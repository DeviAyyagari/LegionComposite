/**
 * Ian Sohl & Xin Tong - 2015
 * Copyright (c) 2015      Los Alamos National Security, LLC
 *                         All rights reserved.
 * Legion Image Composition - Ray Trace Rendering Code
 */

#ifndef RENDER_CU
#define RENDER_CU

#include "cuda.h"
#include "cuda_runtime.h"
#include "cuda_helper.h"
#include "helper_math.h"


#include "composite.h"

typedef float VolumeType;


typedef struct
{
	float4 m[4];
} float4x4; /**< Matrix Holding form */

__device__
float4 mul(const float4x4 &M, const float4 &v){
	/**
	 * Multiply a 4x4 Matrix with a 1x4 vector
	 */
	float4 r;
	r.w = dot(v, M.m[3]);
	r.x = dot(v, M.m[0]);
	r.y = dot(v, M.m[1]);
	r.z = dot(v, M.m[2]);

	return r;
}

__device__
float4 divW(float4 v){
	/**
	 * Divide a 4-vector by it's homogeneous coordinate
	 */
	float invW = 1 / v.w;
	return(make_float4(v.x * invW, v.y * invW, v.z * invW, 1.0f));
}



struct MyRay{
	float3 o;   /**< Origin Point */
	float3 d;   /**< Direction Vector */
}; /**< Ray vector representation */



__device__
int intersectBox(MyRay r, float3 boxmin, float3 boxmax, float *tnear, float *tfar){
	/**
	 * Check if a ray intersects with the data partition
	 */
	// compute intersection of ray with all six bbox planes
	float3 invR = make_float3(1.0f) / r.d;
	float3 tbot = invR * (boxmin - r.o);
	float3 ttop = invR * (boxmax - r.o);

	// re-order intersections to find smallest and largest on each axis
	float3 tmin = fminf(ttop, tbot);
	float3 tmax = fmaxf(ttop, tbot);

	// find the largest tmin and the smallest tmax
	float largest_tmin = fmaxf(fmaxf(tmin.x, tmin.y), fmaxf(tmin.x, tmin.z));
	float smallest_tmax = fminf(fminf(tmax.x, tmax.y), fminf(tmax.x, tmax.z));

	*tnear = largest_tmin;
	*tfar = smallest_tmax;
	
	return smallest_tmax > largest_tmin;
}


__device__
void drawPixel(float* imgPtr, int x, int y, int imageW, float r, float g, float b, float a){
	/**
	 * Populate a Legion region with a particular pixel
	 */
	int writePoint = (y*imageW+x)*4; // X->Y ordering
	imgPtr[writePoint++] = r;
	imgPtr[writePoint++] = g;
	imgPtr[writePoint++] = b;
	imgPtr[writePoint] = a;
}

__device__
float interpolate(float* dataPtr, float3 pos, int3 partitionSize){
	/**
	 * Replicate Texture functionality with a trilinear interpolant
	 */
	int3 originPoint = make_int3(floor(pos.x),floor(pos.y),floor(pos.z)); 		// Find the corner of the box the point is in
	float3 point = pos-make_float3(originPoint.x,originPoint.y,originPoint.z);	// Find the location of the point within the box
	float3 complement = make_float3(1,1,1)-point;								// Compute the distance to the opposite corner
	auto getPoint = [&](int x, int y, int z){									// Lambda function: Get a particular point from volumetric data
	int3 p = originPoint + make_int3(x,y,z);									//		Only works on integer values
		if(p.x>=partitionSize.x ||		// 	Make sure the point is in the array
				p.y>=partitionSize.y || 
				p.z>=partitionSize.z )
			return 0.0f;
		else
			return dataPtr[p.z*partitionSize.y*partitionSize.x+p.y*partitionSize.x+p.x];	// Get the point from legion X->Y->Z
	};
	float sample = 	getPoint(0,0,0) *	complement.x *	complement.y *	complement.z +	// Standard trilinear interpolant
					getPoint(1,0,0) *	point.x		 *	complement.y *	complement.z +
					getPoint(0,1,0) *	complement.x *	point.y		 *	complement.z +
					getPoint(0,0,1) *	complement.x *	complement.y *	point.z 	 +
					getPoint(1,0,1) *	point.x		 *	complement.y *	point.z		 +
					getPoint(0,1,1) *	complement.x *	point.y		 *	point.z		 +
					getPoint(1,1,0) *	point.x		 *	point.y		 *	complement.z +
					getPoint(1,1,1) *	point.x		 *	point.y		 *	point.z;
	return sample;
}


__global__ void
d_render(int imageW, int imageH,
		int3 boxSize, float3 minBound, float3 maxBound,
		float density, float brightness,
		float transferOffset, float transferScale, 
		float* imgPtr,
		float4x4 invPVMMatrix, float* dataPtr)
{
	/**
	 * Kernal renderer for individual ray tracing
	 */

	const int maxSteps = (int)sqrtf(boxSize.x*boxSize.x+boxSize.y*boxSize.y+boxSize.z*boxSize.z);	// The maximum possible number of steps
	const float tstep = 0.1f;				// Distance to step
	const float opacityThreshold = 0.95f;	// Arbitrarily defined alpha cutoff

	const float3 boxMin = make_float3(0, 0, 0); 	// Minimum bounds of data partition
	const float3 boxMax = make_float3(boxSize.x - 1,				// Maximum bound of partition
			boxSize.y - 1,
			boxSize.z - 1);
	
	const float3 boxScale = make_float3(abs(maxBound.x - minBound.x) / boxSize.x,
			abs(maxBound.y - minBound.y) / boxSize.y,
			abs(maxBound.z - minBound.z) / boxSize.z);
	
	uint x = blockIdx.x*blockDim.x + threadIdx.x;	// Current pixel x value
	uint y = blockIdx.y*blockDim.y + threadIdx.y;	// Current pixel y value
	

	if ((x >= imageW) || (y >= imageH)) return;


	drawPixel(imgPtr,x,y,imageW,0.0f,0.0f,0.0f,0.0f);	// Fill pixel with blank

	float u = (x / (float)imageW)*2.0f-1.0f;			// Get the image space coordinates
	float v = (y / (float)imageH)*2.0f-1.0f;

	//unproject eye ray from clip space to object space
	//unproject: http://gamedev.stackexchange.com/questions/8974/how-can-i-convert-a-mouse-click-to-a-ray
	MyRay eyeRay;
	eyeRay.o = make_float3(divW(mul(invPVMMatrix, make_float4(u, v, 2.0f, 1.0f))));
	float3 eyeRay_t = make_float3(divW(mul(invPVMMatrix, make_float4(u, v, -1.0f, 1.0f))));
	eyeRay.d = normalize(eyeRay_t - eyeRay.o);

	// find intersection with box
	float tnear, tfar;
	int hit = intersectBox(eyeRay, minBound, maxBound, &tnear, &tfar);
	float4 cols[] = { 	// Hard-coded transfer function (Fixme)
			make_float4(0.0, 0.5, 0.0, 0.05),
			make_float4(0.0, 0.0, 0.5, 0.05),
			make_float4(0.5, 0.0, 0.0, 0.05),
			make_float4(0.0, 0.0, 0.5, 0.05),
			make_float4(0.0, 0.0, 0.5, 0.05),
	};

	if (hit){
		if (tnear < 0.0f) tnear = 0.0f;     // clamp to near plane
		
//		drawPixel(imgPtr,x,y,imageW,0.5,0.0,0.0,0.1);
//		return;
		// march along ray from front to back, accumulating color
		float4 sum = make_float4(0.0f,0.0f,0.0f,0.0f);
		float t = tnear;
		float3 pos = (eyeRay.o + eyeRay.d*tnear + minBound) / boxScale;
		float3 step = eyeRay.d*tstep;

		for (int i=0; i<maxSteps; i++){
			if(pos.x< boxMax.x && pos.x >= boxMin.x && pos.y< boxMax.y && pos.y >= boxMin.y && pos.z< boxMax.z && pos.z >= boxMin.z){
				float sample = interpolate(dataPtr,pos,boxSize);
				float4 col;
				col = cols[(int)floor(sample*5)];
				col.w *= density;

				// "under" operator for back-to-front blending
				//sum = lerp(sum, col, col.w);

				// pre-multiply alpha
				col.x *= col.w;
				col.y *= col.w;
				col.z *= col.w;
				// "over" operator for front-to-back blending
				sum += col*(1.0f - sum.w);

				// exit early if opaque
				if (sum.w > opacityThreshold)
					break;
				
			}

			t += tstep;
			if (t > tfar) break;
			pos += step;
		}

		sum *= brightness;
		
		drawPixel(imgPtr,x,y,imageW,(float)sum.x,(float)sum.y,(float)sum.z,(float)sum.w);
	}
}


__host__
int iDivUp(int a, int b){
	/**
	 * Integer division with rounding up
	 */
	return (a % b != 0) ? (a / b + 1) : (a / b);
}

__host__
void create_task(const Task *task,
		const std::vector<PhysicalRegion> &regions,
		LegionRuntime::HighLevel::Context ctx, HighLevelRuntime *runtime){
	/**
	 * Image rendering task
	 */
	
	assert(regions.size()==2);

	Image tmpimg = *((Image*)task->args);	// Metadata for current render
	RegionAccessor<AccessorType::Generic, float> dataAccessor = regions[1].get_field_accessor(FID_VAL).typeify<float>(); // Accessor for data
	RegionAccessor<AccessorType::Generic, float> imgAccessor = regions[0].get_field_accessor(FID_VAL).typeify<float>();	// And image
	float density = 0.05f;			// Arbitrary defined constants
	float brightness = 1.0f;		// 	(should be moved into metadata)
	float transferOffset = 0.0f;
	float transferScale = 1.0f;
	int width = tmpimg.width;			// Get total image size
	int height = tmpimg.height;


	float4x4 invPVMMatrix; // Copy over inverse PV Matrix from metadata
	for(int i = 0; i < 4; ++i){
		invPVMMatrix.m[i].x = tmpimg.invPVM[4*i+0];
		invPVMMatrix.m[i].y = tmpimg.invPVM[4*i+1];
		invPVMMatrix.m[i].z = tmpimg.invPVM[4*i+2];
		invPVMMatrix.m[i].w = tmpimg.invPVM[4*i+3];
	}

	dim3 blockSize = dim3(16,16);	// Define kernal execution block size
	dim3 gridSize = dim3(iDivUp(width, blockSize.x), iDivUp(height, blockSize.y)); // Number of pixels per block


	Domain dataDomain = runtime->get_index_space_domain(ctx,regions[1].get_logical_region().get_index_space());
	Rect<1> dataRect = dataDomain.get_rect<1>();	// Get data size domain
	Rect<1> dataSubRect;							// Empty filler rectangle
	ByteOffset dataOffsets[1];						// Byte Offset object
	float* dataPtr = dataAccessor.raw_rect_ptr<1>(dataRect,dataSubRect,dataOffsets); // Get raw framebuffer pointers
	

	Domain imgDomain = runtime->get_index_space_domain(ctx,regions[0].get_logical_region().get_index_space());
	Rect<1> imgRect = imgDomain.get_rect<1>();
	Rect<1> imgSubRect;
	ByteOffset imgOffsets[1];
	float* imgPtr = imgAccessor.raw_rect_ptr<1>(imgRect,imgSubRect,imgOffsets);	// For output image as well

//	int3 lowerBound = make_int3(tmpimg.partition.xmin, tmpimg.partition.ymin, tmpimg.partition.zmin);
//	int3 upperBound = make_int3(tmpimg.partition.xmax,tmpimg.partition.ymax,tmpimg.partition.zmax);


	int3 boxSize = make_int3(tmpimg.partition.datx,tmpimg.partition.daty,tmpimg.partition.datz);
	float3 minBound = make_float3(tmpimg.partition.xmin + 2,tmpimg.partition.ymin + 2,tmpimg.partition.zmin + 2)*3;
	float3 maxBound = make_float3(tmpimg.partition.xmax + 2,tmpimg.partition.ymax + 2,tmpimg.partition.zmax + 2)*3;
//	
//	printf("Min:<%f,%f,%f> Max:<%f,%f,%f>\n",minBound.x,minBound.y,minBound.z,maxBound.x,maxBound.y,maxBound.z);
	
//	float3 minBound = make_float3(0.0,0.0,0.0);
//	float3 maxBound = make_float3(5.0,5.0,5.0);

	d_render<<<gridSize, blockSize>>>(width,height,boxSize,minBound,maxBound,density,brightness,transferOffset,transferScale,imgPtr,invPVMMatrix, dataPtr);

	cudaDeviceSynchronize();
}
#endif