/**
 * Ian Sohl & Xin Tong - 2015
 * Copyright (c) 2015      Los Alamos National Security, LLC
 *                         All rights reserved.
 * Legion Image Composition - Ray Trace Rendering Code
 */

#ifndef RENDER_CU
#define RENDER_CU

//#include "cuda.h"
//#include "cuda_runtime.h"
#include "cuda_helper.h"
#include "helper_math.h"


#include "composite.h"

typedef float VolumeType;


typedef struct
{
	float4 m[4];
} float4x4; /**< Matrix Holding form */

//__device__
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

//__device__
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



//__device__
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
	float largest_tmin = fmaxf(fmaxf(tmin.x, tmin.y), tmin.z);
	float smallest_tmax = fminf(fminf(tmax.x, tmax.y), tmax.z);

	*tnear = largest_tmin;
	*tfar = smallest_tmax;
	
	return smallest_tmax > largest_tmin;
}


//__device__
int intersectBoxAlt(MyRay r, float3 boxmin, float3 boxmax, float *tnear, float *tfar){ 
    float tmin = (boxmin.x - r.o.x) / r.d.x; 
    float tmax = (boxmax.x - r.o.x) / r.d.x; 
 
    if (tmin > tmax){
    	float temp = tmin;
    	tmin = tmax;
    	tmax = temp;
    }
 
    float tymin = (boxmin.y - r.o.y) / r.d.y; 
    float tymax = (boxmax.y - r.o.y) / r.d.y; 
 
    if (tymin > tymax){
    	float temp = tymin;
    	tymin = tymax;
    	tymax = temp;
    } 
 
    if ((tmin > tymax) || (tymin > tmax)) 
        return 0; 
 
    if (tymin > tmin) 
        tmin = tymin; 
 
    if (tymax < tmax) 
        tmax = tymax; 
 
    float tzmin = (boxmin.z - r.o.z) / r.d.z; 
    float tzmax = (boxmax.z - r.o.z) / r.d.z; 
 
    if (tzmin > tzmax){
    	float temp = tzmin;
    	tzmin = tzmax;
    	tzmax = temp;
    }
 
    if ((tmin > tzmax) || (tzmin > tmax)) 
        return false; 
 
    if (tzmin > tmin) 
        tmin = tzmin; 
 
    if (tzmax < tmax) 
        tmax = tzmax; 
 
    *tnear = tmin;
    *tfar = tmax;
    return true; 
} 

//__device__
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

//__device__
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
			return 1.0f;
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


//__global__ void
void d_render(int imageW, int imageH,
		int3 boxSize, float3 minBound, float3 maxBound,
		float density, float brightness,
		float transferOffset, float transferScale, 
		float* imgPtr,
		float4x4 invPVMMatrix, float* dataPtr, uint x, uint y)
{
	/**
	 * Kernal renderer for individual ray tracing
	 */

	const int maxSteps = (int)sqrtf(boxSize.x*boxSize.x+boxSize.y*boxSize.y+boxSize.z*boxSize.z)*10;	// The maximum possible number of steps
	const float tstep = 0.1f;				// Distance to step
	const float opacityThreshold = 0.95f;	// Arbitrarily defined alpha cutoff
	

	//uint x = blockIdx.x*blockDim.x + threadIdx.x;	// Current pixel x value
	//uint y = blockIdx.y*blockDim.y + threadIdx.y;	// Current pixel y value
//printf("anmol: x:%d y:%d blockIdx.x:%d blockDim.x:%d threadIdx.x:%d\n",x,y,blockIdx.x, blockDim.x,threadIdx.x);	
//printf("anmol: x:%d y:%d blockIdx.y:%d blockDim.y:%d threadIdx.y:%d\n",x,y,blockIdx.y, blockDim.y,threadIdx.y);	

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
//	if (x==501 && y == 501){
//		printf("invPVM:\n");
//		printf("\t<%f,%f,%f,%f>\n",invPVMMatrix.m[0].x,invPVMMatrix.m[0].y,invPVMMatrix.m[0].z,invPVMMatrix.m[0].w);
//		printf("\t<%f,%f,%f,%f>\n",invPVMMatrix.m[1].x,invPVMMatrix.m[1].y,invPVMMatrix.m[1].z,invPVMMatrix.m[1].w);
//		printf("\t<%f,%f,%f,%f>\n",invPVMMatrix.m[2].x,invPVMMatrix.m[2].y,invPVMMatrix.m[2].z,invPVMMatrix.m[2].w);
//		printf("\t<%f,%f,%f,%f>\n",invPVMMatrix.m[3].x,invPVMMatrix.m[3].y,invPVMMatrix.m[3].z,invPVMMatrix.m[3].w);
//		printf("x:%d, y:%d, u:%f, v:%f\n",x,y,u,v);
//		printf("EyeRay.o:<%f,%f,%f>  EyeRay.d:<%f,%f,%f>\n",eyeRay.o.x,eyeRay.o.y,eyeRay.o.z,eyeRay.d.x,eyeRay.d.y,eyeRay.d.z);
//	}
	int hit = intersectBoxAlt(eyeRay, minBound, maxBound, &tnear, &tfar);
//	float4 cols[] = { 	// Hard-coded transfer function (Fixme)
//			make_float4(0.0, 0.0, 0.0, 0.00),
//			make_float4(0.0, 0.5, 0.0, 0.05),
//			make_float4(0.0, 0.5, 0.0, 0.05),
//			make_float4(0.0, 0.5, 0.0, 0.05),
//			make_float4(0.0, 0.5, 0.0, 0.05),
//			make_float4(0.0, 0.0, 0.5, 0.05),
//			make_float4(0.0, 0.0, 0.5, 0.05),
//			make_float4(0.0, 0.0, 0.5, 0.05),
//			make_float4(0.0, 0.0, 0.5, 0.05),
//			make_float4(0.0, 0.0, 0.5, 0.05),
//			
//			make_float4(0.5, 0.0, 0.0, 0.00),
//			make_float4(0.5, 0.0, 0.0, 0.05),
//			make_float4(0.5, 0.0, 0.0, 0.05),
//			make_float4(0.5, 0.0, 0.0, 0.05),
//			make_float4(0.5, 0.0, 0.0, 0.05),
//			make_float4(0.0, 0.0, 0.5, 0.05),
//			make_float4(0.0, 0.0, 0.5, 0.05),
//			make_float4(0.0, 0.0, 0.5, 0.05),
//			make_float4(0.0, 0.0, 0.5, 0.05),
//			make_float4(0.0, 0.0, 0.5, 0.05),
//			
//			make_float4(0.0, 0.0, 0.0, 0.05),
//	};

	
	float4 cols[] = {make_float4(0,0,0,0),
			make_float4(0,0,0,0),
			make_float4(0,0,0,0),
			make_float4(0,0,0,0),
			make_float4(0,0,0,0),
			make_float4(0,0,0,0),
			make_float4(0,0,0,0),
			make_float4(0,0,0,0),
			make_float4(0,0,0,0),
			make_float4(0,0,0,0),
			make_float4(0,0,0,0),
			make_float4(0,0,0,0),
			make_float4(0,0,0,0),
			make_float4(0,0,0.363636,0),
			make_float4(0,0.459091,0.459091,0),
			make_float4(0,0.490909,0.490909,0.0318182),
			make_float4(0,0.477273,0.477273,0.0378788),
			make_float4(0,0.513636,0.513636,0.0439394),
			make_float4(0,0.531818,0.531818,0.05),
			make_float4(0,0.545455,0.545455,0.05),
			make_float4(0,0.55,0.55,0.05),
			make_float4(0,0.554545,0.554545,0.05),
			make_float4(0,0.554545,0.554545,0.05),
			make_float4(0,0.559091,0.559091,0.05),
			make_float4(0,0.563636,0.563636,0.05),
			make_float4(0,0.563636,0.563636,0.05),
			make_float4(0,0.568182,0.568182,0.0545455),
			make_float4(0,0.568182,0.568182,0.0545455),
			make_float4(0,0.568182,0.568182,0.0545455),
			make_float4(0,0.568182,0.568182,0.0590909),
			make_float4(0,0.568182,0.568182,0.0636364),
			make_float4(0,0.568182,0.568182,0.0681818),
			make_float4(0.00681818,0.568182,0.568182,0.0727273),
			make_float4(0.00909091,0.568182,0.568182,0.0727273),
			make_float4(0.0106061,0.568182,0.568182,0.0727273),
			make_float4(0.0121212,0.563636,0.563636,0.0727273),
			make_float4(0.0136364,0.559091,0.559091,0.0727273),
			make_float4(0.0181818,0.554545,0.554545,0.0727273),
			make_float4(0.0227273,0.540909,0.540909,0.0727273),
			make_float4(0.0272727,0.540909,0.540909,0.0727273),
			make_float4(0.0340909,0.522727,0.522727,0.0636364),
			make_float4(0.0409091,0.518182,0.518182,0.0590909),
			make_float4(0.0439394,0.513636,0.513636,0.0545455),
			make_float4(0.0469697,0.509091,0.509091,0.0545455),
			make_float4(0.05,0.504545,0.504545,0.05),
			make_float4(0.05,0.495455,0.495455,0.05),
			make_float4(0.0772727,0.486364,0.486364,0.0454545),
			make_float4(0.0954545,0.463636,0.463636,0.0409091),
			make_float4(0.131818,0.195455,0.195455,0.0363636),
			make_float4(0.154545,0.177273,0.177273,0.0363636),
			make_float4(0.181818,0.168182,0.168182,0.0409091),
			make_float4(0.190909,0.15,0.15,0.0727273),
			make_float4(0.204545,0.136364,0.136364,0.0727273),
			make_float4(0.227273,0.122727,0.122727,0.0909091),
			make_float4(0.245455,0.122727,0.122727,0.104545),
			make_float4(0.25,0.118182,0.118182,0.113636),
			make_float4(0.263636,0.104545,0.104545,0.131818),
			make_float4(0.268182,0.1,0.1,0.140909),
			make_float4(0.286364,0.0909091,0.0909091,0.154545),
			make_float4(0.290909,0.0681818,0.0681818,0.168182),
			make_float4(0.3,0.0545455,0.0545455,0.172727),
			make_float4(0.313636,0.0409091,0.0409091,0.186364),
			make_float4(0.318182,0.0136364,0.0136364,0.190909),
			make_float4(0.322727,0.00454545,0.00454545,0.204545),
			make_float4(0.331818,0,0,0.209091),
			make_float4(0.336364,0,0,0.218182),
			make_float4(0.340909,0,0,0.231818),
			make_float4(0.345455,0,0,0.240909),
			make_float4(0.347727,0,0,0.243182),
			make_float4(0.35,0,0,0.245455),
			make_float4(0.356061,0,0,0.25),
			make_float4(0.362121,0,0,0.259091),
			make_float4(0.368182,0,0,0.263636),
			make_float4(0.372727,0,0,0.268182),
			make_float4(0.377273,0,0,0.275),
			make_float4(0.390909,0,0,0.281818),
			make_float4(0.394318,0,0,0.286364),
			make_float4(0.397727,0,0,0.286364),
			make_float4(0.401136,0,0,0.295455),
			make_float4(0.404545,0,0,0.3),
			make_float4(0.409091,0,0,0.3),
			make_float4(0.412121,0,0,0.3),
			make_float4(0.415152,0,0,0.3),
			make_float4(0.418182,0,0,0.3),
			make_float4(0.418182,0,0,0.3),
			make_float4(0.418182,0,0,0.3),
			make_float4(0.418182,0,0,0.3),
			make_float4(0.418182,0,0,0.3),
			make_float4(0.418182,0,0,0.295455),
			make_float4(0.418182,0,0,0.295455),
			make_float4(0.415909,0,0,0.295455),
			make_float4(0.413636,0,0,0.295455),
			make_float4(0.409091,0,0,0.295455),
			make_float4(0.406061,0,0,0.290909),
			make_float4(0.40303,0,0,0.290909),
			make_float4(0.4,0,0,0.286364),
			make_float4(0.39697,0,0,0.281818),
			make_float4(0.393939,0,0,0.277273),
			make_float4(0.390909,0,0,0.263636),
			make_float4(0.386364,0,0,0.259091),
			make_float4(0.381818,0,0,0.254545),
			make_float4(0.377273,0,0,0.245455),
			make_float4(0.370455,0,0,0.236364),
			make_float4(0.363636,0,0,0.231818),
			make_float4(0.359091,0,0,0.222727),
			make_float4(0.354545,0,0,0.195455),
			make_float4(0.347727,0,0,0.181818),
			make_float4(0.340909,0,0,0.168182),
			make_float4(0.336364,0,0,0.159091),
			make_float4(0.331818,0,0,0.15),
			make_float4(0.324242,0,0,0.140909),
			make_float4(0.316667,0,0,0.122727),
			make_float4(0.309091,0,0,0.113636),
			make_float4(0.304545,0,0,0.1),
			make_float4(0.3,0,0,0.0863636),
			make_float4(0.295455,0,0,0.0818182),
			make_float4(0.288636,0,0,0.0727273),
			make_float4(0.281818,0,0,0.0681818),
			make_float4(0.279545,0,0,0.0636364),
			make_float4(0.277273,0,0,0.0590909),
			make_float4(0.270455,0,0,0.0568182),
			make_float4(0.263636,0,0,0.0545455),
			make_float4(0.256818,0,0,0.0545455),
			make_float4(0.25,0,0,0.0545455),
			make_float4(0.245455,0,0,0.0545455),
			make_float4(0.240909,0,0,0.0545455),
			make_float4(0.234848,0,0,0.0545455),
			make_float4(0.228788,0,0,0.0545455),
			make_float4(0.222727,0,0,0.0554545),
			make_float4(0.213636,0,0,0.0563636),
			make_float4(0.204545,0,0,0.0572727),
			make_float4(0.193182,0,0,0.0581818),
			make_float4(0.181818,0,0,0.0590909),
			make_float4(0.118182,0,0,0.06),
			make_float4(0.111364,0,0,0.0609091),
			make_float4(0.104545,0,0,0.0618182),
			make_float4(0.0977273,0,0,0.0627273),
			make_float4(0.0909091,0,0,0.0636364),
			make_float4(0.080303,0,0,0.0645454),
			make_float4(0.069697,0,0,0.0654545),
			make_float4(0.0590909,0,0,0.0663636),
			make_float4(0.0511364,0,0,0.0672727),
			make_float4(0.0431818,0.00454545,0,0.0681818),
			make_float4(0.0352273,0.0727273,0,0.07),
			make_float4(0.0272727,0.109091,0,0.0718182),
			make_float4(0.0204545,0.131818,0,0.0736364),
			make_float4(0.0136364,0.140909,0,0.0754545),
			make_float4(0.00681818,0.163636,0,0.0772727),
			make_float4(0,0.172727,0,0.0784091),
			make_float4(0,0.195455,0,0.0795455),
			make_float4(0,0.213636,0,0.0806818),
			make_float4(0,0.313636,0,0.0818182),
			make_float4(0,0.343182,0,0.0836364),
			make_float4(0,0.372727,0,0.0854545),
			make_float4(0,0.413636,0,0.0872727),
			make_float4(0,0.440909,0,0.0890909),
			make_float4(0,0.459091,0,0.0909091),
			make_float4(0,0.509091,0,0.0954545),
			make_float4(0,0.563636,0,0.1),
			make_float4(0,0.618182,0,0.109091),
			make_float4(0,0.659091,0,0.115152),
			make_float4(0,0.677273,0,0.121212),
			make_float4(0,0.681818,0,0.127273),
			make_float4(0,0.681818,0,0.131818),
			make_float4(0,0.681818,0,0.136364),
			make_float4(0,0.681818,0,0.138636),
			make_float4(0,0.681818,0,0.140909),
			make_float4(0,0.681818,0,0.145455),
			make_float4(0,0.681818,0,0.15),
			make_float4(0,0.681818,0,0.159091),
			make_float4(0,0.680682,0,0.161364),
			make_float4(0,0.679545,0,0.163636),
			make_float4(0,0.678409,0,0.163636),
			make_float4(0,0.677273,0,0.163636),
			make_float4(0,0.668182,0,0.168182),
			make_float4(0,0.659091,0,0.172727),
			make_float4(0,0.65,0,0.172727),
			make_float4(0,0.637879,0,0.172727),
			make_float4(0,0.625758,0,0.181818),
			make_float4(0,0.613636,0,0.186364),
			make_float4(0,0.6,0,0.190909),
			make_float4(0,0.586364,0,0.195455),
			make_float4(0,0.577273,0,0.204545),
			make_float4(0,0.568182,0,0.204545),
			make_float4(0,0.557955,0,0.204545),
			make_float4(0,0.547727,0,0.204545),
			make_float4(0,0.5375,0,0.204545),
			make_float4(0,0.527273,0,0.204545),
			make_float4(0,0.518182,0,0.204545),
			make_float4(0,0.509091,0,0.204545),
			make_float4(0,0.5,0,0.204545),
			make_float4(0,0.489394,0,0.195455),
			make_float4(0,0.478788,0,0.186364),
			make_float4(0,0.468182,0,0.181818),
			make_float4(0,0.454545,0,0.177273),
			make_float4(0,0.445455,0,0.172727),
			make_float4(0,0.436364,0,0.168182),
			make_float4(0,0.430303,0,0.163636),
			make_float4(0,0.424242,0,0.154545),
			make_float4(0,0.418182,0,0.15),
			make_float4(0,0.412121,0,0.145455),
			make_float4(0,0.406061,0,0.142424),
			make_float4(0,0.4,0,0.139394),
			make_float4(0,0.390909,0,0.136364),
			make_float4(0,0.381818,0,0.131818),
			make_float4(0,0.372727,0,0.127273),
			make_float4(0,0.365152,0,0.125),
			make_float4(0,0.357576,0,0.122727),
			make_float4(0,0.35,0,0.120455),
			make_float4(0,0.327273,0,0.118182),
			make_float4(0,0.318182,0,0.115909),
			make_float4(0,0.309091,0,0.113636),
			make_float4(0,0.3,0,0.111364),
			make_float4(0,0.290909,0,0.109091),
			make_float4(0,0.288636,0,0.104545),
			make_float4(0,0.286364,0,0.1),
			make_float4(0,0.272727,0,0.0977273),
			make_float4(0,0.268182,0,0.0909091),
			make_float4(0,0.259091,0,0.0863636),
			make_float4(0,0.254545,0,0.0840909),
			make_float4(0,0.25,0,0.0818182),
			make_float4(0,0.240909,0,0.0818182),
			make_float4(0,0.238636,0,0.0795455),
			make_float4(0,0.236364,0,0.0772727),
			make_float4(0,0.227273,0,0.075),
			make_float4(0,0.222727,0,0.0727273),
			make_float4(0,0.218182,0,0.0727273),
			make_float4(0,0.213636,0,0.0727273),
			make_float4(0,0.204545,0,0.0727273),
			make_float4(0,0.195455,0,0.0727273),
			make_float4(0,0.186364,0,0.0727273),
			make_float4(0,0.177273,0,0.0727273),
			make_float4(0,0.168182,0,0.0727273),
			make_float4(0,0.161364,0,0.0727273),
			make_float4(0,0.154545,0,0.0727273),
			make_float4(0,0.140909,0,0.0727273),
			make_float4(0,0.134091,0,0.0727273),
			make_float4(0,0.127273,0,0.0712121),
			make_float4(0,0.122727,0,0.069697),
			make_float4(0,0.118182,0,0.0681818),
			make_float4(0,0.115909,0,0.0681818),
			make_float4(0,0.113636,0,0.0681818),
			make_float4(0,0.104545,0,0.0681818),
			make_float4(0,0.1,0,0.0681818),
			make_float4(0,0.0954545,0,0.0681818),
			make_float4(0,0.0931818,0,0.0681818),
			make_float4(0,0.0909091,0,0.0681818),
			make_float4(0,0.0863636,0,0.0681818),
			make_float4(0,0.0772727,0,0.0636364),
			make_float4(0,0.075,0,0.0636364),
			make_float4(0,0.0727273,0,0.0590909),
			make_float4(0,0.0704545,0,0.0590909),
			make_float4(0,0.0681818,0,0.0545455),
			make_float4(0,0.0636364,0,0.05),
			make_float4(0,0.0636364,0,0.05),
			make_float4(0,0.0545455,0,0.0136364),
	};
	
	
	if (hit){
		if (tnear < 0.0f) tnear = 0.0f;     // clamp to near plane
		
//		drawPixel(imgPtr,x,y,imageW,0.5,0.0,0.0,0.5);
		// march along ray from front to back, accumulating color
		float4 sum = make_float4(0.0f,0.0f,0.0f,0.0f);
		float t = tnear;
		float3 pos = (eyeRay.o + eyeRay.d*tnear - minBound);
		float3 step = eyeRay.d*tstep;

		for (int i=0; i<maxSteps; i++){
//			if(pos.x< boxMax.x && pos.x >= boxMin.x && pos.y< boxMax.y && pos.y >= boxMin.y && pos.z< boxMax.z && pos.z >= boxMin.z){
				float sample = interpolate(dataPtr,pos,boxSize);
				float4 col;
				col = cols[(int)floor(sample*256)];
//				if(sample>0) col.w = 0.05;
//				if(sample<0.01) col.w = 0;
//				else if(sample < 0.05) col.x = 0.5;
//				else if(sample < 0.2) col.w = 0.0;
//				else if(sample < 0.3) col.y = 0.5;
//				else if(sample < 0.4) col.z = 0.5;
//				else if(sample > 0.99) col.w = 0;
//				else col.w = 0.0;
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
				
//			}
//			sum.x += 0.1

			t += tstep;
			if (t > tfar) break;
			pos += step;
		}

		
		drawPixel(imgPtr,x,y,imageW,(float)sum.x*brightness,(float)sum.y*brightness,(float)sum.z*brightness,(float)sum.w);
	}
}


//__host__
int iDivUp(int a, int b){
	/**
	 * Integer division with rounding up
	 */
	return (a % b != 0) ? (a / b + 1) : (a / b);
}

//__host__
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
	float density = 0.1f;			// Arbitrary defined constants
	float brightness = 0.9f;		// 	(should be moved into metadata)
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

//	dim3 blockSize = dim3(16,16);	// Define kernal execution block size
//	dim3 gridSize = dim3(iDivUp(width, blockSize.x), iDivUp(height, blockSize.y)); // Number of pixels per block


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
	float3 minBound = make_float3(tmpimg.partition.xmin,tmpimg.partition.ymin,tmpimg.partition.zmin);
	float3 maxBound = make_float3(tmpimg.partition.xmax,tmpimg.partition.ymax,tmpimg.partition.zmax);

	printf("Box:<%d,%d,%d>\n",boxSize.x,boxSize.y,boxSize.z);
	printf("Min:<%f,%f,%f> Max:<%f,%f,%f>\n",minBound.x,minBound.y,minBound.z,maxBound.x,maxBound.y,maxBound.z);
	
//	float3 minBound = make_float3(0.0,0.0,0.0);
//	float3 maxBound = make_float3(5.0,5.0,5.0);

//	d_render<<<gridSize, blockSize>>>(width,height,boxSize,minBound,maxBound,density,brightness,transferOffset,transferScale,imgPtr,invPVMMatrix, dataPtr);
for(int i =0;i<width;i++) {
	for(int j=0;j<height;j++)
		d_render(width,height,boxSize,minBound,maxBound,density,brightness,transferOffset,transferScale,imgPtr,invPVMMatrix, dataPtr, i, j);
}
//cudaDeviceSynchronize();
}
#endif
