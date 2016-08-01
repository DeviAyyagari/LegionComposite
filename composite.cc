/**
 * Ian Sohl - 2015
 * Copyright (c) 2015      Los Alamos National Security, LLC
 *                         All rights reserved.
 * Legion Image Composition - Main Code
 */

#include <cstdio>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <stdlib.h>     /* srand, rand */
#include <vector>
#include <time.h>       /* time */
#include "glm/glm/glm.hpp"
#include "glm/glm/vec3.hpp"
#include "glm/glm/vec4.hpp"
#include "glm/glm/mat4x4.hpp"
#include "glm/glm/gtc/type_ptr.hpp"
#include "glm/glm/gtc/matrix_transform.hpp"
#include "glm/glm/gtx/string_cast.hpp"
#include <boost/gil/extension/io/png_io.hpp>
#include "composite.h"
#include "DataMgr.h"
#include <algorithm> 

using namespace LegionRuntime::HighLevel;
using namespace LegionRuntime::Accessor;
using namespace std;

void create_task(const Task *task,
		const std::vector<PhysicalRegion> &regions,
		Context ctx, HighLevelRuntime *runtime);

int sqr(int x){
	/**
	 * Returns an integer square value
	 */
	return x*x;
}

int subCoordTransform(int z, int width, int yoffset){
	/**
	 * Converts a index in total region space to an index in sub-region space
	 */
	return z - yoffset * width;
}



float exclusion(float a){
	/**
	 * Compositing swap function
	 * Computes 1 - value
	 */
	return 1.0-a;
}

float passthrough(float a){
	/**
	 * Compositing swap function
	 * Passes through same value
	 */
	return a;
}

float pass1(float a){
	/**
	 * Compositing swap function
	 * Returns 1.0;
	 */
	return 1.0;
}

float pass0(float a){
	/**
	 * Compositing swap function
	 * Returns 0.0;
	 */
	return 0.0;
}

void cpu_draw_task(const Task *task,
		const std::vector<PhysicalRegion> &regions,
		Context ctx, HighLevelRuntime *runtime){
	Image img = *((Image*)task->args);	// Task metadata
	PhysicalRegion imgPhysicalRegion = regions[0];
	imgPhysicalRegion.wait_until_valid();
	Domain outDomain = runtime->get_index_space_domain(ctx,regions[0].get_logical_region().get_index_space());
	Rect<1> outRect = outDomain.get_rect<1>();			// Get the size of the return image
	RegionAccessor<AccessorType::Generic,float> outputAccessor = regions[0].get_field_accessor(FID_VAL).typeify<float>();
	srand(img.randomseed);
	int primary = rand() % 3;
	float c1 = 0.0;
	float c2 = 0.0;
	float c3 = 0.0;
	if(primary==0) c1 = 1.0;
	if(primary==1) c2 = 1.0;
	if(primary==2) c3 = 1.0;
	float c4 = 1.0;
	for(int y = img.partition.ymin; y < img.partition.ymax; ++y){
		for(int x = img.partition.xmin; x < img.partition.xmax; ++x){
			int p = (y * img.width + x) * 4;
			if(p > outRect.hi.x[0]) assert(false);
			if(p < outRect.lo.x[0]) assert(false);
			outputAccessor.write(DomainPoint::from_point<1>(Point<1>(p)),c1);
			outputAccessor.write(DomainPoint::from_point<1>(Point<1>(p+1)),c2);
			outputAccessor.write(DomainPoint::from_point<1>(Point<1>(p+2)),c3);
			outputAccessor.write(DomainPoint::from_point<1>(Point<1>(p+3)),c4);
		}
	}

}

float* loadRawFile(const char *filename, int nx, int ny, int nz){
    FILE *fp = fopen(filename, "rb");
    size_t size = nx * ny * nz * sizeof(unsigned char);
    int nCells = nx * ny * nz;
    if (!fp){
        fprintf(stderr, "Error opening file '%s'\n", filename);
        return 0;
    }
    unsigned char * temp = (unsigned char *)malloc(size);
    float * data = (float*)malloc(nCells*sizeof(float));
    size_t read = fread(temp, 1, size, fp);
    printf("Read '%s', %zu bytes\n", filename, read);
	fflush(stdout);
    fclose(fp);
    for(int i = 0; i < nCells; i++) {
        data[i] = (float)temp[i] / 255;
    }
    return data;
}

void create_interface_task(const Task *task,
		const std::vector<PhysicalRegion> &regions,
		Context ctx, HighLevelRuntime *runtime){
	Image img = *((Image*)task->args);	// Task metadata
	char filename[100];
	sprintf(filename,"heptane_%d_%d_%d.raw",img.i,img.j,img.k);
	//sprintf(filename,"heptane.raw");
	float *volume = loadRawFile(filename, img.partition.datx, img.partition.daty, img.partition.datz);

	Rect<1> dataBound = Rect<1>(0,img.partition.datx*img.partition.daty*img.partition.datz-1);	// Indexing the region used to hold the data (linearized)
	IndexSpace dataIndexSpace = runtime->create_index_space(ctx, Domain::from_rect<1>(dataBound)); //Create the Index Space (1 index per voxel)
	FieldSpace dataFieldSpace = runtime->create_field_space(ctx);	// Simple field space
	{
		FieldAllocator allocator = runtime->create_field_allocator(ctx,dataFieldSpace);
		allocator.allocate_field(sizeof(float),FID_VAL);			// Only requires one field
	}
	LogicalRegion dataLogicalRegion = runtime->create_logical_region(ctx,dataIndexSpace,dataFieldSpace); // Create the Logical Region
	{																// Populate the region
		RegionRequirement req(dataLogicalRegion, WRITE_DISCARD, EXCLUSIVE, dataLogicalRegion); // Filling requirement
		req.add_field(FID_VAL);
		InlineLauncher dataInlineLauncher(req);						// Inline launchers are simple


		PhysicalRegion dataPhysicalRegion = runtime->map_region(ctx,dataInlineLauncher);	// Map to a physical region
		dataPhysicalRegion.wait_until_valid();						// Should be pretty fast

		RegionAccessor<AccessorType::Generic, float> dataAccessor = dataPhysicalRegion.get_field_accessor(FID_VAL).typeify<float>();
		// The GPU's tested with have much better single precision performance. If this is changed, the renderer needs to be modified, too
		int i = 0;
//		float max = 0;
//		float min = 1000000;
		for(GenericPointInRectIterator<1> pir(dataBound); pir; pir++){	// Step through the data and write to the physical region
//			if(volume[i] > max) max = volume[i];
//			if(volume[i] < min && volume[i]!=0) min = volume[i];
			dataAccessor.write(DomainPoint::from_point<1>(pir.p),volume[i++]); // Same order as data: X->Y->Z
		}
//		cout << "Max: " << max << " and Min: " << min << endl;
		runtime->unmap_region(ctx,dataPhysicalRegion);					// Free up resources
	}
	TaskLauncher loadLauncher(CREATE_TASK_ID, TaskArgument(&img,sizeof(img)));	// Spawn the renderer task
	loadLauncher.add_region_requirement(RegionRequirement(regions[0].get_logical_region(),WRITE_DISCARD,EXCLUSIVE,regions[0].get_logical_region()));
	loadLauncher.add_field(0,FID_VAL);		// Output Image as second region
	loadLauncher.add_region_requirement(RegionRequirement(dataLogicalRegion,READ_ONLY,EXCLUSIVE,dataLogicalRegion));
	loadLauncher.add_field(1,FID_VAL);		// Input Data as third region
	runtime->execute_task(ctx,loadLauncher);	// Launch and terminate render task
	
}

void composite(RegionAccessor<AccessorType::Generic, float> input1, RegionAccessor<AccessorType::Generic, float> input2, RegionAccessor<AccessorType::Generic, float> O, int start, int stop, tripleArgument ta, float (*FA)(float), float (*FB)(float)){

	int linewidth = 4 * ta.co.width;	// Number of indices on one scanline

	int voffset1 = (ta.co1.miny - ta.co.miny);	// Vertical offset of region 1 from the top
	int maxbound1 = linewidth * (ta.co1.maxy - ta.co1.miny+1);	// Vertical size of region 1

	int voffset2 = (ta.co2.miny - ta.co.miny); 	// Vertical offset of region 2 from the top
	int maxbound2 = linewidth * (ta.co2.maxy - ta.co2.miny+1);	// Vertical size of region 2


	Point<1> input1Point(subCoordTransform(start,linewidth,voffset1));	// Scanning current point in region 1
	Point<1> input2Point(subCoordTransform(start,linewidth,voffset2));	// Scanning point for region 2
	Point<1> OPoint(start); 											// Scanning point for output region


	for(int i = start; i < stop;i+=4){	//  Step through all of the data in groups of 4 (RGBA pixels)

		bool range1 = input1Point.x[0] >= 0 && input1Point.x[0] < maxbound1; // Check if region 1 is in range
		bool range2 = input2Point.x[0] >= 0 && input2Point.x[0] < maxbound2; // Check if region 2 is in range

#ifdef ISOSURFACE
		bool transparent1 = true;
		bool transparent2 = true;
		for(int j = 0; j < 3; ++j){
			if(range1 ? input1.read(DomainPoint::from_point<1>(input1Point)) : 0.0!=0.0){
				transparent1 = false;
			}
			if(range2 ? input2.read(DomainPoint::from_point<1>(input2Point)) : 0.0!=0.0){
				transparent2 = false;
			}
			input1Point.x[0]++;
			input2Point.x[0]++;
		}
		float alphaA = 0.0;
		float alphaB = 0.0;
		if(!transparent1){
			alphaA = range1 ? input1.read(DomainPoint::from_point<1>(input1Point)) : 0.0;
		}
		if(!transparent2){
			alphaB = range2 ? input2.read(DomainPoint::from_point<1>(input2Point)) : 0.0;
		}
#else
		input1Point.x[0]+=3;	// Increment by 3 to get the alpha value
		input2Point.x[0]+=3;
		float alphaA = range1 ? input1.read(DomainPoint::from_point<1>(input1Point)) : 0.0; // If in range, get the alpha value
		float alphaB = range2 ? input2.read(DomainPoint::from_point<1>(input2Point)) : 0.0; //		Otherwise alpha is 0.0
#endif
		input1Point.x[0]-=3;	// Step back in place for the red value
		input2Point.x[0]-=3;
		float alphaC = alphaA*FA(alphaB)+alphaB*FB(alphaA); // Compute the output alpha
		if(alphaC!=0){			// If there is a non-zero alpha
			for(int j = 0; j < 3; ++j){	// For each of R, G, B
				float A = range1 ? input1.read(DomainPoint::from_point<1>(input1Point)) : 0.0; // Get the values from the input
				float B = range2 ? input2.read(DomainPoint::from_point<1>(input2Point)) : 0.0;
				O.write(DomainPoint::from_point<1>(OPoint),(A*alphaA*FA(alphaB)+B*alphaB*FB(alphaA))/alphaC); // Compute composite and write
				input1Point.x[0]++; // Step to next point
				input2Point.x[0]++;
				OPoint.x[0]++;
			}
		}
		else{	// If Alpha is zero
			for(int j = 0; j < 3; ++j){
				O.write(DomainPoint::from_point<1>(OPoint),0.0); // Fill RGB with zeros
				OPoint.x[0]++;
			}
			input1Point.x[0]+=3;
			input2Point.x[0]+=3;
		}
		O.write(DomainPoint::from_point<1>(OPoint),alphaC); // Write output alpha
		OPoint.x[0]++;		// Increment for next pixel
		input1Point.x[0]++;
		input2Point.x[0]++;
	}	


	/**
	 * Generic alpha compositing that can be called with adjustable functions for different operations
	 */
	// cout<<"composite called!"; 
	/*Point<1> input1Point(start);
	Point<1> input2Point(start);
	Point<1> OPoint(start); 											// Scanning point for output region


	for(int i = start; i < stop;i+=4){	//  Step through all of the data in groups of 4 (RGBA pixels)
		input1Point.x[0]+=3;	// Increment by 3 to get the alpha value
		input2Point.x[0]+=3;
		float alphaA = input1.read(DomainPoint::from_point<1>(input1Point));
		float alphaB = input2.read(DomainPoint::from_point<1>(input2Point));
		input1Point.x[0]-=3;	// Step back in place for the red value
		input2Point.x[0]-=3;
		float alphaC = alphaA*FA(alphaB)+alphaB*FB(alphaA); // Compute the output alpha
		if(alphaC!=0){			// If there is a non-zero alpha
			for(int j = 0; j < 3; ++j){	// For each of R, G, B
				float A = input1.read(DomainPoint::from_point<1>(input1Point));
				float B = input2.read(DomainPoint::from_point<1>(input2Point));
				O.write(DomainPoint::from_point<1>(OPoint),(A*alphaA*FA(alphaB)+B*alphaB*FB(alphaA))/alphaC); // Compute composite and write
				input1Point.x[0]++; // Step to next point
				input2Point.x[0]++;
				OPoint.x[0]++;
			}
		}
		else{	// If Alpha is zero
			for(int j = 0; j < 3; ++j){
				O.write(DomainPoint::from_point<1>(OPoint),0.0); // Fill RGB with zeros
				OPoint.x[0]++;
			}
			input1Point.x[0]+=3;
			input2Point.x[0]+=3;
		}
		O.write(DomainPoint::from_point<1>(OPoint),alphaC); // Write output alpha
		OPoint.x[0]++;		// Increment for next pixel
		input1Point.x[0]++;
		input2Point.x[0]++;
	}*/
}
/*
void compositeOver(RegionAccessor<AccessorType::Generic, float> input1, RegionAccessor<AccessorType::Generic, float> input2, RegionAccessor<AccessorType::Generic, float> imgO, int start, int stop, compositeArguments co){
	/**
	 *  Alpha 'Over' Compositing
	 */
	//
	//cout<<"compositeOver called!";
//	composite(input1, input2,imgO,start,stop,co,&pass1,&exclusion);
//}

void compositeOver(RegionAccessor<AccessorType::Generic, float> input1, RegionAccessor<AccessorType::Generic, float> input2, RegionAccessor<AccessorType::Generic, float> imgO, int start, int stop, tripleArgument ta){
	/**
	 *  Alpha 'Over' Compositing
	 */
	composite(input1, input2,imgO,start,stop,ta,&pass1,&exclusion);
}

void combine_task(const Task *task,
		const std::vector<PhysicalRegion> &regions,
		Context ctx, HighLevelRuntime *runtime){
	/**
	 * Combining task that actually composites images together
	 */
	//cout<<"\n Inside combine_task";
	assert(regions.size()==3);
	//compositeArguments co = *((compositeArguments*)task->args); // Get metadata properties
		tripleArgument ta = *((tripleArgument*)task->args); // Get metadata properties
	
	Domain outDomain = runtime->get_index_space_domain(ctx,regions[2].get_logical_region().get_index_space());
	Rect<1> outRect = outDomain.get_rect<1>();			// Get the size of the return image
	PhysicalRegion img0 = regions[0];
	img0.wait_until_valid();
	PhysicalRegion img1 = regions[1];
	img1.wait_until_valid();
	PhysicalRegion img2 = regions[2];
	img2.wait_until_valid();
	RegionAccessor<AccessorType::Generic,float> inputAccessor1 = regions[0].get_field_accessor(FID_VAL).typeify<float>();
	RegionAccessor<AccessorType::Generic,float> inputAccessor2 = regions[1].get_field_accessor(FID_VAL).typeify<float>();
	RegionAccessor<AccessorType::Generic,float> outputAccessor = regions[2].get_field_accessor(FID_VAL).typeify<float>();
	compositeOver(inputAccessor1,inputAccessor2,outputAccessor,outRect.lo.x[0],outRect.hi.x[0],ta); // Call the Composite 'Over' version
}



void composite_task(const Task *task,
		const std::vector<PhysicalRegion> &regions,
		Context ctx, HighLevelRuntime *runtime){
	/**
	 * Main recursive 'compositing' task.
	 * This is an inner task with control structures only that don't access data.
	 */
	compositeArguments co = *((compositeArguments*)task->args);	// Task metadata
	int size = (co.width)*(co.maxy-co.miny+1)*4;				// Total image pixel count
	int inputRegionCount = regions.size();
	PhysicalRegion metadataPhysicalRegion = regions[0];
	LogicalRegion metadataLogicalRegion = metadataPhysicalRegion.get_logical_region();
	IndexSpace metadataIndexSpace = metadataLogicalRegion.get_index_space();
	Domain totalDomain = runtime->get_index_space_domain(ctx,metadataIndexSpace);
	Rect<1> totalRect = totalDomain.get_rect<1>();				// Get the size of the metadata region
	cout << "Compositing |";
	for(int i = 0; i < inputRegionCount; ++i){
		if(i <= totalRect.hi.x[0] && i >= totalRect.lo.x[0]) cout << "=";
		else cout << " ";
	}
	cout << endl;

	if(totalRect.lo.x[0]==totalRect.hi.x[0]){					// If the metadata region only has 1 element in it
		int i = totalRect.lo.x[0]+2;
		IndexSpace inputIndex = regions[i].get_logical_region().get_index_space();
		Domain inputDomain = runtime->get_index_space_domain(ctx,inputIndex);
		Rect<1> inputRect = inputDomain.get_rect<1>();
		RegionAccessor<AccessorType::Generic,float> inputAccessor = regions[i].get_field_accessor(FID_VAL).typeify<float>();
		RegionAccessor<AccessorType::Generic,float> imgAccessor = regions[1].get_field_accessor(FID_VAL).typeify<float>();
//		RegionAccessor<AccessorType::Generic,Image> metadataAccessor = regions[0].get_field_accessor(FID_META).typeify<Image>();
//		cout << "Compositing: " << metadataAccessor.read(DomainPoint::from_point<1>(totalRect.lo.x[0])).order << endl;
		for(GenericPointInRectIterator<1>pir(inputRect); pir; pir++){
			imgAccessor.write(DomainPoint::from_point<1>(pir.p),inputAccessor.read(DomainPoint::from_point<1>(pir.p)));
		}
	}
	else{
		FieldSpace metadataOutputField = runtime->create_field_space(ctx); // Setup field spaces for the metadata values
		{
			FieldAllocator allocator = runtime->create_field_allocator(ctx,metadataOutputField);
			allocator.allocate_field(sizeof(Image),FID_META);
		}
		// Create index spaces for splitting the metadata into two (roughly equivalent) halves.
		IndexSpace metadataOutputIndex1 = runtime->create_index_space(ctx, Domain::from_rect<1>(Rect<1>(Point<1>(totalRect.lo.x[0]),Point<1>(static_cast<int>((totalRect.hi.x[0]-totalRect.lo.x[0])/2)+totalRect.lo.x[0]))));
		LogicalRegion metadataOutputLogicalRegion1 = runtime->create_logical_region(ctx,metadataOutputIndex1,metadataOutputField);
		// And bind them to new logical regions
		IndexSpace metadataOutputIndex2 = runtime->create_index_space(ctx, Domain::from_rect<1>(Rect<1>(Point<1>(static_cast<int>((totalRect.hi.x[0]-totalRect.lo.x[0])/2)+totalRect.lo.x[0]+1),Point<1>(totalRect.hi.x[0]))));
		LogicalRegion metadataOutputLogicalRegion2 = runtime->create_logical_region(ctx,metadataOutputIndex2,metadataOutputField);

		compositeArguments co1 = co;	// Create sub-arguments for the two new tasks
		compositeArguments co2 = co;

		// Map and prepare to copy the metadata values
		RegionAccessor<AccessorType::Generic, Image> filenameAccessor = regions[0].get_field_accessor(FID_META).typeify<Image>();

		for(int r = 0; r < 2; ++r){		// Setup for the two new tasks
			RegionRequirement req;
			Domain metadataDomain;
			if(r==0){					// Bind into the first logical region
				req = RegionRequirement(metadataOutputLogicalRegion1,WRITE_DISCARD,EXCLUSIVE,metadataOutputLogicalRegion1);
				metadataDomain = runtime->get_index_space_domain(ctx,metadataOutputIndex1);
			}
			else{						// Or the second region
				req = RegionRequirement(metadataOutputLogicalRegion2,WRITE_DISCARD,EXCLUSIVE,metadataOutputLogicalRegion2);
				metadataDomain = runtime->get_index_space_domain(ctx,metadataOutputIndex2);
			}
			req.add_field(FID_META);	// Need to copy metadata values
			InlineLauncher metadataLauncher(req);
			PhysicalRegion metadataPhysicalRegion = runtime->map_region(ctx,metadataLauncher);
			metadataPhysicalRegion.wait_until_valid();
			RegionAccessor<AccessorType::Generic, Image> accessFilename = metadataPhysicalRegion.get_field_accessor(FID_META).typeify<Image>();
			Rect<1> metadataBound = metadataDomain.get_rect<1>();	// The range of indices of the metadata region

			/// Update the relative image bounds based on data bounds (does nothing for now)
			int miny = co.height;	// Set to opposite extrema
			int maxy = 0;
			for(GenericPointInRectIterator<1> pir(metadataBound); pir; pir++){ // Iterate through images and expand bounds as necessary
				Image tmpimg = filenameAccessor.read(DomainPoint::from_point<1>(pir.p));
				miny = min(tmpimg.ymin,miny);
				maxy = max(tmpimg.ymax,maxy);
				accessFilename.write(DomainPoint::from_point<1>(pir.p),tmpimg);
			}
			if(r==0){	// Define the corresponding sub-arguments
				co1.miny = miny;
				co1.maxy = maxy;
			}
			else{
				co2.miny = miny;
				co2.maxy = maxy;
			}
			runtime->unmap_region(ctx,metadataPhysicalRegion);	// Free up resources
		}


		FieldSpace inputField = runtime->create_field_space(ctx);	// Prepare the image field spaces
		{
			FieldAllocator allocator = runtime->create_field_allocator(ctx,inputField);
			allocator.allocate_field(sizeof(float),FID_VAL);
		}

		Rect<1> inputBound1(Point<1>(0),Point<1>((co.width)*(co1.maxy-co1.miny+1)*4-1)); // Create a new region for the return images
		IndexSpace inputIndex1 = runtime->create_index_space(ctx, Domain::from_rect<1>(inputBound1));

		Rect<1> inputBound2(Point<1>(0),Point<1>((co.width)*(co2.maxy-co2.miny+1)*4-1)); // One for each sub-task
		IndexSpace inputIndex2 = runtime->create_index_space(ctx, Domain::from_rect<1>(inputBound2));

		LogicalRegion imgLogicalRegion1 = runtime->create_logical_region(ctx,inputIndex1,inputField);
		LogicalRegion imgLogicalRegion2 = runtime->create_logical_region(ctx,inputIndex2,inputField);




		TaskLauncher compositeLauncher1(COMPOSITE_TASK_ID, TaskArgument(&co1,sizeof(co1))); // Launch a single task for each half of the remaining data
		compositeLauncher1.add_region_requirement(RegionRequirement(metadataOutputLogicalRegion1,READ_ONLY,EXCLUSIVE,metadataOutputLogicalRegion1));
		compositeLauncher1.add_field(0,FID_META); // Recursively calling composite task, so same arguments
		compositeLauncher1.add_region_requirement(RegionRequirement(imgLogicalRegion1,WRITE_DISCARD,EXCLUSIVE,imgLogicalRegion1));
		compositeLauncher1.add_field(1,FID_VAL);
		for(int i = 2; i < inputRegionCount; ++i){
			compositeLauncher1.add_region_requirement(RegionRequirement(regions[i].get_logical_region(),READ_ONLY,EXCLUSIVE,regions[i].get_logical_region()));
			compositeLauncher1.add_field(i,FID_VAL);
		}
		runtime->execute_task(ctx,compositeLauncher1);

		TaskLauncher compositeLauncher2(COMPOSITE_TASK_ID, TaskArgument(&co2,sizeof(co2))); // Second half of metadata
		compositeLauncher2.add_region_requirement(RegionRequirement(metadataOutputLogicalRegion2,READ_ONLY,EXCLUSIVE,metadataOutputLogicalRegion2));
		compositeLauncher2.add_field(0,FID_META);
		compositeLauncher2.add_region_requirement(RegionRequirement(imgLogicalRegion2,WRITE_DISCARD,EXCLUSIVE,imgLogicalRegion2));
		compositeLauncher2.add_field(1,FID_VAL);
		for(int i = 2; i < inputRegionCount; ++i){
			compositeLauncher2.add_region_requirement(RegionRequirement(regions[i].get_logical_region(),READ_ONLY,EXCLUSIVE,regions[i].get_logical_region()));
			compositeLauncher2.add_field(i,FID_VAL);
		}
		runtime->execute_task(ctx,compositeLauncher2);


		/// Actual Compositing

		tripleArgument ta = {co,co1,co2}; // New argument type that combines all possible metadata for the compositing operation

		LogicalRegion combineOutputRegion = regions[1].get_logical_region(); // Image return region

		const int divisions = 8; // This is an arbitrarily chosen number based on Darwin configuration (Nodes with 32 cores)

		Rect<1> combineColorRect(Point<1>(0),Point<1>(divisions-1)); // Define the number of coloring divisions
		Domain combineColorDomain = Domain::from_rect<1>(combineColorRect);

		DomainColoring combineInputColoring1; // Need three separate colorings for each of the regions
		DomainColoring combineInputColoring2; // 	Since the regions don't necessarily completely overlap
		DomainColoring combineOutputColoring;
		{
			int index = 0;	// Coloring current pixel value
			int linewidth = 4 * co.width; // Size of an image scan line
			int partitionCount = (int)(((int)(size/4)/divisions)*4/linewidth)*linewidth; // Number of pixels per partition

			int voffset1 = (co1.miny - co.miny); // Vertical offset of the first region
			int maxbound1 = linewidth * (co1.maxy - co1.miny+1); // Vertical size of the first region

			int voffset2 = (co2.miny - co.miny); // Vertical offset of the second region
			int maxbound2 = linewidth * (co2.maxy - co2.miny+1); // Vertical size of the second region

			for(int i = 0; i < divisions-1; ++i){ // Step through each color
				Rect<1> subrect(Point<1>(index),Point<1>(index+partitionCount-1)); 	// Find the size of the output region in coloring
				combineOutputColoring[i] = Domain::from_rect<1>(subrect); 			// Get a domain of this area
				int lower1 = subCoordTransform(index,linewidth,voffset1);			// Get the lower bound index for the first input region
				int upper1 = subCoordTransform(index+partitionCount-1,linewidth,voffset1); // And the upper bound index
				if(upper1 < 0 || lower1 >= maxbound1){								// Check if the bounds fit
					Rect<1> subrect1(Point<1>(0),Point<1>(0));						// If not, the coloring doesn't overlap the region
					combineInputColoring1[i] = Domain::from_rect<1>(subrect1);		// Define a hack to indicate that to the composite task
				}
				else{																// Otherwise define the mapping
					lower1 = max(lower1,0);											// Clamp it to the bounds of the region
					upper1 = min(upper1,maxbound1-1);
					Rect<1> subrect1(Point<1>(lower1),Point<1>(upper1-0));			// Find the actual rectangle of the first region
					combineInputColoring1[i] = Domain::from_rect<1>(subrect1);		// And define it in the coloring domain
				}
				int lower2 = subCoordTransform(index,linewidth,voffset2);			// Complete the same process for the second input region
				int upper2 = subCoordTransform(index+partitionCount-1,linewidth,voffset2);
				if(upper2 < 0 || lower2 >= maxbound2){
					Rect<1> subrect2(Point<1>(0),Point<1>(0));
					combineInputColoring2[i] = Domain::from_rect<1>(subrect2);
				}
				else{
					lower2 = max(lower2,0);
					upper2 = min(upper2,maxbound2-1);
					Rect<1> subrect2(Point<1>(lower2),Point<1>(upper2-0));
					combineInputColoring2[i] = Domain::from_rect<1>(subrect2);
				}
				index += partitionCount;											// Keep track of the current index
			}
			Rect<1> subrect(Point<1>(index),Point<1>(size-1));						// In case we have any extra space, fit it all into the last color
			combineOutputColoring[divisions-1] = Domain::from_rect<1>(subrect);

			int lower1 = subCoordTransform(index,linewidth,voffset1);				// Same process, just with everything until the end
			int upper1 = subCoordTransform(size-1,linewidth,voffset1);
			if(upper1 < 0 || lower1 >= maxbound1){
				Rect<1> subrect1(Point<1>(0),Point<1>(0));
				combineInputColoring1[divisions-1] = Domain::from_rect<1>(subrect1);
			}
			else{
				lower1 = max(lower1,0);
				upper1 = min(upper1,maxbound1-1);
				Rect<1> subrect1(Point<1>(lower1),Point<1>(upper1-0));
				combineInputColoring1[divisions-1] = Domain::from_rect<1>(subrect1);
			}
				int lower2 = subCoordTransform(index,linewidth,voffset2);
			int upper2 = subCoordTransform(size-1,linewidth,voffset2);
			if(upper2 < 0 || lower2 >= maxbound2){
				Rect<1> subrect2(Point<1>(0),Point<1>(0));
				combineInputColoring2[divisions-1] = Domain::from_rect<1>(subrect2);
			}
			else{
				lower2 = max(lower2,0);
				upper2 = min(upper2,maxbound2-1);
				Rect<1> subrect2(Point<1>(lower2),Point<1>(upper2+0));
				combineInputColoring2[divisions-1] = Domain::from_rect<1>(subrect2);
			}
		}
		


		/// Pass the individual colorings into Index Partitions
		IndexPartition combineInputIndexPartition1 = runtime->create_index_partition(ctx, inputIndex1, combineColorDomain, combineInputColoring1, false);
		IndexPartition combineInputIndexPartition2 = runtime->create_index_partition(ctx, inputIndex2, combineColorDomain, combineInputColoring2, false);
		IndexPartition combineOutputIndexPartition = runtime->create_index_partition(ctx, combineOutputRegion.get_index_space(), combineColorDomain, combineOutputColoring, true);

		/// And define logical partitions on those index partitions
		LogicalPartition combineInputLogicalPartition1 = runtime->get_logical_partition(ctx, imgLogicalRegion1, combineInputIndexPartition1);
		LogicalPartition combineInputLogicalPartition2 = runtime->get_logical_partition(ctx, imgLogicalRegion2, combineInputIndexPartition2);
		LogicalPartition combineOutputLogicalPartition = runtime->get_logical_partition(ctx, combineOutputRegion, combineOutputIndexPartition);


		ArgumentMap argMap;

		/// Index launcher for the actual combinator tasks
		IndexLauncher combineLauncher(COMBINE_TASK_ID, combineColorDomain,TaskArgument(&ta,sizeof(ta)), argMap);
		combineLauncher.add_region_requirement(RegionRequirement(combineInputLogicalPartition1,0,READ_ONLY,EXCLUSIVE,imgLogicalRegion1));
		combineLauncher.add_field(0,FID_VAL);	// First region is the first input region
		combineLauncher.add_region_requirement(RegionRequirement(combineInputLogicalPartition2,0,READ_ONLY,EXCLUSIVE,imgLogicalRegion2));
		combineLauncher.add_field(1,FID_VAL);	// Second region is the other input region
		combineLauncher.add_region_requirement(RegionRequirement(combineOutputLogicalPartition,0,WRITE_DISCARD,EXCLUSIVE,combineOutputRegion));
		combineLauncher.add_field(2,FID_VAL);	// Third region is the output image region
		runtime->execute_index_space(ctx,combineLauncher);

		/// Clean up created components afterwards
		runtime->destroy_logical_region(ctx,metadataOutputLogicalRegion2);
		runtime->destroy_logical_region(ctx,metadataOutputLogicalRegion1);
		runtime->destroy_field_space(ctx,metadataOutputField);
		runtime->destroy_index_space(ctx,metadataOutputIndex2);
		runtime->destroy_index_space(ctx,metadataOutputIndex1);
		runtime->destroy_logical_region(ctx,imgLogicalRegion2);
		runtime->destroy_logical_region(ctx,imgLogicalRegion1);
		runtime->destroy_index_space(ctx,inputIndex1);
		runtime->destroy_index_space(ctx,inputIndex2);
		runtime->destroy_field_space(ctx,inputField);
		
		

	}

	cout << "Finishing   |";
	for(int i = 0; i < inputRegionCount; ++i){
		if(i <= totalRect.hi.x[0] && i >= totalRect.lo.x[0]) cout << "=";
		else cout << " ";
	}
	cout << endl;
}

void display_task(const Task *task,
		const std::vector<PhysicalRegion> &regions,
		Context ctx, HighLevelRuntime *runtime){
	/**
	 * Task for sending data to Qt Window
	 */
	Image img = *((Image*)task->args); // Load the image metadata
	PhysicalRegion imgPhysicalRegion = regions[0];				// The region that holds the image pixels
	LogicalRegion imgLogicalRegion = imgPhysicalRegion.get_logical_region();
	IndexSpace imgIndexSpace = imgLogicalRegion.get_index_space();
	Domain imgDomain = runtime->get_index_space_domain(ctx,imgIndexSpace);
	Rect<1> imgBound = imgDomain.get_rect<1>();					// Get the size of the pixel data
	RegionAccessor<AccessorType::Generic,float> accessImg = imgPhysicalRegion.get_field_accessor(FID_VAL).typeify<float>();
	imgPhysicalRegion.wait_until_valid();

	//modified by : Vidhi on July 17, 2016
        std::vector< unsigned char> r;  // red
	std::vector< unsigned char> g;  // green
        std::vector< unsigned char> b;  // blue
	std::vector< unsigned char> a;  // alpha
	cout << "Writing to File out.png" << endl;
	int counter = 0;
	for(GenericPointInRectIterator<1> pir(imgBound); pir; pir++){
		float val = accessImg.read(DomainPoint::from_point<1>(pir.p));
	if(counter == 0)
			r.push_back(val*255);//reinterpret_cast<char*>(&val));
		else if(counter == 1)
			g.push_back(val*255);//reinterpret_cast<char*>(&val));
		else if(counter == 2)
			b.push_back(val*255);//reinterpret_cast<char*>(&val));
		else if(counter == 3)
			a.push_back(val*255);//reinterpret_cast<char*>(&val));
		counter++;
		if(counter>3) counter = 0;
	}
	boost::gil::rgba8c_planar_view_t view = boost::gil::planar_rgba_view(img.width, img.height, r.data(), g.data(), b.data(), a.data(), img.width);
	boost::gil::png_write_view("out.png", view);

}

Future setupCombine(Context ctx, HighLevelRuntime *runtime, LogicalRegion input1, LogicalRegion input2, LogicalRegion output, compositeArguments co){
	TaskLauncher combineLauncher(COMBINE_TASK_ID, TaskArgument(&co,sizeof(co)));
	combineLauncher.add_region_requirement(RegionRequirement(input1,READ_ONLY,EXCLUSIVE,input1));
	combineLauncher.add_field(0,FID_VAL);
	combineLauncher.add_region_requirement(RegionRequirement(input2,READ_ONLY,EXCLUSIVE,input2));
	combineLauncher.add_field(1,FID_VAL);
	combineLauncher.add_region_requirement(RegionRequirement(output,WRITE_DISCARD,EXCLUSIVE,output));
	combineLauncher.add_field(2,FID_VAL);
	return runtime->execute_task(ctx,combineLauncher);
}


vector<LogicalRegion> loadRenderCPU(Context ctx, HighLevelRuntime *runtime, int width, int height, Movement mov, IndexSpace imgIndex, FieldSpace imgField){
	vector<LogicalRegion> imgs;
	vector<Future> futures;
	int num_subregions = runtime->get_tunable_value(ctx, SUBREGION_TUNABLE, PARTITIONING_MAPPER_ID);
	int i = 0;
	for(int y = 0; y < height; y+=20){
		for(int x = 0; x < width; x+=20){
			Image img;
			img.width = width;
			img.height = height;
			img.partition = (DataPartition){0,0,0,(float)x,(float)(x+10),(float)y,(float)(y+10),0,0};
			img.order = i++;
			img.core = i % num_subregions;
			img.randomseed = rand();
			LogicalRegion imgLogicalRegion = runtime->create_logical_region(ctx,imgIndex,imgField);
			TaskLauncher loadLauncher(CPU_DRAW_TASK_ID, TaskArgument(&img,sizeof(img)));	// Spawn the renderer task
			loadLauncher.add_region_requirement(RegionRequirement(imgLogicalRegion,READ_WRITE,EXCLUSIVE,imgLogicalRegion));
			loadLauncher.add_field(0,FID_VAL);
			futures.push_back(runtime->execute_task(ctx,loadLauncher));	// Launch and terminate render task

			imgs.push_back(imgLogicalRegion);

		}
	}
	for(unsigned int i = 0; i < futures.size(); ++i){
		futures[i].get_void_result();
	}
	return imgs;
}

vector<LogicalRegion> loadRenderHeptane(Context ctx, HighLevelRuntime *runtime, int width, int height, Movement mov, IndexSpace imgIndex, FieldSpace imgField){
	vector<LogicalRegion> imgs;
	vector<Future> futures;
	int num_subregions = runtime->get_tunable_value(ctx, SUBREGION_TUNABLE, PARTITIONING_MAPPER_ID);
	int t = 0;
	for(int i = 0; i < 4; i++){
		for(int j = 0; j < 4; j++){
			for(int k = 0; k < 4; k++){
				int x = i * 75;
				int y = j * 75;
				int z = k * 75;
				Image img;
				img.width = width;
				img.height = height;
				img.partition = (DataPartition){75,75,75,(float)x,(float)(x+74),(float)y,(float)(y+74),(float)z,(float)(z+74)};
				if(i == 3){ img.partition.xmax = x + 76; img.partition.datx = 77; }
				if(j == 3){ img.partition.ymax = y + 76; img.partition.daty = 77; }
				if(k == 3){ img.partition.zmax = z + 76; img.partition.datz = 77; }
				for(int l = 0; l < 16; ++l)
					img.invPVM[l] = mov.invPVM[l];
				img.i = i;
				img.j = j;
				img.k = k;
				img.order = t++;
				img.core = t % num_subregions;
				img.randomseed = rand();
				LogicalRegion imgLogicalRegion = runtime->create_logical_region(ctx,imgIndex,imgField);
				TaskLauncher loadLauncher(CREATE_INTERFACE_TASK_ID, TaskArgument(&img,sizeof(img)));	// Spawn the renderer task
				loadLauncher.add_region_requirement(RegionRequirement(imgLogicalRegion,READ_WRITE,EXCLUSIVE,imgLogicalRegion));
				loadLauncher.add_field(0,FID_VAL);		// Output Image as second region
				futures.push_back(runtime->execute_task(ctx,loadLauncher));	// Launch and terminate render task
							imgs.push_back(imgLogicalRegion);
				if(t >= num_subregions){
					for(unsigned int i = 0; i < futures.size(); ++i){
						futures[i].get_void_result();
					}
					return imgs;
				}
			}
		}
	}
	for(unsigned int i = 0; i < futures.size(); ++i){
		futures[i].get_void_result();
	}
	return imgs;
}

vector<LogicalRegion> loadRenderHeptane2(Context ctx, HighLevelRuntime *runtime, int width, int height, Movement mov, IndexSpace imgIndex, FieldSpace imgField){
	vector<LogicalRegion> imgs;
	vector<Future> futures;

	Image img;
	img.width = width;
	img.height = height;
	img.partition = (DataPartition){302,302,302,0,301,0,301,0,301};
	for(int j = 0; j < 16; ++j)
		img.invPVM[j] = mov.invPVM[j];
	img.order = 1;
	img.core = 1;
	img.randomseed = rand();
	LogicalRegion imgLogicalRegion = runtime->create_logical_region(ctx,imgIndex,imgField);
	TaskLauncher loadLauncher(CREATE_INTERFACE_TASK_ID, TaskArgument(&img,sizeof(img)));	// Spawn the renderer task
	loadLauncher.add_region_requirement(RegionRequirement(imgLogicalRegion,READ_WRITE,EXCLUSIVE,imgLogicalRegion));
	loadLauncher.add_field(0,FID_VAL);		// Output Image as second region
	futures.push_back(runtime->execute_task(ctx,loadLauncher));	// Launch and terminate render task

	imgs.push_back(imgLogicalRegion);


	for(unsigned int i = 0; i < futures.size(); ++i){
		futures[i].get_void_result();
	}
	return imgs;
}

void top_level_task(const Task *task,
const std::vector<PhysicalRegion> &regions,
		Context ctx, HighLevelRuntime *runtime){

	cout << "Reading data from file..." << endl;

	srand(time(NULL));
	int width = 1000;
	int height = 1000;
	Movement mov = {{141.421, 100., 2000., 2166.42, 0., 141.421, -2828.43, -2651.65, \
			141.421, -100., -2000., -1883.58, 0., 0., 0., 1.},1.0};


	const InputArgs &command_args = HighLevelRuntime::get_input_args();
	if (command_args.argc > 1){
		width = atoi(command_args.argv[1]);
		height = width;
		assert(width >= 0);
	}

	Rect<1> imgBound(Point<1>(0),Point<1>(width*height*4-1));
	IndexSpace imgIndex = runtime->create_index_space(ctx, Domain::from_rect<1>(imgBound));
	FieldSpace imgField = runtime->create_field_space(ctx);
	{
		FieldAllocator allocator = runtime->create_field_allocator(ctx,imgField);
		allocator.allocate_field(sizeof(float),FID_VAL);
	}
	vector<LogicalRegion> imgLogicalRegions = loadRenderHeptane(ctx, runtime, width, height, mov, imgIndex, imgField);

	//int numFiles = 20;							// Choose to create two partitions of the data
	vector<Image> images;						// Array to hold the metadata values in
	int xindex = 0;					
	
	cout<<"\n size of imgLogicalRegions = "<<imgLogicalRegions.size();

	vector <LogicalRegion>  imgLogicalRegions_new;

	for(int i = 0; i < imgLogicalRegions.size(); i++){
	
		PhysicalRegion imgPhysicalRegion = regions[i];
		LogicalRegion imgLogicalRegion = imgLogicalRegions[i];
		IndexSpace imgIndexSpace = imgLogicalRegion.get_index_space();
		Domain imgDomain = runtime->get_index_space_domain(ctx,imgIndexSpace);
		Rect<1> imgBound = imgDomain.get_rect<1>();                                     // Get the size of the pixel data
		RegionAccessor<AccessorType::Generic,float> accessImg = imgPhysicalRegion.get_field_accessor(FID_VAL).typeify<float>();
		imgPhysicalRegion.wait_until_valid();
		std::vector< unsigned char> r;  // red
		std::vector< unsigned char> g;  // green
		std::vector< unsigned char> b;  // blue
		std::vector< unsigned char> a;  // alpha	
	
		//calculate bounding box using a naive approach of traversing all pixels:
		int counter = 0, first_rx = -1,first_ry = -1, first_gx = -1, first_gy = -1, first_bx = -1, first_by = -1,  first_ax = -1, first_ay = -1, last_rx = -1, last_ry = -1, last_gx = -1, last_gy = -1,  last_bx = -1, last_by = -1, last_ax = -1, last_ay = -1;
		int original_counter = 0;
		int box_topleftx, box_toplefty, box_botrightx, box_botrighty;
		for(GenericPointInRectIterator<1> pir(imgBound); pir; pir++){

			float val = accessImg.read(DomainPoint::from_point<1>(pir.p));
			if(counter == 0){
				r.push_back(val*255);
				if((first_rx==-1||first_rx > original_counter%500) && val>0)
				{	first_rx = original_counter%500;}
				if((first_ry==-1 || original_counter/500 <  first_ry) && val > 0)
				{	first_ry = original_counter/500;}
				
				if((last_rx < original_counter%500 || last_rx == -1) &&val>0)
				{	last_rx = original_counter%500;}
				if((last_ry <original_counter/500|| last_ry == -1) && val>0)
				{ 	last_ry = original_counter/500;}
				
			}

			else if(counter == 1){
				g.push_back(val*255);
				if((first_gx==-1||first_gx > original_counter%500) && val>0)
				{       first_gx = original_counter%500;}
				if((first_gy==-1 || original_counter/500 <  first_gy) && val > 0)
				{       first_gy = original_counter/500;}
				
				if((last_gx < original_counter%500 || last_gx == -1) &&val>0)
				{       last_gx = original_counter%500;}
				if((last_gy <original_counter/500|| last_gy == -1) && val>0)
				{       last_gy = original_counter/500;}

			}
			else if(counter == 2){
				b.push_back(val*255);
				if((first_bx==-1||first_bx > original_counter%500) && val>0)
				{       first_bx = original_counter%500;}
				if((first_by==-1 || original_counter/500 <  first_by) && val > 0)
				{       first_by = original_counter/500;}
				
				if((last_bx < original_counter%500 || last_bx == -1) &&val>0)
				{       last_bx = original_counter%500;}
				if((last_by <original_counter/500|| last_by == -1) && val>0)
				{       last_by = original_counter/500;}
			}
			else if(counter == 3){
				a.push_back(val*255);
				
				if((first_ax==-1||first_ax > original_counter%500) && val>0)
				{       first_ax = original_counter%500;}
				if((first_ay==-1 || original_counter/500 <  first_ay) && val > 0)
				{       first_ay = original_counter/500;}
				
				if((last_ax < original_counter%500 || last_ax == -1) &&val>0)
				{       last_ax = original_counter%500;}
				if((last_ay <original_counter/500|| last_ay == -1) && val>0)
				{       last_ay = original_counter/500;}

			}
			counter++;			
			if(counter>3) {
				counter = 0; original_counter++;
			}
		}
		
		//write the individual rendered images (without cropping using the bounding boxes found):
		//cout << "Writing to File out"<<i<<".png" << endl;
		
		//boost::gil::rgba8c_planar_view_t view = boost::gil::planar_rgba_view(width, height, r.data(), g.data(), b.data(), a.data(), width);
		//std::string filename ("out");
		//filename = filename+""+std::to_string(i)+".png";
		//boost::gil::png_write_view(filename, view);
		
		int all_x[] = {first_ax,last_ax,first_rx, last_rx, first_gx, last_gx, first_bx, last_bx}; 
		int all_y[] = {first_ay,last_ay,first_ry, last_ry, first_gy, last_gy, first_by, last_by};

		//coordinates of bounding box:
		box_topleftx = *std::min_element(all_x,all_x+sizeof(all_x)/sizeof(all_x[0]));
		box_toplefty = *std::min_element(all_y,all_y+sizeof(all_y)/sizeof(all_y[0]));
		box_botrightx = *std::max_element(all_x,all_x+sizeof(all_x)/sizeof(all_x[0]));
		box_botrighty = *std::max_element(all_y,all_y+sizeof(all_y)/sizeof(all_y[0]));
		//cout<<"\n Bounding BOX : {("<<*box_topleftx<<","<<*box_toplefty<<") to ("<<*box_botrightx<<","<<*box_botrighty<<")}\n";
		
		//crop the rendered images using the bounding box found : modify the imgLogicalRegion holding the images
		std::vector <unsigned char> rnew;  // red
		std::vector<unsigned char> gnew;  // green
		std::vector<unsigned char> bnew;  // blue
		std::vector< unsigned char> anew;  // alpha       
		
		//new height and width
		int wnew = box_botrightx - box_topleftx +1;
		int hnew = box_botrighty - box_toplefty +1;

		for(int i = 0; i<a.size();i++){
			//check if it lies within the bounding box or not:
			if(i%500 >= box_topleftx && i/500 >= box_toplefty && i%500 <= box_botrightx && i/500 <= box_botrighty)
			{	
				rnew.push_back( r[i]);
				gnew.push_back( g[i]);
				bnew.push_back( b[i]);
				anew.push_back( a[i] );
			}
		
		}

		//write cropped file to verify if the bounding boxes are correctly calculated
		cout<<"\n writing cropped file:"<<endl;
		std::string filename2 ("out");
		filename2 = filename2+"_cropped"+std::to_string(i)+".png";
		boost::gil::rgba8c_planar_view_t view2 = boost::gil::planar_rgba_view(wnew, hnew, rnew.data(), gnew.data(), bnew.data(), anew.data(), wnew);
		boost::gil::png_write_view(filename2, view2);

		//Ian's old (master branch) code copied from here::		

		Image newimg;							// Create a metadata object to hold values
		newimg.width = wnew;					// This data gets sent to the renderer (necessary)
		newimg.height = hnew;					// 		This is total image Width and Height
		for(int j = 0; j < 16; ++j)
			newimg.invPVM[j] = mov.invPVM[j];	// Copy the transformation matrix over
		newimg.xmin = box_topleftx;						// Values for the extent of the render within the image
		newimg.xmax = box_botrightx;					// 		Set to be the entire size for now
		newimg.ymin = box_toplefty;						//		Need to feed partition bounds into the modelview to get these
		newimg.ymax = box_botrighty;
		newimg.partition = (DataPartition)(DataPartition){75,75,75,(float)i,(float)(i+74),(float)i,(float)(i+74),(float)i,(float)(i+74)};//{xindex,i==numFiles-1 ? (int)datx : xindex+xspan+10,0,(int)daty, 0,(int)datz}; // Define the data partitioning
		newimg.order = mov.xdat * (float)xindex;// Feed the partition value into the modelview to get the compositing order
		images.push_back(newimg);				// Add the metadata to the array
		xindex += 75;						// Iterate index values
		imgLogicalRegions_new.push_back(runtime->create_logical_region(ctx,imgIndex,imgField));
		runtime->destroy_logical_region(ctx,imgLogicalRegions[i]);

	}
			
	sort(images.rbegin(),images.rend());		// Sort the metadata in reverse value of order

	cout << "Spawning with <" << images.size() << "> Images" << endl;

	Rect<1> metadataBound(Point<1>(0),Point<1>( imgLogicalRegions_new.size()-1));	// Set up an index space for the metadata
	IndexSpace taskIndex = runtime->create_index_space(ctx, Domain::from_rect<1>(metadataBound));
	FieldSpace metadataField = runtime->create_field_space(ctx);
	{
		FieldAllocator allocator = runtime->create_field_allocator(ctx,metadataField);
		allocator.allocate_field(sizeof(Image),FID_META);	// Classify it with the META field value
	}

	LogicalRegion metadataLogicalRegion = runtime->create_logical_region(ctx, taskIndex, metadataField);


	{	
		// Fill the metadata Logical Region with the previously generated metadata
		RegionRequirement req(metadataLogicalRegion,WRITE_DISCARD,EXCLUSIVE,metadataLogicalRegion);
		req.add_field(FID_META);
		InlineLauncher metadataLauncher(req);

		PhysicalRegion metadataPhysicalRegion = runtime->map_region(ctx,metadataLauncher);
		metadataPhysicalRegion.wait_until_valid();

		RegionAccessor<AccessorType::Generic, Image> accessFilename = metadataPhysicalRegion.get_field_accessor(FID_META).typeify<Image>();
		// Using the 'Image' struct type
		int i = 0;
		for(GenericPointInRectIterator<1> pir(metadataBound); pir; pir++)
			accessFilename.write(DomainPoint::from_point<1>(pir.p),images[i++]);	// Make sure to write in the sorted order

		runtime->unmap_region(ctx,metadataPhysicalRegion);		// Free up resources
	}

	compositeArguments co;
	co.width = width;			// Total image size
	co.height = height;
	co.mov = mov;				// Inverse PV Matrix for tagging image
	co.miny = 0;				// Image possible extent (in Y-Dimension)
	co.maxy = height-1;			// For first level, must be entire image

	TaskLauncher compositeLauncher(COMPOSITE_TASK_ID, TaskArgument(&co,sizeof(co)));
	compositeLauncher.add_region_requirement(RegionRequirement(metadataLogicalRegion,READ_ONLY,EXCLUSIVE,metadataLogicalRegion));
	compositeLauncher.add_field(0,FID_META); 	// Metadata as first region
	LogicalRegion out1 = runtime->create_logical_region(ctx,imgIndex,imgField);
	compositeLauncher.add_region_requirement(RegionRequirement(out1,WRITE_DISCARD,EXCLUSIVE,out1));
	compositeLauncher.add_field(1,FID_VAL);		// Output Image as second region
	for(unsigned int i = 0; i < images.size(); ++i){
		compositeLauncher.add_region_requirement(RegionRequirement(imgLogicalRegions_new[i],READ_ONLY,EXCLUSIVE,imgLogicalRegions_new[i]));
		compositeLauncher.add_field(2+i,FID_VAL);
	}
	runtime->execute_task(ctx,compositeLauncher);
	
/*
	// Build a blanced binary tree for composition
	vector<LogicalRegion> nodes;
	vector<int> cores;
	int num_subregions = runtime->get_tunable_value(ctx, SUBREGION_TUNABLE, PARTITIONING_MAPPER_ID);
	for(unsigned int i = 0; i < imgLogicalRegions.size(); ++i){
		nodes.push_back(imgLogicalRegions[i]);
		cores.push_back(i % num_subregions);
	}
	Future f;
	
	while(nodes.size()>1){
		cout << "Spawning " << nodes.size() << " nodes" << endl;
		vector<LogicalRegion> oldnodes = nodes;
		vector<int> oldcores = cores;
		nodes.clear();
		cores.clear();
		while(oldnodes.size() > 1){
			LogicalRegion output = runtime->create_logical_region(ctx,imgIndex,imgField);
			int primary = rand() % 2;
			co.core = oldcores[primary];
			f = setupCombine(ctx,runtime,oldnodes[0],oldnodes[1],output,co);
			oldnodes.erase(oldnodes.begin(),oldnodes.begin()+2);
			cores.push_back(oldcores[primary]);
			oldcores.erase(oldcores.begin(),oldcores.begin()+2);
			nodes.push_back(output);
		}
		if(oldnodes.size()==1) {
			nodes.push_back(oldnodes[0]);
			cores.push_back(oldcores[0]);
		}
	}
	cout << "Done Spawning" << endl;
	f.get_void_result();*/
	cout << "Done Compositing" << endl;
	
	TaskLauncher displayLauncher(DISPLAY_TASK_ID, TaskArgument(&co,sizeof(co)));	// Spawn a task for sending to Qt
	displayLauncher.add_region_requirement(RegionRequirement(out1,READ_ONLY,EXCLUSIVE,out1));
	displayLauncher.add_field(0,FID_VAL);	// Only needs the image (will map once compositor is done)

	runtime->execute_task(ctx,displayLauncher); // Run the display Task
}

CompositeMapper::CompositeMapper(Machine m, HighLevelRuntime *rt, Processor p) : DefaultMapper(m, rt, p){
	/**
	 * Mapper for the compositor and renderer (will need to be modified for in-situ)
	 */
	stealing_enabled = false;
	const std::set<Processor> &cpu_procs = machine_interface.filter_processors(Processor::LOC_PROC);
	std::copy(cpu_procs.begin(), cpu_procs.end(), std::back_inserter(all_cpus));
}

void CompositeMapper::select_task_options(Task *task){
	/**
	 * Specify properties and location of tasks
	 */
	task->inline_task = false;	// All of these off
	task->spawn_task = false;
	task->map_locally = false;
	task->profile_task = false;
	task->task_priority = 0;	// Can be used to specify some execution order (TO DO)
	if(task->task_id == CREATE_TASK_ID){ // Map the GPU tasks onto the GPU, though
		std::set<Processor> connectedProcs;
		machine.get_local_processors_by_kind(connectedProcs,Processor::TOC_PROC);
		machine.get_shared_processors(task->regions[1].selected_memory,connectedProcs);
		task->target_proc = DefaultMapper::select_random_processor(connectedProcs, Processor::TOC_PROC, machine);
	}
	else if(task->task_id == CPU_DRAW_TASK_ID || task->task_id==CREATE_INTERFACE_TASK_ID){
		Image img = *((Image*)task->args);
		task->target_proc = all_cpus[img.core];
	}
	else if(task->task_id == COMBINE_TASK_ID){
		compositeArguments co = *((compositeArguments*)task->args);
		task->target_proc = all_cpus[co.core];
	}
	else{
			std::set<Processor> all_procs_2;
			machine.get_all_processors(all_procs_2);
			task->target_proc = DefaultMapper::select_random_processor(all_procs_2, Processor::LOC_PROC, machine);
	}
}


bool CompositeMapper::map_task(Task *task){
	/**
	 * Control memory mapping for each task
	 */
	if (task->task_id == CREATE_TASK_ID){ // If running on the GPU
		Memory fb_mem = machine_interface.find_memory_kind(task->target_proc,Memory::GPU_FB_MEM); // Get FrameBuffer Memories
		assert(fb_mem.exists()); // Make sure it is supported
		for (unsigned idx = 0; idx < task->regions.size(); idx++){ 	// Step through all regions
			task->regions[idx].target_ranking.push_back(fb_mem); 	//	and map them to the framebuffer memory
			task->regions[idx].virtual_map = false;
			task->regions[idx].enable_WAR_optimization = war_enabled;
			task->regions[idx].reduction_list = false;
			// Make everything SOA
			task->regions[idx].blocking_factor = task->regions[idx].max_blocking_factor;
		}
	}
	else{
		// Put everything else in the system memory
		Memory sys_mem = machine_interface.find_memory_kind(task->target_proc,Memory::SYSTEM_MEM);
		assert(sys_mem.exists());
		for (unsigned idx = 0; idx < task->regions.size(); idx++)
		{
			task->regions[idx].target_ranking.push_back(sys_mem);
			task->regions[idx].virtual_map = false;
			task->regions[idx].enable_WAR_optimization = war_enabled;
			task->regions[idx].reduction_list = false;
			// Make everything SOA
			task->regions[idx].blocking_factor = task->regions[idx].max_blocking_factor;
		}
	}
	return false;
}

PartitioningMapper::PartitioningMapper(Machine m,
                                       HighLevelRuntime *rt,
                                       Processor p)
  : DefaultMapper(m, rt, p)
{
}

int PartitioningMapper::get_tunable_value(const Task *task,
                                          TunableID tid,
                                          MappingTagID tag)
{
  if (tid == SUBREGION_TUNABLE)
  {
    const std::set<Processor> &cpu_procs = machine_interface.filter_processors(Processor::LOC_PROC);
    return cpu_procs.size();
  }
  // Should never get here
  assert(false);
  return 0;
}



void mapper_registration(Machine machine, HighLevelRuntime *rt, const std::set<Processor> &local_procs){
	/**
	 * Register this mapper for each processor
	 */
	for (std::set<Processor>::const_iterator it = local_procs.begin(); it != local_procs.end(); it++){
		rt->replace_default_mapper(new CompositeMapper(machine, rt, *it), *it); // Step through all processors and create a mapper instance
	    rt->add_mapper(PARTITIONING_MAPPER_ID, new PartitioningMapper(machine, rt, *it), *it);
	}
}


int main(int argc, char **argv){
	//cout<<"\n Inside main!";
	HighLevelRuntime::set_top_level_task_id(TOP_LEVEL_TASK_ID);
	HighLevelRuntime::register_legion_task<top_level_task>(TOP_LEVEL_TASK_ID,
				Processor::LOC_PROC, true/*single*/, false/*index*/,
				AUTO_GENERATE_ID, TaskConfigOptions(), "top_level");
	HighLevelRuntime::register_legion_task<display_task>(DISPLAY_TASK_ID,		// Register Qt Display connection task (Leaf Task)
				Processor::LOC_PROC, true/*single*/, true/*index*/,
				AUTO_GENERATE_ID, TaskConfigOptions(false, false), "display_task");
	HighLevelRuntime::register_legion_task<create_task>(CREATE_TASK_ID,			// Register the GPU render task (Leaf Task, TOC processor)
				Processor::TOC_PROC, true/*single*/, true/*index*/,
				AUTO_GENERATE_ID, TaskConfigOptions(false, false), "create_task");
	HighLevelRuntime::register_legion_task<composite_task>(COMPOSITE_TASK_ID, 	// Register the composite task (Inner Task)
			Processor::LOC_PROC, true/*single*/, true/*index*/,
			AUTO_GENERATE_ID, TaskConfigOptions(false, false), "composite_task");	
	HighLevelRuntime::register_legion_task<combine_task>(COMBINE_TASK_ID,		// Register combination task (Leaf Task)
			Processor::LOC_PROC, true/*single*/, true/*index*/,
			AUTO_GENERATE_ID, TaskConfigOptions(false, false), "combine_task");
	HighLevelRuntime::register_legion_task<create_interface_task>(CREATE_INTERFACE_TASK_ID,		// Register Qt Display connection task (Leaf Task)
				Processor::LOC_PROC, true/*single*/, true/*index*/,
				AUTO_GENERATE_ID, TaskConfigOptions(false, false), "create_interface_task");
	HighLevelRuntime::register_legion_task<cpu_draw_task>(CPU_DRAW_TASK_ID,
				Processor::LOC_PROC, true/*single*/, true/*index*/,
				AUTO_GENERATE_ID, TaskConfigOptions(false, false), "cpu_draw_task");
	HighLevelRuntime::set_registration_callback(mapper_registration);
	return HighLevelRuntime::start(argc, argv);
}
