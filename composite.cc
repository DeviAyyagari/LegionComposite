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
#include <time.h>       /* time */
#include "composite.h"
#include "DataMgr.h"

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
//	sprintf(filename,"heptane_%d_%d_%d.raw",img.i,img.j,img.k);
//	sprintf(filename,"heptane.raw");
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

void composite(RegionAccessor<AccessorType::Generic, float> input1, RegionAccessor<AccessorType::Generic, float> input2, RegionAccessor<AccessorType::Generic, float> O, int start, int stop, compositeArguments co, float (*FA)(float), float (*FB)(float)){
	/**
	 * Generic alpha compositing that can be called with adjustable functions for different operations
	 */
	Point<1> input1Point(start);
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
	}
}

void compositeOver(RegionAccessor<AccessorType::Generic, float> input1, RegionAccessor<AccessorType::Generic, float> input2, RegionAccessor<AccessorType::Generic, float> imgO, int start, int stop, compositeArguments co){
	/**
	 *  Alpha 'Over' Compositing
	 */
	composite(input1, input2,imgO,start,stop,co,&pass1,&exclusion);
}

void combine_task(const Task *task,
		const std::vector<PhysicalRegion> &regions,
		Context ctx, HighLevelRuntime *runtime){
	/**
	 * Combining task that actually composites images together
	 */
	assert(regions.size()==3);
	compositeArguments co = *((compositeArguments*)task->args); // Get metadata properties
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
	compositeOver(inputAccessor1,inputAccessor2,outputAccessor,outRect.lo.x[0],outRect.hi.x[0],co); // Call the Composite 'Over' version
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
	cout << "Writing to File" << endl;
	char filename[50];
	sprintf(filename,"output.raw");
	ofstream oFile(filename, ios::out | ios::binary);
	oFile.write(reinterpret_cast<char*>(&img.width),sizeof(int));
	oFile.write(reinterpret_cast<char*>(&img.height),sizeof(int));
	for(GenericPointInRectIterator<1> pir(imgBound); pir; pir++){
		float val = accessImg.read(DomainPoint::from_point<1>(pir.p));
		oFile.write(reinterpret_cast<char*>(&val),sizeof(float));
	}
	oFile.close();
	cout << "Finished Writing" << endl;
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

	cout << "Done rendering" << endl;


	compositeArguments co;
	co.width = width;			// Total image size
	co.height = height;
	co.mov = mov;				// Inverse PV Matrix for tagging image
	co.miny = 0;				// Image possible extent (in Y-Dimension)
	co.maxy = height-1;			// For first level, must be entire image
	

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
	f.get_void_result();
	cout << "Done Compositing" << endl;
	
	TaskLauncher displayLauncher(DISPLAY_TASK_ID, TaskArgument(&co,sizeof(co)));	// Spawn a task for sending to Qt
	displayLauncher.add_region_requirement(RegionRequirement(nodes[0],READ_ONLY,EXCLUSIVE,nodes[0]));
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
