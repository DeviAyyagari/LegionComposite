/**
 * Ian Sohl - 2015
 * Copyright (c) 2015      Los Alamos National Security, LLC
 *                         All rights reserved.
 * Legion Image Composition - Main Header
 */

#ifndef COMPOSITE_H
#define COMPOSITE_H


#include "legion.h"
#include "default_mapper.h"


using namespace LegionRuntime::HighLevel;
using namespace LegionRuntime::Accessor;


enum TaskIDs{
	TOP_LEVEL_TASK_ID,
	CREATE_TASK_ID,
	DISPLAY_TASK_ID,
	COMBINE_TASK_ID,
	CREATE_INTERFACE_TASK_ID,
	CPU_DRAW_TASK_ID,
};

enum FieldIDs{
	FID_META,
	FID_VAL,
};


struct Movement{
	float invPVM[16]; /**< Inverse PV Matrix for rendering */
	float xdat;		  /**< X[3] value for composition ordering */

	bool operator==( const Movement& rhs ) const {
		/**
		 * Manually check for invPVM equivalence
		 */
		for(int i = 0; i < 16; ++i){
			if(abs(invPVM[i]-rhs.invPVM[i])>0.000001){ // Floating point problems
				return false;
			}
		}
		return true;
	}

	Movement& operator =(const Movement& a){
		/**
		 * Manual assignment
		 */
		for(int i = 0; i < 16; ++i){
			invPVM[i] = a.invPVM[i];
		}
		xdat = a.xdat;
	    return *this;
	}
}; /**< Current data state */

struct compositeArguments{
	int width;
	int height;
	Movement mov;
	int miny;
	int maxy;
	PhaseBarrier barrier;
};

struct DataPartition{
	int datx;
	int daty;
	int datz;
	float xmin;
	float xmax;
	float ymin;
	float ymax;
	float zmin;
	float zmax;
}; /**< Volumetric data partition bounding */

struct Image{
	int width;
	int height;
	float invPVM[16];
	DataPartition partition;
	float order;
	int randomseed;
	char volumeFilename[100];

	bool operator<( const Image& rhs ) const
	{ return order < rhs.order; }
}; /**< Individual image metadata structure (One per leaf render) */


class CompositeMapper : public DefaultMapper {
public:
	CompositeMapper(Machine machine, HighLevelRuntime *rt, Processor local);
public:
	virtual void select_task_options(Task *task);
//	virtual void slice_domain(const Task *task, const Domain &domain, std::vector<DomainSplit> &slices);
	virtual bool map_task(Task *task);
//	virtual bool map_inline(Inline *inline_operation);
protected:
	Processor top_proc;
	std::set<Processor> task_procs;
	std::map<Processor, Memory> all_sysmems;
	std::vector<Processor> all_cpus;
	std::vector<Processor> all_gpus;
};


#endif
