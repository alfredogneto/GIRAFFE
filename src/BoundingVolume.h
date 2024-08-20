#pragma once
#include <math.h>
#include "MatrixFloat.h"

class BoundingVolume
{
public:
	BoundingVolume();
	virtual ~BoundingVolume();
	virtual void SaveConfiguration() = 0;	//Saves in "prev" variables the current variables to serve as next saved configuration
	virtual void Report() = 0;

	float bv_factor;
	float inc_len_factor;

	float factor_kin2;						//Employed in correction of the BV due to kinematics

	bool first_set;

	float x_center[3];						//Central position for linked cells and Verlet evaluation
	float last_center_change[3];			//Last change in center position of the BV (for Verlet evaluations)
	float size;								//Maximum size (ref) for linked cells evaluation

	//Linked list structure - (linked cells)
	BoundingVolume* next;
	BoundingVolume* previous;
	int i_loc, j_loc, k_loc;	//Indexes of the cell where this BV is located in

	//Association with entity:
	char associated_type;	//'P' (particle), 'B' (boundary), 'G' (geometry)
	int associated_ID;		//starts with 1
	int associated_sub_ID;	//starts with 1 (0 if N/A)
};
