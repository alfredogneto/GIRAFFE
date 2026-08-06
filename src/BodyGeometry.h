#pragma once
#include <stdio.h>

class Geometry;
class BoundingVolume;

class BodyGeometry
{
public:
	BodyGeometry();
	~BodyGeometry();

	int n_items;			//number of items
	int* list_items;		//list 
	int number;				//reference number

	bool sequence;			//true if and of the type sequence
	bool list;				//true if and of the type list

	//For the case of sequence
	int initial;
	int increment;

	bool Read(FILE *f);
	void Write(FILE *f);

	BoundingVolume* bv;			//Bounding volume
	Geometry** ptr_geom;		//Pointer to geometries
	float inc_len_factor;
	float max_offset;

	void WriteVTK_XMLRender(FILE *f);
	bool Check();
	void PreCalc();
	void UpdateVariables();
	void UpdateBoundingVolumes();				//Updates bounding volumes
	void SaveLagrange();						//Saves variables

	double mass;								//Mass of the Body - for compute contact damping
};

