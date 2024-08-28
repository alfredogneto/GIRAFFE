#pragma once
#include <stdio.h>

class Geometry;
class BoundingVolume;

class BodyGeometry
{
public:
	BodyGeometry();
	~BodyGeometry();

	int n_items;			//número de items
	int* list_items;		//lista 
	int number;				//número de referência

	bool sequence;			//true se é do tipo sequence
	bool list;				//true se é do tipo list

	//Para o caso de sequence
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
	void SaveLagrange();						//Salva variáveis

	double mass;								//Massa do Body - para computar contact damping
};

