#pragma once
#include <stdio.h>

class Matrix;
class BoundingVolume;

class Boundary
{
public:
	Boundary();
	virtual ~Boundary();

	int node;				//Node of reference
	int number;				//ID of the boundary
	int cs;					//ID of the coordinate system
	int material;			//ID of the material of the boundary (for assignment of leis of interface)

	void Alloc();
	void Free();
	void UpdateVariables();

	//Functions especificas for each type of particle
	virtual bool Read(FILE *f) = 0;
	virtual void Write(FILE *f) = 0;
	virtual void WriteModifyingParameters(FILE *f, int e_material, int e_node, int e_number, int e_cs) = 0;
	virtual bool Check() = 0;								//Checks inconsistencies in the particle for avoid errors of execution

	virtual void PreCalc() = 0;								//Pre-calculation of variables that and done a single time in the start
	virtual void UpdateBoundingVolumes() = 0;				//Updates bounding volumes
	virtual void SaveLagrange() = 0;						//Saves variables
	virtual void WriteVTK_XMLBase(FILE *f) = 0;
	virtual void WriteVTK_XMLRender(FILE *f) = 0;

	//Variables - bounding volumes
	BoundingVolume* bv;												//Bounding volume
	int n_sub_bv;													//Number of sub bounding volumes
	BoundingVolume** sub_bv;										//Sub bounding volumes

	Matrix* I3;
	Matrix* Q0;														//Transformation of coordinates
	//Variables for calculate state current (in the functions of contact)
	Matrix* Qip;
	Matrix* x0ip;
	//Variables for calculate state previous (in the functions of contact)
	Matrix* Qi;
	Matrix* x0i;

	//Explicit
	virtual void InitialEvaluations() = 0;
};

