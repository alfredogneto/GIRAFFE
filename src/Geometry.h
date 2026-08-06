#pragma once
#include <stdio.h>

#define PI 3.1415926535897932384626433832795

class BoundingVolume;
class Matrix;

class Geometry
{
public:
	Geometry();
	virtual ~Geometry();
	void Alloc();
	void Free();
	virtual void WriteVTK_XMLRender(FILE *f) = 0;
	virtual void AllocSpecific() = 0;
	virtual void FreeSpecific() = 0;
	virtual bool Read(FILE *f) = 0;
	virtual void Write(FILE *f) = 0;
	virtual bool Check() = 0;
	virtual void PreCalc() = 0;

	virtual void UpdateVariables() = 0;
	virtual void UpdateBoundingVolumes() = 0;				//Updates bounding volumes
	virtual void SaveLagrange() = 0;						//Saves variables

	virtual void SurfacePosition(Matrix* Gamma, double z, double th, bool next) = 0;	//returns the surface position on the vector Gamma, evaluated at (z,th) on current (i) or next (i+1) configuration
	virtual void NormalExt(double z, double th, Matrix* n, bool next) = 0;							//returns the exterior normal of the surface
	//Degeneration (true or false for first and second convective coordinates - when applicable for surfaces)
	bool deg_u1, deg_u2, deg_u1_u2;
	//Degenerated coordinates (only considered if boolean degeneration values are TRUE)
	double deg_u1_value, deg_u2_value, deg_u1_u2_relation;

	BoundingVolume* bv;		//Bounding volume
	float inc_len_factor;	//increase factor do BV
	float bv_offset;		//offset do BV
	float max[3], min[3];	//coordinates max and min - global CS - to guide the BodyGeometry Axis Aligned BV
	float largest_gnb;		//for auxiliary in the calculation of the 'size' of the bounding volume - Verlet/LinkedCells schemes

	int number;				//ID of the Geometry
	int material;			//ID of the material
	int super_node;			//Super node (when aplicavel)
	int mother_entity;		//ID of the surface mãe usada for degeneration (for cases of hierarquia of contact pairs)
	char* type_name;		//Name of the type of Geometry
	int *nodes;				//List of nodes
	int n_nodes;			//Number of nodes
	int **DOFs;				//Indicates for the indexing of each degree of freedom, 1 or 0, active or inactive
	int nDOFs;				//Number of the GL (global)
	int **GLs;				//Pointer for the GL global used

	double* body_mass;		//Mass of the associated body (pointer)

	//Variables to store DOFs and derivatives
	Matrix* d;
	Matrix* dui;
	Matrix* ddui;
};