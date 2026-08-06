#pragma once
#include <stdio.h>

class Matrix;
class BoundingVolume;

class Particle
{
public:
	Particle();
	virtual ~Particle();
	int number;				//ID of the particle
	int material;			//ID of the material
	int cs;					//ID of the coordinate system of the particle
	char* type_name;		//Name of the type of particle

	int node;				//Node of reference
	int *DOFs;				//Indicates for the indexing of each degree of freedom, 1 or 0, active or inactive for the element in question
	int super_node;			//Super node associated

	//Functions especificas for each type of particle
	virtual void Alloc() = 0;
	virtual void Free() = 0;
	virtual void Mount() = 0;
	virtual void MountGlobal() = 0;
	virtual void MountFieldLoading() = 0;
	virtual void InertialContributions() = 0;
	virtual void UpdateVariables() = 0;
	//Explicit
	virtual void EvaluateExplicit() = 0;
	virtual void EvaluateAccelerations() = 0;
	Matrix* MassMatrix;
	Matrix* invMAA;
	Matrix* MAB;
	Matrix* MBB;
	int* GLA; //connectivity
	int* GLB; //connectivity
	bool allocedMassFragmented;
	Matrix* ResidualVector;
	virtual void InitialEvaluations() = 0;
	void AllocMassFragmented(int nA, int nB);
	void FreeMassFragmented();

	virtual bool Read(FILE *f) = 0;
	virtual void Write(FILE *f) = 0;
	void WriteMonitor(FILE *f, bool first_record, double time);	//Writes in the monitor of the particle

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
	//Transformation of coordinates
	Matrix* Q0;
	//Post processing
	double strain_energy;
	double kinetic_energy;
	double potential_g_energy;
	double angular_momentum[3];
	double angular_momentum_origin[3];
	double angular_momentum_mag;
	double linear_momentum[3];
	double linear_momentum_mag;
};

