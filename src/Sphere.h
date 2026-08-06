#pragma once
#include "Particle.h"

class Sphere :
	public Particle
{
public:
	Sphere();
	~Sphere();

	bool Read(FILE *f);
	void Write(FILE *f);
	void WriteModifyingParameters(FILE *f, int e_material, int e_node, int e_number, int e_cs);
	bool Check();

	//Explicit
	void EvaluateExplicit();
	void EvaluateAccelerations();
	void InitialEvaluations();

	void Alloc();
	void Free();
	void Mount();
	void MountGlobal();
	void MountFieldLoading();
	void InertialContributions();
	void UpdateVariables();
	void PreCalc();								//Pre-calculation of variables that and done a single time in the start
	void UpdateBoundingVolumes();
	void SaveLagrange();						//Saves variables
	void WriteVTK_XMLBase(FILE *f);
	void WriteVTK_XMLRender(FILE *f);

	
	double radius;
	double volume;									//Volume
	double inertia;									//inertia the rotation

	//Rigid body variables
	double mass;													//Mass
	double** Jr;													//Tensor of inertia - formato double**
	double* br;														//Vector br - formato double*
	double** DdT;													//Operator tangent
	double* dT;														//Residual
	double** Ddfield;												//Loading of field linearizado
	double* dfield;													//Residual loading of field
	//Variables cinematicas
	double alphai[3];
	double ud[3];
	double alphad[3];
	double dui[3];
	double ddui[3];
	double omegai[3];
	double domegai[3];
};

