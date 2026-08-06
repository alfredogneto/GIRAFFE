#pragma once
#include "Particle.h"

class MatrixFloat;

class Polyhedron :
	public Particle
{
public:
	Polyhedron();
	~Polyhedron();

	int CADDATA_ID;

	bool Read(FILE *f);
	void Write(FILE *f);
	void WriteModifyingParameters(FILE *f, int e_material, int e_node, int e_number, int e_cs);
	bool Check();

	//Explicit
	void EvaluateExplicit();
	void EvaluateAccelerations();
	void InitialEvaluations();

	bool CheckInside(Matrix& point, int other_particle_ID);
	bool CheckInsideEdge(Matrix& point, int edge_ID, double tol);

	void Alloc();
	void Free();
	void Mount();
	void MountGlobal();
	void MountFieldLoading();
	void InertialContributions();
	void UpdateVariables();
	void PreCalc();
	void UpdateBoundingVolumes();
	void SaveLagrange();
	void WriteVTK_XMLBase(FILE *f);
	void WriteVTK_XMLRender(FILE *f);

	//float bv_factor;		//Controls the size of bounding volumes of edges and vertices
	float inc_len_factor;	//Controls inflation of bounding volumes

	MatrixFloat* x0f;
	MatrixFloat* Q0f;

	//Variables for calculate state current (in the functions of contact)
	Matrix* Qip;
	Matrix* x0ip;
	//Variables for calculate state previous (in the functions of contact)
	Matrix* Qi;
	Matrix* x0i;

	//AceGen Mirror variables
	double **pQi;
	

	//Rigid body variables

	Matrix* mJr;													//Tensor of inertia - formato Matrix
	Matrix* mbr;													//Vector br - formato Matrix

	Matrix* mJrlocal;												//Tensor of inertia - system local
	Matrix* mbrlocal;												//Vector br - system local

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
