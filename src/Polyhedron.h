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

	//Variaveis para calcular estado atual (nas funcoes de contato)
	Matrix* Qip;
	Matrix* x0ip;
	//Variaveis para calcular estado anterior (nas funcoes de contato)
	Matrix* Qi;
	Matrix* x0i;

	//AceGen Mirror variables
	double **pQi;
	

	//Rigid body variables

	Matrix* mJr;													//Tensor de inercia - formato Matrix
	Matrix* mbr;													//Vetor br - formato Matrix

	Matrix* mJrlocal;												//Tensor de inercia - sistema local
	Matrix* mbrlocal;												//Vetor br - sistema local

	double mass;													//Massa
	double** Jr;													//Tensor de inercia - formato double**
	double* br;														//Vetor br - formato double*
	double** DdT;													//Operador tangente
	double* dT;														//Residuo
	double** Ddfield;												//Carregamento de campo linearizado
	double* dfield;													//Residuo carregamento de campo
	//Variaveis cinematicas
	double alphai[3];
	double ud[3];
	double alphad[3];
	double dui[3];
	double ddui[3];
	double omegai[3];
	double domegai[3];
};
