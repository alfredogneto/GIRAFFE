#pragma once
#include "Particle.h"
class NURBSParticle :
	public Particle
{
public:
	NURBSParticle();
	~NURBSParticle();

	int CADDATA_ID;

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
	void PreCalc();
	void UpdateBoundingVolumes();
	void SaveLagrange();
	void WriteVTK_XMLBase(FILE *f);
	void WriteVTK_XMLRender(FILE *f);

	//Vari�veis para calcular estado atual (nas funcoes de contato)
	Matrix* Qip;
	Matrix* x0ip;
	//Vari�veis para calcular estado anterior (nas funcoes de contato)
	Matrix* Qi;
	Matrix* x0i;

	//Rigid body variables
	double mass;													//Massa
	double** Jr;													//Tensor de in�rcia - formato double**
	double* br;														//Vetor br - formato double*
	double** DdT;													//Operador tangente
	double* dT;														//Res�duo
	double** Ddfield;												//Carregamento de campo linearizado
	double* dfield;													//Res�duo carregamento de campo
	//Vari�veis cinem�ticas
	double alphai[3];
	double ud[3];
	double alphad[3];
	double dui[3];
	double ddui[3];
	double omegai[3];
	double domegai[3];
};

