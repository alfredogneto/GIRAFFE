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
	void PreCalc();								//Pr�-c�lculo de vari�veis que � feito uma �nica vez no in�cio
	void UpdateBoundingVolumes();
	void SaveLagrange();						//Salva vari�veis
	void WriteVTK_XMLBase(FILE *f);
	void WriteVTK_XMLRender(FILE *f);

	
	double radius;
	double volume;									//Volume
	double inertia;									//in�rcia � rota��o

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

