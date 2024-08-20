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
	void PreCalc();								//Pré-cálculo de variáveis que é feito uma única vez no início
	void UpdateBoundingVolumes();
	void SaveLagrange();						//Salva variáveis
	void WriteVTK_XMLBase(FILE *f);
	void WriteVTK_XMLRender(FILE *f);

	
	double radius;
	double volume;									//Volume
	double inertia;									//inércia à rotação

	//Rigid body variables
	double mass;													//Massa
	double** Jr;													//Tensor de inércia - formato double**
	double* br;														//Vetor br - formato double*
	double** DdT;													//Operador tangente
	double* dT;														//Resíduo
	double** Ddfield;												//Carregamento de campo linearizado
	double* dfield;													//Resíduo carregamento de campo
	//Variáveis cinemáticas
	double alphai[3];
	double ud[3];
	double alphad[3];
	double dui[3];
	double ddui[3];
	double omegai[3];
	double domegai[3];
};

