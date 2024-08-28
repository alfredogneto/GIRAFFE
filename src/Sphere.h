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
	void PreCalc();								//Pre-calculo de variaveis que e feito uma unica vez no inicio
	void UpdateBoundingVolumes();
	void SaveLagrange();						//Salva variaveis
	void WriteVTK_XMLBase(FILE *f);
	void WriteVTK_XMLRender(FILE *f);

	
	double radius;
	double volume;									//Volume
	double inertia;									//inercia a rotação

	//Rigid body variables
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

