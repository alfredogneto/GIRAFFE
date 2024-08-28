#pragma once
#include <stdio.h>

class Matrix;
class BoundingVolume;

class Particle
{
public:
	Particle();
	virtual ~Particle();
	int number;				//ID da part�cula
	int material;			//ID do material
	int cs;					//ID do sistema de coordenadas da part�cula
	char* type_name;		//Nome do tipo de part�cula

	int node;				//N� de refer�ncia
	int *DOFs;				//Indica para a indexa��o de cada grau de liberdade, 1 ou 0, ativo ou inativo para o elemento em quest�o
	int super_node;			//Super node associated

	//Fun��es espec�ficas para cada tipo de part�cula
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
	void WriteMonitor(FILE *f, bool first_record, double time);	//Escreve no monitor da part�cula

	virtual void WriteModifyingParameters(FILE *f, int e_material, int e_node, int e_number, int e_cs) = 0;
	virtual bool Check() = 0;								//Checa inconsist�ncias na part�cula para evitar erros de execu��o

	virtual void PreCalc() = 0;								//Pr�-c�lculo de vari�veis que � feito uma �nica vez no in�cio
	virtual void UpdateBoundingVolumes() = 0;				//Updates bounding volumes
	virtual void SaveLagrange() = 0;						//Salva vari�veis
	virtual void WriteVTK_XMLBase(FILE *f) = 0;
	virtual void WriteVTK_XMLRender(FILE *f) = 0;

	//Vari�veis - bounding volumes
	BoundingVolume* bv;												//Bounding volume
	int n_sub_bv;													//Number of sub bounding volumes
	BoundingVolume** sub_bv;										//Sub bounding volumes

	

	Matrix* I3;
	//Transforma��o de coordenadas
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

