#pragma once
#include "Geometry.h"
class SECylinder :
	public Geometry
{
public:
	SECylinder();
	~SECylinder();
	void WriteVTK_XMLRender(FILE *f);
	void AllocSpecific();
	void FreeSpecific();
	bool Read(FILE *f);
	void Write(FILE *f);
	bool Check();
	void PreCalc();

	void SurfacePosition(Matrix* Gamma, double z, double th, bool next);
	void NormalExt(double z, double th, Matrix* n, bool next);

	void UpdateVariables();
	void UpdateBoundingVolumes();				//Updates bounding volumes
	void SaveLagrange();						//Salva variáveis

	//Variables - SECylinder
	double *a, *b, *n, *e;
	int csA, csB;
	bool flag_normal_int;

	//Variáveis para calcular estado atual (nas funcoes de bounding volume)
	MatrixFloat* xAf;
	MatrixFloat* xBf;

	//Variáveis para calcular estado atual (nas funcoes de contato)
	Matrix* x_Ai;
	Matrix* x_Ap;
	Matrix* x_Bi;
	Matrix* x_Bp;
	Matrix* Q_Ai; 
	Matrix* Q_Ap;
	Matrix* Q_Bi;
	Matrix* Q_Bp;

	float rad;
	
	Matrix* Q0A;
	Matrix* Q0B;

	Matrix* I3;

	//AceGen Mirror variables
	double *xAi, *xBi;
	double **QAi, **QBi;
};

