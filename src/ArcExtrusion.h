#pragma once
#include "Geometry.h"
class ArcExtrusion :
	public Geometry
{
public:
	ArcExtrusion();
	~ArcExtrusion();

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

	//Variables - ArcExtrusion
	int arc_ID;													//ID do arco a ser extrudado
	int cs;                                                     //ID do sistema de coordenadas para posicionar o arco no espaço
	double* radius;                                             //Raio de curvatura do arco
	Matrix* c_point;											//Centro de curvatura do arco
	Matrix* i_point;											//Ponto inicial do arco
	Matrix* f_point;											//Ponto final do arco
	double* theta_i;
	double* theta_f;

	//Variables - SECylinder - Bounding volume (cylinder)
	float BVradius;
	MatrixFloat* x_local;
	
	//Variáveis para escolha do lado concavo ou convexo
	bool flag_normal_int;									//Flag para indicar uso de normal interior/exterior

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


	Matrix* Q0;
	MatrixFloat* Qf;
	Matrix* I3;

	//AceGen Mirror variables
	double *xAi, *xBi;
	double **QAi, **QBi;
};

