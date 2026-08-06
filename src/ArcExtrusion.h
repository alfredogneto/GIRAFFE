#pragma once
#include "Geometry.h"

class MatrixFloat;

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
	void SaveLagrange();						//Saves variables

	//Variables - ArcExtrusion
	int arc_ID;													//ID of the arc the be extrudado
	int cs;                                                     //ID of the coordinate system for posicionar the arc in the espaço
	double* radius;                                             //Radius of curvature of the arc
	Matrix* c_point;											//Center of curvature of the arc
	Matrix* i_point;											//Point initial of the arc
	Matrix* f_point;											//Point final of the arc
	double* theta_i;
	double* theta_f;

	//Variables - SECylinder - Bounding volume (cylinder)
	float BVradius;
	MatrixFloat* x_local;
	
	//Variables for escolha of the lado concavo or convexo
	bool flag_normal_int;									//Flag for indicate use of normal interior/exterior

	//Variables for calculate state current (in the functions of bounding volume)
	MatrixFloat* xAf;
	MatrixFloat* xBf;

	//Variables for calculate state current (in the functions of contact)
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

