#pragma once
#include "SurfacePair.h"

class FlexibleSECylinder_1;
class FlexibleTriangularSurface_2;
class Matrix;

class FlexibleSECylinder_1_FlexibleTriangularSurface_2 :
	public SurfacePair
{
public:
	FlexibleSECylinder_1_FlexibleTriangularSurface_2();
	~FlexibleSECylinder_1_FlexibleTriangularSurface_2();

	//Chute inicial para coordenadas convectivas do par de superfícies
	void InitialGuess(SSContactData* c_data);
	
	void ContactSS(bool *stick, bool *stickupdated, bool *previouscontact, double* Rc, double** Kc, double** invH, double* convective, double* copy_convective, double* gti, double* gtpupdated, double* epsn, double* epsn0, double* epst, double* cn, double* ct, double* mus, double* mud, double* fn, double* ft);
	double ObjectivePhase1(Matrix& mc);								//Calcula a função objetivo para um conjunto de coordenadas convectivas - Phase 1
	void GradientPhase1(Matrix& mc, Matrix& mGra);					//Calcula o Gradiente da função objetivo - Phase 1
	void HessianPhase1(Matrix& mc, Matrix& mHes);					//Calcula a Hessiana da função objetivo - Phase 1
	
	double Gap(Matrix& mc, bool fixed_normals, Matrix& nA, Matrix& nB);											//Calcula e rotorna o gap (com sinal)
	void GradientGap(Matrix& mc, Matrix& mGra, bool fixed_normals, Matrix& nA, Matrix& nB);						//Calcula o Gradiente do gap
	void HessianGap(Matrix& mc, Matrix& mHes, bool fixed_normals, Matrix& nA, Matrix& nB);						//Calcula a Hessiana do gap
	int VerifyConvectiveRange(Matrix& mc);							//Verifica range de coordenadas convectivas
	void InitializeConvectiveRange();								//Initialize range of validity of convective coordinates

	void EvaluateDOFsVariables();									//Evaluates the necessary DOFs variables to be used as input to AceGen routines
	//AceGen Variables or Pointers
	FlexibleSECylinder_1* surf1;			//Ponteiro para a superfície 1
	double* aA;
	double* bA;
	double* eA;
	double* dA;
	double* duiA;
	double* dduiA;
	bool* normalintA;
	double* xAAi;
	double* xBAi;
	double QAAi[3][3];
	double QBAi[3][3];
	FlexibleTriangularSurface_2* surf2;		//Ponteiro para a superfície 2
	double* dB;
	double* duiB;
	double* dduiB;
	double* xAi;
	double* xBi;
	double* xCi;
	double* xDi;
	double* xEi;
	double* xFi;
};

