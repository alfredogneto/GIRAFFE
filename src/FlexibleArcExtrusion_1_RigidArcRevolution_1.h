#pragma once
#include "SurfacePair.h"
class Matrix;

class FlexibleArcExtrusion_1_RigidArcRevolution_1 :
	public SurfacePair
{
public:
	FlexibleArcExtrusion_1_RigidArcRevolution_1();
	~FlexibleArcExtrusion_1_RigidArcRevolution_1();
	//Chute inicial para coordenadas convectivas do par de superfícies
	void InitialGuess(SSContactData* c_data);
	
	double ObjectivePhase1(Matrix& mc);								//Calcula a função objetivo para um conjunto de coordenadas convectivas - Phase 1
	void GradientPhase1(Matrix& mc, Matrix& mGra);					//Calcula o Gradiente da função objetivo - Phase 1
	void HessianPhase1(Matrix& mc, Matrix& mHes);					//Calcula a Hessiana da função objetivo - Phase 1
	
	double Gap(Matrix& mc, bool fixed_normals, Matrix& nA, Matrix& nB);											//Calcula e rotorna o gap (com sinal)
	void GradientGap(Matrix& mc, Matrix& mGra, bool fixed_normals, Matrix& nA, Matrix& nB);						//Calcula o Gradiente do gap
	void HessianGap(Matrix& mc, Matrix& mHes, bool fixed_normals, Matrix& nA, Matrix& nB);						//Calcula a Hessiana do gap
    int VerifyConvectiveRange(Matrix& mc);							//Verifica range de coordenadas convectivas
	void InitializeConvectiveRange();								//Initialize range of validity of convective coordinates
	bool SpecialLCP(Matrix& solution);

	void ContactSS(bool *stick, bool *stickupdated, bool *previouscontact, double* Rc, double** Kc, double** invH, double* convective, double* copy_convective, double* gti, double* gtpupdated, double* epsn, double* epsn0, double* epst, double* cn, double* ct, double* mus, double* mud, double* fn, double* ft);
};