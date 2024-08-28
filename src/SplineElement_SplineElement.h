#pragma once
#include "SplineElementPair.h"

class SplineElement_SplineElement :
	public SplineElementPair
{
public:
	SplineElement_SplineElement();
	~SplineElement_SplineElement();
	//Chute inicial para coordenadas convectivas do par de superfícies
	void InitialGuess(SPContactData* c_data);
				
	double ObjectivePhase1(Matrix& mc);								//Calcula a função objetivo para um conjunto de coordenadas convectivas - Phase 1
	void GradientPhase1(Matrix& mc, Matrix& mGra);					//Calcula o Gradiente da função objetivo - Phase 1
	void HessianPhase1(Matrix& mc, Matrix& mHes);					//Calcula a Hessiana da função objetivo - Phase 1
	
	int VerifyConvectiveRange(Matrix& mc);							//Verifica range de coordenadas convectivas
	void InitializeConvectiveRange();								//Initialize range of validity of convective coordinates

	void ContactSS(bool *stick, bool *stickupdated, bool *previouscontact, double* Rc, double** Kc, double** invH, double* convective, double* copy_convective, double* gti, double* gtpupdated, double* epsn, double* epsn_n, double* epst, double* cn, double* ct, double* mus, double* mud, double* fn, double* ft);
};