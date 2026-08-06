#pragma once
#include "SplineElementPair.h"

class SplineElement_SplineElement :
	public SplineElementPair
{
public:
	SplineElement_SplineElement();
	~SplineElement_SplineElement();
	//Initial guess for convective coordinates of the pair of surfaces
	void InitialGuess(SPContactData* c_data);
				
	double ObjectivePhase1(const Matrix& mc);								//Calculates the objective function for a set of convective coordinates - Phase 1
	void GradientPhase1(Matrix& mc, Matrix& mGra);					//Calculates the Gradient of the objective function - Phase 1
	void HessianPhase1(Matrix& mc, Matrix& mHes);					//Calculates the Hessian of the objective function - Phase 1
	
	int VerifyConvectiveRange(Matrix& mc);							//Checks range of convective coordinates
	void InitializeConvectiveRange();								//Initialize range of validity of convective coordinates

	void ContactSS(bool *stick, bool *stickupdated, bool *previouscontact, double* Rc, double** Kc, double** invH, double* convective, double* copy_convective, double* gti, double* gtpupdated, double* epsn, double* epsn_n, double* epst, double* cn, double* ct, double* mus, double* mud, double* fn, double* ft);
};