#pragma once
#include "SurfacePair.h"
class Matrix;

class FlexibleSECylinder_1_FlexibleSECylinder_1 :
	public SurfacePair
{
public:
	FlexibleSECylinder_1_FlexibleSECylinder_1();
	~FlexibleSECylinder_1_FlexibleSECylinder_1();
	//Initial guess for convective coordinates of the pair of surfaces
	void InitialGuess(SSContactData* c_data);
				
	double ObjectivePhase1(const Matrix& mc);								//Calculates the objective function for a set of convective coordinates - Phase 1
	void GradientPhase1(Matrix& mc, Matrix& mGra);					//Calculates the Gradient of the objective function - Phase 1
	void HessianPhase1(Matrix& mc, Matrix& mHes);					//Calculates the Hessian of the objective function - Phase 1
	
	double Gap(const Matrix& mc, bool fixed_normals, Matrix& nA, Matrix& nB);											//Calculates and rotorna the gap (with sign)
	void GradientGap(Matrix& mc, Matrix& mGra, bool fixed_normals, Matrix& nA, Matrix& nB);						//Calculates the Gradient of the gap
	void HessianGap(Matrix& mc, Matrix& mHes, bool fixed_normals, Matrix& nA, Matrix& nB);						//Calculates the Hessian of the gap
	int VerifyConvectiveRange(Matrix& mc);							//Checks range of convective coordinates
	void InitializeConvectiveRange();								//Initialize range of validity of convective coordinates

	void ContactSS(bool *stick, bool *stickupdated, bool *previouscontact, double* Rc, double** Kc, double** invH, double* convective, double* copy_convective, double* gti, double* gtpupdated, double* epsn, double* epsn0, double* epst, double* cn, double* ct, double* mus, double* mud, double* fn, double* ft);
};