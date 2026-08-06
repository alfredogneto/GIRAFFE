#pragma once
#include "ContactBodyBody.h"
class ContactArcExtrusionArcRevolution :
	public ContactBodyBody
{
public:
	ContactArcExtrusionArcRevolution();
	~ContactArcExtrusionArcRevolution();

	void Alloc();
	void Free();
	void MountLocalContributions();
	void SetVariables();					//Sets variables for AceGen codes interfaces
	void Report();
	void CompactReport();
	void InitialGuess();

	double ObjectivePhase1(const Matrix& mc);								//Calculates the objective function for a set of convective coordinates - Phase 1
	void GradientPhase1(Matrix& mc, Matrix& mGra);					//Calculates the Gradient of the objective function - Phase 1
	void HessianPhase1(Matrix& mc, Matrix& mHes);					//Calculates the Hessian of the objective function - Phase 1
	double Gap(const Matrix& mc, bool fixed_normals, Matrix& nA, Matrix& nB);											//Calculates and rotorna the gap (with sign)
	void GradientGap(Matrix& mc, Matrix& mGra, bool fixed_normals, Matrix& nA, Matrix& nB);						//Calculates the Gradient of the gap
	void HessianGap(Matrix& mc, Matrix& mHes, bool fixed_normals, Matrix& nA, Matrix& nB);						//Calculates the Hessian of the gap
	int VerifyConvectiveRange(Matrix& mc);							//Checks range of convective coordinates
	void InitializeConvectiveRange();								//Initialize range of validity of convective coordinates

	void PrintAceGenPointers();

	double* cAp;
	double* cBp;
	double* cAi;
	double* cBi;

	//AceGen Pointers - specific
	double* radA;                                             //Radius of curvature of the arc
	double* cpointA;											//Center of curvature of the arc
	double* radB;                                             //Radius of curvature of the arc
	double* cpointB;											//Center of curvature of the arc

	bool* normalintA;
	bool* normalintB;
	double* xAAi;
	double* xBAi;
	double* xABi;
	double** QAAi;
	double** QBAi;
	double** QABi;

	double* gti;
	double* gtpupdated;
	bool *stick;
	bool *stickupdated;
	double** invH;
	bool *interfacelaw0;
};

