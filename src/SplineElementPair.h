#pragma once
#include <stdio.h>

#include "Matrix.h"

class SPContactData;

class SplineElementPair
{
public:
	SplineElementPair();
	virtual ~SplineElementPair();
	//Initial guess for convective coordinates of the pair of surfaces
	virtual void InitialGuess(SPContactData* c_data) = 0;
	//Single function for both contributions
	virtual void ContactSS(bool *stick, bool *stickupdated, bool *previouscontact, double* Rc, double** Kc, double** invH, double* convective, double* copy_convective, double* gti, double* gtpupdated, double* epsn, double* epsn_n, double* epst, double* cn, double* ct, double* mus, double* mud, double* fn, double* ft) = 0;
	
	virtual double ObjectivePhase1(const Matrix& mc) = 0;								//Calculates the objective function for a set of convective coordinates - Phase 1
	virtual void GradientPhase1(Matrix& mc, Matrix& mGra) = 0;					//Calculates the Gradient of the objective function - Phase 1
	virtual void HessianPhase1(Matrix& mc, Matrix& mHes) = 0;					//Calculates the Hessian of the objective function - Phase 1

	virtual int VerifyConvectiveRange(Matrix& mc) = 0;							//Checks range of convective coordinates
	virtual void InitializeConvectiveRange() = 0;								//Initialize range of validity of convective coordinates

	void WriteConvectiveRange();						//Writes convective range in the report
	void PreCalc();																										//Pre calculation
	void Alloc(SPContactData* c_data);																				//Allocates memory
	void Free();																										//Frees memory											
	void DefaultValues();																								//Values Default of tolerancias and other variables

	void EvaluateInvertedHessian(SPContactData* c_data);																//Calculates the inversa of the Hessian
	bool FindMinimumSolution(SPContactData* c_data, Matrix* solution, int &return_info);								//Optimization - determination of minimum
	bool FindMinimumSolutionDegenerated(SPContactData* c_data, Matrix* P_0, Matrix* solution);							//Optimization - determination of minimum

	void BeginStepCheck(SPContactData* c_data);
	bool EndStepCheck(SPContactData* c_data);
	void SolveLCP(SPContactData* c_data);																				//Solves problem local of contact

	int CharacterizeCriticalPoint(Matrix* solution);
	int CharacterizeCriticalPointDegenerated(Matrix* solution, Matrix* P_0, bool print = false);

	//Variables internal
	int spline1_ID;			//ID of the spline 1
	int spline2_ID;			//ID of the spline 2
	int surf1_ID;			//ID do spline element 1
	int surf2_ID;			//ID do spline element 2
	
	bool alloc_control;
	Matrix** cNR1;			//Solution of the phase 1

	double tol_small_1;		//Criterion for number very small - residual != 0
	double tol_eig;			//Criterion for eigenvalue baixo
	double tol_convective;	//Criterion for maximum error in the convective coordinates
	double tol_ascent;
	int seq_number;

	int max_it_1;			//Maximum number of iterations for optimization - phase 1

	int n_pointwise;		//Number of interações pointwise (assigned in the constructor)
		
	Matrix convective_range;//Matrix that stores the ranges of convective coordinates validas for the surfaces
	Matrix convective_max;	//Matrix that stores the values maximos of convective coordinates
	Matrix convective_min;	//Matrix that stores the values minimos of convective coordinates
	double minimum_convective_range;
	//TR report
	FILE **f_TR_report;
	FILE **f_DEG_report;
	void OpenTRReport(int index);
	void OpenDEGReport(int index);
	void InitializeTRReport(int index);
	void InitializeDEGReport(int index);
	bool write_report;
	bool write_report_diverged;

	char name[1000];
	char name_deg[1000];
};

