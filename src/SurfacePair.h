#pragma once
#include <stdio.h>
#include "Matrix.h"

class SSContactData;

class SurfacePair
{
public:
	SurfacePair();
	virtual ~SurfacePair();
	//Initial guess for convective coordinates of the pair of surfaces
	virtual void InitialGuess(SSContactData* c_data) = 0;				
	//Single function for both contributions
	virtual void ContactSS(bool *stick, bool *stickupdated, bool *previouscontact, double* Rc, double** Kc, double** invH, double* convective, double* copy_convective, double* gti, double* gtpupdated, double* epsn, double* epsn0, double* epst, double* cn, double* ct, double* mus, double* mud, double* fn, double* ft) = 0;

	virtual double ObjectivePhase1(const Matrix& mc) = 0;								//Calculates the objective function for a set of convective coordinates - Phase 1
	virtual void GradientPhase1(Matrix& mc, Matrix& mGra) = 0;					//Calculates the Gradient of the objective function - Phase 1
	virtual void HessianPhase1(Matrix& mc, Matrix& mHes) = 0;					//Calculates the Hessian of the objective function - Phase 1
	
	virtual double Gap(const Matrix& mc, bool fixed_normals, Matrix& nA, Matrix& nB) = 0;											//Calculates and rotorna the gap (with sign)
	virtual void GradientGap(Matrix& mc, Matrix& mGra, bool fixed_normals, Matrix& nA, Matrix& nB) = 0;						//Calculates the Gradient of the gap
	virtual void HessianGap(Matrix& mc, Matrix& mHes, bool fixed_normals, Matrix& nA, Matrix& nB) = 0;						//Calculates the Hessian of the gap
	virtual int VerifyConvectiveRange(Matrix& mc) = 0;							//Checks range of convective coordinates
	virtual void InitializeConvectiveRange() = 0;								//Initialize range of validity of convective coordinates
	
	void WriteConvectiveRange();						//Writes convective range in the report
	void PreCalc();																										//Pre calculation
	void Alloc(SSContactData* c_data);																				//Allocates memory
	void Free();																										//Frees memory											
	void DefaultValues();																								//Values Default of tolerancias and other variables
	
	void EvaluateInvertedHessian(SSContactData* c_data);																//Calculates the inversa of the Hessian
	bool FindMinimumSolution(SSContactData* c_data, Matrix* solution, int &return_info);								//Optimization - determination of minimum
	bool FindMinimumSolutionDegenerated(SSContactData* c_data, Matrix* P_0, Matrix* solution);						//Optimization - determination of minimum
	bool FindSaddleSolution(SSContactData* c_data, Matrix* solution, int &return_info, bool return_gap);				//Optimization - determination of saddle
	bool FindSaddleSolutionDegenerated(SSContactData* c_data, Matrix* P_0, Matrix* solution, bool return_gap);		//Optimization - determination of saddle
	bool FindMinimumGapDegenerated(SSContactData* c_data, Matrix* P_0, Matrix* solution, int &return_info, bool fixed_normals, Matrix& nA, Matrix& nB);			//Optimization - determination of minimum of the gap with sign
	
	void BeginStepCheck(SSContactData* c_data);					
	bool EndStepCheck(SSContactData* c_data);
	void SolveLCP(SSContactData* c_data);																				//Solves problem local of contact
	
	int CharacterizeCriticalPoint(Matrix* solution);
	int CharacterizeCriticalPointDegenerated(Matrix* solution, Matrix* P_0, bool print = false);
	void AutomaticDegenerationProcedure();																				//Performs automatic degeneration according to eigenvalues of Hessian
	
	//Variables internal
	int surf1_ID;			//ID of the surface 1
	int surf2_ID;			//ID of the surface 2
	bool inverted;			//Indicates that there is inversão of the types of surface (not if applies when the pair and formado by surfaces of same type)
	bool alloc_control;
	Matrix** cNR1;			//Solution of the phase 1
	Matrix** cNR2;			//Solution of the phase 2
	Matrix** cdeg;			//Solution degenerated in the start of the increment
	
	double tol_small_1;		//Criterion for number very small - residual != 0
	double tol_eig;			//Criterion for eigenvalue baixo
	double tol_convective;	//Criterion for maximum error in the convective coordinates
	double tol_ascent;
	int seq_number;

	int max_it_1;			//Maximum number of iterations for optimization - phase 1
	int max_it_2;			//Maximum number of iterations for optimization - phase 2

	bool* flag_degenerated;	//Flag that indicates ocorrência of degeneration
	int n_pointwise;		//Number of interações pointwise (assigned in the constructor)
	
	double perc;			//Percentual of convective coordinate for consider longe or perto of the range (afeta return value of the function ConvectiveRange)
	Matrix convective_range;//Matrix that stores the ranges of convective coordinates validas for the surfaces
	Matrix convective_max;	//Matrix that stores the values maximos of convective coordinates
	Matrix convective_min;	//Matrix that stores the values minimos of convective coordinates
	double minimum_convective_range;	
	//TR report
	FILE **f_TR_report;
	void OpenTRReport(int index);
	void InitializeTRReport(int index);
	bool write_report;
	bool write_report_diverged;
	bool specialLCP;
	char name[1000];
};

