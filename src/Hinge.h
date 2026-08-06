#pragma once
#include "SpecialConstraint.h"
#include "Matrix.h"

class Hinge :
	public SpecialConstraint
{
public:
	Hinge();
	~Hinge();

	bool Read(FILE *f);		//Reading
	void Write(FILE *f);	//Recording
	void Mount();			//Assembly of the residuals and stiffness tangent
	void MountGlobal();		//Fills the contribution of the element in the matrices global
	void PreCalc();			//Pre-calculation of variables that and done a single time in the start
	void SaveLagrange();	//Saving variables of the configuration converged
	void ActivateDOFs();	//Checks which Lagrange multipliers will be ativados,according to the activation of the GLs of the nodes of the which the special constraint participates
	bool Check();			//Checks inconsistencies in the SC for avoid errors of execution
	void ClearContributions();			//Zeroes matrices and vectors
	void ComputeVelAccel();		//Computes effect of the initial conditions in the nodes of the constraint
	void ComputeInitialGuessDisplacements();
	void WriteVTK_XMLRender(FILE *f);
	void WriteMonitor(FILE *f, bool first_record, double time);	//Writes in the monitor of the SpecialConstraint
	//Variables
	int node_A;
	int node_B;
	int cs;//coordinate system - e3 will be the orientation of the articulação

	Matrix I3;
	Matrix r1;

	Matrix alphaA;
	Matrix alphaB;
	Matrix alphaiA;
	Matrix alphaiB;
	double alpha_escalar_i, g;
	Matrix A;
	Matrix QA;
	Matrix QB;
	Matrix ei3A;
	Matrix ei1B;
	Matrix ei2B;
	Matrix ei1A;
	Matrix ei2A;

	//Stiffness matrix tangent and vector residual
	Matrix* stiffness1;
	Matrix* residual1;
	double** stiffness2;
	double* residual2;
	double* temp_lambda;

	//Spring of torção
	double** stiffness_spring;
	double* residual_spring;
	double thetai, thetad, stiffc;

	//Amortecedor of torção
	double** stiffness_damper;
	double* residual_damper;
	double dampc1, dampc2;
	Matrix omegaiA;
	Matrix omegaiB;
	Matrix domegaiA;
	Matrix domegaiB;

	double temp_v[10000];	//For functions of the AceGen
	//Calculates contributions of the residual and operator tangent - generated in the AceGen
	void EvaluateHingeContribution(double *v, double *residual
		, double **stiffness, double *alphaA, double *alphaB, double *ei3A
		, double *ei1B, double *ei2B, double *lambda);
	//Calculates contributions of the residual and operator tangent - generated in the AceGen (without use SMSD)
	void EvaluateHingeContribution2(double *v, double *residual
		, double **stiffness, double *alphaA, double *alphaB, double *ei3A
		, double *ei1B, double *ei2B, double *lambda);
	//Calculates contributions of the spring of torção (spring)
	void EvaluateTorsionSpring(double *v, double *residual
		, double **stiffness, double *alphaA, double *alphaB, double *ei1A, double *ei1B, double *ei3A
		, double(*stiffc), double(*thetai)
		, double(*thetad));
	//Calculates contributions of the amortecedor of torção
	void EvaluateTorsionDamping(double *v, double *residual
		, double **stiffness, double *alphaA, double *alphaB, double *omegaiA
		, double *omegaiB, double *domegaiA, double *domegaiB, double *ei3A, double
		(*dampc1), double(*dampc2), double(*alpha4), double(*alpha5), double(*alpha6));
};

