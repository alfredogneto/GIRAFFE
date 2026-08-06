#pragma once
#include "SpecialConstraint.h"
#include "Matrix.h"

class TranslationalJoint :
	public SpecialConstraint
{
public:
	TranslationalJoint();
	~TranslationalJoint();

	bool Read(FILE *f);		//Reading
	void Write(FILE *f);	//Recording
	void Mount();			//Assembly of the residuals and stiffness tangent
	void MountGlobal();		//Fills the contribution of the element in the matrices global
	void PreCalc();			//Pre-calculation of variables that and done a single time in the start
	void SaveLagrange();	//Saving variables of the configuration converged
	void ActivateDOFs();	//Checks which Lagrange multipliers will be ativados,according to the activation of the GLs of the nodes of the which the special constraint participates
	bool Check();			//Checks inconsistencies in the SC for avoid errors of execution
	void ClearContributions();			//Zeroes matrices and vectors
	void ComputeVelAccel();					//Computes effect of the initial conditions in the nodes of the constraint
	void ComputeInitialGuessDisplacements();
	void WriteVTK_XMLRender(FILE *f);
	void WriteMonitor(FILE *f, bool first_record, double time);	//Writes in the monitor of the SpecialConstraint
	//Variables
	int node_A;
	int node_B;
	int rot_node;	
	int cs;					//coordinate system - e3 will be the orientation with translation relativa allowed between A and B

	Matrix I3;
	
	Matrix alphaA;
	Matrix alphaiA;
	double alpha_escalar_i, g;
	Matrix A;
	Matrix QA;
	Matrix ei1A;
	Matrix ei2A;

	Matrix uA;
	Matrix uB;

	double** stiffness;
	double* residual;

	double temp_v[1000];	//For functions of the AceGen
	//Calculates contributions of the residual and operator tangent - generated in the AceGen
	void EvaluateTranslationalContribution(double *v, double *residual, double **stiffness, double *uA, double *uB, double *alphaA, double *ei1A, double *ei2A, double *lambda);
};

