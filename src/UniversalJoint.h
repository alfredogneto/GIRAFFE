#pragma once
#include "SpecialConstraint.h"
#include "Matrix.h"

class UniversalJoint :
	public SpecialConstraint
{
public:
	UniversalJoint();
	~UniversalJoint();

	bool Read(FILE *f);		//Reading
	void Write(FILE *f);	//Recording
	void Mount();			//Assembly of the residuals and stiffness tangent
	void MountGlobal();		//Fills the contribution of the element in the matrices global
	void PreCalc();			//Pre-calculation of variables that and done a single time in the start
	void SaveLagrange();	//Saving variables of the configuration converged
	void ActivateDOFs();	//Checks which Lagrange multipliers will be ativados,according to the activation of the GLs of the nodes of the which the special constraint participates
	bool Check();			//Checks inconsistencies in the SC for avoid errors of execution
	void ComputeVelAccel();		//Computes effect of the initial conditions in the nodes of the constraint
	void ComputeInitialGuessDisplacements();
	void WriteVTK_XMLRender(FILE *f);
	void WriteMonitor(FILE *f, bool first_record, double time);	//Writes in the monitor of the SpecialConstraint
	//Variables
	int node_A;
	int node_B;
	int csA;//coordinate system - e3 will be the orientation in that the rotation and transmitida
	int csB;//coordinate system - e3 will be the orientation in that the rotation and transmitida

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
	Matrix ei1A;
	Matrix ei2B;

	//Stiffness matrix tangent and vector residual
	Matrix* stiffness1;
	Matrix* residual1;
	double** stiffness2;
	double* residual2;
	double* temp_lambda;

	double temp_v[1000];	//For functions of the AceGen
	//Calculates contributions of the residual and operator tangent - generated in the AceGen
	void EvaluateUniversalJointContribution(double *v, double *residual
		, double **stiffness, double *alphaA, double *alphaB, double *ei1A
		, double *ei2B, double *lambda);
	//Calculates contributions of the residual and operator tangent - generated in the AceGen - without use SMSD
	void EvaluateUniversalJointContribution2(double *v, double *residual
		, double **stiffness, double *alphaA, double *alphaB, double *ei1A
		, double *ei2B, double *lambda);
};

