#pragma once
#include "Load.h"
#include "Matrix.h"

class NodalFollowerLoad :
	public Load
{
public:
	NodalFollowerLoad();
	~NodalFollowerLoad();

	bool Read(FILE *f);						//Reads input file
	void Write(FILE *f);					//Writes output file
	void WriteVTK_XML(FILE *f);//Writes VTK XML data for post-processing
	void PreCalc();							//Performs pre-calculations of the internal force seguidor
	void UpdateforSolutionStep();			//Updates data necessarios and that be dependentes of DOFs active/inativos - called in the start of each solution step
	void Mount();							//Calculates the internal force and the operator tangent
	bool Check();							//Checking inconsistencies
	void EvaluateExplicit(double t);
	//Variables
	int node_set;
	int cs;
	int n_nodes_copy;

	int* n_nodes_f;							//number of nodes for divisão of forces - 3 components
	int* n_nodes_m;							//number of nodes for divisão of moments - 3 components
	double* mult_f;							//multiplicador for the internal forces of force
	double* mult_m;							//multiplicador for the internal forces of moment

	//Variables used in the functions
	Matrix Q, Xi, alpha;	
	double g,alpha_escalar;
	Matrix I, A;
	Matrix f;
	Matrix m;
	Matrix** q;
	Matrix dqdd;
};
