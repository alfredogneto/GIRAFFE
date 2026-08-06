#pragma once
#include "Load.h"

class NodalLoad :
	public Load
{
public:
	NodalLoad();
	~NodalLoad();
	bool Read(FILE *f);						//Reads input file
	void Write(FILE *f);					//Writes output file
	void WriteVTK_XML(FILE *f);				//Writes VTK XML data for post-processing
	void PreCalc();							//Pre-calculus
	void UpdateforSolutionStep();			//Updates data necessarios and that be dependentes of DOFs active/inativos - called in the start of each solution step
	void Mount();
	void EvaluateExplicit(double t);
	bool Check();							//Checking inconsistencies

	//Variables
	int node_set;
	int cs;

	int* n_nodes_f;							//number of nodes for divisão of forces - 3 components
	int* n_nodes_m;							//number of nodes for divisão of moments - 3 components
	double* mult_f;							//multiplicador for the internal forces of force
	double* mult_m;							//multiplicador for the internal forces of moment
};

