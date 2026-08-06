#pragma once
#include "Load.h"

class ShellLoad :
	public Load
{
public:
	ShellLoad();
	~ShellLoad();

	bool Read(FILE *f);						//Reads input file
	void Write(FILE *f);					//Writes output file
	void WriteVTK_XML(FILE *f);				//Writes VTK XML data for post-processing
	void PreCalc();							//Pre-calculus
	void UpdateforSolutionStep();			//Updates data necessarios and that be dependentes of DOFs active/inativos - called in the start of each solution step
	void Mount();
	bool Check();							//Checking inconsistencies
	void EvaluateExplicit(double t);
	//Variables
	int element_set;						//ElementSet of atuação
	bool area_update;						//flag - area update
};

