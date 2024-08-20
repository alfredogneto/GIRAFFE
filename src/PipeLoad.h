#pragma once
#include "Load.h"

class PipeLoad :
	public Load
{
public:
	PipeLoad();
	~PipeLoad();

	bool Read(FILE *f);						//Reads input file
	void Write(FILE *f);					//Writes output file
	void WriteVTK_XML(FILE *f);				//Writes VTK XML data for post-processing
	void PreCalc();							//Pre-calculus
	void UpdateforSolutionStep();			//Atualiza dados necess�rios e que sejam dependentes de DOFs ativos/inativos - chamado no in�cio de cada solution step
	void Mount();
	bool Check();							//Checking inconsistencies
	void EvaluateExplicit(double t);
	//Variables
	int element_set;						//ElementSet de atua��o
};

