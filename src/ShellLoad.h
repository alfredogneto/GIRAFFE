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
	void UpdateforSolutionStep();			//Atualiza dados necessarios e que sejam dependentes de DOFs ativos/inativos - chamado no inicio de cada solution step
	void Mount();
	bool Check();							//Checking inconsistencies
	void EvaluateExplicit(double t);
	//Variables
	int element_set;						//ElementSet de atuação
	bool area_update;						//flag - area update
};

