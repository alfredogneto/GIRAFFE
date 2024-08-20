#pragma once
#include "Load.h"
class SuperNodalLoad :
	public Load
{
public:
	SuperNodalLoad();
	~SuperNodalLoad();

	int localDOF;							//Local DOF of the SuperNode
	int super_node_ID;						//super node number
	double load_value;						//load value (current)

	bool Read(FILE *f);						//Reads input file
	void Write(FILE *f);					//Writes output file
	void WriteVTK_XML(FILE *f);				//Writes VTK XML data for post-processing
	void PreCalc();							//Pre-calculus
	void UpdateforSolutionStep();			//Atualiza dados necessários e que sejam dependentes de DOFs ativos/inativos - chamado no início de cada solution step
	void Mount();
	bool Check();							//Checking inconsistencies

	void EvaluateExplicit(double t);
};

