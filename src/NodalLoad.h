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
	void UpdateforSolutionStep();			//Atualiza dados necessarios e que sejam dependentes de DOFs ativos/inativos - chamado no inicio de cada solution step
	void Mount();
	void EvaluateExplicit(double t);
	bool Check();							//Checking inconsistencies

	//Variables
	int node_set;
	int cs;

	int* n_nodes_f;							//numero de nós para divisão de forças - 3 componentes
	int* n_nodes_m;							//numero de nós para divisão de momentos - 3 componentes
	double* mult_f;							//multiplicador para os esforços de força
	double* mult_m;							//multiplicador para os esforços de 

	//Matrix Q, Xi, alpha;
	//double g, alpha_escalar;
	//Matrix I, A;
	//Matrix** q;
	//Matrix dqdd;
};

