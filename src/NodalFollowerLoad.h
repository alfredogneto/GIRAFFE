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
	void PreCalc();							//Realiza pré-cálculos do esforço seguidor
	void UpdateforSolutionStep();			//Atualiza dados necessários e que sejam dependentes de DOFs ativos/inativos - chamado no início de cada solution step
	void Mount();							//Calcula o esforço e o operador tangente
	bool Check();							//Checking inconsistencies
	void EvaluateExplicit(double t);
	//Variables
	int node_set;
	int cs;
	int n_nodes_copy;

	int* n_nodes_f;							//número de nós para divisão de forças - 3 componentes
	int* n_nodes_m;							//número de nós para divisão de momentos - 3 componentes
	double* mult_f;							//multiplicador para os esforços de força
	double* mult_m;							//multiplicador para os esforços de momento

	//Variáveis utilizadas nas funções
	Matrix Q, Xi, alpha;	
	double g,alpha_escalar;
	Matrix I, A;
	Matrix f;
	Matrix m;
	Matrix** q;
	Matrix dqdd;
};
