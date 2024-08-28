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
	void PreCalc();							//Realiza pre-calculos do esfor�o seguidor
	void UpdateforSolutionStep();			//Atualiza dados necessarios e que sejam dependentes de DOFs ativos/inativos - chamado no inicio de cada solution step
	void Mount();							//Calcula o esfor�o e o operador tangente
	bool Check();							//Checking inconsistencies
	void EvaluateExplicit(double t);
	//Variables
	int node_set;
	int cs;
	int n_nodes_copy;

	int* n_nodes_f;							//numero de n�s para divis�o de for�as - 3 componentes
	int* n_nodes_m;							//numero de n�s para divis�o de momentos - 3 componentes
	double* mult_f;							//multiplicador para os esfor�os de for�a
	double* mult_m;							//multiplicador para os esfor�os de momento

	//Variaveis utilizadas nas fun��es
	Matrix Q, Xi, alpha;	
	double g,alpha_escalar;
	Matrix I, A;
	Matrix f;
	Matrix m;
	Matrix** q;
	Matrix dqdd;
};
