#pragma once
#include "SpecialConstraint.h"
#include "Matrix.h"

class UniversalJoint :
	public SpecialConstraint
{
public:
	UniversalJoint();
	~UniversalJoint();

	bool Read(FILE *f);		//Leitura
	void Write(FILE *f);	//Gravação
	void Mount();			//Montagem dos residuos e rigidez tangente
	void MountGlobal();		//Preenche a contribuição do elemento nas matrizes globais
	void PreCalc();			//Pre-calculo de variaveis que e feito uma unica vez no inicio
	void SaveLagrange();	//Salvando variaveis da configuração convergida
	void ActivateDOFs();	//Checa quais multiplicadores de lagrange serão ativados,de acordo com a ativação dos GLs dos nós dos quais a special constraint participa
	bool Check();			//Checa inconsistências no SC para evitar erros de execução
	void ComputeVelAccel();		//Computa efeito das condições iniciais nos nós da restrição
	void ComputeInitialGuessDisplacements();
	void WriteVTK_XMLRender(FILE *f);
	void WriteMonitor(FILE *f, bool first_record, double time);	//Escreve no monitor do SpecialConstraint
	//Variaveis
	int node_A;
	int node_B;
	int csA;//coordinate system - e3 sera a orientação em que a rotação e transmitida
	int csB;//coordinate system - e3 sera a orientação em que a rotação e transmitida

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

	//Matriz de rigidez tangente e vetor residuo
	Matrix* stiffness1;
	Matrix* residual1;
	double** stiffness2;
	double* residual2;
	double* temp_lambda;

	double temp_v[1000];	//Para funções do AceGen
	//Calcula contribuições do residuo e operador tangente - gerado no AceGen
	void EvaluateUniversalJointContribution(double *v, double *residual
		, double **stiffness, double *alphaA, double *alphaB, double *ei1A
		, double *ei2B, double *lambda);
	//Calcula contribuições do residuo e operador tangente - gerado no AceGen - sem usar SMSD
	void EvaluateUniversalJointContribution2(double *v, double *residual
		, double **stiffness, double *alphaA, double *alphaB, double *ei1A
		, double *ei2B, double *lambda);
};

