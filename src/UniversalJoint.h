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
	void Write(FILE *f);	//Grava��o
	void Mount();			//Montagem dos residuos e rigidez tangente
	void MountGlobal();		//Preenche a contribui��o do elemento nas matrizes globais
	void PreCalc();			//Pre-calculo de variaveis que e feito uma unica vez no inicio
	void SaveLagrange();	//Salvando variaveis da configura��o convergida
	void ActivateDOFs();	//Checa quais multiplicadores de lagrange ser�o ativados,de acordo com a ativa��o dos GLs dos n�s dos quais a special constraint participa
	bool Check();			//Checa inconsist�ncias no SC para evitar erros de execu��o
	void ComputeVelAccel();		//Computa efeito das condi��es iniciais nos n�s da restri��o
	void ComputeInitialGuessDisplacements();
	void WriteVTK_XMLRender(FILE *f);
	void WriteMonitor(FILE *f, bool first_record, double time);	//Escreve no monitor do SpecialConstraint
	//Variaveis
	int node_A;
	int node_B;
	int csA;//coordinate system - e3 sera a orienta��o em que a rota��o e transmitida
	int csB;//coordinate system - e3 sera a orienta��o em que a rota��o e transmitida

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

	double temp_v[1000];	//Para fun��es do AceGen
	//Calcula contribui��es do residuo e operador tangente - gerado no AceGen
	void EvaluateUniversalJointContribution(double *v, double *residual
		, double **stiffness, double *alphaA, double *alphaB, double *ei1A
		, double *ei2B, double *lambda);
	//Calcula contribui��es do residuo e operador tangente - gerado no AceGen - sem usar SMSD
	void EvaluateUniversalJointContribution2(double *v, double *residual
		, double **stiffness, double *alphaA, double *alphaB, double *ei1A
		, double *ei2B, double *lambda);
};

