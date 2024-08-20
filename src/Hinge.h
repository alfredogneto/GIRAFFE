#pragma once
#include "SpecialConstraint.h"
class Hinge :
	public SpecialConstraint
{
public:
	Hinge();
	~Hinge();

	bool Read(FILE *f);		//Leitura
	void Write(FILE *f);	//Grava��o
	void Mount();			//Montagem dos res�duos e rigidez tangente
	void MountGlobal();		//Preenche a contribui��o do elemento nas matrizes globais
	void PreCalc();			//Pr�-c�lculo de vari�veis que � feito uma �nica vez no in�cio
	void SaveLagrange();	//Salvando vari�veis da configura��o convergida
	void ActivateDOFs();	//Checa quais multiplicadores de lagrange ser�o ativados,de acordo com a ativa��o dos GLs dos n�s dos quais a special constraint participa
	bool Check();			//Checa inconsist�ncias no SC para evitar erros de execu��o
	void ClearContributions();			//Zera matrizes e vetores
	void ComputeVelAccel();		//Computa efeito das condi��es iniciais nos n�s da restri��o
	void ComputeInitialGuessDisplacements();
	void WriteVTK_XMLRender(FILE *f);
	void WriteMonitor(FILE *f, bool first_record, double time);	//Escreve no monitor do SpecialConstraint
	//Vari�veis
	int node_A;
	int node_B;
	int cs;//coordinate system - e3 ser� a orienta��o da articula��o

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
	Matrix ei3A;
	Matrix ei1B;
	Matrix ei2B;
	Matrix ei1A;
	Matrix ei2A;

	//Matriz de rigidez tangente e vetor res�duo
	Matrix* stiffness1;
	Matrix* residual1;
	double** stiffness2;
	double* residual2;
	double* temp_lambda;

	//Mola de tor��o
	double** stiffness_spring;
	double* residual_spring;
	double thetai, thetad, stiffc;

	//Amortecedor de tor��o
	double** stiffness_damper;
	double* residual_damper;
	double dampc1, dampc2;
	Matrix omegaiA;
	Matrix omegaiB;
	Matrix domegaiA;
	Matrix domegaiB;

	double temp_v[10000];	//Para fun��es do AceGen
	//Calcula contribui��es do res�duo e operador tangente - gerado no AceGen
	void EvaluateHingeContribution(double *v, double *residual
		, double **stiffness, double *alphaA, double *alphaB, double *ei3A
		, double *ei1B, double *ei2B, double *lambda);
	//Calcula contribui��es do res�duo e operador tangente - gerado no AceGen (sem usar SMSD)
	void EvaluateHingeContribution2(double *v, double *residual
		, double **stiffness, double *alphaA, double *alphaB, double *ei3A
		, double *ei1B, double *ei2B, double *lambda);
	//Calcula contribui��es da mola de tor��o (spring)
	void EvaluateTorsionSpring(double *v, double *residual
		, double **stiffness, double *alphaA, double *alphaB, double *ei1A, double *ei1B, double *ei3A
		, double(*stiffc), double(*thetai)
		, double(*thetad));
	//Calcula contribui��es do amortecedor de tor��o
	void EvaluateTorsionDamping(double *v, double *residual
		, double **stiffness, double *alphaA, double *alphaB, double *omegaiA
		, double *omegaiB, double *domegaiA, double *domegaiB, double *ei3A, double
		(*dampc1), double(*dampc2), double(*alpha4), double(*alpha5), double(*alpha6));
};

