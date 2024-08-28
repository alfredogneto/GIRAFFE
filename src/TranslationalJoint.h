#pragma once
#include "SpecialConstraint.h"
#include "Matrix.h"

class TranslationalJoint :
	public SpecialConstraint
{
public:
	TranslationalJoint();
	~TranslationalJoint();

	bool Read(FILE *f);		//Leitura
	void Write(FILE *f);	//Grava��o
	void Mount();			//Montagem dos res�duos e rigidez tangente
	void MountGlobal();		//Preenche a contribui��o do elemento nas matrizes globais
	void PreCalc();			//Pr�-c�lculo de vari�veis que � feito uma �nica vez no in�cio
	void SaveLagrange();	//Salvando vari�veis da configura��o convergida
	void ActivateDOFs();	//Checa quais multiplicadores de lagrange ser�o ativados,de acordo com a ativa��o dos GLs dos n�s dos quais a special constraint participa
	bool Check();			//Checa inconsist�ncias no SC para evitar erros de execu��o
	void ClearContributions();			//Zera matrizes e vetores
	void ComputeVelAccel();					//Computa efeito das condi��es iniciais nos n�s da restri��o
	void ComputeInitialGuessDisplacements();
	void WriteVTK_XMLRender(FILE *f);
	void WriteMonitor(FILE *f, bool first_record, double time);	//Escreve no monitor do SpecialConstraint
	//Vari�veis
	int node_A;
	int node_B;
	int rot_node;	
	int cs;					//coordinate system - e3 ser� a orienta��o com transla��o relativa permitida entre A e B

	Matrix I3;
	
	Matrix alphaA;
	Matrix alphaiA;
	double alpha_escalar_i, g;
	Matrix A;
	Matrix QA;
	Matrix ei1A;
	Matrix ei2A;

	Matrix uA;
	Matrix uB;

	double** stiffness;
	double* residual;

	double temp_v[1000];	//Para fun��es do AceGen
	//Calcula contribui��es do res�duo e operador tangente - gerado no AceGen
	void EvaluateTranslationalContribution(double *v, double *residual, double **stiffness, double *uA, double *uB, double *alphaA, double *ei1A, double *ei2A, double *lambda);
};

