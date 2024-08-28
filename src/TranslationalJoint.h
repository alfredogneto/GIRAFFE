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
	void Write(FILE *f);	//Gravação
	void Mount();			//Montagem dos residuos e rigidez tangente
	void MountGlobal();		//Preenche a contribuição do elemento nas matrizes globais
	void PreCalc();			//Pre-calculo de variaveis que e feito uma unica vez no inicio
	void SaveLagrange();	//Salvando variaveis da configuração convergida
	void ActivateDOFs();	//Checa quais multiplicadores de lagrange serão ativados,de acordo com a ativação dos GLs dos nós dos quais a special constraint participa
	bool Check();			//Checa inconsistências no SC para evitar erros de execução
	void ClearContributions();			//Zera matrizes e vetores
	void ComputeVelAccel();					//Computa efeito das condições iniciais nos nós da restrição
	void ComputeInitialGuessDisplacements();
	void WriteVTK_XMLRender(FILE *f);
	void WriteMonitor(FILE *f, bool first_record, double time);	//Escreve no monitor do SpecialConstraint
	//Variaveis
	int node_A;
	int node_B;
	int rot_node;	
	int cs;					//coordinate system - e3 sera a orientação com translação relativa permitida entre A e B

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

	double temp_v[1000];	//Para funções do AceGen
	//Calcula contribuições do residuo e operador tangente - gerado no AceGen
	void EvaluateTranslationalContribution(double *v, double *residual, double **stiffness, double *uA, double *uB, double *alphaA, double *ei1A, double *ei2A, double *lambda);
};

