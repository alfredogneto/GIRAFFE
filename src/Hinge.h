#pragma once
#include "SpecialConstraint.h"
class Hinge :
	public SpecialConstraint
{
public:
	Hinge();
	~Hinge();

	bool Read(FILE *f);		//Leitura
	void Write(FILE *f);	//Gravação
	void Mount();			//Montagem dos resíduos e rigidez tangente
	void MountGlobal();		//Preenche a contribuição do elemento nas matrizes globais
	void PreCalc();			//Pré-cálculo de variáveis que é feito uma única vez no início
	void SaveLagrange();	//Salvando variáveis da configuração convergida
	void ActivateDOFs();	//Checa quais multiplicadores de lagrange serão ativados,de acordo com a ativação dos GLs dos nós dos quais a special constraint participa
	bool Check();			//Checa inconsistências no SC para evitar erros de execução
	void ClearContributions();			//Zera matrizes e vetores
	void ComputeVelAccel();		//Computa efeito das condições iniciais nos nós da restrição
	void ComputeInitialGuessDisplacements();
	void WriteVTK_XMLRender(FILE *f);
	void WriteMonitor(FILE *f, bool first_record, double time);	//Escreve no monitor do SpecialConstraint
	//Variáveis
	int node_A;
	int node_B;
	int cs;//coordinate system - e3 será a orientação da articulação

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

	//Matriz de rigidez tangente e vetor resíduo
	Matrix* stiffness1;
	Matrix* residual1;
	double** stiffness2;
	double* residual2;
	double* temp_lambda;

	//Mola de torção
	double** stiffness_spring;
	double* residual_spring;
	double thetai, thetad, stiffc;

	//Amortecedor de torção
	double** stiffness_damper;
	double* residual_damper;
	double dampc1, dampc2;
	Matrix omegaiA;
	Matrix omegaiB;
	Matrix domegaiA;
	Matrix domegaiB;

	double temp_v[10000];	//Para funções do AceGen
	//Calcula contribuições do resíduo e operador tangente - gerado no AceGen
	void EvaluateHingeContribution(double *v, double *residual
		, double **stiffness, double *alphaA, double *alphaB, double *ei3A
		, double *ei1B, double *ei2B, double *lambda);
	//Calcula contribuições do resíduo e operador tangente - gerado no AceGen (sem usar SMSD)
	void EvaluateHingeContribution2(double *v, double *residual
		, double **stiffness, double *alphaA, double *alphaB, double *ei3A
		, double *ei1B, double *ei2B, double *lambda);
	//Calcula contribuições da mola de torção (spring)
	void EvaluateTorsionSpring(double *v, double *residual
		, double **stiffness, double *alphaA, double *alphaB, double *ei1A, double *ei1B, double *ei3A
		, double(*stiffc), double(*thetai)
		, double(*thetad));
	//Calcula contribuições do amortecedor de torção
	void EvaluateTorsionDamping(double *v, double *residual
		, double **stiffness, double *alphaA, double *alphaB, double *omegaiA
		, double *omegaiB, double *domegaiA, double *domegaiB, double *ei3A, double
		(*dampc1), double(*dampc2), double(*alpha4), double(*alpha5), double(*alpha6));
};

