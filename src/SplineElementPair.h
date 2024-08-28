#pragma once
#include <stdio.h>

#include "Matrix.h"

class SPContactData;

class SplineElementPair
{
public:
	SplineElementPair();
	virtual ~SplineElementPair();
	//Chute inicial para coordenadas convectivas do par de superfícies
	virtual void InitialGuess(SPContactData* c_data) = 0;
	//Single function for both contributions
	virtual void ContactSS(bool *stick, bool *stickupdated, bool *previouscontact, double* Rc, double** Kc, double** invH, double* convective, double* copy_convective, double* gti, double* gtpupdated, double* epsn, double* epsn_n, double* epst, double* cn, double* ct, double* mus, double* mud, double* fn, double* ft) = 0;
	
	virtual double ObjectivePhase1(Matrix& mc) = 0;								//Calcula a função objetivo para um conjunto de coordenadas convectivas - Phase 1
	virtual void GradientPhase1(Matrix& mc, Matrix& mGra) = 0;					//Calcula o Gradiente da função objetivo - Phase 1
	virtual void HessianPhase1(Matrix& mc, Matrix& mHes) = 0;					//Calcula a Hessiana da função objetivo - Phase 1

	virtual int VerifyConvectiveRange(Matrix& mc) = 0;							//Verifica range de coordenadas convectivas
	virtual void InitializeConvectiveRange() = 0;								//Initialize range of validity of convective coordinates

	void WriteConvectiveRange();						//Escreve convective range no report
	void PreCalc();																										//Pré cálculo
	void Alloc(SPContactData* c_data);																				//Aloca memória
	void Free();																										//Libera memória											
	void DefaultValues();																								//Valores Default de tolerâncias e outras variáveis

	void EvaluateInvertedHessian(SPContactData* c_data);																//Calcula a inversa da Hessiana
	bool FindMinimumSolution(SPContactData* c_data, Matrix* solution, int &return_info);								//Otimização - determinação de mínimo
	bool FindMinimumSolutionDegenerated(SPContactData* c_data, Matrix* P_0, Matrix* solution);							//Otimização - determinação de mínimo

	void BeginStepCheck(SPContactData* c_data);
	bool EndStepCheck(SPContactData* c_data);
	void SolveLCP(SPContactData* c_data);																				//Resolve problema local de contato

	int CharacterizeCriticalPoint(Matrix* solution);
	int CharacterizeCriticalPointDegenerated(Matrix* solution, Matrix* P_0, bool print = false);

	//Variáveis internas
	int spline1_ID;			//ID da spline 1
	int spline2_ID;			//ID da spline 2
	int surf1_ID;			//ID do spline element 1
	int surf2_ID;			//ID do spline element 2
	
	bool alloc_control;
	Matrix** cNR1;			//Solução da phase 1

	double tol_small_1;		//Critério para número muito pequeno - resíduo != 0
	double tol_eig;			//Critério para autovalor baixo
	double tol_convective;	//Criterio para maximo erro nas coordenadas convectivas
	double tol_ascent;
	int seq_number;

	int max_it_1;			//Número máximo de iterações para otimização - phase 1

	int n_pointwise;		//Número de interações pointwise (atribuído no construtor)
		
	Matrix convective_range;//Matrix que guarda os ranges de coordenadas convectivas validas para as superficies
	Matrix convective_max;	//Matrix que guarda os valores máximos de coordenadas convectivas
	Matrix convective_min;	//Matrix que guarda os valores mínimos de coordenadas convectivas
	double minimum_convective_range;
	//TR report
	FILE **f_TR_report;
	FILE **f_DEG_report;
	void OpenTRReport(int index);
	void OpenDEGReport(int index);
	void InitializeTRReport(int index);
	void InitializeDEGReport(int index);
	bool write_report;
	bool write_report_diverged;

	char name[1000];
	char name_deg[1000];
};

