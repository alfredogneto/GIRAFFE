#pragma once
#include <stdio.h>

#include "Matrix.h"

class SPContactData;

class SplineElementPair
{
public:
	SplineElementPair();
	virtual ~SplineElementPair();
	//Chute inicial para coordenadas convectivas do par de superf�cies
	virtual void InitialGuess(SPContactData* c_data) = 0;
	//Single function for both contributions
	virtual void ContactSS(bool *stick, bool *stickupdated, bool *previouscontact, double* Rc, double** Kc, double** invH, double* convective, double* copy_convective, double* gti, double* gtpupdated, double* epsn, double* epsn_n, double* epst, double* cn, double* ct, double* mus, double* mud, double* fn, double* ft) = 0;
	
	virtual double ObjectivePhase1(Matrix& mc) = 0;								//Calcula a fun��o objetivo para um conjunto de coordenadas convectivas - Phase 1
	virtual void GradientPhase1(Matrix& mc, Matrix& mGra) = 0;					//Calcula o Gradiente da fun��o objetivo - Phase 1
	virtual void HessianPhase1(Matrix& mc, Matrix& mHes) = 0;					//Calcula a Hessiana da fun��o objetivo - Phase 1

	virtual int VerifyConvectiveRange(Matrix& mc) = 0;							//Verifica range de coordenadas convectivas
	virtual void InitializeConvectiveRange() = 0;								//Initialize range of validity of convective coordinates

	void WriteConvectiveRange();						//Escreve convective range no report
	void PreCalc();																										//Pr� c�lculo
	void Alloc(SPContactData* c_data);																				//Aloca mem�ria
	void Free();																										//Libera mem�ria											
	void DefaultValues();																								//Valores Default de toler�ncias e outras vari�veis

	void EvaluateInvertedHessian(SPContactData* c_data);																//Calcula a inversa da Hessiana
	bool FindMinimumSolution(SPContactData* c_data, Matrix* solution, int &return_info);								//Otimiza��o - determina��o de m�nimo
	bool FindMinimumSolutionDegenerated(SPContactData* c_data, Matrix* P_0, Matrix* solution);							//Otimiza��o - determina��o de m�nimo

	void BeginStepCheck(SPContactData* c_data);
	bool EndStepCheck(SPContactData* c_data);
	void SolveLCP(SPContactData* c_data);																				//Resolve problema local de contato

	int CharacterizeCriticalPoint(Matrix* solution);
	int CharacterizeCriticalPointDegenerated(Matrix* solution, Matrix* P_0, bool print = false);

	//Vari�veis internas
	int spline1_ID;			//ID da spline 1
	int spline2_ID;			//ID da spline 2
	int surf1_ID;			//ID do spline element 1
	int surf2_ID;			//ID do spline element 2
	
	bool alloc_control;
	Matrix** cNR1;			//Solu��o da phase 1

	double tol_small_1;		//Crit�rio para n�mero muito pequeno - res�duo != 0
	double tol_eig;			//Crit�rio para autovalor baixo
	double tol_convective;	//Criterio para maximo erro nas coordenadas convectivas
	double tol_ascent;
	int seq_number;

	int max_it_1;			//N�mero m�ximo de itera��es para otimiza��o - phase 1

	int n_pointwise;		//N�mero de intera��es pointwise (atribu�do no construtor)
		
	Matrix convective_range;//Matrix que guarda os ranges de coordenadas convectivas validas para as superficies
	Matrix convective_max;	//Matrix que guarda os valores m�ximos de coordenadas convectivas
	Matrix convective_min;	//Matrix que guarda os valores m�nimos de coordenadas convectivas
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

