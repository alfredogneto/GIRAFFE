#pragma once
#include <stdio.h>
#include "Matrix.h"

class SSContactData;

class SurfacePair
{
public:
	SurfacePair();
	virtual ~SurfacePair();
	//Chute inicial para coordenadas convectivas do par de superfícies
	virtual void InitialGuess(SSContactData* c_data) = 0;				
	//Single function for both contributions
	virtual void ContactSS(bool *stick, bool *stickupdated, bool *previouscontact, double* Rc, double** Kc, double** invH, double* convective, double* copy_convective, double* gti, double* gtpupdated, double* epsn, double* epsn0, double* epst, double* cn, double* ct, double* mus, double* mud, double* fn, double* ft) = 0;

	virtual double ObjectivePhase1(Matrix& mc) = 0;								//Calcula a função objetivo para um conjunto de coordenadas convectivas - Phase 1
	virtual void GradientPhase1(Matrix& mc, Matrix& mGra) = 0;					//Calcula o Gradiente da função objetivo - Phase 1
	virtual void HessianPhase1(Matrix& mc, Matrix& mHes) = 0;					//Calcula a Hessiana da função objetivo - Phase 1
	
	virtual double Gap(Matrix& mc, bool fixed_normals, Matrix& nA, Matrix& nB) = 0;											//Calcula e rotorna o gap (com sinal)
	virtual void GradientGap(Matrix& mc, Matrix& mGra, bool fixed_normals, Matrix& nA, Matrix& nB) = 0;						//Calcula o Gradiente do gap
	virtual void HessianGap(Matrix& mc, Matrix& mHes, bool fixed_normals, Matrix& nA, Matrix& nB) = 0;						//Calcula a Hessiana do gap
	virtual int VerifyConvectiveRange(Matrix& mc) = 0;							//Verifica range de coordenadas convectivas
	virtual void InitializeConvectiveRange() = 0;								//Initialize range of validity of convective coordinates
	
	void WriteConvectiveRange();						//Escreve convective range no report
	void PreCalc();																										//Pré cálculo
	void Alloc(SSContactData* c_data);																				//Aloca memória
	void Free();																										//Libera memória											
	void DefaultValues();																								//Valores Default de tolerâncias e outras variáveis
	
	void EvaluateInvertedHessian(SSContactData* c_data);																//Calcula a inversa da Hessiana
	bool FindMinimumSolution(SSContactData* c_data, Matrix* solution, int &return_info);								//Otimização - determinação de mínimo
	bool FindMinimumSolutionDegenerated(SSContactData* c_data, Matrix* P_0, Matrix* solution);						//Otimização - determinação de mínimo
	bool FindSaddleSolution(SSContactData* c_data, Matrix* solution, int &return_info, bool return_gap);				//Otimização - determinação de sela
	bool FindSaddleSolutionDegenerated(SSContactData* c_data, Matrix* P_0, Matrix* solution, bool return_gap);		//Otimização - determinação de sela
	bool FindMinimumGapDegenerated(SSContactData* c_data, Matrix* P_0, Matrix* solution, int &return_info, bool fixed_normals, Matrix& nA, Matrix& nB);			//Otimização - determinação de mínimo do gap com sinal
	
	void BeginStepCheck(SSContactData* c_data);					
	bool EndStepCheck(SSContactData* c_data);
	void SolveLCP(SSContactData* c_data);																				//Resolve problema local de contato
	
	int CharacterizeCriticalPoint(Matrix* solution);
	int CharacterizeCriticalPointDegenerated(Matrix* solution, Matrix* P_0, bool print = false);
	void AutomaticDegenerationProcedure();																				//Performs automatic degeneration according to eigenvalues of Hessian
	
	//Variáveis internas
	int surf1_ID;			//ID da superfície 1
	int surf2_ID;			//ID da superfície 2
	bool inverted;			//Indica que há inversão dos tipos de superfície (não se aplica quando o par é formado por superfícies de mesmo tipo)
	bool alloc_control;
	Matrix** cNR1;			//Solução da phase 1
	Matrix** cNR2;			//Solução da phase 2
	Matrix** cdeg;			//Solução degenerada no início do incremento
	
	double tol_small_1;		//Critério para número muito pequeno - resíduo != 0
	double tol_eig;			//Critério para autovalor baixo
	double tol_convective;	//Criterio para maximo erro nas coordenadas convectivas
	double tol_ascent;
	int seq_number;

	int max_it_1;			//Número máximo de iterações para otimização - phase 1
	int max_it_2;			//Número máximo de iterações para otimização - phase 2

	bool* flag_degenerated;	//Flag que indica ocorrência de degeneração
	int n_pointwise;		//Número de interações pointwise (atribuído no construtor)
	
	double perc;			//Percentual de coordenada convectiva para considerar longe ou perto do range (afeta return value da funcao ConvectiveRange)
	Matrix convective_range;//Matrix que guarda os ranges de coordenadas convectivas validas para as superficies
	Matrix convective_max;	//Matrix que guarda os valores máximos de coordenadas convectivas
	Matrix convective_min;	//Matrix que guarda os valores mínimos de coordenadas convectivas
	double minimum_convective_range;	
	//TR report
	FILE **f_TR_report;
	void OpenTRReport(int index);
	void InitializeTRReport(int index);
	bool write_report;
	bool write_report_diverged;
	bool specialLCP;
	char name[1000];
};

