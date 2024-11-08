#pragma once
#include "Geometry.h"

class SSContactData;
class Interface_0;
class Interface_1;

using namespace std;

class ContactBodyBody
{
public:
	ContactBodyBody();
	virtual ~ContactBodyBody();

	int nDOF;
	bool previous_evaluation;			//Flag to evaluate previous contact condition
	Geometry* sA;
	Geometry* sB;

	int index1;				//Body 1 - index
	int index2;				//Body 2 - index
	int sub_index1;			//Body 1 - sub_index
	int sub_index2;			//Body 2 - sub_index
	bool invert;			//Inversion of Body types between index1 and index2
	bool cur_active;		//Flag active/unnactive - current status
	bool prev_active;		//Flag active/unnactive - previous status (of current time-step)

	void PreCalc();
	void ProcessSurfacePair();
	//Explicit
	void FinalProcessSurfacePairsExplicit(double t);
	void EvaluateNormalGap();
	bool HaveErrors();
	void SolveLCP();
	void SurfacePoints();					//Sets GammaA and GammaB with contact positions on both surfaces
	void SolvePreviousContact();

	virtual void Alloc() = 0;
	virtual void Free() = 0;
	virtual void MountLocalContributions() = 0;
	virtual void SetVariables() = 0;					//Sets variables for AceGen codes interfaces
	virtual void Report() = 0;
	virtual void CompactReport() = 0;
	virtual void InitialGuess() = 0;
	
	virtual double ObjectivePhase1(Matrix& mc) = 0;								//Calcula a função objetivo para um conjunto de coordenadas convectivas - Phase 1
	virtual void GradientPhase1(Matrix& mc, Matrix& mGra) = 0;					//Calcula o Gradiente da função objetivo - Phase 1
	virtual void HessianPhase1(Matrix& mc, Matrix& mHes) = 0;					//Calcula a Hessiana da função objetivo - Phase 1

	virtual double Gap(Matrix& mc, bool fixed_normals, Matrix& nA, Matrix& nB) = 0;											//Calcula e rotorna o gap (com sinal)
	virtual void GradientGap(Matrix& mc, Matrix& mGra, bool fixed_normals, Matrix& nA, Matrix& nB) = 0;						//Calcula o Gradiente do gap
	virtual void HessianGap(Matrix& mc, Matrix& mHes, bool fixed_normals, Matrix& nA, Matrix& nB) = 0;						//Calcula a Hessiana do gap
	virtual int VerifyConvectiveRange(Matrix& mc) = 0;							//Verifica range de coordenadas convectivas
	virtual void InitializeConvectiveRange() = 0;								//Initialize range of validity of convective coordinates

	bool FindMinimumSolution(SSContactData* c_data, Matrix* solution, int &return_info);								//Otimização - determinação de minimo
	bool FindMinimumSolutionDegenerated(SSContactData* c_data, Matrix* P_0, Matrix* solution);						//Otimização - determinação de minimo
	bool FindSaddleSolution(SSContactData* c_data, Matrix* solution, int &return_info, bool return_gap);				//Otimização - determinação de sela
	bool FindSaddleSolutionDegenerated(SSContactData* c_data, Matrix* P_0, Matrix* solution, bool return_gap);		//Otimização - determinação de sela
	bool FindMinimumGapDegenerated(SSContactData* c_data, Matrix* P_0, Matrix* solution, int &return_info, bool fixed_normals, Matrix& nA, Matrix& nB);			//Otimização - determinação de minimo do gap com sinal
	int CharacterizeCriticalPoint(Matrix* solution);
	int CharacterizeCriticalPointDegenerated(Matrix* solution, Matrix* P_0, bool print = false);
	void DefaultValues();
	double tol_small_1;		//Criterio para numero muito pequeno - residuo != 0
	double tol_eig;			//Criterio para autovalor baixo
	double tol_convective;	//Criterio para maximo erro nas coordenadas convectivas
	double tol_ascent;
	int max_it_1;			//Numero maximo de iterações para otimização - phase 1
	int max_it_2;			//Numero maximo de iterações para otimização - phase 2

	void MountGlobal();
	void MountGlobalExplicit();
	void MountContact();
	void MountContactExplicit(double t);
	void SaveConfiguration();
	void Clear();
	bool NightOwlContact();
	void WriteVTK_XMLForces(FILE *f);

	//Variables for contact evaluation:
	SSContactData* cd;						//information of the contact between surfaces
	double* Rc;								//Vetor de esforços internos
	double** Kc;							//Matriz de rigidez

	void EvaluateInvertedHessian();			//Calcula a inversa da Hessiana
	
	Matrix* I3;
	Matrix* GammaA;
	Matrix* GammaB;

	//Common AceGenPointers (to all surface pairs)
	double* fn;
	double* ft;
	double* v;
	double* dA;
	double* dB;
	double* duiA;
	double* dduiA;
	double* duiB;
	double* dduiB;
	
	double* me;				//Equivalent body-pair mass (for contact-impact purposes)

	bool eligible;			//Set by the function that creates surface pairs - may change in each iteration of NR
	bool prev_eligible;		//Previous converged eligible

	bool alloc_control;

	double*dA_zero;
	double*dB_zero;

	Matrix* convective_range;//Matrix que guarda os ranges de coordenadas convectivas validas para as superficies
	Matrix* convective_max;	//Matrix que guarda os valores maximos de coordenadas convectivas
	Matrix* convective_min;	//Matrix que guarda os valores minimos de coordenadas convectivas
	double minimum_convective_range;

	//TR report
	FILE *f_TR_report;
	void OpenTRReport();
	void InitializeTRReport();
	bool write_report;
	bool write_report_diverged;
	char name[1000];

	//Interface model parameters
	double(*epsn1);
	double(*n1);
	double(*n2);
	double(*gnb);
	double(*gnbb);
	double(*zetan);
	double(*mus);
	double(*mud);
	double(*epst);
	double(*ct);

	bool interface_0_flag;
	bool interface_1_flag;
	Interface_0* inter_0;
	Interface_1* inter_1;
};

