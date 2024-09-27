#pragma once
#include "SurfacePair.h"

class SSContactData;

class RigidNURBS_1_RigidNURBS_1 :
	public SurfacePair
{
public:
	RigidNURBS_1_RigidNURBS_1();
	~RigidNURBS_1_RigidNURBS_1();

	//Chute inicial para coordenadas convectivas do par de superficies
	void InitialGuess(SSContactData* c_data);
	
	void ContactSS(bool *stick, bool *stickupdated, bool *previouscontact, double* Rc, double** Kc, double** invH, double* convective, double* copy_convective, double* gti, double* gtpupdated, double* epsn, double* epsn0, double* epst, double* cn, double* ct, double* mus, double* mud, double* fn, double* ft);
	double ObjectivePhase1(Matrix& mc);								//Calcula a função objetivo para um conjunto de coordenadas convectivas - Phase 1
	void GradientPhase1(Matrix& mc, Matrix& mGra);					//Calcula o Gradiente da função objetivo - Phase 1
	void HessianPhase1(Matrix& mc, Matrix& mHes);					//Calcula a Hessiana da função objetivo - Phase 1
	
	double Gap(Matrix& mc, bool fixed_normals, Matrix& nA, Matrix& nB);											//Calcula e rotorna o gap (com sinal)
	void GradientGap(Matrix& mc, Matrix& mGra, bool fixed_normals, Matrix& nA, Matrix& nB);						//Calcula o Gradiente do gap
	void HessianGap(Matrix& mc, Matrix& mHes, bool fixed_normals, Matrix& nA, Matrix& nB);						//Calcula a Hessiana do gap
	int VerifyConvectiveRange(Matrix& mc);							//Verifica range de coordenadas convectivas
	void InitializeConvectiveRange();								//Initialize range of validity of convective coordinates

	//Specific for NURBS
	void EvaluateNURBSDerivatives_p(Matrix& mc);						//Evaluates the necessary derivatives to be used as input to AceGen routines
	void EvaluateNURBSDerivatives_i(Matrix& mc);						//Evaluates the necessary derivatives to be used as input to AceGen routines
	void EvaluateNURBSDOFsVariables();									//Evaluates the necessary DOFs variables to be used as input to AceGen routines
	int derivative_order;
	Matrix** dataA;
	Matrix** dataB;
	
	double GAp[3];
	double dGAp[3][2];
	double ddGAp[3][2][2];
	double dddGAp[3][2][2][2];
	double cAp[2];

	double GBp[3];
	double dGBp[3][2];
	double ddGBp[3][2][2];
	double dddGBp[3][2][2][2];
	double cBp[2];

	double GAi[3];
	double dGAi[3][2];
	double cAi[2];

	double GBi[3];
	double dGBi[3][2];
	double cBi[2];

	double xAi[3];
	double uA[3];
	double QAi[3][3];
	double alphaA[3];
	bool* invertnormalA;
	double duiA[3];
	double dduiA[3];
	double dalphaiA[3];
	double ddalphaiA[3];

	double xBi[3];
	double uB[3];
	double QBi[3][3];
	double alphaB[3];
	bool* invertnormalB;
	double duiB[3];
	double dduiB[3];
	double dalphaiB[3];
	double ddalphaiB[3];

	//Last calls of derivatives - improves efficiency
	Matrix last_cp;
	Matrix last_ci;
	bool first_ci;
	bool first_cp;
};

