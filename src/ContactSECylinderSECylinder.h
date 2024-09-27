#pragma once
#include "ContactBodyBody.h"

class Matrix;

class ContactSECylinderSECylinder :
	public ContactBodyBody
{
public:
	ContactSECylinderSECylinder();
	~ContactSECylinderSECylinder();

	void Alloc();
	void Free();
	void MountLocalContributions();
	void SetVariables();					//Sets variables for AceGen codes interfaces
	void Report();
	void CompactReport();
	void InitialGuess();

	double ObjectivePhase1(Matrix& mc);								//Calcula a função objetivo para um conjunto de coordenadas convectivas - Phase 1
	void GradientPhase1(Matrix& mc, Matrix& mGra);					//Calcula o Gradiente da função objetivo - Phase 1
	void HessianPhase1(Matrix& mc, Matrix& mHes);					//Calcula a Hessiana da função objetivo - Phase 1
	double Gap(Matrix& mc, bool fixed_normals, Matrix& nA, Matrix& nB);											//Calcula e rotorna o gap (com sinal)
	void GradientGap(Matrix& mc, Matrix& mGra, bool fixed_normals, Matrix& nA, Matrix& nB);						//Calcula o Gradiente do gap
	void HessianGap(Matrix& mc, Matrix& mHes, bool fixed_normals, Matrix& nA, Matrix& nB);						//Calcula a Hessiana do gap
	int VerifyConvectiveRange(Matrix& mc);							//Verifica range de coordenadas convectivas
	void InitializeConvectiveRange();								//Initialize range of validity of convective coordinates

	void PrintAceGenPointers();

	
	double* cAp;
	double* cBp;
	double* cAi;
	double* cBi;

	//AceGen Pointers - specific
	double* aA;
	double* aB;
	double* bA;
	double* bB;
	double* eA;
	double* eB;
	bool* normalintA;
	bool* normalintB;
	double* xAAi;
	double* xBAi;
	double* xABi;
	double* xBBi;
	double** QAAi;
	double** QBAi;
	double** QABi;
	double** QBBi;

	double* gti;
	double* gtpupdated;
	bool *stick;
	bool *stickupdated;
	double** invH;
	bool *interfacelaw0;
};

