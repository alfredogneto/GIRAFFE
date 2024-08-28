#pragma once
#include "ContactParticleBody.h"

class Polyhedron;
class TriangularFace;
class STLSurface;
class Matrix;

class ContactPolyhedronArcExtrusion :
	public ContactParticleBody
{
public:
	ContactPolyhedronArcExtrusion();
	~ContactPolyhedronArcExtrusion();

	void PreCalc();
	void SolvePreviousContact();
	void Alloc();
	void Free();
	void MountLocalContributions();
	void SetVariables();					//Sets variables for AceGen codes interfaces
	void Report();
	void CompactReport();
	void InitialGuess();
	void SurfacePoints();					//Sets GammaA and GammaB with contact positions on both surfaces
	bool HaveErrors();

	double ObjectivePhase1(Matrix& mc);								//Calcula a função objetivo para um conjunto de coordenadas convectivas - Phase 1
	void GradientPhase1(Matrix& mc, Matrix& mGra);					//Calcula o Gradiente da função objetivo - Phase 1
	void HessianPhase1(Matrix& mc, Matrix& mHes);					//Calcula a Hessiana da função objetivo - Phase 1

	//double Gap(Matrix& mc, bool fixed_normals, Matrix& nA, Matrix& nB);											//Calcula e rotorna o gap (com sinal)
	//void GradientGap(Matrix& mc, Matrix& mGra, bool fixed_normals, Matrix& nA, Matrix& nB);						//Calcula o Gradiente do gap
	//void HessianGap(Matrix& mc, Matrix& mHes, bool fixed_normals, Matrix& nA, Matrix& nB);						//Calcula a Hessiana do gap
	int VerifyConvectiveRange(Matrix& mc);							//Verifica range de coordenadas convectivas
	void InitializeConvectiveRange();								//Initialize range of validity of convective coordinates

	void MountGlobal();
	void MountGlobalExplicit();

	void PrintAceGenPointers();

	//Variables - Particle
	Polyhedron* localpA;
	STLSurface* surfA;
	TriangularFace *faceA;
	
	Matrix* QAp;
	Matrix* xAp;
	
	Matrix* matrixQAi;
	double** QAi;//AceGen mirror
	Matrix* matrixxAi;
	double* xAi;


	double* cAp;
	double* cBp;
	double* cAi;
	double* cBi;

	double* x1A;
	double* x2A;
	double* x3A;

	//Degenerations
	int deg_pointA;
	int deg_curveA;


	//IDs of vertices of STL surfaces to be considered to triangular surface parameterization (following the defined sequence)
	int vertexIDsA[3];
	
	//AceGen Pointers - specific
	double* radB;                                           //Raio de curvatura do arco
	double* cpointB;										//Centro de curvatura do arco

	bool* normalintB;
	double* xABi;
	double* xBBi;
	double** QABi;
	double** QBBi;

	double* gti;
	double* gtpupdated;
	bool *stick;
	bool *stickupdated;
	double** invH;
	bool *interfacelaw0;
};

