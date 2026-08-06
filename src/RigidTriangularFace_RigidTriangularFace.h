#pragma once
#include "SurfacePairGeneralContact.h"


class Polyhedron;
class STLSurface;
class TriangularFace;
class Interface_1;

class RigidTriangularFace_RigidTriangularFace :
	public SurfacePairGeneralContact
{
public:
	RigidTriangularFace_RigidTriangularFace();
	~RigidTriangularFace_RigidTriangularFace();

	void EvaluateNormalGap();
	bool HaveErrors();					//Checks for possible errors
	void SolveLCP();					//Solves LCP for the surface pair
	void SolvePreviousContact();		//Solves the previous contact problem
	void MountLocalContributions();		//Mounts local contributions due to contact occurrence
	void PreCalc();						//PreCalcs vertices for parameterizations and contact degenerations
	void SetVariables();				//Sets variables for AceGen codes interfaces
	void HessianPhase1(Matrix& mHes);
	void SurfacePoints();
	void Report();
	void CompactReport();
	void PredictorTimeStep(double kin);

	void MountLocalContributionsExplicit(double t);
	void SetVariablesExplicit(double t);
	void FinalUpdateExplicit(double t);

	void AllocSpecific();
	void FreeSpecific();
	bool alloc_specific_control;

	//IDs of vertices of STL surfaces to be considered to triangular surface parameterization (following the defined sequence)
	int vertexIDsA[3];
	int vertexIDsB[3];

	bool previous_evaluation;			//Flag to evaluate previous contact condition
	double extra_leng;					//Extra length to consider in parameterization parameters

	//Node ID's
	int node_A;
	int node_B;
	//CAD ID's
	int CAD_AID;
	int CAD_BID;
	//Face ID's
	int faceAID;
	int faceBID;
	//Material ID's (for interface law atribution)
	int material_A;
	int material_B;

	//Pointers
	STLSurface* surfA;
	STLSurface* surfB;
	TriangularFace *faceA;
	TriangularFace *faceB;
	Matrix* nA;
	Matrix* nB;

	double* x1A;
	double* x2A;
	double* x3A;
	double* x1B;
	double* x2B;
	double* x3B;
	double* gti;
	double* gtpupdated;
	bool *stick;
	bool *stickupdated;
	double** invH;

	void PrintAceGenPointers();

	//Pointers
	Matrix* ptrQAp;
	Matrix* ptrQBp;
	Matrix* ptrx0Ap;
	Matrix* ptrx0Bp;

	Matrix* ptrQAi;
	Matrix* ptrQBi;
	Matrix* ptrx0Ai;
	Matrix* ptrx0Bi;

	Matrix* ptrQ0A;
	Matrix* ptrQ0B;

	//Matrix* I3;

	//AceGen variables - alloced in AllocSpecific only on demand
	double* cAp;
	double* cBp;
	double* cAi;
	double* cBi;
	double* xAi;
	double* xBi;
	double* uA;
	double* uB;
	double* alphaA;
	double* alphaB;
	double* duiA;
	double* duiB;
	double* dalphaiA;
	double* dalphaiB;
	double* dduiA;
	double* dduiB;
	double* ddalphaiA;
	double* ddalphaiB;
	double** QAi;
	double** QBi;

	//Explicit
	double** Q0A;
	double** Q0B;
	
	//Contact model parameters
	double(*epsn1);
	double(*n1);
	double(*n2);
	double(*gnb);
	double(*gnbb);
	double(*zetan);
	double(*meq);
	double(*mus);
	double(*mud);
	double(*epst);
	double(*ct);
	Interface_1* inter;

	//Interface work variables
	double* Wnum;
	double* Wteo;
};


