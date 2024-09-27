#pragma once
#include "SurfacePairGeneralContact.h"

class STLSurface;
class TriangularFace;
class Interface_1;
class Matrix;

class FlexibleTriangularFace_RigidTriangularFace :
	public SurfacePairGeneralContact
{
public:
	FlexibleTriangularFace_RigidTriangularFace();
	~FlexibleTriangularFace_RigidTriangularFace();

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

	void AllocSpecificExplicit();
	void FreeSpecificExplicit();
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

	//SuperNode ID
	int super_node_A;
	//Node ID
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

	double* gti;
	double* gtpupdated;
	bool *stick;
	bool *stickupdated;
	double** invH;

	//Pointers
	Matrix** vertices_i_A;
	Matrix** vertices_p_A;
	Matrix* ptrQBp;
	Matrix* ptrx0Bp;
	Matrix* ptrQBi;
	Matrix* ptrx0Bi;

	double* x1B;
	double* x2B;
	double* x3B;

	//AceGen variables - alloced in AllocSpecific only on demand
	double* cAp;
	double* cBp;
	double* cAi;
	double* cBi;
	double* xi1A;
	double* xi2A;
	double* xi3A;
	double* u1A;
	double* u2A;
	double* u3A;
	double* dui1A;
	double* dui2A;
	double* dui3A;
	double* ddui1A;
	double* ddui2A;
	double* ddui3A;

	
	double* xBi;
	double* uB;
	double* alphaB;
	double* duiB;
	double* dalphaiB;
	double* dduiB;
	double* ddalphaiB;
	double** QBi;
	

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
};

