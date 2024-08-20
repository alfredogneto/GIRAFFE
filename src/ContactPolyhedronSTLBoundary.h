#pragma once
#include "ContactParticleBoundary.h"
#include "Polyhedron.h"
#include "STLBoundary.h"
#include "STLSurface.h"
#include "TriangularFace.h"

class ContactPolyhedronSTLBoundary :
	public ContactParticleBoundary
{
public:
	ContactPolyhedronSTLBoundary();
	~ContactPolyhedronSTLBoundary();

	void PreCalc();
	void MountGlobal();
	void MountGlobalExplicit();
	void ProcessContactHierarchy();

	Polyhedron* pA;
	STLBoundary* pB;

	STLSurface* surfA;
	STLSurface* surfB;

	TriangularFace *faceA;
	TriangularFace *faceB;

	Matrix* QAp;
	Matrix* QBp;
	Matrix* x0Ap;
	Matrix* x0Bp;

	Matrix* QAi;
	Matrix* QBi;
	Matrix* x0Ai;
	Matrix* x0Bi;

	Matrix* Q0A;
	Matrix* Q0B;

	//Degenerations
	int deg_pointA;
	int deg_curveA;
	int deg_pointB;
	int deg_curveB;

	void InsertNewContact(int deg_pointA, int deg_curveA, int deg_pointB, int deg_curveB, int faceAID, int faceBID);
	int CreateSurfacePair(int deg_pointA, int deg_curveA, int deg_pointB, int deg_curveB, int faceAID, int faceBID);				//Returns the list index where the contact was created
	void ProcessSurfacePair(int list_index);
	void ProcessSurfacePairs();
	//Explicit
	void FinalProcessSurfacePairsExplicit(double t);
	void FinalProcessSurfacePairExplicit(int list_index, double t);
};

