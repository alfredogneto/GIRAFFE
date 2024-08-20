#pragma once
#include "ContactParticleBoundary.h"
#include "VEMPolyhedron.h"
#include "STLBoundary.h"
#include "STLSurface.h"
#include "TriangularFace.h"

class ContactVEMPolyhedronSTLBoundary :
	public ContactParticleBoundary
{
public:
	ContactVEMPolyhedronSTLBoundary();
	~ContactVEMPolyhedronSTLBoundary();

	void PreCalc();
	void MountGlobal();
	void MountGlobalExplicit();
	void ProcessContactHierarchy();

	VEMPolyhedron* pA;
	STLBoundary* pB;
	STLSurface* surfA;
	STLSurface* surfB;
	TriangularFace *faceA;
	TriangularFace *faceB;

	//Pointers
	Matrix** vertices_i_A;
	Matrix** vertices_p_A;
	Matrix* QBi;
	Matrix* x0Bi;
	Matrix* QBp;
	Matrix* x0Bp;

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

