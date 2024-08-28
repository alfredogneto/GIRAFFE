#include "RigidTriangularFace_RigidTriangularFace.h"

#include "Polyhedron.h"
#include "STLSurface.h"
#include "TriangularFace.h"
#include "Interface_1.h"
#include "SSContactData.h"
#include "ExecutionData.h"
#include "Material.h"
#include "Node.h"
#include "Dynamic.h"
#include "TimeStepControlData.h"


#include"Database.h"
//Variáveis globais
extern
Database db;
#define PI 3.1415926535897932384626433832795

////////////////////////////////////////////////////////////////////
#define _CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include <crtdbg.h>
#ifdef _DEBUG
#define DBG_NEW new ( _NORMAL_BLOCK , __FILE__ , __LINE__ )
// Replace _NORMAL_BLOCK with _CLIENT_BLOCK if you want the
// allocations to be of _CLIENT_BLOCK type
#else
#define DBG_NEW new
#endif
////////////////////////////////////////////////////////////////////

RigidTriangularFace_RigidTriangularFace::RigidTriangularFace_RigidTriangularFace()
{
	index1 = 0;				//index 1
	index2 = 0;				//index 2
	sub_index1 = 0;			//sub index 1
	sub_index2 = 0;			//sub index 2
	active = false;			//default: unnactive

	deg_pointA = 0;	//ID of the point in surface 1 to be considered in the degeneration (zero if no degeneration)
	deg_curveA = 0; //ID of the curve in surface 1 to be considered in the degeneration (zero if no degeneration)
	deg_pointB = 0; //ID of the point in surface 2 to be considered in the degeneration (zero if no degeneration)
	deg_curveB = 0; //ID of the curve in surface 2 to be considered in the degeneration (zero if no degeneration)

	faceAID = 0;	//ID of face A
	faceBID = 0;	//ID of face B
	CAD_AID = 0;	//ID of CAD A
	CAD_BID = 0;	//ID of CAD B
	node_A = 0;		//ID of nó A
	node_B = 0;		//ID of nó B
	material_A = 0;	//ID of material A (for interface law atribution)
	material_B = 0;	//ID of material B (for interface law atribution)

	cd = NULL;
	Rc = NULL;
	Kc = NULL;
	cd = NULL;
	
	vertexIDsA[0] = 0;
	vertexIDsA[1] = 0;
	vertexIDsA[2] = 0;
	vertexIDsB[0] = 0;
	vertexIDsB[1] = 0;
	vertexIDsB[2] = 0;

	surfA = NULL;
	surfB = NULL;
	faceA = NULL;
	faceB = NULL; 

	alloc_specific_control = false;
	

	GammaA = DBG_NEW Matrix(3);
	GammaB = DBG_NEW Matrix(3);

	eligible = false;
	prev_eligible = false;

	previous_evaluation = false;

	extra_leng = 0.0;
	meq = new double;

	Wnum = new double;
	Wteo = new double;
	impactcontrol = new bool;
	*Wnum = 0.0;
	*Wteo = 0.0;
	*impactcontrol = false;
}

RigidTriangularFace_RigidTriangularFace::~RigidTriangularFace_RigidTriangularFace()
{
	SetUnnactive();
	delete GammaA;
	delete GammaB;

	delete meq;
	delete Wnum;
	delete Wteo;
	delete impactcontrol;
}

void RigidTriangularFace_RigidTriangularFace::AllocSpecific()
{
	if (alloc_specific_control == true)
		FreeSpecific();
	if (alloc_specific_control == false)
	{
		cAp = DBG_NEW double[2];
		cBp = DBG_NEW double[2];
		cAi = DBG_NEW double[2];
		cBi = DBG_NEW double[2];
		xAi = DBG_NEW double[3];
		xBi = DBG_NEW double[3];
		uA = DBG_NEW double[3];
		uB = DBG_NEW double[3];
		alphaA = DBG_NEW double[3];
		alphaB = DBG_NEW double[3];
		duiA = DBG_NEW double[3];
		duiB = DBG_NEW double[3];
		dalphaiA = DBG_NEW double[3];
		dalphaiB = DBG_NEW double[3];
		dduiA = DBG_NEW double[3];
		dduiB = DBG_NEW double[3];
		ddalphaiA = DBG_NEW double[3];
		ddalphaiB = DBG_NEW double[3];

		fn = DBG_NEW double[3];
		ft = DBG_NEW double[3];
		vnrel = new double;
		

		QAi = DBG_NEW double*[3];
		for (int i = 0; i < 3; i++)
			QAi[i] = DBG_NEW double[3];
		QBi = DBG_NEW double*[3];
		for (int i = 0; i < 3; i++)
			QBi[i] = DBG_NEW double[3];

		Q0A = DBG_NEW double*[3];
		for (int i = 0; i < 3; i++)
			Q0A[i] = DBG_NEW double[3];
		Q0B = DBG_NEW double*[3];
		for (int i = 0; i < 3; i++)
			Q0B[i] = DBG_NEW double[3];

		Rc = DBG_NEW double[12];
		Kc = DBG_NEW double*[12];
		for (int i = 0; i < 12; i++)
			Kc[i] = DBG_NEW double[12];
		for (int ni = 0; ni < 12; ni++)
			for (int nj = 0; nj < 12; nj++)
				Kc[ni][nj] = 0.0;
		
		alloc_specific_control = true;
	}
	
}

void RigidTriangularFace_RigidTriangularFace::FreeSpecific()
{
	if (alloc_specific_control == true)
	{
		delete[] cAp;
		delete[] cBp;
		delete[] cAi;
		delete[] cBi;
		delete[] xAi;
		delete[] xBi;
		delete[] uA;
		delete[] uB;
		delete[] alphaA;
		delete[] alphaB;
		delete[] duiA;
		delete[] duiB;
		delete[] dalphaiA;
		delete[] dalphaiB;
		delete[] dduiA;
		delete[] dduiB;
		delete[] ddalphaiA;
		delete[] ddalphaiB;

		delete[] fn;
		delete[] ft;
		delete vnrel;

		for (int i = 0; i < 3; i++)
			delete[] QAi[i];
		delete[]QAi;
		for (int i = 0; i < 3; i++)
			delete[] QBi[i];
		delete[]QBi;

		for (int i = 0; i < 3; i++)
			delete[] Q0A[i];
		delete[]Q0A;
		for (int i = 0; i < 3; i++)
			delete[] Q0B[i];
		delete[]Q0B;

		delete[] Rc;
		for (int i = 0; i < 12; i++)
			delete[]Kc[i];
		delete[]Kc;
		 
		alloc_specific_control = false;
	}
}

bool RigidTriangularFace_RigidTriangularFace::HaveErrors()
{
	//cd->Plot();
	if (cd->copy_return_value[0] == 2)
		SolvePreviousContact();
	if (dot(*cd->n[0], *cd->copy_n[0]) < 0.0 )
	{
		if (cd->return_value[0] == 0)
		{
			if (deg_pointA != 0 || deg_pointB != 0)
			{
				//Contact involving a vertex 
				db.myprintf("\nPenetration detected in contact evaluation.");
				if (db.execution_data->print_contact_report)
				{
					//db.myprintf("\nRigidTriangular-RigidTriangular I1: %d I2: %d \nfaceA: %d curveA: %d pointA: %d\nfaceB: %d curveB: %d pointB: %d\n", index1, index2, faceA->ID, deg_curveA, deg_pointA, faceB->ID, deg_curveB, deg_pointB);
					cd->Plot();
					Report();
				}
				return true;
			}
			else
			{
				//Contact not involving a vertex (edge-edge)
				if (cd->copy_return_value[0] == 0)
				{
					db.myprintf("\nPenetration detected in contact evaluation.");
					if (db.execution_data->print_contact_report)
					{
						//db.myprintf("\nRigidTriangular-RigidTriangular I1: %d I2: %d \nfaceA: %d curveA: %d pointA: %d\nfaceB: %d curveB: %d pointB: %d\n", index1, index2, faceA->ID, deg_curveA, deg_pointA, faceB->ID, deg_curveB, deg_pointB);
						cd->Plot();
						Report();
					}
					return true;
				}
				else
					return false;
			}
		}
		else
			return false;
	}
	else
		return false;
}

void RigidTriangularFace_RigidTriangularFace::PreCalc()
{
	
	surfA = static_cast<STLSurface*>(db.cad_data[CAD_AID - 1]);
	surfB = static_cast<STLSurface*>(db.cad_data[CAD_BID - 1]);
	faceA = surfA->faces[faceAID - 1];
	faceB = surfB->faces[faceBID - 1];

	//Case 1 - degeneration of a point in surface A and no degeneration in surface B
	//This point will always be the vertex 1 of the triangular surface A
	if (deg_pointA != 0 && deg_curveA == 0 && deg_pointB == 0 && deg_curveB == 0)
	{
		//Filling info about the vertices IDs - surface A
		if (faceA->verticesIDs[0] == deg_pointA)
		{
			vertexIDsA[0] = faceA->verticesIDs[0];
			vertexIDsA[1] = faceA->verticesIDs[1];
			vertexIDsA[2] = faceA->verticesIDs[2];
		}
		if (faceA->verticesIDs[1] == deg_pointA)
		{
			vertexIDsA[0] = faceA->verticesIDs[1];
			vertexIDsA[1] = faceA->verticesIDs[2];
			vertexIDsA[2] = faceA->verticesIDs[0];
		}
		if (faceA->verticesIDs[2] == deg_pointA)
		{
			vertexIDsA[0] = faceA->verticesIDs[2];
			vertexIDsA[1] = faceA->verticesIDs[0];
			vertexIDsA[2] = faceA->verticesIDs[1];
		}
		//Filling info about the vertices IDs - surface B
		vertexIDsB[0] = faceB->verticesIDs[0];
		vertexIDsB[1] = faceB->verticesIDs[1];
		vertexIDsB[2] = faceB->verticesIDs[2];

		//Degenerated coordinates index (in the canonical basis)
		cd->deg_control[0][0] = true;
		cd->deg_control[0][1] = true;
		cd->deg_control[0][2] = false;
		cd->deg_control[0][3] = false;
	}

	//Case 2 - degeneration of a point in surface B and no degeneration in surface A
	//This point will always be the vertex 1 of the triangular surface B
	if (deg_pointA == 0 && deg_curveA == 0 && deg_pointB != 0 && deg_curveB == 0)
	{
		//Filling info about the vertices IDs - surface B
		if (faceB->verticesIDs[0] == deg_pointB)
		{
			vertexIDsB[0] = faceB->verticesIDs[0];
			vertexIDsB[1] = faceB->verticesIDs[1];
			vertexIDsB[2] = faceB->verticesIDs[2];
		}
		if (faceB->verticesIDs[1] == deg_pointB)
		{
			vertexIDsB[0] = faceB->verticesIDs[1];
			vertexIDsB[1] = faceB->verticesIDs[2];
			vertexIDsB[2] = faceB->verticesIDs[0];
		}
		if (faceB->verticesIDs[2] == deg_pointB)
		{
			vertexIDsB[0] = faceB->verticesIDs[2];
			vertexIDsB[1] = faceB->verticesIDs[0];
			vertexIDsB[2] = faceB->verticesIDs[1];
		}
		//Filling info about the vertices IDs - surface A
		vertexIDsA[0] = faceA->verticesIDs[0];
		vertexIDsA[1] = faceA->verticesIDs[1];
		vertexIDsA[2] = faceA->verticesIDs[2];

		//Degenerated coordinates index (in the canonical basis)
		cd->deg_control[0][0] = false;
		cd->deg_control[0][1] = false;
		cd->deg_control[0][2] = true;
		cd->deg_control[0][3] = true;
	}

	//Case 3 - degeneration of edges in surfaces A and B
	//The edge connecting vertices 1 and 2 will alway be the degenerated one
	if (deg_pointA == 0 && deg_curveA != 0 && deg_pointB == 0 && deg_curveB != 0)
	{
		//Filling info about the vertices IDs - surface A
		int v1A = surfA->edges[deg_curveA - 1].verticesIDs[0];
		int v2A = surfA->edges[deg_curveA - 1].verticesIDs[1];
		if (v1A == faceA->verticesIDs[0])
		{
			if (v2A == faceA->verticesIDs[1])
			{
				vertexIDsA[0] = faceA->verticesIDs[0];
				vertexIDsA[1] = faceA->verticesIDs[1];
				vertexIDsA[2] = faceA->verticesIDs[2];
			}
			else //(v2A == faceA->verticesIDs[2])
			{
				vertexIDsA[0] = faceA->verticesIDs[2];
				vertexIDsA[1] = faceA->verticesIDs[0];
				vertexIDsA[2] = faceA->verticesIDs[1];
			}
		}
		if (v1A == faceA->verticesIDs[1])
		{
			if (v2A == faceA->verticesIDs[2])
			{
				vertexIDsA[0] = faceA->verticesIDs[1];
				vertexIDsA[1] = faceA->verticesIDs[2];
				vertexIDsA[2] = faceA->verticesIDs[0];
			}
			else //(v2A == faceA->verticesIDs[0])
			{
				vertexIDsA[0] = faceA->verticesIDs[0];
				vertexIDsA[1] = faceA->verticesIDs[1];
				vertexIDsA[2] = faceA->verticesIDs[2];
			}
		}
		if (v1A == faceA->verticesIDs[2])
		{
			if (v2A == faceA->verticesIDs[0])
			{
				vertexIDsA[0] = faceA->verticesIDs[2];
				vertexIDsA[1] = faceA->verticesIDs[0];
				vertexIDsA[2] = faceA->verticesIDs[1];
			}
			else //(v2A == faceA->verticesIDs[1])
			{
				vertexIDsA[0] = faceA->verticesIDs[1];
				vertexIDsA[1] = faceA->verticesIDs[2];
				vertexIDsA[2] = faceA->verticesIDs[0];
			}
		}
		//Filling info about the vertices IDs - surface B
		int v1B = surfB->edges[deg_curveB - 1].verticesIDs[0];
		int v2B = surfB->edges[deg_curveB - 1].verticesIDs[1];
		if (v1B == faceB->verticesIDs[0])
		{
			if (v2B == faceB->verticesIDs[1])
			{
				vertexIDsB[0] = faceB->verticesIDs[0];
				vertexIDsB[1] = faceB->verticesIDs[1];
				vertexIDsB[2] = faceB->verticesIDs[2];
			}
			else //(v2B == faceB->verticesIDs[2])
			{
				vertexIDsB[0] = faceB->verticesIDs[2];
				vertexIDsB[1] = faceB->verticesIDs[0];
				vertexIDsB[2] = faceB->verticesIDs[1];
			}
		}
		if (v1B == faceB->verticesIDs[1])
		{
			if (v2B == faceB->verticesIDs[2])
			{
				vertexIDsB[0] = faceB->verticesIDs[1];
				vertexIDsB[1] = faceB->verticesIDs[2];
				vertexIDsB[2] = faceB->verticesIDs[0];
			}
			else //(v2B == faceB->verticesIDs[0])
			{
				vertexIDsB[0] = faceB->verticesIDs[0];
				vertexIDsB[1] = faceB->verticesIDs[1];
				vertexIDsB[2] = faceB->verticesIDs[2];
			}
		}
		if (v1B == faceB->verticesIDs[2])
		{
			if (v2B == faceB->verticesIDs[0])
			{
				vertexIDsB[0] = faceB->verticesIDs[2];
				vertexIDsB[1] = faceB->verticesIDs[0];
				vertexIDsB[2] = faceB->verticesIDs[1];
			}
			else //(v2B == faceB->verticesIDs[1])
			{
				vertexIDsB[0] = faceB->verticesIDs[1];
				vertexIDsB[1] = faceB->verticesIDs[2];
				vertexIDsB[2] = faceB->verticesIDs[0];
			}
		}

		//Degenerated coordinates index (in the canonical basis)
		cd->deg_control[0][0] = false;
		cd->deg_control[0][1] = true;
		cd->deg_control[0][2] = false;
		cd->deg_control[0][3] = true;
	}

	//Case 4 - degeneration of point in surface A and edge on surface B
	//The edge connecting vertices 1 and 2 will alway be the degenerated one
	if (deg_pointA != 0 && deg_curveA == 0 && deg_pointB == 0 && deg_curveB != 0)
	{
		//Filling info about the vertices IDs - surface A
		if (faceA->verticesIDs[0] == deg_pointA)
		{
			vertexIDsA[0] = faceA->verticesIDs[0];
			vertexIDsA[1] = faceA->verticesIDs[1];
			vertexIDsA[2] = faceA->verticesIDs[2];
		}
		if (faceA->verticesIDs[1] == deg_pointA)
		{
			vertexIDsA[0] = faceA->verticesIDs[1];
			vertexIDsA[1] = faceA->verticesIDs[2];
			vertexIDsA[2] = faceA->verticesIDs[0];
		}
		if (faceA->verticesIDs[2] == deg_pointA)
		{
			vertexIDsA[0] = faceA->verticesIDs[2];
			vertexIDsA[1] = faceA->verticesIDs[0];
			vertexIDsA[2] = faceA->verticesIDs[1];
		}

		//Filling info about the vertices IDs - surface B
		int v1B = surfB->edges[deg_curveB - 1].verticesIDs[0];
		int v2B = surfB->edges[deg_curveB - 1].verticesIDs[1];
		if (v1B == faceB->verticesIDs[0])
		{
			if (v2B == faceB->verticesIDs[1])
			{
				vertexIDsB[0] = faceB->verticesIDs[0];
				vertexIDsB[1] = faceB->verticesIDs[1];
				vertexIDsB[2] = faceB->verticesIDs[2];
			}
			else //(v2B == faceB->verticesIDs[2])
			{
				vertexIDsB[0] = faceB->verticesIDs[2];
				vertexIDsB[1] = faceB->verticesIDs[0];
				vertexIDsB[2] = faceB->verticesIDs[1];
			}
		}
		if (v1B == faceB->verticesIDs[1])
		{
			if (v2B == faceB->verticesIDs[2])
			{
				vertexIDsB[0] = faceB->verticesIDs[1];
				vertexIDsB[1] = faceB->verticesIDs[2];
				vertexIDsB[2] = faceB->verticesIDs[0];
			}
			else //(v2B == faceB->verticesIDs[0])
			{
				vertexIDsB[0] = faceB->verticesIDs[0];
				vertexIDsB[1] = faceB->verticesIDs[1];
				vertexIDsB[2] = faceB->verticesIDs[2];
			}
		}
		if (v1B == faceB->verticesIDs[2])
		{
			if (v2B == faceB->verticesIDs[0])
			{
				vertexIDsB[0] = faceB->verticesIDs[2];
				vertexIDsB[1] = faceB->verticesIDs[0];
				vertexIDsB[2] = faceB->verticesIDs[1];
			}
			else //(v2B == faceB->verticesIDs[1])
			{
				vertexIDsB[0] = faceB->verticesIDs[1];
				vertexIDsB[1] = faceB->verticesIDs[2];
				vertexIDsB[2] = faceB->verticesIDs[0];
			}
		}


		//Degenerated coordinates index (in the canonical basis)
		cd->deg_control[0][0] = true;
		cd->deg_control[0][1] = true;
		cd->deg_control[0][2] = false;
		cd->deg_control[0][3] = true;
	}

	//Case 5 - degeneration of point in surface B and edge on surface A
	//The edge connecting vertices 1 and 2 will alway be the degenerated one
	if (deg_pointA == 0 && deg_curveA != 0 && deg_pointB != 0 && deg_curveB == 0)
	{
		//Filling info about the vertices IDs - surface A
		int v1A = surfA->edges[deg_curveA - 1].verticesIDs[0];
		int v2A = surfA->edges[deg_curveA - 1].verticesIDs[1];
		if (v1A == faceA->verticesIDs[0])
		{
			if (v2A == faceA->verticesIDs[1])
			{
				vertexIDsA[0] = faceA->verticesIDs[0];
				vertexIDsA[1] = faceA->verticesIDs[1];
				vertexIDsA[2] = faceA->verticesIDs[2];
			}
			else //(v2A == faceA->verticesIDs[2])
			{
				vertexIDsA[0] = faceA->verticesIDs[2];
				vertexIDsA[1] = faceA->verticesIDs[0];
				vertexIDsA[2] = faceA->verticesIDs[1];
			}
		}
		if (v1A == faceA->verticesIDs[1])
		{
			if (v2A == faceA->verticesIDs[2])
			{
				vertexIDsA[0] = faceA->verticesIDs[1];
				vertexIDsA[1] = faceA->verticesIDs[2];
				vertexIDsA[2] = faceA->verticesIDs[0];
			}
			else //(v2A == faceA->verticesIDs[0])
			{
				vertexIDsA[0] = faceA->verticesIDs[0];
				vertexIDsA[1] = faceA->verticesIDs[1];
				vertexIDsA[2] = faceA->verticesIDs[2];
			}
		}
		if (v1A == faceA->verticesIDs[2])
		{
			if (v2A == faceA->verticesIDs[0])
			{
				vertexIDsA[0] = faceA->verticesIDs[2];
				vertexIDsA[1] = faceA->verticesIDs[0];
				vertexIDsA[2] = faceA->verticesIDs[1];
			}
			else //(v2A == faceA->verticesIDs[1])
			{
				vertexIDsA[0] = faceA->verticesIDs[1];
				vertexIDsA[1] = faceA->verticesIDs[2];
				vertexIDsA[2] = faceA->verticesIDs[0];
			}
		}
		//Filling info about the vertices IDs - surface B
		if (faceB->verticesIDs[0] == deg_pointB)
		{
			vertexIDsB[0] = faceB->verticesIDs[0];
			vertexIDsB[1] = faceB->verticesIDs[1];
			vertexIDsB[2] = faceB->verticesIDs[2];
		}
		if (faceB->verticesIDs[1] == deg_pointB)
		{
			vertexIDsB[0] = faceB->verticesIDs[1];
			vertexIDsB[1] = faceB->verticesIDs[2];
			vertexIDsB[2] = faceB->verticesIDs[0];
		}
		if (faceB->verticesIDs[2] == deg_pointB)
		{
			vertexIDsB[0] = faceB->verticesIDs[2];
			vertexIDsB[1] = faceB->verticesIDs[0];
			vertexIDsB[2] = faceB->verticesIDs[1];
		}

		//Degenerated coordinates index (in the canonical basis)
		cd->deg_control[0][0] = false;
		cd->deg_control[0][1] = true;
		cd->deg_control[0][2] = true;
		cd->deg_control[0][3] = true;
	}
	
	//Case 6 - degeneration of point in surface A and a point on surface B - no convective coordinates are considered in this case (node-node contact)
	if (deg_pointA != 0 && deg_curveA == 0 && deg_pointB != 0 && deg_curveB == 0)
	{
		//Filling info about the vertices IDs - surface A
		if (faceA->verticesIDs[0] == deg_pointA)
		{
			vertexIDsA[0] = faceA->verticesIDs[0];
			vertexIDsA[1] = faceA->verticesIDs[1];
			vertexIDsA[2] = faceA->verticesIDs[2];
		}
		if (faceA->verticesIDs[1] == deg_pointA)
		{
			vertexIDsA[0] = faceA->verticesIDs[1];
			vertexIDsA[1] = faceA->verticesIDs[2];
			vertexIDsA[2] = faceA->verticesIDs[0];
		}
		if (faceA->verticesIDs[2] == deg_pointA)
		{
			vertexIDsA[0] = faceA->verticesIDs[2];
			vertexIDsA[1] = faceA->verticesIDs[0];
			vertexIDsA[2] = faceA->verticesIDs[1];
		}
		//Filling info about the vertices IDs - surface B
		if (faceB->verticesIDs[0] == deg_pointB)
		{
			vertexIDsB[0] = faceB->verticesIDs[0];
			vertexIDsB[1] = faceB->verticesIDs[1];
			vertexIDsB[2] = faceB->verticesIDs[2];
		}
		if (faceB->verticesIDs[1] == deg_pointB)
		{
			vertexIDsB[0] = faceB->verticesIDs[1];
			vertexIDsB[1] = faceB->verticesIDs[2];
			vertexIDsB[2] = faceB->verticesIDs[0];
		}
		if (faceB->verticesIDs[2] == deg_pointB)
		{
			vertexIDsB[0] = faceB->verticesIDs[2];
			vertexIDsB[1] = faceB->verticesIDs[0];
			vertexIDsB[2] = faceB->verticesIDs[1];
		}

		//Degenerated coordinates index (in the canonical basis)
		cd->deg_control[0][0] = true;
		cd->deg_control[0][1] = true;
		cd->deg_control[0][2] = true;
		cd->deg_control[0][3] = true;
	}

	//Case 7 - no degeneration (not treated)
	if (deg_pointA == 0 && deg_curveA == 0 && deg_pointB == 0 && deg_curveB == 0)
	{
		//Filling info about the vertices IDs - surface A
		vertexIDsA[0] = faceA->verticesIDs[0];
		vertexIDsA[1] = faceA->verticesIDs[1];
		vertexIDsA[2] = faceA->verticesIDs[2];

		//Filling info about the vertices IDs - surface B
		vertexIDsB[0] = faceB->verticesIDs[0];
		vertexIDsB[1] = faceB->verticesIDs[1];
		vertexIDsB[2] = faceB->verticesIDs[2];

		//Degenerated coordinates index (in the canonical basis)
		cd->deg_control[0][0] = false;
		cd->deg_control[0][1] = false;
		cd->deg_control[0][2] = false;
		cd->deg_control[0][3] = false;
	}

	//Degeneration basis - canonical basis
	(*cd->P[0])(0, 0) = 1.0;
	(*cd->P[0])(1, 1) = 1.0;
	(*cd->P[0])(2, 2) = 1.0;
	(*cd->P[0])(3, 3) = 1.0;

	//Degenerated coordinates values - in case of degeneration, just fill the desired coordinate value
	cd->copy_deg_coordinates[0][0] = -1.0;
	cd->copy_deg_coordinates[0][1] = -1.0;
	cd->copy_deg_coordinates[0][2] = -1.0;
	cd->copy_deg_coordinates[0][3] = -1.0;

	cd->degenerated[0] = true;
	
	//Fixing convective coordinates
	Matrix temp(4);
	for (int i = 0; i < 4; i++)
		if (cd->deg_control[0][i] == true)
			temp(i, 0) = cd->copy_deg_coordinates[0][i];
	//Transforming into original basis
	temp = (*cd->P[0])*temp;
	for (int i = 0; i < 4; i++)
		cd->convective[0][i] = temp(i, 0);
	//Degenerative operator
	cd->MountDegenerativeOperator();

	x1A = surfA->vertices[vertexIDsA[0] - 1].coord_double->getMatrix();
	x2A = surfA->vertices[vertexIDsA[1] - 1].coord_double->getMatrix();
	x3A = surfA->vertices[vertexIDsA[2] - 1].coord_double->getMatrix();
	x1B = surfB->vertices[vertexIDsB[0] - 1].coord_double->getMatrix();
	x2B = surfB->vertices[vertexIDsB[1] - 1].coord_double->getMatrix();
	x3B = surfB->vertices[vertexIDsB[2] - 1].coord_double->getMatrix();
	gti = cd->copy_g_t[0]->getMatrix();
	gtpupdated = cd->g_t[0]->getMatrix();
	stick = &cd->copy_stick[0];
	stickupdated = &cd->stick[0];
	invH = cd->invHessian[0];

	//Contact Interface law
	bool scape = false;
	for (int i = 0; i < db.number_contactinterfaces && scape == false; i++)
	{
		if ((db.contactinterfaces[i]->material_1 == material_A && db.contactinterfaces[i]->material_2 == material_B) ||
			(db.contactinterfaces[i]->material_1 == material_B && db.contactinterfaces[i]->material_2 == material_A))
		{
			if (typeid(*db.contactinterfaces[i]) == typeid(Interface_1))
			{
				inter = static_cast<Interface_1*>(db.contactinterfaces[i]);
				scape = true;
			}			
		}
	}
	if (scape == false)
		printf("Error in PreCalc of RigidTriangularFace_RigidTriangularFace!\n");

	epsn1 = &inter->epsn1;
	n1 = &inter->n1;
	n2 = &inter->n2;
	gnb = &inter->gnb;
	gnbb = &inter->gnbb;
	zetan = &inter->zetan;
	mus = &inter->mus;
	mud = &inter->mud;
	epst = &inter->epst;
	ct = &inter->ct;

	//Masses (for damping in contact)
	double mA = db.cad_data[CAD_AID - 1]->volume*db.materials[material_A - 1]->rho;
	double mB = db.cad_data[CAD_BID - 1]->volume*db.materials[material_B - 1]->rho;
	double eps = DBL_EPSILON;
	if (mA > eps && mB > eps)
		*meq = (mA*mB) / (mA + mB);
	else
	{
		if (mA > eps)
			*meq = mA;
		else
			*meq = mB;
	}
}

void RigidTriangularFace_RigidTriangularFace::EvaluateNormalGap()
{
	//Cálculo da função gap (escalar)
	SurfacePoints();
	//Gap vetorial
	*cd->g[0] = *GammaA - *GammaB;
	//Normal do contato
	if (norm(*cd->g[0]) != 0.0)
		*cd->n[0] = (1.0 / norm(*cd->g[0]))*(*cd->g[0]);
	else
		zeros(cd->n[0]);
	cd->g_n[0] = +1.0*norm(*cd->g[0]);
	if ((cd->g_n[0] <= *gnb && cd->return_value[0] == 0) /*|| (cd->copy_g_n[0] <= *gnb && cd->copy_return_value[0] == 0)*/)
		eligible = true;
}

void RigidTriangularFace_RigidTriangularFace::SolveLCP()
{
	//Return value:
	//0 - in range
	//4 - out of range
	//1 - problematic return
	double tol_ortho = 1e-12;

	//Case 1 - degeneration of a point in surface A and no degeneration in surface B
	//This point will always be the vertex 1 of the triangular surface A
	if (deg_pointA != 0 && deg_curveA == 0 && deg_pointB == 0 && deg_curveB == 0)
	{
		//Find (zeta,theta) on surface B
		//Surface B
		Matrix x1Bl = *surfB->vertices[vertexIDsB[0] - 1].coord_double;
		Matrix x2Bl = *surfB->vertices[vertexIDsB[1] - 1].coord_double;
		Matrix x3Bl = *surfB->vertices[vertexIDsB[2] - 1].coord_double;
		//Surface A
		Matrix xp = *surfA->vertices[vertexIDsA[0] - 1].coord_double;
		if (previous_evaluation == false)
		{
			x1Bl = *ptrx0Bp + (*ptrQBp) * x1Bl;
			x2Bl = *ptrx0Bp + (*ptrQBp) * x2Bl;
			x3Bl = *ptrx0Bp + (*ptrQBp) * x3Bl;
			xp = *ptrx0Ap + (*ptrQAp) * xp;
		}
		else
		{
			/*x1Bl = *(db.particles[index2]->x0i) + (*(db.particles[index2]->Qi)) * x1Bl;
			x2Bl = *(db.particles[index2]->x0i) + (*(db.particles[index2]->Qi)) * x2Bl;
			x3Bl = *(db.particles[index2]->x0i) + (*(db.particles[index2]->Qi)) * x3Bl;
			xp = *(db.particles[index1]->x0i) + (*(db.particles[index1]->Qi)) * xp;*/

			x1Bl = *ptrx0Bi + (*ptrQBi) * x1Bl;
			x2Bl = *ptrx0Bi + (*ptrQBi) * x2Bl;
			x3Bl = *ptrx0Bi + (*ptrQBi) * x3Bl;
			xp = *ptrx0Ai + (*ptrQAi) * xp;
		}

		Matrix t1 = 0.5 * (x2Bl - x1Bl);
		Matrix t2 = 0.5 * (x3Bl - x1Bl);
		Matrix b = 0.5 * (x2Bl + x3Bl);
		double De = dot(t1, t1) * dot(t2, t2) - dot(t1, t2) * dot(t1, t2);
		double Dz = dot(t2, t2)*dot(xp - b, t1) - dot(t1, t2)*dot(xp - b, t2);
		double Dth = dot(t1, t1)*dot(xp - b, t2) - dot(t1, t2)*dot(xp - b, t1);
		double zB = Dz / De;
		double thB = Dth / De;

		int ret_value;
		if (abs(zB) < (+1.0 + extra_leng) && abs(thB) < (+1.0 + extra_leng) && thB < -(zB - extra_leng))
			ret_value = 0;
		else
			ret_value = 4;

		if (previous_evaluation == false)
		{
			cd->return_value[0] = ret_value;
			cd->convective[0][0] = -1.0;
			cd->convective[0][1] = -1.0;
			cd->convective[0][2] = zB;
			cd->convective[0][3] = thB;
		}
		else
		{
			cd->copy_return_value[0] = ret_value;
			cd->copy_convective[0][0] = -1.0;
			cd->copy_convective[0][1] = -1.0;
			cd->copy_convective[0][2] = zB;
			cd->copy_convective[0][3] = thB;
		}
	}
	
	//Case 2 - degeneration of a point in surface B and no degeneration in surface A
	//This point will always be the vertex 1 of the triangular surface B
	if (deg_pointA == 0 && deg_curveA == 0 && deg_pointB != 0 && deg_curveB == 0)
	{
		//Find (zeta,theta) on surface A

		//Surface A
		Matrix x1Al = *surfA->vertices[vertexIDsA[0] - 1].coord_double;
		Matrix x2Al = *surfA->vertices[vertexIDsA[1] - 1].coord_double;
		Matrix x3Al = *surfA->vertices[vertexIDsA[2] - 1].coord_double;
		//Surface B
		Matrix xp = *surfB->vertices[vertexIDsB[0] - 1].coord_double;

		if (previous_evaluation == false)
		{
			x1Al = *ptrx0Ap + (*ptrQAp) * x1Al;
			x2Al = *ptrx0Ap + (*ptrQAp) * x2Al;
			x3Al = *ptrx0Ap + (*ptrQAp) * x3Al;
			xp = *ptrx0Bp + (*ptrQBp) * xp;
		}
		else
		{
			/*x1Al = *(db.particles[index1]->x0i) + (*(db.particles[index1]->Qi)) * x1Al;
			x2Al = *(db.particles[index1]->x0i) + (*(db.particles[index1]->Qi)) * x2Al;
			x3Al = *(db.particles[index1]->x0i) + (*(db.particles[index1]->Qi)) * x3Al;
			xp = *(db.particles[index2]->x0i) + (*(db.particles[index2]->Qi)) * xp;*/

			x1Al = *ptrx0Ai + (*ptrQAi) * x1Al;
			x2Al = *ptrx0Ai + (*ptrQAi) * x2Al;
			x3Al = *ptrx0Ai + (*ptrQAi) * x3Al;
			xp = *ptrx0Bi + (*ptrQBi) * xp;
		}

		Matrix t1 = 0.5 * (x2Al - x1Al);
		Matrix t2 = 0.5 * (x3Al - x1Al);
		Matrix b = 0.5 * (x2Al + x3Al);
		double De = dot(t1, t1) * dot(t2, t2) - dot(t1, t2) * dot(t1, t2);
		double Dz = dot(t2, t2)*dot(xp - b, t1) - dot(t1, t2)*dot(xp - b, t2);
		double Dth = dot(t1, t1)*dot(xp - b, t2) - dot(t1, t2)*dot(xp - b, t1);
		double zA = Dz / De;
		double thA = Dth / De;

		int ret_value;
		if (abs(zA) < (+1.0 + extra_leng) && abs(thA) < (+1.0 + extra_leng) && thA < -(zA - extra_leng))
			ret_value = 0;
		else
			ret_value = 4;

		if (previous_evaluation == false)
		{
			cd->return_value[0] = ret_value;
			cd->convective[0][0] = zA;
			cd->convective[0][1] = thA;
			cd->convective[0][2] = -1.0;
			cd->convective[0][3] = -1.0;
		}
		else
		{
			cd->copy_return_value[0] = ret_value;
			cd->copy_convective[0][0] = zA;
			cd->copy_convective[0][1] = thA;
			cd->copy_convective[0][2] = -1.0;
			cd->copy_convective[0][3] = -1.0;
		}
	}
	
	//Case 3 - degeneration of edges in surfaces A and B
	//The edge connecting vertices 1 and 2 will alway be the degenerated one
	if (deg_pointA == 0 && deg_curveA != 0 && deg_pointB == 0 && deg_curveB != 0)
	{
		//Minimum distance problem between:
		//edge connecting vertices 1 and 2 of surface A
		//edge connecting vertices 1 and 2 of surface B

		//Find zeta on surface A (theta = -1) and zeta on surface B (theta = -1) -> in both cases the range [-1,1] is within the region of interest

		//Surface A
		Matrix x1Al = *surfA->vertices[vertexIDsA[0] - 1].coord_double;
		Matrix x2Al = *surfA->vertices[vertexIDsA[1] - 1].coord_double;
		//Surface B
		Matrix x1Bl = *surfB->vertices[vertexIDsB[0] - 1].coord_double;
		Matrix x2Bl = *surfB->vertices[vertexIDsB[1] - 1].coord_double;

		if (previous_evaluation == false)
		{
			x1Al = *ptrx0Ap + (*ptrQAp) * x1Al;
			x2Al = *ptrx0Ap + (*ptrQAp) * x2Al;

			x1Bl = *ptrx0Bp + (*ptrQBp) * x1Bl;
			x2Bl = *ptrx0Bp + (*ptrQBp) * x2Bl;
		}
		else
		{
			/*x1Al = *(db.particles[index1]->x0i) + (*(db.particles[index1]->Qi)) * x1Al;
			x2Al = *(db.particles[index1]->x0i) + (*(db.particles[index1]->Qi)) * x2Al;

			x1Bl = *(db.particles[index2]->x0i) + (*(db.particles[index2]->Qi)) * x1Bl;
			x2Bl = *(db.particles[index2]->x0i) + (*(db.particles[index2]->Qi)) * x2Bl;*/

			x1Al = *ptrx0Ai + (*ptrQAi) * x1Al;
			x2Al = *ptrx0Ai + (*ptrQAi) * x2Al;

			x1Bl = *ptrx0Bi + (*ptrQBi) * x1Bl;
			x2Bl = *ptrx0Bi + (*ptrQBi) * x2Bl;
		}

		Matrix bA = 0.5*(x1Al + x2Al);
		Matrix tA = 0.5*(x2Al - x1Al);
		Matrix bB = 0.5*(x1Bl + x2Bl);
		Matrix tB = 0.5*(x2Bl - x1Bl);

		double zeta_A, zeta_B;
		zeta_A = 0.0;
		zeta_B = 0.0;
		int ret_value;
		//Se não houver paralelismo
		if (abs(dot(tA, tA)*dot(tB, tB) - dot(tA, tB)*dot(tA, tB)) > tol_ortho)
		{
			// baseado em Wriggers e Zavarise, 1997
			zeta_A = dot(bA - bB, (1.0 / (dot(tA, tA)*dot(tB, tB) - dot(tA, tB)*dot(tA, tB)))*(tB*dot(tA, tB) - tA * dot(tB, tB)));
			zeta_B = -1.0*dot(bA - bB, (1.0 / (dot(tA, tA)*dot(tB, tB) - dot(tA, tB)*dot(tA, tB)))*(tA*dot(tA, tB) - tB * dot(tA, tA)));

			if ((-1.0 - extra_leng) < zeta_A && zeta_A < (+1.0 + extra_leng) && (-1.0 - extra_leng) < zeta_B && zeta_B < (+1.0 + extra_leng))
				ret_value = 0;
			else
				ret_value = 4;
		}
		else
			ret_value = 4;

		if (previous_evaluation == false)
		{
			cd->return_value[0] = ret_value;
			cd->convective[0][0] = zeta_A;
			cd->convective[0][1] = -1.0;
			cd->convective[0][2] = zeta_B;
			cd->convective[0][3] = -1.0;
		}
		else
		{
			cd->copy_return_value[0] = ret_value;
			cd->copy_convective[0][0] = zeta_A;
			cd->copy_convective[0][1] = -1.0;
			cd->copy_convective[0][2] = zeta_B;
			cd->copy_convective[0][3] = -1.0;
		}
	}
	
	//Case 4 - degeneration of point in surface A and edge on surface B
	//The edge connecting vertices 1 and 2 will alway be the degenerated one
	if (deg_pointA != 0 && deg_curveA == 0 && deg_pointB == 0 && deg_curveB != 0)
	{
		//Surface A - point
		Matrix x1Al = *surfA->vertices[vertexIDsA[0] - 1].coord_double;
		//Surface B - edge
		Matrix x1Bl = *surfB->vertices[vertexIDsB[0] - 1].coord_double;
		Matrix x2Bl = *surfB->vertices[vertexIDsB[1] - 1].coord_double;

		if (previous_evaluation == false)
		{
			x1Al = *ptrx0Ap + (*ptrQAp) * x1Al;
			
			x1Bl = *ptrx0Bp + (*ptrQBp) * x1Bl;
			x2Bl = *ptrx0Bp + (*ptrQBp) * x2Bl;
		}
		else
		{
			/*x1Al = *(db.particles[index1]->x0i) + (*(db.particles[index1]->Qi)) * x1Al;
			
			x1Bl = *(db.particles[index2]->x0i) + (*(db.particles[index2]->Qi)) * x1Bl;
			x2Bl = *(db.particles[index2]->x0i) + (*(db.particles[index2]->Qi)) * x2Bl;*/

			x1Al = *ptrx0Ai + (*ptrQAi) * x1Al;

			x1Bl = *ptrx0Bi + (*ptrQBi) * x1Bl;
			x2Bl = *ptrx0Bi + (*ptrQBi) * x2Bl;
		}

		Matrix b = 0.5*(x1Bl + x2Bl);
		Matrix t = 0.5*(x2Bl - x1Bl);

		double zeta = (1.0/(dot(t,t))) * dot(x1Al - b, t);
		int ret_value;
		
		if (abs(zeta) < (1.0 + extra_leng))
			ret_value = 0;
		else
			ret_value = 4;

		if (previous_evaluation == false)
		{
			cd->return_value[0] = ret_value;
			cd->convective[0][0] = -1.0;
			cd->convective[0][1] = -1.0;
			cd->convective[0][2] = zeta;
			cd->convective[0][3] = -1.0;
		}
		else
		{
			cd->copy_return_value[0] = ret_value;
			cd->copy_convective[0][0] = -1.0;
			cd->copy_convective[0][1] = -1.0;
			cd->copy_convective[0][2] = zeta;
			cd->copy_convective[0][3] = -1.0;
		}
	}

	//Case 5 - degeneration of point in surface B and edge on surface A
	//The edge connecting vertices 1 and 2 will alway be the degenerated one
	if (deg_pointA == 0 && deg_curveA != 0 && deg_pointB != 0 && deg_curveB == 0)
	{
		//Surface B - point
		Matrix x1Bl = *surfB->vertices[vertexIDsB[0] - 1].coord_double;
		//Surface A - edge
		Matrix x1Al = *surfA->vertices[vertexIDsA[0] - 1].coord_double;
		Matrix x2Al = *surfA->vertices[vertexIDsA[1] - 1].coord_double;

		if (previous_evaluation == false)
		{
			x1Bl = *ptrx0Bp + (*ptrQBp) * x1Bl;

			x1Al = *ptrx0Ap + (*ptrQAp) * x1Al;
			x2Al = *ptrx0Ap + (*ptrQAp) * x2Al;
		}
		else
		{
			/*x1Bl = *(db.particles[index2]->x0i) + (*(db.particles[index2]->Qi)) * x1Bl;

			x1Al = *(db.particles[index1]->x0i) + (*(db.particles[index1]->Qi)) * x1Al;
			x2Al = *(db.particles[index1]->x0i) + (*(db.particles[index1]->Qi)) * x2Al;*/

			x1Bl = *ptrx0Bi + (*ptrQBi) * x1Bl;

			x1Al = *ptrx0Ai + (*ptrQAi) * x1Al;
			x2Al = *ptrx0Ai + (*ptrQAi) * x2Al;
		}

		Matrix b = 0.5*(x1Al + x2Al);
		Matrix t = 0.5*(x2Al - x1Al);

		double zeta = (1.0 / (dot(t, t))) * dot(x1Bl - b, t);
		int ret_value;

		if (abs(zeta) < (1.0 + extra_leng))
			ret_value = 0;
		else
			ret_value = 4;

		if (previous_evaluation == false)
		{
			cd->return_value[0] = ret_value;
			cd->convective[0][0] = zeta;
			cd->convective[0][1] = -1.0;
			cd->convective[0][2] = -1.0;
			cd->convective[0][3] = -1.0;
		}
		else
		{
			cd->copy_return_value[0] = ret_value;
			cd->copy_convective[0][0] = zeta;
			cd->copy_convective[0][1] = -1.0;
			cd->copy_convective[0][2] = -1.0;
			cd->copy_convective[0][3] = -1.0;
		}
	}

	//Case 6 - degeneration of point in surface A and a point on surface B
	if (deg_pointA != 0 && deg_curveA == 0 && deg_pointB != 0 && deg_curveB == 0)
	{
		//No LCP is solved in this case
		if (previous_evaluation == false)
		{
			cd->return_value[0] = 0;
			cd->convective[0][0] = -1.0;
			cd->convective[0][1] = -1.0;
			cd->convective[0][2] = -1.0;
			cd->convective[0][3] = -1.0;
		}
		else
		{
			cd->copy_return_value[0] = 0;
			cd->copy_convective[0][0] = -1.0;
			cd->copy_convective[0][1] = -1.0;
			cd->copy_convective[0][2] = -1.0;
			cd->copy_convective[0][3] = -1.0;
		}
	}

	//Case 7 - no degeneration (not treated)
	if (deg_pointA == 0 && deg_curveA == 0 && deg_pointB == 0 && deg_curveB == 0)
	{
		//No LCP is solved in this case
		if (previous_evaluation == false)
		{
			cd->return_value[0] = 4;		//no contact will be addressed
			cd->convective[0][0] = -1.0;
			cd->convective[0][1] = -1.0;
			cd->convective[0][2] = -1.0;
			cd->convective[0][3] = -1.0;
		}
		else
		{
			cd->copy_return_value[0] = 4;	//no contact will be addressed
			cd->copy_convective[0][0] = -1.0;
			cd->copy_convective[0][1] = -1.0;
			cd->copy_convective[0][2] = -1.0;
			cd->copy_convective[0][3] = -1.0;
		}
	}
}

//Solves LCP for the surface pair
void RigidTriangularFace_RigidTriangularFace::SolvePreviousContact()
{
	previous_evaluation = true;
	//LCP
	SolveLCP();
	SurfacePoints();
	//Gap vetorial
	*cd->copy_g[0] = *GammaA - *GammaB;
	//Normal do contato
	if (norm(*cd->copy_g[0]) != 0.0)
		*cd->copy_n[0] = (1.0 / norm(*cd->copy_g[0]))*(*cd->copy_g[0]);
	else
		zeros(cd->copy_n[0]);
	//Gap escalar
	*cd->copy_g_n = norm(*cd->copy_g[0]);
	previous_evaluation = false;
}

void RigidTriangularFace_RigidTriangularFace::SurfacePoints()
{
	double z;
	double th;
	double N1;
	double N2;
	double N3;
	Matrix localx;
	
	//Gamma A
	if (previous_evaluation == false)
	{
		z = cd->convective[0][0];
		th = cd->convective[0][1];
		N1 = -(th / 2.0) - z / 2.0;
		N2 = 0.5 + z / 2.0;
		N3 = 0.5 + th / 2.0;
		localx = N1 * (*surfA->vertices[vertexIDsA[0] - 1].coord_double)
			+ N2 * (*surfA->vertices[vertexIDsA[1] - 1].coord_double)
			+ N3 * (*surfA->vertices[vertexIDsA[2] - 1].coord_double);
		*GammaA = *ptrx0Ap + (*ptrQAp) * localx;
	}
	else
	{
		z = cd->copy_convective[0][0];
		th = cd->copy_convective[0][1];
		N1 = -(th / 2.0) - z / 2.0;
		N2 = 0.5 + z / 2.0;
		N3 = 0.5 + th / 2.0;
		localx = N1 * (*surfA->vertices[vertexIDsA[0] - 1].coord_double)
			+ N2 * (*surfA->vertices[vertexIDsA[1] - 1].coord_double)
			+ N3 * (*surfA->vertices[vertexIDsA[2] - 1].coord_double);
		*GammaA = *ptrx0Ai + (*ptrQAi) * localx;
	}
		
	//Gamma B
	if (previous_evaluation == false)
	{
		z = cd->convective[0][2];
		th = cd->convective[0][3];
		N1 = -(th / 2.0) - z / 2.0;
		N2 = 0.5 + z / 2.0;
		N3 = 0.5 + th / 2.0;
		localx = N1 * (*surfB->vertices[vertexIDsB[0] - 1].coord_double)
			+ N2 * (*surfB->vertices[vertexIDsB[1] - 1].coord_double)
			+ N3 * (*surfB->vertices[vertexIDsB[2] - 1].coord_double);
		*GammaB = *ptrx0Bp + (*ptrQBp) * localx;
	}
	else
	{
		z = cd->copy_convective[0][2];
		th = cd->copy_convective[0][3];
		N1 = -(th / 2.0) - z / 2.0;
		N2 = 0.5 + z / 2.0;
		N3 = 0.5 + th / 2.0;
		localx = N1 * (*surfB->vertices[vertexIDsB[0] - 1].coord_double)
			+ N2 * (*surfB->vertices[vertexIDsB[1] - 1].coord_double)
			+ N3 * (*surfB->vertices[vertexIDsB[2] - 1].coord_double);
		*GammaB = *ptrx0Bi + (*ptrQBi) * localx;
	}
		
}

void RigidTriangularFace_RigidTriangularFace::SetVariables()
{
	for (int i = 0; i < 3; i++)
	{
		xAi[i] = db.nodes[node_A - 1]->copy_coordinates[i];
		xBi[i] = db.nodes[node_B - 1]->copy_coordinates[i];
		uA[i] = db.nodes[node_A - 1]->displacements[i];
		uB[i] = db.nodes[node_B - 1]->displacements[i];
		alphaA[i] = db.nodes[node_A - 1]->displacements[i + 3];
		alphaB[i] = db.nodes[node_B - 1]->displacements[i + 3];

		duiA[i] = db.nodes[node_A - 1]->copy_vel[i];
		duiB[i] = db.nodes[node_B - 1]->copy_vel[i];
		dduiA[i] = db.nodes[node_A - 1]->copy_accel[i];
		dduiB[i] = db.nodes[node_B - 1]->copy_accel[i];
		dalphaiA[i] = db.nodes[node_A - 1]->copy_vel[i + 3];
		dalphaiB[i] = db.nodes[node_B - 1]->copy_vel[i + 3];
		ddalphaiA[i] = db.nodes[node_A - 1]->copy_accel[i + 3];
		ddalphaiB[i] = db.nodes[node_B - 1]->copy_accel[i + 3];
		
		for (int j = 0; j < 3; j++)
		{
			QAi[i][j] = (*ptrQAi)(i, j);
			QBi[i][j] = (*ptrQBi)(i, j);
		}
	}
}

void RigidTriangularFace_RigidTriangularFace::MountLocalContributions()
{
	cAp[0] = cd->convective[0][0];
	cAp[1] = cd->convective[0][1];
	cBp[0] = cd->convective[0][2];
	cBp[1] = cd->convective[0][3];

	cAi[0] = cd->copy_convective[0][0];
	cAi[1] = cd->copy_convective[0][1];
	cBi[0] = cd->copy_convective[0][2];
	cBi[1] = cd->copy_convective[0][3];

	v = DBG_NEW double[10000];
	
	double value = 0.0;
	double* a4;
	double* a5;
	double* a6;
	if (typeid(*db.solution[db.current_solution_number - 1]) == typeid(Dynamic))
	{
		Dynamic* ptr_sol = static_cast<Dynamic*>(db.solution[db.current_solution_number - 1]);
		a4 = &ptr_sol->a4;
		a5 = &ptr_sol->a5;
		a6 = &ptr_sol->a6;
	}
	else
	{
		a4 = &value;
		a5 = &value;
		a6 = &value;
	}
	//Zerando matrizes e vetores
	for (int i = 0; i < 12; i++)
	{
		Rc[i] = 0.0;
		for (int j = 0; j < 12; j++)
			Kc[i][j] = 0.0;
	}

	bool *previouscontact = &prev_eligible;
	
	//Avalia contribuições de contato

//
//#pragma region AceGen
//	double v01; double v010; double v011; double v012; double v013; double v014;
//	double v015; double v016; double v017; double v018; double v019; double v02;
//	double v020; double v021; double v022; double v023; double v024; double v025;
//	double v026; double v027; double v03; double v04; double v05; double v06; double v07;
//	double v08; double v09;
//	int i835, i913, i1159, i1583, i2691, i2928, i3769, i3770, i3771, i3772, i3773, i3774, i3837
//		, i3838, i3839, i3843, i3844, i3845, ii3860, ii3861, ii3862, i3896, i3897, i3898, i3899, i3900
//		, i3901, i4056, i4057, i4058, i4062, i4063, i4064, b673, b674, b706, b728, b729, b730, b731
//		, b737, b738, b747, b748, b764, b765, b766, b767, b783, b800, b801, b813, b840, b841, b997
//		, b998, b1164, b1165, b1166, b1167, b1195, b1199, b1200, b1212, b1213, b1298, b1332, b1935
//		, b1980, b2087, b2088, b2113, b2114, b2137, b2138, b2139, b2156, b2215, b2219, b2223, b2232
//		, b2233, b2286, b2311, b2700, b2701, b2705, b2710, b2711, b3263, b3264, b3281, b3285, b3289
//		, b3294, b3295;
//	v[6318] = -1e0;
//	v[6319] = 0e0;
//	v[6320] = 0e0;
//	v[6321] = 0e0;
//	v[6322] = 0e0;
//	v[6323] = 0e0;
//	v[6324] = 1e0;
//	v[6325] = 0e0;
//	v[6326] = 0e0;
//	v[6327] = 0e0;
//	v[6328] = 0e0;
//	v[6329] = 0e0;
//	v[6330] = 0e0;
//	v[6331] = -1e0;
//	v[6332] = 0e0;
//	v[6333] = 0e0;
//	v[6334] = 0e0;
//	v[6335] = 0e0;
//	v[6336] = 0e0;
//	v[6337] = 1e0;
//	v[6338] = 0e0;
//	v[6339] = 0e0;
//	v[6340] = 0e0;
//	v[6341] = 0e0;
//	v[6342] = 0e0;
//	v[6343] = 0e0;
//	v[6344] = -1e0;
//	v[6345] = 0e0;
//	v[6346] = 0e0;
//	v[6347] = 0e0;
//	v[6348] = 0e0;
//	v[6349] = 0e0;
//	v[6350] = 1e0;
//	v[6351] = 0e0;
//	v[6352] = 0e0;
//	v[6353] = 0e0;
//	v[3782] = -((*a4)*(*ct));
//	v[3767] = (*epsn1)*(*n1);
//	v[5734] = 1e0;
//	v[5735] = 0e0;
//	v[5736] = 0e0;
//	v[5737] = 0e0;
//	v[5738] = 0e0;
//	v[5739] = 0e0;
//	v[5740] = -1e0;
//	v[5741] = 0e0;
//	v[5742] = 0e0;
//	v[5743] = 0e0;
//	v[5744] = 0e0;
//	v[5745] = 0e0;
//	v[5746] = 0e0;
//	v[5747] = 1e0;
//	v[5748] = 0e0;
//	v[5749] = 0e0;
//	v[5750] = 0e0;
//	v[5751] = 0e0;
//	v[5752] = 0e0;
//	v[5753] = -1e0;
//	v[5754] = 0e0;
//	v[5755] = 0e0;
//	v[5756] = 0e0;
//	v[5757] = 0e0;
//	v[5758] = 0e0;
//	v[5759] = 0e0;
//	v[5760] = 1e0;
//	v[5761] = 0e0;
//	v[5762] = 0e0;
//	v[5763] = 0e0;
//	v[5764] = 0e0;
//	v[5765] = 0e0;
//	v[5766] = -1e0;
//	v[5767] = 0e0;
//	v[5768] = 0e0;
//	v[5769] = 0e0;
//	v[725] = (*gnb) - (*gnbb);
//	v[740] = -1e0 + (*n1);
//	v[2115] = -1e0 + v[740];
//	v[745] = -1e0 + (*n2);
//	v[2120] = -1e0 + v[745];
//	v[336] = -x1A[0] / 2e0;
//	v[338] = -x1A[1] / 2e0;
//	v[340] = -x1A[2] / 2e0;
//	v[329] = v[336] + x2A[0] / 2e0;
//	v[330] = v[338] + x2A[1] / 2e0;
//	v[331] = v[340] + x2A[2] / 2e0;
//	v[337] = v[336] + x3A[0] / 2e0;
//	v[339] = v[338] + x3A[1] / 2e0;
//	v[341] = v[340] + x3A[2] / 2e0;
//	v[365] = alphaA[0] / 2e0;
//	v[363] = 2e0*alphaA[0];
//	v[156] = Power(alphaA[0], 2);
//	v[366] = 2e0*alphaA[1];
//	v[5686] = 0e0;
//	v[5687] = 0e0;
//	v[5688] = 0e0;
//	v[5689] = v[363];
//	v[5690] = v[366];
//	v[5691] = 0e0;
//	v[5692] = 0e0;
//	v[5693] = 0e0;
//	v[5694] = 0e0;
//	v[5695] = 0e0;
//	v[5696] = 0e0;
//	v[5697] = 0e0;
//	v[5886] = 0e0;
//	v[5887] = 0e0;
//	v[5888] = 0e0;
//	v[5889] = -v[363];
//	v[5890] = -v[366];
//	v[5891] = 0e0;
//	v[5892] = 0e0;
//	v[5893] = 0e0;
//	v[5894] = 0e0;
//	v[5895] = 0e0;
//	v[5896] = 0e0;
//	v[5897] = 0e0;
//	v[364] = alphaA[1] / 2e0;
//	v[5434] = 0e0;
//	v[5435] = 0e0;
//	v[5436] = 0e0;
//	v[5437] = v[364];
//	v[5438] = v[365];
//	v[5439] = 0e0;
//	v[5440] = 0e0;
//	v[5441] = 0e0;
//	v[5442] = 0e0;
//	v[5443] = 0e0;
//	v[5444] = 0e0;
//	v[5445] = 0e0;
//	v[154] = (alphaA[0] * alphaA[1]) / 2e0;
//	v[149] = Power(alphaA[1], 2);
//	v[415] = -v[149] - v[156];
//	v[393] = alphaA[2] + v[154];
//	v[383] = -alphaA[2] + v[154];
//	v[368] = 2e0*alphaA[2];
//	v[5602] = 0e0;
//	v[5603] = 0e0;
//	v[5604] = 0e0;
//	v[5605] = v[363];
//	v[5606] = 0e0;
//	v[5607] = v[368];
//	v[5608] = 0e0;
//	v[5609] = 0e0;
//	v[5610] = 0e0;
//	v[5611] = 0e0;
//	v[5612] = 0e0;
//	v[5613] = 0e0;
//	v[5518] = 0e0;
//	v[5519] = 0e0;
//	v[5520] = 0e0;
//	v[5521] = 0e0;
//	v[5522] = v[366];
//	v[5523] = v[368];
//	v[5524] = 0e0;
//	v[5525] = 0e0;
//	v[5526] = 0e0;
//	v[5527] = 0e0;
//	v[5528] = 0e0;
//	v[5529] = 0e0;
//	v[5934] = 0e0;
//	v[5935] = 0e0;
//	v[5936] = 0e0;
//	v[5937] = 0e0;
//	v[5938] = -v[366];
//	v[5939] = -v[368];
//	v[5940] = 0e0;
//	v[5941] = 0e0;
//	v[5942] = 0e0;
//	v[5943] = 0e0;
//	v[5944] = 0e0;
//	v[5945] = 0e0;
//	v[5422] = 0e0;
//	v[5423] = 0e0;
//	v[5424] = 0e0;
//	v[5425] = v[363];
//	v[5426] = v[366];
//	v[5427] = v[368];
//	v[5428] = 0e0;
//	v[5429] = 0e0;
//	v[5430] = 0e0;
//	v[5431] = 0e0;
//	v[5432] = 0e0;
//	v[5433] = 0e0;
//	v[367] = alphaA[2] / 2e0;
//	v[5446] = 0e0;
//	v[5447] = 0e0;
//	v[5448] = 0e0;
//	v[5449] = v[367];
//	v[5450] = 0e0;
//	v[5451] = v[365];
//	v[5452] = 0e0;
//	v[5453] = 0e0;
//	v[5454] = 0e0;
//	v[5455] = 0e0;
//	v[5456] = 0e0;
//	v[5457] = 0e0;
//	v[5458] = 0e0;
//	v[5459] = 0e0;
//	v[5460] = 0e0;
//	v[5461] = 0e0;
//	v[5462] = v[367];
//	v[5463] = v[364];
//	v[5464] = 0e0;
//	v[5465] = 0e0;
//	v[5466] = 0e0;
//	v[5467] = 0e0;
//	v[5468] = 0e0;
//	v[5469] = 0e0;
//	v[161] = (alphaA[1] * alphaA[2]) / 2e0;
//	v[410] = alphaA[0] + v[161];
//	v[402] = -alphaA[0] + v[161];
//	v[159] = (alphaA[0] * alphaA[2]) / 2e0;
//	v[406] = -alphaA[1] + v[159];
//	v[389] = alphaA[1] + v[159];
//	v[150] = Power(alphaA[2], 2);
//	v[1087] = 4e0 + v[149] + v[150] + v[156];
//	v[3675] = 1e0 / Power(v[1087], 4);
//	v[1562] = 1e0 / Power(v[1087], 3);
//	v[1634] = -8e0*v[1562] * v[363];
//	v[1632] = 8e0*v[1562] * v[366];
//	v[1629] = -8e0*v[1562] * v[368];
//	v[962] = 1e0 / Power(v[1087], 2);
//	v[4135] = 2e0*v[962];
//	v[4132] = 4e0*v[962];
//	v[3992] = 2e0*v[962];
//	v[3985] = 4e0*v[962];
//	v[397] = -v[150] - v[156];
//	v[378] = -v[149] - v[150];
//	v[377] = 4e0*v[368] * v[962];
//	v[419] = -(v[377] * v[415]) / 2e0;
//	v[376] = -4e0*v[366] * v[962];
//	v[399] = (v[376] * v[397]) / 2e0;
//	v[375] = 4e0*v[363] * v[962];
//	v[379] = -(v[375] * v[378]) / 2e0;
//	v[256] = (*a4)*alphaA[0] + (*a5)*dalphaiA[0] + (*a6)*ddalphaiA[0];
//	v[262] = (*a4)*alphaA[1] + (*a5)*dalphaiA[1] + (*a6)*ddalphaiA[1];
//	v[264] = (*a4)*alphaA[2] + (*a5)*dalphaiA[2] + (*a6)*ddalphaiA[2];
//	v[133] = -cAp[0] / 2e0;
//	v[135] = -cAp[1] / 2e0;
//	v[138] = -cAi[0] / 2e0;
//	v[140] = -cAi[1] / 2e0;
//	v[353] = -x1B[0] / 2e0;
//	v[355] = -x1B[1] / 2e0;
//	v[357] = -x1B[2] / 2e0;
//	v[346] = v[353] + x2B[0] / 2e0;
//	v[347] = v[355] + x2B[1] / 2e0;
//	v[348] = v[357] + x2B[2] / 2e0;
//	v[354] = v[353] + x3B[0] / 2e0;
//	v[356] = v[355] + x3B[1] / 2e0;
//	v[358] = v[357] + x3B[2] / 2e0;
//	v[371] = alphaB[0] / 2e0;
//	v[369] = 2e0*alphaB[0];
//	v[212] = Power(alphaB[0], 2);
//	v[372] = 2e0*alphaB[1];
//	v[5722] = 0e0;
//	v[5723] = 0e0;
//	v[5724] = 0e0;
//	v[5725] = 0e0;
//	v[5726] = 0e0;
//	v[5727] = 0e0;
//	v[5728] = 0e0;
//	v[5729] = 0e0;
//	v[5730] = 0e0;
//	v[5731] = v[369];
//	v[5732] = v[372];
//	v[5733] = 0e0;
//	v[5958] = 0e0;
//	v[5959] = 0e0;
//	v[5960] = 0e0;
//	v[5961] = 0e0;
//	v[5962] = 0e0;
//	v[5963] = 0e0;
//	v[5964] = 0e0;
//	v[5965] = 0e0;
//	v[5966] = 0e0;
//	v[5967] = -v[369];
//	v[5968] = -v[372];
//	v[5969] = 0e0;
//	v[370] = alphaB[1] / 2e0;
//	v[5482] = 0e0;
//	v[5483] = 0e0;
//	v[5484] = 0e0;
//	v[5485] = 0e0;
//	v[5486] = 0e0;
//	v[5487] = 0e0;
//	v[5488] = 0e0;
//	v[5489] = 0e0;
//	v[5490] = 0e0;
//	v[5491] = v[370];
//	v[5492] = v[371];
//	v[5493] = 0e0;
//	v[210] = (alphaB[0] * alphaB[1]) / 2e0;
//	v[205] = Power(alphaB[1], 2);
//	v[514] = -v[205] - v[212];
//	v[492] = alphaB[2] + v[210];
//	v[482] = -alphaB[2] + v[210];
//	v[374] = 2e0*alphaB[2];
//	v[5638] = 0e0;
//	v[5639] = 0e0;
//	v[5640] = 0e0;
//	v[5641] = 0e0;
//	v[5642] = 0e0;
//	v[5643] = 0e0;
//	v[5644] = 0e0;
//	v[5645] = 0e0;
//	v[5646] = 0e0;
//	v[5647] = v[369];
//	v[5648] = 0e0;
//	v[5649] = v[374];
//	v[5554] = 0e0;
//	v[5555] = 0e0;
//	v[5556] = 0e0;
//	v[5557] = 0e0;
//	v[5558] = 0e0;
//	v[5559] = 0e0;
//	v[5560] = 0e0;
//	v[5561] = 0e0;
//	v[5562] = 0e0;
//	v[5563] = 0e0;
//	v[5564] = v[372];
//	v[5565] = v[374];
//	v[6006] = 0e0;
//	v[6007] = 0e0;
//	v[6008] = 0e0;
//	v[6009] = 0e0;
//	v[6010] = 0e0;
//	v[6011] = 0e0;
//	v[6012] = 0e0;
//	v[6013] = 0e0;
//	v[6014] = 0e0;
//	v[6015] = 0e0;
//	v[6016] = -v[372];
//	v[6017] = -v[374];
//	v[5470] = 0e0;
//	v[5471] = 0e0;
//	v[5472] = 0e0;
//	v[5473] = 0e0;
//	v[5474] = 0e0;
//	v[5475] = 0e0;
//	v[5476] = 0e0;
//	v[5477] = 0e0;
//	v[5478] = 0e0;
//	v[5479] = v[369];
//	v[5480] = v[372];
//	v[5481] = v[374];
//	v[373] = alphaB[2] / 2e0;
//	v[5494] = 0e0;
//	v[5495] = 0e0;
//	v[5496] = 0e0;
//	v[5497] = 0e0;
//	v[5498] = 0e0;
//	v[5499] = 0e0;
//	v[5500] = 0e0;
//	v[5501] = 0e0;
//	v[5502] = 0e0;
//	v[5503] = v[373];
//	v[5504] = 0e0;
//	v[5505] = v[371];
//	v[5506] = 0e0;
//	v[5507] = 0e0;
//	v[5508] = 0e0;
//	v[5509] = 0e0;
//	v[5510] = 0e0;
//	v[5511] = 0e0;
//	v[5512] = 0e0;
//	v[5513] = 0e0;
//	v[5514] = 0e0;
//	v[5515] = 0e0;
//	v[5516] = v[373];
//	v[5517] = v[370];
//	v[217] = (alphaB[1] * alphaB[2]) / 2e0;
//	v[509] = alphaB[0] + v[217];
//	v[501] = -alphaB[0] + v[217];
//	v[215] = (alphaB[0] * alphaB[2]) / 2e0;
//	v[505] = -alphaB[1] + v[215];
//	v[488] = alphaB[1] + v[215];
//	v[206] = Power(alphaB[2], 2);
//	v[1097] = 4e0 + v[205] + v[206] + v[212];
//	v[3650] = 1e0 / Power(v[1097], 4);
//	v[1564] = 1e0 / Power(v[1097], 3);
//	v[1660] = -8e0*v[1564] * v[369];
//	v[1658] = 8e0*v[1564] * v[372];
//	v[1655] = -8e0*v[1564] * v[374];
//	v[968] = 1e0 / Power(v[1097], 2);
//	v[4141] = 2e0*v[968];
//	v[4134] = 4e0*v[968];
//	v[3998] = 2e0*v[968];
//	v[3988] = 4e0*v[968];
//	v[496] = -v[206] - v[212];
//	v[477] = -v[205] - v[206];
//	v[476] = 4e0*v[374] * v[968];
//	v[518] = -(v[476] * v[514]) / 2e0;
//	v[475] = -4e0*v[372] * v[968];
//	v[498] = (v[475] * v[496]) / 2e0;
//	v[474] = 4e0*v[369] * v[968];
//	v[478] = -(v[474] * v[477]) / 2e0;
//	v[282] = (*a4)*alphaB[0] + (*a5)*dalphaiB[0] + (*a6)*ddalphaiB[0];
//	v[288] = (*a4)*alphaB[1] + (*a5)*dalphaiB[1] + (*a6)*ddalphaiB[1];
//	v[290] = (*a4)*alphaB[2] + (*a5)*dalphaiB[2] + (*a6)*ddalphaiB[2];
//	v[189] = -cBp[0] / 2e0;
//	v[191] = -cBp[1] / 2e0;
//	v[194] = -cBi[0] / 2e0;
//	v[196] = -cBi[1] / 2e0;
//	v[132] = v[133] + v[135];
//	v[134] = 0.5e0 - v[133];
//	v[136] = 0.5e0 - v[135];
//	v[137] = v[138] + v[140];
//	v[139] = 0.5e0 - v[138];
//	v[141] = 0.5e0 - v[140];
//	v[142] = v[132] * x1A[0] + v[134] * x2A[0] + v[136] * x3A[0];
//	v[143] = v[132] * x1A[1] + v[134] * x2A[1] + v[136] * x3A[1];
//	v[144] = v[132] * x1A[2] + v[134] * x2A[2] + v[136] * x3A[2];
//	v[861] = QAi[0][0] * v[142] + QAi[0][1] * v[143] + QAi[0][2] * v[144];
//	v[859] = QAi[1][0] * v[142] + QAi[1][1] * v[143] + QAi[1][2] * v[144];
//	v[857] = QAi[2][0] * v[142] + QAi[2][1] * v[143] + QAi[2][2] * v[144];
//	v[145] = v[137] * x1A[0] + v[139] * x2A[0] + v[141] * x3A[0];
//	v[146] = v[137] * x1A[1] + v[139] * x2A[1] + v[141] * x3A[1];
//	v[147] = v[137] * x1A[2] + v[139] * x2A[2] + v[141] * x3A[2];
//	v[148] = 4e0 / v[1087];
//	v[3983] = -v[148] / 2e0;
//	v[417] = -(v[148] * v[366]) / 2e0;
//	v[418] = (v[376] * v[415]) / 2e0 + v[417];
//	v[414] = -(v[148] * v[363]) / 2e0;
//	v[416] = v[414] - (v[375] * v[415]) / 2e0;
//	v[411] = v[148] - v[375] * v[410];
//	v[408] = -v[148] + v[376] * v[406];
//	v[403] = -v[148] - v[375] * v[402];
//	v[400] = -(v[148] * v[368]) / 2e0;
//	v[401] = -(v[377] * v[397]) / 2e0 + v[400];
//	v[398] = -(v[375] * v[397]) / 2e0 + v[414];
//	v[396] = v[148] - v[377] * v[393];
//	v[391] = v[148] + v[376] * v[389];
//	v[388] = -(v[148] * v[367]);
//	v[412] = -v[388] + v[376] * v[410];
//	v[457] = QAi[0][2] * v[408] + QAi[1][2] * v[412] + QAi[2][2] * v[418];
//	v[454] = QAi[0][1] * v[408] + QAi[1][1] * v[412] + QAi[2][1] * v[418];
//	v[451] = QAi[0][0] * v[408] + QAi[1][0] * v[412] + QAi[2][0] * v[418];
//	v[472] = v[142] * v[451] + v[143] * v[454] + v[144] * v[457];
//	v[407] = -v[388] - v[375] * v[406];
//	v[456] = QAi[0][2] * v[407] + QAi[1][2] * v[411] + QAi[2][2] * v[416];
//	v[453] = QAi[0][1] * v[407] + QAi[1][1] * v[411] + QAi[2][1] * v[416];
//	v[450] = QAi[0][0] * v[407] + QAi[1][0] * v[411] + QAi[2][0] * v[416];
//	v[471] = v[142] * v[450] + v[143] * v[453] + v[144] * v[456];
//	v[404] = -v[388] + v[376] * v[402];
//	v[390] = -v[388] - v[375] * v[389];
//	v[387] = -v[148] - v[377] * v[383];
//	v[385] = -(v[148] * v[365]);
//	v[409] = -v[385] - v[377] * v[406];
//	v[395] = -v[385] + v[376] * v[393];
//	v[442] = QAi[0][2] * v[395] + QAi[1][2] * v[399] + QAi[2][2] * v[404];
//	v[439] = QAi[0][1] * v[395] + QAi[1][1] * v[399] + QAi[2][1] * v[404];
//	v[436] = QAi[0][0] * v[395] + QAi[1][0] * v[399] + QAi[2][0] * v[404];
//	v[469] = v[142] * v[436] + v[143] * v[439] + v[144] * v[442];
//	v[392] = -v[385] - v[377] * v[389];
//	v[386] = v[376] * v[383] - v[385];
//	v[382] = v[148] * v[364];
//	v[413] = v[382] - v[377] * v[410];
//	v[458] = QAi[0][2] * v[409] + QAi[1][2] * v[413] + QAi[2][2] * v[419];
//	v[455] = QAi[0][1] * v[409] + QAi[1][1] * v[413] + QAi[2][1] * v[419];
//	v[452] = QAi[0][0] * v[409] + QAi[1][0] * v[413] + QAi[2][0] * v[419];
//	v[473] = v[142] * v[452] + v[143] * v[455] + v[144] * v[458];
//	v[405] = v[382] - v[377] * v[402];
//	v[443] = QAi[0][2] * v[396] + QAi[1][2] * v[401] + QAi[2][2] * v[405];
//	v[440] = QAi[0][1] * v[396] + QAi[1][1] * v[401] + QAi[2][1] * v[405];
//	v[437] = QAi[0][0] * v[396] + QAi[1][0] * v[401] + QAi[2][0] * v[405];
//	v[470] = v[142] * v[437] + v[143] * v[440] + v[144] * v[443];
//	v[394] = v[382] - v[375] * v[393];
//	v[441] = QAi[0][2] * v[394] + QAi[1][2] * v[398] + QAi[2][2] * v[403];
//	v[438] = QAi[0][1] * v[394] + QAi[1][1] * v[398] + QAi[2][1] * v[403];
//	v[435] = QAi[0][0] * v[394] + QAi[1][0] * v[398] + QAi[2][0] * v[403];
//	v[468] = v[142] * v[435] + v[143] * v[438] + v[144] * v[441];
//	v[384] = v[382] - v[375] * v[383];
//	v[426] = QAi[0][2] * v[379] + QAi[1][2] * v[384] + QAi[2][2] * v[390];
//	v[423] = QAi[0][1] * v[379] + QAi[1][1] * v[384] + QAi[2][1] * v[390];
//	v[420] = QAi[0][0] * v[379] + QAi[1][0] * v[384] + QAi[2][0] * v[390];
//	v[465] = v[142] * v[420] + v[143] * v[423] + v[144] * v[426];
//	v[381] = -(v[377] * v[378]) / 2e0 + v[400];
//	v[428] = QAi[0][2] * v[381] + QAi[1][2] * v[387] + QAi[2][2] * v[392];
//	v[425] = QAi[0][1] * v[381] + QAi[1][1] * v[387] + QAi[2][1] * v[392];
//	v[422] = QAi[0][0] * v[381] + QAi[1][0] * v[387] + QAi[2][0] * v[392];
//	v[467] = v[142] * v[422] + v[143] * v[425] + v[144] * v[428];
//	v[380] = (v[376] * v[378]) / 2e0 + v[417];
//	v[427] = QAi[0][2] * v[380] + QAi[1][2] * v[386] + QAi[2][2] * v[391];
//	v[424] = QAi[0][1] * v[380] + QAi[1][1] * v[386] + QAi[2][1] * v[391];
//	v[421] = QAi[0][0] * v[380] + QAi[1][0] * v[386] + QAi[2][0] * v[391];
//	v[466] = v[142] * v[421] + v[143] * v[424] + v[144] * v[427];
//	v[253] = (v[148] * v[148]);
//	v[3717] = v[264] / v[253];
//	v[3716] = v[262] / v[253];
//	v[3715] = v[256] / v[253];
//	v[3677] = 1e0 / Power(v[253], 3);
//	v[151] = 1e0 + (v[148] * v[378]) / 2e0;
//	v[4119] = v[151] / v[253];
//	v[1405] = v[151] * v[3715];
//	v[152] = v[148] * v[383];
//	v[3826] = v[152] * v[262];
//	v[1407] = v[152] * v[3716];
//	v[2018] = v[1405] + v[1407];
//	v[153] = v[148] * v[389];
//	v[4157] = v[153] / v[253];
//	v[3825] = v[153] * v[264];
//	v[1408] = v[153] * v[3717];
//	v[2023] = v[1405] + v[1408];
//	v[2012] = v[1407] + v[2023];
//	v[155] = v[148] * v[393];
//	v[1412] = v[155] * v[3715];
//	v[157] = 1e0 + (v[148] * v[397]) / 2e0;
//	v[4118] = v[157] / v[253];
//	v[1401] = v[157] * v[3716];
//	v[2014] = v[1401] + v[1412];
//	v[158] = v[148] * v[402];
//	v[4158] = v[158] / v[253];
//	v[1402] = v[158] * v[3717];
//	v[2024] = v[1401] + v[1402];
//	v[2017] = v[1412] + v[2024];
//	v[160] = v[148] * v[406];
//	v[1415] = v[160] * v[3715];
//	v[162] = v[148] * v[410];
//	v[4153] = v[162] / v[253];
//	v[1397] = v[162] * v[3716];
//	v[163] = 1e0 + (v[148] * v[415]) / 2e0;
//	v[4122] = v[163] / v[253];
//	v[1398] = v[163] * v[3717];
//	v[4120] = v[1398] * v[253];
//	v[2022] = v[1397] + v[1398] + v[1415];
//	v[2019] = -v[1415] + v[2022];
//	v[2013] = -v[1397] + v[2022];
//	v[271] = -(v[382] * v[388]);
//	v[3728] = v[271] - v[375];
//	v[269] = v[385] * v[388];
//	v[3726] = v[269] - v[376];
//	v[260] = -(v[382] * v[385]);
//	v[3723] = v[260] - v[377];
//	v[167] = QAi[0][0] * v[151] + QAi[1][0] * v[152] + QAi[2][0] * v[153];
//	v[168] = QAi[0][1] * v[151] + QAi[1][1] * v[152] + QAi[2][1] * v[153];
//	v[169] = QAi[0][2] * v[151] + QAi[1][2] * v[152] + QAi[2][2] * v[153];
//	v[342] = v[167] * v[337] + v[168] * v[339] + v[169] * v[341];
//	v[332] = v[167] * v[329] + v[168] * v[330] + v[169] * v[331];
//	v[170] = QAi[0][0] * v[155] + QAi[1][0] * v[157] + QAi[2][0] * v[158];
//	v[171] = QAi[0][1] * v[155] + QAi[1][1] * v[157] + QAi[2][1] * v[158];
//	v[172] = QAi[0][2] * v[155] + QAi[1][2] * v[157] + QAi[2][2] * v[158];
//	v[343] = v[170] * v[337] + v[171] * v[339] + v[172] * v[341];
//	v[333] = v[170] * v[329] + v[171] * v[330] + v[172] * v[331];
//	v[173] = QAi[0][0] * v[160] + QAi[1][0] * v[162] + QAi[2][0] * v[163];
//	v[174] = QAi[0][1] * v[160] + QAi[1][1] * v[162] + QAi[2][1] * v[163];
//	v[175] = QAi[0][2] * v[160] + QAi[1][2] * v[162] + QAi[2][2] * v[163];
//	v[344] = v[173] * v[337] + v[174] * v[339] + v[175] * v[341];
//	v[334] = v[173] * v[329] + v[174] * v[330] + v[175] * v[331];
//	v[188] = v[189] + v[191];
//	v[190] = 0.5e0 - v[189];
//	v[192] = 0.5e0 - v[191];
//	v[193] = v[194] + v[196];
//	v[195] = 0.5e0 - v[194];
//	v[197] = 0.5e0 - v[196];
//	v[198] = v[188] * x1B[0] + v[190] * x2B[0] + v[192] * x3B[0];
//	v[199] = v[188] * x1B[1] + v[190] * x2B[1] + v[192] * x3B[1];
//	v[200] = v[188] * x1B[2] + v[190] * x2B[2] + v[192] * x3B[2];
//	v[855] = -(QBi[0][0] * v[198]) - QBi[0][1] * v[199] - QBi[0][2] * v[200];
//	v[853] = -(QBi[1][0] * v[198]) - QBi[1][1] * v[199] - QBi[1][2] * v[200];
//	v[851] = -(QBi[2][0] * v[198]) - QBi[2][1] * v[199] - QBi[2][2] * v[200];
//	v[201] = v[193] * x1B[0] + v[195] * x2B[0] + v[197] * x3B[0];
//	v[202] = v[193] * x1B[1] + v[195] * x2B[1] + v[197] * x3B[1];
//	v[203] = v[193] * x1B[2] + v[195] * x2B[2] + v[197] * x3B[2];
//	v[204] = 4e0 / v[1097];
//	v[3986] = -v[204] / 2e0;
//	v[516] = -(v[204] * v[372]) / 2e0;
//	v[517] = (v[475] * v[514]) / 2e0 + v[516];
//	v[513] = -(v[204] * v[369]) / 2e0;
//	v[515] = v[513] - (v[474] * v[514]) / 2e0;
//	v[510] = v[204] - v[474] * v[509];
//	v[507] = -v[204] + v[475] * v[505];
//	v[502] = -v[204] - v[474] * v[501];
//	v[499] = -(v[204] * v[374]) / 2e0;
//	v[500] = -(v[476] * v[496]) / 2e0 + v[499];
//	v[497] = -(v[474] * v[496]) / 2e0 + v[513];
//	v[495] = v[204] - v[476] * v[492];
//	v[490] = v[204] + v[475] * v[488];
//	v[487] = -(v[204] * v[373]);
//	v[511] = -v[487] + v[475] * v[509];
//	v[556] = QBi[0][2] * v[507] + QBi[1][2] * v[511] + QBi[2][2] * v[517];
//	v[553] = QBi[0][1] * v[507] + QBi[1][1] * v[511] + QBi[2][1] * v[517];
//	v[550] = QBi[0][0] * v[507] + QBi[1][0] * v[511] + QBi[2][0] * v[517];
//	v[571] = v[198] * v[550] + v[199] * v[553] + v[200] * v[556];
//	v[506] = -v[487] - v[474] * v[505];
//	v[555] = QBi[0][2] * v[506] + QBi[1][2] * v[510] + QBi[2][2] * v[515];
//	v[552] = QBi[0][1] * v[506] + QBi[1][1] * v[510] + QBi[2][1] * v[515];
//	v[549] = QBi[0][0] * v[506] + QBi[1][0] * v[510] + QBi[2][0] * v[515];
//	v[570] = v[198] * v[549] + v[199] * v[552] + v[200] * v[555];
//	v[503] = -v[487] + v[475] * v[501];
//	v[489] = -v[487] - v[474] * v[488];
//	v[486] = -v[204] - v[476] * v[482];
//	v[484] = -(v[204] * v[371]);
//	v[508] = -v[484] - v[476] * v[505];
//	v[494] = -v[484] + v[475] * v[492];
//	v[541] = QBi[0][2] * v[494] + QBi[1][2] * v[498] + QBi[2][2] * v[503];
//	v[538] = QBi[0][1] * v[494] + QBi[1][1] * v[498] + QBi[2][1] * v[503];
//	v[535] = QBi[0][0] * v[494] + QBi[1][0] * v[498] + QBi[2][0] * v[503];
//	v[568] = v[198] * v[535] + v[199] * v[538] + v[200] * v[541];
//	v[491] = -v[484] - v[476] * v[488];
//	v[485] = v[475] * v[482] - v[484];
//	v[481] = v[204] * v[370];
//	v[512] = v[481] - v[476] * v[509];
//	v[557] = QBi[0][2] * v[508] + QBi[1][2] * v[512] + QBi[2][2] * v[518];
//	v[554] = QBi[0][1] * v[508] + QBi[1][1] * v[512] + QBi[2][1] * v[518];
//	v[551] = QBi[0][0] * v[508] + QBi[1][0] * v[512] + QBi[2][0] * v[518];
//	v[572] = v[198] * v[551] + v[199] * v[554] + v[200] * v[557];
//	v[504] = v[481] - v[476] * v[501];
//	v[542] = QBi[0][2] * v[495] + QBi[1][2] * v[500] + QBi[2][2] * v[504];
//	v[539] = QBi[0][1] * v[495] + QBi[1][1] * v[500] + QBi[2][1] * v[504];
//	v[536] = QBi[0][0] * v[495] + QBi[1][0] * v[500] + QBi[2][0] * v[504];
//	v[569] = v[198] * v[536] + v[199] * v[539] + v[200] * v[542];
//	v[493] = v[481] - v[474] * v[492];
//	v[540] = QBi[0][2] * v[493] + QBi[1][2] * v[497] + QBi[2][2] * v[502];
//	v[537] = QBi[0][1] * v[493] + QBi[1][1] * v[497] + QBi[2][1] * v[502];
//	v[534] = QBi[0][0] * v[493] + QBi[1][0] * v[497] + QBi[2][0] * v[502];
//	v[567] = v[198] * v[534] + v[199] * v[537] + v[200] * v[540];
//	v[483] = v[481] - v[474] * v[482];
//	v[525] = QBi[0][2] * v[478] + QBi[1][2] * v[483] + QBi[2][2] * v[489];
//	v[522] = QBi[0][1] * v[478] + QBi[1][1] * v[483] + QBi[2][1] * v[489];
//	v[519] = QBi[0][0] * v[478] + QBi[1][0] * v[483] + QBi[2][0] * v[489];
//	v[564] = v[198] * v[519] + v[199] * v[522] + v[200] * v[525];
//	v[582] = -(v[342] * v[564]) - v[343] * v[567] - v[344] * v[570];
//	v[576] = -(v[332] * v[564]) - v[333] * v[567] - v[334] * v[570];
//	v[480] = -(v[476] * v[477]) / 2e0 + v[499];
//	v[527] = QBi[0][2] * v[480] + QBi[1][2] * v[486] + QBi[2][2] * v[491];
//	v[524] = QBi[0][1] * v[480] + QBi[1][1] * v[486] + QBi[2][1] * v[491];
//	v[521] = QBi[0][0] * v[480] + QBi[1][0] * v[486] + QBi[2][0] * v[491];
//	v[566] = v[198] * v[521] + v[199] * v[524] + v[200] * v[527];
//	v[584] = -(v[342] * v[566]) - v[343] * v[569] - v[344] * v[572];
//	v[578] = -(v[332] * v[566]) - v[333] * v[569] - v[334] * v[572];
//	v[479] = (v[475] * v[477]) / 2e0 + v[516];
//	v[526] = QBi[0][2] * v[479] + QBi[1][2] * v[485] + QBi[2][2] * v[490];
//	v[523] = QBi[0][1] * v[479] + QBi[1][1] * v[485] + QBi[2][1] * v[490];
//	v[520] = QBi[0][0] * v[479] + QBi[1][0] * v[485] + QBi[2][0] * v[490];
//	v[565] = v[198] * v[520] + v[199] * v[523] + v[200] * v[526];
//	v[583] = -(v[342] * v[565]) - v[343] * v[568] - v[344] * v[571];
//	v[577] = -(v[332] * v[565]) - v[333] * v[568] - v[334] * v[571];
//	v[279] = (v[204] * v[204]);
//	v[3720] = v[290] / v[279];
//	v[3719] = v[288] / v[279];
//	v[3718] = v[282] / v[279];
//	v[3652] = 1e0 / Power(v[279], 3);
//	v[207] = 1e0 + (v[204] * v[477]) / 2e0;
//	v[4110] = v[207] / v[279];
//	v[1381] = v[207] * v[3718];
//	v[208] = v[204] * v[482];
//	v[3808] = v[208] * v[288];
//	v[1383] = v[208] * v[3719];
//	v[2033] = v[1381] + v[1383];
//	v[209] = v[204] * v[488];
//	v[4169] = v[209] / v[279];
//	v[3807] = v[209] * v[290];
//	v[1384] = v[209] * v[3720];
//	v[2038] = v[1381] + v[1384];
//	v[2027] = v[1383] + v[2038];
//	v[211] = v[204] * v[492];
//	v[1388] = v[211] * v[3718];
//	v[213] = 1e0 + (v[204] * v[496]) / 2e0;
//	v[4109] = v[213] / v[279];
//	v[1377] = v[213] * v[3719];
//	v[2029] = v[1377] + v[1388];
//	v[214] = v[204] * v[501];
//	v[4170] = v[214] / v[279];
//	v[1378] = v[214] * v[3720];
//	v[2039] = v[1377] + v[1378];
//	v[2032] = v[1388] + v[2039];
//	v[216] = v[204] * v[505];
//	v[1391] = v[216] * v[3718];
//	v[218] = v[204] * v[509];
//	v[4165] = v[218] / v[279];
//	v[1373] = v[218] * v[3719];
//	v[219] = 1e0 + (v[204] * v[514]) / 2e0;
//	v[4113] = v[219] / v[279];
//	v[1374] = v[219] * v[3720];
//	v[4111] = v[1374] * v[279];
//	v[2037] = v[1373] + v[1374] + v[1391];
//	v[2034] = -v[1391] + v[2037];
//	v[2028] = -v[1373] + v[2037];
//	v[297] = -(v[481] * v[487]);
//	v[3736] = v[297] - v[474];
//	v[295] = v[484] * v[487];
//	v[3734] = v[295] - v[475];
//	v[286] = -(v[481] * v[484]);
//	v[3731] = v[286] - v[476];
//	v[223] = QBi[0][0] * v[207] + QBi[1][0] * v[208] + QBi[2][0] * v[209];
//	v[224] = QBi[0][1] * v[207] + QBi[1][1] * v[208] + QBi[2][1] * v[209];
//	v[225] = QBi[0][2] * v[207] + QBi[1][2] * v[208] + QBi[2][2] * v[209];
//	v[359] = v[223] * v[354] + v[224] * v[356] + v[225] * v[358];
//	v[349] = v[223] * v[346] + v[224] * v[347] + v[225] * v[348];
//	v[664] = -(invH[3][0] * v[332]) - invH[3][1] * v[342] + invH[3][2] * v[349] + invH[3][3] * v[359];
//	v[649] = -(invH[2][0] * v[332]) - invH[2][1] * v[342] + invH[2][2] * v[349] + invH[2][3] * v[359];
//	v[634] = -(invH[1][0] * v[332]) - invH[1][1] * v[342] + invH[1][2] * v[349] + invH[1][3] * v[359];
//	v[619] = -(invH[0][0] * v[332]) - invH[0][1] * v[342] + invH[0][2] * v[349] + invH[0][3] * v[359];
//	v[226] = QBi[0][0] * v[211] + QBi[1][0] * v[213] + QBi[2][0] * v[214];
//	v[227] = QBi[0][1] * v[211] + QBi[1][1] * v[213] + QBi[2][1] * v[214];
//	v[228] = QBi[0][2] * v[211] + QBi[1][2] * v[213] + QBi[2][2] * v[214];
//	v[360] = v[226] * v[354] + v[227] * v[356] + v[228] * v[358];
//	v[350] = v[226] * v[346] + v[227] * v[347] + v[228] * v[348];
//	v[666] = -(invH[3][0] * v[333]) - invH[3][1] * v[343] + invH[3][2] * v[350] + invH[3][3] * v[360];
//	v[651] = -(invH[2][0] * v[333]) - invH[2][1] * v[343] + invH[2][2] * v[350] + invH[2][3] * v[360];
//	v[636] = -(invH[1][0] * v[333]) - invH[1][1] * v[343] + invH[1][2] * v[350] + invH[1][3] * v[360];
//	v[621] = -(invH[0][0] * v[333]) - invH[0][1] * v[343] + invH[0][2] * v[350] + invH[0][3] * v[360];
//	v[229] = QBi[0][0] * v[216] + QBi[1][0] * v[218] + QBi[2][0] * v[219];
//	v[230] = QBi[0][1] * v[216] + QBi[1][1] * v[218] + QBi[2][1] * v[219];
//	v[231] = QBi[0][2] * v[216] + QBi[1][2] * v[218] + QBi[2][2] * v[219];
//	v[361] = v[229] * v[354] + v[230] * v[356] + v[231] * v[358];
//	v[593] = -(v[359] * v[467]) - v[360] * v[470] - v[361] * v[473];
//	v[592] = -(v[359] * v[466]) - v[360] * v[469] - v[361] * v[472];
//	v[591] = -(v[359] * v[465]) - v[360] * v[468] - v[361] * v[471];
//	v[351] = v[229] * v[346] + v[230] * v[347] + v[231] * v[348];
//	v[668] = -(invH[3][0] * v[334]) - invH[3][1] * v[344] + invH[3][2] * v[351] + invH[3][3] * v[361];
//	v[653] = -(invH[2][0] * v[334]) - invH[2][1] * v[344] + invH[2][2] * v[351] + invH[2][3] * v[361];
//	v[638] = -(invH[1][0] * v[334]) - invH[1][1] * v[344] + invH[1][2] * v[351] + invH[1][3] * v[361];
//	v[623] = -(invH[0][0] * v[334]) - invH[0][1] * v[344] + invH[0][2] * v[351] + invH[0][3] * v[361];
//	v[587] = -(v[349] * v[467]) - v[350] * v[470] - v[351] * v[473];
//	v[586] = -(v[349] * v[466]) - v[350] * v[469] - v[351] * v[472];
//	v[585] = -(v[349] * v[465]) - v[350] * v[468] - v[351] * v[471];
//	v[3762] = uA[0] - uB[0] + xAi[0] - xBi[0];
//	v[3761] = uA[1] - uB[1] + xAi[1] - xBi[1];
//	v[3760] = uA[2] - uB[2] + xAi[2] - xBi[2];
//	v[244] = (*a6)*dduiA[0] + (*a5)*duiA[0] + (*a4)*uA[0];
//	v[245] = (*a6)*dduiA[1] + (*a5)*duiA[1] + (*a4)*uA[1];
//	v[246] = (*a6)*dduiA[2] + (*a5)*duiA[2] + (*a4)*uA[2];
//	v[247] = (*a6)*dduiB[0] + (*a5)*duiB[0] + (*a4)*uB[0];
//	v[3370] = v[244] - v[247];
//	v[248] = (*a6)*dduiB[1] + (*a5)*duiB[1] + (*a4)*uB[1];
//	v[3354] = v[245] - v[248];
//	v[249] = (*a6)*dduiB[2] + (*a5)*duiB[2] + (*a4)*uB[2];
//	v[3335] = v[246] - v[249];
//	v[251] = v[269] + v[376];
//	v[3727] = v[251] / v[253];
//	v[3721] = v[163] * v[251];
//	v[252] = v[260] + v[377];
//	v[3724] = v[252] / v[253];
//	v[3722] = v[157] * v[252];
//	v[254] = v[253] + (v[385] * v[385]);
//	v[3680] = -(v[264] * v[3721]) - v[262] * v[3722] - v[254] * (v[256] + v[3825] + v[3826]);
//	v[2682] = v[254] / v[253];
//	v[1574] = v[3726] / v[253];
//	v[1572] = v[3723] / v[253];
//	v[6246] = 0e0;
//	v[6247] = 0e0;
//	v[6248] = 0e0;
//	v[6249] = 0e0;
//	v[6250] = v[1572];
//	v[6251] = v[1574];
//	v[6252] = 0e0;
//	v[6253] = 0e0;
//	v[6254] = 0e0;
//	v[6255] = 0e0;
//	v[6256] = 0e0;
//	v[6257] = 0e0;
//	v[255] = v[1572] * v[262] + v[1574] * v[264] + v[256] * v[2682];
//	v[257] = v[253] + (v[382] * v[382]);
//	v[3814] = v[155] * v[257] + v[151] * v[3723];
//	v[263] = v[271] + v[375];
//	v[3725] = v[158] * v[257] + v[163] * v[263];
//	v[3679] = -(v[257] * v[262]) - v[264] * v[3725] - v[256] * v[3814];
//	v[2683] = v[257] / v[253];
//	v[1573] = v[3728] / v[253];
//	v[6258] = 0e0;
//	v[6259] = 0e0;
//	v[6260] = 0e0;
//	v[6261] = v[3724];
//	v[6262] = 0e0;
//	v[6263] = v[1573];
//	v[6264] = 0e0;
//	v[6265] = 0e0;
//	v[6266] = 0e0;
//	v[6267] = 0e0;
//	v[6268] = 0e0;
//	v[6269] = 0e0;
//	v[265] = v[1573] * v[264] + v[262] * v[2683] + v[256] * v[3724];
//	v[266] = v[253] + (v[388] * v[388]);
//	v[3815] = v[160] * v[266] + v[151] * v[3726];
//	v[3816] = v[162] * v[266] + v[157] * v[3728];
//	v[3678] = -(v[264] * v[266]) - v[256] * v[3815] - v[262] * v[3816];
//	v[2684] = v[266] / v[253];
//	v[1571] = v[263] / v[253];
//	v[6270] = 0e0;
//	v[6271] = 0e0;
//	v[6272] = 0e0;
//	v[6273] = v[3727];
//	v[6274] = v[1571];
//	v[6275] = 0e0;
//	v[6276] = 0e0;
//	v[6277] = 0e0;
//	v[6278] = 0e0;
//	v[6279] = 0e0;
//	v[6280] = 0e0;
//	v[6281] = 0e0;
//	v[275] = v[1571] * v[262] + v[264] * v[2684] + v[256] * v[3727];
//	v[277] = v[295] + v[475];
//	v[3735] = v[277] / v[279];
//	v[3729] = v[219] * v[277];
//	v[278] = v[286] + v[476];
//	v[3732] = v[278] / v[279];
//	v[3730] = v[213] * v[278];
//	v[280] = v[279] + (v[484] * v[484]);
//	v[3655] = -(v[290] * v[3729]) - v[288] * v[3730] - v[280] * (v[282] + v[3807] + v[3808]);
//	v[2685] = v[280] / v[279];
//	v[1580] = v[3734] / v[279];
//	v[1578] = v[3731] / v[279];
//	v[6282] = 0e0;
//	v[6283] = 0e0;
//	v[6284] = 0e0;
//	v[6285] = 0e0;
//	v[6286] = 0e0;
//	v[6287] = 0e0;
//	v[6288] = 0e0;
//	v[6289] = 0e0;
//	v[6290] = 0e0;
//	v[6291] = 0e0;
//	v[6292] = v[1578];
//	v[6293] = v[1580];
//	v[281] = v[2685] * v[282] + v[1578] * v[288] + v[1580] * v[290];
//	v[283] = v[279] + (v[481] * v[481]);
//	v[3796] = v[211] * v[283] + v[207] * v[3731];
//	v[289] = v[297] + v[474];
//	v[3733] = v[214] * v[283] + v[219] * v[289];
//	v[3654] = -(v[283] * v[288]) - v[290] * v[3733] - v[282] * v[3796];
//	v[2686] = v[283] / v[279];
//	v[1579] = v[3736] / v[279];
//	v[6294] = 0e0;
//	v[6295] = 0e0;
//	v[6296] = 0e0;
//	v[6297] = 0e0;
//	v[6298] = 0e0;
//	v[6299] = 0e0;
//	v[6300] = 0e0;
//	v[6301] = 0e0;
//	v[6302] = 0e0;
//	v[6303] = v[3732];
//	v[6304] = 0e0;
//	v[6305] = v[1579];
//	v[291] = v[2686] * v[288] + v[1579] * v[290] + v[282] * v[3732];
//	v[292] = v[279] + (v[487] * v[487]);
//	v[3797] = v[216] * v[292] + v[207] * v[3734];
//	v[3798] = v[218] * v[292] + v[213] * v[3736];
//	v[3653] = -(v[290] * v[292]) - v[282] * v[3797] - v[288] * v[3798];
//	v[2687] = v[292] / v[279];
//	v[1577] = v[289] / v[279];
//	v[6306] = 0e0;
//	v[6307] = 0e0;
//	v[6308] = 0e0;
//	v[6309] = 0e0;
//	v[6310] = 0e0;
//	v[6311] = 0e0;
//	v[6312] = 0e0;
//	v[6313] = 0e0;
//	v[6314] = 0e0;
//	v[6315] = v[3735];
//	v[6316] = v[1577];
//	v[6317] = 0e0;
//	v[301] = v[1577] * v[288] + v[2687] * v[290] + v[282] * v[3735];
//	v[2100] = v[3370] + v[255] * v[465] + v[265] * v[466] + v[275] * v[467] - v[281] * v[564] - v[291] * v[565]
//		- v[301] * v[566];
//	v[3739] = v[2100] - v[3370];
//	v[2095] = v[3354] + v[255] * v[468] + v[265] * v[469] + v[275] * v[470] - v[281] * v[567] - v[291] * v[568]
//		- v[301] * v[569];
//	v[3738] = v[2095] - v[3354];
//	v[2090] = v[3335] + v[255] * v[471] + v[265] * v[472] + v[275] * v[473] - v[281] * v[570] - v[291] * v[571]
//		- v[301] * v[572];
//	v[3737] = v[2090] - v[3335];
//	v[302] = v[142] * v[167] + v[143] * v[168] + v[144] * v[169] - v[198] * v[223] - v[199] * v[224] - v[200] * v[225] + v[3762];
//	v[303] = v[142] * v[170] + v[143] * v[171] + v[144] * v[172] - v[198] * v[226] - v[199] * v[227] - v[200] * v[228] + v[3761];
//	v[304] = v[142] * v[173] + v[143] * v[174] + v[144] * v[175] - v[198] * v[229] - v[199] * v[230] - v[200] * v[231] + v[3760];
//	v[727] = sqrt((v[302] * v[302]) + (v[303] * v[303]) + (v[304] * v[304]));
//	v[1888] = 1e0 / Power(v[727], 2);
//	v[596] = -(v[302] * (v[354] * v[521] + v[356] * v[524] + v[358] * v[527])) - v[303] * (v[354] * v[536] + v[356] * v[539]
//		+ v[358] * v[542]) - v[304] * (v[354] * v[551] + v[356] * v[554] + v[358] * v[557]) + v[359] * v[566] + v[360] * v[569]
//		+ v[361] * v[572];
//	v[595] = -(v[302] * (v[354] * v[520] + v[356] * v[523] + v[358] * v[526])) - v[303] * (v[354] * v[535] + v[356] * v[538]
//		+ v[358] * v[541]) - v[304] * (v[354] * v[550] + v[356] * v[553] + v[358] * v[556]) + v[359] * v[565] + v[360] * v[568]
//		+ v[361] * v[571];
//	v[594] = -(v[302] * (v[354] * v[519] + v[356] * v[522] + v[358] * v[525])) - v[303] * (v[354] * v[534] + v[356] * v[537]
//		+ v[358] * v[540]) - v[304] * (v[354] * v[549] + v[356] * v[552] + v[358] * v[555]) + v[359] * v[564] + v[360] * v[567]
//		+ v[361] * v[570];
//	v[590] = -(v[302] * (v[346] * v[521] + v[347] * v[524] + v[348] * v[527])) - v[303] * (v[346] * v[536] + v[347] * v[539]
//		+ v[348] * v[542]) - v[304] * (v[346] * v[551] + v[347] * v[554] + v[348] * v[557]) + v[349] * v[566] + v[350] * v[569]
//		+ v[351] * v[572];
//	v[589] = -(v[302] * (v[346] * v[520] + v[347] * v[523] + v[348] * v[526])) - v[303] * (v[346] * v[535] + v[347] * v[538]
//		+ v[348] * v[541]) - v[304] * (v[346] * v[550] + v[347] * v[553] + v[348] * v[556]) + v[349] * v[565] + v[350] * v[568]
//		+ v[351] * v[571];
//	v[588] = -(v[302] * (v[346] * v[519] + v[347] * v[522] + v[348] * v[525])) - v[303] * (v[346] * v[534] + v[347] * v[537]
//		+ v[348] * v[540]) - v[304] * (v[346] * v[549] + v[347] * v[552] + v[348] * v[555]) + v[349] * v[564] + v[350] * v[567]
//		+ v[351] * v[570];
//	v[581] = v[302] * (v[337] * v[422] + v[339] * v[425] + v[341] * v[428]) + v[303] * (v[337] * v[437] + v[339] * v[440]
//		+ v[341] * v[443]) + v[304] * (v[337] * v[452] + v[339] * v[455] + v[341] * v[458]) + v[342] * v[467] + v[343] * v[470]
//		+ v[344] * v[473];
//	v[580] = v[302] * (v[337] * v[421] + v[339] * v[424] + v[341] * v[427]) + v[303] * (v[337] * v[436] + v[339] * v[439]
//		+ v[341] * v[442]) + v[304] * (v[337] * v[451] + v[339] * v[454] + v[341] * v[457]) + v[342] * v[466] + v[343] * v[469]
//		+ v[344] * v[472];
//	v[579] = v[302] * (v[337] * v[420] + v[339] * v[423] + v[341] * v[426]) + v[303] * (v[337] * v[435] + v[339] * v[438]
//		+ v[341] * v[441]) + v[304] * (v[337] * v[450] + v[339] * v[453] + v[341] * v[456]) + v[342] * v[465] + v[343] * v[468]
//		+ v[344] * v[471];
//	v[575] = v[302] * (v[329] * v[422] + v[330] * v[425] + v[331] * v[428]) + v[303] * (v[329] * v[437] + v[330] * v[440]
//		+ v[331] * v[443]) + v[304] * (v[329] * v[452] + v[330] * v[455] + v[331] * v[458]) + v[332] * v[467] + v[333] * v[470]
//		+ v[334] * v[473];
//	v[574] = v[302] * (v[329] * v[421] + v[330] * v[424] + v[331] * v[427]) + v[303] * (v[329] * v[436] + v[330] * v[439]
//		+ v[331] * v[442]) + v[304] * (v[329] * v[451] + v[330] * v[454] + v[331] * v[457]) + v[332] * v[466] + v[333] * v[469]
//		+ v[334] * v[472];
//	v[573] = v[302] * (v[329] * v[420] + v[330] * v[423] + v[331] * v[426]) + v[303] * (v[329] * v[435] + v[330] * v[438]
//		+ v[331] * v[441]) + v[304] * (v[329] * v[450] + v[330] * v[453] + v[331] * v[456]) + v[332] * v[465] + v[333] * v[468]
//		+ v[334] * v[471];
//	v[305] = QAi[0][0] * v[145] + QAi[0][1] * v[146] + QAi[0][2] * v[147] - QBi[0][0] * v[201] - QBi[0][1] * v[202]
//		- QBi[0][2] * v[203] + xAi[0] - xBi[0];
//	v[306] = QAi[1][0] * v[145] + QAi[1][1] * v[146] + QAi[1][2] * v[147] - QBi[1][0] * v[201] - QBi[1][1] * v[202]
//		- QBi[1][2] * v[203] + xAi[1] - xBi[1];
//	v[307] = QAi[2][0] * v[145] + QAi[2][1] * v[146] + QAi[2][2] * v[147] - QBi[2][0] * v[201] - QBi[2][1] * v[202]
//		- QBi[2][2] * v[203] + xAi[2] - xBi[2];
//	if (v[727] > 0.1e-7) { v01 = 1e0 / v[727]; v02 = (-(v01 / v[727])); v03 = (2e0*v01) / Power(v[727], 2); }
//	else {
//		v01 = (12500000e0 / 3e0)*(24e0 - (-0.1e-7 + v[727])*(0.24e10 - 2e0*(-1e0 + 100000000e0*v[727])*(0.2399999997e10
//			- 0.1199999994e18*v[727] - 0.3e17*(v[727] * v[727]))));
//		v02 = (-50000000e0 / 3e0)*(0.3599999994e10 - 0.4799999982e18*v[727] + 0.6e25*Power(v[727], 3)
//			+ 0.1799999982e26*(v[727] * v[727]));
//		v03 = 0.1e17*(799999997e0 - 0.599999994e17*v[727] - 0.3e17*(v[727] * v[727]));
//	};
//	v[314] = v03;
//	v[315] = v02;
//	v[316] = v01;
//	v[317] = v[302] * v[316];
//	v[3352] = (v[245] - v[248])*v[317];
//	v[3333] = (v[246] - v[249])*v[317];
//	v[2103] = v[3333] + v[317] * v[3737];
//	v[2102] = v[3352] + v[317] * v[3738];
//	v[1223] = -(v[317] * v[566]);
//	v[1222] = -(v[317] * v[565]);
//	v[1221] = -(v[317] * v[564]);
//	v[1220] = v[317] * v[467];
//	v[1219] = v[317] * v[466];
//	v[1218] = v[317] * v[465];
//	v[713] = (v[317] * v[317]);
//	v[318] = v[303] * v[316];
//	v[3368] = (v[244] - v[247])*v[318];
//	v[3334] = (v[246] - v[249])*v[318];
//	v[2098] = v[3334] + v[318] * v[3737];
//	v[2096] = v[3368] + v[318] * v[3739];
//	v[1239] = -(v[318] * v[569]);
//	v[2336] = v[1223] + v[1239];
//	v[1237] = -(v[318] * v[568]);
//	v[2340] = v[1222] + v[1237];
//	v[1235] = -(v[318] * v[567]);
//	v[2344] = v[1221] + v[1235];
//	v[1233] = v[318] * v[470];
//	v[2348] = v[1220] + v[1233];
//	v[1231] = v[318] * v[469];
//	v[2352] = v[1219] + v[1231];
//	v[1229] = v[318] * v[468];
//	v[2356] = v[1218] + v[1229];
//	v[3758] = (v[245] - v[248])*v[318];
//	v[2265] = v[244] * v[317] - v[247] * v[317] + v[3758];
//	v[2093] = v[2265] + v[2356] * v[255] + v[2352] * v[265] + v[2348] * v[275] + v[2344] * v[281] + v[2340] * v[291]
//		+ v[2336] * v[301];
//	v[715] = (v[318] * v[318]);
//	v[319] = v[304] * v[316];
//	v[3745] = -(v[319] * v[572]);
//	v[3747] = -v[1239] - v[3745];
//	v[3746] = v[1223] + v[3745];
//	v[3744] = -(v[319] * v[571]);
//	v[3749] = -v[1237] - v[3744];
//	v[3748] = v[1222] + v[3744];
//	v[3743] = -(v[319] * v[570]);
//	v[3751] = -v[1235] - v[3743];
//	v[3750] = v[1221] + v[3743];
//	v[3742] = v[319] * v[473];
//	v[3753] = v[1233] + v[3742];
//	v[3752] = v[1220] + v[3742];
//	v[3741] = v[319] * v[472];
//	v[3755] = v[1231] + v[3741];
//	v[3754] = v[1219] + v[3741];
//	v[3740] = v[319] * v[471];
//	v[3757] = v[1229] + v[3740];
//	v[3756] = v[1218] + v[3740];
//	v[3369] = (v[244] - v[247])*v[319];
//	v[3353] = (v[245] - v[248])*v[319];
//	v[2092] = v[3353] + v[319] * v[3738];
//	v[2091] = v[3369] + v[319] * v[3739];
//	v[1262] = v[318] * v[3756] + v[468] * v[715];
//	v[1261] = v[317] * v[3757] + v[465] * v[713];
//	v[1258] = v[318] * v[3754] + v[469] * v[715];
//	v[1257] = v[317] * v[3755] + v[466] * v[713];
//	v[1254] = v[318] * v[3752] + v[470] * v[715];
//	v[1253] = v[317] * v[3753] + v[467] * v[713];
//	v[1250] = v[318] * v[3750] - v[567] * v[715];
//	v[1249] = -(v[317] * v[3751]) - v[564] * v[713];
//	v[1246] = v[318] * v[3748] - v[568] * v[715];
//	v[1245] = -(v[317] * v[3749]) - v[565] * v[713];
//	v[1242] = v[318] * v[3746] - v[569] * v[715];
//	v[1241] = -(v[317] * v[3747]) - v[566] * v[713];
//	v[2272] = v[246] * v[319] - v[249] * v[319] + v[3758];
//	v[2268] = v[2265] + v[2272] - 2e0*v[3758];
//	v[2101] = v[2272] - v[301] * v[3747] - v[291] * v[3749] - v[281] * v[3751] + v[275] * v[3753] + v[265] * v[3755]
//		+ v[255] * v[3757];
//	v[2097] = v[2268] + v[301] * v[3746] + v[291] * v[3748] + v[281] * v[3750] + v[275] * v[3752] + v[265] * v[3754]
//		+ v[255] * v[3756];
//	v[717] = (v[319] * v[319]);
//	v[1263] = v[2356] * v[319] + v[471] * v[717];
//	v[1259] = v[2352] * v[319] + v[472] * v[717];
//	v[1255] = v[2348] * v[319] + v[473] * v[717];
//	v[1251] = v[2344] * v[319] - v[570] * v[717];
//	v[1247] = v[2340] * v[319] - v[571] * v[717];
//	v[1243] = v[2336] * v[319] - v[572] * v[717];
//	v[320] = sqrt((v[305] * v[305]) + (v[306] * v[306]) + (v[307] * v[307]));
//	if (v[320] > 0.1e-7) { v04 = 1e0 / v[320]; v05 = (-(v04 / v[320])); v06 = (2e0*v04) / Power(v[320], 2); }
//	else {
//		v04 = (12500000e0 / 3e0)*(24e0 - (-0.1e-7 + v[320])*(0.24e10 - 2e0*(-1e0 + 100000000e0*v[320])*(0.2399999997e10
//			- 0.1199999994e18*v[320] - 0.3e17*(v[320] * v[320]))));
//		v05 = (-50000000e0 / 3e0)*(0.3599999994e10 - 0.4799999982e18*v[320] + 0.6e25*Power(v[320], 3)
//			+ 0.1799999982e26*(v[320] * v[320]));
//		v06 = 0.1e17*(799999997e0 - 0.599999994e17*v[320] - 0.3e17*(v[320] * v[320]));
//	};
//	v[325] = v04;
//	v[326] = v[305] * v[325];
//	v[327] = v[306] * v[325];
//	v[328] = v[307] * v[325];
//	v[616] = -(invH[0][0] * v[573]) - invH[0][1] * v[579] - invH[0][2] * v[585] - invH[0][3] * v[591];
//	v[617] = -(invH[0][0] * v[574]) - invH[0][1] * v[580] - invH[0][2] * v[586] - invH[0][3] * v[592];
//	v[618] = -(invH[0][0] * v[575]) - invH[0][1] * v[581] - invH[0][2] * v[587] - invH[0][3] * v[593];
//	v[625] = -(invH[0][0] * v[576]) - invH[0][1] * v[582] - invH[0][2] * v[588] - invH[0][3] * v[594];
//	v[626] = -(invH[0][0] * v[577]) - invH[0][1] * v[583] - invH[0][2] * v[589] - invH[0][3] * v[595];
//	v[627] = -(invH[0][0] * v[578]) - invH[0][1] * v[584] - invH[0][2] * v[590] - invH[0][3] * v[596];
//	v[1179] = -(v[255] * v[616]) - v[265] * v[617] - v[275] * v[618] - v[3370] * v[619] - v[3354] * v[621] - v[3335] * v[623]
//		- v[281] * v[625] - v[291] * v[626] - v[301] * v[627];
//	v[5370] = v[619];
//	v[5371] = v[621];
//	v[5372] = v[623];
//	v[5373] = v[616];
//	v[5374] = v[617];
//	v[5375] = v[618];
//	v[5376] = -v[619];
//	v[5377] = -v[621];
//	v[5378] = -v[623];
//	v[5379] = v[625];
//	v[5380] = v[626];
//	v[5381] = v[627];
//	v[631] = -(invH[1][0] * v[573]) - invH[1][1] * v[579] - invH[1][2] * v[585] - invH[1][3] * v[591];
//	v[632] = -(invH[1][0] * v[574]) - invH[1][1] * v[580] - invH[1][2] * v[586] - invH[1][3] * v[592];
//	v[633] = -(invH[1][0] * v[575]) - invH[1][1] * v[581] - invH[1][2] * v[587] - invH[1][3] * v[593];
//	v[640] = -(invH[1][0] * v[576]) - invH[1][1] * v[582] - invH[1][2] * v[588] - invH[1][3] * v[594];
//	v[641] = -(invH[1][0] * v[577]) - invH[1][1] * v[583] - invH[1][2] * v[589] - invH[1][3] * v[595];
//	v[642] = -(invH[1][0] * v[578]) - invH[1][1] * v[584] - invH[1][2] * v[590] - invH[1][3] * v[596];
//	v[1181] = -(v[255] * v[631]) - v[265] * v[632] - v[275] * v[633] - v[3370] * v[634] - v[3354] * v[636] - v[3335] * v[638]
//		- v[281] * v[640] - v[291] * v[641] - v[301] * v[642];
//	v[5358] = v[634];
//	v[5359] = v[636];
//	v[5360] = v[638];
//	v[5361] = v[631];
//	v[5362] = v[632];
//	v[5363] = v[633];
//	v[5364] = -v[634];
//	v[5365] = -v[636];
//	v[5366] = -v[638];
//	v[5367] = v[640];
//	v[5368] = v[641];
//	v[5369] = v[642];
//	v[646] = -(invH[2][0] * v[573]) - invH[2][1] * v[579] - invH[2][2] * v[585] - invH[2][3] * v[591];
//	v[647] = -(invH[2][0] * v[574]) - invH[2][1] * v[580] - invH[2][2] * v[586] - invH[2][3] * v[592];
//	v[648] = -(invH[2][0] * v[575]) - invH[2][1] * v[581] - invH[2][2] * v[587] - invH[2][3] * v[593];
//	v[655] = -(invH[2][0] * v[576]) - invH[2][1] * v[582] - invH[2][2] * v[588] - invH[2][3] * v[594];
//	v[656] = -(invH[2][0] * v[577]) - invH[2][1] * v[583] - invH[2][2] * v[589] - invH[2][3] * v[595];
//	v[657] = -(invH[2][0] * v[578]) - invH[2][1] * v[584] - invH[2][2] * v[590] - invH[2][3] * v[596];
//	v[1175] = v[255] * v[646] + v[265] * v[647] + v[275] * v[648] + v[3370] * v[649] + v[3354] * v[651] + v[3335] * v[653]
//		+ v[281] * v[655] + v[291] * v[656] + v[301] * v[657];
//	v[5382] = v[649];
//	v[5383] = v[651];
//	v[5384] = v[653];
//	v[5385] = v[646];
//	v[5386] = v[647];
//	v[5387] = v[648];
//	v[5388] = -v[649];
//	v[5389] = -v[651];
//	v[5390] = -v[653];
//	v[5391] = v[655];
//	v[5392] = v[656];
//	v[5393] = v[657];
//	v[1112] = v[334] * v[619] + v[344] * v[634] - v[351] * v[649] - v[361] * v[664];
//	v[1111] = v[333] * v[619] + v[343] * v[634] - v[350] * v[649] - v[360] * v[664];
//	v[1110] = v[332] * v[619] + v[342] * v[634] - v[349] * v[649] - v[359] * v[664];
//	v[1116] = v[334] * v[621] + v[344] * v[636] - v[351] * v[651] - v[361] * v[666];
//	v[1115] = v[333] * v[621] + v[343] * v[636] - v[350] * v[651] - v[360] * v[666];
//	v[1114] = v[332] * v[621] + v[342] * v[636] - v[349] * v[651] - v[359] * v[666];
//	v[1120] = v[334] * v[623] + v[344] * v[638] - v[351] * v[653] - v[361] * v[668];
//	v[1119] = v[333] * v[623] + v[343] * v[638] - v[350] * v[653] - v[360] * v[668];
//	v[1118] = v[332] * v[623] + v[342] * v[638] - v[349] * v[653] - v[359] * v[668];
//	v[661] = -(invH[3][0] * v[573]) - invH[3][1] * v[579] - invH[3][2] * v[585] - invH[3][3] * v[591];
//	v[1124] = v[334] * v[616] + v[344] * v[631] - v[351] * v[646] - v[361] * v[661];
//	v[1123] = v[333] * v[616] + v[343] * v[631] - v[350] * v[646] - v[360] * v[661];
//	v[1122] = v[332] * v[616] + v[342] * v[631] - v[349] * v[646] - v[359] * v[661];
//	v[662] = -(invH[3][0] * v[574]) - invH[3][1] * v[580] - invH[3][2] * v[586] - invH[3][3] * v[592];
//	v[1128] = v[334] * v[617] + v[344] * v[632] - v[351] * v[647] - v[361] * v[662];
//	v[1127] = v[333] * v[617] + v[343] * v[632] - v[350] * v[647] - v[360] * v[662];
//	v[1126] = v[332] * v[617] + v[342] * v[632] - v[349] * v[647] - v[359] * v[662];
//	v[663] = -(invH[3][0] * v[575]) - invH[3][1] * v[581] - invH[3][2] * v[587] - invH[3][3] * v[593];
//	v[1132] = v[334] * v[618] + v[344] * v[633] - v[351] * v[648] - v[361] * v[663];
//	v[1131] = v[333] * v[618] + v[343] * v[633] - v[350] * v[648] - v[360] * v[663];
//	v[1130] = v[332] * v[618] + v[342] * v[633] - v[349] * v[648] - v[359] * v[663];
//	v[670] = -(invH[3][0] * v[576]) - invH[3][1] * v[582] - invH[3][2] * v[588] - invH[3][3] * v[594];
//	v[1148] = v[334] * v[625] + v[344] * v[640] - v[351] * v[655] - v[361] * v[670];
//	v[1147] = v[333] * v[625] + v[343] * v[640] - v[350] * v[655] - v[360] * v[670];
//	v[1146] = v[332] * v[625] + v[342] * v[640] - v[349] * v[655] - v[359] * v[670];
//	v[671] = -(invH[3][0] * v[577]) - invH[3][1] * v[583] - invH[3][2] * v[589] - invH[3][3] * v[595];
//	v[1152] = v[334] * v[626] + v[344] * v[641] - v[351] * v[656] - v[361] * v[671];
//	v[1151] = v[333] * v[626] + v[343] * v[641] - v[350] * v[656] - v[360] * v[671];
//	v[1150] = v[332] * v[626] + v[342] * v[641] - v[349] * v[656] - v[359] * v[671];
//	v[672] = -(invH[3][0] * v[578]) - invH[3][1] * v[584] - invH[3][2] * v[590] - invH[3][3] * v[596];
//	v[1177] = v[255] * v[661] + v[265] * v[662] + v[275] * v[663] + v[3370] * v[664] + v[3354] * v[666] + v[3335] * v[668]
//		+ v[281] * v[670] + v[291] * v[671] + v[301] * v[672];
//	v[1156] = v[334] * v[627] + v[344] * v[642] - v[351] * v[657] - v[361] * v[672];
//	v[1155] = v[333] * v[627] + v[343] * v[642] - v[350] * v[657] - v[360] * v[672];
//	v[1154] = v[332] * v[627] + v[342] * v[642] - v[349] * v[657] - v[359] * v[672];
//	v[5394] = v[664];
//	v[5395] = v[666];
//	v[5396] = v[668];
//	v[5397] = v[661];
//	v[5398] = v[662];
//	v[5399] = v[663];
//	v[5400] = -v[664];
//	v[5401] = -v[666];
//	v[5402] = -v[668];
//	v[5403] = v[670];
//	v[5404] = v[671];
//	v[5405] = v[672];
//	b673 = sqrt(Power(v[318] * v[326] - v[317] * v[327], 2) + Power(-(v[319] * v[326]) + v[317] * v[328], 2) + Power
//	(v[319] * v[327] - v[318] * v[328], 2)) > 0.1e-7;
//	if (b673) {
//		v[675] = v[319] * v[327] - v[318] * v[328];
//		v[676] = -(v[319] * v[326]) + v[317] * v[328];
//		v[677] = v[318] * v[326] - v[317] * v[327];
//		v[678] = sqrt((v[675] * v[675]) + (v[676] * v[676]) + (v[677] * v[677]));
//		v[1942] = 1e0 / Power(v[678], 2);
//		v[1361] = v[678];
//		v[1953] = 1e0 - (v[1361] * v[1361]);
//		v[4225] = 1e0 / Power(v[1953], 0.15e1);
//		v[1948] = 1e0 / sqrt(v[1953]);
//		v[1360] = asin(v[1361]) / 2e0;
//		v[4224] = tan(v[1360]);
//		v[1947] = 1e0 / Power(cos(v[1360]), 2);
//		v[3855] = v[1947] * v[1948];
//		v[680] = 2e0*tan(v[1360]);
//		if (v[678] > 0.1e-7) { v07 = 1e0 / v[678]; v08 = (-(v07 / v[678])); v09 = (2e0*v07) / Power(v[678], 2); }
//		else {
//			v07 = (12500000e0 / 3e0)*(24e0 - (-0.1e-7 + v[678])*(0.24e10 - 2e0*(-1e0 + 100000000e0*v[678])*
//				(0.2399999997e10 - 0.1199999994e18*v[678] - 0.3e17*(v[678] * v[678]))));
//			v08 = (-50000000e0 / 3e0)*(0.3599999994e10 - 0.4799999982e18*v[678] + 0.6e25*Power(v[678], 3)
//				+ 0.1799999982e26*(v[678] * v[678]));
//			v09 = 0.1e17*(799999997e0 - 0.599999994e17*v[678] - 0.3e17*(v[678] * v[678]));
//		};
//		v[684] = v09;
//		v[685] = v08;
//		v[686] = v07;
//		v[4223] = v[680] * v[685] + v[3855] * v[686];
//		v[3759] = v[680] * v[686];
//		v[687] = v[3759] * v[675];
//		v[698] = (v[687] * v[687]);
//		v[688] = v[3759] * v[676];
//		v[696] = (v[687] * v[688]) / 2e0;
//		v[691] = (v[688] * v[688]);
//		v[1344] = -v[691] - v[698];
//		v[689] = v[3759] * v[677];
//		v[1339] = v[689] + v[696];
//		v[1337] = -v[689] + v[696];
//		v[703] = (v[688] * v[689]) / 2e0;
//		v[1343] = v[687] + v[703];
//		v[1341] = -v[687] + v[703];
//		v[701] = (v[687] * v[689]) / 2e0;
//		v[1342] = -v[688] + v[701];
//		v[1338] = v[688] + v[701];
//		v[692] = (v[689] * v[689]);
//		v[1348] = 4e0 + v[691] + v[692] + v[698];
//		v[4226] = 1e0 / Power(v[1348], 3);
//		v[2312] = 1e0 / Power(v[1348], 2);
//		v[1340] = -v[692] - v[698];
//		v[1336] = -v[691] - v[692];
//		v[690] = 4e0 / v[1348];
//		v[693] = 1e0 + (v[1336] * v[690]) / 2e0;
//		v[694] = v[1337] * v[690];
//		v[695] = v[1338] * v[690];
//		v[697] = v[1339] * v[690];
//		v[699] = 1e0 + (v[1340] * v[690]) / 2e0;
//		v[700] = v[1341] * v[690];
//		v[702] = v[1342] * v[690];
//		v[704] = v[1343] * v[690];
//		v[705] = 1e0 + (v[1344] * v[690]) / 2e0;
//	}
//	else {
//		v[693] = 1e0;
//		v[694] = 0e0;
//		v[695] = 0e0;
//		v[697] = 0e0;
//		v[699] = 1e0;
//		v[700] = 0e0;
//		v[702] = 0e0;
//		v[704] = 0e0;
//		v[705] = 1e0;
//	};
//	if ((*previouscontact)) {
//		v[1303] = 1e0 - v[717];
//		v[1301] = 1e0 - v[715];
//		v[1299] = 1e0 - v[713];
//		v[710] = v[145] * v[173] + v[146] * v[174] + v[147] * v[175] - v[201] * v[229] - v[202] * v[230] - v[203] * v[231] + v[3760]
//			+ gti[0] * v[702] + gti[1] * v[704] + gti[2] * v[705];
//		v[3763] = v[319] * v[710];
//		v[709] = v[145] * v[170] + v[146] * v[171] + v[147] * v[172] - v[201] * v[226] - v[202] * v[227] - v[203] * v[228] + v[3761]
//			+ gti[0] * v[697] + gti[1] * v[699] + gti[2] * v[700];
//		v[3765] = v[318] * v[709];
//		v[3789] = v[3763] + v[3765];
//		v[708] = v[145] * v[167] + v[146] * v[168] + v[147] * v[169] - v[201] * v[223] - v[202] * v[224] - v[203] * v[225] + v[3762]
//			+ gti[0] * v[693] + gti[1] * v[694] + gti[2] * v[695];
//		v[3764] = -(v[317] * v[708]);
//		v[3788] = -v[3763] + v[3764];
//		v[3787] = v[3764] - v[3765];
//		v[707] = -(v[317] * v[3789]) + v[1299] * v[708];
//		v[711] = v[318] * v[3788] + v[1301] * v[709];
//		v[712] = v[319] * v[3787] + v[1303] * v[710];
//	}
//	else {
//		v[707] = 0e0;
//		v[711] = 0e0;
//		v[712] = 0e0;
//	};
//	v[714] = v[1261] * v[255] + v[1257] * v[265] + v[1253] * v[275] + v[1249] * v[281] + v[1245] * v[291] + v[1241] * v[301]
//		+ v[2272] * v[317] + v[3370] * v[713];
//	v[716] = v[1262] * v[255] + v[1258] * v[265] + v[1254] * v[275] + v[1250] * v[281] + v[1246] * v[291] + v[1242] * v[301]
//		+ v[2268] * v[318] + v[3354] * v[715];
//	v[718] = v[1263] * v[255] + v[1259] * v[265] + v[1255] * v[275] + v[1251] * v[281] + v[1247] * v[291] + v[1243] * v[301]
//		+ v[2265] * v[319] + v[3335] * v[717];
//	(*vnrel) = sqrt((v[714] * v[714]) + (v[716] * v[716]) + (v[718] * v[718]));
//	v[724] = -((v[3767] * Power(v[725], v[740])) / ((*n2)*Power((*gnbb), v[745])));
//	b728 = v[727] < (*gnb);
//	if (b728) {
//		b730 = v[727] > (*gnbb);
//		if (b730) {
//			v[739] = (*gnb) - v[727];
//			v[4419] = Power(v[739], v[740]);
//			v[4411] = Power(v[739], v[2115]);
//			v[4300] = Power(v[739], v[740]);
//			v[4257] = Power(v[739], -2e0 + v[740]);
//			v[2225] = Power(v[739], v[2115]);
//			v[4184] = Power(v[739], v[740]);
//			v[733] = Power(v[739], (*n1));
//			v[4084] = -((*epsn1)*v[733]);
//			v[3917] = -((*epsn1)*v[733]);
//			v[3775] = -((*epsn1)*v[733]);
//			v[3766] = -((*epsn1)*v[733]);
//			v[732] = v[317] * v[3766];
//			v[734] = v[318] * v[3766];
//			v[735] = v[319] * v[3766];
//		}
//		else {
//			v[4423] = Power(v[727], v[745]);
//			v[4412] = Power(v[727], v[2120]);
//			v[4304] = Power(v[727], v[745]);
//			v[4262] = Power(v[727], -2e0 + v[745]);
//			v[2231] = Power(v[727], v[2120]);
//			v[4188] = Power(v[727], v[745]);
//			v[736] = -((*epsn1)*Power(v[725], (*n1))) + v[724] * (Power((*gnbb), (*n2)) - Power(v[727], (*n2)));
//			v[732] = v[317] * v[736];
//			v[734] = v[318] * v[736];
//			v[735] = v[319] * v[736];
//		};
//	}
//	else {
//		v[732] = 0e0;
//		v[734] = 0e0;
//		v[735] = 0e0;
//	};
//	if (b728) {
//		if (b730) {
//			v[1201] = (*meq)*v[3767];
//			v[1202] = sqrt(v[1201] * Power(v[739], v[740]));
//			v[2224] = (v[1201] * v[740]) / v[1202];
//			v[742] = 2e0*v[1202] * (*zetan);
//			v[741] = v[714] * v[742];
//			v[743] = v[716] * v[742];
//			v[744] = v[718] * v[742];
//		}
//		else {
//			v[1206] = -((*meq)*(*n2)*v[724]);
//			v[1207] = sqrt(v[1206] * Power(v[727], v[745]));
//			v[2230] = (v[1206] * v[745]) / v[1207];
//			v[746] = 2e0*v[1207] * (*zetan);
//			v[741] = v[714] * v[746];
//			v[743] = v[716] * v[746];
//			v[744] = v[718] * v[746];
//		};
//		b747 = v[727] < (*gnb) && v[317] * (v[732] + v[741]) + v[318] * (v[734] + v[743]) + v[319] * (v[735] + v[744]) < 0e0;
//		if (b747) {
//			v[749] = v[741];
//			v[750] = v[743];
//			v[751] = v[744];
//		}
//		else {
//			v[749] = -v[732];
//			v[750] = -v[734];
//			v[751] = -v[735];
//		};
//	}
//	else {
//		v[749] = 0e0;
//		v[750] = 0e0;
//		v[751] = 0e0;
//	};
//	v[752] = v[732] + v[749];
//	v[753] = v[734] + v[750];
//	v[754] = v[735] + v[751];
//	v[2154] = (v[752] * v[752]) + (v[753] * v[753]) + (v[754] * v[754]);
//	v[755] = (*epst)*v[707];
//	v[756] = (*epst)*v[711];
//	v[757] = (*epst)*v[712];
//	v[761] = -((*ct)*(v[1110] * v[244] + v[1114] * v[245] + v[1118] * v[246] - v[1110] * v[247] - v[1114] * v[248]
//		- v[1118] * v[249] + v[1122] * v[255] + v[1126] * v[265] + v[1130] * v[275] + v[1146] * v[281] + v[1150] * v[291]
//		+ v[1154] * v[301])) + v[755];
//	v[762] = -((*ct)*(v[1111] * v[244] + v[1115] * v[245] + v[1119] * v[246] - v[1111] * v[247] - v[1115] * v[248]
//		- v[1119] * v[249] + v[1123] * v[255] + v[1127] * v[265] + v[1131] * v[275] + v[1147] * v[281] + v[1151] * v[291]
//		+ v[1155] * v[301])) + v[756];
//	v[763] = -((*ct)*(v[1112] * v[244] + v[1116] * v[245] + v[1120] * v[246] - v[1112] * v[247] - v[1116] * v[248]
//		- v[1120] * v[249] + v[1124] * v[255] + v[1128] * v[265] + v[1132] * v[275] + v[1148] * v[281] + v[1152] * v[291]
//		+ v[1156] * v[301])) + v[757];
//	v[2151] = (v[761] * v[761]) + (v[762] * v[762]) + (v[763] * v[763]);
//	if (b728) {
//		if ((*stick)) {
//			b766 = sqrt((v[761] * v[761]) + (v[762] * v[762]) + (v[763] * v[763])) <= (*mus)*sqrt((v[752] * v[752]) +
//				(v[753] * v[753]) + (v[754] * v[754]));
//			if (b766) {
//				v[768] = v[761];
//				v[769] = v[762];
//				v[770] = v[763];
//				v[771] = 1e0;
//			}
//			else {
//				v[2153] = 1e0 / sqrt(v[2151]);
//				v[2155] = 1e0 / sqrt(v[2154]);
//				v[781] = sqrt(v[2154]);
//				v[772] = sqrt(v[2151]);
//				if (v[772] > 0.1e-5) { v010 = 1e0 / v[772]; v011 = (-(v010 / v[772])); v012 = (2e0*v010) / Power(v[772], 2); }
//				else {
//					v010 = (24000000e0 - (-1e0 + 1000000e0*v[772])*(71999994e0 - 0.71999982e14*v[772] + 0.6e19*Power(v[772], 3)
//						+ 0.23999982e20*(v[772] * v[772]))) / 24e0;
//					v011 = (-500000e0 / 3e0)*(35999994e0 - 0.47999982e14*v[772] + 0.6e19*Power(v[772], 3) + 0.17999982e20*
//						(v[772] * v[772]));
//					v012 = 0.1e13*(7999997e0 - 0.5999994e13*v[772] - 0.3e13*(v[772] * v[772]));
//				};
//				v[776] = v011;
//				v[777] = v010;
//				v[2152] = (*mud)*v[777] * v[781];
//				v[768] = v[2152] * v[761];
//				v[769] = v[2152] * v[762];
//				v[770] = v[2152] * v[763];
//				v[771] = 0e0;
//			};
//			if (sqrt((v[755] * v[755]) + (v[756] * v[756]) + (v[757] * v[757])) > (*mus)*sqrt((v[752] * v[752]) +
//				(v[753] * v[753]) + (v[754] * v[754]))) {
//				if ((*epst) > 0.1e-5) {
//					v013 = 1e0 / (*epst); v014 = (-(v013 / (*epst))); v015 = (2e0*v013) / Power((*epst), 2
//					);
//				}
//				else {
//					v013 = (24000000e0 - (-1e0 + 1000000e0*(*epst))*(71999994e0 - 0.71999982e14*(*epst) + 0.23999982e20*Power(
//						(*epst), 2) + 0.6e19*Power((*epst), 3))) / 24e0;
//					v014 = (-500000e0 / 3e0)*(35999994e0 - 0.47999982e14*(*epst) + 0.17999982e20*Power((*epst), 2)
//						+ 0.6e19*Power((*epst), 3));
//					v015 = 0.1e13*(7999997e0 - 0.5999994e13*(*epst) - 0.3e13*Power((*epst), 2));
//				};
//				v[790] = sqrt((v[755] * v[755]) + (v[756] * v[756]) + (v[757] * v[757]));
//				if (v[790] > 0.1e-5) { v016 = 1e0 / v[790]; v017 = (-(v016 / v[790])); v018 = (2e0*v016) / Power(v[790], 2); }
//				else {
//					v016 = (24000000e0 - (-1e0 + 1000000e0*v[790])*(71999994e0 - 0.71999982e14*v[790] + 0.6e19*Power(v[790], 3)
//						+ 0.23999982e20*(v[790] * v[790]))) / 24e0;
//					v017 = (-500000e0 / 3e0)*(35999994e0 - 0.47999982e14*v[790] + 0.6e19*Power(v[790], 3) + 0.17999982e20*
//						(v[790] * v[790]));
//					v018 = 0.1e13*(7999997e0 - 0.5999994e13*v[790] - 0.3e13*(v[790] * v[790]));
//				};
//				v[797] = -((*mud)*v013*v016*sqrt(v[2154]));
//				v[796] = v[707] + v[755] * v[797];
//				v[798] = v[711] + v[756] * v[797];
//				v[799] = v[712] + v[757] * v[797];
//			}
//			else {
//				v[796] = 0e0;
//				v[798] = 0e0;
//				v[799] = 0e0;
//			};
//		}
//		else {
//			b800 = sqrt((v[761] * v[761]) + (v[762] * v[762]) + (v[763] * v[763])) <= (*mud)*sqrt((v[752] * v[752]) +
//				(v[753] * v[753]) + (v[754] * v[754]));
//			if (b800) {
//				v[768] = v[761];
//				v[769] = v[762];
//				v[770] = v[763];
//				v[771] = 1e0;
//			}
//			else {
//				v[2161] = 1e0 / sqrt(v[2151]);
//				v[2165] = 1e0 / sqrt(v[2154]);
//				v[811] = sqrt(v[2154]);
//				v[3905] = (*mud)*v[811];
//				v[802] = sqrt(v[2151]);
//				if (v[802] > 0.1e-5) { v019 = 1e0 / v[802]; v020 = (-(v019 / v[802])); v021 = (2e0*v019) / Power(v[802], 2); }
//				else {
//					v019 = (24000000e0 - (-1e0 + 1000000e0*v[802])*(71999994e0 - 0.71999982e14*v[802] + 0.6e19*Power(v[802], 3)
//						+ 0.23999982e20*(v[802] * v[802]))) / 24e0;
//					v020 = (-500000e0 / 3e0)*(35999994e0 - 0.47999982e14*v[802] + 0.6e19*Power(v[802], 3) + 0.17999982e20*
//						(v[802] * v[802]));
//					v021 = 0.1e13*(7999997e0 - 0.5999994e13*v[802] - 0.3e13*(v[802] * v[802]));
//				};
//				v[806] = v020;
//				v[807] = v019;
//				v[3768] = (*mud)*v[807] * v[811];
//				v[768] = v[3768] * v[761];
//				v[769] = v[3768] * v[762];
//				v[770] = v[3768] * v[763];
//				v[771] = 0e0;
//			};
//			if (sqrt((v[755] * v[755]) + (v[756] * v[756]) + (v[757] * v[757])) > (*mud)*sqrt((v[752] * v[752]) +
//				(v[753] * v[753]) + (v[754] * v[754]))) {
//				if ((*epst) > 0.1e-5) {
//					v022 = 1e0 / (*epst); v023 = (-(v022 / (*epst))); v024 = (2e0*v022) / Power((*epst), 2
//					);
//				}
//				else {
//					v022 = (24000000e0 - (-1e0 + 1000000e0*(*epst))*(71999994e0 - 0.71999982e14*(*epst) + 0.23999982e20*Power(
//						(*epst), 2) + 0.6e19*Power((*epst), 3))) / 24e0;
//					v023 = (-500000e0 / 3e0)*(35999994e0 - 0.47999982e14*(*epst) + 0.17999982e20*Power((*epst), 2)
//						+ 0.6e19*Power((*epst), 3));
//					v024 = 0.1e13*(7999997e0 - 0.5999994e13*(*epst) - 0.3e13*Power((*epst), 2));
//				};
//				v[820] = sqrt((v[755] * v[755]) + (v[756] * v[756]) + (v[757] * v[757]));
//				if (v[820] > 0.1e-5) { v025 = 1e0 / v[820]; v026 = (-(v025 / v[820])); v027 = (2e0*v025) / Power(v[820], 2); }
//				else {
//					v025 = (24000000e0 - (-1e0 + 1000000e0*v[820])*(71999994e0 - 0.71999982e14*v[820] + 0.6e19*Power(v[820], 3)
//						+ 0.23999982e20*(v[820] * v[820]))) / 24e0;
//					v026 = (-500000e0 / 3e0)*(35999994e0 - 0.47999982e14*v[820] + 0.6e19*Power(v[820], 3) + 0.17999982e20*
//						(v[820] * v[820]));
//					v027 = 0.1e13*(7999997e0 - 0.5999994e13*v[820] - 0.3e13*(v[820] * v[820]));
//				};
//				v[826] = -((*mud)*v022*v025*sqrt(v[2154]));
//				v[796] = v[707] + v[755] * v[826];
//				v[798] = v[711] + v[756] * v[826];
//				v[799] = v[712] + v[757] * v[826];
//			}
//			else {
//				v[796] = 0e0;
//				v[798] = 0e0;
//				v[799] = 0e0;
//			};
//		};
//	}
//	else {
//		v[768] = 0e0;
//		v[769] = 0e0;
//		v[770] = 0e0;
//	};
//	fn[0] = v[752];
//	fn[1] = v[753];
//	fn[2] = v[754];
//	ft[0] = v[768];
//	ft[1] = v[769];
//	ft[2] = v[770];
//	(*stickupdated) = v[771];
//	gtpupdated[0] = v[707] - v[796];
//	gtpupdated[1] = v[711] - v[798];
//	gtpupdated[2] = v[712] - v[799];
//	b840 = b728;
//	if (b840) {
//		b841 = b730;
//	}
//	else {
//	};
//	v[845] = v[735] * v[851];
//	v[846] = v[735] * v[853];
//	v[900] = v[204] * v[846];
//	v[847] = v[735] * v[855];
//	v[906] = v[204] * v[847];
//	v[848] = v[735] * v[857];
//	v[849] = v[735] * v[859];
//	v[889] = v[148] * v[849];
//	v[850] = v[735] * v[861];
//	v[895] = v[148] * v[850];
//	v[852] = v[734] * v[851];
//	v[901] = v[204] * v[852];
//	v[854] = v[734] * v[853];
//	v[903] = -(v[204] * v[854]) / 2e0;
//	v[856] = v[734] * v[855];
//	v[909] = v[204] * v[856];
//	v[858] = v[734] * v[857];
//	v[890] = v[148] * v[858];
//	v[860] = v[734] * v[859];
//	v[892] = -(v[148] * v[860]) / 2e0;
//	v[862] = v[734] * v[861];
//	v[898] = v[148] * v[862];
//	v[863] = v[732] * v[851];
//	v[907] = v[204] * v[863];
//	v[864] = v[732] * v[853];
//	v[910] = v[204] * v[864];
//	v[865] = v[732] * v[855];
//	v[908] = -(v[204] * v[865]) / 2e0;
//	v[866] = -(v[225] * v[732]) - v[228] * v[734] - v[231] * v[735];
//	v[867] = -(v[224] * v[732]) - v[227] * v[734] - v[230] * v[735];
//	v[868] = -(v[223] * v[732]) - v[226] * v[734] - v[229] * v[735];
//	v[869] = v[732] * v[857];
//	v[896] = v[148] * v[869];
//	v[870] = v[732] * v[859];
//	v[899] = v[148] * v[870];
//	v[871] = v[732] * v[861];
//	v[897] = -(v[148] * v[871]) / 2e0;
//	v[872] = v[169] * v[732] + v[172] * v[734] + v[175] * v[735];
//	v[873] = v[168] * v[732] + v[171] * v[734] + v[174] * v[735];
//	v[874] = v[167] * v[732] + v[170] * v[734] + v[173] * v[735];
//	v[875] = v[900] + v[901];
//	v[1100] = v[875] / 2e0;
//	v[1096] = (v[906] + v[907]) / 2e0;
//	v[877] = v[909] + v[910];
//	v[1095] = v[877] / 2e0;
//	v[878] = (v[514] * v[845]) / 2e0 + v[509] * v[846] + v[505] * v[847] + v[501] * v[852] + (v[496] * v[854]) / 2e0
//		+ v[492] * v[856] + v[488] * v[863] + v[482] * v[864] + (v[477] * v[865]) / 2e0;
//	v[1103] = v[903] + v[908] - 4e0*v[878] * v[968];
//	v[1099] = v[1103] - (v[204] * v[845]) / 2e0 - v[903];
//	v[1094] = v[1099] + v[903] - v[908];
//	v[879] = v[868] * x1B[0] + v[867] * x1B[1] + v[866] * x1B[2];
//	v[880] = v[889] + v[890];
//	v[1090] = v[880] / 2e0;
//	v[1086] = (v[895] + v[896]) / 2e0;
//	v[882] = v[898] + v[899];
//	v[1085] = v[882] / 2e0;
//	v[883] = (v[415] * v[848]) / 2e0 + v[410] * v[849] + v[406] * v[850] + v[402] * v[858] + (v[397] * v[860]) / 2e0
//		+ v[393] * v[862] + v[389] * v[869] + v[383] * v[870] + (v[378] * v[871]) / 2e0;
//	v[1093] = v[892] + v[897] - 4e0*v[883] * v[962];
//	v[1089] = v[1093] - (v[148] * v[848]) / 2e0 - v[892];
//	v[1084] = v[1089] + v[892] - v[897];
//	v[5346] = v[732];
//	v[5347] = v[734];
//	v[5348] = v[735];
//	v[5349] = 2e0*alphaA[0] * v[1084] + alphaA[2] * v[1086] + v[364] * v[882] + v[889] - v[890];
//	v[5350] = 2e0*alphaA[1] * v[1089] + (alphaA[2] * v[880]) / 2e0 + (alphaA[0] * v[882]) / 2e0 - v[895] + v[896];
//	v[5351] = alphaA[0] * v[1086] + 2e0*alphaA[2] * v[1093] + v[364] * v[880] + v[898] - v[899];
//	v[5352] = -v[732];
//	v[5353] = -v[734];
//	v[5354] = -v[735];
//	v[5355] = 2e0*alphaB[0] * v[1094] + alphaB[2] * v[1096] + v[370] * v[877] + v[900] - v[901];
//	v[5356] = 2e0*alphaB[1] * v[1099] + (alphaB[2] * v[875]) / 2e0 + (alphaB[0] * v[877]) / 2e0 - v[906] + v[907];
//	v[5357] = alphaB[0] * v[1096] + 2e0*alphaB[2] * v[1103] + v[370] * v[875] + v[909] - v[910];
//	v[884] = v[874] * x1A[0] + v[873] * x1A[1] + v[872] * x1A[2];
//	v[885] = (-v[879] + v[868] * x3B[0] + v[867] * x3B[1] + v[866] * x3B[2]) / 2e0;
//	v[886] = (-v[879] + v[868] * x2B[0] + v[867] * x2B[1] + v[866] * x2B[2]) / 2e0;
//	v[887] = (-v[884] + v[874] * x2A[0] + v[873] * x2A[1] + v[872] * x2A[2]) / 2e0;
//	v[888] = (-v[884] + v[874] * x3A[0] + v[873] * x3A[1] + v[872] * x3A[2]) / 2e0;
//	for (i835 = 1; i835 <= 12; i835++) {
//		i3774 = (i835 == 10 ? 1 : 0);
//		i3773 = (i835 == 11 ? 1 : 0);
//		i3772 = (i835 == 12 ? 1 : 0);
//		i3771 = (i835 == 4 ? 1 : 0);
//		i3770 = (i835 == 5 ? 1 : 0);
//		i3769 = (i835 == 6 ? 1 : 0);
//		v[924] = v[5393 + i835];
//		v[923] = v[5381 + i835];
//		v[917] = v[5357 + i835];
//		v[916] = v[5369 + i835];
//		v[918] = (-v[916] - v[917]) / 2e0;
//		v[919] = v[5421 + i835];
//		v[964] = -4e0*v[919] * v[962];
//		v[920] = v[5433 + i835];
//		v[921] = v[5445 + i835];
//		v[922] = v[5457 + i835];
//		v[925] = (-v[923] - v[924]) / 2e0;
//		v[926] = v[5469 + i835];
//		v[970] = -4e0*v[926] * v[968];
//		v[927] = v[5481 + i835];
//		v[928] = v[5493 + i835];
//		v[929] = v[5505 + i835];
//		v[930] = v[5517 + i835];
//		v[931] = v[5553 + i835];
//		v[932] = v[5601 + i835];
//		v[933] = v[5637 + i835];
//		v[934] = v[5685 + i835];
//		v[935] = v[5721 + i835];
//		v[936] = (2e0*v[918] * x1A[0] + v[916] * x2A[0] + v[917] * x3A[0]) / 2e0;
//		v[937] = (2e0*v[918] * x1A[1] + v[916] * x2A[1] + v[917] * x3A[1]) / 2e0;
//		v[938] = (2e0*v[918] * x1A[2] + v[916] * x2A[2] + v[917] * x3A[2]) / 2e0;
//		v[939] = -i3769 + v[920];
//		v[941] = i3769 + v[920];
//		v[942] = i3770 + v[921];
//		v[944] = -i3770 + v[921];
//		v[945] = -i3771 + v[922];
//		v[947] = i3771 + v[922];
//		v[948] = (2e0*v[925] * x1B[0] + v[923] * x2B[0] + v[924] * x3B[0]) / 2e0;
//		v[949] = (2e0*v[925] * x1B[1] + v[923] * x2B[1] + v[924] * x3B[1]) / 2e0;
//		v[950] = (2e0*v[925] * x1B[2] + v[923] * x2B[2] + v[924] * x3B[2]) / 2e0;
//		v[951] = -i3772 + v[927];
//		v[953] = i3772 + v[927];
//		v[954] = i3773 + v[928];
//		v[956] = -i3773 + v[928];
//		v[957] = -i3774 + v[929];
//		v[959] = i3774 + v[929];
//		v[961] = -(v[148] * v[930]) / 2e0 + (v[378] * v[964]) / 2e0;
//		v[963] = v[148] * v[939] + v[383] * v[964];
//		v[965] = v[148] * v[942] + v[389] * v[964];
//		v[967] = -(v[204] * v[931]) / 2e0 + (v[477] * v[970]) / 2e0;
//		v[969] = v[204] * v[951] + v[482] * v[970];
//		v[971] = v[204] * v[954] + v[488] * v[970];
//		v[972] = v[5733 + i835] + v[167] * v[936] + v[168] * v[937] + v[169] * v[938] - v[223] * v[948] - v[224] * v[949]
//			- v[225] * v[950] + v[861] * v[961] + v[859] * v[963] + v[857] * v[965] + v[855] * v[967] + v[853] * v[969]
//			+ v[851] * v[971];
//		v[1001] = v[972];
//		v[973] = v[148] * v[941] + v[393] * v[964];
//		v[974] = (-(v[148] * v[932]) + v[397] * v[964]) / 2e0;
//		v[975] = v[148] * v[945] + v[402] * v[964];
//		v[976] = v[204] * v[953] + v[492] * v[970];
//		v[977] = (-(v[204] * v[933]) + v[496] * v[970]) / 2e0;
//		v[978] = v[204] * v[957] + v[501] * v[970];
//		v[979] = v[5745 + i835] + v[170] * v[936] + v[171] * v[937] + v[172] * v[938] - v[226] * v[948] - v[227] * v[949]
//			- v[228] * v[950] + v[861] * v[973] + v[859] * v[974] + v[857] * v[975] + v[855] * v[976] + v[853] * v[977]
//			+ v[851] * v[978];
//		v[1000] = v[979];
//		v[980] = v[148] * v[944] + v[406] * v[964];
//		v[981] = v[732] * v[961] + v[734] * v[973] + v[735] * v[980];
//		v[982] = v[148] * v[947] + v[410] * v[964];
//		v[983] = v[732] * v[963] + v[734] * v[974] + v[735] * v[982];
//		v[984] = (-(v[148] * v[934]) + v[415] * v[964]) / 2e0;
//		v[985] = v[732] * v[965] + v[734] * v[975] + v[735] * v[984];
//		v[986] = v[204] * v[956] + v[505] * v[970];
//		v[987] = v[732] * v[967] + v[734] * v[976] + v[735] * v[986];
//		v[988] = v[204] * v[959] + v[509] * v[970];
//		v[989] = v[732] * v[969] + v[734] * v[977] + v[735] * v[988];
//		v[990] = (-(v[204] * v[935]) + v[514] * v[970]) / 2e0;
//		v[991] = v[5757 + i835] + v[173] * v[936] + v[174] * v[937] + v[175] * v[938] - v[229] * v[948] - v[230] * v[949]
//			- v[231] * v[950] + v[861] * v[980] + v[859] * v[982] + v[857] * v[984] + v[855] * v[986] + v[853] * v[988]
//			+ v[851] * v[990];
//		v[999] = v[991];
//		v[992] = v[732] * v[971] + v[734] * v[978] + v[735] * v[990];
//		v[993] = 0e0;
//		v[994] = 0e0;
//		v[995] = 0e0;
//		v[996] = 0e0;
//		b997 = b728;
//		if (b997) {
//			v[3776] = -(v[1001] * v[317]) - v[1000] * v[318] - v[319] * v[999];
//			b998 = b730;
//			if (b998) {
//				v[995] = v[3775] * v[991];
//				v[991] = 0e0;
//				v[994] = v[3775] * v[979];
//				v[979] = 0e0;
//				v[993] = v[3775] * v[972];
//				v[972] = 0e0;
//				v[996] = -(v[3767] * v[3776] * v[4184]);
//			}
//			else {
//				v[995] = v[736] * v[999];
//				v[991] = 0e0;
//				v[994] = v[1000] * v[736];
//				v[979] = 0e0;
//				v[993] = v[1001] * v[736];
//				v[972] = 0e0;
//				v[996] = (*n2)*v[3776] * v[4188] * v[724];
//			};
//		}
//		else {
//		};
//		v[996] = v[315] * (v[302] * v[993] + v[303] * v[994] + v[304] * v[995]) + v[996];
//		v[3777] = v[996] / v[727];
//		v[1006] = v[304] * v[3777] + v[316] * v[995];
//		v[1008] = v[303] * v[3777] + v[316] * v[994];
//		v[1009] = v[302] * v[3777] + v[316] * v[993];
//		v[1010] = -(v[1006] * v[200]) - v[735] * v[950];
//		v[1011] = -(v[1006] * v[199]) - v[735] * v[949];
//		v[1012] = -(v[1006] * v[198]) - v[735] * v[948];
//		v[1013] = v[1006] * v[144] + v[735] * v[938];
//		v[1014] = v[1006] * v[143] + v[735] * v[937];
//		v[1015] = v[1006] * v[142] + v[735] * v[936];
//		v[1016] = -(v[1008] * v[200]) - v[734] * v[950];
//		v[1017] = -(v[1008] * v[199]) - v[734] * v[949];
//		v[1018] = -(v[1008] * v[198]) - v[734] * v[948];
//		v[1019] = v[1008] * v[144] + v[734] * v[938];
//		v[1020] = v[1008] * v[143] + v[734] * v[937];
//		v[1021] = v[1008] * v[142] + v[734] * v[936];
//		v[1022] = -(v[1009] * v[200]) - v[732] * v[950];
//		v[1023] = -(v[1009] * v[199]) - v[732] * v[949];
//		v[1024] = -(v[1009] * v[198]) - v[732] * v[948];
//		v[1025] = v[1009] * v[144] + v[732] * v[938];
//		v[1026] = v[1009] * v[143] + v[732] * v[937];
//		v[1027] = v[1009] * v[142] + v[732] * v[936];
//		v[1028] = QBi[2][2] * v[1010] + QBi[2][1] * v[1011] + QBi[2][0] * v[1012];
//		v[1029] = QBi[1][2] * v[1010] + QBi[1][1] * v[1011] + QBi[1][0] * v[1012];
//		v[1030] = QBi[0][2] * v[1010] + QBi[0][1] * v[1011] + QBi[0][0] * v[1012];
//		v[1031] = QBi[2][2] * v[1016] + QBi[2][1] * v[1017] + QBi[2][0] * v[1018];
//		v[1032] = QBi[1][2] * v[1016] + QBi[1][1] * v[1017] + QBi[1][0] * v[1018];
//		v[1033] = QBi[0][2] * v[1016] + QBi[0][1] * v[1017] + QBi[0][0] * v[1018];
//		v[1034] = QBi[2][2] * v[1022] + QBi[2][1] * v[1023] + QBi[2][0] * v[1024];
//		v[1035] = QBi[1][2] * v[1022] + QBi[1][1] * v[1023] + QBi[1][0] * v[1024];
//		v[1036] = QBi[0][2] * v[1022] + QBi[0][1] * v[1023] + QBi[0][0] * v[1024];
//		v[1037] = (v[1028] * v[204] + v[845] * v[970]) / 2e0;
//		v[1038] = v[1029] * v[204] + v[846] * v[970];
//		v[1039] = v[1030] * v[204] + v[847] * v[970];
//		v[1040] = v[1031] * v[204] + v[852] * v[970];
//		v[1042] = v[1033] * v[204] + v[856] * v[970];
//		v[1043] = v[1034] * v[204] + v[863] * v[970];
//		v[1044] = v[1035] * v[204] + v[864] * v[970];
//		v[3781] = 8e0*v[1564] * v[878] * v[926] - 4e0*((v[1036] * v[477]) / 2e0 + v[1035] * v[482] + v[1034] * v[488]
//			+ v[1033] * v[492] + (v[1032] * v[496]) / 2e0 + v[1031] * v[501] + v[1030] * v[505] + v[1029] * v[509] + (v[1028] * v[514])
//			/ 2e0 - (v[865] * v[931]) / 2e0 - (v[854] * v[933]) / 2e0 - (v[845] * v[935]) / 2e0 + v[864] * v[951] + v[856] * v[953]
//			+ v[863] * v[954] + v[847] * v[956] + v[852] * v[957] + v[846] * v[959])*v[968];
//		v[1098] = -(v[1032] * v[204]) / 2e0 + v[3781] - (v[854] * v[970]) / 2e0;
//		v[1046] = (v[1036] * v[204] + v[865] * v[970]) / 2e0;
//		v[1047] = -(v[1009] * v[225]) - v[1008] * v[228] - v[1006] * v[231] - QBi[0][2] * v[987] - QBi[1][2] * v[989]
//			- QBi[2][2] * v[992];
//		v[1048] = -(v[1009] * v[224]) - v[1008] * v[227] - v[1006] * v[230] - QBi[0][1] * v[987] - QBi[1][1] * v[989]
//			- QBi[2][1] * v[992];
//		v[1049] = -(v[1009] * v[223]) - v[1008] * v[226] - v[1006] * v[229] - QBi[0][0] * v[987] - QBi[1][0] * v[989]
//			- QBi[2][0] * v[992];
//		v[1050] = v[1049] * x1B[0] + v[1048] * x1B[1] + v[1047] * x1B[2];
//		v[1051] = QAi[2][2] * v[1013] + QAi[2][1] * v[1014] + QAi[2][0] * v[1015];
//		v[1052] = QAi[1][2] * v[1013] + QAi[1][1] * v[1014] + QAi[1][0] * v[1015];
//		v[1053] = QAi[0][2] * v[1013] + QAi[0][1] * v[1014] + QAi[0][0] * v[1015];
//		v[1054] = QAi[2][2] * v[1019] + QAi[2][1] * v[1020] + QAi[2][0] * v[1021];
//		v[1055] = QAi[1][2] * v[1019] + QAi[1][1] * v[1020] + QAi[1][0] * v[1021];
//		v[1056] = QAi[0][2] * v[1019] + QAi[0][1] * v[1020] + QAi[0][0] * v[1021];
//		v[1057] = QAi[2][2] * v[1025] + QAi[2][1] * v[1026] + QAi[2][0] * v[1027];
//		v[1058] = QAi[1][2] * v[1025] + QAi[1][1] * v[1026] + QAi[1][0] * v[1027];
//		v[1059] = QAi[0][2] * v[1025] + QAi[0][1] * v[1026] + QAi[0][0] * v[1027];
//		v[1060] = (v[1051] * v[148] + v[848] * v[964]) / 2e0;
//		v[1061] = v[1052] * v[148] + v[849] * v[964];
//		v[1062] = v[1053] * v[148] + v[850] * v[964];
//		v[1063] = v[1054] * v[148] + v[858] * v[964];
//		v[1065] = v[1056] * v[148] + v[862] * v[964];
//		v[1066] = v[1057] * v[148] + v[869] * v[964];
//		v[1067] = v[1058] * v[148] + v[870] * v[964];
//		v[3779] = 8e0*v[1562] * v[883] * v[919] - 4e0*((v[1059] * v[378]) / 2e0 + v[1058] * v[383] + v[1057] * v[389]
//			+ v[1056] * v[393] + (v[1055] * v[397]) / 2e0 + v[1054] * v[402] + v[1053] * v[406] + v[1052] * v[410] + (v[1051] * v[415])
//			/ 2e0 - (v[871] * v[930]) / 2e0 - (v[860] * v[932]) / 2e0 - (v[848] * v[934]) / 2e0 + v[870] * v[939] + v[862] * v[941]
//			+ v[869] * v[942] + v[850] * v[944] + v[858] * v[945] + v[849] * v[947])*v[962];
//		v[1088] = -(v[1055] * v[148]) / 2e0 + v[3779] - (v[860] * v[964]) / 2e0;
//		v[1069] = (v[1059] * v[148] + v[871] * v[964]) / 2e0;
//		v[1070] = v[1009] * v[169] + v[1008] * v[172] + v[1006] * v[175] + QAi[0][2] * v[981] + QAi[1][2] * v[983]
//			+ QAi[2][2] * v[985];
//		v[1071] = v[1009] * v[168] + v[1008] * v[171] + v[1006] * v[174] + QAi[0][1] * v[981] + QAi[1][1] * v[983]
//			+ QAi[2][1] * v[985];
//		v[1072] = v[1009] * v[167] + v[1008] * v[170] + v[1006] * v[173] + QAi[0][0] * v[981] + QAi[1][0] * v[983]
//			+ QAi[2][0] * v[985];
//		v[1073] = v[1072] * x1A[0] + v[1071] * x1A[1] + v[1070] * x1A[2];
//		v[1074] = (-v[1050] + v[1049] * x3B[0] + v[1048] * x3B[1] + v[1047] * x3B[2]) / 2e0;
//		v[1075] = (-v[1050] + v[1049] * x2B[0] + v[1048] * x2B[1] + v[1047] * x2B[2]) / 2e0;
//		v[3780] = (v[1039] + v[1043]) / 2e0;
//		v[1077] = v[1038] + v[1040];
//		v[1078] = v[1042] + v[1044];
//		v[1079] = (-v[1073] + v[1072] * x2A[0] + v[1071] * x2A[1] + v[1070] * x2A[2]) / 2e0;
//		v[1080] = (-v[1073] + v[1072] * x3A[0] + v[1071] * x3A[1] + v[1070] * x3A[2]) / 2e0;
//		v[3778] = (v[1062] + v[1066]) / 2e0;
//		v[1082] = v[1061] + v[1063];
//		v[1083] = v[1065] + v[1067];
//		v[5830] = 0e0;
//		v[5831] = 0e0;
//		v[5832] = 0e0;
//		v[5833] = 2e0*v[1084];
//		v[5834] = v[1085];
//		v[5835] = v[1086];
//		v[5836] = 0e0;
//		v[5837] = 0e0;
//		v[5838] = 0e0;
//		v[5839] = 0e0;
//		v[5840] = 0e0;
//		v[5841] = 0e0;
//		v[5818] = 0e0;
//		v[5819] = 0e0;
//		v[5820] = 0e0;
//		v[5821] = v[1085];
//		v[5822] = 2e0*v[1089];
//		v[5823] = v[1090];
//		v[5824] = 0e0;
//		v[5825] = 0e0;
//		v[5826] = 0e0;
//		v[5827] = 0e0;
//		v[5828] = 0e0;
//		v[5829] = 0e0;
//		v[5806] = 0e0;
//		v[5807] = 0e0;
//		v[5808] = 0e0;
//		v[5809] = v[1086];
//		v[5810] = v[1090];
//		v[5811] = 2e0*v[1093];
//		v[5812] = 0e0;
//		v[5813] = 0e0;
//		v[5814] = 0e0;
//		v[5815] = 0e0;
//		v[5816] = 0e0;
//		v[5817] = 0e0;
//		v[5770] = 0e0;
//		v[5771] = 0e0;
//		v[5772] = 0e0;
//		v[5773] = 0e0;
//		v[5774] = 0e0;
//		v[5775] = 0e0;
//		v[5776] = 0e0;
//		v[5777] = 0e0;
//		v[5778] = 0e0;
//		v[5779] = 2e0*v[1094];
//		v[5780] = v[1095];
//		v[5781] = v[1096];
//		v[5782] = 0e0;
//		v[5783] = 0e0;
//		v[5784] = 0e0;
//		v[5785] = 0e0;
//		v[5786] = 0e0;
//		v[5787] = 0e0;
//		v[5788] = 0e0;
//		v[5789] = 0e0;
//		v[5790] = 0e0;
//		v[5791] = v[1095];
//		v[5792] = 2e0*v[1099];
//		v[5793] = v[1100];
//		v[5794] = 0e0;
//		v[5795] = 0e0;
//		v[5796] = 0e0;
//		v[5797] = 0e0;
//		v[5798] = 0e0;
//		v[5799] = 0e0;
//		v[5800] = 0e0;
//		v[5801] = 0e0;
//		v[5802] = 0e0;
//		v[5803] = v[1096];
//		v[5804] = v[1100];
//		v[5805] = 2e0*v[1103];
//		v[5842] = v[1009];
//		v[5843] = v[1008];
//		v[5844] = v[1006];
//		v[5845] = v[1061] - v[1063] + 2e0*alphaA[0] * (-v[1060] + v[1088]) + v[1083] * v[364] + alphaA[2] * v[3778] + v[5829
//			+ i835];
//		v[5846] = -v[1062] + v[1066] + (alphaA[2] * v[1082]) / 2e0 + (alphaA[0] * v[1083]) / 2e0 + 2e0*alphaA[1] * (-v[1060]
//			- v[1069] + v[3779]) + v[5817 + i835];
//		v[5847] = v[1065] - v[1067] + 2e0*alphaA[2] * (-v[1069] + v[1088]) + v[1082] * v[364] + alphaA[0] * v[3778] + v[5805
//			+ i835];
//		v[5848] = -v[1009];
//		v[5849] = -v[1008];
//		v[5850] = -v[1006];
//		v[5851] = v[1038] - v[1040] + 2e0*alphaB[0] * (-v[1037] + v[1098]) + v[1078] * v[370] + alphaB[2] * v[3780] + v[5769
//			+ i835];
//		v[5852] = -v[1039] + v[1043] + (alphaB[2] * v[1077]) / 2e0 + (alphaB[0] * v[1078]) / 2e0 + 2e0*alphaB[1] * (-v[1037]
//			- v[1046] + v[3781]) + v[5781 + i835];
//		v[5853] = v[1042] - v[1044] + 2e0*alphaB[2] * (-v[1046] + v[1098]) + v[1077] * v[370] + alphaB[0] * v[3780] + v[5793
//			+ i835];
//		Rc[i835 - 1] += v[5345 + i835] + v[887] * v[916] + v[888] * v[917] + v[886] * v[923] + v[885] * v[924];
//		for (i913 = 1; i913 <= 12; i913++) {
//			Kc[i835 - 1][i913 - 1] += v[1080] * v[5357 + i913] + v[1079] * v[5369 + i913] + v[1075] * v[5381 + i913] + v[1074] * v[5393
//				+ i913] + v[5841 + i913];
//		};/* end for */
//	};/* end for */
//	v[1113] = -(v[1110] * v[768]) - v[1111] * v[769] - v[1112] * v[770];
//	v[1117] = -(v[1114] * v[768]) - v[1115] * v[769] - v[1116] * v[770];
//	v[1121] = -(v[1118] * v[768]) - v[1119] * v[769] - v[1120] * v[770];
//	v[1161] = 0e0;
//	v[1162] = 0e0;
//	v[1163] = 0e0;
//	b1164 = b728;
//	if (b1164) {
//		b1165 = (*stick);
//		if (b1165) {
//			b1166 = b766;
//			if (b1166) {
//				v[1163] = 0e0;
//				v[1162] = 0e0;
//				v[1161] = 0e0;
//			}
//			else {
//			};
//		}
//		else {
//			b1167 = b800;
//			if (b1167) {
//				v[1163] = 0e0;
//				v[1162] = 0e0;
//				v[1161] = 0e0;
//			}
//			else {
//			};
//		};
//	}
//	else {
//	};
//	v[3785] = (*ct)*v[1161];
//	v[3913] = v[1161] * v[3782] - v[768];
//	v[3784] = (*ct)*v[1162];
//	v[3911] = v[1162] * v[3782] - v[769];
//	v[3783] = (*ct)*v[1163];
//	v[3909] = v[1163] * v[3782] - v[770];
//	v[1171] = v[1175] * v[3783];
//	v[1172] = v[1177] * v[3783];
//	v[1173] = v[1179] * v[3783];
//	v[1174] = v[1181] * v[3783];
//	v[1176] = v[1175] * v[3784];
//	v[1178] = v[1177] * v[3784];
//	v[1180] = v[1179] * v[3784];
//	v[1182] = v[1181] * v[3784];
//	v[1183] = v[1175] * v[3785];
//	v[1184] = v[1177] * v[3785];
//	v[1185] = v[1179] * v[3785];
//	v[1186] = v[1181] * v[3785];
//	v[1187] = (*epst)*v[1163];
//	v[2002] = -(v[1187] * v[319]);
//	v[1188] = (*epst)*v[1162];
//	v[2005] = -(v[1188] * v[318]);
//	v[2000] = v[2002] + v[2005];
//	v[1189] = (*epst)*v[1161];
//	v[2006] = -(v[1189] * v[317]);
//	v[2007] = v[2005] + v[2006];
//	v[2003] = v[2002] + v[2006];
//	v[1190] = 0e0;
//	v[1191] = 0e0;
//	v[1192] = 0e0;
//	v[1193] = 0e0;
//	v[1194] = 0e0;
//	b1195 = b728;
//	if (b1195) {
//		v[1196] = 0e0;
//		v[1197] = 0e0;
//		v[1198] = 0e0;
//		b1199 = b747;
//		if (b1199) {
//			v[1198] = 0e0;
//			v[1197] = 0e0;
//			v[1196] = 0e0;
//		}
//		else {
//		};
//		v[1208] = v[1196] * v[714] + v[1197] * v[716] + v[1198] * v[718];
//		v[3786] = v[1208] * (*zetan);
//		b1200 = b730;
//		if (b1200) {
//			v[1192] = v[1198] * v[742];
//			v[1191] = v[1197] * v[742];
//			v[1190] = v[1196] * v[742];
//			v[1194] = (v[1201] * v[3786] * v[740] * Power(v[739], v[2115])) / v[1202];
//		}
//		else {
//			v[1192] = v[1198] * v[746];
//			v[1191] = v[1197] * v[746];
//			v[1190] = v[1196] * v[746];
//			v[1193] = (v[1206] * v[3786] * v[745] * Power(v[727], v[2120])) / v[1207];
//		};
//	}
//	else {
//	};
//	v[2676] = v[1190] * v[317];
//	v[2374] = -(v[2676] * v[301]);
//	v[2372] = -(v[2676] * v[291]);
//	v[2370] = -(v[2676] * v[281]);
//	v[2368] = v[2676] * v[275];
//	v[2366] = v[265] * v[2676];
//	v[2364] = v[255] * v[2676];
//	v[2677] = v[1191] * v[318];
//	v[3884] = v[2676] + v[2677];
//	v[2391] = -(v[2677] * v[301]);
//	v[2388] = -(v[2677] * v[291]);
//	v[2385] = -(v[2677] * v[281]);
//	v[2382] = v[2677] * v[275];
//	v[4356] = v[2368] + v[2382];
//	v[2379] = v[265] * v[2677];
//	v[4357] = v[2366] + v[2379];
//	v[2376] = v[255] * v[2677];
//	v[4358] = v[2364] + v[2376];
//	v[4354] = v[1192] * v[2356] + v[3884] * v[471];
//	v[4351] = v[1192] * v[2352] + v[3884] * v[472];
//	v[4348] = v[1192] * v[2348] + v[3884] * v[473];
//	v[4345] = v[1192] * v[2344] - v[3884] * v[570];
//	v[4342] = v[1192] * v[2340] - v[3884] * v[571];
//	v[4339] = v[1192] * v[2336] - v[3884] * v[572];
//	v[2668] = v[1192] * v[319];
//	v[3887] = v[2668] + v[2677];
//	v[4355] = v[1190] * v[3757] + v[3887] * v[465];
//	v[4352] = v[1190] * v[3755] + v[3887] * v[466];
//	v[4349] = v[1190] * v[3753] + v[3887] * v[467];
//	v[4346] = -(v[1190] * v[3751]) - v[3887] * v[564];
//	v[4343] = -(v[1190] * v[3749]) - v[3887] * v[565];
//	v[4340] = -(v[1190] * v[3747]) - v[3887] * v[566];
//	v[3882] = v[2668] + v[2676];
//	v[4353] = v[1191] * v[3756] + v[3882] * v[468];
//	v[4350] = v[1191] * v[3754] + v[3882] * v[469];
//	v[4347] = v[1191] * v[3752] + v[3882] * v[470];
//	v[4344] = v[1191] * v[3750] - v[3882] * v[567];
//	v[4341] = v[1191] * v[3748] - v[3882] * v[568];
//	v[4338] = v[1191] * v[3746] - v[3882] * v[569];
//	v[2392] = -(v[2668] * v[301]);
//	v[2389] = -(v[2668] * v[291]);
//	v[2386] = -(v[2668] * v[281]);
//	v[2383] = v[2668] * v[275];
//	v[4364] = v[2382] + v[2383];
//	v[4361] = v[2368] + v[2383];
//	v[2380] = v[265] * v[2668];
//	v[4363] = v[2379] + v[2380];
//	v[4360] = v[2366] + v[2380];
//	v[2377] = v[255] * v[2668];
//	v[4362] = v[2376] + v[2377];
//	v[4359] = v[2364] + v[2377];
//	v[2331] = -(v[2676] * v[319]) - v[2677] * v[319] - v[1192] * v[717];
//	v[3991] = (*ct)*(v[1118] * v[1161] + v[1119] * v[1162] + v[1120] * v[1163]) + v[2331];
//	v[2332] = -(v[2668] * v[318]) - v[2676] * v[318] - v[1191] * v[715];
//	v[3990] = (*ct)*(v[1114] * v[1161] + v[1115] * v[1162] + v[1116] * v[1163]) + v[2332];
//	v[2333] = -(v[2668] * v[317]) - v[2677] * v[317] - v[1190] * v[713];
//	v[3989] = (*ct)*(v[1110] * v[1161] + v[1111] * v[1162] + v[1112] * v[1163]) + v[2333];
//	v[1209] = 0e0;
//	v[1210] = 0e0;
//	v[1211] = 0e0;
//	b1212 = b728;
//	if (b1212) {
//		b1213 = b730;
//		if (b1213) {
//			v[1211] = 0e0;
//			v[1210] = 0e0;
//			v[1209] = 0e0;
//			v[1193] = v[1193] - v[1194];
//		}
//		else {
//		};
//	}
//	else {
//	};
//	v[1214] = v[1192] * v[2090];
//	v[1215] = v[1191] * v[2095];
//	v[1211] = v[1211] + v[1192] * v[2093] + v[1191] * v[2098] + v[1190] * v[2103];
//	v[1210] = v[1210] + v[1192] * v[2092] + v[1191] * v[2097] + v[1190] * v[2102];
//	v[1224] = v[1190] * v[2100];
//	v[1209] = v[1209] + v[1192] * v[2091] + v[1191] * v[2096] + v[1190] * v[2101];
//	v[1244] = -((*ct)*(v[1154] * v[1161] + v[1155] * v[1162] + v[1156] * v[1163])) + v[1190] * v[1241] + v[1191] * v[1242]
//		+ v[1192] * v[1243];
//	v[1433] = v[1244] * v[1580];
//	v[1419] = v[1244] * v[2687];
//	v[1248] = -((*ct)*(v[1150] * v[1161] + v[1151] * v[1162] + v[1152] * v[1163])) + v[1190] * v[1245] + v[1191] * v[1246]
//		+ v[1192] * v[1247];
//	v[1426] = v[1248] * v[2686];
//	v[3794] = v[1426] + v[1244] * v[1579];
//	v[3793] = v[1419] + v[1248] * v[1577];
//	v[1252] = -((*ct)*(v[1146] * v[1161] + v[1147] * v[1162] + v[1148] * v[1163])) + v[1190] * v[1249] + v[1191] * v[1250]
//		+ v[1192] * v[1251];
//	v[1431] = v[1252] * v[2685];
//	v[4029] = v[1431] + v[1433];
//	v[3795] = v[1431] + v[1248] * v[1578];
//	v[4019] = v[1433] + v[3795];
//	v[1428] = v[1252] * v[3732];
//	v[4025] = v[1428] + v[3794];
//	v[4021] = v[1426] + v[1428];
//	v[1422] = v[1252] * v[3735];
//	v[4028] = v[1422] + v[3793];
//	v[4020] = v[1419] + v[1422];
//	v[1256] = -((*ct)*(v[1130] * v[1161] + v[1131] * v[1162] + v[1132] * v[1163])) + v[1190] * v[1253] + v[1191] * v[1254]
//		+ v[1192] * v[1255];
//	v[1496] = v[1256] * v[1574];
//	v[1482] = v[1256] * v[2684];
//	v[1260] = -((*ct)*(v[1126] * v[1161] + v[1127] * v[1162] + v[1128] * v[1163])) + v[1190] * v[1257] + v[1191] * v[1258]
//		+ v[1192] * v[1259];
//	v[1489] = v[1260] * v[2683];
//	v[3812] = v[1489] + v[1256] * v[1573];
//	v[3811] = v[1482] + v[1260] * v[1571];
//	v[1264] = -((*ct)*(v[1122] * v[1161] + v[1123] * v[1162] + v[1124] * v[1163])) + v[1190] * v[1261] + v[1191] * v[1262]
//		+ v[1192] * v[1263];
//	v[1494] = v[1264] * v[2682];
//	v[5870] = -v[3989];
//	v[5871] = -v[3990];
//	v[5872] = -v[3991];
//	v[5873] = v[1494] + v[1260] * v[3724] + v[1256] * v[3727];
//	v[5874] = v[1489] + v[1256] * v[1571] + v[1264] * v[1572];
//	v[5875] = v[1482] + v[1260] * v[1573] + v[1264] * v[1574];
//	v[5876] = v[3989];
//	v[5877] = v[3990];
//	v[5878] = v[3991];
//	v[5879] = v[1431] + v[1248] * v[3732] + v[1244] * v[3735];
//	v[5880] = v[1426] + v[1244] * v[1577] + v[1252] * v[1578];
//	v[5881] = v[1419] + v[1248] * v[1579] + v[1252] * v[1580];
//	v[4017] = v[1494] + v[1496];
//	v[3813] = v[1494] + v[1260] * v[1572];
//	v[4007] = v[1496] + v[3813];
//	v[1491] = v[1264] * v[3724];
//	v[4013] = v[1491] + v[3812];
//	v[4009] = v[1489] + v[1491];
//	v[1485] = v[1264] * v[3727];
//	v[4016] = v[1485] + v[3811];
//	v[4008] = v[1482] + v[1485];
//	v[1265] = v[2333] * v[291];
//	v[1268] = v[2333] * v[301];
//	v[1269] = v[2333] * v[281];
//	v[1270] = v[2332] * v[281];
//	v[1273] = v[2332] * v[301];
//	v[1274] = v[2331] * v[301];
//	v[1275] = v[2332] * v[291];
//	v[1278] = v[2331] * v[281];
//	v[1279] = v[2331] * v[291];
//	v[1280] = -(v[2333] * v[265]);
//	v[1281] = -(v[2333] * v[275]);
//	v[1282] = -(v[2333] * v[255]);
//	v[1283] = -(v[2332] * v[255]);
//	v[1284] = -(v[2332] * v[275]);
//	v[1285] = -(v[2331] * v[275]);
//	v[1286] = -(v[2332] * v[265]);
//	v[1287] = -(v[2331] * v[255]);
//	v[1288] = -(v[2331] * v[265]);
//	v[1289] = 0e0;
//	v[1290] = 0e0;
//	v[1291] = 0e0;
//	v[1292] = 0e0;
//	v[1293] = 0e0;
//	v[1294] = 0e0;
//	v[1295] = 0e0;
//	v[1296] = 0e0;
//	v[1297] = 0e0;
//	b1298 = (*previouscontact);
//	if (b1298) {
//		v[1214] = v[1214] - v[1187] * v[710];
//		v[1215] = v[1215] - v[1188] * v[709];
//		v[1300] = v[1189] * v[1299] + v[2000] * v[317];
//		v[1302] = v[1188] * v[1301] + v[2003] * v[318];
//		v[1304] = v[1187] * v[1303] + v[2007] * v[319];
//		v[1211] = v[1211] + v[1187] * v[3787] + v[2007] * v[710];
//		v[1210] = v[1210] + v[1188] * v[3788] + v[2003] * v[709];
//		v[1224] = v[1224] - v[1189] * v[708];
//		v[1209] = v[1209] - v[1189] * v[3789] + v[2000] * v[708];
//		v[1289] = gti[0] * v[1300];
//		v[1290] = gti[1] * v[1300];
//		v[1291] = gti[2] * v[1300];
//		v[1308] = -v[1300];
//		v[1309] = -(v[1300] * v[203]);
//		v[1310] = -(v[1300] * v[202]);
//		v[1311] = -(v[1300] * v[201]);
//		v[1312] = v[1300];
//		v[1313] = v[1300] * v[147];
//		v[1314] = v[1300] * v[146];
//		v[1315] = v[1300] * v[145];
//		v[1292] = gti[0] * v[1302];
//		v[1293] = gti[1] * v[1302];
//		v[1294] = gti[2] * v[1302];
//		v[1316] = -v[1302];
//		v[1317] = -(v[1302] * v[203]);
//		v[1318] = -(v[1302] * v[202]);
//		v[1319] = -(v[1302] * v[201]);
//		v[1320] = v[1302];
//		v[1321] = v[1302] * v[147];
//		v[1322] = v[1302] * v[146];
//		v[1323] = v[1302] * v[145];
//		v[1295] = gti[0] * v[1304];
//		v[1296] = gti[1] * v[1304];
//		v[1297] = gti[2] * v[1304];
//		v[1324] = -v[1304];
//		v[1325] = -(v[1304] * v[203]);
//		v[1326] = -(v[1304] * v[202]);
//		v[1327] = -(v[1304] * v[201]);
//		v[1328] = v[1304];
//		v[1329] = v[1304] * v[147];
//		v[1330] = v[1304] * v[146];
//		v[1331] = v[1304] * v[145];
//	}
//	else {
//		v[1315] = 0e0;
//		v[1314] = 0e0;
//		v[1313] = 0e0;
//		v[1323] = 0e0;
//		v[1322] = 0e0;
//		v[1321] = 0e0;
//		v[1331] = 0e0;
//		v[1330] = 0e0;
//		v[1329] = 0e0;
//		v[1312] = 0e0;
//		v[1320] = 0e0;
//		v[1328] = 0e0;
//		v[1311] = 0e0;
//		v[1310] = 0e0;
//		v[1309] = 0e0;
//		v[1319] = 0e0;
//		v[1318] = 0e0;
//		v[1317] = 0e0;
//		v[1327] = 0e0;
//		v[1326] = 0e0;
//		v[1325] = 0e0;
//		v[1308] = 0e0;
//		v[1316] = 0e0;
//		v[1324] = 0e0;
//	};
//	b1332 = b673;
//	if (b1332) {
//		v[1358] = -(v[1297] * v[690]) / 2e0;
//		v[1357] = -(v[1293] * v[690]) / 2e0;
//		v[1356] = v[1296] * v[690];
//		v[1355] = v[1294] * v[690];
//		v[1351] = v[1295] * v[690];
//		v[1350] = v[1291] * v[690];
//		v[1347] = v[1292] * v[690];
//		v[1346] = v[1290] * v[690];
//		v[1333] = v[1355] + v[1356];
//		v[1334] = v[1350] + v[1351];
//		v[1335] = v[1346] + v[1347];
//		v[1345] = (v[1289] * v[1336]) / 2e0 + v[1290] * v[1337] + v[1291] * v[1338] + v[1292] * v[1339] + (v[1293] * v[1340])
//			/ 2e0 + v[1294] * v[1341] + v[1295] * v[1342] + v[1296] * v[1343] + (v[1297] * v[1344]) / 2e0;
//		v[1966] = (-4e0*v[1345]) / Power(v[1348], 2) + v[1357] + v[1358];
//		v[1965] = -v[1357] + v[1966] - (v[1289] * v[690]) / 2e0;
//		v[1964] = v[1357] - v[1358] + v[1965];
//		v[1349] = (-2e0*v[1346] + 2e0*v[1347] + v[1334] * v[687] + v[1333] * v[688] + 4e0*v[1964] * v[689]) / 2e0;
//		v[1354] = (2e0*v[1350] - 2e0*v[1351] + v[1335] * v[687] + 4e0*v[1965] * v[688] + v[1333] * v[689]) / 2e0;
//		v[1359] = (-2e0*v[1355] + 2e0*v[1356] + 4e0*v[1966] * v[687] + v[1335] * v[688] + v[1334] * v[689]) / 2e0;
//		v[3790] = v[1359] * v[675] + v[1354] * v[676] + v[1349] * v[677];
//		v[1952] = v[3790] * v[686];
//		v[1949] = v[3790] * v[680];
//		v[1362] = v[1949] * v[685] + v[1952] / (Power(cos(v[1360]), 2)*sqrt(v[1953]));
//		v[3859] = v[1362] / v[678];
//		v[3791] = v[1362] / v[678];
//		v[1363] = v[1349] * v[3759] + v[3791] * v[677];
//		v[1365] = v[1354] * v[3759] + v[3791] * v[676];
//		v[1366] = v[1359] * v[3759] + v[3791] * v[675];
//		v[1209] = v[1209] - v[1363] * v[327] + v[1365] * v[328];
//		v[1211] = v[1211] - v[1365] * v[326] + v[1366] * v[327];
//		v[1210] = v[1210] + v[1363] * v[326] - v[1366] * v[328];
//	}
//	else {
//	};
//	v[1211] = v[1211] + 2e0*v[1214] * v[319];
//	v[1210] = v[1210] + 2e0*v[1215] * v[318];
//	v[1209] = v[1209] + 2e0*v[1224] * v[317];
//	v[2394] = v[1209] * v[302] + v[1210] * v[303] + v[1211] * v[304];
//	v[1193] = v[1193] + v[2394] * v[315];
//	v[3792] = v[1193] / v[727];
//	v[1367] = v[1211] * v[316] + v[304] * v[3792];
//	v[1368] = v[1210] * v[316] + v[303] * v[3792];
//	v[1369] = v[1209] * v[316] + v[302] * v[3792];
//	v[1324] = v[1324] - v[1367];
//	v[1328] = v[1328] + v[1367];
//	v[1316] = v[1316] - v[1368];
//	v[1320] = v[1320] + v[1368];
//	v[1308] = v[1308] - v[1369];
//	v[1312] = v[1312] + v[1369];
//	v[1370] = v[1244] * v[2039];
//	v[1371] = v[1244] * v[2038];
//	v[1372] = v[1244] * v[2037];
//	v[1375] = v[1248] * v[2034];
//	v[1376] = v[1244] * v[288] + v[1248] * v[290];
//	v[3799] = v[1376] * v[209];
//	v[1379] = v[1248] * v[2032];
//	v[1380] = v[1370] + v[1375];
//	v[1382] = v[1248] * v[2033] + v[3799] / v[279];
//	v[1385] = v[1252] * v[2027];
//	v[1386] = v[1244] * v[282] + v[1252] * v[290];
//	v[3800] = v[1386] * v[214];
//	v[1387] = v[1248] * v[282] + v[1252] * v[288];
//	v[3801] = v[1387] * v[218];
//	v[4384] = -(v[1244] * v[3653]) - v[1248] * v[3654] - v[1252] * v[3655] + v[3731] * v[3799] + v[278] * v[3800]
//		+ v[277] * v[3801];
//	v[1389] = v[1252] * v[2029] + v[3800] / v[279];
//	v[1390] = v[1382] + v[1389];
//	v[1392] = v[1252] * v[2028] + v[3801] / v[279];
//	v[1393] = v[1371] + v[1392];
//	v[1394] = v[1256] * v[2024];
//	v[1395] = v[1256] * v[2023];
//	v[1396] = v[1256] * v[2022];
//	v[1399] = v[1260] * v[2019];
//	v[1400] = v[1256] * v[262] + v[1260] * v[264];
//	v[3817] = v[1400] * v[153];
//	v[1403] = v[1260] * v[2017];
//	v[1404] = v[1394] + v[1399];
//	v[1406] = v[1260] * v[2018] + v[3817] / v[253];
//	v[1409] = v[1264] * v[2012];
//	v[1410] = v[1256] * v[256] + v[1264] * v[264];
//	v[3818] = v[1410] * v[158];
//	v[1411] = v[1260] * v[256] + v[1264] * v[262];
//	v[3819] = v[1411] * v[162];
//	v[4388] = -(v[1256] * v[3678]) - v[1260] * v[3679] - v[1264] * v[3680] + v[3723] * v[3817] + v[252] * v[3818]
//		+ v[251] * v[3819];
//	v[1413] = v[1264] * v[2014] + v[3818] / v[253];
//	v[1414] = v[1406] + v[1413];
//	v[1416] = v[1264] * v[2013] + v[3819] / v[253];
//	v[1417] = v[1395] + v[1416];
//	v[1325] = v[1325] - v[1367] * v[200] + v[1171] * v[348] + v[1172] * v[358];
//	v[1326] = v[1326] - v[1367] * v[199] + v[1171] * v[347] + v[1172] * v[356];
//	v[1327] = v[1327] - v[1367] * v[198] + v[1171] * v[346] + v[1172] * v[354];
//	v[1418] = QBi[2][2] * v[1325] + QBi[2][1] * v[1326] + QBi[2][0] * v[1327] + v[290] * v[4028];
//	v[1421] = QBi[1][2] * v[1325] + QBi[1][1] * v[1326] + QBi[1][0] * v[1327] + v[1387] * v[3735] + v[288] * v[3793];
//	v[1423] = QBi[0][2] * v[1325] + QBi[0][1] * v[1326] + QBi[0][0] * v[1327] + v[282] * v[4020];
//	v[1317] = v[1317] - v[1368] * v[200] + v[1176] * v[348] + v[1178] * v[358];
//	v[1318] = v[1318] - v[1368] * v[199] + v[1176] * v[347] + v[1178] * v[356];
//	v[1319] = v[1319] - v[1368] * v[198] + v[1176] * v[346] + v[1178] * v[354];
//	v[1424] = QBi[2][2] * v[1317] + QBi[2][1] * v[1318] + QBi[2][0] * v[1319] + v[1386] * v[3732] + v[290] * v[3794];
//	v[1427] = QBi[1][2] * v[1317] + QBi[1][1] * v[1318] + QBi[1][0] * v[1319] + v[288] * v[4025];
//	v[1429] = QBi[0][2] * v[1317] + QBi[0][1] * v[1318] + QBi[0][0] * v[1319] + v[282] * v[4021];
//	v[1309] = v[1309] - v[1369] * v[200] + v[1183] * v[348] + v[1184] * v[358];
//	v[1310] = v[1310] - v[1369] * v[199] + v[1183] * v[347] + v[1184] * v[356];
//	v[1311] = v[1311] - v[1369] * v[198] + v[1183] * v[346] + v[1184] * v[354];
//	v[1430] = QBi[2][2] * v[1309] + QBi[2][1] * v[1310] + QBi[2][0] * v[1311] + v[1376] * v[1578] + v[290] * v[4029];
//	v[1432] = QBi[1][2] * v[1309] + QBi[1][1] * v[1310] + QBi[1][0] * v[1311] + v[288] * v[3795];
//	v[1435] = QBi[0][2] * v[1309] + QBi[0][1] * v[1310] + QBi[0][0] * v[1311] + v[282] * v[4019];
//	v[1436] = -(v[1265] * v[855]);
//	v[1437] = -(v[1265] * v[853]);
//	v[1438] = -(v[1265] * v[851]);
//	v[1439] = -(v[1268] * v[855]);
//	v[1440] = -(v[1268] * v[851]);
//	v[1441] = -(v[1268] * v[853]);
//	v[3831] = -2e0*v[1441];
//	v[1442] = -(v[1269] * v[853]);
//	v[1443] = -(v[1269] * v[851]);
//	v[1444] = -(v[1269] * v[855]);
//	v[1445] = -(v[1270] * v[855]);
//	v[1446] = -(v[1270] * v[853]);
//	v[1447] = -(v[1270] * v[851]);
//	v[3829] = -2e0*v[1447];
//	v[1448] = -(v[1273] * v[851]);
//	v[1449] = -(v[1273] * v[855]);
//	v[3832] = 2e0*v[1449];
//	v[1450] = -(v[1273] * v[853]);
//	v[1451] = -(v[1274] * v[853]);
//	v[1452] = -(v[1274] * v[855]);
//	v[1453] = -(v[1274] * v[851]);
//	v[1454] = v[1442] + v[1445] + v[1448] + v[1451] + 2e0*v[1379] * v[481] - v[1390] * v[484] - v[1380] * v[487];
//	v[1455] = -(v[1275] * v[855]);
//	v[1456] = -(v[1275] * v[851]);
//	v[1457] = -(v[1275] * v[853]);
//	v[1458] = -v[1437] - v[1440] - v[1452] - v[1455] - v[1390] * v[481] + 2e0*v[1385] * v[484] + v[1393] * v[487];
//	v[1459] = v[1432] * v[204] - v[1442] * v[474] + v[1437] * v[475] - v[1441] * v[476];
//	v[1460] = -(v[1278] * v[855]);
//	v[1461] = -(v[1278] * v[853]);
//	v[3830] = 2e0*v[1461];
//	v[1462] = -(v[1278] * v[851]);
//	v[1463] = -(v[1279] * v[853]);
//	v[1464] = -(v[1279] * v[855]);
//	v[1465] = -(v[1279] * v[851]);
//	v[1469] = -v[1443] - v[1456] - v[1460] - v[1463] - v[1380] * v[481] + v[1393] * v[484] + 2e0*v[1372] * v[487];
//	v[3987] = -v[1469] / 2e0;
//	v[1470] = v[1430] * v[204] - v[1443] * v[474] + v[1438] * v[475] - v[1440] * v[476];
//	v[1471] = v[1429] * v[204] - v[1445] * v[474] + v[1455] * v[475] - v[1449] * v[476];
//	v[1472] = v[1439] + v[1450];
//	v[1473] = v[1424] * v[204] - v[1447] * v[474] + v[1456] * v[475] - v[1448] * v[476];
//	v[1474] = v[1423] * v[204] - v[1460] * v[474] + v[1464] * v[475] - v[1452] * v[476];
//	v[1475] = v[1421] * v[204] - v[1461] * v[474] + v[1463] * v[475] - v[1451] * v[476];
//	v[1476] = v[1446] + v[1462];
//	v[1477] = v[1436] + v[1465];
//	v[6714] = 0e0;
//	v[6715] = 0e0;
//	v[6716] = 0e0;
//	v[6717] = 0e0;
//	v[6718] = 0e0;
//	v[6719] = 0e0;
//	v[6720] = 0e0;
//	v[6721] = 0e0;
//	v[6722] = 0e0;
//	v[6723] = -v[1458] / 2e0 - v[1476];
//	v[6724] = (v[1454] - 2e0*v[1477]) / 2e0;
//	v[6725] = -v[1472] + v[3987];
//	v[1478] = 1e0 / Power(v[279], 2);
//	v[4032] = -(v[1478] * v[289]);
//	v[3810] = -(v[1478] * v[288]);
//	v[3809] = -(v[1478] * v[290]);
//	v[3806] = -(v[1478] * v[282]);
//	v[3805] = -(v[1478] * v[277]);
//	v[3804] = -(v[1244] * v[1478]);
//	v[3803] = -(v[1478] * v[278]);
//	v[3802] = -(v[1478] * v[3731]);
//	v[2613] = -(v[1478] * (v[209] * v[280] + v[3729]));
//	v[2612] = -(v[1478] * (v[208] * v[280] + v[3730]));
//	v[2611] = -(v[1478] * v[3796]);
//	v[2610] = -(v[1478] * v[3733]);
//	v[2609] = -(v[1478] * v[3797]);
//	v[2608] = -(v[1478] * v[3798]);
//	v[2607] = -(v[1478] * v[280]);
//	v[3979] = v[2607] * v[282];
//	v[2606] = -(v[1478] * v[283]);
//	v[3978] = v[2606] * v[288];
//	v[3977] = -(v[1478] * v[290] * v[292]);
//	v[2604] = -(v[1478] * v[3799]);
//	v[2603] = -(v[1478] * v[3800]);
//	v[2602] = -(v[1478] * v[3801]);
//	v[2601] = v[1376] * v[3802];
//	v[2600] = v[1386] * v[3803];
//	v[2599] = v[1387] * v[3805];
//	v[2480] = v[1248] * v[3802];
//	v[2479] = v[3734] * v[3804];
//	v[2477] = v[1252] * v[2607];
//	v[4371] = v[2601] + (v[2477] + v[2479])*v[290];
//	v[3973] = v[2477] + v[2480];
//	v[4372] = v[2479] + v[3973];
//	v[2473] = v[1252] * v[3803];
//	v[2470] = v[1248] * v[2606];
//	v[4370] = v[2470] + v[2473];
//	v[2469] = v[3736] * v[3804];
//	v[3972] = v[2469] + v[2470];
//	v[4369] = v[2473] + v[3972];
//	v[4368] = v[2600] + v[290] * v[3972];
//	v[2465] = v[1252] * v[3805];
//	v[2463] = v[1248] * v[4032];
//	v[2462] = v[292] * v[3804];
//	v[4367] = v[2462] + v[2465];
//	v[3971] = v[2462] + v[2463];
//	v[4366] = v[2599] + v[288] * v[3971];
//	v[4365] = v[2465] + v[3971];
//	v[2428] = v[216] * v[3806];
//	v[2425] = v[211] * v[3806];
//	v[2422] = -(v[1478] * v[3807]);
//	v[2421] = -(v[1478] * v[3808]);
//	v[2416] = v[214] * v[3809];
//	v[2415] = v[213] * v[3810];
//	v[3423] = v[2415] + v[2425];
//	v[3415] = v[2416] + v[3423];
//	v[3401] = v[2415] + v[2416];
//	v[2413] = v[207] * v[3806];
//	v[3420] = v[2413] + v[2421] + v[2422];
//	v[3413] = -v[2422] + v[3420];
//	v[3403] = -v[2421] + v[3420];
//	v[2409] = v[219] * v[3809];
//	v[3426] = v[2409] + v[2428];
//	v[2408] = v[218] * v[3810];
//	v[3409] = v[2408] + v[2409];
//	v[3405] = v[2408] + v[3426];
//	v[2036] = v[2609] * v[282] + v[2608] * v[288] + v[3977];
//	v[2031] = v[2611] * v[282] + v[2610] * v[290] + v[3978];
//	v[2026] = v[2612] * v[288] + v[2613] * v[290] + v[3979];
//	v[1859] = v[209] * v[3802];
//	v[1857] = v[214] * v[3803];
//	v[1853] = v[218] * v[3805];
//	v[2598] = v[1372] + v[1379] + v[1385] + v[1387] * v[1853] + v[1386] * v[1857] + v[1376] * v[1859] + v[1252] * v[2026]
//		+ v[1248] * v[2031] + v[1244] * v[2036];
//	v[1479] = (2e0*v[1438] + alphaB[1] * v[1454] - alphaB[0] * v[1458] - 2e0*v[1464] - alphaB[2] * v[1469]
//		+ 4e0*v[204] * v[2598] - v[1476] * v[369] - v[1477] * v[372] - v[1472] * v[374] + v[3829] + v[3830] + v[3831] + v[3832]
//		+ v[1435] * v[477] + 2e0*v[1432] * v[482] + 2e0*v[1430] * v[488] + 2e0*v[1429] * v[492] + v[1427] * v[496]
//		+ 2e0*v[1424] * v[501] + 2e0*v[1423] * v[505] + 2e0*v[1421] * v[509] + v[1418] * v[514]) / 2e0;
//	v[1329] = v[1329] + v[1367] * v[144] + v[1173] * v[331] + v[1174] * v[341];
//	v[1330] = v[1330] + v[1367] * v[143] + v[1173] * v[330] + v[1174] * v[339];
//	v[1331] = v[1331] + v[1367] * v[142] + v[1173] * v[329] + v[1174] * v[337];
//	v[1481] = QAi[2][2] * v[1329] + QAi[2][1] * v[1330] + QAi[2][0] * v[1331] + v[264] * v[4016];
//	v[1484] = QAi[1][2] * v[1329] + QAi[1][1] * v[1330] + QAi[1][0] * v[1331] + v[1411] * v[3727] + v[262] * v[3811];
//	v[1486] = QAi[0][2] * v[1329] + QAi[0][1] * v[1330] + QAi[0][0] * v[1331] + v[256] * v[4008];
//	v[1321] = v[1321] + v[1368] * v[144] + v[1180] * v[331] + v[1182] * v[341];
//	v[1322] = v[1322] + v[1368] * v[143] + v[1180] * v[330] + v[1182] * v[339];
//	v[1323] = v[1323] + v[1368] * v[142] + v[1180] * v[329] + v[1182] * v[337];
//	v[1487] = QAi[2][2] * v[1321] + QAi[2][1] * v[1322] + QAi[2][0] * v[1323] + v[1410] * v[3724] + v[264] * v[3812];
//	v[1490] = QAi[1][2] * v[1321] + QAi[1][1] * v[1322] + QAi[1][0] * v[1323] + v[262] * v[4013];
//	v[1492] = QAi[0][2] * v[1321] + QAi[0][1] * v[1322] + QAi[0][0] * v[1323] + v[256] * v[4009];
//	v[1313] = v[1313] + v[1369] * v[144] + v[1185] * v[331] + v[1186] * v[341];
//	v[1314] = v[1314] + v[1369] * v[143] + v[1185] * v[330] + v[1186] * v[339];
//	v[1315] = v[1315] + v[1369] * v[142] + v[1185] * v[329] + v[1186] * v[337];
//	v[1493] = QAi[2][2] * v[1313] + QAi[2][1] * v[1314] + QAi[2][0] * v[1315] + v[1400] * v[1572] + v[264] * v[4017];
//	v[1495] = QAi[1][2] * v[1313] + QAi[1][1] * v[1314] + QAi[1][0] * v[1315] + v[262] * v[3813];
//	v[1498] = QAi[0][2] * v[1313] + QAi[0][1] * v[1314] + QAi[0][0] * v[1315] + v[256] * v[4007];
//	v[1499] = v[1280] * v[861];
//	v[1500] = v[1280] * v[859];
//	v[1501] = v[1280] * v[857];
//	v[1502] = v[1281] * v[861];
//	v[1503] = v[1281] * v[857];
//	v[1504] = v[1281] * v[859];
//	v[3835] = -2e0*v[1504];
//	v[1505] = v[1282] * v[859];
//	v[1506] = v[1282] * v[857];
//	v[1507] = v[1282] * v[861];
//	v[1508] = v[1283] * v[861];
//	v[1509] = v[1283] * v[859];
//	v[1510] = v[1283] * v[857];
//	v[3833] = -2e0*v[1510];
//	v[1511] = v[1284] * v[857];
//	v[1512] = v[1284] * v[861];
//	v[3836] = 2e0*v[1512];
//	v[1513] = v[1284] * v[859];
//	v[1514] = v[1285] * v[859];
//	v[1515] = v[1285] * v[861];
//	v[1516] = v[1285] * v[857];
//	v[1517] = v[1505] + v[1508] + v[1511] + v[1514] + 2e0*v[1403] * v[382] - v[1414] * v[385] - v[1404] * v[388];
//	v[1518] = v[1286] * v[861];
//	v[1519] = v[1286] * v[857];
//	v[1520] = v[1286] * v[859];
//	v[1521] = -v[1500] - v[1503] - v[1515] - v[1518] - v[1414] * v[382] + 2e0*v[1409] * v[385] + v[1417] * v[388];
//	v[1522] = v[148] * v[1495] - v[1505] * v[375] + v[1500] * v[376] - v[1504] * v[377];
//	v[1523] = v[1287] * v[861];
//	v[1524] = v[1287] * v[859];
//	v[3834] = 2e0*v[1524];
//	v[1525] = v[1287] * v[857];
//	v[1526] = v[1288] * v[859];
//	v[1527] = v[1288] * v[861];
//	v[1528] = v[1288] * v[857];
//	v[1532] = -v[1506] - v[1519] - v[1523] - v[1526] - v[1404] * v[382] + v[1417] * v[385] + 2e0*v[1396] * v[388];
//	v[3984] = -v[1532] / 2e0;
//	v[1533] = v[148] * v[1493] - v[1506] * v[375] + v[1501] * v[376] - v[1503] * v[377];
//	v[1534] = v[148] * v[1492] - v[1508] * v[375] + v[1518] * v[376] - v[1512] * v[377];
//	v[1535] = v[1502] + v[1513];
//	v[1536] = v[148] * v[1487] - v[1510] * v[375] + v[1519] * v[376] - v[1511] * v[377];
//	v[1537] = v[148] * v[1486] - v[1523] * v[375] + v[1527] * v[376] - v[1515] * v[377];
//	v[1538] = v[148] * v[1484] - v[1524] * v[375] + v[1526] * v[376] - v[1514] * v[377];
//	v[1539] = v[1509] + v[1525];
//	v[1540] = v[1499] + v[1528];
//	v[6726] = 0e0;
//	v[6727] = 0e0;
//	v[6728] = 0e0;
//	v[6729] = -v[1521] / 2e0 - v[1539];
//	v[6730] = (v[1517] - 2e0*v[1540]) / 2e0;
//	v[6731] = -v[1535] + v[3984];
//	v[6732] = 0e0;
//	v[6733] = 0e0;
//	v[6734] = 0e0;
//	v[6735] = 0e0;
//	v[6736] = 0e0;
//	v[6737] = 0e0;
//	v[1541] = 1e0 / Power(v[253], 2);
//	v[4034] = -(v[1541] * v[263]);
//	v[3828] = -(v[1541] * v[262]);
//	v[3827] = -(v[1541] * v[264]);
//	v[3824] = -(v[1541] * v[256]);
//	v[3823] = -(v[1541] * v[251]);
//	v[3822] = -(v[1256] * v[1541]);
//	v[3821] = -(v[1541] * v[252]);
//	v[3820] = -(v[1541] * v[3723]);
//	v[2642] = -(v[1541] * (v[153] * v[254] + v[3721]));
//	v[2641] = -(v[1541] * (v[152] * v[254] + v[3722]));
//	v[2640] = -(v[1541] * v[3814]);
//	v[2639] = -(v[1541] * v[3725]);
//	v[2638] = -(v[1541] * v[3815]);
//	v[2637] = -(v[1541] * v[3816]);
//	v[2636] = -(v[1541] * v[254]);
//	v[3982] = v[256] * v[2636];
//	v[2635] = -(v[1541] * v[257]);
//	v[3981] = v[262] * v[2635];
//	v[3980] = -(v[1541] * v[264] * v[266]);
//	v[2633] = -(v[1541] * v[3817]);
//	v[2632] = -(v[1541] * v[3818]);
//	v[2631] = -(v[1541] * v[3819]);
//	v[2630] = v[1400] * v[3820];
//	v[2629] = v[1410] * v[3821];
//	v[2628] = v[1411] * v[3823];
//	v[2544] = v[1260] * v[3820];
//	v[2543] = v[3726] * v[3822];
//	v[2541] = v[1264] * v[2636];
//	v[4379] = v[2630] + (v[2541] + v[2543])*v[264];
//	v[3976] = v[2541] + v[2544];
//	v[4380] = v[2543] + v[3976];
//	v[2537] = v[1264] * v[3821];
//	v[2534] = v[1260] * v[2635];
//	v[4378] = v[2534] + v[2537];
//	v[2533] = v[3728] * v[3822];
//	v[3975] = v[2533] + v[2534];
//	v[4377] = v[2537] + v[3975];
//	v[4376] = v[2629] + v[264] * v[3975];
//	v[2529] = v[1264] * v[3823];
//	v[2527] = v[1260] * v[4034];
//	v[2526] = v[266] * v[3822];
//	v[4375] = v[2526] + v[2529];
//	v[3974] = v[2526] + v[2527];
//	v[4374] = v[2628] + v[262] * v[3974];
//	v[4373] = v[2529] + v[3974];
//	v[2458] = v[160] * v[3824];
//	v[2455] = v[155] * v[3824];
//	v[2452] = -(v[1541] * v[3825]);
//	v[2451] = -(v[1541] * v[3826]);
//	v[2446] = v[158] * v[3827];
//	v[2445] = v[157] * v[3828];
//	v[3453] = v[2445] + v[2455];
//	v[3445] = v[2446] + v[3453];
//	v[3431] = v[2445] + v[2446];
//	v[2443] = v[151] * v[3824];
//	v[3450] = v[2443] + v[2451] + v[2452];
//	v[3443] = -v[2452] + v[3450];
//	v[3433] = -v[2451] + v[3450];
//	v[2439] = v[163] * v[3827];
//	v[3456] = v[2439] + v[2458];
//	v[2438] = v[162] * v[3828];
//	v[3439] = v[2438] + v[2439];
//	v[3435] = v[2438] + v[3456];
//	v[2021] = v[262] * v[2637] + v[256] * v[2638] + v[3980];
//	v[2016] = v[2639] * v[264] + v[256] * v[2640] + v[3981];
//	v[2011] = v[262] * v[2641] + v[264] * v[2642] + v[3982];
//	v[1831] = v[153] * v[3820];
//	v[1829] = v[158] * v[3821];
//	v[1825] = v[162] * v[3823];
//	v[2627] = v[1396] + v[1403] + v[1409] + v[1411] * v[1825] + v[1410] * v[1829] + v[1400] * v[1831] + v[1264] * v[2011]
//		+ v[1260] * v[2016] + v[1256] * v[2021];
//	v[1542] = (2e0*v[1501] + alphaA[1] * v[1517] - alphaA[0] * v[1521] - 2e0*v[1527] - alphaA[2] * v[1532]
//		+ 4e0*v[148] * v[2627] - v[1539] * v[363] - v[1540] * v[366] - v[1535] * v[368] + v[1498] * v[378] + 2e0*v[1495] * v[383]
//		+ v[3833] + v[3834] + v[3835] + v[3836] + 2e0*v[1493] * v[389] + 2e0*v[1492] * v[393] + v[1490] * v[397]
//		+ 2e0*v[1487] * v[402] + 2e0*v[1486] * v[406] + 2e0*v[1484] * v[410] + v[1481] * v[415]) / 2e0;
//	v[1544] = (-2e0*v[1370] + 2e0*v[1375] - v[1444] * v[477] - 2e0*v[1442] * v[482] - 2e0*v[1443] * v[488]
//		- 2e0*v[1445] * v[492] - v[1446] * v[496] + v[3829] * v[501] - 2e0*v[1460] * v[505] - v[3830] * v[509] - v[1462] * v[514])
//		/ 2e0;
//	v[4001] = 8e0*v[1544];
//	v[1546] = -v[1371] + v[1392] + (v[1436] * v[477]) / 2e0 + v[1437] * v[482] + v[1438] * v[488] + v[1455] * v[492] +
//		(v[1457] * v[496]) / 2e0 + v[1456] * v[501] + v[1464] * v[505] + v[1463] * v[509] + (v[1465] * v[514]) / 2e0;
//	v[4000] = 8e0*v[1546];
//	v[1547] = (v[1427] * v[204] - v[1446] * v[474] + v[1457] * v[475] - v[1450] * v[476]) / 2e0;
//	v[1548] = (-2e0*v[1382] + 2e0*v[1389] - v[1439] * v[477] + v[3831] * v[482] - 2e0*v[1440] * v[488] - v[3832] * v[492]
//		- v[1450] * v[496] - 2e0*v[1448] * v[501] - 2e0*v[1452] * v[505] - 2e0*v[1451] * v[509] - v[1453] * v[514]) / 2e0;
//	v[4382] = v[1544] * v[369] - v[1546] * v[372] + v[1548] * v[374];
//	v[3999] = 8e0*v[1548];
//	v[6738] = 0e0;
//	v[6739] = 0e0;
//	v[6740] = 0e0;
//	v[6741] = 0e0;
//	v[6742] = 0e0;
//	v[6743] = 0e0;
//	v[6744] = 0e0;
//	v[6745] = 0e0;
//	v[6746] = 0e0;
//	v[6747] = v[4001];
//	v[6748] = -v[4000];
//	v[6749] = v[3999];
//	v[1565] = -v[1547] + v[1548] * v[1655] + v[1546] * v[1658] + v[1544] * v[1660] - v[1479] * v[3988];
//	v[2660] = v[1565] + (-(v[1435] * v[204]) + v[1444] * v[474] - v[1436] * v[475] + v[1439] * v[476]) / 2e0;
//	v[1549] = (v[1418] * v[204] - v[1462] * v[474] + v[1465] * v[475] - v[1453] * v[476]) / 2e0;
//	v[2658] = v[1547] - v[1549] + v[2660];
//	v[2654] = -v[1549] + v[1565];
//	v[2656] = (v[1470] + v[1474]) / 2e0;
//	v[1551] = v[1473] + v[1475];
//	v[2659] = v[1551] / 2e0;
//	v[1552] = v[1459] + v[1471];
//	v[2655] = v[1552] / 2e0;
//	v[1553] = (-2e0*v[1394] + 2e0*v[1399] - v[1507] * v[378] - 2e0*v[1505] * v[383] - 2e0*v[1506] * v[389]
//		- 2e0*v[1508] * v[393] - v[1509] * v[397] + v[3833] * v[402] - 2e0*v[1523] * v[406] - v[3834] * v[410] - v[1525] * v[415])
//		/ 2e0;
//	v[3995] = 8e0*v[1553];
//	v[1555] = -v[1395] + v[1416] + (v[1499] * v[378]) / 2e0 + v[1500] * v[383] + v[1501] * v[389] + v[1518] * v[393] +
//		(v[1520] * v[397]) / 2e0 + v[1519] * v[402] + v[1527] * v[406] + v[1526] * v[410] + (v[1528] * v[415]) / 2e0;
//	v[3994] = 8e0*v[1555];
//	v[1556] = (v[148] * v[1490] - v[1509] * v[375] + v[1520] * v[376] - v[1513] * v[377]) / 2e0;
//	v[1557] = (-2e0*v[1406] + 2e0*v[1413] - v[1502] * v[378] + v[383] * v[3835] - 2e0*v[1503] * v[389] - v[3836] * v[393]
//		- v[1513] * v[397] - 2e0*v[1511] * v[402] - 2e0*v[1515] * v[406] - 2e0*v[1514] * v[410] - v[1516] * v[415]) / 2e0;
//	v[4386] = v[1553] * v[363] - v[1555] * v[366] + v[1557] * v[368];
//	v[3993] = 8e0*v[1557];
//	v[6786] = 0e0;
//	v[6787] = 0e0;
//	v[6788] = 0e0;
//	v[6789] = v[3995];
//	v[6790] = -v[3994];
//	v[6791] = v[3993];
//	v[6792] = 0e0;
//	v[6793] = 0e0;
//	v[6794] = 0e0;
//	v[6795] = 0e0;
//	v[6796] = 0e0;
//	v[6797] = 0e0;
//	v[1563] = -v[1556] + v[1557] * v[1629] + v[1555] * v[1632] + v[1553] * v[1634] - v[1542] * v[3985];
//	v[2653] = v[1563] + (-(v[148] * v[1498]) + v[1507] * v[375] - v[1499] * v[376] + v[1502] * v[377]) / 2e0;
//	v[1558] = (v[148] * v[1481] - v[1525] * v[375] + v[1528] * v[376] - v[1516] * v[377]) / 2e0;
//	v[2651] = v[1556] - v[1558] + v[2653];
//	v[2647] = -v[1558] + v[1563];
//	v[2649] = (v[1533] + v[1537]) / 2e0;
//	v[1560] = v[1536] + v[1538];
//	v[2652] = v[1560] / 2e0;
//	v[1561] = v[1522] + v[1534];
//	v[5858] = v[1113] + v[1312];
//	v[5859] = v[1117] + v[1320];
//	v[5860] = v[1121] + v[1328];
//	v[5861] = -v[1536] + v[1538] + 2e0*alphaA[0] * v[2647] + alphaA[2] * v[2649] + v[1561] * v[364] + v[1521] * v[3983]
//		+ 2e0*(v[1539] * v[3983] + v[1553] * v[3985]) - v[1122] * v[768] - v[1123] * v[769] - v[1124] * v[770];
//	v[5862] = (v[148] * v[1517]) / 2e0 + v[1533] - v[1537] + (alphaA[2] * v[1560]) / 2e0 + (alphaA[0] * v[1561]) / 2e0
//		+ 2e0*alphaA[1] * v[2651] + 2e0*(v[1540] * v[3983] - v[1555] * v[3985]) - v[1126] * v[768] - v[1127] * v[769]
//		- v[1128] * v[770];
//	v[5863] = -v[1522] + v[1534] + alphaA[0] * v[2649] + 2e0*alphaA[2] * v[2653] + v[1560] * v[364] + v[148] * v[3984] + 2e0*
//		(v[1535] * v[3983] + v[1557] * v[3985]) - v[1130] * v[768] - v[1131] * v[769] - v[1132] * v[770];
//	v[5864] = -v[1113] + v[1308];
//	v[5865] = -v[1117] + v[1316];
//	v[5866] = -v[1121] + v[1324];
//	v[5867] = -v[1473] + v[1475] + 2e0*alphaB[0] * v[2654] + alphaB[2] * v[2656] + v[1552] * v[370] + v[1458] * v[3986]
//		+ 2e0*(v[1476] * v[3986] + v[1544] * v[3988]) - v[1146] * v[768] - v[1147] * v[769] - v[1148] * v[770];
//	v[5868] = v[1470] - v[1474] + (alphaB[2] * v[1551]) / 2e0 + (alphaB[0] * v[1552]) / 2e0 + (v[1454] * v[204]) / 2e0
//		+ 2e0*alphaB[1] * v[2658] + 2e0*(v[1477] * v[3986] - v[1546] * v[3988]) - v[1150] * v[768] - v[1151] * v[769]
//		- v[1152] * v[770];
//	v[5869] = -v[1459] + v[1471] + alphaB[0] * v[2656] + 2e0*alphaB[2] * v[2660] + v[1551] * v[370] + v[204] * v[3987] + 2e0*
//		(v[1472] * v[3986] + v[1548] * v[3988]) - v[1154] * v[768] - v[1155] * v[769] - v[1156] * v[770];
//	v[2648] = v[1561] / 2e0;
//	for (i1159 = 1; i1159 <= 12; i1159++) {
//		i3901 = (i1159 == 9 ? 1 : 0);
//		i3900 = (i1159 == 8 ? 1 : 0);
//		i3899 = (i1159 == 7 ? 1 : 0);
//		i3898 = (i1159 == 3 ? 1 : 0);
//		i3897 = (i1159 == 2 ? 1 : 0);
//		i3896 = (i1159 == 1 ? 1 : 0);
//		ii3862 = i3896 - i3899;
//		ii3861 = i3897 - i3900;
//		ii3860 = i3898 - i3901;
//		i3845 = (i1159 == 11 ? 1 : 0);
//		v[3961] = (*a4)*i3845;
//		v[3848] = i3845 * v[204];
//		i3844 = (i1159 == 10 ? 1 : 0);
//		v[3962] = (*a4)*i3844;
//		v[3847] = -(i3844*v[204]);
//		i3843 = (i1159 == 12 ? 1 : 0);
//		v[3963] = (*a4)*i3843;
//		v[3846] = -(i3843*v[204]);
//		i3839 = (i1159 == 5 ? 1 : 0);
//		v[3966] = (*a4)*i3839;
//		v[3842] = i3839 * v[148];
//		i3838 = (i1159 == 4 ? 1 : 0);
//		v[3967] = (*a4)*i3838;
//		v[3841] = -(i3838*v[148]);
//		i3837 = (i1159 == 6 ? 1 : 0);
//		v[3968] = (*a4)*i3837;
//		v[3840] = -(i3837*v[148]);
//		v[1586] = v[5433 + i1159];
//		v[1588] = v[5457 + i1159];
//		v[1590] = v[5445 + i1159];
//		v[1591] = v[5885 + i1159];
//		v[1593] = v[5421 + i1159];
//		v[1640] = 2e0*v[1593] * v[962];
//		v[1706] = -4e0*v[148] * v[1640];
//		v[3851] = v[1706] * v[253];
//		v[1681] = -2e0*v[1640];
//		v[1594] = v[5933 + i1159];
//		v[1596] = v[5481 + i1159];
//		v[1598] = v[5505 + i1159];
//		v[1600] = v[5493 + i1159];
//		v[1601] = v[5957 + i1159];
//		v[1603] = v[5469 + i1159];
//		v[1666] = 2e0*v[1603] * v[968];
//		v[1780] = -4e0*v[1666] * v[204];
//		v[3852] = v[1780] * v[279];
//		v[1755] = -2e0*v[1666];
//		v[1604] = v[6005 + i1159];
//		v[1623] = i3837 + v[1586];
//		v[1624] = -i3837 + v[1586];
//		v[1625] = i3838 + v[1588];
//		v[1626] = -i3838 + v[1588];
//		v[1627] = -i3839 + v[1590];
//		v[1628] = i3839 + v[1590];
//		v[1630] = v[1593] * v[1629] + 8e0*i3837*v[962];
//		v[1631] = 2e0*alphaA[1] * i3839 - v[1593];
//		v[1633] = v[1593] * v[1632] - 8e0*i3839*v[962];
//		v[1635] = v[1593] * v[1634] + 8e0*i3838*v[962];
//		v[3849] = 2e0*(-(v[1593] * v[376]) / 2e0 - v[3842]);
//		v[1637] = (v[1593] * v[375]) / 2e0 + v[3841];
//		v[1638] = (v[1593] * v[377]) / 2e0 + v[3840];
//		v[1639] = alphaA[2] * v[1640] + v[3840] / 2e0;
//		v[1641] = alphaA[0] * v[1640] + v[3841] / 2e0;
//		v[1642] = -(alphaA[1] * v[1640]) + v[3842] / 2e0;
//		v[1643] = (v[148] * v[1591]) / 2e0 - v[1640] * v[415];
//		v[1742] = v[1643] * v[264];
//		v[1644] = (-(v[1591] * v[377]) - v[1630] * v[415]) / 2e0;
//		v[1645] = (v[148] * v[1631]) / 2e0 - v[1640] * v[397];
//		v[1735] = v[1645] * v[262];
//		v[1646] = (v[1631] * v[376] + v[1633] * v[397]) / 2e0;
//		v[1647] = (v[148] * v[1594]) / 2e0 - v[1640] * v[378];
//		v[1727] = v[1647] * v[256];
//		v[1648] = (-(v[1594] * v[375]) - v[1635] * v[378]) / 2e0;
//		v[1649] = i3843 + v[1596];
//		v[1650] = -i3843 + v[1596];
//		v[1651] = i3844 + v[1598];
//		v[1652] = -i3844 + v[1598];
//		v[1653] = -i3845 + v[1600];
//		v[1654] = i3845 + v[1600];
//		v[1656] = v[1603] * v[1655] + 8e0*i3843*v[968];
//		v[1657] = 2e0*alphaB[1] * i3845 - v[1603];
//		v[1659] = v[1603] * v[1658] - 8e0*i3845*v[968];
//		v[1661] = v[1603] * v[1660] + 8e0*i3844*v[968];
//		v[3850] = 2e0*(-v[3848] - (v[1603] * v[475]) / 2e0);
//		v[1663] = v[3847] + (v[1603] * v[474]) / 2e0;
//		v[1664] = v[3846] + (v[1603] * v[476]) / 2e0;
//		v[1665] = alphaB[2] * v[1666] + v[3846] / 2e0;
//		v[1667] = alphaB[0] * v[1666] + v[3847] / 2e0;
//		v[1668] = -(alphaB[1] * v[1666]) + v[3848] / 2e0;
//		v[1669] = (v[1601] * v[204]) / 2e0 - v[1666] * v[514];
//		v[1816] = v[1669] * v[290];
//		v[1670] = (-(v[1601] * v[476]) - v[1656] * v[514]) / 2e0;
//		v[1671] = (v[1657] * v[204]) / 2e0 - v[1666] * v[496];
//		v[1809] = v[1671] * v[288];
//		v[1672] = (v[1657] * v[475] + v[1659] * v[496]) / 2e0;
//		v[1673] = (v[1604] * v[204]) / 2e0 - v[1666] * v[477];
//		v[1801] = v[1673] * v[282];
//		v[1674] = (-(v[1604] * v[474]) - v[1661] * v[477]) / 2e0;
//		v[1675] = (v[1591] * v[376] + v[3849] + v[1633] * v[415]) / 2e0;
//		v[1676] = (v[1594] * v[376] + v[1633] * v[378] + v[3849]) / 2e0;
//		v[1677] = v[1637] - (v[1591] * v[375]) / 2e0 - (v[1635] * v[415]) / 2e0;
//		v[1678] = v[1637] - (v[1631] * v[375]) / 2e0 - (v[1635] * v[397]) / 2e0;
//		v[1679] = v[1681] - v[1625] * v[375] - v[1635] * v[410];
//		v[1680] = v[148] * v[1625] - 2e0*v[1640] * v[410];
//		v[1682] = -v[1681] + v[1627] * v[376] + v[1633] * v[406];
//		v[1683] = v[148] * v[1627] - 2e0*v[1640] * v[406];
//		v[1744] = v[1683] * v[256];
//		v[1684] = -v[1681] - v[1626] * v[375] - v[1635] * v[402];
//		v[1685] = v[148] * v[1626] - 2e0*v[1640] * v[402];
//		v[1736] = v[1685] * v[264];
//		v[1686] = v[1638] - (v[1631] * v[377]) / 2e0 - (v[1630] * v[397]) / 2e0;
//		v[1687] = v[1638] - (v[1594] * v[377]) / 2e0 - (v[1630] * v[378]) / 2e0;
//		v[1688] = v[1681] - v[1623] * v[377] - v[1630] * v[393];
//		v[1689] = v[148] * v[1623] - 2e0*v[1640] * v[393];
//		v[1690] = v[1681] + v[1628] * v[376] + v[1633] * v[389];
//		v[1691] = v[148] * v[1628] - 2e0*v[1640] * v[389];
//		v[1728] = v[1691] * v[264];
//		v[1692] = -v[1639] + v[1625] * v[376] + v[1633] * v[410];
//		v[1693] = -v[1639] - v[1627] * v[375] - v[1635] * v[406];
//		v[1694] = -v[1639] + v[1626] * v[376] + v[1633] * v[402];
//		v[1695] = -v[1639] - v[1628] * v[375] - v[1635] * v[389];
//		v[1696] = v[1706] + 2e0*v[1639] * v[388];
//		v[1697] = v[1675] * v[857] + v[1692] * v[859] + v[1682] * v[861];
//		v[3869] = -(v[1697] * v[265]);
//		v[1698] = v[1677] * v[857] + v[1679] * v[859] + v[1693] * v[861];
//		v[3870] = -(v[1698] * v[255]);
//		v[1699] = -v[1681] - v[1624] * v[377] - v[1630] * v[383];
//		v[1700] = v[148] * v[1624] - 2e0*v[1640] * v[383];
//		v[1701] = -v[1641] + v[1623] * v[376] + v[1633] * v[393];
//		v[1702] = -v[1641] - v[1627] * v[377] - v[1630] * v[406];
//		v[1703] = -v[1641] - v[1628] * v[377] - v[1630] * v[389];
//		v[1704] = -v[1641] + v[1624] * v[376] + v[1633] * v[383];
//		v[1705] = v[1639] * v[385] + v[1641] * v[388];
//		v[1707] = v[1706] + 2e0*v[1641] * v[385];
//		v[1843] = v[1264] * v[1707];
//		v[1708] = v[1694] * v[857] + v[1646] * v[859] + v[1701] * v[861];
//		v[3875] = -(v[1708] * v[265]);
//		v[1709] = v[1642] - v[1625] * v[377] - v[1630] * v[410];
//		v[1710] = v[1642] - v[1626] * v[377] - v[1630] * v[402];
//		v[1711] = v[1642] - v[1623] * v[375] - v[1635] * v[393];
//		v[1712] = v[1642] - v[1624] * v[375] - v[1635] * v[383];
//		v[1713] = -(v[1641] * v[382]) - v[1642] * v[385];
//		v[1714] = -(v[1639] * v[382]) - v[1642] * v[388];
//		v[1715] = v[1706] + 2e0*v[1642] * v[382];
//		v[1847] = v[1260] * v[1715];
//		v[1716] = v[1644] * v[857] + v[1709] * v[859] + v[1702] * v[861];
//		v[3871] = -(v[1716] * v[275]);
//		v[1717] = v[1710] * v[857] + v[1686] * v[859] + v[1688] * v[861];
//		v[3876] = -(v[1717] * v[275]);
//		v[1718] = v[1684] * v[857] + v[1678] * v[859] + v[1711] * v[861];
//		v[3877] = -(v[1718] * v[255]);
//		v[1719] = v[1695] * v[857] + v[1712] * v[859] + v[1648] * v[861];
//		v[3890] = -(v[1719] * v[255]);
//		v[1720] = v[1703] * v[857] + v[1699] * v[859] + v[1687] * v[861];
//		v[3891] = -(v[1720] * v[275]);
//		v[1721] = v[1285] * v[1644] + v[1288] * v[1675] + v[1287] * v[1677] + v[1283] * v[1684] + v[1280] * v[1690]
//			+ v[1286] * v[1694] + v[1282] * v[1695] + v[1281] * v[1703] + v[1284] * v[1710];
//		v[1722] = v[1286] * v[1646] + v[1283] * v[1678] + v[1287] * v[1679] + v[1284] * v[1686] + v[1288] * v[1692]
//			+ v[1281] * v[1699] + v[1280] * v[1704] + v[1285] * v[1709] + v[1282] * v[1712];
//		v[1723] = v[1690] * v[857] + v[1704] * v[859] + v[1676] * v[861];
//		v[3892] = -(v[1723] * v[265]);
//		v[1724] = v[1282] * v[1648] + v[1280] * v[1676] + v[1288] * v[1682] + v[1281] * v[1687] + v[1284] * v[1688]
//			+ v[1287] * v[1693] + v[1286] * v[1701] + v[1285] * v[1702] + v[1283] * v[1711];
//		v[1725] = v[1727] + v[1700] * v[262];
//		v[1726] = v[1725] + v[1728] + v[3967];
//		v[1729] = v[1727] + v[1728];
//		v[1730] = QAi[0][0] * v[1647] + QAi[2][0] * v[1691] + QAi[1][0] * v[1700];
//		v[1731] = QAi[0][1] * v[1647] + QAi[2][1] * v[1691] + QAi[1][1] * v[1700];
//		v[1732] = QAi[0][2] * v[1647] + QAi[2][2] * v[1691] + QAi[1][2] * v[1700];
//		v[1733] = v[1735] + v[1689] * v[256];
//		v[1734] = v[1733] + v[1736] + v[3966];
//		v[1737] = v[1735] + v[1736];
//		v[1738] = QAi[1][0] * v[1645] + QAi[2][0] * v[1685] + QAi[0][0] * v[1689];
//		v[1739] = QAi[1][1] * v[1645] + QAi[2][1] * v[1685] + QAi[0][1] * v[1689];
//		v[1740] = QAi[1][2] * v[1645] + QAi[2][2] * v[1685] + QAi[0][2] * v[1689];
//		v[1741] = v[1742] + v[1744];
//		v[1743] = v[1742] + v[1680] * v[262];
//		v[1745] = v[1743] + v[1744] + v[3968];
//		v[1746] = QAi[2][0] * v[1643] + QAi[1][0] * v[1680] + QAi[0][0] * v[1683];
//		v[1747] = QAi[2][1] * v[1643] + QAi[1][1] * v[1680] + QAi[0][1] * v[1683];
//		v[1748] = QAi[2][2] * v[1643] + QAi[1][2] * v[1680] + QAi[0][2] * v[1683];
//		v[1749] = (v[3850] + v[1601] * v[475] + v[1659] * v[514]) / 2e0;
//		v[1750] = (v[3850] + v[1604] * v[475] + v[1659] * v[477]) / 2e0;
//		v[1751] = v[1663] - (v[1601] * v[474]) / 2e0 - (v[1661] * v[514]) / 2e0;
//		v[1752] = v[1663] - (v[1657] * v[474]) / 2e0 - (v[1661] * v[496]) / 2e0;
//		v[1753] = v[1755] - v[1651] * v[474] - v[1661] * v[509];
//		v[1754] = v[1651] * v[204] - 2e0*v[1666] * v[509];
//		v[1756] = -v[1755] + v[1653] * v[475] + v[1659] * v[505];
//		v[1757] = v[1653] * v[204] - 2e0*v[1666] * v[505];
//		v[1818] = v[1757] * v[282];
//		v[1758] = -v[1755] - v[1652] * v[474] - v[1661] * v[501];
//		v[1759] = v[1652] * v[204] - 2e0*v[1666] * v[501];
//		v[1810] = v[1759] * v[290];
//		v[1760] = v[1664] - (v[1657] * v[476]) / 2e0 - (v[1656] * v[496]) / 2e0;
//		v[1761] = v[1664] - (v[1604] * v[476]) / 2e0 - (v[1656] * v[477]) / 2e0;
//		v[1762] = v[1755] - v[1649] * v[476] - v[1656] * v[492];
//		v[1763] = v[1649] * v[204] - 2e0*v[1666] * v[492];
//		v[1764] = v[1755] + v[1654] * v[475] + v[1659] * v[488];
//		v[1765] = v[1654] * v[204] - 2e0*v[1666] * v[488];
//		v[1802] = v[1765] * v[290];
//		v[1766] = -v[1665] + v[1651] * v[475] + v[1659] * v[509];
//		v[1767] = -v[1665] - v[1653] * v[474] - v[1661] * v[505];
//		v[1768] = -v[1665] + v[1652] * v[475] + v[1659] * v[501];
//		v[1769] = -v[1665] - v[1654] * v[474] - v[1661] * v[488];
//		v[1770] = v[1780] + 2e0*v[1665] * v[487];
//		v[1771] = -(v[1749] * v[851]) - v[1766] * v[853] - v[1756] * v[855];
//		v[3872] = v[1771] * v[291];
//		v[1772] = -(v[1751] * v[851]) - v[1753] * v[853] - v[1767] * v[855];
//		v[3873] = v[1772] * v[281];
//		v[1773] = -v[1755] - v[1650] * v[476] - v[1656] * v[482];
//		v[1774] = v[1650] * v[204] - 2e0*v[1666] * v[482];
//		v[1775] = -v[1667] + v[1649] * v[475] + v[1659] * v[492];
//		v[1776] = -v[1667] - v[1653] * v[476] - v[1656] * v[505];
//		v[1777] = -v[1667] - v[1654] * v[476] - v[1656] * v[488];
//		v[1778] = -v[1667] + v[1650] * v[475] + v[1659] * v[482];
//		v[1779] = v[1665] * v[484] + v[1667] * v[487];
//		v[1781] = v[1780] + 2e0*v[1667] * v[484];
//		v[1871] = v[1252] * v[1781];
//		v[1782] = -(v[1768] * v[851]) - v[1672] * v[853] - v[1775] * v[855];
//		v[3878] = v[1782] * v[291];
//		v[1783] = v[1668] - v[1651] * v[476] - v[1656] * v[509];
//		v[1784] = v[1668] - v[1652] * v[476] - v[1656] * v[501];
//		v[1785] = v[1668] - v[1649] * v[474] - v[1661] * v[492];
//		v[1786] = v[1668] - v[1650] * v[474] - v[1661] * v[482];
//		v[1787] = -(v[1667] * v[481]) - v[1668] * v[484];
//		v[1788] = -(v[1665] * v[481]) - v[1668] * v[487];
//		v[1789] = v[1780] + 2e0*v[1668] * v[481];
//		v[1875] = v[1248] * v[1789];
//		v[1790] = -(v[1670] * v[851]) - v[1783] * v[853] - v[1776] * v[855];
//		v[3874] = v[1790] * v[301];
//		v[3919] = -v[3869] - v[3870] - v[3871] - v[3872] - v[3873] - v[3874] + (*a4)*v[5757 + i1159];
//		v[1791] = -(v[1784] * v[851]) - v[1760] * v[853] - v[1762] * v[855];
//		v[3879] = v[1791] * v[301];
//		v[1792] = -(v[1758] * v[851]) - v[1752] * v[853] - v[1785] * v[855];
//		v[3880] = v[1792] * v[281];
//		v[3920] = -v[3875] - v[3876] - v[3877] - v[3878] - v[3879] - v[3880] + (*a4)*v[5745 + i1159];
//		v[1793] = -(v[1769] * v[851]) - v[1786] * v[853] - v[1674] * v[855];
//		v[3893] = v[1793] * v[281];
//		v[1794] = -(v[1777] * v[851]) - v[1773] * v[853] - v[1761] * v[855];
//		v[3894] = v[1794] * v[301];
//		v[1795] = -(v[1274] * v[1670]) - v[1279] * v[1749] - v[1278] * v[1751] - v[1270] * v[1758] - v[1265] * v[1764]
//			- v[1275] * v[1768] - v[1269] * v[1769] - v[1268] * v[1777] - v[1273] * v[1784];
//		v[1796] = -(v[1275] * v[1672]) - v[1270] * v[1752] - v[1278] * v[1753] - v[1273] * v[1760] - v[1279] * v[1766]
//			- v[1268] * v[1773] - v[1265] * v[1778] - v[1274] * v[1783] - v[1269] * v[1786];
//		v[1797] = -(v[1764] * v[851]) - v[1778] * v[853] - v[1750] * v[855];
//		v[3895] = v[1797] * v[291];
//		v[3921] = -v[3890] - v[3891] - v[3892] - v[3893] - v[3894] - v[3895] + (*a4)*v[5733 + i1159];
//		v[1798] = -(v[1269] * v[1674]) - v[1265] * v[1750] - v[1279] * v[1756] - v[1268] * v[1761] - v[1273] * v[1762]
//			- v[1278] * v[1767] - v[1275] * v[1775] - v[1274] * v[1776] - v[1270] * v[1785];
//		v[1799] = v[1801] + v[1774] * v[288];
//		v[1800] = v[1799] + v[1802] + v[3962];
//		v[1803] = v[1801] + v[1802];
//		v[1804] = QBi[0][0] * v[1673] + QBi[2][0] * v[1765] + QBi[1][0] * v[1774];
//		v[1805] = QBi[0][1] * v[1673] + QBi[2][1] * v[1765] + QBi[1][1] * v[1774];
//		v[1806] = QBi[0][2] * v[1673] + QBi[2][2] * v[1765] + QBi[1][2] * v[1774];
//		v[1807] = v[1809] + v[1763] * v[282];
//		v[1808] = v[1807] + v[1810] + v[3961];
//		v[1811] = v[1809] + v[1810];
//		v[1812] = QBi[1][0] * v[1671] + QBi[2][0] * v[1759] + QBi[0][0] * v[1763];
//		v[1813] = QBi[1][1] * v[1671] + QBi[2][1] * v[1759] + QBi[0][1] * v[1763];
//		v[1814] = QBi[1][2] * v[1671] + QBi[2][2] * v[1759] + QBi[0][2] * v[1763];
//		v[1815] = v[1816] + v[1818];
//		v[1817] = v[1816] + v[1754] * v[288];
//		v[1819] = v[1817] + v[1818] + v[3963];
//		v[1820] = QBi[2][0] * v[1669] + QBi[1][0] * v[1754] + QBi[0][0] * v[1757];
//		v[1821] = QBi[2][1] * v[1669] + QBi[1][1] * v[1754] + QBi[0][1] * v[1757];
//		v[1822] = QBi[2][2] * v[1669] + QBi[1][2] * v[1754] + QBi[0][2] * v[1757];
//		v[1823] = v[1633] + v[1705];
//		v[1841] = v[1264] * v[1823];
//		v[1824] = -v[1633] + v[1705];
//		v[1826] = (v[162] * v[1823] + v[1680] * v[251] + v[1825] * v[3851]) / v[253];
//		v[1827] = v[1630] + v[1713];
//		v[1849] = v[1264] * v[1827];
//		v[1828] = -v[1630] + v[1713];
//		v[1845] = v[1260] * v[1828];
//		v[1830] = (v[158] * v[1827] + v[1685] * v[252] + v[1829] * v[3851]) / v[253];
//		v[1832] = (v[153] * v[1828] + v[1691] * v[3723] + v[1831] * v[3851]) / v[253];
//		v[1833] = v[1843] + v[1845];
//		v[4014] = v[1833] / v[253];
//		v[1834] = v[1635] + v[1714];
//		v[1839] = v[1260] * v[1834];
//		v[1835] = -v[1635] + v[1714];
//		v[1836] = v[1847] + v[1849];
//		v[4010] = v[1836] / v[253];
//		v[1837] = v[1256] * v[1696] + v[1839] + v[1841];
//		v[4018] = v[1837] / v[253];
//		v[1840] = v[1837] - v[1841];
//		v[1842] = v[1837] - v[1839];
//		v[4011] = v[1842] / v[253];
//		v[1844] = v[1256] * v[1824] + v[1843];
//		v[1846] = v[1844] + v[1845];
//		v[4012] = v[1846] / v[253];
//		v[1848] = v[1256] * v[1835] + v[1847];
//		v[1850] = v[1848] + v[1849];
//		v[4015] = v[1850] / v[253];
//		v[1851] = v[1659] + v[1779];
//		v[1869] = v[1252] * v[1851];
//		v[1852] = -v[1659] + v[1779];
//		v[1854] = (v[1851] * v[218] + v[1754] * v[277] + v[1853] * v[3852]) / v[279];
//		v[1855] = v[1656] + v[1787];
//		v[1877] = v[1252] * v[1855];
//		v[1856] = -v[1656] + v[1787];
//		v[1873] = v[1248] * v[1856];
//		v[1858] = (v[1855] * v[214] + v[1759] * v[278] + v[1857] * v[3852]) / v[279];
//		v[1860] = (v[1856] * v[209] + v[1765] * v[3731] + v[1859] * v[3852]) / v[279];
//		v[1861] = v[1871] + v[1873];
//		v[4026] = v[1861] / v[279];
//		v[1862] = v[1661] + v[1788];
//		v[1867] = v[1248] * v[1862];
//		v[1863] = -v[1661] + v[1788];
//		v[1864] = v[1875] + v[1877];
//		v[4022] = v[1864] / v[279];
//		v[1865] = v[1244] * v[1770] + v[1867] + v[1869];
//		v[4030] = v[1865] / v[279];
//		v[1868] = v[1865] - v[1869];
//		v[1870] = v[1865] - v[1867];
//		v[4023] = v[1870] / v[279];
//		v[1872] = v[1244] * v[1852] + v[1871];
//		v[1874] = v[1872] + v[1873];
//		v[4024] = v[1874] / v[279];
//		v[1876] = v[1244] * v[1863] + v[1875];
//		v[1878] = v[1876] + v[1877];
//		v[4027] = v[1878] / v[279];
//		v[1880] = ii3862 + v[142] * v[1730] + v[143] * v[1731] + v[144] * v[1732] - v[1804] * v[198] - v[1805] * v[199]
//			- v[1806] * v[200];
//		v[1883] = ii3861 + v[142] * v[1738] + v[143] * v[1739] + v[144] * v[1740] - v[1812] * v[198] - v[1813] * v[199]
//			- v[1814] * v[200];
//		v[1886] = ii3860 + v[142] * v[1746] + v[143] * v[1747] + v[144] * v[1748] - v[1820] * v[198] - v[1821] * v[199]
//			- v[1822] * v[200];
//		v[3853] = v[1880] * v[302] + v[1883] * v[303] + v[1886] * v[304];
//		v[1889] = -(v[1193] * v[1888] * v[3853]);
//		v[1890] = v[3853] / v[727];
//		v[3854] = v[1890] * v[315];
//		v[2089] = v[1890];
//		v[1900] = v[1886] * v[316] + v[304] * v[3854];
//		v[1896] = v[1883] * v[316] + v[303] * v[3854];
//		v[1892] = v[1880] * v[316] + v[302] * v[3854];
//		v[1891] = 2e0*v[1892] * v[317];
//		v[1999] = v[1891];
//		v[1893] = 2e0*v[1224] * v[1892];
//		v[1894] = v[1892];
//		v[1895] = 2e0*v[1896] * v[318];
//		v[2001] = v[1895];
//		v[1897] = 2e0*v[1215] * v[1896];
//		v[1898] = v[1896];
//		v[1936] = v[1898];
//		v[1899] = 2e0*v[1900] * v[319];
//		v[2004] = v[1899];
//		v[1901] = 2e0*v[1214] * v[1900];
//		v[1902] = v[1900];
//		v[1938] = v[1902];
//		v[1903] = 0e0;
//		v[1904] = 0e0;
//		v[1905] = 0e0;
//		v[1906] = 0e0;
//		v[1907] = 0e0;
//		v[1908] = 0e0;
//		v[1909] = 0e0;
//		v[1910] = 0e0;
//		v[1911] = 0e0;
//		v[1912] = 0e0;
//		v[1913] = 0e0;
//		v[1914] = 0e0;
//		v[1915] = 0e0;
//		v[1916] = 0e0;
//		v[1917] = 0e0;
//		v[1918] = 0e0;
//		v[1919] = 0e0;
//		v[1920] = 0e0;
//		v[1921] = 0e0;
//		v[1922] = 0e0;
//		v[1923] = 0e0;
//		v[1924] = 0e0;
//		v[1925] = 0e0;
//		v[1926] = 0e0;
//		v[1927] = 0e0;
//		v[1928] = 0e0;
//		v[1929] = 0e0;
//		v[1930] = 0e0;
//		v[1931] = 0e0;
//		v[1932] = 0e0;
//		v[1933] = 0e0;
//		v[1934] = 0e0;
//		b1935 = b673;
//		if (b1935) {
//			v[1937] = v[1902] * v[327] - v[1936] * v[328];
//			v[1939] = -(v[1938] * v[326]) + v[1894] * v[328];
//			v[1940] = v[1936] * v[326] - v[1894] * v[327];
//			v[3858] = v[1359] * v[1937] + v[1354] * v[1939] + v[1349] * v[1940];
//			v[3856] = v[1937] * v[675] + v[1939] * v[676] + v[1940] * v[677];
//			v[1941] = v[3856] / v[678];
//			v[3857] = v[1941] * v[4223];
//			v[1951] = v[1941] * v[3855];
//			v[1906] = -(v[1362] * v[1942] * v[3856]);
//			v[1943] = v[1937] * v[3759] + v[3857] * v[675];
//			v[1962] = 2e0*v[1943] * v[687];
//			v[1945] = v[1939] * v[3759] + v[3857] * v[676];
//			v[1959] = 2e0*v[1945] * v[688];
//			v[1946] = v[1940] * v[3759] + v[3857] * v[677];
//			v[1960] = 2e0*v[1946] * v[689];
//			v[1909] = v[1951] * v[3790] + v[3858] * v[680];
//			v[1908] = v[1941] * v[1949] * v[684];
//			v[1907] = v[1941] * v[3790] * v[685] + v[3858] * v[686];
//			v[1933] = 2e0*v[1951] * v[1952] * v[4224];
//			v[1934] = v[1361] * v[1941] * v[1947] * v[1952] * v[4225];
//			v[1905] = v[1349] * v[3857] + v[1940] * v[3859];
//			v[1904] = v[1354] * v[3857] + v[1939] * v[3859];
//			v[1903] = v[1359] * v[3857] + v[1937] * v[3859];
//			v[1954] = (v[1945] * v[687] + v[1943] * v[688]) / 2e0;
//			v[1955] = v[1959] + v[1962];
//			v[1956] = v[1955] + v[1960];
//			v[1957] = (v[1946] * v[687] + v[1943] * v[689]) / 2e0;
//			v[1958] = (v[1946] * v[688] + v[1945] * v[689]) / 2e0;
//			v[1961] = v[1959] + v[1960];
//			v[1963] = v[1960] + v[1962];
//			v[1912] = (v[1334] * v[1943] + v[1333] * v[1945] + 4e0*v[1946] * v[1964]) / 2e0;
//			v[1911] = (v[1335] * v[1943] + v[1333] * v[1946] + 4e0*v[1945] * v[1965]) / 2e0;
//			v[1910] = (v[1335] * v[1945] + v[1334] * v[1946] + 4e0*v[1943] * v[1966]) / 2e0;
//			v[1967] = -4e0*v[1956] * v[2312];
//			v[1932] = 8e0*v[1345] * v[1956] * v[4226];
//			v[1923] = (v[1289] * v[1967]) / 2e0;
//			v[1927] = (v[1293] * v[1967]) / 2e0;
//			v[1925] = v[1291] * v[1967];
//			v[1929] = v[1295] * v[1967];
//			v[1928] = v[1294] * v[1967];
//			v[1930] = v[1296] * v[1967];
//			v[1924] = v[1290] * v[1967];
//			v[1926] = v[1292] * v[1967];
//			v[1931] = (v[1297] * v[1967]) / 2e0;
//			v[1968] = -v[1946] + v[1954];
//			v[1969] = v[1946] + v[1954];
//			v[1970] = v[1945] + v[1957];
//			v[1971] = -v[1945] + v[1957];
//			v[1972] = -v[1943] + v[1958];
//			v[1973] = v[1943] + v[1958];
//			v[1915] = v[1337] * v[1967] + v[1968] * v[690];
//			v[1917] = v[1339] * v[1967] + v[1969] * v[690];
//			v[1916] = v[1338] * v[1967] + v[1970] * v[690];
//			v[1920] = v[1342] * v[1967] + v[1971] * v[690];
//			v[1914] = (v[1336] * v[1967] - v[1961] * v[690]) / 2e0;
//			v[1919] = v[1341] * v[1967] + v[1972] * v[690];
//			v[1921] = v[1343] * v[1967] + v[1973] * v[690];
//			v[1918] = (v[1340] * v[1967] - v[1963] * v[690]) / 2e0;
//			v[1922] = (v[1344] * v[1967] - v[1955] * v[690]) / 2e0;
//			v[1913] = -(v[1297] * v[1955]) / 2e0 - (v[1289] * v[1961]) / 2e0 - (v[1293] * v[1963]) / 2e0 + v[1290] * v[1968]
//				+ v[1292] * v[1969] + v[1291] * v[1970] + v[1295] * v[1971] + v[1294] * v[1972] + v[1296] * v[1973];
//		}
//		else {
//		};
//		v[1998] = v[1902];
//		v[1992] = v[1898];
//		v[1990] = v[1894];
//		v[1988] = v[1915];
//		v[1987] = v[1916];
//		v[1985] = v[1918];
//		v[1984] = v[1919];
//		v[1982] = v[1921];
//		v[1981] = v[1922];
//		v[1974] = 0e0;
//		v[1975] = 0e0;
//		v[1976] = 0e0;
//		v[1977] = 0e0;
//		v[1978] = 0e0;
//		v[1979] = 0e0;
//		b1980 = (*previouscontact);
//		if (b1980) {
//			v[1996] = v[1189] * v[1990];
//			v[1995] = v[1187] * v[1902];
//			v[1993] = v[1188] * v[1898];
//			v[1922] = 0e0;
//			v[1921] = 0e0;
//			v[1983] = ii3860 + v[145] * v[1746] + v[146] * v[1747] + v[147] * v[1748] + gti[0] * v[1920] + gti[2] * v[1981]
//				+ gti[1] * v[1982] - v[1820] * v[201] - v[1821] * v[202] - v[1822] * v[203];
//			v[3863] = -(v[1983] * v[319]) - v[1998] * v[710];
//			v[1920] = 0e0;
//			v[1919] = 0e0;
//			v[1918] = 0e0;
//			v[1986] = ii3861 + v[145] * v[1738] + v[146] * v[1739] + v[147] * v[1740] + gti[0] * v[1917] + gti[2] * v[1984]
//				+ gti[1] * v[1985] - v[1812] * v[201] - v[1813] * v[202] - v[1814] * v[203];
//			v[3865] = -(v[1986] * v[318]) - v[1992] * v[709];
//			v[1917] = 0e0;
//			v[1916] = 0e0;
//			v[1915] = 0e0;
//			v[1989] = ii3862 + v[145] * v[1730] + v[146] * v[1731] + v[147] * v[1732] + gti[0] * v[1914] + gti[2] * v[1987]
//				+ gti[1] * v[1988] - v[1804] * v[201] - v[1805] * v[202] - v[1806] * v[203];
//			v[3864] = -(v[1989] * v[317]) - v[1990] * v[708];
//			v[1914] = 0e0;
//			v[1991] = v[1993] + v[1996];
//			v[1994] = v[1993] + v[1995];
//			v[1997] = v[1995] + v[1996];
//			v[1979] = v[1187] * v[1983];
//			v[1978] = v[1188] * v[1986];
//			v[1977] = v[1189] * v[1989];
//			v[1974] = -(v[1189] * v[1999]) + v[1990] * v[2000] - v[1994] * v[317];
//			v[1893] = v[1893] + v[1989] * v[2000] + v[1189] * (v[3863] + v[3865]) - v[1994] * v[708];
//			v[1975] = -(v[1188] * v[2001]) + v[1992] * v[2003] - v[1997] * v[318];
//			v[1897] = v[1897] + v[1986] * v[2003] + v[1188] * (v[3863] + v[3864]) - v[1997] * v[709];
//			v[1976] = -(v[1187] * v[2004]) + v[1998] * v[2007] - v[1991] * v[319];
//			v[1901] = v[1901] + v[1983] * v[2007] + v[1187] * (v[3864] + v[3865]) - v[1991] * v[710];
//		}
//		else {
//		};
//		v[2674] = v[1192] * v[1899];
//		v[2667] = v[1191] * v[1895];
//		v[3866] = v[1192] * v[1902];
//		v[2085] = v[1902];
//		v[2082] = v[301] * v[3866];
//		v[2079] = v[291] * v[3866];
//		v[2076] = v[281] * v[3866];
//		v[2073] = v[275] * v[3866];
//		v[2070] = v[265] * v[3866];
//		v[2067] = v[255] * v[3866];
//		v[3867] = v[1191] * v[1898];
//		v[2083] = v[301] * v[3867];
//		v[2080] = v[291] * v[3867];
//		v[2077] = v[281] * v[3867];
//		v[2074] = v[275] * v[3867];
//		v[2071] = v[265] * v[3867];
//		v[2068] = v[255] * v[3867];
//		v[2049] = v[1898];
//		v[2661] = v[1190] * v[1891];
//		v[2041] = v[1894];
//		v[3868] = v[1190] * v[2041];
//		v[2062] = v[301] * v[3868];
//		v[2060] = v[291] * v[3868];
//		v[2058] = v[281] * v[3868];
//		v[2056] = v[275] * v[3868];
//		v[2054] = v[265] * v[3868];
//		v[2052] = v[255] * v[3868];
//		v[2008] = v[3869] + v[3870] + v[3871] + v[3872] + v[3873] + v[3874] + (*a4)*v[6341 + i1159];
//		v[3883] = v[2008] * v[319];
//		v[3889] = v[2085] * v[2090] - v[3883];
//		v[2009] = v[3875] + v[3876] + v[3877] + v[3878] + v[3879] + v[3880] + (*a4)*v[6329 + i1159];
//		v[3886] = -(v[2009] * v[318]);
//		v[3888] = v[2049] * v[2095] + v[3886];
//		v[2010] = v[3890] + v[3891] + v[3892] + v[3893] + v[3894] + v[3895] + (*a4)*v[6317 + i1159];
//		v[3881] = -(v[2010] * v[317]);
//		v[3885] = v[2041] * v[2100] + v[3881];
//		v[2015] = v[1706] * v[2011] + v[1707] * v[2012] + v[1823] * v[2013] + v[1827] * v[2014] + v[1826] * v[262]
//			+ v[1830] * v[264] + v[1726] * v[2682] + v[1733] * v[3724] + v[1741] * v[3727] + (*a4)*v[6245 + i1159];
//		v[2020] = v[1572] * v[1725] + v[1571] * v[1743] + v[1706] * v[2016] + v[1715] * v[2017] + v[1828] * v[2018]
//			+ v[1834] * v[2019] + v[1826] * v[256] + v[1832] * v[264] + v[1734] * v[2683] + (*a4)*v[6257 + i1159];
//		v[2025] = v[1574] * v[1729] + v[1573] * v[1737] + v[1706] * v[2021] + v[1696] * v[2022] + v[1824] * v[2023]
//			+ v[1835] * v[2024] + v[1830] * v[256] + v[1832] * v[262] + v[1745] * v[2684] + (*a4)*v[6269 + i1159];
//		v[2030] = v[1780] * v[2026] + v[1781] * v[2027] + v[1851] * v[2028] + v[1855] * v[2029] + v[1800] * v[2685]
//			+ v[1854] * v[288] + v[1858] * v[290] + v[1807] * v[3732] + v[1815] * v[3735] + (*a4)*v[6281 + i1159];
//		v[2035] = v[1578] * v[1799] + v[1577] * v[1817] + v[1780] * v[2031] + v[1789] * v[2032] + v[1856] * v[2033]
//			+ v[1862] * v[2034] + v[1808] * v[2686] + v[1854] * v[282] + v[1860] * v[290] + (*a4)*v[6293 + i1159];
//		v[2040] = v[1580] * v[1803] + v[1579] * v[1811] + v[1780] * v[2036] + v[1770] * v[2037] + v[1852] * v[2038]
//			+ v[1863] * v[2039] + v[1819] * v[2687] + v[1858] * v[282] + v[1860] * v[288] + (*a4)*v[6305 + i1159];
//		v[2042] = v[3867] + v[3868];
//		v[2043] = v[2052] + v[2068];
//		v[2044] = v[2054] + v[2071];
//		v[2045] = v[2056] + v[2074];
//		v[2046] = v[2058] + v[2077];
//		v[2047] = v[2060] + v[2080];
//		v[2048] = v[2062] + v[2083];
//		v[2051] = v[3866] + v[3868];
//		v[2053] = v[2052] + v[2067];
//		v[2055] = v[2054] + v[2070];
//		v[2057] = v[2056] + v[2073];
//		v[2059] = v[2058] + v[2076];
//		v[2061] = v[2060] + v[2079];
//		v[2063] = v[2062] + v[2082];
//		v[2066] = v[3866] + v[3867];
//		v[2069] = v[2067] + v[2068];
//		v[2072] = v[2070] + v[2071];
//		v[2075] = v[2073] + v[2074];
//		v[2078] = v[2076] + v[2077];
//		v[2081] = v[2079] + v[2080];
//		v[2084] = v[2082] + v[2083];
//		v[2086] = 0e0;
//		b2087 = b728;
//		if (b2087) {
//			b2088 = b730;
//			if (b2088) {
//				v[2086] = -v[2089];
//				v[1894] = 0e0;
//				v[1898] = 0e0;
//				v[1902] = 0e0;
//			}
//			else {
//			};
//		}
//		else {
//		};
//		v[2094] = v[1263] * v[2015] + v[1259] * v[2020] + v[1255] * v[2025] + v[1251] * v[2030] + v[1247] * v[2035]
//			+ v[1243] * v[2040] + v[1899] * v[2090] + v[2041] * v[2091] + v[2049] * v[2092] + v[2085] * v[2093] + v[319] * (v[3881]
//				+ v[3886]) + v[3919] * v[717];
//		v[2099] = v[1262] * v[2015] + v[1258] * v[2020] + v[1254] * v[2025] + v[1250] * v[2030] + v[1246] * v[2035]
//			+ v[1242] * v[2040] + v[1895] * v[2095] + v[2041] * v[2096] + v[2049] * v[2097] + v[2085] * v[2098] - v[318] * (-v[3881]
//				+ v[3883]) + v[3920] * v[715];
//		v[2119] = v[2099];
//		v[1897] = v[1897] - v[2009] * v[3882] + v[1191] * (v[3885] + v[3889]);
//		v[2238] = v[1897];
//		v[1901] = v[1901] - v[2008] * v[3884] + v[1192] * (v[3885] + v[3888]);
//		v[2237] = v[1901];
//		v[1893] = v[1893] - v[2010] * v[3887] + v[1190] * (v[3888] + v[3889]);
//		v[2239] = v[1893];
//		v[2104] = v[1261] * v[2015] + v[1257] * v[2020] + v[1253] * v[2025] + v[1249] * v[2030] + v[1245] * v[2035]
//			+ v[1241] * v[2040] + v[1891] * v[2100] + v[2041] * v[2101] + v[2049] * v[2102] + v[2085] * v[2103] - v[317] * (v[3883]
//				- v[3886]) + v[3921] * v[713];
//		v[2118] = v[2104];
//		v[2105] = 0e0;
//		v[2106] = 0e0;
//		v[2107] = 0e0;
//		v[2108] = 0e0;
//		v[2109] = 0e0;
//		v[2110] = 0e0;
//		v[2111] = 0e0;
//		v[2112] = 0e0;
//		b2113 = b728;
//		if (b2113) {
//			v[2123] = v[1198] * v[2094] + v[1196] * v[2118] + v[1197] * v[2119];
//			b2114 = b730;
//			if (b2114) {
//				v[2117] = v[2086] * v[2224] * (*zetan);
//				v[2116] = v[2117] * v[2225];
//				v[2111] = -((v[1208] * v[2116]) / v[1202]);
//				v[2108] = v[1208] * v[2115] * v[2117] * v[4257];
//				v[2086] = 0e0;
//				v[2104] = 0e0;
//				v[2099] = 0e0;
//				v[2109] = v[2123];
//				v[2094] = 0e0;
//			}
//			else {
//				v[2121] = v[1890] * v[2230] * (*zetan);
//				v[2116] = v[2121] * v[2231];
//				v[2112] = -((v[1208] * v[2116]) / v[1207]);
//				v[1889] = v[1889] + v[1208] * v[2120] * v[2121] * v[4262];
//				v[1890] = 0e0;
//				v[2104] = 0e0;
//				v[2099] = 0e0;
//				v[2110] = v[2123];
//				v[2094] = 0e0;
//			};
//			v[2107] = v[1198] * v[2116];
//			v[2106] = v[1197] * v[2116];
//			v[2105] = v[1196] * v[2116];
//		}
//		else {
//		};
//		v[2229] = v[2105];
//		v[2227] = v[2106];
//		v[2226] = v[2107];
//		v[2124] = (*ct)*(v[1161] * (v[1730] * v[337] + v[1731] * v[339] + v[1732] * v[341]) + v[1162] * (v[1738] * v[337]
//			+ v[1739] * v[339] + v[1740] * v[341]) + v[1163] * (v[1746] * v[337] + v[1747] * v[339] + v[1748] * v[341]));
//		v[2125] = (*ct)*(v[1161] * (v[1730] * v[329] + v[1731] * v[330] + v[1732] * v[331]) + v[1162] * (v[1738] * v[329]
//			+ v[1739] * v[330] + v[1740] * v[331]) + v[1163] * (v[1746] * v[329] + v[1747] * v[330] + v[1748] * v[331]));
//		v[2126] = (*ct)*(v[1161] * (v[1804] * v[354] + v[1805] * v[356] + v[1806] * v[358]) + v[1162] * (v[1812] * v[354]
//			+ v[1813] * v[356] + v[1814] * v[358]) + v[1163] * (v[1820] * v[354] + v[1821] * v[356] + v[1822] * v[358]));
//		v[2127] = (*ct)*(v[1161] * (v[1804] * v[346] + v[1805] * v[347] + v[1806] * v[348]) + v[1162] * (v[1812] * v[346]
//			+ v[1813] * v[347] + v[1814] * v[348]) + v[1163] * (v[1820] * v[346] + v[1821] * v[347] + v[1822] * v[348]));
//		v[2128] = -(i3896*v[1110]) + i3899 * v[1110] - i3897 * v[1114] + i3900 * v[1114] - i3898 * v[1118] + i3901 * v[1118]
//			- i3838 * v[1122] - i3839 * v[1126] - i3837 * v[1130] - i3844 * v[1146] - i3845 * v[1150] - i3843 * v[1154];
//		v[2144] = v[2128];
//		v[2129] = -(i3896*v[1111]) + i3899 * v[1111] - i3897 * v[1115] + i3900 * v[1115] - i3898 * v[1119] + i3901 * v[1119]
//			- i3838 * v[1123] - i3839 * v[1127] - i3837 * v[1131] - i3844 * v[1147] - i3845 * v[1151] - i3843 * v[1155];
//		v[2142] = v[2129];
//		v[2130] = -(i3896*v[1112]) + i3899 * v[1112] - i3897 * v[1116] + i3900 * v[1116] - i3898 * v[1120] + i3901 * v[1120]
//			- i3838 * v[1124] - i3839 * v[1128] - i3837 * v[1132] - i3844 * v[1148] - i3845 * v[1152] - i3843 * v[1156];
//		v[2140] = v[2130];
//		v[2131] = 0e0;
//		v[2132] = 0e0;
//		v[2133] = 0e0;
//		v[2134] = 0e0;
//		v[2135] = 0e0;
//		v[2136] = 0e0;
//		b2137 = b728;
//		if (b2137) {
//			b2138 = (*stick);
//			if (b2138) {
//				b2139 = b766;
//				if (b2139) {
//					v[2136] = v[2130];
//					v[2130] = 0e0;
//					v[2135] = v[2129];
//					v[2129] = 0e0;
//					v[2134] = v[2128];
//					v[2128] = 0e0;
//				}
//				else {
//					v[3902] = (*mud)*(v[2144] * v[761] + v[2142] * v[762] + v[2140] * v[763]);
//					v[2130] = 0e0;
//					v[2129] = 0e0;
//					v[3904] = v[2155] * v[3902] * v[777];
//					v[2128] = 0e0;
//					v[3903] = v[2153] * v[3902] * v[776] * v[781];
//					v[2136] = v[2140] * v[2152] + v[3903] * v[763];
//					v[2135] = v[2142] * v[2152] + v[3903] * v[762];
//					v[2134] = v[2144] * v[2152] + v[3903] * v[761];
//					v[2133] = v[3904] * v[754];
//					v[2132] = v[3904] * v[753];
//					v[2131] = v[3904] * v[752];
//				};
//			}
//			else {
//				b2156 = b800;
//				if (b2156) {
//					v[2136] = v[2140];
//					v[2130] = 0e0;
//					v[2135] = v[2142];
//					v[2129] = 0e0;
//					v[2134] = v[2144];
//					v[2128] = 0e0;
//				}
//				else {
//					v[2162] = v[2144] * v[3905];
//					v[2160] = v[2142] * v[3905];
//					v[2159] = v[2140] * v[3905];
//					v[2130] = 0e0;
//					v[2129] = 0e0;
//					v[3907] = (*mud)*v[2165] * (v[2144] * v[761] + v[2142] * v[762] + v[2140] * v[763])*v[807];
//					v[2128] = 0e0;
//					v[3906] = v[2161] * (v[2162] * v[761] + v[2160] * v[762] + v[2159] * v[763])*v[806];
//					v[2136] = v[3906] * v[763] + v[2159] * v[807];
//					v[2135] = v[3906] * v[762] + v[2160] * v[807];
//					v[2134] = v[3906] * v[761] + v[2162] * v[807];
//					v[2133] = v[3907] * v[754];
//					v[2132] = v[3907] * v[753];
//					v[2131] = v[3907] * v[752];
//				};
//			};
//		}
//		else {
//		};
//		v[3912] = -((*ct)*v[2134]);
//		v[3910] = -((*ct)*v[2135]);
//		v[3908] = -((*ct)*v[2136]);
//		v[2167] = -((*ct)*(v[1163] * v[2040] + v[2136] * v[301])) - i3843 * v[770];
//		v[2168] = -((*ct)*(v[1163] * v[2035] + v[2136] * v[291])) - i3845 * v[770];
//		v[2169] = -((*ct)*(v[1163] * v[2030] + v[2136] * v[281])) - i3844 * v[770];
//		v[2174] = -((*ct)*(v[1163] * v[2025] + v[2136] * v[275])) - i3837 * v[770];
//		v[2175] = -((*ct)*(v[1163] * v[2020] + v[2136] * v[265])) - i3839 * v[770];
//		v[2176] = -((*ct)*(v[1163] * v[2015] + v[2136] * v[255])) - i3838 * v[770];
//		v[3940] = -(v[246] * v[3908]) + v[249] * v[3908] - i3898 * v[3909] + i3901 * v[3909];
//		v[3941] = -(v[245] * v[3908]) + v[248] * v[3908] - i3897 * v[3909] + i3900 * v[3909];
//		v[3942] = -(v[244] * v[3908]) + v[247] * v[3908] - i3896 * v[3909] + i3899 * v[3909];
//		v[2180] = -((*ct)*(v[1162] * v[2040] + v[2135] * v[301])) - i3843 * v[769];
//		v[2181] = -((*ct)*(v[1162] * v[2035] + v[2135] * v[291])) - i3845 * v[769];
//		v[2182] = -((*ct)*(v[1162] * v[2030] + v[2135] * v[281])) - i3844 * v[769];
//		v[2187] = -((*ct)*(v[1162] * v[2025] + v[2135] * v[275])) - i3837 * v[769];
//		v[2188] = -((*ct)*(v[1162] * v[2020] + v[2135] * v[265])) - i3839 * v[769];
//		v[2189] = -((*ct)*(v[1162] * v[2015] + v[2135] * v[255])) - i3838 * v[769];
//		v[3937] = -(v[246] * v[3910]) + v[249] * v[3910] - i3898 * v[3911] + i3901 * v[3911];
//		v[3938] = -(v[245] * v[3910]) + v[248] * v[3910] - i3897 * v[3911] + i3900 * v[3911];
//		v[3939] = -(v[244] * v[3910]) + v[247] * v[3910] - i3896 * v[3911] + i3899 * v[3911];
//		v[2193] = -((*ct)*(v[1161] * v[2040] + v[2134] * v[301])) - i3843 * v[768];
//		v[2194] = -((*ct)*(v[1161] * v[2035] + v[2134] * v[291])) - i3845 * v[768];
//		v[2195] = -((*ct)*(v[1161] * v[2030] + v[2134] * v[281])) - i3844 * v[768];
//		v[2200] = -((*ct)*(v[1161] * v[2025] + v[2134] * v[275])) - i3837 * v[768];
//		v[2201] = -((*ct)*(v[1161] * v[2020] + v[2134] * v[265])) - i3839 * v[768];
//		v[2202] = -((*ct)*(v[1161] * v[2015] + v[2134] * v[255])) - i3838 * v[768];
//		v[3934] = -(v[246] * v[3912]) + v[249] * v[3912] - i3898 * v[3913] + i3901 * v[3913];
//		v[3935] = -(v[245] * v[3912]) + v[248] * v[3912] - i3897 * v[3913] + i3900 * v[3913];
//		v[3936] = -(v[244] * v[3912]) + v[247] * v[3912] - i3896 * v[3913] + i3899 * v[3913];
//		v[2206] = (*epst)*v[2136];
//		v[2207] = (*epst)*v[2135];
//		v[2208] = (*epst)*v[2134];
//		v[2209] = v[2133];
//		v[2210] = v[2133];
//		v[2211] = v[2132];
//		v[2212] = v[2132];
//		v[2213] = v[2131];
//		v[2214] = v[2131];
//		b2215 = b728;
//		if (b2215) {
//			v[2216] = 0e0;
//			v[2217] = 0e0;
//			v[2218] = 0e0;
//			b2219 = b747;
//			if (b2219) {
//				v[2218] = v[2209];
//				v[2209] = 0e0;
//				v[2217] = v[2211];
//				v[2211] = 0e0;
//				v[2216] = v[2213];
//				v[2213] = 0e0;
//			}
//			else {
//				v[2210] = 0e0;
//				v[2209] = 0e0;
//				v[2212] = 0e0;
//				v[2211] = 0e0;
//				v[2214] = 0e0;
//				v[2213] = 0e0;
//			};
//			v[2228] = v[2216] * v[714] + v[2217] * v[716] + v[2218] * v[718];
//			b2223 = b730;
//			if (b2223) {
//				v[2107] = v[2107] + v[2218] * v[742];
//				v[2106] = v[2106] + v[2217] * v[742];
//				v[2109] = v[2109] + v[2228];
//				v[2105] = v[2105] + v[2216] * v[742];
//				v[2111] = v[2111] + 2e0*v[2109] * (*zetan);
//				v[2108] = v[2108] + (v[2111] * v[2224] * v[2225]) / 2e0;
//			}
//			else {
//				v[2107] = v[2226] + v[2218] * v[746];
//				v[2106] = v[2227] + v[2217] * v[746];
//				v[2110] = v[2110] + v[2228];
//				v[2105] = v[2229] + v[2216] * v[746];
//				v[2112] = v[2112] + 2e0*v[2110] * (*zetan);
//				v[1889] = v[1889] + (v[2112] * v[2230] * v[2231]) / 2e0;
//			};
//		}
//		else {
//		};
//		v[3914] = v[2105] * v[317];
//		v[3915] = v[2106] * v[318];
//		v[3916] = v[2107] * v[319];
//		v[4006] = (*ct)*(-(v[1118] * v[2134]) - v[1119] * v[2135] - v[1120] * v[2136]) + v[2674] + v[2042] * v[319]
//			+ v[2085] * v[3884] + v[319] * v[3914] + v[319] * v[3915] - v[2125] * v[623] - v[2124] * v[638] + v[2127] * v[653]
//			+ v[2126] * v[668] + v[2107] * v[717];
//		v[4005] = (*ct)*(-(v[1114] * v[2134]) - v[1115] * v[2135] - v[1116] * v[2136]) + v[2667] + v[2051] * v[318]
//			+ v[2049] * v[3882] + v[318] * v[3914] + v[318] * v[3916] - v[2125] * v[621] - v[2124] * v[636] + v[2127] * v[651]
//			+ v[2126] * v[666] + v[2106] * v[715];
//		v[4004] = (*ct)*(-(v[1110] * v[2134]) - v[1111] * v[2135] - v[1112] * v[2136]) + v[2661] + v[2066] * v[317]
//			+ v[2041] * v[3887] + v[317] * v[3915] + v[317] * v[3916] - v[2125] * v[619] - v[2124] * v[634] + v[2127] * v[649]
//			+ v[2126] * v[664] + v[2105] * v[713];
//		v[2240] = v[1889];
//		v[2236] = v[2214];
//		v[2235] = v[2212];
//		v[2234] = v[2210];
//		b2232 = b728;
//		if (b2232) {
//			v[3918] = -(v[2236] * v[317]) - v[2235] * v[318] - v[2234] * v[319];
//			b2233 = b730;
//			if (b2233) {
//				v[1901] = v[1901] + v[2210] * v[3917];
//				v[2210] = 0e0;
//				v[1897] = v[1897] + v[2212] * v[3917];
//				v[2212] = 0e0;
//				v[1893] = v[1893] + v[2214] * v[3917];
//				v[2214] = 0e0;
//				v[2108] = v[2108] + v[3767] * v[3918] * v[4300];
//				v[1889] = v[1889] - v[2108];
//			}
//			else {
//				v[1901] = v[2237] + v[2234] * v[736];
//				v[2210] = 0e0;
//				v[1897] = v[2238] + v[2235] * v[736];
//				v[2212] = 0e0;
//				v[1893] = v[2239] + v[2236] * v[736];
//				v[2214] = 0e0;
//				v[1889] = v[2240] + (*n2)*v[3918] * v[4304] * v[724];
//			};
//		}
//		else {
//		};
//		v[2241] = v[1192] * v[2040] + v[2107] * v[301];
//		v[3959] = v[2241] * v[319];
//		v[2242] = v[1192] * v[2035] + v[2107] * v[291];
//		v[3957] = v[2242] * v[319];
//		v[2243] = v[1192] * v[2030] + v[2107] * v[281];
//		v[3955] = v[2243] * v[319];
//		v[2244] = v[1192] * v[2025] + v[2107] * v[275];
//		v[3953] = v[2244] * v[319];
//		v[2245] = v[1192] * v[2020] + v[2107] * v[265];
//		v[3951] = v[2245] * v[319];
//		v[2246] = v[1192] * v[2015] + v[2107] * v[255];
//		v[3949] = v[2246] * v[319];
//		v[2247] = v[2107] * v[3335] + v[1192] * v[3919];
//		v[2248] = v[1191] * v[2040] + v[2106] * v[301];
//		v[3960] = v[2248] * v[318];
//		v[2249] = v[1191] * v[2035] + v[2106] * v[291];
//		v[3958] = v[2249] * v[318];
//		v[2250] = v[1191] * v[2030] + v[2106] * v[281];
//		v[3956] = v[2250] * v[318];
//		v[2251] = v[1191] * v[2025] + v[2106] * v[275];
//		v[3954] = v[2251] * v[318];
//		v[2252] = v[1191] * v[2020] + v[2106] * v[265];
//		v[3952] = v[2252] * v[318];
//		v[2253] = v[1191] * v[2015] + v[2106] * v[255];
//		v[3950] = v[2253] * v[318];
//		v[2254] = v[2106] * v[3354] + v[1191] * v[3920];
//		v[2255] = v[1190] * v[2040] + v[2105] * v[301];
//		v[3948] = v[2255] * v[317];
//		v[2256] = v[1190] * v[2035] + v[2105] * v[291];
//		v[3947] = v[2256] * v[317];
//		v[2257] = v[1190] * v[2030] + v[2105] * v[281];
//		v[3946] = v[2257] * v[317];
//		v[2258] = v[1190] * v[2025] + v[2105] * v[275];
//		v[3945] = v[2258] * v[317];
//		v[2259] = v[1190] * v[2020] + v[2105] * v[265];
//		v[3944] = v[2259] * v[317];
//		v[2260] = v[1190] * v[2015] + v[2105] * v[255];
//		v[3943] = v[2260] * v[317];
//		v[1901] = v[1901] + v[2107] * v[2265] + v[2105] * v[3333] + v[2106] * v[3334];
//		v[1897] = v[1897] + v[2106] * v[2268] + v[2105] * v[3352] + v[2107] * v[3353];
//		v[2271] = v[2105] * v[3370] + v[1190] * v[3921];
//		v[1893] = v[1893] + v[2105] * v[2272] + v[2106] * v[3368] + v[2107] * v[3369];
//		v[2277] = 0e0;
//		v[2278] = 0e0;
//		v[2279] = 0e0;
//		v[2280] = 0e0;
//		v[2281] = 0e0;
//		v[2282] = 0e0;
//		v[2283] = 0e0;
//		v[2284] = 0e0;
//		v[2285] = 0e0;
//		b2286 = (*previouscontact);
//		if (b2286) {
//			v[3924] = v[2208] * v[317];
//			v[3923] = v[2207] * v[318];
//			v[3925] = v[3923] + v[3924];
//			v[3922] = v[2206] * v[319];
//			v[3927] = v[3922] + v[3923];
//			v[3926] = v[3922] + v[3924];
//			v[1979] = v[1979] + v[2206] * v[710];
//			v[1978] = v[1978] + v[2207] * v[709];
//			v[1974] = v[1974] + v[1299] * v[2208] - v[317] * v[3927];
//			v[1975] = v[1975] + v[1301] * v[2207] - v[318] * v[3926];
//			v[1976] = v[1976] + v[1303] * v[2206] - v[319] * v[3925];
//			v[1977] = v[1977] + v[2208] * v[708];
//			v[1901] = v[1901] + v[2206] * v[3787] - v[3925] * v[710];
//			v[1897] = v[1897] + v[2207] * v[3788] - v[3926] * v[709];
//			v[1893] = v[1893] - v[2208] * v[3789] - v[3927] * v[708];
//			v[2277] = gti[0] * v[1974];
//			v[2278] = gti[1] * v[1974];
//			v[2279] = gti[2] * v[1974];
//			v[2287] = -v[1974];
//			v[2288] = -(v[1974] * v[203]);
//			v[2289] = -(v[1974] * v[202]);
//			v[2290] = -(v[1974] * v[201]);
//			v[2291] = v[1974];
//			v[2292] = v[147] * v[1974];
//			v[2293] = v[146] * v[1974];
//			v[2294] = v[145] * v[1974];
//			v[2280] = gti[0] * v[1975];
//			v[2281] = gti[1] * v[1975];
//			v[2282] = gti[2] * v[1975];
//			v[2295] = -v[1975];
//			v[2296] = -(v[1975] * v[203]);
//			v[2297] = -(v[1975] * v[202]);
//			v[2298] = -(v[1975] * v[201]);
//			v[2299] = v[1975];
//			v[2300] = v[147] * v[1975];
//			v[2301] = v[146] * v[1975];
//			v[2302] = v[145] * v[1975];
//			v[2283] = gti[0] * v[1976];
//			v[2284] = gti[1] * v[1976];
//			v[2285] = gti[2] * v[1976];
//			v[2303] = -v[1976];
//			v[2304] = -(v[1976] * v[203]);
//			v[2305] = -(v[1976] * v[202]);
//			v[2306] = -(v[1976] * v[201]);
//			v[2307] = v[1976];
//			v[2308] = v[147] * v[1976];
//			v[2309] = v[146] * v[1976];
//			v[2310] = v[145] * v[1976];
//			v[2271] = -v[1977] + v[2271];
//			v[2254] = -v[1978] + v[2254];
//			v[2247] = -v[1979] + v[2247];
//		}
//		else {
//			v[2294] = 0e0;
//			v[2293] = 0e0;
//			v[2292] = 0e0;
//			v[2302] = 0e0;
//			v[2301] = 0e0;
//			v[2300] = 0e0;
//			v[2310] = 0e0;
//			v[2309] = 0e0;
//			v[2308] = 0e0;
//			v[2291] = 0e0;
//			v[2299] = 0e0;
//			v[2307] = 0e0;
//			v[2290] = 0e0;
//			v[2289] = 0e0;
//			v[2288] = 0e0;
//			v[2298] = 0e0;
//			v[2297] = 0e0;
//			v[2296] = 0e0;
//			v[2306] = 0e0;
//			v[2305] = 0e0;
//			v[2304] = 0e0;
//			v[2287] = 0e0;
//			v[2295] = 0e0;
//			v[2303] = 0e0;
//		};
//		b2311 = b673;
//		if (b2311) {
//			v[1931] = v[1931] + (v[2285] * v[690]) / 2e0;
//			v[1930] = v[1930] + v[2284] * v[690];
//			v[1929] = v[1929] + v[2283] * v[690];
//			v[1928] = v[1928] + v[2282] * v[690];
//			v[1927] = v[1927] + (v[2281] * v[690]) / 2e0;
//			v[1926] = v[1926] + v[2280] * v[690];
//			v[1925] = v[1925] + v[2279] * v[690];
//			v[1924] = v[1924] + v[2278] * v[690];
//			v[1913] = v[1913] + (v[1336] * v[2277]) / 2e0 + v[1337] * v[2278] + v[1338] * v[2279] + v[1339] * v[2280] +
//				(v[1340] * v[2281]) / 2e0 + v[1341] * v[2282] + v[1342] * v[2283] + v[1343] * v[2284] + (v[1344] * v[2285]) / 2e0;
//			v[1923] = v[1923] + (v[2277] * v[690]) / 2e0;
//			v[1932] = v[1932] - 4e0*v[1913] * v[2312];
//			v[3928] = v[1923] - v[1932];
//			v[3930] = (v[1925] + v[1929]) / 2e0;
//			v[3929] = (v[1928] + v[1930]) / 2e0;
//			v[1912] = v[1912] - v[1924] + v[1926] + v[3930] * v[687] + v[3929] * v[688] - 2e0*(v[1927] + v[3928])*v[689];
//			v[3931] = (v[1924] + v[1926]) / 2e0;
//			v[1911] = v[1911] + v[1925] - v[1929] + v[3931] * v[687] - 2e0*(v[1931] + v[3928])*v[688] + v[3929] * v[689];
//			v[1910] = v[1910] - v[1928] + v[1930] - 2e0*(v[1927] + v[1931] - v[1932])*v[687] + v[3931] * v[688]
//				+ v[3930] * v[689];
//			v[3932] = v[1910] * v[675] + v[1911] * v[676] + v[1912] * v[677];
//			v[1909] = v[1909] + v[3932] * v[680];
//			v[1907] = v[1907] + v[3932] * v[686];
//			v[1908] = v[1908] + v[1909] * v[685];
//			v[1933] = v[1933] + 2e0*v[1907] * v[1947];
//			v[1934] = v[1934] + (v[1933] * v[1948]) / 2e0;
//			v[1906] = v[1906] + v[1908] + v[1934];
//			v[3933] = v[1906] / v[678];
//			v[1905] = v[1905] + v[1912] * v[3759] + v[3933] * v[677];
//			v[1904] = v[1904] + v[1911] * v[3759] + v[3933] * v[676];
//			v[1903] = v[1903] + v[1910] * v[3759] + v[3933] * v[675];
//			v[1893] = v[1893] - v[1905] * v[327] + v[1904] * v[328];
//			v[1901] = v[1901] - v[1904] * v[326] + v[1903] * v[327];
//			v[1897] = v[1897] + v[1905] * v[326] - v[1903] * v[328];
//		}
//		else {
//		};
//		v[2319] = -(v[2202] * v[646]) - v[2201] * v[647] - v[2200] * v[648] + v[3936] * v[649] + v[3935] * v[651]
//			+ v[3934] * v[653] - v[2195] * v[655] - v[2194] * v[656] - v[2193] * v[657];
//		v[2320] = -(v[2202] * v[661]) - v[2201] * v[662] - v[2200] * v[663] + v[3936] * v[664] + v[3935] * v[666]
//			+ v[3934] * v[668] - v[2195] * v[670] - v[2194] * v[671] - v[2193] * v[672];
//		v[2321] = v[2202] * v[616] + v[2201] * v[617] + v[2200] * v[618] - v[3936] * v[619] - v[3935] * v[621] - v[3934] * v[623]
//			+ v[2195] * v[625] + v[2194] * v[626] + v[2193] * v[627];
//		v[2322] = v[2202] * v[631] + v[2201] * v[632] + v[2200] * v[633] - v[3936] * v[634] - v[3935] * v[636] - v[3934] * v[638]
//			+ v[2195] * v[640] + v[2194] * v[641] + v[2193] * v[642];
//		v[2323] = -(v[2189] * v[646]) - v[2188] * v[647] - v[2187] * v[648] + v[3939] * v[649] + v[3938] * v[651]
//			+ v[3937] * v[653] - v[2182] * v[655] - v[2181] * v[656] - v[2180] * v[657];
//		v[2324] = -(v[2189] * v[661]) - v[2188] * v[662] - v[2187] * v[663] + v[3939] * v[664] + v[3938] * v[666]
//			+ v[3937] * v[668] - v[2182] * v[670] - v[2181] * v[671] - v[2180] * v[672];
//		v[2325] = v[2189] * v[616] + v[2188] * v[617] + v[2187] * v[618] - v[3939] * v[619] - v[3938] * v[621] - v[3937] * v[623]
//			+ v[2182] * v[625] + v[2181] * v[626] + v[2180] * v[627];
//		v[2326] = v[2189] * v[631] + v[2188] * v[632] + v[2187] * v[633] - v[3939] * v[634] - v[3938] * v[636] - v[3937] * v[638]
//			+ v[2182] * v[640] + v[2181] * v[641] + v[2180] * v[642];
//		v[2327] = -(v[2176] * v[646]) - v[2175] * v[647] - v[2174] * v[648] + v[3942] * v[649] + v[3941] * v[651]
//			+ v[3940] * v[653] - v[2169] * v[655] - v[2168] * v[656] - v[2167] * v[657];
//		v[2328] = -(v[2176] * v[661]) - v[2175] * v[662] - v[2174] * v[663] + v[3942] * v[664] + v[3941] * v[666]
//			+ v[3940] * v[668] - v[2169] * v[670] - v[2168] * v[671] - v[2167] * v[672];
//		v[2329] = v[2176] * v[616] + v[2175] * v[617] + v[2174] * v[618] - v[3942] * v[619] - v[3941] * v[621] - v[3940] * v[623]
//			+ v[2169] * v[625] + v[2168] * v[626] + v[2167] * v[627];
//		v[2330] = v[2176] * v[631] + v[2175] * v[632] + v[2174] * v[633] - v[3942] * v[634] - v[3941] * v[636] - v[3940] * v[638]
//			+ v[2169] * v[640] + v[2168] * v[641] + v[2167] * v[642];
//		v[2337] = v[1241] * v[2105] + v[1242] * v[2106] + v[1243] * v[2107] + (*ct)*(-(v[1154] * v[2134]) - v[1155] * v[2135]
//			- v[1156] * v[2136]) + v[1790] * v[2331] + v[1791] * v[2332] + v[1794] * v[2333] + v[2049] * v[4338] + v[2085] * v[4339]
//			+ v[2041] * v[4340] - v[2661] * v[566] - v[2667] * v[569] - v[2674] * v[572] - v[2125] * v[627] - v[2124] * v[642]
//			+ v[2127] * v[657] + v[2126] * v[672];
//		v[2467] = v[2337] * v[3720];
//		v[2341] = v[1245] * v[2105] + v[1246] * v[2106] + v[1247] * v[2107] + (*ct)*(-(v[1150] * v[2134]) - v[1151] * v[2135]
//			- v[1152] * v[2136]) + v[1771] * v[2331] + v[1782] * v[2332] + v[1797] * v[2333] + v[2049] * v[4341] + v[2085] * v[4342]
//			+ v[2041] * v[4343] - v[2661] * v[565] - v[2667] * v[568] - v[2674] * v[571] - v[2125] * v[626] - v[2124] * v[641]
//			+ v[2127] * v[656] + v[2126] * v[671];
//		v[3964] = v[2341] * v[279];
//		v[2471] = v[2341] * v[3719];
//		v[2345] = v[1249] * v[2105] + v[1250] * v[2106] + v[1251] * v[2107] + (*ct)*(-(v[1146] * v[2134]) - v[1147] * v[2135]
//			- v[1148] * v[2136]) + v[1772] * v[2331] + v[1792] * v[2332] + v[1793] * v[2333] + v[2049] * v[4344] + v[2085] * v[4345]
//			+ v[2041] * v[4346] - v[2661] * v[564] - v[2667] * v[567] - v[2674] * v[570] - v[2125] * v[625] - v[2124] * v[640]
//			+ v[2127] * v[655] + v[2126] * v[670];
//		v[3965] = v[2345] * v[279];
//		v[2474] = v[2345] * v[3718];
//		v[2349] = v[1253] * v[2105] + v[1254] * v[2106] + v[1255] * v[2107] + (*ct)*(-(v[1130] * v[2134]) - v[1131] * v[2135]
//			- v[1132] * v[2136]) - v[1716] * v[2331] - v[1717] * v[2332] - v[1720] * v[2333] + v[2049] * v[4347] + v[2085] * v[4348]
//			+ v[2041] * v[4349] + v[2661] * v[467] + v[2667] * v[470] + v[2674] * v[473] - v[2125] * v[618] - v[2124] * v[633]
//			+ v[2127] * v[648] + v[2126] * v[663];
//		v[2531] = v[2349] * v[3717];
//		v[2353] = v[1257] * v[2105] + v[1258] * v[2106] + v[1259] * v[2107] + (*ct)*(-(v[1126] * v[2134]) - v[1127] * v[2135]
//			- v[1128] * v[2136]) - v[1697] * v[2331] - v[1708] * v[2332] - v[1723] * v[2333] + v[2049] * v[4350] + v[2085] * v[4351]
//			+ v[2041] * v[4352] + v[2661] * v[466] + v[2667] * v[469] + v[2674] * v[472] - v[2125] * v[617] - v[2124] * v[632]
//			+ v[2127] * v[647] + v[2126] * v[662];
//		v[3969] = v[2353] * v[253];
//		v[2535] = v[2353] * v[3716];
//		v[2357] = v[1261] * v[2105] + v[1262] * v[2106] + v[1263] * v[2107] + (*ct)*(-(v[1122] * v[2134]) - v[1123] * v[2135]
//			- v[1124] * v[2136]) - v[1698] * v[2331] - v[1718] * v[2332] - v[1719] * v[2333] + v[2049] * v[4353] + v[2085] * v[4354]
//			+ v[2041] * v[4355] + v[2661] * v[465] + v[2667] * v[468] + v[2674] * v[471] - v[2125] * v[616] - v[2124] * v[631]
//			+ v[2127] * v[646] + v[2126] * v[661];
//		v[6846] = v[4004];
//		v[6847] = v[4005];
//		v[6848] = v[4006];
//		v[6849] = v[1260] * v[1826] + v[1256] * v[1830] + v[1706] * (v[2541] + v[1256] * v[2638] + v[1260] * v[2640])
//			+ v[2357] * v[2682] + v[2353] * v[3724] + v[2349] * v[3727] + v[1647] * v[4007] + v[1683] * v[4008] + v[1689] * v[4009]
//			+ v[155] * v[4010] + v[160] * v[4011] + v[151] * v[4012];
//		v[6850] = v[1264] * v[1826] + v[1256] * v[1832] + v[1571] * v[2349] + v[1572] * v[2357] + v[1706] * (v[2534]
//			+ v[1256] * v[2637] + v[1264] * v[2641]) + v[2353] * v[2683] + v[1680] * v[3811] + v[1700] * v[3813] + v[1645] * v[4013]
//			+ v[152] * v[4014] + v[157] * v[4015] + v[1840] * v[4153];
//		v[6851] = v[1264] * v[1830] + v[1260] * v[1832] + v[1573] * v[2353] + v[1574] * v[2357] + v[1706] * (v[2526]
//			+ v[1260] * v[2639] + v[1264] * v[2642]) + v[2349] * v[2684] + v[1685] * v[3812] + v[1643] * v[4016] + v[1691] * v[4017]
//			+ v[163] * v[4018] + v[1844] * v[4157] + v[1848] * v[4158];
//		v[6852] = -v[4004];
//		v[6853] = -v[4005];
//		v[6854] = -v[4006];
//		v[6855] = v[1248] * v[1854] + v[1244] * v[1858] + v[1780] * (v[2477] + v[1244] * v[2609] + v[1248] * v[2611])
//			+ v[2345] * v[2685] + v[2341] * v[3732] + v[2337] * v[3735] + v[1673] * v[4019] + v[1757] * v[4020] + v[1763] * v[4021]
//			+ v[211] * v[4022] + v[216] * v[4023] + v[207] * v[4024];
//		v[6856] = v[1252] * v[1854] + v[1244] * v[1860] + v[1577] * v[2337] + v[1578] * v[2345] + v[1780] * (v[2470]
//			+ v[1244] * v[2608] + v[1252] * v[2612]) + v[2341] * v[2686] + v[1754] * v[3793] + v[1774] * v[3795] + v[1671] * v[4025]
//			+ v[208] * v[4026] + v[213] * v[4027] + v[1868] * v[4165];
//		v[6857] = v[1252] * v[1858] + v[1248] * v[1860] + v[1579] * v[2341] + v[1580] * v[2345] + v[1780] * (v[2462]
//			+ v[1248] * v[2610] + v[1252] * v[2613]) + v[2337] * v[2687] + v[1759] * v[3794] + v[1669] * v[4028] + v[1765] * v[4029]
//			+ v[219] * v[4030] + v[1872] * v[4169] + v[1876] * v[4170];
//		v[3970] = v[2357] * v[253];
//		v[2538] = v[2357] * v[3715];
//		v[2247] = v[2247] + v[2246] * v[471] + v[2245] * v[472] + v[2244] * v[473] - v[2243] * v[570] - v[2242] * v[571]
//			- v[2241] * v[572];
//		v[2358] = v[2085] * v[2374] + v[2085] * v[2391] - v[2674] * v[301] - v[319] * (v[2048] + v[3948] + v[3960])
//			- v[2241] * v[717];
//		v[2359] = v[2085] * v[2372] + v[2085] * v[2388] - v[2674] * v[291] - v[319] * (v[2047] + v[3947] + v[3958])
//			- v[2242] * v[717];
//		v[2360] = v[2085] * v[2370] + v[2085] * v[2385] - v[2674] * v[281] - v[319] * (v[2046] + v[3946] + v[3956])
//			- v[2243] * v[717];
//		v[2361] = v[2674] * v[275] + v[319] * (v[2045] + v[3945] + v[3954]) + v[2085] * v[4356] + v[2244] * v[717];
//		v[2362] = v[265] * v[2674] + v[319] * (v[2044] + v[3944] + v[3952]) + v[2085] * v[4357] + v[2245] * v[717];
//		v[2271] = v[2271] + v[2260] * v[465] + v[2259] * v[466] + v[2258] * v[467] - v[2257] * v[564] - v[2256] * v[565]
//			- v[2255] * v[566];
//		v[1901] = v[1901] + v[2241] * v[2336] + v[2242] * v[2340] + v[2243] * v[2344] + v[2244] * v[2348] + v[2245] * v[2352]
//			+ v[2246] * v[2356] + 2e0*v[2247] * v[319] + v[2042] * v[3335] + v[2043] * v[471] + v[2044] * v[472] + v[2045] * v[473]
//			- v[2046] * v[570] - v[2047] * v[571] - v[2048] * v[572] + v[318] * (v[2253] * v[471] + v[2252] * v[472] + v[2251] * v[473]
//				- v[2250] * v[570] - v[2249] * v[571] - v[2248] * v[572]) + v[317] * (v[2260] * v[471] + v[2259] * v[472]
//					+ v[2258] * v[473] - v[2257] * v[570] - v[2256] * v[571] - v[2255] * v[572]);
//		v[2254] = v[2254] + v[2253] * v[468] + v[2252] * v[469] + v[2251] * v[470] - v[2250] * v[567] - v[2249] * v[568]
//			- v[2248] * v[569];
//		v[2363] = v[255] * v[2674] + v[319] * (v[2043] + v[3943] + v[3950]) + v[2085] * v[4358] + v[2246] * v[717];
//		v[2365] = v[255] * v[2667] + v[318] * (v[2053] + v[3943] + v[3949]) + v[2049] * v[4359] + v[2253] * v[715];
//		v[2367] = v[265] * v[2667] + v[318] * (v[2055] + v[3944] + v[3951]) + v[2049] * v[4360] + v[2252] * v[715];
//		v[2369] = v[2667] * v[275] + v[318] * (v[2057] + v[3945] + v[3953]) + v[2049] * v[4361] + v[2251] * v[715];
//		v[2371] = v[2049] * v[2370] + v[2049] * v[2386] - v[2667] * v[281] - v[318] * (v[2059] + v[3946] + v[3955])
//			- v[2250] * v[715];
//		v[2373] = v[2049] * v[2372] + v[2049] * v[2389] - v[2667] * v[291] - v[318] * (v[2061] + v[3947] + v[3957])
//			- v[2249] * v[715];
//		v[1897] = v[1897] + 2e0*v[2254] * v[318] + v[2051] * v[3354] + v[2248] * v[3746] + v[2249] * v[3748] + v[2250] * v[3750]
//			+ v[2251] * v[3752] + v[2252] * v[3754] + v[2253] * v[3756] + v[2053] * v[468] + v[2055] * v[469] + v[2057] * v[470]
//			- v[2059] * v[567] - v[2061] * v[568] - v[2063] * v[569] + v[319] * (v[2246] * v[468] + v[2245] * v[469] + v[2244] * v[470]
//				- v[2243] * v[567] - v[2242] * v[568] - v[2241] * v[569]) + v[317] * (v[2260] * v[468] + v[2259] * v[469]
//					+ v[2258] * v[470] - v[2257] * v[567] - v[2256] * v[568] - v[2255] * v[569]);
//		v[2375] = v[2049] * v[2374] + v[2049] * v[2392] - v[2667] * v[301] - v[318] * (v[2063] + v[3948] + v[3959])
//			- v[2248] * v[715];
//		v[2378] = v[255] * v[2661] + v[317] * (v[2069] + v[3949] + v[3950]) + v[2041] * v[4362] + v[2260] * v[713];
//		v[2381] = v[265] * v[2661] + v[317] * (v[2072] + v[3951] + v[3952]) + v[2041] * v[4363] + v[2259] * v[713];
//		v[2384] = v[2661] * v[275] + v[317] * (v[2075] + v[3953] + v[3954]) + v[2041] * v[4364] + v[2258] * v[713];
//		v[2387] = v[2041] * v[2385] + v[2041] * v[2386] - v[2661] * v[281] - v[317] * (v[2078] + v[3955] + v[3956])
//			- v[2257] * v[713];
//		v[2390] = v[2041] * v[2388] + v[2041] * v[2389] - v[2661] * v[291] - v[317] * (v[2081] + v[3957] + v[3958])
//			- v[2256] * v[713];
//		v[1893] = v[1893] + 2e0*v[2271] * v[317] + v[2066] * v[3370] - v[2255] * v[3747] - v[2256] * v[3749] - v[2257] * v[3751]
//			+ v[2258] * v[3753] + v[2259] * v[3755] + v[2260] * v[3757] + v[2069] * v[465] + v[2072] * v[466] + v[2075] * v[467]
//			- v[2078] * v[564] - v[2081] * v[565] - v[2084] * v[566] + v[319] * (v[2246] * v[465] + v[2245] * v[466] + v[2244] * v[467]
//				- v[2243] * v[564] - v[2242] * v[565] - v[2241] * v[566]) + v[318] * (v[2253] * v[465] + v[2252] * v[466]
//					+ v[2251] * v[467] - v[2250] * v[564] - v[2249] * v[565] - v[2248] * v[566]);
//		v[2393] = v[2041] * v[2391] + v[2041] * v[2392] - v[2661] * v[301] - v[317] * (v[2084] + v[3959] + v[3960])
//			- v[2255] * v[713];
//		v[1889] = v[1889] + v[2089] * v[2394] * v[314] + (v[1209] * v[1880] + v[1210] * v[1883] + v[1211] * v[1886]
//			+ v[1893] * v[302] + v[1897] * v[303] + v[1901] * v[304])*v[315];
//		v[2396] = v[1901] * v[316] + v[1211] * v[3854] + (v[1193] * v[1886] + v[1889] * v[304]) / v[727];
//		v[2398] = v[1897] * v[316] + v[1210] * v[3854] + (v[1193] * v[1883] + v[1889] * v[303]) / v[727];
//		v[2400] = v[1893] * v[316] + v[1209] * v[3854] + (v[1193] * v[1880] + v[1889] * v[302]) / v[727];
//		v[2303] = v[2303] - v[2396];
//		v[2307] = v[2307] + v[2396];
//		v[2295] = v[2295] - v[2398];
//		v[2299] = v[2299] + v[2398];
//		v[2287] = v[2287] - v[2400];
//		v[2291] = v[2291] + v[2400];
//		v[2401] = v[2337] * v[288] + v[1244] * v[3961];
//		v[2402] = v[2337] * v[282] + v[1244] * v[3962];
//		v[2403] = v[1378] * v[2337] + v[1244] * (v[1811] / v[279] + v[1780] * v[3401]) + v[2401] * v[4109];
//		v[2404] = v[1384] * v[2337] + v[1244] * (v[1803] / v[279] + v[1780] * v[3403]) + v[2402] * v[4110];
//		v[2405] = (v[218] * v[2401] + v[216] * v[2402] + v[1244] * (v[1819] + v[3405] * v[3852]) + v[2337] * v[4111]) / v[279];
//		v[2406] = v[2341] * v[290] + v[1248] * v[3963];
//		v[2407] = v[2341] * v[282] + v[1248] * v[3962];
//		v[2410] = v[1373] * v[2341] + v[1248] * (v[1817] / v[279] + v[1780] * v[3409]) + v[2406] * v[4113];
//		v[2411] = v[2401] + v[2406];
//		v[2412] = v[2403] + v[2410];
//		v[2414] = (v[1376] * v[1765] + v[207] * v[2407] + v[209] * v[2411] + v[2604] * v[3852] + v[1248] * (v[1799]
//			+ v[3413] * v[3852]) + v[1383] * v[3964]) / v[279];
//		v[2417] = (v[214] * v[2406] + v[211] * v[2407] + v[1248] * (v[1808] + v[3415] * v[3852]) + v[1377] * v[3964]) / v[279];
//		v[2418] = v[2345] * v[288] + v[1252] * v[3961];
//		v[2419] = v[2345] * v[290] + v[1252] * v[3963];
//		v[2420] = v[2407] + v[2418];
//		v[2423] = (v[208] * v[2418] + v[209] * v[2419] + v[1252] * (v[1800] + v[3420] * v[3852]) + v[1381] * v[3965]) / v[279];
//		v[2424] = v[2402] + v[2419];
//		v[2426] = (v[1386] * v[1759] + v[213] * v[2418] + v[214] * v[2424] + v[2603] * v[3852] + v[1252] * (v[1807]
//			+ v[3423] * v[3852]) + v[1388] * v[3965]) / v[279];
//		v[2427] = v[2414] + v[2426];
//		v[2429] = (v[1387] * v[1754] + v[219] * v[2419] + v[218] * v[2420] + v[2602] * v[3852] + v[1252] * (v[1815]
//			+ v[3426] * v[3852]) + v[1391] * v[3965]) / v[279];
//		v[2430] = v[2404] + v[2429];
//		v[2431] = v[2349] * v[262] + v[1256] * v[3966];
//		v[2432] = v[2349] * v[256] + v[1256] * v[3967];
//		v[2433] = v[1402] * v[2349] + v[1256] * (v[1737] / v[253] + v[1706] * v[3431]) + v[2431] * v[4118];
//		v[2434] = v[1408] * v[2349] + v[1256] * (v[1729] / v[253] + v[1706] * v[3433]) + v[2432] * v[4119];
//		v[2435] = (v[162] * v[2431] + v[160] * v[2432] + v[1256] * (v[1745] + v[3435] * v[3851]) + v[2349] * v[4120]) / v[253];
//		v[2436] = v[2353] * v[264] + v[1260] * v[3968];
//		v[2437] = v[2353] * v[256] + v[1260] * v[3967];
//		v[2440] = v[1397] * v[2353] + v[1260] * (v[1743] / v[253] + v[1706] * v[3439]) + v[2436] * v[4122];
//		v[2441] = v[2431] + v[2436];
//		v[2442] = v[2433] + v[2440];
//		v[2444] = (v[1400] * v[1691] + v[151] * v[2437] + v[153] * v[2441] + v[2633] * v[3851] + v[1260] * (v[1725]
//			+ v[3443] * v[3851]) + v[1407] * v[3969]) / v[253];
//		v[2447] = (v[158] * v[2436] + v[155] * v[2437] + v[1260] * (v[1734] + v[3445] * v[3851]) + v[1401] * v[3969]) / v[253];
//		v[2448] = v[2357] * v[262] + v[1264] * v[3966];
//		v[2449] = v[2357] * v[264] + v[1264] * v[3968];
//		v[2450] = v[2437] + v[2448];
//		v[2453] = (v[152] * v[2448] + v[153] * v[2449] + v[1264] * (v[1726] + v[3450] * v[3851]) + v[1405] * v[3970]) / v[253];
//		v[2454] = v[2432] + v[2449];
//		v[2456] = (v[1410] * v[1685] + v[157] * v[2448] + v[158] * v[2454] + v[2632] * v[3851] + v[1264] * (v[1733]
//			+ v[3453] * v[3851]) + v[1412] * v[3970]) / v[253];
//		v[2457] = v[2444] + v[2456];
//		v[2459] = (v[1411] * v[1680] + v[163] * v[2449] + v[162] * v[2450] + v[2631] * v[3851] + v[1264] * (v[1741]
//			+ v[3456] * v[3851]) + v[1415] * v[3970]) / v[253];
//		v[2460] = v[2434] + v[2459];
//		v[2304] = v[2304] - v[200] * v[2396] + v[2327] * v[348] + v[2328] * v[358];
//		v[2305] = v[2305] - v[199] * v[2396] + v[2327] * v[347] + v[2328] * v[356];
//		v[2306] = v[2306] - v[198] * v[2396] + v[2327] * v[346] + v[2328] * v[354];
//		v[2296] = v[2296] - v[200] * v[2398] + v[2323] * v[348] + v[2324] * v[358];
//		v[2297] = v[2297] - v[199] * v[2398] + v[2323] * v[347] + v[2324] * v[356];
//		v[2298] = v[2298] - v[198] * v[2398] + v[2323] * v[346] + v[2324] * v[354];
//		v[2288] = v[2288] - v[200] * v[2400] + v[2319] * v[348] + v[2320] * v[358];
//		v[2289] = v[2289] - v[199] * v[2400] + v[2319] * v[347] + v[2320] * v[356];
//		v[2290] = v[2290] - v[198] * v[2400] + v[2319] * v[346] + v[2320] * v[354];
//		v[2461] = QBi[2][2] * v[2304] + QBi[2][1] * v[2305] + QBi[2][0] * v[2306] + v[1577] * v[2406] + v[2467] * v[292]
//			+ v[2419] * v[3735] + v[290] * (v[4030] + v[1780] * v[4365]);
//		v[2464] = QBi[1][2] * v[2304] + QBi[1][1] * v[2305] + QBi[1][0] * v[2306] + v[2401] * v[2687] + (v[1387] * v[1851])
//			/ v[279] + v[2471] * v[289] + v[1868] * v[3719] + v[2420] * v[3735] + v[1780] * v[4366];
//		v[2466] = QBi[0][2] * v[2304] + QBi[0][1] * v[2305] + QBi[0][0] * v[2306] + v[2402] * v[2687] + v[2474] * v[277]
//			+ v[282] * (v[4023] + v[1780] * v[4367]);
//		v[2468] = QBi[2][2] * v[2296] + QBi[2][1] * v[2297] + QBi[2][0] * v[2298] + v[2406] * v[2686] + (v[1386] * v[1855])
//			/ v[279] + v[1876] * v[3720] + v[2424] * v[3732] + v[2467] * v[3736] + v[1780] * v[4368];
//		v[2472] = QBi[1][2] * v[2296] + QBi[1][1] * v[2297] + QBi[1][0] * v[2298] + v[1579] * v[2401] + v[2471] * v[283]
//			+ v[2418] * v[3732] + v[288] * (v[4027] + v[1780] * v[4369]);
//		v[2475] = QBi[0][2] * v[2296] + QBi[0][1] * v[2297] + QBi[0][0] * v[2298] + v[2407] * v[2686] + v[2474] * v[278]
//			+ v[282] * (v[4022] + v[1780] * v[4370]);
//		v[2476] = QBi[2][2] * v[2288] + QBi[2][1] * v[2289] + QBi[2][0] * v[2290] + v[1578] * v[2411] + v[2419] * v[2685] +
//			(v[1376] * v[1856]) / v[279] + v[1872] * v[3720] + v[2467] * v[3734] + v[1780] * v[4371];
//		v[2478] = QBi[1][2] * v[2288] + QBi[1][1] * v[2289] + QBi[1][0] * v[2290] + v[2418] * v[2685] + v[2471] * v[3731]
//			+ v[288] * (v[1780] * v[3973] + v[4026]);
//		v[2481] = QBi[0][2] * v[2288] + QBi[0][1] * v[2289] + QBi[0][0] * v[2290] + v[1580] * v[2402] + v[1578] * v[2407]
//			+ v[2474] * v[280] + v[282] * (v[4024] + v[1780] * v[4372]);
//		v[2482] = -(v[2390] * v[855]);
//		v[2483] = -(v[2390] * v[853]);
//		v[2484] = -(v[2390] * v[851]);
//		v[2485] = -(v[2393] * v[855]);
//		v[2486] = -(v[2393] * v[851]);
//		v[2487] = -(v[2393] * v[853]);
//		v[2488] = -(v[2387] * v[853]);
//		v[2489] = -(v[2387] * v[851]);
//		v[2490] = -(v[2387] * v[855]);
//		v[2491] = -(v[2371] * v[855]);
//		v[2492] = -(v[2371] * v[853]);
//		v[2493] = -(v[2371] * v[851]);
//		v[2494] = -(v[2375] * v[851]);
//		v[2495] = -(v[2375] * v[855]);
//		v[2496] = -(v[2375] * v[853]);
//		v[2497] = -(v[2358] * v[853]);
//		v[2498] = -(v[2358] * v[855]);
//		v[2499] = -(v[2358] * v[851]);
//		v[2500] = -(v[1380] * v[1665]) - v[1390] * v[1667] + 2e0*v[1379] * v[1668] + v[2488] + v[2491] + v[2494] + v[2497]
//			+ 2e0*v[2417] * v[481] - v[2427] * v[484] - v[2412] * v[487];
//		v[2501] = -(v[2373] * v[855]);
//		v[2502] = -(v[2373] * v[851]);
//		v[2503] = -(v[2373] * v[853]);
//		v[2504] = v[1393] * v[1665] + 2e0*v[1385] * v[1667] - v[1390] * v[1668] - v[2483] - v[2486] - v[2498] - v[2501]
//			- v[2427] * v[481] + 2e0*v[2423] * v[484] + v[2430] * v[487];
//		v[2505] = -(v[1441] * v[1656]) + v[1437] * v[1659] - v[1442] * v[1661] - 2e0*v[1432] * v[1666] + v[204] * v[2478]
//			- v[2488] * v[474] + v[2483] * v[475] - v[2487] * v[476];
//		v[2506] = -(v[2360] * v[855]);
//		v[2507] = -(v[2360] * v[853]);
//		v[2508] = -(v[2360] * v[851]);
//		v[2509] = -(v[2359] * v[853]);
//		v[2510] = -(v[2359] * v[855]);
//		v[2511] = -(v[2359] * v[851]);
//		v[2512] = 2e0*v[1372] * v[1665] + v[1393] * v[1667] - v[1380] * v[1668] - v[2489] - v[2502] - v[2506] - v[2509]
//			- v[2412] * v[481] + v[2430] * v[484] + 2e0*v[2405] * v[487];
//		v[2513] = -(v[1440] * v[1656]) + v[1438] * v[1659] - v[1443] * v[1661] - 2e0*v[1430] * v[1666] + v[204] * v[2476]
//			- v[2489] * v[474] + v[2484] * v[475] - v[2486] * v[476];
//		v[2514] = -(v[1449] * v[1656]) + v[1455] * v[1659] - v[1445] * v[1661] - 2e0*v[1429] * v[1666] + v[204] * v[2475]
//			- v[2491] * v[474] + v[2501] * v[475] - v[2495] * v[476];
//		v[2515] = v[2485] + v[2496];
//		v[2516] = -(v[1448] * v[1656]) + v[1456] * v[1659] - v[1447] * v[1661] - 2e0*v[1424] * v[1666] + v[204] * v[2468]
//			- v[2493] * v[474] + v[2502] * v[475] - v[2494] * v[476];
//		v[2517] = -(v[1452] * v[1656]) + v[1464] * v[1659] - v[1460] * v[1661] - 2e0*v[1423] * v[1666] + v[204] * v[2466]
//			- v[2506] * v[474] + v[2510] * v[475] - v[2498] * v[476];
//		v[2518] = -(v[1451] * v[1656]) + v[1463] * v[1659] - v[1461] * v[1661] - 2e0*v[1421] * v[1666] + v[204] * v[2464]
//			- v[2507] * v[474] + v[2509] * v[475] - v[2497] * v[476];
//		v[2519] = v[2492] + v[2508];
//		v[2520] = v[2482] + v[2511];
//		v[2521] = -(QBi[2][2] * v[1795]) - QBi[1][2] * v[1796] - QBi[0][2] * v[1798] - v[1369] * v[1806] - v[1368] * v[1814]
//			- v[1367] * v[1822] - v[231] * v[2396] - v[228] * v[2398] - v[225] * v[2400] + v[2387] * v[525] + v[2390] * v[526]
//			+ v[2393] * v[527] + v[2371] * v[540] + v[2373] * v[541] + v[2375] * v[542] + v[2360] * v[555] + v[2359] * v[556]
//			+ v[2358] * v[557];
//		v[2522] = -(QBi[2][1] * v[1795]) - QBi[1][1] * v[1796] - QBi[0][1] * v[1798] - v[1369] * v[1805] - v[1368] * v[1813]
//			- v[1367] * v[1821] - v[230] * v[2396] - v[227] * v[2398] - v[224] * v[2400] + v[2387] * v[522] + v[2390] * v[523]
//			+ v[2393] * v[524] + v[2371] * v[537] + v[2373] * v[538] + v[2375] * v[539] + v[2360] * v[552] + v[2359] * v[553]
//			+ v[2358] * v[554];
//		v[2523] = -(QBi[2][0] * v[1795]) - QBi[1][0] * v[1796] - QBi[0][0] * v[1798] - v[1369] * v[1804] - v[1368] * v[1812]
//			- v[1367] * v[1820] - v[229] * v[2396] - v[226] * v[2398] - v[223] * v[2400] + v[2387] * v[519] + v[2390] * v[520]
//			+ v[2393] * v[521] + v[2371] * v[534] + v[2373] * v[535] + v[2375] * v[536] + v[2360] * v[549] + v[2359] * v[550]
//			+ v[2358] * v[551];
//		v[2524] = v[2523] * x1B[0] + v[2522] * x1B[1] + v[2521] * x1B[2];
//		v[2308] = v[2308] + v[144] * v[2396] + v[2329] * v[331] + v[2330] * v[341];
//		v[2309] = v[2309] + v[143] * v[2396] + v[2329] * v[330] + v[2330] * v[339];
//		v[2310] = v[2310] + v[142] * v[2396] + v[2329] * v[329] + v[2330] * v[337];
//		v[2300] = v[2300] + v[144] * v[2398] + v[2325] * v[331] + v[2326] * v[341];
//		v[2301] = v[2301] + v[143] * v[2398] + v[2325] * v[330] + v[2326] * v[339];
//		v[2302] = v[2302] + v[142] * v[2398] + v[2325] * v[329] + v[2326] * v[337];
//		v[2292] = v[2292] + v[144] * v[2400] + v[2321] * v[331] + v[2322] * v[341];
//		v[2293] = v[2293] + v[143] * v[2400] + v[2321] * v[330] + v[2322] * v[339];
//		v[2294] = v[2294] + v[142] * v[2400] + v[2321] * v[329] + v[2322] * v[337];
//		v[2525] = QAi[2][2] * v[2308] + QAi[2][1] * v[2309] + QAi[2][0] * v[2310] + v[1571] * v[2436] + v[2531] * v[266]
//			+ v[2449] * v[3727] + v[264] * (v[4018] + v[1706] * v[4373]);
//		v[2528] = QAi[1][2] * v[2308] + QAi[1][1] * v[2309] + QAi[1][0] * v[2310] + (v[1411] * v[1823]) / v[253]
//			+ v[2535] * v[263] + v[2431] * v[2684] + v[1840] * v[3716] + v[2450] * v[3727] + v[1706] * v[4374];
//		v[2530] = QAi[0][2] * v[2308] + QAi[0][1] * v[2309] + QAi[0][0] * v[2310] + v[251] * v[2538] + v[2432] * v[2684]
//			+ v[256] * (v[4011] + v[1706] * v[4375]);
//		v[2532] = QAi[2][2] * v[2300] + QAi[2][1] * v[2301] + QAi[2][0] * v[2302] + (v[1410] * v[1827]) / v[253]
//			+ v[2436] * v[2683] + v[1848] * v[3717] + v[2454] * v[3724] + v[2531] * v[3728] + v[1706] * v[4376];
//		v[2536] = QAi[1][2] * v[2300] + QAi[1][1] * v[2301] + QAi[1][0] * v[2302] + v[1573] * v[2431] + v[2535] * v[257]
//			+ v[2448] * v[3724] + v[262] * (v[4015] + v[1706] * v[4377]);
//		v[2539] = QAi[0][2] * v[2300] + QAi[0][1] * v[2301] + QAi[0][0] * v[2302] + v[252] * v[2538] + v[2437] * v[2683]
//			+ v[256] * (v[4010] + v[1706] * v[4378]);
//		v[2540] = QAi[2][2] * v[2292] + QAi[2][1] * v[2293] + QAi[2][0] * v[2294] + v[1572] * v[2441] + (v[1400] * v[1828])
//			/ v[253] + v[2449] * v[2682] + v[1844] * v[3717] + v[2531] * v[3726] + v[1706] * v[4379];
//		v[2542] = QAi[1][2] * v[2292] + QAi[1][1] * v[2293] + QAi[1][0] * v[2294] + v[2448] * v[2682] + v[2535] * v[3723]
//			+ v[262] * (v[1706] * v[3976] + v[4014]);
//		v[2545] = QAi[0][2] * v[2292] + QAi[0][1] * v[2293] + QAi[0][0] * v[2294] + v[1574] * v[2432] + v[1572] * v[2437]
//			+ v[2538] * v[254] + v[256] * (v[4012] + v[1706] * v[4380]);
//		v[2546] = v[2381] * v[861];
//		v[2547] = v[2381] * v[859];
//		v[2548] = v[2381] * v[857];
//		v[2549] = v[2384] * v[861];
//		v[2550] = v[2384] * v[857];
//		v[2551] = v[2384] * v[859];
//		v[2552] = v[2378] * v[859];
//		v[2553] = v[2378] * v[857];
//		v[2554] = v[2378] * v[861];
//		v[2555] = v[2365] * v[861];
//		v[2556] = v[2365] * v[859];
//		v[2557] = v[2365] * v[857];
//		v[2558] = v[2369] * v[857];
//		v[2559] = v[2369] * v[861];
//		v[2560] = v[2369] * v[859];
//		v[2561] = v[2361] * v[859];
//		v[2562] = v[2361] * v[861];
//		v[2563] = v[2361] * v[857];
//		v[2564] = -(v[1404] * v[1639]) - v[1414] * v[1641] + 2e0*v[1403] * v[1642] + v[2552] + v[2555] + v[2558] + v[2561]
//			+ 2e0*v[2447] * v[382] - v[2457] * v[385] - v[2442] * v[388];
//		v[2565] = v[2367] * v[861];
//		v[2566] = v[2367] * v[857];
//		v[2567] = v[2367] * v[859];
//		v[2568] = v[1417] * v[1639] + 2e0*v[1409] * v[1641] - v[1414] * v[1642] - v[2547] - v[2550] - v[2562] - v[2565]
//			- v[2457] * v[382] + 2e0*v[2453] * v[385] + v[2460] * v[388];
//		v[2569] = -(v[1504] * v[1630]) + v[1500] * v[1633] - v[1505] * v[1635] - 2e0*v[1495] * v[1640] + v[148] * v[2542]
//			- v[2552] * v[375] + v[2547] * v[376] - v[2551] * v[377];
//		v[2570] = v[2363] * v[861];
//		v[2571] = v[2363] * v[859];
//		v[2572] = v[2363] * v[857];
//		v[2573] = v[2362] * v[859];
//		v[2574] = v[2362] * v[861];
//		v[2575] = v[2362] * v[857];
//		v[2576] = 2e0*v[1396] * v[1639] + v[1417] * v[1641] - v[1404] * v[1642] - v[2553] - v[2566] - v[2570] - v[2573]
//			- v[2442] * v[382] + v[2460] * v[385] + 2e0*v[2435] * v[388];
//		v[2577] = -(v[1503] * v[1630]) + v[1501] * v[1633] - v[1506] * v[1635] - 2e0*v[1493] * v[1640] + v[148] * v[2540]
//			- v[2553] * v[375] + v[2548] * v[376] - v[2550] * v[377];
//		v[2578] = -(v[1512] * v[1630]) + v[1518] * v[1633] - v[1508] * v[1635] - 2e0*v[1492] * v[1640] + v[148] * v[2539]
//			- v[2555] * v[375] + v[2565] * v[376] - v[2559] * v[377];
//		v[2579] = v[2549] + v[2560];
//		v[2580] = -(v[1511] * v[1630]) + v[1519] * v[1633] - v[1510] * v[1635] - 2e0*v[1487] * v[1640] + v[148] * v[2532]
//			- v[2557] * v[375] + v[2566] * v[376] - v[2558] * v[377];
//		v[2581] = -(v[1515] * v[1630]) + v[1527] * v[1633] - v[1523] * v[1635] - 2e0*v[1486] * v[1640] + v[148] * v[2530]
//			- v[2570] * v[375] + v[2574] * v[376] - v[2562] * v[377];
//		v[2582] = -(v[1514] * v[1630]) + v[1526] * v[1633] - v[1524] * v[1635] - 2e0*v[1484] * v[1640] + v[148] * v[2528]
//			- v[2571] * v[375] + v[2573] * v[376] - v[2561] * v[377];
//		v[2583] = v[2556] + v[2572];
//		v[2584] = v[2546] + v[2575];
//		v[2585] = QAi[2][2] * v[1721] + QAi[1][2] * v[1722] + QAi[0][2] * v[1724] + v[1369] * v[1732] + v[1368] * v[1740]
//			+ v[1367] * v[1748] + v[175] * v[2396] + v[172] * v[2398] + v[169] * v[2400] + v[2378] * v[426] + v[2381] * v[427]
//			+ v[2384] * v[428] + v[2365] * v[441] + v[2367] * v[442] + v[2369] * v[443] + v[2363] * v[456] + v[2362] * v[457]
//			+ v[2361] * v[458];
//		v[2586] = QAi[2][1] * v[1721] + QAi[1][1] * v[1722] + QAi[0][1] * v[1724] + v[1369] * v[1731] + v[1368] * v[1739]
//			+ v[1367] * v[1747] + v[174] * v[2396] + v[171] * v[2398] + v[168] * v[2400] + v[2378] * v[423] + v[2381] * v[424]
//			+ v[2384] * v[425] + v[2365] * v[438] + v[2367] * v[439] + v[2369] * v[440] + v[2363] * v[453] + v[2362] * v[454]
//			+ v[2361] * v[455];
//		v[2587] = QAi[2][0] * v[1721] + QAi[1][0] * v[1722] + QAi[0][0] * v[1724] + v[1369] * v[1730] + v[1368] * v[1738]
//			+ v[1367] * v[1746] + v[173] * v[2396] + v[170] * v[2398] + v[167] * v[2400] + v[2378] * v[420] + v[2381] * v[421]
//			+ v[2384] * v[422] + v[2365] * v[435] + v[2367] * v[436] + v[2369] * v[437] + v[2363] * v[450] + v[2362] * v[451]
//			+ v[2361] * v[452];
//		v[2588] = v[2587] * x1A[0] + v[2586] * x1A[1] + v[2585] * x1A[2];
//		v[2589] = (-v[2524] + v[2523] * x3B[0] + v[2522] * x3B[1] + v[2521] * x3B[2]) / 2e0;
//		v[2590] = (-v[2524] + v[2523] * x2B[0] + v[2522] * x2B[1] + v[2521] * x2B[2]) / 2e0;
//		v[2591] = (-(v[1462] * v[1601]) - v[1444] * v[1604] - 2e0*v[1445] * v[1649] - 2e0*v[1442] * v[1650]
//			- 2e0*v[1461] * v[1651] - 2e0*v[1447] * v[1652] - 2e0*v[1460] * v[1653] - 2e0*v[1443] * v[1654] - v[1446] * v[1657]
//			- 2e0*v[2403] + 2e0*v[2410] - v[2490] * v[477] - 2e0*v[2488] * v[482] - 2e0*v[2489] * v[488] - 2e0*v[2491] * v[492]
//			- v[2492] * v[496] - 2e0*v[2493] * v[501] - 2e0*v[2506] * v[505] - 2e0*v[2507] * v[509] - v[2508] * v[514]) / 2e0;
//		v[2592] = (-(v[1439] * v[1656]) + v[1436] * v[1659] - v[1444] * v[1661] - 2e0*v[1435] * v[1666] + v[204] * v[2481]
//			- v[2490] * v[474] + v[2482] * v[475] - v[2485] * v[476]) / 2e0;
//		v[2594] = (v[1465] * v[1601] + v[1436] * v[1604] + 2e0*v[1455] * v[1649] + 2e0*v[1437] * v[1650]
//			+ 2e0*v[1463] * v[1651] + 2e0*v[1456] * v[1652] + 2e0*v[1464] * v[1653] + 2e0*v[1438] * v[1654] + v[1457] * v[1657]
//			- 2e0*v[2404] + 2e0*v[2429] + v[2482] * v[477] + 2e0*v[2483] * v[482] + 2e0*v[2484] * v[488] + 2e0*v[2501] * v[492]
//			+ v[2503] * v[496] + 2e0*v[2502] * v[501] + 2e0*v[2510] * v[505] + 2e0*v[2509] * v[509] + v[2511] * v[514]) / 2e0;
//		v[2595] = (-(v[1450] * v[1656]) + v[1457] * v[1659] - v[1446] * v[1661] - 2e0*v[1427] * v[1666] + v[204] * v[2472]
//			- v[2492] * v[474] + v[2503] * v[475] - v[2496] * v[476]) / 2e0;
//		v[2596] = (-(v[1453] * v[1601]) - v[1439] * v[1604] - 2e0*v[1449] * v[1649] - 2e0*v[1441] * v[1650]
//			- 2e0*v[1451] * v[1651] - 2e0*v[1448] * v[1652] - 2e0*v[1452] * v[1653] - 2e0*v[1440] * v[1654] - v[1450] * v[1657]
//			- 2e0*v[2414] + 2e0*v[2426] - v[2485] * v[477] - 2e0*v[2487] * v[482] - 2e0*v[2486] * v[488] - 2e0*v[2495] * v[492]
//			- v[2496] * v[496] - 2e0*v[2494] * v[501] - 2e0*v[2498] * v[505] - 2e0*v[2497] * v[509] - v[2499] * v[514]) / 2e0;
//		v[2657] = v[1658] * v[2594] - v[2595] + v[1655] * v[2596] + 24e0*v[1603] * v[3650] * v[4382] + v[3998] * (-
//			(v[1418] * v[1601]) - v[1435] * v[1604] - 2e0*v[1429] * v[1649] - 2e0*v[1432] * v[1650] - 2e0*v[1421] * v[1651]
//			- 2e0*v[1424] * v[1652] - 2e0*v[1423] * v[1653] - 2e0*v[1430] * v[1654] - v[1427] * v[1657] - 2e0*v[2484]
//			+ 2e0*v[2487] + 2e0*v[2493] - 2e0*v[2495] - alphaB[1] * v[2500] + alphaB[0] * v[2504] - 2e0*v[2507] + 2e0*v[2510]
//			+ alphaB[2] * v[2512] + 8e0*v[1666] * v[2598] + v[2519] * v[369] + v[2520] * v[372] + v[2515] * v[374] - 4e0*v[204] *
//			(v[2405] + v[1868] * v[2408] + v[1865] * v[2409] + v[1859] * v[2411] + v[1874] * v[2413] + v[1878] * v[2415]
//				+ v[1876] * v[2416] + v[2417] + v[1853] * v[2420] + v[1861] * v[2421] + v[1872] * v[2422] + v[2423] + v[1857] * v[2424]
//				+ v[1864] * v[2425] + v[1870] * v[2428] + v[1819] * v[2462] + v[1817] * v[2463] + v[1815] * v[2465] + v[1811] * v[2469]
//				+ v[1808] * v[2470] + v[1807] * v[2473] + v[1800] * v[2477] + v[1803] * v[2479] + v[1799] * v[2480] + v[1754] * v[2599]
//				+ v[1759] * v[2600] + v[1765] * v[2601] + v[1851] * v[2602] + v[1855] * v[2603] + v[1856] * v[2604] + v[2401] * v[2608]
//				+ v[2402] * v[2609] + v[2406] * v[2610] + v[2407] * v[2611] + v[2418] * v[2612] + v[2419] * v[2613] + v[2337] * v[3977]
//				+ v[2341] * v[3978] + v[2345] * v[3979] + 2e0*v[1780] * v[3652] * v[4384]) - v[2481] * v[477] - 2e0*v[2478] * v[482]
//			- 2e0*v[2476] * v[488] - 2e0*v[2475] * v[492] - v[2472] * v[496] - 2e0*v[2468] * v[501] - 2e0*v[2466] * v[505]
//			- 2e0*v[2464] * v[509] - v[2461] * v[514] - 2e0*v[6713 + i1159]) - 2e0*v[1564] * (-4e0*v[1479] * v[1603]
//				+ 4e0*v[2591] * v[369] + v[6737 + i1159]);
//		v[4003] = v[2657] + (v[1453] * v[1656] - v[1465] * v[1659] + v[1462] * v[1661] + 2e0*v[1418] * v[1666]
//			- v[204] * v[2461] + v[2508] * v[474] - v[2511] * v[475] + v[2499] * v[476]) / 2e0;
//		v[4002] = (v[2513] + v[2517]) / 2e0;
//		v[2616] = v[2516] + v[2518];
//		v[2617] = v[2505] + v[2514];
//		v[2618] = (-v[2588] + v[2587] * x2A[0] + v[2586] * x2A[1] + v[2585] * x2A[2]) / 2e0;
//		v[2619] = (-v[2588] + v[2587] * x3A[0] + v[2586] * x3A[1] + v[2585] * x3A[2]) / 2e0;
//		v[2620] = (-(v[1525] * v[1591]) - v[1507] * v[1594] - 2e0*v[1508] * v[1623] - 2e0*v[1505] * v[1624]
//			- 2e0*v[1524] * v[1625] - 2e0*v[1510] * v[1626] - 2e0*v[1523] * v[1627] - 2e0*v[1506] * v[1628] - v[1509] * v[1631]
//			- 2e0*v[2433] + 2e0*v[2440] - v[2554] * v[378] - 2e0*v[2552] * v[383] - 2e0*v[2553] * v[389] - 2e0*v[2555] * v[393]
//			- v[2556] * v[397] - 2e0*v[2557] * v[402] - 2e0*v[2570] * v[406] - 2e0*v[2571] * v[410] - v[2572] * v[415]) / 2e0;
//		v[2621] = (-(v[1502] * v[1630]) + v[1499] * v[1633] - v[1507] * v[1635] - 2e0*v[1498] * v[1640] + v[148] * v[2545]
//			- v[2554] * v[375] + v[2546] * v[376] - v[2549] * v[377]) / 2e0;
//		v[2623] = (v[1528] * v[1591] + v[1499] * v[1594] + 2e0*v[1518] * v[1623] + 2e0*v[1500] * v[1624]
//			+ 2e0*v[1526] * v[1625] + 2e0*v[1519] * v[1626] + 2e0*v[1527] * v[1627] + 2e0*v[1501] * v[1628] + v[1520] * v[1631]
//			- 2e0*v[2434] + 2e0*v[2459] + v[2546] * v[378] + 2e0*v[2547] * v[383] + 2e0*v[2548] * v[389] + 2e0*v[2565] * v[393]
//			+ v[2567] * v[397] + 2e0*v[2566] * v[402] + 2e0*v[2574] * v[406] + 2e0*v[2573] * v[410] + v[2575] * v[415]) / 2e0;
//		v[2624] = (-(v[1513] * v[1630]) + v[1520] * v[1633] - v[1509] * v[1635] - 2e0*v[1490] * v[1640] + v[148] * v[2536]
//			- v[2556] * v[375] + v[2567] * v[376] - v[2560] * v[377]) / 2e0;
//		v[2625] = (-(v[1516] * v[1591]) - v[1502] * v[1594] - 2e0*v[1512] * v[1623] - 2e0*v[1504] * v[1624]
//			- 2e0*v[1514] * v[1625] - 2e0*v[1511] * v[1626] - 2e0*v[1515] * v[1627] - 2e0*v[1503] * v[1628] - v[1513] * v[1631]
//			- 2e0*v[2444] + 2e0*v[2456] - v[2549] * v[378] - 2e0*v[2551] * v[383] - 2e0*v[2550] * v[389] - 2e0*v[2559] * v[393]
//			- v[2560] * v[397] - 2e0*v[2558] * v[402] - 2e0*v[2562] * v[406] - 2e0*v[2561] * v[410] - v[2563] * v[415]) / 2e0;
//		v[2650] = v[1632] * v[2623] - v[2624] + v[1629] * v[2625] + 24e0*v[1593] * v[3675] * v[4386] + v[3992] * (-
//			(v[1481] * v[1591]) - v[1498] * v[1594] - 2e0*v[1492] * v[1623] - 2e0*v[1495] * v[1624] - 2e0*v[1484] * v[1625]
//			- 2e0*v[1487] * v[1626] - 2e0*v[1486] * v[1627] - 2e0*v[1493] * v[1628] - v[1490] * v[1631] - 2e0*v[2548]
//			+ 2e0*v[2551] + 2e0*v[2557] - 2e0*v[2559] - alphaA[1] * v[2564] + alphaA[0] * v[2568] - 2e0*v[2571] + 2e0*v[2574]
//			+ alphaA[2] * v[2576] + 8e0*v[1640] * v[2627] + v[2583] * v[363] + v[2584] * v[366] + v[2579] * v[368] - v[2545] * v[378]
//			- 2e0*v[2542] * v[383] - 2e0*v[2540] * v[389] - 2e0*v[2539] * v[393] - v[2536] * v[397] - 2e0*v[2532] * v[402]
//			- 2e0*v[2530] * v[406] - 2e0*v[2528] * v[410] - v[2525] * v[415] - 4e0*v[148] * (v[2435] + v[1840] * v[2438]
//				+ v[1837] * v[2439] + v[1831] * v[2441] + v[1846] * v[2443] + v[1850] * v[2445] + v[1848] * v[2446] + v[2447]
//				+ v[1825] * v[2450] + v[1833] * v[2451] + v[1844] * v[2452] + v[2453] + v[1829] * v[2454] + v[1836] * v[2455]
//				+ v[1842] * v[2458] + v[1745] * v[2526] + v[1743] * v[2527] + v[1741] * v[2529] + v[1737] * v[2533] + v[1734] * v[2534]
//				+ v[1733] * v[2537] + v[1726] * v[2541] + v[1729] * v[2543] + v[1725] * v[2544] + v[1680] * v[2628] + v[1685] * v[2629]
//				+ v[1691] * v[2630] + v[1823] * v[2631] + v[1827] * v[2632] + v[1828] * v[2633] + v[2431] * v[2637] + v[2432] * v[2638]
//				+ v[2436] * v[2639] + v[2437] * v[2640] + v[2448] * v[2641] + v[2449] * v[2642] + v[2349] * v[3980] + v[2353] * v[3981]
//				+ v[2357] * v[3982] + 2e0*v[1706] * v[3677] * v[4388]) - 2e0*v[6725 + i1159]) - 2e0*v[1562] * (-4e0*v[1542] * v[1593]
//					+ 4e0*v[2620] * v[363] + v[6785 + i1159]);
//		v[3997] = v[2650] + (v[1516] * v[1630] - v[1528] * v[1633] + v[1525] * v[1635] + 2e0*v[1481] * v[1640]
//			- v[148] * v[2525] + v[2572] * v[375] - v[2575] * v[376] + v[2563] * v[377]) / 2e0;
//		v[3996] = (v[2577] + v[2581]) / 2e0;
//		v[2645] = v[2580] + v[2582];
//		v[2646] = v[2569] + v[2578];
//		v[6822] = 0e0;
//		v[6823] = 0e0;
//		v[6824] = 0e0;
//		v[6825] = 2e0*v[2647];
//		v[6826] = v[2648];
//		v[6827] = v[2649];
//		v[6828] = 0e0;
//		v[6829] = 0e0;
//		v[6830] = 0e0;
//		v[6831] = 0e0;
//		v[6832] = 0e0;
//		v[6833] = 0e0;
//		v[6810] = 0e0;
//		v[6811] = 0e0;
//		v[6812] = 0e0;
//		v[6813] = v[2648];
//		v[6814] = 2e0*v[2651];
//		v[6815] = v[2652];
//		v[6816] = 0e0;
//		v[6817] = 0e0;
//		v[6818] = 0e0;
//		v[6819] = 0e0;
//		v[6820] = 0e0;
//		v[6821] = 0e0;
//		v[6798] = 0e0;
//		v[6799] = 0e0;
//		v[6800] = 0e0;
//		v[6801] = v[2649];
//		v[6802] = v[2652];
//		v[6803] = 2e0*v[2653];
//		v[6804] = 0e0;
//		v[6805] = 0e0;
//		v[6806] = 0e0;
//		v[6807] = 0e0;
//		v[6808] = 0e0;
//		v[6809] = 0e0;
//		v[6750] = 0e0;
//		v[6751] = 0e0;
//		v[6752] = 0e0;
//		v[6753] = 0e0;
//		v[6754] = 0e0;
//		v[6755] = 0e0;
//		v[6756] = 0e0;
//		v[6757] = 0e0;
//		v[6758] = 0e0;
//		v[6759] = 2e0*v[2654];
//		v[6760] = v[2655];
//		v[6761] = v[2656];
//		v[6762] = 0e0;
//		v[6763] = 0e0;
//		v[6764] = 0e0;
//		v[6765] = 0e0;
//		v[6766] = 0e0;
//		v[6767] = 0e0;
//		v[6768] = 0e0;
//		v[6769] = 0e0;
//		v[6770] = 0e0;
//		v[6771] = v[2655];
//		v[6772] = 2e0*v[2658];
//		v[6773] = v[2659];
//		v[6774] = 0e0;
//		v[6775] = 0e0;
//		v[6776] = 0e0;
//		v[6777] = 0e0;
//		v[6778] = 0e0;
//		v[6779] = 0e0;
//		v[6780] = 0e0;
//		v[6781] = 0e0;
//		v[6782] = 0e0;
//		v[6783] = v[2656];
//		v[6784] = v[2659];
//		v[6785] = 2e0*v[2660];
//		v[6834] = v[2291];
//		v[6835] = v[2299];
//		v[6836] = v[2307];
//		v[6837] = v[1521] * v[1640] - v[2580] + v[2582] + v[2646] * v[364] + v[2568] * v[3983] + 2e0*(v[2583] * v[3983]
//			+ v[2620] * v[3985] + v[1593] * (v[1539] * v[3992] - v[1562] * v[3995])) + alphaA[2] * v[3996]
//			+ 2e0*alphaA[0] * v[3997] + v[6821 + i1159];
//		v[6838] = -(v[1517] * v[1640]) + v[2577] - v[2581] + (alphaA[2] * v[2645]) / 2e0 + (alphaA[0] * v[2646]) / 2e0
//			- v[2564] * v[3983] + 2e0*alphaA[1] * (-v[2621] + v[2624] + v[3997]) + v[6809 + i1159] + 2e0*(v[2584] * v[3983]
//				+ v[1593] * (v[1540] * v[3992] + v[1562] * v[3994]) - 4e0*v[2623] * v[962]);
//		v[6839] = v[1532] * v[1640] - v[2569] + v[2578] + 2e0*alphaA[2] * (-v[2621] + v[2650]) + v[2645] * v[364]
//			+ v[2576] * v[3983] + 2e0*(v[2579] * v[3983] + v[2625] * v[3985] + v[1593] * (v[1535] * v[3992] - v[1562] * v[3993]))
//			+ alphaA[0] * v[3996] + v[6797 + i1159];
//		v[6840] = v[2287];
//		v[6841] = v[2295];
//		v[6842] = v[2303];
//		v[6843] = v[1458] * v[1666] - v[2516] + v[2518] + v[2617] * v[370] + v[2504] * v[3986] + 2e0*(v[2519] * v[3986]
//			+ v[2591] * v[3988] + v[1603] * (v[1476] * v[3998] - v[1564] * v[4001])) + alphaB[2] * v[4002]
//			+ 2e0*alphaB[0] * v[4003] + v[6749 + i1159];
//		v[6844] = -(v[1454] * v[1666]) + v[2513] - v[2517] + (alphaB[2] * v[2616]) / 2e0 + (alphaB[0] * v[2617]) / 2e0
//			- v[2500] * v[3986] + 2e0*alphaB[1] * (-v[2592] + v[2595] + v[4003]) + v[6761 + i1159] + 2e0*(v[2520] * v[3986]
//				+ v[1603] * (v[1477] * v[3998] + v[1564] * v[4000]) - 4e0*v[2594] * v[968]);
//		v[6845] = v[1469] * v[1666] - v[2505] + v[2514] + 2e0*alphaB[2] * (-v[2592] + v[2657]) + v[2616] * v[370]
//			+ v[2512] * v[3986] + 2e0*(v[2515] * v[3986] + v[2596] * v[3988] + v[1603] * (v[1472] * v[3998] - v[1564] * v[3999]))
//			+ alphaB[0] * v[4002] + v[6773 + i1159];
//		Rc[i1159 - 1] += v[5857 + i1159] + (*a4)*v[5869 + i1159];
//		for (i1583 = 1; i1583 <= 12; i1583++) {
//			Kc[i1159 - 1][i1583 - 1] += v[2619] * v[5357 + i1583] + v[2618] * v[5369 + i1583] + v[2590] * v[5381 + i1583]
//				+ v[2589] * v[5393 + i1583] + v[6833 + i1583] + (*a4)*v[6845 + i1583];
//		};/* end for */
//	};/* end for */
//	v[2696] = 0e0;
//	v[2697] = 0e0;
//	v[2698] = 0e0;
//	v[2699] = 0e0;
//	b2700 = b728;
//	if (b2700) {
//		b2701 = b747;
//		b2705 = b730;
//		if (b2705) {
//			v[2698] = 0e0;
//			v[2697] = 0e0;
//			v[2696] = 0e0;
//			v[2699] = 0e0;
//		}
//		else {
//		};
//	}
//	else {
//	};
//	v[3705] = 2e0*v[2696] * v[317];
//	v[3710] = v[3705] / 2e0;
//	v[3355] = -(v[301] * v[3705]) / 2e0;
//	v[3349] = -(v[291] * v[3705]) / 2e0;
//	v[3346] = -(v[281] * v[3705]) / 2e0;
//	v[3343] = (v[275] * v[3705]) / 2e0;
//	v[3340] = (v[265] * v[3705]) / 2e0;
//	v[3337] = (v[255] * v[3705]) / 2e0;
//	v[3711] = v[2697] * v[318];
//	v[4087] = -v[3710] - v[3711];
//	v[3707] = 2e0*v[3711];
//	v[3356] = -(v[301] * v[3711]);
//	v[3350] = -(v[291] * v[3711]);
//	v[3347] = -(v[281] * v[3711]);
//	v[3344] = v[275] * v[3711];
//	v[3341] = v[265] * v[3711];
//	v[3338] = v[255] * v[3711];
//	v[4440] = v[2698] * (v[2356] + 2e0*v[3740]) - v[4087] * v[471];
//	v[4437] = v[2698] * (v[2352] + 2e0*v[3741]) - v[4087] * v[472];
//	v[4434] = v[2698] * (v[2348] + 2e0*v[3742]) - v[4087] * v[473];
//	v[4431] = v[2698] * (v[2344] + 2e0*v[3743]) + v[4087] * v[570];
//	v[4428] = v[2698] * (v[2340] + 2e0*v[3744]) + v[4087] * v[571];
//	v[4425] = v[2698] * (v[2336] + 2e0*v[3745]) + v[4087] * v[572];
//	v[3327] = 2e0*v[2698] * v[319];
//	v[4413] = v[3327] + v[3710] + v[3711];
//	v[3708] = v[3327] / 2e0;
//	v[4415] = v[3705] + v[3708] + v[3711];
//	v[4414] = v[3707] + v[3708] + v[3710];
//	v[4088] = -v[3708] - v[3711];
//	v[4441] = v[2696] * (2e0*v[1218] + v[3757]) - v[4088] * v[465];
//	v[4438] = v[2696] * (2e0*v[1219] + v[3755]) - v[4088] * v[466];
//	v[4435] = v[2696] * (2e0*v[1220] + v[3753]) - v[4088] * v[467];
//	v[4432] = v[2696] * (2e0*v[1221] - v[3751]) + v[4088] * v[564];
//	v[4429] = v[2696] * (2e0*v[1222] - v[3749]) + v[4088] * v[565];
//	v[4426] = v[2696] * (2e0*v[1223] - v[3747]) + v[4088] * v[566];
//	v[4086] = -v[3708] - v[3710];
//	v[4439] = v[2697] * (2e0*v[1229] + v[3756]) - v[4086] * v[468];
//	v[4436] = v[2697] * (2e0*v[1231] + v[3754]) - v[4086] * v[469];
//	v[4433] = v[2697] * (2e0*v[1233] + v[3752]) - v[4086] * v[470];
//	v[4430] = v[2697] * (2e0*v[1235] + v[3750]) + v[4086] * v[567];
//	v[4427] = v[2697] * (2e0*v[1237] + v[3748]) + v[4086] * v[568];
//	v[4424] = v[2697] * (2e0*v[1239] + v[3746]) + v[4086] * v[569];
//	v[3371] = -(v[301] * v[3327]) / 2e0;
//	v[4459] = 2e0*v[3355] + v[3356] + v[3371];
//	v[4453] = v[3355] + 2e0*v[3356] + v[3371];
//	v[4442] = v[3355] + v[3356] + 2e0*v[3371];
//	v[3366] = -(v[291] * v[3327]) / 2e0;
//	v[4458] = 2e0*v[3349] + v[3350] + v[3366];
//	v[4452] = v[3349] + 2e0*v[3350] + v[3366];
//	v[4443] = v[3349] + v[3350] + 2e0*v[3366];
//	v[3364] = -(v[281] * v[3327]) / 2e0;
//	v[4457] = 2e0*v[3346] + v[3347] + v[3364];
//	v[4451] = v[3346] + 2e0*v[3347] + v[3364];
//	v[4444] = v[3346] + v[3347] + 2e0*v[3364];
//	v[3362] = (v[275] * v[3327]) / 2e0;
//	v[4456] = 2e0*v[3343] + v[3344] + v[3362];
//	v[4450] = v[3343] + 2e0*v[3344] + v[3362];
//	v[4445] = v[3343] + v[3344] + 2e0*v[3362];
//	v[3360] = (v[265] * v[3327]) / 2e0;
//	v[4455] = 2e0*v[3340] + v[3341] + v[3360];
//	v[4449] = v[3340] + 2e0*v[3341] + v[3360];
//	v[4446] = v[3340] + v[3341] + 2e0*v[3360];
//	v[3358] = (v[255] * v[3327]) / 2e0;
//	v[4454] = 2e0*v[3337] + v[3338] + v[3358];
//	v[4448] = v[3337] + 2e0*v[3338] + v[3358];
//	v[4447] = v[3337] + v[3338] + 2e0*v[3358];
//	v[2925] = (v[319] * v[3705]) / 2e0 + v[319] * v[3711] + v[2698] * v[717];
//	v[2924] = (v[318] * v[3327]) / 2e0 + (v[318] * v[3705]) / 2e0 + v[2697] * v[715];
//	v[2923] = (v[317] * v[3327]) / 2e0 + v[317] * v[3711] + v[2696] * v[713];
//	v[2706] = 0e0;
//	v[2707] = 0e0;
//	v[2708] = 0e0;
//	v[2709] = 0e0;
//	b2710 = b728;
//	if (b2710) {
//		b2711 = b730;
//		if (b2711) {
//			v[2708] = 0e0;
//			v[2707] = 0e0;
//			v[2706] = 0e0;
//			v[2709] = -v[2699];
//		}
//		else {
//		};
//	}
//	else {
//	};
//	v[2712] = v[1241] * v[2696] + v[1242] * v[2697] + v[1243] * v[2698];
//	v[4031] = -(v[1478] * v[2712]);
//	v[3477] = v[3734] * v[4031];
//	v[3467] = v[3736] * v[4031];
//	v[3460] = v[292] * v[4031];
//	v[2771] = v[1579] * v[2712];
//	v[2767] = v[1580] * v[2712];
//	v[2757] = v[2687] * v[2712];
//	v[2713] = v[1245] * v[2696] + v[1246] * v[2697] + v[1247] * v[2698];
//	v[3478] = v[2713] * v[3802];
//	v[3468] = v[2606] * v[2713];
//	v[4126] = v[3467] + v[3468];
//	v[3461] = v[2713] * v[4032];
//	v[4125] = v[3460] + v[3461];
//	v[2775] = v[1577] * v[2713];
//	v[4163] = v[2757] + v[2775];
//	v[2761] = v[2686] * v[2713];
//	v[4167] = v[2761] + v[2771];
//	v[2714] = v[1249] * v[2696] + v[1250] * v[2697] + v[1251] * v[2698];
//	v[3475] = v[2607] * v[2714];
//	v[4127] = v[3475] + v[3478];
//	v[4467] = v[3477] + v[4127];
//	v[3471] = v[2714] * v[3803];
//	v[4465] = v[3468] + v[3471];
//	v[4464] = v[3471] + v[4126];
//	v[3463] = v[2714] * v[3805];
//	v[4462] = v[3460] + v[3463];
//	v[4460] = v[3463] + v[4125];
//	v[2765] = v[2685] * v[2714];
//	v[4168] = v[2765] + v[2767];
//	v[4039] = v[1578] * v[2713] + v[2765];
//	v[4159] = v[2767] + v[4039];
//	v[4038] = v[2761] + v[2714] * v[3732];
//	v[4162] = v[2771] + v[4038];
//	v[4037] = v[2757] + v[2714] * v[3735];
//	v[4166] = v[2775] + v[4037];
//	v[2715] = v[1253] * v[2696] + v[1254] * v[2697] + v[1255] * v[2698];
//	v[4033] = -(v[1541] * v[2715]);
//	v[3568] = v[3726] * v[4033];
//	v[3558] = v[3728] * v[4033];
//	v[3551] = v[266] * v[4033];
//	v[2804] = v[1573] * v[2715];
//	v[2800] = v[1574] * v[2715];
//	v[2790] = v[2684] * v[2715];
//	v[2716] = v[1257] * v[2696] + v[1258] * v[2697] + v[1259] * v[2698];
//	v[3569] = v[2716] * v[3820];
//	v[3559] = v[2635] * v[2716];
//	v[4129] = v[3558] + v[3559];
//	v[3552] = v[2716] * v[4034];
//	v[4128] = v[3551] + v[3552];
//	v[2808] = v[1571] * v[2716];
//	v[4151] = v[2790] + v[2808];
//	v[2794] = v[2683] * v[2716];
//	v[4155] = v[2794] + v[2804];
//	v[2717] = v[1261] * v[2696] + v[1262] * v[2697] + v[1263] * v[2698];
//	v[3566] = v[2636] * v[2717];
//	v[4130] = v[3566] + v[3569];
//	v[4475] = v[3568] + v[4130];
//	v[3562] = v[2717] * v[3821];
//	v[4473] = v[3559] + v[3562];
//	v[4472] = v[3562] + v[4129];
//	v[3554] = v[2717] * v[3823];
//	v[4470] = v[3551] + v[3554];
//	v[4468] = v[3554] + v[4128];
//	v[2798] = v[2682] * v[2717];
//	v[6874] = v[2923];
//	v[6875] = v[2924];
//	v[6876] = v[2925];
//	v[6877] = v[2798] + v[2716] * v[3724] + v[2715] * v[3727];
//	v[6878] = v[1571] * v[2715] + v[1572] * v[2717] + v[2794];
//	v[6879] = v[1573] * v[2716] + v[1574] * v[2717] + v[2790];
//	v[6880] = -v[2923];
//	v[6881] = -v[2924];
//	v[6882] = -v[2925];
//	v[6883] = v[2765] + v[2713] * v[3732] + v[2712] * v[3735];
//	v[6884] = v[1577] * v[2712] + v[1578] * v[2714] + v[2761];
//	v[6885] = v[1579] * v[2713] + v[1580] * v[2714] + v[2757];
//	v[4156] = v[2798] + v[2800];
//	v[4045] = v[1572] * v[2716] + v[2798];
//	v[4147] = v[2800] + v[4045];
//	v[4044] = v[2794] + v[2717] * v[3724];
//	v[4150] = v[2804] + v[4044];
//	v[4043] = v[2790] + v[2717] * v[3727];
//	v[4154] = v[2808] + v[4043];
//	v[2718] = -(v[291] * v[2923]);
//	v[2721] = -(v[2923] * v[301]);
//	v[2722] = -(v[281] * v[2923]);
//	v[2723] = -(v[281] * v[2924]);
//	v[2726] = -(v[2924] * v[301]);
//	v[2727] = -(v[2925] * v[301]);
//	v[2728] = -(v[291] * v[2924]);
//	v[2731] = -(v[281] * v[2925]);
//	v[2732] = -(v[291] * v[2925]);
//	v[2733] = v[265] * v[2923];
//	v[2734] = v[275] * v[2923];
//	v[2735] = v[255] * v[2923];
//	v[2736] = v[255] * v[2924];
//	v[2737] = v[275] * v[2924];
//	v[2738] = v[275] * v[2925];
//	v[2739] = v[265] * v[2924];
//	v[2740] = v[255] * v[2925];
//	v[2741] = v[265] * v[2925];
//	v[2708] = v[2103] * v[2696] + v[2098] * v[2697] + v[2093] * v[2698] + v[2708] + v[2090] * v[3327];
//	v[2707] = v[2102] * v[2696] + v[2097] * v[2697] + v[2092] * v[2698] + v[2707] + v[2095] * v[3707];
//	v[2706] = v[2101] * v[2696] + v[2096] * v[2697] + v[2091] * v[2698] + v[2706] + v[2100] * v[3705];
//	v[3373] = v[2706] * v[302] + v[2707] * v[303] + v[2708] * v[304];
//	v[2709] = v[2709] + v[315] * v[3373];
//	v[4035] = v[2709] / v[727];
//	v[2742] = v[2708] * v[316] + v[304] * v[4035] + v[751];
//	v[2743] = v[2707] * v[316] + v[303] * v[4035] + v[750];
//	v[2744] = v[2706] * v[316] + v[302] * v[4035] + v[749];
//	v[2745] = v[2039] * v[2712];
//	v[2746] = v[2038] * v[2712];
//	v[2747] = v[2037] * v[2712];
//	v[2748] = v[2034] * v[2713];
//	v[2749] = v[2712] * v[288] + v[2713] * v[290];
//	v[4036] = v[209] * v[2749];
//	v[3661] = -(v[1478] * v[4036]);
//	v[3658] = v[2749] * v[3802];
//	v[4466] = v[290] * (v[3475] + v[3477]) + v[3658];
//	v[2750] = v[2032] * v[2713];
//	v[2751] = v[2745] + v[2748];
//	v[2752] = v[2033] * v[2713] + v[4036] / v[279];
//	v[2753] = v[2027] * v[2714];
//	v[2754] = v[2712] * v[282] + v[2714] * v[290];
//	v[4040] = v[214] * v[2754];
//	v[3660] = -(v[1478] * v[4040]);
//	v[3657] = v[2754] * v[3803];
//	v[4463] = v[3657] + v[290] * v[4126];
//	v[2755] = v[2713] * v[282] + v[2714] * v[288];
//	v[4041] = v[218] * v[2755];
//	v[4479] = -(v[2712] * v[3653]) - v[2713] * v[3654] - v[2714] * v[3655] + v[3731] * v[4036] + v[278] * v[4040]
//		+ v[277] * v[4041];
//	v[3659] = -(v[1478] * v[4041]);
//	v[3656] = v[2755] * v[3805];
//	v[4461] = v[3656] + v[288] * v[4125];
//	v[3651] = v[2036] * v[2712] + v[2031] * v[2713] + v[2026] * v[2714] + v[2747] + v[1859] * v[2749] + v[2750] + v[2753]
//		+ v[1857] * v[2754] + v[1853] * v[2755];
//	v[2756] = v[290] * v[4166] + v[2742] * v[851];
//	v[2759] = v[282] * v[4037] + v[2742] * v[855];
//	v[2760] = v[288] * v[4162] + v[2743] * v[853];
//	v[2763] = v[282] * v[4038] + v[2743] * v[855];
//	v[2764] = v[1578] * v[2749] + v[290] * v[4168] + v[2744] * v[851];
//	v[2766] = v[288] * v[4039] + v[2744] * v[853];
//	v[2769] = v[282] * v[4159] + v[2744] * v[855];
//	v[2770] = v[2029] * v[2714] + v[4040] / v[279];
//	v[2772] = v[2754] * v[3732] + v[290] * v[4167] + v[2743] * v[851];
//	v[2773] = v[2752] + v[2770];
//	v[2774] = v[2028] * v[2714] + v[4041] / v[279];
//	v[2776] = v[2755] * v[3735] + v[288] * v[4163] + v[2742] * v[853];
//	v[2777] = v[2746] + v[2774];
//	v[2778] = v[2024] * v[2715];
//	v[2779] = v[2023] * v[2715];
//	v[2780] = v[2022] * v[2715];
//	v[2781] = v[2019] * v[2716];
//	v[2782] = v[262] * v[2715] + v[264] * v[2716];
//	v[4042] = v[153] * v[2782];
//	v[3686] = -(v[1541] * v[4042]);
//	v[3683] = v[2782] * v[3820];
//	v[4474] = v[264] * (v[3566] + v[3568]) + v[3683];
//	v[2783] = v[2017] * v[2716];
//	v[2784] = v[2778] + v[2781];
//	v[2785] = v[2018] * v[2716] + v[4042] / v[253];
//	v[2786] = v[2012] * v[2717];
//	v[2787] = v[256] * v[2715] + v[264] * v[2717];
//	v[4046] = v[158] * v[2787];
//	v[3685] = -(v[1541] * v[4046]);
//	v[3682] = v[2787] * v[3821];
//	v[4471] = v[3682] + v[264] * v[4129];
//	v[2788] = v[256] * v[2716] + v[262] * v[2717];
//	v[4047] = v[162] * v[2788];
//	v[4483] = -(v[2715] * v[3678]) - v[2716] * v[3679] - v[2717] * v[3680] + v[3723] * v[4042] + v[252] * v[4046]
//		+ v[251] * v[4047];
//	v[3684] = -(v[1541] * v[4047]);
//	v[3681] = v[2788] * v[3823];
//	v[4469] = v[3681] + v[262] * v[4128];
//	v[3676] = v[2021] * v[2715] + v[2016] * v[2716] + v[2011] * v[2717] + v[2780] + v[1831] * v[2782] + v[2783] + v[2786]
//		+ v[1829] * v[2787] + v[1825] * v[2788];
//	v[2789] = v[264] * v[4154] + v[2742] * v[857];
//	v[2792] = v[256] * v[4043] + v[2742] * v[861];
//	v[2793] = v[262] * v[4150] + v[2743] * v[859];
//	v[2796] = v[256] * v[4044] + v[2743] * v[861];
//	v[2797] = v[1572] * v[2782] + v[264] * v[4156] + v[2744] * v[857];
//	v[2799] = v[262] * v[4045] + v[2744] * v[859];
//	v[2802] = v[256] * v[4147] + v[2744] * v[861];
//	v[2803] = v[2014] * v[2717] + v[4046] / v[253];
//	v[2805] = v[2787] * v[3724] + v[264] * v[4155] + v[2743] * v[857];
//	v[2806] = v[2785] + v[2803];
//	v[2807] = v[2013] * v[2717] + v[4047] / v[253];
//	v[2809] = v[2788] * v[3727] + v[262] * v[4151] + v[2742] * v[859];
//	v[2810] = v[2779] + v[2807];
//	v[2811] = -(v[2718] * v[855]);
//	v[2812] = -(v[2718] * v[853]);
//	v[2813] = -(v[2718] * v[851]);
//	v[2814] = -(v[2721] * v[855]);
//	v[2815] = -(v[2721] * v[851]);
//	v[2816] = -(v[2721] * v[853]);
//	v[4050] = -2e0*v[2816];
//	v[2817] = -(v[2722] * v[853]);
//	v[2818] = -(v[2722] * v[851]);
//	v[2819] = -(v[2722] * v[855]);
//	v[2820] = -(v[2723] * v[855]);
//	v[2821] = -(v[2723] * v[853]);
//	v[2822] = -(v[2723] * v[851]);
//	v[4048] = -2e0*v[2822];
//	v[2823] = -(v[2726] * v[851]);
//	v[2824] = -(v[2726] * v[855]);
//	v[4051] = 2e0*v[2824];
//	v[2825] = -(v[2726] * v[853]);
//	v[2826] = -(v[2727] * v[853]);
//	v[2827] = -(v[2727] * v[855]);
//	v[2828] = -(v[2727] * v[851]);
//	v[2829] = v[2817] + v[2820] + v[2823] + v[2826] + 2e0*v[2750] * v[481] - v[2773] * v[484] - v[2751] * v[487];
//	v[2830] = -(v[2728] * v[855]);
//	v[2831] = -(v[2728] * v[851]);
//	v[2832] = -(v[2728] * v[853]);
//	v[2833] = -v[2812] - v[2815] - v[2827] - v[2830] - v[2773] * v[481] + 2e0*v[2753] * v[484] + v[2777] * v[487];
//	v[2834] = v[204] * v[2766] - v[2817] * v[474] + v[2812] * v[475] - v[2816] * v[476];
//	v[2835] = -(v[2731] * v[855]);
//	v[2836] = -(v[2731] * v[853]);
//	v[4049] = 2e0*v[2836];
//	v[2837] = -(v[2731] * v[851]);
//	v[2838] = -(v[2732] * v[853]);
//	v[2839] = -(v[2732] * v[855]);
//	v[2840] = -(v[2732] * v[851]);
//	v[2841] = -(v[231] * v[2742]) - v[228] * v[2743] - v[225] * v[2744] + v[2722] * v[525] + v[2718] * v[526]
//		+ v[2721] * v[527] + v[2723] * v[540] + v[2728] * v[541] + v[2726] * v[542] + v[2731] * v[555] + v[2732] * v[556]
//		+ v[2727] * v[557];
//	v[2842] = -(v[230] * v[2742]) - v[227] * v[2743] - v[224] * v[2744] + v[2722] * v[522] + v[2718] * v[523]
//		+ v[2721] * v[524] + v[2723] * v[537] + v[2728] * v[538] + v[2726] * v[539] + v[2731] * v[552] + v[2732] * v[553]
//		+ v[2727] * v[554];
//	v[2843] = -(v[229] * v[2742]) - v[226] * v[2743] - v[223] * v[2744] + v[2722] * v[519] + v[2718] * v[520]
//		+ v[2721] * v[521] + v[2723] * v[534] + v[2728] * v[535] + v[2726] * v[536] + v[2731] * v[549] + v[2732] * v[550]
//		+ v[2727] * v[551];
//	v[2844] = -v[2818] - v[2831] - v[2835] - v[2838] - v[2751] * v[481] + v[2777] * v[484] + 2e0*v[2747] * v[487];
//	v[4133] = -v[2844] / 2e0;
//	v[2845] = v[204] * v[2764] - v[2818] * v[474] + v[2813] * v[475] - v[2815] * v[476];
//	v[2846] = v[204] * v[2763] - v[2820] * v[474] + v[2830] * v[475] - v[2824] * v[476];
//	v[2847] = v[2814] + v[2825];
//	v[2848] = v[204] * v[2772] - v[2822] * v[474] + v[2831] * v[475] - v[2823] * v[476];
//	v[2849] = v[204] * v[2759] - v[2835] * v[474] + v[2839] * v[475] - v[2827] * v[476];
//	v[2850] = v[204] * v[2776] - v[2836] * v[474] + v[2838] * v[475] - v[2826] * v[476];
//	v[2851] = v[2821] + v[2837];
//	v[2852] = v[2811] + v[2840];
//	v[7034] = 0e0;
//	v[7035] = 0e0;
//	v[7036] = 0e0;
//	v[7037] = 0e0;
//	v[7038] = 0e0;
//	v[7039] = 0e0;
//	v[7040] = 0e0;
//	v[7041] = 0e0;
//	v[7042] = 0e0;
//	v[7043] = -v[2833] / 2e0 - v[2851];
//	v[7044] = (v[2829] - 2e0*v[2852]) / 2e0;
//	v[7045] = -v[2847] + v[4133];
//	v[2853] = (2e0*v[2813] + alphaB[1] * v[2829] - alphaB[0] * v[2833] - 2e0*v[2839] - alphaB[2] * v[2844]
//		+ 4e0*v[204] * v[3651] - v[2851] * v[369] - v[2852] * v[372] - v[2847] * v[374] + v[4048] + v[4049] + v[4050] + v[4051]
//		+ v[2769] * v[477] + 2e0*v[2766] * v[482] + 2e0*v[2764] * v[488] + 2e0*v[2763] * v[492] + v[2760] * v[496]
//		+ 2e0*v[2772] * v[501] + 2e0*v[2759] * v[505] + 2e0*v[2776] * v[509] + v[2756] * v[514]) / 2e0;
//	v[2854] = v[2843] * x1B[0] + v[2842] * x1B[1] + v[2841] * x1B[2];
//	v[2855] = v[2733] * v[861];
//	v[2856] = v[2733] * v[859];
//	v[2857] = v[2733] * v[857];
//	v[2858] = v[2734] * v[861];
//	v[2859] = v[2734] * v[857];
//	v[2860] = v[2734] * v[859];
//	v[4054] = -2e0*v[2860];
//	v[2861] = v[2735] * v[859];
//	v[2862] = v[2735] * v[857];
//	v[2863] = v[2735] * v[861];
//	v[2864] = v[2736] * v[861];
//	v[2865] = v[2736] * v[859];
//	v[2866] = v[2736] * v[857];
//	v[4052] = -2e0*v[2866];
//	v[2867] = v[2737] * v[857];
//	v[2868] = v[2737] * v[861];
//	v[4055] = 2e0*v[2868];
//	v[2869] = v[2737] * v[859];
//	v[2870] = v[2738] * v[859];
//	v[2871] = v[2738] * v[861];
//	v[2872] = v[2738] * v[857];
//	v[2873] = v[2861] + v[2864] + v[2867] + v[2870] + 2e0*v[2783] * v[382] - v[2806] * v[385] - v[2784] * v[388];
//	v[2874] = v[2739] * v[861];
//	v[2875] = v[2739] * v[857];
//	v[2876] = v[2739] * v[859];
//	v[2877] = -v[2856] - v[2859] - v[2871] - v[2874] - v[2806] * v[382] + 2e0*v[2786] * v[385] + v[2810] * v[388];
//	v[2878] = v[148] * v[2799] - v[2861] * v[375] + v[2856] * v[376] - v[2860] * v[377];
//	v[2879] = v[2740] * v[861];
//	v[2880] = v[2740] * v[859];
//	v[4053] = 2e0*v[2880];
//	v[2881] = v[2740] * v[857];
//	v[2882] = v[2741] * v[859];
//	v[2883] = v[2741] * v[861];
//	v[2884] = v[2741] * v[857];
//	v[2885] = v[175] * v[2742] + v[172] * v[2743] + v[169] * v[2744] + v[2735] * v[426] + v[2733] * v[427] + v[2734] * v[428]
//		+ v[2736] * v[441] + v[2739] * v[442] + v[2737] * v[443] + v[2740] * v[456] + v[2741] * v[457] + v[2738] * v[458];
//	v[2886] = v[174] * v[2742] + v[171] * v[2743] + v[168] * v[2744] + v[2735] * v[423] + v[2733] * v[424] + v[2734] * v[425]
//		+ v[2736] * v[438] + v[2739] * v[439] + v[2737] * v[440] + v[2740] * v[453] + v[2741] * v[454] + v[2738] * v[455];
//	v[2887] = v[173] * v[2742] + v[170] * v[2743] + v[167] * v[2744] + v[2735] * v[420] + v[2733] * v[421] + v[2734] * v[422]
//		+ v[2736] * v[435] + v[2739] * v[436] + v[2737] * v[437] + v[2740] * v[450] + v[2741] * v[451] + v[2738] * v[452];
//	v[2888] = -v[2862] - v[2875] - v[2879] - v[2882] - v[2784] * v[382] + v[2810] * v[385] + 2e0*v[2780] * v[388];
//	v[4131] = -v[2888] / 2e0;
//	v[2889] = v[148] * v[2797] - v[2862] * v[375] + v[2857] * v[376] - v[2859] * v[377];
//	v[2890] = v[148] * v[2796] - v[2864] * v[375] + v[2874] * v[376] - v[2868] * v[377];
//	v[2891] = v[2858] + v[2869];
//	v[2892] = v[148] * v[2805] - v[2866] * v[375] + v[2875] * v[376] - v[2867] * v[377];
//	v[2893] = v[148] * v[2792] - v[2879] * v[375] + v[2883] * v[376] - v[2871] * v[377];
//	v[2894] = v[148] * v[2809] - v[2880] * v[375] + v[2882] * v[376] - v[2870] * v[377];
//	v[2895] = v[2865] + v[2881];
//	v[2896] = v[2855] + v[2884];
//	v[7046] = 0e0;
//	v[7047] = 0e0;
//	v[7048] = 0e0;
//	v[7049] = -v[2877] / 2e0 - v[2895];
//	v[7050] = (v[2873] - 2e0*v[2896]) / 2e0;
//	v[7051] = -v[2891] + v[4131];
//	v[7052] = 0e0;
//	v[7053] = 0e0;
//	v[7054] = 0e0;
//	v[7055] = 0e0;
//	v[7056] = 0e0;
//	v[7057] = 0e0;
//	v[2897] = (2e0*v[2857] + alphaA[1] * v[2873] - alphaA[0] * v[2877] - 2e0*v[2883] - alphaA[2] * v[2888]
//		- v[2895] * v[363] - v[2896] * v[366] + 4e0*v[148] * v[3676] - v[2891] * v[368] + v[2802] * v[378] + 2e0*v[2799] * v[383]
//		+ 2e0*v[2797] * v[389] + 2e0*v[2796] * v[393] + v[2793] * v[397] + 2e0*v[2805] * v[402] + v[4052] + v[4053] + v[4054]
//		+ v[4055] + 2e0*v[2792] * v[406] + 2e0*v[2809] * v[410] + v[2789] * v[415]) / 2e0;
//	v[2898] = v[2887] * x1A[0] + v[2886] * x1A[1] + v[2885] * x1A[2];
//	v[2899] = (-v[2854] + v[2843] * x3B[0] + v[2842] * x3B[1] + v[2841] * x3B[2]) / 2e0;
//	v[2900] = (-v[2854] + v[2843] * x2B[0] + v[2842] * x2B[1] + v[2841] * x2B[2]) / 2e0;
//	v[2901] = (-2e0*v[2745] + 2e0*v[2748] - v[2819] * v[477] - 2e0*v[2817] * v[482] - 2e0*v[2818] * v[488]
//		- 2e0*v[2820] * v[492] - v[2821] * v[496] + v[4048] * v[501] - 2e0*v[2835] * v[505] - v[4049] * v[509] - v[2837] * v[514])
//		/ 2e0;
//	v[4144] = 8e0*v[2901];
//	v[2903] = -v[2746] + v[2774] + (v[2811] * v[477]) / 2e0 + v[2812] * v[482] + v[2813] * v[488] + v[2830] * v[492] +
//		(v[2832] * v[496]) / 2e0 + v[2831] * v[501] + v[2839] * v[505] + v[2838] * v[509] + (v[2840] * v[514]) / 2e0;
//	v[4143] = 8e0*v[2903];
//	v[2904] = (v[204] * v[2760] - v[2821] * v[474] + v[2832] * v[475] - v[2825] * v[476]) / 2e0;
//	v[2905] = (-2e0*v[2752] + 2e0*v[2770] - v[2814] * v[477] + v[4050] * v[482] - 2e0*v[2815] * v[488] - v[4051] * v[492]
//		- v[2825] * v[496] - 2e0*v[2823] * v[501] - 2e0*v[2827] * v[505] - 2e0*v[2826] * v[509] - v[2828] * v[514]) / 2e0;
//	v[4477] = v[2901] * v[369] - v[2903] * v[372] + v[2905] * v[374];
//	v[4142] = 8e0*v[2905];
//	v[7058] = 0e0;
//	v[7059] = 0e0;
//	v[7060] = 0e0;
//	v[7061] = 0e0;
//	v[7062] = 0e0;
//	v[7063] = 0e0;
//	v[7064] = 0e0;
//	v[7065] = 0e0;
//	v[7066] = 0e0;
//	v[7067] = v[4144];
//	v[7068] = -v[4143];
//	v[7069] = v[4142];
//	v[2922] = v[1660] * v[2901] + v[1658] * v[2903] - v[2904] + v[1655] * v[2905] - v[2853] * v[4134];
//	v[3704] = v[2922] + (-(v[204] * v[2769]) + v[2819] * v[474] - v[2811] * v[475] + v[2814] * v[476]) / 2e0;
//	v[2906] = (v[204] * v[2756] - v[2837] * v[474] + v[2840] * v[475] - v[2828] * v[476]) / 2e0;
//	v[3702] = v[2904] - v[2906] + v[3704];
//	v[3698] = -v[2906] + v[2922];
//	v[3700] = (v[2845] + v[2849]) / 2e0;
//	v[2908] = v[2848] + v[2850];
//	v[3703] = v[2908] / 2e0;
//	v[2909] = v[2834] + v[2846];
//	v[3699] = v[2909] / 2e0;
//	v[2910] = (-v[2898] + v[2887] * x2A[0] + v[2886] * x2A[1] + v[2885] * x2A[2]) / 2e0;
//	v[2911] = (-v[2898] + v[2887] * x3A[0] + v[2886] * x3A[1] + v[2885] * x3A[2]) / 2e0;
//	v[2912] = (-2e0*v[2778] + 2e0*v[2781] - v[2863] * v[378] - 2e0*v[2861] * v[383] - 2e0*v[2862] * v[389]
//		- 2e0*v[2864] * v[393] - v[2865] * v[397] + v[402] * v[4052] - 2e0*v[2879] * v[406] - v[4053] * v[410] - v[2881] * v[415])
//		/ 2e0;
//	v[4138] = 8e0*v[2912];
//	v[2914] = -v[2779] + v[2807] + (v[2855] * v[378]) / 2e0 + v[2856] * v[383] + v[2857] * v[389] + v[2874] * v[393] +
//		(v[2876] * v[397]) / 2e0 + v[2875] * v[402] + v[2883] * v[406] + v[2882] * v[410] + (v[2884] * v[415]) / 2e0;
//	v[4137] = 8e0*v[2914];
//	v[2915] = (v[148] * v[2793] - v[2865] * v[375] + v[2876] * v[376] - v[2869] * v[377]) / 2e0;
//	v[2916] = (-2e0*v[2785] + 2e0*v[2803] - v[2858] * v[378] - 2e0*v[2859] * v[389] - v[2869] * v[397]
//		- 2e0*v[2867] * v[402] + v[383] * v[4054] - v[393] * v[4055] - 2e0*v[2871] * v[406] - 2e0*v[2870] * v[410]
//		- v[2872] * v[415]) / 2e0;
//	v[4481] = v[2912] * v[363] - v[2914] * v[366] + v[2916] * v[368];
//	v[4136] = 8e0*v[2916];
//	v[7106] = 0e0;
//	v[7107] = 0e0;
//	v[7108] = 0e0;
//	v[7109] = v[4138];
//	v[7110] = -v[4137];
//	v[7111] = v[4136];
//	v[7112] = 0e0;
//	v[7113] = 0e0;
//	v[7114] = 0e0;
//	v[7115] = 0e0;
//	v[7116] = 0e0;
//	v[7117] = 0e0;
//	v[2921] = v[1634] * v[2912] + v[1632] * v[2914] - v[2915] + v[1629] * v[2916] - v[2897] * v[4132];
//	v[3697] = v[2921] + (-(v[148] * v[2802]) + v[2863] * v[375] - v[2855] * v[376] + v[2858] * v[377]) / 2e0;
//	v[2917] = (v[148] * v[2789] - v[2881] * v[375] + v[2884] * v[376] - v[2872] * v[377]) / 2e0;
//	v[3695] = v[2915] - v[2917] + v[3697];
//	v[3691] = -v[2917] + v[2921];
//	v[3693] = (v[2889] + v[2893]) / 2e0;
//	v[2919] = v[2892] + v[2894];
//	v[3696] = v[2919] / 2e0;
//	v[2920] = v[2878] + v[2890];
//	v[6862] = v[2744];
//	v[6863] = v[2743];
//	v[6864] = v[2742];
//	v[6865] = -v[2892] + v[2894] + v[2920] * v[364] + 2e0*alphaA[0] * v[3691] + alphaA[2] * v[3693] + v[2877] * v[3983]
//		+ 2e0*(v[2895] * v[3983] + v[2912] * v[3985]);
//	v[6866] = (v[148] * v[2873]) / 2e0 + v[2889] - v[2893] + (alphaA[2] * v[2919]) / 2e0 + (alphaA[0] * v[2920]) / 2e0
//		+ 2e0*alphaA[1] * v[3695] + 2e0*(v[2896] * v[3983] - v[2914] * v[3985]);
//	v[6867] = -v[2878] + v[2890] + v[2919] * v[364] + alphaA[0] * v[3693] + 2e0*alphaA[2] * v[3697] + v[148] * v[4131] + 2e0*
//		(v[2891] * v[3983] + v[2916] * v[4132]);
//	v[6868] = -v[2744];
//	v[6869] = -v[2743];
//	v[6870] = -v[2742];
//	v[6871] = -v[2848] + v[2850] + 2e0*alphaB[0] * v[3698] + v[2909] * v[370] + alphaB[2] * v[3700] + v[2833] * v[3986]
//		+ 2e0*(v[2851] * v[3986] + v[2901] * v[3988]);
//	v[6872] = (v[204] * v[2829]) / 2e0 + v[2845] - v[2849] + (alphaB[2] * v[2908]) / 2e0 + (alphaB[0] * v[2909]) / 2e0
//		+ 2e0*alphaB[1] * v[3702] + 2e0*(v[2852] * v[3986] - v[2903] * v[3988]);
//	v[6873] = -v[2834] + v[2846] + v[2908] * v[370] + alphaB[0] * v[3700] + 2e0*alphaB[2] * v[3704] + v[204] * v[4133] + 2e0*
//		(v[2847] * v[3986] + v[2905] * v[4134]);
//	v[3692] = v[2920] / 2e0;
//	for (i2691 = 1; i2691 <= 12; i2691++) {
//		i4064 = (i2691 == 11 ? 1 : 0);
//		v[4107] = (*a4)*i4064;
//		v[4067] = i4064 * v[204];
//		i4063 = (i2691 == 10 ? 1 : 0);
//		v[4108] = (*a4)*i4063;
//		v[4066] = -(i4063*v[204]);
//		i4062 = (i2691 == 12 ? 1 : 0);
//		v[4112] = (*a4)*i4062;
//		v[4065] = -(i4062*v[204]);
//		i4058 = (i2691 == 5 ? 1 : 0);
//		v[4116] = (*a4)*i4058;
//		v[4061] = i4058 * v[148];
//		i4057 = (i2691 == 4 ? 1 : 0);
//		v[4117] = (*a4)*i4057;
//		v[4060] = -(i4057*v[148]);
//		i4056 = (i2691 == 6 ? 1 : 0);
//		v[4121] = (*a4)*i4056;
//		v[4059] = -(i4056*v[148]);
//		v[3172] = v[5757 + i2691];
//		v[3163] = v[5745 + i2691];
//		v[3157] = v[5733 + i2691];
//		v[2950] = v[5381 + i2691];
//		v[2949] = v[5393 + i2691];
//		v[2947] = v[5357 + i2691];
//		v[2946] = v[5369 + i2691];
//		v[2930] = v[5433 + i2691];
//		v[2931] = v[5457 + i2691];
//		v[2932] = v[5445 + i2691];
//		v[2934] = v[5885 + i2691];
//		v[2935] = v[5421 + i2691];
//		v[2972] = 2e0*v[2935] * v[962];
//		v[3041] = -4e0*v[148] * v[2972];
//		v[4070] = v[253] * v[3041];
//		v[3016] = -2e0*v[2972];
//		v[2937] = v[5933 + i2691];
//		v[2938] = v[5481 + i2691];
//		v[2939] = v[5505 + i2691];
//		v[2940] = v[5493 + i2691];
//		v[2942] = v[5957 + i2691];
//		v[2943] = v[5469 + i2691];
//		v[2998] = 2e0*v[2943] * v[968];
//		v[3091] = -4e0*v[204] * v[2998];
//		v[4071] = v[279] * v[3091];
//		v[3066] = -2e0*v[2998];
//		v[2945] = v[6005 + i2691];
//		v[2948] = (-v[2946] - v[2947]) / 2e0;
//		v[2951] = (-v[2949] - v[2950]) / 2e0;
//		v[2955] = i4056 + v[2930];
//		v[2957] = -i4056 + v[2930];
//		v[2958] = i4057 + v[2931];
//		v[2960] = -i4057 + v[2931];
//		v[2961] = -i4058 + v[2932];
//		v[2963] = i4058 + v[2932];
//		v[2964] = v[1629] * v[2935] + 8e0*i4056*v[962];
//		v[2965] = 2e0*alphaA[1] * i4058 - v[2935];
//		v[2966] = v[1632] * v[2935] - 8e0*i4058*v[962];
//		v[2967] = v[1634] * v[2935] + 8e0*i4057*v[962];
//		v[4068] = 2e0*(-(v[2935] * v[376]) / 2e0 - v[4061]);
//		v[2969] = (v[2935] * v[375]) / 2e0 + v[4060];
//		v[2970] = (v[2935] * v[377]) / 2e0 + v[4059];
//		v[2971] = alphaA[2] * v[2972] + v[4059] / 2e0;
//		v[2973] = alphaA[0] * v[2972] + v[4060] / 2e0;
//		v[2974] = -(alphaA[1] * v[2972]) + v[4061] / 2e0;
//		v[2975] = (v[148] * v[2934]) / 2e0 - v[2972] * v[415];
//		v[3126] = v[264] * v[2975];
//		v[2976] = (-(v[2934] * v[377]) - v[2964] * v[415]) / 2e0;
//		v[2977] = (v[148] * v[2965]) / 2e0 - v[2972] * v[397];
//		v[3120] = v[262] * v[2977];
//		v[2978] = (v[2965] * v[376] + v[2966] * v[397]) / 2e0;
//		v[2979] = (v[148] * v[2937]) / 2e0 - v[2972] * v[378];
//		v[3115] = v[256] * v[2979];
//		v[2980] = (-(v[2937] * v[375]) - v[2967] * v[378]) / 2e0;
//		v[2981] = i4062 + v[2938];
//		v[2983] = -i4062 + v[2938];
//		v[2984] = i4063 + v[2939];
//		v[2986] = -i4063 + v[2939];
//		v[2987] = -i4064 + v[2940];
//		v[2989] = i4064 + v[2940];
//		v[2990] = v[1655] * v[2943] + 8e0*i4062*v[968];
//		v[2991] = 2e0*alphaB[1] * i4064 - v[2943];
//		v[2992] = v[1658] * v[2943] - 8e0*i4064*v[968];
//		v[2993] = v[1660] * v[2943] + 8e0*i4063*v[968];
//		v[4069] = 2e0*(-v[4067] - (v[2943] * v[475]) / 2e0);
//		v[2995] = v[4066] + (v[2943] * v[474]) / 2e0;
//		v[2996] = v[4065] + (v[2943] * v[476]) / 2e0;
//		v[2997] = alphaB[2] * v[2998] + v[4065] / 2e0;
//		v[2999] = alphaB[0] * v[2998] + v[4066] / 2e0;
//		v[3000] = -(alphaB[1] * v[2998]) + v[4067] / 2e0;
//		v[3001] = (v[204] * v[2942]) / 2e0 - v[2998] * v[514];
//		v[3175] = v[290] * v[3001];
//		v[3002] = (-(v[2942] * v[476]) - v[2990] * v[514]) / 2e0;
//		v[3003] = (v[204] * v[2991]) / 2e0 - v[2998] * v[496];
//		v[3167] = v[288] * v[3003];
//		v[3004] = (v[2991] * v[475] + v[2992] * v[496]) / 2e0;
//		v[3005] = (v[204] * v[2945]) / 2e0 - v[2998] * v[477];
//		v[3160] = v[282] * v[3005];
//		v[3006] = (-(v[2945] * v[474]) - v[2993] * v[477]) / 2e0;
//		v[3007] = (2e0*v[2948] * x1A[0] + v[2946] * x2A[0] + v[2947] * x3A[0]) / 2e0;
//		v[3008] = (2e0*v[2948] * x1A[1] + v[2946] * x2A[1] + v[2947] * x3A[1]) / 2e0;
//		v[3009] = (2e0*v[2948] * x1A[2] + v[2946] * x2A[2] + v[2947] * x3A[2]) / 2e0;
//		v[3010] = (v[2934] * v[376] + v[4068] + v[2966] * v[415]) / 2e0;
//		v[3011] = (v[2937] * v[376] + v[2966] * v[378] + v[4068]) / 2e0;
//		v[3012] = v[2969] - (v[2934] * v[375]) / 2e0 - (v[2967] * v[415]) / 2e0;
//		v[3013] = v[2969] - (v[2965] * v[375]) / 2e0 - (v[2967] * v[397]) / 2e0;
//		v[3014] = v[3016] - v[2958] * v[375] - v[2967] * v[410];
//		v[3015] = v[148] * v[2958] - 2e0*v[2972] * v[410];
//		v[3017] = -v[3016] + v[2961] * v[376] + v[2966] * v[406];
//		v[3018] = v[148] * v[2961] - 2e0*v[2972] * v[406];
//		v[3128] = v[256] * v[3018];
//		v[3019] = -v[3016] - v[2960] * v[375] - v[2967] * v[402];
//		v[3020] = v[148] * v[2960] - 2e0*v[2972] * v[402];
//		v[3121] = v[264] * v[3020];
//		v[3021] = v[2970] - (v[2965] * v[377]) / 2e0 - (v[2964] * v[397]) / 2e0;
//		v[3022] = v[2970] - (v[2937] * v[377]) / 2e0 - (v[2964] * v[378]) / 2e0;
//		v[3023] = v[3016] - v[2955] * v[377] - v[2964] * v[393];
//		v[3024] = v[148] * v[2955] - 2e0*v[2972] * v[393];
//		v[3025] = v[3016] + v[2963] * v[376] + v[2966] * v[389];
//		v[3026] = v[148] * v[2963] - 2e0*v[2972] * v[389];
//		v[3116] = v[264] * v[3026];
//		v[3027] = -v[2971] + v[2958] * v[376] + v[2966] * v[410];
//		v[3028] = -v[2971] - v[2961] * v[375] - v[2967] * v[406];
//		v[3029] = -v[2971] + v[2960] * v[376] + v[2966] * v[402];
//		v[3030] = -v[2971] - v[2963] * v[375] - v[2967] * v[389];
//		v[3031] = v[3041] + 2e0*v[2971] * v[388];
//		v[3032] = v[3007] * v[451] + v[3008] * v[454] + v[3009] * v[457] + v[3010] * v[857] + v[3027] * v[859] + v[3017] * v[861];
//		v[3033] = v[3007] * v[450] + v[3008] * v[453] + v[3009] * v[456] + v[3012] * v[857] + v[3014] * v[859] + v[3028] * v[861];
//		v[3034] = -v[3016] - v[2957] * v[377] - v[2964] * v[383];
//		v[3035] = v[148] * v[2957] - 2e0*v[2972] * v[383];
//		v[3036] = -v[2973] + v[2955] * v[376] + v[2966] * v[393];
//		v[3037] = -v[2973] - v[2961] * v[377] - v[2964] * v[406];
//		v[3038] = -v[2973] - v[2963] * v[377] - v[2964] * v[389];
//		v[3039] = -v[2973] + v[2957] * v[376] + v[2966] * v[383];
//		v[3040] = v[2971] * v[385] + v[2973] * v[388];
//		v[3042] = v[3041] + 2e0*v[2973] * v[385];
//		v[3142] = v[2717] * v[3042];
//		v[3043] = v[3007] * v[436] + v[3008] * v[439] + v[3009] * v[442] + v[3029] * v[857] + v[2978] * v[859] + v[3036] * v[861];
//		v[3044] = v[2974] - v[2958] * v[377] - v[2964] * v[410];
//		v[3045] = v[2974] - v[2960] * v[377] - v[2964] * v[402];
//		v[3046] = v[2974] - v[2955] * v[375] - v[2967] * v[393];
//		v[3047] = v[2974] - v[2957] * v[375] - v[2967] * v[383];
//		v[3048] = -(v[2973] * v[382]) - v[2974] * v[385];
//		v[3049] = -(v[2971] * v[382]) - v[2974] * v[388];
//		v[3050] = v[3041] + 2e0*v[2974] * v[382];
//		v[3146] = v[2716] * v[3050];
//		v[3051] = v[3007] * v[452] + v[3008] * v[455] + v[3009] * v[458] + v[2976] * v[857] + v[3044] * v[859] + v[3037] * v[861];
//		v[3052] = v[3007] * v[437] + v[3008] * v[440] + v[3009] * v[443] + v[3045] * v[857] + v[3021] * v[859] + v[3023] * v[861];
//		v[3053] = v[3007] * v[435] + v[3008] * v[438] + v[3009] * v[441] + v[3019] * v[857] + v[3013] * v[859] + v[3046] * v[861];
//		v[3054] = v[3007] * v[420] + v[3008] * v[423] + v[3009] * v[426] + v[3030] * v[857] + v[3047] * v[859] + v[2980] * v[861];
//		v[3055] = v[3007] * v[422] + v[3008] * v[425] + v[3009] * v[428] + v[3038] * v[857] + v[3034] * v[859] + v[3022] * v[861];
//		v[3056] = v[3007] * v[421] + v[3008] * v[424] + v[3009] * v[427] + v[3025] * v[857] + v[3039] * v[859] + v[3011] * v[861];
//		v[3057] = (2e0*v[2951] * x1B[0] + v[2950] * x2B[0] + v[2949] * x3B[0]) / 2e0;
//		v[3058] = (2e0*v[2951] * x1B[1] + v[2950] * x2B[1] + v[2949] * x3B[1]) / 2e0;
//		v[3059] = (2e0*v[2951] * x1B[2] + v[2950] * x2B[2] + v[2949] * x3B[2]) / 2e0;
//		v[3060] = (v[4069] + v[2942] * v[475] + v[2992] * v[514]) / 2e0;
//		v[3061] = (v[4069] + v[2945] * v[475] + v[2992] * v[477]) / 2e0;
//		v[3062] = v[2995] - (v[2942] * v[474]) / 2e0 - (v[2993] * v[514]) / 2e0;
//		v[3063] = v[2995] - (v[2991] * v[474]) / 2e0 - (v[2993] * v[496]) / 2e0;
//		v[3064] = v[3066] - v[2984] * v[474] - v[2993] * v[509];
//		v[3065] = v[204] * v[2984] - 2e0*v[2998] * v[509];
//		v[3067] = -v[3066] + v[2987] * v[475] + v[2992] * v[505];
//		v[3068] = v[204] * v[2987] - 2e0*v[2998] * v[505];
//		v[3177] = v[282] * v[3068];
//		v[3069] = -v[3066] - v[2986] * v[474] - v[2993] * v[501];
//		v[3070] = v[204] * v[2986] - 2e0*v[2998] * v[501];
//		v[3168] = v[290] * v[3070];
//		v[3071] = v[2996] - (v[2991] * v[476]) / 2e0 - (v[2990] * v[496]) / 2e0;
//		v[3072] = v[2996] - (v[2945] * v[476]) / 2e0 - (v[2990] * v[477]) / 2e0;
//		v[3073] = v[3066] - v[2981] * v[476] - v[2990] * v[492];
//		v[3074] = v[204] * v[2981] - 2e0*v[2998] * v[492];
//		v[3075] = v[3066] + v[2989] * v[475] + v[2992] * v[488];
//		v[3076] = v[204] * v[2989] - 2e0*v[2998] * v[488];
//		v[3161] = v[290] * v[3076];
//		v[3077] = -v[2997] + v[2984] * v[475] + v[2992] * v[509];
//		v[3078] = -v[2997] - v[2987] * v[474] - v[2993] * v[505];
//		v[3079] = -v[2997] + v[2986] * v[475] + v[2992] * v[501];
//		v[3080] = -v[2997] - v[2989] * v[474] - v[2993] * v[488];
//		v[3081] = v[3091] + 2e0*v[2997] * v[487];
//		v[3082] = v[3057] * v[550] + v[3058] * v[553] + v[3059] * v[556] - v[3060] * v[851] - v[3077] * v[853] - v[3067] * v[855];
//		v[3083] = v[3057] * v[549] + v[3058] * v[552] + v[3059] * v[555] - v[3062] * v[851] - v[3064] * v[853] - v[3078] * v[855];
//		v[3084] = -v[3066] - v[2983] * v[476] - v[2990] * v[482];
//		v[3085] = v[204] * v[2983] - 2e0*v[2998] * v[482];
//		v[3086] = -v[2999] + v[2981] * v[475] + v[2992] * v[492];
//		v[3087] = -v[2999] - v[2987] * v[476] - v[2990] * v[505];
//		v[3088] = -v[2999] - v[2989] * v[476] - v[2990] * v[488];
//		v[3089] = -v[2999] + v[2983] * v[475] + v[2992] * v[482];
//		v[3090] = v[2997] * v[484] + v[2999] * v[487];
//		v[3092] = v[3091] + 2e0*v[2999] * v[484];
//		v[3191] = v[2714] * v[3092];
//		v[3093] = v[3057] * v[535] + v[3058] * v[538] + v[3059] * v[541] - v[3079] * v[851] - v[3004] * v[853] - v[3086] * v[855];
//		v[3094] = v[3000] - v[2984] * v[476] - v[2990] * v[509];
//		v[3095] = v[3000] - v[2986] * v[476] - v[2990] * v[501];
//		v[3096] = v[3000] - v[2981] * v[474] - v[2993] * v[492];
//		v[3097] = v[3000] - v[2983] * v[474] - v[2993] * v[482];
//		v[3098] = -(v[2999] * v[481]) - v[3000] * v[484];
//		v[3099] = -(v[2997] * v[481]) - v[3000] * v[487];
//		v[3100] = v[3091] + 2e0*v[3000] * v[481];
//		v[3195] = v[2713] * v[3100];
//		v[3101] = v[3057] * v[551] + v[3058] * v[554] + v[3059] * v[557] - v[3002] * v[851] - v[3094] * v[853] - v[3087] * v[855];
//		v[3102] = v[3057] * v[536] + v[3058] * v[539] + v[3059] * v[542] - v[3095] * v[851] - v[3071] * v[853] - v[3073] * v[855];
//		v[3103] = v[3057] * v[534] + v[3058] * v[537] + v[3059] * v[540] - v[3069] * v[851] - v[3063] * v[853] - v[3096] * v[855];
//		v[3104] = v[3057] * v[519] + v[3058] * v[522] + v[3059] * v[525] - v[3080] * v[851] - v[3097] * v[853] - v[3006] * v[855];
//		v[3105] = v[3057] * v[521] + v[3058] * v[524] + v[3059] * v[527] - v[3088] * v[851] - v[3084] * v[853] - v[3072] * v[855];
//		v[3106] = v[3057] * v[520] + v[3058] * v[523] + v[3059] * v[526] - v[3075] * v[851] - v[3089] * v[853] - v[3061] * v[855];
//		v[3107] = v[2966] + v[3040];
//		v[3140] = v[2717] * v[3107];
//		v[3108] = -v[2966] + v[3040];
//		v[3109] = (v[251] * v[3015] + v[162] * v[3107] + v[1825] * v[4070]) / v[253];
//		v[3110] = v[2964] + v[3048];
//		v[3148] = v[2717] * v[3110];
//		v[3111] = -v[2964] + v[3048];
//		v[3144] = v[2716] * v[3111];
//		v[3112] = (v[252] * v[3020] + v[158] * v[3110] + v[1829] * v[4070]) / v[253];
//		v[3113] = v[262] * v[3035] + v[3115];
//		v[3114] = v[3113] + v[3116] + v[4117];
//		v[3117] = v[3115] + v[3116];
//		v[3118] = v[256] * v[3024] + v[3120];
//		v[3119] = v[3118] + v[3121] + v[4116];
//		v[3122] = v[3120] + v[3121];
//		v[3123] = v[2743] * v[2977] + v[2739] * v[2978] + v[2736] * v[3013] + v[2740] * v[3014] + v[2742] * v[3015]
//			+ v[2737] * v[3021] + v[2741] * v[3027] + v[2734] * v[3034] + v[2744] * v[3035] + v[2733] * v[3039] + v[2738] * v[3044]
//			+ v[2735] * v[3047];
//		v[3124] = v[2744] * v[2979] + v[2735] * v[2980] + v[2733] * v[3011] + v[2741] * v[3017] + v[2742] * v[3018]
//			+ v[2734] * v[3022] + v[2737] * v[3023] + v[2743] * v[3024] + v[2740] * v[3028] + v[2739] * v[3036] + v[2738] * v[3037]
//			+ v[2736] * v[3046];
//		v[3125] = v[3126] + v[3128];
//		v[3127] = v[262] * v[3015] + v[3126];
//		v[3129] = v[3127] + v[3128] + v[4121];
//		v[3130] = v[2742] * v[2975] + v[2738] * v[2976] + v[2741] * v[3010] + v[2740] * v[3012] + v[2736] * v[3019]
//			+ v[2743] * v[3020] + v[2733] * v[3025] + v[2744] * v[3026] + v[2739] * v[3029] + v[2735] * v[3030] + v[2734] * v[3038]
//			+ v[2737] * v[3045];
//		v[3131] = (v[153] * v[3111] + v[3026] * v[3723] + v[1831] * v[4070]) / v[253];
//		v[3132] = v[3142] + v[3144];
//		v[4152] = v[3132] / v[253];
//		v[3133] = v[2967] + v[3049];
//		v[3138] = v[2716] * v[3133];
//		v[3134] = -v[2967] + v[3049];
//		v[3135] = v[3146] + v[3148];
//		v[4148] = v[3135] / v[253];
//		v[3136] = v[2715] * v[3031] + v[3138] + v[3140];
//		v[3139] = v[3136] - v[3140];
//		v[3141] = v[3136] - v[3138];
//		v[4149] = v[3141] / v[253];
//		v[3143] = v[2715] * v[3108] + v[3142];
//		v[3145] = v[3143] + v[3144];
//		v[3147] = v[2715] * v[3134] + v[3146];
//		v[3149] = v[3147] + v[3148];
//		v[3150] = v[2992] + v[3090];
//		v[3189] = v[2714] * v[3150];
//		v[3151] = -v[2992] + v[3090];
//		v[3152] = (v[277] * v[3065] + v[218] * v[3150] + v[1853] * v[4071]) / v[279];
//		v[3153] = v[2990] + v[3098];
//		v[3197] = v[2714] * v[3153];
//		v[3154] = -v[2990] + v[3098];
//		v[3193] = v[2713] * v[3154];
//		v[3155] = (v[278] * v[3070] + v[214] * v[3153] + v[1857] * v[4071]) / v[279];
//		v[3156] = v[288] * v[3085] + v[3160];
//		v[3158] = v[167] * v[3007] + v[168] * v[3008] + v[169] * v[3009] - v[223] * v[3057] - v[224] * v[3058] - v[225] * v[3059]
//			+ v[3157] + v[3076] * v[851] + v[3085] * v[853] + v[3005] * v[855] + v[3026] * v[857] + v[3035] * v[859]
//			+ v[2979] * v[861];
//		v[3159] = v[3156] + v[3161] + v[4108];
//		v[3162] = v[3160] + v[3161];
//		v[3164] = v[170] * v[3007] + v[171] * v[3008] + v[172] * v[3009] - v[226] * v[3057] - v[227] * v[3058] - v[228] * v[3059]
//			+ v[3163] + v[3070] * v[851] + v[3003] * v[853] + v[3074] * v[855] + v[3020] * v[857] + v[2977] * v[859]
//			+ v[3024] * v[861];
//		v[3165] = v[282] * v[3074] + v[3167];
//		v[3166] = v[3165] + v[3168] + v[4107];
//		v[3169] = v[3167] + v[3168];
//		v[3170] = v[2743] * v[3003] - v[2728] * v[3004] - v[2723] * v[3063] - v[2731] * v[3064] + v[2742] * v[3065]
//			- v[2726] * v[3071] - v[2732] * v[3077] - v[2721] * v[3084] + v[2744] * v[3085] - v[2718] * v[3089] - v[2727] * v[3094]
//			- v[2722] * v[3097];
//		v[3171] = v[2744] * v[3005] - v[2722] * v[3006] - v[2718] * v[3061] - v[2732] * v[3067] + v[2742] * v[3068]
//			- v[2721] * v[3072] - v[2726] * v[3073] + v[2743] * v[3074] - v[2731] * v[3078] - v[2728] * v[3086] - v[2727] * v[3087]
//			- v[2723] * v[3096];
//		v[3173] = v[173] * v[3007] + v[174] * v[3008] + v[175] * v[3009] - v[229] * v[3057] - v[230] * v[3058] - v[231] * v[3059]
//			+ v[3172] + v[3001] * v[851] + v[3065] * v[853] + v[3068] * v[855] + v[2975] * v[857] + v[3015] * v[859]
//			+ v[3018] * v[861];
//		v[4072] = v[302] * v[3158] + v[303] * v[3164] + v[304] * v[3173];
//		v[3374] = v[4072] / v[727];
//		v[3174] = v[3175] + v[3177];
//		v[3176] = v[288] * v[3065] + v[3175];
//		v[3178] = v[3176] + v[3177] + v[4112];
//		v[3179] = v[2742] * v[3001] - v[2727] * v[3002] - v[2732] * v[3060] - v[2731] * v[3062] - v[2723] * v[3069]
//			+ v[2743] * v[3070] - v[2718] * v[3075] + v[2744] * v[3076] - v[2728] * v[3079] - v[2722] * v[3080] - v[2721] * v[3088]
//			- v[2726] * v[3095];
//		v[3180] = (v[209] * v[3154] + v[3076] * v[3731] + v[1859] * v[4071]) / v[279];
//		v[3181] = v[3191] + v[3193];
//		v[4164] = v[3181] / v[279];
//		v[3182] = v[2993] + v[3099];
//		v[3187] = v[2713] * v[3182];
//		v[3183] = -v[2993] + v[3099];
//		v[3184] = v[3195] + v[3197];
//		v[4160] = v[3184] / v[279];
//		v[3185] = v[2712] * v[3081] + v[3187] + v[3189];
//		v[3188] = v[3185] - v[3189];
//		v[3190] = v[3185] - v[3187];
//		v[4161] = v[3190] / v[279];
//		v[3192] = v[2712] * v[3151] + v[3191];
//		v[3194] = v[3192] + v[3193];
//		v[3196] = v[2712] * v[3183] + v[3195];
//		v[3198] = v[3196] + v[3197];
//		v[3199] = v[3158];
//		v[3288] = v[3199];
//		v[3200] = v[3164];
//		v[3287] = v[3200];
//		v[3201] = v[3173];
//		v[3286] = v[3201];
//		v[3202] = -(v[1888] * v[2709] * v[4072]);
//		v[3203] = v[3374];
//		v[4073] = v[315] * v[3203];
//		v[3268] = v[3158] * v[316] + v[302] * v[4073];
//		v[3217] = v[316] * v[3173] + v[304] * v[4073];
//		v[4074] = v[2698] * v[3217];
//		v[3248] = v[301] * v[4074];
//		v[3245] = v[291] * v[4074];
//		v[3242] = v[281] * v[4074];
//		v[3239] = v[275] * v[4074];
//		v[3236] = v[265] * v[4074];
//		v[3233] = v[255] * v[4074];
//		v[3206] = v[316] * v[3164] + v[303] * v[4073];
//		v[4075] = v[2697] * v[3206];
//		v[3247] = v[301] * v[4075];
//		v[3244] = v[291] * v[4075];
//		v[3241] = v[281] * v[4075];
//		v[3238] = v[275] * v[4075];
//		v[3235] = v[265] * v[4075];
//		v[3232] = v[255] * v[4075];
//		v[3204] = v[3268];
//		v[4076] = v[2696] * v[3204];
//		v[3227] = v[301] * v[4076];
//		v[3225] = v[291] * v[4076];
//		v[3223] = v[281] * v[4076];
//		v[3221] = v[275] * v[4076];
//		v[3219] = v[265] * v[4076];
//		v[3216] = v[255] * v[4076];
//		v[3205] = v[4075] + v[4076];
//		v[3207] = v[3216] + v[3232];
//		v[3208] = v[3219] + v[3235];
//		v[3209] = v[3221] + v[3238];
//		v[3210] = v[3223] + v[3241];
//		v[3211] = v[3225] + v[3244];
//		v[3212] = v[3227] + v[3247];
//		v[3215] = v[4074] + v[4076];
//		v[3218] = v[3216] + v[3233];
//		v[3220] = v[3219] + v[3236];
//		v[3222] = v[3221] + v[3239];
//		v[3224] = v[3223] + v[3242];
//		v[3226] = v[3225] + v[3245];
//		v[3228] = v[3227] + v[3248];
//		v[3231] = v[4074] + v[4075];
//		v[3234] = v[3232] + v[3233];
//		v[3237] = v[3235] + v[3236];
//		v[3240] = v[3238] + v[3239];
//		v[3243] = v[3241] + v[3242];
//		v[3246] = v[3244] + v[3245];
//		v[3249] = v[3247] + v[3248];
//		v[3252] = v[2011] * v[3041] + v[2012] * v[3042] + v[2013] * v[3107] + v[262] * v[3109] + v[2014] * v[3110]
//			+ v[264] * v[3112] + v[2682] * v[3114] + v[3118] * v[3724] + v[3125] * v[3727] + (*a4)*v[6245 + i2691];
//		v[3254] = v[2016] * v[3041] + v[2017] * v[3050] + v[256] * v[3109] + v[2018] * v[3111] + v[1572] * v[3113]
//			+ v[2683] * v[3119] + v[1571] * v[3127] + v[264] * v[3131] + v[2019] * v[3133] + (*a4)*v[6257 + i2691];
//		v[3256] = v[2022] * v[3031] + v[2021] * v[3041] + v[2023] * v[3108] + v[256] * v[3112] + v[1574] * v[3117]
//			+ v[1573] * v[3122] + v[2684] * v[3129] + v[262] * v[3131] + v[2024] * v[3134] + (*a4)*v[6269 + i2691];
//		v[3258] = v[2026] * v[3091] + v[2027] * v[3092] + v[2028] * v[3150] + v[288] * v[3152] + v[2029] * v[3153]
//			+ v[290] * v[3155] + v[2685] * v[3159] + v[3165] * v[3732] + v[3174] * v[3735] + (*a4)*v[6281 + i2691];
//		v[3260] = v[2031] * v[3091] + v[2032] * v[3100] + v[282] * v[3152] + v[2033] * v[3154] + v[1578] * v[3156]
//			+ v[2686] * v[3166] + v[1577] * v[3176] + v[290] * v[3180] + v[2034] * v[3182] + (*a4)*v[6293 + i2691];
//		v[3262] = v[2037] * v[3081] + v[2036] * v[3091] + v[2038] * v[3151] + v[282] * v[3155] + v[1580] * v[3162]
//			+ v[1579] * v[3169] + v[2687] * v[3178] + v[288] * v[3180] + v[2039] * v[3183] + (*a4)*v[6305 + i2691];
//		b3263 = b728;
//		if (b3263) {
//			b3264 = b730;
//			if (b3264) {
//				v[3203] = 0e0;
//				v[3204] = 0e0;
//			}
//			else {
//			};
//		}
//		else {
//		};
//		v[3265] = -(v[255] * v[3054]) - v[275] * v[3055] - v[265] * v[3056] + v[281] * v[3104] + v[301] * v[3105]
//			+ v[291] * v[3106] - (*a4)*v[3157];
//		v[4077] = -(v[317] * v[3265]) + v[2100] * v[3268];
//		v[3266] = -(v[265] * v[3043]) - v[275] * v[3052] - v[255] * v[3053] + v[291] * v[3093] + v[301] * v[3102]
//			+ v[281] * v[3103] - (*a4)*v[3163];
//		v[4078] = v[2095] * v[3206] - v[318] * v[3266];
//		v[3267] = -(v[265] * v[3032]) - v[255] * v[3033] - v[275] * v[3051] + v[291] * v[3082] + v[281] * v[3083]
//			+ v[301] * v[3101] - (*a4)*v[3172];
//		v[4079] = v[2090] * v[3217] - v[319] * v[3267];
//		v[3270] = 2e0*v[2095] * v[4075] + v[2697] * (v[4077] + v[4079]) + v[3266] * v[4086];
//		v[3300] = v[3270];
//		v[3272] = 2e0*v[2090] * v[4074] + v[2698] * (v[4077] + v[4078]) + v[3267] * v[4087];
//		v[3299] = v[3272];
//		v[3273] = 2e0*v[2100] * v[4076] + v[2696] * (v[4078] + v[4079]) + v[3265] * v[4088];
//		v[3301] = v[3273];
//		v[3274] = 0e0;
//		v[3275] = 0e0;
//		v[3276] = 0e0;
//		v[3277] = 0e0;
//		v[3278] = 0e0;
//		v[3279] = 0e0;
//		v[3280] = 0e0;
//		b3281 = b728;
//		if (b3281) {
//			v[3282] = 0e0;
//			v[3283] = 0e0;
//			v[3284] = 0e0;
//			b3285 = b747;
//			if (b3285) {
//				v[3284] = v[3201];
//				v[3201] = 0e0;
//				v[3283] = v[3200];
//				v[3200] = 0e0;
//				v[3282] = v[3199];
//				v[3199] = 0e0;
//			}
//			else {
//				v[3279] = -v[3286];
//				v[3201] = 0e0;
//				v[3278] = -v[3287];
//				v[3200] = 0e0;
//				v[3277] = -v[3288];
//				v[3199] = 0e0;
//			};
//			v[4080] = (v[3282] * v[714] + v[3283] * v[716] + v[3284] * v[718])*(*zetan);
//			b3289 = b730;
//			if (b3289) {
//				v[3276] = v[3284] * v[742];
//				v[3275] = v[3283] * v[742];
//				v[3274] = v[3282] * v[742];
//				v[3280] = (v[1201] * v[4080] * v[4411] * v[740]) / v[1202];
//			}
//			else {
//				v[3276] = v[3284] * v[746];
//				v[3275] = v[3283] * v[746];
//				v[3274] = v[3282] * v[746];
//				v[3202] = v[3202] + (v[1206] * v[4080] * v[4412] * v[745]) / v[1207];
//			};
//		}
//		else {
//		};
//		v[4081] = v[317] * v[3274];
//		v[4082] = v[318] * v[3275];
//		v[4083] = v[319] * v[3276];
//		v[3712] = v[319] * (v[3205] + v[4081] + v[4082]) + v[3217] * v[4413] + v[3276] * v[717];
//		v[3709] = v[318] * (v[3215] + v[4081] + v[4083]) + v[3206] * v[4414] + v[3275] * v[715];
//		v[3706] = v[317] * (v[3231] + v[4082] + v[4083]) + v[3268] * v[4415] + v[3274] * v[713];
//		v[3302] = v[3202];
//		v[3298] = v[3277];
//		v[3297] = v[3278];
//		v[3296] = v[3279];
//		b3294 = b728;
//		if (b3294) {
//			v[4085] = -(v[319] * v[3296]) - v[318] * v[3297] - v[317] * v[3298];
//			b3295 = b730;
//			if (b3295) {
//				v[3272] = v[3272] + v[3279] * v[4084];
//				v[3279] = 0e0;
//				v[3270] = v[3270] + v[3278] * v[4084];
//				v[3278] = 0e0;
//				v[3273] = v[3273] + v[3277] * v[4084];
//				v[3277] = 0e0;
//				v[3280] = v[3280] + v[3767] * v[4085] * v[4419];
//				v[3202] = v[3202] - v[3280];
//			}
//			else {
//				v[3272] = v[3299] + v[3296] * v[736];
//				v[3279] = 0e0;
//				v[3270] = v[3300] + v[3297] * v[736];
//				v[3278] = 0e0;
//				v[3273] = v[3301] + v[3298] * v[736];
//				v[3277] = 0e0;
//				v[3202] = v[3302] + (*n2)*v[4085] * v[4423] * v[724];
//			};
//		}
//		else {
//		};
//		v[3303] = v[2698] * v[3262] + v[301] * v[3276];
//		v[4105] = v[319] * v[3303];
//		v[3304] = v[2698] * v[3260] + v[291] * v[3276];
//		v[4103] = v[319] * v[3304];
//		v[3305] = v[2698] * v[3258] + v[281] * v[3276];
//		v[4101] = v[319] * v[3305];
//		v[3306] = v[2698] * v[3256] + v[275] * v[3276];
//		v[4099] = v[319] * v[3306];
//		v[3307] = v[2698] * v[3254] + v[265] * v[3276];
//		v[4097] = v[319] * v[3307];
//		v[3308] = v[2698] * v[3252] + v[255] * v[3276];
//		v[4095] = v[319] * v[3308];
//		v[3309] = v[2697] * v[3262] + v[301] * v[3275];
//		v[4106] = v[318] * v[3309];
//		v[3310] = v[2697] * v[3260] + v[291] * v[3275];
//		v[4104] = v[318] * v[3310];
//		v[3311] = v[2697] * v[3258] + v[281] * v[3275];
//		v[4102] = v[318] * v[3311];
//		v[3312] = v[2697] * v[3256] + v[275] * v[3275];
//		v[4100] = v[318] * v[3312];
//		v[3313] = v[2697] * v[3254] + v[265] * v[3275];
//		v[4098] = v[318] * v[3313];
//		v[3314] = v[2697] * v[3252] + v[255] * v[3275];
//		v[4096] = v[318] * v[3314];
//		v[3315] = v[2696] * v[3262] + v[301] * v[3274];
//		v[4094] = v[317] * v[3315];
//		v[3316] = v[2696] * v[3260] + v[291] * v[3274];
//		v[4093] = v[317] * v[3316];
//		v[3317] = v[2696] * v[3258] + v[281] * v[3274];
//		v[4092] = v[317] * v[3317];
//		v[3318] = v[2696] * v[3256] + v[275] * v[3274];
//		v[4091] = v[317] * v[3318];
//		v[3319] = v[2696] * v[3254] + v[265] * v[3274];
//		v[4090] = v[317] * v[3319];
//		v[3320] = v[2696] * v[3252] + v[255] * v[3274];
//		v[4089] = v[317] * v[3320];
//		v[3321] = -(v[2925] * v[3101]) - v[2924] * v[3102] - v[2923] * v[3105] + v[1241] * v[3274] + v[1242] * v[3275]
//			+ v[1243] * v[3276] + v[3206] * v[4424] + v[3217] * v[4425] + v[3268] * v[4426];
//		v[3465] = v[3321] * v[3720];
//		v[3322] = -(v[2925] * v[3082]) - v[2924] * v[3093] - v[2923] * v[3106] + v[1245] * v[3274] + v[1246] * v[3275]
//			+ v[1247] * v[3276] + v[3206] * v[4427] + v[3217] * v[4428] + v[3268] * v[4429];
//		v[4114] = v[279] * v[3322];
//		v[3469] = v[3322] * v[3719];
//		v[3323] = -(v[2925] * v[3083]) - v[2924] * v[3103] - v[2923] * v[3104] + v[1249] * v[3274] + v[1250] * v[3275]
//			+ v[1251] * v[3276] + v[3206] * v[4430] + v[3217] * v[4431] + v[3268] * v[4432];
//		v[4115] = v[279] * v[3323];
//		v[3472] = v[3323] * v[3718];
//		v[3324] = v[2925] * v[3051] + v[2924] * v[3052] + v[2923] * v[3055] + v[1253] * v[3274] + v[1254] * v[3275]
//			+ v[1255] * v[3276] + v[3206] * v[4433] + v[3217] * v[4434] + v[3268] * v[4435];
//		v[3556] = v[3324] * v[3717];
//		v[3325] = v[2925] * v[3032] + v[2924] * v[3043] + v[2923] * v[3056] + v[1257] * v[3274] + v[1258] * v[3275]
//			+ v[1259] * v[3276] + v[3206] * v[4436] + v[3217] * v[4437] + v[3268] * v[4438];
//		v[4123] = v[253] * v[3325];
//		v[3560] = v[3325] * v[3716];
//		v[3326] = v[2925] * v[3033] + v[2924] * v[3053] + v[2923] * v[3054] + v[1261] * v[3274] + v[1262] * v[3275]
//			+ v[1263] * v[3276] + v[3206] * v[4439] + v[3217] * v[4440] + v[3268] * v[4441];
//		v[7166] = v[3706];
//		v[7167] = v[3709];
//		v[7168] = v[3712];
//		v[7169] = v[2716] * v[3109] + v[2715] * v[3112] + v[2682] * v[3326] + v[3041] * (v[2638] * v[2715] + v[2640] * v[2716]
//			+ v[3566]) + v[3325] * v[3724] + v[3324] * v[3727] + v[3018] * v[4043] + v[3024] * v[4044] + v[3145] * v[4119]
//			+ v[2979] * v[4147] + v[155] * v[4148] + v[160] * v[4149];
//		v[7170] = v[2717] * v[3109] + v[2715] * v[3131] + v[1571] * v[3324] + v[2683] * v[3325] + v[1572] * v[3326] + v[3041] *
//			(v[2637] * v[2715] + v[2641] * v[2717] + v[3559]) + v[3035] * v[4045] + v[3149] * v[4118] + v[2977] * v[4150]
//			+ v[3015] * v[4151] + v[152] * v[4152] + v[3139] * v[4153];
//		v[7171] = v[2717] * v[3112] + v[2716] * v[3131] + v[2684] * v[3324] + v[1573] * v[3325] + v[1574] * v[3326] + v[3041] *
//			(v[2639] * v[2716] + v[2642] * v[2717] + v[3551]) + v[3136] * v[4122] + v[2975] * v[4154] + v[3020] * v[4155]
//			+ v[3026] * v[4156] + v[3143] * v[4157] + v[3147] * v[4158];
//		v[7172] = -v[3706];
//		v[7173] = -v[3709];
//		v[7174] = -v[3712];
//		v[7175] = v[2713] * v[3152] + v[2712] * v[3155] + v[2685] * v[3323] + v[3091] * (v[2609] * v[2712] + v[2611] * v[2713]
//			+ v[3475]) + v[3322] * v[3732] + v[3321] * v[3735] + v[3068] * v[4037] + v[3074] * v[4038] + v[3194] * v[4110]
//			+ v[3005] * v[4159] + v[211] * v[4160] + v[216] * v[4161];
//		v[7176] = v[2714] * v[3152] + v[2712] * v[3180] + v[1577] * v[3321] + v[2686] * v[3322] + v[1578] * v[3323] + v[3091] *
//			(v[2608] * v[2712] + v[2612] * v[2714] + v[3468]) + v[3085] * v[4039] + v[3198] * v[4109] + v[3003] * v[4162]
//			+ v[3065] * v[4163] + v[208] * v[4164] + v[3188] * v[4165];
//		v[7177] = v[2714] * v[3155] + v[2713] * v[3180] + v[2687] * v[3321] + v[1579] * v[3322] + v[1580] * v[3323] + v[3091] *
//			(v[2610] * v[2713] + v[2613] * v[2714] + v[3460]) + v[3185] * v[4113] + v[3001] * v[4166] + v[3070] * v[4167]
//			+ v[3076] * v[4168] + v[3192] * v[4169] + v[3196] * v[4170];
//		v[4124] = v[253] * v[3326];
//		v[3563] = v[3326] * v[3715];
//		v[3328] = -(v[319] * (v[3212] + v[4094] + v[4106])) + v[3217] * v[4442] - v[3303] * v[717];
//		v[3329] = -(v[319] * (v[3211] + v[4093] + v[4104])) + v[3217] * v[4443] - v[3304] * v[717];
//		v[3330] = -(v[319] * (v[3210] + v[4092] + v[4102])) + v[3217] * v[4444] - v[3305] * v[717];
//		v[3331] = v[319] * (v[3209] + v[4091] + v[4100]) + v[3217] * v[4445] + v[3306] * v[717];
//		v[3332] = v[319] * (v[3208] + v[4090] + v[4098]) + v[3217] * v[4446] + v[3307] * v[717];
//		v[3272] = v[3272] + v[2265] * v[3276] + v[2336] * v[3303] + v[2340] * v[3304] + v[2344] * v[3305] + v[2348] * v[3306]
//			+ v[2352] * v[3307] + v[2356] * v[3308] + v[3274] * v[3333] + v[3275] * v[3334] + v[3205] * v[3335] + v[3207] * v[471]
//			+ v[3208] * v[472] + v[3209] * v[473] - v[3210] * v[570] - v[3211] * v[571] - v[3212] * v[572] + 2e0*v[319] * (-
//			(v[2698] * v[3267]) + v[3276] * v[3335] + v[3308] * v[471] + v[3307] * v[472] + v[3306] * v[473] - v[3305] * v[570]
//				- v[3304] * v[571] - v[3303] * v[572]) + v[318] * (v[3314] * v[471] + v[3313] * v[472] + v[3312] * v[473]
//					- v[3311] * v[570] - v[3310] * v[571] - v[3309] * v[572]) + v[317] * (v[3320] * v[471] + v[3319] * v[472]
//						+ v[3318] * v[473] - v[3317] * v[570] - v[3316] * v[571] - v[3315] * v[572]);
//		v[3336] = v[319] * (v[3207] + v[4089] + v[4096]) + v[3217] * v[4447] + v[3308] * v[717];
//		v[3339] = v[318] * (v[3218] + v[4089] + v[4095]) + v[3206] * v[4448] + v[3314] * v[715];
//		v[3342] = v[318] * (v[3220] + v[4090] + v[4097]) + v[3206] * v[4449] + v[3313] * v[715];
//		v[3345] = v[318] * (v[3222] + v[4091] + v[4099]) + v[3206] * v[4450] + v[3312] * v[715];
//		v[3348] = -(v[318] * (v[3224] + v[4092] + v[4101])) + v[3206] * v[4451] - v[3311] * v[715];
//		v[3351] = -(v[318] * (v[3226] + v[4093] + v[4103])) + v[3206] * v[4452] - v[3310] * v[715];
//		v[3270] = v[3270] + v[2268] * v[3275] + v[3274] * v[3352] + v[3276] * v[3353] + v[3215] * v[3354] + v[3309] * v[3746]
//			+ v[3310] * v[3748] + v[3311] * v[3750] + v[3312] * v[3752] + v[3313] * v[3754] + v[3314] * v[3756] + v[3218] * v[468]
//			+ v[3220] * v[469] + v[3222] * v[470] - v[3224] * v[567] - v[3226] * v[568] - v[3228] * v[569] + v[319] * (v[3308] * v[468]
//				+ v[3307] * v[469] + v[3306] * v[470] - v[3305] * v[567] - v[3304] * v[568] - v[3303] * v[569]) + 2e0*v[318] * (-
//				(v[2697] * v[3266]) + v[3275] * v[3354] + v[3314] * v[468] + v[3313] * v[469] + v[3312] * v[470] - v[3311] * v[567]
//					- v[3310] * v[568] - v[3309] * v[569]) + v[317] * (v[3320] * v[468] + v[3319] * v[469] + v[3318] * v[470]
//						- v[3317] * v[567] - v[3316] * v[568] - v[3315] * v[569]);
//		v[3357] = -(v[318] * (v[3228] + v[4094] + v[4105])) + v[3206] * v[4453] - v[3309] * v[715];
//		v[3359] = v[317] * (v[3234] + v[4095] + v[4096]) + v[3268] * v[4454] + v[3320] * v[713];
//		v[3361] = v[317] * (v[3237] + v[4097] + v[4098]) + v[3268] * v[4455] + v[3319] * v[713];
//		v[3363] = v[317] * (v[3240] + v[4099] + v[4100]) + v[3268] * v[4456] + v[3318] * v[713];
//		v[3365] = -(v[317] * (v[3243] + v[4101] + v[4102])) + v[3268] * v[4457] - v[3317] * v[713];
//		v[3367] = -(v[317] * (v[3246] + v[4103] + v[4104])) + v[3268] * v[4458] - v[3316] * v[713];
//		v[3273] = v[3273] + v[2272] * v[3274] + v[3275] * v[3368] + v[3276] * v[3369] + v[3231] * v[3370] - v[3315] * v[3747]
//			- v[3316] * v[3749] - v[3317] * v[3751] + v[3318] * v[3753] + v[3319] * v[3755] + v[3320] * v[3757] + v[3234] * v[465]
//			+ v[3237] * v[466] + v[3240] * v[467] - v[3243] * v[564] - v[3246] * v[565] - v[3249] * v[566] + v[319] * (v[3308] * v[465]
//				+ v[3307] * v[466] + v[3306] * v[467] - v[3305] * v[564] - v[3304] * v[565] - v[3303] * v[566]) + v[318] *
//				(v[3314] * v[465] + v[3313] * v[466] + v[3312] * v[467] - v[3311] * v[564] - v[3310] * v[565] - v[3309] * v[566])
//			+ 2e0*v[317] * (-(v[2696] * v[3265]) + v[3274] * v[3370] + v[3320] * v[465] + v[3319] * v[466] + v[3318] * v[467]
//				- v[3317] * v[564] - v[3316] * v[565] - v[3315] * v[566]);
//		v[3372] = -(v[317] * (v[3249] + v[4105] + v[4106])) + v[3268] * v[4459] - v[3315] * v[713];
//		v[3202] = v[3202] + v[315] * (v[2706] * v[3158] + v[2707] * v[3164] + v[2708] * v[3173] + v[303] * v[3270]
//			+ v[304] * v[3272] + v[302] * v[3273]) + v[314] * v[3373] * v[3374];
//		v[3376] = v[316] * v[3272] + v[2708] * v[4073] + (v[2709] * v[3173] + v[304] * v[3202]) / v[727];
//		v[3378] = v[316] * v[3270] + v[2707] * v[4073] + (v[2709] * v[3164] + v[303] * v[3202]) / v[727];
//		v[3380] = v[316] * v[3273] + v[2706] * v[4073] + (v[2709] * v[3158] + v[302] * v[3202]) / v[727];
//		v[3381] = -(v[2742] * v[3059]) - v[200] * v[3376];
//		v[3382] = -(v[2742] * v[3058]) - v[199] * v[3376];
//		v[3383] = -(v[2742] * v[3057]) - v[198] * v[3376];
//		v[3384] = v[2742] * v[3009] + v[144] * v[3376];
//		v[3385] = v[2742] * v[3008] + v[143] * v[3376];
//		v[3386] = v[2742] * v[3007] + v[142] * v[3376];
//		v[3387] = -(v[2743] * v[3059]) - v[200] * v[3378];
//		v[3388] = -(v[2743] * v[3058]) - v[199] * v[3378];
//		v[3389] = -(v[2743] * v[3057]) - v[198] * v[3378];
//		v[3390] = v[2743] * v[3009] + v[144] * v[3378];
//		v[3391] = v[2743] * v[3008] + v[143] * v[3378];
//		v[3392] = v[2743] * v[3007] + v[142] * v[3378];
//		v[3393] = -(v[2744] * v[3059]) - v[200] * v[3380];
//		v[3394] = -(v[2744] * v[3058]) - v[199] * v[3380];
//		v[3395] = -(v[2744] * v[3057]) - v[198] * v[3380];
//		v[3396] = v[2744] * v[3009] + v[144] * v[3380];
//		v[3397] = v[2744] * v[3008] + v[143] * v[3380];
//		v[3398] = v[2744] * v[3007] + v[142] * v[3380];
//		v[3399] = v[288] * v[3321] + v[2712] * v[4107];
//		v[3400] = v[282] * v[3321] + v[2712] * v[4108];
//		v[3402] = v[1378] * v[3321] + v[2712] * (v[3169] / v[279] + v[3091] * v[3401]) + v[3399] * v[4109];
//		v[3404] = v[1384] * v[3321] + v[2712] * (v[3162] / v[279] + v[3091] * v[3403]) + v[3400] * v[4110];
//		v[3406] = (v[218] * v[3399] + v[216] * v[3400] + v[2712] * (v[3178] + v[3405] * v[4071]) + v[3321] * v[4111]) / v[279];
//		v[3407] = v[290] * v[3322] + v[2713] * v[4112];
//		v[3408] = v[282] * v[3322] + v[2713] * v[4108];
//		v[3410] = v[1373] * v[3322] + v[2713] * (v[3176] / v[279] + v[3091] * v[3409]) + v[3407] * v[4113];
//		v[3411] = v[3399] + v[3407];
//		v[3412] = v[3402] + v[3410];
//		v[3414] = (v[2749] * v[3076] + v[207] * v[3408] + v[209] * v[3411] + v[3661] * v[4071] + v[2713] * (v[3156]
//			+ v[3413] * v[4071]) + v[1383] * v[4114]) / v[279];
//		v[3416] = (v[214] * v[3407] + v[211] * v[3408] + v[2713] * (v[3166] + v[3415] * v[4071]) + v[1377] * v[4114]) / v[279];
//		v[3417] = v[288] * v[3323] + v[2714] * v[4107];
//		v[3418] = v[290] * v[3323] + v[2714] * v[4112];
//		v[3419] = v[3408] + v[3417];
//		v[3421] = (v[208] * v[3417] + v[209] * v[3418] + v[2714] * (v[3159] + v[3420] * v[4071]) + v[1381] * v[4115]) / v[279];
//		v[3422] = v[3400] + v[3418];
//		v[3424] = (v[2754] * v[3070] + v[213] * v[3417] + v[214] * v[3422] + v[3660] * v[4071] + v[2714] * (v[3165]
//			+ v[3423] * v[4071]) + v[1388] * v[4115]) / v[279];
//		v[3425] = v[3414] + v[3424];
//		v[3427] = (v[2755] * v[3065] + v[219] * v[3418] + v[218] * v[3419] + v[3659] * v[4071] + v[2714] * (v[3174]
//			+ v[3426] * v[4071]) + v[1391] * v[4115]) / v[279];
//		v[3428] = v[3404] + v[3427];
//		v[3429] = v[262] * v[3324] + v[2715] * v[4116];
//		v[3430] = v[256] * v[3324] + v[2715] * v[4117];
//		v[3432] = v[1402] * v[3324] + v[2715] * (v[3122] / v[253] + v[3041] * v[3431]) + v[3429] * v[4118];
//		v[3434] = v[1408] * v[3324] + v[2715] * (v[3117] / v[253] + v[3041] * v[3433]) + v[3430] * v[4119];
//		v[3436] = (v[162] * v[3429] + v[160] * v[3430] + v[2715] * (v[3129] + v[3435] * v[4070]) + v[3324] * v[4120]) / v[253];
//		v[3437] = v[264] * v[3325] + v[2716] * v[4121];
//		v[3438] = v[256] * v[3325] + v[2716] * v[4117];
//		v[3440] = v[1397] * v[3325] + v[2716] * (v[3127] / v[253] + v[3041] * v[3439]) + v[3437] * v[4122];
//		v[3441] = v[3429] + v[3437];
//		v[3442] = v[3432] + v[3440];
//		v[3444] = (v[2782] * v[3026] + v[151] * v[3438] + v[153] * v[3441] + v[3686] * v[4070] + v[2716] * (v[3113]
//			+ v[3443] * v[4070]) + v[1407] * v[4123]) / v[253];
//		v[3446] = (v[158] * v[3437] + v[155] * v[3438] + v[2716] * (v[3119] + v[3445] * v[4070]) + v[1401] * v[4123]) / v[253];
//		v[3447] = v[262] * v[3326] + v[2717] * v[4116];
//		v[3448] = v[264] * v[3326] + v[2717] * v[4121];
//		v[3449] = v[3438] + v[3447];
//		v[3451] = (v[152] * v[3447] + v[153] * v[3448] + v[2717] * (v[3114] + v[3450] * v[4070]) + v[1405] * v[4124]) / v[253];
//		v[3452] = v[3430] + v[3448];
//		v[3454] = (v[2787] * v[3020] + v[157] * v[3447] + v[158] * v[3452] + v[3685] * v[4070] + v[2717] * (v[3118]
//			+ v[3453] * v[4070]) + v[1412] * v[4124]) / v[253];
//		v[3455] = v[3444] + v[3454];
//		v[3457] = (v[2788] * v[3015] + v[163] * v[3448] + v[162] * v[3449] + v[3684] * v[4070] + v[2717] * (v[3125]
//			+ v[3456] * v[4070]) + v[1415] * v[4124]) / v[253];
//		v[3458] = v[3434] + v[3457];
//		v[3459] = QBi[2][2] * v[3381] + QBi[2][1] * v[3382] + QBi[2][0] * v[3383] + v[1577] * v[3407] + v[292] * v[3465]
//			+ v[3418] * v[3735] + v[290] * (v[3185] / v[279] + v[3091] * v[4460]);
//		v[3462] = (v[2755] * v[3150]) / v[279] + QBi[1][2] * v[3381] + QBi[1][1] * v[3382] + QBi[1][0] * v[3383]
//			+ v[2687] * v[3399] + v[289] * v[3469] + v[3188] * v[3719] + v[3419] * v[3735] + v[3091] * v[4461];
//		v[3464] = QBi[0][2] * v[3381] + QBi[0][1] * v[3382] + QBi[0][0] * v[3383] + v[2687] * v[3400] + v[277] * v[3472]
//			+ v[282] * (v[4161] + v[3091] * v[4462]);
//		v[3466] = (v[2754] * v[3153]) / v[279] + QBi[2][2] * v[3387] + QBi[2][1] * v[3388] + QBi[2][0] * v[3389]
//			+ v[2686] * v[3407] + v[3196] * v[3720] + v[3422] * v[3732] + v[3465] * v[3736] + v[3091] * v[4463];
//		v[3470] = QBi[1][2] * v[3387] + QBi[1][1] * v[3388] + QBi[1][0] * v[3389] + v[1579] * v[3399] + v[283] * v[3469]
//			+ v[3417] * v[3732] + v[288] * (v[3198] / v[279] + v[3091] * v[4464]);
//		v[3473] = QBi[0][2] * v[3387] + QBi[0][1] * v[3388] + QBi[0][0] * v[3389] + v[2686] * v[3408] + v[278] * v[3472]
//			+ v[282] * (v[4160] + v[3091] * v[4465]);
//		v[3474] = (v[2749] * v[3154]) / v[279] + QBi[2][2] * v[3393] + QBi[2][1] * v[3394] + QBi[2][0] * v[3395]
//			+ v[1578] * v[3411] + v[2685] * v[3418] + v[3192] * v[3720] + v[3465] * v[3734] + v[3091] * v[4466];
//		v[3476] = QBi[1][2] * v[3393] + QBi[1][1] * v[3394] + QBi[1][0] * v[3395] + v[2685] * v[3417] + v[3469] * v[3731]
//			+ v[288] * (v[3091] * v[4127] + v[4164]);
//		v[3479] = QBi[0][2] * v[3393] + QBi[0][1] * v[3394] + QBi[0][0] * v[3395] + v[1580] * v[3400] + v[1578] * v[3408]
//			+ v[280] * v[3472] + v[282] * (v[3194] / v[279] + v[3091] * v[4467]);
//		v[3480] = v[2718] * v[3057] + v[198] * v[3367];
//		v[3481] = v[2718] * v[3058] + v[199] * v[3367];
//		v[3482] = v[2718] * v[3059] + v[200] * v[3367];
//		v[3483] = QBi[0][0] * v[3480] + QBi[0][1] * v[3481] + QBi[0][2] * v[3482];
//		v[3484] = QBi[1][0] * v[3480] + QBi[1][1] * v[3481] + QBi[1][2] * v[3482];
//		v[3485] = QBi[2][0] * v[3480] + QBi[2][1] * v[3481] + QBi[2][2] * v[3482];
//		v[3486] = v[2721] * v[3057] + v[198] * v[3372];
//		v[3487] = v[2721] * v[3058] + v[199] * v[3372];
//		v[3488] = v[2721] * v[3059] + v[200] * v[3372];
//		v[3489] = QBi[0][0] * v[3486] + QBi[0][1] * v[3487] + QBi[0][2] * v[3488];
//		v[3490] = QBi[2][0] * v[3486] + QBi[2][1] * v[3487] + QBi[2][2] * v[3488];
//		v[3491] = QBi[1][0] * v[3486] + QBi[1][1] * v[3487] + QBi[1][2] * v[3488];
//		v[3492] = v[2722] * v[3057] + v[198] * v[3365];
//		v[3493] = v[2722] * v[3058] + v[199] * v[3365];
//		v[3494] = v[2722] * v[3059] + v[200] * v[3365];
//		v[3495] = QBi[1][0] * v[3492] + QBi[1][1] * v[3493] + QBi[1][2] * v[3494];
//		v[3496] = QBi[2][0] * v[3492] + QBi[2][1] * v[3493] + QBi[2][2] * v[3494];
//		v[3497] = QBi[0][0] * v[3492] + QBi[0][1] * v[3493] + QBi[0][2] * v[3494];
//		v[3498] = v[2723] * v[3057] + v[198] * v[3348];
//		v[3499] = v[2723] * v[3058] + v[199] * v[3348];
//		v[3500] = v[2723] * v[3059] + v[200] * v[3348];
//		v[3501] = QBi[0][0] * v[3498] + QBi[0][1] * v[3499] + QBi[0][2] * v[3500];
//		v[3502] = QBi[1][0] * v[3498] + QBi[1][1] * v[3499] + QBi[1][2] * v[3500];
//		v[3503] = QBi[2][0] * v[3498] + QBi[2][1] * v[3499] + QBi[2][2] * v[3500];
//		v[3504] = v[2726] * v[3057] + v[198] * v[3357];
//		v[3505] = v[2726] * v[3058] + v[199] * v[3357];
//		v[3506] = v[2726] * v[3059] + v[200] * v[3357];
//		v[3507] = QBi[2][0] * v[3504] + QBi[2][1] * v[3505] + QBi[2][2] * v[3506];
//		v[3508] = QBi[0][0] * v[3504] + QBi[0][1] * v[3505] + QBi[0][2] * v[3506];
//		v[3509] = QBi[1][0] * v[3504] + QBi[1][1] * v[3505] + QBi[1][2] * v[3506];
//		v[3510] = v[2727] * v[3057] + v[198] * v[3328];
//		v[3511] = v[2727] * v[3058] + v[199] * v[3328];
//		v[3512] = v[2727] * v[3059] + v[200] * v[3328];
//		v[3513] = QBi[1][0] * v[3510] + QBi[1][1] * v[3511] + QBi[1][2] * v[3512];
//		v[3514] = QBi[0][0] * v[3510] + QBi[0][1] * v[3511] + QBi[0][2] * v[3512];
//		v[3515] = QBi[2][0] * v[3510] + QBi[2][1] * v[3511] + QBi[2][2] * v[3512];
//		v[3516] = -(v[2751] * v[2997]) - v[2773] * v[2999] + 2e0*v[2750] * v[3000] + v[3495] + v[3501] + v[3507] + v[3513]
//			+ 2e0*v[3416] * v[481] - v[3425] * v[484] - v[3412] * v[487];
//		v[3517] = v[2728] * v[3057] + v[198] * v[3351];
//		v[3518] = v[2728] * v[3058] + v[199] * v[3351];
//		v[3519] = v[2728] * v[3059] + v[200] * v[3351];
//		v[3520] = QBi[0][0] * v[3517] + QBi[0][1] * v[3518] + QBi[0][2] * v[3519];
//		v[3521] = QBi[2][0] * v[3517] + QBi[2][1] * v[3518] + QBi[2][2] * v[3519];
//		v[3522] = QBi[1][0] * v[3517] + QBi[1][1] * v[3518] + QBi[1][2] * v[3519];
//		v[3523] = v[2777] * v[2997] + 2e0*v[2753] * v[2999] - v[2773] * v[3000] - v[3484] - v[3490] - v[3514] - v[3520]
//			- v[3425] * v[481] + 2e0*v[3421] * v[484] + v[3428] * v[487];
//		v[3524] = -(v[2816] * v[2990]) + v[2812] * v[2992] - v[2817] * v[2993] - 2e0*v[2766] * v[2998] + v[204] * v[3476]
//			- v[3495] * v[474] + v[3484] * v[475] - v[3491] * v[476];
//		v[3525] = v[2731] * v[3057] + v[198] * v[3330];
//		v[3526] = v[2731] * v[3058] + v[199] * v[3330];
//		v[3527] = v[2731] * v[3059] + v[200] * v[3330];
//		v[3528] = QBi[0][0] * v[3525] + QBi[0][1] * v[3526] + QBi[0][2] * v[3527];
//		v[3529] = QBi[1][0] * v[3525] + QBi[1][1] * v[3526] + QBi[1][2] * v[3527];
//		v[3530] = QBi[2][0] * v[3525] + QBi[2][1] * v[3526] + QBi[2][2] * v[3527];
//		v[3531] = v[2732] * v[3057] + v[198] * v[3329];
//		v[3532] = v[2732] * v[3058] + v[199] * v[3329];
//		v[3533] = v[2732] * v[3059] + v[200] * v[3329];
//		v[3534] = QBi[1][0] * v[3531] + QBi[1][1] * v[3532] + QBi[1][2] * v[3533];
//		v[3535] = QBi[0][0] * v[3531] + QBi[0][1] * v[3532] + QBi[0][2] * v[3533];
//		v[3536] = QBi[2][0] * v[3531] + QBi[2][1] * v[3532] + QBi[2][2] * v[3533];
//		v[3537] = 2e0*v[2747] * v[2997] + v[2777] * v[2999] - v[2751] * v[3000] - v[3496] - v[3521] - v[3528] - v[3534]
//			- v[3412] * v[481] + v[3428] * v[484] + 2e0*v[3406] * v[487];
//		v[3538] = -(v[2815] * v[2990]) + v[2813] * v[2992] - v[2818] * v[2993] - 2e0*v[2764] * v[2998] + v[204] * v[3474]
//			- v[3496] * v[474] + v[3485] * v[475] - v[3490] * v[476];
//		v[3539] = -(v[2824] * v[2990]) + v[2830] * v[2992] - v[2820] * v[2993] - 2e0*v[2763] * v[2998] + v[204] * v[3473]
//			- v[3501] * v[474] + v[3520] * v[475] - v[3508] * v[476];
//		v[3540] = v[3489] + v[3509];
//		v[3541] = -(v[2823] * v[2990]) + v[2831] * v[2992] - v[2822] * v[2993] - 2e0*v[2772] * v[2998] + v[204] * v[3466]
//			- v[3503] * v[474] + v[3521] * v[475] - v[3507] * v[476];
//		v[3542] = -(v[2827] * v[2990]) + v[2839] * v[2992] - v[2835] * v[2993] - 2e0*v[2759] * v[2998] + v[204] * v[3464]
//			- v[3528] * v[474] + v[3535] * v[475] - v[3514] * v[476];
//		v[3543] = -(v[2826] * v[2990]) + v[2838] * v[2992] - v[2836] * v[2993] - 2e0*v[2776] * v[2998] + v[204] * v[3462]
//			- v[3529] * v[474] + v[3534] * v[475] - v[3513] * v[476];
//		v[3544] = v[3502] + v[3530];
//		v[3545] = v[3483] + v[3536];
//		v[3546] = -(QBi[1][2] * v[3170]) - QBi[0][2] * v[3171] - QBi[2][2] * v[3179] - v[231] * v[3376] - v[228] * v[3378]
//			- v[225] * v[3380] + v[3365] * v[525] + v[3367] * v[526] + v[3372] * v[527] + v[3348] * v[540] + v[3351] * v[541]
//			+ v[3357] * v[542] + v[3330] * v[555] + v[3329] * v[556] + v[3328] * v[557];
//		v[3547] = -(QBi[1][1] * v[3170]) - QBi[0][1] * v[3171] - QBi[2][1] * v[3179] - v[230] * v[3376] - v[227] * v[3378]
//			- v[224] * v[3380] + v[3365] * v[522] + v[3367] * v[523] + v[3372] * v[524] + v[3348] * v[537] + v[3351] * v[538]
//			+ v[3357] * v[539] + v[3330] * v[552] + v[3329] * v[553] + v[3328] * v[554];
//		v[3548] = -(QBi[1][0] * v[3170]) - QBi[0][0] * v[3171] - QBi[2][0] * v[3179] - v[229] * v[3376] - v[226] * v[3378]
//			- v[223] * v[3380] + v[3365] * v[519] + v[3367] * v[520] + v[3372] * v[521] + v[3348] * v[534] + v[3351] * v[535]
//			+ v[3357] * v[536] + v[3330] * v[549] + v[3329] * v[550] + v[3328] * v[551];
//		v[3549] = v[3548] * x1B[0] + v[3547] * x1B[1] + v[3546] * x1B[2];
//		v[3550] = QAi[2][2] * v[3384] + QAi[2][1] * v[3385] + QAi[2][0] * v[3386] + v[1571] * v[3437] + v[266] * v[3556]
//			+ v[3448] * v[3727] + v[264] * (v[3136] / v[253] + v[3041] * v[4468]);
//		v[3553] = (v[2788] * v[3107]) / v[253] + QAi[1][2] * v[3384] + QAi[1][1] * v[3385] + QAi[1][0] * v[3386]
//			+ v[2684] * v[3429] + v[263] * v[3560] + v[3139] * v[3716] + v[3449] * v[3727] + v[3041] * v[4469];
//		v[3555] = QAi[0][2] * v[3384] + QAi[0][1] * v[3385] + QAi[0][0] * v[3386] + v[2684] * v[3430] + v[251] * v[3563]
//			+ v[256] * (v[4149] + v[3041] * v[4470]);
//		v[3557] = (v[2787] * v[3110]) / v[253] + QAi[2][2] * v[3390] + QAi[2][1] * v[3391] + QAi[2][0] * v[3392]
//			+ v[2683] * v[3437] + v[3147] * v[3717] + v[3452] * v[3724] + v[3556] * v[3728] + v[3041] * v[4471];
//		v[3561] = QAi[1][2] * v[3390] + QAi[1][1] * v[3391] + QAi[1][0] * v[3392] + v[1573] * v[3429] + v[257] * v[3560]
//			+ v[3447] * v[3724] + v[262] * (v[3149] / v[253] + v[3041] * v[4472]);
//		v[3564] = QAi[0][2] * v[3390] + QAi[0][1] * v[3391] + QAi[0][0] * v[3392] + v[2683] * v[3438] + v[252] * v[3563]
//			+ v[256] * (v[4148] + v[3041] * v[4473]);
//		v[3565] = (v[2782] * v[3111]) / v[253] + QAi[2][2] * v[3396] + QAi[2][1] * v[3397] + QAi[2][0] * v[3398]
//			+ v[1572] * v[3441] + v[2682] * v[3448] + v[3143] * v[3717] + v[3556] * v[3726] + v[3041] * v[4474];
//		v[3567] = QAi[1][2] * v[3396] + QAi[1][1] * v[3397] + QAi[1][0] * v[3398] + v[2682] * v[3447] + v[3560] * v[3723]
//			+ v[262] * (v[3041] * v[4130] + v[4152]);
//		v[3570] = QAi[0][2] * v[3396] + QAi[0][1] * v[3397] + QAi[0][0] * v[3398] + v[1574] * v[3430] + v[1572] * v[3438]
//			+ v[254] * v[3563] + v[256] * (v[3145] / v[253] + v[3041] * v[4475]);
//		v[3571] = v[2733] * v[3007] + v[142] * v[3361];
//		v[3572] = v[2733] * v[3008] + v[143] * v[3361];
//		v[3573] = v[2733] * v[3009] + v[144] * v[3361];
//		v[3574] = QAi[0][0] * v[3571] + QAi[0][1] * v[3572] + QAi[0][2] * v[3573];
//		v[3575] = QAi[1][0] * v[3571] + QAi[1][1] * v[3572] + QAi[1][2] * v[3573];
//		v[3576] = QAi[2][0] * v[3571] + QAi[2][1] * v[3572] + QAi[2][2] * v[3573];
//		v[3577] = v[2734] * v[3007] + v[142] * v[3363];
//		v[3578] = v[2734] * v[3008] + v[143] * v[3363];
//		v[3579] = v[2734] * v[3009] + v[144] * v[3363];
//		v[3580] = QAi[0][0] * v[3577] + QAi[0][1] * v[3578] + QAi[0][2] * v[3579];
//		v[3581] = QAi[2][0] * v[3577] + QAi[2][1] * v[3578] + QAi[2][2] * v[3579];
//		v[3582] = QAi[1][0] * v[3577] + QAi[1][1] * v[3578] + QAi[1][2] * v[3579];
//		v[3583] = v[2735] * v[3007] + v[142] * v[3359];
//		v[3584] = v[2735] * v[3008] + v[143] * v[3359];
//		v[3585] = v[2735] * v[3009] + v[144] * v[3359];
//		v[3586] = QAi[1][0] * v[3583] + QAi[1][1] * v[3584] + QAi[1][2] * v[3585];
//		v[3587] = QAi[2][0] * v[3583] + QAi[2][1] * v[3584] + QAi[2][2] * v[3585];
//		v[3588] = QAi[0][0] * v[3583] + QAi[0][1] * v[3584] + QAi[0][2] * v[3585];
//		v[3589] = v[2736] * v[3007] + v[142] * v[3339];
//		v[3590] = v[2736] * v[3008] + v[143] * v[3339];
//		v[3591] = v[2736] * v[3009] + v[144] * v[3339];
//		v[3592] = QAi[0][0] * v[3589] + QAi[0][1] * v[3590] + QAi[0][2] * v[3591];
//		v[3593] = QAi[1][0] * v[3589] + QAi[1][1] * v[3590] + QAi[1][2] * v[3591];
//		v[3594] = QAi[2][0] * v[3589] + QAi[2][1] * v[3590] + QAi[2][2] * v[3591];
//		v[3595] = v[2737] * v[3007] + v[142] * v[3345];
//		v[3596] = v[2737] * v[3008] + v[143] * v[3345];
//		v[3597] = v[2737] * v[3009] + v[144] * v[3345];
//		v[3598] = QAi[2][0] * v[3595] + QAi[2][1] * v[3596] + QAi[2][2] * v[3597];
//		v[3599] = QAi[0][0] * v[3595] + QAi[0][1] * v[3596] + QAi[0][2] * v[3597];
//		v[3600] = QAi[1][0] * v[3595] + QAi[1][1] * v[3596] + QAi[1][2] * v[3597];
//		v[3601] = v[2738] * v[3007] + v[142] * v[3331];
//		v[3602] = v[2738] * v[3008] + v[143] * v[3331];
//		v[3603] = v[2738] * v[3009] + v[144] * v[3331];
//		v[3604] = QAi[1][0] * v[3601] + QAi[1][1] * v[3602] + QAi[1][2] * v[3603];
//		v[3605] = QAi[0][0] * v[3601] + QAi[0][1] * v[3602] + QAi[0][2] * v[3603];
//		v[3606] = QAi[2][0] * v[3601] + QAi[2][1] * v[3602] + QAi[2][2] * v[3603];
//		v[3607] = -(v[2784] * v[2971]) - v[2806] * v[2973] + 2e0*v[2783] * v[2974] + v[3586] + v[3592] + v[3598] + v[3604]
//			+ 2e0*v[3446] * v[382] - v[3455] * v[385] - v[3442] * v[388];
//		v[3608] = v[2739] * v[3007] + v[142] * v[3342];
//		v[3609] = v[2739] * v[3008] + v[143] * v[3342];
//		v[3610] = v[2739] * v[3009] + v[144] * v[3342];
//		v[3611] = QAi[0][0] * v[3608] + QAi[0][1] * v[3609] + QAi[0][2] * v[3610];
//		v[3612] = QAi[2][0] * v[3608] + QAi[2][1] * v[3609] + QAi[2][2] * v[3610];
//		v[3613] = QAi[1][0] * v[3608] + QAi[1][1] * v[3609] + QAi[1][2] * v[3610];
//		v[3614] = v[2810] * v[2971] + 2e0*v[2786] * v[2973] - v[2806] * v[2974] - v[3575] - v[3581] - v[3605] - v[3611]
//			- v[3455] * v[382] + 2e0*v[3451] * v[385] + v[3458] * v[388];
//		v[3615] = -(v[2860] * v[2964]) + v[2856] * v[2966] - v[2861] * v[2967] - 2e0*v[2799] * v[2972] + v[148] * v[3567]
//			- v[3586] * v[375] + v[3575] * v[376] - v[3582] * v[377];
//		v[3616] = v[2740] * v[3007] + v[142] * v[3336];
//		v[3617] = v[2740] * v[3008] + v[143] * v[3336];
//		v[3618] = v[2740] * v[3009] + v[144] * v[3336];
//		v[3619] = QAi[0][0] * v[3616] + QAi[0][1] * v[3617] + QAi[0][2] * v[3618];
//		v[3620] = QAi[1][0] * v[3616] + QAi[1][1] * v[3617] + QAi[1][2] * v[3618];
//		v[3621] = QAi[2][0] * v[3616] + QAi[2][1] * v[3617] + QAi[2][2] * v[3618];
//		v[3622] = v[2741] * v[3007] + v[142] * v[3332];
//		v[3623] = v[2741] * v[3008] + v[143] * v[3332];
//		v[3624] = v[2741] * v[3009] + v[144] * v[3332];
//		v[3625] = QAi[1][0] * v[3622] + QAi[1][1] * v[3623] + QAi[1][2] * v[3624];
//		v[3626] = QAi[0][0] * v[3622] + QAi[0][1] * v[3623] + QAi[0][2] * v[3624];
//		v[3627] = QAi[2][0] * v[3622] + QAi[2][1] * v[3623] + QAi[2][2] * v[3624];
//		v[3628] = 2e0*v[2780] * v[2971] + v[2810] * v[2973] - v[2784] * v[2974] - v[3587] - v[3612] - v[3619] - v[3625]
//			- v[3442] * v[382] + v[3458] * v[385] + 2e0*v[3436] * v[388];
//		v[3629] = -(v[2859] * v[2964]) + v[2857] * v[2966] - v[2862] * v[2967] - 2e0*v[2797] * v[2972] + v[148] * v[3565]
//			- v[3587] * v[375] + v[3576] * v[376] - v[3581] * v[377];
//		v[3630] = -(v[2868] * v[2964]) + v[2874] * v[2966] - v[2864] * v[2967] - 2e0*v[2796] * v[2972] + v[148] * v[3564]
//			- v[3592] * v[375] + v[3611] * v[376] - v[3599] * v[377];
//		v[3631] = v[3580] + v[3600];
//		v[3632] = -(v[2867] * v[2964]) + v[2875] * v[2966] - v[2866] * v[2967] - 2e0*v[2805] * v[2972] + v[148] * v[3557]
//			- v[3594] * v[375] + v[3612] * v[376] - v[3598] * v[377];
//		v[3633] = -(v[2871] * v[2964]) + v[2883] * v[2966] - v[2879] * v[2967] - 2e0*v[2792] * v[2972] + v[148] * v[3555]
//			- v[3619] * v[375] + v[3626] * v[376] - v[3605] * v[377];
//		v[3634] = -(v[2870] * v[2964]) + v[2882] * v[2966] - v[2880] * v[2967] - 2e0*v[2809] * v[2972] + v[148] * v[3553]
//			- v[3620] * v[375] + v[3625] * v[376] - v[3604] * v[377];
//		v[3635] = v[3593] + v[3621];
//		v[3636] = v[3574] + v[3627];
//		v[3637] = QAi[1][2] * v[3123] + QAi[0][2] * v[3124] + QAi[2][2] * v[3130] + v[175] * v[3376] + v[172] * v[3378]
//			+ v[169] * v[3380] + v[3359] * v[426] + v[3361] * v[427] + v[3363] * v[428] + v[3339] * v[441] + v[3342] * v[442]
//			+ v[3345] * v[443] + v[3336] * v[456] + v[3332] * v[457] + v[3331] * v[458];
//		v[3638] = QAi[1][1] * v[3123] + QAi[0][1] * v[3124] + QAi[2][1] * v[3130] + v[174] * v[3376] + v[171] * v[3378]
//			+ v[168] * v[3380] + v[3359] * v[423] + v[3361] * v[424] + v[3363] * v[425] + v[3339] * v[438] + v[3342] * v[439]
//			+ v[3345] * v[440] + v[3336] * v[453] + v[3332] * v[454] + v[3331] * v[455];
//		v[3639] = QAi[1][0] * v[3123] + QAi[0][0] * v[3124] + QAi[2][0] * v[3130] + v[173] * v[3376] + v[170] * v[3378]
//			+ v[167] * v[3380] + v[3359] * v[420] + v[3361] * v[421] + v[3363] * v[422] + v[3339] * v[435] + v[3342] * v[436]
//			+ v[3345] * v[437] + v[3336] * v[450] + v[3332] * v[451] + v[3331] * v[452];
//		v[3640] = v[3639] * x1A[0] + v[3638] * x1A[1] + v[3637] * x1A[2];
//		v[3641] = (-v[3549] + v[3548] * x3B[0] + v[3547] * x3B[1] + v[3546] * x3B[2]) / 2e0;
//		v[3642] = (-v[3549] + v[3548] * x2B[0] + v[3547] * x2B[1] + v[3546] * x2B[2]) / 2e0;
//		v[3643] = (-(v[2837] * v[2942]) - v[2819] * v[2945] - 2e0*v[2820] * v[2981] - 2e0*v[2817] * v[2983]
//			- 2e0*v[2836] * v[2984] - 2e0*v[2822] * v[2986] - 2e0*v[2835] * v[2987] - 2e0*v[2818] * v[2989] - v[2821] * v[2991]
//			- 2e0*v[3402] + 2e0*v[3410] - v[3497] * v[477] - 2e0*v[3495] * v[482] - 2e0*v[3496] * v[488] - 2e0*v[3501] * v[492]
//			- v[3502] * v[496] - 2e0*v[3503] * v[501] - 2e0*v[3528] * v[505] - 2e0*v[3529] * v[509] - v[3530] * v[514]) / 2e0;
//		v[3644] = (-(v[2814] * v[2990]) + v[2811] * v[2992] - v[2819] * v[2993] - 2e0*v[2769] * v[2998] + v[204] * v[3479]
//			- v[3497] * v[474] + v[3483] * v[475] - v[3489] * v[476]) / 2e0;
//		v[3646] = (v[2840] * v[2942] + v[2811] * v[2945] + 2e0*v[2830] * v[2981] + 2e0*v[2812] * v[2983]
//			+ 2e0*v[2838] * v[2984] + 2e0*v[2831] * v[2986] + 2e0*v[2839] * v[2987] + 2e0*v[2813] * v[2989] + v[2832] * v[2991]
//			- 2e0*v[3404] + 2e0*v[3427] + v[3483] * v[477] + 2e0*v[3484] * v[482] + 2e0*v[3485] * v[488] + 2e0*v[3520] * v[492]
//			+ v[3522] * v[496] + 2e0*v[3521] * v[501] + 2e0*v[3535] * v[505] + 2e0*v[3534] * v[509] + v[3536] * v[514]) / 2e0;
//		v[3647] = (-(v[2825] * v[2990]) + v[2832] * v[2992] - v[2821] * v[2993] - 2e0*v[2760] * v[2998] + v[204] * v[3470]
//			- v[3502] * v[474] + v[3522] * v[475] - v[3509] * v[476]) / 2e0;
//		v[3648] = (-(v[2828] * v[2942]) - v[2814] * v[2945] - 2e0*v[2824] * v[2981] - 2e0*v[2816] * v[2983]
//			- 2e0*v[2826] * v[2984] - 2e0*v[2823] * v[2986] - 2e0*v[2827] * v[2987] - 2e0*v[2815] * v[2989] - v[2825] * v[2991]
//			- 2e0*v[3414] + 2e0*v[3424] - v[3489] * v[477] - 2e0*v[3491] * v[482] - 2e0*v[3490] * v[488] - 2e0*v[3508] * v[492]
//			- v[3509] * v[496] - 2e0*v[3507] * v[501] - 2e0*v[3514] * v[505] - 2e0*v[3513] * v[509] - v[3515] * v[514]) / 2e0;
//		v[3701] = v[1658] * v[3646] - v[3647] + v[1655] * v[3648] + 24e0*v[2943] * v[3650] * v[4477] + v[4141] * (-
//			(v[2756] * v[2942]) - v[2769] * v[2945] - 2e0*v[2763] * v[2981] - 2e0*v[2766] * v[2983] - 2e0*v[2776] * v[2984]
//			- 2e0*v[2772] * v[2986] - 2e0*v[2759] * v[2987] - 2e0*v[2764] * v[2989] - v[2760] * v[2991] - 2e0*v[3485]
//			+ 2e0*v[3491] + 2e0*v[3503] - 2e0*v[3508] - alphaB[1] * v[3516] + alphaB[0] * v[3523] - 2e0*v[3529] + 2e0*v[3535]
//			+ alphaB[2] * v[3537] + 8e0*v[2998] * v[3651] + v[3544] * v[369] + v[3545] * v[372] + v[3540] * v[374] - 4e0*v[204] *
//			(v[2421] * v[3181] + v[2425] * v[3184] + v[2409] * v[3185] + v[2408] * v[3188] + v[2428] * v[3190] + v[2422] * v[3192]
//				+ v[2413] * v[3194] + v[2416] * v[3196] + v[2415] * v[3198] + v[2608] * v[3399] + v[2609] * v[3400] + v[3406]
//				+ v[2610] * v[3407] + v[2611] * v[3408] + v[1859] * v[3411] + v[3416] + v[2612] * v[3417] + v[2613] * v[3418]
//				+ v[1853] * v[3419] + v[3421] + v[1857] * v[3422] + v[3178] * v[3460] + v[3176] * v[3461] + v[3174] * v[3463]
//				+ v[3169] * v[3467] + v[3166] * v[3468] + v[3165] * v[3471] + v[3159] * v[3475] + v[3162] * v[3477] + v[3156] * v[3478]
//				+ v[3065] * v[3656] + v[3070] * v[3657] + v[3076] * v[3658] + v[3150] * v[3659] + v[3153] * v[3660] + v[3154] * v[3661]
//				+ v[3321] * v[3977] + v[3322] * v[3978] + v[3323] * v[3979] + 2e0*v[3091] * v[3652] * v[4479]) - v[3479] * v[477]
//			- 2e0*v[3476] * v[482] - 2e0*v[3474] * v[488] - 2e0*v[3473] * v[492] - v[3470] * v[496] - 2e0*v[3466] * v[501]
//			- 2e0*v[3464] * v[505] - 2e0*v[3462] * v[509] - v[3459] * v[514] - 2e0*v[7033 + i2691]) - 2e0*v[1564] * (
//				-4e0*v[2853] * v[2943] + 4e0*v[3643] * v[369] + v[7057 + i2691]);
//		v[4146] = v[3701] + (v[2828] * v[2990] - v[2840] * v[2992] + v[2837] * v[2993] + 2e0*v[2756] * v[2998]
//			- v[204] * v[3459] + v[3530] * v[474] - v[3536] * v[475] + v[3515] * v[476]) / 2e0;
//		v[4145] = (v[3538] + v[3542]) / 2e0;
//		v[3664] = v[3541] + v[3543];
//		v[3665] = v[3524] + v[3539];
//		v[3666] = (-v[3640] + v[3639] * x2A[0] + v[3638] * x2A[1] + v[3637] * x2A[2]) / 2e0;
//		v[3667] = (-v[3640] + v[3639] * x3A[0] + v[3638] * x3A[1] + v[3637] * x3A[2]) / 2e0;
//		v[3668] = (-(v[2881] * v[2934]) - v[2863] * v[2937] - 2e0*v[2864] * v[2955] - 2e0*v[2861] * v[2957]
//			- 2e0*v[2880] * v[2958] - 2e0*v[2866] * v[2960] - 2e0*v[2879] * v[2961] - 2e0*v[2862] * v[2963] - v[2865] * v[2965]
//			- 2e0*v[3432] + 2e0*v[3440] - v[3588] * v[378] - 2e0*v[3586] * v[383] - 2e0*v[3587] * v[389] - 2e0*v[3592] * v[393]
//			- v[3593] * v[397] - 2e0*v[3594] * v[402] - 2e0*v[3619] * v[406] - 2e0*v[3620] * v[410] - v[3621] * v[415]) / 2e0;
//		v[3669] = (-(v[2858] * v[2964]) + v[2855] * v[2966] - v[2863] * v[2967] - 2e0*v[2802] * v[2972] + v[148] * v[3570]
//			- v[3588] * v[375] + v[3574] * v[376] - v[3580] * v[377]) / 2e0;
//		v[3671] = (v[2884] * v[2934] + v[2855] * v[2937] + 2e0*v[2874] * v[2955] + 2e0*v[2856] * v[2957]
//			+ 2e0*v[2882] * v[2958] + 2e0*v[2875] * v[2960] + 2e0*v[2883] * v[2961] + 2e0*v[2857] * v[2963] + v[2876] * v[2965]
//			- 2e0*v[3434] + 2e0*v[3457] + v[3574] * v[378] + 2e0*v[3575] * v[383] + 2e0*v[3576] * v[389] + 2e0*v[3611] * v[393]
//			+ v[3613] * v[397] + 2e0*v[3612] * v[402] + 2e0*v[3626] * v[406] + 2e0*v[3625] * v[410] + v[3627] * v[415]) / 2e0;
//		v[3672] = (-(v[2869] * v[2964]) + v[2876] * v[2966] - v[2865] * v[2967] - 2e0*v[2793] * v[2972] + v[148] * v[3561]
//			- v[3593] * v[375] + v[3613] * v[376] - v[3600] * v[377]) / 2e0;
//		v[3673] = (-(v[2872] * v[2934]) - v[2858] * v[2937] - 2e0*v[2868] * v[2955] - 2e0*v[2860] * v[2957]
//			- 2e0*v[2870] * v[2958] - 2e0*v[2867] * v[2960] - 2e0*v[2871] * v[2961] - 2e0*v[2859] * v[2963] - v[2869] * v[2965]
//			- 2e0*v[3444] + 2e0*v[3454] - v[3580] * v[378] - 2e0*v[3582] * v[383] - 2e0*v[3581] * v[389] - 2e0*v[3599] * v[393]
//			- v[3600] * v[397] - 2e0*v[3598] * v[402] - 2e0*v[3605] * v[406] - 2e0*v[3604] * v[410] - v[3606] * v[415]) / 2e0;
//		v[3694] = v[1632] * v[3671] - v[3672] + v[1629] * v[3673] + 24e0*v[2935] * v[3675] * v[4481] + v[4135] * (-
//			(v[2789] * v[2934]) - v[2802] * v[2937] - 2e0*v[2796] * v[2955] - 2e0*v[2799] * v[2957] - 2e0*v[2809] * v[2958]
//			- 2e0*v[2805] * v[2960] - 2e0*v[2792] * v[2961] - 2e0*v[2797] * v[2963] - v[2793] * v[2965] - 2e0*v[3576]
//			+ 2e0*v[3582] + 2e0*v[3594] - 2e0*v[3599] - alphaA[1] * v[3607] + alphaA[0] * v[3614] - 2e0*v[3620] + 2e0*v[3626]
//			+ alphaA[2] * v[3628] + v[363] * v[3635] + v[3636] * v[366] + 8e0*v[2972] * v[3676] + v[3631] * v[368] - v[3570] * v[378]
//			- 2e0*v[3567] * v[383] - 2e0*v[3565] * v[389] - 2e0*v[3564] * v[393] - v[3561] * v[397] - 2e0*v[3557] * v[402]
//			- 2e0*v[3555] * v[406] - 2e0*v[3553] * v[410] - v[3550] * v[415] - 4e0*v[148] * (v[2451] * v[3132] + v[2455] * v[3135]
//				+ v[2439] * v[3136] + v[2438] * v[3139] + v[2458] * v[3141] + v[2452] * v[3143] + v[2443] * v[3145] + v[2446] * v[3147]
//				+ v[2445] * v[3149] + v[2637] * v[3429] + v[2638] * v[3430] + v[3436] + v[2639] * v[3437] + v[2640] * v[3438]
//				+ v[1831] * v[3441] + v[3446] + v[2641] * v[3447] + v[2642] * v[3448] + v[1825] * v[3449] + v[3451] + v[1829] * v[3452]
//				+ v[3129] * v[3551] + v[3127] * v[3552] + v[3125] * v[3554] + v[3122] * v[3558] + v[3119] * v[3559] + v[3118] * v[3562]
//				+ v[3114] * v[3566] + v[3117] * v[3568] + v[3113] * v[3569] + v[3015] * v[3681] + v[3020] * v[3682] + v[3026] * v[3683]
//				+ v[3107] * v[3684] + v[3110] * v[3685] + v[3111] * v[3686] + v[3324] * v[3980] + v[3325] * v[3981] + v[3326] * v[3982]
//				+ 2e0*v[3041] * v[3677] * v[4483]) - 2e0*v[7045 + i2691]) - 2e0*v[1562] * (-4e0*v[2897] * v[2935]
//					+ 4e0*v[363] * v[3668] + v[7105 + i2691]);
//		v[4140] = v[3694] + (v[2872] * v[2964] - v[2884] * v[2966] + v[2881] * v[2967] + 2e0*v[2789] * v[2972]
//			- v[148] * v[3550] + v[3621] * v[375] - v[3627] * v[376] + v[3606] * v[377]) / 2e0;
//		v[4139] = (v[3629] + v[3633]) / 2e0;
//		v[3689] = v[3632] + v[3634];
//		v[3690] = v[3615] + v[3630];
//		v[7142] = 0e0;
//		v[7143] = 0e0;
//		v[7144] = 0e0;
//		v[7145] = 2e0*v[3691];
//		v[7146] = v[3692];
//		v[7147] = v[3693];
//		v[7148] = 0e0;
//		v[7149] = 0e0;
//		v[7150] = 0e0;
//		v[7151] = 0e0;
//		v[7152] = 0e0;
//		v[7153] = 0e0;
//		v[7130] = 0e0;
//		v[7131] = 0e0;
//		v[7132] = 0e0;
//		v[7133] = v[3692];
//		v[7134] = 2e0*v[3695];
//		v[7135] = v[3696];
//		v[7136] = 0e0;
//		v[7137] = 0e0;
//		v[7138] = 0e0;
//		v[7139] = 0e0;
//		v[7140] = 0e0;
//		v[7141] = 0e0;
//		v[7118] = 0e0;
//		v[7119] = 0e0;
//		v[7120] = 0e0;
//		v[7121] = v[3693];
//		v[7122] = v[3696];
//		v[7123] = 2e0*v[3697];
//		v[7124] = 0e0;
//		v[7125] = 0e0;
//		v[7126] = 0e0;
//		v[7127] = 0e0;
//		v[7128] = 0e0;
//		v[7129] = 0e0;
//		v[7070] = 0e0;
//		v[7071] = 0e0;
//		v[7072] = 0e0;
//		v[7073] = 0e0;
//		v[7074] = 0e0;
//		v[7075] = 0e0;
//		v[7076] = 0e0;
//		v[7077] = 0e0;
//		v[7078] = 0e0;
//		v[7079] = 2e0*v[3698];
//		v[7080] = v[3699];
//		v[7081] = v[3700];
//		v[7082] = 0e0;
//		v[7083] = 0e0;
//		v[7084] = 0e0;
//		v[7085] = 0e0;
//		v[7086] = 0e0;
//		v[7087] = 0e0;
//		v[7088] = 0e0;
//		v[7089] = 0e0;
//		v[7090] = 0e0;
//		v[7091] = v[3699];
//		v[7092] = 2e0*v[3702];
//		v[7093] = v[3703];
//		v[7094] = 0e0;
//		v[7095] = 0e0;
//		v[7096] = 0e0;
//		v[7097] = 0e0;
//		v[7098] = 0e0;
//		v[7099] = 0e0;
//		v[7100] = 0e0;
//		v[7101] = 0e0;
//		v[7102] = 0e0;
//		v[7103] = v[3700];
//		v[7104] = v[3703];
//		v[7105] = 2e0*v[3704];
//		v[7154] = v[3380];
//		v[7155] = v[3378];
//		v[7156] = v[3376];
//		v[7157] = v[2877] * v[2972] - v[3632] + v[3634] + v[364] * v[3690] + v[3614] * v[3983] + 2e0*(v[3635] * v[3983]
//			+ v[3668] * v[3985] + v[2935] * (v[2895] * v[3992] - v[1562] * v[4138])) + alphaA[2] * v[4139]
//			+ 2e0*alphaA[0] * v[4140] + v[7141 + i2691];
//		v[7158] = -(v[2873] * v[2972]) + v[3629] - v[3633] + (alphaA[2] * v[3689]) / 2e0 + (alphaA[0] * v[3690]) / 2e0
//			- v[3607] * v[3983] + 2e0*alphaA[1] * (-v[3669] + v[3672] + v[4140]) + v[7129 + i2691] + 2e0*(v[3636] * v[3983]
//				+ v[2935] * (v[2896] * v[3992] + v[1562] * v[4137]) - 4e0*v[3671] * v[962]);
//		v[7159] = v[2888] * v[2972] - v[3615] + v[3630] + v[364] * v[3689] + 2e0*alphaA[2] * (-v[3669] + v[3694])
//			+ v[3628] * v[3983] + 2e0*(v[3631] * v[3983] + v[3673] * v[3985] + v[2935] * (v[2891] * v[4135] - v[1562] * v[4136]))
//			+ alphaA[0] * v[4139] + v[7117 + i2691];
//		v[7160] = -v[3380];
//		v[7161] = -v[3378];
//		v[7162] = -v[3376];
//		v[7163] = v[2833] * v[2998] - v[3541] + v[3543] + v[3665] * v[370] + v[3523] * v[3986] + 2e0*(v[3544] * v[3986]
//			+ v[3643] * v[3988] + v[2943] * (v[2851] * v[3998] - v[1564] * v[4144])) + alphaB[2] * v[4145]
//			+ 2e0*alphaB[0] * v[4146] + v[7069 + i2691];
//		v[7164] = -(v[2829] * v[2998]) + v[3538] - v[3542] + (alphaB[2] * v[3664]) / 2e0 + (alphaB[0] * v[3665]) / 2e0
//			- v[3516] * v[3986] + 2e0*alphaB[1] * (-v[3644] + v[3647] + v[4146]) + v[7081 + i2691] + 2e0*(v[3545] * v[3986]
//				+ v[2943] * (v[2852] * v[3998] + v[1564] * v[4143]) - 4e0*v[3646] * v[968]);
//		v[7165] = v[2844] * v[2998] - v[3524] + v[3539] + v[3664] * v[370] + 2e0*alphaB[2] * (-v[3644] + v[3701])
//			+ v[3537] * v[3986] + 2e0*(v[3540] * v[3986] + v[3648] * v[3988] + v[2943] * (v[2847] * v[4141] - v[1564] * v[4142]))
//			+ alphaB[0] * v[4145] + v[7093 + i2691];
//		Rc[i2691 - 1] += v[2910] * v[2946] + v[2911] * v[2947] + v[2899] * v[2949] + v[2900] * v[2950] + v[6861 + i2691] + (*a4
//			)*v[6873 + i2691];
//		for (i2928 = 1; i2928 <= 12; i2928++) {
//			Kc[i2691 - 1][i2928 - 1] += v[3667] * v[5357 + i2928] + v[3666] * v[5369 + i2928] + v[3642] * v[5381 + i2928]
//				+ v[3641] * v[5393 + i2928] + v[7153 + i2928] + (*a4)*v[7165 + i2928];
//		};/* end for */
//	};/* end for */
//#pragma endregion


#pragma region AceGen
double v01; double v010; double v011; double v012; double v013; double v014;
double v015; double v016; double v017; double v018; double v019; double v02;
double v020; double v021; double v022; double v023; double v024; double v025;
double v026; double v027; double v03; double v04; double v05; double v06; double v07;
double v08; double v09;
int i840, i933, i1194, i1799, i3212, i3510, i4599, i4600, i4601, i4602, i4603, i4604, i4618
, i4619, i4620, b29, b30, b329, b678, b733, b735
, b752, b771, b805, b845, b846, b850, b1023, b1024
, b1035, b1199, b1200, b1201, b1202, b1269, b1273, b1274, b1291, b1292, b1446, b1480, b1533
, b2109, b2151, b2197, b2519, b2520, b2529, b2530, b2554, b2555, b2556, b2572, b2667, b2671
, b2675, b2686, b2687, b2788, b2813, b2869, b3221, b3222, b3226, b3231, b3232, b3317, b3789
, b3971, b3972, b3980, b3984, b3988, b3998, b3999, b4085;
v[7654] = 0e0;
v[7655] = 0e0;
v[7656] = -1e0;
v[7657] = 0e0;
v[7658] = 0e0;
v[7659] = 0e0;
v[7660] = 0e0;
v[7661] = 0e0;
v[7662] = 1e0;
v[7663] = 0e0;
v[7664] = 0e0;
v[7665] = 0e0;
v[7642] = 0e0;
v[7643] = -1e0;
v[7644] = 0e0;
v[7645] = 0e0;
v[7646] = 0e0;
v[7647] = 0e0;
v[7648] = 0e0;
v[7649] = 1e0;
v[7650] = 0e0;
v[7651] = 0e0;
v[7652] = 0e0;
v[7653] = 0e0;
v[7630] = -1e0;
v[7631] = 0e0;
v[7632] = 0e0;
v[7633] = 0e0;
v[7634] = 0e0;
v[7635] = 0e0;
v[7636] = 1e0;
v[7637] = 0e0;
v[7638] = 0e0;
v[7639] = 0e0;
v[7640] = 0e0;
v[7641] = 0e0;
v[6558] = 0e0;
v[6559] = 0e0;
v[6560] = 1e0;
v[6561] = 0e0;
v[6562] = 0e0;
v[6563] = 0e0;
v[6564] = 0e0;
v[6565] = 0e0;
v[6566] = -1e0;
v[6567] = 0e0;
v[6568] = 0e0;
v[6569] = 0e0;
v[6546] = 0e0;
v[6547] = 1e0;
v[6548] = 0e0;
v[6549] = 0e0;
v[6550] = 0e0;
v[6551] = 0e0;
v[6552] = 0e0;
v[6553] = -1e0;
v[6554] = 0e0;
v[6555] = 0e0;
v[6556] = 0e0;
v[6557] = 0e0;
v[6534] = 1e0;
v[6535] = 0e0;
v[6536] = 0e0;
v[6537] = 0e0;
v[6538] = 0e0;
v[6539] = 0e0;
v[6540] = -1e0;
v[6541] = 0e0;
v[6542] = 0e0;
v[6543] = 0e0;
v[6544] = 0e0;
v[6545] = 0e0;
v[1] = gti[0];
v[2] = gti[1];
v[3] = gti[2];
v[4] = (*epst);
v[5] = (*mus);
v[6] = (*mud);
v[7] = (*zetan);
v[8] = (*meq);
v[9] = (*ct);
v[10] = (*a4);
v[11] = (*a5);
v[12] = (*a6);
v[13] = invH[0][0];
v[14] = invH[0][1];
v[15] = invH[0][2];
v[16] = invH[0][3];
v[17] = invH[1][0];
v[18] = invH[1][1];
v[19] = invH[1][2];
v[20] = invH[1][3];
v[21] = invH[2][0];
v[22] = invH[2][1];
v[23] = invH[2][2];
v[24] = invH[2][3];
v[25] = invH[3][0];
v[26] = invH[3][1];
v[27] = invH[3][2];
v[28] = invH[3][3];
b29 = (*stick);
b30 = (*previouscontact);
v[31] = (*epsn1);
v[32] = (*gnb);
v[33] = (*gnbb);
v[730] = v[32] - v[33];
v[34] = (*n1);
v[4485] = v[31] * v[34];
v[745] = -1e0 + v[34];
v[2531] = -1e0 + v[745];
v[35] = (*n2);
v[750] = -1e0 + v[35];
v[2534] = -1e0 + v[750];
v[37] = x1A[0];
v[341] = -v[37] / 2e0;
v[38] = x1A[1];
v[343] = -v[38] / 2e0;
v[39] = x1A[2];
v[345] = -v[39] / 2e0;
v[40] = x2A[0];
v[334] = v[341] + v[40] / 2e0;
v[41] = x2A[1];
v[335] = v[343] + v[41] / 2e0;
v[42] = x2A[2];
v[336] = v[345] + v[42] / 2e0;
v[43] = x3A[0];
v[342] = v[341] + v[43] / 2e0;
v[44] = x3A[1];
v[344] = v[343] + v[44] / 2e0;
v[45] = x3A[2];
v[346] = v[345] + v[45] / 2e0;
v[46] = xAi[0];
v[47] = xAi[1];
v[48] = xAi[2];
v[49] = QAi[0][0];
v[50] = QAi[0][1];
v[51] = QAi[0][2];
v[52] = QAi[1][0];
v[53] = QAi[1][1];
v[54] = QAi[1][2];
v[55] = QAi[2][0];
v[56] = QAi[2][1];
v[57] = QAi[2][2];
v[58] = uA[0];
v[59] = uA[1];
v[60] = uA[2];
v[61] = alphaA[0];
v[370] = v[61] / 2e0;
v[368] = 2e0*v[61];
v[4914] = 4e0*v[368];
v[4731] = -v[368] / 2e0;
v[156] = (v[61] * v[61]);
v[62] = alphaA[1];
v[371] = 2e0*v[62];
v[6486] = 0e0;
v[6487] = 0e0;
v[6488] = 0e0;
v[6489] = v[368];
v[6490] = v[371];
v[6491] = 0e0;
v[6492] = 0e0;
v[6493] = 0e0;
v[6494] = 0e0;
v[6495] = 0e0;
v[6496] = 0e0;
v[6497] = 0e0;
v[4732] = -v[371] / 2e0;
v[6746] = 0e0;
v[6747] = 0e0;
v[6748] = 0e0;
v[6749] = -v[368];
v[6750] = -v[371];
v[6751] = 0e0;
v[6752] = 0e0;
v[6753] = 0e0;
v[6754] = 0e0;
v[6755] = 0e0;
v[6756] = 0e0;
v[6757] = 0e0;
v[369] = v[62] / 2e0;
v[6234] = 0e0;
v[6235] = 0e0;
v[6236] = 0e0;
v[6237] = v[369];
v[6238] = v[370];
v[6239] = 0e0;
v[6240] = 0e0;
v[6241] = 0e0;
v[6242] = 0e0;
v[6243] = 0e0;
v[6244] = 0e0;
v[6245] = 0e0;
v[154] = v[369] * v[61];
v[149] = (v[62] * v[62]);
v[420] = -v[149] - v[156];
v[4558] = v[420] / 2e0;
v[63] = alphaA[2];
v[398] = v[154] + v[63];
v[4570] = -2e0*v[398];
v[388] = v[154] - v[63];
v[4568] = -2e0*v[388];
v[373] = 2e0*v[63];
v[6402] = 0e0;
v[6403] = 0e0;
v[6404] = 0e0;
v[6405] = v[368];
v[6406] = 0e0;
v[6407] = v[373];
v[6408] = 0e0;
v[6409] = 0e0;
v[6410] = 0e0;
v[6411] = 0e0;
v[6412] = 0e0;
v[6413] = 0e0;
v[6318] = 0e0;
v[6319] = 0e0;
v[6320] = 0e0;
v[6321] = 0e0;
v[6322] = v[371];
v[6323] = v[373];
v[6324] = 0e0;
v[6325] = 0e0;
v[6326] = 0e0;
v[6327] = 0e0;
v[6328] = 0e0;
v[6329] = 0e0;
v[4730] = -v[373] / 2e0;
v[6794] = 0e0;
v[6795] = 0e0;
v[6796] = 0e0;
v[6797] = 0e0;
v[6798] = -v[371];
v[6799] = -v[373];
v[6800] = 0e0;
v[6801] = 0e0;
v[6802] = 0e0;
v[6803] = 0e0;
v[6804] = 0e0;
v[6805] = 0e0;
v[6222] = 0e0;
v[6223] = 0e0;
v[6224] = 0e0;
v[6225] = v[368];
v[6226] = v[371];
v[6227] = v[373];
v[6228] = 0e0;
v[6229] = 0e0;
v[6230] = 0e0;
v[6231] = 0e0;
v[6232] = 0e0;
v[6233] = 0e0;
v[372] = v[63] / 2e0;
v[6246] = 0e0;
v[6247] = 0e0;
v[6248] = 0e0;
v[6249] = v[372];
v[6250] = 0e0;
v[6251] = v[370];
v[6252] = 0e0;
v[6253] = 0e0;
v[6254] = 0e0;
v[6255] = 0e0;
v[6256] = 0e0;
v[6257] = 0e0;
v[6258] = 0e0;
v[6259] = 0e0;
v[6260] = 0e0;
v[6261] = 0e0;
v[6262] = v[372];
v[6263] = v[369];
v[6264] = 0e0;
v[6265] = 0e0;
v[6266] = 0e0;
v[6267] = 0e0;
v[6268] = 0e0;
v[6269] = 0e0;
v[161] = v[372] * v[62];
v[415] = v[161] + v[61];
v[4571] = -2e0*v[415];
v[407] = v[161] - v[61];
v[4569] = -2e0*v[407];
v[159] = v[372] * v[61];
v[411] = v[159] - v[62];
v[4572] = -2e0*v[411];
v[394] = v[159] + v[62];
v[4567] = -2e0*v[394];
v[150] = (v[63] * v[63]);
v[1126] = 4e0 + v[149] + v[150] + v[156];
v[4915] = 24e0 / Power(v[1126], 4);
v[1782] = 1e0 / Power(v[1126], 3);
v[4913] = -2e0*v[1782];
v[4425] = 8e0*v[1782];
v[1844] = -(v[368] * v[4425]);
v[1842] = v[371] * v[4425];
v[1839] = -(v[373] * v[4425]);
v[982] = 1e0 / (v[1126] * v[1126]);
v[4426] = 4e0*v[982];
v[4742] = 2e0*v[982];
v[4587] = 8e0*v[982];
v[402] = -v[150] - v[156];
v[4559] = v[402] / 2e0;
v[383] = -v[149] - v[150];
v[4560] = v[383] / 2e0;
v[382] = v[373] * v[4426];
v[4432] = -v[382] / 2e0;
v[424] = v[420] * v[4432];
v[381] = -(v[371] * v[4426]);
v[4429] = v[381] / 2e0;
v[404] = v[402] * v[4429];
v[380] = v[368] * v[4426];
v[4431] = -v[380] / 2e0;
v[384] = v[383] * v[4431];
v[256] = dalphaiA[0] * v[11] + ddalphaiA[0] * v[12] + v[10] * v[61];
v[262] = dalphaiA[1] * v[11] + ddalphaiA[1] * v[12] + v[10] * v[62];
v[264] = dalphaiA[2] * v[11] + ddalphaiA[2] * v[12] + v[10] * v[63];
v[133] = -cAp[0] / 2e0;
v[135] = -cAp[1] / 2e0;
v[138] = -cAi[0] / 2e0;
v[140] = -cAi[1] / 2e0;
v[80] = x1B[0];
v[358] = -v[80] / 2e0;
v[81] = x1B[1];
v[360] = -v[81] / 2e0;
v[82] = x1B[2];
v[362] = -v[82] / 2e0;
v[83] = x2B[0];
v[351] = v[358] + v[83] / 2e0;
v[84] = x2B[1];
v[352] = v[360] + v[84] / 2e0;
v[85] = x2B[2];
v[353] = v[362] + v[85] / 2e0;
v[86] = x3B[0];
v[359] = v[358] + v[86] / 2e0;
v[87] = x3B[1];
v[361] = v[360] + v[87] / 2e0;
v[88] = x3B[2];
v[363] = v[362] + v[88] / 2e0;
v[89] = xBi[0];
v[90] = xBi[1];
v[91] = xBi[2];
v[92] = QBi[0][0];
v[93] = QBi[0][1];
v[94] = QBi[0][2];
v[95] = QBi[1][0];
v[96] = QBi[1][1];
v[97] = QBi[1][2];
v[98] = QBi[2][0];
v[99] = QBi[2][1];
v[100] = QBi[2][2];
v[101] = uB[0];
v[102] = uB[1];
v[103] = uB[2];
v[104] = alphaB[0];
v[376] = v[104] / 2e0;
v[374] = 2e0*v[104];
v[4904] = 4e0*v[374];
v[4719] = -v[374] / 2e0;
v[212] = (v[104] * v[104]);
v[105] = alphaB[1];
v[377] = 2e0*v[105];
v[6522] = 0e0;
v[6523] = 0e0;
v[6524] = 0e0;
v[6525] = 0e0;
v[6526] = 0e0;
v[6527] = 0e0;
v[6528] = 0e0;
v[6529] = 0e0;
v[6530] = 0e0;
v[6531] = v[374];
v[6532] = v[377];
v[6533] = 0e0;
v[4720] = -v[377] / 2e0;
v[6818] = 0e0;
v[6819] = 0e0;
v[6820] = 0e0;
v[6821] = 0e0;
v[6822] = 0e0;
v[6823] = 0e0;
v[6824] = 0e0;
v[6825] = 0e0;
v[6826] = 0e0;
v[6827] = -v[374];
v[6828] = -v[377];
v[6829] = 0e0;
v[375] = v[105] / 2e0;
v[6282] = 0e0;
v[6283] = 0e0;
v[6284] = 0e0;
v[6285] = 0e0;
v[6286] = 0e0;
v[6287] = 0e0;
v[6288] = 0e0;
v[6289] = 0e0;
v[6290] = 0e0;
v[6291] = v[375];
v[6292] = v[376];
v[6293] = 0e0;
v[210] = v[104] * v[375];
v[205] = (v[105] * v[105]);
v[519] = -v[205] - v[212];
v[4537] = v[519] / 2e0;
v[106] = alphaB[2];
v[497] = v[106] + v[210];
v[4564] = -2e0*v[497];
v[487] = -v[106] + v[210];
v[4562] = -2e0*v[487];
v[379] = 2e0*v[106];
v[6438] = 0e0;
v[6439] = 0e0;
v[6440] = 0e0;
v[6441] = 0e0;
v[6442] = 0e0;
v[6443] = 0e0;
v[6444] = 0e0;
v[6445] = 0e0;
v[6446] = 0e0;
v[6447] = v[374];
v[6448] = 0e0;
v[6449] = v[379];
v[6354] = 0e0;
v[6355] = 0e0;
v[6356] = 0e0;
v[6357] = 0e0;
v[6358] = 0e0;
v[6359] = 0e0;
v[6360] = 0e0;
v[6361] = 0e0;
v[6362] = 0e0;
v[6363] = 0e0;
v[6364] = v[377];
v[6365] = v[379];
v[4718] = -v[379] / 2e0;
v[6866] = 0e0;
v[6867] = 0e0;
v[6868] = 0e0;
v[6869] = 0e0;
v[6870] = 0e0;
v[6871] = 0e0;
v[6872] = 0e0;
v[6873] = 0e0;
v[6874] = 0e0;
v[6875] = 0e0;
v[6876] = -v[377];
v[6877] = -v[379];
v[6270] = 0e0;
v[6271] = 0e0;
v[6272] = 0e0;
v[6273] = 0e0;
v[6274] = 0e0;
v[6275] = 0e0;
v[6276] = 0e0;
v[6277] = 0e0;
v[6278] = 0e0;
v[6279] = v[374];
v[6280] = v[377];
v[6281] = v[379];
v[378] = v[106] / 2e0;
v[6294] = 0e0;
v[6295] = 0e0;
v[6296] = 0e0;
v[6297] = 0e0;
v[6298] = 0e0;
v[6299] = 0e0;
v[6300] = 0e0;
v[6301] = 0e0;
v[6302] = 0e0;
v[6303] = v[378];
v[6304] = 0e0;
v[6305] = v[376];
v[6306] = 0e0;
v[6307] = 0e0;
v[6308] = 0e0;
v[6309] = 0e0;
v[6310] = 0e0;
v[6311] = 0e0;
v[6312] = 0e0;
v[6313] = 0e0;
v[6314] = 0e0;
v[6315] = 0e0;
v[6316] = v[378];
v[6317] = v[375];
v[217] = v[105] * v[378];
v[514] = v[104] + v[217];
v[4565] = -2e0*v[514];
v[506] = -v[104] + v[217];
v[4563] = -2e0*v[506];
v[215] = v[104] * v[378];
v[510] = -v[105] + v[215];
v[4566] = -2e0*v[510];
v[493] = v[105] + v[215];
v[4561] = -2e0*v[493];
v[206] = (v[106] * v[106]);
v[1133] = 4e0 + v[205] + v[206] + v[212];
v[4905] = 24e0 / Power(v[1133], 4);
v[1790] = 1e0 / Power(v[1133], 3);
v[4903] = -2e0*v[1790];
v[4427] = 8e0*v[1790];
v[1870] = -(v[374] * v[4427]);
v[1868] = v[377] * v[4427];
v[1865] = -(v[379] * v[4427]);
v[988] = 1e0 / (v[1133] * v[1133]);
v[4428] = 4e0*v[988];
v[4756] = 2e0*v[988];
v[4588] = 8e0*v[988];
v[501] = -v[206] - v[212];
v[4538] = v[501] / 2e0;
v[482] = -v[205] - v[206];
v[4539] = v[482] / 2e0;
v[481] = v[379] * v[4428];
v[4439] = -v[481] / 2e0;
v[523] = v[4439] * v[519];
v[480] = -(v[377] * v[4428]);
v[4436] = v[480] / 2e0;
v[503] = v[4436] * v[501];
v[479] = v[374] * v[4428];
v[4438] = -v[479] / 2e0;
v[483] = v[4438] * v[482];
v[282] = v[10] * v[104] + dalphaiB[0] * v[11] + ddalphaiB[0] * v[12];
v[288] = v[10] * v[105] + dalphaiB[1] * v[11] + ddalphaiB[1] * v[12];
v[290] = v[10] * v[106] + dalphaiB[2] * v[11] + ddalphaiB[2] * v[12];
v[189] = -cBp[0] / 2e0;
v[191] = -cBp[1] / 2e0;
v[194] = -cBi[0] / 2e0;
v[196] = -cBi[1] / 2e0;
v[132] = v[133] + v[135];
v[134] = 0.5e0 - v[133];
v[136] = 0.5e0 - v[135];
v[137] = v[138] + v[140];
v[139] = 0.5e0 - v[138];
v[141] = 0.5e0 - v[140];
v[142] = v[132] * v[37] + v[134] * v[40] + v[136] * v[43];
v[143] = v[132] * v[38] + v[134] * v[41] + v[136] * v[44];
v[144] = v[132] * v[39] + v[134] * v[42] + v[136] * v[45];
v[881] = v[142] * v[49] + v[143] * v[50] + v[144] * v[51];
v[879] = v[142] * v[52] + v[143] * v[53] + v[144] * v[54];
v[877] = v[142] * v[55] + v[143] * v[56] + v[144] * v[57];
v[145] = v[137] * v[37] + v[139] * v[40] + v[141] * v[43];
v[146] = v[137] * v[38] + v[139] * v[41] + v[141] * v[44];
v[147] = v[137] * v[39] + v[139] * v[42] + v[141] * v[45];
v[148] = 4e0 / v[1126];
v[4585] = 2e0*v[148];
v[4430] = -v[148] / 2e0;
v[422] = v[371] * v[4430];
v[423] = v[422] + v[420] * v[4429];
v[419] = v[368] * v[4430];
v[421] = v[419] + v[420] * v[4431];
v[416] = v[148] - v[380] * v[415];
v[413] = -v[148] + v[381] * v[411];
v[408] = -v[148] - v[380] * v[407];
v[405] = v[373] * v[4430];
v[406] = v[405] + v[402] * v[4432];
v[403] = v[419] + v[402] * v[4431];
v[401] = v[148] - v[382] * v[398];
v[396] = v[148] + v[381] * v[394];
v[393] = -(v[148] * v[372]);
v[4590] = 2e0*v[393];
v[417] = -v[393] + v[381] * v[415];
v[462] = v[413] * v[51] + v[417] * v[54] + v[423] * v[57];
v[459] = v[413] * v[50] + v[417] * v[53] + v[423] * v[56];
v[456] = v[413] * v[49] + v[417] * v[52] + v[423] * v[55];
v[477] = v[142] * v[456] + v[143] * v[459] + v[144] * v[462];
v[412] = -v[393] - v[380] * v[411];
v[461] = v[412] * v[51] + v[416] * v[54] + v[421] * v[57];
v[458] = v[412] * v[50] + v[416] * v[53] + v[421] * v[56];
v[455] = v[412] * v[49] + v[416] * v[52] + v[421] * v[55];
v[476] = v[142] * v[455] + v[143] * v[458] + v[144] * v[461];
v[409] = -v[393] + v[381] * v[407];
v[395] = -v[393] - v[380] * v[394];
v[392] = -v[148] - v[382] * v[388];
v[390] = -(v[148] * v[370]);
v[4591] = 2e0*v[390];
v[414] = -v[390] - v[382] * v[411];
v[400] = -v[390] + v[381] * v[398];
v[447] = v[400] * v[51] + v[404] * v[54] + v[409] * v[57];
v[444] = v[400] * v[50] + v[404] * v[53] + v[409] * v[56];
v[441] = v[400] * v[49] + v[404] * v[52] + v[409] * v[55];
v[474] = v[142] * v[441] + v[143] * v[444] + v[144] * v[447];
v[397] = -v[390] - v[382] * v[394];
v[391] = v[381] * v[388] - v[390];
v[387] = v[148] * v[369];
v[4592] = 2e0*v[387];
v[418] = v[387] - v[382] * v[415];
v[463] = v[414] * v[51] + v[418] * v[54] + v[424] * v[57];
v[460] = v[414] * v[50] + v[418] * v[53] + v[424] * v[56];
v[457] = v[414] * v[49] + v[418] * v[52] + v[424] * v[55];
v[478] = v[142] * v[457] + v[143] * v[460] + v[144] * v[463];
v[410] = v[387] - v[382] * v[407];
v[448] = v[401] * v[51] + v[406] * v[54] + v[410] * v[57];
v[445] = v[401] * v[50] + v[406] * v[53] + v[410] * v[56];
v[442] = v[401] * v[49] + v[406] * v[52] + v[410] * v[55];
v[475] = v[142] * v[442] + v[143] * v[445] + v[144] * v[448];
v[399] = v[387] - v[380] * v[398];
v[446] = v[399] * v[51] + v[403] * v[54] + v[408] * v[57];
v[443] = v[399] * v[50] + v[403] * v[53] + v[408] * v[56];
v[440] = v[399] * v[49] + v[403] * v[52] + v[408] * v[55];
v[473] = v[142] * v[440] + v[143] * v[443] + v[144] * v[446];
v[389] = v[387] - v[380] * v[388];
v[431] = v[384] * v[51] + v[389] * v[54] + v[395] * v[57];
v[428] = v[384] * v[50] + v[389] * v[53] + v[395] * v[56];
v[425] = v[384] * v[49] + v[389] * v[52] + v[395] * v[55];
v[470] = v[142] * v[425] + v[143] * v[428] + v[144] * v[431];
v[386] = v[405] + v[383] * v[4432];
v[433] = v[386] * v[51] + v[392] * v[54] + v[397] * v[57];
v[430] = v[386] * v[50] + v[392] * v[53] + v[397] * v[56];
v[427] = v[386] * v[49] + v[392] * v[52] + v[397] * v[55];
v[472] = v[142] * v[427] + v[143] * v[430] + v[144] * v[433];
v[385] = v[422] + v[383] * v[4429];
v[432] = v[385] * v[51] + v[391] * v[54] + v[396] * v[57];
v[429] = v[385] * v[50] + v[391] * v[53] + v[396] * v[56];
v[426] = v[385] * v[49] + v[391] * v[52] + v[396] * v[55];
v[471] = v[142] * v[426] + v[143] * v[429] + v[144] * v[432];
v[253] = (v[148] * v[148]);
v[4435] = v[264] / v[253];
v[4434] = v[262] / v[253];
v[4433] = v[256] / v[253];
v[4916] = 2e0 / Power(v[253], 3);
v[151] = 1e0 - v[383] * v[4430];
v[4881] = v[151] / v[253];
v[1581] = v[151] * v[4433];
v[152] = v[148] * v[388];
v[4555] = v[152] * v[262];
v[1583] = v[152] * v[4434];
v[2254] = v[1581] + v[1583];
v[153] = v[148] * v[394];
v[4929] = v[153] / v[253];
v[4554] = v[153] * v[264];
v[1584] = v[153] * v[4435];
v[2259] = v[1581] + v[1584];
v[2248] = v[1583] + v[2259];
v[155] = v[148] * v[398];
v[1588] = v[155] * v[4433];
v[157] = 1e0 - v[402] * v[4430];
v[4880] = v[157] / v[253];
v[1577] = v[157] * v[4434];
v[2250] = v[1577] + v[1588];
v[158] = v[148] * v[407];
v[4930] = v[158] / v[253];
v[1578] = v[158] * v[4435];
v[2260] = v[1577] + v[1578];
v[2253] = v[1588] + v[2260];
v[160] = v[148] * v[411];
v[1591] = v[160] * v[4433];
v[162] = v[148] * v[415];
v[4925] = v[162] / v[253];
v[1573] = v[162] * v[4434];
v[163] = 1e0 - v[420] * v[4430];
v[4884] = v[163] / v[253];
v[1574] = v[163] * v[4435];
v[4882] = v[1574] * v[253];
v[2258] = v[1573] + v[1574] + v[1591];
v[2255] = -v[1591] + v[2258];
v[2249] = -v[1573] + v[2258];
v[271] = -(v[387] * v[393]);
v[4450] = v[271] - v[380];
v[269] = v[390] * v[393];
v[4448] = v[269] - v[381];
v[260] = -(v[387] * v[390]);
v[4445] = v[260] - v[382];
v[167] = v[151] * v[49] + v[152] * v[52] + v[153] * v[55];
v[168] = v[151] * v[50] + v[152] * v[53] + v[153] * v[56];
v[169] = v[151] * v[51] + v[152] * v[54] + v[153] * v[57];
v[347] = v[167] * v[342] + v[168] * v[344] + v[169] * v[346];
v[337] = v[167] * v[334] + v[168] * v[335] + v[169] * v[336];
v[170] = v[155] * v[49] + v[157] * v[52] + v[158] * v[55];
v[171] = v[155] * v[50] + v[157] * v[53] + v[158] * v[56];
v[172] = v[155] * v[51] + v[157] * v[54] + v[158] * v[57];
v[348] = v[170] * v[342] + v[171] * v[344] + v[172] * v[346];
v[338] = v[170] * v[334] + v[171] * v[335] + v[172] * v[336];
v[173] = v[160] * v[49] + v[162] * v[52] + v[163] * v[55];
v[174] = v[160] * v[50] + v[162] * v[53] + v[163] * v[56];
v[175] = v[160] * v[51] + v[162] * v[54] + v[163] * v[57];
v[349] = v[173] * v[342] + v[174] * v[344] + v[175] * v[346];
v[339] = v[173] * v[334] + v[174] * v[335] + v[175] * v[336];
v[188] = v[189] + v[191];
v[190] = 0.5e0 - v[189];
v[192] = 0.5e0 - v[191];
v[193] = v[194] + v[196];
v[195] = 0.5e0 - v[194];
v[197] = 0.5e0 - v[196];
v[198] = v[188] * v[80] + v[190] * v[83] + v[192] * v[86];
v[199] = v[188] * v[81] + v[190] * v[84] + v[192] * v[87];
v[200] = v[188] * v[82] + v[190] * v[85] + v[192] * v[88];
v[875] = -(v[198] * v[92]) - v[199] * v[93] - v[200] * v[94];
v[873] = -(v[198] * v[95]) - v[199] * v[96] - v[200] * v[97];
v[871] = -(v[100] * v[200]) - v[198] * v[98] - v[199] * v[99];
v[201] = v[193] * v[80] + v[195] * v[83] + v[197] * v[86];
v[202] = v[193] * v[81] + v[195] * v[84] + v[197] * v[87];
v[203] = v[193] * v[82] + v[195] * v[85] + v[197] * v[88];
v[204] = 4e0 / v[1133];
v[4586] = 2e0*v[204];
v[4437] = -v[204] / 2e0;
v[521] = v[377] * v[4437];
v[522] = v[4436] * v[519] + v[521];
v[518] = v[374] * v[4437];
v[520] = v[518] + v[4438] * v[519];
v[515] = v[204] - v[479] * v[514];
v[512] = -v[204] + v[480] * v[510];
v[507] = -v[204] - v[479] * v[506];
v[504] = v[379] * v[4437];
v[505] = v[4439] * v[501] + v[504];
v[502] = v[4438] * v[501] + v[518];
v[500] = v[204] - v[481] * v[497];
v[495] = v[204] + v[480] * v[493];
v[492] = -(v[204] * v[378]);
v[4594] = 2e0*v[492];
v[516] = -v[492] + v[480] * v[514];
v[561] = v[100] * v[522] + v[512] * v[94] + v[516] * v[97];
v[558] = v[512] * v[93] + v[516] * v[96] + v[522] * v[99];
v[555] = v[512] * v[92] + v[516] * v[95] + v[522] * v[98];
v[576] = v[198] * v[555] + v[199] * v[558] + v[200] * v[561];
v[511] = -v[492] - v[479] * v[510];
v[560] = v[100] * v[520] + v[511] * v[94] + v[515] * v[97];
v[557] = v[511] * v[93] + v[515] * v[96] + v[520] * v[99];
v[554] = v[511] * v[92] + v[515] * v[95] + v[520] * v[98];
v[575] = v[198] * v[554] + v[199] * v[557] + v[200] * v[560];
v[508] = -v[492] + v[480] * v[506];
v[494] = -v[492] - v[479] * v[493];
v[491] = -v[204] - v[481] * v[487];
v[489] = -(v[204] * v[376]);
v[4595] = 2e0*v[489];
v[513] = -v[489] - v[481] * v[510];
v[499] = -v[489] + v[480] * v[497];
v[546] = v[100] * v[508] + v[499] * v[94] + v[503] * v[97];
v[543] = v[499] * v[93] + v[503] * v[96] + v[508] * v[99];
v[540] = v[499] * v[92] + v[503] * v[95] + v[508] * v[98];
v[573] = v[198] * v[540] + v[199] * v[543] + v[200] * v[546];
v[496] = -v[489] - v[481] * v[493];
v[490] = v[480] * v[487] - v[489];
v[486] = v[204] * v[375];
v[4596] = 2e0*v[486];
v[517] = v[486] - v[481] * v[514];
v[562] = v[100] * v[523] + v[513] * v[94] + v[517] * v[97];
v[559] = v[513] * v[93] + v[517] * v[96] + v[523] * v[99];
v[556] = v[513] * v[92] + v[517] * v[95] + v[523] * v[98];
v[577] = v[198] * v[556] + v[199] * v[559] + v[200] * v[562];
v[509] = v[486] - v[481] * v[506];
v[547] = v[100] * v[509] + v[500] * v[94] + v[505] * v[97];
v[544] = v[500] * v[93] + v[505] * v[96] + v[509] * v[99];
v[541] = v[500] * v[92] + v[505] * v[95] + v[509] * v[98];
v[574] = v[198] * v[541] + v[199] * v[544] + v[200] * v[547];
v[498] = v[486] - v[479] * v[497];
v[545] = v[100] * v[507] + v[498] * v[94] + v[502] * v[97];
v[542] = v[498] * v[93] + v[502] * v[96] + v[507] * v[99];
v[539] = v[498] * v[92] + v[502] * v[95] + v[507] * v[98];
v[572] = v[198] * v[539] + v[199] * v[542] + v[200] * v[545];
v[488] = v[486] - v[479] * v[487];
v[530] = v[100] * v[494] + v[483] * v[94] + v[488] * v[97];
v[527] = v[483] * v[93] + v[488] * v[96] + v[494] * v[99];
v[524] = v[483] * v[92] + v[488] * v[95] + v[494] * v[98];
v[569] = v[198] * v[524] + v[199] * v[527] + v[200] * v[530];
v[587] = -(v[347] * v[569]) - v[348] * v[572] - v[349] * v[575];
v[581] = -(v[337] * v[569]) - v[338] * v[572] - v[339] * v[575];
v[485] = v[4439] * v[482] + v[504];
v[532] = v[100] * v[496] + v[485] * v[94] + v[491] * v[97];
v[529] = v[485] * v[93] + v[491] * v[96] + v[496] * v[99];
v[526] = v[485] * v[92] + v[491] * v[95] + v[496] * v[98];
v[571] = v[198] * v[526] + v[199] * v[529] + v[200] * v[532];
v[589] = -(v[347] * v[571]) - v[348] * v[574] - v[349] * v[577];
v[583] = -(v[337] * v[571]) - v[338] * v[574] - v[339] * v[577];
v[484] = v[4436] * v[482] + v[521];
v[531] = v[100] * v[495] + v[484] * v[94] + v[490] * v[97];
v[528] = v[484] * v[93] + v[490] * v[96] + v[495] * v[99];
v[525] = v[484] * v[92] + v[490] * v[95] + v[495] * v[98];
v[570] = v[198] * v[525] + v[199] * v[528] + v[200] * v[531];
v[588] = -(v[347] * v[570]) - v[348] * v[573] - v[349] * v[576];
v[582] = -(v[337] * v[570]) - v[338] * v[573] - v[339] * v[576];
v[279] = (v[204] * v[204]);
v[4442] = v[290] / v[279];
v[4441] = v[288] / v[279];
v[4440] = v[282] / v[279];
v[4906] = 2e0 / Power(v[279], 3);
v[207] = 1e0 - v[4437] * v[482];
v[4872] = v[207] / v[279];
v[1557] = v[207] * v[4440];
v[208] = v[204] * v[487];
v[4534] = v[208] * v[288];
v[1559] = v[208] * v[4441];
v[2269] = v[1557] + v[1559];
v[209] = v[204] * v[493];
v[4941] = v[209] / v[279];
v[4533] = v[209] * v[290];
v[1560] = v[209] * v[4442];
v[2274] = v[1557] + v[1560];
v[2263] = v[1559] + v[2274];
v[211] = v[204] * v[497];
v[1564] = v[211] * v[4440];
v[213] = 1e0 - v[4437] * v[501];
v[4871] = v[213] / v[279];
v[1553] = v[213] * v[4441];
v[2265] = v[1553] + v[1564];
v[214] = v[204] * v[506];
v[4942] = v[214] / v[279];
v[1554] = v[214] * v[4442];
v[2275] = v[1553] + v[1554];
v[2268] = v[1564] + v[2275];
v[216] = v[204] * v[510];
v[1567] = v[216] * v[4440];
v[218] = v[204] * v[514];
v[4937] = v[218] / v[279];
v[1549] = v[218] * v[4441];
v[219] = 1e0 - v[4437] * v[519];
v[4875] = v[219] / v[279];
v[1550] = v[219] * v[4442];
v[4873] = v[1550] * v[279];
v[2273] = v[1549] + v[1550] + v[1567];
v[2270] = -v[1567] + v[2273];
v[2264] = -v[1549] + v[2273];
v[297] = -(v[486] * v[492]);
v[4458] = v[297] - v[479];
v[295] = v[489] * v[492];
v[4456] = v[295] - v[480];
v[286] = -(v[486] * v[489]);
v[4453] = v[286] - v[481];
v[223] = v[207] * v[92] + v[208] * v[95] + v[209] * v[98];
v[224] = v[207] * v[93] + v[208] * v[96] + v[209] * v[99];
v[225] = v[100] * v[209] + v[207] * v[94] + v[208] * v[97];
v[364] = v[223] * v[359] + v[224] * v[361] + v[225] * v[363];
v[354] = v[223] * v[351] + v[224] * v[352] + v[225] * v[353];
v[669] = -(v[25] * v[337]) - v[26] * v[347] + v[27] * v[354] + v[28] * v[364];
v[654] = -(v[21] * v[337]) - v[22] * v[347] + v[23] * v[354] + v[24] * v[364];
v[639] = -(v[17] * v[337]) - v[18] * v[347] + v[19] * v[354] + v[20] * v[364];
v[624] = -(v[13] * v[337]) - v[14] * v[347] + v[15] * v[354] + v[16] * v[364];
v[226] = v[211] * v[92] + v[213] * v[95] + v[214] * v[98];
v[227] = v[211] * v[93] + v[213] * v[96] + v[214] * v[99];
v[228] = v[100] * v[214] + v[211] * v[94] + v[213] * v[97];
v[365] = v[226] * v[359] + v[227] * v[361] + v[228] * v[363];
v[355] = v[226] * v[351] + v[227] * v[352] + v[228] * v[353];
v[671] = -(v[25] * v[338]) - v[26] * v[348] + v[27] * v[355] + v[28] * v[365];
v[656] = -(v[21] * v[338]) - v[22] * v[348] + v[23] * v[355] + v[24] * v[365];
v[641] = -(v[17] * v[338]) - v[18] * v[348] + v[19] * v[355] + v[20] * v[365];
v[626] = -(v[13] * v[338]) - v[14] * v[348] + v[15] * v[355] + v[16] * v[365];
v[229] = v[216] * v[92] + v[218] * v[95] + v[219] * v[98];
v[230] = v[216] * v[93] + v[218] * v[96] + v[219] * v[99];
v[231] = v[100] * v[219] + v[216] * v[94] + v[218] * v[97];
v[366] = v[229] * v[359] + v[230] * v[361] + v[231] * v[363];
v[598] = -(v[364] * v[472]) - v[365] * v[475] - v[366] * v[478];
v[597] = -(v[364] * v[471]) - v[365] * v[474] - v[366] * v[477];
v[596] = -(v[364] * v[470]) - v[365] * v[473] - v[366] * v[476];
v[356] = v[229] * v[351] + v[230] * v[352] + v[231] * v[353];
v[673] = -(v[25] * v[339]) - v[26] * v[349] + v[27] * v[356] + v[28] * v[366];
v[658] = -(v[21] * v[339]) - v[22] * v[349] + v[23] * v[356] + v[24] * v[366];
v[643] = -(v[17] * v[339]) - v[18] * v[349] + v[19] * v[356] + v[20] * v[366];
v[628] = -(v[13] * v[339]) - v[14] * v[349] + v[15] * v[356] + v[16] * v[366];
v[592] = -(v[354] * v[472]) - v[355] * v[475] - v[356] * v[478];
v[591] = -(v[354] * v[471]) - v[355] * v[474] - v[356] * v[477];
v[590] = -(v[354] * v[470]) - v[355] * v[473] - v[356] * v[476];
v[4480] = -v[101] + v[46] + v[58] - v[89];
v[4479] = -v[102] + v[47] + v[59] - v[90];
v[4478] = -v[103] + v[48] + v[60] - v[91];
v[244] = duiA[0] * v[11] + dduiA[0] * v[12] + v[10] * v[58];
v[245] = duiA[1] * v[11] + dduiA[1] * v[12] + v[10] * v[59];
v[246] = duiA[2] * v[11] + dduiA[2] * v[12] + v[10] * v[60];
v[247] = v[10] * v[101] + duiB[0] * v[11] + dduiB[0] * v[12];
v[3283] = v[244] - v[247];
v[248] = v[10] * v[102] + duiB[1] * v[11] + dduiB[1] * v[12];
v[3262] = v[245] - v[248];
v[249] = v[10] * v[103] + duiB[2] * v[11] + dduiB[2] * v[12];
v[3241] = v[246] - v[249];
v[251] = v[269] + v[381];
v[4449] = v[251] / v[253];
v[4443] = v[163] * v[251];
v[252] = v[260] + v[382];
v[4446] = v[252] / v[253];
v[4444] = v[157] * v[252];
v[254] = v[253] + (v[390] * v[390]);
v[4401] = -(v[264] * v[4443]) - v[262] * v[4444] - v[254] * (v[256] + v[4554] + v[4555]);
v[3195] = v[254] / v[253];
v[1787] = v[4448] / v[253];
v[1785] = v[4445] / v[253];
v[7666] = 0e0;
v[7667] = 0e0;
v[7668] = 0e0;
v[7669] = 0e0;
v[7670] = v[1785];
v[7671] = v[1787];
v[7672] = 0e0;
v[7673] = 0e0;
v[7674] = 0e0;
v[7675] = 0e0;
v[7676] = 0e0;
v[7677] = 0e0;
v[7118] = 0e0;
v[7119] = 0e0;
v[7120] = 0e0;
v[7121] = 0e0;
v[7122] = v[10] * v[1785];
v[7123] = v[10] * v[1787];
v[7124] = 0e0;
v[7125] = 0e0;
v[7126] = 0e0;
v[7127] = 0e0;
v[7128] = 0e0;
v[7129] = 0e0;
v[255] = v[1785] * v[262] + v[1787] * v[264] + v[256] * v[3195];
v[257] = v[253] + (v[387] * v[387]);
v[4543] = v[155] * v[257] + v[151] * v[4445];
v[263] = v[271] + v[380];
v[4447] = v[158] * v[257] + v[163] * v[263];
v[4400] = -(v[257] * v[262]) - v[264] * v[4447] - v[256] * v[4543];
v[3198] = v[257] / v[253];
v[1786] = v[4450] / v[253];
v[7678] = 0e0;
v[7679] = 0e0;
v[7680] = 0e0;
v[7681] = v[4446];
v[7682] = 0e0;
v[7683] = v[1786];
v[7684] = 0e0;
v[7685] = 0e0;
v[7686] = 0e0;
v[7687] = 0e0;
v[7688] = 0e0;
v[7689] = 0e0;
v[7142] = 0e0;
v[7143] = 0e0;
v[7144] = 0e0;
v[7145] = v[10] * v[4446];
v[7146] = 0e0;
v[7147] = v[10] * v[1786];
v[7148] = 0e0;
v[7149] = 0e0;
v[7150] = 0e0;
v[7151] = 0e0;
v[7152] = 0e0;
v[7153] = 0e0;
v[265] = v[1786] * v[264] + v[262] * v[3198] + v[256] * v[4446];
v[266] = v[253] + (v[393] * v[393]);
v[4544] = v[160] * v[266] + v[151] * v[4448];
v[4545] = v[162] * v[266] + v[157] * v[4450];
v[4399] = -(v[264] * v[266]) - v[256] * v[4544] - v[262] * v[4545];
v[3200] = v[266] / v[253];
v[1784] = v[263] / v[253];
v[7690] = 0e0;
v[7691] = 0e0;
v[7692] = 0e0;
v[7693] = v[4449];
v[7694] = v[1784];
v[7695] = 0e0;
v[7696] = 0e0;
v[7697] = 0e0;
v[7698] = 0e0;
v[7699] = 0e0;
v[7700] = 0e0;
v[7701] = 0e0;
v[7166] = 0e0;
v[7167] = 0e0;
v[7168] = 0e0;
v[7169] = v[10] * v[4449];
v[7170] = v[10] * v[1784];
v[7171] = 0e0;
v[7172] = 0e0;
v[7173] = 0e0;
v[7174] = 0e0;
v[7175] = 0e0;
v[7176] = 0e0;
v[7177] = 0e0;
v[275] = v[1784] * v[262] + v[264] * v[3200] + v[256] * v[4449];
v[277] = v[295] + v[480];
v[4457] = v[277] / v[279];
v[4451] = v[219] * v[277];
v[278] = v[286] + v[481];
v[4454] = v[278] / v[279];
v[4452] = v[213] * v[278];
v[280] = v[279] + (v[489] * v[489]);
v[4376] = -(v[290] * v[4451]) - v[288] * v[4452] - v[280] * (v[282] + v[4533] + v[4534]);
v[3202] = v[280] / v[279];
v[1795] = v[4456] / v[279];
v[1793] = v[4453] / v[279];
v[7702] = 0e0;
v[7703] = 0e0;
v[7704] = 0e0;
v[7705] = 0e0;
v[7706] = 0e0;
v[7707] = 0e0;
v[7708] = 0e0;
v[7709] = 0e0;
v[7710] = 0e0;
v[7711] = 0e0;
v[7712] = v[1793];
v[7713] = v[1795];
v[7190] = 0e0;
v[7191] = 0e0;
v[7192] = 0e0;
v[7193] = 0e0;
v[7194] = 0e0;
v[7195] = 0e0;
v[7196] = 0e0;
v[7197] = 0e0;
v[7198] = 0e0;
v[7199] = 0e0;
v[7200] = v[10] * v[1793];
v[7201] = v[10] * v[1795];
v[281] = v[1793] * v[288] + v[1795] * v[290] + v[282] * v[3202];
v[283] = v[279] + (v[486] * v[486]);
v[4522] = v[211] * v[283] + v[207] * v[4453];
v[289] = v[297] + v[479];
v[4455] = v[214] * v[283] + v[219] * v[289];
v[4375] = -(v[283] * v[288]) - v[290] * v[4455] - v[282] * v[4522];
v[3205] = v[283] / v[279];
v[1794] = v[4458] / v[279];
v[7714] = 0e0;
v[7715] = 0e0;
v[7716] = 0e0;
v[7717] = 0e0;
v[7718] = 0e0;
v[7719] = 0e0;
v[7720] = 0e0;
v[7721] = 0e0;
v[7722] = 0e0;
v[7723] = v[4454];
v[7724] = 0e0;
v[7725] = v[1794];
v[7214] = 0e0;
v[7215] = 0e0;
v[7216] = 0e0;
v[7217] = 0e0;
v[7218] = 0e0;
v[7219] = 0e0;
v[7220] = 0e0;
v[7221] = 0e0;
v[7222] = 0e0;
v[7223] = v[10] * v[4454];
v[7224] = 0e0;
v[7225] = v[10] * v[1794];
v[291] = v[1794] * v[290] + v[288] * v[3205] + v[282] * v[4454];
v[292] = v[279] + (v[492] * v[492]);
v[4523] = v[216] * v[292] + v[207] * v[4456];
v[4524] = v[218] * v[292] + v[213] * v[4458];
v[4374] = -(v[290] * v[292]) - v[282] * v[4523] - v[288] * v[4524];
v[3207] = v[292] / v[279];
v[1792] = v[289] / v[279];
v[7726] = 0e0;
v[7727] = 0e0;
v[7728] = 0e0;
v[7729] = 0e0;
v[7730] = 0e0;
v[7731] = 0e0;
v[7732] = 0e0;
v[7733] = 0e0;
v[7734] = 0e0;
v[7735] = v[4457];
v[7736] = v[1792];
v[7737] = 0e0;
v[7238] = 0e0;
v[7239] = 0e0;
v[7240] = 0e0;
v[7241] = 0e0;
v[7242] = 0e0;
v[7243] = 0e0;
v[7244] = 0e0;
v[7245] = 0e0;
v[7246] = 0e0;
v[7247] = v[10] * v[4457];
v[7248] = v[10] * v[1792];
v[7249] = 0e0;
v[301] = v[1792] * v[288] + v[290] * v[3207] + v[282] * v[4457];
v[2344] = v[3283] + v[255] * v[470] + v[265] * v[471] + v[275] * v[472] - v[281] * v[569] - v[291] * v[570]
- v[301] * v[571];
v[4459] = v[2344] - v[3283];
v[2299] = v[3262] + v[255] * v[473] + v[265] * v[474] + v[275] * v[475] - v[281] * v[572] - v[291] * v[573]
- v[301] * v[574];
v[4460] = v[2299] - v[3262];
v[2279] = v[3241] + v[255] * v[476] + v[265] * v[477] + v[275] * v[478] - v[281] * v[575] - v[291] * v[576]
- v[301] * v[577];
v[4461] = v[2279] - v[3241];
v[302] = v[142] * v[167] + v[143] * v[168] + v[144] * v[169] - v[198] * v[223] - v[199] * v[224] - v[200] * v[225] + v[4480];
v[303] = v[142] * v[170] + v[143] * v[171] + v[144] * v[172] - v[198] * v[226] - v[199] * v[227] - v[200] * v[228] + v[4479];
v[304] = v[142] * v[173] + v[143] * v[174] + v[144] * v[175] - v[198] * v[229] - v[199] * v[230] - v[200] * v[231] + v[4478];
v[732] = sqrt((v[302] * v[302]) + (v[303] * v[303]) + (v[304] * v[304]));
v[1014] = 1e0 / (v[732] * v[732]);
v[601] = -(v[302] * (v[359] * v[526] + v[361] * v[529] + v[363] * v[532])) - v[303] * (v[359] * v[541] + v[361] * v[544]
	+ v[363] * v[547]) - v[304] * (v[359] * v[556] + v[361] * v[559] + v[363] * v[562]) + v[364] * v[571] + v[365] * v[574]
	+ v[366] * v[577];
v[600] = -(v[302] * (v[359] * v[525] + v[361] * v[528] + v[363] * v[531])) - v[303] * (v[359] * v[540] + v[361] * v[543]
	+ v[363] * v[546]) - v[304] * (v[359] * v[555] + v[361] * v[558] + v[363] * v[561]) + v[364] * v[570] + v[365] * v[573]
	+ v[366] * v[576];
v[599] = -(v[302] * (v[359] * v[524] + v[361] * v[527] + v[363] * v[530])) - v[303] * (v[359] * v[539] + v[361] * v[542]
	+ v[363] * v[545]) - v[304] * (v[359] * v[554] + v[361] * v[557] + v[363] * v[560]) + v[364] * v[569] + v[365] * v[572]
	+ v[366] * v[575];
v[595] = -(v[302] * (v[351] * v[526] + v[352] * v[529] + v[353] * v[532])) - v[303] * (v[351] * v[541] + v[352] * v[544]
	+ v[353] * v[547]) - v[304] * (v[351] * v[556] + v[352] * v[559] + v[353] * v[562]) + v[354] * v[571] + v[355] * v[574]
	+ v[356] * v[577];
v[594] = -(v[302] * (v[351] * v[525] + v[352] * v[528] + v[353] * v[531])) - v[303] * (v[351] * v[540] + v[352] * v[543]
	+ v[353] * v[546]) - v[304] * (v[351] * v[555] + v[352] * v[558] + v[353] * v[561]) + v[354] * v[570] + v[355] * v[573]
	+ v[356] * v[576];
v[593] = -(v[302] * (v[351] * v[524] + v[352] * v[527] + v[353] * v[530])) - v[303] * (v[351] * v[539] + v[352] * v[542]
	+ v[353] * v[545]) - v[304] * (v[351] * v[554] + v[352] * v[557] + v[353] * v[560]) + v[354] * v[569] + v[355] * v[572]
	+ v[356] * v[575];
v[586] = v[302] * (v[342] * v[427] + v[344] * v[430] + v[346] * v[433]) + v[303] * (v[342] * v[442] + v[344] * v[445]
	+ v[346] * v[448]) + v[304] * (v[342] * v[457] + v[344] * v[460] + v[346] * v[463]) + v[347] * v[472] + v[348] * v[475]
	+ v[349] * v[478];
v[585] = v[302] * (v[342] * v[426] + v[344] * v[429] + v[346] * v[432]) + v[303] * (v[342] * v[441] + v[344] * v[444]
	+ v[346] * v[447]) + v[304] * (v[342] * v[456] + v[344] * v[459] + v[346] * v[462]) + v[347] * v[471] + v[348] * v[474]
	+ v[349] * v[477];
v[584] = v[302] * (v[342] * v[425] + v[344] * v[428] + v[346] * v[431]) + v[303] * (v[342] * v[440] + v[344] * v[443]
	+ v[346] * v[446]) + v[304] * (v[342] * v[455] + v[344] * v[458] + v[346] * v[461]) + v[347] * v[470] + v[348] * v[473]
	+ v[349] * v[476];
v[580] = v[302] * (v[334] * v[427] + v[335] * v[430] + v[336] * v[433]) + v[303] * (v[334] * v[442] + v[335] * v[445]
	+ v[336] * v[448]) + v[304] * (v[334] * v[457] + v[335] * v[460] + v[336] * v[463]) + v[337] * v[472] + v[338] * v[475]
	+ v[339] * v[478];
v[579] = v[302] * (v[334] * v[426] + v[335] * v[429] + v[336] * v[432]) + v[303] * (v[334] * v[441] + v[335] * v[444]
	+ v[336] * v[447]) + v[304] * (v[334] * v[456] + v[335] * v[459] + v[336] * v[462]) + v[337] * v[471] + v[338] * v[474]
	+ v[339] * v[477];
v[578] = v[302] * (v[334] * v[425] + v[335] * v[428] + v[336] * v[431]) + v[303] * (v[334] * v[440] + v[335] * v[443]
	+ v[336] * v[446]) + v[304] * (v[334] * v[455] + v[335] * v[458] + v[336] * v[461]) + v[337] * v[470] + v[338] * v[473]
	+ v[339] * v[476];
v[305] = v[46] + v[145] * v[49] + v[146] * v[50] + v[147] * v[51] - v[89] - v[201] * v[92] - v[202] * v[93] - v[203] * v[94];
v[306] = v[47] + v[145] * v[52] + v[146] * v[53] + v[147] * v[54] - v[90] - v[201] * v[95] - v[202] * v[96] - v[203] * v[97];
v[307] = -(v[100] * v[203]) + v[48] + v[145] * v[55] + v[146] * v[56] + v[147] * v[57] - v[91] - v[201] * v[98]
- v[202] * v[99];
if (v[732] > 0.1e-7) { v01 = 1e0 / v[732]; v02 = (-(v01 / v[732])); v03 = (2e0*v01) / (v[732] * v[732]); }
else {
	v01 = (12500000e0 / 3e0)*(24e0 - (-0.1e-7 + v[732])*(0.24e10 - 2e0*(-1e0 + 100000000e0*v[732])*(0.2399999997e10
		- 0.1199999994e18*v[732] - 0.3e17*(v[732] * v[732]))));
	v02 = (-50000000e0 / 3e0)*(0.3599999994e10 - 0.4799999982e18*v[732] + 0.6e25*Power(v[732], 3)
		+ 0.1799999982e26*(v[732] * v[732]));
	v03 = 0.1e17*(799999997e0 - 0.599999994e17*v[732] - 0.3e17*(v[732] * v[732]));
};
v[314] = v03;
v[315] = v02;
v[316] = v01;
v[317] = v[302] * v[316];
v[318] = v[303] * v[316];
v[319] = v[304] * v[316];
v[320] = sqrt((v[305] * v[305]) + (v[306] * v[306]) + (v[307] * v[307]));
if (v[320] > 0.1e-7) { v04 = 1e0 / v[320]; v05 = (-(v04 / v[320])); v06 = (2e0*v04) / (v[320] * v[320]); }
else {
	v04 = (12500000e0 / 3e0)*(24e0 - (-0.1e-7 + v[320])*(0.24e10 - 2e0*(-1e0 + 100000000e0*v[320])*(0.2399999997e10
		- 0.1199999994e18*v[320] - 0.3e17*(v[320] * v[320]))));
	v05 = (-50000000e0 / 3e0)*(0.3599999994e10 - 0.4799999982e18*v[320] + 0.6e25*Power(v[320], 3)
		+ 0.1799999982e26*(v[320] * v[320]));
	v06 = 0.1e17*(799999997e0 - 0.599999994e17*v[320] - 0.3e17*(v[320] * v[320]));
};
v[325] = v04;
v[326] = v[305] * v[325];
v[327] = v[306] * v[325];
v[328] = v[307] * v[325];
b329 = v[317] * v[326] + v[318] * v[327] + v[319] * v[328] < 0e0;
if (b329) {
	v[331] = -v[317];
	v[332] = -v[318];
	v[333] = -v[319];
}
else {
	v[331] = v[317];
	v[332] = v[318];
	v[333] = v[319];
};
v[4606] = 2e0*v[333];
v[3240] = (v[245] - v[248])*v[333];
v[3239] = (v[244] - v[247])*v[333];
v[2502] = v[3239] + v[333] * v[4459];
v[2486] = v[3240] + v[333] * v[4460];
v[1397] = -(v[333] * v[577]);
v[1395] = -(v[333] * v[576]);
v[1393] = -(v[333] * v[575]);
v[1391] = v[333] * v[478];
v[1389] = v[333] * v[477];
v[1387] = v[333] * v[476];
v[4474] = (v[246] - v[249])*v[333];
v[722] = (v[333] * v[333]);
v[4607] = 2e0*v[332];
v[3264] = (v[246] - v[249])*v[332];
v[3261] = (v[244] - v[247])*v[332];
v[2419] = v[3261] + v[332] * v[4459];
v[2359] = v[3264] + v[332] * v[4461];
v[1396] = -(v[332] * v[574]);
v[2332] = v[1396] + v[1397];
v[1394] = -(v[332] * v[573]);
v[2334] = v[1394] + v[1395];
v[1392] = -(v[332] * v[572]);
v[2336] = v[1392] + v[1393];
v[1390] = v[332] * v[475];
v[2338] = v[1390] + v[1391];
v[1388] = v[332] * v[474];
v[2340] = v[1388] + v[1389];
v[1386] = v[332] * v[473];
v[2342] = v[1386] + v[1387];
v[2759] = (v[245] - v[248])*v[332] + v[4474];
v[2315] = v[2342] * v[255] + v[2340] * v[265] + v[2338] * v[275] + v[2759] + v[2336] * v[281] + v[2334] * v[291]
+ v[2332] * v[301];
v[720] = (v[332] * v[332]);
v[4608] = 2e0*v[331];
v[4472] = v[331] * v[470];
v[4473] = v[1387] + v[4472];
v[4470] = v[331] * v[471];
v[4471] = v[1389] + v[4470];
v[4468] = v[331] * v[472];
v[4469] = v[1391] + v[4468];
v[4466] = -(v[331] * v[569]);
v[4467] = v[1393] + v[4466];
v[4464] = -(v[331] * v[570]);
v[4465] = v[1395] + v[4464];
v[4462] = -(v[331] * v[571]);
v[4463] = v[1397] + v[4462];
v[3286] = (v[246] - v[249])*v[331];
v[3285] = (v[245] - v[248])*v[331];
v[2298] = v[3285] + v[331] * v[4460];
v[2277] = v[3286] + v[331] * v[4461];
v[1352] = v[332] * v[4473] + v[473] * v[720];
v[1350] = v[332] * v[4471] + v[474] * v[720];
v[1348] = v[332] * v[4469] + v[475] * v[720];
v[1346] = v[332] * v[4467] - v[572] * v[720];
v[1344] = v[332] * v[4465] - v[573] * v[720];
v[1342] = v[332] * v[4463] - v[574] * v[720];
v[5072] = v[2332] + 2e0*v[4462];
v[2460] = v[1396] + v[4462];
v[5048] = 2e0*v[1397] + v[2460];
v[5060] = 2e0*v[1396] + v[4463];
v[5071] = v[2334] + 2e0*v[4464];
v[2462] = v[1394] + v[4464];
v[5047] = 2e0*v[1395] + v[2462];
v[5059] = 2e0*v[1394] + v[4465];
v[5070] = v[2336] + 2e0*v[4466];
v[2464] = v[1392] + v[4466];
v[5046] = 2e0*v[1393] + v[2464];
v[5058] = 2e0*v[1392] + v[4467];
v[5069] = v[2338] + 2e0*v[4468];
v[2466] = v[1390] + v[4468];
v[5045] = 2e0*v[1391] + v[2466];
v[5057] = 2e0*v[1390] + v[4469];
v[5068] = v[2340] + 2e0*v[4470];
v[2468] = v[1388] + v[4470];
v[5044] = 2e0*v[1389] + v[2468];
v[5056] = 2e0*v[1388] + v[4471];
v[5067] = v[2342] + 2e0*v[4472];
v[2470] = v[1386] + v[4472];
v[5043] = 2e0*v[1387] + v[2470];
v[5055] = 2e0*v[1386] + v[4473];
v[2734] = (v[244] - v[247])*v[331] + v[4474];
v[2709] = v[2734] + v[2759] - 2e0*v[4474];
v[2435] = v[2470] * v[255] + v[2468] * v[265] + v[2709] + v[2466] * v[275] + v[2464] * v[281] + v[2462] * v[291]
+ v[2460] * v[301];
v[2376] = v[2734] + v[301] * v[4463] + v[291] * v[4465] + v[281] * v[4467] + v[275] * v[4469] + v[265] * v[4471]
+ v[255] * v[4473];
v[1304] = v[2470] * v[333] + v[476] * v[722];
v[1302] = v[2468] * v[333] + v[477] * v[722];
v[1300] = v[2466] * v[333] + v[478] * v[722];
v[1298] = v[2464] * v[333] - v[575] * v[722];
v[1296] = v[2462] * v[333] - v[576] * v[722];
v[1294] = v[2460] * v[333] - v[577] * v[722];
v[718] = (v[331] * v[331]);
v[1408] = v[2342] * v[331] + v[470] * v[718];
v[1406] = v[2340] * v[331] + v[471] * v[718];
v[1404] = v[2338] * v[331] + v[472] * v[718];
v[1402] = v[2336] * v[331] - v[569] * v[718];
v[1400] = v[2334] * v[331] - v[570] * v[718];
v[1398] = v[2332] * v[331] - v[571] * v[718];
v[621] = -(v[13] * v[578]) - v[14] * v[584] - v[15] * v[590] - v[16] * v[596];
v[622] = -(v[13] * v[579]) - v[14] * v[585] - v[15] * v[591] - v[16] * v[597];
v[623] = -(v[13] * v[580]) - v[14] * v[586] - v[15] * v[592] - v[16] * v[598];
v[630] = -(v[13] * v[581]) - v[14] * v[587] - v[15] * v[593] - v[16] * v[599];
v[631] = -(v[13] * v[582]) - v[14] * v[588] - v[15] * v[594] - v[16] * v[600];
v[632] = -(v[13] * v[583]) - v[14] * v[589] - v[15] * v[595] - v[16] * v[601];
v[1240] = -(v[255] * v[621]) - v[265] * v[622] - v[275] * v[623] - v[3283] * v[624] - v[3262] * v[626] - v[3241] * v[628]
- v[281] * v[630] - v[291] * v[631] - v[301] * v[632];
v[6158] = v[624];
v[6159] = v[626];
v[6160] = v[628];
v[6161] = v[621];
v[6162] = v[622];
v[6163] = v[623];
v[6164] = -v[624];
v[6165] = -v[626];
v[6166] = -v[628];
v[6167] = v[630];
v[6168] = v[631];
v[6169] = v[632];
v[636] = -(v[17] * v[578]) - v[18] * v[584] - v[19] * v[590] - v[20] * v[596];
v[637] = -(v[17] * v[579]) - v[18] * v[585] - v[19] * v[591] - v[20] * v[597];
v[638] = -(v[17] * v[580]) - v[18] * v[586] - v[19] * v[592] - v[20] * v[598];
v[645] = -(v[17] * v[581]) - v[18] * v[587] - v[19] * v[593] - v[20] * v[599];
v[646] = -(v[17] * v[582]) - v[18] * v[588] - v[19] * v[594] - v[20] * v[600];
v[647] = -(v[17] * v[583]) - v[18] * v[589] - v[19] * v[595] - v[20] * v[601];
v[1242] = -(v[255] * v[636]) - v[265] * v[637] - v[275] * v[638] - v[3283] * v[639] - v[3262] * v[641] - v[3241] * v[643]
- v[281] * v[645] - v[291] * v[646] - v[301] * v[647];
v[6170] = v[639];
v[6171] = v[641];
v[6172] = v[643];
v[6173] = v[636];
v[6174] = v[637];
v[6175] = v[638];
v[6176] = -v[639];
v[6177] = -v[641];
v[6178] = -v[643];
v[6179] = v[645];
v[6180] = v[646];
v[6181] = v[647];
v[651] = -(v[21] * v[578]) - v[22] * v[584] - v[23] * v[590] - v[24] * v[596];
v[652] = -(v[21] * v[579]) - v[22] * v[585] - v[23] * v[591] - v[24] * v[597];
v[653] = -(v[21] * v[580]) - v[22] * v[586] - v[23] * v[592] - v[24] * v[598];
v[660] = -(v[21] * v[581]) - v[22] * v[587] - v[23] * v[593] - v[24] * v[599];
v[661] = -(v[21] * v[582]) - v[22] * v[588] - v[23] * v[594] - v[24] * v[600];
v[662] = -(v[21] * v[583]) - v[22] * v[589] - v[23] * v[595] - v[24] * v[601];
v[1236] = v[255] * v[651] + v[265] * v[652] + v[275] * v[653] + v[3283] * v[654] + v[3262] * v[656] + v[3241] * v[658]
+ v[281] * v[660] + v[291] * v[661] + v[301] * v[662];
v[6182] = v[654];
v[6183] = v[656];
v[6184] = v[658];
v[6185] = v[651];
v[6186] = v[652];
v[6187] = v[653];
v[6188] = -v[654];
v[6189] = -v[656];
v[6190] = -v[658];
v[6191] = v[660];
v[6192] = v[661];
v[6193] = v[662];
v[1147] = v[339] * v[624] + v[349] * v[639] - v[356] * v[654] - v[366] * v[669];
v[1146] = v[338] * v[624] + v[348] * v[639] - v[355] * v[654] - v[365] * v[669];
v[1145] = v[337] * v[624] + v[347] * v[639] - v[354] * v[654] - v[364] * v[669];
v[1151] = v[339] * v[626] + v[349] * v[641] - v[356] * v[656] - v[366] * v[671];
v[1150] = v[338] * v[626] + v[348] * v[641] - v[355] * v[656] - v[365] * v[671];
v[1149] = v[337] * v[626] + v[347] * v[641] - v[354] * v[656] - v[364] * v[671];
v[1155] = v[339] * v[628] + v[349] * v[643] - v[356] * v[658] - v[366] * v[673];
v[1154] = v[338] * v[628] + v[348] * v[643] - v[355] * v[658] - v[365] * v[673];
v[1153] = v[337] * v[628] + v[347] * v[643] - v[354] * v[658] - v[364] * v[673];
v[666] = -(v[25] * v[578]) - v[26] * v[584] - v[27] * v[590] - v[28] * v[596];
v[1159] = v[339] * v[621] + v[349] * v[636] - v[356] * v[651] - v[366] * v[666];
v[1158] = v[338] * v[621] + v[348] * v[636] - v[355] * v[651] - v[365] * v[666];
v[1157] = v[337] * v[621] + v[347] * v[636] - v[354] * v[651] - v[364] * v[666];
v[667] = -(v[25] * v[579]) - v[26] * v[585] - v[27] * v[591] - v[28] * v[597];
v[1163] = v[339] * v[622] + v[349] * v[637] - v[356] * v[652] - v[366] * v[667];
v[1162] = v[338] * v[622] + v[348] * v[637] - v[355] * v[652] - v[365] * v[667];
v[1161] = v[337] * v[622] + v[347] * v[637] - v[354] * v[652] - v[364] * v[667];
v[668] = -(v[25] * v[580]) - v[26] * v[586] - v[27] * v[592] - v[28] * v[598];
v[1167] = v[339] * v[623] + v[349] * v[638] - v[356] * v[653] - v[366] * v[668];
v[1166] = v[338] * v[623] + v[348] * v[638] - v[355] * v[653] - v[365] * v[668];
v[1165] = v[337] * v[623] + v[347] * v[638] - v[354] * v[653] - v[364] * v[668];
v[675] = -(v[25] * v[581]) - v[26] * v[587] - v[27] * v[593] - v[28] * v[599];
v[1183] = v[339] * v[630] + v[349] * v[645] - v[356] * v[660] - v[366] * v[675];
v[1182] = v[338] * v[630] + v[348] * v[645] - v[355] * v[660] - v[365] * v[675];
v[1181] = v[337] * v[630] + v[347] * v[645] - v[354] * v[660] - v[364] * v[675];
v[676] = -(v[25] * v[582]) - v[26] * v[588] - v[27] * v[594] - v[28] * v[600];
v[1187] = v[339] * v[631] + v[349] * v[646] - v[356] * v[661] - v[366] * v[676];
v[1186] = v[338] * v[631] + v[348] * v[646] - v[355] * v[661] - v[365] * v[676];
v[1185] = v[337] * v[631] + v[347] * v[646] - v[354] * v[661] - v[364] * v[676];
v[677] = -(v[25] * v[583]) - v[26] * v[589] - v[27] * v[595] - v[28] * v[601];
v[1238] = v[255] * v[666] + v[265] * v[667] + v[275] * v[668] + v[3283] * v[669] + v[3262] * v[671] + v[3241] * v[673]
+ v[281] * v[675] + v[291] * v[676] + v[301] * v[677];
v[1191] = v[339] * v[632] + v[349] * v[647] - v[356] * v[662] - v[366] * v[677];
v[1190] = v[338] * v[632] + v[348] * v[647] - v[355] * v[662] - v[365] * v[677];
v[1189] = v[337] * v[632] + v[347] * v[647] - v[354] * v[662] - v[364] * v[677];
v[6194] = v[669];
v[6195] = v[671];
v[6196] = v[673];
v[6197] = v[666];
v[6198] = v[667];
v[6199] = v[668];
v[6200] = -v[669];
v[6201] = -v[671];
v[6202] = -v[673];
v[6203] = v[675];
v[6204] = v[676];
v[6205] = v[677];
b678 = sqrt(Power(-(v[327] * v[331]) + v[326] * v[332], 2) + Power(v[328] * v[331] - v[326] * v[333], 2) + Power(-
(v[328] * v[332]) + v[327] * v[333], 2)) > 0.1e-7;
if (b678) {
	v[680] = -(v[328] * v[332]) + v[327] * v[333];
	v[681] = v[328] * v[331] - v[326] * v[333];
	v[682] = -(v[327] * v[331]) + v[326] * v[332];
	v[683] = sqrt((v[680] * v[680]) + (v[681] * v[681]) + (v[682] * v[682]));
	v[2159] = 1e0 / (v[683] * v[683]);
	v[1527] = v[683];
	v[2170] = 1e0 - (v[1527] * v[1527]);
	v[4987] = 1e0 / Power(v[2170], 0.15e1);
	v[2165] = 1e0 / sqrt(v[2170]);
	v[1526] = asin(v[1527]) / 2e0;
	v[2164] = 1e0 / Power(cos(v[1526]), 2);
	v[4609] = v[2164] * v[2165];
	v[685] = 2e0*tan(v[1526]);
	if (v[683] > 0.1e-7) { v07 = 1e0 / v[683]; v08 = (-(v07 / v[683])); v09 = (2e0*v07) / (v[683] * v[683]); }
	else {
		v07 = (12500000e0 / 3e0)*(24e0 - (-0.1e-7 + v[683])*(0.24e10 - 2e0*(-1e0 + 100000000e0*v[683])*
			(0.2399999997e10 - 0.1199999994e18*v[683] - 0.3e17*(v[683] * v[683]))));
		v08 = (-50000000e0 / 3e0)*(0.3599999994e10 - 0.4799999982e18*v[683] + 0.6e25*Power(v[683], 3)
			+ 0.1799999982e26*(v[683] * v[683]));
		v09 = 0.1e17*(799999997e0 - 0.599999994e17*v[683] - 0.3e17*(v[683] * v[683]));
	};
	v[689] = v09;
	v[690] = v08;
	v[691] = v07;
	v[4986] = v[685] * v[690] + v[4609] * v[691];
	v[4475] = v[685] * v[691];
	v[692] = v[4475] * v[680];
	v[4658] = 2e0*v[692];
	v[4515] = v[692] / 2e0;
	v[703] = (v[692] * v[692]);
	v[693] = v[4475] * v[681];
	v[4476] = v[693] / 2e0;
	v[701] = v[4476] * v[692];
	v[696] = (v[693] * v[693]);
	v[1481] = -v[696] - v[703];
	v[694] = v[4475] * v[682];
	v[4655] = 2e0*v[694];
	v[1511] = -v[694] + v[701];
	v[1502] = v[694] + v[701];
	v[708] = v[4476] * v[694];
	v[1493] = -v[692] + v[708];
	v[1485] = v[692] + v[708];
	v[706] = v[4515] * v[694];
	v[1506] = v[693] + v[706];
	v[1489] = -v[693] + v[706];
	v[697] = (v[694] * v[694]);
	v[1521] = 4e0 + v[696] + v[697] + v[703];
	v[4988] = 1e0 / Power(v[1521], 3);
	v[4654] = -4e0 / (v[1521] * v[1521]);
	v[1516] = -v[696] - v[697];
	v[1498] = -v[697] - v[703];
	v[695] = 4e0 / v[1521];
	v[4477] = v[695] / 2e0;
	v[698] = 1e0 + v[1516] * v[4477];
	v[699] = v[1511] * v[695];
	v[700] = v[1506] * v[695];
	v[702] = v[1502] * v[695];
	v[704] = 1e0 + v[1498] * v[4477];
	v[705] = v[1493] * v[695];
	v[707] = v[1489] * v[695];
	v[709] = v[1485] * v[695];
	v[710] = 1e0 + v[1481] * v[4477];
}
else {
	v[698] = 1e0;
	v[699] = 0e0;
	v[700] = 0e0;
	v[702] = 0e0;
	v[704] = 1e0;
	v[705] = 0e0;
	v[707] = 0e0;
	v[709] = 0e0;
	v[710] = 1e0;
};
if (b30) {
	v[1451] = 1e0 - v[722];
	v[1449] = 1e0 - v[720];
	v[1447] = 1e0 - v[718];
	v[715] = v[145] * v[173] + v[146] * v[174] + v[147] * v[175] - v[201] * v[229] - v[202] * v[230] - v[203] * v[231] + v[4478]
		+ v[1] * v[707] + v[2] * v[709] + v[3] * v[710];
	v[4481] = v[333] * v[715];
	v[714] = v[145] * v[170] + v[146] * v[171] + v[147] * v[172] - v[201] * v[226] - v[202] * v[227] - v[203] * v[228] + v[4479]
		+ v[1] * v[702] + v[2] * v[704] + v[3] * v[705];
	v[4483] = v[332] * v[714];
	v[4512] = v[4481] + v[4483];
	v[713] = v[145] * v[167] + v[146] * v[168] + v[147] * v[169] - v[201] * v[223] - v[202] * v[224] - v[203] * v[225] + v[4480]
		+ v[1] * v[698] + v[2] * v[699] + v[3] * v[700];
	v[4482] = -(v[331] * v[713]);
	v[4514] = v[4482] - v[4483];
	v[4513] = -v[4481] + v[4482];
	v[712] = -(v[331] * v[4512]) + v[1447] * v[713];
	v[716] = v[332] * v[4513] + v[1449] * v[714];
	v[717] = v[333] * v[4514] + v[1451] * v[715];
}
else {
	v[712] = 0e0;
	v[716] = 0e0;
	v[717] = 0e0;
};
v[719] = v[1408] * v[255] + v[1406] * v[265] + v[1404] * v[275] + v[1402] * v[281] + v[1400] * v[291] + v[1398] * v[301]
+ v[2759] * v[331] + v[3283] * v[718];
v[721] = v[1352] * v[255] + v[1350] * v[265] + v[1348] * v[275] + v[1346] * v[281] + v[1344] * v[291] + v[1342] * v[301]
+ v[2734] * v[332] + v[3262] * v[720];
v[723] = v[1304] * v[255] + v[1302] * v[265] + v[1300] * v[275] + v[1298] * v[281] + v[1296] * v[291] + v[1294] * v[301]
+ v[2709] * v[333] + v[3241] * v[722];
(*vnrel) = sqrt((v[719] * v[719]) + (v[721] * v[721]) + (v[723] * v[723]));
v[729] = -((v[4485] * Power(v[730], v[745])) / (v[35] * Power(v[33], v[750])));
v[4493] = -(v[35] * v[729]);
b733 = v[732] < v[32];
if (b733) {
	b735 = v[732] > v[33];
	if (b735) {
		v[744] = v[32] - v[732];
		v[4484] = -(v[31] * Power(v[744], v[34]));
		v[737] = v[331] * v[4484];
		v[739] = v[332] * v[4484];
		v[740] = v[333] * v[4484];
	}
	else {
		v[741] = -(v[31] * Power(v[730], v[34])) + v[729] * (Power(v[33], v[35]) - Power(v[732], v[35]));
		v[737] = v[331] * v[741];
		v[739] = v[332] * v[741];
		v[740] = v[333] * v[741];
	};
}
else {
	v[737] = 0e0;
	v[739] = 0e0;
	v[740] = 0e0;
};
if (b733) {
	v[4486] = 2e0*v[7];
	if (b735) {
		v[1278] = v[4485] * v[8];
		v[1279] = sqrt(v[1278] * Power(v[744], v[745]));
		v[747] = v[1279] * v[4486];
		v[746] = v[719] * v[747];
		v[748] = v[721] * v[747];
		v[749] = v[723] * v[747];
	}
	else {
		v[1286] = v[4493] * v[8];
		v[1287] = sqrt(v[1286] * Power(v[732], v[750]));
		v[751] = v[1287] * v[4486];
		v[746] = v[719] * v[751];
		v[748] = v[721] * v[751];
		v[749] = v[723] * v[751];
	};
	b752 = v[732] < v[32] && v[331] * (v[737] + v[746]) + v[332] * (v[739] + v[748]) + v[333] * (v[740] + v[749]) < 0e0;
	if (b752) {
		v[754] = v[746];
		v[755] = v[748];
		v[756] = v[749];
	}
	else {
		v[754] = -v[737];
		v[755] = -v[739];
		v[756] = -v[740];
	};
}
else {
	v[754] = 0e0;
	v[755] = 0e0;
	v[756] = 0e0;
};
v[757] = v[737] + v[754];
v[758] = v[739] + v[755];
v[759] = v[740] + v[756];
v[2570] = (v[757] * v[757]) + (v[758] * v[758]) + (v[759] * v[759]);
v[760] = v[4] * v[712];
v[761] = v[4] * v[716];
v[762] = v[4] * v[717];
v[766] = v[760] - (v[1145] * (v[244] - v[247]) + v[1149] * (v[245] - v[248]) + v[1153] * (v[246] - v[249])
	+ v[1157] * v[255] + v[1161] * v[265] + v[1165] * v[275] + v[1181] * v[281] + v[1185] * v[291] + v[1189] * v[301])*v[9];
v[767] = v[761] - (v[1146] * (v[244] - v[247]) + v[1150] * (v[245] - v[248]) + v[1154] * (v[246] - v[249])
	+ v[1158] * v[255] + v[1162] * v[265] + v[1166] * v[275] + v[1182] * v[281] + v[1186] * v[291] + v[1190] * v[301])*v[9];
v[768] = v[762] - (v[1147] * (v[244] - v[247]) + v[1151] * (v[245] - v[248]) + v[1155] * (v[246] - v[249])
	+ v[1159] * v[255] + v[1163] * v[265] + v[1167] * v[275] + v[1183] * v[281] + v[1187] * v[291] + v[1191] * v[301])*v[9];
v[2566] = (v[766] * v[766]) + (v[767] * v[767]) + (v[768] * v[768]);
if (b733) {
	if (b29) {
		b771 = sqrt((v[766] * v[766]) + (v[767] * v[767]) + (v[768] * v[768])) <= v[5] * sqrt((v[757] * v[757]) +
			(v[758] * v[758]) + (v[759] * v[759]));
		if (b771) {
			v[773] = v[766];
			v[774] = v[767];
			v[775] = v[768];
			v[776] = 1e0;
		}
		else {
			v[4487] = v[6] * sqrt(v[2570]);
			v[777] = sqrt(v[2566]);
			if (v[777] > 0.1e-5) { v010 = 1e0 / v[777]; v011 = (-(v010 / v[777])); v012 = (2e0*v010) / (v[777] * v[777]); }
			else {
				v010 = (24000000e0 - (-1e0 + 1000000e0*v[777])*(71999994e0 - 0.71999982e14*v[777] + 0.6e19*Power(v[777], 3)
					+ 0.23999982e20*(v[777] * v[777]))) / 24e0;
				v011 = (-500000e0 / 3e0)*(35999994e0 - 0.47999982e14*v[777] + 0.6e19*Power(v[777], 3) + 0.17999982e20*
					(v[777] * v[777]));
				v012 = 0.1e13*(7999997e0 - 0.5999994e13*v[777] - 0.3e13*(v[777] * v[777]));
			};
			v[781] = v011;
			v[782] = v010;
			v[783] = v[766] * v[782];
			v[784] = v[767] * v[782];
			v[785] = v[768] * v[782];
			v[773] = v[4487] * v[783];
			v[774] = v[4487] * v[784];
			v[775] = v[4487] * v[785];
			v[776] = 0e0;
		};
		if (sqrt((v[760] * v[760]) + (v[761] * v[761]) + (v[762] * v[762])) > v[5] * sqrt((v[757] * v[757]) + (v[758] * v[758]
			) + (v[759] * v[759]))) {
			if (v[4] > 0.1e-5) { v013 = 1e0 / v[4]; v014 = (-(v013 / v[4])); v015 = (2e0*v013) / (v[4] * v[4]); }
			else {
				v013 = (24000000e0 - (-1e0 + 1000000e0*v[4])*(71999994e0 - 0.71999982e14*v[4] + 0.6e19*Power(v[4], 3)
					+ 0.23999982e20*(v[4] * v[4]))) / 24e0;
				v014 = (-500000e0 / 3e0)*(35999994e0 - 0.47999982e14*v[4] + 0.6e19*Power(v[4], 3) + 0.17999982e20*
					(v[4] * v[4]));
				v015 = 0.1e13*(7999997e0 - 0.5999994e13*v[4] - 0.3e13*(v[4] * v[4]));
			};
			v[795] = sqrt((v[760] * v[760]) + (v[761] * v[761]) + (v[762] * v[762]));
			if (v[795] > 0.1e-5) { v016 = 1e0 / v[795]; v017 = (-(v016 / v[795])); v018 = (2e0*v016) / (v[795] * v[795]); }
			else {
				v016 = (24000000e0 - (-1e0 + 1000000e0*v[795])*(71999994e0 - 0.71999982e14*v[795] + 0.6e19*Power(v[795], 3)
					+ 0.23999982e20*(v[795] * v[795]))) / 24e0;
				v017 = (-500000e0 / 3e0)*(35999994e0 - 0.47999982e14*v[795] + 0.6e19*Power(v[795], 3) + 0.17999982e20*
					(v[795] * v[795]));
				v018 = 0.1e13*(7999997e0 - 0.5999994e13*v[795] - 0.3e13*(v[795] * v[795]));
			};
			v[802] = -(v013*v016*v[6] * sqrt(v[2570]));
			v[801] = v[712] + v[760] * v[802];
			v[803] = v[716] + v[761] * v[802];
			v[804] = v[717] + v[762] * v[802];
		}
		else {
			v[801] = 0e0;
			v[803] = 0e0;
			v[804] = 0e0;
		};
	}
	else {
		b805 = sqrt((v[766] * v[766]) + (v[767] * v[767]) + (v[768] * v[768])) <= v[6] * sqrt((v[757] * v[757]) +
			(v[758] * v[758]) + (v[759] * v[759]));
		if (b805) {
			v[773] = v[766];
			v[774] = v[767];
			v[775] = v[768];
			v[776] = 1e0;
		}
		else {
			v[816] = sqrt(v[2570]);
			v[4488] = v[6] * v[816];
			v[807] = sqrt(v[2566]);
			if (v[807] > 0.1e-5) { v019 = 1e0 / v[807]; v020 = (-(v019 / v[807])); v021 = (2e0*v019) / (v[807] * v[807]); }
			else {
				v019 = (24000000e0 - (-1e0 + 1000000e0*v[807])*(71999994e0 - 0.71999982e14*v[807] + 0.6e19*Power(v[807], 3)
					+ 0.23999982e20*(v[807] * v[807]))) / 24e0;
				v020 = (-500000e0 / 3e0)*(35999994e0 - 0.47999982e14*v[807] + 0.6e19*Power(v[807], 3) + 0.17999982e20*
					(v[807] * v[807]));
				v021 = 0.1e13*(7999997e0 - 0.5999994e13*v[807] - 0.3e13*(v[807] * v[807]));
			};
			v[811] = v020;
			v[812] = v019;
			v[813] = v[766] * v[812];
			v[814] = v[767] * v[812];
			v[815] = v[768] * v[812];
			v[773] = v[4488] * v[813];
			v[774] = v[4488] * v[814];
			v[775] = v[4488] * v[815];
			v[776] = 0e0;
		};
		if (sqrt((v[760] * v[760]) + (v[761] * v[761]) + (v[762] * v[762])) > v[6] * sqrt((v[757] * v[757]) + (v[758] * v[758]
			) + (v[759] * v[759]))) {
			if (v[4] > 0.1e-5) { v022 = 1e0 / v[4]; v023 = (-(v022 / v[4])); v024 = (2e0*v022) / (v[4] * v[4]); }
			else {
				v022 = (24000000e0 - (-1e0 + 1000000e0*v[4])*(71999994e0 - 0.71999982e14*v[4] + 0.6e19*Power(v[4], 3)
					+ 0.23999982e20*(v[4] * v[4]))) / 24e0;
				v023 = (-500000e0 / 3e0)*(35999994e0 - 0.47999982e14*v[4] + 0.6e19*Power(v[4], 3) + 0.17999982e20*
					(v[4] * v[4]));
				v024 = 0.1e13*(7999997e0 - 0.5999994e13*v[4] - 0.3e13*(v[4] * v[4]));
			};
			v[825] = sqrt((v[760] * v[760]) + (v[761] * v[761]) + (v[762] * v[762]));
			if (v[825] > 0.1e-5) { v025 = 1e0 / v[825]; v026 = (-(v025 / v[825])); v027 = (2e0*v025) / (v[825] * v[825]); }
			else {
				v025 = (24000000e0 - (-1e0 + 1000000e0*v[825])*(71999994e0 - 0.71999982e14*v[825] + 0.6e19*Power(v[825], 3)
					+ 0.23999982e20*(v[825] * v[825]))) / 24e0;
				v026 = (-500000e0 / 3e0)*(35999994e0 - 0.47999982e14*v[825] + 0.6e19*Power(v[825], 3) + 0.17999982e20*
					(v[825] * v[825]));
				v027 = 0.1e13*(7999997e0 - 0.5999994e13*v[825] - 0.3e13*(v[825] * v[825]));
			};
			v[831] = -(v022*v025*v[6] * sqrt(v[2570]));
			v[801] = v[712] + v[760] * v[831];
			v[803] = v[716] + v[761] * v[831];
			v[804] = v[717] + v[762] * v[831];
		}
		else {
			v[801] = 0e0;
			v[803] = 0e0;
			v[804] = 0e0;
		};
	};
}
else {
	v[773] = 0e0;
	v[774] = 0e0;
	v[775] = 0e0;
};
v[4583] = -(v[1155] * v[775]);
v[4579] = -(v[1151] * v[775]);
v[4575] = -(v[1147] * v[775]);
v[4582] = -(v[1154] * v[774]);
v[4578] = -(v[1150] * v[774]);
v[4574] = -(v[1146] * v[774]);
v[4581] = -(v[1153] * v[773]);
v[4577] = -(v[1149] * v[773]);
v[4573] = -(v[1145] * v[773]);
fn[0] = v[757];
fn[1] = v[758];
fn[2] = v[759];
ft[0] = v[773];
ft[1] = v[774];
ft[2] = v[775];
(*stickupdated) = v[776];
gtpupdated[0] = v[712] - v[801];
gtpupdated[1] = v[716] - v[803];
gtpupdated[2] = v[717] - v[804];
b845 = b733;
if (b845) {
	b846 = b735;
}
else {
};
b850 = b329;
if (b850) {
	v[851] = 0e0;
	v[852] = 0e0;
	v[853] = 0e0;
}
else {
	v[853] = 0e0;
	v[852] = 0e0;
	v[851] = 0e0;
};
v[858] = v[304] * v[851] + v[303] * v[852] + v[302] * v[853];
v[860] = v[315] * v[858];
v[4489] = v[860] / v[732];
v[861] = v[304] * v[4489] + v[740] + v[316] * v[851];
v[863] = v[303] * v[4489] + v[739] + v[316] * v[852];
v[864] = v[302] * v[4489] + v[737] + v[316] * v[853];
v[865] = v[861] * v[871];
v[4499] = v[865] / 2e0;
v[866] = v[861] * v[873];
v[920] = v[204] * v[866];
v[867] = v[861] * v[875];
v[926] = v[204] * v[867];
v[868] = v[861] * v[877];
v[4502] = v[868] / 2e0;
v[869] = v[861] * v[879];
v[909] = v[148] * v[869];
v[870] = v[861] * v[881];
v[915] = v[148] * v[870];
v[872] = v[863] * v[871];
v[921] = v[204] * v[872];
v[874] = v[863] * v[873];
v[4498] = v[874] / 2e0;
v[923] = v[4437] * v[874];
v[876] = v[863] * v[875];
v[929] = v[204] * v[876];
v[878] = v[863] * v[877];
v[910] = v[148] * v[878];
v[880] = v[863] * v[879];
v[4501] = v[880] / 2e0;
v[912] = v[4430] * v[880];
v[882] = v[863] * v[881];
v[918] = v[148] * v[882];
v[883] = v[864] * v[871];
v[927] = v[204] * v[883];
v[884] = v[864] * v[873];
v[930] = v[204] * v[884];
v[885] = v[864] * v[875];
v[4497] = v[885] / 2e0;
v[928] = v[4437] * v[885];
v[886] = -(v[231] * v[861]) - v[228] * v[863] - v[225] * v[864];
v[887] = -(v[230] * v[861]) - v[227] * v[863] - v[224] * v[864];
v[888] = -(v[229] * v[861]) - v[226] * v[863] - v[223] * v[864];
v[889] = v[864] * v[877];
v[916] = v[148] * v[889];
v[890] = v[864] * v[879];
v[919] = v[148] * v[890];
v[891] = v[864] * v[881];
v[4500] = v[891] / 2e0;
v[917] = v[4430] * v[891];
v[892] = v[175] * v[861] + v[172] * v[863] + v[169] * v[864];
v[893] = v[174] * v[861] + v[171] * v[863] + v[168] * v[864];
v[894] = v[173] * v[861] + v[170] * v[863] + v[167] * v[864];
v[895] = v[920] + v[921];
v[896] = v[926] + v[927];
v[897] = v[929] + v[930];
v[898] = v[4497] * v[482] + v[4498] * v[501] + v[4499] * v[519] + v[514] * v[866] + v[510] * v[867] + v[506] * v[872]
+ v[497] * v[876] + v[493] * v[883] + v[487] * v[884];
v[1138] = -(v[4428] * v[898]) + v[923] + v[928];
v[1135] = v[1138] + v[4437] * v[865] - v[923];
v[1132] = v[1135] + v[923] - v[928];
v[899] = v[82] * v[886] + v[81] * v[887] + v[80] * v[888];
v[900] = v[909] + v[910];
v[901] = v[915] + v[916];
v[902] = v[918] + v[919];
v[903] = v[383] * v[4500] + v[402] * v[4501] + v[420] * v[4502] + v[415] * v[869] + v[411] * v[870] + v[407] * v[878]
+ v[398] * v[882] + v[394] * v[889] + v[388] * v[890];
v[1131] = -(v[4426] * v[903]) + v[912] + v[917];
v[1128] = v[1131] + v[4430] * v[868] - v[912];
v[1125] = v[1128] + v[912] - v[917];
v[6146] = v[864];
v[6147] = v[863];
v[6148] = v[861];
v[6149] = v[1125] * v[368] + v[372] * v[901] + v[369] * v[902] + v[909] - v[910];
v[6150] = v[1128] * v[371] + v[372] * v[900] + v[370] * v[902] - v[915] + v[916];
v[6151] = v[1131] * v[373] + v[369] * v[900] + v[370] * v[901] + v[918] - v[919];
v[6152] = -v[864];
v[6153] = -v[863];
v[6154] = -v[861];
v[6155] = v[1132] * v[374] + v[378] * v[896] + v[375] * v[897] + v[920] - v[921];
v[6156] = v[1135] * v[377] + v[378] * v[895] + v[376] * v[897] - v[926] + v[927];
v[6157] = v[1138] * v[379] + v[375] * v[895] + v[376] * v[896] + v[929] - v[930];
v[904] = v[39] * v[892] + v[38] * v[893] + v[37] * v[894];
v[905] = (v[88] * v[886] + v[87] * v[887] + v[86] * v[888] - v[899]) / 2e0;
v[906] = (v[85] * v[886] + v[84] * v[887] + v[83] * v[888] - v[899]) / 2e0;
v[907] = (v[45] * v[892] + v[44] * v[893] + v[43] * v[894] - v[904]) / 2e0;
v[908] = (v[42] * v[892] + v[41] * v[893] + v[40] * v[894] - v[904]) / 2e0;
for (i840 = 1; i840 <= 12; i840++) {
	v[978] = (i840 == 10 ? 1 : 0);
	v[975] = (i840 == 11 ? 1 : 0);
	v[972] = (i840 == 12 ? 1 : 0);
	v[966] = (i840 == 4 ? 1 : 0);
	v[963] = (i840 == 5 ? 1 : 0);
	v[960] = (i840 == 6 ? 1 : 0);
	v[944] = v[6193 + i840];
	v[943] = v[6181 + i840];
	v[937] = v[6169 + i840];
	v[936] = v[6157 + i840];
	v[4490] = -v[936] - v[937];
	v[939] = v[6221 + i840];
	v[984] = -(v[4426] * v[939]);
	v[940] = v[6233 + i840];
	v[941] = v[6245 + i840];
	v[942] = v[6257 + i840];
	v[4491] = -v[943] - v[944];
	v[946] = v[6269 + i840];
	v[990] = -(v[4428] * v[946]);
	v[947] = v[6281 + i840];
	v[948] = v[6293 + i840];
	v[949] = v[6305 + i840];
	v[950] = v[6317 + i840];
	v[951] = v[6353 + i840];
	v[952] = v[6401 + i840];
	v[953] = v[6437 + i840];
	v[954] = v[6485 + i840];
	v[955] = v[6521 + i840];
	v[956] = (v[37] * v[4490] + v[40] * v[936] + v[43] * v[937]) / 2e0;
	v[957] = (v[38] * v[4490] + v[41] * v[936] + v[44] * v[937]) / 2e0;
	v[958] = (v[39] * v[4490] + v[42] * v[936] + v[45] * v[937]) / 2e0;
	v[959] = v[940] - v[960];
	v[961] = v[940] + v[960];
	v[962] = v[941] + v[963];
	v[964] = v[941] - v[963];
	v[965] = v[942] - v[966];
	v[967] = v[942] + v[966];
	v[968] = (v[4491] * v[80] + v[83] * v[943] + v[86] * v[944]) / 2e0;
	v[969] = (v[4491] * v[81] + v[84] * v[943] + v[87] * v[944]) / 2e0;
	v[970] = (v[4491] * v[82] + v[85] * v[943] + v[88] * v[944]) / 2e0;
	v[971] = v[947] - v[972];
	v[973] = v[947] + v[972];
	v[974] = v[948] + v[975];
	v[976] = v[948] - v[975];
	v[977] = v[949] - v[978];
	v[979] = v[949] + v[978];
	v[981] = -(v[383] * v[4742] * v[939]) + v[4430] * v[950];
	v[983] = v[148] * v[959] + v[388] * v[984];
	v[985] = v[148] * v[962] + v[394] * v[984];
	v[987] = -(v[4756] * v[482] * v[946]) + v[4437] * v[951];
	v[989] = v[204] * v[971] + v[487] * v[990];
	v[991] = v[204] * v[974] + v[493] * v[990];
	v[992] = v[6533 + i840] + v[167] * v[956] + v[168] * v[957] + v[169] * v[958] - v[223] * v[968] - v[224] * v[969]
		- v[225] * v[970] + v[881] * v[981] + v[879] * v[983] + v[877] * v[985] + v[875] * v[987] + v[873] * v[989]
		+ v[871] * v[991];
	v[993] = v[148] * v[961] + v[398] * v[984];
	v[994] = (-(v[148] * v[952]) + v[402] * v[984]) / 2e0;
	v[995] = v[148] * v[965] + v[407] * v[984];
	v[996] = v[204] * v[973] + v[497] * v[990];
	v[997] = (-(v[204] * v[953]) + v[501] * v[990]) / 2e0;
	v[998] = v[204] * v[977] + v[506] * v[990];
	v[999] = v[6545 + i840] + v[170] * v[956] + v[171] * v[957] + v[172] * v[958] - v[226] * v[968] - v[227] * v[969]
		- v[228] * v[970] + v[881] * v[993] + v[879] * v[994] + v[877] * v[995] + v[875] * v[996] + v[873] * v[997]
		+ v[871] * v[998];
	v[1000] = v[148] * v[964] + v[411] * v[984];
	v[1001] = v[1000] * v[861] + v[864] * v[981] + v[863] * v[993];
	v[1002] = v[148] * v[967] + v[415] * v[984];
	v[1003] = v[1002] * v[861] + v[864] * v[983] + v[863] * v[994];
	v[1004] = (-(v[148] * v[954]) + v[420] * v[984]) / 2e0;
	v[1005] = v[1004] * v[861] + v[864] * v[985] + v[863] * v[995];
	v[1006] = v[204] * v[976] + v[510] * v[990];
	v[1007] = v[1006] * v[861] + v[864] * v[987] + v[863] * v[996];
	v[1008] = v[204] * v[979] + v[514] * v[990];
	v[1009] = v[1008] * v[861] + v[864] * v[989] + v[863] * v[997];
	v[1010] = (-(v[204] * v[955]) + v[519] * v[990]) / 2e0;
	v[1011] = v[6557 + i840] + v[1010] * v[871] + v[1008] * v[873] + v[1006] * v[875] + v[1004] * v[877] + v[1002] * v[879]
		+ v[1000] * v[881] + v[173] * v[956] + v[174] * v[957] + v[175] * v[958] - v[229] * v[968] - v[230] * v[969]
		- v[231] * v[970];
	v[4492] = v[1011] * v[304] + v[302] * v[992] + v[303] * v[999];
	v[1012] = v[1010] * v[861] + v[864] * v[991] + v[863] * v[998];
	v[1013] = v[4492] / v[732];
	v[1015] = -(v[1014] * v[4492] * v[860]);
	v[1034] = v[1015];
	v[1016] = v[1013] * v[315];
	v[1017] = v[992];
	v[1032] = v[1017];
	v[1018] = v[999];
	v[1030] = v[1018];
	v[1019] = v[1011];
	v[1028] = v[1019];
	v[1020] = 0e0;
	v[1021] = 0e0;
	v[1022] = 0e0;
	b1023 = b733;
	if (b1023) {
		v[4496] = v[1017] * v[331];
		v[4495] = v[1018] * v[332];
		v[4494] = v[1019] * v[333];
		b1024 = b735;
		if (b1024) {
			v[1022] = v[1019] * v[4484];
			v[1019] = 0e0;
			v[1021] = v[1018] * v[4484];
			v[1018] = 0e0;
			v[1020] = v[1017] * v[4484];
			v[1017] = 0e0;
			v[1015] = v[1015] - v[31] * v[34] * (-v[4494] - v[4495] - v[4496])*Power(v[744], v[745]);
		}
		else {
			v[1022] = v[1028] * v[741];
			v[1019] = 0e0;
			v[1021] = v[1030] * v[741];
			v[1018] = 0e0;
			v[1020] = v[1032] * v[741];
			v[1017] = 0e0;
			v[1015] = v[1034] + v[4493] * (v[4494] + v[4495] + v[4496])*Power(v[732], v[750]);
		};
	}
	else {
	};
	b1035 = b329;
	if (b1035) {
		v[1036] = -v[1022];
		v[1037] = -v[1021];
		v[1038] = -v[1020];
	}
	else {
		v[1036] = v[1022];
		v[1037] = v[1021];
		v[1038] = v[1020];
	};
	v[1015] = v[1015] + v[1013] * v[314] * v[858] + v[315] * (v[1038] * v[302] + v[1037] * v[303] + v[1036] * v[304]
		+ v[1011] * v[851] + v[853] * v[992] + v[852] * v[999]);
	v[1046] = v[1036] * v[316] + v[1016] * v[851] + (v[1015] * v[304] + v[1011] * v[860]) / v[732];
	v[1048] = v[1037] * v[316] + v[1016] * v[852] + (v[1015] * v[303] + v[860] * v[999]) / v[732];
	v[1050] = v[1038] * v[316] + v[1016] * v[853] + (v[1015] * v[302] + v[860] * v[992]) / v[732];
	v[1051] = -(v[1046] * v[200]) - v[861] * v[970];
	v[1052] = -(v[1046] * v[199]) - v[861] * v[969];
	v[1053] = -(v[1046] * v[198]) - v[861] * v[968];
	v[1054] = v[1046] * v[144] + v[861] * v[958];
	v[1055] = v[1046] * v[143] + v[861] * v[957];
	v[1056] = v[1046] * v[142] + v[861] * v[956];
	v[1057] = -(v[1048] * v[200]) - v[863] * v[970];
	v[1058] = -(v[1048] * v[199]) - v[863] * v[969];
	v[1059] = -(v[1048] * v[198]) - v[863] * v[968];
	v[1060] = v[1048] * v[144] + v[863] * v[958];
	v[1061] = v[1048] * v[143] + v[863] * v[957];
	v[1062] = v[1048] * v[142] + v[863] * v[956];
	v[1063] = -(v[1050] * v[200]) - v[864] * v[970];
	v[1064] = -(v[1050] * v[199]) - v[864] * v[969];
	v[1065] = -(v[1050] * v[198]) - v[864] * v[968];
	v[1066] = v[1050] * v[144] + v[864] * v[958];
	v[1067] = v[1050] * v[143] + v[864] * v[957];
	v[1068] = v[1050] * v[142] + v[864] * v[956];
	v[1069] = v[100] * v[1051] + v[1053] * v[98] + v[1052] * v[99];
	v[1070] = v[1053] * v[95] + v[1052] * v[96] + v[1051] * v[97];
	v[1071] = v[1053] * v[92] + v[1052] * v[93] + v[1051] * v[94];
	v[1072] = v[100] * v[1057] + v[1059] * v[98] + v[1058] * v[99];
	v[1073] = v[1059] * v[95] + v[1058] * v[96] + v[1057] * v[97];
	v[1074] = v[1059] * v[92] + v[1058] * v[93] + v[1057] * v[94];
	v[1075] = v[100] * v[1063] + v[1065] * v[98] + v[1064] * v[99];
	v[1076] = v[1065] * v[95] + v[1064] * v[96] + v[1063] * v[97];
	v[1077] = v[1065] * v[92] + v[1064] * v[93] + v[1063] * v[94];
	v[1078] = (v[1069] * v[204] + v[865] * v[990]) / 2e0;
	v[1079] = v[1070] * v[204] + v[866] * v[990];
	v[1080] = v[1071] * v[204] + v[867] * v[990];
	v[1081] = v[1072] * v[204] + v[872] * v[990];
	v[1083] = v[1074] * v[204] + v[876] * v[990];
	v[1084] = v[1075] * v[204] + v[883] * v[990];
	v[1085] = v[1076] * v[204] + v[884] * v[990];
	v[4504] = v[4427] * v[898] * v[946] - v[4428] * (v[1069] * v[4537] + v[1073] * v[4538] + v[1077] * v[4539]
		+ v[1076] * v[487] + v[1075] * v[493] + v[1074] * v[497] + v[1072] * v[506] + v[1071] * v[510] + v[1070] * v[514]
		- v[4497] * v[951] - v[4498] * v[953] - v[4499] * v[955] + v[884] * v[971] + v[876] * v[973] + v[883] * v[974]
		+ v[867] * v[976] + v[872] * v[977] + v[866] * v[979]);
	v[1134] = v[1073] * v[4437] + v[4504] - v[4498] * v[990];
	v[1087] = (v[1077] * v[204] + v[885] * v[990]) / 2e0;
	v[1088] = -(v[100] * v[1012]) - v[1050] * v[225] - v[1048] * v[228] - v[1046] * v[231] - v[1007] * v[94]
		- v[1009] * v[97];
	v[1089] = -(v[1050] * v[224]) - v[1048] * v[227] - v[1046] * v[230] - v[1007] * v[93] - v[1009] * v[96] - v[1012] * v[99];
	v[1090] = -(v[1050] * v[223]) - v[1048] * v[226] - v[1046] * v[229] - v[1007] * v[92] - v[1009] * v[95] - v[1012] * v[98];
	v[1091] = v[1090] * v[80] + v[1089] * v[81] + v[1088] * v[82];
	v[1092] = v[1056] * v[55] + v[1055] * v[56] + v[1054] * v[57];
	v[1093] = v[1056] * v[52] + v[1055] * v[53] + v[1054] * v[54];
	v[1094] = v[1056] * v[49] + v[1055] * v[50] + v[1054] * v[51];
	v[1095] = v[1062] * v[55] + v[1061] * v[56] + v[1060] * v[57];
	v[1096] = v[1062] * v[52] + v[1061] * v[53] + v[1060] * v[54];
	v[1097] = v[1062] * v[49] + v[1061] * v[50] + v[1060] * v[51];
	v[1098] = v[1068] * v[55] + v[1067] * v[56] + v[1066] * v[57];
	v[1099] = v[1068] * v[52] + v[1067] * v[53] + v[1066] * v[54];
	v[1100] = v[1068] * v[49] + v[1067] * v[50] + v[1066] * v[51];
	v[1101] = (v[1092] * v[148] + v[868] * v[984]) / 2e0;
	v[1102] = v[1093] * v[148] + v[869] * v[984];
	v[1103] = v[1094] * v[148] + v[870] * v[984];
	v[1104] = v[1095] * v[148] + v[878] * v[984];
	v[1106] = v[1097] * v[148] + v[882] * v[984];
	v[1107] = v[1098] * v[148] + v[889] * v[984];
	v[1108] = v[1099] * v[148] + v[890] * v[984];
	v[4503] = v[4425] * v[903] * v[939] - v[4426] * (v[1099] * v[388] + v[1098] * v[394] + v[1097] * v[398] + v[1095] * v[407]
		+ v[1094] * v[411] + v[1093] * v[415] + v[1092] * v[4558] + v[1096] * v[4559] + v[1100] * v[4560] - v[4500] * v[950]
		- v[4501] * v[952] - v[4502] * v[954] + v[890] * v[959] + v[882] * v[961] + v[889] * v[962] + v[870] * v[964]
		+ v[878] * v[965] + v[869] * v[967]);
	v[1127] = v[1096] * v[4430] + v[4503] - v[4501] * v[984];
	v[1110] = (v[1100] * v[148] + v[891] * v[984]) / 2e0;
	v[1111] = v[1050] * v[169] + v[1048] * v[172] + v[1046] * v[175] + v[1001] * v[51] + v[1003] * v[54] + v[1005] * v[57];
	v[1112] = v[1050] * v[168] + v[1048] * v[171] + v[1046] * v[174] + v[1001] * v[50] + v[1003] * v[53] + v[1005] * v[56];
	v[1113] = v[1050] * v[167] + v[1048] * v[170] + v[1046] * v[173] + v[1001] * v[49] + v[1003] * v[52] + v[1005] * v[55];
	v[1114] = v[1113] * v[37] + v[1112] * v[38] + v[1111] * v[39];
	v[1115] = (-v[1091] + v[1090] * v[86] + v[1089] * v[87] + v[1088] * v[88]) / 2e0;
	v[1116] = (-v[1091] + v[1090] * v[83] + v[1089] * v[84] + v[1088] * v[85]) / 2e0;
	v[1117] = v[1080] + v[1084];
	v[1118] = v[1079] + v[1081];
	v[1119] = v[1083] + v[1085];
	v[1120] = (-v[1114] + v[1113] * v[43] + v[1112] * v[44] + v[1111] * v[45]) / 2e0;
	v[1121] = (-v[1114] + v[1113] * v[40] + v[1112] * v[41] + v[1111] * v[42]) / 2e0;
	v[1122] = v[1103] + v[1107];
	v[1123] = v[1102] + v[1104];
	v[1124] = v[1106] + v[1108];
	v[6690] = 0e0;
	v[6691] = 0e0;
	v[6692] = 0e0;
	v[6693] = 0e0;
	v[6694] = v[902];
	v[6695] = v[901];
	v[6696] = 0e0;
	v[6697] = 0e0;
	v[6698] = 0e0;
	v[6699] = 0e0;
	v[6700] = 0e0;
	v[6701] = 0e0;
	v[6666] = 0e0;
	v[6667] = 0e0;
	v[6668] = 0e0;
	v[6669] = v[902];
	v[6670] = 0e0;
	v[6671] = v[900];
	v[6672] = 0e0;
	v[6673] = 0e0;
	v[6674] = 0e0;
	v[6675] = 0e0;
	v[6676] = 0e0;
	v[6677] = 0e0;
	v[6642] = 0e0;
	v[6643] = 0e0;
	v[6644] = 0e0;
	v[6645] = v[901];
	v[6646] = v[900];
	v[6647] = 0e0;
	v[6648] = 0e0;
	v[6649] = 0e0;
	v[6650] = 0e0;
	v[6651] = 0e0;
	v[6652] = 0e0;
	v[6653] = 0e0;
	v[6594] = 0e0;
	v[6595] = 0e0;
	v[6596] = 0e0;
	v[6597] = 0e0;
	v[6598] = 0e0;
	v[6599] = 0e0;
	v[6600] = 0e0;
	v[6601] = 0e0;
	v[6602] = 0e0;
	v[6603] = 0e0;
	v[6604] = v[897];
	v[6605] = v[896];
	v[6582] = 0e0;
	v[6583] = 0e0;
	v[6584] = 0e0;
	v[6585] = 0e0;
	v[6586] = 0e0;
	v[6587] = 0e0;
	v[6588] = 0e0;
	v[6589] = 0e0;
	v[6590] = 0e0;
	v[6591] = v[897];
	v[6592] = 0e0;
	v[6593] = v[895];
	v[6570] = 0e0;
	v[6571] = 0e0;
	v[6572] = 0e0;
	v[6573] = 0e0;
	v[6574] = 0e0;
	v[6575] = 0e0;
	v[6576] = 0e0;
	v[6577] = 0e0;
	v[6578] = 0e0;
	v[6579] = v[896];
	v[6580] = v[895];
	v[6581] = 0e0;
	v[6714] = v[1050];
	v[6715] = v[1048];
	v[6716] = v[1046];
	v[6717] = v[1102] - v[1104] + (-v[1101] + v[1127])*v[368] + v[1124] * v[369] + v[1122] * v[372] + v[6689 + i840] / 2e0
		+ 2e0*v[1125] * v[966];
	v[6718] = -v[1103] + v[1107] + v[1124] * v[370] + v[1123] * v[372] + v[371] * (-v[1101] - v[1110] + v[4503]) + v[6665
		+ i840] / 2e0 + 2e0*v[1128] * v[963];
	v[6719] = v[1106] - v[1108] + v[1123] * v[369] + v[1122] * v[370] + (-v[1110] + v[1127])*v[373] + v[6641 + i840] / 2e0
		+ 2e0*v[1131] * v[960];
	v[6720] = -v[1050];
	v[6721] = -v[1048];
	v[6722] = -v[1046];
	v[6723] = v[1079] - v[1081] + (-v[1078] + v[1134])*v[374] + v[1119] * v[375] + v[1117] * v[378] + v[6593 + i840] / 2e0
		+ 2e0*v[1132] * v[978];
	v[6724] = -v[1080] + v[1084] + v[1119] * v[376] + v[1118] * v[378] + v[377] * (-v[1078] - v[1087] + v[4504]) + v[6581
		+ i840] / 2e0 + 2e0*v[1135] * v[975];
	v[6725] = v[1083] - v[1085] + v[1118] * v[375] + v[1117] * v[376] + (-v[1087] + v[1134])*v[379] + v[6569 + i840] / 2e0
		+ 2e0*v[1138] * v[972];
	Rc[i840 - 1] += v[6145 + i840] + v[908] * v[936] + v[907] * v[937] + v[906] * v[943] + v[905] * v[944];
	for (i933 = 1; i933 <= 12; i933++) {
		Kc[i840 - 1][i933 - 1] += v[1121] * v[6157 + i933] + v[1120] * v[6169 + i933] + v[1116] * v[6181 + i933] + v[1115] * v[6193
			+ i933] + v[6713 + i933];
	};/* end for */
};/* end for */
v[1196] = 0e0;
v[1197] = 0e0;
v[1198] = 0e0;
b1199 = b733;
if (b1199) {
	b1200 = b29;
	if (b1200) {
		b1201 = b771;
		if (b1201) {
			v[1198] = 0e0;
			v[1197] = 0e0;
			v[1196] = 0e0;
		}
		else {
		};
	}
	else {
		b1202 = b805;
		if (b1202) {
			v[1198] = 0e0;
			v[1197] = 0e0;
			v[1196] = 0e0;
		}
		else {
		};
	};
}
else {
};
v[4507] = v[1196] * v[9];
v[4506] = v[1197] * v[9];
v[4505] = v[1198] * v[9];
v[1219] = v[1236] * v[4505];
v[1220] = v[1238] * v[4505];
v[1221] = v[1240] * v[4505];
v[1222] = v[1242] * v[4505];
v[1237] = v[1236] * v[4506];
v[1239] = v[1238] * v[4506];
v[1241] = v[1240] * v[4506];
v[1243] = v[1242] * v[4506];
v[1257] = v[1236] * v[4507];
v[1258] = v[1238] * v[4507];
v[1259] = v[1240] * v[4507];
v[1260] = v[1242] * v[4507];
v[1261] = v[1198] * v[4];
v[2239] = -(v[1261] * v[333]);
v[1262] = v[1197] * v[4];
v[2240] = -(v[1262] * v[332]);
v[2242] = v[2239] + v[2240];
v[1263] = v[1196] * v[4];
v[2235] = -(v[1263] * v[331]);
v[2243] = v[2235] + v[2239];
v[2241] = v[2235] + v[2240];
v[1264] = 0e0;
v[1265] = 0e0;
v[1266] = 0e0;
v[1267] = 0e0;
v[1268] = 0e0;
b1269 = b733;
if (b1269) {
	v[1270] = 0e0;
	v[1271] = 0e0;
	v[1272] = 0e0;
	b1273 = b752;
	if (b1273) {
		v[1272] = 0e0;
		v[1271] = 0e0;
		v[1270] = 0e0;
	}
	else {
	};
	v[1277] = v[1270] * v[719] + v[1271] * v[721] + v[1272] * v[723];
	v[4508] = v[1277] * v[7];
	b1274 = b735;
	if (b1274) {
		v[1266] = v[1272] * v[747];
		v[1265] = v[1271] * v[747];
		v[1264] = v[1270] * v[747];
		v[1268] = (v[1278] * v[4508] * v[745] * Power(v[744], v[2531])) / v[1279];
	}
	else {
		v[1266] = v[1272] * v[751];
		v[1265] = v[1271] * v[751];
		v[1264] = v[1270] * v[751];
		v[1267] = (v[1286] * v[4508] * v[750] * Power(v[732], v[2534])) / v[1287];
	};
}
else {
};
v[4509] = -(v[1264] * v[331]);
v[4510] = -(v[1265] * v[332]);
v[4511] = -(v[1266] * v[333]);
v[2283] = v[331] * (v[4510] + v[4511]) - v[1264] * v[718];
v[4576] = v[10] * (v[2283] + (v[1145] * v[1196] + v[1146] * v[1197] + v[1147] * v[1198])*v[9]);
v[2282] = v[332] * (v[4509] + v[4511]) - v[1265] * v[720];
v[4580] = v[10] * (v[2282] + (v[1149] * v[1196] + v[1150] * v[1197] + v[1151] * v[1198])*v[9]);
v[2281] = v[333] * (v[4509] + v[4510]) - v[1266] * v[722];
v[4584] = v[10] * (v[2281] + (v[1153] * v[1196] + v[1154] * v[1197] + v[1155] * v[1198])*v[9]);
v[1288] = 0e0;
v[1289] = 0e0;
v[1290] = 0e0;
b1291 = b733;
if (b1291) {
	b1292 = b735;
	if (b1292) {
		v[1290] = 0e0;
		v[1289] = 0e0;
		v[1288] = 0e0;
		v[1267] = v[1267] - v[1268];
	}
	else {
	};
}
else {
};
v[1288] = v[1288] + v[1266] * v[2502];
v[1289] = v[1289] + v[1266] * v[2486];
v[1293] = v[1266] * v[2279];
v[1290] = v[1290] + v[1266] * v[2435];
v[1288] = v[1288] + v[1265] * v[2419];
v[1333] = v[1265] * v[2299];
v[1289] = v[1289] + v[1265] * v[2376];
v[1290] = v[1290] + v[1265] * v[2359];
v[1381] = v[1264] * v[2344];
v[1288] = v[1288] + v[1264] * v[2315];
v[1289] = v[1289] + v[1264] * v[2298];
v[1290] = v[1290] + v[1264] * v[2277];
v[1399] = v[1266] * v[1294] + v[1265] * v[1342] + v[1264] * v[1398] - (v[1189] * v[1196] + v[1190] * v[1197]
	+ v[1191] * v[1198])*v[9];
v[1624] = v[1399] * v[1795];
v[1598] = v[1399] * v[3207];
v[1401] = v[1266] * v[1296] + v[1265] * v[1344] + v[1264] * v[1400] - (v[1185] * v[1196] + v[1186] * v[1197]
	+ v[1187] * v[1198])*v[9];
v[1611] = v[1401] * v[3205];
v[4520] = v[1611] + v[1399] * v[1794];
v[4519] = v[1598] + v[1401] * v[1792];
v[1403] = v[1266] * v[1298] + v[1265] * v[1346] + v[1264] * v[1402] - (v[1181] * v[1196] + v[1182] * v[1197]
	+ v[1183] * v[1198])*v[9];
v[1622] = v[1403] * v[3202];
v[4762] = v[1622] + v[1624];
v[4521] = v[1622] + v[1401] * v[1793];
v[4750] = v[1624] + v[4521];
v[1613] = v[1403] * v[4454];
v[4758] = v[1613] + v[4520];
v[4752] = v[1611] + v[1613];
v[1601] = v[1403] * v[4457];
v[4761] = v[1601] + v[4519];
v[4751] = v[1598] + v[1601];
v[1405] = v[1266] * v[1300] + v[1265] * v[1348] + v[1264] * v[1404] - (v[1165] * v[1196] + v[1166] * v[1197]
	+ v[1167] * v[1198])*v[9];
v[1708] = v[1405] * v[1787];
v[1682] = v[1405] * v[3200];
v[1407] = v[1266] * v[1302] + v[1265] * v[1350] + v[1264] * v[1406] - (v[1161] * v[1196] + v[1162] * v[1197]
	+ v[1163] * v[1198])*v[9];
v[1695] = v[1407] * v[3198];
v[4541] = v[1695] + v[1405] * v[1786];
v[4540] = v[1682] + v[1407] * v[1784];
v[1409] = v[1266] * v[1304] + v[1265] * v[1352] + v[1264] * v[1408] - (v[1157] * v[1196] + v[1158] * v[1197]
	+ v[1159] * v[1198])*v[9];
v[1706] = v[1409] * v[3195];
v[4748] = v[1706] + v[1708];
v[4542] = v[1706] + v[1407] * v[1785];
v[4736] = v[1708] + v[4542];
v[1697] = v[1409] * v[4446];
v[4744] = v[1697] + v[4541];
v[4738] = v[1695] + v[1697];
v[1685] = v[1409] * v[4449];
v[4747] = v[1685] + v[4540];
v[4737] = v[1682] + v[1685];
v[1419] = v[2283] * v[291];
v[1420] = v[2283] * v[301];
v[1421] = v[2283] * v[281];
v[1422] = v[2282] * v[281];
v[1423] = v[2282] * v[301];
v[1424] = v[2281] * v[301];
v[1425] = v[2282] * v[291];
v[1426] = v[2281] * v[281];
v[1427] = v[2281] * v[291];
v[1428] = -(v[2283] * v[265]);
v[1429] = -(v[2283] * v[275]);
v[1430] = -(v[2283] * v[255]);
v[1431] = -(v[2282] * v[255]);
v[1432] = -(v[2282] * v[275]);
v[1433] = -(v[2281] * v[275]);
v[1434] = -(v[2282] * v[265]);
v[1435] = -(v[2281] * v[255]);
v[1436] = -(v[2281] * v[265]);
v[1437] = 0e0;
v[1438] = 0e0;
v[1439] = 0e0;
v[1440] = 0e0;
v[1441] = 0e0;
v[1442] = 0e0;
v[1443] = 0e0;
v[1444] = 0e0;
v[1445] = 0e0;
b1446 = b30;
if (b1446) {
	v[1293] = v[1293] - v[1261] * v[715];
	v[1333] = v[1333] - v[1262] * v[714];
	v[1448] = v[1263] * v[1447] + v[2242] * v[331];
	v[1450] = v[1262] * v[1449] + v[2243] * v[332];
	v[1452] = v[1261] * v[1451] + v[2241] * v[333];
	v[1381] = v[1381] - v[1263] * v[713];
	v[1288] = v[1288] - v[1263] * v[4512] + v[2242] * v[713];
	v[1289] = v[1289] + v[1262] * v[4513] + v[2243] * v[714];
	v[1290] = v[1290] + v[1261] * v[4514] + v[2241] * v[715];
	v[1437] = v[1] * v[1448];
	v[1438] = v[1448] * v[2];
	v[1439] = v[1448] * v[3];
	v[1456] = -v[1448];
	v[1457] = -(v[1448] * v[203]);
	v[1458] = -(v[1448] * v[202]);
	v[1459] = -(v[1448] * v[201]);
	v[1460] = v[1448];
	v[1461] = v[1448] * v[147];
	v[1462] = v[1448] * v[146];
	v[1463] = v[1448] * v[145];
	v[1440] = v[1] * v[1450];
	v[1441] = v[1450] * v[2];
	v[1442] = v[1450] * v[3];
	v[1464] = -v[1450];
	v[1465] = -(v[1450] * v[203]);
	v[1466] = -(v[1450] * v[202]);
	v[1467] = -(v[1450] * v[201]);
	v[1468] = v[1450];
	v[1469] = v[1450] * v[147];
	v[1470] = v[1450] * v[146];
	v[1471] = v[145] * v[1450];
	v[1443] = v[1] * v[1452];
	v[1444] = v[1452] * v[2];
	v[1445] = v[1452] * v[3];
	v[1472] = -v[1452];
	v[1473] = -(v[1452] * v[203]);
	v[1474] = -(v[1452] * v[202]);
	v[1475] = -(v[1452] * v[201]);
	v[1476] = v[1452];
	v[1477] = v[1452] * v[147];
	v[1478] = v[1452] * v[146];
	v[1479] = v[145] * v[1452];
}
else {
	v[1463] = 0e0;
	v[1462] = 0e0;
	v[1461] = 0e0;
	v[1471] = 0e0;
	v[1470] = 0e0;
	v[1469] = 0e0;
	v[1479] = 0e0;
	v[1478] = 0e0;
	v[1477] = 0e0;
	v[1460] = 0e0;
	v[1468] = 0e0;
	v[1476] = 0e0;
	v[1459] = 0e0;
	v[1458] = 0e0;
	v[1457] = 0e0;
	v[1467] = 0e0;
	v[1466] = 0e0;
	v[1465] = 0e0;
	v[1475] = 0e0;
	v[1474] = 0e0;
	v[1473] = 0e0;
	v[1456] = 0e0;
	v[1464] = 0e0;
	v[1472] = 0e0;
};
v[4615] = v[1437] / 2e0;
v[4616] = v[1441] / 2e0;
v[4617] = v[1445] / 2e0;
b1480 = b678;
if (b1480) {
	v[1514] = -(v[1438] * v[695]);
	v[1509] = v[1439] * v[695];
	v[1496] = v[1442] * v[695];
	v[1483] = -(v[1445] * v[4477]);
	v[1487] = v[1444] * v[695];
	v[1491] = v[1443] * v[695];
	v[1495] = v[1487] + v[1496];
	v[1500] = -(v[1441] * v[4477]);
	v[1504] = v[1440] * v[695];
	v[1508] = v[1491] + v[1509];
	v[1515] = v[1504] - v[1514];
	v[1517] = v[1444] * v[1485] + v[1443] * v[1489] + v[1442] * v[1493] + v[1440] * v[1502] + v[1439] * v[1506]
		+ v[1438] * v[1511] + v[1516] * v[4615] + v[1498] * v[4616] + v[1481] * v[4617];
	v[2183] = v[1483] + v[1500] - (4e0*v[1517]) / (v[1521] * v[1521]);
	v[4614] = 4e0*v[2183];
	v[2181] = -v[1483] + v[2183] - v[1437] * v[4477];
	v[4613] = 4e0*(v[1483] - v[1500] + v[2181]);
	v[1522] = v[1504] + v[1514] + v[1495] * v[4476] + v[1508] * v[4515] + 2e0*v[2181] * v[694];
	v[1524] = (-2e0*v[1491] + 2e0*v[1509] + v[1515] * v[692] + v[4613] * v[693] + v[1495] * v[694]) / 2e0;
	v[1525] = (2e0*v[1487] - 2e0*v[1496] + v[4614] * v[692] + v[1515] * v[693] + v[1508] * v[694]) / 2e0;
	v[4516] = v[1525] * v[680] + v[1524] * v[681] + v[1522] * v[682];
	v[2169] = v[4516] * v[691];
	v[2166] = v[4516] * v[685];
	v[1528] = v[2166] * v[690] + v[2169] / (Power(cos(v[1526]), 2)*sqrt(v[2170]));
	v[4517] = v[1528] / v[683];
	v[1529] = v[1522] * v[4475] + v[4517] * v[682];
	v[1531] = v[1524] * v[4475] + v[4517] * v[681];
	v[1532] = v[1525] * v[4475] + v[4517] * v[680];
	v[1288] = v[1288] - v[1529] * v[327] + v[1531] * v[328];
	v[1289] = v[1289] + v[1529] * v[326] - v[1532] * v[328];
	v[1290] = v[1290] - v[1531] * v[326] + v[1532] * v[327];
}
else {
};
v[1288] = v[1288] + v[1381] * v[4608];
v[1289] = v[1289] + v[1333] * v[4607];
v[1290] = v[1290] + v[1293] * v[4606];
b1533 = b329;
if (b1533) {
	v[1534] = -v[1290];
	v[1535] = -v[1289];
	v[1536] = -v[1288];
}
else {
	v[1534] = v[1290];
	v[1535] = v[1289];
	v[1536] = v[1288];
};
v[1541] = v[1536] * v[302] + v[1535] * v[303] + v[1534] * v[304];
v[1267] = v[1267] + v[1541] * v[315];
v[4518] = v[1267] / v[732];
v[1543] = v[1534] * v[316] + v[304] * v[4518];
v[1544] = v[1535] * v[316] + v[303] * v[4518];
v[1545] = v[1536] * v[316] + v[302] * v[4518];
v[1472] = v[1472] - v[1543];
v[1476] = v[1476] + v[1543];
v[1464] = v[1464] - v[1544];
v[1468] = v[1468] + v[1544];
v[1456] = v[1456] - v[1545];
v[1460] = v[1460] + v[1545];
v[1546] = v[1399] * v[2275];
v[1547] = v[1399] * v[2274];
v[1548] = v[1399] * v[2273];
v[1551] = v[1401] * v[2270];
v[1552] = v[1399] * v[288] + v[1401] * v[290];
v[4525] = v[1552] * v[209];
v[1555] = v[1401] * v[2268];
v[1556] = v[1546] + v[1551];
v[1558] = v[1401] * v[2269] + v[4525] / v[279];
v[1561] = v[1403] * v[2263];
v[1562] = v[1399] * v[282] + v[1403] * v[290];
v[4526] = v[1562] * v[214];
v[1563] = v[1401] * v[282] + v[1403] * v[288];
v[4527] = v[1563] * v[218];
v[5034] = -(v[1399] * v[4374]) - v[1401] * v[4375] - v[1403] * v[4376] + v[4453] * v[4525] + v[278] * v[4526]
+ v[277] * v[4527];
v[1565] = v[1403] * v[2265] + v[4526] / v[279];
v[1566] = v[1558] + v[1565];
v[1568] = v[1403] * v[2264] + v[4527] / v[279];
v[1569] = v[1547] + v[1568];
v[1570] = v[1405] * v[2260];
v[1571] = v[1405] * v[2259];
v[1572] = v[1405] * v[2258];
v[1575] = v[1407] * v[2255];
v[1576] = v[1405] * v[262] + v[1407] * v[264];
v[4546] = v[153] * v[1576];
v[1579] = v[1407] * v[2253];
v[1580] = v[1570] + v[1575];
v[1582] = v[1407] * v[2254] + v[4546] / v[253];
v[1585] = v[1409] * v[2248];
v[1586] = v[1405] * v[256] + v[1409] * v[264];
v[4547] = v[158] * v[1586];
v[1587] = v[1407] * v[256] + v[1409] * v[262];
v[4548] = v[1587] * v[162];
v[5038] = -(v[1405] * v[4399]) - v[1407] * v[4400] - v[1409] * v[4401] + v[4445] * v[4546] + v[252] * v[4547]
+ v[251] * v[4548];
v[1589] = v[1409] * v[2250] + v[4547] / v[253];
v[1590] = v[1582] + v[1589];
v[1592] = v[1409] * v[2249] + v[4548] / v[253];
v[1593] = v[1571] + v[1592];
v[1473] = v[1473] - v[1543] * v[200] + v[1219] * v[353] + v[1220] * v[363];
v[1474] = v[1474] - v[1543] * v[199] + v[1219] * v[352] + v[1220] * v[361];
v[1475] = v[1475] - v[1543] * v[198] + v[1219] * v[351] + v[1220] * v[359];
v[1606] = v[100] * v[1473] + v[290] * v[4761] + v[1475] * v[98] + v[1474] * v[99];
v[1607] = v[1563] * v[4457] + v[288] * v[4519] + v[1475] * v[95] + v[1474] * v[96] + v[1473] * v[97];
v[1608] = v[282] * v[4751] + v[1475] * v[92] + v[1474] * v[93] + v[1473] * v[94];
v[1465] = v[1465] - v[1544] * v[200] + v[1237] * v[353] + v[1239] * v[363];
v[1466] = v[1466] - v[1544] * v[199] + v[1237] * v[352] + v[1239] * v[361];
v[1467] = v[1467] - v[1544] * v[198] + v[1237] * v[351] + v[1239] * v[359];
v[1618] = v[100] * v[1465] + v[1562] * v[4454] + v[290] * v[4520] + v[1467] * v[98] + v[1466] * v[99];
v[1619] = v[288] * v[4758] + v[1467] * v[95] + v[1466] * v[96] + v[1465] * v[97];
v[1620] = v[282] * v[4752] + v[1467] * v[92] + v[1466] * v[93] + v[1465] * v[94];
v[1457] = v[1457] - v[1545] * v[200] + v[1257] * v[353] + v[1258] * v[363];
v[1458] = v[1458] - v[1545] * v[199] + v[1257] * v[352] + v[1258] * v[361];
v[1459] = v[1459] - v[1545] * v[198] + v[1257] * v[351] + v[1258] * v[359];
v[1630] = v[100] * v[1457] + v[1552] * v[1793] + v[290] * v[4762] + v[1459] * v[98] + v[1458] * v[99];
v[1631] = v[288] * v[4521] + v[1459] * v[95] + v[1458] * v[96] + v[1457] * v[97];
v[1632] = v[282] * v[4750] + v[1459] * v[92] + v[1458] * v[93] + v[1457] * v[94];
v[1633] = -(v[1419] * v[875]);
v[1634] = -(v[1419] * v[873]);
v[1635] = -(v[1419] * v[871]);
v[1636] = -(v[1420] * v[875]);
v[1637] = -(v[1420] * v[871]);
v[1638] = -(v[1420] * v[873]);
v[1639] = -(v[1421] * v[873]);
v[1640] = -(v[1421] * v[871]);
v[1641] = -(v[1421] * v[875]);
v[1642] = -(v[1422] * v[875]);
v[1643] = -(v[1422] * v[873]);
v[1644] = -(v[1422] * v[871]);
v[1645] = -(v[1423] * v[871]);
v[1646] = -(v[1423] * v[875]);
v[1647] = -(v[1423] * v[873]);
v[1648] = -(v[1424] * v[873]);
v[1649] = -(v[1424] * v[875]);
v[1650] = -(v[1424] * v[871]);
v[1651] = v[1639] + v[1642] + v[1645] + v[1648] + v[1555] * v[4596] - v[1566] * v[489] - v[1556] * v[492];
v[1652] = -(v[1425] * v[875]);
v[1653] = -(v[1425] * v[871]);
v[1654] = -(v[1425] * v[873]);
v[1655] = -v[1634] - v[1637] - v[1649] - v[1652] + v[1561] * v[4595] - v[1566] * v[486] + v[1569] * v[492];
v[1656] = v[1631] * v[204] - v[1639] * v[479] + v[1634] * v[480] - v[1638] * v[481];
v[1657] = -(v[1426] * v[875]);
v[1658] = -(v[1426] * v[873]);
v[1659] = -(v[1426] * v[871]);
v[1660] = -(v[1427] * v[873]);
v[1661] = -(v[1427] * v[875]);
v[1662] = -(v[1427] * v[871]);
v[1666] = -v[1640] - v[1653] - v[1657] - v[1660] + v[1548] * v[4594] - v[1556] * v[486] + v[1569] * v[489];
v[1667] = v[1630] * v[204] - v[1640] * v[479] + v[1635] * v[480] - v[1637] * v[481];
v[1668] = v[1620] * v[204] - v[1642] * v[479] + v[1652] * v[480] - v[1646] * v[481];
v[1669] = v[1636] + v[1647];
v[1670] = v[1618] * v[204] - v[1644] * v[479] + v[1653] * v[480] - v[1645] * v[481];
v[1671] = v[1608] * v[204] - v[1657] * v[479] + v[1661] * v[480] - v[1649] * v[481];
v[1672] = v[1607] * v[204] - v[1658] * v[479] + v[1660] * v[480] - v[1648] * v[481];
v[1673] = v[1643] + v[1659];
v[1674] = v[1633] + v[1662];
v[7394] = 0e0;
v[7395] = 0e0;
v[7396] = 0e0;
v[7397] = 0e0;
v[7398] = 0e0;
v[7399] = 0e0;
v[7400] = 0e0;
v[7401] = 0e0;
v[7402] = 0e0;
v[7403] = -v[1655] / 2e0 - v[1673];
v[7404] = v[1651] / 2e0 - v[1674];
v[7405] = -v[1666] / 2e0 - v[1669];
v[1675] = 1e0 / (v[279] * v[279]);
v[4765] = -(v[1675] * v[289]);
v[4536] = -(v[1675] * v[288]);
v[4535] = -(v[1675] * v[290]);
v[4532] = -(v[1675] * v[282]);
v[4531] = -(v[1675] * v[277]);
v[4530] = -(v[1399] * v[1675]);
v[4529] = -(v[1675] * v[278]);
v[4528] = -(v[1675] * v[4453]);
v[3161] = -(v[1675] * (v[209] * v[280] + v[4451]));
v[3160] = -(v[1675] * (v[208] * v[280] + v[4452]));
v[3159] = -(v[1675] * v[4522]);
v[3158] = -(v[1675] * v[4455]);
v[3157] = -(v[1675] * v[4523]);
v[3156] = -(v[1675] * v[4524]);
v[3155] = -(v[1675] * v[280]);
v[4723] = v[282] * v[3155];
v[3154] = -(v[1675] * v[283]);
v[4722] = v[288] * v[3154];
v[4721] = -(v[1675] * v[290] * v[292]);
v[3152] = -(v[1675] * v[4525]);
v[3151] = -(v[1675] * v[4526]);
v[3150] = -(v[1675] * v[4527]);
v[3149] = v[1552] * v[4528];
v[3148] = v[1562] * v[4529];
v[3147] = v[1563] * v[4531];
v[2980] = v[1401] * v[4528];
v[2979] = v[4456] * v[4530];
v[2977] = v[1403] * v[3155];
v[4703] = v[2977] + v[2980];
v[5027] = v[2979] + v[4703];
v[2967] = v[1403] * v[4529];
v[2964] = v[1401] * v[3154];
v[4702] = v[2964] + v[2967];
v[2963] = v[4458] * v[4530];
v[5026] = v[2963] + v[4702];
v[2953] = v[1403] * v[4531];
v[2951] = v[1401] * v[4765];
v[2950] = v[292] * v[4530];
v[4700] = v[2950] + v[2953];
v[5025] = v[2951] + v[4700];
v[2913] = v[216] * v[4532];
v[2910] = v[211] * v[4532];
v[2907] = -(v[1675] * v[4533]);
v[2906] = -(v[1675] * v[4534]);
v[2901] = v[214] * v[4535];
v[2900] = v[213] * v[4536];
v[4144] = v[2900] + v[2910];
v[4136] = v[2901] + v[4144];
v[4122] = v[2900] + v[2901];
v[2898] = v[207] * v[4532];
v[4141] = v[2898] + v[2906] + v[2907];
v[4134] = -v[2907] + v[4141];
v[4124] = -v[2906] + v[4141];
v[2894] = v[219] * v[4535];
v[4147] = v[2894] + v[2913];
v[2893] = v[218] * v[4536];
v[4130] = v[2893] + v[2894];
v[4126] = v[2893] + v[4147];
v[2272] = v[288] * v[3156] + v[282] * v[3157] + v[4721];
v[2267] = v[290] * v[3158] + v[282] * v[3159] + v[4722];
v[2262] = v[288] * v[3160] + v[290] * v[3161] + v[4723];
v[2073] = v[209] * v[4528];
v[2071] = v[214] * v[4529];
v[2067] = v[218] * v[4531];
v[3146] = v[1548] + v[1555] + v[1561] + v[1563] * v[2067] + v[1562] * v[2071] + v[1552] * v[2073] + v[1403] * v[2262]
+ v[1401] * v[2267] + v[1399] * v[2272];
v[1676] = v[1635] - v[1638] - v[1644] + v[1646] + v[1658] - v[1661] + v[1651] * v[375] - v[1655] * v[376]
- v[1666] * v[378] + v[1606] * v[4537] + v[1619] * v[4538] + v[1632] * v[4539] + v[3146] * v[4586] + v[1669] * v[4718]
+ v[1673] * v[4719] + v[1674] * v[4720] + v[1631] * v[487] + v[1630] * v[493] + v[1620] * v[497] + v[1618] * v[506]
+ v[1608] * v[510] + v[1607] * v[514];
v[1477] = v[1477] + v[144] * v[1543] + v[1221] * v[336] + v[1222] * v[346];
v[1478] = v[1478] + v[143] * v[1543] + v[1221] * v[335] + v[1222] * v[344];
v[1479] = v[1479] + v[142] * v[1543] + v[1221] * v[334] + v[1222] * v[342];
v[1690] = v[264] * v[4747] + v[1479] * v[55] + v[1478] * v[56] + v[1477] * v[57];
v[1691] = v[1587] * v[4449] + v[262] * v[4540] + v[1479] * v[52] + v[1478] * v[53] + v[1477] * v[54];
v[1692] = v[256] * v[4737] + v[1479] * v[49] + v[1478] * v[50] + v[1477] * v[51];
v[1469] = v[1469] + v[144] * v[1544] + v[1241] * v[336] + v[1243] * v[346];
v[1470] = v[1470] + v[143] * v[1544] + v[1241] * v[335] + v[1243] * v[344];
v[1471] = v[1471] + v[142] * v[1544] + v[1241] * v[334] + v[1243] * v[342];
v[1702] = v[1586] * v[4446] + v[264] * v[4541] + v[1471] * v[55] + v[1470] * v[56] + v[1469] * v[57];
v[1703] = v[262] * v[4744] + v[1471] * v[52] + v[1470] * v[53] + v[1469] * v[54];
v[1704] = v[256] * v[4738] + v[1471] * v[49] + v[1470] * v[50] + v[1469] * v[51];
v[1461] = v[1461] + v[144] * v[1545] + v[1259] * v[336] + v[1260] * v[346];
v[1462] = v[1462] + v[143] * v[1545] + v[1259] * v[335] + v[1260] * v[344];
v[1463] = v[1463] + v[142] * v[1545] + v[1259] * v[334] + v[1260] * v[342];
v[1714] = v[1576] * v[1785] + v[264] * v[4748] + v[1463] * v[55] + v[1462] * v[56] + v[1461] * v[57];
v[1715] = v[262] * v[4542] + v[1463] * v[52] + v[1462] * v[53] + v[1461] * v[54];
v[1716] = v[256] * v[4736] + v[1463] * v[49] + v[1462] * v[50] + v[1461] * v[51];
v[1717] = v[1428] * v[881];
v[1718] = v[1428] * v[879];
v[1719] = v[1428] * v[877];
v[1720] = v[1429] * v[881];
v[1721] = v[1429] * v[877];
v[1722] = v[1429] * v[879];
v[1723] = v[1430] * v[879];
v[1724] = v[1430] * v[877];
v[1725] = v[1430] * v[881];
v[1726] = v[1431] * v[881];
v[1727] = v[1431] * v[879];
v[1728] = v[1431] * v[877];
v[1729] = v[1432] * v[877];
v[1730] = v[1432] * v[881];
v[1731] = v[1432] * v[879];
v[1732] = v[1433] * v[879];
v[1733] = v[1433] * v[881];
v[1734] = v[1433] * v[877];
v[1735] = v[1723] + v[1726] + v[1729] + v[1732] - v[1590] * v[390] - v[1580] * v[393] + v[1579] * v[4592];
v[1736] = v[1434] * v[881];
v[1737] = v[1434] * v[877];
v[1738] = v[1434] * v[879];
v[1739] = -v[1718] - v[1721] - v[1733] - v[1736] - v[1590] * v[387] + v[1593] * v[393] + v[1585] * v[4591];
v[1740] = v[148] * v[1715] - v[1723] * v[380] + v[1718] * v[381] - v[1722] * v[382];
v[1741] = v[1435] * v[881];
v[1742] = v[1435] * v[879];
v[1743] = v[1435] * v[877];
v[1744] = v[1436] * v[879];
v[1745] = v[1436] * v[881];
v[1746] = v[1436] * v[877];
v[1750] = -v[1724] - v[1737] - v[1741] - v[1744] - v[1580] * v[387] + v[1593] * v[390] + v[1572] * v[4590];
v[1751] = v[148] * v[1714] - v[1724] * v[380] + v[1719] * v[381] - v[1721] * v[382];
v[1752] = v[148] * v[1704] - v[1726] * v[380] + v[1736] * v[381] - v[1730] * v[382];
v[1753] = v[1720] + v[1731];
v[1754] = v[148] * v[1702] - v[1728] * v[380] + v[1737] * v[381] - v[1729] * v[382];
v[1755] = v[148] * v[1692] - v[1741] * v[380] + v[1745] * v[381] - v[1733] * v[382];
v[1756] = v[148] * v[1691] - v[1742] * v[380] + v[1744] * v[381] - v[1732] * v[382];
v[1757] = v[1727] + v[1743];
v[1758] = v[1717] + v[1746];
v[7406] = 0e0;
v[7407] = 0e0;
v[7408] = 0e0;
v[7409] = -v[1739] / 2e0 - v[1757];
v[7410] = v[1735] / 2e0 - v[1758];
v[7411] = -v[1750] / 2e0 - v[1753];
v[7412] = 0e0;
v[7413] = 0e0;
v[7414] = 0e0;
v[7415] = 0e0;
v[7416] = 0e0;
v[7417] = 0e0;
v[1759] = 1e0 / (v[253] * v[253]);
v[4767] = -(v[1759] * v[263]);
v[4557] = -(v[1759] * v[262]);
v[4556] = -(v[1759] * v[264]);
v[4553] = -(v[1759] * v[256]);
v[4552] = -(v[1759] * v[251]);
v[4551] = -(v[1405] * v[1759]);
v[4550] = -(v[1759] * v[252]);
v[4549] = -(v[1759] * v[4445]);
v[3190] = -(v[1759] * (v[153] * v[254] + v[4443]));
v[3189] = -(v[1759] * (v[152] * v[254] + v[4444]));
v[3188] = -(v[1759] * v[4543]);
v[3187] = -(v[1759] * v[4447]);
v[3186] = -(v[1759] * v[4544]);
v[3185] = -(v[1759] * v[4545]);
v[3184] = -(v[1759] * v[254]);
v[4735] = v[256] * v[3184];
v[3183] = -(v[1759] * v[257]);
v[4734] = v[262] * v[3183];
v[4733] = -(v[1759] * v[264] * v[266]);
v[3181] = -(v[1759] * v[4546]);
v[3180] = -(v[1759] * v[4547]);
v[3179] = -(v[1759] * v[4548]);
v[3178] = v[1576] * v[4549];
v[3177] = v[1586] * v[4550];
v[3176] = v[1587] * v[4552];
v[3077] = v[1407] * v[4549];
v[3076] = v[4448] * v[4551];
v[3074] = v[1409] * v[3184];
v[4711] = v[3074] + v[3077];
v[5030] = v[3076] + v[4711];
v[3064] = v[1409] * v[4550];
v[3061] = v[1407] * v[3183];
v[4710] = v[3061] + v[3064];
v[3060] = v[4450] * v[4551];
v[5029] = v[3060] + v[4710];
v[3050] = v[1409] * v[4552];
v[3048] = v[1407] * v[4767];
v[3047] = v[266] * v[4551];
v[4708] = v[3047] + v[3050];
v[5028] = v[3048] + v[4708];
v[2943] = v[160] * v[4553];
v[2940] = v[155] * v[4553];
v[2937] = -(v[1759] * v[4554]);
v[2936] = -(v[1759] * v[4555]);
v[2931] = v[158] * v[4556];
v[2930] = v[157] * v[4557];
v[4174] = v[2930] + v[2940];
v[4166] = v[2931] + v[4174];
v[4152] = v[2930] + v[2931];
v[2928] = v[151] * v[4553];
v[4171] = v[2928] + v[2936] + v[2937];
v[4164] = -v[2937] + v[4171];
v[4154] = -v[2936] + v[4171];
v[2924] = v[163] * v[4556];
v[4177] = v[2924] + v[2943];
v[2923] = v[162] * v[4557];
v[4160] = v[2923] + v[2924];
v[4156] = v[2923] + v[4177];
v[2257] = v[262] * v[3185] + v[256] * v[3186] + v[4733];
v[2252] = v[264] * v[3187] + v[256] * v[3188] + v[4734];
v[2247] = v[262] * v[3189] + v[264] * v[3190] + v[4735];
v[2045] = v[153] * v[4549];
v[2043] = v[158] * v[4550];
v[2039] = v[162] * v[4552];
v[3175] = v[1572] + v[1579] + v[1585] + v[1587] * v[2039] + v[1586] * v[2043] + v[1576] * v[2045] + v[1409] * v[2247]
+ v[1407] * v[2252] + v[1405] * v[2257];
v[1760] = v[1719] - v[1722] - v[1728] + v[1730] + v[1742] - v[1745] + v[1735] * v[369] - v[1739] * v[370]
- v[1750] * v[372] + v[1715] * v[388] + v[1714] * v[394] + v[1704] * v[398] + v[1702] * v[407] + v[1692] * v[411]
+ v[1691] * v[415] + v[1690] * v[4558] + v[1703] * v[4559] + v[1716] * v[4560] + v[3175] * v[4585] + v[1753] * v[4730]
+ v[1757] * v[4731] + v[1758] * v[4732];
v[1762] = (-2e0*v[1546] + 2e0*v[1551] + v[1640] * v[4561] + v[1639] * v[4562] + v[1644] * v[4563] + v[1642] * v[4564]
	+ v[1658] * v[4565] + v[1657] * v[4566] - v[1641] * v[482] - v[1643] * v[501] - v[1659] * v[519]) / 2e0;
v[1764] = -v[1547] + v[1568] + v[1662] * v[4537] + v[1654] * v[4538] + v[1633] * v[4539] + v[1634] * v[487]
+ v[1635] * v[493] + v[1652] * v[497] + v[1653] * v[506] + v[1661] * v[510] + v[1660] * v[514];
v[1765] = (v[1619] * v[204] - v[1643] * v[479] + v[1654] * v[480] - v[1647] * v[481]) / 2e0;
v[1766] = (-2e0*v[1558] + 2e0*v[1565] + v[1637] * v[4561] + v[1638] * v[4562] + v[1645] * v[4563] + v[1646] * v[4564]
	+ v[1648] * v[4565] + v[1649] * v[4566] - v[1636] * v[482] - v[1647] * v[501] - v[1650] * v[519]) / 2e0;
v[5032] = v[1762] * v[374] - v[1764] * v[377] + v[1766] * v[379];
v[7418] = 0e0;
v[7419] = 0e0;
v[7420] = 0e0;
v[7421] = 0e0;
v[7422] = 0e0;
v[7423] = 0e0;
v[7424] = 0e0;
v[7425] = 0e0;
v[7426] = 0e0;
v[7427] = 8e0*v[1762];
v[7428] = -8e0*v[1764];
v[7429] = 8e0*v[1766];
v[1791] = -v[1765] + v[1766] * v[1865] + v[1764] * v[1868] + v[1762] * v[1870] - v[1676] * v[4428];
v[3208] = v[1791] + (-(v[1632] * v[204]) + v[1641] * v[479] - v[1633] * v[480] + v[1636] * v[481]) / 2e0;
v[1767] = (v[1606] * v[204] - v[1659] * v[479] + v[1662] * v[480] - v[1650] * v[481]) / 2e0;
v[3206] = v[1765] - v[1767] + v[3208];
v[3203] = -v[1767] + v[1791];
v[1768] = v[1667] + v[1671];
v[1769] = v[1670] + v[1672];
v[1770] = v[1656] + v[1668];
v[1771] = (-2e0*v[1570] + 2e0*v[1575] - v[1725] * v[383] - v[1727] * v[402] - v[1743] * v[420] + v[1724] * v[4567]
	+ v[1723] * v[4568] + v[1728] * v[4569] + v[1726] * v[4570] + v[1742] * v[4571] + v[1741] * v[4572]) / 2e0;
v[1773] = -v[1571] + v[1592] + v[1718] * v[388] + v[1719] * v[394] + v[1736] * v[398] + v[1737] * v[407] + v[1745] * v[411]
+ v[1744] * v[415] + v[1746] * v[4558] + v[1738] * v[4559] + v[1717] * v[4560];
v[1774] = (v[148] * v[1703] - v[1727] * v[380] + v[1738] * v[381] - v[1731] * v[382]) / 2e0;
v[1775] = (-2e0*v[1582] + 2e0*v[1589] - v[1720] * v[383] - v[1731] * v[402] - v[1734] * v[420] + v[1721] * v[4567]
	+ v[1722] * v[4568] + v[1729] * v[4569] + v[1730] * v[4570] + v[1732] * v[4571] + v[1733] * v[4572]) / 2e0;
v[5036] = v[1771] * v[368] - v[1773] * v[371] + v[1775] * v[373];
v[7502] = 0e0;
v[7503] = 0e0;
v[7504] = 0e0;
v[7505] = 8e0*v[1771];
v[7506] = -8e0*v[1773];
v[7507] = 8e0*v[1775];
v[7508] = 0e0;
v[7509] = 0e0;
v[7510] = 0e0;
v[7511] = 0e0;
v[7512] = 0e0;
v[7513] = 0e0;
v[1783] = -v[1774] + v[1775] * v[1839] + v[1773] * v[1842] + v[1771] * v[1844] - v[1760] * v[4426];
v[3201] = v[1783] + (-(v[148] * v[1716]) + v[1725] * v[380] - v[1717] * v[381] + v[1720] * v[382]) / 2e0;
v[1776] = (v[148] * v[1690] - v[1743] * v[380] + v[1746] * v[381] - v[1734] * v[382]) / 2e0;
v[3199] = v[1774] - v[1776] + v[3201];
v[3196] = -v[1776] + v[1783];
v[1777] = v[1751] + v[1755];
v[1778] = v[1754] + v[1756];
v[1779] = v[1740] + v[1752];
v[6730] = v[1460] + v[4573] + v[4574] + v[4575] - v[4576];
v[6731] = v[1468] + v[4577] + v[4578] + v[4579] - v[4580];
v[6732] = v[1476] + v[4581] + v[4582] + v[4583] - v[4584];
v[6733] = -v[1754] + v[1756] + v[3196] * v[368] + v[1779] * v[369] + v[1777] * v[372] + v[1739] * v[4430] + 2e0*
(v[1771] * v[4426] + v[1757] * v[4430]) + v[10] * (v[1706] + v[1407] * v[4446] + v[1405] * v[4449]) - v[1157] * v[773]
- v[1158] * v[774] - v[1159] * v[775];
v[6734] = v[1751] - v[1755] + v[10] * (v[1695] + v[1405] * v[1784] + v[1409] * v[1785]) + v[1779] * v[370]
+ v[3199] * v[371] + v[1778] * v[372] - v[1735] * v[4430] + 2e0*(-(v[1773] * v[4426]) + v[1758] * v[4430])
- v[1161] * v[773] - v[1162] * v[774] - v[1163] * v[775];
v[6735] = -v[1740] + v[1752] + v[10] * (v[1682] + v[1407] * v[1786] + v[1409] * v[1787]) + v[1778] * v[369]
+ v[1777] * v[370] + v[3201] * v[373] + v[1750] * v[4430] + 2e0*(v[1775] * v[4426] + v[1753] * v[4430]) - v[1165] * v[773]
- v[1166] * v[774] - v[1167] * v[775];
v[6736] = v[1456] - v[4573] - v[4574] - v[4575] + v[4576];
v[6737] = v[1464] - v[4577] - v[4578] - v[4579] + v[4580];
v[6738] = v[1472] - v[4581] - v[4582] - v[4583] + v[4584];
v[6739] = -v[1670] + v[1672] + v[3203] * v[374] + v[1770] * v[375] + v[1768] * v[378] + v[1655] * v[4437] + 2e0*
(v[1762] * v[4428] + v[1673] * v[4437]) + v[10] * (v[1622] + v[1401] * v[4454] + v[1399] * v[4457]) - v[1181] * v[773]
- v[1182] * v[774] - v[1183] * v[775];
v[6740] = v[1667] - v[1671] + v[10] * (v[1611] + v[1399] * v[1792] + v[1403] * v[1793]) + v[1770] * v[376]
+ v[3206] * v[377] + v[1769] * v[378] - v[1651] * v[4437] + 2e0*(-(v[1764] * v[4428]) + v[1674] * v[4437])
- v[1185] * v[773] - v[1186] * v[774] - v[1187] * v[775];
v[6741] = -v[1656] + v[1668] + v[10] * (v[1598] + v[1401] * v[1794] + v[1403] * v[1795]) + v[1769] * v[375]
+ v[1768] * v[376] + v[3208] * v[379] + v[1666] * v[4437] + 2e0*(v[1766] * v[4428] + v[1669] * v[4437]) - v[1189] * v[773]
- v[1190] * v[774] - v[1191] * v[775];
for (i1194 = 1; i1194 <= 12; i1194++) {
	i4604 = (i1194 == 9 ? 1 : 0);
	i4603 = (i1194 == 3 ? 1 : 0);
	i4618 = i4603 - i4604;
	i4602 = (i1194 == 8 ? 1 : 0);
	i4601 = (i1194 == 2 ? 1 : 0);
	i4619 = i4601 - i4602;
	i4600 = (i1194 == 7 ? 1 : 0);
	i4599 = (i1194 == 1 ? 1 : 0);
	i4620 = i4599 - i4600;
	v[1802] = v[6233 + i1194];
	v[1804] = v[6257 + i1194];
	v[1806] = v[6245 + i1194];
	v[1807] = v[6745 + i1194];
	v[1809] = v[6221 + i1194];
	v[1850] = v[1809] * v[4426];
	v[1915] = -(v[1850] * v[4585]);
	v[4709] = v[1915] * v[264];
	v[4707] = v[1915] * v[262];
	v[4597] = v[1915] * v[253];
	v[1810] = v[6793 + i1194];
	v[1812] = v[6281 + i1194];
	v[1814] = v[6305 + i1194];
	v[1816] = v[6293 + i1194];
	v[1817] = v[6817 + i1194];
	v[1819] = v[6269 + i1194];
	v[1876] = v[1819] * v[4428];
	v[1991] = -(v[1876] * v[4586]);
	v[4701] = v[1991] * v[290];
	v[4699] = v[1991] * v[288];
	v[4598] = v[1991] * v[279];
	v[1820] = v[6865 + i1194];
	v[1827] = (i1194 == 12 ? 1 : 0);
	v[4691] = v[10] * v[1827];
	v[1828] = (i1194 == 11 ? 1 : 0);
	v[4689] = v[10] * v[1828];
	v[1829] = (i1194 == 10 ? 1 : 0);
	v[4690] = v[10] * v[1829];
	v[1830] = (i1194 == 6 ? 1 : 0);
	v[4696] = v[10] * v[1830];
	v[1831] = (i1194 == 5 ? 1 : 0);
	v[4694] = v[10] * v[1831];
	v[1832] = (i1194 == 4 ? 1 : 0);
	v[4695] = v[10] * v[1832];
	v[1833] = v[1802] + v[1830];
	v[4724] = 2e0*v[1833];
	v[1834] = v[1802] - v[1830];
	v[4725] = 2e0*v[1834];
	v[1835] = v[1804] + v[1832];
	v[4726] = 2e0*v[1835];
	v[1836] = v[1804] - v[1832];
	v[4727] = 2e0*v[1836];
	v[1837] = v[1806] - v[1831];
	v[4728] = 2e0*v[1837];
	v[1838] = v[1806] + v[1831];
	v[4729] = 2e0*v[1838];
	v[1840] = v[1809] * v[1839] + v[1830] * v[4587];
	v[1841] = -v[1809] + v[1831] * v[371];
	v[1843] = v[1809] * v[1842] - v[1831] * v[4587];
	v[1845] = v[1809] * v[1844] + v[1832] * v[4587];
	v[4589] = 2e0*(-(v[148] * v[1831]) - v[1809] * v[4429]);
	v[1847] = -(v[148] * v[1832]) - v[1809] * v[4431];
	v[1848] = -(v[148] * v[1830]) - v[1809] * v[4432];
	v[1849] = v[1850] * v[372] + v[1830] * v[4430];
	v[1851] = v[1850] * v[370] + v[1832] * v[4430];
	v[1852] = -(v[1850] * v[369]) - v[1831] * v[4430];
	v[1853] = (v[148] * v[1807] - v[1850] * v[420]) / 2e0;
	v[1951] = v[1853] * v[264];
	v[1854] = (-(v[1807] * v[382]) - v[1840] * v[420]) / 2e0;
	v[1855] = (v[148] * v[1841] - v[1850] * v[402]) / 2e0;
	v[1944] = v[1855] * v[262];
	v[1856] = (v[1841] * v[381] + v[1843] * v[402]) / 2e0;
	v[1857] = (v[148] * v[1810] - v[1850] * v[383]) / 2e0;
	v[1936] = v[1857] * v[256];
	v[1858] = (-(v[1810] * v[380]) - v[1845] * v[383]) / 2e0;
	v[1859] = v[1812] + v[1827];
	v[4712] = 2e0*v[1859];
	v[1860] = v[1812] - v[1827];
	v[4713] = 2e0*v[1860];
	v[1861] = v[1814] + v[1829];
	v[4714] = 2e0*v[1861];
	v[1862] = v[1814] - v[1829];
	v[4715] = 2e0*v[1862];
	v[1863] = v[1816] - v[1828];
	v[4716] = 2e0*v[1863];
	v[1864] = v[1816] + v[1828];
	v[4717] = 2e0*v[1864];
	v[1866] = v[1819] * v[1865] + v[1827] * v[4588];
	v[1867] = -v[1819] + v[1828] * v[377];
	v[1869] = v[1819] * v[1868] - v[1828] * v[4588];
	v[1871] = v[1819] * v[1870] + v[1829] * v[4588];
	v[4593] = 2e0*(-(v[1828] * v[204]) - v[1819] * v[4436]);
	v[1873] = -(v[1829] * v[204]) - v[1819] * v[4438];
	v[1874] = -(v[1827] * v[204]) - v[1819] * v[4439];
	v[1875] = v[1876] * v[378] + v[1827] * v[4437];
	v[1877] = v[1876] * v[376] + v[1829] * v[4437];
	v[1878] = -(v[1876] * v[375]) - v[1828] * v[4437];
	v[1879] = (v[1817] * v[204] - v[1876] * v[519]) / 2e0;
	v[2027] = v[1879] * v[290];
	v[1880] = (-(v[1817] * v[481]) - v[1866] * v[519]) / 2e0;
	v[1881] = (v[1867] * v[204] - v[1876] * v[501]) / 2e0;
	v[2020] = v[1881] * v[288];
	v[1882] = (v[1867] * v[480] + v[1869] * v[501]) / 2e0;
	v[1883] = (v[1820] * v[204] - v[1876] * v[482]) / 2e0;
	v[2012] = v[1883] * v[282];
	v[1884] = (-(v[1820] * v[479]) - v[1871] * v[482]) / 2e0;
	v[1885] = (v[1807] * v[381] + v[1843] * v[420] + v[4589]) / 2e0;
	v[1886] = (v[1810] * v[381] + v[1843] * v[383] + v[4589]) / 2e0;
	v[1887] = v[1847] + v[1807] * v[4431] - v[1845] * v[4558];
	v[1888] = v[1847] + v[1841] * v[4431] - v[1845] * v[4559];
	v[1889] = -v[1850] - v[1835] * v[380] - v[1845] * v[415];
	v[1890] = v[148] * v[1835] - v[1850] * v[415];
	v[1891] = v[1850] + v[1837] * v[381] + v[1843] * v[411];
	v[1892] = v[148] * v[1837] - v[1850] * v[411];
	v[1953] = v[1892] * v[256];
	v[1893] = v[1850] - v[1836] * v[380] - v[1845] * v[407];
	v[1894] = v[148] * v[1836] - v[1850] * v[407];
	v[1945] = v[1894] * v[264];
	v[1895] = v[1848] + v[1841] * v[4432] - v[1840] * v[4559];
	v[1896] = v[1848] + v[1810] * v[4432] - v[1840] * v[4560];
	v[1897] = -v[1850] - v[1833] * v[382] - v[1840] * v[398];
	v[1898] = v[148] * v[1833] - v[1850] * v[398];
	v[1899] = -v[1850] + v[1838] * v[381] + v[1843] * v[394];
	v[1900] = v[148] * v[1838] - v[1850] * v[394];
	v[1937] = v[1900] * v[264];
	v[1901] = -v[1849] + v[1835] * v[381] + v[1843] * v[415];
	v[1902] = -v[1849] - v[1837] * v[380] - v[1845] * v[411];
	v[1903] = -v[1849] + v[1836] * v[381] + v[1843] * v[407];
	v[1904] = -v[1849] - v[1838] * v[380] - v[1845] * v[394];
	v[1905] = v[1915] + v[1849] * v[4590];
	v[1906] = v[1885] * v[877] + v[1901] * v[879] + v[1891] * v[881];
	v[1907] = v[1887] * v[877] + v[1889] * v[879] + v[1902] * v[881];
	v[1908] = v[1850] - v[1834] * v[382] - v[1840] * v[388];
	v[1909] = v[148] * v[1834] - v[1850] * v[388];
	v[1910] = -v[1851] + v[1833] * v[381] + v[1843] * v[398];
	v[1911] = -v[1851] - v[1837] * v[382] - v[1840] * v[411];
	v[1912] = -v[1851] - v[1838] * v[382] - v[1840] * v[394];
	v[1913] = -v[1851] + v[1834] * v[381] + v[1843] * v[388];
	v[1914] = v[1849] * v[390] + v[1851] * v[393];
	v[1916] = v[1915] + v[1851] * v[4591];
	v[2057] = v[1409] * v[1916];
	v[1917] = v[1903] * v[877] + v[1856] * v[879] + v[1910] * v[881];
	v[1918] = v[1852] - v[1835] * v[382] - v[1840] * v[415];
	v[1919] = v[1852] - v[1836] * v[382] - v[1840] * v[407];
	v[1920] = v[1852] - v[1833] * v[380] - v[1845] * v[398];
	v[1921] = v[1852] - v[1834] * v[380] - v[1845] * v[388];
	v[1922] = -(v[1851] * v[387]) - v[1852] * v[390];
	v[1923] = -(v[1849] * v[387]) - v[1852] * v[393];
	v[1924] = v[1915] + v[1852] * v[4592];
	v[2061] = v[1407] * v[1924];
	v[1925] = v[1854] * v[877] + v[1918] * v[879] + v[1911] * v[881];
	v[1926] = v[1919] * v[877] + v[1895] * v[879] + v[1897] * v[881];
	v[1927] = v[1893] * v[877] + v[1888] * v[879] + v[1920] * v[881];
	v[1928] = v[1904] * v[877] + v[1921] * v[879] + v[1858] * v[881];
	v[1929] = v[1912] * v[877] + v[1908] * v[879] + v[1896] * v[881];
	v[1930] = v[1433] * v[1854] + v[1436] * v[1885] + v[1435] * v[1887] + v[1431] * v[1893] + v[1428] * v[1899]
		+ v[1434] * v[1903] + v[1430] * v[1904] + v[1429] * v[1912] + v[1432] * v[1919];
	v[1931] = v[1434] * v[1856] + v[1431] * v[1888] + v[1435] * v[1889] + v[1432] * v[1895] + v[1436] * v[1901]
		+ v[1429] * v[1908] + v[1428] * v[1913] + v[1433] * v[1918] + v[1430] * v[1921];
	v[1932] = v[1899] * v[877] + v[1913] * v[879] + v[1886] * v[881];
	v[1933] = v[1430] * v[1858] + v[1428] * v[1886] + v[1436] * v[1891] + v[1429] * v[1896] + v[1432] * v[1897]
		+ v[1435] * v[1902] + v[1434] * v[1910] + v[1433] * v[1911] + v[1431] * v[1920];
	v[1934] = v[1936] + v[1909] * v[262];
	v[1935] = v[1934] + v[1937] + v[4695];
	v[1938] = v[1936] + v[1937];
	v[1939] = v[1857] * v[49] + v[1909] * v[52] + v[1900] * v[55];
	v[1940] = v[1857] * v[50] + v[1909] * v[53] + v[1900] * v[56];
	v[1941] = v[1857] * v[51] + v[1909] * v[54] + v[1900] * v[57];
	v[1942] = v[1944] + v[1898] * v[256];
	v[1943] = v[1942] + v[1945] + v[4694];
	v[1946] = v[1944] + v[1945];
	v[1947] = v[1898] * v[49] + v[1855] * v[52] + v[1894] * v[55];
	v[1948] = v[1898] * v[50] + v[1855] * v[53] + v[1894] * v[56];
	v[1949] = v[1898] * v[51] + v[1855] * v[54] + v[1894] * v[57];
	v[1950] = v[1951] + v[1953];
	v[1952] = v[1951] + v[1890] * v[262];
	v[1954] = v[1952] + v[1953] + v[4696];
	v[1955] = v[1892] * v[49] + v[1890] * v[52] + v[1853] * v[55];
	v[1956] = v[1892] * v[50] + v[1890] * v[53] + v[1853] * v[56];
	v[1957] = v[1892] * v[51] + v[1890] * v[54] + v[1853] * v[57];
	v[1958] = i4599 * v[10];
	v[1959] = i4601 * v[10];
	v[1960] = i4603 * v[10];
	v[1961] = (v[4593] + v[1817] * v[480] + v[1869] * v[519]) / 2e0;
	v[1962] = (v[4593] + v[1820] * v[480] + v[1869] * v[482]) / 2e0;
	v[1963] = v[1873] + v[1817] * v[4438] - v[1871] * v[4537];
	v[1964] = v[1873] + v[1867] * v[4438] - v[1871] * v[4538];
	v[1965] = -v[1876] - v[1861] * v[479] - v[1871] * v[514];
	v[1966] = v[1861] * v[204] - v[1876] * v[514];
	v[1967] = v[1876] + v[1863] * v[480] + v[1869] * v[510];
	v[1968] = v[1863] * v[204] - v[1876] * v[510];
	v[2029] = v[1968] * v[282];
	v[1969] = v[1876] - v[1862] * v[479] - v[1871] * v[506];
	v[1970] = v[1862] * v[204] - v[1876] * v[506];
	v[2021] = v[1970] * v[290];
	v[1971] = v[1874] + v[1867] * v[4439] - v[1866] * v[4538];
	v[1972] = v[1874] + v[1820] * v[4439] - v[1866] * v[4539];
	v[1973] = -v[1876] - v[1859] * v[481] - v[1866] * v[497];
	v[1974] = v[1859] * v[204] - v[1876] * v[497];
	v[1975] = -v[1876] + v[1864] * v[480] + v[1869] * v[493];
	v[1976] = v[1864] * v[204] - v[1876] * v[493];
	v[2013] = v[1976] * v[290];
	v[1977] = -v[1875] + v[1861] * v[480] + v[1869] * v[514];
	v[1978] = -v[1875] - v[1863] * v[479] - v[1871] * v[510];
	v[1979] = -v[1875] + v[1862] * v[480] + v[1869] * v[506];
	v[1980] = -v[1875] - v[1864] * v[479] - v[1871] * v[493];
	v[1981] = v[1991] + v[1875] * v[4594];
	v[1982] = -(v[1961] * v[871]) - v[1977] * v[873] - v[1967] * v[875];
	v[1983] = -(v[1963] * v[871]) - v[1965] * v[873] - v[1978] * v[875];
	v[1984] = v[1876] - v[1860] * v[481] - v[1866] * v[487];
	v[1985] = v[1860] * v[204] - v[1876] * v[487];
	v[1986] = -v[1877] + v[1859] * v[480] + v[1869] * v[497];
	v[1987] = -v[1877] - v[1863] * v[481] - v[1866] * v[510];
	v[1988] = -v[1877] - v[1864] * v[481] - v[1866] * v[493];
	v[1989] = -v[1877] + v[1860] * v[480] + v[1869] * v[487];
	v[1990] = v[1875] * v[489] + v[1877] * v[492];
	v[1992] = v[1991] + v[1877] * v[4595];
	v[2085] = v[1403] * v[1992];
	v[1993] = -(v[1979] * v[871]) - v[1882] * v[873] - v[1986] * v[875];
	v[1994] = v[1878] - v[1861] * v[481] - v[1866] * v[514];
	v[1995] = v[1878] - v[1862] * v[481] - v[1866] * v[506];
	v[1996] = v[1878] - v[1859] * v[479] - v[1871] * v[497];
	v[1997] = v[1878] - v[1860] * v[479] - v[1871] * v[487];
	v[1998] = -(v[1877] * v[486]) - v[1878] * v[489];
	v[1999] = -(v[1875] * v[486]) - v[1878] * v[492];
	v[2000] = v[1991] + v[1878] * v[4596];
	v[2089] = v[1401] * v[2000];
	v[2001] = -(v[1880] * v[871]) - v[1994] * v[873] - v[1987] * v[875];
	v[2002] = -(v[1995] * v[871]) - v[1971] * v[873] - v[1973] * v[875];
	v[2003] = -(v[1969] * v[871]) - v[1964] * v[873] - v[1996] * v[875];
	v[2004] = -(v[1980] * v[871]) - v[1997] * v[873] - v[1884] * v[875];
	v[2005] = -(v[1988] * v[871]) - v[1984] * v[873] - v[1972] * v[875];
	v[2006] = -(v[1424] * v[1880]) - v[1427] * v[1961] - v[1426] * v[1963] - v[1422] * v[1969] - v[1419] * v[1975]
		- v[1425] * v[1979] - v[1421] * v[1980] - v[1420] * v[1988] - v[1423] * v[1995];
	v[2007] = -(v[1425] * v[1882]) - v[1422] * v[1964] - v[1426] * v[1965] - v[1423] * v[1971] - v[1427] * v[1977]
		- v[1420] * v[1984] - v[1419] * v[1989] - v[1424] * v[1994] - v[1421] * v[1997];
	v[2008] = -(v[1975] * v[871]) - v[1989] * v[873] - v[1962] * v[875];
	v[2009] = -(v[1421] * v[1884]) - v[1419] * v[1962] - v[1427] * v[1967] - v[1420] * v[1972] - v[1423] * v[1973]
		- v[1426] * v[1978] - v[1425] * v[1986] - v[1424] * v[1987] - v[1422] * v[1996];
	v[2010] = v[2012] + v[1985] * v[288];
	v[2011] = v[2010] + v[2013] + v[4690];
	v[2014] = v[2012] + v[2013];
	v[2015] = v[1883] * v[92] + v[1985] * v[95] + v[1976] * v[98];
	v[2016] = v[1883] * v[93] + v[1985] * v[96] + v[1976] * v[99];
	v[2017] = v[100] * v[1976] + v[1883] * v[94] + v[1985] * v[97];
	v[2018] = v[2020] + v[1974] * v[282];
	v[2019] = v[2018] + v[2021] + v[4689];
	v[2022] = v[2020] + v[2021];
	v[2023] = v[1974] * v[92] + v[1881] * v[95] + v[1970] * v[98];
	v[2024] = v[1974] * v[93] + v[1881] * v[96] + v[1970] * v[99];
	v[2025] = v[100] * v[1970] + v[1974] * v[94] + v[1881] * v[97];
	v[2026] = v[2027] + v[2029];
	v[2028] = v[2027] + v[1966] * v[288];
	v[2030] = v[2028] + v[2029] + v[4691];
	v[2031] = v[1968] * v[92] + v[1966] * v[95] + v[1879] * v[98];
	v[2032] = v[1968] * v[93] + v[1966] * v[96] + v[1879] * v[99];
	v[2033] = v[100] * v[1879] + v[1968] * v[94] + v[1966] * v[97];
	v[2034] = i4600 * v[10];
	v[2035] = i4602 * v[10];
	v[2036] = i4604 * v[10];
	v[2037] = v[1843] + v[1914];
	v[2055] = v[1409] * v[2037];
	v[2038] = -v[1843] + v[1914];
	v[2040] = (v[162] * v[2037] + v[1890] * v[251] + v[2039] * v[4597]) / v[253];
	v[2041] = v[1840] + v[1922];
	v[2063] = v[1409] * v[2041];
	v[2042] = -v[1840] + v[1922];
	v[2059] = v[1407] * v[2042];
	v[2044] = (v[158] * v[2041] + v[1894] * v[252] + v[2043] * v[4597]) / v[253];
	v[2046] = (v[153] * v[2042] + v[1900] * v[4445] + v[2045] * v[4597]) / v[253];
	v[2047] = v[2057] + v[2059];
	v[4745] = v[2047] / v[253];
	v[2048] = v[1845] + v[1923];
	v[2053] = v[1407] * v[2048];
	v[2049] = -v[1845] + v[1923];
	v[2050] = v[2061] + v[2063];
	v[4739] = v[2050] / v[253];
	v[2051] = v[1405] * v[1905] + v[2053] + v[2055];
	v[4749] = v[2051] / v[253];
	v[2054] = v[2051] - v[2055];
	v[2056] = v[2051] - v[2053];
	v[4740] = v[2056] / v[253];
	v[2058] = v[1405] * v[2038] + v[2057];
	v[2060] = v[2058] + v[2059];
	v[4741] = v[2060] / v[253];
	v[2062] = v[1405] * v[2049] + v[2061];
	v[2064] = v[2062] + v[2063];
	v[4746] = v[2064] / v[253];
	v[2065] = v[1869] + v[1990];
	v[2083] = v[1403] * v[2065];
	v[2066] = -v[1869] + v[1990];
	v[2068] = (v[2065] * v[218] + v[1966] * v[277] + v[2067] * v[4598]) / v[279];
	v[2069] = v[1866] + v[1998];
	v[2091] = v[1403] * v[2069];
	v[2070] = -v[1866] + v[1998];
	v[2087] = v[1401] * v[2070];
	v[2072] = (v[2069] * v[214] + v[1970] * v[278] + v[2071] * v[4598]) / v[279];
	v[2074] = (v[2070] * v[209] + v[1976] * v[4453] + v[2073] * v[4598]) / v[279];
	v[2075] = v[2085] + v[2087];
	v[4759] = v[2075] / v[279];
	v[2076] = v[1871] + v[1999];
	v[2081] = v[1401] * v[2076];
	v[2077] = -v[1871] + v[1999];
	v[2078] = v[2089] + v[2091];
	v[4753] = v[2078] / v[279];
	v[2079] = v[1399] * v[1981] + v[2081] + v[2083];
	v[4763] = v[2079] / v[279];
	v[2082] = v[2079] - v[2083];
	v[2084] = v[2079] - v[2081];
	v[4754] = v[2084] / v[279];
	v[2086] = v[1399] * v[2066] + v[2085];
	v[2088] = v[2086] + v[2087];
	v[4755] = v[2088] / v[279];
	v[2090] = v[1399] * v[2077] + v[2089];
	v[2092] = v[2090] + v[2091];
	v[4760] = v[2092] / v[279];
	v[2094] = i4620 + v[142] * v[1939] + v[143] * v[1940] + v[144] * v[1941] - v[198] * v[2015] - v[199] * v[2016]
		- v[200] * v[2017];
	v[2097] = i4619 + v[142] * v[1947] + v[143] * v[1948] + v[144] * v[1949] - v[198] * v[2023] - v[199] * v[2024]
		- v[200] * v[2025];
	v[2100] = i4618 + v[142] * v[1955] + v[143] * v[1956] + v[144] * v[1957] - v[198] * v[2031] - v[199] * v[2032]
		- v[200] * v[2033];
	v[4605] = v[2094] * v[302] + v[2097] * v[303] + v[2100] * v[304];
	v[2104] = v[4605] / v[732];
	v[2102] = -(v[1014] * v[1267] * v[4605]);
	v[2103] = v[2104] * v[315];
	v[2105] = v[2104];
	v[2106] = v[2103] * v[302] + v[2094] * v[316];
	v[2107] = v[2103] * v[303] + v[2097] * v[316];
	v[2108] = v[2103] * v[304] + v[2100] * v[316];
	b2109 = b329;
	if (b2109) {
		v[2110] = -v[2106];
		v[2111] = -v[2107];
		v[2112] = -v[2108];
	}
	else {
		v[2110] = v[2106];
		v[2111] = v[2107];
		v[2112] = v[2108];
	};
	v[2113] = v[2112] * v[4606];
	v[2114] = 2e0*v[1293] * v[2112];
	v[2115] = v[2111] * v[4607];
	v[2116] = 2e0*v[1333] * v[2111];
	v[2117] = v[2110] * v[4608];
	v[2118] = 2e0*v[1381] * v[2110];
	v[2119] = 0e0;
	v[2120] = 0e0;
	v[2121] = 0e0;
	v[2122] = 0e0;
	v[2123] = 0e0;
	v[2124] = 0e0;
	v[2125] = 0e0;
	v[2126] = 0e0;
	v[2127] = 0e0;
	v[2128] = 0e0;
	v[2129] = 0e0;
	v[2130] = 0e0;
	v[2131] = 0e0;
	v[2132] = 0e0;
	v[2133] = 0e0;
	v[2134] = 0e0;
	v[2135] = 0e0;
	v[2136] = 0e0;
	v[2137] = 0e0;
	v[2138] = 0e0;
	v[2139] = 0e0;
	v[2140] = 0e0;
	v[2141] = 0e0;
	v[2142] = 0e0;
	v[2143] = 0e0;
	v[2144] = 0e0;
	v[2145] = 0e0;
	v[2146] = 0e0;
	v[2147] = 0e0;
	v[2148] = 0e0;
	v[2149] = 0e0;
	v[2150] = 0e0;
	b2151 = b678;
	if (b2151) {
		v[2152] = v[2112] * v[327];
		v[2153] = -(v[2112] * v[326]);
		v[2154] = v[2152] - v[2111] * v[328];
		v[2155] = v[2111] * v[326];
		v[2156] = v[2153] + v[2110] * v[328];
		v[2157] = v[2155] - v[2110] * v[327];
		v[4612] = v[1525] * v[2154] + v[1524] * v[2156] + v[1522] * v[2157];
		v[4610] = v[2154] * v[680] + v[2156] * v[681] + v[2157] * v[682];
		v[2158] = v[4610] / v[683];
		v[4611] = v[2158] * v[4986];
		v[2168] = v[2158] * v[4609];
		v[2122] = -(v[1528] * v[2159] * v[4610]);
		v[2160] = v[2154] * v[4475] + v[4611] * v[680];
		v[2175] = v[2160] * v[4658];
		v[2162] = v[2156] * v[4475] + v[4611] * v[681];
		v[2179] = 2e0*v[2162] * v[693];
		v[2163] = v[2157] * v[4475] + v[4611] * v[682];
		v[2176] = v[2163] * v[4655];
		v[2125] = v[2168] * v[4516] + v[4612] * v[685];
		v[2124] = v[2158] * v[2166] * v[689];
		v[2123] = v[2158] * v[4516] * v[690] + v[4612] * v[691];
		v[2149] = v[2168] * v[2169] * v[685];
		v[2150] = v[1527] * v[2158] * v[2164] * v[2169] * v[4987];
		v[2121] = v[2157] * v[4517] + v[1522] * v[4611];
		v[2120] = v[2156] * v[4517] + v[1524] * v[4611];
		v[2119] = v[2154] * v[4517] + v[1525] * v[4611];
		v[2171] = (v[2162] * v[692] + v[2160] * v[693]) / 2e0;
		v[2172] = v[2175] + v[2179];
		v[2173] = v[2172] + v[2176];
		v[2174] = (v[2163] * v[692] + v[2160] * v[694]) / 2e0;
		v[2177] = v[2175] + v[2176];
		v[2178] = (v[2163] * v[693] + v[2162] * v[694]) / 2e0;
		v[2180] = v[2176] + v[2179];
		v[2128] = (v[1508] * v[2160] + v[1495] * v[2162] + 4e0*v[2163] * v[2181]) / 2e0;
		v[2127] = (v[1515] * v[2160] + v[1495] * v[2163] + v[2162] * v[4613]) / 2e0;
		v[2126] = (v[1515] * v[2162] + v[1508] * v[2163] + v[2160] * v[4614]) / 2e0;
		v[2184] = v[2173] * v[4654];
		v[2148] = 8e0*v[1517] * v[2173] * v[4988];
		v[2147] = v[2184] * v[4615];
		v[2185] = v[2163] + v[2171];
		v[2186] = v[2163] - v[2171];
		v[2146] = v[1438] * v[2184];
		v[2187] = -v[2162] + v[2174];
		v[2188] = v[2162] + v[2174];
		v[2145] = v[1439] * v[2184];
		v[2133] = v[1502] * v[2184] + v[2185] * v[695];
		v[2144] = v[1440] * v[2184];
		v[2134] = (v[1498] * v[2184] - v[2177] * v[695]) / 2e0;
		v[2143] = v[2184] * v[4616];
		v[2189] = v[2160] + v[2178];
		v[2190] = -v[2160] + v[2178];
		v[2142] = v[1442] * v[2184];
		v[2136] = v[1489] * v[2184] + v[2187] * v[695];
		v[2141] = v[1443] * v[2184];
		v[2137] = v[1485] * v[2184] + v[2189] * v[695];
		v[2140] = v[1444] * v[2184];
		v[2138] = (v[1481] * v[2184] - v[2172] * v[695]) / 2e0;
		v[2139] = v[2184] * v[4617];
		v[2135] = v[1493] * v[2184] + v[2190] * v[695];
		v[2132] = v[1506] * v[2184] + v[2188] * v[695];
		v[2131] = v[1511] * v[2184] - v[2186] * v[695];
		v[2130] = (v[1516] * v[2184] - v[2180] * v[695]) / 2e0;
		v[2129] = v[1440] * v[2185] - v[1438] * v[2186] + v[1443] * v[2187] + v[1439] * v[2188] + v[1444] * v[2189]
			+ v[1442] * v[2190] - v[2180] * v[4615] - v[2177] * v[4616] - v[2172] * v[4617];
	}
	else {
	};
	v[2191] = 0e0;
	v[2192] = 0e0;
	v[2193] = 0e0;
	v[2194] = 0e0;
	v[2195] = 0e0;
	v[2196] = 0e0;
	b2197 = b30;
	if (b2197) {
		v[4623] = -(v[2110] * v[713]);
		v[4622] = -(v[2111] * v[714]);
		v[4621] = -(v[2112] * v[715]);
		v[2237] = v[1263] * v[2110];
		v[2206] = i4618 + v[145] * v[1955] + v[146] * v[1956] + v[147] * v[1957] - v[201] * v[2031] - v[202] * v[2032]
			- v[203] * v[2033] + v[2138] * v[3];
		v[2138] = 0e0;
		v[2207] = v[2] * v[2137] + v[2206];
		v[2137] = 0e0;
		v[2208] = v[1] * v[2136] + v[2207];
		v[4625] = -(v[2208] * v[333]);
		v[2136] = 0e0;
		v[2217] = i4619 + v[145] * v[1947] + v[146] * v[1948] + v[147] * v[1949] - v[201] * v[2023] - v[202] * v[2024]
			- v[2025] * v[203] + v[2135] * v[3];
		v[2135] = 0e0;
		v[2218] = v[2] * v[2134] + v[2217];
		v[2134] = 0e0;
		v[2219] = v[1] * v[2133] + v[2218];
		v[4624] = -(v[2219] * v[332]);
		v[2133] = 0e0;
		v[2228] = i4620 + v[145] * v[1939] + v[146] * v[1940] + v[147] * v[1941] - v[201] * v[2015] - v[2016] * v[202]
			- v[2017] * v[203] + v[2132] * v[3];
		v[2132] = 0e0;
		v[2229] = v[2] * v[2131] + v[2228];
		v[2131] = 0e0;
		v[2230] = v[1] * v[2130] + v[2229];
		v[4626] = -(v[2230] * v[331]);
		v[2130] = 0e0;
		v[2231] = v[1261] * v[2112];
		v[2193] = v[2112] * v[2241];
		v[2118] = v[2118] + v[1263] * v[4621];
		v[2116] = v[2116] + v[1262] * v[4621];
		v[2233] = v[1262] * v[2111];
		v[2234] = v[2231] + v[2233];
		v[2192] = v[2111] * v[2243];
		v[2118] = v[2118] + v[1263] * v[4622];
		v[2114] = v[2114] + v[1261] * v[4622];
		v[2236] = v[2233] + v[2237];
		v[2238] = v[2231] + v[2237];
		v[2191] = v[2110] * v[2242];
		v[2116] = v[2116] + v[1262] * v[4623];
		v[2114] = v[2114] + v[1261] * v[4623];
		v[2191] = -(v[1263] * v[2117]) + v[2191];
		v[2196] = v[1261] * v[2208];
		v[2195] = v[1262] * v[2219];
		v[2194] = v[1263] * v[2230];
		v[2192] = -(v[1262] * v[2115]) + v[2192];
		v[2193] = -(v[1261] * v[2113]) + v[2193];
		v[2193] = v[2193] - v[2236] * v[333];
		v[2114] = v[2114] + v[2208] * v[2241] + v[1261] * (v[4624] + v[4626]) - v[2236] * v[715];
		v[2191] = v[2191] - v[2234] * v[331];
		v[2118] = v[2118] + v[2230] * v[2242] + v[1263] * (v[4624] + v[4625]) - v[2234] * v[713];
		v[2192] = v[2192] - v[2238] * v[332];
		v[2116] = v[2116] + v[2219] * v[2243] + v[1262] * (v[4625] + v[4626]) - v[2238] * v[714];
	}
	else {
	};
	v[4630] = v[2110] * v[2344];
	v[4629] = v[2111] * v[2299];
	v[2301] = v[2111] * v[4509];
	v[4628] = v[1266] * v[2112];
	v[4627] = v[2112] * v[2279];
	v[2284] = v[2112] * v[4509];
	v[2244] = -v[1960] + v[2036] - v[1907] * v[255] - v[1906] * v[265] - v[1925] * v[275] + v[1983] * v[281]
		+ v[1982] * v[291] + v[2001] * v[301];
	v[4635] = v[2244] * v[333];
	v[2245] = -v[1959] + v[2035] - v[1927] * v[255] - v[1917] * v[265] - v[1926] * v[275] + v[2003] * v[281]
		+ v[1993] * v[291] + v[2002] * v[301];
	v[4636] = -(v[2245] * v[332]);
	v[4633] = v[1265] * v[2244] + v[1266] * v[2245];
	v[2246] = -v[1958] + v[2034] - v[1928] * v[255] - v[1932] * v[265] - v[1929] * v[275] + v[2004] * v[281]
		+ v[2008] * v[291] + v[2005] * v[301];
	v[4634] = v[1264] * v[2245] + v[1265] * v[2246];
	v[4632] = v[1264] * v[2244] + v[1266] * v[2246];
	v[4631] = -(v[2246] * v[331]);
	v[2251] = v[1915] * v[2247] + v[1916] * v[2248] + v[2037] * v[2249] + v[2041] * v[2250] + v[2040] * v[262]
		+ v[2044] * v[264] + v[1935] * v[3195] + v[1942] * v[4446] + v[1950] * v[4449] + v[7117 + i1194];
	v[2256] = v[1785] * v[1934] + v[1784] * v[1952] + v[1915] * v[2252] + v[1924] * v[2253] + v[2042] * v[2254]
		+ v[2048] * v[2255] + v[2040] * v[256] + v[2046] * v[264] + v[1943] * v[3198] + v[7141 + i1194];
	v[2261] = v[1787] * v[1938] + v[1786] * v[1946] + v[1915] * v[2257] + v[1905] * v[2258] + v[2038] * v[2259]
		+ v[2049] * v[2260] + v[2044] * v[256] + v[2046] * v[262] + v[1954] * v[3200] + v[7165 + i1194];
	v[2266] = v[1991] * v[2262] + v[1992] * v[2263] + v[2065] * v[2264] + v[2069] * v[2265] + v[2068] * v[288]
		+ v[2072] * v[290] + v[2011] * v[3202] + v[2018] * v[4454] + v[2026] * v[4457] + v[7189 + i1194];
	v[2271] = v[1793] * v[2010] + v[1792] * v[2028] + v[1991] * v[2267] + v[2000] * v[2268] + v[2070] * v[2269]
		+ v[2076] * v[2270] + v[2068] * v[282] + v[2074] * v[290] + v[2019] * v[3205] + v[7213 + i1194];
	v[2276] = v[1795] * v[2014] + v[1794] * v[2022] + v[1991] * v[2272] + v[1981] * v[2273] + v[2066] * v[2274]
		+ v[2077] * v[2275] + v[2072] * v[282] + v[2074] * v[288] + v[2030] * v[3207] + v[7237 + i1194];
	v[2278] = v[1408] * v[2251] + v[1406] * v[2256] + v[1404] * v[2261] + v[1402] * v[2266] + v[1400] * v[2271]
		+ v[1398] * v[2276] + v[2112] * v[2277];
	v[2118] = v[2118] + v[1264] * v[4627];
	v[2362] = v[2112] * v[4510];
	v[2278] = v[2278] + v[2111] * v[2298];
	v[2118] = v[2118] + v[1264] * v[4629];
	v[2278] = v[2278] + v[2110] * v[2315];
	v[2316] = v[1264] * v[2110];
	v[2318] = v[2316] * v[255];
	v[2319] = v[2316] * v[265];
	v[2320] = v[2316] * v[275];
	v[2321] = v[2316] * v[281];
	v[2322] = v[2316] * v[291];
	v[2323] = v[2316] * v[301];
	v[2421] = v[2110] * v[4510];
	v[2278] = v[2278] + v[2117] * v[2344];
	v[2360] = v[1352] * v[2251] + v[1350] * v[2256] + v[1348] * v[2261] + v[1346] * v[2266] + v[1344] * v[2271]
		+ v[1342] * v[2276] + v[2112] * v[2359];
	v[2116] = v[2116] + v[1265] * v[4627];
	v[2458] = v[301] * v[4628];
	v[2456] = v[291] * v[4628];
	v[2454] = v[281] * v[4628];
	v[2452] = v[275] * v[4628];
	v[2450] = v[265] * v[4628];
	v[2448] = v[255] * v[4628];
	v[2360] = v[2360] + v[2111] * v[2376];
	v[2377] = v[1265] * v[2111];
	v[2379] = v[2377] * v[255];
	v[2380] = v[2377] * v[265];
	v[2381] = v[2377] * v[275];
	v[2382] = v[2377] * v[281];
	v[2383] = v[2377] * v[291];
	v[2384] = v[2377] * v[301];
	v[2385] = v[2316] + v[2377];
	v[2387] = v[2318] + v[2379];
	v[2388] = v[2319] + v[2380];
	v[2389] = v[2320] + v[2381];
	v[2390] = v[2321] + v[2382];
	v[2391] = v[2322] + v[2383];
	v[2392] = v[2323] + v[2384];
	v[2488] = v[2111] * v[4511];
	v[2360] = v[2115] * v[2299] + v[2360];
	v[2360] = v[2360] + v[2110] * v[2419];
	v[2116] = v[2116] + v[1265] * v[4630];
	v[2504] = v[2110] * v[4511];
	v[2436] = v[1304] * v[2251] + v[1302] * v[2256] + v[1300] * v[2261] + v[1298] * v[2266] + v[1296] * v[2271]
		+ v[1294] * v[2276] + v[2112] * v[2435];
	v[2437] = v[2377] + v[4628];
	v[2439] = v[2379] + v[2448];
	v[2440] = v[2380] + v[2450];
	v[2441] = v[2381] + v[2452];
	v[2442] = v[2382] + v[2454];
	v[2443] = v[2383] + v[2456];
	v[2444] = v[2384] + v[2458];
	v[2446] = v[2316] + v[4628];
	v[2449] = v[2318] + v[2448];
	v[2451] = v[2319] + v[2450];
	v[2453] = v[2320] + v[2452];
	v[2455] = v[2321] + v[2454];
	v[2457] = v[2322] + v[2456];
	v[2459] = v[2323] + v[2458];
	v[2436] = v[2113] * v[2279] + v[2436];
	v[2478] = -(v[1266] * v[2113]) + v[2284] + v[2362];
	v[2436] = v[2436] + v[2111] * v[2486];
	v[2114] = v[2114] + v[1266] * v[4629];
	v[2494] = -(v[1265] * v[2115]) + v[2301] + v[2488];
	v[2436] = v[2436] + v[2110] * v[2502];
	v[2114] = v[2114] + v[1266] * v[4630];
	v[2510] = -(v[1264] * v[2117]) + v[2421] + v[2504];
	v[2518] = 0e0;
	b2519 = b733;
	if (b2519) {
		b2520 = b735;
		if (b2520) {
			v[2518] = -v[2105];
			v[2110] = 0e0;
			v[2111] = 0e0;
			v[2112] = 0e0;
		}
		else {
		};
	}
	else {
	};
	v[2436] = v[2436] + v[333] * (v[4631] + v[4636]) - v[2244] * v[722];
	v[2538] = v[2436];
	v[2360] = v[2360] - v[332] * (-v[4631] + v[4635]) - v[2245] * v[720];
	v[2537] = v[2360];
	v[2114] = v[2114] - v[331] * v[4632] - v[332] * v[4633];
	v[2692] = v[2114];
	v[2118] = v[2118] - v[333] * v[4632] - v[332] * v[4634];
	v[2696] = v[2118];
	v[2116] = v[2116] - v[333] * v[4633] - v[331] * v[4634];
	v[2694] = v[2116];
	v[2278] = v[2278] - v[331] * (v[4635] - v[4636]) - v[2246] * v[718];
	v[2536] = v[2278];
	v[2521] = 0e0;
	v[2522] = 0e0;
	v[2523] = 0e0;
	v[2524] = 0e0;
	v[2525] = 0e0;
	v[2526] = 0e0;
	v[2527] = 0e0;
	v[2528] = 0e0;
	b2529 = b733;
	if (b2529) {
		v[4637] = v[1270] * v[2536];
		v[2539] = v[1271] * v[2537] + v[1272] * v[2538] + v[4637];
		b2530 = b735;
		if (b2530) {
			v[2677] = Power(v[744], v[2531]);
			v[2676] = (v[1278] * v[745]) / v[1279];
			v[2533] = v[2518] * v[2676] * v[7];
			v[2532] = v[2533] * v[2677];
			v[2527] = -((v[1277] * v[2532]) / v[1279]);
			v[2524] = v[1277] * v[2531] * v[2533] * Power(v[744], -2e0 + v[745]);
			v[2518] = 0e0;
			v[2278] = 0e0;
			v[2360] = 0e0;
			v[2525] = v[2539];
			v[2436] = 0e0;
		}
		else {
			v[2685] = Power(v[732], v[2534]);
			v[2684] = (v[1286] * v[750]) / v[1287];
			v[2535] = v[2105] * v[2684] * v[7];
			v[2532] = v[2535] * v[2685];
			v[2528] = -((v[1277] * v[2532]) / v[1287]);
			v[2102] = v[2102] + v[1277] * v[2534] * v[2535] * Power(v[732], -2e0 + v[750]);
			v[2105] = 0e0;
			v[2278] = 0e0;
			v[2360] = 0e0;
			v[2526] = v[2539];
			v[2436] = 0e0;
		};
		v[2521] = v[1270] * v[2532];
		v[2522] = v[1271] * v[2532];
		v[2523] = v[1272] * v[2532];
	}
	else {
	};
	v[2683] = v[2521];
	v[2681] = v[2522];
	v[2679] = v[2523];
	v[2541] = (v[1196] * (v[1939] * v[342] + v[1940] * v[344] + v[1941] * v[346]) + v[1197] * (v[1947] * v[342]
		+ v[1948] * v[344] + v[1949] * v[346]) + v[1198] * (v[1955] * v[342] + v[1956] * v[344] + v[1957] * v[346]))*v[9];
	v[2542] = (v[1196] * (v[1939] * v[334] + v[1940] * v[335] + v[1941] * v[336]) + v[1197] * (v[1947] * v[334]
		+ v[1948] * v[335] + v[1949] * v[336]) + v[1198] * (v[1955] * v[334] + v[1956] * v[335] + v[1957] * v[336]))*v[9];
	v[2543] = (v[1196] * (v[2015] * v[359] + v[2016] * v[361] + v[2017] * v[363]) + v[1197] * (v[2023] * v[359]
		+ v[2024] * v[361] + v[2025] * v[363]) + v[1198] * (v[2031] * v[359] + v[2032] * v[361] + v[2033] * v[363]))*v[9];
	v[2544] = (v[1196] * (v[2015] * v[351] + v[2016] * v[352] + v[2017] * v[353]) + v[1197] * (v[2023] * v[351]
		+ v[2024] * v[352] + v[2025] * v[353]) + v[1198] * (v[2031] * v[351] + v[2032] * v[352] + v[2033] * v[353]))*v[9];
	v[2545] = (-i4599 + i4600)*v[1145] + (-i4601 + i4602)*v[1149] + (-i4603 + i4604)*v[1153] - v[1189] * v[1827]
		- v[1185] * v[1828] - v[1181] * v[1829] - v[1165] * v[1830] - v[1161] * v[1831] - v[1157] * v[1832];
	v[2561] = v[2545];
	v[2546] = (-i4599 + i4600)*v[1146] + (-i4601 + i4602)*v[1150] + (-i4603 + i4604)*v[1154] - v[1190] * v[1827]
		- v[1186] * v[1828] - v[1182] * v[1829] - v[1166] * v[1830] - v[1162] * v[1831] - v[1158] * v[1832];
	v[2559] = v[2546];
	v[2547] = (-i4599 + i4600)*v[1147] + (-i4601 + i4602)*v[1151] + (-i4603 + i4604)*v[1155] - v[1191] * v[1827]
		- v[1187] * v[1828] - v[1183] * v[1829] - v[1167] * v[1830] - v[1163] * v[1831] - v[1159] * v[1832];
	v[2557] = v[2547];
	v[2548] = 0e0;
	v[2549] = 0e0;
	v[2550] = 0e0;
	v[2551] = 0e0;
	v[2552] = 0e0;
	v[2553] = 0e0;
	b2554 = b733;
	if (b2554) {
		b2555 = b29;
		if (b2555) {
			b2556 = b771;
			if (b2556) {
				v[2553] = v[2547];
				v[2547] = 0e0;
				v[2552] = v[2546];
				v[2546] = 0e0;
				v[2551] = v[2545];
				v[2545] = 0e0;
			}
			else {
				v[2569] = v[2561] * v[4487];
				v[2567] = v[2559] * v[4487];
				v[2565] = v[2557] * v[4487];
				v[2547] = 0e0;
				v[2546] = 0e0;
				v[4639] = (v[6] * (v[2561] * v[783] + v[2559] * v[784] + v[2557] * v[785])) / sqrt(v[2570]);
				v[2545] = 0e0;
				v[4638] = ((v[2569] * v[766] + v[2567] * v[767] + v[2565] * v[768])*v[781]) / sqrt(v[2566]);
				v[2553] = v[4638] * v[768] + v[2565] * v[782];
				v[2552] = v[4638] * v[767] + v[2567] * v[782];
				v[2551] = v[4638] * v[766] + v[2569] * v[782];
				v[2550] = v[4639] * v[759];
				v[2549] = v[4639] * v[758];
				v[2548] = v[4639] * v[757];
			};
		}
		else {
			b2572 = b805;
			if (b2572) {
				v[2553] = v[2557];
				v[2547] = 0e0;
				v[2552] = v[2559];
				v[2546] = 0e0;
				v[2551] = v[2561];
				v[2545] = 0e0;
			}
			else {
				v[2578] = v[2561] * v[6] * v[816];
				v[2576] = v[2559] * v[6] * v[816];
				v[2575] = v[2557] * v[6] * v[816];
				v[2547] = 0e0;
				v[2546] = 0e0;
				v[4641] = (v[6] * (v[2561] * v[813] + v[2559] * v[814] + v[2557] * v[815])) / sqrt(v[2570]);
				v[2545] = 0e0;
				v[4640] = ((v[2578] * v[766] + v[2576] * v[767] + v[2575] * v[768])*v[811]) / sqrt(v[2566]);
				v[2553] = v[4640] * v[768] + v[2575] * v[812];
				v[2552] = v[4640] * v[767] + v[2576] * v[812];
				v[2551] = v[4640] * v[766] + v[2578] * v[812];
				v[2550] = v[4641] * v[759];
				v[2549] = v[4641] * v[758];
				v[2548] = v[4641] * v[757];
			};
		};
	}
	else {
	};
	v[2584] = -(v[1827] * v[775]) - (v[1198] * v[2276] + v[2553] * v[301])*v[9];
	v[2585] = -(v[1828] * v[775]) - (v[1198] * v[2271] + v[2553] * v[291])*v[9];
	v[2586] = -(v[1829] * v[775]) - (v[1198] * v[2266] + v[2553] * v[281])*v[9];
	v[2590] = -(v[1830] * v[775]) - (v[1198] * v[2261] + v[2553] * v[275])*v[9];
	v[2591] = -(v[1831] * v[775]) - (v[1198] * v[2256] + v[2553] * v[265])*v[9];
	v[2592] = -(v[1832] * v[775]) - (v[1198] * v[2251] + v[255] * v[2553])*v[9];
	v[4668] = (i4603 - i4604)*v[775] + (v[1198] * (v[1960] - v[2036]) + (v[246] - v[249])*v[2553])*v[9];
	v[4669] = (i4601 - i4602)*v[775] + (v[1198] * (v[1959] - v[2035]) + (v[245] - v[248])*v[2553])*v[9];
	v[4670] = (i4599 - i4600)*v[775] + (v[1198] * (v[1958] - v[2034]) + (v[244] - v[247])*v[2553])*v[9];
	v[2609] = -(v[1827] * v[774]) - (v[1197] * v[2276] + v[2552] * v[301])*v[9];
	v[2610] = -(v[1828] * v[774]) - (v[1197] * v[2271] + v[2552] * v[291])*v[9];
	v[2611] = -(v[1829] * v[774]) - (v[1197] * v[2266] + v[2552] * v[281])*v[9];
	v[2615] = -(v[1830] * v[774]) - (v[1197] * v[2261] + v[2552] * v[275])*v[9];
	v[2616] = -(v[1831] * v[774]) - (v[1197] * v[2256] + v[2552] * v[265])*v[9];
	v[2617] = -(v[1832] * v[774]) - (v[1197] * v[2251] + v[255] * v[2552])*v[9];
	v[4665] = (i4603 - i4604)*v[774] + (v[1197] * (v[1960] - v[2036]) + (v[246] - v[249])*v[2552])*v[9];
	v[4666] = (i4601 - i4602)*v[774] + (v[1197] * (v[1959] - v[2035]) + (v[245] - v[248])*v[2552])*v[9];
	v[4667] = (i4599 - i4600)*v[774] + (v[1197] * (v[1958] - v[2034]) + (v[244] - v[247])*v[2552])*v[9];
	v[2634] = -(v[1827] * v[773]) - (v[1196] * v[2276] + v[2551] * v[301])*v[9];
	v[2635] = -(v[1828] * v[773]) - (v[1196] * v[2271] + v[2551] * v[291])*v[9];
	v[2636] = -(v[1829] * v[773]) - (v[1196] * v[2266] + v[2551] * v[281])*v[9];
	v[2640] = -(v[1830] * v[773]) - (v[1196] * v[2261] + v[2551] * v[275])*v[9];
	v[2641] = -(v[1831] * v[773]) - (v[1196] * v[2256] + v[2551] * v[265])*v[9];
	v[2642] = -(v[1832] * v[773]) - (v[1196] * v[2251] + v[255] * v[2551])*v[9];
	v[4662] = (i4603 - i4604)*v[773] + (v[1196] * (v[1960] - v[2036]) + (v[246] - v[249])*v[2551])*v[9];
	v[4663] = (i4601 - i4602)*v[773] + (v[1196] * (v[1959] - v[2035]) + (v[245] - v[248])*v[2551])*v[9];
	v[4664] = (i4599 - i4600)*v[773] + (v[1196] * (v[1958] - v[2034]) + (v[244] - v[247])*v[2551])*v[9];
	v[2658] = v[2553] * v[4];
	v[2659] = v[2552] * v[4];
	v[2660] = v[2551] * v[4];
	v[2661] = v[2550];
	v[2662] = v[2550];
	v[2663] = v[2549];
	v[2664] = v[2549];
	v[2665] = v[2548];
	v[2666] = v[2548];
	b2667 = b733;
	if (b2667) {
		v[2668] = 0e0;
		v[2669] = 0e0;
		v[2670] = 0e0;
		b2671 = b752;
		if (b2671) {
			v[2670] = v[2661];
			v[2661] = 0e0;
			v[2669] = v[2663];
			v[2663] = 0e0;
			v[2668] = v[2665];
			v[2665] = 0e0;
		}
		else {
			v[2662] = 0e0;
			v[2661] = 0e0;
			v[2664] = 0e0;
			v[2663] = 0e0;
			v[2666] = 0e0;
			v[2665] = 0e0;
		};
		v[2682] = v[2668] * v[719];
		v[2680] = v[2669] * v[721];
		v[2678] = v[2670] * v[723];
		b2675 = b735;
		if (b2675) {
			v[2525] = v[2525] + v[2678];
			v[2523] = v[2523] + v[2670] * v[747];
			v[2525] = v[2525] + v[2680];
			v[2522] = v[2522] + v[2669] * v[747];
			v[2525] = v[2525] + v[2682];
			v[2521] = v[2521] + v[2668] * v[747];
			v[2527] = v[2527] + v[2525] * v[4486];
			v[2524] = v[2524] + (v[2527] * v[2676] * v[2677]) / 2e0;
		}
		else {
			v[2526] = v[2526] + v[2678];
			v[2523] = v[2679] + v[2670] * v[751];
			v[2526] = v[2526] + v[2680];
			v[2522] = v[2681] + v[2669] * v[751];
			v[2526] = v[2526] + v[2682];
			v[2521] = v[2683] + v[2668] * v[751];
			v[2528] = v[2528] + v[2526] * v[4486];
			v[2102] = v[2102] + (v[2528] * v[2684] * v[2685]) / 2e0;
		};
	}
	else {
	};
	v[4642] = -(v[2521] * v[331]);
	v[4643] = -(v[2522] * v[332]);
	v[4644] = -(v[2523] * v[333]);
	v[4706] = v[10] * (v[2510] + v[331] * (-v[2437] + v[4643] + v[4644]) + v[2542] * v[624] + v[2541] * v[639]
		- v[2544] * v[654] - v[2543] * v[669] - v[2521] * v[718] + (v[1145] * v[2551] + v[1146] * v[2552] + v[1147] * v[2553]
			)*v[9]);
	v[4705] = v[10] * (v[2494] + v[332] * (-v[2446] + v[4642] + v[4644]) + v[2542] * v[626] + v[2541] * v[641]
		- v[2544] * v[656] - v[2543] * v[671] - v[2522] * v[720] + (v[1149] * v[2551] + v[1150] * v[2552] + v[1151] * v[2553]
			)*v[9]);
	v[4704] = v[10] * (v[2478] + v[333] * (-v[2385] + v[4642] + v[4643]) + v[2542] * v[628] + v[2541] * v[643]
		- v[2544] * v[658] - v[2543] * v[673] - v[2523] * v[722] + (v[1153] * v[2551] + v[1154] * v[2552] + v[1155] * v[2553]
			)*v[9]);
	v[2697] = v[2102];
	v[2695] = v[2666];
	v[2693] = v[2664];
	v[2691] = v[2662];
	b2686 = b733;
	if (b2686) {
		v[4647] = v[2666] * v[331];
		v[4646] = v[2664] * v[332];
		v[4645] = v[2662] * v[333];
		b2687 = b735;
		if (b2687) {
			v[2114] = v[2114] + v[2662] * v[4484];
			v[2662] = 0e0;
			v[2116] = v[2116] + v[2664] * v[4484];
			v[2664] = 0e0;
			v[2118] = v[2118] + v[2666] * v[4484];
			v[2666] = 0e0;
			v[2524] = v[2524] + v[31] * v[34] * (-v[4645] - v[4646] - v[4647])*Power(v[744], v[745]);
			v[2102] = v[2102] - v[2524];
		}
		else {
			v[2114] = v[2692] + v[2691] * v[741];
			v[2662] = 0e0;
			v[2116] = v[2694] + v[2693] * v[741];
			v[2664] = 0e0;
			v[2118] = v[2696] + v[2695] * v[741];
			v[2666] = 0e0;
			v[2102] = v[2697] + v[4493] * (v[4645] + v[4646] + v[4647])*Power(v[732], v[750]);
		};
	}
	else {
	};
	v[2698] = v[1266] * v[2276] + v[2523] * v[301];
	v[4676] = v[2698] * v[333];
	v[2699] = v[1266] * v[2271] + v[2523] * v[291];
	v[4675] = v[2699] * v[333];
	v[2700] = v[1266] * v[2266] + v[2523] * v[281];
	v[4674] = v[2700] * v[333];
	v[2701] = v[1266] * v[2261] + v[2523] * v[275];
	v[4673] = v[2701] * v[333];
	v[2702] = v[1266] * v[2256] + v[2523] * v[265];
	v[4672] = v[2702] * v[333];
	v[2703] = v[1266] * v[2251] + v[2523] * v[255];
	v[4671] = v[2703] * v[333];
	v[2118] = v[2118] + v[2523] * v[3239];
	v[2116] = v[2116] + v[2523] * v[3240];
	v[2708] = -(v[1266] * v[2244]) + v[2523] * v[3241];
	v[2114] = v[2114] + v[2523] * v[2709];
	v[2725] = v[1265] * v[2276] + v[2522] * v[301];
	v[4687] = v[2725] * v[332];
	v[2726] = v[1265] * v[2271] + v[2522] * v[291];
	v[4685] = v[2726] * v[332];
	v[2727] = v[1265] * v[2266] + v[2522] * v[281];
	v[4683] = v[2727] * v[332];
	v[2728] = v[1265] * v[2261] + v[2522] * v[275];
	v[4681] = v[2728] * v[332];
	v[2729] = v[1265] * v[2256] + v[2522] * v[265];
	v[4679] = v[2729] * v[332];
	v[2730] = v[1265] * v[2251] + v[2522] * v[255];
	v[4677] = v[2730] * v[332];
	v[2118] = v[2118] + v[2522] * v[3261];
	v[2733] = -(v[1265] * v[2245]) + v[2522] * v[3262];
	v[2116] = v[2116] + v[2522] * v[2734];
	v[2114] = v[2114] + v[2522] * v[3264];
	v[2752] = v[1264] * v[2276] + v[2521] * v[301];
	v[4688] = v[2752] * v[331];
	v[2753] = v[1264] * v[2271] + v[2521] * v[291];
	v[4686] = v[2753] * v[331];
	v[2754] = v[1264] * v[2266] + v[2521] * v[281];
	v[4684] = v[2754] * v[331];
	v[2755] = v[1264] * v[2261] + v[2521] * v[275];
	v[4682] = v[2755] * v[331];
	v[2756] = v[1264] * v[2256] + v[2521] * v[265];
	v[4680] = v[2756] * v[331];
	v[2757] = v[1264] * v[2251] + v[2521] * v[255];
	v[4678] = v[2757] * v[331];
	v[2758] = -(v[1264] * v[2246]) + v[2521] * v[3283];
	v[2118] = v[2118] + v[2521] * v[2759];
	v[2116] = v[2116] + v[2521] * v[3285];
	v[2114] = v[2114] + v[2521] * v[3286];
	v[2779] = 0e0;
	v[2780] = 0e0;
	v[2781] = 0e0;
	v[2782] = 0e0;
	v[2783] = 0e0;
	v[2784] = 0e0;
	v[2785] = 0e0;
	v[2786] = 0e0;
	v[2787] = 0e0;
	b2788 = b30;
	if (b2788) {
		v[4650] = v[2660] * v[331];
		v[4649] = v[2659] * v[332];
		v[4653] = v[4649] + v[4650];
		v[4648] = v[2658] * v[333];
		v[4652] = v[4648] + v[4650];
		v[4651] = v[4648] + v[4649];
		v[2196] = v[2196] + v[2658] * v[715];
		v[2195] = v[2195] + v[2659] * v[714];
		v[2191] = v[2191] + v[1447] * v[2660] - v[331] * v[4651];
		v[2192] = v[2192] + v[1449] * v[2659] - v[332] * v[4652];
		v[2193] = v[2193] + v[1451] * v[2658] - v[333] * v[4653];
		v[2194] = v[2194] + v[2660] * v[713];
		v[2118] = v[2118] - v[2660] * v[4512] - v[4651] * v[713];
		v[2116] = v[2116] + v[2659] * v[4513] - v[4652] * v[714];
		v[2114] = v[2114] + v[2658] * v[4514] - v[4653] * v[715];
		v[2779] = v[1] * v[2191];
		v[2780] = v[2] * v[2191];
		v[2781] = v[2191] * v[3];
		v[2789] = -v[2191];
		v[2790] = -(v[203] * v[2191]);
		v[2791] = -(v[202] * v[2191]);
		v[2792] = -(v[201] * v[2191]);
		v[2793] = v[2191];
		v[2794] = v[147] * v[2191];
		v[2795] = v[146] * v[2191];
		v[2796] = v[145] * v[2191];
		v[2782] = v[1] * v[2192];
		v[2783] = v[2] * v[2192];
		v[2784] = v[2192] * v[3];
		v[2797] = -v[2192];
		v[2798] = -(v[203] * v[2192]);
		v[2799] = -(v[202] * v[2192]);
		v[2800] = -(v[201] * v[2192]);
		v[2801] = v[2192];
		v[2802] = v[147] * v[2192];
		v[2803] = v[146] * v[2192];
		v[2804] = v[145] * v[2192];
		v[2785] = v[1] * v[2193];
		v[2786] = v[2] * v[2193];
		v[2787] = v[2193] * v[3];
		v[2805] = -v[2193];
		v[2806] = -(v[203] * v[2193]);
		v[2807] = -(v[202] * v[2193]);
		v[2808] = -(v[201] * v[2193]);
		v[2809] = v[2193];
		v[2810] = v[147] * v[2193];
		v[2811] = v[146] * v[2193];
		v[2812] = v[145] * v[2193];
		v[2758] = -v[2194] + v[2758];
		v[2733] = -v[2195] + v[2733];
		v[2708] = -v[2196] + v[2708];
	}
	else {
		v[2796] = 0e0;
		v[2795] = 0e0;
		v[2794] = 0e0;
		v[2804] = 0e0;
		v[2803] = 0e0;
		v[2802] = 0e0;
		v[2812] = 0e0;
		v[2811] = 0e0;
		v[2810] = 0e0;
		v[2793] = 0e0;
		v[2801] = 0e0;
		v[2809] = 0e0;
		v[2792] = 0e0;
		v[2791] = 0e0;
		v[2790] = 0e0;
		v[2800] = 0e0;
		v[2799] = 0e0;
		v[2798] = 0e0;
		v[2808] = 0e0;
		v[2807] = 0e0;
		v[2806] = 0e0;
		v[2789] = 0e0;
		v[2797] = 0e0;
		v[2805] = 0e0;
	};
	b2813 = b678;
	if (b2813) {
		v[2129] = v[2129] + (v[1481] * v[2787]) / 2e0;
		v[2139] = v[2139] + v[2787] * v[4477];
		v[2129] = v[2129] + v[1485] * v[2786];
		v[2140] = v[2140] + v[2786] * v[695];
		v[2129] = v[2129] + v[1489] * v[2785];
		v[2141] = v[2141] + v[2785] * v[695];
		v[2129] = v[2129] + v[1493] * v[2784];
		v[2142] = v[2142] + v[2784] * v[695];
		v[2129] = v[2129] + (v[1498] * v[2783]) / 2e0;
		v[2143] = v[2143] + v[2783] * v[4477];
		v[2129] = v[2129] + v[1502] * v[2782];
		v[2144] = v[2144] + v[2782] * v[695];
		v[2129] = v[2129] + v[1506] * v[2781];
		v[2145] = v[2145] + v[2781] * v[695];
		v[2129] = v[2129] + v[1511] * v[2780];
		v[2146] = v[2146] + v[2780] * v[695];
		v[2129] = v[2129] + (v[1516] * v[2779]) / 2e0;
		v[2147] = v[2147] + v[2779] * v[4477];
		v[2148] = v[2148] + v[2129] * v[4654];
		v[4659] = -v[2143] + v[2148];
		v[2127] = v[2127] - v[2141];
		v[2823] = v[2141] + v[2145];
		v[2127] = v[2127] + v[2145];
		v[2126] = v[2126] + v[2140] + (v[2823] * v[694]) / 2e0;
		v[2825] = v[2140] + v[2142];
		v[2126] = v[2126] - v[2142];
		v[2128] = v[2128] + v[2144] + v[2825] * v[4476] + v[2823] * v[4515] + v[4655] * (-v[2147] + v[4659]);
		v[2128] = v[2128] - v[2146];
		v[4656] = v[2128] * v[682];
		v[2827] = v[2144] + v[2146];
		v[2125] = v[2125] + v[4656] * v[685];
		v[2123] = v[2123] + v[4656] * v[691];
		v[2121] = v[2121] + v[2128] * v[4475];
		v[2127] = v[2127] + (v[2827] * v[692] - 4e0*(v[2139] + v[2147] - v[2148])*v[693] + v[2825] * v[694]) / 2e0;
		v[4657] = v[2127] * v[681];
		v[2125] = v[2125] + v[4657] * v[685];
		v[2123] = v[2123] + v[4657] * v[691];
		v[2120] = v[2120] + v[2127] * v[4475];
		v[2126] = v[2126] + v[2827] * v[4476] + v[4658] * (-v[2139] + v[4659]);
		v[4660] = v[2126] * v[680];
		v[2125] = v[2125] + v[4660] * v[685];
		v[2123] = v[2123] + v[4660] * v[691];
		v[2119] = v[2119] + v[2126] * v[4475];
		v[2124] = v[2124] + v[2125] * v[690];
		v[2122] = v[2122] + v[2124];
		v[2149] = v[2149] + 2e0*v[2123] * v[2164];
		v[2150] = v[2150] + (v[2149] * v[2165]) / 2e0;
		v[2122] = v[2122] + v[2150];
		v[4661] = v[2122] / v[683];
		v[2121] = v[2121] + v[4661] * v[682];
		v[2120] = v[2120] + v[4661] * v[681];
		v[2119] = v[2119] + v[4661] * v[680];
		v[2118] = v[2118] - v[2121] * v[327];
		v[2116] = v[2116] + v[2121] * v[326];
		v[2118] = v[2118] + v[2120] * v[328];
		v[2114] = v[2114] - v[2120] * v[326];
		v[2116] = v[2116] - v[2119] * v[328];
		v[2114] = v[2114] + v[2119] * v[327];
	}
	else {
	};
	v[2833] = -(v[2642] * v[651]) - v[2641] * v[652] - v[2640] * v[653] + v[4664] * v[654] + v[4663] * v[656]
		+ v[4662] * v[658] - v[2636] * v[660] - v[2635] * v[661] - v[2634] * v[662];
	v[2834] = -(v[2642] * v[666]) - v[2641] * v[667] - v[2640] * v[668] + v[4664] * v[669] + v[4663] * v[671]
		+ v[4662] * v[673] - v[2636] * v[675] - v[2635] * v[676] - v[2634] * v[677];
	v[2835] = v[2642] * v[621] + v[2641] * v[622] + v[2640] * v[623] - v[4664] * v[624] - v[4663] * v[626] - v[4662] * v[628]
		+ v[2636] * v[630] + v[2635] * v[631] + v[2634] * v[632];
	v[2836] = v[2642] * v[636] + v[2641] * v[637] + v[2640] * v[638] - v[4664] * v[639] - v[4663] * v[641] - v[4662] * v[643]
		+ v[2636] * v[645] + v[2635] * v[646] + v[2634] * v[647];
	v[2837] = -(v[2617] * v[651]) - v[2616] * v[652] - v[2615] * v[653] + v[4667] * v[654] + v[4666] * v[656]
		+ v[4665] * v[658] - v[2611] * v[660] - v[2610] * v[661] - v[2609] * v[662];
	v[2838] = -(v[2617] * v[666]) - v[2616] * v[667] - v[2615] * v[668] + v[4667] * v[669] + v[4666] * v[671]
		+ v[4665] * v[673] - v[2611] * v[675] - v[2610] * v[676] - v[2609] * v[677];
	v[2839] = v[2617] * v[621] + v[2616] * v[622] + v[2615] * v[623] - v[4667] * v[624] - v[4666] * v[626] - v[4665] * v[628]
		+ v[2611] * v[630] + v[2610] * v[631] + v[2609] * v[632];
	v[2840] = v[2617] * v[636] + v[2616] * v[637] + v[2615] * v[638] - v[4667] * v[639] - v[4666] * v[641] - v[4665] * v[643]
		+ v[2611] * v[645] + v[2610] * v[646] + v[2609] * v[647];
	v[2841] = -(v[2592] * v[651]) - v[2591] * v[652] - v[2590] * v[653] + v[4670] * v[654] + v[4669] * v[656]
		+ v[4668] * v[658] - v[2586] * v[660] - v[2585] * v[661] - v[2584] * v[662];
	v[2842] = -(v[2592] * v[666]) - v[2591] * v[667] - v[2590] * v[668] + v[4670] * v[669] + v[4669] * v[671]
		+ v[4668] * v[673] - v[2586] * v[675] - v[2585] * v[676] - v[2584] * v[677];
	v[2843] = v[2592] * v[621] + v[2591] * v[622] + v[2590] * v[623] - v[4670] * v[624] - v[4669] * v[626] - v[4668] * v[628]
		+ v[2586] * v[630] + v[2585] * v[631] + v[2584] * v[632];
	v[2844] = v[2592] * v[636] + v[2591] * v[637] + v[2590] * v[638] - v[4670] * v[639] - v[4669] * v[641] - v[4668] * v[643]
		+ v[2586] * v[645] + v[2585] * v[646] + v[2584] * v[647];
	v[2845] = v[2001] * v[2281] + v[2002] * v[2282] + v[2005] * v[2283] + v[2316] * v[2332] + v[1398] * v[2521]
		+ v[1342] * v[2522] + v[1294] * v[2523] + v[2377] * v[4463] + v[2460] * v[4628] + v[2510] * v[571] + v[2494] * v[574]
		+ v[2478] * v[577] - v[2542] * v[632] - v[2541] * v[647] + v[2544] * v[662] + v[2543] * v[677] + (-(v[1189] * v[2551])
			- v[1190] * v[2552] - v[1191] * v[2553])*v[9];
	v[2961] = v[2845] * v[4442];
	v[2846] = v[1982] * v[2281] + v[1993] * v[2282] + v[2008] * v[2283] + v[2316] * v[2334] + v[1400] * v[2521]
		+ v[1344] * v[2522] + v[1296] * v[2523] + v[2377] * v[4465] + v[2462] * v[4628] + v[2510] * v[570] + v[2494] * v[573]
		+ v[2478] * v[576] - v[2542] * v[631] - v[2541] * v[646] + v[2544] * v[661] + v[2543] * v[676] + (-(v[1185] * v[2551])
			- v[1186] * v[2552] - v[1187] * v[2553])*v[9];
	v[4692] = v[279] * v[2846];
	v[2965] = v[2846] * v[4441];
	v[2847] = v[1983] * v[2281] + v[2003] * v[2282] + v[2004] * v[2283] + v[2316] * v[2336] + v[1402] * v[2521]
		+ v[1346] * v[2522] + v[1298] * v[2523] + v[2377] * v[4467] + v[2464] * v[4628] + v[2510] * v[569] + v[2494] * v[572]
		+ v[2478] * v[575] - v[2542] * v[630] - v[2541] * v[645] + v[2544] * v[660] + v[2543] * v[675] + (-(v[1181] * v[2551])
			- v[1182] * v[2552] - v[1183] * v[2553])*v[9];
	v[4693] = v[279] * v[2847];
	v[2968] = v[2847] * v[4440];
	v[2848] = -(v[1925] * v[2281]) - v[1926] * v[2282] - v[1929] * v[2283] + v[2316] * v[2338] + v[1404] * v[2521]
		+ v[1348] * v[2522] + v[1300] * v[2523] + v[2377] * v[4469] + v[2466] * v[4628] - v[2510] * v[472] - v[2494] * v[475]
		- v[2478] * v[478] - v[2542] * v[623] - v[2541] * v[638] + v[2544] * v[653] + v[2543] * v[668] + (-(v[1165] * v[2551])
			- v[1166] * v[2552] - v[1167] * v[2553])*v[9];
	v[3058] = v[2848] * v[4435];
	v[2849] = -(v[1906] * v[2281]) - v[1917] * v[2282] - v[1932] * v[2283] + v[2316] * v[2340] + v[1406] * v[2521]
		+ v[1350] * v[2522] + v[1302] * v[2523] + v[2377] * v[4471] + v[2468] * v[4628] - v[2510] * v[471] - v[2494] * v[474]
		- v[2478] * v[477] - v[2542] * v[622] - v[2541] * v[637] + v[2544] * v[652] + v[2543] * v[667] + (-(v[1161] * v[2551])
			- v[1162] * v[2552] - v[1163] * v[2553])*v[9];
	v[4697] = v[253] * v[2849];
	v[3062] = v[2849] * v[4434];
	v[2850] = -(v[1907] * v[2281]) - v[1927] * v[2282] - v[1928] * v[2283] + v[2316] * v[2342] + v[1408] * v[2521]
		+ v[1352] * v[2522] + v[1304] * v[2523] + v[2377] * v[4473] + v[2470] * v[4628] - v[2510] * v[470] - v[2494] * v[473]
		- v[2478] * v[476] - v[2542] * v[621] - v[2541] * v[636] + v[2544] * v[651] + v[2543] * v[666] + (-(v[1157] * v[2551])
			- v[1158] * v[2552] - v[1159] * v[2553])*v[9];
	v[4698] = v[253] * v[2850];
	v[3065] = v[2850] * v[4433];
	v[2758] = v[2758] + v[2757] * v[470] + v[2756] * v[471] + v[2755] * v[472] - v[2754] * v[569] - v[2753] * v[570]
		- v[2752] * v[571];
	v[2118] = v[2118] + v[2332] * v[2752] + v[2334] * v[2753] + v[2336] * v[2754] + v[2338] * v[2755] + v[2340] * v[2756]
		+ v[2342] * v[2757] + v[2758] * v[4608];
	v[2708] = v[2708] + v[2703] * v[476] + v[2702] * v[477] + v[2701] * v[478] - v[2700] * v[575] - v[2699] * v[576]
		- v[2698] * v[577];
	v[2851] = v[2510] * v[301] - v[331] * (v[2444] + v[4676] + v[4687]) - v[2752] * v[718];
	v[2852] = v[2510] * v[291] - v[331] * (v[2443] + v[4675] + v[4685]) - v[2753] * v[718];
	v[2853] = v[2510] * v[281] - v[331] * (v[2442] + v[4674] + v[4683]) - v[2754] * v[718];
	v[2854] = -(v[2510] * v[275]) + v[331] * (v[2441] + v[4673] + v[4681]) + v[2755] * v[718];
	v[2855] = -(v[2510] * v[265]) + v[331] * (v[2440] + v[4672] + v[4679]) + v[2756] * v[718];
	v[2118] = v[2118] + v[2437] * v[3283] + v[2439] * v[470] + v[2440] * v[471] + v[2441] * v[472] - v[2442] * v[569]
		- v[2443] * v[570] - v[2444] * v[571] + v[333] * (v[2703] * v[470] + v[2702] * v[471] + v[2701] * v[472] - v[2700] * v[569]
			- v[2699] * v[570] - v[2698] * v[571]) + v[332] * (v[2730] * v[470] + v[2729] * v[471] + v[2728] * v[472]
				- v[2727] * v[569] - v[2726] * v[570] - v[2725] * v[571]);
	v[2733] = v[2733] + v[2730] * v[473] + v[2729] * v[474] + v[2728] * v[475] - v[2727] * v[572] - v[2726] * v[573]
		- v[2725] * v[574];
	v[2856] = -(v[2510] * v[255]) + v[331] * (v[2439] + v[4671] + v[4677]) + v[2757] * v[718];
	v[2116] = v[2116] + v[2725] * v[4463] + v[2726] * v[4465] + v[2727] * v[4467] + v[2728] * v[4469] + v[2729] * v[4471]
		+ v[2730] * v[4473] + v[2733] * v[4607] + v[333] * (v[2703] * v[473] + v[2702] * v[474] + v[2701] * v[475]
			- v[2700] * v[572] - v[2699] * v[573] - v[2698] * v[574]) + v[331] * (v[2757] * v[473] + v[2756] * v[474]
				+ v[2755] * v[475] - v[2754] * v[572] - v[2753] * v[573] - v[2752] * v[574]);
	v[2857] = -(v[2494] * v[255]) + v[332] * (v[2449] + v[4671] + v[4678]) + v[2730] * v[720];
	v[2858] = -(v[2494] * v[265]) + v[332] * (v[2451] + v[4672] + v[4680]) + v[2729] * v[720];
	v[2859] = -(v[2494] * v[275]) + v[332] * (v[2453] + v[4673] + v[4682]) + v[2728] * v[720];
	v[2860] = v[2494] * v[281] - v[332] * (v[2455] + v[4674] + v[4684]) - v[2727] * v[720];
	v[2861] = v[2494] * v[291] - v[332] * (v[2457] + v[4675] + v[4686]) - v[2726] * v[720];
	v[2116] = v[2116] + v[2446] * v[3262] + v[2449] * v[473] + v[2451] * v[474] + v[2453] * v[475] - v[2455] * v[572]
		- v[2457] * v[573] - v[2459] * v[574];
	v[2862] = v[2494] * v[301] - v[332] * (v[2459] + v[4676] + v[4688]) - v[2725] * v[720];
	v[2114] = v[2114] + v[2460] * v[2698] + v[2462] * v[2699] + v[2464] * v[2700] + v[2466] * v[2701] + v[2468] * v[2702]
		+ v[2470] * v[2703] + v[2708] * v[4606] + v[332] * (v[2730] * v[476] + v[2729] * v[477] + v[2728] * v[478]
			- v[2727] * v[575] - v[2726] * v[576] - v[2725] * v[577]) + v[331] * (v[2757] * v[476] + v[2756] * v[477]
				+ v[2755] * v[478] - v[2754] * v[575] - v[2753] * v[576] - v[2752] * v[577]);
	v[2863] = -(v[2478] * v[255]) + v[333] * (v[2387] + v[4677] + v[4678]) + v[2703] * v[722];
	v[2864] = -(v[2478] * v[265]) + v[333] * (v[2388] + v[4679] + v[4680]) + v[2702] * v[722];
	v[2865] = -(v[2478] * v[275]) + v[333] * (v[2389] + v[4681] + v[4682]) + v[2701] * v[722];
	v[2866] = v[2478] * v[281] - v[333] * (v[2390] + v[4683] + v[4684]) - v[2700] * v[722];
	v[2867] = v[2478] * v[291] - v[333] * (v[2391] + v[4685] + v[4686]) - v[2699] * v[722];
	v[2114] = v[2114] + v[2385] * v[3241] + v[2387] * v[476] + v[2388] * v[477] + v[2389] * v[478] - v[2390] * v[575]
		- v[2391] * v[576] - v[2392] * v[577];
	v[2868] = v[2478] * v[301] - v[333] * (v[2392] + v[4687] + v[4688]) - v[2698] * v[722];
	b2869 = b329;
	if (b2869) {
		v[2870] = -v[2114];
		v[2871] = -v[2116];
		v[2872] = -v[2118];
	}
	else {
		v[2870] = v[2114];
		v[2871] = v[2116];
		v[2872] = v[2118];
	};
	v[2102] = v[2102] + v[1541] * v[2104] * v[314] + (v[1536] * v[2094] + v[1535] * v[2097] + v[1534] * v[2100]
		+ v[2872] * v[302] + v[2871] * v[303] + v[2870] * v[304])*v[315];
	v[2881] = v[1534] * v[2103] + v[2870] * v[316] + (v[1267] * v[2100] + v[2102] * v[304]) / v[732];
	v[2883] = v[1535] * v[2103] + v[2871] * v[316] + (v[1267] * v[2097] + v[2102] * v[303]) / v[732];
	v[2885] = v[1536] * v[2103] + v[2872] * v[316] + (v[1267] * v[2094] + v[2102] * v[302]) / v[732];
	v[2805] = v[2805] - v[2881];
	v[2809] = v[2809] + v[2881];
	v[2797] = v[2797] - v[2883];
	v[2801] = v[2801] + v[2883];
	v[2789] = v[2789] - v[2885];
	v[2793] = v[2793] + v[2885];
	v[2886] = v[2845] * v[288] + v[1399] * v[4689];
	v[2887] = v[282] * v[2845] + v[1399] * v[4690];
	v[2888] = v[1554] * v[2845] + v[1399] * (v[2022] / v[279] + v[1991] * v[4122]) + v[2886] * v[4871];
	v[2889] = v[1560] * v[2845] + v[1399] * (v[2014] / v[279] + v[1991] * v[4124]) + v[2887] * v[4872];
	v[2890] = (v[218] * v[2886] + v[216] * v[2887] + v[1399] * (v[2030] + v[4126] * v[4598]) + v[2845] * v[4873]) / v[279];
	v[2891] = v[2846] * v[290] + v[1401] * v[4691];
	v[2892] = v[282] * v[2846] + v[1401] * v[4690];
	v[2895] = v[1549] * v[2846] + v[1401] * (v[2028] / v[279] + v[1991] * v[4130]) + v[2891] * v[4875];
	v[2896] = v[2886] + v[2891];
	v[2897] = v[2888] + v[2895];
	v[2899] = (v[1552] * v[1976] + v[207] * v[2892] + v[209] * v[2896] + v[3152] * v[4598] + v[1401] * (v[2010]
		+ v[4134] * v[4598]) + v[1559] * v[4692]) / v[279];
	v[2902] = (v[214] * v[2891] + v[211] * v[2892] + v[1401] * (v[2019] + v[4136] * v[4598]) + v[1553] * v[4692]) / v[279];
	v[2903] = v[2847] * v[288] + v[1403] * v[4689];
	v[2904] = v[2847] * v[290] + v[1403] * v[4691];
	v[2905] = v[2892] + v[2903];
	v[2908] = (v[208] * v[2903] + v[209] * v[2904] + v[1403] * (v[2011] + v[4141] * v[4598]) + v[1557] * v[4693]) / v[279];
	v[2909] = v[2887] + v[2904];
	v[2911] = (v[1562] * v[1970] + v[213] * v[2903] + v[214] * v[2909] + v[3151] * v[4598] + v[1403] * (v[2018]
		+ v[4144] * v[4598]) + v[1564] * v[4693]) / v[279];
	v[2912] = v[2899] + v[2911];
	v[2914] = (v[1563] * v[1966] + v[219] * v[2904] + v[218] * v[2905] + v[3150] * v[4598] + v[1403] * (v[2026]
		+ v[4147] * v[4598]) + v[1567] * v[4693]) / v[279];
	v[2915] = v[2889] + v[2914];
	v[2916] = v[262] * v[2848] + v[1405] * v[4694];
	v[2917] = v[256] * v[2848] + v[1405] * v[4695];
	v[2918] = v[1578] * v[2848] + v[1405] * (v[1946] / v[253] + v[1915] * v[4152]) + v[2916] * v[4880];
	v[2919] = v[1584] * v[2848] + v[1405] * (v[1938] / v[253] + v[1915] * v[4154]) + v[2917] * v[4881];
	v[2920] = (v[162] * v[2916] + v[160] * v[2917] + v[1405] * (v[1954] + v[4156] * v[4597]) + v[2848] * v[4882]) / v[253];
	v[2921] = v[264] * v[2849] + v[1407] * v[4696];
	v[2922] = v[256] * v[2849] + v[1407] * v[4695];
	v[2925] = v[1573] * v[2849] + v[1407] * (v[1952] / v[253] + v[1915] * v[4160]) + v[2921] * v[4884];
	v[2926] = v[2916] + v[2921];
	v[2927] = v[2918] + v[2925];
	v[2929] = (v[1576] * v[1900] + v[151] * v[2922] + v[153] * v[2926] + v[3181] * v[4597] + v[1407] * (v[1934]
		+ v[4164] * v[4597]) + v[1583] * v[4697]) / v[253];
	v[2932] = (v[158] * v[2921] + v[155] * v[2922] + v[1407] * (v[1943] + v[4166] * v[4597]) + v[1577] * v[4697]) / v[253];
	v[2933] = v[262] * v[2850] + v[1409] * v[4694];
	v[2934] = v[264] * v[2850] + v[1409] * v[4696];
	v[2935] = v[2922] + v[2933];
	v[2938] = (v[152] * v[2933] + v[153] * v[2934] + v[1409] * (v[1935] + v[4171] * v[4597]) + v[1581] * v[4698]) / v[253];
	v[2939] = v[2917] + v[2934];
	v[2941] = (v[1586] * v[1894] + v[157] * v[2933] + v[158] * v[2939] + v[3180] * v[4597] + v[1409] * (v[1942]
		+ v[4174] * v[4597]) + v[1588] * v[4698]) / v[253];
	v[2942] = v[2929] + v[2941];
	v[2944] = (v[1587] * v[1890] + v[163] * v[2934] + v[162] * v[2935] + v[3179] * v[4597] + v[1409] * (v[1950]
		+ v[4177] * v[4597]) + v[1591] * v[4698]) / v[253];
	v[2945] = v[2919] + v[2944];
	v[2806] = v[2806] - v[200] * v[2881] + v[2841] * v[353] + v[2842] * v[363];
	v[2807] = v[2807] - v[199] * v[2881] + v[2841] * v[352] + v[2842] * v[361];
	v[2808] = v[2808] - v[198] * v[2881] + v[2841] * v[351] + v[2842] * v[359];
	v[2798] = v[2798] - v[200] * v[2883] + v[2837] * v[353] + v[2838] * v[363];
	v[2799] = v[2799] - v[199] * v[2883] + v[2837] * v[352] + v[2838] * v[361];
	v[2800] = v[2800] - v[198] * v[2883] + v[2837] * v[351] + v[2838] * v[359];
	v[2790] = v[2790] - v[200] * v[2885] + v[2833] * v[353] + v[2834] * v[363];
	v[2791] = v[2791] - v[199] * v[2885] + v[2833] * v[352] + v[2834] * v[361];
	v[2792] = v[2792] - v[198] * v[2885] + v[2833] * v[351] + v[2834] * v[359];
	v[2988] = v[100] * v[2806] + v[1792] * v[2891] + v[292] * v[2961] + v[2904] * v[4457] + v[290] * (v[4763]
		+ v[1991] * v[5025]) + v[2808] * v[98] + v[2807] * v[99];
	v[2989] = (v[1563] * v[2065] + v[2082] * v[288] + v[277] * v[2905] + v[2886] * v[292]) / v[279] + v[289] * v[2965]
		+ v[1991] * v[3147] + v[2950] * v[4699] + v[2951] * v[4699] + v[2808] * v[95] + v[2807] * v[96] + v[2806] * v[97];
	v[2990] = v[277] * v[2968] + v[2887] * v[3207] + v[282] * (v[1991] * v[4700] + v[4754]) + v[2808] * v[92]
		+ v[2807] * v[93] + v[2806] * v[94];
	v[2991] = v[100] * v[2798] + (v[1562] * v[2069] + v[283] * v[2891] + v[2090] * v[290] + v[278] * v[2909]) / v[279]
		+ v[1991] * v[3148] + v[2961] * v[4458] + (v[2963] + v[2964])*v[4701] + v[2800] * v[98] + v[2799] * v[99];
	v[2992] = v[1794] * v[2886] + v[283] * v[2965] + v[2903] * v[4454] + v[288] * (v[4760] + v[1991] * v[5026])
		+ v[2800] * v[95] + v[2799] * v[96] + v[2798] * v[97];
	v[2993] = v[278] * v[2968] + v[2892] * v[3205] + v[282] * (v[1991] * v[4702] + v[4753]) + v[2800] * v[92]
		+ v[2799] * v[93] + v[2798] * v[94];
	v[2994] = v[100] * v[2790] + v[1991] * v[3149] + (v[1552] * v[2070] + v[2086] * v[290] + v[280] * v[2904]
		+ v[2896] * v[4453]) / v[279] + v[2961] * v[4456] + (v[2977] + v[2979])*v[4701] + v[2792] * v[98] + v[2791] * v[99];
	v[2995] = v[2903] * v[3202] + v[2965] * v[4453] + v[288] * (v[1991] * v[4703] + v[4759]) + v[2792] * v[95]
		+ v[2791] * v[96] + v[2790] * v[97];
	v[2996] = v[1795] * v[2887] + v[1793] * v[2892] + v[280] * v[2968] + v[282] * (v[4755] + v[1991] * v[5027])
		+ v[2792] * v[92] + v[2791] * v[93] + v[2790] * v[94];
	v[2997] = -(v[2852] * v[875]);
	v[2998] = -(v[2852] * v[873]);
	v[2999] = -(v[2852] * v[871]);
	v[3000] = -(v[2851] * v[875]);
	v[3001] = -(v[2851] * v[871]);
	v[3002] = -(v[2851] * v[873]);
	v[3003] = -(v[2853] * v[873]);
	v[3004] = -(v[2853] * v[871]);
	v[3005] = -(v[2853] * v[875]);
	v[3006] = -(v[2860] * v[875]);
	v[3007] = -(v[2860] * v[873]);
	v[3008] = -(v[2860] * v[871]);
	v[3009] = -(v[2862] * v[871]);
	v[3010] = -(v[2862] * v[875]);
	v[3011] = -(v[2862] * v[873]);
	v[3012] = -(v[2868] * v[873]);
	v[3013] = -(v[2868] * v[875]);
	v[3014] = -(v[2868] * v[871]);
	v[3015] = -(v[1556] * v[1875]) - v[1566] * v[1877] + 2e0*v[1555] * v[1878] + v[3003] + v[3006] + v[3009] + v[3012]
		+ v[2902] * v[4596] - v[2912] * v[489] - v[2897] * v[492];
	v[3016] = -(v[2861] * v[875]);
	v[3017] = -(v[2861] * v[871]);
	v[3018] = -(v[2861] * v[873]);
	v[3019] = v[1569] * v[1875] + 2e0*v[1561] * v[1877] - v[1566] * v[1878] - v[2998] - v[3001] - v[3013] - v[3016]
		+ v[2908] * v[4595] - v[2912] * v[486] + v[2915] * v[492];
	v[3020] = -(v[1638] * v[1866]) + v[1634] * v[1869] - v[1639] * v[1871] - v[1631] * v[1876] + v[204] * v[2995]
		- v[3003] * v[479] + v[2998] * v[480] - v[3002] * v[481];
	v[3021] = -(v[2866] * v[875]);
	v[3022] = -(v[2866] * v[873]);
	v[3023] = -(v[2866] * v[871]);
	v[3024] = -(v[2867] * v[873]);
	v[3025] = -(v[2867] * v[875]);
	v[3026] = -(v[2867] * v[871]);
	v[3027] = 2e0*v[1548] * v[1875] + v[1569] * v[1877] - v[1556] * v[1878] - v[3004] - v[3017] - v[3021] - v[3024]
		+ v[2890] * v[4594] - v[2897] * v[486] + v[2915] * v[489];
	v[3028] = -(v[1637] * v[1866]) + v[1635] * v[1869] - v[1640] * v[1871] - v[1630] * v[1876] + v[204] * v[2994]
		- v[3004] * v[479] + v[2999] * v[480] - v[3001] * v[481];
	v[3029] = -(v[1646] * v[1866]) + v[1652] * v[1869] - v[1642] * v[1871] - v[1620] * v[1876] + v[204] * v[2993]
		- v[3006] * v[479] + v[3016] * v[480] - v[3010] * v[481];
	v[3030] = v[3000] + v[3011];
	v[3031] = -(v[1645] * v[1866]) + v[1653] * v[1869] - v[1644] * v[1871] - v[1618] * v[1876] + v[204] * v[2991]
		- v[3008] * v[479] + v[3017] * v[480] - v[3009] * v[481];
	v[3032] = -(v[1649] * v[1866]) + v[1661] * v[1869] - v[1657] * v[1871] - v[1608] * v[1876] + v[204] * v[2990]
		- v[3021] * v[479] + v[3025] * v[480] - v[3013] * v[481];
	v[3033] = -(v[1648] * v[1866]) + v[1660] * v[1869] - v[1658] * v[1871] - v[1607] * v[1876] + v[204] * v[2989]
		- v[3022] * v[479] + v[3024] * v[480] - v[3012] * v[481];
	v[3034] = v[3007] + v[3023];
	v[3035] = v[2997] + v[3026];
	v[3036] = -(v[100] * v[2006]) - v[1545] * v[2017] - v[1544] * v[2025] - v[1543] * v[2033] - v[231] * v[2881]
		- v[228] * v[2883] - v[225] * v[2885] + v[2853] * v[530] + v[2852] * v[531] + v[2851] * v[532] + v[2860] * v[545]
		+ v[2861] * v[546] + v[2862] * v[547] + v[2866] * v[560] + v[2867] * v[561] + v[2868] * v[562] - v[2009] * v[94]
		- v[2007] * v[97];
	v[3037] = -(v[1545] * v[2016]) - v[1544] * v[2024] - v[1543] * v[2032] - v[230] * v[2881] - v[227] * v[2883]
		- v[224] * v[2885] + v[2853] * v[527] + v[2852] * v[528] + v[2851] * v[529] + v[2860] * v[542] + v[2861] * v[543]
		+ v[2862] * v[544] + v[2866] * v[557] + v[2867] * v[558] + v[2868] * v[559] - v[2009] * v[93] - v[2007] * v[96]
		- v[2006] * v[99];
	v[3038] = -(v[1545] * v[2015]) - v[1544] * v[2023] - v[1543] * v[2031] - v[229] * v[2881] - v[226] * v[2883]
		- v[223] * v[2885] + v[2853] * v[524] + v[2852] * v[525] + v[2851] * v[526] + v[2860] * v[539] + v[2861] * v[540]
		+ v[2862] * v[541] + v[2866] * v[554] + v[2867] * v[555] + v[2868] * v[556] - v[2009] * v[92] - v[2007] * v[95]
		- v[2006] * v[98];
	v[3039] = v[3038] * v[80] + v[3037] * v[81] + v[3036] * v[82];
	v[2810] = v[2810] + v[144] * v[2881] + v[2843] * v[336] + v[2844] * v[346];
	v[2811] = v[2811] + v[143] * v[2881] + v[2843] * v[335] + v[2844] * v[344];
	v[2812] = v[2812] + v[142] * v[2881] + v[2843] * v[334] + v[2844] * v[342];
	v[2802] = v[2802] + v[144] * v[2883] + v[2839] * v[336] + v[2840] * v[346];
	v[2803] = v[2803] + v[143] * v[2883] + v[2839] * v[335] + v[2840] * v[344];
	v[2804] = v[2804] + v[142] * v[2883] + v[2839] * v[334] + v[2840] * v[342];
	v[2794] = v[2794] + v[144] * v[2885] + v[2835] * v[336] + v[2836] * v[346];
	v[2795] = v[2795] + v[143] * v[2885] + v[2835] * v[335] + v[2836] * v[344];
	v[2796] = v[2796] + v[142] * v[2885] + v[2835] * v[334] + v[2836] * v[342];
	v[3085] = v[1784] * v[2921] + v[266] * v[3058] + v[2934] * v[4449] + v[264] * (v[4749] + v[1915] * v[5028])
		+ v[2812] * v[55] + v[2811] * v[56] + v[2810] * v[57];
	v[3086] = (v[1587] * v[2037] + v[2054] * v[262] + v[266] * v[2916] + v[251] * v[2935]) / v[253] + v[263] * v[3062]
		+ v[1915] * v[3176] + v[3047] * v[4707] + v[3048] * v[4707] + v[2812] * v[52] + v[2811] * v[53] + v[2810] * v[54];
	v[3087] = v[251] * v[3065] + v[2917] * v[3200] + v[256] * (v[1915] * v[4708] + v[4740]) + v[2812] * v[49]
		+ v[2811] * v[50] + v[2810] * v[51];
	v[3088] = (v[1586] * v[2041] + v[2062] * v[264] + v[257] * v[2921] + v[252] * v[2939]) / v[253] + v[1915] * v[3177]
		+ v[3058] * v[4450] + (v[3060] + v[3061])*v[4709] + v[2804] * v[55] + v[2803] * v[56] + v[2802] * v[57];
	v[3089] = v[1786] * v[2916] + v[257] * v[3062] + v[2933] * v[4446] + v[262] * (v[4746] + v[1915] * v[5029])
		+ v[2804] * v[52] + v[2803] * v[53] + v[2802] * v[54];
	v[3090] = v[252] * v[3065] + v[2922] * v[3198] + v[256] * (v[1915] * v[4710] + v[4739]) + v[2804] * v[49]
		+ v[2803] * v[50] + v[2802] * v[51];
	v[3091] = v[1915] * v[3178] + (v[1576] * v[2042] + v[2058] * v[264] + v[254] * v[2934] + v[2926] * v[4445]) / v[253]
		+ v[3058] * v[4448] + (v[3074] + v[3076])*v[4709] + v[2796] * v[55] + v[2795] * v[56] + v[2794] * v[57];
	v[3092] = v[2933] * v[3195] + v[3062] * v[4445] + v[262] * (v[1915] * v[4711] + v[4745]) + v[2796] * v[52]
		+ v[2795] * v[53] + v[2794] * v[54];
	v[3093] = v[1787] * v[2917] + v[1785] * v[2922] + v[254] * v[3065] + v[2796] * v[49] + v[2795] * v[50] + v[256] * (v[4741]
		+ v[1915] * v[5030]) + v[2794] * v[51];
	v[3094] = v[2855] * v[881];
	v[3095] = v[2855] * v[879];
	v[3096] = v[2855] * v[877];
	v[3097] = v[2854] * v[881];
	v[3098] = v[2854] * v[877];
	v[3099] = v[2854] * v[879];
	v[3100] = v[2856] * v[879];
	v[3101] = v[2856] * v[877];
	v[3102] = v[2856] * v[881];
	v[3103] = v[2857] * v[881];
	v[3104] = v[2857] * v[879];
	v[3105] = v[2857] * v[877];
	v[3106] = v[2859] * v[877];
	v[3107] = v[2859] * v[881];
	v[3108] = v[2859] * v[879];
	v[3109] = v[2865] * v[879];
	v[3110] = v[2865] * v[881];
	v[3111] = v[2865] * v[877];
	v[3112] = -(v[1580] * v[1849]) - v[1590] * v[1851] + 2e0*v[1579] * v[1852] + v[3100] + v[3103] + v[3106] + v[3109]
		- v[2942] * v[390] - v[2927] * v[393] + v[2932] * v[4592];
	v[3113] = v[2858] * v[881];
	v[3114] = v[2858] * v[877];
	v[3115] = v[2858] * v[879];
	v[3116] = v[1593] * v[1849] + 2e0*v[1585] * v[1851] - v[1590] * v[1852] - v[3095] - v[3098] - v[3110] - v[3113]
		- v[2942] * v[387] + v[2945] * v[393] + v[2938] * v[4591];
	v[3117] = -(v[1722] * v[1840]) + v[1718] * v[1843] - v[1723] * v[1845] - v[1715] * v[1850] + v[148] * v[3092]
		- v[3100] * v[380] + v[3095] * v[381] - v[3099] * v[382];
	v[3118] = v[2863] * v[881];
	v[3119] = v[2863] * v[879];
	v[3120] = v[2863] * v[877];
	v[3121] = v[2864] * v[879];
	v[3122] = v[2864] * v[881];
	v[3123] = v[2864] * v[877];
	v[3124] = 2e0*v[1572] * v[1849] + v[1593] * v[1851] - v[1580] * v[1852] - v[3101] - v[3114] - v[3118] - v[3121]
		- v[2927] * v[387] + v[2945] * v[390] + v[2920] * v[4590];
	v[3125] = -(v[1721] * v[1840]) + v[1719] * v[1843] - v[1724] * v[1845] - v[1714] * v[1850] + v[148] * v[3091]
		- v[3101] * v[380] + v[3096] * v[381] - v[3098] * v[382];
	v[3126] = -(v[1730] * v[1840]) + v[1736] * v[1843] - v[1726] * v[1845] - v[1704] * v[1850] + v[148] * v[3090]
		- v[3103] * v[380] + v[3113] * v[381] - v[3107] * v[382];
	v[3127] = v[3097] + v[3108];
	v[3128] = -(v[1729] * v[1840]) + v[1737] * v[1843] - v[1728] * v[1845] - v[1702] * v[1850] + v[148] * v[3088]
		- v[3105] * v[380] + v[3114] * v[381] - v[3106] * v[382];
	v[3129] = -(v[1733] * v[1840]) + v[1745] * v[1843] - v[1741] * v[1845] - v[1692] * v[1850] + v[148] * v[3087]
		- v[3118] * v[380] + v[3122] * v[381] - v[3110] * v[382];
	v[3130] = -(v[1732] * v[1840]) + v[1744] * v[1843] - v[1742] * v[1845] - v[1691] * v[1850] + v[148] * v[3086]
		- v[3119] * v[380] + v[3121] * v[381] - v[3109] * v[382];
	v[3131] = v[3104] + v[3120];
	v[3132] = v[3094] + v[3123];
	v[3133] = v[1545] * v[1941] + v[1544] * v[1949] + v[1543] * v[1957] + v[175] * v[2881] + v[172] * v[2883]
		+ v[169] * v[2885] + v[2856] * v[431] + v[2855] * v[432] + v[2854] * v[433] + v[2857] * v[446] + v[2858] * v[447]
		+ v[2859] * v[448] + v[2863] * v[461] + v[2864] * v[462] + v[2865] * v[463] + v[1933] * v[51] + v[1931] * v[54]
		+ v[1930] * v[57];
	v[3134] = v[1545] * v[1940] + v[1544] * v[1948] + v[1543] * v[1956] + v[174] * v[2881] + v[171] * v[2883]
		+ v[168] * v[2885] + v[2856] * v[428] + v[2855] * v[429] + v[2854] * v[430] + v[2857] * v[443] + v[2858] * v[444]
		+ v[2859] * v[445] + v[2863] * v[458] + v[2864] * v[459] + v[2865] * v[460] + v[1933] * v[50] + v[1931] * v[53]
		+ v[1930] * v[56];
	v[3135] = v[1545] * v[1939] + v[1544] * v[1947] + v[1543] * v[1955] + v[173] * v[2881] + v[170] * v[2883]
		+ v[167] * v[2885] + v[2856] * v[425] + v[2855] * v[426] + v[2854] * v[427] + v[2857] * v[440] + v[2858] * v[441]
		+ v[2859] * v[442] + v[2863] * v[455] + v[2864] * v[456] + v[2865] * v[457] + v[1933] * v[49] + v[1931] * v[52]
		+ v[1930] * v[55];
	v[3136] = v[3135] * v[37] + v[3134] * v[38] + v[3133] * v[39];
	v[3137] = (-v[3039] + v[3038] * v[86] + v[3037] * v[87] + v[3036] * v[88]) / 2e0;
	v[3138] = (-v[3039] + v[3038] * v[83] + v[3037] * v[84] + v[3036] * v[85]) / 2e0;
	v[3139] = (-(v[1659] * v[1817]) - v[1641] * v[1820] - v[1643] * v[1867] - 2e0*v[2888] + 2e0*v[2895]
		+ v[3004] * v[4561] + v[3003] * v[4562] + v[3008] * v[4563] + v[3006] * v[4564] + v[3022] * v[4565] + v[3021] * v[4566]
		- v[1642] * v[4712] - v[1639] * v[4713] - v[1658] * v[4714] - v[1644] * v[4715] - v[1657] * v[4716] - v[1640] * v[4717]
		- v[3005] * v[482] - v[3007] * v[501] - v[3023] * v[519]) / 2e0;
	v[3140] = (-(v[1636] * v[1866]) + v[1633] * v[1869] - v[1641] * v[1871] - v[1632] * v[1876] + v[204] * v[2996]
		- v[3005] * v[479] + v[2997] * v[480] - v[3000] * v[481]) / 2e0;
	v[3142] = (v[1662] * v[1817] + v[1633] * v[1820] + v[1654] * v[1867] - 2e0*v[2889] + 2e0*v[2914] - v[2999] * v[4561]
		- v[2998] * v[4562] - v[3017] * v[4563] - v[3016] * v[4564] - v[3024] * v[4565] - v[3025] * v[4566] + v[1652] * v[4712]
		+ v[1634] * v[4713] + v[1660] * v[4714] + v[1653] * v[4715] + v[1661] * v[4716] + v[1635] * v[4717] + v[2997] * v[482]
		+ v[3018] * v[501] + v[3026] * v[519]) / 2e0;
	v[3143] = (-(v[1647] * v[1866]) + v[1654] * v[1869] - v[1643] * v[1871] - v[1619] * v[1876] + v[204] * v[2992]
		- v[3007] * v[479] + v[3018] * v[480] - v[3011] * v[481]) / 2e0;
	v[3144] = (-(v[1650] * v[1817]) - v[1636] * v[1820] - v[1647] * v[1867] - 2e0*v[2899] + 2e0*v[2911]
		+ v[3001] * v[4561] + v[3002] * v[4562] + v[3009] * v[4563] + v[3010] * v[4564] + v[3012] * v[4565] + v[3013] * v[4566]
		- v[1646] * v[4712] - v[1638] * v[4713] - v[1648] * v[4714] - v[1645] * v[4715] - v[1649] * v[4716] - v[1637] * v[4717]
		- v[3000] * v[482] - v[3011] * v[501] - v[3014] * v[519]) / 2e0;
	v[3204] = v[1868] * v[3142] - v[3143] + v[1865] * v[3144] + v[1819] * v[4905] * v[5032] - v[4428] * ((v[1606] * v[1817])
		/ 2e0 + (v[1632] * v[1820]) / 2e0 + v[1620] * v[1859] + v[1631] * v[1860] + v[1607] * v[1861] + v[1618] * v[1862]
		+ v[1608] * v[1863] + v[1630] * v[1864] + (v[1619] * v[1867]) / 2e0 + v[2999] - v[3002] - v[3008] + v[3010] + v[3022]
		- v[3025] - 2e0*v[1876] * v[3146] + v[3015] * v[375] - v[3019] * v[376] - v[3027] * v[378] + v[2988] * v[4537]
		+ v[2992] * v[4538] + v[2996] * v[4539] + v[3030] * v[4718] + v[3034] * v[4719] + v[3035] * v[4720] + v[2995] * v[487]
		+ v[2994] * v[493] + v[2993] * v[497] + v[4586] * (v[2890] + v[2082] * v[2893] + v[2079] * v[2894] + v[2073] * v[2896]
			+ v[2088] * v[2898] + v[2092] * v[2900] + v[2090] * v[2901] + v[2902] + v[2067] * v[2905] + v[2075] * v[2906]
			+ v[2086] * v[2907] + v[2908] + v[2071] * v[2909] + v[2078] * v[2910] + v[2084] * v[2913] + v[2030] * v[2950]
			+ v[2028] * v[2951] + v[2026] * v[2953] + v[2022] * v[2963] + v[2019] * v[2964] + v[2018] * v[2967] + v[2011] * v[2977]
			+ v[2014] * v[2979] + v[2010] * v[2980] + v[1966] * v[3147] + v[1970] * v[3148] + v[1976] * v[3149] + v[2065] * v[3150]
			+ v[2069] * v[3151] + v[2070] * v[3152] + v[2886] * v[3156] + v[2887] * v[3157] + v[2891] * v[3158] + v[2892] * v[3159]
			+ v[2903] * v[3160] + v[2904] * v[3161] + v[2845] * v[4721] + v[2846] * v[4722] + v[2847] * v[4723]
			+ v[1991] * v[4906] * v[5034]) + v[2991] * v[506] + v[2990] * v[510] + v[2989] * v[514] + v[7393 + i1194]) + v[4903] * (
				-4e0*v[1676] * v[1819] + v[3139] * v[4904] + v[7417 + i1194]);
	v[4757] = v[3204] + (v[1650] * v[1866] - v[1662] * v[1869] + v[1659] * v[1871] + v[1606] * v[1876] - v[204] * v[2988]
		+ v[3023] * v[479] - v[3026] * v[480] + v[3014] * v[481]) / 2e0;
	v[3163] = v[3028] + v[3032];
	v[3164] = v[3031] + v[3033];
	v[3165] = v[3020] + v[3029];
	v[3166] = (-v[3136] + v[3135] * v[43] + v[3134] * v[44] + v[3133] * v[45]) / 2e0;
	v[3167] = (-v[3136] + v[3135] * v[40] + v[3134] * v[41] + v[3133] * v[42]) / 2e0;
	v[3168] = (-(v[1743] * v[1807]) - v[1725] * v[1810] - v[1727] * v[1841] - 2e0*v[2918] + 2e0*v[2925] - v[3102] * v[383]
		- v[3104] * v[402] - v[3120] * v[420] + v[3101] * v[4567] + v[3100] * v[4568] + v[3105] * v[4569] + v[3103] * v[4570]
		+ v[3119] * v[4571] + v[3118] * v[4572] - v[1726] * v[4724] - v[1723] * v[4725] - v[1742] * v[4726] - v[1728] * v[4727]
		- v[1741] * v[4728] - v[1724] * v[4729]) / 2e0;
	v[3169] = (-(v[1720] * v[1840]) + v[1717] * v[1843] - v[1725] * v[1845] - v[1716] * v[1850] + v[148] * v[3093]
		- v[3102] * v[380] + v[3094] * v[381] - v[3097] * v[382]) / 2e0;
	v[3171] = (v[1746] * v[1807] + v[1717] * v[1810] + v[1738] * v[1841] - 2e0*v[2919] + 2e0*v[2944] + v[3094] * v[383]
		+ v[3115] * v[402] + v[3123] * v[420] - v[3096] * v[4567] - v[3095] * v[4568] - v[3114] * v[4569] - v[3113] * v[4570]
		- v[3121] * v[4571] - v[3122] * v[4572] + v[1736] * v[4724] + v[1718] * v[4725] + v[1744] * v[4726] + v[1737] * v[4727]
		+ v[1745] * v[4728] + v[1719] * v[4729]) / 2e0;
	v[3172] = (-(v[1731] * v[1840]) + v[1738] * v[1843] - v[1727] * v[1845] - v[1703] * v[1850] + v[148] * v[3089]
		- v[3104] * v[380] + v[3115] * v[381] - v[3108] * v[382]) / 2e0;
	v[3173] = (-(v[1734] * v[1807]) - v[1720] * v[1810] - v[1731] * v[1841] - 2e0*v[2929] + 2e0*v[2941] - v[3097] * v[383]
		- v[3108] * v[402] - v[3111] * v[420] + v[3098] * v[4567] + v[3099] * v[4568] + v[3106] * v[4569] + v[3107] * v[4570]
		+ v[3109] * v[4571] + v[3110] * v[4572] - v[1730] * v[4724] - v[1722] * v[4725] - v[1732] * v[4726] - v[1729] * v[4727]
		- v[1733] * v[4728] - v[1721] * v[4729]) / 2e0;
	v[3197] = v[1842] * v[3171] - v[3172] + v[1839] * v[3173] + v[1809] * v[4915] * v[5036] - v[4426] * ((v[1690] * v[1807])
		/ 2e0 + (v[1716] * v[1810]) / 2e0 + v[1704] * v[1833] + v[1715] * v[1834] + v[1691] * v[1835] + v[1702] * v[1836]
		+ v[1692] * v[1837] + v[1714] * v[1838] + (v[1703] * v[1841]) / 2e0 + v[3096] - v[3099] - v[3105] + v[3107] + v[3119]
		- v[3122] - 2e0*v[1850] * v[3175] + v[3112] * v[369] - v[3116] * v[370] - v[3124] * v[372] + v[3092] * v[388]
		+ v[3091] * v[394] + v[3090] * v[398] + v[3088] * v[407] + v[3087] * v[411] + v[3086] * v[415] + v[3085] * v[4558]
		+ v[3089] * v[4559] + v[3093] * v[4560] + v[3127] * v[4730] + v[3131] * v[4731] + v[3132] * v[4732] + v[4585] * (v[2920]
			+ v[2054] * v[2923] + v[2051] * v[2924] + v[2045] * v[2926] + v[2060] * v[2928] + v[2064] * v[2930] + v[2062] * v[2931]
			+ v[2932] + v[2039] * v[2935] + v[2047] * v[2936] + v[2058] * v[2937] + v[2938] + v[2043] * v[2939] + v[2050] * v[2940]
			+ v[2056] * v[2943] + v[1954] * v[3047] + v[1952] * v[3048] + v[1950] * v[3050] + v[1946] * v[3060] + v[1943] * v[3061]
			+ v[1942] * v[3064] + v[1935] * v[3074] + v[1938] * v[3076] + v[1934] * v[3077] + v[1890] * v[3176] + v[1894] * v[3177]
			+ v[1900] * v[3178] + v[2037] * v[3179] + v[2041] * v[3180] + v[2042] * v[3181] + v[2916] * v[3185] + v[2917] * v[3186]
			+ v[2921] * v[3187] + v[2922] * v[3188] + v[2933] * v[3189] + v[2934] * v[3190] + v[2848] * v[4733] + v[2849] * v[4734]
			+ v[2850] * v[4735] + v[1915] * v[4916] * v[5038]) + v[7405 + i1194]) + v[4913] * (-4e0*v[1760] * v[1809]
				+ v[3168] * v[4914] + v[7501 + i1194]);
	v[4743] = v[3197] + (v[1734] * v[1840] - v[1746] * v[1843] + v[1743] * v[1845] + v[1690] * v[1850] - v[148] * v[3085]
		+ v[3120] * v[380] - v[3123] * v[381] + v[3111] * v[382]) / 2e0;
	v[3192] = v[3125] + v[3129];
	v[3193] = v[3128] + v[3130];
	v[3194] = v[3117] + v[3126];
	v[7574] = 0e0;
	v[7575] = 0e0;
	v[7576] = 0e0;
	v[7577] = 0e0;
	v[7578] = v[1779];
	v[7579] = v[1777];
	v[7580] = 0e0;
	v[7581] = 0e0;
	v[7582] = 0e0;
	v[7583] = 0e0;
	v[7584] = 0e0;
	v[7585] = 0e0;
	v[7538] = 0e0;
	v[7539] = 0e0;
	v[7540] = 0e0;
	v[7541] = v[1779];
	v[7542] = 0e0;
	v[7543] = v[1778];
	v[7544] = 0e0;
	v[7545] = 0e0;
	v[7546] = 0e0;
	v[7547] = 0e0;
	v[7548] = 0e0;
	v[7549] = 0e0;
	v[7514] = 0e0;
	v[7515] = 0e0;
	v[7516] = 0e0;
	v[7517] = v[1777];
	v[7518] = v[1778];
	v[7519] = 0e0;
	v[7520] = 0e0;
	v[7521] = 0e0;
	v[7522] = 0e0;
	v[7523] = 0e0;
	v[7524] = 0e0;
	v[7525] = 0e0;
	v[7490] = 0e0;
	v[7491] = 0e0;
	v[7492] = 0e0;
	v[7493] = 0e0;
	v[7494] = 0e0;
	v[7495] = 0e0;
	v[7496] = 0e0;
	v[7497] = 0e0;
	v[7498] = 0e0;
	v[7499] = 0e0;
	v[7500] = v[1770];
	v[7501] = v[1768];
	v[7454] = 0e0;
	v[7455] = 0e0;
	v[7456] = 0e0;
	v[7457] = 0e0;
	v[7458] = 0e0;
	v[7459] = 0e0;
	v[7460] = 0e0;
	v[7461] = 0e0;
	v[7462] = 0e0;
	v[7463] = v[1770];
	v[7464] = 0e0;
	v[7465] = v[1769];
	v[7430] = 0e0;
	v[7431] = 0e0;
	v[7432] = 0e0;
	v[7433] = 0e0;
	v[7434] = 0e0;
	v[7435] = 0e0;
	v[7436] = 0e0;
	v[7437] = 0e0;
	v[7438] = 0e0;
	v[7439] = v[1768];
	v[7440] = v[1769];
	v[7441] = 0e0;
	v[7586] = v[2793] - v[4706];
	v[7587] = v[2801] - v[4705];
	v[7588] = v[2809] - v[4704];
	v[7589] = -v[3128] + v[3130] + v[3194] * v[369] + v[3192] * v[372] + v[10] * (v[1407] * v[2040] + v[1405] * v[2044]
		+ v[1915] * (v[3074] + v[1405] * v[3186] + v[1407] * v[3188]) + v[2850] * v[3195] + v[2849] * v[4446] + v[2848] * v[4449]
		+ v[1857] * v[4736] + v[1892] * v[4737] + v[1898] * v[4738] + v[155] * v[4739] + v[160] * v[4740] + v[151] * v[4741])
		+ 2e0*(v[1832] * v[3196] + v[3168] * v[4426] + v[3131] * v[4430] + v[1809] * (-(v[1771] * v[4425]) + v[1757] * v[4742])
			) + v[368] * v[4743] + (v[1739] * v[1850] - v[148] * v[3116] + v[7573 + i1194]) / 2e0;
	v[7590] = v[3125] - v[3129] + v[3194] * v[370] + v[3193] * v[372] + 2e0*(v[1831] * v[3199] - v[3171] * v[4426]
		+ v[3132] * v[4430] + v[1809] * (v[1773] * v[4425] + v[1758] * v[4742])) + v[371] * (-v[3169] + v[3172] + v[4743])
		+ v[10] * (v[1409] * v[2040] + v[1405] * v[2046] + v[1784] * v[2848] + v[1785] * v[2850] + v[1915] * (v[3061]
			+ v[1405] * v[3185] + v[1409] * v[3189]) + v[2849] * v[3198] + v[1890] * v[4540] + v[1909] * v[4542] + v[1855] * v[4744]
			+ v[152] * v[4745] + v[157] * v[4746] + v[2054] * v[4925]) + (-(v[1735] * v[1850]) + v[148] * v[3112] + v[7537 + i1194])
		/ 2e0;
	v[7591] = -v[3117] + v[3126] + v[3193] * v[369] + v[3192] * v[370] + (-v[3169] + v[3197])*v[373] + 2e0*
		(v[1830] * v[3201] + v[3173] * v[4426] + v[3127] * v[4430] + v[1809] * (-(v[1775] * v[4425]) + v[1753] * v[4742]))
		+ v[10] * (v[1409] * v[2044] + v[1407] * v[2046] + v[1786] * v[2849] + v[1787] * v[2850] + v[1915] * (v[3047]
			+ v[1407] * v[3187] + v[1409] * v[3190]) + v[2848] * v[3200] + v[1894] * v[4541] + v[1853] * v[4747] + v[1900] * v[4748]
			+ v[163] * v[4749] + v[2058] * v[4929] + v[2062] * v[4930]) + (v[1750] * v[1850] - v[148] * v[3124] + v[7513 + i1194])
		/ 2e0;
	v[7592] = v[2789] + v[4706];
	v[7593] = v[2797] + v[4705];
	v[7594] = v[2805] + v[4704];
	v[7595] = -v[3031] + v[3033] + v[3165] * v[375] + v[3163] * v[378] + v[10] * (v[1401] * v[2068] + v[1399] * v[2072]
		+ v[1991] * (v[2977] + v[1399] * v[3157] + v[1401] * v[3159]) + v[2847] * v[3202] + v[2846] * v[4454] + v[2845] * v[4457]
		+ v[1883] * v[4750] + v[1968] * v[4751] + v[1974] * v[4752] + v[211] * v[4753] + v[216] * v[4754] + v[207] * v[4755])
		+ 2e0*(v[1829] * v[3203] + v[3139] * v[4428] + v[3034] * v[4437] + v[1819] * (-(v[1762] * v[4427]) + v[1673] * v[4756])
			) + v[374] * v[4757] + (v[1655] * v[1876] - v[204] * v[3019] + v[7489 + i1194]) / 2e0;
	v[7596] = v[3028] - v[3032] + v[3165] * v[376] + v[3164] * v[378] + 2e0*(v[1828] * v[3206] - v[3142] * v[4428]
		+ v[3035] * v[4437] + v[1819] * (v[1764] * v[4427] + v[1674] * v[4756])) + v[377] * (-v[3140] + v[3143] + v[4757])
		+ v[10] * (v[1403] * v[2068] + v[1399] * v[2074] + v[1792] * v[2845] + v[1793] * v[2847] + v[1991] * (v[2964]
			+ v[1399] * v[3156] + v[1403] * v[3160]) + v[2846] * v[3205] + v[1966] * v[4519] + v[1985] * v[4521] + v[1881] * v[4758]
			+ v[208] * v[4759] + v[213] * v[4760] + v[2082] * v[4937]) + (-(v[1651] * v[1876]) + v[204] * v[3015] + v[7453 + i1194])
		/ 2e0;
	v[7597] = -v[3020] + v[3029] + v[3164] * v[375] + v[3163] * v[376] + (-v[3140] + v[3204])*v[379] + 2e0*
		(v[1827] * v[3208] + v[3144] * v[4428] + v[3030] * v[4437] + v[1819] * (-(v[1766] * v[4427]) + v[1669] * v[4756]))
		+ v[10] * (v[1403] * v[2072] + v[1401] * v[2074] + v[1794] * v[2846] + v[1795] * v[2847] + v[1991] * (v[2950]
			+ v[1401] * v[3158] + v[1403] * v[3161]) + v[2845] * v[3207] + v[1970] * v[4520] + v[1879] * v[4761] + v[1976] * v[4762]
			+ v[219] * v[4763] + v[2086] * v[4941] + v[2090] * v[4942]) + (v[1666] * v[1876] - v[204] * v[3027] + v[7429 + i1194])
		/ 2e0;
	Rc[i1194 - 1] += v[6729 + i1194];
	for (i1799 = 1; i1799 <= 12; i1799++) {
		Kc[i1194 - 1][i1799 - 1] += v[3167] * v[6157 + i1799] + v[3166] * v[6169 + i1799] + v[3138] * v[6181 + i1799]
			+ v[3137] * v[6193 + i1799] + v[7585 + i1799];
	};/* end for */
};/* end for */
v[3217] = 0e0;
v[3218] = 0e0;
v[3219] = 0e0;
v[3220] = 0e0;
b3221 = b733;
if (b3221) {
	b3222 = b752;
	b3226 = b735;
	if (b3226) {
		v[3219] = 0e0;
		v[3218] = 0e0;
		v[3217] = 0e0;
		v[3220] = 0e0;
	}
	else {
	};
}
else {
};
v[4768] = -(v[3217] * v[331]);
v[4769] = -(v[3218] * v[332]);
v[4770] = -(v[3219] * v[333]);
v[3227] = 0e0;
v[3228] = 0e0;
v[3229] = 0e0;
v[3230] = 0e0;
b3231 = b733;
if (b3231) {
	b3232 = b735;
	if (b3232) {
		v[3229] = 0e0;
		v[3228] = 0e0;
		v[3227] = 0e0;
		v[3230] = -v[3220];
	}
	else {
	};
}
else {
};
v[3233] = v[301] * v[3219];
v[3927] = -(v[3233] * v[333]);
v[3234] = v[291] * v[3219];
v[3924] = -(v[3234] * v[333]);
v[3235] = v[281] * v[3219];
v[3930] = -(v[3235] * v[333]);
v[3236] = v[275] * v[3219];
v[3936] = v[3236] * v[333];
v[3237] = v[265] * v[3219];
v[3933] = v[3237] * v[333];
v[3238] = v[255] * v[3219];
v[4814] = v[3238] * v[470] + v[3237] * v[471] + v[3236] * v[472] - v[3235] * v[569] - v[3234] * v[570] - v[3233] * v[571];
v[4801] = v[3238] * v[473] + v[3237] * v[474] + v[3236] * v[475] - v[3235] * v[572] - v[3234] * v[573] - v[3233] * v[574];
v[3903] = v[3238] * v[333];
v[3227] = v[3227] + v[3219] * v[3239];
v[3228] = v[3228] + v[3219] * v[3240];
v[3826] = v[3219] * v[3241] + v[3238] * v[476] + v[3237] * v[477] + v[3236] * v[478] - v[3235] * v[575] - v[3234] * v[576]
- v[3233] * v[577];
v[3229] = v[2709] * v[3219] + v[3229];
v[3255] = v[301] * v[3218];
v[3928] = -(v[3255] * v[332]);
v[4818] = -v[3927] - v[3928];
v[3256] = v[291] * v[3218];
v[3925] = -(v[3256] * v[332]);
v[4817] = -v[3924] - v[3925];
v[3257] = v[281] * v[3218];
v[3931] = -(v[3257] * v[332]);
v[4819] = -v[3930] - v[3931];
v[3258] = v[275] * v[3218];
v[3937] = v[3258] * v[332];
v[4816] = v[3936] + v[3937];
v[3259] = v[265] * v[3218];
v[3934] = v[3259] * v[332];
v[4815] = v[3933] + v[3934];
v[3260] = v[255] * v[3218];
v[4813] = v[3260] * v[470] + v[3259] * v[471] + v[3258] * v[472] - v[3257] * v[569] - v[3256] * v[570] - v[3255] * v[571];
v[4794] = v[3260] * v[476] + v[3259] * v[477] + v[3258] * v[478] - v[3257] * v[575] - v[3256] * v[576] - v[3255] * v[577];
v[3904] = v[3260] * v[332];
v[4811] = v[3903] + v[3904];
v[3227] = v[3227] + v[3218] * v[3261];
v[3859] = v[3218] * v[3262] + v[3260] * v[473] + v[3259] * v[474] + v[3258] * v[475] - v[3257] * v[572] - v[3256] * v[573]
- v[3255] * v[574];
v[3228] = v[2734] * v[3218] + v[3228];
v[3229] = v[3229] + v[3218] * v[3264];
v[3277] = v[301] * v[3217];
v[3874] = -(v[3277] * v[331]);
v[5074] = 2e0*v[3874] + v[3927] + v[3928];
v[5062] = v[3874] + v[3927] - v[3255] * v[4607];
v[5049] = v[3874] + v[3928] - v[3233] * v[4606];
v[4827] = -v[3874] - v[3928];
v[4822] = -v[3874] - v[3927];
v[3278] = v[291] * v[3217];
v[3876] = -(v[3278] * v[331]);
v[5073] = 2e0*v[3876] + v[3924] + v[3925];
v[5063] = v[3876] + v[3924] - v[3256] * v[4607];
v[5051] = v[3876] + v[3925] - v[3234] * v[4606];
v[4826] = -v[3876] - v[3925];
v[4821] = -v[3876] - v[3924];
v[3279] = v[281] * v[3217];
v[3872] = -(v[3279] * v[331]);
v[5075] = 2e0*v[3872] + v[3930] + v[3931];
v[5061] = v[3872] + v[3930] - v[3257] * v[4607];
v[5050] = v[3872] + v[3931] - v[3235] * v[4606];
v[4828] = -v[3872] - v[3931];
v[4823] = -v[3872] - v[3930];
v[3280] = v[275] * v[3217];
v[3880] = v[3280] * v[331];
v[5077] = 2e0*v[3880] + v[4816];
v[4803] = v[3880] + v[3936];
v[5065] = 2e0*v[3937] + v[4803];
v[4795] = v[3880] + v[3937];
v[5052] = 2e0*v[3936] + v[4795];
v[3281] = v[265] * v[3217];
v[3882] = v[3281] * v[331];
v[5076] = 2e0*v[3882] + v[4815];
v[4804] = v[3882] + v[3933];
v[5066] = 2e0*v[3934] + v[4804];
v[4797] = v[3882] + v[3934];
v[5054] = 2e0*v[3933] + v[4797];
v[3282] = v[255] * v[3217];
v[4800] = v[3282] * v[473] + v[3281] * v[474] + v[3280] * v[475] - v[3279] * v[572] - v[3278] * v[573] - v[3277] * v[574];
v[4793] = v[3282] * v[476] + v[3281] * v[477] + v[3280] * v[478] - v[3279] * v[575] - v[3278] * v[576] - v[3277] * v[577];
v[3878] = v[3282] * v[331];
v[5078] = 2e0*v[3878] + v[4811];
v[4802] = v[3878] + v[3903];
v[5064] = 2e0*v[3904] + v[4802];
v[4796] = v[3878] + v[3904];
v[5053] = 2e0*v[3903] + v[4796];
v[3905] = v[3217] * v[3283] + v[3282] * v[470] + v[3281] * v[471] + v[3280] * v[472] - v[3279] * v[569] - v[3278] * v[570]
- v[3277] * v[571];
v[3227] = v[2759] * v[3217] + v[3227];
v[3228] = v[3228] + v[3217] * v[3285];
v[3229] = v[3229] + v[3217] * v[3286];
v[3287] = v[1398] * v[3217] + v[1342] * v[3218] + v[1294] * v[3219];
v[4764] = -(v[1675] * v[3287]);
v[4198] = v[4456] * v[4764];
v[4188] = v[4458] * v[4764];
v[4181] = v[292] * v[4764];
v[3356] = v[1794] * v[3287];
v[3352] = v[1795] * v[3287];
v[3342] = v[3207] * v[3287];
v[3288] = v[1400] * v[3217] + v[1344] * v[3218] + v[1296] * v[3219];
v[4199] = v[3288] * v[4528];
v[4189] = v[3154] * v[3288];
v[4182] = v[3288] * v[4765];
v[3360] = v[1792] * v[3288];
v[4935] = v[3342] + v[3360];
v[3346] = v[3205] * v[3288];
v[4939] = v[3346] + v[3356];
v[3289] = v[1402] * v[3217] + v[1346] * v[3218] + v[1298] * v[3219];
v[4196] = v[3155] * v[3289];
v[4891] = v[4196] + v[4199];
v[5094] = v[4198] + v[4891];
v[4192] = v[3289] * v[4529];
v[4890] = v[4189] + v[4192];
v[5093] = v[4188] + v[4890];
v[4184] = v[3289] * v[4531];
v[4888] = v[4181] + v[4184];
v[5092] = v[4182] + v[4888];
v[3350] = v[3202] * v[3289];
v[4940] = v[3350] + v[3352];
v[4775] = v[1793] * v[3288] + v[3350];
v[4931] = v[3352] + v[4775];
v[4774] = v[3346] + v[3289] * v[4454];
v[4934] = v[3356] + v[4774];
v[4773] = v[3342] + v[3289] * v[4457];
v[4938] = v[3360] + v[4773];
v[3290] = v[1404] * v[3217] + v[1348] * v[3218] + v[1300] * v[3219];
v[4766] = -(v[1759] * v[3290]);
v[4289] = v[4448] * v[4766];
v[4279] = v[4450] * v[4766];
v[4272] = v[266] * v[4766];
v[3389] = v[1786] * v[3290];
v[3385] = v[1787] * v[3290];
v[3375] = v[3200] * v[3290];
v[3291] = v[1406] * v[3217] + v[1350] * v[3218] + v[1302] * v[3219];
v[4290] = v[3291] * v[4549];
v[4280] = v[3183] * v[3291];
v[4273] = v[3291] * v[4767];
v[3393] = v[1784] * v[3291];
v[4923] = v[3375] + v[3393];
v[3379] = v[3198] * v[3291];
v[4927] = v[3379] + v[3389];
v[3292] = v[1408] * v[3217] + v[1352] * v[3218] + v[1304] * v[3219];
v[4287] = v[3184] * v[3292];
v[4896] = v[4287] + v[4290];
v[5097] = v[4289] + v[4896];
v[4283] = v[3292] * v[4550];
v[4895] = v[4280] + v[4283];
v[5096] = v[4279] + v[4895];
v[4275] = v[3292] * v[4552];
v[4893] = v[4272] + v[4275];
v[5095] = v[4273] + v[4893];
v[3383] = v[3195] * v[3292];
v[4928] = v[3383] + v[3385];
v[4781] = v[1785] * v[3291] + v[3383];
v[4919] = v[3385] + v[4781];
v[4780] = v[3379] + v[3292] * v[4446];
v[4922] = v[3389] + v[4780];
v[4779] = v[3375] + v[3292] * v[4449];
v[4926] = v[3393] + v[4779];
v[3293] = v[333] * (v[4768] + v[4769]) - v[3219] * v[722];
v[3294] = v[332] * (v[4768] + v[4770]) - v[3218] * v[720];
v[3295] = -(v[331] * (-v[4769] - v[4770])) - v[3217] * v[718];
v[7614] = -v[3295];
v[7615] = -v[3294];
v[7616] = -v[3293];
v[7617] = v[3383] + v[3291] * v[4446] + v[3290] * v[4449];
v[7618] = v[1784] * v[3290] + v[1785] * v[3292] + v[3379];
v[7619] = v[1786] * v[3291] + v[1787] * v[3292] + v[3375];
v[7620] = v[3295];
v[7621] = v[3294];
v[7622] = v[3293];
v[7623] = v[3350] + v[3288] * v[4454] + v[3287] * v[4457];
v[7624] = v[1792] * v[3287] + v[1793] * v[3289] + v[3346];
v[7625] = v[1794] * v[3288] + v[1795] * v[3289] + v[3342];
v[3299] = -(v[331] * v[4818]) - v[3277] * v[718];
v[3300] = -(v[332] * v[4822]) - v[3255] * v[720];
v[3301] = -(v[333] * v[4827]) - v[3233] * v[722];
v[3302] = -(v[331] * v[4817]) - v[3278] * v[718];
v[3303] = -(v[332] * v[4821]) - v[3256] * v[720];
v[3304] = -(v[333] * v[4826]) - v[3234] * v[722];
v[3305] = -(v[331] * v[4819]) - v[3279] * v[718];
v[3306] = -(v[332] * v[4823]) - v[3257] * v[720];
v[3307] = -(v[333] * v[4828]) - v[3235] * v[722];
v[3308] = v[331] * v[4816] + v[3280] * v[718];
v[3309] = v[332] * v[4803] + v[3258] * v[720];
v[3310] = v[333] * v[4795] + v[3236] * v[722];
v[3311] = v[331] * v[4815] + v[3281] * v[718];
v[3312] = v[332] * v[4804] + v[3259] * v[720];
v[3313] = v[333] * v[4797] + v[3237] * v[722];
v[3227] = v[3227] + v[2332] * v[3277] + v[2334] * v[3278] + v[2336] * v[3279] + v[2338] * v[3280] + v[2340] * v[3281]
+ v[2342] * v[3282] + v[3905] * v[4608] + v[332] * v[4813] + v[333] * v[4814];
v[3314] = v[331] * v[4811] + v[3282] * v[718];
v[3315] = v[332] * v[4802] + v[3260] * v[720];
v[3316] = v[333] * v[4796] + v[3238] * v[722];
v[3228] = v[3228] + v[3255] * v[4463] + v[3256] * v[4465] + v[3257] * v[4467] + v[3258] * v[4469] + v[3259] * v[4471]
+ v[3260] * v[4473] + v[3859] * v[4607] + v[331] * v[4800] + v[333] * v[4801];
v[3229] = v[3229] + v[2460] * v[3233] + v[2462] * v[3234] + v[2464] * v[3235] + v[2466] * v[3236] + v[2468] * v[3237]
+ v[2470] * v[3238] + v[3826] * v[4606] + v[331] * v[4793] + v[332] * v[4794];
b3317 = b329;
if (b3317) {
	v[3318] = -v[3229];
	v[3319] = -v[3228];
	v[3320] = -v[3227];
}
else {
	v[3318] = v[3229];
	v[3319] = v[3228];
	v[3320] = v[3227];
};
v[3325] = v[304] * v[3318] + v[303] * v[3319] + v[302] * v[3320];
v[3230] = v[3230] + v[315] * v[3325];
v[4771] = v[3230] / v[732];
v[3327] = v[316] * v[3318] + v[304] * v[4771] + v[756];
v[3328] = v[316] * v[3319] + v[303] * v[4771] + v[755];
v[3329] = v[316] * v[3320] + v[302] * v[4771] + v[754];
v[3330] = v[2275] * v[3287];
v[3331] = v[2274] * v[3287];
v[3332] = v[2273] * v[3287];
v[3333] = v[2270] * v[3288];
v[3334] = v[288] * v[3287] + v[290] * v[3288];
v[4772] = v[209] * v[3334];
v[4382] = -(v[1675] * v[4772]);
v[4379] = v[3334] * v[4528];
v[3335] = v[2268] * v[3288];
v[3336] = v[3330] + v[3333];
v[3337] = v[2269] * v[3288] + v[4772] / v[279];
v[3338] = v[2263] * v[3289];
v[3339] = v[282] * v[3287] + v[290] * v[3289];
v[4776] = v[214] * v[3339];
v[4381] = -(v[1675] * v[4776]);
v[4378] = v[3339] * v[4529];
v[3340] = v[282] * v[3288] + v[288] * v[3289];
v[4777] = v[218] * v[3340];
v[5101] = -(v[3287] * v[4374]) - v[3288] * v[4375] - v[3289] * v[4376] + v[4453] * v[4772] + v[278] * v[4776]
+ v[277] * v[4777];
v[4380] = -(v[1675] * v[4777]);
v[4377] = v[3340] * v[4531];
v[4372] = v[2272] * v[3287] + v[2267] * v[3288] + v[2262] * v[3289] + v[3332] + v[2073] * v[3334] + v[3335] + v[3338]
+ v[2071] * v[3339] + v[2067] * v[3340];
v[3341] = v[290] * v[4938] + v[3327] * v[871];
v[3344] = v[282] * v[4773] + v[3327] * v[875];
v[3345] = v[288] * v[4934] + v[3328] * v[873];
v[3348] = v[282] * v[4774] + v[3328] * v[875];
v[3349] = v[1793] * v[3334] + v[290] * v[4940] + v[3329] * v[871];
v[3351] = v[288] * v[4775] + v[3329] * v[873];
v[3354] = v[282] * v[4931] + v[3329] * v[875];
v[3355] = v[2265] * v[3289] + v[4776] / v[279];
v[3357] = v[3339] * v[4454] + v[290] * v[4939] + v[3328] * v[871];
v[3358] = v[3337] + v[3355];
v[3359] = v[2264] * v[3289] + v[4777] / v[279];
v[3361] = v[3340] * v[4457] + v[288] * v[4935] + v[3327] * v[873];
v[3362] = v[3331] + v[3359];
v[3363] = v[2260] * v[3290];
v[3364] = v[2259] * v[3290];
v[3365] = v[2258] * v[3290];
v[3366] = v[2255] * v[3291];
v[3367] = v[262] * v[3290] + v[264] * v[3291];
v[4778] = v[153] * v[3367];
v[4407] = -(v[1759] * v[4778]);
v[4404] = v[3367] * v[4549];
v[3368] = v[2253] * v[3291];
v[3369] = v[3363] + v[3366];
v[3370] = v[2254] * v[3291] + v[4778] / v[253];
v[3371] = v[2248] * v[3292];
v[3372] = v[256] * v[3290] + v[264] * v[3292];
v[4782] = v[158] * v[3372];
v[4406] = -(v[1759] * v[4782]);
v[4403] = v[3372] * v[4550];
v[3373] = v[256] * v[3291] + v[262] * v[3292];
v[4783] = v[162] * v[3373];
v[5105] = -(v[3290] * v[4399]) - v[3291] * v[4400] - v[3292] * v[4401] + v[4445] * v[4778] + v[252] * v[4782]
+ v[251] * v[4783];
v[4405] = -(v[1759] * v[4783]);
v[4402] = v[3373] * v[4552];
v[4397] = v[2257] * v[3290] + v[2252] * v[3291] + v[2247] * v[3292] + v[3365] + v[2045] * v[3367] + v[3368] + v[3371]
+ v[2043] * v[3372] + v[2039] * v[3373];
v[3374] = v[264] * v[4926] + v[3327] * v[877];
v[3377] = v[256] * v[4779] + v[3327] * v[881];
v[3378] = v[262] * v[4922] + v[3328] * v[879];
v[3381] = v[256] * v[4780] + v[3328] * v[881];
v[3382] = v[1785] * v[3367] + v[264] * v[4928] + v[3329] * v[877];
v[3384] = v[262] * v[4781] + v[3329] * v[879];
v[3387] = v[256] * v[4919] + v[3329] * v[881];
v[3388] = v[2250] * v[3292] + v[4782] / v[253];
v[3390] = v[3372] * v[4446] + v[264] * v[4927] + v[3328] * v[877];
v[3391] = v[3370] + v[3388];
v[3392] = v[2249] * v[3292] + v[4783] / v[253];
v[3394] = v[3373] * v[4449] + v[262] * v[4923] + v[3327] * v[879];
v[3395] = v[3364] + v[3392];
v[3396] = -(v[3302] * v[875]);
v[3397] = -(v[3302] * v[873]);
v[3398] = -(v[3302] * v[871]);
v[3399] = -(v[3299] * v[875]);
v[3400] = -(v[3299] * v[871]);
v[3401] = -(v[3299] * v[873]);
v[3402] = -(v[3305] * v[873]);
v[3403] = -(v[3305] * v[871]);
v[3404] = -(v[3305] * v[875]);
v[3405] = -(v[3306] * v[875]);
v[3406] = -(v[3306] * v[873]);
v[3407] = -(v[3306] * v[871]);
v[3408] = -(v[3300] * v[871]);
v[3409] = -(v[3300] * v[875]);
v[3410] = -(v[3300] * v[873]);
v[3411] = -(v[3301] * v[873]);
v[3412] = -(v[3301] * v[875]);
v[3413] = -(v[3301] * v[871]);
v[3414] = v[3402] + v[3405] + v[3408] + v[3411] + v[3335] * v[4596] - v[3358] * v[489] - v[3336] * v[492];
v[3415] = -(v[3303] * v[875]);
v[3416] = -(v[3303] * v[871]);
v[3417] = -(v[3303] * v[873]);
v[3418] = -v[3397] - v[3400] - v[3412] - v[3415] + v[3338] * v[4595] - v[3358] * v[486] + v[3362] * v[492];
v[3419] = v[204] * v[3351] - v[3402] * v[479] + v[3397] * v[480] - v[3401] * v[481];
v[3420] = -(v[3307] * v[875]);
v[3421] = -(v[3307] * v[873]);
v[3422] = -(v[3307] * v[871]);
v[3423] = -(v[3304] * v[873]);
v[3424] = -(v[3304] * v[875]);
v[3425] = -(v[3304] * v[871]);
v[3426] = -(v[231] * v[3327]) - v[228] * v[3328] - v[225] * v[3329] + v[3305] * v[530] + v[3302] * v[531]
+ v[3299] * v[532] + v[3306] * v[545] + v[3303] * v[546] + v[3300] * v[547] + v[3307] * v[560] + v[3304] * v[561]
+ v[3301] * v[562];
v[3427] = -(v[230] * v[3327]) - v[227] * v[3328] - v[224] * v[3329] + v[3305] * v[527] + v[3302] * v[528]
+ v[3299] * v[529] + v[3306] * v[542] + v[3303] * v[543] + v[3300] * v[544] + v[3307] * v[557] + v[3304] * v[558]
+ v[3301] * v[559];
v[3428] = -(v[229] * v[3327]) - v[226] * v[3328] - v[223] * v[3329] + v[3305] * v[524] + v[3302] * v[525]
+ v[3299] * v[526] + v[3306] * v[539] + v[3303] * v[540] + v[3300] * v[541] + v[3307] * v[554] + v[3304] * v[555]
+ v[3301] * v[556];
v[3429] = -v[3403] - v[3416] - v[3420] - v[3423] + v[3332] * v[4594] - v[3336] * v[486] + v[3362] * v[489];
v[3430] = v[204] * v[3349] - v[3403] * v[479] + v[3398] * v[480] - v[3400] * v[481];
v[3431] = v[204] * v[3348] - v[3405] * v[479] + v[3415] * v[480] - v[3409] * v[481];
v[3432] = v[3399] + v[3410];
v[3433] = v[204] * v[3357] - v[3407] * v[479] + v[3416] * v[480] - v[3408] * v[481];
v[3434] = v[204] * v[3344] - v[3420] * v[479] + v[3424] * v[480] - v[3412] * v[481];
v[3435] = v[204] * v[3361] - v[3421] * v[479] + v[3423] * v[480] - v[3411] * v[481];
v[3436] = v[3406] + v[3422];
v[3437] = v[3396] + v[3425];
v[7882] = 0e0;
v[7883] = 0e0;
v[7884] = 0e0;
v[7885] = 0e0;
v[7886] = 0e0;
v[7887] = 0e0;
v[7888] = 0e0;
v[7889] = 0e0;
v[7890] = 0e0;
v[7891] = -v[3418] / 2e0 - v[3436];
v[7892] = v[3414] / 2e0 - v[3437];
v[7893] = -v[3429] / 2e0 - v[3432];
v[3438] = v[3398] - v[3401] - v[3407] + v[3409] + v[3421] - v[3424] + v[3414] * v[375] - v[3418] * v[376]
- v[3429] * v[378] + v[3341] * v[4537] + v[3345] * v[4538] + v[3354] * v[4539] + v[4372] * v[4586] + v[3432] * v[4718]
+ v[3436] * v[4719] + v[3437] * v[4720] + v[3351] * v[487] + v[3349] * v[493] + v[3348] * v[497] + v[3357] * v[506]
+ v[3344] * v[510] + v[3361] * v[514];
v[3439] = v[3428] * v[80] + v[3427] * v[81] + v[3426] * v[82];
v[3440] = v[3311] * v[881];
v[3441] = v[3311] * v[879];
v[3442] = v[3311] * v[877];
v[3443] = v[3308] * v[881];
v[3444] = v[3308] * v[877];
v[3445] = v[3308] * v[879];
v[3446] = v[3314] * v[879];
v[3447] = v[3314] * v[877];
v[3448] = v[3314] * v[881];
v[3449] = v[3315] * v[881];
v[3450] = v[3315] * v[879];
v[3451] = v[3315] * v[877];
v[3452] = v[3309] * v[877];
v[3453] = v[3309] * v[881];
v[3454] = v[3309] * v[879];
v[3455] = v[3310] * v[879];
v[3456] = v[3310] * v[881];
v[3457] = v[3310] * v[877];
v[3458] = v[3446] + v[3449] + v[3452] + v[3455] - v[3391] * v[390] - v[3369] * v[393] + v[3368] * v[4592];
v[3459] = v[3312] * v[881];
v[3460] = v[3312] * v[877];
v[3461] = v[3312] * v[879];
v[3462] = -v[3441] - v[3444] - v[3456] - v[3459] - v[3391] * v[387] + v[3395] * v[393] + v[3371] * v[4591];
v[3463] = v[148] * v[3384] - v[3446] * v[380] + v[3441] * v[381] - v[3445] * v[382];
v[3464] = v[3316] * v[881];
v[3465] = v[3316] * v[879];
v[3466] = v[3316] * v[877];
v[3467] = v[3313] * v[879];
v[3468] = v[3313] * v[881];
v[3469] = v[3313] * v[877];
v[3470] = v[175] * v[3327] + v[172] * v[3328] + v[169] * v[3329] + v[3314] * v[431] + v[3311] * v[432] + v[3308] * v[433]
+ v[3315] * v[446] + v[3312] * v[447] + v[3309] * v[448] + v[3316] * v[461] + v[3313] * v[462] + v[3310] * v[463];
v[3471] = v[174] * v[3327] + v[171] * v[3328] + v[168] * v[3329] + v[3314] * v[428] + v[3311] * v[429] + v[3308] * v[430]
+ v[3315] * v[443] + v[3312] * v[444] + v[3309] * v[445] + v[3316] * v[458] + v[3313] * v[459] + v[3310] * v[460];
v[3472] = v[173] * v[3327] + v[170] * v[3328] + v[167] * v[3329] + v[3314] * v[425] + v[3311] * v[426] + v[3308] * v[427]
+ v[3315] * v[440] + v[3312] * v[441] + v[3309] * v[442] + v[3316] * v[455] + v[3313] * v[456] + v[3310] * v[457];
v[3473] = -v[3447] - v[3460] - v[3464] - v[3467] - v[3369] * v[387] + v[3395] * v[390] + v[3365] * v[4590];
v[3474] = v[148] * v[3382] - v[3447] * v[380] + v[3442] * v[381] - v[3444] * v[382];
v[3475] = v[148] * v[3381] - v[3449] * v[380] + v[3459] * v[381] - v[3453] * v[382];
v[3476] = v[3443] + v[3454];
v[3477] = v[148] * v[3390] - v[3451] * v[380] + v[3460] * v[381] - v[3452] * v[382];
v[3478] = v[148] * v[3377] - v[3464] * v[380] + v[3468] * v[381] - v[3456] * v[382];
v[3479] = v[148] * v[3394] - v[3465] * v[380] + v[3467] * v[381] - v[3455] * v[382];
v[3480] = v[3450] + v[3466];
v[3481] = v[3440] + v[3469];
v[7894] = 0e0;
v[7895] = 0e0;
v[7896] = 0e0;
v[7897] = -v[3462] / 2e0 - v[3480];
v[7898] = v[3458] / 2e0 - v[3481];
v[7899] = -v[3473] / 2e0 - v[3476];
v[7900] = 0e0;
v[7901] = 0e0;
v[7902] = 0e0;
v[7903] = 0e0;
v[7904] = 0e0;
v[7905] = 0e0;
v[3482] = v[3442] - v[3445] - v[3451] + v[3453] + v[3465] - v[3468] + v[3458] * v[369] - v[3462] * v[370]
- v[3473] * v[372] + v[3384] * v[388] + v[3382] * v[394] + v[3381] * v[398] + v[3390] * v[407] + v[3377] * v[411]
+ v[3394] * v[415] + v[3374] * v[4558] + v[3378] * v[4559] + v[3387] * v[4560] + v[4397] * v[4585] + v[3476] * v[4730]
+ v[3480] * v[4731] + v[3481] * v[4732];
v[3483] = v[3472] * v[37] + v[3471] * v[38] + v[3470] * v[39];
v[3484] = (-v[3439] + v[3428] * v[86] + v[3427] * v[87] + v[3426] * v[88]) / 2e0;
v[3485] = (-v[3439] + v[3428] * v[83] + v[3427] * v[84] + v[3426] * v[85]) / 2e0;
v[3486] = (-2e0*v[3330] + 2e0*v[3333] + v[3403] * v[4561] + v[3402] * v[4562] + v[3407] * v[4563] + v[3405] * v[4564]
	+ v[3421] * v[4565] + v[3420] * v[4566] - v[3404] * v[482] - v[3406] * v[501] - v[3422] * v[519]) / 2e0;
v[3488] = -v[3331] + v[3359] + v[3425] * v[4537] + v[3417] * v[4538] + v[3396] * v[4539] + v[3397] * v[487]
+ v[3398] * v[493] + v[3415] * v[497] + v[3416] * v[506] + v[3424] * v[510] + v[3423] * v[514];
v[3489] = (v[204] * v[3345] - v[3406] * v[479] + v[3417] * v[480] - v[3410] * v[481]) / 2e0;
v[3490] = (-2e0*v[3337] + 2e0*v[3355] + v[3400] * v[4561] + v[3401] * v[4562] + v[3408] * v[4563] + v[3409] * v[4564]
	+ v[3411] * v[4565] + v[3412] * v[4566] - v[3399] * v[482] - v[3410] * v[501] - v[3413] * v[519]) / 2e0;
v[5099] = v[3486] * v[374] - v[3488] * v[377] + v[3490] * v[379];
v[7906] = 0e0;
v[7907] = 0e0;
v[7908] = 0e0;
v[7909] = 0e0;
v[7910] = 0e0;
v[7911] = 0e0;
v[7912] = 0e0;
v[7913] = 0e0;
v[7914] = 0e0;
v[7915] = 8e0*v[3486];
v[7916] = -8e0*v[3488];
v[7917] = 8e0*v[3490];
v[3507] = v[1870] * v[3486] + v[1868] * v[3488] - v[3489] + v[1865] * v[3490] - v[3438] * v[4428];
v[4419] = v[3507] + (-(v[204] * v[3354]) + v[3404] * v[479] - v[3396] * v[480] + v[3399] * v[481]) / 2e0;
v[3491] = (v[204] * v[3341] - v[3422] * v[479] + v[3425] * v[480] - v[3413] * v[481]) / 2e0;
v[4418] = v[3489] - v[3491] + v[4419];
v[4416] = -v[3491] + v[3507];
v[3492] = v[3430] + v[3434];
v[3493] = v[3433] + v[3435];
v[3494] = v[3419] + v[3431];
v[3495] = (-v[3483] + v[3472] * v[43] + v[3471] * v[44] + v[3470] * v[45]) / 2e0;
v[3496] = (-v[3483] + v[3472] * v[40] + v[3471] * v[41] + v[3470] * v[42]) / 2e0;
v[3497] = (-2e0*v[3363] + 2e0*v[3366] - v[3448] * v[383] - v[3450] * v[402] - v[3466] * v[420] + v[3447] * v[4567]
	+ v[3446] * v[4568] + v[3451] * v[4569] + v[3449] * v[4570] + v[3465] * v[4571] + v[3464] * v[4572]) / 2e0;
v[3499] = -v[3364] + v[3392] + v[3441] * v[388] + v[3442] * v[394] + v[3459] * v[398] + v[3460] * v[407] + v[3468] * v[411]
+ v[3467] * v[415] + v[3469] * v[4558] + v[3461] * v[4559] + v[3440] * v[4560];
v[3500] = (v[148] * v[3378] - v[3450] * v[380] + v[3461] * v[381] - v[3454] * v[382]) / 2e0;
v[3501] = (-2e0*v[3370] + 2e0*v[3388] - v[3443] * v[383] - v[3454] * v[402] - v[3457] * v[420] + v[3444] * v[4567]
	+ v[3445] * v[4568] + v[3452] * v[4569] + v[3453] * v[4570] + v[3455] * v[4571] + v[3456] * v[4572]) / 2e0;
v[5103] = v[3497] * v[368] - v[3499] * v[371] + v[3501] * v[373];
v[7990] = 0e0;
v[7991] = 0e0;
v[7992] = 0e0;
v[7993] = 8e0*v[3497];
v[7994] = -8e0*v[3499];
v[7995] = 8e0*v[3501];
v[7996] = 0e0;
v[7997] = 0e0;
v[7998] = 0e0;
v[7999] = 0e0;
v[8000] = 0e0;
v[8001] = 0e0;
v[3506] = v[1844] * v[3497] + v[1842] * v[3499] - v[3500] + v[1839] * v[3501] - v[3482] * v[4426];
v[4415] = v[3506] + (-(v[148] * v[3387]) + v[3448] * v[380] - v[3440] * v[381] + v[3443] * v[382]) / 2e0;
v[3502] = (v[148] * v[3374] - v[3466] * v[380] + v[3469] * v[381] - v[3457] * v[382]) / 2e0;
v[4414] = v[3500] - v[3502] + v[4415];
v[4412] = -v[3502] + v[3506];
v[3503] = v[3474] + v[3478];
v[3504] = v[3477] + v[3479];
v[3505] = v[3463] + v[3475];
v[7602] = v[3329];
v[7603] = v[3328];
v[7604] = v[3327];
v[7605] = -v[3477] + v[3479] + v[3505] * v[369] + v[3503] * v[372] + v[368] * v[4412] + v[3462] * v[4430] + 2e0*
(v[3497] * v[4426] + v[3480] * v[4430]);
v[7606] = v[3474] - v[3478] + v[3505] * v[370] + v[3504] * v[372] + v[371] * v[4414] - v[3458] * v[4430] + 2e0*(-
(v[3499] * v[4426]) + v[3481] * v[4430]);
v[7607] = -v[3463] + v[3475] + v[3504] * v[369] + v[3503] * v[370] + v[373] * v[4415] + v[3473] * v[4430] + 2e0*
(v[3501] * v[4426] + v[3476] * v[4430]);
v[7608] = -v[3329];
v[7609] = -v[3328];
v[7610] = -v[3327];
v[7611] = -v[3433] + v[3435] + v[3494] * v[375] + v[3492] * v[378] + v[374] * v[4416] + v[3418] * v[4437] + 2e0*
(v[3486] * v[4428] + v[3436] * v[4437]);
v[7612] = v[3430] - v[3434] + v[3494] * v[376] + v[3493] * v[378] + v[377] * v[4418] - v[3414] * v[4437] + 2e0*(-
(v[3488] * v[4428]) + v[3437] * v[4437]);
v[7613] = -v[3419] + v[3431] + v[3493] * v[375] + v[3492] * v[376] + v[379] * v[4419] + v[3429] * v[4437] + 2e0*
(v[3490] * v[4428] + v[3432] * v[4437]);
for (i3212 = 1; i3212 <= 12; i3212++) {
	v[3570] = (i3212 == 11 ? 1 : 0);
	v[4869] = v[10] * v[3570];
	v[3567] = (i3212 == 10 ? 1 : 0);
	v[4870] = v[10] * v[3567];
	v[3564] = (i3212 == 12 ? 1 : 0);
	v[4874] = v[10] * v[3564];
	v[3544] = (i3212 == 5 ? 1 : 0);
	v[4878] = v[10] * v[3544];
	v[3541] = (i3212 == 4 ? 1 : 0);
	v[4879] = v[10] * v[3541];
	v[3538] = (i3212 == 6 ? 1 : 0);
	v[4883] = v[10] * v[3538];
	v[3532] = v[6181 + i3212];
	v[3531] = v[6193 + i3212];
	v[3529] = v[6157 + i3212];
	v[3528] = v[6169 + i3212];
	v[3512] = v[6233 + i3212];
	v[3513] = v[6257 + i3212];
	v[3514] = v[6245 + i3212];
	v[3516] = v[6745 + i3212];
	v[3517] = v[6221 + i3212];
	v[3554] = v[3517] * v[4426];
	v[3622] = -(v[3554] * v[4585]);
	v[4894] = v[264] * v[3622];
	v[4892] = v[262] * v[3622];
	v[4788] = v[253] * v[3622];
	v[3519] = v[6793 + i3212];
	v[3520] = v[6281 + i3212];
	v[3521] = v[6305 + i3212];
	v[3522] = v[6293 + i3212];
	v[3524] = v[6817 + i3212];
	v[3525] = v[6269 + i3212];
	v[3580] = v[3525] * v[4428];
	v[3671] = -(v[3580] * v[4586]);
	v[4889] = v[290] * v[3671];
	v[4887] = v[288] * v[3671];
	v[4789] = v[279] * v[3671];
	v[3527] = v[6865 + i3212];
	v[4784] = -v[3528] - v[3529];
	v[4786] = -v[3531] - v[3532];
	v[3534] = v[10] * v[7629 + i3212];
	v[3535] = v[10] * v[7641 + i3212];
	v[4820] = -(v[3218] * v[3534]) - v[3217] * v[3535];
	v[3536] = v[10] * v[7653 + i3212];
	v[4825] = -(v[3219] * v[3534]) - v[3217] * v[3536];
	v[4824] = -(v[3219] * v[3535]) - v[3218] * v[3536];
	v[3537] = v[3512] + v[3538];
	v[4907] = 2e0*v[3537];
	v[3539] = v[3512] - v[3538];
	v[4908] = 2e0*v[3539];
	v[3540] = v[3513] + v[3541];
	v[4909] = 2e0*v[3540];
	v[3542] = v[3513] - v[3541];
	v[4910] = 2e0*v[3542];
	v[3543] = v[3514] - v[3544];
	v[4911] = 2e0*v[3543];
	v[3545] = v[3514] + v[3544];
	v[4912] = 2e0*v[3545];
	v[3546] = v[1839] * v[3517] + v[3538] * v[4587];
	v[3547] = -v[3517] + v[3544] * v[371];
	v[3548] = v[1842] * v[3517] - v[3544] * v[4587];
	v[3549] = v[1844] * v[3517] + v[3541] * v[4587];
	v[4785] = 2e0*(-(v[148] * v[3544]) - v[3517] * v[4429]);
	v[3551] = -(v[148] * v[3541]) - v[3517] * v[4431];
	v[3552] = -(v[148] * v[3538]) - v[3517] * v[4432];
	v[3553] = v[3554] * v[372] + v[3538] * v[4430];
	v[3555] = v[3554] * v[370] + v[3541] * v[4430];
	v[3556] = -(v[3554] * v[369]) - v[3544] * v[4430];
	v[3557] = (v[148] * v[3516] - v[3554] * v[420]) / 2e0;
	v[3706] = v[264] * v[3557];
	v[3558] = (-(v[3516] * v[382]) - v[3546] * v[420]) / 2e0;
	v[3559] = (v[148] * v[3547] - v[3554] * v[402]) / 2e0;
	v[3700] = v[262] * v[3559];
	v[3560] = (v[3547] * v[381] + v[3548] * v[402]) / 2e0;
	v[3561] = (v[148] * v[3519] - v[3554] * v[383]) / 2e0;
	v[3695] = v[256] * v[3561];
	v[3562] = (-(v[3519] * v[380]) - v[3549] * v[383]) / 2e0;
	v[3563] = v[3520] + v[3564];
	v[4897] = 2e0*v[3563];
	v[3565] = v[3520] - v[3564];
	v[4898] = 2e0*v[3565];
	v[3566] = v[3521] + v[3567];
	v[4899] = 2e0*v[3566];
	v[3568] = v[3521] - v[3567];
	v[4900] = 2e0*v[3568];
	v[3569] = v[3522] - v[3570];
	v[4901] = 2e0*v[3569];
	v[3571] = v[3522] + v[3570];
	v[4902] = 2e0*v[3571];
	v[3572] = v[1865] * v[3525] + v[3564] * v[4588];
	v[3573] = -v[3525] + v[3570] * v[377];
	v[3574] = v[1868] * v[3525] - v[3570] * v[4588];
	v[3575] = v[1870] * v[3525] + v[3567] * v[4588];
	v[4787] = 2e0*(-(v[204] * v[3570]) - v[3525] * v[4436]);
	v[3577] = -(v[204] * v[3567]) - v[3525] * v[4438];
	v[3578] = -(v[204] * v[3564]) - v[3525] * v[4439];
	v[3579] = v[3580] * v[378] + v[3564] * v[4437];
	v[3581] = v[3580] * v[376] + v[3567] * v[4437];
	v[3582] = -(v[3580] * v[375]) - v[3570] * v[4437];
	v[3583] = (v[204] * v[3524] - v[3580] * v[519]) / 2e0;
	v[3755] = v[290] * v[3583];
	v[3584] = (-(v[3524] * v[481]) - v[3572] * v[519]) / 2e0;
	v[3585] = (v[204] * v[3573] - v[3580] * v[501]) / 2e0;
	v[3747] = v[288] * v[3585];
	v[3586] = (v[3573] * v[480] + v[3574] * v[501]) / 2e0;
	v[3587] = (v[204] * v[3527] - v[3580] * v[482]) / 2e0;
	v[3740] = v[282] * v[3587];
	v[3588] = (-(v[3527] * v[479]) - v[3575] * v[482]) / 2e0;
	v[3589] = (v[3529] * v[40] + v[3528] * v[43] + v[37] * v[4784]) / 2e0;
	v[3590] = (v[3529] * v[41] + v[3528] * v[44] + v[38] * v[4784]) / 2e0;
	v[3591] = (v[3529] * v[42] + v[3528] * v[45] + v[39] * v[4784]) / 2e0;
	v[3592] = (v[3516] * v[381] + v[3548] * v[420] + v[4785]) / 2e0;
	v[3593] = (v[3519] * v[381] + v[3548] * v[383] + v[4785]) / 2e0;
	v[3594] = v[3551] + v[3516] * v[4431] - v[3549] * v[4558];
	v[3595] = v[3551] + v[3547] * v[4431] - v[3549] * v[4559];
	v[3596] = -v[3554] - v[3540] * v[380] - v[3549] * v[415];
	v[3597] = v[148] * v[3540] - v[3554] * v[415];
	v[3598] = v[3554] + v[3543] * v[381] + v[3548] * v[411];
	v[3599] = v[148] * v[3543] - v[3554] * v[411];
	v[3708] = v[256] * v[3599];
	v[3600] = v[3554] - v[3542] * v[380] - v[3549] * v[407];
	v[3601] = v[148] * v[3542] - v[3554] * v[407];
	v[3701] = v[264] * v[3601];
	v[3602] = v[3552] + v[3547] * v[4432] - v[3546] * v[4559];
	v[3603] = v[3552] + v[3519] * v[4432] - v[3546] * v[4560];
	v[3604] = -v[3554] - v[3537] * v[382] - v[3546] * v[398];
	v[3605] = v[148] * v[3537] - v[3554] * v[398];
	v[3606] = -v[3554] + v[3545] * v[381] + v[3548] * v[394];
	v[3607] = v[148] * v[3545] - v[3554] * v[394];
	v[3696] = v[264] * v[3607];
	v[3608] = -v[3553] + v[3540] * v[381] + v[3548] * v[415];
	v[3609] = -v[3553] - v[3543] * v[380] - v[3549] * v[411];
	v[3610] = -v[3553] + v[3542] * v[381] + v[3548] * v[407];
	v[3611] = -v[3553] - v[3545] * v[380] - v[3549] * v[394];
	v[3612] = v[3622] + v[3553] * v[4590];
	v[3613] = v[3589] * v[456] + v[3590] * v[459] + v[3591] * v[462] + v[3592] * v[877] + v[3608] * v[879] + v[3598] * v[881];
	v[4845] = v[333] * v[3613];
	v[3614] = v[3589] * v[455] + v[3590] * v[458] + v[3591] * v[461] + v[3594] * v[877] + v[3596] * v[879] + v[3609] * v[881];
	v[4805] = v[333] * v[3614];
	v[3615] = v[3554] - v[3539] * v[382] - v[3546] * v[388];
	v[3616] = v[148] * v[3539] - v[3554] * v[388];
	v[3617] = -v[3555] + v[3537] * v[381] + v[3548] * v[398];
	v[3618] = -v[3555] - v[3543] * v[382] - v[3546] * v[411];
	v[3619] = -v[3555] - v[3545] * v[382] - v[3546] * v[394];
	v[3620] = -v[3555] + v[3539] * v[381] + v[3548] * v[388];
	v[3621] = v[3553] * v[390] + v[3555] * v[393];
	v[3623] = v[3622] + v[3555] * v[4591];
	v[3722] = v[3292] * v[3623];
	v[3624] = v[3589] * v[441] + v[3590] * v[444] + v[3591] * v[447] + v[3610] * v[877] + v[3560] * v[879] + v[3617] * v[881];
	v[4846] = v[332] * v[3624];
	v[3625] = v[3556] - v[3540] * v[382] - v[3546] * v[415];
	v[3626] = v[3556] - v[3542] * v[382] - v[3546] * v[407];
	v[3627] = v[3556] - v[3537] * v[380] - v[3549] * v[398];
	v[3628] = v[3556] - v[3539] * v[380] - v[3549] * v[388];
	v[3629] = -(v[3555] * v[387]) - v[3556] * v[390];
	v[3630] = -(v[3553] * v[387]) - v[3556] * v[393];
	v[3631] = v[3622] + v[3556] * v[4592];
	v[3726] = v[3291] * v[3631];
	v[3632] = v[3589] * v[457] + v[3590] * v[460] + v[3591] * v[463] + v[3558] * v[877] + v[3625] * v[879] + v[3618] * v[881];
	v[4842] = v[333] * v[3632];
	v[3633] = v[3589] * v[442] + v[3590] * v[445] + v[3591] * v[448] + v[3626] * v[877] + v[3602] * v[879] + v[3604] * v[881];
	v[4843] = v[332] * v[3633];
	v[3634] = v[3589] * v[440] + v[3590] * v[443] + v[3591] * v[446] + v[3600] * v[877] + v[3595] * v[879] + v[3627] * v[881];
	v[4808] = v[332] * v[3634];
	v[4810] = v[4805] + v[4808];
	v[3635] = v[3589] * v[425] + v[3590] * v[428] + v[3591] * v[431] + v[3611] * v[877] + v[3628] * v[879] + v[3562] * v[881];
	v[4812] = v[331] * v[3635];
	v[3636] = v[3589] * v[427] + v[3590] * v[430] + v[3591] * v[433] + v[3619] * v[877] + v[3615] * v[879] + v[3603] * v[881];
	v[4844] = v[331] * v[3636];
	v[3637] = v[3589] * v[426] + v[3590] * v[429] + v[3591] * v[432] + v[3606] * v[877] + v[3620] * v[879] + v[3593] * v[881];
	v[4847] = v[331] * v[3637];
	v[3638] = (v[4786] * v[80] + v[3532] * v[83] + v[3531] * v[86]) / 2e0;
	v[3639] = (v[4786] * v[81] + v[3532] * v[84] + v[3531] * v[87]) / 2e0;
	v[3640] = (v[4786] * v[82] + v[3532] * v[85] + v[3531] * v[88]) / 2e0;
	v[3641] = (v[4787] + v[3524] * v[480] + v[3574] * v[519]) / 2e0;
	v[3642] = (v[4787] + v[3527] * v[480] + v[3574] * v[482]) / 2e0;
	v[3643] = v[3577] + v[3524] * v[4438] - v[3575] * v[4537];
	v[3644] = v[3577] + v[3573] * v[4438] - v[3575] * v[4538];
	v[3645] = -v[3580] - v[3566] * v[479] - v[3575] * v[514];
	v[3646] = v[204] * v[3566] - v[3580] * v[514];
	v[3647] = v[3580] + v[3569] * v[480] + v[3574] * v[510];
	v[3648] = v[204] * v[3569] - v[3580] * v[510];
	v[3757] = v[282] * v[3648];
	v[3649] = v[3580] - v[3568] * v[479] - v[3575] * v[506];
	v[3650] = v[204] * v[3568] - v[3580] * v[506];
	v[3748] = v[290] * v[3650];
	v[3651] = v[3578] + v[3573] * v[4439] - v[3572] * v[4538];
	v[3652] = v[3578] + v[3527] * v[4439] - v[3572] * v[4539];
	v[3653] = -v[3580] - v[3563] * v[481] - v[3572] * v[497];
	v[3654] = v[204] * v[3563] - v[3580] * v[497];
	v[3655] = -v[3580] + v[3571] * v[480] + v[3574] * v[493];
	v[3656] = v[204] * v[3571] - v[3580] * v[493];
	v[3741] = v[290] * v[3656];
	v[3657] = -v[3579] + v[3566] * v[480] + v[3574] * v[514];
	v[3658] = -v[3579] - v[3569] * v[479] - v[3575] * v[510];
	v[3659] = -v[3579] + v[3568] * v[480] + v[3574] * v[506];
	v[3660] = -v[3579] - v[3571] * v[479] - v[3575] * v[493];
	v[3661] = v[3671] + v[3579] * v[4594];
	v[3662] = v[3638] * v[555] + v[3639] * v[558] + v[3640] * v[561] - v[3641] * v[871] - v[3657] * v[873] - v[3647] * v[875];
	v[4836] = -(v[333] * v[3662]);
	v[3663] = v[3638] * v[554] + v[3639] * v[557] + v[3640] * v[560] - v[3643] * v[871] - v[3645] * v[873] - v[3658] * v[875];
	v[4839] = -(v[333] * v[3663]);
	v[3664] = v[3580] - v[3565] * v[481] - v[3572] * v[487];
	v[3665] = v[204] * v[3565] - v[3580] * v[487];
	v[3666] = -v[3581] + v[3563] * v[480] + v[3574] * v[497];
	v[3667] = -v[3581] - v[3569] * v[481] - v[3572] * v[510];
	v[3668] = -v[3581] - v[3571] * v[481] - v[3572] * v[493];
	v[3669] = -v[3581] + v[3565] * v[480] + v[3574] * v[487];
	v[3670] = v[3579] * v[489] + v[3581] * v[492];
	v[3672] = v[3671] + v[3581] * v[4595];
	v[3771] = v[3289] * v[3672];
	v[3673] = v[3638] * v[540] + v[3639] * v[543] + v[3640] * v[546] - v[3659] * v[871] - v[3586] * v[873] - v[3666] * v[875];
	v[4837] = -(v[332] * v[3673]);
	v[3674] = v[3582] - v[3566] * v[481] - v[3572] * v[514];
	v[3675] = v[3582] - v[3568] * v[481] - v[3572] * v[506];
	v[3676] = v[3582] - v[3563] * v[479] - v[3575] * v[497];
	v[3677] = v[3582] - v[3565] * v[479] - v[3575] * v[487];
	v[3678] = -(v[3581] * v[486]) - v[3582] * v[489];
	v[3679] = -(v[3579] * v[486]) - v[3582] * v[492];
	v[3680] = v[3671] + v[3582] * v[4596];
	v[3775] = v[3288] * v[3680];
	v[3681] = v[3638] * v[556] + v[3639] * v[559] + v[3640] * v[562] - v[3584] * v[871] - v[3674] * v[873] - v[3667] * v[875];
	v[4833] = -(v[333] * v[3681]);
	v[3682] = v[3638] * v[541] + v[3639] * v[544] + v[3640] * v[547] - v[3675] * v[871] - v[3651] * v[873] - v[3653] * v[875];
	v[4834] = -(v[332] * v[3682]);
	v[3683] = v[3638] * v[539] + v[3639] * v[542] + v[3640] * v[545] - v[3649] * v[871] - v[3644] * v[873] - v[3676] * v[875];
	v[4840] = -(v[332] * v[3683]);
	v[3684] = v[3638] * v[524] + v[3639] * v[527] + v[3640] * v[530] - v[3660] * v[871] - v[3677] * v[873] - v[3588] * v[875];
	v[4841] = -(v[331] * v[3684]);
	v[3685] = v[3638] * v[526] + v[3639] * v[529] + v[3640] * v[532] - v[3668] * v[871] - v[3664] * v[873] - v[3652] * v[875];
	v[4835] = -(v[331] * v[3685]);
	v[3686] = v[3638] * v[525] + v[3639] * v[528] + v[3640] * v[531] - v[3655] * v[871] - v[3669] * v[873] - v[3642] * v[875];
	v[4838] = -(v[331] * v[3686]);
	v[3687] = v[3548] + v[3621];
	v[3720] = v[3292] * v[3687];
	v[3688] = -v[3548] + v[3621];
	v[3689] = (v[251] * v[3597] + v[162] * v[3687] + v[2039] * v[4788]) / v[253];
	v[3690] = v[3546] + v[3629];
	v[3728] = v[3292] * v[3690];
	v[3691] = -v[3546] + v[3629];
	v[3724] = v[3291] * v[3691];
	v[3692] = (v[252] * v[3601] + v[158] * v[3690] + v[2043] * v[4788]) / v[253];
	v[3693] = v[262] * v[3616] + v[3695];
	v[3694] = v[3693] + v[3696] + v[4879];
	v[3697] = v[3695] + v[3696];
	v[3698] = v[256] * v[3605] + v[3700];
	v[3699] = v[3698] + v[3701] + v[4878];
	v[3702] = v[3700] + v[3701];
	v[3703] = v[3328] * v[3559] + v[3312] * v[3560] + v[3315] * v[3595] + v[3316] * v[3596] + v[3327] * v[3597]
		+ v[3309] * v[3602] + v[3313] * v[3608] + v[3308] * v[3615] + v[3329] * v[3616] + v[3311] * v[3620] + v[3310] * v[3625]
		+ v[3314] * v[3628];
	v[3704] = v[3329] * v[3561] + v[3314] * v[3562] + v[3311] * v[3593] + v[3313] * v[3598] + v[3327] * v[3599]
		+ v[3308] * v[3603] + v[3309] * v[3604] + v[3328] * v[3605] + v[3316] * v[3609] + v[3312] * v[3617] + v[3310] * v[3618]
		+ v[3315] * v[3627];
	v[3705] = v[3706] + v[3708];
	v[3707] = v[262] * v[3597] + v[3706];
	v[3709] = v[3707] + v[3708] + v[4883];
	v[3710] = v[3327] * v[3557] + v[3310] * v[3558] + v[3313] * v[3592] + v[3316] * v[3594] + v[3315] * v[3600]
		+ v[3328] * v[3601] + v[3311] * v[3606] + v[3329] * v[3607] + v[3312] * v[3610] + v[3314] * v[3611] + v[3308] * v[3619]
		+ v[3309] * v[3626];
	v[3711] = (v[153] * v[3691] + v[3607] * v[4445] + v[2045] * v[4788]) / v[253];
	v[3712] = v[3722] + v[3724];
	v[4924] = v[3712] / v[253];
	v[3713] = v[3549] + v[3630];
	v[3718] = v[3291] * v[3713];
	v[3714] = -v[3549] + v[3630];
	v[3715] = v[3726] + v[3728];
	v[4920] = v[3715] / v[253];
	v[3716] = v[3290] * v[3612] + v[3718] + v[3720];
	v[3719] = v[3716] - v[3720];
	v[3721] = v[3716] - v[3718];
	v[4921] = v[3721] / v[253];
	v[3723] = v[3290] * v[3688] + v[3722];
	v[3725] = v[3723] + v[3724];
	v[3727] = v[3290] * v[3714] + v[3726];
	v[3729] = v[3727] + v[3728];
	v[3730] = v[3574] + v[3670];
	v[3769] = v[3289] * v[3730];
	v[3731] = -v[3574] + v[3670];
	v[3732] = (v[277] * v[3646] + v[218] * v[3730] + v[2067] * v[4789]) / v[279];
	v[3733] = v[3572] + v[3678];
	v[3777] = v[3289] * v[3733];
	v[3734] = -v[3572] + v[3678];
	v[3773] = v[3288] * v[3734];
	v[3735] = (v[278] * v[3650] + v[214] * v[3733] + v[2071] * v[4789]) / v[279];
	v[3736] = v[288] * v[3665] + v[3740];
	v[3738] = v[167] * v[3589] + v[168] * v[3590] + v[169] * v[3591] - v[223] * v[3638] - v[224] * v[3639] - v[225] * v[3640]
		+ v[6533 + i3212] + v[3656] * v[871] + v[3665] * v[873] + v[3587] * v[875] + v[3607] * v[877] + v[3616] * v[879]
		+ v[3561] * v[881];
	v[3739] = v[3736] + v[3741] + v[4870];
	v[3742] = v[3740] + v[3741];
	v[3744] = v[170] * v[3589] + v[171] * v[3590] + v[172] * v[3591] - v[226] * v[3638] - v[227] * v[3639] - v[228] * v[3640]
		+ v[6545 + i3212] + v[3650] * v[871] + v[3585] * v[873] + v[3654] * v[875] + v[3601] * v[877] + v[3559] * v[879]
		+ v[3605] * v[881];
	v[3745] = v[282] * v[3654] + v[3747];
	v[3746] = v[3745] + v[3748] + v[4869];
	v[3749] = v[3747] + v[3748];
	v[3750] = v[3328] * v[3585] - v[3303] * v[3586] - v[3306] * v[3644] - v[3307] * v[3645] + v[3327] * v[3646]
		- v[3300] * v[3651] - v[3304] * v[3657] - v[3299] * v[3664] + v[3329] * v[3665] - v[3302] * v[3669] - v[3301] * v[3674]
		- v[3305] * v[3677];
	v[3751] = v[3329] * v[3587] - v[3305] * v[3588] - v[3302] * v[3642] - v[3304] * v[3647] + v[3327] * v[3648]
		- v[3299] * v[3652] - v[3300] * v[3653] + v[3328] * v[3654] - v[3307] * v[3658] - v[3303] * v[3666] - v[3301] * v[3667]
		- v[3306] * v[3676];
	v[3753] = v[173] * v[3589] + v[174] * v[3590] + v[175] * v[3591] - v[229] * v[3638] - v[230] * v[3639] - v[231] * v[3640]
		+ v[6557 + i3212] + v[3583] * v[871] + v[3646] * v[873] + v[3648] * v[875] + v[3557] * v[877] + v[3597] * v[879]
		+ v[3599] * v[881];
	v[4790] = v[302] * v[3738] + v[303] * v[3744] + v[304] * v[3753];
	v[3781] = v[4790] / v[732];
	v[3754] = v[3755] + v[3757];
	v[3756] = v[288] * v[3646] + v[3755];
	v[3758] = v[3756] + v[3757] + v[4874];
	v[3759] = v[3327] * v[3583] - v[3301] * v[3584] - v[3304] * v[3641] - v[3307] * v[3643] - v[3306] * v[3649]
		+ v[3328] * v[3650] - v[3302] * v[3655] + v[3329] * v[3656] - v[3303] * v[3659] - v[3305] * v[3660] - v[3299] * v[3668]
		- v[3300] * v[3675];
	v[3760] = (v[209] * v[3734] + v[3656] * v[4453] + v[2073] * v[4789]) / v[279];
	v[3761] = v[3771] + v[3773];
	v[4936] = v[3761] / v[279];
	v[3762] = v[3575] + v[3679];
	v[3767] = v[3288] * v[3762];
	v[3763] = -v[3575] + v[3679];
	v[3764] = v[3775] + v[3777];
	v[4932] = v[3764] / v[279];
	v[3765] = v[3287] * v[3661] + v[3767] + v[3769];
	v[3768] = v[3765] - v[3769];
	v[3770] = v[3765] - v[3767];
	v[4933] = v[3770] / v[279];
	v[3772] = v[3287] * v[3731] + v[3771];
	v[3774] = v[3772] + v[3773];
	v[3776] = v[3287] * v[3763] + v[3775];
	v[3778] = v[3776] + v[3777];
	v[3779] = -(v[1014] * v[3230] * v[4790]);
	v[3780] = v[315] * v[3781];
	v[3783] = v[3738];
	v[3987] = v[3783];
	v[3784] = v[316] * v[3738] + v[302] * v[3780];
	v[3785] = v[3744];
	v[3986] = v[3785];
	v[3786] = v[316] * v[3744] + v[303] * v[3780];
	v[3787] = v[3753];
	v[3985] = v[3787];
	v[3788] = v[316] * v[3753] + v[304] * v[3780];
	b3789 = b329;
	if (b3789) {
		v[3790] = -v[3784];
		v[3791] = -v[3786];
		v[3792] = -v[3788];
	}
	else {
		v[3790] = v[3784];
		v[3791] = v[3786];
		v[3792] = v[3788];
	};
	v[4809] = v[333] * v[3790];
	v[4807] = v[3790] * v[470] + v[4812];
	v[4806] = v[332] * v[3790];
	v[3922] = v[3277] * v[3790];
	v[3920] = v[3278] * v[3790];
	v[3918] = v[3279] * v[3790];
	v[3916] = v[3280] * v[3790];
	v[3914] = v[3281] * v[3790];
	v[3912] = v[3282] * v[3790];
	v[4799] = v[333] * v[3791];
	v[4798] = v[331] * v[3791];
	v[3870] = v[3255] * v[3791];
	v[3868] = v[3256] * v[3791];
	v[3866] = v[3257] * v[3791];
	v[3864] = v[3258] * v[3791];
	v[3862] = v[3259] * v[3791];
	v[3860] = v[3260] * v[3791];
	v[4792] = v[332] * v[3792];
	v[4791] = v[331] * v[3792];
	v[3805] = v[3792] * v[4606];
	v[3806] = v[3792] * v[5043];
	v[3807] = v[3792] * v[5044];
	v[3808] = v[3792] * v[5045];
	v[3809] = v[3792] * v[5046];
	v[3810] = v[3792] * v[5047];
	v[3811] = v[3792] * v[5048];
	v[3812] = v[3238] * v[3792];
	v[3813] = v[3237] * v[3792];
	v[3814] = v[3236] * v[3792];
	v[3815] = v[3235] * v[3792];
	v[3816] = v[3234] * v[3792];
	v[3817] = v[3233] * v[3792];
	v[3818] = v[3792] * v[4793];
	v[3825] = v[3792] * v[4794];
	v[3827] = 2e0*v[3792] * v[3826];
	v[3828] = v[3792] * v[5049];
	v[3829] = v[3792] * v[5050];
	v[3830] = v[3792] * v[5051];
	v[3831] = v[3792] * v[5052];
	v[3832] = v[3792] * v[5053];
	v[3833] = v[3792] * v[5054];
	v[3840] = v[3791] * v[4607];
	v[3841] = v[476] * v[4792] + v[3791] * v[5055];
	v[3842] = v[477] * v[4792] + v[3791] * v[5056];
	v[3843] = v[478] * v[4792] + v[3791] * v[5057];
	v[3844] = v[3791] * v[5058] - v[4792] * v[575];
	v[3845] = v[3791] * v[5059] - v[4792] * v[576];
	v[3846] = v[3791] * v[5060] - v[4792] * v[577];
	v[3853] = v[3812] + v[3860];
	v[3854] = v[3813] + v[3862];
	v[3855] = v[3814] + v[3864];
	v[3856] = v[3815] + v[3866];
	v[3857] = v[3816] + v[3868];
	v[3858] = v[3817] + v[3870];
	v[3818] = v[3818] + v[3791] * v[4800];
	v[3825] = v[3825] + 2e0*v[3791] * v[3859];
	v[3827] = v[3827] + v[3791] * v[4801];
	v[3873] = v[3791] * v[5061];
	v[3875] = v[3791] * v[5062];
	v[3877] = v[3791] * v[5063];
	v[3879] = v[3791] * v[5064];
	v[3881] = v[3791] * v[5065];
	v[3883] = v[3791] * v[5066];
	v[3884] = v[3790] * v[4608];
	v[3885] = v[476] * v[4791] + v[473] * v[4798] + v[331] * v[4810] + v[3790] * v[5067] + v[3635] * v[718];
	v[3886] = v[477] * v[4791] + v[474] * v[4798] + v[3790] * v[5068];
	v[3887] = v[478] * v[4791] + v[475] * v[4798] + v[3790] * v[5069];
	v[3888] = v[3790] * v[5070] - v[4798] * v[572] - v[4791] * v[575];
	v[3889] = v[3790] * v[5071] - v[4798] * v[573] - v[4791] * v[576];
	v[3890] = v[3790] * v[5072] - v[4798] * v[574] - v[4791] * v[577];
	v[3818] = v[3818] + 2e0*v[3790] * v[3905] + v[3282] * v[4810] + v[3635] * v[4811];
	v[3906] = v[3812] + v[3912];
	v[3907] = v[3813] + v[3914];
	v[3908] = v[3814] + v[3916];
	v[3909] = v[3815] + v[3918];
	v[3910] = v[3816] + v[3920];
	v[3911] = v[3817] + v[3922];
	v[3825] = v[3825] + v[3634] * v[4802] + v[3260] * (v[4805] + v[4812]) + v[3790] * v[4813];
	v[3913] = v[3860] + v[3912];
	v[3915] = v[3862] + v[3914];
	v[3917] = v[3864] + v[3916];
	v[3919] = v[3866] + v[3918];
	v[3921] = v[3868] + v[3920];
	v[3923] = v[3870] + v[3922];
	v[3827] = v[3827] + v[3614] * v[4796] + v[3238] * (v[4808] + v[4812]) + v[3790] * v[4814];
	v[3926] = v[3790] * v[5073];
	v[3929] = v[3790] * v[5074];
	v[3932] = v[3790] * v[5075];
	v[3935] = v[3790] * v[5076];
	v[3938] = v[3790] * v[5077];
	v[3939] = v[3790] * v[5078];
	v[3940] = v[2247] * v[3622] + v[2248] * v[3623] + v[2249] * v[3687] + v[262] * v[3689] + v[2250] * v[3690]
		+ v[264] * v[3692] + v[3195] * v[3694] + v[3698] * v[4446] + v[3705] * v[4449] + v[10] * v[7665 + i3212];
	v[3941] = v[2252] * v[3622] + v[2253] * v[3631] + v[256] * v[3689] + v[2254] * v[3691] + v[1785] * v[3693]
		+ v[3198] * v[3699] + v[1784] * v[3707] + v[264] * v[3711] + v[2255] * v[3713] + v[10] * v[7677 + i3212];
	v[3942] = v[2258] * v[3612] + v[2257] * v[3622] + v[2259] * v[3688] + v[256] * v[3692] + v[1787] * v[3697]
		+ v[1786] * v[3702] + v[3200] * v[3709] + v[262] * v[3711] + v[2260] * v[3714] + v[10] * v[7689 + i3212];
	v[3943] = v[2262] * v[3671] + v[2263] * v[3672] + v[2264] * v[3730] + v[288] * v[3732] + v[2265] * v[3733]
		+ v[290] * v[3735] + v[3202] * v[3739] + v[3745] * v[4454] + v[3754] * v[4457] + v[10] * v[7701 + i3212];
	v[3944] = v[2267] * v[3671] + v[2268] * v[3680] + v[282] * v[3732] + v[2269] * v[3734] + v[1793] * v[3736]
		+ v[3205] * v[3746] + v[1792] * v[3756] + v[290] * v[3760] + v[2270] * v[3762] + v[10] * v[7713 + i3212];
	v[3945] = v[2273] * v[3661] + v[2272] * v[3671] + v[2274] * v[3731] + v[282] * v[3735] + v[1795] * v[3742]
		+ v[1794] * v[3749] + v[3207] * v[3758] + v[288] * v[3760] + v[2275] * v[3763] + v[10] * v[7725 + i3212];
	v[3946] = v[3217] * v[3792];
	v[3948] = v[3217] * v[3791];
	v[3950] = v[3217] * v[3790];
	v[3954] = v[3218] * v[3792];
	v[3964] = v[3219] * v[3792];
	v[3956] = v[3218] * v[3791];
	v[3958] = v[3950] + v[3956];
	v[3960] = v[3218] * v[3790];
	v[3818] = v[3818] + v[3637] * v[4815] + v[3636] * v[4816] - v[3686] * v[4817] - v[3685] * v[4818] - v[3684] * v[4819]
		+ v[332] * (v[3281] * v[3624] + v[3280] * v[3633] - v[3278] * v[3673] - v[3277] * v[3682] - v[3279] * v[3683] + v[4820])
		+ v[333] * (v[3281] * v[3613] + v[3280] * v[3632] - v[3278] * v[3662] - v[3279] * v[3663] - v[3277] * v[3681] + v[4825]
			);
	v[4008] = v[3818];
	v[3825] = v[3825] + v[3633] * v[4803] + v[3624] * v[4804] + v[331] * (v[3258] * v[3636] + v[3259] * v[3637]
		- v[3257] * v[3684] - v[3255] * v[3685] - v[3256] * v[3686] + v[4820]) - v[3673] * v[4821] - v[3682] * v[4822]
		- v[3683] * v[4823] + v[333] * (v[3259] * v[3613] + v[3258] * v[3632] - v[3256] * v[3662] - v[3257] * v[3663]
			- v[3255] * v[3681] + v[4824]);
	v[4006] = v[3825];
	v[3827] = v[3827] + v[3632] * v[4795] + v[3613] * v[4797] + v[332] * (v[3237] * v[3624] + v[3236] * v[3633]
		- v[3234] * v[3673] - v[3233] * v[3682] - v[3235] * v[3683] + v[4824]) + v[331] * (v[3236] * v[3636] + v[3237] * v[3637]
			- v[3235] * v[3684] - v[3233] * v[3685] - v[3234] * v[3686] + v[4825]) - v[3662] * v[4826] - v[3681] * v[4827]
		- v[3663] * v[4828];
	v[4004] = v[3827];
	v[3962] = v[3956] + v[3964];
	v[3965] = v[3950] + v[3964];
	v[3967] = v[3219] * v[3791];
	v[3969] = v[3219] * v[3790];
	b3971 = b733;
	if (b3971) {
		b3972 = b735;
		if (b3972) {
			v[3790] = 0e0;
			v[3791] = 0e0;
			v[3792] = 0e0;
		}
		else {
		};
	}
	else {
	};
	v[3973] = 0e0;
	v[3974] = 0e0;
	v[3975] = 0e0;
	v[3976] = 0e0;
	v[3977] = 0e0;
	v[3978] = 0e0;
	v[3979] = 0e0;
	b3980 = b733;
	if (b3980) {
		v[3981] = 0e0;
		v[3982] = 0e0;
		v[3983] = 0e0;
		b3984 = b752;
		if (b3984) {
			v[3983] = v[3787];
			v[3787] = 0e0;
			v[3982] = v[3785];
			v[3785] = 0e0;
			v[3981] = v[3783];
			v[3783] = 0e0;
		}
		else {
			v[3978] = -v[3985];
			v[3787] = 0e0;
			v[3977] = -v[3986];
			v[3785] = 0e0;
			v[3976] = -v[3987];
			v[3783] = 0e0;
		};
		v[4829] = v[7] * (v[3981] * v[719] + v[3982] * v[721] + v[3983] * v[723]);
		b3988 = b735;
		if (b3988) {
			v[3975] = v[3983] * v[747];
			v[3974] = v[3982] * v[747];
			v[3973] = v[3981] * v[747];
			v[3979] = (v[1278] * v[4829] * v[745] * Power(v[744], v[2531])) / v[1279];
		}
		else {
			v[3975] = v[3983] * v[751];
			v[3974] = v[3982] * v[751];
			v[3973] = v[3981] * v[751];
			v[3779] = v[3779] + (v[1286] * v[4829] * v[750] * Power(v[732], v[2534])) / v[1287];
		};
	}
	else {
	};
	v[4848] = v[331] * v[3973];
	v[4849] = v[332] * v[3974];
	v[4850] = v[333] * v[3975];
	v[4009] = v[3779];
	v[4007] = v[3976];
	v[4005] = v[3977];
	v[4003] = v[3978];
	b3998 = b733;
	if (b3998) {
		v[4832] = v[331] * v[3976];
		v[4831] = v[332] * v[3977];
		v[4830] = v[333] * v[3978];
		b3999 = b735;
		if (b3999) {
			v[3827] = v[3827] + v[3978] * v[4484];
			v[3978] = 0e0;
			v[3825] = v[3825] + v[3977] * v[4484];
			v[3977] = 0e0;
			v[3818] = v[3818] + v[3976] * v[4484];
			v[3976] = 0e0;
			v[3979] = v[3979] + v[31] * v[34] * (-v[4830] - v[4831] - v[4832])*Power(v[744], v[745]);
			v[3779] = v[3779] - v[3979];
		}
		else {
			v[3827] = v[4004] + v[4003] * v[741];
			v[3978] = 0e0;
			v[3825] = v[4006] + v[4005] * v[741];
			v[3977] = 0e0;
			v[3818] = v[4008] + v[4007] * v[741];
			v[3976] = 0e0;
			v[3779] = v[4009] + v[4493] * (v[4830] + v[4831] + v[4832])*Power(v[732], v[750]);
		};
	}
	else {
	};
	v[4010] = v[3219] * v[3945] + v[301] * v[3975];
	v[4856] = v[333] * v[4010];
	v[4011] = v[3219] * v[3944] + v[291] * v[3975];
	v[4855] = v[333] * v[4011];
	v[4012] = v[3219] * v[3943] + v[281] * v[3975];
	v[4854] = v[333] * v[4012];
	v[4013] = v[3219] * v[3942] + v[275] * v[3975];
	v[4853] = v[333] * v[4013];
	v[4014] = v[3219] * v[3941] + v[265] * v[3975];
	v[4852] = v[333] * v[4014];
	v[4015] = v[3219] * v[3940] + v[255] * v[3975];
	v[4851] = v[333] * v[4015];
	v[3818] = v[3818] + v[3239] * v[3975];
	v[3825] = v[3825] + v[3240] * v[3975];
	v[3827] = v[3827] + v[2709] * v[3975];
	v[4029] = v[3218] * v[3945] + v[301] * v[3974];
	v[4867] = v[332] * v[4029];
	v[4030] = v[3218] * v[3944] + v[291] * v[3974];
	v[4865] = v[332] * v[4030];
	v[4031] = v[3218] * v[3943] + v[281] * v[3974];
	v[4863] = v[332] * v[4031];
	v[4032] = v[3218] * v[3942] + v[275] * v[3974];
	v[4861] = v[332] * v[4032];
	v[4033] = v[3218] * v[3941] + v[265] * v[3974];
	v[4859] = v[332] * v[4033];
	v[4034] = v[3218] * v[3940] + v[255] * v[3974];
	v[4857] = v[332] * v[4034];
	v[3818] = v[3818] + v[3261] * v[3974];
	v[3825] = v[3825] + v[2734] * v[3974];
	v[3827] = v[3827] + v[3264] * v[3974];
	v[4048] = v[3217] * v[3945] + v[301] * v[3973];
	v[4868] = v[331] * v[4048];
	v[4049] = v[3217] * v[3944] + v[291] * v[3973];
	v[4866] = v[331] * v[4049];
	v[4050] = v[3217] * v[3943] + v[281] * v[3973];
	v[4864] = v[331] * v[4050];
	v[4051] = v[3217] * v[3942] + v[275] * v[3973];
	v[4862] = v[331] * v[4051];
	v[4052] = v[3217] * v[3941] + v[265] * v[3973];
	v[4860] = v[331] * v[4052];
	v[4053] = v[3217] * v[3940] + v[255] * v[3973];
	v[4858] = v[331] * v[4053];
	v[3818] = v[3818] + v[2759] * v[3973];
	v[3825] = v[3825] + v[3285] * v[3973];
	v[3827] = v[3827] + v[3286] * v[3973];
	v[4055] = v[1398] * v[3973] + v[1342] * v[3974] + v[1294] * v[3975] + v[3217] * (v[3890] + v[331] * (v[4833] + v[4834])
		- v[3685] * v[718]) + v[3218] * (v[3846] + v[332] * (v[4833] + v[4835]) - v[4806] * v[571] - v[3682] * v[720]) + v[3219] *
		(v[3811] + v[333] * (v[4834] + v[4835]) - v[4809] * v[571] - v[4799] * v[574] - v[3681] * v[722]);
	v[4186] = v[4055] * v[4442];
	v[4056] = v[1400] * v[3973] + v[1344] * v[3974] + v[1296] * v[3975] + v[3217] * (v[3889] + v[331] * (v[4836] + v[4837])
		- v[3686] * v[718]) + v[3218] * (v[3845] + v[332] * (v[4836] + v[4838]) - v[4806] * v[570] - v[3673] * v[720]) + v[3219] *
		(v[3810] + v[333] * (v[4837] + v[4838]) - v[4809] * v[570] - v[4799] * v[573] - v[3662] * v[722]);
	v[4876] = v[279] * v[4056];
	v[4190] = v[4056] * v[4441];
	v[4057] = v[1402] * v[3973] + v[1346] * v[3974] + v[1298] * v[3975] + v[3217] * (v[3888] + v[331] * (v[4839] + v[4840])
		- v[3684] * v[718]) + v[3218] * (v[3844] + v[332] * (v[4839] + v[4841]) - v[4806] * v[569] - v[3683] * v[720]) + v[3219] *
		(v[3809] + v[333] * (v[4840] + v[4841]) - v[4809] * v[569] - v[4799] * v[572] - v[3663] * v[722]);
	v[4877] = v[279] * v[4057];
	v[4193] = v[4057] * v[4440];
	v[4058] = v[1404] * v[3973] + v[1348] * v[3974] + v[1300] * v[3975] + v[3217] * (v[3887] + v[331] * (v[4842] + v[4843])
		+ v[3636] * v[718]) + v[3218] * (v[3843] + v[472] * v[4806] + v[332] * (v[4842] + v[4844]) + v[3633] * v[720]) + v[3219] *
		(v[3808] + v[475] * v[4799] + v[472] * v[4809] + v[333] * (v[4843] + v[4844]) + v[3632] * v[722]);
	v[4277] = v[4058] * v[4435];
	v[4059] = v[1406] * v[3973] + v[1350] * v[3974] + v[1302] * v[3975] + v[3217] * (v[3886] + v[331] * (v[4845] + v[4846])
		+ v[3637] * v[718]) + v[3218] * (v[3842] + v[471] * v[4806] + v[332] * (v[4845] + v[4847]) + v[3624] * v[720]) + v[3219] *
		(v[3807] + v[474] * v[4799] + v[471] * v[4809] + v[333] * (v[4846] + v[4847]) + v[3613] * v[722]);
	v[4885] = v[253] * v[4059];
	v[4281] = v[4059] * v[4434];
	v[4060] = v[3217] * v[3885] + v[1408] * v[3973] + v[1352] * v[3974] + v[1304] * v[3975] + v[3218] * (v[3841] + v[332] *
		(v[4805] + v[4807]) + v[3634] * v[720]) + v[3219] * (v[3806] + v[473] * v[4799] + v[333] * (v[4807] + v[4808])
			+ v[3614] * v[722]);
	v[4886] = v[253] * v[4060];
	v[4284] = v[4060] * v[4433];
	v[4422] = v[3219] * v[3805] + v[331] * v[3946] + v[332] * v[3954] + v[333] * (v[3958] + v[4848] + v[4849])
		+ v[3975] * v[722];
	v[4421] = v[3218] * v[3840] + v[331] * v[3948] + v[333] * v[3967] + v[332] * (v[3965] + v[4848] + v[4850])
		+ v[3974] * v[720];
	v[4420] = v[3217] * v[3884] + v[332] * v[3960] + v[333] * v[3969] + v[331] * (v[3962] + v[4849] + v[4850])
		+ v[3973] * v[718];
	v[8086] = v[4420];
	v[8087] = v[4421];
	v[8088] = v[4422];
	v[8089] = v[3291] * v[3689] + v[3290] * v[3692] + v[3195] * v[4060] + v[3622] * (v[3186] * v[3290] + v[3188] * v[3291]
		+ v[4287]) + v[4059] * v[4446] + v[4058] * v[4449] + v[3599] * v[4779] + v[3605] * v[4780] + v[3725] * v[4881]
		+ v[3561] * v[4919] + v[155] * v[4920] + v[160] * v[4921];
	v[8090] = v[3292] * v[3689] + v[3290] * v[3711] + v[1784] * v[4058] + v[3198] * v[4059] + v[1785] * v[4060] + v[3622] *
		(v[3185] * v[3290] + v[3189] * v[3292] + v[4280]) + v[3616] * v[4781] + v[3729] * v[4880] + v[3559] * v[4922]
		+ v[3597] * v[4923] + v[152] * v[4924] + v[3719] * v[4925];
	v[8091] = v[3292] * v[3692] + v[3291] * v[3711] + v[3200] * v[4058] + v[1786] * v[4059] + v[1787] * v[4060] + v[3622] *
		(v[3187] * v[3291] + v[3190] * v[3292] + v[4272]) + v[3716] * v[4884] + v[3557] * v[4926] + v[3601] * v[4927]
		+ v[3607] * v[4928] + v[3723] * v[4929] + v[3727] * v[4930];
	v[8092] = -v[4420];
	v[8093] = -v[4421];
	v[8094] = -v[4422];
	v[8095] = v[3288] * v[3732] + v[3287] * v[3735] + v[3202] * v[4057] + v[3671] * (v[3157] * v[3287] + v[3159] * v[3288]
		+ v[4196]) + v[4056] * v[4454] + v[4055] * v[4457] + v[3648] * v[4773] + v[3654] * v[4774] + v[3774] * v[4872]
		+ v[3587] * v[4931] + v[211] * v[4932] + v[216] * v[4933];
	v[8096] = v[3289] * v[3732] + v[3287] * v[3760] + v[1792] * v[4055] + v[3205] * v[4056] + v[1793] * v[4057] + v[3671] *
		(v[3156] * v[3287] + v[3160] * v[3289] + v[4189]) + v[3665] * v[4775] + v[3778] * v[4871] + v[3585] * v[4934]
		+ v[3646] * v[4935] + v[208] * v[4936] + v[3768] * v[4937];
	v[8097] = v[3289] * v[3735] + v[3288] * v[3760] + v[3207] * v[4055] + v[1794] * v[4056] + v[1795] * v[4057] + v[3671] *
		(v[3158] * v[3288] + v[3161] * v[3289] + v[4181]) + v[3765] * v[4875] + v[3583] * v[4938] + v[3650] * v[4939]
		+ v[3656] * v[4940] + v[3772] * v[4941] + v[3776] * v[4942];
	v[4067] = v[3929] - v[331] * (v[3858] + v[4856] + v[4867]) - v[4048] * v[718];
	v[4068] = v[3926] - v[331] * (v[3857] + v[4855] + v[4865]) - v[4049] * v[718];
	v[4069] = v[3932] - v[331] * (v[3856] + v[4854] + v[4863]) - v[4050] * v[718];
	v[4070] = v[3938] + v[331] * (v[3855] + v[4853] + v[4861]) + v[4051] * v[718];
	v[4071] = v[3935] + v[331] * (v[3854] + v[4852] + v[4859]) + v[4052] * v[718];
	v[4072] = v[3939] + v[331] * (v[3853] + v[4851] + v[4857]) + v[4053] * v[718];
	v[3818] = v[3818] + v[3241] * v[3946] + v[3262] * v[3948] + v[3283] * v[3962] + v[2332] * v[4048] + v[2334] * v[4049]
		+ v[2336] * v[4050] + v[2338] * v[4051] + v[2340] * v[4052] + v[2342] * v[4053] + v[3853] * v[470] + v[3854] * v[471]
		+ v[3855] * v[472] - v[3856] * v[569] - v[3857] * v[570] - v[3858] * v[571] + v[333] * (v[4015] * v[470] + v[4014] * v[471]
			+ v[4013] * v[472] - v[4012] * v[569] - v[4011] * v[570] - v[4010] * v[571]) + v[332] * (v[4034] * v[470]
				+ v[4033] * v[471] + v[4032] * v[472] - v[4031] * v[569] - v[4030] * v[570] - v[4029] * v[571]) + v[4608] * (-
				(v[3217] * v[3534]) + v[3282] * v[3635] + v[3280] * v[3636] + v[3281] * v[3637] - v[3279] * v[3684] - v[3277] * v[3685]
					- v[3278] * v[3686] + v[3283] * v[3973] + v[4053] * v[470] + v[4052] * v[471] + v[4051] * v[472] - v[4050] * v[569]
					- v[4049] * v[570] - v[4048] * v[571]);
	v[4073] = v[3879] + v[332] * (v[3906] + v[4851] + v[4858]) + v[4034] * v[720];
	v[4074] = v[3883] + v[332] * (v[3907] + v[4852] + v[4860]) + v[4033] * v[720];
	v[4075] = v[3881] + v[332] * (v[3908] + v[4853] + v[4862]) + v[4032] * v[720];
	v[4076] = v[3873] - v[332] * (v[3909] + v[4854] + v[4864]) - v[4031] * v[720];
	v[4077] = v[3877] - v[332] * (v[3910] + v[4855] + v[4866]) - v[4030] * v[720];
	v[4078] = v[3875] - v[332] * (v[3911] + v[4856] + v[4868]) - v[4029] * v[720];
	v[3825] = v[3825] + v[3241] * v[3954] + v[3283] * v[3960] + v[3262] * v[3965] + v[4029] * v[4463] + v[4030] * v[4465]
		+ v[4031] * v[4467] + v[4032] * v[4469] + v[4033] * v[4471] + v[4034] * v[4473] + v[3906] * v[473] + v[3907] * v[474]
		+ v[3908] * v[475] - v[3909] * v[572] - v[3910] * v[573] - v[3911] * v[574] + v[333] * (v[4015] * v[473] + v[4014] * v[474]
			+ v[4013] * v[475] - v[4012] * v[572] - v[4011] * v[573] - v[4010] * v[574]) + v[4607] * (-(v[3218] * v[3535])
				+ v[3259] * v[3624] + v[3258] * v[3633] + v[3260] * v[3634] - v[3256] * v[3673] - v[3255] * v[3682] - v[3257] * v[3683]
				+ v[3262] * v[3974] + v[4034] * v[473] + v[4033] * v[474] + v[4032] * v[475] - v[4031] * v[572] - v[4030] * v[573]
				- v[4029] * v[574]) + v[331] * (v[4053] * v[473] + v[4052] * v[474] + v[4051] * v[475] - v[4050] * v[572]
					- v[4049] * v[573] - v[4048] * v[574]);
	v[4079] = v[3832] + v[333] * (v[3913] + v[4857] + v[4858]) + v[4015] * v[722];
	v[4080] = v[3833] + v[333] * (v[3915] + v[4859] + v[4860]) + v[4014] * v[722];
	v[4081] = v[3831] + v[333] * (v[3917] + v[4861] + v[4862]) + v[4013] * v[722];
	v[4082] = v[3829] - v[333] * (v[3919] + v[4863] + v[4864]) - v[4012] * v[722];
	v[4083] = v[3830] - v[333] * (v[3921] + v[4865] + v[4866]) - v[4011] * v[722];
	v[4084] = v[3828] - v[333] * (v[3923] + v[4867] + v[4868]) - v[4010] * v[722];
	v[3827] = v[3827] + v[3241] * v[3958] + v[3262] * v[3967] + v[3283] * v[3969] + v[2460] * v[4010] + v[2462] * v[4011]
		+ v[2464] * v[4012] + v[2466] * v[4013] + v[2468] * v[4014] + v[2470] * v[4015] + v[3913] * v[476] + v[3915] * v[477]
		+ v[3917] * v[478] - v[3919] * v[575] - v[3921] * v[576] - v[3923] * v[577] + v[4606] * (-(v[3219] * v[3536])
			+ v[3237] * v[3613] + v[3238] * v[3614] + v[3236] * v[3632] - v[3234] * v[3662] - v[3235] * v[3663] - v[3233] * v[3681]
			+ v[3241] * v[3975] + v[4015] * v[476] + v[4014] * v[477] + v[4013] * v[478] - v[4012] * v[575] - v[4011] * v[576]
			- v[4010] * v[577]) + v[332] * (v[4034] * v[476] + v[4033] * v[477] + v[4032] * v[478] - v[4031] * v[575]
				- v[4030] * v[576] - v[4029] * v[577]) + v[331] * (v[4053] * v[476] + v[4052] * v[477] + v[4051] * v[478]
					- v[4050] * v[575] - v[4049] * v[576] - v[4048] * v[577]);
	b4085 = b329;
	if (b4085) {
		v[4086] = -v[3827];
		v[4087] = -v[3825];
		v[4088] = -v[3818];
	}
	else {
		v[4086] = v[3827];
		v[4087] = v[3825];
		v[4088] = v[3818];
	};
	v[3779] = v[3779] + v[314] * v[3325] * v[3781] + v[315] * (v[3320] * v[3738] + v[3319] * v[3744] + v[3318] * v[3753]
		+ v[304] * v[4086] + v[303] * v[4087] + v[302] * v[4088]);
	v[4097] = v[3318] * v[3780] + v[316] * v[4086] + (v[3230] * v[3753] + v[304] * v[3779]) / v[732];
	v[4099] = v[3319] * v[3780] + v[316] * v[4087] + (v[3230] * v[3744] + v[303] * v[3779]) / v[732];
	v[4101] = v[3320] * v[3780] + v[316] * v[4088] + (v[3230] * v[3738] + v[302] * v[3779]) / v[732];
	v[4102] = -(v[3327] * v[3640]) - v[200] * v[4097];
	v[4103] = -(v[3327] * v[3639]) - v[199] * v[4097];
	v[4104] = -(v[3327] * v[3638]) - v[198] * v[4097];
	v[4105] = v[3327] * v[3591] + v[144] * v[4097];
	v[4106] = v[3327] * v[3590] + v[143] * v[4097];
	v[4107] = v[3327] * v[3589] + v[142] * v[4097];
	v[4108] = -(v[3328] * v[3640]) - v[200] * v[4099];
	v[4109] = -(v[3328] * v[3639]) - v[199] * v[4099];
	v[4110] = -(v[3328] * v[3638]) - v[198] * v[4099];
	v[4111] = v[3328] * v[3591] + v[144] * v[4099];
	v[4112] = v[3328] * v[3590] + v[143] * v[4099];
	v[4113] = v[3328] * v[3589] + v[142] * v[4099];
	v[4114] = -(v[3329] * v[3640]) - v[200] * v[4101];
	v[4115] = -(v[3329] * v[3639]) - v[199] * v[4101];
	v[4116] = -(v[3329] * v[3638]) - v[198] * v[4101];
	v[4117] = v[3329] * v[3591] + v[144] * v[4101];
	v[4118] = v[3329] * v[3590] + v[143] * v[4101];
	v[4119] = v[3329] * v[3589] + v[142] * v[4101];
	v[4120] = v[288] * v[4055] + v[3287] * v[4869];
	v[4121] = v[282] * v[4055] + v[3287] * v[4870];
	v[4123] = v[1554] * v[4055] + v[3287] * (v[3749] / v[279] + v[3671] * v[4122]) + v[4120] * v[4871];
	v[4125] = v[1560] * v[4055] + v[3287] * (v[3742] / v[279] + v[3671] * v[4124]) + v[4121] * v[4872];
	v[4127] = (v[218] * v[4120] + v[216] * v[4121] + v[3287] * (v[3758] + v[4126] * v[4789]) + v[4055] * v[4873]) / v[279];
	v[4128] = v[290] * v[4056] + v[3288] * v[4874];
	v[4129] = v[282] * v[4056] + v[3288] * v[4870];
	v[4131] = v[1549] * v[4056] + v[3288] * (v[3756] / v[279] + v[3671] * v[4130]) + v[4128] * v[4875];
	v[4132] = v[4120] + v[4128];
	v[4133] = v[4123] + v[4131];
	v[4135] = (v[3334] * v[3656] + v[207] * v[4129] + v[209] * v[4132] + v[4382] * v[4789] + v[3288] * (v[3736]
		+ v[4134] * v[4789]) + v[1559] * v[4876]) / v[279];
	v[4137] = (v[214] * v[4128] + v[211] * v[4129] + v[3288] * (v[3746] + v[4136] * v[4789]) + v[1553] * v[4876]) / v[279];
	v[4138] = v[288] * v[4057] + v[3289] * v[4869];
	v[4139] = v[290] * v[4057] + v[3289] * v[4874];
	v[4140] = v[4129] + v[4138];
	v[4142] = (v[208] * v[4138] + v[209] * v[4139] + v[3289] * (v[3739] + v[4141] * v[4789]) + v[1557] * v[4877]) / v[279];
	v[4143] = v[4121] + v[4139];
	v[4145] = (v[3339] * v[3650] + v[213] * v[4138] + v[214] * v[4143] + v[4381] * v[4789] + v[3289] * (v[3745]
		+ v[4144] * v[4789]) + v[1564] * v[4877]) / v[279];
	v[4146] = v[4135] + v[4145];
	v[4148] = (v[3340] * v[3646] + v[219] * v[4139] + v[218] * v[4140] + v[4380] * v[4789] + v[3289] * (v[3754]
		+ v[4147] * v[4789]) + v[1567] * v[4877]) / v[279];
	v[4149] = v[4125] + v[4148];
	v[4150] = v[262] * v[4058] + v[3290] * v[4878];
	v[4151] = v[256] * v[4058] + v[3290] * v[4879];
	v[4153] = v[1578] * v[4058] + v[3290] * (v[3702] / v[253] + v[3622] * v[4152]) + v[4150] * v[4880];
	v[4155] = v[1584] * v[4058] + v[3290] * (v[3697] / v[253] + v[3622] * v[4154]) + v[4151] * v[4881];
	v[4157] = (v[162] * v[4150] + v[160] * v[4151] + v[3290] * (v[3709] + v[4156] * v[4788]) + v[4058] * v[4882]) / v[253];
	v[4158] = v[264] * v[4059] + v[3291] * v[4883];
	v[4159] = v[256] * v[4059] + v[3291] * v[4879];
	v[4161] = v[1573] * v[4059] + v[3291] * (v[3707] / v[253] + v[3622] * v[4160]) + v[4158] * v[4884];
	v[4162] = v[4150] + v[4158];
	v[4163] = v[4153] + v[4161];
	v[4165] = (v[3367] * v[3607] + v[151] * v[4159] + v[153] * v[4162] + v[4407] * v[4788] + v[3291] * (v[3693]
		+ v[4164] * v[4788]) + v[1583] * v[4885]) / v[253];
	v[4167] = (v[158] * v[4158] + v[155] * v[4159] + v[3291] * (v[3699] + v[4166] * v[4788]) + v[1577] * v[4885]) / v[253];
	v[4168] = v[262] * v[4060] + v[3292] * v[4878];
	v[4169] = v[264] * v[4060] + v[3292] * v[4883];
	v[4170] = v[4159] + v[4168];
	v[4172] = (v[152] * v[4168] + v[153] * v[4169] + v[3292] * (v[3694] + v[4171] * v[4788]) + v[1581] * v[4886]) / v[253];
	v[4173] = v[4151] + v[4169];
	v[4175] = (v[3372] * v[3601] + v[157] * v[4168] + v[158] * v[4173] + v[4406] * v[4788] + v[3292] * (v[3698]
		+ v[4174] * v[4788]) + v[1588] * v[4886]) / v[253];
	v[4176] = v[4165] + v[4175];
	v[4178] = (v[3373] * v[3597] + v[163] * v[4169] + v[162] * v[4170] + v[4405] * v[4788] + v[3292] * (v[3705]
		+ v[4177] * v[4788]) + v[1591] * v[4886]) / v[253];
	v[4179] = v[4155] + v[4178];
	v[4180] = v[100] * v[4102] + v[1792] * v[4128] + v[292] * v[4186] + v[4139] * v[4457] + v[290] * (v[3765] / v[279]
		+ v[3671] * v[5092]) + v[4104] * v[98] + v[4103] * v[99];
	v[4183] = (v[3340] * v[3730] + v[288] * v[3768] + v[292] * v[4120] + v[277] * v[4140]) / v[279] + v[289] * v[4190]
		+ v[3671] * v[4377] + v[4181] * v[4887] + v[4182] * v[4887] + v[4104] * v[95] + v[4103] * v[96] + v[4102] * v[97];
	v[4185] = v[3207] * v[4121] + v[277] * v[4193] + v[282] * (v[3671] * v[4888] + v[4933]) + v[4104] * v[92]
		+ v[4103] * v[93] + v[4102] * v[94];
	v[4187] = v[100] * v[4108] + (v[3339] * v[3733] + v[290] * v[3776] + v[283] * v[4128] + v[278] * v[4143]) / v[279]
		+ v[3671] * v[4378] + v[4186] * v[4458] + (v[4188] + v[4189])*v[4889] + v[4110] * v[98] + v[4109] * v[99];
	v[4191] = v[1794] * v[4120] + v[283] * v[4190] + v[4138] * v[4454] + v[288] * (v[3778] / v[279] + v[3671] * v[5093])
		+ v[4110] * v[95] + v[4109] * v[96] + v[4108] * v[97];
	v[4194] = v[3205] * v[4129] + v[278] * v[4193] + v[282] * (v[3671] * v[4890] + v[4932]) + v[4110] * v[92]
		+ v[4109] * v[93] + v[4108] * v[94];
	v[4195] = v[100] * v[4114] + v[3671] * v[4379] + (v[3334] * v[3734] + v[290] * v[3772] + v[280] * v[4139]
		+ v[4132] * v[4453]) / v[279] + v[4186] * v[4456] + (v[4196] + v[4198])*v[4889] + v[4116] * v[98] + v[4115] * v[99];
	v[4197] = v[3202] * v[4138] + v[4190] * v[4453] + v[288] * (v[3671] * v[4891] + v[4936]) + v[4116] * v[95]
		+ v[4115] * v[96] + v[4114] * v[97];
	v[4200] = v[1795] * v[4121] + v[1793] * v[4129] + v[280] * v[4193] + v[282] * (v[3774] / v[279] + v[3671] * v[5094])
		+ v[4116] * v[92] + v[4115] * v[93] + v[4114] * v[94];
	v[4201] = v[3302] * v[3638] + v[198] * v[4068];
	v[4202] = v[3302] * v[3639] + v[199] * v[4068];
	v[4203] = v[3302] * v[3640] + v[200] * v[4068];
	v[4204] = v[4201] * v[92] + v[4202] * v[93] + v[4203] * v[94];
	v[4205] = v[4201] * v[95] + v[4202] * v[96] + v[4203] * v[97];
	v[4206] = v[100] * v[4203] + v[4201] * v[98] + v[4202] * v[99];
	v[4207] = v[3299] * v[3638] + v[198] * v[4067];
	v[4208] = v[3299] * v[3639] + v[199] * v[4067];
	v[4209] = v[3299] * v[3640] + v[200] * v[4067];
	v[4210] = v[4207] * v[92] + v[4208] * v[93] + v[4209] * v[94];
	v[4211] = v[100] * v[4209] + v[4207] * v[98] + v[4208] * v[99];
	v[4212] = v[4207] * v[95] + v[4208] * v[96] + v[4209] * v[97];
	v[4213] = v[3305] * v[3638] + v[198] * v[4069];
	v[4214] = v[3305] * v[3639] + v[199] * v[4069];
	v[4215] = v[3305] * v[3640] + v[200] * v[4069];
	v[4216] = v[4213] * v[95] + v[4214] * v[96] + v[4215] * v[97];
	v[4217] = v[100] * v[4215] + v[4213] * v[98] + v[4214] * v[99];
	v[4218] = v[4213] * v[92] + v[4214] * v[93] + v[4215] * v[94];
	v[4219] = v[3306] * v[3638] + v[198] * v[4076];
	v[4220] = v[3306] * v[3639] + v[199] * v[4076];
	v[4221] = v[3306] * v[3640] + v[200] * v[4076];
	v[4222] = v[4219] * v[92] + v[4220] * v[93] + v[4221] * v[94];
	v[4223] = v[4219] * v[95] + v[4220] * v[96] + v[4221] * v[97];
	v[4224] = v[100] * v[4221] + v[4219] * v[98] + v[4220] * v[99];
	v[4225] = v[3300] * v[3638] + v[198] * v[4078];
	v[4226] = v[3300] * v[3639] + v[199] * v[4078];
	v[4227] = v[3300] * v[3640] + v[200] * v[4078];
	v[4228] = v[100] * v[4227] + v[4225] * v[98] + v[4226] * v[99];
	v[4229] = v[4225] * v[92] + v[4226] * v[93] + v[4227] * v[94];
	v[4230] = v[4225] * v[95] + v[4226] * v[96] + v[4227] * v[97];
	v[4231] = v[3301] * v[3638] + v[198] * v[4084];
	v[4232] = v[3301] * v[3639] + v[199] * v[4084];
	v[4233] = v[3301] * v[3640] + v[200] * v[4084];
	v[4234] = v[4231] * v[95] + v[4232] * v[96] + v[4233] * v[97];
	v[4235] = v[4231] * v[92] + v[4232] * v[93] + v[4233] * v[94];
	v[4236] = v[100] * v[4233] + v[4231] * v[98] + v[4232] * v[99];
	v[4237] = -(v[3336] * v[3579]) - v[3358] * v[3581] + 2e0*v[3335] * v[3582] + v[4216] + v[4222] + v[4228] + v[4234]
		+ v[4137] * v[4596] - v[4146] * v[489] - v[4133] * v[492];
	v[4238] = v[3303] * v[3638] + v[198] * v[4077];
	v[4239] = v[3303] * v[3639] + v[199] * v[4077];
	v[4240] = v[3303] * v[3640] + v[200] * v[4077];
	v[4241] = v[4238] * v[92] + v[4239] * v[93] + v[4240] * v[94];
	v[4242] = v[100] * v[4240] + v[4238] * v[98] + v[4239] * v[99];
	v[4243] = v[4238] * v[95] + v[4239] * v[96] + v[4240] * v[97];
	v[4244] = v[3362] * v[3579] + 2e0*v[3338] * v[3581] - v[3358] * v[3582] - v[4205] - v[4211] - v[4235] - v[4241]
		+ v[4142] * v[4595] - v[4146] * v[486] + v[4149] * v[492];
	v[4245] = -(v[3401] * v[3572]) + v[3397] * v[3574] - v[3402] * v[3575] - v[3351] * v[3580] + v[204] * v[4197]
		- v[4216] * v[479] + v[4205] * v[480] - v[4212] * v[481];
	v[4246] = v[3307] * v[3638] + v[198] * v[4082];
	v[4247] = v[3307] * v[3639] + v[199] * v[4082];
	v[4248] = v[3307] * v[3640] + v[200] * v[4082];
	v[4249] = v[4246] * v[92] + v[4247] * v[93] + v[4248] * v[94];
	v[4250] = v[4246] * v[95] + v[4247] * v[96] + v[4248] * v[97];
	v[4251] = v[100] * v[4248] + v[4246] * v[98] + v[4247] * v[99];
	v[4252] = v[3304] * v[3638] + v[198] * v[4083];
	v[4253] = v[3304] * v[3639] + v[199] * v[4083];
	v[4254] = v[3304] * v[3640] + v[200] * v[4083];
	v[4255] = v[4252] * v[95] + v[4253] * v[96] + v[4254] * v[97];
	v[4256] = v[4252] * v[92] + v[4253] * v[93] + v[4254] * v[94];
	v[4257] = v[100] * v[4254] + v[4252] * v[98] + v[4253] * v[99];
	v[4258] = 2e0*v[3332] * v[3579] + v[3362] * v[3581] - v[3336] * v[3582] - v[4217] - v[4242] - v[4249] - v[4255]
		+ v[4127] * v[4594] - v[4133] * v[486] + v[4149] * v[489];
	v[4259] = -(v[3400] * v[3572]) + v[3398] * v[3574] - v[3403] * v[3575] - v[3349] * v[3580] + v[204] * v[4195]
		- v[4217] * v[479] + v[4206] * v[480] - v[4211] * v[481];
	v[4260] = -(v[3409] * v[3572]) + v[3415] * v[3574] - v[3405] * v[3575] - v[3348] * v[3580] + v[204] * v[4194]
		- v[4222] * v[479] + v[4241] * v[480] - v[4229] * v[481];
	v[4261] = v[4210] + v[4230];
	v[4262] = -(v[3408] * v[3572]) + v[3416] * v[3574] - v[3407] * v[3575] - v[3357] * v[3580] + v[204] * v[4187]
		- v[4224] * v[479] + v[4242] * v[480] - v[4228] * v[481];
	v[4263] = -(v[3412] * v[3572]) + v[3424] * v[3574] - v[3420] * v[3575] - v[3344] * v[3580] + v[204] * v[4185]
		- v[4249] * v[479] + v[4256] * v[480] - v[4235] * v[481];
	v[4264] = -(v[3411] * v[3572]) + v[3423] * v[3574] - v[3421] * v[3575] - v[3361] * v[3580] + v[204] * v[4183]
		- v[4250] * v[479] + v[4255] * v[480] - v[4234] * v[481];
	v[4265] = v[4223] + v[4251];
	v[4266] = v[4204] + v[4257];
	v[4267] = -(v[100] * v[3759]) - v[231] * v[4097] - v[228] * v[4099] - v[225] * v[4101] + v[4069] * v[530]
		+ v[4068] * v[531] + v[4067] * v[532] + v[4076] * v[545] + v[4077] * v[546] + v[4078] * v[547] + v[4082] * v[560]
		+ v[4083] * v[561] + v[4084] * v[562] - v[3751] * v[94] - v[3750] * v[97];
	v[4268] = -(v[230] * v[4097]) - v[227] * v[4099] - v[224] * v[4101] + v[4069] * v[527] + v[4068] * v[528]
		+ v[4067] * v[529] + v[4076] * v[542] + v[4077] * v[543] + v[4078] * v[544] + v[4082] * v[557] + v[4083] * v[558]
		+ v[4084] * v[559] - v[3751] * v[93] - v[3750] * v[96] - v[3759] * v[99];
	v[4269] = -(v[229] * v[4097]) - v[226] * v[4099] - v[223] * v[4101] + v[4069] * v[524] + v[4068] * v[525]
		+ v[4067] * v[526] + v[4076] * v[539] + v[4077] * v[540] + v[4078] * v[541] + v[4082] * v[554] + v[4083] * v[555]
		+ v[4084] * v[556] - v[3751] * v[92] - v[3750] * v[95] - v[3759] * v[98];
	v[4270] = v[4269] * v[80] + v[4268] * v[81] + v[4267] * v[82];
	v[4271] = v[1784] * v[4158] + v[266] * v[4277] + v[4169] * v[4449] + v[264] * (v[3716] / v[253] + v[3622] * v[5095])
		+ v[4107] * v[55] + v[4106] * v[56] + v[4105] * v[57];
	v[4274] = (v[3373] * v[3687] + v[262] * v[3719] + v[266] * v[4150] + v[251] * v[4170]) / v[253] + v[263] * v[4281]
		+ v[3622] * v[4402] + v[4272] * v[4892] + v[4273] * v[4892] + v[4107] * v[52] + v[4106] * v[53] + v[4105] * v[54];
	v[4276] = v[3200] * v[4151] + v[251] * v[4284] + v[4107] * v[49] + v[256] * (v[3622] * v[4893] + v[4921])
		+ v[4106] * v[50] + v[4105] * v[51];
	v[4278] = (v[3372] * v[3690] + v[264] * v[3727] + v[257] * v[4158] + v[252] * v[4173]) / v[253] + v[3622] * v[4403]
		+ v[4277] * v[4450] + (v[4279] + v[4280])*v[4894] + v[4113] * v[55] + v[4112] * v[56] + v[4111] * v[57];
	v[4282] = v[1786] * v[4150] + v[257] * v[4281] + v[4168] * v[4446] + v[262] * (v[3729] / v[253] + v[3622] * v[5096])
		+ v[4113] * v[52] + v[4112] * v[53] + v[4111] * v[54];
	v[4285] = v[3198] * v[4159] + v[252] * v[4284] + v[4113] * v[49] + v[256] * (v[3622] * v[4895] + v[4920])
		+ v[4112] * v[50] + v[4111] * v[51];
	v[4286] = v[3622] * v[4404] + (v[3367] * v[3691] + v[264] * v[3723] + v[254] * v[4169] + v[4162] * v[4445]) / v[253]
		+ v[4277] * v[4448] + (v[4287] + v[4289])*v[4894] + v[4119] * v[55] + v[4118] * v[56] + v[4117] * v[57];
	v[4288] = v[3195] * v[4168] + v[4281] * v[4445] + v[262] * (v[3622] * v[4896] + v[4924]) + v[4119] * v[52]
		+ v[4118] * v[53] + v[4117] * v[54];
	v[4291] = v[1787] * v[4151] + v[1785] * v[4159] + v[254] * v[4284] + v[4119] * v[49] + v[4118] * v[50] + v[256] * (v[3725]
		/ v[253] + v[3622] * v[5097]) + v[4117] * v[51];
	v[4292] = v[3311] * v[3589] + v[142] * v[4071];
	v[4293] = v[3311] * v[3590] + v[143] * v[4071];
	v[4294] = v[3311] * v[3591] + v[144] * v[4071];
	v[4295] = v[4292] * v[49] + v[4293] * v[50] + v[4294] * v[51];
	v[4296] = v[4292] * v[52] + v[4293] * v[53] + v[4294] * v[54];
	v[4297] = v[4292] * v[55] + v[4293] * v[56] + v[4294] * v[57];
	v[4298] = v[3308] * v[3589] + v[142] * v[4070];
	v[4299] = v[3308] * v[3590] + v[143] * v[4070];
	v[4300] = v[3308] * v[3591] + v[144] * v[4070];
	v[4301] = v[4298] * v[49] + v[4299] * v[50] + v[4300] * v[51];
	v[4302] = v[4298] * v[55] + v[4299] * v[56] + v[4300] * v[57];
	v[4303] = v[4298] * v[52] + v[4299] * v[53] + v[4300] * v[54];
	v[4304] = v[3314] * v[3589] + v[142] * v[4072];
	v[4305] = v[3314] * v[3590] + v[143] * v[4072];
	v[4306] = v[3314] * v[3591] + v[144] * v[4072];
	v[4307] = v[4304] * v[52] + v[4305] * v[53] + v[4306] * v[54];
	v[4308] = v[4304] * v[55] + v[4305] * v[56] + v[4306] * v[57];
	v[4309] = v[4304] * v[49] + v[4305] * v[50] + v[4306] * v[51];
	v[4310] = v[3315] * v[3589] + v[142] * v[4073];
	v[4311] = v[3315] * v[3590] + v[143] * v[4073];
	v[4312] = v[3315] * v[3591] + v[144] * v[4073];
	v[4313] = v[4310] * v[49] + v[4311] * v[50] + v[4312] * v[51];
	v[4314] = v[4310] * v[52] + v[4311] * v[53] + v[4312] * v[54];
	v[4315] = v[4310] * v[55] + v[4311] * v[56] + v[4312] * v[57];
	v[4316] = v[3309] * v[3589] + v[142] * v[4075];
	v[4317] = v[3309] * v[3590] + v[143] * v[4075];
	v[4318] = v[3309] * v[3591] + v[144] * v[4075];
	v[4319] = v[4316] * v[55] + v[4317] * v[56] + v[4318] * v[57];
	v[4320] = v[4316] * v[49] + v[4317] * v[50] + v[4318] * v[51];
	v[4321] = v[4316] * v[52] + v[4317] * v[53] + v[4318] * v[54];
	v[4322] = v[3310] * v[3589] + v[142] * v[4081];
	v[4323] = v[3310] * v[3590] + v[143] * v[4081];
	v[4324] = v[3310] * v[3591] + v[144] * v[4081];
	v[4325] = v[4322] * v[52] + v[4323] * v[53] + v[4324] * v[54];
	v[4326] = v[4322] * v[49] + v[4323] * v[50] + v[4324] * v[51];
	v[4327] = v[4322] * v[55] + v[4323] * v[56] + v[4324] * v[57];
	v[4328] = -(v[3369] * v[3553]) - v[3391] * v[3555] + 2e0*v[3368] * v[3556] - v[393] * v[4163] - v[390] * v[4176]
		+ v[4307] + v[4313] + v[4319] + v[4325] + v[4167] * v[4592];
	v[4329] = v[3312] * v[3589] + v[142] * v[4074];
	v[4330] = v[3312] * v[3590] + v[143] * v[4074];
	v[4331] = v[3312] * v[3591] + v[144] * v[4074];
	v[4332] = v[4329] * v[49] + v[4330] * v[50] + v[4331] * v[51];
	v[4333] = v[4329] * v[55] + v[4330] * v[56] + v[4331] * v[57];
	v[4334] = v[4329] * v[52] + v[4330] * v[53] + v[4331] * v[54];
	v[4335] = v[3395] * v[3553] + 2e0*v[3371] * v[3555] - v[3391] * v[3556] - v[387] * v[4176] + v[393] * v[4179] - v[4296]
		- v[4302] - v[4326] - v[4332] + v[4172] * v[4591];
	v[4336] = -(v[3445] * v[3546]) + v[3441] * v[3548] - v[3446] * v[3549] - v[3384] * v[3554] + v[148] * v[4288]
		+ v[381] * v[4296] - v[382] * v[4303] - v[380] * v[4307];
	v[4337] = v[3316] * v[3589] + v[142] * v[4079];
	v[4338] = v[3316] * v[3590] + v[143] * v[4079];
	v[4339] = v[3316] * v[3591] + v[144] * v[4079];
	v[4340] = v[4337] * v[49] + v[4338] * v[50] + v[4339] * v[51];
	v[4341] = v[4337] * v[52] + v[4338] * v[53] + v[4339] * v[54];
	v[4342] = v[4337] * v[55] + v[4338] * v[56] + v[4339] * v[57];
	v[4343] = v[3313] * v[3589] + v[142] * v[4080];
	v[4344] = v[3313] * v[3590] + v[143] * v[4080];
	v[4345] = v[3313] * v[3591] + v[144] * v[4080];
	v[4346] = v[4343] * v[52] + v[4344] * v[53] + v[4345] * v[54];
	v[4347] = v[4343] * v[49] + v[4344] * v[50] + v[4345] * v[51];
	v[4348] = v[4343] * v[55] + v[4344] * v[56] + v[4345] * v[57];
	v[4349] = 2e0*v[3365] * v[3553] + v[3395] * v[3555] - v[3369] * v[3556] - v[387] * v[4163] + v[390] * v[4179] - v[4308]
		- v[4333] - v[4340] - v[4346] + v[4157] * v[4590];
	v[4350] = -(v[3444] * v[3546]) + v[3442] * v[3548] - v[3447] * v[3549] - v[3382] * v[3554] + v[148] * v[4286]
		+ v[381] * v[4297] - v[382] * v[4302] - v[380] * v[4308];
	v[4351] = -(v[3453] * v[3546]) + v[3459] * v[3548] - v[3449] * v[3549] - v[3381] * v[3554] + v[148] * v[4285]
		- v[380] * v[4313] - v[382] * v[4320] + v[381] * v[4332];
	v[4352] = v[4301] + v[4321];
	v[4353] = -(v[3452] * v[3546]) + v[3460] * v[3548] - v[3451] * v[3549] - v[3390] * v[3554] + v[148] * v[4278]
		- v[380] * v[4315] - v[382] * v[4319] + v[381] * v[4333];
	v[4354] = -(v[3456] * v[3546]) + v[3468] * v[3548] - v[3464] * v[3549] - v[3377] * v[3554] + v[148] * v[4276]
		- v[382] * v[4326] - v[380] * v[4340] + v[381] * v[4347];
	v[4355] = -(v[3455] * v[3546]) + v[3467] * v[3548] - v[3465] * v[3549] - v[3394] * v[3554] + v[148] * v[4274]
		- v[382] * v[4325] - v[380] * v[4341] + v[381] * v[4346];
	v[4356] = v[4314] + v[4342];
	v[4357] = v[4295] + v[4348];
	v[4358] = v[175] * v[4097] + v[172] * v[4099] + v[169] * v[4101] + v[4072] * v[431] + v[4071] * v[432] + v[4070] * v[433]
		+ v[4073] * v[446] + v[4074] * v[447] + v[4075] * v[448] + v[4079] * v[461] + v[4080] * v[462] + v[4081] * v[463]
		+ v[3704] * v[51] + v[3703] * v[54] + v[3710] * v[57];
	v[4359] = v[174] * v[4097] + v[171] * v[4099] + v[168] * v[4101] + v[4072] * v[428] + v[4071] * v[429] + v[4070] * v[430]
		+ v[4073] * v[443] + v[4074] * v[444] + v[4075] * v[445] + v[4079] * v[458] + v[4080] * v[459] + v[4081] * v[460]
		+ v[3704] * v[50] + v[3703] * v[53] + v[3710] * v[56];
	v[4360] = v[173] * v[4097] + v[170] * v[4099] + v[167] * v[4101] + v[4072] * v[425] + v[4071] * v[426] + v[4070] * v[427]
		+ v[4073] * v[440] + v[4074] * v[441] + v[4075] * v[442] + v[4079] * v[455] + v[4080] * v[456] + v[4081] * v[457]
		+ v[3704] * v[49] + v[3703] * v[52] + v[3710] * v[55];
	v[4361] = v[39] * v[4358] + v[38] * v[4359] + v[37] * v[4360];
	v[4362] = (-v[4270] + v[4269] * v[86] + v[4268] * v[87] + v[4267] * v[88]) / 2e0;
	v[4363] = (-v[4270] + v[4269] * v[83] + v[4268] * v[84] + v[4267] * v[85]) / 2e0;
	v[4364] = (-(v[3422] * v[3524]) - v[3404] * v[3527] - v[3406] * v[3573] - 2e0*v[4123] + 2e0*v[4131]
		+ v[4217] * v[4561] + v[4216] * v[4562] + v[4224] * v[4563] + v[4222] * v[4564] + v[4250] * v[4565] + v[4249] * v[4566]
		- v[4218] * v[482] - v[3405] * v[4897] - v[3402] * v[4898] - v[3421] * v[4899] - v[3407] * v[4900] - v[3420] * v[4901]
		- v[3403] * v[4902] - v[4223] * v[501] - v[4251] * v[519]) / 2e0;
	v[4365] = (-(v[3399] * v[3572]) + v[3396] * v[3574] - v[3404] * v[3575] - v[3354] * v[3580] + v[204] * v[4200]
		- v[4218] * v[479] + v[4204] * v[480] - v[4210] * v[481]) / 2e0;
	v[4367] = (v[3425] * v[3524] + v[3396] * v[3527] + v[3417] * v[3573] - 2e0*v[4125] + 2e0*v[4148] - v[4206] * v[4561]
		- v[4205] * v[4562] - v[4242] * v[4563] - v[4241] * v[4564] - v[4255] * v[4565] - v[4256] * v[4566] + v[4204] * v[482]
		+ v[3415] * v[4897] + v[3397] * v[4898] + v[3423] * v[4899] + v[3416] * v[4900] + v[3424] * v[4901] + v[3398] * v[4902]
		+ v[4243] * v[501] + v[4257] * v[519]) / 2e0;
	v[4368] = (-(v[3410] * v[3572]) + v[3417] * v[3574] - v[3406] * v[3575] - v[3345] * v[3580] + v[204] * v[4191]
		- v[4223] * v[479] + v[4243] * v[480] - v[4230] * v[481]) / 2e0;
	v[4369] = (-(v[3413] * v[3524]) - v[3399] * v[3527] - v[3410] * v[3573] - 2e0*v[4135] + 2e0*v[4145]
		+ v[4211] * v[4561] + v[4212] * v[4562] + v[4228] * v[4563] + v[4229] * v[4564] + v[4234] * v[4565] + v[4235] * v[4566]
		- v[4210] * v[482] - v[3409] * v[4897] - v[3401] * v[4898] - v[3411] * v[4899] - v[3408] * v[4900] - v[3412] * v[4901]
		- v[3400] * v[4902] - v[4230] * v[501] - v[4236] * v[519]) / 2e0;
	v[4417] = v[1868] * v[4367] - v[4368] + v[1865] * v[4369] + v[3525] * v[4905] * v[5099] - v[4428] * ((v[3341] * v[3524])
		/ 2e0 + (v[3354] * v[3527]) / 2e0 + v[3348] * v[3563] + v[3351] * v[3565] + v[3361] * v[3566] + v[3357] * v[3568]
		+ v[3344] * v[3569] + v[3349] * v[3571] + (v[3345] * v[3573]) / 2e0 + v[4206] - v[4212] - v[4224] + v[4229]
		+ v[375] * v[4237] - v[376] * v[4244] + v[4250] - v[4256] - v[378] * v[4258] - 2e0*v[3580] * v[4372] + v[4180] * v[4537]
		+ v[4191] * v[4538] + v[4200] * v[4539] + v[4261] * v[4718] + v[4265] * v[4719] + v[4266] * v[4720] + v[4197] * v[487]
		+ v[4195] * v[493] + v[4194] * v[497] + v[4187] * v[506] + v[4185] * v[510] + v[4586] * (v[2906] * v[3761]
			+ v[2910] * v[3764] + v[2894] * v[3765] + v[2893] * v[3768] + v[2913] * v[3770] + v[2907] * v[3772] + v[2898] * v[3774]
			+ v[2901] * v[3776] + v[2900] * v[3778] + v[3156] * v[4120] + v[3157] * v[4121] + v[4127] + v[3158] * v[4128]
			+ v[3159] * v[4129] + v[2073] * v[4132] + v[4137] + v[3160] * v[4138] + v[3161] * v[4139] + v[2067] * v[4140] + v[4142]
			+ v[2071] * v[4143] + v[3758] * v[4181] + v[3756] * v[4182] + v[3754] * v[4184] + v[3749] * v[4188] + v[3746] * v[4189]
			+ v[3745] * v[4192] + v[3739] * v[4196] + v[3742] * v[4198] + v[3736] * v[4199] + v[3646] * v[4377] + v[3650] * v[4378]
			+ v[3656] * v[4379] + v[3730] * v[4380] + v[3733] * v[4381] + v[3734] * v[4382] + v[4055] * v[4721] + v[4056] * v[4722]
			+ v[4057] * v[4723] + v[3671] * v[4906] * v[5101]) + v[4183] * v[514] + v[7881 + i3212]) + v[4903] * (
				-4e0*v[3438] * v[3525] + v[4364] * v[4904] + v[7905 + i3212]);
	v[4918] = v[4417] + (v[3413] * v[3572] - v[3425] * v[3574] + v[3422] * v[3575] + v[3341] * v[3580] - v[204] * v[4180]
		+ v[4251] * v[479] - v[4257] * v[480] + v[4236] * v[481]) / 2e0;
	v[4384] = v[4259] + v[4263];
	v[4385] = v[4262] + v[4264];
	v[4386] = v[4245] + v[4260];
	v[4387] = (v[43] * v[4360] - v[4361] + v[4359] * v[44] + v[4358] * v[45]) / 2e0;
	v[4388] = (v[42] * v[4358] + v[41] * v[4359] + v[40] * v[4360] - v[4361]) / 2e0;
	v[4389] = (-(v[3466] * v[3516]) - v[3448] * v[3519] - v[3450] * v[3547] - 2e0*v[4153] + 2e0*v[4161] - v[383] * v[4309]
		- v[402] * v[4314] - v[420] * v[4342] + v[4308] * v[4567] + v[4307] * v[4568] + v[4315] * v[4569] + v[4313] * v[4570]
		+ v[4341] * v[4571] + v[4340] * v[4572] - v[3449] * v[4907] - v[3446] * v[4908] - v[3465] * v[4909] - v[3451] * v[4910]
		- v[3464] * v[4911] - v[3447] * v[4912]) / 2e0;
	v[4390] = (-(v[3443] * v[3546]) + v[3440] * v[3548] - v[3448] * v[3549] - v[3387] * v[3554] + v[148] * v[4291]
		+ v[381] * v[4295] - v[382] * v[4301] - v[380] * v[4309]) / 2e0;
	v[4392] = (v[3469] * v[3516] + v[3440] * v[3519] + v[3461] * v[3547] - 2e0*v[4155] + 2e0*v[4178] + v[383] * v[4295]
		+ v[402] * v[4334] + v[420] * v[4348] - v[4297] * v[4567] - v[4296] * v[4568] - v[4333] * v[4569] - v[4332] * v[4570]
		- v[4346] * v[4571] - v[4347] * v[4572] + v[3459] * v[4907] + v[3441] * v[4908] + v[3467] * v[4909] + v[3460] * v[4910]
		+ v[3468] * v[4911] + v[3442] * v[4912]) / 2e0;
	v[4393] = (-(v[3454] * v[3546]) + v[3461] * v[3548] - v[3450] * v[3549] - v[3378] * v[3554] + v[148] * v[4282]
		- v[380] * v[4314] - v[382] * v[4321] + v[381] * v[4334]) / 2e0;
	v[4394] = (-(v[3457] * v[3516]) - v[3443] * v[3519] - v[3454] * v[3547] - 2e0*v[4165] + 2e0*v[4175] - v[383] * v[4301]
		- v[402] * v[4321] - v[420] * v[4327] + v[4302] * v[4567] + v[4303] * v[4568] + v[4319] * v[4569] + v[4320] * v[4570]
		+ v[4325] * v[4571] + v[4326] * v[4572] - v[3453] * v[4907] - v[3445] * v[4908] - v[3455] * v[4909] - v[3452] * v[4910]
		- v[3456] * v[4911] - v[3444] * v[4912]) / 2e0;
	v[4413] = v[1842] * v[4392] - v[4393] + v[1839] * v[4394] + v[3517] * v[4915] * v[5103] - v[4426] * ((v[3374] * v[3516])
		/ 2e0 + (v[3387] * v[3519]) / 2e0 + v[3381] * v[3537] + v[3384] * v[3539] + v[3394] * v[3540] + v[3390] * v[3542]
		+ v[3377] * v[3543] + v[3382] * v[3545] + (v[3378] * v[3547]) / 2e0 + v[415] * v[4274] + v[411] * v[4276]
		+ v[407] * v[4278] + v[398] * v[4285] + v[394] * v[4286] + v[388] * v[4288] + v[4297] - v[4303] - v[4315] + v[4320]
		+ v[369] * v[4328] - v[370] * v[4335] + v[4341] - v[4347] - v[372] * v[4349] - 2e0*v[3554] * v[4397] + v[4271] * v[4558]
		+ v[4282] * v[4559] + v[4291] * v[4560] + v[4352] * v[4730] + v[4356] * v[4731] + v[4357] * v[4732] + v[4585] *
		(v[2936] * v[3712] + v[2940] * v[3715] + v[2924] * v[3716] + v[2923] * v[3719] + v[2943] * v[3721] + v[2937] * v[3723]
			+ v[2928] * v[3725] + v[2931] * v[3727] + v[2930] * v[3729] + v[3185] * v[4150] + v[3186] * v[4151] + v[4157]
			+ v[3187] * v[4158] + v[3188] * v[4159] + v[2045] * v[4162] + v[4167] + v[3189] * v[4168] + v[3190] * v[4169]
			+ v[2039] * v[4170] + v[4172] + v[2043] * v[4173] + v[3709] * v[4272] + v[3707] * v[4273] + v[3705] * v[4275]
			+ v[3702] * v[4279] + v[3699] * v[4280] + v[3698] * v[4283] + v[3694] * v[4287] + v[3697] * v[4289] + v[3693] * v[4290]
			+ v[3597] * v[4402] + v[3601] * v[4403] + v[3607] * v[4404] + v[3687] * v[4405] + v[3690] * v[4406] + v[3691] * v[4407]
			+ v[4058] * v[4733] + v[4059] * v[4734] + v[4060] * v[4735] + v[3622] * v[4916] * v[5105]) + v[7893 + i3212]) + v[4913] * (
				-4e0*v[3482] * v[3517] + v[4389] * v[4914] + v[7989 + i3212]);
	v[4917] = (v[3457] * v[3546] - v[3469] * v[3548] + v[3466] * v[3549] + v[3374] * v[3554] - v[148] * v[4271]
		+ v[382] * v[4327] + v[380] * v[4342] - v[381] * v[4348]) / 2e0 + v[4413];
	v[4409] = v[4350] + v[4354];
	v[4410] = v[4353] + v[4355];
	v[4411] = v[4336] + v[4351];
	v[8062] = 0e0;
	v[8063] = 0e0;
	v[8064] = 0e0;
	v[8065] = 0e0;
	v[8066] = v[3505];
	v[8067] = v[3503];
	v[8068] = 0e0;
	v[8069] = 0e0;
	v[8070] = 0e0;
	v[8071] = 0e0;
	v[8072] = 0e0;
	v[8073] = 0e0;
	v[8026] = 0e0;
	v[8027] = 0e0;
	v[8028] = 0e0;
	v[8029] = v[3505];
	v[8030] = 0e0;
	v[8031] = v[3504];
	v[8032] = 0e0;
	v[8033] = 0e0;
	v[8034] = 0e0;
	v[8035] = 0e0;
	v[8036] = 0e0;
	v[8037] = 0e0;
	v[8002] = 0e0;
	v[8003] = 0e0;
	v[8004] = 0e0;
	v[8005] = v[3503];
	v[8006] = v[3504];
	v[8007] = 0e0;
	v[8008] = 0e0;
	v[8009] = 0e0;
	v[8010] = 0e0;
	v[8011] = 0e0;
	v[8012] = 0e0;
	v[8013] = 0e0;
	v[7978] = 0e0;
	v[7979] = 0e0;
	v[7980] = 0e0;
	v[7981] = 0e0;
	v[7982] = 0e0;
	v[7983] = 0e0;
	v[7984] = 0e0;
	v[7985] = 0e0;
	v[7986] = 0e0;
	v[7987] = 0e0;
	v[7988] = v[3494];
	v[7989] = v[3492];
	v[7942] = 0e0;
	v[7943] = 0e0;
	v[7944] = 0e0;
	v[7945] = 0e0;
	v[7946] = 0e0;
	v[7947] = 0e0;
	v[7948] = 0e0;
	v[7949] = 0e0;
	v[7950] = 0e0;
	v[7951] = v[3494];
	v[7952] = 0e0;
	v[7953] = v[3493];
	v[7918] = 0e0;
	v[7919] = 0e0;
	v[7920] = 0e0;
	v[7921] = 0e0;
	v[7922] = 0e0;
	v[7923] = 0e0;
	v[7924] = 0e0;
	v[7925] = 0e0;
	v[7926] = 0e0;
	v[7927] = v[3492];
	v[7928] = v[3493];
	v[7929] = 0e0;
	v[8074] = v[4101];
	v[8075] = v[4099];
	v[8076] = v[4097];
	v[8077] = -v[4353] + v[4355] + v[372] * v[4409] + v[369] * v[4411] + 2e0*(v[3541] * v[4412] + v[4389] * v[4426]
		+ v[4356] * v[4430] + v[3517] * (-(v[3497] * v[4425]) + v[3480] * v[4742])) + v[368] * v[4917] + (v[3462] * v[3554]
			- v[148] * v[4335] + v[8061 + i3212]) / 2e0;
	v[8078] = v[4350] - v[4354] + v[372] * v[4410] + v[370] * v[4411] + 2e0*(v[3544] * v[4414] - v[4392] * v[4426]
		+ v[4357] * v[4430] + v[3517] * (v[3499] * v[4425] + v[3481] * v[4742])) + v[371] * (-v[4390] + v[4393] + v[4917]) + (-
		(v[3458] * v[3554]) + v[148] * v[4328] + v[8025 + i3212]) / 2e0;
	v[8079] = -v[4336] + v[4351] + v[370] * v[4409] + v[369] * v[4410] + v[373] * (-v[4390] + v[4413]) + 2e0*
		(v[3538] * v[4415] + v[4394] * v[4426] + v[4352] * v[4430] + v[3517] * (-(v[3501] * v[4425]) + v[3476] * v[4742])) +
		(v[3473] * v[3554] - v[148] * v[4349] + v[8001 + i3212]) / 2e0;
	v[8080] = -v[4101];
	v[8081] = -v[4099];
	v[8082] = -v[4097];
	v[8083] = -v[4262] + v[4264] + v[378] * v[4384] + v[375] * v[4386] + 2e0*(v[3567] * v[4416] + v[4364] * v[4428]
		+ v[4265] * v[4437] + v[3525] * (-(v[3486] * v[4427]) + v[3436] * v[4756])) + v[374] * v[4918] + (v[3418] * v[3580]
			- v[204] * v[4244] + v[7977 + i3212]) / 2e0;
	v[8084] = v[4259] - v[4263] + v[378] * v[4385] + v[376] * v[4386] + 2e0*(v[3570] * v[4418] - v[4367] * v[4428]
		+ v[4266] * v[4437] + v[3525] * (v[3488] * v[4427] + v[3437] * v[4756])) + v[377] * (-v[4365] + v[4368] + v[4918]) + (-
		(v[3414] * v[3580]) + v[204] * v[4237] + v[7941 + i3212]) / 2e0;
	v[8085] = -v[4245] + v[4260] + v[376] * v[4384] + v[375] * v[4385] + v[379] * (-v[4365] + v[4417]) + 2e0*
		(v[3564] * v[4419] + v[4369] * v[4428] + v[4261] * v[4437] + v[3525] * (-(v[3490] * v[4427]) + v[3432] * v[4756])) +
		(v[3429] * v[3580] - v[204] * v[4258] + v[7917 + i3212]) / 2e0;
	Rc[i3212 - 1] += v[3495] * v[3528] + v[3496] * v[3529] + v[3484] * v[3531] + v[3485] * v[3532] + v[7601 + i3212]
		+ v[10] * v[7613 + i3212];
	for (i3510 = 1; i3510 <= 12; i3510++) {
		Kc[i3212 - 1][i3510 - 1] += v[4388] * v[6157 + i3510] + v[4387] * v[6169 + i3510] + v[4363] * v[6181 + i3510]
			+ v[4362] * v[6193 + i3510] + v[8073 + i3510] + v[10] * v[8085 + i3510];
	};/* end for */
};/* end for */
#pragma endregion

	delete[] v;
}

void RigidTriangularFace_RigidTriangularFace::HessianPhase1(Matrix& mHes)
{
	double Hes[4][4];
	v = DBG_NEW double[1000];
#pragma region AceGen
	int i01; int i02;
	v[180] = -x1A[0] / 2e0;
	v[183] = -x1A[1] / 2e0;
	v[186] = -x1A[2] / 2e0;
	v[179] = v[180] + x2A[0] / 2e0;
	v[182] = v[183] + x2A[1] / 2e0;
	v[185] = v[186] + x2A[2] / 2e0;
	v[181] = v[180] + x3A[0] / 2e0;
	v[184] = v[183] + x3A[1] / 2e0;
	v[187] = v[186] + x3A[2] / 2e0;
	v[110] = Power(alphaA[0], 2);
	v[108] = (alphaA[0] * alphaA[1]) / 2e0;
	v[103] = Power(alphaA[1], 2);
	v[115] = (alphaA[1] * alphaA[2]) / 2e0;
	v[113] = (alphaA[0] * alphaA[2]) / 2e0;
	v[104] = Power(alphaA[2], 2);
	v[263] = v[103] + v[104];
	v[195] = -x1B[0] / 2e0;
	v[198] = -x1B[1] / 2e0;
	v[201] = -x1B[2] / 2e0;
	v[194] = v[195] + x2B[0] / 2e0;
	v[197] = v[198] + x2B[1] / 2e0;
	v[200] = v[201] + x2B[2] / 2e0;
	v[196] = v[195] + x3B[0] / 2e0;
	v[199] = v[198] + x3B[1] / 2e0;
	v[202] = v[201] + x3B[2] / 2e0;
	v[152] = Power(alphaB[0], 2);
	v[150] = (alphaB[0] * alphaB[1]) / 2e0;
	v[145] = Power(alphaB[1], 2);
	v[157] = (alphaB[1] * alphaB[2]) / 2e0;
	v[155] = (alphaB[0] * alphaB[2]) / 2e0;
	v[146] = Power(alphaB[2], 2);
	v[264] = v[145] + v[146];
	v[102] = 4e0 / (4e0 + v[110] + v[263]);
	v[105] = 1e0 - (v[102] * v[263]) / 2e0;
	v[106] = v[102] * (-alphaA[2] + v[108]);
	v[107] = v[102] * (alphaA[1] + v[113]);
	v[109] = v[102] * (alphaA[2] + v[108]);
	v[111] = 1e0 - (v[102] * (v[104] + v[110])) / 2e0;
	v[112] = v[102] * (-alphaA[0] + v[115]);
	v[114] = v[102] * (-alphaA[1] + v[113]);
	v[116] = v[102] * (alphaA[0] + v[115]);
	v[117] = 1e0 - (v[102] * (v[103] + v[110])) / 2e0;
	v[121] = QAi[0][0] * v[105] + QAi[1][0] * v[106] + QAi[2][0] * v[107];
	v[122] = QAi[0][1] * v[105] + QAi[1][1] * v[106] + QAi[2][1] * v[107];
	v[123] = QAi[0][2] * v[105] + QAi[1][2] * v[106] + QAi[2][2] * v[107];
	v[189] = v[121] * v[181] + v[122] * v[184] + v[123] * v[187];
	v[265] = 2e0*v[189];
	v[188] = v[121] * v[179] + v[122] * v[182] + v[123] * v[185];
	v[124] = QAi[0][0] * v[109] + QAi[1][0] * v[111] + QAi[2][0] * v[112];
	v[125] = QAi[0][1] * v[109] + QAi[1][1] * v[111] + QAi[2][1] * v[112];
	v[126] = QAi[0][2] * v[109] + QAi[1][2] * v[111] + QAi[2][2] * v[112];
	v[191] = v[124] * v[181] + v[125] * v[184] + v[126] * v[187];
	v[266] = 2e0*v[191];
	v[190] = v[124] * v[179] + v[125] * v[182] + v[126] * v[185];
	v[127] = QAi[0][0] * v[114] + QAi[1][0] * v[116] + QAi[2][0] * v[117];
	v[128] = QAi[0][1] * v[114] + QAi[1][1] * v[116] + QAi[2][1] * v[117];
	v[129] = QAi[0][2] * v[114] + QAi[1][2] * v[116] + QAi[2][2] * v[117];
	v[193] = v[127] * v[181] + v[128] * v[184] + v[129] * v[187];
	v[267] = 2e0*v[193];
	v[192] = v[127] * v[179] + v[128] * v[182] + v[129] * v[185];
	v[144] = 4e0 / (4e0 + v[152] + v[264]);
	v[147] = 1e0 - (v[144] * v[264]) / 2e0;
	v[148] = v[144] * (-alphaB[2] + v[150]);
	v[149] = v[144] * (alphaB[1] + v[155]);
	v[151] = v[144] * (alphaB[2] + v[150]);
	v[153] = 1e0 - (v[144] * (v[146] + v[152])) / 2e0;
	v[154] = v[144] * (-alphaB[0] + v[157]);
	v[156] = v[144] * (-alphaB[1] + v[155]);
	v[158] = v[144] * (alphaB[0] + v[157]);
	v[159] = 1e0 - (v[144] * (v[145] + v[152])) / 2e0;
	v[163] = QBi[0][0] * v[147] + QBi[1][0] * v[148] + QBi[2][0] * v[149];
	v[164] = QBi[0][1] * v[147] + QBi[1][1] * v[148] + QBi[2][1] * v[149];
	v[165] = QBi[0][2] * v[147] + QBi[1][2] * v[148] + QBi[2][2] * v[149];
	v[204] = v[163] * v[196] + v[164] * v[199] + v[165] * v[202];
	v[203] = v[163] * v[194] + v[164] * v[197] + v[165] * v[200];
	v[268] = 2e0*v[203];
	v[166] = QBi[0][0] * v[151] + QBi[1][0] * v[153] + QBi[2][0] * v[154];
	v[167] = QBi[0][1] * v[151] + QBi[1][1] * v[153] + QBi[2][1] * v[154];
	v[168] = QBi[0][2] * v[151] + QBi[1][2] * v[153] + QBi[2][2] * v[154];
	v[206] = v[166] * v[196] + v[167] * v[199] + v[168] * v[202];
	v[205] = v[166] * v[194] + v[167] * v[197] + v[168] * v[200];
	v[269] = 2e0*v[205];
	v[169] = QBi[0][0] * v[156] + QBi[1][0] * v[158] + QBi[2][0] * v[159];
	v[170] = QBi[0][1] * v[156] + QBi[1][1] * v[158] + QBi[2][1] * v[159];
	v[171] = QBi[0][2] * v[156] + QBi[1][2] * v[158] + QBi[2][2] * v[159];
	v[208] = v[169] * v[196] + v[170] * v[199] + v[171] * v[202];
	v[207] = v[169] * v[194] + v[170] * v[197] + v[171] * v[200];
	v[270] = 2e0*v[207];
	Hes[0][0] = 1e0*((v[188] * v[188]) + (v[190] * v[190]) + (v[192] * v[192]));
	Hes[0][1] = 0.5e0*(v[188] * v[265] + v[190] * v[266] + v[192] * v[267]);
	Hes[0][2] = 0.5e0*(-(v[188] * v[268]) - v[190] * v[269] - v[192] * v[270]);
	Hes[0][3] = -1e0*(v[188] * v[204] + v[190] * v[206] + v[192] * v[208]);
	Hes[1][1] = 1e0*((v[189] * v[189]) + (v[191] * v[191]) + (v[193] * v[193]));
	Hes[1][2] = -1e0*(v[189] * v[203] + v[191] * v[205] + v[193] * v[207]);
	Hes[1][3] = 0.5e0*(-(v[204] * v[265]) - v[206] * v[266] - v[208] * v[267]);
	Hes[2][2] = 1e0*((v[203] * v[203]) + (v[205] * v[205]) + (v[207] * v[207]));
	Hes[2][3] = 0.5e0*(v[204] * v[268] + v[206] * v[269] + v[208] * v[270]);
	Hes[3][3] = 1e0*((v[204] * v[204]) + (v[206] * v[206]) + (v[208] * v[208]));
	for (i01 = 1; i01 < 4; i01++) {
		for (i02 = 0; i02 < i01; i02++) {
			Hes[i01][i02] = Hes[i02][i01];
		}
	};
#pragma endregion

	delete[] v;
	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
			mHes(i, j) = Hes[i][j];
}

void RigidTriangularFace_RigidTriangularFace::Report()
{
	if (active)
	{
		if (db.execution_data->print_contact_report)
		{
			if (prev_eligible == false && eligible == true)
			{
				db.myprintf("\n---------- New (active) contact ---------- I1: %d I2: %d \nfaceA: %d curveA: %d pointA: %d\nfaceB: %d curveB: %d pointB: %d\n", index1, index2, faceA->ID, deg_curveA, deg_pointA, faceB->ID, deg_curveB, deg_pointB);
				cd->Plot();
			}
			if (prev_eligible == true && eligible == false)
			{
				db.myprintf("\n---------- Old (unnactive) contact ---------- I1: %d I2: %d \nfaceA: %d curveA: %d pointA: %d\nfaceB: %d curveB: %d pointB: %d\n", index1, index2, faceA->ID, deg_curveA, deg_pointA, faceB->ID, deg_curveB, deg_pointB);
				cd->Plot();
			}
		}
		if (db.execution_data->print_contact_report)
		{
			if (cd->g_n[0] < *gnbb)
				db.myprintf("---------- Barrier activated! ---------- %.1f %c of contact layer is active.\n", 100.0*(1.0 - cd->g_n[0] / (*gnb)), 37);
			//cd->Plot();
			//cd->copy_g_t[0]->print();
			//cd->g_t[0]->print();
		}
	}
}

void RigidTriangularFace_RigidTriangularFace::CompactReport()
{
	db.myprintf("Eligible %d\n", (int)eligible);
	db.myprintf("Deg point A: %d\n", deg_pointA);
	db.myprintf("Deg point B: %d\n", deg_pointB);
	db.myprintf("Deg curve A: %d\n", deg_curveA);
	db.myprintf("Deg curve B: %d\n", deg_curveB);
	db.myprintf("Gap %.6e\n", cd->g_n[0]);
}

void RigidTriangularFace_RigidTriangularFace::PredictorTimeStep(double kin)
{
	//Predictor - time step control
	if (eligible == true)
	{
		//Work estimated to be stored in contact spring
		double work = inter->IntegralElasticForce(*gnb, 0.9*(*gnb));
		//if (kin / (-work) > 1.0)
		{
			//Sets zero to counter impact
			td->steps_count_impact = 0;
			double stif = inter->EvaluateElasticStiffness(cd->g_n[0]);
			//db.myprintf("m: %lf stif %lf\n", *meq,stif);
			td->time_step_impact = (1.0 / td->n_steps_impact)* PI * sqrt((*meq) / stif);
		}
		//else
		//	td->time_step_impact = db.solution[db.current_solution_number - 1]->end_time;	//valor alto
		
	}
	else
	{
		//The configuration is already saved when this funcion is called. 
		//Kinematic data to evaluate the relative speed
		Matrix vOA(3);
		Matrix vOB(3);
		Matrix dalphaA(3);
		Matrix dalphaB(3);
		Matrix P_A(3);
		Matrix P_B(3);

		//Vectors "(P-O)" - global system
		P_A = *GammaA - *ptrx0Ap;
		P_B = *GammaB - *ptrx0Bp;

		for (int i = 0; i < 3; i++)
		{
			vOA(i, 0) = db.nodes[node_A - 1]->copy_vel[i];
			vOB(i, 0) = db.nodes[node_B - 1]->copy_vel[i];
			dalphaA(i, 0) = db.nodes[node_A - 1]->copy_vel[i + 3];
			dalphaB(i, 0) = db.nodes[node_B - 1]->copy_vel[i + 3];
		}

		//Obs: Xi é I3 nesse caso - acaba de salvar a configuração e o alpha avaliado em torno desse valor é nulo.
		Matrix omega_A =  dalphaA;
		Matrix omega_B =  dalphaB;

		Matrix VP_A = vOA + cross(omega_A, P_A);
		Matrix VP_B = vOB + cross(omega_B, P_B);

		double vrel = dot(VP_A - VP_B, *cd->n[0]);
		double deltaS = (cd->g_n[0] - (*gnb));
		if (vrel < 0)
			td->time_step_impact = (deltaS) / (-vrel);
		else
			td->time_step_impact = db.solution[db.current_solution_number - 1]->end_time;	//valor alto
	}
}

void RigidTriangularFace_RigidTriangularFace::PrintAceGenPointers()
{
	
	ptrx0Ai->print();
	ptrx0Bi->print();
	ptrQAi->print();
	ptrQBi->print();
	ptrQ0A->print();
	ptrQ0B->print();

	db.PrintPtr(cAp, 2);
	db.PrintPtr(cBp, 2);
	db.PrintPtr(cAi, 2);
	db.PrintPtr(cBi, 2);

	db.PrintPtr(xAi, 3);
	db.PrintPtr(xBi, 3);
	db.PrintPtr(uA, 3);
	db.PrintPtr(uB, 3);

	db.PrintPtr(alphaA, 3);
	db.PrintPtr(alphaB, 3);
	db.PrintPtr(duiA, 3);
	db.PrintPtr(duiB, 3);

	db.PrintPtr(dalphaiA, 3);
	db.PrintPtr(dalphaiB, 3);
	
	db.PrintPtr(QAi, 3, 3);
	db.PrintPtr(QBi, 3, 3);

	db.PrintPtr(Q0A, 3, 3);
	db.PrintPtr(Q0B, 3, 3);

	db.PrintPtr(invH, 4, 4);

	

	cd->Plot();
}

void RigidTriangularFace_RigidTriangularFace::MountLocalContributionsExplicit(double t)
{
	
	cAp[0] = cd->convective[0][0];
	cAp[1] = cd->convective[0][1];
	cBp[0] = cd->convective[0][2];
	cBp[1] = cd->convective[0][3];

	cAi[0] = cd->copy_convective[0][0];
	cAi[1] = cd->copy_convective[0][1];
	cBi[0] = cd->copy_convective[0][2];
	cBi[1] = cd->copy_convective[0][3];
	
	v = DBG_NEW double[766];

	bool *previouscontact = &prev_eligible;
	
	double v01; double v010; double v011; double v012; double v013; double v014;
	double v015; double v016; double v017; double v018; double v019; double v02;
	double v020; double v021; double v022; double v023; double v024; double v025;
	double v026; double v027; double v03; double v04; double v05; double v06; double v07;
	double v08; double v09;
	int b352, b354;
	v[1] = gti[0];
	v[2] = gti[1];
	v[3] = gti[2];
	v[4] = (*epst);
	v[5] = (*mus);
	v[6] = (*mud);
	v[8] = (*meq);
	v[9] = (*ct);
	v[12] = (*epsn1);
	v[13] = (*gnb);
	v[14] = (*gnbb);
	v[349] = v[13] - v[14];
	v[15] = (*n1);
	v[536] = v[12] * v[15];
	v[364] = -1e0 + v[15];
	v[16] = (*n2);
	v[369] = -1e0 + v[16];
	v[17] = x1A[0];
	v[18] = x1A[1];
	v[19] = x1A[2];
	v[20] = x2A[0];
	v[21] = x2A[1];
	v[22] = x2A[2];
	v[23] = x3A[0];
	v[24] = x3A[1];
	v[25] = x3A[2];
	v[26] = xAi[0];
	v[27] = xAi[1];
	v[28] = xAi[2];
	v[29] = QAi[0][0];
	v[30] = QAi[0][1];
	v[31] = QAi[0][2];
	v[32] = QAi[1][0];
	v[33] = QAi[1][1];
	v[34] = QAi[1][2];
	v[35] = QAi[2][0];
	v[36] = QAi[2][1];
	v[37] = QAi[2][2];
	v[38] = Q0A[0][0];
	v[39] = Q0A[0][1];
	v[40] = Q0A[0][2];
	v[156] = v[31] * v[38] + v[34] * v[39] + v[37] * v[40];
	v[155] = v[30] * v[38] + v[33] * v[39] + v[36] * v[40];
	v[154] = v[29] * v[38] + v[32] * v[39] + v[35] * v[40];
	v[41] = Q0A[1][0];
	v[42] = Q0A[1][1];
	v[43] = Q0A[1][2];
	v[162] = v[31] * v[41] + v[34] * v[42] + v[37] * v[43];
	v[161] = v[30] * v[41] + v[33] * v[42] + v[36] * v[43];
	v[160] = v[29] * v[41] + v[32] * v[42] + v[35] * v[43];
	v[44] = Q0A[2][0];
	v[45] = Q0A[2][1];
	v[46] = Q0A[2][2];
	v[168] = v[31] * v[44] + v[34] * v[45] + v[37] * v[46];
	v[167] = v[30] * v[44] + v[33] * v[45] + v[36] * v[46];
	v[166] = v[29] * v[44] + v[32] * v[45] + v[35] * v[46];
	v[50] = alphaA[0];
	v[142] = (v[50] * v[50]);
	v[51] = alphaA[1];
	v[513] = v[51] / 2e0;
	v[140] = v[50] * v[513];
	v[135] = (v[51] * v[51]);
	v[52] = alphaA[2];
	v[147] = v[513] * v[52];
	v[145] = (v[50] * v[52]) / 2e0;
	v[136] = (v[52] * v[52]);
	v[515] = v[135] + v[136];
	v[56] = dalphaiA[0];
	v[57] = dalphaiA[1];
	v[58] = dalphaiA[2];
	v[119] = -0.5e0*cAp[0];
	v[121] = -0.5e0*cAp[1];
	v[124] = -0.5e0*cAi[0];
	v[126] = -0.5e0*cAi[1];
	v[63] = x1B[0];
	v[64] = x1B[1];
	v[65] = x1B[2];
	v[66] = x2B[0];
	v[67] = x2B[1];
	v[68] = x2B[2];
	v[69] = x3B[0];
	v[70] = x3B[1];
	v[71] = x3B[2];
	v[72] = xBi[0];
	v[73] = xBi[1];
	v[74] = xBi[2];
	v[75] = QBi[0][0];
	v[76] = QBi[0][1];
	v[77] = QBi[0][2];
	v[78] = QBi[1][0];
	v[79] = QBi[1][1];
	v[80] = QBi[1][2];
	v[81] = QBi[2][0];
	v[82] = QBi[2][1];
	v[83] = QBi[2][2];
	v[84] = Q0B[0][0];
	v[85] = Q0B[0][1];
	v[86] = Q0B[0][2];
	v[221] = v[77] * v[84] + v[80] * v[85] + v[83] * v[86];
	v[220] = v[76] * v[84] + v[79] * v[85] + v[82] * v[86];
	v[219] = v[75] * v[84] + v[78] * v[85] + v[81] * v[86];
	v[87] = Q0B[1][0];
	v[88] = Q0B[1][1];
	v[89] = Q0B[1][2];
	v[227] = v[77] * v[87] + v[80] * v[88] + v[83] * v[89];
	v[226] = v[76] * v[87] + v[79] * v[88] + v[82] * v[89];
	v[225] = v[75] * v[87] + v[78] * v[88] + v[81] * v[89];
	v[90] = Q0B[2][0];
	v[91] = Q0B[2][1];
	v[92] = Q0B[2][2];
	v[233] = v[77] * v[90] + v[80] * v[91] + v[83] * v[92];
	v[232] = v[76] * v[90] + v[79] * v[91] + v[82] * v[92];
	v[231] = v[75] * v[90] + v[78] * v[91] + v[81] * v[92];
	v[96] = alphaB[0];
	v[207] = (v[96] * v[96]);
	v[97] = alphaB[1];
	v[514] = v[97] / 2e0;
	v[205] = v[514] * v[96];
	v[200] = (v[97] * v[97]);
	v[98] = alphaB[2];
	v[212] = v[514] * v[98];
	v[210] = (v[96] * v[98]) / 2e0;
	v[201] = (v[98] * v[98]);
	v[517] = v[200] + v[201];
	v[102] = dalphaiB[0];
	v[103] = dalphaiB[1];
	v[104] = dalphaiB[2];
	v[184] = -0.5e0*cBp[0];
	v[186] = -0.5e0*cBp[1];
	v[189] = -0.5e0*cBi[0];
	v[191] = -0.5e0*cBi[1];
	v[118] = v[119] + v[121];
	v[120] = 0.5e0 - v[119];
	v[122] = 0.5e0 - v[121];
	v[123] = v[124] + v[126];
	v[125] = 0.5e0 - v[124];
	v[127] = 0.5e0 - v[126];
	v[128] = v[118] * v[17] + v[120] * v[20] + v[122] * v[23];
	v[129] = v[118] * v[18] + v[120] * v[21] + v[122] * v[24];
	v[130] = v[118] * v[19] + v[120] * v[22] + v[122] * v[25];
	v[131] = v[123] * v[17] + v[125] * v[20] + v[127] * v[23];
	v[132] = v[123] * v[18] + v[125] * v[21] + v[127] * v[24];
	v[133] = v[123] * v[19] + v[125] * v[22] + v[127] * v[25];
	v[134] = 4e0 / (4e0 + v[142] + v[515]);
	v[516] = -0.5e0*v[134];
	v[137] = 1e0 + v[515] * v[516];
	v[138] = v[134] * (v[140] - v[52]);
	v[139] = v[134] * (v[145] + v[51]);
	v[141] = v[134] * (v[140] + v[52]);
	v[143] = 1e0 + (v[136] + v[142])*v[516];
	v[144] = v[134] * (v[147] - v[50]);
	v[146] = v[134] * (v[145] - v[51]);
	v[148] = v[134] * (v[147] + v[50]);
	v[149] = 1e0 - (-v[135] - v[142])*v[516];
	v[150] = -(v[516] * v[52]);
	v[151] = -(v[134] * v[513]);
	v[152] = -(v[50] * v[516]);
	v[153] = v[137] * v[154] + v[141] * v[155] + v[146] * v[156];
	v[157] = v[138] * v[154] + v[143] * v[155] + v[148] * v[156];
	v[158] = v[139] * v[154] + v[144] * v[155] + v[149] * v[156];
	v[159] = v[137] * v[160] + v[141] * v[161] + v[146] * v[162];
	v[163] = v[138] * v[160] + v[143] * v[161] + v[148] * v[162];
	v[164] = v[139] * v[160] + v[144] * v[161] + v[149] * v[162];
	v[165] = v[137] * v[166] + v[141] * v[167] + v[146] * v[168];
	v[169] = v[138] * v[166] + v[143] * v[167] + v[148] * v[168];
	v[170] = v[139] * v[166] + v[144] * v[167] + v[149] * v[168];
	v[183] = v[184] + v[186];
	v[185] = 0.5e0 - v[184];
	v[187] = 0.5e0 - v[186];
	v[188] = v[189] + v[191];
	v[190] = 0.5e0 - v[189];
	v[192] = 0.5e0 - v[191];
	v[193] = v[183] * v[63] + v[185] * v[66] + v[187] * v[69];
	v[194] = v[183] * v[64] + v[185] * v[67] + v[187] * v[70];
	v[195] = v[183] * v[65] + v[185] * v[68] + v[187] * v[71];
	v[196] = v[188] * v[63] + v[190] * v[66] + v[192] * v[69];
	v[197] = v[188] * v[64] + v[190] * v[67] + v[192] * v[70];
	v[198] = v[188] * v[65] + v[190] * v[68] + v[192] * v[71];
	v[199] = 4e0 / (4e0 + v[207] + v[517]);
	v[518] = -0.5e0*v[199];
	v[202] = 1e0 + v[517] * v[518];
	v[203] = v[199] * (v[205] - v[98]);
	v[204] = v[199] * (v[210] + v[97]);
	v[206] = v[199] * (v[205] + v[98]);
	v[208] = 1e0 + (v[201] + v[207])*v[518];
	v[209] = v[199] * (v[212] - v[96]);
	v[211] = v[199] * (v[210] - v[97]);
	v[213] = v[199] * (v[212] + v[96]);
	v[214] = 1e0 - (-v[200] - v[207])*v[518];
	v[215] = -(v[518] * v[98]);
	v[216] = -(v[199] * v[514]);
	v[217] = -(v[518] * v[96]);
	v[218] = v[202] * v[219] + v[206] * v[220] + v[211] * v[221];
	v[222] = v[203] * v[219] + v[208] * v[220] + v[213] * v[221];
	v[223] = v[204] * v[219] + v[209] * v[220] + v[214] * v[221];
	v[224] = v[202] * v[225] + v[206] * v[226] + v[211] * v[227];
	v[228] = v[203] * v[225] + v[208] * v[226] + v[213] * v[227];
	v[229] = v[204] * v[225] + v[209] * v[226] + v[214] * v[227];
	v[230] = v[202] * v[231] + v[206] * v[232] + v[211] * v[233];
	v[234] = v[203] * v[231] + v[208] * v[232] + v[213] * v[233];
	v[235] = v[204] * v[231] + v[209] * v[232] + v[214] * v[233];
	v[525] = uA[0] - uB[0] + v[26] - v[72];
	v[524] = uA[1] - uB[1] + v[27] - v[73];
	v[523] = uA[2] - uB[2] + v[28] - v[74];
	v[248] = v[128] * v[153] + v[129] * v[157] + v[130] * v[158] - v[193] * v[218] - v[194] * v[222] - v[195] * v[223] + v[525];
	v[249] = v[128] * v[159] + v[129] * v[163] + v[130] * v[164] - v[193] * v[224] - v[194] * v[228] - v[195] * v[229] + v[524];
	v[250] = v[128] * v[165] + v[129] * v[169] + v[130] * v[170] - v[193] * v[230] - v[194] * v[234] - v[195] * v[235] + v[523];
	v[351] = sqrt((v[248] * v[248]) + (v[249] * v[249]) + (v[250] * v[250]));
	v[251] = v[131] * v[154] + v[132] * v[155] + v[133] * v[156] - v[196] * v[219] - v[197] * v[220] - v[198] * v[221] + v[26]
		- v[72];
	v[252] = v[131] * v[160] + v[132] * v[161] + v[133] * v[162] - v[196] * v[225] - v[197] * v[226] - v[198] * v[227] + v[27]
		- v[73];
	v[253] = v[131] * v[166] + v[132] * v[167] + v[133] * v[168] - v[196] * v[231] - v[197] * v[232] - v[198] * v[233] + v[28]
		- v[74];
	if (v[351] > 0.1e-7) { v01 = 1e0 / v[351]; v02 = (-(v01 / v[351])); v03 = (2e0*v01) / (v[351] * v[351]); }
	else {
		v01 = (12500000e0 / 3e0)*(24e0 - (-0.1e-7 + v[351])*(0.24e10 - 2e0*(-1e0 + 100000000e0*v[351])*(0.2399999997e10
			- 0.1199999994e18*v[351] - 0.3e17*(v[351] * v[351]))));
		v02 = (-50000000e0 / 3e0)*(0.3599999994e10 - 0.4799999982e18*v[351] + 0.6e25*Power(v[351], 3)
			+ 0.1799999982e26*(v[351] * v[351]));
		v03 = 0.1e17*(799999997e0 - 0.599999994e17*v[351] - 0.3e17*(v[351] * v[351]));
	};
	v[262] = v01;
	v[263] = v[248] * v[262];
	v[264] = v[249] * v[262];
	v[265] = v[250] * v[262];
	v[266] = sqrt((v[251] * v[251]) + (v[252] * v[252]) + (v[253] * v[253]));
	if (v[266] > 0.1e-7) { v04 = 1e0 / v[266]; v05 = (-(v04 / v[266])); v06 = (2e0*v04) / (v[266] * v[266]); }
	else {
		v04 = (12500000e0 / 3e0)*(24e0 - (-0.1e-7 + v[266])*(0.24e10 - 2e0*(-1e0 + 100000000e0*v[266])*(0.2399999997e10
			- 0.1199999994e18*v[266] - 0.3e17*(v[266] * v[266]))));
		v05 = (-50000000e0 / 3e0)*(0.3599999994e10 - 0.4799999982e18*v[266] + 0.6e25*Power(v[266], 3)
			+ 0.1799999982e26*(v[266] * v[266]));
		v06 = 0.1e17*(799999997e0 - 0.599999994e17*v[266] - 0.3e17*(v[266] * v[266]));
	};
	v[271] = v04;
	v[272] = v[251] * v[271];
	v[273] = v[252] * v[271];
	v[274] = v[253] * v[271];
	if (v[263] * v[272] + v[264] * v[273] + v[265] * v[274] < 0e0) {
		v[277] = -v[263];
		v[278] = -v[264];
		v[279] = -v[265];
	}
	else {
		v[277] = v[263];
		v[278] = v[264];
		v[279] = v[265];
	};
	v[534] = (v[279] * v[279]);
	v[346] = 1e0 - v[534];
	v[533] = (v[278] * v[278]);
	v[341] = 1e0 - v[533];
	v[532] = (v[277] * v[277]);
	v[336] = 1e0 - v[532];
	if (sqrt(Power(-(v[273] * v[277]) + v[272] * v[278], 2) + Power(v[274] * v[277] - v[272] * v[279], 2) + Power(-
		(v[274] * v[278]) + v[273] * v[279], 2)) > 0.1e-7) {
		v[282] = -(v[274] * v[278]) + v[273] * v[279];
		v[283] = v[274] * v[277] - v[272] * v[279];
		v[284] = -(v[273] * v[277]) + v[272] * v[278];
		v[285] = sqrt((v[282] * v[282]) + (v[283] * v[283]) + (v[284] * v[284]));
		if (v[285] > 0.1e-7) { v07 = 1e0 / v[285]; v08 = (-(v07 / v[285])); v09 = (2e0*v07) / (v[285] * v[285]); }
		else {
			v07 = (12500000e0 / 3e0)*(24e0 - (-0.1e-7 + v[285])*(0.24e10 - 2e0*(-1e0 + 100000000e0*v[285])*
				(0.2399999997e10 - 0.1199999994e18*v[285] - 0.3e17*(v[285] * v[285]))));
			v08 = (-50000000e0 / 3e0)*(0.3599999994e10 - 0.4799999982e18*v[285] + 0.6e25*Power(v[285], 3)
				+ 0.1799999982e26*(v[285] * v[285]));
			v09 = 0.1e17*(799999997e0 - 0.599999994e17*v[285] - 0.3e17*(v[285] * v[285]));
		};
		v[519] = 2e0*v07*tan(asin(v[285]) / 2e0);
		v[294] = v[282] * v[519];
		v[305] = (v[294] * v[294]);
		v[295] = v[283] * v[519];
		v[520] = v[295] / 2e0;
		v[303] = v[294] * v[520];
		v[298] = (v[295] * v[295]);
		v[296] = v[284] * v[519];
		v[310] = v[296] * v[520];
		v[308] = (v[294] * v[296]) / 2e0;
		v[299] = (v[296] * v[296]);
		v[521] = v[298] + v[299];
		v[297] = 4e0 / (4e0 + v[305] + v[521]);
		v[522] = -0.5e0*v[297];
		v[300] = 1e0 + v[521] * v[522];
		v[301] = v[297] * (-v[296] + v[303]);
		v[302] = v[297] * (v[295] + v[308]);
		v[304] = v[297] * (v[296] + v[303]);
		v[306] = 1e0 + (v[299] + v[305])*v[522];
		v[307] = v[297] * (-v[294] + v[310]);
		v[309] = v[297] * (-v[295] + v[308]);
		v[311] = v[297] * (v[294] + v[310]);
		v[312] = 1e0 - (-v[298] - v[305])*v[522];
	}
	else {
		v[300] = 1e0;
		v[301] = 0e0;
		v[302] = 0e0;
		v[304] = 0e0;
		v[306] = 1e0;
		v[307] = 0e0;
		v[309] = 0e0;
		v[311] = 0e0;
		v[312] = 1e0;
	};
	if ((*previouscontact)) {
		v[317] = v[131] * v[165] + v[132] * v[169] + v[133] * v[170] - v[196] * v[230] - v[197] * v[234] - v[198] * v[235]
			+ v[1] * v[309] + v[2] * v[311] + v[3] * v[312] + v[523];
		v[526] = v[279] * v[317];
		v[316] = v[131] * v[159] + v[132] * v[163] + v[133] * v[164] - v[196] * v[224] - v[197] * v[228] - v[198] * v[229]
			+ v[1] * v[304] + v[2] * v[306] + v[3] * v[307] + v[524];
		v[528] = v[278] * v[316];
		v[315] = v[131] * v[153] + v[132] * v[157] + v[133] * v[158] - v[196] * v[218] - v[197] * v[222] - v[198] * v[223]
			+ v[1] * v[300] + v[2] * v[301] + v[3] * v[302] + v[525];
		v[527] = -(v[277] * v[315]);
		v[314] = v[315] * v[336] - v[277] * (v[526] + v[528]);
		v[318] = v[316] * v[341] + v[278] * (-v[526] + v[527]);
		v[319] = v[317] * v[346] + v[279] * (v[527] - v[528]);
	}
	else {
		v[314] = 0e0;
		v[318] = 0e0;
		v[319] = 0e0;
	};
	v[320] = v[134] * v[56] + v[150] * v[57] + v[151] * v[58];
	v[321] = -(v[150] * v[56]) + v[134] * v[57] + v[152] * v[58];
	v[322] = -(v[151] * v[56]) - v[152] * v[57] + v[134] * v[58];
	v[323] = v[102] * v[199] + v[103] * v[215] + v[104] * v[216];
	v[324] = v[103] * v[199] - v[102] * v[215] + v[104] * v[217];
	v[325] = v[104] * v[199] - v[102] * v[216] - v[103] * v[217];
	v[326] = duiA[0] - duiB[0] + v[130] * (-(v[157] * v[320]) + v[153] * v[321]) + v[129] * (v[158] * v[320] - v[153] * v[322])
		+ v[128] * (-(v[158] * v[321]) + v[157] * v[322]) + v[195] * (v[222] * v[323] - v[218] * v[324]) - v[194] * (v[223] * v[323]
			- v[218] * v[325]) + v[193] * (v[223] * v[324] - v[222] * v[325]);
	v[529] = v[277] * v[326];
	v[327] = duiA[1] - duiB[1] + v[130] * (-(v[163] * v[320]) + v[159] * v[321]) + v[129] * (v[164] * v[320] - v[159] * v[322])
		+ v[128] * (-(v[164] * v[321]) + v[163] * v[322]) + v[195] * (v[228] * v[323] - v[224] * v[324]) - v[194] * (v[229] * v[323]
			- v[224] * v[325]) + v[193] * (v[229] * v[324] - v[228] * v[325]);
	v[530] = v[278] * v[327];
	v[539] = v[279] * (v[529] + v[530]);
	v[337] = v[277] * v[530];
	v[328] = duiA[2] - duiB[2] + v[130] * (-(v[169] * v[320]) + v[165] * v[321]) + v[129] * (v[170] * v[320] - v[165] * v[322])
		+ v[128] * (-(v[170] * v[321]) + v[169] * v[322]) + v[195] * (v[234] * v[323] - v[230] * v[324]) - v[194] * (v[235] * v[323]
			- v[230] * v[325]) + v[193] * (v[235] * v[324] - v[234] * v[325]);
	v[531] = v[279] * v[328];
	v[538] = v[278] * (v[529] + v[531]);
	v[338] = v[277] * v[531];
	v[330] = v[337] + v[338] + v[326] * v[532];
	v[332] = v[327] * v[533] + v[538];
	v[334] = v[328] * v[534] + v[539];
	v[348] = -((v[536] * Power(v[349], v[364])) / (v[16] * Power(v[14], v[369])));
	b352 = v[351] < v[13];
	if (b352) {
		b354 = v[351] > v[14];
		if (b354) {
			v[363] = v[13] - v[351];
			v[535] = -(v[12] * Power(v[363], v[15]));
			v[356] = v[277] * v[535];
			v[358] = v[278] * v[535];
			v[359] = v[279] * v[535];
		}
		else {
			v[360] = -(v[12] * Power(v[349], v[15])) + v[348] * (Power(v[14], v[16]) - Power(v[351], v[16]));
			v[356] = v[277] * v[360];
			v[358] = v[278] * v[360];
			v[359] = v[279] * v[360];
		};
	}
	else {
		v[356] = 0e0;
		v[358] = 0e0;
		v[359] = 0e0;
	};
	if (b352) {
		v[537] = 2e0*(*zetan);
		if (b354) {
			v[366] = v[537] * sqrt(v[536] * v[8] * Power(v[363], v[364]));
			v[365] = v[330] * v[366];
			v[367] = v[332] * v[366];
			v[368] = v[334] * v[366];
		}
		else {
			v[370] = v[537] * sqrt(-(v[16] * v[348] * v[8] * Power(v[351], v[369])));
			v[365] = v[330] * v[370];
			v[367] = v[332] * v[370];
			v[368] = v[334] * v[370];
		};
		if (v[351] < v[13] && v[277] * (v[356] + v[365]) + v[278] * (v[358] + v[367]) + v[279] * (v[359] + v[368]) < 0e0) {
			v[373] = v[365];
			v[374] = v[367];
			v[375] = v[368];
		}
		else {
			v[373] = -v[356];
			v[374] = -v[358];
			v[375] = -v[359];
		};
	}
	else {
		v[373] = 0e0;
		v[374] = 0e0;
		v[375] = 0e0;
	};
	v[376] = v[356] + v[373];
	v[377] = v[358] + v[374];
	v[378] = v[359] + v[375];
	v[379] = v[314] * v[4];
	v[380] = v[318] * v[4];
	v[381] = v[319] * v[4];
	v[385] = v[379] + (v[326] * v[336] - v[337] - v[338])*v[9];
	v[386] = v[380] - (-(v[327] * v[341]) + v[538])*v[9];
	v[387] = v[381] - (-(v[328] * v[346]) + v[539])*v[9];
	if (b352) {
		if ((*stick)) {
			if (sqrt((v[385] * v[385]) + (v[386] * v[386]) + (v[387] * v[387])) <= v[5] * sqrt((v[376] * v[376]) +
				(v[377] * v[377]) + (v[378] * v[378]))) {
				v[392] = v[385];
				v[393] = v[386];
				v[394] = v[387];
				v[395] = 1e0;
			}
			else {
				v[396] = sqrt((v[385] * v[385]) + (v[386] * v[386]) + (v[387] * v[387]));
				if (v[396] > 0.1e-5) { v010 = 1e0 / v[396]; v011 = (-(v010 / v[396])); v012 = (2e0*v010) / (v[396] * v[396]); }
				else {
					v010 = (24000000e0 - (-1e0 + 1000000e0*v[396])*(71999994e0 - 0.71999982e14*v[396] + 0.6e19*Power(v[396], 3)
						+ 0.23999982e20*(v[396] * v[396]))) / 24e0;
					v011 = (-500000e0 / 3e0)*(35999994e0 - 0.47999982e14*v[396] + 0.6e19*Power(v[396], 3) + 0.17999982e20*
						(v[396] * v[396]));
					v012 = 0.1e13*(7999997e0 - 0.5999994e13*v[396] - 0.3e13*(v[396] * v[396]));
				};
				v[540] = v010 * v[6] * sqrt((v[376] * v[376]) + (v[377] * v[377]) + (v[378] * v[378]));
				v[392] = v[385] * v[540];
				v[393] = v[386] * v[540];
				v[394] = v[387] * v[540];
				v[395] = 0e0;
			};
			if (sqrt((v[379] * v[379]) + (v[380] * v[380]) + (v[381] * v[381])) > v[5] * sqrt((v[376] * v[376]) + (v[377] * v[377]
				) + (v[378] * v[378]))) {
				if (v[4] > 0.1e-5) { v013 = 1e0 / v[4]; v014 = (-(v013 / v[4])); v015 = (2e0*v013) / (v[4] * v[4]); }
				else {
					v013 = (24000000e0 - (-1e0 + 1000000e0*v[4])*(71999994e0 - 0.71999982e14*v[4] + 0.6e19*Power(v[4], 3)
						+ 0.23999982e20*(v[4] * v[4]))) / 24e0;
					v014 = (-500000e0 / 3e0)*(35999994e0 - 0.47999982e14*v[4] + 0.6e19*Power(v[4], 3) + 0.17999982e20*
						(v[4] * v[4]));
					v015 = 0.1e13*(7999997e0 - 0.5999994e13*v[4] - 0.3e13*(v[4] * v[4]));
				};
				v[414] = sqrt((v[379] * v[379]) + (v[380] * v[380]) + (v[381] * v[381]));
				if (v[414] > 0.1e-5) { v016 = 1e0 / v[414]; v017 = (-(v016 / v[414])); v018 = (2e0*v016) / (v[414] * v[414]); }
				else {
					v016 = (24000000e0 - (-1e0 + 1000000e0*v[414])*(71999994e0 - 0.71999982e14*v[414] + 0.6e19*Power(v[414], 3)
						+ 0.23999982e20*(v[414] * v[414]))) / 24e0;
					v017 = (-500000e0 / 3e0)*(35999994e0 - 0.47999982e14*v[414] + 0.6e19*Power(v[414], 3) + 0.17999982e20*
						(v[414] * v[414]));
					v018 = 0.1e13*(7999997e0 - 0.5999994e13*v[414] - 0.3e13*(v[414] * v[414]));
				};
				v[421] = -(v013*v016*v[6] * sqrt((v[376] * v[376]) + (v[377] * v[377]) + (v[378] * v[378])));
				v[420] = v[314] + v[379] * v[421];
				v[422] = v[318] + v[380] * v[421];
				v[423] = v[319] + v[381] * v[421];
			}
			else {
				v[420] = 0e0;
				v[422] = 0e0;
				v[423] = 0e0;
			};
		}
		else {
			if (sqrt((v[385] * v[385]) + (v[386] * v[386]) + (v[387] * v[387])) <= v[6] * sqrt((v[376] * v[376]) +
				(v[377] * v[377]) + (v[378] * v[378]))) {
				v[392] = v[385];
				v[393] = v[386];
				v[394] = v[387];
				v[395] = 1e0;
			}
			else {
				v[426] = sqrt((v[385] * v[385]) + (v[386] * v[386]) + (v[387] * v[387]));
				if (v[426] > 0.1e-5) { v019 = 1e0 / v[426]; v020 = (-(v019 / v[426])); v021 = (2e0*v019) / (v[426] * v[426]); }
				else {
					v019 = (24000000e0 - (-1e0 + 1000000e0*v[426])*(71999994e0 - 0.71999982e14*v[426] + 0.6e19*Power(v[426], 3)
						+ 0.23999982e20*(v[426] * v[426]))) / 24e0;
					v020 = (-500000e0 / 3e0)*(35999994e0 - 0.47999982e14*v[426] + 0.6e19*Power(v[426], 3) + 0.17999982e20*
						(v[426] * v[426]));
					v021 = 0.1e13*(7999997e0 - 0.5999994e13*v[426] - 0.3e13*(v[426] * v[426]));
				};
				v[541] = v019 * v[6] * sqrt((v[376] * v[376]) + (v[377] * v[377]) + (v[378] * v[378]));
				v[392] = v[385] * v[541];
				v[393] = v[386] * v[541];
				v[394] = v[387] * v[541];
				v[395] = 0e0;
			};
			if (sqrt((v[379] * v[379]) + (v[380] * v[380]) + (v[381] * v[381])) > v[6] * sqrt((v[376] * v[376]) + (v[377] * v[377]
				) + (v[378] * v[378]))) {
				if (v[4] > 0.1e-5) { v022 = 1e0 / v[4]; v023 = (-(v022 / v[4])); v024 = (2e0*v022) / (v[4] * v[4]); }
				else {
					v022 = (24000000e0 - (-1e0 + 1000000e0*v[4])*(71999994e0 - 0.71999982e14*v[4] + 0.6e19*Power(v[4], 3)
						+ 0.23999982e20*(v[4] * v[4]))) / 24e0;
					v023 = (-500000e0 / 3e0)*(35999994e0 - 0.47999982e14*v[4] + 0.6e19*Power(v[4], 3) + 0.17999982e20*
						(v[4] * v[4]));
					v024 = 0.1e13*(7999997e0 - 0.5999994e13*v[4] - 0.3e13*(v[4] * v[4]));
				};
				v[444] = sqrt((v[379] * v[379]) + (v[380] * v[380]) + (v[381] * v[381]));
				if (v[444] > 0.1e-5) { v025 = 1e0 / v[444]; v026 = (-(v025 / v[444])); v027 = (2e0*v025) / (v[444] * v[444]); }
				else {
					v025 = (24000000e0 - (-1e0 + 1000000e0*v[444])*(71999994e0 - 0.71999982e14*v[444] + 0.6e19*Power(v[444], 3)
						+ 0.23999982e20*(v[444] * v[444]))) / 24e0;
					v026 = (-500000e0 / 3e0)*(35999994e0 - 0.47999982e14*v[444] + 0.6e19*Power(v[444], 3) + 0.17999982e20*
						(v[444] * v[444]));
					v027 = 0.1e13*(7999997e0 - 0.5999994e13*v[444] - 0.3e13*(v[444] * v[444]));
				};
				v[450] = -(v022*v025*v[6] * sqrt((v[376] * v[376]) + (v[377] * v[377]) + (v[378] * v[378])));
				v[420] = v[314] + v[379] * v[450];
				v[422] = v[318] + v[380] * v[450];
				v[423] = v[319] + v[381] * v[450];
			}
			else {
				v[420] = 0e0;
				v[422] = 0e0;
				v[423] = 0e0;
			};
		};
	}
	else {
		v[392] = 0e0;
		v[393] = 0e0;
		v[394] = 0e0;
	};
	v[544] = v[378] + v[394];
	v[543] = v[377] + v[393];
	v[542] = v[376] + v[392];
	fn[0] = v[376];
	fn[1] = v[377];
	fn[2] = v[378];
	ft[0] = v[392];
	ft[1] = v[393];
	ft[2] = v[394];
	(*stickupdated) = v[395];
	gtpupdated[0] = v[314] - v[420];
	gtpupdated[1] = v[318] - v[422];
	gtpupdated[2] = v[319] - v[423];
	Rc[0] = v[153] * v[542] + v[159] * v[543] + v[165] * v[544];
	Rc[1] = v[157] * v[542] + v[163] * v[543] + v[169] * v[544];
	Rc[2] = v[158] * v[542] + v[164] * v[543] + v[170] * v[544];
	Rc[3] = (-(v[130] * v[157]) + v[129] * v[158])*v[542] + (-(v[130] * v[163]) + v[129] * v[164])*v[543] + (-
		(v[130] * v[169]) + v[129] * v[170])*v[544];
	Rc[4] = (v[130] * v[153] - v[128] * v[158])*v[542] + (v[130] * v[159] - v[128] * v[164])*v[543] + (v[130] * v[165]
		- v[128] * v[170])*v[544];
	Rc[5] = (-(v[129] * v[153]) + v[128] * v[157])*v[542] + (-(v[129] * v[159]) + v[128] * v[163])*v[543] + (-
		(v[129] * v[165]) + v[128] * v[169])*v[544];
	Rc[6] = -(v[218] * v[542]) - v[224] * v[543] - v[230] * v[544];
	Rc[7] = -(v[222] * v[542]) - v[228] * v[543] - v[234] * v[544];
	Rc[8] = -(v[223] * v[542]) - v[229] * v[543] - v[235] * v[544];
	Rc[9] = -((-(v[195] * v[222]) + v[194] * v[223])*v[542]) - (-(v[195] * v[228]) + v[194] * v[229])*v[543] - (-
		(v[195] * v[234]) + v[194] * v[235])*v[544];
	Rc[10] = -((v[195] * v[218] - v[193] * v[223])*v[542]) - (v[195] * v[224] - v[193] * v[229])*v[543] - (v[195] * v[230]
		- v[193] * v[235])*v[544];
	Rc[11] = -((-(v[194] * v[218]) + v[193] * v[222])*v[542]) - (-(v[194] * v[224]) + v[193] * v[228])*v[543] - (-
		(v[194] * v[230]) + v[193] * v[234])*v[544];


	delete[]v;
}
void RigidTriangularFace_RigidTriangularFace::SetVariablesExplicit(double t)
{
	for (int i = 0; i < 3; i++)
	{
		xAi[i] = db.nodes[node_A - 1]->copy_coordinates[i];
		xBi[i] = db.nodes[node_B - 1]->copy_coordinates[i];
		
		uA[i] = (*db.nodes[node_A - 1]->u)(i, 0);
		uB[i] = (*db.nodes[node_B - 1]->u)(i, 0);
		
		alphaA[i] = (*db.nodes[node_A - 1]->alpha)(i, 0);
		alphaB[i] = (*db.nodes[node_B - 1]->alpha)(i, 0);

		duiA[i] = (*db.nodes[node_A - 1]->du)(i, 0);
		duiB[i] = (*db.nodes[node_B - 1]->du)(i, 0);

		dalphaiA[i] = (*db.nodes[node_A - 1]->omega)(i, 0);
		dalphaiB[i] = (*db.nodes[node_B - 1]->omega)(i, 0);
		

		for (int j = 0; j < 3; j++)
		{
			QAi[i][j] = (*db.nodes[node_A - 1]->Q)(i, j);
			QBi[i][j] = (*db.nodes[node_B - 1]->Q)(i, j);
		}
		
		for (int j = 0; j < 3; j++)
		{
			Q0A[i][j] = (*ptrQ0A)(i, j);
			Q0B[i][j] = (*ptrQ0B)(i, j);
		}
	}
}

void RigidTriangularFace_RigidTriangularFace::FinalUpdateExplicit(double t)
{
	cAp[0] = cd->convective[0][0];
	cAp[1] = cd->convective[0][1];
	cBp[0] = cd->convective[0][2];
	cBp[1] = cd->convective[0][3];

	cAi[0] = cd->copy_convective[0][0];
	cAi[1] = cd->copy_convective[0][1];
	cBi[0] = cd->copy_convective[0][2];
	cBi[1] = cd->copy_convective[0][3];

	v = DBG_NEW double[704];

	bool* previouscontact = &prev_eligible;

	double v01; double v010; double v011; double v012; double v013; double v014;
	double v015; double v016; double v017; double v018; double v019; double v02;
	double v020; double v021; double v022; double v023; double v024; double v025;
	double v026; double v027; double v03; double v04; double v05; double v06; double v07;
	double v08; double v09;
	int b352, b354;
	v[1] = gti[0];
	v[2] = gti[1];
	v[3] = gti[2];
	v[4] = (*epst);
	v[5] = (*mus);
	v[6] = (*mud);
	v[8] = (*meq);
	v[9] = (*ct);
	v[12] = (*epsn1);
	v[13] = (*gnb);
	v[14] = (*gnbb);
	v[349] = v[13] - v[14];
	v[15] = (*n1);
	v[536] = v[12] * v[15];
	v[364] = -1e0 + v[15];
	v[16] = (*n2);
	v[369] = -1e0 + v[16];
	v[17] = x1A[0];
	v[18] = x1A[1];
	v[19] = x1A[2];
	v[20] = x2A[0];
	v[21] = x2A[1];
	v[22] = x2A[2];
	v[23] = x3A[0];
	v[24] = x3A[1];
	v[25] = x3A[2];
	v[26] = xAi[0];
	v[27] = xAi[1];
	v[28] = xAi[2];
	v[29] = QAi[0][0];
	v[30] = QAi[0][1];
	v[31] = QAi[0][2];
	v[32] = QAi[1][0];
	v[33] = QAi[1][1];
	v[34] = QAi[1][2];
	v[35] = QAi[2][0];
	v[36] = QAi[2][1];
	v[37] = QAi[2][2];
	v[38] = Q0A[0][0];
	v[39] = Q0A[0][1];
	v[40] = Q0A[0][2];
	v[156] = v[31] * v[38] + v[34] * v[39] + v[37] * v[40];
	v[155] = v[30] * v[38] + v[33] * v[39] + v[36] * v[40];
	v[154] = v[29] * v[38] + v[32] * v[39] + v[35] * v[40];
	v[41] = Q0A[1][0];
	v[42] = Q0A[1][1];
	v[43] = Q0A[1][2];
	v[162] = v[31] * v[41] + v[34] * v[42] + v[37] * v[43];
	v[161] = v[30] * v[41] + v[33] * v[42] + v[36] * v[43];
	v[160] = v[29] * v[41] + v[32] * v[42] + v[35] * v[43];
	v[44] = Q0A[2][0];
	v[45] = Q0A[2][1];
	v[46] = Q0A[2][2];
	v[168] = v[31] * v[44] + v[34] * v[45] + v[37] * v[46];
	v[167] = v[30] * v[44] + v[33] * v[45] + v[36] * v[46];
	v[166] = v[29] * v[44] + v[32] * v[45] + v[35] * v[46];
	v[50] = alphaA[0];
	v[142] = (v[50] * v[50]);
	v[51] = alphaA[1];
	v[513] = v[51] / 2e0;
	v[140] = v[50] * v[513];
	v[135] = (v[51] * v[51]);
	v[52] = alphaA[2];
	v[147] = v[513] * v[52];
	v[145] = (v[50] * v[52]) / 2e0;
	v[136] = (v[52] * v[52]);
	v[515] = v[135] + v[136];
	v[56] = dalphaiA[0];
	v[57] = dalphaiA[1];
	v[58] = dalphaiA[2];
	v[119] = -0.5e0 * cAp[0];
	v[121] = -0.5e0 * cAp[1];
	v[124] = -0.5e0 * cAi[0];
	v[126] = -0.5e0 * cAi[1];
	v[63] = x1B[0];
	v[64] = x1B[1];
	v[65] = x1B[2];
	v[66] = x2B[0];
	v[67] = x2B[1];
	v[68] = x2B[2];
	v[69] = x3B[0];
	v[70] = x3B[1];
	v[71] = x3B[2];
	v[72] = xBi[0];
	v[73] = xBi[1];
	v[74] = xBi[2];
	v[75] = QBi[0][0];
	v[76] = QBi[0][1];
	v[77] = QBi[0][2];
	v[78] = QBi[1][0];
	v[79] = QBi[1][1];
	v[80] = QBi[1][2];
	v[81] = QBi[2][0];
	v[82] = QBi[2][1];
	v[83] = QBi[2][2];
	v[84] = Q0B[0][0];
	v[85] = Q0B[0][1];
	v[86] = Q0B[0][2];
	v[221] = v[77] * v[84] + v[80] * v[85] + v[83] * v[86];
	v[220] = v[76] * v[84] + v[79] * v[85] + v[82] * v[86];
	v[219] = v[75] * v[84] + v[78] * v[85] + v[81] * v[86];
	v[87] = Q0B[1][0];
	v[88] = Q0B[1][1];
	v[89] = Q0B[1][2];
	v[227] = v[77] * v[87] + v[80] * v[88] + v[83] * v[89];
	v[226] = v[76] * v[87] + v[79] * v[88] + v[82] * v[89];
	v[225] = v[75] * v[87] + v[78] * v[88] + v[81] * v[89];
	v[90] = Q0B[2][0];
	v[91] = Q0B[2][1];
	v[92] = Q0B[2][2];
	v[233] = v[77] * v[90] + v[80] * v[91] + v[83] * v[92];
	v[232] = v[76] * v[90] + v[79] * v[91] + v[82] * v[92];
	v[231] = v[75] * v[90] + v[78] * v[91] + v[81] * v[92];
	v[96] = alphaB[0];
	v[207] = (v[96] * v[96]);
	v[97] = alphaB[1];
	v[514] = v[97] / 2e0;
	v[205] = v[514] * v[96];
	v[200] = (v[97] * v[97]);
	v[98] = alphaB[2];
	v[212] = v[514] * v[98];
	v[210] = (v[96] * v[98]) / 2e0;
	v[201] = (v[98] * v[98]);
	v[517] = v[200] + v[201];
	v[102] = dalphaiB[0];
	v[103] = dalphaiB[1];
	v[104] = dalphaiB[2];
	v[184] = -0.5e0 * cBp[0];
	v[186] = -0.5e0 * cBp[1];
	v[189] = -0.5e0 * cBi[0];
	v[191] = -0.5e0 * cBi[1];
	v[118] = v[119] + v[121];
	v[120] = 0.5e0 - v[119];
	v[122] = 0.5e0 - v[121];
	v[123] = v[124] + v[126];
	v[125] = 0.5e0 - v[124];
	v[127] = 0.5e0 - v[126];
	v[128] = v[118] * v[17] + v[120] * v[20] + v[122] * v[23];
	v[129] = v[118] * v[18] + v[120] * v[21] + v[122] * v[24];
	v[130] = v[118] * v[19] + v[120] * v[22] + v[122] * v[25];
	v[131] = v[123] * v[17] + v[125] * v[20] + v[127] * v[23];
	v[132] = v[123] * v[18] + v[125] * v[21] + v[127] * v[24];
	v[133] = v[123] * v[19] + v[125] * v[22] + v[127] * v[25];
	v[134] = 4e0 / (4e0 + v[142] + v[515]);
	v[516] = -0.5e0 * v[134];
	v[137] = 1e0 + v[515] * v[516];
	v[138] = v[134] * (v[140] - v[52]);
	v[139] = v[134] * (v[145] + v[51]);
	v[141] = v[134] * (v[140] + v[52]);
	v[143] = 1e0 + (v[136] + v[142]) * v[516];
	v[144] = v[134] * (v[147] - v[50]);
	v[146] = v[134] * (v[145] - v[51]);
	v[148] = v[134] * (v[147] + v[50]);
	v[149] = 1e0 - (-v[135] - v[142]) * v[516];
	v[150] = -(v[516] * v[52]);
	v[151] = -(v[134] * v[513]);
	v[152] = -(v[50] * v[516]);
	v[153] = v[137] * v[154] + v[141] * v[155] + v[146] * v[156];
	v[157] = v[138] * v[154] + v[143] * v[155] + v[148] * v[156];
	v[158] = v[139] * v[154] + v[144] * v[155] + v[149] * v[156];
	v[159] = v[137] * v[160] + v[141] * v[161] + v[146] * v[162];
	v[163] = v[138] * v[160] + v[143] * v[161] + v[148] * v[162];
	v[164] = v[139] * v[160] + v[144] * v[161] + v[149] * v[162];
	v[165] = v[137] * v[166] + v[141] * v[167] + v[146] * v[168];
	v[169] = v[138] * v[166] + v[143] * v[167] + v[148] * v[168];
	v[170] = v[139] * v[166] + v[144] * v[167] + v[149] * v[168];
	v[183] = v[184] + v[186];
	v[185] = 0.5e0 - v[184];
	v[187] = 0.5e0 - v[186];
	v[188] = v[189] + v[191];
	v[190] = 0.5e0 - v[189];
	v[192] = 0.5e0 - v[191];
	v[193] = v[183] * v[63] + v[185] * v[66] + v[187] * v[69];
	v[194] = v[183] * v[64] + v[185] * v[67] + v[187] * v[70];
	v[195] = v[183] * v[65] + v[185] * v[68] + v[187] * v[71];
	v[196] = v[188] * v[63] + v[190] * v[66] + v[192] * v[69];
	v[197] = v[188] * v[64] + v[190] * v[67] + v[192] * v[70];
	v[198] = v[188] * v[65] + v[190] * v[68] + v[192] * v[71];
	v[199] = 4e0 / (4e0 + v[207] + v[517]);
	v[518] = -0.5e0 * v[199];
	v[202] = 1e0 + v[517] * v[518];
	v[203] = v[199] * (v[205] - v[98]);
	v[204] = v[199] * (v[210] + v[97]);
	v[206] = v[199] * (v[205] + v[98]);
	v[208] = 1e0 + (v[201] + v[207]) * v[518];
	v[209] = v[199] * (v[212] - v[96]);
	v[211] = v[199] * (v[210] - v[97]);
	v[213] = v[199] * (v[212] + v[96]);
	v[214] = 1e0 - (-v[200] - v[207]) * v[518];
	v[215] = -(v[518] * v[98]);
	v[216] = -(v[199] * v[514]);
	v[217] = -(v[518] * v[96]);
	v[218] = v[202] * v[219] + v[206] * v[220] + v[211] * v[221];
	v[222] = v[203] * v[219] + v[208] * v[220] + v[213] * v[221];
	v[223] = v[204] * v[219] + v[209] * v[220] + v[214] * v[221];
	v[224] = v[202] * v[225] + v[206] * v[226] + v[211] * v[227];
	v[228] = v[203] * v[225] + v[208] * v[226] + v[213] * v[227];
	v[229] = v[204] * v[225] + v[209] * v[226] + v[214] * v[227];
	v[230] = v[202] * v[231] + v[206] * v[232] + v[211] * v[233];
	v[234] = v[203] * v[231] + v[208] * v[232] + v[213] * v[233];
	v[235] = v[204] * v[231] + v[209] * v[232] + v[214] * v[233];
	v[525] = uA[0] - uB[0] + v[26] - v[72];
	v[524] = uA[1] - uB[1] + v[27] - v[73];
	v[523] = uA[2] - uB[2] + v[28] - v[74];
	v[248] = v[128] * v[153] + v[129] * v[157] + v[130] * v[158] - v[193] * v[218] - v[194] * v[222] - v[195] * v[223] + v[525];
	v[249] = v[128] * v[159] + v[129] * v[163] + v[130] * v[164] - v[193] * v[224] - v[194] * v[228] - v[195] * v[229] + v[524];
	v[250] = v[128] * v[165] + v[129] * v[169] + v[130] * v[170] - v[193] * v[230] - v[194] * v[234] - v[195] * v[235] + v[523];
	v[351] = sqrt((v[248] * v[248]) + (v[249] * v[249]) + (v[250] * v[250]));
	v[251] = v[131] * v[154] + v[132] * v[155] + v[133] * v[156] - v[196] * v[219] - v[197] * v[220] - v[198] * v[221] + v[26]
		- v[72];
	v[252] = v[131] * v[160] + v[132] * v[161] + v[133] * v[162] - v[196] * v[225] - v[197] * v[226] - v[198] * v[227] + v[27]
		- v[73];
	v[253] = v[131] * v[166] + v[132] * v[167] + v[133] * v[168] - v[196] * v[231] - v[197] * v[232] - v[198] * v[233] + v[28]
		- v[74];
	if (v[351] > 0.1e-7) { v01 = 1e0 / v[351]; v02 = (-(v01 / v[351])); v03 = (2e0 * v01) / (v[351] * v[351]); }
	else {
		v01 = (12500000e0 / 3e0) * (24e0 - (-0.1e-7 + v[351]) * (0.24e10 - 2e0 * (-1e0 + 100000000e0 * v[351]) * (0.2399999997e10
			- 0.1199999994e18 * v[351] - 0.3e17 * (v[351] * v[351]))));
		v02 = (-50000000e0 / 3e0) * (0.3599999994e10 - 0.4799999982e18 * v[351] + 0.6e25 * Power(v[351], 3)
			+ 0.1799999982e26 * (v[351] * v[351]));
		v03 = 0.1e17 * (799999997e0 - 0.599999994e17 * v[351] - 0.3e17 * (v[351] * v[351]));
	};
	v[262] = v01;
	v[263] = v[248] * v[262];
	v[264] = v[249] * v[262];
	v[265] = v[250] * v[262];
	v[266] = sqrt((v[251] * v[251]) + (v[252] * v[252]) + (v[253] * v[253]));
	if (v[266] > 0.1e-7) { v04 = 1e0 / v[266]; v05 = (-(v04 / v[266])); v06 = (2e0 * v04) / (v[266] * v[266]); }
	else {
		v04 = (12500000e0 / 3e0) * (24e0 - (-0.1e-7 + v[266]) * (0.24e10 - 2e0 * (-1e0 + 100000000e0 * v[266]) * (0.2399999997e10
			- 0.1199999994e18 * v[266] - 0.3e17 * (v[266] * v[266]))));
		v05 = (-50000000e0 / 3e0) * (0.3599999994e10 - 0.4799999982e18 * v[266] + 0.6e25 * Power(v[266], 3)
			+ 0.1799999982e26 * (v[266] * v[266]));
		v06 = 0.1e17 * (799999997e0 - 0.599999994e17 * v[266] - 0.3e17 * (v[266] * v[266]));
	};
	v[271] = v04;
	v[272] = v[251] * v[271];
	v[273] = v[252] * v[271];
	v[274] = v[253] * v[271];
	if (v[263] * v[272] + v[264] * v[273] + v[265] * v[274] < 0e0) {
		v[277] = -v[263];
		v[278] = -v[264];
		v[279] = -v[265];
	}
	else {
		v[277] = v[263];
		v[278] = v[264];
		v[279] = v[265];
	};
	v[534] = (v[279] * v[279]);
	v[346] = 1e0 - v[534];
	v[533] = (v[278] * v[278]);
	v[341] = 1e0 - v[533];
	v[532] = (v[277] * v[277]);
	v[336] = 1e0 - v[532];
	if (sqrt(Power(-(v[273] * v[277]) + v[272] * v[278], 2) + Power(v[274] * v[277] - v[272] * v[279], 2) + Power(-
		(v[274] * v[278]) + v[273] * v[279], 2)) > 0.1e-7) {
		v[282] = -(v[274] * v[278]) + v[273] * v[279];
		v[283] = v[274] * v[277] - v[272] * v[279];
		v[284] = -(v[273] * v[277]) + v[272] * v[278];
		v[285] = sqrt((v[282] * v[282]) + (v[283] * v[283]) + (v[284] * v[284]));
		if (v[285] > 0.1e-7) { v07 = 1e0 / v[285]; v08 = (-(v07 / v[285])); v09 = (2e0 * v07) / (v[285] * v[285]); }
		else {
			v07 = (12500000e0 / 3e0) * (24e0 - (-0.1e-7 + v[285]) * (0.24e10 - 2e0 * (-1e0 + 100000000e0 * v[285]) *
				(0.2399999997e10 - 0.1199999994e18 * v[285] - 0.3e17 * (v[285] * v[285]))));
			v08 = (-50000000e0 / 3e0) * (0.3599999994e10 - 0.4799999982e18 * v[285] + 0.6e25 * Power(v[285], 3)
				+ 0.1799999982e26 * (v[285] * v[285]));
			v09 = 0.1e17 * (799999997e0 - 0.599999994e17 * v[285] - 0.3e17 * (v[285] * v[285]));
		};
		v[519] = 2e0 * v07 * tan(asin(v[285]) / 2e0);
		v[294] = v[282] * v[519];
		v[305] = (v[294] * v[294]);
		v[295] = v[283] * v[519];
		v[520] = v[295] / 2e0;
		v[303] = v[294] * v[520];
		v[298] = (v[295] * v[295]);
		v[296] = v[284] * v[519];
		v[310] = v[296] * v[520];
		v[308] = (v[294] * v[296]) / 2e0;
		v[299] = (v[296] * v[296]);
		v[521] = v[298] + v[299];
		v[297] = 4e0 / (4e0 + v[305] + v[521]);
		v[522] = -0.5e0 * v[297];
		v[300] = 1e0 + v[521] * v[522];
		v[301] = v[297] * (-v[296] + v[303]);
		v[302] = v[297] * (v[295] + v[308]);
		v[304] = v[297] * (v[296] + v[303]);
		v[306] = 1e0 + (v[299] + v[305]) * v[522];
		v[307] = v[297] * (-v[294] + v[310]);
		v[309] = v[297] * (-v[295] + v[308]);
		v[311] = v[297] * (v[294] + v[310]);
		v[312] = 1e0 - (-v[298] - v[305]) * v[522];
	}
	else {
		v[300] = 1e0;
		v[301] = 0e0;
		v[302] = 0e0;
		v[304] = 0e0;
		v[306] = 1e0;
		v[307] = 0e0;
		v[309] = 0e0;
		v[311] = 0e0;
		v[312] = 1e0;
	};
	if ((*previouscontact)) {
		v[317] = v[131] * v[165] + v[132] * v[169] + v[133] * v[170] - v[196] * v[230] - v[197] * v[234] - v[198] * v[235]
			+ v[1] * v[309] + v[2] * v[311] + v[3] * v[312] + v[523];
		v[526] = v[279] * v[317];
		v[316] = v[131] * v[159] + v[132] * v[163] + v[133] * v[164] - v[196] * v[224] - v[197] * v[228] - v[198] * v[229]
			+ v[1] * v[304] + v[2] * v[306] + v[3] * v[307] + v[524];
		v[528] = v[278] * v[316];
		v[315] = v[131] * v[153] + v[132] * v[157] + v[133] * v[158] - v[196] * v[218] - v[197] * v[222] - v[198] * v[223]
			+ v[1] * v[300] + v[2] * v[301] + v[3] * v[302] + v[525];
		v[527] = -(v[277] * v[315]);
		v[314] = v[315] * v[336] - v[277] * (v[526] + v[528]);
		v[318] = v[316] * v[341] + v[278] * (-v[526] + v[527]);
		v[319] = v[317] * v[346] + v[279] * (v[527] - v[528]);
	}
	else {
		v[314] = 0e0;
		v[318] = 0e0;
		v[319] = 0e0;
	};
	v[320] = v[134] * v[56] + v[150] * v[57] + v[151] * v[58];
	v[321] = -(v[150] * v[56]) + v[134] * v[57] + v[152] * v[58];
	v[322] = -(v[151] * v[56]) - v[152] * v[57] + v[134] * v[58];
	v[323] = v[102] * v[199] + v[103] * v[215] + v[104] * v[216];
	v[324] = v[103] * v[199] - v[102] * v[215] + v[104] * v[217];
	v[325] = v[104] * v[199] - v[102] * v[216] - v[103] * v[217];
	v[326] = duiA[0] - duiB[0] + v[130] * (-(v[157] * v[320]) + v[153] * v[321]) + v[129] * (v[158] * v[320] - v[153] * v[322])
		+ v[128] * (-(v[158] * v[321]) + v[157] * v[322]) + v[195] * (v[222] * v[323] - v[218] * v[324]) - v[194] * (v[223] * v[323]
			- v[218] * v[325]) + v[193] * (v[223] * v[324] - v[222] * v[325]);
	v[529] = v[277] * v[326];
	v[327] = duiA[1] - duiB[1] + v[130] * (-(v[163] * v[320]) + v[159] * v[321]) + v[129] * (v[164] * v[320] - v[159] * v[322])
		+ v[128] * (-(v[164] * v[321]) + v[163] * v[322]) + v[195] * (v[228] * v[323] - v[224] * v[324]) - v[194] * (v[229] * v[323]
			- v[224] * v[325]) + v[193] * (v[229] * v[324] - v[228] * v[325]);
	v[530] = v[278] * v[327];
	v[539] = v[279] * (v[529] + v[530]);
	v[337] = v[277] * v[530];
	v[328] = duiA[2] - duiB[2] + v[130] * (-(v[169] * v[320]) + v[165] * v[321]) + v[129] * (v[170] * v[320] - v[165] * v[322])
		+ v[128] * (-(v[170] * v[321]) + v[169] * v[322]) + v[195] * (v[234] * v[323] - v[230] * v[324]) - v[194] * (v[235] * v[323]
			- v[230] * v[325]) + v[193] * (v[235] * v[324] - v[234] * v[325]);
	v[531] = v[279] * v[328];
	v[538] = v[278] * (v[529] + v[531]);
	v[338] = v[277] * v[531];
	v[330] = v[337] + v[338] + v[326] * v[532];
	v[332] = v[327] * v[533] + v[538];
	v[334] = v[328] * v[534] + v[539];
	v[348] = -((v[536] * Power(v[349], v[364])) / (v[16] * Power(v[14], v[369])));
	b352 = v[351] < v[13];
	if (b352) {
		b354 = v[351] > v[14];
		if (b354) {
			v[363] = v[13] - v[351];
			v[535] = -(v[12] * Power(v[363], v[15]));
			v[356] = v[277] * v[535];
			v[358] = v[278] * v[535];
			v[359] = v[279] * v[535];
		}
		else {
			v[360] = -(v[12] * Power(v[349], v[15])) + v[348] * (Power(v[14], v[16]) - Power(v[351], v[16]));
			v[356] = v[277] * v[360];
			v[358] = v[278] * v[360];
			v[359] = v[279] * v[360];
		};
	}
	else {
		v[356] = 0e0;
		v[358] = 0e0;
		v[359] = 0e0;
	};
	if (b352) {
		v[537] = 2e0 * (*zetan);
		if (b354) {
			v[366] = v[537] * sqrt(v[536] * v[8] * Power(v[363], v[364]));
			v[365] = v[330] * v[366];
			v[367] = v[332] * v[366];
			v[368] = v[334] * v[366];
		}
		else {
			v[370] = v[537] * sqrt(-(v[16] * v[348] * v[8] * Power(v[351], v[369])));
			v[365] = v[330] * v[370];
			v[367] = v[332] * v[370];
			v[368] = v[334] * v[370];
		};
		if (v[351] < v[13] && v[277] * (v[356] + v[365]) + v[278] * (v[358] + v[367]) + v[279] * (v[359] + v[368]) < 0e0) {
			v[373] = v[365];
			v[374] = v[367];
			v[375] = v[368];
		}
		else {
			v[373] = -v[356];
			v[374] = -v[358];
			v[375] = -v[359];
		};
	}
	else {
		v[373] = 0e0;
		v[374] = 0e0;
		v[375] = 0e0;
	};
	v[376] = v[356] + v[373];
	v[377] = v[358] + v[374];
	v[378] = v[359] + v[375];
	v[379] = v[314] * v[4];
	v[380] = v[318] * v[4];
	v[381] = v[319] * v[4];
	v[385] = v[379] + (v[326] * v[336] - v[337] - v[338]) * v[9];
	v[386] = v[380] - (-(v[327] * v[341]) + v[538]) * v[9];
	v[387] = v[381] - (-(v[328] * v[346]) + v[539]) * v[9];
	if (b352) {
		if ((*stick)) {
			if (sqrt((v[385] * v[385]) + (v[386] * v[386]) + (v[387] * v[387])) <= v[5] * sqrt((v[376] * v[376]) +
				(v[377] * v[377]) + (v[378] * v[378]))) {
				v[392] = v[385];
				v[393] = v[386];
				v[394] = v[387];
				v[395] = 1e0;
			}
			else {
				v[396] = sqrt((v[385] * v[385]) + (v[386] * v[386]) + (v[387] * v[387]));
				if (v[396] > 0.1e-5) { v010 = 1e0 / v[396]; v011 = (-(v010 / v[396])); v012 = (2e0 * v010) / (v[396] * v[396]); }
				else {
					v010 = (24000000e0 - (-1e0 + 1000000e0 * v[396]) * (71999994e0 - 0.71999982e14 * v[396] + 0.6e19 * Power(v[396], 3)
						+ 0.23999982e20 * (v[396] * v[396]))) / 24e0;
					v011 = (-500000e0 / 3e0) * (35999994e0 - 0.47999982e14 * v[396] + 0.6e19 * Power(v[396], 3) + 0.17999982e20 *
						(v[396] * v[396]));
					v012 = 0.1e13 * (7999997e0 - 0.5999994e13 * v[396] - 0.3e13 * (v[396] * v[396]));
				};
				v[540] = v010 * v[6] * sqrt((v[376] * v[376]) + (v[377] * v[377]) + (v[378] * v[378]));
				v[392] = v[385] * v[540];
				v[393] = v[386] * v[540];
				v[394] = v[387] * v[540];
				v[395] = 0e0;
			};
			if (sqrt((v[379] * v[379]) + (v[380] * v[380]) + (v[381] * v[381])) > v[5] * sqrt((v[376] * v[376]) + (v[377] * v[377]
				) + (v[378] * v[378]))) {
				if (v[4] > 0.1e-5) { v013 = 1e0 / v[4]; v014 = (-(v013 / v[4])); v015 = (2e0 * v013) / (v[4] * v[4]); }
				else {
					v013 = (24000000e0 - (-1e0 + 1000000e0 * v[4]) * (71999994e0 - 0.71999982e14 * v[4] + 0.6e19 * Power(v[4], 3)
						+ 0.23999982e20 * (v[4] * v[4]))) / 24e0;
					v014 = (-500000e0 / 3e0) * (35999994e0 - 0.47999982e14 * v[4] + 0.6e19 * Power(v[4], 3) + 0.17999982e20 *
						(v[4] * v[4]));
					v015 = 0.1e13 * (7999997e0 - 0.5999994e13 * v[4] - 0.3e13 * (v[4] * v[4]));
				};
				v[414] = sqrt((v[379] * v[379]) + (v[380] * v[380]) + (v[381] * v[381]));
				if (v[414] > 0.1e-5) { v016 = 1e0 / v[414]; v017 = (-(v016 / v[414])); v018 = (2e0 * v016) / (v[414] * v[414]); }
				else {
					v016 = (24000000e0 - (-1e0 + 1000000e0 * v[414]) * (71999994e0 - 0.71999982e14 * v[414] + 0.6e19 * Power(v[414], 3)
						+ 0.23999982e20 * (v[414] * v[414]))) / 24e0;
					v017 = (-500000e0 / 3e0) * (35999994e0 - 0.47999982e14 * v[414] + 0.6e19 * Power(v[414], 3) + 0.17999982e20 *
						(v[414] * v[414]));
					v018 = 0.1e13 * (7999997e0 - 0.5999994e13 * v[414] - 0.3e13 * (v[414] * v[414]));
				};
				v[421] = -(v013 * v016 * v[6] * sqrt((v[376] * v[376]) + (v[377] * v[377]) + (v[378] * v[378])));
				v[420] = v[314] + v[379] * v[421];
				v[422] = v[318] + v[380] * v[421];
				v[423] = v[319] + v[381] * v[421];
			}
			else {
				v[420] = 0e0;
				v[422] = 0e0;
				v[423] = 0e0;
			};
		}
		else {
			if (sqrt((v[385] * v[385]) + (v[386] * v[386]) + (v[387] * v[387])) <= v[6] * sqrt((v[376] * v[376]) +
				(v[377] * v[377]) + (v[378] * v[378]))) {
				v[392] = v[385];
				v[393] = v[386];
				v[394] = v[387];
				v[395] = 1e0;
			}
			else {
				v[426] = sqrt((v[385] * v[385]) + (v[386] * v[386]) + (v[387] * v[387]));
				if (v[426] > 0.1e-5) { v019 = 1e0 / v[426]; v020 = (-(v019 / v[426])); v021 = (2e0 * v019) / (v[426] * v[426]); }
				else {
					v019 = (24000000e0 - (-1e0 + 1000000e0 * v[426]) * (71999994e0 - 0.71999982e14 * v[426] + 0.6e19 * Power(v[426], 3)
						+ 0.23999982e20 * (v[426] * v[426]))) / 24e0;
					v020 = (-500000e0 / 3e0) * (35999994e0 - 0.47999982e14 * v[426] + 0.6e19 * Power(v[426], 3) + 0.17999982e20 *
						(v[426] * v[426]));
					v021 = 0.1e13 * (7999997e0 - 0.5999994e13 * v[426] - 0.3e13 * (v[426] * v[426]));
				};
				v[541] = v019 * v[6] * sqrt((v[376] * v[376]) + (v[377] * v[377]) + (v[378] * v[378]));
				v[392] = v[385] * v[541];
				v[393] = v[386] * v[541];
				v[394] = v[387] * v[541];
				v[395] = 0e0;
			};
			if (sqrt((v[379] * v[379]) + (v[380] * v[380]) + (v[381] * v[381])) > v[6] * sqrt((v[376] * v[376]) + (v[377] * v[377]
				) + (v[378] * v[378]))) {
				if (v[4] > 0.1e-5) { v022 = 1e0 / v[4]; v023 = (-(v022 / v[4])); v024 = (2e0 * v022) / (v[4] * v[4]); }
				else {
					v022 = (24000000e0 - (-1e0 + 1000000e0 * v[4]) * (71999994e0 - 0.71999982e14 * v[4] + 0.6e19 * Power(v[4], 3)
						+ 0.23999982e20 * (v[4] * v[4]))) / 24e0;
					v023 = (-500000e0 / 3e0) * (35999994e0 - 0.47999982e14 * v[4] + 0.6e19 * Power(v[4], 3) + 0.17999982e20 *
						(v[4] * v[4]));
					v024 = 0.1e13 * (7999997e0 - 0.5999994e13 * v[4] - 0.3e13 * (v[4] * v[4]));
				};
				v[444] = sqrt((v[379] * v[379]) + (v[380] * v[380]) + (v[381] * v[381]));
				if (v[444] > 0.1e-5) { v025 = 1e0 / v[444]; v026 = (-(v025 / v[444])); v027 = (2e0 * v025) / (v[444] * v[444]); }
				else {
					v025 = (24000000e0 - (-1e0 + 1000000e0 * v[444]) * (71999994e0 - 0.71999982e14 * v[444] + 0.6e19 * Power(v[444], 3)
						+ 0.23999982e20 * (v[444] * v[444]))) / 24e0;
					v026 = (-500000e0 / 3e0) * (35999994e0 - 0.47999982e14 * v[444] + 0.6e19 * Power(v[444], 3) + 0.17999982e20 *
						(v[444] * v[444]));
					v027 = 0.1e13 * (7999997e0 - 0.5999994e13 * v[444] - 0.3e13 * (v[444] * v[444]));
				};
				v[450] = -(v022 * v025 * v[6] * sqrt((v[376] * v[376]) + (v[377] * v[377]) + (v[378] * v[378])));
				v[420] = v[314] + v[379] * v[450];
				v[422] = v[318] + v[380] * v[450];
				v[423] = v[319] + v[381] * v[450];
			}
			else {
				v[420] = 0e0;
				v[422] = 0e0;
				v[423] = 0e0;
			};
		};
	}
	else {
		v[392] = 0e0;
		v[393] = 0e0;
		v[394] = 0e0;
	};
	v[544] = v[378] + v[394];
	v[543] = v[377] + v[393];
	v[542] = v[376] + v[392];
	fn[0] = v[376];
	fn[1] = v[377];
	fn[2] = v[378];
	ft[0] = v[392];
	ft[1] = v[393];
	ft[2] = v[394];
	(*stickupdated) = v[395];
	gtpupdated[0] = v[314] - v[420];
	gtpupdated[1] = v[318] - v[422];
	gtpupdated[2] = v[319] - v[423];
	/*Rc[0] = v[153] * v[542] + v[159] * v[543] + v[165] * v[544];
	Rc[1] = v[157] * v[542] + v[163] * v[543] + v[169] * v[544];
	Rc[2] = v[158] * v[542] + v[164] * v[543] + v[170] * v[544];
	Rc[3] = (-(v[130] * v[157]) + v[129] * v[158]) * v[542] + (-(v[130] * v[163]) + v[129] * v[164]) * v[543] + (-
		(v[130] * v[169]) + v[129] * v[170]) * v[544];
	Rc[4] = (v[130] * v[153] - v[128] * v[158]) * v[542] + (v[130] * v[159] - v[128] * v[164]) * v[543] + (v[130] * v[165]
		- v[128] * v[170]) * v[544];
	Rc[5] = (-(v[129] * v[153]) + v[128] * v[157]) * v[542] + (-(v[129] * v[159]) + v[128] * v[163]) * v[543] + (-
		(v[129] * v[165]) + v[128] * v[169]) * v[544];
	Rc[6] = -(v[218] * v[542]) - v[224] * v[543] - v[230] * v[544];
	Rc[7] = -(v[222] * v[542]) - v[228] * v[543] - v[234] * v[544];
	Rc[8] = -(v[223] * v[542]) - v[229] * v[543] - v[235] * v[544];
	Rc[9] = -((-(v[195] * v[222]) + v[194] * v[223]) * v[542]) - (-(v[195] * v[228]) + v[194] * v[229]) * v[543] - (-
		(v[195] * v[234]) + v[194] * v[235]) * v[544];
	Rc[10] = -((v[195] * v[218] - v[193] * v[223]) * v[542]) - (v[195] * v[224] - v[193] * v[229]) * v[543] - (v[195] * v[230]
		- v[193] * v[235]) * v[544];
	Rc[11] = -((-(v[194] * v[218]) + v[193] * v[222]) * v[542]) - (-(v[194] * v[224]) + v[193] * v[228]) * v[543] - (-
		(v[194] * v[230]) + v[193] * v[234]) * v[544];*/

	delete[]v;
}