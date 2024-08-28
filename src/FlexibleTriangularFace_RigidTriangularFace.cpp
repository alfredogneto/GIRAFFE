#include "FlexibleTriangularFace_RigidTriangularFace.h"

#include "VEMPolyhedron.h"
#include "STLSurface.h"
#include "TriangularFace.h"
#include "Interface_1.h"
#include "Matrix.h"
#include "SSContactData.h"
#include "ExecutionData.h"
#include "Material.h"
#include "SuperNode.h"
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


FlexibleTriangularFace_RigidTriangularFace::FlexibleTriangularFace_RigidTriangularFace()
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
	super_node_A = 0;		//ID of super node A
	node_B = 0;				//ID of node B
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
}


FlexibleTriangularFace_RigidTriangularFace::~FlexibleTriangularFace_RigidTriangularFace()
{
	SetUnnactive();
	delete GammaA;
	delete GammaB;

	delete meq;
	
}
void FlexibleTriangularFace_RigidTriangularFace::AllocSpecific()
{
	if (alloc_specific_control == false)
	{
		cAp = DBG_NEW double[2];
		cBp = DBG_NEW double[2];
		cAi = DBG_NEW double[2];
		cBi = DBG_NEW double[2];
		xi1A = DBG_NEW double[3];
		xi2A = DBG_NEW double[3];
		xi3A = DBG_NEW double[3];
		u1A = DBG_NEW double[3];
		u2A = DBG_NEW double[3];
		u3A = DBG_NEW double[3];
		dui1A = DBG_NEW double[3];
		dui2A = DBG_NEW double[3];
		dui3A = DBG_NEW double[3];
		ddui1A = DBG_NEW double[3];
		ddui2A = DBG_NEW double[3];
		ddui3A = DBG_NEW double[3];

		
		xBi = DBG_NEW double[3];
		uB = DBG_NEW double[3];
		alphaB = DBG_NEW double[3];
		duiB = DBG_NEW double[3];
		dalphaiB = DBG_NEW double[3];
		dduiB = DBG_NEW double[3];
		ddalphaiB = DBG_NEW double[3];
		QBi = DBG_NEW double*[3];
		for (int i = 0; i < 3; i++)
			QBi[i] = DBG_NEW double[3];
		

		fn = DBG_NEW double[3];
		ft = DBG_NEW double[3];

		Rc = DBG_NEW double[15];
		Kc = DBG_NEW double*[15];

		for (int i = 0; i < 15; i++)
			Kc[i] = DBG_NEW double[15];
		for (int ni = 0; ni < 15; ni++)
			for (int nj = 0; nj < 15; nj++)
				Kc[ni][nj] = 0.0;

		alloc_specific_control = true;
		vnrel = new double;
	}
}

void FlexibleTriangularFace_RigidTriangularFace::FreeSpecific()
{
	if (alloc_specific_control == true)
	{
		delete[] cAp;
		delete[] cBp;
		delete[] cAi;
		delete[] cBi;

		delete[] xi1A;
		delete[] xi2A;
		delete[] xi3A;
		delete[] u1A;
		delete[] u2A;
		delete[] u3A;
		delete[] dui1A;
		delete[] dui2A;
		delete[] dui3A;
		delete[] ddui1A;
		delete[] ddui2A;
		delete[] ddui3A;

		delete[] xBi;
		delete[] uB;
		delete[] alphaB;
		delete[] duiB;
		delete[] dalphaiB;
		delete[] dduiB;
		delete[] ddalphaiB;
		for (int i = 0; i < 3; i++)
			delete[] QBi[i];
		delete[]QBi;

		delete[] fn;
		delete[] ft;

		delete[] Rc;
		for (int i = 0; i < 15; i++)
			delete[]Kc[i];
		delete[]Kc;

		alloc_specific_control = false;
		delete vnrel;
	}
}

bool FlexibleTriangularFace_RigidTriangularFace::HaveErrors()
{
	if (cd->copy_return_value[0] == 2)
		SolvePreviousContact();
	if (dot(*cd->n[0], *cd->copy_n[0]) < 0.0)
	{
		if (cd->return_value[0] == 0)
		{
			if (deg_pointA != 0 || deg_pointB != 0)
			{
				//Contact involving a vertex 
				db.myprintf("\nPenetration detected in contact evaluation.");
				if (db.execution_data->print_contact_report)
				{
					db.myprintf("\nFlexibleTriangular-RigidTriangular I1: %d I2: %d \nfaceA: %d curveA: %d pointA: %d\nfaceB: %d curveB: %d pointB: %d\n", index1, index2, faceA->ID, deg_curveA, deg_pointA, faceB->ID, deg_curveB, deg_pointB);
					cd->Plot();
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
						db.myprintf("\nFlexibleTriangular-RigidTriangular I1: %d I2: %d \nfaceA: %d curveA: %d pointA: %d\nfaceB: %d curveB: %d pointB: %d\n", index1, index2, faceA->ID, deg_curveA, deg_pointA, faceB->ID, deg_curveB, deg_pointB);
						cd->Plot();
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

void FlexibleTriangularFace_RigidTriangularFace::PreCalc()
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


void FlexibleTriangularFace_RigidTriangularFace::EvaluateNormalGap()
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


void FlexibleTriangularFace_RigidTriangularFace::SolveLCP()
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
		Matrix xp;
		if (previous_evaluation == false)
		{
			x1Bl = *ptrx0Bp + (*ptrQBp) * x1Bl;
			x2Bl = *ptrx0Bp + (*ptrQBp) * x2Bl;
			x3Bl = *ptrx0Bp + (*ptrQBp) * x3Bl;
			xp = *vertices_p_A[vertexIDsA[0] - 1];
		}
		else
		{
			x1Bl = *ptrx0Bi + (*ptrQBi) * x1Bl;
			x2Bl = *ptrx0Bi + (*ptrQBi) * x2Bl;
			x3Bl = *ptrx0Bi + (*ptrQBi) * x3Bl;
			xp = *vertices_i_A[vertexIDsA[0] - 1];
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
		Matrix x1Al;
		Matrix x2Al;
		Matrix x3Al;
		//Surface B
		Matrix xp = *surfB->vertices[vertexIDsB[0] - 1].coord_double;

		if (previous_evaluation == false)
		{
			x1Al = *vertices_p_A[vertexIDsA[0] - 1];
			x2Al = *vertices_p_A[vertexIDsA[1] - 1];
			x3Al = *vertices_p_A[vertexIDsA[2] - 1];
			xp = *ptrx0Bp + (*ptrQBp) * xp;
		}
		else
		{
			x1Al = *vertices_i_A[vertexIDsA[0] - 1];
			x2Al = *vertices_i_A[vertexIDsA[1] - 1];
			x3Al = *vertices_i_A[vertexIDsA[2] - 1];
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
		Matrix x1Al;
		Matrix x2Al;
		//Surface B
		Matrix x1Bl = *surfB->vertices[vertexIDsB[0] - 1].coord_double;
		Matrix x2Bl = *surfB->vertices[vertexIDsB[1] - 1].coord_double;

		if (previous_evaluation == false)
		{
			x1Al = *vertices_p_A[vertexIDsA[0] - 1];
			x2Al = *vertices_p_A[vertexIDsA[1] - 1];
			x1Bl = *ptrx0Bp + (*ptrQBp) * x1Bl;
			x2Bl = *ptrx0Bp + (*ptrQBp) * x2Bl;
		}
		else
		{
			x1Al = *vertices_i_A[vertexIDsA[0] - 1];
			x2Al = *vertices_i_A[vertexIDsA[1] - 1];
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
		Matrix x1Al;
		//Surface B - edge
		Matrix x1Bl = *surfB->vertices[vertexIDsB[0] - 1].coord_double;
		Matrix x2Bl = *surfB->vertices[vertexIDsB[1] - 1].coord_double;

		if (previous_evaluation == false)
		{
			x1Al = *vertices_p_A[vertexIDsA[0] - 1];
			x1Bl = *ptrx0Bp + (*ptrQBp) * x1Bl;
			x2Bl = *ptrx0Bp + (*ptrQBp) * x2Bl;
		}
		else
		{
			x1Al = *vertices_i_A[vertexIDsA[0] - 1];
			x1Bl = *ptrx0Bi + (*ptrQBi) * x1Bl;
			x2Bl = *ptrx0Bi + (*ptrQBi) * x2Bl;
		}

		Matrix b = 0.5*(x1Bl + x2Bl);
		Matrix t = 0.5*(x2Bl - x1Bl);

		double zeta = (1.0 / (dot(t, t))) * dot(x1Al - b, t);
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
		Matrix x1Al;
		Matrix x2Al;

		if (previous_evaluation == false)
		{
			x1Bl = *ptrx0Bp + (*ptrQBp) * x1Bl;
			x1Al = *vertices_p_A[vertexIDsA[0] - 1];
			x2Al = *vertices_p_A[vertexIDsA[1] - 1];
		}
		else
		{
			x1Bl = *ptrx0Bi + (*ptrQBi) * x1Bl;
			x1Al = *vertices_i_A[vertexIDsA[0] - 1];
			x2Al = *vertices_i_A[vertexIDsA[1] - 1];
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
void FlexibleTriangularFace_RigidTriangularFace::SolvePreviousContact()
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

void FlexibleTriangularFace_RigidTriangularFace::SurfacePoints()
{
	double z;
	double th;
	double N1;
	double N2;
	double N3;
	//Gamma A
	if (previous_evaluation == false)
	{
		z = cd->convective[0][0];
		th = cd->convective[0][1];
		N1 = -(th / 2.0) - z / 2.0;
		N2 = 0.5 + z / 2.0;
		N3 = 0.5 + th / 2.0;
		*GammaA = N1 * (*vertices_p_A[vertexIDsA[0] - 1])
			+ N2 * (*vertices_p_A[vertexIDsA[1] - 1])
			+ N3 * (*vertices_p_A[vertexIDsA[2] - 1]);
	}
	else
	{
		z = cd->copy_convective[0][0];
		th = cd->copy_convective[0][1];
		N1 = -(th / 2.0) - z / 2.0;
		N2 = 0.5 + z / 2.0;
		N3 = 0.5 + th / 2.0;
		*GammaA = N1 * (*vertices_i_A[vertexIDsA[0] - 1])
			+ N2 * (*vertices_i_A[vertexIDsA[1] - 1])
			+ N3 * (*vertices_i_A[vertexIDsA[2] - 1]);
	}

	//Gamma B
	Matrix localx;
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


void FlexibleTriangularFace_RigidTriangularFace::SetVariables()
{
	for (int i = 0; i < 3; i++)
	{
		xi1A[i] = db.super_nodes[super_node_A - 1]->copy_coordinates[(vertexIDsA[0] - 1) * 3 + i];
		xi2A[i] = db.super_nodes[super_node_A - 1]->copy_coordinates[(vertexIDsA[1] - 1) * 3 + i];
		xi3A[i] = db.super_nodes[super_node_A - 1]->copy_coordinates[(vertexIDsA[2] - 1) * 3 + i];

		u1A[i] = db.super_nodes[super_node_A - 1]->displacements[(vertexIDsA[0] - 1) * 3 + i];
		u2A[i] = db.super_nodes[super_node_A - 1]->displacements[(vertexIDsA[1] - 1) * 3 + i];
		u3A[i] = db.super_nodes[super_node_A - 1]->displacements[(vertexIDsA[2] - 1) * 3 + i];

		dui1A[i] = db.super_nodes[super_node_A - 1]->copy_vel[(vertexIDsA[0] - 1) * 3 + i];
		dui2A[i] = db.super_nodes[super_node_A - 1]->copy_vel[(vertexIDsA[1] - 1) * 3 + i];
		dui3A[i] = db.super_nodes[super_node_A - 1]->copy_vel[(vertexIDsA[2] - 1) * 3 + i];

		ddui1A[i] = db.super_nodes[super_node_A - 1]->copy_accel[(vertexIDsA[0] - 1) * 3 + i];
		ddui2A[i] = db.super_nodes[super_node_A - 1]->copy_accel[(vertexIDsA[1] - 1) * 3 + i];
		ddui3A[i] = db.super_nodes[super_node_A - 1]->copy_accel[(vertexIDsA[2] - 1) * 3 + i];

		
		xBi[i] = db.nodes[node_B - 1]->copy_coordinates[i];
		uB[i] = db.nodes[node_B - 1]->displacements[i];
		alphaB[i] = db.nodes[node_B - 1]->displacements[i + 3];
		duiB[i] = db.nodes[node_B - 1]->copy_vel[i];
		dduiB[i] = db.nodes[node_B - 1]->copy_accel[i];
		dalphaiB[i] = db.nodes[node_B - 1]->copy_vel[i + 3];
		ddalphaiB[i] = db.nodes[node_B - 1]->copy_accel[i + 3];
		for (int j = 0; j < 3; j++)
			QBi[i][j] = (*ptrQBi)(i, j);
	}
}


void FlexibleTriangularFace_RigidTriangularFace::MountLocalContributions()
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
	for (int i = 0; i < 15; i++)
	{
		Rc[i] = 0.0;
		for (int j = 0; j < 15; j++)
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
//	int i705, i753, i933, i1225, i1998, i2139, i2660, i2661, i2662, i2710, i2711, i2712, i2718
//		, i2719, i2720, i2729, i2730, i2731, i2732, i2733, i2734, i2735, i2736, i2737, i2864, i2865
//		, i2866, b537, b538, b570, b598, b599, b600, b601, b607, b608, b617, b618, b634, b635, b636
//		, b637, b653, b670, b671, b683, b710, b711, b803, b804, b938, b939, b940, b941, b969, b973
//		, b974, b986, b987, b1049, b1080, b1448, b1493, b1571, b1572, b1600, b1601, b1627, b1628
//		, b1629, b1646, b1714, b1718, b1722, b1731, b1732, b1766, b1788, b2007, b2008, b2012, b2017
//		, b2018, b2335, b2336, b2351, b2355, b2359, b2364, b2365;
//	v[2671] = -((*a4)*(*ct));
//	v[2658] = (*epsn1)*(*n1);
//	v[595] = (*gnb) - (*gnbb);
//	v[610] = -1e0 + (*n1);
//	v[1602] = -1e0 + v[610];
//	v[615] = -1e0 + (*n2);
//	v[1607] = -1e0 + v[615];
//	v[144] = u1A[0] + xi1A[0];
//	v[281] = -v[144] / 2e0;
//	v[148] = u1A[1] + xi1A[1];
//	v[283] = -v[148] / 2e0;
//	v[152] = u1A[2] + xi1A[2];
//	v[285] = -v[152] / 2e0;
//	v[145] = u2A[0] + xi2A[0];
//	v[277] = v[145] / 2e0 + v[281];
//	v[149] = u2A[1] + xi2A[1];
//	v[278] = v[149] / 2e0 + v[283];
//	v[153] = u2A[2] + xi2A[2];
//	v[279] = v[153] / 2e0 + v[285];
//	v[146] = u3A[0] + xi3A[0];
//	v[282] = v[146] / 2e0 + v[281];
//	v[150] = u3A[1] + xi3A[1];
//	v[284] = v[150] / 2e0 + v[283];
//	v[154] = u3A[2] + xi3A[2];
//	v[286] = v[154] / 2e0 + v[285];
//	v[129] = -cAp[0] / 2e0;
//	v[131] = -cAp[1] / 2e0;
//	v[134] = -cAi[0] / 2e0;
//	v[136] = -cAi[1] / 2e0;
//	v[295] = -x1B[0] / 2e0;
//	v[297] = -x1B[1] / 2e0;
//	v[299] = -x1B[2] / 2e0;
//	v[288] = v[295] + x2B[0] / 2e0;
//	v[289] = v[297] + x2B[1] / 2e0;
//	v[290] = v[299] + x2B[2] / 2e0;
//	v[296] = v[295] + x3B[0] / 2e0;
//	v[298] = v[297] + x3B[1] / 2e0;
//	v[300] = v[299] + x3B[2] / 2e0;
//	v[307] = alphaB[0] / 2e0;
//	v[305] = 2e0*alphaB[0];
//	v[180] = Power(alphaB[0], 2);
//	v[308] = 2e0*alphaB[1];
//	v[4365] = 0e0;
//	v[4366] = 0e0;
//	v[4367] = 0e0;
//	v[4368] = 0e0;
//	v[4369] = 0e0;
//	v[4370] = 0e0;
//	v[4371] = 0e0;
//	v[4372] = 0e0;
//	v[4373] = 0e0;
//	v[4374] = 0e0;
//	v[4375] = 0e0;
//	v[4376] = 0e0;
//	v[4377] = v[305];
//	v[4378] = v[308];
//	v[4379] = 0e0;
//	v[4703] = 0e0;
//	v[4704] = 0e0;
//	v[4705] = 0e0;
//	v[4706] = 0e0;
//	v[4707] = 0e0;
//	v[4708] = 0e0;
//	v[4709] = 0e0;
//	v[4710] = 0e0;
//	v[4711] = 0e0;
//	v[4712] = 0e0;
//	v[4713] = 0e0;
//	v[4714] = 0e0;
//	v[4715] = -v[305];
//	v[4716] = -v[308];
//	v[4717] = 0e0;
//	v[306] = alphaB[1] / 2e0;
//	v[4200] = 0e0;
//	v[4201] = 0e0;
//	v[4202] = 0e0;
//	v[4203] = 0e0;
//	v[4204] = 0e0;
//	v[4205] = 0e0;
//	v[4206] = 0e0;
//	v[4207] = 0e0;
//	v[4208] = 0e0;
//	v[4209] = 0e0;
//	v[4210] = 0e0;
//	v[4211] = 0e0;
//	v[4212] = v[306];
//	v[4213] = v[307];
//	v[4214] = 0e0;
//	v[178] = (alphaB[0] * alphaB[1]) / 2e0;
//	v[173] = Power(alphaB[1], 2);
//	v[351] = -v[173] - v[180];
//	v[329] = alphaB[2] + v[178];
//	v[319] = -alphaB[2] + v[178];
//	v[310] = 2e0*alphaB[2];
//	v[4305] = 0e0;
//	v[4306] = 0e0;
//	v[4307] = 0e0;
//	v[4308] = 0e0;
//	v[4309] = 0e0;
//	v[4310] = 0e0;
//	v[4311] = 0e0;
//	v[4312] = 0e0;
//	v[4313] = 0e0;
//	v[4314] = 0e0;
//	v[4315] = 0e0;
//	v[4316] = 0e0;
//	v[4317] = v[305];
//	v[4318] = 0e0;
//	v[4319] = v[310];
//	v[4245] = 0e0;
//	v[4246] = 0e0;
//	v[4247] = 0e0;
//	v[4248] = 0e0;
//	v[4249] = 0e0;
//	v[4250] = 0e0;
//	v[4251] = 0e0;
//	v[4252] = 0e0;
//	v[4253] = 0e0;
//	v[4254] = 0e0;
//	v[4255] = 0e0;
//	v[4256] = 0e0;
//	v[4257] = 0e0;
//	v[4258] = v[308];
//	v[4259] = v[310];
//	v[4763] = 0e0;
//	v[4764] = 0e0;
//	v[4765] = 0e0;
//	v[4766] = 0e0;
//	v[4767] = 0e0;
//	v[4768] = 0e0;
//	v[4769] = 0e0;
//	v[4770] = 0e0;
//	v[4771] = 0e0;
//	v[4772] = 0e0;
//	v[4773] = 0e0;
//	v[4774] = 0e0;
//	v[4775] = 0e0;
//	v[4776] = -v[308];
//	v[4777] = -v[310];
//	v[4185] = 0e0;
//	v[4186] = 0e0;
//	v[4187] = 0e0;
//	v[4188] = 0e0;
//	v[4189] = 0e0;
//	v[4190] = 0e0;
//	v[4191] = 0e0;
//	v[4192] = 0e0;
//	v[4193] = 0e0;
//	v[4194] = 0e0;
//	v[4195] = 0e0;
//	v[4196] = 0e0;
//	v[4197] = v[305];
//	v[4198] = v[308];
//	v[4199] = v[310];
//	v[309] = alphaB[2] / 2e0;
//	v[4215] = 0e0;
//	v[4216] = 0e0;
//	v[4217] = 0e0;
//	v[4218] = 0e0;
//	v[4219] = 0e0;
//	v[4220] = 0e0;
//	v[4221] = 0e0;
//	v[4222] = 0e0;
//	v[4223] = 0e0;
//	v[4224] = 0e0;
//	v[4225] = 0e0;
//	v[4226] = 0e0;
//	v[4227] = v[309];
//	v[4228] = 0e0;
//	v[4229] = v[307];
//	v[4230] = 0e0;
//	v[4231] = 0e0;
//	v[4232] = 0e0;
//	v[4233] = 0e0;
//	v[4234] = 0e0;
//	v[4235] = 0e0;
//	v[4236] = 0e0;
//	v[4237] = 0e0;
//	v[4238] = 0e0;
//	v[4239] = 0e0;
//	v[4240] = 0e0;
//	v[4241] = 0e0;
//	v[4242] = 0e0;
//	v[4243] = v[309];
//	v[4244] = v[306];
//	v[185] = (alphaB[1] * alphaB[2]) / 2e0;
//	v[346] = alphaB[0] + v[185];
//	v[338] = -alphaB[0] + v[185];
//	v[183] = (alphaB[0] * alphaB[2]) / 2e0;
//	v[342] = -alphaB[1] + v[183];
//	v[325] = alphaB[1] + v[183];
//	v[174] = Power(alphaB[2], 2);
//	v[859] = 4e0 + v[173] + v[174] + v[180];
//	v[2576] = 1e0 / Power(v[859], 4);
//	v[1215] = 1e0 / Power(v[859], 3);
//	v[1272] = -8e0*v[1215] * v[305];
//	v[1270] = 8e0*v[1215] * v[308];
//	v[1267] = -8e0*v[1215] * v[310];
//	v[783] = 1e0 / Power(v[859], 2);
//	v[2935] = 2e0*v[783];
//	v[2934] = 4e0*v[783];
//	v[2827] = 2e0*v[783];
//	v[2826] = 4e0*v[783];
//	v[333] = -v[174] - v[180];
//	v[314] = -v[173] - v[174];
//	v[313] = 4e0*v[310] * v[783];
//	v[355] = -(v[313] * v[351]) / 2e0;
//	v[312] = -4e0*v[308] * v[783];
//	v[335] = (v[312] * v[333]) / 2e0;
//	v[311] = 4e0*v[305] * v[783];
//	v[315] = -(v[311] * v[314]) / 2e0;
//	v[230] = (*a4)*alphaB[0] + (*a5)*dalphaiB[0] + (*a6)*ddalphaiB[0];
//	v[236] = (*a4)*alphaB[1] + (*a5)*dalphaiB[1] + (*a6)*ddalphaiB[1];
//	v[238] = (*a4)*alphaB[2] + (*a5)*dalphaiB[2] + (*a6)*ddalphaiB[2];
//	v[157] = -cBp[0] / 2e0;
//	v[159] = -cBp[1] / 2e0;
//	v[162] = -cBi[0] / 2e0;
//	v[164] = -cBi[1] / 2e0;
//	v[128] = v[129] + v[131];
//	v[130] = 0.5e0 - v[129];
//	v[430] = v[130] * v[286];
//	v[429] = v[130] * v[284];
//	v[428] = v[130] * v[282];
//	v[132] = 0.5e0 - v[131];
//	v[4988] = -v[128];
//	v[4989] = 0e0;
//	v[4990] = 0e0;
//	v[4991] = -v[130];
//	v[4992] = 0e0;
//	v[4993] = 0e0;
//	v[4994] = -v[132];
//	v[4995] = 0e0;
//	v[4996] = 0e0;
//	v[4997] = 1e0;
//	v[4998] = 0e0;
//	v[4999] = 0e0;
//	v[5000] = 0e0;
//	v[5001] = 0e0;
//	v[5002] = 0e0;
//	v[5003] = 0e0;
//	v[5004] = -v[128];
//	v[5005] = 0e0;
//	v[5006] = 0e0;
//	v[5007] = -v[130];
//	v[5008] = 0e0;
//	v[5009] = 0e0;
//	v[5010] = -v[132];
//	v[5011] = 0e0;
//	v[5012] = 0e0;
//	v[5013] = 1e0;
//	v[5014] = 0e0;
//	v[5015] = 0e0;
//	v[5016] = 0e0;
//	v[5017] = 0e0;
//	v[4410] = 0e0;
//	v[4411] = 0e0;
//	v[4412] = v[128];
//	v[4413] = 0e0;
//	v[4414] = 0e0;
//	v[4415] = v[130];
//	v[4416] = 0e0;
//	v[4417] = 0e0;
//	v[4418] = v[132];
//	v[4419] = 0e0;
//	v[4420] = 0e0;
//	v[4421] = -1e0;
//	v[4422] = 0e0;
//	v[4423] = 0e0;
//	v[4424] = 0e0;
//	v[4395] = 0e0;
//	v[4396] = v[128];
//	v[4397] = 0e0;
//	v[4398] = 0e0;
//	v[4399] = v[130];
//	v[4400] = 0e0;
//	v[4401] = 0e0;
//	v[4402] = v[132];
//	v[4403] = 0e0;
//	v[4404] = 0e0;
//	v[4405] = -1e0;
//	v[4406] = 0e0;
//	v[4407] = 0e0;
//	v[4408] = 0e0;
//	v[4409] = 0e0;
//	v[4380] = v[128];
//	v[4381] = 0e0;
//	v[4382] = 0e0;
//	v[4383] = v[130];
//	v[4384] = 0e0;
//	v[4385] = 0e0;
//	v[4386] = v[132];
//	v[4387] = 0e0;
//	v[4388] = 0e0;
//	v[4389] = -1e0;
//	v[4390] = 0e0;
//	v[4391] = 0e0;
//	v[4392] = 0e0;
//	v[4393] = 0e0;
//	v[4394] = 0e0;
//	v[421] = v[132] * v[279];
//	v[420] = v[132] * v[278];
//	v[419] = v[132] * v[277];
//	v[133] = v[134] + v[136];
//	v[135] = 0.5e0 - v[134];
//	v[137] = 0.5e0 - v[136];
//	v[156] = v[157] + v[159];
//	v[158] = 0.5e0 - v[157];
//	v[160] = 0.5e0 - v[159];
//	v[161] = v[162] + v[164];
//	v[163] = 0.5e0 - v[162];
//	v[165] = 0.5e0 - v[164];
//	v[166] = v[156] * x1B[0] + v[158] * x2B[0] + v[160] * x3B[0];
//	v[167] = v[156] * x1B[1] + v[158] * x2B[1] + v[160] * x3B[1];
//	v[168] = v[156] * x1B[2] + v[158] * x2B[2] + v[160] * x3B[2];
//	v[722] = -(QBi[0][0] * v[166]) - QBi[0][1] * v[167] - QBi[0][2] * v[168];
//	v[720] = -(QBi[1][0] * v[166]) - QBi[1][1] * v[167] - QBi[1][2] * v[168];
//	v[718] = -(QBi[2][0] * v[166]) - QBi[2][1] * v[167] - QBi[2][2] * v[168];
//	v[169] = v[161] * x1B[0] + v[163] * x2B[0] + v[165] * x3B[0];
//	v[170] = v[161] * x1B[1] + v[163] * x2B[1] + v[165] * x3B[1];
//	v[171] = v[161] * x1B[2] + v[163] * x2B[2] + v[165] * x3B[2];
//	v[172] = 4e0 / v[859];
//	v[2824] = -v[172] / 2e0;
//	v[353] = -(v[172] * v[308]) / 2e0;
//	v[354] = (v[312] * v[351]) / 2e0 + v[353];
//	v[350] = -(v[172] * v[305]) / 2e0;
//	v[352] = v[350] - (v[311] * v[351]) / 2e0;
//	v[347] = v[172] - v[311] * v[346];
//	v[344] = -v[172] + v[312] * v[342];
//	v[339] = -v[172] - v[311] * v[338];
//	v[336] = -(v[172] * v[310]) / 2e0;
//	v[337] = -(v[313] * v[333]) / 2e0 + v[336];
//	v[334] = -(v[311] * v[333]) / 2e0 + v[350];
//	v[332] = v[172] - v[313] * v[329];
//	v[327] = v[172] + v[312] * v[325];
//	v[324] = -(v[172] * v[309]);
//	v[348] = -v[324] + v[312] * v[346];
//	v[393] = QBi[0][2] * v[344] + QBi[1][2] * v[348] + QBi[2][2] * v[354];
//	v[390] = QBi[0][1] * v[344] + QBi[1][1] * v[348] + QBi[2][1] * v[354];
//	v[387] = QBi[0][0] * v[344] + QBi[1][0] * v[348] + QBi[2][0] * v[354];
//	v[408] = v[166] * v[387] + v[167] * v[390] + v[168] * v[393];
//	v[343] = -v[324] - v[311] * v[342];
//	v[392] = QBi[0][2] * v[343] + QBi[1][2] * v[347] + QBi[2][2] * v[352];
//	v[389] = QBi[0][1] * v[343] + QBi[1][1] * v[347] + QBi[2][1] * v[352];
//	v[386] = QBi[0][0] * v[343] + QBi[1][0] * v[347] + QBi[2][0] * v[352];
//	v[407] = v[166] * v[386] + v[167] * v[389] + v[168] * v[392];
//	v[340] = -v[324] + v[312] * v[338];
//	v[326] = -v[324] - v[311] * v[325];
//	v[323] = -v[172] - v[313] * v[319];
//	v[321] = -(v[172] * v[307]);
//	v[345] = -v[321] - v[313] * v[342];
//	v[331] = -v[321] + v[312] * v[329];
//	v[378] = QBi[0][2] * v[331] + QBi[1][2] * v[335] + QBi[2][2] * v[340];
//	v[375] = QBi[0][1] * v[331] + QBi[1][1] * v[335] + QBi[2][1] * v[340];
//	v[372] = QBi[0][0] * v[331] + QBi[1][0] * v[335] + QBi[2][0] * v[340];
//	v[405] = v[166] * v[372] + v[167] * v[375] + v[168] * v[378];
//	v[328] = -v[321] - v[313] * v[325];
//	v[322] = v[312] * v[319] - v[321];
//	v[318] = v[172] * v[306];
//	v[349] = v[318] - v[313] * v[346];
//	v[394] = QBi[0][2] * v[345] + QBi[1][2] * v[349] + QBi[2][2] * v[355];
//	v[391] = QBi[0][1] * v[345] + QBi[1][1] * v[349] + QBi[2][1] * v[355];
//	v[388] = QBi[0][0] * v[345] + QBi[1][0] * v[349] + QBi[2][0] * v[355];
//	v[409] = v[166] * v[388] + v[167] * v[391] + v[168] * v[394];
//	v[341] = v[318] - v[313] * v[338];
//	v[379] = QBi[0][2] * v[332] + QBi[1][2] * v[337] + QBi[2][2] * v[341];
//	v[376] = QBi[0][1] * v[332] + QBi[1][1] * v[337] + QBi[2][1] * v[341];
//	v[373] = QBi[0][0] * v[332] + QBi[1][0] * v[337] + QBi[2][0] * v[341];
//	v[406] = v[166] * v[373] + v[167] * v[376] + v[168] * v[379];
//	v[330] = v[318] - v[311] * v[329];
//	v[377] = QBi[0][2] * v[330] + QBi[1][2] * v[334] + QBi[2][2] * v[339];
//	v[374] = QBi[0][1] * v[330] + QBi[1][1] * v[334] + QBi[2][1] * v[339];
//	v[371] = QBi[0][0] * v[330] + QBi[1][0] * v[334] + QBi[2][0] * v[339];
//	v[404] = v[166] * v[371] + v[167] * v[374] + v[168] * v[377];
//	v[320] = v[318] - v[311] * v[319];
//	v[362] = QBi[0][2] * v[315] + QBi[1][2] * v[320] + QBi[2][2] * v[326];
//	v[359] = QBi[0][1] * v[315] + QBi[1][1] * v[320] + QBi[2][1] * v[326];
//	v[356] = QBi[0][0] * v[315] + QBi[1][0] * v[320] + QBi[2][0] * v[326];
//	v[401] = v[166] * v[356] + v[167] * v[359] + v[168] * v[362];
//	v[434] = -(v[282] * v[401]) - v[284] * v[404] - v[286] * v[407];
//	v[422] = -(v[277] * v[401]) - v[278] * v[404] - v[279] * v[407];
//	v[317] = -(v[313] * v[314]) / 2e0 + v[336];
//	v[364] = QBi[0][2] * v[317] + QBi[1][2] * v[323] + QBi[2][2] * v[328];
//	v[361] = QBi[0][1] * v[317] + QBi[1][1] * v[323] + QBi[2][1] * v[328];
//	v[358] = QBi[0][0] * v[317] + QBi[1][0] * v[323] + QBi[2][0] * v[328];
//	v[403] = v[166] * v[358] + v[167] * v[361] + v[168] * v[364];
//	v[436] = -(v[282] * v[403]) - v[284] * v[406] - v[286] * v[409];
//	v[424] = -(v[277] * v[403]) - v[278] * v[406] - v[279] * v[409];
//	v[316] = (v[312] * v[314]) / 2e0 + v[353];
//	v[363] = QBi[0][2] * v[316] + QBi[1][2] * v[322] + QBi[2][2] * v[327];
//	v[360] = QBi[0][1] * v[316] + QBi[1][1] * v[322] + QBi[2][1] * v[327];
//	v[357] = QBi[0][0] * v[316] + QBi[1][0] * v[322] + QBi[2][0] * v[327];
//	v[402] = v[166] * v[357] + v[167] * v[360] + v[168] * v[363];
//	v[435] = -(v[282] * v[402]) - v[284] * v[405] - v[286] * v[408];
//	v[423] = -(v[277] * v[402]) - v[278] * v[405] - v[279] * v[408];
//	v[227] = (v[172] * v[172]);
//	v[2626] = v[238] / v[227];
//	v[2625] = v[236] / v[227];
//	v[2624] = v[230] / v[227];
//	v[2578] = 1e0 / Power(v[227], 3);
//	v[175] = 1e0 + (v[172] * v[314]) / 2e0;
//	v[2917] = v[175] / v[227];
//	v[1130] = v[175] * v[2624];
//	v[176] = v[172] * v[319];
//	v[2697] = v[176] * v[236];
//	v[1132] = v[176] * v[2625];
//	v[1535] = v[1130] + v[1132];
//	v[177] = v[172] * v[325];
//	v[2958] = v[177] / v[227];
//	v[2696] = v[177] * v[238];
//	v[1133] = v[177] * v[2626];
//	v[1540] = v[1130] + v[1133];
//	v[1529] = v[1132] + v[1540];
//	v[179] = v[172] * v[329];
//	v[1137] = v[179] * v[2624];
//	v[181] = 1e0 + (v[172] * v[333]) / 2e0;
//	v[2916] = v[181] / v[227];
//	v[1126] = v[181] * v[2625];
//	v[1531] = v[1126] + v[1137];
//	v[182] = v[172] * v[338];
//	v[2959] = v[182] / v[227];
//	v[1127] = v[182] * v[2626];
//	v[1541] = v[1126] + v[1127];
//	v[1534] = v[1137] + v[1541];
//	v[184] = v[172] * v[342];
//	v[1140] = v[184] * v[2624];
//	v[186] = v[172] * v[346];
//	v[2954] = v[186] / v[227];
//	v[1122] = v[186] * v[2625];
//	v[187] = 1e0 + (v[172] * v[351]) / 2e0;
//	v[2920] = v[187] / v[227];
//	v[1123] = v[187] * v[2626];
//	v[2918] = v[1123] * v[227];
//	v[1539] = v[1122] + v[1123] + v[1140];
//	v[1536] = -v[1140] + v[1539];
//	v[1530] = -v[1122] + v[1539];
//	v[245] = -(v[318] * v[324]);
//	v[2634] = v[245] - v[311];
//	v[243] = v[321] * v[324];
//	v[2632] = v[243] - v[312];
//	v[234] = -(v[318] * v[321]);
//	v[2629] = v[234] - v[313];
//	v[191] = QBi[0][0] * v[175] + QBi[1][0] * v[176] + QBi[2][0] * v[177];
//	v[192] = QBi[0][1] * v[175] + QBi[1][1] * v[176] + QBi[2][1] * v[177];
//	v[193] = QBi[0][2] * v[175] + QBi[1][2] * v[176] + QBi[2][2] * v[177];
//	v[301] = v[191] * v[296] + v[192] * v[298] + v[193] * v[300];
//	v[455] = -(v[132] * v[301]);
//	v[452] = -(v[130] * v[301]);
//	v[449] = -(v[128] * v[301]);
//	v[291] = v[191] * v[288] + v[192] * v[289] + v[193] * v[290];
//	v[443] = -(v[132] * v[291]);
//	v[440] = -(v[130] * v[291]);
//	v[437] = -(v[128] * v[291]);
//	v[194] = QBi[0][0] * v[179] + QBi[1][0] * v[181] + QBi[2][0] * v[182];
//	v[195] = QBi[0][1] * v[179] + QBi[1][1] * v[181] + QBi[2][1] * v[182];
//	v[196] = QBi[0][2] * v[179] + QBi[1][2] * v[181] + QBi[2][2] * v[182];
//	v[302] = v[194] * v[296] + v[195] * v[298] + v[196] * v[300];
//	v[456] = -(v[132] * v[302]);
//	v[453] = -(v[130] * v[302]);
//	v[450] = -(v[128] * v[302]);
//	v[292] = v[194] * v[288] + v[195] * v[289] + v[196] * v[290];
//	v[444] = -(v[132] * v[292]);
//	v[441] = -(v[130] * v[292]);
//	v[438] = -(v[128] * v[292]);
//	v[197] = QBi[0][0] * v[184] + QBi[1][0] * v[186] + QBi[2][0] * v[187];
//	v[198] = QBi[0][1] * v[184] + QBi[1][1] * v[186] + QBi[2][1] * v[187];
//	v[199] = QBi[0][2] * v[184] + QBi[1][2] * v[186] + QBi[2][2] * v[187];
//	v[303] = v[197] * v[296] + v[198] * v[298] + v[199] * v[300];
//	v[457] = -(v[132] * v[303]);
//	v[454] = -(v[130] * v[303]);
//	v[451] = -(v[128] * v[303]);
//	v[293] = v[197] * v[288] + v[198] * v[289] + v[199] * v[290];
//	v[445] = -(v[132] * v[293]);
//	v[442] = -(v[130] * v[293]);
//	v[439] = -(v[128] * v[293]);
//	v[200] = uB[0] + xBi[0];
//	v[201] = uB[1] + xBi[1];
//	v[202] = uB[2] + xBi[2];
//	v[212] = (*a6)*ddui1A[0] + (*a5)*dui1A[0] + (*a4)*u1A[0];
//	v[2637] = v[128] * v[212];
//	v[213] = (*a6)*ddui1A[1] + (*a5)*dui1A[1] + (*a4)*u1A[1];
//	v[2640] = v[128] * v[213];
//	v[214] = (*a6)*ddui1A[2] + (*a5)*dui1A[2] + (*a4)*u1A[2];
//	v[215] = (*a6)*ddui2A[0] + (*a5)*dui2A[0] + (*a4)*u2A[0];
//	v[2636] = v[130] * v[215];
//	v[216] = (*a6)*ddui2A[1] + (*a5)*dui2A[1] + (*a4)*u2A[1];
//	v[2639] = v[130] * v[216];
//	v[217] = (*a6)*ddui2A[2] + (*a5)*dui2A[2] + (*a4)*u2A[2];
//	v[218] = (*a6)*ddui3A[0] + (*a5)*dui3A[0] + (*a4)*u3A[0];
//	v[2635] = v[132] * v[218];
//	v[219] = (*a6)*ddui3A[1] + (*a5)*dui3A[1] + (*a4)*u3A[1];
//	v[2638] = v[132] * v[219];
//	v[220] = (*a6)*ddui3A[2] + (*a5)*dui3A[2] + (*a4)*u3A[2];
//	v[1411] = v[128] * v[214] + v[130] * v[217] + v[132] * v[220];
//	v[221] = (*a6)*dduiB[0] + (*a5)*duiB[0] + (*a4)*uB[0];
//	v[222] = (*a6)*dduiB[1] + (*a5)*duiB[1] + (*a4)*uB[1];
//	v[223] = (*a6)*dduiB[2] + (*a5)*duiB[2] + (*a4)*uB[2];
//	v[2398] = v[1411] - v[223];
//	v[225] = v[243] + v[312];
//	v[2633] = v[225] / v[227];
//	v[2627] = v[187] * v[225];
//	v[226] = v[234] + v[313];
//	v[2630] = v[226] / v[227];
//	v[2628] = v[181] * v[226];
//	v[228] = v[227] + (v[321] * v[321]);
//	v[2581] = -(v[238] * v[2627]) - v[236] * v[2628] - v[228] * (v[230] + v[2696] + v[2697]);
//	v[1992] = v[228] / v[227];
//	v[1222] = v[2632] / v[227];
//	v[1220] = v[2629] / v[227];
//	v[4928] = 0e0;
//	v[4929] = 0e0;
//	v[4930] = 0e0;
//	v[4931] = 0e0;
//	v[4932] = 0e0;
//	v[4933] = 0e0;
//	v[4934] = 0e0;
//	v[4935] = 0e0;
//	v[4936] = 0e0;
//	v[4937] = 0e0;
//	v[4938] = 0e0;
//	v[4939] = 0e0;
//	v[4940] = 0e0;
//	v[4941] = v[1220];
//	v[4942] = v[1222];
//	v[229] = v[1992] * v[230] + v[1220] * v[236] + v[1222] * v[238];
//	v[231] = v[227] + (v[318] * v[318]);
//	v[2685] = v[179] * v[231] + v[175] * v[2629];
//	v[237] = v[245] + v[311];
//	v[2631] = v[182] * v[231] + v[187] * v[237];
//	v[2580] = -(v[231] * v[236]) - v[238] * v[2631] - v[230] * v[2685];
//	v[1993] = v[231] / v[227];
//	v[1221] = v[2634] / v[227];
//	v[4943] = 0e0;
//	v[4944] = 0e0;
//	v[4945] = 0e0;
//	v[4946] = 0e0;
//	v[4947] = 0e0;
//	v[4948] = 0e0;
//	v[4949] = 0e0;
//	v[4950] = 0e0;
//	v[4951] = 0e0;
//	v[4952] = 0e0;
//	v[4953] = 0e0;
//	v[4954] = 0e0;
//	v[4955] = v[2630];
//	v[4956] = 0e0;
//	v[4957] = v[1221];
//	v[239] = v[1993] * v[236] + v[1221] * v[238] + v[230] * v[2630];
//	v[240] = v[227] + (v[324] * v[324]);
//	v[2686] = v[184] * v[240] + v[175] * v[2632];
//	v[2687] = v[186] * v[240] + v[181] * v[2634];
//	v[2579] = -(v[238] * v[240]) - v[230] * v[2686] - v[236] * v[2687];
//	v[1994] = v[240] / v[227];
//	v[1219] = v[237] / v[227];
//	v[4958] = 0e0;
//	v[4959] = 0e0;
//	v[4960] = 0e0;
//	v[4961] = 0e0;
//	v[4962] = 0e0;
//	v[4963] = 0e0;
//	v[4964] = 0e0;
//	v[4965] = 0e0;
//	v[4966] = 0e0;
//	v[4967] = 0e0;
//	v[4968] = 0e0;
//	v[4969] = 0e0;
//	v[4970] = v[2633];
//	v[4971] = v[1219];
//	v[4972] = 0e0;
//	v[249] = v[1219] * v[236] + v[1994] * v[238] + v[230] * v[2633];
//	v[1590] = -v[222] - v[229] * v[404] - v[239] * v[405] - v[249] * v[406];
//	v[1584] = -v[221] - v[229] * v[401] - v[239] * v[402] - v[249] * v[403];
//	v[1585] = v[1584] + v[2635] + v[2636] + v[2637];
//	v[2416] = -v[1584] + v[1585] - v[221];
//	v[2905] = v[221] + v[2416];
//	v[1579] = v[1590] + v[2638] + v[2639] + v[2640];
//	v[2407] = v[1579] - v[1590] - v[222];
//	v[2894] = v[222] + v[2407];
//	v[1574] = v[2398] - v[229] * v[407] - v[239] * v[408] - v[249] * v[409];
//	v[1591] = -v[1411] + v[1574];
//	v[250] = v[128] * v[144] + v[130] * v[145] + v[132] * v[146] - v[166] * v[191] - v[167] * v[192] - v[168] * v[193] - v[200];
//	v[413] = -v[250] / 2e0;
//	v[431] = v[132] * v[282] - v[413];
//	v[425] = v[128] * v[282] + v[413];
//	v[414] = v[130] * v[277] - v[413];
//	v[410] = v[128] * v[277] + v[413];
//	v[251] = v[128] * v[148] + v[130] * v[149] + v[132] * v[150] - v[166] * v[194] - v[167] * v[195] - v[168] * v[196] - v[201];
//	v[415] = -v[251] / 2e0;
//	v[432] = v[132] * v[284] - v[415];
//	v[426] = v[128] * v[284] + v[415];
//	v[416] = v[130] * v[278] - v[415];
//	v[411] = v[128] * v[278] + v[415];
//	v[252] = v[128] * v[152] + v[130] * v[153] + v[132] * v[154] - v[166] * v[197] - v[167] * v[198] - v[168] * v[199] - v[202];
//	v[597] = sqrt((v[250] * v[250]) + (v[251] * v[251]) + (v[252] * v[252]));
//	v[1395] = 1e0 / Power(v[597], 2);
//	v[460] = -(v[250] * (v[296] * v[358] + v[298] * v[361] + v[300] * v[364])) - v[251] * (v[296] * v[373] + v[298] * v[376]
//		+ v[300] * v[379]) - v[252] * (v[296] * v[388] + v[298] * v[391] + v[300] * v[394]) + v[301] * v[403] + v[302] * v[406]
//		+ v[303] * v[409];
//	v[459] = -(v[250] * (v[296] * v[357] + v[298] * v[360] + v[300] * v[363])) - v[251] * (v[296] * v[372] + v[298] * v[375]
//		+ v[300] * v[378]) - v[252] * (v[296] * v[387] + v[298] * v[390] + v[300] * v[393]) + v[301] * v[402] + v[302] * v[405]
//		+ v[303] * v[408];
//	v[458] = -(v[250] * (v[296] * v[356] + v[298] * v[359] + v[300] * v[362])) - v[251] * (v[296] * v[371] + v[298] * v[374]
//		+ v[300] * v[377]) - v[252] * (v[296] * v[386] + v[298] * v[389] + v[300] * v[392]) + v[301] * v[401] + v[302] * v[404]
//		+ v[303] * v[407];
//	v[448] = -(v[250] * (v[288] * v[358] + v[289] * v[361] + v[290] * v[364])) - v[251] * (v[288] * v[373] + v[289] * v[376]
//		+ v[290] * v[379]) - v[252] * (v[288] * v[388] + v[289] * v[391] + v[290] * v[394]) + v[291] * v[403] + v[292] * v[406]
//		+ v[293] * v[409];
//	v[447] = -(v[250] * (v[288] * v[357] + v[289] * v[360] + v[290] * v[363])) - v[251] * (v[288] * v[372] + v[289] * v[375]
//		+ v[290] * v[378]) - v[252] * (v[288] * v[387] + v[289] * v[390] + v[290] * v[393]) + v[291] * v[402] + v[292] * v[405]
//		+ v[293] * v[408];
//	v[446] = -(v[250] * (v[288] * v[356] + v[289] * v[359] + v[290] * v[362])) - v[251] * (v[288] * v[371] + v[289] * v[374]
//		+ v[290] * v[377]) - v[252] * (v[288] * v[386] + v[289] * v[389] + v[290] * v[392]) + v[291] * v[401] + v[292] * v[404]
//		+ v[293] * v[407];
//	v[417] = -v[252] / 2e0;
//	v[433] = v[132] * v[286] - v[417];
//	v[427] = v[128] * v[286] + v[417];
//	v[418] = v[130] * v[279] - v[417];
//	v[412] = v[128] * v[279] + v[417];
//	v[253] = -(QBi[0][0] * v[169]) - QBi[0][1] * v[170] - QBi[0][2] * v[171] - xBi[0] + v[133] * xi1A[0] + v[135] * xi2A[0]
//		+ v[137] * xi3A[0];
//	v[254] = -(QBi[1][0] * v[169]) - QBi[1][1] * v[170] - QBi[1][2] * v[171] - xBi[1] + v[133] * xi1A[1] + v[135] * xi2A[1]
//		+ v[137] * xi3A[1];
//	v[255] = -(QBi[2][0] * v[169]) - QBi[2][1] * v[170] - QBi[2][2] * v[171] - xBi[2] + v[133] * xi1A[2] + v[135] * xi2A[2]
//		+ v[137] * xi3A[2];
//	if (v[597] > 0.1e-7) { v01 = 1e0 / v[597]; v02 = (-(v01 / v[597])); v03 = (2e0*v01) / Power(v[597], 2); }
//	else {
//		v01 = (12500000e0 / 3e0)*(24e0 - (-0.1e-7 + v[597])*(0.24e10 - 2e0*(-1e0 + 100000000e0*v[597])*(0.2399999997e10
//			- 0.1199999994e18*v[597] - 0.3e17*(v[597] * v[597]))));
//		v02 = (-50000000e0 / 3e0)*(0.3599999994e10 - 0.4799999982e18*v[597] + 0.6e25*Power(v[597], 3)
//			+ 0.1799999982e26*(v[597] * v[597]));
//		v03 = 0.1e17*(799999997e0 - 0.599999994e17*v[597] - 0.3e17*(v[597] * v[597]));
//	};
//	v[262] = v03;
//	v[263] = v02;
//	v[264] = v01;
//	v[265] = v[250] * v[264];
//	v[1588] = v[1591] * v[265];
//	v[1587] = v[1590] * v[265];
//	v[1003] = -(v[265] * v[403]);
//	v[1002] = -(v[265] * v[402]);
//	v[1001] = -(v[265] * v[401]);
//	v[1000] = -(v[221] * v[265]);
//	v[577] = (v[265] * v[265]);
//	v[266] = v[251] * v[264];
//	v[2881] = v[218] * v[265] + v[219] * v[266];
//	v[2880] = v[215] * v[265] + v[216] * v[266];
//	v[2879] = v[212] * v[265] + v[213] * v[266];
//	v[2287] = (v[265] * v[266]) / 2e0;
//	v[1582] = v[1591] * v[266];
//	v[1580] = v[1584] * v[266];
//	v[1011] = -(v[266] * v[406]);
//	v[1813] = v[1003] + v[1011];
//	v[1009] = -(v[266] * v[405]);
//	v[1817] = v[1002] + v[1009];
//	v[1007] = -(v[266] * v[404]);
//	v[1821] = v[1001] + v[1007];
//	v[1006] = -(v[222] * v[266]);
//	v[1752] = v[1000] + v[1006] + v[2635] * v[265] + v[2636] * v[265] + v[2637] * v[265] + v[2638] * v[266] + v[2639] * v[266]
//		+ v[2640] * v[266];
//	v[1577] = v[1752] + v[1821] * v[229] + v[1817] * v[239] + v[1813] * v[249];
//	v[582] = (v[266] * v[266]);
//	v[581] = 2e0*v[132] * v[2287];
//	v[580] = 2e0*v[130] * v[2287];
//	v[579] = 2e0*v[128] * v[2287];
//	v[5048] = v[128] * v[577];
//	v[5049] = v[579];
//	v[5050] = 0e0;
//	v[5051] = v[130] * v[577];
//	v[5052] = v[580];
//	v[5053] = 0e0;
//	v[5054] = v[132] * v[577];
//	v[5055] = v[581];
//	v[5056] = 0e0;
//	v[5057] = -v[577];
//	v[5058] = 0e0;
//	v[5059] = 0e0;
//	v[5060] = 0e0;
//	v[5061] = 0e0;
//	v[5062] = 0e0;
//	v[5033] = v[579];
//	v[5034] = v[128] * v[582];
//	v[5035] = 0e0;
//	v[5036] = v[580];
//	v[5037] = v[130] * v[582];
//	v[5038] = 0e0;
//	v[5039] = v[581];
//	v[5040] = v[132] * v[582];
//	v[5041] = 0e0;
//	v[5042] = 0e0;
//	v[5043] = -v[582];
//	v[5044] = 0e0;
//	v[5045] = 0e0;
//	v[5046] = 0e0;
//	v[5047] = 0e0;
//	v[267] = v[252] * v[264];
//	v[2652] = v[128] * v[267];
//	v[2651] = v[130] * v[267];
//	v[2650] = v[132] * v[267];
//	v[2643] = -(v[267] * v[409]);
//	v[2645] = -v[1011] - v[2643];
//	v[2644] = v[1003] + v[2643];
//	v[2642] = -(v[267] * v[408]);
//	v[2647] = -v[1009] - v[2642];
//	v[2646] = v[1002] + v[2642];
//	v[2641] = -(v[267] * v[407]);
//	v[2649] = -v[1007] - v[2641];
//	v[2648] = v[1001] + v[2641];
//	v[4973] = 0e0;
//	v[4974] = 0e0;
//	v[4975] = v[2652];
//	v[4976] = 0e0;
//	v[4977] = 0e0;
//	v[4978] = v[2651];
//	v[4979] = 0e0;
//	v[4980] = 0e0;
//	v[4981] = v[2650];
//	v[4982] = 0e0;
//	v[4983] = 0e0;
//	v[4984] = 0e0;
//	v[4985] = 0e0;
//	v[4986] = 0e0;
//	v[4987] = 0e0;
//	v[1576] = v[1579] * v[267];
//	v[1575] = v[1585] * v[267];
//	v[1022] = v[2648] * v[266] - v[404] * v[582];
//	v[1021] = -(v[2649] * v[265]) - v[401] * v[577];
//	v[1018] = v[2646] * v[266] - v[405] * v[582];
//	v[1017] = -(v[2647] * v[265]) - v[402] * v[577];
//	v[1014] = v[2644] * v[266] - v[406] * v[582];
//	v[1013] = -(v[2645] * v[265]) - v[403] * v[577];
//	v[587] = (v[267] * v[267]);
//	v[5018] = 0e0;
//	v[5019] = 0e0;
//	v[5020] = v[128] * v[587];
//	v[5021] = 0e0;
//	v[5022] = 0e0;
//	v[5023] = v[130] * v[587];
//	v[5024] = 0e0;
//	v[5025] = 0e0;
//	v[5026] = v[132] * v[587];
//	v[5027] = 0e0;
//	v[5028] = 0e0;
//	v[5029] = -v[587];
//	v[5030] = 0e0;
//	v[5031] = 0e0;
//	v[5032] = 0e0;
//	v[1023] = v[1821] * v[267] - v[407] * v[587];
//	v[1019] = v[1817] * v[267] - v[408] * v[587];
//	v[1015] = v[1813] * v[267] - v[409] * v[587];
//	v[1005] = v[220] * v[2650] + v[217] * v[2651] + v[214] * v[2652] - v[223] * v[267];
//	v[2415] = v[1005] + v[1006];
//	v[2406] = v[1000] + v[1005];
//	v[1586] = v[2415] - v[249] * v[2645] - v[239] * v[2647] - v[229] * v[2649];
//	v[1581] = v[2406] + v[249] * v[2644] + v[239] * v[2646] + v[229] * v[2648];
//	v[268] = sqrt((v[253] * v[253]) + (v[254] * v[254]) + (v[255] * v[255]));
//	if (v[268] > 0.1e-7) { v04 = 1e0 / v[268]; v05 = (-(v04 / v[268])); v06 = (2e0*v04) / Power(v[268], 2); }
//	else {
//		v04 = (12500000e0 / 3e0)*(24e0 - (-0.1e-7 + v[268])*(0.24e10 - 2e0*(-1e0 + 100000000e0*v[268])*(0.2399999997e10
//			- 0.1199999994e18*v[268] - 0.3e17*(v[268] * v[268]))));
//		v05 = (-50000000e0 / 3e0)*(0.3599999994e10 - 0.4799999982e18*v[268] + 0.6e25*Power(v[268], 3)
//			+ 0.1799999982e26*(v[268] * v[268]));
//		v06 = 0.1e17*(799999997e0 - 0.599999994e17*v[268] - 0.3e17*(v[268] * v[268]));
//	};
//	v[273] = v04;
//	v[274] = v[253] * v[273];
//	v[275] = v[254] * v[273];
//	v[276] = v[255] * v[273];
//	v[477] = -(invH[0][0] * v[410]) - invH[0][1] * v[425] - invH[0][2] * v[437] - invH[0][3] * v[449];
//	v[478] = -(invH[0][0] * v[411]) - invH[0][1] * v[426] - invH[0][2] * v[438] - invH[0][3] * v[450];
//	v[479] = -(invH[0][0] * v[412]) - invH[0][1] * v[427] - invH[0][2] * v[439] - invH[0][3] * v[451];
//	v[480] = -(invH[0][0] * v[414]) - invH[0][1] * v[428] - invH[0][2] * v[440] - invH[0][3] * v[452];
//	v[481] = -(invH[0][0] * v[416]) - invH[0][1] * v[429] - invH[0][2] * v[441] - invH[0][3] * v[453];
//	v[482] = -(invH[0][0] * v[418]) - invH[0][1] * v[430] - invH[0][2] * v[442] - invH[0][3] * v[454];
//	v[483] = -(invH[0][0] * v[419]) - invH[0][1] * v[431] - invH[0][2] * v[443] - invH[0][3] * v[455];
//	v[484] = -(invH[0][0] * v[420]) - invH[0][1] * v[432] - invH[0][2] * v[444] - invH[0][3] * v[456];
//	v[485] = -(invH[0][0] * v[421]) - invH[0][1] * v[433] - invH[0][2] * v[445] - invH[0][3] * v[457];
//	v[486] = invH[0][0] * v[277] + invH[0][1] * v[282] - invH[0][2] * v[291] - invH[0][3] * v[301];
//	v[487] = invH[0][0] * v[278] + invH[0][1] * v[284] - invH[0][2] * v[292] - invH[0][3] * v[302];
//	v[488] = invH[0][0] * v[279] + invH[0][1] * v[286] - invH[0][2] * v[293] - invH[0][3] * v[303];
//	v[489] = -(invH[0][0] * v[422]) - invH[0][1] * v[434] - invH[0][2] * v[446] - invH[0][3] * v[458];
//	v[490] = -(invH[0][0] * v[423]) - invH[0][1] * v[435] - invH[0][2] * v[447] - invH[0][3] * v[459];
//	v[491] = -(invH[0][0] * v[424]) - invH[0][1] * v[436] - invH[0][2] * v[448] - invH[0][3] * v[460];
//	v[955] = -(v[212] * v[477]) - v[213] * v[478] - v[214] * v[479] - v[215] * v[480] - v[216] * v[481] - v[217] * v[482]
//		- v[218] * v[483] - v[219] * v[484] - v[220] * v[485] - v[221] * v[486] - v[222] * v[487] - v[223] * v[488] - v[229] * v[489]
//		- v[239] * v[490] - v[249] * v[491];
//	v[4121] = v[477];
//	v[4122] = v[478];
//	v[4123] = v[479];
//	v[4124] = v[480];
//	v[4125] = v[481];
//	v[4126] = v[482];
//	v[4127] = v[483];
//	v[4128] = v[484];
//	v[4129] = v[485];
//	v[4130] = v[486];
//	v[4131] = v[487];
//	v[4132] = v[488];
//	v[4133] = v[489];
//	v[4134] = v[490];
//	v[4135] = v[491];
//	v[492] = -(invH[1][0] * v[410]) - invH[1][1] * v[425] - invH[1][2] * v[437] - invH[1][3] * v[449];
//	v[493] = -(invH[1][0] * v[411]) - invH[1][1] * v[426] - invH[1][2] * v[438] - invH[1][3] * v[450];
//	v[494] = -(invH[1][0] * v[412]) - invH[1][1] * v[427] - invH[1][2] * v[439] - invH[1][3] * v[451];
//	v[495] = -(invH[1][0] * v[414]) - invH[1][1] * v[428] - invH[1][2] * v[440] - invH[1][3] * v[452];
//	v[496] = -(invH[1][0] * v[416]) - invH[1][1] * v[429] - invH[1][2] * v[441] - invH[1][3] * v[453];
//	v[497] = -(invH[1][0] * v[418]) - invH[1][1] * v[430] - invH[1][2] * v[442] - invH[1][3] * v[454];
//	v[498] = -(invH[1][0] * v[419]) - invH[1][1] * v[431] - invH[1][2] * v[443] - invH[1][3] * v[455];
//	v[499] = -(invH[1][0] * v[420]) - invH[1][1] * v[432] - invH[1][2] * v[444] - invH[1][3] * v[456];
//	v[500] = -(invH[1][0] * v[421]) - invH[1][1] * v[433] - invH[1][2] * v[445] - invH[1][3] * v[457];
//	v[501] = invH[1][0] * v[277] + invH[1][1] * v[282] - invH[1][2] * v[291] - invH[1][3] * v[301];
//	v[502] = invH[1][0] * v[278] + invH[1][1] * v[284] - invH[1][2] * v[292] - invH[1][3] * v[302];
//	v[503] = invH[1][0] * v[279] + invH[1][1] * v[286] - invH[1][2] * v[293] - invH[1][3] * v[303];
//	v[504] = -(invH[1][0] * v[422]) - invH[1][1] * v[434] - invH[1][2] * v[446] - invH[1][3] * v[458];
//	v[505] = -(invH[1][0] * v[423]) - invH[1][1] * v[435] - invH[1][2] * v[447] - invH[1][3] * v[459];
//	v[506] = -(invH[1][0] * v[424]) - invH[1][1] * v[436] - invH[1][2] * v[448] - invH[1][3] * v[460];
//	v[953] = -(v[212] * v[492]) - v[213] * v[493] - v[214] * v[494] - v[215] * v[495] - v[216] * v[496] - v[217] * v[497]
//		- v[218] * v[498] - v[219] * v[499] - v[220] * v[500] - v[221] * v[501] - v[222] * v[502] - v[223] * v[503] - v[229] * v[504]
//		- v[239] * v[505] - v[249] * v[506];
//	v[4106] = v[492];
//	v[4107] = v[493];
//	v[4108] = v[494];
//	v[4109] = v[495];
//	v[4110] = v[496];
//	v[4111] = v[497];
//	v[4112] = v[498];
//	v[4113] = v[499];
//	v[4114] = v[500];
//	v[4115] = v[501];
//	v[4116] = v[502];
//	v[4117] = v[503];
//	v[4118] = v[504];
//	v[4119] = v[505];
//	v[4120] = v[506];
//	v[507] = -(invH[2][0] * v[410]) - invH[2][1] * v[425] - invH[2][2] * v[437] - invH[2][3] * v[449];
//	v[508] = -(invH[2][0] * v[411]) - invH[2][1] * v[426] - invH[2][2] * v[438] - invH[2][3] * v[450];
//	v[509] = -(invH[2][0] * v[412]) - invH[2][1] * v[427] - invH[2][2] * v[439] - invH[2][3] * v[451];
//	v[510] = -(invH[2][0] * v[414]) - invH[2][1] * v[428] - invH[2][2] * v[440] - invH[2][3] * v[452];
//	v[511] = -(invH[2][0] * v[416]) - invH[2][1] * v[429] - invH[2][2] * v[441] - invH[2][3] * v[453];
//	v[512] = -(invH[2][0] * v[418]) - invH[2][1] * v[430] - invH[2][2] * v[442] - invH[2][3] * v[454];
//	v[513] = -(invH[2][0] * v[419]) - invH[2][1] * v[431] - invH[2][2] * v[443] - invH[2][3] * v[455];
//	v[514] = -(invH[2][0] * v[420]) - invH[2][1] * v[432] - invH[2][2] * v[444] - invH[2][3] * v[456];
//	v[515] = -(invH[2][0] * v[421]) - invH[2][1] * v[433] - invH[2][2] * v[445] - invH[2][3] * v[457];
//	v[516] = invH[2][0] * v[277] + invH[2][1] * v[282] - invH[2][2] * v[291] - invH[2][3] * v[301];
//	v[517] = invH[2][0] * v[278] + invH[2][1] * v[284] - invH[2][2] * v[292] - invH[2][3] * v[302];
//	v[518] = invH[2][0] * v[279] + invH[2][1] * v[286] - invH[2][2] * v[293] - invH[2][3] * v[303];
//	v[519] = -(invH[2][0] * v[422]) - invH[2][1] * v[434] - invH[2][2] * v[446] - invH[2][3] * v[458];
//	v[520] = -(invH[2][0] * v[423]) - invH[2][1] * v[435] - invH[2][2] * v[447] - invH[2][3] * v[459];
//	v[521] = -(invH[2][0] * v[424]) - invH[2][1] * v[436] - invH[2][2] * v[448] - invH[2][3] * v[460];
//	v[949] = v[212] * v[507] + v[213] * v[508] + v[214] * v[509] + v[215] * v[510] + v[216] * v[511] + v[217] * v[512]
//		+ v[218] * v[513] + v[219] * v[514] + v[220] * v[515] + v[221] * v[516] + v[222] * v[517] + v[223] * v[518] + v[229] * v[519]
//		+ v[239] * v[520] + v[249] * v[521];
//	v[4136] = v[507];
//	v[4137] = v[508];
//	v[4138] = v[509];
//	v[4139] = v[510];
//	v[4140] = v[511];
//	v[4141] = v[512];
//	v[4142] = v[513];
//	v[4143] = v[514];
//	v[4144] = v[515];
//	v[4145] = v[516];
//	v[4146] = v[517];
//	v[4147] = v[518];
//	v[4148] = v[519];
//	v[4149] = v[520];
//	v[4150] = v[521];
//	v[522] = -(invH[3][0] * v[410]) - invH[3][1] * v[425] - invH[3][2] * v[437] - invH[3][3] * v[449];
//	v[874] = v[279] * v[477] + v[286] * v[492] - v[293] * v[507] - v[303] * v[522];
//	v[873] = v[278] * v[477] + v[284] * v[492] - v[292] * v[507] - v[302] * v[522];
//	v[872] = v[277] * v[477] + v[282] * v[492] - v[291] * v[507] - v[301] * v[522];
//	v[523] = -(invH[3][0] * v[411]) - invH[3][1] * v[426] - invH[3][2] * v[438] - invH[3][3] * v[450];
//	v[878] = v[279] * v[478] + v[286] * v[493] - v[293] * v[508] - v[303] * v[523];
//	v[877] = v[278] * v[478] + v[284] * v[493] - v[292] * v[508] - v[302] * v[523];
//	v[876] = v[277] * v[478] + v[282] * v[493] - v[291] * v[508] - v[301] * v[523];
//	v[524] = -(invH[3][0] * v[412]) - invH[3][1] * v[427] - invH[3][2] * v[439] - invH[3][3] * v[451];
//	v[882] = v[279] * v[479] + v[286] * v[494] - v[293] * v[509] - v[303] * v[524];
//	v[881] = v[278] * v[479] + v[284] * v[494] - v[292] * v[509] - v[302] * v[524];
//	v[880] = v[277] * v[479] + v[282] * v[494] - v[291] * v[509] - v[301] * v[524];
//	v[525] = -(invH[3][0] * v[414]) - invH[3][1] * v[428] - invH[3][2] * v[440] - invH[3][3] * v[452];
//	v[886] = v[279] * v[480] + v[286] * v[495] - v[293] * v[510] - v[303] * v[525];
//	v[885] = v[278] * v[480] + v[284] * v[495] - v[292] * v[510] - v[302] * v[525];
//	v[884] = v[277] * v[480] + v[282] * v[495] - v[291] * v[510] - v[301] * v[525];
//	v[526] = -(invH[3][0] * v[416]) - invH[3][1] * v[429] - invH[3][2] * v[441] - invH[3][3] * v[453];
//	v[890] = v[279] * v[481] + v[286] * v[496] - v[293] * v[511] - v[303] * v[526];
//	v[889] = v[278] * v[481] + v[284] * v[496] - v[292] * v[511] - v[302] * v[526];
//	v[888] = v[277] * v[481] + v[282] * v[496] - v[291] * v[511] - v[301] * v[526];
//	v[527] = -(invH[3][0] * v[418]) - invH[3][1] * v[430] - invH[3][2] * v[442] - invH[3][3] * v[454];
//	v[894] = v[279] * v[482] + v[286] * v[497] - v[293] * v[512] - v[303] * v[527];
//	v[893] = v[278] * v[482] + v[284] * v[497] - v[292] * v[512] - v[302] * v[527];
//	v[892] = v[277] * v[482] + v[282] * v[497] - v[291] * v[512] - v[301] * v[527];
//	v[528] = -(invH[3][0] * v[419]) - invH[3][1] * v[431] - invH[3][2] * v[443] - invH[3][3] * v[455];
//	v[898] = v[279] * v[483] + v[286] * v[498] - v[293] * v[513] - v[303] * v[528];
//	v[897] = v[278] * v[483] + v[284] * v[498] - v[292] * v[513] - v[302] * v[528];
//	v[896] = v[277] * v[483] + v[282] * v[498] - v[291] * v[513] - v[301] * v[528];
//	v[529] = -(invH[3][0] * v[420]) - invH[3][1] * v[432] - invH[3][2] * v[444] - invH[3][3] * v[456];
//	v[902] = v[279] * v[484] + v[286] * v[499] - v[293] * v[514] - v[303] * v[529];
//	v[901] = v[278] * v[484] + v[284] * v[499] - v[292] * v[514] - v[302] * v[529];
//	v[900] = v[277] * v[484] + v[282] * v[499] - v[291] * v[514] - v[301] * v[529];
//	v[530] = -(invH[3][0] * v[421]) - invH[3][1] * v[433] - invH[3][2] * v[445] - invH[3][3] * v[457];
//	v[906] = v[279] * v[485] + v[286] * v[500] - v[293] * v[515] - v[303] * v[530];
//	v[905] = v[278] * v[485] + v[284] * v[500] - v[292] * v[515] - v[302] * v[530];
//	v[904] = v[277] * v[485] + v[282] * v[500] - v[291] * v[515] - v[301] * v[530];
//	v[531] = invH[3][0] * v[277] + invH[3][1] * v[282] - invH[3][2] * v[291] - invH[3][3] * v[301];
//	v[910] = v[279] * v[486] + v[286] * v[501] - v[293] * v[516] - v[303] * v[531];
//	v[909] = v[278] * v[486] + v[284] * v[501] - v[292] * v[516] - v[302] * v[531];
//	v[908] = v[277] * v[486] + v[282] * v[501] - v[291] * v[516] - v[301] * v[531];
//	v[532] = invH[3][0] * v[278] + invH[3][1] * v[284] - invH[3][2] * v[292] - invH[3][3] * v[302];
//	v[914] = v[279] * v[487] + v[286] * v[502] - v[293] * v[517] - v[303] * v[532];
//	v[913] = v[278] * v[487] + v[284] * v[502] - v[292] * v[517] - v[302] * v[532];
//	v[912] = v[277] * v[487] + v[282] * v[502] - v[291] * v[517] - v[301] * v[532];
//	v[533] = invH[3][0] * v[279] + invH[3][1] * v[286] - invH[3][2] * v[293] - invH[3][3] * v[303];
//	v[918] = v[279] * v[488] + v[286] * v[503] - v[293] * v[518] - v[303] * v[533];
//	v[917] = v[278] * v[488] + v[284] * v[503] - v[292] * v[518] - v[302] * v[533];
//	v[916] = v[277] * v[488] + v[282] * v[503] - v[291] * v[518] - v[301] * v[533];
//	v[534] = -(invH[3][0] * v[422]) - invH[3][1] * v[434] - invH[3][2] * v[446] - invH[3][3] * v[458];
//	v[922] = v[279] * v[489] + v[286] * v[504] - v[293] * v[519] - v[303] * v[534];
//	v[921] = v[278] * v[489] + v[284] * v[504] - v[292] * v[519] - v[302] * v[534];
//	v[920] = v[277] * v[489] + v[282] * v[504] - v[291] * v[519] - v[301] * v[534];
//	v[535] = -(invH[3][0] * v[423]) - invH[3][1] * v[435] - invH[3][2] * v[447] - invH[3][3] * v[459];
//	v[926] = v[279] * v[490] + v[286] * v[505] - v[293] * v[520] - v[303] * v[535];
//	v[925] = v[278] * v[490] + v[284] * v[505] - v[292] * v[520] - v[302] * v[535];
//	v[924] = v[277] * v[490] + v[282] * v[505] - v[291] * v[520] - v[301] * v[535];
//	v[536] = -(invH[3][0] * v[424]) - invH[3][1] * v[436] - invH[3][2] * v[448] - invH[3][3] * v[460];
//	v[951] = v[212] * v[522] + v[213] * v[523] + v[214] * v[524] + v[215] * v[525] + v[216] * v[526] + v[217] * v[527]
//		+ v[218] * v[528] + v[219] * v[529] + v[220] * v[530] + v[221] * v[531] + v[222] * v[532] + v[223] * v[533] + v[229] * v[534]
//		+ v[239] * v[535] + v[249] * v[536];
//	v[930] = v[279] * v[491] + v[286] * v[506] - v[293] * v[521] - v[303] * v[536];
//	v[929] = v[278] * v[491] + v[284] * v[506] - v[292] * v[521] - v[302] * v[536];
//	v[928] = v[277] * v[491] + v[282] * v[506] - v[291] * v[521] - v[301] * v[536];
//	v[4151] = v[522];
//	v[4152] = v[523];
//	v[4153] = v[524];
//	v[4154] = v[525];
//	v[4155] = v[526];
//	v[4156] = v[527];
//	v[4157] = v[528];
//	v[4158] = v[529];
//	v[4159] = v[530];
//	v[4160] = v[531];
//	v[4161] = v[532];
//	v[4162] = v[533];
//	v[4163] = v[534];
//	v[4164] = v[535];
//	v[4165] = v[536];
//	b537 = sqrt(Power(v[266] * v[274] - v[265] * v[275], 2) + Power(-(v[267] * v[274]) + v[265] * v[276], 2) + Power
//	(v[267] * v[275] - v[266] * v[276], 2)) > 0.1e-7;
//	if (b537) {
//		v[539] = v[267] * v[275] - v[266] * v[276];
//		v[540] = -(v[267] * v[274]) + v[265] * v[276];
//		v[541] = v[266] * v[274] - v[265] * v[275];
//		v[542] = sqrt((v[539] * v[539]) + (v[540] * v[540]) + (v[541] * v[541]));
//		v[1455] = 1e0 / Power(v[542], 2);
//		v[1109] = v[542];
//		v[1466] = 1e0 - (v[1109] * v[1109]);
//		v[3015] = 1e0 / Power(v[1466], 0.15e1);
//		v[1461] = 1e0 / sqrt(v[1466]);
//		v[1108] = asin(v[1109]) / 2e0;
//		v[3014] = tan(v[1108]);
//		v[1460] = 1e0 / Power(cos(v[1108]), 2);
//		v[2724] = v[1460] * v[1461];
//		v[544] = 2e0*tan(v[1108]);
//		if (v[542] > 0.1e-7) { v07 = 1e0 / v[542]; v08 = (-(v07 / v[542])); v09 = (2e0*v07) / Power(v[542], 2); }
//		else {
//			v07 = (12500000e0 / 3e0)*(24e0 - (-0.1e-7 + v[542])*(0.24e10 - 2e0*(-1e0 + 100000000e0*v[542])*
//				(0.2399999997e10 - 0.1199999994e18*v[542] - 0.3e17*(v[542] * v[542]))));
//			v08 = (-50000000e0 / 3e0)*(0.3599999994e10 - 0.4799999982e18*v[542] + 0.6e25*Power(v[542], 3)
//				+ 0.1799999982e26*(v[542] * v[542]));
//			v09 = 0.1e17*(799999997e0 - 0.599999994e17*v[542] - 0.3e17*(v[542] * v[542]));
//		};
//		v[548] = v09;
//		v[549] = v08;
//		v[550] = v07;
//		v[3013] = v[544] * v[549] + v[2724] * v[550];
//		v[2653] = v[544] * v[550];
//		v[551] = v[2653] * v[539];
//		v[562] = (v[551] * v[551]);
//		v[552] = v[2653] * v[540];
//		v[560] = (v[551] * v[552]) / 2e0;
//		v[555] = (v[552] * v[552]);
//		v[1092] = -v[555] - v[562];
//		v[553] = v[2653] * v[541];
//		v[1087] = v[553] + v[560];
//		v[1085] = -v[553] + v[560];
//		v[567] = (v[552] * v[553]) / 2e0;
//		v[1091] = v[551] + v[567];
//		v[1089] = -v[551] + v[567];
//		v[565] = (v[551] * v[553]) / 2e0;
//		v[1090] = -v[552] + v[565];
//		v[1086] = v[552] + v[565];
//		v[556] = (v[553] * v[553]);
//		v[1096] = 4e0 + v[555] + v[556] + v[562];
//		v[3016] = 1e0 / Power(v[1096], 3);
//		v[1789] = 1e0 / Power(v[1096], 2);
//		v[1088] = -v[556] - v[562];
//		v[1084] = -v[555] - v[556];
//		v[554] = 4e0 / v[1096];
//		v[557] = 1e0 + (v[1084] * v[554]) / 2e0;
//		v[558] = v[1085] * v[554];
//		v[559] = v[1086] * v[554];
//		v[561] = v[1087] * v[554];
//		v[563] = 1e0 + (v[1088] * v[554]) / 2e0;
//		v[564] = v[1089] * v[554];
//		v[566] = v[1090] * v[554];
//		v[568] = v[1091] * v[554];
//		v[569] = 1e0 + (v[1092] * v[554]) / 2e0;
//	}
//	else {
//		v[557] = 1e0;
//		v[558] = 0e0;
//		v[559] = 0e0;
//		v[561] = 0e0;
//		v[563] = 1e0;
//		v[564] = 0e0;
//		v[566] = 0e0;
//		v[568] = 0e0;
//		v[569] = 1e0;
//	};
//	if ((*previouscontact)) {
//		v[1054] = 1e0 - v[587];
//		v[1052] = 1e0 - v[582];
//		v[1050] = 1e0 - v[577];
//		v[574] = v[133] * v[152] + v[135] * v[153] + v[137] * v[154] - v[169] * v[197] - v[170] * v[198] - v[171] * v[199] - v[202]
//			+ gti[0] * v[566] + gti[1] * v[568] + gti[2] * v[569];
//		v[2654] = v[267] * v[574];
//		v[573] = v[133] * v[148] + v[135] * v[149] + v[137] * v[150] - v[169] * v[194] - v[170] * v[195] - v[171] * v[196] - v[201]
//			+ gti[0] * v[561] + gti[1] * v[563] + gti[2] * v[564];
//		v[2656] = v[266] * v[573];
//		v[2678] = v[2654] + v[2656];
//		v[572] = v[133] * v[144] + v[135] * v[145] + v[137] * v[146] - v[169] * v[191] - v[170] * v[192] - v[171] * v[193] - v[200]
//			+ gti[0] * v[557] + gti[1] * v[558] + gti[2] * v[559];
//		v[2655] = -(v[265] * v[572]);
//		v[2677] = -v[2654] + v[2655];
//		v[2676] = v[2655] - v[2656];
//		v[571] = -(v[265] * v[2678]) + v[1050] * v[572];
//		v[575] = v[266] * v[2677] + v[1052] * v[573];
//		v[576] = v[267] * v[2676] + v[1054] * v[574];
//	}
//	else {
//		v[571] = 0e0;
//		v[575] = 0e0;
//		v[576] = 0e0;
//	};
//	v[578] = v[1021] * v[229] + v[1017] * v[239] + v[1013] * v[249] + v[2415] * v[265] + v[2416] * v[577] + v[213] * v[579]
//		+ v[216] * v[580] + v[219] * v[581];
//	v[586] = v[1022] * v[229] + v[1018] * v[239] + v[1014] * v[249] + v[2406] * v[266] + v[212] * v[579] + v[215] * v[580]
//		+ v[218] * v[581] + v[2407] * v[582];
//	v[588] = v[1023] * v[229] + v[1019] * v[239] + v[1015] * v[249] + v[1752] * v[267] + v[2398] * v[587];
//	(*vnrel) = sqrt((v[578] * v[578]) + (v[586] * v[586]) + (v[588] * v[588]));
//	v[594] = -((v[2658] * Power(v[595], v[610])) / ((*n2)*Power((*gnbb), v[615])));
//	b598 = v[597] < (*gnb);
//	if (b598) {
//		b600 = v[597] > (*gnbb);
//		if (b600) {
//			v[609] = (*gnb) - v[597];
//			v[3188] = Power(v[609], v[610]);
//			v[3183] = Power(v[609], v[1602]);
//			v[3090] = Power(v[609], v[610]);
//			v[3047] = Power(v[609], -2e0 + v[610]);
//			v[1724] = Power(v[609], v[1602]);
//			v[2970] = Power(v[609], v[610]);
//			v[603] = Power(v[609], (*n1));
//			v[2884] = -((*epsn1)*v[603]);
//			v[2767] = -((*epsn1)*v[603]);
//			v[2663] = -((*epsn1)*v[603]);
//			v[2657] = -((*epsn1)*v[603]);
//			v[602] = v[265] * v[2657];
//			v[604] = v[2657] * v[266];
//			v[605] = v[2657] * v[267];
//		}
//		else {
//			v[3192] = Power(v[597], v[615]);
//			v[3184] = Power(v[597], v[1607]);
//			v[3094] = Power(v[597], v[615]);
//			v[3052] = Power(v[597], -2e0 + v[615]);
//			v[1730] = Power(v[597], v[1607]);
//			v[2974] = Power(v[597], v[615]);
//			v[606] = -((*epsn1)*Power(v[595], (*n1))) + v[594] * (Power((*gnbb), (*n2)) - Power(v[597], (*n2)));
//			v[602] = v[265] * v[606];
//			v[604] = v[266] * v[606];
//			v[605] = v[267] * v[606];
//		};
//	}
//	else {
//		v[602] = 0e0;
//		v[604] = 0e0;
//		v[605] = 0e0;
//	};
//	v[2668] = v[605] / 2e0;
//	v[2667] = v[604] / 2e0;
//	v[4425] = 0e0;
//	v[4426] = 0e0;
//	v[4427] = 0e0;
//	v[4428] = 0e0;
//	v[4429] = 0e0;
//	v[4430] = 0e0;
//	v[4431] = v[602];
//	v[4432] = v[604];
//	v[4433] = v[605];
//	v[4434] = 0e0;
//	v[4435] = 0e0;
//	v[4436] = 0e0;
//	v[4437] = 0e0;
//	v[4438] = 0e0;
//	v[4439] = 0e0;
//	v[4440] = 0e0;
//	v[4441] = 0e0;
//	v[4442] = 0e0;
//	v[4443] = v[602];
//	v[4444] = v[604];
//	v[4445] = v[605];
//	v[4446] = 0e0;
//	v[4447] = 0e0;
//	v[4448] = 0e0;
//	v[4449] = 0e0;
//	v[4450] = 0e0;
//	v[4451] = 0e0;
//	v[4452] = 0e0;
//	v[4453] = 0e0;
//	v[4454] = 0e0;
//	v[4455] = v[602];
//	v[4456] = v[604];
//	v[4457] = v[605];
//	v[4458] = 0e0;
//	v[4459] = 0e0;
//	v[4460] = 0e0;
//	v[4461] = 0e0;
//	v[4462] = 0e0;
//	v[4463] = 0e0;
//	v[4464] = 0e0;
//	v[4465] = 0e0;
//	v[4466] = 0e0;
//	v[4467] = 0e0;
//	v[4468] = 0e0;
//	v[4469] = 0e0;
//	v[2666] = v[602] / 2e0;
//	if (b598) {
//		if (b600) {
//			v[975] = (*meq)*v[2658];
//			v[976] = sqrt(v[975] * Power(v[609], v[610]));
//			v[1723] = (v[610] * v[975]) / v[976];
//			v[612] = 2e0*v[976] * (*zetan);
//			v[611] = v[578] * v[612];
//			v[613] = v[586] * v[612];
//			v[614] = v[588] * v[612];
//		}
//		else {
//			v[980] = -((*meq)*(*n2)*v[594]);
//			v[981] = sqrt(v[980] * Power(v[597], v[615]));
//			v[1729] = (v[615] * v[980]) / v[981];
//			v[616] = 2e0*v[981] * (*zetan);
//			v[611] = v[578] * v[616];
//			v[613] = v[586] * v[616];
//			v[614] = v[588] * v[616];
//		};
//		b617 = v[597] < (*gnb) && v[265] * (v[602] + v[611]) + v[266] * (v[604] + v[613]) + v[267] * (v[605] + v[614]) < 0e0;
//		if (b617) {
//			v[619] = v[611];
//			v[620] = v[613];
//			v[621] = v[614];
//		}
//		else {
//			v[619] = -v[602];
//			v[620] = -v[604];
//			v[621] = -v[605];
//		};
//	}
//	else {
//		v[619] = 0e0;
//		v[620] = 0e0;
//		v[621] = 0e0;
//	};
//	v[622] = v[602] + v[619];
//	v[623] = v[604] + v[620];
//	v[624] = v[605] + v[621];
//	v[1644] = (v[622] * v[622]) + (v[623] * v[623]) + (v[624] * v[624]);
//	v[625] = (*epst)*v[571];
//	v[626] = (*epst)*v[575];
//	v[627] = (*epst)*v[576];
//	v[631] = v[625] - (*ct)*(v[212] * v[872] + v[213] * v[876] + v[214] * v[880] + v[215] * v[884] + v[216] * v[888]
//		+ v[217] * v[892] + v[218] * v[896] + v[219] * v[900] + v[220] * v[904] + v[221] * v[908] + v[222] * v[912] + v[223] * v[916]
//		+ v[229] * v[920] + v[239] * v[924] + v[249] * v[928]);
//	v[632] = v[626] - (*ct)*(v[212] * v[873] + v[213] * v[877] + v[214] * v[881] + v[215] * v[885] + v[216] * v[889]
//		+ v[217] * v[893] + v[218] * v[897] + v[219] * v[901] + v[220] * v[905] + v[221] * v[909] + v[222] * v[913] + v[223] * v[917]
//		+ v[229] * v[921] + v[239] * v[925] + v[249] * v[929]);
//	v[633] = v[627] - (*ct)*(v[212] * v[874] + v[213] * v[878] + v[214] * v[882] + v[215] * v[886] + v[216] * v[890]
//		+ v[217] * v[894] + v[218] * v[898] + v[219] * v[902] + v[220] * v[906] + v[221] * v[910] + v[222] * v[914] + v[223] * v[918]
//		+ v[229] * v[922] + v[239] * v[926] + v[249] * v[930]);
//	v[1641] = (v[631] * v[631]) + (v[632] * v[632]) + (v[633] * v[633]);
//	if (b598) {
//		if ((*stick)) {
//			b636 = sqrt((v[631] * v[631]) + (v[632] * v[632]) + (v[633] * v[633])) <= (*mus)*sqrt((v[622] * v[622]) +
//				(v[623] * v[623]) + (v[624] * v[624]));
//			if (b636) {
//				v[638] = v[631];
//				v[639] = v[632];
//				v[640] = v[633];
//				v[641] = 1e0;
//			}
//			else {
//				v[1643] = 1e0 / sqrt(v[1641]);
//				v[1645] = 1e0 / sqrt(v[1644]);
//				v[651] = sqrt(v[1644]);
//				v[642] = sqrt(v[1641]);
//				if (v[642] > 0.1e-5) { v010 = 1e0 / v[642]; v011 = (-(v010 / v[642])); v012 = (2e0*v010) / Power(v[642], 2); }
//				else {
//					v010 = (24000000e0 - (-1e0 + 1000000e0*v[642])*(71999994e0 - 0.71999982e14*v[642] + 0.6e19*Power(v[642], 3)
//						+ 0.23999982e20*(v[642] * v[642]))) / 24e0;
//					v011 = (-500000e0 / 3e0)*(35999994e0 - 0.47999982e14*v[642] + 0.6e19*Power(v[642], 3) + 0.17999982e20*
//						(v[642] * v[642]));
//					v012 = 0.1e13*(7999997e0 - 0.5999994e13*v[642] - 0.3e13*(v[642] * v[642]));
//				};
//				v[646] = v011;
//				v[647] = v010;
//				v[1642] = (*mud)*v[647] * v[651];
//				v[638] = v[1642] * v[631];
//				v[639] = v[1642] * v[632];
//				v[640] = v[1642] * v[633];
//				v[641] = 0e0;
//			};
//			if (sqrt((v[625] * v[625]) + (v[626] * v[626]) + (v[627] * v[627])) > (*mus)*sqrt((v[622] * v[622]) +
//				(v[623] * v[623]) + (v[624] * v[624]))) {
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
//				v[660] = sqrt((v[625] * v[625]) + (v[626] * v[626]) + (v[627] * v[627]));
//				if (v[660] > 0.1e-5) { v016 = 1e0 / v[660]; v017 = (-(v016 / v[660])); v018 = (2e0*v016) / Power(v[660], 2); }
//				else {
//					v016 = (24000000e0 - (-1e0 + 1000000e0*v[660])*(71999994e0 - 0.71999982e14*v[660] + 0.6e19*Power(v[660], 3)
//						+ 0.23999982e20*(v[660] * v[660]))) / 24e0;
//					v017 = (-500000e0 / 3e0)*(35999994e0 - 0.47999982e14*v[660] + 0.6e19*Power(v[660], 3) + 0.17999982e20*
//						(v[660] * v[660]));
//					v018 = 0.1e13*(7999997e0 - 0.5999994e13*v[660] - 0.3e13*(v[660] * v[660]));
//				};
//				v[667] = -((*mud)*v013*v016*sqrt(v[1644]));
//				v[666] = v[571] + v[625] * v[667];
//				v[668] = v[575] + v[626] * v[667];
//				v[669] = v[576] + v[627] * v[667];
//			}
//			else {
//				v[666] = 0e0;
//				v[668] = 0e0;
//				v[669] = 0e0;
//			};
//		}
//		else {
//			b670 = sqrt((v[631] * v[631]) + (v[632] * v[632]) + (v[633] * v[633])) <= (*mud)*sqrt((v[622] * v[622]) +
//				(v[623] * v[623]) + (v[624] * v[624]));
//			if (b670) {
//				v[638] = v[631];
//				v[639] = v[632];
//				v[640] = v[633];
//				v[641] = 1e0;
//			}
//			else {
//				v[1651] = 1e0 / sqrt(v[1641]);
//				v[1655] = 1e0 / sqrt(v[1644]);
//				v[681] = sqrt(v[1644]);
//				v[2758] = (*mud)*v[681];
//				v[672] = sqrt(v[1641]);
//				if (v[672] > 0.1e-5) { v019 = 1e0 / v[672]; v020 = (-(v019 / v[672])); v021 = (2e0*v019) / Power(v[672], 2); }
//				else {
//					v019 = (24000000e0 - (-1e0 + 1000000e0*v[672])*(71999994e0 - 0.71999982e14*v[672] + 0.6e19*Power(v[672], 3)
//						+ 0.23999982e20*(v[672] * v[672]))) / 24e0;
//					v020 = (-500000e0 / 3e0)*(35999994e0 - 0.47999982e14*v[672] + 0.6e19*Power(v[672], 3) + 0.17999982e20*
//						(v[672] * v[672]));
//					v021 = 0.1e13*(7999997e0 - 0.5999994e13*v[672] - 0.3e13*(v[672] * v[672]));
//				};
//				v[676] = v020;
//				v[677] = v019;
//				v[2659] = (*mud)*v[677] * v[681];
//				v[638] = v[2659] * v[631];
//				v[639] = v[2659] * v[632];
//				v[640] = v[2659] * v[633];
//				v[641] = 0e0;
//			};
//			if (sqrt((v[625] * v[625]) + (v[626] * v[626]) + (v[627] * v[627])) > (*mud)*sqrt((v[622] * v[622]) +
//				(v[623] * v[623]) + (v[624] * v[624]))) {
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
//				v[690] = sqrt((v[625] * v[625]) + (v[626] * v[626]) + (v[627] * v[627]));
//				if (v[690] > 0.1e-5) { v025 = 1e0 / v[690]; v026 = (-(v025 / v[690])); v027 = (2e0*v025) / Power(v[690], 2); }
//				else {
//					v025 = (24000000e0 - (-1e0 + 1000000e0*v[690])*(71999994e0 - 0.71999982e14*v[690] + 0.6e19*Power(v[690], 3)
//						+ 0.23999982e20*(v[690] * v[690]))) / 24e0;
//					v026 = (-500000e0 / 3e0)*(35999994e0 - 0.47999982e14*v[690] + 0.6e19*Power(v[690], 3) + 0.17999982e20*
//						(v[690] * v[690]));
//					v027 = 0.1e13*(7999997e0 - 0.5999994e13*v[690] - 0.3e13*(v[690] * v[690]));
//				};
//				v[696] = -((*mud)*v022*v025*sqrt(v[1644]));
//				v[666] = v[571] + v[625] * v[696];
//				v[668] = v[575] + v[626] * v[696];
//				v[669] = v[576] + v[627] * v[696];
//			}
//			else {
//				v[666] = 0e0;
//				v[668] = 0e0;
//				v[669] = 0e0;
//			};
//		};
//	}
//	else {
//		v[638] = 0e0;
//		v[639] = 0e0;
//		v[640] = 0e0;
//	};
//	fn[0] = v[622];
//	fn[1] = v[623];
//	fn[2] = v[624];
//	ft[0] = v[638];
//	ft[1] = v[639];
//	ft[2] = v[640];
//	(*stickupdated) = v[641];
//	gtpupdated[0] = v[571] - v[666];
//	gtpupdated[1] = v[575] - v[668];
//	gtpupdated[2] = v[576] - v[669];
//	b710 = b598;
//	if (b710) {
//		b711 = b600;
//	}
//	else {
//	};
//	v[715] = v[605] * v[718];
//	v[716] = v[605] * v[720];
//	v[740] = v[172] * v[716];
//	v[717] = v[605] * v[722];
//	v[746] = v[172] * v[717];
//	v[719] = v[604] * v[718];
//	v[741] = v[172] * v[719];
//	v[721] = v[604] * v[720];
//	v[743] = -(v[172] * v[721]) / 2e0;
//	v[723] = v[604] * v[722];
//	v[749] = v[172] * v[723];
//	v[724] = v[602] * v[718];
//	v[747] = v[172] * v[724];
//	v[725] = v[602] * v[720];
//	v[750] = v[172] * v[725];
//	v[726] = v[602] * v[722];
//	v[748] = -(v[172] * v[726]) / 2e0;
//	v[727] = -(v[193] * v[602]) - v[196] * v[604] - v[199] * v[605];
//	v[728] = -(v[192] * v[602]) - v[195] * v[604] - v[198] * v[605];
//	v[729] = -(v[191] * v[602]) - v[194] * v[604] - v[197] * v[605];
//	v[730] = v[144] * v[602] + v[148] * v[604] + v[152] * v[605];
//	v[731] = v[740] + v[741];
//	v[862] = v[731] / 2e0;
//	v[858] = (v[746] + v[747]) / 2e0;
//	v[733] = v[749] + v[750];
//	v[857] = v[733] / 2e0;
//	v[734] = (v[351] * v[715]) / 2e0 + v[346] * v[716] + v[342] * v[717] + v[338] * v[719] + (v[333] * v[721]) / 2e0
//		+ v[329] * v[723] + v[325] * v[724] + v[319] * v[725] + (v[314] * v[726]) / 2e0;
//	v[865] = v[743] + v[748] - 4e0*v[734] * v[783];
//	v[861] = -(v[172] * v[715]) / 2e0 - v[743] + v[865];
//	v[856] = v[743] - v[748] + v[861];
//	v[4091] = v[128] * v[602];
//	v[4092] = v[128] * v[604];
//	v[4093] = v[128] * v[605];
//	v[4094] = v[130] * v[602];
//	v[4095] = v[130] * v[604];
//	v[4096] = v[130] * v[605];
//	v[4097] = v[132] * v[602];
//	v[4098] = v[132] * v[604];
//	v[4099] = v[132] * v[605];
//	v[4100] = -v[602];
//	v[4101] = -v[604];
//	v[4102] = -v[605];
//	v[4103] = v[306] * v[733] + v[740] - v[741] + 2e0*alphaB[0] * v[856] + alphaB[2] * v[858];
//	v[4104] = (alphaB[2] * v[731]) / 2e0 + (alphaB[0] * v[733]) / 2e0 - v[746] + v[747] + 2e0*alphaB[1] * v[861];
//	v[4105] = v[306] * v[731] + v[749] - v[750] + alphaB[0] * v[858] + 2e0*alphaB[2] * v[865];
//	v[735] = v[729] * x1B[0] + v[728] * x1B[1] + v[727] * x1B[2];
//	v[736] = (-v[735] + v[729] * x3B[0] + v[728] * x3B[1] + v[727] * x3B[2]) / 2e0;
//	v[737] = (-v[735] + v[729] * x2B[0] + v[728] * x2B[1] + v[727] * x2B[2]) / 2e0;
//	v[738] = (v[145] * v[602] + v[149] * v[604] + v[153] * v[605] - v[730]) / 2e0;
//	v[739] = (v[146] * v[602] + v[150] * v[604] + v[154] * v[605] - v[730]) / 2e0;
//	for (i705 = 1; i705 <= 15; i705++) {
//		i2662 = (i705 == 13 ? 1 : 0);
//		i2661 = (i705 == 14 ? 1 : 0);
//		i2660 = (i705 == 15 ? 1 : 0);
//		v[759] = v[4150 + i705];
//		v[758] = v[4135 + i705];
//		v[757] = v[4120 + i705];
//		v[756] = v[4105 + i705];
//		v[760] = (-v[758] - v[759]) / 2e0;
//		v[761] = v[4184 + i705];
//		v[785] = -4e0*v[761] * v[783];
//		v[762] = v[4199 + i705];
//		v[763] = v[4214 + i705];
//		v[764] = v[4229 + i705];
//		v[765] = v[4244 + i705];
//		v[766] = v[4304 + i705];
//		v[767] = v[4364 + i705];
//		v[768] = (-v[756] - v[757]) / 2e0;
//		v[769] = (2e0*v[760] * x1B[0] + v[758] * x2B[0] + v[759] * x3B[0]) / 2e0;
//		v[770] = (2e0*v[760] * x1B[1] + v[758] * x2B[1] + v[759] * x3B[1]) / 2e0;
//		v[771] = (2e0*v[760] * x1B[2] + v[758] * x2B[2] + v[759] * x3B[2]) / 2e0;
//		v[772] = -i2660 + v[762];
//		v[774] = i2660 + v[762];
//		v[775] = i2661 + v[763];
//		v[777] = -i2661 + v[763];
//		v[778] = -i2662 + v[764];
//		v[780] = i2662 + v[764];
//		v[782] = -(v[172] * v[765]) / 2e0 + (v[314] * v[785]) / 2e0;
//		v[784] = v[172] * v[772] + v[319] * v[785];
//		v[786] = v[172] * v[775] + v[325] * v[785];
//		v[787] = v[4379 + i705] + (v[146] * v[756]) / 2e0 + (v[145] * v[757]) / 2e0 + v[144] * v[768] - v[191] * v[769]
//			- v[192] * v[770] - v[193] * v[771] + v[722] * v[782] + v[720] * v[784] + v[718] * v[786];
//		v[807] = v[787];
//		v[788] = v[172] * v[774] + v[329] * v[785];
//		v[789] = (-(v[172] * v[766]) + v[333] * v[785]) / 2e0;
//		v[790] = v[172] * v[778] + v[338] * v[785];
//		v[791] = v[4394 + i705] + (v[150] * v[756]) / 2e0 + (v[149] * v[757]) / 2e0 + v[148] * v[768] - v[194] * v[769]
//			- v[195] * v[770] - v[196] * v[771] + v[722] * v[788] + v[720] * v[789] + v[718] * v[790];
//		v[806] = v[791];
//		v[792] = v[172] * v[777] + v[342] * v[785];
//		v[793] = v[602] * v[782] + v[604] * v[788] + v[605] * v[792];
//		v[794] = v[172] * v[780] + v[346] * v[785];
//		v[795] = v[602] * v[784] + v[604] * v[789] + v[605] * v[794];
//		v[796] = (-(v[172] * v[767]) + v[351] * v[785]) / 2e0;
//		v[797] = v[4409 + i705] + (v[154] * v[756]) / 2e0 + (v[153] * v[757]) / 2e0 + v[152] * v[768] - v[197] * v[769]
//			- v[198] * v[770] - v[199] * v[771] + v[722] * v[792] + v[720] * v[794] + v[718] * v[796];
//		v[805] = v[797];
//		v[798] = v[602] * v[786] + v[604] * v[790] + v[605] * v[796];
//		v[799] = 0e0;
//		v[800] = 0e0;
//		v[801] = 0e0;
//		v[802] = 0e0;
//		b803 = b598;
//		if (b803) {
//			v[2664] = -(v[267] * v[805]) - v[266] * v[806] - v[265] * v[807];
//			b804 = b600;
//			if (b804) {
//				v[801] = v[2663] * v[797];
//				v[797] = 0e0;
//				v[800] = v[2663] * v[791];
//				v[791] = 0e0;
//				v[799] = v[2663] * v[787];
//				v[787] = 0e0;
//				v[802] = -(v[2658] * v[2664] * v[2970]);
//			}
//			else {
//				v[801] = v[606] * v[805];
//				v[797] = 0e0;
//				v[800] = v[606] * v[806];
//				v[791] = 0e0;
//				v[799] = v[606] * v[807];
//				v[787] = 0e0;
//				v[802] = (*n2)*v[2664] * v[2974] * v[594];
//			};
//		}
//		else {
//		};
//		v[802] = v[263] * (v[250] * v[799] + v[251] * v[800] + v[252] * v[801]) + v[802];
//		v[2665] = v[802] / v[597];
//		v[812] = v[252] * v[2665] + v[264] * v[801];
//		v[814] = v[251] * v[2665] + v[264] * v[800];
//		v[815] = v[250] * v[2665] + v[264] * v[799];
//		v[816] = -(v[605] * v[771]) - v[168] * v[812];
//		v[817] = -(v[605] * v[770]) - v[167] * v[812];
//		v[818] = -(v[605] * v[769]) - v[166] * v[812];
//		v[819] = -(v[604] * v[771]) - v[168] * v[814];
//		v[820] = -(v[604] * v[770]) - v[167] * v[814];
//		v[821] = -(v[604] * v[769]) - v[166] * v[814];
//		v[822] = -(v[602] * v[771]) - v[168] * v[815];
//		v[823] = -(v[602] * v[770]) - v[167] * v[815];
//		v[824] = -(v[602] * v[769]) - v[166] * v[815];
//		v[825] = v[4454 + i705] + v[152] * v[812] + v[148] * v[814] + v[144] * v[815];
//		v[826] = QBi[2][2] * v[816] + QBi[2][1] * v[817] + QBi[2][0] * v[818];
//		v[827] = QBi[1][2] * v[816] + QBi[1][1] * v[817] + QBi[1][0] * v[818];
//		v[828] = QBi[0][2] * v[816] + QBi[0][1] * v[817] + QBi[0][0] * v[818];
//		v[829] = QBi[2][2] * v[819] + QBi[2][1] * v[820] + QBi[2][0] * v[821];
//		v[830] = QBi[1][2] * v[819] + QBi[1][1] * v[820] + QBi[1][0] * v[821];
//		v[831] = QBi[0][2] * v[819] + QBi[0][1] * v[820] + QBi[0][0] * v[821];
//		v[832] = QBi[2][2] * v[822] + QBi[2][1] * v[823] + QBi[2][0] * v[824];
//		v[833] = QBi[1][2] * v[822] + QBi[1][1] * v[823] + QBi[1][0] * v[824];
//		v[834] = QBi[0][2] * v[822] + QBi[0][1] * v[823] + QBi[0][0] * v[824];
//		v[835] = (v[715] * v[785] + v[172] * v[826]) / 2e0;
//		v[836] = v[716] * v[785] + v[172] * v[827];
//		v[837] = v[717] * v[785] + v[172] * v[828];
//		v[838] = v[719] * v[785] + v[172] * v[829];
//		v[840] = v[723] * v[785] + v[172] * v[831];
//		v[841] = v[724] * v[785] + v[172] * v[832];
//		v[842] = v[725] * v[785] + v[172] * v[833];
//		v[2670] = 8e0*v[1215] * v[734] * v[761] - 4e0*v[783] * (-(v[726] * v[765]) / 2e0 - (v[721] * v[766]) / 2e0 -
//			(v[715] * v[767]) / 2e0 + v[725] * v[772] + v[723] * v[774] + v[724] * v[775] + v[717] * v[777] + v[719] * v[778]
//			+ v[716] * v[780] + (v[351] * v[826]) / 2e0 + v[346] * v[827] + v[342] * v[828] + v[338] * v[829] + (v[333] * v[830]) / 2e0
//			+ v[329] * v[831] + v[325] * v[832] + v[319] * v[833] + (v[314] * v[834]) / 2e0);
//		v[860] = v[2670] - (v[721] * v[785]) / 2e0 - (v[172] * v[830]) / 2e0;
//		v[844] = (v[726] * v[785] + v[172] * v[834]) / 2e0;
//		v[845] = -(QBi[0][2] * v[793]) - QBi[1][2] * v[795] - QBi[2][2] * v[798] - v[199] * v[812] - v[196] * v[814]
//			- v[193] * v[815];
//		v[846] = -(QBi[0][1] * v[793]) - QBi[1][1] * v[795] - QBi[2][1] * v[798] - v[198] * v[812] - v[195] * v[814]
//			- v[192] * v[815];
//		v[847] = -(QBi[0][0] * v[793]) - QBi[1][0] * v[795] - QBi[2][0] * v[798] - v[197] * v[812] - v[194] * v[814]
//			- v[191] * v[815];
//		v[848] = v[847] * x1B[0] + v[846] * x1B[1] + v[845] * x1B[2];
//		v[849] = (-v[848] + v[847] * x3B[0] + v[846] * x3B[1] + v[845] * x3B[2]) / 2e0;
//		v[850] = (-v[848] + v[847] * x2B[0] + v[846] * x2B[1] + v[845] * x2B[2]) / 2e0;
//		v[2669] = (v[837] + v[841]) / 2e0;
//		v[852] = v[836] + v[838];
//		v[853] = v[840] + v[842];
//		v[4470] = 0e0;
//		v[4471] = 0e0;
//		v[4472] = 0e0;
//		v[4473] = 0e0;
//		v[4474] = 0e0;
//		v[4475] = 0e0;
//		v[4476] = 0e0;
//		v[4477] = 0e0;
//		v[4478] = 0e0;
//		v[4479] = 0e0;
//		v[4480] = 0e0;
//		v[4481] = 0e0;
//		v[4482] = 2e0*v[856];
//		v[4483] = v[857];
//		v[4484] = v[858];
//		v[4485] = 0e0;
//		v[4486] = 0e0;
//		v[4487] = 0e0;
//		v[4488] = 0e0;
//		v[4489] = 0e0;
//		v[4490] = 0e0;
//		v[4491] = 0e0;
//		v[4492] = 0e0;
//		v[4493] = 0e0;
//		v[4494] = 0e0;
//		v[4495] = 0e0;
//		v[4496] = 0e0;
//		v[4497] = v[857];
//		v[4498] = 2e0*v[861];
//		v[4499] = v[862];
//		v[4500] = 0e0;
//		v[4501] = 0e0;
//		v[4502] = 0e0;
//		v[4503] = 0e0;
//		v[4504] = 0e0;
//		v[4505] = 0e0;
//		v[4506] = 0e0;
//		v[4507] = 0e0;
//		v[4508] = 0e0;
//		v[4509] = 0e0;
//		v[4510] = 0e0;
//		v[4511] = 0e0;
//		v[4512] = v[858];
//		v[4513] = v[862];
//		v[4514] = 2e0*v[865];
//		v[4515] = v[602] * v[768] + v[128] * v[815];
//		v[4516] = v[604] * v[768] + v[128] * v[814];
//		v[4517] = v[605] * v[768] + v[128] * v[812];
//		v[4518] = v[2666] * v[757] + v[130] * v[815];
//		v[4519] = v[2667] * v[757] + v[130] * v[814];
//		v[4520] = v[2668] * v[757] + v[130] * v[812];
//		v[4521] = v[2666] * v[756] + v[132] * v[815];
//		v[4522] = v[2667] * v[756] + v[132] * v[814];
//		v[4523] = v[2668] * v[756] + v[132] * v[812];
//		v[4524] = -v[815];
//		v[4525] = -v[814];
//		v[4526] = -v[812];
//		v[4527] = alphaB[2] * v[2669] + v[4469 + i705] + v[836] - v[838] + v[306] * v[853] + 2e0*alphaB[0] * (-v[835] + v[860]);
//		v[4528] = v[4484 + i705] - v[837] + v[841] + 2e0*alphaB[1] * (v[2670] - v[835] - v[844]) + (alphaB[2] * v[852]) / 2e0 +
//			(alphaB[0] * v[853]) / 2e0;
//		v[4529] = alphaB[0] * v[2669] + v[4499 + i705] + v[840] - v[842] + v[306] * v[852] + 2e0*alphaB[2] * (-v[844] + v[860]);
//		v[854] = (v[4439 + i705] + v[153] * v[812] + v[149] * v[814] + v[145] * v[815] - v[825]) / 2e0;
//		v[855] = (v[4424 + i705] + v[154] * v[812] + v[150] * v[814] + v[146] * v[815] - v[825]) / 2e0;
//		Rc[i705 - 1] += v[4090 + i705] + v[739] * v[756] + v[738] * v[757] + v[737] * v[758] + v[736] * v[759];
//		for (i753 = 1; i753 <= 15; i753++) {
//			Kc[i705 - 1][i753 - 1] += v[4514 + i753] + v[4150 + i753] * v[849] + v[4135 + i753] * v[850] + v[4120 + i753] * v[854]
//				+ v[4105 + i753] * v[855];
//		};/* end for */
//	};/* end for */
//	v[935] = 0e0;
//	v[936] = 0e0;
//	v[937] = 0e0;
//	b938 = b598;
//	if (b938) {
//		b939 = (*stick);
//		if (b939) {
//			b940 = b636;
//			if (b940) {
//				v[937] = 0e0;
//				v[936] = 0e0;
//				v[935] = 0e0;
//			}
//			else {
//			};
//		}
//		else {
//			b941 = b670;
//			if (b941) {
//				v[937] = 0e0;
//				v[936] = 0e0;
//				v[935] = 0e0;
//			}
//			else {
//			};
//		};
//	}
//	else {
//	};
//	v[2674] = (*ct)*v[935];
//	v[2766] = -v[638] + v[2671] * v[935];
//	v[2673] = (*ct)*v[936];
//	v[2764] = -v[639] + v[2671] * v[936];
//	v[2672] = (*ct)*v[937];
//	v[2762] = -v[640] + v[2671] * v[937];
//	v[945] = v[2672] * v[949];
//	v[946] = v[2672] * v[951];
//	v[2704] = (v[2672] * v[953]) / 2e0;
//	v[2705] = (v[2672] * v[955]) / 2e0;
//	v[950] = v[2673] * v[949];
//	v[952] = v[2673] * v[951];
//	v[2706] = (v[2673] * v[953]) / 2e0;
//	v[2707] = (v[2673] * v[955]) / 2e0;
//	v[957] = v[2674] * v[949];
//	v[958] = v[2674] * v[951];
//	v[2708] = (v[2674] * v[953]) / 2e0;
//	v[2709] = (v[2674] * v[955]) / 2e0;
//	v[961] = (*epst)*v[937];
//	v[1515] = -(v[267] * v[961]);
//	v[962] = (*epst)*v[936];
//	v[1518] = -(v[266] * v[962]);
//	v[1513] = v[1515] + v[1518];
//	v[963] = (*epst)*v[935];
//	v[1519] = -(v[265] * v[963]);
//	v[1520] = v[1518] + v[1519];
//	v[1516] = v[1515] + v[1519];
//	v[964] = 0e0;
//	v[965] = 0e0;
//	v[966] = 0e0;
//	v[967] = 0e0;
//	v[968] = 0e0;
//	b969 = b598;
//	if (b969) {
//		v[970] = 0e0;
//		v[971] = 0e0;
//		v[972] = 0e0;
//		b973 = b617;
//		if (b973) {
//			v[972] = 0e0;
//			v[971] = 0e0;
//			v[970] = 0e0;
//		}
//		else {
//		};
//		v[982] = v[578] * v[970] + v[586] * v[971] + v[588] * v[972];
//		v[2675] = v[982] * (*zetan);
//		b974 = b600;
//		if (b974) {
//			v[966] = v[612] * v[972];
//			v[965] = v[612] * v[971];
//			v[964] = v[612] * v[970];
//			v[968] = (v[2675] * v[610] * v[975] * Power(v[609], v[1602])) / v[976];
//		}
//		else {
//			v[966] = v[616] * v[972];
//			v[965] = v[616] * v[971];
//			v[964] = v[616] * v[970];
//			v[967] = (v[2675] * v[615] * v[980] * Power(v[597], v[1607])) / v[981];
//		};
//	}
//	else {
//	};
//	v[5693] = v[128] * v[964];
//	v[5694] = 0e0;
//	v[5695] = 0e0;
//	v[5696] = v[130] * v[964];
//	v[5697] = 0e0;
//	v[5698] = 0e0;
//	v[5699] = v[132] * v[964];
//	v[5700] = 0e0;
//	v[5701] = 0e0;
//	v[5702] = -v[964];
//	v[5703] = 0e0;
//	v[5704] = 0e0;
//	v[5705] = 0e0;
//	v[5706] = 0e0;
//	v[5707] = 0e0;
//	v[1991] = v[265] * v[964];
//	v[1849] = -(v[577] * v[964]);
//	v[1831] = -(v[1991] * v[249]);
//	v[1829] = -(v[1991] * v[239]);
//	v[1827] = -(v[1991] * v[229]);
//	v[5663] = 0e0;
//	v[5664] = 0e0;
//	v[5665] = 0e0;
//	v[5666] = 0e0;
//	v[5667] = 0e0;
//	v[5668] = 0e0;
//	v[5669] = v[965];
//	v[5670] = v[964];
//	v[5671] = 0e0;
//	v[5672] = 0e0;
//	v[5673] = 0e0;
//	v[5674] = 0e0;
//	v[5675] = 0e0;
//	v[5676] = 0e0;
//	v[5677] = 0e0;
//	v[5648] = 0e0;
//	v[5649] = 0e0;
//	v[5650] = 0e0;
//	v[5651] = v[965];
//	v[5652] = v[964];
//	v[5653] = 0e0;
//	v[5654] = 0e0;
//	v[5655] = 0e0;
//	v[5656] = 0e0;
//	v[5657] = 0e0;
//	v[5658] = 0e0;
//	v[5659] = 0e0;
//	v[5660] = 0e0;
//	v[5661] = 0e0;
//	v[5662] = 0e0;
//	v[5633] = v[965];
//	v[5634] = v[964];
//	v[5635] = 0e0;
//	v[5636] = 0e0;
//	v[5637] = 0e0;
//	v[5638] = 0e0;
//	v[5639] = 0e0;
//	v[5640] = 0e0;
//	v[5641] = 0e0;
//	v[5642] = 0e0;
//	v[5643] = 0e0;
//	v[5644] = 0e0;
//	v[5645] = 0e0;
//	v[5646] = 0e0;
//	v[5647] = 0e0;
//	v[5678] = 0e0;
//	v[5679] = v[128] * v[965];
//	v[5680] = 0e0;
//	v[5681] = 0e0;
//	v[5682] = v[130] * v[965];
//	v[5683] = 0e0;
//	v[5684] = 0e0;
//	v[5685] = v[132] * v[965];
//	v[5686] = 0e0;
//	v[5687] = 0e0;
//	v[5688] = -v[965];
//	v[5689] = 0e0;
//	v[5690] = 0e0;
//	v[5691] = 0e0;
//	v[5692] = 0e0;
//	v[1989] = v[266] * v[965];
//	v[1850] = -(v[582] * v[965]);
//	v[1839] = -(v[1989] * v[249]);
//	v[1836] = -(v[1989] * v[239]);
//	v[1833] = -(v[1989] * v[229]);
//	v[5603] = 0e0;
//	v[5604] = 0e0;
//	v[5605] = v[128] * v[966];
//	v[5606] = 0e0;
//	v[5607] = 0e0;
//	v[5608] = v[130] * v[966];
//	v[5609] = 0e0;
//	v[5610] = 0e0;
//	v[5611] = v[132] * v[966];
//	v[5612] = 0e0;
//	v[5613] = 0e0;
//	v[5614] = -v[966];
//	v[5615] = 0e0;
//	v[5616] = 0e0;
//	v[5617] = 0e0;
//	v[1990] = -(v[267] * v[966]);
//	v[2789] = v[1990] - v[1991];
//	v[3139] = v[2789] * v[404] + v[2648] * v[965];
//	v[3136] = v[2789] * v[405] + v[2646] * v[965];
//	v[3133] = v[2789] * v[406] + v[2644] * v[965];
//	v[2788] = -v[1989] + v[1990];
//	v[3138] = v[2788] * v[401] - v[2649] * v[964];
//	v[3135] = v[2788] * v[402] - v[2647] * v[964];
//	v[3132] = v[2788] * v[403] - v[2645] * v[964];
//	v[1840] = v[1990] * v[249];
//	v[1837] = v[1990] * v[239];
//	v[1834] = v[1990] * v[229];
//	v[1808] = -(v[1989] * v[267]) - v[1991] * v[267] - v[587] * v[966];
//	v[1032] = v[1990] * v[266];
//	v[1977] = -v[1032] - v[1850];
//	v[1809] = v[1032] + v[1850] - v[1991] * v[266];
//	v[1027] = v[1990] * v[265];
//	v[1976] = -v[1027] - v[1849];
//	v[5813] = 0e0;
//	v[5814] = 0e0;
//	v[5815] = 0e0;
//	v[5816] = 0e0;
//	v[5817] = 0e0;
//	v[5818] = 0e0;
//	v[5819] = v[1976];
//	v[5820] = v[1977];
//	v[5821] = -v[1808];
//	v[5822] = 0e0;
//	v[5823] = 0e0;
//	v[5824] = 0e0;
//	v[5825] = 0e0;
//	v[5826] = 0e0;
//	v[5827] = 0e0;
//	v[5828] = 0e0;
//	v[5829] = 0e0;
//	v[5830] = 0e0;
//	v[5831] = v[1976];
//	v[5832] = v[1977];
//	v[5833] = -v[1808];
//	v[5834] = 0e0;
//	v[5835] = 0e0;
//	v[5836] = 0e0;
//	v[5837] = 0e0;
//	v[5838] = 0e0;
//	v[5839] = 0e0;
//	v[5840] = 0e0;
//	v[5841] = 0e0;
//	v[5842] = 0e0;
//	v[5843] = v[1976];
//	v[5844] = v[1977];
//	v[5845] = -v[1808];
//	v[5846] = 0e0;
//	v[5847] = 0e0;
//	v[5848] = 0e0;
//	v[5849] = 0e0;
//	v[5850] = 0e0;
//	v[5851] = 0e0;
//	v[5852] = 0e0;
//	v[5853] = 0e0;
//	v[5854] = 0e0;
//	v[5855] = 0e0;
//	v[5856] = 0e0;
//	v[5857] = 0e0;
//	v[1810] = v[1027] + v[1849] - v[1989] * v[265];
//	v[983] = 0e0;
//	v[984] = 0e0;
//	v[985] = 0e0;
//	b986 = b598;
//	if (b986) {
//		b987 = b600;
//		if (b987) {
//			v[985] = 0e0;
//			v[984] = 0e0;
//			v[983] = 0e0;
//			v[967] = v[967] - v[968];
//		}
//		else {
//		};
//	}
//	else {
//	};
//	v[988] = v[1574] * v[966];
//	v[989] = v[1579] * v[965];
//	v[990] = v[1989] + v[1991];
//	v[3137] = v[1821] * v[966] - v[407] * v[990];
//	v[3134] = v[1817] * v[966] - v[408] * v[990];
//	v[3131] = v[1813] * v[966] - v[409] * v[990];
//	v[5618] = 0e0;
//	v[5619] = 0e0;
//	v[5620] = v[128] * v[990];
//	v[5621] = 0e0;
//	v[5622] = 0e0;
//	v[5623] = v[130] * v[990];
//	v[5624] = 0e0;
//	v[5625] = 0e0;
//	v[5626] = v[132] * v[990];
//	v[5627] = 0e0;
//	v[5628] = 0e0;
//	v[5629] = 0e0;
//	v[5630] = 0e0;
//	v[5631] = 0e0;
//	v[5632] = 0e0;
//	v[985] = v[1588] * v[964] + v[1582] * v[965] + v[1577] * v[966] + v[985];
//	v[997] = v[213] * v[964] + v[212] * v[965];
//	v[998] = v[216] * v[964] + v[215] * v[965];
//	v[999] = v[219] * v[964] + v[218] * v[965];
//	v[1407] = v[128] * v[997] + v[130] * v[998] + v[132] * v[999];
//	v[984] = v[1587] * v[964] + v[1581] * v[965] + v[1576] * v[966] + v[984];
//	v[1004] = v[1585] * v[964];
//	v[983] = v[1586] * v[964] + v[1580] * v[965] + v[1575] * v[966] + v[983];
//	v[1016] = -((*ct)*(v[928] * v[935] + v[929] * v[936] + v[930] * v[937])) + v[1013] * v[964] + v[1014] * v[965]
//		+ v[1015] * v[966];
//	v[1158] = v[1016] * v[1222];
//	v[1144] = v[1016] * v[1994];
//	v[1020] = -((*ct)*(v[924] * v[935] + v[925] * v[936] + v[926] * v[937])) + v[1017] * v[964] + v[1018] * v[965]
//		+ v[1019] * v[966];
//	v[1151] = v[1020] * v[1993];
//	v[2683] = v[1151] + v[1016] * v[1221];
//	v[2682] = v[1144] + v[1020] * v[1219];
//	v[1024] = -((*ct)*(v[920] * v[935] + v[921] * v[936] + v[922] * v[937])) + v[1021] * v[964] + v[1022] * v[965]
//		+ v[1023] * v[966];
//	v[1156] = v[1024] * v[1992];
//	v[4549] = v[128] * v[1976] + (*ct)*(-(v[872] * v[935]) - v[873] * v[936] - v[874] * v[937]) + v[579] * v[965];
//	v[4550] = v[128] * v[1977] + (*ct)*(-(v[876] * v[935]) - v[877] * v[936] - v[878] * v[937]) + v[579] * v[964];
//	v[4551] = -(v[128] * v[1808]) + (*ct)*(-(v[880] * v[935]) - v[881] * v[936] - v[882] * v[937]);
//	v[4552] = v[130] * v[1976] + (*ct)*(-(v[884] * v[935]) - v[885] * v[936] - v[886] * v[937]) + v[580] * v[965];
//	v[4553] = v[130] * v[1977] + (*ct)*(-(v[888] * v[935]) - v[889] * v[936] - v[890] * v[937]) + v[580] * v[964];
//	v[4554] = -(v[130] * v[1808]) + (*ct)*(-(v[892] * v[935]) - v[893] * v[936] - v[894] * v[937]);
//	v[4555] = v[132] * v[1976] + (*ct)*(-(v[896] * v[935]) - v[897] * v[936] - v[898] * v[937]) + v[581] * v[965];
//	v[4556] = v[132] * v[1977] + (*ct)*(-(v[900] * v[935]) - v[901] * v[936] - v[902] * v[937]) + v[581] * v[964];
//	v[4557] = -(v[132] * v[1808]) + (*ct)*(-(v[904] * v[935]) - v[905] * v[936] - v[906] * v[937]);
//	v[4558] = v[1810] + (*ct)*(-(v[908] * v[935]) - v[909] * v[936] - v[910] * v[937]);
//	v[4559] = v[1809] + (*ct)*(-(v[912] * v[935]) - v[913] * v[936] - v[914] * v[937]);
//	v[4560] = v[1808] + (*ct)*(-(v[916] * v[935]) - v[917] * v[936] - v[918] * v[937]);
//	v[4561] = v[1156] + v[1020] * v[2630] + v[1016] * v[2633];
//	v[4562] = v[1151] + v[1016] * v[1219] + v[1024] * v[1220];
//	v[4563] = v[1144] + v[1020] * v[1221] + v[1024] * v[1222];
//	v[2848] = v[1156] + v[1158];
//	v[2684] = v[1156] + v[1020] * v[1220];
//	v[2838] = v[1158] + v[2684];
//	v[1153] = v[1024] * v[2630];
//	v[2844] = v[1153] + v[2683];
//	v[2840] = v[1151] + v[1153];
//	v[1147] = v[1024] * v[2633];
//	v[2847] = v[1147] + v[2682];
//	v[2839] = v[1144] + v[1147];
//	v[1025] = v[1810] * v[239];
//	v[1028] = v[1810] * v[249];
//	v[1029] = v[1810] * v[229];
//	v[1030] = v[1809] * v[229];
//	v[1033] = v[1809] * v[249];
//	v[1034] = v[1808] * v[249];
//	v[1035] = v[1809] * v[239];
//	v[1038] = v[1808] * v[229];
//	v[1039] = v[1808] * v[239];
//	v[1040] = 0e0;
//	v[1041] = 0e0;
//	v[1042] = 0e0;
//	v[1043] = 0e0;
//	v[1044] = 0e0;
//	v[1045] = 0e0;
//	v[1046] = 0e0;
//	v[1047] = 0e0;
//	v[1048] = 0e0;
//	b1049 = (*previouscontact);
//	if (b1049) {
//		v[988] = -(v[574] * v[961]) + v[988];
//		v[989] = -(v[573] * v[962]) + v[989];
//		v[1051] = v[1513] * v[265] + v[1050] * v[963];
//		v[1053] = v[1516] * v[266] + v[1052] * v[962];
//		v[1055] = v[1520] * v[267] + v[1054] * v[961];
//		v[985] = v[1520] * v[574] + v[2676] * v[961] + v[985];
//		v[984] = v[1516] * v[573] + v[2677] * v[962] + v[984];
//		v[1004] = v[1004] - v[572] * v[963];
//		v[983] = v[1513] * v[572] - v[2678] * v[963] + v[983];
//		v[1040] = gti[0] * v[1051];
//		v[1041] = gti[1] * v[1051];
//		v[1042] = gti[2] * v[1051];
//		v[1059] = -v[1051];
//		v[1060] = -(v[1051] * v[171]);
//		v[1061] = -(v[1051] * v[170]);
//		v[1062] = -(v[1051] * v[169]);
//		v[1063] = v[1051] * v[137];
//		v[1064] = v[1051] * v[135];
//		v[1065] = v[1051] * v[133];
//		v[1043] = gti[0] * v[1053];
//		v[1044] = gti[1] * v[1053];
//		v[1045] = gti[2] * v[1053];
//		v[1066] = -v[1053];
//		v[1067] = -(v[1053] * v[171]);
//		v[1068] = -(v[1053] * v[170]);
//		v[1069] = -(v[1053] * v[169]);
//		v[1070] = v[1053] * v[137];
//		v[1071] = v[1053] * v[135];
//		v[1072] = v[1053] * v[133];
//		v[1046] = gti[0] * v[1055];
//		v[1047] = gti[1] * v[1055];
//		v[1048] = gti[2] * v[1055];
//		v[1073] = -v[1055];
//		v[1074] = -(v[1055] * v[171]);
//		v[1075] = -(v[1055] * v[170]);
//		v[1076] = -(v[1055] * v[169]);
//		v[1077] = v[1055] * v[137];
//		v[1078] = v[1055] * v[135];
//		v[1079] = v[1055] * v[133];
//	}
//	else {
//		v[1065] = 0e0;
//		v[1064] = 0e0;
//		v[1063] = 0e0;
//		v[1072] = 0e0;
//		v[1071] = 0e0;
//		v[1070] = 0e0;
//		v[1079] = 0e0;
//		v[1078] = 0e0;
//		v[1077] = 0e0;
//		v[1062] = 0e0;
//		v[1061] = 0e0;
//		v[1060] = 0e0;
//		v[1069] = 0e0;
//		v[1068] = 0e0;
//		v[1067] = 0e0;
//		v[1076] = 0e0;
//		v[1075] = 0e0;
//		v[1074] = 0e0;
//		v[1059] = 0e0;
//		v[1066] = 0e0;
//		v[1073] = 0e0;
//	};
//	b1080 = b537;
//	if (b1080) {
//		v[1106] = -(v[1048] * v[554]) / 2e0;
//		v[1105] = -(v[1044] * v[554]) / 2e0;
//		v[1104] = v[1047] * v[554];
//		v[1103] = v[1045] * v[554];
//		v[1099] = v[1046] * v[554];
//		v[1098] = v[1042] * v[554];
//		v[1095] = v[1043] * v[554];
//		v[1094] = v[1041] * v[554];
//		v[1081] = v[1103] + v[1104];
//		v[1082] = v[1098] + v[1099];
//		v[1083] = v[1094] + v[1095];
//		v[1093] = (v[1040] * v[1084]) / 2e0 + v[1041] * v[1085] + v[1042] * v[1086] + v[1043] * v[1087] + (v[1044] * v[1088])
//			/ 2e0 + v[1045] * v[1089] + v[1046] * v[1090] + v[1047] * v[1091] + (v[1048] * v[1092]) / 2e0;
//		v[1479] = (-4e0*v[1093]) / Power(v[1096], 2) + v[1105] + v[1106];
//		v[1478] = -v[1105] + v[1479] - (v[1040] * v[554]) / 2e0;
//		v[1477] = v[1105] - v[1106] + v[1478];
//		v[1097] = (-2e0*v[1094] + 2e0*v[1095] + v[1082] * v[551] + v[1081] * v[552] + 4e0*v[1477] * v[553]) / 2e0;
//		v[1102] = (2e0*v[1098] - 2e0*v[1099] + v[1083] * v[551] + 4e0*v[1478] * v[552] + v[1081] * v[553]) / 2e0;
//		v[1107] = (-2e0*v[1103] + 2e0*v[1104] + 4e0*v[1479] * v[551] + v[1083] * v[552] + v[1082] * v[553]) / 2e0;
//		v[2679] = v[1107] * v[539] + v[1102] * v[540] + v[1097] * v[541];
//		v[1465] = v[2679] * v[550];
//		v[1462] = v[2679] * v[544];
//		v[1110] = v[1462] * v[549] + v[1465] / (Power(cos(v[1108]), 2)*sqrt(v[1466]));
//		v[2728] = v[1110] / v[542];
//		v[2680] = v[1110] / v[542];
//		v[1111] = v[1097] * v[2653] + v[2680] * v[541];
//		v[1113] = v[1102] * v[2653] + v[2680] * v[540];
//		v[1114] = v[1107] * v[2653] + v[2680] * v[539];
//		v[983] = -(v[1111] * v[275]) + v[1113] * v[276] + v[983];
//		v[985] = -(v[1113] * v[274]) + v[1114] * v[275] + v[985];
//		v[984] = v[1111] * v[274] - v[1114] * v[276] + v[984];
//	}
//	else {
//	};
//	v[985] = v[985] + 2e0*v[267] * v[988] + v[1411] * v[990];
//	v[984] = v[1407] * v[265] + v[984] + 2e0*v[266] * v[989];
//	v[983] = 2e0*v[1004] * v[265] + v[1407] * v[266] + v[983];
//	v[1842] = v[250] * v[983] + v[251] * v[984] + v[252] * v[985];
//	v[967] = v[1842] * v[263] + v[967];
//	v[2681] = v[967] / v[597];
//	v[1115] = v[252] * v[2681] + v[264] * v[985];
//	v[1116] = v[251] * v[2681] + v[264] * v[984];
//	v[1117] = v[250] * v[2681] + v[264] * v[983];
//	v[1073] = v[1073] - v[1115];
//	v[1066] = v[1066] - v[1116];
//	v[1059] = v[1059] - v[1117];
//	v[1119] = v[1016] * v[1541];
//	v[1120] = v[1016] * v[1540];
//	v[1121] = v[1016] * v[1539];
//	v[1124] = v[1020] * v[1536];
//	v[1125] = v[1016] * v[236] + v[1020] * v[238];
//	v[2688] = v[1125] * v[177];
//	v[1128] = v[1020] * v[1534];
//	v[1129] = v[1119] + v[1124];
//	v[1131] = v[1020] * v[1535] + v[2688] / v[227];
//	v[1134] = v[1024] * v[1529];
//	v[1135] = v[1016] * v[230] + v[1024] * v[238];
//	v[2689] = v[1135] * v[182];
//	v[1136] = v[1020] * v[230] + v[1024] * v[236];
//	v[2690] = v[1136] * v[186];
//	v[3152] = -(v[1016] * v[2579]) - v[1020] * v[2580] - v[1024] * v[2581] + v[2629] * v[2688] + v[226] * v[2689]
//		+ v[225] * v[2690];
//	v[1138] = v[1024] * v[1531] + v[2689] / v[227];
//	v[1139] = v[1131] + v[1138];
//	v[1141] = v[1024] * v[1530] + v[2690] / v[227];
//	v[1142] = v[1120] + v[1141];
//	v[1074] = v[1074] - v[1115] * v[168] + v[290] * v[945] + v[300] * v[946];
//	v[1075] = v[1075] - v[1115] * v[167] + v[289] * v[945] + v[298] * v[946];
//	v[1076] = v[1076] - v[1115] * v[166] + v[288] * v[945] + v[296] * v[946];
//	v[1143] = QBi[2][2] * v[1074] + QBi[2][1] * v[1075] + QBi[2][0] * v[1076] + v[238] * v[2847];
//	v[1146] = QBi[1][2] * v[1074] + QBi[1][1] * v[1075] + QBi[1][0] * v[1076] + v[1136] * v[2633] + v[236] * v[2682];
//	v[1148] = QBi[0][2] * v[1074] + QBi[0][1] * v[1075] + QBi[0][0] * v[1076] + v[230] * v[2839];
//	v[1067] = v[1067] - v[1116] * v[168] + v[290] * v[950] + v[300] * v[952];
//	v[1068] = v[1068] - v[1116] * v[167] + v[289] * v[950] + v[298] * v[952];
//	v[1069] = v[1069] - v[1116] * v[166] + v[288] * v[950] + v[296] * v[952];
//	v[1149] = QBi[2][2] * v[1067] + QBi[2][1] * v[1068] + QBi[2][0] * v[1069] + v[1135] * v[2630] + v[238] * v[2683];
//	v[1152] = QBi[1][2] * v[1067] + QBi[1][1] * v[1068] + QBi[1][0] * v[1069] + v[236] * v[2844];
//	v[1154] = QBi[0][2] * v[1067] + QBi[0][1] * v[1068] + QBi[0][0] * v[1069] + v[230] * v[2840];
//	v[1060] = v[1060] - v[1117] * v[168] + v[290] * v[957] + v[300] * v[958];
//	v[1061] = v[1061] - v[1117] * v[167] + v[289] * v[957] + v[298] * v[958];
//	v[1062] = v[1062] - v[1117] * v[166] + v[288] * v[957] + v[296] * v[958];
//	v[1155] = QBi[2][2] * v[1060] + QBi[2][1] * v[1061] + QBi[2][0] * v[1062] + v[1125] * v[1220] + v[238] * v[2848];
//	v[1157] = QBi[1][2] * v[1060] + QBi[1][1] * v[1061] + QBi[1][0] * v[1062] + v[236] * v[2684];
//	v[1160] = QBi[0][2] * v[1060] + QBi[0][1] * v[1061] + QBi[0][0] * v[1062] + v[230] * v[2838];
//	v[1161] = -(v[1025] * v[722]);
//	v[1162] = -(v[1025] * v[720]);
//	v[1163] = -(v[1025] * v[718]);
//	v[1164] = -(v[1028] * v[722]);
//	v[1165] = -(v[1028] * v[718]);
//	v[1166] = -(v[1028] * v[720]);
//	v[2702] = -2e0*v[1166];
//	v[1167] = -(v[1029] * v[720]);
//	v[1168] = -(v[1029] * v[718]);
//	v[1169] = -(v[1029] * v[722]);
//	v[1170] = -(v[1030] * v[722]);
//	v[1171] = -(v[1030] * v[720]);
//	v[1172] = -(v[1030] * v[718]);
//	v[2700] = -2e0*v[1172];
//	v[1173] = -(v[1033] * v[718]);
//	v[1174] = -(v[1033] * v[722]);
//	v[2703] = 2e0*v[1174];
//	v[1175] = -(v[1033] * v[720]);
//	v[1176] = -(v[1034] * v[720]);
//	v[1177] = -(v[1034] * v[722]);
//	v[1178] = -(v[1034] * v[718]);
//	v[1179] = v[1167] + v[1170] + v[1173] + v[1176] + 2e0*v[1128] * v[318] - v[1139] * v[321] - v[1129] * v[324];
//	v[1180] = -(v[1035] * v[722]);
//	v[1181] = -(v[1035] * v[718]);
//	v[1182] = -(v[1035] * v[720]);
//	v[1183] = -v[1162] - v[1165] - v[1177] - v[1180] - v[1139] * v[318] + 2e0*v[1134] * v[321] + v[1142] * v[324];
//	v[1184] = v[1157] * v[172] - v[1167] * v[311] + v[1162] * v[312] - v[1166] * v[313];
//	v[1185] = -(v[1038] * v[722]);
//	v[1186] = -(v[1038] * v[720]);
//	v[2701] = 2e0*v[1186];
//	v[1187] = -(v[1038] * v[718]);
//	v[1188] = -(v[1039] * v[720]);
//	v[1189] = -(v[1039] * v[722]);
//	v[1190] = -(v[1039] * v[718]);
//	v[1194] = -v[1168] - v[1181] - v[1185] - v[1188] - v[1129] * v[318] + v[1142] * v[321] + 2e0*v[1121] * v[324];
//	v[2825] = -v[1194] / 2e0;
//	v[1195] = v[1155] * v[172] - v[1168] * v[311] + v[1163] * v[312] - v[1165] * v[313];
//	v[1196] = v[1154] * v[172] - v[1170] * v[311] + v[1180] * v[312] - v[1174] * v[313];
//	v[1197] = v[1164] + v[1175];
//	v[1198] = v[1149] * v[172] - v[1172] * v[311] + v[1181] * v[312] - v[1173] * v[313];
//	v[1199] = v[1148] * v[172] - v[1185] * v[311] + v[1189] * v[312] - v[1177] * v[313];
//	v[1200] = v[1146] * v[172] - v[1186] * v[311] + v[1188] * v[312] - v[1176] * v[313];
//	v[1201] = v[1171] + v[1187];
//	v[1202] = v[1161] + v[1190];
//	v[5798] = 0e0;
//	v[5799] = 0e0;
//	v[5800] = 0e0;
//	v[5801] = 0e0;
//	v[5802] = 0e0;
//	v[5803] = 0e0;
//	v[5804] = 0e0;
//	v[5805] = 0e0;
//	v[5806] = 0e0;
//	v[5807] = 0e0;
//	v[5808] = 0e0;
//	v[5809] = 0e0;
//	v[5810] = -v[1183] / 2e0 - v[1201];
//	v[5811] = (v[1179] - 2e0*v[1202]) / 2e0;
//	v[5812] = -v[1197] + v[2825];
//	v[1203] = 1e0 / Power(v[227], 2);
//	v[2852] = -(v[1203] * v[237]);
//	v[2699] = -(v[1203] * v[236]);
//	v[2698] = -(v[1203] * v[238]);
//	v[2695] = -(v[1203] * v[230]);
//	v[2694] = -(v[1203] * v[225]);
//	v[2693] = -(v[1016] * v[1203]);
//	v[2692] = -(v[1203] * v[226]);
//	v[2691] = -(v[1203] * v[2629]);
//	v[1971] = -(v[1203] * (v[177] * v[228] + v[2627]));
//	v[1970] = -(v[1203] * (v[176] * v[228] + v[2628]));
//	v[1969] = -(v[1203] * v[2685]);
//	v[1968] = -(v[1203] * v[2631]);
//	v[1967] = -(v[1203] * v[2686]);
//	v[1966] = -(v[1203] * v[2687]);
//	v[1965] = -(v[1203] * v[228]);
//	v[2814] = v[1965] * v[230];
//	v[1964] = -(v[1203] * v[231]);
//	v[2813] = v[1964] * v[236];
//	v[2812] = -(v[1203] * v[238] * v[240]);
//	v[1962] = -(v[1203] * v[2688]);
//	v[1961] = -(v[1203] * v[2689]);
//	v[1960] = -(v[1203] * v[2690]);
//	v[1959] = v[1125] * v[2691];
//	v[1958] = v[1135] * v[2692];
//	v[1957] = v[1136] * v[2694];
//	v[1902] = v[1020] * v[2691];
//	v[1901] = v[2632] * v[2693];
//	v[1899] = v[1024] * v[1965];
//	v[3147] = v[1959] + (v[1899] + v[1901])*v[238];
//	v[2811] = v[1899] + v[1902];
//	v[3148] = v[1901] + v[2811];
//	v[1895] = v[1024] * v[2692];
//	v[1892] = v[1020] * v[1964];
//	v[3146] = v[1892] + v[1895];
//	v[1891] = v[2634] * v[2693];
//	v[2810] = v[1891] + v[1892];
//	v[3145] = v[1895] + v[2810];
//	v[3144] = v[1958] + v[238] * v[2810];
//	v[1887] = v[1024] * v[2694];
//	v[1885] = v[1020] * v[2852];
//	v[1884] = v[240] * v[2693];
//	v[3143] = v[1884] + v[1887];
//	v[2809] = v[1884] + v[1885];
//	v[3142] = v[1957] + v[236] * v[2809];
//	v[3141] = v[1887] + v[2809];
//	v[1880] = v[184] * v[2695];
//	v[1877] = v[179] * v[2695];
//	v[1874] = -(v[1203] * v[2696]);
//	v[1873] = -(v[1203] * v[2697]);
//	v[1868] = v[182] * v[2698];
//	v[1867] = v[181] * v[2699];
//	v[2470] = v[1867] + v[1877];
//	v[2462] = v[1868] + v[2470];
//	v[2448] = v[1867] + v[1868];
//	v[1865] = v[175] * v[2695];
//	v[2467] = v[1865] + v[1873] + v[1874];
//	v[2460] = -v[1874] + v[2467];
//	v[2450] = -v[1873] + v[2467];
//	v[1861] = v[187] * v[2698];
//	v[2473] = v[1861] + v[1880];
//	v[1860] = v[186] * v[2699];
//	v[2456] = v[1860] + v[1861];
//	v[2452] = v[1860] + v[2473];
//	v[1538] = v[1967] * v[230] + v[1966] * v[236] + v[2812];
//	v[1533] = v[1969] * v[230] + v[1968] * v[238] + v[2813];
//	v[1528] = v[1970] * v[236] + v[1971] * v[238] + v[2814];
//	v[1369] = v[177] * v[2691];
//	v[1367] = v[182] * v[2692];
//	v[1363] = v[186] * v[2694];
//	v[1956] = v[1121] + v[1128] + v[1134] + v[1136] * v[1363] + v[1135] * v[1367] + v[1125] * v[1369] + v[1024] * v[1528]
//		+ v[1020] * v[1533] + v[1016] * v[1538];
//	v[1204] = (2e0*v[1163] + alphaB[1] * v[1179] - alphaB[0] * v[1183] - 2e0*v[1189] - alphaB[2] * v[1194]
//		+ 4e0*v[172] * v[1956] + v[2700] + v[2701] + v[2702] + v[2703] - v[1201] * v[305] - v[1202] * v[308] - v[1197] * v[310]
//		+ v[1160] * v[314] + 2e0*v[1157] * v[319] + 2e0*v[1155] * v[325] + 2e0*v[1154] * v[329] + v[1152] * v[333]
//		+ 2e0*v[1149] * v[338] + 2e0*v[1148] * v[342] + 2e0*v[1146] * v[346] + v[1143] * v[351]) / 2e0;
//	v[1206] = (-2e0*v[1119] + 2e0*v[1124] - v[1169] * v[314] - 2e0*v[1167] * v[319] - 2e0*v[1168] * v[325]
//		- 2e0*v[1170] * v[329] - v[1171] * v[333] + v[2700] * v[338] - 2e0*v[1185] * v[342] - v[2701] * v[346] - v[1187] * v[351])
//		/ 2e0;
//	v[2830] = 8e0*v[1206];
//	v[1208] = -v[1120] + v[1141] + (v[1161] * v[314]) / 2e0 + v[1162] * v[319] + v[1163] * v[325] + v[1180] * v[329] +
//		(v[1182] * v[333]) / 2e0 + v[1181] * v[338] + v[1189] * v[342] + v[1188] * v[346] + (v[1190] * v[351]) / 2e0;
//	v[2829] = 8e0*v[1208];
//	v[1209] = (v[1152] * v[172] - v[1171] * v[311] + v[1182] * v[312] - v[1175] * v[313]) / 2e0;
//	v[1210] = (-2e0*v[1131] + 2e0*v[1138] - v[1164] * v[314] + v[2702] * v[319] - 2e0*v[1165] * v[325] - v[2703] * v[329]
//		- v[1175] * v[333] - 2e0*v[1173] * v[338] - 2e0*v[1177] * v[342] - 2e0*v[1176] * v[346] - v[1178] * v[351]) / 2e0;
//	v[3150] = v[1206] * v[305] - v[1208] * v[308] + v[1210] * v[310];
//	v[2828] = 8e0*v[1210];
//	v[5858] = 0e0;
//	v[5859] = 0e0;
//	v[5860] = 0e0;
//	v[5861] = 0e0;
//	v[5862] = 0e0;
//	v[5863] = 0e0;
//	v[5864] = 0e0;
//	v[5865] = 0e0;
//	v[5866] = 0e0;
//	v[5867] = 0e0;
//	v[5868] = 0e0;
//	v[5869] = 0e0;
//	v[5870] = v[2830];
//	v[5871] = -v[2829];
//	v[5872] = v[2828];
//	v[1216] = -v[1209] + v[1210] * v[1267] + v[1208] * v[1270] + v[1206] * v[1272] - v[1204] * v[2826];
//	v[1988] = v[1216] + (-(v[1160] * v[172]) + v[1169] * v[311] - v[1161] * v[312] + v[1164] * v[313]) / 2e0;
//	v[1211] = (v[1143] * v[172] - v[1187] * v[311] + v[1190] * v[312] - v[1178] * v[313]) / 2e0;
//	v[1986] = v[1209] - v[1211] + v[1988];
//	v[1982] = -v[1211] + v[1216];
//	v[1984] = (v[1195] + v[1199]) / 2e0;
//	v[1213] = v[1198] + v[1200];
//	v[1987] = v[1213] / 2e0;
//	v[1214] = v[1184] + v[1196];
//	v[1983] = v[1214] / 2e0;
//	v[1077] = v[1077] + v[1115] * v[132] + v[2704];
//	v[1070] = v[1070] + v[1116] * v[132] + v[2706];
//	v[1063] = v[1063] + v[1117] * v[132] + v[2708];
//	v[1078] = v[1078] + v[1115] * v[130] + v[2705];
//	v[1071] = v[1071] + v[1116] * v[130] + v[2707];
//	v[1064] = v[1064] + v[1117] * v[130] + v[2709];
//	v[1079] = v[1079] + v[1115] * v[128] - v[2704] - v[2705];
//	v[1072] = v[1072] + v[1116] * v[128] - v[2706] - v[2707];
//	v[1065] = v[1065] + v[1117] * v[128] - v[2708] - v[2709];
//	v[4534] = v[1065] - v[638] * v[872] - v[639] * v[873] - v[640] * v[874];
//	v[4535] = v[1072] - v[638] * v[876] - v[639] * v[877] - v[640] * v[878];
//	v[4536] = v[1079] - v[638] * v[880] - v[639] * v[881] - v[640] * v[882];
//	v[4537] = v[1064] - v[638] * v[884] - v[639] * v[885] - v[640] * v[886];
//	v[4538] = v[1071] - v[638] * v[888] - v[639] * v[889] - v[640] * v[890];
//	v[4539] = v[1078] - v[638] * v[892] - v[639] * v[893] - v[640] * v[894];
//	v[4540] = v[1063] - v[638] * v[896] - v[639] * v[897] - v[640] * v[898];
//	v[4541] = v[1070] - v[638] * v[900] - v[639] * v[901] - v[640] * v[902];
//	v[4542] = v[1077] - v[638] * v[904] - v[639] * v[905] - v[640] * v[906];
//	v[4543] = v[1059] - v[638] * v[908] - v[639] * v[909] - v[640] * v[910];
//	v[4544] = v[1066] - v[638] * v[912] - v[639] * v[913] - v[640] * v[914];
//	v[4545] = v[1073] - v[638] * v[916] - v[639] * v[917] - v[640] * v[918];
//	v[4546] = -v[1198] + v[1200] + 2e0*alphaB[0] * v[1982] + alphaB[2] * v[1984] + v[1183] * v[2824] + 2e0*
//		(v[1201] * v[2824] + v[1206] * v[2826]) + v[1214] * v[306] - v[638] * v[920] - v[639] * v[921] - v[640] * v[922];
//	v[4547] = v[1195] - v[1199] + (alphaB[2] * v[1213]) / 2e0 + (alphaB[0] * v[1214]) / 2e0 + (v[1179] * v[172]) / 2e0
//		+ 2e0*alphaB[1] * v[1986] + 2e0*(v[1202] * v[2824] - v[1208] * v[2826]) - v[638] * v[924] - v[639] * v[925]
//		- v[640] * v[926];
//	v[4548] = -v[1184] + v[1196] + alphaB[0] * v[1984] + 2e0*alphaB[2] * v[1988] + v[172] * v[2825] + 2e0*
//		(v[1197] * v[2824] + v[1210] * v[2826]) + v[1213] * v[306] - v[638] * v[928] - v[639] * v[929] - v[640] * v[930];
//	for (i933 = 1; i933 <= 15; i933++) {
//		i2737 = (i933 == 7 ? 1 : 0);
//		i2736 = (i933 == 4 ? 1 : 0);
//		i2735 = (i933 == 10 ? 1 : 0);
//		i2734 = (i933 == 8 ? 1 : 0);
//		i2733 = (i933 == 5 ? 1 : 0);
//		i2732 = (i933 == 11 ? 1 : 0);
//		i2731 = (i933 == 9 ? 1 : 0);
//		i2730 = (i933 == 6 ? 1 : 0);
//		i2729 = (i933 == 12 ? 1 : 0);
//		i2720 = (i933 == 3 ? 1 : 0);
//		i2719 = (i933 == 2 ? 1 : 0);
//		i2718 = (i933 == 1 ? 1 : 0);
//		i2712 = (i933 == 14 ? 1 : 0);
//		v[2804] = (*a4)*i2712;
//		v[2715] = i2712 * v[172];
//		i2711 = (i933 == 13 ? 1 : 0);
//		v[2805] = (*a4)*i2711;
//		v[2714] = -(i2711*v[172]);
//		i2710 = (i933 == 15 ? 1 : 0);
//		v[2806] = (*a4)*i2710;
//		v[2713] = -(i2710*v[172]);
//		v[1228] = v[4199 + i933];
//		v[1230] = v[4229 + i933];
//		v[1232] = v[4214 + i933];
//		v[1233] = v[4702 + i933];
//		v[1235] = v[4184 + i933];
//		v[1278] = 2e0*v[1235] * v[783];
//		v[1318] = -4e0*v[1278] * v[172];
//		v[2717] = v[1318] * v[227];
//		v[1293] = -2e0*v[1278];
//		v[1236] = v[4762 + i933];
//		v[1612] = -i2718 / 2e0;
//		v[1613] = -i2719 / 2e0;
//		v[1614] = -i2720 / 2e0;
//		v[1261] = i2710 + v[1228];
//		v[1262] = -i2710 + v[1228];
//		v[1263] = i2711 + v[1230];
//		v[1264] = -i2711 + v[1230];
//		v[1265] = -i2712 + v[1232];
//		v[1266] = i2712 + v[1232];
//		v[1268] = v[1235] * v[1267] + 8e0*i2710*v[783];
//		v[1269] = 2e0*alphaB[1] * i2712 - v[1235];
//		v[1271] = v[1235] * v[1270] - 8e0*i2712*v[783];
//		v[1273] = v[1235] * v[1272] + 8e0*i2711*v[783];
//		v[2716] = 2e0*(-v[2715] - (v[1235] * v[312]) / 2e0);
//		v[1275] = v[2714] + (v[1235] * v[311]) / 2e0;
//		v[1276] = v[2713] + (v[1235] * v[313]) / 2e0;
//		v[1277] = alphaB[2] * v[1278] + v[2713] / 2e0;
//		v[1279] = alphaB[0] * v[1278] + v[2714] / 2e0;
//		v[1280] = -(alphaB[1] * v[1278]) + v[2715] / 2e0;
//		v[1281] = (v[1233] * v[172]) / 2e0 - v[1278] * v[351];
//		v[1354] = v[1281] * v[238];
//		v[1282] = (-(v[1233] * v[313]) - v[1268] * v[351]) / 2e0;
//		v[1283] = (v[1269] * v[172]) / 2e0 - v[1278] * v[333];
//		v[1347] = v[1283] * v[236];
//		v[1284] = (v[1269] * v[312] + v[1271] * v[333]) / 2e0;
//		v[1285] = (v[1236] * v[172]) / 2e0 - v[1278] * v[314];
//		v[1339] = v[1285] * v[230];
//		v[1286] = (-(v[1236] * v[311]) - v[1273] * v[314]) / 2e0;
//		v[1287] = (v[2716] + v[1233] * v[312] + v[1271] * v[351]) / 2e0;
//		v[1288] = (v[2716] + v[1236] * v[312] + v[1271] * v[314]) / 2e0;
//		v[1289] = v[1275] - (v[1233] * v[311]) / 2e0 - (v[1273] * v[351]) / 2e0;
//		v[1290] = v[1275] - (v[1269] * v[311]) / 2e0 - (v[1273] * v[333]) / 2e0;
//		v[1291] = v[1293] - v[1263] * v[311] - v[1273] * v[346];
//		v[1292] = v[1263] * v[172] - 2e0*v[1278] * v[346];
//		v[1294] = -v[1293] + v[1265] * v[312] + v[1271] * v[342];
//		v[1295] = v[1265] * v[172] - 2e0*v[1278] * v[342];
//		v[1356] = v[1295] * v[230];
//		v[1296] = -v[1293] - v[1264] * v[311] - v[1273] * v[338];
//		v[1297] = v[1264] * v[172] - 2e0*v[1278] * v[338];
//		v[1348] = v[1297] * v[238];
//		v[1298] = v[1276] - (v[1269] * v[313]) / 2e0 - (v[1268] * v[333]) / 2e0;
//		v[1299] = v[1276] - (v[1236] * v[313]) / 2e0 - (v[1268] * v[314]) / 2e0;
//		v[1300] = v[1293] - v[1261] * v[313] - v[1268] * v[329];
//		v[1301] = v[1261] * v[172] - 2e0*v[1278] * v[329];
//		v[1302] = v[1293] + v[1266] * v[312] + v[1271] * v[325];
//		v[1303] = v[1266] * v[172] - 2e0*v[1278] * v[325];
//		v[1340] = v[1303] * v[238];
//		v[1304] = -v[1277] + v[1263] * v[312] + v[1271] * v[346];
//		v[1305] = -v[1277] - v[1265] * v[311] - v[1273] * v[342];
//		v[1306] = -v[1277] + v[1264] * v[312] + v[1271] * v[338];
//		v[1307] = -v[1277] - v[1266] * v[311] - v[1273] * v[325];
//		v[1308] = v[1318] + 2e0*v[1277] * v[324];
//		v[1309] = -(v[1287] * v[718]) - v[1304] * v[720] - v[1294] * v[722];
//		v[2744] = v[1309] * v[239];
//		v[1310] = -(v[1289] * v[718]) - v[1291] * v[720] - v[1305] * v[722];
//		v[2745] = v[1310] * v[229];
//		v[1311] = -v[1293] - v[1262] * v[313] - v[1268] * v[319];
//		v[1312] = v[1262] * v[172] - 2e0*v[1278] * v[319];
//		v[1313] = -v[1279] + v[1261] * v[312] + v[1271] * v[329];
//		v[1314] = -v[1279] - v[1265] * v[313] - v[1268] * v[342];
//		v[1315] = -v[1279] - v[1266] * v[313] - v[1268] * v[325];
//		v[1316] = -v[1279] + v[1262] * v[312] + v[1271] * v[319];
//		v[1317] = v[1277] * v[321] + v[1279] * v[324];
//		v[1319] = v[1318] + 2e0*v[1279] * v[321];
//		v[1381] = v[1024] * v[1319];
//		v[1320] = -(v[1306] * v[718]) - v[1284] * v[720] - v[1313] * v[722];
//		v[1321] = v[1280] - v[1263] * v[313] - v[1268] * v[346];
//		v[1322] = v[1280] - v[1264] * v[313] - v[1268] * v[338];
//		v[1323] = v[1280] - v[1261] * v[311] - v[1273] * v[329];
//		v[1324] = v[1280] - v[1262] * v[311] - v[1273] * v[319];
//		v[1325] = -(v[1279] * v[318]) - v[1280] * v[321];
//		v[1326] = -(v[1277] * v[318]) - v[1280] * v[324];
//		v[1327] = v[1318] + 2e0*v[1280] * v[318];
//		v[1385] = v[1020] * v[1327];
//		v[1328] = -(v[1282] * v[718]) - v[1321] * v[720] - v[1314] * v[722];
//		v[2746] = v[1328] * v[249];
//		v[2769] = -v[2744] - v[2745] - v[2746];
//		v[1329] = -(v[1322] * v[718]) - v[1298] * v[720] - v[1300] * v[722];
//		v[1330] = -(v[1296] * v[718]) - v[1290] * v[720] - v[1323] * v[722];
//		v[1523] = v[1330] * v[229] + v[1320] * v[239] + v[1329] * v[249];
//		v[1331] = -(v[1307] * v[718]) - v[1324] * v[720] - v[1286] * v[722];
//		v[1332] = -(v[1315] * v[718]) - v[1311] * v[720] - v[1299] * v[722];
//		v[1333] = -(v[1034] * v[1282]) - v[1039] * v[1287] - v[1038] * v[1289] - v[1030] * v[1296] - v[1025] * v[1302]
//			- v[1035] * v[1306] - v[1029] * v[1307] - v[1028] * v[1315] - v[1033] * v[1322];
//		v[1334] = -(v[1035] * v[1284]) - v[1030] * v[1290] - v[1038] * v[1291] - v[1033] * v[1298] - v[1039] * v[1304]
//			- v[1028] * v[1311] - v[1025] * v[1316] - v[1034] * v[1321] - v[1029] * v[1324];
//		v[1335] = -(v[1302] * v[718]) - v[1316] * v[720] - v[1288] * v[722];
//		v[1526] = v[1331] * v[229] + v[1335] * v[239] + v[1332] * v[249];
//		v[1336] = -(v[1029] * v[1286]) - v[1025] * v[1288] - v[1039] * v[1294] - v[1028] * v[1299] - v[1033] * v[1300]
//			- v[1038] * v[1305] - v[1035] * v[1313] - v[1034] * v[1314] - v[1030] * v[1323];
//		v[1337] = v[1339] + v[1312] * v[236];
//		v[1338] = v[1337] + v[1340] + v[2805];
//		v[1341] = v[1339] + v[1340];
//		v[1342] = QBi[0][0] * v[1285] + QBi[2][0] * v[1303] + QBi[1][0] * v[1312];
//		v[1343] = QBi[0][1] * v[1285] + QBi[2][1] * v[1303] + QBi[1][1] * v[1312];
//		v[1344] = QBi[0][2] * v[1285] + QBi[2][2] * v[1303] + QBi[1][2] * v[1312];
//		v[1345] = v[1347] + v[1301] * v[230];
//		v[1346] = v[1345] + v[1348] + v[2804];
//		v[1349] = v[1347] + v[1348];
//		v[1350] = QBi[1][0] * v[1283] + QBi[2][0] * v[1297] + QBi[0][0] * v[1301];
//		v[1351] = QBi[1][1] * v[1283] + QBi[2][1] * v[1297] + QBi[0][1] * v[1301];
//		v[1352] = QBi[1][2] * v[1283] + QBi[2][2] * v[1297] + QBi[0][2] * v[1301];
//		v[1353] = v[1354] + v[1356];
//		v[1355] = v[1354] + v[1292] * v[236];
//		v[1357] = v[1355] + v[1356] + v[2806];
//		v[1358] = QBi[2][0] * v[1281] + QBi[1][0] * v[1292] + QBi[0][0] * v[1295];
//		v[2752] = v[1342] * v[935] + v[1350] * v[936] + v[1358] * v[937];
//		v[1359] = QBi[2][1] * v[1281] + QBi[1][1] * v[1292] + QBi[0][1] * v[1295];
//		v[2753] = v[1343] * v[935] + v[1351] * v[936] + v[1359] * v[937];
//		v[1360] = QBi[2][2] * v[1281] + QBi[1][2] * v[1292] + QBi[0][2] * v[1295];
//		v[2754] = v[1344] * v[935] + v[1352] * v[936] + v[1360] * v[937];
//		v[1361] = v[1271] + v[1317];
//		v[1379] = v[1024] * v[1361];
//		v[1362] = -v[1271] + v[1317];
//		v[1364] = (v[1361] * v[186] + v[1292] * v[225] + v[1363] * v[2717]) / v[227];
//		v[1365] = v[1268] + v[1325];
//		v[1387] = v[1024] * v[1365];
//		v[1366] = -v[1268] + v[1325];
//		v[1383] = v[1020] * v[1366];
//		v[1368] = (v[1365] * v[182] + v[1297] * v[226] + v[1367] * v[2717]) / v[227];
//		v[1370] = (v[1366] * v[177] + v[1303] * v[2629] + v[1369] * v[2717]) / v[227];
//		v[1371] = v[1381] + v[1383];
//		v[2845] = v[1371] / v[227];
//		v[1372] = v[1273] + v[1326];
//		v[1377] = v[1020] * v[1372];
//		v[1373] = -v[1273] + v[1326];
//		v[1374] = v[1385] + v[1387];
//		v[2841] = v[1374] / v[227];
//		v[1375] = v[1016] * v[1308] + v[1377] + v[1379];
//		v[2849] = v[1375] / v[227];
//		v[1378] = v[1375] - v[1379];
//		v[1380] = v[1375] - v[1377];
//		v[2842] = v[1380] / v[227];
//		v[1382] = v[1016] * v[1362] + v[1381];
//		v[1384] = v[1382] + v[1383];
//		v[2843] = v[1384] / v[227];
//		v[1386] = v[1016] * v[1373] + v[1385];
//		v[1388] = v[1386] + v[1387];
//		v[2846] = v[1388] / v[227];
//		v[1389] = -i2735 + i2718 * v[128] + i2736 * v[130] + i2737 * v[132] - v[1342] * v[166] - v[1343] * v[167]
//			- v[1344] * v[168];
//		v[1391] = -i2732 + i2719 * v[128] + i2733 * v[130] + i2734 * v[132] - v[1350] * v[166] - v[1351] * v[167]
//			- v[1352] * v[168];
//		v[1393] = -i2729 + i2720 * v[128] + i2730 * v[130] + i2731 * v[132] - v[1358] * v[166] - v[1359] * v[167]
//			- v[1360] * v[168];
//		v[2721] = v[1389] * v[250] + v[1391] * v[251] + v[1393] * v[252];
//		v[1396] = -(v[1395] * v[2721] * v[967]);
//		v[1397] = v[2721] / v[597];
//		v[2722] = v[1397] * v[263];
//		v[1573] = v[1397];
//		v[1412] = v[1393] * v[264] + v[252] * v[2722];
//		v[1402] = v[1391] * v[264] + v[251] * v[2722];
//		v[2803] = v[1402] * v[265];
//		v[1399] = v[1389] * v[264] + v[250] * v[2722];
//		v[1398] = 2e0*v[1399] * v[265];
//		v[1512] = v[1398];
//		v[1400] = v[1399];
//		v[2723] = v[1400] * v[266] + v[2803];
//		v[1401] = 2e0*v[1402] * v[266];
//		v[1514] = v[1401];
//		v[1403] = v[132] * v[2723];
//		v[1404] = v[130] * v[2723];
//		v[1405] = v[128] * v[2723];
//		v[1406] = v[1400] * v[1407] + 2e0*v[1402] * v[989];
//		v[1408] = 2e0*v[1004] * v[1400] + v[1402] * v[1407];
//		v[1409] = v[1402];
//		v[1449] = v[1409];
//		v[1410] = 2e0*v[1412] * v[267];
//		v[1517] = v[1410];
//		v[1414] = (*a4)*v[5617 + i933] + 2e0*v[1412] * v[988];
//		v[1415] = v[1412];
//		v[1451] = v[1415];
//		v[1416] = 0e0;
//		v[1417] = 0e0;
//		v[1418] = 0e0;
//		v[1419] = 0e0;
//		v[1420] = 0e0;
//		v[1421] = 0e0;
//		v[1422] = 0e0;
//		v[1423] = 0e0;
//		v[1424] = 0e0;
//		v[1425] = 0e0;
//		v[1426] = 0e0;
//		v[1427] = 0e0;
//		v[1428] = 0e0;
//		v[1429] = 0e0;
//		v[1430] = 0e0;
//		v[1431] = 0e0;
//		v[1432] = 0e0;
//		v[1433] = 0e0;
//		v[1434] = 0e0;
//		v[1435] = 0e0;
//		v[1436] = 0e0;
//		v[1437] = 0e0;
//		v[1438] = 0e0;
//		v[1439] = 0e0;
//		v[1440] = 0e0;
//		v[1441] = 0e0;
//		v[1442] = 0e0;
//		v[1443] = 0e0;
//		v[1444] = 0e0;
//		v[1445] = 0e0;
//		v[1446] = 0e0;
//		v[1447] = 0e0;
//		b1448 = b537;
//		if (b1448) {
//			v[1450] = v[1415] * v[275] - v[1449] * v[276];
//			v[1452] = -(v[1451] * v[274]) + v[1400] * v[276];
//			v[1453] = v[1449] * v[274] - v[1400] * v[275];
//			v[2727] = v[1107] * v[1450] + v[1102] * v[1452] + v[1097] * v[1453];
//			v[2725] = v[1450] * v[539] + v[1452] * v[540] + v[1453] * v[541];
//			v[1454] = v[2725] / v[542];
//			v[2726] = v[1454] * v[3013];
//			v[1464] = v[1454] * v[2724];
//			v[1419] = -(v[1110] * v[1455] * v[2725]);
//			v[1456] = v[1450] * v[2653] + v[2726] * v[539];
//			v[1475] = 2e0*v[1456] * v[551];
//			v[1458] = v[1452] * v[2653] + v[2726] * v[540];
//			v[1472] = 2e0*v[1458] * v[552];
//			v[1459] = v[1453] * v[2653] + v[2726] * v[541];
//			v[1473] = 2e0*v[1459] * v[553];
//			v[1422] = v[1464] * v[2679] + v[2727] * v[544];
//			v[1421] = v[1454] * v[1462] * v[548];
//			v[1420] = v[1454] * v[2679] * v[549] + v[2727] * v[550];
//			v[1446] = 2e0*v[1464] * v[1465] * v[3014];
//			v[1447] = v[1109] * v[1454] * v[1460] * v[1465] * v[3015];
//			v[1418] = v[1097] * v[2726] + v[1453] * v[2728];
//			v[1417] = v[1102] * v[2726] + v[1452] * v[2728];
//			v[1416] = v[1107] * v[2726] + v[1450] * v[2728];
//			v[1467] = (v[1458] * v[551] + v[1456] * v[552]) / 2e0;
//			v[1468] = v[1472] + v[1475];
//			v[1469] = v[1468] + v[1473];
//			v[1470] = (v[1459] * v[551] + v[1456] * v[553]) / 2e0;
//			v[1471] = (v[1459] * v[552] + v[1458] * v[553]) / 2e0;
//			v[1474] = v[1472] + v[1473];
//			v[1476] = v[1473] + v[1475];
//			v[1425] = (v[1082] * v[1456] + v[1081] * v[1458] + 4e0*v[1459] * v[1477]) / 2e0;
//			v[1424] = (v[1083] * v[1456] + v[1081] * v[1459] + 4e0*v[1458] * v[1478]) / 2e0;
//			v[1423] = (v[1083] * v[1458] + v[1082] * v[1459] + 4e0*v[1456] * v[1479]) / 2e0;
//			v[1480] = -4e0*v[1469] * v[1789];
//			v[1445] = 8e0*v[1093] * v[1469] * v[3016];
//			v[1436] = (v[1040] * v[1480]) / 2e0;
//			v[1440] = (v[1044] * v[1480]) / 2e0;
//			v[1438] = v[1042] * v[1480];
//			v[1442] = v[1046] * v[1480];
//			v[1441] = v[1045] * v[1480];
//			v[1443] = v[1047] * v[1480];
//			v[1437] = v[1041] * v[1480];
//			v[1439] = v[1043] * v[1480];
//			v[1444] = (v[1048] * v[1480]) / 2e0;
//			v[1481] = -v[1459] + v[1467];
//			v[1482] = v[1459] + v[1467];
//			v[1483] = v[1458] + v[1470];
//			v[1484] = -v[1458] + v[1470];
//			v[1485] = -v[1456] + v[1471];
//			v[1486] = v[1456] + v[1471];
//			v[1428] = v[1085] * v[1480] + v[1481] * v[554];
//			v[1430] = v[1087] * v[1480] + v[1482] * v[554];
//			v[1429] = v[1086] * v[1480] + v[1483] * v[554];
//			v[1433] = v[1090] * v[1480] + v[1484] * v[554];
//			v[1427] = (v[1084] * v[1480] - v[1474] * v[554]) / 2e0;
//			v[1432] = v[1089] * v[1480] + v[1485] * v[554];
//			v[1434] = v[1091] * v[1480] + v[1486] * v[554];
//			v[1431] = (v[1088] * v[1480] - v[1476] * v[554]) / 2e0;
//			v[1435] = (v[1092] * v[1480] - v[1468] * v[554]) / 2e0;
//			v[1426] = -(v[1048] * v[1468]) / 2e0 - (v[1040] * v[1474]) / 2e0 - (v[1044] * v[1476]) / 2e0 + v[1041] * v[1481]
//				+ v[1043] * v[1482] + v[1042] * v[1483] + v[1046] * v[1484] + v[1045] * v[1485] + v[1047] * v[1486];
//		}
//		else {
//		};
//		v[1511] = v[1415];
//		v[1505] = v[1409];
//		v[1503] = v[1400];
//		v[1501] = v[1428];
//		v[1500] = v[1429];
//		v[1498] = v[1431];
//		v[1497] = v[1432];
//		v[1495] = v[1434];
//		v[1494] = v[1435];
//		v[1487] = 0e0;
//		v[1488] = 0e0;
//		v[1489] = 0e0;
//		v[1490] = 0e0;
//		v[1491] = 0e0;
//		v[1492] = 0e0;
//		b1493 = (*previouscontact);
//		if (b1493) {
//			v[1509] = v[1503] * v[963];
//			v[1508] = v[1415] * v[961];
//			v[1506] = v[1409] * v[962];
//			v[1435] = 0e0;
//			v[1434] = 0e0;
//			v[1496] = -i2729 + i2720 * v[133] + i2730 * v[135] + i2731 * v[137] + gti[0] * v[1433] + gti[2] * v[1494]
//				+ gti[1] * v[1495] - v[1358] * v[169] - v[1359] * v[170] - v[1360] * v[171];
//			v[2738] = -(v[1496] * v[267]) - v[1511] * v[574];
//			v[1433] = 0e0;
//			v[1432] = 0e0;
//			v[1431] = 0e0;
//			v[1499] = -i2732 + i2719 * v[133] + i2733 * v[135] + i2734 * v[137] + gti[0] * v[1430] + gti[2] * v[1497]
//				+ gti[1] * v[1498] - v[1350] * v[169] - v[1351] * v[170] - v[1352] * v[171];
//			v[2740] = -(v[1499] * v[266]) - v[1505] * v[573];
//			v[1430] = 0e0;
//			v[1429] = 0e0;
//			v[1428] = 0e0;
//			v[1502] = -i2735 + i2718 * v[133] + i2736 * v[135] + i2737 * v[137] + gti[0] * v[1427] + gti[2] * v[1500]
//				+ gti[1] * v[1501] - v[1342] * v[169] - v[1343] * v[170] - v[1344] * v[171];
//			v[2739] = -(v[1502] * v[265]) - v[1503] * v[572];
//			v[1427] = 0e0;
//			v[1504] = v[1506] + v[1509];
//			v[1507] = v[1506] + v[1508];
//			v[1510] = v[1508] + v[1509];
//			v[1492] = v[1496] * v[961];
//			v[1491] = v[1499] * v[962];
//			v[1490] = v[1502] * v[963];
//			v[1487] = v[1503] * v[1513] - v[1507] * v[265] - v[1512] * v[963];
//			v[1408] = v[1408] + v[1502] * v[1513] - v[1507] * v[572] + (v[2738] + v[2740])*v[963];
//			v[1488] = v[1505] * v[1516] - v[1510] * v[266] - v[1514] * v[962];
//			v[1406] = v[1406] + v[1499] * v[1516] - v[1510] * v[573] + (v[2738] + v[2739])*v[962];
//			v[1489] = v[1511] * v[1520] - v[1504] * v[267] - v[1517] * v[961];
//			v[1414] = v[1414] + v[1496] * v[1520] - v[1504] * v[574] + (v[2739] + v[2740])*v[961];
//		}
//		else {
//		};
//		v[2787] = -(v[1410] * v[966]);
//		v[2786] = -(v[1401] * v[965]);
//		v[2741] = v[1415] * v[966];
//		v[1569] = v[1415];
//		v[2751] = v[1569] * v[1591];
//		v[1566] = v[249] * v[2741];
//		v[1563] = v[239] * v[2741];
//		v[1560] = v[229] * v[2741];
//		v[2742] = v[1409] * v[965];
//		v[1567] = v[249] * v[2742];
//		v[1564] = v[239] * v[2742];
//		v[1561] = v[229] * v[2742];
//		v[1548] = v[1409];
//		v[2785] = -(v[1398] * v[964]);
//		v[1543] = v[1400];
//		v[2743] = v[1543] * v[964];
//		v[1555] = v[249] * v[2743];
//		v[1553] = v[239] * v[2743];
//		v[1551] = v[229] * v[2743];
//		v[1521] = (*a4)*i2729 + v[2744] + v[2745] + v[2746];
//		v[2748] = v[1411] * v[1412] - v[1521] * v[267] + (*a4)*v[4972 + i933];
//		v[1522] = v[1523] + (*a4)*v[5002 + i933];
//		v[1524] = (*a4)*i2732 + v[1523];
//		v[2750] = -(v[1524] * v[266]) + v[2748];
//		v[1525] = v[1526] + (*a4)*v[4987 + i933];
//		v[2749] = -(v[1525] * v[265]) - v[1522] * v[266];
//		v[1527] = (*a4)*i2735 + v[1526];
//		v[2747] = -(v[1527] * v[265]) + v[2748];
//		v[1532] = v[1318] * v[1528] + v[1319] * v[1529] + v[1361] * v[1530] + v[1365] * v[1531] + v[1338] * v[1992]
//			+ v[1364] * v[236] + v[1368] * v[238] + v[1345] * v[2630] + v[1353] * v[2633] + (*a4)*v[4927 + i933];
//		v[1537] = v[1220] * v[1337] + v[1219] * v[1355] + v[1318] * v[1533] + v[1327] * v[1534] + v[1366] * v[1535]
//			+ v[1372] * v[1536] + v[1346] * v[1993] + v[1364] * v[230] + v[1370] * v[238] + (*a4)*v[4942 + i933];
//		v[1542] = v[1222] * v[1341] + v[1221] * v[1349] + v[1318] * v[1538] + v[1308] * v[1539] + v[1362] * v[1540]
//			+ v[1373] * v[1541] + v[1357] * v[1994] + v[1368] * v[230] + v[1370] * v[236] + (*a4)*v[4957 + i933];
//		v[1544] = v[2742] + v[2743];
//		v[1550] = v[2741] + v[2743];
//		v[1552] = v[1551] + v[1560];
//		v[1554] = v[1553] + v[1563];
//		v[1556] = v[1555] + v[1566];
//		v[1559] = v[2741] + v[2742];
//		v[1562] = v[1560] + v[1561];
//		v[1565] = v[1563] + v[1564];
//		v[1568] = v[1566] + v[1567];
//		v[1570] = 0e0;
//		b1571 = b598;
//		if (b1571) {
//			b1572 = b600;
//			if (b1572) {
//				v[1570] = -v[1573];
//				v[1400] = 0e0;
//				v[1409] = 0e0;
//				v[1415] = 0e0;
//			}
//			else {
//			};
//		}
//		else {
//		};
//		v[1578] = v[1023] * v[1532] + v[1019] * v[1537] + v[1015] * v[1542] + v[1410] * v[1574] + v[1543] * v[1575]
//			+ v[1548] * v[1576] + v[1569] * v[1577] + v[267] * v[2749] + (*a4)*v[5017 + i933] + v[2769] * v[587];
//		v[1583] = v[1022] * v[1532] + v[1018] * v[1537] + v[1014] * v[1542] + v[1401] * v[1579] + v[1543] * v[1580]
//			+ v[1548] * v[1581] + v[1569] * v[1582] + v[1405] * v[212] + v[1404] * v[215] + v[1403] * v[218] + v[266] * v[2747] + (*a4
//				)*v[5032 + i933] - v[1523] * v[582];
//		v[1606] = v[1583];
//		v[1406] = v[1406] + v[1522] * v[1990] - v[1524] * v[1991] + (v[1543] * v[1584] + v[2747] + v[2751])*v[965];
//		v[1737] = v[1406];
//		v[1589] = v[1021] * v[1532] + v[1017] * v[1537] + v[1013] * v[1542] + v[1398] * v[1585] + v[1543] * v[1586]
//			+ v[1548] * v[1587] + v[1569] * v[1588] + v[1405] * v[213] + v[1404] * v[216] + v[1403] * v[219] + v[265] * v[2750] + (*a4
//				)*v[5047 + i933] - v[1526] * v[577];
//		v[1605] = v[1589];
//		v[1414] = v[1414] + (v[1548] * v[1579] + v[1543] * v[1585] + v[2749])*v[966] - v[1521] * v[990];
//		v[1736] = v[1414];
//		v[1408] = v[1408] - v[1527] * v[1989] + v[1525] * v[1990] + (v[1548] * v[1590] + v[2750] + v[2751])*v[964];
//		v[1738] = v[1408];
//		v[1592] = 0e0;
//		v[1593] = 0e0;
//		v[1594] = 0e0;
//		v[1595] = 0e0;
//		v[1596] = 0e0;
//		v[1597] = 0e0;
//		v[1598] = 0e0;
//		v[1599] = 0e0;
//		b1600 = b598;
//		if (b1600) {
//			v[1610] = v[1605] * v[970] + v[1606] * v[971] + v[1578] * v[972];
//			b1601 = b600;
//			if (b1601) {
//				v[1604] = v[1570] * v[1723] * (*zetan);
//				v[1603] = v[1604] * v[1724];
//				v[1598] = -((v[1603] * v[982]) / v[976]);
//				v[1595] = v[1602] * v[1604] * v[3047] * v[982];
//				v[1570] = 0e0;
//				v[1589] = 0e0;
//				v[1583] = 0e0;
//				v[1596] = v[1610];
//				v[1578] = 0e0;
//			}
//			else {
//				v[1608] = v[1397] * v[1729] * (*zetan);
//				v[1603] = v[1608] * v[1730];
//				v[1599] = -((v[1603] * v[982]) / v[981]);
//				v[1396] = v[1396] + v[1607] * v[1608] * v[3052] * v[982];
//				v[1397] = 0e0;
//				v[1589] = 0e0;
//				v[1583] = 0e0;
//				v[1597] = v[1610];
//				v[1578] = 0e0;
//			};
//			v[1594] = v[1603] * v[972];
//			v[1593] = v[1603] * v[971];
//			v[1592] = v[1603] * v[970];
//		}
//		else {
//		};
//		v[1728] = v[1592];
//		v[1726] = v[1593];
//		v[1725] = v[1594];
//		v[1611] = (*ct)*((i2736 / 2e0 + v[1612])*v[935] + (i2733 / 2e0 + v[1613])*v[936] + (i2730 / 2e0 + v[1614])*v[937]);
//		v[1615] = (*ct)*((i2737 / 2e0 + v[1612])*v[935] + (i2734 / 2e0 + v[1613])*v[936] + (i2731 / 2e0 + v[1614])*v[937]);
//		v[1616] = (*ct)*(v[2752] * v[296] + v[2753] * v[298] + v[2754] * v[300]);
//		v[1617] = (*ct)*(v[2752] * v[288] + v[2753] * v[289] + v[2754] * v[290]);
//		v[1618] = -(i2718*v[872]) - i2719 * v[876] - i2720 * v[880] - i2736 * v[884] - i2733 * v[888] - i2730 * v[892]
//			- i2737 * v[896] - i2734 * v[900] - i2731 * v[904] - i2735 * v[908] - i2732 * v[912] - i2729 * v[916] - i2711 * v[920]
//			- i2712 * v[924] - i2710 * v[928];
//		v[1634] = v[1618];
//		v[1619] = -(i2718*v[873]) - i2719 * v[877] - i2720 * v[881] - i2736 * v[885] - i2733 * v[889] - i2730 * v[893]
//			- i2737 * v[897] - i2734 * v[901] - i2731 * v[905] - i2735 * v[909] - i2732 * v[913] - i2729 * v[917] - i2711 * v[921]
//			- i2712 * v[925] - i2710 * v[929];
//		v[1632] = v[1619];
//		v[1620] = -(i2718*v[874]) - i2719 * v[878] - i2720 * v[882] - i2736 * v[886] - i2733 * v[890] - i2730 * v[894]
//			- i2737 * v[898] - i2734 * v[902] - i2731 * v[906] - i2735 * v[910] - i2732 * v[914] - i2729 * v[918] - i2711 * v[922]
//			- i2712 * v[926] - i2710 * v[930];
//		v[1630] = v[1620];
//		v[1621] = 0e0;
//		v[1622] = 0e0;
//		v[1623] = 0e0;
//		v[1624] = 0e0;
//		v[1625] = 0e0;
//		v[1626] = 0e0;
//		b1627 = b598;
//		if (b1627) {
//			b1628 = (*stick);
//			if (b1628) {
//				b1629 = b636;
//				if (b1629) {
//					v[1626] = v[1620];
//					v[1620] = 0e0;
//					v[1625] = v[1619];
//					v[1619] = 0e0;
//					v[1624] = v[1618];
//					v[1618] = 0e0;
//				}
//				else {
//					v[2755] = (*mud)*(v[1634] * v[631] + v[1632] * v[632] + v[1630] * v[633]);
//					v[1620] = 0e0;
//					v[1619] = 0e0;
//					v[2757] = v[1645] * v[2755] * v[647];
//					v[1618] = 0e0;
//					v[2756] = v[1643] * v[2755] * v[646] * v[651];
//					v[1626] = v[1630] * v[1642] + v[2756] * v[633];
//					v[1625] = v[1632] * v[1642] + v[2756] * v[632];
//					v[1624] = v[1634] * v[1642] + v[2756] * v[631];
//					v[1623] = v[2757] * v[624];
//					v[1622] = v[2757] * v[623];
//					v[1621] = v[2757] * v[622];
//				};
//			}
//			else {
//				b1646 = b670;
//				if (b1646) {
//					v[1626] = v[1630];
//					v[1620] = 0e0;
//					v[1625] = v[1632];
//					v[1619] = 0e0;
//					v[1624] = v[1634];
//					v[1618] = 0e0;
//				}
//				else {
//					v[1652] = v[1634] * v[2758];
//					v[1650] = v[1632] * v[2758];
//					v[1649] = v[1630] * v[2758];
//					v[1620] = 0e0;
//					v[1619] = 0e0;
//					v[2760] = (*mud)*v[1655] * (v[1634] * v[631] + v[1632] * v[632] + v[1630] * v[633])*v[677];
//					v[1618] = 0e0;
//					v[2759] = v[1651] * (v[1652] * v[631] + v[1650] * v[632] + v[1649] * v[633])*v[676];
//					v[1626] = v[2759] * v[633] + v[1649] * v[677];
//					v[1625] = v[2759] * v[632] + v[1650] * v[677];
//					v[1624] = v[2759] * v[631] + v[1652] * v[677];
//					v[1623] = v[2760] * v[624];
//					v[1622] = v[2760] * v[623];
//					v[1621] = v[2760] * v[622];
//				};
//			};
//		}
//		else {
//		};
//		v[2765] = -((*ct)*v[1624]);
//		v[2763] = -((*ct)*v[1625]);
//		v[2761] = -((*ct)*v[1626]);
//		v[1657] = -(i2710*v[640]) - (*ct)*(v[1626] * v[249] + v[1542] * v[937]);
//		v[1658] = -(i2712*v[640]) - (*ct)*(v[1626] * v[239] + v[1537] * v[937]);
//		v[1659] = -(i2711*v[640]) - (*ct)*(v[1626] * v[229] + v[1532] * v[937]);
//		v[1660] = v[223] * v[2761] + i2729 * v[2762];
//		v[1662] = v[222] * v[2761] + i2732 * v[2762];
//		v[1663] = v[221] * v[2761] + i2735 * v[2762];
//		v[1664] = v[220] * v[2761] + i2731 * v[2762];
//		v[1665] = v[219] * v[2761] + i2734 * v[2762];
//		v[1666] = v[218] * v[2761] + i2737 * v[2762];
//		v[1667] = v[217] * v[2761] + i2730 * v[2762];
//		v[1668] = v[216] * v[2761] + i2733 * v[2762];
//		v[1669] = v[215] * v[2761] + i2736 * v[2762];
//		v[1670] = v[214] * v[2761] + i2720 * v[2762];
//		v[1671] = v[213] * v[2761] + i2719 * v[2762];
//		v[1672] = v[212] * v[2761] + i2718 * v[2762];
//		v[1673] = -(i2710*v[639]) - (*ct)*(v[1625] * v[249] + v[1542] * v[936]);
//		v[1674] = -(i2712*v[639]) - (*ct)*(v[1625] * v[239] + v[1537] * v[936]);
//		v[1675] = -(i2711*v[639]) - (*ct)*(v[1625] * v[229] + v[1532] * v[936]);
//		v[1676] = v[223] * v[2763] + i2729 * v[2764];
//		v[1678] = v[222] * v[2763] + i2732 * v[2764];
//		v[1679] = v[221] * v[2763] + i2735 * v[2764];
//		v[1680] = v[220] * v[2763] + i2731 * v[2764];
//		v[1681] = v[219] * v[2763] + i2734 * v[2764];
//		v[1682] = v[218] * v[2763] + i2737 * v[2764];
//		v[1683] = v[217] * v[2763] + i2730 * v[2764];
//		v[1684] = v[216] * v[2763] + i2733 * v[2764];
//		v[1685] = v[215] * v[2763] + i2736 * v[2764];
//		v[1686] = v[214] * v[2763] + i2720 * v[2764];
//		v[1687] = v[213] * v[2763] + i2719 * v[2764];
//		v[1688] = v[212] * v[2763] + i2718 * v[2764];
//		v[1689] = -(i2710*v[638]) - (*ct)*(v[1624] * v[249] + v[1542] * v[935]);
//		v[1690] = -(i2712*v[638]) - (*ct)*(v[1624] * v[239] + v[1537] * v[935]);
//		v[1691] = -(i2711*v[638]) - (*ct)*(v[1624] * v[229] + v[1532] * v[935]);
//		v[1692] = v[223] * v[2765] + i2729 * v[2766];
//		v[1694] = v[222] * v[2765] + i2732 * v[2766];
//		v[1695] = v[221] * v[2765] + i2735 * v[2766];
//		v[1696] = v[220] * v[2765] + i2731 * v[2766];
//		v[1697] = v[219] * v[2765] + i2734 * v[2766];
//		v[1698] = v[218] * v[2765] + i2737 * v[2766];
//		v[1699] = v[217] * v[2765] + i2730 * v[2766];
//		v[1700] = v[216] * v[2765] + i2733 * v[2766];
//		v[1701] = v[215] * v[2765] + i2736 * v[2766];
//		v[1702] = v[214] * v[2765] + i2720 * v[2766];
//		v[1703] = v[213] * v[2765] + i2719 * v[2766];
//		v[1704] = v[212] * v[2765] + i2718 * v[2766];
//		v[1705] = (*epst)*v[1626];
//		v[1706] = (*epst)*v[1625];
//		v[1707] = (*epst)*v[1624];
//		v[1708] = v[1623];
//		v[1709] = v[1623];
//		v[1710] = v[1622];
//		v[1711] = v[1622];
//		v[1712] = v[1621];
//		v[1713] = v[1621];
//		b1714 = b598;
//		if (b1714) {
//			v[1715] = 0e0;
//			v[1716] = 0e0;
//			v[1717] = 0e0;
//			b1718 = b617;
//			if (b1718) {
//				v[1717] = v[1708];
//				v[1708] = 0e0;
//				v[1716] = v[1710];
//				v[1710] = 0e0;
//				v[1715] = v[1712];
//				v[1712] = 0e0;
//			}
//			else {
//				v[1709] = 0e0;
//				v[1708] = 0e0;
//				v[1711] = 0e0;
//				v[1710] = 0e0;
//				v[1713] = 0e0;
//				v[1712] = 0e0;
//			};
//			v[1727] = v[1715] * v[578] + v[1716] * v[586] + v[1717] * v[588];
//			b1722 = b600;
//			if (b1722) {
//				v[1594] = v[1594] + v[1717] * v[612];
//				v[1593] = v[1593] + v[1716] * v[612];
//				v[1596] = v[1596] + v[1727];
//				v[1592] = v[1592] + v[1715] * v[612];
//				v[1598] = v[1598] + 2e0*v[1596] * (*zetan);
//				v[1595] = v[1595] + (v[1598] * v[1723] * v[1724]) / 2e0;
//			}
//			else {
//				v[1594] = v[1725] + v[1717] * v[616];
//				v[1593] = v[1726] + v[1716] * v[616];
//				v[1597] = v[1597] + v[1727];
//				v[1592] = v[1728] + v[1715] * v[616];
//				v[1599] = v[1599] + 2e0*v[1597] * (*zetan);
//				v[1396] = v[1396] + (v[1599] * v[1729] * v[1730]) / 2e0;
//			};
//		}
//		else {
//		};
//		v[2834] = v[1592] * v[577];
//		v[2770] = v[1592] * v[265];
//		v[2835] = v[1593] * v[582];
//		v[2817] = -(v[1548] * v[1990]) - v[2786] + v[2835];
//		v[2772] = v[1593] * v[266];
//		v[2837] = v[1594] * v[587];
//		v[2771] = v[1594] * v[267];
//		v[1739] = v[1396];
//		v[1735] = v[1713];
//		v[1734] = v[1711];
//		v[1733] = v[1709];
//		b1731 = b598;
//		if (b1731) {
//			v[2768] = -(v[1735] * v[265]) - v[1734] * v[266] - v[1733] * v[267];
//			b1732 = b600;
//			if (b1732) {
//				v[1414] = v[1414] + v[1709] * v[2767];
//				v[1709] = 0e0;
//				v[1406] = v[1406] + v[1711] * v[2767];
//				v[1711] = 0e0;
//				v[1408] = v[1408] + v[1713] * v[2767];
//				v[1713] = 0e0;
//				v[1595] = v[1595] + v[2658] * v[2768] * v[3090];
//				v[1396] = v[1396] - v[1595];
//			}
//			else {
//				v[1414] = v[1736] + v[1733] * v[606];
//				v[1709] = 0e0;
//				v[1406] = v[1737] + v[1734] * v[606];
//				v[1711] = 0e0;
//				v[1408] = v[1738] + v[1735] * v[606];
//				v[1713] = 0e0;
//				v[1396] = v[1739] + (*n2)*v[2768] * v[3094] * v[594];
//			};
//		}
//		else {
//		};
//		v[1740] = v[1594] * v[249] + v[1542] * v[966];
//		v[2801] = v[1740] * v[267];
//		v[1741] = v[1594] * v[239] + v[1537] * v[966];
//		v[2798] = v[1741] * v[267];
//		v[1742] = v[1594] * v[229] + v[1532] * v[966];
//		v[2796] = v[1742] * v[267];
//		v[1743] = v[1594] * v[2398] + (*a4)*v[5602 + i933] + v[2769] * v[966];
//		v[1744] = v[2741] + v[2771];
//		v[2833] = v[1744] * v[266] + v[2817];
//		v[2816] = -(v[1543] * v[1990]) + v[1744] * v[265] - v[2785] + v[2834];
//		v[1745] = v[1593] * v[249] + v[1542] * v[965];
//		v[2802] = v[1745] * v[266];
//		v[1746] = v[1593] * v[239] + v[1537] * v[965];
//		v[2799] = v[1746] * v[266];
//		v[1747] = v[1593] * v[229] + v[1532] * v[965];
//		v[2797] = -(v[1747] * v[266]);
//		v[1748] = v[1593] * v[2407] + (*a4)*v[5677 + i933] - v[1523] * v[965];
//		v[1749] = v[1592] * v[249] + v[1542] * v[964];
//		v[2795] = v[1749] * v[265];
//		v[2790] = v[1555] + v[1567] + v[2795] + v[2802];
//		v[1750] = v[1592] * v[239] + v[1537] * v[964];
//		v[2794] = v[1750] * v[265];
//		v[2791] = v[1553] + v[1564] + v[2794] + v[2799];
//		v[1751] = v[1592] * v[229] + v[1532] * v[964];
//		v[2793] = -(v[1751] * v[265]);
//		v[2792] = -v[1551] - v[1561] + v[2793] + v[2797];
//		v[1414] = v[1414] + v[1594] * v[1752] - v[223] * (v[2770] + v[2772]);
//		v[1753] = v[1593] * v[212] + v[1592] * v[213] + (*a4)*v[5632 + i933];
//		v[1754] = v[1593] * v[215] + v[1592] * v[216] + (*a4)*v[5647 + i933];
//		v[1755] = v[1593] * v[218] + v[1592] * v[219] + (*a4)*v[5662 + i933];
//		v[2800] = v[128] * v[1753] + v[130] * v[1754] + v[132] * v[1755];
//		v[1406] = v[1406] + v[1593] * v[2406] - v[222] * (v[2770] + v[2771]);
//		v[1756] = v[1592] * v[2416] + (*a4)*v[5692 + i933] - v[1526] * v[964];
//		v[1408] = v[1408] + v[1592] * v[2415] - v[221] * (v[2771] + v[2772]);
//		v[1757] = 0e0;
//		v[1758] = 0e0;
//		v[1759] = 0e0;
//		v[1760] = 0e0;
//		v[1761] = 0e0;
//		v[1762] = 0e0;
//		v[1763] = 0e0;
//		v[1764] = 0e0;
//		v[1765] = 0e0;
//		b1766 = (*previouscontact);
//		if (b1766) {
//			v[2775] = v[1707] * v[265];
//			v[2774] = v[1706] * v[266];
//			v[2776] = v[2774] + v[2775];
//			v[2773] = v[1705] * v[267];
//			v[2778] = v[2773] + v[2774];
//			v[2777] = v[2773] + v[2775];
//			v[1492] = v[1492] + v[1705] * v[574];
//			v[1491] = v[1491] + v[1706] * v[573];
//			v[1487] = v[1487] + v[1050] * v[1707] - v[265] * v[2778];
//			v[1488] = v[1488] + v[1052] * v[1706] - v[266] * v[2777];
//			v[1489] = v[1489] + v[1054] * v[1705] - v[267] * v[2776];
//			v[1490] = v[1490] + v[1707] * v[572];
//			v[1414] = v[1414] + v[1705] * v[2676] - v[2776] * v[574];
//			v[1406] = v[1406] + v[1706] * v[2677] - v[2777] * v[573];
//			v[1408] = v[1408] - v[1707] * v[2678] - v[2778] * v[572];
//			v[1757] = gti[0] * v[1487];
//			v[1758] = gti[1] * v[1487];
//			v[1759] = gti[2] * v[1487];
//			v[1767] = -v[1487];
//			v[1768] = -(v[1487] * v[171]);
//			v[1769] = -(v[1487] * v[170]);
//			v[1770] = -(v[1487] * v[169]);
//			v[1771] = v[137] * v[1487];
//			v[1772] = v[135] * v[1487];
//			v[1773] = v[133] * v[1487];
//			v[1760] = gti[0] * v[1488];
//			v[1761] = gti[1] * v[1488];
//			v[1762] = gti[2] * v[1488];
//			v[1774] = -v[1488];
//			v[1775] = -(v[1488] * v[171]);
//			v[1776] = -(v[1488] * v[170]);
//			v[1777] = -(v[1488] * v[169]);
//			v[1778] = v[137] * v[1488];
//			v[1779] = v[135] * v[1488];
//			v[1780] = v[133] * v[1488];
//			v[1763] = gti[0] * v[1489];
//			v[1764] = gti[1] * v[1489];
//			v[1765] = gti[2] * v[1489];
//			v[1781] = -v[1489];
//			v[1782] = -(v[1489] * v[171]);
//			v[1783] = -(v[1489] * v[170]);
//			v[1784] = -(v[1489] * v[169]);
//			v[1785] = v[137] * v[1489];
//			v[1786] = v[135] * v[1489];
//			v[1787] = v[133] * v[1489];
//			v[1756] = -v[1490] + v[1756];
//			v[1748] = -v[1491] + v[1748];
//			v[1743] = -v[1492] + v[1743];
//		}
//		else {
//			v[1773] = 0e0;
//			v[1772] = 0e0;
//			v[1771] = 0e0;
//			v[1780] = 0e0;
//			v[1779] = 0e0;
//			v[1778] = 0e0;
//			v[1787] = 0e0;
//			v[1786] = 0e0;
//			v[1785] = 0e0;
//			v[1770] = 0e0;
//			v[1769] = 0e0;
//			v[1768] = 0e0;
//			v[1777] = 0e0;
//			v[1776] = 0e0;
//			v[1775] = 0e0;
//			v[1784] = 0e0;
//			v[1783] = 0e0;
//			v[1782] = 0e0;
//			v[1767] = 0e0;
//			v[1774] = 0e0;
//			v[1781] = 0e0;
//		};
//		b1788 = b537;
//		if (b1788) {
//			v[1444] = v[1444] + (v[1765] * v[554]) / 2e0;
//			v[1443] = v[1443] + v[1764] * v[554];
//			v[1442] = v[1442] + v[1763] * v[554];
//			v[1441] = v[1441] + v[1762] * v[554];
//			v[1440] = v[1440] + (v[1761] * v[554]) / 2e0;
//			v[1439] = v[1439] + v[1760] * v[554];
//			v[1438] = v[1438] + v[1759] * v[554];
//			v[1437] = v[1437] + v[1758] * v[554];
//			v[1426] = v[1426] + (v[1084] * v[1757]) / 2e0 + v[1085] * v[1758] + v[1086] * v[1759] + v[1087] * v[1760] +
//				(v[1088] * v[1761]) / 2e0 + v[1089] * v[1762] + v[1090] * v[1763] + v[1091] * v[1764] + (v[1092] * v[1765]) / 2e0;
//			v[1436] = v[1436] + (v[1757] * v[554]) / 2e0;
//			v[1445] = v[1445] - 4e0*v[1426] * v[1789];
//			v[2779] = v[1436] - v[1445];
//			v[2781] = (v[1438] + v[1442]) / 2e0;
//			v[2780] = (v[1441] + v[1443]) / 2e0;
//			v[1425] = v[1425] - v[1437] + v[1439] + v[2781] * v[551] + v[2780] * v[552] - 2e0*(v[1440] + v[2779])*v[553];
//			v[2782] = (v[1437] + v[1439]) / 2e0;
//			v[1424] = v[1424] + v[1438] - v[1442] + v[2782] * v[551] - 2e0*(v[1444] + v[2779])*v[552] + v[2780] * v[553];
//			v[1423] = v[1423] - v[1441] + v[1443] - 2e0*(v[1440] + v[1444] - v[1445])*v[551] + v[2782] * v[552]
//				+ v[2781] * v[553];
//			v[2783] = v[1423] * v[539] + v[1424] * v[540] + v[1425] * v[541];
//			v[1422] = v[1422] + v[2783] * v[544];
//			v[1420] = v[1420] + v[2783] * v[550];
//			v[1421] = v[1421] + v[1422] * v[549];
//			v[1446] = v[1446] + 2e0*v[1420] * v[1460];
//			v[1447] = v[1447] + (v[1446] * v[1461]) / 2e0;
//			v[1419] = v[1419] + v[1421] + v[1447];
//			v[2784] = v[1419] / v[542];
//			v[1418] = v[1418] + v[1425] * v[2653] + v[2784] * v[541];
//			v[1417] = v[1417] + v[1424] * v[2653] + v[2784] * v[540];
//			v[1416] = v[1416] + v[1423] * v[2653] + v[2784] * v[539];
//			v[1408] = v[1408] - v[1418] * v[275] + v[1417] * v[276];
//			v[1414] = v[1414] - v[1417] * v[274] + v[1416] * v[275];
//			v[1406] = v[1406] + v[1418] * v[274] - v[1416] * v[276];
//		}
//		else {
//		};
//		v[1796] = -(v[1704] * v[507]) - v[1703] * v[508] - v[1702] * v[509] - v[1701] * v[510] - v[1700] * v[511]
//			- v[1699] * v[512] - v[1698] * v[513] - v[1697] * v[514] - v[1696] * v[515] - v[1695] * v[516] - v[1694] * v[517]
//			- v[1692] * v[518] - v[1691] * v[519] - v[1690] * v[520] - v[1689] * v[521];
//		v[1797] = -(v[1704] * v[522]) - v[1703] * v[523] - v[1702] * v[524] - v[1701] * v[525] - v[1700] * v[526]
//			- v[1699] * v[527] - v[1698] * v[528] - v[1697] * v[529] - v[1696] * v[530] - v[1695] * v[531] - v[1694] * v[532]
//			- v[1692] * v[533] - v[1691] * v[534] - v[1690] * v[535] - v[1689] * v[536];
//		v[2822] = (v[1704] * v[492] + v[1703] * v[493] + v[1702] * v[494] + v[1701] * v[495] + v[1700] * v[496] + v[1699] * v[497]
//			+ v[1698] * v[498] + v[1697] * v[499] + v[1696] * v[500] + v[1695] * v[501] + v[1694] * v[502] + v[1692] * v[503]
//			+ v[1691] * v[504] + v[1690] * v[505] + v[1689] * v[506]) / 2e0;
//		v[2823] = (v[1704] * v[477] + v[1703] * v[478] + v[1702] * v[479] + v[1701] * v[480] + v[1700] * v[481] + v[1699] * v[482]
//			+ v[1698] * v[483] + v[1697] * v[484] + v[1696] * v[485] + v[1695] * v[486] + v[1694] * v[487] + v[1692] * v[488]
//			+ v[1691] * v[489] + v[1690] * v[490] + v[1689] * v[491]) / 2e0;
//		v[1800] = -(v[1688] * v[507]) - v[1687] * v[508] - v[1686] * v[509] - v[1685] * v[510] - v[1684] * v[511]
//			- v[1683] * v[512] - v[1682] * v[513] - v[1681] * v[514] - v[1680] * v[515] - v[1679] * v[516] - v[1678] * v[517]
//			- v[1676] * v[518] - v[1675] * v[519] - v[1674] * v[520] - v[1673] * v[521];
//		v[1801] = -(v[1688] * v[522]) - v[1687] * v[523] - v[1686] * v[524] - v[1685] * v[525] - v[1684] * v[526]
//			- v[1683] * v[527] - v[1682] * v[528] - v[1681] * v[529] - v[1680] * v[530] - v[1679] * v[531] - v[1678] * v[532]
//			- v[1676] * v[533] - v[1675] * v[534] - v[1674] * v[535] - v[1673] * v[536];
//		v[2820] = (v[1688] * v[492] + v[1687] * v[493] + v[1686] * v[494] + v[1685] * v[495] + v[1684] * v[496] + v[1683] * v[497]
//			+ v[1682] * v[498] + v[1681] * v[499] + v[1680] * v[500] + v[1679] * v[501] + v[1678] * v[502] + v[1676] * v[503]
//			+ v[1675] * v[504] + v[1674] * v[505] + v[1673] * v[506]) / 2e0;
//		v[2821] = (v[1688] * v[477] + v[1687] * v[478] + v[1686] * v[479] + v[1685] * v[480] + v[1684] * v[481] + v[1683] * v[482]
//			+ v[1682] * v[483] + v[1681] * v[484] + v[1680] * v[485] + v[1679] * v[486] + v[1678] * v[487] + v[1676] * v[488]
//			+ v[1675] * v[489] + v[1674] * v[490] + v[1673] * v[491]) / 2e0;
//		v[1804] = -(v[1672] * v[507]) - v[1671] * v[508] - v[1670] * v[509] - v[1669] * v[510] - v[1668] * v[511]
//			- v[1667] * v[512] - v[1666] * v[513] - v[1665] * v[514] - v[1664] * v[515] - v[1663] * v[516] - v[1662] * v[517]
//			- v[1660] * v[518] - v[1659] * v[519] - v[1658] * v[520] - v[1657] * v[521];
//		v[1805] = -(v[1672] * v[522]) - v[1671] * v[523] - v[1670] * v[524] - v[1669] * v[525] - v[1668] * v[526]
//			- v[1667] * v[527] - v[1666] * v[528] - v[1665] * v[529] - v[1664] * v[530] - v[1663] * v[531] - v[1662] * v[532]
//			- v[1660] * v[533] - v[1659] * v[534] - v[1658] * v[535] - v[1657] * v[536];
//		v[2818] = (v[1672] * v[492] + v[1671] * v[493] + v[1670] * v[494] + v[1669] * v[495] + v[1668] * v[496] + v[1667] * v[497]
//			+ v[1666] * v[498] + v[1665] * v[499] + v[1664] * v[500] + v[1663] * v[501] + v[1662] * v[502] + v[1660] * v[503]
//			+ v[1659] * v[504] + v[1658] * v[505] + v[1657] * v[506]) / 2e0;
//		v[2819] = (v[1672] * v[477] + v[1671] * v[478] + v[1670] * v[479] + v[1669] * v[480] + v[1668] * v[481] + v[1667] * v[482]
//			+ v[1666] * v[483] + v[1665] * v[484] + v[1664] * v[485] + v[1663] * v[486] + v[1662] * v[487] + v[1660] * v[488]
//			+ v[1659] * v[489] + v[1658] * v[490] + v[1657] * v[491]) / 2e0;
//		v[1814] = v[1013] * v[1592] + v[1014] * v[1593] + v[1015] * v[1594] + v[1328] * v[1808] + v[1329] * v[1809]
//			+ v[1332] * v[1810] + v[1569] * v[3131] + v[1543] * v[3132] + v[1548] * v[3133] + v[2785] * v[403] + v[2786] * v[406]
//			+ v[2787] * v[409] - v[1611] * v[491] - v[1615] * v[506] + v[1617] * v[521] + v[1616] * v[536] + (*ct)*(-
//			(v[1624] * v[928]) - v[1625] * v[929] - v[1626] * v[930]);
//		v[1889] = v[1814] * v[2626];
//		v[1818] = v[1017] * v[1592] + v[1018] * v[1593] + v[1019] * v[1594] + v[1309] * v[1808] + v[1320] * v[1809]
//			+ v[1335] * v[1810] + v[1569] * v[3134] + v[1543] * v[3135] + v[1548] * v[3136] + v[2785] * v[402] + v[2786] * v[405]
//			+ v[2787] * v[408] - v[1611] * v[490] - v[1615] * v[505] + v[1617] * v[520] + v[1616] * v[535] + (*ct)*(-
//			(v[1624] * v[924]) - v[1625] * v[925] - v[1626] * v[926]);
//		v[2807] = v[1818] * v[227];
//		v[1893] = v[1818] * v[2625];
//		v[1822] = v[1021] * v[1592] + v[1022] * v[1593] + v[1023] * v[1594] + v[1310] * v[1808] + v[1330] * v[1809]
//			+ v[1331] * v[1810] + v[1569] * v[3137] + v[1543] * v[3138] + v[1548] * v[3139] + v[2785] * v[401] + v[2786] * v[404]
//			+ v[2787] * v[407] - v[1611] * v[489] - v[1615] * v[504] + v[1617] * v[519] + v[1616] * v[534] + (*ct)*(-
//			(v[1624] * v[920]) - v[1625] * v[921] - v[1626] * v[922]);
//		v[2808] = v[1822] * v[227];
//		v[1896] = v[1822] * v[2624];
//		v[1823] = v[1544] + v[2770] + v[2772];
//		v[2836] = v[1823] * v[267];
//		v[2815] = -v[2787] + v[2836] + v[2837] + v[1451] * v[990];
//		v[5933] = v[128] * v[2816] - v[1611] * v[477] - v[1615] * v[492] + v[1617] * v[507] + v[1616] * v[522] + v[1593] * v[579] +
//			(*ct)*(-(v[1624] * v[872]) - v[1625] * v[873] - v[1626] * v[874]) + v[1405] * v[965];
//		v[5934] = v[128] * v[2833] - v[1611] * v[478] - v[1615] * v[493] + v[1617] * v[508] + v[1616] * v[523] + v[1592] * v[579] +
//			(*ct)*(-(v[1624] * v[876]) - v[1625] * v[877] - v[1626] * v[878]) + v[1405] * v[964];
//		v[5935] = v[128] * v[2815] - v[1611] * v[479] - v[1615] * v[494] + v[1617] * v[509] + v[1616] * v[524] + (*ct)*(-
//			(v[1624] * v[880]) - v[1625] * v[881] - v[1626] * v[882]);
//		v[5936] = v[130] * v[2816] - v[1611] * v[480] - v[1615] * v[495] + v[1617] * v[510] + v[1616] * v[525] + v[1593] * v[580] +
//			(*ct)*(-(v[1624] * v[884]) - v[1625] * v[885] - v[1626] * v[886]) + v[1404] * v[965];
//		v[5937] = v[130] * v[2833] - v[1611] * v[481] - v[1615] * v[496] + v[1617] * v[511] + v[1616] * v[526] + v[1592] * v[580] +
//			(*ct)*(-(v[1624] * v[888]) - v[1625] * v[889] - v[1626] * v[890]) + v[1404] * v[964];
//		v[5938] = v[130] * v[2815] - v[1611] * v[482] - v[1615] * v[497] + v[1617] * v[512] + v[1616] * v[527] + (*ct)*(-
//			(v[1624] * v[892]) - v[1625] * v[893] - v[1626] * v[894]);
//		v[5939] = v[132] * v[2816] - v[1611] * v[483] - v[1615] * v[498] + v[1617] * v[513] + v[1616] * v[528] + v[1593] * v[581] +
//			(*ct)*(-(v[1624] * v[896]) - v[1625] * v[897] - v[1626] * v[898]) + v[1403] * v[965];
//		v[5940] = v[132] * v[2833] - v[1611] * v[484] - v[1615] * v[499] + v[1617] * v[514] + v[1616] * v[529] + v[1592] * v[581] +
//			(*ct)*(-(v[1624] * v[900]) - v[1625] * v[901] - v[1626] * v[902]) + v[1403] * v[964];
//		v[5941] = v[132] * v[2815] - v[1611] * v[485] - v[1615] * v[500] + v[1617] * v[515] + v[1616] * v[530] + (*ct)*(-
//			(v[1624] * v[904]) - v[1625] * v[905] - v[1626] * v[906]);
//		v[5942] = v[265] * (-v[1559] - v[2771] - v[2772]) + v[2785] + v[1543] * v[2788] - v[2834] - v[1611] * v[486]
//			- v[1615] * v[501] + v[1617] * v[516] + v[1616] * v[531] + (*ct)*(-(v[1624] * v[908]) - v[1625] * v[909]
//				- v[1626] * v[910]);
//		v[5943] = v[266] * (-v[1550] - v[2770] - v[2771]) + v[2786] + v[1548] * v[2789] - v[2835] - v[1611] * v[487]
//			- v[1615] * v[502] + v[1617] * v[517] + v[1616] * v[532] + (*ct)*(-(v[1624] * v[912]) - v[1625] * v[913]
//				- v[1626] * v[914]);
//		v[5944] = v[2787] - v[2836] - v[2837] - v[1611] * v[488] - v[1615] * v[503] + v[1617] * v[518] + v[1616] * v[533] + (*ct)*
//			(-(v[1624] * v[916]) - v[1625] * v[917] - v[1626] * v[918]) - v[1569] * v[990];
//		v[5945] = v[1020] * v[1364] + v[1016] * v[1368] + v[1318] * (v[1899] + v[1016] * v[1967] + v[1020] * v[1969])
//			+ v[1822] * v[1992] + v[1818] * v[2630] + v[1814] * v[2633] + v[1285] * v[2838] + v[1295] * v[2839] + v[1301] * v[2840]
//			+ v[179] * v[2841] + v[184] * v[2842] + v[175] * v[2843];
//		v[5946] = v[1024] * v[1364] + v[1016] * v[1370] + v[1219] * v[1814] + v[1220] * v[1822] + v[1318] * (v[1892]
//			+ v[1016] * v[1966] + v[1024] * v[1970]) + v[1818] * v[1993] + v[1292] * v[2682] + v[1312] * v[2684] + v[1283] * v[2844]
//			+ v[176] * v[2845] + v[181] * v[2846] + v[1378] * v[2954];
//		v[5947] = v[1024] * v[1368] + v[1020] * v[1370] + v[1221] * v[1818] + v[1222] * v[1822] + v[1318] * (v[1884]
//			+ v[1020] * v[1968] + v[1024] * v[1971]) + v[1814] * v[1994] + v[1297] * v[2683] + v[1281] * v[2847] + v[1303] * v[2848]
//			+ v[187] * v[2849] + v[1382] * v[2958] + v[1386] * v[2959];
//		v[1743] = v[1743] - v[1742] * v[407] - v[1741] * v[408] - v[1740] * v[409];
//		v[1824] = v[1569] * v[1831] + v[1569] * v[1839] + v[249] * v[2787] - v[267] * v[2790] - v[1740] * v[587];
//		v[1825] = v[1569] * v[1829] + v[1569] * v[1836] + v[239] * v[2787] - v[267] * v[2791] - v[1741] * v[587];
//		v[1756] = v[1756] - v[1751] * v[401] - v[1750] * v[402] - v[1749] * v[403];
//		v[1414] = v[1414] + v[1740] * v[1813] + v[1741] * v[1817] + v[1742] * v[1821] + v[1411] * v[1823] - v[1544] * v[223]
//			+ 2e0*v[1743] * v[267] + v[2792] * v[407] - v[2791] * v[408] - v[2790] * v[409];
//		v[1748] = v[1748] - v[1747] * v[404] - v[1746] * v[405] - v[1745] * v[406];
//		v[1826] = v[1569] * v[1827] + v[1569] * v[1833] + v[229] * v[2787] + v[267] * v[2792] - v[1742] * v[587];
//		v[1828] = v[1548] * v[1827] + v[1548] * v[1834] + v[229] * v[2786] - v[266] * (v[1552] - v[2793] + v[2796])
//			- v[1747] * v[582];
//		v[1830] = v[1548] * v[1829] + v[1548] * v[1837] + v[239] * v[2786] - v[266] * (v[1554] + v[2794] + v[2798])
//			- v[1746] * v[582];
//		v[1406] = v[1406] - v[1550] * v[222] + v[1745] * v[2644] + v[1746] * v[2646] + v[1747] * v[2648] + 2e0*v[1748] * v[266]
//			+ v[1744] * v[2894] - v[1552] * v[404] - v[1554] * v[405] - v[1556] * v[406] + v[267] * (-(v[1742] * v[404])
//				- v[1741] * v[405] - v[1740] * v[406]) + v[265] * (v[2800] - v[1751] * v[404] - v[1750] * v[405] - v[1749] * v[406]);
//		v[1832] = v[1548] * v[1831] + v[1548] * v[1840] + v[249] * v[2786] - v[266] * (v[1556] + v[2795] + v[2801])
//			- v[1745] * v[582];
//		v[1835] = v[1543] * v[1833] + v[1543] * v[1834] + v[229] * v[2785] - v[265] * (v[1562] + v[2796] - v[2797])
//			- v[1751] * v[577];
//		v[1838] = v[1543] * v[1836] + v[1543] * v[1837] + v[239] * v[2785] - v[265] * (v[1565] + v[2798] + v[2799])
//			- v[1750] * v[577];
//		v[1408] = v[1408] - v[1559] * v[221] - v[1749] * v[2645] - v[1750] * v[2647] - v[1751] * v[2649] + 2e0*v[1756] * v[265]
//			+ v[1744] * v[2905] - v[1562] * v[401] - v[1565] * v[402] - v[1568] * v[403] + v[267] * (-(v[1742] * v[401])
//				- v[1741] * v[402] - v[1740] * v[403]) + v[266] * (v[2800] - v[1747] * v[401] - v[1746] * v[402] - v[1745] * v[403]);
//		v[1841] = v[1543] * v[1839] + v[1543] * v[1840] + v[249] * v[2785] - v[265] * (v[1568] + v[2801] + v[2802])
//			- v[1749] * v[577];
//		v[1396] = v[1396] + v[1573] * v[1842] * v[262] + v[263] * (v[1408] * v[250] + v[1406] * v[251] + v[1414] * v[252]
//			+ v[1389] * v[983] + v[1391] * v[984] + v[1393] * v[985]);
//		v[1844] = v[1414] * v[264] + (v[1396] * v[252] + v[1393] * v[967]) / v[597] + v[2722] * v[985];
//		v[1846] = v[1406] * v[264] + (v[1396] * v[251] + v[1391] * v[967]) / v[597] + v[2722] * v[984];
//		v[1848] = v[1408] * v[264] + (v[1396] * v[250] + v[1389] * v[967]) / v[597] + v[2722] * v[983];
//		v[1781] = v[1781] - v[1844];
//		v[1774] = v[1774] - v[1846];
//		v[1767] = v[1767] - v[1848];
//		v[1852] = i2720 * v[1115] + i2719 * v[1116] + i2718 * v[1117] + v[152] * v[1844] + v[148] * v[1846] + v[144] * v[1848]
//			+ v[214] * v[2815] + v[212] * v[2816] + v[213] * v[2817] + (*a4)*v[5842 + i933] + v[2803] * v[997] + v[266] *
//			(v[1744] * v[213] + v[1753] * v[265] + v[1399] * v[997]);
//		v[1853] = v[1814] * v[236] + v[1016] * v[2804];
//		v[1854] = v[1814] * v[230] + v[1016] * v[2805];
//		v[1855] = v[1127] * v[1814] + v[1016] * (v[1349] / v[227] + v[1318] * v[2448]) + v[1853] * v[2916];
//		v[1856] = v[1133] * v[1814] + v[1016] * (v[1341] / v[227] + v[1318] * v[2450]) + v[1854] * v[2917];
//		v[1857] = (v[184] * v[1854] + v[1853] * v[186] + v[1016] * (v[1357] + v[2452] * v[2717]) + v[1814] * v[2918]) / v[227];
//		v[1858] = v[1818] * v[238] + v[1020] * v[2806];
//		v[1859] = v[1818] * v[230] + v[1020] * v[2805];
//		v[1862] = v[1122] * v[1818] + v[1020] * (v[1355] / v[227] + v[1318] * v[2456]) + v[1858] * v[2920];
//		v[1863] = v[1853] + v[1858];
//		v[1864] = v[1855] + v[1862];
//		v[1866] = (v[1125] * v[1303] + v[175] * v[1859] + v[177] * v[1863] + v[1962] * v[2717] + v[1020] * (v[1337]
//			+ v[2460] * v[2717]) + v[1132] * v[2807]) / v[227];
//		v[1869] = (v[182] * v[1858] + v[179] * v[1859] + v[1020] * (v[1346] + v[2462] * v[2717]) + v[1126] * v[2807]) / v[227];
//		v[1870] = v[1822] * v[236] + v[1024] * v[2804];
//		v[1871] = v[1822] * v[238] + v[1024] * v[2806];
//		v[1872] = v[1859] + v[1870];
//		v[1875] = (v[176] * v[1870] + v[177] * v[1871] + v[1024] * (v[1338] + v[2467] * v[2717]) + v[1130] * v[2808]) / v[227];
//		v[1876] = v[1854] + v[1871];
//		v[1878] = (v[1135] * v[1297] + v[181] * v[1870] + v[182] * v[1876] + v[1961] * v[2717] + v[1024] * (v[1345]
//			+ v[2470] * v[2717]) + v[1137] * v[2808]) / v[227];
//		v[1879] = v[1866] + v[1878];
//		v[1881] = (v[1136] * v[1292] + v[187] * v[1871] + v[186] * v[1872] + v[1960] * v[2717] + v[1024] * (v[1353]
//			+ v[2473] * v[2717]) + v[1140] * v[2808]) / v[227];
//		v[1882] = v[1856] + v[1881];
//		v[1782] = v[1782] - v[168] * v[1844] + v[1804] * v[290] + v[1805] * v[300];
//		v[1783] = v[1783] - v[167] * v[1844] + v[1804] * v[289] + v[1805] * v[298];
//		v[1784] = v[1784] - v[166] * v[1844] + v[1804] * v[288] + v[1805] * v[296];
//		v[1775] = v[1775] - v[168] * v[1846] + v[1800] * v[290] + v[1801] * v[300];
//		v[1776] = v[1776] - v[167] * v[1846] + v[1800] * v[289] + v[1801] * v[298];
//		v[1777] = v[1777] - v[166] * v[1846] + v[1800] * v[288] + v[1801] * v[296];
//		v[1768] = v[1768] - v[168] * v[1848] + v[1796] * v[290] + v[1797] * v[300];
//		v[1769] = v[1769] - v[167] * v[1848] + v[1796] * v[289] + v[1797] * v[298];
//		v[1770] = v[1770] - v[166] * v[1848] + v[1796] * v[288] + v[1797] * v[296];
//		v[1883] = QBi[2][2] * v[1782] + QBi[2][1] * v[1783] + QBi[2][0] * v[1784] + v[1219] * v[1858] + v[1889] * v[240]
//			+ v[1871] * v[2633] + v[238] * (v[2849] + v[1318] * v[3141]);
//		v[1886] = QBi[1][2] * v[1782] + QBi[1][1] * v[1783] + QBi[1][0] * v[1784] + v[1853] * v[1994] + (v[1136] * v[1361])
//			/ v[227] + v[1893] * v[237] + v[1378] * v[2625] + v[1872] * v[2633] + v[1318] * v[3142];
//		v[1888] = QBi[0][2] * v[1782] + QBi[0][1] * v[1783] + QBi[0][0] * v[1784] + v[1854] * v[1994] + v[1896] * v[225]
//			+ v[230] * (v[2842] + v[1318] * v[3143]);
//		v[1890] = QBi[2][2] * v[1775] + QBi[2][1] * v[1776] + QBi[2][0] * v[1777] + v[1858] * v[1993] + (v[1135] * v[1365])
//			/ v[227] + v[1386] * v[2626] + v[1876] * v[2630] + v[1889] * v[2634] + v[1318] * v[3144];
//		v[1894] = QBi[1][2] * v[1775] + QBi[1][1] * v[1776] + QBi[1][0] * v[1777] + v[1221] * v[1853] + v[1893] * v[231]
//			+ v[1870] * v[2630] + v[236] * (v[2846] + v[1318] * v[3145]);
//		v[1897] = QBi[0][2] * v[1775] + QBi[0][1] * v[1776] + QBi[0][0] * v[1777] + v[1859] * v[1993] + v[1896] * v[226]
//			+ v[230] * (v[2841] + v[1318] * v[3146]);
//		v[1898] = QBi[2][2] * v[1768] + QBi[2][1] * v[1769] + QBi[2][0] * v[1770] + v[1220] * v[1863] + v[1871] * v[1992] +
//			(v[1125] * v[1366]) / v[227] + v[1382] * v[2626] + v[1889] * v[2632] + v[1318] * v[3147];
//		v[1900] = QBi[1][2] * v[1768] + QBi[1][1] * v[1769] + QBi[1][0] * v[1770] + v[1870] * v[1992] + v[1893] * v[2629]
//			+ v[236] * (v[1318] * v[2811] + v[2845]);
//		v[1903] = QBi[0][2] * v[1768] + QBi[0][1] * v[1769] + QBi[0][0] * v[1770] + v[1222] * v[1854] + v[1220] * v[1859]
//			+ v[1896] * v[228] + v[230] * (v[2843] + v[1318] * v[3148]);
//		v[1904] = -(v[1838] * v[722]);
//		v[1905] = -(v[1838] * v[720]);
//		v[1906] = -(v[1838] * v[718]);
//		v[1907] = -(v[1841] * v[722]);
//		v[1908] = -(v[1841] * v[718]);
//		v[1909] = -(v[1841] * v[720]);
//		v[1910] = -(v[1835] * v[720]);
//		v[1911] = -(v[1835] * v[718]);
//		v[1912] = -(v[1835] * v[722]);
//		v[1913] = -(v[1828] * v[722]);
//		v[1914] = -(v[1828] * v[720]);
//		v[1915] = -(v[1828] * v[718]);
//		v[1916] = -(v[1832] * v[718]);
//		v[1917] = -(v[1832] * v[722]);
//		v[1918] = -(v[1832] * v[720]);
//		v[1919] = -(v[1824] * v[720]);
//		v[1920] = -(v[1824] * v[722]);
//		v[1921] = -(v[1824] * v[718]);
//		v[1922] = -(v[1129] * v[1277]) - v[1139] * v[1279] + 2e0*v[1128] * v[1280] + v[1910] + v[1913] + v[1916] + v[1919]
//			+ 2e0*v[1869] * v[318] - v[1879] * v[321] - v[1864] * v[324];
//		v[1923] = -(v[1830] * v[722]);
//		v[1924] = -(v[1830] * v[718]);
//		v[1925] = -(v[1830] * v[720]);
//		v[1926] = v[1142] * v[1277] + 2e0*v[1134] * v[1279] - v[1139] * v[1280] - v[1905] - v[1908] - v[1920] - v[1923]
//			- v[1879] * v[318] + 2e0*v[1875] * v[321] + v[1882] * v[324];
//		v[1927] = -(v[1166] * v[1268]) + v[1162] * v[1271] - v[1167] * v[1273] - 2e0*v[1157] * v[1278] + v[172] * v[1900]
//			- v[1910] * v[311] + v[1905] * v[312] - v[1909] * v[313];
//		v[1928] = -(v[1826] * v[722]);
//		v[1929] = -(v[1826] * v[720]);
//		v[1930] = -(v[1826] * v[718]);
//		v[1931] = -(v[1825] * v[720]);
//		v[1932] = -(v[1825] * v[722]);
//		v[1933] = -(v[1825] * v[718]);
//		v[1934] = 2e0*v[1121] * v[1277] + v[1142] * v[1279] - v[1129] * v[1280] - v[1911] - v[1924] - v[1928] - v[1931]
//			- v[1864] * v[318] + v[1882] * v[321] + 2e0*v[1857] * v[324];
//		v[1935] = -(v[1165] * v[1268]) + v[1163] * v[1271] - v[1168] * v[1273] - 2e0*v[1155] * v[1278] + v[172] * v[1898]
//			- v[1911] * v[311] + v[1906] * v[312] - v[1908] * v[313];
//		v[1936] = -(v[1174] * v[1268]) + v[1180] * v[1271] - v[1170] * v[1273] - 2e0*v[1154] * v[1278] + v[172] * v[1897]
//			- v[1913] * v[311] + v[1923] * v[312] - v[1917] * v[313];
//		v[1937] = v[1907] + v[1918];
//		v[1938] = -(v[1173] * v[1268]) + v[1181] * v[1271] - v[1172] * v[1273] - 2e0*v[1149] * v[1278] + v[172] * v[1890]
//			- v[1915] * v[311] + v[1924] * v[312] - v[1916] * v[313];
//		v[1939] = -(v[1177] * v[1268]) + v[1189] * v[1271] - v[1185] * v[1273] - 2e0*v[1148] * v[1278] + v[172] * v[1888]
//			- v[1928] * v[311] + v[1932] * v[312] - v[1920] * v[313];
//		v[1940] = -(v[1176] * v[1268]) + v[1188] * v[1271] - v[1186] * v[1273] - 2e0*v[1146] * v[1278] + v[172] * v[1886]
//			- v[1929] * v[311] + v[1931] * v[312] - v[1919] * v[313];
//		v[1941] = v[1914] + v[1930];
//		v[1942] = v[1904] + v[1933];
//		v[1943] = -(QBi[2][2] * v[1333]) - QBi[1][2] * v[1334] - QBi[0][2] * v[1336] - v[1117] * v[1344] - v[1116] * v[1352]
//			- v[1115] * v[1360] - v[1848] * v[193] - v[1846] * v[196] - v[1844] * v[199] + v[1835] * v[362] + v[1838] * v[363]
//			+ v[1841] * v[364] + v[1828] * v[377] + v[1830] * v[378] + v[1832] * v[379] + v[1826] * v[392] + v[1825] * v[393]
//			+ v[1824] * v[394];
//		v[1944] = -(QBi[2][1] * v[1333]) - QBi[1][1] * v[1334] - QBi[0][1] * v[1336] - v[1117] * v[1343] - v[1116] * v[1351]
//			- v[1115] * v[1359] - v[1848] * v[192] - v[1846] * v[195] - v[1844] * v[198] + v[1835] * v[359] + v[1838] * v[360]
//			+ v[1841] * v[361] + v[1828] * v[374] + v[1830] * v[375] + v[1832] * v[376] + v[1826] * v[389] + v[1825] * v[390]
//			+ v[1824] * v[391];
//		v[1945] = -(QBi[2][0] * v[1333]) - QBi[1][0] * v[1334] - QBi[0][0] * v[1336] - v[1117] * v[1342] - v[1116] * v[1350]
//			- v[1115] * v[1358] - v[1848] * v[191] - v[1846] * v[194] - v[1844] * v[197] + v[1835] * v[356] + v[1838] * v[357]
//			+ v[1841] * v[358] + v[1828] * v[371] + v[1830] * v[372] + v[1832] * v[373] + v[1826] * v[386] + v[1825] * v[387]
//			+ v[1824] * v[388];
//		v[1946] = v[1945] * x1B[0] + v[1944] * x1B[1] + v[1943] * x1B[2];
//		v[1947] = (-v[1946] + v[1945] * x3B[0] + v[1944] * x3B[1] + v[1943] * x3B[2]) / 2e0;
//		v[1948] = (-v[1946] + v[1945] * x2B[0] + v[1944] * x2B[1] + v[1943] * x2B[2]) / 2e0;
//		v[1949] = (-(v[1187] * v[1233]) - v[1169] * v[1236] - 2e0*v[1170] * v[1261] - 2e0*v[1167] * v[1262]
//			- 2e0*v[1186] * v[1263] - 2e0*v[1172] * v[1264] - 2e0*v[1185] * v[1265] - 2e0*v[1168] * v[1266] - v[1171] * v[1269]
//			- 2e0*v[1855] + 2e0*v[1862] - v[1912] * v[314] - 2e0*v[1910] * v[319] - 2e0*v[1911] * v[325] - 2e0*v[1913] * v[329]
//			- v[1914] * v[333] - 2e0*v[1915] * v[338] - 2e0*v[1928] * v[342] - 2e0*v[1929] * v[346] - v[1930] * v[351]) / 2e0;
//		v[1950] = (-(v[1164] * v[1268]) + v[1161] * v[1271] - v[1169] * v[1273] - 2e0*v[1160] * v[1278] + v[172] * v[1903]
//			- v[1912] * v[311] + v[1904] * v[312] - v[1907] * v[313]) / 2e0;
//		v[1952] = (v[1190] * v[1233] + v[1161] * v[1236] + 2e0*v[1180] * v[1261] + 2e0*v[1162] * v[1262]
//			+ 2e0*v[1188] * v[1263] + 2e0*v[1181] * v[1264] + 2e0*v[1189] * v[1265] + 2e0*v[1163] * v[1266] + v[1182] * v[1269]
//			- 2e0*v[1856] + 2e0*v[1881] + v[1904] * v[314] + 2e0*v[1905] * v[319] + 2e0*v[1906] * v[325] + 2e0*v[1923] * v[329]
//			+ v[1925] * v[333] + 2e0*v[1924] * v[338] + 2e0*v[1932] * v[342] + 2e0*v[1931] * v[346] + v[1933] * v[351]) / 2e0;
//		v[1953] = (-(v[1175] * v[1268]) + v[1182] * v[1271] - v[1171] * v[1273] - 2e0*v[1152] * v[1278] + v[172] * v[1894]
//			- v[1914] * v[311] + v[1925] * v[312] - v[1918] * v[313]) / 2e0;
//		v[1954] = (-(v[1178] * v[1233]) - v[1164] * v[1236] - 2e0*v[1174] * v[1261] - 2e0*v[1166] * v[1262]
//			- 2e0*v[1176] * v[1263] - 2e0*v[1173] * v[1264] - 2e0*v[1177] * v[1265] - 2e0*v[1165] * v[1266] - v[1175] * v[1269]
//			- 2e0*v[1866] + 2e0*v[1878] - v[1907] * v[314] - 2e0*v[1909] * v[319] - 2e0*v[1908] * v[325] - 2e0*v[1917] * v[329]
//			- v[1918] * v[333] - 2e0*v[1916] * v[338] - 2e0*v[1920] * v[342] - 2e0*v[1919] * v[346] - v[1921] * v[351]) / 2e0;
//		v[1985] = v[1270] * v[1952] - v[1953] + v[1267] * v[1954] + 24e0*v[1235] * v[2576] * v[3150] + v[2827] * (-
//			(v[1143] * v[1233]) - v[1160] * v[1236] - 2e0*v[1154] * v[1261] - 2e0*v[1157] * v[1262] - 2e0*v[1146] * v[1263]
//			- 2e0*v[1149] * v[1264] - 2e0*v[1148] * v[1265] - 2e0*v[1155] * v[1266] - v[1152] * v[1269] - 2e0*v[1906]
//			+ 2e0*v[1909] + 2e0*v[1915] - 2e0*v[1917] - alphaB[1] * v[1922] + alphaB[0] * v[1926] - 2e0*v[1929] + 2e0*v[1932]
//			+ alphaB[2] * v[1934] + 8e0*v[1278] * v[1956] + v[1941] * v[305] + v[1942] * v[308] + v[1937] * v[310] - v[1903] * v[314]
//			- 4e0*v[172] * (v[1857] + v[1378] * v[1860] + v[1375] * v[1861] + v[1369] * v[1863] + v[1384] * v[1865]
//				+ v[1388] * v[1867] + v[1386] * v[1868] + v[1869] + v[1363] * v[1872] + v[1371] * v[1873] + v[1382] * v[1874] + v[1875]
//				+ v[1367] * v[1876] + v[1374] * v[1877] + v[1380] * v[1880] + v[1357] * v[1884] + v[1355] * v[1885] + v[1353] * v[1887]
//				+ v[1349] * v[1891] + v[1346] * v[1892] + v[1345] * v[1895] + v[1338] * v[1899] + v[1341] * v[1901] + v[1337] * v[1902]
//				+ v[1292] * v[1957] + v[1297] * v[1958] + v[1303] * v[1959] + v[1361] * v[1960] + v[1365] * v[1961] + v[1366] * v[1962]
//				+ v[1853] * v[1966] + v[1854] * v[1967] + v[1858] * v[1968] + v[1859] * v[1969] + v[1870] * v[1970] + v[1871] * v[1971]
//				+ v[1814] * v[2812] + v[1818] * v[2813] + v[1822] * v[2814] + 2e0*v[1318] * v[2578] * v[3152]) - 2e0*v[1900] * v[319]
//			- 2e0*v[1898] * v[325] - 2e0*v[1897] * v[329] - v[1894] * v[333] - 2e0*v[1890] * v[338] - 2e0*v[1888] * v[342]
//			- 2e0*v[1886] * v[346] - v[1883] * v[351] - 2e0*v[5797 + i933]) - 2e0*v[1215] * (-4e0*v[1204] * v[1235]
//				+ 4e0*v[1949] * v[305] + v[5857 + i933]);
//		v[2832] = v[1985] + (v[1178] * v[1268] - v[1190] * v[1271] + v[1187] * v[1273] + 2e0*v[1143] * v[1278]
//			- v[172] * v[1883] + v[1930] * v[311] - v[1933] * v[312] + v[1921] * v[313]) / 2e0;
//		v[2831] = (v[1935] + v[1939]) / 2e0;
//		v[1974] = v[1938] + v[1940];
//		v[1975] = v[1927] + v[1936];
//		v[1980] = (i2730*v[1115] + i2733 * v[1116] + i2736 * v[1117] + v[153] * v[1844] + v[149] * v[1846] + v[145] * v[1848]
//			- v[1852] + v[217] * v[2815] + v[215] * v[2816] + v[216] * v[2817] + (*a4)*v[5827 + i933] + v[2803] * v[998] + v[266] *
//			(v[1744] * v[216] + v[1754] * v[265] + v[1399] * v[998])) / 2e0;
//		v[1981] = (i2731*v[1115] + i2734 * v[1116] + i2737 * v[1117] + v[154] * v[1844] + v[150] * v[1846] + v[146] * v[1848]
//			- v[1852] + v[220] * v[2815] + v[218] * v[2816] + v[219] * v[2817] + (*a4)*v[5812 + i933] + v[2803] * v[999] + v[266] *
//			(v[1744] * v[219] + v[1755] * v[265] + v[1399] * v[999])) / 2e0;
//		v[1785] = v[1785] + v[132] * v[1844] + v[2818];
//		v[1778] = v[1778] + v[132] * v[1846] + v[2820];
//		v[1771] = v[1771] + v[132] * v[1848] + v[2822];
//		v[1786] = v[1786] + v[130] * v[1844] + v[2819];
//		v[1779] = v[1779] + v[130] * v[1846] + v[2821];
//		v[1772] = v[1772] + v[130] * v[1848] + v[2823];
//		v[1787] = v[1787] + v[128] * v[1844] - v[2818] - v[2819];
//		v[1780] = v[1780] + v[128] * v[1846] - v[2820] - v[2821];
//		v[1773] = v[1773] + v[128] * v[1848] - v[2822] - v[2823];
//		v[5873] = 0e0;
//		v[5874] = 0e0;
//		v[5875] = 0e0;
//		v[5876] = 0e0;
//		v[5877] = 0e0;
//		v[5878] = 0e0;
//		v[5879] = 0e0;
//		v[5880] = 0e0;
//		v[5881] = 0e0;
//		v[5882] = 0e0;
//		v[5883] = 0e0;
//		v[5884] = 0e0;
//		v[5885] = 2e0*v[1982];
//		v[5886] = v[1983];
//		v[5887] = v[1984];
//		v[5888] = 0e0;
//		v[5889] = 0e0;
//		v[5890] = 0e0;
//		v[5891] = 0e0;
//		v[5892] = 0e0;
//		v[5893] = 0e0;
//		v[5894] = 0e0;
//		v[5895] = 0e0;
//		v[5896] = 0e0;
//		v[5897] = 0e0;
//		v[5898] = 0e0;
//		v[5899] = 0e0;
//		v[5900] = v[1983];
//		v[5901] = 2e0*v[1986];
//		v[5902] = v[1987];
//		v[5903] = 0e0;
//		v[5904] = 0e0;
//		v[5905] = 0e0;
//		v[5906] = 0e0;
//		v[5907] = 0e0;
//		v[5908] = 0e0;
//		v[5909] = 0e0;
//		v[5910] = 0e0;
//		v[5911] = 0e0;
//		v[5912] = 0e0;
//		v[5913] = 0e0;
//		v[5914] = 0e0;
//		v[5915] = v[1984];
//		v[5916] = v[1987];
//		v[5917] = 2e0*v[1988];
//		v[5918] = v[1773];
//		v[5919] = v[1780];
//		v[5920] = v[1787];
//		v[5921] = v[1772];
//		v[5922] = v[1779];
//		v[5923] = v[1786];
//		v[5924] = v[1771];
//		v[5925] = v[1778];
//		v[5926] = v[1785];
//		v[5927] = v[1767];
//		v[5928] = v[1774];
//		v[5929] = v[1781];
//		v[5930] = v[1183] * v[1278] - v[1938] + v[1940] + v[1926] * v[2824] + 2e0*(v[1941] * v[2824] + v[1949] * v[2826]
//			+ v[1235] * (v[1201] * v[2827] - v[1215] * v[2830])) + alphaB[2] * v[2831] + 2e0*alphaB[0] * v[2832] + v[1975] * v[306]
//			+ v[5872 + i933];
//		v[5931] = -(v[1179] * v[1278]) + v[1935] - v[1939] + (alphaB[2] * v[1974]) / 2e0 + (alphaB[0] * v[1975]) / 2e0
//			- v[1922] * v[2824] + 2e0*alphaB[1] * (-v[1950] + v[1953] + v[2832]) + v[5887 + i933] + 2e0*(v[1942] * v[2824]
//				+ v[1235] * (v[1202] * v[2827] + v[1215] * v[2829]) - 4e0*v[1952] * v[783]);
//		v[5932] = v[1194] * v[1278] - v[1927] + v[1936] + 2e0*alphaB[2] * (-v[1950] + v[1985]) + v[1934] * v[2824] + 2e0*
//			(v[1937] * v[2824] + v[1954] * v[2826] + v[1235] * (v[1197] * v[2827] - v[1215] * v[2828])) + alphaB[0] * v[2831]
//			+ v[1974] * v[306] + v[5902 + i933];
//		Rc[i933 - 1] += v[4533 + i933] + (*a4)*v[4548 + i933];
//		for (i1225 = 1; i1225 <= 15; i1225++) {
//			Kc[i933 - 1][i1225 - 1] += v[1981] * v[4105 + i1225] + v[1980] * v[4120 + i1225] + v[1948] * v[4135 + i1225]
//				+ v[1947] * v[4150 + i1225] + v[5917 + i1225] + (*a4)*v[5932 + i1225];
//		};/* end for */
//	};/* end for */
//	v[2003] = 0e0;
//	v[2004] = 0e0;
//	v[2005] = 0e0;
//	v[2006] = 0e0;
//	b2007 = b598;
//	if (b2007) {
//		b2008 = b617;
//		b2012 = b600;
//		if (b2012) {
//			v[2005] = 0e0;
//			v[2004] = 0e0;
//			v[2003] = 0e0;
//			v[2006] = 0e0;
//		}
//		else {
//		};
//	}
//	else {
//	};
//	v[6076] = v[128] * v[2003];
//	v[6077] = 0e0;
//	v[6078] = 0e0;
//	v[6079] = v[130] * v[2003];
//	v[6080] = 0e0;
//	v[6081] = 0e0;
//	v[6082] = v[132] * v[2003];
//	v[6083] = 0e0;
//	v[6084] = 0e0;
//	v[6085] = -v[2003];
//	v[6086] = 0e0;
//	v[6087] = 0e0;
//	v[6088] = 0e0;
//	v[6089] = 0e0;
//	v[6090] = 0e0;
//	v[2621] = v[2003] * v[265];
//	v[2438] = 2e0*v[2621];
//	v[2436] = -(v[2003] * v[577]);
//	v[2408] = -(v[249] * v[2621]);
//	v[2403] = -(v[239] * v[2621]);
//	v[2400] = -(v[229] * v[2621]);
//	v[6061] = 0e0;
//	v[6062] = v[128] * v[2004];
//	v[6063] = 0e0;
//	v[6064] = 0e0;
//	v[6065] = v[130] * v[2004];
//	v[6066] = 0e0;
//	v[6067] = 0e0;
//	v[6068] = v[132] * v[2004];
//	v[6069] = 0e0;
//	v[6070] = 0e0;
//	v[6071] = -v[2004];
//	v[6072] = 0e0;
//	v[6073] = 0e0;
//	v[6074] = 0e0;
//	v[6075] = 0e0;
//	v[6046] = 0e0;
//	v[6047] = 0e0;
//	v[6048] = 0e0;
//	v[6049] = 0e0;
//	v[6050] = 0e0;
//	v[6051] = 0e0;
//	v[6052] = v[2004];
//	v[6053] = v[2003];
//	v[6054] = 0e0;
//	v[6055] = 0e0;
//	v[6056] = 0e0;
//	v[6057] = 0e0;
//	v[6058] = 0e0;
//	v[6059] = 0e0;
//	v[6060] = 0e0;
//	v[6031] = 0e0;
//	v[6032] = 0e0;
//	v[6033] = 0e0;
//	v[6034] = v[2004];
//	v[6035] = v[2003];
//	v[6036] = 0e0;
//	v[6037] = 0e0;
//	v[6038] = 0e0;
//	v[6039] = 0e0;
//	v[6040] = 0e0;
//	v[6041] = 0e0;
//	v[6042] = 0e0;
//	v[6043] = 0e0;
//	v[6044] = 0e0;
//	v[6045] = 0e0;
//	v[6016] = v[2004];
//	v[6017] = v[2003];
//	v[6018] = 0e0;
//	v[6019] = 0e0;
//	v[6020] = 0e0;
//	v[6021] = 0e0;
//	v[6022] = 0e0;
//	v[6023] = 0e0;
//	v[6024] = 0e0;
//	v[6025] = 0e0;
//	v[6026] = 0e0;
//	v[6027] = 0e0;
//	v[6028] = 0e0;
//	v[6029] = 0e0;
//	v[6030] = 0e0;
//	v[2619] = v[2004] * v[266];
//	v[2441] = 2e0*v[2619];
//	v[2437] = -(v[2004] * v[582]);
//	v[2409] = -(v[249] * v[2619]);
//	v[2404] = -(v[239] * v[2619]);
//	v[2401] = -(v[229] * v[2619]);
//	v[5986] = 0e0;
//	v[5987] = 0e0;
//	v[5988] = v[128] * v[2005];
//	v[5989] = 0e0;
//	v[5990] = 0e0;
//	v[5991] = v[130] * v[2005];
//	v[5992] = 0e0;
//	v[5993] = 0e0;
//	v[5994] = v[132] * v[2005];
//	v[5995] = 0e0;
//	v[5996] = 0e0;
//	v[5997] = -v[2005];
//	v[5998] = 0e0;
//	v[5999] = 0e0;
//	v[6000] = 0e0;
//	v[2850] = v[2005] * v[267];
//	v[2887] = -v[2619] - v[2850];
//	v[3204] = v[2003] * (2e0*v[1001] - v[2649]) + v[2887] * v[401];
//	v[3201] = v[2003] * (2e0*v[1002] - v[2647]) + v[2887] * v[402];
//	v[3198] = v[2003] * (2e0*v[1003] - v[2645]) + v[2887] * v[403];
//	v[2886] = -v[2621] - v[2850];
//	v[3203] = v[2004] * (2e0*v[1007] + v[2648]) + v[2886] * v[404];
//	v[3200] = v[2004] * (2e0*v[1009] + v[2646]) + v[2886] * v[405];
//	v[3197] = v[2004] * (2e0*v[1011] + v[2644]) + v[2886] * v[406];
//	v[2393] = 2e0*v[2850];
//	v[2617] = (v[132] * v[2393]) / 2e0;
//	v[2615] = (v[130] * v[2393]) / 2e0;
//	v[2613] = (v[128] * v[2393]) / 2e0;
//	v[2417] = -(v[2393] * v[249]) / 2e0;
//	v[3216] = 2e0*v[2408] + v[2409] + v[2417];
//	v[3212] = v[2408] + 2e0*v[2409] + v[2417];
//	v[3205] = v[2408] + v[2409] + 2e0*v[2417];
//	v[2413] = -(v[239] * v[2393]) / 2e0;
//	v[3214] = 2e0*v[2403] + v[2404] + v[2413];
//	v[3210] = v[2403] + 2e0*v[2404] + v[2413];
//	v[3206] = v[2403] + v[2404] + 2e0*v[2413];
//	v[2411] = -(v[229] * v[2393]) / 2e0;
//	v[3213] = 2e0*v[2400] + v[2401] + v[2411];
//	v[3209] = v[2400] + 2e0*v[2401] + v[2411];
//	v[3208] = v[2400] + v[2401] + 2e0*v[2411];
//	v[2386] = -(v[2619] * v[267]) - v[2621] * v[267] - v[2005] * v[587];
//	v[2033] = -(v[266] * v[2850]);
//	v[2593] = -v[2033] - v[2437];
//	v[2387] = v[2033] + v[2437] - v[2621] * v[266];
//	v[2028] = -(v[265] * v[2850]);
//	v[2592] = -v[2028] - v[2436];
//	v[6196] = 0e0;
//	v[6197] = 0e0;
//	v[6198] = 0e0;
//	v[6199] = 0e0;
//	v[6200] = 0e0;
//	v[6201] = 0e0;
//	v[6202] = v[2592];
//	v[6203] = v[2593];
//	v[6204] = -v[2386];
//	v[6205] = 0e0;
//	v[6206] = 0e0;
//	v[6207] = 0e0;
//	v[6208] = 0e0;
//	v[6209] = 0e0;
//	v[6210] = 0e0;
//	v[6226] = 0e0;
//	v[6227] = 0e0;
//	v[6228] = 0e0;
//	v[6229] = v[2592];
//	v[6230] = v[2593];
//	v[6231] = -v[2386];
//	v[6232] = 0e0;
//	v[6233] = 0e0;
//	v[6234] = 0e0;
//	v[6235] = 0e0;
//	v[6236] = 0e0;
//	v[6237] = 0e0;
//	v[6238] = 0e0;
//	v[6239] = 0e0;
//	v[6240] = 0e0;
//	v[6256] = v[2592];
//	v[6257] = v[2593];
//	v[6258] = -v[2386];
//	v[6259] = 0e0;
//	v[6260] = 0e0;
//	v[6261] = 0e0;
//	v[6262] = 0e0;
//	v[6263] = 0e0;
//	v[6264] = 0e0;
//	v[6265] = 0e0;
//	v[6266] = 0e0;
//	v[6267] = 0e0;
//	v[6268] = 0e0;
//	v[6269] = 0e0;
//	v[6270] = 0e0;
//	v[2388] = v[2028] + v[2436] - v[2619] * v[265];
//	v[2013] = 0e0;
//	v[2014] = 0e0;
//	v[2015] = 0e0;
//	v[2016] = 0e0;
//	b2017 = b598;
//	if (b2017) {
//		b2018 = b600;
//		if (b2018) {
//			v[2015] = 0e0;
//			v[2014] = 0e0;
//			v[2013] = 0e0;
//			v[2016] = -v[2006];
//		}
//		else {
//		};
//	}
//	else {
//	};
//	v[2019] = v[2619] + v[2621];
//	v[3202] = v[2005] * (v[1821] + 2e0*v[2641]) - v[2019] * v[407];
//	v[3199] = v[2005] * (v[1817] + 2e0*v[2642]) - v[2019] * v[408];
//	v[3196] = v[2005] * (v[1813] + 2e0*v[2643]) - v[2019] * v[409];
//	v[2618] = v[132] * v[2019];
//	v[2616] = v[130] * v[2019];
//	v[2614] = v[128] * v[2019];
//	v[6001] = 0e0;
//	v[6002] = 0e0;
//	v[6003] = v[2614];
//	v[6004] = 0e0;
//	v[6005] = 0e0;
//	v[6006] = v[2616];
//	v[6007] = 0e0;
//	v[6008] = 0e0;
//	v[6009] = v[2618];
//	v[6010] = 0e0;
//	v[6011] = 0e0;
//	v[6012] = 0e0;
//	v[6013] = 0e0;
//	v[6014] = 0e0;
//	v[6015] = 0e0;
//	v[2604] = v[2019] * v[220];
//	v[3241] = v[220] * v[2393] + v[2604];
//	v[3164] = v[2604] + v[2005] * v[2881];
//	v[2598] = v[2019] * v[217];
//	v[3236] = v[217] * v[2393] + v[2598];
//	v[3165] = v[2598] + v[2005] * v[2880];
//	v[2444] = v[2019] * v[214];
//	v[3219] = v[214] * v[2393] + v[2444];
//	v[3163] = v[2444] + v[2005] * v[2879];
//	v[2020] = v[2004] * v[212] + v[2003] * v[213];
//	v[2913] = v[2020] * v[265] + v[213] * v[2850];
//	v[3221] = v[213] * v[2441] + v[2913];
//	v[2440] = v[2020] * v[266];
//	v[2912] = v[2440] + v[212] * v[2850];
//	v[3220] = v[212] * v[2438] + v[2912];
//	v[2021] = v[2004] * v[215] + v[2003] * v[216];
//	v[2930] = v[2021] * v[265] + v[216] * v[2850];
//	v[3238] = v[216] * v[2441] + v[2930];
//	v[2595] = v[2021] * v[266];
//	v[2929] = v[2595] + v[215] * v[2850];
//	v[3237] = v[215] * v[2438] + v[2929];
//	v[2022] = v[2004] * v[218] + v[2003] * v[219];
//	v[2932] = v[2022] * v[265] + v[219] * v[2850];
//	v[3243] = v[219] * v[2441] + v[2932];
//	v[2601] = v[2022] * v[266];
//	v[2931] = v[2601] + v[218] * v[2850];
//	v[3242] = v[218] * v[2438] + v[2931];
//	v[2342] = v[128] * v[2020] + v[130] * v[2021] + v[132] * v[2022];
//	v[3166] = v[1590] * v[2003] + v[2342];
//	v[3161] = v[1584] * v[2004] + v[2342];
//	v[2023] = v[1013] * v[2003] + v[1014] * v[2004] + v[1015] * v[2005];
//	v[2851] = -(v[1203] * v[2023]);
//	v[2494] = v[2632] * v[2851];
//	v[2484] = v[2634] * v[2851];
//	v[2477] = v[240] * v[2851];
//	v[2072] = v[1221] * v[2023];
//	v[2068] = v[1222] * v[2023];
//	v[2058] = v[1994] * v[2023];
//	v[2024] = v[1017] * v[2003] + v[1018] * v[2004] + v[1019] * v[2005];
//	v[2495] = v[2024] * v[2691];
//	v[2485] = v[1964] * v[2024];
//	v[2924] = v[2484] + v[2485];
//	v[2478] = v[2024] * v[2852];
//	v[2923] = v[2477] + v[2478];
//	v[2076] = v[1219] * v[2024];
//	v[2952] = v[2058] + v[2076];
//	v[2062] = v[1993] * v[2024];
//	v[2956] = v[2062] + v[2072];
//	v[2025] = v[1021] * v[2003] + v[1022] * v[2004] + v[1023] * v[2005];
//	v[2492] = v[1965] * v[2025];
//	v[2925] = v[2492] + v[2495];
//	v[3229] = v[2494] + v[2925];
//	v[2488] = v[2025] * v[2692];
//	v[3227] = v[2485] + v[2488];
//	v[3226] = v[2488] + v[2924];
//	v[2480] = v[2025] * v[2694];
//	v[3224] = v[2477] + v[2480];
//	v[3222] = v[2480] + v[2923];
//	v[2066] = v[1992] * v[2025];
//	v[5967] = v[128] * v[2592] + v[2004] * v[579];
//	v[5968] = v[128] * v[2593] + v[2003] * v[579];
//	v[5969] = -(v[128] * v[2386]);
//	v[5970] = v[130] * v[2592] + v[2004] * v[580];
//	v[5971] = v[130] * v[2593] + v[2003] * v[580];
//	v[5972] = -(v[130] * v[2386]);
//	v[5973] = v[132] * v[2592] + v[2004] * v[581];
//	v[5974] = v[132] * v[2593] + v[2003] * v[581];
//	v[5975] = -(v[132] * v[2386]);
//	v[5976] = v[2388];
//	v[5977] = v[2387];
//	v[5978] = v[2386];
//	v[5979] = v[2066] + v[2024] * v[2630] + v[2023] * v[2633];
//	v[5980] = v[1219] * v[2023] + v[1220] * v[2025] + v[2062];
//	v[5981] = v[1221] * v[2024] + v[1222] * v[2025] + v[2058];
//	v[2957] = v[2066] + v[2068];
//	v[2857] = v[1220] * v[2024] + v[2066];
//	v[2948] = v[2068] + v[2857];
//	v[2856] = v[2062] + v[2025] * v[2630];
//	v[2951] = v[2072] + v[2856];
//	v[2855] = v[2058] + v[2025] * v[2633];
//	v[2955] = v[2076] + v[2855];
//	v[2026] = v[2388] * v[239];
//	v[2029] = v[2388] * v[249];
//	v[2030] = v[229] * v[2388];
//	v[2031] = v[229] * v[2387];
//	v[2034] = v[2387] * v[249];
//	v[2035] = v[2386] * v[249];
//	v[2036] = v[2387] * v[239];
//	v[2039] = v[229] * v[2386];
//	v[2040] = v[2386] * v[239];
//	v[2015] = v[1588] * v[2003] + v[1582] * v[2004] + v[1577] * v[2005] + v[2015] + v[1411] * v[2019] + v[1574] * v[2393];
//	v[2014] = v[1587] * v[2003] + v[1581] * v[2004] + v[1576] * v[2005] + v[2014] + v[1579] * v[2441] + v[2342] * v[265];
//	v[2013] = v[1586] * v[2003] + v[1580] * v[2004] + v[1575] * v[2005] + v[2013] + v[1585] * v[2438] + v[2342] * v[266];
//	v[2419] = v[2013] * v[250] + v[2014] * v[251] + v[2015] * v[252];
//	v[2016] = v[2016] + v[2419] * v[263];
//	v[2853] = v[2016] / v[597];
//	v[2041] = v[2015] * v[264] + v[252] * v[2853] + v[621];
//	v[2042] = v[2014] * v[264] + v[251] * v[2853] + v[620];
//	v[2043] = v[2013] * v[264] + v[250] * v[2853] + v[619];
//	v[6211] = 0e0;
//	v[6212] = 0e0;
//	v[6213] = 0e0;
//	v[6214] = 0e0;
//	v[6215] = 0e0;
//	v[6216] = 0e0;
//	v[6217] = v[2043];
//	v[6218] = v[2042];
//	v[6219] = v[2041];
//	v[6220] = 0e0;
//	v[6221] = 0e0;
//	v[6222] = 0e0;
//	v[6223] = 0e0;
//	v[6224] = 0e0;
//	v[6225] = 0e0;
//	v[6241] = 0e0;
//	v[6242] = 0e0;
//	v[6243] = 0e0;
//	v[6244] = v[2043];
//	v[6245] = v[2042];
//	v[6246] = v[2041];
//	v[6247] = 0e0;
//	v[6248] = 0e0;
//	v[6249] = 0e0;
//	v[6250] = 0e0;
//	v[6251] = 0e0;
//	v[6252] = 0e0;
//	v[6253] = 0e0;
//	v[6254] = 0e0;
//	v[6255] = 0e0;
//	v[6271] = v[2043];
//	v[6272] = v[2042];
//	v[6273] = v[2041];
//	v[6274] = 0e0;
//	v[6275] = 0e0;
//	v[6276] = 0e0;
//	v[6277] = 0e0;
//	v[6278] = 0e0;
//	v[6279] = 0e0;
//	v[6280] = 0e0;
//	v[6281] = 0e0;
//	v[6282] = 0e0;
//	v[6283] = 0e0;
//	v[6284] = 0e0;
//	v[6285] = 0e0;
//	v[2045] = v[152] * v[2041] + v[148] * v[2042] + v[144] * v[2043] - v[212] * v[2436] - v[213] * v[2437] + v[2440] * v[265]
//		+ v[2444] * v[267] + v[2005] * (v[267] * v[2879] + v[214] * v[587]);
//	v[2046] = v[1541] * v[2023];
//	v[2047] = v[1540] * v[2023];
//	v[2048] = v[1539] * v[2023];
//	v[2049] = v[1536] * v[2024];
//	v[2050] = v[2023] * v[236] + v[2024] * v[238];
//	v[2854] = v[177] * v[2050];
//	v[2587] = -(v[1203] * v[2854]);
//	v[2584] = v[2050] * v[2691];
//	v[3228] = v[238] * (v[2492] + v[2494]) + v[2584];
//	v[2051] = v[1534] * v[2024];
//	v[2052] = v[2046] + v[2049];
//	v[2053] = v[1535] * v[2024] + v[2854] / v[227];
//	v[2054] = v[1529] * v[2025];
//	v[2055] = v[2023] * v[230] + v[2025] * v[238];
//	v[2858] = v[182] * v[2055];
//	v[2586] = -(v[1203] * v[2858]);
//	v[2583] = v[2055] * v[2692];
//	v[3225] = v[2583] + v[238] * v[2924];
//	v[2056] = v[2024] * v[230] + v[2025] * v[236];
//	v[2859] = v[186] * v[2056];
//	v[3233] = -(v[2023] * v[2579]) - v[2024] * v[2580] - v[2025] * v[2581] + v[2629] * v[2854] + v[226] * v[2858]
//		+ v[225] * v[2859];
//	v[2585] = -(v[1203] * v[2859]);
//	v[2582] = v[2056] * v[2694];
//	v[3223] = v[2582] + v[236] * v[2923];
//	v[2577] = v[1538] * v[2023] + v[1533] * v[2024] + v[1528] * v[2025] + v[2048] + v[1369] * v[2050] + v[2051] + v[2054]
//		+ v[1367] * v[2055] + v[1363] * v[2056];
//	v[2057] = v[238] * v[2955] + v[2041] * v[718];
//	v[2060] = v[230] * v[2855] + v[2041] * v[722];
//	v[2061] = v[236] * v[2951] + v[2042] * v[720];
//	v[2064] = v[230] * v[2856] + v[2042] * v[722];
//	v[2065] = v[1220] * v[2050] + v[238] * v[2957] + v[2043] * v[718];
//	v[2067] = v[236] * v[2857] + v[2043] * v[720];
//	v[2070] = v[230] * v[2948] + v[2043] * v[722];
//	v[2071] = v[1531] * v[2025] + v[2858] / v[227];
//	v[2073] = v[2055] * v[2630] + v[238] * v[2956] + v[2042] * v[718];
//	v[2074] = v[2053] + v[2071];
//	v[2075] = v[1530] * v[2025] + v[2859] / v[227];
//	v[2077] = v[2056] * v[2633] + v[236] * v[2952] + v[2041] * v[720];
//	v[2078] = v[2047] + v[2075];
//	v[2079] = -(v[2026] * v[722]);
//	v[2080] = -(v[2026] * v[720]);
//	v[2081] = -(v[2026] * v[718]);
//	v[2082] = -(v[2029] * v[722]);
//	v[2083] = -(v[2029] * v[718]);
//	v[2084] = -(v[2029] * v[720]);
//	v[2862] = -2e0*v[2084];
//	v[2085] = -(v[2030] * v[720]);
//	v[2086] = -(v[2030] * v[718]);
//	v[2087] = -(v[2030] * v[722]);
//	v[2088] = -(v[2031] * v[722]);
//	v[2089] = -(v[2031] * v[720]);
//	v[2090] = -(v[2031] * v[718]);
//	v[2860] = -2e0*v[2090];
//	v[2091] = -(v[2034] * v[718]);
//	v[2092] = -(v[2034] * v[722]);
//	v[2863] = 2e0*v[2092];
//	v[2093] = -(v[2034] * v[720]);
//	v[2094] = -(v[2035] * v[720]);
//	v[2095] = -(v[2035] * v[722]);
//	v[2096] = -(v[2035] * v[718]);
//	v[2097] = v[2085] + v[2088] + v[2091] + v[2094] + 2e0*v[2051] * v[318] - v[2074] * v[321] - v[2052] * v[324];
//	v[2098] = -(v[2036] * v[722]);
//	v[2099] = -(v[2036] * v[718]);
//	v[2100] = -(v[2036] * v[720]);
//	v[2101] = -v[2080] - v[2083] - v[2095] - v[2098] - v[2074] * v[318] + 2e0*v[2054] * v[321] + v[2078] * v[324];
//	v[2102] = v[172] * v[2067] - v[2085] * v[311] + v[2080] * v[312] - v[2084] * v[313];
//	v[2103] = -(v[2039] * v[722]);
//	v[2104] = -(v[2039] * v[720]);
//	v[2861] = 2e0*v[2104];
//	v[2105] = -(v[2039] * v[718]);
//	v[2106] = -(v[2040] * v[720]);
//	v[2107] = -(v[2040] * v[722]);
//	v[2108] = -(v[2040] * v[718]);
//	v[2109] = -(v[199] * v[2041]) - v[196] * v[2042] - v[193] * v[2043] + v[2030] * v[362] + v[2026] * v[363]
//		+ v[2029] * v[364] + v[2031] * v[377] + v[2036] * v[378] + v[2034] * v[379] + v[2039] * v[392] + v[2040] * v[393]
//		+ v[2035] * v[394];
//	v[2110] = -(v[198] * v[2041]) - v[195] * v[2042] - v[192] * v[2043] + v[2030] * v[359] + v[2026] * v[360]
//		+ v[2029] * v[361] + v[2031] * v[374] + v[2036] * v[375] + v[2034] * v[376] + v[2039] * v[389] + v[2040] * v[390]
//		+ v[2035] * v[391];
//	v[2111] = -(v[197] * v[2041]) - v[194] * v[2042] - v[191] * v[2043] + v[2030] * v[356] + v[2026] * v[357]
//		+ v[2029] * v[358] + v[2031] * v[371] + v[2036] * v[372] + v[2034] * v[373] + v[2039] * v[386] + v[2040] * v[387]
//		+ v[2035] * v[388];
//	v[2112] = -v[2086] - v[2099] - v[2103] - v[2106] - v[2052] * v[318] + v[2078] * v[321] + 2e0*v[2048] * v[324];
//	v[2933] = -v[2112] / 2e0;
//	v[2113] = v[172] * v[2065] - v[2086] * v[311] + v[2081] * v[312] - v[2083] * v[313];
//	v[2114] = v[172] * v[2064] - v[2088] * v[311] + v[2098] * v[312] - v[2092] * v[313];
//	v[2115] = v[2082] + v[2093];
//	v[2116] = v[172] * v[2073] - v[2090] * v[311] + v[2099] * v[312] - v[2091] * v[313];
//	v[2117] = v[172] * v[2060] - v[2103] * v[311] + v[2107] * v[312] - v[2095] * v[313];
//	v[2118] = v[172] * v[2077] - v[2104] * v[311] + v[2106] * v[312] - v[2094] * v[313];
//	v[2119] = v[2089] + v[2105];
//	v[2120] = v[2079] + v[2108];
//	v[6181] = 0e0;
//	v[6182] = 0e0;
//	v[6183] = 0e0;
//	v[6184] = 0e0;
//	v[6185] = 0e0;
//	v[6186] = 0e0;
//	v[6187] = 0e0;
//	v[6188] = 0e0;
//	v[6189] = 0e0;
//	v[6190] = 0e0;
//	v[6191] = 0e0;
//	v[6192] = 0e0;
//	v[6193] = -v[2101] / 2e0 - v[2119];
//	v[6194] = (v[2097] - 2e0*v[2120]) / 2e0;
//	v[6195] = -v[2115] + v[2933];
//	v[2121] = (2e0*v[2081] + alphaB[1] * v[2097] - alphaB[0] * v[2101] - 2e0*v[2107] - alphaB[2] * v[2112]
//		+ 4e0*v[172] * v[2577] + v[2860] + v[2861] + v[2862] + v[2863] - v[2119] * v[305] - v[2120] * v[308] - v[2115] * v[310]
//		+ v[2070] * v[314] + 2e0*v[2067] * v[319] + 2e0*v[2065] * v[325] + 2e0*v[2064] * v[329] + v[2061] * v[333]
//		+ 2e0*v[2073] * v[338] + 2e0*v[2060] * v[342] + 2e0*v[2077] * v[346] + v[2057] * v[351]) / 2e0;
//	v[2122] = v[2111] * x1B[0] + v[2110] * x1B[1] + v[2109] * x1B[2];
//	v[2123] = (-v[2122] + v[2111] * x3B[0] + v[2110] * x3B[1] + v[2109] * x3B[2]) / 2e0;
//	v[2124] = (-v[2122] + v[2111] * x2B[0] + v[2110] * x2B[1] + v[2109] * x2B[2]) / 2e0;
//	v[2125] = (-2e0*v[2046] + 2e0*v[2049] - v[2087] * v[314] - 2e0*v[2085] * v[319] - 2e0*v[2086] * v[325]
//		- 2e0*v[2088] * v[329] - v[2089] * v[333] + v[2860] * v[338] - 2e0*v[2103] * v[342] - v[2861] * v[346] - v[2105] * v[351])
//		/ 2e0;
//	v[2938] = 8e0*v[2125];
//	v[2127] = -v[2047] + v[2075] + (v[2079] * v[314]) / 2e0 + v[2080] * v[319] + v[2081] * v[325] + v[2098] * v[329] +
//		(v[2100] * v[333]) / 2e0 + v[2099] * v[338] + v[2107] * v[342] + v[2106] * v[346] + (v[2108] * v[351]) / 2e0;
//	v[2937] = 8e0*v[2127];
//	v[2128] = (v[172] * v[2061] - v[2089] * v[311] + v[2100] * v[312] - v[2093] * v[313]) / 2e0;
//	v[2129] = (-2e0*v[2053] + 2e0*v[2071] - v[2082] * v[314] + v[2862] * v[319] - 2e0*v[2083] * v[325] - v[2863] * v[329]
//		- v[2093] * v[333] - 2e0*v[2091] * v[338] - 2e0*v[2095] * v[342] - 2e0*v[2094] * v[346] - v[2096] * v[351]) / 2e0;
//	v[3231] = v[2125] * v[305] - v[2127] * v[308] + v[2129] * v[310];
//	v[2936] = 8e0*v[2129];
//	v[6286] = 0e0;
//	v[6287] = 0e0;
//	v[6288] = 0e0;
//	v[6289] = 0e0;
//	v[6290] = 0e0;
//	v[6291] = 0e0;
//	v[6292] = 0e0;
//	v[6293] = 0e0;
//	v[6294] = 0e0;
//	v[6295] = 0e0;
//	v[6296] = 0e0;
//	v[6297] = 0e0;
//	v[6298] = v[2938];
//	v[6299] = -v[2937];
//	v[6300] = v[2936];
//	v[2136] = v[1272] * v[2125] + v[1270] * v[2127] - v[2128] + v[1267] * v[2129] - v[2121] * v[2934];
//	v[2612] = v[2136] + (-(v[172] * v[2070]) + v[2087] * v[311] - v[2079] * v[312] + v[2082] * v[313]) / 2e0;
//	v[2130] = (v[172] * v[2057] - v[2105] * v[311] + v[2108] * v[312] - v[2096] * v[313]) / 2e0;
//	v[2610] = v[2128] - v[2130] + v[2612];
//	v[2606] = -v[2130] + v[2136];
//	v[2608] = (v[2113] + v[2117]) / 2e0;
//	v[2132] = v[2116] + v[2118];
//	v[2611] = v[2132] / 2e0;
//	v[2133] = v[2102] + v[2114];
//	v[5952] = v[128] * v[2043];
//	v[5953] = v[128] * v[2042];
//	v[5954] = v[128] * v[2041];
//	v[5955] = v[130] * v[2043];
//	v[5956] = v[130] * v[2042];
//	v[5957] = v[130] * v[2041];
//	v[5958] = v[132] * v[2043];
//	v[5959] = v[132] * v[2042];
//	v[5960] = v[132] * v[2041];
//	v[5961] = -v[2043];
//	v[5962] = -v[2042];
//	v[5963] = -v[2041];
//	v[5964] = -v[2116] + v[2118] + 2e0*alphaB[0] * v[2606] + alphaB[2] * v[2608] + v[2101] * v[2824] + 2e0*
//		(v[2119] * v[2824] + v[2125] * v[2826]) + v[2133] * v[306];
//	v[5965] = (v[172] * v[2097]) / 2e0 + v[2113] - v[2117] + (alphaB[2] * v[2132]) / 2e0 + (alphaB[0] * v[2133]) / 2e0
//		+ 2e0*alphaB[1] * v[2610] + 2e0*(v[2120] * v[2824] - v[2127] * v[2826]);
//	v[5966] = -v[2102] + v[2114] + alphaB[0] * v[2608] + 2e0*alphaB[2] * v[2612] + v[172] * v[2933] + 2e0*
//		(v[2115] * v[2824] + v[2129] * v[2934]) + v[2132] * v[306];
//	v[2607] = v[2133] / 2e0;
//	v[2134] = (v[153] * v[2041] + v[149] * v[2042] + v[145] * v[2043] - v[2045] - v[215] * v[2436] - v[216] * v[2437]
//		+ v[2595] * v[265] + v[2598] * v[267] + v[2005] * (v[267] * v[2880] + v[217] * v[587])) / 2e0;
//	v[2135] = (v[154] * v[2041] + v[150] * v[2042] + v[146] * v[2043] - v[2045] - v[218] * v[2436] - v[219] * v[2437]
//		+ v[2601] * v[265] + v[2604] * v[267] + v[2005] * (v[267] * v[2881] + v[220] * v[587])) / 2e0;
//	for (i1998 = 1; i1998 <= 15; i1998++) {
//		i2866 = (i1998 == 14 ? 1 : 0);
//		v[2914] = (*a4)*i2866;
//		v[2869] = i2866 * v[172];
//		i2865 = (i1998 == 13 ? 1 : 0);
//		v[2915] = (*a4)*i2865;
//		v[2868] = -(i2865*v[172]);
//		i2864 = (i1998 == 15 ? 1 : 0);
//		v[2919] = (*a4)*i2864;
//		v[2867] = -(i2864*v[172]);
//		v[2152] = v[4135 + i1998];
//		v[2151] = v[4150 + i1998];
//		v[2142] = v[4120 + i1998];
//		v[2896] = v[2142] / 2e0;
//		v[2141] = v[4105 + i1998];
//		v[2895] = v[2141] / 2e0;
//		v[2143] = v[4199 + i1998];
//		v[2144] = v[4229 + i1998];
//		v[2145] = v[4214 + i1998];
//		v[2147] = v[4702 + i1998];
//		v[2148] = v[4184 + i1998];
//		v[2172] = 2e0*v[2148] * v[783];
//		v[2215] = -4e0*v[172] * v[2172];
//		v[2871] = v[2215] * v[227];
//		v[2190] = -2e0*v[2172];
//		v[2150] = v[4762 + i1998];
//		v[2153] = (-v[2151] - v[2152]) / 2e0;
//		v[2154] = (-v[2141] - v[2142]) / 2e0;
//		v[2888] = v[214] * v[2154] + (v[2142] * v[217]) / 2e0 + (v[2141] * v[220]) / 2e0;
//		v[2155] = i2864 + v[2143];
//		v[2157] = -i2864 + v[2143];
//		v[2158] = i2865 + v[2144];
//		v[2160] = -i2865 + v[2144];
//		v[2161] = -i2866 + v[2145];
//		v[2163] = i2866 + v[2145];
//		v[2164] = v[1267] * v[2148] + 8e0*i2864*v[783];
//		v[2165] = 2e0*alphaB[1] * i2866 - v[2148];
//		v[2166] = v[1270] * v[2148] - 8e0*i2866*v[783];
//		v[2167] = v[1272] * v[2148] + 8e0*i2865*v[783];
//		v[2870] = 2e0*(-v[2869] - (v[2148] * v[312]) / 2e0);
//		v[2169] = v[2868] + (v[2148] * v[311]) / 2e0;
//		v[2170] = v[2867] + (v[2148] * v[313]) / 2e0;
//		v[2171] = alphaB[2] * v[2172] + v[2867] / 2e0;
//		v[2173] = alphaB[0] * v[2172] + v[2868] / 2e0;
//		v[2174] = -(alphaB[1] * v[2172]) + v[2869] / 2e0;
//		v[2175] = (v[172] * v[2147]) / 2e0 - v[2172] * v[351];
//		v[2250] = v[2175] * v[238];
//		v[2176] = (-(v[2147] * v[313]) - v[2164] * v[351]) / 2e0;
//		v[2177] = (v[172] * v[2165]) / 2e0 - v[2172] * v[333];
//		v[2244] = v[2177] * v[236];
//		v[2178] = (v[2165] * v[312] + v[2166] * v[333]) / 2e0;
//		v[2179] = (v[172] * v[2150]) / 2e0 - v[2172] * v[314];
//		v[2239] = v[2179] * v[230];
//		v[2180] = (-(v[2150] * v[311]) - v[2167] * v[314]) / 2e0;
//		v[2181] = (2e0*v[2153] * x1B[0] + v[2152] * x2B[0] + v[2151] * x3B[0]) / 2e0;
//		v[2182] = (2e0*v[2153] * x1B[1] + v[2152] * x2B[1] + v[2151] * x3B[1]) / 2e0;
//		v[2183] = (2e0*v[2153] * x1B[2] + v[2152] * x2B[2] + v[2151] * x3B[2]) / 2e0;
//		v[2184] = (v[2870] + v[2147] * v[312] + v[2166] * v[351]) / 2e0;
//		v[2185] = (v[2870] + v[2150] * v[312] + v[2166] * v[314]) / 2e0;
//		v[2186] = v[2169] - (v[2147] * v[311]) / 2e0 - (v[2167] * v[351]) / 2e0;
//		v[2187] = v[2169] - (v[2165] * v[311]) / 2e0 - (v[2167] * v[333]) / 2e0;
//		v[2188] = v[2190] - v[2158] * v[311] - v[2167] * v[346];
//		v[2189] = v[172] * v[2158] - 2e0*v[2172] * v[346];
//		v[2191] = -v[2190] + v[2161] * v[312] + v[2166] * v[342];
//		v[2192] = v[172] * v[2161] - 2e0*v[2172] * v[342];
//		v[2252] = v[2192] * v[230];
//		v[2193] = -v[2190] - v[2160] * v[311] - v[2167] * v[338];
//		v[2194] = v[172] * v[2160] - 2e0*v[2172] * v[338];
//		v[2245] = v[2194] * v[238];
//		v[2195] = v[2170] - (v[2165] * v[313]) / 2e0 - (v[2164] * v[333]) / 2e0;
//		v[2196] = v[2170] - (v[2150] * v[313]) / 2e0 - (v[2164] * v[314]) / 2e0;
//		v[2197] = v[2190] - v[2155] * v[313] - v[2164] * v[329];
//		v[2198] = v[172] * v[2155] - 2e0*v[2172] * v[329];
//		v[2199] = v[2190] + v[2163] * v[312] + v[2166] * v[325];
//		v[2200] = v[172] * v[2163] - 2e0*v[2172] * v[325];
//		v[2240] = v[2200] * v[238];
//		v[2201] = -v[2171] + v[2158] * v[312] + v[2166] * v[346];
//		v[2202] = -v[2171] - v[2161] * v[311] - v[2167] * v[342];
//		v[2203] = -v[2171] + v[2160] * v[312] + v[2166] * v[338];
//		v[2204] = -v[2171] - v[2163] * v[311] - v[2167] * v[325];
//		v[2205] = v[2215] + 2e0*v[2171] * v[324];
//		v[2206] = v[2181] * v[387] + v[2182] * v[390] + v[2183] * v[393] - v[2184] * v[718] - v[2201] * v[720] - v[2191] * v[722];
//		v[2889] = v[2206] * v[239];
//		v[2207] = v[2181] * v[386] + v[2182] * v[389] + v[2183] * v[392] - v[2186] * v[718] - v[2188] * v[720] - v[2202] * v[722];
//		v[2890] = v[2207] * v[229];
//		v[2208] = -v[2190] - v[2157] * v[313] - v[2164] * v[319];
//		v[2209] = v[172] * v[2157] - 2e0*v[2172] * v[319];
//		v[2210] = -v[2173] + v[2155] * v[312] + v[2166] * v[329];
//		v[2211] = -v[2173] - v[2161] * v[313] - v[2164] * v[342];
//		v[2212] = -v[2173] - v[2163] * v[313] - v[2164] * v[325];
//		v[2213] = -v[2173] + v[2157] * v[312] + v[2166] * v[319];
//		v[2214] = v[2171] * v[321] + v[2173] * v[324];
//		v[2216] = v[2215] + 2e0*v[2173] * v[321];
//		v[2266] = v[2025] * v[2216];
//		v[2217] = v[2181] * v[372] + v[2182] * v[375] + v[2183] * v[378] - v[2203] * v[718] - v[2178] * v[720] - v[2210] * v[722];
//		v[2897] = v[2217] * v[239];
//		v[2218] = v[2174] - v[2158] * v[313] - v[2164] * v[346];
//		v[2219] = v[2174] - v[2160] * v[313] - v[2164] * v[338];
//		v[2220] = v[2174] - v[2155] * v[311] - v[2167] * v[329];
//		v[2221] = v[2174] - v[2157] * v[311] - v[2167] * v[319];
//		v[2222] = -(v[2173] * v[318]) - v[2174] * v[321];
//		v[2223] = -(v[2171] * v[318]) - v[2174] * v[324];
//		v[2224] = v[2215] + 2e0*v[2174] * v[318];
//		v[2270] = v[2024] * v[2224];
//		v[2225] = v[2181] * v[388] + v[2182] * v[391] + v[2183] * v[394] - v[2176] * v[718] - v[2218] * v[720] - v[2211] * v[722];
//		v[2891] = v[2225] * v[249];
//		v[2226] = v[2181] * v[373] + v[2182] * v[376] + v[2183] * v[379] - v[2219] * v[718] - v[2195] * v[720] - v[2197] * v[722];
//		v[2898] = v[2226] * v[249];
//		v[2227] = v[2181] * v[371] + v[2182] * v[374] + v[2183] * v[377] - v[2193] * v[718] - v[2187] * v[720] - v[2220] * v[722];
//		v[2899] = v[2227] * v[229];
//		v[2323] = v[2897] + v[2898] + v[2899];
//		v[2228] = v[2181] * v[356] + v[2182] * v[359] + v[2183] * v[362] - v[2204] * v[718] - v[2221] * v[720] - v[2180] * v[722];
//		v[2907] = v[2228] * v[229];
//		v[2229] = v[2181] * v[358] + v[2182] * v[361] + v[2183] * v[364] - v[2212] * v[718] - v[2208] * v[720] - v[2196] * v[722];
//		v[2908] = v[2229] * v[249];
//		v[2230] = v[2181] * v[357] + v[2182] * v[360] + v[2183] * v[363] - v[2199] * v[718] - v[2213] * v[720] - v[2185] * v[722];
//		v[2909] = v[2230] * v[239];
//		v[2327] = v[2907] + v[2908] + v[2909];
//		v[2231] = v[2166] + v[2214];
//		v[2264] = v[2025] * v[2231];
//		v[2232] = -v[2166] + v[2214];
//		v[2233] = (v[186] * v[2231] + v[2189] * v[225] + v[1363] * v[2871]) / v[227];
//		v[2234] = v[2164] + v[2222];
//		v[2272] = v[2025] * v[2234];
//		v[2235] = -v[2164] + v[2222];
//		v[2268] = v[2024] * v[2235];
//		v[2236] = (v[182] * v[2234] + v[2194] * v[226] + v[1367] * v[2871]) / v[227];
//		v[2237] = v[2239] + v[2209] * v[236];
//		v[2238] = v[2237] + v[2240] + v[2915];
//		v[2241] = v[2239] + v[2240];
//		v[2242] = v[2244] + v[2198] * v[230];
//		v[2243] = v[2242] + v[2245] + v[2914];
//		v[2246] = v[2244] + v[2245];
//		v[2247] = v[2042] * v[2177] - v[2036] * v[2178] - v[2031] * v[2187] - v[2039] * v[2188] + v[2041] * v[2189]
//			- v[2034] * v[2195] - v[2040] * v[2201] - v[2029] * v[2208] + v[2043] * v[2209] - v[2026] * v[2213] - v[2035] * v[2218]
//			- v[2030] * v[2221];
//		v[2248] = v[2043] * v[2179] - v[2030] * v[2180] - v[2026] * v[2185] - v[2040] * v[2191] + v[2041] * v[2192]
//			- v[2029] * v[2196] - v[2034] * v[2197] + v[2042] * v[2198] - v[2039] * v[2202] - v[2036] * v[2210] - v[2035] * v[2211]
//			- v[2031] * v[2220];
//		v[2249] = v[2250] + v[2252];
//		v[2251] = v[2250] + v[2189] * v[236];
//		v[2253] = v[2251] + v[2252] + v[2919];
//		v[2254] = v[2041] * v[2175] - v[2035] * v[2176] - v[2040] * v[2184] - v[2039] * v[2186] - v[2031] * v[2193]
//			+ v[2042] * v[2194] - v[2026] * v[2199] + v[2043] * v[2200] - v[2036] * v[2203] - v[2030] * v[2204] - v[2029] * v[2212]
//			- v[2034] * v[2219];
//		v[2255] = (v[177] * v[2235] + v[2200] * v[2629] + v[1369] * v[2871]) / v[227];
//		v[2256] = v[2266] + v[2268];
//		v[2953] = v[2256] / v[227];
//		v[2257] = v[2167] + v[2223];
//		v[2262] = v[2024] * v[2257];
//		v[2258] = -v[2167] + v[2223];
//		v[2259] = v[2270] + v[2272];
//		v[2949] = v[2259] / v[227];
//		v[2260] = v[2023] * v[2205] + v[2262] + v[2264];
//		v[2263] = v[2260] - v[2264];
//		v[2265] = v[2260] - v[2262];
//		v[2950] = v[2265] / v[227];
//		v[2267] = v[2023] * v[2232] + v[2266];
//		v[2269] = v[2267] + v[2268];
//		v[2271] = v[2023] * v[2258] + v[2270];
//		v[2273] = v[2271] + v[2272];
//		v[2275] = (v[146] * v[2141]) / 2e0 + (v[145] * v[2142]) / 2e0 + v[144] * v[2154] - v[191] * v[2181] - v[192] * v[2182]
//			- v[193] * v[2183] + v[4379 + i1998] + v[2200] * v[718] + v[2209] * v[720] + v[2179] * v[722];
//		v[2277] = (v[150] * v[2141]) / 2e0 + (v[149] * v[2142]) / 2e0 + v[148] * v[2154] - v[194] * v[2181] - v[195] * v[2182]
//			- v[196] * v[2183] + v[4394 + i1998] + v[2194] * v[718] + v[2177] * v[720] + v[2198] * v[722];
//		v[2279] = (v[154] * v[2141]) / 2e0 + (v[153] * v[2142]) / 2e0 + v[152] * v[2154] - v[197] * v[2181] - v[198] * v[2182]
//			- v[199] * v[2183] + v[4409 + i1998] + v[2175] * v[718] + v[2189] * v[720] + v[2192] * v[722];
//		v[2872] = v[2275] * v[250] + v[2277] * v[251] + v[2279] * v[252];
//		v[2420] = v[2872] / v[597];
//		v[2280] = v[2275];
//		v[2358] = v[2280];
//		v[2281] = v[2277];
//		v[2357] = v[2281];
//		v[2282] = v[2279];
//		v[2356] = v[2282];
//		v[2283] = -(v[1395] * v[2016] * v[2872]);
//		v[2284] = v[2420];
//		v[2873] = v[2284] * v[263];
//		v[2337] = v[2275] * v[264] + v[250] * v[2873];
//		v[2299] = v[2279] * v[264] + v[252] * v[2873];
//		v[2874] = v[2005] * v[2299];
//		v[2317] = v[249] * v[2874];
//		v[2314] = v[239] * v[2874];
//		v[2311] = v[229] * v[2874];
//		v[2288] = v[2277] * v[264] + v[251] * v[2873];
//		v[2875] = v[2004] * v[2288];
//		v[2316] = v[249] * v[2875];
//		v[2313] = v[239] * v[2875];
//		v[2310] = v[229] * v[2875];
//		v[2285] = v[2337];
//		v[2878] = v[2285] * v[266];
//		v[2877] = v[2288] * v[265] + v[2878];
//		v[2876] = v[2003] * v[2285];
//		v[2305] = v[249] * v[2876];
//		v[2303] = v[239] * v[2876];
//		v[2301] = v[229] * v[2876];
//		v[2286] = v[2141] * v[2287] + v[132] * v[2877];
//		v[2289] = v[2142] * v[2287] + v[130] * v[2877];
//		v[2290] = v[265] * (v[128] * v[2288] + v[2154] * v[266]) + v[128] * v[2878];
//		v[2291] = v[2875] + v[2876];
//		v[2292] = v[2301] + v[2310];
//		v[2293] = v[2303] + v[2313];
//		v[2294] = v[2305] + v[2316];
//		v[2300] = v[2874] + v[2876];
//		v[2302] = v[2301] + v[2311];
//		v[2304] = v[2303] + v[2314];
//		v[2306] = v[2305] + v[2317];
//		v[2309] = v[2874] + v[2875];
//		v[2312] = v[2310] + v[2311];
//		v[2315] = v[2313] + v[2314];
//		v[2318] = v[2316] + v[2317];
//		v[2320] = (i1998 == 12 ? (*a4) : 0e0) + v[2889] + v[2890] + v[2891];
//		v[2882] = v[1411] * v[2299] + v[1591] * v[2299] - v[2320] * v[267] + v[267] * v[2888] + (*a4)*v[4972 + i1998];
//		v[2322] = v[2323] + (*a4)*v[5002 + i1998];
//		v[2324] = (i1998 == 11 ? (*a4) : 0e0) + v[2323];
//		v[2326] = v[2327] + (*a4)*v[4987 + i1998];
//		v[2328] = (i1998 == 10 ? (*a4) : 0e0) + v[2327];
//		v[2330] = v[1528] * v[2215] + v[1529] * v[2216] + v[1530] * v[2231] + v[1531] * v[2234] + v[1992] * v[2238]
//			+ v[2233] * v[236] + v[2236] * v[238] + v[2242] * v[2630] + v[2249] * v[2633] + (*a4)*v[4927 + i1998];
//		v[2332] = v[1533] * v[2215] + v[1534] * v[2224] + v[1535] * v[2235] + v[1220] * v[2237] + v[1993] * v[2243]
//			+ v[1219] * v[2251] + v[1536] * v[2257] + v[2233] * v[230] + v[2255] * v[238] + (*a4)*v[4942 + i1998];
//		v[2334] = v[1539] * v[2205] + v[1538] * v[2215] + v[1540] * v[2232] + v[1222] * v[2241] + v[1221] * v[2246]
//			+ v[1994] * v[2253] + v[1541] * v[2258] + v[2236] * v[230] + v[2255] * v[236] + (*a4)*v[4957 + i1998];
//		b2335 = b598;
//		if (b2335) {
//			b2336 = b600;
//			if (b2336) {
//				v[2284] = 0e0;
//				v[2285] = 0e0;
//			}
//			else {
//			};
//		}
//		else {
//		};
//		v[2339] = -(v[2324] * v[2621]) - v[2322] * v[2850] + 2e0*v[1579] * v[2875] + v[2004] * (-(v[2328] * v[265]) + v[2882])
//			+ v[2154] * v[2913] + (v[2142] * v[2930]) / 2e0 + (v[2141] * v[2932]) / 2e0 + v[2337] * v[3161];
//		v[2370] = v[2339];
//		v[2341] = -(v[2019] * v[2320]) + v[2005] * (v[1579] * v[2288] + v[1585] * v[2337] - v[2326] * v[265] - v[2322] * v[266])
//			+ 2e0*v[1574] * v[2874] + v[2154] * v[3163] + v[2895] * v[3164] + v[2896] * v[3165] + (*a4)*v[6000 + i1998];
//		v[2369] = v[2341];
//		v[2343] = -(v[2328] * v[2619]) - v[2326] * v[2850] + 2e0*v[1585] * v[2876] + v[2003] * (-(v[2324] * v[266]) + v[2882])
//			+ v[2154] * v[2912] + (v[2142] * v[2929]) / 2e0 + (v[2141] * v[2931]) / 2e0 + v[2288] * v[3166];
//		v[2371] = v[2343];
//		v[2344] = 0e0;
//		v[2345] = 0e0;
//		v[2346] = 0e0;
//		v[2347] = 0e0;
//		v[2348] = 0e0;
//		v[2349] = 0e0;
//		v[2350] = 0e0;
//		b2351 = b598;
//		if (b2351) {
//			v[2352] = 0e0;
//			v[2353] = 0e0;
//			v[2354] = 0e0;
//			b2355 = b617;
//			if (b2355) {
//				v[2354] = v[2282];
//				v[2282] = 0e0;
//				v[2353] = v[2281];
//				v[2281] = 0e0;
//				v[2352] = v[2280];
//				v[2280] = 0e0;
//			}
//			else {
//				v[2349] = -v[2356];
//				v[2282] = 0e0;
//				v[2348] = -v[2357];
//				v[2281] = 0e0;
//				v[2347] = -v[2358];
//				v[2280] = 0e0;
//			};
//			v[2883] = (v[2352] * v[578] + v[2353] * v[586] + v[2354] * v[588])*(*zetan);
//			b2359 = b600;
//			if (b2359) {
//				v[2346] = v[2354] * v[612];
//				v[2345] = v[2353] * v[612];
//				v[2344] = v[2352] * v[612];
//				v[2350] = (v[2883] * v[3183] * v[610] * v[975]) / v[976];
//			}
//			else {
//				v[2346] = v[2354] * v[616];
//				v[2345] = v[2353] * v[616];
//				v[2344] = v[2352] * v[616];
//				v[2283] = v[2283] + (v[2883] * v[3184] * v[615] * v[980]) / v[981];
//			};
//		}
//		else {
//		};
//		v[2945] = v[2344] * v[265];
//		v[2926] = v[2344] * v[577];
//		v[2943] = v[2345] * v[266];
//		v[2927] = v[2345] * v[582];
//		v[2947] = v[2346] * v[587];
//		v[2944] = v[2346] * v[267];
//		v[2372] = v[2283];
//		v[2368] = v[2347];
//		v[2367] = v[2348];
//		v[2366] = v[2349];
//		b2364 = b598;
//		if (b2364) {
//			v[2885] = -(v[2368] * v[265]) - v[2367] * v[266] - v[2366] * v[267];
//			b2365 = b600;
//			if (b2365) {
//				v[2341] = v[2341] + v[2349] * v[2884];
//				v[2349] = 0e0;
//				v[2339] = v[2339] + v[2348] * v[2884];
//				v[2348] = 0e0;
//				v[2343] = v[2343] + v[2347] * v[2884];
//				v[2347] = 0e0;
//				v[2350] = v[2350] + v[2658] * v[2885] * v[3188];
//				v[2283] = v[2283] - v[2350];
//			}
//			else {
//				v[2341] = v[2369] + v[2366] * v[606];
//				v[2349] = 0e0;
//				v[2339] = v[2370] + v[2367] * v[606];
//				v[2348] = 0e0;
//				v[2343] = v[2371] + v[2368] * v[606];
//				v[2347] = 0e0;
//				v[2283] = v[2372] + (*n2)*v[2885] * v[3192] * v[594];
//			};
//		}
//		else {
//		};
//		v[2373] = v[2005] * v[2334] + v[2346] * v[249];
//		v[2910] = v[2373] * v[267];
//		v[2374] = v[2005] * v[2332] + v[2346] * v[239];
//		v[2903] = v[2374] * v[267];
//		v[2375] = v[2005] * v[2330] + v[229] * v[2346];
//		v[2901] = v[2375] * v[267];
//		v[2376] = v[2874] + v[2944];
//		v[2942] = v[2376] * v[266] + v[2927];
//		v[2941] = v[2376] * v[265] + v[2926];
//		v[2377] = v[2004] * v[2334] + v[2345] * v[249];
//		v[2911] = v[2377] * v[266];
//		v[2378] = v[2004] * v[2332] + v[2345] * v[239];
//		v[2904] = v[2378] * v[266];
//		v[2379] = v[2004] * v[2330] + v[229] * v[2345];
//		v[2902] = v[2379] * v[266];
//		v[2380] = v[2003] * v[2334] + v[2344] * v[249];
//		v[2900] = v[2380] * v[265];
//		v[2381] = v[2003] * v[2332] + v[2344] * v[239];
//		v[2893] = v[2381] * v[265];
//		v[2382] = v[2003] * v[2330] + v[229] * v[2344];
//		v[2892] = v[2382] * v[265];
//		v[2383] = v[213] * v[2344] + v[212] * v[2345] + (*a4)*v[6015 + i1998];
//		v[2384] = v[216] * v[2344] + v[215] * v[2345] + (*a4)*v[6030 + i1998];
//		v[2385] = v[219] * v[2344] + v[218] * v[2345] + (*a4)*v[6045 + i1998];
//		v[2906] = v[128] * v[2383] + v[130] * v[2384] + v[132] * v[2385];
//		v[2389] = v[1013] * v[2344] + v[1014] * v[2345] + v[1015] * v[2346] + v[2225] * v[2386] + v[2226] * v[2387]
//			+ v[2229] * v[2388] + v[2299] * v[3196] + v[2288] * v[3197] + v[2337] * v[3198];
//		v[2482] = v[2389] * v[2626];
//		v[2390] = v[1017] * v[2344] + v[1018] * v[2345] + v[1019] * v[2346] + v[2206] * v[2386] + v[2217] * v[2387]
//			+ v[2230] * v[2388] + v[2299] * v[3199] + v[2288] * v[3200] + v[2337] * v[3201];
//		v[2921] = v[227] * v[2390];
//		v[2486] = v[2390] * v[2625];
//		v[2391] = v[1021] * v[2344] + v[1022] * v[2345] + v[1023] * v[2346] + v[2207] * v[2386] + v[2227] * v[2387]
//			+ v[2228] * v[2388] + v[2299] * v[3202] + v[2288] * v[3203] + v[2337] * v[3204];
//		v[2922] = v[227] * v[2391];
//		v[2489] = v[2391] * v[2624];
//		v[2392] = v[2291] + v[2943] + v[2945];
//		v[2946] = v[2392] * v[267];
//		v[2928] = v[2946] + v[2947];
//		v[6361] = v[2004] * v[2290] + v[2154] * v[2592] + v[2337] * (v[128] * v[2438] + v[2613]) + v[128] * v[2941]
//			+ v[2345] * v[579];
//		v[6362] = v[2003] * v[2290] + v[2154] * v[2593] + v[2288] * (v[128] * v[2441] + v[2613]) + v[128] * v[2942]
//			+ v[2344] * v[579];
//		v[6363] = -(v[2154] * v[2386]) + v[2299] * (2e0*v[2613] + v[2614]) + v[128] * v[2928];
//		v[6364] = v[2004] * v[2289] + v[2337] * (v[130] * v[2438] + v[2615]) + v[2592] * v[2896] + v[130] * v[2941]
//			+ v[2345] * v[580];
//		v[6365] = v[2003] * v[2289] + v[2288] * (v[130] * v[2441] + v[2615]) + v[2593] * v[2896] + v[130] * v[2942]
//			+ v[2344] * v[580];
//		v[6366] = v[2299] * (2e0*v[2615] + v[2616]) - v[2386] * v[2896] + v[130] * v[2928];
//		v[6367] = v[2004] * v[2286] + v[2337] * (v[132] * v[2438] + v[2617]) + v[2592] * v[2895] + v[132] * v[2941]
//			+ v[2345] * v[581];
//		v[6368] = v[2003] * v[2286] + v[2288] * (v[132] * v[2441] + v[2617]) + v[2593] * v[2895] + v[132] * v[2942]
//			+ v[2344] * v[581];
//		v[6369] = v[2299] * (2e0*v[2617] + v[2618]) - v[2386] * v[2895] + v[132] * v[2928];
//		v[6370] = v[2337] * (-v[2438] + v[2887]) - v[2926] + v[265] * (-v[2309] - v[2943] - v[2944]);
//		v[6371] = v[2288] * (-v[2441] + v[2886]) - v[2927] + v[266] * (-v[2300] - v[2944] - v[2945]);
//		v[6372] = v[2299] * (-v[2393] - v[2619] - v[2621]) - v[2946] - v[2947];
//		v[6373] = v[2024] * v[2233] + v[2023] * v[2236] + v[1992] * v[2391] + v[2215] * (v[1967] * v[2023] + v[1969] * v[2024]
//			+ v[2492]) + v[2390] * v[2630] + v[2389] * v[2633] + v[2192] * v[2855] + v[2198] * v[2856] + v[2269] * v[2917]
//			+ v[2179] * v[2948] + v[179] * v[2949] + v[184] * v[2950];
//		v[6374] = v[2025] * v[2233] + v[2023] * v[2255] + v[1219] * v[2389] + v[1993] * v[2390] + v[1220] * v[2391] + v[2215] *
//			(v[1966] * v[2023] + v[1970] * v[2025] + v[2485]) + v[2209] * v[2857] + v[2273] * v[2916] + v[2177] * v[2951]
//			+ v[2189] * v[2952] + v[176] * v[2953] + v[2263] * v[2954];
//		v[6375] = v[2025] * v[2236] + v[2024] * v[2255] + v[1994] * v[2389] + v[1221] * v[2390] + v[1222] * v[2391] + v[2215] *
//			(v[1968] * v[2024] + v[1971] * v[2025] + v[2477]) + v[2260] * v[2920] + v[2175] * v[2955] + v[2194] * v[2956]
//			+ v[2200] * v[2957] + v[2267] * v[2958] + v[2271] * v[2959];
//		v[2394] = -(v[267] * (v[2294] + v[2900] + v[2911])) + v[2299] * v[3205] - v[2373] * v[587];
//		v[2395] = -(v[267] * (v[2293] + v[2893] + v[2904])) + v[2299] * v[3206] - v[2374] * v[587];
//		v[2341] = -(v[223] * v[2291]) + v[2341] + v[1752] * v[2346] + v[1813] * v[2373] + v[1817] * v[2374] + v[1821] * v[2375]
//			+ v[1411] * v[2392] - v[2292] * v[407] - v[2293] * v[408] - v[2294] * v[409] + v[266] * (-(v[223] * v[2345])
//				- v[2379] * v[407] - v[2378] * v[408] - v[2377] * v[409]) + v[265] * (-(v[223] * v[2344]) - v[2382] * v[407]
//					- v[2381] * v[408] - v[2380] * v[409]) + 2e0*v[267] * (v[2346] * v[2398] + v[2005] * (v[2888] - v[2889] - v[2890]
//						- v[2891]) - v[2375] * v[407] - v[2374] * v[408] - v[2373] * v[409] + (*a4)*v[5985 + i1998]);
//		v[2399] = -(v[267] * (v[2292] + v[2892] + v[2902])) + v[2299] * v[3208] - v[2375] * v[587];
//		v[2402] = -(v[266] * (v[2302] + v[2892] + v[2901])) + v[2288] * v[3209] - v[2379] * v[582];
//		v[2405] = -(v[266] * (v[2304] + v[2893] + v[2903])) + v[2288] * v[3210] - v[2378] * v[582];
//		v[2339] = -(v[222] * v[2300]) + v[2339] + v[2345] * v[2406] + v[2377] * v[2644] + v[2378] * v[2646] + v[2379] * v[2648]
//			+ v[2376] * v[2894] - v[2302] * v[404] - v[2304] * v[405] - v[2306] * v[406] + v[267] * (-(v[222] * v[2346])
//				- v[2375] * v[404] - v[2374] * v[405] - v[2373] * v[406]) + v[265] * (-(v[222] * v[2344]) + v[2906] - v[2382] * v[404]
//					- v[2381] * v[405] - v[2380] * v[406]) + 2e0*v[266] * (v[2345] * v[2407] + v[2004] * (v[213] * v[2154] + v[219] * v[2895]
//						+ v[216] * v[2896] - v[2897] - v[2898] - v[2899]) - v[2379] * v[404] - v[2378] * v[405] - v[2377] * v[406] + (*a4)*v[6060
//						+ i1998]);
//		v[2410] = -(v[266] * (v[2306] + v[2900] + v[2910])) + v[2288] * v[3212] - v[2377] * v[582];
//		v[2412] = -(v[265] * (v[2312] + v[2901] + v[2902])) + v[2337] * v[3213] - v[2382] * v[577];
//		v[2414] = -(v[265] * (v[2315] + v[2903] + v[2904])) + v[2337] * v[3214] - v[2381] * v[577];
//		v[2343] = -(v[221] * v[2309]) + v[2343] + v[2344] * v[2415] - v[2380] * v[2645] - v[2381] * v[2647] - v[2382] * v[2649]
//			+ v[2376] * v[2905] - v[2312] * v[401] - v[2315] * v[402] - v[2318] * v[403] + v[267] * (-(v[221] * v[2346])
//				- v[2375] * v[401] - v[2374] * v[402] - v[2373] * v[403]) + v[266] * (-(v[221] * v[2345]) + v[2906] - v[2379] * v[401]
//					- v[2378] * v[402] - v[2377] * v[403]) + 2e0*v[265] * (v[2344] * v[2416] + v[2003] * ((v[2142] * v[215]) / 2e0
//						+ v[212] * v[2154] + (v[2141] * v[218]) / 2e0 - v[2907] - v[2908] - v[2909]) - v[2382] * v[401] - v[2381] * v[402]
//						- v[2380] * v[403] + (*a4)*v[6075 + i1998]);
//		v[2418] = -(v[265] * (v[2318] + v[2910] + v[2911])) + v[2337] * v[3216] - v[2380] * v[577];
//		v[2283] = v[2283] + v[2419] * v[2420] * v[262] + (v[2013] * v[2275] + v[2014] * v[2277] + v[2015] * v[2279]
//			+ v[2343] * v[250] + v[2339] * v[251] + v[2341] * v[252])*v[263];
//		v[2422] = v[2341] * v[264] + v[2015] * v[2873] + (v[2016] * v[2279] + v[2283] * v[252]) / v[597];
//		v[2424] = v[2339] * v[264] + v[2014] * v[2873] + (v[2016] * v[2277] + v[2283] * v[251]) / v[597];
//		v[2426] = v[2343] * v[264] + v[2013] * v[2873] + (v[2016] * v[2275] + v[2283] * v[250]) / v[597];
//		v[2427] = -(v[2041] * v[2183]) - v[168] * v[2422];
//		v[2428] = -(v[2041] * v[2182]) - v[167] * v[2422];
//		v[2429] = -(v[2041] * v[2181]) - v[166] * v[2422];
//		v[2430] = -(v[2042] * v[2183]) - v[168] * v[2424];
//		v[2431] = -(v[2042] * v[2182]) - v[167] * v[2424];
//		v[2432] = -(v[2042] * v[2181]) - v[166] * v[2424];
//		v[2433] = -(v[2043] * v[2183]) - v[168] * v[2426];
//		v[2434] = -(v[2043] * v[2182]) - v[167] * v[2426];
//		v[2435] = -(v[2043] * v[2181]) - v[166] * v[2426];
//		v[2445] = 2e0*v[2287] * v[2383] + v[152] * v[2422] + v[148] * v[2424] + v[144] * v[2426] + v[2376] * v[2879]
//			+ v[212] * v[2926] + v[213] * v[2927] + v[214] * v[2928] + v[2299] * v[3219] + v[2337] * v[3220] + v[2288] * v[3221] + (*a4
//				)*v[6255 + i1998] + v[6270 + i1998];
//		v[2446] = v[236] * v[2389] + v[2023] * v[2914];
//		v[2447] = v[230] * v[2389] + v[2023] * v[2915];
//		v[2449] = v[1127] * v[2389] + v[2023] * (v[2246] / v[227] + v[2215] * v[2448]) + v[2446] * v[2916];
//		v[2451] = v[1133] * v[2389] + v[2023] * (v[2241] / v[227] + v[2215] * v[2450]) + v[2447] * v[2917];
//		v[2453] = (v[186] * v[2446] + v[184] * v[2447] + v[2023] * (v[2253] + v[2452] * v[2871]) + v[2389] * v[2918]) / v[227];
//		v[2454] = v[238] * v[2390] + v[2024] * v[2919];
//		v[2455] = v[230] * v[2390] + v[2024] * v[2915];
//		v[2457] = v[1122] * v[2390] + v[2024] * (v[2251] / v[227] + v[2215] * v[2456]) + v[2454] * v[2920];
//		v[2458] = v[2446] + v[2454];
//		v[2459] = v[2449] + v[2457];
//		v[2461] = (v[2050] * v[2200] + v[175] * v[2455] + v[177] * v[2458] + v[2587] * v[2871] + v[2024] * (v[2237]
//			+ v[2460] * v[2871]) + v[1132] * v[2921]) / v[227];
//		v[2463] = (v[182] * v[2454] + v[179] * v[2455] + v[2024] * (v[2243] + v[2462] * v[2871]) + v[1126] * v[2921]) / v[227];
//		v[2464] = v[236] * v[2391] + v[2025] * v[2914];
//		v[2465] = v[238] * v[2391] + v[2025] * v[2919];
//		v[2466] = v[2455] + v[2464];
//		v[2468] = (v[176] * v[2464] + v[177] * v[2465] + v[2025] * (v[2238] + v[2467] * v[2871]) + v[1130] * v[2922]) / v[227];
//		v[2469] = v[2447] + v[2465];
//		v[2471] = (v[2055] * v[2194] + v[181] * v[2464] + v[182] * v[2469] + v[2586] * v[2871] + v[2025] * (v[2242]
//			+ v[2470] * v[2871]) + v[1137] * v[2922]) / v[227];
//		v[2472] = v[2461] + v[2471];
//		v[2474] = (v[2056] * v[2189] + v[187] * v[2465] + v[186] * v[2466] + v[2585] * v[2871] + v[2025] * (v[2249]
//			+ v[2473] * v[2871]) + v[1140] * v[2922]) / v[227];
//		v[2475] = v[2451] + v[2474];
//		v[2476] = QBi[2][2] * v[2427] + QBi[2][1] * v[2428] + QBi[2][0] * v[2429] + v[1219] * v[2454] + v[240] * v[2482]
//			+ v[2465] * v[2633] + v[238] * (v[2260] / v[227] + v[2215] * v[3222]);
//		v[2479] = (v[2056] * v[2231]) / v[227] + QBi[1][2] * v[2427] + QBi[1][1] * v[2428] + QBi[1][0] * v[2429]
//			+ v[1994] * v[2446] + v[237] * v[2486] + v[2263] * v[2625] + v[2466] * v[2633] + v[2215] * v[3223];
//		v[2481] = QBi[0][2] * v[2427] + QBi[0][1] * v[2428] + QBi[0][0] * v[2429] + v[1994] * v[2447] + v[225] * v[2489]
//			+ v[230] * (v[2950] + v[2215] * v[3224]);
//		v[2483] = (v[2055] * v[2234]) / v[227] + QBi[2][2] * v[2430] + QBi[2][1] * v[2431] + QBi[2][0] * v[2432]
//			+ v[1993] * v[2454] + v[2271] * v[2626] + v[2469] * v[2630] + v[2482] * v[2634] + v[2215] * v[3225];
//		v[2487] = QBi[1][2] * v[2430] + QBi[1][1] * v[2431] + QBi[1][0] * v[2432] + v[1221] * v[2446] + v[231] * v[2486]
//			+ v[2464] * v[2630] + v[236] * (v[2273] / v[227] + v[2215] * v[3226]);
//		v[2490] = QBi[0][2] * v[2430] + QBi[0][1] * v[2431] + QBi[0][0] * v[2432] + v[1993] * v[2455] + v[226] * v[2489]
//			+ v[230] * (v[2949] + v[2215] * v[3227]);
//		v[2491] = (v[2050] * v[2235]) / v[227] + QBi[2][2] * v[2433] + QBi[2][1] * v[2434] + QBi[2][0] * v[2435]
//			+ v[1220] * v[2458] + v[1992] * v[2465] + v[2267] * v[2626] + v[2482] * v[2632] + v[2215] * v[3228];
//		v[2493] = QBi[1][2] * v[2433] + QBi[1][1] * v[2434] + QBi[1][0] * v[2435] + v[1992] * v[2464] + v[2486] * v[2629]
//			+ v[236] * (v[2215] * v[2925] + v[2953]);
//		v[2496] = QBi[0][2] * v[2433] + QBi[0][1] * v[2434] + QBi[0][0] * v[2435] + v[1222] * v[2447] + v[1220] * v[2455]
//			+ v[228] * v[2489] + v[230] * (v[2269] / v[227] + v[2215] * v[3229]);
//		v[2497] = v[2026] * v[2181] + v[166] * v[2414];
//		v[2498] = v[2026] * v[2182] + v[167] * v[2414];
//		v[2499] = v[2026] * v[2183] + v[168] * v[2414];
//		v[2500] = QBi[0][0] * v[2497] + QBi[0][1] * v[2498] + QBi[0][2] * v[2499];
//		v[2501] = QBi[1][0] * v[2497] + QBi[1][1] * v[2498] + QBi[1][2] * v[2499];
//		v[2502] = QBi[2][0] * v[2497] + QBi[2][1] * v[2498] + QBi[2][2] * v[2499];
//		v[2503] = v[2029] * v[2181] + v[166] * v[2418];
//		v[2504] = v[2029] * v[2182] + v[167] * v[2418];
//		v[2505] = v[2029] * v[2183] + v[168] * v[2418];
//		v[2506] = QBi[0][0] * v[2503] + QBi[0][1] * v[2504] + QBi[0][2] * v[2505];
//		v[2507] = QBi[2][0] * v[2503] + QBi[2][1] * v[2504] + QBi[2][2] * v[2505];
//		v[2508] = QBi[1][0] * v[2503] + QBi[1][1] * v[2504] + QBi[1][2] * v[2505];
//		v[2509] = v[2030] * v[2181] + v[166] * v[2412];
//		v[2510] = v[2030] * v[2182] + v[167] * v[2412];
//		v[2511] = v[2030] * v[2183] + v[168] * v[2412];
//		v[2512] = QBi[1][0] * v[2509] + QBi[1][1] * v[2510] + QBi[1][2] * v[2511];
//		v[2513] = QBi[2][0] * v[2509] + QBi[2][1] * v[2510] + QBi[2][2] * v[2511];
//		v[2514] = QBi[0][0] * v[2509] + QBi[0][1] * v[2510] + QBi[0][2] * v[2511];
//		v[2515] = v[2031] * v[2181] + v[166] * v[2402];
//		v[2516] = v[2031] * v[2182] + v[167] * v[2402];
//		v[2517] = v[2031] * v[2183] + v[168] * v[2402];
//		v[2518] = QBi[0][0] * v[2515] + QBi[0][1] * v[2516] + QBi[0][2] * v[2517];
//		v[2519] = QBi[1][0] * v[2515] + QBi[1][1] * v[2516] + QBi[1][2] * v[2517];
//		v[2520] = QBi[2][0] * v[2515] + QBi[2][1] * v[2516] + QBi[2][2] * v[2517];
//		v[2521] = v[2034] * v[2181] + v[166] * v[2410];
//		v[2522] = v[2034] * v[2182] + v[167] * v[2410];
//		v[2523] = v[2034] * v[2183] + v[168] * v[2410];
//		v[2524] = QBi[2][0] * v[2521] + QBi[2][1] * v[2522] + QBi[2][2] * v[2523];
//		v[2525] = QBi[0][0] * v[2521] + QBi[0][1] * v[2522] + QBi[0][2] * v[2523];
//		v[2526] = QBi[1][0] * v[2521] + QBi[1][1] * v[2522] + QBi[1][2] * v[2523];
//		v[2527] = v[2035] * v[2181] + v[166] * v[2394];
//		v[2528] = v[2035] * v[2182] + v[167] * v[2394];
//		v[2529] = v[2035] * v[2183] + v[168] * v[2394];
//		v[2530] = QBi[1][0] * v[2527] + QBi[1][1] * v[2528] + QBi[1][2] * v[2529];
//		v[2531] = QBi[0][0] * v[2527] + QBi[0][1] * v[2528] + QBi[0][2] * v[2529];
//		v[2532] = QBi[2][0] * v[2527] + QBi[2][1] * v[2528] + QBi[2][2] * v[2529];
//		v[2533] = -(v[2052] * v[2171]) - v[2074] * v[2173] + 2e0*v[2051] * v[2174] + v[2512] + v[2518] + v[2524] + v[2530]
//			+ 2e0*v[2463] * v[318] - v[2472] * v[321] - v[2459] * v[324];
//		v[2534] = v[2036] * v[2181] + v[166] * v[2405];
//		v[2535] = v[2036] * v[2182] + v[167] * v[2405];
//		v[2536] = v[2036] * v[2183] + v[168] * v[2405];
//		v[2537] = QBi[0][0] * v[2534] + QBi[0][1] * v[2535] + QBi[0][2] * v[2536];
//		v[2538] = QBi[2][0] * v[2534] + QBi[2][1] * v[2535] + QBi[2][2] * v[2536];
//		v[2539] = QBi[1][0] * v[2534] + QBi[1][1] * v[2535] + QBi[1][2] * v[2536];
//		v[2540] = v[2078] * v[2171] + 2e0*v[2054] * v[2173] - v[2074] * v[2174] - v[2501] - v[2507] - v[2531] - v[2537]
//			- v[2472] * v[318] + 2e0*v[2468] * v[321] + v[2475] * v[324];
//		v[2541] = -(v[2084] * v[2164]) + v[2080] * v[2166] - v[2085] * v[2167] - 2e0*v[2067] * v[2172] + v[172] * v[2493]
//			- v[2512] * v[311] + v[2501] * v[312] - v[2508] * v[313];
//		v[2542] = v[2039] * v[2181] + v[166] * v[2399];
//		v[2543] = v[2039] * v[2182] + v[167] * v[2399];
//		v[2544] = v[2039] * v[2183] + v[168] * v[2399];
//		v[2545] = QBi[0][0] * v[2542] + QBi[0][1] * v[2543] + QBi[0][2] * v[2544];
//		v[2546] = QBi[1][0] * v[2542] + QBi[1][1] * v[2543] + QBi[1][2] * v[2544];
//		v[2547] = QBi[2][0] * v[2542] + QBi[2][1] * v[2543] + QBi[2][2] * v[2544];
//		v[2548] = v[2040] * v[2181] + v[166] * v[2395];
//		v[2549] = v[2040] * v[2182] + v[167] * v[2395];
//		v[2550] = v[2040] * v[2183] + v[168] * v[2395];
//		v[2551] = QBi[1][0] * v[2548] + QBi[1][1] * v[2549] + QBi[1][2] * v[2550];
//		v[2552] = QBi[0][0] * v[2548] + QBi[0][1] * v[2549] + QBi[0][2] * v[2550];
//		v[2553] = QBi[2][0] * v[2548] + QBi[2][1] * v[2549] + QBi[2][2] * v[2550];
//		v[2554] = 2e0*v[2048] * v[2171] + v[2078] * v[2173] - v[2052] * v[2174] - v[2513] - v[2538] - v[2545] - v[2551]
//			- v[2459] * v[318] + v[2475] * v[321] + 2e0*v[2453] * v[324];
//		v[2555] = -(v[2083] * v[2164]) + v[2081] * v[2166] - v[2086] * v[2167] - 2e0*v[2065] * v[2172] + v[172] * v[2491]
//			- v[2513] * v[311] + v[2502] * v[312] - v[2507] * v[313];
//		v[2556] = -(v[2092] * v[2164]) + v[2098] * v[2166] - v[2088] * v[2167] - 2e0*v[2064] * v[2172] + v[172] * v[2490]
//			- v[2518] * v[311] + v[2537] * v[312] - v[2525] * v[313];
//		v[2557] = v[2506] + v[2526];
//		v[2558] = -(v[2091] * v[2164]) + v[2099] * v[2166] - v[2090] * v[2167] - 2e0*v[2073] * v[2172] + v[172] * v[2483]
//			- v[2520] * v[311] + v[2538] * v[312] - v[2524] * v[313];
//		v[2559] = -(v[2095] * v[2164]) + v[2107] * v[2166] - v[2103] * v[2167] - 2e0*v[2060] * v[2172] + v[172] * v[2481]
//			- v[2545] * v[311] + v[2552] * v[312] - v[2531] * v[313];
//		v[2560] = -(v[2094] * v[2164]) + v[2106] * v[2166] - v[2104] * v[2167] - 2e0*v[2077] * v[2172] + v[172] * v[2479]
//			- v[2546] * v[311] + v[2551] * v[312] - v[2530] * v[313];
//		v[2561] = v[2519] + v[2547];
//		v[2562] = v[2500] + v[2553];
//		v[2563] = -(QBi[1][2] * v[2247]) - QBi[0][2] * v[2248] - QBi[2][2] * v[2254] - v[199] * v[2422] - v[196] * v[2424]
//			- v[193] * v[2426] + v[2412] * v[362] + v[2414] * v[363] + v[2418] * v[364] + v[2402] * v[377] + v[2405] * v[378]
//			+ v[2410] * v[379] + v[2399] * v[392] + v[2395] * v[393] + v[2394] * v[394];
//		v[2564] = -(QBi[1][1] * v[2247]) - QBi[0][1] * v[2248] - QBi[2][1] * v[2254] - v[198] * v[2422] - v[195] * v[2424]
//			- v[192] * v[2426] + v[2412] * v[359] + v[2414] * v[360] + v[2418] * v[361] + v[2402] * v[374] + v[2405] * v[375]
//			+ v[2410] * v[376] + v[2399] * v[389] + v[2395] * v[390] + v[2394] * v[391];
//		v[2565] = -(QBi[1][0] * v[2247]) - QBi[0][0] * v[2248] - QBi[2][0] * v[2254] - v[197] * v[2422] - v[194] * v[2424]
//			- v[191] * v[2426] + v[2412] * v[356] + v[2414] * v[357] + v[2418] * v[358] + v[2402] * v[371] + v[2405] * v[372]
//			+ v[2410] * v[373] + v[2399] * v[386] + v[2395] * v[387] + v[2394] * v[388];
//		v[2566] = v[2565] * x1B[0] + v[2564] * x1B[1] + v[2563] * x1B[2];
//		v[2567] = (-v[2566] + v[2565] * x3B[0] + v[2564] * x3B[1] + v[2563] * x3B[2]) / 2e0;
//		v[2568] = (-v[2566] + v[2565] * x2B[0] + v[2564] * x2B[1] + v[2563] * x2B[2]) / 2e0;
//		v[2569] = (-(v[2105] * v[2147]) - v[2087] * v[2150] - 2e0*v[2088] * v[2155] - 2e0*v[2085] * v[2157]
//			- 2e0*v[2104] * v[2158] - 2e0*v[2090] * v[2160] - 2e0*v[2103] * v[2161] - 2e0*v[2086] * v[2163] - v[2089] * v[2165]
//			- 2e0*v[2449] + 2e0*v[2457] - v[2514] * v[314] - 2e0*v[2512] * v[319] - 2e0*v[2513] * v[325] - 2e0*v[2518] * v[329]
//			- v[2519] * v[333] - 2e0*v[2520] * v[338] - 2e0*v[2545] * v[342] - 2e0*v[2546] * v[346] - v[2547] * v[351]) / 2e0;
//		v[2570] = (-(v[2082] * v[2164]) + v[2079] * v[2166] - v[2087] * v[2167] - 2e0*v[2070] * v[2172] + v[172] * v[2496]
//			- v[2514] * v[311] + v[2500] * v[312] - v[2506] * v[313]) / 2e0;
//		v[2572] = (v[2108] * v[2147] + v[2079] * v[2150] + 2e0*v[2098] * v[2155] + 2e0*v[2080] * v[2157]
//			+ 2e0*v[2106] * v[2158] + 2e0*v[2099] * v[2160] + 2e0*v[2107] * v[2161] + 2e0*v[2081] * v[2163] + v[2100] * v[2165]
//			- 2e0*v[2451] + 2e0*v[2474] + v[2500] * v[314] + 2e0*v[2501] * v[319] + 2e0*v[2502] * v[325] + 2e0*v[2537] * v[329]
//			+ v[2539] * v[333] + 2e0*v[2538] * v[338] + 2e0*v[2552] * v[342] + 2e0*v[2551] * v[346] + v[2553] * v[351]) / 2e0;
//		v[2573] = (-(v[2093] * v[2164]) + v[2100] * v[2166] - v[2089] * v[2167] - 2e0*v[2061] * v[2172] + v[172] * v[2487]
//			- v[2519] * v[311] + v[2539] * v[312] - v[2526] * v[313]) / 2e0;
//		v[2574] = (-(v[2096] * v[2147]) - v[2082] * v[2150] - 2e0*v[2092] * v[2155] - 2e0*v[2084] * v[2157]
//			- 2e0*v[2094] * v[2158] - 2e0*v[2091] * v[2160] - 2e0*v[2095] * v[2161] - 2e0*v[2083] * v[2163] - v[2093] * v[2165]
//			- 2e0*v[2461] + 2e0*v[2471] - v[2506] * v[314] - 2e0*v[2508] * v[319] - 2e0*v[2507] * v[325] - 2e0*v[2525] * v[329]
//			- v[2526] * v[333] - 2e0*v[2524] * v[338] - 2e0*v[2531] * v[342] - 2e0*v[2530] * v[346] - v[2532] * v[351]) / 2e0;
//		v[2609] = v[1270] * v[2572] - v[2573] + v[1267] * v[2574] + 24e0*v[2148] * v[2576] * v[3231] + v[2935] * (-
//			(v[2057] * v[2147]) - v[2070] * v[2150] - 2e0*v[2064] * v[2155] - 2e0*v[2067] * v[2157] - 2e0*v[2077] * v[2158]
//			- 2e0*v[2073] * v[2160] - 2e0*v[2060] * v[2161] - 2e0*v[2065] * v[2163] - v[2061] * v[2165] - 2e0*v[2502]
//			+ 2e0*v[2508] + 2e0*v[2520] - 2e0*v[2525] - alphaB[1] * v[2533] + alphaB[0] * v[2540] - 2e0*v[2546] + 2e0*v[2552]
//			+ alphaB[2] * v[2554] + 8e0*v[2172] * v[2577] + v[2561] * v[305] + v[2562] * v[308] + v[2557] * v[310] - v[2496] * v[314]
//			- 2e0*v[2493] * v[319] - 4e0*v[172] * (v[1873] * v[2256] + v[1877] * v[2259] + v[1861] * v[2260] + v[1860] * v[2263]
//				+ v[1880] * v[2265] + v[1874] * v[2267] + v[1865] * v[2269] + v[1868] * v[2271] + v[1867] * v[2273] + v[1966] * v[2446]
//				+ v[1967] * v[2447] + v[2453] + v[1968] * v[2454] + v[1969] * v[2455] + v[1369] * v[2458] + v[2463] + v[1970] * v[2464]
//				+ v[1971] * v[2465] + v[1363] * v[2466] + v[2468] + v[1367] * v[2469] + v[2253] * v[2477] + v[2251] * v[2478]
//				+ v[2249] * v[2480] + v[2246] * v[2484] + v[2243] * v[2485] + v[2242] * v[2488] + v[2238] * v[2492] + v[2241] * v[2494]
//				+ v[2237] * v[2495] + v[2189] * v[2582] + v[2194] * v[2583] + v[2200] * v[2584] + v[2231] * v[2585] + v[2234] * v[2586]
//				+ v[2235] * v[2587] + v[2389] * v[2812] + v[2390] * v[2813] + v[2391] * v[2814] + 2e0*v[2215] * v[2578] * v[3233])
//			- 2e0*v[2491] * v[325] - 2e0*v[2490] * v[329] - v[2487] * v[333] - 2e0*v[2483] * v[338] - 2e0*v[2481] * v[342]
//			- 2e0*v[2479] * v[346] - v[2476] * v[351] - 2e0*v[6180 + i1998]) - 2e0*v[1215] * (-4e0*v[2121] * v[2148]
//				+ 4e0*v[2569] * v[305] + v[6285 + i1998]);
//		v[2940] = v[2609] + (v[2096] * v[2164] - v[2108] * v[2166] + v[2105] * v[2167] + 2e0*v[2057] * v[2172]
//			- v[172] * v[2476] + v[2547] * v[311] - v[2553] * v[312] + v[2532] * v[313]) / 2e0;
//		v[2939] = (v[2555] + v[2559]) / 2e0;
//		v[2590] = v[2558] + v[2560];
//		v[2591] = v[2541] + v[2556];
//		v[6301] = 0e0;
//		v[6302] = 0e0;
//		v[6303] = 0e0;
//		v[6304] = 0e0;
//		v[6305] = 0e0;
//		v[6306] = 0e0;
//		v[6307] = 0e0;
//		v[6308] = 0e0;
//		v[6309] = 0e0;
//		v[6310] = 0e0;
//		v[6311] = 0e0;
//		v[6312] = 0e0;
//		v[6313] = 2e0*v[2606];
//		v[6314] = v[2607];
//		v[6315] = v[2608];
//		v[6316] = 0e0;
//		v[6317] = 0e0;
//		v[6318] = 0e0;
//		v[6319] = 0e0;
//		v[6320] = 0e0;
//		v[6321] = 0e0;
//		v[6322] = 0e0;
//		v[6323] = 0e0;
//		v[6324] = 0e0;
//		v[6325] = 0e0;
//		v[6326] = 0e0;
//		v[6327] = 0e0;
//		v[6328] = v[2607];
//		v[6329] = 2e0*v[2610];
//		v[6330] = v[2611];
//		v[6331] = 0e0;
//		v[6332] = 0e0;
//		v[6333] = 0e0;
//		v[6334] = 0e0;
//		v[6335] = 0e0;
//		v[6336] = 0e0;
//		v[6337] = 0e0;
//		v[6338] = 0e0;
//		v[6339] = 0e0;
//		v[6340] = 0e0;
//		v[6341] = 0e0;
//		v[6342] = 0e0;
//		v[6343] = v[2608];
//		v[6344] = v[2611];
//		v[6345] = 2e0*v[2612];
//		v[6346] = v[2043] * v[2154] + v[128] * v[2426];
//		v[6347] = v[2042] * v[2154] + v[128] * v[2424];
//		v[6348] = v[2041] * v[2154] + v[128] * v[2422];
//		v[6349] = v[130] * v[2426] + v[2043] * v[2896];
//		v[6350] = v[130] * v[2424] + v[2042] * v[2896];
//		v[6351] = v[130] * v[2422] + v[2041] * v[2896];
//		v[6352] = v[132] * v[2426] + v[2043] * v[2895];
//		v[6353] = v[132] * v[2424] + v[2042] * v[2895];
//		v[6354] = v[132] * v[2422] + v[2041] * v[2895];
//		v[6355] = -v[2426];
//		v[6356] = -v[2424];
//		v[6357] = -v[2422];
//		v[6358] = v[2101] * v[2172] - v[2558] + v[2560] + v[2540] * v[2824] + 2e0*(v[2561] * v[2824] + v[2569] * v[2826]
//			+ v[2148] * (v[2119] * v[2827] - v[1215] * v[2938])) + alphaB[2] * v[2939] + 2e0*alphaB[0] * v[2940] + v[2591] * v[306]
//			+ v[6300 + i1998];
//		v[6359] = -(v[2097] * v[2172]) + v[2555] - v[2559] + (alphaB[2] * v[2590]) / 2e0 + (alphaB[0] * v[2591]) / 2e0
//			- v[2533] * v[2824] + 2e0*alphaB[1] * (-v[2570] + v[2573] + v[2940]) + v[6315 + i1998] + 2e0*(v[2562] * v[2824]
//				+ v[2148] * (v[2120] * v[2827] + v[1215] * v[2937]) - 4e0*v[2572] * v[783]);
//		v[6360] = v[2112] * v[2172] - v[2541] + v[2556] + 2e0*alphaB[2] * (-v[2570] + v[2609]) + v[2554] * v[2824] + 2e0*
//			(v[2557] * v[2824] + v[2574] * v[2826] + v[2148] * (v[2115] * v[2935] - v[1215] * v[2936])) + alphaB[0] * v[2939]
//			+ v[2590] * v[306] + v[6330 + i1998];
//		v[2599] = (2e0*v[2287] * v[2384] + v[153] * v[2422] + v[149] * v[2424] + v[145] * v[2426] - v[2445] + v[2376] * v[2880]
//			+ v[215] * v[2926] + v[216] * v[2927] + v[217] * v[2928] + v[2299] * v[3236] + v[2337] * v[3237] + v[2288] * v[3238] + (*a4
//				)*v[6225 + i1998] + v[6240 + i1998]) / 2e0;
//		v[2605] = (2e0*v[2287] * v[2385] + v[154] * v[2422] + v[150] * v[2424] + v[146] * v[2426] - v[2445] + v[2376] * v[2881]
//			+ v[218] * v[2926] + v[219] * v[2927] + v[220] * v[2928] + v[2299] * v[3241] + v[2337] * v[3242] + v[2288] * v[3243] + (*a4
//				)*v[6195 + i1998] + v[6210 + i1998]) / 2e0;
//		Rc[i1998 - 1] += v[2135] * v[2141] + v[2134] * v[2142] + v[2123] * v[2151] + v[2124] * v[2152] + v[5951 + i1998] + (*a4
//			)*v[5966 + i1998];
//		for (i2139 = 1; i2139 <= 15; i2139++) {
//			Kc[i1998 - 1][i2139 - 1] += v[2605] * v[4105 + i2139] + v[2599] * v[4120 + i2139] + v[2568] * v[4135 + i2139]
//				+ v[2567] * v[4150 + i2139] + v[6345 + i2139] + (*a4)*v[6360 + i2139];
//		};/* end for */
//	};/* end for */
//#pragma endregion

#pragma region
double v01; double v010; double v011; double v012; double v013; double v014;
double v015; double v016; double v017; double v018; double v019; double v02;
double v020; double v021; double v022; double v023; double v024; double v025;
double v026; double v027; double v03; double v04; double v05; double v06; double v07;
double v08; double v09;
int i710, i773, i971, i1444, i2518, i2733, i3444, i3445, i3446, i3447, i3448, i3449, i3456
, i3457, i3458, i3459, i3460, i3461, i3462, i3463, i3464, b29, b30, b277, b542
, b603, b605, b622, b641, b675
, b715, b716, b720, b829, b830, b841, b976, b977, b978, b979, b1055, b1059, b1060
, b1077, b1078, b1209, b1240, b1293, b1632, b1695, b1741, b1975, b1976, b1985, b1986, b2013
, b2014, b2015, b2031, b2144, b2148, b2152, b2163, b2164, b2258, b2280, b2325, b2527, b2528
, b2532, b2537, b2538, b2626, b2894, b3025, b3026, b3034, b3038, b3042, b3052, b3053
, b3148;
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
v[600] = v[32] - v[33];
v[34] = (*n1);
v[3369] = v[31] * v[34];
v[615] = -1e0 + v[34];
v[1987] = -1e0 + v[615];
v[35] = (*n2);
v[620] = -1e0 + v[35];
v[1990] = -1e0 + v[620];
v[36] = xi1A[0];
v[37] = xi1A[1];
v[38] = xi1A[2];
v[39] = xi2A[0];
v[40] = xi2A[1];
v[41] = xi2A[2];
v[42] = xi3A[0];
v[43] = xi3A[1];
v[44] = xi3A[2];
v[45] = u1A[0];
v[144] = v[36] + v[45];
v[286] = -v[144] / 2e0;
v[46] = u1A[1];
v[148] = v[37] + v[46];
v[288] = -v[148] / 2e0;
v[47] = u1A[2];
v[152] = v[38] + v[47];
v[290] = -v[152] / 2e0;
v[48] = u2A[0];
v[145] = v[39] + v[48];
v[3375] = v[145] / 2e0;
v[282] = v[286] + v[3375];
v[49] = u2A[1];
v[149] = v[40] + v[49];
v[3377] = v[149] / 2e0;
v[283] = v[288] + v[3377];
v[50] = u2A[2];
v[153] = v[41] + v[50];
v[3379] = v[153] / 2e0;
v[284] = v[290] + v[3379];
v[51] = u3A[0];
v[146] = v[42] + v[51];
v[3376] = v[146] / 2e0;
v[287] = v[286] + v[3376];
v[52] = u3A[1];
v[150] = v[43] + v[52];
v[3378] = v[150] / 2e0;
v[289] = v[288] + v[3378];
v[53] = u3A[2];
v[154] = v[44] + v[53];
v[3380] = v[154] / 2e0;
v[291] = v[290] + v[3380];
v[129] = -cAp[0] / 2e0;
v[131] = -cAp[1] / 2e0;
v[134] = -cAi[0] / 2e0;
v[136] = -cAi[1] / 2e0;
v[76] = x1B[0];
v[300] = -v[76] / 2e0;
v[77] = x1B[1];
v[302] = -v[77] / 2e0;
v[78] = x1B[2];
v[304] = -v[78] / 2e0;
v[79] = x2B[0];
v[293] = v[300] + v[79] / 2e0;
v[80] = x2B[1];
v[294] = v[302] + v[80] / 2e0;
v[81] = x2B[2];
v[295] = v[304] + v[81] / 2e0;
v[82] = x3B[0];
v[301] = v[300] + v[82] / 2e0;
v[83] = x3B[1];
v[303] = v[302] + v[83] / 2e0;
v[84] = x3B[2];
v[305] = v[304] + v[84] / 2e0;
v[85] = xBi[0];
v[86] = xBi[1];
v[87] = xBi[2];
v[88] = QBi[0][0];
v[89] = QBi[0][1];
v[90] = QBi[0][2];
v[91] = QBi[1][0];
v[92] = QBi[1][1];
v[93] = QBi[1][2];
v[94] = QBi[2][0];
v[95] = QBi[2][1];
v[96] = QBi[2][2];
v[97] = uB[0];
v[98] = uB[1];
v[99] = uB[2];
v[100] = alphaB[0];
v[312] = v[100] / 2e0;
v[310] = 2e0*v[100];
v[3709] = 4e0*v[310];
v[3563] = -v[310] / 2e0;
v[180] = (v[100] * v[100]);
v[101] = alphaB[1];
v[313] = 2e0*v[101];
v[5178] = 0e0;
v[5179] = 0e0;
v[5180] = 0e0;
v[5181] = 0e0;
v[5182] = 0e0;
v[5183] = 0e0;
v[5184] = 0e0;
v[5185] = 0e0;
v[5186] = 0e0;
v[5187] = 0e0;
v[5188] = 0e0;
v[5189] = 0e0;
v[5190] = v[310];
v[5191] = v[313];
v[5192] = 0e0;
v[3564] = -v[313] / 2e0;
v[5546] = 0e0;
v[5547] = 0e0;
v[5548] = 0e0;
v[5549] = 0e0;
v[5550] = 0e0;
v[5551] = 0e0;
v[5552] = 0e0;
v[5553] = 0e0;
v[5554] = 0e0;
v[5555] = 0e0;
v[5556] = 0e0;
v[5557] = 0e0;
v[5558] = -v[310];
v[5559] = -v[313];
v[5560] = 0e0;
v[311] = v[101] / 2e0;
v[5013] = 0e0;
v[5014] = 0e0;
v[5015] = 0e0;
v[5016] = 0e0;
v[5017] = 0e0;
v[5018] = 0e0;
v[5019] = 0e0;
v[5020] = 0e0;
v[5021] = 0e0;
v[5022] = 0e0;
v[5023] = 0e0;
v[5024] = 0e0;
v[5025] = v[311];
v[5026] = v[312];
v[5027] = 0e0;
v[178] = v[100] * v[311];
v[173] = (v[101] * v[101]);
v[356] = -v[173] - v[180];
v[3426] = v[356] / 2e0;
v[102] = alphaB[2];
v[334] = v[102] + v[178];
v[3432] = -2e0*v[334];
v[324] = -v[102] + v[178];
v[3430] = -2e0*v[324];
v[315] = 2e0*v[102];
v[5118] = 0e0;
v[5119] = 0e0;
v[5120] = 0e0;
v[5121] = 0e0;
v[5122] = 0e0;
v[5123] = 0e0;
v[5124] = 0e0;
v[5125] = 0e0;
v[5126] = 0e0;
v[5127] = 0e0;
v[5128] = 0e0;
v[5129] = 0e0;
v[5130] = v[310];
v[5131] = 0e0;
v[5132] = v[315];
v[5058] = 0e0;
v[5059] = 0e0;
v[5060] = 0e0;
v[5061] = 0e0;
v[5062] = 0e0;
v[5063] = 0e0;
v[5064] = 0e0;
v[5065] = 0e0;
v[5066] = 0e0;
v[5067] = 0e0;
v[5068] = 0e0;
v[5069] = 0e0;
v[5070] = 0e0;
v[5071] = v[313];
v[5072] = v[315];
v[3562] = -v[315] / 2e0;
v[5606] = 0e0;
v[5607] = 0e0;
v[5608] = 0e0;
v[5609] = 0e0;
v[5610] = 0e0;
v[5611] = 0e0;
v[5612] = 0e0;
v[5613] = 0e0;
v[5614] = 0e0;
v[5615] = 0e0;
v[5616] = 0e0;
v[5617] = 0e0;
v[5618] = 0e0;
v[5619] = -v[313];
v[5620] = -v[315];
v[4998] = 0e0;
v[4999] = 0e0;
v[5000] = 0e0;
v[5001] = 0e0;
v[5002] = 0e0;
v[5003] = 0e0;
v[5004] = 0e0;
v[5005] = 0e0;
v[5006] = 0e0;
v[5007] = 0e0;
v[5008] = 0e0;
v[5009] = 0e0;
v[5010] = v[310];
v[5011] = v[313];
v[5012] = v[315];
v[314] = v[102] / 2e0;
v[5028] = 0e0;
v[5029] = 0e0;
v[5030] = 0e0;
v[5031] = 0e0;
v[5032] = 0e0;
v[5033] = 0e0;
v[5034] = 0e0;
v[5035] = 0e0;
v[5036] = 0e0;
v[5037] = 0e0;
v[5038] = 0e0;
v[5039] = 0e0;
v[5040] = v[314];
v[5041] = 0e0;
v[5042] = v[312];
v[5043] = 0e0;
v[5044] = 0e0;
v[5045] = 0e0;
v[5046] = 0e0;
v[5047] = 0e0;
v[5048] = 0e0;
v[5049] = 0e0;
v[5050] = 0e0;
v[5051] = 0e0;
v[5052] = 0e0;
v[5053] = 0e0;
v[5054] = 0e0;
v[5055] = 0e0;
v[5056] = v[314];
v[5057] = v[311];
v[185] = v[101] * v[314];
v[351] = v[100] + v[185];
v[3433] = -2e0*v[351];
v[343] = -v[100] + v[185];
v[3431] = -2e0*v[343];
v[183] = v[100] * v[314];
v[347] = -v[101] + v[183];
v[3434] = -2e0*v[347];
v[330] = v[101] + v[183];
v[3429] = -2e0*v[330];
v[174] = (v[102] * v[102]);
v[898] = 4e0 + v[173] + v[174] + v[180];
v[3710] = 24e0 / Power(v[898], 4);
v[1435] = 1e0 / Power(v[898], 3);
v[3708] = -2e0*v[1435];
v[3584] = -4e0*v[1435];
v[3329] = 8e0*v[1435];
v[1500] = -(v[310] * v[3329]);
v[1498] = v[313] * v[3329];
v[1495] = -(v[315] * v[3329]);
v[803] = 1e0 / (v[898] * v[898]);
v[3330] = 4e0*v[803];
v[3715] = 2e0*v[803];
v[3450] = 8e0*v[803];
v[338] = -v[174] - v[180];
v[3427] = v[338] / 2e0;
v[319] = -v[173] - v[174];
v[3428] = v[319] / 2e0;
v[318] = v[315] * v[3330];
v[3334] = -v[318] / 2e0;
v[360] = v[3334] * v[356];
v[317] = -(v[313] * v[3330]);
v[3331] = v[317] / 2e0;
v[340] = v[3331] * v[338];
v[316] = v[310] * v[3330];
v[3333] = -v[316] / 2e0;
v[320] = v[319] * v[3333];
v[230] = v[10] * v[100] + dalphaiB[0] * v[11] + ddalphaiB[0] * v[12];
v[236] = v[10] * v[101] + dalphaiB[1] * v[11] + ddalphaiB[1] * v[12];
v[238] = v[10] * v[102] + dalphaiB[2] * v[11] + ddalphaiB[2] * v[12];
v[157] = -cBp[0] / 2e0;
v[159] = -cBp[1] / 2e0;
v[162] = -cBi[0] / 2e0;
v[164] = -cBi[1] / 2e0;
v[128] = v[129] + v[131];
v[130] = 0.5e0 - v[129];
v[435] = v[130] * v[291];
v[434] = v[130] * v[289];
v[433] = v[130] * v[287];
v[132] = 0.5e0 - v[131];
v[6214] = 0e0;
v[6215] = 0e0;
v[6216] = -v[128];
v[6217] = 0e0;
v[6218] = 0e0;
v[6219] = -v[130];
v[6220] = 0e0;
v[6221] = 0e0;
v[6222] = -v[132];
v[6223] = 0e0;
v[6224] = 0e0;
v[6225] = 0e0;
v[6226] = 0e0;
v[6227] = 0e0;
v[6228] = 0e0;
v[6199] = 0e0;
v[6200] = -v[128];
v[6201] = 0e0;
v[6202] = 0e0;
v[6203] = -v[130];
v[6204] = 0e0;
v[6205] = 0e0;
v[6206] = -v[132];
v[6207] = 0e0;
v[6208] = 0e0;
v[6209] = 0e0;
v[6210] = 0e0;
v[6211] = 0e0;
v[6212] = 0e0;
v[6213] = 0e0;
v[6184] = -v[128];
v[6185] = 0e0;
v[6186] = 0e0;
v[6187] = -v[130];
v[6188] = 0e0;
v[6189] = 0e0;
v[6190] = -v[132];
v[6191] = 0e0;
v[6192] = 0e0;
v[6193] = 0e0;
v[6194] = 0e0;
v[6195] = 0e0;
v[6196] = 0e0;
v[6197] = 0e0;
v[6198] = 0e0;
v[5223] = 0e0;
v[5224] = 0e0;
v[5225] = v[128];
v[5226] = 0e0;
v[5227] = 0e0;
v[5228] = v[130];
v[5229] = 0e0;
v[5230] = 0e0;
v[5231] = v[132];
v[5232] = 0e0;
v[5233] = 0e0;
v[5234] = -1e0;
v[5235] = 0e0;
v[5236] = 0e0;
v[5237] = 0e0;
v[5208] = 0e0;
v[5209] = v[128];
v[5210] = 0e0;
v[5211] = 0e0;
v[5212] = v[130];
v[5213] = 0e0;
v[5214] = 0e0;
v[5215] = v[132];
v[5216] = 0e0;
v[5217] = 0e0;
v[5218] = -1e0;
v[5219] = 0e0;
v[5220] = 0e0;
v[5221] = 0e0;
v[5222] = 0e0;
v[5193] = v[128];
v[5194] = 0e0;
v[5195] = 0e0;
v[5196] = v[130];
v[5197] = 0e0;
v[5198] = 0e0;
v[5199] = v[132];
v[5200] = 0e0;
v[5201] = 0e0;
v[5202] = -1e0;
v[5203] = 0e0;
v[5204] = 0e0;
v[5205] = 0e0;
v[5206] = 0e0;
v[5207] = 0e0;
v[426] = v[132] * v[284];
v[425] = v[132] * v[283];
v[424] = v[132] * v[282];
v[133] = v[134] + v[136];
v[135] = 0.5e0 - v[134];
v[137] = 0.5e0 - v[136];
v[156] = v[157] + v[159];
v[158] = 0.5e0 - v[157];
v[160] = 0.5e0 - v[159];
v[161] = v[162] + v[164];
v[163] = 0.5e0 - v[162];
v[165] = 0.5e0 - v[164];
v[166] = v[156] * v[76] + v[158] * v[79] + v[160] * v[82];
v[167] = v[156] * v[77] + v[158] * v[80] + v[160] * v[83];
v[168] = v[156] * v[78] + v[158] * v[81] + v[160] * v[84];
v[742] = -(v[166] * v[88]) - v[167] * v[89] - v[168] * v[90];
v[740] = -(v[166] * v[91]) - v[167] * v[92] - v[168] * v[93];
v[738] = -(v[166] * v[94]) - v[167] * v[95] - v[168] * v[96];
v[169] = v[161] * v[76] + v[163] * v[79] + v[165] * v[82];
v[170] = v[161] * v[77] + v[163] * v[80] + v[165] * v[83];
v[171] = v[161] * v[78] + v[163] * v[81] + v[165] * v[84];
v[172] = 4e0 / v[898];
v[3443] = 2e0*v[172];
v[3332] = -v[172] / 2e0;
v[358] = v[313] * v[3332];
v[359] = v[3331] * v[356] + v[358];
v[355] = v[310] * v[3332];
v[357] = v[355] + v[3333] * v[356];
v[352] = v[172] - v[316] * v[351];
v[349] = -v[172] + v[317] * v[347];
v[344] = -v[172] - v[316] * v[343];
v[341] = v[315] * v[3332];
v[342] = v[3334] * v[338] + v[341];
v[339] = v[3333] * v[338] + v[355];
v[337] = v[172] - v[318] * v[334];
v[332] = v[172] + v[317] * v[330];
v[329] = -(v[172] * v[314]);
v[3452] = 2e0*v[329];
v[353] = -v[329] + v[317] * v[351];
v[398] = v[349] * v[90] + v[353] * v[93] + v[359] * v[96];
v[395] = v[349] * v[89] + v[353] * v[92] + v[359] * v[95];
v[392] = v[349] * v[88] + v[353] * v[91] + v[359] * v[94];
v[413] = v[166] * v[392] + v[167] * v[395] + v[168] * v[398];
v[348] = -v[329] - v[316] * v[347];
v[397] = v[348] * v[90] + v[352] * v[93] + v[357] * v[96];
v[394] = v[348] * v[89] + v[352] * v[92] + v[357] * v[95];
v[391] = v[348] * v[88] + v[352] * v[91] + v[357] * v[94];
v[412] = v[166] * v[391] + v[167] * v[394] + v[168] * v[397];
v[345] = -v[329] + v[317] * v[343];
v[331] = -v[329] - v[316] * v[330];
v[328] = -v[172] - v[318] * v[324];
v[326] = -(v[172] * v[312]);
v[3453] = 2e0*v[326];
v[350] = -v[326] - v[318] * v[347];
v[336] = -v[326] + v[317] * v[334];
v[383] = v[336] * v[90] + v[340] * v[93] + v[345] * v[96];
v[380] = v[336] * v[89] + v[340] * v[92] + v[345] * v[95];
v[377] = v[336] * v[88] + v[340] * v[91] + v[345] * v[94];
v[410] = v[166] * v[377] + v[167] * v[380] + v[168] * v[383];
v[333] = -v[326] - v[318] * v[330];
v[327] = v[317] * v[324] - v[326];
v[323] = v[172] * v[311];
v[3454] = 2e0*v[323];
v[354] = v[323] - v[318] * v[351];
v[399] = v[350] * v[90] + v[354] * v[93] + v[360] * v[96];
v[396] = v[350] * v[89] + v[354] * v[92] + v[360] * v[95];
v[393] = v[350] * v[88] + v[354] * v[91] + v[360] * v[94];
v[414] = v[166] * v[393] + v[167] * v[396] + v[168] * v[399];
v[346] = v[323] - v[318] * v[343];
v[384] = v[337] * v[90] + v[342] * v[93] + v[346] * v[96];
v[381] = v[337] * v[89] + v[342] * v[92] + v[346] * v[95];
v[378] = v[337] * v[88] + v[342] * v[91] + v[346] * v[94];
v[411] = v[166] * v[378] + v[167] * v[381] + v[168] * v[384];
v[335] = v[323] - v[316] * v[334];
v[382] = v[335] * v[90] + v[339] * v[93] + v[344] * v[96];
v[379] = v[335] * v[89] + v[339] * v[92] + v[344] * v[95];
v[376] = v[335] * v[88] + v[339] * v[91] + v[344] * v[94];
v[409] = v[166] * v[376] + v[167] * v[379] + v[168] * v[382];
v[325] = v[323] - v[316] * v[324];
v[367] = v[320] * v[90] + v[325] * v[93] + v[331] * v[96];
v[364] = v[320] * v[89] + v[325] * v[92] + v[331] * v[95];
v[361] = v[320] * v[88] + v[325] * v[91] + v[331] * v[94];
v[406] = v[166] * v[361] + v[167] * v[364] + v[168] * v[367];
v[439] = -(v[287] * v[406]) - v[289] * v[409] - v[291] * v[412];
v[427] = -(v[282] * v[406]) - v[283] * v[409] - v[284] * v[412];
v[322] = v[319] * v[3334] + v[341];
v[369] = v[322] * v[90] + v[328] * v[93] + v[333] * v[96];
v[366] = v[322] * v[89] + v[328] * v[92] + v[333] * v[95];
v[363] = v[322] * v[88] + v[328] * v[91] + v[333] * v[94];
v[408] = v[166] * v[363] + v[167] * v[366] + v[168] * v[369];
v[441] = -(v[287] * v[408]) - v[289] * v[411] - v[291] * v[414];
v[429] = -(v[282] * v[408]) - v[283] * v[411] - v[284] * v[414];
v[321] = v[319] * v[3331] + v[358];
v[368] = v[321] * v[90] + v[327] * v[93] + v[332] * v[96];
v[365] = v[321] * v[89] + v[327] * v[92] + v[332] * v[95];
v[362] = v[321] * v[88] + v[327] * v[91] + v[332] * v[94];
v[407] = v[166] * v[362] + v[167] * v[365] + v[168] * v[368];
v[440] = -(v[287] * v[407]) - v[289] * v[410] - v[291] * v[413];
v[428] = -(v[282] * v[407]) - v[283] * v[410] - v[284] * v[413];
v[227] = (v[172] * v[172]);
v[3588] = v[10] / v[227];
v[3337] = v[238] / v[227];
v[3336] = v[236] / v[227];
v[3335] = v[230] / v[227];
v[3711] = 2e0 / Power(v[227], 3);
v[175] = 1e0 - v[319] * v[3332];
v[3691] = v[175] / v[227];
v[1318] = v[175] * v[3335];
v[176] = v[172] * v[324];
v[3423] = v[176] * v[236];
v[1320] = v[176] * v[3336];
v[1794] = v[1318] + v[1320];
v[177] = v[172] * v[330];
v[3422] = v[177] * v[238];
v[1321] = v[177] * v[3337];
v[1799] = v[1318] + v[1321];
v[1788] = v[1320] + v[1799];
v[179] = v[172] * v[334];
v[1325] = v[179] * v[3335];
v[181] = 1e0 - v[3332] * v[338];
v[3690] = v[181] / v[227];
v[1314] = v[181] * v[3336];
v[1790] = v[1314] + v[1325];
v[182] = v[172] * v[343];
v[1315] = v[182] * v[3337];
v[1800] = v[1314] + v[1315];
v[1793] = v[1325] + v[1800];
v[184] = v[172] * v[347];
v[1328] = v[184] * v[3335];
v[186] = v[172] * v[351];
v[1310] = v[186] * v[3336];
v[187] = 1e0 - v[3332] * v[356];
v[3694] = v[187] / v[227];
v[1311] = v[187] * v[3337];
v[3692] = v[1311] * v[227];
v[1798] = v[1310] + v[1311] + v[1328];
v[1795] = -v[1328] + v[1798];
v[1789] = -v[1310] + v[1798];
v[245] = -(v[323] * v[329]);
v[3345] = v[245] - v[316];
v[243] = v[326] * v[329];
v[3343] = v[243] - v[317];
v[234] = -(v[323] * v[326]);
v[3340] = v[234] - v[318];
v[191] = v[175] * v[88] + v[176] * v[91] + v[177] * v[94];
v[192] = v[175] * v[89] + v[176] * v[92] + v[177] * v[95];
v[193] = v[175] * v[90] + v[176] * v[93] + v[177] * v[96];
v[306] = v[191] * v[301] + v[192] * v[303] + v[193] * v[305];
v[460] = -(v[132] * v[306]);
v[457] = -(v[130] * v[306]);
v[454] = -(v[128] * v[306]);
v[296] = v[191] * v[293] + v[192] * v[294] + v[193] * v[295];
v[448] = -(v[132] * v[296]);
v[445] = -(v[130] * v[296]);
v[442] = -(v[128] * v[296]);
v[194] = v[179] * v[88] + v[181] * v[91] + v[182] * v[94];
v[195] = v[179] * v[89] + v[181] * v[92] + v[182] * v[95];
v[196] = v[179] * v[90] + v[181] * v[93] + v[182] * v[96];
v[307] = v[194] * v[301] + v[195] * v[303] + v[196] * v[305];
v[461] = -(v[132] * v[307]);
v[458] = -(v[130] * v[307]);
v[455] = -(v[128] * v[307]);
v[297] = v[194] * v[293] + v[195] * v[294] + v[196] * v[295];
v[449] = -(v[132] * v[297]);
v[446] = -(v[130] * v[297]);
v[443] = -(v[128] * v[297]);
v[197] = v[184] * v[88] + v[186] * v[91] + v[187] * v[94];
v[198] = v[184] * v[89] + v[186] * v[92] + v[187] * v[95];
v[199] = v[184] * v[90] + v[186] * v[93] + v[187] * v[96];
v[308] = v[197] * v[301] + v[198] * v[303] + v[199] * v[305];
v[462] = -(v[132] * v[308]);
v[459] = -(v[130] * v[308]);
v[456] = -(v[128] * v[308]);
v[298] = v[197] * v[293] + v[198] * v[294] + v[199] * v[295];
v[450] = -(v[132] * v[298]);
v[447] = -(v[130] * v[298]);
v[444] = -(v[128] * v[298]);
v[200] = v[85] + v[97];
v[201] = v[86] + v[98];
v[202] = v[87] + v[99];
v[212] = dui1A[0] * v[11] + ddui1A[0] * v[12] + v[10] * v[45];
v[3351] = v[128] * v[212];
v[213] = dui1A[1] * v[11] + ddui1A[1] * v[12] + v[10] * v[46];
v[3348] = v[128] * v[213];
v[214] = dui1A[2] * v[11] + ddui1A[2] * v[12] + v[10] * v[47];
v[3354] = v[128] * v[214];
v[215] = dui2A[0] * v[11] + ddui2A[0] * v[12] + v[10] * v[48];
v[3350] = v[130] * v[215];
v[216] = dui2A[1] * v[11] + ddui2A[1] * v[12] + v[10] * v[49];
v[3347] = v[130] * v[216];
v[217] = dui2A[2] * v[11] + ddui2A[2] * v[12] + v[10] * v[50];
v[3353] = v[130] * v[217];
v[218] = dui3A[0] * v[11] + ddui3A[0] * v[12] + v[10] * v[51];
v[3349] = v[132] * v[218];
v[219] = dui3A[1] * v[11] + ddui3A[1] * v[12] + v[10] * v[52];
v[3346] = v[132] * v[219];
v[220] = dui3A[2] * v[11] + ddui3A[2] * v[12] + v[10] * v[53];
v[3352] = v[132] * v[220];
v[1637] = v[3352] + v[3353] + v[3354];
v[221] = duiB[0] * v[11] + dduiB[0] * v[12] + v[10] * v[97];
v[222] = duiB[1] * v[11] + dduiB[1] * v[12] + v[10] * v[98];
v[223] = duiB[2] * v[11] + dduiB[2] * v[12] + v[10] * v[99];
v[2542] = v[1637] - v[223];
v[225] = v[243] + v[317];
v[3344] = v[225] / v[227];
v[3338] = v[187] * v[225];
v[226] = v[234] + v[318];
v[3341] = v[226] / v[227];
v[3339] = v[181] * v[226];
v[228] = v[227] + (v[326] * v[326]);
v[3310] = -(v[238] * v[3338]) - v[236] * v[3339] - v[228] * (v[230] + v[3422] + v[3423]);
v[2507] = v[228] / v[227];
v[1440] = v[3343] / v[227];
v[1438] = v[3340] / v[227];
v[6124] = 0e0;
v[6125] = 0e0;
v[6126] = 0e0;
v[6127] = 0e0;
v[6128] = 0e0;
v[6129] = 0e0;
v[6130] = 0e0;
v[6131] = 0e0;
v[6132] = 0e0;
v[6133] = 0e0;
v[6134] = 0e0;
v[6135] = 0e0;
v[6136] = 0e0;
v[6137] = v[1438];
v[6138] = v[1440];
v[5786] = 0e0;
v[5787] = 0e0;
v[5788] = 0e0;
v[5789] = 0e0;
v[5790] = 0e0;
v[5791] = 0e0;
v[5792] = 0e0;
v[5793] = 0e0;
v[5794] = 0e0;
v[5795] = 0e0;
v[5796] = 0e0;
v[5797] = 0e0;
v[5798] = 0e0;
v[5799] = v[10] * v[1438];
v[5800] = v[10] * v[1440];
v[229] = v[1438] * v[236] + v[1440] * v[238] + v[230] * v[2507];
v[231] = v[227] + (v[323] * v[323]);
v[3411] = v[179] * v[231] + v[175] * v[3340];
v[237] = v[245] + v[316];
v[3342] = v[182] * v[231] + v[187] * v[237];
v[3309] = -(v[231] * v[236]) - v[238] * v[3342] - v[230] * v[3411];
v[2510] = v[231] / v[227];
v[1439] = v[3345] / v[227];
v[6139] = 0e0;
v[6140] = 0e0;
v[6141] = 0e0;
v[6142] = 0e0;
v[6143] = 0e0;
v[6144] = 0e0;
v[6145] = 0e0;
v[6146] = 0e0;
v[6147] = 0e0;
v[6148] = 0e0;
v[6149] = 0e0;
v[6150] = 0e0;
v[6151] = v[3341];
v[6152] = 0e0;
v[6153] = v[1439];
v[5816] = 0e0;
v[5817] = 0e0;
v[5818] = 0e0;
v[5819] = 0e0;
v[5820] = 0e0;
v[5821] = 0e0;
v[5822] = 0e0;
v[5823] = 0e0;
v[5824] = 0e0;
v[5825] = 0e0;
v[5826] = 0e0;
v[5827] = 0e0;
v[5828] = v[10] * v[3341];
v[5829] = 0e0;
v[5830] = v[10] * v[1439];
v[239] = v[1439] * v[238] + v[236] * v[2510] + v[230] * v[3341];
v[240] = v[227] + (v[329] * v[329]);
v[3412] = v[184] * v[240] + v[175] * v[3343];
v[3413] = v[186] * v[240] + v[181] * v[3345];
v[3308] = -(v[238] * v[240]) - v[230] * v[3412] - v[236] * v[3413];
v[2512] = v[240] / v[227];
v[1437] = v[237] / v[227];
v[6154] = 0e0;
v[6155] = 0e0;
v[6156] = 0e0;
v[6157] = 0e0;
v[6158] = 0e0;
v[6159] = 0e0;
v[6160] = 0e0;
v[6161] = 0e0;
v[6162] = 0e0;
v[6163] = 0e0;
v[6164] = 0e0;
v[6165] = 0e0;
v[6166] = v[3344];
v[6167] = v[1437];
v[6168] = 0e0;
v[5846] = 0e0;
v[5847] = 0e0;
v[5848] = 0e0;
v[5849] = 0e0;
v[5850] = 0e0;
v[5851] = 0e0;
v[5852] = 0e0;
v[5853] = 0e0;
v[5854] = 0e0;
v[5855] = 0e0;
v[5856] = 0e0;
v[5857] = 0e0;
v[5858] = v[10] * v[3344];
v[5859] = v[10] * v[1437];
v[5860] = 0e0;
v[249] = v[1437] * v[236] + v[238] * v[2512] + v[230] * v[3344];
v[1901] = -v[221] - v[229] * v[406] - v[239] * v[407] - v[249] * v[408];
v[1838] = v[1901] + v[3349] + v[3350] + v[3351];
v[2907] = v[1838] - v[1901];
v[2593] = -v[221] + v[2907];
v[1814] = -v[222] - v[229] * v[409] - v[239] * v[410] - v[249] * v[411];
v[1886] = v[1814] + v[3346] + v[3347] + v[3348];
v[2905] = -v[1814] + v[1886];
v[2569] = -v[222] + v[2905];
v[1804] = -v[223] - v[229] * v[412] - v[239] * v[413] - v[249] * v[414];
v[1938] = v[1637] + v[1804];
v[250] = v[128] * v[144] + v[130] * v[145] + v[132] * v[146] - v[166] * v[191] - v[167] * v[192] - v[168] * v[193] - v[200];
v[418] = -v[250] / 2e0;
v[436] = v[132] * v[287] - v[418];
v[430] = v[128] * v[287] + v[418];
v[419] = v[130] * v[282] - v[418];
v[415] = v[128] * v[282] + v[418];
v[251] = v[128] * v[148] + v[130] * v[149] + v[132] * v[150] - v[166] * v[194] - v[167] * v[195] - v[168] * v[196] - v[201];
v[420] = -v[251] / 2e0;
v[437] = v[132] * v[289] - v[420];
v[431] = v[128] * v[289] + v[420];
v[421] = v[130] * v[283] - v[420];
v[416] = v[128] * v[283] + v[420];
v[252] = v[128] * v[152] + v[130] * v[153] + v[132] * v[154] - v[166] * v[197] - v[167] * v[198] - v[168] * v[199] - v[202];
v[602] = sqrt((v[250] * v[250]) + (v[251] * v[251]) + (v[252] * v[252]));
v[820] = 1e0 / (v[602] * v[602]);
v[465] = -(v[250] * (v[301] * v[363] + v[303] * v[366] + v[305] * v[369])) - v[251] * (v[301] * v[378] + v[303] * v[381]
	+ v[305] * v[384]) - v[252] * (v[301] * v[393] + v[303] * v[396] + v[305] * v[399]) + v[306] * v[408] + v[307] * v[411]
	+ v[308] * v[414];
v[464] = -(v[250] * (v[301] * v[362] + v[303] * v[365] + v[305] * v[368])) - v[251] * (v[301] * v[377] + v[303] * v[380]
	+ v[305] * v[383]) - v[252] * (v[301] * v[392] + v[303] * v[395] + v[305] * v[398]) + v[306] * v[407] + v[307] * v[410]
	+ v[308] * v[413];
v[463] = -(v[250] * (v[301] * v[361] + v[303] * v[364] + v[305] * v[367])) - v[251] * (v[301] * v[376] + v[303] * v[379]
	+ v[305] * v[382]) - v[252] * (v[301] * v[391] + v[303] * v[394] + v[305] * v[397]) + v[306] * v[406] + v[307] * v[409]
	+ v[308] * v[412];
v[453] = -(v[250] * (v[293] * v[363] + v[294] * v[366] + v[295] * v[369])) - v[251] * (v[293] * v[378] + v[294] * v[381]
	+ v[295] * v[384]) - v[252] * (v[293] * v[393] + v[294] * v[396] + v[295] * v[399]) + v[296] * v[408] + v[297] * v[411]
	+ v[298] * v[414];
v[452] = -(v[250] * (v[293] * v[362] + v[294] * v[365] + v[295] * v[368])) - v[251] * (v[293] * v[377] + v[294] * v[380]
	+ v[295] * v[383]) - v[252] * (v[293] * v[392] + v[294] * v[395] + v[295] * v[398]) + v[296] * v[407] + v[297] * v[410]
	+ v[298] * v[413];
v[451] = -(v[250] * (v[293] * v[361] + v[294] * v[364] + v[295] * v[367])) - v[251] * (v[293] * v[376] + v[294] * v[379]
	+ v[295] * v[382]) - v[252] * (v[293] * v[391] + v[294] * v[394] + v[295] * v[397]) + v[296] * v[406] + v[297] * v[409]
	+ v[298] * v[412];
v[422] = -v[252] / 2e0;
v[438] = v[132] * v[291] - v[422];
v[432] = v[128] * v[291] + v[422];
v[423] = v[130] * v[284] - v[422];
v[417] = v[128] * v[284] + v[422];
v[253] = v[133] * v[36] + v[135] * v[39] + v[137] * v[42] - v[85] - v[169] * v[88] - v[170] * v[89] - v[171] * v[90];
v[254] = v[133] * v[37] + v[135] * v[40] + v[137] * v[43] - v[86] - v[169] * v[91] - v[170] * v[92] - v[171] * v[93];
v[255] = v[133] * v[38] + v[135] * v[41] + v[137] * v[44] - v[87] - v[169] * v[94] - v[170] * v[95] - v[171] * v[96];
if (v[602] > 0.1e-7) { v01 = 1e0 / v[602]; v02 = (-(v01 / v[602])); v03 = (2e0*v01) / (v[602] * v[602]); }
else {
	v01 = (12500000e0 / 3e0)*(24e0 - (-0.1e-7 + v[602])*(0.24e10 - 2e0*(-1e0 + 100000000e0*v[602])*(0.2399999997e10
		- 0.1199999994e18*v[602] - 0.3e17*(v[602] * v[602]))));
	v02 = (-50000000e0 / 3e0)*(0.3599999994e10 - 0.4799999982e18*v[602] + 0.6e25*Power(v[602], 3)
		+ 0.1799999982e26*(v[602] * v[602]));
	v03 = 0.1e17*(799999997e0 - 0.599999994e17*v[602] - 0.3e17*(v[602] * v[602]));
};
v[262] = v03;
v[263] = v02;
v[264] = v01;
v[265] = v[250] * v[264];
v[266] = v[251] * v[264];
v[267] = v[252] * v[264];
v[268] = sqrt((v[253] * v[253]) + (v[254] * v[254]) + (v[255] * v[255]));
if (v[268] > 0.1e-7) { v04 = 1e0 / v[268]; v05 = (-(v04 / v[268])); v06 = (2e0*v04) / (v[268] * v[268]); }
else {
	v04 = (12500000e0 / 3e0)*(24e0 - (-0.1e-7 + v[268])*(0.24e10 - 2e0*(-1e0 + 100000000e0*v[268])*(0.2399999997e10
		- 0.1199999994e18*v[268] - 0.3e17*(v[268] * v[268]))));
	v05 = (-50000000e0 / 3e0)*(0.3599999994e10 - 0.4799999982e18*v[268] + 0.6e25*Power(v[268], 3)
		+ 0.1799999982e26*(v[268] * v[268]));
	v06 = 0.1e17*(799999997e0 - 0.599999994e17*v[268] - 0.3e17*(v[268] * v[268]));
};
v[273] = v04;
v[274] = v[253] * v[273];
v[275] = v[254] * v[273];
v[276] = v[255] * v[273];
b277 = v[265] * v[274] + v[266] * v[275] + v[267] * v[276] < 0e0;
if (b277) {
	v[279] = -v[265];
	v[280] = -v[266];
	v[281] = -v[267];
}
else {
	v[279] = v[265];
	v[280] = v[266];
	v[281] = v[267];
};
v[3714] = v[128] * v[281];
v[3713] = v[130] * v[281];
v[3712] = v[132] * v[281];
v[6169] = 0e0;
v[6170] = 0e0;
v[6171] = v[3714];
v[6172] = 0e0;
v[6173] = 0e0;
v[6174] = v[3713];
v[6175] = 0e0;
v[6176] = 0e0;
v[6177] = v[3712];
v[6178] = 0e0;
v[6179] = 0e0;
v[6180] = 0e0;
v[6181] = 0e0;
v[6182] = 0e0;
v[6183] = 0e0;
v[3466] = 2e0*v[281];
v[1167] = -(v[281] * v[414]);
v[1165] = -(v[281] * v[413]);
v[1163] = -(v[281] * v[412]);
v[2179] = v[281] * (-v[222] + v[3346] + v[3347] + v[3348]);
v[1952] = v[2179] + (v[1814] + v[222])*v[281];
v[2178] = v[281] * (-v[221] + v[3349] + v[3350] + v[3351]);
v[1963] = v[2178] + (v[1901] + v[221])*v[281];
v[592] = (v[281] * v[281]);
v[1160] = v[281] * (-v[223] + v[3352] + v[3353] + v[3354]);
v[3469] = 2e0*v[280];
v[1900] = v[1901] * v[280];
v[1858] = v[1804] * v[280];
v[1166] = -(v[280] * v[411]);
v[1832] = v[1166] + v[1167];
v[1164] = -(v[280] * v[410]);
v[1834] = v[1164] + v[1165];
v[1162] = -(v[280] * v[409]);
v[1836] = v[1162] + v[1163];
v[1161] = -(v[222] * v[280]);
v[2598] = v[1160] + v[1161];
v[1823] = v[1836] * v[229] + v[1834] * v[239] + v[1832] * v[249] + v[2598];
v[587] = (v[280] * v[280]);
v[3814] = v[215] * v[279] + v[216] * v[280];
v[3813] = v[218] * v[279] + v[219] * v[280];
v[3805] = v[212] * v[279] + v[213] * v[280];
v[3472] = 2e0*v[279];
v[3397] = v[279] * v[281];
v[3361] = v[279] * v[280];
v[3359] = -(v[279] * v[406]);
v[3360] = v[1163] + v[3359];
v[3357] = -(v[279] * v[407]);
v[3358] = v[1165] + v[3357];
v[3355] = -(v[279] * v[408]);
v[3356] = v[1167] + v[3355];
v[1813] = v[1814] * v[279];
v[1802] = v[1804] * v[279];
v[1133] = v[280] * v[3360] - v[409] * v[587];
v[1131] = v[280] * v[3358] - v[410] * v[587];
v[1129] = v[280] * v[3356] - v[411] * v[587];
v[3846] = v[1832] + 2e0*v[3355];
v[1923] = v[1166] + v[3355];
v[3824] = 2e0*v[1167] + v[1923];
v[3840] = 2e0*v[1166] + v[3356];
v[3845] = v[1834] + 2e0*v[3357];
v[1925] = v[1164] + v[3357];
v[3823] = 2e0*v[1165] + v[1925];
v[3839] = 2e0*v[1164] + v[3358];
v[3844] = v[1836] + 2e0*v[3359];
v[1927] = v[1162] + v[3359];
v[3822] = 2e0*v[1163] + v[1927];
v[3838] = 2e0*v[1162] + v[3360];
v[1122] = -(v[221] * v[279]);
v[2571] = v[1122] + v[1160];
v[2550] = v[1122] + v[1161];
v[1909] = v[1927] * v[229] + v[1925] * v[239] + v[1923] * v[249] + v[2550] + v[280] * v[2905] + v[279] * v[2907];
v[1871] = v[2571] + v[249] * v[3356] + v[239] * v[3358] + v[229] * v[3360];
v[1090] = v[1927] * v[281] - v[412] * v[592];
v[1088] = v[1925] * v[281] - v[413] * v[592];
v[1086] = v[1923] * v[281] - v[414] * v[592];
v[586] = v[132] * v[3361];
v[585] = v[130] * v[3361];
v[584] = v[128] * v[3361];
v[582] = (v[279] * v[279]);
v[1175] = v[1836] * v[279] - v[406] * v[582];
v[1173] = v[1834] * v[279] - v[407] * v[582];
v[1171] = v[1832] * v[279] - v[408] * v[582];
v[482] = -(v[13] * v[415]) - v[14] * v[430] - v[15] * v[442] - v[16] * v[454];
v[483] = -(v[13] * v[416]) - v[14] * v[431] - v[15] * v[443] - v[16] * v[455];
v[484] = -(v[13] * v[417]) - v[14] * v[432] - v[15] * v[444] - v[16] * v[456];
v[485] = -(v[13] * v[419]) - v[14] * v[433] - v[15] * v[445] - v[16] * v[457];
v[486] = -(v[13] * v[421]) - v[14] * v[434] - v[15] * v[446] - v[16] * v[458];
v[487] = -(v[13] * v[423]) - v[14] * v[435] - v[15] * v[447] - v[16] * v[459];
v[488] = -(v[13] * v[424]) - v[14] * v[436] - v[15] * v[448] - v[16] * v[460];
v[489] = -(v[13] * v[425]) - v[14] * v[437] - v[15] * v[449] - v[16] * v[461];
v[490] = -(v[13] * v[426]) - v[14] * v[438] - v[15] * v[450] - v[16] * v[462];
v[491] = v[13] * v[282] + v[14] * v[287] - v[15] * v[296] - v[16] * v[306];
v[492] = v[13] * v[283] + v[14] * v[289] - v[15] * v[297] - v[16] * v[307];
v[493] = v[13] * v[284] + v[14] * v[291] - v[15] * v[298] - v[16] * v[308];
v[494] = -(v[13] * v[427]) - v[14] * v[439] - v[15] * v[451] - v[16] * v[463];
v[495] = -(v[13] * v[428]) - v[14] * v[440] - v[15] * v[452] - v[16] * v[464];
v[496] = -(v[13] * v[429]) - v[14] * v[441] - v[15] * v[453] - v[16] * v[465];
v[1025] = -(v[212] * v[482]) - v[213] * v[483] - v[214] * v[484] - v[215] * v[485] - v[216] * v[486] - v[217] * v[487]
- v[218] * v[488] - v[219] * v[489] - v[220] * v[490] - v[221] * v[491] - v[222] * v[492] - v[223] * v[493] - v[229] * v[494]
- v[239] * v[495] - v[249] * v[496];
v[4919] = v[482];
v[4920] = v[483];
v[4921] = v[484];
v[4922] = v[485];
v[4923] = v[486];
v[4924] = v[487];
v[4925] = v[488];
v[4926] = v[489];
v[4927] = v[490];
v[4928] = v[491];
v[4929] = v[492];
v[4930] = v[493];
v[4931] = v[494];
v[4932] = v[495];
v[4933] = v[496];
v[497] = -(v[17] * v[415]) - v[18] * v[430] - v[19] * v[442] - v[20] * v[454];
v[498] = -(v[17] * v[416]) - v[18] * v[431] - v[19] * v[443] - v[20] * v[455];
v[499] = -(v[17] * v[417]) - v[18] * v[432] - v[19] * v[444] - v[20] * v[456];
v[500] = -(v[17] * v[419]) - v[18] * v[433] - v[19] * v[445] - v[20] * v[457];
v[501] = -(v[17] * v[421]) - v[18] * v[434] - v[19] * v[446] - v[20] * v[458];
v[502] = -(v[17] * v[423]) - v[18] * v[435] - v[19] * v[447] - v[20] * v[459];
v[503] = -(v[17] * v[424]) - v[18] * v[436] - v[19] * v[448] - v[20] * v[460];
v[504] = -(v[17] * v[425]) - v[18] * v[437] - v[19] * v[449] - v[20] * v[461];
v[505] = -(v[17] * v[426]) - v[18] * v[438] - v[19] * v[450] - v[20] * v[462];
v[506] = v[17] * v[282] + v[18] * v[287] - v[19] * v[296] - v[20] * v[306];
v[507] = v[17] * v[283] + v[18] * v[289] - v[19] * v[297] - v[20] * v[307];
v[508] = v[17] * v[284] + v[18] * v[291] - v[19] * v[298] - v[20] * v[308];
v[509] = -(v[17] * v[427]) - v[18] * v[439] - v[19] * v[451] - v[20] * v[463];
v[510] = -(v[17] * v[428]) - v[18] * v[440] - v[19] * v[452] - v[20] * v[464];
v[511] = -(v[17] * v[429]) - v[18] * v[441] - v[19] * v[453] - v[20] * v[465];
v[1023] = -(v[212] * v[497]) - v[213] * v[498] - v[214] * v[499] - v[215] * v[500] - v[216] * v[501] - v[217] * v[502]
- v[218] * v[503] - v[219] * v[504] - v[220] * v[505] - v[221] * v[506] - v[222] * v[507] - v[223] * v[508] - v[229] * v[509]
- v[239] * v[510] - v[249] * v[511];
v[4934] = v[497];
v[4935] = v[498];
v[4936] = v[499];
v[4937] = v[500];
v[4938] = v[501];
v[4939] = v[502];
v[4940] = v[503];
v[4941] = v[504];
v[4942] = v[505];
v[4943] = v[506];
v[4944] = v[507];
v[4945] = v[508];
v[4946] = v[509];
v[4947] = v[510];
v[4948] = v[511];
v[512] = -(v[21] * v[415]) - v[22] * v[430] - v[23] * v[442] - v[24] * v[454];
v[513] = -(v[21] * v[416]) - v[22] * v[431] - v[23] * v[443] - v[24] * v[455];
v[514] = -(v[21] * v[417]) - v[22] * v[432] - v[23] * v[444] - v[24] * v[456];
v[515] = -(v[21] * v[419]) - v[22] * v[433] - v[23] * v[445] - v[24] * v[457];
v[516] = -(v[21] * v[421]) - v[22] * v[434] - v[23] * v[446] - v[24] * v[458];
v[517] = -(v[21] * v[423]) - v[22] * v[435] - v[23] * v[447] - v[24] * v[459];
v[518] = -(v[21] * v[424]) - v[22] * v[436] - v[23] * v[448] - v[24] * v[460];
v[519] = -(v[21] * v[425]) - v[22] * v[437] - v[23] * v[449] - v[24] * v[461];
v[520] = -(v[21] * v[426]) - v[22] * v[438] - v[23] * v[450] - v[24] * v[462];
v[521] = v[21] * v[282] + v[22] * v[287] - v[23] * v[296] - v[24] * v[306];
v[522] = v[21] * v[283] + v[22] * v[289] - v[23] * v[297] - v[24] * v[307];
v[523] = v[21] * v[284] + v[22] * v[291] - v[23] * v[298] - v[24] * v[308];
v[524] = -(v[21] * v[427]) - v[22] * v[439] - v[23] * v[451] - v[24] * v[463];
v[525] = -(v[21] * v[428]) - v[22] * v[440] - v[23] * v[452] - v[24] * v[464];
v[526] = -(v[21] * v[429]) - v[22] * v[441] - v[23] * v[453] - v[24] * v[465];
v[1019] = v[212] * v[512] + v[213] * v[513] + v[214] * v[514] + v[215] * v[515] + v[216] * v[516] + v[217] * v[517]
+ v[218] * v[518] + v[219] * v[519] + v[220] * v[520] + v[221] * v[521] + v[222] * v[522] + v[223] * v[523] + v[229] * v[524]
+ v[239] * v[525] + v[249] * v[526];
v[4949] = v[512];
v[4950] = v[513];
v[4951] = v[514];
v[4952] = v[515];
v[4953] = v[516];
v[4954] = v[517];
v[4955] = v[518];
v[4956] = v[519];
v[4957] = v[520];
v[4958] = v[521];
v[4959] = v[522];
v[4960] = v[523];
v[4961] = v[524];
v[4962] = v[525];
v[4963] = v[526];
v[527] = -(v[25] * v[415]) - v[26] * v[430] - v[27] * v[442] - v[28] * v[454];
v[912] = v[284] * v[482] + v[291] * v[497] - v[298] * v[512] - v[308] * v[527];
v[911] = v[283] * v[482] + v[289] * v[497] - v[297] * v[512] - v[307] * v[527];
v[910] = v[282] * v[482] + v[287] * v[497] - v[296] * v[512] - v[306] * v[527];
v[528] = -(v[25] * v[416]) - v[26] * v[431] - v[27] * v[443] - v[28] * v[455];
v[916] = v[284] * v[483] + v[291] * v[498] - v[298] * v[513] - v[308] * v[528];
v[915] = v[283] * v[483] + v[289] * v[498] - v[297] * v[513] - v[307] * v[528];
v[914] = v[282] * v[483] + v[287] * v[498] - v[296] * v[513] - v[306] * v[528];
v[529] = -(v[25] * v[417]) - v[26] * v[432] - v[27] * v[444] - v[28] * v[456];
v[920] = v[284] * v[484] + v[291] * v[499] - v[298] * v[514] - v[308] * v[529];
v[919] = v[283] * v[484] + v[289] * v[499] - v[297] * v[514] - v[307] * v[529];
v[918] = v[282] * v[484] + v[287] * v[499] - v[296] * v[514] - v[306] * v[529];
v[530] = -(v[25] * v[419]) - v[26] * v[433] - v[27] * v[445] - v[28] * v[457];
v[924] = v[284] * v[485] + v[291] * v[500] - v[298] * v[515] - v[308] * v[530];
v[923] = v[283] * v[485] + v[289] * v[500] - v[297] * v[515] - v[307] * v[530];
v[922] = v[282] * v[485] + v[287] * v[500] - v[296] * v[515] - v[306] * v[530];
v[531] = -(v[25] * v[421]) - v[26] * v[434] - v[27] * v[446] - v[28] * v[458];
v[928] = v[284] * v[486] + v[291] * v[501] - v[298] * v[516] - v[308] * v[531];
v[927] = v[283] * v[486] + v[289] * v[501] - v[297] * v[516] - v[307] * v[531];
v[926] = v[282] * v[486] + v[287] * v[501] - v[296] * v[516] - v[306] * v[531];
v[532] = -(v[25] * v[423]) - v[26] * v[435] - v[27] * v[447] - v[28] * v[459];
v[932] = v[284] * v[487] + v[291] * v[502] - v[298] * v[517] - v[308] * v[532];
v[931] = v[283] * v[487] + v[289] * v[502] - v[297] * v[517] - v[307] * v[532];
v[930] = v[282] * v[487] + v[287] * v[502] - v[296] * v[517] - v[306] * v[532];
v[533] = -(v[25] * v[424]) - v[26] * v[436] - v[27] * v[448] - v[28] * v[460];
v[936] = v[284] * v[488] + v[291] * v[503] - v[298] * v[518] - v[308] * v[533];
v[935] = v[283] * v[488] + v[289] * v[503] - v[297] * v[518] - v[307] * v[533];
v[934] = v[282] * v[488] + v[287] * v[503] - v[296] * v[518] - v[306] * v[533];
v[534] = -(v[25] * v[425]) - v[26] * v[437] - v[27] * v[449] - v[28] * v[461];
v[940] = v[284] * v[489] + v[291] * v[504] - v[298] * v[519] - v[308] * v[534];
v[939] = v[283] * v[489] + v[289] * v[504] - v[297] * v[519] - v[307] * v[534];
v[938] = v[282] * v[489] + v[287] * v[504] - v[296] * v[519] - v[306] * v[534];
v[535] = -(v[25] * v[426]) - v[26] * v[438] - v[27] * v[450] - v[28] * v[462];
v[944] = v[284] * v[490] + v[291] * v[505] - v[298] * v[520] - v[308] * v[535];
v[943] = v[283] * v[490] + v[289] * v[505] - v[297] * v[520] - v[307] * v[535];
v[942] = v[282] * v[490] + v[287] * v[505] - v[296] * v[520] - v[306] * v[535];
v[536] = v[25] * v[282] + v[26] * v[287] - v[27] * v[296] - v[28] * v[306];
v[948] = v[284] * v[491] + v[291] * v[506] - v[298] * v[521] - v[308] * v[536];
v[947] = v[283] * v[491] + v[289] * v[506] - v[297] * v[521] - v[307] * v[536];
v[946] = v[282] * v[491] + v[287] * v[506] - v[296] * v[521] - v[306] * v[536];
v[537] = v[25] * v[283] + v[26] * v[289] - v[27] * v[297] - v[28] * v[307];
v[952] = v[284] * v[492] + v[291] * v[507] - v[298] * v[522] - v[308] * v[537];
v[951] = v[283] * v[492] + v[289] * v[507] - v[297] * v[522] - v[307] * v[537];
v[950] = v[282] * v[492] + v[287] * v[507] - v[296] * v[522] - v[306] * v[537];
v[538] = v[25] * v[284] + v[26] * v[291] - v[27] * v[298] - v[28] * v[308];
v[956] = v[284] * v[493] + v[291] * v[508] - v[298] * v[523] - v[308] * v[538];
v[955] = v[283] * v[493] + v[289] * v[508] - v[297] * v[523] - v[307] * v[538];
v[954] = v[282] * v[493] + v[287] * v[508] - v[296] * v[523] - v[306] * v[538];
v[539] = -(v[25] * v[427]) - v[26] * v[439] - v[27] * v[451] - v[28] * v[463];
v[960] = v[284] * v[494] + v[291] * v[509] - v[298] * v[524] - v[308] * v[539];
v[959] = v[283] * v[494] + v[289] * v[509] - v[297] * v[524] - v[307] * v[539];
v[958] = v[282] * v[494] + v[287] * v[509] - v[296] * v[524] - v[306] * v[539];
v[540] = -(v[25] * v[428]) - v[26] * v[440] - v[27] * v[452] - v[28] * v[464];
v[964] = v[284] * v[495] + v[291] * v[510] - v[298] * v[525] - v[308] * v[540];
v[963] = v[283] * v[495] + v[289] * v[510] - v[297] * v[525] - v[307] * v[540];
v[962] = v[282] * v[495] + v[287] * v[510] - v[296] * v[525] - v[306] * v[540];
v[541] = -(v[25] * v[429]) - v[26] * v[441] - v[27] * v[453] - v[28] * v[465];
v[1021] = v[212] * v[527] + v[213] * v[528] + v[214] * v[529] + v[215] * v[530] + v[216] * v[531] + v[217] * v[532]
+ v[218] * v[533] + v[219] * v[534] + v[220] * v[535] + v[221] * v[536] + v[222] * v[537] + v[223] * v[538] + v[229] * v[539]
+ v[239] * v[540] + v[249] * v[541];
v[968] = v[284] * v[496] + v[291] * v[511] - v[298] * v[526] - v[308] * v[541];
v[967] = v[283] * v[496] + v[289] * v[511] - v[297] * v[526] - v[307] * v[541];
v[966] = v[282] * v[496] + v[287] * v[511] - v[296] * v[526] - v[306] * v[541];
v[4964] = v[527];
v[4965] = v[528];
v[4966] = v[529];
v[4967] = v[530];
v[4968] = v[531];
v[4969] = v[532];
v[4970] = v[533];
v[4971] = v[534];
v[4972] = v[535];
v[4973] = v[536];
v[4974] = v[537];
v[4975] = v[538];
v[4976] = v[539];
v[4977] = v[540];
v[4978] = v[541];
b542 = sqrt(Power(-(v[275] * v[279]) + v[274] * v[280], 2) + Power(v[276] * v[279] - v[274] * v[281], 2) + Power(-
(v[276] * v[280]) + v[275] * v[281], 2)) > 0.1e-7;
if (b542) {
	v[544] = -(v[276] * v[280]) + v[275] * v[281];
	v[545] = v[276] * v[279] - v[274] * v[281];
	v[546] = -(v[275] * v[279]) + v[274] * v[280];
	v[547] = sqrt((v[544] * v[544]) + (v[545] * v[545]) + (v[546] * v[546]));
	v[1703] = 1e0 / (v[547] * v[547]);
	v[1287] = v[547];
	v[1714] = 1e0 - (v[1287] * v[1287]);
	v[3770] = 1e0 / Power(v[1714], 0.15e1);
	v[1709] = 1e0 / sqrt(v[1714]);
	v[1286] = asin(v[1287]) / 2e0;
	v[1708] = 1e0 / Power(cos(v[1286]), 2);
	v[3473] = v[1708] * v[1709];
	v[549] = 2e0*tan(v[1286]);
	if (v[547] > 0.1e-7) { v07 = 1e0 / v[547]; v08 = (-(v07 / v[547])); v09 = (2e0*v07) / (v[547] * v[547]); }
	else {
		v07 = (12500000e0 / 3e0)*(24e0 - (-0.1e-7 + v[547])*(0.24e10 - 2e0*(-1e0 + 100000000e0*v[547])*
			(0.2399999997e10 - 0.1199999994e18*v[547] - 0.3e17*(v[547] * v[547]))));
		v08 = (-50000000e0 / 3e0)*(0.3599999994e10 - 0.4799999982e18*v[547] + 0.6e25*Power(v[547], 3)
			+ 0.1799999982e26*(v[547] * v[547]));
		v09 = 0.1e17*(799999997e0 - 0.599999994e17*v[547] - 0.3e17*(v[547] * v[547]));
	};
	v[553] = v09;
	v[554] = v08;
	v[555] = v07;
	v[3769] = v[549] * v[554] + v[3473] * v[555];
	v[3362] = v[549] * v[555];
	v[556] = v[3362] * v[544];
	v[3529] = 2e0*v[556];
	v[3404] = v[556] / 2e0;
	v[567] = (v[556] * v[556]);
	v[557] = v[3362] * v[545];
	v[3363] = v[557] / 2e0;
	v[565] = v[3363] * v[556];
	v[560] = (v[557] * v[557]);
	v[1241] = -v[560] - v[567];
	v[558] = v[3362] * v[546];
	v[3526] = 2e0*v[558];
	v[1271] = -v[558] + v[565];
	v[1262] = v[558] + v[565];
	v[572] = v[3363] * v[558];
	v[1253] = -v[556] + v[572];
	v[1245] = v[556] + v[572];
	v[570] = v[3404] * v[558];
	v[1266] = v[557] + v[570];
	v[1249] = -v[557] + v[570];
	v[561] = (v[558] * v[558]);
	v[1281] = 4e0 + v[560] + v[561] + v[567];
	v[3771] = 1e0 / Power(v[1281], 3);
	v[3525] = -4e0 / (v[1281] * v[1281]);
	v[1276] = -v[560] - v[561];
	v[1258] = -v[561] - v[567];
	v[559] = 4e0 / v[1281];
	v[3364] = v[559] / 2e0;
	v[562] = 1e0 + v[1276] * v[3364];
	v[563] = v[1271] * v[559];
	v[564] = v[1266] * v[559];
	v[566] = v[1262] * v[559];
	v[568] = 1e0 + v[1258] * v[3364];
	v[569] = v[1253] * v[559];
	v[571] = v[1249] * v[559];
	v[573] = v[1245] * v[559];
	v[574] = 1e0 + v[1241] * v[3364];
}
else {
	v[562] = 1e0;
	v[563] = 0e0;
	v[564] = 0e0;
	v[566] = 0e0;
	v[568] = 1e0;
	v[569] = 0e0;
	v[571] = 0e0;
	v[573] = 0e0;
	v[574] = 1e0;
};
if (b30) {
	v[1214] = 1e0 - v[592];
	v[1212] = 1e0 - v[587];
	v[1210] = 1e0 - v[582];
	v[579] = v[133] * v[152] + v[135] * v[153] + v[137] * v[154] - v[169] * v[197] - v[170] * v[198] - v[171] * v[199] - v[202]
		+ v[1] * v[571] + v[2] * v[573] + v[3] * v[574];
	v[3365] = v[281] * v[579];
	v[578] = v[133] * v[148] + v[135] * v[149] + v[137] * v[150] - v[169] * v[194] - v[170] * v[195] - v[171] * v[196] - v[201]
		+ v[1] * v[566] + v[2] * v[568] + v[3] * v[569];
	v[3367] = v[280] * v[578];
	v[3401] = v[3365] + v[3367];
	v[577] = v[133] * v[144] + v[135] * v[145] + v[137] * v[146] - v[169] * v[191] - v[170] * v[192] - v[171] * v[193] - v[200]
		+ v[1] * v[562] + v[2] * v[563] + v[3] * v[564];
	v[3366] = -(v[279] * v[577]);
	v[3403] = v[3366] - v[3367];
	v[3402] = -v[3365] + v[3366];
	v[576] = -(v[279] * v[3401]) + v[1210] * v[577];
	v[580] = v[280] * v[3402] + v[1212] * v[578];
	v[581] = v[281] * v[3403] + v[1214] * v[579];
}
else {
	v[576] = 0e0;
	v[580] = 0e0;
	v[581] = 0e0;
};
v[583] = v[1175] * v[229] + v[1173] * v[239] + v[1171] * v[249] + v[2598] * v[279] + v[2593] * v[582] + v[213] * v[584]
+ v[216] * v[585] + v[219] * v[586];
v[591] = v[1133] * v[229] + v[1131] * v[239] + v[1129] * v[249] + v[2571] * v[280] + v[212] * v[584] + v[215] * v[585]
+ v[218] * v[586] + v[2569] * v[587];
v[593] = v[1090] * v[229] + v[1088] * v[239] + v[1086] * v[249] + v[2178] * v[279] + v[2179] * v[280] + v[2542] * v[592];
(*vnrel) = sqrt((v[583] * v[583]) + (v[591] * v[591]) + (v[593] * v[593]));
v[599] = -((v[3369] * Power(v[600], v[615])) / (v[35] * Power(v[33], v[620])));
v[3382] = -(v[35] * v[599]);
b603 = v[602] < v[32];
if (b603) {
	b605 = v[602] > v[33];
	if (b605) {
		v[614] = v[32] - v[602];
		v[3368] = -(v[31] * Power(v[614], v[34]));
		v[607] = v[279] * v[3368];
		v[609] = v[280] * v[3368];
		v[610] = v[281] * v[3368];
	}
	else {
		v[611] = -(v[31] * Power(v[600], v[34])) + v[599] * (Power(v[33], v[35]) - Power(v[602], v[35]));
		v[607] = v[279] * v[611];
		v[609] = v[280] * v[611];
		v[610] = v[281] * v[611];
	};
}
else {
	v[607] = 0e0;
	v[609] = 0e0;
	v[610] = 0e0;
};
if (b603) {
	v[3370] = 2e0*v[7];
	if (b605) {
		v[1064] = v[3369] * v[8];
		v[1065] = sqrt(v[1064] * Power(v[614], v[615]));
		v[617] = v[1065] * v[3370];
		v[616] = v[583] * v[617];
		v[618] = v[591] * v[617];
		v[619] = v[593] * v[617];
	}
	else {
		v[1072] = v[3382] * v[8];
		v[1073] = sqrt(v[1072] * Power(v[602], v[620]));
		v[621] = v[1073] * v[3370];
		v[616] = v[583] * v[621];
		v[618] = v[591] * v[621];
		v[619] = v[593] * v[621];
	};
	b622 = v[602] < v[32] && v[279] * (v[607] + v[616]) + v[280] * (v[609] + v[618]) + v[281] * (v[610] + v[619]) < 0e0;
	if (b622) {
		v[624] = v[616];
		v[625] = v[618];
		v[626] = v[619];
	}
	else {
		v[624] = -v[607];
		v[625] = -v[609];
		v[626] = -v[610];
	};
}
else {
	v[624] = 0e0;
	v[625] = 0e0;
	v[626] = 0e0;
};
v[627] = v[607] + v[624];
v[628] = v[609] + v[625];
v[629] = v[610] + v[626];
v[2029] = (v[627] * v[627]) + (v[628] * v[628]) + (v[629] * v[629]);
v[630] = v[4] * v[576];
v[631] = v[4] * v[580];
v[632] = v[4] * v[581];
v[636] = v[630] - v[9] * (v[212] * v[910] + v[213] * v[914] + v[214] * v[918] + v[215] * v[922] + v[216] * v[926]
	+ v[217] * v[930] + v[218] * v[934] + v[219] * v[938] + v[220] * v[942] + v[221] * v[946] + v[222] * v[950] + v[223] * v[954]
	+ v[229] * v[958] + v[239] * v[962] + v[249] * v[966]);
v[637] = v[631] - v[9] * (v[212] * v[911] + v[213] * v[915] + v[214] * v[919] + v[215] * v[923] + v[216] * v[927]
	+ v[217] * v[931] + v[218] * v[935] + v[219] * v[939] + v[220] * v[943] + v[221] * v[947] + v[222] * v[951] + v[223] * v[955]
	+ v[229] * v[959] + v[239] * v[963] + v[249] * v[967]);
v[638] = v[632] - v[9] * (v[212] * v[912] + v[213] * v[916] + v[214] * v[920] + v[215] * v[924] + v[216] * v[928]
	+ v[217] * v[932] + v[218] * v[936] + v[219] * v[940] + v[220] * v[944] + v[221] * v[948] + v[222] * v[952] + v[223] * v[956]
	+ v[229] * v[960] + v[239] * v[964] + v[249] * v[968]);
v[2025] = (v[636] * v[636]) + (v[637] * v[637]) + (v[638] * v[638]);
if (b603) {
	if (b29) {
		b641 = sqrt((v[636] * v[636]) + (v[637] * v[637]) + (v[638] * v[638])) <= v[5] * sqrt((v[627] * v[627]) +
			(v[628] * v[628]) + (v[629] * v[629]));
		if (b641) {
			v[643] = v[636];
			v[644] = v[637];
			v[645] = v[638];
			v[646] = 1e0;
		}
		else {
			v[3371] = v[6] * sqrt(v[2029]);
			v[647] = sqrt(v[2025]);
			if (v[647] > 0.1e-5) { v010 = 1e0 / v[647]; v011 = (-(v010 / v[647])); v012 = (2e0*v010) / (v[647] * v[647]); }
			else {
				v010 = (24000000e0 - (-1e0 + 1000000e0*v[647])*(71999994e0 - 0.71999982e14*v[647] + 0.6e19*Power(v[647], 3)
					+ 0.23999982e20*(v[647] * v[647]))) / 24e0;
				v011 = (-500000e0 / 3e0)*(35999994e0 - 0.47999982e14*v[647] + 0.6e19*Power(v[647], 3) + 0.17999982e20*
					(v[647] * v[647]));
				v012 = 0.1e13*(7999997e0 - 0.5999994e13*v[647] - 0.3e13*(v[647] * v[647]));
			};
			v[651] = v011;
			v[652] = v010;
			v[653] = v[636] * v[652];
			v[654] = v[637] * v[652];
			v[655] = v[638] * v[652];
			v[643] = v[3371] * v[653];
			v[644] = v[3371] * v[654];
			v[645] = v[3371] * v[655];
			v[646] = 0e0;
		};
		if (sqrt((v[630] * v[630]) + (v[631] * v[631]) + (v[632] * v[632])) > v[5] * sqrt((v[627] * v[627]) + (v[628] * v[628]
			) + (v[629] * v[629]))) {
			if (v[4] > 0.1e-5) { v013 = 1e0 / v[4]; v014 = (-(v013 / v[4])); v015 = (2e0*v013) / (v[4] * v[4]); }
			else {
				v013 = (24000000e0 - (-1e0 + 1000000e0*v[4])*(71999994e0 - 0.71999982e14*v[4] + 0.6e19*Power(v[4], 3)
					+ 0.23999982e20*(v[4] * v[4]))) / 24e0;
				v014 = (-500000e0 / 3e0)*(35999994e0 - 0.47999982e14*v[4] + 0.6e19*Power(v[4], 3) + 0.17999982e20*
					(v[4] * v[4]));
				v015 = 0.1e13*(7999997e0 - 0.5999994e13*v[4] - 0.3e13*(v[4] * v[4]));
			};
			v[665] = sqrt((v[630] * v[630]) + (v[631] * v[631]) + (v[632] * v[632]));
			if (v[665] > 0.1e-5) { v016 = 1e0 / v[665]; v017 = (-(v016 / v[665])); v018 = (2e0*v016) / (v[665] * v[665]); }
			else {
				v016 = (24000000e0 - (-1e0 + 1000000e0*v[665])*(71999994e0 - 0.71999982e14*v[665] + 0.6e19*Power(v[665], 3)
					+ 0.23999982e20*(v[665] * v[665]))) / 24e0;
				v017 = (-500000e0 / 3e0)*(35999994e0 - 0.47999982e14*v[665] + 0.6e19*Power(v[665], 3) + 0.17999982e20*
					(v[665] * v[665]));
				v018 = 0.1e13*(7999997e0 - 0.5999994e13*v[665] - 0.3e13*(v[665] * v[665]));
			};
			v[672] = -(v013*v016*v[6] * sqrt(v[2029]));
			v[671] = v[576] + v[630] * v[672];
			v[673] = v[580] + v[631] * v[672];
			v[674] = v[581] + v[632] * v[672];
		}
		else {
			v[671] = 0e0;
			v[673] = 0e0;
			v[674] = 0e0;
		};
	}
	else {
		b675 = sqrt((v[636] * v[636]) + (v[637] * v[637]) + (v[638] * v[638])) <= v[6] * sqrt((v[627] * v[627]) +
			(v[628] * v[628]) + (v[629] * v[629]));
		if (b675) {
			v[643] = v[636];
			v[644] = v[637];
			v[645] = v[638];
			v[646] = 1e0;
		}
		else {
			v[686] = sqrt(v[2029]);
			v[3372] = v[6] * v[686];
			v[677] = sqrt(v[2025]);
			if (v[677] > 0.1e-5) { v019 = 1e0 / v[677]; v020 = (-(v019 / v[677])); v021 = (2e0*v019) / (v[677] * v[677]); }
			else {
				v019 = (24000000e0 - (-1e0 + 1000000e0*v[677])*(71999994e0 - 0.71999982e14*v[677] + 0.6e19*Power(v[677], 3)
					+ 0.23999982e20*(v[677] * v[677]))) / 24e0;
				v020 = (-500000e0 / 3e0)*(35999994e0 - 0.47999982e14*v[677] + 0.6e19*Power(v[677], 3) + 0.17999982e20*
					(v[677] * v[677]));
				v021 = 0.1e13*(7999997e0 - 0.5999994e13*v[677] - 0.3e13*(v[677] * v[677]));
			};
			v[681] = v020;
			v[682] = v019;
			v[683] = v[636] * v[682];
			v[684] = v[637] * v[682];
			v[685] = v[638] * v[682];
			v[643] = v[3372] * v[683];
			v[644] = v[3372] * v[684];
			v[645] = v[3372] * v[685];
			v[646] = 0e0;
		};
		if (sqrt((v[630] * v[630]) + (v[631] * v[631]) + (v[632] * v[632])) > v[6] * sqrt((v[627] * v[627]) + (v[628] * v[628]
			) + (v[629] * v[629]))) {
			if (v[4] > 0.1e-5) { v022 = 1e0 / v[4]; v023 = (-(v022 / v[4])); v024 = (2e0*v022) / (v[4] * v[4]); }
			else {
				v022 = (24000000e0 - (-1e0 + 1000000e0*v[4])*(71999994e0 - 0.71999982e14*v[4] + 0.6e19*Power(v[4], 3)
					+ 0.23999982e20*(v[4] * v[4]))) / 24e0;
				v023 = (-500000e0 / 3e0)*(35999994e0 - 0.47999982e14*v[4] + 0.6e19*Power(v[4], 3) + 0.17999982e20*
					(v[4] * v[4]));
				v024 = 0.1e13*(7999997e0 - 0.5999994e13*v[4] - 0.3e13*(v[4] * v[4]));
			};
			v[695] = sqrt((v[630] * v[630]) + (v[631] * v[631]) + (v[632] * v[632]));
			if (v[695] > 0.1e-5) { v025 = 1e0 / v[695]; v026 = (-(v025 / v[695])); v027 = (2e0*v025) / (v[695] * v[695]); }
			else {
				v025 = (24000000e0 - (-1e0 + 1000000e0*v[695])*(71999994e0 - 0.71999982e14*v[695] + 0.6e19*Power(v[695], 3)
					+ 0.23999982e20*(v[695] * v[695]))) / 24e0;
				v026 = (-500000e0 / 3e0)*(35999994e0 - 0.47999982e14*v[695] + 0.6e19*Power(v[695], 3) + 0.17999982e20*
					(v[695] * v[695]));
				v027 = 0.1e13*(7999997e0 - 0.5999994e13*v[695] - 0.3e13*(v[695] * v[695]));
			};
			v[701] = -(v022*v025*v[6] * sqrt(v[2029]));
			v[671] = v[576] + v[630] * v[701];
			v[673] = v[580] + v[631] * v[701];
			v[674] = v[581] + v[632] * v[701];
		}
		else {
			v[671] = 0e0;
			v[673] = 0e0;
			v[674] = 0e0;
		};
	};
}
else {
	v[643] = 0e0;
	v[644] = 0e0;
	v[645] = 0e0;
};
fn[0] = v[627];
fn[1] = v[628];
fn[2] = v[629];
ft[0] = v[643];
ft[1] = v[644];
ft[2] = v[645];
(*stickupdated) = v[646];
gtpupdated[0] = v[576] - v[671];
gtpupdated[1] = v[580] - v[673];
gtpupdated[2] = v[581] - v[674];
b715 = b603;
if (b715) {
	b716 = b605;
}
else {
};
b720 = b277;
if (b720) {
	v[721] = 0e0;
	v[722] = 0e0;
	v[723] = 0e0;
}
else {
	v[723] = 0e0;
	v[722] = 0e0;
	v[721] = 0e0;
};
v[728] = v[252] * v[721] + v[251] * v[722] + v[250] * v[723];
v[730] = v[263] * v[728];
v[3373] = v[730] / v[602];
v[731] = v[252] * v[3373] + v[610] + v[264] * v[721];
v[3391] = v[731] / 2e0;
v[733] = v[251] * v[3373] + v[609] + v[264] * v[722];
v[3390] = v[733] / 2e0;
v[734] = v[250] * v[3373] + v[607] + v[264] * v[723];
v[5253] = 0e0;
v[5254] = 0e0;
v[5255] = 0e0;
v[5256] = v[734];
v[5257] = v[733];
v[5258] = v[731];
v[5259] = 0e0;
v[5260] = 0e0;
v[5261] = 0e0;
v[5262] = 0e0;
v[5263] = 0e0;
v[5264] = 0e0;
v[5265] = 0e0;
v[5266] = 0e0;
v[5267] = 0e0;
v[5238] = 0e0;
v[5239] = 0e0;
v[5240] = 0e0;
v[5241] = 0e0;
v[5242] = 0e0;
v[5243] = 0e0;
v[5244] = v[734];
v[5245] = v[733];
v[5246] = v[731];
v[5247] = 0e0;
v[5248] = 0e0;
v[5249] = 0e0;
v[5250] = 0e0;
v[5251] = 0e0;
v[5252] = 0e0;
v[5268] = v[734];
v[5269] = v[733];
v[5270] = v[731];
v[5271] = 0e0;
v[5272] = 0e0;
v[5273] = 0e0;
v[5274] = 0e0;
v[5275] = 0e0;
v[5276] = 0e0;
v[5277] = 0e0;
v[5278] = 0e0;
v[5279] = 0e0;
v[5280] = 0e0;
v[5281] = 0e0;
v[5282] = 0e0;
v[3389] = v[734] / 2e0;
v[735] = v[731] * v[738];
v[3388] = v[735] / 2e0;
v[736] = v[731] * v[740];
v[760] = v[172] * v[736];
v[737] = v[731] * v[742];
v[766] = v[172] * v[737];
v[739] = v[733] * v[738];
v[761] = v[172] * v[739];
v[741] = v[733] * v[740];
v[3387] = v[741] / 2e0;
v[763] = v[3332] * v[741];
v[743] = v[733] * v[742];
v[769] = v[172] * v[743];
v[744] = v[734] * v[738];
v[767] = v[172] * v[744];
v[745] = v[734] * v[740];
v[770] = v[172] * v[745];
v[746] = v[734] * v[742];
v[3386] = v[746] / 2e0;
v[768] = v[3332] * v[746];
v[747] = -(v[199] * v[731]) - v[196] * v[733] - v[193] * v[734];
v[748] = -(v[198] * v[731]) - v[195] * v[733] - v[192] * v[734];
v[749] = -(v[197] * v[731]) - v[194] * v[733] - v[191] * v[734];
v[750] = v[152] * v[731] + v[148] * v[733] + v[144] * v[734];
v[751] = v[760] + v[761];
v[752] = v[766] + v[767];
v[753] = v[769] + v[770];
v[754] = v[319] * v[3386] + v[338] * v[3387] + v[3388] * v[356] + v[351] * v[736] + v[347] * v[737] + v[343] * v[739]
+ v[334] * v[743] + v[330] * v[744] + v[324] * v[745];
v[903] = -(v[3330] * v[754]) + v[763] + v[768];
v[900] = v[3332] * v[735] - v[763] + v[903];
v[897] = v[763] - v[768] + v[900];
v[4904] = v[128] * v[734];
v[4905] = v[128] * v[733];
v[4906] = v[128] * v[731];
v[4907] = v[130] * v[734];
v[4908] = v[130] * v[733];
v[4909] = v[130] * v[731];
v[4910] = v[132] * v[734];
v[4911] = v[132] * v[733];
v[4912] = v[132] * v[731];
v[4913] = -v[734];
v[4914] = -v[733];
v[4915] = -v[731];
v[4916] = v[314] * v[752] + v[311] * v[753] + v[760] - v[761] + v[310] * v[897];
v[4917] = v[314] * v[751] + v[312] * v[753] - v[766] + v[767] + v[313] * v[900];
v[4918] = v[311] * v[751] + v[312] * v[752] + v[769] - v[770] + v[315] * v[903];
v[755] = v[749] * v[76] + v[748] * v[77] + v[747] * v[78];
v[756] = (-v[755] + v[749] * v[82] + v[748] * v[83] + v[747] * v[84]) / 2e0;
v[757] = (-v[755] + v[749] * v[79] + v[748] * v[80] + v[747] * v[81]) / 2e0;
v[758] = (v[154] * v[731] + v[150] * v[733] + v[146] * v[734] - v[750]) / 2e0;
v[759] = (v[153] * v[731] + v[149] * v[733] + v[145] * v[734] - v[750]) / 2e0;
for (i710 = 1; i710 <= 15; i710++) {
	v[799] = (i710 == 13 ? 1 : 0);
	v[796] = (i710 == 14 ? 1 : 0);
	v[793] = (i710 == 15 ? 1 : 0);
	v[784] = v[4933 + i710];
	v[783] = v[4918 + i710];
	v[777] = v[4963 + i710];
	v[776] = v[4948 + i710];
	v[3374] = -v[776] - v[777];
	v[779] = v[4997 + i710];
	v[805] = -(v[3330] * v[779]);
	v[780] = v[5012 + i710];
	v[781] = v[5027 + i710];
	v[782] = v[5042 + i710];
	v[785] = (-v[783] - v[784]) / 2e0;
	v[786] = v[5057 + i710];
	v[787] = v[5117 + i710];
	v[788] = v[5177 + i710];
	v[789] = (v[3374] * v[76] + v[776] * v[79] + v[777] * v[82]) / 2e0;
	v[790] = (v[3374] * v[77] + v[776] * v[80] + v[777] * v[83]) / 2e0;
	v[791] = (v[3374] * v[78] + v[776] * v[81] + v[777] * v[84]) / 2e0;
	v[792] = v[780] - v[793];
	v[794] = v[780] + v[793];
	v[795] = v[781] + v[796];
	v[797] = v[781] - v[796];
	v[798] = v[782] - v[799];
	v[800] = v[782] + v[799];
	v[802] = -(v[319] * v[3715] * v[779]) + v[3332] * v[786];
	v[804] = v[172] * v[792] + v[324] * v[805];
	v[806] = v[172] * v[795] + v[330] * v[805];
	v[807] = v[5192 + i710] + v[3375] * v[783] + v[3376] * v[784] + v[144] * v[785] - v[191] * v[789] - v[192] * v[790]
		- v[193] * v[791] + v[742] * v[802] + v[740] * v[804] + v[738] * v[806];
	v[808] = v[172] * v[794] + v[334] * v[805];
	v[809] = (-(v[172] * v[787]) + v[338] * v[805]) / 2e0;
	v[810] = v[172] * v[798] + v[343] * v[805];
	v[811] = v[5207 + i710] + v[3377] * v[783] + v[3378] * v[784] + v[148] * v[785] - v[194] * v[789] - v[195] * v[790]
		- v[196] * v[791] + v[742] * v[808] + v[740] * v[809] + v[738] * v[810];
	v[812] = v[172] * v[797] + v[347] * v[805];
	v[813] = v[734] * v[802] + v[733] * v[808] + v[731] * v[812];
	v[814] = v[172] * v[800] + v[351] * v[805];
	v[815] = v[734] * v[804] + v[733] * v[809] + v[731] * v[814];
	v[816] = (-(v[172] * v[788]) + v[356] * v[805]) / 2e0;
	v[817] = v[5222 + i710] + v[3379] * v[783] + v[3380] * v[784] + v[152] * v[785] - v[197] * v[789] - v[198] * v[790]
		- v[199] * v[791] + v[742] * v[812] + v[740] * v[814] + v[738] * v[816];
	v[3381] = v[250] * v[807] + v[251] * v[811] + v[252] * v[817];
	v[818] = v[734] * v[806] + v[733] * v[810] + v[731] * v[816];
	v[819] = v[3381] / v[602];
	v[821] = -(v[3381] * v[730] * v[820]);
	v[840] = v[821];
	v[822] = v[263] * v[819];
	v[823] = v[807];
	v[838] = v[823];
	v[824] = v[811];
	v[836] = v[824];
	v[825] = v[817];
	v[834] = v[825];
	v[826] = 0e0;
	v[827] = 0e0;
	v[828] = 0e0;
	b829 = b603;
	if (b829) {
		v[3385] = v[279] * v[823];
		v[3384] = v[280] * v[824];
		v[3383] = v[281] * v[825];
		b830 = b605;
		if (b830) {
			v[828] = v[3368] * v[825];
			v[825] = 0e0;
			v[827] = v[3368] * v[824];
			v[824] = 0e0;
			v[826] = v[3368] * v[823];
			v[823] = 0e0;
			v[821] = v[821] - v[31] * (-v[3383] - v[3384] - v[3385])*v[34] * Power(v[614], v[615]);
		}
		else {
			v[828] = v[611] * v[834];
			v[825] = 0e0;
			v[827] = v[611] * v[836];
			v[824] = 0e0;
			v[826] = v[611] * v[838];
			v[823] = 0e0;
			v[821] = v[840] + v[3382] * (v[3383] + v[3384] + v[3385])*Power(v[602], v[620]);
		};
	}
	else {
	};
	b841 = b277;
	if (b841) {
		v[842] = -v[828];
		v[843] = -v[827];
		v[844] = -v[826];
	}
	else {
		v[842] = v[828];
		v[843] = v[827];
		v[844] = v[826];
	};
	v[821] = v[262] * v[728] * v[819] + v[821] + v[263] * (v[723] * v[807] + v[722] * v[811] + v[721] * v[817] + v[252] * v[842]
		+ v[251] * v[843] + v[250] * v[844]);
	v[852] = (v[730] * v[817] + v[252] * v[821]) / v[602] + v[721] * v[822] + v[264] * v[842];
	v[854] = (v[730] * v[811] + v[251] * v[821]) / v[602] + v[722] * v[822] + v[264] * v[843];
	v[856] = (v[730] * v[807] + v[250] * v[821]) / v[602] + v[723] * v[822] + v[264] * v[844];
	v[857] = -(v[731] * v[791]) - v[168] * v[852];
	v[858] = -(v[731] * v[790]) - v[167] * v[852];
	v[859] = -(v[731] * v[789]) - v[166] * v[852];
	v[860] = -(v[733] * v[791]) - v[168] * v[854];
	v[861] = -(v[733] * v[790]) - v[167] * v[854];
	v[862] = -(v[733] * v[789]) - v[166] * v[854];
	v[863] = -(v[734] * v[791]) - v[168] * v[856];
	v[864] = -(v[734] * v[790]) - v[167] * v[856];
	v[865] = -(v[734] * v[789]) - v[166] * v[856];
	v[866] = v[5267 + i710] + v[152] * v[852] + v[148] * v[854] + v[144] * v[856];
	v[867] = v[859] * v[94] + v[858] * v[95] + v[857] * v[96];
	v[868] = v[859] * v[91] + v[858] * v[92] + v[857] * v[93];
	v[869] = v[859] * v[88] + v[858] * v[89] + v[857] * v[90];
	v[870] = v[862] * v[94] + v[861] * v[95] + v[860] * v[96];
	v[871] = v[862] * v[91] + v[861] * v[92] + v[860] * v[93];
	v[872] = v[862] * v[88] + v[861] * v[89] + v[860] * v[90];
	v[873] = v[865] * v[94] + v[864] * v[95] + v[863] * v[96];
	v[874] = v[865] * v[91] + v[864] * v[92] + v[863] * v[93];
	v[875] = v[865] * v[88] + v[864] * v[89] + v[863] * v[90];
	v[876] = (v[735] * v[805] + v[172] * v[867]) / 2e0;
	v[877] = v[736] * v[805] + v[172] * v[868];
	v[878] = v[737] * v[805] + v[172] * v[869];
	v[879] = v[739] * v[805] + v[172] * v[870];
	v[881] = v[743] * v[805] + v[172] * v[872];
	v[882] = v[744] * v[805] + v[172] * v[873];
	v[883] = v[745] * v[805] + v[172] * v[874];
	v[3392] = v[3329] * v[754] * v[779] - v[3330] * (-(v[3386] * v[786]) - v[3387] * v[787] - v[3388] * v[788]
		+ v[745] * v[792] + v[743] * v[794] + v[744] * v[795] + v[737] * v[797] + v[739] * v[798] + v[736] * v[800]
		+ v[3426] * v[867] + v[351] * v[868] + v[347] * v[869] + v[343] * v[870] + v[3427] * v[871] + v[334] * v[872]
		+ v[330] * v[873] + v[324] * v[874] + v[3428] * v[875]);
	v[899] = v[3392] - v[3387] * v[805] + v[3332] * v[871];
	v[885] = (v[746] * v[805] + v[172] * v[875]) / 2e0;
	v[886] = -(v[199] * v[852]) - v[196] * v[854] - v[193] * v[856] - v[813] * v[90] - v[815] * v[93] - v[818] * v[96];
	v[887] = -(v[198] * v[852]) - v[195] * v[854] - v[192] * v[856] - v[813] * v[89] - v[815] * v[92] - v[818] * v[95];
	v[888] = -(v[197] * v[852]) - v[194] * v[854] - v[191] * v[856] - v[813] * v[88] - v[815] * v[91] - v[818] * v[94];
	v[889] = v[78] * v[886] + v[77] * v[887] + v[76] * v[888];
	v[890] = (v[84] * v[886] + v[83] * v[887] + v[82] * v[888] - v[889]) / 2e0;
	v[891] = (v[81] * v[886] + v[80] * v[887] + v[79] * v[888] - v[889]) / 2e0;
	v[892] = v[878] + v[882];
	v[893] = v[877] + v[879];
	v[894] = v[881] + v[883];
	v[5313] = 0e0;
	v[5314] = 0e0;
	v[5315] = 0e0;
	v[5316] = 0e0;
	v[5317] = 0e0;
	v[5318] = 0e0;
	v[5319] = 0e0;
	v[5320] = 0e0;
	v[5321] = 0e0;
	v[5322] = 0e0;
	v[5323] = 0e0;
	v[5324] = 0e0;
	v[5325] = 0e0;
	v[5326] = v[753];
	v[5327] = v[752];
	v[5298] = 0e0;
	v[5299] = 0e0;
	v[5300] = 0e0;
	v[5301] = 0e0;
	v[5302] = 0e0;
	v[5303] = 0e0;
	v[5304] = 0e0;
	v[5305] = 0e0;
	v[5306] = 0e0;
	v[5307] = 0e0;
	v[5308] = 0e0;
	v[5309] = 0e0;
	v[5310] = v[753];
	v[5311] = 0e0;
	v[5312] = v[751];
	v[5283] = 0e0;
	v[5284] = 0e0;
	v[5285] = 0e0;
	v[5286] = 0e0;
	v[5287] = 0e0;
	v[5288] = 0e0;
	v[5289] = 0e0;
	v[5290] = 0e0;
	v[5291] = 0e0;
	v[5292] = 0e0;
	v[5293] = 0e0;
	v[5294] = 0e0;
	v[5295] = v[752];
	v[5296] = v[751];
	v[5297] = 0e0;
	v[5373] = v[734] * v[785] + v[128] * v[856];
	v[5374] = v[733] * v[785] + v[128] * v[854];
	v[5375] = v[731] * v[785] + v[128] * v[852];
	v[5376] = v[3389] * v[783] + v[130] * v[856];
	v[5377] = v[3390] * v[783] + v[130] * v[854];
	v[5378] = v[3391] * v[783] + v[130] * v[852];
	v[5379] = v[3389] * v[784] + v[132] * v[856];
	v[5380] = v[3390] * v[784] + v[132] * v[854];
	v[5381] = v[3391] * v[784] + v[132] * v[852];
	v[5382] = -v[856];
	v[5383] = -v[854];
	v[5384] = -v[852];
	v[5385] = v[5312 + i710] / 2e0 + v[877] - v[879] + v[314] * v[892] + v[311] * v[894] + 2e0*v[799] * v[897] + v[310] * (
		-v[876] + v[899]);
	v[5386] = v[5297 + i710] / 2e0 - v[878] + v[882] + v[313] * (v[3392] - v[876] - v[885]) + v[314] * v[893] + v[312] * v[894]
		+ 2e0*v[796] * v[900];
	v[5387] = v[5282 + i710] / 2e0 + v[881] - v[883] + v[312] * v[892] + v[311] * v[893] + v[315] * (-v[885] + v[899])
		+ 2e0*v[793] * v[903];
	v[895] = (v[5237 + i710] + v[154] * v[852] + v[150] * v[854] + v[146] * v[856] - v[866]) / 2e0;
	v[896] = (v[5252 + i710] + v[153] * v[852] + v[149] * v[854] + v[145] * v[856] - v[866]) / 2e0;
	Rc[i710 - 1] += v[4903 + i710] + v[757] * v[776] + v[756] * v[777] + v[759] * v[783] + v[758] * v[784];
	for (i773 = 1; i773 <= 15; i773++) {
		Kc[i710 - 1][i773 - 1] += v[5372 + i773] + v[4963 + i773] * v[890] + v[4948 + i773] * v[891] + v[4933 + i773] * v[895]
			+ v[4918 + i773] * v[896];
	};/* end for */
};/* end for */
v[973] = 0e0;
v[974] = 0e0;
v[975] = 0e0;
b976 = b603;
if (b976) {
	b977 = b29;
	if (b977) {
		b978 = b641;
		if (b978) {
			v[975] = 0e0;
			v[974] = 0e0;
			v[973] = 0e0;
		}
		else {
		};
	}
	else {
		b979 = b675;
		if (b979) {
			v[975] = 0e0;
			v[974] = 0e0;
			v[973] = 0e0;
		}
		else {
		};
	};
}
else {
};
v[3395] = v[9] * v[973];
v[3394] = v[9] * v[974];
v[3393] = v[9] * v[975];
v[999] = v[1019] * v[3393];
v[1000] = v[1021] * v[3393];
v[3435] = (v[1023] * v[3393]) / 2e0;
v[3436] = (v[1025] * v[3393]) / 2e0;
v[1020] = v[1019] * v[3394];
v[1022] = v[1021] * v[3394];
v[3437] = (v[1023] * v[3394]) / 2e0;
v[3438] = (v[1025] * v[3394]) / 2e0;
v[1043] = v[1019] * v[3395];
v[1044] = v[1021] * v[3395];
v[3439] = (v[1023] * v[3395]) / 2e0;
v[3440] = (v[1025] * v[3395]) / 2e0;
v[1047] = v[4] * v[975];
v[1780] = -(v[1047] * v[281]);
v[1048] = v[4] * v[974];
v[1781] = -(v[1048] * v[280]);
v[1783] = v[1780] + v[1781];
v[1049] = v[4] * v[973];
v[1776] = -(v[1049] * v[279]);
v[1784] = v[1776] + v[1780];
v[1782] = v[1776] + v[1781];
v[1050] = 0e0;
v[1051] = 0e0;
v[1052] = 0e0;
v[1053] = 0e0;
v[1054] = 0e0;
b1055 = b603;
if (b1055) {
	v[1056] = 0e0;
	v[1057] = 0e0;
	v[1058] = 0e0;
	b1059 = b622;
	if (b1059) {
		v[1058] = 0e0;
		v[1057] = 0e0;
		v[1056] = 0e0;
	}
	else {
	};
	v[1063] = v[1056] * v[583] + v[1057] * v[591] + v[1058] * v[593];
	v[3396] = v[1063] * v[7];
	b1060 = b605;
	if (b1060) {
		v[1052] = v[1058] * v[617];
		v[1051] = v[1057] * v[617];
		v[1050] = v[1056] * v[617];
		v[1054] = (v[1064] * v[3396] * v[615] * Power(v[614], v[1987])) / v[1065];
	}
	else {
		v[1052] = v[1058] * v[621];
		v[1051] = v[1057] * v[621];
		v[1050] = v[1056] * v[621];
		v[1053] = (v[1072] * v[3396] * v[620] * Power(v[602], v[1990])) / v[1073];
	};
}
else {
};
v[3488] = v[1050] * v[279];
v[3442] = v[1050] * v[582];
v[1189] = -(v[1050] * v[3361]);
v[1860] = -(v[1051] * v[587]);
v[3500] = v[1052] * v[281];
v[1098] = v[1052] * v[3397];
v[1096] = v[280] * v[3500];
v[3399] = v[1096] - v[1189] - v[1860];
v[1863] = -v[1096] + v[1860];
v[1074] = 0e0;
v[1075] = 0e0;
v[1076] = 0e0;
b1077 = b603;
if (b1077) {
	b1078 = b605;
	if (b1078) {
		v[1076] = 0e0;
		v[1075] = 0e0;
		v[1074] = 0e0;
		v[1053] = v[1053] - v[1054];
	}
	else {
	};
}
else {
};
v[1074] = v[1074] + v[1052] * v[1963];
v[1075] = v[1075] + v[1052] * v[1952];
v[1085] = v[1052] * v[1938];
v[1076] = v[1076] + v[1052] * v[1909];
v[1074] = v[1074] + v[1051] * v[1900];
v[1121] = v[1051] * v[1886];
v[1075] = v[1075] + v[1051] * v[1871];
v[1126] = v[1051] * v[280];
v[3492] = -(v[1126] * v[281]) - v[1050] * v[3397];
v[3400] = v[3492] - v[1052] * v[592];
v[1861] = -(v[1126] * v[279]);
v[3398] = v[1098] - v[1861] + v[3442];
v[3441] = v[1861] + v[3398];
v[1076] = v[1076] + v[1051] * v[1858];
v[1156] = v[1050] * v[1838];
v[1157] = v[1051] * v[212] + v[1050] * v[213];
v[1158] = v[1051] * v[215] + v[1050] * v[216];
v[1159] = v[1051] * v[218] + v[1050] * v[219];
v[1658] = v[1157] * v[128] + v[1158] * v[130] + v[1159] * v[132];
v[1074] = v[1074] + v[1050] * v[1823];
v[1075] = v[1075] + v[1050] * v[1813];
v[1168] = v[1126] + v[3488];
v[1076] = v[1076] + v[1050] * v[1802];
v[1172] = v[1052] * v[1086] + v[1051] * v[1129] + v[1050] * v[1171] - v[9] * (v[966] * v[973] + v[967] * v[974]
	+ v[968] * v[975]);
v[1361] = v[1172] * v[1440];
v[1335] = v[1172] * v[2512];
v[1174] = v[1052] * v[1088] + v[1051] * v[1131] + v[1050] * v[1173] - v[9] * (v[962] * v[973] + v[963] * v[974]
	+ v[964] * v[975]);
v[1348] = v[1174] * v[2510];
v[3409] = v[1348] + v[1172] * v[1439];
v[3408] = v[1335] + v[1174] * v[1437];
v[1176] = v[1052] * v[1090] + v[1051] * v[1133] + v[1050] * v[1175] - v[9] * (v[958] * v[973] + v[959] * v[974]
	+ v[960] * v[975]);
v[1359] = v[1176] * v[2507];
v[3595] = v[1359] + v[1361];
v[3410] = v[1359] + v[1174] * v[1438];
v[3589] = v[1361] + v[3410];
v[1350] = v[1176] * v[3341];
v[3593] = v[1350] + v[3409];
v[3591] = v[1348] + v[1350];
v[1338] = v[1176] * v[3344];
v[3594] = v[1338] + v[3408];
v[3590] = v[1335] + v[1338];
v[1186] = -(v[239] * v[3398]);
v[1187] = -(v[249] * v[3398]);
v[1188] = -(v[229] * v[3398]);
v[1190] = -(v[229] * v[3399]);
v[1191] = -(v[249] * v[3399]);
v[1193] = v[249] * v[3400];
v[1194] = -(v[239] * v[3399]);
v[1195] = v[229] * v[3400];
v[1196] = v[239] * v[3400];
v[1200] = 0e0;
v[1201] = 0e0;
v[1202] = 0e0;
v[1203] = 0e0;
v[1204] = 0e0;
v[1205] = 0e0;
v[1206] = 0e0;
v[1207] = 0e0;
v[1208] = 0e0;
b1209 = b30;
if (b1209) {
	v[1085] = v[1085] - v[1047] * v[579];
	v[1121] = v[1121] - v[1048] * v[578];
	v[1211] = v[1049] * v[1210] + v[1783] * v[279];
	v[1213] = v[1048] * v[1212] + v[1784] * v[280];
	v[1215] = v[1047] * v[1214] + v[1782] * v[281];
	v[1156] = v[1156] - v[1049] * v[577];
	v[1074] = v[1074] - v[1049] * v[3401] + v[1783] * v[577];
	v[1075] = v[1075] + v[1048] * v[3402] + v[1784] * v[578];
	v[1076] = v[1076] + v[1047] * v[3403] + v[1782] * v[579];
	v[1200] = v[1] * v[1211];
	v[1201] = v[1211] * v[2];
	v[1202] = v[1211] * v[3];
	v[1219] = -v[1211];
	v[1220] = -(v[1211] * v[171]);
	v[1221] = -(v[1211] * v[170]);
	v[1222] = -(v[1211] * v[169]);
	v[1223] = v[1211] * v[137];
	v[1224] = v[1211] * v[135];
	v[1225] = v[1211] * v[133];
	v[1203] = v[1] * v[1213];
	v[1204] = v[1213] * v[2];
	v[1205] = v[1213] * v[3];
	v[1226] = -v[1213];
	v[1227] = -(v[1213] * v[171]);
	v[1228] = -(v[1213] * v[170]);
	v[1229] = -(v[1213] * v[169]);
	v[1230] = v[1213] * v[137];
	v[1231] = v[1213] * v[135];
	v[1232] = v[1213] * v[133];
	v[1206] = v[1] * v[1215];
	v[1207] = v[1215] * v[2];
	v[1208] = v[1215] * v[3];
	v[1233] = -v[1215];
	v[1234] = -(v[1215] * v[171]);
	v[1235] = -(v[1215] * v[170]);
	v[1236] = -(v[1215] * v[169]);
	v[1237] = v[1215] * v[137];
	v[1238] = v[1215] * v[135];
	v[1239] = v[1215] * v[133];
}
else {
	v[1225] = 0e0;
	v[1224] = 0e0;
	v[1223] = 0e0;
	v[1232] = 0e0;
	v[1231] = 0e0;
	v[1230] = 0e0;
	v[1239] = 0e0;
	v[1238] = 0e0;
	v[1237] = 0e0;
	v[1222] = 0e0;
	v[1221] = 0e0;
	v[1220] = 0e0;
	v[1229] = 0e0;
	v[1228] = 0e0;
	v[1227] = 0e0;
	v[1236] = 0e0;
	v[1235] = 0e0;
	v[1234] = 0e0;
	v[1219] = 0e0;
	v[1226] = 0e0;
	v[1233] = 0e0;
};
v[3479] = v[1200] / 2e0;
v[3480] = v[1204] / 2e0;
v[3481] = v[1208] / 2e0;
b1240 = b542;
if (b1240) {
	v[1274] = -(v[1201] * v[559]);
	v[1269] = v[1202] * v[559];
	v[1256] = v[1205] * v[559];
	v[1243] = -(v[1208] * v[3364]);
	v[1247] = v[1207] * v[559];
	v[1251] = v[1206] * v[559];
	v[1255] = v[1247] + v[1256];
	v[1260] = -(v[1204] * v[3364]);
	v[1264] = v[1203] * v[559];
	v[1268] = v[1251] + v[1269];
	v[1275] = v[1264] - v[1274];
	v[1277] = v[1207] * v[1245] + v[1206] * v[1249] + v[1205] * v[1253] + v[1203] * v[1262] + v[1202] * v[1266]
		+ v[1201] * v[1271] + v[1276] * v[3479] + v[1258] * v[3480] + v[1241] * v[3481];
	v[1727] = v[1243] + v[1260] - (4e0*v[1277]) / (v[1281] * v[1281]);
	v[3478] = 4e0*v[1727];
	v[1725] = -v[1243] + v[1727] - v[1200] * v[3364];
	v[3477] = 4e0*(v[1243] - v[1260] + v[1725]);
	v[1282] = v[1264] + v[1274] + v[1255] * v[3363] + v[1268] * v[3404] + 2e0*v[1725] * v[558];
	v[1284] = (-2e0*v[1251] + 2e0*v[1269] + v[1275] * v[556] + v[3477] * v[557] + v[1255] * v[558]) / 2e0;
	v[1285] = (2e0*v[1247] - 2e0*v[1256] + v[3478] * v[556] + v[1275] * v[557] + v[1268] * v[558]) / 2e0;
	v[3405] = v[1285] * v[544] + v[1284] * v[545] + v[1282] * v[546];
	v[1713] = v[3405] * v[555];
	v[1710] = v[3405] * v[549];
	v[1288] = v[1710] * v[554] + v[1713] / (Power(cos(v[1286]), 2)*sqrt(v[1714]));
	v[3406] = v[1288] / v[547];
	v[1289] = v[1282] * v[3362] + v[3406] * v[546];
	v[1291] = v[1284] * v[3362] + v[3406] * v[545];
	v[1292] = v[1285] * v[3362] + v[3406] * v[544];
	v[1074] = v[1074] - v[1289] * v[275] + v[1291] * v[276];
	v[1075] = v[1075] + v[1289] * v[274] - v[1292] * v[276];
	v[1076] = v[1076] - v[1291] * v[274] + v[1292] * v[275];
}
else {
};
v[1074] = v[1074] + v[1156] * v[3472];
v[1074] = v[1074] + v[1658] * v[280];
v[1075] = v[1075] + v[1658] * v[279] + v[1121] * v[3469];
v[1076] = v[1076] + v[1168] * v[1637] + v[1085] * v[3466];
b1293 = b277;
if (b1293) {
	v[1294] = -v[1076];
	v[1295] = -v[1075];
	v[1296] = -v[1074];
}
else {
	v[1294] = v[1076];
	v[1295] = v[1075];
	v[1296] = v[1074];
};
v[1301] = v[1296] * v[250] + v[1295] * v[251] + v[1294] * v[252];
v[1053] = v[1053] + v[1301] * v[263];
v[3407] = v[1053] / v[602];
v[1303] = v[1294] * v[264] + v[252] * v[3407];
v[1304] = v[1295] * v[264] + v[251] * v[3407];
v[1305] = v[1296] * v[264] + v[250] * v[3407];
v[1233] = v[1233] - v[1303];
v[1226] = v[1226] - v[1304];
v[1219] = v[1219] - v[1305];
v[1307] = v[1172] * v[1800];
v[1308] = v[1172] * v[1799];
v[1309] = v[1172] * v[1798];
v[1312] = v[1174] * v[1795];
v[1313] = v[1172] * v[236] + v[1174] * v[238];
v[3414] = v[1313] * v[177];
v[1316] = v[1174] * v[1793];
v[1317] = v[1307] + v[1312];
v[1319] = v[1174] * v[1794] + v[3414] / v[227];
v[1322] = v[1176] * v[1788];
v[1323] = v[1172] * v[230] + v[1176] * v[238];
v[3415] = v[1323] * v[182];
v[1324] = v[1174] * v[230] + v[1176] * v[236];
v[3416] = v[1324] * v[186];
v[3812] = -(v[1172] * v[3308]) - v[1174] * v[3309] - v[1176] * v[3310] + v[3340] * v[3414] + v[226] * v[3415]
+ v[225] * v[3416];
v[1326] = v[1176] * v[1790] + v[3415] / v[227];
v[1327] = v[1319] + v[1326];
v[1329] = v[1176] * v[1789] + v[3416] / v[227];
v[1330] = v[1308] + v[1329];
v[1234] = v[1234] - v[1303] * v[168] + v[1000] * v[305] + v[295] * v[999];
v[1235] = v[1235] - v[1303] * v[167] + v[1000] * v[303] + v[294] * v[999];
v[1236] = v[1236] - v[1303] * v[166] + v[1000] * v[301] + v[293] * v[999];
v[1343] = v[238] * v[3594] + v[1236] * v[94] + v[1235] * v[95] + v[1234] * v[96];
v[1344] = v[1324] * v[3344] + v[236] * v[3408] + v[1236] * v[91] + v[1235] * v[92] + v[1234] * v[93];
v[1345] = v[230] * v[3590] + v[1236] * v[88] + v[1235] * v[89] + v[1234] * v[90];
v[1227] = v[1227] - v[1304] * v[168] + v[1020] * v[295] + v[1022] * v[305];
v[1228] = v[1228] - v[1304] * v[167] + v[1020] * v[294] + v[1022] * v[303];
v[1229] = v[1229] - v[1304] * v[166] + v[1020] * v[293] + v[1022] * v[301];
v[1355] = v[1323] * v[3341] + v[238] * v[3409] + v[1229] * v[94] + v[1228] * v[95] + v[1227] * v[96];
v[1356] = v[236] * v[3593] + v[1229] * v[91] + v[1228] * v[92] + v[1227] * v[93];
v[1357] = v[230] * v[3591] + v[1229] * v[88] + v[1228] * v[89] + v[1227] * v[90];
v[1220] = v[1220] - v[1305] * v[168] + v[1043] * v[295] + v[1044] * v[305];
v[1221] = v[1221] - v[1305] * v[167] + v[1043] * v[294] + v[1044] * v[303];
v[1222] = v[1222] - v[1305] * v[166] + v[1043] * v[293] + v[1044] * v[301];
v[1367] = v[1313] * v[1438] + v[238] * v[3595] + v[1222] * v[94] + v[1221] * v[95] + v[1220] * v[96];
v[1368] = v[236] * v[3410] + v[1222] * v[91] + v[1221] * v[92] + v[1220] * v[93];
v[1369] = v[230] * v[3589] + v[1222] * v[88] + v[1221] * v[89] + v[1220] * v[90];
v[1370] = -(v[1186] * v[742]);
v[1371] = -(v[1186] * v[740]);
v[1372] = -(v[1186] * v[738]);
v[1373] = -(v[1187] * v[742]);
v[1374] = -(v[1187] * v[738]);
v[1375] = -(v[1187] * v[740]);
v[1376] = -(v[1188] * v[740]);
v[1377] = -(v[1188] * v[738]);
v[1378] = -(v[1188] * v[742]);
v[1379] = -(v[1190] * v[742]);
v[1380] = -(v[1190] * v[740]);
v[1381] = -(v[1190] * v[738]);
v[1382] = -(v[1191] * v[738]);
v[1383] = -(v[1191] * v[742]);
v[1384] = -(v[1191] * v[740]);
v[1385] = -(v[1193] * v[740]);
v[1386] = -(v[1193] * v[742]);
v[1387] = -(v[1193] * v[738]);
v[1388] = v[1376] + v[1379] + v[1382] + v[1385] - v[1327] * v[326] - v[1317] * v[329] + v[1316] * v[3454];
v[3586] = v[1388] / 2e0;
v[1389] = -(v[1194] * v[742]);
v[1390] = -(v[1194] * v[738]);
v[1391] = -(v[1194] * v[740]);
v[1392] = -v[1371] - v[1374] - v[1386] - v[1389] - v[1327] * v[323] + v[1330] * v[329] + v[1322] * v[3453];
v[3587] = v[1392] / 2e0;
v[1393] = v[1368] * v[172] - v[1376] * v[316] + v[1371] * v[317] - v[1375] * v[318];
v[1394] = -(v[1195] * v[742]);
v[1395] = -(v[1195] * v[740]);
v[1396] = -(v[1195] * v[738]);
v[1397] = -(v[1196] * v[740]);
v[1398] = -(v[1196] * v[742]);
v[1399] = -(v[1196] * v[738]);
v[1403] = -v[1377] - v[1390] - v[1394] - v[1397] - v[1317] * v[323] + v[1330] * v[326] + v[1309] * v[3452];
v[3585] = v[1403] / 2e0;
v[1404] = v[1367] * v[172] - v[1377] * v[316] + v[1372] * v[317] - v[1374] * v[318];
v[1405] = v[1357] * v[172] - v[1379] * v[316] + v[1389] * v[317] - v[1383] * v[318];
v[1406] = v[1373] + v[1384];
v[1407] = v[1355] * v[172] - v[1381] * v[316] + v[1390] * v[317] - v[1382] * v[318];
v[1408] = v[1345] * v[172] - v[1394] * v[316] + v[1398] * v[317] - v[1386] * v[318];
v[1409] = v[1344] * v[172] - v[1395] * v[316] + v[1397] * v[317] - v[1385] * v[318];
v[1410] = v[1380] + v[1396];
v[1411] = v[1370] + v[1399];
v[5951] = 0e0;
v[5952] = 0e0;
v[5953] = 0e0;
v[5954] = 0e0;
v[5955] = 0e0;
v[5956] = 0e0;
v[5957] = 0e0;
v[5958] = 0e0;
v[5959] = 0e0;
v[5960] = 0e0;
v[5961] = 0e0;
v[5962] = 0e0;
v[5963] = -v[1410] - v[3587];
v[5964] = -v[1411] + v[3586];
v[5965] = -v[1406] - v[3585];
v[1412] = 1e0 / (v[227] * v[227]);
v[3599] = -(v[1412] * v[237]);
v[3425] = -(v[1412] * v[236]);
v[3424] = -(v[1412] * v[238]);
v[3421] = -(v[1412] * v[230]);
v[3420] = -(v[1412] * v[225]);
v[3419] = -(v[1172] * v[1412]);
v[3418] = -(v[1412] * v[226]);
v[3417] = -(v[1412] * v[3340]);
v[2491] = -(v[1412] * (v[177] * v[228] + v[3338]));
v[2490] = -(v[1412] * (v[176] * v[228] + v[3339]));
v[2489] = -(v[1412] * v[3411]);
v[2488] = -(v[1412] * v[3342]);
v[2487] = -(v[1412] * v[3412]);
v[2486] = -(v[1412] * v[3413]);
v[2485] = -(v[1412] * v[228]);
v[3567] = v[230] * v[2485];
v[2484] = -(v[1412] * v[231]);
v[3566] = v[236] * v[2484];
v[3565] = -(v[1412] * v[238] * v[240]);
v[2482] = -(v[1412] * v[3414]);
v[2481] = -(v[1412] * v[3415]);
v[2480] = -(v[1412] * v[3416]);
v[2479] = v[1313] * v[3417];
v[2478] = v[1323] * v[3418];
v[2477] = v[1324] * v[3420];
v[2407] = v[1174] * v[3417];
v[2406] = v[3343] * v[3419];
v[2404] = v[1176] * v[2485];
v[3555] = v[2404] + v[2407];
v[3808] = v[2406] + v[3555];
v[2394] = v[1176] * v[3418];
v[2391] = v[1174] * v[2484];
v[3554] = v[2391] + v[2394];
v[2390] = v[3345] * v[3419];
v[3807] = v[2390] + v[3554];
v[2380] = v[1176] * v[3420];
v[2378] = v[1174] * v[3599];
v[2377] = v[240] * v[3419];
v[3552] = v[2377] + v[2380];
v[3806] = v[2378] + v[3552];
v[2370] = v[184] * v[3421];
v[2367] = v[179] * v[3421];
v[2364] = -(v[1412] * v[3422]);
v[2363] = -(v[1412] * v[3423]);
v[2358] = v[182] * v[3424];
v[2357] = v[181] * v[3425];
v[3199] = v[2357] + v[2367];
v[3191] = v[2358] + v[3199];
v[3177] = v[2357] + v[2358];
v[2355] = v[175] * v[3421];
v[3196] = v[2355] + v[2363] + v[2364];
v[3189] = -v[2364] + v[3196];
v[3179] = -v[2363] + v[3196];
v[2351] = v[187] * v[3424];
v[3202] = v[2351] + v[2370];
v[2350] = v[186] * v[3425];
v[3185] = v[2350] + v[2351];
v[3181] = v[2350] + v[3202];
v[1797] = v[236] * v[2486] + v[230] * v[2487] + v[3565];
v[1792] = v[238] * v[2488] + v[230] * v[2489] + v[3566];
v[1787] = v[236] * v[2490] + v[238] * v[2491] + v[3567];
v[1599] = v[177] * v[3417];
v[1597] = v[182] * v[3418];
v[1593] = v[186] * v[3420];
v[2476] = v[1309] + v[1316] + v[1322] + v[1324] * v[1593] + v[1323] * v[1597] + v[1313] * v[1599] + v[1176] * v[1787]
+ v[1174] * v[1792] + v[1172] * v[1797];
v[1413] = v[1372] - v[1375] - v[1381] + v[1383] + v[1395] - v[1398] + v[1388] * v[311] - v[1392] * v[312]
- v[1403] * v[314] + v[1368] * v[324] + v[1367] * v[330] + v[1357] * v[334] + v[1343] * v[3426] + v[1356] * v[3427]
+ v[1369] * v[3428] + v[1355] * v[343] + v[2476] * v[3443] + v[1345] * v[347] + v[1344] * v[351] + v[1406] * v[3562]
+ v[1410] * v[3563] + v[1411] * v[3564];
v[1415] = (-2e0*v[1307] + 2e0*v[1312] - v[1378] * v[319] - v[1380] * v[338] + v[1377] * v[3429] + v[1376] * v[3430]
	+ v[1381] * v[3431] + v[1379] * v[3432] + v[1395] * v[3433] + v[1394] * v[3434] - v[1396] * v[356]) / 2e0;
v[3816] = v[1415] * v[3584] + v[1410] * v[803];
v[1417] = -v[1308] + v[1329] + v[1371] * v[324] + v[1372] * v[330] + v[1389] * v[334] + v[1399] * v[3426]
+ v[1391] * v[3427] + v[1370] * v[3428] + v[1390] * v[343] + v[1398] * v[347] + v[1397] * v[351];
v[3818] = -(v[1417] * v[3584]) + v[1411] * v[803];
v[1418] = (v[1356] * v[172] - v[1380] * v[316] + v[1391] * v[317] - v[1384] * v[318]) / 2e0;
v[1419] = (-2e0*v[1319] + 2e0*v[1326] - v[1373] * v[319] - v[1384] * v[338] + v[1374] * v[3429] + v[1375] * v[3430]
	+ v[1382] * v[3431] + v[1383] * v[3432] + v[1385] * v[3433] + v[1386] * v[3434] - v[1387] * v[356]) / 2e0;
v[3820] = v[1419] * v[3584] + v[1406] * v[803];
v[3810] = v[1415] * v[310] - v[1417] * v[313] + v[1419] * v[315];
v[5966] = 0e0;
v[5967] = 0e0;
v[5968] = 0e0;
v[5969] = 0e0;
v[5970] = 0e0;
v[5971] = 0e0;
v[5972] = 0e0;
v[5973] = 0e0;
v[5974] = 0e0;
v[5975] = 0e0;
v[5976] = 0e0;
v[5977] = 0e0;
v[5978] = 8e0*v[1415];
v[5979] = -8e0*v[1417];
v[5980] = 8e0*v[1419];
v[1436] = -v[1418] + v[1419] * v[1495] + v[1417] * v[1498] + v[1415] * v[1500] - v[1413] * v[3330];
v[2513] = v[1436] + (-(v[1369] * v[172]) + v[1378] * v[316] - v[1370] * v[317] + v[1373] * v[318]) / 2e0;
v[1420] = (v[1343] * v[172] - v[1396] * v[316] + v[1399] * v[317] - v[1387] * v[318]) / 2e0;
v[2511] = v[1418] - v[1420] + v[2513];
v[2508] = -v[1420] + v[1436];
v[1421] = v[1404] + v[1408];
v[1422] = v[1407] + v[1409];
v[5981] = 0e0;
v[5982] = 0e0;
v[5983] = 0e0;
v[5984] = 0e0;
v[5985] = 0e0;
v[5986] = 0e0;
v[5987] = 0e0;
v[5988] = 0e0;
v[5989] = 0e0;
v[5990] = 0e0;
v[5991] = 0e0;
v[5992] = 0e0;
v[5993] = v[1421];
v[5994] = v[1422];
v[5995] = 0e0;
v[1423] = v[1393] + v[1405];
v[6011] = 0e0;
v[6012] = 0e0;
v[6013] = 0e0;
v[6014] = 0e0;
v[6015] = 0e0;
v[6016] = 0e0;
v[6017] = 0e0;
v[6018] = 0e0;
v[6019] = 0e0;
v[6020] = 0e0;
v[6021] = 0e0;
v[6022] = 0e0;
v[6023] = v[1423];
v[6024] = 0e0;
v[6025] = v[1422];
v[6056] = 0e0;
v[6057] = 0e0;
v[6058] = 0e0;
v[6059] = 0e0;
v[6060] = 0e0;
v[6061] = 0e0;
v[6062] = 0e0;
v[6063] = 0e0;
v[6064] = 0e0;
v[6065] = 0e0;
v[6066] = 0e0;
v[6067] = 0e0;
v[6068] = 0e0;
v[6069] = v[1423];
v[6070] = v[1421];
v[1237] = v[1237] + v[1303] * v[132] + v[3435];
v[1230] = v[1230] + v[1304] * v[132] + v[3437];
v[1223] = v[1223] + v[1305] * v[132] + v[3439];
v[1238] = v[1238] + v[130] * v[1303] + v[3436];
v[1231] = v[1231] + v[130] * v[1304] + v[3438];
v[1224] = v[1224] + v[130] * v[1305] + v[3440];
v[1239] = v[1239] + v[128] * v[1303] - v[3435] - v[3436];
v[1232] = v[1232] + v[128] * v[1304] - v[3437] - v[3438];
v[1225] = v[1225] + v[128] * v[1305] - v[3439] - v[3440];
v[5392] = v[1225] - v[643] * v[910] - v[644] * v[911] - v[645] * v[912] + v[10] * (v[128] * v[3441] + v[1051] * v[584]
	+ v[9] * (-(v[910] * v[973]) - v[911] * v[974] - v[912] * v[975]));
v[5393] = v[1232] - v[643] * v[914] - v[644] * v[915] - v[645] * v[916] + v[10] * (-(v[128] * v[1863]) + v[1050] * v[584]
	+ v[9] * (-(v[914] * v[973]) - v[915] * v[974] - v[916] * v[975]));
v[5394] = v[1239] - v[643] * v[918] - v[644] * v[919] - v[645] * v[920] + v[10] * (-(v[128] * v[3400]) + v[9] * (-
(v[918] * v[973]) - v[919] * v[974] - v[920] * v[975]));
v[5395] = v[1224] - v[643] * v[922] - v[644] * v[923] - v[645] * v[924] + v[10] * (v[130] * v[3441] + v[1051] * v[585]
	+ v[9] * (-(v[922] * v[973]) - v[923] * v[974] - v[924] * v[975]));
v[5396] = v[1231] - v[643] * v[926] - v[644] * v[927] - v[645] * v[928] + v[10] * (-(v[130] * v[1863]) + v[1050] * v[585]
	+ v[9] * (-(v[926] * v[973]) - v[927] * v[974] - v[928] * v[975]));
v[5397] = v[1238] - v[643] * v[930] - v[644] * v[931] - v[645] * v[932] + v[10] * (-(v[130] * v[3400]) + v[9] * (-
(v[930] * v[973]) - v[931] * v[974] - v[932] * v[975]));
v[5398] = v[1223] - v[643] * v[934] - v[644] * v[935] - v[645] * v[936] + v[10] * (v[132] * v[3441] + v[1051] * v[586]
	+ v[9] * (-(v[934] * v[973]) - v[935] * v[974] - v[936] * v[975]));
v[5399] = v[1230] - v[643] * v[938] - v[644] * v[939] - v[645] * v[940] + v[10] * (-(v[132] * v[1863]) + v[1050] * v[586]
	+ v[9] * (-(v[938] * v[973]) - v[939] * v[974] - v[940] * v[975]));
v[5400] = v[1237] - v[643] * v[942] - v[644] * v[943] - v[645] * v[944] + v[10] * (-(v[132] * v[3400]) + v[9] * (-
(v[942] * v[973]) - v[943] * v[974] - v[944] * v[975]));
v[5401] = v[1219] - v[643] * v[946] - v[644] * v[947] - v[645] * v[948] + v[10] * (-v[1098] + v[1861] - v[3442] + v[9] * (-
(v[946] * v[973]) - v[947] * v[974] - v[948] * v[975]));
v[5402] = v[1226] - v[643] * v[950] - v[644] * v[951] - v[645] * v[952] + v[10] * (v[1189] + v[1863] + v[9] * (-
(v[950] * v[973]) - v[951] * v[974] - v[952] * v[975]));
v[5403] = v[1233] - v[643] * v[954] - v[644] * v[955] - v[645] * v[956] + v[10] * (v[3400] + v[9] * (-(v[954] * v[973])
	- v[955] * v[974] - v[956] * v[975]));
v[5404] = -v[1407] + v[1409] + v[2508] * v[310] + v[1423] * v[311] + v[1421] * v[314] + v[1392] * v[3332] + 2e0*
(v[1415] * v[3330] + v[1410] * v[3332]) + v[10] * (v[1359] + v[1174] * v[3341] + v[1172] * v[3344]) - v[643] * v[958]
- v[644] * v[959] - v[645] * v[960];
v[5405] = v[1404] - v[1408] + v[10] * (v[1348] + v[1172] * v[1437] + v[1176] * v[1438]) + v[1423] * v[312]
+ v[2511] * v[313] + v[1422] * v[314] - v[1388] * v[3332] + 2e0*(-(v[1417] * v[3330]) + v[1411] * v[3332])
- v[643] * v[962] - v[644] * v[963] - v[645] * v[964];
v[5406] = -v[1393] + v[1405] + v[10] * (v[1335] + v[1174] * v[1439] + v[1176] * v[1440]) + v[1422] * v[311]
+ v[1421] * v[312] + v[2513] * v[315] + v[1403] * v[3332] + 2e0*(v[1419] * v[3330] + v[1406] * v[3332]) - v[643] * v[966]
- v[644] * v[967] - v[645] * v[968];
for (i971 = 1; i971 <= 15; i971++) {
	i3464 = (i971 == 9 ? 1 : 0);
	i3463 = (i971 == 6 ? 1 : 0);
	i3462 = (i971 == 12 ? 1 : 0);
	i3461 = (i971 == 8 ? 1 : 0);
	i3460 = (i971 == 5 ? 1 : 0);
	i3459 = (i971 == 11 ? 1 : 0);
	i3458 = (i971 == 7 ? 1 : 0);
	i3457 = (i971 == 4 ? 1 : 0);
	i3456 = (i971 == 10 ? 1 : 0);
	i3449 = (i971 == 14 ? 1 : 0);
	v[3546] = i3449 * v[10];
	i3448 = (i971 == 13 ? 1 : 0);
	v[3547] = i3448 * v[10];
	i3447 = (i971 == 15 ? 1 : 0);
	v[3548] = i3447 * v[10];
	i3446 = (i971 == 3 ? 1 : 0);
	i3445 = (i971 == 2 ? 1 : 0);
	i3444 = (i971 == 1 ? 1 : 0);
	v[1998] = -i3444 / 2e0;
	v[1999] = -i3445 / 2e0;
	v[2000] = -i3446 / 2e0;
	v[1456] = v[5012 + i971];
	v[1458] = v[5042 + i971];
	v[1460] = v[5027 + i971];
	v[1461] = v[5545 + i971];
	v[1463] = v[4997 + i971];
	v[3583] = 4e0*v[1463];
	v[1506] = v[1463] * v[3330];
	v[1545] = -(v[1506] * v[3443]);
	v[3553] = v[1545] * v[238];
	v[3551] = v[1545] * v[236];
	v[3455] = v[1545] * v[227];
	v[1464] = v[5605 + i971];
	v[1471] = i3444 * v[10];
	v[1473] = i3445 * v[10];
	v[1475] = i3446 * v[10];
	v[1477] = i3457 * v[10];
	v[1479] = i3460 * v[10];
	v[1481] = i3463 * v[10];
	v[1483] = i3458 * v[10];
	v[1485] = i3461 * v[10];
	v[1487] = i3464 * v[10];
	v[3467] = v[128] * v[1475] + v[130] * v[1481] + v[132] * v[1487];
	v[1489] = i3447 + v[1456];
	v[3556] = 2e0*v[1489];
	v[1490] = -i3447 + v[1456];
	v[3557] = 2e0*v[1490];
	v[1491] = i3448 + v[1458];
	v[3558] = 2e0*v[1491];
	v[1492] = -i3448 + v[1458];
	v[3559] = 2e0*v[1492];
	v[1493] = -i3449 + v[1460];
	v[3560] = 2e0*v[1493];
	v[1494] = i3449 + v[1460];
	v[3561] = 2e0*v[1494];
	v[1496] = v[1463] * v[1495] + i3447 * v[3450];
	v[1497] = -v[1463] + i3449 * v[313];
	v[1499] = v[1463] * v[1498] - i3449 * v[3450];
	v[1501] = v[1463] * v[1500] + i3448 * v[3450];
	v[3451] = 2e0*(-(i3449*v[172]) - v[1463] * v[3331]);
	v[1503] = -(i3448*v[172]) - v[1463] * v[3333];
	v[1504] = -(i3447*v[172]) - v[1463] * v[3334];
	v[1505] = v[1506] * v[314] + i3447 * v[3332];
	v[1507] = v[1506] * v[312] + i3448 * v[3332];
	v[1508] = -(v[1506] * v[311]) - i3449 * v[3332];
	v[1509] = (v[1461] * v[172] - v[1506] * v[356]) / 2e0;
	v[1581] = v[1509] * v[238];
	v[1510] = (-(v[1461] * v[318]) - v[1496] * v[356]) / 2e0;
	v[1511] = (v[1497] * v[172] - v[1506] * v[338]) / 2e0;
	v[1574] = v[1511] * v[236];
	v[1512] = (v[1497] * v[317] + v[1499] * v[338]) / 2e0;
	v[1513] = (v[1464] * v[172] - v[1506] * v[319]) / 2e0;
	v[1566] = v[1513] * v[230];
	v[1514] = (-(v[1464] * v[316]) - v[1501] * v[319]) / 2e0;
	v[1515] = (v[1461] * v[317] + v[3451] + v[1499] * v[356]) / 2e0;
	v[1516] = (v[1464] * v[317] + v[1499] * v[319] + v[3451]) / 2e0;
	v[1517] = v[1503] + v[1461] * v[3333] - v[1501] * v[3426];
	v[1518] = v[1503] + v[1497] * v[3333] - v[1501] * v[3427];
	v[1519] = -v[1506] - v[1491] * v[316] - v[1501] * v[351];
	v[1520] = v[1491] * v[172] - v[1506] * v[351];
	v[1521] = v[1506] + v[1493] * v[317] + v[1499] * v[347];
	v[1522] = v[1493] * v[172] - v[1506] * v[347];
	v[1583] = v[1522] * v[230];
	v[1523] = v[1506] - v[1492] * v[316] - v[1501] * v[343];
	v[1524] = v[1492] * v[172] - v[1506] * v[343];
	v[1575] = v[1524] * v[238];
	v[1525] = v[1504] + v[1497] * v[3334] - v[1496] * v[3427];
	v[1526] = v[1504] + v[1464] * v[3334] - v[1496] * v[3428];
	v[1527] = -v[1506] - v[1489] * v[318] - v[1496] * v[334];
	v[1528] = v[1489] * v[172] - v[1506] * v[334];
	v[1529] = -v[1506] + v[1494] * v[317] + v[1499] * v[330];
	v[1530] = v[1494] * v[172] - v[1506] * v[330];
	v[1567] = v[1530] * v[238];
	v[1531] = -v[1505] + v[1491] * v[317] + v[1499] * v[351];
	v[1532] = -v[1505] - v[1493] * v[316] - v[1501] * v[347];
	v[1533] = -v[1505] + v[1492] * v[317] + v[1499] * v[343];
	v[1534] = -v[1505] - v[1494] * v[316] - v[1501] * v[330];
	v[1535] = v[1545] + v[1505] * v[3452];
	v[1536] = -(v[1515] * v[738]) - v[1531] * v[740] - v[1521] * v[742];
	v[3497] = v[1536] * v[239];
	v[1537] = -(v[1517] * v[738]) - v[1519] * v[740] - v[1532] * v[742];
	v[3498] = v[1537] * v[229];
	v[1538] = v[1506] - v[1490] * v[318] - v[1496] * v[324];
	v[1539] = v[1490] * v[172] - v[1506] * v[324];
	v[1540] = -v[1507] + v[1489] * v[317] + v[1499] * v[334];
	v[1541] = -v[1507] - v[1493] * v[318] - v[1496] * v[347];
	v[1542] = -v[1507] - v[1494] * v[318] - v[1496] * v[330];
	v[1543] = -v[1507] + v[1490] * v[317] + v[1499] * v[324];
	v[1544] = v[1505] * v[326] + v[1507] * v[329];
	v[1546] = v[1545] + v[1507] * v[3453];
	v[1611] = v[1176] * v[1546];
	v[1547] = -(v[1533] * v[738]) - v[1512] * v[740] - v[1540] * v[742];
	v[1853] = v[1547] * v[239];
	v[1548] = v[1508] - v[1491] * v[318] - v[1496] * v[351];
	v[1549] = v[1508] - v[1492] * v[318] - v[1496] * v[343];
	v[1550] = v[1508] - v[1489] * v[316] - v[1501] * v[334];
	v[1551] = v[1508] - v[1490] * v[316] - v[1501] * v[324];
	v[1552] = -(v[1507] * v[323]) - v[1508] * v[326];
	v[1553] = -(v[1505] * v[323]) - v[1508] * v[329];
	v[1554] = v[1545] + v[1508] * v[3454];
	v[1615] = v[1174] * v[1554];
	v[1555] = -(v[1510] * v[738]) - v[1548] * v[740] - v[1541] * v[742];
	v[3499] = v[1555] * v[249];
	v[1556] = -(v[1549] * v[738]) - v[1525] * v[740] - v[1527] * v[742];
	v[1854] = v[1556] * v[249];
	v[1557] = -(v[1523] * v[738]) - v[1518] * v[740] - v[1550] * v[742];
	v[1855] = v[1557] * v[229];
	v[1558] = -(v[1534] * v[738]) - v[1551] * v[740] - v[1514] * v[742];
	v[1559] = -(v[1542] * v[738]) - v[1538] * v[740] - v[1526] * v[742];
	v[1560] = -(v[1193] * v[1510]) - v[1196] * v[1515] - v[1195] * v[1517] - v[1190] * v[1523] - v[1186] * v[1529]
		- v[1194] * v[1533] - v[1188] * v[1534] - v[1187] * v[1542] - v[1191] * v[1549];
	v[1561] = -(v[1194] * v[1512]) - v[1190] * v[1518] - v[1195] * v[1519] - v[1191] * v[1525] - v[1196] * v[1531]
		- v[1187] * v[1538] - v[1186] * v[1543] - v[1193] * v[1548] - v[1188] * v[1551];
	v[1562] = -(v[1529] * v[738]) - v[1543] * v[740] - v[1516] * v[742];
	v[1563] = -(v[1188] * v[1514]) - v[1186] * v[1516] - v[1196] * v[1521] - v[1187] * v[1526] - v[1191] * v[1527]
		- v[1195] * v[1532] - v[1194] * v[1540] - v[1193] * v[1541] - v[1190] * v[1550];
	v[1564] = v[1566] + v[1539] * v[236];
	v[1565] = v[1564] + v[1567] + v[3547];
	v[1568] = v[1566] + v[1567];
	v[1569] = v[1513] * v[88] + v[1539] * v[91] + v[1530] * v[94];
	v[1570] = v[1513] * v[89] + v[1539] * v[92] + v[1530] * v[95];
	v[1571] = v[1513] * v[90] + v[1539] * v[93] + v[1530] * v[96];
	v[1572] = v[1574] + v[1528] * v[230];
	v[1573] = v[1572] + v[1575] + v[3546];
	v[1576] = v[1574] + v[1575];
	v[1577] = v[1528] * v[88] + v[1511] * v[91] + v[1524] * v[94];
	v[1578] = v[1528] * v[89] + v[1511] * v[92] + v[1524] * v[95];
	v[1579] = v[1528] * v[90] + v[1511] * v[93] + v[1524] * v[96];
	v[1580] = v[1581] + v[1583];
	v[1582] = v[1581] + v[1520] * v[236];
	v[1584] = v[1582] + v[1583] + v[3548];
	v[1585] = v[1522] * v[88] + v[1520] * v[91] + v[1509] * v[94];
	v[3506] = v[1569] * v[973] + v[1577] * v[974] + v[1585] * v[975];
	v[1586] = v[1522] * v[89] + v[1520] * v[92] + v[1509] * v[95];
	v[3507] = v[1570] * v[973] + v[1578] * v[974] + v[1586] * v[975];
	v[1587] = v[1522] * v[90] + v[1520] * v[93] + v[1509] * v[96];
	v[3508] = v[1571] * v[973] + v[1579] * v[974] + v[1587] * v[975];
	v[1588] = i3456 * v[10];
	v[3490] = -v[1588] - v[1558] * v[229] - v[1562] * v[239] - v[1559] * v[249];
	v[3489] = v[128] * v[1471] + v[130] * v[1477] + v[132] * v[1483] + v[3490];
	v[1589] = i3459 * v[10];
	v[1590] = i3462 * v[10];
	v[3516] = -v[1590] + v[3467] - v[3497] - v[3498] - v[3499];
	v[1591] = v[1499] + v[1544];
	v[1609] = v[1176] * v[1591];
	v[1592] = -v[1499] + v[1544];
	v[1594] = (v[1591] * v[186] + v[1520] * v[225] + v[1593] * v[3455]) / v[227];
	v[1595] = v[1496] + v[1552];
	v[1617] = v[1176] * v[1595];
	v[1596] = -v[1496] + v[1552];
	v[1613] = v[1174] * v[1596];
	v[1598] = (v[1595] * v[182] + v[1524] * v[226] + v[1597] * v[3455]) / v[227];
	v[1600] = (v[1596] * v[177] + v[1530] * v[3340] + v[1599] * v[3455]) / v[227];
	v[1601] = v[1611] + v[1613];
	v[1602] = v[1501] + v[1553];
	v[1607] = v[1174] * v[1602];
	v[1603] = -v[1501] + v[1553];
	v[1604] = v[1615] + v[1617];
	v[1605] = v[1172] * v[1535] + v[1607] + v[1609];
	v[1608] = v[1605] - v[1609];
	v[1610] = v[1605] - v[1607];
	v[1612] = v[1172] * v[1592] + v[1611];
	v[1614] = v[1612] + v[1613];
	v[1616] = v[1172] * v[1603] + v[1615];
	v[1618] = v[1616] + v[1617];
	v[1619] = -i3456 + i3444 * v[128] + i3457 * v[130] + i3458 * v[132] - v[1569] * v[166] - v[1570] * v[167]
		- v[1571] * v[168];
	v[1621] = -i3459 + i3445 * v[128] + i3460 * v[130] + i3461 * v[132] - v[1577] * v[166] - v[1578] * v[167]
		- v[1579] * v[168];
	v[1623] = -i3462 + i3446 * v[128] + i3463 * v[130] + i3464 * v[132] - v[1585] * v[166] - v[1586] * v[167]
		- v[1587] * v[168];
	v[3465] = v[1619] * v[250] + v[1621] * v[251] + v[1623] * v[252];
	v[1627] = v[3465] / v[602];
	v[1625] = -(v[1053] * v[3465] * v[820]);
	v[1626] = v[1627] * v[263];
	v[1628] = v[1627];
	v[1629] = v[1626] * v[250] + v[1619] * v[264];
	v[1630] = v[1626] * v[251] + v[1621] * v[264];
	v[1631] = v[1626] * v[252] + v[1623] * v[264];
	b1632 = b277;
	if (b1632) {
		v[1633] = -v[1629];
		v[1634] = -v[1630];
		v[1635] = -v[1631];
	}
	else {
		v[1633] = v[1629];
		v[1634] = v[1630];
		v[1635] = v[1631];
	};
	v[3471] = v[1633] * v[280];
	v[3470] = v[1634] * v[279];
	v[3468] = v[1168] * v[1635];
	v[1636] = v[1635] * v[3466];
	v[1638] = v[1635] * v[1637] + v[281] * v[3467];
	v[1639] = 2e0*v[1085] * v[1635] + v[1168] * v[3467];
	v[1643] = i3464 * v[1303] + i3461 * v[1304] + i3458 * v[1305] + v[1168] * (v[1635] * v[220] + v[1487] * v[281]);
	v[1644] = i3463 * v[1303] + i3460 * v[1304] + i3457 * v[1305] + v[1168] * (v[1635] * v[217] + v[1481] * v[281]);
	v[1645] = i3446 * v[1303] + i3445 * v[1304] + i3444 * v[1305] + v[1168] * (v[1635] * v[214] + v[1475] * v[281]);
	v[1646] = v[1634] * v[3469];
	v[1650] = v[1634] * v[1658];
	v[1651] = 2e0*v[1121] * v[1634];
	v[1655] = v[132] * (v[3470] + v[3471]);
	v[1656] = v[130] * (v[3470] + v[3471]);
	v[1657] = v[128] * (v[3470] + v[3471]);
	v[1651] = v[1651] + v[1633] * v[1658];
	v[1662] = v[1633] * v[3472];
	v[1650] = 2e0*v[1156] * v[1633] + v[1650];
	v[1663] = 0e0;
	v[1664] = 0e0;
	v[1665] = 0e0;
	v[1666] = 0e0;
	v[1667] = 0e0;
	v[1668] = 0e0;
	v[1669] = 0e0;
	v[1670] = 0e0;
	v[1671] = 0e0;
	v[1672] = 0e0;
	v[1673] = 0e0;
	v[1674] = 0e0;
	v[1675] = 0e0;
	v[1676] = 0e0;
	v[1677] = 0e0;
	v[1678] = 0e0;
	v[1679] = 0e0;
	v[1680] = 0e0;
	v[1681] = 0e0;
	v[1682] = 0e0;
	v[1683] = 0e0;
	v[1684] = 0e0;
	v[1685] = 0e0;
	v[1686] = 0e0;
	v[1687] = 0e0;
	v[1688] = 0e0;
	v[1689] = 0e0;
	v[1690] = 0e0;
	v[1691] = 0e0;
	v[1692] = 0e0;
	v[1693] = 0e0;
	v[1694] = 0e0;
	b1695 = b542;
	if (b1695) {
		v[1696] = v[1635] * v[275];
		v[1697] = -(v[1635] * v[274]);
		v[1698] = v[1696] - v[1634] * v[276];
		v[1699] = v[1634] * v[274];
		v[1700] = v[1697] + v[1633] * v[276];
		v[1701] = v[1699] - v[1633] * v[275];
		v[3476] = v[1285] * v[1698] + v[1284] * v[1700] + v[1282] * v[1701];
		v[3474] = v[1698] * v[544] + v[1700] * v[545] + v[1701] * v[546];
		v[1702] = v[3474] / v[547];
		v[3475] = v[1702] * v[3769];
		v[1712] = v[1702] * v[3473];
		v[1666] = -(v[1288] * v[1703] * v[3474]);
		v[1704] = v[1698] * v[3362] + v[3475] * v[544];
		v[1719] = v[1704] * v[3529];
		v[1706] = v[1700] * v[3362] + v[3475] * v[545];
		v[1723] = 2e0*v[1706] * v[557];
		v[1707] = v[1701] * v[3362] + v[3475] * v[546];
		v[1720] = v[1707] * v[3526];
		v[1669] = v[1712] * v[3405] + v[3476] * v[549];
		v[1668] = v[1702] * v[1710] * v[553];
		v[1667] = v[1702] * v[3405] * v[554] + v[3476] * v[555];
		v[1693] = v[1712] * v[1713] * v[549];
		v[1694] = v[1287] * v[1702] * v[1708] * v[1713] * v[3770];
		v[1665] = v[1701] * v[3406] + v[1282] * v[3475];
		v[1664] = v[1700] * v[3406] + v[1284] * v[3475];
		v[1663] = v[1698] * v[3406] + v[1285] * v[3475];
		v[1715] = (v[1706] * v[556] + v[1704] * v[557]) / 2e0;
		v[1716] = v[1719] + v[1723];
		v[1717] = v[1716] + v[1720];
		v[1718] = (v[1707] * v[556] + v[1704] * v[558]) / 2e0;
		v[1721] = v[1719] + v[1720];
		v[1722] = (v[1707] * v[557] + v[1706] * v[558]) / 2e0;
		v[1724] = v[1720] + v[1723];
		v[1672] = (v[1268] * v[1704] + v[1255] * v[1706] + 4e0*v[1707] * v[1725]) / 2e0;
		v[1671] = (v[1275] * v[1704] + v[1255] * v[1707] + v[1706] * v[3477]) / 2e0;
		v[1670] = (v[1275] * v[1706] + v[1268] * v[1707] + v[1704] * v[3478]) / 2e0;
		v[1728] = v[1717] * v[3525];
		v[1692] = 8e0*v[1277] * v[1717] * v[3771];
		v[1691] = v[1728] * v[3479];
		v[1729] = v[1707] + v[1715];
		v[1730] = v[1707] - v[1715];
		v[1690] = v[1201] * v[1728];
		v[1731] = -v[1706] + v[1718];
		v[1732] = v[1706] + v[1718];
		v[1689] = v[1202] * v[1728];
		v[1677] = v[1262] * v[1728] + v[1729] * v[559];
		v[1688] = v[1203] * v[1728];
		v[1678] = (v[1258] * v[1728] - v[1721] * v[559]) / 2e0;
		v[1687] = v[1728] * v[3480];
		v[1733] = v[1704] + v[1722];
		v[1734] = -v[1704] + v[1722];
		v[1686] = v[1205] * v[1728];
		v[1680] = v[1249] * v[1728] + v[1731] * v[559];
		v[1685] = v[1206] * v[1728];
		v[1681] = v[1245] * v[1728] + v[1733] * v[559];
		v[1684] = v[1207] * v[1728];
		v[1682] = (v[1241] * v[1728] - v[1716] * v[559]) / 2e0;
		v[1683] = v[1728] * v[3481];
		v[1679] = v[1253] * v[1728] + v[1734] * v[559];
		v[1676] = v[1266] * v[1728] + v[1732] * v[559];
		v[1675] = v[1271] * v[1728] - v[1730] * v[559];
		v[1674] = (v[1276] * v[1728] - v[1724] * v[559]) / 2e0;
		v[1673] = v[1203] * v[1729] - v[1201] * v[1730] + v[1206] * v[1731] + v[1202] * v[1732] + v[1207] * v[1733]
			+ v[1205] * v[1734] - v[1724] * v[3479] - v[1721] * v[3480] - v[1716] * v[3481];
	}
	else {
	};
	v[1735] = 0e0;
	v[1736] = 0e0;
	v[1737] = 0e0;
	v[1738] = 0e0;
	v[1739] = 0e0;
	v[1740] = 0e0;
	b1741 = b30;
	if (b1741) {
		v[3484] = -(v[1633] * v[577]);
		v[3483] = -(v[1634] * v[578]);
		v[3482] = -(v[1635] * v[579]);
		v[1778] = v[1049] * v[1633];
		v[1749] = -i3462 + i3446 * v[133] + i3463 * v[135] + i3464 * v[137] - v[1585] * v[169] - v[1586] * v[170]
			- v[1587] * v[171] + v[1682] * v[3];
		v[1682] = 0e0;
		v[1750] = v[1749] + v[1681] * v[2];
		v[1681] = 0e0;
		v[1751] = v[1] * v[1680] + v[1750];
		v[3486] = -(v[1751] * v[281]);
		v[1680] = 0e0;
		v[1759] = -i3459 + i3445 * v[133] + i3460 * v[135] + i3461 * v[137] - v[1577] * v[169] - v[1578] * v[170]
			- v[1579] * v[171] + v[1679] * v[3];
		v[1679] = 0e0;
		v[1760] = v[1759] + v[1678] * v[2];
		v[1678] = 0e0;
		v[1761] = v[1] * v[1677] + v[1760];
		v[3485] = -(v[1761] * v[280]);
		v[1677] = 0e0;
		v[1769] = -i3456 + i3444 * v[133] + i3457 * v[135] + i3458 * v[137] - v[1569] * v[169] - v[1570] * v[170]
			- v[1571] * v[171] + v[1676] * v[3];
		v[1676] = 0e0;
		v[1770] = v[1769] + v[1675] * v[2];
		v[1675] = 0e0;
		v[1771] = v[1] * v[1674] + v[1770];
		v[3487] = -(v[1771] * v[279]);
		v[1674] = 0e0;
		v[1772] = v[1047] * v[1635];
		v[1737] = v[1635] * v[1782];
		v[1650] = v[1650] + v[1049] * v[3482];
		v[1651] = v[1651] + v[1048] * v[3482];
		v[1774] = v[1048] * v[1634];
		v[1775] = v[1772] + v[1774];
		v[1736] = v[1634] * v[1784];
		v[1650] = v[1650] + v[1049] * v[3483];
		v[1639] = v[1639] + v[1047] * v[3483];
		v[1777] = v[1774] + v[1778];
		v[1779] = v[1772] + v[1778];
		v[1735] = v[1633] * v[1783];
		v[1651] = v[1651] + v[1048] * v[3484];
		v[1639] = v[1639] + v[1047] * v[3484];
		v[1735] = -(v[1049] * v[1662]) + v[1735];
		v[1740] = v[1047] * v[1751];
		v[1739] = v[1048] * v[1761];
		v[1738] = v[1049] * v[1771];
		v[1736] = -(v[1048] * v[1646]) + v[1736];
		v[1737] = -(v[1047] * v[1636]) + v[1737];
		v[1737] = v[1737] - v[1777] * v[281];
		v[1639] = v[1639] + v[1751] * v[1782] + v[1047] * (v[3485] + v[3487]) - v[1777] * v[579];
		v[1735] = v[1735] - v[1775] * v[279];
		v[1650] = v[1650] + v[1771] * v[1783] + v[1049] * (v[3485] + v[3486]) - v[1775] * v[577];
		v[1736] = v[1736] - v[1779] * v[280];
		v[1651] = v[1651] + v[1761] * v[1784] + v[1048] * (v[3486] + v[3487]) - v[1779] * v[578];
	}
	else {
	};
	v[3541] = -(v[1052] * v[1636]);
	v[3496] = -(v[1051] * v[1646]);
	v[3495] = -(v[1126] * v[1633]);
	v[1816] = -(v[1634] * v[3488]);
	v[3494] = v[1052] * v[1635];
	v[3491] = v[1635] * v[1804];
	v[1806] = -(v[1635] * v[3488]);
	v[3493] = -(v[1126] * v[1635]) + v[1806];
	v[1785] = v[1590] + v[3497] + v[3498] + v[3499];
	v[3501] = -(v[1785] * v[281]);
	v[1786] = v[1589] + v[1853] + v[1854] + v[1855];
	v[3502] = v[1786] * v[280];
	v[1791] = v[1545] * v[1787] + v[1546] * v[1788] + v[1591] * v[1789] + v[1595] * v[1790] + v[1594] * v[236]
		+ v[1598] * v[238] + v[1565] * v[2507] + v[1572] * v[3341] + v[1580] * v[3344] + v[5785 + i971];
	v[1796] = v[1438] * v[1564] + v[1437] * v[1582] + v[1545] * v[1792] + v[1554] * v[1793] + v[1596] * v[1794]
		+ v[1602] * v[1795] + v[1594] * v[230] + v[1600] * v[238] + v[1573] * v[2510] + v[5815 + i971];
	v[1801] = v[1440] * v[1568] + v[1439] * v[1576] + v[1545] * v[1797] + v[1535] * v[1798] + v[1592] * v[1799]
		+ v[1603] * v[1800] + v[1598] * v[230] + v[1600] * v[236] + v[1584] * v[2512] + v[5845 + i971];
	v[1803] = v[1175] * v[1791] + v[1173] * v[1796] + v[1171] * v[1801] + v[1635] * v[1802] + v[3489] * v[582]
		+ v[1473] * v[584] + v[1479] * v[585] + v[1485] * v[586];
	v[1650] = v[1650] + v[1050] * v[3491];
	v[1803] = v[1803] + v[1634] * v[1813] + v[1638] * v[279];
	v[1650] = v[1650] + v[1050] * (v[1638] + v[1634] * v[1814]);
	v[1803] = v[1803] + v[1633] * v[1823];
	v[1824] = v[1050] * v[1633];
	v[1825] = v[1824] * v[229];
	v[1826] = v[1824] * v[239];
	v[1827] = v[1824] * v[249];
	v[1803] = v[1803] + v[1662] * v[1838] + v[1657] * v[213] + v[1656] * v[216] + v[1655] * v[219];
	v[1842] = -(v[1050] * v[1662]);
	v[3503] = -(v[1050] * v[1785]) + v[1052] * v[3489];
	v[1856] = v[128] * v[1473] + v[130] * v[1479] + v[132] * v[1485] - v[1589] - v[1853] - v[1854] - v[1855];
	v[3504] = v[1856] * v[280];
	v[1857] = v[1638] + v[279] * v[3490] + v[3501];
	v[1859] = v[1133] * v[1791] + v[1131] * v[1796] + v[1129] * v[1801] + v[1635] * v[1858] + v[1471] * v[584]
		+ v[1477] * v[585] + v[1483] * v[586] + v[1856] * v[587];
	v[1651] = v[1651] + v[1051] * v[3491];
	v[1921] = v[249] * v[3494];
	v[1919] = v[239] * v[3494];
	v[1917] = v[229] * v[3494];
	v[1859] = v[1859] + v[1634] * v[1871] + v[1857] * v[280];
	v[1872] = v[1051] * v[1634];
	v[1873] = v[1872] * v[229];
	v[1874] = v[1872] * v[239];
	v[1875] = v[1872] * v[249];
	v[1876] = v[1824] + v[1872];
	v[1877] = v[1825] + v[1873];
	v[1878] = v[1826] + v[1874];
	v[1879] = v[1827] + v[1875];
	v[1859] = v[1859] + v[1646] * v[1886];
	v[1859] = v[1859] + v[1633] * v[1900];
	v[1651] = v[1651] + v[1051] * (v[1857] + v[1633] * v[1901]);
	v[1905] = v[1842] + v[3495];
	v[1859] = v[1859] + v[1657] * v[212] + v[1656] * v[215] + v[1655] * v[218];
	v[1993] = v[1859];
	v[1910] = v[1090] * v[1791] + v[1088] * v[1796] + v[1086] * v[1801] + v[1635] * v[1909] + v[3516] * v[592];
	v[1911] = v[1872] + v[3494];
	v[1912] = v[1873] + v[1917];
	v[1913] = v[1874] + v[1919];
	v[1914] = v[1875] + v[1921];
	v[1650] = v[1650] + v[1126] * v[3490] + v[2907] * v[3494];
	v[1916] = v[1824] + v[3494];
	v[1918] = v[1825] + v[1917];
	v[1920] = v[1826] + v[1919];
	v[1922] = v[1827] + v[1921];
	v[1651] = v[1651] + v[2905] * v[3494];
	v[1910] = v[1910] + v[1636] * v[1938];
	v[1910] = v[1910] + v[1634] * v[1952];
	v[1953] = v[1052] * v[1634];
	v[3582] = -(v[1953] * v[281]);
	v[3536] = v[1816] + v[3496] + v[3582];
	v[1639] = v[1639] - v[1126] * v[1785] + v[1814] * v[1953];
	v[1910] = v[1910] + v[1633] * v[1963];
	v[1964] = v[1052] * v[1633];
	v[3533] = -(v[1964] * v[281]);
	v[1639] = v[1639] + v[1901] * v[1964];
	v[1974] = 0e0;
	b1975 = b603;
	if (b1975) {
		b1976 = b605;
		if (b1976) {
			v[1974] = -v[1628];
			v[1633] = 0e0;
			v[1634] = 0e0;
			v[1635] = 0e0;
		}
		else {
		};
	}
	else {
	};
	v[1910] = v[1910] + v[281] * (v[279] * v[3489] + v[3504]);
	v[1994] = v[1910];
	v[1651] = v[1651] - v[1786] * v[3488] + v[1856] * v[3500];
	v[2171] = v[1651];
	v[1803] = v[1803] - v[279] * (-v[3501] + v[3502]);
	v[1992] = v[1803];
	v[1650] = v[1650] - v[1050] * v[3502] + v[281] * v[3503];
	v[2173] = v[1650];
	v[1639] = v[1639] + v[279] * v[3503] + v[1052] * v[3504];
	v[2169] = v[1639];
	v[1977] = 0e0;
	v[1978] = 0e0;
	v[1979] = 0e0;
	v[1980] = 0e0;
	v[1981] = 0e0;
	v[1982] = 0e0;
	v[1983] = 0e0;
	v[1984] = 0e0;
	b1985 = b603;
	if (b1985) {
		v[3505] = v[1056] * v[1992];
		v[1995] = v[1057] * v[1993] + v[1058] * v[1994] + v[3505];
		b1986 = b605;
		if (b1986) {
			v[2154] = Power(v[614], v[1987]);
			v[2153] = (v[1064] * v[615]) / v[1065];
			v[1989] = v[1974] * v[2153] * v[7];
			v[1988] = v[1989] * v[2154];
			v[1983] = -((v[1063] * v[1988]) / v[1065]);
			v[1980] = v[1063] * v[1987] * v[1989] * Power(v[614], -2e0 + v[615]);
			v[1974] = 0e0;
			v[1803] = 0e0;
			v[1859] = 0e0;
			v[1981] = v[1995];
			v[1910] = 0e0;
		}
		else {
			v[2162] = Power(v[602], v[1990]);
			v[2161] = (v[1072] * v[620]) / v[1073];
			v[1991] = v[1628] * v[2161] * v[7];
			v[1988] = v[1991] * v[2162];
			v[1984] = -((v[1063] * v[1988]) / v[1073]);
			v[1625] = v[1625] + v[1063] * v[1990] * v[1991] * Power(v[602], -2e0 + v[620]);
			v[1628] = 0e0;
			v[1803] = 0e0;
			v[1859] = 0e0;
			v[1982] = v[1995];
			v[1910] = 0e0;
		};
		v[1977] = v[1056] * v[1988];
		v[1978] = v[1057] * v[1988];
		v[1979] = v[1058] * v[1988];
	}
	else {
	};
	v[2160] = v[1977];
	v[2158] = v[1978];
	v[2156] = v[1979];
	v[1997] = v[9] * ((i3457 / 2e0 + v[1998])*v[973] + (i3460 / 2e0 + v[1999])*v[974] + (i3463 / 2e0 + v[2000])*v[975]);
	v[2001] = v[9] * ((i3458 / 2e0 + v[1998])*v[973] + (i3461 / 2e0 + v[1999])*v[974] + (i3464 / 2e0 + v[2000])*v[975]);
	v[2002] = (v[301] * v[3506] + v[303] * v[3507] + v[305] * v[3508])*v[9];
	v[2003] = (v[293] * v[3506] + v[294] * v[3507] + v[295] * v[3508])*v[9];
	v[2004] = -(i3444*v[910]) - i3445 * v[914] - i3446 * v[918] - i3457 * v[922] - i3460 * v[926] - i3463 * v[930]
		- i3458 * v[934] - i3461 * v[938] - i3464 * v[942] - i3456 * v[946] - i3459 * v[950] - i3462 * v[954] - i3448 * v[958]
		- i3449 * v[962] - i3447 * v[966];
	v[2020] = v[2004];
	v[2005] = -(i3444*v[911]) - i3445 * v[915] - i3446 * v[919] - i3457 * v[923] - i3460 * v[927] - i3463 * v[931]
		- i3458 * v[935] - i3461 * v[939] - i3464 * v[943] - i3456 * v[947] - i3459 * v[951] - i3462 * v[955] - i3448 * v[959]
		- i3449 * v[963] - i3447 * v[967];
	v[2018] = v[2005];
	v[2006] = -(i3444*v[912]) - i3445 * v[916] - i3446 * v[920] - i3457 * v[924] - i3460 * v[928] - i3463 * v[932]
		- i3458 * v[936] - i3461 * v[940] - i3464 * v[944] - i3456 * v[948] - i3459 * v[952] - i3462 * v[956] - i3448 * v[960]
		- i3449 * v[964] - i3447 * v[968];
	v[2016] = v[2006];
	v[2007] = 0e0;
	v[2008] = 0e0;
	v[2009] = 0e0;
	v[2010] = 0e0;
	v[2011] = 0e0;
	v[2012] = 0e0;
	b2013 = b603;
	if (b2013) {
		b2014 = b29;
		if (b2014) {
			b2015 = b641;
			if (b2015) {
				v[2012] = v[2006];
				v[2006] = 0e0;
				v[2011] = v[2005];
				v[2005] = 0e0;
				v[2010] = v[2004];
				v[2004] = 0e0;
			}
			else {
				v[2028] = v[2020] * v[3371];
				v[2026] = v[2018] * v[3371];
				v[2024] = v[2016] * v[3371];
				v[2006] = 0e0;
				v[2005] = 0e0;
				v[3510] = (v[6] * (v[2020] * v[653] + v[2018] * v[654] + v[2016] * v[655])) / sqrt(v[2029]);
				v[2004] = 0e0;
				v[3509] = ((v[2028] * v[636] + v[2026] * v[637] + v[2024] * v[638])*v[651]) / sqrt(v[2025]);
				v[2012] = v[3509] * v[638] + v[2024] * v[652];
				v[2011] = v[3509] * v[637] + v[2026] * v[652];
				v[2010] = v[3509] * v[636] + v[2028] * v[652];
				v[2009] = v[3510] * v[629];
				v[2008] = v[3510] * v[628];
				v[2007] = v[3510] * v[627];
			};
		}
		else {
			b2031 = b675;
			if (b2031) {
				v[2012] = v[2016];
				v[2006] = 0e0;
				v[2011] = v[2018];
				v[2005] = 0e0;
				v[2010] = v[2020];
				v[2004] = 0e0;
			}
			else {
				v[2037] = v[2020] * v[6] * v[686];
				v[2035] = v[2018] * v[6] * v[686];
				v[2034] = v[2016] * v[6] * v[686];
				v[2006] = 0e0;
				v[2005] = 0e0;
				v[3512] = (v[6] * (v[2020] * v[683] + v[2018] * v[684] + v[2016] * v[685])) / sqrt(v[2029]);
				v[2004] = 0e0;
				v[3511] = ((v[2037] * v[636] + v[2035] * v[637] + v[2034] * v[638])*v[681]) / sqrt(v[2025]);
				v[2012] = v[3511] * v[638] + v[2034] * v[682];
				v[2011] = v[3511] * v[637] + v[2035] * v[682];
				v[2010] = v[3511] * v[636] + v[2037] * v[682];
				v[2009] = v[3512] * v[629];
				v[2008] = v[3512] * v[628];
				v[2007] = v[3512] * v[627];
			};
		};
	}
	else {
	};
	v[2043] = -(i3447*v[645]) - v[9] * (v[2012] * v[249] + v[1801] * v[975]);
	v[2044] = -(i3449*v[645]) - v[9] * (v[2012] * v[239] + v[1796] * v[975]);
	v[2045] = -(i3448*v[645]) - v[9] * (v[2012] * v[229] + v[1791] * v[975]);
	v[2046] = -(i3462*v[645]) - v[9] * (v[2012] * v[223] + v[1590] * v[975]);
	v[2047] = -(i3459*v[645]) - v[9] * (v[2012] * v[222] + v[1589] * v[975]);
	v[2048] = -(i3456*v[645]) - v[9] * (v[2012] * v[221] + v[1588] * v[975]);
	v[2049] = -(i3464*v[645]) - v[9] * (v[2012] * v[220] + v[1487] * v[975]);
	v[2050] = -(i3461*v[645]) - v[9] * (v[2012] * v[219] + v[1485] * v[975]);
	v[2051] = -(i3458*v[645]) - v[9] * (v[2012] * v[218] + v[1483] * v[975]);
	v[2052] = -(i3463*v[645]) - v[9] * (v[2012] * v[217] + v[1481] * v[975]);
	v[2053] = -(i3460*v[645]) - v[9] * (v[2012] * v[216] + v[1479] * v[975]);
	v[2054] = -(i3457*v[645]) - v[9] * (v[2012] * v[215] + v[1477] * v[975]);
	v[2055] = -(i3446*v[645]) - v[9] * (v[2012] * v[214] + v[1475] * v[975]);
	v[2056] = -(i3445*v[645]) - v[9] * (v[2012] * v[213] + v[1473] * v[975]);
	v[2057] = -(i3444*v[645]) - v[9] * (v[2012] * v[212] + v[1471] * v[975]);
	v[2074] = -(i3447*v[644]) - v[9] * (v[2011] * v[249] + v[1801] * v[974]);
	v[2075] = -(i3449*v[644]) - v[9] * (v[2011] * v[239] + v[1796] * v[974]);
	v[2076] = -(i3448*v[644]) - v[9] * (v[2011] * v[229] + v[1791] * v[974]);
	v[2077] = -(i3462*v[644]) - v[9] * (v[2011] * v[223] + v[1590] * v[974]);
	v[2078] = -(i3459*v[644]) - v[9] * (v[2011] * v[222] + v[1589] * v[974]);
	v[2079] = -(i3456*v[644]) - v[9] * (v[2011] * v[221] + v[1588] * v[974]);
	v[2080] = -(i3464*v[644]) - v[9] * (v[2011] * v[220] + v[1487] * v[974]);
	v[2081] = -(i3461*v[644]) - v[9] * (v[2011] * v[219] + v[1485] * v[974]);
	v[2082] = -(i3458*v[644]) - v[9] * (v[2011] * v[218] + v[1483] * v[974]);
	v[2083] = -(i3463*v[644]) - v[9] * (v[2011] * v[217] + v[1481] * v[974]);
	v[2084] = -(i3460*v[644]) - v[9] * (v[2011] * v[216] + v[1479] * v[974]);
	v[2085] = -(i3457*v[644]) - v[9] * (v[2011] * v[215] + v[1477] * v[974]);
	v[2086] = -(i3446*v[644]) - v[9] * (v[2011] * v[214] + v[1475] * v[974]);
	v[2087] = -(i3445*v[644]) - v[9] * (v[2011] * v[213] + v[1473] * v[974]);
	v[2088] = -(i3444*v[644]) - v[9] * (v[2011] * v[212] + v[1471] * v[974]);
	v[2105] = -(i3447*v[643]) - v[9] * (v[2010] * v[249] + v[1801] * v[973]);
	v[2106] = -(i3449*v[643]) - v[9] * (v[2010] * v[239] + v[1796] * v[973]);
	v[2107] = -(i3448*v[643]) - v[9] * (v[2010] * v[229] + v[1791] * v[973]);
	v[2108] = -(i3462*v[643]) - v[9] * (v[2010] * v[223] + v[1590] * v[973]);
	v[2109] = -(i3459*v[643]) - v[9] * (v[2010] * v[222] + v[1589] * v[973]);
	v[2110] = -(i3456*v[643]) - v[9] * (v[2010] * v[221] + v[1588] * v[973]);
	v[2111] = -(i3464*v[643]) - v[9] * (v[2010] * v[220] + v[1487] * v[973]);
	v[2112] = -(i3461*v[643]) - v[9] * (v[2010] * v[219] + v[1485] * v[973]);
	v[2113] = -(i3458*v[643]) - v[9] * (v[2010] * v[218] + v[1483] * v[973]);
	v[2114] = -(i3463*v[643]) - v[9] * (v[2010] * v[217] + v[1481] * v[973]);
	v[2115] = -(i3460*v[643]) - v[9] * (v[2010] * v[216] + v[1479] * v[973]);
	v[2116] = -(i3457*v[643]) - v[9] * (v[2010] * v[215] + v[1477] * v[973]);
	v[2117] = -(i3446*v[643]) - v[9] * (v[2010] * v[214] + v[1475] * v[973]);
	v[2118] = -(i3445*v[643]) - v[9] * (v[2010] * v[213] + v[1473] * v[973]);
	v[2119] = -(i3444*v[643]) - v[9] * (v[2010] * v[212] + v[1471] * v[973]);
	v[2135] = v[2012] * v[4];
	v[2136] = v[2011] * v[4];
	v[2137] = v[2010] * v[4];
	v[2138] = v[2009];
	v[2139] = v[2009];
	v[2140] = v[2008];
	v[2141] = v[2008];
	v[2142] = v[2007];
	v[2143] = v[2007];
	b2144 = b603;
	if (b2144) {
		v[2145] = 0e0;
		v[2146] = 0e0;
		v[2147] = 0e0;
		b2148 = b622;
		if (b2148) {
			v[2147] = v[2138];
			v[2138] = 0e0;
			v[2146] = v[2140];
			v[2140] = 0e0;
			v[2145] = v[2142];
			v[2142] = 0e0;
		}
		else {
			v[2139] = 0e0;
			v[2138] = 0e0;
			v[2141] = 0e0;
			v[2140] = 0e0;
			v[2143] = 0e0;
			v[2142] = 0e0;
		};
		v[2159] = v[2145] * v[583];
		v[2157] = v[2146] * v[591];
		v[2155] = v[2147] * v[593];
		b2152 = b605;
		if (b2152) {
			v[1981] = v[1981] + v[2155];
			v[1979] = v[1979] + v[2147] * v[617];
			v[1981] = v[1981] + v[2157];
			v[1978] = v[1978] + v[2146] * v[617];
			v[1981] = v[1981] + v[2159];
			v[1977] = v[1977] + v[2145] * v[617];
			v[1983] = v[1983] + v[1981] * v[3370];
			v[1980] = v[1980] + (v[1983] * v[2153] * v[2154]) / 2e0;
		}
		else {
			v[1982] = v[1982] + v[2155];
			v[1979] = v[2156] + v[2147] * v[621];
			v[1982] = v[1982] + v[2157];
			v[1978] = v[2158] + v[2146] * v[621];
			v[1982] = v[1982] + v[2159];
			v[1977] = v[2160] + v[2145] * v[621];
			v[1984] = v[1984] + v[1982] * v[3370];
			v[1625] = v[1625] + (v[1984] * v[2161] * v[2162]) / 2e0;
		};
	}
	else {
	};
	v[3580] = v[1977] * v[582];
	v[3568] = -v[1842] + v[3580];
	v[3518] = -(v[1977] * v[279]);
	v[3569] = -v[3496] + v[1978] * v[587];
	v[3517] = v[1978] * v[280];
	v[3581] = v[1979] * v[281];
	v[3570] = v[1979] * v[592];
	v[2174] = v[1625];
	v[2172] = v[2143];
	v[2170] = v[2141];
	v[2168] = v[2139];
	b2163 = b603;
	if (b2163) {
		v[3515] = v[2143] * v[279];
		v[3514] = v[2141] * v[280];
		v[3513] = v[2139] * v[281];
		b2164 = b605;
		if (b2164) {
			v[1639] = v[1639] + v[2139] * v[3368];
			v[2139] = 0e0;
			v[1651] = v[1651] + v[2141] * v[3368];
			v[2141] = 0e0;
			v[1650] = v[1650] + v[2143] * v[3368];
			v[2143] = 0e0;
			v[1980] = v[1980] + v[31] * v[34] * (-v[3513] - v[3514] - v[3515])*Power(v[614], v[615]);
			v[1625] = v[1625] - v[1980];
		}
		else {
			v[1639] = v[2169] + v[2168] * v[611];
			v[2139] = 0e0;
			v[1651] = v[2171] + v[2170] * v[611];
			v[2141] = 0e0;
			v[1650] = v[2173] + v[2172] * v[611];
			v[2143] = 0e0;
			v[1625] = v[2174] + v[3382] * (v[3513] + v[3514] + v[3515])*Power(v[602], v[620]);
		};
	}
	else {
	};
	v[2175] = v[1052] * v[1801] + v[1979] * v[249];
	v[3538] = v[2175] * v[281];
	v[2176] = v[1052] * v[1796] + v[1979] * v[239];
	v[3537] = v[2176] * v[281];
	v[2177] = v[1052] * v[1791] + v[1979] * v[229];
	v[3535] = v[2177] * v[281];
	v[1650] = v[1650] + v[1979] * v[2178];
	v[1651] = v[1651] + v[1979] * v[2179];
	v[2180] = v[1979] * v[2542] + v[1052] * v[3516];
	v[2181] = v[1964] + v[1979] * v[279];
	v[3577] = v[2181] * v[281] + v[279] * v[3494] + v[3568];
	v[2184] = v[1953] + v[1979] * v[280];
	v[3578] = v[2184] * v[281] + v[280] * v[3494] + v[3569];
	v[1639] = v[1639] + v[1979] * v[2550];
	v[2202] = v[1051] * v[1801] + v[1978] * v[249];
	v[3544] = v[2202] * v[280];
	v[2203] = v[1051] * v[1796] + v[1978] * v[239];
	v[3542] = v[2203] * v[280];
	v[2204] = v[1051] * v[1791] + v[1978] * v[229];
	v[3539] = v[2204] * v[280];
	v[1650] = v[1650] - v[221] * v[3517];
	v[2205] = v[1051] * v[1856] + v[1978] * v[2569];
	v[1651] = v[1651] + v[1978] * v[2571];
	v[1639] = v[1639] - v[223] * v[3517];
	v[2224] = v[1050] * v[1801] + v[1977] * v[249];
	v[3545] = v[2224] * v[279];
	v[2225] = v[1050] * v[1796] + v[1977] * v[239];
	v[3543] = v[2225] * v[279];
	v[2226] = v[1050] * v[1791] + v[1977] * v[229];
	v[3540] = v[2226] * v[279];
	v[2227] = v[1977] * v[2593] + v[1050] * v[3489];
	v[2228] = v[1051] * v[1471] + v[1050] * v[1473] + v[1978] * v[212] + v[1977] * v[213];
	v[2229] = v[1051] * v[1477] + v[1050] * v[1479] + v[1978] * v[215] + v[1977] * v[216];
	v[2230] = v[1051] * v[1483] + v[1050] * v[1485] + v[1978] * v[218] + v[1977] * v[219];
	v[3534] = v[128] * v[2228] + v[130] * v[2229] + v[132] * v[2230];
	v[1650] = v[1650] + v[1977] * v[2598];
	v[1651] = v[1651] + v[222] * v[3518];
	v[1639] = v[1639] + v[223] * v[3518];
	v[2249] = 0e0;
	v[2250] = 0e0;
	v[2251] = 0e0;
	v[2252] = 0e0;
	v[2253] = 0e0;
	v[2254] = 0e0;
	v[2255] = 0e0;
	v[2256] = 0e0;
	v[2257] = 0e0;
	b2258 = b30;
	if (b2258) {
		v[3521] = v[2137] * v[279];
		v[3520] = v[2136] * v[280];
		v[3524] = v[3520] + v[3521];
		v[3519] = v[2135] * v[281];
		v[3523] = v[3519] + v[3521];
		v[3522] = v[3519] + v[3520];
		v[1740] = v[1740] + v[2135] * v[579];
		v[1739] = v[1739] + v[2136] * v[578];
		v[1735] = v[1735] + v[1210] * v[2137] - v[279] * v[3522];
		v[1736] = v[1736] + v[1212] * v[2136] - v[280] * v[3523];
		v[1737] = v[1737] + v[1214] * v[2135] - v[281] * v[3524];
		v[1738] = v[1738] + v[2137] * v[577];
		v[1650] = v[1650] - v[2137] * v[3401] - v[3522] * v[577];
		v[1651] = v[1651] + v[2136] * v[3402] - v[3523] * v[578];
		v[1639] = v[1639] + v[2135] * v[3403] - v[3524] * v[579];
		v[2249] = v[1] * v[1735];
		v[2250] = v[1735] * v[2];
		v[2251] = v[1735] * v[3];
		v[2259] = -v[1735];
		v[2260] = -(v[171] * v[1735]);
		v[2261] = -(v[170] * v[1735]);
		v[2262] = -(v[169] * v[1735]);
		v[2263] = v[137] * v[1735];
		v[2264] = v[135] * v[1735];
		v[2265] = v[133] * v[1735];
		v[2252] = v[1] * v[1736];
		v[2253] = v[1736] * v[2];
		v[2254] = v[1736] * v[3];
		v[2266] = -v[1736];
		v[2267] = -(v[171] * v[1736]);
		v[2268] = -(v[170] * v[1736]);
		v[2269] = -(v[169] * v[1736]);
		v[2270] = v[137] * v[1736];
		v[2271] = v[135] * v[1736];
		v[2272] = v[133] * v[1736];
		v[2255] = v[1] * v[1737];
		v[2256] = v[1737] * v[2];
		v[2257] = v[1737] * v[3];
		v[2273] = -v[1737];
		v[2274] = -(v[171] * v[1737]);
		v[2275] = -(v[170] * v[1737]);
		v[2276] = -(v[169] * v[1737]);
		v[2277] = v[137] * v[1737];
		v[2278] = v[135] * v[1737];
		v[2279] = v[133] * v[1737];
		v[2227] = -v[1738] + v[2227];
		v[2205] = -v[1739] + v[2205];
		v[2180] = -v[1740] + v[2180];
	}
	else {
		v[2265] = 0e0;
		v[2264] = 0e0;
		v[2263] = 0e0;
		v[2272] = 0e0;
		v[2271] = 0e0;
		v[2270] = 0e0;
		v[2279] = 0e0;
		v[2278] = 0e0;
		v[2277] = 0e0;
		v[2262] = 0e0;
		v[2261] = 0e0;
		v[2260] = 0e0;
		v[2269] = 0e0;
		v[2268] = 0e0;
		v[2267] = 0e0;
		v[2276] = 0e0;
		v[2275] = 0e0;
		v[2274] = 0e0;
		v[2259] = 0e0;
		v[2266] = 0e0;
		v[2273] = 0e0;
	};
	b2280 = b542;
	if (b2280) {
		v[1673] = v[1673] + (v[1241] * v[2257]) / 2e0;
		v[1683] = v[1683] + v[2257] * v[3364];
		v[1673] = v[1673] + v[1245] * v[2256];
		v[1684] = v[1684] + v[2256] * v[559];
		v[1673] = v[1673] + v[1249] * v[2255];
		v[1685] = v[1685] + v[2255] * v[559];
		v[1673] = v[1673] + v[1253] * v[2254];
		v[1686] = v[1686] + v[2254] * v[559];
		v[1673] = v[1673] + (v[1258] * v[2253]) / 2e0;
		v[1687] = v[1687] + v[2253] * v[3364];
		v[1673] = v[1673] + v[1262] * v[2252];
		v[1688] = v[1688] + v[2252] * v[559];
		v[1673] = v[1673] + v[1266] * v[2251];
		v[1689] = v[1689] + v[2251] * v[559];
		v[1673] = v[1673] + v[1271] * v[2250];
		v[1690] = v[1690] + v[2250] * v[559];
		v[1673] = v[1673] + (v[1276] * v[2249]) / 2e0;
		v[1691] = v[1691] + v[2249] * v[3364];
		v[1692] = v[1692] + v[1673] * v[3525];
		v[3530] = -v[1687] + v[1692];
		v[1671] = v[1671] - v[1685];
		v[2290] = v[1685] + v[1689];
		v[1671] = v[1671] + v[1689];
		v[1670] = v[1670] + v[1684] + (v[2290] * v[558]) / 2e0;
		v[2292] = v[1684] + v[1686];
		v[1670] = v[1670] - v[1686];
		v[1672] = v[1672] + v[1688] + v[2292] * v[3363] + v[2290] * v[3404] + v[3526] * (-v[1691] + v[3530]);
		v[1672] = v[1672] - v[1690];
		v[3527] = v[1672] * v[546];
		v[2294] = v[1688] + v[1690];
		v[1669] = v[1669] + v[3527] * v[549];
		v[1667] = v[1667] + v[3527] * v[555];
		v[1665] = v[1665] + v[1672] * v[3362];
		v[1671] = v[1671] + (v[2294] * v[556] - 4e0*(v[1683] + v[1691] - v[1692])*v[557] + v[2292] * v[558]) / 2e0;
		v[3528] = v[1671] * v[545];
		v[1669] = v[1669] + v[3528] * v[549];
		v[1667] = v[1667] + v[3528] * v[555];
		v[1664] = v[1664] + v[1671] * v[3362];
		v[1670] = v[1670] + v[2294] * v[3363] + v[3529] * (-v[1683] + v[3530]);
		v[3531] = v[1670] * v[544];
		v[1669] = v[1669] + v[3531] * v[549];
		v[1667] = v[1667] + v[3531] * v[555];
		v[1663] = v[1663] + v[1670] * v[3362];
		v[1668] = v[1668] + v[1669] * v[554];
		v[1666] = v[1666] + v[1668];
		v[1693] = v[1693] + 2e0*v[1667] * v[1708];
		v[1694] = v[1694] + (v[1693] * v[1709]) / 2e0;
		v[1666] = v[1666] + v[1694];
		v[3532] = v[1666] / v[547];
		v[1665] = v[1665] + v[3532] * v[546];
		v[1664] = v[1664] + v[3532] * v[545];
		v[1663] = v[1663] + v[3532] * v[544];
		v[1650] = v[1650] - v[1665] * v[275];
		v[1651] = v[1651] + v[1665] * v[274];
		v[1650] = v[1650] + v[1664] * v[276];
		v[1639] = v[1639] - v[1664] * v[274];
		v[1651] = v[1651] - v[1663] * v[276];
		v[1639] = v[1639] + v[1663] * v[275];
	}
	else {
	};
	v[2300] = -(v[2119] * v[512]) - v[2118] * v[513] - v[2117] * v[514] - v[2116] * v[515] - v[2115] * v[516]
		- v[2114] * v[517] - v[2113] * v[518] - v[2112] * v[519] - v[2111] * v[520] - v[2110] * v[521] - v[2109] * v[522]
		- v[2108] * v[523] - v[2107] * v[524] - v[2106] * v[525] - v[2105] * v[526];
	v[2301] = -(v[2119] * v[527]) - v[2118] * v[528] - v[2117] * v[529] - v[2116] * v[530] - v[2115] * v[531]
		- v[2114] * v[532] - v[2113] * v[533] - v[2112] * v[534] - v[2111] * v[535] - v[2110] * v[536] - v[2109] * v[537]
		- v[2108] * v[538] - v[2107] * v[539] - v[2106] * v[540] - v[2105] * v[541];
	v[3575] = (v[2119] * v[497] + v[2118] * v[498] + v[2117] * v[499] + v[2116] * v[500] + v[2115] * v[501] + v[2114] * v[502]
		+ v[2113] * v[503] + v[2112] * v[504] + v[2111] * v[505] + v[2110] * v[506] + v[2109] * v[507] + v[2108] * v[508]
		+ v[2107] * v[509] + v[2106] * v[510] + v[2105] * v[511]) / 2e0;
	v[3576] = (v[2119] * v[482] + v[2118] * v[483] + v[2117] * v[484] + v[2116] * v[485] + v[2115] * v[486] + v[2114] * v[487]
		+ v[2113] * v[488] + v[2112] * v[489] + v[2111] * v[490] + v[2110] * v[491] + v[2109] * v[492] + v[2108] * v[493]
		+ v[2107] * v[494] + v[2106] * v[495] + v[2105] * v[496]) / 2e0;
	v[2304] = -(v[2088] * v[512]) - v[2087] * v[513] - v[2086] * v[514] - v[2085] * v[515] - v[2084] * v[516]
		- v[2083] * v[517] - v[2082] * v[518] - v[2081] * v[519] - v[2080] * v[520] - v[2079] * v[521] - v[2078] * v[522]
		- v[2077] * v[523] - v[2076] * v[524] - v[2075] * v[525] - v[2074] * v[526];
	v[2305] = -(v[2088] * v[527]) - v[2087] * v[528] - v[2086] * v[529] - v[2085] * v[530] - v[2084] * v[531]
		- v[2083] * v[532] - v[2082] * v[533] - v[2081] * v[534] - v[2080] * v[535] - v[2079] * v[536] - v[2078] * v[537]
		- v[2077] * v[538] - v[2076] * v[539] - v[2075] * v[540] - v[2074] * v[541];
	v[3573] = (v[2088] * v[497] + v[2087] * v[498] + v[2086] * v[499] + v[2085] * v[500] + v[2084] * v[501] + v[2083] * v[502]
		+ v[2082] * v[503] + v[2081] * v[504] + v[2080] * v[505] + v[2079] * v[506] + v[2078] * v[507] + v[2077] * v[508]
		+ v[2076] * v[509] + v[2075] * v[510] + v[2074] * v[511]) / 2e0;
	v[3574] = (v[2088] * v[482] + v[2087] * v[483] + v[2086] * v[484] + v[2085] * v[485] + v[2084] * v[486] + v[2083] * v[487]
		+ v[2082] * v[488] + v[2081] * v[489] + v[2080] * v[490] + v[2079] * v[491] + v[2078] * v[492] + v[2077] * v[493]
		+ v[2076] * v[494] + v[2075] * v[495] + v[2074] * v[496]) / 2e0;
	v[2308] = -(v[2057] * v[512]) - v[2056] * v[513] - v[2055] * v[514] - v[2054] * v[515] - v[2053] * v[516]
		- v[2052] * v[517] - v[2051] * v[518] - v[2050] * v[519] - v[2049] * v[520] - v[2048] * v[521] - v[2047] * v[522]
		- v[2046] * v[523] - v[2045] * v[524] - v[2044] * v[525] - v[2043] * v[526];
	v[2309] = -(v[2057] * v[527]) - v[2056] * v[528] - v[2055] * v[529] - v[2054] * v[530] - v[2053] * v[531]
		- v[2052] * v[532] - v[2051] * v[533] - v[2050] * v[534] - v[2049] * v[535] - v[2048] * v[536] - v[2047] * v[537]
		- v[2046] * v[538] - v[2045] * v[539] - v[2044] * v[540] - v[2043] * v[541];
	v[3571] = (v[2057] * v[497] + v[2056] * v[498] + v[2055] * v[499] + v[2054] * v[500] + v[2053] * v[501] + v[2052] * v[502]
		+ v[2051] * v[503] + v[2050] * v[504] + v[2049] * v[505] + v[2048] * v[506] + v[2047] * v[507] + v[2046] * v[508]
		+ v[2045] * v[509] + v[2044] * v[510] + v[2043] * v[511]) / 2e0;
	v[3572] = (v[2057] * v[482] + v[2056] * v[483] + v[2055] * v[484] + v[2054] * v[485] + v[2053] * v[486] + v[2052] * v[487]
		+ v[2051] * v[488] + v[2050] * v[489] + v[2049] * v[490] + v[2048] * v[491] + v[2047] * v[492] + v[2046] * v[493]
		+ v[2045] * v[494] + v[2044] * v[495] + v[2043] * v[496]) / 2e0;
	v[2312] = v[1824] * v[1832] + v[1171] * v[1977] + v[1129] * v[1978] + v[1086] * v[1979] + v[1872] * v[3356]
		- v[1559] * v[3398] - v[1556] * v[3399] + v[1555] * v[3492] + v[1923] * v[3494] + (v[1842] + v[3495])*v[408] + (v[1816]
			+ v[3496])*v[411] + v[281] * (-(v[1964] * v[408]) - v[1953] * v[411]) + v[3493] * v[414] - v[1997] * v[496]
		- v[2001] * v[511] + v[2003] * v[526] + v[2002] * v[541] + v[1052] * (-(v[1636] * v[414]) - v[1555] * v[592]) + v[9] * (-
		(v[2010] * v[966]) - v[2011] * v[967] - v[2012] * v[968]);
	v[2388] = v[2312] * v[3337];
	v[2313] = v[1824] * v[1834] + v[1173] * v[1977] + v[1131] * v[1978] + v[1088] * v[1979] + v[1872] * v[3358]
		- v[1562] * v[3398] - v[1547] * v[3399] + v[1536] * v[3492] + v[1925] * v[3494] + (v[1842] + v[3495])*v[407] + (v[1816]
			+ v[3496])*v[410] + v[281] * (-(v[1964] * v[407]) - v[1953] * v[410]) + v[3493] * v[413] - v[1997] * v[495]
		- v[2001] * v[510] + v[2003] * v[525] + v[2002] * v[540] + v[1052] * (-(v[1636] * v[413]) - v[1536] * v[592]) + v[9] * (-
		(v[2010] * v[962]) - v[2011] * v[963] - v[2012] * v[964]);
	v[3549] = v[227] * v[2313];
	v[2392] = v[2313] * v[3336];
	v[2314] = v[1824] * v[1836] + v[1175] * v[1977] + v[1133] * v[1978] + v[1090] * v[1979] + v[1872] * v[3360]
		- v[1558] * v[3398] - v[1557] * v[3399] + v[1537] * v[3492] + v[1927] * v[3494] + (v[1842] + v[3495])*v[406] + (v[1816]
			+ v[3496])*v[409] + v[281] * (-(v[1964] * v[406]) - v[1953] * v[409]) + v[3493] * v[412] - v[1997] * v[494]
		- v[2001] * v[509] + v[2003] * v[524] + v[2002] * v[539] + v[1052] * (-(v[1636] * v[412]) - v[1537] * v[592]) + v[9] * (-
		(v[2010] * v[958]) - v[2011] * v[959] - v[2012] * v[960]);
	v[3550] = v[227] * v[2314];
	v[2395] = v[2314] * v[3335];
	v[2227] = v[2227] - v[2226] * v[406] - v[2225] * v[407] - v[2224] * v[408];
	v[1650] = v[1650] + v[1832] * v[2224] + v[1834] * v[2225] + v[1836] * v[2226] + v[2227] * v[3472];
	v[2180] = v[2180] - v[2177] * v[412] - v[2176] * v[413] - v[2175] * v[414];
	v[2315] = v[249] * (v[1905] + v[3533]) - v[279] * (v[1914] + v[3538] + v[3544]) - v[2224] * v[582];
	v[2316] = v[239] * (v[1905] + v[3533]) - v[279] * (v[1913] + v[3537] + v[3542]) - v[2225] * v[582];
	v[1650] = v[1650] - v[1911] * v[221] - v[1912] * v[406] - v[1913] * v[407] - v[1914] * v[408] + v[281] * (-
		(v[2177] * v[406]) - v[2176] * v[407] - v[2175] * v[408]) + v[280] * (v[3534] - v[2204] * v[406] - v[2203] * v[407]
			- v[2202] * v[408]);
	v[2205] = v[2205] - v[2204] * v[409] - v[2203] * v[410] - v[2202] * v[411];
	v[2317] = v[229] * (v[1905] + v[3533]) - v[279] * (v[1912] + v[3535] + v[3539]) - v[2226] * v[582];
	v[1651] = v[1651] + v[2202] * v[3356] + v[2203] * v[3358] + v[2204] * v[3360] + v[2205] * v[3469] + v[281] * (-
		(v[2177] * v[409]) - v[2176] * v[410] - v[2175] * v[411]) + v[279] * (v[3534] - v[2226] * v[409] - v[2225] * v[410]
			- v[2224] * v[411]);
	v[2318] = v[229] * v[3536] - v[280] * (v[1918] + v[3535] + v[3540]) - v[2204] * v[587];
	v[2319] = v[239] * v[3536] - v[280] * (v[1920] + v[3537] + v[3543]) - v[2203] * v[587];
	v[1651] = v[1651] - v[1916] * v[222] - v[1918] * v[409] - v[1920] * v[410] - v[1922] * v[411];
	v[2320] = v[249] * v[3536] - v[280] * (v[1922] + v[3538] + v[3545]) - v[2202] * v[587];
	v[2321] = v[1876] + v[3517] - v[3518];
	v[3579] = v[2321] * v[281] - v[3541] + v[3570];
	v[1639] = v[1639] + v[1923] * v[2175] + v[1925] * v[2176] + v[1927] * v[2177] - v[1876] * v[223] + v[1637] * v[2321]
		+ v[2180] * v[3466] + v[280] * (-(v[2204] * v[412]) - v[2203] * v[413] - v[2202] * v[414]) + v[279] * (-(v[2226] * v[412]
			) - v[2225] * v[413] - v[2224] * v[414]);
	v[2322] = v[281] * (-v[1877] - v[3539] - v[3540]) + v[229] * (v[3493] + v[3541]) - v[2177] * v[592];
	v[2323] = v[239] * (v[3493] + v[3541]) + v[281] * (-v[1878] - v[3542] - v[3543]) - v[2176] * v[592];
	v[1639] = v[1639] + v[2184] * v[2905] + v[2181] * v[2907] - v[1877] * v[412] - v[1878] * v[413] - v[1879] * v[414];
	v[2324] = v[249] * (v[3493] + v[3541]) + v[281] * (-v[1879] - v[3544] - v[3545]) - v[2175] * v[592];
	b2325 = b277;
	if (b2325) {
		v[2326] = -v[1639];
		v[2327] = -v[1651];
		v[2328] = -v[1650];
	}
	else {
		v[2326] = v[1639];
		v[2327] = v[1651];
		v[2328] = v[1650];
	};
	v[1625] = v[1625] + v[1301] * v[1627] * v[262] + (v[1296] * v[1619] + v[1295] * v[1621] + v[1294] * v[1623]
		+ v[2328] * v[250] + v[2327] * v[251] + v[2326] * v[252])*v[263];
	v[2337] = v[1294] * v[1626] + v[2326] * v[264] + (v[1053] * v[1623] + v[1625] * v[252]) / v[602];
	v[2339] = v[1295] * v[1626] + v[2327] * v[264] + (v[1053] * v[1621] + v[1625] * v[251]) / v[602];
	v[2341] = v[1296] * v[1626] + v[2328] * v[264] + (v[1053] * v[1619] + v[1625] * v[250]) / v[602];
	v[2273] = v[2273] - v[2337];
	v[2266] = v[2266] - v[2339];
	v[2259] = v[2259] - v[2341];
	v[2342] = v[1645] - v[1473] * v[1863] + v[152] * v[2337] + v[148] * v[2339] + v[144] * v[2341] + (v[212] * v[2181]
		+ v[213] * v[2184] + v[214] * v[2321])*v[281] + v[2228] * v[3361] + v[1471] * v[3441] + v[1157] * (v[3470] + v[3471])
		+ v[212] * v[3568] + v[213] * v[3569] + v[214] * v[3570] + v[3494] * v[3805] + v[1052] * (v[1636] * v[214]
			+ v[1475] * v[592]);
	v[2343] = v[2312] * v[236] + v[1172] * v[3546];
	v[2344] = v[230] * v[2312] + v[1172] * v[3547];
	v[2345] = v[1315] * v[2312] + v[1172] * (v[1576] / v[227] + v[1545] * v[3177]) + v[2343] * v[3690];
	v[2346] = v[1321] * v[2312] + v[1172] * (v[1568] / v[227] + v[1545] * v[3179]) + v[2344] * v[3691];
	v[2347] = (v[186] * v[2343] + v[184] * v[2344] + v[1172] * (v[1584] + v[3181] * v[3455]) + v[2312] * v[3692]) / v[227];
	v[2348] = v[2313] * v[238] + v[1174] * v[3548];
	v[2349] = v[230] * v[2313] + v[1174] * v[3547];
	v[2352] = v[1310] * v[2313] + v[1174] * (v[1582] / v[227] + v[1545] * v[3185]) + v[2348] * v[3694];
	v[2353] = v[2343] + v[2348];
	v[2354] = v[2345] + v[2352];
	v[2356] = (v[1313] * v[1530] + v[175] * v[2349] + v[177] * v[2353] + v[2482] * v[3455] + v[1174] * (v[1564]
		+ v[3189] * v[3455]) + v[1320] * v[3549]) / v[227];
	v[2359] = (v[182] * v[2348] + v[179] * v[2349] + v[1174] * (v[1573] + v[3191] * v[3455]) + v[1314] * v[3549]) / v[227];
	v[2360] = v[2314] * v[236] + v[1176] * v[3546];
	v[2361] = v[2314] * v[238] + v[1176] * v[3548];
	v[2362] = v[2349] + v[2360];
	v[2365] = (v[176] * v[2360] + v[177] * v[2361] + v[1176] * (v[1565] + v[3196] * v[3455]) + v[1318] * v[3550]) / v[227];
	v[2366] = v[2344] + v[2361];
	v[2368] = (v[1323] * v[1524] + v[181] * v[2360] + v[182] * v[2366] + v[2481] * v[3455] + v[1176] * (v[1572]
		+ v[3199] * v[3455]) + v[1325] * v[3550]) / v[227];
	v[2369] = v[2356] + v[2368];
	v[2371] = (v[1324] * v[1520] + v[187] * v[2361] + v[186] * v[2362] + v[2480] * v[3455] + v[1176] * (v[1580]
		+ v[3202] * v[3455]) + v[1328] * v[3550]) / v[227];
	v[2372] = v[2346] + v[2371];
	v[2274] = v[2274] - v[168] * v[2337] + v[2308] * v[295] + v[2309] * v[305];
	v[2275] = v[2275] - v[167] * v[2337] + v[2308] * v[294] + v[2309] * v[303];
	v[2276] = v[2276] - v[166] * v[2337] + v[2308] * v[293] + v[2309] * v[301];
	v[2267] = v[2267] - v[168] * v[2339] + v[2304] * v[295] + v[2305] * v[305];
	v[2268] = v[2268] - v[167] * v[2339] + v[2304] * v[294] + v[2305] * v[303];
	v[2269] = v[2269] - v[166] * v[2339] + v[2304] * v[293] + v[2305] * v[301];
	v[2260] = v[2260] - v[168] * v[2341] + v[2300] * v[295] + v[2301] * v[305];
	v[2261] = v[2261] - v[167] * v[2341] + v[2300] * v[294] + v[2301] * v[303];
	v[2262] = v[2262] - v[166] * v[2341] + v[2300] * v[293] + v[2301] * v[301];
	v[2415] = v[1437] * v[2348] + v[2388] * v[240] + v[2361] * v[3344] + v[238] * (v[1605] / v[227] + v[1545] * v[3806])
		+ v[2276] * v[94] + v[2275] * v[95] + v[2274] * v[96];
	v[2416] = v[237] * v[2392] + (v[1324] * v[1591] + v[1608] * v[236] + v[225] * v[2362] + v[2343] * v[240]) / v[227]
		+ v[1545] * v[2477] + v[2377] * v[3551] + v[2378] * v[3551] + v[2276] * v[91] + v[2275] * v[92] + v[2274] * v[93];
	v[2417] = v[225] * v[2395] + v[2344] * v[2512] + v[230] * (v[1610] / v[227] + v[1545] * v[3552]) + v[2276] * v[88]
		+ v[2275] * v[89] + v[2274] * v[90];
	v[2418] = (v[1323] * v[1595] + v[231] * v[2348] + v[226] * v[2366] + v[1616] * v[238]) / v[227] + v[1545] * v[2478]
		+ v[2388] * v[3345] + (v[2390] + v[2391])*v[3553] + v[2269] * v[94] + v[2268] * v[95] + v[2267] * v[96];
	v[2419] = v[1439] * v[2343] + v[231] * v[2392] + v[2360] * v[3341] + v[236] * (v[1618] / v[227] + v[1545] * v[3807])
		+ v[2269] * v[91] + v[2268] * v[92] + v[2267] * v[93];
	v[2420] = v[226] * v[2395] + v[2349] * v[2510] + v[230] * (v[1604] / v[227] + v[1545] * v[3554]) + v[2269] * v[88]
		+ v[2268] * v[89] + v[2267] * v[90];
	v[2421] = v[1545] * v[2479] + (v[1313] * v[1596] + v[228] * v[2361] + v[1612] * v[238] + v[2353] * v[3340]) / v[227]
		+ v[2388] * v[3343] + (v[2404] + v[2406])*v[3553] + v[2262] * v[94] + v[2261] * v[95] + v[2260] * v[96];
	v[2422] = v[2360] * v[2507] + v[2392] * v[3340] + v[236] * (v[1601] / v[227] + v[1545] * v[3555]) + v[2262] * v[91]
		+ v[2261] * v[92] + v[2260] * v[93];
	v[2423] = v[1440] * v[2344] + v[1438] * v[2349] + v[228] * v[2395] + v[230] * (v[1614] / v[227] + v[1545] * v[3808])
		+ v[2262] * v[88] + v[2261] * v[89] + v[2260] * v[90];
	v[2424] = -(v[2316] * v[742]);
	v[2425] = -(v[2316] * v[740]);
	v[2426] = -(v[2316] * v[738]);
	v[2427] = -(v[2315] * v[742]);
	v[2428] = -(v[2315] * v[738]);
	v[2429] = -(v[2315] * v[740]);
	v[2430] = -(v[2317] * v[740]);
	v[2431] = -(v[2317] * v[738]);
	v[2432] = -(v[2317] * v[742]);
	v[2433] = -(v[2318] * v[742]);
	v[2434] = -(v[2318] * v[740]);
	v[2435] = -(v[2318] * v[738]);
	v[2436] = -(v[2320] * v[738]);
	v[2437] = -(v[2320] * v[742]);
	v[2438] = -(v[2320] * v[740]);
	v[2439] = -(v[2324] * v[740]);
	v[2440] = -(v[2324] * v[742]);
	v[2441] = -(v[2324] * v[738]);
	v[2442] = -(v[1317] * v[1505]) - v[1327] * v[1507] + 2e0*v[1316] * v[1508] + v[2430] + v[2433] + v[2436] + v[2439]
		- v[2369] * v[326] - v[2354] * v[329] + v[2359] * v[3454];
	v[2443] = -(v[2319] * v[742]);
	v[2444] = -(v[2319] * v[738]);
	v[2445] = -(v[2319] * v[740]);
	v[2446] = v[1330] * v[1505] + 2e0*v[1322] * v[1507] - v[1327] * v[1508] - v[2425] - v[2428] - v[2440] - v[2443]
		- v[2369] * v[323] + v[2372] * v[329] + v[2365] * v[3453];
	v[2447] = -(v[1375] * v[1496]) + v[1371] * v[1499] - v[1376] * v[1501] - v[1368] * v[1506] + v[172] * v[2422]
		- v[2430] * v[316] + v[2425] * v[317] - v[2429] * v[318];
	v[2448] = -(v[2322] * v[742]);
	v[2449] = -(v[2322] * v[740]);
	v[2450] = -(v[2322] * v[738]);
	v[2451] = -(v[2323] * v[740]);
	v[2452] = -(v[2323] * v[742]);
	v[2453] = -(v[2323] * v[738]);
	v[2454] = 2e0*v[1309] * v[1505] + v[1330] * v[1507] - v[1317] * v[1508] - v[2431] - v[2444] - v[2448] - v[2451]
		- v[2354] * v[323] + v[2372] * v[326] + v[2347] * v[3452];
	v[2455] = -(v[1374] * v[1496]) + v[1372] * v[1499] - v[1377] * v[1501] - v[1367] * v[1506] + v[172] * v[2421]
		- v[2431] * v[316] + v[2426] * v[317] - v[2428] * v[318];
	v[2456] = -(v[1383] * v[1496]) + v[1389] * v[1499] - v[1379] * v[1501] - v[1357] * v[1506] + v[172] * v[2420]
		- v[2433] * v[316] + v[2443] * v[317] - v[2437] * v[318];
	v[2457] = v[2427] + v[2438];
	v[2458] = -(v[1382] * v[1496]) + v[1390] * v[1499] - v[1381] * v[1501] - v[1355] * v[1506] + v[172] * v[2418]
		- v[2435] * v[316] + v[2444] * v[317] - v[2436] * v[318];
	v[2459] = -(v[1386] * v[1496]) + v[1398] * v[1499] - v[1394] * v[1501] - v[1345] * v[1506] + v[172] * v[2417]
		- v[2448] * v[316] + v[2452] * v[317] - v[2440] * v[318];
	v[2460] = -(v[1385] * v[1496]) + v[1397] * v[1499] - v[1395] * v[1501] - v[1344] * v[1506] + v[172] * v[2416]
		- v[2449] * v[316] + v[2451] * v[317] - v[2439] * v[318];
	v[2461] = v[2434] + v[2450];
	v[2462] = v[2424] + v[2453];
	v[2463] = -(v[1305] * v[1571]) - v[1304] * v[1579] - v[1303] * v[1587] - v[199] * v[2337] - v[196] * v[2339]
		- v[193] * v[2341] + v[2317] * v[367] + v[2316] * v[368] + v[2315] * v[369] + v[2318] * v[382] + v[2319] * v[383]
		+ v[2320] * v[384] + v[2322] * v[397] + v[2323] * v[398] + v[2324] * v[399] - v[1563] * v[90] - v[1561] * v[93]
		- v[1560] * v[96];
	v[2464] = -(v[1305] * v[1570]) - v[1304] * v[1578] - v[1303] * v[1586] - v[198] * v[2337] - v[195] * v[2339]
		- v[192] * v[2341] + v[2317] * v[364] + v[2316] * v[365] + v[2315] * v[366] + v[2318] * v[379] + v[2319] * v[380]
		+ v[2320] * v[381] + v[2322] * v[394] + v[2323] * v[395] + v[2324] * v[396] - v[1563] * v[89] - v[1561] * v[92]
		- v[1560] * v[95];
	v[2465] = -(v[1305] * v[1569]) - v[1304] * v[1577] - v[1303] * v[1585] - v[197] * v[2337] - v[194] * v[2339]
		- v[191] * v[2341] + v[2317] * v[361] + v[2316] * v[362] + v[2315] * v[363] + v[2318] * v[376] + v[2319] * v[377]
		+ v[2320] * v[378] + v[2322] * v[391] + v[2323] * v[392] + v[2324] * v[393] - v[1563] * v[88] - v[1561] * v[91]
		- v[1560] * v[94];
	v[2466] = v[2465] * v[76] + v[2464] * v[77] + v[2463] * v[78];
	v[2467] = (-v[2466] + v[2465] * v[82] + v[2464] * v[83] + v[2463] * v[84]) / 2e0;
	v[2468] = (-v[2466] + v[2465] * v[79] + v[2464] * v[80] + v[2463] * v[81]) / 2e0;
	v[2469] = (-(v[1396] * v[1461]) - v[1378] * v[1464] - v[1380] * v[1497] - 2e0*v[2345] + 2e0*v[2352] - v[2432] * v[319]
		- v[2434] * v[338] + v[2431] * v[3429] + v[2430] * v[3430] + v[2435] * v[3431] + v[2433] * v[3432] + v[2449] * v[3433]
		+ v[2448] * v[3434] - v[1379] * v[3556] - v[1376] * v[3557] - v[1395] * v[3558] - v[1381] * v[3559] - v[2450] * v[356]
		- v[1394] * v[3560] - v[1377] * v[3561]) / 2e0;
	v[2470] = (-(v[1373] * v[1496]) + v[1370] * v[1499] - v[1378] * v[1501] - v[1369] * v[1506] + v[172] * v[2423]
		- v[2432] * v[316] + v[2424] * v[317] - v[2427] * v[318]) / 2e0;
	v[2472] = (v[1399] * v[1461] + v[1370] * v[1464] + v[1391] * v[1497] - 2e0*v[2346] + 2e0*v[2371] + v[2424] * v[319]
		+ v[2445] * v[338] - v[2426] * v[3429] - v[2425] * v[3430] - v[2444] * v[3431] - v[2443] * v[3432] - v[2451] * v[3433]
		- v[2452] * v[3434] + v[1389] * v[3556] + v[1371] * v[3557] + v[1397] * v[3558] + v[1390] * v[3559] + v[2453] * v[356]
		+ v[1398] * v[3560] + v[1372] * v[3561]) / 2e0;
	v[2473] = (-(v[1384] * v[1496]) + v[1391] * v[1499] - v[1380] * v[1501] - v[1356] * v[1506] + v[172] * v[2419]
		- v[2434] * v[316] + v[2445] * v[317] - v[2438] * v[318]) / 2e0;
	v[2474] = (-(v[1387] * v[1461]) - v[1373] * v[1464] - v[1384] * v[1497] - 2e0*v[2356] + 2e0*v[2368] - v[2427] * v[319]
		- v[2438] * v[338] + v[2428] * v[3429] + v[2429] * v[3430] + v[2436] * v[3431] + v[2437] * v[3432] + v[2439] * v[3433]
		+ v[2440] * v[3434] - v[1383] * v[3556] - v[1375] * v[3557] - v[1385] * v[3558] - v[1382] * v[3559] - v[2441] * v[356]
		- v[1386] * v[3560] - v[1374] * v[3561]) / 2e0;
	v[2509] = v[1498] * v[2472] - v[2473] + v[1495] * v[2474] + v[1463] * v[3710] * v[3810] - v[3330] * ((v[1343] * v[1461])
		/ 2e0 + (v[1369] * v[1464]) / 2e0 + v[1357] * v[1489] + v[1368] * v[1490] + v[1344] * v[1491] + v[1355] * v[1492]
		+ v[1345] * v[1493] + v[1367] * v[1494] + (v[1356] * v[1497]) / 2e0 + v[2426] - v[2429] - v[2435] + v[2437] + v[2449]
		- v[2452] - 2e0*v[1506] * v[2476] + v[2442] * v[311] - v[2446] * v[312] - v[2454] * v[314] + v[2422] * v[324]
		+ v[2421] * v[330] + v[2420] * v[334] + v[2415] * v[3426] + v[2419] * v[3427] + v[2423] * v[3428] + v[2418] * v[343]
		+ v[2417] * v[347] + v[2416] * v[351] + v[2457] * v[3562] + v[2461] * v[3563] + v[2462] * v[3564] + v[3443] * (v[2347]
			+ v[1608] * v[2350] + v[1605] * v[2351] + v[1599] * v[2353] + v[1614] * v[2355] + v[1618] * v[2357] + v[1616] * v[2358]
			+ v[2359] + v[1593] * v[2362] + v[1601] * v[2363] + v[1612] * v[2364] + v[2365] + v[1597] * v[2366] + v[1604] * v[2367]
			+ v[1610] * v[2370] + v[1584] * v[2377] + v[1582] * v[2378] + v[1580] * v[2380] + v[1576] * v[2390] + v[1573] * v[2391]
			+ v[1572] * v[2394] + v[1565] * v[2404] + v[1568] * v[2406] + v[1564] * v[2407] + v[1520] * v[2477] + v[1524] * v[2478]
			+ v[1530] * v[2479] + v[1591] * v[2480] + v[1595] * v[2481] + v[1596] * v[2482] + v[2343] * v[2486] + v[2344] * v[2487]
			+ v[2348] * v[2488] + v[2349] * v[2489] + v[2360] * v[2490] + v[2361] * v[2491] + v[2312] * v[3565] + v[2313] * v[3566]
			+ v[2314] * v[3567] + v[1545] * v[3711] * v[3812]) + v[5950 + i971]) + v[3708] * (-(v[1413] * v[3583])
				+ v[2469] * v[3709] + v[5965 + i971]);
	v[3592] = v[2509] + (v[1387] * v[1496] - v[1399] * v[1499] + v[1396] * v[1501] + v[1343] * v[1506] - v[172] * v[2415]
		+ v[2450] * v[316] - v[2453] * v[317] + v[2441] * v[318]) / 2e0;
	v[2493] = v[2455] + v[2459];
	v[2494] = v[2458] + v[2460];
	v[2495] = v[2447] + v[2456];
	v[2496] = (v[1643] - v[1485] * v[1863] + v[154] * v[2337] + v[150] * v[2339] + v[146] * v[2341] - v[2342] +
		(v[218] * v[2181] + v[2184] * v[219] + v[220] * v[2321])*v[281] + v[2230] * v[3361] + v[1483] * v[3441] + v[1159] *
		(v[3470] + v[3471]) + v[218] * v[3568] + v[219] * v[3569] + v[220] * v[3570] + v[3494] * v[3813] + v[1052] *
		(v[1636] * v[220] + v[1487] * v[592])) / 2e0;
	v[2497] = (v[1644] - v[1479] * v[1863] + v[153] * v[2337] + v[149] * v[2339] + v[145] * v[2341] - v[2342] +
		(v[215] * v[2181] + v[216] * v[2184] + v[217] * v[2321])*v[281] + v[2229] * v[3361] + v[1477] * v[3441] + v[1158] *
		(v[3470] + v[3471]) + v[215] * v[3568] + v[216] * v[3569] + v[217] * v[3570] + v[3494] * v[3814] + v[1052] *
		(v[1636] * v[217] + v[1481] * v[592])) / 2e0;
	v[2277] = v[2277] + v[132] * v[2337] + v[3571];
	v[2270] = v[2270] + v[132] * v[2339] + v[3573];
	v[2263] = v[2263] + v[132] * v[2341] + v[3575];
	v[2278] = v[2278] + v[130] * v[2337] + v[3572];
	v[2271] = v[2271] + v[130] * v[2339] + v[3574];
	v[2264] = v[2264] + v[130] * v[2341] + v[3576];
	v[2279] = v[2279] + v[128] * v[2337] - v[3571] - v[3572];
	v[2272] = v[2272] + v[128] * v[2339] - v[3573] - v[3574];
	v[2265] = v[2265] + v[128] * v[2341] - v[3575] - v[3576];
	v[6071] = v[2265] + v[2497] * v[482] + v[2496] * v[497] + v[2468] * v[512] + v[2467] * v[527] + v[10] * (v[1051] * v[1657]
		+ v[128] * v[3577] - v[1997] * v[482] - v[2001] * v[497] + v[2003] * v[512] + v[2002] * v[527] + v[1978] * v[584] - v[9] *
		(v[2010] * v[910] + v[2011] * v[911] + v[2012] * v[912]));
	v[6072] = v[2272] + v[2497] * v[483] + v[2496] * v[498] + v[2468] * v[513] + v[2467] * v[528] + v[10] * (v[1050] * v[1657]
		+ v[128] * v[3578] - v[1997] * v[483] - v[2001] * v[498] + v[2003] * v[513] + v[2002] * v[528] + v[1977] * v[584] - v[9] *
		(v[2010] * v[914] + v[2011] * v[915] + v[2012] * v[916]));
	v[6073] = v[2279] + v[2497] * v[484] + v[2496] * v[499] + v[2468] * v[514] + v[2467] * v[529] + v[10] * (v[128] * (v[3468]
		+ v[3579]) - v[1997] * v[484] - v[2001] * v[499] + v[2003] * v[514] + v[2002] * v[529] - v[9] * (v[2010] * v[918]
			+ v[2011] * v[919] + v[2012] * v[920]));
	v[6074] = v[2264] + v[2497] * v[485] + v[2496] * v[500] + v[2468] * v[515] + v[2467] * v[530] + v[10] * (v[1051] * v[1656]
		+ v[130] * v[3577] - v[1997] * v[485] - v[2001] * v[500] + v[2003] * v[515] + v[2002] * v[530] + v[1978] * v[585] - v[9] *
		(v[2010] * v[922] + v[2011] * v[923] + v[2012] * v[924]));
	v[6075] = v[2271] + v[2497] * v[486] + v[2496] * v[501] + v[2468] * v[516] + v[2467] * v[531] + v[10] * (v[1050] * v[1656]
		+ v[130] * v[3578] - v[1997] * v[486] - v[2001] * v[501] + v[2003] * v[516] + v[2002] * v[531] + v[1977] * v[585] - v[9] *
		(v[2010] * v[926] + v[2011] * v[927] + v[2012] * v[928]));
	v[6076] = v[2278] + v[2497] * v[487] + v[2496] * v[502] + v[2468] * v[517] + v[2467] * v[532] + v[10] * (v[130] * (v[3468]
		+ v[3579]) - v[1997] * v[487] - v[2001] * v[502] + v[2003] * v[517] + v[2002] * v[532] - v[9] * (v[2010] * v[930]
			+ v[2011] * v[931] + v[2012] * v[932]));
	v[6077] = v[2263] + v[2497] * v[488] + v[2496] * v[503] + v[2468] * v[518] + v[2467] * v[533] + v[10] * (v[1051] * v[1655]
		+ v[132] * v[3577] - v[1997] * v[488] - v[2001] * v[503] + v[2003] * v[518] + v[2002] * v[533] + v[1978] * v[586] - v[9] *
		(v[2010] * v[934] + v[2011] * v[935] + v[2012] * v[936]));
	v[6078] = v[2270] + v[2497] * v[489] + v[2496] * v[504] + v[2468] * v[519] + v[2467] * v[534] + v[10] * (v[1050] * v[1655]
		+ v[132] * v[3578] - v[1997] * v[489] - v[2001] * v[504] + v[2003] * v[519] + v[2002] * v[534] + v[1977] * v[586] - v[9] *
		(v[2010] * v[938] + v[2011] * v[939] + v[2012] * v[940]));
	v[6079] = v[2277] + v[2497] * v[490] + v[2496] * v[505] + v[2468] * v[520] + v[2467] * v[535] + v[10] * (v[132] * (v[3468]
		+ v[3579]) - v[1997] * v[490] - v[2001] * v[505] + v[2003] * v[520] + v[2002] * v[535] - v[9] * (v[2010] * v[942]
			+ v[2011] * v[943] + v[2012] * v[944]));
	v[6080] = v[2259] + v[2497] * v[491] + v[2496] * v[506] + v[2468] * v[521] + v[2467] * v[536] + v[10] * (v[1905] + v[3533]
		- v[3580] - v[279] * (v[1911] + v[3517] + v[3581]) - v[1997] * v[491] - v[2001] * v[506] + v[2003] * v[521]
		+ v[2002] * v[536] - v[9] * (v[2010] * v[946] + v[2011] * v[947] + v[2012] * v[948]));
	v[6081] = v[2266] + v[2497] * v[492] + v[2496] * v[507] + v[2468] * v[522] + v[2467] * v[537] - v[10] * (-v[1816]
		+ v[3569] + v[280] * (v[1916] - v[3518] + v[3581]) - v[3582] + v[1997] * v[492] + v[2001] * v[507] - v[2003] * v[522]
		- v[2002] * v[537] + v[9] * (v[2010] * v[950] + v[2011] * v[951] + v[2012] * v[952]));
	v[6082] = v[2273] + v[2497] * v[493] + v[2496] * v[508] + v[2468] * v[523] + v[2467] * v[538] - v[10] * (-v[3493]
		+ v[3579] + v[1997] * v[493] + v[2001] * v[508] - v[2003] * v[523] - v[2002] * v[538] + v[9] * (v[2010] * v[954]
			+ v[2011] * v[955] + v[2012] * v[956]));
	v[6083] = -v[2458] + v[2460] + 2e0*i3448*v[2508] + v[2495] * v[311] + v[2493] * v[314] + (v[2446] + 2e0*v[2461]
		)*v[3332] + v[2469] * v[3450] + v[1506] * v[3587] + (v[1614] * v[175] + v[1604] * v[179] + v[1610] * v[184])*v[3588]
		+ v[10] * (v[1545] * v[2404] + v[1172] * (v[1598] + v[1545] * v[2487]) + v[1174] * (v[1594] + v[1545] * v[2489])
			+ v[2314] * v[2507] + v[2313] * v[3341] + v[2312] * v[3344] + v[1513] * v[3589] + v[1522] * v[3590] + v[1528] * v[3591])
		+ v[310] * v[3592] + v[3583] * v[3816] + v[2497] * v[494] + v[2496] * v[509] + v[2468] * v[524] + v[2467] * v[539] + v[6055
		+ i971] / 2e0;
	v[6084] = v[2455] - v[2459] + 2e0*i3449*v[2511] + v[2495] * v[312] + v[2494] * v[314] - (v[2442] - 2e0*v[2462]
		)*v[3332] - v[2472] * v[3450] - v[1506] * v[3586] + (v[1601] * v[176] + v[1618] * v[181] + v[1608] * v[186])*v[3588]
		+ v[313] * (-v[2470] + v[2473] + v[3592]) + v[10] * (v[1437] * v[2312] + v[1438] * v[2314] + v[1545] * v[2391] + v[1172] *
		(v[1600] + v[1545] * v[2486]) + v[1176] * (v[1594] + v[1545] * v[2490]) + v[2313] * v[2510] + v[1520] * v[3408]
			+ v[1539] * v[3410] + v[1511] * v[3593]) + v[3583] * v[3818] + v[2497] * v[495] + v[2496] * v[510] + v[2468] * v[525]
		+ v[2467] * v[540] + v[6010 + i971] / 2e0;
	v[6085] = -v[2447] + v[2456] + 2e0*i3447*v[2513] + v[2494] * v[311] + v[2493] * v[312] + (-v[2470] + v[2509])*v[315]
		+ (v[2454] + 2e0*v[2457])*v[3332] + v[2474] * v[3450] + v[1506] * v[3585] + (v[1612] * v[177] + v[1616] * v[182]
			+ v[1605] * v[187])*v[3588] + v[10] * (v[1439] * v[2313] + v[1440] * v[2314] + v[1545] * v[2377] + v[1174] * (v[1600]
				+ v[1545] * v[2488]) + v[1176] * (v[1598] + v[1545] * v[2491]) + v[2312] * v[2512] + v[1524] * v[3409]
				+ v[1509] * v[3594] + v[1530] * v[3595]) + v[3583] * v[3820] + v[2497] * v[496] + v[2496] * v[511] + v[2468] * v[526]
		+ v[2467] * v[541] + v[5980 + i971] / 2e0;
	Rc[i971 - 1] += v[5391 + i971];
	for (i1444 = 1; i1444 <= 15; i1444++) {
		Kc[i971 - 1][i1444 - 1] += v[6070 + i1444];
	};/* end for */
};/* end for */
v[2523] = 0e0;
v[2524] = 0e0;
v[2525] = 0e0;
v[2526] = 0e0;
b2527 = b603;
if (b2527) {
	b2528 = b622;
	b2532 = b605;
	if (b2532) {
		v[2525] = 0e0;
		v[2524] = 0e0;
		v[2523] = 0e0;
		v[2526] = 0e0;
	}
	else {
	};
}
else {
};
v[3601] = v[2523] * v[582];
v[3597] = v[2523] * v[279];
v[3600] = v[2524] * v[587];
v[3596] = v[2524] * v[280];
v[2533] = 0e0;
v[2534] = 0e0;
v[2535] = 0e0;
v[2536] = 0e0;
b2537 = b603;
if (b2537) {
	b2538 = b605;
	if (b2538) {
		v[2535] = 0e0;
		v[2534] = 0e0;
		v[2533] = 0e0;
		v[2536] = -v[2526];
	}
	else {
	};
}
else {
};
v[2539] = v[249] * v[2525];
v[2999] = -(v[2539] * v[281]);
v[2540] = v[239] * v[2525];
v[2996] = -(v[2540] * v[281]);
v[2541] = v[229] * v[2525];
v[3638] = -(v[2541] * v[406]) - v[2540] * v[407] - v[2539] * v[408];
v[3627] = -(v[2541] * v[409]) - v[2540] * v[410] - v[2539] * v[411];
v[2983] = -(v[2541] * v[281]);
v[2533] = v[2178] * v[2525] + v[2533];
v[2534] = v[2179] * v[2525] + v[2534];
v[2921] = v[2525] * v[2542] - v[2541] * v[412] - v[2540] * v[413] - v[2539] * v[414];
v[2544] = v[2525] * v[279];
v[2547] = v[2525] * v[280];
v[2535] = v[2535] + v[2525] * v[2550];
v[2554] = -(v[2525] * v[592]);
v[2555] = -(v[2547] * v[281]);
v[3609] = -v[2555] + v[3600];
v[2556] = -(v[2544] * v[281]);
v[3610] = -v[2556] + v[3601];
v[2566] = v[249] * v[2524];
v[3000] = -(v[2566] * v[280]);
v[2567] = v[239] * v[2524];
v[2997] = -(v[2567] * v[280]);
v[2568] = v[229] * v[2524];
v[3619] = -(v[2568] * v[412]) - v[2567] * v[413] - v[2566] * v[414];
v[2984] = -(v[2568] * v[280]);
v[3635] = -v[2983] - v[2984];
v[2533] = v[2533] - v[221] * v[3596];
v[2954] = v[2524] * v[2569] - v[2568] * v[409] - v[2567] * v[410] - v[2566] * v[411];
v[2534] = v[2534] + v[2524] * v[2571];
v[2535] = v[2535] - v[223] * v[3596];
v[2590] = v[249] * v[2523];
v[2963] = -(v[2590] * v[279]);
v[3848] = 2e0*v[2963] + v[2999] + v[3000];
v[3842] = v[2963] + v[2999] - v[2566] * v[3469];
v[3826] = v[2963] + v[3000] - v[2539] * v[3466];
v[3646] = -v[2963] - v[3000];
v[3643] = -v[2963] - v[2999];
v[2591] = v[239] * v[2523];
v[2965] = -(v[2591] * v[279]);
v[3847] = 2e0*v[2965] + v[2996] + v[2997];
v[3843] = v[2965] + v[2996] - v[2567] * v[3469];
v[3828] = v[2965] + v[2997] - v[2540] * v[3466];
v[3645] = -v[2965] - v[2997];
v[3642] = -v[2965] - v[2996];
v[2592] = v[229] * v[2523];
v[3617] = -(v[2592] * v[412]) - v[2591] * v[413] - v[2590] * v[414];
v[2961] = -(v[2592] * v[279]);
v[3849] = 2e0*v[2961] + v[2983] + v[2984];
v[3841] = v[2961] + v[2983] - v[2568] * v[3469];
v[3827] = v[2961] + v[2984] - v[2541] * v[3466];
v[3639] = -v[2961] - v[2984];
v[3636] = -v[2961] - v[2983];
v[2985] = v[2523] * v[2593] - v[2592] * v[406] - v[2591] * v[407] - v[2590] * v[408];
v[2595] = v[213] * v[2523] + v[212] * v[2524];
v[2596] = v[216] * v[2523] + v[215] * v[2524];
v[2597] = v[219] * v[2523] + v[218] * v[2524];
v[2989] = v[128] * v[2595] + v[130] * v[2596] + v[132] * v[2597];
v[3637] = v[2989] - v[2568] * v[406] - v[2567] * v[407] - v[2566] * v[408];
v[3626] = v[2989] - v[2592] * v[409] - v[2591] * v[410] - v[2590] * v[411];
v[2533] = v[2533] + v[2523] * v[2598];
v[2534] = v[2534] - v[222] * v[3597];
v[2599] = v[3596] + v[3597];
v[3837] = v[212] * v[2544] + v[213] * v[2547] + v[214] * v[2599];
v[3834] = v[215] * v[2544] + v[216] * v[2547] + v[217] * v[2599];
v[3831] = v[218] * v[2544] + v[219] * v[2547] + v[220] * v[2599];
v[6229] = 0e0;
v[6230] = 0e0;
v[6231] = v[128] * v[2599];
v[6232] = 0e0;
v[6233] = 0e0;
v[6234] = v[130] * v[2599];
v[6235] = 0e0;
v[6236] = 0e0;
v[6237] = v[132] * v[2599];
v[6238] = 0e0;
v[6239] = 0e0;
v[6240] = 0e0;
v[6241] = 0e0;
v[6242] = 0e0;
v[6243] = 0e0;
v[2535] = v[2535] - v[223] * v[3597];
v[2602] = v[1171] * v[2523] + v[1129] * v[2524] + v[1086] * v[2525];
v[3598] = -(v[1412] * v[2602]);
v[3223] = v[3343] * v[3598];
v[3213] = v[3345] * v[3598];
v[3206] = v[240] * v[3598];
v[2666] = v[1439] * v[2602];
v[2662] = v[1440] * v[2602];
v[2652] = v[2512] * v[2602];
v[2603] = v[1173] * v[2523] + v[1131] * v[2524] + v[1088] * v[2525];
v[3224] = v[2603] * v[3417];
v[3214] = v[2484] * v[2603];
v[3207] = v[2603] * v[3599];
v[2670] = v[1437] * v[2603];
v[3721] = v[2652] + v[2670];
v[2656] = v[2510] * v[2603];
v[3724] = v[2656] + v[2666];
v[2604] = v[1175] * v[2523] + v[1133] * v[2524] + v[1090] * v[2525];
v[3221] = v[2485] * v[2604];
v[3701] = v[3221] + v[3224];
v[3865] = v[3223] + v[3701];
v[3217] = v[2604] * v[3418];
v[3700] = v[3214] + v[3217];
v[3864] = v[3213] + v[3700];
v[3209] = v[2604] * v[3420];
v[3698] = v[3206] + v[3209];
v[3863] = v[3207] + v[3698];
v[2660] = v[2507] * v[2604];
v[3725] = v[2660] + v[2662];
v[3606] = v[1438] * v[2603] + v[2660];
v[3717] = v[2662] + v[3606];
v[3605] = v[2656] + v[2604] * v[3341];
v[3720] = v[2666] + v[3605];
v[3604] = v[2652] + v[2604] * v[3344];
v[3723] = v[2670] + v[3604];
v[2605] = v[2554] - v[2599] * v[281];
v[6409] = -v[2556];
v[6410] = -v[2555];
v[6411] = -v[2605];
v[6412] = 0e0;
v[6413] = 0e0;
v[6414] = 0e0;
v[6415] = 0e0;
v[6416] = 0e0;
v[6417] = 0e0;
v[6418] = 0e0;
v[6419] = 0e0;
v[6420] = 0e0;
v[6421] = 0e0;
v[6422] = 0e0;
v[6423] = 0e0;
v[6379] = 0e0;
v[6380] = 0e0;
v[6381] = 0e0;
v[6382] = -v[2556];
v[6383] = -v[2555];
v[6384] = -v[2605];
v[6385] = 0e0;
v[6386] = 0e0;
v[6387] = 0e0;
v[6388] = 0e0;
v[6389] = 0e0;
v[6390] = 0e0;
v[6391] = 0e0;
v[6392] = 0e0;
v[6393] = 0e0;
v[6349] = 0e0;
v[6350] = 0e0;
v[6351] = 0e0;
v[6352] = 0e0;
v[6353] = 0e0;
v[6354] = 0e0;
v[6355] = -v[2556];
v[6356] = -v[2555];
v[6357] = -v[2605];
v[6358] = 0e0;
v[6359] = 0e0;
v[6360] = 0e0;
v[6361] = 0e0;
v[6362] = 0e0;
v[6363] = 0e0;
v[6105] = v[128] * (-v[2556] + v[3601]) + v[2524] * v[584];
v[6106] = v[128] * (-v[2555] + v[3600]) + v[2523] * v[584];
v[6107] = -(v[128] * v[2554]) + v[2599] * v[3714];
v[6108] = v[130] * (-v[2556] + v[3601]) + v[2524] * v[585];
v[6109] = v[130] * (-v[2555] + v[3600]) + v[2523] * v[585];
v[6110] = -(v[130] * v[2554]) + v[2599] * v[3713];
v[6111] = v[132] * (-v[2556] + v[3601]) + v[2524] * v[586];
v[6112] = v[132] * (-v[2555] + v[3600]) + v[2523] * v[586];
v[6113] = -(v[132] * v[2554]) + v[2599] * v[3712];
v[6114] = v[2556] - v[279] * v[3596] - v[3601];
v[6115] = v[2555] - v[280] * v[3597] - v[3600];
v[6116] = v[2605];
v[6117] = v[2660] + v[2603] * v[3341] + v[2602] * v[3344];
v[6118] = v[1437] * v[2602] + v[1438] * v[2604] + v[2656];
v[6119] = v[1439] * v[2603] + v[1440] * v[2604] + v[2652];
v[2617] = -(v[279] * (-v[2999] - v[3000])) - v[2590] * v[582];
v[2618] = -(v[280] * v[3643]) - v[2566] * v[587];
v[2619] = -(v[281] * v[3646]) - v[2539] * v[592];
v[2620] = -(v[279] * (-v[2996] - v[2997])) - v[2591] * v[582];
v[2621] = -(v[280] * v[3642]) - v[2567] * v[587];
v[2622] = -(v[281] * v[3645]) - v[2540] * v[592];
v[2533] = v[2533] + v[1832] * v[2590] + v[1834] * v[2591] + v[1836] * v[2592] + v[2985] * v[3472] + v[280] * v[3637]
+ v[281] * v[3638];
v[2623] = -(v[279] * v[3635]) - v[2592] * v[582];
v[2624] = -(v[280] * v[3636]) - v[2568] * v[587];
v[2625] = -(v[281] * v[3639]) - v[2541] * v[592];
v[2534] = v[2534] + v[2566] * v[3356] + v[2567] * v[3358] + v[2568] * v[3360] + v[2954] * v[3469] + v[279] * v[3626]
+ v[281] * v[3627];
v[2535] = v[2535] + v[1923] * v[2539] + v[1925] * v[2540] + v[1927] * v[2541] + v[1637] * v[2599] + v[2547] * v[2905]
+ v[2544] * v[2907] + v[2921] * v[3466] + v[279] * v[3617] + v[280] * v[3619];
b2626 = b277;
if (b2626) {
	v[2627] = -v[2535];
	v[2628] = -v[2534];
	v[2629] = -v[2533];
}
else {
	v[2627] = v[2535];
	v[2628] = v[2534];
	v[2629] = v[2533];
};
v[2634] = v[252] * v[2627] + v[251] * v[2628] + v[250] * v[2629];
v[2536] = v[2536] + v[263] * v[2634];
v[3602] = v[2536] / v[602];
v[2636] = v[2627] * v[264] + v[252] * v[3602] + v[626];
v[2637] = v[2628] * v[264] + v[251] * v[3602] + v[625];
v[2638] = v[2629] * v[264] + v[250] * v[3602] + v[624];
v[6424] = v[2638];
v[6425] = v[2637];
v[6426] = v[2636];
v[6427] = 0e0;
v[6428] = 0e0;
v[6429] = 0e0;
v[6430] = 0e0;
v[6431] = 0e0;
v[6432] = 0e0;
v[6433] = 0e0;
v[6434] = 0e0;
v[6435] = 0e0;
v[6436] = 0e0;
v[6437] = 0e0;
v[6438] = 0e0;
v[6394] = 0e0;
v[6395] = 0e0;
v[6396] = 0e0;
v[6397] = v[2638];
v[6398] = v[2637];
v[6399] = v[2636];
v[6400] = 0e0;
v[6401] = 0e0;
v[6402] = 0e0;
v[6403] = 0e0;
v[6404] = 0e0;
v[6405] = 0e0;
v[6406] = 0e0;
v[6407] = 0e0;
v[6408] = 0e0;
v[6364] = 0e0;
v[6365] = 0e0;
v[6366] = 0e0;
v[6367] = 0e0;
v[6368] = 0e0;
v[6369] = 0e0;
v[6370] = v[2638];
v[6371] = v[2637];
v[6372] = v[2636];
v[6373] = 0e0;
v[6374] = 0e0;
v[6375] = 0e0;
v[6376] = 0e0;
v[6377] = 0e0;
v[6378] = 0e0;
v[2639] = -(v[214] * v[2605]) + v[152] * v[2636] + v[148] * v[2637] + v[144] * v[2638] + v[2595] * v[3361]
+ v[213] * v[3609] + v[212] * v[3610];
v[2640] = v[1800] * v[2602];
v[2641] = v[1799] * v[2602];
v[2642] = v[1798] * v[2602];
v[2643] = v[1795] * v[2603];
v[2644] = v[236] * v[2602] + v[238] * v[2603];
v[3603] = v[177] * v[2644];
v[3316] = -(v[1412] * v[3603]);
v[3313] = v[2644] * v[3417];
v[2645] = v[1793] * v[2603];
v[2646] = v[2640] + v[2643];
v[2647] = v[1794] * v[2603] + v[3603] / v[227];
v[2648] = v[1788] * v[2604];
v[2649] = v[230] * v[2602] + v[238] * v[2604];
v[3607] = v[182] * v[2649];
v[3315] = -(v[1412] * v[3607]);
v[3312] = v[2649] * v[3418];
v[2650] = v[230] * v[2603] + v[236] * v[2604];
v[3608] = v[186] * v[2650];
v[3869] = -(v[2602] * v[3308]) - v[2603] * v[3309] - v[2604] * v[3310] + v[3340] * v[3603] + v[226] * v[3607]
+ v[225] * v[3608];
v[3314] = -(v[1412] * v[3608]);
v[3311] = v[2650] * v[3420];
v[3306] = v[1797] * v[2602] + v[1792] * v[2603] + v[1787] * v[2604] + v[2642] + v[1599] * v[2644] + v[2645] + v[2648]
+ v[1597] * v[2649] + v[1593] * v[2650];
v[2651] = v[238] * v[3723] + v[2636] * v[738];
v[2654] = v[230] * v[3604] + v[2636] * v[742];
v[2655] = v[236] * v[3720] + v[2637] * v[740];
v[2658] = v[230] * v[3605] + v[2637] * v[742];
v[2659] = v[1438] * v[2644] + v[238] * v[3725] + v[2638] * v[738];
v[2661] = v[236] * v[3606] + v[2638] * v[740];
v[2664] = v[230] * v[3717] + v[2638] * v[742];
v[2665] = v[1790] * v[2604] + v[3607] / v[227];
v[2667] = v[2649] * v[3341] + v[238] * v[3724] + v[2637] * v[738];
v[2668] = v[2647] + v[2665];
v[2669] = v[1789] * v[2604] + v[3608] / v[227];
v[2671] = v[2650] * v[3344] + v[236] * v[3721] + v[2636] * v[740];
v[2672] = v[2641] + v[2669];
v[2673] = -(v[2620] * v[742]);
v[2674] = -(v[2620] * v[740]);
v[2675] = -(v[2620] * v[738]);
v[2676] = -(v[2617] * v[742]);
v[2677] = -(v[2617] * v[738]);
v[2678] = -(v[2617] * v[740]);
v[2679] = -(v[2623] * v[740]);
v[2680] = -(v[2623] * v[738]);
v[2681] = -(v[2623] * v[742]);
v[2682] = -(v[2624] * v[742]);
v[2683] = -(v[2624] * v[740]);
v[2684] = -(v[2624] * v[738]);
v[2685] = -(v[2618] * v[738]);
v[2686] = -(v[2618] * v[742]);
v[2687] = -(v[2618] * v[740]);
v[2688] = -(v[2619] * v[740]);
v[2689] = -(v[2619] * v[742]);
v[2690] = -(v[2619] * v[738]);
v[2691] = v[2679] + v[2682] + v[2685] + v[2688] - v[2668] * v[326] - v[2646] * v[329] + v[2645] * v[3454];
v[2692] = -(v[2621] * v[742]);
v[2693] = -(v[2621] * v[738]);
v[2694] = -(v[2621] * v[740]);
v[2695] = -v[2674] - v[2677] - v[2689] - v[2692] - v[2668] * v[323] + v[2672] * v[329] + v[2648] * v[3453];
v[2696] = v[172] * v[2661] - v[2679] * v[316] + v[2674] * v[317] - v[2678] * v[318];
v[2697] = -(v[2625] * v[742]);
v[2698] = -(v[2625] * v[740]);
v[2699] = -(v[2625] * v[738]);
v[2700] = -(v[2622] * v[740]);
v[2701] = -(v[2622] * v[742]);
v[2702] = -(v[2622] * v[738]);
v[2703] = -(v[199] * v[2636]) - v[196] * v[2637] - v[193] * v[2638] + v[2623] * v[367] + v[2620] * v[368]
+ v[2617] * v[369] + v[2624] * v[382] + v[2621] * v[383] + v[2618] * v[384] + v[2625] * v[397] + v[2622] * v[398]
+ v[2619] * v[399];
v[2704] = -(v[198] * v[2636]) - v[195] * v[2637] - v[192] * v[2638] + v[2623] * v[364] + v[2620] * v[365]
+ v[2617] * v[366] + v[2624] * v[379] + v[2621] * v[380] + v[2618] * v[381] + v[2625] * v[394] + v[2622] * v[395]
+ v[2619] * v[396];
v[2705] = -(v[197] * v[2636]) - v[194] * v[2637] - v[191] * v[2638] + v[2623] * v[361] + v[2620] * v[362]
+ v[2617] * v[363] + v[2624] * v[376] + v[2621] * v[377] + v[2618] * v[378] + v[2625] * v[391] + v[2622] * v[392]
+ v[2619] * v[393];
v[2706] = -v[2680] - v[2693] - v[2697] - v[2700] - v[2646] * v[323] + v[2672] * v[326] + v[2642] * v[3452];
v[2707] = v[172] * v[2659] - v[2680] * v[316] + v[2675] * v[317] - v[2677] * v[318];
v[2708] = v[172] * v[2658] - v[2682] * v[316] + v[2692] * v[317] - v[2686] * v[318];
v[2709] = v[2676] + v[2687];
v[2710] = v[172] * v[2667] - v[2684] * v[316] + v[2693] * v[317] - v[2685] * v[318];
v[2711] = v[172] * v[2654] - v[2697] * v[316] + v[2701] * v[317] - v[2689] * v[318];
v[2712] = v[172] * v[2671] - v[2698] * v[316] + v[2700] * v[317] - v[2688] * v[318];
v[2713] = v[2683] + v[2699];
v[2714] = v[2673] + v[2702];
v[6334] = 0e0;
v[6335] = 0e0;
v[6336] = 0e0;
v[6337] = 0e0;
v[6338] = 0e0;
v[6339] = 0e0;
v[6340] = 0e0;
v[6341] = 0e0;
v[6342] = 0e0;
v[6343] = 0e0;
v[6344] = 0e0;
v[6345] = 0e0;
v[6346] = -v[2695] / 2e0 - v[2713];
v[6347] = v[2691] / 2e0 - v[2714];
v[6348] = -v[2706] / 2e0 - v[2709];
v[2715] = v[2675] - v[2678] - v[2684] + v[2686] + v[2698] - v[2701] + v[2691] * v[311] - v[2695] * v[312]
- v[2706] * v[314] + v[2661] * v[324] + v[2659] * v[330] + v[2658] * v[334] + v[2651] * v[3426] + v[2655] * v[3427]
+ v[2664] * v[3428] + v[2667] * v[343] + v[3306] * v[3443] + v[2654] * v[347] + v[2671] * v[351] + v[2709] * v[3562]
+ v[2713] * v[3563] + v[2714] * v[3564];
v[2716] = v[2705] * v[76] + v[2704] * v[77] + v[2703] * v[78];
v[2717] = (-v[2716] + v[2705] * v[82] + v[2704] * v[83] + v[2703] * v[84]) / 2e0;
v[2718] = (-v[2716] + v[2705] * v[79] + v[2704] * v[80] + v[2703] * v[81]) / 2e0;
v[2719] = (-2e0*v[2640] + 2e0*v[2643] - v[2681] * v[319] - v[2683] * v[338] + v[2680] * v[3429] + v[2679] * v[3430]
	+ v[2684] * v[3431] + v[2682] * v[3432] + v[2698] * v[3433] + v[2697] * v[3434] - v[2699] * v[356]) / 2e0;
v[2721] = -v[2641] + v[2669] + v[2674] * v[324] + v[2675] * v[330] + v[2692] * v[334] + v[2702] * v[3426]
+ v[2694] * v[3427] + v[2673] * v[3428] + v[2693] * v[343] + v[2701] * v[347] + v[2700] * v[351];
v[2722] = (v[172] * v[2655] - v[2683] * v[316] + v[2694] * v[317] - v[2687] * v[318]) / 2e0;
v[2723] = (-2e0*v[2647] + 2e0*v[2665] - v[2676] * v[319] - v[2687] * v[338] + v[2677] * v[3429] + v[2678] * v[3430]
	+ v[2685] * v[3431] + v[2686] * v[3432] + v[2688] * v[3433] + v[2689] * v[3434] - v[2690] * v[356]) / 2e0;
v[3867] = v[2719] * v[310] - v[2721] * v[313] + v[2723] * v[315];
v[6439] = 0e0;
v[6440] = 0e0;
v[6441] = 0e0;
v[6442] = 0e0;
v[6443] = 0e0;
v[6444] = 0e0;
v[6445] = 0e0;
v[6446] = 0e0;
v[6447] = 0e0;
v[6448] = 0e0;
v[6449] = 0e0;
v[6450] = 0e0;
v[6451] = 8e0*v[2719];
v[6452] = -8e0*v[2721];
v[6453] = 8e0*v[2723];
v[2730] = v[1500] * v[2719] + v[1498] * v[2721] - v[2722] + v[1495] * v[2723] - v[2715] * v[3330];
v[3326] = v[2730] + (-(v[172] * v[2664]) + v[2681] * v[316] - v[2673] * v[317] + v[2676] * v[318]) / 2e0;
v[2724] = (v[172] * v[2651] - v[2699] * v[316] + v[2702] * v[317] - v[2690] * v[318]) / 2e0;
v[3325] = v[2722] - v[2724] + v[3326];
v[3323] = -v[2724] + v[2730];
v[2725] = v[2707] + v[2711];
v[2726] = v[2710] + v[2712];
v[2727] = v[2696] + v[2708];
v[6090] = v[128] * v[2638];
v[6091] = v[128] * v[2637];
v[6092] = v[128] * v[2636];
v[6093] = v[130] * v[2638];
v[6094] = v[130] * v[2637];
v[6095] = v[130] * v[2636];
v[6096] = v[132] * v[2638];
v[6097] = v[132] * v[2637];
v[6098] = v[132] * v[2636];
v[6099] = -v[2638];
v[6100] = -v[2637];
v[6101] = -v[2636];
v[6102] = -v[2710] + v[2712] + v[2727] * v[311] + v[2725] * v[314] + v[310] * v[3323] + v[2695] * v[3332] + 2e0*
(v[2719] * v[3330] + v[2713] * v[3332]);
v[6103] = v[2707] - v[2711] + v[2727] * v[312] + v[2726] * v[314] + v[313] * v[3325] - v[2691] * v[3332] + 2e0*(-
(v[2721] * v[3330]) + v[2714] * v[3332]);
v[6104] = -v[2696] + v[2708] + v[2726] * v[311] + v[2725] * v[312] + v[315] * v[3326] + v[2706] * v[3332] + 2e0*
(v[2723] * v[3330] + v[2709] * v[3332]);
v[2728] = (-(v[220] * v[2605]) + v[154] * v[2636] + v[150] * v[2637] + v[146] * v[2638] - v[2639] + v[2597] * v[3361]
	+ v[219] * v[3609] + v[218] * v[3610]) / 2e0;
v[2729] = (-(v[217] * v[2605]) + v[153] * v[2636] + v[149] * v[2637] + v[145] * v[2638] - v[2639] + v[2596] * v[3361]
	+ v[216] * v[3609] + v[215] * v[3610]) / 2e0;
for (i2518 = 1; i2518 <= 15; i2518++) {
	v[2767] = (i2518 == 14 ? 1 : 0);
	v[3688] = v[10] * v[2767];
	v[2764] = (i2518 == 13 ? 1 : 0);
	v[3689] = v[10] * v[2764];
	v[2761] = (i2518 == 15 ? 1 : 0);
	v[3693] = v[10] * v[2761];
	v[2746] = v[4948 + i2518];
	v[2745] = v[4963 + i2518];
	v[2736] = v[4933 + i2518];
	v[2759] = -v[2736] / 2e0;
	v[3680] = -(v[220] * v[2759]);
	v[3671] = v[219] * v[2759];
	v[3666] = v[218] * v[2759];
	v[2735] = v[4918 + i2518];
	v[2737] = v[5012 + i2518];
	v[2738] = v[5042 + i2518];
	v[2739] = v[5027 + i2518];
	v[2741] = v[5545 + i2518];
	v[2742] = v[4997 + i2518];
	v[2777] = v[2742] * v[3330];
	v[2819] = -(v[2777] * v[3443]);
	v[3699] = v[238] * v[2819];
	v[3697] = v[236] * v[2819];
	v[3613] = v[227] * v[2819];
	v[2744] = v[5605 + i2518];
	v[3611] = -v[2745] - v[2746];
	v[2748] = (i2518 == 1 ? v[10] : 0e0);
	v[2749] = (i2518 == 2 ? v[10] : 0e0);
	v[2750] = (i2518 == 4 ? v[10] : 0e0);
	v[2751] = (i2518 == 5 ? v[10] : 0e0);
	v[2752] = (i2518 == 7 ? v[10] : 0e0);
	v[2753] = (i2518 == 8 ? v[10] : 0e0);
	v[2754] = (i2518 == 10 ? v[10] : 0e0);
	v[2755] = (i2518 == 11 ? v[10] : 0e0);
	v[3641] = v[2523] * v[2755];
	v[2756] = (i2518 == 12 ? v[10] : 0e0);
	v[3644] = v[2523] * v[2756];
	v[2757] = v[2735] / 2e0;
	v[3678] = v[217] * v[2757];
	v[3669] = -(v[216] * v[2757]);
	v[3664] = -(v[215] * v[2757]);
	v[2758] = -v[2757] + v[2759];
	v[3679] = v[214] * v[2758];
	v[3670] = -(v[213] * v[2758]);
	v[3665] = -(v[212] * v[2758]);
	v[3620] = v[3678] + v[3679] + v[3680];
	v[3618] = v[2596] * v[2757] + v[2595] * v[2758] - v[2597] * v[2759];
	v[2760] = v[2737] + v[2761];
	v[3702] = 2e0*v[2760];
	v[2762] = v[2737] - v[2761];
	v[3703] = 2e0*v[2762];
	v[2763] = v[2738] + v[2764];
	v[3704] = 2e0*v[2763];
	v[2765] = v[2738] - v[2764];
	v[3705] = 2e0*v[2765];
	v[2766] = v[2739] - v[2767];
	v[3706] = 2e0*v[2766];
	v[2768] = v[2739] + v[2767];
	v[3707] = 2e0*v[2768];
	v[2769] = v[1495] * v[2742] + v[2761] * v[3450];
	v[2770] = -v[2742] + v[2767] * v[313];
	v[2771] = v[1498] * v[2742] - v[2767] * v[3450];
	v[2772] = v[1500] * v[2742] + v[2764] * v[3450];
	v[3612] = 2e0*(-(v[172] * v[2767]) - v[2742] * v[3331]);
	v[2774] = -(v[172] * v[2764]) - v[2742] * v[3333];
	v[2775] = -(v[172] * v[2761]) - v[2742] * v[3334];
	v[2776] = v[2777] * v[314] + v[2761] * v[3332];
	v[2778] = v[2777] * v[312] + v[2764] * v[3332];
	v[2779] = -(v[2777] * v[311]) - v[2767] * v[3332];
	v[2780] = (v[172] * v[2741] - v[2777] * v[356]) / 2e0;
	v[2854] = v[238] * v[2780];
	v[2781] = (-(v[2741] * v[318]) - v[2769] * v[356]) / 2e0;
	v[2782] = (v[172] * v[2770] - v[2777] * v[338]) / 2e0;
	v[2848] = v[236] * v[2782];
	v[2783] = (v[2770] * v[317] + v[2771] * v[338]) / 2e0;
	v[2784] = (v[172] * v[2744] - v[2777] * v[319]) / 2e0;
	v[2843] = v[230] * v[2784];
	v[2785] = (-(v[2744] * v[316]) - v[2772] * v[319]) / 2e0;
	v[2786] = (v[3611] * v[76] + v[2746] * v[79] + v[2745] * v[82]) / 2e0;
	v[2787] = (v[3611] * v[77] + v[2746] * v[80] + v[2745] * v[83]) / 2e0;
	v[2788] = (v[3611] * v[78] + v[2746] * v[81] + v[2745] * v[84]) / 2e0;
	v[2789] = (v[2741] * v[317] + v[2771] * v[356] + v[3612]) / 2e0;
	v[2790] = (v[2744] * v[317] + v[2771] * v[319] + v[3612]) / 2e0;
	v[2791] = v[2774] + v[2741] * v[3333] - v[2772] * v[3426];
	v[2792] = v[2774] + v[2770] * v[3333] - v[2772] * v[3427];
	v[2793] = -v[2777] - v[2763] * v[316] - v[2772] * v[351];
	v[2794] = v[172] * v[2763] - v[2777] * v[351];
	v[2795] = v[2777] + v[2766] * v[317] + v[2771] * v[347];
	v[2796] = v[172] * v[2766] - v[2777] * v[347];
	v[2856] = v[230] * v[2796];
	v[2797] = v[2777] - v[2765] * v[316] - v[2772] * v[343];
	v[2798] = v[172] * v[2765] - v[2777] * v[343];
	v[2849] = v[238] * v[2798];
	v[2799] = v[2775] + v[2770] * v[3334] - v[2769] * v[3427];
	v[2800] = v[2775] + v[2744] * v[3334] - v[2769] * v[3428];
	v[2801] = -v[2777] - v[2760] * v[318] - v[2769] * v[334];
	v[2802] = v[172] * v[2760] - v[2777] * v[334];
	v[2803] = -v[2777] + v[2768] * v[317] + v[2771] * v[330];
	v[2804] = v[172] * v[2768] - v[2777] * v[330];
	v[2844] = v[238] * v[2804];
	v[2805] = -v[2776] + v[2763] * v[317] + v[2771] * v[351];
	v[2806] = -v[2776] - v[2766] * v[316] - v[2772] * v[347];
	v[2807] = -v[2776] + v[2765] * v[317] + v[2771] * v[343];
	v[2808] = -v[2776] - v[2768] * v[316] - v[2772] * v[330];
	v[2809] = v[2819] + v[2776] * v[3452];
	v[2810] = v[2786] * v[392] + v[2787] * v[395] + v[2788] * v[398] - v[2789] * v[738] - v[2805] * v[740] - v[2795] * v[742];
	v[3657] = -(v[281] * v[2810]);
	v[2811] = v[2786] * v[391] + v[2787] * v[394] + v[2788] * v[397] - v[2791] * v[738] - v[2793] * v[740] - v[2806] * v[742];
	v[3629] = -(v[281] * v[2811]);
	v[2812] = v[2777] - v[2762] * v[318] - v[2769] * v[324];
	v[2813] = v[172] * v[2762] - v[2777] * v[324];
	v[2814] = -v[2778] + v[2760] * v[317] + v[2771] * v[334];
	v[2815] = -v[2778] - v[2766] * v[318] - v[2769] * v[347];
	v[2816] = -v[2778] - v[2768] * v[318] - v[2769] * v[330];
	v[2817] = -v[2778] + v[2762] * v[317] + v[2771] * v[324];
	v[2818] = v[2776] * v[326] + v[2778] * v[329];
	v[2820] = v[2819] + v[2778] * v[3453];
	v[2870] = v[2604] * v[2820];
	v[2821] = v[2786] * v[377] + v[2787] * v[380] + v[2788] * v[383] - v[2807] * v[738] - v[2783] * v[740] - v[2814] * v[742];
	v[3658] = -(v[280] * v[2821]);
	v[2822] = v[2779] - v[2763] * v[318] - v[2769] * v[351];
	v[2823] = v[2779] - v[2765] * v[318] - v[2769] * v[343];
	v[2824] = v[2779] - v[2760] * v[316] - v[2772] * v[334];
	v[2825] = v[2779] - v[2762] * v[316] - v[2772] * v[324];
	v[2826] = -(v[2778] * v[323]) - v[2779] * v[326];
	v[2827] = -(v[2776] * v[323]) - v[2779] * v[329];
	v[2828] = v[2819] + v[2779] * v[3454];
	v[2874] = v[2603] * v[2828];
	v[2829] = v[2786] * v[393] + v[2787] * v[396] + v[2788] * v[399] - v[2781] * v[738] - v[2822] * v[740] - v[2815] * v[742];
	v[3654] = -(v[281] * v[2829]);
	v[2830] = v[2786] * v[378] + v[2787] * v[381] + v[2788] * v[384] - v[2823] * v[738] - v[2799] * v[740] - v[2801] * v[742];
	v[3655] = -(v[280] * v[2830]);
	v[2831] = v[2786] * v[376] + v[2787] * v[379] + v[2788] * v[382] - v[2797] * v[738] - v[2792] * v[740] - v[2824] * v[742];
	v[3630] = -(v[280] * v[2831]);
	v[3634] = v[3629] + v[3630];
	v[2832] = v[2786] * v[361] + v[2787] * v[364] + v[2788] * v[367] - v[2808] * v[738] - v[2825] * v[740] - v[2785] * v[742];
	v[3631] = v[279] * v[2832];
	v[2833] = v[2786] * v[363] + v[2787] * v[366] + v[2788] * v[369] - v[2816] * v[738] - v[2812] * v[740] - v[2800] * v[742];
	v[3656] = -(v[279] * v[2833]);
	v[2834] = v[2786] * v[362] + v[2787] * v[365] + v[2788] * v[368] - v[2803] * v[738] - v[2817] * v[740] - v[2790] * v[742];
	v[3659] = -(v[279] * v[2834]);
	v[2835] = v[2771] + v[2818];
	v[2868] = v[2604] * v[2835];
	v[2836] = -v[2771] + v[2818];
	v[2837] = (v[225] * v[2794] + v[186] * v[2835] + v[1593] * v[3613]) / v[227];
	v[2838] = v[2769] + v[2826];
	v[2876] = v[2604] * v[2838];
	v[2839] = -v[2769] + v[2826];
	v[2872] = v[2603] * v[2839];
	v[2840] = (v[226] * v[2798] + v[182] * v[2838] + v[1597] * v[3613]) / v[227];
	v[2841] = v[236] * v[2813] + v[2843];
	v[2842] = v[2841] + v[2844] + v[3689];
	v[2845] = v[2843] + v[2844];
	v[2846] = v[230] * v[2802] + v[2848];
	v[2847] = v[2846] + v[2849] + v[3688];
	v[2850] = v[2848] + v[2849];
	v[2851] = v[2637] * v[2782] - v[2621] * v[2783] - v[2624] * v[2792] - v[2625] * v[2793] + v[2636] * v[2794]
		- v[2618] * v[2799] - v[2622] * v[2805] - v[2617] * v[2812] + v[2638] * v[2813] - v[2620] * v[2817] - v[2619] * v[2822]
		- v[2623] * v[2825];
	v[2852] = v[2638] * v[2784] - v[2623] * v[2785] - v[2620] * v[2790] - v[2622] * v[2795] + v[2636] * v[2796]
		- v[2617] * v[2800] - v[2618] * v[2801] + v[2637] * v[2802] - v[2625] * v[2806] - v[2621] * v[2814] - v[2619] * v[2815]
		- v[2624] * v[2824];
	v[2853] = v[2854] + v[2856];
	v[2855] = v[236] * v[2794] + v[2854];
	v[2857] = v[2855] + v[2856] + v[3693];
	v[2858] = v[2636] * v[2780] - v[2619] * v[2781] - v[2622] * v[2789] - v[2625] * v[2791] - v[2624] * v[2797]
		+ v[2637] * v[2798] - v[2620] * v[2803] + v[2638] * v[2804] - v[2621] * v[2807] - v[2623] * v[2808] - v[2617] * v[2816]
		- v[2618] * v[2823];
	v[2859] = (v[177] * v[2839] + v[2804] * v[3340] + v[1599] * v[3613]) / v[227];
	v[2860] = v[2870] + v[2872];
	v[3722] = v[2860] / v[227];
	v[2861] = v[2772] + v[2827];
	v[2866] = v[2603] * v[2861];
	v[2862] = -v[2772] + v[2827];
	v[2863] = v[2874] + v[2876];
	v[3718] = v[2863] / v[227];
	v[2864] = v[2602] * v[2809] + v[2866] + v[2868];
	v[2867] = v[2864] - v[2868];
	v[2869] = v[2864] - v[2866];
	v[3719] = v[2869] / v[227];
	v[2871] = v[2602] * v[2836] + v[2870];
	v[2873] = v[2871] + v[2872];
	v[2875] = v[2602] * v[2862] + v[2874];
	v[2877] = v[2875] + v[2876];
	v[2879] = v[145] * v[2757] + v[144] * v[2758] - v[146] * v[2759] - v[191] * v[2786] - v[192] * v[2787] - v[193] * v[2788]
		+ v[5192 + i2518] + v[2804] * v[738] + v[2813] * v[740] + v[2784] * v[742];
	v[2881] = v[149] * v[2757] + v[148] * v[2758] - v[150] * v[2759] - v[194] * v[2786] - v[195] * v[2787] - v[196] * v[2788]
		+ v[5207 + i2518] + v[2798] * v[738] + v[2782] * v[740] + v[2802] * v[742];
	v[2883] = v[153] * v[2757] + v[152] * v[2758] - v[154] * v[2759] - v[197] * v[2786] - v[198] * v[2787] - v[199] * v[2788]
		+ v[5222 + i2518] + v[2780] * v[738] + v[2794] * v[740] + v[2796] * v[742];
	v[3614] = v[250] * v[2879] + v[251] * v[2881] + v[252] * v[2883];
	v[2886] = v[3614] / v[602];
	v[2884] = -(v[2536] * v[3614] * v[820]);
	v[2885] = v[263] * v[2886];
	v[2888] = v[2879];
	v[3041] = v[2888];
	v[2889] = v[264] * v[2879] + v[250] * v[2885];
	v[2890] = v[2881];
	v[3040] = v[2890];
	v[2891] = v[264] * v[2881] + v[251] * v[2885];
	v[2892] = v[2883];
	v[3039] = v[2892];
	v[2893] = v[264] * v[2883] + v[252] * v[2885];
	b2894 = b277;
	if (b2894) {
		v[2895] = -v[2889];
		v[2896] = -v[2891];
		v[2897] = -v[2893];
	}
	else {
		v[2895] = v[2889];
		v[2896] = v[2891];
		v[2897] = v[2893];
	};
	v[3633] = -(v[281] * v[2895]);
	v[3632] = v[2895] * v[406];
	v[3628] = v[280] * v[2895];
	v[2994] = v[2590] * v[2895];
	v[2992] = v[2591] * v[2895];
	v[2990] = v[2592] * v[2895];
	v[3625] = -(v[281] * v[2896]);
	v[3624] = -(v[279] * v[2896]);
	v[2959] = v[2566] * v[2896];
	v[2957] = v[2567] * v[2896];
	v[2955] = v[2568] * v[2896];
	v[3640] = -(v[223] * v[2897]);
	v[3623] = v[128] * v[2897];
	v[3622] = v[130] * v[2897];
	v[3621] = v[132] * v[2897];
	v[3616] = -(v[280] * v[2897]);
	v[3615] = -(v[279] * v[2897]);
	v[2898] = v[1637] * v[2897] + v[281] * v[3620] + v[10] * v[6168 + i2518];
	v[2906] = v[2897] * v[2905];
	v[2908] = v[2897] * v[2907];
	v[2909] = v[2897] * v[3466];
	v[3660] = v[2525] * v[2909];
	v[2910] = v[2897] * v[3822];
	v[2911] = v[2897] * v[3823];
	v[2912] = v[2897] * v[3824];
	v[2913] = v[2541] * v[2897];
	v[2914] = v[2540] * v[2897];
	v[2915] = v[2539] * v[2897];
	v[2916] = v[2897] * v[3617] + v[280] * v[3618];
	v[2920] = v[279] * v[3618] + v[2897] * v[3619];
	v[2922] = 2e0*v[2897] * v[2921] + v[2599] * v[3620] + v[10] * v[6228 + i2518];
	v[2932] = v[2897] * v[3826];
	v[2933] = v[2897] * v[3827];
	v[2934] = v[2897] * v[3828];
	v[2935] = v[2897] * v[3831] + v[10] * v[6348 + i2518] + v[6363 + i2518];
	v[2936] = v[2897] * v[3834] + v[10] * v[6378 + i2518] + v[6393 + i2518];
	v[2937] = v[2897] * v[3837] + v[10] * v[6408 + i2518] + v[6423 + i2518];
	v[2938] = v[279] * (-(v[2759] * v[280]) + v[132] * v[2896]);
	v[2939] = v[279] * (v[2757] * v[280] + v[130] * v[2896]);
	v[2940] = v[279] * (v[2758] * v[280] + v[128] * v[2896]);
	v[2944] = v[2896] * v[3469];
	v[3662] = v[2524] * v[2944];
	v[2945] = v[2896] * v[3838] + v[3616] * v[412];
	v[2946] = v[2896] * v[3839] + v[3616] * v[413];
	v[2947] = v[2896] * v[3840] + v[3616] * v[414];
	v[2951] = v[2913] + v[2955];
	v[2952] = v[2914] + v[2957];
	v[2953] = v[2915] + v[2959];
	v[2916] = v[2916] + v[2896] * v[3626];
	v[2920] = v[2920] + 2e0*v[2896] * v[2954];
	v[2922] = v[2922] + v[2896] * v[3627];
	v[2962] = v[2896] * v[3841];
	v[2964] = v[2896] * v[3842];
	v[2966] = v[2896] * v[3843];
	v[2970] = v[2938] + v[132] * v[3628];
	v[2971] = v[2939] + v[130] * v[3628];
	v[2972] = v[2940] + v[128] * v[3628];
	v[2973] = v[2895] * v[3472];
	v[3663] = v[2523] * v[2973];
	v[2974] = v[279] * v[3634] + v[2895] * v[3844] + v[3624] * v[409] + v[3615] * v[412] - v[2832] * v[582];
	v[2975] = v[2895] * v[3845] + v[3624] * v[410] + v[3615] * v[413];
	v[2976] = v[2895] * v[3846] + v[3624] * v[411] + v[3615] * v[414];
	v[2916] = v[2916] + 2e0*v[2895] * v[2985] + v[2592] * v[3634] - v[2832] * v[3635];
	v[2986] = v[2913] + v[2990];
	v[2987] = v[2914] + v[2992];
	v[2988] = v[2915] + v[2994];
	v[2920] = v[2920] + v[2568] * (v[3629] - v[3631]) - v[2831] * v[3636] + v[2895] * v[3637];
	v[2922] = v[2922] + v[2541] * v[3630] - v[2541] * v[3631] + v[2895] * v[3638] - v[2811] * v[3639];
	v[2998] = v[2895] * v[3847];
	v[3001] = v[2895] * v[3848];
	v[3002] = v[2895] * v[3849];
	v[3006] = v[1787] * v[2819] + v[1788] * v[2820] + v[1789] * v[2835] + v[236] * v[2837] + v[1790] * v[2838]
		+ v[238] * v[2840] + v[2507] * v[2842] + v[2846] * v[3341] + v[2853] * v[3344] + v[10] * v[6123 + i2518];
	v[3007] = v[1792] * v[2819] + v[1793] * v[2828] + v[230] * v[2837] + v[1794] * v[2839] + v[1438] * v[2841]
		+ v[2510] * v[2847] + v[1437] * v[2855] + v[238] * v[2859] + v[1795] * v[2861] + v[10] * v[6138 + i2518];
	v[3008] = v[1798] * v[2809] + v[1797] * v[2819] + v[1799] * v[2836] + v[230] * v[2840] + v[1440] * v[2845]
		+ v[1439] * v[2850] + v[2512] * v[2857] + v[236] * v[2859] + v[1800] * v[2862] + v[10] * v[6153 + i2518];
	v[2916] = v[2916] + v[2834] * v[2996] + v[2834] * v[2997] + v[2833] * v[2999] + v[2833] * v[3000] + v[2523] * v[3640]
		- v[280] * (v[2591] * v[2821] + v[2590] * v[2830] + v[3641]) - v[281] * (v[2591] * v[2810] + v[2590] * v[2829] + v[3644]
			);
	v[3009] = -(v[2897] * v[3597]);
	v[2916] = v[2523] * (-(v[222] * v[2896]) + v[2898]) + v[2916];
	v[3010] = -(v[2896] * v[3597]);
	v[3011] = v[2523] * v[2895];
	v[3013] = v[2754] + v[3664] + v[3665] + v[3666] + v[10] * v[6183 + i2518];
	v[3014] = v[2755] + v[3669] + v[3670] + v[3671] + v[10] * v[6198 + i2518];
	v[3016] = v[3009] - v[2897] * v[3596];
	v[3021] = v[2525] * v[2897];
	v[3017] = v[2524] * v[2896];
	v[3018] = v[3011] + v[3017];
	v[2920] = v[281] * (-(v[2567] * v[2810]) - v[2566] * v[2829]) + v[2920] + v[2524] * (-(v[2754] * v[279])
		- v[2756] * v[281] - v[221] * v[2895] + v[2898] + v[3640]) + v[279] * (-(v[2566] * v[2833]) - v[2567] * v[2834] - v[3641]
			) - v[2821] * v[3642] - v[2830] * v[3643];
	v[3019] = -(v[2895] * v[3596]) - v[3663];
	v[2922] = v[280] * (-(v[2540] * v[2821]) - v[2539] * v[2830]) + v[2922] - v[2544] * v[3013] - v[2547] * v[3014]
		- v[2756] * v[3596] + v[279] * (-(v[2539] * v[2833]) - v[2540] * v[2834] - v[3644]) - v[2810] * v[3645]
		- v[2829] * v[3646];
	v[3058] = v[2922];
	v[3020] = v[3017] + v[3021];
	v[3022] = v[3011] + v[3021];
	v[2920] = v[2920] + v[2525] * (v[2906] - v[281] * v[3014]);
	v[3060] = v[2920];
	v[2916] = v[2916] + v[2525] * (v[2908] - v[281] * v[3013]) - v[2754] * v[3596];
	v[3062] = v[2916];
	v[3023] = v[2525] * v[2896];
	v[3024] = v[2525] * v[2895];
	b3025 = b603;
	if (b3025) {
		b3026 = b605;
		if (b3026) {
			v[2895] = 0e0;
			v[2896] = 0e0;
			v[2897] = 0e0;
		}
		else {
		};
	}
	else {
	};
	v[3027] = 0e0;
	v[3028] = 0e0;
	v[3029] = 0e0;
	v[3030] = 0e0;
	v[3031] = 0e0;
	v[3032] = 0e0;
	v[3033] = 0e0;
	b3034 = b603;
	if (b3034) {
		v[3035] = 0e0;
		v[3036] = 0e0;
		v[3037] = 0e0;
		b3038 = b622;
		if (b3038) {
			v[3037] = v[2892];
			v[2892] = 0e0;
			v[3036] = v[2890];
			v[2890] = 0e0;
			v[3035] = v[2888];
			v[2888] = 0e0;
		}
		else {
			v[3032] = -v[3039];
			v[2892] = 0e0;
			v[3031] = -v[3040];
			v[2890] = 0e0;
			v[3030] = -v[3041];
			v[2888] = 0e0;
		};
		v[3647] = (v[3035] * v[583] + v[3036] * v[591] + v[3037] * v[593])*v[7];
		b3042 = b605;
		if (b3042) {
			v[3029] = v[3037] * v[617];
			v[3028] = v[3036] * v[617];
			v[3027] = v[3035] * v[617];
			v[3033] = (v[1064] * v[3647] * v[615] * Power(v[614], v[1987])) / v[1065];
		}
		else {
			v[3029] = v[3037] * v[621];
			v[3028] = v[3036] * v[621];
			v[3027] = v[3035] * v[621];
			v[2884] = v[2884] + (v[1072] * v[3647] * v[620] * Power(v[602], v[1990])) / v[1073];
		};
	}
	else {
	};
	v[3686] = v[3027] * v[582];
	v[3653] = -(v[279] * v[3027]);
	v[3687] = v[3028] * v[587];
	v[3652] = v[280] * v[3028];
	v[3661] = v[3029] * v[592];
	v[3651] = v[3660] + v[3661];
	v[3063] = v[2884];
	v[3061] = v[3030];
	v[3059] = v[3031];
	v[3057] = v[3032];
	b3052 = b603;
	if (b3052) {
		v[3650] = v[279] * v[3030];
		v[3649] = v[280] * v[3031];
		v[3648] = v[281] * v[3032];
		b3053 = b605;
		if (b3053) {
			v[2922] = v[2922] + v[3032] * v[3368];
			v[3032] = 0e0;
			v[2920] = v[2920] + v[3031] * v[3368];
			v[3031] = 0e0;
			v[2916] = v[2916] + v[3030] * v[3368];
			v[3030] = 0e0;
			v[3033] = v[3033] + v[31] * v[34] * (-v[3648] - v[3649] - v[3650])*Power(v[614], v[615]);
			v[2884] = v[2884] - v[3033];
		}
		else {
			v[2922] = v[3058] + v[3057] * v[611];
			v[3032] = 0e0;
			v[2920] = v[3060] + v[3059] * v[611];
			v[3031] = 0e0;
			v[2916] = v[3062] + v[3061] * v[611];
			v[3030] = 0e0;
			v[2884] = v[3063] + v[3382] * (v[3648] + v[3649] + v[3650])*Power(v[602], v[620]);
		};
	}
	else {
	};
	v[3064] = v[2525] * v[3008] + v[249] * v[3029];
	v[3673] = v[281] * v[3064];
	v[3065] = v[2525] * v[3007] + v[239] * v[3029];
	v[3668] = v[281] * v[3065];
	v[3066] = v[2525] * v[3006] + v[229] * v[3029];
	v[3667] = v[281] * v[3066];
	v[2916] = v[2916] + v[2178] * v[3029];
	v[2920] = v[2920] + v[2179] * v[3029];
	v[3068] = v[279] * v[3029];
	v[3071] = v[280] * v[3029];
	v[2922] = v[2922] + v[2550] * v[3029];
	v[3089] = v[2524] * v[3008] + v[249] * v[3028];
	v[3681] = v[280] * v[3089];
	v[3090] = v[2524] * v[3007] + v[239] * v[3028];
	v[3676] = v[280] * v[3090];
	v[3091] = v[2524] * v[3006] + v[229] * v[3028];
	v[3674] = v[280] * v[3091];
	v[2916] = v[2916] - v[221] * v[3652];
	v[2920] = v[2920] + v[2571] * v[3028];
	v[2922] = v[2922] - v[223] * v[3652];
	v[3111] = v[2523] * v[3008] + v[249] * v[3027];
	v[3682] = v[279] * v[3111];
	v[3685] = -v[2959] - v[2994] - v[3681] - v[3682];
	v[3112] = v[2523] * v[3007] + v[239] * v[3027];
	v[3677] = v[279] * v[3112];
	v[3683] = v[2957] + v[2992] + v[3676] + v[3677];
	v[3113] = v[2523] * v[3006] + v[229] * v[3027];
	v[3675] = v[279] * v[3113];
	v[3684] = v[2955] + v[2990] + v[3674] + v[3675];
	v[3115] = v[2524] * v[2748] + v[2523] * v[2749] + v[213] * v[3027] + v[212] * v[3028];
	v[3116] = v[2524] * v[2750] + v[2523] * v[2751] + v[216] * v[3027] + v[215] * v[3028];
	v[3117] = v[2524] * v[2752] + v[2523] * v[2753] + v[219] * v[3027] + v[218] * v[3028];
	v[3672] = v[128] * v[3115] + v[130] * v[3116] + v[132] * v[3117];
	v[2916] = v[2916] + v[2598] * v[3027];
	v[2920] = v[2920] + v[222] * v[3653];
	v[2922] = v[2922] + v[223] * v[3653];
	v[3121] = v[1171] * v[3027] + v[1129] * v[3028] + v[1086] * v[3029] + v[2523] * (v[2976] + v[279] * (v[3654] + v[3655])
		- v[2833] * v[582]) + v[2524] * (v[2947] + v[280] * (v[3654] + v[3656]) - v[3628] * v[408] - v[2830] * v[587]) + v[2525] *
		(v[2912] + v[281] * (v[3655] + v[3656]) + v[3633] * v[408] + v[3625] * v[411] - v[2829] * v[592]);
	v[3211] = v[3121] * v[3337];
	v[3122] = v[1173] * v[3027] + v[1131] * v[3028] + v[1088] * v[3029] + v[2523] * (v[2975] + v[279] * (v[3657] + v[3658])
		- v[2834] * v[582]) + v[2524] * (v[2946] + v[280] * (v[3657] + v[3659]) - v[3628] * v[407] - v[2821] * v[587]) + v[2525] *
		(v[2911] + v[281] * (v[3658] + v[3659]) + v[3633] * v[407] + v[3625] * v[410] - v[2810] * v[592]);
	v[3695] = v[227] * v[3122];
	v[3215] = v[3122] * v[3336];
	v[3123] = v[2523] * v[2974] + v[1175] * v[3027] + v[1133] * v[3028] + v[1090] * v[3029] + v[2524] * (v[2945] - v[280] * (
		-v[3629] + v[3631] + v[3632]) - v[2831] * v[587]) + v[2525] * (v[2910] + v[281] * (v[3630] - v[3631] - v[3632])
			+ v[3625] * v[409] - v[2811] * v[592]);
	v[3696] = v[227] * v[3123];
	v[3218] = v[3123] * v[3335];
	v[3136] = v[3001] - v[279] * (v[2953] + v[3673] + v[3681]) - v[3111] * v[582];
	v[3137] = v[2998] - v[279] * (v[2952] + v[3668] + v[3676]) - v[3112] * v[582];
	v[2916] = v[2916] - v[221] * v[3020] + v[1832] * v[3111] + v[1834] * v[3112] + v[1836] * v[3113] - v[2951] * v[406]
		- v[2952] * v[407] - v[2953] * v[408] + v[281] * (-(v[3066] * v[406]) - v[3065] * v[407] - v[3064] * v[408]) + v[280] *
		(v[3672] - v[3091] * v[406] - v[3090] * v[407] - v[3089] * v[408]) + v[3472] * (-(v[2592] * v[2832]) - v[2590] * v[2833]
			- v[2591] * v[2834] + v[2593] * v[3027] + v[2523] * (v[128] * v[2748] + v[130] * v[2750] + v[132] * v[2752] - v[2754]
				- v[3664] - v[3665] - v[3666]) - v[3113] * v[406] - v[3112] * v[407] - v[3111] * v[408]);
	v[3138] = v[3002] - v[279] * (v[2951] + v[3667] + v[3674]) - v[3113] * v[582];
	v[3139] = v[2962] - v[280] * (v[2986] + v[3667] + v[3675]) - v[3091] * v[587];
	v[3140] = v[2966] - v[280] * (v[2987] + v[3668] + v[3677]) - v[3090] * v[587];
	v[2920] = v[2920] - v[222] * v[3022] + v[3089] * v[3356] + v[3090] * v[3358] + v[3091] * v[3360] - v[2986] * v[409]
		- v[2987] * v[410] - v[2988] * v[411] + v[281] * (-(v[3066] * v[409]) - v[3065] * v[410] - v[3064] * v[411]) + v[3469] * (-
		(v[2567] * v[2821]) - v[2566] * v[2830] - v[2568] * v[2831] + v[2569] * v[3028] + v[2524] * (v[128] * v[2749]
			+ v[130] * v[2751] + v[132] * v[2753] - v[2755] - v[3669] - v[3670] - v[3671]) - v[3091] * v[409] - v[3090] * v[410]
			- v[3089] * v[411]) + v[279] * (v[3672] - v[3113] * v[409] - v[3112] * v[410] - v[3111] * v[411]);
	v[3141] = v[2964] - v[280] * (v[2988] + v[3673] + v[3682]) - v[3089] * v[587];
	v[3142] = v[3018] + v[3652] - v[3653];
	v[3143] = v[3024] + v[3068];
	v[3144] = v[3023] + v[3071];
	v[6559] = -(v[2556] * v[2758]) + v[2524] * v[2972] + v[2544] * v[3623] + v[128] * v[3663] + v[3143] * v[3714] +
		(v[2523] * v[2758] + v[128] * v[3027])*v[582] + v[3028] * v[584];
	v[6560] = -(v[2555] * v[2758]) + v[2523] * v[2972] + v[2547] * v[3623] + v[128] * v[3662] + v[3144] * v[3714]
		+ v[3027] * v[584] + (v[2524] * v[2758] + v[128] * v[3028])*v[587];
	v[6561] = -(v[2554] * v[2758]) + v[2599] * (v[2758] * v[281] + v[3623]) + v[128] * v[3651] + v[3142] * v[3714];
	v[6562] = -(v[2556] * v[2757]) + v[2524] * v[2971] + v[2544] * v[3622] + v[130] * v[3663] + v[3143] * v[3713] +
		(v[2523] * v[2757] + v[130] * v[3027])*v[582] + v[3028] * v[585];
	v[6563] = -(v[2555] * v[2757]) + v[2523] * v[2971] + v[2547] * v[3622] + v[130] * v[3662] + v[3144] * v[3713]
		+ v[3027] * v[585] + (v[2524] * v[2757] + v[130] * v[3028])*v[587];
	v[6564] = -(v[2554] * v[2757]) + v[2599] * (v[2757] * v[281] + v[3622]) + v[130] * v[3651] + v[3142] * v[3713];
	v[6565] = v[2556] * v[2759] + v[2524] * v[2970] + v[2544] * v[3621] + v[132] * v[3663] + v[3143] * v[3712] + (-
		(v[2523] * v[2759]) + v[132] * v[3027])*v[582] + v[3028] * v[586];
	v[6566] = v[2555] * v[2759] + v[2523] * v[2970] + v[2547] * v[3621] + v[132] * v[3662] + v[3144] * v[3712]
		+ v[3027] * v[586] + (-(v[2524] * v[2759]) + v[132] * v[3028])*v[587];
	v[6567] = v[2554] * v[2759] + v[2599] * (-(v[2759] * v[281]) + v[3621]) + v[132] * v[3651] + v[3142] * v[3712];
	v[6568] = v[3019] + v[281] * (-v[3024] - v[3068]) + v[279] * (-v[3020] - v[3652]) - v[3686];
	v[6569] = v[3010] - v[280] * v[3022] + v[281] * (-v[3023] - v[3071]) - v[3027] * v[3361] - v[3662] - v[3687];
	v[6570] = v[3016] + v[281] * (-v[3018] - v[3652] + v[3653]) - v[3660] - v[3661];
	v[6571] = v[2603] * v[2837] + v[2602] * v[2840] + v[2507] * v[3123] + v[2819] * (v[2487] * v[2602] + v[2489] * v[2603]
		+ v[3221]) + v[3122] * v[3341] + v[3121] * v[3344] + v[2796] * v[3604] + v[2802] * v[3605] + v[2873] * v[3691]
		+ v[2784] * v[3717] + v[179] * v[3718] + v[184] * v[3719];
	v[6572] = v[2604] * v[2837] + v[2602] * v[2859] + (v[186] * v[2867]) / v[227] + v[1437] * v[3121] + v[2510] * v[3122]
		+ v[1438] * v[3123] + v[2819] * (v[2486] * v[2602] + v[2490] * v[2604] + v[3214]) + v[2813] * v[3606] + v[2877] * v[3690]
		+ v[2782] * v[3720] + v[2794] * v[3721] + v[176] * v[3722];
	v[6573] = v[2604] * v[2840] + v[2603] * v[2859] + (v[177] * v[2871]) / v[227] + (v[182] * v[2875]) / v[227]
		+ v[2512] * v[3121] + v[1439] * v[3122] + v[1440] * v[3123] + v[2819] * (v[2488] * v[2603] + v[2491] * v[2604] + v[3206])
		+ v[2864] * v[3694] + v[2780] * v[3723] + v[2798] * v[3724] + v[2804] * v[3725];
	v[3145] = v[2933] - v[281] * v[3684] - v[3066] * v[592];
	v[3146] = v[2934] - v[281] * v[3683] - v[3065] * v[592];
	v[2922] = v[2922] - v[223] * v[3018] - v[222] * v[3023] - v[221] * v[3024] + v[1923] * v[3064] + v[1925] * v[3065]
		+ v[1927] * v[3066] + v[1637] * v[3142] + v[2907] * v[3143] + v[2905] * v[3144] - v[3684] * v[412] - v[3683] * v[413]
		+ v[3685] * v[414] + v[3466] * (-(v[2540] * v[2810]) - v[2541] * v[2811] - v[2539] * v[2829] + v[2542] * v[3029]
			- v[3066] * v[412] - v[3065] * v[413] - v[3064] * v[414] - v[2525] * (v[2756] - v[3678] - v[3679] - v[3680]
				+ v[10] * v[6213 + i2518]));
	v[3147] = v[2932] + v[281] * v[3685] - v[3064] * v[592];
	b3148 = b277;
	if (b3148) {
		v[3149] = -v[2922];
		v[3150] = -v[2920];
		v[3151] = -v[2916];
	}
	else {
		v[3149] = v[2922];
		v[3150] = v[2920];
		v[3151] = v[2916];
	};
	v[2884] = v[2884] + v[262] * v[2634] * v[2886] + v[263] * (v[2629] * v[2879] + v[2628] * v[2881] + v[2627] * v[2883]
		+ v[252] * v[3149] + v[251] * v[3150] + v[250] * v[3151]);
	v[3160] = v[2627] * v[2885] + v[264] * v[3149] + (v[2536] * v[2883] + v[252] * v[2884]) / v[602];
	v[3162] = v[2628] * v[2885] + v[264] * v[3150] + (v[2536] * v[2881] + v[251] * v[2884]) / v[602];
	v[3164] = v[2629] * v[2885] + v[264] * v[3151] + (v[2536] * v[2879] + v[250] * v[2884]) / v[602];
	v[3165] = -(v[2636] * v[2788]) - v[168] * v[3160];
	v[3166] = -(v[2636] * v[2787]) - v[167] * v[3160];
	v[3167] = -(v[2636] * v[2786]) - v[166] * v[3160];
	v[3168] = -(v[2637] * v[2788]) - v[168] * v[3162];
	v[3169] = -(v[2637] * v[2787]) - v[167] * v[3162];
	v[3170] = -(v[2637] * v[2786]) - v[166] * v[3162];
	v[3171] = -(v[2638] * v[2788]) - v[168] * v[3164];
	v[3172] = -(v[2638] * v[2787]) - v[167] * v[3164];
	v[3173] = -(v[2638] * v[2786]) - v[166] * v[3164];
	v[3174] = v[2937] + v[281] * (v[214] * v[3142] + v[212] * v[3143] + v[213] * v[3144]) + v[152] * v[3160]
		+ v[148] * v[3162] + v[144] * v[3164] + v[3115] * v[3361] + v[2595] * (-v[3624] + v[3628]) + v[214] * v[3651]
		+ v[212] * v[3686] + v[213] * v[3687] + v[2523] * (v[212] * v[2973] + v[2748] * v[582]) + v[2524] * (v[213] * v[2944]
			+ v[2749] * v[587]);
	v[3175] = v[236] * v[3121] + v[2602] * v[3688];
	v[3176] = v[230] * v[3121] + v[2602] * v[3689];
	v[3178] = v[1315] * v[3121] + v[2602] * (v[2850] / v[227] + v[2819] * v[3177]) + v[3175] * v[3690];
	v[3180] = v[1321] * v[3121] + v[2602] * (v[2845] / v[227] + v[2819] * v[3179]) + v[3176] * v[3691];
	v[3182] = (v[186] * v[3175] + v[184] * v[3176] + v[2602] * (v[2857] + v[3181] * v[3613]) + v[3121] * v[3692]) / v[227];
	v[3183] = v[238] * v[3122] + v[2603] * v[3693];
	v[3184] = v[230] * v[3122] + v[2603] * v[3689];
	v[3186] = v[1310] * v[3122] + v[2603] * (v[2855] / v[227] + v[2819] * v[3185]) + v[3183] * v[3694];
	v[3187] = v[3175] + v[3183];
	v[3188] = v[3178] + v[3186];
	v[3190] = (v[2644] * v[2804] + v[175] * v[3184] + v[177] * v[3187] + v[3316] * v[3613] + v[2603] * (v[2841]
		+ v[3189] * v[3613]) + v[1320] * v[3695]) / v[227];
	v[3192] = (v[182] * v[3183] + v[179] * v[3184] + v[2603] * (v[2847] + v[3191] * v[3613]) + v[1314] * v[3695]) / v[227];
	v[3193] = v[236] * v[3123] + v[2604] * v[3688];
	v[3194] = v[238] * v[3123] + v[2604] * v[3693];
	v[3195] = v[3184] + v[3193];
	v[3197] = (v[176] * v[3193] + v[177] * v[3194] + v[2604] * (v[2842] + v[3196] * v[3613]) + v[1318] * v[3696]) / v[227];
	v[3198] = v[3176] + v[3194];
	v[3200] = (v[2649] * v[2798] + v[181] * v[3193] + v[182] * v[3198] + v[3315] * v[3613] + v[2604] * (v[2846]
		+ v[3199] * v[3613]) + v[1325] * v[3696]) / v[227];
	v[3201] = v[3190] + v[3200];
	v[3203] = (v[2650] * v[2794] + v[187] * v[3194] + v[186] * v[3195] + v[3314] * v[3613] + v[2604] * (v[2853]
		+ v[3202] * v[3613]) + v[1328] * v[3696]) / v[227];
	v[3204] = v[3180] + v[3203];
	v[3205] = v[1437] * v[3183] + v[240] * v[3211] + v[3194] * v[3344] + v[238] * (v[2864] / v[227] + v[2819] * v[3863])
		+ v[3167] * v[94] + v[3166] * v[95] + v[3165] * v[96];
	v[3208] = (v[2650] * v[2835] + v[236] * v[2867] + v[240] * v[3175] + v[225] * v[3195]) / v[227] + v[237] * v[3215]
		+ v[2819] * v[3311] + v[3206] * v[3697] + v[3207] * v[3697] + v[3167] * v[91] + v[3166] * v[92] + v[3165] * v[93];
	v[3210] = v[2512] * v[3176] + v[225] * v[3218] + v[230] * (v[2819] * v[3698] + v[3719]) + v[3167] * v[88]
		+ v[3166] * v[89] + v[3165] * v[90];
	v[3212] = (v[2649] * v[2838] + v[238] * v[2875] + v[231] * v[3183] + v[226] * v[3198]) / v[227] + v[2819] * v[3312]
		+ v[3211] * v[3345] + (v[3213] + v[3214])*v[3699] + v[3170] * v[94] + v[3169] * v[95] + v[3168] * v[96];
	v[3216] = v[1439] * v[3175] + v[231] * v[3215] + v[3193] * v[3341] + v[236] * (v[2877] / v[227] + v[2819] * v[3864])
		+ v[3170] * v[91] + v[3169] * v[92] + v[3168] * v[93];
	v[3219] = v[2510] * v[3184] + v[226] * v[3218] + v[230] * (v[2819] * v[3700] + v[3718]) + v[3170] * v[88]
		+ v[3169] * v[89] + v[3168] * v[90];
	v[3220] = v[2819] * v[3313] + (v[2644] * v[2839] + v[238] * v[2871] + v[228] * v[3194] + v[3187] * v[3340]) / v[227]
		+ v[3211] * v[3343] + (v[3221] + v[3223])*v[3699] + v[3173] * v[94] + v[3172] * v[95] + v[3171] * v[96];
	v[3222] = v[2507] * v[3193] + v[3215] * v[3340] + v[236] * (v[2819] * v[3701] + v[3722]) + v[3173] * v[91]
		+ v[3172] * v[92] + v[3171] * v[93];
	v[3225] = v[1440] * v[3176] + v[1438] * v[3184] + v[228] * v[3218] + v[230] * (v[2873] / v[227] + v[2819] * v[3865])
		+ v[3173] * v[88] + v[3172] * v[89] + v[3171] * v[90];
	v[3226] = v[2620] * v[2786] + v[166] * v[3137];
	v[3227] = v[2620] * v[2787] + v[167] * v[3137];
	v[3228] = v[2620] * v[2788] + v[168] * v[3137];
	v[3229] = v[3226] * v[88] + v[3227] * v[89] + v[3228] * v[90];
	v[3230] = v[3226] * v[91] + v[3227] * v[92] + v[3228] * v[93];
	v[3231] = v[3226] * v[94] + v[3227] * v[95] + v[3228] * v[96];
	v[3232] = v[2617] * v[2786] + v[166] * v[3136];
	v[3233] = v[2617] * v[2787] + v[167] * v[3136];
	v[3234] = v[2617] * v[2788] + v[168] * v[3136];
	v[3235] = v[3232] * v[88] + v[3233] * v[89] + v[3234] * v[90];
	v[3236] = v[3232] * v[94] + v[3233] * v[95] + v[3234] * v[96];
	v[3237] = v[3232] * v[91] + v[3233] * v[92] + v[3234] * v[93];
	v[3238] = v[2623] * v[2786] + v[166] * v[3138];
	v[3239] = v[2623] * v[2787] + v[167] * v[3138];
	v[3240] = v[2623] * v[2788] + v[168] * v[3138];
	v[3241] = v[3238] * v[91] + v[3239] * v[92] + v[3240] * v[93];
	v[3242] = v[3238] * v[94] + v[3239] * v[95] + v[3240] * v[96];
	v[3243] = v[3238] * v[88] + v[3239] * v[89] + v[3240] * v[90];
	v[3244] = v[2624] * v[2786] + v[166] * v[3139];
	v[3245] = v[2624] * v[2787] + v[167] * v[3139];
	v[3246] = v[2624] * v[2788] + v[168] * v[3139];
	v[3247] = v[3244] * v[88] + v[3245] * v[89] + v[3246] * v[90];
	v[3248] = v[3244] * v[91] + v[3245] * v[92] + v[3246] * v[93];
	v[3249] = v[3244] * v[94] + v[3245] * v[95] + v[3246] * v[96];
	v[3250] = v[2618] * v[2786] + v[166] * v[3141];
	v[3251] = v[2618] * v[2787] + v[167] * v[3141];
	v[3252] = v[2618] * v[2788] + v[168] * v[3141];
	v[3253] = v[3250] * v[94] + v[3251] * v[95] + v[3252] * v[96];
	v[3254] = v[3250] * v[88] + v[3251] * v[89] + v[3252] * v[90];
	v[3255] = v[3250] * v[91] + v[3251] * v[92] + v[3252] * v[93];
	v[3256] = v[2619] * v[2786] + v[166] * v[3147];
	v[3257] = v[2619] * v[2787] + v[167] * v[3147];
	v[3258] = v[2619] * v[2788] + v[168] * v[3147];
	v[3259] = v[3256] * v[91] + v[3257] * v[92] + v[3258] * v[93];
	v[3260] = v[3256] * v[88] + v[3257] * v[89] + v[3258] * v[90];
	v[3261] = v[3256] * v[94] + v[3257] * v[95] + v[3258] * v[96];
	v[3262] = -(v[2646] * v[2776]) - v[2668] * v[2778] + 2e0*v[2645] * v[2779] + v[3241] + v[3247] + v[3253] + v[3259]
		- v[3201] * v[326] - v[3188] * v[329] + v[3192] * v[3454];
	v[3263] = v[2621] * v[2786] + v[166] * v[3140];
	v[3264] = v[2621] * v[2787] + v[167] * v[3140];
	v[3265] = v[2621] * v[2788] + v[168] * v[3140];
	v[3266] = v[3263] * v[88] + v[3264] * v[89] + v[3265] * v[90];
	v[3267] = v[3263] * v[94] + v[3264] * v[95] + v[3265] * v[96];
	v[3268] = v[3263] * v[91] + v[3264] * v[92] + v[3265] * v[93];
	v[3269] = v[2672] * v[2776] + 2e0*v[2648] * v[2778] - v[2668] * v[2779] - v[3201] * v[323] - v[3230] - v[3236] - v[3260]
		- v[3266] + v[3204] * v[329] + v[3197] * v[3453];
	v[3270] = -(v[2678] * v[2769]) + v[2674] * v[2771] - v[2679] * v[2772] - v[2661] * v[2777] + v[172] * v[3222]
		+ v[317] * v[3230] - v[318] * v[3237] - v[316] * v[3241];
	v[3271] = v[2625] * v[2786] + v[166] * v[3145];
	v[3272] = v[2625] * v[2787] + v[167] * v[3145];
	v[3273] = v[2625] * v[2788] + v[168] * v[3145];
	v[3274] = v[3271] * v[88] + v[3272] * v[89] + v[3273] * v[90];
	v[3275] = v[3271] * v[91] + v[3272] * v[92] + v[3273] * v[93];
	v[3276] = v[3271] * v[94] + v[3272] * v[95] + v[3273] * v[96];
	v[3277] = v[2622] * v[2786] + v[166] * v[3146];
	v[3278] = v[2622] * v[2787] + v[167] * v[3146];
	v[3279] = v[2622] * v[2788] + v[168] * v[3146];
	v[3280] = v[3277] * v[91] + v[3278] * v[92] + v[3279] * v[93];
	v[3281] = v[3277] * v[88] + v[3278] * v[89] + v[3279] * v[90];
	v[3282] = v[3277] * v[94] + v[3278] * v[95] + v[3279] * v[96];
	v[3283] = 2e0*v[2642] * v[2776] + v[2672] * v[2778] - v[2646] * v[2779] - v[3188] * v[323] - v[3242] + v[3204] * v[326]
		- v[3267] - v[3274] - v[3280] + v[3182] * v[3452];
	v[3284] = -(v[2677] * v[2769]) + v[2675] * v[2771] - v[2680] * v[2772] - v[2659] * v[2777] + v[172] * v[3220]
		+ v[317] * v[3231] - v[318] * v[3236] - v[316] * v[3242];
	v[3285] = -(v[2686] * v[2769]) + v[2692] * v[2771] - v[2682] * v[2772] - v[2658] * v[2777] + v[172] * v[3219]
		- v[316] * v[3247] - v[318] * v[3254] + v[317] * v[3266];
	v[3286] = v[3235] + v[3255];
	v[3287] = -(v[2685] * v[2769]) + v[2693] * v[2771] - v[2684] * v[2772] - v[2667] * v[2777] + v[172] * v[3212]
		- v[316] * v[3249] - v[318] * v[3253] + v[317] * v[3267];
	v[3288] = -(v[2689] * v[2769]) + v[2701] * v[2771] - v[2697] * v[2772] - v[2654] * v[2777] + v[172] * v[3210]
		- v[318] * v[3260] - v[316] * v[3274] + v[317] * v[3281];
	v[3289] = -(v[2688] * v[2769]) + v[2700] * v[2771] - v[2698] * v[2772] - v[2671] * v[2777] + v[172] * v[3208]
		- v[318] * v[3259] - v[316] * v[3275] + v[317] * v[3280];
	v[3290] = v[3248] + v[3276];
	v[3291] = v[3229] + v[3282];
	v[3292] = -(v[199] * v[3160]) - v[196] * v[3162] - v[193] * v[3164] + v[3138] * v[367] + v[3137] * v[368]
		+ v[3136] * v[369] + v[3139] * v[382] + v[3140] * v[383] + v[3141] * v[384] + v[3145] * v[397] + v[3146] * v[398]
		+ v[3147] * v[399] - v[2852] * v[90] - v[2851] * v[93] - v[2858] * v[96];
	v[3293] = -(v[198] * v[3160]) - v[195] * v[3162] - v[192] * v[3164] + v[3138] * v[364] + v[3137] * v[365]
		+ v[3136] * v[366] + v[3139] * v[379] + v[3140] * v[380] + v[3141] * v[381] + v[3145] * v[394] + v[3146] * v[395]
		+ v[3147] * v[396] - v[2852] * v[89] - v[2851] * v[92] - v[2858] * v[95];
	v[3294] = -(v[197] * v[3160]) - v[194] * v[3162] - v[191] * v[3164] + v[3138] * v[361] + v[3137] * v[362]
		+ v[3136] * v[363] + v[3139] * v[376] + v[3140] * v[377] + v[3141] * v[378] + v[3145] * v[391] + v[3146] * v[392]
		+ v[3147] * v[393] - v[2852] * v[88] - v[2851] * v[91] - v[2858] * v[94];
	v[3295] = v[3294] * v[76] + v[3293] * v[77] + v[3292] * v[78];
	v[3296] = (-v[3295] + v[3294] * v[82] + v[3293] * v[83] + v[3292] * v[84]) / 2e0;
	v[3297] = (-v[3295] + v[3294] * v[79] + v[3293] * v[80] + v[3292] * v[81]) / 2e0;
	v[3298] = (-(v[2699] * v[2741]) - v[2681] * v[2744] - v[2683] * v[2770] - 2e0*v[3178] + 2e0*v[3186] - v[319] * v[3243]
		- v[3248] * v[338] + v[3242] * v[3429] + v[3241] * v[3430] + v[3249] * v[3431] + v[3247] * v[3432] + v[3275] * v[3433]
		+ v[3274] * v[3434] - v[3276] * v[356] - v[2682] * v[3702] - v[2679] * v[3703] - v[2698] * v[3704] - v[2684] * v[3705]
		- v[2697] * v[3706] - v[2680] * v[3707]) / 2e0;
	v[3299] = (-(v[2676] * v[2769]) + v[2673] * v[2771] - v[2681] * v[2772] - v[2664] * v[2777] + v[172] * v[3225]
		+ v[317] * v[3229] - v[318] * v[3235] - v[316] * v[3243]) / 2e0;
	v[3301] = (v[2702] * v[2741] + v[2673] * v[2744] + v[2694] * v[2770] - 2e0*v[3180] + 2e0*v[3203] + v[319] * v[3229]
		+ v[3268] * v[338] - v[3231] * v[3429] - v[3230] * v[3430] - v[3267] * v[3431] - v[3266] * v[3432] - v[3280] * v[3433]
		- v[3281] * v[3434] + v[3282] * v[356] + v[2692] * v[3702] + v[2674] * v[3703] + v[2700] * v[3704] + v[2693] * v[3705]
		+ v[2701] * v[3706] + v[2675] * v[3707]) / 2e0;
	v[3302] = (-(v[2687] * v[2769]) + v[2694] * v[2771] - v[2683] * v[2772] - v[2655] * v[2777] + v[172] * v[3216]
		- v[316] * v[3248] - v[318] * v[3255] + v[317] * v[3268]) / 2e0;
	v[3303] = (-(v[2690] * v[2741]) - v[2676] * v[2744] - v[2687] * v[2770] - 2e0*v[3190] + 2e0*v[3200] - v[319] * v[3235]
		- v[3255] * v[338] + v[3236] * v[3429] + v[3237] * v[3430] + v[3253] * v[3431] + v[3254] * v[3432] + v[3259] * v[3433]
		+ v[3260] * v[3434] - v[3261] * v[356] - v[2686] * v[3702] - v[2678] * v[3703] - v[2688] * v[3704] - v[2685] * v[3705]
		- v[2689] * v[3706] - v[2677] * v[3707]) / 2e0;
	v[3324] = v[1498] * v[3301] - v[3302] + v[1495] * v[3303] + v[2742] * v[3710] * v[3867] - v[3330] * ((v[2651] * v[2741])
		/ 2e0 + (v[2664] * v[2744]) / 2e0 + v[2658] * v[2760] + v[2661] * v[2762] + v[2671] * v[2763] + v[2667] * v[2765]
		+ v[2654] * v[2766] + v[2659] * v[2768] + (v[2655] * v[2770]) / 2e0 + v[3231] - v[3237] + v[3222] * v[324] - v[3249]
		+ v[3254] + v[311] * v[3262] - v[312] * v[3269] + v[3275] - v[3281] - v[314] * v[3283] + v[3220] * v[330]
		- 2e0*v[2777] * v[3306] + v[3219] * v[334] + v[3205] * v[3426] + v[3216] * v[3427] + v[3225] * v[3428] + v[3212] * v[343]
		+ v[3210] * v[347] + v[3208] * v[351] + v[3286] * v[3562] + v[3290] * v[3563] + v[3291] * v[3564] + v[3443] *
		(v[2363] * v[2860] + v[2367] * v[2863] + v[2351] * v[2864] + v[2350] * v[2867] + v[2370] * v[2869] + v[2364] * v[2871]
			+ v[2355] * v[2873] + v[2358] * v[2875] + v[2357] * v[2877] + v[2486] * v[3175] + v[2487] * v[3176] + v[3182]
			+ v[2488] * v[3183] + v[2489] * v[3184] + v[1599] * v[3187] + v[3192] + v[2490] * v[3193] + v[2491] * v[3194]
			+ v[1593] * v[3195] + v[3197] + v[1597] * v[3198] + v[2857] * v[3206] + v[2855] * v[3207] + v[2853] * v[3209]
			+ v[2850] * v[3213] + v[2847] * v[3214] + v[2846] * v[3217] + v[2842] * v[3221] + v[2845] * v[3223] + v[2841] * v[3224]
			+ v[2794] * v[3311] + v[2798] * v[3312] + v[2804] * v[3313] + v[2835] * v[3314] + v[2838] * v[3315] + v[2839] * v[3316]
			+ v[3121] * v[3565] + v[3122] * v[3566] + v[3123] * v[3567] + v[2819] * v[3711] * v[3869]) + v[6333 + i2518]) + v[3708] * (
				-4e0*v[2715] * v[2742] + v[3298] * v[3709] + v[6438 + i2518]);
	v[3716] = (v[2690] * v[2769] - v[2702] * v[2771] + v[2699] * v[2772] + v[2651] * v[2777] - v[172] * v[3205]
		+ v[318] * v[3261] + v[316] * v[3276] - v[317] * v[3282]) / 2e0 + v[3324];
	v[3318] = v[3284] + v[3288];
	v[3319] = v[3287] + v[3289];
	v[3320] = v[3270] + v[3285];
	v[6529] = 0e0;
	v[6530] = 0e0;
	v[6531] = 0e0;
	v[6532] = 0e0;
	v[6533] = 0e0;
	v[6534] = 0e0;
	v[6535] = 0e0;
	v[6536] = 0e0;
	v[6537] = 0e0;
	v[6538] = 0e0;
	v[6539] = 0e0;
	v[6540] = 0e0;
	v[6541] = 0e0;
	v[6542] = v[2727];
	v[6543] = v[2725];
	v[6484] = 0e0;
	v[6485] = 0e0;
	v[6486] = 0e0;
	v[6487] = 0e0;
	v[6488] = 0e0;
	v[6489] = 0e0;
	v[6490] = 0e0;
	v[6491] = 0e0;
	v[6492] = 0e0;
	v[6493] = 0e0;
	v[6494] = 0e0;
	v[6495] = 0e0;
	v[6496] = v[2727];
	v[6497] = 0e0;
	v[6498] = v[2726];
	v[6454] = 0e0;
	v[6455] = 0e0;
	v[6456] = 0e0;
	v[6457] = 0e0;
	v[6458] = 0e0;
	v[6459] = 0e0;
	v[6460] = 0e0;
	v[6461] = 0e0;
	v[6462] = 0e0;
	v[6463] = 0e0;
	v[6464] = 0e0;
	v[6465] = 0e0;
	v[6466] = v[2725];
	v[6467] = v[2726];
	v[6468] = 0e0;
	v[6544] = v[2638] * v[2758] + v[128] * v[3164];
	v[6545] = v[2637] * v[2758] + v[128] * v[3162];
	v[6546] = v[2636] * v[2758] + v[128] * v[3160];
	v[6547] = v[2638] * v[2757] + v[130] * v[3164];
	v[6548] = v[2637] * v[2757] + v[130] * v[3162];
	v[6549] = v[2636] * v[2757] + v[130] * v[3160];
	v[6550] = -(v[2638] * v[2759]) + v[132] * v[3164];
	v[6551] = -(v[2637] * v[2759]) + v[132] * v[3162];
	v[6552] = -(v[2636] * v[2759]) + v[132] * v[3160];
	v[6553] = -v[3164];
	v[6554] = -v[3162];
	v[6555] = -v[3160];
	v[6556] = -v[3287] + v[3289] + v[314] * v[3318] + v[311] * v[3320] + 2e0*(v[2764] * v[3323] + v[3298] * v[3330]
		+ v[3290] * v[3332] + v[2742] * (-(v[2719] * v[3329]) + v[2713] * v[3715])) + v[310] * v[3716] + (v[2695] * v[2777]
			- v[172] * v[3269] + v[6528 + i2518]) / 2e0;
	v[6557] = v[3284] - v[3288] + v[314] * v[3319] + v[312] * v[3320] + 2e0*(v[2767] * v[3325] - v[3301] * v[3330]
		+ v[3291] * v[3332] + v[2742] * (v[2721] * v[3329] + v[2714] * v[3715])) + v[313] * (-v[3299] + v[3302] + v[3716]) + (-
		(v[2691] * v[2777]) + v[172] * v[3262] + v[6483 + i2518]) / 2e0;
	v[6558] = -v[3270] + v[3285] + v[312] * v[3318] + v[311] * v[3319] + v[315] * (-v[3299] + v[3324]) + 2e0*
		(v[2761] * v[3326] + v[3303] * v[3330] + v[3286] * v[3332] + v[2742] * (-(v[2723] * v[3329]) + v[2709] * v[3715])) +
		(v[2706] * v[2777] - v[172] * v[3283] + v[6453 + i2518]) / 2e0;
	v[3321] = (v[2935] + v[281] * (v[220] * v[3142] + v[218] * v[3143] + v[219] * v[3144]) + v[154] * v[3160]
		+ v[150] * v[3162] + v[146] * v[3164] - v[3174] + v[3117] * v[3361] + v[2597] * (-v[3624] + v[3628]) + v[220] * v[3651]
		+ v[218] * v[3686] + v[219] * v[3687] + v[2523] * (v[218] * v[2973] + v[2752] * v[582]) + v[2524] * (v[219] * v[2944]
			+ v[2753] * v[587])) / 2e0;
	v[3322] = (v[2936] + v[281] * (v[217] * v[3142] + v[215] * v[3143] + v[216] * v[3144]) + v[153] * v[3160]
		+ v[149] * v[3162] + v[145] * v[3164] - v[3174] + v[3116] * v[3361] + v[2596] * (-v[3624] + v[3628]) + v[217] * v[3651]
		+ v[215] * v[3686] + v[216] * v[3687] + v[2523] * (v[215] * v[2973] + v[2750] * v[582]) + v[2524] * (v[216] * v[2944]
			+ v[2751] * v[587])) / 2e0;
	Rc[i2518 - 1] += v[2729] * v[2735] + v[2728] * v[2736] + v[2717] * v[2745] + v[2718] * v[2746] + v[6089 + i2518]
		+ v[10] * v[6104 + i2518];
	for (i2733 = 1; i2733 <= 15; i2733++) {
		Kc[i2518 - 1][i2733 - 1] += v[3322] * v[4918 + i2733] + v[3321] * v[4933 + i2733] + v[3297] * v[4948 + i2733]
			+ v[3296] * v[4963 + i2733] + v[6543 + i2733] + v[10] * v[6558 + i2733];
	};/* end for */
};/* end for */
#pragma endregion


	delete[] v;
}

void FlexibleTriangularFace_RigidTriangularFace::HessianPhase1(Matrix& mHes)
{
	double Hes[4][4];
	v = DBG_NEW double[500];

#pragma region AceGen
	int i01; int i02;
	v[135] = (-u1A[0] - xi1A[0]) / 2e0;
	v[141] = (-u1A[1] - xi1A[1]) / 2e0;
	v[147] = (-u1A[2] - xi1A[2]) / 2e0;
	v[134] = v[135] + (u2A[0] + xi2A[0]) / 2e0;
	v[140] = v[141] + (u2A[1] + xi2A[1]) / 2e0;
	v[146] = v[147] + (u2A[2] + xi2A[2]) / 2e0;
	v[137] = v[135] + (u3A[0] + xi3A[0]) / 2e0;
	v[185] = 2e0*v[137];
	v[143] = v[141] + (u3A[1] + xi3A[1]) / 2e0;
	v[186] = 2e0*v[143];
	v[149] = v[147] + (u3A[2] + xi3A[2]) / 2e0;
	v[187] = 2e0*v[149];
	v[151] = -x1B[0] / 2e0;
	v[154] = -x1B[1] / 2e0;
	v[157] = -x1B[2] / 2e0;
	v[150] = v[151] + x2B[0] / 2e0;
	v[153] = v[154] + x2B[1] / 2e0;
	v[156] = v[157] + x2B[2] / 2e0;
	v[152] = v[151] + x3B[0] / 2e0;
	v[155] = v[154] + x3B[1] / 2e0;
	v[158] = v[157] + x3B[2] / 2e0;
	v[105] = Power(alphaB[0], 2);
	v[103] = (alphaB[0] * alphaB[1]) / 2e0;
	v[98] = Power(alphaB[1], 2);
	v[110] = (alphaB[1] * alphaB[2]) / 2e0;
	v[108] = (alphaB[0] * alphaB[2]) / 2e0;
	v[99] = Power(alphaB[2], 2);
	v[184] = v[98] + v[99];
	v[97] = 4e0 / (4e0 + v[105] + v[184]);
	v[100] = 1e0 - (v[184] * v[97]) / 2e0;
	v[101] = (-alphaB[2] + v[103])*v[97];
	v[102] = (alphaB[1] + v[108])*v[97];
	v[104] = (alphaB[2] + v[103])*v[97];
	v[106] = 1e0 - (v[97] * (v[105] + v[99])) / 2e0;
	v[107] = (-alphaB[0] + v[110])*v[97];
	v[109] = (-alphaB[1] + v[108])*v[97];
	v[111] = (alphaB[0] + v[110])*v[97];
	v[112] = 1e0 - (v[97] * (v[105] + v[98])) / 2e0;
	v[116] = QBi[0][0] * v[100] + QBi[1][0] * v[101] + QBi[2][0] * v[102];
	v[117] = QBi[0][1] * v[100] + QBi[1][1] * v[101] + QBi[2][1] * v[102];
	v[118] = QBi[0][2] * v[100] + QBi[1][2] * v[101] + QBi[2][2] * v[102];
	v[160] = v[116] * v[152] + v[117] * v[155] + v[118] * v[158];
	v[159] = v[116] * v[150] + v[117] * v[153] + v[118] * v[156];
	v[188] = 2e0*v[159];
	v[119] = QBi[0][0] * v[104] + QBi[1][0] * v[106] + QBi[2][0] * v[107];
	v[120] = QBi[0][1] * v[104] + QBi[1][1] * v[106] + QBi[2][1] * v[107];
	v[121] = QBi[0][2] * v[104] + QBi[1][2] * v[106] + QBi[2][2] * v[107];
	v[162] = v[119] * v[152] + v[120] * v[155] + v[121] * v[158];
	v[161] = v[119] * v[150] + v[120] * v[153] + v[121] * v[156];
	v[189] = 2e0*v[161];
	v[122] = QBi[0][0] * v[109] + QBi[1][0] * v[111] + QBi[2][0] * v[112];
	v[123] = QBi[0][1] * v[109] + QBi[1][1] * v[111] + QBi[2][1] * v[112];
	v[124] = QBi[0][2] * v[109] + QBi[1][2] * v[111] + QBi[2][2] * v[112];
	v[164] = v[122] * v[152] + v[123] * v[155] + v[124] * v[158];
	v[163] = v[122] * v[150] + v[123] * v[153] + v[124] * v[156];
	v[190] = 2e0*v[163];
	Hes[0][0] = 1e0*((v[134] * v[134]) + (v[140] * v[140]) + (v[146] * v[146]));
	Hes[0][1] = 0.5e0*(v[134] * v[185] + v[140] * v[186] + v[146] * v[187]);
	Hes[0][2] = 0.5e0*(-(v[134] * v[188]) - v[140] * v[189] - v[146] * v[190]);
	Hes[0][3] = -1e0*(v[134] * v[160] + v[140] * v[162] + v[146] * v[164]);
	Hes[1][1] = 1e0*((v[137] * v[137]) + (v[143] * v[143]) + (v[149] * v[149]));
	Hes[1][2] = -1e0*(v[137] * v[159] + v[143] * v[161] + v[149] * v[163]);
	Hes[1][3] = 0.5e0*(-(v[160] * v[185]) - v[162] * v[186] - v[164] * v[187]);
	Hes[2][2] = 1e0*((v[159] * v[159]) + (v[161] * v[161]) + (v[163] * v[163]));
	Hes[2][3] = 0.5e0*(v[160] * v[188] + v[162] * v[189] + v[164] * v[190]);
	Hes[3][3] = 1e0*((v[160] * v[160]) + (v[162] * v[162]) + (v[164] * v[164]));
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

void FlexibleTriangularFace_RigidTriangularFace::Report()
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
			if (cd->g_n[0] < *gnbb && eligible == true)
				db.myprintf("---------- Barrier activated! ---------- %.1f %c of contact layer is active.\n", 100.0*(1.0 - cd->g_n[0] / (*gnb)), 37);
				db.myprintf("\n---------I1: %d I2: %d \nfaceA: %d curveA: %d pointA: %d\nfaceB: %d curveB: %d pointB: %d\n", index1, index2, faceA->ID, deg_curveA, deg_pointA, faceB->ID, deg_curveB, deg_pointB);

		}
	}
}

void FlexibleTriangularFace_RigidTriangularFace::CompactReport()
{
	db.myprintf("Eligible %d\n", (int)eligible);
	db.myprintf("Deg point A: %d\n", deg_pointA);
	db.myprintf("Deg point B: %d\n", deg_pointB);
	db.myprintf("Deg curve A: %d\n", deg_curveA);
	db.myprintf("Deg curve B: %d\n", deg_curveB);
	db.myprintf("Gap %.6e\n", cd->g_n[0]);
}

void FlexibleTriangularFace_RigidTriangularFace::PredictorTimeStep(double kin)
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
		//Gamma A
		double za = cd->convective[0][0];
		double tha = cd->convective[0][1];
		double N1A = -(tha / 2.0) - za / 2.0;
		double N2A = 0.5 + za / 2.0;
		double N3A = 0.5 + tha / 2.0;
		
		Matrix VP_A(3);
		Matrix VP_B(3);
		Matrix vOB(3);
		Matrix dalphaB(3);
		Matrix P_B(3);

		//Vectors "(P-O)" - global system
		P_B = *GammaB - *ptrx0Bp;
		for (int i = 0; i < 3; i++)
		{
			vOB(i, 0) = db.nodes[node_B - 1]->copy_vel[i];
			dalphaB(i, 0) = db.nodes[node_B - 1]->copy_vel[i + 3];
		}
		//Obs: Xi é I3 nesse caso - acaba de salvar a configuração e o alpha avaliado em torno desse valor é nulo.
		Matrix omega_B = dalphaB;
		VP_B = vOB + cross(omega_B, P_B);

		for (int i = 0; i < 3; i++)
			VP_A(i, 0) = N1A * dui1A[i] + N2A * dui2A[i] + N3A * dui3A[i];
			
		double vrel = dot(VP_A - VP_B, *cd->n[0]);
		double deltaS = (cd->g_n[0] - (*gnb));
		if (vrel < 0)
			td->time_step_impact = (deltaS) / (-vrel);
		else
			td->time_step_impact = db.solution[db.current_solution_number - 1]->end_time;	//valor alto
	}
}
void FlexibleTriangularFace_RigidTriangularFace::AllocSpecificExplicit()
{
	//TODO-Explicit
}
void FlexibleTriangularFace_RigidTriangularFace::FreeSpecificExplicit()
{
	//TODO-Explicit
}
void FlexibleTriangularFace_RigidTriangularFace::MountLocalContributionsExplicit(double t)
{
	//TODO-Explicit
}
void FlexibleTriangularFace_RigidTriangularFace::SetVariablesExplicit(double t)
{
	//TODO-Explicit
}

void FlexibleTriangularFace_RigidTriangularFace::FinalUpdateExplicit(double t)
{
	//TODO-Explicit
}