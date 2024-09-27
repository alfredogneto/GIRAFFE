#include "FlexibleTriangularFace_FlexibleTriangularFace.h"

#include "SSContactData.h"
#include "ExecutionData.h"
#include "SuperNode.h"
#include "STLSurface.h"
#include "TriangularFace.h"
#include "Interface_1.h"
#include "VEMPolyhedron.h"
#include "Dynamic.h"
#include "Material.h"
#include "TimeStepControlData.h"
#include "Matrix.h"

#include "Database.h"
//Variaveis globais
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

FlexibleTriangularFace_FlexibleTriangularFace::FlexibleTriangularFace_FlexibleTriangularFace()
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
	super_node_B = 0;		//ID of super node B
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


FlexibleTriangularFace_FlexibleTriangularFace::~FlexibleTriangularFace_FlexibleTriangularFace()
{
	SetUnnactive();
	delete GammaA;
	delete GammaB;

	delete meq;
	
}

void FlexibleTriangularFace_FlexibleTriangularFace::AllocSpecific()
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
		xi1B = DBG_NEW double[3];
		xi2B = DBG_NEW double[3];
		xi3B = DBG_NEW double[3];
		u1A = DBG_NEW double[3];
		u2A = DBG_NEW double[3];
		u3A = DBG_NEW double[3];
		u1B = DBG_NEW double[3];
		u2B = DBG_NEW double[3];
		u3B = DBG_NEW double[3];
		dui1A = DBG_NEW double[3];
		dui2A = DBG_NEW double[3];
		dui3A = DBG_NEW double[3];
		dui1B = DBG_NEW double[3];
		dui2B = DBG_NEW double[3];
		dui3B = DBG_NEW double[3];
		ddui1A = DBG_NEW double[3];
		ddui2A = DBG_NEW double[3];
		ddui3A = DBG_NEW double[3];
		ddui1B = DBG_NEW double[3];
		ddui2B = DBG_NEW double[3];
		ddui3B = DBG_NEW double[3];

		fn = DBG_NEW double[3];
		ft = DBG_NEW double[3];
		
		Rc = DBG_NEW double[18];
		Kc = DBG_NEW double*[18];

		for (int i = 0; i < 18; i++)
			Kc[i] = DBG_NEW double[18];
		for (int ni = 0; ni < 18; ni++)
			for (int nj = 0; nj < 18; nj++)
				Kc[ni][nj] = 0.0;

		alloc_specific_control = true;
		vnrel = new double;
	}

}

void FlexibleTriangularFace_FlexibleTriangularFace::FreeSpecific()
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
		delete[] xi1B;
		delete[] xi2B;
		delete[] xi3B;
		delete[] u1A;
		delete[] u2A;
		delete[] u3A;
		delete[] u1B;
		delete[] u2B;
		delete[] u3B;
		delete[] dui1A;
		delete[] dui2A;
		delete[] dui3A;
		delete[] dui1B;
		delete[] dui2B;
		delete[] dui3B;
		delete[] ddui1A;
		delete[] ddui2A;
		delete[] ddui3A;
		delete[] ddui1B;
		delete[] ddui2B;
		delete[] ddui3B;

		delete[] fn;
		delete[] ft;
		
		delete[] Rc;
		for (int i = 0; i < 18; i++)
			delete[]Kc[i];
		delete[]Kc;

		alloc_specific_control = false;
		delete vnrel;
	}
}

bool FlexibleTriangularFace_FlexibleTriangularFace::HaveErrors()
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
					db.myprintf("\nFlexibleTriangular-FlexibleTriangular I1: %d I2: %d \nfaceA: %d curveA: %d pointA: %d\nfaceB: %d curveB: %d pointB: %d\n", index1, index2, faceA->ID, deg_curveA, deg_pointA, faceB->ID, deg_curveB, deg_pointB);
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
						db.myprintf("\nFlexibleTriangular-FlexibleTriangular I1: %d I2: %d \nfaceA: %d curveA: %d pointA: %d\nfaceB: %d curveB: %d pointB: %d\n", index1, index2, faceA->ID, deg_curveA, deg_pointA, faceB->ID, deg_curveB, deg_pointB);
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

void FlexibleTriangularFace_FlexibleTriangularFace::PreCalc()
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

void FlexibleTriangularFace_FlexibleTriangularFace::EvaluateNormalGap()
{
	//Calculo da função gap (escalar)
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


void FlexibleTriangularFace_FlexibleTriangularFace::SolveLCP()
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
		Matrix x1Bl;
		Matrix x2Bl;
		Matrix x3Bl;
		//Surface A
		Matrix xp;
		if (previous_evaluation == false)
		{
			x1Bl = *vertices_p_B[vertexIDsB[0] - 1];
			x2Bl = *vertices_p_B[vertexIDsB[1] - 1];
			x3Bl = *vertices_p_B[vertexIDsB[2] - 1];
			xp = *vertices_p_A[vertexIDsA[0] - 1];
		}
		else
		{
			x1Bl = *vertices_i_B[vertexIDsB[0] - 1];
			x2Bl = *vertices_i_B[vertexIDsB[1] - 1];
			x3Bl = *vertices_i_B[vertexIDsB[2] - 1];
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
		Matrix xp;

		if (previous_evaluation == false)
		{
			x1Al = *vertices_p_A[vertexIDsA[0] - 1];
			x2Al = *vertices_p_A[vertexIDsA[1] - 1];
			x3Al = *vertices_p_A[vertexIDsA[2] - 1];
			xp = *vertices_p_B[vertexIDsB[0] - 1];
		}
		else
		{
			x1Al = *vertices_i_A[vertexIDsA[0] - 1];
			x2Al = *vertices_i_A[vertexIDsA[1] - 1];
			x3Al = *vertices_i_A[vertexIDsA[2] - 1];
			xp = *vertices_i_B[vertexIDsB[0] - 1];
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
		Matrix x1Bl;
		Matrix x2Bl;

		if (previous_evaluation == false)
		{
			x1Al = *vertices_p_A[vertexIDsA[0] - 1];
			x2Al = *vertices_p_A[vertexIDsA[1] - 1];
			x1Bl = *vertices_p_B[vertexIDsB[0] - 1];
			x2Bl = *vertices_p_B[vertexIDsB[1] - 1];
		}
		else
		{
			x1Al = *vertices_i_A[vertexIDsA[0] - 1];
			x2Al = *vertices_i_A[vertexIDsA[1] - 1];
			x1Bl = *vertices_i_B[vertexIDsB[0] - 1];
			x2Bl = *vertices_i_B[vertexIDsB[1] - 1];
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
		Matrix x1Bl;
		Matrix x2Bl;

		if (previous_evaluation == false)
		{
			x1Al = *vertices_p_A[vertexIDsA[0] - 1];
			x1Bl = *vertices_p_B[vertexIDsB[0] - 1];
			x2Bl = *vertices_p_B[vertexIDsB[1] - 1];
		}
		else
		{
			x1Al = *vertices_i_A[vertexIDsA[0] - 1];
			x1Bl = *vertices_i_B[vertexIDsB[0] - 1];
			x2Bl = *vertices_i_B[vertexIDsB[1] - 1];
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
		Matrix x1Bl;
		//Surface A - edge
		Matrix x1Al;
		Matrix x2Al;

		if (previous_evaluation == false)
		{
			x1Bl = *vertices_p_B[vertexIDsB[0] - 1];
			x1Al = *vertices_p_A[vertexIDsA[0] - 1];
			x2Al = *vertices_p_A[vertexIDsA[1] - 1];
		}
		else
		{
			x1Bl = *vertices_i_B[vertexIDsB[0] - 1];
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
void FlexibleTriangularFace_FlexibleTriangularFace::SolvePreviousContact()
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

void FlexibleTriangularFace_FlexibleTriangularFace::SurfacePoints()
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
	if (previous_evaluation == false)
	{
		z = cd->convective[0][2];
		th = cd->convective[0][3];
		N1 = -(th / 2.0) - z / 2.0;
		N2 = 0.5 + z / 2.0;
		N3 = 0.5 + th / 2.0;
		*GammaB = N1 * (*vertices_p_B[vertexIDsB[0] - 1])
			+ N2 * (*vertices_p_B[vertexIDsB[1] - 1])
			+ N3 * (*vertices_p_B[vertexIDsB[2] - 1]);
	}
	else 
	{
		z = cd->copy_convective[0][2];
		th = cd->copy_convective[0][3];
		N1 = -(th / 2.0) - z / 2.0;
		N2 = 0.5 + z / 2.0;
		N3 = 0.5 + th / 2.0;
		*GammaB = N1 * (*vertices_i_B[vertexIDsB[0] - 1])
			+ N2 * (*vertices_i_B[vertexIDsB[1] - 1])
			+ N3 * (*vertices_i_B[vertexIDsB[2] - 1]);
	}
		
}

void FlexibleTriangularFace_FlexibleTriangularFace::SetVariables()
{
	for (int i = 0; i < 3; i++)
	{
		xi1A[i] = db.super_nodes[super_node_A - 1]->copy_coordinates[(vertexIDsA[0] - 1) * 3 + i];
		xi2A[i] = db.super_nodes[super_node_A - 1]->copy_coordinates[(vertexIDsA[1] - 1) * 3 + i];
		xi3A[i] = db.super_nodes[super_node_A - 1]->copy_coordinates[(vertexIDsA[2] - 1) * 3 + i];

		xi1B[i] = db.super_nodes[super_node_B - 1]->copy_coordinates[(vertexIDsB[0] - 1) * 3 + i];
		xi2B[i] = db.super_nodes[super_node_B - 1]->copy_coordinates[(vertexIDsB[1] - 1) * 3 + i];
		xi3B[i] = db.super_nodes[super_node_B - 1]->copy_coordinates[(vertexIDsB[2] - 1) * 3 + i];

		u1A[i] = db.super_nodes[super_node_A - 1]->displacements[(vertexIDsA[0] - 1) * 3 + i];
		u2A[i] = db.super_nodes[super_node_A - 1]->displacements[(vertexIDsA[1] - 1) * 3 + i];
		u3A[i] = db.super_nodes[super_node_A - 1]->displacements[(vertexIDsA[2] - 1) * 3 + i];

		u1B[i] = db.super_nodes[super_node_B - 1]->displacements[(vertexIDsB[0] - 1) * 3 + i];
		u2B[i] = db.super_nodes[super_node_B - 1]->displacements[(vertexIDsB[1] - 1) * 3 + i];
		u3B[i] = db.super_nodes[super_node_B - 1]->displacements[(vertexIDsB[2] - 1) * 3 + i];

		dui1A[i] = db.super_nodes[super_node_A - 1]->copy_vel[(vertexIDsA[0] - 1) * 3 + i];
		dui2A[i] = db.super_nodes[super_node_A - 1]->copy_vel[(vertexIDsA[1] - 1) * 3 + i];
		dui3A[i] = db.super_nodes[super_node_A - 1]->copy_vel[(vertexIDsA[2] - 1) * 3 + i];

		dui1B[i] = db.super_nodes[super_node_B - 1]->copy_vel[(vertexIDsB[0] - 1) * 3 + i];
		dui2B[i] = db.super_nodes[super_node_B - 1]->copy_vel[(vertexIDsB[1] - 1) * 3 + i];
		dui3B[i] = db.super_nodes[super_node_B - 1]->copy_vel[(vertexIDsB[2] - 1) * 3 + i];

		ddui1A[i] = db.super_nodes[super_node_A - 1]->copy_accel[(vertexIDsA[0] - 1) * 3 + i];
		ddui2A[i] = db.super_nodes[super_node_A - 1]->copy_accel[(vertexIDsA[1] - 1) * 3 + i];
		ddui3A[i] = db.super_nodes[super_node_A - 1]->copy_accel[(vertexIDsA[2] - 1) * 3 + i];

		ddui1B[i] = db.super_nodes[super_node_B - 1]->copy_accel[(vertexIDsB[0] - 1) * 3 + i];
		ddui2B[i] = db.super_nodes[super_node_B - 1]->copy_accel[(vertexIDsB[1] - 1) * 3 + i];
		ddui3B[i] = db.super_nodes[super_node_B - 1]->copy_accel[(vertexIDsB[2] - 1) * 3 + i];
	}
}

void FlexibleTriangularFace_FlexibleTriangularFace::MountLocalContributions()
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
	for (int i = 0; i < 18; i++)
	{
		Rc[i] = 0.0;
		for (int j = 0; j < 18; j++)
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
//	int i585, i603, i717, i878, i1319, i1360, i1576, i1577, i1578, i1579, i1580, i1581, i1591
//		, i1592, i1593, i1594, i1595, i1596, i1597, i1598, i1599, i1600, i1601, i1602, b411, b412
//		, b444, b478, b479, b480, b481, b487, b488, b497, b498, b514, b515, b516, b517, b533, b550
//		, b551, b563, b590, b591, b619, b620, b722, b723, b724, b725, b753, b757, b758, b770, b771
//		, b804, b832, b983, b1028, b1069, b1070, b1080, b1081, b1110, b1111, b1112, b1129, b1203
//		, b1207, b1211, b1220, b1221, b1248, b1267, b1328, b1329, b1333, b1338, b1339, b1407, b1408
//		, b1416, b1420, b1424, b1429, b1430;
//	v[1553] = (*a4)*(*ct);
//	v[1545] = (*epsn1)*(*n1);
//	v[475] = (*gnb) - (*gnbb);
//	v[490] = -1e0 + (*n1);
//	v[1082] = -1e0 + v[490];
//	v[495] = -1e0 + (*n2);
//	v[1087] = -1e0 + v[495];
//	v[141] = u1A[0] + xi1A[0];
//	v[230] = -v[141] / 2e0;
//	v[145] = u1A[1] + xi1A[1];
//	v[232] = -v[145] / 2e0;
//	v[149] = u1A[2] + xi1A[2];
//	v[234] = -v[149] / 2e0;
//	v[142] = u2A[0] + xi2A[0];
//	v[226] = v[142] / 2e0 + v[230];
//	v[146] = u2A[1] + xi2A[1];
//	v[227] = v[146] / 2e0 + v[232];
//	v[150] = u2A[2] + xi2A[2];
//	v[228] = v[150] / 2e0 + v[234];
//	v[143] = u3A[0] + xi3A[0];
//	v[231] = v[143] / 2e0 + v[230];
//	v[147] = u3A[1] + xi3A[1];
//	v[233] = v[147] / 2e0 + v[232];
//	v[151] = u3A[2] + xi3A[2];
//	v[235] = v[151] / 2e0 + v[234];
//	v[126] = -cAp[0] / 2e0;
//	v[128] = -cAp[1] / 2e0;
//	v[131] = -cAi[0] / 2e0;
//	v[133] = -cAi[1] / 2e0;
//	v[169] = u1B[0] + xi1B[0];
//	v[241] = -v[169] / 2e0;
//	v[173] = u1B[1] + xi1B[1];
//	v[243] = -v[173] / 2e0;
//	v[177] = u1B[2] + xi1B[2];
//	v[245] = -v[177] / 2e0;
//	v[170] = u2B[0] + xi2B[0];
//	v[237] = v[170] / 2e0 + v[241];
//	v[174] = u2B[1] + xi2B[1];
//	v[238] = v[174] / 2e0 + v[243];
//	v[178] = u2B[2] + xi2B[2];
//	v[239] = v[178] / 2e0 + v[245];
//	v[171] = u3B[0] + xi3B[0];
//	v[242] = v[171] / 2e0 + v[241];
//	v[175] = u3B[1] + xi3B[1];
//	v[244] = v[175] / 2e0 + v[243];
//	v[179] = u3B[2] + xi3B[2];
//	v[246] = v[179] / 2e0 + v[245];
//	v[154] = -cBp[0] / 2e0;
//	v[156] = -cBp[1] / 2e0;
//	v[159] = -cBi[0] / 2e0;
//	v[161] = -cBi[1] / 2e0;
//	v[125] = v[126] + v[128];
//	v[307] = -(v[125] * v[246]);
//	v[306] = -(v[125] * v[244]);
//	v[305] = -(v[125] * v[242]);
//	v[289] = -(v[125] * v[239]);
//	v[288] = -(v[125] * v[238]);
//	v[287] = -(v[125] * v[237]);
//	v[127] = 0.5e0 - v[126];
//	v[310] = -(v[127] * v[246]);
//	v[309] = -(v[127] * v[244]);
//	v[308] = -(v[127] * v[242]);
//	v[292] = -(v[127] * v[239]);
//	v[291] = -(v[127] * v[238]);
//	v[290] = -(v[127] * v[237]);
//	v[274] = v[127] * v[235];
//	v[273] = v[127] * v[233];
//	v[272] = v[127] * v[231];
//	v[129] = 0.5e0 - v[128];
//	v[313] = -(v[129] * v[246]);
//	v[312] = -(v[129] * v[244]);
//	v[311] = -(v[129] * v[242]);
//	v[295] = -(v[129] * v[239]);
//	v[294] = -(v[129] * v[238]);
//	v[293] = -(v[129] * v[237]);
//	v[259] = v[129] * v[228];
//	v[258] = v[129] * v[227];
//	v[257] = v[129] * v[226];
//	v[130] = v[131] + v[133];
//	v[132] = 0.5e0 - v[131];
//	v[134] = 0.5e0 - v[133];
//	v[153] = v[154] + v[156];
//	v[280] = -(v[153] * v[235]);
//	v[279] = -(v[153] * v[233]);
//	v[278] = -(v[153] * v[231]);
//	v[262] = -(v[153] * v[228]);
//	v[261] = -(v[153] * v[227]);
//	v[260] = -(v[153] * v[226]);
//	v[155] = 0.5e0 - v[154];
//	v[319] = v[155] * v[246];
//	v[318] = v[155] * v[244];
//	v[317] = v[155] * v[242];
//	v[283] = -(v[155] * v[235]);
//	v[282] = -(v[155] * v[233]);
//	v[281] = -(v[155] * v[231]);
//	v[265] = -(v[155] * v[228]);
//	v[264] = -(v[155] * v[227]);
//	v[263] = -(v[155] * v[226]);
//	v[157] = 0.5e0 - v[156];
//	v[2935] = 0e0;
//	v[2936] = v[125];
//	v[2937] = 0e0;
//	v[2938] = 0e0;
//	v[2939] = v[127];
//	v[2940] = 0e0;
//	v[2941] = 0e0;
//	v[2942] = v[129];
//	v[2943] = 0e0;
//	v[2944] = 0e0;
//	v[2945] = -v[153];
//	v[2946] = 0e0;
//	v[2947] = 0e0;
//	v[2948] = -v[155];
//	v[2949] = 0e0;
//	v[2950] = 0e0;
//	v[2951] = -v[157];
//	v[2952] = 0e0;
//	v[2917] = v[125];
//	v[2918] = 0e0;
//	v[2919] = 0e0;
//	v[2920] = v[127];
//	v[2921] = 0e0;
//	v[2922] = 0e0;
//	v[2923] = v[129];
//	v[2924] = 0e0;
//	v[2925] = 0e0;
//	v[2926] = -v[153];
//	v[2927] = 0e0;
//	v[2928] = 0e0;
//	v[2929] = -v[155];
//	v[2930] = 0e0;
//	v[2931] = 0e0;
//	v[2932] = -v[157];
//	v[2933] = 0e0;
//	v[2934] = 0e0;
//	v[2953] = 0e0;
//	v[2954] = 0e0;
//	v[2955] = v[125];
//	v[2956] = 0e0;
//	v[2957] = 0e0;
//	v[2958] = v[127];
//	v[2959] = 0e0;
//	v[2960] = 0e0;
//	v[2961] = v[129];
//	v[2962] = 0e0;
//	v[2963] = 0e0;
//	v[2964] = -v[153];
//	v[2965] = 0e0;
//	v[2966] = 0e0;
//	v[2967] = -v[155];
//	v[2968] = 0e0;
//	v[2969] = 0e0;
//	v[2970] = -v[157];
//	v[304] = v[157] * v[239];
//	v[303] = v[157] * v[238];
//	v[302] = v[157] * v[237];
//	v[286] = -(v[157] * v[235]);
//	v[285] = -(v[157] * v[233]);
//	v[284] = -(v[157] * v[231]);
//	v[268] = -(v[157] * v[228]);
//	v[267] = -(v[157] * v[227]);
//	v[266] = -(v[157] * v[226]);
//	v[158] = v[159] + v[161];
//	v[160] = 0.5e0 - v[159];
//	v[162] = 0.5e0 - v[161];
//	v[181] = (*a6)*ddui1A[0] + (*a5)*dui1A[0] + (*a4)*u1A[0];
//	v[1527] = v[125] * v[181];
//	v[182] = (*a6)*ddui1A[1] + (*a5)*dui1A[1] + (*a4)*u1A[1];
//	v[1533] = v[125] * v[182];
//	v[183] = (*a6)*ddui1A[2] + (*a5)*dui1A[2] + (*a4)*u1A[2];
//	v[1539] = v[125] * v[183];
//	v[184] = (*a6)*ddui2A[0] + (*a5)*dui2A[0] + (*a4)*u2A[0];
//	v[1526] = v[127] * v[184];
//	v[185] = (*a6)*ddui2A[1] + (*a5)*dui2A[1] + (*a4)*u2A[1];
//	v[1532] = v[127] * v[185];
//	v[186] = (*a6)*ddui2A[2] + (*a5)*dui2A[2] + (*a4)*u2A[2];
//	v[1538] = v[127] * v[186];
//	v[187] = (*a6)*ddui3A[0] + (*a5)*dui3A[0] + (*a4)*u3A[0];
//	v[1525] = v[129] * v[187];
//	v[188] = (*a6)*ddui3A[1] + (*a5)*dui3A[1] + (*a4)*u3A[1];
//	v[1531] = v[129] * v[188];
//	v[189] = (*a6)*ddui3A[2] + (*a5)*dui3A[2] + (*a4)*u3A[2];
//	v[1537] = v[129] * v[189];
//	v[190] = (*a6)*ddui1B[0] + (*a5)*dui1B[0] + (*a4)*u1B[0];
//	v[1524] = -(v[153] * v[190]);
//	v[191] = (*a6)*ddui1B[1] + (*a5)*dui1B[1] + (*a4)*u1B[1];
//	v[1530] = -(v[153] * v[191]);
//	v[192] = (*a6)*ddui1B[2] + (*a5)*dui1B[2] + (*a4)*u1B[2];
//	v[1536] = -(v[153] * v[192]);
//	v[193] = (*a6)*ddui2B[0] + (*a5)*dui2B[0] + (*a4)*u2B[0];
//	v[1523] = -(v[155] * v[193]);
//	v[194] = (*a6)*ddui2B[1] + (*a5)*dui2B[1] + (*a4)*u2B[1];
//	v[1529] = -(v[155] * v[194]);
//	v[195] = (*a6)*ddui2B[2] + (*a5)*dui2B[2] + (*a4)*u2B[2];
//	v[1535] = -(v[155] * v[195]);
//	v[196] = (*a6)*ddui3B[0] + (*a5)*dui3B[0] + (*a4)*u3B[0];
//	v[1522] = -(v[157] * v[196]);
//	v[1057] = v[1522] + v[1523] + v[1524] + v[1525] + v[1526] + v[1527];
//	v[197] = (*a6)*ddui3B[1] + (*a5)*dui3B[1] + (*a4)*u3B[1];
//	v[1528] = -(v[157] * v[197]);
//	v[1060] = v[1528] + v[1529] + v[1530] + v[1531] + v[1532] + v[1533];
//	v[198] = (*a6)*ddui3B[2] + (*a5)*dui3B[2] + (*a4)*u3B[2];
//	v[1534] = -(v[157] * v[198]);
//	v[945] = v[1534] + v[1535] + v[1536] + v[1537] + v[1538] + v[1539];
//	v[199] = v[125] * v[141] + v[127] * v[142] + v[129] * v[143] - v[153] * v[169] - v[155] * v[170] - v[157] * v[171];
//	v[251] = -v[199] / 2e0;
//	v[320] = v[157] * v[242] + v[251];
//	v[314] = v[153] * v[242] - v[251];
//	v[299] = v[155] * v[237] + v[251];
//	v[296] = v[153] * v[237] - v[251];
//	v[275] = v[129] * v[231] - v[251];
//	v[269] = v[125] * v[231] + v[251];
//	v[252] = v[127] * v[226] - v[251];
//	v[248] = v[125] * v[226] + v[251];
//	v[200] = v[125] * v[145] + v[127] * v[146] + v[129] * v[147] - v[153] * v[173] - v[155] * v[174] - v[157] * v[175];
//	v[253] = -v[200] / 2e0;
//	v[321] = v[157] * v[244] + v[253];
//	v[315] = v[153] * v[244] - v[253];
//	v[300] = v[155] * v[238] + v[253];
//	v[297] = v[153] * v[238] - v[253];
//	v[276] = v[129] * v[233] - v[253];
//	v[270] = v[125] * v[233] + v[253];
//	v[254] = v[127] * v[227] - v[253];
//	v[249] = v[125] * v[227] + v[253];
//	v[201] = v[125] * v[149] + v[127] * v[150] + v[129] * v[151] - v[153] * v[177] - v[155] * v[178] - v[157] * v[179];
//	v[477] = sqrt((v[199] * v[199]) + (v[200] * v[200]) + (v[201] * v[201]));
//	v[919] = 1e0 / Power(v[477], 2);
//	v[255] = -v[201] / 2e0;
//	v[322] = v[157] * v[246] + v[255];
//	v[316] = v[153] * v[246] - v[255];
//	v[301] = v[155] * v[239] + v[255];
//	v[298] = v[153] * v[239] - v[255];
//	v[277] = v[129] * v[235] - v[255];
//	v[271] = v[125] * v[235] + v[255];
//	v[256] = v[127] * v[228] - v[255];
//	v[250] = v[125] * v[228] + v[255];
//	v[202] = v[130] * xi1A[0] - v[158] * xi1B[0] + v[132] * xi2A[0] - v[160] * xi2B[0] + v[134] * xi3A[0] - v[162] * xi3B[0];
//	v[203] = v[130] * xi1A[1] - v[158] * xi1B[1] + v[132] * xi2A[1] - v[160] * xi2B[1] + v[134] * xi3A[1] - v[162] * xi3B[1];
//	v[204] = v[130] * xi1A[2] - v[158] * xi1B[2] + v[132] * xi2A[2] - v[160] * xi2B[2] + v[134] * xi3A[2] - v[162] * xi3B[2];
//	if (v[477] > 0.1e-7) { v01 = 1e0 / v[477]; v02 = (-(v01 / v[477])); v03 = (2e0*v01) / Power(v[477], 2); }
//	else {
//		v01 = (12500000e0 / 3e0)*(24e0 - (-0.1e-7 + v[477])*(0.24e10 - 2e0*(-1e0 + 100000000e0*v[477])*(0.2399999997e10
//			- 0.1199999994e18*v[477] - 0.3e17*(v[477] * v[477]))));
//		v02 = (-50000000e0 / 3e0)*(0.3599999994e10 - 0.4799999982e18*v[477] + 0.6e25*Power(v[477], 3)
//			+ 0.1799999982e26*(v[477] * v[477]));
//		v03 = 0.1e17*(799999997e0 - 0.599999994e17*v[477] - 0.3e17*(v[477] * v[477]));
//	};
//	v[211] = v03;
//	v[212] = v02;
//	v[213] = v01;
//	v[214] = v[199] * v[213];
//	v[451] = (v[214] * v[214]);
//	v[215] = v[200] * v[213];
//	v[3501] = v[125] * v[214];
//	v[3502] = v[125] * v[215];
//	v[3503] = 0e0;
//	v[3504] = v[127] * v[214];
//	v[3505] = v[127] * v[215];
//	v[3506] = 0e0;
//	v[3507] = v[129] * v[214];
//	v[3508] = v[129] * v[215];
//	v[3509] = 0e0;
//	v[3510] = -(v[153] * v[214]);
//	v[3511] = -(v[153] * v[215]);
//	v[3512] = 0e0;
//	v[3513] = -(v[155] * v[214]);
//	v[3514] = -(v[155] * v[215]);
//	v[3515] = 0e0;
//	v[3516] = -(v[157] * v[214]);
//	v[3517] = -(v[157] * v[215]);
//	v[3518] = 0e0;
//	v[1686] = v[193] * v[214] + v[194] * v[215];
//	v[1683] = v[196] * v[214] + v[197] * v[215];
//	v[1675] = v[190] * v[214] + v[191] * v[215];
//	v[1669] = v[187] * v[214] + v[188] * v[215];
//	v[1668] = v[184] * v[214] + v[185] * v[215];
//	v[1667] = v[181] * v[214] + v[182] * v[215];
//	v[1380] = -(v[214] * v[215]) / 2e0;
//	v[1066] = v[1522] * v[214] + v[1523] * v[214] + v[1524] * v[214] + v[1525] * v[214] + v[1526] * v[214] + v[1527] * v[214]
//		+ v[1528] * v[215] + v[1529] * v[215] + v[1530] * v[215] + v[1531] * v[215] + v[1532] * v[215] + v[1533] * v[215];
//	v[459] = (v[215] * v[215]);
//	v[458] = 2e0*v[1380] * v[157];
//	v[457] = 2e0*v[1380] * v[155];
//	v[456] = 2e0*v[1380] * v[153];
//	v[455] = -2e0*v[129] * v[1380];
//	v[454] = -2e0*v[127] * v[1380];
//	v[453] = -2e0*v[125] * v[1380];
//	v[3519] = v[453];
//	v[3520] = v[125] * v[459];
//	v[3521] = 0e0;
//	v[3522] = v[454];
//	v[3523] = v[127] * v[459];
//	v[3524] = 0e0;
//	v[3525] = v[455];
//	v[3526] = v[129] * v[459];
//	v[3527] = 0e0;
//	v[3528] = v[456];
//	v[3529] = -(v[153] * v[459]);
//	v[3530] = 0e0;
//	v[3531] = v[457];
//	v[3532] = -(v[155] * v[459]);
//	v[3533] = 0e0;
//	v[3534] = v[458];
//	v[3535] = -(v[157] * v[459]);
//	v[3536] = 0e0;
//	v[3537] = v[125] * v[451];
//	v[3538] = v[453];
//	v[3539] = 0e0;
//	v[3540] = v[127] * v[451];
//	v[3541] = v[454];
//	v[3542] = 0e0;
//	v[3543] = v[129] * v[451];
//	v[3544] = v[455];
//	v[3545] = 0e0;
//	v[3546] = -(v[153] * v[451]);
//	v[3547] = v[456];
//	v[3548] = 0e0;
//	v[3549] = -(v[155] * v[451]);
//	v[3550] = v[457];
//	v[3551] = 0e0;
//	v[3552] = -(v[157] * v[451]);
//	v[3553] = v[458];
//	v[3554] = 0e0;
//	v[216] = v[201] * v[213];
//	v[1661] = v[215] * v[216];
//	v[1660] = v[214] * v[216];
//	v[1585] = (*a4)*v[216];
//	v[1065] = v[1060] * v[216];
//	v[1064] = v[1057] * v[216];
//	v[467] = (v[216] * v[216]);
//	v[3483] = 0e0;
//	v[3484] = 0e0;
//	v[3485] = v[125] * v[467];
//	v[3486] = 0e0;
//	v[3487] = 0e0;
//	v[3488] = v[127] * v[467];
//	v[3489] = 0e0;
//	v[3490] = 0e0;
//	v[3491] = v[129] * v[467];
//	v[3492] = 0e0;
//	v[3493] = 0e0;
//	v[3494] = -(v[153] * v[467]);
//	v[3495] = 0e0;
//	v[3496] = 0e0;
//	v[3497] = -(v[155] * v[467]);
//	v[3498] = 0e0;
//	v[3499] = 0e0;
//	v[3500] = -(v[157] * v[467]);
//	v[794] = (v[1534] + v[1535] + v[1536] + v[1537] + v[1538] + v[1539])*v[216];
//	v[217] = sqrt((v[202] * v[202]) + (v[203] * v[203]) + (v[204] * v[204]));
//	if (v[217] > 0.1e-7) { v04 = 1e0 / v[217]; v05 = (-(v04 / v[217])); v06 = (2e0*v04) / Power(v[217], 2); }
//	else {
//		v04 = (12500000e0 / 3e0)*(24e0 - (-0.1e-7 + v[217])*(0.24e10 - 2e0*(-1e0 + 100000000e0*v[217])*(0.2399999997e10
//			- 0.1199999994e18*v[217] - 0.3e17*(v[217] * v[217]))));
//		v05 = (-50000000e0 / 3e0)*(0.3599999994e10 - 0.4799999982e18*v[217] + 0.6e25*Power(v[217], 3)
//			+ 0.1799999982e26*(v[217] * v[217]));
//		v06 = 0.1e17*(799999997e0 - 0.599999994e17*v[217] - 0.3e17*(v[217] * v[217]));
//	};
//	v[222] = v04;
//	v[223] = v[202] * v[222];
//	v[224] = v[203] * v[222];
//	v[225] = v[204] * v[222];
//	v[339] = -(invH[0][0] * v[248]) - invH[0][1] * v[269] - invH[0][2] * v[287] - invH[0][3] * v[305];
//	v[340] = -(invH[0][0] * v[249]) - invH[0][1] * v[270] - invH[0][2] * v[288] - invH[0][3] * v[306];
//	v[341] = -(invH[0][0] * v[250]) - invH[0][1] * v[271] - invH[0][2] * v[289] - invH[0][3] * v[307];
//	v[342] = -(invH[0][0] * v[252]) - invH[0][1] * v[272] - invH[0][2] * v[290] - invH[0][3] * v[308];
//	v[343] = -(invH[0][0] * v[254]) - invH[0][1] * v[273] - invH[0][2] * v[291] - invH[0][3] * v[309];
//	v[344] = -(invH[0][0] * v[256]) - invH[0][1] * v[274] - invH[0][2] * v[292] - invH[0][3] * v[310];
//	v[345] = -(invH[0][0] * v[257]) - invH[0][1] * v[275] - invH[0][2] * v[293] - invH[0][3] * v[311];
//	v[346] = -(invH[0][0] * v[258]) - invH[0][1] * v[276] - invH[0][2] * v[294] - invH[0][3] * v[312];
//	v[347] = -(invH[0][0] * v[259]) - invH[0][1] * v[277] - invH[0][2] * v[295] - invH[0][3] * v[313];
//	v[348] = -(invH[0][0] * v[260]) - invH[0][1] * v[278] - invH[0][2] * v[296] - invH[0][3] * v[314];
//	v[349] = -(invH[0][0] * v[261]) - invH[0][1] * v[279] - invH[0][2] * v[297] - invH[0][3] * v[315];
//	v[350] = -(invH[0][0] * v[262]) - invH[0][1] * v[280] - invH[0][2] * v[298] - invH[0][3] * v[316];
//	v[351] = -(invH[0][0] * v[263]) - invH[0][1] * v[281] - invH[0][2] * v[299] - invH[0][3] * v[317];
//	v[352] = -(invH[0][0] * v[264]) - invH[0][1] * v[282] - invH[0][2] * v[300] - invH[0][3] * v[318];
//	v[353] = -(invH[0][0] * v[265]) - invH[0][1] * v[283] - invH[0][2] * v[301] - invH[0][3] * v[319];
//	v[354] = -(invH[0][0] * v[266]) - invH[0][1] * v[284] - invH[0][2] * v[302] - invH[0][3] * v[320];
//	v[355] = -(invH[0][0] * v[267]) - invH[0][1] * v[285] - invH[0][2] * v[303] - invH[0][3] * v[321];
//	v[356] = -(invH[0][0] * v[268]) - invH[0][1] * v[286] - invH[0][2] * v[304] - invH[0][3] * v[322];
//	v[739] = -(v[181] * v[339]) - v[182] * v[340] - v[183] * v[341] - v[184] * v[342] - v[185] * v[343] - v[186] * v[344]
//		- v[187] * v[345] - v[188] * v[346] - v[189] * v[347] - v[190] * v[348] - v[191] * v[349] - v[192] * v[350] - v[193] * v[351]
//		- v[194] * v[352] - v[195] * v[353] - v[196] * v[354] - v[197] * v[355] - v[198] * v[356];
//	v[2841] = v[339];
//	v[2842] = v[340];
//	v[2843] = v[341];
//	v[2844] = v[342];
//	v[2845] = v[343];
//	v[2846] = v[344];
//	v[2847] = v[345];
//	v[2848] = v[346];
//	v[2849] = v[347];
//	v[2850] = v[348];
//	v[2851] = v[349];
//	v[2852] = v[350];
//	v[2853] = v[351];
//	v[2854] = v[352];
//	v[2855] = v[353];
//	v[2856] = v[354];
//	v[2857] = v[355];
//	v[2858] = v[356];
//	v[357] = -(invH[1][0] * v[248]) - invH[1][1] * v[269] - invH[1][2] * v[287] - invH[1][3] * v[305];
//	v[358] = -(invH[1][0] * v[249]) - invH[1][1] * v[270] - invH[1][2] * v[288] - invH[1][3] * v[306];
//	v[359] = -(invH[1][0] * v[250]) - invH[1][1] * v[271] - invH[1][2] * v[289] - invH[1][3] * v[307];
//	v[360] = -(invH[1][0] * v[252]) - invH[1][1] * v[272] - invH[1][2] * v[290] - invH[1][3] * v[308];
//	v[361] = -(invH[1][0] * v[254]) - invH[1][1] * v[273] - invH[1][2] * v[291] - invH[1][3] * v[309];
//	v[362] = -(invH[1][0] * v[256]) - invH[1][1] * v[274] - invH[1][2] * v[292] - invH[1][3] * v[310];
//	v[363] = -(invH[1][0] * v[257]) - invH[1][1] * v[275] - invH[1][2] * v[293] - invH[1][3] * v[311];
//	v[364] = -(invH[1][0] * v[258]) - invH[1][1] * v[276] - invH[1][2] * v[294] - invH[1][3] * v[312];
//	v[365] = -(invH[1][0] * v[259]) - invH[1][1] * v[277] - invH[1][2] * v[295] - invH[1][3] * v[313];
//	v[366] = -(invH[1][0] * v[260]) - invH[1][1] * v[278] - invH[1][2] * v[296] - invH[1][3] * v[314];
//	v[367] = -(invH[1][0] * v[261]) - invH[1][1] * v[279] - invH[1][2] * v[297] - invH[1][3] * v[315];
//	v[368] = -(invH[1][0] * v[262]) - invH[1][1] * v[280] - invH[1][2] * v[298] - invH[1][3] * v[316];
//	v[369] = -(invH[1][0] * v[263]) - invH[1][1] * v[281] - invH[1][2] * v[299] - invH[1][3] * v[317];
//	v[370] = -(invH[1][0] * v[264]) - invH[1][1] * v[282] - invH[1][2] * v[300] - invH[1][3] * v[318];
//	v[371] = -(invH[1][0] * v[265]) - invH[1][1] * v[283] - invH[1][2] * v[301] - invH[1][3] * v[319];
//	v[372] = -(invH[1][0] * v[266]) - invH[1][1] * v[284] - invH[1][2] * v[302] - invH[1][3] * v[320];
//	v[373] = -(invH[1][0] * v[267]) - invH[1][1] * v[285] - invH[1][2] * v[303] - invH[1][3] * v[321];
//	v[374] = -(invH[1][0] * v[268]) - invH[1][1] * v[286] - invH[1][2] * v[304] - invH[1][3] * v[322];
//	v[737] = -(v[181] * v[357]) - v[182] * v[358] - v[183] * v[359] - v[184] * v[360] - v[185] * v[361] - v[186] * v[362]
//		- v[187] * v[363] - v[188] * v[364] - v[189] * v[365] - v[190] * v[366] - v[191] * v[367] - v[192] * v[368] - v[193] * v[369]
//		- v[194] * v[370] - v[195] * v[371] - v[196] * v[372] - v[197] * v[373] - v[198] * v[374];
//	v[2823] = v[357];
//	v[2824] = v[358];
//	v[2825] = v[359];
//	v[2826] = v[360];
//	v[2827] = v[361];
//	v[2828] = v[362];
//	v[2829] = v[363];
//	v[2830] = v[364];
//	v[2831] = v[365];
//	v[2832] = v[366];
//	v[2833] = v[367];
//	v[2834] = v[368];
//	v[2835] = v[369];
//	v[2836] = v[370];
//	v[2837] = v[371];
//	v[2838] = v[372];
//	v[2839] = v[373];
//	v[2840] = v[374];
//	v[375] = -(invH[2][0] * v[248]) - invH[2][1] * v[269] - invH[2][2] * v[287] - invH[2][3] * v[305];
//	v[376] = -(invH[2][0] * v[249]) - invH[2][1] * v[270] - invH[2][2] * v[288] - invH[2][3] * v[306];
//	v[377] = -(invH[2][0] * v[250]) - invH[2][1] * v[271] - invH[2][2] * v[289] - invH[2][3] * v[307];
//	v[378] = -(invH[2][0] * v[252]) - invH[2][1] * v[272] - invH[2][2] * v[290] - invH[2][3] * v[308];
//	v[379] = -(invH[2][0] * v[254]) - invH[2][1] * v[273] - invH[2][2] * v[291] - invH[2][3] * v[309];
//	v[380] = -(invH[2][0] * v[256]) - invH[2][1] * v[274] - invH[2][2] * v[292] - invH[2][3] * v[310];
//	v[381] = -(invH[2][0] * v[257]) - invH[2][1] * v[275] - invH[2][2] * v[293] - invH[2][3] * v[311];
//	v[382] = -(invH[2][0] * v[258]) - invH[2][1] * v[276] - invH[2][2] * v[294] - invH[2][3] * v[312];
//	v[383] = -(invH[2][0] * v[259]) - invH[2][1] * v[277] - invH[2][2] * v[295] - invH[2][3] * v[313];
//	v[384] = -(invH[2][0] * v[260]) - invH[2][1] * v[278] - invH[2][2] * v[296] - invH[2][3] * v[314];
//	v[385] = -(invH[2][0] * v[261]) - invH[2][1] * v[279] - invH[2][2] * v[297] - invH[2][3] * v[315];
//	v[386] = -(invH[2][0] * v[262]) - invH[2][1] * v[280] - invH[2][2] * v[298] - invH[2][3] * v[316];
//	v[387] = -(invH[2][0] * v[263]) - invH[2][1] * v[281] - invH[2][2] * v[299] - invH[2][3] * v[317];
//	v[388] = -(invH[2][0] * v[264]) - invH[2][1] * v[282] - invH[2][2] * v[300] - invH[2][3] * v[318];
//	v[389] = -(invH[2][0] * v[265]) - invH[2][1] * v[283] - invH[2][2] * v[301] - invH[2][3] * v[319];
//	v[390] = -(invH[2][0] * v[266]) - invH[2][1] * v[284] - invH[2][2] * v[302] - invH[2][3] * v[320];
//	v[391] = -(invH[2][0] * v[267]) - invH[2][1] * v[285] - invH[2][2] * v[303] - invH[2][3] * v[321];
//	v[392] = -(invH[2][0] * v[268]) - invH[2][1] * v[286] - invH[2][2] * v[304] - invH[2][3] * v[322];
//	v[735] = v[181] * v[375] + v[182] * v[376] + v[183] * v[377] + v[184] * v[378] + v[185] * v[379] + v[186] * v[380]
//		+ v[187] * v[381] + v[188] * v[382] + v[189] * v[383] + v[190] * v[384] + v[191] * v[385] + v[192] * v[386] + v[193] * v[387]
//		+ v[194] * v[388] + v[195] * v[389] + v[196] * v[390] + v[197] * v[391] + v[198] * v[392];
//	v[2859] = v[375];
//	v[2860] = v[376];
//	v[2861] = v[377];
//	v[2862] = v[378];
//	v[2863] = v[379];
//	v[2864] = v[380];
//	v[2865] = v[381];
//	v[2866] = v[382];
//	v[2867] = v[383];
//	v[2868] = v[384];
//	v[2869] = v[385];
//	v[2870] = v[386];
//	v[2871] = v[387];
//	v[2872] = v[388];
//	v[2873] = v[389];
//	v[2874] = v[390];
//	v[2875] = v[391];
//	v[2876] = v[392];
//	v[393] = -(invH[3][0] * v[248]) - invH[3][1] * v[269] - invH[3][2] * v[287] - invH[3][3] * v[305];
//	v[646] = v[228] * v[339] + v[235] * v[357] - v[239] * v[375] - v[246] * v[393];
//	v[645] = v[227] * v[339] + v[233] * v[357] - v[238] * v[375] - v[244] * v[393];
//	v[644] = v[226] * v[339] + v[231] * v[357] - v[237] * v[375] - v[242] * v[393];
//	v[394] = -(invH[3][0] * v[249]) - invH[3][1] * v[270] - invH[3][2] * v[288] - invH[3][3] * v[306];
//	v[650] = v[228] * v[340] + v[235] * v[358] - v[239] * v[376] - v[246] * v[394];
//	v[649] = v[227] * v[340] + v[233] * v[358] - v[238] * v[376] - v[244] * v[394];
//	v[648] = v[226] * v[340] + v[231] * v[358] - v[237] * v[376] - v[242] * v[394];
//	v[395] = -(invH[3][0] * v[250]) - invH[3][1] * v[271] - invH[3][2] * v[289] - invH[3][3] * v[307];
//	v[654] = v[228] * v[341] + v[235] * v[359] - v[239] * v[377] - v[246] * v[395];
//	v[653] = v[227] * v[341] + v[233] * v[359] - v[238] * v[377] - v[244] * v[395];
//	v[652] = v[226] * v[341] + v[231] * v[359] - v[237] * v[377] - v[242] * v[395];
//	v[396] = -(invH[3][0] * v[252]) - invH[3][1] * v[272] - invH[3][2] * v[290] - invH[3][3] * v[308];
//	v[658] = v[228] * v[342] + v[235] * v[360] - v[239] * v[378] - v[246] * v[396];
//	v[657] = v[227] * v[342] + v[233] * v[360] - v[238] * v[378] - v[244] * v[396];
//	v[656] = v[226] * v[342] + v[231] * v[360] - v[237] * v[378] - v[242] * v[396];
//	v[397] = -(invH[3][0] * v[254]) - invH[3][1] * v[273] - invH[3][2] * v[291] - invH[3][3] * v[309];
//	v[662] = v[228] * v[343] + v[235] * v[361] - v[239] * v[379] - v[246] * v[397];
//	v[661] = v[227] * v[343] + v[233] * v[361] - v[238] * v[379] - v[244] * v[397];
//	v[660] = v[226] * v[343] + v[231] * v[361] - v[237] * v[379] - v[242] * v[397];
//	v[398] = -(invH[3][0] * v[256]) - invH[3][1] * v[274] - invH[3][2] * v[292] - invH[3][3] * v[310];
//	v[666] = v[228] * v[344] + v[235] * v[362] - v[239] * v[380] - v[246] * v[398];
//	v[665] = v[227] * v[344] + v[233] * v[362] - v[238] * v[380] - v[244] * v[398];
//	v[664] = v[226] * v[344] + v[231] * v[362] - v[237] * v[380] - v[242] * v[398];
//	v[399] = -(invH[3][0] * v[257]) - invH[3][1] * v[275] - invH[3][2] * v[293] - invH[3][3] * v[311];
//	v[670] = v[228] * v[345] + v[235] * v[363] - v[239] * v[381] - v[246] * v[399];
//	v[669] = v[227] * v[345] + v[233] * v[363] - v[238] * v[381] - v[244] * v[399];
//	v[668] = v[226] * v[345] + v[231] * v[363] - v[237] * v[381] - v[242] * v[399];
//	v[400] = -(invH[3][0] * v[258]) - invH[3][1] * v[276] - invH[3][2] * v[294] - invH[3][3] * v[312];
//	v[674] = v[228] * v[346] + v[235] * v[364] - v[239] * v[382] - v[246] * v[400];
//	v[673] = v[227] * v[346] + v[233] * v[364] - v[238] * v[382] - v[244] * v[400];
//	v[672] = v[226] * v[346] + v[231] * v[364] - v[237] * v[382] - v[242] * v[400];
//	v[401] = -(invH[3][0] * v[259]) - invH[3][1] * v[277] - invH[3][2] * v[295] - invH[3][3] * v[313];
//	v[678] = v[228] * v[347] + v[235] * v[365] - v[239] * v[383] - v[246] * v[401];
//	v[677] = v[227] * v[347] + v[233] * v[365] - v[238] * v[383] - v[244] * v[401];
//	v[676] = v[226] * v[347] + v[231] * v[365] - v[237] * v[383] - v[242] * v[401];
//	v[402] = -(invH[3][0] * v[260]) - invH[3][1] * v[278] - invH[3][2] * v[296] - invH[3][3] * v[314];
//	v[682] = v[228] * v[348] + v[235] * v[366] - v[239] * v[384] - v[246] * v[402];
//	v[681] = v[227] * v[348] + v[233] * v[366] - v[238] * v[384] - v[244] * v[402];
//	v[680] = v[226] * v[348] + v[231] * v[366] - v[237] * v[384] - v[242] * v[402];
//	v[403] = -(invH[3][0] * v[261]) - invH[3][1] * v[279] - invH[3][2] * v[297] - invH[3][3] * v[315];
//	v[686] = v[228] * v[349] + v[235] * v[367] - v[239] * v[385] - v[246] * v[403];
//	v[685] = v[227] * v[349] + v[233] * v[367] - v[238] * v[385] - v[244] * v[403];
//	v[684] = v[226] * v[349] + v[231] * v[367] - v[237] * v[385] - v[242] * v[403];
//	v[404] = -(invH[3][0] * v[262]) - invH[3][1] * v[280] - invH[3][2] * v[298] - invH[3][3] * v[316];
//	v[690] = v[228] * v[350] + v[235] * v[368] - v[239] * v[386] - v[246] * v[404];
//	v[689] = v[227] * v[350] + v[233] * v[368] - v[238] * v[386] - v[244] * v[404];
//	v[688] = v[226] * v[350] + v[231] * v[368] - v[237] * v[386] - v[242] * v[404];
//	v[405] = -(invH[3][0] * v[263]) - invH[3][1] * v[281] - invH[3][2] * v[299] - invH[3][3] * v[317];
//	v[694] = v[228] * v[351] + v[235] * v[369] - v[239] * v[387] - v[246] * v[405];
//	v[693] = v[227] * v[351] + v[233] * v[369] - v[238] * v[387] - v[244] * v[405];
//	v[692] = v[226] * v[351] + v[231] * v[369] - v[237] * v[387] - v[242] * v[405];
//	v[406] = -(invH[3][0] * v[264]) - invH[3][1] * v[282] - invH[3][2] * v[300] - invH[3][3] * v[318];
//	v[698] = v[228] * v[352] + v[235] * v[370] - v[239] * v[388] - v[246] * v[406];
//	v[697] = v[227] * v[352] + v[233] * v[370] - v[238] * v[388] - v[244] * v[406];
//	v[696] = v[226] * v[352] + v[231] * v[370] - v[237] * v[388] - v[242] * v[406];
//	v[407] = -(invH[3][0] * v[265]) - invH[3][1] * v[283] - invH[3][2] * v[301] - invH[3][3] * v[319];
//	v[702] = v[228] * v[353] + v[235] * v[371] - v[239] * v[389] - v[246] * v[407];
//	v[701] = v[227] * v[353] + v[233] * v[371] - v[238] * v[389] - v[244] * v[407];
//	v[700] = v[226] * v[353] + v[231] * v[371] - v[237] * v[389] - v[242] * v[407];
//	v[408] = -(invH[3][0] * v[266]) - invH[3][1] * v[284] - invH[3][2] * v[302] - invH[3][3] * v[320];
//	v[706] = v[228] * v[354] + v[235] * v[372] - v[239] * v[390] - v[246] * v[408];
//	v[705] = v[227] * v[354] + v[233] * v[372] - v[238] * v[390] - v[244] * v[408];
//	v[704] = v[226] * v[354] + v[231] * v[372] - v[237] * v[390] - v[242] * v[408];
//	v[409] = -(invH[3][0] * v[267]) - invH[3][1] * v[285] - invH[3][2] * v[303] - invH[3][3] * v[321];
//	v[710] = v[228] * v[355] + v[235] * v[373] - v[239] * v[391] - v[246] * v[409];
//	v[709] = v[227] * v[355] + v[233] * v[373] - v[238] * v[391] - v[244] * v[409];
//	v[708] = v[226] * v[355] + v[231] * v[373] - v[237] * v[391] - v[242] * v[409];
//	v[410] = -(invH[3][0] * v[268]) - invH[3][1] * v[286] - invH[3][2] * v[304] - invH[3][3] * v[322];
//	v[733] = v[181] * v[393] + v[182] * v[394] + v[183] * v[395] + v[184] * v[396] + v[185] * v[397] + v[186] * v[398]
//		+ v[187] * v[399] + v[188] * v[400] + v[189] * v[401] + v[190] * v[402] + v[191] * v[403] + v[192] * v[404] + v[193] * v[405]
//		+ v[194] * v[406] + v[195] * v[407] + v[196] * v[408] + v[197] * v[409] + v[198] * v[410];
//	v[714] = v[228] * v[356] + v[235] * v[374] - v[239] * v[392] - v[246] * v[410];
//	v[4005] = v[646];
//	v[4006] = v[650];
//	v[4007] = v[654];
//	v[4008] = v[658];
//	v[4009] = v[662];
//	v[4010] = v[666];
//	v[4011] = v[670];
//	v[4012] = v[674];
//	v[4013] = v[678];
//	v[4014] = v[682];
//	v[4015] = v[686];
//	v[4016] = v[690];
//	v[4017] = v[694];
//	v[4018] = v[698];
//	v[4019] = v[702];
//	v[4020] = v[706];
//	v[4021] = v[710];
//	v[4022] = v[714];
//	v[713] = v[227] * v[356] + v[233] * v[374] - v[238] * v[392] - v[244] * v[410];
//	v[4023] = v[645];
//	v[4024] = v[649];
//	v[4025] = v[653];
//	v[4026] = v[657];
//	v[4027] = v[661];
//	v[4028] = v[665];
//	v[4029] = v[669];
//	v[4030] = v[673];
//	v[4031] = v[677];
//	v[4032] = v[681];
//	v[4033] = v[685];
//	v[4034] = v[689];
//	v[4035] = v[693];
//	v[4036] = v[697];
//	v[4037] = v[701];
//	v[4038] = v[705];
//	v[4039] = v[709];
//	v[4040] = v[713];
//	v[712] = v[226] * v[356] + v[231] * v[374] - v[237] * v[392] - v[242] * v[410];
//	v[4041] = v[644];
//	v[4042] = v[648];
//	v[4043] = v[652];
//	v[4044] = v[656];
//	v[4045] = v[660];
//	v[4046] = v[664];
//	v[4047] = v[668];
//	v[4048] = v[672];
//	v[4049] = v[676];
//	v[4050] = v[680];
//	v[4051] = v[684];
//	v[4052] = v[688];
//	v[4053] = v[692];
//	v[4054] = v[696];
//	v[4055] = v[700];
//	v[4056] = v[704];
//	v[4057] = v[708];
//	v[4058] = v[712];
//	v[2877] = v[393];
//	v[2878] = v[394];
//	v[2879] = v[395];
//	v[2880] = v[396];
//	v[2881] = v[397];
//	v[2882] = v[398];
//	v[2883] = v[399];
//	v[2884] = v[400];
//	v[2885] = v[401];
//	v[2886] = v[402];
//	v[2887] = v[403];
//	v[2888] = v[404];
//	v[2889] = v[405];
//	v[2890] = v[406];
//	v[2891] = v[407];
//	v[2892] = v[408];
//	v[2893] = v[409];
//	v[2894] = v[410];
//	b411 = sqrt(Power(v[215] * v[223] - v[214] * v[224], 2) + Power(-(v[216] * v[223]) + v[214] * v[225], 2) + Power
//	(v[216] * v[224] - v[215] * v[225], 2)) > 0.1e-7;
//	if (b411) {
//		v[413] = v[216] * v[224] - v[215] * v[225];
//		v[414] = -(v[216] * v[223]) + v[214] * v[225];
//		v[415] = v[215] * v[223] - v[214] * v[224];
//		v[416] = sqrt((v[413] * v[413]) + (v[414] * v[414]) + (v[415] * v[415]));
//		v[990] = 1e0 / Power(v[416], 2);
//		v[861] = v[416];
//		v[1001] = 1e0 - (v[861] * v[861]);
//		v[1756] = 1e0 / Power(v[1001], 0.15e1);
//		v[996] = 1e0 / sqrt(v[1001]);
//		v[860] = asin(v[861]) / 2e0;
//		v[1755] = tan(v[860]);
//		v[995] = 1e0 / Power(cos(v[860]), 2);
//		v[1586] = v[995] * v[996];
//		v[418] = 2e0*tan(v[860]);
//		if (v[416] > 0.1e-7) { v07 = 1e0 / v[416]; v08 = (-(v07 / v[416])); v09 = (2e0*v07) / Power(v[416], 2); }
//		else {
//			v07 = (12500000e0 / 3e0)*(24e0 - (-0.1e-7 + v[416])*(0.24e10 - 2e0*(-1e0 + 100000000e0*v[416])*
//				(0.2399999997e10 - 0.1199999994e18*v[416] - 0.3e17*(v[416] * v[416]))));
//			v08 = (-50000000e0 / 3e0)*(0.3599999994e10 - 0.4799999982e18*v[416] + 0.6e25*Power(v[416], 3)
//				+ 0.1799999982e26*(v[416] * v[416]));
//			v09 = 0.1e17*(799999997e0 - 0.599999994e17*v[416] - 0.3e17*(v[416] * v[416]));
//		};
//		v[422] = v09;
//		v[423] = v08;
//		v[424] = v07;
//		v[1754] = v[418] * v[423] + v[1586] * v[424];
//		v[1540] = v[418] * v[424];
//		v[425] = v[1540] * v[413];
//		v[436] = (v[425] * v[425]);
//		v[426] = v[1540] * v[414];
//		v[434] = (v[425] * v[426]) / 2e0;
//		v[429] = (v[426] * v[426]);
//		v[844] = -v[429] - v[436];
//		v[427] = v[1540] * v[415];
//		v[839] = v[427] + v[434];
//		v[837] = -v[427] + v[434];
//		v[441] = (v[426] * v[427]) / 2e0;
//		v[843] = v[425] + v[441];
//		v[841] = -v[425] + v[441];
//		v[439] = (v[425] * v[427]) / 2e0;
//		v[842] = -v[426] + v[439];
//		v[838] = v[426] + v[439];
//		v[430] = (v[427] * v[427]);
//		v[848] = 4e0 + v[429] + v[430] + v[436];
//		v[1757] = 1e0 / Power(v[848], 3);
//		v[1268] = 1e0 / Power(v[848], 2);
//		v[840] = -v[430] - v[436];
//		v[836] = -v[429] - v[430];
//		v[428] = 4e0 / v[848];
//		v[431] = 1e0 + (v[428] * v[836]) / 2e0;
//		v[432] = v[428] * v[837];
//		v[433] = v[428] * v[838];
//		v[435] = v[428] * v[839];
//		v[437] = 1e0 + (v[428] * v[840]) / 2e0;
//		v[438] = v[428] * v[841];
//		v[440] = v[428] * v[842];
//		v[442] = v[428] * v[843];
//		v[443] = 1e0 + (v[428] * v[844]) / 2e0;
//	}
//	else {
//		v[431] = 1e0;
//		v[432] = 0e0;
//		v[433] = 0e0;
//		v[435] = 0e0;
//		v[437] = 1e0;
//		v[438] = 0e0;
//		v[440] = 0e0;
//		v[442] = 0e0;
//		v[443] = 1e0;
//	};
//	if ((*previouscontact)) {
//		v[809] = 1e0 - v[467];
//		v[807] = 1e0 - v[459];
//		v[805] = 1e0 - v[451];
//		v[448] = v[130] * v[149] + v[132] * v[150] + v[134] * v[151] - v[158] * v[177] - v[160] * v[178] - v[162] * v[179]
//			+ gti[0] * v[440] + gti[1] * v[442] + gti[2] * v[443];
//		v[1541] = v[216] * v[448];
//		v[447] = v[130] * v[145] + v[132] * v[146] + v[134] * v[147] - v[158] * v[173] - v[160] * v[174] - v[162] * v[175]
//			+ gti[0] * v[435] + gti[1] * v[437] + gti[2] * v[438];
//		v[1543] = v[215] * v[447];
//		v[1560] = v[1541] + v[1543];
//		v[446] = v[130] * v[141] + v[132] * v[142] + v[134] * v[143] - v[158] * v[169] - v[160] * v[170] - v[162] * v[171]
//			+ gti[0] * v[431] + gti[1] * v[432] + gti[2] * v[433];
//		v[1542] = -(v[214] * v[446]);
//		v[1559] = -v[1541] + v[1542];
//		v[1558] = v[1542] - v[1543];
//		v[445] = -(v[1560] * v[214]) + v[446] * v[805];
//		v[449] = v[1559] * v[215] + v[447] * v[807];
//		v[450] = v[1558] * v[216] + v[448] * v[809];
//	}
//	else {
//		v[445] = 0e0;
//		v[449] = 0e0;
//		v[450] = 0e0;
//	};
//	v[452] = v[1057] * v[451] + v[182] * v[453] + v[185] * v[454] + v[188] * v[455] + v[191] * v[456] + v[194] * v[457]
//		+ v[197] * v[458] + v[214] * v[794];
//	v[466] = v[181] * v[453] + v[184] * v[454] + v[187] * v[455] + v[190] * v[456] + v[193] * v[457] + v[196] * v[458]
//		+ v[1060] * v[459] + v[215] * v[794];
//	v[468] = v[1066] * v[216] + v[467] * v[945];
//	(*vnrel) = sqrt((v[452] * v[452]) + (v[466] * v[466]) + (v[468] * v[468]));
//	v[474] = -((v[1545] * Power(v[475], v[490])) / ((*n2)*Power((*gnbb), v[495])));
//	b478 = v[477] < (*gnb);
//	if (b478) {
//		b480 = v[477] > (*gnbb);
//		if (b480) {
//			v[489] = (*gnb) - v[477];
//			v[1926] = Power(v[489], v[490]);
//			v[1921] = Power(v[489], v[1082]);
//			v[1832] = Power(v[489], v[490]);
//			v[1789] = Power(v[489], -2e0 + v[490]);
//			v[1213] = Power(v[489], v[1082]);
//			v[1706] = Power(v[489], v[490]);
//			v[483] = Power(v[489], (*n1));
//			v[1671] = -((*epsn1)*v[483]);
//			v[1616] = -((*epsn1)*v[483]);
//			v[1547] = -((*epsn1)*v[483]);
//			v[1544] = -((*epsn1)*v[483]);
//			v[482] = v[1544] * v[214];
//			v[484] = v[1544] * v[215];
//			v[485] = v[1544] * v[216];
//		}
//		else {
//			v[1930] = Power(v[477], v[495]);
//			v[1922] = Power(v[477], v[1087]);
//			v[1836] = Power(v[477], v[495]);
//			v[1794] = Power(v[477], -2e0 + v[495]);
//			v[1219] = Power(v[477], v[1087]);
//			v[1710] = Power(v[477], v[495]);
//			v[486] = -((*epsn1)*Power(v[475], (*n1))) + v[474] * (Power((*gnbb), (*n2)) - Power(v[477], (*n2)));
//			v[482] = v[214] * v[486];
//			v[484] = v[215] * v[486];
//			v[485] = v[216] * v[486];
//		};
//	}
//	else {
//		v[482] = 0e0;
//		v[484] = 0e0;
//		v[485] = 0e0;
//	};
//	v[1552] = v[485] / 2e0;
//	v[1551] = v[484] / 2e0;
//	v[2805] = v[125] * v[482];
//	v[2806] = v[125] * v[484];
//	v[2807] = v[125] * v[485];
//	v[2808] = v[127] * v[482];
//	v[2809] = v[127] * v[484];
//	v[2810] = v[127] * v[485];
//	v[2811] = v[129] * v[482];
//	v[2812] = v[129] * v[484];
//	v[2813] = v[129] * v[485];
//	v[2814] = -(v[153] * v[482]);
//	v[2815] = -(v[153] * v[484]);
//	v[2816] = -(v[153] * v[485]);
//	v[2817] = -(v[155] * v[482]);
//	v[2818] = -(v[155] * v[484]);
//	v[2819] = -(v[155] * v[485]);
//	v[2820] = -(v[157] * v[482]);
//	v[2821] = -(v[157] * v[484]);
//	v[2822] = -(v[157] * v[485]);
//	v[3025] = 0e0;
//	v[3026] = 0e0;
//	v[3027] = 0e0;
//	v[3028] = 0e0;
//	v[3029] = 0e0;
//	v[3030] = 0e0;
//	v[3031] = v[482];
//	v[3032] = v[484];
//	v[3033] = v[485];
//	v[3034] = 0e0;
//	v[3035] = 0e0;
//	v[3036] = 0e0;
//	v[3037] = 0e0;
//	v[3038] = 0e0;
//	v[3039] = 0e0;
//	v[3040] = 0e0;
//	v[3041] = 0e0;
//	v[3042] = 0e0;
//	v[3043] = 0e0;
//	v[3044] = 0e0;
//	v[3045] = 0e0;
//	v[3046] = v[482];
//	v[3047] = v[484];
//	v[3048] = v[485];
//	v[3049] = 0e0;
//	v[3050] = 0e0;
//	v[3051] = 0e0;
//	v[3052] = 0e0;
//	v[3053] = 0e0;
//	v[3054] = 0e0;
//	v[3055] = 0e0;
//	v[3056] = 0e0;
//	v[3057] = 0e0;
//	v[3058] = 0e0;
//	v[3059] = 0e0;
//	v[3060] = 0e0;
//	v[2989] = 0e0;
//	v[2990] = 0e0;
//	v[2991] = 0e0;
//	v[2992] = 0e0;
//	v[2993] = 0e0;
//	v[2994] = 0e0;
//	v[2995] = 0e0;
//	v[2996] = 0e0;
//	v[2997] = 0e0;
//	v[2998] = 0e0;
//	v[2999] = 0e0;
//	v[3000] = 0e0;
//	v[3001] = -v[482];
//	v[3002] = -v[484];
//	v[3003] = -v[485];
//	v[3004] = 0e0;
//	v[3005] = 0e0;
//	v[3006] = 0e0;
//	v[2971] = 0e0;
//	v[2972] = 0e0;
//	v[2973] = 0e0;
//	v[2974] = 0e0;
//	v[2975] = 0e0;
//	v[2976] = 0e0;
//	v[2977] = 0e0;
//	v[2978] = 0e0;
//	v[2979] = 0e0;
//	v[2980] = 0e0;
//	v[2981] = 0e0;
//	v[2982] = 0e0;
//	v[2983] = 0e0;
//	v[2984] = 0e0;
//	v[2985] = 0e0;
//	v[2986] = -v[482];
//	v[2987] = -v[484];
//	v[2988] = -v[485];
//	v[3061] = v[482];
//	v[3062] = v[484];
//	v[3063] = v[485];
//	v[3064] = 0e0;
//	v[3065] = 0e0;
//	v[3066] = 0e0;
//	v[3067] = 0e0;
//	v[3068] = 0e0;
//	v[3069] = 0e0;
//	v[3070] = 0e0;
//	v[3071] = 0e0;
//	v[3072] = 0e0;
//	v[3073] = 0e0;
//	v[3074] = 0e0;
//	v[3075] = 0e0;
//	v[3076] = 0e0;
//	v[3077] = 0e0;
//	v[3078] = 0e0;
//	v[3007] = 0e0;
//	v[3008] = 0e0;
//	v[3009] = 0e0;
//	v[3010] = 0e0;
//	v[3011] = 0e0;
//	v[3012] = 0e0;
//	v[3013] = 0e0;
//	v[3014] = 0e0;
//	v[3015] = 0e0;
//	v[3016] = -v[482];
//	v[3017] = -v[484];
//	v[3018] = -v[485];
//	v[3019] = 0e0;
//	v[3020] = 0e0;
//	v[3021] = 0e0;
//	v[3022] = 0e0;
//	v[3023] = 0e0;
//	v[3024] = 0e0;
//	v[1550] = v[482] / 2e0;
//	if (b478) {
//		if (b480) {
//			v[759] = (*meq)*v[1545];
//			v[760] = sqrt(v[759] * Power(v[489], v[490]));
//			v[1212] = (v[490] * v[759]) / v[760];
//			v[492] = 2e0*v[760] * (*zetan);
//			v[491] = v[452] * v[492];
//			v[493] = v[466] * v[492];
//			v[494] = v[468] * v[492];
//		}
//		else {
//			v[764] = -((*meq)*(*n2)*v[474]);
//			v[765] = sqrt(v[764] * Power(v[477], v[495]));
//			v[1218] = (v[495] * v[764]) / v[765];
//			v[496] = 2e0*v[765] * (*zetan);
//			v[491] = v[452] * v[496];
//			v[493] = v[466] * v[496];
//			v[494] = v[468] * v[496];
//		};
//		b497 = v[477] < (*gnb) && v[214] * (v[482] + v[491]) + v[215] * (v[484] + v[493]) + v[216] * (v[485] + v[494]) < 0e0;
//		if (b497) {
//			v[499] = v[491];
//			v[500] = v[493];
//			v[501] = v[494];
//		}
//		else {
//			v[499] = -v[482];
//			v[500] = -v[484];
//			v[501] = -v[485];
//		};
//	}
//	else {
//		v[499] = 0e0;
//		v[500] = 0e0;
//		v[501] = 0e0;
//	};
//	v[502] = v[482] + v[499];
//	v[503] = v[484] + v[500];
//	v[504] = v[485] + v[501];
//	v[1127] = (v[502] * v[502]) + (v[503] * v[503]) + (v[504] * v[504]);
//	v[505] = (*epst)*v[445];
//	v[506] = (*epst)*v[449];
//	v[507] = (*epst)*v[450];
//	v[511] = v[505] - (*ct)*(v[181] * v[644] + v[182] * v[648] + v[183] * v[652] + v[184] * v[656] + v[185] * v[660]
//		+ v[186] * v[664] + v[187] * v[668] + v[188] * v[672] + v[189] * v[676] + v[190] * v[680] + v[191] * v[684] + v[192] * v[688]
//		+ v[193] * v[692] + v[194] * v[696] + v[195] * v[700] + v[196] * v[704] + v[197] * v[708] + v[198] * v[712]);
//	v[512] = v[506] - (*ct)*(v[181] * v[645] + v[182] * v[649] + v[183] * v[653] + v[184] * v[657] + v[185] * v[661]
//		+ v[186] * v[665] + v[187] * v[669] + v[188] * v[673] + v[189] * v[677] + v[190] * v[681] + v[191] * v[685] + v[192] * v[689]
//		+ v[193] * v[693] + v[194] * v[697] + v[195] * v[701] + v[196] * v[705] + v[197] * v[709] + v[198] * v[713]);
//	v[513] = v[507] - (*ct)*(v[181] * v[646] + v[182] * v[650] + v[183] * v[654] + v[184] * v[658] + v[185] * v[662]
//		+ v[186] * v[666] + v[187] * v[670] + v[188] * v[674] + v[189] * v[678] + v[190] * v[682] + v[191] * v[686] + v[192] * v[690]
//		+ v[193] * v[694] + v[194] * v[698] + v[195] * v[702] + v[196] * v[706] + v[197] * v[710] + v[198] * v[714]);
//	v[1124] = (v[511] * v[511]) + (v[512] * v[512]) + (v[513] * v[513]);
//	if (b478) {
//		if ((*stick)) {
//			b516 = sqrt((v[511] * v[511]) + (v[512] * v[512]) + (v[513] * v[513])) <= (*mus)*sqrt((v[502] * v[502]) +
//				(v[503] * v[503]) + (v[504] * v[504]));
//			if (b516) {
//				v[518] = v[511];
//				v[519] = v[512];
//				v[520] = v[513];
//				v[521] = 1e0;
//			}
//			else {
//				v[1126] = 1e0 / sqrt(v[1124]);
//				v[1128] = 1e0 / sqrt(v[1127]);
//				v[531] = sqrt(v[1127]);
//				v[522] = sqrt(v[1124]);
//				if (v[522] > 0.1e-5) { v010 = 1e0 / v[522]; v011 = (-(v010 / v[522])); v012 = (2e0*v010) / Power(v[522], 2); }
//				else {
//					v010 = (24000000e0 - (-1e0 + 1000000e0*v[522])*(71999994e0 - 0.71999982e14*v[522] + 0.6e19*Power(v[522], 3)
//						+ 0.23999982e20*(v[522] * v[522]))) / 24e0;
//					v011 = (-500000e0 / 3e0)*(35999994e0 - 0.47999982e14*v[522] + 0.6e19*Power(v[522], 3) + 0.17999982e20*
//						(v[522] * v[522]));
//					v012 = 0.1e13*(7999997e0 - 0.5999994e13*v[522] - 0.3e13*(v[522] * v[522]));
//				};
//				v[526] = v011;
//				v[527] = v010;
//				v[1125] = (*mud)*v[527] * v[531];
//				v[518] = v[1125] * v[511];
//				v[519] = v[1125] * v[512];
//				v[520] = v[1125] * v[513];
//				v[521] = 0e0;
//			};
//			if (sqrt((v[505] * v[505]) + (v[506] * v[506]) + (v[507] * v[507])) > (*mus)*sqrt((v[502] * v[502]) +
//				(v[503] * v[503]) + (v[504] * v[504]))) {
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
//				v[540] = sqrt((v[505] * v[505]) + (v[506] * v[506]) + (v[507] * v[507]));
//				if (v[540] > 0.1e-5) { v016 = 1e0 / v[540]; v017 = (-(v016 / v[540])); v018 = (2e0*v016) / Power(v[540], 2); }
//				else {
//					v016 = (24000000e0 - (-1e0 + 1000000e0*v[540])*(71999994e0 - 0.71999982e14*v[540] + 0.6e19*Power(v[540], 3)
//						+ 0.23999982e20*(v[540] * v[540]))) / 24e0;
//					v017 = (-500000e0 / 3e0)*(35999994e0 - 0.47999982e14*v[540] + 0.6e19*Power(v[540], 3) + 0.17999982e20*
//						(v[540] * v[540]));
//					v018 = 0.1e13*(7999997e0 - 0.5999994e13*v[540] - 0.3e13*(v[540] * v[540]));
//				};
//				v[547] = -((*mud)*v013*v016*sqrt(v[1127]));
//				v[546] = v[445] + v[505] * v[547];
//				v[548] = v[449] + v[506] * v[547];
//				v[549] = v[450] + v[507] * v[547];
//			}
//			else {
//				v[546] = 0e0;
//				v[548] = 0e0;
//				v[549] = 0e0;
//			};
//		}
//		else {
//			b550 = sqrt((v[511] * v[511]) + (v[512] * v[512]) + (v[513] * v[513])) <= (*mud)*sqrt((v[502] * v[502]) +
//				(v[503] * v[503]) + (v[504] * v[504]));
//			if (b550) {
//				v[518] = v[511];
//				v[519] = v[512];
//				v[520] = v[513];
//				v[521] = 1e0;
//			}
//			else {
//				v[1134] = 1e0 / sqrt(v[1124]);
//				v[1138] = 1e0 / sqrt(v[1127]);
//				v[561] = sqrt(v[1127]);
//				v[1610] = (*mud)*v[561];
//				v[552] = sqrt(v[1124]);
//				if (v[552] > 0.1e-5) { v019 = 1e0 / v[552]; v020 = (-(v019 / v[552])); v021 = (2e0*v019) / Power(v[552], 2); }
//				else {
//					v019 = (24000000e0 - (-1e0 + 1000000e0*v[552])*(71999994e0 - 0.71999982e14*v[552] + 0.6e19*Power(v[552], 3)
//						+ 0.23999982e20*(v[552] * v[552]))) / 24e0;
//					v020 = (-500000e0 / 3e0)*(35999994e0 - 0.47999982e14*v[552] + 0.6e19*Power(v[552], 3) + 0.17999982e20*
//						(v[552] * v[552]));
//					v021 = 0.1e13*(7999997e0 - 0.5999994e13*v[552] - 0.3e13*(v[552] * v[552]));
//				};
//				v[556] = v020;
//				v[557] = v019;
//				v[1546] = (*mud)*v[557] * v[561];
//				v[518] = v[1546] * v[511];
//				v[519] = v[1546] * v[512];
//				v[520] = v[1546] * v[513];
//				v[521] = 0e0;
//			};
//			if (sqrt((v[505] * v[505]) + (v[506] * v[506]) + (v[507] * v[507])) > (*mud)*sqrt((v[502] * v[502]) +
//				(v[503] * v[503]) + (v[504] * v[504]))) {
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
//				v[570] = sqrt((v[505] * v[505]) + (v[506] * v[506]) + (v[507] * v[507]));
//				if (v[570] > 0.1e-5) { v025 = 1e0 / v[570]; v026 = (-(v025 / v[570])); v027 = (2e0*v025) / Power(v[570], 2); }
//				else {
//					v025 = (24000000e0 - (-1e0 + 1000000e0*v[570])*(71999994e0 - 0.71999982e14*v[570] + 0.6e19*Power(v[570], 3)
//						+ 0.23999982e20*(v[570] * v[570]))) / 24e0;
//					v026 = (-500000e0 / 3e0)*(35999994e0 - 0.47999982e14*v[570] + 0.6e19*Power(v[570], 3) + 0.17999982e20*
//						(v[570] * v[570]));
//					v027 = 0.1e13*(7999997e0 - 0.5999994e13*v[570] - 0.3e13*(v[570] * v[570]));
//				};
//				v[576] = -((*mud)*v022*v025*sqrt(v[1127]));
//				v[546] = v[445] + v[505] * v[576];
//				v[548] = v[449] + v[506] * v[576];
//				v[549] = v[450] + v[507] * v[576];
//			}
//			else {
//				v[546] = 0e0;
//				v[548] = 0e0;
//				v[549] = 0e0;
//			};
//		};
//	}
//	else {
//		v[518] = 0e0;
//		v[519] = 0e0;
//		v[520] = 0e0;
//	};
//	fn[0] = v[502];
//	fn[1] = v[503];
//	fn[2] = v[504];
//	ft[0] = v[518];
//	ft[1] = v[519];
//	ft[2] = v[520];
//	(*stickupdated) = v[521];
//	gtpupdated[0] = v[445] - v[546];
//	gtpupdated[1] = v[449] - v[548];
//	gtpupdated[2] = v[450] - v[549];
//	b590 = b478;
//	if (b590) {
//		b591 = b480;
//	}
//	else {
//	};
//	v[595] = -(v[169] * v[482]) - v[173] * v[484] - v[177] * v[485];
//	v[596] = v[141] * v[482] + v[145] * v[484] + v[149] * v[485];
//	v[597] = (-(v[171] * v[482]) - v[175] * v[484] - v[179] * v[485] - v[595]) / 2e0;
//	v[598] = (-(v[170] * v[482]) - v[174] * v[484] - v[178] * v[485] - v[595]) / 2e0;
//	v[599] = (v[142] * v[482] + v[146] * v[484] + v[150] * v[485] - v[596]) / 2e0;
//	v[600] = (v[143] * v[482] + v[147] * v[484] + v[151] * v[485] - v[596]) / 2e0;
//	for (i585 = 1; i585 <= 18; i585++) {
//		v[609] = v[2876 + i585];
//		v[608] = v[2858 + i585];
//		v[607] = v[2840 + i585];
//		v[606] = v[2822 + i585];
//		v[610] = (-v[606] - v[607]) / 2e0;
//		v[611] = (-v[608] - v[609]) / 2e0;
//		v[612] = (2e0*v[2916 + i585] + v[143] * v[606] + v[142] * v[607] - v[170] * v[608] - v[171] * v[609] + 2e0*v[141] * v[610]
//			- 2e0*v[169] * v[611]) / 2e0;
//		v[623] = v[612];
//		v[613] = (2e0*v[2934 + i585] + v[147] * v[606] + v[146] * v[607] - v[174] * v[608] - v[175] * v[609] + 2e0*v[145] * v[610]
//			- 2e0*v[173] * v[611]) / 2e0;
//		v[622] = v[613];
//		v[614] = (2e0*v[2952 + i585] + v[151] * v[606] + v[150] * v[607] - v[178] * v[608] - v[179] * v[609] + 2e0*v[149] * v[610]
//			- 2e0*v[177] * v[611]) / 2e0;
//		v[621] = v[614];
//		v[615] = 0e0;
//		v[616] = 0e0;
//		v[617] = 0e0;
//		v[618] = 0e0;
//		b619 = b478;
//		if (b619) {
//			v[1548] = -(v[216] * v[621]) - v[215] * v[622] - v[214] * v[623];
//			b620 = b480;
//			if (b620) {
//				v[617] = v[1547] * v[614];
//				v[614] = 0e0;
//				v[616] = v[1547] * v[613];
//				v[613] = 0e0;
//				v[615] = v[1547] * v[612];
//				v[612] = 0e0;
//				v[618] = -(v[1545] * v[1548] * v[1706]);
//			}
//			else {
//				v[617] = v[486] * v[621];
//				v[614] = 0e0;
//				v[616] = v[486] * v[622];
//				v[613] = 0e0;
//				v[615] = v[486] * v[623];
//				v[612] = 0e0;
//				v[618] = (*n2)*v[1548] * v[1710] * v[474];
//			};
//		}
//		else {
//		};
//		v[618] = v[212] * (v[199] * v[615] + v[200] * v[616] + v[201] * v[617]) + v[618];
//		v[1549] = v[618] / v[477];
//		v[628] = v[1549] * v[201] + v[213] * v[617];
//		v[630] = v[1549] * v[200] + v[213] * v[616];
//		v[631] = v[1549] * v[199] + v[213] * v[615];
//		v[3079] = v[482] * v[610] + v[125] * v[631];
//		v[3080] = v[484] * v[610] + v[125] * v[630];
//		v[3081] = v[485] * v[610] + v[125] * v[628];
//		v[3082] = v[1550] * v[607] + v[127] * v[631];
//		v[3083] = v[1551] * v[607] + v[127] * v[630];
//		v[3084] = v[1552] * v[607] + v[127] * v[628];
//		v[3085] = v[1550] * v[606] + v[129] * v[631];
//		v[3086] = v[1551] * v[606] + v[129] * v[630];
//		v[3087] = v[1552] * v[606] + v[129] * v[628];
//		v[3088] = -(v[482] * v[611]) - v[153] * v[631];
//		v[3089] = -(v[484] * v[611]) - v[153] * v[630];
//		v[3090] = -(v[485] * v[611]) - v[153] * v[628];
//		v[3091] = -(v[1550] * v[608]) - v[155] * v[631];
//		v[3092] = -(v[1551] * v[608]) - v[155] * v[630];
//		v[3093] = -(v[1552] * v[608]) - v[155] * v[628];
//		v[3094] = -(v[1550] * v[609]) - v[157] * v[631];
//		v[3095] = -(v[1551] * v[609]) - v[157] * v[630];
//		v[3096] = -(v[1552] * v[609]) - v[157] * v[628];
//		v[632] = v[3006 + i585] - v[177] * v[628] - v[173] * v[630] - v[169] * v[631];
//		v[633] = v[3060 + i585] + v[149] * v[628] + v[145] * v[630] + v[141] * v[631];
//		v[634] = (v[2970 + i585] - v[179] * v[628] - v[175] * v[630] - v[171] * v[631] - v[632]) / 2e0;
//		v[635] = (v[2988 + i585] - v[178] * v[628] - v[174] * v[630] - v[170] * v[631] - v[632]) / 2e0;
//		v[636] = (v[3042 + i585] + v[150] * v[628] + v[146] * v[630] + v[142] * v[631] - v[633]) / 2e0;
//		v[637] = (v[3024 + i585] + v[151] * v[628] + v[147] * v[630] + v[143] * v[631] - v[633]) / 2e0;
//		Rc[i585 - 1] += v[2804 + i585] + v[600] * v[606] + v[599] * v[607] + v[598] * v[608] + v[597] * v[609];
//		for (i603 = 1; i603 <= 18; i603++) {
//			Kc[i585 - 1][i603 - 1] += v[3078 + i603] + v[2876 + i603] * v[634] + v[2858 + i603] * v[635] + v[2840 + i603] * v[636]
//				+ v[2822 + i603] * v[637];
//		};/* end for */
//	};/* end for */
//	v[719] = 0e0;
//	v[720] = 0e0;
//	v[721] = 0e0;
//	b722 = b478;
//	if (b722) {
//		b723 = (*stick);
//		if (b723) {
//			b724 = b516;
//			if (b724) {
//				v[721] = 0e0;
//				v[720] = 0e0;
//				v[719] = 0e0;
//			}
//			else {
//			};
//		}
//		else {
//			b725 = b550;
//			if (b725) {
//				v[721] = 0e0;
//				v[720] = 0e0;
//				v[719] = 0e0;
//			}
//			else {
//			};
//		};
//	}
//	else {
//	};
//	v[1556] = (*ct)*v[719];
//	v[1276] = v[1553] * v[719];
//	v[1555] = (*ct)*v[720];
//	v[1281] = v[1553] * v[720];
//	v[1554] = (*ct)*v[721];
//	v[1286] = v[1553] * v[721];
//	v[1564] = (v[1554] * v[733]) / 2e0;
//	v[1565] = (v[1554] * v[735]) / 2e0;
//	v[1570] = (v[1554] * v[737]) / 2e0;
//	v[1571] = (v[1554] * v[739]) / 2e0;
//	v[1566] = (v[1555] * v[733]) / 2e0;
//	v[1567] = (v[1555] * v[735]) / 2e0;
//	v[1572] = (v[1555] * v[737]) / 2e0;
//	v[1573] = (v[1555] * v[739]) / 2e0;
//	v[1568] = (v[1556] * v[733]) / 2e0;
//	v[1569] = (v[1556] * v[735]) / 2e0;
//	v[1574] = (v[1556] * v[737]) / 2e0;
//	v[1575] = (v[1556] * v[739]) / 2e0;
//	v[745] = (*epst)*v[721];
//	v[1050] = -(v[216] * v[745]);
//	v[746] = (*epst)*v[720];
//	v[1053] = -(v[215] * v[746]);
//	v[1048] = v[1050] + v[1053];
//	v[747] = (*epst)*v[719];
//	v[1054] = -(v[214] * v[747]);
//	v[1055] = v[1053] + v[1054];
//	v[1051] = v[1050] + v[1054];
//	v[748] = 0e0;
//	v[749] = 0e0;
//	v[750] = 0e0;
//	v[751] = 0e0;
//	v[752] = 0e0;
//	b753 = b478;
//	if (b753) {
//		v[754] = 0e0;
//		v[755] = 0e0;
//		v[756] = 0e0;
//		b757 = b497;
//		if (b757) {
//			v[756] = 0e0;
//			v[755] = 0e0;
//			v[754] = 0e0;
//		}
//		else {
//		};
//		v[766] = v[452] * v[754] + v[466] * v[755] + v[468] * v[756];
//		v[1557] = v[766] * (*zetan);
//		b758 = b480;
//		if (b758) {
//			v[750] = v[492] * v[756];
//			v[749] = v[492] * v[755];
//			v[748] = v[492] * v[754];
//			v[752] = (v[1557] * v[490] * v[759] * Power(v[489], v[1082])) / v[760];
//		}
//		else {
//			v[750] = v[496] * v[756];
//			v[749] = v[496] * v[755];
//			v[748] = v[496] * v[754];
//			v[751] = (v[1557] * v[495] * v[764] * Power(v[477], v[1087])) / v[765];
//		};
//	}
//	else {
//	};
//	v[3717] = v[125] * v[748];
//	v[3718] = 0e0;
//	v[3719] = 0e0;
//	v[3720] = v[127] * v[748];
//	v[3721] = 0e0;
//	v[3722] = 0e0;
//	v[3723] = v[129] * v[748];
//	v[3724] = 0e0;
//	v[3725] = 0e0;
//	v[3726] = -(v[153] * v[748]);
//	v[3727] = 0e0;
//	v[3728] = 0e0;
//	v[3729] = -(v[155] * v[748]);
//	v[3730] = 0e0;
//	v[3731] = 0e0;
//	v[3732] = -(v[157] * v[748]);
//	v[3733] = 0e0;
//	v[3734] = 0e0;
//	v[1300] = -(v[451] * v[748]);
//	v[3663] = 0e0;
//	v[3664] = 0e0;
//	v[3665] = 0e0;
//	v[3666] = 0e0;
//	v[3667] = 0e0;
//	v[3668] = 0e0;
//	v[3669] = 0e0;
//	v[3670] = 0e0;
//	v[3671] = 0e0;
//	v[3672] = 0e0;
//	v[3673] = 0e0;
//	v[3674] = 0e0;
//	v[3675] = 0e0;
//	v[3676] = 0e0;
//	v[3677] = 0e0;
//	v[3678] = v[749];
//	v[3679] = v[748];
//	v[3680] = 0e0;
//	v[3645] = 0e0;
//	v[3646] = 0e0;
//	v[3647] = 0e0;
//	v[3648] = 0e0;
//	v[3649] = 0e0;
//	v[3650] = 0e0;
//	v[3651] = 0e0;
//	v[3652] = 0e0;
//	v[3653] = 0e0;
//	v[3654] = 0e0;
//	v[3655] = 0e0;
//	v[3656] = 0e0;
//	v[3657] = v[749];
//	v[3658] = v[748];
//	v[3659] = 0e0;
//	v[3660] = 0e0;
//	v[3661] = 0e0;
//	v[3662] = 0e0;
//	v[3627] = 0e0;
//	v[3628] = 0e0;
//	v[3629] = 0e0;
//	v[3630] = 0e0;
//	v[3631] = 0e0;
//	v[3632] = 0e0;
//	v[3633] = 0e0;
//	v[3634] = 0e0;
//	v[3635] = 0e0;
//	v[3636] = v[749];
//	v[3637] = v[748];
//	v[3638] = 0e0;
//	v[3639] = 0e0;
//	v[3640] = 0e0;
//	v[3641] = 0e0;
//	v[3642] = 0e0;
//	v[3643] = 0e0;
//	v[3644] = 0e0;
//	v[3609] = 0e0;
//	v[3610] = 0e0;
//	v[3611] = 0e0;
//	v[3612] = 0e0;
//	v[3613] = 0e0;
//	v[3614] = 0e0;
//	v[3615] = v[749];
//	v[3616] = v[748];
//	v[3617] = 0e0;
//	v[3618] = 0e0;
//	v[3619] = 0e0;
//	v[3620] = 0e0;
//	v[3621] = 0e0;
//	v[3622] = 0e0;
//	v[3623] = 0e0;
//	v[3624] = 0e0;
//	v[3625] = 0e0;
//	v[3626] = 0e0;
//	v[3591] = 0e0;
//	v[3592] = 0e0;
//	v[3593] = 0e0;
//	v[3594] = v[749];
//	v[3595] = v[748];
//	v[3596] = 0e0;
//	v[3597] = 0e0;
//	v[3598] = 0e0;
//	v[3599] = 0e0;
//	v[3600] = 0e0;
//	v[3601] = 0e0;
//	v[3602] = 0e0;
//	v[3603] = 0e0;
//	v[3604] = 0e0;
//	v[3605] = 0e0;
//	v[3606] = 0e0;
//	v[3607] = 0e0;
//	v[3608] = 0e0;
//	v[3573] = v[749];
//	v[3574] = v[748];
//	v[3575] = 0e0;
//	v[3576] = 0e0;
//	v[3577] = 0e0;
//	v[3578] = 0e0;
//	v[3579] = 0e0;
//	v[3580] = 0e0;
//	v[3581] = 0e0;
//	v[3582] = 0e0;
//	v[3583] = 0e0;
//	v[3584] = 0e0;
//	v[3585] = 0e0;
//	v[3586] = 0e0;
//	v[3587] = 0e0;
//	v[3588] = 0e0;
//	v[3589] = 0e0;
//	v[3590] = 0e0;
//	v[3681] = 0e0;
//	v[3682] = v[125] * v[749];
//	v[3683] = 0e0;
//	v[3684] = 0e0;
//	v[3685] = v[127] * v[749];
//	v[3686] = 0e0;
//	v[3687] = 0e0;
//	v[3688] = v[129] * v[749];
//	v[3689] = 0e0;
//	v[3690] = 0e0;
//	v[3691] = -(v[153] * v[749]);
//	v[3692] = 0e0;
//	v[3693] = 0e0;
//	v[3694] = -(v[155] * v[749]);
//	v[3695] = 0e0;
//	v[3696] = 0e0;
//	v[3697] = -(v[157] * v[749]);
//	v[3698] = 0e0;
//	v[1301] = -(v[459] * v[749]);
//	v[1606] = v[216] * v[750];
//	v[1304] = -(v[215] * v[750]);
//	v[1303] = -(v[214] * v[750]);
//	v[1302] = -(v[467] * v[750]);
//	v[3119] = -(v[125] * v[1300]) + v[453] * v[749];
//	v[3120] = -(v[125] * v[1301]) + v[453] * v[748];
//	v[3121] = -(v[125] * v[1302]);
//	v[3122] = -(v[127] * v[1300]) + v[454] * v[749];
//	v[3123] = -(v[127] * v[1301]) + v[454] * v[748];
//	v[3124] = -(v[127] * v[1302]);
//	v[3125] = -(v[129] * v[1300]) + v[455] * v[749];
//	v[3126] = -(v[129] * v[1301]) + v[455] * v[748];
//	v[3127] = -(v[129] * v[1302]);
//	v[3128] = v[1300] * v[153] + v[456] * v[749];
//	v[3129] = v[1301] * v[153] + v[456] * v[748];
//	v[3130] = v[1302] * v[153];
//	v[3131] = v[1300] * v[155] + v[457] * v[749];
//	v[3132] = v[1301] * v[155] + v[457] * v[748];
//	v[3133] = v[1302] * v[155];
//	v[3134] = v[1300] * v[157] + v[458] * v[749];
//	v[3135] = v[1301] * v[157] + v[458] * v[748];
//	v[3136] = v[1302] * v[157];
//	v[3861] = 0e0;
//	v[3862] = 0e0;
//	v[3863] = 0e0;
//	v[3864] = 0e0;
//	v[3865] = 0e0;
//	v[3866] = 0e0;
//	v[3867] = -v[1300];
//	v[3868] = -v[1301];
//	v[3869] = -v[1302];
//	v[3870] = 0e0;
//	v[3871] = 0e0;
//	v[3872] = 0e0;
//	v[3873] = 0e0;
//	v[3874] = 0e0;
//	v[3875] = 0e0;
//	v[3876] = 0e0;
//	v[3877] = 0e0;
//	v[3878] = 0e0;
//	v[3897] = 0e0;
//	v[3898] = 0e0;
//	v[3899] = 0e0;
//	v[3900] = -v[1300];
//	v[3901] = -v[1301];
//	v[3902] = -v[1302];
//	v[3903] = 0e0;
//	v[3904] = 0e0;
//	v[3905] = 0e0;
//	v[3906] = 0e0;
//	v[3907] = 0e0;
//	v[3908] = 0e0;
//	v[3909] = 0e0;
//	v[3910] = 0e0;
//	v[3911] = 0e0;
//	v[3912] = 0e0;
//	v[3913] = 0e0;
//	v[3914] = 0e0;
//	v[3789] = 0e0;
//	v[3790] = 0e0;
//	v[3791] = 0e0;
//	v[3792] = 0e0;
//	v[3793] = 0e0;
//	v[3794] = 0e0;
//	v[3795] = 0e0;
//	v[3796] = 0e0;
//	v[3797] = 0e0;
//	v[3798] = 0e0;
//	v[3799] = 0e0;
//	v[3800] = 0e0;
//	v[3801] = v[1300];
//	v[3802] = v[1301];
//	v[3803] = v[1302];
//	v[3804] = 0e0;
//	v[3805] = 0e0;
//	v[3806] = 0e0;
//	v[3753] = 0e0;
//	v[3754] = 0e0;
//	v[3755] = 0e0;
//	v[3756] = 0e0;
//	v[3757] = 0e0;
//	v[3758] = 0e0;
//	v[3759] = 0e0;
//	v[3760] = 0e0;
//	v[3761] = 0e0;
//	v[3762] = 0e0;
//	v[3763] = 0e0;
//	v[3764] = 0e0;
//	v[3765] = 0e0;
//	v[3766] = 0e0;
//	v[3767] = 0e0;
//	v[3768] = v[1300];
//	v[3769] = v[1301];
//	v[3770] = v[1302];
//	v[3933] = -v[1300];
//	v[3934] = -v[1301];
//	v[3935] = -v[1302];
//	v[3936] = 0e0;
//	v[3937] = 0e0;
//	v[3938] = 0e0;
//	v[3939] = 0e0;
//	v[3940] = 0e0;
//	v[3941] = 0e0;
//	v[3942] = 0e0;
//	v[3943] = 0e0;
//	v[3944] = 0e0;
//	v[3945] = 0e0;
//	v[3946] = 0e0;
//	v[3947] = 0e0;
//	v[3948] = 0e0;
//	v[3949] = 0e0;
//	v[3950] = 0e0;
//	v[3825] = 0e0;
//	v[3826] = 0e0;
//	v[3827] = 0e0;
//	v[3828] = 0e0;
//	v[3829] = 0e0;
//	v[3830] = 0e0;
//	v[3831] = 0e0;
//	v[3832] = 0e0;
//	v[3833] = 0e0;
//	v[3834] = v[1300];
//	v[3835] = v[1301];
//	v[3836] = v[1302];
//	v[3837] = 0e0;
//	v[3838] = 0e0;
//	v[3839] = 0e0;
//	v[3840] = 0e0;
//	v[3841] = 0e0;
//	v[3842] = 0e0;
//	v[939] = -(v[157] * v[750]);
//	v[938] = -(v[155] * v[750]);
//	v[937] = -(v[153] * v[750]);
//	v[936] = v[129] * v[750];
//	v[935] = v[127] * v[750];
//	v[934] = v[125] * v[750];
//	v[3555] = 0e0;
//	v[3556] = 0e0;
//	v[3557] = v[934];
//	v[3558] = 0e0;
//	v[3559] = 0e0;
//	v[3560] = v[935];
//	v[3561] = 0e0;
//	v[3562] = 0e0;
//	v[3563] = v[936];
//	v[3564] = 0e0;
//	v[3565] = 0e0;
//	v[3566] = v[937];
//	v[3567] = 0e0;
//	v[3568] = 0e0;
//	v[3569] = v[938];
//	v[3570] = 0e0;
//	v[3571] = 0e0;
//	v[3572] = v[939];
//	v[3735] = v[934];
//	v[3736] = 0e0;
//	v[3737] = 0e0;
//	v[3738] = v[935];
//	v[3739] = 0e0;
//	v[3740] = 0e0;
//	v[3741] = v[936];
//	v[3742] = 0e0;
//	v[3743] = 0e0;
//	v[3744] = v[937];
//	v[3745] = 0e0;
//	v[3746] = 0e0;
//	v[3747] = v[938];
//	v[3748] = 0e0;
//	v[3749] = 0e0;
//	v[3750] = v[939];
//	v[3751] = 0e0;
//	v[3752] = 0e0;
//	v[3699] = 0e0;
//	v[3700] = v[934];
//	v[3701] = 0e0;
//	v[3702] = 0e0;
//	v[3703] = v[935];
//	v[3704] = 0e0;
//	v[3705] = 0e0;
//	v[3706] = v[936];
//	v[3707] = 0e0;
//	v[3708] = 0e0;
//	v[3709] = v[937];
//	v[3710] = 0e0;
//	v[3711] = 0e0;
//	v[3712] = v[938];
//	v[3713] = 0e0;
//	v[3714] = 0e0;
//	v[3715] = v[939];
//	v[3716] = 0e0;
//	v[767] = 0e0;
//	v[768] = 0e0;
//	v[769] = 0e0;
//	b770 = b478;
//	if (b770) {
//		b771 = b480;
//		if (b771) {
//			v[769] = 0e0;
//			v[768] = 0e0;
//			v[767] = 0e0;
//			v[751] = v[751] - v[752];
//		}
//		else {
//		};
//	}
//	else {
//	};
//	v[772] = v[750] * v[945];
//	v[769] = v[1066] * v[750] + v[769];
//	v[785] = v[1060] * v[749];
//	v[768] = v[1065] * v[750] + v[768] + v[749] * v[794];
//	v[786] = v[214] * v[748] + v[215] * v[749];
//	v[3879] = 0e0;
//	v[3880] = 0e0;
//	v[3881] = 0e0;
//	v[3882] = 0e0;
//	v[3883] = 0e0;
//	v[3884] = 0e0;
//	v[3885] = -v[1303];
//	v[3886] = -v[1304];
//	v[3887] = v[786];
//	v[3888] = 0e0;
//	v[3889] = 0e0;
//	v[3890] = 0e0;
//	v[3891] = 0e0;
//	v[3892] = 0e0;
//	v[3893] = 0e0;
//	v[3894] = 0e0;
//	v[3895] = 0e0;
//	v[3896] = 0e0;
//	v[3915] = 0e0;
//	v[3916] = 0e0;
//	v[3917] = 0e0;
//	v[3918] = -v[1303];
//	v[3919] = -v[1304];
//	v[3920] = v[786];
//	v[3921] = 0e0;
//	v[3922] = 0e0;
//	v[3923] = 0e0;
//	v[3924] = 0e0;
//	v[3925] = 0e0;
//	v[3926] = 0e0;
//	v[3927] = 0e0;
//	v[3928] = 0e0;
//	v[3929] = 0e0;
//	v[3930] = 0e0;
//	v[3931] = 0e0;
//	v[3932] = 0e0;
//	v[3807] = 0e0;
//	v[3808] = 0e0;
//	v[3809] = 0e0;
//	v[3810] = 0e0;
//	v[3811] = 0e0;
//	v[3812] = 0e0;
//	v[3813] = 0e0;
//	v[3814] = 0e0;
//	v[3815] = 0e0;
//	v[3816] = 0e0;
//	v[3817] = 0e0;
//	v[3818] = 0e0;
//	v[3819] = v[1303];
//	v[3820] = v[1304];
//	v[3821] = -v[786];
//	v[3822] = 0e0;
//	v[3823] = 0e0;
//	v[3824] = 0e0;
//	v[3771] = 0e0;
//	v[3772] = 0e0;
//	v[3773] = 0e0;
//	v[3774] = 0e0;
//	v[3775] = 0e0;
//	v[3776] = 0e0;
//	v[3777] = 0e0;
//	v[3778] = 0e0;
//	v[3779] = 0e0;
//	v[3780] = 0e0;
//	v[3781] = 0e0;
//	v[3782] = 0e0;
//	v[3783] = 0e0;
//	v[3784] = 0e0;
//	v[3785] = 0e0;
//	v[3786] = v[1303];
//	v[3787] = v[1304];
//	v[3788] = -v[786];
//	v[3951] = -v[1303];
//	v[3952] = -v[1304];
//	v[3953] = v[786];
//	v[3954] = 0e0;
//	v[3955] = 0e0;
//	v[3956] = 0e0;
//	v[3957] = 0e0;
//	v[3958] = 0e0;
//	v[3959] = 0e0;
//	v[3960] = 0e0;
//	v[3961] = 0e0;
//	v[3962] = 0e0;
//	v[3963] = 0e0;
//	v[3964] = 0e0;
//	v[3965] = 0e0;
//	v[3966] = 0e0;
//	v[3967] = 0e0;
//	v[3968] = 0e0;
//	v[3843] = 0e0;
//	v[3844] = 0e0;
//	v[3845] = 0e0;
//	v[3846] = 0e0;
//	v[3847] = 0e0;
//	v[3848] = 0e0;
//	v[3849] = 0e0;
//	v[3850] = 0e0;
//	v[3851] = 0e0;
//	v[3852] = v[1303];
//	v[3853] = v[1304];
//	v[3854] = -v[786];
//	v[3855] = 0e0;
//	v[3856] = 0e0;
//	v[3857] = 0e0;
//	v[3858] = 0e0;
//	v[3859] = 0e0;
//	v[3860] = 0e0;
//	v[3137] = -(v[125] * v[1303]);
//	v[3138] = -(v[125] * v[1304]);
//	v[3139] = v[125] * v[786];
//	v[3140] = -(v[127] * v[1303]);
//	v[3141] = -(v[127] * v[1304]);
//	v[3142] = v[127] * v[786];
//	v[3143] = -(v[129] * v[1303]);
//	v[3144] = -(v[129] * v[1304]);
//	v[3145] = v[129] * v[786];
//	v[3146] = v[1303] * v[153];
//	v[3147] = v[1304] * v[153];
//	v[3148] = -(v[153] * v[786]);
//	v[3149] = v[1303] * v[155];
//	v[3150] = v[1304] * v[155];
//	v[3151] = -(v[155] * v[786]);
//	v[3152] = v[1303] * v[157];
//	v[3153] = v[1304] * v[157];
//	v[3154] = -(v[157] * v[786]);
//	v[787] = v[182] * v[748] + v[181] * v[749];
//	v[788] = v[185] * v[748] + v[184] * v[749];
//	v[789] = v[188] * v[748] + v[187] * v[749];
//	v[790] = v[191] * v[748] + v[190] * v[749];
//	v[791] = v[194] * v[748] + v[193] * v[749];
//	v[792] = v[197] * v[748] + v[196] * v[749];
//	v[940] = v[125] * v[787] + v[127] * v[788] + v[129] * v[789] - v[153] * v[790] - v[155] * v[791] - v[157] * v[792];
//	v[793] = v[1057] * v[748];
//	v[767] = v[1064] * v[750] + v[767] + v[748] * v[794];
//	v[795] = 0e0;
//	v[796] = 0e0;
//	v[797] = 0e0;
//	v[798] = 0e0;
//	v[799] = 0e0;
//	v[800] = 0e0;
//	v[801] = 0e0;
//	v[802] = 0e0;
//	v[803] = 0e0;
//	b804 = (*previouscontact);
//	if (b804) {
//		v[772] = -(v[448] * v[745]) + v[772];
//		v[785] = -(v[447] * v[746]) + v[785];
//		v[806] = v[1048] * v[214] + v[747] * v[805];
//		v[808] = v[1051] * v[215] + v[746] * v[807];
//		v[810] = v[1055] * v[216] + v[745] * v[809];
//		v[769] = v[1055] * v[448] + v[1558] * v[745] + v[769];
//		v[768] = v[1051] * v[447] + v[1559] * v[746] + v[768];
//		v[793] = -(v[446] * v[747]) + v[793];
//		v[767] = v[1048] * v[446] - v[1560] * v[747] + v[767];
//		v[795] = gti[0] * v[806];
//		v[796] = gti[1] * v[806];
//		v[797] = gti[2] * v[806];
//		v[814] = -(v[162] * v[806]);
//		v[815] = -(v[160] * v[806]);
//		v[816] = -(v[158] * v[806]);
//		v[817] = v[134] * v[806];
//		v[818] = v[132] * v[806];
//		v[819] = v[130] * v[806];
//		v[798] = gti[0] * v[808];
//		v[799] = gti[1] * v[808];
//		v[800] = gti[2] * v[808];
//		v[820] = -(v[162] * v[808]);
//		v[821] = -(v[160] * v[808]);
//		v[822] = -(v[158] * v[808]);
//		v[823] = v[134] * v[808];
//		v[824] = v[132] * v[808];
//		v[825] = v[130] * v[808];
//		v[801] = gti[0] * v[810];
//		v[802] = gti[1] * v[810];
//		v[803] = gti[2] * v[810];
//		v[826] = -(v[162] * v[810]);
//		v[827] = -(v[160] * v[810]);
//		v[828] = -(v[158] * v[810]);
//		v[829] = v[134] * v[810];
//		v[830] = v[132] * v[810];
//		v[831] = v[130] * v[810];
//	}
//	else {
//		v[819] = 0e0;
//		v[818] = 0e0;
//		v[817] = 0e0;
//		v[825] = 0e0;
//		v[824] = 0e0;
//		v[823] = 0e0;
//		v[831] = 0e0;
//		v[830] = 0e0;
//		v[829] = 0e0;
//		v[816] = 0e0;
//		v[815] = 0e0;
//		v[814] = 0e0;
//		v[822] = 0e0;
//		v[821] = 0e0;
//		v[820] = 0e0;
//		v[828] = 0e0;
//		v[827] = 0e0;
//		v[826] = 0e0;
//	};
//	b832 = b411;
//	if (b832) {
//		v[858] = -(v[428] * v[803]) / 2e0;
//		v[857] = -(v[428] * v[799]) / 2e0;
//		v[856] = v[428] * v[802];
//		v[855] = v[428] * v[800];
//		v[851] = v[428] * v[801];
//		v[850] = v[428] * v[797];
//		v[847] = v[428] * v[798];
//		v[846] = v[428] * v[796];
//		v[833] = v[855] + v[856];
//		v[834] = v[850] + v[851];
//		v[835] = v[846] + v[847];
//		v[845] = (v[795] * v[836]) / 2e0 + v[796] * v[837] + v[797] * v[838] + v[798] * v[839] + (v[799] * v[840]) / 2e0
//			+ v[800] * v[841] + v[801] * v[842] + v[802] * v[843] + (v[803] * v[844]) / 2e0;
//		v[1014] = (-4e0*v[845]) / Power(v[848], 2) + v[857] + v[858];
//		v[1013] = v[1014] - (v[428] * v[795]) / 2e0 - v[857];
//		v[1012] = v[1013] + v[857] - v[858];
//		v[849] = (4e0*v[1012] * v[427] + v[426] * v[833] + v[425] * v[834] - 2e0*v[846] + 2e0*v[847]) / 2e0;
//		v[854] = (4e0*v[1013] * v[426] + v[427] * v[833] + v[425] * v[835] + 2e0*v[850] - 2e0*v[851]) / 2e0;
//		v[859] = (4e0*v[1014] * v[425] + v[427] * v[834] + v[426] * v[835] - 2e0*v[855] + 2e0*v[856]) / 2e0;
//		v[1561] = v[415] * v[849] + v[414] * v[854] + v[413] * v[859];
//		v[1000] = v[1561] * v[424];
//		v[997] = v[1561] * v[418];
//		v[862] = v[423] * v[997] + v[1000] / (Power(cos(v[860]), 2)*sqrt(v[1001]));
//		v[1590] = v[862] / v[416];
//		v[1562] = v[862] / v[416];
//		v[863] = v[1562] * v[415] + v[1540] * v[849];
//		v[865] = v[1562] * v[414] + v[1540] * v[854];
//		v[866] = v[1562] * v[413] + v[1540] * v[859];
//		v[767] = v[767] - v[224] * v[863] + v[225] * v[865];
//		v[769] = v[769] - v[223] * v[865] + v[224] * v[866];
//		v[768] = v[768] + v[223] * v[863] - v[225] * v[866];
//	}
//	else {
//	};
//	v[769] = v[769] + 2e0*v[216] * v[772] + v[786] * v[945];
//	v[768] = v[768] + 2e0*v[215] * v[785] + v[214] * v[940];
//	v[767] = v[767] + 2e0*v[214] * v[793] + v[215] * v[940];
//	v[1291] = v[199] * v[767] + v[200] * v[768] + v[201] * v[769];
//	v[751] = v[1291] * v[212] + v[751];
//	v[1563] = v[751] / v[477];
//	v[867] = v[1563] * v[201] + v[213] * v[769];
//	v[868] = v[1563] * v[200] + v[213] * v[768];
//	v[869] = v[1563] * v[199] + v[213] * v[767];
//	v[826] = v[1564] + v[826] - v[157] * v[867];
//	v[820] = v[1566] + v[820] - v[157] * v[868];
//	v[814] = v[1568] + v[814] - v[157] * v[869];
//	v[827] = v[1565] + v[827] - v[155] * v[867];
//	v[821] = v[1567] + v[821] - v[155] * v[868];
//	v[815] = v[1569] + v[815] - v[155] * v[869];
//	v[828] = -v[1564] - v[1565] + v[828] - v[153] * v[867];
//	v[822] = -v[1566] - v[1567] + v[822] - v[153] * v[868];
//	v[816] = -v[1568] - v[1569] + v[816] - v[153] * v[869];
//	v[829] = v[1570] + v[829] + v[129] * v[867];
//	v[823] = v[1572] + v[823] + v[129] * v[868];
//	v[817] = v[1574] + v[817] + v[129] * v[869];
//	v[830] = v[1571] + v[830] + v[127] * v[867];
//	v[824] = v[1573] + v[824] + v[127] * v[868];
//	v[818] = v[1575] + v[818] + v[127] * v[869];
//	v[831] = -v[1570] - v[1571] + v[831] + v[125] * v[867];
//	v[825] = -v[1572] - v[1573] + v[825] + v[125] * v[868];
//	v[819] = -v[1574] - v[1575] + v[819] + v[125] * v[869];
//	v[3101] = -(v[518] * v[644]) - v[519] * v[645] - v[520] * v[646] + v[819];
//	v[3102] = -(v[518] * v[648]) - v[519] * v[649] - v[520] * v[650] + v[825];
//	v[3103] = -(v[518] * v[652]) - v[519] * v[653] - v[520] * v[654] + v[831];
//	v[3104] = -(v[518] * v[656]) - v[519] * v[657] - v[520] * v[658] + v[818];
//	v[3105] = -(v[518] * v[660]) - v[519] * v[661] - v[520] * v[662] + v[824];
//	v[3106] = -(v[518] * v[664]) - v[519] * v[665] - v[520] * v[666] + v[830];
//	v[3107] = -(v[518] * v[668]) - v[519] * v[669] - v[520] * v[670] + v[817];
//	v[3108] = -(v[518] * v[672]) - v[519] * v[673] - v[520] * v[674] + v[823];
//	v[3109] = -(v[518] * v[676]) - v[519] * v[677] - v[520] * v[678] + v[829];
//	v[3110] = -(v[518] * v[680]) - v[519] * v[681] - v[520] * v[682] + v[816];
//	v[3111] = -(v[518] * v[684]) - v[519] * v[685] - v[520] * v[686] + v[822];
//	v[3112] = -(v[518] * v[688]) - v[519] * v[689] - v[520] * v[690] + v[828];
//	v[3113] = -(v[518] * v[692]) - v[519] * v[693] - v[520] * v[694] + v[815];
//	v[3114] = -(v[518] * v[696]) - v[519] * v[697] - v[520] * v[698] + v[821];
//	v[3115] = -(v[518] * v[700]) - v[519] * v[701] - v[520] * v[702] + v[827];
//	v[3116] = -(v[518] * v[704]) - v[519] * v[705] - v[520] * v[706] + v[814];
//	v[3117] = -(v[518] * v[708]) - v[519] * v[709] - v[520] * v[710] + v[820];
//	v[3118] = -(v[518] * v[712]) - v[519] * v[713] - v[520] * v[714] + v[826];
//	for (i717 = 1; i717 <= 18; i717++) {
//		i1602 = (i717 == 16 ? 1 : 0);
//		i1601 = (i717 == 13 ? 1 : 0);
//		i1600 = (i717 == 7 ? 1 : 0);
//		i1599 = (i717 == 4 ? 1 : 0);
//		i1598 = (i717 == 17 ? 1 : 0);
//		i1597 = (i717 == 14 ? 1 : 0);
//		i1596 = (i717 == 8 ? 1 : 0);
//		i1595 = (i717 == 5 ? 1 : 0);
//		i1594 = (i717 == 18 ? 1 : 0);
//		i1593 = (i717 == 15 ? 1 : 0);
//		i1592 = (i717 == 9 ? 1 : 0);
//		i1591 = (i717 == 6 ? 1 : 0);
//		i1581 = (i717 == 12 ? 1 : 0);
//		i1580 = (i717 == 3 ? 1 : 0);
//		i1579 = (i717 == 11 ? 1 : 0);
//		i1578 = (i717 == 2 ? 1 : 0);
//		i1577 = (i717 == 10 ? 1 : 0);
//		i1576 = (i717 == 1 ? 1 : 0);
//		v[948] = v[3136 + i717];
//		v[872] = v[2840 + i717];
//		v[873] = v[2822 + i717];
//		v[874] = v[2858 + i717];
//		v[875] = v[2876 + i717];
//		v[1092] = -i1576 / 2e0;
//		v[1093] = -i1578 / 2e0;
//		v[1094] = -i1580 / 2e0;
//		v[1097] = -i1577 / 2e0;
//		v[1098] = -i1579 / 2e0;
//		v[1099] = -i1581 / 2e0;
//		v[913] = i1576 * v[125] + i1599 * v[127] + i1600 * v[129] - i1577 * v[153] - i1601 * v[155] - i1602 * v[157];
//		v[915] = i1578 * v[125] + i1595 * v[127] + i1596 * v[129] - i1579 * v[153] - i1597 * v[155] - i1598 * v[157];
//		v[917] = i1580 * v[125] + i1591 * v[127] + i1592 * v[129] - i1581 * v[153] - i1593 * v[155] - i1594 * v[157];
//		v[1582] = v[199] * v[913] + v[200] * v[915] + v[201] * v[917];
//		v[920] = -(v[1582] * v[751] * v[919]);
//		v[921] = v[1582] / v[477];
//		v[1583] = v[212] * v[921];
//		v[1071] = v[921];
//		v[946] = v[1583] * v[201] + v[213] * v[917];
//		v[926] = v[1583] * v[200] + v[213] * v[915];
//		v[1631] = v[214] * v[926];
//		v[923] = v[1583] * v[199] + v[213] * v[913];
//		v[922] = 2e0*v[214] * v[923];
//		v[1047] = v[922];
//		v[924] = v[923];
//		v[1584] = v[1631] + v[215] * v[924];
//		v[925] = 2e0*v[215] * v[926];
//		v[1049] = v[925];
//		v[927] = -(v[157] * v[1584]);
//		v[928] = -(v[155] * v[1584]);
//		v[929] = -(v[153] * v[1584]);
//		v[930] = v[129] * v[1584];
//		v[931] = v[127] * v[1584];
//		v[932] = v[125] * v[1584];
//		v[933] = v[1585] * v[3698 + i717] + 2e0*v[785] * v[926] + v[924] * v[940];
//		v[941] = v[1585] * v[3734 + i717] + 2e0*v[793] * v[924] + v[926] * v[940];
//		v[942] = v[926];
//		v[984] = v[942];
//		v[943] = 2e0*v[216] * v[946];
//		v[1052] = v[943];
//		v[947] = v[1585] * v[2952 + i717] + v[945] * v[946];
//		v[949] = 2e0*v[772] * v[946] + (*a4)*v[948];
//		v[950] = v[946];
//		v[986] = v[950];
//		v[1633] = v[786] * v[986];
//		v[951] = 0e0;
//		v[952] = 0e0;
//		v[953] = 0e0;
//		v[954] = 0e0;
//		v[955] = 0e0;
//		v[956] = 0e0;
//		v[957] = 0e0;
//		v[958] = 0e0;
//		v[959] = 0e0;
//		v[960] = 0e0;
//		v[961] = 0e0;
//		v[962] = 0e0;
//		v[963] = 0e0;
//		v[964] = 0e0;
//		v[965] = 0e0;
//		v[966] = 0e0;
//		v[967] = 0e0;
//		v[968] = 0e0;
//		v[969] = 0e0;
//		v[970] = 0e0;
//		v[971] = 0e0;
//		v[972] = 0e0;
//		v[973] = 0e0;
//		v[974] = 0e0;
//		v[975] = 0e0;
//		v[976] = 0e0;
//		v[977] = 0e0;
//		v[978] = 0e0;
//		v[979] = 0e0;
//		v[980] = 0e0;
//		v[981] = 0e0;
//		v[982] = 0e0;
//		b983 = b411;
//		if (b983) {
//			v[985] = v[224] * v[950] - v[225] * v[984];
//			v[987] = v[225] * v[924] - v[223] * v[986];
//			v[988] = -(v[224] * v[924]) + v[223] * v[984];
//			v[1589] = v[859] * v[985] + v[854] * v[987] + v[849] * v[988];
//			v[1587] = v[413] * v[985] + v[414] * v[987] + v[415] * v[988];
//			v[989] = v[1587] / v[416];
//			v[1588] = v[1754] * v[989];
//			v[999] = v[1586] * v[989];
//			v[954] = -(v[1587] * v[862] * v[990]);
//			v[991] = v[1588] * v[413] + v[1540] * v[985];
//			v[1010] = 2e0*v[425] * v[991];
//			v[993] = v[1588] * v[414] + v[1540] * v[987];
//			v[1007] = 2e0*v[426] * v[993];
//			v[994] = v[1588] * v[415] + v[1540] * v[988];
//			v[1008] = 2e0*v[427] * v[994];
//			v[957] = v[1589] * v[418] + v[1561] * v[999];
//			v[956] = v[422] * v[989] * v[997];
//			v[955] = v[1589] * v[424] + v[1561] * v[423] * v[989];
//			v[981] = 2e0*v[1000] * v[1755] * v[999];
//			v[982] = v[1000] * v[1756] * v[861] * v[989] * v[995];
//			v[953] = v[1588] * v[849] + v[1590] * v[988];
//			v[952] = v[1588] * v[854] + v[1590] * v[987];
//			v[951] = v[1588] * v[859] + v[1590] * v[985];
//			v[1002] = (v[426] * v[991] + v[425] * v[993]) / 2e0;
//			v[1003] = v[1007] + v[1010];
//			v[1004] = v[1003] + v[1008];
//			v[1005] = (v[427] * v[991] + v[425] * v[994]) / 2e0;
//			v[1006] = (v[427] * v[993] + v[426] * v[994]) / 2e0;
//			v[1009] = v[1007] + v[1008];
//			v[1011] = v[1008] + v[1010];
//			v[960] = (v[834] * v[991] + v[833] * v[993] + 4e0*v[1012] * v[994]) / 2e0;
//			v[959] = (v[835] * v[991] + 4e0*v[1013] * v[993] + v[833] * v[994]) / 2e0;
//			v[958] = (4e0*v[1014] * v[991] + v[835] * v[993] + v[834] * v[994]) / 2e0;
//			v[1015] = -4e0*v[1004] * v[1268];
//			v[980] = 8e0*v[1004] * v[1757] * v[845];
//			v[971] = (v[1015] * v[795]) / 2e0;
//			v[975] = (v[1015] * v[799]) / 2e0;
//			v[973] = v[1015] * v[797];
//			v[977] = v[1015] * v[801];
//			v[976] = v[1015] * v[800];
//			v[978] = v[1015] * v[802];
//			v[972] = v[1015] * v[796];
//			v[974] = v[1015] * v[798];
//			v[979] = (v[1015] * v[803]) / 2e0;
//			v[1016] = v[1002] - v[994];
//			v[1017] = v[1002] + v[994];
//			v[1018] = v[1005] + v[993];
//			v[1019] = v[1005] - v[993];
//			v[1020] = v[1006] - v[991];
//			v[1021] = v[1006] + v[991];
//			v[963] = v[1016] * v[428] + v[1015] * v[837];
//			v[965] = v[1017] * v[428] + v[1015] * v[839];
//			v[964] = v[1018] * v[428] + v[1015] * v[838];
//			v[968] = v[1019] * v[428] + v[1015] * v[842];
//			v[962] = (-(v[1009] * v[428]) + v[1015] * v[836]) / 2e0;
//			v[967] = v[1020] * v[428] + v[1015] * v[841];
//			v[969] = v[1021] * v[428] + v[1015] * v[843];
//			v[966] = (-(v[1011] * v[428]) + v[1015] * v[840]) / 2e0;
//			v[970] = (-(v[1003] * v[428]) + v[1015] * v[844]) / 2e0;
//			v[961] = -(v[1009] * v[795]) / 2e0 + v[1016] * v[796] + v[1018] * v[797] + v[1017] * v[798] - (v[1011] * v[799]) / 2e0
//				+ v[1020] * v[800] + v[1019] * v[801] + v[1021] * v[802] - (v[1003] * v[803]) / 2e0;
//		}
//		else {
//		};
//		v[1046] = v[950];
//		v[1040] = v[942];
//		v[1038] = v[924];
//		v[1036] = v[963];
//		v[1035] = v[964];
//		v[1033] = v[966];
//		v[1032] = v[967];
//		v[1030] = v[969];
//		v[1029] = v[970];
//		v[1022] = 0e0;
//		v[1023] = 0e0;
//		v[1024] = 0e0;
//		v[1025] = 0e0;
//		v[1026] = 0e0;
//		v[1027] = 0e0;
//		b1028 = (*previouscontact);
//		if (b1028) {
//			v[1044] = v[1038] * v[747];
//			v[1043] = v[745] * v[950];
//			v[1041] = v[746] * v[942];
//			v[970] = 0e0;
//			v[969] = 0e0;
//			v[1031] = gti[2] * v[1029] + gti[1] * v[1030] + i1580 * v[130] + i1591 * v[132] + i1592 * v[134] - i1581 * v[158]
//				- i1593 * v[160] - i1594 * v[162] + gti[0] * v[968];
//			v[1603] = -(v[1031] * v[216]) - v[1046] * v[448];
//			v[968] = 0e0;
//			v[967] = 0e0;
//			v[966] = 0e0;
//			v[1034] = gti[2] * v[1032] + gti[1] * v[1033] + i1578 * v[130] + i1595 * v[132] + i1596 * v[134] - i1579 * v[158]
//				- i1597 * v[160] - i1598 * v[162] + gti[0] * v[965];
//			v[1605] = -(v[1034] * v[215]) - v[1040] * v[447];
//			v[965] = 0e0;
//			v[964] = 0e0;
//			v[963] = 0e0;
//			v[1037] = gti[2] * v[1035] + gti[1] * v[1036] + i1576 * v[130] + i1599 * v[132] + i1600 * v[134] - i1577 * v[158]
//				- i1601 * v[160] - i1602 * v[162] + gti[0] * v[962];
//			v[1604] = -(v[1037] * v[214]) - v[1038] * v[446];
//			v[962] = 0e0;
//			v[1039] = v[1041] + v[1044];
//			v[1042] = v[1041] + v[1043];
//			v[1045] = v[1043] + v[1044];
//			v[1027] = v[1031] * v[745];
//			v[1026] = v[1034] * v[746];
//			v[1025] = v[1037] * v[747];
//			v[1022] = v[1038] * v[1048] - v[1042] * v[214] - v[1047] * v[747];
//			v[941] = v[1037] * v[1048] - v[1042] * v[446] + (v[1603] + v[1605])*v[747] + v[941];
//			v[1023] = v[1040] * v[1051] - v[1045] * v[215] - v[1049] * v[746];
//			v[933] = v[1034] * v[1051] - v[1045] * v[447] + (v[1603] + v[1604])*v[746] + v[933];
//			v[1024] = v[1046] * v[1055] - v[1039] * v[216] - v[1052] * v[745];
//			v[949] = v[1031] * v[1055] - v[1039] * v[448] + (v[1604] + v[1605])*v[745] + v[949];
//		}
//		else {
//		};
//		v[1632] = v[750] * v[943];
//		v[1063] = v[950];
//		v[1637] = v[749] * v[925];
//		v[1061] = v[942];
//		v[1306] = v[1061] * v[1606];
//		v[1056] = v[924];
//		v[1305] = v[1056] * v[1606];
//		v[1058] = (*a4)*v[3536 + i717] + v[1056] * v[794] + v[1057] * v[922] + v[197] * v[927] + v[194] * v[928] + v[191] * v[929]
//			+ v[188] * v[930] + v[185] * v[931] + v[182] * v[932] + v[214] * v[947];
//		v[1085] = v[1058];
//		v[933] = v[933] + v[749] * v[947];
//		v[1226] = v[933];
//		v[941] = v[941] + v[748] * v[947];
//		v[1227] = v[941];
//		v[1059] = v[1056] * v[748] + v[749] * v[942];
//		v[949] = v[750] * (v[1056] * v[1057] + v[1060] * v[942]) + v[949];
//		v[1225] = v[949];
//		v[1062] = (*a4)*v[3518 + i717] + v[1061] * v[794] + v[1060] * v[925] + v[196] * v[927] + v[193] * v[928] + v[190] * v[929]
//			+ v[187] * v[930] + v[184] * v[931] + v[181] * v[932] + v[215] * v[947];
//		v[1086] = v[1062];
//		v[1067] = v[1056] * v[1064] + v[1061] * v[1065] + v[1063] * v[1066] + (*a4)*(v[3482 + i717] + v[216] * v[3500 + i717])
//			+ v[943] * v[945];
//		v[1068] = 0e0;
//		b1069 = b478;
//		if (b1069) {
//			b1070 = b480;
//			if (b1070) {
//				v[1068] = -v[1071];
//				v[924] = 0e0;
//				v[942] = 0e0;
//				v[950] = 0e0;
//			}
//			else {
//			};
//		}
//		else {
//		};
//		v[1072] = 0e0;
//		v[1073] = 0e0;
//		v[1074] = 0e0;
//		v[1075] = 0e0;
//		v[1076] = 0e0;
//		v[1077] = 0e0;
//		v[1078] = 0e0;
//		v[1079] = 0e0;
//		b1080 = b478;
//		if (b1080) {
//			v[1090] = v[1085] * v[754] + v[1086] * v[755] + v[1067] * v[756];
//			b1081 = b480;
//			if (b1081) {
//				v[1084] = v[1068] * v[1212] * (*zetan);
//				v[1083] = v[1084] * v[1213];
//				v[1078] = -((v[1083] * v[766]) / v[760]);
//				v[1075] = v[1082] * v[1084] * v[1789] * v[766];
//				v[1068] = 0e0;
//				v[1058] = 0e0;
//				v[1062] = 0e0;
//				v[1076] = v[1090];
//				v[1067] = 0e0;
//			}
//			else {
//				v[1088] = v[1218] * v[921] * (*zetan);
//				v[1083] = v[1088] * v[1219];
//				v[1079] = -((v[1083] * v[766]) / v[765]);
//				v[920] = v[1087] * v[1088] * v[1794] * v[766] + v[920];
//				v[921] = 0e0;
//				v[1058] = 0e0;
//				v[1062] = 0e0;
//				v[1077] = v[1090];
//				v[1067] = 0e0;
//			};
//			v[1074] = v[1083] * v[756];
//			v[1073] = v[1083] * v[755];
//			v[1072] = v[1083] * v[754];
//		}
//		else {
//		};
//		v[1217] = v[1072];
//		v[1215] = v[1073];
//		v[1214] = v[1074];
//		v[1091] = (*ct)*((i1599 / 2e0 + v[1092])*v[719] + (i1595 / 2e0 + v[1093])*v[720] + (i1591 / 2e0 + v[1094])*v[721]);
//		v[1095] = (*ct)*((i1600 / 2e0 + v[1092])*v[719] + (i1596 / 2e0 + v[1093])*v[720] + (i1592 / 2e0 + v[1094])*v[721]);
//		v[1096] = (*ct)*((i1601 / 2e0 + v[1097])*v[719] + (i1597 / 2e0 + v[1098])*v[720] + (i1593 / 2e0 + v[1099])*v[721]);
//		v[1100] = (*ct)*((i1602 / 2e0 + v[1097])*v[719] + (i1598 / 2e0 + v[1098])*v[720] + (i1594 / 2e0 + v[1099])*v[721]);
//		v[1101] = -(i1576*v[644]) - i1578 * v[648] - i1580 * v[652] - i1599 * v[656] - i1595 * v[660] - i1591 * v[664]
//			- i1600 * v[668] - i1596 * v[672] - i1592 * v[676] - i1577 * v[680] - i1579 * v[684] - i1581 * v[688] - i1601 * v[692]
//			- i1597 * v[696] - i1593 * v[700] - i1602 * v[704] - i1598 * v[708] - i1594 * v[712];
//		v[1117] = v[1101];
//		v[1102] = -(i1576*v[645]) - i1578 * v[649] - i1580 * v[653] - i1599 * v[657] - i1595 * v[661] - i1591 * v[665]
//			- i1600 * v[669] - i1596 * v[673] - i1592 * v[677] - i1577 * v[681] - i1579 * v[685] - i1581 * v[689] - i1601 * v[693]
//			- i1597 * v[697] - i1593 * v[701] - i1602 * v[705] - i1598 * v[709] - i1594 * v[713];
//		v[1115] = v[1102];
//		v[1103] = -(i1576*v[646]) - i1578 * v[650] - i1580 * v[654] - i1599 * v[658] - i1595 * v[662] - i1591 * v[666]
//			- i1600 * v[670] - i1596 * v[674] - i1592 * v[678] - i1577 * v[682] - i1579 * v[686] - i1581 * v[690] - i1601 * v[694]
//			- i1597 * v[698] - i1593 * v[702] - i1602 * v[706] - i1598 * v[710] - i1594 * v[714];
//		v[1113] = v[1103];
//		v[1104] = 0e0;
//		v[1105] = 0e0;
//		v[1106] = 0e0;
//		v[1107] = 0e0;
//		v[1108] = 0e0;
//		v[1109] = 0e0;
//		b1110 = b478;
//		if (b1110) {
//			b1111 = (*stick);
//			if (b1111) {
//				b1112 = b516;
//				if (b1112) {
//					v[1109] = v[1103];
//					v[1103] = 0e0;
//					v[1108] = v[1102];
//					v[1102] = 0e0;
//					v[1107] = v[1101];
//					v[1101] = 0e0;
//				}
//				else {
//					v[1607] = (*mud)*(v[1117] * v[511] + v[1115] * v[512] + v[1113] * v[513]);
//					v[1103] = 0e0;
//					v[1102] = 0e0;
//					v[1609] = v[1128] * v[1607] * v[527];
//					v[1101] = 0e0;
//					v[1608] = v[1126] * v[1607] * v[526] * v[531];
//					v[1109] = v[1113] * v[1125] + v[1608] * v[513];
//					v[1108] = v[1115] * v[1125] + v[1608] * v[512];
//					v[1107] = v[1117] * v[1125] + v[1608] * v[511];
//					v[1106] = v[1609] * v[504];
//					v[1105] = v[1609] * v[503];
//					v[1104] = v[1609] * v[502];
//				};
//			}
//			else {
//				b1129 = b550;
//				if (b1129) {
//					v[1109] = v[1113];
//					v[1103] = 0e0;
//					v[1108] = v[1115];
//					v[1102] = 0e0;
//					v[1107] = v[1117];
//					v[1101] = 0e0;
//				}
//				else {
//					v[1135] = v[1117] * v[1610];
//					v[1133] = v[1115] * v[1610];
//					v[1132] = v[1113] * v[1610];
//					v[1103] = 0e0;
//					v[1102] = 0e0;
//					v[1612] = (*mud)*v[1138] * (v[1117] * v[511] + v[1115] * v[512] + v[1113] * v[513])*v[557];
//					v[1101] = 0e0;
//					v[1611] = v[1134] * (v[1135] * v[511] + v[1133] * v[512] + v[1132] * v[513])*v[556];
//					v[1109] = v[1611] * v[513] + v[1132] * v[557];
//					v[1108] = v[1611] * v[512] + v[1133] * v[557];
//					v[1107] = v[1611] * v[511] + v[1135] * v[557];
//					v[1106] = v[1612] * v[504];
//					v[1105] = v[1612] * v[503];
//					v[1104] = v[1612] * v[502];
//				};
//			};
//		}
//		else {
//		};
//		v[1615] = -((*ct)*v[1107]);
//		v[1614] = -((*ct)*v[1108]);
//		v[1613] = -((*ct)*v[1109]);
//		v[1140] = v[1613] * v[198] - i1594 * v[520];
//		v[1141] = v[1613] * v[197] - i1598 * v[520];
//		v[1142] = v[1613] * v[196] - i1602 * v[520];
//		v[1143] = v[1613] * v[195] - i1593 * v[520];
//		v[1144] = v[1613] * v[194] - i1597 * v[520];
//		v[1145] = v[1613] * v[193] - i1601 * v[520];
//		v[1146] = v[1613] * v[192] - i1581 * v[520];
//		v[1147] = v[1613] * v[191] - i1579 * v[520];
//		v[1148] = v[1613] * v[190] - i1577 * v[520];
//		v[1149] = v[1613] * v[189] - i1592 * v[520];
//		v[1150] = v[1613] * v[188] - i1596 * v[520];
//		v[1151] = v[1613] * v[187] - i1600 * v[520];
//		v[1152] = v[1613] * v[186] - i1591 * v[520];
//		v[1153] = v[1613] * v[185] - i1595 * v[520];
//		v[1154] = v[1613] * v[184] - i1599 * v[520];
//		v[1155] = v[1613] * v[183] - i1580 * v[520];
//		v[1156] = v[1613] * v[182] - i1578 * v[520];
//		v[1157] = v[1613] * v[181] - i1576 * v[520];
//		v[1158] = v[1614] * v[198] - i1594 * v[519];
//		v[1159] = v[1614] * v[197] - i1598 * v[519];
//		v[1160] = v[1614] * v[196] - i1602 * v[519];
//		v[1161] = v[1614] * v[195] - i1593 * v[519];
//		v[1162] = v[1614] * v[194] - i1597 * v[519];
//		v[1163] = v[1614] * v[193] - i1601 * v[519];
//		v[1164] = v[1614] * v[192] - i1581 * v[519];
//		v[1165] = v[1614] * v[191] - i1579 * v[519];
//		v[1166] = v[1614] * v[190] - i1577 * v[519];
//		v[1167] = v[1614] * v[189] - i1592 * v[519];
//		v[1168] = v[1614] * v[188] - i1596 * v[519];
//		v[1169] = v[1614] * v[187] - i1600 * v[519];
//		v[1170] = v[1614] * v[186] - i1591 * v[519];
//		v[1171] = v[1614] * v[185] - i1595 * v[519];
//		v[1172] = v[1614] * v[184] - i1599 * v[519];
//		v[1173] = v[1614] * v[183] - i1580 * v[519];
//		v[1174] = v[1614] * v[182] - i1578 * v[519];
//		v[1175] = v[1614] * v[181] - i1576 * v[519];
//		v[1176] = v[1615] * v[198] - i1594 * v[518];
//		v[1177] = v[1615] * v[197] - i1598 * v[518];
//		v[1178] = v[1615] * v[196] - i1602 * v[518];
//		v[1179] = v[1615] * v[195] - i1593 * v[518];
//		v[1180] = v[1615] * v[194] - i1597 * v[518];
//		v[1181] = v[1615] * v[193] - i1601 * v[518];
//		v[1182] = v[1615] * v[192] - i1581 * v[518];
//		v[1183] = v[1615] * v[191] - i1579 * v[518];
//		v[1184] = v[1615] * v[190] - i1577 * v[518];
//		v[1185] = v[1615] * v[189] - i1592 * v[518];
//		v[1186] = v[1615] * v[188] - i1596 * v[518];
//		v[1187] = v[1615] * v[187] - i1600 * v[518];
//		v[1188] = v[1615] * v[186] - i1591 * v[518];
//		v[1189] = v[1615] * v[185] - i1595 * v[518];
//		v[1190] = v[1615] * v[184] - i1599 * v[518];
//		v[1191] = v[1615] * v[183] - i1580 * v[518];
//		v[1192] = v[1615] * v[182] - i1578 * v[518];
//		v[1193] = v[1615] * v[181] - i1576 * v[518];
//		v[1194] = (*epst)*v[1109];
//		v[1195] = (*epst)*v[1108];
//		v[1196] = (*epst)*v[1107];
//		v[1197] = v[1106];
//		v[1198] = v[1106];
//		v[1199] = v[1105];
//		v[1200] = v[1105];
//		v[1201] = v[1104];
//		v[1202] = v[1104];
//		b1203 = b478;
//		if (b1203) {
//			v[1204] = 0e0;
//			v[1205] = 0e0;
//			v[1206] = 0e0;
//			b1207 = b497;
//			if (b1207) {
//				v[1206] = v[1197];
//				v[1197] = 0e0;
//				v[1205] = v[1199];
//				v[1199] = 0e0;
//				v[1204] = v[1201];
//				v[1201] = 0e0;
//			}
//			else {
//				v[1198] = 0e0;
//				v[1197] = 0e0;
//				v[1200] = 0e0;
//				v[1199] = 0e0;
//				v[1202] = 0e0;
//				v[1201] = 0e0;
//			};
//			v[1216] = v[1204] * v[452] + v[1205] * v[466] + v[1206] * v[468];
//			b1211 = b480;
//			if (b1211) {
//				v[1074] = v[1074] + v[1206] * v[492];
//				v[1073] = v[1073] + v[1205] * v[492];
//				v[1076] = v[1076] + v[1216];
//				v[1072] = v[1072] + v[1204] * v[492];
//				v[1078] = v[1078] + 2e0*v[1076] * (*zetan);
//				v[1075] = v[1075] + (v[1078] * v[1212] * v[1213]) / 2e0;
//			}
//			else {
//				v[1074] = v[1214] + v[1206] * v[496];
//				v[1073] = v[1215] + v[1205] * v[496];
//				v[1077] = v[1077] + v[1216];
//				v[1072] = v[1217] + v[1204] * v[496];
//				v[1079] = v[1079] + 2e0*v[1077] * (*zetan);
//				v[920] = (v[1079] * v[1218] * v[1219]) / 2e0 + v[920];
//			};
//		}
//		else {
//		};
//		v[1638] = v[1073] * v[459];
//		v[1640] = v[1306] + v[1637] + v[1638];
//		v[1634] = v[1074] * v[467];
//		v[1657] = v[1633] + v[1634];
//		v[1639] = v[1632] + v[1633] + v[1634];
//		v[1228] = v[920];
//		v[1224] = v[1202];
//		v[1223] = v[1200];
//		v[1222] = v[1198];
//		b1220 = b478;
//		if (b1220) {
//			v[1617] = -(v[1224] * v[214]) - v[1223] * v[215] - v[1222] * v[216];
//			b1221 = b480;
//			if (b1221) {
//				v[949] = v[1198] * v[1616] + v[949];
//				v[1198] = 0e0;
//				v[933] = v[1200] * v[1616] + v[933];
//				v[1200] = 0e0;
//				v[941] = v[1202] * v[1616] + v[941];
//				v[1202] = 0e0;
//				v[1075] = v[1075] + v[1545] * v[1617] * v[1832];
//				v[920] = -v[1075] + v[920];
//			}
//			else {
//				v[949] = v[1225] + v[1222] * v[486];
//				v[1198] = 0e0;
//				v[933] = v[1226] + v[1223] * v[486];
//				v[1200] = 0e0;
//				v[941] = v[1227] + v[1224] * v[486];
//				v[1202] = 0e0;
//				v[920] = v[1228] + (*n2)*v[1617] * v[1836] * v[474];
//			};
//		}
//		else {
//		};
//		v[1229] = (*a4)*v[3554 + i717] + v[1074] * v[945];
//		v[949] = v[1066] * v[1074] + v[949];
//		v[1230] = v[1074] * v[216] + v[1063] * v[750];
//		v[1642] = v[1306] + v[1637] + v[1638] + v[1230] * v[215];
//		v[1656] = -v[1306] + v[1642];
//		v[1636] = v[1305] + v[1230] * v[214] + v[1072] * v[451] + v[748] * v[922];
//		v[1655] = -v[1305] + v[1636];
//		v[3987] = -(v[125] * v[1655]) - v[1073] * v[453] - v[749] * v[932];
//		v[3988] = -(v[125] * v[1656]) - v[1072] * v[453] - v[748] * v[932];
//		v[3989] = -(v[125] * v[1657]) - v[934] * v[943];
//		v[3990] = -(v[127] * v[1655]) - v[1073] * v[454] - v[749] * v[931];
//		v[3991] = -(v[127] * v[1656]) - v[1072] * v[454] - v[748] * v[931];
//		v[3992] = -(v[127] * v[1657]) - v[935] * v[943];
//		v[3993] = -(v[129] * v[1655]) - v[1073] * v[455] - v[749] * v[930];
//		v[3994] = -(v[129] * v[1656]) - v[1072] * v[455] - v[748] * v[930];
//		v[3995] = -(v[129] * v[1657]) - v[936] * v[943];
//		v[3996] = v[153] * v[1655] - v[1073] * v[456] - v[749] * v[929];
//		v[3997] = v[153] * v[1656] - v[1072] * v[456] - v[748] * v[929];
//		v[3998] = v[153] * v[1657] - v[937] * v[943];
//		v[3999] = v[155] * v[1655] - v[1073] * v[457] - v[749] * v[928];
//		v[4000] = v[155] * v[1656] - v[1072] * v[457] - v[748] * v[928];
//		v[4001] = v[155] * v[1657] - v[938] * v[943];
//		v[4002] = v[157] * v[1655] - v[1073] * v[458] - v[749] * v[927];
//		v[4003] = v[157] * v[1656] - v[1072] * v[458] - v[748] * v[927];
//		v[4004] = v[157] * v[1657] - v[939] * v[943];
//		v[1231] = v[1060] * v[1073] + (*a4)*v[3680 + i717];
//		v[933] = v[1073] * v[794] + v[933];
//		v[1232] = v[1073] * v[181] + v[1072] * v[182] + (*a4)*v[3572 + i717];
//		v[1233] = v[1073] * v[184] + v[1072] * v[185] + (*a4)*v[3590 + i717];
//		v[1234] = v[1073] * v[187] + v[1072] * v[188] + (*a4)*v[3608 + i717];
//		v[1235] = v[1073] * v[190] + v[1072] * v[191] + (*a4)*v[3626 + i717];
//		v[1236] = v[1073] * v[193] + v[1072] * v[194] + (*a4)*v[3644 + i717];
//		v[1237] = v[1073] * v[196] + v[1072] * v[197] + (*a4)*v[3662 + i717];
//		v[1630] = v[1232] * v[125] + v[1233] * v[127] + v[1234] * v[129] - v[1235] * v[153] - v[1236] * v[155] - v[1237] * v[157];
//		v[1238] = v[1057] * v[1072] + (*a4)*v[3716 + i717];
//		v[941] = v[1072] * v[794] + v[941];
//		v[1239] = 0e0;
//		v[1240] = 0e0;
//		v[1241] = 0e0;
//		v[1242] = 0e0;
//		v[1243] = 0e0;
//		v[1244] = 0e0;
//		v[1245] = 0e0;
//		v[1246] = 0e0;
//		v[1247] = 0e0;
//		b1248 = (*previouscontact);
//		if (b1248) {
//			v[1620] = v[1196] * v[214];
//			v[1619] = v[1195] * v[215];
//			v[1621] = v[1619] + v[1620];
//			v[1618] = v[1194] * v[216];
//			v[1623] = v[1618] + v[1619];
//			v[1622] = v[1618] + v[1620];
//			v[1027] = v[1027] + v[1194] * v[448];
//			v[1026] = v[1026] + v[1195] * v[447];
//			v[1022] = v[1022] - v[1623] * v[214] + v[1196] * v[805];
//			v[1023] = v[1023] - v[1622] * v[215] + v[1195] * v[807];
//			v[1024] = v[1024] - v[1621] * v[216] + v[1194] * v[809];
//			v[1025] = v[1025] + v[1196] * v[446];
//			v[949] = v[1194] * v[1558] - v[1621] * v[448] + v[949];
//			v[933] = v[1195] * v[1559] - v[1622] * v[447] + v[933];
//			v[941] = -(v[1196] * v[1560]) - v[1623] * v[446] + v[941];
//			v[1239] = gti[0] * v[1022];
//			v[1240] = gti[1] * v[1022];
//			v[1241] = gti[2] * v[1022];
//			v[1249] = -(v[1022] * v[162]);
//			v[1250] = -(v[1022] * v[160]);
//			v[1251] = -(v[1022] * v[158]);
//			v[1252] = v[1022] * v[134];
//			v[1253] = v[1022] * v[132];
//			v[1254] = v[1022] * v[130];
//			v[1242] = gti[0] * v[1023];
//			v[1243] = gti[1] * v[1023];
//			v[1244] = gti[2] * v[1023];
//			v[1255] = -(v[1023] * v[162]);
//			v[1256] = -(v[1023] * v[160]);
//			v[1257] = -(v[1023] * v[158]);
//			v[1258] = v[1023] * v[134];
//			v[1259] = v[1023] * v[132];
//			v[1260] = v[1023] * v[130];
//			v[1245] = gti[0] * v[1024];
//			v[1246] = gti[1] * v[1024];
//			v[1247] = gti[2] * v[1024];
//			v[1261] = -(v[1024] * v[162]);
//			v[1262] = -(v[1024] * v[160]);
//			v[1263] = -(v[1024] * v[158]);
//			v[1264] = v[1024] * v[134];
//			v[1265] = v[1024] * v[132];
//			v[1266] = v[1024] * v[130];
//			v[1238] = -v[1025] + v[1238];
//			v[1231] = -v[1026] + v[1231];
//			v[1229] = -v[1027] + v[1229];
//		}
//		else {
//			v[1254] = 0e0;
//			v[1253] = 0e0;
//			v[1252] = 0e0;
//			v[1260] = 0e0;
//			v[1259] = 0e0;
//			v[1258] = 0e0;
//			v[1266] = 0e0;
//			v[1265] = 0e0;
//			v[1264] = 0e0;
//			v[1251] = 0e0;
//			v[1250] = 0e0;
//			v[1249] = 0e0;
//			v[1257] = 0e0;
//			v[1256] = 0e0;
//			v[1255] = 0e0;
//			v[1263] = 0e0;
//			v[1262] = 0e0;
//			v[1261] = 0e0;
//		};
//		b1267 = b411;
//		if (b1267) {
//			v[979] = (v[1247] * v[428]) / 2e0 + v[979];
//			v[978] = v[1246] * v[428] + v[978];
//			v[977] = v[1245] * v[428] + v[977];
//			v[976] = v[1244] * v[428] + v[976];
//			v[975] = (v[1243] * v[428]) / 2e0 + v[975];
//			v[974] = v[1242] * v[428] + v[974];
//			v[973] = v[1241] * v[428] + v[973];
//			v[972] = v[1240] * v[428] + v[972];
//			v[961] = (v[1239] * v[836]) / 2e0 + v[1240] * v[837] + v[1241] * v[838] + v[1242] * v[839] + (v[1243] * v[840]) / 2e0
//				+ v[1244] * v[841] + v[1245] * v[842] + v[1246] * v[843] + (v[1247] * v[844]) / 2e0 + v[961];
//			v[971] = (v[1239] * v[428]) / 2e0 + v[971];
//			v[980] = -4e0*v[1268] * v[961] + v[980];
//			v[1624] = v[971] - v[980];
//			v[1626] = (v[973] + v[977]) / 2e0;
//			v[1625] = (v[976] + v[978]) / 2e0;
//			v[960] = v[1626] * v[425] + v[1625] * v[426] + v[960] - v[972] + v[974] - 2e0*v[427] * (v[1624] + v[975]);
//			v[1627] = (v[972] + v[974]) / 2e0;
//			v[959] = v[1627] * v[425] + v[1625] * v[427] + v[959] + v[973] - v[977] - 2e0*v[426] * (v[1624] + v[979]);
//			v[958] = v[1627] * v[426] + v[1626] * v[427] + v[958] - v[976] + v[978] - 2e0*v[425] * (v[975] + v[979] - v[980]);
//			v[1628] = v[413] * v[958] + v[414] * v[959] + v[415] * v[960];
//			v[957] = v[1628] * v[418] + v[957];
//			v[955] = v[1628] * v[424] + v[955];
//			v[956] = v[956] + v[423] * v[957];
//			v[981] = v[981] + 2e0*v[955] * v[995];
//			v[982] = v[982] + (v[981] * v[996]) / 2e0;
//			v[954] = v[954] + v[956] + v[982];
//			v[1629] = v[954] / v[416];
//			v[953] = v[1629] * v[415] + v[953] + v[1540] * v[960];
//			v[952] = v[1629] * v[414] + v[952] + v[1540] * v[959];
//			v[951] = v[1629] * v[413] + v[951] + v[1540] * v[958];
//			v[941] = v[941] + v[225] * v[952] - v[224] * v[953];
//			v[949] = v[949] + v[224] * v[951] - v[223] * v[952];
//			v[933] = v[933] - v[225] * v[951] + v[223] * v[953];
//		}
//		else {
//		};
//		v[1647] = (-(v[1193] * v[393]) - v[1192] * v[394] - v[1191] * v[395] - v[1190] * v[396] - v[1189] * v[397]
//			- v[1188] * v[398] - v[1187] * v[399] - v[1186] * v[400] - v[1185] * v[401] - v[1184] * v[402] - v[1183] * v[403]
//			- v[1182] * v[404] - v[1181] * v[405] - v[1180] * v[406] - v[1179] * v[407] - v[1178] * v[408] - v[1177] * v[409]
//			- v[1176] * v[410] + v[1276] * v[875]) / 2e0;
//		v[1648] = (-(v[1193] * v[375]) - v[1192] * v[376] - v[1191] * v[377] - v[1190] * v[378] - v[1189] * v[379]
//			- v[1188] * v[380] - v[1187] * v[381] - v[1186] * v[382] - v[1185] * v[383] - v[1184] * v[384] - v[1183] * v[385]
//			- v[1182] * v[386] - v[1181] * v[387] - v[1180] * v[388] - v[1179] * v[389] - v[1178] * v[390] - v[1177] * v[391]
//			- v[1176] * v[392] + v[1276] * v[874]) / 2e0;
//		v[1653] = (v[1193] * v[357] + v[1192] * v[358] + v[1191] * v[359] + v[1190] * v[360] + v[1189] * v[361] + v[1188] * v[362]
//			+ v[1187] * v[363] + v[1186] * v[364] + v[1185] * v[365] + v[1184] * v[366] + v[1183] * v[367] + v[1182] * v[368]
//			+ v[1181] * v[369] + v[1180] * v[370] + v[1179] * v[371] + v[1178] * v[372] + v[1177] * v[373] + v[1176] * v[374]
//			- v[1276] * v[873]) / 2e0;
//		v[1654] = (v[1193] * v[339] + v[1192] * v[340] + v[1191] * v[341] + v[1190] * v[342] + v[1189] * v[343] + v[1188] * v[344]
//			+ v[1187] * v[345] + v[1186] * v[346] + v[1185] * v[347] + v[1184] * v[348] + v[1183] * v[349] + v[1182] * v[350]
//			+ v[1181] * v[351] + v[1180] * v[352] + v[1179] * v[353] + v[1178] * v[354] + v[1177] * v[355] + v[1176] * v[356]
//			- v[1276] * v[872]) / 2e0;
//		v[1645] = (-(v[1175] * v[393]) - v[1174] * v[394] - v[1173] * v[395] - v[1172] * v[396] - v[1171] * v[397]
//			- v[1170] * v[398] - v[1169] * v[399] - v[1168] * v[400] - v[1167] * v[401] - v[1166] * v[402] - v[1165] * v[403]
//			- v[1164] * v[404] - v[1163] * v[405] - v[1162] * v[406] - v[1161] * v[407] - v[1160] * v[408] - v[1159] * v[409]
//			- v[1158] * v[410] + v[1281] * v[875]) / 2e0;
//		v[1646] = (-(v[1175] * v[375]) - v[1174] * v[376] - v[1173] * v[377] - v[1172] * v[378] - v[1171] * v[379]
//			- v[1170] * v[380] - v[1169] * v[381] - v[1168] * v[382] - v[1167] * v[383] - v[1166] * v[384] - v[1165] * v[385]
//			- v[1164] * v[386] - v[1163] * v[387] - v[1162] * v[388] - v[1161] * v[389] - v[1160] * v[390] - v[1159] * v[391]
//			- v[1158] * v[392] + v[1281] * v[874]) / 2e0;
//		v[1651] = (v[1175] * v[357] + v[1174] * v[358] + v[1173] * v[359] + v[1172] * v[360] + v[1171] * v[361] + v[1170] * v[362]
//			+ v[1169] * v[363] + v[1168] * v[364] + v[1167] * v[365] + v[1166] * v[366] + v[1165] * v[367] + v[1164] * v[368]
//			+ v[1163] * v[369] + v[1162] * v[370] + v[1161] * v[371] + v[1160] * v[372] + v[1159] * v[373] + v[1158] * v[374]
//			- v[1281] * v[873]) / 2e0;
//		v[1652] = (v[1175] * v[339] + v[1174] * v[340] + v[1173] * v[341] + v[1172] * v[342] + v[1171] * v[343] + v[1170] * v[344]
//			+ v[1169] * v[345] + v[1168] * v[346] + v[1167] * v[347] + v[1166] * v[348] + v[1165] * v[349] + v[1164] * v[350]
//			+ v[1163] * v[351] + v[1162] * v[352] + v[1161] * v[353] + v[1160] * v[354] + v[1159] * v[355] + v[1158] * v[356]
//			- v[1281] * v[872]) / 2e0;
//		v[1643] = (-(v[1157] * v[393]) - v[1156] * v[394] - v[1155] * v[395] - v[1154] * v[396] - v[1153] * v[397]
//			- v[1152] * v[398] - v[1151] * v[399] - v[1150] * v[400] - v[1149] * v[401] - v[1148] * v[402] - v[1147] * v[403]
//			- v[1146] * v[404] - v[1145] * v[405] - v[1144] * v[406] - v[1143] * v[407] - v[1142] * v[408] - v[1141] * v[409]
//			- v[1140] * v[410] + v[1286] * v[875]) / 2e0;
//		v[1644] = (-(v[1157] * v[375]) - v[1156] * v[376] - v[1155] * v[377] - v[1154] * v[378] - v[1153] * v[379]
//			- v[1152] * v[380] - v[1151] * v[381] - v[1150] * v[382] - v[1149] * v[383] - v[1148] * v[384] - v[1147] * v[385]
//			- v[1146] * v[386] - v[1145] * v[387] - v[1144] * v[388] - v[1143] * v[389] - v[1142] * v[390] - v[1141] * v[391]
//			- v[1140] * v[392] + v[1286] * v[874]) / 2e0;
//		v[1649] = (v[1157] * v[357] + v[1156] * v[358] + v[1155] * v[359] + v[1154] * v[360] + v[1153] * v[361] + v[1152] * v[362]
//			+ v[1151] * v[363] + v[1150] * v[364] + v[1149] * v[365] + v[1148] * v[366] + v[1147] * v[367] + v[1146] * v[368]
//			+ v[1145] * v[369] + v[1144] * v[370] + v[1143] * v[371] + v[1142] * v[372] + v[1141] * v[373] + v[1140] * v[374]
//			- v[1286] * v[873]) / 2e0;
//		v[1650] = (v[1157] * v[339] + v[1156] * v[340] + v[1155] * v[341] + v[1154] * v[342] + v[1153] * v[343] + v[1152] * v[344]
//			+ v[1151] * v[345] + v[1150] * v[346] + v[1149] * v[347] + v[1148] * v[348] + v[1147] * v[349] + v[1146] * v[350]
//			+ v[1145] * v[351] + v[1144] * v[352] + v[1143] * v[353] + v[1142] * v[354] + v[1141] * v[355] + v[1140] * v[356]
//			- v[1286] * v[872]) / 2e0;
//		v[1290] = v[1059] + v[1072] * v[214] + v[1073] * v[215];
//		v[4059] = -(v[1056] * v[934]);
//		v[4060] = -(v[1061] * v[934]);
//		v[4061] = -(v[125] * v[1290]);
//		v[4062] = -(v[1056] * v[935]);
//		v[4063] = -(v[1061] * v[935]);
//		v[4064] = -(v[127] * v[1290]);
//		v[4065] = -(v[1056] * v[936]);
//		v[4066] = -(v[1061] * v[936]);
//		v[4067] = -(v[129] * v[1290]);
//		v[4068] = -(v[1056] * v[937]);
//		v[4069] = -(v[1061] * v[937]);
//		v[4070] = v[1290] * v[153];
//		v[4071] = -(v[1056] * v[938]);
//		v[4072] = -(v[1061] * v[938]);
//		v[4073] = v[1290] * v[155];
//		v[4074] = -(v[1056] * v[939]);
//		v[4075] = -(v[1061] * v[939]);
//		v[4076] = v[1290] * v[157];
//		v[1635] = v[1290] * v[216];
//		v[1641] = v[1632] + v[1635] + v[1657];
//		v[949] = 2e0*v[1229] * v[216] + v[1290] * v[945] + v[949];
//		v[933] = v[1060] * v[1230] + v[1630] * v[214] + 2e0*v[1231] * v[215] + v[933];
//		v[941] = v[1057] * v[1230] + 2e0*v[1238] * v[214] + v[1630] * v[215] + v[941];
//		v[920] = v[1071] * v[1291] * v[211] + v[920] + v[212] * (v[767] * v[913] + v[768] * v[915] + v[769] * v[917]
//			+ v[200] * v[933] + v[199] * v[941] + v[201] * v[949]);
//		v[1293] = v[1583] * v[769] + (v[751] * v[917] + v[201] * v[920]) / v[477] + v[213] * v[949];
//		v[1295] = v[1583] * v[768] + (v[751] * v[915] + v[200] * v[920]) / v[477] + v[213] * v[933];
//		v[1297] = v[1583] * v[767] + (v[751] * v[913] + v[199] * v[920]) / v[477] + v[213] * v[941];
//		v[1299] = -(v[1297] * v[169]) - v[1295] * v[173] - v[1293] * v[177] - v[1636] * v[190] - v[1640] * v[191]
//			- v[1635] * v[192] - v[1639] * v[192] + (*a4)*(v[3824 + i717] + v[216] * v[3842 + i717]) - v[1631] * v[790]
//			- i1581 * v[867] - i1579 * v[868] - i1577 * v[869] - v[215] * (v[1230] * v[191] + v[1235] * v[214] + v[790] * v[923]);
//		v[1307] = -2e0*v[1232] * v[1380] + v[1297] * v[141] + v[1295] * v[145] + v[1293] * v[149] + v[1636] * v[181]
//			+ v[1642] * v[182] + v[1641] * v[183] + (*a4)*(v[3932 + i717] + v[216] * v[3950 + i717]) + v[1584] * v[787]
//			+ i1580 * v[867] + i1578 * v[868] + i1576 * v[869];
//		v[1308] = (-v[1299] - v[1297] * v[171] - v[1295] * v[175] - v[1293] * v[179] - v[1636] * v[196] - v[1640] * v[197]
//			- v[1635] * v[198] - v[1639] * v[198] + (*a4)*(v[3752 + i717] + v[216] * v[3770 + i717]) - v[1631] * v[792]
//			- i1594 * v[867] - i1598 * v[868] - i1602 * v[869] - v[215] * (v[1230] * v[197] + v[1237] * v[214] + v[792] * v[923]))
//			/ 2e0;
//		v[1309] = (-v[1299] - v[1297] * v[170] - v[1295] * v[174] - v[1293] * v[178] - v[1636] * v[193] - v[1640] * v[194]
//			- v[1635] * v[195] - v[1639] * v[195] + (*a4)*(v[3788 + i717] + v[216] * v[3806 + i717]) - v[1631] * v[791]
//			- i1593 * v[867] - i1597 * v[868] - i1601 * v[869] - v[215] * (v[1230] * v[194] + v[1236] * v[214] + v[791] * v[923]))
//			/ 2e0;
//		v[1261] = v[1261] - v[1293] * v[157] + v[1643];
//		v[1310] = (-v[1307] - 2e0*v[1233] * v[1380] + v[1297] * v[142] + v[1295] * v[146] + v[1293] * v[150] + v[1636] * v[184]
//			+ v[1642] * v[185] + v[1641] * v[186] + (*a4)*(v[3896 + i717] + v[216] * v[3914 + i717]) + v[1584] * v[788]
//			+ i1591 * v[867] + i1595 * v[868] + i1599 * v[869]) / 2e0;
//		v[1255] = v[1255] - v[1295] * v[157] + v[1645];
//		v[1249] = v[1249] - v[1297] * v[157] + v[1647];
//		v[1311] = (-v[1307] - 2e0*v[1234] * v[1380] + v[1297] * v[143] + v[1295] * v[147] + v[1293] * v[151] + v[1636] * v[187]
//			+ v[1642] * v[188] + v[1641] * v[189] + (*a4)*(v[3860 + i717] + v[216] * v[3878 + i717]) + v[1584] * v[789]
//			+ i1592 * v[867] + i1596 * v[868] + i1600 * v[869]) / 2e0;
//		v[1262] = v[1262] - v[1293] * v[155] + v[1644];
//		v[1256] = v[1256] - v[1295] * v[155] + v[1646];
//		v[1250] = v[1250] - v[1297] * v[155] + v[1648];
//		v[1263] = v[1263] - v[1293] * v[153] - v[1643] - v[1644];
//		v[1257] = v[1257] - v[1295] * v[153] - v[1645] - v[1646];
//		v[1251] = v[1251] - v[1297] * v[153] - v[1647] - v[1648];
//		v[1264] = v[1264] + v[129] * v[1293] + v[1649];
//		v[1258] = v[1258] + v[129] * v[1295] + v[1651];
//		v[1252] = v[1252] + v[129] * v[1297] + v[1653];
//		v[1265] = v[1265] + v[127] * v[1293] + v[1650];
//		v[1259] = v[1259] + v[127] * v[1295] + v[1652];
//		v[1253] = v[1253] + v[127] * v[1297] + v[1654];
//		v[1266] = v[1266] + v[125] * v[1293] - v[1649] - v[1650];
//		v[1260] = v[1260] + v[125] * v[1295] - v[1651] - v[1652];
//		v[1254] = v[1254] + v[125] * v[1297] - v[1653] - v[1654];
//		v[3969] = v[1254];
//		v[3970] = v[1260];
//		v[3971] = v[1266];
//		v[3972] = v[1253];
//		v[3973] = v[1259];
//		v[3974] = v[1265];
//		v[3975] = v[1252];
//		v[3976] = v[1258];
//		v[3977] = v[1264];
//		v[3978] = v[1251];
//		v[3979] = v[1257];
//		v[3980] = v[1263];
//		v[3981] = v[1250];
//		v[3982] = v[1256];
//		v[3983] = v[1262];
//		v[3984] = v[1249];
//		v[3985] = v[1255];
//		v[3986] = v[1261];
//		Rc[i717 - 1] += v[3100 + i717] + (*a4)*(v[3118 + i717] + (*ct)*(v[719] * (-(v[226] * v[872]) - v[231] * v[873]
//			+ v[237] * v[874] + v[242] * v[875]) + v[720] * (-(v[227] * v[872]) - v[233] * v[873] + v[238] * v[874] + v[244] * v[875])
//			+ v[721] * (-(v[228] * v[872]) - v[235] * v[873] + v[239] * v[874] + v[246] * v[875])) + v[216] * v[948]);
//		for (i878 = 1; i878 <= 18; i878++) {
//			v[1312] = v[2840 + i878];
//			v[1313] = v[2822 + i878];
//			v[1314] = v[2858 + i878];
//			v[1315] = v[2876 + i878];
//			Kc[i717 - 1][i878 - 1] += v[1310] * v[1312] + v[1311] * v[1313] + v[1309] * v[1314] + v[1308] * v[1315] + v[3968 + i878] +
//				(*a4)*(-(v[1091] * v[1312]) - v[1095] * v[1313] + v[1096] * v[1314] + v[1100] * v[1315] - v[3986 + i878] + (*ct)*(-
//				(v[1109] * v[4004 + i878]) - v[1108] * v[4022 + i878] - v[1107] * v[4040 + i878]) - v[216] * v[4058 + i878]);
//		};/* end for */
//	};/* end for */
//	v[1324] = 0e0;
//	v[1325] = 0e0;
//	v[1326] = 0e0;
//	v[1327] = 0e0;
//	b1328 = b478;
//	if (b1328) {
//		b1329 = b497;
//		b1333 = b480;
//		if (b1333) {
//			v[1326] = 0e0;
//			v[1325] = 0e0;
//			v[1324] = 0e0;
//			v[1327] = 0e0;
//		}
//		else {
//		};
//	}
//	else {
//	};
//	v[4301] = v[125] * v[1324];
//	v[4302] = 0e0;
//	v[4303] = 0e0;
//	v[4304] = v[127] * v[1324];
//	v[4305] = 0e0;
//	v[4306] = 0e0;
//	v[4307] = v[129] * v[1324];
//	v[4308] = 0e0;
//	v[4309] = 0e0;
//	v[4310] = -(v[1324] * v[153]);
//	v[4311] = 0e0;
//	v[4312] = 0e0;
//	v[4313] = -(v[1324] * v[155]);
//	v[4314] = 0e0;
//	v[4315] = 0e0;
//	v[4316] = -(v[1324] * v[157]);
//	v[4317] = 0e0;
//	v[4318] = 0e0;
//	v[1472] = 2e0*v[1324] * v[214];
//	v[1467] = -(v[1324] * v[451]);
//	v[1404] = 2e0*v[1057] * v[1324];
//	v[4265] = 0e0;
//	v[4266] = v[125] * v[1325];
//	v[4267] = 0e0;
//	v[4268] = 0e0;
//	v[4269] = v[127] * v[1325];
//	v[4270] = 0e0;
//	v[4271] = 0e0;
//	v[4272] = v[129] * v[1325];
//	v[4273] = 0e0;
//	v[4274] = 0e0;
//	v[4275] = -(v[1325] * v[153]);
//	v[4276] = 0e0;
//	v[4277] = 0e0;
//	v[4278] = -(v[1325] * v[155]);
//	v[4279] = 0e0;
//	v[4280] = 0e0;
//	v[4281] = -(v[1325] * v[157]);
//	v[4282] = 0e0;
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
//	v[4258] = 0e0;
//	v[4259] = 0e0;
//	v[4260] = 0e0;
//	v[4261] = 0e0;
//	v[4262] = v[1325];
//	v[4263] = v[1324];
//	v[4264] = 0e0;
//	v[4229] = 0e0;
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
//	v[4241] = v[1325];
//	v[4242] = v[1324];
//	v[4243] = 0e0;
//	v[4244] = 0e0;
//	v[4245] = 0e0;
//	v[4246] = 0e0;
//	v[4211] = 0e0;
//	v[4212] = 0e0;
//	v[4213] = 0e0;
//	v[4214] = 0e0;
//	v[4215] = 0e0;
//	v[4216] = 0e0;
//	v[4217] = 0e0;
//	v[4218] = 0e0;
//	v[4219] = 0e0;
//	v[4220] = v[1325];
//	v[4221] = v[1324];
//	v[4222] = 0e0;
//	v[4223] = 0e0;
//	v[4224] = 0e0;
//	v[4225] = 0e0;
//	v[4226] = 0e0;
//	v[4227] = 0e0;
//	v[4228] = 0e0;
//	v[4193] = 0e0;
//	v[4194] = 0e0;
//	v[4195] = 0e0;
//	v[4196] = 0e0;
//	v[4197] = 0e0;
//	v[4198] = 0e0;
//	v[4199] = v[1325];
//	v[4200] = v[1324];
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
//	v[4175] = 0e0;
//	v[4176] = 0e0;
//	v[4177] = 0e0;
//	v[4178] = v[1325];
//	v[4179] = v[1324];
//	v[4180] = 0e0;
//	v[4181] = 0e0;
//	v[4182] = 0e0;
//	v[4183] = 0e0;
//	v[4184] = 0e0;
//	v[4185] = 0e0;
//	v[4186] = 0e0;
//	v[4187] = 0e0;
//	v[4188] = 0e0;
//	v[4189] = 0e0;
//	v[4190] = 0e0;
//	v[4191] = 0e0;
//	v[4192] = 0e0;
//	v[4157] = v[1325];
//	v[4158] = v[1324];
//	v[4159] = 0e0;
//	v[4160] = 0e0;
//	v[4161] = 0e0;
//	v[4162] = 0e0;
//	v[4163] = 0e0;
//	v[4164] = 0e0;
//	v[4165] = 0e0;
//	v[4166] = 0e0;
//	v[4167] = 0e0;
//	v[4168] = 0e0;
//	v[4169] = 0e0;
//	v[4170] = 0e0;
//	v[4171] = 0e0;
//	v[4172] = 0e0;
//	v[4173] = 0e0;
//	v[4174] = 0e0;
//	v[1475] = 2e0*v[1325] * v[215];
//	v[1468] = -(v[1325] * v[459]);
//	v[1396] = 2e0*v[1060] * v[1325];
//	v[1666] = 2e0*v[1326];
//	v[1658] = v[1326] * v[216];
//	v[1478] = 2e0*v[1658];
//	v[1518] = -(v[1478] * v[157]) / 2e0;
//	v[1516] = -(v[1478] * v[155]) / 2e0;
//	v[1514] = -(v[1478] * v[153]) / 2e0;
//	v[1512] = (v[129] * v[1478]) / 2e0;
//	v[1510] = (v[127] * v[1478]) / 2e0;
//	v[1507] = -v[1468] + (v[1478] * v[215]) / 2e0;
//	v[1506] = (v[125] * v[1478]) / 2e0;
//	v[1505] = -v[1467] + (v[1478] * v[214]) / 2e0;
//	v[1471] = -(v[1326] * v[215]);
//	v[1470] = -(v[1326] * v[214]);
//	v[1469] = -(v[1326] * v[467]);
//	v[4099] = -(v[125] * v[1467]) + v[1325] * v[453];
//	v[4100] = -(v[125] * v[1468]) + v[1324] * v[453];
//	v[4101] = -(v[125] * v[1469]);
//	v[4102] = -(v[127] * v[1467]) + v[1325] * v[454];
//	v[4103] = -(v[127] * v[1468]) + v[1324] * v[454];
//	v[4104] = -(v[127] * v[1469]);
//	v[4105] = -(v[129] * v[1467]) + v[1325] * v[455];
//	v[4106] = -(v[129] * v[1468]) + v[1324] * v[455];
//	v[4107] = -(v[129] * v[1469]);
//	v[4108] = v[1467] * v[153] + v[1325] * v[456];
//	v[4109] = v[1468] * v[153] + v[1324] * v[456];
//	v[4110] = v[1469] * v[153];
//	v[4111] = v[1467] * v[155] + v[1325] * v[457];
//	v[4112] = v[1468] * v[155] + v[1324] * v[457];
//	v[4113] = v[1469] * v[155];
//	v[4114] = v[1467] * v[157] + v[1325] * v[458];
//	v[4115] = v[1468] * v[157] + v[1324] * v[458];
//	v[4116] = v[1469] * v[157];
//	v[4499] = 0e0;
//	v[4500] = 0e0;
//	v[4501] = 0e0;
//	v[4502] = 0e0;
//	v[4503] = 0e0;
//	v[4504] = 0e0;
//	v[4505] = -v[1467];
//	v[4506] = -v[1468];
//	v[4507] = -v[1469];
//	v[4508] = 0e0;
//	v[4509] = 0e0;
//	v[4510] = 0e0;
//	v[4511] = 0e0;
//	v[4512] = 0e0;
//	v[4513] = 0e0;
//	v[4514] = 0e0;
//	v[4515] = 0e0;
//	v[4516] = 0e0;
//	v[4553] = 0e0;
//	v[4554] = 0e0;
//	v[4555] = 0e0;
//	v[4556] = -v[1467];
//	v[4557] = -v[1468];
//	v[4558] = -v[1469];
//	v[4559] = 0e0;
//	v[4560] = 0e0;
//	v[4561] = 0e0;
//	v[4562] = 0e0;
//	v[4563] = 0e0;
//	v[4564] = 0e0;
//	v[4565] = 0e0;
//	v[4566] = 0e0;
//	v[4567] = 0e0;
//	v[4568] = 0e0;
//	v[4569] = 0e0;
//	v[4570] = 0e0;
//	v[4391] = 0e0;
//	v[4392] = 0e0;
//	v[4393] = 0e0;
//	v[4394] = 0e0;
//	v[4395] = 0e0;
//	v[4396] = 0e0;
//	v[4397] = 0e0;
//	v[4398] = 0e0;
//	v[4399] = 0e0;
//	v[4400] = 0e0;
//	v[4401] = 0e0;
//	v[4402] = 0e0;
//	v[4403] = v[1467];
//	v[4404] = v[1468];
//	v[4405] = v[1469];
//	v[4406] = 0e0;
//	v[4407] = 0e0;
//	v[4408] = 0e0;
//	v[4337] = 0e0;
//	v[4338] = 0e0;
//	v[4339] = 0e0;
//	v[4340] = 0e0;
//	v[4341] = 0e0;
//	v[4342] = 0e0;
//	v[4343] = 0e0;
//	v[4344] = 0e0;
//	v[4345] = 0e0;
//	v[4346] = 0e0;
//	v[4347] = 0e0;
//	v[4348] = 0e0;
//	v[4349] = 0e0;
//	v[4350] = 0e0;
//	v[4351] = 0e0;
//	v[4352] = v[1467];
//	v[4353] = v[1468];
//	v[4354] = v[1469];
//	v[4607] = -v[1467];
//	v[4608] = -v[1468];
//	v[4609] = -v[1469];
//	v[4610] = 0e0;
//	v[4611] = 0e0;
//	v[4612] = 0e0;
//	v[4613] = 0e0;
//	v[4614] = 0e0;
//	v[4615] = 0e0;
//	v[4616] = 0e0;
//	v[4617] = 0e0;
//	v[4618] = 0e0;
//	v[4619] = 0e0;
//	v[4620] = 0e0;
//	v[4621] = 0e0;
//	v[4622] = 0e0;
//	v[4623] = 0e0;
//	v[4624] = 0e0;
//	v[4445] = 0e0;
//	v[4446] = 0e0;
//	v[4447] = 0e0;
//	v[4448] = 0e0;
//	v[4449] = 0e0;
//	v[4450] = 0e0;
//	v[4451] = 0e0;
//	v[4452] = 0e0;
//	v[4453] = 0e0;
//	v[4454] = v[1467];
//	v[4455] = v[1468];
//	v[4456] = v[1469];
//	v[4457] = 0e0;
//	v[4458] = 0e0;
//	v[4459] = 0e0;
//	v[4460] = 0e0;
//	v[4461] = 0e0;
//	v[4462] = 0e0;
//	v[1403] = -(v[1326] * v[157]);
//	v[1402] = -(v[1326] * v[155]);
//	v[1401] = -(v[1326] * v[153]);
//	v[1400] = v[129] * v[1326];
//	v[1399] = v[127] * v[1326];
//	v[1398] = v[125] * v[1326];
//	v[4139] = 0e0;
//	v[4140] = 0e0;
//	v[4141] = v[1398];
//	v[4142] = 0e0;
//	v[4143] = 0e0;
//	v[4144] = v[1399];
//	v[4145] = 0e0;
//	v[4146] = 0e0;
//	v[4147] = v[1400];
//	v[4148] = 0e0;
//	v[4149] = 0e0;
//	v[4150] = v[1401];
//	v[4151] = 0e0;
//	v[4152] = 0e0;
//	v[4153] = v[1402];
//	v[4154] = 0e0;
//	v[4155] = 0e0;
//	v[4156] = v[1403];
//	v[4319] = v[1398];
//	v[4320] = 0e0;
//	v[4321] = 0e0;
//	v[4322] = v[1399];
//	v[4323] = 0e0;
//	v[4324] = 0e0;
//	v[4325] = v[1400];
//	v[4326] = 0e0;
//	v[4327] = 0e0;
//	v[4328] = v[1401];
//	v[4329] = 0e0;
//	v[4330] = 0e0;
//	v[4331] = v[1402];
//	v[4332] = 0e0;
//	v[4333] = 0e0;
//	v[4334] = v[1403];
//	v[4335] = 0e0;
//	v[4336] = 0e0;
//	v[4283] = 0e0;
//	v[4284] = v[1398];
//	v[4285] = 0e0;
//	v[4286] = 0e0;
//	v[4287] = v[1399];
//	v[4288] = 0e0;
//	v[4289] = 0e0;
//	v[4290] = v[1400];
//	v[4291] = 0e0;
//	v[4292] = 0e0;
//	v[4293] = v[1401];
//	v[4294] = 0e0;
//	v[4295] = 0e0;
//	v[4296] = v[1402];
//	v[4297] = 0e0;
//	v[4298] = 0e0;
//	v[4299] = v[1403];
//	v[4300] = 0e0;
//	v[1392] = 2e0*v[1326] * v[945];
//	v[1334] = 0e0;
//	v[1335] = 0e0;
//	v[1336] = 0e0;
//	v[1337] = 0e0;
//	b1338 = b478;
//	if (b1338) {
//		b1339 = b480;
//		if (b1339) {
//			v[1336] = 0e0;
//			v[1335] = 0e0;
//			v[1334] = 0e0;
//			v[1337] = -v[1327];
//		}
//		else {
//		};
//	}
//	else {
//	};
//	v[1340] = (v[1472] + v[1475]) / 2e0;
//	v[4517] = 0e0;
//	v[4518] = 0e0;
//	v[4519] = 0e0;
//	v[4520] = 0e0;
//	v[4521] = 0e0;
//	v[4522] = 0e0;
//	v[4523] = -v[1470];
//	v[4524] = -v[1471];
//	v[4525] = v[1340];
//	v[4526] = 0e0;
//	v[4527] = 0e0;
//	v[4528] = 0e0;
//	v[4529] = 0e0;
//	v[4530] = 0e0;
//	v[4531] = 0e0;
//	v[4532] = 0e0;
//	v[4533] = 0e0;
//	v[4534] = 0e0;
//	v[4571] = 0e0;
//	v[4572] = 0e0;
//	v[4573] = 0e0;
//	v[4574] = -v[1470];
//	v[4575] = -v[1471];
//	v[4576] = v[1340];
//	v[4577] = 0e0;
//	v[4578] = 0e0;
//	v[4579] = 0e0;
//	v[4580] = 0e0;
//	v[4581] = 0e0;
//	v[4582] = 0e0;
//	v[4583] = 0e0;
//	v[4584] = 0e0;
//	v[4585] = 0e0;
//	v[4586] = 0e0;
//	v[4587] = 0e0;
//	v[4588] = 0e0;
//	v[4409] = 0e0;
//	v[4410] = 0e0;
//	v[4411] = 0e0;
//	v[4412] = 0e0;
//	v[4413] = 0e0;
//	v[4414] = 0e0;
//	v[4415] = 0e0;
//	v[4416] = 0e0;
//	v[4417] = 0e0;
//	v[4418] = 0e0;
//	v[4419] = 0e0;
//	v[4420] = 0e0;
//	v[4421] = v[1470];
//	v[4422] = v[1471];
//	v[4423] = -v[1340];
//	v[4424] = 0e0;
//	v[4425] = 0e0;
//	v[4426] = 0e0;
//	v[4355] = 0e0;
//	v[4356] = 0e0;
//	v[4357] = 0e0;
//	v[4358] = 0e0;
//	v[4359] = 0e0;
//	v[4360] = 0e0;
//	v[4361] = 0e0;
//	v[4362] = 0e0;
//	v[4363] = 0e0;
//	v[4364] = 0e0;
//	v[4365] = 0e0;
//	v[4366] = 0e0;
//	v[4367] = 0e0;
//	v[4368] = 0e0;
//	v[4369] = 0e0;
//	v[4370] = v[1470];
//	v[4371] = v[1471];
//	v[4372] = -v[1340];
//	v[4625] = -v[1470];
//	v[4626] = -v[1471];
//	v[4627] = v[1340];
//	v[4628] = 0e0;
//	v[4629] = 0e0;
//	v[4630] = 0e0;
//	v[4631] = 0e0;
//	v[4632] = 0e0;
//	v[4633] = 0e0;
//	v[4634] = 0e0;
//	v[4635] = 0e0;
//	v[4636] = 0e0;
//	v[4637] = 0e0;
//	v[4638] = 0e0;
//	v[4639] = 0e0;
//	v[4640] = 0e0;
//	v[4641] = 0e0;
//	v[4642] = 0e0;
//	v[4463] = 0e0;
//	v[4464] = 0e0;
//	v[4465] = 0e0;
//	v[4466] = 0e0;
//	v[4467] = 0e0;
//	v[4468] = 0e0;
//	v[4469] = 0e0;
//	v[4470] = 0e0;
//	v[4471] = 0e0;
//	v[4472] = v[1470];
//	v[4473] = v[1471];
//	v[4474] = -v[1340];
//	v[4475] = 0e0;
//	v[4476] = 0e0;
//	v[4477] = 0e0;
//	v[4478] = 0e0;
//	v[4479] = 0e0;
//	v[4480] = 0e0;
//	v[1519] = -(v[1340] * v[157]);
//	v[1517] = -(v[1340] * v[155]);
//	v[1515] = -(v[1340] * v[153]);
//	v[1513] = v[129] * v[1340];
//	v[1511] = v[127] * v[1340];
//	v[1509] = v[125] * v[1340];
//	v[4117] = -(v[125] * v[1470]);
//	v[4118] = -(v[125] * v[1471]);
//	v[4119] = v[1509];
//	v[4120] = -(v[127] * v[1470]);
//	v[4121] = -(v[127] * v[1471]);
//	v[4122] = v[1511];
//	v[4123] = -(v[129] * v[1470]);
//	v[4124] = -(v[129] * v[1471]);
//	v[4125] = v[1513];
//	v[4126] = v[1470] * v[153];
//	v[4127] = v[1471] * v[153];
//	v[4128] = v[1515];
//	v[4129] = v[1470] * v[155];
//	v[4130] = v[1471] * v[155];
//	v[4131] = v[1517];
//	v[4132] = v[1470] * v[157];
//	v[4133] = v[1471] * v[157];
//	v[4134] = v[1519];
//	v[1508] = -v[1469] + v[1340] * v[216];
//	v[1503] = v[1340] * v[189];
//	v[1973] = v[1503] + v[1478] * v[189];
//	v[1896] = v[1503] + v[1326] * v[1669];
//	v[1497] = v[1340] * v[186];
//	v[1967] = v[1497] + v[1478] * v[186];
//	v[1895] = v[1497] + v[1326] * v[1668];
//	v[1491] = -(v[1340] * v[195]);
//	v[1961] = v[1491] - v[1478] * v[195];
//	v[1899] = v[1491] - v[1326] * v[1686];
//	v[1485] = -(v[1340] * v[198]);
//	v[1955] = v[1485] - v[1478] * v[198];
//	v[1898] = v[1485] - v[1326] * v[1683];
//	v[1479] = v[1340] * v[183];
//	v[1949] = v[1479] + v[1478] * v[183];
//	v[1894] = v[1479] + v[1326] * v[1667];
//	v[1465] = -(v[1340] * v[192]);
//	v[1943] = v[1465] - v[1478] * v[192];
//	v[1897] = v[1465] - v[1326] * v[1675];
//	v[1341] = v[1325] * v[181] + v[1324] * v[182];
//	v[1682] = v[1658] * v[182] + v[1341] * v[214];
//	v[1951] = v[1682] + v[1475] * v[182];
//	v[1474] = v[1341] * v[215];
//	v[1681] = v[1474] + v[1658] * v[181];
//	v[1950] = v[1681] + v[1472] * v[181];
//	v[1342] = v[1325] * v[184] + v[1324] * v[185];
//	v[1690] = v[1658] * v[185] + v[1342] * v[214];
//	v[1969] = v[1690] + v[1475] * v[185];
//	v[1494] = v[1342] * v[215];
//	v[1689] = v[1494] + v[1658] * v[184];
//	v[1968] = v[1689] + v[1472] * v[184];
//	v[1343] = v[1325] * v[187] + v[1324] * v[188];
//	v[1692] = v[1658] * v[188] + v[1343] * v[214];
//	v[1975] = v[1692] + v[1475] * v[188];
//	v[1500] = v[1343] * v[215];
//	v[1691] = v[1500] + v[1658] * v[187];
//	v[1974] = v[1691] + v[1472] * v[187];
//	v[1344] = v[1325] * v[190] + v[1324] * v[191];
//	v[1677] = -(v[1658] * v[191]) - v[1344] * v[214];
//	v[1945] = v[1677] - v[1475] * v[191];
//	v[1460] = -(v[1344] * v[215]);
//	v[1676] = v[1460] - v[1658] * v[190];
//	v[1944] = v[1676] - v[1472] * v[190];
//	v[1345] = v[1325] * v[193] + v[1324] * v[194];
//	v[1688] = -(v[1658] * v[194]) - v[1345] * v[214];
//	v[1963] = v[1688] - v[1475] * v[194];
//	v[1488] = -(v[1345] * v[215]);
//	v[1687] = v[1488] - v[1658] * v[193];
//	v[1962] = v[1687] - v[1472] * v[193];
//	v[1346] = v[1325] * v[196] + v[1324] * v[197];
//	v[1685] = -(v[1658] * v[197]) - v[1346] * v[214];
//	v[1957] = v[1685] - v[1475] * v[197];
//	v[1482] = -(v[1346] * v[215]);
//	v[1684] = v[1482] - v[1658] * v[196];
//	v[1956] = v[1684] - v[1472] * v[196];
//	v[1405] = v[125] * v[1341] + v[127] * v[1342] + v[129] * v[1343] - v[1344] * v[153] - v[1345] * v[155] - v[1346] * v[157];
//	v[1336] = v[1066] * v[1326] + v[1336] + v[1392] * v[216] + v[1340] * v[945];
//	v[1335] = v[1065] * v[1326] + v[1335] + v[1405] * v[214] + v[1396] * v[215] + v[1325] * v[794];
//	v[1334] = v[1064] * v[1326] + v[1334] + v[1404] * v[214] + v[1405] * v[215] + v[1324] * v[794];
//	v[1451] = v[1334] * v[199] + v[1335] * v[200] + v[1336] * v[201];
//	v[1337] = v[1337] + v[1451] * v[212];
//	v[1659] = v[1337] / v[477];
//	v[1347] = v[1659] * v[201] + v[1336] * v[213] + v[501];
//	v[1348] = v[1659] * v[200] + v[1335] * v[213] + v[500];
//	v[1349] = v[1659] * v[199] + v[1334] * v[213] + v[499];
//	v[4081] = v[125] * v[1349];
//	v[4082] = v[125] * v[1348];
//	v[4083] = v[125] * v[1347];
//	v[4084] = v[127] * v[1349];
//	v[4085] = v[127] * v[1348];
//	v[4086] = v[127] * v[1347];
//	v[4087] = v[129] * v[1349];
//	v[4088] = v[129] * v[1348];
//	v[4089] = v[129] * v[1347];
//	v[4090] = -(v[1349] * v[153]);
//	v[4091] = -(v[1348] * v[153]);
//	v[4092] = -(v[1347] * v[153]);
//	v[4093] = -(v[1349] * v[155]);
//	v[4094] = -(v[1348] * v[155]);
//	v[4095] = -(v[1347] * v[155]);
//	v[4096] = -(v[1349] * v[157]);
//	v[4097] = -(v[1348] * v[157]);
//	v[4098] = -(v[1347] * v[157]);
//	v[4535] = 0e0;
//	v[4536] = 0e0;
//	v[4537] = 0e0;
//	v[4538] = 0e0;
//	v[4539] = 0e0;
//	v[4540] = 0e0;
//	v[4541] = v[1349];
//	v[4542] = v[1348];
//	v[4543] = v[1347];
//	v[4544] = 0e0;
//	v[4545] = 0e0;
//	v[4546] = 0e0;
//	v[4547] = 0e0;
//	v[4548] = 0e0;
//	v[4549] = 0e0;
//	v[4550] = 0e0;
//	v[4551] = 0e0;
//	v[4552] = 0e0;
//	v[4589] = 0e0;
//	v[4590] = 0e0;
//	v[4591] = 0e0;
//	v[4592] = v[1349];
//	v[4593] = v[1348];
//	v[4594] = v[1347];
//	v[4595] = 0e0;
//	v[4596] = 0e0;
//	v[4597] = 0e0;
//	v[4598] = 0e0;
//	v[4599] = 0e0;
//	v[4600] = 0e0;
//	v[4601] = 0e0;
//	v[4602] = 0e0;
//	v[4603] = 0e0;
//	v[4604] = 0e0;
//	v[4605] = 0e0;
//	v[4606] = 0e0;
//	v[4427] = 0e0;
//	v[4428] = 0e0;
//	v[4429] = 0e0;
//	v[4430] = 0e0;
//	v[4431] = 0e0;
//	v[4432] = 0e0;
//	v[4433] = 0e0;
//	v[4434] = 0e0;
//	v[4435] = 0e0;
//	v[4436] = 0e0;
//	v[4437] = 0e0;
//	v[4438] = 0e0;
//	v[4439] = -v[1349];
//	v[4440] = -v[1348];
//	v[4441] = -v[1347];
//	v[4442] = 0e0;
//	v[4443] = 0e0;
//	v[4444] = 0e0;
//	v[4373] = 0e0;
//	v[4374] = 0e0;
//	v[4375] = 0e0;
//	v[4376] = 0e0;
//	v[4377] = 0e0;
//	v[4378] = 0e0;
//	v[4379] = 0e0;
//	v[4380] = 0e0;
//	v[4381] = 0e0;
//	v[4382] = 0e0;
//	v[4383] = 0e0;
//	v[4384] = 0e0;
//	v[4385] = 0e0;
//	v[4386] = 0e0;
//	v[4387] = 0e0;
//	v[4388] = -v[1349];
//	v[4389] = -v[1348];
//	v[4390] = -v[1347];
//	v[4643] = v[1349];
//	v[4644] = v[1348];
//	v[4645] = v[1347];
//	v[4646] = 0e0;
//	v[4647] = 0e0;
//	v[4648] = 0e0;
//	v[4649] = 0e0;
//	v[4650] = 0e0;
//	v[4651] = 0e0;
//	v[4652] = 0e0;
//	v[4653] = 0e0;
//	v[4654] = 0e0;
//	v[4655] = 0e0;
//	v[4656] = 0e0;
//	v[4657] = 0e0;
//	v[4658] = 0e0;
//	v[4659] = 0e0;
//	v[4660] = 0e0;
//	v[4481] = 0e0;
//	v[4482] = 0e0;
//	v[4483] = 0e0;
//	v[4484] = 0e0;
//	v[4485] = 0e0;
//	v[4486] = 0e0;
//	v[4487] = 0e0;
//	v[4488] = 0e0;
//	v[4489] = 0e0;
//	v[4490] = -v[1349];
//	v[4491] = -v[1348];
//	v[4492] = -v[1347];
//	v[4493] = 0e0;
//	v[4494] = 0e0;
//	v[4495] = 0e0;
//	v[4496] = 0e0;
//	v[4497] = 0e0;
//	v[4498] = 0e0;
//	v[1351] = -(v[1349] * v[169]) - v[1348] * v[173] - v[1347] * v[177] + v[1467] * v[190] + v[1468] * v[191]
//		+ v[1460] * v[214] + v[1465] * v[216] - v[1326] * (v[1660] * v[190] + v[1661] * v[191] + v[192] * v[467]);
//	v[1353] = v[1349] * v[141] + v[1348] * v[145] + v[1347] * v[149] - v[1467] * v[181] - v[1468] * v[182] + v[1474] * v[214]
//		+ v[1479] * v[216] + v[1326] * (v[1667] * v[216] + v[183] * v[467]);
//	v[1354] = (-v[1351] - v[1349] * v[171] - v[1348] * v[175] - v[1347] * v[179] + v[1467] * v[196] + v[1468] * v[197]
//		+ v[1482] * v[214] + v[1485] * v[216] - v[1326] * (v[1660] * v[196] + v[1661] * v[197] + v[198] * v[467])) / 2e0;
//	v[1355] = (-v[1351] - v[1349] * v[170] - v[1348] * v[174] - v[1347] * v[178] + v[1467] * v[193] + v[1468] * v[194]
//		+ v[1488] * v[214] + v[1491] * v[216] - v[1326] * (v[1660] * v[193] + v[1661] * v[194] + v[195] * v[467])) / 2e0;
//	v[1356] = (-v[1353] + v[1349] * v[142] + v[1348] * v[146] + v[1347] * v[150] - v[1467] * v[184] - v[1468] * v[185]
//		+ v[1494] * v[214] + v[1497] * v[216] + v[1326] * (v[1668] * v[216] + v[186] * v[467])) / 2e0;
//	v[1357] = (-v[1353] + v[1349] * v[143] + v[1348] * v[147] + v[1347] * v[151] - v[1467] * v[187] - v[1468] * v[188]
//		+ v[1500] * v[214] + v[1503] * v[216] + v[1326] * (v[1669] * v[216] + v[189] * v[467])) / 2e0;
//	for (i1319 = 1; i1319 <= 18; i1319++) {
//		v[1391] = v[4116 + i1319];
//		v[1389] = v[2952 + i1319];
//		v[1365] = v[2876 + i1319];
//		v[1696] = -v[1365] / 2e0;
//		v[1364] = v[2858 + i1319];
//		v[1695] = -v[1364] / 2e0;
//		v[1363] = v[2840 + i1319];
//		v[1693] = v[1363] / 2e0;
//		v[1362] = v[2822 + i1319];
//		v[1694] = v[1362] / 2e0;
//		v[1366] = (-v[1362] - v[1363]) / 2e0;
//		v[1367] = (-v[1364] - v[1365]) / 2e0;
//		v[1673] = v[1366] * v[183] + (v[1363] * v[186]) / 2e0 + (v[1362] * v[189]) / 2e0 - v[1367] * v[192] - (v[1364] * v[195])
//			/ 2e0 - (v[1365] * v[198]) / 2e0;
//		v[1369] = (2e0*v[1366] * v[141] + v[1363] * v[142] + v[1362] * v[143] - 2e0*v[1367] * v[169] - v[1364] * v[170]
//			- v[1365] * v[171] + 2e0*v[2916 + i1319]) / 2e0;
//		v[1371] = (2e0*v[1366] * v[145] + v[1363] * v[146] + v[1362] * v[147] - 2e0*v[1367] * v[173] - v[1364] * v[174]
//			- v[1365] * v[175] + 2e0*v[2934 + i1319]) / 2e0;
//		v[1372] = (2e0*v[1389] + 2e0*v[1366] * v[149] + v[1363] * v[150] + v[1362] * v[151] - 2e0*v[1367] * v[177]
//			- v[1364] * v[178] - v[1365] * v[179]) / 2e0;
//		v[1662] = v[1369] * v[199] + v[1371] * v[200] + v[1372] * v[201];
//		v[1452] = v[1662] / v[477];
//		v[1373] = v[1369];
//		v[1423] = v[1373];
//		v[1374] = v[1371];
//		v[1422] = v[1374];
//		v[1375] = v[1372];
//		v[1421] = v[1375];
//		v[1376] = -(v[1337] * v[1662] * v[919]);
//		v[1377] = v[1452];
//		v[1663] = v[1377] * v[212];
//		v[1461] = v[1663] * v[199] + v[1369] * v[213];
//		v[1393] = v[1663] * v[201] + v[1372] * v[213];
//		v[1381] = v[1663] * v[200] + v[1371] * v[213];
//		v[1378] = v[1461];
//		v[1665] = v[1378] * v[215];
//		v[1664] = v[1665] + v[1381] * v[214];
//		v[1379] = v[1365] * v[1380] - v[157] * v[1664];
//		v[1382] = v[1364] * v[1380] - v[155] * v[1664];
//		v[1383] = -(v[153] * v[1665]) - v[214] * (v[1381] * v[153] + v[1367] * v[215]);
//		v[1384] = -(v[1362] * v[1380]) + v[129] * v[1664];
//		v[1385] = -(v[1363] * v[1380]) + v[127] * v[1664];
//		v[1386] = v[125] * v[1665] + v[214] * (v[125] * v[1381] + v[1366] * v[215]);
//		v[1387] = v[1324] * v[1378] + v[1325] * v[1381];
//		v[1388] = v[1381];
//		v[1390] = ((*a4)*v[1389] + v[1673])*v[216] + v[1393] * v[945];
//		v[1394] = (2e0*(*a4)*v[1391] + 2e0*v[1392] * v[1393] + v[1057] * v[1378] * v[1666] + v[1060] * v[1388] * v[1666]
//			+ 2e0*v[1366] * v[1894] + v[1363] * v[1895] + v[1362] * v[1896] + 2e0*v[1367] * v[1897] + v[1365] * v[1898]
//			+ v[1364] * v[1899]) / 2e0;
//		v[1434] = v[1394];
//		v[1397] = v[1325] * v[1390] + v[1388] * v[1396] + v[1378] * v[1405] + v[1367] * v[1677] + v[1366] * v[1682] +
//			(v[1365] * v[1685]) / 2e0 + (v[1364] * v[1688]) / 2e0 + (v[1363] * v[1690]) / 2e0 + (v[1362] * v[1692]) / 2e0
//			+ v[1585] * v[4282 + i1319];
//		v[1435] = v[1397];
//		v[1406] = v[1324] * v[1390] + v[1378] * v[1404] + v[1388] * v[1405] + v[1367] * v[1676] + v[1366] * v[1681] +
//			(v[1365] * v[1684]) / 2e0 + (v[1364] * v[1687]) / 2e0 + (v[1363] * v[1689]) / 2e0 + (v[1362] * v[1691]) / 2e0
//			+ v[1585] * v[4318 + i1319];
//		v[1436] = v[1406];
//		b1407 = b478;
//		if (b1407) {
//			b1408 = b480;
//			if (b1408) {
//				v[1377] = 0e0;
//				v[1378] = 0e0;
//				v[1388] = 0e0;
//			}
//			else {
//			};
//		}
//		else {
//		};
//		v[1409] = 0e0;
//		v[1410] = 0e0;
//		v[1411] = 0e0;
//		v[1412] = 0e0;
//		v[1413] = 0e0;
//		v[1414] = 0e0;
//		v[1415] = 0e0;
//		b1416 = b478;
//		if (b1416) {
//			v[1417] = 0e0;
//			v[1418] = 0e0;
//			v[1419] = 0e0;
//			b1420 = b497;
//			if (b1420) {
//				v[1419] = v[1375];
//				v[1375] = 0e0;
//				v[1418] = v[1374];
//				v[1374] = 0e0;
//				v[1417] = v[1373];
//				v[1373] = 0e0;
//			}
//			else {
//				v[1414] = -v[1421];
//				v[1375] = 0e0;
//				v[1413] = -v[1422];
//				v[1374] = 0e0;
//				v[1412] = -v[1423];
//				v[1373] = 0e0;
//			};
//			v[1670] = (v[1417] * v[452] + v[1418] * v[466] + v[1419] * v[468])*(*zetan);
//			b1424 = b480;
//			if (b1424) {
//				v[1411] = v[1419] * v[492];
//				v[1410] = v[1418] * v[492];
//				v[1409] = v[1417] * v[492];
//				v[1415] = (v[1670] * v[1921] * v[490] * v[759]) / v[760];
//			}
//			else {
//				v[1411] = v[1419] * v[496];
//				v[1410] = v[1418] * v[496];
//				v[1409] = v[1417] * v[496];
//				v[1376] = v[1376] + (v[1670] * v[1922] * v[495] * v[764]) / v[765];
//			};
//		}
//		else {
//		};
//		v[1678] = v[1409] * v[451];
//		v[1679] = v[1410] * v[459];
//		v[1437] = v[1376];
//		v[1433] = v[1412];
//		v[1432] = v[1413];
//		v[1431] = v[1414];
//		b1429 = b478;
//		if (b1429) {
//			v[1672] = -(v[1433] * v[214]) - v[1432] * v[215] - v[1431] * v[216];
//			b1430 = b480;
//			if (b1430) {
//				v[1394] = v[1394] + v[1414] * v[1671];
//				v[1414] = 0e0;
//				v[1397] = v[1397] + v[1413] * v[1671];
//				v[1413] = 0e0;
//				v[1406] = v[1406] + v[1412] * v[1671];
//				v[1412] = 0e0;
//				v[1415] = v[1415] + v[1545] * v[1672] * v[1926];
//				v[1376] = v[1376] - v[1415];
//			}
//			else {
//				v[1394] = v[1434] + v[1431] * v[486];
//				v[1414] = 0e0;
//				v[1397] = v[1435] + v[1432] * v[486];
//				v[1413] = 0e0;
//				v[1406] = v[1436] + v[1433] * v[486];
//				v[1412] = 0e0;
//				v[1376] = v[1437] + (*n2)*v[1672] * v[1930] * v[474];
//			};
//		}
//		else {
//		};
//		v[1439] = v[1326] * v[1393] + v[1411] * v[216];
//		v[1698] = v[1679] + v[1439] * v[215];
//		v[1697] = v[1678] + v[1439] * v[214];
//		v[1440] = v[1410] * v[181] + v[1409] * v[182] + (*a4)*v[4156 + i1319];
//		v[1441] = v[1410] * v[184] + v[1409] * v[185] + (*a4)*v[4174 + i1319];
//		v[1442] = v[1410] * v[187] + v[1409] * v[188] + (*a4)*v[4192 + i1319];
//		v[1443] = v[1410] * v[190] + v[1409] * v[191] + (*a4)*v[4210 + i1319];
//		v[1444] = v[1410] * v[193] + v[1409] * v[194] + (*a4)*v[4228 + i1319];
//		v[1445] = v[1410] * v[196] + v[1409] * v[197] + (*a4)*v[4246 + i1319];
//		v[1674] = v[125] * v[1440] + v[127] * v[1441] + v[129] * v[1442] - v[1443] * v[153] - v[1444] * v[155] - v[1445] * v[157];
//		v[1446] = v[1387] + v[1409] * v[214] + v[1410] * v[215];
//		v[1680] = -(v[1446] * v[216]) - v[1411] * v[467];
//		v[4679] = v[1325] * v[1386] + v[1366] * v[1505] + v[1461] * (v[125] * v[1472] + v[1506]) + v[125] * v[1697]
//			+ v[1410] * v[453];
//		v[4680] = v[1324] * v[1386] + v[1381] * (v[125] * v[1475] + v[1506]) + v[1366] * v[1507] + v[125] * v[1698]
//			+ v[1409] * v[453];
//		v[4681] = v[1366] * v[1508] + v[1393] * (2e0*v[1506] + v[1509]) - v[125] * v[1680];
//		v[4682] = v[1325] * v[1385] + v[1461] * (v[127] * v[1472] + v[1510]) + v[1505] * v[1693] + v[127] * v[1697]
//			+ v[1410] * v[454];
//		v[4683] = v[1324] * v[1385] + v[1381] * (v[127] * v[1475] + v[1510]) + v[1507] * v[1693] + v[127] * v[1698]
//			+ v[1409] * v[454];
//		v[4684] = v[1393] * (2e0*v[1510] + v[1511]) - v[127] * v[1680] + v[1508] * v[1693];
//		v[4685] = v[1325] * v[1384] + v[1461] * (v[129] * v[1472] + v[1512]) + v[1505] * v[1694] + v[129] * v[1697]
//			+ v[1410] * v[455];
//		v[4686] = v[1324] * v[1384] + v[1381] * (v[129] * v[1475] + v[1512]) + v[1507] * v[1694] + v[129] * v[1698]
//			+ v[1409] * v[455];
//		v[4687] = v[1393] * (2e0*v[1512] + v[1513]) - v[129] * v[1680] + v[1508] * v[1694];
//		v[4688] = v[1325] * v[1383] - v[1367] * v[1505] + v[1461] * (v[1514] - v[1472] * v[153]) - v[153] * v[1697]
//			+ v[1410] * v[456];
//		v[4689] = v[1324] * v[1383] - v[1367] * v[1507] + v[1381] * (v[1514] - v[1475] * v[153]) - v[153] * v[1698]
//			+ v[1409] * v[456];
//		v[4690] = -(v[1367] * v[1508]) + v[1393] * (2e0*v[1514] + v[1515]) + v[153] * v[1680];
//		v[4691] = v[1325] * v[1382] + v[1461] * (v[1516] - v[1472] * v[155]) + v[1505] * v[1695] - v[155] * v[1697]
//			+ v[1410] * v[457];
//		v[4692] = v[1324] * v[1382] + v[1381] * (v[1516] - v[1475] * v[155]) + v[1507] * v[1695] - v[155] * v[1698]
//			+ v[1409] * v[457];
//		v[4693] = v[1393] * (2e0*v[1516] + v[1517]) + v[155] * v[1680] + v[1508] * v[1695];
//		v[4694] = v[1325] * v[1379] + v[1461] * (v[1518] - v[1472] * v[157]) + v[1505] * v[1696] - v[157] * v[1697]
//			+ v[1410] * v[458];
//		v[4695] = v[1324] * v[1379] + v[1381] * (v[1518] - v[1475] * v[157]) + v[1507] * v[1696] - v[157] * v[1698]
//			+ v[1409] * v[458];
//		v[4696] = v[1393] * (2e0*v[1518] + v[1519]) + v[157] * v[1680] + v[1508] * v[1696];
//		v[1394] = v[1394] + v[1066] * v[1411] + v[1446] * v[945] + 2e0*v[216] * (v[1326] * v[1673] + (*a4)*v[4138 + i1319]
//			+ v[1411] * v[945]);
//		v[1397] = v[1397] + v[1060] * v[1439] + v[1674] * v[214] + 2e0*v[215] * (v[1060] * v[1410] + v[1325] * (v[1366] * v[182]
//			+ (v[1363] * v[185]) / 2e0 + (v[1362] * v[188]) / 2e0 - v[1367] * v[191] - (v[1364] * v[194]) / 2e0 - (v[1365] * v[197])
//			/ 2e0) + (*a4)*v[4264 + i1319]) + v[1410] * v[794];
//		v[1406] = v[1406] + v[1057] * v[1439] + v[1674] * v[215] + 2e0*v[214] * (v[1057] * v[1409] + v[1324] * (v[1366] * v[181]
//			+ (v[1363] * v[184]) / 2e0 + (v[1362] * v[187]) / 2e0 - v[1367] * v[190] - (v[1364] * v[193]) / 2e0 - (v[1365] * v[196])
//			/ 2e0) + (*a4)*v[4300 + i1319]) + v[1409] * v[794];
//		v[1376] = v[1376] + v[1451] * v[1452] * v[211] + (v[1334] * v[1369] + v[1335] * v[1371] + v[1336] * v[1372]
//			+ v[1406] * v[199] + v[1397] * v[200] + v[1394] * v[201])*v[212];
//		v[1454] = v[1336] * v[1663] + v[1394] * v[213] + (v[1337] * v[1372] + v[1376] * v[201]) / v[477];
//		v[1456] = v[1335] * v[1663] + v[1397] * v[213] + (v[1337] * v[1371] + v[1376] * v[200]) / v[477];
//		v[1458] = v[1334] * v[1663] + v[1406] * v[213] + (v[1337] * v[1369] + v[1376] * v[199]) / v[477];
//		v[4661] = v[1349] * v[1366] + v[125] * v[1458];
//		v[4662] = v[1348] * v[1366] + v[125] * v[1456];
//		v[4663] = v[1347] * v[1366] + v[125] * v[1454];
//		v[4664] = v[127] * v[1458] + v[1349] * v[1693];
//		v[4665] = v[127] * v[1456] + v[1348] * v[1693];
//		v[4666] = v[127] * v[1454] + v[1347] * v[1693];
//		v[4667] = v[129] * v[1458] + v[1349] * v[1694];
//		v[4668] = v[129] * v[1456] + v[1348] * v[1694];
//		v[4669] = v[129] * v[1454] + v[1347] * v[1694];
//		v[4670] = -(v[1349] * v[1367]) - v[1458] * v[153];
//		v[4671] = -(v[1348] * v[1367]) - v[1456] * v[153];
//		v[4672] = -(v[1347] * v[1367]) - v[1454] * v[153];
//		v[4673] = -(v[1458] * v[155]) + v[1349] * v[1695];
//		v[4674] = -(v[1456] * v[155]) + v[1348] * v[1695];
//		v[4675] = -(v[1454] * v[155]) + v[1347] * v[1695];
//		v[4676] = -(v[1458] * v[157]) + v[1349] * v[1696];
//		v[4677] = -(v[1456] * v[157]) + v[1348] * v[1696];
//		v[4678] = -(v[1454] * v[157]) + v[1347] * v[1696];
//		v[1466] = 2e0*v[1380] * v[1443] - v[1439] * v[1675] - v[1458] * v[169] - v[1456] * v[173] - v[1454] * v[177]
//			- v[1678] * v[190] - v[1679] * v[191] + v[1680] * v[192] + v[1393] * v[1943] + v[1461] * v[1944] + v[1381] * v[1945] + (*a4
//				)*(v[4444 + i1319] + v[216] * v[4462 + i1319]) + v[4480 + i1319];
//		v[1480] = -2e0*v[1380] * v[1440] + v[145] * v[1456] + v[141] * v[1458] + v[1454] * v[149] + v[1439] * v[1667]
//			+ v[1678] * v[181] + v[1679] * v[182] - v[1680] * v[183] + v[1393] * v[1949] + v[1461] * v[1950] + v[1381] * v[1951] + (*a4
//				)*(v[4606 + i1319] + v[216] * v[4624 + i1319]) + v[4642 + i1319];
//		v[1486] = (2e0*v[1380] * v[1445] - v[1466] - v[1439] * v[1683] - v[1458] * v[171] - v[1456] * v[175] - v[1454] * v[179]
//			+ v[1393] * v[1955] + v[1461] * v[1956] + v[1381] * v[1957] - v[1678] * v[196] - v[1679] * v[197] + v[1680] * v[198] + (*a4
//				)*(v[4336 + i1319] + v[216] * v[4354 + i1319]) + v[4372 + i1319]) / 2e0;
//		v[1492] = (2e0*v[1380] * v[1444] - v[1466] - v[1439] * v[1686] - v[1458] * v[170] - v[1456] * v[174] - v[1454] * v[178]
//			- v[1678] * v[193] - v[1679] * v[194] + v[1680] * v[195] + v[1393] * v[1961] + v[1461] * v[1962] + v[1381] * v[1963] + (*a4
//				)*(v[4390 + i1319] + v[216] * v[4408 + i1319]) + v[4426 + i1319]) / 2e0;
//		v[1498] = (-2e0*v[1380] * v[1441] + v[142] * v[1458] + v[1456] * v[146] - v[1480] + v[1454] * v[150] + v[1439] * v[1668]
//			+ v[1678] * v[184] + v[1679] * v[185] - v[1680] * v[186] + v[1393] * v[1967] + v[1461] * v[1968] + v[1381] * v[1969] + (*a4
//				)*(v[4552 + i1319] + v[216] * v[4570 + i1319]) + v[4588 + i1319]) / 2e0;
//		v[1504] = (-2e0*v[1380] * v[1442] + v[143] * v[1458] + v[1456] * v[147] - v[1480] + v[1454] * v[151] + v[1439] * v[1669]
//			+ v[1678] * v[187] + v[1679] * v[188] - v[1680] * v[189] + v[1393] * v[1973] + v[1461] * v[1974] + v[1381] * v[1975] + (*a4
//				)*(v[4498 + i1319] + v[216] * v[4516 + i1319]) + v[4534 + i1319]) / 2e0;
//		Rc[i1319 - 1] += v[1357] * v[1362] + v[1356] * v[1363] + v[1355] * v[1364] + v[1354] * v[1365] + v[4080 + i1319] + (*a4)*
//			(v[1391] * v[216] + v[4098 + i1319]);
//		for (i1360 = 1; i1360 <= 18; i1360++) {
//			Kc[i1319 - 1][i1360 - 1] += v[1504] * v[2822 + i1360] + v[1498] * v[2840 + i1360] + v[1492] * v[2858 + i1360]
//				+ v[1486] * v[2876 + i1360] + v[4660 + i1360] + (*a4)*v[4678 + i1360];
//		};/* end for */
//	};/* end for */
//#pragma endregion

#pragma region
double v01; double v010; double v011; double v012; double v013; double v014;
double v015; double v016; double v017; double v018; double v019; double v02;
double v020; double v021; double v022; double v023; double v024; double v025;
double v026; double v027; double v03; double v04; double v05; double v06; double v07;
double v08; double v09;
int i590, i623, i758, i1107, i1844, i1973, i2325, i2326, i2327, i2328, i2329, i2330, i2331
, i2332, i2333, i2334, i2335, i2336, i2337, i2338, i2339, i2340, i2341, i2342, b29, b30, b226
, b416, b483, b485, b502, b521
, b555, b595, b596, b600, b645, b646, b657, b763, b764, b765, b766
, b851, b855, b856, b873, b874, b990, b1018, b1071, b1173, b1254, b1300, b1430, b1431, b1440
, b1441, b1471, b1472, b1473, b1489, b1620, b1624, b1628, b1639, b1640, b1747, b1766, b1799
, b1853, b1854, b1858, b1863, b1864, b1952, b2019, b2105, b2106, b2114, b2118, b2122, b2132
, b2133, b2232;
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
v[480] = v[32] - v[33];
v[34] = (*n1);
v[2283] = v[31] * v[34];
v[495] = -1e0 + v[34];
v[1442] = -1e0 + v[495];
v[35] = (*n2);
v[500] = -1e0 + v[35];
v[1445] = -1e0 + v[500];
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
v[141] = v[36] + v[45];
v[235] = -v[141] / 2e0;
v[46] = u1A[1];
v[145] = v[37] + v[46];
v[237] = -v[145] / 2e0;
v[47] = u1A[2];
v[149] = v[38] + v[47];
v[239] = -v[149] / 2e0;
v[48] = u2A[0];
v[142] = v[39] + v[48];
v[231] = v[142] / 2e0 + v[235];
v[49] = u2A[1];
v[146] = v[40] + v[49];
v[232] = v[146] / 2e0 + v[237];
v[50] = u2A[2];
v[150] = v[41] + v[50];
v[233] = v[150] / 2e0 + v[239];
v[51] = u3A[0];
v[143] = v[42] + v[51];
v[236] = v[143] / 2e0 + v[235];
v[52] = u3A[1];
v[147] = v[43] + v[52];
v[238] = v[147] / 2e0 + v[237];
v[53] = u3A[2];
v[151] = v[44] + v[53];
v[240] = v[151] / 2e0 + v[239];
v[126] = -cAp[0] / 2e0;
v[128] = -cAp[1] / 2e0;
v[131] = -cAi[0] / 2e0;
v[133] = -cAi[1] / 2e0;
v[76] = xi1B[0];
v[77] = xi1B[1];
v[78] = xi1B[2];
v[79] = xi2B[0];
v[80] = xi2B[1];
v[81] = xi2B[2];
v[82] = xi3B[0];
v[83] = xi3B[1];
v[84] = xi3B[2];
v[85] = u1B[0];
v[169] = v[76] + v[85];
v[246] = -v[169] / 2e0;
v[86] = u1B[1];
v[173] = v[77] + v[86];
v[248] = -v[173] / 2e0;
v[87] = u1B[2];
v[177] = v[78] + v[87];
v[250] = -v[177] / 2e0;
v[88] = u2B[0];
v[170] = v[79] + v[88];
v[242] = v[170] / 2e0 + v[246];
v[89] = u2B[1];
v[174] = v[80] + v[89];
v[243] = v[174] / 2e0 + v[248];
v[90] = u2B[2];
v[178] = v[81] + v[90];
v[244] = v[178] / 2e0 + v[250];
v[91] = u3B[0];
v[171] = v[82] + v[91];
v[247] = v[171] / 2e0 + v[246];
v[92] = u3B[1];
v[175] = v[83] + v[92];
v[249] = v[175] / 2e0 + v[248];
v[93] = u3B[2];
v[179] = v[84] + v[93];
v[251] = v[179] / 2e0 + v[250];
v[154] = -cBp[0] / 2e0;
v[156] = -cBp[1] / 2e0;
v[159] = -cBi[0] / 2e0;
v[161] = -cBi[1] / 2e0;
v[125] = v[126] + v[128];
v[312] = -(v[125] * v[251]);
v[311] = -(v[125] * v[249]);
v[310] = -(v[125] * v[247]);
v[294] = -(v[125] * v[244]);
v[293] = -(v[125] * v[243]);
v[292] = -(v[125] * v[242]);
v[127] = 0.5e0 - v[126];
v[315] = -(v[127] * v[251]);
v[314] = -(v[127] * v[249]);
v[313] = -(v[127] * v[247]);
v[297] = -(v[127] * v[244]);
v[296] = -(v[127] * v[243]);
v[295] = -(v[127] * v[242]);
v[279] = v[127] * v[240];
v[278] = v[127] * v[238];
v[277] = v[127] * v[236];
v[129] = 0.5e0 - v[128];
v[318] = -(v[129] * v[251]);
v[317] = -(v[129] * v[249]);
v[316] = -(v[129] * v[247]);
v[300] = -(v[129] * v[244]);
v[299] = -(v[129] * v[243]);
v[298] = -(v[129] * v[242]);
v[264] = v[129] * v[233];
v[263] = v[129] * v[232];
v[262] = v[129] * v[231];
v[130] = v[131] + v[133];
v[132] = 0.5e0 - v[131];
v[134] = 0.5e0 - v[133];
v[153] = v[154] + v[156];
v[285] = -(v[153] * v[240]);
v[284] = -(v[153] * v[238]);
v[283] = -(v[153] * v[236]);
v[267] = -(v[153] * v[233]);
v[266] = -(v[153] * v[232]);
v[265] = -(v[153] * v[231]);
v[155] = 0.5e0 - v[154];
v[324] = v[155] * v[251];
v[323] = v[155] * v[249];
v[322] = v[155] * v[247];
v[288] = -(v[155] * v[240]);
v[287] = -(v[155] * v[238]);
v[286] = -(v[155] * v[236]);
v[270] = -(v[155] * v[233]);
v[269] = -(v[155] * v[232]);
v[268] = -(v[155] * v[231]);
v[157] = 0.5e0 - v[156];
v[3661] = 0e0;
v[3662] = 0e0;
v[3663] = v[125];
v[3664] = 0e0;
v[3665] = 0e0;
v[3666] = v[127];
v[3667] = 0e0;
v[3668] = 0e0;
v[3669] = v[129];
v[3670] = 0e0;
v[3671] = 0e0;
v[3672] = -v[153];
v[3673] = 0e0;
v[3674] = 0e0;
v[3675] = -v[155];
v[3676] = 0e0;
v[3677] = 0e0;
v[3678] = -v[157];
v[3643] = 0e0;
v[3644] = v[125];
v[3645] = 0e0;
v[3646] = 0e0;
v[3647] = v[127];
v[3648] = 0e0;
v[3649] = 0e0;
v[3650] = v[129];
v[3651] = 0e0;
v[3652] = 0e0;
v[3653] = -v[153];
v[3654] = 0e0;
v[3655] = 0e0;
v[3656] = -v[155];
v[3657] = 0e0;
v[3658] = 0e0;
v[3659] = -v[157];
v[3660] = 0e0;
v[3625] = v[125];
v[3626] = 0e0;
v[3627] = 0e0;
v[3628] = v[127];
v[3629] = 0e0;
v[3630] = 0e0;
v[3631] = v[129];
v[3632] = 0e0;
v[3633] = 0e0;
v[3634] = -v[153];
v[3635] = 0e0;
v[3636] = 0e0;
v[3637] = -v[155];
v[3638] = 0e0;
v[3639] = 0e0;
v[3640] = -v[157];
v[3641] = 0e0;
v[3642] = 0e0;
v[309] = v[157] * v[244];
v[308] = v[157] * v[243];
v[307] = v[157] * v[242];
v[291] = -(v[157] * v[240]);
v[290] = -(v[157] * v[238]);
v[289] = -(v[157] * v[236]);
v[273] = -(v[157] * v[233]);
v[272] = -(v[157] * v[232]);
v[271] = -(v[157] * v[231]);
v[158] = v[159] + v[161];
v[160] = 0.5e0 - v[159];
v[162] = 0.5e0 - v[161];
v[181] = dui1A[0] * v[11] + ddui1A[0] * v[12] + v[10] * v[45];
v[2268] = v[125] * v[181];
v[182] = dui1A[1] * v[11] + ddui1A[1] * v[12] + v[10] * v[46];
v[2262] = v[125] * v[182];
v[183] = dui1A[2] * v[11] + ddui1A[2] * v[12] + v[10] * v[47];
v[2274] = v[125] * v[183];
v[184] = dui2A[0] * v[11] + ddui2A[0] * v[12] + v[10] * v[48];
v[2267] = v[127] * v[184];
v[185] = dui2A[1] * v[11] + ddui2A[1] * v[12] + v[10] * v[49];
v[2261] = v[127] * v[185];
v[186] = dui2A[2] * v[11] + ddui2A[2] * v[12] + v[10] * v[50];
v[2273] = v[127] * v[186];
v[187] = dui3A[0] * v[11] + ddui3A[0] * v[12] + v[10] * v[51];
v[2266] = v[129] * v[187];
v[188] = dui3A[1] * v[11] + ddui3A[1] * v[12] + v[10] * v[52];
v[2260] = v[129] * v[188];
v[189] = dui3A[2] * v[11] + ddui3A[2] * v[12] + v[10] * v[53];
v[2272] = v[129] * v[189];
v[190] = dui1B[0] * v[11] + ddui1B[0] * v[12] + v[10] * v[85];
v[2265] = -(v[153] * v[190]);
v[191] = dui1B[1] * v[11] + ddui1B[1] * v[12] + v[10] * v[86];
v[2259] = -(v[153] * v[191]);
v[192] = dui1B[2] * v[11] + ddui1B[2] * v[12] + v[10] * v[87];
v[2271] = -(v[153] * v[192]);
v[193] = dui2B[0] * v[11] + ddui2B[0] * v[12] + v[10] * v[88];
v[2264] = -(v[155] * v[193]);
v[194] = dui2B[1] * v[11] + ddui2B[1] * v[12] + v[10] * v[89];
v[2258] = -(v[155] * v[194]);
v[195] = dui2B[2] * v[11] + ddui2B[2] * v[12] + v[10] * v[90];
v[2270] = -(v[155] * v[195]);
v[196] = dui3B[0] * v[11] + ddui3B[0] * v[12] + v[10] * v[91];
v[2263] = -(v[157] * v[196]);
v[1343] = v[2263] + v[2264] + v[2265] + v[2266] + v[2267] + v[2268];
v[197] = dui3B[1] * v[11] + ddui3B[1] * v[12] + v[10] * v[92];
v[2257] = -(v[157] * v[197]);
v[1358] = v[2257] + v[2258] + v[2259] + v[2260] + v[2261] + v[2262];
v[198] = dui3B[2] * v[11] + ddui3B[2] * v[12] + v[10] * v[93];
v[2269] = -(v[157] * v[198]);
v[1178] = v[2269] + v[2270] + v[2271] + v[2272] + v[2273] + v[2274];
v[199] = v[125] * v[141] + v[127] * v[142] + v[129] * v[143] - v[153] * v[169] - v[155] * v[170] - v[157] * v[171];
v[256] = -v[199] / 2e0;
v[325] = v[157] * v[247] + v[256];
v[319] = v[153] * v[247] - v[256];
v[304] = v[155] * v[242] + v[256];
v[301] = v[153] * v[242] - v[256];
v[280] = v[129] * v[236] - v[256];
v[274] = v[125] * v[236] + v[256];
v[257] = v[127] * v[231] - v[256];
v[253] = v[125] * v[231] + v[256];
v[200] = v[125] * v[145] + v[127] * v[146] + v[129] * v[147] - v[153] * v[173] - v[155] * v[174] - v[157] * v[175];
v[258] = -v[200] / 2e0;
v[326] = v[157] * v[249] + v[258];
v[320] = v[153] * v[249] - v[258];
v[305] = v[155] * v[243] + v[258];
v[302] = v[153] * v[243] - v[258];
v[281] = v[129] * v[238] - v[258];
v[275] = v[125] * v[238] + v[258];
v[259] = v[127] * v[232] - v[258];
v[254] = v[125] * v[232] + v[258];
v[201] = v[125] * v[149] + v[127] * v[150] + v[129] * v[151] - v[153] * v[177] - v[155] * v[178] - v[157] * v[179];
v[482] = sqrt((v[199] * v[199]) + (v[200] * v[200]) + (v[201] * v[201]));
v[636] = 1e0 / (v[482] * v[482]);
v[260] = -v[201] / 2e0;
v[327] = v[157] * v[251] + v[260];
v[321] = v[153] * v[251] - v[260];
v[306] = v[155] * v[244] + v[260];
v[303] = v[153] * v[244] - v[260];
v[282] = v[129] * v[240] - v[260];
v[276] = v[125] * v[240] + v[260];
v[261] = v[127] * v[233] - v[260];
v[255] = v[125] * v[233] + v[260];
v[202] = v[130] * v[36] + v[132] * v[39] + v[134] * v[42] - v[158] * v[76] - v[160] * v[79] - v[162] * v[82];
v[203] = v[130] * v[37] + v[132] * v[40] + v[134] * v[43] - v[158] * v[77] - v[160] * v[80] - v[162] * v[83];
v[204] = v[130] * v[38] + v[132] * v[41] + v[134] * v[44] - v[158] * v[78] - v[160] * v[81] - v[162] * v[84];
if (v[482] > 0.1e-7) { v01 = 1e0 / v[482]; v02 = (-(v01 / v[482])); v03 = (2e0*v01) / (v[482] * v[482]); }
else {
	v01 = (12500000e0 / 3e0)*(24e0 - (-0.1e-7 + v[482])*(0.24e10 - 2e0*(-1e0 + 100000000e0*v[482])*(0.2399999997e10
		- 0.1199999994e18*v[482] - 0.3e17*(v[482] * v[482]))));
	v02 = (-50000000e0 / 3e0)*(0.3599999994e10 - 0.4799999982e18*v[482] + 0.6e25*Power(v[482], 3)
		+ 0.1799999982e26*(v[482] * v[482]));
	v03 = 0.1e17*(799999997e0 - 0.599999994e17*v[482] - 0.3e17*(v[482] * v[482]));
};
v[211] = v03;
v[212] = v02;
v[213] = v01;
v[214] = v[199] * v[213];
v[215] = v[200] * v[213];
v[216] = v[201] * v[213];
v[217] = sqrt((v[202] * v[202]) + (v[203] * v[203]) + (v[204] * v[204]));
if (v[217] > 0.1e-7) { v04 = 1e0 / v[217]; v05 = (-(v04 / v[217])); v06 = (2e0*v04) / (v[217] * v[217]); }
else {
	v04 = (12500000e0 / 3e0)*(24e0 - (-0.1e-7 + v[217])*(0.24e10 - 2e0*(-1e0 + 100000000e0*v[217])*(0.2399999997e10
		- 0.1199999994e18*v[217] - 0.3e17*(v[217] * v[217]))));
	v05 = (-50000000e0 / 3e0)*(0.3599999994e10 - 0.4799999982e18*v[217] + 0.6e25*Power(v[217], 3)
		+ 0.1799999982e26*(v[217] * v[217]));
	v06 = 0.1e17*(799999997e0 - 0.599999994e17*v[217] - 0.3e17*(v[217] * v[217]));
};
v[222] = v04;
v[223] = v[202] * v[222];
v[224] = v[203] * v[222];
v[225] = v[204] * v[222];
b226 = v[214] * v[223] + v[215] * v[224] + v[216] * v[225] < 0e0;
if (b226) {
	v[228] = -v[214];
	v[229] = -v[215];
	v[230] = -v[216];
}
else {
	v[228] = v[214];
	v[229] = v[215];
	v[230] = v[216];
};
v[2426] = v[10] * v[230];
v[2344] = 2e0*v[230];
v[1415] = (v[2257] + v[2258] + v[2259] + v[2260] + v[2261] + v[2262])*v[230];
v[1422] = (v[2263] + v[2264] + v[2265] + v[2266] + v[2267] + v[2268])*v[230];
v[472] = (v[230] * v[230]);
v[956] = (v[2269] + v[2270] + v[2271] + v[2272] + v[2273] + v[2274])*v[230];
v[2347] = 2e0*v[229];
v[464] = (v[229] * v[229]);
v[2501] = v[181] * v[228] + v[182] * v[229];
v[2500] = v[184] * v[228] + v[185] * v[229];
v[2499] = v[187] * v[228] + v[188] * v[229];
v[2498] = -(v[190] * v[228]) - v[191] * v[229];
v[2497] = -(v[193] * v[228]) - v[194] * v[229];
v[2496] = -(v[196] * v[228]) - v[197] * v[229];
v[2350] = 2e0*v[228];
v[2275] = -(v[228] * v[229]);
v[1377] = v[1343] * v[228] + v[1358] * v[229];
v[463] = v[157] * v[2275];
v[462] = v[155] * v[2275];
v[461] = v[153] * v[2275];
v[460] = -(v[129] * v[2275]);
v[459] = -(v[127] * v[2275]);
v[458] = -(v[125] * v[2275]);
v[456] = (v[228] * v[228]);
v[344] = -(v[13] * v[253]) - v[14] * v[274] - v[15] * v[292] - v[16] * v[310];
v[345] = -(v[13] * v[254]) - v[14] * v[275] - v[15] * v[293] - v[16] * v[311];
v[346] = -(v[13] * v[255]) - v[14] * v[276] - v[15] * v[294] - v[16] * v[312];
v[347] = -(v[13] * v[257]) - v[14] * v[277] - v[15] * v[295] - v[16] * v[313];
v[348] = -(v[13] * v[259]) - v[14] * v[278] - v[15] * v[296] - v[16] * v[314];
v[349] = -(v[13] * v[261]) - v[14] * v[279] - v[15] * v[297] - v[16] * v[315];
v[350] = -(v[13] * v[262]) - v[14] * v[280] - v[15] * v[298] - v[16] * v[316];
v[351] = -(v[13] * v[263]) - v[14] * v[281] - v[15] * v[299] - v[16] * v[317];
v[352] = -(v[13] * v[264]) - v[14] * v[282] - v[15] * v[300] - v[16] * v[318];
v[353] = -(v[13] * v[265]) - v[14] * v[283] - v[15] * v[301] - v[16] * v[319];
v[354] = -(v[13] * v[266]) - v[14] * v[284] - v[15] * v[302] - v[16] * v[320];
v[355] = -(v[13] * v[267]) - v[14] * v[285] - v[15] * v[303] - v[16] * v[321];
v[356] = -(v[13] * v[268]) - v[14] * v[286] - v[15] * v[304] - v[16] * v[322];
v[357] = -(v[13] * v[269]) - v[14] * v[287] - v[15] * v[305] - v[16] * v[323];
v[358] = -(v[13] * v[270]) - v[14] * v[288] - v[15] * v[306] - v[16] * v[324];
v[359] = -(v[13] * v[271]) - v[14] * v[289] - v[15] * v[307] - v[16] * v[325];
v[360] = -(v[13] * v[272]) - v[14] * v[290] - v[15] * v[308] - v[16] * v[326];
v[361] = -(v[13] * v[273]) - v[14] * v[291] - v[15] * v[309] - v[16] * v[327];
v[818] = -(v[181] * v[344]) - v[182] * v[345] - v[183] * v[346] - v[184] * v[347] - v[185] * v[348] - v[186] * v[349]
- v[187] * v[350] - v[188] * v[351] - v[189] * v[352] - v[190] * v[353] - v[191] * v[354] - v[192] * v[355] - v[193] * v[356]
- v[194] * v[357] - v[195] * v[358] - v[196] * v[359] - v[197] * v[360] - v[198] * v[361];
v[3531] = v[344];
v[3532] = v[345];
v[3533] = v[346];
v[3534] = v[347];
v[3535] = v[348];
v[3536] = v[349];
v[3537] = v[350];
v[3538] = v[351];
v[3539] = v[352];
v[3540] = v[353];
v[3541] = v[354];
v[3542] = v[355];
v[3543] = v[356];
v[3544] = v[357];
v[3545] = v[358];
v[3546] = v[359];
v[3547] = v[360];
v[3548] = v[361];
v[362] = -(v[17] * v[253]) - v[18] * v[274] - v[19] * v[292] - v[20] * v[310];
v[363] = -(v[17] * v[254]) - v[18] * v[275] - v[19] * v[293] - v[20] * v[311];
v[364] = -(v[17] * v[255]) - v[18] * v[276] - v[19] * v[294] - v[20] * v[312];
v[365] = -(v[17] * v[257]) - v[18] * v[277] - v[19] * v[295] - v[20] * v[313];
v[366] = -(v[17] * v[259]) - v[18] * v[278] - v[19] * v[296] - v[20] * v[314];
v[367] = -(v[17] * v[261]) - v[18] * v[279] - v[19] * v[297] - v[20] * v[315];
v[368] = -(v[17] * v[262]) - v[18] * v[280] - v[19] * v[298] - v[20] * v[316];
v[369] = -(v[17] * v[263]) - v[18] * v[281] - v[19] * v[299] - v[20] * v[317];
v[370] = -(v[17] * v[264]) - v[18] * v[282] - v[19] * v[300] - v[20] * v[318];
v[371] = -(v[17] * v[265]) - v[18] * v[283] - v[19] * v[301] - v[20] * v[319];
v[372] = -(v[17] * v[266]) - v[18] * v[284] - v[19] * v[302] - v[20] * v[320];
v[373] = -(v[17] * v[267]) - v[18] * v[285] - v[19] * v[303] - v[20] * v[321];
v[374] = -(v[17] * v[268]) - v[18] * v[286] - v[19] * v[304] - v[20] * v[322];
v[375] = -(v[17] * v[269]) - v[18] * v[287] - v[19] * v[305] - v[20] * v[323];
v[376] = -(v[17] * v[270]) - v[18] * v[288] - v[19] * v[306] - v[20] * v[324];
v[377] = -(v[17] * v[271]) - v[18] * v[289] - v[19] * v[307] - v[20] * v[325];
v[378] = -(v[17] * v[272]) - v[18] * v[290] - v[19] * v[308] - v[20] * v[326];
v[379] = -(v[17] * v[273]) - v[18] * v[291] - v[19] * v[309] - v[20] * v[327];
v[816] = -(v[181] * v[362]) - v[182] * v[363] - v[183] * v[364] - v[184] * v[365] - v[185] * v[366] - v[186] * v[367]
- v[187] * v[368] - v[188] * v[369] - v[189] * v[370] - v[190] * v[371] - v[191] * v[372] - v[192] * v[373] - v[193] * v[374]
- v[194] * v[375] - v[195] * v[376] - v[196] * v[377] - v[197] * v[378] - v[198] * v[379];
v[3549] = v[362];
v[3550] = v[363];
v[3551] = v[364];
v[3552] = v[365];
v[3553] = v[366];
v[3554] = v[367];
v[3555] = v[368];
v[3556] = v[369];
v[3557] = v[370];
v[3558] = v[371];
v[3559] = v[372];
v[3560] = v[373];
v[3561] = v[374];
v[3562] = v[375];
v[3563] = v[376];
v[3564] = v[377];
v[3565] = v[378];
v[3566] = v[379];
v[380] = -(v[21] * v[253]) - v[22] * v[274] - v[23] * v[292] - v[24] * v[310];
v[381] = -(v[21] * v[254]) - v[22] * v[275] - v[23] * v[293] - v[24] * v[311];
v[382] = -(v[21] * v[255]) - v[22] * v[276] - v[23] * v[294] - v[24] * v[312];
v[383] = -(v[21] * v[257]) - v[22] * v[277] - v[23] * v[295] - v[24] * v[313];
v[384] = -(v[21] * v[259]) - v[22] * v[278] - v[23] * v[296] - v[24] * v[314];
v[385] = -(v[21] * v[261]) - v[22] * v[279] - v[23] * v[297] - v[24] * v[315];
v[386] = -(v[21] * v[262]) - v[22] * v[280] - v[23] * v[298] - v[24] * v[316];
v[387] = -(v[21] * v[263]) - v[22] * v[281] - v[23] * v[299] - v[24] * v[317];
v[388] = -(v[21] * v[264]) - v[22] * v[282] - v[23] * v[300] - v[24] * v[318];
v[389] = -(v[21] * v[265]) - v[22] * v[283] - v[23] * v[301] - v[24] * v[319];
v[390] = -(v[21] * v[266]) - v[22] * v[284] - v[23] * v[302] - v[24] * v[320];
v[391] = -(v[21] * v[267]) - v[22] * v[285] - v[23] * v[303] - v[24] * v[321];
v[392] = -(v[21] * v[268]) - v[22] * v[286] - v[23] * v[304] - v[24] * v[322];
v[393] = -(v[21] * v[269]) - v[22] * v[287] - v[23] * v[305] - v[24] * v[323];
v[394] = -(v[21] * v[270]) - v[22] * v[288] - v[23] * v[306] - v[24] * v[324];
v[395] = -(v[21] * v[271]) - v[22] * v[289] - v[23] * v[307] - v[24] * v[325];
v[396] = -(v[21] * v[272]) - v[22] * v[290] - v[23] * v[308] - v[24] * v[326];
v[397] = -(v[21] * v[273]) - v[22] * v[291] - v[23] * v[309] - v[24] * v[327];
v[814] = v[181] * v[380] + v[182] * v[381] + v[183] * v[382] + v[184] * v[383] + v[185] * v[384] + v[186] * v[385]
+ v[187] * v[386] + v[188] * v[387] + v[189] * v[388] + v[190] * v[389] + v[191] * v[390] + v[192] * v[391] + v[193] * v[392]
+ v[194] * v[393] + v[195] * v[394] + v[196] * v[395] + v[197] * v[396] + v[198] * v[397];
v[3567] = v[380];
v[3568] = v[381];
v[3569] = v[382];
v[3570] = v[383];
v[3571] = v[384];
v[3572] = v[385];
v[3573] = v[386];
v[3574] = v[387];
v[3575] = v[388];
v[3576] = v[389];
v[3577] = v[390];
v[3578] = v[391];
v[3579] = v[392];
v[3580] = v[393];
v[3581] = v[394];
v[3582] = v[395];
v[3583] = v[396];
v[3584] = v[397];
v[398] = -(v[25] * v[253]) - v[26] * v[274] - v[27] * v[292] - v[28] * v[310];
v[687] = v[233] * v[344] + v[240] * v[362] - v[244] * v[380] - v[251] * v[398];
v[686] = v[232] * v[344] + v[238] * v[362] - v[243] * v[380] - v[249] * v[398];
v[685] = v[231] * v[344] + v[236] * v[362] - v[242] * v[380] - v[247] * v[398];
v[399] = -(v[25] * v[254]) - v[26] * v[275] - v[27] * v[293] - v[28] * v[311];
v[691] = v[233] * v[345] + v[240] * v[363] - v[244] * v[381] - v[251] * v[399];
v[690] = v[232] * v[345] + v[238] * v[363] - v[243] * v[381] - v[249] * v[399];
v[689] = v[231] * v[345] + v[236] * v[363] - v[242] * v[381] - v[247] * v[399];
v[400] = -(v[25] * v[255]) - v[26] * v[276] - v[27] * v[294] - v[28] * v[312];
v[695] = v[233] * v[346] + v[240] * v[364] - v[244] * v[382] - v[251] * v[400];
v[694] = v[232] * v[346] + v[238] * v[364] - v[243] * v[382] - v[249] * v[400];
v[693] = v[231] * v[346] + v[236] * v[364] - v[242] * v[382] - v[247] * v[400];
v[401] = -(v[25] * v[257]) - v[26] * v[277] - v[27] * v[295] - v[28] * v[313];
v[699] = v[233] * v[347] + v[240] * v[365] - v[244] * v[383] - v[251] * v[401];
v[698] = v[232] * v[347] + v[238] * v[365] - v[243] * v[383] - v[249] * v[401];
v[697] = v[231] * v[347] + v[236] * v[365] - v[242] * v[383] - v[247] * v[401];
v[402] = -(v[25] * v[259]) - v[26] * v[278] - v[27] * v[296] - v[28] * v[314];
v[703] = v[233] * v[348] + v[240] * v[366] - v[244] * v[384] - v[251] * v[402];
v[702] = v[232] * v[348] + v[238] * v[366] - v[243] * v[384] - v[249] * v[402];
v[701] = v[231] * v[348] + v[236] * v[366] - v[242] * v[384] - v[247] * v[402];
v[403] = -(v[25] * v[261]) - v[26] * v[279] - v[27] * v[297] - v[28] * v[315];
v[707] = v[233] * v[349] + v[240] * v[367] - v[244] * v[385] - v[251] * v[403];
v[706] = v[232] * v[349] + v[238] * v[367] - v[243] * v[385] - v[249] * v[403];
v[705] = v[231] * v[349] + v[236] * v[367] - v[242] * v[385] - v[247] * v[403];
v[404] = -(v[25] * v[262]) - v[26] * v[280] - v[27] * v[298] - v[28] * v[316];
v[711] = v[233] * v[350] + v[240] * v[368] - v[244] * v[386] - v[251] * v[404];
v[710] = v[232] * v[350] + v[238] * v[368] - v[243] * v[386] - v[249] * v[404];
v[709] = v[231] * v[350] + v[236] * v[368] - v[242] * v[386] - v[247] * v[404];
v[405] = -(v[25] * v[263]) - v[26] * v[281] - v[27] * v[299] - v[28] * v[317];
v[715] = v[233] * v[351] + v[240] * v[369] - v[244] * v[387] - v[251] * v[405];
v[714] = v[232] * v[351] + v[238] * v[369] - v[243] * v[387] - v[249] * v[405];
v[713] = v[231] * v[351] + v[236] * v[369] - v[242] * v[387] - v[247] * v[405];
v[406] = -(v[25] * v[264]) - v[26] * v[282] - v[27] * v[300] - v[28] * v[318];
v[719] = v[233] * v[352] + v[240] * v[370] - v[244] * v[388] - v[251] * v[406];
v[718] = v[232] * v[352] + v[238] * v[370] - v[243] * v[388] - v[249] * v[406];
v[717] = v[231] * v[352] + v[236] * v[370] - v[242] * v[388] - v[247] * v[406];
v[407] = -(v[25] * v[265]) - v[26] * v[283] - v[27] * v[301] - v[28] * v[319];
v[723] = v[233] * v[353] + v[240] * v[371] - v[244] * v[389] - v[251] * v[407];
v[722] = v[232] * v[353] + v[238] * v[371] - v[243] * v[389] - v[249] * v[407];
v[721] = v[231] * v[353] + v[236] * v[371] - v[242] * v[389] - v[247] * v[407];
v[408] = -(v[25] * v[266]) - v[26] * v[284] - v[27] * v[302] - v[28] * v[320];
v[727] = v[233] * v[354] + v[240] * v[372] - v[244] * v[390] - v[251] * v[408];
v[726] = v[232] * v[354] + v[238] * v[372] - v[243] * v[390] - v[249] * v[408];
v[725] = v[231] * v[354] + v[236] * v[372] - v[242] * v[390] - v[247] * v[408];
v[409] = -(v[25] * v[267]) - v[26] * v[285] - v[27] * v[303] - v[28] * v[321];
v[731] = v[233] * v[355] + v[240] * v[373] - v[244] * v[391] - v[251] * v[409];
v[730] = v[232] * v[355] + v[238] * v[373] - v[243] * v[391] - v[249] * v[409];
v[729] = v[231] * v[355] + v[236] * v[373] - v[242] * v[391] - v[247] * v[409];
v[410] = -(v[25] * v[268]) - v[26] * v[286] - v[27] * v[304] - v[28] * v[322];
v[735] = v[233] * v[356] + v[240] * v[374] - v[244] * v[392] - v[251] * v[410];
v[734] = v[232] * v[356] + v[238] * v[374] - v[243] * v[392] - v[249] * v[410];
v[733] = v[231] * v[356] + v[236] * v[374] - v[242] * v[392] - v[247] * v[410];
v[411] = -(v[25] * v[269]) - v[26] * v[287] - v[27] * v[305] - v[28] * v[323];
v[739] = v[233] * v[357] + v[240] * v[375] - v[244] * v[393] - v[251] * v[411];
v[738] = v[232] * v[357] + v[238] * v[375] - v[243] * v[393] - v[249] * v[411];
v[737] = v[231] * v[357] + v[236] * v[375] - v[242] * v[393] - v[247] * v[411];
v[412] = -(v[25] * v[270]) - v[26] * v[288] - v[27] * v[306] - v[28] * v[324];
v[743] = v[233] * v[358] + v[240] * v[376] - v[244] * v[394] - v[251] * v[412];
v[742] = v[232] * v[358] + v[238] * v[376] - v[243] * v[394] - v[249] * v[412];
v[741] = v[231] * v[358] + v[236] * v[376] - v[242] * v[394] - v[247] * v[412];
v[413] = -(v[25] * v[271]) - v[26] * v[289] - v[27] * v[307] - v[28] * v[325];
v[747] = v[233] * v[359] + v[240] * v[377] - v[244] * v[395] - v[251] * v[413];
v[746] = v[232] * v[359] + v[238] * v[377] - v[243] * v[395] - v[249] * v[413];
v[745] = v[231] * v[359] + v[236] * v[377] - v[242] * v[395] - v[247] * v[413];
v[414] = -(v[25] * v[272]) - v[26] * v[290] - v[27] * v[308] - v[28] * v[326];
v[751] = v[233] * v[360] + v[240] * v[378] - v[244] * v[396] - v[251] * v[414];
v[750] = v[232] * v[360] + v[238] * v[378] - v[243] * v[396] - v[249] * v[414];
v[749] = v[231] * v[360] + v[236] * v[378] - v[242] * v[396] - v[247] * v[414];
v[415] = -(v[25] * v[273]) - v[26] * v[291] - v[27] * v[309] - v[28] * v[327];
v[812] = v[181] * v[398] + v[182] * v[399] + v[183] * v[400] + v[184] * v[401] + v[185] * v[402] + v[186] * v[403]
+ v[187] * v[404] + v[188] * v[405] + v[189] * v[406] + v[190] * v[407] + v[191] * v[408] + v[192] * v[409] + v[193] * v[410]
+ v[194] * v[411] + v[195] * v[412] + v[196] * v[413] + v[197] * v[414] + v[198] * v[415];
v[755] = v[233] * v[361] + v[240] * v[379] - v[244] * v[397] - v[251] * v[415];
v[754] = v[232] * v[361] + v[238] * v[379] - v[243] * v[397] - v[249] * v[415];
v[753] = v[231] * v[361] + v[236] * v[379] - v[242] * v[397] - v[247] * v[415];
v[3585] = v[398];
v[3586] = v[399];
v[3587] = v[400];
v[3588] = v[401];
v[3589] = v[402];
v[3590] = v[403];
v[3591] = v[404];
v[3592] = v[405];
v[3593] = v[406];
v[3594] = v[407];
v[3595] = v[408];
v[3596] = v[409];
v[3597] = v[410];
v[3598] = v[411];
v[3599] = v[412];
v[3600] = v[413];
v[3601] = v[414];
v[3602] = v[415];
b416 = sqrt(Power(-(v[224] * v[228]) + v[223] * v[229], 2) + Power(v[225] * v[228] - v[223] * v[230], 2) + Power(-
(v[225] * v[229]) + v[224] * v[230], 2)) > 0.1e-7;
if (b416) {
	v[418] = -(v[225] * v[229]) + v[224] * v[230];
	v[419] = v[225] * v[228] - v[223] * v[230];
	v[420] = -(v[224] * v[228]) + v[223] * v[229];
	v[421] = sqrt((v[418] * v[418]) + (v[419] * v[419]) + (v[420] * v[420]));
	v[1262] = 1e0 / (v[421] * v[421]);
	v[1065] = v[421];
	v[1273] = 1e0 - (v[1065] * v[1065]);
	v[2488] = 1e0 / Power(v[1273], 0.15e1);
	v[1268] = 1e0 / sqrt(v[1273]);
	v[1064] = asin(v[1065]) / 2e0;
	v[1267] = 1e0 / Power(cos(v[1064]), 2);
	v[2351] = v[1267] * v[1268];
	v[423] = 2e0*tan(v[1064]);
	if (v[421] > 0.1e-7) { v07 = 1e0 / v[421]; v08 = (-(v07 / v[421])); v09 = (2e0*v07) / (v[421] * v[421]); }
	else {
		v07 = (12500000e0 / 3e0)*(24e0 - (-0.1e-7 + v[421])*(0.24e10 - 2e0*(-1e0 + 100000000e0*v[421])*
			(0.2399999997e10 - 0.1199999994e18*v[421] - 0.3e17*(v[421] * v[421]))));
		v08 = (-50000000e0 / 3e0)*(0.3599999994e10 - 0.4799999982e18*v[421] + 0.6e25*Power(v[421], 3)
			+ 0.1799999982e26*(v[421] * v[421]));
		v09 = 0.1e17*(799999997e0 - 0.599999994e17*v[421] - 0.3e17*(v[421] * v[421]));
	};
	v[427] = v09;
	v[428] = v08;
	v[429] = v07;
	v[2487] = v[423] * v[428] + v[2351] * v[429];
	v[2276] = v[423] * v[429];
	v[430] = v[2276] * v[418];
	v[2388] = 2e0*v[430];
	v[2306] = v[430] / 2e0;
	v[441] = (v[430] * v[430]);
	v[431] = v[2276] * v[419];
	v[2277] = v[431] / 2e0;
	v[439] = v[2277] * v[430];
	v[434] = (v[431] * v[431]);
	v[1019] = -v[434] - v[441];
	v[432] = v[2276] * v[420];
	v[2385] = 2e0*v[432];
	v[1049] = -v[432] + v[439];
	v[1040] = v[432] + v[439];
	v[446] = v[2277] * v[432];
	v[1031] = -v[430] + v[446];
	v[1023] = v[430] + v[446];
	v[444] = v[2306] * v[432];
	v[1044] = v[431] + v[444];
	v[1027] = -v[431] + v[444];
	v[435] = (v[432] * v[432]);
	v[1059] = 4e0 + v[434] + v[435] + v[441];
	v[2489] = 1e0 / Power(v[1059], 3);
	v[2384] = -4e0 / (v[1059] * v[1059]);
	v[1054] = -v[434] - v[435];
	v[1036] = -v[435] - v[441];
	v[433] = 4e0 / v[1059];
	v[2278] = v[433] / 2e0;
	v[436] = 1e0 + v[1054] * v[2278];
	v[437] = v[1049] * v[433];
	v[438] = v[1044] * v[433];
	v[440] = v[1040] * v[433];
	v[442] = 1e0 + v[1036] * v[2278];
	v[443] = v[1031] * v[433];
	v[445] = v[1027] * v[433];
	v[447] = v[1023] * v[433];
	v[448] = 1e0 + v[1019] * v[2278];
}
else {
	v[436] = 1e0;
	v[437] = 0e0;
	v[438] = 0e0;
	v[440] = 0e0;
	v[442] = 1e0;
	v[443] = 0e0;
	v[445] = 0e0;
	v[447] = 0e0;
	v[448] = 1e0;
};
if (b30) {
	v[995] = 1e0 - v[472];
	v[993] = 1e0 - v[464];
	v[991] = 1e0 - v[456];
	v[453] = v[130] * v[149] + v[132] * v[150] + v[134] * v[151] - v[158] * v[177] - v[160] * v[178] - v[162] * v[179]
		+ v[1] * v[445] + v[2] * v[447] + v[3] * v[448];
	v[2279] = v[230] * v[453];
	v[452] = v[130] * v[145] + v[132] * v[146] + v[134] * v[147] - v[158] * v[173] - v[160] * v[174] - v[162] * v[175]
		+ v[1] * v[440] + v[2] * v[442] + v[3] * v[443];
	v[2281] = v[229] * v[452];
	v[2303] = v[2279] + v[2281];
	v[451] = v[130] * v[141] + v[132] * v[142] + v[134] * v[143] - v[158] * v[169] - v[160] * v[170] - v[162] * v[171]
		+ v[1] * v[436] + v[2] * v[437] + v[3] * v[438];
	v[2280] = -(v[228] * v[451]);
	v[2305] = v[2280] - v[2281];
	v[2304] = -v[2279] + v[2280];
	v[450] = -(v[228] * v[2303]) + v[451] * v[991];
	v[454] = v[229] * v[2304] + v[452] * v[993];
	v[455] = v[230] * v[2305] + v[453] * v[995];
}
else {
	v[450] = 0e0;
	v[454] = 0e0;
	v[455] = 0e0;
};
v[457] = v[1343] * v[456] + v[182] * v[458] + v[185] * v[459] + v[188] * v[460] + v[191] * v[461] + v[194] * v[462]
+ v[197] * v[463] + v[228] * v[956];
v[471] = v[181] * v[458] + v[184] * v[459] + v[187] * v[460] + v[190] * v[461] + v[193] * v[462] + v[196] * v[463]
+ v[1358] * v[464] + v[229] * v[956];
v[473] = v[1422] * v[228] + v[1415] * v[229] + v[1178] * v[472];
(*vnrel) = sqrt((v[457] * v[457]) + (v[471] * v[471]) + (v[473] * v[473]));
v[479] = -((v[2283] * Power(v[480], v[495])) / (v[35] * Power(v[33], v[500])));
v[2291] = -(v[35] * v[479]);
b483 = v[482] < v[32];
if (b483) {
	b485 = v[482] > v[33];
	if (b485) {
		v[494] = v[32] - v[482];
		v[2282] = -(v[31] * Power(v[494], v[34]));
		v[487] = v[228] * v[2282];
		v[489] = v[2282] * v[229];
		v[490] = v[2282] * v[230];
	}
	else {
		v[491] = -(v[31] * Power(v[480], v[34])) + v[479] * (Power(v[33], v[35]) - Power(v[482], v[35]));
		v[487] = v[228] * v[491];
		v[489] = v[229] * v[491];
		v[490] = v[230] * v[491];
	};
}
else {
	v[487] = 0e0;
	v[489] = 0e0;
	v[490] = 0e0;
};
if (b483) {
	v[2284] = 2e0*v[7];
	if (b485) {
		v[860] = v[2283] * v[8];
		v[861] = sqrt(v[860] * Power(v[494], v[495]));
		v[497] = v[2284] * v[861];
		v[496] = v[457] * v[497];
		v[498] = v[471] * v[497];
		v[499] = v[473] * v[497];
	}
	else {
		v[868] = v[2291] * v[8];
		v[869] = sqrt(v[868] * Power(v[482], v[500]));
		v[501] = v[2284] * v[869];
		v[496] = v[457] * v[501];
		v[498] = v[471] * v[501];
		v[499] = v[473] * v[501];
	};
	b502 = v[482] < v[32] && v[228] * (v[487] + v[496]) + v[229] * (v[489] + v[498]) + v[230] * (v[490] + v[499]) < 0e0;
	if (b502) {
		v[504] = v[496];
		v[505] = v[498];
		v[506] = v[499];
	}
	else {
		v[504] = -v[487];
		v[505] = -v[489];
		v[506] = -v[490];
	};
}
else {
	v[504] = 0e0;
	v[505] = 0e0;
	v[506] = 0e0;
};
v[507] = v[487] + v[504];
v[508] = v[489] + v[505];
v[509] = v[490] + v[506];
v[1487] = (v[507] * v[507]) + (v[508] * v[508]) + (v[509] * v[509]);
v[510] = v[4] * v[450];
v[511] = v[4] * v[454];
v[512] = v[4] * v[455];
v[516] = v[510] - (v[181] * v[685] + v[182] * v[689] + v[183] * v[693] + v[184] * v[697] + v[185] * v[701] + v[186] * v[705]
	+ v[187] * v[709] + v[188] * v[713] + v[189] * v[717] + v[190] * v[721] + v[191] * v[725] + v[192] * v[729] + v[193] * v[733]
	+ v[194] * v[737] + v[195] * v[741] + v[196] * v[745] + v[197] * v[749] + v[198] * v[753])*v[9];
v[517] = v[511] - (v[181] * v[686] + v[182] * v[690] + v[183] * v[694] + v[184] * v[698] + v[185] * v[702] + v[186] * v[706]
	+ v[187] * v[710] + v[188] * v[714] + v[189] * v[718] + v[190] * v[722] + v[191] * v[726] + v[192] * v[730] + v[193] * v[734]
	+ v[194] * v[738] + v[195] * v[742] + v[196] * v[746] + v[197] * v[750] + v[198] * v[754])*v[9];
v[518] = v[512] - (v[181] * v[687] + v[182] * v[691] + v[183] * v[695] + v[184] * v[699] + v[185] * v[703] + v[186] * v[707]
	+ v[187] * v[711] + v[188] * v[715] + v[189] * v[719] + v[190] * v[723] + v[191] * v[727] + v[192] * v[731] + v[193] * v[735]
	+ v[194] * v[739] + v[195] * v[743] + v[196] * v[747] + v[197] * v[751] + v[198] * v[755])*v[9];
v[1483] = (v[516] * v[516]) + (v[517] * v[517]) + (v[518] * v[518]);
if (b483) {
	if (b29) {
		b521 = sqrt((v[516] * v[516]) + (v[517] * v[517]) + (v[518] * v[518])) <= v[5] * sqrt((v[507] * v[507]) +
			(v[508] * v[508]) + (v[509] * v[509]));
		if (b521) {
			v[523] = v[516];
			v[524] = v[517];
			v[525] = v[518];
			v[526] = 1e0;
		}
		else {
			v[2285] = v[6] * sqrt(v[1487]);
			v[527] = sqrt(v[1483]);
			if (v[527] > 0.1e-5) { v010 = 1e0 / v[527]; v011 = (-(v010 / v[527])); v012 = (2e0*v010) / (v[527] * v[527]); }
			else {
				v010 = (24000000e0 - (-1e0 + 1000000e0*v[527])*(71999994e0 - 0.71999982e14*v[527] + 0.6e19*Power(v[527], 3)
					+ 0.23999982e20*(v[527] * v[527]))) / 24e0;
				v011 = (-500000e0 / 3e0)*(35999994e0 - 0.47999982e14*v[527] + 0.6e19*Power(v[527], 3) + 0.17999982e20*
					(v[527] * v[527]));
				v012 = 0.1e13*(7999997e0 - 0.5999994e13*v[527] - 0.3e13*(v[527] * v[527]));
			};
			v[531] = v011;
			v[532] = v010;
			v[533] = v[516] * v[532];
			v[534] = v[517] * v[532];
			v[535] = v[518] * v[532];
			v[523] = v[2285] * v[533];
			v[524] = v[2285] * v[534];
			v[525] = v[2285] * v[535];
			v[526] = 0e0;
		};
		if (sqrt((v[510] * v[510]) + (v[511] * v[511]) + (v[512] * v[512])) > v[5] * sqrt((v[507] * v[507]) + (v[508] * v[508]
			) + (v[509] * v[509]))) {
			if (v[4] > 0.1e-5) { v013 = 1e0 / v[4]; v014 = (-(v013 / v[4])); v015 = (2e0*v013) / (v[4] * v[4]); }
			else {
				v013 = (24000000e0 - (-1e0 + 1000000e0*v[4])*(71999994e0 - 0.71999982e14*v[4] + 0.6e19*Power(v[4], 3)
					+ 0.23999982e20*(v[4] * v[4]))) / 24e0;
				v014 = (-500000e0 / 3e0)*(35999994e0 - 0.47999982e14*v[4] + 0.6e19*Power(v[4], 3) + 0.17999982e20*
					(v[4] * v[4]));
				v015 = 0.1e13*(7999997e0 - 0.5999994e13*v[4] - 0.3e13*(v[4] * v[4]));
			};
			v[545] = sqrt((v[510] * v[510]) + (v[511] * v[511]) + (v[512] * v[512]));
			if (v[545] > 0.1e-5) { v016 = 1e0 / v[545]; v017 = (-(v016 / v[545])); v018 = (2e0*v016) / (v[545] * v[545]); }
			else {
				v016 = (24000000e0 - (-1e0 + 1000000e0*v[545])*(71999994e0 - 0.71999982e14*v[545] + 0.6e19*Power(v[545], 3)
					+ 0.23999982e20*(v[545] * v[545]))) / 24e0;
				v017 = (-500000e0 / 3e0)*(35999994e0 - 0.47999982e14*v[545] + 0.6e19*Power(v[545], 3) + 0.17999982e20*
					(v[545] * v[545]));
				v018 = 0.1e13*(7999997e0 - 0.5999994e13*v[545] - 0.3e13*(v[545] * v[545]));
			};
			v[552] = -(v013*v016*v[6] * sqrt(v[1487]));
			v[551] = v[450] + v[510] * v[552];
			v[553] = v[454] + v[511] * v[552];
			v[554] = v[455] + v[512] * v[552];
		}
		else {
			v[551] = 0e0;
			v[553] = 0e0;
			v[554] = 0e0;
		};
	}
	else {
		b555 = sqrt((v[516] * v[516]) + (v[517] * v[517]) + (v[518] * v[518])) <= v[6] * sqrt((v[507] * v[507]) +
			(v[508] * v[508]) + (v[509] * v[509]));
		if (b555) {
			v[523] = v[516];
			v[524] = v[517];
			v[525] = v[518];
			v[526] = 1e0;
		}
		else {
			v[566] = sqrt(v[1487]);
			v[2286] = v[566] * v[6];
			v[557] = sqrt(v[1483]);
			if (v[557] > 0.1e-5) { v019 = 1e0 / v[557]; v020 = (-(v019 / v[557])); v021 = (2e0*v019) / (v[557] * v[557]); }
			else {
				v019 = (24000000e0 - (-1e0 + 1000000e0*v[557])*(71999994e0 - 0.71999982e14*v[557] + 0.6e19*Power(v[557], 3)
					+ 0.23999982e20*(v[557] * v[557]))) / 24e0;
				v020 = (-500000e0 / 3e0)*(35999994e0 - 0.47999982e14*v[557] + 0.6e19*Power(v[557], 3) + 0.17999982e20*
					(v[557] * v[557]));
				v021 = 0.1e13*(7999997e0 - 0.5999994e13*v[557] - 0.3e13*(v[557] * v[557]));
			};
			v[561] = v020;
			v[562] = v019;
			v[563] = v[516] * v[562];
			v[564] = v[517] * v[562];
			v[565] = v[518] * v[562];
			v[523] = v[2286] * v[563];
			v[524] = v[2286] * v[564];
			v[525] = v[2286] * v[565];
			v[526] = 0e0;
		};
		if (sqrt((v[510] * v[510]) + (v[511] * v[511]) + (v[512] * v[512])) > v[6] * sqrt((v[507] * v[507]) + (v[508] * v[508]
			) + (v[509] * v[509]))) {
			if (v[4] > 0.1e-5) { v022 = 1e0 / v[4]; v023 = (-(v022 / v[4])); v024 = (2e0*v022) / (v[4] * v[4]); }
			else {
				v022 = (24000000e0 - (-1e0 + 1000000e0*v[4])*(71999994e0 - 0.71999982e14*v[4] + 0.6e19*Power(v[4], 3)
					+ 0.23999982e20*(v[4] * v[4]))) / 24e0;
				v023 = (-500000e0 / 3e0)*(35999994e0 - 0.47999982e14*v[4] + 0.6e19*Power(v[4], 3) + 0.17999982e20*
					(v[4] * v[4]));
				v024 = 0.1e13*(7999997e0 - 0.5999994e13*v[4] - 0.3e13*(v[4] * v[4]));
			};
			v[575] = sqrt((v[510] * v[510]) + (v[511] * v[511]) + (v[512] * v[512]));
			if (v[575] > 0.1e-5) { v025 = 1e0 / v[575]; v026 = (-(v025 / v[575])); v027 = (2e0*v025) / (v[575] * v[575]); }
			else {
				v025 = (24000000e0 - (-1e0 + 1000000e0*v[575])*(71999994e0 - 0.71999982e14*v[575] + 0.6e19*Power(v[575], 3)
					+ 0.23999982e20*(v[575] * v[575]))) / 24e0;
				v026 = (-500000e0 / 3e0)*(35999994e0 - 0.47999982e14*v[575] + 0.6e19*Power(v[575], 3) + 0.17999982e20*
					(v[575] * v[575]));
				v027 = 0.1e13*(7999997e0 - 0.5999994e13*v[575] - 0.3e13*(v[575] * v[575]));
			};
			v[581] = -(v022*v025*v[6] * sqrt(v[1487]));
			v[551] = v[450] + v[510] * v[581];
			v[553] = v[454] + v[511] * v[581];
			v[554] = v[455] + v[512] * v[581];
		}
		else {
			v[551] = 0e0;
			v[553] = 0e0;
			v[554] = 0e0;
		};
	};
}
else {
	v[523] = 0e0;
	v[524] = 0e0;
	v[525] = 0e0;
};
fn[0] = v[507];
fn[1] = v[508];
fn[2] = v[509];
ft[0] = v[523];
ft[1] = v[524];
ft[2] = v[525];
(*stickupdated) = v[526];
gtpupdated[0] = v[450] - v[551];
gtpupdated[1] = v[454] - v[553];
gtpupdated[2] = v[455] - v[554];
b595 = b483;
if (b595) {
	b596 = b485;
}
else {
};
b600 = b226;
if (b600) {
	v[601] = 0e0;
	v[602] = 0e0;
	v[603] = 0e0;
}
else {
	v[603] = 0e0;
	v[602] = 0e0;
	v[601] = 0e0;
};
v[608] = v[201] * v[601] + v[200] * v[602] + v[199] * v[603];
v[610] = v[212] * v[608];
v[2287] = v[610] / v[482];
v[611] = v[201] * v[2287] + v[490] + v[213] * v[601];
v[2297] = v[611] / 2e0;
v[613] = v[200] * v[2287] + v[489] + v[213] * v[602];
v[2296] = v[613] / 2e0;
v[614] = v[199] * v[2287] + v[487] + v[213] * v[603];
v[3513] = v[125] * v[614];
v[3514] = v[125] * v[613];
v[3515] = v[125] * v[611];
v[3516] = v[127] * v[614];
v[3517] = v[127] * v[613];
v[3518] = v[127] * v[611];
v[3519] = v[129] * v[614];
v[3520] = v[129] * v[613];
v[3521] = v[129] * v[611];
v[3522] = -(v[153] * v[614]);
v[3523] = -(v[153] * v[613]);
v[3524] = -(v[153] * v[611]);
v[3525] = -(v[155] * v[614]);
v[3526] = -(v[155] * v[613]);
v[3527] = -(v[155] * v[611]);
v[3528] = -(v[157] * v[614]);
v[3529] = -(v[157] * v[613]);
v[3530] = -(v[157] * v[611]);
v[3751] = 0e0;
v[3752] = 0e0;
v[3753] = 0e0;
v[3754] = v[614];
v[3755] = v[613];
v[3756] = v[611];
v[3757] = 0e0;
v[3758] = 0e0;
v[3759] = 0e0;
v[3760] = 0e0;
v[3761] = 0e0;
v[3762] = 0e0;
v[3763] = 0e0;
v[3764] = 0e0;
v[3765] = 0e0;
v[3766] = 0e0;
v[3767] = 0e0;
v[3768] = 0e0;
v[3733] = 0e0;
v[3734] = 0e0;
v[3735] = 0e0;
v[3736] = 0e0;
v[3737] = 0e0;
v[3738] = 0e0;
v[3739] = v[614];
v[3740] = v[613];
v[3741] = v[611];
v[3742] = 0e0;
v[3743] = 0e0;
v[3744] = 0e0;
v[3745] = 0e0;
v[3746] = 0e0;
v[3747] = 0e0;
v[3748] = 0e0;
v[3749] = 0e0;
v[3750] = 0e0;
v[3697] = 0e0;
v[3698] = 0e0;
v[3699] = 0e0;
v[3700] = 0e0;
v[3701] = 0e0;
v[3702] = 0e0;
v[3703] = 0e0;
v[3704] = 0e0;
v[3705] = 0e0;
v[3706] = 0e0;
v[3707] = 0e0;
v[3708] = 0e0;
v[3709] = -v[614];
v[3710] = -v[613];
v[3711] = -v[611];
v[3712] = 0e0;
v[3713] = 0e0;
v[3714] = 0e0;
v[3679] = 0e0;
v[3680] = 0e0;
v[3681] = 0e0;
v[3682] = 0e0;
v[3683] = 0e0;
v[3684] = 0e0;
v[3685] = 0e0;
v[3686] = 0e0;
v[3687] = 0e0;
v[3688] = 0e0;
v[3689] = 0e0;
v[3690] = 0e0;
v[3691] = 0e0;
v[3692] = 0e0;
v[3693] = 0e0;
v[3694] = -v[614];
v[3695] = -v[613];
v[3696] = -v[611];
v[3769] = v[614];
v[3770] = v[613];
v[3771] = v[611];
v[3772] = 0e0;
v[3773] = 0e0;
v[3774] = 0e0;
v[3775] = 0e0;
v[3776] = 0e0;
v[3777] = 0e0;
v[3778] = 0e0;
v[3779] = 0e0;
v[3780] = 0e0;
v[3781] = 0e0;
v[3782] = 0e0;
v[3783] = 0e0;
v[3784] = 0e0;
v[3785] = 0e0;
v[3786] = 0e0;
v[3715] = 0e0;
v[3716] = 0e0;
v[3717] = 0e0;
v[3718] = 0e0;
v[3719] = 0e0;
v[3720] = 0e0;
v[3721] = 0e0;
v[3722] = 0e0;
v[3723] = 0e0;
v[3724] = -v[614];
v[3725] = -v[613];
v[3726] = -v[611];
v[3727] = 0e0;
v[3728] = 0e0;
v[3729] = 0e0;
v[3730] = 0e0;
v[3731] = 0e0;
v[3732] = 0e0;
v[2295] = v[614] / 2e0;
v[615] = -(v[177] * v[611]) - v[173] * v[613] - v[169] * v[614];
v[616] = v[149] * v[611] + v[145] * v[613] + v[141] * v[614];
v[617] = (-(v[179] * v[611]) - v[175] * v[613] - v[171] * v[614] - v[615]) / 2e0;
v[618] = (-(v[178] * v[611]) - v[174] * v[613] - v[170] * v[614] - v[615]) / 2e0;
v[619] = (v[151] * v[611] + v[147] * v[613] + v[143] * v[614] - v[616]) / 2e0;
v[620] = (v[150] * v[611] + v[146] * v[613] + v[142] * v[614] - v[616]) / 2e0;
for (i590 = 1; i590 <= 18; i590++) {
	v[630] = v[3584 + i590];
	v[629] = v[3566 + i590];
	v[627] = v[3548 + i590];
	v[626] = v[3530 + i590];
	v[628] = (-v[626] - v[627]) / 2e0;
	v[2288] = 2e0*v[628];
	v[631] = (-v[629] - v[630]) / 2e0;
	v[2289] = -2e0*v[631];
	v[632] = (v[141] * v[2288] + v[169] * v[2289] + 2e0*v[3624 + i590] + v[142] * v[626] + v[143] * v[627] - v[170] * v[629]
		- v[171] * v[630]) / 2e0;
	v[633] = (v[145] * v[2288] + v[173] * v[2289] + 2e0*v[3642 + i590] + v[146] * v[626] + v[147] * v[627] - v[174] * v[629]
		- v[175] * v[630]) / 2e0;
	v[634] = (v[149] * v[2288] + v[177] * v[2289] + 2e0*v[3660 + i590] + v[150] * v[626] + v[151] * v[627] - v[178] * v[629]
		- v[179] * v[630]) / 2e0;
	v[2290] = v[199] * v[632] + v[200] * v[633] + v[201] * v[634];
	v[635] = v[2290] / v[482];
	v[637] = -(v[2290] * v[610] * v[636]);
	v[656] = v[637];
	v[638] = v[212] * v[635];
	v[639] = v[632];
	v[654] = v[639];
	v[640] = v[633];
	v[652] = v[640];
	v[641] = v[634];
	v[650] = v[641];
	v[642] = 0e0;
	v[643] = 0e0;
	v[644] = 0e0;
	b645 = b483;
	if (b645) {
		v[2294] = v[228] * v[639];
		v[2293] = v[229] * v[640];
		v[2292] = v[230] * v[641];
		b646 = b485;
		if (b646) {
			v[644] = v[2282] * v[641];
			v[641] = 0e0;
			v[643] = v[2282] * v[640];
			v[640] = 0e0;
			v[642] = v[2282] * v[639];
			v[639] = 0e0;
			v[637] = v[637] - (-v[2292] - v[2293] - v[2294])*v[31] * v[34] * Power(v[494], v[495]);
		}
		else {
			v[644] = v[491] * v[650];
			v[641] = 0e0;
			v[643] = v[491] * v[652];
			v[640] = 0e0;
			v[642] = v[491] * v[654];
			v[639] = 0e0;
			v[637] = v[656] + v[2291] * (v[2292] + v[2293] + v[2294])*Power(v[482], v[500]);
		};
	}
	else {
	};
	b657 = b226;
	if (b657) {
		v[658] = -v[644];
		v[659] = -v[643];
		v[660] = -v[642];
	}
	else {
		v[658] = v[644];
		v[659] = v[643];
		v[660] = v[642];
	};
	v[637] = v[211] * v[608] * v[635] + v[637] + v[212] * (v[603] * v[632] + v[602] * v[633] + v[601] * v[634] + v[201] * v[658]
		+ v[200] * v[659] + v[199] * v[660]);
	v[668] = (v[610] * v[634] + v[201] * v[637]) / v[482] + v[601] * v[638] + v[213] * v[658];
	v[670] = (v[610] * v[633] + v[200] * v[637]) / v[482] + v[602] * v[638] + v[213] * v[659];
	v[672] = (v[610] * v[632] + v[199] * v[637]) / v[482] + v[603] * v[638] + v[213] * v[660];
	v[3787] = v[614] * v[628] + v[125] * v[672];
	v[3788] = v[613] * v[628] + v[125] * v[670];
	v[3789] = v[611] * v[628] + v[125] * v[668];
	v[3790] = v[2295] * v[626] + v[127] * v[672];
	v[3791] = v[2296] * v[626] + v[127] * v[670];
	v[3792] = v[2297] * v[626] + v[127] * v[668];
	v[3793] = v[2295] * v[627] + v[129] * v[672];
	v[3794] = v[2296] * v[627] + v[129] * v[670];
	v[3795] = v[2297] * v[627] + v[129] * v[668];
	v[3796] = -(v[614] * v[631]) - v[153] * v[672];
	v[3797] = -(v[613] * v[631]) - v[153] * v[670];
	v[3798] = -(v[611] * v[631]) - v[153] * v[668];
	v[3799] = -(v[2295] * v[629]) - v[155] * v[672];
	v[3800] = -(v[2296] * v[629]) - v[155] * v[670];
	v[3801] = -(v[2297] * v[629]) - v[155] * v[668];
	v[3802] = -(v[2295] * v[630]) - v[157] * v[672];
	v[3803] = -(v[2296] * v[630]) - v[157] * v[670];
	v[3804] = -(v[2297] * v[630]) - v[157] * v[668];
	v[673] = v[3714 + i590] - v[177] * v[668] - v[173] * v[670] - v[169] * v[672];
	v[674] = v[3768 + i590] + v[149] * v[668] + v[145] * v[670] + v[141] * v[672];
	v[675] = (v[3678 + i590] - v[179] * v[668] - v[175] * v[670] - v[171] * v[672] - v[673]) / 2e0;
	v[676] = (v[3696 + i590] - v[178] * v[668] - v[174] * v[670] - v[170] * v[672] - v[673]) / 2e0;
	v[677] = (v[3732 + i590] + v[151] * v[668] + v[147] * v[670] + v[143] * v[672] - v[674]) / 2e0;
	v[678] = (v[3750 + i590] + v[150] * v[668] + v[146] * v[670] + v[142] * v[672] - v[674]) / 2e0;
	Rc[i590 - 1] += v[3512 + i590] + v[620] * v[626] + v[619] * v[627] + v[618] * v[629] + v[617] * v[630];
	for (i623 = 1; i623 <= 18; i623++) {
		Kc[i590 - 1][i623 - 1] += v[3786 + i623] + v[3584 + i623] * v[675] + v[3566 + i623] * v[676] + v[3548 + i623] * v[677]
			+ v[3530 + i623] * v[678];
	};/* end for */
};/* end for */
v[760] = 0e0;
v[761] = 0e0;
v[762] = 0e0;
b763 = b483;
if (b763) {
	b764 = b29;
	if (b764) {
		b765 = b521;
		if (b765) {
			v[762] = 0e0;
			v[761] = 0e0;
			v[760] = 0e0;
		}
		else {
		};
	}
	else {
		b766 = b555;
		if (b766) {
			v[762] = 0e0;
			v[761] = 0e0;
			v[760] = 0e0;
		}
		else {
		};
	};
}
else {
};
v[2300] = v[760] * v[9];
v[2299] = v[761] * v[9];
v[2298] = v[762] * v[9];
v[2310] = (v[2298] * v[812]) / 2e0;
v[2311] = (v[2298] * v[814]) / 2e0;
v[2316] = (v[2298] * v[816]) / 2e0;
v[2317] = (v[2298] * v[818]) / 2e0;
v[2312] = (v[2299] * v[812]) / 2e0;
v[2313] = (v[2299] * v[814]) / 2e0;
v[2318] = (v[2299] * v[816]) / 2e0;
v[2319] = (v[2299] * v[818]) / 2e0;
v[2314] = (v[2300] * v[812]) / 2e0;
v[2315] = (v[2300] * v[814]) / 2e0;
v[2320] = (v[2300] * v[816]) / 2e0;
v[2321] = (v[2300] * v[818]) / 2e0;
v[843] = v[4] * v[762];
v[1336] = -(v[230] * v[843]);
v[844] = v[4] * v[761];
v[1337] = -(v[229] * v[844]);
v[1339] = v[1336] + v[1337];
v[845] = v[4] * v[760];
v[1332] = -(v[228] * v[845]);
v[1340] = v[1332] + v[1336];
v[1338] = v[1332] + v[1337];
v[846] = 0e0;
v[847] = 0e0;
v[848] = 0e0;
v[849] = 0e0;
v[850] = 0e0;
b851 = b483;
if (b851) {
	v[852] = 0e0;
	v[853] = 0e0;
	v[854] = 0e0;
	b855 = b502;
	if (b855) {
		v[854] = 0e0;
		v[853] = 0e0;
		v[852] = 0e0;
	}
	else {
	};
	v[859] = v[457] * v[852] + v[471] * v[853] + v[473] * v[854];
	v[2301] = v[7] * v[859];
	b856 = b485;
	if (b856) {
		v[848] = v[497] * v[854];
		v[847] = v[497] * v[853];
		v[846] = v[497] * v[852];
		v[850] = (v[2301] * v[495] * v[860] * Power(v[494], v[1442])) / v[861];
	}
	else {
		v[848] = v[501] * v[854];
		v[847] = v[501] * v[853];
		v[846] = v[501] * v[852];
		v[849] = (v[2301] * v[500] * v[868] * Power(v[482], v[1445])) / v[869];
	};
}
else {
};
v[2322] = v[456] * v[846];
v[2323] = v[464] * v[847];
v[2302] = v[230] * v[848];
v[910] = v[125] * v[2302];
v[906] = v[127] * v[2302];
v[902] = v[129] * v[2302];
v[898] = -(v[153] * v[2302]);
v[894] = -(v[155] * v[2302]);
v[890] = -(v[157] * v[2302]);
v[870] = 0e0;
v[871] = 0e0;
v[872] = 0e0;
b873 = b483;
if (b873) {
	b874 = b485;
	if (b874) {
		v[872] = 0e0;
		v[871] = 0e0;
		v[870] = 0e0;
		v[849] = v[849] - v[850];
	}
	else {
	};
}
else {
};
v[870] = v[1422] * v[848] + v[870];
v[871] = v[1415] * v[848] + v[871];
v[887] = v[1178] * v[848];
v[872] = v[1377] * v[848] + v[872];
v[924] = v[1358] * v[847];
v[871] = v[871] + v[847] * v[956];
v[949] = v[1343] * v[846];
v[950] = v[182] * v[846] + v[181] * v[847];
v[951] = v[185] * v[846] + v[184] * v[847];
v[952] = v[188] * v[846] + v[187] * v[847];
v[953] = v[191] * v[846] + v[190] * v[847];
v[954] = v[194] * v[846] + v[193] * v[847];
v[955] = v[197] * v[846] + v[196] * v[847];
v[1214] = v[125] * v[950] + v[127] * v[951] + v[129] * v[952] - v[153] * v[953] - v[155] * v[954] - v[157] * v[955];
v[870] = v[870] + v[846] * v[956];
v[957] = v[228] * v[846] + v[229] * v[847];
v[2324] = v[472] * v[848] + v[230] * v[957];
v[981] = 0e0;
v[982] = 0e0;
v[983] = 0e0;
v[984] = 0e0;
v[985] = 0e0;
v[986] = 0e0;
v[987] = 0e0;
v[988] = 0e0;
v[989] = 0e0;
b990 = b30;
if (b990) {
	v[887] = -(v[453] * v[843]) + v[887];
	v[924] = -(v[452] * v[844]) + v[924];
	v[992] = v[1339] * v[228] + v[845] * v[991];
	v[994] = v[1340] * v[229] + v[844] * v[993];
	v[996] = v[1338] * v[230] + v[843] * v[995];
	v[949] = -(v[451] * v[845]) + v[949];
	v[870] = v[1339] * v[451] - v[2303] * v[845] + v[870];
	v[871] = v[1340] * v[452] + v[2304] * v[844] + v[871];
	v[872] = v[1338] * v[453] + v[2305] * v[843] + v[872];
	v[981] = v[1] * v[992];
	v[982] = v[2] * v[992];
	v[983] = v[3] * v[992];
	v[1000] = -(v[162] * v[992]);
	v[1001] = -(v[160] * v[992]);
	v[1002] = -(v[158] * v[992]);
	v[1003] = v[134] * v[992];
	v[1004] = v[132] * v[992];
	v[1005] = v[130] * v[992];
	v[984] = v[1] * v[994];
	v[985] = v[2] * v[994];
	v[986] = v[3] * v[994];
	v[1006] = -(v[162] * v[994]);
	v[1007] = -(v[160] * v[994]);
	v[1008] = -(v[158] * v[994]);
	v[1009] = v[134] * v[994];
	v[1010] = v[132] * v[994];
	v[1011] = v[130] * v[994];
	v[987] = v[1] * v[996];
	v[988] = v[2] * v[996];
	v[989] = v[3] * v[996];
	v[1012] = -(v[162] * v[996]);
	v[1013] = -(v[160] * v[996]);
	v[1014] = -(v[158] * v[996]);
	v[1015] = v[134] * v[996];
	v[1016] = v[132] * v[996];
	v[1017] = v[130] * v[996];
}
else {
	v[1005] = 0e0;
	v[1004] = 0e0;
	v[1003] = 0e0;
	v[1011] = 0e0;
	v[1010] = 0e0;
	v[1009] = 0e0;
	v[1017] = 0e0;
	v[1016] = 0e0;
	v[1015] = 0e0;
	v[1002] = 0e0;
	v[1001] = 0e0;
	v[1000] = 0e0;
	v[1008] = 0e0;
	v[1007] = 0e0;
	v[1006] = 0e0;
	v[1014] = 0e0;
	v[1013] = 0e0;
	v[1012] = 0e0;
};
v[2357] = v[981] / 2e0;
v[2358] = v[985] / 2e0;
v[2359] = v[989] / 2e0;
b1018 = b416;
if (b1018) {
	v[1052] = -(v[433] * v[982]);
	v[1047] = v[433] * v[983];
	v[1034] = v[433] * v[986];
	v[1021] = -(v[2278] * v[989]);
	v[1025] = v[433] * v[988];
	v[1029] = v[433] * v[987];
	v[1033] = v[1025] + v[1034];
	v[1038] = -(v[2278] * v[985]);
	v[1042] = v[433] * v[984];
	v[1046] = v[1029] + v[1047];
	v[1053] = v[1042] - v[1052];
	v[1055] = v[1054] * v[2357] + v[1036] * v[2358] + v[1019] * v[2359] + v[1049] * v[982] + v[1044] * v[983]
		+ v[1040] * v[984] + v[1031] * v[986] + v[1027] * v[987] + v[1023] * v[988];
	v[1286] = v[1021] + v[1038] - (4e0*v[1055]) / (v[1059] * v[1059]);
	v[2356] = 4e0*v[1286];
	v[1284] = -v[1021] + v[1286] - v[2278] * v[981];
	v[2355] = 4e0*(v[1021] - v[1038] + v[1284]);
	v[1060] = v[1042] + v[1052] + v[1033] * v[2277] + v[1046] * v[2306] + 2e0*v[1284] * v[432];
	v[1062] = (-2e0*v[1029] + 2e0*v[1047] + v[1053] * v[430] + v[2355] * v[431] + v[1033] * v[432]) / 2e0;
	v[1063] = (2e0*v[1025] - 2e0*v[1034] + v[2356] * v[430] + v[1053] * v[431] + v[1046] * v[432]) / 2e0;
	v[2307] = v[1063] * v[418] + v[1062] * v[419] + v[1060] * v[420];
	v[1272] = v[2307] * v[429];
	v[1269] = v[2307] * v[423];
	v[1066] = v[1269] * v[428] + v[1272] / (Power(cos(v[1064]), 2)*sqrt(v[1273]));
	v[2308] = v[1066] / v[421];
	v[1067] = v[1060] * v[2276] + v[2308] * v[420];
	v[1069] = v[1062] * v[2276] + v[2308] * v[419];
	v[1070] = v[1063] * v[2276] + v[2308] * v[418];
	v[870] = -(v[1067] * v[224]) + v[1069] * v[225] + v[870];
	v[871] = v[1067] * v[223] - v[1070] * v[225] + v[871];
	v[872] = -(v[1069] * v[223]) + v[1070] * v[224] + v[872];
}
else {
};
v[870] = v[870] + v[2350] * v[949];
v[870] = v[1214] * v[229] + v[870];
v[871] = v[1214] * v[228] + v[871] + v[2347] * v[924];
v[872] = v[872] + v[2344] * v[887] + v[1178] * v[957];
b1071 = b226;
if (b1071) {
	v[1072] = -v[872];
	v[1073] = -v[871];
	v[1074] = -v[870];
}
else {
	v[1072] = v[872];
	v[1073] = v[871];
	v[1074] = v[870];
};
v[1079] = v[1074] * v[199] + v[1073] * v[200] + v[1072] * v[201];
v[849] = v[1079] * v[212] + v[849];
v[2309] = v[849] / v[482];
v[1081] = v[1072] * v[213] + v[201] * v[2309];
v[1082] = v[1073] * v[213] + v[200] * v[2309];
v[1083] = v[1074] * v[213] + v[199] * v[2309];
v[1012] = v[1012] - v[1081] * v[157] + v[2310];
v[1006] = v[1006] - v[1082] * v[157] + v[2312];
v[1000] = v[1000] - v[1083] * v[157] + v[2314];
v[1013] = v[1013] - v[1081] * v[155] + v[2311];
v[1007] = v[1007] - v[1082] * v[155] + v[2313];
v[1001] = v[1001] - v[1083] * v[155] + v[2315];
v[1014] = v[1014] - v[1081] * v[153] - v[2310] - v[2311];
v[1008] = v[1008] - v[1082] * v[153] - v[2312] - v[2313];
v[1002] = v[1002] - v[1083] * v[153] - v[2314] - v[2315];
v[1015] = v[1015] + v[1081] * v[129] + v[2316];
v[1009] = v[1009] + v[1082] * v[129] + v[2318];
v[1003] = v[1003] + v[1083] * v[129] + v[2320];
v[1016] = v[1016] + v[1081] * v[127] + v[2317];
v[1010] = v[1010] + v[1082] * v[127] + v[2319];
v[1004] = v[1004] + v[1083] * v[127] + v[2321];
v[1017] = v[1017] + v[1081] * v[125] - v[2316] - v[2317];
v[1011] = v[1011] + v[1082] * v[125] - v[2318] - v[2319];
v[1005] = v[1005] + v[1083] * v[125] - v[2320] - v[2321];
v[3809] = v[1005] - v[523] * v[685] - v[524] * v[686] - v[525] * v[687] + v[10] * (v[125] * v[2322] + v[458] * v[847] + (-
(v[685] * v[760]) - v[686] * v[761] - v[687] * v[762])*v[9] + v[228] * v[910]);
v[3810] = v[1011] - v[523] * v[689] - v[524] * v[690] - v[525] * v[691] + v[10] * (v[125] * v[2323] + v[458] * v[846] + (-
(v[689] * v[760]) - v[690] * v[761] - v[691] * v[762])*v[9] + v[229] * v[910]);
v[3811] = v[1017] - v[523] * v[693] - v[524] * v[694] - v[525] * v[695] + v[10] * (v[125] * v[2324] + (-(v[693] * v[760])
	- v[694] * v[761] - v[695] * v[762])*v[9]);
v[3812] = v[1004] - v[523] * v[697] - v[524] * v[698] - v[525] * v[699] + v[10] * (v[127] * v[2322] + v[459] * v[847] + (-
(v[697] * v[760]) - v[698] * v[761] - v[699] * v[762])*v[9] + v[228] * v[906]);
v[3813] = v[1010] - v[523] * v[701] - v[524] * v[702] - v[525] * v[703] + v[10] * (v[127] * v[2323] + v[459] * v[846] + (-
(v[701] * v[760]) - v[702] * v[761] - v[703] * v[762])*v[9] + v[229] * v[906]);
v[3814] = v[1016] - v[523] * v[705] - v[524] * v[706] - v[525] * v[707] + v[10] * (v[127] * v[2324] + (-(v[705] * v[760])
	- v[706] * v[761] - v[707] * v[762])*v[9]);
v[3815] = v[1003] - v[523] * v[709] - v[524] * v[710] - v[525] * v[711] + v[10] * (v[129] * v[2322] + v[460] * v[847] + (-
(v[709] * v[760]) - v[710] * v[761] - v[711] * v[762])*v[9] + v[228] * v[902]);
v[3816] = v[1009] - v[523] * v[713] - v[524] * v[714] - v[525] * v[715] + v[10] * (v[129] * v[2323] + v[460] * v[846] + (-
(v[713] * v[760]) - v[714] * v[761] - v[715] * v[762])*v[9] + v[229] * v[902]);
v[3817] = v[1015] - v[523] * v[717] - v[524] * v[718] - v[525] * v[719] + v[10] * (v[129] * v[2324] + (-(v[717] * v[760])
	- v[718] * v[761] - v[719] * v[762])*v[9]);
v[3818] = v[1002] - v[523] * v[721] - v[524] * v[722] - v[525] * v[723] + v[10] * (-(v[153] * v[2322]) + v[461] * v[847]
	+ v[228] * v[898] + (-(v[721] * v[760]) - v[722] * v[761] - v[723] * v[762])*v[9]);
v[3819] = v[1008] - v[523] * v[725] - v[524] * v[726] - v[525] * v[727] + v[10] * (-(v[153] * v[2323]) + v[461] * v[846]
	+ v[229] * v[898] + (-(v[725] * v[760]) - v[726] * v[761] - v[727] * v[762])*v[9]);
v[3820] = v[1014] - v[523] * v[729] - v[524] * v[730] - v[525] * v[731] + v[10] * (-(v[153] * v[2324]) + (-(v[729] * v[760]
	) - v[730] * v[761] - v[731] * v[762])*v[9]);
v[3821] = v[1001] - v[523] * v[733] - v[524] * v[734] - v[525] * v[735] + v[10] * (-(v[155] * v[2322]) + v[462] * v[847]
	+ v[228] * v[894] + (-(v[733] * v[760]) - v[734] * v[761] - v[735] * v[762])*v[9]);
v[3822] = v[1007] - v[523] * v[737] - v[524] * v[738] - v[525] * v[739] + v[10] * (-(v[155] * v[2323]) + v[462] * v[846]
	+ v[229] * v[894] + (-(v[737] * v[760]) - v[738] * v[761] - v[739] * v[762])*v[9]);
v[3823] = v[1013] - v[523] * v[741] - v[524] * v[742] - v[525] * v[743] + v[10] * (-(v[155] * v[2324]) + (-(v[741] * v[760]
	) - v[742] * v[761] - v[743] * v[762])*v[9]);
v[3824] = v[1000] - v[523] * v[745] - v[524] * v[746] - v[525] * v[747] + v[10] * (-(v[157] * v[2322]) + v[463] * v[847]
	+ v[228] * v[890] + (-(v[745] * v[760]) - v[746] * v[761] - v[747] * v[762])*v[9]);
v[3825] = v[1006] - v[523] * v[749] - v[524] * v[750] - v[525] * v[751] + v[10] * (-(v[157] * v[2323]) + v[463] * v[846]
	+ v[229] * v[890] + (-(v[749] * v[760]) - v[750] * v[761] - v[751] * v[762])*v[9]);
v[3826] = v[1012] - v[523] * v[753] - v[524] * v[754] - v[525] * v[755] + v[10] * (-(v[157] * v[2324]) + (-(v[753] * v[760]
	) - v[754] * v[761] - v[755] * v[762])*v[9]);
for (i758 = 1; i758 <= 18; i758++) {
	i2342 = (i758 == 18 ? 1 : 0);
	i2341 = (i758 == 15 ? 1 : 0);
	i2340 = (i758 == 9 ? 1 : 0);
	i2339 = (i758 == 6 ? 1 : 0);
	i2338 = (i758 == 17 ? 1 : 0);
	i2337 = (i758 == 14 ? 1 : 0);
	i2336 = (i758 == 8 ? 1 : 0);
	i2335 = (i758 == 5 ? 1 : 0);
	i2334 = (i758 == 16 ? 1 : 0);
	i2333 = (i758 == 13 ? 1 : 0);
	i2332 = (i758 == 7 ? 1 : 0);
	i2331 = (i758 == 4 ? 1 : 0);
	i2330 = (i758 == 12 ? 1 : 0);
	i2329 = (i758 == 11 ? 1 : 0);
	i2328 = (i758 == 10 ? 1 : 0);
	i2327 = (i758 == 3 ? 1 : 0);
	i2326 = (i758 == 2 ? 1 : 0);
	i2325 = (i758 == 1 ? 1 : 0);
	v[1453] = -i2325 / 2e0;
	v[1454] = -i2326 / 2e0;
	v[1455] = -i2327 / 2e0;
	v[1458] = -i2328 / 2e0;
	v[1459] = -i2329 / 2e0;
	v[1460] = -i2330 / 2e0;
	v[1127] = i2325 * v[10];
	v[1129] = i2326 * v[10];
	v[1131] = i2327 * v[10];
	v[1133] = i2331 * v[10];
	v[1135] = i2335 * v[10];
	v[1137] = i2339 * v[10];
	v[1139] = i2332 * v[10];
	v[1141] = i2336 * v[10];
	v[1143] = i2340 * v[10];
	v[1145] = i2328 * v[10];
	v[1147] = i2329 * v[10];
	v[1149] = i2330 * v[10];
	v[1151] = i2333 * v[10];
	v[1153] = i2337 * v[10];
	v[1155] = i2341 * v[10];
	v[1157] = i2334 * v[10];
	v[2377] = v[1127] * v[125] + v[1133] * v[127] + v[1139] * v[129] - v[1145] * v[153] - v[1151] * v[155] - v[1157] * v[157];
	v[1158] = i2325 * v[125] + i2331 * v[127] + i2332 * v[129] - i2328 * v[153] - i2333 * v[155] - i2334 * v[157];
	v[1160] = i2338 * v[10];
	v[2376] = v[1129] * v[125] + v[1135] * v[127] + v[1141] * v[129] - v[1147] * v[153] - v[1153] * v[155] - v[1160] * v[157];
	v[1161] = i2326 * v[125] + i2335 * v[127] + i2336 * v[129] - i2329 * v[153] - i2337 * v[155] - i2338 * v[157];
	v[1163] = i2342 * v[10];
	v[2345] = v[1131] * v[125] + v[1137] * v[127] + v[1143] * v[129] - v[1149] * v[153] - v[1155] * v[155] - v[1163] * v[157];
	v[1164] = i2327 * v[125] + i2339 * v[127] + i2340 * v[129] - i2330 * v[153] - i2341 * v[155] - i2342 * v[157];
	v[2343] = v[1158] * v[199] + v[1161] * v[200] + v[1164] * v[201];
	v[1168] = v[2343] / v[482];
	v[1166] = -(v[2343] * v[636] * v[849]);
	v[1167] = v[1168] * v[212];
	v[1169] = v[1168];
	v[1170] = v[1167] * v[199] + v[1158] * v[213];
	v[1171] = v[1167] * v[200] + v[1161] * v[213];
	v[1172] = v[1167] * v[201] + v[1164] * v[213];
	b1173 = b226;
	if (b1173) {
		v[1174] = -v[1170];
		v[1175] = -v[1171];
		v[1176] = -v[1172];
	}
	else {
		v[1174] = v[1170];
		v[1175] = v[1171];
		v[1176] = v[1172];
	};
	v[2349] = -(v[1174] * v[229]);
	v[2348] = -(v[1175] * v[228]);
	v[2346] = -(v[1176] * v[957]);
	v[1177] = v[1176] * v[2344];
	v[1179] = v[1176] * v[1178] + v[230] * v[2345];
	v[1180] = 2e0*v[1176] * v[887] + v[2345] * v[957];
	v[1187] = -(i2342*v[1081]) - i2338 * v[1082] - i2334 * v[1083] + (-(v[1176] * v[198]) - v[1163] * v[230])*v[957];
	v[1188] = -(i2341*v[1081]) - i2337 * v[1082] - i2333 * v[1083] + (-(v[1176] * v[195]) - v[1155] * v[230])*v[957];
	v[1189] = -(i2330*v[1081]) - i2329 * v[1082] - i2328 * v[1083] + (-(v[1176] * v[192]) - v[1149] * v[230])*v[957];
	v[1190] = i2340 * v[1081] + i2336 * v[1082] + i2332 * v[1083] + (v[1176] * v[189] + v[1143] * v[230])*v[957];
	v[1191] = i2339 * v[1081] + i2335 * v[1082] + i2331 * v[1083] + (v[1176] * v[186] + v[1137] * v[230])*v[957];
	v[1192] = i2327 * v[1081] + i2326 * v[1082] + i2325 * v[1083] + (v[1176] * v[183] + v[1131] * v[230])*v[957];
	v[1193] = v[1175] * v[2347];
	v[1200] = v[1175] * v[1214];
	v[1201] = 2e0*v[1175] * v[924];
	v[1208] = v[157] * (v[2348] + v[2349]);
	v[1209] = v[155] * (v[2348] + v[2349]);
	v[1210] = v[153] * (v[2348] + v[2349]);
	v[1211] = v[129] * (-v[2348] - v[2349]);
	v[1212] = v[127] * (-v[2348] - v[2349]);
	v[1213] = v[125] * (-v[2348] - v[2349]);
	v[1201] = v[1201] + v[1174] * v[1214];
	v[1221] = v[1174] * v[2350];
	v[1200] = v[1200] + 2e0*v[1174] * v[949];
	v[1222] = 0e0;
	v[1223] = 0e0;
	v[1224] = 0e0;
	v[1225] = 0e0;
	v[1226] = 0e0;
	v[1227] = 0e0;
	v[1228] = 0e0;
	v[1229] = 0e0;
	v[1230] = 0e0;
	v[1231] = 0e0;
	v[1232] = 0e0;
	v[1233] = 0e0;
	v[1234] = 0e0;
	v[1235] = 0e0;
	v[1236] = 0e0;
	v[1237] = 0e0;
	v[1238] = 0e0;
	v[1239] = 0e0;
	v[1240] = 0e0;
	v[1241] = 0e0;
	v[1242] = 0e0;
	v[1243] = 0e0;
	v[1244] = 0e0;
	v[1245] = 0e0;
	v[1246] = 0e0;
	v[1247] = 0e0;
	v[1248] = 0e0;
	v[1249] = 0e0;
	v[1250] = 0e0;
	v[1251] = 0e0;
	v[1252] = 0e0;
	v[1253] = 0e0;
	b1254 = b416;
	if (b1254) {
		v[1255] = v[1176] * v[224];
		v[1256] = -(v[1176] * v[223]);
		v[1257] = v[1255] - v[1175] * v[225];
		v[1258] = v[1175] * v[223];
		v[1259] = v[1256] + v[1174] * v[225];
		v[1260] = v[1258] - v[1174] * v[224];
		v[2354] = v[1063] * v[1257] + v[1062] * v[1259] + v[1060] * v[1260];
		v[2352] = v[1257] * v[418] + v[1259] * v[419] + v[1260] * v[420];
		v[1261] = v[2352] / v[421];
		v[2353] = v[1261] * v[2487];
		v[1271] = v[1261] * v[2351];
		v[1225] = -(v[1066] * v[1262] * v[2352]);
		v[1263] = v[1257] * v[2276] + v[2353] * v[418];
		v[1278] = v[1263] * v[2388];
		v[1265] = v[1259] * v[2276] + v[2353] * v[419];
		v[1282] = 2e0*v[1265] * v[431];
		v[1266] = v[1260] * v[2276] + v[2353] * v[420];
		v[1279] = v[1266] * v[2385];
		v[1228] = v[1271] * v[2307] + v[2354] * v[423];
		v[1227] = v[1261] * v[1269] * v[427];
		v[1226] = v[1261] * v[2307] * v[428] + v[2354] * v[429];
		v[1252] = v[1271] * v[1272] * v[423];
		v[1253] = v[1065] * v[1261] * v[1267] * v[1272] * v[2488];
		v[1224] = v[1260] * v[2308] + v[1060] * v[2353];
		v[1223] = v[1259] * v[2308] + v[1062] * v[2353];
		v[1222] = v[1257] * v[2308] + v[1063] * v[2353];
		v[1274] = (v[1265] * v[430] + v[1263] * v[431]) / 2e0;
		v[1275] = v[1278] + v[1282];
		v[1276] = v[1275] + v[1279];
		v[1277] = (v[1266] * v[430] + v[1263] * v[432]) / 2e0;
		v[1280] = v[1278] + v[1279];
		v[1281] = (v[1266] * v[431] + v[1265] * v[432]) / 2e0;
		v[1283] = v[1279] + v[1282];
		v[1231] = (v[1046] * v[1263] + v[1033] * v[1265] + 4e0*v[1266] * v[1284]) / 2e0;
		v[1230] = (v[1053] * v[1263] + v[1033] * v[1266] + v[1265] * v[2355]) / 2e0;
		v[1229] = (v[1053] * v[1265] + v[1046] * v[1266] + v[1263] * v[2356]) / 2e0;
		v[1287] = v[1276] * v[2384];
		v[1251] = 8e0*v[1055] * v[1276] * v[2489];
		v[1250] = v[1287] * v[2357];
		v[1288] = v[1266] + v[1274];
		v[1289] = v[1266] - v[1274];
		v[1249] = v[1287] * v[982];
		v[1290] = -v[1265] + v[1277];
		v[1291] = v[1265] + v[1277];
		v[1248] = v[1287] * v[983];
		v[1236] = v[1040] * v[1287] + v[1288] * v[433];
		v[1247] = v[1287] * v[984];
		v[1237] = (v[1036] * v[1287] - v[1280] * v[433]) / 2e0;
		v[1246] = v[1287] * v[2358];
		v[1292] = v[1263] + v[1281];
		v[1293] = -v[1263] + v[1281];
		v[1245] = v[1287] * v[986];
		v[1239] = v[1027] * v[1287] + v[1290] * v[433];
		v[1244] = v[1287] * v[987];
		v[1240] = v[1023] * v[1287] + v[1292] * v[433];
		v[1243] = v[1287] * v[988];
		v[1241] = (v[1019] * v[1287] - v[1275] * v[433]) / 2e0;
		v[1242] = v[1287] * v[2359];
		v[1238] = v[1031] * v[1287] + v[1293] * v[433];
		v[1235] = v[1044] * v[1287] + v[1291] * v[433];
		v[1234] = v[1049] * v[1287] - v[1289] * v[433];
		v[1233] = (v[1054] * v[1287] - v[1283] * v[433]) / 2e0;
		v[1232] = -(v[1283] * v[2357]) - v[1280] * v[2358] - v[1275] * v[2359] - v[1289] * v[982] + v[1291] * v[983]
			+ v[1288] * v[984] + v[1293] * v[986] + v[1290] * v[987] + v[1292] * v[988];
	}
	else {
	};
	v[1294] = 0e0;
	v[1295] = 0e0;
	v[1296] = 0e0;
	v[1297] = 0e0;
	v[1298] = 0e0;
	v[1299] = 0e0;
	b1300 = b30;
	if (b1300) {
		v[2362] = -(v[1174] * v[451]);
		v[2361] = -(v[1175] * v[452]);
		v[2360] = -(v[1176] * v[453]);
		v[1334] = v[1174] * v[845];
		v[1307] = i2327 * v[130] + i2339 * v[132] + i2340 * v[134] - i2330 * v[158] - i2341 * v[160] - i2342 * v[162]
			+ v[1241] * v[3];
		v[1241] = 0e0;
		v[1308] = v[1307] + v[1240] * v[2];
		v[1240] = 0e0;
		v[1309] = v[1] * v[1239] + v[1308];
		v[2364] = -(v[1309] * v[230]);
		v[1239] = 0e0;
		v[1316] = i2326 * v[130] + i2335 * v[132] + i2336 * v[134] - i2329 * v[158] - i2337 * v[160] - i2338 * v[162]
			+ v[1238] * v[3];
		v[1238] = 0e0;
		v[1317] = v[1316] + v[1237] * v[2];
		v[1237] = 0e0;
		v[1318] = v[1] * v[1236] + v[1317];
		v[2363] = -(v[1318] * v[229]);
		v[1236] = 0e0;
		v[1325] = i2325 * v[130] + i2331 * v[132] + i2332 * v[134] - i2328 * v[158] - i2333 * v[160] - i2334 * v[162]
			+ v[1235] * v[3];
		v[1235] = 0e0;
		v[1326] = v[1325] + v[1234] * v[2];
		v[1234] = 0e0;
		v[1327] = v[1] * v[1233] + v[1326];
		v[2365] = -(v[1327] * v[228]);
		v[1233] = 0e0;
		v[1328] = v[1176] * v[843];
		v[1296] = v[1176] * v[1338];
		v[1200] = v[1200] + v[2360] * v[845];
		v[1201] = v[1201] + v[2360] * v[844];
		v[1330] = v[1175] * v[844];
		v[1331] = v[1328] + v[1330];
		v[1295] = v[1175] * v[1340];
		v[1200] = v[1200] + v[2361] * v[845];
		v[1180] = v[1180] + v[2361] * v[843];
		v[1333] = v[1330] + v[1334];
		v[1335] = v[1328] + v[1334];
		v[1294] = v[1174] * v[1339];
		v[1201] = v[1201] + v[2362] * v[844];
		v[1180] = v[1180] + v[2362] * v[843];
		v[1294] = v[1294] - v[1221] * v[845];
		v[1299] = v[1309] * v[843];
		v[1298] = v[1318] * v[844];
		v[1297] = v[1327] * v[845];
		v[1295] = v[1295] - v[1193] * v[844];
		v[1296] = v[1296] - v[1177] * v[843];
		v[1296] = v[1296] - v[1333] * v[230];
		v[1180] = v[1180] + v[1309] * v[1338] - v[1333] * v[453] + (v[2363] + v[2365])*v[843];
		v[1294] = v[1294] - v[1331] * v[228];
		v[1200] = v[1200] + v[1327] * v[1339] - v[1331] * v[451] + (v[2363] + v[2364])*v[845];
		v[1295] = v[1295] - v[1335] * v[229];
		v[1201] = v[1201] + v[1318] * v[1340] - v[1335] * v[452] + (v[2364] + v[2365])*v[844];
	}
	else {
	};
	v[2366] = v[1176] * v[848];
	v[1395] = v[125] * v[2366];
	v[1392] = v[127] * v[2366];
	v[1389] = v[129] * v[2366];
	v[1386] = -(v[153] * v[2366]);
	v[1383] = -(v[155] * v[2366]);
	v[1380] = -(v[157] * v[2366]);
	v[1341] = v[1179] * v[228] + v[2377] * v[456] + v[1129] * v[458] + v[1135] * v[459] + v[1141] * v[460] + v[1147] * v[461]
		+ v[1153] * v[462] + v[1160] * v[463] + v[1174] * v[956];
	v[1342] = v[1174] * v[846];
	v[1341] = v[1341] + v[1221] * v[1343] + v[1213] * v[182] + v[1212] * v[185] + v[1211] * v[188] + v[1210] * v[191]
		+ v[1209] * v[194] + v[1208] * v[197];
	v[1447] = v[1341];
	v[1356] = v[1179] * v[229] + v[1127] * v[458] + v[1133] * v[459] + v[1139] * v[460] + v[1145] * v[461] + v[1151] * v[462]
		+ v[1157] * v[463] + v[2376] * v[464] + v[1175] * v[956];
	v[1357] = v[1342] + v[1175] * v[847];
	v[1356] = v[1356] + v[1193] * v[1358];
	v[1356] = v[1356] + v[1213] * v[181] + v[1212] * v[184] + v[1211] * v[187] + v[1210] * v[190] + v[1209] * v[193]
		+ v[1208] * v[196];
	v[1448] = v[1356];
	v[1371] = v[1127] * v[228] + v[1129] * v[229];
	v[1372] = v[1133] * v[228] + v[1135] * v[229];
	v[1373] = v[1139] * v[228] + v[1141] * v[229];
	v[1374] = v[1145] * v[228] + v[1147] * v[229];
	v[1375] = v[1151] * v[228] + v[1153] * v[229];
	v[1376] = v[1157] * v[228] + v[1160] * v[229];
	v[2367] = v[125] * v[1371] + v[127] * v[1372] + v[129] * v[1373] - v[1374] * v[153] - v[1375] * v[155] - v[1376] * v[157];
	v[1378] = v[1176] * v[1377] + v[2345] * v[472];
	v[1200] = v[1200] + v[1343] * v[2366] + v[1179] * v[846] + v[1157] * v[890] + v[1151] * v[894] + v[1145] * v[898]
		+ v[1139] * v[902] + v[1133] * v[906] + v[1127] * v[910];
	v[1649] = v[1200];
	v[1201] = v[1201] + v[1358] * v[2366] + v[1179] * v[847] + v[1160] * v[890] + v[1153] * v[894] + v[1147] * v[898]
		+ v[1141] * v[902] + v[1135] * v[906] + v[1129] * v[910];
	v[1647] = v[1201];
	v[1397] = v[1187] + (-(v[1221] * v[196]) - v[1157] * v[456])*v[846] + (-(v[1193] * v[197]) - v[1160] * v[464]
		)*v[847] + (v[1176] * v[2496] - v[1163] * v[472])*v[848] + (v[2348] + v[2349])*v[955];
	v[1398] = v[1188] + (-(v[1221] * v[193]) - v[1151] * v[456])*v[846] + (-(v[1193] * v[194]) - v[1153] * v[464]
		)*v[847] + (v[1176] * v[2497] - v[1155] * v[472])*v[848] + (v[2348] + v[2349])*v[954];
	v[1399] = v[1189] + (-(v[1221] * v[190]) - v[1145] * v[456])*v[846] + (-(v[1193] * v[191]) - v[1147] * v[464]
		)*v[847] + (v[1176] * v[2498] - v[1149] * v[472])*v[848] + (v[2348] + v[2349])*v[953];
	v[1400] = v[1190] + (v[1221] * v[187] + v[1139] * v[456])*v[846] + (v[1193] * v[188] + v[1141] * v[464])*v[847] +
		(v[1176] * v[2499] + v[1143] * v[472])*v[848] + (-v[2348] - v[2349])*v[952];
	v[1401] = v[1191] + (v[1221] * v[184] + v[1133] * v[456])*v[846] + (v[1193] * v[185] + v[1135] * v[464])*v[847] +
		(v[1176] * v[2500] + v[1137] * v[472])*v[848] + (-v[2348] - v[2349])*v[951];
	v[1402] = v[1192] + (v[1221] * v[181] + v[1127] * v[456])*v[846] + (v[1193] * v[182] + v[1129] * v[464])*v[847] +
		(v[1176] * v[2501] + v[1131] * v[472])*v[848] + (-v[2348] - v[2349])*v[950];
	v[1378] = v[1177] * v[1178] + v[1378];
	v[1378] = v[1378] + v[1175] * v[1415];
	v[1416] = v[1175] * v[848];
	v[1378] = v[1378] + v[1174] * v[1422];
	v[1423] = v[1174] * v[848];
	v[1429] = 0e0;
	b1430 = b483;
	if (b1430) {
		b1431 = b485;
		if (b1431) {
			v[1429] = -v[1169];
			v[1174] = 0e0;
			v[1175] = 0e0;
			v[1176] = 0e0;
		}
		else {
		};
	}
	else {
	};
	v[1378] = v[1378] + v[230] * v[2367];
	v[1449] = v[1378];
	v[1180] = v[1180] + v[2367] * v[848];
	v[1645] = v[1180];
	v[1432] = 0e0;
	v[1433] = 0e0;
	v[1434] = 0e0;
	v[1435] = 0e0;
	v[1436] = 0e0;
	v[1437] = 0e0;
	v[1438] = 0e0;
	v[1439] = 0e0;
	b1440 = b483;
	if (b1440) {
		v[2368] = v[1447] * v[852];
		v[1450] = v[2368] + v[1448] * v[853] + v[1449] * v[854];
		b1441 = b485;
		if (b1441) {
			v[1630] = Power(v[494], v[1442]);
			v[1629] = (v[495] * v[860]) / v[861];
			v[1444] = v[1429] * v[1629] * v[7];
			v[1443] = v[1444] * v[1630];
			v[1438] = -((v[1443] * v[859]) / v[861]);
			v[1435] = v[1442] * v[1444] * v[859] * Power(v[494], -2e0 + v[495]);
			v[1429] = 0e0;
			v[1341] = 0e0;
			v[1356] = 0e0;
			v[1436] = v[1450];
			v[1378] = 0e0;
		}
		else {
			v[1638] = Power(v[482], v[1445]);
			v[1637] = (v[500] * v[868]) / v[869];
			v[1446] = v[1169] * v[1637] * v[7];
			v[1443] = v[1446] * v[1638];
			v[1439] = -((v[1443] * v[859]) / v[869]);
			v[1166] = v[1166] + v[1445] * v[1446] * v[859] * Power(v[482], -2e0 + v[500]);
			v[1169] = 0e0;
			v[1341] = 0e0;
			v[1356] = 0e0;
			v[1437] = v[1450];
			v[1378] = 0e0;
		};
		v[1432] = v[1443] * v[852];
		v[1433] = v[1443] * v[853];
		v[1434] = v[1443] * v[854];
	}
	else {
	};
	v[1636] = v[1432];
	v[1634] = v[1433];
	v[1632] = v[1434];
	v[1452] = ((i2331 / 2e0 + v[1453])*v[760] + (i2335 / 2e0 + v[1454])*v[761] + (i2339 / 2e0 + v[1455])*v[762])*v[9];
	v[1456] = ((i2332 / 2e0 + v[1453])*v[760] + (i2336 / 2e0 + v[1454])*v[761] + (i2340 / 2e0 + v[1455])*v[762])*v[9];
	v[1457] = ((i2333 / 2e0 + v[1458])*v[760] + (i2337 / 2e0 + v[1459])*v[761] + (i2341 / 2e0 + v[1460])*v[762])*v[9];
	v[1461] = ((i2334 / 2e0 + v[1458])*v[760] + (i2338 / 2e0 + v[1459])*v[761] + (i2342 / 2e0 + v[1460])*v[762])*v[9];
	v[1462] = -(i2325*v[685]) - i2326 * v[689] - i2327 * v[693] - i2331 * v[697] - i2335 * v[701] - i2339 * v[705]
		- i2332 * v[709] - i2336 * v[713] - i2340 * v[717] - i2328 * v[721] - i2329 * v[725] - i2330 * v[729] - i2333 * v[733]
		- i2337 * v[737] - i2341 * v[741] - i2334 * v[745] - i2338 * v[749] - i2342 * v[753];
	v[1478] = v[1462];
	v[1463] = -(i2325*v[686]) - i2326 * v[690] - i2327 * v[694] - i2331 * v[698] - i2335 * v[702] - i2339 * v[706]
		- i2332 * v[710] - i2336 * v[714] - i2340 * v[718] - i2328 * v[722] - i2329 * v[726] - i2330 * v[730] - i2333 * v[734]
		- i2337 * v[738] - i2341 * v[742] - i2334 * v[746] - i2338 * v[750] - i2342 * v[754];
	v[1476] = v[1463];
	v[1464] = -(i2325*v[687]) - i2326 * v[691] - i2327 * v[695] - i2331 * v[699] - i2335 * v[703] - i2339 * v[707]
		- i2332 * v[711] - i2336 * v[715] - i2340 * v[719] - i2328 * v[723] - i2329 * v[727] - i2330 * v[731] - i2333 * v[735]
		- i2337 * v[739] - i2341 * v[743] - i2334 * v[747] - i2338 * v[751] - i2342 * v[755];
	v[1474] = v[1464];
	v[1465] = 0e0;
	v[1466] = 0e0;
	v[1467] = 0e0;
	v[1468] = 0e0;
	v[1469] = 0e0;
	v[1470] = 0e0;
	b1471 = b483;
	if (b1471) {
		b1472 = b29;
		if (b1472) {
			b1473 = b521;
			if (b1473) {
				v[1470] = v[1464];
				v[1464] = 0e0;
				v[1469] = v[1463];
				v[1463] = 0e0;
				v[1468] = v[1462];
				v[1462] = 0e0;
			}
			else {
				v[1486] = v[1478] * v[2285];
				v[1484] = v[1476] * v[2285];
				v[1482] = v[1474] * v[2285];
				v[1464] = 0e0;
				v[1463] = 0e0;
				v[2370] = ((v[1478] * v[533] + v[1476] * v[534] + v[1474] * v[535])*v[6]) / sqrt(v[1487]);
				v[1462] = 0e0;
				v[2369] = ((v[1486] * v[516] + v[1484] * v[517] + v[1482] * v[518])*v[531]) / sqrt(v[1483]);
				v[1470] = v[2369] * v[518] + v[1482] * v[532];
				v[1469] = v[2369] * v[517] + v[1484] * v[532];
				v[1468] = v[2369] * v[516] + v[1486] * v[532];
				v[1467] = v[2370] * v[509];
				v[1466] = v[2370] * v[508];
				v[1465] = v[2370] * v[507];
			};
		}
		else {
			b1489 = b555;
			if (b1489) {
				v[1470] = v[1474];
				v[1464] = 0e0;
				v[1469] = v[1476];
				v[1463] = 0e0;
				v[1468] = v[1478];
				v[1462] = 0e0;
			}
			else {
				v[1495] = v[1478] * v[566] * v[6];
				v[1493] = v[1476] * v[566] * v[6];
				v[1492] = v[1474] * v[566] * v[6];
				v[1464] = 0e0;
				v[1463] = 0e0;
				v[2372] = ((v[1478] * v[563] + v[1476] * v[564] + v[1474] * v[565])*v[6]) / sqrt(v[1487]);
				v[1462] = 0e0;
				v[2371] = ((v[1495] * v[516] + v[1493] * v[517] + v[1492] * v[518])*v[561]) / sqrt(v[1483]);
				v[1470] = v[2371] * v[518] + v[1492] * v[562];
				v[1469] = v[2371] * v[517] + v[1493] * v[562];
				v[1468] = v[2371] * v[516] + v[1495] * v[562];
				v[1467] = v[2372] * v[509];
				v[1466] = v[2372] * v[508];
				v[1465] = v[2372] * v[507];
			};
		};
	}
	else {
	};
	v[1501] = -(i2342*v[525]) - (v[1470] * v[198] + v[1163] * v[762])*v[9];
	v[1502] = -(i2338*v[525]) - (v[1470] * v[197] + v[1160] * v[762])*v[9];
	v[1503] = -(i2334*v[525]) - (v[1470] * v[196] + v[1157] * v[762])*v[9];
	v[1504] = -(i2341*v[525]) - (v[1470] * v[195] + v[1155] * v[762])*v[9];
	v[1505] = -(i2337*v[525]) - (v[1470] * v[194] + v[1153] * v[762])*v[9];
	v[1506] = -(i2333*v[525]) - (v[1470] * v[193] + v[1151] * v[762])*v[9];
	v[1507] = -(i2330*v[525]) - (v[1470] * v[192] + v[1149] * v[762])*v[9];
	v[1508] = -(i2329*v[525]) - (v[1470] * v[191] + v[1147] * v[762])*v[9];
	v[1509] = -(i2328*v[525]) - (v[1470] * v[190] + v[1145] * v[762])*v[9];
	v[1510] = -(i2340*v[525]) - (v[1470] * v[189] + v[1143] * v[762])*v[9];
	v[1511] = -(i2336*v[525]) - (v[1470] * v[188] + v[1141] * v[762])*v[9];
	v[1512] = -(i2332*v[525]) - (v[1470] * v[187] + v[1139] * v[762])*v[9];
	v[1513] = -(i2339*v[525]) - (v[1470] * v[186] + v[1137] * v[762])*v[9];
	v[1514] = -(i2335*v[525]) - (v[1470] * v[185] + v[1135] * v[762])*v[9];
	v[1515] = -(i2331*v[525]) - (v[1470] * v[184] + v[1133] * v[762])*v[9];
	v[1516] = -(i2327*v[525]) - (v[1470] * v[183] + v[1131] * v[762])*v[9];
	v[1517] = -(i2326*v[525]) - (v[1470] * v[182] + v[1129] * v[762])*v[9];
	v[1518] = -(i2325*v[525]) - (v[1470] * v[181] + v[1127] * v[762])*v[9];
	v[1538] = -(i2342*v[524]) - (v[1469] * v[198] + v[1163] * v[761])*v[9];
	v[1539] = -(i2338*v[524]) - (v[1469] * v[197] + v[1160] * v[761])*v[9];
	v[1540] = -(i2334*v[524]) - (v[1469] * v[196] + v[1157] * v[761])*v[9];
	v[1541] = -(i2341*v[524]) - (v[1469] * v[195] + v[1155] * v[761])*v[9];
	v[1542] = -(i2337*v[524]) - (v[1469] * v[194] + v[1153] * v[761])*v[9];
	v[1543] = -(i2333*v[524]) - (v[1469] * v[193] + v[1151] * v[761])*v[9];
	v[1544] = -(i2330*v[524]) - (v[1469] * v[192] + v[1149] * v[761])*v[9];
	v[1545] = -(i2329*v[524]) - (v[1469] * v[191] + v[1147] * v[761])*v[9];
	v[1546] = -(i2328*v[524]) - (v[1469] * v[190] + v[1145] * v[761])*v[9];
	v[1547] = -(i2340*v[524]) - (v[1469] * v[189] + v[1143] * v[761])*v[9];
	v[1548] = -(i2336*v[524]) - (v[1469] * v[188] + v[1141] * v[761])*v[9];
	v[1549] = -(i2332*v[524]) - (v[1469] * v[187] + v[1139] * v[761])*v[9];
	v[1550] = -(i2339*v[524]) - (v[1469] * v[186] + v[1137] * v[761])*v[9];
	v[1551] = -(i2335*v[524]) - (v[1469] * v[185] + v[1135] * v[761])*v[9];
	v[1552] = -(i2331*v[524]) - (v[1469] * v[184] + v[1133] * v[761])*v[9];
	v[1553] = -(i2327*v[524]) - (v[1469] * v[183] + v[1131] * v[761])*v[9];
	v[1554] = -(i2326*v[524]) - (v[1469] * v[182] + v[1129] * v[761])*v[9];
	v[1555] = -(i2325*v[524]) - (v[1469] * v[181] + v[1127] * v[761])*v[9];
	v[1575] = -(i2342*v[523]) - (v[1468] * v[198] + v[1163] * v[760])*v[9];
	v[1576] = -(i2338*v[523]) - (v[1468] * v[197] + v[1160] * v[760])*v[9];
	v[1577] = -(i2334*v[523]) - (v[1468] * v[196] + v[1157] * v[760])*v[9];
	v[1578] = -(i2341*v[523]) - (v[1468] * v[195] + v[1155] * v[760])*v[9];
	v[1579] = -(i2337*v[523]) - (v[1468] * v[194] + v[1153] * v[760])*v[9];
	v[1580] = -(i2333*v[523]) - (v[1468] * v[193] + v[1151] * v[760])*v[9];
	v[1581] = -(i2330*v[523]) - (v[1468] * v[192] + v[1149] * v[760])*v[9];
	v[1582] = -(i2329*v[523]) - (v[1468] * v[191] + v[1147] * v[760])*v[9];
	v[1583] = -(i2328*v[523]) - (v[1468] * v[190] + v[1145] * v[760])*v[9];
	v[1584] = -(i2340*v[523]) - (v[1468] * v[189] + v[1143] * v[760])*v[9];
	v[1585] = -(i2336*v[523]) - (v[1468] * v[188] + v[1141] * v[760])*v[9];
	v[1586] = -(i2332*v[523]) - (v[1468] * v[187] + v[1139] * v[760])*v[9];
	v[1587] = -(i2339*v[523]) - (v[1468] * v[186] + v[1137] * v[760])*v[9];
	v[1588] = -(i2335*v[523]) - (v[1468] * v[185] + v[1135] * v[760])*v[9];
	v[1589] = -(i2331*v[523]) - (v[1468] * v[184] + v[1133] * v[760])*v[9];
	v[1590] = -(i2327*v[523]) - (v[1468] * v[183] + v[1131] * v[760])*v[9];
	v[1591] = -(i2326*v[523]) - (v[1468] * v[182] + v[1129] * v[760])*v[9];
	v[1592] = -(i2325*v[523]) - (v[1468] * v[181] + v[1127] * v[760])*v[9];
	v[1611] = v[1470] * v[4];
	v[1612] = v[1469] * v[4];
	v[1613] = v[1468] * v[4];
	v[1614] = v[1467];
	v[1615] = v[1467];
	v[1616] = v[1466];
	v[1617] = v[1466];
	v[1618] = v[1465];
	v[1619] = v[1465];
	b1620 = b483;
	if (b1620) {
		v[1621] = 0e0;
		v[1622] = 0e0;
		v[1623] = 0e0;
		b1624 = b502;
		if (b1624) {
			v[1623] = v[1614];
			v[1614] = 0e0;
			v[1622] = v[1616];
			v[1616] = 0e0;
			v[1621] = v[1618];
			v[1618] = 0e0;
		}
		else {
			v[1615] = 0e0;
			v[1614] = 0e0;
			v[1617] = 0e0;
			v[1616] = 0e0;
			v[1619] = 0e0;
			v[1618] = 0e0;
		};
		v[1635] = v[1621] * v[457];
		v[1633] = v[1622] * v[471];
		v[1631] = v[1623] * v[473];
		b1628 = b485;
		if (b1628) {
			v[1436] = v[1436] + v[1631];
			v[1434] = v[1434] + v[1623] * v[497];
			v[1436] = v[1436] + v[1633];
			v[1433] = v[1433] + v[1622] * v[497];
			v[1436] = v[1436] + v[1635];
			v[1432] = v[1432] + v[1621] * v[497];
			v[1438] = v[1438] + v[1436] * v[2284];
			v[1435] = v[1435] + (v[1438] * v[1629] * v[1630]) / 2e0;
		}
		else {
			v[1437] = v[1437] + v[1631];
			v[1434] = v[1632] + v[1623] * v[501];
			v[1437] = v[1437] + v[1633];
			v[1433] = v[1634] + v[1622] * v[501];
			v[1437] = v[1437] + v[1635];
			v[1432] = v[1636] + v[1621] * v[501];
			v[1439] = v[1439] + v[1437] * v[2284];
			v[1166] = v[1166] + (v[1439] * v[1637] * v[1638]) / 2e0;
		};
	}
	else {
	};
	v[2393] = v[1432] * v[456];
	v[2394] = v[1433] * v[464];
	v[1650] = v[1166];
	v[1648] = v[1619];
	v[1646] = v[1617];
	v[1644] = v[1615];
	b1639 = b483;
	if (b1639) {
		v[2375] = v[1619] * v[228];
		v[2374] = v[1617] * v[229];
		v[2373] = v[1615] * v[230];
		b1640 = b485;
		if (b1640) {
			v[1180] = v[1180] + v[1615] * v[2282];
			v[1615] = 0e0;
			v[1201] = v[1201] + v[1617] * v[2282];
			v[1617] = 0e0;
			v[1200] = v[1200] + v[1619] * v[2282];
			v[1619] = 0e0;
			v[1435] = v[1435] + (-v[2373] - v[2374] - v[2375])*v[31] * v[34] * Power(v[494], v[495]);
			v[1166] = v[1166] - v[1435];
		}
		else {
			v[1180] = v[1645] + v[1644] * v[491];
			v[1615] = 0e0;
			v[1201] = v[1647] + v[1646] * v[491];
			v[1617] = 0e0;
			v[1200] = v[1649] + v[1648] * v[491];
			v[1619] = 0e0;
			v[1166] = v[1650] + v[2291] * (v[2373] + v[2374] + v[2375])*Power(v[482], v[500]);
		};
	}
	else {
	};
	v[1200] = v[1200] + v[1422] * v[1434];
	v[1201] = v[1201] + v[1415] * v[1434];
	v[1651] = v[1178] * v[1434] + v[2345] * v[848];
	v[1652] = v[1423] + v[1434] * v[228];
	v[2408] = v[1652] * v[230] + v[2393] + v[1221] * v[846];
	v[1658] = v[1416] + v[1434] * v[229];
	v[2409] = v[1658] * v[230] + v[2394] + v[1193] * v[847];
	v[1682] = v[1358] * v[1433] + v[2376] * v[847];
	v[1201] = v[1201] + v[1433] * v[956];
	v[1707] = v[1343] * v[1432] + v[2377] * v[846];
	v[1708] = v[1433] * v[181] + v[1432] * v[182] + v[1129] * v[846] + v[1127] * v[847];
	v[1709] = v[1433] * v[184] + v[1432] * v[185] + v[1135] * v[846] + v[1133] * v[847];
	v[1710] = v[1433] * v[187] + v[1432] * v[188] + v[1141] * v[846] + v[1139] * v[847];
	v[1711] = v[1433] * v[190] + v[1432] * v[191] + v[1147] * v[846] + v[1145] * v[847];
	v[1712] = v[1433] * v[193] + v[1432] * v[194] + v[1153] * v[846] + v[1151] * v[847];
	v[1713] = v[1433] * v[196] + v[1432] * v[197] + v[1160] * v[846] + v[1157] * v[847];
	v[2392] = v[125] * v[1708] + v[127] * v[1709] + v[129] * v[1710] - v[153] * v[1711] - v[155] * v[1712] - v[157] * v[1713];
	v[1200] = v[1200] + v[1432] * v[956];
	v[1738] = 0e0;
	v[1739] = 0e0;
	v[1740] = 0e0;
	v[1741] = 0e0;
	v[1742] = 0e0;
	v[1743] = 0e0;
	v[1744] = 0e0;
	v[1745] = 0e0;
	v[1746] = 0e0;
	b1747 = b30;
	if (b1747) {
		v[2380] = v[1613] * v[228];
		v[2379] = v[1612] * v[229];
		v[2383] = v[2379] + v[2380];
		v[2378] = v[1611] * v[230];
		v[2382] = v[2378] + v[2380];
		v[2381] = v[2378] + v[2379];
		v[1299] = v[1299] + v[1611] * v[453];
		v[1298] = v[1298] + v[1612] * v[452];
		v[1294] = v[1294] - v[228] * v[2381] + v[1613] * v[991];
		v[1295] = v[1295] - v[229] * v[2382] + v[1612] * v[993];
		v[1296] = v[1296] - v[230] * v[2383] + v[1611] * v[995];
		v[1297] = v[1297] + v[1613] * v[451];
		v[1200] = v[1200] - v[1613] * v[2303] - v[2381] * v[451];
		v[1201] = v[1201] + v[1612] * v[2304] - v[2382] * v[452];
		v[1180] = v[1180] + v[1611] * v[2305] - v[2383] * v[453];
		v[1738] = v[1] * v[1294];
		v[1739] = v[1294] * v[2];
		v[1740] = v[1294] * v[3];
		v[1748] = -(v[1294] * v[162]);
		v[1749] = -(v[1294] * v[160]);
		v[1750] = -(v[1294] * v[158]);
		v[1751] = v[1294] * v[134];
		v[1752] = v[1294] * v[132];
		v[1753] = v[1294] * v[130];
		v[1741] = v[1] * v[1295];
		v[1742] = v[1295] * v[2];
		v[1743] = v[1295] * v[3];
		v[1754] = -(v[1295] * v[162]);
		v[1755] = -(v[1295] * v[160]);
		v[1756] = -(v[1295] * v[158]);
		v[1757] = v[1295] * v[134];
		v[1758] = v[1295] * v[132];
		v[1759] = v[1295] * v[130];
		v[1744] = v[1] * v[1296];
		v[1745] = v[1296] * v[2];
		v[1746] = v[1296] * v[3];
		v[1760] = -(v[1296] * v[162]);
		v[1761] = -(v[1296] * v[160]);
		v[1762] = -(v[1296] * v[158]);
		v[1763] = v[1296] * v[134];
		v[1764] = v[1296] * v[132];
		v[1765] = v[1296] * v[130];
		v[1707] = -v[1297] + v[1707];
		v[1682] = -v[1298] + v[1682];
		v[1651] = -v[1299] + v[1651];
	}
	else {
		v[1753] = 0e0;
		v[1752] = 0e0;
		v[1751] = 0e0;
		v[1759] = 0e0;
		v[1758] = 0e0;
		v[1757] = 0e0;
		v[1765] = 0e0;
		v[1764] = 0e0;
		v[1763] = 0e0;
		v[1750] = 0e0;
		v[1749] = 0e0;
		v[1748] = 0e0;
		v[1756] = 0e0;
		v[1755] = 0e0;
		v[1754] = 0e0;
		v[1762] = 0e0;
		v[1761] = 0e0;
		v[1760] = 0e0;
	};
	b1766 = b416;
	if (b1766) {
		v[1232] = v[1232] + (v[1019] * v[1746]) / 2e0;
		v[1242] = v[1242] + v[1746] * v[2278];
		v[1232] = v[1232] + v[1023] * v[1745];
		v[1243] = v[1243] + v[1745] * v[433];
		v[1232] = v[1232] + v[1027] * v[1744];
		v[1244] = v[1244] + v[1744] * v[433];
		v[1232] = v[1232] + v[1031] * v[1743];
		v[1245] = v[1245] + v[1743] * v[433];
		v[1232] = v[1232] + (v[1036] * v[1742]) / 2e0;
		v[1246] = v[1246] + v[1742] * v[2278];
		v[1232] = v[1232] + v[1040] * v[1741];
		v[1247] = v[1247] + v[1741] * v[433];
		v[1232] = v[1232] + v[1044] * v[1740];
		v[1248] = v[1248] + v[1740] * v[433];
		v[1232] = v[1232] + v[1049] * v[1739];
		v[1249] = v[1249] + v[1739] * v[433];
		v[1232] = v[1232] + (v[1054] * v[1738]) / 2e0;
		v[1250] = v[1250] + v[1738] * v[2278];
		v[1251] = v[1251] + v[1232] * v[2384];
		v[2389] = -v[1246] + v[1251];
		v[1230] = v[1230] - v[1244];
		v[1776] = v[1244] + v[1248];
		v[1230] = v[1230] + v[1248];
		v[1229] = v[1229] + v[1243] + (v[1776] * v[432]) / 2e0;
		v[1778] = v[1243] + v[1245];
		v[1229] = v[1229] - v[1245];
		v[1231] = v[1231] + v[1247] + v[1778] * v[2277] + v[1776] * v[2306] + v[2385] * (-v[1250] + v[2389]);
		v[1231] = v[1231] - v[1249];
		v[2386] = v[1231] * v[420];
		v[1780] = v[1247] + v[1249];
		v[1228] = v[1228] + v[2386] * v[423];
		v[1226] = v[1226] + v[2386] * v[429];
		v[1224] = v[1224] + v[1231] * v[2276];
		v[1230] = v[1230] + (v[1780] * v[430] - 4e0*(v[1242] + v[1250] - v[1251])*v[431] + v[1778] * v[432]) / 2e0;
		v[2387] = v[1230] * v[419];
		v[1228] = v[1228] + v[2387] * v[423];
		v[1226] = v[1226] + v[2387] * v[429];
		v[1223] = v[1223] + v[1230] * v[2276];
		v[1229] = v[1229] + v[1780] * v[2277] + v[2388] * (-v[1242] + v[2389]);
		v[2390] = v[1229] * v[418];
		v[1228] = v[1228] + v[2390] * v[423];
		v[1226] = v[1226] + v[2390] * v[429];
		v[1222] = v[1222] + v[1229] * v[2276];
		v[1227] = v[1227] + v[1228] * v[428];
		v[1225] = v[1225] + v[1227];
		v[1252] = v[1252] + 2e0*v[1226] * v[1267];
		v[1253] = v[1253] + (v[1252] * v[1268]) / 2e0;
		v[1225] = v[1225] + v[1253];
		v[2391] = v[1225] / v[421];
		v[1224] = v[1224] + v[2391] * v[420];
		v[1223] = v[1223] + v[2391] * v[419];
		v[1222] = v[1222] + v[2391] * v[418];
		v[1200] = v[1200] - v[1224] * v[224];
		v[1201] = v[1201] + v[1224] * v[223];
		v[1200] = v[1200] + v[1223] * v[225];
		v[1180] = v[1180] - v[1223] * v[223];
		v[1201] = v[1201] - v[1222] * v[225];
		v[1180] = v[1180] + v[1222] * v[224];
	}
	else {
	};
	v[2400] = (-(v[1592] * v[398]) - v[1591] * v[399] - v[1590] * v[400] - v[1589] * v[401] - v[1588] * v[402]
		- v[1587] * v[403] - v[1586] * v[404] - v[1585] * v[405] - v[1584] * v[406] - v[1583] * v[407] - v[1582] * v[408]
		- v[1581] * v[409] - v[1580] * v[410] - v[1579] * v[411] - v[1578] * v[412] - v[1577] * v[413] - v[1576] * v[414]
		- v[1575] * v[415]) / 2e0;
	v[2401] = (-(v[1592] * v[380]) - v[1591] * v[381] - v[1590] * v[382] - v[1589] * v[383] - v[1588] * v[384]
		- v[1587] * v[385] - v[1586] * v[386] - v[1585] * v[387] - v[1584] * v[388] - v[1583] * v[389] - v[1582] * v[390]
		- v[1581] * v[391] - v[1580] * v[392] - v[1579] * v[393] - v[1578] * v[394] - v[1577] * v[395] - v[1576] * v[396]
		- v[1575] * v[397]) / 2e0;
	v[2406] = (v[1592] * v[362] + v[1591] * v[363] + v[1590] * v[364] + v[1589] * v[365] + v[1588] * v[366] + v[1587] * v[367]
		+ v[1586] * v[368] + v[1585] * v[369] + v[1584] * v[370] + v[1583] * v[371] + v[1582] * v[372] + v[1581] * v[373]
		+ v[1580] * v[374] + v[1579] * v[375] + v[1578] * v[376] + v[1577] * v[377] + v[1576] * v[378] + v[1575] * v[379]) / 2e0;
	v[2407] = (v[1592] * v[344] + v[1591] * v[345] + v[1590] * v[346] + v[1589] * v[347] + v[1588] * v[348] + v[1587] * v[349]
		+ v[1586] * v[350] + v[1585] * v[351] + v[1584] * v[352] + v[1583] * v[353] + v[1582] * v[354] + v[1581] * v[355]
		+ v[1580] * v[356] + v[1579] * v[357] + v[1578] * v[358] + v[1577] * v[359] + v[1576] * v[360] + v[1575] * v[361]) / 2e0;
	v[2398] = (-(v[1555] * v[398]) - v[1554] * v[399] - v[1553] * v[400] - v[1552] * v[401] - v[1551] * v[402]
		- v[1550] * v[403] - v[1549] * v[404] - v[1548] * v[405] - v[1547] * v[406] - v[1546] * v[407] - v[1545] * v[408]
		- v[1544] * v[409] - v[1543] * v[410] - v[1542] * v[411] - v[1541] * v[412] - v[1540] * v[413] - v[1539] * v[414]
		- v[1538] * v[415]) / 2e0;
	v[2399] = (-(v[1555] * v[380]) - v[1554] * v[381] - v[1553] * v[382] - v[1552] * v[383] - v[1551] * v[384]
		- v[1550] * v[385] - v[1549] * v[386] - v[1548] * v[387] - v[1547] * v[388] - v[1546] * v[389] - v[1545] * v[390]
		- v[1544] * v[391] - v[1543] * v[392] - v[1542] * v[393] - v[1541] * v[394] - v[1540] * v[395] - v[1539] * v[396]
		- v[1538] * v[397]) / 2e0;
	v[2404] = (v[1555] * v[362] + v[1554] * v[363] + v[1553] * v[364] + v[1552] * v[365] + v[1551] * v[366] + v[1550] * v[367]
		+ v[1549] * v[368] + v[1548] * v[369] + v[1547] * v[370] + v[1546] * v[371] + v[1545] * v[372] + v[1544] * v[373]
		+ v[1543] * v[374] + v[1542] * v[375] + v[1541] * v[376] + v[1540] * v[377] + v[1539] * v[378] + v[1538] * v[379]) / 2e0;
	v[2405] = (v[1555] * v[344] + v[1554] * v[345] + v[1553] * v[346] + v[1552] * v[347] + v[1551] * v[348] + v[1550] * v[349]
		+ v[1549] * v[350] + v[1548] * v[351] + v[1547] * v[352] + v[1546] * v[353] + v[1545] * v[354] + v[1544] * v[355]
		+ v[1543] * v[356] + v[1542] * v[357] + v[1541] * v[358] + v[1540] * v[359] + v[1539] * v[360] + v[1538] * v[361]) / 2e0;
	v[2396] = (-(v[1518] * v[398]) - v[1517] * v[399] - v[1516] * v[400] - v[1515] * v[401] - v[1514] * v[402]
		- v[1513] * v[403] - v[1512] * v[404] - v[1511] * v[405] - v[1510] * v[406] - v[1509] * v[407] - v[1508] * v[408]
		- v[1507] * v[409] - v[1506] * v[410] - v[1505] * v[411] - v[1504] * v[412] - v[1503] * v[413] - v[1502] * v[414]
		- v[1501] * v[415]) / 2e0;
	v[2397] = (-(v[1518] * v[380]) - v[1517] * v[381] - v[1516] * v[382] - v[1515] * v[383] - v[1514] * v[384]
		- v[1513] * v[385] - v[1512] * v[386] - v[1511] * v[387] - v[1510] * v[388] - v[1509] * v[389] - v[1508] * v[390]
		- v[1507] * v[391] - v[1506] * v[392] - v[1505] * v[393] - v[1504] * v[394] - v[1503] * v[395] - v[1502] * v[396]
		- v[1501] * v[397]) / 2e0;
	v[2402] = (v[1518] * v[362] + v[1517] * v[363] + v[1516] * v[364] + v[1515] * v[365] + v[1514] * v[366] + v[1513] * v[367]
		+ v[1512] * v[368] + v[1511] * v[369] + v[1510] * v[370] + v[1509] * v[371] + v[1508] * v[372] + v[1507] * v[373]
		+ v[1506] * v[374] + v[1505] * v[375] + v[1504] * v[376] + v[1503] * v[377] + v[1502] * v[378] + v[1501] * v[379]) / 2e0;
	v[2403] = (v[1518] * v[344] + v[1517] * v[345] + v[1516] * v[346] + v[1515] * v[347] + v[1514] * v[348] + v[1513] * v[349]
		+ v[1512] * v[350] + v[1511] * v[351] + v[1510] * v[352] + v[1509] * v[353] + v[1508] * v[354] + v[1507] * v[355]
		+ v[1506] * v[356] + v[1505] * v[357] + v[1504] * v[358] + v[1503] * v[359] + v[1502] * v[360] + v[1501] * v[361]) / 2e0;
	v[1200] = v[1200] + v[1707] * v[2350];
	v[1200] = v[1200] + v[229] * v[2392];
	v[1201] = v[1201] + v[1682] * v[2347] + v[228] * v[2392];
	v[1798] = v[1357] + v[1432] * v[228] + v[1433] * v[229];
	v[2395] = v[1798] * v[230] + v[1434] * v[472] + v[1177] * v[848];
	v[1180] = v[1180] + v[1178] * v[1798] + v[1651] * v[2344];
	v[1180] = v[1180] + v[1343] * v[1652] + v[1358] * v[1658];
	b1799 = b226;
	if (b1799) {
		v[1800] = -v[1180];
		v[1801] = -v[1201];
		v[1802] = -v[1200];
	}
	else {
		v[1800] = v[1180];
		v[1801] = v[1201];
		v[1802] = v[1200];
	};
	v[1166] = v[1166] + v[1079] * v[1168] * v[211] + (v[1074] * v[1158] + v[1073] * v[1161] + v[1072] * v[1164]
		+ v[1802] * v[199] + v[1801] * v[200] + v[1800] * v[201])*v[212];
	v[1811] = v[1072] * v[1167] + v[1800] * v[213] + (v[1166] * v[201] + v[1164] * v[849]) / v[482];
	v[1813] = v[1073] * v[1167] + v[1801] * v[213] + (v[1166] * v[200] + v[1161] * v[849]) / v[482];
	v[1815] = v[1074] * v[1167] + v[1802] * v[213] + (v[1166] * v[199] + v[1158] * v[849]) / v[482];
	v[1816] = v[1399] - v[177] * v[1811] - v[173] * v[1813] - v[169] * v[1815] + v[1711] * v[2275] - v[190] * v[2393]
		- v[191] * v[2394] - v[192] * v[2395] - v[230] * (v[1652] * v[190] + v[1658] * v[191] + v[1374] * v[848]);
	v[1817] = v[1402] + v[149] * v[1811] + v[145] * v[1813] + v[141] * v[1815] - v[1708] * v[2275] + v[181] * v[2393]
		+ v[182] * v[2394] + v[183] * v[2395] + v[230] * (v[1652] * v[181] + v[1658] * v[182] + v[1371] * v[848]);
	v[1818] = (v[1397] - v[179] * v[1811] - v[175] * v[1813] - v[171] * v[1815] - v[1816] + v[1713] * v[2275]
		- v[196] * v[2393] - v[197] * v[2394] - v[198] * v[2395] - v[230] * (v[1652] * v[196] + v[1658] * v[197] + v[1376] * v[848]
			)) / 2e0;
	v[1819] = (v[1398] - v[178] * v[1811] - v[174] * v[1813] - v[170] * v[1815] - v[1816] + v[1712] * v[2275]
		- v[193] * v[2393] - v[194] * v[2394] - v[195] * v[2395] - v[230] * (v[1652] * v[193] + v[1658] * v[194] + v[1375] * v[848]
			)) / 2e0;
	v[1760] = v[1760] - v[157] * v[1811] + v[2396];
	v[1754] = v[1754] - v[157] * v[1813] + v[2398];
	v[1748] = v[1748] - v[157] * v[1815] + v[2400];
	v[1761] = v[1761] - v[155] * v[1811] + v[2397];
	v[1755] = v[1755] - v[155] * v[1813] + v[2399];
	v[1749] = v[1749] - v[155] * v[1815] + v[2401];
	v[1762] = v[1762] - v[153] * v[1811] - v[2396] - v[2397];
	v[1756] = v[1756] - v[153] * v[1813] - v[2398] - v[2399];
	v[1750] = v[1750] - v[153] * v[1815] - v[2400] - v[2401];
	v[1829] = (v[1400] + v[151] * v[1811] + v[147] * v[1813] + v[143] * v[1815] - v[1817] - v[1710] * v[2275]
		+ v[187] * v[2393] + v[188] * v[2394] + v[189] * v[2395] + v[230] * (v[1652] * v[187] + v[1658] * v[188] + v[1373] * v[848]
			)) / 2e0;
	v[1830] = (v[1401] + v[150] * v[1811] + v[146] * v[1813] + v[142] * v[1815] - v[1817] - v[1709] * v[2275]
		+ v[184] * v[2393] + v[185] * v[2394] + v[186] * v[2395] + v[230] * (v[1652] * v[184] + v[1658] * v[185] + v[1372] * v[848]
			)) / 2e0;
	v[1763] = v[1763] + v[129] * v[1811] + v[2402];
	v[1757] = v[1757] + v[129] * v[1813] + v[2404];
	v[1751] = v[1751] + v[129] * v[1815] + v[2406];
	v[1764] = v[1764] + v[127] * v[1811] + v[2403];
	v[1758] = v[1758] + v[127] * v[1813] + v[2405];
	v[1752] = v[1752] + v[127] * v[1815] + v[2407];
	v[1765] = v[1765] + v[125] * v[1811] - v[2402] - v[2403];
	v[1759] = v[1759] + v[125] * v[1813] - v[2404] - v[2405];
	v[1753] = v[1753] + v[125] * v[1815] - v[2406] - v[2407];
	v[4155] = v[1753] + v[1830] * v[344] + v[1829] * v[362] + v[1819] * v[380] + v[1818] * v[398] + v[10] * (v[1395] * v[228]
		+ v[125] * v[2408] - v[1452] * v[344] - v[1456] * v[362] + v[1457] * v[380] + v[1461] * v[398] + v[1433] * v[458]
		+ v[1213] * v[847] + (-(v[1468] * v[685]) - v[1469] * v[686] - v[1470] * v[687])*v[9]);
	v[4156] = v[1759] + v[1830] * v[345] + v[1829] * v[363] + v[1819] * v[381] + v[1818] * v[399] + v[10] * (v[1395] * v[229]
		+ v[125] * v[2409] - v[1452] * v[345] - v[1456] * v[363] + v[1457] * v[381] + v[1461] * v[399] + v[1432] * v[458]
		+ v[1213] * v[846] + (-(v[1468] * v[689]) - v[1469] * v[690] - v[1470] * v[691])*v[9]);
	v[4157] = v[1765] + v[1830] * v[346] + v[1829] * v[364] + v[1819] * v[382] + v[1818] * v[400] + v[10] * (v[125] * (
		-v[2346] + v[2395]) - v[1452] * v[346] - v[1456] * v[364] + v[1457] * v[382] + v[1461] * v[400] + (-(v[1468] * v[693])
			- v[1469] * v[694] - v[1470] * v[695])*v[9]);
	v[4158] = v[1752] + v[1830] * v[347] + v[1829] * v[365] + v[1819] * v[383] + v[1818] * v[401] + v[10] * (v[1392] * v[228]
		+ v[127] * v[2408] - v[1452] * v[347] - v[1456] * v[365] + v[1457] * v[383] + v[1461] * v[401] + v[1433] * v[459]
		+ v[1212] * v[847] + (-(v[1468] * v[697]) - v[1469] * v[698] - v[1470] * v[699])*v[9]);
	v[4159] = v[1758] + v[1830] * v[348] + v[1829] * v[366] + v[1819] * v[384] + v[1818] * v[402] + v[10] * (v[1392] * v[229]
		+ v[127] * v[2409] - v[1452] * v[348] - v[1456] * v[366] + v[1457] * v[384] + v[1461] * v[402] + v[1432] * v[459]
		+ v[1212] * v[846] + (-(v[1468] * v[701]) - v[1469] * v[702] - v[1470] * v[703])*v[9]);
	v[4160] = v[1764] + v[1830] * v[349] + v[1829] * v[367] + v[1819] * v[385] + v[1818] * v[403] + v[10] * (v[127] * (
		-v[2346] + v[2395]) - v[1452] * v[349] - v[1456] * v[367] + v[1457] * v[385] + v[1461] * v[403] + (-(v[1468] * v[705])
			- v[1469] * v[706] - v[1470] * v[707])*v[9]);
	v[4161] = v[1751] + v[1830] * v[350] + v[1829] * v[368] + v[1819] * v[386] + v[1818] * v[404] + v[10] * (v[1389] * v[228]
		+ v[129] * v[2408] - v[1452] * v[350] - v[1456] * v[368] + v[1457] * v[386] + v[1461] * v[404] + v[1433] * v[460]
		+ v[1211] * v[847] + (-(v[1468] * v[709]) - v[1469] * v[710] - v[1470] * v[711])*v[9]);
	v[4162] = v[1757] + v[1830] * v[351] + v[1829] * v[369] + v[1819] * v[387] + v[1818] * v[405] + v[10] * (v[1389] * v[229]
		+ v[129] * v[2409] - v[1452] * v[351] - v[1456] * v[369] + v[1457] * v[387] + v[1461] * v[405] + v[1432] * v[460]
		+ v[1211] * v[846] + (-(v[1468] * v[713]) - v[1469] * v[714] - v[1470] * v[715])*v[9]);
	v[4163] = v[1763] + v[1830] * v[352] + v[1829] * v[370] + v[1819] * v[388] + v[1818] * v[406] + v[10] * (v[129] * (
		-v[2346] + v[2395]) - v[1452] * v[352] - v[1456] * v[370] + v[1457] * v[388] + v[1461] * v[406] + (-(v[1468] * v[717])
			- v[1469] * v[718] - v[1470] * v[719])*v[9]);
	v[4164] = v[1750] + v[1830] * v[353] + v[1829] * v[371] + v[1819] * v[389] + v[1818] * v[407] + v[10] * (v[1386] * v[228]
		- v[153] * v[2408] - v[1452] * v[353] - v[1456] * v[371] + v[1457] * v[389] + v[1461] * v[407] + v[1433] * v[461]
		+ v[1210] * v[847] + (-(v[1468] * v[721]) - v[1469] * v[722] - v[1470] * v[723])*v[9]);
	v[4165] = v[1756] + v[1830] * v[354] + v[1829] * v[372] + v[1819] * v[390] + v[1818] * v[408] + v[10] * (v[1386] * v[229]
		- v[153] * v[2409] - v[1452] * v[354] - v[1456] * v[372] + v[1457] * v[390] + v[1461] * v[408] + v[1432] * v[461]
		+ v[1210] * v[846] + (-(v[1468] * v[725]) - v[1469] * v[726] - v[1470] * v[727])*v[9]);
	v[4166] = v[1762] + v[1830] * v[355] + v[1829] * v[373] + v[1819] * v[391] + v[1818] * v[409] + v[10] * (v[153] * (v[2346]
		- v[2395]) - v[1452] * v[355] - v[1456] * v[373] + v[1457] * v[391] + v[1461] * v[409] + (-(v[1468] * v[729])
			- v[1469] * v[730] - v[1470] * v[731])*v[9]);
	v[4167] = v[1749] + v[1830] * v[356] + v[1829] * v[374] + v[1819] * v[392] + v[1818] * v[410] + v[10] * (v[1383] * v[228]
		- v[155] * v[2408] - v[1452] * v[356] - v[1456] * v[374] + v[1457] * v[392] + v[1461] * v[410] + v[1433] * v[462]
		+ v[1209] * v[847] + (-(v[1468] * v[733]) - v[1469] * v[734] - v[1470] * v[735])*v[9]);
	v[4168] = v[1755] + v[1830] * v[357] + v[1829] * v[375] + v[1819] * v[393] + v[1818] * v[411] + v[10] * (v[1383] * v[229]
		- v[155] * v[2409] - v[1452] * v[357] - v[1456] * v[375] + v[1457] * v[393] + v[1461] * v[411] + v[1432] * v[462]
		+ v[1209] * v[846] + (-(v[1468] * v[737]) - v[1469] * v[738] - v[1470] * v[739])*v[9]);
	v[4169] = v[1761] + v[1830] * v[358] + v[1829] * v[376] + v[1819] * v[394] + v[1818] * v[412] + v[10] * (v[155] * (v[2346]
		- v[2395]) - v[1452] * v[358] - v[1456] * v[376] + v[1457] * v[394] + v[1461] * v[412] + (-(v[1468] * v[741])
			- v[1469] * v[742] - v[1470] * v[743])*v[9]);
	v[4170] = v[1748] + v[1830] * v[359] + v[1829] * v[377] + v[1819] * v[395] + v[1818] * v[413] + v[10] * (v[1380] * v[228]
		- v[157] * v[2408] - v[1452] * v[359] - v[1456] * v[377] + v[1457] * v[395] + v[1461] * v[413] + v[1433] * v[463]
		+ v[1208] * v[847] + (-(v[1468] * v[745]) - v[1469] * v[746] - v[1470] * v[747])*v[9]);
	v[4171] = v[1754] + v[1830] * v[360] + v[1829] * v[378] + v[1819] * v[396] + v[1818] * v[414] + v[10] * (v[1380] * v[229]
		- v[157] * v[2409] - v[1452] * v[360] - v[1456] * v[378] + v[1457] * v[396] + v[1461] * v[414] + v[1432] * v[463]
		+ v[1208] * v[846] + (-(v[1468] * v[749]) - v[1469] * v[750] - v[1470] * v[751])*v[9]);
	v[4172] = v[1760] + v[1830] * v[361] + v[1829] * v[379] + v[1819] * v[397] + v[1818] * v[415] + v[10] * (v[157] * (v[2346]
		- v[2395]) - v[1452] * v[361] - v[1456] * v[379] + v[1457] * v[397] + v[1461] * v[415] + (-(v[1468] * v[753])
			- v[1469] * v[754] - v[1470] * v[755])*v[9]);
	Rc[i758 - 1] += v[3808 + i758];
	for (i1107 = 1; i1107 <= 18; i1107++) {
		Kc[i758 - 1][i1107 - 1] += v[4154 + i1107];
	};/* end for */
};/* end for */
v[1849] = 0e0;
v[1850] = 0e0;
v[1851] = 0e0;
v[1852] = 0e0;
b1853 = b483;
if (b1853) {
	b1854 = b502;
	b1858 = b485;
	if (b1858) {
		v[1851] = 0e0;
		v[1850] = 0e0;
		v[1849] = 0e0;
		v[1852] = 0e0;
	}
	else {
	};
}
else {
};
v[2412] = -(v[1849] * v[456]);
v[2411] = -(v[1850] * v[464]);
v[2410] = -(v[1851] * v[472]);
v[1859] = 0e0;
v[1860] = 0e0;
v[1861] = 0e0;
v[1862] = 0e0;
b1863 = b483;
if (b1863) {
	b1864 = b485;
	if (b1864) {
		v[1861] = 0e0;
		v[1860] = 0e0;
		v[1859] = 0e0;
		v[1862] = -v[1852];
	}
	else {
	};
}
else {
};
v[1859] = v[1422] * v[1851] + v[1859];
v[1860] = v[1415] * v[1851] + v[1860];
v[1865] = v[1178] * v[1851];
v[1866] = v[1851] * v[228];
v[2415] = v[1866] * v[230] - v[2412];
v[1872] = v[1851] * v[229];
v[2416] = v[1872] * v[230] - v[2411];
v[1896] = v[1358] * v[1850];
v[1860] = v[1860] + v[1850] * v[956];
v[1921] = v[1343] * v[1849];
v[1922] = v[182] * v[1849] + v[181] * v[1850];
v[1923] = v[1849] * v[185] + v[184] * v[1850];
v[1924] = v[1850] * v[187] + v[1849] * v[188];
v[1925] = v[1850] * v[190] + v[1849] * v[191];
v[1926] = v[1850] * v[193] + v[1849] * v[194];
v[1927] = v[1850] * v[196] + v[1849] * v[197];
v[2084] = v[125] * v[1922] + v[127] * v[1923] + v[129] * v[1924] - v[153] * v[1925] - v[155] * v[1926] - v[157] * v[1927];
v[1859] = v[1859] + v[1849] * v[956];
v[1928] = v[1849] * v[228] + v[1850] * v[229];
v[4415] = v[1866];
v[4416] = v[1872];
v[4417] = v[1928];
v[4418] = 0e0;
v[4419] = 0e0;
v[4420] = 0e0;
v[4421] = 0e0;
v[4422] = 0e0;
v[4423] = 0e0;
v[4424] = 0e0;
v[4425] = 0e0;
v[4426] = 0e0;
v[4427] = 0e0;
v[4428] = 0e0;
v[4429] = 0e0;
v[4430] = 0e0;
v[4431] = 0e0;
v[4432] = 0e0;
v[4379] = 0e0;
v[4380] = 0e0;
v[4381] = 0e0;
v[4382] = v[1866];
v[4383] = v[1872];
v[4384] = v[1928];
v[4385] = 0e0;
v[4386] = 0e0;
v[4387] = 0e0;
v[4388] = 0e0;
v[4389] = 0e0;
v[4390] = 0e0;
v[4391] = 0e0;
v[4392] = 0e0;
v[4393] = 0e0;
v[4394] = 0e0;
v[4395] = 0e0;
v[4396] = 0e0;
v[4343] = 0e0;
v[4344] = 0e0;
v[4345] = 0e0;
v[4346] = 0e0;
v[4347] = 0e0;
v[4348] = 0e0;
v[4349] = v[1866];
v[4350] = v[1872];
v[4351] = v[1928];
v[4352] = 0e0;
v[4353] = 0e0;
v[4354] = 0e0;
v[4355] = 0e0;
v[4356] = 0e0;
v[4357] = 0e0;
v[4358] = 0e0;
v[4359] = 0e0;
v[4360] = 0e0;
v[4307] = 0e0;
v[4308] = 0e0;
v[4309] = 0e0;
v[4310] = 0e0;
v[4311] = 0e0;
v[4312] = 0e0;
v[4313] = 0e0;
v[4314] = 0e0;
v[4315] = 0e0;
v[4316] = -v[1866];
v[4317] = -v[1872];
v[4318] = -v[1928];
v[4319] = 0e0;
v[4320] = 0e0;
v[4321] = 0e0;
v[4322] = 0e0;
v[4323] = 0e0;
v[4324] = 0e0;
v[4271] = 0e0;
v[4272] = 0e0;
v[4273] = 0e0;
v[4274] = 0e0;
v[4275] = 0e0;
v[4276] = 0e0;
v[4277] = 0e0;
v[4278] = 0e0;
v[4279] = 0e0;
v[4280] = 0e0;
v[4281] = 0e0;
v[4282] = 0e0;
v[4283] = -v[1866];
v[4284] = -v[1872];
v[4285] = -v[1928];
v[4286] = 0e0;
v[4287] = 0e0;
v[4288] = 0e0;
v[4235] = 0e0;
v[4236] = 0e0;
v[4237] = 0e0;
v[4238] = 0e0;
v[4239] = 0e0;
v[4240] = 0e0;
v[4241] = 0e0;
v[4242] = 0e0;
v[4243] = 0e0;
v[4244] = 0e0;
v[4245] = 0e0;
v[4246] = 0e0;
v[4247] = 0e0;
v[4248] = 0e0;
v[4249] = 0e0;
v[4250] = -v[1866];
v[4251] = -v[1872];
v[4252] = -v[1928];
v[4213] = v[125] * v[1866];
v[4214] = v[125] * v[1872];
v[4215] = v[125] * v[1928];
v[4216] = v[127] * v[1866];
v[4217] = v[127] * v[1872];
v[4218] = v[127] * v[1928];
v[4219] = v[129] * v[1866];
v[4220] = v[129] * v[1872];
v[4221] = v[129] * v[1928];
v[4222] = -(v[153] * v[1866]);
v[4223] = -(v[153] * v[1872]);
v[4224] = -(v[153] * v[1928]);
v[4225] = -(v[155] * v[1866]);
v[4226] = -(v[155] * v[1872]);
v[4227] = -(v[155] * v[1928]);
v[4228] = -(v[157] * v[1866]);
v[4229] = -(v[157] * v[1872]);
v[4230] = -(v[157] * v[1928]);
v[2418] = -(v[1928] * v[195]);
v[2417] = -(v[1928] * v[198]);
v[2414] = -(v[192] * v[1928]);
v[2060] = v[181] * v[1866] + v[182] * v[1872] + v[183] * v[1928];
v[2058] = v[184] * v[1866] + v[185] * v[1872] + v[186] * v[1928];
v[2056] = v[1866] * v[187] + v[1872] * v[188] + v[189] * v[1928];
v[2054] = -(v[1866] * v[190]) - v[1872] * v[191] + v[2414];
v[2052] = -(v[1866] * v[193]) - v[1872] * v[194] + v[2418];
v[2050] = -(v[1866] * v[196]) - v[1872] * v[197] + v[2417];
v[4195] = -(v[125] * v[2412]) + v[1850] * v[458];
v[4196] = -(v[125] * v[2411]) + v[1849] * v[458];
v[4197] = -(v[125] * v[2410]);
v[4198] = -(v[127] * v[2412]) + v[1850] * v[459];
v[4199] = -(v[127] * v[2411]) + v[1849] * v[459];
v[4200] = -(v[127] * v[2410]);
v[4201] = -(v[129] * v[2412]) + v[1850] * v[460];
v[4202] = -(v[129] * v[2411]) + v[1849] * v[460];
v[4203] = -(v[129] * v[2410]);
v[4204] = v[153] * v[2412] + v[1850] * v[461];
v[4205] = v[153] * v[2411] + v[1849] * v[461];
v[4206] = v[153] * v[2410];
v[4207] = v[155] * v[2412] + v[1850] * v[462];
v[4208] = v[155] * v[2411] + v[1849] * v[462];
v[4209] = v[155] * v[2410];
v[4210] = v[157] * v[2412] + v[1850] * v[463];
v[4211] = v[157] * v[2411] + v[1849] * v[463];
v[4212] = v[157] * v[2410];
v[1859] = v[1859] + v[2084] * v[229] + v[1921] * v[2350];
v[1860] = v[1860] + v[2084] * v[228] + v[1896] * v[2347];
v[1861] = v[1861] + v[1343] * v[1866] + v[1358] * v[1872] + v[1178] * v[1928] + v[1865] * v[2344];
b1952 = b226;
if (b1952) {
	v[1953] = -v[1861];
	v[1954] = -v[1860];
	v[1955] = -v[1859];
}
else {
	v[1953] = v[1861];
	v[1954] = v[1860];
	v[1955] = v[1859];
};
v[1960] = v[1955] * v[199] + v[1954] * v[200] + v[1953] * v[201];
v[1862] = v[1862] + v[1960] * v[212];
v[2413] = v[1862] / v[482];
v[1962] = v[1953] * v[213] + v[201] * v[2413] + v[506];
v[1963] = v[1954] * v[213] + v[200] * v[2413] + v[505];
v[1964] = v[1955] * v[213] + v[199] * v[2413] + v[504];
v[4177] = v[125] * v[1964];
v[4178] = v[125] * v[1963];
v[4179] = v[125] * v[1962];
v[4180] = v[127] * v[1964];
v[4181] = v[127] * v[1963];
v[4182] = v[127] * v[1962];
v[4183] = v[129] * v[1964];
v[4184] = v[129] * v[1963];
v[4185] = v[129] * v[1962];
v[4186] = -(v[153] * v[1964]);
v[4187] = -(v[153] * v[1963]);
v[4188] = -(v[153] * v[1962]);
v[4189] = -(v[155] * v[1964]);
v[4190] = -(v[155] * v[1963]);
v[4191] = -(v[155] * v[1962]);
v[4192] = -(v[157] * v[1964]);
v[4193] = -(v[157] * v[1963]);
v[4194] = -(v[157] * v[1962]);
v[4433] = v[1964];
v[4434] = v[1963];
v[4435] = v[1962];
v[4436] = 0e0;
v[4437] = 0e0;
v[4438] = 0e0;
v[4439] = 0e0;
v[4440] = 0e0;
v[4441] = 0e0;
v[4442] = 0e0;
v[4443] = 0e0;
v[4444] = 0e0;
v[4445] = 0e0;
v[4446] = 0e0;
v[4447] = 0e0;
v[4448] = 0e0;
v[4449] = 0e0;
v[4450] = 0e0;
v[4397] = 0e0;
v[4398] = 0e0;
v[4399] = 0e0;
v[4400] = v[1964];
v[4401] = v[1963];
v[4402] = v[1962];
v[4403] = 0e0;
v[4404] = 0e0;
v[4405] = 0e0;
v[4406] = 0e0;
v[4407] = 0e0;
v[4408] = 0e0;
v[4409] = 0e0;
v[4410] = 0e0;
v[4411] = 0e0;
v[4412] = 0e0;
v[4413] = 0e0;
v[4414] = 0e0;
v[4361] = 0e0;
v[4362] = 0e0;
v[4363] = 0e0;
v[4364] = 0e0;
v[4365] = 0e0;
v[4366] = 0e0;
v[4367] = v[1964];
v[4368] = v[1963];
v[4369] = v[1962];
v[4370] = 0e0;
v[4371] = 0e0;
v[4372] = 0e0;
v[4373] = 0e0;
v[4374] = 0e0;
v[4375] = 0e0;
v[4376] = 0e0;
v[4377] = 0e0;
v[4378] = 0e0;
v[4325] = 0e0;
v[4326] = 0e0;
v[4327] = 0e0;
v[4328] = 0e0;
v[4329] = 0e0;
v[4330] = 0e0;
v[4331] = 0e0;
v[4332] = 0e0;
v[4333] = 0e0;
v[4334] = -v[1964];
v[4335] = -v[1963];
v[4336] = -v[1962];
v[4337] = 0e0;
v[4338] = 0e0;
v[4339] = 0e0;
v[4340] = 0e0;
v[4341] = 0e0;
v[4342] = 0e0;
v[4289] = 0e0;
v[4290] = 0e0;
v[4291] = 0e0;
v[4292] = 0e0;
v[4293] = 0e0;
v[4294] = 0e0;
v[4295] = 0e0;
v[4296] = 0e0;
v[4297] = 0e0;
v[4298] = 0e0;
v[4299] = 0e0;
v[4300] = 0e0;
v[4301] = -v[1964];
v[4302] = -v[1963];
v[4303] = -v[1962];
v[4304] = 0e0;
v[4305] = 0e0;
v[4306] = 0e0;
v[4253] = 0e0;
v[4254] = 0e0;
v[4255] = 0e0;
v[4256] = 0e0;
v[4257] = 0e0;
v[4258] = 0e0;
v[4259] = 0e0;
v[4260] = 0e0;
v[4261] = 0e0;
v[4262] = 0e0;
v[4263] = 0e0;
v[4264] = 0e0;
v[4265] = 0e0;
v[4266] = 0e0;
v[4267] = 0e0;
v[4268] = -v[1964];
v[4269] = -v[1963];
v[4270] = -v[1962];
v[1965] = -(v[177] * v[1962]) - v[173] * v[1963] - v[169] * v[1964] + v[1925] * v[2275] + v[192] * v[2410]
+ v[230] * v[2414] - v[190] * v[2415] - v[191] * v[2416];
v[1966] = v[149] * v[1962] + v[145] * v[1963] + v[141] * v[1964] - v[1922] * v[2275] + v[2060] * v[230] - v[183] * v[2410]
- v[182] * v[2411] - v[181] * v[2412];
v[1967] = (-(v[179] * v[1962]) - v[175] * v[1963] - v[171] * v[1964] - v[1965] + v[1927] * v[2275] + v[198] * v[2410]
	- v[196] * v[2415] - v[197] * v[2416] + v[230] * v[2417]) / 2e0;
v[1968] = (-(v[178] * v[1962]) - v[174] * v[1963] - v[170] * v[1964] - v[1965] + v[1926] * v[2275] + v[195] * v[2410]
	- v[193] * v[2415] - v[194] * v[2416] + v[230] * v[2418]) / 2e0;
v[1969] = (v[151] * v[1962] + v[147] * v[1963] + v[143] * v[1964] - v[1966] - v[1924] * v[2275] + v[2056] * v[230]
	- v[189] * v[2410] - v[188] * v[2411] - v[187] * v[2412]) / 2e0;
v[1970] = (v[150] * v[1962] + v[146] * v[1963] + v[142] * v[1964] - v[1966] - v[1923] * v[2275] + v[2058] * v[230]
	- v[186] * v[2410] - v[185] * v[2411] - v[184] * v[2412]) / 2e0;
for (i1844 = 1; i1844 <= 18; i1844++) {
	v[2030] = v[4212 + i1844];
	v[2027] = v[3624 + i1844];
	v[2025] = v[3642 + i1844];
	v[2023] = v[3660 + i1844];
	v[1978] = v[3584 + i1844];
	v[2002] = -v[1978] / 2e0;
	v[1977] = v[3566 + i1844];
	v[1976] = v[3548 + i1844];
	v[1999] = -v[1976] / 2e0;
	v[1975] = v[3530 + i1844];
	v[1979] = (i1844 == 1 ? v[10] : 0e0);
	v[1980] = (i1844 == 2 ? v[10] : 0e0);
	v[1981] = (i1844 == 4 ? v[10] : 0e0);
	v[1982] = (i1844 == 5 ? v[10] : 0e0);
	v[1983] = (i1844 == 7 ? v[10] : 0e0);
	v[1984] = (i1844 == 8 ? v[10] : 0e0);
	v[1985] = (i1844 == 10 ? v[10] : 0e0);
	v[1986] = (i1844 == 11 ? v[10] : 0e0);
	v[1987] = (i1844 == 13 ? v[10] : 0e0);
	v[1988] = (i1844 == 14 ? v[10] : 0e0);
	v[1989] = (i1844 == 16 ? v[10] : 0e0);
	v[1990] = (i1844 == 17 ? v[10] : 0e0);
	v[1991] = (i1844 == 3 ? v[10] : 0e0);
	v[1992] = (i1844 == 6 ? v[10] : 0e0);
	v[1993] = (i1844 == 9 ? v[10] : 0e0);
	v[1994] = (i1844 == 12 ? v[10] : 0e0);
	v[1995] = (i1844 == 15 ? v[10] : 0e0);
	v[1996] = (i1844 == 18 ? v[10] : 0e0);
	v[1997] = v[1975] / 2e0;
	v[1998] = -v[1997] + v[1999];
	v[2000] = v[1977] / 2e0;
	v[2001] = -v[2000] + v[2002];
	v[2440] = v[186] * v[1997] + v[183] * v[1998] - v[189] * v[1999] - v[195] * v[2000] - v[192] * v[2001] + v[198] * v[2002];
	v[2438] = v[185] * v[1997] + v[182] * v[1998] - v[188] * v[1999] - v[194] * v[2000] - v[191] * v[2001] + v[197] * v[2002];
	v[2437] = v[184] * v[1997] + v[181] * v[1998] - v[187] * v[1999] - v[193] * v[2000] - v[190] * v[2001] + v[196] * v[2002];
	v[2427] = v[1923] * v[1997] + v[1922] * v[1998] - v[1924] * v[1999] - v[1926] * v[2000] - v[1925] * v[2001]
		+ v[1927] * v[2002];
	v[2004] = v[142] * v[1997] + v[141] * v[1998] - v[143] * v[1999] - v[170] * v[2000] - v[169] * v[2001] + v[171] * v[2002]
		+ v[2027];
	v[2006] = v[146] * v[1997] + v[145] * v[1998] - v[147] * v[1999] - v[174] * v[2000] - v[173] * v[2001] + v[175] * v[2002]
		+ v[2025];
	v[2008] = v[150] * v[1997] + v[149] * v[1998] - v[151] * v[1999] - v[178] * v[2000] - v[177] * v[2001] + v[179] * v[2002]
		+ v[2023];
	v[2419] = v[199] * v[2004] + v[200] * v[2006] + v[2008] * v[201];
	v[2011] = v[2419] / v[482];
	v[2009] = -(v[1862] * v[2419] * v[636]);
	v[2010] = v[2011] * v[212];
	v[2013] = v[2004];
	v[2121] = v[2013];
	v[2014] = v[199] * v[2010] + v[2004] * v[213];
	v[2015] = v[2006];
	v[2120] = v[2015];
	v[2016] = v[200] * v[2010] + v[2006] * v[213];
	v[2017] = v[2008];
	v[2119] = v[2017];
	v[2018] = v[201] * v[2010] + v[2008] * v[213];
	b2019 = b226;
	if (b2019) {
		v[2020] = -v[2014];
		v[2021] = -v[2016];
		v[2022] = -v[2018];
	}
	else {
		v[2020] = v[2014];
		v[2021] = v[2016];
		v[2022] = v[2018];
	};
	v[2429] = -(v[2020] * v[229]);
	v[2428] = -(v[2021] * v[228]);
	v[2425] = v[125] * v[2022] + v[1998] * v[230];
	v[2424] = v[127] * v[2022] + v[1997] * v[230];
	v[2423] = v[129] * v[2022] - v[1999] * v[230];
	v[2422] = v[153] * v[2022] + v[2001] * v[230];
	v[2421] = v[155] * v[2022] + v[2000] * v[230];
	v[2420] = -(v[157] * v[2022]) + v[2002] * v[230];
	v[2024] = v[1178] * v[2022] + v[230] * (v[10] * v[2023] + v[2440]);
	v[2026] = v[1358] * v[2022] + v[230] * (v[10] * v[2025] + v[2438]);
	v[2028] = v[1343] * v[2022] + v[230] * (v[10] * v[2027] + v[2437]);
	v[2029] = v[2022] * v[2344];
	v[2434] = -(v[1851] * v[2029]);
	v[2031] = 2e0*v[1865] * v[2022] + v[10] * v[2030] - v[2002] * v[2050] + v[2000] * v[2052] + v[2001] * v[2054]
		- v[1999] * v[2056] + v[1997] * v[2058] + v[1998] * v[2060];
	v[2138] = v[2031];
	v[2051] = v[2022] * v[2050] + v[2426] * v[4234 + i1844] + v[4252 + i1844];
	v[2053] = v[2022] * v[2052] + v[2426] * v[4270 + i1844] + v[4288 + i1844];
	v[2055] = v[2022] * v[2054] + v[2426] * v[4306 + i1844] + v[4324 + i1844];
	v[2057] = v[2022] * v[2056] + v[2426] * v[4342 + i1844] + v[4360 + i1844];
	v[2059] = v[2022] * v[2058] + v[2426] * v[4378 + i1844] + v[4396 + i1844];
	v[2061] = v[2022] * v[2060] + v[2426] * v[4414 + i1844] + v[4432 + i1844];
	v[2062] = v[228] * (-(v[157] * v[2021]) + v[2002] * v[229]);
	v[2063] = -(v[228] * (v[155] * v[2021] + v[2000] * v[229]));
	v[2064] = -(v[228] * (v[153] * v[2021] + v[2001] * v[229]));
	v[2065] = v[228] * (v[129] * v[2021] - v[1999] * v[229]);
	v[2066] = v[228] * (v[127] * v[2021] + v[1997] * v[229]);
	v[2067] = v[228] * (v[125] * v[2021] + v[1998] * v[229]);
	v[2068] = v[2021] * v[2347];
	v[2435] = -(v[1850] * v[2068]);
	v[2069] = v[2021] * v[2084] + v[229] * v[2427];
	v[2070] = 2e0*v[1896] * v[2021] + v[228] * v[2427];
	v[2077] = v[2062] + v[157] * v[2429];
	v[2078] = v[2063] + v[155] * v[2429];
	v[2079] = v[2064] + v[153] * v[2429];
	v[2080] = v[2065] - v[129] * v[2429];
	v[2081] = v[2066] - v[127] * v[2429];
	v[2082] = v[2067] - v[125] * v[2429];
	v[2083] = v[2020] * v[2350];
	v[2436] = -(v[1849] * v[2083]);
	v[2069] = 2e0*v[1921] * v[2020] + v[2069];
	v[2070] = v[2070] + v[2020] * v[2084];
	v[2091] = v[1849] * v[2020];
	v[2092] = v[1850] * v[2021] + v[2091];
	v[2070] = v[1850] * v[2024] + v[1851] * v[2026] + v[2070];
	v[2140] = v[2070];
	v[2069] = v[1849] * v[2024] + v[1851] * v[2028] + v[2069];
	v[2142] = v[2069];
	v[2093] = v[1851] * v[2021];
	v[2099] = v[1851] * v[2020];
	b2105 = b483;
	if (b2105) {
		b2106 = b485;
		if (b2106) {
			v[2020] = 0e0;
			v[2021] = 0e0;
			v[2022] = 0e0;
		}
		else {
		};
	}
	else {
	};
	v[2107] = 0e0;
	v[2108] = 0e0;
	v[2109] = 0e0;
	v[2110] = 0e0;
	v[2111] = 0e0;
	v[2112] = 0e0;
	v[2113] = 0e0;
	b2114 = b483;
	if (b2114) {
		v[2115] = 0e0;
		v[2116] = 0e0;
		v[2117] = 0e0;
		b2118 = b502;
		if (b2118) {
			v[2117] = v[2017];
			v[2017] = 0e0;
			v[2116] = v[2015];
			v[2015] = 0e0;
			v[2115] = v[2013];
			v[2013] = 0e0;
		}
		else {
			v[2112] = -v[2119];
			v[2017] = 0e0;
			v[2111] = -v[2120];
			v[2015] = 0e0;
			v[2110] = -v[2121];
			v[2013] = 0e0;
		};
		v[2430] = (v[2115] * v[457] + v[2116] * v[471] + v[2117] * v[473])*v[7];
		b2122 = b485;
		if (b2122) {
			v[2109] = v[2117] * v[497];
			v[2108] = v[2116] * v[497];
			v[2107] = v[2115] * v[497];
			v[2113] = (v[2430] * v[495] * v[860] * Power(v[494], v[1442])) / v[861];
		}
		else {
			v[2109] = v[2117] * v[501];
			v[2108] = v[2116] * v[501];
			v[2107] = v[2115] * v[501];
			v[2009] = v[2009] + (v[2430] * v[500] * v[868] * Power(v[482], v[1445])) / v[869];
		};
	}
	else {
	};
	v[2441] = v[2107] * v[456];
	v[2442] = v[2108] * v[464];
	v[2443] = v[2109] * v[472];
	v[2143] = v[2009];
	v[2141] = v[2110];
	v[2139] = v[2111];
	v[2137] = v[2112];
	b2132 = b483;
	if (b2132) {
		v[2433] = v[2110] * v[228];
		v[2432] = v[2111] * v[229];
		v[2431] = v[2112] * v[230];
		b2133 = b485;
		if (b2133) {
			v[2031] = v[2031] + v[2112] * v[2282];
			v[2112] = 0e0;
			v[2070] = v[2070] + v[2111] * v[2282];
			v[2111] = 0e0;
			v[2069] = v[2069] + v[2110] * v[2282];
			v[2110] = 0e0;
			v[2113] = v[2113] + (-v[2431] - v[2432] - v[2433])*v[31] * v[34] * Power(v[494], v[495]);
			v[2009] = v[2009] - v[2113];
		}
		else {
			v[2031] = v[2138] + v[2137] * v[491];
			v[2112] = 0e0;
			v[2070] = v[2140] + v[2139] * v[491];
			v[2111] = 0e0;
			v[2069] = v[2142] + v[2141] * v[491];
			v[2110] = 0e0;
			v[2009] = v[2143] + v[2291] * (v[2431] + v[2432] + v[2433])*Power(v[482], v[500]);
		};
	}
	else {
	};
	v[2069] = v[2069] + v[1422] * v[2109];
	v[2070] = v[2070] + v[1415] * v[2109];
	v[2145] = v[2099] + v[2109] * v[228];
	v[2151] = v[2093] + v[2109] * v[229];
	v[2070] = v[2070] + v[2108] * v[956];
	v[2201] = v[1850] * v[1979] + v[1849] * v[1980] + v[182] * v[2107] + v[181] * v[2108];
	v[2202] = v[1850] * v[1981] + v[1849] * v[1982] + v[185] * v[2107] + v[184] * v[2108];
	v[2203] = v[1850] * v[1983] + v[1849] * v[1984] + v[188] * v[2107] + v[187] * v[2108];
	v[2204] = v[1850] * v[1985] + v[1849] * v[1986] + v[191] * v[2107] + v[190] * v[2108];
	v[2205] = v[1850] * v[1987] + v[1849] * v[1988] + v[194] * v[2107] + v[193] * v[2108];
	v[2206] = v[1850] * v[1989] + v[1849] * v[1990] + v[197] * v[2107] + v[196] * v[2108];
	v[2439] = v[125] * v[2201] + v[127] * v[2202] + v[129] * v[2203] - v[153] * v[2204] - v[155] * v[2205] - v[157] * v[2206];
	v[2069] = v[2069] + v[2107] * v[956];
	v[4469] = v[1850] * v[2082] + v[1866] * v[2425] - v[125] * v[2436] + (v[1849] * v[1998] + v[125] * v[2107])*v[456]
		+ v[2108] * v[458];
	v[4470] = v[1849] * v[2082] + v[1872] * v[2425] - v[125] * v[2435] + v[2107] * v[458] + (v[1850] * v[1998]
		+ v[125] * v[2108])*v[464];
	v[4471] = v[1928] * v[2425] - v[125] * v[2434] + (v[1851] * v[1998] + v[125] * v[2109])*v[472];
	v[4472] = v[1850] * v[2081] + v[1866] * v[2424] - v[127] * v[2436] + (v[1849] * v[1997] + v[127] * v[2107])*v[456]
		+ v[2108] * v[459];
	v[4473] = v[1849] * v[2081] + v[1872] * v[2424] - v[127] * v[2435] + v[2107] * v[459] + (v[1850] * v[1997]
		+ v[127] * v[2108])*v[464];
	v[4474] = v[1928] * v[2424] - v[127] * v[2434] + (v[1851] * v[1997] + v[127] * v[2109])*v[472];
	v[4475] = v[1850] * v[2080] + v[1866] * v[2423] - v[129] * v[2436] + (-(v[1849] * v[1999]) + v[129] * v[2107])*v[456]
		+ v[2108] * v[460];
	v[4476] = v[1849] * v[2080] + v[1872] * v[2423] - v[129] * v[2435] + v[2107] * v[460] + (-(v[1850] * v[1999])
		+ v[129] * v[2108])*v[464];
	v[4477] = v[1928] * v[2423] - v[129] * v[2434] + (-(v[1851] * v[1999]) + v[129] * v[2109])*v[472];
	v[4478] = v[1850] * v[2079] - v[1866] * v[2422] + v[153] * v[2436] - (v[1849] * v[2001] + v[153] * v[2107])*v[456]
		+ v[2108] * v[461];
	v[4479] = v[1849] * v[2079] - v[1872] * v[2422] + v[153] * v[2435] + v[2107] * v[461] - (v[1850] * v[2001]
		+ v[153] * v[2108])*v[464];
	v[4480] = -(v[1928] * v[2422]) + v[153] * v[2434] - (v[1851] * v[2001] + v[153] * v[2109])*v[472];
	v[4481] = v[1850] * v[2078] - v[1866] * v[2421] + v[155] * v[2436] - (v[1849] * v[2000] + v[155] * v[2107])*v[456]
		+ v[2108] * v[462];
	v[4482] = v[1849] * v[2078] - v[1872] * v[2421] + v[155] * v[2435] + v[2107] * v[462] - (v[1850] * v[2000]
		+ v[155] * v[2108])*v[464];
	v[4483] = -(v[1928] * v[2421]) + v[155] * v[2434] - (v[1851] * v[2000] + v[155] * v[2109])*v[472];
	v[4484] = v[1850] * v[2077] + v[1866] * v[2420] + v[157] * v[2436] + (v[1849] * v[2002] - v[157] * v[2107])*v[456]
		+ v[2108] * v[463];
	v[4485] = v[1849] * v[2077] + v[1872] * v[2420] + v[157] * v[2435] + v[2107] * v[463] + (v[1850] * v[2002]
		- v[157] * v[2108])*v[464];
	v[4486] = v[1928] * v[2420] + v[157] * v[2434] + (v[1851] * v[2002] - v[157] * v[2109])*v[472];
	v[2069] = v[2069] + v[2350] * (v[1343] * v[2107] + v[1849] * (v[125] * v[1979] + v[127] * v[1981] + v[129] * v[1983]
		- v[153] * v[1985] - v[155] * v[1987] - v[157] * v[1989] + v[2437])) + v[229] * v[2439];
	v[2070] = v[2070] + v[2347] * (v[1358] * v[2108] + v[1850] * (v[125] * v[1980] + v[127] * v[1982] + v[129] * v[1984]
		- v[153] * v[1986] - v[155] * v[1988] - v[157] * v[1990] + v[2438])) + v[228] * v[2439];
	v[2231] = v[2092] + v[2107] * v[228] + v[2108] * v[229];
	v[4487] = v[125] * v[2145];
	v[4488] = v[125] * v[2151];
	v[4489] = v[125] * v[2231];
	v[4490] = v[127] * v[2145];
	v[4491] = v[127] * v[2151];
	v[4492] = v[127] * v[2231];
	v[4493] = v[129] * v[2145];
	v[4494] = v[129] * v[2151];
	v[4495] = v[129] * v[2231];
	v[4496] = -(v[153] * v[2145]);
	v[4497] = -(v[153] * v[2151]);
	v[4498] = -(v[153] * v[2231]);
	v[4499] = -(v[155] * v[2145]);
	v[4500] = -(v[155] * v[2151]);
	v[4501] = -(v[155] * v[2231]);
	v[4502] = -(v[157] * v[2145]);
	v[4503] = -(v[157] * v[2151]);
	v[4504] = -(v[157] * v[2231]);
	v[2031] = v[2031] + v[1343] * v[2145] + v[1358] * v[2151] + v[1178] * v[2231] + v[2344] * (v[1178] * v[2109] + v[1851] *
		(v[125] * v[1991] + v[127] * v[1992] + v[129] * v[1993] - v[153] * v[1994] - v[155] * v[1995] - v[157] * v[1996] + v[2440]
			));
	b2232 = b226;
	if (b2232) {
		v[2233] = -v[2031];
		v[2234] = -v[2070];
		v[2235] = -v[2069];
	}
	else {
		v[2233] = v[2031];
		v[2234] = v[2070];
		v[2235] = v[2069];
	};
	v[2009] = v[2009] + v[1960] * v[2011] * v[211] + v[212] * (v[1955] * v[2004] + v[1954] * v[2006] + v[1953] * v[2008]
		+ v[201] * v[2233] + v[200] * v[2234] + v[199] * v[2235]);
	v[2244] = v[1953] * v[2010] + v[213] * v[2233] + (v[1862] * v[2008] + v[2009] * v[201]) / v[482];
	v[2246] = v[1954] * v[2010] + v[213] * v[2234] + (v[1862] * v[2006] + v[200] * v[2009]) / v[482];
	v[2248] = v[1955] * v[2010] + v[213] * v[2235] + (v[1862] * v[2004] + v[199] * v[2009]) / v[482];
	v[4451] = v[1964] * v[1998] + v[125] * v[2248];
	v[4452] = v[1963] * v[1998] + v[125] * v[2246];
	v[4453] = v[1962] * v[1998] + v[125] * v[2244];
	v[4454] = v[1964] * v[1997] + v[127] * v[2248];
	v[4455] = v[1963] * v[1997] + v[127] * v[2246];
	v[4456] = v[1962] * v[1997] + v[127] * v[2244];
	v[4457] = -(v[1964] * v[1999]) + v[129] * v[2248];
	v[4458] = -(v[1963] * v[1999]) + v[129] * v[2246];
	v[4459] = -(v[1962] * v[1999]) + v[129] * v[2244];
	v[4460] = -(v[1964] * v[2001]) - v[153] * v[2248];
	v[4461] = -(v[1963] * v[2001]) - v[153] * v[2246];
	v[4462] = -(v[1962] * v[2001]) - v[153] * v[2244];
	v[4463] = -(v[1964] * v[2000]) - v[155] * v[2248];
	v[4464] = -(v[1963] * v[2000]) - v[155] * v[2246];
	v[4465] = -(v[1962] * v[2000]) - v[155] * v[2244];
	v[4466] = v[1964] * v[2002] - v[157] * v[2248];
	v[4467] = v[1963] * v[2002] - v[157] * v[2246];
	v[4468] = v[1962] * v[2002] - v[157] * v[2244];
	v[2249] = v[2055] - v[177] * v[2244] - v[173] * v[2246] - v[169] * v[2248] + v[2204] * v[2275] + (-(v[190] * v[2145])
		- v[191] * v[2151] - v[192] * v[2231])*v[230] + v[1925] * (v[2428] + v[2429]) - v[190] * v[2441] - v[191] * v[2442]
		- v[192] * v[2443] + v[1849] * (-(v[190] * v[2083]) - v[1985] * v[456]) + v[1850] * (-(v[191] * v[2068])
			- v[1986] * v[464]) + v[1851] * (-(v[192] * v[2029]) - v[1994] * v[472]);
	v[2250] = v[2061] + v[149] * v[2244] + v[145] * v[2246] + v[141] * v[2248] - v[2201] * v[2275] + (v[181] * v[2145]
		+ v[182] * v[2151] + v[183] * v[2231])*v[230] + v[1922] * (-v[2428] - v[2429]) + v[181] * v[2441] + v[182] * v[2442]
		+ v[183] * v[2443] + v[1849] * (v[181] * v[2083] + v[1979] * v[456]) + v[1850] * (v[182] * v[2068] + v[1980] * v[464])
		+ v[1851] * (v[183] * v[2029] + v[1991] * v[472]);
	v[2251] = (v[2051] - v[179] * v[2244] - v[175] * v[2246] - v[171] * v[2248] - v[2249] + v[2206] * v[2275] + (-
		(v[196] * v[2145]) - v[197] * v[2151] - v[198] * v[2231])*v[230] + v[1927] * (v[2428] + v[2429]) - v[196] * v[2441]
		- v[197] * v[2442] - v[198] * v[2443] + v[1849] * (-(v[196] * v[2083]) - v[1989] * v[456]) + v[1850] * (-
		(v[197] * v[2068]) - v[1990] * v[464]) + v[1851] * (-(v[198] * v[2029]) - v[1996] * v[472])) / 2e0;
	v[2252] = (v[2053] - v[178] * v[2244] - v[174] * v[2246] - v[170] * v[2248] - v[2249] + v[2205] * v[2275] + (-
		(v[193] * v[2145]) - v[194] * v[2151] - v[195] * v[2231])*v[230] + v[1926] * (v[2428] + v[2429]) - v[193] * v[2441]
		- v[194] * v[2442] - v[195] * v[2443] + v[1849] * (-(v[193] * v[2083]) - v[1987] * v[456]) + v[1850] * (-
		(v[194] * v[2068]) - v[1988] * v[464]) + v[1851] * (-(v[195] * v[2029]) - v[1995] * v[472])) / 2e0;
	v[2253] = (v[2057] + v[151] * v[2244] + v[147] * v[2246] + v[143] * v[2248] - v[2250] - v[2203] * v[2275] +
		(v[187] * v[2145] + v[188] * v[2151] + v[189] * v[2231])*v[230] + v[1924] * (-v[2428] - v[2429]) + v[187] * v[2441]
		+ v[188] * v[2442] + v[189] * v[2443] + v[1849] * (v[187] * v[2083] + v[1983] * v[456]) + v[1850] * (v[188] * v[2068]
			+ v[1984] * v[464]) + v[1851] * (v[189] * v[2029] + v[1993] * v[472])) / 2e0;
	v[2254] = (v[2059] + v[150] * v[2244] + v[146] * v[2246] + v[142] * v[2248] - v[2250] - v[2202] * v[2275] +
		(v[184] * v[2145] + v[185] * v[2151] + v[186] * v[2231])*v[230] + v[1923] * (-v[2428] - v[2429]) + v[184] * v[2441]
		+ v[185] * v[2442] + v[186] * v[2443] + v[1849] * (v[184] * v[2083] + v[1981] * v[456]) + v[1850] * (v[185] * v[2068]
			+ v[1982] * v[464]) + v[1851] * (v[186] * v[2029] + v[1992] * v[472])) / 2e0;
	Rc[i1844 - 1] += v[1970] * v[1975] + v[1969] * v[1976] + v[1968] * v[1977] + v[1967] * v[1978] + v[4176 + i1844] + v[10] *
		(v[2030] * v[230] + v[4194 + i1844]);
	for (i1973 = 1; i1973 <= 18; i1973++) {
		Kc[i1844 - 1][i1973 - 1] += v[2254] * v[3530 + i1973] + v[2253] * v[3548 + i1973] + v[2252] * v[3566 + i1973]
			+ v[2251] * v[3584 + i1973] + v[4450 + i1973] + v[10] * (v[4468 + i1973] + v[230] * v[4486 + i1973]);
	};/* end for */
};/* end for */
#pragma endregion

	delete[] v;
}

void FlexibleTriangularFace_FlexibleTriangularFace::HessianPhase1(Matrix& mHes)
{
	double Hes[4][4];
	v = DBG_NEW double[500];


#pragma region AceGen
	int i01; int i02;
	v[83] = (-u1A[0] - xi1A[0]) / 2e0;
	v[89] = (-u1A[1] - xi1A[1]) / 2e0;
	v[95] = (-u1A[2] - xi1A[2]) / 2e0;
	v[82] = v[83] + (u2A[0] + xi2A[0]) / 2e0;
	v[88] = v[89] + (u2A[1] + xi2A[1]) / 2e0;
	v[94] = v[95] + (u2A[2] + xi2A[2]) / 2e0;
	v[85] = v[83] + (u3A[0] + xi3A[0]) / 2e0;
	v[135] = 2e0*v[85];
	v[91] = v[89] + (u3A[1] + xi3A[1]) / 2e0;
	v[136] = 2e0*v[91];
	v[97] = v[95] + (u3A[2] + xi3A[2]) / 2e0;
	v[137] = 2e0*v[97];
	v[101] = (-u1B[0] - xi1B[0]) / 2e0;
	v[107] = (-u1B[1] - xi1B[1]) / 2e0;
	v[113] = (-u1B[2] - xi1B[2]) / 2e0;
	v[100] = v[101] + (u2B[0] + xi2B[0]) / 2e0;
	v[138] = 2e0*v[100];
	v[106] = v[107] + (u2B[1] + xi2B[1]) / 2e0;
	v[139] = 2e0*v[106];
	v[112] = v[113] + (u2B[2] + xi2B[2]) / 2e0;
	v[140] = 2e0*v[112];
	v[103] = v[101] + (u3B[0] + xi3B[0]) / 2e0;
	v[109] = v[107] + (u3B[1] + xi3B[1]) / 2e0;
	v[115] = v[113] + (u3B[2] + xi3B[2]) / 2e0;
	Hes[0][0] = 1e0*((v[82] * v[82]) + (v[88] * v[88]) + (v[94] * v[94]));
	Hes[0][1] = 0.5e0*(v[135] * v[82] + v[136] * v[88] + v[137] * v[94]);
	Hes[0][2] = 0.5e0*(-(v[138] * v[82]) - v[139] * v[88] - v[140] * v[94]);
	Hes[0][3] = -1e0*(v[103] * v[82] + v[109] * v[88] + v[115] * v[94]);
	Hes[1][1] = 1e0*((v[85] * v[85]) + (v[91] * v[91]) + (v[97] * v[97]));
	Hes[1][2] = -1e0*(v[100] * v[85] + v[106] * v[91] + v[112] * v[97]);
	Hes[1][3] = 0.5e0*(-(v[103] * v[135]) - v[109] * v[136] - v[115] * v[137]);
	Hes[2][2] = 1e0*((v[100] * v[100]) + (v[106] * v[106]) + (v[112] * v[112]));
	Hes[2][3] = 0.5e0*(v[103] * v[138] + v[109] * v[139] + v[115] * v[140]);
	Hes[3][3] = 1e0*((v[103] * v[103]) + (v[109] * v[109]) + (v[115] * v[115]));
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

void FlexibleTriangularFace_FlexibleTriangularFace::Report()
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
		}
	}
}

void FlexibleTriangularFace_FlexibleTriangularFace::CompactReport()
{
	db.myprintf("Eligible %d\n", (int)eligible);
	db.myprintf("Deg point A: %d\n", deg_pointA);
	db.myprintf("Deg point B: %d\n", deg_pointB);
	db.myprintf("Deg curve A: %d\n", deg_curveA);
	db.myprintf("Deg curve B: %d\n", deg_curveB);
	db.myprintf("Gap %.6e\n", cd->g_n[0]);
}

void FlexibleTriangularFace_FlexibleTriangularFace::PredictorTimeStep(double kin)
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
		//Gamma B
		double zb = cd->convective[0][2];
		double thb = cd->convective[0][3];
		double N1B = -(thb / 2.0) - zb / 2.0;
		double N2B = 0.5 + zb / 2.0;
		double N3B = 0.5 + thb / 2.0;

		Matrix VP_A(3);
		Matrix VP_B(3);
		for (int i = 0; i < 3; i++)
		{
			VP_A(i, 0) = N1A * dui1A[i] + N2A * dui2A[i] + N3A * dui3A[i];
			VP_B(i, 0) = N1B * dui1B[i] + N2B * dui2B[i] + N3B * dui3B[i];
		}
		
		double vrel = dot(VP_A - VP_B, *cd->n[0]);
		double deltaS = (cd->g_n[0] - (*gnb));
		if (vrel < 0)
			td->time_step_impact = (deltaS) / (-vrel);
		else
			td->time_step_impact = db.solution[db.current_solution_number - 1]->end_time;	//valor alto
	}
}
void FlexibleTriangularFace_FlexibleTriangularFace::AllocSpecificExplicit()
{
	//TODO-Explicit
}
void FlexibleTriangularFace_FlexibleTriangularFace::FreeSpecificExplicit()
{
	//TODO-Explicit
}
void FlexibleTriangularFace_FlexibleTriangularFace::MountLocalContributionsExplicit(double t)
{
	//TODO-Explicit
}
void FlexibleTriangularFace_FlexibleTriangularFace::SetVariablesExplicit(double t)
{
	//TODO-Explicit
}
void FlexibleTriangularFace_FlexibleTriangularFace::FinalUpdateExplicit(double t)
{
	//TODO-Explicit
}