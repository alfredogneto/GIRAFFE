#include "ContactVEMPolyhedronVEMPolyhedron.h"

#include "VEMPolyhedron.h"
#include "STLSurface.h"
#include "TriangularFace.h"
#include "Matrix.h"
#include "FlexibleTriangularFace_FlexibleTriangularFace.h"
#include "SurfacePairGeneralContact.h"
#include "GeneralContactSearch.h"
#include "SSContactData.h"
#include "SuperNode.h"

#include "Database.h"
//Variaveis globais
extern
Database db;


ContactVEMPolyhedronVEMPolyhedron::ContactVEMPolyhedronVEMPolyhedron()
{
	index1 = 0;				//Particle 1 - index
	index2 = 0;				//Particle 2 - index
	sub_index1 = 0;			//Particle 1 - sub_index
	sub_index2 = 0;			//Particle 2 - sub_index

	prev_active = false;
	cur_active = false;

	I3 = new Matrix(3, 3);
	(*I3)(0, 0) = 1.0;
	(*I3)(1, 1) = 1.0;
	(*I3)(2, 2) = 1.0;

	contact_pairs.clear();

	//Degenerations
	deg_pointA = 0;
	deg_curveA = 0;
	deg_pointB = 0;
	deg_curveB = 0;
}


ContactVEMPolyhedronVEMPolyhedron::~ContactVEMPolyhedronVEMPolyhedron()
{
	for (int i = 0; i < contact_pairs.size(); i++)
		delete contact_pairs[i];
	contact_pairs.clear();

	//Cleaning
	delete I3;
}


void ContactVEMPolyhedronVEMPolyhedron::PreCalc()
{
	pA = static_cast<VEMPolyhedron*>(db.particles[index1]);
	pB = static_cast<VEMPolyhedron*>(db.particles[index2]);
	surfA = static_cast<STLSurface*>(db.cad_data[pA->CADDATA_ID - 1]);
	surfB = static_cast<STLSurface*>(db.cad_data[pB->CADDATA_ID - 1]);

	//Determining the faces to be used - according to the sequence of geometric entities in the polyhedron:
	//1 - faces
	//2 - edges
	//3 - vertices
	//Sub bounding volumes (faces+edges+vertices)
	int n_facesA = surfA->n_faces;
	int n_edgesA = (int)surfA->edges.size();
	int n_verticesA = (int)surfA->vertices.size();

	if (sub_index1 < n_facesA)
	{
		//No degeneration on face A
		deg_pointA = 0;
		deg_curveA = 0;
		//Face associated
		faceA = surfA->faces[sub_index1];
	}
	else
	{
		if (sub_index1 < (n_facesA + n_edgesA))
		{
			//Degeneration on face A - edge on face A is used in contact model
			deg_pointA = 0;
			deg_curveA = sub_index1 - n_facesA + 1;		//ID starts with 1
			//Face associated
			faceA = surfA->faces[surfA->edges[deg_curveA - 1].faceIDs[0] - 1];
		}
		else
		{
			//Degeneration on face A - point on face A is used in contact model
			deg_pointA = sub_index1 - n_facesA - n_edgesA + 1;		//ID starts with 1
			deg_curveA = 0;
			//Face associated
			faceA = surfA->faces[surfA->vertices[deg_pointA - 1].faceIDs[0] - 1];
		}
	}

	//Sub bounding volumes (faces+edges+vertices)
	int n_facesB = surfB->n_faces;
	int n_edgesB = (int)surfB->edges.size();
	int n_verticesB = (int)surfB->vertices.size();
	if (sub_index2 < n_facesB)
	{
		//No degeneration on face B
		deg_pointB = 0;
		deg_curveB = 0;
		//Face associated
		faceB = surfB->faces[sub_index2];
	}
	else
	{
		if (sub_index2 < (n_facesB + n_edgesB))
		{
			//Degeneration on face B - edge on face B is used in contact model
			deg_pointB = 0;
			deg_curveB = sub_index2 - n_facesB + 1;		//ID starts with 1
			//Face associated
			faceB = surfB->faces[surfB->edges[deg_curveB - 1].faceIDs[0] - 1];
		}
		else
		{
			//Degeneration on face B - point on face B is used in contact model
			deg_pointB = sub_index2 - n_facesB - n_edgesB + 1;		//ID starts with 1
			deg_curveB = 0;
			//Face associated
			faceB = surfB->faces[surfB->vertices[deg_pointB - 1].faceIDs[0] - 1];
		}
	}

	//pointers to evaluate DOFs of the particle
	vertices_i_A = pA->vertices_i;
	vertices_i_B = pB->vertices_i;
	vertices_p_A = pA->vertices_p;
	vertices_p_B = pB->vertices_p;

	////Tests to exit with no creation of a contact pair: edge-edge
	//if (deg_curveA != 0 && surfA->edges[deg_curveA - 1].concave_indicator == 0)
	//	return;
	//if (deg_curveB != 0 && surfB->edges[deg_curveB - 1].concave_indicator == 0)
	//	return;
	//Not covered cases:
	if (
		(deg_pointA == 0 && deg_curveA == 0 && deg_pointB == 0 && deg_curveB == 0) ||
		(deg_pointA == 0 && deg_curveA != 0 && deg_pointB == 0 && deg_curveB == 0) ||
		(deg_pointA == 0 && deg_curveA == 0 && deg_pointB == 0 && deg_curveB != 0)
		)
	{
		return;
	}

	CreateSurfacePair(deg_pointA, deg_curveA, deg_pointB, deg_curveB, faceA->ID, faceB->ID);
}

void ContactVEMPolyhedronVEMPolyhedron::InsertNewContact(int deg_pointA, int deg_curveA, int deg_pointB, int deg_curveB, int faceAID, int faceBID)
{
	FlexibleTriangularFace_FlexibleTriangularFace* temp = new FlexibleTriangularFace_FlexibleTriangularFace();
	temp->index1 = index1;
	temp->index2 = index2;
	temp->sub_index1 = sub_index1;
	temp->sub_index2 = sub_index2;

	temp->deg_pointA = deg_pointA;
	temp->deg_curveA = deg_curveA;
	temp->deg_pointB = deg_pointB;
	temp->deg_curveB = deg_curveB;
	temp->faceAID = faceAID;
	temp->faceBID = faceBID;
	//Assigning pointers - not creating copies, just pointing
	temp->vertices_i_A = vertices_i_A;
	temp->vertices_i_B = vertices_i_B;
	temp->vertices_p_A = vertices_p_A;
	temp->vertices_p_B = vertices_p_B;
	temp->super_node_A = pA->super_node;
	temp->super_node_B = pB->super_node;
	temp->material_A = pA->material;
	temp->material_B = pB->material;
	temp->CAD_AID = pA->CADDATA_ID;
	temp->CAD_BID = pB->CADDATA_ID;

	temp->SetActive();
	contact_pairs.push_back(temp);
}

int ContactVEMPolyhedronVEMPolyhedron::CreateSurfacePair(int deg_pointA, int deg_curveA, int deg_pointB, int deg_curveB, int faceAID, int faceBID)
{
	int list_index = -1;
	for (int i = 0; i < contact_pairs.size(); i++)
	{
		FlexibleTriangularFace_FlexibleTriangularFace* ptr = static_cast<FlexibleTriangularFace_FlexibleTriangularFace*>(contact_pairs[i]);
		if (ptr->deg_pointA == deg_pointA && ptr->deg_pointB == deg_pointB &&
			ptr->deg_curveA == deg_curveA && ptr->deg_curveB == deg_curveB && ptr->faceAID == faceAID && ptr->faceBID == faceBID)
		{
			contact_pairs[i]->SetActive();
			list_index = i;
		}
	}
	if (list_index == -1)
	{
		InsertNewContact(deg_pointA, deg_curveA, deg_pointB, deg_curveB, faceAID, faceBID);
		list_index = (int)contact_pairs.size() - 1;
	}
	return list_index;
}
void ContactVEMPolyhedronVEMPolyhedron::ProcessSurfacePair(int list_index)
{
	contact_pairs[list_index]->eligible = false;
	contact_pairs[list_index]->SolveLCP();
	contact_pairs[list_index]->EvaluateNormalGap();

}
void ContactVEMPolyhedronVEMPolyhedron::ProcessSurfacePairs()
{
	for (int i = 0; i < (int)contact_pairs.size(); i++)
		ProcessSurfacePair(i);
}

void ContactVEMPolyhedronVEMPolyhedron::FinalProcessSurfacePairExplicit(int list_index, double t)
{
	//Copies
	bool prev_elig = contact_pairs[list_index]->eligible;
	double prev_conv[4];
	for (int i = 0; i < 4; i++)
		prev_conv[i] = contact_pairs[list_index]->cd->convective[0][i];
	int prev_ret_value = contact_pairs[list_index]->cd->return_value[0];
	double prev_g_n = contact_pairs[list_index]->cd->g_n[0];

	contact_pairs[list_index]->SolveLCP();
	contact_pairs[list_index]->EvaluateNormalGap();//here the normal of contact is updated, which ensures a successful test for penetration in HaveErrors function

	//Returning values modified in SolveLCP and EvaluateNormalGap routines (except gap and normal)
	contact_pairs[list_index]->eligible = prev_elig;
	for (int i = 0; i < 4; i++)
		contact_pairs[list_index]->cd->convective[0][i] = prev_conv[i];
	contact_pairs[list_index]->cd->return_value[0] = prev_ret_value;
	contact_pairs[list_index]->cd->g_n[0] = prev_g_n;

	if (contact_pairs[list_index]->eligible)
		contact_pairs[list_index]->FinalUpdateExplicit(t);
}
void ContactVEMPolyhedronVEMPolyhedron::FinalProcessSurfacePairsExplicit(double t)
{
	for (int i = 0; i < (int)contact_pairs.size(); i++)
		FinalProcessSurfacePairExplicit(i,t);
}

void ContactVEMPolyhedronVEMPolyhedron::MountGlobalExplicit()
{
	//TODO-Explicit
}

void ContactVEMPolyhedronVEMPolyhedron::MountGlobal()
{
	//Variaveis temporarias para salvar a indexação global dos graus de liberdade a serem setados na matriz de rigidez global
	int GL_global_1 = 0;
	int GL_global_2 = 0;
	double anterior = 0;
	for (int c = 0; c < contact_pairs.size(); c++)
	{
		FlexibleTriangularFace_FlexibleTriangularFace* ptr = static_cast<FlexibleTriangularFace_FlexibleTriangularFace*>(contact_pairs[c]);
		if (contact_pairs[c]->GetActive())
		{
			if (contact_pairs[c]->eligible)
			{
				for (int i = 0; i < 18; i++)
				{
					if (i < 3)
						GL_global_1 = db.super_nodes[db.particles[index1]->super_node - 1]->DOFs[(ptr->vertexIDsA[0] - 1) * 3 + (i)];
					if (i >= 3 && i < 6)
						GL_global_1 = db.super_nodes[db.particles[index1]->super_node - 1]->DOFs[(ptr->vertexIDsA[1] - 1) * 3 + (i - 3)];
					if (i >= 6 && i < 9)
						GL_global_1 = db.super_nodes[db.particles[index1]->super_node - 1]->DOFs[(ptr->vertexIDsA[2] - 1) * 3 + (i - 6)];
					if (i >= 9 && i < 12)
						GL_global_1 = db.super_nodes[db.particles[index2]->super_node - 1]->DOFs[(ptr->vertexIDsB[0] - 1) * 3 + (i - 9)];
					if (i >= 12 && i < 15)
						GL_global_1 = db.super_nodes[db.particles[index2]->super_node - 1]->DOFs[(ptr->vertexIDsB[1] - 1) * 3 + (i - 12)];
					if (i >= 15)
						GL_global_1 = db.super_nodes[db.particles[index2]->super_node - 1]->DOFs[(ptr->vertexIDsB[2] - 1) * 3 + (i - 15)];

					//Caso o grau de liberdade seja livre:
					if (GL_global_1 > 0)
					{
						anterior = db.global_P_A(GL_global_1 - 1, 0);
						db.global_P_A(GL_global_1 - 1, 0) = anterior + ptr->Rc[i];
						anterior = db.global_I_A(GL_global_1 - 1, 0);
						db.global_I_A(GL_global_1 - 1, 0) = anterior + ptr->Rc[i];
						////Convergence criteria
						//anterior = db.global_ABS_P_A(GL_global_1 - 1, 0);
						//db.global_ABS_P_A(GL_global_1 - 1, 0) = anterior + abs(ptr->Rc[i]);
						//db.global_COUNT_ABS_P_A(GL_global_1 - 1, 0) = db.global_COUNT_ABS_P_A(GL_global_1 - 1, 0) + 1.0;
					}
					if (GL_global_1 < 0)
					{
						anterior = db.global_P_B(-GL_global_1 - 1, 0);
						db.global_P_B(-GL_global_1 - 1, 0) = anterior + ptr->Rc[i];
						////Convergence criteria
						//anterior = db.global_ABS_P_B(-GL_global_1 - 1, 0);
						//db.global_ABS_P_B(-GL_global_1 - 1, 0) = anterior + abs(ptr->Rc[i]);
						//db.global_COUNT_ABS_P_B(-GL_global_1 - 1, 0) = db.global_COUNT_ABS_P_B(-GL_global_1 - 1, 0) + 1.0;
					}
					for (int j = 0; j < 18; j++)
					{
						if (j < 3)
							GL_global_2 = db.super_nodes[db.particles[index1]->super_node - 1]->DOFs[(ptr->vertexIDsA[0] - 1) * 3 + (j)];
						if (j >= 3 && j < 6)
							GL_global_2 = db.super_nodes[db.particles[index1]->super_node - 1]->DOFs[(ptr->vertexIDsA[1] - 1) * 3 + (j - 3)];
						if (j >= 6 && j < 9)
							GL_global_2 = db.super_nodes[db.particles[index1]->super_node - 1]->DOFs[(ptr->vertexIDsA[2] - 1) * 3 + (j - 6)];
						if (j >= 9 && j < 12)
							GL_global_2 = db.super_nodes[db.particles[index2]->super_node - 1]->DOFs[(ptr->vertexIDsB[0] - 1) * 3 + (j - 9)];
						if (j >= 12 && j < 15)
							GL_global_2 = db.super_nodes[db.particles[index2]->super_node - 1]->DOFs[(ptr->vertexIDsB[1] - 1) * 3 + (j - 12)];
						if (j >= 15)
							GL_global_2 = db.super_nodes[db.particles[index2]->super_node - 1]->DOFs[(ptr->vertexIDsB[2] - 1) * 3 + (j - 15)];

						//Caso os graus de liberdade sejam ambos livres (Matriz Kaa)
						if (GL_global_1 > 0 && GL_global_2 > 0)
							db.global_stiffness_AA.setValue(GL_global_1 - 1, GL_global_2 - 1, ptr->Kc[i][j]);
						//Caso os graus de liberdade sejam ambos fixos (Matriz Kbb)
						if (GL_global_1 < 0 && GL_global_2 < 0)
							db.global_stiffness_BB.setValue(-GL_global_1 - 1, -GL_global_2 - 1, ptr->Kc[i][j]);
						//Caso os graus de liberdade sejam livre e fixo (Matriz Kab)
						if (GL_global_1 > 0 && GL_global_2 < 0)
							db.global_stiffness_AB.setValue(GL_global_1 - 1, -GL_global_2 - 1, ptr->Kc[i][j]);
						//Caso os graus de liberdade sejam fixo e livre (Matriz Kba)
						if (GL_global_1 < 0 && GL_global_2 > 0)
							db.global_stiffness_BA.setValue(-GL_global_1 - 1, GL_global_2 - 1, ptr->Kc[i][j]);
					}
				}
			}
		}
	}
}

void ContactVEMPolyhedronVEMPolyhedron::ProcessContactHierarchy()
{
	for (int c = 0; c < (int)contact_pairs.size(); c++)
	{
		//If the pair is eligible, checks its hierarchy
		if (cur_active && contact_pairs[c]->eligible)
		{
			//Indexes:	particle 1: index1
			//			particle 2: index2
			//Pointers
			VEMPolyhedron* pA = static_cast<VEMPolyhedron*>(db.particles[index1]);
			VEMPolyhedron* pB = static_cast<VEMPolyhedron*>(db.particles[index2]);
			STLSurface* surfA;
			STLSurface* surfB;
			//Surfaces A and B (pointers)
			surfA = static_cast<STLSurface*>(db.cad_data[pA->CADDATA_ID - 1]);
			surfB = static_cast<STLSurface*>(db.cad_data[pB->CADDATA_ID - 1]);

			//Edge in A and point in B (vertex-edge contact) - check hierarchy
			if (deg_curveA != 0 && deg_pointB != 0)
			{
				//Faces associated with the edge
				int face1ID;
				int face2ID;
				if (surfA->edges[deg_curveA - 1].faceIDs.size() == 1)
				{
					face1ID = surfA->edges[deg_curveA - 1].faceIDs[0];
					face2ID = -1; //single face associated
				}
				else
				{
					face1ID = surfA->edges[deg_curveA - 1].faceIDs[0];
					face2ID = surfA->edges[deg_curveA - 1].faceIDs[1];
				}
				//Convex edge in A
				if (surfA->edges[deg_curveA - 1].concave_indicator <= 0)
				{
					/////////////////////////////////////////////////////////////////////////////
					//Searching for contacts associated with face1ID and face2ID - subcycling
					bool checked = false;
					for (int subcont = 0; subcont < db.gcs->contactPP_list[index1].size() && checked == false; subcont++)
					{
						//Check if are the same particles, but not the same pair 
						if (db.gcs->contactPP_list[index1][subcont]->index2 == index2 &&
							(db.gcs->contactPP_list[index1][subcont]->sub_index1 != sub_index1 ||
								db.gcs->contactPP_list[index1][subcont]->sub_index2 != sub_index2))
						{
							for (int subc = 0; subc < (int)db.gcs->contactPP_list[index1][subcont]->contact_pairs.size() && checked == false; subc++)
							{
								//If the pair is active and eligible
								if (db.gcs->contactPP_list[index1][subcont]->cur_active && db.gcs->contactPP_list[index1][subcont]->contact_pairs[subc]->eligible)
								{
									ContactVEMPolyhedronVEMPolyhedron* subptr = static_cast<ContactVEMPolyhedronVEMPolyhedron*>(db.gcs->contactPP_list[index1][subcont]);
									//Test: subptr->faceA->ID == face1ID  & face-point contact
									if (subptr->deg_curveA == 0 && subptr->deg_pointA == 0 && subptr->faceA->ID == face1ID && subptr->deg_pointB == deg_pointB)
									{
										contact_pairs[c]->eligible = false;
										checked = true;
									}
									//Test: subptr->faceA->ID == face2ID  & face-point contact
									if (subptr->deg_curveA == 0 && subptr->deg_pointA == 0 && subptr->faceA->ID == face2ID && subptr->deg_pointB == deg_pointB)
									{
										contact_pairs[c]->eligible = false;
										checked = true;
									}
								}
							}
						}
					}
					///////////////////////////////////////////////////////////////////////////////
				}
				//Concave edge in A
				else
				{
					/////////////////////////////////////////////////////////////////////////////
					//Searching for contacts associated with face1ID and face2ID - subcycling
					bool checked = false;
					bool face1found = false;
					bool face2found = false;
					for (int subcont = 0; subcont < db.gcs->contactPP_list[index1].size() && checked == false; subcont++)
					{
						//Check if are the same particles, but not the same pair (count != subcount)

						//Check if are the same particles, but not the same pair 
						if (db.gcs->contactPP_list[index1][subcont]->index2 == index2 &&
							(db.gcs->contactPP_list[index1][subcont]->sub_index1 != sub_index1 ||
								db.gcs->contactPP_list[index1][subcont]->sub_index2 != sub_index2))
						{
							for (int subc = 0; subc < (int)db.gcs->contactPP_list[index1][subcont]->contact_pairs.size() && checked == false; subc++)
							{
								//If the pair is active and eligible
								if (db.gcs->contactPP_list[index1][subcont]->cur_active && db.gcs->contactPP_list[index1][subcont]->contact_pairs[subc]->eligible)
								{
									ContactVEMPolyhedronVEMPolyhedron* subptr = static_cast<ContactVEMPolyhedronVEMPolyhedron*>(db.gcs->contactPP_list[index1][subcont]);
									//Test: subptr->faceA->ID == face1ID  & face-point contact
									if (subptr->deg_curveA == 0 && subptr->deg_pointA == 0 && subptr->faceA->ID == face1ID && subptr->deg_pointB == deg_pointB)
										face1found = true;
									//Test: subptr->faceA->ID == face2ID  & face-point contact
									if (subptr->deg_curveA == 0 && subptr->deg_pointA == 0 && subptr->faceA->ID == face2ID && subptr->deg_pointB == deg_pointB)
										face2found = true;
									if (face1found && face2found)
									{
										contact_pairs[c]->eligible = false;
										checked = true;
									}
								}
							}
						}
					}
					///////////////////////////////////////////////////////////////////////////////
				}

				//Test related to edge-edge contact - for concave edges
				/////////////////////////////////////////////////////////////////////////////
				//Searching - subcycling
				bool checked = false;
				int cc_ed = 0;		//counter
				for (int subcont = 0; subcont < db.gcs->contactPP_list[index1].size() && checked == false; subcont++)
				{
					//Check if are the same particles, but not the same pair 
					if (db.gcs->contactPP_list[index1][subcont]->index2 == index2 &&
						(db.gcs->contactPP_list[index1][subcont]->sub_index1 != sub_index1 ||
							db.gcs->contactPP_list[index1][subcont]->sub_index2 != sub_index2))
					{
						for (int subc = 0; subc < (int)db.gcs->contactPP_list[index1][subcont]->contact_pairs.size() && checked == false; subc++)
						{
							//If the pair is active and eligible
							if (db.gcs->contactPP_list[index1][subcont]->cur_active && db.gcs->contactPP_list[index1][subcont]->contact_pairs[subc]->eligible)
							{
								ContactVEMPolyhedronVEMPolyhedron* subptr = static_cast<ContactVEMPolyhedronVEMPolyhedron*>(db.gcs->contactPP_list[index1][subcont]);

								//Edges associated with vertex B
								for (int edge_count = 0; edge_count < (int)surfB->vertices[deg_pointB - 1].edgeIDs.size(); edge_count++)
								{
									if (subptr->deg_curveA == deg_curveA && subptr->deg_curveB == surfB->vertices[deg_pointB - 1].edgeIDs[edge_count])
										cc_ed++;
								}
								if (cc_ed >= 2)
								{
									contact_pairs[c]->eligible = false;
									checked = true;
								}
							}
						}
					}
				}
				/////////////////////////////////////////////////////////////////////////////
			}
			//Edge in B and point in A (vertex-edge contact) - check hierarchy
			if (deg_curveB != 0 && deg_pointA != 0)
			{
				//Faces associated with the edgebool single_face;
				int face1ID;
				int face2ID;
				if (surfB->edges[deg_curveB - 1].faceIDs.size() == 1)
				{
					face1ID = surfB->edges[deg_curveB - 1].faceIDs[0];
					face2ID = -1;
				}
				else
				{
					face1ID = surfB->edges[deg_curveB - 1].faceIDs[0];
					face2ID = surfB->edges[deg_curveB - 1].faceIDs[1];
				}
				//Convex edge in B
				if (surfB->edges[deg_curveB - 1].concave_indicator <= 0)
				{
					/////////////////////////////////////////////////////////////////////////////
					//Searching for contacts associated with face1ID and face2ID - subcycling
					bool checked = false;
					for (int subcont = 0; subcont < db.gcs->contactPP_list[index1].size() && checked == false; subcont++)
					{
						//Check if are the same particles, but not the same pair 
						if (db.gcs->contactPP_list[index1][subcont]->index2 == index2 &&
							(db.gcs->contactPP_list[index1][subcont]->sub_index1 != sub_index1 ||
								db.gcs->contactPP_list[index1][subcont]->sub_index2 != sub_index2))
						{
							for (int subc = 0; subc < (int)db.gcs->contactPP_list[index1][subcont]->contact_pairs.size() && checked == false; subc++)
							{
								//If the pair is active and eligible
								if (db.gcs->contactPP_list[index1][subcont]->cur_active && db.gcs->contactPP_list[index1][subcont]->contact_pairs[subc]->eligible)
								{
									ContactVEMPolyhedronVEMPolyhedron* subptr = static_cast<ContactVEMPolyhedronVEMPolyhedron*>(db.gcs->contactPP_list[index1][subcont]);
									//Test: subptr->faceB->ID == face1ID  & face-point contact
									if (subptr->deg_curveB == 0 && subptr->deg_pointB == 0 && subptr->faceB->ID == face1ID && subptr->deg_pointA == deg_pointA)
									{
										contact_pairs[c]->eligible = false;
										checked = true;
									}
									//Test: subptr->faceB->ID == face2ID  & face-point contact
									if (subptr->deg_curveB == 0 && subptr->deg_pointB == 0 && subptr->faceB->ID == face2ID && subptr->deg_pointA == deg_pointA)
									{
										contact_pairs[c]->eligible = false;
										checked = true;
									}
								}
							}
						}
					}
					///////////////////////////////////////////////////////////////////////////////
				}

				//Concave edge in B
				else
				{
					/////////////////////////////////////////////////////////////////////////////
					//Searching for contacts associated with face1ID and face2ID - subcycling
					bool checked = false;
					bool face1found = false;
					bool face2found = false;
					for (int subcont = 0; subcont < db.gcs->contactPP_list[index1].size() && checked == false; subcont++)
					{
						//Check if are the same particles, but not the same pair 
						if (db.gcs->contactPP_list[index1][subcont]->index2 == index2 &&
							(db.gcs->contactPP_list[index1][subcont]->sub_index1 != sub_index1 ||
								db.gcs->contactPP_list[index1][subcont]->sub_index2 != sub_index2))
						{
							for (int subc = 0; subc < (int)db.gcs->contactPP_list[index1][subcont]->contact_pairs.size() && checked == false; subc++)
							{
								//If the pair is active and eligible
								if (db.gcs->contactPP_list[index1][subcont]->cur_active && db.gcs->contactPP_list[index1][subcont]->contact_pairs[subc]->eligible)
								{
									ContactVEMPolyhedronVEMPolyhedron* subptr = static_cast<ContactVEMPolyhedronVEMPolyhedron*>(db.gcs->contactPP_list[index1][subcont]);
									//Test: subptr->faceB->ID == face1ID  & face-point contact
									if (subptr->deg_curveB == 0 && subptr->deg_pointB == 0 && subptr->faceB->ID == face1ID && subptr->deg_pointA == deg_pointA)
										face1found = true;
									//Test: subptr->faceB->ID == face2ID  & face-point contact
									if (subptr->deg_curveB == 0 && subptr->deg_pointB == 0 && subptr->faceB->ID == face2ID && subptr->deg_pointA == deg_pointA)
										face2found = true;
									if (face1found && face2found)
									{
										contact_pairs[c]->eligible = false;
										checked = true;
									}
								}
							}
						}
					}
					///////////////////////////////////////////////////////////////////////////////
				}

				//Test related to edge-edge contact - for concave edges
				/////////////////////////////////////////////////////////////////////////////
				//Searching - subcycling
				bool checked = false;
				int cc_ed = 0;		//counter
				for (int subcont = 0; subcont < db.gcs->contactPP_list[index1].size() && checked == false; subcont++)
				{
					//Check if are the same particles, but not the same pair 
					if (db.gcs->contactPP_list[index1][subcont]->index2 == index2 &&
						(db.gcs->contactPP_list[index1][subcont]->sub_index1 != sub_index1 ||
							db.gcs->contactPP_list[index1][subcont]->sub_index2 != sub_index2))
					{
						for (int subc = 0; subc < (int)db.gcs->contactPP_list[index1][subcont]->contact_pairs.size() && checked == false; subc++)
						{
							//If the pair is active and eligible
							if (db.gcs->contactPP_list[index1][subcont]->cur_active && db.gcs->contactPP_list[index1][subcont]->contact_pairs[subc]->eligible)
							{
								ContactVEMPolyhedronVEMPolyhedron* subptr = static_cast<ContactVEMPolyhedronVEMPolyhedron*>(db.gcs->contactPP_list[index1][subcont]);

								//Edges associated with vertex A
								for (int edge_count = 0; edge_count < (int)surfA->vertices[deg_pointA - 1].edgeIDs.size(); edge_count++)
								{
									if (subptr->deg_curveB == deg_curveB && subptr->deg_curveA == surfA->vertices[deg_pointA - 1].edgeIDs[edge_count])
										cc_ed++;
								}
								if (cc_ed >= 2)
								{
									contact_pairs[c]->eligible = false;
									checked = true;
								}
							}
						}
					}
				}
				/////////////////////////////////////////////////////////////////////////////
			}

			//Vertex in A and vertex in B (vertex-vertex contact) - check hierarchy
			if (deg_pointA != 0 && deg_pointB != 0)
			{
				/////////////////////////////////////////////////////////////////////////////
				//Searching - subcycling
				bool checked = false;
				for (int subcont = 0; subcont < db.gcs->contactPP_list[index1].size() && checked == false; subcont++)
				{
					//Check if are the same particles, but not the same pair 
					if (db.gcs->contactPP_list[index1][subcont]->index2 == index2 &&
						(db.gcs->contactPP_list[index1][subcont]->sub_index1 != sub_index1 ||
							db.gcs->contactPP_list[index1][subcont]->sub_index2 != sub_index2))
					{
						for (int subc = 0; subc < (int)db.gcs->contactPP_list[index1][subcont]->contact_pairs.size() && checked == false; subc++)
						{
							//If the pair is active and eligible
							if (db.gcs->contactPP_list[index1][subcont]->cur_active && db.gcs->contactPP_list[index1][subcont]->contact_pairs[subc]->eligible)
							{
								ContactVEMPolyhedronVEMPolyhedron* subptr = static_cast<ContactVEMPolyhedronVEMPolyhedron*>(db.gcs->contactPP_list[index1][subcont]);

								//Edges associated with vertex A
								//surfA->vertices[ptr->deg_pointA - 1].edgeIDs.size();
								for (int edge_count = 0; edge_count < (int)surfA->vertices[deg_pointA - 1].edgeIDs.size(); edge_count++)
								{
									if (subptr->deg_curveA == surfA->vertices[deg_pointA - 1].edgeIDs[edge_count]
										&& subptr->deg_pointB == deg_pointB)
									{
										contact_pairs[c]->eligible = false;
										checked = true;
									}
								}
								//Faces associated with vertex A
								//surfA->vertices[ptr->deg_pointA - 1].faceIDs.size();
								for (int face_count = 0; face_count < (int)surfA->vertices[deg_pointA - 1].faceIDs.size(); face_count++)
								{
									if (subptr->deg_curveA == 0 && subptr->deg_pointA == 0 && subptr->faceA->ID == surfA->vertices[deg_pointA - 1].faceIDs[face_count]
										&& subptr->deg_pointB == deg_pointB)
									{
										contact_pairs[c]->eligible = false;
										checked = true;
									}
								}

								//Edges associated with vertex B
								//surfB->vertices[ptr->deg_pointB - 1].edgeIDs.size();
								for (int edge_count = 0; edge_count < (int)surfB->vertices[deg_pointB - 1].edgeIDs.size(); edge_count++)
								{
									if (subptr->deg_curveB == surfB->vertices[deg_pointB - 1].edgeIDs[edge_count]
										&& subptr->deg_pointA == deg_pointA)
									{
										contact_pairs[c]->eligible = false;
										checked = true;
									}
								}
								//Faces associated with vertex B
								//surfB->vertices[ptr->deg_pointB - 1].faceIDs.size();
								for (int face_count = 0; face_count < (int)surfB->vertices[deg_pointB - 1].faceIDs.size(); face_count++)
								{
									if (subptr->deg_curveB == 0 && subptr->deg_pointB == 0 && subptr->faceB->ID == surfB->vertices[deg_pointB - 1].faceIDs[face_count]
										&& subptr->deg_pointA == deg_pointA)
									{
										contact_pairs[c]->eligible = false;
										checked = true;
									}
								}
							}
						}
					}
				}
				///////////////////////////////////////////////////////////////////////////////
			}
		}
	}
}