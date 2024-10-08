#include "ContactPolyhedronPolyhedron.h"

#include "Polyhedron.h"
#include "STLSurface.h"
#include "TriangularFace.h"
#include "Matrix.h"
#include "RigidTriangularFace_RigidTriangularFace.h"
#include "Node.h"
#include "GeneralContactSearch.h"
#include "SurfacePairGeneralContact.h"

#include "Database.h"
//Variaveis globais
extern
Database db;

ContactPolyhedronPolyhedron::ContactPolyhedronPolyhedron()
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

ContactPolyhedronPolyhedron::~ContactPolyhedronPolyhedron()
{
	for (int i = 0; i < contact_pairs.size(); i++)
		delete contact_pairs[i];
	contact_pairs.clear();

	//Cleaning
	delete I3;
}

void ContactPolyhedronPolyhedron::PreCalc()
{
	pA = static_cast<Polyhedron*>(db.particles[index1]);
	pB = static_cast<Polyhedron*>(db.particles[index2]);
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


	//pointers to evaluate rigid body motion of the particle
	QAp = pA->Qip;
	QBp = pB->Qip;
	x0Ap = pA->x0ip;
	x0Bp = pB->x0ip;

	QAi = pA->Qi;
	QBi = pB->Qi;
	x0Ai = pA->x0i;
	x0Bi = pB->x0i;

	Q0A = pA->Q0;
	Q0B = pB->Q0;

	//Tests to exit with no creation of a contact pair: edge-edge
	if (deg_curveA != 0 && surfA->edges[deg_curveA - 1].concave_indicator == 0)
		return;
	if (deg_curveB != 0 && surfB->edges[deg_curveB - 1].concave_indicator == 0)
		return;
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

void ContactPolyhedronPolyhedron::InsertNewContact(int deg_pointA, int deg_curveA, int deg_pointB, int deg_curveB, int faceAID, int faceBID)
{
	RigidTriangularFace_RigidTriangularFace* temp = new RigidTriangularFace_RigidTriangularFace();
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
	temp->ptrx0Ai = x0Ai;
	temp->ptrx0Bi = x0Bi;
	temp->ptrQAi = QAi;
	temp->ptrQBi = QBi;
	temp->ptrx0Ap = x0Ap;
	temp->ptrx0Bp = x0Bp;
	temp->ptrQAp = QAp;
	temp->ptrQBp = QBp;
	temp->node_A = pA->node;
	temp->node_B = pB->node;
	temp->material_A = pA->material;
	temp->material_B = pB->material;
	temp->CAD_AID = pA->CADDATA_ID;
	temp->CAD_BID = pB->CADDATA_ID;
	temp->ptrQ0A = Q0A;
	temp->ptrQ0B = Q0B;

	//temp->I3 = I3;
	temp->SetActive();
	contact_pairs.push_back(temp);
}

int ContactPolyhedronPolyhedron::CreateSurfacePair(int deg_pointA, int deg_curveA, int deg_pointB, int deg_curveB, int faceAID, int faceBID)
{
	int list_index = -1;
	for (int i = 0; i < contact_pairs.size(); i++)
	{
		RigidTriangularFace_RigidTriangularFace* ptr = static_cast<RigidTriangularFace_RigidTriangularFace*>(contact_pairs[i]);
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
	//ProcessSurfacePair(list_index);
	
	return list_index;
}
void ContactPolyhedronPolyhedron::ProcessSurfacePair(int list_index)
{
	contact_pairs[list_index]->eligible = false;
	contact_pairs[list_index]->SolveLCP();
	contact_pairs[list_index]->EvaluateNormalGap();
		
}
void ContactPolyhedronPolyhedron::ProcessSurfacePairs()
{
	for (int i = 0; i < (int)contact_pairs.size(); i++)
		ProcessSurfacePair(i);
}


void ContactPolyhedronPolyhedron::FinalProcessSurfacePairExplicit(int list_index, double t)
{
	contact_pairs[list_index]->eligible = false;
	contact_pairs[list_index]->SolveLCP();
	contact_pairs[list_index]->EvaluateNormalGap();//here the normal of contact is updated, which ensures a successful test for penetration in HaveErrors function
	if (contact_pairs[list_index]->eligible == true)
	{
		contact_pairs[list_index]->Alloc();
		contact_pairs[list_index]->SetVariablesExplicit(t);		//Sets variables for next evaluations
		contact_pairs[list_index]->FinalUpdateExplicit(t);
	}
}
void ContactPolyhedronPolyhedron::FinalProcessSurfacePairsExplicit(double t)
{
	for (int i = 0; i < (int)contact_pairs.size(); i++)
		FinalProcessSurfacePairExplicit(i,t);
}

void ContactPolyhedronPolyhedron::MountGlobalExplicit()
{
	//Variaveis temporarias para salvar a indexa��o global dos graus de liberdade a serem setados na matriz de rigidez global
	//Obs: aqui os contatos s�o subtraidos pois no explicito ja esta isolada a acelera��o e aqui entram como se fosse for�as externas, sendo necessaria a invers�o do sinal
	int GL_global_1 = 0;
	double anterior = 0;
	for (int c = 0; c < contact_pairs.size(); c++)
	{
		if (contact_pairs[c]->GetActive())
		{
			if (contact_pairs[c]->eligible)
			{
				for (int i = 0; i < 12; i++)
				{
					if (i < 6)
						GL_global_1 = db.nodes[db.particles[index1]->node - 1]->GLs[i];
					else
						GL_global_1 = db.nodes[db.particles[index2]->node - 1]->GLs[i - 6];
					//Caso o grau de liberdade seja livre:
					if (GL_global_1 > 0)
					{
						anterior = db.global_P_A(GL_global_1 - 1, 0);
						db.global_P_A(GL_global_1 - 1, 0) = anterior - contact_pairs[c]->Rc[i];
						anterior = db.global_I_A(GL_global_1 - 1, 0);
						db.global_I_A(GL_global_1 - 1, 0) = anterior - contact_pairs[c]->Rc[i];
					}
					else
					{
						if (GL_global_1 != 0)//se o GL e ativo
						{
							anterior = db.global_P_B(-GL_global_1 - 1, 0);
							db.global_P_B(-GL_global_1 - 1, 0) = anterior - contact_pairs[c]->Rc[i];
						}
					}
				}
			}
		}
	}
}

void ContactPolyhedronPolyhedron::MountGlobal()
{
	//Variaveis temporarias para salvar a indexa��o global dos graus de liberdade a serem setados na matriz de rigidez global
	int GL_global_1 = 0;
	int GL_global_2 = 0;
	double anterior = 0;
	for (int c = 0; c < contact_pairs.size(); c++)
	{
		if (contact_pairs[c]->GetActive())
		{
			if (contact_pairs[c]->eligible)
			{
				for (int i = 0; i < 12; i++)
				{
					if (i < 6)
						GL_global_1 = db.nodes[db.particles[index1]->node - 1]->GLs[i];
					else
						GL_global_1 = db.nodes[db.particles[index2]->node - 1]->GLs[i - 6];
					//Caso o grau de liberdade seja livre:
					if (GL_global_1 > 0)
					{
						anterior = db.global_P_A(GL_global_1 - 1, 0);
						db.global_P_A(GL_global_1 - 1, 0) = anterior + contact_pairs[c]->Rc[i];
						anterior = db.global_I_A(GL_global_1 - 1, 0);
						db.global_I_A(GL_global_1 - 1, 0) = anterior + contact_pairs[c]->Rc[i];
					}
					else
					{
						if (GL_global_1 != 0)//se o GL e ativo
						{
							anterior = db.global_P_B(-GL_global_1 - 1, 0);
							db.global_P_B(-GL_global_1 - 1, 0) = anterior + contact_pairs[c]->Rc[i];
						}
					}
					for (int j = 0; j < 12; j++)
					{
						if (j < 6)
							GL_global_2 = db.nodes[db.particles[index1]->node - 1]->GLs[j];
						else
							GL_global_2 = db.nodes[db.particles[index2]->node - 1]->GLs[j - 6];
						//Caso os graus de liberdade sejam ambos livres (Matriz Kaa)
						if (GL_global_1 > 0 && GL_global_2 > 0)
							db.global_stiffness_AA.setValue(GL_global_1 - 1, GL_global_2 - 1, contact_pairs[c]->Kc[i][j]);
						//Caso os graus de liberdade sejam ambos fixos (Matriz Kbb)
						if (GL_global_1 < 0 && GL_global_2 < 0)
							db.global_stiffness_BB.setValue(-GL_global_1 - 1, -GL_global_2 - 1, contact_pairs[c]->Kc[i][j]);
						//Caso os graus de liberdade sejam livre e fixo (Matriz Kab)
						if (GL_global_1 > 0 && GL_global_2 < 0)
							db.global_stiffness_AB.setValue(GL_global_1 - 1, -GL_global_2 - 1, contact_pairs[c]->Kc[i][j]);
						//Caso os graus de liberdade sejam fixo e livre (Matriz Kba)
						if (GL_global_1 < 0 && GL_global_2 > 0)
							db.global_stiffness_BA.setValue(-GL_global_1 - 1, GL_global_2 - 1, contact_pairs[c]->Kc[i][j]);
					}
				}
			}
		}
	}
}

void ContactPolyhedronPolyhedron::ProcessContactHierarchy()
{
	for (int c = 0; c < (int)contact_pairs.size(); c++)
	{
		//If the pair is eligible, checks its hierarchy
		if (cur_active && contact_pairs[c]->eligible)
		{
			//Indexes:	particle 1: index1
			//			particle 2: index2
			//Pointers
			Polyhedron* pA = static_cast<Polyhedron*>(db.particles[index1]);
			Polyhedron* pB = static_cast<Polyhedron*>(db.particles[index2]);
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
									ContactPolyhedronPolyhedron* subptr = static_cast<ContactPolyhedronPolyhedron*>(db.gcs->contactPP_list[index1][subcont]);
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
									ContactPolyhedronPolyhedron* subptr = static_cast<ContactPolyhedronPolyhedron*>(db.gcs->contactPP_list[index1][subcont]);
									//Test: subptr->faceA->ID == face1ID  & face-point contact
									if (subptr->deg_curveA == 0 && subptr->deg_pointA == 0 && subptr->faceA->ID == face1ID && subptr->deg_pointB ==deg_pointB)
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
								ContactPolyhedronPolyhedron* subptr = static_cast<ContactPolyhedronPolyhedron*>(db.gcs->contactPP_list[index1][subcont]);

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
									ContactPolyhedronPolyhedron* subptr = static_cast<ContactPolyhedronPolyhedron*>(db.gcs->contactPP_list[index1][subcont]);
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
									ContactPolyhedronPolyhedron* subptr = static_cast<ContactPolyhedronPolyhedron*>(db.gcs->contactPP_list[index1][subcont]);
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
								ContactPolyhedronPolyhedron* subptr = static_cast<ContactPolyhedronPolyhedron*>(db.gcs->contactPP_list[index1][subcont]);

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
								ContactPolyhedronPolyhedron* subptr = static_cast<ContactPolyhedronPolyhedron*>(db.gcs->contactPP_list[index1][subcont]);

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