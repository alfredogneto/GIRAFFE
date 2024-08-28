#include "STLBoundary.h"

#include "BoundingSphere.h"
#include "CADData.h"
#include "STLSurface.h"
#include "Node.h"
#include "CoordinateSystem.h"
#include "GeneralContactSearch.h"
#include "BoundingTriangularBox.h"
#include "BoundingCylinder.h"
#include "Particle.h"
#include "VEMPolyhedron.h"
#include "Encoding.h"
#include "TriangularFace.h"
#include "MatrixFloat.h"
#include"Database.h"
//Variáveis globais
extern
Database db;


STLBoundary::STLBoundary()
{
	Alloc();
	bv = new BoundingSphere();
	CADDATA_ID = 0;
	bv_factor = 0.0f;
	inc_len_factor = 0.1f;

	x0f = new MatrixFloat(3);
	Q0f = new MatrixFloat(3, 3);
}


STLBoundary::~STLBoundary()
{
	Free();

	delete bv;
	if (sub_bv != NULL)
	{
		for (int i = 0; i < n_sub_bv; i++)
		{
			delete sub_bv[i];
		}
		delete[] sub_bv;
	}
	delete x0f;
	delete Q0f;
}

bool STLBoundary::Check()
{
	if (node > db.number_nodes)
		return false;
	if (cs > db.number_CS)
		return false;
	if (CADDATA_ID > db.number_cad_data)
		return false;
	if (typeid(*db.cad_data[CADDATA_ID - 1]) != typeid(STLSurface))
		return false;
	return true;
}

bool STLBoundary::Read(FILE *f)
{
	char s[1000];
	fscanf(f, "%s", s);
	number = atoi(s);

	fscanf(f, "%s", s);
	if (!strcmp(s, "Mat"))
	{
		fscanf(f, "%s", s);
		material = atoi(s);
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "CS"))
	{
		fscanf(f, "%s", s);
		cs = atoi(s);
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "CADData"))
	{
		fscanf(f, "%s", s);
		CADDATA_ID = atoi(s);
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "Node"))
	{
		fscanf(f, "%s", s);
		node = atoi(s);
	}
	else
		return false;
	return true;
}

void STLBoundary::Write(FILE *f)
{
	fprintf(f, "STLBoundary\t%d\tMat\t%d\tCS\t%d\tCADData\t%d\tNode\t%d\n",
		number,
		material,
		cs,
		CADDATA_ID,
		node
	);
}

void STLBoundary::WriteModifyingParameters(FILE *f, int e_material, int e_node, int e_number, int e_cs)
{
	fprintf(f, "STLBoundary\t%d\tMat\t%d\tCS\t%d\tCADData\t%d\tNode\t%d\n",
		e_number,
		e_material,
		e_cs,
		CADDATA_ID,
		e_node
	);
}

//Explicit
void STLBoundary::InitialEvaluations()
{
	//Seta variáveis nodais associadas com o tratamento que essa partícula vai empregar
	*db.nodes[node - 1]->Q0 = *Q0;
	db.nodes[node - 1]->flag_material_description = true;
	db.nodes[node - 1]->flag_pseudo_moment = false;
}

void STLBoundary::PreCalc()
{
	////////////////////////////////////////////////////////////////
	//Matriz para transformação de coordenadas  local-global
	//vlocal = (sistema do CAD)
	// Q0*vlocal = vglobal
	////////////////////////////////////////////////////////////////
	(*Q0)(0, 0) = (*db.CS[cs - 1]->E1)(0, 0);
	(*Q0)(1, 0) = (*db.CS[cs - 1]->E1)(1, 0);
	(*Q0)(2, 0) = (*db.CS[cs - 1]->E1)(2, 0);
	(*Q0)(0, 1) = (*db.CS[cs - 1]->E2)(0, 0);
	(*Q0)(1, 1) = (*db.CS[cs - 1]->E2)(1, 0);
	(*Q0)(2, 1) = (*db.CS[cs - 1]->E2)(2, 0);
	(*Q0)(0, 2) = (*db.CS[cs - 1]->E3)(0, 0);
	(*Q0)(1, 2) = (*db.CS[cs - 1]->E3)(1, 0);
	(*Q0)(2, 2) = (*db.CS[cs - 1]->E3)(2, 0);

	//////////////////////////////////////Bounding Volumes////////////////////////////////////
	if (db.gcs_exist)
		inc_len_factor = db.gcs->inc_len_factor;
	//PreCalc do CAD é chamado antes do PreCalc das partículas
	BoundingSphere* ptr_bv1 = static_cast<BoundingSphere*>(bv);
	STLSurface* ptr_cad = static_cast<STLSurface*>(db.cad_data[CADDATA_ID - 1]);
	//Initial settings of the spherical bounding volume
	ptr_bv1->radius = (1.0f + inc_len_factor) * ptr_cad->radius;
	ptr_bv1->ref_radius = (1.0f + inc_len_factor) * ptr_cad->radius;
	ptr_bv1->size = ptr_bv1->radius;
	//Associated entity:
	ptr_bv1->associated_type = 'B';
	ptr_bv1->associated_ID = number;
	ptr_bv1->associated_sub_ID = 0;
	//Sub bounding volumes (faces+edges+vertices)
	n_sub_bv = ptr_cad->n_faces + (int)ptr_cad->edges.size() + (int)ptr_cad->vertices.size();
	sub_bv = new BoundingVolume*[n_sub_bv];
	//BV on faces
	BoundingTriangularBox* ptr_bv2;
	for (int i = 0; i < ptr_cad->n_faces; i++)
	{
		sub_bv[i] = new BoundingTriangularBox();
		ptr_bv2 = static_cast<BoundingTriangularBox*>(sub_bv[i]);
		ptr_bv2->thickness = inc_len_factor * ptr_cad->faces[i]->radius;
		ptr_bv2->ref_thickness = inc_len_factor * ptr_cad->faces[i]->radius;
		ptr_bv2->inc_len_factor = inc_len_factor;
		float templen = (1.0f + inc_len_factor)*ptr_cad->faces[i]->radius;
		ptr_bv2->size = sqrt(templen *templen + ptr_bv2->thickness*ptr_bv2->thickness);
		//printf("Size %.6e\n", ptr_bv2->size);
		//Associated entity:
		ptr_bv2->associated_type = 'B';
		ptr_bv2->associated_ID = number;
		ptr_bv2->associated_sub_ID = i + 1;
	}
	//BV on edges
	BoundingCylinder* ptr_bv3;
	for (int i = ptr_cad->n_faces; i < ptr_cad->n_faces + (int)ptr_cad->edges.size(); i++)
	{
		sub_bv[i] = new BoundingCylinder();
		//Determining the radius of the cylinder (smallest radius related to a neighboring face)
		float temp_radius = FLT_MAX;
		for (int j = 0; j < ptr_cad->edges[i - ptr_cad->n_faces].faceIDs.size(); j++)
		{
			if (ptr_cad->faces[ptr_cad->edges[i - ptr_cad->n_faces].faceIDs[j] - 1]->radius < temp_radius)
				temp_radius = ptr_cad->faces[ptr_cad->edges[i - ptr_cad->n_faces].faceIDs[j] - 1]->radius;
		}
		ptr_bv3 = static_cast<BoundingCylinder*>(sub_bv[i]);
		ptr_bv3->radius = inc_len_factor * temp_radius;
		ptr_bv3->ref_radius = inc_len_factor * temp_radius;
		ptr_bv3->inc_len_factor = inc_len_factor;
		float templen = 0.5f*(1.0f + inc_len_factor)*ptr_cad->edges[i - ptr_cad->n_faces].length + ptr_bv3->radius;
		ptr_bv3->size = sqrt(templen*templen + ptr_bv3->radius*ptr_bv3->radius);
		//Associated entity:
		ptr_bv3->associated_type = 'B';
		ptr_bv3->associated_ID = number;
		ptr_bv3->associated_sub_ID = i + 1;
	}
	//BV on vertices
	BoundingSphere* ptr_bv4;
	for (int i = ptr_cad->n_faces + (int)ptr_cad->edges.size(); i < ptr_cad->n_faces + (int)ptr_cad->edges.size() + (int)ptr_cad->vertices.size(); i++)
	{
		sub_bv[i] = new BoundingSphere();
		//Determining the radius of the cylinder (smallest radius related to a neighboring face)
		float temp_radius = FLT_MAX;
		for (int j = 0; j < ptr_cad->vertices[i - ptr_cad->n_faces - (int)ptr_cad->edges.size()].faceIDs.size(); j++)
		{
			if (ptr_cad->faces[ptr_cad->vertices[i - ptr_cad->n_faces - (int)ptr_cad->edges.size()].faceIDs[j] - 1]->radius < temp_radius)
				temp_radius = ptr_cad->faces[ptr_cad->vertices[i - ptr_cad->n_faces - (int)ptr_cad->edges.size()].faceIDs[j] - 1]->radius;
		}
		ptr_bv4 = static_cast<BoundingSphere*>(sub_bv[i]);
		ptr_bv4->radius = inc_len_factor * temp_radius;
		ptr_bv4->ref_radius = inc_len_factor * temp_radius;
		ptr_bv4->size = ptr_bv4->radius;
		//Associated entity:
		ptr_bv4->associated_type = 'B';
		ptr_bv4->associated_ID = number;
		ptr_bv4->associated_sub_ID = i + 1;
	}
	UpdateVariables();
	UpdateBoundingVolumes();
	SaveLagrange();
}

void STLBoundary::UpdateBoundingVolumes()
{
	//Evaluating x0f and Q0f

	//Conversão da matriz de rotação para single precision
	*x0f = *x0ip;
	*Q0f = *Qip;

	//Updating the main bounding volume
	BoundingSphere* ptr_bv1 = static_cast<BoundingSphere*>(bv);
	*ptr_bv1->center = *x0f;

	//Updating the sub bounding volumes
	STLSurface* ptr_cad = static_cast<STLSurface*>(db.cad_data[CADDATA_ID - 1]);
	//BV on faces
	BoundingTriangularBox* ptr_bv2;
	for (int i = 0; i < ptr_cad->n_faces; i++)
	{
		ptr_bv2 = static_cast<BoundingTriangularBox*>(sub_bv[i]);
		*ptr_bv2->x0 = *x0f + (*Q0f) * (*ptr_cad->vertices[ptr_cad->faces[i]->verticesIDs[0] - 1].coord_float);
		*ptr_bv2->x1 = *x0f + (*Q0f) * (*ptr_cad->vertices[ptr_cad->faces[i]->verticesIDs[1] - 1].coord_float);
		*ptr_bv2->x2 = *x0f + (*Q0f) * (*ptr_cad->vertices[ptr_cad->faces[i]->verticesIDs[2] - 1].coord_float);
	}
	//BV on edges
	BoundingCylinder* ptr_bv3;
	for (int i = ptr_cad->n_faces; i < ptr_cad->n_faces + (int)ptr_cad->edges.size(); i++)
	{
		ptr_bv3 = static_cast<BoundingCylinder*>(sub_bv[i]);
		*ptr_bv3->xb = *x0f + (*Q0f) * (*ptr_cad->vertices[ptr_cad->edges[i - ptr_cad->n_faces].verticesIDs[0] - 1].coord_float);
		*ptr_bv3->xt = *x0f + (*Q0f) * (*ptr_cad->vertices[ptr_cad->edges[i - ptr_cad->n_faces].verticesIDs[1] - 1].coord_float);
		//ptr_bv3->Report();
	}
	//BV on vertices
	BoundingSphere* ptr_bv4;
	for (int i = ptr_cad->n_faces + (int)ptr_cad->edges.size(); i < ptr_cad->n_faces + (int)ptr_cad->edges.size() + (int)ptr_cad->vertices.size(); i++)
	{
		ptr_bv4 = static_cast<BoundingSphere*>(sub_bv[i]);
		*ptr_bv4->center = *x0f + (*Q0f) * (*ptr_cad->vertices[i - ptr_cad->n_faces - (int)ptr_cad->edges.size()].coord_float);
	}
}

void STLBoundary::SaveLagrange()
{
	//Do not call UpdateVaribles because node data are updated first!
	//Particle variables
	*x0i = *x0ip;
	*Qi = *Qip;

	//Saving bounding volumes
	bv->SaveConfiguration();
	for (int i = 0; i < n_sub_bv; i++)
		sub_bv[i]->SaveConfiguration();
}
void STLBoundary::WriteVTK_XMLBase(FILE *f)
{

}

void STLBoundary::WriteVTK_XMLRender(FILE *f)
{
	Matrix pos(3);
	pos(0, 0) = db.nodes[node - 1]->copy_coordinates[0];
	pos(1, 0) = db.nodes[node - 1]->copy_coordinates[1];
	pos(2, 0) = db.nodes[node - 1]->copy_coordinates[2];
	
	if (db.nodes[node - 1]->flag_material_description)
		db.cad_data[CADDATA_ID - 1]->WriteVTK_XMLRender(f, pos, (*Q0) * (*db.nodes[node - 1]->Q), number);
	else
		db.cad_data[CADDATA_ID - 1]->WriteVTK_XMLRender(f, pos, (*db.nodes[node - 1]->Q)*(*Q0), number);
}
