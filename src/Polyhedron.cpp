#include "Polyhedron.h"

#include "BoundingSphere.h"
#include "CADData.h"
#include "MatrixFloat.h"
#include "STLSurface.h"
#include "CoordinateSystem.h"
#include "Material.h"
#include "GeneralContactSearch.h"
#include "Interface_1.h"
#include "BoundingTriangularBox.h"
#include "BoundingCylinder.h"
#include "Environment.h"
#include "Node.h"
#include "ContactParticleParticle.h"
#include "Dynamic.h"
#include "TriangularFace.h"

//extern
//FILE *fdebug;

#include"Database.h"
//Variaveis globais
extern
Database db;

Polyhedron::Polyhedron()
{
	Alloc();
	bv = new BoundingSphere();
	n_sub_bv = 0;
	CADDATA_ID = 0;
	sprintf(type_name, "Polyhedron");
	//bv_factor = 0.1f;
	strain_energy = 0.0;
	potential_g_energy = 0.0;
	inc_len_factor = 0.1f;

	x0f = new MatrixFloat(3);
	Q0f = new MatrixFloat(3, 3);

	pQi = new double*[3];
	for (int i = 0; i < 3; i++)
		pQi[i] = new double[3];

	allocedMassFragmented = false;

	GLA = NULL;
	GLB = NULL;
}

Polyhedron::~Polyhedron()
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
	
	for (int i = 0; i < 3; i++)
		delete pQi[i];
	delete[] pQi;
}

bool Polyhedron::Check()
{
	if (node > db.number_nodes)
		return false;
	if (cs > db.number_CS)
		return false;
	if (material > db.number_materials)
		return false;
	if (CADDATA_ID > db.number_cad_data)
		return false;
	if (typeid(*db.cad_data[CADDATA_ID - 1]) != typeid(STLSurface))
		return false;
	return true;
}

bool Polyhedron::Read(FILE *f)
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

void Polyhedron::Write(FILE *f)
{
	fprintf(f, "Polyhedron\t%d\tMat\t%d\tCS\t%d\tCADData\t%d\tNode\t%d\n",
		number,
		material,
		cs,
		CADDATA_ID,
		node
	);
}

void Polyhedron::WriteModifyingParameters(FILE *f, int e_material, int e_node, int e_number, int e_cs)
{
	fprintf(f, "Polyhedron\t%d\tMat\t%d\tCS\t%d\tCADData\t%d\tNode\t%d\n",
		e_number,
		e_material,
		e_cs,
		CADDATA_ID,
		e_node
	);
}

void Polyhedron::PreCalc()
{
	////////////////////////////////////////////////////////////////
	//Matriz para transforma��o de coordenadas  local-global
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

	//Calculando valores do br. Feito apenas uma vez no inicio. 
	//A origem do CAD coincide com o Polo, ent�o G-O = G
	//Como G e dado no sistema do CAD, ent�o G precisa ser transformado
	Matrix local_b = *db.cad_data[CADDATA_ID - 1]->G;
	
	(*mbr) = (*Q0) * local_b;
	br[0] = (*mbr)(0, 0);
	br[1] = (*mbr)(1, 0);
	br[2] = (*mbr)(2, 0);

	//Atualizando J_O para o sistema global e copiando valores para Jr. Feito apenas uma vez no inicio
	*mJr = db.materials[material - 1]->rho * (*Q0) * (*db.cad_data[CADDATA_ID - 1]->J_O)*transp(*Q0);
	mJr->MatrixToPtr(Jr, 3);
	mass = db.cad_data[CADDATA_ID - 1]->volume * db.materials[material - 1]->rho;

	//Explicit
	*mbrlocal = *db.cad_data[CADDATA_ID - 1]->G;
	*mJrlocal = (db.materials[material - 1]->rho) * (*db.cad_data[CADDATA_ID - 1]->J_O);
	
	//Skew(b) - including the mass
	Matrix mB = mass * skew(*mbrlocal);

	//Mass matrix evaluation - explicit
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			(*MassMatrix)(i, j) = mass * (*I3)(i, j);
			(*MassMatrix)(i + 0, j + 3) = -1.0 * mB(i, j);
			(*MassMatrix)(i + 3, j + 0) = mB(i, j);
			(*MassMatrix)(i + 3, j + 3) = (*mJrlocal)(i, j);
		}
	}

	//////////////////////////////////////Bounding Volumes////////////////////////////////////
	if (db.gcs_exist)
		inc_len_factor = db.gcs->inc_len_factor;

	//Computing the largest gnb available, from interface laws (to avoid smaller BVs of small geometric entities)
	double largest_gnb = 0.0;
	for (int i = 0; i < db.number_contactinterfaces; i++)
	{
		Interface_1* ptr = static_cast<Interface_1*>(db.contactinterfaces[i]);
		if (ptr->gnb > largest_gnb)
			largest_gnb = ptr->gnb;
	}
	//PreCalc do CAD e chamado antes do PreCalc das particulas
	BoundingSphere* ptr_bv1 = static_cast<BoundingSphere*>(bv);
	STLSurface* ptr_cad = static_cast<STLSurface*>(db.cad_data[CADDATA_ID - 1]);
	//Initial settings of the spherical bounding volume
	ptr_bv1->radius = (1.0f + inc_len_factor) * ptr_cad->radius;
	ptr_bv1->ref_radius = (1.0f + inc_len_factor) * ptr_cad->radius;
	ptr_bv1->size = ptr_bv1->radius;
	//fprintf(fdebug, "Polyhedron %d\t%.6e\t", number,ptr_bv1->size);
	//Associated entity:
	ptr_bv1->associated_type = 'P';
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
		//Test of small radius
		float eff_radius = ptr_cad->faces[i]->radius;
		if (eff_radius < (1.0f + inc_len_factor)*(float)largest_gnb)
			eff_radius = (1.0f + inc_len_factor)*(float)largest_gnb;
		//end of Test of small radius
		ptr_bv2->thickness = inc_len_factor * eff_radius;
		ptr_bv2->ref_thickness = inc_len_factor * eff_radius;
		ptr_bv2->inc_len_factor = inc_len_factor;
		float templen = (1.0f + inc_len_factor)*eff_radius;
		ptr_bv2->size = sqrt(templen *templen + ptr_bv2->thickness*ptr_bv2->thickness);
		//fprintf(fdebug, "%.6e\t", ptr_bv2->size);
		//Associated entity:
		ptr_bv2->associated_type = 'P';
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
		//Test of small radius
		if (temp_radius < (1.0f + inc_len_factor)*(float)largest_gnb)
			temp_radius = (1.0f + inc_len_factor)*(float)largest_gnb;
		//end of Test of small radius
		ptr_bv3 = static_cast<BoundingCylinder*>(sub_bv[i]);
		ptr_bv3->radius = inc_len_factor * temp_radius;
		ptr_bv3->ref_radius = inc_len_factor * temp_radius;
		ptr_bv3->inc_len_factor = inc_len_factor;
		float templen = 0.5f*(1.0f + inc_len_factor)*ptr_cad->edges[i - ptr_cad->n_faces].length + ptr_bv3->radius;
		ptr_bv3->size = sqrt(templen*templen + ptr_bv3->radius*ptr_bv3->radius);
		//fprintf(fdebug, "%.6e\t", ptr_bv3->size);
		//Associated entity:
		ptr_bv3->associated_type = 'P';
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
		//Test of small radius
		if (temp_radius < (1.0f + inc_len_factor)*(float)largest_gnb)
			temp_radius = (1.0f + inc_len_factor)*(float)largest_gnb;
		//end of Test of small radius

		ptr_bv4 = static_cast<BoundingSphere*>(sub_bv[i]);
		ptr_bv4->radius = inc_len_factor * temp_radius;
		ptr_bv4->ref_radius = inc_len_factor * temp_radius;
		ptr_bv4->size = ptr_bv4->radius;
		//fprintf(fdebug, "%.6e\t", ptr_bv4->size);
		//Associated entity:
		ptr_bv4->associated_type = 'P';
		ptr_bv4->associated_ID = number;
		ptr_bv4->associated_sub_ID = i + 1;
	}
	//fprintf(fdebug, "\n");
	UpdateVariables();
	UpdateBoundingVolumes();
	SaveLagrange();
}

void Polyhedron::UpdateBoundingVolumes()
{
	//Evaluating x0f and Q0f
	
	//Convers�o da matriz de rota��o para single precision
	*x0f = *x0ip;
	*Q0f = *Qip;
	 
	//Updating the main bounding volume
	BoundingSphere* ptr_bv1 = static_cast<BoundingSphere*>(bv);
	*ptr_bv1->center = *x0f;
	//ptr_bv1->Report();
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
		//ptr_bv2->Report();
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
		//ptr_bv4->Report();
	}
}

void Polyhedron::SaveLagrange()
{
	//Do not call UpdateVaribles because node data are updated first!
	//Particle variables
	*x0i = *x0ip;
	*Qi = *Qip;

	//Saving bounding volumes
	bv->SaveConfiguration();
	for (int i = 0; i < n_sub_bv; i++)
		sub_bv[i]->SaveConfiguration();

	///////////////////////Potential gravitational energy//////////////////////////////////////////
	Matrix center(3);	//Updated center of mass of the body
	Matrix g(3);		//Vetor gravidade
	if (db.environment_exist == true)
	{
		//Se existe campo gravitacional
		if (db.environment->g_exist == true)
		{
			//Gravity field
			double l_factor = db.environment->bool_g.GetLinearFactorAtCurrentTime();
			g(0, 0) = l_factor * db.environment->G(0, 0);
			g(1, 0) = l_factor * db.environment->G(1, 0);
			g(2, 0) = l_factor * db.environment->G(2, 0);
			//Center of mass - approximated value
			center = *x0ip + (*Qip) * (*db.cad_data[CADDATA_ID - 1]->G);
		}
	}
	potential_g_energy = - mass * dot(g, center);

	Qi->MatrixToPtr(pQi, 3);
}
void Polyhedron::WriteVTK_XMLBase(FILE *f)
{

}

void Polyhedron::WriteVTK_XMLRender(FILE *f)
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

bool Polyhedron::CheckInside(Matrix& point, int other_particle_ID)
{
	STLSurface* ptr_cad = static_cast<STLSurface*>(db.cad_data[CADDATA_ID - 1]);
	Matrix normal(3);
	Matrix localt1(3);
	Matrix localt2(3);
	int ID0, ID1, ID2;
	Matrix c1(3);
	Matrix centroid(3);
	bool return_value = true;		//inside default
	bool convex = true;				//convex default
	int local_index1, local_index2;
	int face_index;

	//Check is based only on active contacts - already accounting for geometric proximity - avoid searching on all triangles of the particle
	if (db.gcs_exist == true)
	{
		if (number < other_particle_ID)
		{
			local_index1 = number - 1;
			local_index2 = other_particle_ID - 1;
		}
		else
		{
			local_index1 = other_particle_ID - 1;
			local_index2 = number - 1;
		}

		//test convex/concave
		for (int i = 0; i < db.gcs->contactPP_list[local_index1].size() && return_value == true; i++)
		{
			if (db.gcs->contactPP_list[local_index1][i]->index2 == local_index2)
			{
				if (number < other_particle_ID)
					face_index = db.gcs->contactPP_list[local_index1][i]->sub_index1;
				else
					face_index = db.gcs->contactPP_list[local_index1][i]->sub_index2;

				if (ptr_cad->faces[face_index]->has_concave_edge)
					convex = false;
			}
		}
		if (convex)
		{
			//loop on contact pairs associated with the particle number
			
			for (int i = 0; i < db.gcs->contactPP_list[local_index1].size() && return_value == true; i++)
			{
				if (db.gcs->contactPP_list[local_index1][i]->index2 == local_index2)
				{
					if (number < other_particle_ID)
						face_index = db.gcs->contactPP_list[local_index1][i]->sub_index1;
					else
						face_index = db.gcs->contactPP_list[local_index1][i]->sub_index2;
					//Centroid
					centroid(0, 0) = (double)(*ptr_cad->faces[face_index]->centroid)(0, 0);
					centroid(1, 0) = (double)(*ptr_cad->faces[face_index]->centroid)(1, 0);
					centroid(2, 0) = (double)(*ptr_cad->faces[face_index]->centroid)(2, 0);
					//Coordinates (point - centroid)
					c1 = point - *x0ip - (*Qip)*centroid;
					//Normal
					ID0 = ptr_cad->faces[face_index]->verticesIDs[0];
					ID1 = ptr_cad->faces[face_index]->verticesIDs[1];
					ID2 = ptr_cad->faces[face_index]->verticesIDs[2];
					localt1 = *ptr_cad->vertices[ID1 - 1].coord_double - *ptr_cad->vertices[ID0 - 1].coord_double;
					localt2 = *ptr_cad->vertices[ID2 - 1].coord_double - *ptr_cad->vertices[ID0 - 1].coord_double;
					normal = 1.0 / norm(cross(localt1, localt2))*cross(localt1, localt2);
					normal = (*Qip) * normal;
					
					if (dot(c1, normal) > DBL_EPSILON)
						return_value = false;
				}
			}
		}
		else
		{
			//TODO - based on Jordan Theorem
			return_value = false;
		}
		return return_value;
	}
	else
		return return_value;

}

bool Polyhedron::CheckInsideEdge(Matrix& point, int edge_ID, double tol)
{
	STLSurface* ptr_cad = static_cast<STLSurface*>(db.cad_data[CADDATA_ID - 1]);

	//Checando se e uma aresta c�ncava. Se assim o for, retorna false
	if (ptr_cad->edges[edge_ID - 1].concave_indicator == 1)
		return false;
	else
	{
		Matrix n1(3);
		Matrix n2(3);
		Matrix t(3);
		Matrix O(3);
		//Vertices da aresta em quest�o
		O = *x0ip + (*Qip)*(*ptr_cad->vertices[ptr_cad->edges[edge_ID - 1].verticesIDs[0] - 1].coord_double);
		n1 = (*Qip)*(*ptr_cad->faces[ptr_cad->edges[edge_ID - 1].faceIDs[0] - 1]->normal);
		n2 = (*Qip)*(*ptr_cad->faces[ptr_cad->edges[edge_ID - 1].faceIDs[1] - 1]->normal);
		t = (1.0/norm(cross(n1, n2)))*cross(n1, n2);
		//Condi��o 1
		if (dot(point - O, n1) < tol && dot(point - O, cross(t, n1)) < tol)
			return true;
		//Condi��o 2
		if (dot(point - O, n2) < tol && dot(point - O, cross(n2, t)) < tol)
			return true;
		//Condi��es n�o atendidas (nem uma nem outra - ponto fora)
		return false;
	}
}

void Polyhedron::Alloc()
{
	n_sub_bv = 0;
	sub_bv = NULL;

	super_node = 0;
	material = 0;
	node = 0;
	number = 0;
	cs = 0;

	DOFs = new int[db.number_GLs_node];
	type_name = new char[20];//Nome do tipo da particula

	//Rotina para ativar os GLS do n� da particula
	for (int j = 0; j < db.number_GLs_node; j++)
	{
		DOFs[j] = 0;
	}
	//Transla��o
	DOFs[0] = 1;
	DOFs[1] = 1;
	DOFs[2] = 1;
	//Rota��o
	DOFs[3] = 1;
	DOFs[4] = 1;
	DOFs[5] = 1;

	//Variaveis internas
	mass = 0.0;

	I3 = new Matrix(3, 3);
	(*I3)(0, 0) = 1.0;
	(*I3)(1, 1) = 1.0;
	(*I3)(2, 2) = 1.0;

	mJr = new Matrix(3, 3);
	mbr = new Matrix(3);
	MassMatrix = new Matrix(6, 6);
	ResidualVector = new Matrix(6);

	Q0 = new Matrix(3, 3);
	Qip = new Matrix(3, 3);
	x0ip = new Matrix(3);

	Qi = new Matrix(3, 3);
	x0i = new Matrix(3);

	DdT = new double*[6];
	for (int i = 0; i < 6; i++)
		DdT[i] = new double[6];
	dT = new double[6];
	Jr = new double*[3];
	for (int i = 0; i < 3; i++)
		Jr[i] = new double[3];
	br = new double[3];
	Ddfield = new double*[6];
	for (int i = 0; i < 6; i++)
		Ddfield[i] = new double[6];
	dfield = new double[6];

	for (int i = 0; i < 6; i++)
	{
		dfield[i] = 0.0;
		dT[i] = 0.0;
		for (int j = 0; j < 6; j++)
		{
			Ddfield[i][j] = 0.0;
			DdT[i][j] = 0.0;
		}
	}

	mJrlocal = new Matrix(3, 3);
	mbrlocal = new Matrix(3);

	kinetic_energy = 0.0;

	angular_momentum[0] = 0.0;
	angular_momentum[1] = 0.0;
	angular_momentum[2] = 0.0;
	angular_momentum_mag = 0.0;

	angular_momentum_origin[0] = 0.0;
	angular_momentum_origin[1] = 0.0;
	angular_momentum_origin[2] = 0.0;

	linear_momentum[0] = 0.0;
	linear_momentum[1] = 0.0;
	linear_momentum[2] = 0.0;
	linear_momentum_mag = 0.0;
}

void Polyhedron::Free()
{
	delete[] DOFs;
	delete[] type_name;
	delete Q0;
	delete Qip;
	delete x0ip;
	delete Qi;
	delete x0i;
	for (int i = 0; i < 6; i++)
		delete[]DdT[i];
	delete[] DdT;
	delete[] dT;

	for (int i = 0; i < 3; i++)
		delete[]Jr[i];
	delete[]Jr;
	delete[]br;
	for (int i = 0; i < 6; i++)
		delete[]Ddfield[i];
	delete[] Ddfield;
	delete[] dfield;
	delete mJr;
	delete mbr;
	delete mJrlocal;
	delete mbrlocal;
	delete MassMatrix;
	delete ResidualVector;
	FreeMassFragmented();
}

void Polyhedron::UpdateVariables()
{
	//Rotation
	Matrix alpha(3);
	alpha(0, 0) = db.nodes[node - 1]->displacements[3];
	alpha(1, 0) = db.nodes[node - 1]->displacements[4];
	alpha(2, 0) = db.nodes[node - 1]->displacements[5];
	double alpha_escalar = norm(alpha);
	Matrix A = skew(alpha);
	double g = 4.0 / (4.0 + alpha_escalar * alpha_escalar);
	Matrix QdA = *I3 + g * (A + 0.5*((A)*(A)));
	if (db.nodes[node-1]->flag_material_description)
		*Qip = (*Q0) * (*db.nodes[node - 1]->Q) * QdA;
	else
		*Qip = QdA * (*db.nodes[node - 1]->Q) * (*Q0);
	//Position
	(*x0ip)(0, 0) = db.nodes[node - 1]->copy_coordinates[0] + db.nodes[node - 1]->displacements[0];
	(*x0ip)(1, 0) = db.nodes[node - 1]->copy_coordinates[1] + db.nodes[node - 1]->displacements[1];
	(*x0ip)(2, 0) = db.nodes[node - 1]->copy_coordinates[2] + db.nodes[node - 1]->displacements[2];
}

void Polyhedron::Mount()
{
	//Zera matrizes diversas
	for (int i = 0; i < 6; i++)
	{
		dfield[i] = 0.0;
		dT[i] = 0.0;
		for (int j = 0; j < 6; j++)
		{
			Ddfield[i][j] = 0.0;
			DdT[i][j] = 0.0;
		}
	}
	//Valores de variaveis cinematicas - atribui��o
	for (int i = 0; i < 3; i++)
	{
		alphai[i] = db.nodes[node - 1]->copy_coordinates[i + 3];
		ud[i] = db.nodes[node - 1]->displacements[i];
		alphad[i] = db.nodes[node - 1]->displacements[i + 3];
		dui[i] = db.nodes[node - 1]->copy_vel[i];
		ddui[i] = db.nodes[node - 1]->copy_accel[i];
		omegai[i] = db.nodes[node - 1]->copy_vel[i + 3];
		domegai[i] = db.nodes[node - 1]->copy_accel[i + 3];
	}
	MountFieldLoading();
	if (db.psy_coupling_exist == false)
		InertialContributions();
}

//Explicit
void Polyhedron::InitialEvaluations()
{
	//Seta variaveis nodais associadas com o tratamento que essa particula vai empregar
	*db.nodes[node - 1]->Q0 = *Q0;
	db.nodes[node - 1]->flag_material_description = true;
	db.nodes[node - 1]->flag_pseudo_moment = false;

	AllocMassFragmented(db.nodes[node - 1]->n_GL_free, db.nodes[node - 1]->n_GL_fixed);


	//Connectivity (local)
	int countA = 0;
	int countB = 0;
	for (int i = 0; i < 6; i++)
	{
		if (db.nodes[node - 1]->GLs[i] > 0)
		{
			GLA[countA] = i;
			countA++;
		}
			
		else
		{
			GLB[countB] = i;
			countB++;
		}
	}

	//Distribui��o da matriz de massa 6x6 (fragmenta��o)
	for (int i = 0; i < db.nodes[node - 1]->n_GL_free; i++)
		for (int j = 0; j < db.nodes[node - 1]->n_GL_free; j++)
			(*invMAA)(i, j) = (*MassMatrix)(GLA[i], GLA[j]);
	for (int i = 0; i < db.nodes[node - 1]->n_GL_free; i++)
		for (int j = 0; j < db.nodes[node - 1]->n_GL_fixed; j++)
			(*MAB)(i, j) = (*MassMatrix)(GLA[i], GLB[j]);
	for (int i = 0; i < db.nodes[node - 1]->n_GL_fixed; i++)
		for (int j = 0; j < db.nodes[node - 1]->n_GL_fixed; j++)
			(*MBB)(i, j) = (*MassMatrix)(GLB[i], GLB[j]);
	//Inversa da matriz MAA
	*invMAA = invert(*invMAA);
}

void Polyhedron::EvaluateExplicit()
{
	//Node pointer
	Node* nd = db.nodes[node - 1];
	Matrix mBr = skew(*mbrlocal);
	//Auxiliary variables
	Matrix A = skew(*nd->alpha);		
	double g = 4.0 / (4.0 + norm(*nd->alpha)*norm(*nd->alpha));
	Matrix Qdelta = *I3 + g * (A + 0.5*(A)*(A));
	Matrix Qp = (*nd->Q0) * (*nd->Q) * Qdelta;
	Matrix Qpt = transp(Qp);
	//Residual vector contributions
	Matrix force_i = - 1.0 * mass * skew(*nd->omega) * (cross(*nd->omega, *mbrlocal));
	Matrix moment_i = - 1.0 * skew(*nd->omega) * ((*mJrlocal) * (*nd->omega));
	
	//Peso pr�prio
	Matrix force_w;
	Matrix moment_w;

	//Variaveis para calculo de steps
	double l_factor;
	if (db.environment_exist == true)
	{
		//Se existe campo gravitacional
		if (db.environment->g_exist == true)
		{
			l_factor = db.environment->bool_g.GetLinearFactorAtCurrentTime();
			force_w = l_factor * mass * Qpt * db.environment->G;
			moment_w =  l_factor * mass * cross(*mbrlocal ,Qpt * db.environment->G);
		}
	}
	
	Matrix force = force_i + force_w;
	Matrix moment = moment_i + moment_w;

	//Residual vector evaluation
	for (int i = 0; i < 3; i++)
	{
		(*ResidualVector)(i, 0) = force(i, 0);
		(*ResidualVector)(i + 3, 0) = moment(i, 0);
	}

	//Espalhamento - Residual vector
	//Variaveis temporarias para salvar a indexa��o global dos graus de liberdade a serem setados na matriz de rigidez global
	int GL_global = 0;
	double anterior = 0;
	for (int i = 0; i < 6; i++)
	{
		GL_global = db.nodes[node - 1]->GLs[i];
		//Caso o grau de liberdade seja livre:
		if (GL_global > 0)
		{
			anterior = db.global_P_A(GL_global - 1, 0);
			db.global_P_A(GL_global - 1, 0) = anterior + (*ResidualVector)(i, 0);
			db.global_I_A(GL_global - 1, 0) = anterior + (*ResidualVector)(i, 0);
		}
		else
		{
			anterior = db.global_P_B(-GL_global - 1, 0);
			db.global_P_B(-GL_global - 1, 0) = anterior + (*ResidualVector)(i, 0);
		}
	}

	//Post-processing
	kinetic_energy = 0.5 * mass * dot(*nd->du, *nd->du) + 0.5 * dot(*nd->omega, (*mJrlocal)*(*nd->omega)) + mass * dot(cross(*nd->omega, *mbrlocal), *nd->du);
	Matrix angular_mom = ((*mJrlocal) + mass * mBr * mBr)*(*nd->omega);
	angular_momentum[0] = angular_mom(0, 0);
	angular_momentum[1] = angular_mom(1, 0);
	angular_momentum[2] = angular_mom(2, 0);
	angular_momentum_mag = norm(angular_mom);
}

void Polyhedron::EvaluateAccelerations()
{
	Matrix residualA(db.nodes[node - 1]->n_GL_free);
	Matrix residualB(db.nodes[node - 1]->n_GL_fixed);

	Matrix accelA(db.nodes[node - 1]->n_GL_free);
	Matrix accelB(db.nodes[node - 1]->n_GL_fixed);
	Matrix accel(6);

	//Tomando as informa��es do vetor global e copiando em residualA:
	int GL_global = 0;
	for (int i = 0; i < db.nodes[node - 1]->n_GL_free; i++)
	{
		GL_global = db.nodes[node-1]->GLs[GLA[i]];
		residualA(i, 0) = db.global_P_A(GL_global - 1, 0);
	}
	//Tomando as informa��es de accelB:
	for (int i = 0; i < db.nodes[node - 1]->n_GL_fixed; i++)
		accelB(i, 0) = db.nodes[node - 1]->accel[GLB[i]];
	
	//Calculando accelA:
	if (db.nodes[node - 1]->n_GL_free != 0 && db.nodes[node - 1]->n_GL_fixed != 0)
		accelA = (*invMAA)*((residualA) - (*MAB) * accelB);
	else
	{
		if (db.nodes[node - 1]->n_GL_fixed == 0)
			accelA = (*invMAA)*(residualA);
	}
	//Calculando residualB:
	if (db.nodes[node - 1]->n_GL_free != 0 && db.nodes[node - 1]->n_GL_fixed != 0)
		residualB = transp(*MAB) * accelA + (*MBB) * accelB;
	else
	{
		if (db.nodes[node - 1]->n_GL_free == 0)
			residualB = (*MBB) * accelB;
	}
	
	//Calculando o vetor accel(6):
	for (int i = 0; i < db.nodes[node - 1]->n_GL_free; i++)
		accel(GLA[i], 0) = accelA(i, 0);
	for (int i = 0; i < db.nodes[node - 1]->n_GL_fixed; i++)
		accel(GLB[i], 0) = accelB(i, 0);
	//Salvando accelera��es nos n�s
	for (int k = 0; k < 3; k++)
	{
		GL_global = db.nodes[node - 1]->GLs[k];
		if (GL_global > 0)
			(*db.nodes[node - 1]->ddu)(k, 0) = accel(k, 0);
	}

	for (int k = 0; k < 3; k++)
	{
		GL_global = db.nodes[node - 1]->GLs[k + 3];
		if (GL_global > 0)
			(*db.nodes[node - 1]->domega)(k, 0) = accel(k + 3, 0);
	}
	//Salvando rea��es vinculares em P_B
	for (int i = 0; i < db.nodes[node - 1]->n_GL_fixed; i++)
	{
		GL_global = db.nodes[node - 1]->GLs[GLB[i]];
		double anterior = db.global_P_B(-GL_global - 1, 0);
		db.global_P_B(-GL_global - 1, 0) = anterior + residualB(i, 0);
	}
}

void Polyhedron::MountGlobal()
{
	//Variaveis temporarias para salvar a indexa��o global dos graus de liberdade a serem setados na matriz de rigidez global
	int GL_global_1 = 0;
	int GL_global_2 = 0;
	double anterior = 0;
	for (int i = 0; i < 6; i++)
	{
		//////////////MONTAGEM DO VETOR DE ESFOR�OS DESBALANCEADOS//////////////////
		//Toma do vetor de GL globais, a indexa��o de cada grau de liberdade global
		GL_global_1 = db.nodes[node - 1]->GLs[i];

		//Caso o grau de liberdade seja livre:
		if (GL_global_1 > 0)
		{
			anterior = db.global_P_A(GL_global_1 - 1, 0);
			db.global_P_A(GL_global_1 - 1, 0) = anterior + dT[i] - dfield[i];
			anterior = db.global_I_A(GL_global_1 - 1, 0);
			db.global_I_A(GL_global_1 - 1, 0) = anterior + dT[i] - dfield[i];
		}
		else
		{
			anterior = db.global_P_B(-GL_global_1 - 1, 0);
			db.global_P_B(-GL_global_1 - 1, 0) = anterior + dT[i] - dfield[i];
		}
		for (int j = 0; j < 6; j++)
		{
			//////////////////////MONTAGEM DA MATRIZ DE RIGIDEZ/////////////////////////
			//Toma do vetor de GL globais, a indexa��o de cada grau de liberdade global
			GL_global_2 = db.nodes[node - 1]->GLs[j];
			//Caso os graus de liberdade sejam ambos livres (Matriz Kaa)
			if (GL_global_1 > 0 && GL_global_2 > 0)
				db.global_stiffness_AA.setValue(GL_global_1 - 1, GL_global_2 - 1, (DdT[i][j] - Ddfield[i][j]));
			//Caso os graus de liberdade sejam ambos fixos (Matriz Kbb)
			if (GL_global_1 < 0 && GL_global_2 < 0)
				db.global_stiffness_BB.setValue(-GL_global_1 - 1, -GL_global_2 - 1, (DdT[i][j] - Ddfield[i][j]));
			//Caso os graus de liberdade sejam livre e fixo (Matriz Kab)
			if (GL_global_1 > 0 && GL_global_2 < 0)
				db.global_stiffness_AB.setValue(GL_global_1 - 1, -GL_global_2 - 1, (DdT[i][j] - Ddfield[i][j]));
			//Caso os graus de liberdade sejam fixo e livre (Matriz Kba)
			if (GL_global_1 < 0 && GL_global_2 > 0)
				db.global_stiffness_BA.setValue(-GL_global_1 - 1, GL_global_2 - 1, (DdT[i][j] - Ddfield[i][j]));
		}
	}
}

void Polyhedron::MountFieldLoading()
{
	//Ponteiro para o valor da massa
	double* m = &mass;
	double v[1000];													//Variavel necessaria para o AceGen		
	//Variaveis para calculo de steps
	double l_factor;
	bool in_gcs_domain = true;
	if (db.environment_exist == true)
	{
		//Se existe campo gravitacional
		if (db.environment->g_exist == true)
		{
			//if (db.gcs_exist)
			//{
			//	//If the user has set boundings for contacts, consider here the same boundings for the gravitational field
			//	if (db.gcs->user_set_boundings)
			//	{
			//		double center[3];
			//		//teest is done considering the last converged position
			//		center[0] = db.nodes[node - 1]->copy_coordinates[0];
			//		center[1] = db.nodes[node - 1]->copy_coordinates[1];
			//		center[2] = db.nodes[node - 1]->copy_coordinates[2];


			//		if (center[0] < db.gcs->min_x || center[0] > db.gcs->max_x ||
			//			center[1] < db.gcs->min_y || center[1] > db.gcs->max_y ||
			//			center[2] < db.gcs->min_z || center[2] > db.gcs->max_z)
			//		{
			//			in_gcs_domain = false;
			//		}

			//	}

			//}
			//includes the weight only it the particle lies inside the domain of interest, defined by gcs
			if (in_gcs_domain)
			{
				l_factor = db.environment->bool_g.GetLinearFactorAtCurrentTime();
				double g[3];					//Vetor gravidade
				g[0] = l_factor * db.environment->G(0, 0);
				g[1] = l_factor * db.environment->G(1, 0);
				g[2] = l_factor * db.environment->G(2, 0);

				//AceGen
				int i01; int i02;
				v[257] = 0.5e0*alphad[2];
				v[255] = Power(alphad[2], 2);
				v[254] = alphad[0] * v[257];
				v[253] = 0.5e0*alphad[1];
				v[252] = 2e0*alphad[2];
				v[251] = Power(alphad[1], 2);
				v[250] = alphad[0] * v[253];
				v[249] = 2e0*alphad[1];
				v[248] = Power(alphad[0], 2);
				v[247] = 2e0*alphad[0];
				v[246] = Power(alphai[2], 2);
				v[245] = 0.5e0*alphai[0] * alphai[2];
				v[244] = 0.5e0*alphai[1];
				v[243] = Power(alphai[1], 2);
				v[261] = v[243] + v[246];
				v[242] = alphai[0] * v[244];
				v[241] = Power(alphai[0], 2);
				v[97] = alphai[2] * v[244];
				v[172] = -v[248] - v[251];
				v[150] = alphad[2] + v[250];
				v[140] = -alphad[2] + v[250];
				v[81] = alphad[2] * v[253];
				v[167] = alphad[0] + v[81];
				v[159] = -alphad[0] + v[81];
				v[163] = -alphad[1] + v[254];
				v[146] = alphad[1] + v[254];
				v[154] = -v[248] - v[255];
				v[135] = -v[251] - v[255];
				v[130] = 4e0 + v[248] + v[251] + v[255];
				v[132] = 1e0 / Power(v[130], 2);
				v[256] = -4e0*v[132];
				v[134] = v[252] * v[256];
				v[260] = 0.5e0*v[134];
				v[176] = v[172] * v[260];
				v[133] = v[249] * v[256];
				v[259] = 0.5e0*v[133];
				v[156] = v[154] * v[259];
				v[131] = v[247] * v[256];
				v[258] = 0.5e0*v[131];
				v[180] = v[131] * v[253];
				v[177] = -(v[131] * v[257]);
				v[136] = v[135] * v[258];
				v[68] = 4e0 / v[130];
				v[181] = -0.5e0*v[68];
				v[183] = v[181] - alphad[0] * v[258];
				v[174] = v[181] * v[249];
				v[175] = v[174] + v[172] * v[259];
				v[171] = v[181] * v[247];
				v[173] = v[171] + v[172] * v[258];
				v[168] = v[131] * v[167] + v[68];
				v[165] = v[133] * v[163] - v[68];
				v[160] = v[131] * v[159] - v[68];
				v[157] = v[181] * v[252];
				v[158] = v[157] + v[154] * v[260];
				v[155] = v[171] + v[154] * v[258];
				v[153] = v[134] * v[150] + v[68];
				v[148] = v[133] * v[146] + v[68];
				v[145] = alphad[2] * v[181];
				v[169] = -v[145] + v[133] * v[167];
				v[164] = -v[145] + v[131] * v[163];
				v[161] = -v[145] + v[133] * v[159];
				v[147] = -v[145] + v[131] * v[146];
				v[144] = v[134] * v[140] - v[68];
				v[142] = alphad[0] * v[181];
				v[166] = -v[142] + v[134] * v[163];
				v[152] = -v[142] + v[133] * v[150];
				v[149] = -v[142] + v[134] * v[146];
				v[143] = v[133] * v[140] - v[142];
				v[139] = -(alphad[1] * v[181]);
				v[170] = v[139] + v[134] * v[167];
				v[162] = v[139] + v[134] * v[159];
				v[151] = v[139] + v[131] * v[150];
				v[141] = v[139] + v[131] * v[140];
				v[138] = v[157] + v[135] * v[260];
				v[137] = v[174] + v[135] * v[259];
				v[71] = 1e0 - v[135] * v[181];
				v[72] = v[140] * v[68];
				v[73] = v[146] * v[68];
				v[75] = v[150] * v[68];
				v[77] = 1e0 - v[154] * v[181];
				v[78] = v[159] * v[68];
				v[80] = v[163] * v[68];
				v[82] = v[167] * v[68];
				v[83] = 1e0 - v[172] * v[181];
				v[84] = 4e0 / (4e0 + v[241] + v[261]);
				v[262] = -0.5e0*v[84];
				v[87] = 1e0 + v[261] * v[262];
				v[88] = (-alphai[2] + v[242])*v[84];
				v[89] = (alphai[1] + v[245])*v[84];
				v[91] = (alphai[2] + v[242])*v[84];
				v[93] = 1e0 + (v[241] + v[246])*v[262];
				v[94] = v[84] * (-alphai[0] + v[97]);
				v[96] = (-alphai[1] + v[245])*v[84];
				v[98] = v[84] * (alphai[0] + v[97]);
				v[99] = 1e0 + (v[241] + v[243])*v[262];
				v[222] = br[0] * (v[166] * v[87] + v[170] * v[91] + v[176] * v[96]) + br[1] * (v[166] * v[88] + v[170] * v[93] + v[176] * v[98]
					) + br[2] * (v[166] * v[89] + v[170] * v[94] + v[176] * v[99]);
				v[218] = br[0] * (v[165] * v[87] + v[169] * v[91] + v[175] * v[96]) + br[1] * (v[165] * v[88] + v[169] * v[93] + v[175] * v[98]
					) + br[2] * (v[165] * v[89] + v[169] * v[94] + v[175] * v[99]);
				v[214] = br[0] * (v[164] * v[87] + v[168] * v[91] + v[173] * v[96]) + br[1] * (v[164] * v[88] + v[168] * v[93] + v[173] * v[98]
					) + br[2] * (v[164] * v[89] + v[168] * v[94] + v[173] * v[99]);
				v[207] = br[0] * (v[153] * v[87] + v[158] * v[91] + v[162] * v[96]) + br[1] * (v[153] * v[88] + v[158] * v[93] + v[162] * v[98]
					) + br[2] * (v[153] * v[89] + v[158] * v[94] + v[162] * v[99]);
				v[203] = br[0] * (v[152] * v[87] + v[156] * v[91] + v[161] * v[96]) + br[1] * (v[152] * v[88] + v[156] * v[93] + v[161] * v[98]
					) + br[2] * (v[152] * v[89] + v[156] * v[94] + v[161] * v[99]);
				v[224] = (*m)*(g[2] * v[203] - g[1] * v[218]);
				v[199] = br[0] * (v[151] * v[87] + v[155] * v[91] + v[160] * v[96]) + br[1] * (v[151] * v[88] + v[155] * v[93] + v[160] * v[98]
					) + br[2] * (v[151] * v[89] + v[155] * v[94] + v[160] * v[99]);
				v[223] = (*m)*(g[2] * v[199] - g[1] * v[214]);
				v[195] = br[0] * (v[138] * v[87] + v[144] * v[91] + v[149] * v[96]) + br[1] * (v[138] * v[88] + v[144] * v[93] + v[149] * v[98]
					) + br[2] * (v[138] * v[89] + v[144] * v[94] + v[149] * v[99]);
				v[191] = br[0] * (v[137] * v[87] + v[143] * v[91] + v[148] * v[96]) + br[1] * (v[137] * v[88] + v[143] * v[93] + v[148] * v[98]
					) + br[2] * (v[137] * v[89] + v[143] * v[94] + v[148] * v[99]);
				v[227] = (*m)*(-(g[2] * v[191]) + g[0] * v[218]);
				v[209] = (*m)*(g[1] * v[191] - g[0] * v[203]);
				v[187] = br[0] * (v[136] * v[87] + v[141] * v[91] + v[147] * v[96]) + br[1] * (v[136] * v[88] + v[141] * v[93] + v[147] * v[98]
					) + br[2] * (v[136] * v[89] + v[141] * v[94] + v[147] * v[99]);
				v[226] = (*m)*(-(g[2] * v[187]) + g[0] * v[214]);
				v[208] = (*m)*(g[1] * v[187] - g[0] * v[199]);
				v[112] = br[0] * (v[71] * v[87] + v[72] * v[91] + v[73] * v[96]) + br[1] * (v[71] * v[88] + v[72] * v[93] + v[73] * v[98])
					+ br[2] * (v[71] * v[89] + v[72] * v[94] + v[73] * v[99]);
				v[113] = br[0] * (v[75] * v[87] + v[77] * v[91] + v[78] * v[96]) + br[1] * (v[75] * v[88] + v[77] * v[93] + v[78] * v[98])
					+ br[2] * (v[75] * v[89] + v[77] * v[94] + v[78] * v[99]);
				v[119] = (*m)*(g[1] * v[112] - g[0] * v[113]);
				v[232] = -(v[119] * v[180]);
				v[114] = br[0] * (v[80] * v[87] + v[82] * v[91] + v[83] * v[96]) + br[1] * (v[80] * v[88] + v[82] * v[93] + v[83] * v[98])
					+ br[2] * (v[80] * v[89] + v[82] * v[94] + v[83] * v[99]);
				v[121] = (*m)*(g[2] * v[113] - g[1] * v[114]);
				v[236] = -(alphad[2] * v[121] * v[259]);
				v[120] = (*m)*(-(g[2] * v[112]) + g[0] * v[114]);
				v[235] = -(v[120] * v[177]);
				dfield[0] = g[0] * (*m);
				dfield[1] = g[1] * (*m);
				dfield[2] = g[2] * (*m);
				dfield[3] = -(v[119] * v[139]) - v[120] * v[145] + v[121] * v[68];
				dfield[4] = -(v[119] * v[142]) + v[121] * v[145] + v[120] * v[68];
				dfield[5] = v[121] * v[139] + v[120] * v[142] + v[119] * v[68];
				Ddfield[0][0] = 0e0;
				Ddfield[0][1] = 0e0;
				Ddfield[0][2] = 0e0;
				Ddfield[0][3] = 0e0;
				Ddfield[0][4] = 0e0;
				Ddfield[0][5] = 0e0;
				Ddfield[1][1] = 0e0;
				Ddfield[1][2] = 0e0;
				Ddfield[1][3] = 0e0;
				Ddfield[1][4] = 0e0;
				Ddfield[1][5] = 0e0;
				Ddfield[2][2] = 0e0;
				Ddfield[2][3] = 0e0;
				Ddfield[2][4] = 0e0;
				Ddfield[2][5] = 0e0;
				Ddfield[3][3] = v[121] * v[131] - v[139] * v[208] - v[145] * v[226] + v[232] + v[235] + v[223] * v[68];
				Ddfield[3][4] = v[120] * v[131] + v[121] * v[177] - v[119] * v[183] - v[142] * v[208] + v[145] * v[223] + v[226] * v[68];
				Ddfield[3][5] = v[119] * v[131] + v[121] * v[180] + v[120] * v[183] + v[139] * v[223] + v[142] * v[226] + v[208] * v[68];
				Ddfield[4][4] = v[120] * v[133] - v[142] * v[209] + v[145] * v[224] - v[232] + v[236] + v[227] * v[68];
				Ddfield[4][5] = v[119] * v[133] - v[120] * v[180] + v[139] * v[224] + v[142] * v[227] + v[121] * (-v[181]
					+ alphad[1] * v[259]) + v[209] * v[68];
				Ddfield[5][5] = v[119] * v[134] - v[235] - v[236] + (*m)*(v[142] * (-(g[2] * v[195]) + g[0] * v[222]) + v[139] *
					(g[2] * v[207] - g[1] * v[222]) + (g[1] * v[195] - g[0] * v[207])*v[68]);
				for (i01 = 1; i01 < 6; i01++) {
					for (i02 = 0; i02 < i01; i02++) {
						Ddfield[i01][i02] = Ddfield[i02][i01];
					}
				};
			}

		}
	}

	//Outras contribui��es de esfor�os de campo, colocar aqui...
}

//Calcula contribui��es de inercia
void Polyhedron::InertialContributions()
{
	double value = 0.0;
	double* a1;
	double* a2;
	double* a3;
	double* a4;
	double* a5;
	double* a6;
	double xpolei[3];
	xpolei[0] = db.nodes[node - 1]->copy_coordinates[0];
	xpolei[1] = db.nodes[node - 1]->copy_coordinates[1];
	xpolei[2] = db.nodes[node - 1]->copy_coordinates[2];

	if (typeid(*db.solution[db.current_solution_number - 1]) == typeid(Dynamic))
	{
		Dynamic* ptr_sol = static_cast<Dynamic*>(db.solution[db.current_solution_number - 1]);
		//Ponteiros para parametros do Newmark
		a1 = &ptr_sol->a1;
		a2 = &ptr_sol->a2;
		a3 = &ptr_sol->a3;
		a4 = &ptr_sol->a4;
		a5 = &ptr_sol->a5;
		a6 = &ptr_sol->a6;
	}
	else
	{
		a1 = &value;
		a2 = &value;
		a3 = &value;
		a4 = &value;
		a5 = &value;
		a6 = &value;
	}

	//Ponteiro para o valor da massa
	double* m = &mass;
	double v[1000];													//Variavel necessaria para o AceGen		
	//AceGen
	v[134] = Power(alphai[0], 2);
	v[132] = 0.5e0*alphai[0] * alphai[1];
	v[127] = Power(alphai[1], 2);
	v[139] = 0.5e0*alphai[1] * alphai[2];
	v[137] = 0.5e0*alphai[0] * alphai[2];
	v[128] = Power(alphai[2], 2);
	v[523] = v[127] + v[128];
	v[214] = 0.5e0*alphad[0];
	v[212] = 2e0*alphad[0];
	v[118] = Power(alphad[0], 2);
	v[215] = 2e0*alphad[1];
	v[213] = 0.5e0*alphad[1];
	v[116] = 0.5e0*alphad[0] * alphad[1];
	v[111] = Power(alphad[1], 2);
	v[260] = -v[111] - v[118];
	v[238] = alphad[2] + v[116];
	v[228] = -alphad[2] + v[116];
	v[217] = 2e0*alphad[2];
	v[216] = 0.5e0*alphad[2];
	v[123] = 0.5e0*alphad[1] * alphad[2];
	v[255] = alphad[0] + v[123];
	v[247] = -alphad[0] + v[123];
	v[121] = 0.5e0*alphad[0] * alphad[2];
	v[251] = -alphad[1] + v[121];
	v[234] = alphad[1] + v[121];
	v[112] = Power(alphad[2], 2);
	v[242] = -v[112] - v[118];
	v[223] = -v[111] - v[112];
	v[218] = 4e0 + v[111] + v[112] + v[118];
	v[220] = 1e0 / Power(v[218], 2);
	v[222] = -4e0*v[217] * v[220];
	v[264] = 0.5e0*v[222] * v[260];
	v[221] = -4e0*v[215] * v[220];
	v[360] = -(v[214] * v[221]);
	v[356] = -(v[216] * v[221]);
	v[244] = 0.5e0*v[221] * v[242];
	v[219] = -4e0*v[212] * v[220];
	v[355] = -(v[216] * v[219]);
	v[224] = 0.5e0*v[219] * v[223];
	v[191] = (*a1)*alphad[0] - (*a3)*domegai[0] - (*a2)*omegai[0];
	v[185] = (*a4)*alphad[0] + (*a6)*domegai[0] + (*a5)*omegai[0];
	v[192] = (*a1)*alphad[1] - (*a3)*domegai[1] - (*a2)*omegai[1];
	v[186] = (*a4)*alphad[1] + (*a6)*domegai[1] + (*a5)*omegai[1];
	v[193] = (*a1)*alphad[2] - (*a3)*domegai[2] - (*a2)*omegai[2];
	v[187] = (*a4)*alphad[2] + (*a6)*domegai[2] + (*a5)*omegai[2];
	v[211] = (*a1)*(*m);
	v[110] = 4e0 / v[218];
	v[358] = -0.5e0*v[110];
	v[361] = v[213] * v[221] - v[358];
	v[359] = -(v[214] * v[219]) + v[358];
	v[357] = -(v[216] * v[222]) + v[358];
	v[262] = -0.5e0*v[110] * v[215];
	v[263] = 0.5e0*v[221] * v[260] + v[262];
	v[259] = -0.5e0*v[110] * v[212];
	v[261] = v[259] + 0.5e0*v[219] * v[260];
	v[256] = v[110] + v[219] * v[255];
	v[253] = -v[110] + v[221] * v[251];
	v[248] = -v[110] + v[219] * v[247];
	v[245] = -0.5e0*v[110] * v[217];
	v[246] = 0.5e0*v[222] * v[242] + v[245];
	v[243] = 0.5e0*v[219] * v[242] + v[259];
	v[241] = v[110] + v[222] * v[238];
	v[236] = v[110] + v[221] * v[234];
	v[233] = -(v[110] * v[216]);
	v[257] = -v[233] + v[221] * v[255];
	v[252] = -v[233] + v[219] * v[251];
	v[249] = -v[233] + v[221] * v[247];
	v[235] = -v[233] + v[219] * v[234];
	v[232] = -v[110] + v[222] * v[228];
	v[230] = -(v[110] * v[214]);
	v[254] = -v[230] + v[222] * v[251];
	v[240] = -v[230] + v[221] * v[238];
	v[237] = -v[230] + v[222] * v[234];
	v[231] = v[221] * v[228] - v[230];
	v[227] = v[110] * v[213];
	v[258] = v[227] + v[222] * v[255];
	v[250] = v[227] + v[222] * v[247];
	v[239] = v[227] + v[219] * v[238];
	v[229] = v[227] + v[219] * v[228];
	v[226] = 0.5e0*v[222] * v[223] + v[245];
	v[225] = 0.5e0*v[221] * v[223] + v[262];
	v[113] = 1e0 + 0.5e0*v[110] * v[223];
	v[319] = (*a1)*v[113] + v[191] * v[224] + v[192] * v[229] + v[193] * v[235];
	v[301] = (*a4)*v[113] + v[185] * v[224] + v[186] * v[229] + v[187] * v[235];
	v[114] = v[110] * v[228];
	v[320] = (*a1)*v[114] + v[191] * v[225] + v[192] * v[231] + v[193] * v[236];
	v[302] = (*a4)*v[114] + v[185] * v[225] + v[186] * v[231] + v[187] * v[236];
	v[115] = v[110] * v[234];
	v[321] = (*a1)*v[115] + v[191] * v[226] + v[192] * v[232] + v[193] * v[237];
	v[303] = (*a4)*v[115] + v[185] * v[226] + v[186] * v[232] + v[187] * v[237];
	v[117] = v[110] * v[238];
	v[322] = (*a1)*v[117] + v[191] * v[239] + v[192] * v[243] + v[193] * v[248];
	v[304] = (*a4)*v[117] + v[185] * v[239] + v[186] * v[243] + v[187] * v[248];
	v[119] = 1e0 + 0.5e0*v[110] * v[242];
	v[323] = (*a1)*v[119] + v[191] * v[240] + v[192] * v[244] + v[193] * v[249];
	v[305] = (*a4)*v[119] + v[185] * v[240] + v[186] * v[244] + v[187] * v[249];
	v[120] = v[110] * v[247];
	v[324] = (*a1)*v[120] + v[191] * v[241] + v[192] * v[246] + v[193] * v[250];
	v[306] = (*a4)*v[120] + v[185] * v[241] + v[186] * v[246] + v[187] * v[250];
	v[122] = v[110] * v[251];
	v[325] = (*a1)*v[122] + v[191] * v[252] + v[192] * v[256] + v[193] * v[261];
	v[310] = (*a4)*v[122] + v[185] * v[252] + v[186] * v[256] + v[187] * v[261];
	v[124] = v[110] * v[255];
	v[326] = (*a1)*v[124] + v[191] * v[253] + v[192] * v[257] + v[193] * v[263];
	v[311] = (*a4)*v[124] + v[185] * v[253] + v[186] * v[257] + v[187] * v[263];
	v[125] = 1e0 + 0.5e0*v[110] * v[260];
	v[327] = (*a1)*v[125] + v[191] * v[254] + v[192] * v[258] + v[193] * v[264];
	v[312] = (*a4)*v[125] + v[185] * v[254] + v[186] * v[258] + v[187] * v[264];
	v[126] = 4e0 / (4e0 + v[134] + v[523]);
	v[129] = 1e0 - 0.5e0*v[126] * v[523];
	v[130] = v[126] * (-alphai[2] + v[132]);
	v[131] = v[126] * (alphai[1] + v[137]);
	v[133] = v[126] * (alphai[2] + v[132]);
	v[135] = 1e0 - 0.5e0*v[126] * (v[128] + v[134]);
	v[136] = v[126] * (-alphai[0] + v[139]);
	v[138] = v[126] * (-alphai[1] + v[137]);
	v[285] = v[129] * v[254] + v[133] * v[258] + v[138] * v[264];
	v[284] = v[129] * v[253] + v[133] * v[257] + v[138] * v[263];
	v[283] = v[129] * v[252] + v[133] * v[256] + v[138] * v[261];
	v[276] = v[129] * v[241] + v[133] * v[246] + v[138] * v[250];
	v[275] = v[129] * v[240] + v[133] * v[244] + v[138] * v[249];
	v[274] = v[129] * v[239] + v[133] * v[243] + v[138] * v[248];
	v[267] = v[129] * v[226] + v[133] * v[232] + v[138] * v[237];
	v[266] = v[129] * v[225] + v[133] * v[231] + v[138] * v[236];
	v[265] = v[129] * v[224] + v[133] * v[229] + v[138] * v[235];
	v[140] = v[126] * (alphai[0] + v[139]);
	v[288] = v[130] * v[254] + v[135] * v[258] + v[140] * v[264];
	v[287] = v[130] * v[253] + v[135] * v[257] + v[140] * v[263];
	v[286] = v[130] * v[252] + v[135] * v[256] + v[140] * v[261];
	v[279] = v[130] * v[241] + v[135] * v[246] + v[140] * v[250];
	v[278] = v[130] * v[240] + v[135] * v[244] + v[140] * v[249];
	v[277] = v[130] * v[239] + v[135] * v[243] + v[140] * v[248];
	v[270] = v[130] * v[226] + v[135] * v[232] + v[140] * v[237];
	v[269] = v[130] * v[225] + v[135] * v[231] + v[140] * v[236];
	v[268] = v[130] * v[224] + v[135] * v[229] + v[140] * v[235];
	v[141] = 1e0 - 0.5e0*v[126] * (v[127] + v[134]);
	v[291] = v[131] * v[254] + v[136] * v[258] + v[141] * v[264];
	v[388] = Jr[0][0] * v[285] + Jr[1][0] * v[288] + Jr[2][0] * v[291];
	v[385] = Jr[0][1] * v[285] + Jr[1][1] * v[288] + Jr[2][1] * v[291];
	v[382] = Jr[0][2] * v[285] + Jr[1][2] * v[288] + Jr[2][2] * v[291];
	v[300] = br[0] * v[285] + br[1] * v[288] + br[2] * v[291];
	v[290] = v[131] * v[253] + v[136] * v[257] + v[141] * v[263];
	v[387] = Jr[0][0] * v[284] + Jr[1][0] * v[287] + Jr[2][0] * v[290];
	v[384] = Jr[0][1] * v[284] + Jr[1][1] * v[287] + Jr[2][1] * v[290];
	v[381] = Jr[0][2] * v[284] + Jr[1][2] * v[287] + Jr[2][2] * v[290];
	v[299] = br[0] * v[284] + br[1] * v[287] + br[2] * v[290];
	v[289] = v[131] * v[252] + v[136] * v[256] + v[141] * v[261];
	v[386] = Jr[0][0] * v[283] + Jr[1][0] * v[286] + Jr[2][0] * v[289];
	v[383] = Jr[0][1] * v[283] + Jr[1][1] * v[286] + Jr[2][1] * v[289];
	v[380] = Jr[0][2] * v[283] + Jr[1][2] * v[286] + Jr[2][2] * v[289];
	v[298] = br[0] * v[283] + br[1] * v[286] + br[2] * v[289];
	v[282] = v[131] * v[241] + v[136] * v[246] + v[141] * v[250];
	v[379] = Jr[0][0] * v[276] + Jr[1][0] * v[279] + Jr[2][0] * v[282];
	v[376] = Jr[0][1] * v[276] + Jr[1][1] * v[279] + Jr[2][1] * v[282];
	v[373] = Jr[0][2] * v[276] + Jr[1][2] * v[279] + Jr[2][2] * v[282];
	v[297] = br[0] * v[276] + br[1] * v[279] + br[2] * v[282];
	v[281] = v[131] * v[240] + v[136] * v[244] + v[141] * v[249];
	v[378] = Jr[0][0] * v[275] + Jr[1][0] * v[278] + Jr[2][0] * v[281];
	v[375] = Jr[0][1] * v[275] + Jr[1][1] * v[278] + Jr[2][1] * v[281];
	v[372] = Jr[0][2] * v[275] + Jr[1][2] * v[278] + Jr[2][2] * v[281];
	v[296] = br[0] * v[275] + br[1] * v[278] + br[2] * v[281];
	v[280] = v[131] * v[239] + v[136] * v[243] + v[141] * v[248];
	v[377] = Jr[0][0] * v[274] + Jr[1][0] * v[277] + Jr[2][0] * v[280];
	v[374] = Jr[0][1] * v[274] + Jr[1][1] * v[277] + Jr[2][1] * v[280];
	v[371] = Jr[0][2] * v[274] + Jr[1][2] * v[277] + Jr[2][2] * v[280];
	v[295] = br[0] * v[274] + br[1] * v[277] + br[2] * v[280];
	v[273] = v[131] * v[226] + v[136] * v[232] + v[141] * v[237];
	v[370] = Jr[0][0] * v[267] + Jr[1][0] * v[270] + Jr[2][0] * v[273];
	v[367] = Jr[0][1] * v[267] + Jr[1][1] * v[270] + Jr[2][1] * v[273];
	v[364] = Jr[0][2] * v[267] + Jr[1][2] * v[270] + Jr[2][2] * v[273];
	v[294] = br[0] * v[267] + br[1] * v[270] + br[2] * v[273];
	v[272] = v[131] * v[225] + v[136] * v[231] + v[141] * v[236];
	v[369] = Jr[0][0] * v[266] + Jr[1][0] * v[269] + Jr[2][0] * v[272];
	v[366] = Jr[0][1] * v[266] + Jr[1][1] * v[269] + Jr[2][1] * v[272];
	v[363] = Jr[0][2] * v[266] + Jr[1][2] * v[269] + Jr[2][2] * v[272];
	v[293] = br[0] * v[266] + br[1] * v[269] + br[2] * v[272];
	v[271] = v[131] * v[224] + v[136] * v[229] + v[141] * v[235];
	v[368] = Jr[0][0] * v[265] + Jr[1][0] * v[268] + Jr[2][0] * v[271];
	v[365] = Jr[0][1] * v[265] + Jr[1][1] * v[268] + Jr[2][1] * v[271];
	v[362] = Jr[0][2] * v[265] + Jr[1][2] * v[268] + Jr[2][2] * v[271];
	v[292] = br[0] * v[265] + br[1] * v[268] + br[2] * v[271];
	v[142] = v[113] * v[129] + v[114] * v[133] + v[115] * v[138];
	v[143] = v[113] * v[130] + v[114] * v[135] + v[115] * v[140];
	v[144] = v[113] * v[131] + v[114] * v[136] + v[115] * v[141];
	v[160] = Jr[0][2] * v[142] + Jr[1][2] * v[143] + Jr[2][2] * v[144];
	v[159] = Jr[0][1] * v[142] + Jr[1][1] * v[143] + Jr[2][1] * v[144];
	v[158] = Jr[0][0] * v[142] + Jr[1][0] * v[143] + Jr[2][0] * v[144];
	v[391] = v[158] * v[267] + v[159] * v[270] + v[160] * v[273] + v[144] * v[364] + v[143] * v[367] + v[142] * v[370];
	v[390] = v[158] * v[266] + v[159] * v[269] + v[160] * v[272] + v[144] * v[363] + v[143] * v[366] + v[142] * v[369];
	v[389] = v[158] * v[265] + v[159] * v[268] + v[160] * v[271] + v[144] * v[362] + v[143] * v[365] + v[142] * v[368];
	v[145] = v[117] * v[129] + v[119] * v[133] + v[120] * v[138];
	v[146] = v[117] * v[130] + v[119] * v[135] + v[120] * v[140];
	v[147] = v[117] * v[131] + v[119] * v[136] + v[120] * v[141];
	v[394] = v[158] * v[276] + v[159] * v[279] + v[160] * v[282] + v[147] * v[364] + v[146] * v[367] + v[145] * v[370];
	v[393] = v[158] * v[275] + v[159] * v[278] + v[160] * v[281] + v[147] * v[363] + v[146] * v[366] + v[145] * v[369];
	v[392] = v[158] * v[274] + v[159] * v[277] + v[160] * v[280] + v[147] * v[362] + v[146] * v[365] + v[145] * v[368];
	v[166] = Jr[0][2] * v[145] + Jr[1][2] * v[146] + Jr[2][2] * v[147];
	v[165] = Jr[0][1] * v[145] + Jr[1][1] * v[146] + Jr[2][1] * v[147];
	v[164] = Jr[0][0] * v[145] + Jr[1][0] * v[146] + Jr[2][0] * v[147];
	v[403] = v[164] * v[276] + v[165] * v[279] + v[166] * v[282] + v[147] * v[373] + v[146] * v[376] + v[145] * v[379];
	v[402] = v[164] * v[275] + v[165] * v[278] + v[166] * v[281] + v[147] * v[372] + v[146] * v[375] + v[145] * v[378];
	v[401] = v[164] * v[274] + v[165] * v[277] + v[166] * v[280] + v[147] * v[371] + v[146] * v[374] + v[145] * v[377];
	v[400] = v[164] * v[267] + v[165] * v[270] + v[166] * v[273] + v[144] * v[373] + v[143] * v[376] + v[142] * v[379];
	v[399] = v[164] * v[266] + v[165] * v[269] + v[166] * v[272] + v[144] * v[372] + v[143] * v[375] + v[142] * v[378];
	v[398] = v[164] * v[265] + v[165] * v[268] + v[166] * v[271] + v[144] * v[371] + v[143] * v[374] + v[142] * v[377];
	v[148] = v[122] * v[129] + v[124] * v[133] + v[125] * v[138];
	v[149] = v[122] * v[130] + v[124] * v[135] + v[125] * v[140];
	v[150] = v[122] * v[131] + v[124] * v[136] + v[125] * v[141];
	v[406] = v[164] * v[285] + v[165] * v[288] + v[166] * v[291] + v[150] * v[373] + v[149] * v[376] + v[148] * v[379];
	v[405] = v[164] * v[284] + v[165] * v[287] + v[166] * v[290] + v[150] * v[372] + v[149] * v[375] + v[148] * v[378];
	v[404] = v[164] * v[283] + v[165] * v[286] + v[166] * v[289] + v[150] * v[371] + v[149] * v[374] + v[148] * v[377];
	v[397] = v[158] * v[285] + v[159] * v[288] + v[160] * v[291] + v[150] * v[364] + v[149] * v[367] + v[148] * v[370];
	v[396] = v[158] * v[284] + v[159] * v[287] + v[160] * v[290] + v[150] * v[363] + v[149] * v[366] + v[148] * v[369];
	v[395] = v[158] * v[283] + v[159] * v[286] + v[160] * v[289] + v[150] * v[362] + v[149] * v[365] + v[148] * v[368];
	v[172] = Jr[0][2] * v[148] + Jr[1][2] * v[149] + Jr[2][2] * v[150];
	v[171] = Jr[0][1] * v[148] + Jr[1][1] * v[149] + Jr[2][1] * v[150];
	v[170] = Jr[0][0] * v[148] + Jr[1][0] * v[149] + Jr[2][0] * v[150];
	v[415] = v[170] * v[285] + v[171] * v[288] + v[172] * v[291] + v[150] * v[382] + v[149] * v[385] + v[148] * v[388];
	v[414] = v[170] * v[284] + v[171] * v[287] + v[172] * v[290] + v[150] * v[381] + v[149] * v[384] + v[148] * v[387];
	v[413] = v[170] * v[283] + v[171] * v[286] + v[172] * v[289] + v[150] * v[380] + v[149] * v[383] + v[148] * v[386];
	v[412] = v[170] * v[276] + v[171] * v[279] + v[172] * v[282] + v[147] * v[382] + v[146] * v[385] + v[145] * v[388];
	v[411] = v[170] * v[275] + v[171] * v[278] + v[172] * v[281] + v[147] * v[381] + v[146] * v[384] + v[145] * v[387];
	v[410] = v[170] * v[274] + v[171] * v[277] + v[172] * v[280] + v[147] * v[380] + v[146] * v[383] + v[145] * v[386];
	v[409] = v[170] * v[267] + v[171] * v[270] + v[172] * v[273] + v[144] * v[382] + v[143] * v[385] + v[142] * v[388];
	v[408] = v[170] * v[266] + v[171] * v[269] + v[172] * v[272] + v[144] * v[381] + v[143] * v[384] + v[142] * v[387];
	v[407] = v[170] * v[265] + v[171] * v[268] + v[172] * v[271] + v[144] * v[380] + v[143] * v[383] + v[142] * v[386];
	v[151] = ud[0] + xpolei[0];
	v[152] = ud[1] + xpolei[1];
	v[153] = ud[2] + xpolei[2];
	v[157] = v[142] * v[158] + v[143] * v[159] + v[144] * v[160];
	v[471] = v[157] * v[356];
	v[161] = v[145] * v[158] + v[146] * v[159] + v[147] * v[160];
	v[473] = v[161] * v[356];
	v[162] = v[148] * v[158] + v[149] * v[159] + v[150] * v[160];
	v[475] = v[162] * v[356];
	v[163] = v[142] * v[164] + v[143] * v[165] + v[144] * v[166];
	v[470] = -(v[163] * v[355]);
	v[167] = v[145] * v[164] + v[146] * v[165] + v[147] * v[166];
	v[472] = -(v[167] * v[355]);
	v[168] = v[148] * v[164] + v[149] * v[165] + v[150] * v[166];
	v[474] = -(v[168] * v[355]);
	v[169] = v[142] * v[170] + v[143] * v[171] + v[144] * v[172];
	v[461] = v[110] * v[169] + v[157] * v[227] + v[163] * v[230];
	v[455] = v[169] * v[360];
	v[449] = v[110] * v[163] - v[169] * v[230] + v[157] * v[233];
	v[443] = v[110] * v[157] - v[169] * v[227] - v[163] * v[233];
	v[173] = v[145] * v[170] + v[146] * v[171] + v[147] * v[172];
	v[462] = v[110] * v[173] + v[161] * v[227] + v[167] * v[230];
	v[456] = v[173] * v[360];
	v[450] = v[110] * v[167] - v[173] * v[230] + v[161] * v[233];
	v[444] = v[110] * v[161] - v[173] * v[227] - v[167] * v[233];
	v[174] = v[148] * v[170] + v[149] * v[171] + v[150] * v[172];
	v[463] = v[110] * v[174] + v[162] * v[227] + v[168] * v[230];
	v[457] = v[174] * v[360];
	v[451] = v[110] * v[168] - v[174] * v[230] + v[162] * v[233];
	v[445] = v[110] * v[162] - v[174] * v[227] - v[168] * v[233];
	v[175] = br[0] * v[142] + br[1] * v[143] + br[2] * v[144];
	v[509] = -v[151] - v[175];
	v[499] = (v[175] * v[175]);
	v[338] = (*a1)*v[175];
	v[353] = -(v[230] * v[338]);
	v[350] = -(v[110] * v[338]);
	v[176] = br[0] * v[145] + br[1] * v[146] + br[2] * v[147];
	v[524] = (*m)*v[176];
	v[508] = v[152] + v[176];
	v[500] = (v[176] * v[176]);
	v[494] = v[175] * v[524];
	v[337] = -((*a1)*v[176]);
	v[352] = -(v[227] * v[337]);
	v[348] = -(v[110] * v[337]);
	v[177] = br[0] * v[148] + br[1] * v[149] + br[2] * v[150];
	v[505] = -v[153] - v[177];
	v[498] = v[177] * v[524];
	v[497] = (*m)*v[175] * v[177];
	v[495] = (v[177] * v[177]);
	v[339] = -((*a1)*v[177]);
	v[345] = v[233] * v[339];
	v[343] = v[110] * v[339];
	v[178] = (*a6)*ddui[0] + (*a5)*dui[0] + (*a4)*ud[0];
	v[179] = (*a6)*ddui[1] + (*a5)*dui[1] + (*a4)*ud[1];
	v[180] = (*a6)*ddui[2] + (*a5)*dui[2] + (*a4)*ud[2];
	v[181] = -((*a3)*ddui[0]) - (*a2)*dui[0] + (*a1)*ud[0];
	v[182] = -((*a3)*ddui[1]) - (*a2)*dui[1] + (*a1)*ud[1];
	v[417] = v[182] * v[293] - v[181] * v[296];
	v[416] = v[182] * v[292] - v[181] * v[295];
	v[203] = -(v[176] * v[181]) + v[175] * v[182];
	v[458] = v[203] * v[360];
	v[183] = -((*a3)*ddui[2]) - (*a2)*dui[2] + (*a1)*ud[2];
	v[423] = -(v[183] * v[293]) + v[181] * v[299];
	v[422] = -(v[183] * v[292]) + v[181] * v[298];
	v[420] = v[183] * v[296] - v[182] * v[299];
	v[419] = v[183] * v[295] - v[182] * v[298];
	v[205] = -(v[177] * v[182]) + v[176] * v[183];
	v[477] = v[205] * v[356];
	v[204] = v[177] * v[181] - v[175] * v[183];
	v[476] = -(v[204] * v[355]);
	v[466] = (*m)*(v[203] * v[221] + v[204] * v[360] + v[205] * v[361] + v[110] * v[417] + v[227] * v[420] + v[230] * v[423]);
	v[464] = (*m)*(v[203] * v[219] + v[204] * v[359] - v[205] * v[360] + v[110] * v[416] + v[227] * v[419] + v[230] * v[422]);
	v[452] = (*m)*(v[204] * v[219] + v[205] * v[355] - v[203] * v[359] - v[230] * v[416] + v[233] * v[419] + v[110] * v[422]);
	v[184] = v[113] * v[185] + v[114] * v[186] + v[115] * v[187];
	v[484] = -(v[169] * v[184]);
	v[482] = v[163] * v[184];
	v[188] = v[117] * v[185] + v[119] * v[186] + v[120] * v[187];
	v[485] = v[173] * v[188];
	v[480] = -(v[161] * v[188]);
	v[439] = v[168] * v[184] - v[162] * v[188];
	v[438] = v[167] * v[184] + v[480];
	v[437] = -(v[157] * v[188]) + v[482];
	v[309] = -(v[188] * v[294]) + v[184] * v[297] + v[176] * v[303] - v[175] * v[306];
	v[308] = -(v[188] * v[293]) + v[184] * v[296] + v[176] * v[302] - v[175] * v[305];
	v[307] = -(v[188] * v[292]) + v[184] * v[295] + v[176] * v[301] - v[175] * v[304];
	v[197] = v[176] * v[184] - v[175] * v[188];
	v[189] = v[122] * v[185] + v[124] * v[186] + v[125] * v[187];
	v[507] = v[174] * v[189] - v[484] + v[485];
	v[483] = -(v[168] * v[189]);
	v[504] = v[167] * v[188] + v[482] - v[483];
	v[481] = v[162] * v[189];
	v[502] = v[157] * v[184] - v[480] + v[481];
	v[442] = v[184] * (v[163] * v[303] - v[157] * v[306] - v[188] * v[391] + v[184] * v[400]) + v[188] * (v[167] * v[303]
		- v[161] * v[306] - v[188] * v[394] + v[184] * v[403]) + v[189] * (v[168] * v[303] - v[162] * v[306] - v[188] * v[397]
			+ v[184] * v[406]) + v[303] * v[437] + v[306] * v[438] + v[312] * v[439];
	v[441] = v[184] * (v[163] * v[302] - v[157] * v[305] - v[188] * v[390] + v[184] * v[399]) + v[188] * (v[167] * v[302]
		- v[161] * v[305] - v[188] * v[393] + v[184] * v[402]) + v[189] * (v[168] * v[302] - v[162] * v[305] - v[188] * v[396]
			+ v[184] * v[405]) + v[302] * v[437] + v[305] * v[438] + v[311] * v[439];
	v[440] = v[184] * (v[163] * v[301] - v[157] * v[304] - v[188] * v[389] + v[184] * v[398]) + v[188] * (v[167] * v[301]
		- v[161] * v[304] - v[188] * v[392] + v[184] * v[401]) + v[189] * (v[168] * v[301] - v[162] * v[304] - v[188] * v[395]
			+ v[184] * v[404]) + v[301] * v[437] + v[304] * v[438] + v[310] * v[439];
	v[433] = -(v[174] * v[184]) + v[481];
	v[432] = -(v[173] * v[184]) + v[161] * v[189];
	v[431] = v[157] * v[189] + v[484];
	v[436] = v[184] * (-(v[169] * v[303]) + v[157] * v[312] + v[189] * v[391] - v[184] * v[409]) + v[188] * (-(v[173] * v[303])
		+ v[161] * v[312] + v[189] * v[394] - v[184] * v[412]) + v[189] * (-(v[174] * v[303]) + v[162] * v[312] + v[189] * v[397]
			- v[184] * v[415]) + v[303] * v[431] + v[306] * v[432] + v[312] * v[433];
	v[435] = v[184] * (-(v[169] * v[302]) + v[157] * v[311] + v[189] * v[390] - v[184] * v[408]) + v[188] * (-(v[173] * v[302])
		+ v[161] * v[311] + v[189] * v[393] - v[184] * v[411]) + v[189] * (-(v[174] * v[302]) + v[162] * v[311] + v[189] * v[396]
			- v[184] * v[414]) + v[302] * v[431] + v[305] * v[432] + v[311] * v[433];
	v[434] = v[184] * (-(v[169] * v[301]) + v[157] * v[310] + v[189] * v[389] - v[184] * v[407]) + v[188] * (-(v[173] * v[301])
		+ v[161] * v[310] + v[189] * v[392] - v[184] * v[410]) + v[189] * (-(v[174] * v[301]) + v[162] * v[310] + v[189] * v[395]
			- v[184] * v[413]) + v[301] * v[431] + v[304] * v[432] + v[310] * v[433];
	v[427] = v[174] * v[188] + v[483];
	v[426] = -(v[167] * v[189]) + v[485];
	v[425] = v[169] * v[188] - v[163] * v[189];
	v[430] = v[184] * (v[169] * v[306] - v[163] * v[312] - v[189] * v[400] + v[188] * v[409]) + v[188] * (v[173] * v[306]
		- v[167] * v[312] - v[189] * v[403] + v[188] * v[412]) + v[189] * (v[174] * v[306] - v[168] * v[312] - v[189] * v[406]
			+ v[188] * v[415]) + v[303] * v[425] + v[306] * v[426] + v[312] * v[427];
	v[429] = v[184] * (v[169] * v[305] - v[163] * v[311] - v[189] * v[399] + v[188] * v[408]) + v[188] * (v[173] * v[305]
		- v[167] * v[311] - v[189] * v[402] + v[188] * v[411]) + v[189] * (v[174] * v[305] - v[168] * v[311] - v[189] * v[405]
			+ v[188] * v[414]) + v[302] * v[425] + v[305] * v[426] + v[311] * v[427];
	v[428] = v[184] * (v[169] * v[304] - v[163] * v[310] - v[189] * v[398] + v[188] * v[407]) + v[188] * (v[173] * v[304]
		- v[167] * v[310] - v[189] * v[401] + v[188] * v[410]) + v[189] * (v[174] * v[304] - v[168] * v[310] - v[189] * v[404]
			+ v[188] * v[413]) + v[301] * v[425] + v[304] * v[426] + v[310] * v[427];
	v[318] = v[189] * v[294] - v[184] * v[300] - v[177] * v[303] + v[175] * v[312];
	v[317] = v[189] * v[293] - v[184] * v[299] - v[177] * v[302] + v[175] * v[311];
	v[316] = v[189] * v[292] - v[184] * v[298] - v[177] * v[301] + v[175] * v[310];
	v[315] = -(v[189] * v[297]) + v[188] * v[300] + v[177] * v[306] - v[176] * v[312];
	v[314] = -(v[189] * v[296]) + v[188] * v[299] + v[177] * v[305] - v[176] * v[311];
	v[313] = -(v[189] * v[295]) + v[188] * v[298] + v[177] * v[304] - v[176] * v[310];
	v[208] = v[184] * v[425] + v[188] * v[426] + v[189] * v[427];
	v[469] = v[208] * v[356];
	v[207] = v[184] * v[431] + v[188] * v[432] + v[189] * v[433];
	v[468] = -(v[207] * v[355]);
	v[206] = v[184] * v[437] + v[188] * v[438] + v[189] * v[439];
	v[454] = v[206] * v[360];
	v[200] = v[177] * v[188] - v[176] * v[189];
	v[199] = -(v[177] * v[184]) + v[175] * v[189];
	v[190] = v[113] * v[191] + v[114] * v[192] + v[115] * v[193];
	v[194] = v[117] * v[191] + v[119] * v[192] + v[120] * v[193];
	v[195] = v[122] * v[191] + v[124] * v[192] + v[125] * v[193];
	v[479] = 0.5e0*(*m)*((v[178] * v[178]) + (v[179] * v[179]) + (v[180] * v[180]));
	v[486] = 0.5e0*(v[184] * v[502] + v[188] * v[504] + v[189] * v[507]);
	v[487] = (*m)*(v[177] * (-(v[179] * v[184]) + v[178] * v[188]) + v[176] * (v[180] * v[184] - v[178] * v[189]) + v[175] * (-
		(v[180] * v[188]) + v[179] * v[189]));
	v[489] = (*m)*(v[178] + v[200]);
	v[490] = (*m)*(v[179] + v[199]);
	v[491] = (*m)*(v[180] + v[197]);
	v[493] = v[188] * (v[161] + v[494]) + v[189] * (v[162] + v[497]) + v[184] * (v[157] - (*m)*(v[495] + v[500]));
	v[496] = v[184] * (v[163] + v[494]) + v[189] * (v[168] + v[498]) + v[188] * (v[167] - (*m)*(v[495] + v[499]));
	v[501] = v[184] * (v[169] + v[497]) + v[188] * (v[173] + v[498]) + v[189] * (v[174] - (*m)*(v[499] + v[500]));
	dT[0] = (*m)*(v[181] + v[177] * v[194] - v[176] * v[195] + v[188] * v[197] - v[189] * v[199]);
	dT[1] = (*m)*(v[182] - v[177] * v[190] + v[175] * v[195] - v[184] * v[197] + v[189] * v[200]);
	dT[2] = (*m)*(v[183] + v[176] * v[190] - v[175] * v[194] + v[184] * v[199] - v[188] * v[200]);
	dT[3] = v[110] * v[208] - v[206] * v[227] - v[207] * v[233] + (*m)*(v[110] * v[205] - v[203] * v[227] - v[204] * v[233])
		+ v[190] * v[443] + v[194] * v[444] + v[195] * v[445];
	dT[4] = v[110] * v[207] - v[206] * v[230] + v[208] * v[233] + (*m)*(v[110] * v[204] - v[203] * v[230] + v[205] * v[233])
		+ v[190] * v[449] + v[194] * v[450] + v[195] * v[451];
	dT[5] = v[110] * v[206] + v[208] * v[227] + v[207] * v[230] + (*m)*(v[110] * v[203] + v[205] * v[227] + v[204] * v[230])
		+ v[190] * v[461] + v[194] * v[462] + v[195] * v[463];
	DdT[0][0] = v[211];
	DdT[0][1] = 0e0;
	DdT[0][2] = 0e0;
	DdT[0][3] = (*m)*(-(v[195] * v[295]) + v[194] * v[298] + v[197] * v[304] + v[188] * v[307] - v[199] * v[310]
		- v[189] * v[316] + v[177] * v[322] - v[176] * v[325]);
	DdT[0][4] = (*m)*(-(v[195] * v[296]) + v[194] * v[299] + v[197] * v[305] + v[188] * v[308] - v[199] * v[311]
		- v[189] * v[317] + v[177] * v[323] - v[176] * v[326]);
	DdT[0][5] = (*m)*(-(v[195] * v[297]) + v[194] * v[300] + v[197] * v[306] + v[188] * v[309] - v[199] * v[312]
		- v[189] * v[318] + v[177] * v[324] - v[176] * v[327]);
	DdT[1][0] = 0e0;
	DdT[1][1] = v[211];
	DdT[1][2] = 0e0;
	DdT[1][3] = (*m)*(v[195] * v[292] - v[190] * v[298] - v[197] * v[301] - v[184] * v[307] + v[200] * v[310] + v[189] * v[313]
		- v[177] * v[319] + v[175] * v[325]);
	DdT[1][4] = (*m)*(v[195] * v[293] - v[190] * v[299] - v[197] * v[302] - v[184] * v[308] + v[200] * v[311] + v[189] * v[314]
		- v[177] * v[320] + v[175] * v[326]);
	DdT[1][5] = (*m)*(v[195] * v[294] - v[190] * v[300] - v[197] * v[303] - v[184] * v[309] + v[200] * v[312] + v[189] * v[315]
		- v[177] * v[321] + v[175] * v[327]);
	DdT[2][0] = 0e0;
	DdT[2][1] = 0e0;
	DdT[2][2] = v[211];
	DdT[2][3] = (*m)*(-(v[194] * v[292]) + v[190] * v[295] + v[199] * v[301] - v[200] * v[304] - v[188] * v[313]
		+ v[184] * v[316] + v[176] * v[319] - v[175] * v[322]);
	DdT[2][4] = (*m)*(-(v[194] * v[293]) + v[190] * v[296] + v[199] * v[302] - v[200] * v[305] - v[188] * v[314]
		+ v[184] * v[317] + v[176] * v[320] - v[175] * v[323]);
	DdT[2][5] = (*m)*(-(v[194] * v[294]) + v[190] * v[297] + v[199] * v[303] - v[200] * v[306] - v[188] * v[315]
		+ v[184] * v[318] + v[176] * v[321] - v[175] * v[324]);
	DdT[3][0] = (*m)*(v[345] + v[352]);
	DdT[3][1] = (*m)*(-(v[227] * v[338]) + v[343]);
	DdT[3][2] = (*m)*(v[233] * v[338] + v[348]);
	DdT[3][3] = v[208] * v[219] + v[110] * v[428] - v[233] * v[434] - v[227] * v[440] + v[319] * v[443] + v[322] * v[444]
		+ v[325] * v[445] + v[454] + v[468] + v[190] * (v[157] * v[219] + v[110] * v[389] - v[233] * v[398] - v[227] * v[407] + v[455]
			+ v[470]) + v[194] * (v[161] * v[219] + v[110] * v[392] - v[233] * v[401] - v[227] * v[410] + v[456] + v[472]) + v[195] *
			(v[162] * v[219] + v[110] * v[395] - v[233] * v[404] - v[227] * v[413] + v[457] + v[474]) + (*m)*(v[205] * v[219]
				- v[227] * v[416] + v[110] * v[419] - v[233] * v[422] + v[458] + v[476]);
	DdT[3][4] = v[208] * v[221] - v[207] * v[356] - v[206] * v[361] + v[190] * (v[157] * v[221] - v[163] * v[356]
		- v[169] * v[361] + v[110] * v[390] - v[233] * v[399] - v[227] * v[408]) + v[194] * (v[161] * v[221] - v[167] * v[356]
			- v[173] * v[361] + v[110] * v[393] - v[233] * v[402] - v[227] * v[411]) + v[195] * (v[162] * v[221] - v[168] * v[356]
				- v[174] * v[361] + v[110] * v[396] - v[233] * v[405] - v[227] * v[414]) + v[110] * v[429] - v[233] * v[435] - v[227] * v[441]
		+ v[320] * v[443] + v[323] * v[444] + v[326] * v[445] + v[452];
	DdT[3][5] = v[208] * v[222] + v[206] * v[356] - v[207] * v[357] + v[190] * (v[157] * v[222] + v[169] * v[356]
		- v[163] * v[357] + v[110] * v[391] - v[233] * v[400] - v[227] * v[409]) + v[194] * (v[161] * v[222] + v[173] * v[356]
			- v[167] * v[357] + v[110] * v[394] - v[233] * v[403] - v[227] * v[412]) + v[195] * (v[162] * v[222] + v[174] * v[356]
				- v[168] * v[357] + v[110] * v[397] - v[233] * v[406] - v[227] * v[415]) + v[110] * v[430] - v[233] * v[436] - v[227] * v[442]
		+ v[321] * v[443] + v[324] * v[444] + v[327] * v[445] + v[464];
	DdT[4][0] = (*m)*(-(v[230] * v[337]) - v[343]);
	DdT[4][1] = (*m)*(v[345] + v[353]);
	DdT[4][2] = (*m)*(-(v[233] * v[337]) + v[350]);
	DdT[4][3] = v[207] * v[219] + v[208] * v[355] - v[206] * v[359] + v[190] * (v[163] * v[219] + v[157] * v[355]
		- v[169] * v[359] + v[233] * v[389] + v[110] * v[398] - v[230] * v[407]) + v[194] * (v[167] * v[219] + v[161] * v[355]
			- v[173] * v[359] + v[233] * v[392] + v[110] * v[401] - v[230] * v[410]) + v[195] * (v[168] * v[219] + v[162] * v[355]
				- v[174] * v[359] + v[233] * v[395] + v[110] * v[404] - v[230] * v[413]) + v[233] * v[428] + v[110] * v[434] - v[230] * v[440]
		+ v[319] * v[449] + v[322] * v[450] + v[325] * v[451] + v[452];
	DdT[4][4] = v[207] * v[221] + v[233] * v[429] + v[110] * v[435] - v[230] * v[441] + v[320] * v[449] + v[323] * v[450]
		+ v[326] * v[451] - v[454] + v[469] + v[190] * (v[163] * v[221] + v[233] * v[390] + v[110] * v[399] - v[230] * v[408] - v[455]
			+ v[471]) + v[194] * (v[167] * v[221] + v[233] * v[393] + v[110] * v[402] - v[230] * v[411] - v[456] + v[473]) + v[195] *
			(v[168] * v[221] + v[233] * v[396] + v[110] * v[405] - v[230] * v[414] - v[457] + v[475]) + (*m)*(v[204] * v[221]
				- v[230] * v[417] + v[233] * v[420] + v[110] * v[423] - v[458] + v[477]);
	DdT[4][5] = v[207] * v[222] - v[206] * v[355] + v[208] * v[357] + v[190] * (v[163] * v[222] - v[169] * v[355]
		+ v[157] * v[357] + v[233] * v[391] + v[110] * v[400] - v[230] * v[409]) + v[194] * (v[167] * v[222] - v[173] * v[355]
			+ v[161] * v[357] + v[233] * v[394] + v[110] * v[403] - v[230] * v[412]) + v[195] * (v[168] * v[222] - v[174] * v[355]
				+ v[162] * v[357] + v[233] * v[397] + v[110] * v[406] - v[230] * v[415]) + v[233] * v[430] + v[110] * v[436] - v[230] * v[442]
		+ v[321] * v[449] + v[324] * v[450] + v[327] * v[451] + v[466];
	DdT[5][0] = (*m)*(-(v[230] * v[339]) - v[348]);
	DdT[5][1] = (*m)*(v[227] * v[339] - v[350]);
	DdT[5][2] = (*m)*(v[352] + v[353]);
	DdT[5][3] = v[206] * v[219] + v[207] * v[359] - v[208] * v[360] + v[190] * (v[169] * v[219] + v[163] * v[359]
		- v[157] * v[360] + v[227] * v[389] + v[230] * v[398] + v[110] * v[407]) + v[194] * (v[173] * v[219] + v[167] * v[359]
			- v[161] * v[360] + v[227] * v[392] + v[230] * v[401] + v[110] * v[410]) + v[195] * (v[174] * v[219] + v[168] * v[359]
				- v[162] * v[360] + v[227] * v[395] + v[230] * v[404] + v[110] * v[413]) + v[227] * v[428] + v[230] * v[434] + v[110] * v[440]
		+ v[319] * v[461] + v[322] * v[462] + v[325] * v[463] + v[464];
	DdT[5][4] = v[206] * v[221] + v[207] * v[360] + v[208] * v[361] + v[190] * (v[169] * v[221] + v[163] * v[360]
		+ v[157] * v[361] + v[227] * v[390] + v[230] * v[399] + v[110] * v[408]) + v[194] * (v[173] * v[221] + v[167] * v[360]
			+ v[161] * v[361] + v[227] * v[393] + v[230] * v[402] + v[110] * v[411]) + v[195] * (v[174] * v[221] + v[168] * v[360]
				+ v[162] * v[361] + v[227] * v[396] + v[230] * v[405] + v[110] * v[414]) + v[227] * v[429] + v[230] * v[435] + v[110] * v[441]
		+ v[320] * v[461] + v[323] * v[462] + v[326] * v[463] + v[466];
	DdT[5][5] = v[206] * v[222] + v[227] * v[430] + v[230] * v[436] + v[110] * v[442] + v[321] * v[461] + v[324] * v[462]
		+ v[327] * v[463] - v[468] - v[469] + v[190] * (v[169] * v[222] + v[227] * v[391] + v[230] * v[400] + v[110] * v[409] - v[470]
			- v[471]) + v[194] * (v[173] * v[222] + v[227] * v[394] + v[230] * v[403] + v[110] * v[412] - v[472] - v[473]) + v[195] *
			(v[174] * v[222] + v[227] * v[397] + v[230] * v[406] + v[110] * v[415] - v[474] - v[475]) + (*m)*(v[203] * v[222] + v[110] *
		(v[182] * v[294] - v[181] * v[297]) + v[230] * (-(v[183] * v[294]) + v[181] * v[300]) + v[227] * (v[183] * v[297]
			- v[182] * v[300]) - v[476] - v[477]);
	//(*T1) = v[479];
	//(*T2) = v[486];
	//(*T3) = v[487];
	//(*T) = v[479] + v[486] + v[487];
	//L[0] = v[489];
	//L[1] = v[490];
	//L[2] = v[491];
	//(*magL) = sqrt((v[489] * v[489]) + (v[490] * v[490]) + (v[491] * v[491]));
	//HG[0] = v[493];
	//HG[1] = v[496];
	//HG[2] = v[501];
	//HO[0] = v[502] + (*m)*(v[152] * v[197] - v[153] * v[199] + v[179] * v[505] + v[180] * v[508]);
	//HO[1] = v[504] + (*m)*(-(v[151] * v[197]) + v[153] * v[200] - v[178] * v[505] + v[180] * v[509]);
	//HO[2] = v[507] + (*m)*(v[151] * v[199] - v[152] * v[200] - v[178] * v[508] - v[179] * v[509]);
	//(*magHG) = sqrt((v[493] * v[493]) + (v[496] * v[496]) + (v[501] * v[501]));

	//Kinetic Energy
	kinetic_energy = v[479] + v[486] + v[487];
	//Linear momentum
	linear_momentum[0] = v[489];
	linear_momentum[1] = v[490];
	linear_momentum[2] = v[491];
	linear_momentum_mag = sqrt((v[489] * v[489]) + (v[490] * v[490]) + (v[491] * v[491]));
	//Angular momentum
	angular_momentum[0] = v[493];
	angular_momentum[1] = v[496];
	angular_momentum[2] = v[501];
	angular_momentum_mag = sqrt((v[493] * v[493]) + (v[496] * v[496]) + (v[501] * v[501]));
	angular_momentum_origin[0] = v[502] + (*m)*(v[152] * v[197] - v[153] * v[199] + v[179] * v[505] + v[180] * v[508]);
	angular_momentum_origin[1] = v[504] + (*m)*(-(v[151] * v[197]) + v[153] * v[200] - v[178] * v[505] + v[180] * v[509]);
	angular_momentum_origin[2] = v[507] + (*m)*(v[151] * v[199] - v[152] * v[200] - v[178] * v[508] - v[179] * v[509]);
}
