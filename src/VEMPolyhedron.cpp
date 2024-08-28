#include "VEMPolyhedron.h"

#include "Matrix.h"
#include "MatrixFloat.h"
#include "BoundingSphere.h"
#include "CADData.h"
#include "STLSurface.h"
#include "CoordinateSystem.h"
#include "SuperNode.h"
#include "BoundingTriangularBox.h"
#include "BoundingCylinder.h"
#include "Hooke.h"
#include "TriangularFace.h"
#include "GeneralContactSearch.h"
#include "Environment.h"
#include "Interface_1.h"
#include "Encoding.h"
#include "Dynamic.h"
#include"Database.h"
//Variáveis globais
extern
Database db;

VEMPolyhedron::VEMPolyhedron()
{
	material = 0;
	super_node = 0;
	number = 0;
	node = 0;
	cs = 0;
	kinetic_energy = 0.0;
	strain_energy = 0.0;
	potential_g_energy = 0.0;
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
	n_vertices = 0;
	n_sub_bv = 0;
	n_sub_bv = 0;
	CADDATA_ID = 0;
	//bv_factor = 0.1f;
	inc_len_factor = 0.1f;
	deformed_volume = 0.0;

	//Load factor
	l_factor = 0.0;

	sub_bv = NULL;
	vertices_i = NULL;
	vertices_p = NULL;

	/*count_local_abs_load = NULL;
	local_abs_load = NULL;*/
	local_residual = NULL;
	local_e_load = NULL;
	local_stiffness = NULL;
	local_mass = NULL;
	local_damping = NULL;
	
	deformed_barycenter = new MatrixFloat(3);
	deformed_faces_radii = NULL;
	deformed_edges_lengths = NULL;
	deformed_volume = 0.0f;

	domainData = NULL;														//Domain data
	elementData = NULL;														//Element data
	
	type_name = new char[20];//Nome do tipo da partícula
	sprintf(type_name, "VEMPolyhedron");

	I3 = new Matrix(3, 3);
	(*I3)(0, 0) = 1.0;
	(*I3)(1, 1) = 1.0;
	(*I3)(2, 2) = 1.0;
	Q0 = new Matrix(3, 3);
	bv = new BoundingSphere();
	x0i = new Matrix(3);
	x0p = new Matrix(3);

	//control of specific allocation
	allocced = false;

	//AceGen
	Xref = NULL;
	Xcur = NULL;
	displacement = NULL;
	triangles = NULL;
	tetrahedrons = NULL;
	ht = NULL;
	hp = NULL;
	residual = NULL;
	residualload = NULL;
	elementPostprocessing = NULL;

	tangent = NULL;
	massmatrix = NULL;
	//Stabilization type: 1-FEM 2-standard
	stab_flag = 1;
	//Stabilization values (default)
	stab_factor = 0.5;
	stab_factor_mass = 1.0;
	stab_factor_load = 1.0;
	//Lumped mass flag
	lumped_mass_flag = false;

	//Damping
	damping_flag = false;
	first_damping_evaluation = true;
	damping_factor = 0.0;
	alpha = 0.0;
	beta = 0.0;
}

VEMPolyhedron::~VEMPolyhedron()
{
	Free();
	delete[] type_name;
	delete I3;
	delete Q0;
	delete bv;
	delete x0i;
	delete x0p;
	delete deformed_barycenter;
	if (sub_bv != NULL)
	{
		for (int i = 0; i < n_sub_bv; i++)
		{
			delete sub_bv[i];
		}
		delete[] sub_bv;
	}
}

bool VEMPolyhedron::Check()
{
	if (super_node > db.number_super_nodes)
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

bool VEMPolyhedron::Read(FILE *f)
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
	if (!strcmp(s, "SuperNode"))
	{
		fscanf(f, "%s", s);
		super_node = atoi(s);
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "StabFlag"))
	{
		fscanf(f, "%s", s);
		stab_flag = atoi(s);
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "StabBetaStiff"))
	{
		fscanf(f, "%s", s);
		stab_factor = atof(s);
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "StabBetaMassLoad"))
	{
		fscanf(f, "%s", s);
		stab_factor_mass = atof(s);
		stab_factor_load = atof(s);
	}
	else
	{
		if (!strcmp(s, "StabBetaMass"))
		{
			fscanf(f, "%s", s);
			stab_factor_mass = atof(s);
		}
		else
			return false;
		fscanf(f, "%s", s);
		if (!strcmp(s, "StabBetaLoad"))
		{
			fscanf(f, "%s", s);
			stab_factor_load = atof(s);
		}
		else
			return false;
	}

	fscanf(f, "%s", s);
	if (!strcmp(s, "LumpedMass"))
	{
		fscanf(f, "%s", s);
		lumped_mass_flag = (bool)atoi(s);
	}
	else
		return false;

	//Salva a posição (stream)
	fpos_t pos;
	fgetpos(f, &pos);
	fscanf(f, "%s", s);
	if (!strcmp(s, "RayleighDamping"))
	{
		damping_flag = true;
		fscanf(f, "%s", s);
		if (!strcmp(s, "Alpha"))
		{
			fscanf(f, "%s", s);
			alpha = atof(s);
		}
		else
			return false;

		fscanf(f, "%s", s);
		if (!strcmp(s, "Beta"))
		{
			fscanf(f, "%s", s);
			beta = atof(s);
		}
		else
			return false;
	}
	else
	{
		fsetpos(f, &pos);
		damping_flag = false;
		alpha = 0.0;
		beta = 0.0;
	}

	return true;
}

void VEMPolyhedron::Write(FILE *f)
{
	fprintf(f, "VEMPolyhedron\t%d\tMat\t%d\tCS\t%d\tCADData\t%d\tSuperNode\t%d\tStabFlag\t%d\tStabBetaStiff\t%.6e\tStabBetaMassLoad\t%.6e\tLumpedMass\t%d\tRayleighDamping\tAlpha\t%.6e\tBeta\t%.6e\n",
		number,
		material,
		cs,
		CADDATA_ID,
		super_node,
		stab_flag,
		stab_factor,
		stab_factor_mass,
		(int)lumped_mass_flag,
		alpha,
		beta
	);
}
void VEMPolyhedron::WriteModifyingParameters(FILE *f, int e_material, int e_node, int e_number, int e_cs)
{
	//TODO (used in saving configuration)
}

void VEMPolyhedron::PreCalc()
{
	////////////////////////////////////////////////////////////////
	//Matriz para transformação de coordenadas  local-global
	// vlocal = (sistema do CAD)
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

	//CAD pointer
	STLSurface* ptr_cad = static_cast<STLSurface*>(db.cad_data[CADDATA_ID - 1]);
	n_vertices = (int)ptr_cad->vertices.size();
	n_faces = (int)ptr_cad->n_faces;
	n_edges = (int)ptr_cad->edges.size();
	//Allocs variables of VEMPolyhedron class
	Alloc();
	//Setting number of DOFs on the corresponding super node
	db.super_nodes[super_node - 1]->Alloc(n_vertices * 3, 0);
	//Filling coordinates on the super node associated
	Matrix xv(3);
	(*x0i)(0, 0) = db.super_nodes[super_node - 1]->super_node_coordinates[0];
	(*x0i)(1, 0) = db.super_nodes[super_node - 1]->super_node_coordinates[1];
	(*x0i)(2, 0) = db.super_nodes[super_node - 1]->super_node_coordinates[2];
	for (int i = 0; i < n_vertices; i++)
	{
		xv(0, 0) = (*ptr_cad->vertices[i].coord_double)(0, 0);
		xv(1, 0) = (*ptr_cad->vertices[i].coord_double)(1, 0);
		xv(2, 0) = (*ptr_cad->vertices[i].coord_double)(2, 0);
		xv = (*x0i) + (*Q0) * xv;	//to global CS
		db.super_nodes[super_node - 1]->ref_coordinates[3 * i] = xv(0, 0);
		db.super_nodes[super_node - 1]->ref_coordinates[3 * i + 1] = xv(1, 0);
		db.super_nodes[super_node - 1]->ref_coordinates[3 * i + 2] = xv(2, 0);

		db.super_nodes[super_node - 1]->copy_coordinates[3 * i] = xv(0, 0);
		db.super_nodes[super_node - 1]->copy_coordinates[3 * i + 1] = xv(1, 0);
		db.super_nodes[super_node - 1]->copy_coordinates[3 * i + 2] = xv(2, 0);

		//Local variables
		*vertices_i[i] = xv;
		*vertices_p[i] = xv;
	}
	*x0p = *x0i;
	UpdateVariables();

	//Allocation of subBV's
	//Sub bounding volumes (faces+edges+vertices)
	n_sub_bv = ptr_cad->n_faces + (int)ptr_cad->edges.size() + (int)ptr_cad->vertices.size();
	sub_bv = new BoundingVolume*[n_sub_bv];
	//BV on faces
	for (int i = 0; i < ptr_cad->n_faces; i++)
		sub_bv[i] = new BoundingTriangularBox();
	//BV on edges
	for (int i = ptr_cad->n_faces; i < ptr_cad->n_faces + (int)ptr_cad->edges.size(); i++)
		sub_bv[i] = new BoundingCylinder();
	//BV on vertices
	for (int i = ptr_cad->n_faces + (int)ptr_cad->edges.size(); i < ptr_cad->n_faces + (int)ptr_cad->edges.size() + (int)ptr_cad->vertices.size(); i++)
		sub_bv[i] = new BoundingSphere();
		
	UpdateGeometricEntities();//BV's initial setup
	UpdateBoundingVolumes();
	SaveLagrange();
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	//Mass
	mass = ptr_cad->volume*db.materials[material - 1]->rho;
	//////////////AceGen variables///////////////////
	Xref = db.super_nodes[super_node - 1]->ref_coordinates;			//Reference coordinates
	Xcur = db.super_nodes[super_node - 1]->copy_coordinates;		//Current coordinates
	displacement = db.super_nodes[super_node - 1]->displacements;	//Displacements
	
	elementData[0] = n_vertices;
	elementData[1] = ptr_cad->n_faces;
	elementData[2] = (double)ptr_cad->tetras.size();
	for (int i = 0; i < 26; i++)
		elementPostprocessing[i] = 0.0;

	/*
	{"E -elastic modulus", "\[Nu] -Poisson ratio", "rho - specific mass", \
	"qX body force X", "qY body force Y", "qZ body force Z", "\[Beta]SE \
	-beta (VE-FE ratio) for Solid static", "\[Beta]mass -beta (VE-FE \
	ratio) for Solid mass", "\[Beta]load -beta (VE-FE ratio) for Solid \
	load", "st -stabilisation type: 1-FE, 2-standard"}
	*/
	if (typeid(*db.materials[material - 1]) == typeid(Hooke))
	{
		Hooke* ptr = static_cast<Hooke*>(db.materials[material - 1]);
		domainData[0] = ptr->E;
		domainData[1] = ptr->nu;
		domainData[2] = ptr->rho;
		//Stabilization
		domainData[6] = stab_factor;		//Stabilization factor
		domainData[7] = stab_factor_mass;	//Stabilization factor
		domainData[8] = stab_factor_load;	//Stabilization factor
		domainData[9] = stab_flag;			//Stabilization type: 1-FEM 2-standard
	}
	
	residual = local_residual;
	tangent = local_stiffness;
	massmatrix = local_mass;
	residualload = local_e_load;

	triangles = new int[ptr_cad->n_faces * 3];
	for (int f = 0; f < ptr_cad->n_faces; f++)
	{
		triangles[f * 3 + 0] = ptr_cad->faces[f]->verticesIDs[0];
		triangles[f * 3 + 1] = ptr_cad->faces[f]->verticesIDs[1];
		triangles[f * 3 + 2] = ptr_cad->faces[f]->verticesIDs[2];
	}
	tetrahedrons = new int[(int)ptr_cad->tetras.size() * 4];
	for (int t = 0; t < (int)ptr_cad->tetras.size(); t++)
	{
		tetrahedrons[t * 4 + 0] = ptr_cad->tetras[t].verticesIDs[0];
		tetrahedrons[t * 4 + 1] = ptr_cad->tetras[t].verticesIDs[1];
		tetrahedrons[t * 4 + 2] = ptr_cad->tetras[t].verticesIDs[2];
		tetrahedrons[t * 4 + 3] = ptr_cad->tetras[t].verticesIDs[3];
	}
}

void VEMPolyhedron::UpdateBoundingVolumes()
{
	UpdateGeometricEntities();//BV's initial setup
	//////////////////////////////////////Bounding Volumes////////////////////////////////////
	if (db.gcs_exist)
		inc_len_factor = db.gcs->inc_len_factor;
	STLSurface* ptr_cad = static_cast<STLSurface*>(db.cad_data[CADDATA_ID - 1]);

	BoundingSphere* ptr_bv1 = static_cast<BoundingSphere*>(bv);
	BoundingTriangularBox* ptr_bv2;
	BoundingCylinder* ptr_bv3;
	BoundingSphere* ptr_bv4;
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//Updating the main bounding volume
	*ptr_bv1->center = *x0p;
	//ptr_bv1->Report();
	//Updating the sub bounding volumes
	
	//BV on faces
	for (int i = 0; i < ptr_cad->n_faces; i++)
	{
		ptr_bv2 = static_cast<BoundingTriangularBox*>(sub_bv[i]);
		*ptr_bv2->x0 = *vertices_p[ptr_cad->faces[i]->verticesIDs[0] - 1];
		*ptr_bv2->x1 = *vertices_p[ptr_cad->faces[i]->verticesIDs[1] - 1];
		*ptr_bv2->x2 = *vertices_p[ptr_cad->faces[i]->verticesIDs[2] - 1];
		//ptr_bv2->Report();
	}
	//BV on edges
	for (int i = ptr_cad->n_faces; i < ptr_cad->n_faces + (int)ptr_cad->edges.size(); i++)
	{
		ptr_bv3 = static_cast<BoundingCylinder*>(sub_bv[i]);
		*ptr_bv3->xb = *vertices_p[ptr_cad->edges[i - ptr_cad->n_faces].verticesIDs[0] - 1];
		*ptr_bv3->xt = *vertices_p[ptr_cad->edges[i - ptr_cad->n_faces].verticesIDs[1] - 1];
		//ptr_bv3->Report();
	}
	//BV on vertices
	for (int i = ptr_cad->n_faces + (int)ptr_cad->edges.size(); i < ptr_cad->n_faces + (int)ptr_cad->edges.size() + (int)ptr_cad->vertices.size(); i++)
	{
		ptr_bv4 = static_cast<BoundingSphere*>(sub_bv[i]);
		*ptr_bv4->center = *vertices_p[i - ptr_cad->n_faces - (int)ptr_cad->edges.size()];
		//ptr_bv4->Report();
	}
}

void VEMPolyhedron::SaveLagrange()
{
	PostProcessing(); 
	for (int i=0;i<n_vertices;i++)
		*vertices_i[i] = *vertices_p[i];
	*x0i = *x0p;

	//Saving bounding volumes
	bv->SaveConfiguration();
	for (int i = 0; i < n_sub_bv; i++)
		sub_bv[i]->SaveConfiguration();

	//CAD pointer
	STLSurface* ptr_cad = static_cast<STLSurface*>(db.cad_data[CADDATA_ID - 1]);
	//Kinetic energy 
	kinetic_energy = 0.0;
	//Strain energy
	strain_energy = elementPostprocessing[28];
	/*linear_momentum[0]= elementPostprocessing[26];
	linear_momentum[1] = elementPostprocessing[27];
	linear_momentum[2] = elementPostprocessing[28];*/
	if (lumped_mass_flag)
	{
		for (int i = 0; i < n_vertices; i++)
		{
			kinetic_energy += 0.5*mass * ptr_cad->vertice_factors[i] *
				(db.super_nodes[super_node - 1]->copy_vel[i * 3 + 0] * db.super_nodes[super_node - 1]->copy_vel[i * 3 + 0] +
					db.super_nodes[super_node - 1]->copy_vel[i * 3 + 1] * db.super_nodes[super_node - 1]->copy_vel[i * 3 + 1] +
					db.super_nodes[super_node - 1]->copy_vel[i * 3 + 2] * db.super_nodes[super_node - 1]->copy_vel[i * 3 + 2]);
			
		}
	}
	else
	{
		Matrix temp_vel(3 * n_vertices);
		Matrix temp_mass(3 * n_vertices, 3 * n_vertices);
		for (int i = 0; i < 3 * n_vertices; i++)
		{
			temp_vel(i, 0) = db.super_nodes[super_node - 1]->copy_vel[i];
			for (int j = 0; j < 3 * n_vertices; j++)
				temp_mass(i, j) = local_mass[i][j];
		}
		kinetic_energy = 0.5*(transp(temp_vel)*temp_mass*temp_vel)(0,0);
	}

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
			////Center of mass
			center = *x0p;
		}
	}
	//printf("\n%.6f\n", deformed_volume);
	//printf("%.6f  %.6f  %.6f\n\n", (*deformed_barycenter)(0, 0), (*deformed_barycenter)(1, 0), (*deformed_barycenter)(2, 0));
	potential_g_energy = -mass * dot(g, center);

	first_damping_evaluation = true;
}
void VEMPolyhedron::WriteVTK_XMLBase(FILE *f)
{
	//Not used
}

void VEMPolyhedron::UpdateVolumeAndBarycenter()
{
	//CAD pointer
	STLSurface* ptr_cad = static_cast<STLSurface*>(db.cad_data[CADDATA_ID - 1]);

	///////////VOLUME////////////
	deformed_volume = 0.0;
	MatrixFloat a(3), b(3), ab(3), bc(3), ca(3), c(3), n(3);
	for (int i = 0; i < ptr_cad->n_faces; i++)
	{
		a = *vertices_p[ptr_cad->faces[i]->verticesIDs[0] - 1];
		b = *vertices_p[ptr_cad->faces[i]->verticesIDs[1] - 1];
		c = *vertices_p[ptr_cad->faces[i]->verticesIDs[2] - 1];
		n = cross(b - a, c - a);
		deformed_volume += dot(a, n) / 6;
	}
	///////////CENTROID/////////////
	zeros(deformed_barycenter);
	for (int i = 0; i < ptr_cad->n_faces; i++)
	{
		a = *vertices_p[ptr_cad->faces[i]->verticesIDs[0] - 1];
		b = *vertices_p[ptr_cad->faces[i]->verticesIDs[1] - 1];
		c = *vertices_p[ptr_cad->faces[i]->verticesIDs[2] - 1];
		n = cross(b - a, c - a);
		ab = a + b;
		bc = b + c;
		ca = c + a;
		for (int j = 0; j < 3; j++)
		{
			(*deformed_barycenter)(j, 0) = (*deformed_barycenter)(j, 0) + n(j, 0)*(ab(j, 0)*ab(j, 0) + bc(j, 0)*bc(j, 0) + ca(j, 0)*ca(j, 0)) / (24 * 2 * deformed_volume);
		}
	}
}


void VEMPolyhedron::UpdateGeometricEntities()
{
	//CAD pointer
	STLSurface* ptr_cad = static_cast<STLSurface*>(db.cad_data[CADDATA_ID - 1]);
	//////////Sphere radius///////////
	deformed_radius = 0.0;
	float temp_len = 0.0;
	for (int i = 0; i < n_vertices; i ++)
	{
		temp_len = sqrt(dot((*vertices_p[i]) + (-1)*(*deformed_barycenter), (*vertices_p[i]) + (-1)*(*deformed_barycenter)));
		if (temp_len > deformed_radius)
			deformed_radius = temp_len;
	}
	///////////Faces radii////////////
	MatrixFloat temp_centroid(3);
	for (int i = 0; i < n_faces; i++)
	{
		temp_centroid = 0.33333333333333 * (*vertices_p[ptr_cad->faces[i]->verticesIDs[0] - 1]
			+ *vertices_p[ptr_cad->faces[i]->verticesIDs[1] - 1]
			+ *vertices_p[ptr_cad->faces[i]->verticesIDs[2] - 1]);
		deformed_faces_radii[i] = 0.0;
		float temp_len = 0.0;
		temp_len = sqrt(dot(*vertices_p[ptr_cad->faces[i]->verticesIDs[0] - 1] + (-1)* temp_centroid, *vertices_p[ptr_cad->faces[i]->verticesIDs[0] - 1] + (-1)* temp_centroid));
		if (temp_len > deformed_faces_radii[i])
			deformed_faces_radii[i] = temp_len;
		temp_len = sqrt(dot(*vertices_p[ptr_cad->faces[i]->verticesIDs[1] - 1] + (-1)* temp_centroid, *vertices_p[ptr_cad->faces[i]->verticesIDs[1] - 1] + (-1)* temp_centroid));
		if (temp_len > deformed_faces_radii[i])
			deformed_faces_radii[i] = temp_len;
		temp_len = sqrt(dot(*vertices_p[ptr_cad->faces[i]->verticesIDs[2] - 1] + (-1)* temp_centroid, *vertices_p[ptr_cad->faces[i]->verticesIDs[2] - 1] + (-1)* temp_centroid));
		if (temp_len > deformed_faces_radii[i])
			deformed_faces_radii[i] = temp_len;
	}
	///////////Edges lengths////////////
	for (int i = 0; i < n_edges; i++)
	{
		deformed_edges_lengths[i] = (float)sqrt(dot(*vertices_p[ptr_cad->edges[i].verticesIDs[0] - 1] - *vertices_p[ptr_cad->edges[i].verticesIDs[1] - 1],
			*vertices_p[ptr_cad->edges[i].verticesIDs[0] - 1] - *vertices_p[ptr_cad->edges[i].verticesIDs[1] - 1]));
	}

	//Computing the largest gnb available, from interface laws (to avoid smaller BVs of small geometric entities)
	double largest_gnb = 0.0;
	for (int i = 0; i < db.number_contactinterfaces; i++)
	{
		Interface_1* ptr = static_cast<Interface_1*>(db.contactinterfaces[i]);
		if (ptr->gnb > largest_gnb)
			largest_gnb = ptr->gnb;
	}

	BoundingSphere* ptr_bv1 = static_cast<BoundingSphere*>(bv);
	//Initial settings of the spherical bounding volume
	ptr_bv1->radius = (1.0f + inc_len_factor) * deformed_radius;
	ptr_bv1->size = ptr_bv1->radius;
	//Associated entity:
	ptr_bv1->associated_type = 'P';
	ptr_bv1->associated_ID = number;
	ptr_bv1->associated_sub_ID = 0;
	
	//BV on faces
	BoundingTriangularBox* ptr_bv2;
	for (int i = 0; i < ptr_cad->n_faces; i++)
	{
		ptr_bv2 = static_cast<BoundingTriangularBox*>(sub_bv[i]);
		//Test of small radius
		float eff_radius = deformed_faces_radii[i];
		if (eff_radius < (1.0f + inc_len_factor)*(float)largest_gnb)
			eff_radius = (1.0f + inc_len_factor)*(float)largest_gnb;
		//end of Test of small radius
		ptr_bv2->thickness = inc_len_factor * eff_radius;
		ptr_bv2->inc_len_factor = inc_len_factor;
		float templen = (1.0f + inc_len_factor)*eff_radius;
		ptr_bv2->size = sqrt(templen *templen + ptr_bv2->thickness*ptr_bv2->thickness);
		//Associated entity:
		ptr_bv2->associated_type = 'P';
		ptr_bv2->associated_ID = number;
		ptr_bv2->associated_sub_ID = i + 1;
	}
	//BV on edges
	BoundingCylinder* ptr_bv3;
	for (int i = ptr_cad->n_faces; i < ptr_cad->n_faces + (int)ptr_cad->edges.size(); i++)
	{
		//Determining the radius of the cylinder (smallest radius related to a neighboring face)
		float temp_radius = FLT_MAX;
		for (int j = 0; j < ptr_cad->edges[i - ptr_cad->n_faces].faceIDs.size(); j++)
		{
			if (deformed_faces_radii[ptr_cad->edges[i - ptr_cad->n_faces].faceIDs[j] - 1] < temp_radius)
				temp_radius = deformed_faces_radii[ptr_cad->edges[i - ptr_cad->n_faces].faceIDs[j] - 1];
		}
		//Test of small radius
		if (temp_radius < (1.0f + inc_len_factor)*(float)largest_gnb)
			temp_radius = (1.0f + inc_len_factor)*(float)largest_gnb;
		//end of Test of small radius
		ptr_bv3 = static_cast<BoundingCylinder*>(sub_bv[i]);
		ptr_bv3->radius = inc_len_factor * temp_radius;
		ptr_bv3->inc_len_factor = inc_len_factor;
		float templen = 0.5f*(1.0f + inc_len_factor)*deformed_edges_lengths[i - ptr_cad->n_faces] + ptr_bv3->radius;
		ptr_bv3->size = sqrt(templen*templen + ptr_bv3->radius*ptr_bv3->radius);
		//Associated entity:
		ptr_bv3->associated_type = 'P';
		ptr_bv3->associated_ID = number;
		ptr_bv3->associated_sub_ID = i + 1;

	}
	//BV on vertices
	BoundingSphere* ptr_bv4;
	for (int i = ptr_cad->n_faces + (int)ptr_cad->edges.size(); i < ptr_cad->n_faces + (int)ptr_cad->edges.size() + (int)ptr_cad->vertices.size(); i++)
	{
		//Determining the radius of the cylinder (smallest radius related to a neighboring face)
		float temp_radius = FLT_MAX;
		for (int j = 0; j < ptr_cad->vertices[i - ptr_cad->n_faces - (int)ptr_cad->edges.size()].faceIDs.size(); j++)
		{
			if (deformed_faces_radii[ptr_cad->vertices[i - ptr_cad->n_faces - (int)ptr_cad->edges.size()].faceIDs[j] - 1] < temp_radius)
				temp_radius = deformed_faces_radii[ptr_cad->vertices[i - ptr_cad->n_faces - (int)ptr_cad->edges.size()].faceIDs[j] - 1];
		}
		//Test of small radius
		if (temp_radius < (1.0f + inc_len_factor)*(float)largest_gnb)
			temp_radius = (1.0f + inc_len_factor)*(float)largest_gnb;
		//end of Test of small radius
		ptr_bv4 = static_cast<BoundingSphere*>(sub_bv[i]);
		ptr_bv4->radius = inc_len_factor * temp_radius;
		ptr_bv4->size = ptr_bv4->radius;
		//Associated entity:
		ptr_bv4->associated_type = 'P';
		ptr_bv4->associated_ID = number;
		ptr_bv4->associated_sub_ID = i + 1;
	}
}

void VEMPolyhedron::WriteVTK_XMLRender(FILE *f)
{
	//PostProcessing();
	//db.PrintPtr(elementPostprocessing, 26);

	STLSurface* ptr_cad = static_cast<STLSurface*>(db.cad_data[CADDATA_ID - 1]);
	std::vector<int> vint;
	
	//Abre piece
	fprintf(f, "\n<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n<Points>\n<DataArray type = \"Float32\" NumberOfComponents = \"3\" format = \"binary\">\n", n_vertices, ptr_cad->n_faces);
	std::vector<float> coord2;
	std::vector<float> float_vector;

	for (int i = 0; i < n_vertices; i++)
	{
		coord2.push_back((float)(db.super_nodes[super_node - 1]->copy_coordinates[3 * i + 0]));
		coord2.push_back((float)(db.super_nodes[super_node - 1]->copy_coordinates[3 * i + 1]));
		coord2.push_back((float)(db.super_nodes[super_node - 1]->copy_coordinates[3 * i + 2]));
	}

	fprintf(f, encodeData<float>(coord2).c_str());
	fprintf(f, "\n</DataArray>\n</Points>");

	fprintf(f, "\n<Cells><DataArray type=\"Int32\" Name=\"connectivity\" format=\"binary\">\n");
	for (int i = 0; i < ptr_cad->n_faces; i++)
	{
		vint.push_back(ptr_cad->faces[i]->verticesIDs[0] - 1);
		vint.push_back(ptr_cad->faces[i]->verticesIDs[1] - 1);
		vint.push_back(ptr_cad->faces[i]->verticesIDs[2] - 1);
	}

	fprintf(f, encodeData<int>(vint).c_str());

	vint.clear();
	fprintf(f, "\n</DataArray>\n<DataArray type=\"Int32\" Name=\"offsets\" format=\"binary\">\n");
	int cur_offset = 0;
	for (int i = 0; i < ptr_cad->n_faces; i++)
	{
		cur_offset = cur_offset + 3;
		vint.push_back(cur_offset);
	}
	fprintf(f, encodeData<int>(vint).c_str());

	vint.clear();
	fprintf(f, "\n</DataArray>\n<DataArray type=\"Int32\" Name=\"types\" format=\"binary\">\n");
	for (int i = 0; i < ptr_cad->n_faces; i++)
	{
		vint.push_back(5);
	}
	fprintf(f, encodeData<int>(vint).c_str());

	vint.clear();
	fprintf(f, "</DataArray>\n</Cells>\n");

	//Opens PointData
	fprintf(f, "\t\t\t<PointData Vectors = \"Displacement Velocity Acceleration\">\n");
	float_vector.clear();
	//Opens DataArray
	fprintf(f, "\t\t\t\t<DataArray Name = \"Displacement\" type = \"Float32\" NumberOfComponents= \"3\" format=\"binary\">\n");
	for (int n = 0; n < (int)ptr_cad->vertices.size(); n++)
	{
		float_vector.push_back((float)(db.super_nodes[super_node - 1]->copy_coordinates[3 * n + 0] - db.super_nodes[super_node - 1]->ref_coordinates[3 * n + 0]));
		float_vector.push_back((float)(db.super_nodes[super_node - 1]->copy_coordinates[3 * n + 1] - db.super_nodes[super_node - 1]->ref_coordinates[3 * n + 1]));
		float_vector.push_back((float)(db.super_nodes[super_node - 1]->copy_coordinates[3 * n + 2] - db.super_nodes[super_node - 1]->ref_coordinates[3 * n + 2]));

	}
	fprintf(f, encodeData<float>(float_vector).c_str());
	fprintf(f, "\n");
	//Closes DataArray
	fprintf(f, "\t\t\t\t</DataArray>\n");
	float_vector.clear();
	//Opens DataArray
	fprintf(f, "\t\t\t\t<DataArray Name = \"Velocity\" type = \"Float32\" NumberOfComponents= \"3\" format=\"binary\">\n");
	for (int n = 0; n < (int)ptr_cad->vertices.size(); n++)
	{
		float_vector.push_back((float)(db.super_nodes[super_node - 1]->copy_vel[3 * n + 0]));
		float_vector.push_back((float)(db.super_nodes[super_node - 1]->copy_vel[3 * n + 1]));
		float_vector.push_back((float)(db.super_nodes[super_node - 1]->copy_vel[3 * n + 2]));
	}
	fprintf(f, encodeData<float>(float_vector).c_str());
	fprintf(f, "\n");
	//Closes DataArray
	fprintf(f, "\t\t\t\t</DataArray>\n");
	float_vector.clear();
	//Opens DataArray
	fprintf(f, "\t\t\t\t<DataArray Name = \"Acceleration\" type = \"Float32\" NumberOfComponents= \"3\" format=\"binary\">\n");
	for (int n = 0; n < (int)ptr_cad->vertices.size(); n++)
	{
		float_vector.push_back((float)(db.super_nodes[super_node - 1]->copy_accel[3 * n + 0]));
		float_vector.push_back((float)(db.super_nodes[super_node - 1]->copy_accel[3 * n + 1]));
		float_vector.push_back((float)(db.super_nodes[super_node - 1]->copy_accel[3 * n + 2]));
	}
	fprintf(f, encodeData<float>(float_vector).c_str());
	fprintf(f, "\n");
	//Closes DataArray
	fprintf(f, "\t\t\t\t</DataArray>\n");
	//Closes PointData
	fprintf(f, "\t\t\t</PointData>\n");



	int res_element = 26;
	//Opens CellData
	fprintf(f, "<CellData FieldData = \"ElementData\">\n");
	//Opens DataArray
	fprintf(f, "<DataArray Name = \"ElementResults\" type = \"Float32\" NumberOfComponents=\"%d\" format=\"binary\">\n", res_element);
	//Imprime os resultados do elemento
	float_vector.clear();
	for (int i = 0; i < ptr_cad->n_faces; i++)
	{
		//For all the results are the same
		for (int r = 0; r < res_element; r++)
			float_vector.push_back((float)elementPostprocessing[r]);
	}
	
	fprintf(f, encodeData(float_vector).c_str());
	fprintf(f, "\n");
	//Closes DataArray
	fprintf(f, "</DataArray>\n");

	//Opens DataArray
	fprintf(f, "\t\t\t\t<DataArray Name=\"ParticleData\" type=\"Int32\" NumberOfComponents=\"%d\" format=\"binary\">\n", 4);
	//Imprime os dados da partícula
	vint.clear();
	for (int i = 0; i < ptr_cad->n_faces; i++)
	{
		vint.push_back(number);
		vint.push_back(super_node);
		vint.push_back(CADDATA_ID);
		vint.push_back(material);
	}
	fprintf(f, encodeData<int>(vint).c_str());
	fprintf(f, "\n");
	//Closes DataArray
	fprintf(f, "\t\t\t\t</DataArray>\n");
	//Closes CellData
	fprintf(f, "</CellData>\n");
	fprintf(f, "</Piece>\n");
}

void VEMPolyhedron::Alloc()
{
	Free();
	//Local variables for kinematics
	vertices_i = new Matrix*[n_vertices];
	vertices_p = new Matrix*[n_vertices];
	for (int i = 0; i < n_vertices; i++)
	{
		vertices_i[i] = new Matrix(3);
		vertices_p[i] = new Matrix(3);
	}
	deformed_faces_radii = new float[n_faces];
	for (int i = 0; i < n_faces; i++)
		deformed_faces_radii[i] = 0.0f;
	deformed_edges_lengths = new float[n_edges];
	for (int i = 0; i < n_edges; i++)
		deformed_edges_lengths[i] = 0.0f;
	/*count_local_abs_load = new int[3 * n_vertices];
	local_abs_load = new double[3 * n_vertices];*/
	local_residual = new double[3 * n_vertices];
	local_e_load = new double[3 * n_vertices];
	local_stiffness = new double*[3 * n_vertices];
	local_mass = new double*[3 * n_vertices];
	local_damping = new double*[3 * n_vertices];
	for (int i = 0; i < 3 * n_vertices; i++)
	{
		local_stiffness[i] = new double[3 * n_vertices];
		local_mass[i] = new double[3 * n_vertices];
		local_damping[i] = new double[3 * n_vertices];
		/*count_local_abs_load[i] = 0;
		local_abs_load[i] = 0.0;*/
		local_residual[i] = 0.0;
		local_e_load[i] = 0.0;
		for (int j = 0; j < 3 * n_vertices; j++)
		{
			local_stiffness[i][j] = 0.0;
			local_mass[i][j] = 0.0;
			local_damping[i][j] = 0.0;
		}
	}
	domainData = new double[12];
	//Filling values
	for (int i = 0; i < 12; i++)
		domainData[i] = 0.0;
	elementData = new double[20];
	//Filling values
	for (int i = 0; i < 20; i++)
		elementData[i] = 0.0;
	elementPostprocessing = new double[29];
	//Filling values
	for (int i = 0; i < 29; i++)
		elementPostprocessing[i] = 0.0;

	//flag
	allocced = true;
}

void VEMPolyhedron::Free()
{
	if (allocced == true)
	{
		if (vertices_i != NULL)
		{
			for (int i = 0; i < n_vertices; i++)
				delete vertices_i[i];
			delete[]vertices_i;
		}
		if (vertices_p != NULL)
		{
			for (int i = 0; i < n_vertices; i++)
				delete vertices_p[i];
			delete[]vertices_p;
		}
		/*if (count_local_abs_load != NULL)
		{
			delete[]count_local_abs_load;
		}
		if (local_abs_load != NULL)
		{
			delete[]local_abs_load;
		}*/
		if (local_residual != NULL)
		{
			delete[]local_residual;
		}
		if (local_e_load != NULL)
		{
			delete[]local_e_load;
		}
		if (local_stiffness != NULL)
		{
			for (int i = 0; i < 3 * n_vertices; i++)
				delete[]local_stiffness[i];
			delete[]local_stiffness;
		}
		if (local_mass != NULL)
		{
			for (int i = 0; i < 3 * n_vertices; i++)
				delete[]local_mass[i];
			delete[]local_mass;
		}
		if (local_damping != NULL)
		{
			for (int i = 0; i < 3 * n_vertices; i++)
				delete[]local_damping[i];
			delete[]local_damping;
		}
		if (triangles != NULL)
			delete[]triangles;
		if (tetrahedrons != NULL)
			delete[]tetrahedrons;
		if (domainData != NULL)
			delete []domainData;
		if (elementData != NULL)
			delete []elementData;
		if (elementPostprocessing != NULL)
			delete[]elementPostprocessing;
		if (deformed_faces_radii != NULL)
			delete[]deformed_faces_radii;
		if (deformed_edges_lengths != NULL)
			delete[]deformed_edges_lengths;
		//flag
		allocced = false;
	}
}

void VEMPolyhedron::UpdateVariables()
{
	for (int i = 0; i < n_vertices; i++)
	{
		//Local variables
		(*vertices_p[i])(0, 0) = db.super_nodes[super_node - 1]->copy_coordinates[3 * i] + db.super_nodes[super_node - 1]->displacements[3 * i];
		(*vertices_p[i])(1, 0) = db.super_nodes[super_node - 1]->copy_coordinates[3 * i + 1] + db.super_nodes[super_node - 1]->displacements[3 * i + 1];
		(*vertices_p[i])(2, 0) = db.super_nodes[super_node - 1]->copy_coordinates[3 * i + 2] + db.super_nodes[super_node - 1]->displacements[3 * i + 2];
	}
	//Updating the position of the centroid
	UpdateVolumeAndBarycenter();
	(*x0p)(0, 0) = (double)(*deformed_barycenter)(0, 0);
	(*x0p)(1, 0) = (double)(*deformed_barycenter)(1, 0);
	(*x0p)(2, 0) = (double)(*deformed_barycenter)(2, 0);
}

void VEMPolyhedron::MountGlobal()
{
	//Variáveis temporárias para salvar a indexação global dos graus de liberdade a serem setados na matriz de rigidez global
	int GL_global_1 = 0;
	int GL_global_2 = 0;
	double anterior = 0;
	for (int i = 0; i < db.super_nodes[super_node-1]->n_DOFs; i++)
	{
		GL_global_1 = db.super_nodes[super_node - 1]->DOFs[i];
		//Caso o grau de liberdade seja livre:
		if (GL_global_1 > 0)
		{
			anterior = db.global_P_A(GL_global_1 - 1, 0);
			db.global_P_A(GL_global_1 - 1, 0) = anterior + local_residual[i];
			anterior = db.global_I_A(GL_global_1 - 1, 0);
			db.global_I_A(GL_global_1 - 1, 0) = anterior + local_residual[i];
			//Convergence criteria
			/*anterior = db.global_ABS_P_A(GL_global_1 - 1, 0);
			db.global_ABS_P_A(GL_global_1 - 1, 0) = anterior + local_abs_load[i];
			db.global_COUNT_ABS_P_A(GL_global_1 - 1, 0) = db.global_COUNT_ABS_P_A(GL_global_1 - 1, 0) + count_local_abs_load[i];*/
		}
		//Caso o grau de liberdade seja fixo:
		if (GL_global_1 < 0)
		{
			anterior = db.global_P_B(-GL_global_1 - 1, 0);
			db.global_P_B(-GL_global_1 - 1, 0) = anterior + local_residual[i];
			//Convergence criteria
			/*anterior = db.global_ABS_P_B(-GL_global_1 - 1, 0);
			db.global_ABS_P_B(-GL_global_1 - 1, 0) = anterior + local_abs_load[i];
			db.global_COUNT_ABS_P_B(-GL_global_1 - 1, 0) = db.global_COUNT_ABS_P_B(-GL_global_1 - 1, 0) + count_local_abs_load[i];*/
		}
		for (int j = 0; j < db.super_nodes[super_node - 1]->n_DOFs; j++)
		{
			GL_global_2 = db.super_nodes[super_node - 1]->DOFs[j];
			
			//Caso os graus de liberdade sejam ambos livres (Matriz Kaa)
			if (GL_global_1 > 0 && GL_global_2 > 0)
				db.global_stiffness_AA.setValue(GL_global_1 - 1, GL_global_2 - 1, local_stiffness[i][j]);
			//Caso os graus de liberdade sejam ambos fixos (Matriz Kbb)
			if (GL_global_1 < 0 && GL_global_2 < 0)
				db.global_stiffness_BB.setValue(-GL_global_1 - 1, -GL_global_2 - 1, local_stiffness[i][j]);
			//Caso os graus de liberdade sejam livre e fixo (Matriz Kab)
			if (GL_global_1 > 0 && GL_global_2 < 0)
				db.global_stiffness_AB.setValue(GL_global_1 - 1, -GL_global_2 - 1, local_stiffness[i][j]);
			//Caso os graus de liberdade sejam fixo e livre (Matriz Kba)
			if (GL_global_1 < 0 && GL_global_2 > 0)
				db.global_stiffness_BA.setValue(-GL_global_1 - 1, GL_global_2 - 1, local_stiffness[i][j]);
		}
	}
}

void VEMPolyhedron::Mount()
{
	//elementData[0] = n_vertices;
	//elementData[1] = ptr_cad->n_faces;

	double *v;
	v = new double[3305 + 51 * (int)elementData[0] + 16 * (int)elementData[1] + 3 * (int)elementData[0] * (int)elementData[1] + 18 * (int)elementData[1] * (int)elementData[1]];
	for (int i = 0; i < 3 * n_vertices; i++)
	{
		/*count_local_abs_load[i] = 0;
		local_abs_load[i] = 0.0;*/
		local_residual[i] = 0.0;
		local_e_load[i] = 0.0;
		for (int j = 0; j < 3 * n_vertices; j++)
		{
			local_stiffness[i][j] = 0.0;
			local_mass[i][j] = 0.0;
		}
	}

	//Particle weight
	
	bool in_gcs_domain = true;
	l_factor = 0.0;
	if (db.environment_exist == true)
	{
		//Se existe campo gravitacional
		if (db.environment->g_exist == true)
		{
			if (db.gcs_exist)
			{
				//If the user has set boundings for contacts, consider here the same boundings for the gravitational field
				if (db.gcs->user_set_boundings)
				{
					double center[3];
					//teest is done considering the last converged position
					center[0] = db.super_nodes[super_node - 1]->super_node_coordinates[0];
					center[1] = db.super_nodes[super_node - 1]->super_node_coordinates[1];
					center[2] = db.super_nodes[super_node - 1]->super_node_coordinates[2];


					if (center[0] < db.gcs->min_x || center[0] > db.gcs->max_x ||
						center[1] < db.gcs->min_y || center[1] > db.gcs->max_y ||
						center[2] < db.gcs->min_z || center[2] > db.gcs->max_z)
					{
						in_gcs_domain = false;
					}

				}

			}
			//includes the weight only it the particle lies inside the domain of interest, defined by gcs
			if (in_gcs_domain)
			{
				l_factor = db.environment->bool_g.GetLinearFactorAtCurrentTime();
				domainData[3] = l_factor * db.environment->G(0, 0)*db.materials[material - 1]->rho;
				domainData[4] = l_factor * db.environment->G(1, 0)*db.materials[material - 1]->rho;
				domainData[5] = l_factor * db.environment->G(2, 0)*db.materials[material - 1]->rho;

			}
			else
			{
				domainData[3] = 0.0;
				domainData[4] = 0.0;
				domainData[5] = 0.0;
			}
		}
	}
	else
	{
		domainData[3] = 0.0;
		domainData[4] = 0.0;
		domainData[5] = 0.0;
	}

#pragma region
	int i2, i4, i6, i7, i9, i30, i32, i35, i37, i39, i47, i49, i53, i57, i61, i64, i66, i68, i70, i72
		, i74, i76, i78, i80, i237, i242, i257, i259, i260, i261, i264, i267, i279, i281, i283
		, i284, i286, i288, i290, i291, i293, i295, i297, i299, i301
		, i302, i304, i305, i307, i309, i311, i314, i316, i318, i328, i329, i333
		, i337, i340, i341, i345, i349, i352, i353, i357, i361, i364, i462, i464, i465, i478, i479
		, i481, i482, i484, i486, i487, i488, i490, i495, i497, i501, i505, i509, i548, i655, i694
		, i724, i727, i730, i761, i763, i764, i777, i778, i780, i781, i783, i785
		, i786, i787, i789, i790, i794, i795, i796, i799, i801, i803, i816, i818, i820, i822, i823
		, i824, i828, i832, i835, i836, i840, i844, i847, i848, i852, i856, i859, i1000, i1002, i1003
		, i1016, i1017, i1019, i1020, i1022, i1024, i1025, i1026, i1028, i1029, i1033, i1034, i1035
		, i1039, i1041, i1043, i1052, i1053, i1054, i1055, i1057, i1058, i1059, i1060, i1061, i1062
		, i1063, i1064, i1066, i1067, i1068, i1069, i1164, i1180, i1182, i1189, i1191, i1196, i1198
		, i1200, i1271, i1273, i1275, i1277, i1302, i1304, i1311, i1312, i1314, i1315, i1317, i1320
		, i1322, i1325, i1327, i1329, i1356, i1362, i1364, i1377, i1390, i1391, i1393, i1394, i1396
		, i1398, i1399, i1400, i1402, i1403, i1407, i1408, i1409, i1412, i1414, i1416, i1429, i1431
		, i1433, i1435, i1436, i1437, i1441, i1445, i1448, i1449, i1453, i1457, i1460, i1461, i1465
		, i1469, i1472, i1497, i1499, i1512, i1525, i1526, i1528, i1529, i1531, i1533, i1534, i1535
		, i1537, i1538, i1542, i1543, i1544, i1869, b238, i0;
	i0 = 247;
	i2 = (int)(elementData[0]);
	i1869 = 3 * i2;
	i462 = i0;
	i0 = i0 + i1869;
	i316 = i0;
	i0 = i0 + i1869;
	i257 = i0;
	i0 = i0 + i1869;
	i4 = (int)(elementData[1]);
	i6 = (int)(elementData[2]);
	i7 = i0;
	i0 = i0 + i1869;
	for (i9 = 1; i9 <= i1869; i9++) {
		v[2909 + i7 + i9] = displacement[i9 - 1];
	};/* end for */
	v[10] = 0e0;
	v[11] = 0e0;
	v[12] = 0e0;
	v[13] = 0e0;
	v[14] = 0e0;
	v[15] = 0e0;
	v[16] = 0e0;
	v[17] = 0e0;
	v[18] = 0e0;
	v[19] = 0e0;
	v[20] = 0e0;
	v[21] = 0e0;
	v[22] = 0e0;
	v[23] = 0e0;
	v[24] = 0e0;
	v[25] = 0e0;
	v[26] = 0e0;
	i724 = i0;
	i0 = i0 + i2;
	i727 = i0;
	i0 = i0 + i2;
	i730 = i0;
	i0 = i0 + i2;
	for (i30 = 1; i30 <= i2; i30++) {
		i32 = 3 * (-1 + i30);
		i39 = 3 + i32;
		v[2909 + i30 + i724] = (double)(i39);
		i37 = 2 + i32;
		v[2909 + i30 + i727] = (double)(i37);
		i35 = 1 + i32;
		v[2909 + i30 + i730] = (double)(i35);
		v[31] = Xref[i35 - 1];
		v[33] = Xref[i37 - 1];
		v[34] = Xref[i39 - 1];
		v[10] = v[10] + v[31];
		v[11] = v[11] + v[33];
		v[12] = v[12] + v[34];
		v[13] = v[13] + v[2909 + i35 + i7] - v[31] + Xcur[i35 - 1];
		v[14] = v[14] + v[2909 + i37 + i7] - v[33] + Xcur[i37 - 1];
		v[15] = v[15] + v[2909 + i39 + i7] - v[34] + Xcur[i39 - 1];
	};/* end for */
	i261 = i0;
	i0 = i0 + i4;
	i264 = i0;
	i0 = i0 + i4;
	i267 = i0;
	i0 = i0 + i4;
	i279 = i0;
	i0 = i0 + i4;
	i281 = i0;
	i0 = i0 + i4;
	i284 = i0;
	i0 = i0 + i4;
	i288 = i0;
	i0 = i0 + i4;
	i291 = i0;
	i0 = i0 + i4;
	i295 = i0;
	i0 = i0 + i4;
	i299 = i0;
	i0 = i0 + i4;
	i302 = i0;
	i0 = i0 + i4;
	i305 = i0;
	i0 = i0 + i4;
	i309 = i0;
	i0 = i0 + i4;
	for (i47 = 1; i47 <= i4; i47++) {
		i49 = 3 * (-1 + i47);
		i53 = 3 * (-1 + triangles[i49]);
		i68 = 3 + i53;
		v[2909 + i302 + i47] = (double)(i68);
		i66 = 2 + i53;
		v[2909 + i305 + i47] = (double)(i66);
		i64 = 1 + i53;
		v[2909 + i309 + i47] = (double)(i64);
		i57 = 3 * (-1 + triangles[1 + i49]);
		i74 = 3 + i57;
		v[2909 + i291 + i47] = (double)(i74);
		i72 = 2 + i57;
		v[2909 + i295 + i47] = (double)(i72);
		i70 = 1 + i57;
		v[2909 + i299 + i47] = (double)(i70);
		i61 = 3 * (-1 + triangles[2 + i49]);
		i80 = 3 + i61;
		v[2909 + i281 + i47] = (double)(i80);
		i78 = 2 + i61;
		v[2909 + i284 + i47] = (double)(i78);
		i76 = 1 + i61;
		v[2909 + i288 + i47] = (double)(i76);
		v[52] = Xref[i64 - 1];
		v[54] = Xref[i66 - 1];
		v[55] = Xref[i68 - 1];
		v[56] = Xref[i70 - 1];
		v[58] = Xref[i72 - 1];
		v[59] = Xref[i74 - 1];
		v[60] = Xref[i76 - 1];
		v[109] = v[56] - v[60];
		v[108] = v[52] - v[60];
		v[62] = Xref[i78 - 1];
		v[111] = v[58] - v[62];
		v[110] = v[54] - v[62];
		v[63] = Xref[i80 - 1];
		v[113] = v[59] - v[63];
		v[112] = v[55] - v[63];
		v[114] = -(v[111] * v[112]) + v[110] * v[113];
		v[115] = v[109] * v[112] - v[108] * v[113];
		v[116] = -(v[109] * v[110]) + v[108] * v[111];
		v[119] = sqrt((v[114] * v[114]) + (v[115] * v[115]) + (v[116] * v[116]));
		v[1870] = v[119] / 6e0;
		v[2909 + i279 + i47] = v[119];
		v[140] = v[1870] * (v[2909 + i68 + i7] + v[2909 + i7 + i74] + v[2909 + i7 + i80] - v[55] - v[59] - v[63] + Xcur[i68 - 1]
			+ Xcur[i74 - 1] + Xcur[i80 - 1]);
		v[136] = v[1870] * (v[2909 + i66 + i7] + v[2909 + i7 + i72] + v[2909 + i7 + i78] - v[54] - v[58] - v[62] + Xcur[i66 - 1]
			+ Xcur[i72 - 1] + Xcur[i78 - 1]);
		v[132] = v[1870] * (v[2909 + i64 + i7] + v[2909 + i7 + i70] + v[2909 + i7 + i76] - v[52] - v[56] - v[60] + Xcur[i64 - 1]
			+ Xcur[i70 - 1] + Xcur[i76 - 1]);
		v[126] = v[114] / v[119];
		v[2909 + i267 + i47] = v[126];
		v[127] = v[115] / v[119];
		v[2909 + i264 + i47] = v[127];
		v[128] = v[116] / v[119];
		v[2909 + i261 + i47] = v[128];
		v[26] = v[119] / 2e0 + v[26];
		v[25] = v[25] + (v[119] * (v[126] * v[52] + v[127] * v[54] + v[128] * v[55] + v[126] * v[56] + v[127] * v[58] + v[128] * v[59]
			+ v[126] * v[60] + v[127] * v[62] + v[128] * v[63])) / 18e0;
		v[16] = v[126] * v[132] + v[16];
		v[17] = v[127] * v[132] + v[17];
		v[18] = v[128] * v[132] + v[18];
		v[19] = v[126] * v[136] + v[19];
		v[20] = v[127] * v[136] + v[20];
		v[21] = v[128] * v[136] + v[21];
		v[22] = v[126] * v[140] + v[22];
		v[23] = v[127] * v[140] + v[23];
		v[24] = v[128] * v[140] + v[24];
	};/* end for */
	v[143] = v[16] / v[25];
	v[144] = v[17] / v[25];
	v[145] = v[18] / v[25];
	v[179] = -((-v[13] + v[10] * v[143] + v[11] * v[144] + v[12] * v[145]) / i2);
	v[146] = v[19] / v[25];
	v[147] = v[20] / v[25];
	v[148] = v[21] / v[25];
	v[181] = -((-v[14] + v[10] * v[146] + v[11] * v[147] + v[12] * v[148]) / i2);
	v[149] = v[22] / v[25];
	v[150] = v[23] / v[25];
	v[151] = v[24] / v[25];
	v[183] = -((v[10] * v[149] - v[15] + v[11] * v[150] + v[12] * v[151]) / i2);
	elementData[3] = v[179];
	elementData[4] = v[181];
	elementData[5] = v[183];
	elementData[6] = v[143];
	elementData[7] = v[146];
	elementData[8] = v[149];
	elementData[9] = v[144];
	elementData[10] = v[147];
	elementData[11] = v[150];
	elementData[12] = v[145];
	elementData[13] = v[148];
	elementData[14] = v[151];
	elementData[15] = elementData[15];
	elementData[16] = elementData[16];
	elementData[17] = elementData[17];
	elementData[18] = v[25];
	elementData[19] = v[26];
	v[208] = 1e0 + v[143];
	v[212] = 1e0 + v[147];
	v[440] = -(v[144] * v[146]) + v[208] * v[212];
	v[431] = v[145] * v[146] - v[148] * v[208];
	v[422] = v[144] * v[148] - v[145] * v[212];
	v[413] = v[144] * v[149] - v[150] * v[208];
	v[386] = v[146] * v[150] - v[149] * v[212];
	v[216] = 1e0 + v[151];
	v[404] = -(v[145] * v[149]) + v[208] * v[216];
	v[395] = v[145] * v[150] - v[144] * v[216];
	v[377] = v[148] * v[149] - v[146] * v[216];
	v[368] = -(v[148] * v[150]) + v[212] * v[216];
	v[244] = (v[144] * v[144]) + (v[145] * v[145]) + (v[146] * v[146]) + (v[148] * v[148]) + (v[149] * v[149]) +
		(v[150] * v[150]) + (v[208] * v[208]) + (v[212] * v[212]) + (v[216] * v[216]);
	v[223] = v[208] * v[368] + v[144] * v[377] + v[145] * v[386];
	v[1872] = -1e0 + (v[223] * v[223]);
	v[450] = 1e0 / Power(v[223], 0.16666666666666669e1);
	v[246] = 1e0 / Power(v[223], 0.6666666666666666e0);
	v[224] = domainData[0];
	v[225] = domainData[1];
	v[226] = domainData[2];
	v[227] = domainData[3];
	v[228] = domainData[4];
	v[229] = domainData[5];
	v[230] = v[224] / (2e0 + 2e0*v[225]);
	v[1887] = v[230] * v[244];
	v[231] = v[224] / (3e0 - 6e0*v[225]);
	v[1886] = v[231] / 2e0;
	v[1982] = ((5e0 / 9e0)*v[1887]) / Power(v[223], 0.26666666666666666e1) + v[1886] * (1e0 + 1e0 / (v[223] * v[223]));
	v[233] = domainData[6];
	v[234] = domainData[7];
	v[235] = domainData[8];
	i237 = (int)(domainData[9]);
	b238 = i237 == 1;
	if (b238) {
		v[240] = 1e0 - v[233];
	}
	else {
		v[240] = 1e0;
	};
	v[1871] = v[240] * v[25];
	v[449] = v[1871] * v[230];
	v[452] = (-2e0 / 3e0)*v[449] * v[450];
	v[248] = v[246] * v[449];
	v[241] = (v[1871] * (2e0*v[230] * (-3e0 + v[244] * v[246]) + v[231] * (v[1872] - 2e0*log(v[223])))) / 4e0;
	v[245] = v[1871] * ((v[1872] * v[1886]) / v[223] - (v[1887] * v[450]) / 3e0);
	v[1876] = v[245] / v[25];
	v[247] = (v[216] * v[248] + v[245] * v[440]) / v[25];
	v[249] = (v[150] * v[248] + v[245] * v[431]) / v[25];
	v[250] = (v[149] * v[248] + v[245] * v[422]) / v[25];
	v[251] = (v[148] * v[248] + v[245] * v[413]) / v[25];
	v[252] = (v[212] * v[248] + v[245] * v[404]) / v[25];
	v[253] = (v[146] * v[248] + v[245] * v[395]) / v[25];
	v[254] = (v[145] * v[248] + v[245] * v[386]) / v[25];
	v[255] = (v[144] * v[248] + v[245] * v[377]) / v[25];
	v[256] = (v[208] * v[248] + v[245] * v[368]) / v[25];
	for (i259 = 1; i259 <= i1869; i259++) {
		v[2909 + i257 + i259] = 0e0;
	};/* end for */
	i329 = i0;
	i0 = i0 + i4;
	i341 = i0;
	i0 = i0 + i4;
	i353 = i0;
	i0 = i0 + i4;
	for (i260 = i4; i260 >= 1; i260 = -1 + i260) {
		v[308] = v[2909 + i260 + i279] / 6e0;
		v[2909 + i260 + i329] = v[308];
		v[2909 + i260 + i341] = v[308];
		v[2909 + i260 + i353] = v[308];
		v[274] = v[2909 + i260 + i267];
		v[272] = v[2909 + i260 + i264];
		v[270] = v[2909 + i260 + i261];
		v[263] = v[247] * v[270];
		v[266] = v[263] + v[249] * v[272];
		v[269] = v[266] + v[250] * v[274];
		v[1873] = v[269] * v[308];
		v[271] = v[251] * v[270];
		v[273] = v[271] + v[252] * v[272];
		v[275] = v[273] + v[253] * v[274];
		v[1874] = v[275] * v[308];
		v[276] = v[254] * v[270];
		v[277] = v[255] * v[272] + v[276];
		v[278] = v[256] * v[274] + v[277];
		v[1875] = v[278] * v[308];
		i283 = (int)(v[2909 + i260 + i281]);
		v[2909 + i257 + i283] = v[1873] + v[2909 + i257 + i283];
		i286 = (int)(v[2909 + i260 + i284]);
		v[2909 + i257 + i286] = v[1874] + v[2909 + i257 + i286];
		i290 = (int)(v[2909 + i260 + i288]);
		v[2909 + i257 + i290] = v[1875] + v[2909 + i257 + i290];
		i293 = (int)(v[2909 + i260 + i291]);
		v[2909 + i257 + i293] = v[1873] + v[2909 + i257 + i293];
		i297 = (int)(v[2909 + i260 + i295]);
		v[2909 + i257 + i297] = v[1874] + v[2909 + i257 + i297];
		i301 = (int)(v[2909 + i260 + i299]);
		v[2909 + i257 + i301] = v[1875] + v[2909 + i257 + i301];
		i304 = (int)(v[2909 + i260 + i302]);
		v[2909 + i257 + i304] = v[1873] + v[2909 + i257 + i304];
		i307 = (int)(v[2909 + i260 + i305]);
		v[2909 + i257 + i307] = v[1874] + v[2909 + i257 + i307];
		i311 = (int)(v[2909 + i260 + i309]);
		v[2909 + i257 + i311] = v[1875] + v[2909 + i257 + i311];
	};/* end for */
	v[247] = 0e0;
	v[249] = 0e0;
	v[250] = 0e0;
	v[251] = 0e0;
	v[252] = 0e0;
	v[253] = 0e0;
	v[254] = 0e0;
	v[255] = 0e0;
	v[256] = 0e0;
	for (i242 = 1; i242 <= i1869; i242++) {
		for (i318 = 1; i318 <= i1869; i318++) {
			v[2909 + i316 + i318] = 0e0;
		};/* end for */
		v[2909 + i242 + i316] = 1e0 + v[2909 + i242 + i316];
		v[319] = 0e0;
		v[320] = 0e0;
		v[321] = 0e0;
		v[322] = 0e0;
		v[323] = 0e0;
		v[324] = 0e0;
		v[325] = 0e0;
		v[326] = 0e0;
		v[327] = 0e0;
		for (i328 = 1; i328 <= i4; i328++) {
			v[367] = v[2909 + i261 + i328];
			v[366] = v[2909 + i264 + i328];
			v[365] = v[2909 + i267 + i328];
			i364 = (int)(v[2909 + i281 + i328]);
			i361 = (int)(v[2909 + i284 + i328]);
			v[358] = v[2909 + i328 + i353];
			i357 = (int)(v[2909 + i288 + i328]);
			i352 = (int)(v[2909 + i291 + i328]);
			i349 = (int)(v[2909 + i295 + i328]);
			v[346] = v[2909 + i328 + i341];
			i345 = (int)(v[2909 + i299 + i328]);
			i340 = (int)(v[2909 + i302 + i328]);
			i337 = (int)(v[2909 + i305 + i328]);
			v[334] = v[2909 + i328 + i329];
			i333 = (int)(v[2909 + i309 + i328]);
			v[332] = v[2909 + i316 + i333];
			v[2909 + i316 + i333] = v[332];
			v[336] = v[2909 + i316 + i337];
			v[2909 + i316 + i337] = v[336];
			v[339] = v[2909 + i316 + i340];
			v[2909 + i316 + i340] = v[339];
			v[344] = v[2909 + i316 + i345];
			v[2909 + i316 + i345] = v[344];
			v[348] = v[2909 + i316 + i349];
			v[2909 + i316 + i349] = v[348];
			v[351] = v[2909 + i316 + i352];
			v[2909 + i316 + i352] = v[351];
			v[356] = v[2909 + i316 + i357];
			v[355] = v[332] * v[334] + v[344] * v[346] + v[356] * v[358];
			v[2909 + i316 + i357] = v[356];
			v[360] = v[2909 + i316 + i361];
			v[359] = v[334] * v[336] + v[346] * v[348] + v[358] * v[360];
			v[2909 + i316 + i361] = v[360];
			v[363] = v[2909 + i316 + i364];
			v[362] = v[334] * v[339] + v[346] * v[351] + v[358] * v[363];
			v[2909 + i316 + i364] = v[363];
			v[319] = v[319] + v[355] * v[365];
			v[320] = v[320] + v[355] * v[366];
			v[321] = v[321] + v[355] * v[367];
			v[322] = v[322] + v[359] * v[365];
			v[323] = v[323] + v[359] * v[366];
			v[324] = v[324] + v[359] * v[367];
			v[325] = v[325] + v[362] * v[365];
			v[326] = v[326] + v[362] * v[366];
			v[327] = v[327] + v[362] * v[367];
		};/* end for */
		v[1885] = v[327] / v[25];
		v[445] = v[1876] * v[327];
		v[1884] = v[326] / v[25];
		v[436] = v[1876] * v[326];
		v[1883] = v[325] / v[25];
		v[427] = v[1876] * v[325];
		v[1882] = v[324] / v[25];
		v[417] = v[1876] * v[324];
		v[1881] = v[323] / v[25];
		v[408] = v[1876] * v[323];
		v[1880] = v[322] / v[25];
		v[399] = v[1876] * v[322];
		v[1879] = v[321] / v[25];
		v[390] = v[1876] * v[321];
		v[1878] = v[320] / v[25];
		v[381] = v[1876] * v[320];
		v[1877] = v[319] / v[25];
		v[372] = v[1876] * v[319];
		v[319] = 0e0;
		v[320] = 0e0;
		v[321] = 0e0;
		v[322] = 0e0;
		v[323] = 0e0;
		v[324] = 0e0;
		v[325] = 0e0;
		v[326] = 0e0;
		v[441] = v[1877] * v[368] + v[1878] * v[377] + v[1879] * v[386] + v[1880] * v[395] + v[1881] * v[404] + v[1882] * v[413]
			+ v[1883] * v[422] + v[1884] * v[431] + v[1885] * v[440];
		v[1888] = v[441] * v[452];
		v[327] = 0e0;
		v[451] = v[1871] * v[1982] * v[441] + (v[144] * v[1878] + v[145] * v[1879] + v[146] * v[1880] + v[148] * v[1882]
			+ v[149] * v[1883] + v[150] * v[1884] + v[1877] * v[208] + v[1881] * v[212] + v[1885] * v[216])*v[452];
		v[453] = (v[1888] * v[216] + v[1885] * v[248] + v[212] * v[372] - v[146] * v[381] - v[144] * v[399] + v[208] * v[408]
			+ v[440] * v[451]) / v[25];
		v[454] = (v[150] * v[1888] + v[1884] * v[248] - v[148] * v[372] + v[146] * v[390] + v[145] * v[399] - v[208] * v[417]
			+ v[431] * v[451]) / v[25];
		v[455] = (v[149] * v[1888] + v[1883] * v[248] + v[148] * v[381] - v[212] * v[390] - v[145] * v[408] + v[144] * v[417]
			+ v[422] * v[451]) / v[25];
		v[456] = (v[148] * v[1888] + v[1882] * v[248] - v[150] * v[372] + v[149] * v[381] + v[144] * v[427] - v[208] * v[436]
			+ v[413] * v[451]) / v[25];
		v[457] = (v[1888] * v[212] + v[1881] * v[248] + v[216] * v[372] - v[149] * v[390] - v[145] * v[427] + v[208] * v[445]
			+ v[404] * v[451]) / v[25];
		v[458] = (v[146] * v[1888] + v[1880] * v[248] - v[216] * v[381] + v[150] * v[390] + v[145] * v[436] - v[144] * v[445]
			+ v[395] * v[451]) / v[25];
		v[459] = (v[145] * v[1888] + v[1879] * v[248] + v[150] * v[399] - v[149] * v[408] - v[212] * v[427] + v[146] * v[436]
			+ v[386] * v[451]) / v[25];
		v[460] = (v[144] * v[1888] + v[1878] * v[248] - v[216] * v[399] + v[149] * v[417] + v[148] * v[427] - v[146] * v[445]
			+ v[377] * v[451]) / v[25];
		v[461] = (v[1888] * v[208] + v[1877] * v[248] + v[216] * v[408] - v[150] * v[417] - v[148] * v[436] + v[212] * v[445]
			+ v[368] * v[451]) / v[25];
		for (i464 = 1; i464 <= i1869; i464++) {
			v[2909 + i462 + i464] = 0e0;
		};/* end for */
		for (i465 = i4; i465 >= 1; i465 = -1 + i465) {
			v[489] = v[2909 + i279 + i465] / 6e0;
			v[473] = v[2909 + i267 + i465];
			v[471] = v[2909 + i264 + i465];
			v[469] = v[2909 + i261 + i465];
			v[466] = v[453] * v[469];
			v[467] = v[466] + v[454] * v[471];
			v[468] = v[467] + v[455] * v[473];
			v[1889] = v[468] * v[489];
			v[470] = v[456] * v[469];
			v[472] = v[470] + v[457] * v[471];
			v[474] = v[472] + v[458] * v[473];
			v[1890] = v[474] * v[489];
			v[475] = v[459] * v[469];
			v[476] = v[460] * v[471] + v[475];
			v[477] = v[461] * v[473] + v[476];
			v[1891] = v[477] * v[489];
			i478 = (int)(v[2909 + i281 + i465]);
			v[2909 + i462 + i478] = v[1889] + v[2909 + i462 + i478];
			i479 = (int)(v[2909 + i284 + i465]);
			v[2909 + i462 + i479] = v[1890] + v[2909 + i462 + i479];
			i481 = (int)(v[2909 + i288 + i465]);
			v[2909 + i462 + i481] = v[1891] + v[2909 + i462 + i481];
			i482 = (int)(v[2909 + i291 + i465]);
			v[2909 + i462 + i482] = v[1889] + v[2909 + i462 + i482];
			i484 = (int)(v[2909 + i295 + i465]);
			v[2909 + i462 + i484] = v[1890] + v[2909 + i462 + i484];
			i486 = (int)(v[2909 + i299 + i465]);
			v[2909 + i462 + i486] = v[1891] + v[2909 + i462 + i486];
			i487 = (int)(v[2909 + i302 + i465]);
			v[2909 + i462 + i487] = v[1889] + v[2909 + i462 + i487];
			i488 = (int)(v[2909 + i305 + i465]);
			v[2909 + i462 + i488] = v[1890] + v[2909 + i462 + i488];
			i490 = (int)(v[2909 + i309 + i465]);
			v[2909 + i462 + i490] = v[1891] + v[2909 + i462 + i490];
		};/* end for */
		v[453] = 0e0;
		v[454] = 0e0;
		v[455] = 0e0;
		v[456] = 0e0;
		v[457] = 0e0;
		v[458] = 0e0;
		v[459] = 0e0;
		v[460] = 0e0;
		v[461] = 0e0;
		residual[i242 - 1] += v[2909 + i242 + i257];
		for (i314 = i242; i314 <= i1869; i314++) {
			tangent[i242 - 1][i314 - 1] += v[2909 + i314 + i462];
		};/* end for */
	};/* end for */
	if (v[234] != 1e0) {
		v[1930] = v[226] / 2e0;
		v[1901] = 2e0 / i2;
		i1000 = i0;
		i0 = i0 + i1869;
		i801 = i0;
		i0 = i0 + i1869;
		i761 = i0;
		i0 = i0 + i1869;
		v[734] = (b238 ? 1e0 - v[235] : 1e0);
		v[657] = (b238 ? 1e0 - v[234] : 1e0);
		for (i495 = 1; i495 <= i4; i495++) {
			i497 = 3 * (-1 + i495);
			i501 = 3 * (-1 + triangles[i497]);
			i505 = 3 * (-1 + triangles[1 + i497]);
			i509 = 3 * (-1 + triangles[2 + i497]);
			v[500] = Xref[i501];
			v[502] = Xref[1 + i501];
			v[503] = Xref[2 + i501];
			v[504] = Xref[i505];
			v[506] = Xref[1 + i505];
			v[507] = Xref[2 + i505];
			v[508] = Xref[i509];
			v[558] = v[504] - v[508];
			v[557] = v[500] - v[508];
			v[510] = Xref[1 + i509];
			v[560] = v[506] - v[510];
			v[559] = v[502] - v[510];
			v[565] = -(v[558] * v[559]) + v[557] * v[560];
			v[511] = Xref[2 + i509];
			v[562] = v[507] - v[511];
			v[561] = v[503] - v[511];
			v[564] = v[558] * v[561] - v[557] * v[562];
			v[563] = -(v[560] * v[561]) + v[559] * v[562];
			v[568] = sqrt((v[563] * v[563]) + (v[564] * v[564]) + (v[565] * v[565]));
			v[577] = v[565] / v[568];
			v[576] = v[564] / v[568];
			v[575] = v[563] / v[568];
			for (i548 = 1; i548 <= 6; i548++) {
				v[2965] = 0.4459484909159649e0;
				v[2966] = 0.10810301816807023e0;
				v[2967] = 0.4459484909159649e0;
				v[2968] = 0.9157621350977074e-1;
				v[2969] = 0.8168475729804585e0;
				v[2970] = 0.9157621350977074e-1;
				v[549] = v[2964 + i548];
				v[2971] = 0.4459484909159649e0;
				v[2972] = 0.4459484909159649e0;
				v[2973] = 0.10810301816807023e0;
				v[2974] = 0.9157621350977074e-1;
				v[2975] = 0.9157621350977074e-1;
				v[2976] = 0.8168475729804585e0;
				v[550] = v[2970 + i548];
				v[553] = 1e0 - v[549] - v[550];
				v[2983] = 0.11169079483900574e0;
				v[2984] = 0.11169079483900574e0;
				v[2985] = 0.11169079483900574e0;
				v[2986] = 0.5497587182766094e-1;
				v[2987] = 0.5497587182766094e-1;
				v[2988] = 0.5497587182766094e-1;
				v[578] = v[2982 + i548] * v[568];
				v[1892] = v[578] / 6e0;
				v[749] = -(v[1892] * v[734]);
				v[744] = v[229] * v[749];
				v[736] = 2e0*v[749];
				v[1910] = v[227] * v[736];
				v[1906] = v[228] * v[736];
				v[1904] = v[229] * v[736];
				v[670] = v[1892] * v[657];
				v[668] = (2e0 / 3e0)*v[578] * v[657];
				v[1896] = v[668] / v[25];
				v[1897] = v[1896] / 6e0;
				v[579] = v[500] * v[549] + v[504] * v[550] + v[508] * v[553];
				v[742] = v[575] * v[579];
				v[580] = v[502] * v[549] + v[506] * v[550] + v[510] * v[553];
				v[1934] = v[226] * v[580];
				v[740] = v[576] * v[580];
				v[581] = v[503] * v[549] + v[507] * v[550] + v[511] * v[553];
				v[1932] = v[226] * v[581];
				v[748] = v[581] * v[736] * (v[740] + v[742]);
				v[743] = v[577] * v[581];
				v[1902] = v[740] + v[743];
				v[1905] = (v[1902] + v[742]) / i2;
				v[1899] = v[579] * v[743];
				v[879] = v[742] + v[743];
				v[751] = v[580] * v[736] * v[879];
				v[594] = v[1930] * v[579];
				v[1898] = 6e0*v[575] * v[594];
				v[1894] = 2e0*v[594];
				v[599] = (v[579] * v[579]);
				v[1931] = v[575] * v[599];
				v[1895] = v[599] / 2e0;
				v[753] = v[1895] * v[575];
				v[1907] = v[1902] * v[579] + v[753];
				v[682] = 2e0*v[670] * v[753];
				v[613] = (v[580] * v[580]);
				v[1923] = v[594] * v[613];
				v[916] = (v[576] * v[613]) / 2e0;
				v[1917] = 2e0*v[916];
				v[1909] = v[751] + v[736] * v[916];
				v[616] = v[1894] * v[581];
				v[1933] = v[580] * v[616];
				v[621] = (v[581] * v[581]);
				v[1893] = v[577] * v[621];
				v[892] = v[1893] * v[594];
				v[885] = v[1893] / 2e0;
				v[1908] = v[748] + v[1893] * v[749];
				v[901] = v[577] * v[616] + v[1894] * v[740];
				v[912] = v[1896] * v[226] * (v[1895] * v[1902] + (v[575] * Power(v[579], 3)) / 6e0);
				v[923] = v[1897] * (v[1898] * v[613] + v[226] * (v[576] * Power(v[580], 3) + 3e0*v[613] * v[743]));
				v[934] = v[1897] * (v[1898] * v[621] + v[226] * (v[577] * Power(v[581], 3) + 3e0*v[621] * v[740]));
				v[658] = v[226] * v[670];
				v[947] = v[1901] * v[658] * v[916];
				v[659] = v[668] / 2e0;
				v[1915] = v[659] / v[25];
				v[1900] = v[580] * v[659];
				v[949] = v[1931] * v[658] + v[659] * v[901];
				v[905] = v[1899] * v[659] + v[682];
				v[898] = v[575] * v[616] * v[659];
				v[945] = v[1893] * v[658] + v[898];
				v[678] = v[1900] * v[226];
				v[1916] = v[678] / v[25];
				v[925] = v[1899] * v[678] + v[1894] * v[659] * v[916];
				v[674] = v[1900] * v[616];
				v[938] = v[576] * v[674] + v[659] * v[892];
				v[673] = v[1932] * v[659];
				v[936] = v[575] * v[674] + v[673] * v[916];
				v[932] = v[673] * v[740] + v[898];
				v[927] = (v[575] * v[673]) / v[25];
				v[943] = v[1901] * (v[1902] * v[658] + (v[1898] * v[659]) / 6e0);
				i824 = i0;
				i0 = i0 + i4;
				i836 = i0;
				i0 = i0 + i4;
				i848 = i0;
				i0 = i0 + i4;
				for (i694 = i4; i694 >= 1; i694 = -1 + i694) {
					v[1903] = v[2909 + i279 + i694] / 6e0;
					v[2909 + i694 + i824] = v[1903];
					v[2909 + i694 + i836] = v[1903];
					v[2909 + i694 + i848] = v[1903];
				};/* end for */
				v[739] = v[1904] * v[1905];
				v[741] = (-(v[12] * v[739]) + v[1893] * v[744] + v[229] * v[748]) / v[25];
				v[745] = (-(v[11] * v[739]) + v[1917] * v[744] + v[229] * v[751]) / v[25];
				v[746] = (v[1904] * v[1907] - v[10] * v[739]) / v[25];
				v[747] = v[1905] * v[1906];
				v[750] = (v[1908] * v[228] - v[12] * v[747]) / v[25];
				v[752] = (v[1909] * v[228] - v[11] * v[747]) / v[25];
				v[754] = (v[1906] * v[1907] - v[10] * v[747]) / v[25];
				v[755] = v[1905] * v[1910];
				v[757] = (v[1908] * v[227] - v[12] * v[755]) / v[25];
				v[759] = (v[1909] * v[227] - v[11] * v[755]) / v[25];
				v[760] = (v[1907] * v[1910] - v[10] * v[755]) / v[25];
				for (i763 = 1; i763 <= i1869; i763++) {
					v[2909 + i761 + i763] = 0e0;
				};/* end for */
				for (i764 = i4; i764 >= 1; i764 = -1 + i764) {
					v[788] = v[2909 + i279 + i764] / 6e0;
					v[772] = v[2909 + i267 + i764];
					v[770] = v[2909 + i264 + i764];
					v[768] = v[2909 + i261 + i764];
					v[765] = v[741] * v[768];
					v[766] = v[765] + v[745] * v[770];
					v[767] = v[766] + v[746] * v[772];
					v[1911] = v[767] * v[788];
					v[769] = v[750] * v[768];
					v[771] = v[769] + v[752] * v[770];
					v[773] = v[771] + v[754] * v[772];
					v[1912] = v[773] * v[788];
					v[774] = v[757] * v[768];
					v[775] = v[759] * v[770] + v[774];
					v[776] = v[760] * v[772] + v[775];
					v[1913] = v[776] * v[788];
					i777 = (int)(v[2909 + i281 + i764]);
					v[2909 + i761 + i777] = v[1911] + v[2909 + i761 + i777];
					i778 = (int)(v[2909 + i284 + i764]);
					v[2909 + i761 + i778] = v[1912] + v[2909 + i761 + i778];
					i780 = (int)(v[2909 + i288 + i764]);
					v[2909 + i761 + i780] = v[1913] + v[2909 + i761 + i780];
					i781 = (int)(v[2909 + i291 + i764]);
					v[2909 + i761 + i781] = v[1911] + v[2909 + i761 + i781];
					i783 = (int)(v[2909 + i295 + i764]);
					v[2909 + i761 + i783] = v[1912] + v[2909 + i761 + i783];
					i785 = (int)(v[2909 + i299 + i764]);
					v[2909 + i761 + i785] = v[1913] + v[2909 + i761 + i785];
					i786 = (int)(v[2909 + i302 + i764]);
					v[2909 + i761 + i786] = v[1911] + v[2909 + i761 + i786];
					i787 = (int)(v[2909 + i305 + i764]);
					v[2909 + i761 + i787] = v[1912] + v[2909 + i761 + i787];
					i789 = (int)(v[2909 + i309 + i764]);
					v[2909 + i761 + i789] = v[1913] + v[2909 + i761 + i789];
				};/* end for */
				for (i790 = i2; i790 >= 1; i790 = -1 + i790) {
					v[791] = v[739];
					v[792] = v[747];
					v[793] = v[755];
					i794 = (int)(v[2909 + i724 + i790]);
					v[2909 + i761 + i794] = v[2909 + i761 + i794] + v[791];
					i795 = (int)(v[2909 + i727 + i790]);
					v[2909 + i761 + i795] = v[2909 + i761 + i795] + v[792];
					i796 = (int)(v[2909 + i730 + i790]);
					v[2909 + i761 + i796] = v[2909 + i761 + i796] + v[793];
				};/* end for */
				v[741] = 0e0;
				v[745] = 0e0;
				v[746] = 0e0;
				v[750] = 0e0;
				v[752] = 0e0;
				v[754] = 0e0;
				v[757] = 0e0;
				v[759] = 0e0;
				v[760] = 0e0;
				v[739] = 0e0;
				v[747] = 0e0;
				v[755] = 0e0;
				for (i655 = 1; i655 <= i1869; i655++) {
					for (i803 = 1; i803 <= i1869; i803++) {
						v[2909 + i801 + i803] = 0e0;
					};/* end for */
					v[2909 + i655 + i801] = 1e0 + v[2909 + i655 + i801];
					v[804] = 0e0;
					v[805] = 0e0;
					v[806] = 0e0;
					v[807] = 0e0;
					v[808] = 0e0;
					v[809] = 0e0;
					v[810] = 0e0;
					v[811] = 0e0;
					v[812] = 0e0;
					v[813] = 0e0;
					v[814] = 0e0;
					v[815] = 0e0;
					for (i816 = 1; i816 <= i2; i816++) {
						i822 = (int)(v[2909 + i724 + i816]);
						i820 = (int)(v[2909 + i727 + i816]);
						i818 = (int)(v[2909 + i730 + i816]);
						v[817] = v[2909 + i801 + i818];
						v[2909 + i801 + i818] = v[817];
						v[819] = v[2909 + i801 + i820];
						v[2909 + i801 + i820] = v[819];
						v[821] = v[2909 + i801 + i822];
						v[2909 + i801 + i822] = v[821];
						v[804] = v[804] + v[817];
						v[805] = v[805] + v[819];
						v[806] = v[806] + v[821];
					};/* end for */
					for (i823 = 1; i823 <= i4; i823++) {
						v[862] = v[2909 + i261 + i823];
						v[861] = v[2909 + i264 + i823];
						v[860] = v[2909 + i267 + i823];
						i859 = (int)(v[2909 + i281 + i823]);
						i856 = (int)(v[2909 + i284 + i823]);
						v[853] = v[2909 + i823 + i848];
						i852 = (int)(v[2909 + i288 + i823]);
						i847 = (int)(v[2909 + i291 + i823]);
						i844 = (int)(v[2909 + i295 + i823]);
						v[841] = v[2909 + i823 + i836];
						i840 = (int)(v[2909 + i299 + i823]);
						i835 = (int)(v[2909 + i302 + i823]);
						i832 = (int)(v[2909 + i305 + i823]);
						v[829] = v[2909 + i823 + i824];
						i828 = (int)(v[2909 + i309 + i823]);
						v[827] = v[2909 + i801 + i828];
						v[2909 + i801 + i828] = v[827];
						v[831] = v[2909 + i801 + i832];
						v[2909 + i801 + i832] = v[831];
						v[834] = v[2909 + i801 + i835];
						v[2909 + i801 + i835] = v[834];
						v[839] = v[2909 + i801 + i840];
						v[2909 + i801 + i840] = v[839];
						v[843] = v[2909 + i801 + i844];
						v[2909 + i801 + i844] = v[843];
						v[846] = v[2909 + i801 + i847];
						v[2909 + i801 + i847] = v[846];
						v[851] = v[2909 + i801 + i852];
						v[850] = v[827] * v[829] + v[839] * v[841] + v[851] * v[853];
						v[2909 + i801 + i852] = v[851];
						v[855] = v[2909 + i801 + i856];
						v[854] = v[829] * v[831] + v[841] * v[843] + v[853] * v[855];
						v[2909 + i801 + i856] = v[855];
						v[858] = v[2909 + i801 + i859];
						v[857] = v[829] * v[834] + v[841] * v[846] + v[853] * v[858];
						v[2909 + i801 + i859] = v[858];
						v[807] = v[807] + v[850] * v[860];
						v[808] = v[808] + v[850] * v[861];
						v[809] = v[809] + v[850] * v[862];
						v[810] = v[810] + v[854] * v[860];
						v[811] = v[811] + v[854] * v[861];
						v[812] = v[812] + v[854] * v[862];
						v[813] = v[813] + v[857] * v[860];
						v[814] = v[814] + v[857] * v[861];
						v[815] = v[815] + v[857] * v[862];
					};/* end for */
					v[1928] = v[815] / v[25];
					v[972] = v[1915] * v[815];
					v[1914] = v[814] / v[25];
					v[965] = v[1914] * v[678];
					v[961] = v[1914] * v[659];
					v[1927] = v[813] / v[25];
					v[954] = v[1915] * v[813];
					v[1925] = v[812] / v[25];
					v[930] = v[1915] * v[812];
					v[1924] = v[811] / v[25];
					v[920] = v[1916] * v[811];
					v[915] = v[1915] * v[811];
					v[1922] = v[810] / v[25];
					v[907] = v[1915] * v[810];
					v[1920] = v[809] / v[25];
					v[887] = v[1915] * v[809];
					v[1919] = v[808] / v[25];
					v[878] = v[1916] * v[808];
					v[873] = v[1915] * v[808];
					v[876] = v[1917] * v[873];
					v[1918] = v[807] / v[25];
					v[866] = v[1915] * v[807];
					v[804] = -(v[10] * v[1918]) + v[804];
					v[863] = v[576] * v[866];
					v[871] = v[807] * v[912];
					v[807] = 0e0;
					v[804] = -(v[11] * v[1919]) + v[804];
					v[882] = v[1923] * v[863] + v[808] * v[923];
					v[808] = 0e0;
					v[804] = -(v[12] * v[1920]) + v[804];
					v[1921] = v[804] / i2;
					v[895] = v[1921] * v[659];
					v[884] = v[575] * v[878] + v[809] * v[927];
					v[890] = v[878] * v[885] + v[866] * v[892] + v[809] * v[934];
					v[809] = 0e0;
					v[897] = v[878] * v[879] + v[866] * v[901] + v[1920] * v[932] + v[804] * v[943];
					v[900] = v[882] + v[1920] * v[936] + v[804] * v[947];
					v[804] = 0e0;
					v[805] = -(v[10] * v[1922]) + v[805];
					v[903] = v[576] * v[907];
					v[913] = v[810] * v[912];
					v[810] = 0e0;
					v[805] = -(v[11] * v[1924]) + v[805];
					v[924] = v[1923] * v[903] + v[811] * v[923];
					v[811] = 0e0;
					v[805] = -(v[12] * v[1925]) + v[805];
					v[1926] = v[805] / i2;
					v[941] = v[1926] * v[659];
					v[928] = v[575] * v[920] + v[812] * v[927];
					v[935] = v[892] * v[907] + v[885] * v[920] + v[812] * v[934];
					v[812] = 0e0;
					v[944] = v[901] * v[907] + v[879] * v[920] + v[1925] * v[932] + v[805] * v[943];
					v[948] = v[924] + v[1925] * v[936] + v[805] * v[947];
					v[805] = 0e0;
					v[806] = -(v[10] * v[1927]) + v[806];
					v[951] = v[576] * v[954];
					v[959] = v[813] * v[912];
					v[813] = 0e0;
					v[806] = -(v[11] * v[1914]) + v[806];
					v[968] = v[814] * v[923] + v[1923] * v[951];
					v[814] = 0e0;
					v[806] = -(v[12] * v[1928]) + v[806];
					v[1929] = v[806] / i2;
					v[979] = v[1929] * v[659];
					v[970] = v[815] * v[927] + v[575] * v[965];
					v[975] = v[815] * v[934] + v[892] * v[954] + v[885] * v[965];
					v[815] = 0e0;
					v[981] = v[1928] * v[932] + v[806] * v[943] + v[901] * v[954] + v[879] * v[965];
					v[983] = v[1928] * v[936] + v[806] * v[947] + v[968];
					v[806] = 0e0;
					v[988] = (v[1930] * (v[1931] * v[954] + v[1917] * v[961] + v[1893] * v[972]) + v[981]) / i2;
					v[989] = (v[1929] * v[945] + v[1933] * (v[951] + v[575] * v[961]) + v[975] + v[1932] * (v[1927] * v[682]
						+ v[916] * v[961] + v[740] * v[979]) - v[12] * v[988]) / v[25];
					v[990] = (v[1934] * (v[1927] * v[905] + v[885] * v[972] + v[879] * v[979]) + v[983] - v[11] * v[988]) / v[25];
					v[991] = (v[1914] * v[925] + v[1928] * v[938] + v[1929] * v[949] + v[959] + v[1895] * v[970] - v[10] * v[988]) / v[25];
					v[992] = (v[1930] * (v[1931] * v[907] + v[1917] * v[915] + v[1893] * v[930]) + v[944]) / i2;
					v[993] = (v[1933] * (v[903] + v[575] * v[915]) + v[935] + v[1932] * (v[1922] * v[682] + v[915] * v[916]
						+ v[740] * v[941]) + v[1926] * v[945] - v[12] * v[992]) / v[25];
					v[994] = (v[1934] * (v[1922] * v[905] + v[885] * v[930] + v[879] * v[941]) + v[948] - v[11] * v[992]) / v[25];
					v[995] = (v[913] + v[1924] * v[925] + v[1895] * v[928] + v[1925] * v[938] + v[1926] * v[949] - v[10] * v[992]) / v[25];
					v[996] = (v[1930] * (v[1931] * v[866] + v[876] + v[1893] * v[887]) + v[897]) / i2;
					v[997] = (v[1933] * (v[863] + v[575] * v[873]) + v[890] + v[1932] * (v[1918] * v[682] + v[876] / 2e0 + v[740] * v[895])
						+ v[1921] * v[945] - v[12] * v[996]) / v[25];
					v[998] = (v[900] + v[1934] * (v[885] * v[887] + v[879] * v[895] + v[1918] * v[905]) - v[11] * v[996]) / v[25];
					v[999] = (v[871] + v[1895] * v[884] + v[1919] * v[925] + v[1920] * v[938] + v[1921] * v[949] - v[10] * v[996]) / v[25];
					for (i1002 = 1; i1002 <= i1869; i1002++) {
						v[2909 + i1000 + i1002] = 0e0;
					};/* end for */
					for (i1003 = i4; i1003 >= 1; i1003 = -1 + i1003) {
						v[1027] = v[2909 + i1003 + i279] / 6e0;
						v[1011] = v[2909 + i1003 + i267];
						v[1009] = v[2909 + i1003 + i264];
						v[1007] = v[2909 + i1003 + i261];
						v[1004] = v[1007] * v[989];
						v[1005] = v[1004] + v[1009] * v[990];
						v[1006] = v[1005] + v[1011] * v[991];
						v[1935] = v[1006] * v[1027];
						v[1008] = v[1007] * v[993];
						v[1010] = v[1008] + v[1009] * v[994];
						v[1012] = v[1010] + v[1011] * v[995];
						v[1936] = v[1012] * v[1027];
						v[1013] = v[1007] * v[997];
						v[1014] = v[1013] + v[1009] * v[998];
						v[1015] = v[1014] + v[1011] * v[999];
						v[1937] = v[1015] * v[1027];
						i1016 = (int)(v[2909 + i1003 + i281]);
						v[2909 + i1000 + i1016] = v[1935] + v[2909 + i1000 + i1016];
						i1017 = (int)(v[2909 + i1003 + i284]);
						v[2909 + i1000 + i1017] = v[1936] + v[2909 + i1000 + i1017];
						i1019 = (int)(v[2909 + i1003 + i288]);
						v[2909 + i1000 + i1019] = v[1937] + v[2909 + i1000 + i1019];
						i1020 = (int)(v[2909 + i1003 + i291]);
						v[2909 + i1000 + i1020] = v[1935] + v[2909 + i1000 + i1020];
						i1022 = (int)(v[2909 + i1003 + i295]);
						v[2909 + i1000 + i1022] = v[1936] + v[2909 + i1000 + i1022];
						i1024 = (int)(v[2909 + i1003 + i299]);
						v[2909 + i1000 + i1024] = v[1937] + v[2909 + i1000 + i1024];
						i1025 = (int)(v[2909 + i1003 + i302]);
						v[2909 + i1000 + i1025] = v[1935] + v[2909 + i1000 + i1025];
						i1026 = (int)(v[2909 + i1003 + i305]);
						v[2909 + i1000 + i1026] = v[1936] + v[2909 + i1000 + i1026];
						i1028 = (int)(v[2909 + i1003 + i309]);
						v[2909 + i1000 + i1028] = v[1937] + v[2909 + i1000 + i1028];
					};/* end for */
					for (i1029 = i2; i1029 >= 1; i1029 = -1 + i1029) {
						v[1030] = v[988];
						v[1031] = v[992];
						v[1032] = v[996];
						i1033 = (int)(v[2909 + i1029 + i724]);
						v[2909 + i1000 + i1033] = v[1030] + v[2909 + i1000 + i1033];
						i1034 = (int)(v[2909 + i1029 + i727]);
						v[2909 + i1000 + i1034] = v[1031] + v[2909 + i1000 + i1034];
						i1035 = (int)(v[2909 + i1029 + i730]);
						v[2909 + i1000 + i1035] = v[1032] + v[2909 + i1000 + i1035];
					};/* end for */
					v[989] = 0e0;
					v[990] = 0e0;
					v[991] = 0e0;
					v[993] = 0e0;
					v[994] = 0e0;
					v[995] = 0e0;
					v[997] = 0e0;
					v[998] = 0e0;
					v[999] = 0e0;
					v[988] = 0e0;
					v[992] = 0e0;
					v[996] = 0e0;
					residualload[i655 - 1] += v[2909 + i655 + i761];
					for (i799 = i655; i799 <= i1869; i799++) {
						massmatrix[i655 - 1][i799 - 1] += v[2909 + i1000 + i799];
					};/* end for */
				};/* end for */
			};/* end for */
		};/* end for */
	}
	else {
	};
	if (b238) {
		i1302 = i0;
		i0 = i0 + i1869;
		i1275 = i0;
		i0 = i0 + i1869;
		i1271 = i0;
		i0 = i0 + i1869;
		i1198 = i0;
		i0 = i0 + i1869;
		i1189 = i0;
		i0 = i0 + i1869;
		i1180 = i0;
		i0 = i0 + i1869;
		for (i1039 = 1; i1039 <= i6; i1039++) {
			i1043 = 4 * (-1 + i1039);
			i1063 = 3 * (-1 + tetrahedrons[3 + i1043]);
			i1069 = 3 + i1063;
			v[1084] = Xref[i1069 - 1];
			v[1126] = -v[1084] + v[2909 + i1069 + i7] + Xcur[i1069 - 1];
			i1064 = 2 + i1063;
			v[1079] = Xref[i1064 - 1];
			v[1119] = -v[1079] + v[2909 + i1064 + i7] + Xcur[i1064 - 1];
			i1055 = 1 + i1063;
			v[1074] = Xref[i1055 - 1];
			v[1112] = -v[1074] + v[2909 + i1055 + i7] + Xcur[i1055 - 1];
			i1061 = 3 * (-1 + tetrahedrons[2 + i1043]);
			i1068 = 3 + i1061;
			v[1083] = Xref[i1068 - 1];
			v[1125] = -v[1083] + v[2909 + i1068 + i7] + Xcur[i1068 - 1];
			v[1094] = v[1083] - v[1084];
			i1062 = 2 + i1061;
			v[1078] = Xref[i1062 - 1];
			v[1118] = -v[1078] + v[2909 + i1062 + i7] + Xcur[i1062 - 1];
			v[1091] = v[1078] - v[1079];
			i1054 = 1 + i1061;
			v[1073] = Xref[i1054 - 1];
			v[1111] = -v[1073] + v[2909 + i1054 + i7] + Xcur[i1054 - 1];
			v[1088] = v[1073] - v[1074];
			i1059 = 3 * (-1 + tetrahedrons[1 + i1043]);
			i1067 = 3 + i1059;
			v[1082] = Xref[i1067 - 1];
			v[1124] = -v[1082] + v[2909 + i1067 + i7] + Xcur[i1067 - 1];
			v[1093] = v[1082] - v[1084];
			i1060 = 2 + i1059;
			v[1077] = Xref[i1060 - 1];
			v[1117] = -v[1077] + v[2909 + i1060 + i7] + Xcur[i1060 - 1];
			v[1090] = v[1077] - v[1079];
			v[1941] = -(v[1091] * v[1093]) + v[1090] * v[1094];
			i1053 = 1 + i1059;
			v[1072] = Xref[i1053 - 1];
			v[1110] = -v[1072] + v[2909 + i1053 + i7] + Xcur[i1053 - 1];
			v[1087] = v[1072] - v[1074];
			i1057 = 3 * (-1 + tetrahedrons[i1043]);
			i1066 = 3 + i1057;
			v[1081] = Xref[i1066 - 1];
			v[1123] = -v[1081] + v[2909 + i1066 + i7] + Xcur[i1066 - 1];
			v[1092] = v[1081] - v[1084];
			i1058 = 2 + i1057;
			v[1076] = Xref[i1058 - 1];
			v[1116] = -v[1076] + v[2909 + i1058 + i7] + Xcur[i1058 - 1];
			v[1089] = v[1076] - v[1079];
			v[1940] = v[1091] * v[1092] - v[1089] * v[1094];
			v[1939] = -(v[1090] * v[1092]) + v[1089] * v[1093];
			i1052 = 1 + i1057;
			v[1071] = Xref[i1052 - 1];
			v[1109] = -v[1071] + v[2909 + i1052 + i7] + Xcur[i1052 - 1];
			v[1086] = v[1071] - v[1074];
			v[1105] = v[1088] * v[1939] + v[1087] * v[1940] + v[1086] * v[1941];
			v[1938] = v[1105] / 24e0;
			v[1192] = -(v[1938] * v[235]);
			v[1187] = v[1938] * v[226] * v[234];
			v[1166] = v[1938] * v[233];
			v[1942] = v[1166] * v[230];
			v[1104] = (-(v[1087] * v[1089]) + v[1086] * v[1090]) / v[1105];
			v[1103] = (v[1087] * v[1092] - v[1086] * v[1093]) / v[1105];
			v[1102] = v[1939] / v[1105];
			v[1101] = (v[1088] * v[1089] - v[1086] * v[1091]) / v[1105];
			v[1100] = (-(v[1088] * v[1092]) + v[1086] * v[1094]) / v[1105];
			v[1099] = v[1940] / v[1105];
			v[1098] = (-(v[1088] * v[1090]) + v[1087] * v[1091]) / v[1105];
			v[1108] = -v[1098] - v[1101] - v[1104];
			v[1152] = 1e0 + v[1098] * v[1123] + v[1101] * v[1124] + v[1104] * v[1125] + v[1108] * v[1126];
			v[1149] = v[1098] * v[1116] + v[1101] * v[1117] + v[1104] * v[1118] + v[1108] * v[1119];
			v[1146] = v[1098] * v[1109] + v[1101] * v[1110] + v[1104] * v[1111] + v[1108] * v[1112];
			v[1097] = (v[1088] * v[1093] - v[1087] * v[1094]) / v[1105];
			v[1107] = -v[1097] - v[1100] - v[1103];
			v[1151] = v[1097] * v[1123] + v[1100] * v[1124] + v[1103] * v[1125] + v[1107] * v[1126];
			v[1148] = 1e0 + v[1097] * v[1116] + v[1100] * v[1117] + v[1103] * v[1118] + v[1107] * v[1119];
			v[1257] = -(v[1149] * v[1151]) + v[1148] * v[1152];
			v[1145] = v[1097] * v[1109] + v[1100] * v[1110] + v[1103] * v[1111] + v[1107] * v[1112];
			v[1254] = v[1146] * v[1151] - v[1145] * v[1152];
			v[1251] = -(v[1146] * v[1148]) + v[1145] * v[1149];
			v[1096] = v[1941] / v[1105];
			v[1106] = -v[1096] - v[1099] - v[1102];
			v[1150] = v[1096] * v[1123] + v[1099] * v[1124] + v[1102] * v[1125] + v[1106] * v[1126];
			v[1147] = v[1096] * v[1116] + v[1099] * v[1117] + v[1102] * v[1118] + v[1106] * v[1119];
			v[1256] = v[1149] * v[1150] - v[1147] * v[1152];
			v[1255] = -(v[1148] * v[1150]) + v[1147] * v[1151];
			v[1144] = 1e0 + v[1096] * v[1109] + v[1099] * v[1110] + v[1102] * v[1111] + v[1106] * v[1112];
			v[1253] = -(v[1146] * v[1150]) + v[1144] * v[1152];
			v[1252] = v[1145] * v[1150] - v[1144] * v[1151];
			v[1250] = v[1146] * v[1147] - v[1144] * v[1149];
			v[1249] = -(v[1145] * v[1147]) + v[1144] * v[1148];
			v[1159] = v[1146] * v[1255] + v[1145] * v[1256] + v[1144] * v[1257];
			v[1259] = 1e0 / Power(v[1159], 0.16666666666666669e1);
			v[1261] = (-2e0 / 3e0)*v[1259] * v[1942];
			v[1171] = v[1942] / Power(v[1159], 0.6666666666666666e0);
			v[1947] = ((v[1144] * v[1144]) + (v[1145] * v[1145]) + (v[1146] * v[1146]) + (v[1147] * v[1147]) + (v[1148] * v[1148]
				) + (v[1149] * v[1149]) + (v[1150] * v[1150]) + (v[1151] * v[1151]) + (v[1152] * v[1152]))*v[230];
			v[1992] = (1e0 + 1e0 / (v[1159] * v[1159]))*v[1886] + ((5e0 / 9e0)*v[1947]) / Power(v[1159]
				, 0.26666666666666666e1);
			v[1168] = v[1166] * (((-1e0 + (v[1159] * v[1159]))*v[1886]) / v[1159] - (v[1259] * v[1947]) / 3e0);
			v[1170] = v[1144] * v[1171] + v[1168] * v[1257];
			v[1172] = v[1147] * v[1171] + v[1168] * v[1254];
			v[1173] = v[1150] * v[1171] + v[1168] * v[1251];
			v[1174] = v[1145] * v[1171] + v[1168] * v[1256];
			v[1175] = v[1148] * v[1171] + v[1168] * v[1253];
			v[1176] = v[1151] * v[1171] + v[1168] * v[1250];
			v[1177] = v[1146] * v[1171] + v[1168] * v[1255];
			v[1178] = v[1149] * v[1171] + v[1168] * v[1252];
			v[1179] = v[1152] * v[1171] + v[1168] * v[1249];
			for (i1182 = 1; i1182 <= i1869; i1182++) {
				v[2909 + i1180 + i1182] = 0e0;
			};/* end for */
			v[2909 + i1052 + i1180] = v[1096] * v[1170] + v[1097] * v[1174] + v[1098] * v[1177] + v[2909 + i1052 + i1180];
			v[2909 + i1058 + i1180] = v[1096] * v[1172] + v[1097] * v[1175] + v[1098] * v[1178] + v[2909 + i1058 + i1180];
			v[2909 + i1066 + i1180] = v[1096] * v[1173] + v[1097] * v[1176] + v[1098] * v[1179] + v[2909 + i1066 + i1180];
			v[2909 + i1053 + i1180] = v[1099] * v[1170] + v[1100] * v[1174] + v[1101] * v[1177] + v[2909 + i1053 + i1180];
			v[2909 + i1060 + i1180] = v[1099] * v[1172] + v[1100] * v[1175] + v[1101] * v[1178] + v[2909 + i1060 + i1180];
			v[2909 + i1067 + i1180] = v[1099] * v[1173] + v[1100] * v[1176] + v[1101] * v[1179] + v[2909 + i1067 + i1180];
			v[2909 + i1054 + i1180] = v[1102] * v[1170] + v[1103] * v[1174] + v[1104] * v[1177] + v[2909 + i1054 + i1180];
			v[2909 + i1062 + i1180] = v[1102] * v[1172] + v[1103] * v[1175] + v[1104] * v[1178] + v[2909 + i1062 + i1180];
			v[2909 + i1068 + i1180] = v[1102] * v[1173] + v[1103] * v[1176] + v[1104] * v[1179] + v[2909 + i1068 + i1180];
			v[2909 + i1055 + i1180] = v[1106] * v[1170] + v[1107] * v[1174] + v[1108] * v[1177] + v[2909 + i1055 + i1180];
			v[2909 + i1064 + i1180] = v[1106] * v[1172] + v[1107] * v[1175] + v[1108] * v[1178] + v[2909 + i1064 + i1180];
			v[2909 + i1069 + i1180] = v[1106] * v[1173] + v[1107] * v[1176] + v[1108] * v[1179] + v[2909 + i1069 + i1180];
			for (i1041 = 1; i1041 <= 4; i1041++) {
				v[3053] = 0.585410196624969e0;
				v[3054] = 0.138196601125011e0;
				v[3055] = 0.138196601125011e0;
				v[3056] = 0.138196601125011e0;
				v[1047] = v[3052 + i1041];
				v[1943] = v[1047] * v[1192];
				v[3057] = 0.138196601125011e0;
				v[3058] = 0.585410196624969e0;
				v[3059] = 0.138196601125011e0;
				v[3060] = 0.138196601125011e0;
				v[1048] = v[3056 + i1041];
				v[1944] = v[1048] * v[1192];
				v[3061] = 0.138196601125011e0;
				v[3062] = 0.138196601125011e0;
				v[3063] = 0.585410196624969e0;
				v[3064] = 0.138196601125011e0;
				v[1049] = v[3060 + i1041];
				v[1945] = v[1049] * v[1192];
				v[1051] = 1e0 - v[1047] - v[1048] - v[1049];
				v[1946] = v[1051] * v[1192];
				for (i1191 = 1; i1191 <= i1869; i1191++) {
					v[2909 + i1189 + i1191] = 0e0;
				};/* end for */
				v[2909 + i1052 + i1189] = v[1943] * v[227] + v[2909 + i1052 + i1189];
				v[2909 + i1058 + i1189] = v[1943] * v[228] + v[2909 + i1058 + i1189];
				v[2909 + i1066 + i1189] = v[1943] * v[229] + v[2909 + i1066 + i1189];
				v[2909 + i1053 + i1189] = v[1944] * v[227] + v[2909 + i1053 + i1189];
				v[2909 + i1060 + i1189] = v[1944] * v[228] + v[2909 + i1060 + i1189];
				v[2909 + i1067 + i1189] = v[1944] * v[229] + v[2909 + i1067 + i1189];
				v[2909 + i1054 + i1189] = v[1945] * v[227] + v[2909 + i1054 + i1189];
				v[2909 + i1062 + i1189] = v[1945] * v[228] + v[2909 + i1062 + i1189];
				v[2909 + i1068 + i1189] = v[1945] * v[229] + v[2909 + i1068 + i1189];
				v[2909 + i1055 + i1189] = v[1946] * v[227] + v[2909 + i1055 + i1189];
				v[2909 + i1064 + i1189] = v[1946] * v[228] + v[2909 + i1064 + i1189];
				v[2909 + i1069 + i1189] = v[1946] * v[229] + v[2909 + i1069 + i1189];
				for (i1164 = 1; i1164 <= i1869; i1164++) {
					for (i1200 = 1; i1200 <= i1869; i1200++) {
						v[2909 + i1198 + i1200] = 0e0;
					};/* end for */
					v[2909 + i1164 + i1198] = 1e0 + v[2909 + i1164 + i1198];
					v[1202] = v[2909 + i1069 + i1198];
					v[2909 + i1069 + i1198] = v[1202];
					v[1206] = v[2909 + i1064 + i1198];
					v[2909 + i1064 + i1198] = v[1206];
					v[1210] = v[2909 + i1055 + i1198];
					v[2909 + i1055 + i1198] = v[1210];
					v[1214] = v[2909 + i1068 + i1198];
					v[2909 + i1068 + i1198] = v[1214];
					v[1218] = v[2909 + i1062 + i1198];
					v[2909 + i1062 + i1198] = v[1218];
					v[1222] = v[2909 + i1054 + i1198];
					v[2909 + i1054 + i1198] = v[1222];
					v[1226] = v[2909 + i1067 + i1198];
					v[2909 + i1067 + i1198] = v[1226];
					v[1230] = v[2909 + i1060 + i1198];
					v[2909 + i1060 + i1198] = v[1230];
					v[1234] = v[2909 + i1053 + i1198];
					v[2909 + i1053 + i1198] = v[1234];
					v[1238] = v[2909 + i1066 + i1198];
					v[1237] = v[1108] * v[1202] + v[1104] * v[1214] + v[1101] * v[1226] + v[1098] * v[1238];
					v[1239] = v[1107] * v[1202] + v[1103] * v[1214] + v[1100] * v[1226] + v[1097] * v[1238];
					v[1240] = v[1106] * v[1202] + v[1102] * v[1214] + v[1099] * v[1226] + v[1096] * v[1238];
					v[2909 + i1066 + i1198] = v[1238];
					v[1242] = v[2909 + i1058 + i1198];
					v[1241] = v[1108] * v[1206] + v[1104] * v[1218] + v[1101] * v[1230] + v[1098] * v[1242];
					v[1243] = v[1107] * v[1206] + v[1103] * v[1218] + v[1100] * v[1230] + v[1097] * v[1242];
					v[1244] = v[1106] * v[1206] + v[1102] * v[1218] + v[1099] * v[1230] + v[1096] * v[1242];
					v[2909 + i1058 + i1198] = v[1242];
					v[1246] = v[2909 + i1052 + i1198];
					v[1245] = v[1108] * v[1210] + v[1104] * v[1222] + v[1101] * v[1234] + v[1098] * v[1246];
					v[1247] = v[1107] * v[1210] + v[1103] * v[1222] + v[1100] * v[1234] + v[1097] * v[1246];
					v[1248] = v[1106] * v[1210] + v[1102] * v[1222] + v[1099] * v[1234] + v[1096] * v[1246];
					v[2909 + i1052 + i1198] = v[1246];
					v[1258] = v[1237] * v[1249] + v[1239] * v[1250] + v[1240] * v[1251] + v[1241] * v[1252] + v[1243] * v[1253]
						+ v[1244] * v[1254] + v[1245] * v[1255] + v[1247] * v[1256] + v[1248] * v[1257];
					v[1948] = v[1258] * v[1261];
					v[1260] = (v[1152] * v[1237] + v[1151] * v[1239] + v[1150] * v[1240] + v[1149] * v[1241] + v[1148] * v[1243]
						+ v[1147] * v[1244] + v[1146] * v[1245] + v[1145] * v[1247] + v[1144] * v[1248])*v[1261]
						+ v[1166] * v[1258] * v[1992];
					v[1262] = v[1168] * (v[1148] * v[1237] - v[1149] * v[1239] - v[1151] * v[1241] + v[1152] * v[1243])
						+ v[1171] * v[1248] + v[1257] * v[1260] + v[1144] * v[1948];
					v[1263] = v[1171] * v[1244] + v[1168] * (-(v[1145] * v[1237]) + v[1146] * v[1239] + v[1151] * v[1245]
						- v[1152] * v[1247]) + v[1254] * v[1260] + v[1147] * v[1948];
					v[1264] = v[1171] * v[1240] + v[1168] * (v[1145] * v[1241] - v[1146] * v[1243] - v[1148] * v[1245]
						+ v[1149] * v[1247]) + v[1251] * v[1260] + v[1150] * v[1948];
					v[1265] = v[1168] * (-(v[1147] * v[1237]) + v[1149] * v[1240] + v[1150] * v[1241] - v[1152] * v[1244])
						+ v[1171] * v[1247] + v[1256] * v[1260] + v[1145] * v[1948];
					v[1266] = v[1171] * v[1243] + v[1168] * (v[1144] * v[1237] - v[1146] * v[1240] - v[1150] * v[1245]
						+ v[1152] * v[1248]) + v[1253] * v[1260] + v[1148] * v[1948];
					v[1267] = v[1171] * v[1239] + v[1168] * (-(v[1144] * v[1241]) + v[1146] * v[1244] + v[1147] * v[1245]
						- v[1149] * v[1248]) + v[1250] * v[1260] + v[1151] * v[1948];
					v[1268] = v[1168] * (v[1147] * v[1239] - v[1148] * v[1240] - v[1150] * v[1243] + v[1151] * v[1244])
						+ v[1171] * v[1245] + v[1255] * v[1260] + v[1146] * v[1948];
					v[1269] = v[1171] * v[1241] + v[1168] * (-(v[1144] * v[1239]) + v[1145] * v[1240] + v[1150] * v[1247]
						- v[1151] * v[1248]) + v[1252] * v[1260] + v[1149] * v[1948];
					v[1270] = v[1171] * v[1237] + v[1168] * (v[1144] * v[1243] - v[1145] * v[1244] - v[1147] * v[1247]
						+ v[1148] * v[1248]) + v[1249] * v[1260] + v[1152] * v[1948];
					for (i1273 = 1; i1273 <= i1869; i1273++) {
						v[2909 + i1271 + i1273] = 0e0;
					};/* end for */
					v[2909 + i1052 + i1271] = v[1096] * v[1262] + v[1097] * v[1265] + v[1098] * v[1268] + v[2909 + i1052 + i1271];
					v[2909 + i1058 + i1271] = v[1096] * v[1263] + v[1097] * v[1266] + v[1098] * v[1269] + v[2909 + i1058 + i1271];
					v[2909 + i1066 + i1271] = v[1096] * v[1264] + v[1097] * v[1267] + v[1098] * v[1270] + v[2909 + i1066 + i1271];
					v[2909 + i1053 + i1271] = v[1099] * v[1262] + v[1100] * v[1265] + v[1101] * v[1268] + v[2909 + i1053 + i1271];
					v[2909 + i1060 + i1271] = v[1099] * v[1263] + v[1100] * v[1266] + v[1101] * v[1269] + v[2909 + i1060 + i1271];
					v[2909 + i1067 + i1271] = v[1099] * v[1264] + v[1100] * v[1267] + v[1101] * v[1270] + v[2909 + i1067 + i1271];
					v[2909 + i1054 + i1271] = v[1102] * v[1262] + v[1103] * v[1265] + v[1104] * v[1268] + v[2909 + i1054 + i1271];
					v[2909 + i1062 + i1271] = v[1102] * v[1263] + v[1103] * v[1266] + v[1104] * v[1269] + v[2909 + i1062 + i1271];
					v[2909 + i1068 + i1271] = v[1102] * v[1264] + v[1103] * v[1267] + v[1104] * v[1270] + v[2909 + i1068 + i1271];
					v[2909 + i1055 + i1271] = v[1106] * v[1262] + v[1107] * v[1265] + v[1108] * v[1268] + v[2909 + i1055 + i1271];
					v[2909 + i1064 + i1271] = v[1106] * v[1263] + v[1107] * v[1266] + v[1108] * v[1269] + v[2909 + i1064 + i1271];
					v[2909 + i1069 + i1271] = v[1106] * v[1264] + v[1107] * v[1267] + v[1108] * v[1270] + v[2909 + i1069 + i1271];
					for (i1277 = 1; i1277 <= i1869; i1277++) {
						v[2909 + i1275 + i1277] = 0e0;
					};/* end for */
					v[2909 + i1164 + i1275] = 1e0 + v[2909 + i1164 + i1275];
					v[1279] = v[2909 + i1069 + i1275];
					v[2909 + i1069 + i1275] = v[1279];
					v[1281] = v[2909 + i1064 + i1275];
					v[2909 + i1064 + i1275] = v[1281];
					v[1283] = v[2909 + i1055 + i1275];
					v[2909 + i1055 + i1275] = v[1283];
					v[1285] = v[2909 + i1068 + i1275];
					v[2909 + i1068 + i1275] = v[1285];
					v[1287] = v[2909 + i1062 + i1275];
					v[2909 + i1062 + i1275] = v[1287];
					v[1289] = v[2909 + i1054 + i1275];
					v[2909 + i1054 + i1275] = v[1289];
					v[1291] = v[2909 + i1067 + i1275];
					v[2909 + i1067 + i1275] = v[1291];
					v[1293] = v[2909 + i1060 + i1275];
					v[2909 + i1060 + i1275] = v[1293];
					v[1295] = v[2909 + i1053 + i1275];
					v[2909 + i1053 + i1275] = v[1295];
					v[1297] = v[2909 + i1066 + i1275];
					v[1296] = v[1187] * (v[1051] * v[1279] + v[1049] * v[1285] + v[1048] * v[1291] + v[1047] * v[1297]);
					v[2909 + i1066 + i1275] = v[1297];
					v[1299] = v[2909 + i1058 + i1275];
					v[1298] = v[1187] * (v[1051] * v[1281] + v[1049] * v[1287] + v[1048] * v[1293] + v[1047] * v[1299]);
					v[2909 + i1058 + i1275] = v[1299];
					v[1301] = v[2909 + i1052 + i1275];
					v[1300] = v[1187] * (v[1051] * v[1283] + v[1049] * v[1289] + v[1048] * v[1295] + v[1047] * v[1301]);
					v[2909 + i1052 + i1275] = v[1301];
					for (i1304 = 1; i1304 <= i1869; i1304++) {
						v[2909 + i1302 + i1304] = 0e0;
					};/* end for */
					v[2909 + i1052 + i1302] = v[1047] * v[1300] + v[2909 + i1052 + i1302];
					v[2909 + i1058 + i1302] = v[1047] * v[1298] + v[2909 + i1058 + i1302];
					v[2909 + i1066 + i1302] = v[1047] * v[1296] + v[2909 + i1066 + i1302];
					v[2909 + i1053 + i1302] = v[1048] * v[1300] + v[2909 + i1053 + i1302];
					v[2909 + i1060 + i1302] = v[1048] * v[1298] + v[2909 + i1060 + i1302];
					v[2909 + i1067 + i1302] = v[1048] * v[1296] + v[2909 + i1067 + i1302];
					v[2909 + i1054 + i1302] = v[1049] * v[1300] + v[2909 + i1054 + i1302];
					v[2909 + i1062 + i1302] = v[1049] * v[1298] + v[2909 + i1062 + i1302];
					v[2909 + i1068 + i1302] = v[1049] * v[1296] + v[2909 + i1068 + i1302];
					v[2909 + i1055 + i1302] = v[1051] * v[1300] + v[2909 + i1055 + i1302];
					v[2909 + i1064 + i1302] = v[1051] * v[1298] + v[2909 + i1064 + i1302];
					v[2909 + i1069 + i1302] = v[1051] * v[1296] + v[2909 + i1069 + i1302];
					residual[i1164 - 1] += v[2909 + i1164 + i1180];
					residualload[i1164 - 1] += v[2909 + i1164 + i1189];
					for (i1196 = i1164; i1196 <= i1869; i1196++) {
						tangent[i1164 - 1][i1196 - 1] += v[2909 + i1196 + i1271];
						massmatrix[i1164 - 1][i1196 - 1] += v[2909 + i1196 + i1302];
					};/* end for */
				};/* end for */
			};/* end for */
		};/* end for */
	}
	else {
	};
	if (i237 == 2) {
		i1497 = i0;
		i0 = i0 + i1869;
		i1414 = i0;
		i0 = i0 + i1869;
		i1362 = i0;
		i0 = i0 + i1869;
		v[1310] = 0e0;
		for (i1311 = 1; i1311 <= i2; i1311++) {
			i1314 = 3 * (-1 + i1311);
			i1312 = 1 + i1314;
			i1315 = 2 + i1314;
			i1317 = 3 + i1314;
			v[1310] = Power(tangent[i1312 - 1][i1312 - 1], 2) + Power(tangent[i1315 - 1][i1315 - 1], 2) + Power(tangent[i1317
				- 1][i1317 - 1], 2) + v[1310];
		};/* end for */
		v[1359] = -2e0*v[233] * sqrt(v[1310] / (3e0*i2));
		for (i1320 = 1; i1320 <= i2; i1320++) {
			i1322 = 3 * (-1 + i1320);
			i1329 = 3 + i1322;
			i1327 = 2 + i1322;
			i1325 = 1 + i1322;
			v[1321] = Xref[i1325 - 1];
			v[1323] = Xref[i1327 - 1];
			v[1324] = Xref[i1329 - 1];
			v[1358] = v[1359] * (-v[1324] - v[1321] * v[149] - v[1323] * v[150] - v[1324] * v[151] - v[183] + v[2909 + i1329 + i7]
				+ Xcur[i1329 - 1]);
			v[1360] = v[1359] * (-v[1323] - v[1321] * v[146] - v[1323] * v[147] - v[1324] * v[148] - v[181] + v[2909 + i1327 + i7]
				+ Xcur[i1327 - 1]);
			v[1361] = v[1359] * (-v[1321] - v[1321] * v[143] - v[1323] * v[144] - v[1324] * v[145] - v[179] + v[2909 + i1325 + i7]
				+ Xcur[i1325 - 1]);
			for (i1364 = 1; i1364 <= i1869; i1364++) {
				v[2909 + i1362 + i1364] = 0e0;
			};/* end for */
			v[2909 + i1329 + i1362] = -v[1358] + v[2909 + i1329 + i1362];
			v[2909 + i1327 + i1362] = -v[1360] + v[2909 + i1327 + i1362];
			v[2909 + i1325 + i1362] = -v[1361] + v[2909 + i1325 + i1362];
			v[1365] = v[1358] / i2;
			v[1366] = (v[1324] * v[1358] - v[12] * v[1365]) / v[25];
			v[1367] = (v[1323] * v[1358] - v[11] * v[1365]) / v[25];
			v[1368] = (v[1321] * v[1358] - v[10] * v[1365]) / v[25];
			v[1369] = v[1360] / i2;
			v[1370] = (v[1324] * v[1360] - v[12] * v[1369]) / v[25];
			v[1371] = (v[1323] * v[1360] - v[11] * v[1369]) / v[25];
			v[1372] = (v[1321] * v[1360] - v[10] * v[1369]) / v[25];
			v[1373] = v[1361] / i2;
			v[1374] = (v[1324] * v[1361] - v[12] * v[1373]) / v[25];
			v[1375] = (v[1323] * v[1361] - v[11] * v[1373]) / v[25];
			v[1376] = (v[1321] * v[1361] - v[10] * v[1373]) / v[25];
			i1437 = i0;
			i0 = i0 + i4;
			i1449 = i0;
			i0 = i0 + i4;
			i1461 = i0;
			i0 = i0 + i4;
			for (i1377 = i4; i1377 >= 1; i1377 = -1 + i1377) {
				v[1401] = v[2909 + i1377 + i279] / 6e0;
				v[2909 + i1377 + i1437] = v[1401];
				v[2909 + i1377 + i1449] = v[1401];
				v[2909 + i1377 + i1461] = v[1401];
				v[1385] = v[2909 + i1377 + i267];
				v[1383] = v[2909 + i1377 + i264];
				v[1381] = v[2909 + i1377 + i261];
				v[1378] = v[1366] * v[1381];
				v[1379] = v[1378] + v[1367] * v[1383];
				v[1380] = v[1379] + v[1368] * v[1385];
				v[1949] = v[1380] * v[1401];
				v[1382] = v[1370] * v[1381];
				v[1384] = v[1382] + v[1371] * v[1383];
				v[1386] = v[1384] + v[1372] * v[1385];
				v[1950] = v[1386] * v[1401];
				v[1387] = v[1374] * v[1381];
				v[1388] = v[1375] * v[1383] + v[1387];
				v[1389] = v[1376] * v[1385] + v[1388];
				v[1951] = v[1389] * v[1401];
				i1390 = (int)(v[2909 + i1377 + i281]);
				v[2909 + i1362 + i1390] = v[1949] + v[2909 + i1362 + i1390];
				i1391 = (int)(v[2909 + i1377 + i284]);
				v[2909 + i1362 + i1391] = v[1950] + v[2909 + i1362 + i1391];
				i1393 = (int)(v[2909 + i1377 + i288]);
				v[2909 + i1362 + i1393] = v[1951] + v[2909 + i1362 + i1393];
				i1394 = (int)(v[2909 + i1377 + i291]);
				v[2909 + i1362 + i1394] = v[1949] + v[2909 + i1362 + i1394];
				i1396 = (int)(v[2909 + i1377 + i295]);
				v[2909 + i1362 + i1396] = v[1950] + v[2909 + i1362 + i1396];
				i1398 = (int)(v[2909 + i1377 + i299]);
				v[2909 + i1362 + i1398] = v[1951] + v[2909 + i1362 + i1398];
				i1399 = (int)(v[2909 + i1377 + i302]);
				v[2909 + i1362 + i1399] = v[1949] + v[2909 + i1362 + i1399];
				i1400 = (int)(v[2909 + i1377 + i305]);
				v[2909 + i1362 + i1400] = v[1950] + v[2909 + i1362 + i1400];
				i1402 = (int)(v[2909 + i1377 + i309]);
				v[2909 + i1362 + i1402] = v[1951] + v[2909 + i1362 + i1402];
			};/* end for */
			for (i1403 = i2; i1403 >= 1; i1403 = -1 + i1403) {
				v[1404] = v[1365];
				v[1405] = v[1369];
				v[1406] = v[1373];
				i1407 = (int)(v[2909 + i1403 + i724]);
				v[2909 + i1362 + i1407] = v[1404] + v[2909 + i1362 + i1407];
				i1408 = (int)(v[2909 + i1403 + i727]);
				v[2909 + i1362 + i1408] = v[1405] + v[2909 + i1362 + i1408];
				i1409 = (int)(v[2909 + i1403 + i730]);
				v[2909 + i1362 + i1409] = v[1406] + v[2909 + i1362 + i1409];
			};/* end for */
			v[1366] = 0e0;
			v[1367] = 0e0;
			v[1368] = 0e0;
			v[1370] = 0e0;
			v[1371] = 0e0;
			v[1372] = 0e0;
			v[1374] = 0e0;
			v[1375] = 0e0;
			v[1376] = 0e0;
			v[1365] = 0e0;
			v[1369] = 0e0;
			v[1373] = 0e0;
			for (i1356 = 1; i1356 <= i1869; i1356++) {
				for (i1416 = 1; i1416 <= i1869; i1416++) {
					v[2909 + i1414 + i1416] = 0e0;
				};/* end for */
				v[2909 + i1356 + i1414] = 1e0 + v[2909 + i1356 + i1414];
				v[1417] = 0e0;
				v[1418] = 0e0;
				v[1419] = 0e0;
				v[1420] = 0e0;
				v[1421] = 0e0;
				v[1422] = 0e0;
				v[1423] = 0e0;
				v[1424] = 0e0;
				v[1425] = 0e0;
				v[1426] = 0e0;
				v[1427] = 0e0;
				v[1428] = 0e0;
				for (i1429 = 1; i1429 <= i2; i1429++) {
					i1435 = (int)(v[2909 + i1429 + i724]);
					i1433 = (int)(v[2909 + i1429 + i727]);
					i1431 = (int)(v[2909 + i1429 + i730]);
					v[1430] = v[2909 + i1414 + i1431];
					v[2909 + i1414 + i1431] = v[1430];
					v[1432] = v[2909 + i1414 + i1433];
					v[2909 + i1414 + i1433] = v[1432];
					v[1434] = v[2909 + i1414 + i1435];
					v[2909 + i1414 + i1435] = v[1434];
					v[1417] = v[1417] + v[1430];
					v[1418] = v[1418] + v[1432];
					v[1419] = v[1419] + v[1434];
				};/* end for */
				for (i1436 = 1; i1436 <= i4; i1436++) {
					v[1475] = v[2909 + i1436 + i261];
					v[1474] = v[2909 + i1436 + i264];
					v[1473] = v[2909 + i1436 + i267];
					i1472 = (int)(v[2909 + i1436 + i281]);
					i1469 = (int)(v[2909 + i1436 + i284]);
					v[1466] = v[2909 + i1436 + i1461];
					i1465 = (int)(v[2909 + i1436 + i288]);
					i1460 = (int)(v[2909 + i1436 + i291]);
					i1457 = (int)(v[2909 + i1436 + i295]);
					v[1454] = v[2909 + i1436 + i1449];
					i1453 = (int)(v[2909 + i1436 + i299]);
					i1448 = (int)(v[2909 + i1436 + i302]);
					i1445 = (int)(v[2909 + i1436 + i305]);
					v[1442] = v[2909 + i1436 + i1437];
					i1441 = (int)(v[2909 + i1436 + i309]);
					v[1440] = v[2909 + i1414 + i1441];
					v[2909 + i1414 + i1441] = v[1440];
					v[1444] = v[2909 + i1414 + i1445];
					v[2909 + i1414 + i1445] = v[1444];
					v[1447] = v[2909 + i1414 + i1448];
					v[2909 + i1414 + i1448] = v[1447];
					v[1452] = v[2909 + i1414 + i1453];
					v[2909 + i1414 + i1453] = v[1452];
					v[1456] = v[2909 + i1414 + i1457];
					v[2909 + i1414 + i1457] = v[1456];
					v[1459] = v[2909 + i1414 + i1460];
					v[2909 + i1414 + i1460] = v[1459];
					v[1464] = v[2909 + i1414 + i1465];
					v[1463] = v[1440] * v[1442] + v[1452] * v[1454] + v[1464] * v[1466];
					v[2909 + i1414 + i1465] = v[1464];
					v[1468] = v[2909 + i1414 + i1469];
					v[1467] = v[1442] * v[1444] + v[1454] * v[1456] + v[1466] * v[1468];
					v[2909 + i1414 + i1469] = v[1468];
					v[1471] = v[2909 + i1414 + i1472];
					v[1470] = v[1442] * v[1447] + v[1454] * v[1459] + v[1466] * v[1471];
					v[2909 + i1414 + i1472] = v[1471];
					v[1420] = v[1420] + v[1463] * v[1473];
					v[1421] = v[1421] + v[1463] * v[1474];
					v[1422] = v[1422] + v[1463] * v[1475];
					v[1423] = v[1423] + v[1467] * v[1473];
					v[1424] = v[1424] + v[1467] * v[1474];
					v[1425] = v[1425] + v[1467] * v[1475];
					v[1426] = v[1426] + v[1470] * v[1473];
					v[1427] = v[1427] + v[1470] * v[1474];
					v[1428] = v[1428] + v[1470] * v[1475];
				};/* end for */
				v[1489] = v[2909 + i1325 + i1414];
				v[1960] = v[1428] / v[25];
				v[1959] = v[1427] / v[25];
				v[1958] = v[1426] / v[25];
				v[1957] = v[1425] / v[25];
				v[1956] = v[1424] / v[25];
				v[1955] = v[1423] / v[25];
				v[1954] = v[1422] / v[25];
				v[1953] = v[1421] / v[25];
				v[1952] = v[1420] / v[25];
				v[1417] = v[1417] - v[10] * v[1952];
				v[1420] = 0e0;
				v[1417] = v[1417] - v[11] * v[1953];
				v[1421] = 0e0;
				v[1417] = v[1417] - v[12] * v[1954];
				v[1422] = 0e0;
				v[1479] = v[1417] / i2 + v[1321] * v[1952] + v[1323] * v[1953] + v[1324] * v[1954];
				v[1417] = 0e0;
				v[1418] = v[1418] - v[10] * v[1955];
				v[1423] = 0e0;
				v[1418] = v[1418] - v[11] * v[1956];
				v[1424] = 0e0;
				v[1418] = v[1418] - v[12] * v[1957];
				v[1425] = 0e0;
				v[1483] = v[1418] / i2 + v[1321] * v[1955] + v[1323] * v[1956] + v[1324] * v[1957];
				v[1418] = 0e0;
				v[1419] = v[1419] - v[10] * v[1958];
				v[1426] = 0e0;
				v[1419] = v[1419] - v[11] * v[1959];
				v[1427] = 0e0;
				v[1419] = v[1419] - v[12] * v[1960];
				v[1428] = 0e0;
				v[1487] = v[1419] / i2 + v[1321] * v[1958] + v[1323] * v[1959] + v[1324] * v[1960];
				v[1419] = 0e0;
				v[2909 + i1325 + i1414] = v[1489];
				v[1491] = v[2909 + i1327 + i1414];
				v[2909 + i1327 + i1414] = v[1491];
				v[1493] = v[2909 + i1329 + i1414];
				v[2909 + i1329 + i1414] = v[1493];
				v[1494] = -(v[1359] * (v[1479] - v[1489]));
				v[1495] = -(v[1359] * (v[1483] - v[1491]));
				v[1496] = -(v[1359] * (v[1487] - v[1493]));
				for (i1499 = 1; i1499 <= i1869; i1499++) {
					v[2909 + i1497 + i1499] = 0e0;
				};/* end for */
				v[2909 + i1329 + i1497] = -v[1496] + v[2909 + i1329 + i1497];
				v[2909 + i1327 + i1497] = -v[1495] + v[2909 + i1327 + i1497];
				v[2909 + i1325 + i1497] = -v[1494] + v[2909 + i1325 + i1497];
				v[1500] = v[1496] / i2;
				v[1501] = (v[1324] * v[1496] - v[12] * v[1500]) / v[25];
				v[1502] = (v[1323] * v[1496] - v[11] * v[1500]) / v[25];
				v[1503] = (v[1321] * v[1496] - v[10] * v[1500]) / v[25];
				v[1504] = v[1495] / i2;
				v[1505] = (v[1324] * v[1495] - v[12] * v[1504]) / v[25];
				v[1506] = (v[1323] * v[1495] - v[11] * v[1504]) / v[25];
				v[1507] = (v[1321] * v[1495] - v[10] * v[1504]) / v[25];
				v[1508] = v[1494] / i2;
				v[1509] = (v[1324] * v[1494] - v[12] * v[1508]) / v[25];
				v[1510] = (v[1323] * v[1494] - v[11] * v[1508]) / v[25];
				v[1511] = (v[1321] * v[1494] - v[10] * v[1508]) / v[25];
				for (i1512 = i4; i1512 >= 1; i1512 = -1 + i1512) {
					v[1536] = v[2909 + i1512 + i279] / 6e0;
					v[1520] = v[2909 + i1512 + i267];
					v[1518] = v[2909 + i1512 + i264];
					v[1516] = v[2909 + i1512 + i261];
					v[1513] = v[1501] * v[1516];
					v[1514] = v[1513] + v[1502] * v[1518];
					v[1515] = v[1514] + v[1503] * v[1520];
					v[1961] = v[1515] * v[1536];
					v[1517] = v[1505] * v[1516];
					v[1519] = v[1517] + v[1506] * v[1518];
					v[1521] = v[1519] + v[1507] * v[1520];
					v[1962] = v[1521] * v[1536];
					v[1522] = v[1509] * v[1516];
					v[1523] = v[1510] * v[1518] + v[1522];
					v[1524] = v[1511] * v[1520] + v[1523];
					v[1963] = v[1524] * v[1536];
					i1525 = (int)(v[2909 + i1512 + i281]);
					v[2909 + i1497 + i1525] = v[1961] + v[2909 + i1497 + i1525];
					i1526 = (int)(v[2909 + i1512 + i284]);
					v[2909 + i1497 + i1526] = v[1962] + v[2909 + i1497 + i1526];
					i1528 = (int)(v[2909 + i1512 + i288]);
					v[2909 + i1497 + i1528] = v[1963] + v[2909 + i1497 + i1528];
					i1529 = (int)(v[2909 + i1512 + i291]);
					v[2909 + i1497 + i1529] = v[1961] + v[2909 + i1497 + i1529];
					i1531 = (int)(v[2909 + i1512 + i295]);
					v[2909 + i1497 + i1531] = v[1962] + v[2909 + i1497 + i1531];
					i1533 = (int)(v[2909 + i1512 + i299]);
					v[2909 + i1497 + i1533] = v[1963] + v[2909 + i1497 + i1533];
					i1534 = (int)(v[2909 + i1512 + i302]);
					v[2909 + i1497 + i1534] = v[1961] + v[2909 + i1497 + i1534];
					i1535 = (int)(v[2909 + i1512 + i305]);
					v[2909 + i1497 + i1535] = v[1962] + v[2909 + i1497 + i1535];
					i1537 = (int)(v[2909 + i1512 + i309]);
					v[2909 + i1497 + i1537] = v[1963] + v[2909 + i1497 + i1537];
				};/* end for */
				for (i1538 = i2; i1538 >= 1; i1538 = -1 + i1538) {
					v[1539] = v[1500];
					v[1540] = v[1504];
					v[1541] = v[1508];
					i1542 = (int)(v[2909 + i1538 + i724]);
					v[2909 + i1497 + i1542] = v[1539] + v[2909 + i1497 + i1542];
					i1543 = (int)(v[2909 + i1538 + i727]);
					v[2909 + i1497 + i1543] = v[1540] + v[2909 + i1497 + i1543];
					i1544 = (int)(v[2909 + i1538 + i730]);
					v[2909 + i1497 + i1544] = v[1541] + v[2909 + i1497 + i1544];
				};/* end for */
				v[1501] = 0e0;
				v[1502] = 0e0;
				v[1503] = 0e0;
				v[1505] = 0e0;
				v[1506] = 0e0;
				v[1507] = 0e0;
				v[1509] = 0e0;
				v[1510] = 0e0;
				v[1511] = 0e0;
				v[1500] = 0e0;
				v[1504] = 0e0;
				v[1508] = 0e0;
				residual[i1356 - 1] += v[2909 + i1356 + i1362];
				for (i1412 = i1356; i1412 <= i1869; i1412++) {
					tangent[i1356 - 1][i1412 - 1] += v[2909 + i1412 + i1497];
				};/* end for */
			};/* end for */
		};/* end for */
	}
	else {
	};
#pragma endregion


	//Filling the other triangular region of the symmetric matrix
	for (int i = 0; i < 3 * n_vertices; i++)
	{
		for (int j = i + 1; j < 3 * n_vertices; j++)
		{
			tangent[j][i] = tangent[i][j];
		}
	}
	for (int i = 0; i < 3 * n_vertices; i++)
	{
		for (int j = i + 1; j < 3 * n_vertices; j++)
		{
			massmatrix[j][i] = massmatrix[i][j];
		}
	}
	
	//db.PrintPtr(elementData,20);
	//db.PrintPtr(domainData, 4);
	//db.PrintPtr(Xref, (int)ptr_cad->vertices.size() * 3);
	//db.PrintPtr(Xcur, (int)ptr_cad->vertices.size() * 3);
	//db.PrintPtr(displacement, (int)ptr_cad->vertices.size() * 3);
	//db.PrintPtr((double*)triangles, 3 * (int)ptr_cad->vertices.size());
	//db.PrintPtr(residual, 3 * n_vertices);
	//db.PrintPtr(local_stiffness, 3 * n_vertices, 3 * n_vertices);

	//db.PrintPtr(residualload, 3 * n_vertices);
	//db.PrintPtr(massmatrix, 3 * n_vertices, 3 * n_vertices);

	
	DampingContributions();//damping is called first, because it employs stiffness and mass just evaluated
	InertialContributions();

	////Updating local_abs_load:
	//for (int i = 0; i < 3 * n_vertices; i++)
	//{
	//	count_local_abs_load[i]++;
	//	local_abs_load[i] += abs(local_residual[i]);
	//}

	MountFieldLoading();

	delete[] v;
}

//Explicit
void VEMPolyhedron::EvaluateExplicit()
{

}


void VEMPolyhedron::EvaluateAccelerations()
{

}

void VEMPolyhedron::MountFieldLoading()
{
	if (lumped_mass_flag)
	{
		if (l_factor != 0.0)
		{
			//CAD pointer
			STLSurface* ptr_cad = static_cast<STLSurface*>(db.cad_data[CADDATA_ID - 1]);
			for (int i = 0; i < n_vertices; i++)
			{
				local_residual[i * 3 + 0] -= mass * ptr_cad->vertice_factors[i] * l_factor * db.environment->G(0, 0);
				local_residual[i * 3 + 1] -= mass * ptr_cad->vertice_factors[i] * l_factor * db.environment->G(1, 0);
				local_residual[i * 3 + 2] -= mass * ptr_cad->vertice_factors[i] * l_factor * db.environment->G(2, 0);
				//Updating local_abs_load:
				/*local_abs_load[i * 3 + 0] += abs(mass * ptr_cad->vertice_factors[i] * l_factor * db.environment->G(0, 0));
				local_abs_load[i * 3 + 1] += abs(mass * ptr_cad->vertice_factors[i] * l_factor * db.environment->G(1, 0));
				local_abs_load[i * 3 + 2] += abs(mass * ptr_cad->vertice_factors[i] * l_factor * db.environment->G(2, 0));
				count_local_abs_load[i * 3 + 0]++;
				count_local_abs_load[i * 3 + 1]++;
				count_local_abs_load[i * 3 + 2]++;*/
			}
		}	
	}
	else
	{
		//local_e_load already with negative sign, from AceGen code
		for (int i = 0; i < 3 * n_vertices; i++)
			local_residual[i] += local_e_load[i];
		////Updating local_abs_load:
		//for (int i = 0; i < 3 * n_vertices; i++)
		//{
		//	count_local_abs_load[i]++;
		//	local_abs_load[i] += abs(local_residual[i]);
		//}
	}
	
}

//Calcula contribuições de inércia
void VEMPolyhedron::InertialContributions()
{
	//CAD pointer
	STLSurface* ptr_cad = static_cast<STLSurface*>(db.cad_data[CADDATA_ID - 1]);
	//Contributions - Newmark
	double a1;
	if (typeid(*db.solution[db.current_solution_number - 1]) == typeid(Dynamic))
	{
		Dynamic* ptr_sol = static_cast<Dynamic*>(db.solution[db.current_solution_number - 1]);
		a1 = ptr_sol->a1;
	}
	else
		a1 = 0.0;

	//Lumped mass matrix - direct insertion on residual and tangent
	if (lumped_mass_flag)
	{
		for (int i = 0; i < n_vertices; i++)
		{
			local_residual[i * 3 + 0] += mass * ptr_cad->vertice_factors[i] * (db.super_nodes[super_node - 1]->accel[i * 3 + 0]);
			local_residual[i * 3 + 1] += mass * ptr_cad->vertice_factors[i] * (db.super_nodes[super_node - 1]->accel[i * 3 + 1]);
			local_residual[i * 3 + 2] += mass * ptr_cad->vertice_factors[i] * (db.super_nodes[super_node - 1]->accel[i * 3 + 2]);
			local_stiffness[i * 3 + 0][i * 3 + 0] += a1 * mass * ptr_cad->vertice_factors[i];
			local_stiffness[i * 3 + 1][i * 3 + 1] += a1 * mass * ptr_cad->vertice_factors[i];
			local_stiffness[i * 3 + 2][i * 3 + 2] += a1 * mass * ptr_cad->vertice_factors[i];
		}
	}
	else
	{
		//Contribution to the residual and tangent (already evaluated in Mount() function)
		//Residual - inertial contribution
		Matrix temp(3 * n_vertices);
		for (int i = 0; i < 3 * n_vertices; i++)
			for (int j = 0; j < 3 * n_vertices; j++)
				temp(i, 0) = temp(i, 0) + local_mass[i][j] * db.super_nodes[super_node - 1]->accel[j];
		for (int i = 0; i < 3 * n_vertices; i++)
		{
			local_residual[i] += temp(i,0);
			for (int j = 0; j < 3 * n_vertices; j++)
				local_stiffness[i][j] += a1 * local_mass[i][j];
		}
		////Updating local_abs_load:
		//for (int i = 0; i < 3 * n_vertices; i++)
		//{
		//	count_local_abs_load[i]++;
		//	local_abs_load[i] += abs(temp(i, 0));
		//}
	}
}

void VEMPolyhedron::DampingContributions()
{
	//Contributions - Newmark
	
	if (typeid(*db.solution[db.current_solution_number - 1]) == typeid(Dynamic))
	{
		double a4;
		Dynamic* ptr_sol = static_cast<Dynamic*>(db.solution[db.current_solution_number - 1]);
		a4 = ptr_sol->a4;

		if (damping_flag)
		{
			if (first_damping_evaluation)
			{
				for (int i = 0; i < 3 * n_vertices; i++)
				{
					for (int j = 0; j < 3 * n_vertices; j++)
						local_damping[i][j] = alpha * local_mass[i][j] + beta * local_stiffness[i][j];
				}
				first_damping_evaluation = false;
			}
			//Contribution to the residual and tangent
			Matrix temp(3 * n_vertices);
			for (int i = 0; i < 3 * n_vertices; i++)
				for (int j = 0; j < 3 * n_vertices; j++)
					temp(i, 0) = temp(i, 0) + local_damping[i][j] * db.super_nodes[super_node - 1]->vel[j];
			for (int i = 0; i < 3 * n_vertices; i++)
			{
				local_residual[i] += temp(i, 0);
				for (int j = 0; j < 3 * n_vertices; j++)
					local_stiffness[i][j] += a4 * local_damping[i][j];
			}
			////Updating local_abs_load:
			//for (int i = 0; i < 3 * n_vertices; i++)
			//{
			//	count_local_abs_load[i]++;
			//	local_abs_load[i] += abs(temp(i, 0));
			//}
		}
	}
}

void VEMPolyhedron::PostProcessing()
{
	/*
	{{1, "Sxx" -> AceGen`$V[1660, 1]}, {2, 
  "Sxy" -> AceGen`$V[1661, 1]}, {3, "Sxz" -> AceGen`$V[1662, 1]}, {4, 
  "Syx" -> AceGen`$V[1663, 1]}, {5, "Syy" -> AceGen`$V[1664, 1]}, {6, 
  "Syz" -> AceGen`$V[1665, 1]}, {7, "Szx" -> AceGen`$V[1666, 1]}, {8, 
  "Szy" -> AceGen`$V[1667, 1]}, {9, "Szz" -> AceGen`$V[1668, 1]}, {10,
   "p" -> AceGen`$V[1669, 1]}, {11, "Exx" -> AceGen`$V[1642, 1]}, {12,
   "Exy" -> AceGen`$V[1643, 1]}, {13, 
  "Exz" -> AceGen`$V[1644, 1]}, {14, 
  "Eyx" -> AceGen`$V[1643, 1]}, {15, 
  "Eyy" -> AceGen`$V[1645, 1]}, {16, 
  "Eyz" -> AceGen`$V[1646, 1]}, {17, 
  "Ezx" -> AceGen`$V[1644, 1]}, {18, 
  "Ezy" -> AceGen`$V[1646, 1]}, {19, 
  "Ezz" -> AceGen`$V[1647, 1]}, {20, 
  "Mises stress" -> AceGen`$V[1673, 1]}, {21, 
  "Specific strain energy" -> AceGen`$V[1641, 1]}, {22, 
  "xxC" -> AceGen`$V[1568, 1]}, {23, 
  "xyC" -> AceGen`$V[1569, 1]}, {24, 
  "xzC" -> AceGen`$V[1570, 1]}, {25, "V" -> AceGen`$V[1571, 1]}, {26, 
  "A" -> AceGen`$V[1572, 1]}, {27, 
  "WSEproj" -> AceGen`$V[1674, 1]}, {28, "WSESt" -> 0}, {29, 
  "WSETot" -> AceGen`$V[1674, 1]}}
	*/
	double v[36000];
	int i1546, i1550, i1551, i1553, i1709, i1711, i1713, i1722, i1723, i1724, i1725, i1727
		, i1728, i1729, i1730, i1731, i1732, i1733, i1734, i1736, i1737, i1738, i1739, i2042
		, i0;
	i0 = 32;
	i1546 = (int)(elementData[0]);
	i2042 = 3 * i1546;
	i1550 = (int)(elementData[2]);
	i1551 = i0;
	i0 = i0 + i2042;
	for (i1553 = 1; i1553 <= i2042; i1553++) {
		v[2958 + i1551 + i1553] = displacement[i1553 - 1];
	};/* end for */
	v[1558] = elementData[7];
	v[1559] = elementData[8];
	v[1560] = elementData[9];
	v[1562] = elementData[11];
	v[1563] = elementData[12];
	v[1564] = elementData[13];
	v[1569] = elementData[18];
	v[1586] = 1e0 + elementData[6];
	v[1587] = 1e0 + elementData[10];
	v[2048] = v[1560] * v[1564] - v[1563] * v[1587];
	v[1588] = 1e0 + elementData[14];
	v[2047] = v[1562] * v[1563] - v[1560] * v[1588];
	v[2045] = -(v[1562] * v[1564]) + v[1587] * v[1588];
	v[2043] = v[1586] * v[2045] + v[1558] * v[2047] + v[1559] * v[2048];
	v[1624] = (v[1558] * v[1558]) + (v[1559] * v[1559]) + (v[1586] * v[1586]);
	v[1627] = (v[1560] * v[1560]) + (v[1562] * v[1562]) + (v[1587] * v[1587]);
	v[1629] = (v[1563] * v[1563]) + (v[1564] * v[1564]) + (v[1588] * v[1588]);
	v[1646] = v[1624] + v[1627] + v[1629];
	v[2044] = -1e0 + (v[2043] * v[2043]);
	v[1648] = 1e0 / Power(v[2043], 0.6666666666666666e0);
	v[1631] = domainData[0];
	v[1632] = domainData[1];
	v[1637] = v[1631] / (2e0 + 2e0*v[1632]);
	v[2046] = v[1637] * v[1648];
	v[1638] = v[1631] / (3e0 - 6e0*v[1632]);
	v[2049] = (v[1638] * v[2044]) / (2e0*v[2043]);
	v[1639] = (v[1637] * (-3e0 + v[1646] * v[1648])) / 2e0 + (v[1638] * (v[2044] - 2e0*log(v[2043]))) / 4e0;
	v[1641] = (v[1559] * v[1562] + v[1560] * v[1586] + v[1558] * v[1587]) / 2e0;
	v[1642] = (v[1558] * v[1564] + v[1563] * v[1586] + v[1559] * v[1588]) / 2e0;
	v[1644] = (v[1560] * v[1563] + v[1564] * v[1587] + v[1562] * v[1588]) / 2e0;
	v[1647] = -(v[1637] * v[1646]) / (3e0*Power(v[2043], 0.16666666666666669e1)) + v[2049];
	v[1649] = v[1647] * v[2045] + v[1586] * v[2046];
	v[1650] = (v[1559] * v[1564] - v[1558] * v[1588])*v[1647] + v[1560] * v[2046];
	v[1651] = (v[1558] * v[1562] - v[1559] * v[1587])*v[1647] + v[1563] * v[2046];
	v[1652] = v[1558] * v[2046] + v[1647] * v[2047];
	v[1653] = (-(v[1559] * v[1563]) + v[1586] * v[1588])*v[1647] + v[1587] * v[2046];
	v[1654] = (v[1559] * v[1560] - v[1562] * v[1586])*v[1647] + v[1564] * v[2046];
	v[1658] = (v[1586] * v[1649] + v[1560] * v[1650] + v[1563] * v[1651]) / v[2043];
	v[1659] = (v[1558] * v[1649] + v[1587] * v[1650] + v[1564] * v[1651]) / v[2043];
	v[1660] = (v[1559] * v[1649] + v[1562] * v[1650] + v[1588] * v[1651]) / v[2043];
	v[1662] = (v[1558] * v[1652] + v[1587] * v[1653] + v[1564] * v[1654]) / v[2043];
	v[1663] = (v[1559] * v[1652] + v[1562] * v[1653] + v[1588] * v[1654]) / v[2043];
	v[1666] = (v[1562] * ((v[1558] * v[1563] - v[1564] * v[1586])*v[1647] + v[1562] * v[2046]) + v[1588] * ((-
		(v[1558] * v[1560]) + v[1586] * v[1587])*v[1647] + v[1588] * v[2046]) + v[1559] * (v[1559] * v[2046]
			+ v[1647] * v[2048])) / v[2043];
	v[1672] = v[1569] * v[1639];
	elementPostprocessing[0] = v[1658];
	elementPostprocessing[1] = v[1659];
	elementPostprocessing[2] = v[1660];
	elementPostprocessing[3] = v[1659];
	elementPostprocessing[4] = v[1662];
	elementPostprocessing[5] = v[1663];
	elementPostprocessing[6] = v[1660];
	elementPostprocessing[7] = v[1663];
	elementPostprocessing[8] = v[1666];
	elementPostprocessing[9] = v[2049];
	elementPostprocessing[10] = (-1e0 + v[1624]) / 2e0;
	elementPostprocessing[11] = v[1641];
	elementPostprocessing[12] = v[1642];
	elementPostprocessing[13] = v[1641];
	elementPostprocessing[14] = (-1e0 + v[1627]) / 2e0;
	elementPostprocessing[15] = v[1644];
	elementPostprocessing[16] = v[1642];
	elementPostprocessing[17] = v[1644];
	elementPostprocessing[18] = (-1e0 + v[1629]) / 2e0;
	elementPostprocessing[19] = sqrt(0.15e1*(2e0*(v[1659] * v[1659]) + 2e0*(v[1660] * v[1660]) + 2e0*
		(v[1663] * v[1663]) + Power(v[1658] - v[2049], 2) + Power(v[1662] - v[2049], 2) + Power(v[1666] - v[2049], 2)));
	elementPostprocessing[20] = v[1639];
	elementPostprocessing[21] = elementData[15];
	elementPostprocessing[22] = elementData[16];
	elementPostprocessing[23] = elementData[17];
	elementPostprocessing[24] = v[1569];
	elementPostprocessing[25] = elementData[19];
	elementPostprocessing[26] = v[1672];
	elementPostprocessing[27] = 0e0;
	elementPostprocessing[28] = v[1672];
	v[1703] = domainData[6];
	if ((int)(domainData[9]) == 1) {
		v[1708] = 0e0;
		for (i1709 = 1; i1709 <= i1550; i1709++) {
			i1713 = 4 * (-1 + i1709);
			i1733 = 3 * (-1 + tetrahedrons[3 + i1713]);
			i1739 = 3 + i1733;
			v[1754] = Xref[i1739 - 1];
			v[1796] = -v[1754] + v[2958 + i1551 + i1739] + Xcur[i1739 - 1];
			i1734 = 2 + i1733;
			v[1749] = Xref[i1734 - 1];
			v[1789] = -v[1749] + v[2958 + i1551 + i1734] + Xcur[i1734 - 1];
			i1725 = 1 + i1733;
			v[1744] = Xref[i1725 - 1];
			v[1782] = -v[1744] + v[2958 + i1551 + i1725] + Xcur[i1725 - 1];
			i1731 = 3 * (-1 + tetrahedrons[2 + i1713]);
			i1738 = 3 + i1731;
			v[1753] = Xref[i1738 - 1];
			v[1795] = -v[1753] + v[2958 + i1551 + i1738] + Xcur[i1738 - 1];
			v[1764] = v[1753] - v[1754];
			i1732 = 2 + i1731;
			v[1748] = Xref[i1732 - 1];
			v[1788] = -v[1748] + v[2958 + i1551 + i1732] + Xcur[i1732 - 1];
			v[1761] = v[1748] - v[1749];
			i1724 = 1 + i1731;
			v[1743] = Xref[i1724 - 1];
			v[1781] = -v[1743] + v[2958 + i1551 + i1724] + Xcur[i1724 - 1];
			v[1758] = v[1743] - v[1744];
			i1729 = 3 * (-1 + tetrahedrons[1 + i1713]);
			i1737 = 3 + i1729;
			v[1752] = Xref[i1737 - 1];
			v[1794] = -v[1752] + v[2958 + i1551 + i1737] + Xcur[i1737 - 1];
			v[1763] = v[1752] - v[1754];
			i1730 = 2 + i1729;
			v[1747] = Xref[i1730 - 1];
			v[1787] = -v[1747] + v[2958 + i1551 + i1730] + Xcur[i1730 - 1];
			v[1760] = v[1747] - v[1749];
			v[2052] = -(v[1761] * v[1763]) + v[1760] * v[1764];
			i1723 = 1 + i1729;
			v[1742] = Xref[i1723 - 1];
			v[1780] = -v[1742] + v[2958 + i1551 + i1723] + Xcur[i1723 - 1];
			v[1757] = v[1742] - v[1744];
			i1727 = 3 * (-1 + tetrahedrons[i1713]);
			i1736 = 3 + i1727;
			v[1751] = Xref[i1736 - 1];
			v[1793] = -v[1751] + v[2958 + i1551 + i1736] + Xcur[i1736 - 1];
			v[1762] = v[1751] - v[1754];
			i1728 = 2 + i1727;
			v[1746] = Xref[i1728 - 1];
			v[1786] = -v[1746] + v[2958 + i1551 + i1728] + Xcur[i1728 - 1];
			v[1759] = v[1746] - v[1749];
			v[2051] = v[1761] * v[1762] - v[1759] * v[1764];
			v[2050] = -(v[1760] * v[1762]) + v[1759] * v[1763];
			i1722 = 1 + i1727;
			v[1741] = Xref[i1722 - 1];
			v[1779] = -v[1741] + v[2958 + i1551 + i1722] + Xcur[i1722 - 1];
			v[1756] = v[1741] - v[1744];
			v[1775] = v[1758] * v[2050] + v[1757] * v[2051] + v[1756] * v[2052];
			v[1774] = (-(v[1757] * v[1759]) + v[1756] * v[1760]) / v[1775];
			v[1773] = (v[1757] * v[1762] - v[1756] * v[1763]) / v[1775];
			v[1772] = v[2050] / v[1775];
			v[1771] = (v[1758] * v[1759] - v[1756] * v[1761]) / v[1775];
			v[1770] = (-(v[1758] * v[1762]) + v[1756] * v[1764]) / v[1775];
			v[1769] = v[2051] / v[1775];
			v[1768] = (-(v[1758] * v[1760]) + v[1757] * v[1761]) / v[1775];
			v[1778] = -v[1768] - v[1771] - v[1774];
			v[1822] = 1e0 + v[1768] * v[1793] + v[1771] * v[1794] + v[1774] * v[1795] + v[1778] * v[1796];
			v[1819] = v[1768] * v[1786] + v[1771] * v[1787] + v[1774] * v[1788] + v[1778] * v[1789];
			v[1816] = v[1768] * v[1779] + v[1771] * v[1780] + v[1774] * v[1781] + v[1778] * v[1782];
			v[1767] = (v[1758] * v[1763] - v[1757] * v[1764]) / v[1775];
			v[1777] = -v[1767] - v[1770] - v[1773];
			v[1821] = v[1767] * v[1793] + v[1770] * v[1794] + v[1773] * v[1795] + v[1777] * v[1796];
			v[1818] = 1e0 + v[1767] * v[1786] + v[1770] * v[1787] + v[1773] * v[1788] + v[1777] * v[1789];
			v[1815] = v[1767] * v[1779] + v[1770] * v[1780] + v[1773] * v[1781] + v[1777] * v[1782];
			v[1766] = v[2052] / v[1775];
			v[1776] = -v[1766] - v[1769] - v[1772];
			v[1820] = v[1766] * v[1793] + v[1769] * v[1794] + v[1772] * v[1795] + v[1776] * v[1796];
			v[1817] = v[1766] * v[1786] + v[1769] * v[1787] + v[1772] * v[1788] + v[1776] * v[1789];
			v[1814] = 1e0 + v[1766] * v[1779] + v[1769] * v[1780] + v[1772] * v[1781] + v[1776] * v[1782];
			v[1829] = v[1816] * (-(v[1818] * v[1820]) + v[1817] * v[1821]) + v[1815] * (v[1819] * v[1820] - v[1817] * v[1822])
				+ v[1814] * (-(v[1819] * v[1821]) + v[1818] * v[1822]);
			v[1831] = (v[1775] * (2e0*v[1637] * (-3e0 + ((v[1814] * v[1814]) + (v[1815] * v[1815]) + (v[1816] * v[1816]) +
				(v[1817] * v[1817]) + (v[1818] * v[1818]) + (v[1819] * v[1819]) + (v[1820] * v[1820]) + (v[1821] * v[1821]) +
				(v[1822] * v[1822])) / Power(v[1829], 0.6666666666666666e0)) + v[1638] * (-1e0 + (v[1829] * v[1829]) - 2e0*log
				(v[1829])))) / 96e0;
			for (i1711 = 1; i1711 <= 4; i1711++) {
				v[1708] = v[1708] + v[1831];
			};/* end for */
		};/* end for */
		elementPostprocessing[27] = v[1708];
		elementPostprocessing[28] = v[1672] * (1e0 - v[1703]) + v[1703] * v[1708];
	}
	else {
	};
}

void VEMPolyhedron::InitialEvaluations()
{

}