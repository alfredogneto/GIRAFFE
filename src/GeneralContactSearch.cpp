#include "GeneralContactSearch.h"
#include "Database.h"
#include "CollisionDetection.h"

//Variáveis globais
extern
Database db;

//FILE* ftest;

//FILE *fdebug = fopen("debug.txt", "w");

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

GeneralContactSearch::GeneralContactSearch()
{
	method = 0;
	contactPP_list = NULL;
	contactPB_list = NULL;
	contactBOBO_list = NULL;
	contactPBO_list = NULL;
	cleanup_count = 0;
	cleanup_freq = 10;

	
	inc_len_factor = 0.1f;
	
	n_monitoring_PP = 0;
	n_active_PP = 0;
	n_monitoring_PB = 0;
	n_active_PB = 0;
	n_monitoring_BOBO = 0;
	n_active_BOBO = 0;
	n_monitoring_PBO = 0;
	n_active_PBO = 0;

	plot_solution_times = true;
	duration_verlet = 0;
	duration_linkedcells = 0;
	duration_collision_detection = 0;
	duration_mount_contact = 0;
	global_n_collisiondetection = 0;

	//Domain of interest
	max_x = 0.0f;
	min_x = 0.0f;
	max_y = 0.0f;
	min_y = 0.0f;
	max_z = 0.0f;
	min_z = 0.0f;
	user_set_boundings = false;

	cells_P_sub = NULL;
	cells_B_sub = NULL;
	cells_BO_sub = NULL;
	//ftest = fopen("test.txt", "w");
}

GeneralContactSearch::~GeneralContactSearch()
{
	if (contactPP_list != NULL)
	{
		for (int i = 0; i < db.number_particles; i++)
		{
			for (int j = 0; j < contactPP_list[i].size(); j++)
				delete contactPP_list[i][j];
			contactPP_list[i].clear();
		}
		delete[] contactPP_list;
	}

	if (contactPB_list != NULL)
	{
		for (int i = 0; i < db.number_particles; i++)
		{
			for (int j = 0; j < contactPB_list[i].size(); j++)
				delete contactPB_list[i][j];
			contactPB_list[i].clear();
		}
		delete[] contactPB_list;
	}

	if (contactBOBO_list != NULL)
	{
		for (int i = 0; i < db.number_body_geometries; i++)
		{
			for (int j = 0; j < contactBOBO_list[i].size(); j++)
				delete contactBOBO_list[i][j];
			contactBOBO_list[i].clear();
		}
		delete[] contactBOBO_list;
	}

	if (contactPBO_list != NULL)
	{
		for (int i = 0; i < db.number_particles; i++)
		{
			for (int j = 0; j < contactPBO_list[i].size(); j++)
				delete contactPBO_list[i][j];
			contactPBO_list[i].clear();
		}
		delete[] contactPBO_list;
	}

	if (cells_P_sub != NULL)
		delete[] cells_P_sub;
	if (cells_B_sub != NULL)
		delete[] cells_B_sub;
	if (cells_BO_sub != NULL)
		delete[] cells_BO_sub;
	//fclose(ftest);
}

//Pre calculation (allocation,etc.)
void GeneralContactSearch::PreCalc()
{
	//PreCalc - largest sizes of BV's and sub-BV's
	float largest_P = 0.0f;
	float largest_Psub = 0.0f;
	float largest_B = 0.0f;
	float largest_Bsub = 0.0f;
	float largest_BO = 0.0f;
	float largest_BOsub = 0.0f;

	//Setting maximum values for particles and sub particles
	for (int i = 0; i < db.number_particles; i++)
	{
		if (db.particles[i]->bv->size > largest_P)
			largest_P = db.particles[i]->bv->size;
		for (int ii = 0; ii < db.particles[i]->n_sub_bv; ii++)
		{
			if (db.particles[i]->sub_bv[ii]->size > largest_Psub)
				largest_Psub = db.particles[i]->sub_bv[ii]->size;
		}
	}

	//Setting maximum values for boundaries and sub boundaries
	for (int i = 0; i < db.number_boundaries; i++)
	{
		if (db.boundaries[i]->bv->size > largest_B)
			largest_B = db.boundaries[i]->bv->size;
		for (int ii = 0; ii < db.boundaries[i]->n_sub_bv; ii++)
		{
			if (db.boundaries[i]->sub_bv[ii]->size > largest_Bsub)
				largest_Bsub = db.boundaries[i]->sub_bv[ii]->size;
		}
	}

	//Setting maximum values for body geometries and sub body geometries
	for (int i = 0; i < db.number_body_geometries; i++)
	{
		if (db.body_geometries[i]->bv->size > largest_BO)
			largest_BO = db.body_geometries[i]->bv->size;
		for (int ii = 0; ii < db.body_geometries[i]->n_items; ii++)
		{
			if (db.body_geometries[i]->ptr_geom[ii]->bv->size > largest_BOsub)
				largest_BOsub = db.body_geometries[i]->ptr_geom[ii]->bv->size;
		}
	}

	float largest_P_BO = largest_P;
	if (largest_P_BO < largest_BO)
		largest_P_BO = largest_BO;

	//Domain of interest
	if (user_set_boundings == false)
	{
		//Setting boundaries - using nodal information of overall system at reference configuration
		for (int i = 0; i < db.number_nodes; i++)
		{
			if ((float)db.nodes[i]->ref_coordinates[0] > max_x)
				max_x = (float)db.nodes[i]->ref_coordinates[0];
			if ((float)db.nodes[i]->ref_coordinates[0] < min_x)
				min_x = (float)db.nodes[i]->ref_coordinates[0];
			if ((float)db.nodes[i]->ref_coordinates[1] > max_y)
				max_y = (float)db.nodes[i]->ref_coordinates[1];
			if ((float)db.nodes[i]->ref_coordinates[1] < min_y)
				min_y = (float)db.nodes[i]->ref_coordinates[1];
			if ((float)db.nodes[i]->ref_coordinates[2] > max_z)
				max_z = (float)db.nodes[i]->ref_coordinates[2];
			if ((float)db.nodes[i]->ref_coordinates[2] < min_z)
				min_z = (float)db.nodes[i]->ref_coordinates[2];
		}
		max_x += largest_P_BO;
		max_y += largest_P_BO;
		max_z += largest_P_BO;
		min_x -= largest_P_BO;
		min_y -= largest_P_BO;
		min_z -= largest_P_BO;
	}
	

	//Allocation - particle-particle contact
	if (db.number_particles != 0)
	{
		contactPP_list = DBG_NEW vector<ContactParticleParticle*>[db.number_particles];
		for (int i = 0; i < db.number_particles; i++)
		{
			contactPP_list[i].clear();
			contactPP_list[i].reserve(10);
		}
	}

	//Allocation - particle-boundary contact
	if (db.number_particles != 0)
	{
		contactPB_list = DBG_NEW vector<ContactParticleBoundary*>[db.number_particles];
		for (int i = 0; i < db.number_particles; i++)
		{
			contactPB_list[i].clear();
			contactPB_list[i].reserve(10);
		}
	}

	//Allocation - body-body contact
	if (db.number_body_geometries != 0)
	{
		contactBOBO_list = DBG_NEW vector<ContactBodyBody*>[db.number_body_geometries];
		for (int i = 0; i < db.number_body_geometries; i++)
		{
			contactBOBO_list[i].clear();
			contactBOBO_list[i].reserve(10);
		}
	}

	//Allocation - particle-body contact
	if (db.number_particles != 0)
	{
		contactPBO_list = DBG_NEW vector<ContactParticleBody*>[db.number_particles];
		for (int i = 0; i < db.number_particles; i++)
		{
			contactPBO_list[i].clear();
			contactPBO_list[i].reserve(10);
		}
	}

	//LinkedCells ou VerletLinkedCells
	if (method == 2 || method == 4)
	{
		//Maximum expected size for the grid construction
		float largest_PP_PB_PBO = largest_P;
		if (largest_PP_PB_PBO < largest_B)
			largest_PP_PB_PBO = largest_B;
		if (largest_PP_PB_PBO < largest_BO)
			largest_PP_PB_PBO = largest_BO;
		//Maximum expected size for the sub grid construction
		float largest_PP_PB_PBO_sub = largest_Psub;
		if (largest_PP_PB_PBO_sub < largest_Bsub)
			largest_PP_PB_PBO_sub = largest_Bsub;
		if (largest_PP_PB_PBO_sub < largest_BOsub)
			largest_PP_PB_PBO_sub = largest_BOsub;

		//Allocation
		cells_P_sub = new LinkedCells[db.number_particles];
		cells_B_sub = new LinkedCells[db.number_boundaries];
		cells_BO_sub = new LinkedCells[db.number_body_geometries];

		//Particles
		cells_P.PreCalc(largest_PP_PB_PBO);
		//Boundaries
		cells_B.PreCalc(largest_PP_PB_PBO);
		//Body geometries
		cells_BO.PreCalc(largest_PP_PB_PBO);
		//Sub particles
		for (int i = 0; i < db.number_particles; i++)
			cells_P_sub[i].PreCalc(largest_PP_PB_PBO_sub);
		//Sub boundaries
		for (int i = 0; i < db.number_boundaries; i++)
			cells_B_sub[i].PreCalc(largest_PP_PB_PBO_sub);
		//Body Geometries
		for (int i = 0; i < db.number_body_geometries; i++)
			cells_BO_sub[i].PreCalc(largest_PP_PB_PBO_sub);
	}

	//Verlet ou VerletLinkedCells
	if (method == 3 || method == 4)
	{
		verlet_table.AllocTables();
		if (method == 4)
			verlet_table.combine_linked_cells = true;
		else
			verlet_table.combine_linked_cells = false;
	}

}

bool GeneralContactSearch::Read(FILE *f)
{
	char s[1000];
	
	fscanf(f, "%s", s);
	if (!strcmp(s, "Method"))
	{
		fscanf(f, "%s", s);
		bool method_set = false;
		if (!strcmp(s, "AlltoAll"))
		{
			method = 1;
			method_set = true;
		}
		if (!strcmp(s, "LinkedCells"))
		{
			method = 2;
			method_set = true;
		}
		if (!strcmp(s, "Verlet"))
		{
			method = 3;
			method_set = true;
			verlet_table.Read(f);
		}
		if (!strcmp(s, "VerletLinkedCells"))
		{
			method = 4;
			method_set = true;
			verlet_table.Read(f);
		}
		
		if (method_set == false)
			return false;
	}
	else
		return false;


	//Salva a posição (stream)
	fpos_t pos;
	fgetpos(f, &pos);
	fscanf(f, "%s", s);
	if (!strcmp(s, "BVFactor"))
	{
		fscanf(f, "%s", s);
		
		inc_len_factor = (float)atof(s);
	}
	else
		fsetpos(f, &pos);

	//Salva a posição (stream)
	fgetpos(f, &pos);
	fscanf(f, "%s", s);
	if (!strcmp(s, "Domain"))
	{
		user_set_boundings = true;
		fscanf(f, "%s", s);
		min_x = (float)atof(s);
		fscanf(f, "%s", s);
		max_x = (float)atof(s);
		fscanf(f, "%s", s);
		min_y = (float)atof(s);
		fscanf(f, "%s", s);
		max_y = (float)atof(s);
		fscanf(f, "%s", s);
		min_z = (float)atof(s);
		fscanf(f, "%s", s);
		max_z = (float)atof(s);
	}
	else
		fsetpos(f, &pos);

	//Salva a posição (stream)
	fgetpos(f, &pos);
	fscanf(f, "%s", s);
	if (!strcmp(s, "BoolTable"))
		bool_table.Read(f);
	else
	{
		fsetpos(f, &pos);
		bool_table.SetDefault(true);
	}

	return true;
}
void GeneralContactSearch::Write(FILE *f)
{
	char s[30];
	switch (method)
	{
		case 1:
			sprintf(s, "AlltoAll");
			break;
		case 2:
			sprintf(s, "LinkedCells");
			break;
		case 3:
			sprintf(s, "Verlet");
			break;
		case 4:
			sprintf(s, "VerletLinkedCells");
			break;
		default:
			sprintf(s, "");
			break;
	}
	
	fprintf(f, "\nGeneralContactSearch\tMethod\t%s\t", s);
	if (method == 3 || method == 4)
		verlet_table.Write(f);
	fprintf(f,"BVFactor\t%.6e\t", inc_len_factor);
	//if (user_set_boundings)
	{
		fprintf(f, "\nDomain\n%.6e\t%.6e\n%.6e\t%.6e\n%.6e\t%.6e\n", min_x, max_x, min_y, max_y, min_z, max_z);
	}
	bool_table.Write(f);
	fprintf(f, "\n");
}

void GeneralContactSearch::ContactSearch()
{
	//Methods:
	//1 - All to All contact search
	//2 - Linked Cell contact search
	//3 - Verlet contact search

	if (plot_solution_times)
		db.myprintf("\nCONTACT SEARCH TIME CONTROL:\n");

	switch (method)
	{
	case 1:
		AlltoAll();
		break;
	case 2:
		LinkedCellsMethod();
		break;
	case 3:
		VerletMethod();
		break;
	case 4:
		VerletMethod();
		break;
	}

} 

void GeneralContactSearch::AlltoAll()
{
	//Solution time evaluation
	high_resolution_clock::time_point t_begin = high_resolution_clock::now();
	
	//All particle-particle contacts are set as unnactivated
	for (int i = 0; i < db.number_particles; i++)
		for (int cont = 0; cont < contactPP_list[i].size(); cont++)
			contactPP_list[i][cont]->cur_active = false;

	//All particle-boundary contacts are set as unnactivated
	for (int i = 0; i < db.number_particles; i++)
		for (int cont = 0; cont < contactPB_list[i].size(); cont++)
			contactPB_list[i][cont]->cur_active = false;

	//All body-body contacts are set as unnactivated
	for (int i = 0; i < db.number_body_geometries; i++)
		for (int cont = 0; cont < contactBOBO_list[i].size(); cont++)
			contactBOBO_list[i][cont]->cur_active = false;

	//All particle-body contacts are set as unnactivated
	for (int i = 0; i < db.number_particles; i++)
		for (int cont = 0; cont < contactPBO_list[i].size(); cont++)
			contactPBO_list[i][cont]->cur_active = false;
	
	int n_collisiondetection = 0;
#pragma omp parallel
	{
#pragma omp for schedule(dynamic)  reduction(+:n_collisiondetection)
		//All to all search - particle - particle
		for (int i = 0; i < db.number_particles; i++)
		{
			for (int j = i + 1; j < db.number_particles; j++)
			{
				n_collisiondetection++;
				//If collision is detected, continue the search procedure
				if (CollisionDetection(db.particles[i]->bv, db.particles[j]->bv))
				{
					//Searching subdivisions of particles
					for (int ii = 0; ii < db.particles[i]->n_sub_bv; ii++)
					{
						for (int jj = 0; jj < db.particles[j]->n_sub_bv; jj++)
						{
							n_collisiondetection++;
							if (CollisionDetection(db.particles[i]->sub_bv[ii], db.particles[j]->sub_bv[jj]))
							{
								//Search for already-existing contact
								int index = CheckIfContactParticleParticleExists(i, j, ii, jj);
								//New contact detected - creation of a new object
								if (index == -1)
								{
									ContactParticleParticle* pair;
									if (typeid(*db.particles[i]) == typeid(Polyhedron) && typeid(*db.particles[j]) == typeid(Polyhedron))
									{
										pair = new ContactPolyhedronPolyhedron();
										pair->index1 = i;
										pair->index2 = j;
										pair->sub_index1 = ii;
										pair->sub_index2 = jj;
										pair->cur_active = true;
										pair->PreCalc();
										contactPP_list[i].push_back(pair);
									}
									if (typeid(*db.particles[i]) == typeid(VEMPolyhedron) && typeid(*db.particles[j]) == typeid(VEMPolyhedron))
									{
										pair = new ContactVEMPolyhedronVEMPolyhedron();
										pair->index1 = i;
										pair->index2 = j;
										pair->sub_index1 = ii;
										pair->sub_index2 = jj;
										pair->cur_active = true;
										pair->PreCalc();
										contactPP_list[i].push_back(pair);
									}
									//outros tipos de pares de colisão entre diferentes tipos de partículas
								}
								//Previously existing contact is made active
								else
								{
									contactPP_list[i][index]->cur_active = true;
								}
									
							}

						}
					}
				}
			}
		}
	}

#pragma omp parallel
	{
#pragma omp for schedule(dynamic)  reduction(+:n_collisiondetection)
		//All to all search - particle - boundary
		for (int i = 0; i < db.number_particles; i++)
		{
			for (int j = 0; j < db.number_boundaries; j++)
			{
				n_collisiondetection++;
				//If collision is detected, continue the search procedure
				if (CollisionDetection(db.particles[i]->bv, db.boundaries[j]->bv))
				{
					//Searching subdivisions of particles
					for (int ii = 0; ii < db.particles[i]->n_sub_bv; ii++)
					{
						for (int jj = 0; jj < db.boundaries[j]->n_sub_bv; jj++)
						{
							n_collisiondetection++;
							if (CollisionDetection(db.particles[i]->sub_bv[ii], db.boundaries[j]->sub_bv[jj]))
							{
								//Search for already-existing contact
								int index = CheckIfContactParticleBoundaryExists(i, j, ii, jj);
								//New contact detected - creation of a new object
								if (index == -1)
								{
									ContactParticleBoundary* pair;
									if (typeid(*db.particles[i]) == typeid(Polyhedron) && typeid(*db.boundaries[j]) == typeid(STLBoundary))
									{
										pair = new ContactPolyhedronSTLBoundary();
										pair->index1 = i;
										pair->index2 = j;
										pair->sub_index1 = ii;
										pair->sub_index2 = jj;
										pair->cur_active = true;
										pair->PreCalc();
										contactPB_list[i].push_back(pair);
									}
									if (typeid(*db.particles[i]) == typeid(VEMPolyhedron) && typeid(*db.boundaries[j]) == typeid(STLBoundary))
									{
										pair = new ContactVEMPolyhedronSTLBoundary();
										pair->index1 = i;
										pair->index2 = j;
										pair->sub_index1 = ii;
										pair->sub_index2 = jj;
										pair->cur_active = true;
										pair->PreCalc();
										contactPB_list[i].push_back(pair);
									}
									//outros tipos de pares de colisão entre diferentes tipos de partículas e contornos
								}
								//Previously existing contact is made active
								else
								{
									contactPB_list[i][index]->cur_active = true;
								}

							}

						}
					}
				}
			}
		}
	}

#pragma omp parallel
	{
#pragma omp for schedule(dynamic)  reduction(+:n_collisiondetection)
		//All to all search - body - body
		for (int i = 0; i < db.number_body_geometries; i++)
		{
			for (int j = i + 1; j < db.number_body_geometries; j++)
			{
				n_collisiondetection++;
				//If collision is detected, continue the search procedure
				if (CollisionDetection(db.body_geometries[i]->bv, db.body_geometries[j]->bv))
				{
					//Searching subdivisions of particles
					for (int ii = 0; ii < db.body_geometries[i]->n_items; ii++)
					{
						for (int jj = 0; jj < db.body_geometries[j]->n_items; jj++)
						{
							n_collisiondetection++;
							if (CollisionDetection(db.body_geometries[i]->ptr_geom[ii]->bv, db.body_geometries[j]->ptr_geom[jj]->bv))
							{
								//Search for already-existing contact
								int index = CheckIfContactBodyBodyExists(i, j, ii, jj);
								//New contact detected - creation of a new object
								if (index == -1)
								{
									ContactBodyBody* pair;
									if (typeid(*db.body_geometries[i]->ptr_geom[ii]) == typeid(SECylinder) && typeid(*db.body_geometries[j]->ptr_geom[jj]) == typeid(SECylinder))
									{
										pair = new ContactSECylinderSECylinder();
										pair->index1 = i;
										pair->index2 = j;
										pair->sub_index1 = ii;
										pair->sub_index2 = jj;
										pair->cur_active = true;
										pair->PreCalc();
										contactBOBO_list[i].push_back(pair);
									}
									if (typeid(*db.body_geometries[i]->ptr_geom[ii]) == typeid(ArcExtrusion) && typeid(*db.body_geometries[j]->ptr_geom[jj]) == typeid(ArcRevolution))
									{
										pair = new ContactArcExtrusionArcRevolution();
										pair->index1 = i;
										pair->index2 = j;
										pair->sub_index1 = ii;
										pair->sub_index2 = jj;
										pair->cur_active = true;
										pair->PreCalc();
										contactBOBO_list[i].push_back(pair);
									}
									if (typeid(*db.body_geometries[i]->ptr_geom[ii]) == typeid(ArcRevolution) && typeid(*db.body_geometries[j]->ptr_geom[jj]) == typeid(ArcExtrusion))
									{
										pair = new ContactArcExtrusionArcRevolution();
										pair->invert = true;
										pair->index1 = i;
										pair->index2 = j;
										pair->sub_index1 = ii;
										pair->sub_index2 = jj;
										pair->cur_active = true;
										pair->PreCalc();
										contactBOBO_list[i].push_back(pair);
									}
									//outros tipos de pares de colisão entre diferentes tipos de geometria de bodies
								}
								//Previously existing contact is made active
								else
								{
									contactBOBO_list[i][index]->cur_active = true;
								}

							}

						}
					}
				}
			}
		}
	}

#pragma omp parallel
	{
#pragma omp for schedule(dynamic)  reduction(+:n_collisiondetection)
		//All to all search - particle - body
		for (int i = 0; i < db.number_particles; i++)
		{
			for (int j = 0; j < db.number_body_geometries; j++)
			{
				n_collisiondetection++;
				//If collision is detected, continue the search procedure
				if (CollisionDetection(db.particles[i]->bv, db.body_geometries[j]->bv))
				{
					//Searching subdivisions of particles
					for (int ii = 0; ii < db.particles[i]->n_sub_bv; ii++)
					{
						for (int jj = 0; jj < db.body_geometries[j]->n_items; jj++)
						{
							n_collisiondetection++;
							if (CollisionDetection(db.particles[i]->sub_bv[ii], db.body_geometries[j]->ptr_geom[jj]->bv))
							{
								//Search for already-existing contact
								int index = CheckIfContactParticleBodyExists(i, j, ii, jj);
								//New contact detected - creation of a new object
								if (index == -1)
								{
									ContactParticleBody* pair;
									if (typeid(*db.particles[i]) == typeid(Polyhedron) && typeid(*db.body_geometries[j]->ptr_geom[jj]) == typeid(ArcExtrusion))
									{
										pair = new ContactPolyhedronArcExtrusion();
										pair->index1 = i;
										pair->index2 = j;
										pair->sub_index1 = ii;
										pair->sub_index2 = jj;
										pair->cur_active = true;
										pair->PreCalc();
										contactPBO_list[i].push_back(pair);
									}
									
									//outros tipos de pares de colisão entre diferentes tipos de partículas e bodies
								}
								//Previously existing contact is made active
								else
								{
									contactPBO_list[i][index]->cur_active = true;
								}
							}

						}
					}
				}
			}
		}
	}
	global_n_collisiondetection = n_collisiondetection;
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(high_resolution_clock::now() - t_begin).count();
	if (plot_solution_times)
		db.myprintf("Main loops total: \t%.6f sec\n",duration / 1e6);
	duration_collision_detection = duration;

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//To implement (other combinations e.g.: surface - particle, surface set - surface set, etc). Create new search lists and new contact vectors
}
 
void GeneralContactSearch::VerletMethod()
{
	//Solution time evaluation
	high_resolution_clock::time_point t_begin = high_resolution_clock::now();

	//All particle-particle contacts are set as unnactivated
	for (int i = 0; i < db.number_particles; i++)
		for (int cont = 0; cont < contactPP_list[i].size(); cont++)
			contactPP_list[i][cont]->cur_active = false;

	//All particle-boundary contacts are set as unnactivated
	for (int i = 0; i < db.number_particles; i++)
		for (int cont = 0; cont < contactPB_list[i].size(); cont++)
			contactPB_list[i][cont]->cur_active = false;

	//All body-body contacts are set as unnactivated
	for (int i = 0; i < db.number_body_geometries; i++)
		for (int cont = 0; cont < contactBOBO_list[i].size(); cont++)
			contactBOBO_list[i][cont]->cur_active = false;

	//All body-body contacts are set as unnactivated
	for (int i = 0; i < db.number_particles; i++)
		for (int cont = 0; cont < contactPBO_list[i].size(); cont++)
			contactPBO_list[i][cont]->cur_active = false;

	if (verlet_table.combine_linked_cells)
	{
		//Solution time evaluation
		high_resolution_clock::time_point t_partial = high_resolution_clock::now();
		UpdateCells();
		auto duration = std::chrono::duration_cast<std::chrono::microseconds>(high_resolution_clock::now() - t_partial).count();
		if (plot_solution_times)
			db.myprintf("Refresh Cells time: \t%.6f sec\n", duration / 1e6);
		duration_linkedcells = duration;
	}
	

	if (verlet_table.acc_steps == verlet_table.n_sampling)
	{
		//Solution time evaluation
		high_resolution_clock::time_point t_partial = high_resolution_clock::now();
		verlet_table.RefreshVerletTables();
		verlet_table.SetSamplingRate();
		verlet_table.acc_steps = 0;
		auto duration = std::chrono::duration_cast<std::chrono::microseconds>(high_resolution_clock::now() - t_partial).count();
		if (plot_solution_times)
		{
			db.myprintf("Refresh Verlet tables time: \t%.6f sec.\n", duration / 1e6);
			if (verlet_table.n_sampling != 1)
				db.myprintf("Verlet Sampling : %d\n", verlet_table.n_sampling);
		}
		duration_verlet = duration;

		int n_collisiondetection = 0;
#pragma omp parallel
		{
#pragma omp for schedule(dynamic)  reduction(+:n_collisiondetection)
			for (int i = 0; i < db.number_particles; i++)
			{
				for (int cont = 0; cont < (int)verlet_table.neighborhood_PP[i].size(); cont++)
				{
					int pi = i;
					int pj = verlet_table.neighborhood_PP[i][cont][0];
					n_collisiondetection++;
					if (CollisionDetection(db.particles[pi]->bv, db.particles[pj]->bv))
					{
						int pii = verlet_table.neighborhood_PP[i][cont][1];
						int pjj = verlet_table.neighborhood_PP[i][cont][2];
						n_collisiondetection++;
						if (CollisionDetection(db.particles[pi]->sub_bv[pii], db.particles[pj]->sub_bv[pjj]))
						{
							//Search for already-existing contact
							int index = CheckIfContactParticleParticleExists(pi, pj, pii, pjj);
							//New contact detected - creation of a new object
							if (index == -1)
							{
								ContactParticleParticle* pair;
								if (typeid(*db.particles[pi]) == typeid(Polyhedron) && typeid(*db.particles[pj]) == typeid(Polyhedron))
								{
									pair = new ContactPolyhedronPolyhedron();
									pair->index1 = pi;
									pair->index2 = pj;
									pair->sub_index1 = pii;
									pair->sub_index2 = pjj;
									pair->cur_active = true;
									pair->PreCalc();
									contactPP_list[i].push_back(pair);
								}
								if (typeid(*db.particles[pi]) == typeid(VEMPolyhedron) && typeid(*db.particles[pj]) == typeid(VEMPolyhedron))
								{
									pair = new ContactVEMPolyhedronVEMPolyhedron();
									pair->index1 = pi;
									pair->index2 = pj;
									pair->sub_index1 = pii;
									pair->sub_index2 = pjj;
									pair->cur_active = true;
									pair->PreCalc();
									contactPP_list[i].push_back(pair);
								}
							}
							//Previously existing contact is made active
							else
								contactPP_list[i][index]->cur_active = true;
						}
					}
				}
			}
		}

#pragma omp parallel
		{
#pragma omp for schedule(dynamic)  reduction(+:n_collisiondetection)
			for (int i = 0; i < db.number_particles; i++)
			{
				for (int cont = 0; cont < (int)verlet_table.neighborhood_PB[i].size(); cont++)
				{
					int pi = i;
					int bj = verlet_table.neighborhood_PB[i][cont][0];
					n_collisiondetection++;
					if (CollisionDetection(db.particles[pi]->bv, db.boundaries[bj]->bv))
					{
						int pii = verlet_table.neighborhood_PB[i][cont][1];
						int bjj = verlet_table.neighborhood_PB[i][cont][2];
						n_collisiondetection++;
						if (CollisionDetection(db.particles[pi]->sub_bv[pii], db.boundaries[bj]->sub_bv[bjj]))
						{
							//Search for already-existing contact
							int index = CheckIfContactParticleBoundaryExists(pi, bj, pii, bjj);
							//New contact detected - creation of a new object
							if (index == -1)
							{
								ContactParticleBoundary* pair;
								if (typeid(*db.particles[pi]) == typeid(Polyhedron) && typeid(*db.boundaries[bj]) == typeid(STLBoundary))
								{
									pair = new ContactPolyhedronSTLBoundary();
									pair->index1 = pi;
									pair->index2 = bj;
									pair->sub_index1 = pii;
									pair->sub_index2 = bjj;
									pair->cur_active = true;
									pair->PreCalc();
									contactPB_list[i].push_back(pair);
								}
								if (typeid(*db.particles[pi]) == typeid(VEMPolyhedron) && typeid(*db.boundaries[bj]) == typeid(STLBoundary))
								{
									pair = new ContactVEMPolyhedronSTLBoundary();
									pair->index1 = pi;
									pair->index2 = bj;
									pair->sub_index1 = pii;
									pair->sub_index2 = bjj;
									pair->cur_active = true;
									pair->PreCalc();
									contactPB_list[i].push_back(pair);
								}
							}
							//Previously existing contact is made active
							else
								contactPB_list[i][index]->cur_active = true;
						}
					}
				}
			}
		}

#pragma omp parallel
		{
#pragma omp for schedule(dynamic)  reduction(+:n_collisiondetection)
			for (int i = 0; i < db.number_body_geometries; i++)
			{
				for (int cont = 0; cont < (int)verlet_table.neighborhood_BOBO[i].size(); cont++)
				{
					int pi = i;
					int pj = verlet_table.neighborhood_BOBO[i][cont][0];
					n_collisiondetection++;
					if (CollisionDetection(db.body_geometries[pi]->bv, db.body_geometries[pj]->bv))
					{
						int pii = verlet_table.neighborhood_BOBO[i][cont][1];
						int pjj = verlet_table.neighborhood_BOBO[i][cont][2];
						n_collisiondetection++;
						if (CollisionDetection(db.body_geometries[pi]->ptr_geom[pii]->bv, db.body_geometries[pj]->ptr_geom[pjj]->bv))
						{
							//Search for already-existing contact
							int index = CheckIfContactBodyBodyExists(pi, pj, pii, pjj);
							//New contact detected - creation of a new object
							if (index == -1)
							{
								ContactBodyBody* pair;
								if (typeid(*db.body_geometries[pi]->ptr_geom[pii]) == typeid(SECylinder) && typeid(*db.body_geometries[pj]->ptr_geom[pjj]) == typeid(SECylinder))
								{
									pair = new ContactSECylinderSECylinder();
									pair->index1 = pi;
									pair->index2 = pj;
									pair->sub_index1 = pii;
									pair->sub_index2 = pjj;
									pair->cur_active = true;
									pair->PreCalc();
									contactBOBO_list[i].push_back(pair);
								}
								if (typeid(*db.body_geometries[pi]->ptr_geom[pii]) == typeid(ArcExtrusion) && typeid(*db.body_geometries[pj]->ptr_geom[pjj]) == typeid(ArcRevolution))
								{
									pair = new ContactArcExtrusionArcRevolution();
									pair->index1 = pi;
									pair->index2 = pj;
									pair->sub_index1 = pii;
									pair->sub_index2 = pjj;
									pair->cur_active = true;
									pair->PreCalc();
									contactBOBO_list[i].push_back(pair);
								}
								if (typeid(*db.body_geometries[pi]->ptr_geom[pii]) == typeid(ArcRevolution) && typeid(*db.body_geometries[pj]->ptr_geom[pjj]) == typeid(ArcExtrusion))
								{
									pair = new ContactArcExtrusionArcRevolution();
									pair->invert = true;
									pair->index1 = pi;
									pair->index2 = pj;
									pair->sub_index1 = pii;
									pair->sub_index2 = pjj;
									pair->cur_active = true;
									pair->PreCalc();
									contactBOBO_list[i].push_back(pair);
								}
							}
							//Previously existing contact is made active
							else
								contactBOBO_list[i][index]->cur_active = true;
						}
					}
				}
			}
		}

#pragma omp parallel
		{
#pragma omp for schedule(dynamic)  reduction(+:n_collisiondetection)
			for (int i = 0; i < db.number_particles; i++)
			{
				for (int cont = 0; cont < (int)verlet_table.neighborhood_PBO[i].size(); cont++)
				{
					int pi = i;
					int bj = verlet_table.neighborhood_PBO[i][cont][0];
					n_collisiondetection++;
					if (CollisionDetection(db.particles[pi]->bv, db.body_geometries[bj]->bv))
					{
						int pii = verlet_table.neighborhood_PBO[i][cont][1];
						int bjj = verlet_table.neighborhood_PBO[i][cont][2];
						n_collisiondetection++;
						if (CollisionDetection(db.particles[pi]->sub_bv[pii], db.body_geometries[bj]->ptr_geom[bjj]->bv))
						{
							//Search for already-existing contact
							int index = CheckIfContactParticleBodyExists(pi, bj, pii, bjj);
							//New contact detected - creation of a new object
							if (index == -1)
							{
								ContactParticleBody* pair;
								if (typeid(*db.particles[pi]) == typeid(Polyhedron) && typeid(*db.body_geometries[bj]->ptr_geom[bjj]) == typeid(ArcExtrusion))
								{
									pair = new ContactPolyhedronArcExtrusion();
									pair->index1 = pi;
									pair->index2 = bj;
									pair->sub_index1 = pii;
									pair->sub_index2 = bjj;
									pair->cur_active = true;
									pair->PreCalc();
									contactPBO_list[i].push_back(pair);
								}

								//outros tipos de pares de colisão entre diferentes tipos de partículas e bodies
							}
							//Previously existing contact is made active
							else
							{
								contactPBO_list[i][index]->cur_active = true;
							}
							
						}
					}
				}
			}
		}
		global_n_collisiondetection = n_collisiondetection;
	}
	else
	{
		//Particle-particle
		for (int i = 0; i < db.number_particles; i++)
			for (int ii = 0; ii < contactPP_list[i].size(); ii++)
				contactPP_list[i][ii]->cur_active = contactPP_list[i][ii]->prev_active;
		//Particle-boundary
		for (int i = 0; i < db.number_particles; i++)
			for (int ii = 0; ii < contactPB_list[i].size(); ii++)
				contactPB_list[i][ii]->cur_active = contactPB_list[i][ii]->prev_active;
		//Particle-body
		for (int i = 0; i < db.number_particles; i++)
			for (int ii = 0; ii < contactPBO_list[i].size(); ii++)
				contactPBO_list[i][ii]->cur_active = contactPBO_list[i][ii]->prev_active;
		//Body-body
		for (int i = 0; i < db.number_body_geometries; i++)
			for (int ii = 0; ii < contactBOBO_list[i].size(); ii++)
				contactBOBO_list[i][ii]->cur_active = contactBOBO_list[i][ii]->prev_active;
		verlet_table.acc_steps++;
	}
	
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(high_resolution_clock::now() - t_begin).count();
	if (plot_solution_times)
		db.myprintf("Main loops total: \t%.6f sec\n", duration / 1e6);
	duration_collision_detection = duration;

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//To implement (other combinations e.g.: surface - particle, surface set - surface set, etc). Create new search lists and new contact vectors
}

//Updates teh cells to be used with LinkedCellsMethod or Verlet combined with LinkedCells 
void GeneralContactSearch::UpdateCells()
{
	//Emptying cells
	cells_P.EmptyCells();
	cells_B.EmptyCells();
	cells_BO.EmptyCells();
	//Sub particles
	for (int i = 0; i < db.number_particles; i++)
		cells_P_sub[i].EmptyCells();
	//Sub boundaries
	for (int i = 0; i < db.number_boundaries; i++)
		cells_B_sub[i].EmptyCells();
	//Sub bodies
	for (int i = 0; i < db.number_body_geometries; i++)
		cells_BO_sub[i].EmptyCells();

	//Updating the cells - particles
	for (int i = 0; i < db.number_particles; i++)
		cells_P.InsertBoundingVolume(db.particles[i]->bv);
	//Updating the cells - boundaries
	for (int i = 0; i < db.number_boundaries; i++)
		cells_B.InsertBoundingVolume(db.boundaries[i]->bv);
	//Updating the cells - bodies
	for (int i = 0; i < db.number_body_geometries; i++)
		cells_BO.InsertBoundingVolume(db.body_geometries[i]->bv);
	
	//Updating the cells - sub particles
	for (int i = 0; i < db.number_particles; i++)
		for (int ii = 0; ii < db.particles[i]->n_sub_bv; ii++)
			cells_P_sub[i].InsertBoundingVolume(db.particles[i]->sub_bv[ii]);
	//Updating the cells - sub boundaries		
	for (int i = 0; i < db.number_boundaries; i++)
		for (int ii = 0; ii < db.boundaries[i]->n_sub_bv; ii++)
			cells_B_sub[i].InsertBoundingVolume(db.boundaries[i]->sub_bv[ii]);
	//Updating the cells - sub bodies		
	for (int i = 0; i < db.number_body_geometries; i++)
		for (int ii = 0; ii < db.body_geometries[i]->n_items; ii++)
			cells_BO_sub[i].InsertBoundingVolume(db.body_geometries[i]->ptr_geom[ii]->bv);
}

void GeneralContactSearch::LinkedCellsMethod()
{
	//Solution time evaluation
	high_resolution_clock::time_point t_begin = high_resolution_clock::now();

	//All particle-particle contacts are set as unnactivated
	for (int i = 0; i < db.number_particles; i++)
		for (int cont = 0; cont < contactPP_list[i].size(); cont++)
			contactPP_list[i][cont]->cur_active = false;

	//All particle-boundary contacts are set as unnactivated
	for (int i = 0; i < db.number_particles; i++)
		for (int cont = 0; cont < contactPB_list[i].size(); cont++)
			contactPB_list[i][cont]->cur_active = false;

	//All body-body contacts are set as unnactivated
	for (int i = 0; i < db.number_body_geometries; i++)
		for (int cont = 0; cont < contactBOBO_list[i].size(); cont++)
			contactBOBO_list[i][cont]->cur_active = false;
	
	//Solution time evaluation
	high_resolution_clock::time_point t_partial = high_resolution_clock::now();
	UpdateCells();
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(high_resolution_clock::now() - t_partial).count();
	if (plot_solution_times)
		db.myprintf("Refresh Cells time: \t%.6f sec\n", duration / 1e6);
	duration_linkedcells = duration;
	
	int n_collisiondetection = 0;
#pragma omp parallel
	{
#pragma omp for schedule(dynamic)  reduction(+:n_collisiondetection)
		//LinkedCells algorithm - PP contact search
		for (int i = 0; i < db.number_particles; i++)
		{
			//Pointers defined inside the loop to avoid paralellization problems
			BoundingVolume* ptr_bv;
			BoundingVolume* ptr_subbv;
			//Cell of the particle BV
			int i_loc = db.particles[i]->bv->i_loc;
			int j_loc = db.particles[i]->bv->j_loc;
			int k_loc = db.particles[i]->bv->k_loc;
			//Loops to run on 26 neighbor cells + own central cell
			for (int cell_i = (i_loc - 1); cell_i <= (i_loc + 1); cell_i++)
			{
				for (int cell_j = (j_loc - 1); cell_j <= (j_loc + 1); cell_j++)
				{
					for (int cell_k = (k_loc - 1); cell_k <= (k_loc + 1); cell_k++)
					{
						ptr_bv = cells_P.GetCellListRoot(cell_i, cell_j, cell_k);
						for (int ob = 0; ob < cells_P.GetNObjects(cell_i, cell_j, cell_k); ob++)
						{
							//Considers only if the ID of the found particle is larger than the current one
							//This avoids self-contact occurrence and doubling contact interactions
							if (ptr_bv->associated_ID > db.particles[i]->number)
							{
								n_collisiondetection++;
								//If collision is detected, continue the search procedure to the subgrid level
								if (CollisionDetection(db.particles[i]->bv, ptr_bv))
								{
									//Searching subdivisions between particles:
									for (int ii = 0; ii < db.particles[i]->n_sub_bv; ii++)
									{
										//Cell of the particle subBV
										int i_subloc = db.particles[i]->sub_bv[ii]->i_loc;
										int j_subloc = db.particles[i]->sub_bv[ii]->j_loc;
										int k_subloc = db.particles[i]->sub_bv[ii]->k_loc;
										//Loops to run on 26 neighbor cells + own central cell
										for (int subcell_i = (i_subloc - 1); subcell_i <= (i_subloc + 1); subcell_i++)
										{
											for (int subcell_j = (j_subloc - 1); subcell_j <= (j_subloc + 1); subcell_j++)
											{
												for (int subcell_k = (k_subloc - 1); subcell_k <= (k_subloc + 1); subcell_k++)
												{
													ptr_subbv = cells_P_sub[ptr_bv->associated_ID - 1].GetCellListRoot(subcell_i, subcell_j, subcell_k);
													
													for (int subob = 0; subob < cells_P_sub[ptr_bv->associated_ID - 1].GetNObjects(subcell_i, subcell_j, subcell_k); subob++)
													{
														n_collisiondetection++;
														//If collision is detected, create the contact pair (or activate a previously created one)
														if (CollisionDetection(db.particles[i]->sub_bv[ii], ptr_subbv))
														{
															int ID = ptr_subbv->associated_ID;
															int subID = ptr_subbv->associated_sub_ID;
															//Search for already-existing contact
															int index = CheckIfContactParticleParticleExists(i, ID - 1, ii, subID - 1);
															//New contact detected - creation of a new object
															if (index == -1)
															{
																ContactParticleParticle* pair;
																if (typeid(*db.particles[i]) == typeid(Polyhedron) && typeid(*db.particles[ID - 1]) == typeid(Polyhedron))
																{
																	pair = new ContactPolyhedronPolyhedron();
																	pair->index1 = i;
																	pair->index2 = ID - 1;
																	pair->sub_index1 = ii;
																	pair->sub_index2 = subID - 1;
																	pair->cur_active = true;
																	pair->PreCalc();
																	contactPP_list[i].push_back(pair);
																}
																if (typeid(*db.particles[i]) == typeid(VEMPolyhedron) && typeid(*db.particles[ID - 1]) == typeid(VEMPolyhedron))
																{
																	pair = new ContactVEMPolyhedronVEMPolyhedron();
																	pair->index1 = i;
																	pair->index2 = ID - 1;
																	pair->sub_index1 = ii;
																	pair->sub_index2 = subID - 1;
																	pair->cur_active = true;
																	pair->PreCalc();
																	contactPP_list[i].push_back(pair);
																}
																//outros tipos de pares de colisão entre diferentes tipos de partículas
															}
															//Previously existing contact is made active
															else
															{
																contactPP_list[i][index]->cur_active = true;
															}
														}
														//Atualização do ponteiro. Obs: o último apontara NULL, mas não entrará no loop (sai antes)
														ptr_subbv = ptr_subbv->next;
													}
												}
											}
										}
									}
								}
							}

							//Atualização do ponteiro. Obs: o último apontara NULL, mas não entrará no loop (sai antes)
							ptr_bv = ptr_bv->next;
						}
					}
				}
			}
		}
	}

#pragma omp parallel
	{
#pragma omp for schedule(dynamic)  reduction(+:n_collisiondetection)
		//LinkedCells algorithm - PB contact search
		for (int i = 0; i < db.number_particles; i++)
		{
			//Pointers defined inside the loop to avoid paralellization problems
			BoundingVolume* ptr_bv;
			BoundingVolume* ptr_subbv;
			//Cell of the particle BV
			int i_loc = db.particles[i]->bv->i_loc;
			int j_loc = db.particles[i]->bv->j_loc;
			int k_loc = db.particles[i]->bv->k_loc;
			//Loops to run on 26 neighbor cells + own central cell
			for (int cell_i = (i_loc - 1); cell_i <= (i_loc + 1); cell_i++)
			{
				for (int cell_j = (j_loc - 1); cell_j <= (j_loc + 1); cell_j++)
				{
					for (int cell_k = (k_loc - 1); cell_k <= (k_loc + 1); cell_k++)
					{
						ptr_bv = cells_B.GetCellListRoot(cell_i, cell_j, cell_k);
						for (int ob = 0; ob < cells_B.GetNObjects(cell_i, cell_j, cell_k); ob++)
						{
							n_collisiondetection++;
							//If collision is detected, continue the search procedure to the subgrid level
							if (CollisionDetection(db.particles[i]->bv, ptr_bv))
							{
								//Searching subdivisions between particles/boundaries:
								for (int ii = 0; ii < db.particles[i]->n_sub_bv; ii++)
								{
									//Cell of the particle subBV 
									int i_subloc = db.particles[i]->sub_bv[ii]->i_loc;
									int j_subloc = db.particles[i]->sub_bv[ii]->j_loc;
									int k_subloc = db.particles[i]->sub_bv[ii]->k_loc;
									//Loops to run on 26 neighbor cells + own central cell
									for (int subcell_i = (i_subloc - 1); subcell_i <= (i_subloc + 1); subcell_i++)
									{
										for (int subcell_j = (j_subloc - 1); subcell_j <= (j_subloc + 1); subcell_j++)
										{
											for (int subcell_k = (k_subloc - 1); subcell_k <= (k_subloc + 1); subcell_k++)
											{
												ptr_subbv = cells_B_sub[ptr_bv->associated_ID - 1].GetCellListRoot(subcell_i, subcell_j, subcell_k);
												for (int subob = 0; subob < cells_B_sub[ptr_bv->associated_ID - 1].GetNObjects(subcell_i, subcell_j, subcell_k); subob++)
												{
													n_collisiondetection++;
													//If collision is detected, create the contact pair (or activate a previously created one)
													if (CollisionDetection(db.particles[i]->sub_bv[ii], ptr_subbv))
													{
														int ID = ptr_subbv->associated_ID;
														int subID = ptr_subbv->associated_sub_ID;
														//Search for already-existing contact
														int index = CheckIfContactParticleBoundaryExists(i, ID - 1, ii, subID - 1);
														//New contact detected - creation of a new object
														if (index == -1)
														{
															ContactParticleBoundary* pair;
															if (typeid(*db.particles[i]) == typeid(Polyhedron) && typeid(*db.boundaries[ID - 1]) == typeid(STLBoundary))
															{
																pair = new ContactPolyhedronSTLBoundary();
																pair->index1 = i;
																pair->index2 = ID - 1;
																pair->sub_index1 = ii;
																pair->sub_index2 = subID - 1;
																pair->cur_active = true;
																pair->PreCalc();
																contactPB_list[i].push_back(pair);
															}
															if (typeid(*db.particles[i]) == typeid(VEMPolyhedron) && typeid(*db.boundaries[ID - 1]) == typeid(STLBoundary))
															{
																pair = new ContactVEMPolyhedronSTLBoundary();
																pair->index1 = i;
																pair->index2 = ID - 1;
																pair->sub_index1 = ii;
																pair->sub_index2 = subID - 1;
																pair->cur_active = true;
																pair->PreCalc();
																contactPB_list[i].push_back(pair);
															}
															//outros tipos de pares de colisão entre diferentes tipos de partículas e contornos
														}
														//Previously existing contact is made active
														else
														{
															contactPB_list[i][index]->cur_active = true;
														}
													}
													//Atualização do ponteiro. Obs: o último apontara NULL, mas não entrará no loop (sai antes)
													ptr_subbv = ptr_subbv->next;
												}
											}
										}
									}
								}
							}
							//Atualização do ponteiro. Obs: o último apontara NULL, mas não entrará no loop (sai antes)
							ptr_bv = ptr_bv->next;
						}
					}
				}
			}
		}
	}


#pragma omp parallel
	{
#pragma omp for schedule(dynamic)  reduction(+:n_collisiondetection)
		//LinkedCells algorithm - BOBO contact search
		for (int i = 0; i < db.number_body_geometries; i++)
		{
			//Pointers defined inside the loop to avoid paralellization problems
			BoundingVolume* ptr_bv;
			BoundingVolume* ptr_subbv;
			//Cell of the particle BV
			int i_loc = db.body_geometries[i]->bv->i_loc;
			int j_loc = db.body_geometries[i]->bv->j_loc;
			int k_loc = db.body_geometries[i]->bv->k_loc;
			//Loops to run on 26 neighbor cells + own central cell
			for (int cell_i = (i_loc - 1); cell_i <= (i_loc + 1); cell_i++)
			{
				for (int cell_j = (j_loc - 1); cell_j <= (j_loc + 1); cell_j++)
				{
					for (int cell_k = (k_loc - 1); cell_k <= (k_loc + 1); cell_k++)
					{
						ptr_bv = cells_BO.GetCellListRoot(cell_i, cell_j, cell_k);
						for (int ob = 0; ob < cells_BO.GetNObjects(cell_i, cell_j, cell_k); ob++)
						{
							//Considers only if the ID of the found particle is larger than the current one
							//This avoids self-contact occurrence and doubling contact interactions
							if (ptr_bv->associated_ID > db.body_geometries[i]->number)
							{
								n_collisiondetection++;
								//If collision is detected, continue the search procedure to the subgrid level
								if (CollisionDetection(db.body_geometries[i]->bv, ptr_bv))
								{
									//Searching subdivisions between particles:
									for (int ii = 0; ii < db.body_geometries[i]->n_items; ii++)
									{
										//Cell of the particle subBV
										int i_subloc = db.body_geometries[i]->ptr_geom[ii]->bv->i_loc;
										int j_subloc = db.body_geometries[i]->ptr_geom[ii]->bv->j_loc;
										int k_subloc = db.body_geometries[i]->ptr_geom[ii]->bv->k_loc;
										//Loops to run on 26 neighbor cells + own central cell
										for (int subcell_i = (i_subloc - 1); subcell_i <= (i_subloc + 1); subcell_i++)
										{
											for (int subcell_j = (j_subloc - 1); subcell_j <= (j_subloc + 1); subcell_j++)
											{
												for (int subcell_k = (k_subloc - 1); subcell_k <= (k_subloc + 1); subcell_k++)
												{
													ptr_subbv = cells_BO_sub[ptr_bv->associated_ID - 1].GetCellListRoot(subcell_i, subcell_j, subcell_k);

													for (int subob = 0; subob < cells_BO_sub[ptr_bv->associated_ID - 1].GetNObjects(subcell_i, subcell_j, subcell_k); subob++)
													{
														n_collisiondetection++;
														//If collision is detected, create the contact pair (or activate a previously created one)
														if (CollisionDetection(db.body_geometries[i]->ptr_geom[ii]->bv, ptr_subbv))
														{
															int ID = ptr_subbv->associated_ID;//one-based - set on PreCalc of BodyGeometry
															int subID = ptr_subbv->associated_sub_ID;//zero-based - set on PreCalc of BodyGeometry
															//Search for already-existing contact
															int index = CheckIfContactBodyBodyExists(i, ID - 1, ii, subID);
															//New contact detected - creation of a new object
															if (index == -1)
															{
																ContactBodyBody* pair;
																if (typeid(*db.body_geometries[i]->ptr_geom[ii]) == typeid(SECylinder) && typeid(*db.body_geometries[ID - 1]->ptr_geom[subID]) == typeid(SECylinder))
																{
																	pair = new ContactSECylinderSECylinder();
																	pair->index1 = i;
																	pair->index2 = ID - 1;
																	pair->sub_index1 = ii;
																	pair->sub_index2 = subID;
																	pair->cur_active = true;
																	pair->PreCalc();
																	contactBOBO_list[i].push_back(pair);
																}
																if (typeid(*db.body_geometries[i]->ptr_geom[ii]) == typeid(ArcExtrusion) && typeid(*db.body_geometries[ID - 1]->ptr_geom[subID]) == typeid(ArcRevolution))
																{
																	pair = new ContactArcExtrusionArcRevolution();
																	pair->index1 = i;
																	pair->index2 = ID - 1;
																	pair->sub_index1 = ii;
																	pair->sub_index2 = subID;
																	pair->cur_active = true;
																	pair->PreCalc();
																	contactBOBO_list[i].push_back(pair);
																}
																if (typeid(*db.body_geometries[i]->ptr_geom[ii]) == typeid(ArcRevolution) && typeid(*db.body_geometries[ID - 1]->ptr_geom[subID]) == typeid(ArcExtrusion))
																{
																	pair = new ContactArcExtrusionArcRevolution();
																	pair->invert = true;
																	pair->index1 = i;
																	pair->index2 = ID - 1;
																	pair->sub_index1 = ii;
																	pair->sub_index2 = subID;
																	pair->cur_active = true;
																	pair->PreCalc();
																	contactBOBO_list[i].push_back(pair);
																}
																
																//outros tipos de pares de colisão entre diferentes tipos de partículas
															}
															//Previously existing contact is made active
															else
															{
																contactBOBO_list[i][index]->cur_active = true;
															}
														}
														//Atualização do ponteiro. Obs: o último apontara NULL, mas não entrará no loop (sai antes)
														ptr_subbv = ptr_subbv->next;
													}
												}
											}
										}
									}
								}
							}

							//Atualização do ponteiro. Obs: o último apontara NULL, mas não entrará no loop (sai antes)
							ptr_bv = ptr_bv->next;
						}
					}
				}
			}
		}
	}


#pragma omp parallel
	{
#pragma omp for schedule(dynamic)  reduction(+:n_collisiondetection)
		//LinkedCells algorithm - PBO contact search
		for (int i = 0; i < db.number_particles; i++)
		{
			//Pointers defined inside the loop to avoid paralellization problems
			BoundingVolume* ptr_bv;
			BoundingVolume* ptr_subbv;
			//Cell of the particle BV
			int i_loc = db.particles[i]->bv->i_loc;
			int j_loc = db.particles[i]->bv->j_loc;
			int k_loc = db.particles[i]->bv->k_loc;
			//Loops to run on 26 neighbor cells + own central cell
			for (int cell_i = (i_loc - 1); cell_i <= (i_loc + 1); cell_i++)
			{
				for (int cell_j = (j_loc - 1); cell_j <= (j_loc + 1); cell_j++)
				{
					for (int cell_k = (k_loc - 1); cell_k <= (k_loc + 1); cell_k++)
					{
						ptr_bv = cells_BO.GetCellListRoot(cell_i, cell_j, cell_k);
						for (int ob = 0; ob < cells_BO.GetNObjects(cell_i, cell_j, cell_k); ob++)
						{
							n_collisiondetection++;
							//If collision is detected, continue the search procedure to the subgrid level
							if (CollisionDetection(db.particles[i]->bv, ptr_bv))
							{
								//Searching subdivisions between particles/boundaries:
								for (int ii = 0; ii < db.particles[i]->n_sub_bv; ii++)
								{
									//Cell of the particle subBV 
									int i_subloc = db.particles[i]->sub_bv[ii]->i_loc;
									int j_subloc = db.particles[i]->sub_bv[ii]->j_loc;
									int k_subloc = db.particles[i]->sub_bv[ii]->k_loc;
									//Loops to run on 26 neighbor cells + own central cell
									for (int subcell_i = (i_subloc - 1); subcell_i <= (i_subloc + 1); subcell_i++)
									{
										for (int subcell_j = (j_subloc - 1); subcell_j <= (j_subloc + 1); subcell_j++)
										{
											for (int subcell_k = (k_subloc - 1); subcell_k <= (k_subloc + 1); subcell_k++)
											{
												ptr_subbv = cells_BO_sub[ptr_bv->associated_ID - 1].GetCellListRoot(subcell_i, subcell_j, subcell_k);
												for (int subob = 0; subob < cells_BO_sub[ptr_bv->associated_ID - 1].GetNObjects(subcell_i, subcell_j, subcell_k); subob++)
												{
													n_collisiondetection++;
													//If collision is detected, create the contact pair (or activate a previously created one)
													if (CollisionDetection(db.particles[i]->sub_bv[ii], ptr_subbv))
													{
														int ID = ptr_subbv->associated_ID;
														int subID = ptr_subbv->associated_sub_ID;
														//Search for already-existing contact
														int index = CheckIfContactParticleBodyExists(i, ID - 1, ii, subID - 1);
														//New contact detected - creation of a new object
														if (index == -1)
														{
															ContactParticleBody* pair;
															if (typeid(*db.particles[i]) == typeid(Polyhedron) && typeid(*db.body_geometries[ID - 1]->ptr_geom[subID]) == typeid(ArcExtrusion))
															{
																pair = new ContactPolyhedronArcExtrusion();
																pair->index1 = i;
																pair->index2 = ID - 1;
																pair->sub_index1 = ii;
																pair->sub_index2 = subID;
																pair->cur_active = true;
																pair->PreCalc();
																contactPBO_list[i].push_back(pair);
															}

															//outros tipos de pares de colisão entre diferentes tipos de partículas e bodies
														}
														//Previously existing contact is made active
														else
														{
															contactPBO_list[i][index]->cur_active = true;
														}
													}
													//Atualização do ponteiro. Obs: o último apontara NULL, mas não entrará no loop (sai antes)
													ptr_subbv = ptr_subbv->next;
												}
											}
										}
									}
								}
							}
							//Atualização do ponteiro. Obs: o último apontara NULL, mas não entrará no loop (sai antes)
							ptr_bv = ptr_bv->next;
						}
					}
				}
			}
		}
	}

	global_n_collisiondetection = n_collisiondetection;
	duration = std::chrono::duration_cast<std::chrono::microseconds>(high_resolution_clock::now() - t_begin).count();
	if (plot_solution_times)
		db.myprintf("Main loops total: \t%.6f sec\n", duration / 1e6);
	duration_collision_detection = duration;
}


int GeneralContactSearch::CheckIfContactParticleParticleExists(int i, int j, int ii, int jj)
{
	for (int cont = 0; cont < contactPP_list[i].size(); cont++)
	{
		if (contactPP_list[i][cont]->index2 == j && contactPP_list[i][cont]->sub_index1 == ii && contactPP_list[i][cont]->sub_index2 == jj)
			return cont;
	}
	return -1;
}

int GeneralContactSearch::CheckIfContactParticleBoundaryExists(int i, int j, int ii, int jj)
{
	for (int cont = 0; cont < contactPB_list[i].size(); cont++)
	{
		if (contactPB_list[i][cont]->index2 == j && contactPB_list[i][cont]->sub_index1 == ii && contactPB_list[i][cont]->sub_index2 == jj)
			return cont;
	}
	return -1;
}

int GeneralContactSearch::CheckIfContactBodyBodyExists(int i, int j, int ii, int jj)
{
	for (int cont = 0; cont < contactBOBO_list[i].size(); cont++)
	{
		if (contactBOBO_list[i][cont]->index2 == j && contactBOBO_list[i][cont]->sub_index1 == ii && contactBOBO_list[i][cont]->sub_index2 == jj)
			return cont;
	}
	return -1;
}

int GeneralContactSearch::CheckIfContactParticleBodyExists(int i, int j, int ii, int jj)
{
	for (int cont = 0; cont < contactPBO_list[i].size(); cont++)
	{
		if (contactPBO_list[i][cont]->index2 == j && contactPBO_list[i][cont]->sub_index1 == ii && contactPBO_list[i][cont]->sub_index2 == jj)
			return cont;
	}
	return -1;
}

void GeneralContactSearch::UpdateBoundingVolumes()
{
#pragma omp parallel
	{
		//Percorre partículas e atualiza os bounding volumes
#pragma omp for schedule(dynamic)
		for (int i = 0; i < db.number_particles; i++)
			db.particles[i]->UpdateBoundingVolumes();
	}
#pragma omp parallel
	{
		//Percorre boundaries e atualiza os bounding volumes
#pragma omp for schedule(dynamic)
		for (int i = 0; i < db.number_boundaries; i++)
			db.boundaries[i]->UpdateBoundingVolumes();
	}

#pragma omp parallel
	{
		//Percorre boundaries e atualiza os bounding volumes
#pragma omp for schedule(dynamic)
		for (int i = 0; i < db.number_body_geometries; i++)
			db.body_geometries[i]->UpdateBoundingVolumes();
	}

	//Se necessário inserir novas atualizações para superfícies, etc.
}

bool GeneralContactSearch::HaveErrors()
{
	db.gcs->ReportContact();

	//Checks for errors for active contact-contact pairs
	for (int i = 0; i < db.number_particles; i++)
		for (int cont = 0; cont < contactPP_list[i].size(); cont++)
			if (contactPP_list[i][cont]->cur_active == true)
			{
				if (contactPP_list[i][cont]->HaveErrors())
					return true;
				//db.myprintf("\n%d, %d\n", (int)contactPP_list[i][cont]->cur_active,(int)contactPP_list[i][cont]->prev_active);
				//db.myprintf("\nGeneral Contact Search has found an error in a commputed contact! Rolling over back on time...\n");
				/*if (i == 56)
				{
					if (contactPP_list[i][cont]->index2 == 290)
					{
						db.myprintf("HERE-56-290\n");
						
						for (int pi = 0; pi < contactPP_list[i][cont]->contact_pairs.size(); pi++)
						{
							if (contactPP_list[i][cont]->contact_pairs[pi]->GetActive())
							{
								db.myprintf("\nContact ---------- I1: %d I2: %d curveA: %d pointA: %d\nfaceB: %d curveB: %d pointB: %d\n",
									contactPP_list[i][cont]->contact_pairs[pi]->index1, contactPP_list[i][cont]->contact_pairs[pi]->index2,
									contactPP_list[i][cont]->contact_pairs[pi]->deg_curveA, contactPP_list[i][cont]->contact_pairs[pi]->deg_pointA,
									contactPP_list[i][cont]->contact_pairs[pi]->deg_curveB, contactPP_list[i][cont]->contact_pairs[pi]->deg_pointB);
								contactPP_list[i][cont]->contact_pairs[pi]->cd->Plot();
							}
						}
					}
					
					
				}*/
			}
	//Checks for errors for active contact-contact pairs
	for (int i = 0; i < db.number_particles; i++)
		for (int cont = 0; cont < contactPB_list[i].size(); cont++)
			if (contactPB_list[i][cont]->cur_active == true)
			{
				if (contactPB_list[i][cont]->HaveErrors())
				{
					//db.myprintf("\nError on contact: %d %d %d %d\n", contactPB_list[i][cont]->index1, contactPB_list[i][cont]->index2, 
					//	contactPB_list[i][cont]->sub_index1, contactPB_list[i][cont]->sub_index2);
					////Bounding volumes data
					//db.myprintf("\nBV's:\n");
					//db.particles[contactPB_list[i][cont]->index1]->bv->Report();
					//db.boundaries[contactPB_list[i][cont]->index2]->bv->Report();
					//db.particles[contactPB_list[i][cont]->index1]->sub_bv[contactPB_list[i][cont]->sub_index1]->Report();
					//db.boundaries[contactPB_list[i][cont]->index2]->sub_bv[contactPB_list[i][cont]->sub_index2]->Report();
					////Distance between BV's
					//float center_i[3];
					//center_i[0] = db.particles[contactPB_list[i][cont]->index1]->bv->x_center[0];
					//center_i[1] = db.particles[contactPB_list[i][cont]->index1]->bv->x_center[1];
					//center_i[2] = db.particles[contactPB_list[i][cont]->index1]->bv->x_center[2];
					//float center_j[3];
					//center_j[0] = db.boundaries[contactPB_list[i][cont]->index2]->bv->x_center[0];
					//center_j[1] = db.boundaries[contactPB_list[i][cont]->index2]->bv->x_center[1];
					//center_j[2] = db.boundaries[contactPB_list[i][cont]->index2]->bv->x_center[2];
					//float dij_SQ = (center_i[0] - center_j[0])*(center_i[0] - center_j[0]) +
					//	(center_i[1] - center_j[1])*(center_i[1] - center_j[1]) +
					//	(center_i[2] - center_j[2])*(center_i[2] - center_j[2]);
					//float size_i = db.particles[contactPB_list[i][cont]->index1]->bv->size;
					//float size_j = db.boundaries[contactPB_list[i][cont]->index2]->bv->size;
					//float cuttoff_dij =  (size_i + size_j)*(size_i + size_j);
					//db.myprintf("\ndij_SQ\t%.6e\tcuttoff_dij\t%.6e\n", dij_SQ, cuttoff_dij);
					////Distance between subBV's
					//float center_ii[3];
					//center_ii[0] = db.particles[contactPB_list[i][cont]->index1]->sub_bv[contactPB_list[i][cont]->sub_index1]->x_center[0];
					//center_ii[1] = db.particles[contactPB_list[i][cont]->index1]->sub_bv[contactPB_list[i][cont]->sub_index1]->x_center[1];
					//center_ii[2] = db.particles[contactPB_list[i][cont]->index1]->sub_bv[contactPB_list[i][cont]->sub_index1]->x_center[2];
					//float center_jj[3];
					//center_jj[0] = db.boundaries[contactPB_list[i][cont]->index2]->sub_bv[contactPB_list[i][cont]->sub_index2]->x_center[0];
					//center_jj[1] = db.boundaries[contactPB_list[i][cont]->index2]->sub_bv[contactPB_list[i][cont]->sub_index2]->x_center[1];
					//center_jj[2] = db.boundaries[contactPB_list[i][cont]->index2]->sub_bv[contactPB_list[i][cont]->sub_index2]->x_center[2];
					//float diijj_SQ = (center_ii[0] - center_jj[0])*(center_ii[0] - center_jj[0]) +
					//	(center_ii[1] - center_jj[1])*(center_ii[1] - center_jj[1]) +
					//	(center_ii[2] - center_jj[2])*(center_ii[2] - center_jj[2]);
					//float size_ii = db.particles[contactPB_list[i][cont]->index1]->sub_bv[contactPB_list[i][cont]->sub_index1]->size;
					//float size_jj = db.boundaries[contactPB_list[i][cont]->index2]->sub_bv[contactPB_list[i][cont]->sub_index2]->size;
					//float cuttoff_diijj = (size_ii + size_jj)*(size_ii + size_jj);
					//db.myprintf("\ndij_SQ\t%.6e\tcuttoff_dij\t%.6e\n", diijj_SQ, cuttoff_diijj);
					return true;
				}
			}

	//Checks for errors for active contact-contact pairs
	for (int i = 0; i < db.number_body_geometries; i++)
		for (int cont = 0; cont < contactBOBO_list[i].size(); cont++)
			if (contactBOBO_list[i][cont]->cur_active == true)
			{
				if (contactBOBO_list[i][cont]->HaveErrors())
					return true;
			}

	//Checks for errors for active contact-contact pairs
	for (int i = 0; i < db.number_particles; i++)
		for (int cont = 0; cont < contactPBO_list[i].size(); cont++)
			if (contactPBO_list[i][cont]->cur_active == true)
			{
				if (contactPBO_list[i][cont]->HaveErrors())
					return true;
			}
#pragma omp parallel
	{
		//updates variables associated with each particle (centroid and rotation tensor)
#pragma omp for schedule(dynamic)
		for (int i = 0; i < db.number_particles; i++)
			db.particles[i]->UpdateVariables();
	}

#pragma omp parallel
	{
		//updates variables associated with each boundary (centroid and rotation tensor)
#pragma omp for schedule(dynamic)
		for (int i = 0; i < db.number_boundaries; i++)
			db.boundaries[i]->UpdateVariables();
	}

#pragma omp parallel
	{
		//updates variables associated with each body (centroid and rotation tensor)
#pragma omp for schedule(dynamic)
		for (int i = 0; i < db.number_body_geometries; i++)
			db.body_geometries[i]->UpdateVariables();
	}

	//All particle-particle contacts current status are saved into a previous status for further comparison with the new check
	for (int i = 0; i < db.number_particles; i++)
		for (int cont = 0; cont < contactPP_list[i].size(); cont++)
			contactPP_list[i][cont]->prev_active = contactPP_list[i][cont]->cur_active;

	//All particle-boundary contacts current status are saved into a previous status for further comparison with the new check
	for (int i = 0; i < db.number_particles; i++)
		for (int cont = 0; cont < contactPB_list[i].size(); cont++)
			contactPB_list[i][cont]->prev_active = contactPB_list[i][cont]->cur_active;

	//All body-body contacts current status are saved into a previous status for further comparison with the new check
	for (int i = 0; i < db.number_body_geometries; i++)
		for (int cont = 0; cont < contactBOBO_list[i].size(); cont++)
			contactBOBO_list[i][cont]->prev_active = contactBOBO_list[i][cont]->cur_active;

	//All particle-body contacts current status are saved into a previous status for further comparison with the new check
	for (int i = 0; i < db.number_particles; i++)
		for (int cont = 0; cont < contactPBO_list[i].size(); cont++)
			contactPBO_list[i][cont]->prev_active = contactPBO_list[i][cont]->cur_active;

	UpdateBoundingVolumes();
	ContactSearch();
	
	//Check for non considered contacts (that should have to be considered)
	
	//particle-particle
	for (int i = 0; i < db.number_particles; i++)
	{
		for (int cont = 0; cont < contactPP_list[i].size(); cont++)
		{
			if (contactPP_list[i][cont]->cur_active == true && contactPP_list[i][cont]->prev_active == false)
			{
				if (contactPP_list[i][cont]->NightOwlContact())
				{
					db.myprintf("\nNight owl contact detected! Rolling over back on time...\n");
					return true;
				}
			}
		}
	}

	//particle-boundary
	for (int i = 0; i < db.number_particles; i++)
	{
		for (int cont = 0; cont < contactPB_list[i].size(); cont++)
		{
			if (contactPB_list[i][cont]->cur_active == true && contactPB_list[i][cont]->prev_active == false)
			{
				if (contactPB_list[i][cont]->NightOwlContact())
				{
					db.myprintf("\nNight owl contact detected! Rolling over back on time...\n");
					return true;
				}
			}
		}
	}

	//body-body
	for (int i = 0; i < db.number_body_geometries; i++)
	{
		for (int cont = 0; cont < contactBOBO_list[i].size(); cont++)
		{
			if (contactBOBO_list[i][cont]->cur_active == true && contactBOBO_list[i][cont]->prev_active == false)
			{
				if (contactBOBO_list[i][cont]->NightOwlContact())
				{
					db.myprintf("\nNight owl contact detected! Rolling over back on time...\n");
					return true;
				}
			}
		}
	}

	//particle-body
	for (int i = 0; i < db.number_particles; i++)
	{
		for (int cont = 0; cont < contactPBO_list[i].size(); cont++)
		{
			if (contactPBO_list[i][cont]->cur_active == true && contactPBO_list[i][cont]->prev_active == false)
			{
				if (contactPBO_list[i][cont]->NightOwlContact())
				{
					db.myprintf("\nNight owl contact detected! Rolling over back on time...\n");
					return true;
				}
			}
		}
	}

	return false;
}
void GeneralContactSearch::SaveConfiguration()
{
	////Result test
	////TESTE - integral de interface
	//double integral1 = 0.0;
	//double integral2 = 0.0;
	//double integral3 = 0.0;
	//double integral4 = 0.0;
	//for (int i = 0; i < db.gcs->contactPB_list[0].size(); i++)
	//{
	//	if (db.gcs->contactPB_list[0][i]->contact_pairs.size() != 0)
	//	{
	//		RigidTriangularFace_RigidTriangularFace* ptr = static_cast<RigidTriangularFace_RigidTriangularFace*>(db.gcs->contactPB_list[0][i]->contact_pairs[0]);
	//		integral1 += ptr->inter->IntegralElasticForce(ptr->cd->copy_g_n[0], ptr->cd->g_n[0]);
	//		integral2 += ptr->inter->IntegralTrapElasticForce(ptr->cd->copy_g_n[0], ptr->cd->g_n[0]);
	//		integral3 += *ptr->Wnum;
	//		integral4 += *ptr->Wteo;
	//		
	//	}

	//}
	////Informações a serem salvas
	//fprintf(ftest, "%.6e\t", db.last_converged_time+db.current_time_step);
	//fprintf(ftest, "%.6e\t%.6e\t%.6e\t%.6e\t%.6e\n", integral1, integral2, integral2 - integral1,integral3,integral4);

#pragma omp parallel
	{
		//particle-particle
#pragma omp for schedule(dynamic)
		for (int i = 0; i < db.number_particles; i++)
		{
			for (int cont = 0; cont < contactPP_list[i].size(); cont++)
			{
				if (contactPP_list[i][cont]->cur_active == true)
					contactPP_list[i][cont]->SaveConfiguration();
				else
				{
					for (int ii = 0; ii < (int)contactPP_list[i][cont]->contact_pairs.size(); ii++)
					{
						contactPP_list[i][cont]->contact_pairs[ii]->cd->copy_return_value[0] = 2;
						//contactPP_list[i][cont]->contact_pairs[ii]->cd->return_value[0] = 2;
					}
						
				}
			}
		}
	}
#pragma omp parallel
	{
		//particle-boundary
#pragma omp for schedule(dynamic)
		for (int i = 0; i < db.number_particles; i++)
		{
			for (int cont = 0; cont < contactPB_list[i].size(); cont++)
			{
				if (contactPB_list[i][cont]->cur_active == true)
					contactPB_list[i][cont]->SaveConfiguration();
				else
				{
					for (int ii = 0; ii < (int)contactPB_list[i][cont]->contact_pairs.size(); ii++)
					{
						contactPB_list[i][cont]->contact_pairs[ii]->cd->copy_return_value[0] = 2;
						//contactPB_list[i][cont]->contact_pairs[ii]->cd->return_value[0] = 2;
					}
				}
			}
		}
	}

#pragma omp parallel
	{
		//body-body
#pragma omp for schedule(dynamic)
		for (int i = 0; i < db.number_body_geometries; i++)
		{
			for (int cont = 0; cont < contactBOBO_list[i].size(); cont++)
			{
				if (contactBOBO_list[i][cont]->cur_active == true)
					contactBOBO_list[i][cont]->SaveConfiguration();
				else
				{
					contactBOBO_list[i][cont]->cd->copy_return_value[0] = 2;
					contactBOBO_list[i][cont]->cd->return_value[0] = 2;
				}
			}
		}
	}

#pragma omp parallel
	{
		//particle-body
#pragma omp for schedule(dynamic)
		for (int i = 0; i < db.number_particles; i++)
		{
			for (int cont = 0; cont < contactPBO_list[i].size(); cont++)
			{
				if (contactPBO_list[i][cont]->cur_active == true)
					contactPBO_list[i][cont]->SaveConfiguration();
				else
				{
					contactPBO_list[i][cont]->cd->copy_return_value[0] = 2;
					contactPBO_list[i][cont]->cd->return_value[0] = 2;
				}
			}
		}
	}

	if (cleanup_count == cleanup_freq)
	{
		int cleancount = 0;

		//particle-particle
		for (int i = 0; i < db.number_particles; i++)
		{
			for (int cont = 0; cont < contactPP_list[i].size(); cont++)
			{
				if (contactPP_list[i][cont]->cur_active == false)
				{
					cleancount++;
					delete contactPP_list[i][cont];
					contactPP_list[i].erase(contactPP_list[i].begin() + cont);
				}
			}
			/*contactPP_list[i].shrink_to_fit();*/
		}
		//particle-boundary
		for (int i = 0; i < db.number_particles; i++)
		{
			for (int cont = 0; cont < contactPB_list[i].size(); cont++)
			{
				if (contactPB_list[i][cont]->cur_active == false)
				{
					cleancount++;
					delete contactPB_list[i][cont];
					contactPB_list[i].erase(contactPB_list[i].begin() + cont);
				}
			}
			/*contactPB_list[i].shrink_to_fit();*/
		}
		//body-body
		for (int i = 0; i < db.number_body_geometries; i++)
		{
			for (int cont = 0; cont < contactBOBO_list[i].size(); cont++)
			{
				if (contactBOBO_list[i][cont]->cur_active == false)
				{
					cleancount++;
					delete contactBOBO_list[i][cont];
					contactBOBO_list[i].erase(contactBOBO_list[i].begin() + cont);
				}
			}
			/*contactBOBO_list[i].shrink_to_fit();*/
		}
		//particle-body
		for (int i = 0; i < db.number_particles; i++)
		{
			for (int cont = 0; cont < contactPBO_list[i].size(); cont++)
			{
				if (contactPBO_list[i][cont]->cur_active == false)
				{
					cleancount++;
					delete contactPBO_list[i][cont];
					contactPBO_list[i].erase(contactPBO_list[i].begin() + cont);
				}
			}
			/*contactPB_list[i].shrink_to_fit();*/
		}
		db.myprintf("\nGeneral Contact Search. Cleaned contacts: %d\n", cleancount);
		cleanup_count = 0;
	}
	else
		cleanup_count++;
}
void GeneralContactSearch::MountContacts()
{
	//Solution time evaluation
	high_resolution_clock::time_point t_begin = high_resolution_clock::now();

#pragma omp parallel
	{
		//updates variables associated with each particle (centroid and rotation tensor)
#pragma omp for schedule(dynamic)
		for (int i = 0; i < db.number_particles; i++)
			db.particles[i]->UpdateVariables();
	}
	
#pragma omp parallel
	{
		//updates variables associated with each boundary (centroid and rotation tensor)
#pragma omp for schedule(dynamic)
		for (int i = 0; i < db.number_boundaries; i++)
			db.boundaries[i]->UpdateVariables();
	} 

#pragma omp parallel
	{
		//updates variables associated with each body 
#pragma omp for schedule(dynamic)
		for (int i = 0; i < db.number_body_geometries; i++)
			db.body_geometries[i]->UpdateVariables();
	}

#pragma omp parallel
	{
		//particle-particle - creates surface pairs (initial geometric information, LCP and normal Gap)
#pragma omp for schedule(dynamic)
		for (int i = 0; i < db.number_particles; i++)
		{
			for (int cont = 0; cont < contactPP_list[i].size(); cont++)
				if (contactPP_list[i][cont]->cur_active)
					contactPP_list[i][cont]->ProcessSurfacePairs();
		}
	}
	ProcessContactHierarchy();
#pragma omp parallel
	{
		//particle-boundary - creates surface pairs (initial geometric information, LCP and normal Gap)
#pragma omp for schedule(dynamic)
		for (int i = 0; i < db.number_particles; i++)
		{
			for (int cont = 0; cont < contactPB_list[i].size(); cont++)
				if (contactPB_list[i][cont]->cur_active)
					contactPB_list[i][cont]->ProcessSurfacePairs();
		}
	}
#pragma omp parallel
	{
		//body-body - creates surface pairs (initial geometric information, LCP and normal Gap)
#pragma omp for schedule(dynamic)
		for (int i = 0; i < db.number_body_geometries; i++)
		{
			for (int cont = 0; cont < contactBOBO_list[i].size(); cont++)
				if (contactBOBO_list[i][cont]->cur_active)
					contactBOBO_list[i][cont]->ProcessSurfacePair();
		}
	}
#pragma omp parallel
	{
		//particle-body - creates surface pairs (initial geometric information, LCP and normal Gap)
#pragma omp for schedule(dynamic)
		for (int i = 0; i < db.number_particles; i++)
		{
			for (int cont = 0; cont < contactPBO_list[i].size(); cont++)
				if (contactPBO_list[i][cont]->cur_active)
					contactPBO_list[i][cont]->ProcessSurfacePair();
		}
	}

#pragma omp parallel
	{
		//particle-particle
#pragma omp for schedule(dynamic)
		for (int i = 0; i < db.number_particles; i++)
			for (int cont = 0; cont < contactPP_list[i].size(); cont++)
				if (contactPP_list[i][cont]->cur_active)
					contactPP_list[i][cont]->MountContacts();
	}

#pragma omp parallel
	{
		//particle-boundary
#pragma omp for schedule(dynamic)
		for (int i = 0; i < db.number_particles; i++)
			for (int cont = 0; cont < contactPB_list[i].size(); cont++)
				if (contactPB_list[i][cont]->cur_active)
					contactPB_list[i][cont]->MountContacts();
	}

#pragma omp parallel
	{
		//body-body
#pragma omp for schedule(dynamic)
		for (int i = 0; i < db.number_body_geometries; i++)
			for (int cont = 0; cont < contactBOBO_list[i].size(); cont++)
				if (contactBOBO_list[i][cont]->cur_active)
					contactBOBO_list[i][cont]->MountContact();
	}

#pragma omp parallel
	{
		//particle-body
#pragma omp for schedule(dynamic)
		for (int i = 0; i < db.number_particles; i++)
			for (int cont = 0; cont < contactPBO_list[i].size(); cont++)
				if (contactPBO_list[i][cont]->cur_active)
					contactPBO_list[i][cont]->MountContact();
	}

	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(high_resolution_clock::now() - t_begin).count();
	duration_mount_contact = duration;
}

void GeneralContactSearch::FinalUpdateContactsExplicit(double t)
{
	//Solution time evaluation
	high_resolution_clock::time_point t_begin = high_resolution_clock::now();

#pragma omp parallel
	{
		//updates variables associated with each particle (centroid and rotation tensor)
#pragma omp for schedule(dynamic)
		for (int i = 0; i < db.number_particles; i++)
			db.particles[i]->UpdateVariables();
	}

#pragma omp parallel
	{
		//updates variables associated with each boundary (centroid and rotation tensor)
#pragma omp for schedule(dynamic)
		for (int i = 0; i < db.number_boundaries; i++)
			db.boundaries[i]->UpdateVariables();
	}

#pragma omp parallel
	{
		//updates variables associated with each body 
#pragma omp for schedule(dynamic)
		for (int i = 0; i < db.number_body_geometries; i++)
			db.body_geometries[i]->UpdateVariables();
	}


#pragma omp parallel
	{
		//particle-particle - creates surface pairs (initial geometric information, LCP and normal Gap)
#pragma omp for schedule(dynamic)
		for (int i = 0; i < db.number_particles; i++)
		{
			for (int cont = 0; cont < contactPP_list[i].size(); cont++)
				if (contactPP_list[i][cont]->cur_active)
					contactPP_list[i][cont]->FinalProcessSurfacePairsExplicit(t);
		}
	}
	ProcessContactHierarchy();
#pragma omp parallel
	{
		//particle-boundary - creates surface pairs (initial geometric information, LCP and normal Gap)
#pragma omp for schedule(dynamic)
		for (int i = 0; i < db.number_particles; i++)
		{
			for (int cont = 0; cont < contactPB_list[i].size(); cont++)
				if (contactPB_list[i][cont]->cur_active)
					contactPB_list[i][cont]->FinalProcessSurfacePairsExplicit(t);
		}
	}
#pragma omp parallel
	{
		//body-body - creates surface pairs (initial geometric information, LCP and normal Gap)
#pragma omp for schedule(dynamic)
		for (int i = 0; i < db.number_body_geometries; i++)
		{
			for (int cont = 0; cont < contactBOBO_list[i].size(); cont++)
				if (contactBOBO_list[i][cont]->cur_active)
					contactBOBO_list[i][cont]->FinalProcessSurfacePairsExplicit(t);
		}
	}
#pragma omp parallel
	{
		//particle-body - creates surface pairs (initial geometric information, LCP and normal Gap)
#pragma omp for schedule(dynamic)
		for (int i = 0; i < db.number_particles; i++)
		{
			for (int cont = 0; cont < contactPBO_list[i].size(); cont++)
				if (contactPBO_list[i][cont]->cur_active)
					contactPBO_list[i][cont]->FinalProcessSurfacePairsExplicit(t);
		}
	}
}

void GeneralContactSearch::MountContactsExplicit(double t)
{
	//Solution time evaluation
	high_resolution_clock::time_point t_begin = high_resolution_clock::now();

#pragma omp parallel
	{
		//updates variables associated with each particle (centroid and rotation tensor)
#pragma omp for schedule(dynamic)
		for (int i = 0; i < db.number_particles; i++)
			db.particles[i]->UpdateVariables();
	}

#pragma omp parallel
	{
		//updates variables associated with each boundary (centroid and rotation tensor)
#pragma omp for schedule(dynamic)
		for (int i = 0; i < db.number_boundaries; i++)
			db.boundaries[i]->UpdateVariables();
	}

#pragma omp parallel
	{
		//updates variables associated with each body 
#pragma omp for schedule(dynamic)
		for (int i = 0; i < db.number_body_geometries; i++)
			db.body_geometries[i]->UpdateVariables();
	}

#pragma omp parallel
	{
		//particle-particle - creates surface pairs (initial geometric information, LCP and normal Gap)
#pragma omp for schedule(dynamic)
		for (int i = 0; i < db.number_particles; i++)
		{
			for (int cont = 0; cont < contactPP_list[i].size(); cont++)
				if (contactPP_list[i][cont]->cur_active)
					contactPP_list[i][cont]->ProcessSurfacePairs();
		}
	}
	ProcessContactHierarchy();
#pragma omp parallel
	{
		//particle-boundary - creates surface pairs (initial geometric information, LCP and normal Gap)
#pragma omp for schedule(dynamic)
		for (int i = 0; i < db.number_particles; i++)
		{
			for (int cont = 0; cont < contactPB_list[i].size(); cont++)
				if (contactPB_list[i][cont]->cur_active)
					contactPB_list[i][cont]->ProcessSurfacePairs();
		}
	}
#pragma omp parallel
	{
		//body-body - creates surface pairs (initial geometric information, LCP and normal Gap)
#pragma omp for schedule(dynamic)
		for (int i = 0; i < db.number_body_geometries; i++)
		{
			for (int cont = 0; cont < contactBOBO_list[i].size(); cont++)
				if (contactBOBO_list[i][cont]->cur_active)
					contactBOBO_list[i][cont]->ProcessSurfacePair();
		}
	}
#pragma omp parallel
	{
		//particle-body - creates surface pairs (initial geometric information, LCP and normal Gap)
#pragma omp for schedule(dynamic)
		for (int i = 0; i < db.number_particles; i++)
		{
			for (int cont = 0; cont < contactPBO_list[i].size(); cont++)
				if (contactPBO_list[i][cont]->cur_active)
					contactPBO_list[i][cont]->ProcessSurfacePair();
		}
	}

#pragma omp parallel
	{
		//particle-particle
#pragma omp for schedule(dynamic)
		for (int i = 0; i < db.number_particles; i++)
			for (int cont = 0; cont < contactPP_list[i].size(); cont++)
				if (contactPP_list[i][cont]->cur_active)
					contactPP_list[i][cont]->MountContactsExplicit(t);
	}

#pragma omp parallel
	{
		//particle-boundary
#pragma omp for schedule(dynamic)
		for (int i = 0; i < db.number_particles; i++)
			for (int cont = 0; cont < contactPB_list[i].size(); cont++)
				if (contactPB_list[i][cont]->cur_active)
					contactPB_list[i][cont]->MountContactsExplicit(t);
	}

#pragma omp parallel
	{
		//body-body
#pragma omp for schedule(dynamic)
		for (int i = 0; i < db.number_body_geometries; i++)
			for (int cont = 0; cont < contactBOBO_list[i].size(); cont++)
				if (contactBOBO_list[i][cont]->cur_active)
					contactBOBO_list[i][cont]->MountContactExplicit(t);
	}

#pragma omp parallel
	{
		//particle-body
#pragma omp for schedule(dynamic)
		for (int i = 0; i < db.number_particles; i++)
			for (int cont = 0; cont < contactPBO_list[i].size(); cont++)
				if (contactPBO_list[i][cont]->cur_active)
					contactPBO_list[i][cont]->MountContactExplicit(t);
	}

	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(high_resolution_clock::now() - t_begin).count();
	duration_mount_contact = duration;
}

void GeneralContactSearch::MountContactsGlobal()
{
	//particle-particle
	for (int i = 0; i < db.number_particles; i++)
		for (int cont = 0; cont < contactPP_list[i].size(); cont++)
			if (contactPP_list[i][cont]->cur_active)
				contactPP_list[i][cont]->MountGlobal();

	//particle-boundary
	for (int i = 0; i < db.number_particles; i++)
		for (int cont = 0; cont < contactPB_list[i].size(); cont++)
			if (contactPB_list[i][cont]->cur_active)
				contactPB_list[i][cont]->MountGlobal();

	//body-body
	for (int i = 0; i < db.number_body_geometries; i++)
		for (int cont = 0; cont < contactBOBO_list[i].size(); cont++)
			if (contactBOBO_list[i][cont]->cur_active)
				contactBOBO_list[i][cont]->MountGlobal();

	//particle-body
	for (int i = 0; i < db.number_particles; i++)
		for (int cont = 0; cont < contactPBO_list[i].size(); cont++)
			if (contactPBO_list[i][cont]->cur_active)
				contactPBO_list[i][cont]->MountGlobal();

}

void GeneralContactSearch::MountContactsGlobalExplicit()
{
	//particle-particle
	for (int i = 0; i < db.number_particles; i++)
		for (int cont = 0; cont < contactPP_list[i].size(); cont++)
			if (contactPP_list[i][cont]->cur_active)
				contactPP_list[i][cont]->MountGlobalExplicit();

	//particle-boundary
	for (int i = 0; i < db.number_particles; i++)
		for (int cont = 0; cont < contactPB_list[i].size(); cont++)
			if (contactPB_list[i][cont]->cur_active)
				contactPB_list[i][cont]->MountGlobalExplicit();

	//body-body
	for (int i = 0; i < db.number_body_geometries; i++)
		for (int cont = 0; cont < contactBOBO_list[i].size(); cont++)
			if (contactBOBO_list[i][cont]->cur_active)
				contactBOBO_list[i][cont]->MountGlobalExplicit();

	//particle-body
	for (int i = 0; i < db.number_particles; i++)
		for (int cont = 0; cont < contactPBO_list[i].size(); cont++)
			if (contactPBO_list[i][cont]->cur_active)
				contactPBO_list[i][cont]->MountGlobalExplicit();
}

void GeneralContactSearch::SolutionStepInitialCheck()
{
	UpdateBoundingVolumes();
	ContactSearch();
}

//Saves a text file with contact data (debugging purposes)
void GeneralContactSearch::ReportContact()
{
	//double kin = 0.0;
	n_monitoring_PP = 0;
	n_active_PP = 0;
	for (int i = 0; i < db.number_particles; i++)
	{
		//kin += db.particles[i]->kinetic_energy;
		for (int cont = 0; cont < contactPP_list[i].size(); cont++)
		{
			if (contactPP_list[i][cont]->cur_active)
			{
				for (int cp = 0; cp < contactPP_list[i][cont]->contact_pairs.size(); cp++)
				{
					n_monitoring_PP++;
					contactPP_list[i][cont]->contact_pairs[cp]->Report();
					if (contactPP_list[i][cont]->contact_pairs[cp]->eligible == true)
					{
						n_active_PP++;
						//contactPP_list[i][cont]->contact_pairs[cp]->Report();
					}
						
				}
			}
		}
	}

	db.myprintf("\nNumber of monitored contacts (P-P):\t%d\n", n_monitoring_PP);
	db.myprintf("\nNumber of active contacts (P-P):\t%d\n", n_active_PP);

	n_monitoring_PB = 0;
	n_active_PB = 0;
	for (int i = 0; i < db.number_particles; i++)
	{
		for (int cont = 0; cont < contactPB_list[i].size(); cont++)
		{
			if (contactPB_list[i][cont]->cur_active)
			{
				for (int cp = 0; cp < contactPB_list[i][cont]->contact_pairs.size(); cp++)
				{
					n_monitoring_PB++;
					contactPB_list[i][cont]->contact_pairs[cp]->Report();
					if (contactPB_list[i][cont]->contact_pairs[cp]->eligible == true)
					{
						n_active_PB++;
						//contactPB_list[i][cont]->contact_pairs[cp]->Report();
					}

				}
			}
		}
	}

	db.myprintf("\nNumber of monitored contacts (P-B):\t%d\n", n_monitoring_PB);
	db.myprintf("\nNumber of active contacts (P-B):\t%d\n", n_active_PB);

	n_monitoring_BOBO = 0;
	n_active_BOBO = 0;
	for (int i = 0; i < db.number_body_geometries; i++)
	{
		for (int cont = 0; cont < contactBOBO_list[i].size(); cont++)
		{
			if (contactBOBO_list[i][cont]->cur_active)
			{
				n_monitoring_BOBO++;
				contactBOBO_list[i][cont]->Report();
				if (contactBOBO_list[i][cont]->eligible == true)
					n_active_BOBO++;	
			}
		}
	}

	db.myprintf("\nNumber of monitored contacts (BO-BO):\t%d\n", n_monitoring_BOBO);
	db.myprintf("\nNumber of active contacts (BO-BO):\t%d\n", n_active_BOBO);

	n_monitoring_PBO = 0;
	n_active_PBO = 0;
	for (int i = 0; i < db.number_particles; i++)
	{
		for (int cont = 0; cont < contactPBO_list[i].size(); cont++)
		{
			if (contactPBO_list[i][cont]->cur_active)
			{
				n_monitoring_PBO++;
				contactPBO_list[i][cont]->Report();
				if (contactPBO_list[i][cont]->eligible == true)
					n_active_PBO++;
			}
		}
	}

	db.myprintf("\nNumber of monitored contacts (P-BO):\t%d\n", n_monitoring_PBO);
	db.myprintf("\nNumber of active contacts (P-BO):\t%d\n", n_active_PBO);
}

void GeneralContactSearch::WriteVTK_XMLForces(FILE *f)
{
	for (int i = 0; i < db.number_particles; i++)
	{
		for (int cont = 0; cont < contactPP_list[i].size(); cont++)
		{
			if (contactPP_list[i][cont]->cur_active)
			{
				contactPP_list[i][cont]->WriteVTK_XMLForces(f);
			}
		}
	}

	for (int i = 0; i < db.number_particles; i++)
	{
		for (int cont = 0; cont < contactPB_list[i].size(); cont++)
		{
			if (contactPB_list[i][cont]->cur_active)
			{
				contactPB_list[i][cont]->WriteVTK_XMLForces(f);
			}
		}
	}

	for (int i = 0; i < db.number_body_geometries; i++)
	{
		for (int cont = 0; cont < contactBOBO_list[i].size(); cont++)
		{
			if (contactBOBO_list[i][cont]->cur_active)
			{
				contactBOBO_list[i][cont]->WriteVTK_XMLForces(f);
			}
		}
	}

	for (int i = 0; i < db.number_particles; i++)
	{
		for (int cont = 0; cont < contactPBO_list[i].size(); cont++)
		{
			if (contactPBO_list[i][cont]->cur_active)
			{
				contactPBO_list[i][cont]->WriteVTK_XMLForces(f);
			}
		}
	}
}

//Checa falta de dados de interfaces de contato
bool GeneralContactSearch::Check()
{
	//TODO
	//Verificar a existencia de todos as interfaces de contato necessarias para as interacoes entre particulas.
	return true;
}

//Returns the maximum time step required for contact well-resolution purposes
double GeneralContactSearch::TimeStepControl()
{
	double temp_step = db.solution[db.current_solution_number-1]->end_time;		//valor alto
	for (int i = 0; i < db.number_particles; i++)
	{
		for (int cont = 0; cont < contactPP_list[i].size(); cont++)
		{
			if (contactPP_list[i][cont]->cur_active)
			{
				//Verificação - energia cinética
				double kin1 = db.particles[contactPP_list[i][cont]->index1]->kinetic_energy;
				double kin2 = db.particles[contactPP_list[i][cont]->index2]->kinetic_energy;
				double step = contactPP_list[i][cont]->TimeStepControl(kin1 + kin2);
				if (step < temp_step)
					temp_step = step;
			}
		}
	}
	for (int i = 0; i < db.number_particles; i++)
	{
		for (int cont = 0; cont < contactPB_list[i].size(); cont++)
		{
			if (contactPB_list[i][cont]->cur_active)
			{
				//Verificação - energia cinética
				double kin1 = db.particles[contactPB_list[i][cont]->index1]->kinetic_energy;
				double step = contactPB_list[i][cont]->TimeStepControl(kin1);
				if (step < temp_step)
					temp_step = step;
			}
		}
	}
	return temp_step;
}

void GeneralContactSearch::ProcessContactHierarchy()
{
	//Testing hierarchy for PP (particle-particle)
	for (int i = 0; i < db.number_particles; i++)
	{
		for (int cont = 0; cont < contactPP_list[i].size(); cont++)
		{
			contactPP_list[i][cont]->ProcessContactHierarchy();
		}
	}
	//Testing hierarchy for PB (particle-boundary)
	for (int i = 0; i < db.number_particles; i++)
	{
		for (int cont = 0; cont < contactPB_list[i].size(); cont++)
		{
			contactPB_list[i][cont]->ProcessContactHierarchy();
		}
	}
}