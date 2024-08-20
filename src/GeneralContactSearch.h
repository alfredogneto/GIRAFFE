#pragma once
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "Matrix.h"
#include "BoolTable.h"
#include "ContactParticleParticle.h"
#include "ContactParticleBoundary.h"
#include "ContactBodyBody.h"
#include "ContactParticleBody.h"

//Particle pairs:
#include "ContactPolyhedronPolyhedron.h"
#include "ContactVEMPolyhedronVEMPolyhedron.h"
//Particle-boundary pairs:
#include "ContactPolyhedronSTLBoundary.h"
#include "ContactVEMPolyhedronSTLBoundary.h"
//Geometry pairs:
#include "ContactSECylinderSECylinder.h"
#include "ContactArcExtrusionArcRevolution.h"
//Particle-body pairs:
#include "ContactPolyhedronArcExtrusion.h"

#include "LinkedCells.h"
#include "Verlet.h"

#include <vector>
#include <array>
#include <string>
#include <ctype.h>

#include <cstdio>

#include <omp.h>
#include <process.h>
#include <chrono>

using namespace std;

class GeneralContactSearch
{
public:
	int method;
	vector<ContactParticleParticle*>* contactPP_list;
	vector<ContactParticleBoundary*>* contactPB_list;
	vector<ContactBodyBody*>* contactBOBO_list;
	vector<ContactParticleBody*>* contactPBO_list;
	
	BoolTable bool_table;							//Bool table que ativa ou desativa na seq. de solutions

	GeneralContactSearch();
	~GeneralContactSearch();

	void PreCalc();									//Pre calculation (allocation,etc.)
	bool Read(FILE *f);								//Reads data
	void Write(FILE *f);							//Writes data 
	bool Check();									//Checa falta de dados de interfaces de contato
	void UpdateBoundingVolumes();					//Updates all the bounding volumes of particles, surfaces, etc.
	bool HaveErrors();								//Checks for untreated potential collisions for rolling back in time
	void SaveConfiguration();						//Saves a converged configuration data
	void MountContacts();							//Mounts all the contact local contributions detected during the global check
	void MountContactsExplicit(double t);			//Mounts all the contact loads (explicit)
	void MountContactsGlobalExplicit();				//Mounts all the contact global contributions detected during the global check (explicit)
	void FinalUpdateContactsExplicit(double t);		//Updates gap, normal and another quantities for the final check on the explicit method
	void MountContactsGlobal();						//Mounts all the contact global contributions detected during the global check
	void SolutionStepInitialCheck();				//Computes the bounding volumes on the initial of a solution step (first computing)
	void ContactSearch();							//Overall collision detection
	void AlltoAll();								//Checks all to all entities for collision detection
	void LinkedCellsMethod();						//Uses the linked cell algorithm to enhance the collision detection search in space
	void UpdateCells();								//Updates teh cells to be used with LinkedCellsMethod or Verlet combined with LinkedCells 
	void VerletMethod();							//Uses the verlet algorithm to enhance the collision detection search in space
	void ReportContact();							//Saves a text file with contact data (debugging purposes)
	void WriteVTK_XMLForces(FILE *f);				//writes forces on a Paraview post-processing file
	double TimeStepControl();						//Returns the maximum time step required for contact well-resolution purposes

	void ProcessContactHierarchy();
	//Cleaning variables
	int cleanup_count;								
	int cleanup_freq;
	//Checks is there is already a contact created for the given indexes
	int CheckIfContactParticleParticleExists(int i, int j, int ii, int jj);
	int CheckIfContactParticleBoundaryExists(int i, int j, int ii, int jj);
	int CheckIfContactBodyBodyExists(int i, int j, int ii, int jj);
	int CheckIfContactParticleBodyExists(int i, int j, int ii, int jj);
	//return values:
	//"-1" if no contact exists
	//"index corresponding to the contact" in the vector contact_list[particle index]

	//Bounding Volume global setting factors
	float inc_len_factor;
	
	//////////////////////////////////////////////////////////////////////////////////////////////
	//LinkedCells
	LinkedCells cells_P;				//Cells for particles
	LinkedCells cells_B;				//Cells for boundaries
	LinkedCells cells_BO;				//Cells for bodies
	LinkedCells* cells_P_sub;			//Cells for sub-particles
	LinkedCells* cells_B_sub;			//Cells for sub-boundaries
	LinkedCells* cells_BO_sub;			//Cells for sub-bodies
	//////////////////////////////////////////////////////////////////////////////////////////////
	//Verlet
	Verlet verlet_table;
	//////////////////////////////////////////////////////////////////////////////////////////////
	//Valid domain boundings
	bool user_set_boundings;			//Boolean flag to indicate that user input bounding limits
	float max_x;						//max x
	float min_x;						//min x
	float max_y;						//max y
	float min_y;						//min y
	float max_z;						//max z
	float min_z;						//min z


	//Post-processing variables
	int n_active_PP;
	int n_active_PB;
	int n_active_BOBO;
	int n_active_PBO;
	int n_monitoring_PP;
	int n_monitoring_PB;
	int n_monitoring_BOBO;
	int n_monitoring_PBO;

	//Variables for measuring processing times
	bool plot_solution_times;
	long long duration_verlet;
	long long duration_linkedcells;
	long long duration_collision_detection;
	long long duration_mount_contact;
	long long global_n_collisiondetection;
};

