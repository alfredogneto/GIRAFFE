#pragma once
#include <stdio.h>
#include <vector>
class Matrix;

/////////////////////////////////////////////////////////////////
//IDs of the elements of the Giraffe								/////
///	1 - Beam_1												/////
///	2 - Pipe_1												/////
///	3 - Shell_1												/////
///	4 - Mass_1												/////	
///	5 - SpringDashpot_1										/////	
///	6 - RigidBody_1											/////
///	7 - Solid_1												/////
/////////////////////////////////////////////////////////////////

class Element
{
public:
	int material;	//ID of the material
	int section;	//ID of the section - or pipe_sec (if for element pipe)
	int *nodes;		//Nodes global - conectividade
	int number;		//ID of the element
	int n_nodes;	//Number of nodes of the element
	int cs;			//ID of the coordinate system of the element
	char* type_name;//Name of the type of the element
	int VTK_type;
	int *VTK_nodes;//Indexing for convert numbering of the formato giraffe for the formato of the celula equivalente of the paraview

	int **DOFs;		//Indicates for the indexing of each degree of freedom, 1 or 0, active or inactive for the element in question

	int nDOFs;
	int nDOFs_free;
	int nDOFs_fixed;

	//Variables for check of the function Band
	int temp_band_free;
	int temp_band_fixed;

	int lowest_free_global_DOF;
	int highest_free_global_DOF;
	int lowest_fixed_global_DOF;
	int highest_fixed_global_DOF;
	bool free_marked;
	bool fixed_marked;

	double strain_energy;
	double kinetic_energy;
	double potential_gravitational_energy;

	Element();
	virtual ~Element();
	virtual bool Check() = 0;				//Checks inconsistencies in the element for avoid errors of execution
	virtual bool Read(FILE *f) = 0;
	virtual void Write(FILE *f) = 0;
	virtual void WriteResults(FILE *f) = 0;	//Writes file of results
	virtual void WriteMonitor(FILE *f, bool first_record, double time) = 0;	//Writes in the monitor of the element
	virtual void WriteVTK_XMLBase(std::vector<float> *float_vector) = 0;
	virtual void WriteVTK_XMLRender(FILE *f) = 0;
	
	virtual void Mount() = 0;									//Assembles elements
	virtual void MountMass() = 0;								//Assembles the mass matrix
	virtual void MountMassModal() = 0;							//Assembles the mass matrix for performing of the analysis modal
	virtual void MountDampingModal() = 0;						//Assembles the damping matrix for performing of the analysis modal
	virtual void MountElementLoads() = 0;						//Assembles loadings associated to the element
	virtual void MountDamping(bool update_rayleigh) = 0;		//Assembles the damping matrix
	virtual void TransformMatrix() = 0;							//Assembles transformation matrix of coordinates
	virtual void MountGlobal() = 0;								//Fills the contribution of the element in the matrices global

	virtual void SaveLagrange() = 0;	//Saves variables in the gauss points useful for description lagrangiana updated
	void Band(int* band_fixed, int* band_free);		//Calculates the band generated in the matrix global through the element

	virtual void PreCalc() = 0;		//Pre-calculation of variables that and done a single time in the start
	virtual void MountDyn() = 0;	//Assemblies - Newmark
	virtual void MountDynModal() = 0;	//Assemblies for analysis modal - insertion of the mass matrix and damping in the stiffness matrix for subsequent assembly global
	virtual void Zeros() = 0;		//Limpa the matrices internal of the element
};

