#include "SuperNode.h"

SuperNode::SuperNode()
{
	super_node_coordinates[0] = 0.0;
	super_node_coordinates[1] = 0.0;
	super_node_coordinates[2] = 0.0;

	n_displacement_DOFs = 0;
	n_temperature_DOFs = 0;
	n_DOFs = 0;

	ID = 0;

	ref_coordinates = NULL;
	copy_coordinates = NULL;
	displacements = NULL;
	vel = NULL;
	copy_vel = NULL;
	accel = NULL;
	copy_accel = NULL;

	temperatures = NULL;
	copy_temperatures = NULL;
	ref_temperatures = NULL;

	constraints = NULL;
	DOFs = NULL;

	alloced = false;

	n_DOFs = 0;
}

SuperNode::~SuperNode()
{
	Free();
}

//Allocates the memory for the super node
void SuperNode::Alloc(int e_displacement_DOFs, int e_temperature_DOFs)
{
	Free();
	if (alloced == false)
	{
		n_displacement_DOFs = e_displacement_DOFs;
		n_temperature_DOFs = e_temperature_DOFs;
		n_DOFs = e_displacement_DOFs + e_temperature_DOFs;

		ref_coordinates = new double[e_displacement_DOFs];
		copy_coordinates = new double[e_displacement_DOFs];
		displacements = new double[e_displacement_DOFs];
		vel = new double[e_displacement_DOFs];
		copy_vel = new double[e_displacement_DOFs];
		accel = new double[e_displacement_DOFs];
		copy_accel = new double[e_displacement_DOFs];
		
		temperatures = new double[e_temperature_DOFs];
		copy_temperatures = new double[e_temperature_DOFs];
		ref_temperatures = new double[e_temperature_DOFs];

		constraints = new bool[n_DOFs];
		DOFs = new int[n_DOFs];
		alloced = true;
		
		for (int i = 0; i < e_displacement_DOFs; i++)
		{
			ref_coordinates[i] = 0.0;
			copy_coordinates[i] = 0.0;
			displacements[i] = 0.0;
			vel[i] = 0.0;
			copy_vel[i] = 0.0;
			accel[i] = 0.0;
			copy_accel[i] = 0.0;
		}

		for (int i = 0; i < e_temperature_DOFs; i++)
		{
			temperatures[i] = 0.0;
			copy_temperatures[i] = 0.0;
			ref_temperatures[i] = 0.0;
		}

		for (int i = 0; i < n_DOFs; i++)
		{
			constraints[i] = false;	//default (free DOF)
			DOFs[i] = 0;			//Number to be set by the function:  SetGlobalDOFs (solution.cpp)
		}
	}
}
//Frees the memory for the super node
void SuperNode::Free()
{
	if (alloced == true)
	{
		delete[]ref_coordinates;
		delete[]copy_coordinates;
		delete[]displacements;
		delete[]vel;
		delete[]copy_vel;
		delete[]accel;
		delete[]copy_accel;

		delete[]temperatures;
		delete[]copy_temperatures;
		delete[]ref_temperatures;

		delete[]constraints;
		delete[]DOFs;

		alloced = false;
		n_DOFs = 0;
		n_displacement_DOFs = 0;
		n_temperature_DOFs = 0;
	}
}
//Read function
bool SuperNode::Read(FILE *f)
{
	char s[1000];
	fscanf(f, "%s", s);
	//Verifica a palavra chave "SuperNode"
	if (!strcmp(s, "SuperNode"))
	{
		fscanf(f, "%s", s);
		ID = atoi(s);
	}
	else
		return false;
	//Leitura das coordenadas
	fscanf(f, "%s", s);
	super_node_coordinates[0] = atof(s);
	fscanf(f, "%s", s);
	super_node_coordinates[1] = atof(s);
	fscanf(f, "%s", s);
	super_node_coordinates[2] = atof(s);
	
	return true;
}
//Write function
void SuperNode::Write(FILE *f)
{
	fprintf(f, "SuperNode\t%d\t%.12e\t%.12e\t%.12e\n", ID,
		super_node_coordinates[0],
		super_node_coordinates[1],
		super_node_coordinates[2]);
}
//Saves a converged configuration
void SuperNode::SaveConfiguration()
{
	//Displacements
	for (int i = 0; i < n_displacement_DOFs; i++)
	{
		copy_coordinates[i] += displacements[i];
		copy_accel[i] = accel[i];
		copy_vel[i] = vel[i];
		displacements[i] = 0.0;
	}
	//Temperatures
	for (int i = 0; i < n_temperature_DOFs; i++)
	{
		copy_temperatures[i] = temperatures[i];
	}
}
