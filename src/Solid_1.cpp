#include "Solid_1.h"
#include "Beam_1.h"
#include "Encoding.h"	

#include"Database.h"
//Global variables
extern
Database db;

Solid_1::Solid_1()
{
	strain_energy = 0.0;
	kinetic_energy = 0.0;
	potential_gravitational_energy = 0.0;

	//VTK_type = 24;
	nDOFs = 24;
	material = 0;
	section = 0;
	n_nodes = 8;
	number = 0;
	nodes = new int[n_nodes];
	VTK_nodes = new int[n_nodes];
	VTK_nodes[0] = 0;
	VTK_nodes[1] = 1;
	VTK_nodes[2] = 2;
	VTK_nodes[3] = 3;
	VTK_nodes[4] = 4;
	VTK_nodes[5] = 5;
	VTK_nodes[6] = 6;
	VTK_nodes[7] = 7;
	
	DOFs = new int *[n_nodes];
	for (int i = 0; i < n_nodes;i++)
		DOFs[i] = new int[db.number_GLs_node];

	//Routine for activate the GLS of each node of the element
	for (int i = 0; i < n_nodes; i++)
	{
		for (int j = 0; j < db.number_GLs_node; j++)
		{
			DOFs[i][j] = 0;
		}
		DOFs[i][0] = 1;
		DOFs[i][1] = 1;
		DOFs[i][2] = 1;
	}
}

Solid_1::~Solid_1()
{
	delete[] nodes;
	delete[] VTK_nodes;
	if (DOFs != NULL)
	{
		for (int i = 0; i < n_nodes; i++)
			delete[] DOFs[i];
		delete[] DOFs;
	}
}

//Checks inconsistencies in the element for avoid errors of execution
bool Solid_1::Check()
{
	return true;
}

bool Solid_1::Read(FILE *f)
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
	if (!strcmp(s, "Nodes"))
	{
		for (int n = 0; n < n_nodes; n++)
		{
			fscanf(f, "%s", s);
			nodes[n] = atoi(s);
		}
	}
	else
		return false;
	return true;
}

void Solid_1::Write(FILE *f)
{
	fprintf(f, "Solid_1\t%d\tMat\t%d\tCS\t%d\tNodes\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
		number,
		material,
		cs,
		nodes[0],
		nodes[1],
		nodes[2],
		nodes[3], 
		nodes[4], 
		nodes[5], 
		nodes[6], 
		nodes[7] );
}
//Writes file of results
void Solid_1::WriteResults(FILE *f)
{

}

void Solid_1::WriteVTK_XMLBase(std::vector<float> *float_vector)
{
	//Prints the results of the element
	int res_element = 0;
	//Prints values null for results that not fazem sentido for this element type
	for (int i = res_element; i < db.n_element_results; i++)
		float_vector->push_back(0.0);
}

void Solid_1::WriteVTK_XMLRender(FILE *f)
{
	//DOES NOTHING
}

//Writes in the monitor of the element//Writes in the monitor of the element
void Solid_1::WriteMonitor(FILE *f, bool first_record, double time)
{

}

//Assembles elements
void Solid_1::Mount()
{

}
//Assembles loadings associated to the element
void Solid_1::MountElementLoads()
{

}
//Assembles transformation matrix of coordinates
void Solid_1::TransformMatrix()
{

}
//Fills the contribution of the element in the matrices global
void Solid_1::MountGlobal()
{

}
//Saves variables in the gauss points useful for description lagrangiana updated
void Solid_1::SaveLagrange()
{

}
//Pre-calculation of variables that and done a single time in the start
void Solid_1::PreCalc()
{

}

//Assembles the mass matrix
void Solid_1::MountMass()
{

}

//Assembles the mass matrix
void Solid_1::MountMassModal()
{
	Zeros();
}

//Assembles the damping matrix for performing of the analysis modal
void Solid_1::MountDampingModal()
{
	Zeros();
}

//Assembles the damping matrix
void Solid_1::MountDamping(bool update_rayleigh)
{

}

//Assemblies - Newmark
void Solid_1::MountDyn()
{

}

//Assemblies for analysis modal - insertion of the mass matrix and damping in the stiffness matrix for subsequent assembly global
void Solid_1::MountDynModal()
{
	
}

//Zeroes matrices local of the element
void Solid_1::Zeros()
{
	kinetic_energy = 0.0;
	strain_energy = 0.0;
	potential_gravitational_energy = 0.0;
}