#include "InitialCondition.h"
#include"Database.h"

//Variáveis globais
extern
Database db;

InitialCondition::InitialCondition()
{
	number = 0;
	node = 0;
	super_node = 0;
	du = Matrix(3);
	omega = Matrix(3);
	solution = 1;
}

InitialCondition::~InitialCondition()
{

}
bool InitialCondition::Read(FILE *f)
{
	char s[1000];
	fscanf(f, "%s", s);
	//Verifica a palavra chave "InitialCondition"
	if (!strcmp(s, "InitialCondition"))
	{
		fscanf(f, "%s", s);
		number = atoi(s);
	}
	else
		return false;

	fscanf(f, "%s", s);
	//Verifica a palavra chave "Node"
	if (!strcmp(s, "Node"))
	{
		fscanf(f, "%s", s);
		node = atoi(s);
	}
	else
	{
		if (!strcmp(s, "SuperNode"))
		{
			fscanf(f, "%s", s);
			super_node = atoi(s);
		}
		else
			return false;
	}

	fscanf(f, "%s", s);
	if (!strcmp(s, "DU"))
	{
		fscanf(f, "%s", s);
		du(0, 0) = atof(s);
		fscanf(f, "%s", s);
		du(1, 0) = atof(s);
		fscanf(f, "%s", s);
		du(2, 0) = atof(s);
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "OMEGA"))
	{
		fscanf(f, "%s", s);
		omega(0, 0) = atof(s);
		fscanf(f, "%s", s);
		omega(1, 0) = atof(s);
		fscanf(f, "%s", s);
		omega(2, 0) = atof(s);
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "SolutionStep"))
	{
		fscanf(f, "%s", s);
		solution = atoi(s);
	}
	else
		return false;

	return true;
}
void InitialCondition::Write(FILE *f)
{
	if (node != 0)
	fprintf(f, "InitialCondition\t%d\tNode\t%d\tDU\t%.6e\t%.6e\t%.6e\tOMEGA\t%.6e\t%.6e\t%.6e\tSolutionStep\t%d\n",
		number, node, du(0,0), du(1,0), du(2,0),omega(0,0),omega(1,0),omega(2,0),solution);
	if (super_node != 0)
		fprintf(f, "InitialCondition\t%d\tSuperNode\t%d\tDU\t%.6e\t%.6e\t%.6e\tOMEGA\t%.6e\t%.6e\t%.6e\tSolutionStep\t%d\n",
			number, node, du(0, 0), du(1, 0), du(2, 0), omega(0, 0), omega(1, 0), omega(2, 0), solution);
}

bool InitialCondition::Check()
{
	if (node != 0 && node > db.number_nodes)
		return false;
	if (super_node != 0 && super_node > db.number_super_nodes)
		return false;
	return true;
}

void InitialCondition::ComputeInitialCondition()
{
	if (db.current_solution_number == solution)
	{
		if (node != 0)
		{
			if (db.nodes[node - 1]->flag_material_description == false)
			{
				for (int index = 0; index < 3; index++)
				{
					db.nodes[node - 1]->copy_vel[index] = du(index, 0);
					db.nodes[node - 1]->copy_vel[index + 3] = omega(index, 0);
					db.nodes[node - 1]->vel[index] = du(index, 0);
					db.nodes[node - 1]->vel[index + 3] = omega(index, 0);
					//Explicit
					(*db.nodes[node - 1]->du)(index, 0) = du(index, 0);
					(*db.nodes[node - 1]->omega)(index, 0) = omega(index, 0);
				}
			}
			else
			{
				Matrix transf_du = transp(*db.nodes[node - 1]->Q0) * du;
				Matrix transf_omega = transp(*db.nodes[node - 1]->Q0) * omega;

				for (int index = 0; index < 3; index++)
				{
					db.nodes[node - 1]->copy_vel[index] = transf_du(index, 0);
					db.nodes[node - 1]->copy_vel[index + 3] = transf_omega(index, 0);
					db.nodes[node - 1]->vel[index] = transf_du(index, 0);
					db.nodes[node - 1]->vel[index + 3] = transf_omega(index, 0);
					//Explicit
					(*db.nodes[node - 1]->du)(index, 0) = transf_du(index, 0);
					(*db.nodes[node - 1]->omega)(index, 0) = transf_omega(index, 0);
				}
			}
			
		}
		if (super_node != 0)
		{
			for (int index = 0; index < db.super_nodes[super_node - 1]->n_displacement_DOFs/3; index++)
			{
				db.super_nodes[super_node - 1]->copy_vel[index * 3 + 0] = du(0, 0);
				db.super_nodes[super_node - 1]->copy_vel[index * 3 + 1] = du(1, 0);
				db.super_nodes[super_node - 1]->copy_vel[index * 3 + 2] = du(2, 0);
			}
		}
	}
}