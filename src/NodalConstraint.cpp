#include "NodalConstraint.h"

#include "Node.h"
#include "NodeSet.h"


#include "IO.h"

NodalConstraint::NodalConstraint()
{
	number = 0;
	//setting default values for constraints (false) - in case of not reading
	UX_table.SetDefault(false);
	UY_table.SetDefault(false);
	UZ_table.SetDefault(false);
	ROTX_table.SetDefault(false);
	ROTY_table.SetDefault(false);
	ROTZ_table.SetDefault(false);
	node_set = 0;
	single_node = 0;
}

NodalConstraint::~NodalConstraint()
{
}

//Reads input file
bool NodalConstraint::Read(FILE *f)
{
	char s[1000];
	fscanf(f, "%s", s);
	number = atoi(s);
	fscanf(f, "%s", s);
	if (!strcmp(s, "NodeSet"))
	{
		fscanf(f, "%s", s);
		node_set = atoi(s);
	}
	else
		return false;
	bool DOF_OK = false;
	bool flag_continue = true;
	fpos_t pos;
	while (flag_continue == true)
	{
		DOF_OK = false;
		TryComment(f);
		fgetpos(f, &pos);
		if (fscanf(f, "%s", s) == EOF)
			return true;
		if (!strcmp(s, "UX"))
		{
			fscanf(f, "%s", s);
			if (!strcmp(s, "BoolTable"))
				UX_table.Read(f);
			else
				return false;
			DOF_OK = true;
		}
		if (!strcmp(s, "UY"))
		{
			fscanf(f, "%s", s);
			if (!strcmp(s, "BoolTable"))
				UY_table.Read(f);
			else
				return false;
			DOF_OK = true;
		}
		if (!strcmp(s, "UZ"))
		{
			fscanf(f, "%s", s);
			if (!strcmp(s, "BoolTable"))
				UZ_table.Read(f);
			else
				return false;
			DOF_OK = true;
		}
		if (!strcmp(s, "ROTX"))
		{
			fscanf(f, "%s", s);
			if (!strcmp(s, "BoolTable"))
				ROTX_table.Read(f);
			else
				return false;
			DOF_OK = true;
		}
		if (!strcmp(s, "ROTY"))
		{
			fscanf(f, "%s", s);
			if (!strcmp(s, "BoolTable"))
				ROTY_table.Read(f);
			else
				return false;
			DOF_OK = true;
		}
		if (!strcmp(s, "ROTZ"))
		{
			fscanf(f, "%s", s);
			if (!strcmp(s, "BoolTable"))
				ROTZ_table.Read(f);
			else
				return false;
			DOF_OK = true;
		}
		if (DOF_OK == true)
			flag_continue = true;
		else
		{
			fsetpos(f, &pos);
			flag_continue = false;
		}
	}
	return true;
}

//Checking inconsistencies
bool NodalConstraint::Check()
{
	if (node_set > db.number_node_sets)
		return false;
	return true;
}

//Writes output file
void NodalConstraint::Write(FILE *f)
{
	fprintf(f, "NodalConstraint\t%d\tNodeSet\t%d\n", number, node_set);
	fprintf(f, "UX\t");
	UX_table.Write(f);
	fprintf(f, "UY\t");
	UY_table.Write(f);
	fprintf(f, "UZ\t");
	UZ_table.Write(f);
	fprintf(f, "ROTX\t");
	ROTX_table.Write(f);
	fprintf(f, "ROTY\t");
	ROTY_table.Write(f);
	fprintf(f, "ROTZ\t");
	ROTZ_table.Write(f);
}
//Writes VTK XML data for post-processing
void NodalConstraint::WriteVTK_XML(FILE *f)
{
	
}
//Pre-calculus
void NodalConstraint::PreCalc()
{

}

void NodalConstraint::Mount()
{
	int node;
	//For all nodes of the node set - Computes DOFs information to nodes constraint variable
	for (int index = 0; index < db.node_sets[node_set - 1]->n_nodes; index++)
	{
		//Node number
		node = db.node_sets[node_set - 1]->node_list[index];
		/*db.nodes[node - 1]->constraints[0] = (int)UX_table.GetAt(db.current_solution_number - 1);
		db.nodes[node - 1]->constraints[1] = (int)UY_table.GetAt(db.current_solution_number - 1);
		db.nodes[node - 1]->constraints[2] = (int)UZ_table.GetAt(db.current_solution_number - 1);
		db.nodes[node - 1]->constraints[3] = (int)ROTX_table.GetAt(db.current_solution_number - 1);
		db.nodes[node - 1]->constraints[4] = (int)ROTY_table.GetAt(db.current_solution_number - 1);
		db.nodes[node - 1]->constraints[5] = (int)ROTZ_table.GetAt(db.current_solution_number - 1);*/
		if (db.nodes[node - 1]->constraints[0] + (int)UX_table.GetAt(db.current_solution_number - 1) >= 1) db.nodes[node - 1]->constraints[0] = 1;
		if (db.nodes[node - 1]->constraints[1] + (int)UY_table.GetAt(db.current_solution_number - 1) >= 1) db.nodes[node - 1]->constraints[1] = 1;
		if (db.nodes[node - 1]->constraints[2] + (int)UZ_table.GetAt(db.current_solution_number - 1) >= 1) db.nodes[node - 1]->constraints[2] = 1;
		if (db.nodes[node - 1]->constraints[3] + (int)ROTX_table.GetAt(db.current_solution_number - 1) >= 1) db.nodes[node - 1]->constraints[3] = 1;
		if (db.nodes[node - 1]->constraints[4] + (int)ROTY_table.GetAt(db.current_solution_number - 1) >= 1) db.nodes[node - 1]->constraints[4] = 1;
		if (db.nodes[node - 1]->constraints[5] + (int)ROTZ_table.GetAt(db.current_solution_number - 1) >= 1) db.nodes[node - 1]->constraints[5] = 1;
	}
}

void NodalConstraint::MountSingleNodeConstraint()
{
	if (db.nodes[single_node - 1]->constraints[0] + (int)UX_table.GetAt(db.current_solution_number - 1) >= 1) db.nodes[single_node - 1]->constraints[0] = 1;
	if (db.nodes[single_node - 1]->constraints[1] + (int)UY_table.GetAt(db.current_solution_number - 1) >= 1) db.nodes[single_node - 1]->constraints[1] = 1;
	if (db.nodes[single_node - 1]->constraints[2] + (int)UZ_table.GetAt(db.current_solution_number - 1) >= 1) db.nodes[single_node - 1]->constraints[2] = 1;
	if (db.nodes[single_node - 1]->constraints[3] + (int)ROTX_table.GetAt(db.current_solution_number - 1) >= 1) db.nodes[single_node - 1]->constraints[3] = 1;
	if (db.nodes[single_node - 1]->constraints[4] + (int)ROTY_table.GetAt(db.current_solution_number - 1) >= 1) db.nodes[single_node - 1]->constraints[4] = 1;
	if (db.nodes[single_node - 1]->constraints[5] + (int)ROTZ_table.GetAt(db.current_solution_number - 1) >= 1) db.nodes[single_node - 1]->constraints[5] = 1;

}
