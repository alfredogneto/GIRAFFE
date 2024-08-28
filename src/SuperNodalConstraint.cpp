#include "SuperNodalConstraint.h"

#include "SuperNode.h"
#include "SuperNodeSet.h"
#include "IO.h"

SuperNodalConstraint::SuperNodalConstraint()
{
	number = 0;
	//setting default values for booltable (false) - in case of not reading
	bool_table.SetDefault(false);
	super_node = 0;
	super_node_set = 0;
	all = true;
	DOF_table.clear();
}


SuperNodalConstraint::~SuperNodalConstraint()
{
}

//Reads input file
bool SuperNodalConstraint::Read(FILE *f)
{
	char s[1000];
	fscanf(f, "%s", s);
	number = atoi(s);
	fscanf(f, "%s", s);
	if (!strcmp(s, "SuperNode"))
	{
		fscanf(f, "%s", s);
		super_node = atoi(s);
	}
	else
	{
		if (!strcmp(s, "SuperNodeSet"))
		{
			fscanf(f, "%s", s);
			super_node_set = atoi(s);
		}
		else
		{
			return false;
		}
	}

	
	TryComment(f);
	//Reading - if EOF, return error on reading (expected data not found)
	if (fscanf(f, "%s", s) == EOF)
		return false;
	//Case ALL 
	if (!strcmp(s, "ALL"))
	{
		all = true;
		fscanf(f, "%s", s);
		if (!strcmp(s, "BoolTable"))
			bool_table.Read(f);
		else
			return false;
	}
	else
	{
		if (!strcmp(s, "SELECTED"))
		{
			all = false;
			//Reads and returns an error if EOF is detected (data with DOFs or a BoolTable were expected)
			if (fscanf(f, "%s", s) == EOF)
				return false;
			while (strcmp(s, "BoolTable"))
			{
				//Attributes a number to the list
				DOF_table.push_back((atoi(s)));
				//Reads the next and returns an error if EOF is detected (data with DOFs or a BoolTable were expected)
				if (fscanf(f, "%s", s) == EOF)
					return false;
			}
			bool_table.Read(f);
		}
		else
			return false;
	}
	return true;
}

//Checking inconsistencies
bool SuperNodalConstraint::Check()
{
	if (super_node != 0)
		if (super_node > db.number_super_nodes)
			return false;
	if (super_node_set != 0)
		if (super_node_set > db.number_super_node_sets)
			return false;
	return true;
}

//Writes output file
void SuperNodalConstraint::Write(FILE *f)
{
	fprintf(f, "SuperNodalConstraint\t%d\t", number);
	if (super_node != 0)
		fprintf(f, "SuperNode\t%d\n", super_node);
	if (super_node_set != 0)
		fprintf(f, "SuperNodeSet\t%d\n", super_node_set);

	if (all == true)
	{
		fprintf(f, "ALL\t");
		bool_table.Write(f);
	}
	else
	{
		fprintf(f, "SELECTED\t");
		//Print DOFs vector
		for (int i = 0; i < (int)DOF_table.size(); i++)
			fprintf(f, "%d\t", DOF_table[i]);
		bool_table.Write(f);
	}
}
//Writes VTK XML data for post-processing
void SuperNodalConstraint::WriteVTK_XML(FILE *f)
{

}
//Pre-calculus
void SuperNodalConstraint::PreCalc()
{

}

void SuperNodalConstraint::Mount()
{
	//If in the current load step the constraint is active
	if (bool_table.GetAt(db.current_solution_number - 1))
	{
		if (super_node != 0)
		{
			if (all == true)
			{
				//Rolls over all DOFs and set them as constrained (true)
				for (int i = 0; i < db.super_nodes[super_node - 1]->n_DOFs; i++)
				{
					db.super_nodes[super_node - 1]->constraints[i] = true;
				}
			}
			else
			{
				//Rolls over the table of constrained DOFs
				for (int i = 0; i < (int)DOF_table.size(); i++)
				{
					//If the DOF is within the range of DOFs of the super node, contrain it
					if (DOF_table[i] <= db.super_nodes[super_node - 1]->n_DOFs)
						db.super_nodes[super_node - 1]->constraints[DOF_table[i] - 1] = true;
				}
			}
		}
		if (super_node_set != 0)
		{
			for (int sn = 0; sn < db.super_node_sets[super_node_set - 1]->n_super_nodes; sn++)
			{
				int temp_super_node = db.super_node_sets[super_node_set - 1]->super_node_list[sn];

				if (all == true)
				{
					//Rolls over all DOFs and set them as constrained (true)
					for (int i = 0; i < db.super_nodes[temp_super_node - 1]->n_DOFs; i++)
					{
						db.super_nodes[temp_super_node - 1]->constraints[i] = true;
					}
				}
				else
				{
					//Rolls over the table of constrained DOFs
					for (int i = 0; i < (int)DOF_table.size(); i++)
					{
						//If the DOF is within the range of DOFs of the super node, contrain it
						if (DOF_table[i] <= db.super_nodes[temp_super_node - 1]->n_DOFs)
							db.super_nodes[temp_super_node - 1]->constraints[DOF_table[i] - 1] = true;
					}
				}
			}
		}
	}
}
