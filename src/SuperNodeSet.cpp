#include "SuperNodeSet.h"

#include"Database.h"
//Global variables
extern
Database db;

SuperNodeSet::SuperNodeSet()
{
	n_super_nodes = 0;
	super_node_list = NULL;

	sequence = false;
	list = false;
	initial = 0;
	increment = 0;
}

SuperNodeSet::~SuperNodeSet()
{
	if (super_node_list != NULL)
	{
		delete[] super_node_list;
	}
}

bool SuperNodeSet::Read(FILE *f)
{
	char s[1000];
	fscanf(f, "%s", s);
	//Checks the keyword "Set"
	if (!strcmp(s, "SuperNodeSet"))
	{
		fscanf(f, "%s", s);
		number = atoi(s);
	}
	else
		return false;
	fscanf(f, "%s", s);
	//Checks the keyword "SuperNodes"
	if (!strcmp(s, "SuperNodes"))
	{
		fscanf(f, "%s", s);
		n_super_nodes = atoi(s);
		//Allocation of the vector of nodes
		super_node_list = new int[n_super_nodes];
	}
	else
		return false;
	//Two possibilities of reading:
	//1 - List
	//2 - Sequence
	fscanf(f, "%s", s);
	if (!strcmp(s, "List"))
	{
		list = true;
		for (int i = 0; i < n_super_nodes; i++)
		{
			fscanf(f, "%s", s);//Reading of the number of the node
			super_node_list[i] = atoi(s);
		}
	}
	else
	{
		if (!strcmp(s, "Sequence"))
		{
			sequence = true;
			fscanf(f, "%s", s);
			if (!strcmp(s, "Initial"))
			{
				fscanf(f, "%s", s);
				initial = atoi(s);
			}
			else
				return false;
			fscanf(f, "%s", s);
			if (!strcmp(s, "Increment"))
			{
				fscanf(f, "%s", s);
				increment = atoi(s);
			}
			else
				return false;
			//Geração of the list of nodes
			for (int i = 0; i < n_super_nodes; i++)
			{
				super_node_list[i] = initial + i * increment;
			}
		}
		else
			return false;
	}
	//If atingiu this point, sign of reading correct of tudo: returns true
	return true;
}

void SuperNodeSet::Write(FILE *f)
{
	fprintf(f, "SuperNodeSet\t%d\tNodes\t%d\tList\t", number, n_super_nodes);
	for (int i = 0; i < n_super_nodes; i++)
		fprintf(f, "%d\t", super_node_list[i]);
	fprintf(f, "\n");
}

void SuperNodeSet::WriteMonitor(FILE *f, bool first_record, double time)
{
	//TO DO
}

