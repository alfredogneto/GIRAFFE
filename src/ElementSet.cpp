#include "ElementSet.h"

#include"Database.h"
//Global variables
extern
Database db;

ElementSet::ElementSet()
{
	n_el = 0;
	el_list = NULL;

	sequence = false;
	list = false;
	initial = 0;
	increment = 0;
}

ElementSet::~ElementSet()
{
	if (el_list != NULL)
	{
		delete[] el_list;
	}
}

bool ElementSet::Read(FILE *f)
{
	char s[1000];
	fscanf(f, "%s", s);
	//Checks the keyword "ElementSet"
	if (!strcmp(s, "ElementSet"))
	{
		fscanf(f, "%s", s);
		number = atoi(s);
	}
	else
		return false;
	fscanf(f, "%s", s);
	//Checks the keyword "Elements"
	if (!strcmp(s, "Elements"))
	{
		fscanf(f, "%s", s);
		n_el = atoi(s);
		//Allocation of the vector of nodes
		el_list = new int[n_el];
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
		for (int i = 0; i < n_el; i++)
		{
			fscanf(f, "%s", s);//Reading of the number of the element
			el_list[i] = atoi(s);
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
			//Geração of the list of elements
			for (int i = 0; i < n_el; i++)
			{
				el_list[i] = initial + i*increment;
			}
		}
		else
			return false;
	}
	//If atingiu this point, sign of reading correct of tudo: returns true
	return true;
}

void ElementSet::Write(FILE *f)
{
	fprintf(f, "ElementSet\t%d\tSurfaces\t%d\tList\t", number, n_el);
	for (int i = 0; i < n_el; i++)
		fprintf(f, "%d\t", el_list[i]);
	fprintf(f, "\n");
}
