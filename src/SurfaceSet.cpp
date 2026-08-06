#include "SurfaceSet.h"
#include <string>


#include"Database.h"
//Global variables
extern
Database db;

SurfaceSet::SurfaceSet()
{
	n_surf = 0;
	surf_list = NULL;

	sequence = false;
	list = false;
	initial = 0;
	increment = 0;
}

SurfaceSet::~SurfaceSet()
{
	if (surf_list != NULL)
	{
		delete[] surf_list;
	}
}

bool SurfaceSet::Read(FILE *f)
{
	char s[1000];
	fscanf(f, "%s", s);
	//Checks the keyword "SurfaceSet"
	if (!strcmp(s, "SurfaceSet"))
	{
		fscanf(f, "%s", s);
		number = atoi(s);
	}
	else
		return false;
	fscanf(f, "%s", s);
	//Checks the keyword "Surfaces"
	if (!strcmp(s, "Surfaces"))
	{
		fscanf(f, "%s", s);
		n_surf = atoi(s);
		//Allocation of the vector of nodes
		surf_list = new int[n_surf];
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
		for (int i = 0; i < n_surf; i++)
		{
			fscanf(f, "%s", s);//Reading of the number of the node
			surf_list[i] = atoi(s);
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
			for (int i = 0; i < n_surf; i++)
			{
				surf_list[i] = initial + i*increment;
			}
		}
		else
			return false;
	}
	//If atingiu this point, sign of reading correct of tudo: returns true
	return true;
}

void SurfaceSet::Write(FILE *f)
{
	fprintf(f, "SurfaceSet\t%d\tSurfaces\t%d\tList\t", number, n_surf);
	for (int i = 0; i < n_surf; i++)
		fprintf(f, "%d\t", surf_list[i]);
	fprintf(f, "\n");
}
