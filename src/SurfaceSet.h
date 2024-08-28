#pragma once
#include <stdio.h>

class SurfaceSet
{
public:
	SurfaceSet();
	~SurfaceSet();

	int n_surf;			//número de superfícies
	int* surf_list;		//lista de superfícies
	int number;			//número de referência

	
	bool sequence;		//true se é do tipo sequence
	bool list;			//true se é do tipo list

	//Para o caso de sequence
	int initial;
	int increment;

	bool Read(FILE *f);
	void Write(FILE *f);
};

