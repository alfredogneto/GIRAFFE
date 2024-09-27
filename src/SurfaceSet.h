#pragma once
#include <stdio.h>

class SurfaceSet
{
public:
	SurfaceSet();
	~SurfaceSet();

	int n_surf;			//numero de superficies
	int* surf_list;		//lista de superficies
	int number;			//numero de referência

	
	bool sequence;		//true se e do tipo sequence
	bool list;			//true se e do tipo list

	//Para o caso de sequence
	int initial;
	int increment;

	bool Read(FILE *f);
	void Write(FILE *f);
};

