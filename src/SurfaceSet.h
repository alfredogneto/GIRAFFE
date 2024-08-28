#pragma once
#include <stdio.h>

class SurfaceSet
{
public:
	SurfaceSet();
	~SurfaceSet();

	int n_surf;			//n�mero de superf�cies
	int* surf_list;		//lista de superf�cies
	int number;			//n�mero de refer�ncia

	
	bool sequence;		//true se � do tipo sequence
	bool list;			//true se � do tipo list

	//Para o caso de sequence
	int initial;
	int increment;

	bool Read(FILE *f);
	void Write(FILE *f);
};

