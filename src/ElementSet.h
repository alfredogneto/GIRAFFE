#pragma once
#include <stdio.h>

class ElementSet
{
public:
	ElementSet();
	~ElementSet();

	int n_el;			//numero de elementos
	int* el_list;		//lista de superficies
	int number;			//numero de referência

	
	bool sequence;		//true se e do tipo sequence
	bool list;			//true se e do tipo list

	//Para o caso de sequence
	int initial;
	int increment;

	bool Read(FILE *f);
	void Write(FILE *f);
};

