#pragma once
#include <stdio.h>

class ElementSet
{
public:
	ElementSet();
	~ElementSet();

	int n_el;			//n�mero de elementos
	int* el_list;		//lista de superf�cies
	int number;			//n�mero de refer�ncia

	
	bool sequence;		//true se � do tipo sequence
	bool list;			//true se � do tipo list

	//Para o caso de sequence
	int initial;
	int increment;

	bool Read(FILE *f);
	void Write(FILE *f);
};

