#pragma once
#include <stdio.h>

class ElementSet
{
public:
	ElementSet();
	~ElementSet();

	int n_el;			//number of elements
	int* el_list;		//list of surfaces
	int number;			//reference number

	
	bool sequence;		//true if and of the type sequence
	bool list;			//true if and of the type list

	//For the case of sequence
	int initial;
	int increment;

	bool Read(FILE *f);
	void Write(FILE *f);
};

