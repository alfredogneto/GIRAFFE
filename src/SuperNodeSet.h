#pragma once
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "Matrix.h"
#include <math.h>
class SuperNodeSet
{
public:
	SuperNodeSet();
	~SuperNodeSet();

	int n_super_nodes;			//number of super nodes
	int* super_node_list;		//list of super nodes
	int number;					//reference number

	bool sequence;		//true if and of the type sequence
	bool list;			//true if and of the type list

	//For the case of sequence
	int initial;
	int increment;

	bool Read(FILE *f);
	void Write(FILE *f);
	void WriteMonitor(FILE *f, bool first_record, double time);
};

