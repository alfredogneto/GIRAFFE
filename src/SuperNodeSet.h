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

	int n_super_nodes;			//n�mero de super nodes
	int* super_node_list;		//lista de super nodes
	int number;					//n�mero de refer�ncia

	bool sequence;		//true se � do tipo sequence
	bool list;			//true se � do tipo list

	//Para o caso de sequence
	int initial;
	int increment;

	bool Read(FILE *f);
	void Write(FILE *f);
	void WriteMonitor(FILE *f, bool first_record, double time);
};

