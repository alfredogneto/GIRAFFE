#pragma once
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "Matrix.h"
#include <math.h>
class NodeSet
{
public:
	NodeSet();
	~NodeSet();

	int n_nodes;		//número de nós
	int* node_list;		//lista de nós
	int number;			//número de referência

	
	bool sequence;		//true se é do tipo sequence
	bool list;			//true se é do tipo list

	//Para o caso de sequence
	int initial;
	int increment;

	bool Read(FILE *f);
	void Write(FILE *f);
	void WriteMonitor(FILE *f, bool first_record, double time);
};

