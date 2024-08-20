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

	int n_nodes;		//n�mero de n�s
	int* node_list;		//lista de n�s
	int number;			//n�mero de refer�ncia

	
	bool sequence;		//true se � do tipo sequence
	bool list;			//true se � do tipo list

	//Para o caso de sequence
	int initial;
	int increment;

	bool Read(FILE *f);
	void Write(FILE *f);
	void WriteMonitor(FILE *f, bool first_record, double time);
};

