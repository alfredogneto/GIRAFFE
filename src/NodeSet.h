#pragma once
#include <stdio.h>

class NodeSet
{
public:
	NodeSet();
	~NodeSet();

	int n_nodes;		//numero de nós
	int* node_list;		//lista de nós
	int number;			//numero de referência

	
	bool sequence;		//true se e do tipo sequence
	bool list;			//true se e do tipo list

	//Para o caso de sequence
	int initial;
	int increment;

	bool Read(FILE *f);
	void Write(FILE *f);
	void WriteMonitor(FILE *f, bool first_record, double time);
};

