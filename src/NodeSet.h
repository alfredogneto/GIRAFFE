#pragma once
#include <stdio.h>

class NodeSet
{
public:
	NodeSet();
	~NodeSet();

	int n_nodes;		//number of nodes
	int* node_list;		//list of nodes
	int number;			//reference number

	
	bool sequence;		//true if and of the type sequence
	bool list;			//true if and of the type list

	//For the case of sequence
	int initial;
	int increment;

	bool Read(FILE *f);
	void Write(FILE *f);
	void WriteMonitor(FILE *f, bool first_record, double time);
};

