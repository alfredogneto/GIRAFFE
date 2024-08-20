#pragma once
#include "Constraint.h"
#include "BoolTable.h"

class NodalConstraint :
	public Constraint
{
public:
	NodalConstraint();
	~NodalConstraint();

	bool Read(FILE *f);						//Reads input file
	void Write(FILE *f);					//Writes output file
	void WriteVTK_XML(FILE *f);				//Writes VTK XML data for post-processing
	void PreCalc();							//Pre-calculus
	void Mount();
	void MountSingleNodeConstraint();
	bool Check();							//Checking inconsistencies

	//Variables
	int node_set;

	//BoolTables for each kind of DOF
	BoolTable UX_table;
	BoolTable UY_table;
	BoolTable UZ_table;
	BoolTable ROTX_table;
	BoolTable ROTY_table;
	BoolTable ROTZ_table;

	//single node
	int single_node;
};

