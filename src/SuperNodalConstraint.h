#pragma once
#include "BoolTable.h"
#include "Constraint.h"

using namespace std;

class SuperNodalConstraint :
	public Constraint
{
public:
	SuperNodalConstraint();
	~SuperNodalConstraint();

	bool Read(FILE *f);						//Reads input file
	void Write(FILE *f);					//Writes output file
	void WriteVTK_XML(FILE *f);				//Writes VTK XML data for post-processing
	void PreCalc();							//Pre-calculus
	void Mount();
	bool Check();							//Checking inconsistencies

	//Variables
	int super_node;							//super node ID
	int super_node_set;						//super node set
	bool all;								//indicates that ALL DOFs has to be constrained
	
	vector<int> DOF_table;					//In case of selected DOFs contains the ones chosen by the user to be constrained
	BoolTable bool_table;					//Controls the activation of not of the constraints in solution steps
};

