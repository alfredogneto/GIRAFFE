#pragma once
#include "Displacement.h"
#include "BoolTable.h"
#include "Matrix.h"

class Table;
class MathCode;

//using namespace exprtk;

class NodalDisplacement :
	public Displacement
{
public:
	NodalDisplacement();
	~NodalDisplacement();
	bool Read(FILE *f);						//Reads input file
	void Write(FILE *f);					//Writes output file
	void WriteVTK_XML(FILE *f);				//Writes VTK XML data for post-processing
	void PreCalc();							//Pre-calculus
	void Mount();
	bool Check();							//Checking inconsistencies
	double GetValueAt(double t, int position);
	//Variables
	int node_set;
	int cs;
	char angular_parameters[100];
	Matrix Q;
	Matrix I;
	void MountEulerVector(Matrix* euler, Matrix* rodrigues, int temp_node, double time);

	void EvaluateExplicit(double t);

	void MountSingleNodeDisplacement();
	
	MathCode* mcode;									//To work with a mathematical expression to generate a code
	Table* table;										//Table with displacement data
	int n_times;										//NTimes (of table)
	int n_values;										//Number of values (columns) on table

	//Variables - definition to a single node and global CS
	int single_node;

	BoolTable bool_table;				//Bool table que ativa ou desativa na seq. de solutions
};

