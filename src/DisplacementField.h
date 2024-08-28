#pragma once
#include "Displacement.h"
#include "Matrix.h"


class DisplacementField :
	public Displacement
{
public:
	DisplacementField();
	~DisplacementField();

	bool Read(FILE *f);						//Reads input file
	void Write(FILE *f);					//Writes output file
	void WriteVTK_XML(FILE *f);				//Writes VTK XML data for post-processing
	void PreCalc();							//Pre-calculus
	void Mount();
	bool Check();							//Checking inconsistencies
	void Alloc();							//Alocação dinamica de variaveis
	void Free();							//Desalocação dinamica de variaveis

	//Variables
	bool alloced;							//Controle de alocação de memória
	int solution_step;						//Solution step to prescribe displacement field
	int cs;									//Coordinate system for input data
	int *nodes;								//List of nodes to apply the displacement field
	double** displacements;					//List of displacements to be applied
	int n_nodes;							//NNodes
	int n_values;							//Number of values (columns) on table

	void EvaluateExplicit(double t);

	char angular_parameters[100];
	Matrix Q;
	Matrix I;
	void MountEulerVector(Matrix* euler, Matrix* rodrigues, int temp_node, int index_node, double factor);
};

