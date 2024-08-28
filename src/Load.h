#pragma once
#include <stdio.h>

class Table;
class MathCode;

class Load
{
public:
	Load();
	virtual ~Load();
	int number;											//ID
	virtual bool Read(FILE *f) = 0;						//Reads input file
	virtual void Write(FILE *f) = 0;					//Writes output file
	virtual void WriteVTK_XML(FILE *f) = 0;				//Writes VTK XML data for post-processing
	virtual void PreCalc() = 0;							//Pre-calculus - chamado no inicio - somente uma vez
	virtual void UpdateforSolutionStep() = 0;			//Atualiza dados necessarios e que sejam dependentes de DOFs ativos/inativos - chamado no inicio de cada solution step
	virtual void Mount() = 0;							//Evaluates the load and its tangent operator (when applicable). 
														//Additionally, when applicable inserts contributions to the global vector/matrix
	virtual bool Check() = 0;							//Checking inconsistencies

	virtual void EvaluateExplicit(double t) = 0;

	double GetValueAt(double t, int position);

protected:
	MathCode* mcode;									//To work with a mathematical expression to generate a code
	Table* table;										//Table with load data
	int n_times;										//NTimes (of table)
	int n_values;										//Number of values (columns) on table
};

