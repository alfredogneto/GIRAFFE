#pragma once
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

class Constraint
{
public:
	Constraint();
	virtual ~Constraint();

	virtual bool Read(FILE *f) = 0;						//Reads input file
	virtual void Write(FILE *f) = 0;					//Writes output file
	virtual void WriteVTK_XML(FILE *f) = 0;				//Writes VTK XML data for post-processing
	virtual void PreCalc() = 0;							//Pre-calculus
	virtual void Mount() = 0;							//Computes it to nodal constraints variables
	virtual bool Check() = 0;							//Checking inconsistencies

	int number;											//ID
};

