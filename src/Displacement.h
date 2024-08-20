#pragma once
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "Table.h"
#include "Matrix.h"
#include "MathCode.h"

class Displacement
{
public:
	Displacement();
	virtual ~Displacement();

	virtual bool Read(FILE *f) = 0;						//Reads input file
	virtual void Write(FILE *f) = 0;					//Writes output file
	virtual void WriteVTK_XML(FILE *f) = 0;				//Writes VTK XML data for post-processing
	virtual void PreCalc() = 0;							//Pre-calculus
	virtual void Mount() = 0;							//Evaluates the displacement and incluces in global vector
	virtual bool Check() = 0;							//Checking inconsistencies

	virtual void EvaluateExplicit(double t) = 0;
	
	int number;											//ID
	
};

