#pragma once
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "Matrix.h"

class AnalyticalSurface
{
public:
	AnalyticalSurface();
	virtual ~AnalyticalSurface();

	int number;
	virtual void PreCalc() = 0;
	virtual bool Read(FILE *f) = 0;
	virtual void Write(FILE *f) = 0;
	virtual void WriteVTK_XMLRender(FILE *f) = 0;
	virtual double Distance(Matrix &O) = 0;					//Calculates the distance between the Analytical Surface and a given point in space
	virtual void N_O(Matrix *O, Matrix *Normal) = 0;		//Calculates and returns the normal direction of the closest point on analytical suface to point O in space
	virtual void T1_O(Matrix *O, Matrix *T1) = 0;			//Calculates and returns the T1 direction of the closest point on analytical suface to point O in space
	virtual void T2_O(Matrix *O, Matrix *T2) = 0;			//Calculates and returns the T2 direction of the closest point on analytical suface to point O in space
	
};

