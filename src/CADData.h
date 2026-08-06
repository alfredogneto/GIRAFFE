#pragma once
#include <stdio.h>

class Matrix;

class CADData
{
public:
	CADData() {}
	virtual ~CADData() {}
	int number;
	char file[200];																//File name CAD
	virtual bool Read(FILE *f) = 0;												//Reading of the input file
	virtual void Write(FILE *f) = 0;											//File writing of output
	virtual void PreCalc() = 0;													//PreCalc
	virtual void WriteVTK_XMLRender(FILE *f, const Matrix& pos, const Matrix& rot, int number) = 0;		//Plots CAD Data

	//Geometric evaluation functions
	virtual void EvaluateVolume() = 0;
	virtual void EvaluateCentroid() = 0;
	virtual void EvaluateInertiaTensor() = 0;
	virtual void EvaluateRadius() = 0;

	Matrix* J_O;																//Tensor of inertia in relação the origin of the CS local (CAD)
	Matrix* G;																	//Position of the baricentro in relação the origin of the CS local (CAD)
	double volume;																//volume

	float radius;																//radius encompassing the whole surface (w/r to the origin)
};	