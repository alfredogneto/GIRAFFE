#pragma once
#include "AnalyticalSurface.h"

class Plane :
	public AnalyticalSurface
{
public:
	Plane();
	~Plane();

	Matrix* N;	//normal
	Matrix* T1;
	Matrix* T2;
	Matrix* P;	//Point pertencente to the plane
	void PreCalc();
	double len;					//Size of the bounding box (for pós-processamento)
	bool Read(FILE *f);
	void Write(FILE *f);
	double Distance(Matrix &O);	//Calculates the distance between the Analytical Surface and a given point in space
	void N_O(Matrix *O, Matrix *Normal);		//Calculates and returns the normal direction of the closest point on analytical suface to point O in space
	void T1_O(Matrix *O, Matrix *Tang1);			//Calculates and returns the T1 direction of the closest point on analytical suface to point O in space
	void T2_O(Matrix *O, Matrix *Tang2);			//Calculates and returns the T2 direction of the closest point on analytical suface to point O in space
	void WriteVTK_XMLRender(FILE *f);
	
};

