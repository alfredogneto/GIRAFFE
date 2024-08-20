#pragma once
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "Matrix.h"
#include <vector>
#include "MatrixFloat.h"

using namespace std;

class TriangularFace
{
public:
	TriangularFace();
	~TriangularFace();

	//IDs of vertices of the triangle
	int CAD_ID;
	int verticesIDs[3];
	double vertice_weight[3];
	int edgesIDs[3];
	int ID;												//ID da face

	MatrixFloat* centroid;								//Coordenadas centroide (float)
	float radius;										//Radius
	double area;

	bool has_concave_edge;
	Matrix* normal;
	//Functions
	void PreCalc();										//PreCalc
	void EvaluateCentroid();							//Evaluates the triangle centroid
	void EvaluateRadius();								//Evaluates the radius enconpassing the triangle (centered in centroid)
	void Print(FILE *f);
};