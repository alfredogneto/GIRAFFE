#pragma once
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "Matrix.h"
#include <vector>
#include "MatrixFloat.h"

using namespace std;
class Tetrahedron
{
public:
	Tetrahedron();
	~Tetrahedron();
	void Print(FILE *f);
	//IDs of vertices of the triangle
	int CAD_ID;
	int verticesIDs[4];
	int edgesIDs[6];
	int ID;												//ID do tetraedro
};

