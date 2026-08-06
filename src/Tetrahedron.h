#pragma once
#include <stdio.h>

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

