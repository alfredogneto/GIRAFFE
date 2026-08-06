#pragma once
#include <stdio.h>
#include <vector>

class Matrix;

using namespace std;

class Edge
{
public:
	Edge();
	~Edge();
	Edge(const Edge &copied);								//Copy constructor
	Edge& operator=(const Edge& copied);
	void Print(FILE *f);
	int ID;
	int verticesIDs[2];

	int concave_indicator;

	int CAD_ID;					//associated CAD ID
	void PreCalc();
	void EvaluateLength();
	float length;

	vector<int> faceIDs;
	//Pointers to face normals
	Matrix* n1;
	Matrix* n2;
};

bool operator == (const Edge &e1, const Edge &e2);				//Equality check with tol_equal

