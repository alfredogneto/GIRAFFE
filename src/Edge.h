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
	Edge(const Edge &copied);								//Construtor de cópia
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

bool operator == (Edge &e1, Edge &e2);				//Verificação de igualdade com tol_equal

