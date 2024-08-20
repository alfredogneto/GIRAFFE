#pragma once
#include "Matrix.h"
#include "MatrixFloat.h"

#include <vector>

using namespace std;

class Vertex
{
public:
	Vertex();
	~Vertex();
	Vertex(const Vertex &copied);						//Construtor de cópia
	float tol;											//Tolerance on coordinates
	int ID;
	MatrixFloat* coord_float;
	Matrix* coord_double;
	Vertex& operator=(const Vertex& copied);
	void Print(FILE* f);
	vector<int> faceIDs;
	vector<int> edgeIDs;
};

bool operator == (Vertex &v1, Vertex &v2);				//Verificação de igualdade com tol_equal

