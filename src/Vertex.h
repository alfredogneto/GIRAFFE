#pragma once
#include <vector>

class Matrix;
class MatrixFloat;

using namespace std;

class Vertex
{
public:
	Vertex();
	~Vertex();
	Vertex(const Vertex &copied);						//Copy constructor
	float tol;											//Tolerance on coordinates
	int ID;
	MatrixFloat* coord_float;
	Matrix* coord_double;
	Vertex& operator=(const Vertex& copied);
	void Print(FILE* f);
	vector<int> faceIDs;
	vector<int> edgeIDs;
};

bool operator == (const Vertex &v1, const Vertex &v2);				//Equality check with tol_equal

