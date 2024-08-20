#include "Vertex.h"

Vertex::Vertex()
{
	coord_float = new MatrixFloat(3);
	coord_double = new Matrix(3);
	tol = (float)1e-6;
	faceIDs.clear();
	edgeIDs.clear();
}

//Construtor de cópia
Vertex::Vertex(const Vertex &copied)
{
	coord_float = new MatrixFloat(3);
	coord_double = new Matrix(3);
	*coord_float = *copied.coord_float;
	*coord_double = *copied.coord_double;
	tol = copied.tol;
	ID = copied.ID;
	
	faceIDs.clear();
	for (int i = 0; i < copied.faceIDs.size(); i++)
		faceIDs.push_back(copied.faceIDs[i]);

	edgeIDs.clear();
	for (int i = 0; i < copied.edgeIDs.size(); i++)
		edgeIDs.push_back(copied.edgeIDs[i]);
}

Vertex& Vertex::operator=(const Vertex& copied) 
{
	*coord_float = *copied.coord_float;
	*coord_double = *copied.coord_double;
	tol = copied.tol;
	ID = copied.ID;

	faceIDs.clear();
	for (int i = 0; i < copied.faceIDs.size(); i++)
		faceIDs.push_back(copied.faceIDs[i]);

	return *this;
}

Vertex::~Vertex()
{
	delete coord_float;
	delete coord_double;
	faceIDs.clear();
	edgeIDs.clear();
}

//Verificação de igualdade
bool operator == (Vertex &v1, Vertex &v2)
{
	float tol_equal;
	if (v1.tol > v2.tol)
		tol_equal = v1.tol;
	else
		tol_equal = v2.tol;

	float dist = sqrt(dot(*v1.coord_float - *v2.coord_float, *v1.coord_float - *v2.coord_float));
	if (dist <= tol_equal)
		return true;
	else
		return false;
}
void Vertex::Print(FILE* f)
{
	fprintf(f,"Vertex %d\t%.6f\t%.6f\t%.6f\n", ID - 1, (*coord_float)(0, 0), (*coord_float)(1, 0), (*coord_float)(2, 0));
}
