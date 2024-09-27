#include "Edge.h"

#include "STLSurface.h"
#include "Matrix.h"
#include "MatrixFloat.h"
#include"Database.h"
//Variaveis globais
extern
Database db;

Edge::Edge()
{
	CAD_ID = 0;
	ID = 0;
	verticesIDs[0] = 0;
	verticesIDs[1] = 0;
	concave_indicator = 0;
	faceIDs.clear();
	length = 0.0f;
}

//Construtor de cópia
Edge::Edge(const Edge &copied)
{
	ID = copied.ID;
	verticesIDs[0] = copied.verticesIDs[0];
	verticesIDs[1] = copied.verticesIDs[1];

	concave_indicator = copied.concave_indicator;

	faceIDs.clear();
	for (int i = 0; i < copied.faceIDs.size(); i++)
		faceIDs.push_back(copied.faceIDs[i]);
}

Edge& Edge::operator=(const Edge& copied)
{
	ID = copied.ID;
	verticesIDs[0] = copied.verticesIDs[0];
	verticesIDs[1] = copied.verticesIDs[1];

	concave_indicator = copied.concave_indicator;

	faceIDs.clear();
	for (int i = 0; i < copied.faceIDs.size(); i++)
		faceIDs.push_back(copied.faceIDs[i]);

	return *this;
}

Edge::~Edge()
{
	faceIDs.clear();
}

//Verificação de igualdade
bool operator == (Edge &e1, Edge &e2)
{
	if (e1.verticesIDs[0] == e2.verticesIDs[0] && e1.verticesIDs[1] == e2.verticesIDs[1])
		return true;
	if (e1.verticesIDs[0] == e2.verticesIDs[1] && e1.verticesIDs[1] == e2.verticesIDs[0])
		return true;
	return false;
}
void Edge::Print(FILE* f)
{
	if (concave_indicator > 0)
		fprintf(f,"Edge %d\t%d\t%d\tConcave\n", ID - 1, verticesIDs[0] - 1,verticesIDs[1] - 1);
	else
	{
		if (concave_indicator == 0)
			fprintf(f, "Edge %d\t%d\t%d\tFlat\n", ID - 1, verticesIDs[0] - 1, verticesIDs[1] - 1);
		else
			fprintf(f, "Edge %d\t%d\t%d\tConvex\n", ID - 1, verticesIDs[0] - 1, verticesIDs[1] - 1);
	}
		
}

void Edge::PreCalc()
{
	EvaluateLength();
}
void Edge::EvaluateLength()
{
	STLSurface* ptr = static_cast<STLSurface*>(db.cad_data[CAD_ID - 1]);
	//Evaluates the lenght of the edge
	length = sqrt(dot(	*ptr->vertices[verticesIDs[1] - 1].coord_float - *ptr->vertices[verticesIDs[0] - 1].coord_float, 
						*ptr->vertices[verticesIDs[1] - 1].coord_float - *ptr->vertices[verticesIDs[0] - 1].coord_float));
	
}
