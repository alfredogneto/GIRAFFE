#include "TriangularFace.h"
#include "STLSurface.h"
#include"Database.h"
//Variáveis globais
extern
Database db;

TriangularFace::TriangularFace()
{
	radius = 0;
	centroid = new MatrixFloat(3);
	ID = 0;
	has_concave_edge = false;
	normal = new Matrix(3);
	vertice_weight[0] = 0.0;
	vertice_weight[1] = 0.0;
	vertice_weight[2] = 0.0;
}

//PreCalc
void TriangularFace::PreCalc()
{
	EvaluateCentroid();
	EvaluateRadius();

	//Normal
	STLSurface* ptr = static_cast<STLSurface*>(db.cad_data[CAD_ID - 1]);
	Matrix e1 = *ptr->vertices[verticesIDs[1] - 1].coord_double - *ptr->vertices[verticesIDs[0] - 1].coord_double;
	Matrix e2 = *ptr->vertices[verticesIDs[2] - 1].coord_double - *ptr->vertices[verticesIDs[0] - 1].coord_double;
	*normal = (1.0 / norm(cross(e1, e2)))*cross(e1, e2);

	//Area
	area = 0.5 * norm(cross(e1, e2));

	double temp;
	//Ponderações angulares dos vértices (de acordo com o ângulo interno do triângulo)
	//vertice [0]
	e1 = *ptr->vertices[verticesIDs[1] - 1].coord_double - *ptr->vertices[verticesIDs[0] - 1].coord_double;
	e2 = *ptr->vertices[verticesIDs[2] - 1].coord_double - *ptr->vertices[verticesIDs[0] - 1].coord_double;
	temp = asin(norm(cross(e1, e2)) / (norm(e1)*norm(e2)));
	if (dot(e1, e2) < 0.0)
		temp = PI - temp;
	vertice_weight[0] = temp / PI;
	//vertice [1]
	e1 = *ptr->vertices[verticesIDs[0] - 1].coord_double - *ptr->vertices[verticesIDs[1] - 1].coord_double;
	e2 = *ptr->vertices[verticesIDs[2] - 1].coord_double - *ptr->vertices[verticesIDs[1] - 1].coord_double;
	temp = asin(norm(cross(e1, e2)) / (norm(e1)*norm(e2)));
	if (dot(e1, e2) < 0.0)
		temp = PI - temp;
	vertice_weight[1] = temp / PI;
	//vertice [2]
	e1 = *ptr->vertices[verticesIDs[0] - 1].coord_double - *ptr->vertices[verticesIDs[2] - 1].coord_double;
	e2 = *ptr->vertices[verticesIDs[1] - 1].coord_double - *ptr->vertices[verticesIDs[2] - 1].coord_double;
	temp = asin(norm(cross(e1, e2)) / (norm(e1)*norm(e2)));
	if (dot(e1, e2) < 0.0)
		temp = PI - temp;
	vertice_weight[2] = temp / PI;
}


TriangularFace::~TriangularFace()
{
	delete centroid;
	delete normal;
}

//Evaluates the triangle centroid
void TriangularFace::EvaluateCentroid()
{
	
	STLSurface* ptr = static_cast<STLSurface*>(db.cad_data[CAD_ID - 1]);
	*centroid = (float)0.3333333333333333333333333333*(*ptr->vertices[verticesIDs[0] - 1].coord_float + *ptr->vertices[verticesIDs[1] - 1].coord_float + *ptr->vertices[verticesIDs[2] - 1].coord_float);
}

//Evaluates the radius enconpassing the triangle (centered in centroid)
void TriangularFace::EvaluateRadius()
{
	STLSurface* ptr = static_cast<STLSurface*>(db.cad_data[CAD_ID - 1]);
	//Evaluates the triangular face radius (w/r to the centroid)
	radius = 0.0;
	float temp_len = 0.0;
	temp_len = sqrt(dot(*ptr->vertices[verticesIDs[0] - 1].coord_float - *centroid, *ptr->vertices[verticesIDs[0] - 1].coord_float - *centroid));
	if (temp_len > radius)
		radius = temp_len;
	temp_len = sqrt(dot(*ptr->vertices[verticesIDs[1] - 1].coord_float - *centroid, *ptr->vertices[verticesIDs[1] - 1].coord_float - *centroid));
	if (temp_len > radius)
		radius = temp_len;
	temp_len = sqrt(dot(*ptr->vertices[verticesIDs[2] - 1].coord_float - *centroid, *ptr->vertices[verticesIDs[2] - 1].coord_float - *centroid));
	if (temp_len > radius)
		radius = temp_len;
}

void TriangularFace::Print(FILE *f)
{
	fprintf(f,"Face %d\t%d\t%d\t%d\n", ID - 1, verticesIDs[0] - 1, verticesIDs[1] - 1, verticesIDs[2] - 1);
}
