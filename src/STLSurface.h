#pragma once
#include <stdio.h>
#include <vector>
#include "CADData.h"

#include "Vertex.h"
#include "Edge.h"
#include "Tetrahedron.h"

class TriangularFace;
class MatrixFloat;
using namespace std;

class STLSurface :
	public CADData
{
public:
	STLSurface();
	~STLSurface();

	vector<float> coord_float;										//Coordinates of the points of the file stl (CAD)
	vector<double> coord_double;									//Coordinates of the points of the file stl (CAD)
	int n_CAD_points;												//Number of points of the file stl (CAD)

	bool Read(FILE *f);												//Reading of the input file
	void Write(FILE *f);											//File writing of output
	void PreCalc();													//PreCalc

	bool ReadCADFile();												//File reading of CAD
	bool ReadMeshFile();											//File reading of mesh (optional)
	void WriteVTK_XMLRender(FILE *f, const Matrix& pos, const Matrix& rot, int number);		//Plots CAD Data
	void WriteVTK_XMLMesh(MatrixFloat& pos, MatrixFloat& rot);		//Plots CAD Data

	void CreateVerticesEdges();
	void MergeVertices();
	void MergeEdges();
	void OrganizeNumbering();
	void SetNumberingInfo();
	void ReplaceVertexID(int index, int newID);
	void ReplaceEdgeID(int index, int newID);
	void MarkConcaveEdges();
	void PointNormalEdges();
	void GenerateTetraMesh();
	void EvaluateVerticeFactors();

	void PrintSurfaceReport();

	int GetVertexAssociatedwithBothEdges(int edge1, int edge2);

	//Geometric evaluation functions
	void EvaluateVolume();
	void EvaluateCentroid();
	void EvaluateInertiaTensor();
	void EvaluateRadius();

	//Specific variables - stl
	int n_faces;				//number of faces
	TriangularFace** faces;		//faces vector
	vector<Vertex> vertices;
	vector<Edge> edges;
	vector<Tetrahedron> tetras;

	double total_ref_area;			//particle reference surface area
	double* vertice_factors;		//factors of total mass for each vertex to compound mass matrix

	bool mesh_available;
};