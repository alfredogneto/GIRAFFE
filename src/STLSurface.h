#pragma once
#include "CADData.h"
#include "TriangularFace.h"
#include "Vertex.h"
#include "Edge.h"
#include "Tetrahedron.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

using namespace std;

class STLSurface :
	public CADData
{
public:
	STLSurface();
	~STLSurface();

	vector<float> coord_float;										//Coordenadas dos pontos do arquivo stl (CAD)
	vector<double> coord_double;									//Coordenadas dos pontos do arquivo stl (CAD)
	int n_CAD_points;												//Número de pontos do arquivo stl (CAD)

	bool Read(FILE *f);												//Leitura do arquivo de entrada
	void Write(FILE *f);											//Escrita do arquivo de saída
	void PreCalc();													//PreCalc

	bool ReadCADFile();												//Leitura do arquivo de CAD
	bool ReadMeshFile();											//Leitura do arquivo de malha (opcional)
	void WriteVTK_XMLRender(FILE *f, Matrix& pos, Matrix& rot, int number);		//Plots CAD Data
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