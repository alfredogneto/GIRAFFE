#pragma once
#include <stdio.h>
#include <vector>
#include "CADData.h"
#include <Eigen/dense>

class NURBSSurface;
class MatrixFloat;

class NURBSMultipatchSurface :
	public CADData
{
public:
	NURBSMultipatchSurface();
	~NURBSMultipatchSurface();

	//Dados NURBS
	std::vector<NURBSSurface*> patches; // vector with the NURBS patches that define the whole surface
	int n_patches; // number of patches

	bool Read(FILE *f);												//Leitura do arquivo de entrada
	void Write(FILE *f);											//Escrita do arquivo de saida
	void PreCalc();													//PreCalc

	//Geometric evaluation functions
	void EvaluateVolume();
	void EvaluateCentroid();
	void EvaluateInertiaTensor();
	void EvaluateRadius();
	//Marina
	bool ReadNurbsData();


	bool ReadCADFile();												//Leitura do arquivo de CAD
	void WriteVTK_XMLRender(FILE *f, Matrix& pos, Matrix& rot, int number);
};
