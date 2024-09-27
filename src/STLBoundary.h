#pragma once
#include "Boundary.h"
class MatrixFloat;

class STLBoundary :
	public Boundary
{
public:
	STLBoundary();
	~STLBoundary();

	int CADDATA_ID;

	bool Read(FILE *f);
	void Write(FILE *f);
	void WriteModifyingParameters(FILE *f, int e_material, int e_node, int e_number, int e_cs);
	bool Check();

	void PreCalc();
	void UpdateBoundingVolumes();
	void SaveLagrange();
	void WriteVTK_XMLBase(FILE *f);
	void WriteVTK_XMLRender(FILE *f);

	float bv_factor;		//Controls the size of bounding volumes of edges and vertices
	float inc_len_factor;	//Controls inflation of bounding volumes

	MatrixFloat* x0f;
	MatrixFloat* Q0f;

	//Explicit
	void InitialEvaluations();
};

