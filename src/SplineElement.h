#pragma once
#include <stdio.h>

#include "Box.h"

class Matrix;

class SplineElement
{
public:
	SplineElement();
	~SplineElement();

	void PreCalc();

	bool ReadCommon(FILE *f);

	//Functions
	void FillNodes();														//Updates the variables internal of the surface, for grab info of the pilot node for use subsequent with position updated
	void SaveConfiguration();												//Saves vectors of configuration converged
	void CenterPoint(Matrix* center, double* radii);						//Returns coordinates global of the point central of the surface the be used for calculations coarse of its location (pinball)
	bool Check();															//Checks inconsistencies for avoid errors of execution
	void SetMinMaxRange();
	void SplinePoint(double & zeta, Matrix& point);					//Obtains point of the spline
	void UpdateBox();

	int *nodes;		//Nodes global - conectividade
	int n_nodes;	//Number of nodes of the surface
	int **DOFs;		//Indicates for the indexing of each degree of freedom, 1 or 0, active or inactive for the element in question
	int nDOFs;		//Number of GL
	//int VTK_type;	//Type of celula for VTK
	int *VTK_nodes;	//Indexing for convert numbering of the formato giraffe for the formato of the celula equivalente of the paraview
	int **GLs;		//Pointer for the GL global used in the surface

	Box box;

	//Bool indicators of coordinates
	bool entered_u1;
	double u1_min, u1_max, u1_range;
	//Divisions on degeneration
	int div1, div2;

	double* knot_element;

	double* radius;

	Matrix* x_Ai;
	Matrix* x_Bi;
	Matrix* x_Ci;

	Matrix* x_Ap;
	Matrix* x_Bp;
	Matrix* x_Cp;

	Matrix* d;
	Matrix* dui;
	Matrix* ddui;

	Matrix* I3;																//Identidade of order 3

};

