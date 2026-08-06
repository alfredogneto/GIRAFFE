#pragma once
#include <stdio.h>

class Matrix;
class SplineElement;

class Spline
{
public:
	Spline();
	~Spline();

	bool Read(FILE *f);
	void Write(FILE *f);	

	void PreCalc();								//Calculations realizados a single time
	void WriteVTK_XML_SplineMesh(FILE *f);		//Mesh of the spline junto to the paraview
	void WriteVTK_XML_SplineRender(FILE *f);	//Mesh renderizada of the spline junto to the paraview
	void CalculateSpline();						//Calculation of points of the spline
	void CalculateSplineTangentNormal();		//Calculation of tangent vector and normal the spline
	void SaveConfiguration();					//Updates the variables internal of the surface, for grab info of the pilot node for use subsequent with position updated
	bool Check();								//Checks inconsistencies for avoid errors of execution
	
	//Variables internal of the spline		
	bool alloc = false;			//Control of allocation of the spline for destructor
	int number;					//Reference number
	double radius;				//Radius that defines the surface external the spline (offset)
	int nodeset;				//Idenificador of the nodeset that defines the points of control of the spline	
	int size_nodeset;			//Number of nodes definidos in the nodeset
	int size_sp_nodes;			//Number of nodes used for descrever the spline only for visualização	
	int size_sp_elements;		//Number of elements (trechos) of spline

	int* nodeset_list;			//List of nodes of the nodeset	
	double* knot;				//Knot vector calculated according to the quantidade of nodes of the nodeset
	double* sp0;				//Coefficients for assemble the first iteration of the spline (p=0)
	double* sp1;				//Coefficients for assemble the second iteration of the spline (p=1)
	double* sp1_dd;				//Coefficients for assemble the second iteration of the spline (p=1) for second derivative
	double** sp2;				//Coefficients for assemble the third iteration of the spline (p=2)
	double** sp2_d;				//Coefficients for assemble the primeria derivative of the spline quadrada
	double** sp2_dd;			//Coefficients for assemble the second derivative of the spline quadrada	
	int* sp_elements_list;		//List of elements (trechos) of spline
	
	Matrix** x_sp_Ai;			//Points of the spline
	Matrix** x_sp_d;			//Points of the first derivative
	Matrix** x_sp_dd;			//Points of the second derivative
	Matrix** x_sp_tangent;		//Vector of tangentes for plotting
	Matrix** x_sp_normal;		//Vector of normal for plortagem

	Matrix** x_Ai;				//Points of control of the spline

	SplineElement** sp_element;	//List of elements of spline

};

