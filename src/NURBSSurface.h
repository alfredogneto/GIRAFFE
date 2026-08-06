#pragma once
#include "CADData.h"

class NURBSSurface :
	public CADData
{
public:
	NURBSSurface();
	~NURBSSurface();

	
	//Data NURBS
	int U_dim;
	int U_order;
	int V_dim;
	int V_order;
	double* U_knot_vector;
	double* V_knot_vector;
	double** weights;
	Matrix** control_points;
	Matrix** Pw;	//control points in 4D (already with pesos)
	double** Bin;	//Binomial coefficients

	bool Read(FILE *f);												//Reading of the input file
	void Write(FILE *f);											//File writing of output
	void PreCalc();													//PreCalc

	//Geometric evaluation functions
	void EvaluateVolume();
	void EvaluateCentroid();
	void EvaluateInertiaTensor();
	void EvaluateRadius();


	bool ReadCADFile();												//File reading of CAD
	int FindSpan(int &n, int &p, const double &u, double* U);
	void BasisFunctions(int &i, double &u, int &p, double* U, double* N);
	void DersBasisFunctions(int &i, const double &u, int &p, int &n, double* U, double** ders);
	void NURBSPoint(double &u, double &v, Matrix &point);
	void NURBSDerivatives(const double &uc, const double &vc, Matrix** &Skl, int &d);
	void WriteVTK_XMLRender(FILE *f, const Matrix& pos, const Matrix& rot, int number);
};

