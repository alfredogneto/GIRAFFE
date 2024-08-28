#pragma once
#include "CADData.h"

class NURBSSurface :
	public CADData
{
public:
	NURBSSurface();
	~NURBSSurface();

	
	//Dados NURBS
	int U_dim;
	int U_order;
	int V_dim;
	int V_order;
	double* U_knot_vector;
	double* V_knot_vector;
	double** weights;
	Matrix** control_points;
	Matrix** Pw;	//control points em 4D (ja com pesos)
	double** Bin;	//Binomial coefficients

	bool Read(FILE *f);												//Leitura do arquivo de entrada
	void Write(FILE *f);											//Escrita do arquivo de saida
	void PreCalc();													//PreCalc

	//Geometric evaluation functions
	void EvaluateVolume();
	void EvaluateCentroid();
	void EvaluateInertiaTensor();
	void EvaluateRadius();


	bool ReadCADFile();												//Leitura do arquivo de CAD
	int FindSpan(int &n, int &p, double &u, double* U);
	void BasisFunctions(int &i, double &u, int &p, double* U, double* N);
	void DersBasisFunctions(int &i, double &u, int &p, int &n, double* U, double** ders);
	void NURBSPoint(double &u, double &v, Matrix &point);
	void NURBSDerivatives(double &uc, double &vc, Matrix** &Skl, int &d);
	void WriteVTK_XMLRender(FILE *f, Matrix& pos, Matrix& rot, int number);
};

