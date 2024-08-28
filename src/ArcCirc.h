#pragma once
#include <stdio.h>
#include "Matrix.h"
#include "MatrixFloat.h"


class ArcCirc
{
public:
	ArcCirc();
	~ArcCirc();
	bool Read(FILE *f);					//Leitura
	void Write(FILE *f);				//Gravação
	void PreCalc();						//Pré-cálculo

	int number;							//ID do arco
	Matrix i_point;					//ponto de início do arco - plano xy
	Matrix f_point;					//ponto de fim do arco - plano xy
	Matrix c_point;					//ponto de centro de curvatura do arco - plano xy

	//variáveis internas
	double theta_i;
	double theta_f;
	double radius;

	void Center(Matrix & center);
	void BoundingRadiusAndCenter(MatrixFloat& center, float* e_radius);
	void BoundingRectangle(float* x_min, float* x_max, float* y_min, float* y_max);
	bool InsideArc(double theta);
	double AngularRange();
	double CentralTheta();
};

