#pragma once
#include <stdio.h>
#include "Matrix.h"
#include "MatrixFloat.h"


class ArcCirc
{
public:
	ArcCirc();
	~ArcCirc();
	bool Read(FILE *f);					//Reading
	void Write(FILE *f);				//Recording
	void PreCalc();						//Pre-calculation

	int number;							//ID of the arc
	Matrix i_point;					//point of start of the arc - plane xy
	Matrix f_point;					//point of end of the arc - plane xy
	Matrix c_point;					//point of center of curvature of the arc - plane xy

	//variables internal
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

