#include "ArcCirc.h"
//#include <stdlib.h>
#include <string.h>
#include <math.h>
#define PI 3.1415926535897932384626433832795

ArcCirc::ArcCirc()
{
	i_point = Matrix(2);
	f_point = Matrix(2);		
	c_point = Matrix(2);
	theta_i = 0.0;
	theta_f = 0.0;
	radius = 0.0;
}


ArcCirc::~ArcCirc()
{
}

//Leitura
bool ArcCirc::Read(FILE *f)
{
	char s[1000];
	fscanf(f, "%s", s);
	//Verifica a palavra chave "Arc"
	if (!strcmp(s, "Arc"))
	{
		fscanf(f, "%s", s);
		number = atoi(s);
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "InitialPoint"))
	{
		fscanf(f, "%s", s);
		i_point(0, 0) = atof(s);
		fscanf(f, "%s", s);
		i_point(1, 0) = atof(s);
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "FinalPoint"))
	{
		fscanf(f, "%s", s);
		f_point(0, 0) = atof(s);
		fscanf(f, "%s", s);
		f_point(1, 0) = atof(s);
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "CenterPoint"))
	{
		fscanf(f, "%s", s);
		c_point(0, 0) = atof(s);
		fscanf(f, "%s", s);
		c_point(1, 0) = atof(s);
	}
	else
		return false;
	return true;
}

//Gravação
void ArcCirc::Write(FILE *f)
{
	fprintf(f, "Arc\t%d\tInitialPoint\t%.6e\t%.6e\tFinalPoint\t%.6e\t%.6e\tCenterPoint\t%.6e\t%.6e\n", number, i_point(0, 0), i_point(1, 0), f_point(0, 0), f_point(1, 0), c_point(0, 0), c_point(1, 0));
}

//Pré-cálculo
void ArcCirc::PreCalc()
{
	theta_i = ArcReduction(atan2(i_point(1, 0) - c_point(1, 0), i_point(0, 0) - c_point(0, 0)));
	theta_f = ArcReduction(atan2(f_point(1, 0) - c_point(1, 0), f_point(0, 0) - c_point(0, 0)));
	radius = sqrt((i_point(1, 0) - c_point(1, 0))*(i_point(1, 0) - c_point(1, 0)) + (i_point(0, 0) - c_point(0, 0))*(i_point(0, 0) - c_point(0, 0)));
	//printf("thetai %.6e thetaf %.6e radius %.6e\n", theta_i, theta_f, radius);
}

void ArcCirc::BoundingRadiusAndCenter(MatrixFloat& center, float* e_radius)
{
	//Evaluates the center and the radius, to be used to produce bounding cylinders evaluation
	if (AngularRange() >= PI)
	{
		center(0, 0) = (float)c_point(0, 0);
		center(1, 0) = (float)c_point(1, 0);
		center(2, 0) = 0.0f;
		*e_radius = (float)radius;
	}
	else
	{
		center(0, 0) = (float)(0.5*(i_point(0, 0) + f_point(0, 0)));
		center(1, 0) = (float)(0.5*(i_point(1, 0) + f_point(1, 0)));
		center(2, 0) = 0.0f;
		*e_radius = (float)(0.5 * norm(f_point - i_point));
	}
}

void ArcCirc::Center(Matrix & center)
{
	//Evaluates the center with double precision
	if (AngularRange() >= PI)
	{
		center(0, 0) = c_point(0, 0);
		center(1, 0) = c_point(1, 0);
		center(2, 0) = 0.0;
	}
	else
	{
		center(0, 0) = (0.5*(i_point(0, 0) + f_point(0, 0)));
		center(1, 0) = (0.5*(i_point(1, 0) + f_point(1, 0)));
		center(2, 0) = 0.0;
	}
}


void ArcCirc::BoundingRectangle(float* x_min, float* x_max, float* y_min, float* y_max)
{
	//Evaluates bounding rectangle based on arc locus
	double temp;
	//Initial point:
	temp = radius * cos(theta_i);
	double temp_x_min = temp;
	double temp_x_max = temp;
	temp = radius * sin(theta_i);
	double temp_y_min = temp;
	double temp_y_max = temp;
	//Final point:
	temp = radius * cos(theta_f);
	if (temp < temp_x_min)
		temp_x_min = temp;
	if (temp > temp_x_max)
		temp_x_max = temp;
	temp = radius * sin(theta_f);
	if (temp < temp_y_min)
		temp_y_min = temp;
	if (temp > temp_y_max)
		temp_y_max = temp;
	// 0 rad
	if (InsideArc(0.0))
	{
		temp = radius;		//radius * cos (0.0)
		if (temp < temp_x_min)
			temp_x_min = temp;
		if (temp > temp_x_max)
			temp_x_max = temp;
		temp = 0.0;			//radius * sin (0.0)
		if (temp < temp_y_min)
			temp_y_min = temp;
		if (temp > temp_y_max)
			temp_y_max = temp;
	}
	// PI/2 rad
	if (InsideArc(0.5 * PI))
	{
		temp = 0.0;			//radius * cos (PI/2)
		if (temp < temp_x_min)
			temp_x_min = temp;
		if (temp > temp_x_max)
			temp_x_max = temp;
		temp = radius;		//radius * sin (PI/2)
		if (temp < temp_y_min)
			temp_y_min = temp;
		if (temp > temp_y_max)
			temp_y_max = temp;
	}
	// PI rad
	if (InsideArc(PI))
	{
		temp = -radius;		//radius * cos (PI)
		if (temp < temp_x_min)
			temp_x_min = temp;
		if (temp > temp_x_max)
			temp_x_max = temp;
		temp = 0.0;			//radius * sin (PI)
		if (temp < temp_y_min)
			temp_y_min = temp;
		if (temp > temp_y_max)
			temp_y_max = temp;
	}
	// 3 PI / 2 rad
	if (InsideArc(1.5 * PI))
	{
		temp = 0.0;			//radius * cos (3 PI / 2)
		if (temp < temp_x_min)
			temp_x_min = temp;
		if (temp > temp_x_max)
			temp_x_max = temp;
		temp = -radius;		//radius * sin (3 PI / 2)
		if (temp < temp_y_min)
			temp_y_min = temp;
		if (temp > temp_y_max)
			temp_y_max = temp;
	}
	//Final attributions:
	*x_min = (float)(temp_x_min + c_point(0, 0));
	*x_max = (float)(temp_x_max + c_point(0, 0));
	*y_min = (float)(temp_y_min + c_point(1, 0));
	*y_max = (float)(temp_y_max + c_point(1, 0));
}

bool ArcCirc::InsideArc(double theta)
{
	//Redução de theta ao range [-PI,PI]
	double theta_mod = ArcReduction(theta);
	if (theta_i <= theta_f)
	{
		if (theta_mod >= theta_i && theta_mod <= theta_f)
			return true;
		else
			return false;
	}
	else
	{
		if (theta_mod <= theta_f || theta_mod >= theta_i)
			return true;
		else
			return false;
	}
}

double ArcCirc::AngularRange()
{
	if (theta_f >= theta_i)
		return abs(theta_f - theta_i);
	else
		return (2 * PI - abs(theta_f - theta_i));
}

double ArcCirc::CentralTheta()
{
	if (theta_i <= theta_f)
	{
		return (0.5 * (theta_i + theta_f));
	}
	else
	{
		return (0.5 * (theta_i + theta_f) + PI);
	}
}
