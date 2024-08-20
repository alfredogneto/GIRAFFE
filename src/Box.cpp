#include "Box.h"


Box::Box()
{
	for (int i = 0; i < 8; i++)
	{
		x[i] = 0;
		y[i] = 0;
		z[i] = 0;
	}
	xmax = 0;
	ymax = 0;
	zmax = 0;
	xmin = 0;
	ymin = 0;
	zmin = 0;

	center[0] = 0;
	center[1] = 0;
	center[2] = 0;
}

Box::~Box()
{
}

double Box::Get_Diagonal_Size()
{
	double ret = sqrt((xmin - xmax)*(xmin - xmax) + (ymin - ymax)*(ymin - ymax) + (zmin - zmax)*(zmin - zmax));
	return ret;
}

//seta vértices, centro  e prepara retorno de máximos e mínimos
void Box::SetVertices(double* e_x, double* e_y, double* e_z)
{
	center[0] = 0;
	center[1] = 0;
	center[2] = 0;
	for (int i = 0; i < 8; i++)
	{
		x[i] = e_x[i];
		y[i] = e_y[i];
		z[i] = e_z[i];
		center[0] += x[i];
		center[1] += y[i];
		center[2] += z[i];
	}
	center[0] = center[0] / 8;
	center[1] = center[1] / 8;
	center[2] = center[2] / 8;
	//Calcula os valores máximos e deixa-os prontos para serem exportados, quando requeridos
	xmax = x[0];
	xmin = x[0];
	ymax = y[0];
	ymin = y[0];
	zmax = z[0];
	zmin = z[0];
	for (int i = 1; i < 8; i++)
	{
		if (x[i] > xmax)
			xmax = x[i];
		if (x[i] < xmin)
			xmin = x[i];

		if (y[i] > ymax)
			ymax = y[i];
		if (y[i] < ymin)
			ymin = y[i];

		if (z[i] > zmax)
			zmax = z[i];
		if (z[i] < zmin)
			zmin = z[i];
	}
	 

}

//returns x max
double Box::Get_x_Max(double inflation_factor)
{
	return (center[0] + (xmax-center[0])*inflation_factor);
}
//returns y max
double Box::Get_y_Max(double inflation_factor)
{
	return (center[1] + (ymax - center[1])*inflation_factor);
}
//returns z max
double Box::Get_z_Max(double inflation_factor)
{
	return (center[2] + (zmax - center[2])*inflation_factor);
}
//returns x min
double Box::Get_x_Min(double inflation_factor)
{
	return (center[0] + (xmin - center[0])*inflation_factor);
}
//returns y min
double Box::Get_y_Min(double inflation_factor)
{
	return (center[1] + (ymin - center[1])*inflation_factor);
}
//returns z min
double Box::Get_z_Min(double inflation_factor)
{
	return (center[2] + (zmin - center[2])*inflation_factor);
}


bool BoxOverlap(Box& b1, Box & b2, double inflation_factor)
{
	if (b1.Get_x_Max(inflation_factor) < b2.Get_x_Min(inflation_factor)) return false;
	if (b1.Get_x_Min(inflation_factor) > b2.Get_x_Max(inflation_factor)) return false;

	if (b1.Get_y_Max(inflation_factor) < b2.Get_y_Min(inflation_factor)) return false;
	if (b1.Get_y_Min(inflation_factor) > b2.Get_y_Max(inflation_factor)) return false;

	if (b1.Get_z_Max(inflation_factor) < b2.Get_z_Min(inflation_factor)) return false;
	if (b1.Get_z_Min(inflation_factor) > b2.Get_z_Max(inflation_factor)) return false;

	return true; // boxes present overlap
}