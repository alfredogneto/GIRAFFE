#pragma once

//This class contains data about a bounding box of a given 3D object
class Box
{
public:
	Box();
	~Box();
	
	double Get_x_Max(double inflation_factor);		//returns x max
	double Get_y_Max(double inflation_factor);		//returns y max
	double Get_z_Max(double inflation_factor);		//returns z max

	double Get_x_Min(double inflation_factor);		//returns x min
	double Get_y_Min(double inflation_factor);		//returns y min
	double Get_z_Min(double inflation_factor);		//returns z min

	double Get_Diagonal_Size();

	void SetVertices(double* e_x, double* e_y, double* e_z);					//seta vertices e prepara retorno de maximos e minimos

protected:
	//Coordinates of 8 vertices of box
	double x[8];
	double y[8];
	double z[8];

	double center[3];	//Coordenadas do centro geometrico dos 8 vertices do bounding box

	double xmax, ymax, zmax, xmin, ymin, zmin;
};

bool BoxOverlap(Box& b1, Box & b2, double inflation_factor);			//Checa overlap entre dois boxes

