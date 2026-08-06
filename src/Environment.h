#pragma once
#include <stdio.h>

#include "Matrix.h"
#include "BoolTable.h"

class Environment
{
public:
	Environment(void);
	~Environment(void);
	void flush();															//Deletes alocações of the memory
	bool Read(FILE *f);
	void Write(FILE *f);

	//Variables booleanas for controlar the that exists of data ambientais
	bool g_exist;
	bool ocean_data_exist;
	bool wind_data_exist;

	//BoolTables
	BoolTable bool_g;																//Sequencial of passos for the weight próprio
	BoolTable bool_current;															//Sequencial of passos for the corrente
	
	//Gravity data
	Matrix G;																//Field gravitacional
	
	//Sea current data
	void SetNCurrentPoints(int value);										//Sets the number of points of data of correnteza maritima
	void SetCurrentData(int line, double depth, double speed, double alpha);	//Sets in the linha "line" the variables referring to the correnteza
	Matrix VelocityAt(double H);											//Performs interpolation linear and returns the value of the velocity in certa profundidade
	double GetDepthAtIndex(int index);										//Returns the Depth in the due index								
	double GetSpeedAtIndex(int index);										//Returns the Speed in the due index
	double GetAngleAtIndex(int index);										//Returns the Angle in the due index
	double rho_fluid;														//Mass specifies of the agua
	Matrix surface_position;												//Position of the lamina d'agua in the ambiente
	int n_current_points;													//Number of points of data of the correnteza maritima
	double** current;														//Matrix that stores the data referring to the correnteza maritima
	
	//WindData
	Matrix TimeAt(double time);												//Performs interpolation linear in the time - returns the matrix with all the data already interpolados in the time
	Matrix WindVelocityAt(Matrix pos, double time);							//Performs interpolation linear in the espaço and time - returns vector velocity in the point of interest
	double GetTimeAt(int index);											//Returns the Time in the due index
	double GetWindSpeedAt(int index);										//Returns the Wind Speed in the due index
	double GetDeltaAt(int index);											//Returns the Delta in the due index
	double GetWindVertAt(int index);										//Returns the Wind Vert in the due index
	double GetHorizShearAt(int index);										//Returns the Shear Linear horizontal in the due index
	double GetVertPLShearAt(int index);										//Returns the Power Law Shear vertical in the due index
	double GetVertLinShearAt(int index);									//Returns the Shear Linear vertical in the due index
	double GetGustSpeedAt(int index);										//Returns the Gust Speed in the due index

	void SetWindData(int line, double time, double wind_speed, double delta, double wind_vert, double HSHR, double VSHR, double Lin_VSHR, double gust_speed);		//Sets in the linha "line" the variables referring to to the wind

	/*
	wind[line][0] = time;
	wind[line][1] = wind_speed;
	wind[line][2] = delta*PI / 180;
	wind[line][3] = wind_vert;
	wind[line][4] = HSHR;
	wind[line][5] = VSHR;
	wind[line][6] = Lin_VSHR;
	wind[line][7] = gust_speed;
	
	*/
	int cs;																	//ID of the coordinate system for define the wind
	
	double VHor(double proj, double zlocal, double zhub, double wind_speed, double HSHR, double Lin_VSHR, double VSHR, double Vgust);	//Performs calculation of the velocity horizontal
	double rho_air;															//Mass specifies of the ar
	Matrix reference_position;												//Position of the incidência of the wind - system global
	Matrix transform3;														//Transformation matrix of coordinates 3x3 - local-global
	bool transform3_calculated;												//Variable boolean that indicates that the transformation matrix of coordinates already was calculated
	
	double** wind;															//Matrix that stores the data referring to to the wind (all the data of reading - for each instant of time)
	
	int n_wind_points;														//Number of points of data of the wind (number of instantes of time)
	void SetNWindPoints(int value);											//Sets the number of points of data of the ar (number of instantes - allocation of memory)
	
};
