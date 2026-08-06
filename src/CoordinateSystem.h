#pragma once
#include <stdio.h>

class Matrix;

class CoordinateSystem
{
public:
	CoordinateSystem();
	~CoordinateSystem();

	int number;		//ID do CS
	Matrix* E1;		//Versores
	Matrix* E2;
	Matrix* E3;

	Matrix* Q; //transformation matrix of coordinates

	//////////////Transformation of coordinates///////////////
	//                     UL = Q*UG                       //
	//    UL -> vector with coordinates in the system local     //
	//    UG -> vector with coordinates in the system global    //
	/////////////////////////////////////////////////////////
	bool Read(FILE *f);
	void Write(FILE *f);
};

