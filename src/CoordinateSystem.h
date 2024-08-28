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

	Matrix* Q; //matriz de transformação de coordenadas

	//////////////Transformação de coordenadas///////////////
	//                     UL = Q*UG                       //
	//    UL -> vetor com coordenadas no sistema local     //
	//    UG -> vetor com coordenadas no sistema global    //
	/////////////////////////////////////////////////////////
	bool Read(FILE *f);
	void Write(FILE *f);
};

