#pragma once
#include <stdio.h>
#include "Matrix.h"

class Point
{
public:
	Point();
	~Point();
	//[0] X - coordenada do ponto
	//[1] Y - coordenada do ponto
	//[2] Z - coordenada do ponto
	Matrix coordinates;			//Coordenadas do ponto
	int number;	//Numero do ponto
	bool Read(FILE *f);
	void Write(FILE *f);
};

