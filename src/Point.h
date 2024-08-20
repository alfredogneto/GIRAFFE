#pragma once
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "Matrix.h"
#include <math.h>
class Point
{
public:
	Point();
	~Point();
	//[0] X - coordenada do ponto
	//[1] Y - coordenada do ponto
	//[2] Z - coordenada do ponto
	Matrix coordinates;			//Coordenadas do ponto
	int number;	//Número do ponto
	bool Read(FILE *f);
	void Write(FILE *f);
};

