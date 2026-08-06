#pragma once
#include <stdio.h>
#include "Matrix.h"

class Point
{
public:
	Point();
	~Point();
	//[0] X - coordinate of the point
	//[1] Y - coordinate of the point
	//[2] Z - coordinate of the point
	Matrix coordinates;			//Coordinates of the point
	int number;	//Number of the point
	bool Read(FILE *f);
	void Write(FILE *f);
};

