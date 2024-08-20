#pragma once
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "Matrix.h"

class SurfaceRegion
{
public:
	SurfaceRegion();
	~SurfaceRegion();

	int number;
	int n_elements;
	int* elements;

	bool Read(FILE *f);
	void Write(FILE *f);
};

