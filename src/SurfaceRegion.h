#pragma once
#include <stdio.h>

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

