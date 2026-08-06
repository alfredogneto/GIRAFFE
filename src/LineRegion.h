#pragma once
#include <stdio.h>

class LineRegion
{
public:
	LineRegion();
	~LineRegion();

	int number;
	int n_elements;
	int* elements;

	bool Read(FILE *f);
	void Write(FILE *f);
};

