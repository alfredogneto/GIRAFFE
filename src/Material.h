#pragma once
#include <stdio.h>

class Material
{
public:
	Material() {}
	virtual ~Material() {}

	int number;		//ID of the material
	double rho;		//Mass specifies [M/L3]

	virtual bool Read(FILE *f) = 0;
	virtual void Write(FILE *f) = 0;
};

