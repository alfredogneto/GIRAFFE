#pragma once
#include <stdio.h>

class Material
{
public:
	Material() {}
	virtual ~Material() {}

	int number;		//ID do material
	double rho;		//Massa especifica [M/L3]

	virtual bool Read(FILE *f) = 0;
	virtual void Write(FILE *f) = 0;
};

