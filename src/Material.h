#pragma once
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
class Material
{
public:
	Material();
	virtual ~Material();

	int number;		//ID do material
	double rho;		//Massa específica [M/L3]

	virtual bool Read(FILE *f) = 0;
	virtual void Write(FILE *f) = 0;
};

