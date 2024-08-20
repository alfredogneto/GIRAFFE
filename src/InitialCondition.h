#pragma once
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "Matrix.h"

class InitialCondition
{
public:
	InitialCondition();
	~InitialCondition();
	bool Read(FILE *f);
	void Write(FILE *f);
	bool Check();
	int number;
	int node;
	int super_node;
	Matrix du;
	Matrix omega;

	int solution;						//ID da solution para aplicar o IC
	
	void ComputeInitialCondition();
};

