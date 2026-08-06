#pragma once
#include <stdio.h>
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

	int solution;						//ID of the solution for apply the IC
	
	void ComputeInitialCondition();
};

