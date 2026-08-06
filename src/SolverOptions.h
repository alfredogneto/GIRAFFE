#pragma once
#include <stdio.h>

class SolverOptions
{
public:
	SolverOptions();
	~SolverOptions();
	int processors;						//Number of processors - parallel processing - PARDISO
	int solver;							//Solver type - 0 diretor 1 iterativo - PARDISO
	void PreCalc();						//Setting parallel processing and solver options
	bool Read(FILE *f);
	void Write(FILE *f);
};

