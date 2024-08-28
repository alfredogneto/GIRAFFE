#pragma once
#include <stdio.h>

class SolverOptions
{
public:
	SolverOptions();
	~SolverOptions();
	int processors;						//Número de processadores - processamento paralelo - PARDISO
	int solver;							//Tipo de solver - 0 diretor 1 iterativo - PARDISO
	void PreCalc();						//Setting parallel processing and solver options
	bool Read(FILE *f);
	void Write(FILE *f);
};

