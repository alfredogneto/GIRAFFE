#pragma once
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "Matrix.h"
#include <math.h>
#include "BoolTable.h"
#include "Solution.h"

class ConcomitantSolution
{
public:
	ConcomitantSolution();
	~ConcomitantSolution();
	void StartConcomitantSolution();
	void UpdateConcomitantSolution(double time);
	void EndConcomitantSolution();

	bool Read(FILE *f);					//Leitura
	void Write(FILE *f);				//Grava��o
	bool Check();						//Checking inconsistencies

	int sample;

protected:
	//vari�veis internas
	BoolTable	bool_concomitant;
	Solution* sol;
	FILE *f_output;
};

