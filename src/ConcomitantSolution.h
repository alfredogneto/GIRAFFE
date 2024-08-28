#pragma once
#include <stdio.h>
#include "BoolTable.h"

class Solution;

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

