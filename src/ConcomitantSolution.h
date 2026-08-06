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

	bool Read(FILE *f);					//Reading
	void Write(FILE *f);				//Recording
	bool Check();						//Checking inconsistencies

	int sample;

protected:
	//variables internal
	BoolTable	bool_concomitant;
	Solution* sol;
	FILE *f_output;
};

