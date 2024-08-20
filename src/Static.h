#pragma once
#include "Solution.h"

class Static :
	public Solution
{
public:
	Static();
	~Static();
	bool Read(FILE *f);									//Reads input file
	void Write(FILE *f);								//Writes output file
	bool Solve();										//Solves solution routine
	void WriteResults(double time_value);				//Escreve resultados de acordo com uma amostragem especificada
	
	double i_time_step;
	double max_time_step;
	double min_time_step;
	int max_it;
	int min_it;
	int conv_increase;
	double inc_factor;
	int sample;

	//Internal control variables
	int file_index;
};

