#pragma once
#include "Solution.h"

class Dynamic :
	public Solution
{
public:
	Dynamic();
	~Dynamic();
	bool Read(FILE *f);			//Reads input file
	void Write(FILE *f);		//Writes output file
	bool Solve();				//Solves solution routine
	void UpdateDyn();			//Atualizações - Newmark
	void WriteResults(double time_value);		//Escreve resultados de acordo com uma amostragem especificada

	//Parameters
	double i_time_step;
	double max_time_step;
	double min_time_step;
	int max_it;
	int min_it;
	int conv_increase;
	double inc_factor;
	int sample;
	//Damping
	double alpha;
	double beta;
	int update;
	//Newmark
	double gamma_new;
	double beta_new;
	double a1, a2, a3, a4, a5, a6;
	void CalculateNewmarkCoeff(double time_step);

	//Flag to indicate existence of contact-inpact control
	bool contact_impact_control;
	int n_steps_impact;

	//Internal control variables
	int file_index;

	bool zero_IC_flag;		//Zeros the ICs from previous step - dafult false
};

