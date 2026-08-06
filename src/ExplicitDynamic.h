#pragma once
#include "Solution.h"

class Matrix;

class ExplicitDynamic :
	public Solution
{
public:
	ExplicitDynamic();
	~ExplicitDynamic();
	bool Read(FILE *f);			//Reads input file
	void Write(FILE *f);		//Writes output file
	bool Solve();				//Solves solution routine
	
	void WriteResults();		//Writes results according to a amostragem specified
	void SetGlobalSizeExplicit();				//Sets tamanhos of contributions global - explicito
	void MountExplicit();
	void Euler();
	void RungeKutta4();
	void MountLoadsExplicit(double t);
	void MountContactsExplicit(double t);
	void MontSpecialConstraintsExplicit(double t);
	void MountDisplacementsExplicit(double t);
	void InitialEvaluations();
	void MountGlobalExplicit();
	void FinalUpdateContactsExplicit(double t);

	//Parameters
	double i_time_step;
	double max_time_step;
	double min_time_step;
	char* method;	//Name of the type of the method
	
	int sample;
	//Damping
	double alpha;
	double beta;
	int update;
	

	//Internal control variables
	int file_index;
	double time;
	double time_step;

	bool zero_IC_flag;		//Zeros the ICs from previous step - dafult false

	Matrix* I3;
};

