#pragma once
#include "Solution.h"
class ExplicitDynamic :
	public Solution
{
public:
	ExplicitDynamic();
	~ExplicitDynamic();
	bool Read(FILE *f);			//Reads input file
	void Write(FILE *f);		//Writes output file
	bool Solve();				//Solves solution routine
	
	void WriteResults();		//Escreve resultados de acordo com uma amostragem especificada
	void SetGlobalSizeExplicit();				//Seta tamanhos de contribuições globais - explícito
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
	char* method;	//Nome do tipo do método
	
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

