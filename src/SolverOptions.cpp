#include "SolverOptions.h"
#ifdef GIRAFFE_USE_MKL
#include <mkl.h>
#endif
#include <omp.h>

#include"Database.h"
//Global variables
extern
Database db;

SolverOptions::SolverOptions()
{
	processors = 1;						//Number of processors - parallel processing - PARDISO
	solver = 0;							//Solver type - 0 diretor 1 iterativo - PARDISO
}


SolverOptions::~SolverOptions()
{
}

bool SolverOptions::Read(FILE *f)
{
	char s[1000];
	fscanf(f, "%s", s);
	if (!strcmp(s, "Processors"))
	{
		fscanf(f, "%s", s);
		processors = atoi(s);
	}
	else
		return false;
	fscanf(f, "%s", s);
	if (!strcmp(s, "LinSys"))
	{
		fscanf(f, "%s", s);
		if (!strcmp(s, "Direct"))
		{
			solver = 0;
		}
		else
		{
			if (!strcmp(s, "Iterative"))
			{
				solver = 1;
			}
			else
				return false;
		}
	}
	else
		return false;
	return true;
}
void SolverOptions::Write(FILE *f)
{
	char str[1000];
	if (solver == 0)
		strcpy(str, "Direct");
	if (solver == 1)
		strcpy(str, "Iterative");
	fprintf(f, "SolverOptions\tProcessors\t%d\tLinSys\t%s\n", processors,str);
}

void SolverOptions::PreCalc()
{
#ifdef GIRAFFE_USE_MKL
	if (processors != 1)
		MKL_Set_Num_Threads(processors - 1);
	else
		MKL_Set_Num_Threads(processors);
#endif
	omp_set_num_threads(processors);
}
