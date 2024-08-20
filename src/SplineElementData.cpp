#include "SplineElementData.h"

SplineElementData::SplineElementData(int n_solutions)
{
	n_sol = n_solutions;
	G_p = new Matrix*[n_solutions];
	for (int i = 0; i < n_solutions; i++)
		G_p[i] = new Matrix(3);
}

SplineElementData::~SplineElementData()
{
	if (G_p != NULL)
	{
		for (int i = 0; i < n_sol; i++)
			delete G_p[i];
		delete[]G_p;
	}	
}