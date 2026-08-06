#include "ElasticPlasticIsoHardening.h"
#include <string>


ElasticPlasticIsoHardening::ElasticPlasticIsoHardening()
{
}


ElasticPlasticIsoHardening::~ElasticPlasticIsoHardening()
{
}


bool ElasticPlasticIsoHardening::Read(FILE *f)
{
	char s[1000];
	fscanf(f, "%s", s);
	number = atoi(s);

	fscanf(f, "%s", s);
	if (!strcmp(s, "E"))
	{
		fscanf(f, "%s", s);
		E = atof(s);
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "Nu"))
	{
		fscanf(f, "%s", s);
		nu = atof(s);
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "Rho"))
	{
		fscanf(f, "%s", s);
		rho = atof(s);
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "H"))
	{
		fscanf(f, "%s", s);
		H = atof(s);
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "YieldingStrength"))
	{
		fscanf(f, "%s", s);
		sigma_y_0 = atof(s);
	}
	else
		return false;

	return true;
}

void ElasticPlasticIsoHardening::Write(FILE *f)
{
	fprintf(f, "ElasticPlasticIsoHardening\t%d\tE\t%.6e\tNu\t%.6e\tRho\t%.6e\tH\t%.6e\tYieldingStrength\t%.6e\n",
		number,
		E,
		nu,
		rho,
		H,
		sigma_y_0);
}