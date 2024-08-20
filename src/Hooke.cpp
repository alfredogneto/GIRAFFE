#include "Hooke.h"


Hooke::Hooke()
{
	number = 0;
	E = 0;
	nu = 0;
}


Hooke::~Hooke()
{
}

bool Hooke::Read(FILE *f)
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
	return true;
}

void Hooke::Write(FILE *f)
{
	fprintf(f, "Hooke\t%d\tE\t%.6e\tNu\t%.6e\tRho\t%.6e\n",
		number,
		E,
		nu,
		rho);
}
