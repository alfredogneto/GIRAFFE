#include "Orthotropic.h"



Orthotropic::Orthotropic()
{
	number = 0;
	E1 = 0;
	E2 = 0;
	G12 = 0;
	G23 = 0;
	nu12 = 0;
}


Orthotropic::~Orthotropic()
{
}

bool Orthotropic::Read(FILE *f) 
{
	char s[1000];
	fscanf(f, "%s", s);
	number = atoi(s);

	fscanf(f, "%s", s);
	if (!strcmp(s, "E1"))
	{
		fscanf(f, "%s", s);
		E1 = atof(s);
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "E2"))
	{
		fscanf(f, "%s", s);
		E2 = atof(s);
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "G12"))
	{
		fscanf(f, "%s", s);
		G12 = atof(s);
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "G23"))
	{
		fscanf(f, "%s", s);
		G23 = atof(s);
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "Nu12"))
	{
		fscanf(f, "%s", s);
		nu12 = atof(s);
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

void Orthotropic::Write(FILE *f)
{
	fprintf(f, "Orthotropic\t%d\tE1\t%.6e\tE2\t%.6e\tG12\t%.6e\tG23\t%.6e\tNu12\t%.6e\tRho\t%.6e\n",
		number,
		E1,
		E2,
		G12,
		G23,
		nu12,
		rho);
}
