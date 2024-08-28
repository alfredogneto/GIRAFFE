#include "PipeSection.h"
#include <string>

#include "SolidSection.h"


#define PI 3.1415926535897932384626433832795

PipeSection::PipeSection()
{
	sec_details = new SolidSection();
}

PipeSection::~PipeSection()
{
	delete sec_details;
}

bool PipeSection::Read(FILE *f)
{
	char s[1000];
	fscanf(f, "%s", s);
	//Verifica a palavra chave "PS"
	if (!strcmp(s, "PS"))
	{
		fscanf(f, "%s", s);
		number = atoi(s);
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "EA"))
	{
		fscanf(f, "%s", s);
		EA = atof(s);
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "EI"))
	{
		fscanf(f, "%s", s);
		EI = atof(s);
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "GJ"))
	{
		fscanf(f, "%s", s);
		GJ = atof(s);
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "GA"))
	{
		fscanf(f, "%s", s);
		GA = atof(s);
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "Rho"))
	{
		fscanf(f, "%s", s);
		Rho = atof(s);
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "CDt"))
	{
		fscanf(f, "%s", s);
		CDt = atof(s);
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "CDn"))
	{
		fscanf(f, "%s", s);
		CDn = atof(s);
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "CAt"))
	{
		fscanf(f, "%s", s);
		CAt = atof(s);
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "CAn"))
	{
		fscanf(f, "%s", s);
		CAn = atof(s);
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "De"))
	{
		fscanf(f, "%s", s);
		De = atof(s);
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "Di"))
	{
		fscanf(f, "%s", s);
		Di = atof(s);
	}
	else
		return false;
	return true;
}

void PipeSection::Write(FILE *f)
{
	fprintf(f, "PS\t%d\tEA\t%.6e\tEI\t%.6e\tGJ\t%.6e\tGA\t%.6e\tRho\t%.6e\tCDt\t%.6e\tCDn\t%.6e\tCAt\t%.6e\tCAn\t%.6e\tDe\t%.6e\tDi\t%.6e\n", number, EA, EI, GJ, GA, Rho, CDt, CDn, CAt, CAn, De, Di);
}

void PipeSection::PreCalc()
{
	int n_circ = 24;
	sec_details->Alloc(n_circ);
	double r = De / 2;
	double theta = 0;
	for (int index = 0; index < n_circ; index++)
	{
		theta = (index * 2 * PI) / n_circ;
		sec_details->points[index][0] = r*cos(theta);
		sec_details->points[index][1] = r*sin(theta);
	}
}
