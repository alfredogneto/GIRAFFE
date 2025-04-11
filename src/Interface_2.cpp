#include "Interface_2.h"
#include"Database.h"
#include <math.h>
//Variáveis globais
extern
Database db;


Interface_2::Interface_2()
{
	epsn = 0.0;
	cn = 0.0;
	mus = 0.0;
	mud = 0.0;
	epst = 0.0;
	ct = 0.0;

	material_1 = 0;
	material_2 = 0;
}


Interface_2::~Interface_2()
{
}

//Checa inconsistências no elemento para evitar erros de execução
bool Interface_2::Check()
{
	if (material_1 > db.number_materials)
		return false;
	if (material_2 > db.number_materials)
		return false;
	return true;
}

bool Interface_2::Read(FILE * f)
{
	char s[1000];
	fscanf(f, "%s", s);
	number = atoi(s);

	//Materials
	fscanf(f, "%s", s);
	if (!strcmp(s, "Materials"))
	{
		fscanf(f, "%s", s);
		material_1 = atoi(s);
		fscanf(f, "%s", s);
		material_2 = atoi(s);
	}
	else
		return false;

	//Salva a posição (stream)
	fpos_t pos;
	fgetpos(f, &pos);

	fscanf(f, "%s", s);
	if (!strcmp(s, "MU"))
	{
		fscanf(f, "%s", s);
		mus = atof(s);
		mud = atof(s);
	}
	else
	{
		fsetpos(f, &pos);
		fscanf(f, "%s", s);
		if (!strcmp(s, "MUS"))
		{
			fscanf(f, "%s", s);
			mus = atof(s);
		}
		else
			return false;
		fscanf(f, "%s", s);
		if (!strcmp(s, "MUD"))
		{
			fscanf(f, "%s", s);
			mud = atof(s);
		}
		else
			return false;
	}

	fscanf(f, "%s", s);
	if (!strcmp(s, "EPN"))
	{
		fscanf(f, "%s", s);
		epsn = atof(s);
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "CN"))
	{
		fscanf(f, "%s", s);
		cn = atof(s);
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "EPT"))
	{
		fscanf(f, "%s", s);
		epst = atof(s);
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "CT"))
	{
		fscanf(f, "%s", s);
		ct = atof(s);
	}
	else
		return false;

	return true;
}

void Interface_2::Write(FILE *f)
{
	fprintf(f, "Interface_1\t%d\tMaterials\t%d\t%d\nMUS\t%.6e\nMUD\t%.6e\nEPN\t%.6e\nCN\t%.6e\nEPT\t%.6e\nCT\t%.6e\nFactor\t%.6e\n",
		number, material_1, material_2, mus, mud, epsn, cn, epst, ct);
}

//Pré-cálculo de variáveis que é feito uma única vez no início
void Interface_2::PreCalc()
{

}

double Interface_2::EvaluateElasticForce(double gap)
{
	return 0.0;
}

double Interface_2::EvaluateElasticStiffness(double gap)
{
	return 0.0;
}

double Interface_2::IntegralElasticForce(double gap_i, double gap_f)
{
	return 0.0;
}

double Interface_2::IntegralTrapElasticForce(double gap_i, double gap_f)
{
	return 0.0;
}
