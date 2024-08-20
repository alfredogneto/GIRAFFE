#include "Interface_0.h"
#include"Database.h"
#include <math.h>
//Variáveis globais
extern
Database db;

Interface_0::Interface_0()
{
	epsn1 = 0.0;
	n1 = 0.0;
	zetan = 0.0;

	mus = 0.0;
	mud = 0.0;
	epst = 0.0;
	ct = 0.0;

	//Not used (but created)
	n2 = 0.0;
	gnb = 0.0;
	gnbb = 0.0;

	material_1 = 0;
	material_2 = 0;
}

Interface_0::~Interface_0()
{

}

//Checa inconsistências no elemento para evitar erros de execução
bool Interface_0::Check()
{
	if (material_1 > db.number_materials)
		return false;
	if (material_2 > db.number_materials)
		return false;
	return true;
}

bool Interface_0::Read(FILE *f)
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
	if (!strcmp(s, "EPN1"))
	{
		fscanf(f, "%s", s);
		epsn1 = atof(s);
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "N1"))
	{
		fscanf(f, "%s", s);
		n1 = atof(s);
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "ZetaN"))
	{
		fscanf(f, "%s", s);
		zetan = atof(s);
	}
	else
		return false;

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

void Interface_0::Write(FILE *f)
{
	fprintf(f, "Interface_0\t%d\tMaterials\t%d\t%d\nEPN1\t%.6e\nN1\t%.6e\nZetaN\t%.6e\nMUS\t%.6e\nMUD\t%.6e\nEPT\t%.6e\nCT\t%.6e\n",
		number, material_1, material_2, epsn1, n1, zetan, mus, mud, epst, ct);
}

//Pré-cálculo de variáveis que é feito uma única vez no início
void Interface_0::PreCalc()
{

}