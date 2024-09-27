#include "Interface_1.h"
#include"Database.h"
#include <math.h>
//Variaveis globais
extern
Database db;

Interface_1::Interface_1()
{
	epsn1 = 0.0;
	n1 = 0.0;
	gnb = 0.0;
	gnbb = 0.0;
	factor = 0.0;
	n2 = 0.0;
	zetan = 0.0;

	epsn2 = 0.0;
	c2 = 0.0;

	mus = 0.0;
	mud = 0.0;
	epst = 0.0;
	ct = 0.0;
	
	material_1 = 0;
	material_2 = 0;
}

Interface_1::~Interface_1()
{

}

//Checa inconsistências no elemento para evitar erros de execução
bool Interface_1::Check()
{
	if (material_1 > db.number_materials)
		return false;
	if (material_2 > db.number_materials)
		return false;
	return true;
}

bool Interface_1::Read(FILE *f)
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
	if (!strcmp(s, "N2"))
	{
		fscanf(f, "%s", s);
		n2 = atof(s);
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "GNB"))
	{
		fscanf(f, "%s", s);
		gnb = atof(s);
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "Factor"))
	{
		fscanf(f, "%s", s);
		factor = atof(s);
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

void Interface_1::Write(FILE *f)
{
	fprintf(f, "Interface_1\t%d\tMaterials\t%d\t%d\nEPN1\t%.6e\nN1\t%.6e\nN2\t%.6e\nGNB\t%.6e\nFactor\t%.6e\nZetaN\t%.6e\nMUS\t%.6e\nMUD\t%.6e\nEPT\t%.6e\nCT\t%.6e\n", 
		number, material_1, material_2, epsn1, n1, n2, gnb, factor, zetan, mus, mud, epst, ct);
}

//Pre-calculo de variaveis que e feito uma unica vez no inicio
void Interface_1::PreCalc()
{
	gnbb = factor * gnb;
	epsn2 = (-epsn1 * n1 * pow(gnb - gnbb, n1 - 1)) / (n2 *  pow(gnbb, n2 - 1));
	c2 = epsn1 * pow(gnb - gnbb, n1) - epsn2 * pow(gnbb, n2);
}

double Interface_1::EvaluateElasticForce(double gap)
{
	/*
	\[Epsilon]2 = (-\[Epsilon]1 n1 (gb - gbb)^(n1 - 1))/(n2 gbb^(n2 - 1));
	c2 = \[Epsilon]1 (gb - gbb)^n1 - \[Epsilon]2 gbb^n2;
	fn1 = \[Epsilon]1 (gb - g)^n1;
	fn2 = \[Epsilon]2 (g)^n2 + c2;
	*/
	double fn = 0.0;
	
	if (gap > gnb)
	{
		fn = 0.0;
	}
	else
	{
		if (gap > gnbb)
		{
			fn = epsn1 * pow(gnb - gap, n1);
		}
		else
		{
			fn = epsn2 * pow(gap,n2) + c2;
		}
	}
	return fn;
}
double Interface_1::EvaluateElasticStiffness(double gap)
{
	double dfndgap = 0.0;
	//Returns a positive value (since dfndgap is negative), return its opposite value
	if (gap > gnb)
	{
		dfndgap = 0.0;
	}
	else
	{
		if (gap > gnbb)
		{
			dfndgap = epsn1 * n1 * pow(gnb - gap, n1 - 1.0);
		}
		else
		{
			dfndgap = - epsn2 * n2 * pow(gap, n2 - 1.0);
		}
	}
	return dfndgap;
}

double Interface_1::IntegralElasticForce(double gap_i, double gap_f)
{
	double result = 0.0;
	double multiplier;
	double g_i, g_f;
	if (gap_i > gap_f)
	{
		g_i = gap_i;
		g_f = gap_f;
		multiplier = +1.0;
	}
		
	else
	{
		g_i = gap_f;
		g_f = gap_i;
		multiplier = -1.0;
	}
		
	//-(epsn1 / (n1 + 1.0))*pow(gnb - gap_f, n1 + 1.0) + (epsn1 / (n1 + 1.0))*pow(gnb - gap_i, n1 + 1.0);
	//(epsn2 / (n2 + 1.0))*pow(gap_f, n2 + 1.0) + c2 * gap_f - (epsn2 / (n2 + 1.0))*pow(gap_i, n2 + 1.0) - c2 * gap_i;
	if (g_i > gnb)
	{
		if (g_f > gnb)
		{
			result = 0.0;
		}
		else
		{
			if (g_f > gnbb)
				result = -(epsn1 / (n1 + 1.0))*pow(gnb - g_f, n1 + 1.0);
			else
				result = -(epsn1 / (n1 + 1.0))*pow(gnb - gnbb, n1 + 1.0) +
				(epsn2 / (n2 + 1.0))*pow(g_f, n2 + 1.0) + c2 * g_f - (epsn2 / (n2 + 1.0))*pow(gnbb, n2 + 1.0) - c2 * gnbb;
		}
		
	}
	else
	{
		if (g_i > gnbb)
		{
			if (g_f > gnbb)
				result = -(epsn1 / (n1 + 1.0))*pow(gnb - g_f, n1 + 1.0) + (epsn1 / (n1 + 1.0))*pow(gnb - g_i, n1 + 1.0);
			else
				result = -(epsn1 / (n1 + 1.0))*pow(gnb - gnbb, n1 + 1.0) + (epsn1 / (n1 + 1.0))*pow(gnb - g_i, n1 + 1.0) +
				(epsn2 / (n2 + 1.0))*pow(g_f, n2 + 1.0) + c2 * g_f - (epsn2 / (n2 + 1.0))*pow(gnbb, n2 + 1.0) - c2 * gnbb;
		}
		else
		{
			result = (epsn2 / (n2 + 1.0))*pow(g_f, n2 + 1.0) + c2 * g_f - (epsn2 / (n2 + 1.0))*pow(g_i, n2 + 1.0) - c2 * g_i;
		}
	}
	return result* multiplier;
}
double Interface_1::IntegralTrapElasticForce(double gap_i, double gap_f)
{
	return 0.5*(EvaluateElasticForce(gap_i) + EvaluateElasticForce(gap_f))*(gap_f - gap_i);	
	//return (1.0/6.0)*(EvaluateElasticForce(gap_i) + EvaluateElasticForce(gap_f)+4.0*EvaluateElasticForce((gap_i+gap_f)/2.0))*(gap_f - gap_i);
}
