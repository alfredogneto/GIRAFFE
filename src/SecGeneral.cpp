#include "SecGeneral.h"


SecGeneral::SecGeneral()
{
	//aloca, mas não utiliza section details
	sec_details = new SolidSection();
}


SecGeneral::~SecGeneral()
{
	delete sec_details;
}

bool SecGeneral::Read(FILE *f)
{
	char s[1000];
	fscanf(f, "%s", s);
	number = atoi(s);
	fscanf(f, "%s", s);
	if (!strcmp(s, "A"))
	{
		fscanf(f, "%s", s);
		A = atof(s);
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "I11"))
	{
		fscanf(f, "%s", s);
		I11 = atof(s);
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "I22"))
	{
		fscanf(f, "%s", s);
		I22 = atof(s);
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "I12"))
	{
		fscanf(f, "%s", s);
		I12 = atof(s);
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "JT"))
	{
		fscanf(f, "%s", s);
		It = atof(s);
	}
	else
		return false;

	//Salva a posição (stream)
	fpos_t pos;
	fgetpos(f, &pos);
	fscanf(f, "%s", s);
	if (!strcmp(s, "AD"))
	{
		fscanf(f, "%s", s);
		aerodynamicdataID = atoi(s);

		fscanf(f, "%s", s);
		if (!strcmp(s, "AC"))
		{
			fscanf(f, "%s", s);
			AC(0, 0) = atof(s);
			fscanf(f, "%s", s);
			AC(1, 0) = atof(s);
		}
		else
			return false;

		fscanf(f, "%s", s);
		if (!strcmp(s, "AeroLength"))
		{
			fscanf(f, "%s", s);
			aero_length = atof(s);
		}
		else
			return false;
	}
	else
		fsetpos(f, &pos);	//volta à posição anterior

	return true;
}

void SecGeneral::Write(FILE *f)
{
	fprintf(f, "General\t%d\tA\t%.6e\tI11\t%.6e\tI22\t%.6e\tI12\t%.6e\tJT\t%.6e\tAD\t%d\n", number, A, I11, I22, I12, It, aerodynamicdataID);
}
void SecGeneral::PreCalc()
{
	I33 = I11 + I22;
}