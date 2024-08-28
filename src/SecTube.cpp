#include "SecTube.h"
#include <string>

#include "SolidSection.h"

#define PI 3.1415926535897932384626433832795

SecTube::SecTube()
{
	sec_details = new SolidSection();

	AC = Matrix(2);
	aero_length = 0;
}


SecTube::~SecTube()
{
	delete sec_details;
}

bool SecTube::Read(FILE *f)
{
	char s[1000];
	fscanf(f, "%s", s);
	number = atoi(s);
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

void SecTube::Write(FILE *f)
{
	fprintf(f, "Tube\t%d\tDe\t%.6e\tDi\t%.6e\tAD\t%d\n", number, De, Di, aerodynamicdataID);
}

void SecTube::PreCalc()
{
	SolidSection* ptr_sd = static_cast<SolidSection*>(sec_details);
	A = (PI/4.0)*(De*De-Di*Di);
	I11 = (PI / 64.0)*(De*De*De*De - Di*Di*Di*Di);
	I22 = (PI / 64.0)*(De*De*De*De - Di*Di*Di*Di);
	I12 = 0.0;
	I33 = (PI / 32.0)*(De*De*De*De - Di*Di*Di*Di);
	It = I33;

	int n_circ = 24;
	ptr_sd->Alloc(n_circ);
	double r = De/2;
	double theta = 0;
	for (int index = 0; index < n_circ; index++)
	{
		theta = (index * 2 * PI) / n_circ;
		ptr_sd->points[index][0] = r*cos(theta);
		ptr_sd->points[index][1] = r*sin(theta);
	}
}

