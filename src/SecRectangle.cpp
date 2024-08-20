#include "SecRectangle.h"


SecRectangle::SecRectangle()
{
	sec_details = new SolidSection();

	AC = Matrix(2);
	aero_length = 0;
}

SecRectangle::~SecRectangle()
{
	delete sec_details;
}

bool SecRectangle::Read(FILE *f)
{
	char s[1000];
	fscanf(f, "%s", s);
	number = atoi(s);
	fscanf(f, "%s", s);
	if (!strcmp(s, "B"))
	{
		fscanf(f, "%s", s);
		b = atof(s);
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "H"))
	{
		fscanf(f, "%s", s);
		h = atof(s);
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

void SecRectangle::Write(FILE *f)
{
	fprintf(f, "Rectangle\t%d\tB\t%.6e\tH\t%.6e\n", number, b,h);
}

void SecRectangle::PreCalc()
{
	SolidSection* ptr_sd = static_cast<SolidSection*>(sec_details);
	A = b*h;
	I11 = b*h*h*h / 12.0;
	I22 = h*b*b*b / 12.0;
	I12 = 0.0;
	I33 = I11 + I22;
	//Momento de torção - Timoshenko - solução da torção de Saint-Venant - Theory of Elasticity - pág. ~278 eq 160
	double temp = 0;
	for (int n = 1; n < 22; n = n + 2)
	{ 
		temp += (1.0 / (pow((double)n,5)))*tanh(n*PI*h / (2 * b));
	}
	It = (1.0 / 3.0)*b*b*b*h*(1.0 - 192.0*b*temp / (pow(PI,5)*h));

	ptr_sd->Alloc(4);
	ptr_sd->points[0][0] = -b / 2;
	ptr_sd->points[0][1] = +h / 2;
	ptr_sd->points[1][0] = +b / 2;
	ptr_sd->points[1][1] = +h / 2;
	ptr_sd->points[2][0] = +b / 2;
	ptr_sd->points[2][1] = -h / 2;
	ptr_sd->points[3][0] = -b / 2;
	ptr_sd->points[3][1] = -h / 2;
}