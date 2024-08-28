#include "CoordinateSystem.h"

#include <string>

#include "Matrix.h"

CoordinateSystem::CoordinateSystem()
{
	E1 = new Matrix(3);
	E2 = new Matrix(3);
	E3 = new Matrix(3);
	Q = new Matrix(3, 3);
}

CoordinateSystem::~CoordinateSystem()
{
	delete E1;
	delete E2;
	delete E3;
	delete Q;
}

bool CoordinateSystem::Read(FILE *f)
{
	char s[1000];
	fscanf(f, "%s", s);
	//Verifica a palavra chave "CS"
	if (!strcmp(s, "CS"))
	{
		fscanf(f, "%s", s);
		number = atoi(s);
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "E1"))
	{
		fscanf(f, "%s", s);
		(*E1)(0, 0) = atof(s);
		fscanf(f, "%s", s);
		(*E1)(1, 0) = atof(s);
		fscanf(f, "%s", s);
		(*E1)(2, 0) = atof(s);
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "E3"))
	{
		fscanf(f, "%s", s);
		(*E3)(0, 0) = atof(s);
		fscanf(f, "%s", s);
		(*E3)(1, 0) = atof(s);
		fscanf(f, "%s", s);
		(*E3)(2, 0) = atof(s);
	}
	else
		return false;
	double tolortho = 1e-8;
	if (abs(dot(*E3, *E1)) >= tolortho)
	{
		printf("Non-orthogonal vectors set to E1 and E3 from Coordinate System %d.", number);
		return false;
	}

	*E2 = cross(*E3, *E1);

	//Tratamento aos versores - verificação e normalização unitaria
	double norm1 = norm(*E1);
	double norm2 = norm(*E2);
	double norm3 = norm(*E3);
	if (norm1 != 1.0)
		(*E1) = (1.0 / norm1)*(*E1);
	if (norm2 != 1.0)
		(*E2) = (1.0 / norm2)*(*E2);
	if (norm3 != 1.0)
		(*E3) = (1.0 / norm3)*(*E3);

	//Montagem da matrix de transformação de coordenadas
	(*Q)(0, 0) = (*E1)(0, 0);
	(*Q)(0, 1) = (*E1)(1, 0);
	(*Q)(0, 2) = (*E1)(2, 0);

	(*Q)(1, 0) = (*E2)(0, 0);
	(*Q)(1, 1) = (*E2)(1, 0);
	(*Q)(1, 2) = (*E2)(2, 0);

	(*Q)(2, 0) = (*E3)(0, 0);
	(*Q)(2, 1) = (*E3)(1, 0);
	(*Q)(2, 2) = (*E3)(2, 0);


	return true;
}

void CoordinateSystem::Write(FILE *f)
{
	fprintf(f, "CS\t%d\tE1\t%.6e\t%.6e\t%.6e\tE3\t%.6e\t%.6e\t%.6e\n", number, (*E1)(0, 0), (*E1)(1, 0), (*E1)(2, 0), (*E3)(0, 0), (*E3)(1, 0), (*E3)(2, 0));
}
