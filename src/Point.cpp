#include "Point.h"
#include <string>

#include"Database.h"
//Variáveis globais
extern
Database db;

Point::Point()
{
	coordinates = Matrix(3);
}

Point::~Point()
{
	
}

bool Point::Read(FILE *f)
{
	char s[1000];
	fscanf(f, "%s", s);
	//Verifica a palavra chave "Point"
	if (!strcmp(s, "Point"))
	{
		fscanf(f, "%s", s);
		number = atoi(s);
	}
	else
		return false;
	
	//Leitura das coordenadas
	fscanf(f, "%s", s);
	coordinates(0, 0) = atof(s);

	fscanf(f, "%s", s);
	coordinates(1, 0) = atof(s);

	fscanf(f, "%s", s);
	coordinates(2, 0) = atof(s);

	return true;
}

void Point::Write(FILE *f)
{					
	fprintf(f, "Point\t%d\t%.6e\t%.6e\t%.6e\n", number, coordinates(0, 0), coordinates(1, 0), coordinates(2, 0));
}
