#pragma once
#include "Section.h"

class SecSuperEllipse :
	public Section
{
public:
	SecSuperEllipse();
	~SecSuperEllipse();

	bool Read(FILE *f);
	void Write(FILE *f);
	void PreCalc();
	double a, b, n;		//Eixos e expoentes
	int div_mesh_a;

	double MDF_SaintVenantSE();

	double Beta(double p1, double p2);
};

