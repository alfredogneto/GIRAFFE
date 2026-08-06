#pragma once
#include <stdio.h>

class SolidSection;

class PipeSection
{
public:
	PipeSection();
	~PipeSection();

	int number;										//ID of the section
	double EA;										//EA -  stiffness axial
	double EI;										//EI -  stiffness flexional 					
	double GJ;										//GJ -  stiffness the torção
	double GA;										//GA -  stiffness to the cisalhamento
	double Rho;										//Mass by unit of length
	double CDt;										//Coefficient of arrasto in the direction tangential
	double CDn;										//Coefficient of arrasto in the normal direction
	double CAt;										//Coefficient of mass adicional in the direction tangential
	double CAn;										//Coefficient of mass adicional in the normal direction
	double De;										//Diametro external
	double Di;										//Diametro internal

	bool Read(FILE *f);
	void Write(FILE *f);

	SolidSection *sec_details;
	void PreCalc();
};

