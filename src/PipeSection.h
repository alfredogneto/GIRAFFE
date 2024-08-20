#pragma once
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "SolidSection.h"
#include <math.h>
#define PI 3.1415926535897932384626433832795

class PipeSection
{
public:
	PipeSection();
	~PipeSection();

	int number;										//ID da seção
	double EA;										//EA -  rigidez axial
	double EI;										//EI -  rigidez flexional 					
	double GJ;										//GJ -  rigidez à torção
	double GA;										//GA -  rigidez ao cisalhamento
	double Rho;										//Massa por unidade de comprimento
	double CDt;										//Coeficiente de arrasto na direção tangencial
	double CDn;										//Coeficiente de arrasto na direção normal
	double CAt;										//Coeficiente de massa adicional na direção tangencial
	double CAn;										//Coeficiente de massa adicional na direção normal
	double De;										//Diâmetro externo
	double Di;										//Diâmetro interno

	bool Read(FILE *f);
	void Write(FILE *f);

	SolidSection *sec_details;
	void PreCalc();
};

