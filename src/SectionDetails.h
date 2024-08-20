#pragma once
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "Beam_1.h"

class SectionDetails
{
public:
	SectionDetails();
	virtual ~SectionDetails();
	virtual bool Read(FILE *f) = 0;
	virtual void Write(FILE *f) = 0;
	virtual void WriteVTK_XMLRender(FILE *f, Beam_1* elem) = 0;
	//Variáveis para todos os SD
	char section_type[100];
	double axis_position[2];						//Coordenadas x,y do eixo. A posição em que o eixo está em relação aos pontos inseridos. Deve ser condizente com os dados calculados
	int number;
};

