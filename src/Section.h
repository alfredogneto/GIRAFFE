#pragma once
#include <stdio.h>
#include "Matrix.h"

class SectionDetails;

class Section
{
public:
	Section() {}
	virtual ~Section() {}

	int number;							//ID da seção
	double A, I11, I22, I12, I33, It;	//propriedades da seção
	//A - área da ST  [L^2]
	//I11 - momento de inércia em torno do eixo 1 [L^4]
	//I22 - momento de inércia em torno do eixo 2 [L^4]
	//I12 - produto de inércia em relação aos eixos 1 e 2 [L^4]
	//I33 - momento de inércia em relação ao eixo 3 [L^4]
	//It - momento de torção  [L^4]

	int aerodynamicdataID = 0;				//ID de propriedades aerodinâmicas associadas à seção transversal em questão
	Matrix AC;								//Centro aerodinâmico (no plano e1 x e2 da ST)
	double aero_length;						//Comprimento de referência para cálculo aerodinâmico (corda)
	SectionDetails *sec_details;

	virtual bool Read(FILE *f) = 0;
	virtual void Write(FILE *f) = 0;
	virtual void PreCalc() = 0;
};

