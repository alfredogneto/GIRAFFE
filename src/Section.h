#pragma once
#include <stdio.h>
#include "Matrix.h"

class SectionDetails;

class Section
{
public:
	Section() {}
	virtual ~Section() {}

	int number;							//ID da se��o
	double A, I11, I22, I12, I33, It;	//propriedades da se��o
	//A - �rea da ST  [L^2]
	//I11 - momento de in�rcia em torno do eixo 1 [L^4]
	//I22 - momento de in�rcia em torno do eixo 2 [L^4]
	//I12 - produto de in�rcia em rela��o aos eixos 1 e 2 [L^4]
	//I33 - momento de in�rcia em rela��o ao eixo 3 [L^4]
	//It - momento de tor��o  [L^4]

	int aerodynamicdataID = 0;				//ID de propriedades aerodin�micas associadas � se��o transversal em quest�o
	Matrix AC;								//Centro aerodin�mico (no plano e1 x e2 da ST)
	double aero_length;						//Comprimento de refer�ncia para c�lculo aerodin�mico (corda)
	SectionDetails *sec_details;

	virtual bool Read(FILE *f) = 0;
	virtual void Write(FILE *f) = 0;
	virtual void PreCalc() = 0;
};

