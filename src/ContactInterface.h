#pragma once
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "Matrix.h"
#include <vector>

class ContactInterface
{
public:
	ContactInterface();
	virtual ~ContactInterface();

	virtual bool Check() = 0;				//Checa inconsist�ncias para evitar erros de execu��o
	virtual bool Read(FILE *f) = 0;
	virtual void Write(FILE *f) = 0;
	virtual void PreCalc() = 0;				//Pr�-c�lculo de vari�veis que � feito uma �nica vez no in�cio


	int number;			//ID
	int material_1;		//Material 1 ID
	int material_2;		//Material 2 ID
};

