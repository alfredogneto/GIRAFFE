#pragma once
#include <stdio.h>


class ContactInterface
{
public:
	ContactInterface() {}
	virtual ~ContactInterface() {}

	virtual bool Check() = 0;				//Checa inconsistências para evitar erros de execução
	virtual bool Read(FILE *f) = 0;
	virtual void Write(FILE *f) = 0;
	virtual void PreCalc() = 0;				//Pre-calculo de variaveis que e feito uma unica vez no inicio


	int number;			//ID
	int material_1;		//Material 1 ID
	int material_2;		//Material 2 ID
};

