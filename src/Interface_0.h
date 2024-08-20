#pragma once
#include "ContactInterface.h"
class Interface_0 :
	public ContactInterface
{
public:
	Interface_0();
	~Interface_0();

	bool Check();				//Checa inconsist�ncias para evitar erros de execu��o
	bool Read(FILE *f);
	void Write(FILE *f);
	void PreCalc();				//Pr�-c�lculo de vari�veis que � feito uma �nica vez no in�cio

	//Normal contact
	double epsn1;
	double n1;
	double zetan;

	//Tangential contact
	double epst;
	double mus;
	double mud;
	double ct;

	//Not used (but created)
	double n2;
	double gnb;
	double gnbb;
};

