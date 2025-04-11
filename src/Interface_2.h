#pragma once
#include "ContactInterface.h"
#include "Interface_2.h"

//Marina
class Interface_2 :
	public ContactInterface
{
public:
	Interface_2();
	~Interface_2();

	bool Check();				//Checa inconsist�ncias no elemento para evitar erros de execu��o
	bool Read(FILE *f);
	void Write(FILE *f);
	void PreCalc();				//Pr�-c�lculo de vari�veis que � feito uma �nica vez no in�cio

	double EvaluateElasticForce(double gap);
	double EvaluateElasticStiffness(double gap);
	double IntegralElasticForce(double gap_i, double gap_f);
	double IntegralTrapElasticForce(double gap_i, double gap_f);

	//Normal contact
	double epsn;
	double cn;

	//Tangential contact
	double epst;
	double mus;
	double mud;
	double ct;
};

