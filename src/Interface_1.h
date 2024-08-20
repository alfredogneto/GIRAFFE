#pragma once
#include "ContactInterface.h"

class Interface_1 :
	public ContactInterface
{
public:
	Interface_1();
	~Interface_1();

	bool Check();				//Checa inconsistências para evitar erros de execução
	bool Read(FILE *f);
	void Write(FILE *f);
	void PreCalc();				//Pré-cálculo de variáveis que é feito uma única vez no início

	double EvaluateElasticForce(double gap);
	double EvaluateElasticStiffness(double gap);
	double IntegralElasticForce(double gap_i, double gap_f);
	double IntegralTrapElasticForce(double gap_i, double gap_f);

	//Normal contact
	double epsn1;
	double n1;
	double gnb;
	double gnbb;
	double factor;
	double n2;
	double zetan;

	double epsn2;
	double c2;

	//Tangential contact
	double epst;
	double mus;
	double mud;
	double ct;
};

