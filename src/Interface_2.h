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

	bool Check();				//Checa inconsistências no elemento para evitar erros de execução
	bool Read(FILE *f);
	void Write(FILE *f);
	void PreCalc();				//Pré-cálculo de variáveis que é feito uma única vez no início

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

