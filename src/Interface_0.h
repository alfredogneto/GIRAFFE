#pragma once
#include "ContactInterface.h"
class Interface_0 :
	public ContactInterface
{
public:
	Interface_0();
	~Interface_0();

	bool Check();				//Checks inconsistencies for avoid errors of execution
	bool Read(FILE *f);
	void Write(FILE *f);
	void PreCalc();				//Pre-calculation of variables that and done a single time in the start

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

