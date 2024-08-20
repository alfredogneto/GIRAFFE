#pragma once
#include "Section.h"
class SecUserDefined :
	public Section
{
public:
	SecUserDefined();
	~SecUserDefined();

	bool Read(FILE *f);
	void Write(FILE *f);
	void PreCalc();

	//Variables
	double GA;
	double EA;
	double ES1;
	double ES2;
	double EI11;
	double EI22;
	double EI12;
	double GS1;
	double GS2;
	double GS1S;
	double GS2S;
	double GJT;
	double J11;
	double J22;
	double J12;
	
	Matrix SC;
	Matrix BC;
	double Rho;
	int sec_details_ID;

};