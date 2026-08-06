#pragma once
#include "Section.h"

class SecRectangle :
	public Section
{
public:
	SecRectangle();
	~SecRectangle();

	bool Read(FILE *f);
	void Write(FILE *f);
	void PreCalc();
	double b, h;		//Base, altura
};

