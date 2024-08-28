#pragma once
#include "Section.h"

class SecTube :
	public Section
{
public:
	SecTube();
	~SecTube();

	bool Read(FILE *f);
	void Write(FILE *f);
	void PreCalc();
	double De, Di;		//diametros
};

