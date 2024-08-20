#pragma once
#include "Section.h"
class SecGeneral :
	public Section
{
public:
	SecGeneral();
	~SecGeneral();

	bool Read(FILE *f);
	void Write(FILE *f);
	void PreCalc();
};

