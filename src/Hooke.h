#pragma once
#include "Material.h"
class Hooke :
	public Material
{
public:
	Hooke();
	~Hooke();
	double E, nu;	//Propriedades do material (Hooke)

	bool Read(FILE *f);
	void Write(FILE *f);
};

