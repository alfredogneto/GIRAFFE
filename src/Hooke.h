#pragma once
#include "Material.h"
class Hooke :
	public Material
{
public:
	Hooke();
	~Hooke();
	double E, nu;	//Propriedades of the material (Hooke)

	bool Read(FILE *f);
	void Write(FILE *f);
};

