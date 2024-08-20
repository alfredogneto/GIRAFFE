#pragma once
#include "Material.h"
class Orthotropic: 
	public Material
{
public:
	Orthotropic();
	~Orthotropic();
	double E1, E2, G12, G23, nu12;	//Propriedades do material ortotropico

	bool Read(FILE *f);
	void Write(FILE *f);
};

