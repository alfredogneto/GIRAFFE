#pragma once
#include "Material.h"
class ElasticPlasticIsoHardening:
	public Material
{
public:
	ElasticPlasticIsoHardening();
	~ElasticPlasticIsoHardening();

	double E, H, nu, rho, sigma_y_0;					//material properties

	bool Read(FILE *f);
	void Write(FILE *f);

};

