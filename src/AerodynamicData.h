#pragma once
#include "Matrix.h"
class Table;

//This class contains a structures to save aerodynamic data
class AerodynamicData
{
public:
	AerodynamicData();
	~AerodynamicData();

	int number;				//ID
	Matrix ref_position;	//Reference position for aerodynamic forces
	double ref_length;		//Comprimento de refer�ncia (corda)
	Table* CL;				//Lift coeff. vs. alpha
	Table* CD;				//Drag coeff. vs. alpha
	Table* CM;				//Pitching coeff. vs. alpha

	bool Read(FILE *f);
	void Write(FILE *f);
};

