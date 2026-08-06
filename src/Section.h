#pragma once
#include <stdio.h>
#include "Matrix.h"

class SectionDetails;

class Section
{
public:
	Section() {}
	virtual ~Section() {}

	int number;							//ID of the section
	double A, I11, I22, I12, I33, It;	//propriedades of the section
	//A - area of the ST  [L^2]
	//I11 - moment of inertia in torno of the axis 1 [L^4]
	//I22 - moment of inertia in torno of the axis 2 [L^4]
	//I12 - product of inertia in relação to the eixos 1 and 2 [L^4]
	//I33 - moment of inertia in relação to the axis 3 [L^4]
	//It - moment of torção  [L^4]

	int aerodynamicdataID = 0;				//ID of propriedades aerodinamicas associated the section transversal in question
	Matrix AC;								//Center aerodinamico (in the plane e1 x e2 of the ST)
	double aero_length;						//Length of reference for calculation aerodinamico (corda)
	SectionDetails *sec_details;

	virtual bool Read(FILE *f) = 0;
	virtual void Write(FILE *f) = 0;
	virtual void PreCalc() = 0;
};

