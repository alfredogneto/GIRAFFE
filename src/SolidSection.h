#pragma once
#include "SectionDetails.h"

class SolidSection :
	public SectionDetails
{
public:
	SolidSection();
	~SolidSection();
	SolidSection(SolidSection &copied);
	void Alloc(int e_points);											//Performs allocation of the tables envolvidas
	bool Read(FILE *f);
	void Write(FILE *f);
	void WriteVTK_XMLRender(FILE *f, Beam_1* elem);

	int n_points;
	double** points;													//Contem coordinates of the points that formam the ST
	
};


