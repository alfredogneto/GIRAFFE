#pragma once
#include "SectionDetails.h"

class MultiCellSection :
	public SectionDetails
{
public:
	MultiCellSection();
	~MultiCellSection();
	MultiCellSection(MultiCellSection &copied);
	void Alloc(int e_points, int e_webs);											//Performs allocation of the tables envolvidas
	bool Read(FILE *f);
	void Write(FILE *f);
	void WriteVTK_XMLRender(FILE *f, Beam_1* elem);
	int n_points;
	int n_webs;
	double** points;				//Contem coordinates of the points that formam the ST
	int** webs;						//Contem indices that indicam numerações of the webs

};

