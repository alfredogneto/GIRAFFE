#pragma once
#include "SectionDetails.h"
class SolidSection :
	public SectionDetails
{
public:
	SolidSection();
	~SolidSection();
	SolidSection(SolidSection &copied);
	void Alloc(int e_points);											//Realiza alocação das tabelas envolvidas
	bool Read(FILE *f);
	void Write(FILE *f);
	void WriteVTK_XMLRender(FILE *f, Beam_1* elem);

	int n_points;
	double** points;													//Contém coordenadas dos pontos que formam a ST
	
};


