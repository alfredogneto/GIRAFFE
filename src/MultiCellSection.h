#pragma once
#include "SectionDetails.h"

class MultiCellSection :
	public SectionDetails
{
public:
	MultiCellSection();
	~MultiCellSection();
	MultiCellSection(MultiCellSection &copied);
	void Alloc(int e_points, int e_webs);											//Realiza alocação das tabelas envolvidas
	bool Read(FILE *f);
	void Write(FILE *f);
	void WriteVTK_XMLRender(FILE *f, Beam_1* elem);
	int n_points;
	int n_webs;
	double** points;				//Contém coordenadas dos pontos que formam a ST
	int** webs;						//Contém índices que indicam numerações das webs

};

