#pragma once
#include "SectionDetails.h"

class MultiCellSection :
	public SectionDetails
{
public:
	MultiCellSection();
	~MultiCellSection();
	MultiCellSection(MultiCellSection &copied);
	void Alloc(int e_points, int e_webs);											//Realiza aloca��o das tabelas envolvidas
	bool Read(FILE *f);
	void Write(FILE *f);
	void WriteVTK_XMLRender(FILE *f, Beam_1* elem);
	int n_points;
	int n_webs;
	double** points;				//Cont�m coordenadas dos pontos que formam a ST
	int** webs;						//Cont�m �ndices que indicam numera��es das webs

};

