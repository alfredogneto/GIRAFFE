#pragma once
#include "Section.h"
class SecHelicalFiber :
	public Section
{
public:
	SecHelicalFiber();
	~SecHelicalFiber();

	bool Read(FILE *f);
	void Write(FILE *f);
	void PreCalc();

	double R;		//helix radius
	double r;		//cross-section radius
	double alpha;	//helix angle
	int nt;			//number of turns
	int nh;			//number of helices

	//Fills contents of pointers to D, Mr and Jr (needed for Beam_1)
	void ComputeStiffnessMass(Matrix& D, Matrix& Mr, Matrix& Jr, int material_ID);
};

