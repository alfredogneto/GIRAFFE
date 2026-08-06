#pragma once
#include <stdio.h>

class Beam_1;

class SectionDetails
{
public:
	SectionDetails() {}
	virtual ~SectionDetails() {}
	virtual bool Read(FILE *f) = 0;
	virtual void Write(FILE *f) = 0;
	virtual void WriteVTK_XMLRender(FILE *f, Beam_1* elem) = 0;
	//Variables for all the SD
	char section_type[100];
	double axis_position[2];						//Coordinates x,y of the axis. A position in that the axis this in relação to the points inseridos. Must be condizente with the data calculated
	int number;
};

