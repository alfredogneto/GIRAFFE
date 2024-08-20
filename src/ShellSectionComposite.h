#pragma once
#include "ShellSection.h"

struct Lamina
{
	int id;
	double thickness;
	double angle;
};

class ShellSectionComposite:
	public ShellSection
{
public:
	ShellSectionComposite();
	~ShellSectionComposite();
	bool Read(FILE *f);
	void Write(FILE *f);

	//Composite Variables
	int n_laminas;	
	Lamina* db_laminas;
};

