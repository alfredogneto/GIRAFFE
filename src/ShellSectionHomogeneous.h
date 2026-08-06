#pragma once
#include "ShellSection.h"

class ShellSectionHomogeneous:
	public ShellSection
{
public:
	ShellSectionHomogeneous();
	~ShellSectionHomogeneous();
	bool Read(FILE *f);
	void Write(FILE *f);
};

