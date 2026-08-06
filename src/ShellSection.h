#pragma once
#include <stdio.h>

class ShellSection
{
public:
	ShellSection() {}
	virtual ~ShellSection() {}

	virtual bool Read(FILE *f) = 0;
	virtual void Write(FILE *f) = 0;

	int number;
	double thickness;
};

