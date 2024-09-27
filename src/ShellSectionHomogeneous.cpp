#include "ShellSectionHomogeneous.h"
#include <string>

ShellSectionHomogeneous::ShellSectionHomogeneous()
{
	thickness = 0;
	number = 0;
}


ShellSectionHomogeneous::~ShellSectionHomogeneous()
{
}

bool ShellSectionHomogeneous::Read(FILE *f)
{
	char s[1000];
	fscanf(f, "%s", s);
	number = atoi(s);
	fscanf(f, "%s", s);
	if (!strcmp(s, "Thickness"))
	{
		fscanf(f, "%s", s);
		thickness = atof(s);
	}
	else
		return false;
	return true;
}

void ShellSectionHomogeneous::Write(FILE *f)
{
	fprintf(f, "Homogeneous\t%d\tThickness\t%.6e\n", number, thickness);
}
