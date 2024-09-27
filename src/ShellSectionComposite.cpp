#include "ShellSectionComposite.h"
#include <string>

ShellSectionComposite::ShellSectionComposite()
{
	number = 0;
	thickness = 0;
	n_laminas = 0;
}

ShellSectionComposite::~ShellSectionComposite()
{
	delete[] db_laminas;
}

bool ShellSectionComposite::Read(FILE *f)
{
	char s[1000];
	fscanf(f, "%s", s);
	number = atoi(s);

	fscanf(f, "%s", s);
	if (!strcmp(s, "Laminas"))
	{
		fscanf(f, "%s", s);
		n_laminas = atoi(s);
	}

	db_laminas = new Lamina[n_laminas];

	for (int i = 0; i < n_laminas; i++)
	{
		fscanf(f, "%s", s);
		db_laminas[i].id = atoi(s);

		fscanf(f, "%s", s);
		db_laminas[i].thickness = atof(s);

		fscanf(f, "%s", s);
		db_laminas[i].angle = atof(s);
	}
	return true;
}

void ShellSectionComposite::Write(FILE *f)
{
	fprintf(f, "Composite\t%d\tLaminas\t%d\n", number, n_laminas);
	//Escrita da tabela
	for (int i = 0; i < n_laminas; i++)
	{
		fprintf(f, "%d\t%.6e\t%.6e\n", db_laminas[i].id, db_laminas[i].thickness, db_laminas[i].angle);
	}
}
