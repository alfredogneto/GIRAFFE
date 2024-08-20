#include "SurfaceRegion.h"


SurfaceRegion::SurfaceRegion()
{
	number = 0;
	n_elements = 0;
	elements = NULL;
}


SurfaceRegion::~SurfaceRegion()
{
	if (elements != NULL)
		delete[]elements;
}

bool SurfaceRegion::Read(FILE *f)
{
	char s[1000];
	fscanf(f, "%s", s);
	if (!strcmp(s, "SR"))
	{
		fscanf(f, "%s", s);
		number = atoi(s);
		fscanf(f, "%s", s);
		if (!strcmp(s, "NElements"))
		{
			fscanf(f, "%s", s);
			n_elements = atoi(s);

			//Alocação do vetor de elementos
			elements = new int[n_elements];
		}
		else
			return false;

		fscanf(f, "%s", s);
		if (!strcmp(s, "Elements"))
		{
			for (int i = 0; i < n_elements; i++)
			{
				fscanf(f, "%s", s);
				elements[i] = atoi(s);
			}
		}
		else
			return false;
	}
	else
		return false;
	return true;
}

void SurfaceRegion::Write(FILE *f)
{
	fprintf(f, "SR\t%d\tNElements\t%d\tElements",
		number, n_elements);
	for (int i = 0; i < n_elements; i++)
		fprintf(f, "\t%d", elements[i]);
	fprintf(f, "\n");
}