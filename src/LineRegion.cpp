#include "LineRegion.h"
#include <string>

LineRegion::LineRegion()
{
	number = 0;
	n_elements = 0;
	elements = NULL;
}


LineRegion::~LineRegion()
{
	if (elements != NULL)
		delete []elements;
}

bool LineRegion::Read(FILE *f)
{
	char s[1000];
	fscanf(f, "%s", s);
	if (!strcmp(s, "LR"))
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

void LineRegion::Write(FILE *f)
{
	fprintf(f, "LR\t%d\tNElements\t%d\tElements",
		number, n_elements);
	for (int i = 0; i < n_elements;i++)
		fprintf(f, "\t%d",elements[i]);
	fprintf(f, "\n");
}
