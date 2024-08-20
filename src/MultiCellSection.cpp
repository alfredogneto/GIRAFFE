#include "MultiCellSection.h"

MultiCellSection::MultiCellSection()
{
	n_points = 0;;
	n_webs = 0;
	points = NULL;
	webs = NULL;
	number = 0;				//ID
	axis_position[0] = 0.0;	//Ref position
	axis_position[1] = 0.0;	//Ref position
	sprintf(section_type, "MultiCellSection");
}
//Construtor de cópia
MultiCellSection::MultiCellSection(MultiCellSection &copied)
{
	n_points = 0;;
	n_webs = 0;
	points = NULL;
	webs = NULL;
	number = 0;				//ID
	axis_position[0] = 0.0;	//Ref position
	axis_position[1] = 0.0;	//Ref position
	sprintf(section_type, "MultiCellSection");

	Alloc(copied.n_points, copied.n_webs);
	//Copia os valores dos pontos
	for (long i = 0; i < copied.n_points; i++)
	{
		points[i][0] = copied.points[i][0];
		points[i][1] = copied.points[i][1];
	}
	//Copia os valores das webs
	for (long i = 0; i < copied.n_webs; i++)
	{
		webs[i][0] = copied.webs[i][0];
		webs[i][1] = copied.webs[i][1];
	}
	strcpy(section_type, copied.section_type);
	axis_position[0] = copied.axis_position[0];
	axis_position[1] = copied.axis_position[1];
}

void MultiCellSection::Alloc(int e_points, int e_webs)
{
	//Desaloca matrizes
	if (points != NULL)
	{
		for (int i = 0; i < n_points; i++)
			delete points[i];
		delete[]points;
	}
	if (webs != NULL)
	{
		for (int i = 0; i < n_webs; i++)
			delete webs[i];
		delete[]webs;
	}
	n_points = e_points;
	n_webs = e_webs;
	//Aloca matrizes
	points = new double*[n_points];
	for (int i = 0; i < n_points; i++)
		points[i] = new double[2];
	webs = new int*[n_webs];
	for (int i = 0; i < n_webs; i++)
		webs[i] = new int[2];
}

MultiCellSection::~MultiCellSection()
{
	//Desaloca matrizes
	if (points != NULL)
	{
		for (int i = 0; i < n_points; i++)
			delete points[i];
		delete[]points;
	}
	if (webs != NULL)
	{
		for (int i = 0; i < n_webs; i++)
			delete webs[i];
		delete[]webs;
	}
}
bool MultiCellSection::Read(FILE *f)
{
	char s[1000];
	fscanf(f, "%s", s);
	number = atoi(s);
	fscanf(f, "%s", s);
	if (!strcmp(s, "AxisPosition"))
	{
		fscanf(f, "%s", s);
		axis_position[0] = atof(s);
		fscanf(f, "%s", s);
		axis_position[1] = atof(s);
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "NPoints"))
	{
		fscanf(f, "%s", s);
		n_points = atoi(s);
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "NWebs"))
	{
		fscanf(f, "%s", s);
		n_webs = atoi(s);
	}
	else
		return false;

	Alloc(n_points, n_webs);
	//Leitura dos pontos
	double x, y;
	int index;
	for (int p = 0; p < n_points; p++)
	{
		fscanf(f, "%s", s);
		if (!strcmp(s, "Point"))
		{
			fscanf(f, "%s", s);
			index = atoi(s);
			fscanf(f, "%s", s);
			x = atof(s);
			fscanf(f, "%s", s);
			y = atof(s);
			points[p][0] = x;
			points[p][1] = y;
		}
		else
			return false;
	}
	//Leitura dos webs
	int p1, p2;
	for (int web = 0; web < n_webs; web++)
	{
		fscanf(f, "%s", s);
		if (!strcmp(s, "Web"))
		{
			fscanf(f, "%s", s);
			index = atoi(s);
			fscanf(f, "%s", s);
			p1 = atoi(s);
			fscanf(f, "%s", s);
			p2 = atoi(s);
			webs[web][0] = p1;
			webs[web][1] = p2;
		}
		else
			return false;
	}
	return true;
}
void MultiCellSection::Write(FILE *f)
{
	fprintf(f, section_type);
	fprintf(f, "\t%d\t", number);
	fprintf(f, "\tAxisPosition\t%.6f\t%.6f\t", axis_position[0], axis_position[1]);
	fprintf(f, "NPoints\t%d\tNWebs\t%d\n", n_points, n_webs);
	for (int p = 0; p < n_points; p++)
	{
		fprintf(f, "Point\t%d\t%.6f\t%.6f\n", p + 1, points[p][0], points[p][1]);
	}
	for (int web = 0; web < n_webs; web++)
	{
		fprintf(f, "Web\t%d\t%d\t%d\n", web + 1, webs[web][0], webs[web][1]);
	}
}

void MultiCellSection::WriteVTK_XMLRender(FILE *f, Beam_1* elem)
{
	//DOES NOTHING
}