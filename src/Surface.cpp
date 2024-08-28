#include "Surface.h"
#include <string>

#include "Matrix.h"
#include "NSContactData.h"

Surface::Surface()
{
}


Surface::~Surface()
{
}

bool Surface::ReadCommon(FILE *f)
{
	char s[1000];

	//Salva a posição (stream)
	fpos_t pos;
	fgetpos(f, &pos);
	fscanf(f, "%s", s);
	if (!strcmp(s, "Degeneration"))
	{
		//Salva a posição (stream)
		fgetpos(f, &pos);
		fscanf(f, "%s", s);
		if (!strcmp(s, "Coord1"))
		{
			fscanf(f, "%s", s);
			degeneration[0] = true;
			entered_u1 = true;
			deg_coordinates[0][0][0] = atof(s);
		}
		else
			fsetpos(f, &pos);

		//Salva a posição (stream)
		fgetpos(f, &pos);
		fscanf(f, "%s", s);
		if (!strcmp(s, "Coord2"))
		{
			fscanf(f, "%s", s);
			degeneration[1] = true;
			entered_u2 = true;
			deg_coordinates[0][0][1] = atof(s);
		}
		else
			fsetpos(f, &pos);

		//Salva a posição (stream)
		fgetpos(f, &pos);
		fscanf(f, "%s", s);
		if (!strcmp(s, "Div1"))
		{
			fscanf(f, "%s", s);
			degeneration[0] = true;
			div1 = atoi(s);
		}
		else
			fsetpos(f, &pos);
		//Salva a posição (stream)
		fgetpos(f, &pos);
		fscanf(f, "%s", s);
		if (!strcmp(s, "Div2"))
		{
			fscanf(f, "%s", s);
			degeneration[1] = true;
			div2 = atoi(s);
		}
		else
			fsetpos(f, &pos);
	}
	else
		fsetpos(f, &pos);

	return true;
}

void Surface::AllocDegeneration()
{
	if (alloced_degeneration == false)
	{
		deg_coordinates = new double**[div1];
		for (int i = 0; i < div1; i++)
		{
			deg_coordinates[i] = new double*[div2];
			for (int j = 0; j < div2; j++)
			{
				deg_coordinates[i][j] = new double[2];
			}
		}
			
		for (int i = 0; i < div1; i++)
		{
			for (int j = 0; j < div2; j++)
			{
				deg_coordinates[i][j][0] = 0.0;
				deg_coordinates[i][j][1] = 0.0;
			}
		}

		alloced_div1 = div1;
		alloced_div2 = div2;
		alloced_degeneration = true;
	}
}

void Surface::FreeDegeneration()
{
	if (alloced_degeneration == true)
	{
		for (int i = 0; i < alloced_div1; i++)
		{
			for (int j = 0; j < alloced_div2; j++)
			{
				delete[] deg_coordinates[i][j];
			}
			delete[] deg_coordinates[i];
		}
		delete[] deg_coordinates;
		deg_coordinates = NULL;
		alloced_degeneration = false;
	}
}

void Surface::InitializeDegeneration()
{
	//Degeneration
	degeneration[0] = false;
	degeneration[1] = false;
	div1 = 1;
	div2 = 1;
	
	deg_coordinates = NULL;
	entered_u1 = false;
	entered_u2 = false;
	alloced_degeneration = false;
	alloced_div1 = 0;
	alloced_div2 = 0;
	AllocDegeneration();
}

void Surface::DegenerationPreCalc()
{
	//Alocação de acordo com nível de degeneração
	double copy_u1, copy_u2;
	if (entered_u1)
	{
		copy_u1 = deg_coordinates[0][0][0];
		div1 = 1;
	}
		
	if (entered_u2)
	{
		copy_u2 = deg_coordinates[0][0][1];
		div2 = 1;
	}
		
	FreeDegeneration();
	AllocDegeneration();
	SetMinMaxRange();

	//Preenche valores de coordenadas degeneradas
	double step_u1 = (u1_max - u1_min) / div1;
	double step_u2 = (u2_max - u2_min) / div2;
	double cur_u1 = u1_min + 0.5 * step_u1;
	double cur_u2 = u2_min + 0.5 * step_u2;
	
	for (int i = 0; i < div1; i++)
	{
		cur_u2 = u2_min + 0.5 * step_u2;
		for (int j = 0; j < div2; j++)
		{
			deg_coordinates[i][j][0] = cur_u1;
			deg_coordinates[i][j][1] = cur_u2;
			cur_u2 += step_u2;
		}
		cur_u1 += step_u1;
	}
	//Correção em caso de entrada manual de coordenadas
	if (entered_u1)
	{
		for (int j = 0; j < div2; j++)
			deg_coordinates[0][j][0] = copy_u1;
	}
	if (entered_u2)
	{
		for (int i = 0; i < div1; i++)
			deg_coordinates[i][0][1] = copy_u2;
	}
}
