#include "AerodynamicData.h"
#include <string>
#include "Table.h"

AerodynamicData::AerodynamicData()
{
	ref_position = Matrix(2);
	ref_length = 0.0;
	CL = NULL;
	CD = NULL;
	CM = NULL;
}


AerodynamicData::~AerodynamicData()
{
	if (CL != NULL)
		delete CL;
	if (CD != NULL)
		delete CD;
	if (CM != NULL)
		delete CM;
}

bool AerodynamicData::Read(FILE *f)
{
	char s[1000];
	fscanf(f, "%s", s);
	if (!strcmp(s, "AD"))
	{
		fscanf(f, "%s", s);
		number = atoi(s);
	}
	else
		return false;
	fscanf(f, "%s", s);
	if (!strcmp(s, "DragCoefficientData"))
	{
		int np = 0;
		fscanf(f, "%s", s);
		if (!strcmp(s, "NPoints"))
		{
			fscanf(f, "%s", s);
			np = atoi(s);
		}
		else
			return false;
		//Saves the position (stream)
		fpos_t pos;
		fgetpos(f, &pos);
		fscanf(f, "%s", s);
		if (strcmp(s, "(Alpha/CD)"))
			fsetpos(f, &pos);	//turn the position previous
		//Allocation and reading of the table
		CD = new Table(np, 1);
		CD->Read(f);
	}
	else
		return false;
	fscanf(f, "%s", s);
	if (!strcmp(s, "LiftCoefficientData"))
	{
		int np = 0;
		fscanf(f, "%s", s);
		if (!strcmp(s, "NPoints"))
		{
			fscanf(f, "%s", s);
			np = atoi(s);
		}
		else
			return false;
		//Saves the position (stream)
		fpos_t pos;
		fgetpos(f, &pos);
		fscanf(f, "%s", s);
		if (strcmp(s, "(Alpha/CL)"))
			fsetpos(f, &pos);	//turn the position previous
		//Allocation and reading of the table
		CL = new Table(np, 1);
		CL->Read(f);
	}
	else
		return false;
	fscanf(f, "%s", s);
	if (!strcmp(s, "PitchCoefficientData"))
	{
		int np = 0;
		fscanf(f, "%s", s);
		if (!strcmp(s, "NPoints"))
		{
			fscanf(f, "%s", s);
			np = atoi(s);
		}
		else
			return false;
		//Saves the position (stream)
		fpos_t pos;
		fgetpos(f, &pos);
		fscanf(f, "%s", s);
		if (strcmp(s, "(Alpha/CM)"))
			fsetpos(f, &pos);	//turn the position previous
		//Allocation and reading of the table
		CM = new Table(np, 1);
		CM->Read(f);
	}
	else
		return false;

	return true;
}
void AerodynamicData::Write(FILE *f)
{
	fprintf(f, "AD\t%d\n", number);
	fprintf(f, "DragCoefficientData\tNPoints\t%d\t(Alpha/CD)\n", CD->n_times);
	CD->Write(f);
	fprintf(f, "LiftCoefficientData\tNPoints\t%d\t(Alpha/CL)\n", CL->n_times);
	CL->Write(f);
	fprintf(f, "PitchCoefficientData\tNPoints\t%d\t(Alpha/CM)\n", CM->n_times);
	CM->Write(f);
}

