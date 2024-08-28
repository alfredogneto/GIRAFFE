#include "ModalParameters.h"
#include <string>

ModalParameters::ModalParameters()
{
	sample = 1;
	number_modes = 1;
	number_frames = 12;
	tolerance = 0.0;
	export_matrices = true;
	compute_eigenvectors = true;
	concomitant_modal = false;
}

ModalParameters::~ModalParameters()
{
}
bool ModalParameters::Read(FILE *f)
{
	char s[1000];
	fscanf(f, "%s", s);
	if (!strcmp(s, "ExportMatrices"))
	{
		fscanf(f, "%s", s);
		if (atoi(s) != 0 && atoi(s) != 1)//Se não for 0 ou 1
			return false;
		if (atoi(s) == 1)
			export_matrices = true;
		else
			export_matrices = false;
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "NumberModes"))
	{
		fscanf(f, "%s", s);
		number_modes = atoi(s);
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "Tolerance"))
	{
		fscanf(f, "%s", s);
		tolerance = atof(s);
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "ComputeEigenvectors"))
	{
		fscanf(f, "%s", s);
		if (atoi(s) != 0 && atoi(s) != 1)//Se não for 0 ou 1
			return false;
		if (atoi(s) == 1)
			compute_eigenvectors = true;
		else
			compute_eigenvectors = false;
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "NumberFrames"))
	{
		fscanf(f, "%s", s);
		number_frames = atoi(s);
	}
	else
		return false;
	
	//Palavras chaves opcionais (analise modal concomitante)
	
	fpos_t pos;
	fgetpos(f, &pos);
	fscanf(f, "%s", s);
	if (!strcmp(s, "ConcomitantModal"))
	{
		fscanf(f, "%s", s);
		if (atoi(s) != 0 && atoi(s) != 1)//Se não for 0 ou 1
			return false;
		if (atoi(s) == 1)
			concomitant_modal = true;
		else
			concomitant_modal = false;
	}
	else
		fsetpos(f, &pos);
	
	fgetpos(f, &pos);
	fscanf(f, "%s", s);
	if (!strcmp(s, "Sample"))
	{
		fscanf(f, "%s", s);
		sample = atoi(s);
	}
	else
		fsetpos(f, &pos);

	return true;
}
void ModalParameters::Write(FILE *f)
{
	fprintf(f, "ExportMatrices\t%d\tNumberModes\t%d\tTolerance\t%.6e\tComputeEigenvectors\t%dNumberFrames%dConcomitantModal\t%d\tSample\t%d\n",
		export_matrices, number_modes, tolerance, compute_eigenvectors, number_frames, (int)concomitant_modal, sample);
}
