#include "ConcomitantSolution.h"
#include <direct.h>

#include "Modal.h"

#include "Database.h"
//Variáveis globais
extern
Database db;

ConcomitantSolution::ConcomitantSolution()
{
	bool_concomitant.SetDefault(true);
	sample = 1;
	sol = NULL;
	f_output = NULL;
}


ConcomitantSolution::~ConcomitantSolution()
{
	if (sol != NULL)
		delete sol;
}

//Leitura
bool ConcomitantSolution::Read(FILE *f)
{
	char s[1000];
	fscanf(f, "%s", s);
	if (!strcmp(s, "Sample"))
	{
		fscanf(f, "%s", s);
		sample = atoi(s);
	}
	else
		return false;
	fscanf(f, "%s", s);
	if (!strcmp(s, "BoolTable"))
		bool_concomitant.Read(f);
	else
		return false;
	//Reading concomitant analysis data
	fscanf(f, "%s", s);
	if (!strcmp(s, "Modal"))
		sol = new Modal();
	else
		return false;
	Modal* ptr = static_cast<Modal*>(sol);
	fscanf(f, "%s", s);
	if (!strcmp(s, "NumberModes"))
	{
		fscanf(f, "%s", s);
		ptr->number_modes = atoi(s);
	}
	else
		return false;
	fscanf(f, "%s", s);
	if (!strcmp(s, "Tolerance"))
	{
		fscanf(f, "%s", s);
		ptr->tolerance = atof(s);
	}
	else
		return false;
	//Setting aditional parameters
	ptr->compute_eigenvectors = false;
	ptr->export_matrices = false;
	ptr->concomitant_solution = true;
	return true;
}

//Gravação
void ConcomitantSolution::Write(FILE *f)
{
	Modal* ptr = static_cast<Modal*>(sol);
	fprintf(f, "\nConcomitantSolution\tSample\t%d\t", sample);
	bool_concomitant.Write(f);
	fprintf(f, "Modal\tNumberModes\t%d\tTolerance\t%.6e\n", ptr->number_modes, ptr->tolerance);
}

//Checking inconsistencies
bool ConcomitantSolution::Check()
{
	return true;
}
void ConcomitantSolution::StartConcomitantSolution()
{
	Modal* ptr = static_cast<Modal*>(sol);
	//Abre os arquivos dos nós
	char name[200];
	strcpy(name, db.folder_name);
	strcat(name, "post/");
	_mkdir(name);
	strcat(name, "concomitant_solution/");
	_mkdir(name);
	strcat(name, "concomitant_solution.txt");
	f_output = fopen(name, "w");
	//Cabeçalho
	fprintf(f_output, "TIME\t");
	for (int i = 0; i < ptr->number_modes; i++)
		fprintf(f_output, "Re_%d\tIm_%d\t", i + 1, i + 1);
	fprintf(f_output, "\n");
}

void ConcomitantSolution::UpdateConcomitantSolution(double time)
{
	//Somente atualiza se for ativo para a atual solution
	if (bool_concomitant.GetAt(db.current_solution_number - 1))
	{
		Modal* ptr = static_cast<Modal*>(sol);
		bool return_variable = ptr->Solve();
		if (return_variable == true)
		{
			//Writing results
			fprintf(f_output, "%.6f\t", time);
			for (int mode = 0; mode < ptr->number_modes; mode++)
				fprintf(f_output, "%.6e\t%.6e\t", ptr->eig(mode, 0), ptr->eig(mode, 1));
			fprintf(f_output, "\n");
		}
		else
			db.myprintf("\nConcomitant solution did not present results!\n");
	}
}

void ConcomitantSolution::EndConcomitantSolution()
{
	fclose(f_output);
}
