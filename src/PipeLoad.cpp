#include "PipeLoad.h"

#include "IO.h"
#include "Table.h"
#include "MathCode.h"
#include "ElementSet.h"
#include "Element.h"
#include "Pipe_1.h"


PipeLoad::PipeLoad()
{
	number = 0;
	table = NULL;
	mcode = NULL;
	n_times = 0;
	n_values = 4;
	//(0,0) P0I
	//(1,0) P0E
	//(2,0) RhoI
	//(3,0) RhoE
}


PipeLoad::~PipeLoad()
{
	if (table != NULL)
		delete table;
	if (mcode != NULL)
		delete mcode;
}

//Pre-calculus
void PipeLoad::PreCalc()
{

}

//Atualiza dados necessários e que sejam dependentes de DOFs ativos/inativos - chamado no início de cada solution step
void PipeLoad::UpdateforSolutionStep()
{

}

bool PipeLoad::Read(FILE *f)
{
	char s[1000];
	fscanf(f, "%s", s);
	number = atoi(s);

	//Leitura do ElementSet
	fscanf(f, "%s", s);
	if (!strcmp(s, "ElementSet"))
	{
		fscanf(f, "%s", s);
		element_set = atoi(s);
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "NTimes"))
	{
		fscanf(f, "%s", s);
		n_times = atoi(s);
		//Allocating table
		table = new Table(n_times, n_values);
		TryComment(f);
		//Reads table data
		if (!table->Read(f))
			return false;

	}
	else
	{
		if (!strcmp(s, "MathCode"))
		{
			//Allocating MathCode
			mcode = new MathCode(n_values);
			//Reads MathCode
			if (!mcode->Read(f))
				return false;
		}
		else
			return false;
	}
	return true;
}

//Checking inconsistencies
bool PipeLoad::Check()
{
	if (element_set > db.number_element_sets)
		return false;
	int element = 0;
	for (int index = 0; index < db.element_sets[element_set - 1]->n_el; index++)
	{
		element = db.element_sets[element_set - 1]->el_list[index];
		if (typeid(*db.elements[element - 1]) != typeid(Pipe_1))
		{
			db.myprintf("Element number %d assigned to Load number %d.\nIt is not of the kind 'Pipe_1'.\n", element, number);
			return false;
		}	
	}
	return true;
}

void PipeLoad::Write(FILE *f)
{
	fprintf(f, "PipeLoad\t%d\tElementSet\t%d\t", number, element_set);
	if (table != NULL)
	{
		fprintf(f, "NTimes\t%d\n", n_times);
		table->Write(f);
	}
	if (mcode != NULL)
	{
		fprintf(f, "MathCode\n");
		mcode->Write(f);
	}
}

void PipeLoad::Mount()
{
	//Percorre elementos
	int element;
	for (int index = 0; index < db.element_sets[element_set - 1]->n_el; index++)
	{
		element = db.element_sets[element_set - 1]->el_list[index];
		if (typeid(*db.elements[element - 1]) == typeid(Pipe_1))
		{
			Pipe_1* ptr = static_cast<Pipe_1*>(db.elements[element - 1]);
			ptr->MountPipeSpecialLoads(number);
		}
		else
			db.myprintf("Element number %d assigned to Load number %d was ignored!\n", element, number);
	}
}


void PipeLoad::EvaluateExplicit(double t)
{

}

//Writes VTK XML data for post-processing
void PipeLoad::WriteVTK_XML(FILE *f)
{

}
