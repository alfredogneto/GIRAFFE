#include "ShellLoad.h"
#include "IO.h"
#include "Table.h"
#include "MathCode.h"
#include "ElementSet.h"
#include "Element.h"
#include "Shell_1.h"

ShellLoad::ShellLoad()
{
	area_update = false;
	number = 0;
	table = NULL;
	mcode = NULL;
	n_times = 0;
	n_values = 1;
}


ShellLoad::~ShellLoad()
{
	if (table != NULL)
		delete table;
	if (mcode != NULL)
		delete mcode;
}

bool ShellLoad::Read(FILE *f)
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

	//Leitura do AreaUpdate
	fscanf(f, "%s", s);
	if (!strcmp(s, "AreaUpdate"))
	{
		fscanf(f, "%s", s);
		if (atoi(s) != 0 && atoi(s) != 1)//Se não for 0 ou 1
			return false;
		if (atoi(s) == 1)
			area_update = true;
		else
			area_update = false;
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
bool ShellLoad::Check()
{
	if (element_set > db.number_element_sets)
		return false;
	int element = 0;
	for (int index = 0; index < db.element_sets[element_set - 1]->n_el; index++)
	{
		element = db.element_sets[element_set - 1]->el_list[index];
		if (typeid(*db.elements[element - 1]) != typeid(Shell_1))
		{
			db.myprintf("Element number %d assigned to Load number %d.\nIt is not of the kind 'Shell_1'.\n", element, number);
			return false;
		}
	}
	return true;
}

void ShellLoad::Write(FILE *f)
{
	fprintf(f, "ShellLoad\t%d\tElementSet\t%d\tAreaUpdate\t%d\t", number, element_set, area_update);
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

//Pre-calculus
void ShellLoad::PreCalc()
{

}

//Atualiza dados necessários e que sejam dependentes de DOFs ativos/inativos - chamado no início de cada solution step
void ShellLoad::UpdateforSolutionStep()
{

}

void ShellLoad::Mount()
{
	//Percorre elementos
	int element;
	for (int index = 0; index < db.element_sets[element_set - 1]->n_el; index++)
	{
		element = db.element_sets[element_set - 1]->el_list[index];
		if (typeid(*db.elements[element - 1]) == typeid(Shell_1))
		{
			Shell_1* ptr = static_cast<Shell_1*>(db.elements[element - 1]);
			ptr->MountShellSpecialLoads(number);
		}
		else
			db.myprintf("Element number %d assigned to Load number %d was ignored!\n", element, number);
	}
}


void ShellLoad::EvaluateExplicit(double t)
{

}


//Writes VTK XML data for post-processing
void ShellLoad::WriteVTK_XML(FILE *f)
{

}
