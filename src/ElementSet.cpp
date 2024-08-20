#include "ElementSet.h"

#include"Database.h"
//Variáveis globais
extern
Database db;

ElementSet::ElementSet()
{
	n_el = 0;
	el_list = NULL;

	sequence = false;
	list = false;
	initial = 0;
	increment = 0;
}

ElementSet::~ElementSet()
{
	if (el_list != NULL)
	{
		delete[] el_list;
	}
}

bool ElementSet::Read(FILE *f)
{
	char s[1000];
	fscanf(f, "%s", s);
	//Verifica a palavra chave "ElementSet"
	if (!strcmp(s, "ElementSet"))
	{
		fscanf(f, "%s", s);
		number = atoi(s);
	}
	else
		return false;
	fscanf(f, "%s", s);
	//Verifica a palavra chave "Elements"
	if (!strcmp(s, "Elements"))
	{
		fscanf(f, "%s", s);
		n_el = atoi(s);
		//Alocação do vetor de nós
		el_list = new int[n_el];
	}
	else
		return false;
	//Duas possibilidades de leitura:
	//1 - List
	//2 - Sequence
	fscanf(f, "%s", s);
	if (!strcmp(s, "List"))
	{
		list = true;
		for (int i = 0; i < n_el; i++)
		{
			fscanf(f, "%s", s);//Leitura do número do elemento
			el_list[i] = atoi(s);
		}
	}
	else
	{
		if (!strcmp(s, "Sequence"))
		{
			sequence = true;
			fscanf(f, "%s", s);
			if (!strcmp(s, "Initial"))
			{
				fscanf(f, "%s", s);
				initial = atoi(s);
			}
			else
				return false;
			fscanf(f, "%s", s);
			if (!strcmp(s, "Increment"))
			{
				fscanf(f, "%s", s);
				increment = atoi(s);
			}
			else
				return false;
			//Geração da lista de elementos
			for (int i = 0; i < n_el; i++)
			{
				el_list[i] = initial + i*increment;
			}
		}
		else
			return false;
	}
	//Se atingiu esse ponto, sinal de leitura correta de tudo: retorna true
	return true;
}

void ElementSet::Write(FILE *f)
{
	fprintf(f, "ElementSet\t%d\tSurfaces\t%d\tList\t", number, n_el);
	for (int i = 0; i < n_el; i++)
		fprintf(f, "%d\t", el_list[i]);
	fprintf(f, "\n");
}
