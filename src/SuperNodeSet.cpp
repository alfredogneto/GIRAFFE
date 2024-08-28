#include "SuperNodeSet.h"

#include"Database.h"
//Variaveis globais
extern
Database db;

SuperNodeSet::SuperNodeSet()
{
	n_super_nodes = 0;
	super_node_list = NULL;

	sequence = false;
	list = false;
	initial = 0;
	increment = 0;
}

SuperNodeSet::~SuperNodeSet()
{
	if (super_node_list != NULL)
	{
		delete[] super_node_list;
	}
}

bool SuperNodeSet::Read(FILE *f)
{
	char s[1000];
	fscanf(f, "%s", s);
	//Verifica a palavra chave "Set"
	if (!strcmp(s, "SuperNodeSet"))
	{
		fscanf(f, "%s", s);
		number = atoi(s);
	}
	else
		return false;
	fscanf(f, "%s", s);
	//Verifica a palavra chave "SuperNodes"
	if (!strcmp(s, "SuperNodes"))
	{
		fscanf(f, "%s", s);
		n_super_nodes = atoi(s);
		//Alocação do vetor de nós
		super_node_list = new int[n_super_nodes];
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
		for (int i = 0; i < n_super_nodes; i++)
		{
			fscanf(f, "%s", s);//Leitura do numero do nó
			super_node_list[i] = atoi(s);
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
			//Geração da lista de nós
			for (int i = 0; i < n_super_nodes; i++)
			{
				super_node_list[i] = initial + i * increment;
			}
		}
		else
			return false;
	}
	//Se atingiu esse ponto, sinal de leitura correta de tudo: retorna true
	return true;
}

void SuperNodeSet::Write(FILE *f)
{
	fprintf(f, "SuperNodeSet\t%d\tNodes\t%d\tList\t", number, n_super_nodes);
	for (int i = 0; i < n_super_nodes; i++)
		fprintf(f, "%d\t", super_node_list[i]);
	fprintf(f, "\n");
}

void SuperNodeSet::WriteMonitor(FILE *f, bool first_record, double time)
{
	//TO DO
}

