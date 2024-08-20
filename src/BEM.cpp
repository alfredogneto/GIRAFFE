#include "BEM.h"
#include"Database.h"

//Variáveis globais
extern
Database db;

BEM::BEM()
{
	B = 0;
	R = 0;
	Rhub = 0;
	tol_bem = 1e-6;
	CS_rotor = 0;
	node_rotor = 0;
}


BEM::~BEM()
{
}

//Leitura
bool BEM::Read(FILE *f)
{
	char s[1000];
	fscanf(f, "%s", s);
	//Verifica a palavra chave "NumberBlades"
	if (!strcmp(s, "NumberBlades"))
	{
		fscanf(f, "%s", s);
		B = atoi(s);
	}
	else
		return false;

	fscanf(f, "%s", s);
	//Verifica a palavra chave "RotorRadius"
	if (!strcmp(s, "RotorRadius"))
	{
		fscanf(f, "%s", s);
		R = atof(s);
	}
	else
		return false;

	fscanf(f, "%s", s);
	//Verifica a palavra chave "HubRadius"
	if (!strcmp(s, "HubRadius"))
	{
		fscanf(f, "%s", s);
		Rhub = atof(s);
	}
	else
		return false;

	fscanf(f, "%s", s);
	//Verifica a palavra chave "RotorCS"
	if (!strcmp(s, "RotorCS"))
	{
		fscanf(f, "%s", s);
		CS_rotor = atoi(s);
	}
	else
		return false;

	fscanf(f, "%s", s);
	//Verifica a palavra chave "RotorNode"
	if (!strcmp(s, "RotorNode"))
	{
		fscanf(f, "%s", s);
		node_rotor = atoi(s);
	}
	else
		return false;

	return true;
}

//Checking inconsistencies
bool BEM::Check()
{
	if (CS_rotor > db.number_CS)
		return false;
	if (node_rotor > db.number_nodes)
		return false;
	return true;
}

//Gravação
void BEM::Write(FILE *f)
{
	fprintf(f, "BEM\nNumberBlades\t%d\nRotorRadius\t%.6e\nHubRadius\t%.6e\tRotorCS\t%d\tRotorNode\t%d\n", B,R,Rhub, CS_rotor, node_rotor);
}
