#pragma once
#include <stdio.h>

class BEM
{
public:
	BEM();
	~BEM();

	bool Read(FILE *f);					//Leitura
	void Write(FILE *f);				//Gravação
	bool Check();						//Checking inconsistencies
	//variaveis internas
	int B;		//Numero de pas
	double R;		//Raio do rotor
	double Rhub;	//Raio do hub
	int CS_rotor;	//Sistema de coordenadas para descrever o plano do rotor (plano xz)
	int node_rotor;	//Nó associado a posição do rotor

	double tol_bem;	//tolerancia do BEM
};

