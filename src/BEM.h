#pragma once
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "Matrix.h"
#include <math.h>
class BEM
{
public:
	BEM();
	~BEM();

	bool Read(FILE *f);					//Leitura
	void Write(FILE *f);				//Grava��o
	bool Check();						//Checking inconsistencies
	//vari�veis internas
	int B;		//N�mero de p�s
	double R;		//Raio do rotor
	double Rhub;	//Raio do hub
	int CS_rotor;	//Sistema de coordenadas para descrever o plano do rotor (plano xz)
	int node_rotor;	//N� associado � posi��o do rotor

	double tol_bem;	//toler�ncia do BEM
};

