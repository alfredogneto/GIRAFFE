#pragma once
#include "Matrix.h"
#include <stdio.h>
#include <math.h>

class LagrangeSaveLRAS
{
public:
	LagrangeSaveLRAS(int e_elements);
	~LagrangeSaveLRAS();

	int **contact_status;								//Indica true se houver contato, e false se n�o houver. Para controle do ponto de in�cio de contato
	Matrix*** x0prev;									//Ponto anterior de contato
	double **gt1s;										//Deslizamentos acumulados nas dire��es 1 e 2
	double **gt2s;										//Deslizamentos acumulados nas dire��es 1 e 2
	double **gn;										//gn
	double **Fx;										//For�a total de contato em x
	double **Fy;										//For�a total de contato em y
	double **Fz;										//For�a total de contato em z

	int n_elements;
};

