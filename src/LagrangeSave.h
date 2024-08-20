#pragma once
#include "Matrix.h"
#include <stdio.h>
#include <math.h>

class LagrangeSave
{
public:
	LagrangeSave(void);
	~LagrangeSave(void);

	//Variáveis salvas nos pontos de Gauss do elementos de viga (para serem utilizadas posteriormente sem necessidade de interpolação
	Matrix** Q_i;					//Tensor rotação no instante i
	Matrix** alpha_i;				//Vetor rotação no instante i 
	Matrix** d_alpha_i;				//Derivada temporal do vetor rotação no instante i
	Matrix** theta_i;				//Vetor rotação de Euler no instante i
	Matrix** kappa_i_ref;			//Vetor curvatura generalizada no instante i
	Matrix** u_i;					//Vetor deslocamento do eixo da barra no instante i
	Matrix** dz_i;					//Vetor derivada do deslocamento no instante i
	
	LagrangeSave &operator = (LagrangeSave const &save1);
};
