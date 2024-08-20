#pragma once
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "Matrix.h"

class ModalParameters
{
public:
	ModalParameters();
	~ModalParameters();
	bool Read(FILE *f);
	void Write(FILE *f);

	int sample;													//Amostragem para realiza��o do problema de autovalor e salvar matrizes
	double tolerance;											//Toler�ncia para extra��o de autovalores - Arpack
	int number_frames;											//N�mero de frames para p�s-processar os modos (vtk)
	int number_modes;											//N�mero de modos a serem calculados
	bool export_matrices;										//Exportar ou n�o as matrizes para c�lculo externo ao Giraffe
	bool compute_eigenvectors;									//Salva ou n�o os autovetores

	bool concomitant_modal;										//Flag que indica se deve ou n�o ser realizada an�lise modal concomitante � est�tica/din�mica
};

