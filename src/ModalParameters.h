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

	int sample;													//Amostragem para realização do problema de autovalor e salvar matrizes
	double tolerance;											//Tolerância para extração de autovalores - Arpack
	int number_frames;											//Número de frames para pós-processar os modos (vtk)
	int number_modes;											//Número de modos a serem calculados
	bool export_matrices;										//Exportar ou não as matrizes para cálculo externo ao Giraffe
	bool compute_eigenvectors;									//Salva ou não os autovetores

	bool concomitant_modal;										//Flag que indica se deve ou não ser realizada análise modal concomitante à estática/dinâmica
};

