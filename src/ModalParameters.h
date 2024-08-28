#pragma once
#include <stdio.h>

class ModalParameters
{
public:
	ModalParameters();
	~ModalParameters();
	bool Read(FILE *f);
	void Write(FILE *f);

	int sample;													//Amostragem para realização do problema de autovalor e salvar matrizes
	double tolerance;											//Tolerancia para extração de autovalores - Arpack
	int number_frames;											//Numero de frames para pós-processar os modos (vtk)
	int number_modes;											//Numero de modos a serem calculados
	bool export_matrices;										//Exportar ou não as matrizes para calculo externo ao Giraffe
	bool compute_eigenvectors;									//Salva ou não os autovetores

	bool concomitant_modal;										//Flag que indica se deve ou não ser realizada analise modal concomitante a estatica/dinamica
};

