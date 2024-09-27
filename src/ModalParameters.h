#pragma once
#include <stdio.h>

class ModalParameters
{
public:
	ModalParameters();
	~ModalParameters();
	bool Read(FILE *f);
	void Write(FILE *f);

	int sample;													//Amostragem para realiza��o do problema de autovalor e salvar matrizes
	double tolerance;											//Tolerancia para extra��o de autovalores - Arpack
	int number_frames;											//Numero de frames para p�s-processar os modos (vtk)
	int number_modes;											//Numero de modos a serem calculados
	bool export_matrices;										//Exportar ou n�o as matrizes para calculo externo ao Giraffe
	bool compute_eigenvectors;									//Salva ou n�o os autovetores

	bool concomitant_modal;										//Flag que indica se deve ou n�o ser realizada analise modal concomitante a estatica/dinamica
};

