#pragma once
#include <stdio.h>

class ModalParameters
{
public:
	ModalParameters();
	~ModalParameters();
	bool Read(FILE *f);
	void Write(FILE *f);

	int sample;													//Amostragem for performing of the eigenvalue problem and save matrices
	double tolerance;											//Tolerance for extração of eigenvalues - Arpack
	int number_frames;											//Number of frames for pós-processar the modos (vtk)
	int number_modes;											//Number of modos the to be calculated
	bool export_matrices;										//Exportar or not the matrices for calculation external to the Giraffe
	bool compute_eigenvectors;									//Saves or not the eigenvectors

	bool concomitant_modal;										//Flag that indicates if must or not be realizada analysis modal concomitante the estatica/dynamics
};

