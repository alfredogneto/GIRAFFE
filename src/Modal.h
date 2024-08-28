#pragma once
#include "Solution.h"
#include "Matrix.h"

class Modal:
	public Solution
{
public:
	Modal();
	~Modal();

	bool Read(FILE *f);											//Leitura de arquivo de entrada
	void Write(FILE *f);										//Grava��o de arquivo de sa�da
	bool Solve();												//Solves solution routine
	void CreateOutputFolder();									//Creates output folder
	void WriteResults(double time_value);						//Escreve resultados
	void WriteMatrices();										//Escreve as matrizes na pasta modal
	void WriteModes();											//Escreve arquivos com os modos de vibrar
	void WriteVTK_XML();										//Escreve modos de vibrar de todos os elementos

	//Functions to evaluate mode displacements
	void ComputeModalDisplacement(int node, int DOF, int mode, double* Re, double* Im);	//Calcula o valor do deslocamento modal no n� e grau de liberdade em quest�o
	double ComputeModeNorm(int mode);													//Computa a norma do autovetor
	
	//Input variables
	bool export_matrices;
	int number_modes;
	double tolerance;
	bool compute_eigenvectors;
	int number_frames;

	//Internal variables
	char copy_name[1000];										//Endereco do diretorio para salvar arquivos
	Matrix eig;													//Matrix com autovalores complexos
	Matrix eigenvectors;										//Matrix com autovetores complexos
	double mode_factor;											//Fator multiplicativo para o modo de vibrar
	bool concomitant_solution;									//Flag que indica se � an�lise concomitante
	
	//Special functions for matrices
	void ZerosLocalMatrices();				//Zera matrizes dos elementos (para uso na an�lise modal)
	void MountMassModal();					//Monta matriz de massa para an�lise modal
	void MountDampingModal();				//Monta matriz de amortecimento para an�lise modal
	void MountDynModal();					//Montagens para analise Modal

};

