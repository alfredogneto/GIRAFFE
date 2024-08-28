#pragma once
#include "Contact.h"
#include <stdio.h>

#include "Matrix.h"

class GeneralPLR :
	public Contact
{
public:
	GeneralPLR();
	~GeneralPLR();

	bool Check();				//Checa inconsist�ncias no elemento para evitar erros de execu��o
	void Mount();													//Monta contatos
	void MountGlobal();												//Preenche a contribui��o do contato nas matrizes globais
	bool Read(FILE *f);
	void Write(FILE *f);
	void WriteResults(FILE *f);										//Escreve arquivo de resultados
	void WriteMonitor(FILE *f, bool first_record, double time);		//Escreve no monitor do contato
	void PreCalc();													//Pr�-c�lculo de vari�veis que � feito uma �nica vez no in�cio
	void PinballCheck();											//Checks proximity between each beam from LR to AS
	void BeginStepCheck();											//checagem inicial do contato  - in�cio de cada incremento
	void SaveLagrange();											//Salva vari�veis para descri��o lagrangiana atualizada
	bool HaveErrors();												//Retorna 1 - h� algum erro, mesmo que tenha convergido. Ex: penetra��o excessiva
	void Band(int* band_fixed, int* band_free);						//Calcula a banda gerada na matriz global pelo contato
	void MountDyn();												//Montagens - Newmark
	void Alloc();													//Aloca estrutura de matrizes para endere�ar vetores e matrizes devidos ao contato
	void AllocSpecific(int i, int j);								//Aloca matrizes para contato espec�fico
	void FreeSpecific(int i, int j);								//Aloca matrizes para contato espec�fico
	void WriteVTK_XMLRender(FILE *f);
	void WriteVTK_XMLForces(FILE *f);
	//Vari�veis internas
	double mu, epn, ept, pinball;
	double c;														//Coeficiente de dissipa��o de energia
	Matrix z1, z2;													//Posi��es centrais das part�culas no PinballSearch
	int n_particles;												//N�mero de part�culas
	int n_elements;													//Num�ro de elementos dentro do line region
	int n_LR;														//N�mero (�ndice) do line region associado ao contato
	Matrix ****c_stiffness;											//Matriz de rigidez
	Matrix ****c_damping;											//Matriz de amortecimento
	Matrix ****c_loading;											//Vetor de esfor�os externos
	bool** alloc_control;											//Controle de aloca��o din�mica
	Matrix I3;
};

