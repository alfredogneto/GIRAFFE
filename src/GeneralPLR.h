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

	bool Check();				//Checa inconsistências no elemento para evitar erros de execução
	void Mount();													//Monta contatos
	void MountGlobal();												//Preenche a contribuição do contato nas matrizes globais
	bool Read(FILE *f);
	void Write(FILE *f);
	void WriteResults(FILE *f);										//Escreve arquivo de resultados
	void WriteMonitor(FILE *f, bool first_record, double time);		//Escreve no monitor do contato
	void PreCalc();													//Pré-cálculo de variáveis que é feito uma única vez no início
	void PinballCheck();											//Checks proximity between each beam from LR to AS
	void BeginStepCheck();											//checagem inicial do contato  - início de cada incremento
	void SaveLagrange();											//Salva variáveis para descrição lagrangiana atualizada
	bool HaveErrors();												//Retorna 1 - há algum erro, mesmo que tenha convergido. Ex: penetração excessiva
	void Band(int* band_fixed, int* band_free);						//Calcula a banda gerada na matriz global pelo contato
	void MountDyn();												//Montagens - Newmark
	void Alloc();													//Aloca estrutura de matrizes para endereçar vetores e matrizes devidos ao contato
	void AllocSpecific(int i, int j);								//Aloca matrizes para contato específico
	void FreeSpecific(int i, int j);								//Aloca matrizes para contato específico
	void WriteVTK_XMLRender(FILE *f);
	void WriteVTK_XMLForces(FILE *f);
	//Variáveis internas
	double mu, epn, ept, pinball;
	double c;														//Coeficiente de dissipação de energia
	Matrix z1, z2;													//Posições centrais das partículas no PinballSearch
	int n_particles;												//Número de partículas
	int n_elements;													//Numéro de elementos dentro do line region
	int n_LR;														//Número (índice) do line region associado ao contato
	Matrix ****c_stiffness;											//Matriz de rigidez
	Matrix ****c_damping;											//Matriz de amortecimento
	Matrix ****c_loading;											//Vetor de esforços externos
	bool** alloc_control;											//Controle de alocação dinâmica
	Matrix I3;
};

