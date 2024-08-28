#pragma once
#include <stdio.h>
#include "BoolTable.h"

class Contact
{
public:
	Contact() {}
	virtual ~Contact() {}

	char* type_name;											//Nome do tipo do elemento
	int number;													//Numero 
	bool** activate;											//Based on pinball check, contains 1 or 0 info to activate/not the contact
	bool* typeOK1;												//Contains 1 or 0, based on verification if the type of the underlying element is compatible with the kind of contact
	bool* typeOK2;												//Contains 1 or 0, based on verification if the type of the underlying element is compatible with the kind of contact
	BoolTable bool_table;										//Bool table que ativa ou desativa na seq. de solutions
	virtual void Mount() = 0;									//Monta contatos
	virtual void MountGlobal() = 0;								//Preenche a contribuição do contato nas matrizes globais
	virtual void Band(int* band_fixed, int* band_free) = 0;		//Calcula a banda gerada na matriz global pelo contato

	virtual bool Read(FILE *f) = 0;
	virtual void Write(FILE *f) = 0;
	virtual void WriteResults(FILE *f) = 0;										//Escreve arquivo de resultados
	virtual void WriteMonitor(FILE *f, bool first_record, double time) = 0;		//Escreve no monitor do contato
	virtual void WriteVTK_XMLRender(FILE *f) = 0;
	virtual void WriteVTK_XMLForces(FILE *f) = 0;
	virtual void PreCalc() = 0;													//Pre-calculo de variaveis que e feito uma unica vez no inicio
	virtual void PinballCheck() = 0;											//Checks proximity between each beam from LR to AS
	virtual void BeginStepCheck() = 0;											//checagem inicial do contato  - inicio de cada incremento
	virtual void SaveLagrange() = 0;											//Salva variaveis para descrição lagrangiana atualizada
	virtual bool HaveErrors() = 0;												//Retorna 1 - ha algum erro, mesmo que tenha convergido. Ex: penetração excessiva
	virtual void MountDyn() = 0;												//Montagens - Newmark
	virtual bool Check() = 0;													//Checa inconsistências no elemento para evitar erros de execução
};

