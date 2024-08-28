#pragma once

#include "BoolTable.h"

class SpecialConstraint
{
public:
	SpecialConstraint() {}
	virtual ~SpecialConstraint() {}

	virtual bool Read(FILE *f) = 0;		//Leitura
	virtual void Write(FILE *f) = 0;	//Gravação
	virtual void Mount() = 0;			//Montagem dos residuos e rigidez tangente
	virtual void MountGlobal() = 0;		//Preenche a contribuição do elemento nas matrizes globais
	virtual void PreCalc() = 0;			//Pre-calculo de variaveis que e feito uma unica vez no inicio
	virtual void SaveLagrange() = 0;	//Salvando variaveis da configuração convergida
	virtual void ActivateDOFs() = 0;	//Realiza ativação dos GLs dos nós dos quais a special constraint participa
	virtual bool Check() = 0;			//Checa inconsistências no SC para evitar erros de execução
	virtual void ComputeVelAccel() = 0;	//Computa efeito das condições iniciais nos nós da restrição
	virtual void WriteVTK_XMLRender(FILE *f) = 0;
	virtual void ComputeInitialGuessDisplacements() = 0;
	virtual void WriteMonitor(FILE *f, bool first_record, double time) = 0;	//Escreve no monitor do SpecialConstraint
	int number;							//ID da Special Constraint

	int* active_lambda;					//Controle da ativação ou não dos multiplicadores de lagrange, de acordo com a ativação ou não dos GLs dos nós associados ao vinculo					
	//1 - GL ativo
	//0 - GL inativo
	int* GLs;							//Numeros dos GLs dos multiplicadores de lagrange
	int n_GL;							//Numero de multiplicadores de lagrange
	double* lambda;						//Multiplicadores de lagrange do vinculo
	double* copy_lambda;				//Multiplicadores de lagrange do vinculo (ultima configuração convergida)
	BoolTable bool_table;				//Bool table que ativa ou desativa na seq. de solutions
};

