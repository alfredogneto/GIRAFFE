#pragma once
#pragma once
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "Matrix.h"
#include "BoolTable.h"
#define PI 3.1415926535897932384626433832795

class SpecialConstraint
{
public:
	SpecialConstraint();
	virtual ~SpecialConstraint();

	virtual bool Read(FILE *f) = 0;		//Leitura
	virtual void Write(FILE *f) = 0;	//Grava��o
	virtual void Mount() = 0;			//Montagem dos res�duos e rigidez tangente
	virtual void MountGlobal() = 0;		//Preenche a contribui��o do elemento nas matrizes globais
	virtual void PreCalc() = 0;			//Pr�-c�lculo de vari�veis que � feito uma �nica vez no in�cio
	virtual void SaveLagrange() = 0;	//Salvando vari�veis da configura��o convergida
	virtual void ActivateDOFs() = 0;	//Realiza ativa��o dos GLs dos n�s dos quais a special constraint participa
	virtual bool Check() = 0;			//Checa inconsist�ncias no SC para evitar erros de execu��o
	virtual void ComputeVelAccel() = 0;	//Computa efeito das condi��es iniciais nos n�s da restri��o
	virtual void WriteVTK_XMLRender(FILE *f) = 0;
	virtual void ComputeInitialGuessDisplacements() = 0;
	virtual void WriteMonitor(FILE *f, bool first_record, double time) = 0;	//Escreve no monitor do SpecialConstraint
	int number;							//ID da Special Constraint

	int* active_lambda;					//Controle da ativa��o ou n�o dos multiplicadores de lagrange, de acordo com a ativa��o ou n�o dos GLs dos n�s associados ao v�nculo					
	//1 - GL ativo
	//0 - GL inativo
	int* GLs;							//N�meros dos GLs dos multiplicadores de lagrange
	int n_GL;							//N�mero de multiplicadores de lagrange
	double* lambda;						//Multiplicadores de lagrange do v�nculo
	double* copy_lambda;				//Multiplicadores de lagrange do v�nculo (�ltima configura��o convergida)
	BoolTable bool_table;				//Bool table que ativa ou desativa na seq. de solutions
};

