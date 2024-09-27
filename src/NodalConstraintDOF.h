#pragma once
#include "SpecialConstraint.h"
#include "BoolTable.h"
#include "Matrix.h"

class NodalConstraintDOF:
	public SpecialConstraint
{
public:
	NodalConstraintDOF();
	~NodalConstraintDOF();

	bool Read(FILE* f);		//Leitura
	void Write(FILE* f);	//Grava��o
	void Mount();			//Montagem dos residuos e rigidez tangente
	void MountGlobal();		//Preenche a contribui��o do elemento nas matrizes globais
	void PreCalc();			//Pre-calculo de variaveis que e feito uma unica vez no inicio
	void SaveLagrange();	//Salvando variaveis da configura��o convergida
	void ActivateDOFs();	//Checa quais multiplicadores de lagrange ser�o ativados,de acordo com a ativa��o dos GLs dos n�s dos quais a special constraint participa
	bool Check();			//Checa inconsist�ncias no SC para evitar erros de execu��o
	void ComputeVelAccel();		//Computa efeito das condi��es iniciais nos n�s da restri��o
	void ComputeInitialGuessDisplacements();
	void WriteVTK_XMLRender(FILE* f);
	void WriteMonitor(FILE *f, bool first_record, double time);	//Escreve no monitor do SpecialConstraint
	//Variaveis
	int node_A;
	int node_B;

	//BoolTables for each kind of DOF
	BoolTable UX_table;
	BoolTable UY_table;
	BoolTable UZ_table;
	BoolTable ROTX_table;
	BoolTable ROTY_table;
	BoolTable ROTZ_table;

	Matrix I6;
	//Matriz de rigidez tangente e vetor residuo
	Matrix* stiffness;
	Matrix* residual;
	Matrix r;

};

