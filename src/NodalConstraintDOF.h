#pragma once
#include "SpecialConstraint.h"
#include "BoolTable.h"

class NodalConstraintDOF:
	public SpecialConstraint
{
public:
	NodalConstraintDOF();
	~NodalConstraintDOF();

	bool Read(FILE* f);		//Leitura
	void Write(FILE* f);	//Gravação
	void Mount();			//Montagem dos resíduos e rigidez tangente
	void MountGlobal();		//Preenche a contribuição do elemento nas matrizes globais
	void PreCalc();			//Pré-cálculo de variáveis que é feito uma única vez no início
	void SaveLagrange();	//Salvando variáveis da configuração convergida
	void ActivateDOFs();	//Checa quais multiplicadores de lagrange serão ativados,de acordo com a ativação dos GLs dos nós dos quais a special constraint participa
	bool Check();			//Checa inconsistências no SC para evitar erros de execução
	void ComputeVelAccel();		//Computa efeito das condições iniciais nos nós da restrição
	void ComputeInitialGuessDisplacements();
	void WriteVTK_XMLRender(FILE* f);
	void WriteMonitor(FILE *f, bool first_record, double time);	//Escreve no monitor do SpecialConstraint
	//Variáveis
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
	//Matriz de rigidez tangente e vetor resíduo
	Matrix* stiffness;
	Matrix* residual;
	Matrix r;

};

