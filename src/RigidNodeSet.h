#pragma once
#include "SpecialConstraint.h"
#include "NodeSet.h"
#include "SameRotation.h"

class RigidNodeSet :
	public SpecialConstraint
{
public:
	RigidNodeSet();
	~RigidNodeSet();

	bool Read(FILE *f);		//Leitura
	void Write(FILE *f);	//Grava��o
	void Mount();			//Montagem dos res�duos e rigidez tangente
	void MountGlobal();		//Preenche a contribui��o do elemento nas matrizes globais
	void PreCalc();			//Pr�-c�lculo de vari�veis que � feito uma �nica vez no in�cio
	void SaveLagrange();	//Salvando vari�veis da configura��o convergida
	void ActivateDOFs();	//Checa quais multiplicadores de lagrange ser�o ativados,de acordo com a ativa��o dos GLs dos n�s dos quais a special constraint participa
	void Alloc();			//Aloca matrizes
	void EvaluateM2M3(double* v, double** M2, double** M3, double** Qp, double* alphap, double* di, double* lamb);
	bool Check();			//Checa inconsist�ncias no SC para evitar erros de execu��o
	void ComputeVelAccel();		//Computa efeito das condi��es iniciais nos n�s da restri��o
	void ComputeInitialGuessDisplacements();
	void WriteVTK_XMLRender(FILE *f);
	void WriteMonitor(FILE *f, bool first_record, double time);	//Escreve no monitor do SpecialConstraint
	//Vari�veis
	int pilot_node;			//N� piloto do corpo r�gido
	int node_set_ID;		//N�mero do node set vinculado ao n� piloto
	int n_nodes_set;		//N�mero de n�s do set vinculado ao n� piloto
	//Flag - ativa��o de GL de rota��o
	bool flag_rotation;
	
	Matrix I3;
	//Matriz de rigidez tangente e vetor res�duo
	Matrix** stiffness1;
	Matrix** residual1;
	Matrix** stiffness2;
	Matrix** residual2;
	Matrix** distancei;

	double v[600];
};

