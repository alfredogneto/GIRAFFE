#pragma once
#include "SpecialConstraint.h"

#include "Matrix.h"
class RigidNodeSet :
	public SpecialConstraint
{
public:
	RigidNodeSet();
	~RigidNodeSet();

	bool Read(FILE *f);		//Leitura
	void Write(FILE *f);	//Gravação
	void Mount();			//Montagem dos residuos e rigidez tangente
	void MountGlobal();		//Preenche a contribuição do elemento nas matrizes globais
	void PreCalc();			//Pre-calculo de variaveis que e feito uma unica vez no inicio
	void SaveLagrange();	//Salvando variaveis da configuração convergida
	void ActivateDOFs();	//Checa quais multiplicadores de lagrange serão ativados,de acordo com a ativação dos GLs dos nós dos quais a special constraint participa
	void Alloc();			//Aloca matrizes
	void EvaluateM2M3(double* v, double** M2, double** M3, double** Qp, double* alphap, double* di, double* lamb);
	bool Check();			//Checa inconsistências no SC para evitar erros de execução
	void ComputeVelAccel();		//Computa efeito das condições iniciais nos nós da restrição
	void ComputeInitialGuessDisplacements();
	void WriteVTK_XMLRender(FILE *f);
	void WriteMonitor(FILE *f, bool first_record, double time);	//Escreve no monitor do SpecialConstraint
	//Variaveis
	int pilot_node;			//Nó piloto do corpo rigido
	int node_set_ID;		//Numero do node set vinculado ao nó piloto
	int n_nodes_set;		//Numero de nós do set vinculado ao nó piloto
	//Flag - ativação de GL de rotação
	bool flag_rotation;
	
	Matrix I3;
	//Matriz de rigidez tangente e vetor residuo
	Matrix** stiffness1;
	Matrix** residual1;
	Matrix** stiffness2;
	Matrix** residual2;
	Matrix** distancei;

	double v[600];
};

