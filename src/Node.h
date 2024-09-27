#pragma once
#include <stdio.h>
#include "Matrix.h"

class Node
{
public:
	Node(int e_nGL);
	~Node();
	//[0] X - coordenada do nó
	//[1] Y - coordenada do nó
	//[2] Z - coordenada do nó
	//[3] RX - coordenada do nó
	//[4] RY - coordenada do nó
	//[5] RZ - coordenada do nó

	double* ref_coordinates;	//Coordenadas do nó na configuração de referência
	double* copy_coordinates;	//Coordenadas do nó na configuração de cópia (ultima convergida)
	double* copy_rot_euler;		//Coordenadas de rotação (Euler) do nó na configuração de cópia (ultima convergida)
	//double* rot_axes;
	double* displacements;		//Deslocamentos do nó em relação a ultima cópia de coordenadas
	double* vel;				//Velocidades
	double* copy_vel;			//Velocidades (cópia da ultima convergida)
	double* accel;				//Acelerações
	double* copy_accel;			//Acelerações (cópia da ultima convergida)

	int* constraints;
	//[0] Constraint em X - restrição do nó
	//[1] Constraint em Y - restrição do nó
	//[2] Constraint em Z - restrição do nó
	//[3] Constraint em RX - restrição do nó
	//[4] Constraint em RY - restrição do nó
	//[5] Constraint em RZ - restrição do nó

	int* GLs;
	//Possui numerações dos graus de liberdade do nó:
	//Positivo - GL livre
	//Negativo - GL prescrito
	//Não ha valor ZERO, tanto os GL livres, como prescritos iniciam sua numeração do numero UM
	/*
	[0] X  -  GL do nó
	[1] Y  -  GL do nó
	[2] Z  -  GL do nó
	[3] RX  -  GL do nó
	[4] RY  -  GL do nó
	[5] RZ  -  GL do nó
	[6]... - outros graus de liberdade possiveis
	*/
	int* active_GL;
	//Possui informações da utilização do GL
	/*
	1 - ativo
	0 - inativo
	*/
	int n_GL_free;
	int n_GL_fixed;

	int number;	//Numero do nó - diretamente relacionado com a indexação do vetor de nós do database - informação redundante para facilitar possiveis verificações

	bool Read(FILE *f);
	void Write(FILE *f);
	void WriteUpdated(FILE *f);
	void WriteInitialConditions(FILE *f);
	void WriteResults(FILE *f);
	void WriteVTK(FILE *f);
	void WriteMonitor(FILE *f, bool first_record, double time);
	void SaveConfiguration();
	void ZeroIC();
	Matrix InvXi(Matrix& alpha);

	int nGL;	//Numero de graus de liberdade presentes no nó

	//Variaveis internas para conversão de rotação e calculos do monitor
	double theta;
	Matrix* rot_euler;
	Matrix* rot_rodrigues;
	Matrix* load;
	Matrix* moment;
	Matrix* force;
	Matrix* A;
	double g;
	Matrix* I;
	Matrix* Xi;
	Matrix* Xi_T_inv;
	Matrix* Q;
	Matrix* I3;
	bool flag_material_description;
	bool flag_pseudo_moment;
	Matrix* Q0; //Descrição material (transformação de coordenadas inicial)


	double alpha_escalar;
	double theta_escalar;
	Matrix* alpha_v;
	Matrix* theta_v;

	//Variaveis locais
	Matrix alpha_1;
	Matrix alpha_2;
	Matrix alpha_3;

	Matrix theta3;
	double theta_escalar3;

	//Explicit
	Matrix* copy_du;
	Matrix* copy_omega;
	Matrix* u;
	Matrix* du;
	Matrix* ddu;
	Matrix* alpha;
	Matrix* omega;
	Matrix* domega;
};

