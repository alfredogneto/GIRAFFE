#pragma once
#include <stdio.h>
#include "Pipe_1.h"
#include "Matrix.h"
#include "MonitorNodesUserDefined.h"

#ifdef I
#undef I
#endif

class Node
{
public:
	Node(int e_nGL);
	~Node();
	//[0] X - coordenada do n�
	//[1] Y - coordenada do n�
	//[2] Z - coordenada do n�
	//[3] RX - coordenada do n�
	//[4] RY - coordenada do n�
	//[5] RZ - coordenada do n�

	double* ref_coordinates;	//Coordenadas do n� na configura��o de refer�ncia
	double* copy_coordinates;	//Coordenadas do n� na configura��o de c�pia (ultima convergida)
	double* copy_rot_euler;		//Coordenadas de rota��o (Euler) do n� na configura��o de c�pia (ultima convergida)
	//double* rot_axes;
	double* displacements;		//Deslocamentos do n� em rela��o a ultima c�pia de coordenadas
	double* vel;				//Velocidades
	double* copy_vel;			//Velocidades (c�pia da ultima convergida)
	double* accel;				//Acelera��es
	double* copy_accel;			//Acelera��es (c�pia da ultima convergida)

	int* constraints;
	//[0] Constraint em X - restri��o do n�
	//[1] Constraint em Y - restri��o do n�
	//[2] Constraint em Z - restri��o do n�
	//[3] Constraint em RX - restri��o do n�
	//[4] Constraint em RY - restri��o do n�
	//[5] Constraint em RZ - restri��o do n�

	int* GLs;
	//Possui numera��es dos graus de liberdade do n�:
	//Positivo - GL livre
	//Negativo - GL prescrito
	//N�o ha valor ZERO, tanto os GL livres, como prescritos iniciam sua numera��o do numero UM
	/*
	[0] X  -  GL do n�
	[1] Y  -  GL do n�
	[2] Z  -  GL do n�
	[3] RX  -  GL do n�
	[4] RY  -  GL do n�
	[5] RZ  -  GL do n�
	[6]... - outros graus de liberdade possiveis
	*/
	int* active_GL;
	//Possui informa��es da utiliza��o do GL
	/*
	1 - ativo
	0 - inativo
	*/
	int n_GL_free;
	int n_GL_fixed;

	int number;	//Numero do n� - diretamente relacionado com a indexa��o do vetor de n�s do database - informa��o redundante para facilitar possiveis verifica��es

	bool Read(FILE *f);
	void PreCalc();
	void Write(FILE *f);
	void WriteUpdated(FILE *f);
	void WriteInitialConditions(FILE *f);
	void WriteResults(FILE *f);
	void WriteVTK(FILE *f);
	void WriteMonitor(FILE *f, bool first_record, double time);
	void WriteMonitorUserDef(FILE *f, bool first_record, double time, const UserDefMonitorParams& params);
	void SaveConfiguration();
	void ZeroIC();
	Matrix InvXi(Matrix& alpha);


	int nGL;	//Numero de graus de liberdade presentes no n�

	//Variaveis internas para convers�o de rota��o e calculos do monitor
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
	Matrix* Q0; //Descri��o material (transforma��o de coordenadas inicial)


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

	std::vector<Pipe_1*> attached_pipes;
};

