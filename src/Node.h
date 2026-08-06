#pragma once
#include <stdio.h>
#include "Matrix.h"

class Node
{
public:
	Node(int e_nGL);
	~Node();
	//[0] X - coordinate of the node
	//[1] Y - coordinate of the node
	//[2] Z - coordinate of the node
	//[3] RX - coordinate of the node
	//[4] RY - coordinate of the node
	//[5] RZ - coordinate of the node

	double* ref_coordinates;	//Coordinates of the node in the configuration of reference
	double* copy_coordinates;	//Coordinates of the node in the configuration of copies (last converged)
	double* copy_rot_euler;		//Coordinates of rotation (Euler) of the node in the configuration of copies (last converged)
	//double* rot_axes;
	double* displacements;		//Displacements of the node in relação the last copies of coordinates
	double* vel;				//Velocities
	double* copy_vel;			//Velocities (copies of the last converged)
	double* accel;				//Acelerações
	double* copy_accel;			//Acelerações (copies of the last converged)

	int* constraints;
	//[0] Constraint in X - constraint of the node
	//[1] Constraint in Y - constraint of the node
	//[2] Constraint in Z - constraint of the node
	//[3] Constraint in RX - constraint of the node
	//[4] Constraint in RY - constraint of the node
	//[5] Constraint in RZ - constraint of the node

	int* GLs;
	//Has numerações of the degrees of freedom of the node:
	//Positive - GL free
	//Negative - GL prescrito
	//Not there is value ZERO, tanto the GL free, as prescritos iniciam its numbering of the number A
	/*
	[0] X  -  GL of the node
	[1] Y  -  GL of the node
	[2] Z  -  GL of the node
	[3] RX  -  GL of the node
	[4] RY  -  GL of the node
	[5] RZ  -  GL of the node
	[6]... - other degrees of freedom possiveis
	*/
	int* active_GL;
	//Has information of the utilização of the GL
	/*
	1 - active
	0 - inactive
	*/
	int n_GL_free;
	int n_GL_fixed;

	int number;	//Number of the node - diretamente relacionado with the indexing of the vector of nodes of the database - information redundante for facilitar possiveis verificações

	bool Read(FILE *f);
	void Write(FILE *f);
	void WriteUpdated(FILE *f);
	void WriteInitialConditions(FILE *f);
	void WriteResults(FILE *f);
	void WriteVTK(FILE *f);
	void WriteMonitor(FILE *f, bool first_record, double time);
	void SaveConfiguration();
	void ZeroIC();
	Matrix InvXi(const Matrix& alpha);

	int nGL;	//Number of degrees of freedom presentes in the node

	//Variables internal for conversion of rotation and calculations of the monitor
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
	Matrix* Q0; //Description material (transformation of coordinates initial)


	double alpha_escalar;
	double theta_escalar;
	Matrix* alpha_v;
	Matrix* theta_v;

	//Variables local
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

