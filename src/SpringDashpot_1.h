#pragma once
#include "Element.h"
#include "Matrix.h"

class SpringDashpot_1 :
	public Element
{
public:
	SpringDashpot_1();
	~SpringDashpot_1();
	bool Check();				//Checks inconsistencies in the element for avoid errors of execution
	bool Read(FILE *f);
	void Write(FILE *f);
	void WriteResults(FILE *f);//Writes file of results
	void WriteMonitor(FILE *f, bool first_record, double time);	//Writes in the monitor of the element
	void WriteVTK_XMLBase(std::vector<float> *float_vector);
	void WriteVTK_XMLRender(FILE *f);
	void Mount();			//Assembles elements
	void MountElementLoads();			//Assembles loadings associated to the element
	void MountMass();					//Assembles the mass matrix
	void MountMassModal();				//Assembles the mass matrix for performing of the analysis modal
	void MountDamping(bool update_rayleigh);				//Assembles the damping matrix
	void MountDampingModal();			//Assembles the damping matrix for performing of the analysis modal
	void MountDyn();					//Assemblies - Newmark
	void MountDynModal();				//Assemblies for analysis modal - insertion of the mass matrix and damping in the stiffness matrix for subsequent assembly global
	void TransformMatrix();	//Assembles transformation matrix of coordinates
	void MountGlobal();		//Fills the contribution of the element in the matrices global

	void SaveLagrange();	//Saves variables in the gauss points useful for description lagrangiana updated
	void PreCalc();		//Pre-calculation of variables that and done a single time in the start
	void Zeros();			//Zeroes matrices local of the element

	//Variables of the element
	double k;	//Stiffness
	double c;	//Damping
	double initial_distance;
	double gn;		//gap normal
	Matrix z1, z2;	//Positions of the nodes
	Matrix xd1, xd2;//Velocities of the nodes
	Matrix c_stiffness;									//Stiffness matrix
	Matrix c_damping;									//Damping matrix
	Matrix c_damping_modal;									//Damping matrix
	Matrix c_stiffness_force;							//Vector of internal forces elasticos
	Matrix c_damping_force;								//Vector of internal forces of damping
	Matrix I3;
	Matrix z1z2;	//distance current between nodes
	Matrix n;		//normal direction of the spring
	
	Matrix non;
	double f;
	Matrix C1;
	Matrix C2;
	Matrix C3;

	//Pós-processamento
	double elastic_force;
	double damping_force;

	Matrix last_n;	//last normal direction converged of the spring
	bool first_evaluation;
};

