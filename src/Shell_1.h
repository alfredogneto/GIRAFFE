#pragma once
#include "Element.h"
#include "Matrix.h"

class Shell_1 :
	public Element
{
public:
	Shell_1();
	~Shell_1();
	bool Check();				//Checks inconsistencies in the element for avoid errors of execution
	bool Read(FILE *f);
	void Write(FILE *f);
	void WriteResults(FILE *f);									//Writes file of results
	void WriteMonitor(FILE *f, bool first_record, double time);	//Writes in the monitor of the element
	void WriteVTK_XMLBase(std::vector<float> *float_vector);
	void WriteVTK_XMLRender(FILE *f);
	void Mount();												//Assembles elements
	void MountElementLoads();									//Assembles loadings associated to the element
	void MountMass();											//Assembles the mass matrix
	void MountMassModal();										//Assembles the mass matrix for performing of the analysis modal
	void MountDampingModal();									//Assembles the damping matrix for performing of the analysis modal
	void MountDamping(bool update_rayleigh);					//Assembles the damping matrix
	void MountDyn();											//Assemblies - Newmark
	void MountDynModal();										//Assemblies for analysis modal - insertion of the mass matrix and damping in the stiffness matrix for subsequent assembly global
	void TransformMatrix();										//Assembles transformation matrix of coordinates
	void MountGlobal();											//Fills the contribution of the element in the matrices global
	void MountFieldLoads();										//Assembles loadings of field (ex: weight próprio)
	void MountShellSpecialLoads(int l_number);					//Assembles loadings of pressão in the shell

	void SaveLagrange();										//Saves variables in the gauss points useful for description lagrangiana updated
	void PreCalc();												//Pre-calculation of variables that and done a single time in the start
	void Zeros();												//Zeroes matrices

	//Functions of way and its derivatives
	double* N1u;												
	double* N2u;
	double* N3u;
	double* N4u;
	double* N5u;
	double* N6u;
	double* N4a;
	double* N5a;
	double* N6a;

	double* N1u_x1;
	double* N2u_x1;
	double* N3u_x1;
	double* N4u_x1;
	double* N5u_x1;
	double* N6u_x1;
	double* N4a_x1;
	double* N5a_x1;
	double* N6a_x1;

	double* N1u_x2;
	double* N2u_x2;
	double* N3u_x2;
	double* N4u_x2;
	double* N5u_x2;
	double* N6u_x2;
	double* N4a_x2;
	double* N5a_x2;
	double* N6a_x2;

	Matrix** N;							//6x27//Matrix of the functions of way
	Matrix** deltaN;					//15x27//Matrix deltaN (with functions of way and its derivatives) 

	//////////////////////////////// Functions of way for integration of degree 4 (6 gauss points)
	//Functions of way
	double* N1u4;
	double* N2u4;
	double* N3u4;
	double* N4u4;
	double* N5u4;
	double* N6u4;
	double* N4a4;
	double* N5a4;
	double* N6a4;

	double* N1u4_x1;
	double* N2u4_x1;
	double* N3u4_x1;
	double* N4u4_x1;
	double* N5u4_x1;
	double* N6u4_x1;

	double* N1u4_x2;
	double* N2u4_x2;
	double* N3u4_x2;
	double* N4u4_x2;
	double* N5u4_x2;
	double* N6u4_x2;

	Matrix** N4;							//6x27//Matrix of the functions of way
	Matrix** N4_u_x1;						//3x18//Matrix of the functions of way
	Matrix** N4_u_x2;						//3x18//Matrix of the functions of way
	double* w4;								//pesos
	Matrix** alpha_i4;
	////////////////////////////////////
	
	//Variables for descrever the cinematica
	Matrix** alpha_delta;
	Matrix** alpha_delta_x1;
	Matrix** alpha_delta_x2;
	Matrix** u_delta;
	Matrix** u_delta_x1;
	Matrix** u_delta_x2;
	Matrix* e1r;	//Direction 1 of the shell (global)
	Matrix* e2r;	//Direction 2	of the shell (global)
	Matrix* e3r;	//Direction 2	of the shell (global)
	Matrix* e1rlocal;	//Direction 1 of the shell (local)
	Matrix* e2rlocal;	//Direction 2	of the shell (local)
	Matrix* e3rlocal;	//Direction 3	of the shell (local)

	//Variables the to be saved - Lag. Updated
	Matrix** Q_i;
	Matrix** Q_delta;
	Matrix** Xi_delta;
	Matrix** u_i;
	Matrix** alpha_i;
	Matrix** kappa_r1_i;
	Matrix** kappa_r2_i;
	Matrix** kappa_r1_delta;
	Matrix** kappa_r2_delta;
	Matrix** z_x1_i;
	Matrix** z_x2_i;

	Matrix* stiffness;								//Stiffness matrix
	Matrix* mass;									//Mass matrix
	Matrix* mass_modal;								//Mass matrix for analysis modal
	Matrix* damping;								//Damping matrix
	Matrix* damping_modal;							//Damping matrix for analysis modal
	Matrix* damping_loading;						//Internal forces of damping
	Matrix* rayleigh_damping;						//Damping matrix initial of the problem
	Matrix* i_loading;								//Vector of internal forces internal
	Matrix* inertial_loading;						//Vector of internal forces of inertia
	Matrix* P_loading;								//Vector of internal forces unbalanced
	Matrix* e_loading;								//Vector of internal forces external

	Matrix* transform;								//Transformation matrix of coordinates
	Matrix* transform3;								//Transformation matrix of coordinates 3x3

	double lambda, mu;								//Lame
	double stiff_drill;								//Stiffness to the drilling
	double area;									//area of the element
	double alpha1;									//Pesos of integration

	double coef1, coef2, coef3;						//Coefficients for use in the mass matrix

	Matrix I3;

	Matrix** eta_r1;			//strain in er1
	Matrix** eta_r2;			//strain in er2
	Matrix** kappa_r1;			//rot/comp in er1
	Matrix** kappa_r2;			//rot/comp in er2
	Matrix** n_r1;				//force in er1
	Matrix** n_r2;				//force in er2
	Matrix** m_r1;				//moment in er1
	Matrix** m_r2;				//moment in er2

	Matrix** eta_r1_global;			//strain in er1
	Matrix** eta_r2_global;			//strain in er2
	Matrix** kappa_r1_global;		//rot/comp in er1
	Matrix** kappa_r2_global;		//rot/comp in er2
	Matrix** n_r1_global;			//force in er1
	Matrix** n_r2_global;			//force in er2
	Matrix** m_r1_global;			//moment in er1
	Matrix** m_r2_global;			//moment in er2

	//Variables internal for use in the dynamics
	Matrix** alpha_dot;
	Matrix** Xi_dot;
	Matrix** Mip;
	Matrix** Jip;
	Matrix** M;
	Matrix** Md1;
	Matrix v_ipp;//Estimativa of the velocity in the instant subsequent


	//Matrices for transformation of the results
	Matrix n_1;
	Matrix n_2;
	Matrix m_1;
	Matrix m_2;

	Matrix e_1;
	Matrix e_2;
	Matrix k_1;
	Matrix k_2;

	//Calculates the contributions inerciais for the analysis dynamics - inclui all the contributions for the way fraca and for the operator tangent
	void EvaluateInertialContributions(double* v, double(*a1)
		, double(*a2), double(*a3), double(*a4), double(*a5), double(*a6)
		, double* alphai, double* alphad, double* ui, double* ud, double* omegai
		, double* domegai, double* dui, double* ddui, double* e3r, double(*coef1)
		, double(*coef2), double* dT, double** DdT);

	void EvaluateMassModal(double* v, double* alphai, double* coef1, double* coef2, double* coef3, double* e3r, double** matrixm);

	//Variables for function generated in the AceGen
	double temp_v[2000];				//temporary variable for calculations internal
	Matrix *DdT;
	Matrix *dT;
	double** pDdT;
	double* tempkin;

	//Variables for calculation of steps
	int t1, t2;
	double load_multiplier, l_factor, mult;

	//COMPOSITE MATERIAL
	Matrix* D_comp;				//Matrix Constitutiva Compósito
	double thick_comp;			//Espessura of the Material Compósito
	double rho_comp;			//Rho of the material compósito
	double stiff_drill_comp1;	//Rigidiz to the drilling in the direction 1
	double stiff_drill_comp2;	//Stiffness to the drilling in the direction 2
	bool composite_shell;		//Checks if the shell and of material compósito


};