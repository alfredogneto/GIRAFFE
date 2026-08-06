#pragma once
#include "Element.h"

class LagrangeSave;

class Beam_1 :
	public Element
{
public:
	Beam_1(void);
	~Beam_1(void);
	bool Check();						//Checks inconsistencies in the element for avoid errors of execution
	bool Read(FILE *f);
	void Write(FILE *f);
	void WriteResults(FILE *f);			//Writes file of results
	void WriteVTK_XMLBase(std::vector<float> *float_vector);
	void WriteVTK_XMLRender(FILE *f);
	void WriteMonitor(FILE *f, bool first_record, double time);	//Writes in the monitor of the element
	void Mount();						//Assembles elements
	void MountElementLoads();			//Assembles loadings associated to the element
	void MountMass();					//Assembles the mass matrix
	void MountMassModal();				//Assembles the mass matrix for performing of the analysis modal
	void MountDampingModal();			//Assembles the damping matrix for performing of the analysis modal
	void MountDamping(bool update_rayleigh);				//Assembles the damping matrix
	void MountDyn() ;					//Assemblies - Newmark
	void MountDynModal();				//Assemblies for analysis modal - insertion of the mass matrix and damping in the stiffness matrix for subsequent assembly global
	void MountFieldLoading();			//Assembles loadings of field
	void TransformMatrix();				//Assembles transformation matrix of coordinates
	void MountGlobal();					//Fills the contribution of the element in the matrices global
	void SaveLagrange();				//Saves variables in the gauss points useful for description lagrangiana updated
	void PreCalc();						//Pre-calculation of variables that and done a single time in the start
	double CalculateLength();			//Calculates the length indeformado
	void Zeros();						//Zeroes some matrices used in the calculations
	void EvaluateMassModal(double* v, double* alphai, double** Jr, double** Mr, double* br, double** matrixm);

	//Calculates the contributions inerciais for the analysis dynamics - inclui all the contributions for the way fraca and for the operator tangent
	void EvaluateInertialContributions(double* v, double(*a1)
		, double(*a2), double(*a3), double(*a4), double(*a5), double(*a6)
		, double* alphai, double* alphad, double* ui, double* ud, double* omegai
		, double* domegai, double* dui, double* ddui, double** Jr, double** Mr
		, double* br, double* dT, double** DdT);
	//Variables for funçõa generated in the AceGen
	double temp_v[2000];				//temporary variable for calculations internal
	double** pJr;						//Pointer double** - conversion of matrix
	double** pMr;						//Pointer double** - conversion of matrix
	Matrix *br;
	Matrix *DdT;
	Matrix *dT;
	double** pDdT;
	double* tempkin;
	
	double alpha1;						//Weight of the method of quadratura Gaussiana
	double length;						//Length indeformado of the element
	double jacobian;					//Jacobian 
	double alpha_escalar_delta;
	double g;
	double* N1;							//Functions of way and its derivatives
	double* N2;
	double* N3;
	double* dN1;
	double* dN2;
	double* dN3;	
	double* csi;						//Coordinate natural of the element isoparametrico in each gauss point
	//Variables saved in each gauss point for subsequent pós-processamento and facilidade to the lidar with Lag. Updated
	Matrix** N;							//Matrix of the functions of way
	Matrix** deltaN;					//Matrix with the derivatives of the functions of way and with the functions of way
	Matrix** alpha_delta;	
	Matrix** d_alpha_delta;
	Matrix** u_delta;
	Matrix** d_u_delta;
	Matrix** A_delta;
	Matrix** Q_delta;
	Matrix** Xi_delta;
	Matrix** d_A_delta;
	Matrix** d_Xi_delta;
	Matrix** d_z;
	Matrix** d_Z;
	Matrix** eta_r;
	Matrix** kappa_r;
	Matrix** epsilon_r;
	Matrix** sigma_r;
	Matrix** n_r;
	Matrix** m_r;
	Matrix** n;
	Matrix** m;
	Matrix* D;							
	Matrix* I3;
	Matrix* e3r;
	Matrix* B1;
	Matrix* Qtransp;
	Matrix* B2;
	Matrix* B2temp;
	Matrix* stiffness;								//Stiffness matrix
	Matrix* constitutive_stiffness;					//Stiffness matrix constitutiva
	Matrix* geometric_stiffness;					//Geometric stiffness matrix
	Matrix* loading_stiffness;						//Stiffness matrix of loading
	Matrix* mass;									//Mass matrix
	Matrix* mass_modal;								//Mass matrix for analysis modal
	Matrix* damping;								//Damping matrix
	Matrix* damping_modal;							//Damping matrix for analysis modal
	Matrix* rayleigh_damping;						//Damping matrix initial of the problem
	Matrix* i_loading;								//Vector of internal forces internal
	Matrix* e_loading;								//Vector of internal forces external
	Matrix* P_loading;								//Vector internal force unbalanced
	Matrix* inertial_loading;						//Vector of internal forces inerciais
	Matrix* damping_loading;						//Vector of internal forces of damping
	Matrix* transform;								//Transformation matrix of coordinates
	Matrix* transform3;								//Transformation matrix of coordinates 3x3
	Matrix* V_alpha_dz_n;
	Matrix* V_alpha_m;
	Matrix* d_V_dalpha_apha_m;
	Matrix* G_d_u_alpha;
	Matrix* G_d_u_alpha_transp;
	Matrix* G_alpha_alpha;
	Matrix* G_alpha_d_alpha;
	Matrix* G_alpha_d_alpha_transp;
	Matrix* B;
	Matrix* G;
	LagrangeSave* lag_save;							//For save the variables devidas when there is convergence (access through of the function SaveLagrange())

	//Variables internal for use in the dynamics
	Matrix** alpha_dot;
	Matrix** Xi_dot;
	Matrix* Mr;
	Matrix* Jr;
	Matrix** Mip;
	Matrix** Jip;
	Matrix** M;
	Matrix** Md1;

	//Variables for calculation of steps
	int t1, t2;
	double load_multiplier, l_factor, mult;

	Matrix* e3rg;
	double T0;	//PreTension
};

