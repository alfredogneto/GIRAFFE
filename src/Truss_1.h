#pragma once
#include "Element.h"
#include "Matrix.h"

class Truss_1 :
	public Element
{
public:
	Truss_1();
	~Truss_1();

	bool Check();									//Checks inconsistencies in the element for avoid errors of execution
	bool Read(FILE *f);
	void Write(FILE *f);
	void WriteResults(FILE *f);						//Writes file of results
	void WriteMonitor(FILE *f, bool first_record, double time);	//Writes in the monitor of the element
	void WriteVTK_XMLBase(std::vector<float> *float_vector);
	void WriteVTK_XMLRender(FILE *f);

	void Mount();									//Assembles elements
	void MountMass();								//Assembles the mass matrix
	void MountAddedMass();							//Assembles the mass matrix adicional
	void MountAddedMassModal();						//Assembles the mass matrix adicional for analysis modal
	void MountMassModal();							//Assembles the mass matrix for performing of the analysis modal
	void MountDampingModal();						//Assembles the damping matrix for performing of the analysis modal
	void MountElementLoads();						//Assembles loadings associated to the element
	void MountDamping(bool update_rayleigh);		//Assembles the damping matrix
	void TransformMatrix();							//Assembles transformation matrix of coordinates
	void MountGlobal();								//Fills the contribution of the element in the matrices global
	void MountSeaCurrentLoading();					//Assembles loading of correnteza/damping hidro (Morison)
	void SaveLagrange();							//Saves variables in the gauss points useful for description lagrangiana updated
	
	void PreCalc();									//Pre-calculation of variables that and done a single time in the start
	void MountDyn();								//Assemblies - Newmark
	void MountDynModal();							//Assemblies for analysis modal - insertion of the mass matrix and damping in the stiffness matrix for subsequent assembly global
	void Zeros();									//Limpa the matrices internal of the element

	Matrix I3;
	double T;										//Normal force
	double L;										//Undeformed length
	double l;										//Deformed length
	double tau;										//Kirschhoff stress
	double A;										//cross sectional undeformed area
	double a;										//cross sectional deformed area
	double Vol;										//Undeformed volume
	double E, H, nu, rho, tau_y_0;					//material properties - local copies from data.materials
	double Ahydro;									//area for hydraulic evaluations purposes (buoyancy, Morison)
	double T0;	//PreTension
	
	Matrix c_internal_loads;						//internal loads contribution
	Matrix c_inertial_loads;						//inertial loads contribution
	Matrix c_damping_loads;							//damping loads contribution
	Matrix c_external_loads;						//external loads contribution
	
	Matrix c_damping_matrix;						//damping matrix
	Matrix c_stiffness_matrix;						//stiffness matrix
	Matrix c_mass_matrix;							//mass matrix
	Matrix c_external_loads_stiffness;				//Stiffness of loadings external
		
	Matrix c_rayleigh_damping;						//rayleigh damping matrix
	Matrix c_mass_modal;							//mass matrix for modal analysis
	Matrix normal;									//normal direction

	//Plasticity history variables
	double epsb;									//accumulated strain (to rule the hardening)
	double epsp;									//plastic strain
	double l_p;										//Plastic deformed length
	//stored copies (converged)
	double epsb_i;
	double epsp_i;
	double l_p_i;									

	//Returns the value of the Yieding function for given tau
	double YieldingFunction(double tau, double epsb);
	//Returns the sign of a number
	double sign(double number);

	//Flags - materials
	bool flag_plastic;
	bool flag_hydro;
	int pipe_sec;
	double rho_len;

	//Variables for calculation of steps
	int t1, t2;
	double load_multiplier, l_factor, mult;

	double* N1;		//Functions of way
	double* N2;
	Matrix** N;		//Matrix of the functions of way

	//Variables for functions geradas in the AceGen
	double C1t;							//auxiliary variables for calculation of the Morison
	double C1n;
	double rho_adt;						//auxiliary variables for calculation of the mass adicional
	double rho_adn;
	double temp_v[2000];				//temporary variable for calculations internal
	void EvaluateMorisonContributions(double* v, double(*a4), double
		(*a5), double(*a6), double(*C1t), double(*C1n), double* dUinf, double* xA, double* xB
		, double* uA, double* uB, double* ud, double* dui, double* ddui
		, double* force, double** stiffness1, double** stiffness2);
	void EvaluateAddedMassContributions(double* v, double(*a1), double(*a2), double(*a3)
		, double(*rhoadt), double(*rhoadn), double* xA, double* xB
		, double* uA, double* uB, double* ud, double* dui, double* ddui
		, double* force, double** stiffness1, double** stiffness2);
};		