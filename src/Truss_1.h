#pragma once
#include "Element.h"
class Truss_1 :
	public Element
{
public:
	Truss_1();
	~Truss_1();

	bool Check();									//Checa inconsist�ncias no elemento para evitar erros de execu��o
	bool Read(FILE *f);
	void Write(FILE *f);
	void WriteResults(FILE *f);						//Escreve arquivo de resultados
	void WriteMonitor(FILE *f, bool first_record, double time);	//Escreve no monitor do elemento
	void WriteVTK_XMLBase(std::vector<float> *float_vector);
	void WriteVTK_XMLRender(FILE *f);

	void Mount();									//Monta elementos
	void MountMass();								//Monta a matriz de massa
	void MountAddedMass();							//Monta a matriz de massa adicional
	void MountAddedMassModal();						//Monta a matriz de massa adicional para an�lise modal
	void MountMassModal();							//Monta a matriz de massa para realiza��o da an�lise modal
	void MountDampingModal();						//Monta a matriz de amortecimento para realiza��o da an�lise modal
	void MountElementLoads();						//Monta carregamentos associados ao elemento
	void MountDamping(bool update_rayleigh);		//Monta a matriz de amortecimento
	void TransformMatrix();							//Monta matriz de transforma��o de coordenadas
	void MountGlobal();								//Preenche a contribui��o do elemento nas matrizes globais
	void MountSeaCurrentLoading();					//Monta carregamento de correnteza/amortecimento hidro (Morison)
	void SaveLagrange();							//Salva vari�veis nos pontos de Gauss �teis para descri��o lagrangiana atualizada
	
	void PreCalc();									//Pr�-c�lculo de vari�veis que � feito uma �nica vez no in�cio
	void MountDyn();								//Montagens - Newmark
	void MountDynModal();							//Montagens para an�lise modal - inser��o da matriz de massa e amortecimento na matriz de rigidez para posterior montagem global
	void Zeros();									//Limpa as matrizes internas do elemento

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
	Matrix c_external_loads_stiffness;				//Rigidez de carregamentos externos
		
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

	//Vari�veis para c�lculo de steps
	int t1, t2;
	double load_multiplier, l_factor, mult;

	double* N1;		//Fun��es de forma
	double* N2;
	Matrix** N;		//Matriz das fun��es de forma

	//Vari�veis para fun��es geradas no AceGen
	double C1t;							//vari�veis auxiliares para c�lculo do Morison
	double C1n;
	double rho_adt;						//vari�veis auxiliares para c�lculo da massa adicional
	double rho_adn;
	double temp_v[2000];				//vari�vel tempor�ria para c�lculos internos
	void EvaluateMorisonContributions(double* v, double(*a4), double
		(*a5), double(*a6), double(*C1t), double(*C1n), double* dUinf, double* xA, double* xB
		, double* uA, double* uB, double* ud, double* dui, double* ddui
		, double* force, double** stiffness1, double** stiffness2);
	void EvaluateAddedMassContributions(double* v, double(*a1), double(*a2), double(*a3)
		, double(*rhoadt), double(*rhoadn), double* xA, double* xB
		, double* uA, double* uB, double* ud, double* dui, double* ddui
		, double* force, double** stiffness1, double** stiffness2);
};		