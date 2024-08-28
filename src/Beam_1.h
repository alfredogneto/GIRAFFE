#pragma once
#include "Element.h"

class LagrangeSave;

class Beam_1 :
	public Element
{
public:
	Beam_1(void);
	~Beam_1(void);
	bool Check();						//Checa inconsist�ncias no elemento para evitar erros de execu��o
	bool Read(FILE *f);
	void Write(FILE *f);
	void WriteResults(FILE *f);			//Escreve arquivo de resultados
	void WriteVTK_XMLBase(std::vector<float> *float_vector);
	void WriteVTK_XMLRender(FILE *f);
	void WriteMonitor(FILE *f, bool first_record, double time);	//Escreve no monitor do elemento
	void Mount();						//Monta elementos
	void MountElementLoads();			//Monta carregamentos associados ao elemento
	void MountMass();					//Monta a matriz de massa
	void MountMassModal();				//Monta a matriz de massa para realiza��o da an�lise modal
	void MountDampingModal();			//Monta a matriz de amortecimento para realiza��o da an�lise modal
	void MountDamping(bool update_rayleigh);				//Monta a matriz de amortecimento
	void MountDyn() ;					//Montagens - Newmark
	void MountDynModal();				//Montagens para an�lise modal - inser��o da matriz de massa e amortecimento na matriz de rigidez para posterior montagem global
	void MountFieldLoading();			//Monta carregamentos de campo
	void TransformMatrix();				//Monta matriz de transforma��o de coordenadas
	void MountGlobal();					//Preenche a contribui��o do elemento nas matrizes globais
	void SaveLagrange();				//Salva vari�veis nos pontos de Gauss �teis para descri��o lagrangiana atualizada
	void PreCalc();						//Pr�-c�lculo de vari�veis que � feito uma �nica vez no in�cio
	double CalculateLength();			//Calcula o comprimento indeformado
	void Zeros();						//Zera algumas matrizes utilizadas nos c�lculos
	void EvaluateMassModal(double* v, double* alphai, double** Jr, double** Mr, double* br, double** matrixm);

	//Calcula as contribui��es inerciais para a an�lise din�mica - inclui todas as contribui��es para a forma fraca e para o operador tangente
	void EvaluateInertialContributions(double* v, double(*a1)
		, double(*a2), double(*a3), double(*a4), double(*a5), double(*a6)
		, double* alphai, double* alphad, double* ui, double* ud, double* omegai
		, double* domegai, double* dui, double* ddui, double** Jr, double** Mr
		, double* br, double* dT, double** DdT);
	//Vari�veis para fun��a gerada no AceGen
	double temp_v[2000];				//vari�vel tempor�ria para c�lculos internos
	double** pJr;						//Ponteiro double** - convers�o de matriz
	double** pMr;						//Ponteiro double** - convers�o de matriz
	Matrix *br;
	Matrix *DdT;
	Matrix *dT;
	double** pDdT;
	double* tempkin;
	
	double alpha1;						//Peso do m�todo de quadratura Gaussiana
	double length;						//Comprimento indeformado do elemento
	double jacobian;					//Jacobiano 
	double alpha_escalar_delta;
	double g;
	double* N1;							//Fun��es de forma e suas derivadas
	double* N2;
	double* N3;
	double* dN1;
	double* dN2;
	double* dN3;	
	double* csi;						//Coordenada natural do elemento isoparam�trico em cada ponto de Gauss
	//Vari�veis salvas em cada ponto de Gauss para posterior p�s-processamento e facilidade ao lidar com Lag. Atualizado
	Matrix** N;							//Matriz das fun��es de forma
	Matrix** deltaN;					//Matriz com as derivadas das fun��es de forma e com as fun��es de forma
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
	Matrix* stiffness;								//Matriz de rigidez
	Matrix* constitutive_stiffness;					//Matriz de rigidez constitutiva
	Matrix* geometric_stiffness;					//Matriz de rigidez geom�trica
	Matrix* loading_stiffness;						//Matriz de rigidez de carregamento
	Matrix* mass;									//Matriz de massa
	Matrix* mass_modal;								//Matriz de massa para an�lise modal
	Matrix* damping;								//Matriz de amortecimento
	Matrix* damping_modal;							//Matriz de amortecimento para an�lise modal
	Matrix* rayleigh_damping;						//Matriz de amortecimento inicial do problema
	Matrix* i_loading;								//Vetor de esfor�os internos
	Matrix* e_loading;								//Vetor de esfor�os externos
	Matrix* P_loading;								//Vetor esfor�o desbalanceado
	Matrix* inertial_loading;						//Vetor de esfor�os inerciais
	Matrix* damping_loading;						//Vetor de esfor�os de amortecimento
	Matrix* transform;								//Matriz de transforma��o de coordenadas
	Matrix* transform3;								//Matriz de transforma��o de coordenadas 3x3
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
	LagrangeSave* lag_save;							//Para salvar as vari�veis devidas quando h� converg�ncia (acesso atrav�s da fun��o SaveLagrange())

	//Vari�veis internas para uso na din�mica
	Matrix** alpha_dot;
	Matrix** Xi_dot;
	Matrix* Mr;
	Matrix* Jr;
	Matrix** Mip;
	Matrix** Jip;
	Matrix** M;
	Matrix** Md1;

	//Vari�veis para c�lculo de steps
	int t1, t2;
	double load_multiplier, l_factor, mult;

	Matrix* e3rg;
	double T0;	//PreTension
};

