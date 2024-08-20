#pragma once
#include "Element.h"
#include "LagrangeSave.h"

class Pipe_1 :
	public Element
{
public:
	Pipe_1();
	~Pipe_1();
	bool Check();				//Checa inconsist�ncias no elemento para evitar erros de execu��o
	bool Read(FILE *f);
	void Write(FILE *f);
	void WriteResults(FILE *f);//Escreve arquivo de resultados
	void WriteVTK_XMLBase(std::vector<float> *float_vector);
	void WriteVTK_XMLRender(FILE *f);
	void WriteMonitor(FILE *f, bool first_record, double time);	//Escreve no monitor do elemento
	void Mount();						//Monta elementos
	void MountElementLoads();			//Monta carregamentos associados ao elemento
	void MountMass();					//Monta a matriz de massa
	void MountMassModal();				//Monta a matriz de massa para realiza��o da an�lise modal
	void MountDampingModal();			//Monta a matriz de amortecimento para realiza��o da an�lise modal
	void MountDamping(bool update_rayleigh);	//Monta a matriz de amortecimento
	void MountDyn();					//Montagens - Newmark
	void MountDynModal();				//Montagens para an�lise modal - inser��o da matriz de massa e amortecimento na matriz de rigidez para posterior montagem global
	void MountFieldLoading();			//Monta carregamentos de campo
	void MountSeaCurrentLoading();		//Monta a parte est�tica dos esfor�os de Morison - correnteza mar�tima
	void MountPipeSpecialLoads(int l_number);//Monta carregamentos de press�o interna/externa no tubo
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

	void EvaluateMorisonContributions(double* v, double(*a4), double
		(*a5), double(*a6), double(*C1t), double(*C1n), double* dU, double* e3r
		, double* alphai, double* alphad, double* ud, double* dui, double* ddui
		, double* force, double** stiffness);
	double rho_adt;	//massa adicional
	double rho_adn;	//massa adicional
	//Vari�veis para fun��a gerada no AceGen
	double temp_v[2000];				//vari�vel tempor�ria para c�lculos internos
	double** pJr;						//Ponteiro double** - convers�o de matriz
	double** pMr;						//Ponteiro double** - convers�o de matriz
	Matrix *br;
	Matrix *DdT;
	Matrix *dT;
	double** pDdT;
	Matrix* inertial_loading;			//Vetor de esfor�os internos
	Matrix* morison_loading;
	double** pL;						//Loading stiffness (ponteiro)

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
	Matrix* damping;								//Matriz de amortecimento
	Matrix* mass_modal;									//Matriz de massa
	Matrix* damping_modal;								//Matriz de amortecimento
	Matrix* damping_loading;						//Matriz de amortecimento
	Matrix* rayleigh_damping;						//Matriz de amortecimento inicial do problema
	Matrix* i_loading;								//Vetor de esfor�os internos
	Matrix* e_loading;								//Vetor de esfor�os externos
	Matrix* P_loading;								//Vetor esfor�o desbalanceado
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

	//Vari�veis internas para o esfor�o de correnteza mar�tima
	Matrix** e3ip;
	Matrix** zi;
	Matrix** vel;
	Matrix** velr;
	Matrix** element_vel;
	Matrix** ut;
	Matrix** un;
	Matrix** d_e3_d_alpha;
	Matrix** Lt;
	Matrix** Ln;
	Matrix** L_u_alpha;
	Matrix** f_current;
	Matrix** L;
	Matrix* e3rg;
	double signal_t;
	double Cdt;
	double Cdn;
	double Aext;
	double rho_f;
	double depth;

	Matrix** t_e;
	Matrix** n_e;
	Matrix** vtr;
	Matrix** vnr;
	Matrix** Mdt;
	Matrix** Mdn;
	Matrix** Md2;
	double Un_;
	double un_;
	double Ut_;
	double ut_;
	double C1t;
	double C1n;

	//Vari�veis internas para o esfor�o pipe load
	Matrix** kip;
	Matrix** temp_f;
	Matrix** temp_m;
	Matrix** temp_l;
	Matrix** Kip;
	Matrix** E3ip;
	Matrix** UpsilonN;
	Matrix* K1ua;
	Matrix* K1aa;
	Matrix* K2ua;
	Matrix* K2au;
	Matrix* Kext;
	Matrix* O1;
	double p0i;
	double p0e;
	double rhoi;
	double rhoe;
	double Aint;

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
	//Fun��es Pipe:
	double De();								//Retorna o di�metro externo
};

