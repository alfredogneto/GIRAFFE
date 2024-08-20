#pragma once
#include "Element.h"
#include "LagrangeSave.h"

class Pipe_1 :
	public Element
{
public:
	Pipe_1();
	~Pipe_1();
	bool Check();				//Checa inconsistências no elemento para evitar erros de execução
	bool Read(FILE *f);
	void Write(FILE *f);
	void WriteResults(FILE *f);//Escreve arquivo de resultados
	void WriteVTK_XMLBase(std::vector<float> *float_vector);
	void WriteVTK_XMLRender(FILE *f);
	void WriteMonitor(FILE *f, bool first_record, double time);	//Escreve no monitor do elemento
	void Mount();						//Monta elementos
	void MountElementLoads();			//Monta carregamentos associados ao elemento
	void MountMass();					//Monta a matriz de massa
	void MountMassModal();				//Monta a matriz de massa para realização da análise modal
	void MountDampingModal();			//Monta a matriz de amortecimento para realização da análise modal
	void MountDamping(bool update_rayleigh);	//Monta a matriz de amortecimento
	void MountDyn();					//Montagens - Newmark
	void MountDynModal();				//Montagens para análise modal - inserção da matriz de massa e amortecimento na matriz de rigidez para posterior montagem global
	void MountFieldLoading();			//Monta carregamentos de campo
	void MountSeaCurrentLoading();		//Monta a parte estática dos esforços de Morison - correnteza marítima
	void MountPipeSpecialLoads(int l_number);//Monta carregamentos de pressão interna/externa no tubo
	void TransformMatrix();				//Monta matriz de transformação de coordenadas
	void MountGlobal();					//Preenche a contribuição do elemento nas matrizes globais
	void SaveLagrange();				//Salva variáveis nos pontos de Gauss úteis para descrição lagrangiana atualizada
	void PreCalc();						//Pré-cálculo de variáveis que é feito uma única vez no início
	double CalculateLength();			//Calcula o comprimento indeformado
	void Zeros();						//Zera algumas matrizes utilizadas nos cálculos

	void EvaluateMassModal(double* v, double* alphai, double** Jr, double** Mr, double* br, double** matrixm);
	
	//Calcula as contribuições inerciais para a análise dinâmica - inclui todas as contribuições para a forma fraca e para o operador tangente
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
	//Variáveis para funçõa gerada no AceGen
	double temp_v[2000];				//variável temporária para cálculos internos
	double** pJr;						//Ponteiro double** - conversão de matriz
	double** pMr;						//Ponteiro double** - conversão de matriz
	Matrix *br;
	Matrix *DdT;
	Matrix *dT;
	double** pDdT;
	Matrix* inertial_loading;			//Vetor de esforços internos
	Matrix* morison_loading;
	double** pL;						//Loading stiffness (ponteiro)

	double alpha1;						//Peso do método de quadratura Gaussiana
	double length;						//Comprimento indeformado do elemento
	double jacobian;					//Jacobiano 
	double alpha_escalar_delta;
	double g;
	double* N1;							//Funções de forma e suas derivadas
	double* N2;
	double* N3;
	double* dN1;
	double* dN2;
	double* dN3;
	double* csi;						//Coordenada natural do elemento isoparamétrico em cada ponto de Gauss
	//Variáveis salvas em cada ponto de Gauss para posterior pós-processamento e facilidade ao lidar com Lag. Atualizado
	Matrix** N;							//Matriz das funções de forma
	Matrix** deltaN;					//Matriz com as derivadas das funções de forma e com as funções de forma
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
	Matrix* geometric_stiffness;					//Matriz de rigidez geométrica
	Matrix* loading_stiffness;						//Matriz de rigidez de carregamento
	Matrix* mass;									//Matriz de massa
	Matrix* damping;								//Matriz de amortecimento
	Matrix* mass_modal;									//Matriz de massa
	Matrix* damping_modal;								//Matriz de amortecimento
	Matrix* damping_loading;						//Matriz de amortecimento
	Matrix* rayleigh_damping;						//Matriz de amortecimento inicial do problema
	Matrix* i_loading;								//Vetor de esforços internos
	Matrix* e_loading;								//Vetor de esforços externos
	Matrix* P_loading;								//Vetor esforço desbalanceado
	Matrix* transform;								//Matriz de transformação de coordenadas
	Matrix* transform3;								//Matriz de transformação de coordenadas 3x3
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
	LagrangeSave* lag_save;							//Para salvar as variáveis devidas quando há convergência (acesso através da função SaveLagrange())

	//Variáveis internas para o esforço de correnteza marítima
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

	//Variáveis internas para o esforço pipe load
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

	//Variáveis internas para uso na dinâmica
	Matrix** alpha_dot;
	Matrix** Xi_dot;
	Matrix* Mr;
	Matrix* Jr;
	Matrix** Mip;
	Matrix** Jip;
	Matrix** M;
	Matrix** Md1;

	//Variáveis para cálculo de steps
	int t1, t2;
	double load_multiplier, l_factor, mult;
	//Funções Pipe:
	double De();								//Retorna o diâmetro externo
};

