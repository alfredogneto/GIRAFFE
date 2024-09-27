#pragma once
#include "Element.h"
#include "Matrix.h"

class Shell_1 :
	public Element
{
public:
	Shell_1();
	~Shell_1();
	bool Check();				//Checa inconsist�ncias no elemento para evitar erros de execu��o
	bool Read(FILE *f);
	void Write(FILE *f);
	void WriteResults(FILE *f);									//Escreve arquivo de resultados
	void WriteMonitor(FILE *f, bool first_record, double time);	//Escreve no monitor do elemento
	void WriteVTK_XMLBase(std::vector<float> *float_vector);
	void WriteVTK_XMLRender(FILE *f);
	void Mount();												//Monta elementos
	void MountElementLoads();									//Monta carregamentos associados ao elemento
	void MountMass();											//Monta a matriz de massa
	void MountMassModal();										//Monta a matriz de massa para realiza��o da analise modal
	void MountDampingModal();									//Monta a matriz de amortecimento para realiza��o da analise modal
	void MountDamping(bool update_rayleigh);					//Monta a matriz de amortecimento
	void MountDyn();											//Montagens - Newmark
	void MountDynModal();										//Montagens para analise modal - inser��o da matriz de massa e amortecimento na matriz de rigidez para posterior montagem global
	void TransformMatrix();										//Monta matriz de transforma��o de coordenadas
	void MountGlobal();											//Preenche a contribui��o do elemento nas matrizes globais
	void MountFieldLoads();										//Monta carregamentos de campo (ex: peso pr�prio)
	void MountShellSpecialLoads(int l_number);					//Monta carregamentos de press�o na casca

	void SaveLagrange();										//Salva variaveis nos pontos de Gauss uteis para descri��o lagrangiana atualizada
	void PreCalc();												//Pre-calculo de variaveis que e feito uma unica vez no inicio
	void Zeros();												//Zera matrizes

	//Fun��es de forma e suas derivadas
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

	Matrix** N;							//6x27//Matriz das fun��es de forma
	Matrix** deltaN;					//15x27//Matriz deltaN (com fun��es de forma e suas derivadas) 

	//////////////////////////////// Fun��es de forma para integra��o de grau 4 (6 pontos de Gauss)
	//Fun��es de forma
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

	Matrix** N4;							//6x27//Matriz das fun��es de forma
	Matrix** N4_u_x1;						//3x18//Matriz das fun��es de forma
	Matrix** N4_u_x2;						//3x18//Matriz das fun��es de forma
	double* w4;								//pesos
	Matrix** alpha_i4;
	////////////////////////////////////
	
	//Variaveis para descrever a cinematica
	Matrix** alpha_delta;
	Matrix** alpha_delta_x1;
	Matrix** alpha_delta_x2;
	Matrix** u_delta;
	Matrix** u_delta_x1;
	Matrix** u_delta_x2;
	Matrix* e1r;	//Dire��o 1 da casca (global)
	Matrix* e2r;	//Dire��o 2	da casca (global)
	Matrix* e3r;	//Dire��o 2	da casca (global)
	Matrix* e1rlocal;	//Dire��o 1 da casca (local)
	Matrix* e2rlocal;	//Dire��o 2	da casca (local)
	Matrix* e3rlocal;	//Dire��o 3	da casca (local)

	//Variaveis a serem salvas - Lag. Atualizado
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

	Matrix* stiffness;								//Matriz de rigidez
	Matrix* mass;									//Matriz de massa
	Matrix* mass_modal;								//Matriz de massa para analise modal
	Matrix* damping;								//Matriz de amortecimento
	Matrix* damping_modal;							//Matriz de amortecimento para analise modal
	Matrix* damping_loading;						//Esfor�os de amortecimento
	Matrix* rayleigh_damping;						//Matriz de amortecimento inicial do problema
	Matrix* i_loading;								//Vetor de esfor�os internos
	Matrix* inertial_loading;						//Vetor de esfor�os de inercia
	Matrix* P_loading;								//Vetor de esfor�os desbalanceados
	Matrix* e_loading;								//Vetor de esfor�os externos

	Matrix* transform;								//Matriz de transforma��o de coordenadas
	Matrix* transform3;								//Matriz de transforma��o de coordenadas 3x3

	double lambda, mu;								//Lame
	double stiff_drill;								//Rigidez ao drilling
	double area;									//area do elemento
	double alpha1;									//Pesos de integra��o

	double coef1, coef2, coef3;						//Coeficientes para uso na matriz de massa

	Matrix I3;

	Matrix** eta_r1;			//deforma��o em er1
	Matrix** eta_r2;			//deforma��o em er2
	Matrix** kappa_r1;			//rot/comp em er1
	Matrix** kappa_r2;			//rot/comp em er2
	Matrix** n_r1;				//for�a em er1
	Matrix** n_r2;				//for�a em er2
	Matrix** m_r1;				//momento em er1
	Matrix** m_r2;				//momento em er2

	Matrix** eta_r1_global;			//deforma��o em er1
	Matrix** eta_r2_global;			//deforma��o em er2
	Matrix** kappa_r1_global;		//rot/comp em er1
	Matrix** kappa_r2_global;		//rot/comp em er2
	Matrix** n_r1_global;			//for�a em er1
	Matrix** n_r2_global;			//for�a em er2
	Matrix** m_r1_global;			//momento em er1
	Matrix** m_r2_global;			//momento em er2

	//Variaveis internas para uso na dinamica
	Matrix** alpha_dot;
	Matrix** Xi_dot;
	Matrix** Mip;
	Matrix** Jip;
	Matrix** M;
	Matrix** Md1;
	Matrix v_ipp;//Estimativa da velocidade no instante posterior


	//Matrizes para transforma��o dos resultados
	Matrix n_1;
	Matrix n_2;
	Matrix m_1;
	Matrix m_2;

	Matrix e_1;
	Matrix e_2;
	Matrix k_1;
	Matrix k_2;

	//Calcula as contribui��es inerciais para a analise dinamica - inclui todas as contribui��es para a forma fraca e para o operador tangente
	void EvaluateInertialContributions(double* v, double(*a1)
		, double(*a2), double(*a3), double(*a4), double(*a5), double(*a6)
		, double* alphai, double* alphad, double* ui, double* ud, double* omegai
		, double* domegai, double* dui, double* ddui, double* e3r, double(*coef1)
		, double(*coef2), double* dT, double** DdT);

	void EvaluateMassModal(double* v, double* alphai, double* coef1, double* coef2, double* coef3, double* e3r, double** matrixm);

	//Variaveis para fun��o gerada no AceGen
	double temp_v[2000];				//variavel temporaria para calculos internos
	Matrix *DdT;
	Matrix *dT;
	double** pDdT;
	double* tempkin;

	//Variaveis para calculo de steps
	int t1, t2;
	double load_multiplier, l_factor, mult;

	//COMPOSITE MATERIAL
	Matrix* D_comp;				//Matriz Constitutiva Comp�sito
	double thick_comp;			//Espessura do Material Comp�sito
	double rho_comp;			//Rho do material comp�sito
	double stiff_drill_comp1;	//Rigidiz ao drilling na dire��o 1
	double stiff_drill_comp2;	//Rigidez ao drilling na dire��o 2
	bool composite_shell;		//Verifica se a casca e de material comp�sito


};