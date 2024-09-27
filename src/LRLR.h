#pragma once
#include "Contact.h"
class Matrix;

////////Contato do tipo LRLR - "Line region to line region"///////////
/////Desenvolvido para trabalhar em conjunto com elementos Pipe_1/////
class LRLR :
	public Contact
{
public:
	bool Check();				//Checa inconsistências no elemento para evitar erros de execução
	LRLR();
	~LRLR();
	int n_LR1;
	int n_LR2;
	double mu;
	double ept;
	double epn;
	double pinball;
	double sum_Fx;												//Variaveis que armazenam somas vetoriais dos esforços (monitor)
	double sum_Fy;
	double sum_Fz;
	bool Read(FILE *f);
	void Write(FILE *f);
	void WriteResults(FILE *f);									//Escreve arquivo de resultados
	void WriteMonitor(FILE *f, bool first_record, double time);	//Escreve no monitor do contato
	void WriteVTK_XMLRender(FILE *f);
	void WriteVTK_XMLForces(FILE *f);
	void SaveLagrange();										//Salva variaveis para descrição lagrangiana atualizada
	void Mount();
	void MountGlobal();											//Preenche a contribuição do contato nas matrizes globais
	void Band(int* band_fixed, int* band_free);					//Calcula a banda gerada na matriz global pelo contato
	void PreCalc();												//Pre-calculo de variaveis que e feito uma unica vez no inicio
	void PinballCheck();										//Checks proximity between each beam from LR to LR
	void BeginStepCheck();											//checagem inicial do contato  - inicio de cada incremento
	void Alloc(int e_elements1, int e_elements2);				//Aloca na memória as variaveis que dependem do numero de elementos
	bool HaveErrors();											//Retorna 0 - não ha cruzamento 1 - ha cruzamento
	void MountDyn();												//Montagens - Newmark
	//Funções internas
	void FillNodes(int e_element1, int e_element2);				//Preenche as variaveis dos nós com valores atualizados
	int FindMinimumParameters(int i, int j);					//Determina csi1 e csi2 do ponto de minima distancia
	void EvaluateParameters(int i, int j);						//Calcula parametros em função de csi_1 e csi_2
	void CalculateLengthsAndTangents();							//Calcula o comprimento dos elementos //Calcula os vetores t1 e t2 dos elementos
	void PlotContactStatus(int i, int j);						//Plota caracteristicas do contato, independente de ter convergido ou não
	//Variaveis internas:
	int n_elements1;											//Numero de elementos da line region correspondente ao contato
	int n_elements2;											//Numero de elementos da line region correspondente ao contato
	int temp_element1;
	int temp_element2;
	int temp_node;
	double error;
	double tol_NR;												//Tolerência NR
	int max_it;													//Numero maximo de iterações
	double tol_ortho;											//Tolerancia de erro a ortogonalidade
	double dot_test;
	double gt1, gt2;
	double l1, l2;
	double gte1, gte2;
	double Fat_max;
	double t1t2;
	double tt1;
	double tt2;
	double m1;
	double m2;
	double Fat_try;
	double Fat;
	double Fat1;
	double Fat2;
	double delta_lambda, delta_lambda_1, delta_lambda_2;
	double d;
	double N_1_1, N_3_1, dN_1_1, dN_3_1, ddN_1_1, ddN_3_1;		//Variaveis relacionadas as funções de forma do elemento 1
	double N_1_2, N_3_2, dN_1_2, dN_3_2, ddN_1_2, ddN_3_2;		//Variaveis relacionadas as funções de forma do elemento 2
	double D_;
	double dot11;
	double dot22;
	double dot12;
	Matrix* mean_position1;										//Posição da media aritmetica dos nós de cada elemento - calculado para realizar o pinball search
	Matrix* mean_position2;										//Posição da media aritmetica dos nós de cada elemento - calculado para realizar o pinball search
	Matrix* gte;
	Matrix* node_1_1;											//Coordenadas do nó 1 - elemento 1
	Matrix* node_3_1;											//Coordenadas do nó 3 - elemento 1
	Matrix* node_1_2;											//Coordenadas do nó 1 - elemento 2
	Matrix* node_3_2;											//Coordenadas do nó 3 - elemento 2
	Matrix* N1;
	Matrix* dN1;
	Matrix* ddN1;
	Matrix* N2;
	Matrix* dN2;
	Matrix* ddN2;
	Matrix* x1;													//Posições nodais em um unico vetor:
	Matrix* x2;													//Posições nodais em um unico vetor:
	Matrix* z1;
	Matrix* z2;
	Matrix* t1ext;
	Matrix* t2ext;
	Matrix* I361;
	Matrix* I362;
	Matrix* Ntio;
	Matrix* G1;
	Matrix* G1b;
	Matrix* G2;
	Matrix* G2b;
	Matrix* G3;
	Matrix* G3b;
	Matrix* aux1;
	Matrix* aux2;
	Matrix* G4;
	Matrix* G4b;
	Matrix* aux3;
	Matrix* G5;
	Matrix* aux4;
	Matrix* aux5;
	Matrix* aux6;
	Matrix* aux7;
	Matrix* Q1;
	Matrix* Q2;
	Matrix* Q3;
	Matrix* Q4;
	Matrix* G6;
	Matrix* G6b;
	Matrix* G7;
	Matrix* G7b;
	Matrix* I3;													//Identidade de ordem 3
	Matrix* n;
	Matrix* dz1;		
	Matrix* dz2;		
	Matrix* ddz1;	
	Matrix* ddz2;	
	Matrix* N_ext;
	Matrix* dN_ext;
	Matrix* d1;
	Matrix* d2;
	Matrix* E1;
	Matrix* E2;
	Matrix* ST;													//Contribuições da parte do atrito ao operador tangente
	Matrix* STb;												//Contribuições da parte do atrito ao operador tangente
	Matrix* R;													//Residuo
	Matrix* A;													//Matriz Jacobiana - que tb. e a matriz A, uma vez convergido o NR
	Matrix* B;
	Matrix* C;
	Matrix* D;
	Matrix* E;
	Matrix* F;
	Matrix* G;
	Matrix* t1;													//Direções tangenciais dos elementos 1 e 2
	Matrix* t2;
	Matrix* delta_csi;
	//Variaveis dependentes do numero de elementos
	Matrix ***last_dif_pos;
	Matrix ***c_stiffness;										//Matriz de rigidez
	Matrix ***c_loading;										//Vetor de esforços externos
	Matrix*** z2z1;
	bool** first_mount;
	bool** flag_cross;
	bool** sticking;
	double** last_g_n;
	double** csi_1_0;
	double** csi_2_0;
	double** copy_csi_1_converged;
	double** copy_csi_2_converged;
	double** gt1s;
	double** gt2s;
	double** g_t1s_temp;
	double** g_t2s_temp;
	int** return_value;
	int** last_return_value;
	double** g_n;
	double* L1;													//Comprimento
	double* r1;													//Raio
	double* L2;													//Comprimento
	double* r2;													//Raio
	double** csi_1;
	double** csi_2;
	double** factor_1;
	double** factor_2;
	//Variaveis para verificação da função Band
	int temp_band_free;
	int temp_band_fixed;
	int lowest_free_global_DOF;
	int highest_free_global_DOF;
	int lowest_fixed_global_DOF;
	int highest_fixed_global_DOF;
	bool free_marked;
	bool fixed_marked;
};

