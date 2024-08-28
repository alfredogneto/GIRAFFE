#pragma once
#include "Element.h"

class RigidBody_1 :
	public Element
{
public:
	RigidBody_1();
	~RigidBody_1();
	bool Check();													//Checa inconsistências no elemento para evitar erros de execução
	bool Read(FILE *f);												//Leitura do arquivo de entrada
	void Write(FILE *f);											//Escrita do arquivo de saída
	void WriteResults(FILE *f);										//Escreve arquivo de resultados
	void WriteMonitor(FILE *f, bool first_record, double time);		//Escreve no monitor do elemento
	void WriteVTK_XMLBase(std::vector<float> *float_vector);		//Escrita VTK XML - usado no arquivo base
	void WriteVTK_XMLRender(FILE *f);								//Escrita VTK XML - usado no arquivo render
	void Mount();													//Monta contribuições do elemento para rigidez e resíduo
	void MountElementLoads();										//Monta carregamentos associados ao elemento
	void MountMass();												//Monta a matriz de massa
	void MountMassModal();											//Monta a matriz de massa para realização da análise modal
	void MountDampingModal();										//Monta a matriz de amortecimento para realização da análise modal
	void MountDamping(bool update_rayleigh);						//Monta a matriz de amortecimento
	void MountDyn();												//Montagens - Newmark
	void MountDynModal();											//Montagens para análise modal - inserção da matriz de massa e amortecimento na matriz de rigidez para posterior montagem global
	void TransformMatrix();											//Monta matriz de transformação de coordenadas
	void MountGlobal();												//Preenche a contribuição do elemento nas matrizes globais
	void Zeros();													//Zera algumas matrizes utilizadas nos cálculos
	void SaveLagrange();											//Salva dados para configuração convergida
	void PreCalc();													//Pré-cálculo de variáveis que é feito uma única vez no início
	void InertialContributions();				                    //Calcula contribuições de inércia do Rigid Body
	void MountFieldLoading();											//Calcula a contribuição dos esforços de campo (aceleração da gravidade)

	//Variáveis do elemento
	int RB_data_ID;													//ID dos parâmetros do Rigid Body
	double v[1000];													//Variável necessária para o AceGen			
	double** Jr;													//Tensor de inércia - formato double**
	double* br;														//Vetor br - formato double*
	double** DdT;													//Operador tangente
	double* dT;														//Resíduo
	double** Ddfield;												//Carregamento de campo linearizado
	double* dfield;													//Resíduo carregamento de campo
	
	double*	T1;														//Energia cinética translacional
	double*	T2;														//Energia cinética rotacional
	double*	T3;														//Energia cinética acoplamento
	double*	T;														//Energia cinética
	double* L;														//Momentum Linear 
	double* magL;													//Momentum Linear 
	double* HG;													    //Momentum Angular around Barycenter
	double* magHG;												    //Momentum Angular around Barycenter

	
	//Variáveis cinemáticas - declaração
	double alphai[3];
	double ud[3];
	double alphad[3];
	double dui[3];
	double ddui[3];
	double omegai[3];
	double domegai[3];
		
};

