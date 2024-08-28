#pragma once
#include "Element.h"

class RigidBody_1 :
	public Element
{
public:
	RigidBody_1();
	~RigidBody_1();
	bool Check();													//Checa inconsist�ncias no elemento para evitar erros de execu��o
	bool Read(FILE *f);												//Leitura do arquivo de entrada
	void Write(FILE *f);											//Escrita do arquivo de saida
	void WriteResults(FILE *f);										//Escreve arquivo de resultados
	void WriteMonitor(FILE *f, bool first_record, double time);		//Escreve no monitor do elemento
	void WriteVTK_XMLBase(std::vector<float> *float_vector);		//Escrita VTK XML - usado no arquivo base
	void WriteVTK_XMLRender(FILE *f);								//Escrita VTK XML - usado no arquivo render
	void Mount();													//Monta contribui��es do elemento para rigidez e residuo
	void MountElementLoads();										//Monta carregamentos associados ao elemento
	void MountMass();												//Monta a matriz de massa
	void MountMassModal();											//Monta a matriz de massa para realiza��o da analise modal
	void MountDampingModal();										//Monta a matriz de amortecimento para realiza��o da analise modal
	void MountDamping(bool update_rayleigh);						//Monta a matriz de amortecimento
	void MountDyn();												//Montagens - Newmark
	void MountDynModal();											//Montagens para analise modal - inser��o da matriz de massa e amortecimento na matriz de rigidez para posterior montagem global
	void TransformMatrix();											//Monta matriz de transforma��o de coordenadas
	void MountGlobal();												//Preenche a contribui��o do elemento nas matrizes globais
	void Zeros();													//Zera algumas matrizes utilizadas nos calculos
	void SaveLagrange();											//Salva dados para configura��o convergida
	void PreCalc();													//Pre-calculo de variaveis que e feito uma unica vez no inicio
	void InertialContributions();				                    //Calcula contribui��es de inercia do Rigid Body
	void MountFieldLoading();											//Calcula a contribui��o dos esfor�os de campo (acelera��o da gravidade)

	//Variaveis do elemento
	int RB_data_ID;													//ID dos parametros do Rigid Body
	double v[1000];													//Variavel necessaria para o AceGen			
	double** Jr;													//Tensor de inercia - formato double**
	double* br;														//Vetor br - formato double*
	double** DdT;													//Operador tangente
	double* dT;														//Residuo
	double** Ddfield;												//Carregamento de campo linearizado
	double* dfield;													//Residuo carregamento de campo
	
	double*	T1;														//Energia cinetica translacional
	double*	T2;														//Energia cinetica rotacional
	double*	T3;														//Energia cinetica acoplamento
	double*	T;														//Energia cinetica
	double* L;														//Momentum Linear 
	double* magL;													//Momentum Linear 
	double* HG;													    //Momentum Angular around Barycenter
	double* magHG;												    //Momentum Angular around Barycenter

	
	//Variaveis cinematicas - declara��o
	double alphai[3];
	double ud[3];
	double alphad[3];
	double dui[3];
	double ddui[3];
	double omegai[3];
	double domegai[3];
		
};

