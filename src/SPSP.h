#pragma once
#include "Contact.h"
#include "SPContactData.h"
#include "SplineElementPair.h"
//Tipos de pares de superfície:
#include "SplineElement_SplineElement.h"

#include <typeinfo>

////////Contato do tipo SPSP - "Spline to Spline"///////////
class SPSP :
	public Contact
{
public:
	SPSP();
	~SPSP();
	int n_SP1;		//Spline ID 1
	int n_SP2;		//Spline ID 2
	double mus;		//Coef. de atrito estático
	double mud;		//Coef. de atrito dinâmico
	double ept;		//penalidade tangencial
	double epn;		//penalidade normal
	double epn_n;	//potência da penalidade normal
	double ct;		//coef. de dissipação tangencial
	double cn;		//coef. de dissipação normal
	double pinball;	//pinball radius
	bool write_report;
	bool write_report_diverged;
	
	bool Check();												//Checa inconsistências no elemento para evitar erros de execução
	bool Read(FILE *f);											//Leitura do arquivo de entrada
	void Write(FILE *f);										//Gravação do arquivo de saída
	void WriteResults(FILE *f);									//Escreve arquivo de resultados
	void WriteMonitor(FILE *f, bool first_record, double time);	//Escreve no monitor do contato
	void WriteVTK_XMLRender(FILE *f);							//Plota superfícies de contato
	void WriteVTK_XMLForces(FILE *f);							//Plota forças de contato
	void SaveLagrange();										//Salva variáveis para descrição lagrangiana atualizada
	bool HaveErrors();											//Retorna 1 - há algum erro, mesmo que tenha convergido. Ex: penetração excessiva
	void MountDyn();											//Montagens referentes à parte dependente de velocidades - dinâmica
	void Mount();												//Montagem - contribuições locais de contato
	void MountGlobal();								//Preenche a contribuição do contato nas matrizes globais
	void Band(int* band_fixed, int* band_free);		//Calcula a banda gerada na matriz global pelo contato
	void PreCalc();									//Pré-cálculo de variáveis que é feito uma única vez no início
	void PinballCheck();							//Checks proximity between Surface to Surface
	void Alloc();									//Aloca na memória
	void SetPairs();								//Atribui pares de superfícies, de acordo com varredura feita entre os surface sets 1 e 2										
	void BeginStepCheck();											//checagem inicial do contato  - início de cada incremento

	void AllocSpecific(int surface1_index, int surface2_index, int sol_index);				//Alocação específica para o contato surface1, surface2
	void FreeSpecific(int surface1_index, int surface2_index, int sol_index);				//Liberação específica para o contato surface1, surface2
	void EvaluateNormalGap(int surf1_index, int surf2_index, int sol_index);				//Calcula matriz com gap normal entre surface1, surface2
	void ReportContact(int surface1_index, int surface2_index, int sol_index);				//Imprime na tela informações do par de contato surface1, surface2
	
	
	//Variáveis internas
	int number_surfaces1;														//Número de superfícies do SurfaceSet 1
	int number_surfaces2;														//Número de superfícies do SurfaceSet 2
	bool*** alloc_control;														//Controle de alocação dinâmica
	SPContactData*** cd;														//Dados das soluções de contato
	int number_pointwise_interactions;			//Número máximo de interações pontuais por par superfície-superfície
	Matrix**** i_loading;						//Vetor de esforços internos
	double***** contact_stiffness;				//Matriz de rigidez
	int * DOFs_surfaces1;						//Número de graus de liberdade que cada uma das superfícies no Surface Set possui - influência no vetor de resíduos e na regidez tangente
	int * DOFs_surfaces2;						//Número de graus de liberdade que cada uma das superfícies no Surface Set possui - influência no vetor de resíduos e na regidez tangente			
	Matrix* I3;									//Matriz identidade de ordem 3
	SplineElementPair*** surf_pair;				//Pares de superfícies para atribuições de funções específicas para cada par

	//Saving contact forces (for post-processing)
	Matrix**** fn;
	Matrix**** ft;

	//Tolerance for contact detection purposes (to avoid zero gaps at the begining of the simulation, for example)
	double pen_tol = 1e-14;
	bool specialLCP;
};

