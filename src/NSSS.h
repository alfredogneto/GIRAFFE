#pragma once
#include "Contact.h"
#include "NSContactData.h"
#include <typeinfo>

////////Contato do tipo NSSS - "Node Set to Surface Set"///////////
class NSSS:
	public Contact
{
public:
	NSSS();
	~NSSS();
	int n_NS;		//Node set ID
	int n_SS;		//Surface set ID
	double mu;		//Coef. de atrito
	double ept;		//penalidade tangencial
	double epn;		//penalidade normal
	double cn;		//coef. de dissipação normal
	double ct;		//coef. de dissipação tangencial
	double pinball;	//pinball radius
	double radius;

	bool Check();				//Checa inconsistências no elemento para evitar erros de execução
	bool Read(FILE *f);
	void Write(FILE *f);
	void WriteResults(FILE *f);									//Escreve arquivo de resultados
	void WriteMonitor(FILE *f, bool first_record, double time);	//Escreve no monitor do contato
	void WriteVTK_XMLRender(FILE *f);							//Plota superfícies de contato
	void WriteVTK_XMLForces(FILE *f);							//Plota forças de contato
	void SaveLagrange();										//Salva variáveis para descrição lagrangiana atualizada
	bool HaveErrors();											//Retorna 1 - há algum erro, mesmo que tenha convergido. Ex: penetração excessiva
	void MountDyn();											//Montagens - Newmark
	void Mount();
	void MountGlobal();								//Preenche a contribuição do contato nas matrizes globais
	void Band(int* band_fixed, int* band_free);		//Calcula a banda gerada na matriz global pelo contato
	void PreCalc();									//Pré-cálculo de variáveis que é feito uma única vez no início
	void PinballCheck();							//Checks proximity between each beam from LR to AS
	void Alloc();									//Aloca na memória
	void BeginStepCheck();											//checagem inicial do contato  - início de cada incremento
	
	void AllocSpecific(int node_index, int surface_index, int sol_index);				//Alocação específica para o contato node_index, surface_index
	void FreeSpecific (int node_index, int surface_index, int sol_index);				//Liberação específica para o contato node_index, surface_index
	void EvaluateTangentialGap(int node_index, int surf_index, int sol_index);			//Calcula matriz com gap tangencial

	void ReportContact(int node_index, int surf_index, int sol_index);					//Imprime na tela informações do par de contato
	int number_nodes;								//Número de nós do NodeSet
	int number_surfaces;							//Número de superfícies do SurfaceSet

	bool*** alloc_control;							//Controle de alocação dinâmica
	
	//Variáveis internas
	int number_pointwise_interactions;			//Número máximo de interações pontuais por par nó-superfície
	NSContactData*** cd;						//Dados das soluções de contato
	Matrix*** xS_p;
	Matrix*** QS;
	Matrix*** alphaS;
	Matrix*** uS;
	Matrix**** i_loading;						//Vetor de esforços internos
	double***** contact_stiffness;				//Matriz de rigidez
	
	int * DOFs_surfaces;						//Número de graus de liberdade que cada uma das superfícies no Surface Set possui - influência no vetor de resíduos e na regidez tangente			

	Matrix* I3;

	
};

