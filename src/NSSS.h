#pragma once
#include "Contact.h"

class Matrix;
class NSContactData;

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
	double cn;		//coef. de dissipa��o normal
	double ct;		//coef. de dissipa��o tangencial
	double pinball;	//pinball radius
	double radius;

	bool Check();				//Checa inconsist�ncias no elemento para evitar erros de execu��o
	bool Read(FILE *f);
	void Write(FILE *f);
	void WriteResults(FILE *f);									//Escreve arquivo de resultados
	void WriteMonitor(FILE *f, bool first_record, double time);	//Escreve no monitor do contato
	void WriteVTK_XMLRender(FILE *f);							//Plota superficies de contato
	void WriteVTK_XMLForces(FILE *f);							//Plota for�as de contato
	void SaveLagrange();										//Salva variaveis para descri��o lagrangiana atualizada
	bool HaveErrors();											//Retorna 1 - ha algum erro, mesmo que tenha convergido. Ex: penetra��o excessiva
	void MountDyn();											//Montagens - Newmark
	void Mount();
	void MountGlobal();								//Preenche a contribui��o do contato nas matrizes globais
	void Band(int* band_fixed, int* band_free);		//Calcula a banda gerada na matriz global pelo contato
	void PreCalc();									//Pre-calculo de variaveis que e feito uma unica vez no inicio
	void PinballCheck();							//Checks proximity between each beam from LR to AS
	void Alloc();									//Aloca na mem�ria
	void BeginStepCheck();											//checagem inicial do contato  - inicio de cada incremento
	
	void AllocSpecific(int node_index, int surface_index, int sol_index);				//Aloca��o especifica para o contato node_index, surface_index
	void FreeSpecific (int node_index, int surface_index, int sol_index);				//Libera��o especifica para o contato node_index, surface_index
	void EvaluateTangentialGap(int node_index, int surf_index, int sol_index);			//Calcula matriz com gap tangencial

	void ReportContact(int node_index, int surf_index, int sol_index);					//Imprime na tela informa��es do par de contato
	int number_nodes;								//Numero de n�s do NodeSet
	int number_surfaces;							//Numero de superficies do SurfaceSet

	bool*** alloc_control;							//Controle de aloca��o dinamica
	
	//Variaveis internas
	int number_pointwise_interactions;			//Numero maximo de intera��es pontuais por par n�-superficie
	NSContactData*** cd;						//Dados das solu��es de contato
	Matrix*** xS_p;
	Matrix*** QS;
	Matrix*** alphaS;
	Matrix*** uS;
	Matrix**** i_loading;						//Vetor de esfor�os internos
	double***** contact_stiffness;				//Matriz de rigidez
	
	int * DOFs_surfaces;						//Numero de graus de liberdade que cada uma das superficies no Surface Set possui - influ�ncia no vetor de residuos e na regidez tangente			

	Matrix* I3;

	
};

