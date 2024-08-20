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
	double cn;		//coef. de dissipa��o normal
	double ct;		//coef. de dissipa��o tangencial
	double pinball;	//pinball radius
	double radius;

	bool Check();				//Checa inconsist�ncias no elemento para evitar erros de execu��o
	bool Read(FILE *f);
	void Write(FILE *f);
	void WriteResults(FILE *f);									//Escreve arquivo de resultados
	void WriteMonitor(FILE *f, bool first_record, double time);	//Escreve no monitor do contato
	void WriteVTK_XMLRender(FILE *f);							//Plota superf�cies de contato
	void WriteVTK_XMLForces(FILE *f);							//Plota for�as de contato
	void SaveLagrange();										//Salva vari�veis para descri��o lagrangiana atualizada
	bool HaveErrors();											//Retorna 1 - h� algum erro, mesmo que tenha convergido. Ex: penetra��o excessiva
	void MountDyn();											//Montagens - Newmark
	void Mount();
	void MountGlobal();								//Preenche a contribui��o do contato nas matrizes globais
	void Band(int* band_fixed, int* band_free);		//Calcula a banda gerada na matriz global pelo contato
	void PreCalc();									//Pr�-c�lculo de vari�veis que � feito uma �nica vez no in�cio
	void PinballCheck();							//Checks proximity between each beam from LR to AS
	void Alloc();									//Aloca na mem�ria
	void BeginStepCheck();											//checagem inicial do contato  - in�cio de cada incremento
	
	void AllocSpecific(int node_index, int surface_index, int sol_index);				//Aloca��o espec�fica para o contato node_index, surface_index
	void FreeSpecific (int node_index, int surface_index, int sol_index);				//Libera��o espec�fica para o contato node_index, surface_index
	void EvaluateTangentialGap(int node_index, int surf_index, int sol_index);			//Calcula matriz com gap tangencial

	void ReportContact(int node_index, int surf_index, int sol_index);					//Imprime na tela informa��es do par de contato
	int number_nodes;								//N�mero de n�s do NodeSet
	int number_surfaces;							//N�mero de superf�cies do SurfaceSet

	bool*** alloc_control;							//Controle de aloca��o din�mica
	
	//Vari�veis internas
	int number_pointwise_interactions;			//N�mero m�ximo de intera��es pontuais por par n�-superf�cie
	NSContactData*** cd;						//Dados das solu��es de contato
	Matrix*** xS_p;
	Matrix*** QS;
	Matrix*** alphaS;
	Matrix*** uS;
	Matrix**** i_loading;						//Vetor de esfor�os internos
	double***** contact_stiffness;				//Matriz de rigidez
	
	int * DOFs_surfaces;						//N�mero de graus de liberdade que cada uma das superf�cies no Surface Set possui - influ�ncia no vetor de res�duos e na regidez tangente			

	Matrix* I3;

	
};

