#pragma once
#include "Contact.h"
#include "SPContactData.h"
#include "SplineElementPair.h"
//Tipos de pares de superf�cie:
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
	double mus;		//Coef. de atrito est�tico
	double mud;		//Coef. de atrito din�mico
	double ept;		//penalidade tangencial
	double epn;		//penalidade normal
	double epn_n;	//pot�ncia da penalidade normal
	double ct;		//coef. de dissipa��o tangencial
	double cn;		//coef. de dissipa��o normal
	double pinball;	//pinball radius
	bool write_report;
	bool write_report_diverged;
	
	bool Check();												//Checa inconsist�ncias no elemento para evitar erros de execu��o
	bool Read(FILE *f);											//Leitura do arquivo de entrada
	void Write(FILE *f);										//Grava��o do arquivo de sa�da
	void WriteResults(FILE *f);									//Escreve arquivo de resultados
	void WriteMonitor(FILE *f, bool first_record, double time);	//Escreve no monitor do contato
	void WriteVTK_XMLRender(FILE *f);							//Plota superf�cies de contato
	void WriteVTK_XMLForces(FILE *f);							//Plota for�as de contato
	void SaveLagrange();										//Salva vari�veis para descri��o lagrangiana atualizada
	bool HaveErrors();											//Retorna 1 - h� algum erro, mesmo que tenha convergido. Ex: penetra��o excessiva
	void MountDyn();											//Montagens referentes � parte dependente de velocidades - din�mica
	void Mount();												//Montagem - contribui��es locais de contato
	void MountGlobal();								//Preenche a contribui��o do contato nas matrizes globais
	void Band(int* band_fixed, int* band_free);		//Calcula a banda gerada na matriz global pelo contato
	void PreCalc();									//Pr�-c�lculo de vari�veis que � feito uma �nica vez no in�cio
	void PinballCheck();							//Checks proximity between Surface to Surface
	void Alloc();									//Aloca na mem�ria
	void SetPairs();								//Atribui pares de superf�cies, de acordo com varredura feita entre os surface sets 1 e 2										
	void BeginStepCheck();											//checagem inicial do contato  - in�cio de cada incremento

	void AllocSpecific(int surface1_index, int surface2_index, int sol_index);				//Aloca��o espec�fica para o contato surface1, surface2
	void FreeSpecific(int surface1_index, int surface2_index, int sol_index);				//Libera��o espec�fica para o contato surface1, surface2
	void EvaluateNormalGap(int surf1_index, int surf2_index, int sol_index);				//Calcula matriz com gap normal entre surface1, surface2
	void ReportContact(int surface1_index, int surface2_index, int sol_index);				//Imprime na tela informa��es do par de contato surface1, surface2
	
	
	//Vari�veis internas
	int number_surfaces1;														//N�mero de superf�cies do SurfaceSet 1
	int number_surfaces2;														//N�mero de superf�cies do SurfaceSet 2
	bool*** alloc_control;														//Controle de aloca��o din�mica
	SPContactData*** cd;														//Dados das solu��es de contato
	int number_pointwise_interactions;			//N�mero m�ximo de intera��es pontuais por par superf�cie-superf�cie
	Matrix**** i_loading;						//Vetor de esfor�os internos
	double***** contact_stiffness;				//Matriz de rigidez
	int * DOFs_surfaces1;						//N�mero de graus de liberdade que cada uma das superf�cies no Surface Set possui - influ�ncia no vetor de res�duos e na regidez tangente
	int * DOFs_surfaces2;						//N�mero de graus de liberdade que cada uma das superf�cies no Surface Set possui - influ�ncia no vetor de res�duos e na regidez tangente			
	Matrix* I3;									//Matriz identidade de ordem 3
	SplineElementPair*** surf_pair;				//Pares de superf�cies para atribui��es de fun��es espec�ficas para cada par

	//Saving contact forces (for post-processing)
	Matrix**** fn;
	Matrix**** ft;

	//Tolerance for contact detection purposes (to avoid zero gaps at the begining of the simulation, for example)
	double pen_tol = 1e-14;
	bool specialLCP;
};

