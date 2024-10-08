#pragma once
#include "Contact.h"

class SSContactData;
class SurfacePair;
class Matrix;

#include <typeinfo>

////////Contato do tipo SSSS - "Surface Set to Surface Set"///////////
class SSSS :
	public Contact
{
public:
	SSSS();
	~SSSS();
	int n_SS1;		//Surface set ID 1
	int n_SS2;		//Surface set ID 2
	double mus;		//Coef. de atrito estatico
	double mud;		//Coef. de atrito dinamico
	double ept;		//penalidade tangencial
	double epn;		//penalidade normal
	double epn0;	//expoente da lei normal de interface (opcional)
	double ct;		//coef. de dissipa��o tangencial
	double cn;		//coef. de dissipa��o normal
	double pinball;	//pinball radius
	bool write_report;
	bool write_report_diverged;
	
	bool Check();												//Checa inconsist�ncias no elemento para evitar erros de execu��o
	bool Read(FILE *f);											//Leitura do arquivo de entrada
	void Write(FILE *f);										//Grava��o do arquivo de saida
	void WriteResults(FILE *f);									//Escreve arquivo de resultados
	void WriteMonitor(FILE *f, bool first_record, double time);	//Escreve no monitor do contato
	void WriteVTK_XMLRender(FILE *f);							//Plota superficies de contato
	void WriteVTK_XMLForces(FILE *f);							//Plota for�as de contato
	void SaveLagrange();										//Salva variaveis para descri��o lagrangiana atualizada
	bool HaveErrors();											//Retorna 1 - ha algum erro, mesmo que tenha convergido. Ex: penetra��o excessiva
	void MountDyn();											//Montagens referentes a parte dependente de velocidades - dinamica
	void Mount();												//Montagem - contribui��es locais de contato
	void MountGlobal();								//Preenche a contribui��o do contato nas matrizes globais
	void Band(int* band_fixed, int* band_free);		//Calcula a banda gerada na matriz global pelo contato
	void PreCalc();									//Pre-calculo de variaveis que e feito uma unica vez no inicio
	void PinballCheck();							//Checks proximity between Surface to Surface
	void Alloc();									//Aloca na mem�ria
	void SetPairs();								//Atribui pares de superficies, de acordo com varredura feita entre os surface sets 1 e 2										
	void BeginStepCheck();											//checagem inicial do contato  - inicio de cada incremento

	void AllocSpecific(int surface1_index, int surface2_index, int sol_index);				//Aloca��o especifica para o contato surface1, surface2
	void FreeSpecific(int surface1_index, int surface2_index, int sol_index);				//Libera��o especifica para o contato surface1, surface2
	void EvaluateTangentialGap(int surf1_index, int surf2_index, int sol_index);		//Calcula matriz com gap tangencial entre surface1, surface2
	void EvaluateNormalGap(int surf1_index, int surf2_index, int sol_index);			//Calcula matriz com gap normal entre surface1, surface2
	void ReportContact(int surface1_index, int surface2_index, int sol_index);				//Imprime na tela informa��es do par de contato surface1, surface2
	
	
	//Variaveis internas
	int number_surfaces1;														//Numero de superficies do SurfaceSet 1
	int number_surfaces2;														//Numero de superficies do SurfaceSet 2
	bool*** alloc_control;														//Controle de aloca��o dinamica
	SSContactData*** cd;														//Dados das solu��es de contato
	int number_pointwise_interactions;			//Numero maximo de intera��es pontuais por par superficie-superficie
	Matrix**** i_loading;						//Vetor de esfor�os internos
	double***** contact_stiffness;				//Matriz de rigidez
	int * DOFs_surfaces1;						//Numero de graus de liberdade que cada uma das superficies no Surface Set possui - influ�ncia no vetor de residuos e na regidez tangente
	int * DOFs_surfaces2;						//Numero de graus de liberdade que cada uma das superficies no Surface Set possui - influ�ncia no vetor de residuos e na regidez tangente			
	Matrix* I3;									//Matriz identidade de ordem 3
	SurfacePair*** surf_pair;					//Pares de superficies para atribui��es de fun��es especificas para cada par

	//Saving contact forces (for post-processing)
	Matrix**** fn;
	Matrix**** ft;

	//Tolerance for contact detection purposes (to avoid zero gaps at the begining of the simulation, for example)
	double pen_tol = 1e-14;
	bool specialLCP;
};

