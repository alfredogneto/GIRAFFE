#pragma once
#include <stdio.h>

class ConvergenceCriteria
{
public:
	ConvergenceCriteria();
	~ConvergenceCriteria();

	//Fun��es para estabelecer e checar converg�ncia
	void EstablishResidualCriteria();
	bool CheckResidualConvergence();
	bool CheckGLConvergence();

	bool Read(FILE *f);			//Leitura
	void Write(FILE *f);		//Grava��o

	double force_tol;			//Fator aplicado � norma do vetor de for�as do sistema para estabelecer crit�rio de converg�ncia
	double moment_tol;			//Fator aplicado � norma do vetor de momentos do sistema para estabelecer crit�rio de converg�ncia
	
	double force_ref;			//Valor de refer�ncia para avalia��o do crit�rio de parada de for�as
	int force_n_times;			//N�mero de valores utilizado no c�lculo da m�dia do valor de refer�ncia
	double moment_ref;			//Valor de refer�ncia para avalia��o do crit�rio de parada de momentos
	int moment_n_times;			//N�mero de valores utilizado no c�lculo da m�dia do valor de refer�ncia


	double force_min;			//Valor de for�a considerado pequeno - usado quando h� somente for�as nulas - valor dimensional
	double moment_min;			//Valor de momento considerado pequeno - usado quando h� somente momentos nulos - valor dimensional
	double constraint_min;		//Valor de erro de constraint, considerado pequeno - sempre � utilizado quando h� joints

	double disp_tol;			//Fator aplicado � norma do vetor de deslocamentos incrementais do sistema para estabelecer crit�rio de converg�ncia
	double rot_tol;				//Fator aplicado � norma do vetor de rota��es incrementais do sistema para estabelecer crit�rio de converg�ncia
	double lag_tol;				//Fator aplicado � norma do vetor de multiplicadores de lagrange do sistema para estabelecer crit�rio de converg�ncia

	double disp_min;			//Valor de deslocamento considerado pequeno - usado quando h� somente deslocamentos nulos - valor dimensional
	double rot_min;				//Valor de rota��o considerado pequeno - usado quando h� somente rota��es nulos - valor dimensional
	double lag_min;				//Valor de multiplicador de lagrange considerado pequeno - usado quando h� somente multiplicadores de lagrange nulos - valor dimensional

	double divergence_ref;		//Valor bastante elevado que indica diverg�ncia do modelo
	bool diverged;				//Vari�vel booleada que indica diverg�ncia

	double disp_criterion;		//Crit�rio de converg�ncia para norma de deslocamento - calculado automaticamente
	double rot_criterion;		//Crit�rio de converg�ncia para norma de rota��o - calculado automaticamente
	double lag_criterion;		//Crit�rio de converg�ncia para norma de multiplicador de lagrange - calculado automaticamente
	double force_criterion;		//Crit�rio de converg�ncia para norma de for�a - calculado automaticamente
	double moment_criterion;	//Crit�rio de converg�ncia para norma de momento - calculado automaticamente	
	double constraint_criterion;//Crit�rio de converg�ncia para norma de constraint - calculado automaticamente

	bool NaNDetector(double &value);

	void PlotDivergenceReport();//Imprime na tela informa��es para auxiliar o porqu� da diverg�ncia

	//N�s ou SpecialConstraint associados ao m�ximo erro
	int super_node_disp;
	int node_disp;
	int node_rot;
	int sc_lag;
	//N�s ou SpecialConstraint associados ao m�ximo erro
	int super_node_force;
	int node_force;
	int node_moment;
	int sc_constraint;

	int n_conv_evaluated;	//Contador do n�mero de vezes que o crit�rio foi estabelecido.

};

