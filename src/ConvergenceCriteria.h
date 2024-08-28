#pragma once
#include <stdio.h>

class ConvergenceCriteria
{
public:
	ConvergenceCriteria();
	~ConvergenceCriteria();

	//Funções para estabelecer e checar convergência
	void EstablishResidualCriteria();
	bool CheckResidualConvergence();
	bool CheckGLConvergence();

	bool Read(FILE *f);			//Leitura
	void Write(FILE *f);		//Gravação

	double force_tol;			//Fator aplicado à norma do vetor de forças do sistema para estabelecer critério de convergência
	double moment_tol;			//Fator aplicado à norma do vetor de momentos do sistema para estabelecer critério de convergência
	
	double force_ref;			//Valor de referência para avaliação do critério de parada de forças
	int force_n_times;			//Número de valores utilizado no cálculo da média do valor de referência
	double moment_ref;			//Valor de referência para avaliação do critério de parada de momentos
	int moment_n_times;			//Número de valores utilizado no cálculo da média do valor de referência


	double force_min;			//Valor de força considerado pequeno - usado quando há somente forças nulas - valor dimensional
	double moment_min;			//Valor de momento considerado pequeno - usado quando há somente momentos nulos - valor dimensional
	double constraint_min;		//Valor de erro de constraint, considerado pequeno - sempre é utilizado quando há joints

	double disp_tol;			//Fator aplicado à norma do vetor de deslocamentos incrementais do sistema para estabelecer critério de convergência
	double rot_tol;				//Fator aplicado à norma do vetor de rotações incrementais do sistema para estabelecer critério de convergência
	double lag_tol;				//Fator aplicado à norma do vetor de multiplicadores de lagrange do sistema para estabelecer critério de convergência

	double disp_min;			//Valor de deslocamento considerado pequeno - usado quando há somente deslocamentos nulos - valor dimensional
	double rot_min;				//Valor de rotação considerado pequeno - usado quando há somente rotações nulos - valor dimensional
	double lag_min;				//Valor de multiplicador de lagrange considerado pequeno - usado quando há somente multiplicadores de lagrange nulos - valor dimensional

	double divergence_ref;		//Valor bastante elevado que indica divergência do modelo
	bool diverged;				//Variável booleada que indica divergência

	double disp_criterion;		//Critério de convergência para norma de deslocamento - calculado automaticamente
	double rot_criterion;		//Critério de convergência para norma de rotação - calculado automaticamente
	double lag_criterion;		//Critério de convergência para norma de multiplicador de lagrange - calculado automaticamente
	double force_criterion;		//Critério de convergência para norma de força - calculado automaticamente
	double moment_criterion;	//Critério de convergência para norma de momento - calculado automaticamente	
	double constraint_criterion;//Critério de convergência para norma de constraint - calculado automaticamente

	bool NaNDetector(double &value);

	void PlotDivergenceReport();//Imprime na tela informações para auxiliar o porquê da divergência

	//Nós ou SpecialConstraint associados ao máximo erro
	int super_node_disp;
	int node_disp;
	int node_rot;
	int sc_lag;
	//Nós ou SpecialConstraint associados ao máximo erro
	int super_node_force;
	int node_force;
	int node_moment;
	int sc_constraint;

	int n_conv_evaluated;	//Contador do número de vezes que o critério foi estabelecido.

};

