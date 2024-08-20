#include "Matrix.h"

#pragma once
class TimeStepControlData
{
public:
	TimeStepControlData();
	~TimeStepControlData();
	double time_step_impact;				//Time step a ser empregado quando da ocorrência de impacto
	int n_steps_impact;						//Número de time-steps a serem empregados na resolução do impacto
	double prediction_impact;				//Tempo previsto para início do impacto
	int steps_count_impact;					//Contagem do número de steps transcorridos desde o início da interação de impacto
};	

