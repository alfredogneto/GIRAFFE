#include "Matrix.h"

#pragma once
class TimeStepControlData
{
public:
	TimeStepControlData();
	~TimeStepControlData();
	double time_step_impact;				//Time step a ser empregado quando da ocorr�ncia de impacto
	int n_steps_impact;						//N�mero de time-steps a serem empregados na resolu��o do impacto
	double prediction_impact;				//Tempo previsto para in�cio do impacto
	int steps_count_impact;					//Contagem do n�mero de steps transcorridos desde o in�cio da intera��o de impacto
};	

