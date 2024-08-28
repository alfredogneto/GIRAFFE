#pragma once

class TimeStepControlData
{
public:
	TimeStepControlData() {}
	~TimeStepControlData() {}
	double time_step_impact = 0.0;				//Time step a ser empregado quando da ocorr�ncia de impacto
	int n_steps_impact = 0;						//N�mero de time-steps a serem empregados na resolu��o do impacto
	double prediction_impact = 0.0;				//Tempo previsto para in�cio do impacto
	int steps_count_impact = 0;					//Contagem do n�mero de steps transcorridos desde o in�cio da intera��o de impacto
};	

