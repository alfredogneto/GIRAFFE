#pragma once

class TimeStepControlData
{
public:
	TimeStepControlData() {}
	~TimeStepControlData() {}
	double time_step_impact = 0.0;				//Time step a ser empregado quando da ocorr�ncia de impacto
	int n_steps_impact = 0;						//Numero de time-steps a serem empregados na resolu��o do impacto
	double prediction_impact = 0.0;				//Tempo previsto para inicio do impacto
	int steps_count_impact = 0;					//Contagem do numero de steps transcorridos desde o inicio da intera��o de impacto
};	

