#pragma once

class TimeStepControlData
{
public:
	TimeStepControlData() {}
	~TimeStepControlData() {}
	double time_step_impact = 0.0;				//Time step a ser empregado quando da ocorrência de impacto
	int n_steps_impact = 0;						//Número de time-steps a serem empregados na resolução do impacto
	double prediction_impact = 0.0;				//Tempo previsto para início do impacto
	int steps_count_impact = 0;					//Contagem do número de steps transcorridos desde o início da interação de impacto
};	

