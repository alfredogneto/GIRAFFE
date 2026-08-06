#pragma once

class TimeStepControlData
{
public:
	TimeStepControlData() {}
	~TimeStepControlData() {}
	double time_step_impact = 0.0;				//Time step the be empregado when of the ocorrência of impacto
	int n_steps_impact = 0;						//Number of time-steps the to be empregados in the resolução of the impacto
	double prediction_impact = 0.0;				//Time previsto for start of the impacto
	int steps_count_impact = 0;					//Contagem of the number of steps transcorridos desde the start of the interação of impacto
};	

