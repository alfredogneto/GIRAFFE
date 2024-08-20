#include "TimeStepControlData.h"
#include "stdio.h"
#include "math.h"

#include"Database.h"
//Variáveis globais
extern
Database db;


////////////////////////////////////////////////////////////////////
#define _CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include <crtdbg.h>
#ifdef _DEBUG
#define DBG_NEW new ( _NORMAL_BLOCK , __FILE__ , __LINE__ )
// Replace _NORMAL_BLOCK with _CLIENT_BLOCK if you want the
// allocations to be of _CLIENT_BLOCK type
#else
#define DBG_NEW new
#endif
////////////////////////////////////////////////////////////////////

TimeStepControlData::TimeStepControlData()
{
	//Input
	n_steps_impact = 0;					//Número de time-steps a serem empregados na resolução do impacto
	//To be evaluated
	time_step_impact = 0.0;					//Estimativa de duração de impacto
	prediction_impact = 0.0;				//Tempo previsto para início do impacto
	//Update
	steps_count_impact = 0;					//Contagem do número de steps transcorridos desde o início da interação de impacto
}

TimeStepControlData::~TimeStepControlData()
{
	
}
