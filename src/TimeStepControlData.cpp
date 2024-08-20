#include "TimeStepControlData.h"
#include "stdio.h"
#include "math.h"

#include"Database.h"
//Vari�veis globais
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
	n_steps_impact = 0;					//N�mero de time-steps a serem empregados na resolu��o do impacto
	//To be evaluated
	time_step_impact = 0.0;					//Estimativa de dura��o de impacto
	prediction_impact = 0.0;				//Tempo previsto para in�cio do impacto
	//Update
	steps_count_impact = 0;					//Contagem do n�mero de steps transcorridos desde o in�cio da intera��o de impacto
}

TimeStepControlData::~TimeStepControlData()
{
	
}
