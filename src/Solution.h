#pragma once
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "Matrix.h"
#include <vector>
#include <direct.h>
#include <omp.h>
#include <chrono>
#include <iostream>
#include <process.h>

using namespace std;
using namespace std::chrono;

class Solution
{
public:
	//Tracking solution - index and start/end times
	int solution_number;					//Solution number
	double start_time;						//Start time
	double end_time;						//End time

	Solution();
	virtual ~Solution();
	virtual bool Read(FILE *f) = 0;			//Reads input file
	virtual void Write(FILE *f) = 0;		//Writes output file
	virtual bool Solve() = 0;				//Solves solution routine

	void SetGlobalDOFs();					//Realiza a numera��o dos graus de liberdade globais
	void DOFsActive();						//Checa e aponta os graus de liberdade ativos e inativos de acordo com os elementos, special constraints e constraints
	void MountLocal();						//Montagem dos elementos e part�culas (informa��es locais)
	void MountElementLoads();				//Montagem dos carregamentos de campo em elementos (informa��es locais)
	void MountContacts();					//Montagem dos contatos (informa��es locais)
	void MountSpecialConstraints();			//Montagem das special constraints (informa��es locais)
	void MountGlobal();						//Espalhamento das informa��es locais nas matrizes/vetores globais
	void MountDisplacements();				//Inclui info de deslocamentos impostos nas matrizes/vetores globais
	void MountLoads();						//Inclui info de carregamentos impostos nas matrizes/vetores globais
	void UpdateDisps();						//Atualiza os deslocamentos nodais e multiplicadores de Lagrange (Newton-Raphson)
	void Zeros();							//Zera vari�veis iterativas para chute inicial nulo
	void ZerosVelAccel();					//Zera velocidades e acelera��es para realiza��o de an�lise est�tica (seguida da din�mica)
	void Clear();							//Limpa as matrizes globais
	void SaveConfiguration();				//Salva configura��o convergida
	void RestoreConfiguration();			//Restaura a �ltima configura��o que convergiu
	void SetGlobalSize();					//Seta o tamanho da matriz de rigidez global
	void PinballCheck();					//Varre contatos e verifica pinball
	void BeginStepCheck();					//checagem inicial de in�cio de step para contatos
	void MountMass();						//Monta matriz de massa e esforcos inerciais
	void MountDamping(bool update_rayleigh);//Monta matriz de amortecimento
	void MountDyn();						//Composicao de rigidez e esforcos para analise dinamica (informa��es locais)
	void ComputeInitialConditions(bool zero_ICs);		//Computa as condi��es inciais nodais impostas
	bool HaveErrors();						//Verifica erros que impedem avan�a da an�lise, mesmo em caso de converg�ncia
	void MountSparse();						//Monta matrizes esparsas
};

