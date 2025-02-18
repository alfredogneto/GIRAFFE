#pragma once
#include <stdio.h>

using namespace std;

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

	void SetGlobalDOFs();					//Realiza a numeração dos graus de liberdade globais
	void DOFsActive();						//Checa e aponta os graus de liberdade ativos e inativos de acordo com os elementos, special constraints e constraints
	void MountLocal();						//Montagem dos elementos e particulas (informações locais)
	void MountElementLoads();				//Montagem dos carregamentos de campo em elementos (informações locais)
	void MountContacts();					//Montagem dos contatos (informações locais)
	void MountSpecialConstraints();			//Montagem das special constraints (informações locais)
	void MountGlobal();						//Espalhamento das informações locais nas matrizes/vetores globais
	void MountDisplacements();				//Inclui info de deslocamentos impostos nas matrizes/vetores globais
	void MountLoads();						//Inclui info de carregamentos impostos nas matrizes/vetores globais
	void UpdateDisps();						//Atualiza os deslocamentos nodais e multiplicadores de Lagrange (Newton-Raphson)
	void Zeros();							//Zera variaveis iterativas para chute inicial nulo
	void ZerosVelAccel();					//Zera velocidades e acelerações para realização de analise estatica (seguida da dinamica)
	void Clear();							//Limpa as matrizes globais
	void SaveConfiguration();				//Salva configuração convergida
	void RestoreConfiguration();			//Restaura a ultima configuração que convergiu
	void SetGlobalSize();					//Seta o tamanho da matriz de rigidez global
	void PinballCheck();					//Varre contatos e verifica pinball
	void BeginStepCheck();					//checagem inicial de inicio de step para contatos
	void MountMass();						//Monta matriz de massa e esforcos inerciais
	void MountDamping(bool update_rayleigh);//Monta matriz de amortecimento
	void MountDyn();						//Composicao de rigidez e esforcos para analise dinamica (informações locais)
	void ComputeInitialConditions(bool zero_ICs);		//Computa as condições inciais nodais impostas
	bool HaveErrors();						//Verifica erros que impedem avança da analise, mesmo em caso de convergência
	void MountSparse();						//Monta matrizes esparsas
};

