#pragma once
#include <stdio.h>
#include <vector>
class Matrix;

/////////////////////////////////////////////////////////////////
//IDs dos elementos do Giraffe								/////
///	1 - Beam_1												/////
///	2 - Pipe_1												/////
///	3 - Shell_1												/////
///	4 - Mass_1												/////	
///	5 - SpringDashpot_1										/////	
///	6 - RigidBody_1											/////
///	7 - Solid_1												/////
/////////////////////////////////////////////////////////////////

class Element
{
public:
	int material;	//ID do material
	int section;	//ID da se��o - ou pipe_sec (se for elemento pipe)
	int *nodes;		//N�s globais - conectividade
	int number;		//ID do elemento
	int n_nodes;	//Numero de n�s do elemento
	int cs;			//ID do sistema de coordenadas do elemento
	char* type_name;//Nome do tipo do elemento
	int VTK_type;
	int *VTK_nodes;//Indexa��o para converter numera��o do formato giraffe para o formato da celula equivalente do paraview

	int **DOFs;		//Indica para a indexa��o de cada grau de liberdade, 1 ou 0, ativo ou inativo para o elemento em quest�o

	int nDOFs;
	int nDOFs_free;
	int nDOFs_fixed;

	//Variaveis para verifica��o da fun��o Band
	int temp_band_free;
	int temp_band_fixed;

	int lowest_free_global_DOF;
	int highest_free_global_DOF;
	int lowest_fixed_global_DOF;
	int highest_fixed_global_DOF;
	bool free_marked;
	bool fixed_marked;

	double strain_energy;
	double kinetic_energy;
	double potential_gravitational_energy;

	Element();
	virtual ~Element();
	virtual bool Check() = 0;				//Checa inconsist�ncias no elemento para evitar erros de execu��o
	virtual bool Read(FILE *f) = 0;
	virtual void Write(FILE *f) = 0;
	virtual void WriteResults(FILE *f) = 0;	//Escreve arquivo de resultados
	virtual void WriteMonitor(FILE *f, bool first_record, double time) = 0;	//Escreve no monitor do elemento
	virtual void WriteVTK_XMLBase(std::vector<float> *float_vector) = 0;
	virtual void WriteVTK_XMLRender(FILE *f) = 0;
	
	virtual void Mount() = 0;									//Monta elementos
	virtual void MountMass() = 0;								//Monta a matriz de massa
	virtual void MountMassModal() = 0;							//Monta a matriz de massa para realiza��o da analise modal
	virtual void MountDampingModal() = 0;						//Monta a matriz de amortecimento para realiza��o da analise modal
	virtual void MountElementLoads() = 0;						//Monta carregamentos associados ao elemento
	virtual void MountDamping(bool update_rayleigh) = 0;		//Monta a matriz de amortecimento
	virtual void TransformMatrix() = 0;							//Monta matriz de transforma��o de coordenadas
	virtual void MountGlobal() = 0;								//Preenche a contribui��o do elemento nas matrizes globais

	virtual void SaveLagrange() = 0;	//Salva variaveis nos pontos de Gauss uteis para descri��o lagrangiana atualizada
	void Band(int* band_fixed, int* band_free);		//Calcula a banda gerada na matriz global pelo elemento

	virtual void PreCalc() = 0;		//Pre-calculo de variaveis que e feito uma unica vez no inicio
	virtual void MountDyn() = 0;	//Montagens - Newmark
	virtual void MountDynModal() = 0;	//Montagens para analise modal - inser��o da matriz de massa e amortecimento na matriz de rigidez para posterior montagem global
	virtual void Zeros() = 0;		//Limpa as matrizes internas do elemento
};

