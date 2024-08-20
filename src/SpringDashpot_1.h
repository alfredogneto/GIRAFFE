#pragma once
#include "Element.h"
class SpringDashpot_1 :
	public Element
{
public:
	SpringDashpot_1();
	~SpringDashpot_1();
	bool Check();				//Checa inconsistências no elemento para evitar erros de execução
	bool Read(FILE *f);
	void Write(FILE *f);
	void WriteResults(FILE *f);//Escreve arquivo de resultados
	void WriteMonitor(FILE *f, bool first_record, double time);	//Escreve no monitor do elemento
	void WriteVTK_XMLBase(std::vector<float> *float_vector);
	void WriteVTK_XMLRender(FILE *f);
	void Mount();			//Monta elementos
	void MountElementLoads();			//Monta carregamentos associados ao elemento
	void MountMass();					//Monta a matriz de massa
	void MountMassModal();				//Monta a matriz de massa para realização da análise modal
	void MountDamping(bool update_rayleigh);				//Monta a matriz de amortecimento
	void MountDampingModal();			//Monta a matriz de amortecimento para realização da análise modal
	void MountDyn();					//Montagens - Newmark
	void MountDynModal();				//Montagens para análise modal - inserção da matriz de massa e amortecimento na matriz de rigidez para posterior montagem global
	void TransformMatrix();	//Monta matriz de transformação de coordenadas
	void MountGlobal();		//Preenche a contribuição do elemento nas matrizes globais

	void SaveLagrange();	//Salva variáveis nos pontos de Gauss úteis para descrição lagrangiana atualizada
	void PreCalc();		//Pré-cálculo de variáveis que é feito uma única vez no início
	void Zeros();			//Zera matrizes locais do elemento

	//Variáveis do elemento
	double k;	//Stiffness
	double c;	//Damping
	double initial_distance;
	double gn;		//gap normal
	Matrix z1, z2;	//Posições dos nós
	Matrix xd1, xd2;//Velocidades dos nós
	Matrix c_stiffness;									//Matriz de rigidez
	Matrix c_damping;									//Matriz de amortecimento
	Matrix c_damping_modal;									//Matriz de amortecimento
	Matrix c_stiffness_force;							//Vetor de esforços elásticos
	Matrix c_damping_force;								//Vetor de esforços de amortecimento
	Matrix I3;
	Matrix z1z2;	//distância atual entre nós
	Matrix n;		//direção normal da mola
	
	Matrix non;
	double f;
	Matrix C1;
	Matrix C2;
	Matrix C3;

	//Pós-processamento
	double elastic_force;
	double damping_force;

	Matrix last_n;	//última direção normal convergida da mola
	bool first_evaluation;
};

