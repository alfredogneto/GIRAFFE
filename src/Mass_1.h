#pragma once
#include "Element.h"
#include "Matrix.h"

class Mass_1 :
	public Element
{
public:
	Mass_1();
	~Mass_1();
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
	void MountMassModal();				//Monta a matriz de massa para realização da analise modal
	void MountDampingModal();			//Monta a matriz de amortecimento para realização da analise modal
	void MountDamping(bool update_rayleigh);				//Monta a matriz de amortecimento
	void MountDyn();					//Montagens - Newmark
	void MountDynModal();				//Montagens para analise modal - inserção da matriz de massa e amortecimento na matriz de rigidez para posterior montagem global
	void TransformMatrix();	//Monta matriz de transformação de coordenadas
	void MountGlobal();		//Preenche a contribuição do elemento nas matrizes globais

	void SaveLagrange();	//Salva variaveis nos pontos de Gauss uteis para descrição lagrangiana atualizada
	void PreCalc();		//Pre-calculo de variaveis que e feito uma unica vez no inicio

	//Variaveis do elemento
	double m;	//massa
	Matrix c_mass;									//Matriz de massa
	Matrix c_loading;								//Vetor de esforços 
	Matrix I3;
	void Zeros();			//Zera matrizes locais do elemento
	double r;	//raio para representação grafica
	//Variaveis para calculo de steps
	int t1, t2;
	double load_multiplier, l_factor, mult;

};

