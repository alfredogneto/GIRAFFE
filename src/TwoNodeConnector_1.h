#pragma once
#include "Element.h"

class TwoNodeConnector_1 :
	public Element
{
public:
	TwoNodeConnector_1();
	~TwoNodeConnector_1();

	bool Check();												//Checa inconsistências no elemento para evitar erros de execução
	bool Read(FILE *f);
	void Write(FILE *f);
	void WriteResults(FILE *f);									//Escreve arquivo de resultados
	void WriteMonitor(FILE *f, bool first_record, double time);	//Escreve no monitor do elemento
	void WriteVTK_XMLBase(std::vector<float> *float_vector);
	void WriteVTK_XMLRender(FILE *f);
	void Mount();												//Monta elemento
	void MountElementLoads();									//Monta carregamentos associados ao elemento
	void MountMass();											//Monta a matriz de massa
	void MountMassModal();										//Monta a matriz de massa para realização da análise modal
	void MountDampingModal();									//Monta a matriz de amortecimento para realização da análise modal
	void MountDamping(bool update_rayleigh);					//Monta a matriz de amortecimento
	void MountDyn();											//Montagens - Newmark
	void MountDynModal();										//Montagens para análise modal - inserção da matriz de massa e amortecimento na matriz de rigidez para posterior montagem global
	void TransformMatrix();										//Monta matriz de transformação de coordenadas
	void MountGlobal();											//Preenche a contribuição do elemento nas matrizes globais
	void Zeros();												//Zera algumas matrizes utilizadas nos cálculos

	void SaveLagrange();										//Salva variáveis nos pontos de Gauss úteis para descrição lagrangiana atualizada
	void PreCalc();												//Pré-cálculo de variáveis que é feito uma única vez no início

	//Variáveis do elemento
	Matrix* K;													//Matriz de rigidez 12x12
	Matrix* C;													//Matriz de amortecimento 12x12
	Matrix* Q;													//Transformação de coordenadas 12x12

	Matrix* K_m;												//Matriz de rigidez 6x6
	Matrix* C_m;												//Matriz de amortecimento 6x6
	Matrix* Q_m;												//Transformação de coordenadas 6x6

	Matrix* f;													//Força para compor global
	Matrix* fd;													//Força - parte do amortecimento
	Matrix* Kt;													//Operador tangente para compor global
	Matrix* Ct;													//Operador tangente - parte do amortecimento
};