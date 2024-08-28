#pragma once
#include "Element.h"

class TwoNodeConnector_1 :
	public Element
{
public:
	TwoNodeConnector_1();
	~TwoNodeConnector_1();

	bool Check();												//Checa inconsist�ncias no elemento para evitar erros de execu��o
	bool Read(FILE *f);
	void Write(FILE *f);
	void WriteResults(FILE *f);									//Escreve arquivo de resultados
	void WriteMonitor(FILE *f, bool first_record, double time);	//Escreve no monitor do elemento
	void WriteVTK_XMLBase(std::vector<float> *float_vector);
	void WriteVTK_XMLRender(FILE *f);
	void Mount();												//Monta elemento
	void MountElementLoads();									//Monta carregamentos associados ao elemento
	void MountMass();											//Monta a matriz de massa
	void MountMassModal();										//Monta a matriz de massa para realiza��o da an�lise modal
	void MountDampingModal();									//Monta a matriz de amortecimento para realiza��o da an�lise modal
	void MountDamping(bool update_rayleigh);					//Monta a matriz de amortecimento
	void MountDyn();											//Montagens - Newmark
	void MountDynModal();										//Montagens para an�lise modal - inser��o da matriz de massa e amortecimento na matriz de rigidez para posterior montagem global
	void TransformMatrix();										//Monta matriz de transforma��o de coordenadas
	void MountGlobal();											//Preenche a contribui��o do elemento nas matrizes globais
	void Zeros();												//Zera algumas matrizes utilizadas nos c�lculos

	void SaveLagrange();										//Salva vari�veis nos pontos de Gauss �teis para descri��o lagrangiana atualizada
	void PreCalc();												//Pr�-c�lculo de vari�veis que � feito uma �nica vez no in�cio

	//Vari�veis do elemento
	Matrix* K;													//Matriz de rigidez 12x12
	Matrix* C;													//Matriz de amortecimento 12x12
	Matrix* Q;													//Transforma��o de coordenadas 12x12

	Matrix* K_m;												//Matriz de rigidez 6x6
	Matrix* C_m;												//Matriz de amortecimento 6x6
	Matrix* Q_m;												//Transforma��o de coordenadas 6x6

	Matrix* f;													//For�a para compor global
	Matrix* fd;													//For�a - parte do amortecimento
	Matrix* Kt;													//Operador tangente para compor global
	Matrix* Ct;													//Operador tangente - parte do amortecimento
};