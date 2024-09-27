#pragma once
#include <stdio.h>

#include "Box.h"

class Matrix;

class SplineElement
{
public:
	SplineElement();
	~SplineElement();

	void PreCalc();

	bool ReadCommon(FILE *f);

	//Funções
	void FillNodes();														//Atualiza as variaveis internas da superficie, para pegarem info do pilot node para uso posterior com posição atualizada
	void SaveConfiguration();												//Salva vetores de configuração convergida
	void CenterPoint(Matrix* center, double* radii);						//Retorna coordenadas globais do ponto central da superficie a ser utilizado para calculos grosseiros de sua localização (pinball)
	bool Check();															//Checa inconsistências para evitar erros de execução
	void SetMinMaxRange();
	void SplinePoint(double & zeta, Matrix& point);					//Obtem ponto da spline
	void UpdateBox();

	int *nodes;		//Nós globais - conectividade
	int n_nodes;	//Numero de nós da superficie
	int **DOFs;		//Indica para a indexação de cada grau de liberdade, 1 ou 0, ativo ou inativo para o elemento em questão
	int nDOFs;		//Numero de GL
	//int VTK_type;	//Tipo de celula para VTK
	int *VTK_nodes;	//Indexação para converter numeração do formato giraffe para o formato da celula equivalente do paraview
	int **GLs;		//Ponteiro para os GL globais utilizados na superficie

	Box box;

	//Bool indicators of coordinates
	bool entered_u1;
	double u1_min, u1_max, u1_range;
	//Divisions on degeneration
	int div1, div2;

	double* knot_element;

	double* radius;

	Matrix* x_Ai;
	Matrix* x_Bi;
	Matrix* x_Ci;

	Matrix* x_Ap;
	Matrix* x_Bp;
	Matrix* x_Cp;

	Matrix* d;
	Matrix* dui;
	Matrix* ddui;

	Matrix* I3;																//Identidade de ordem 3

};

