#pragma once
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "Matrix.h"
#include "Box.h"

class SplineElement
{
public:
	SplineElement();
	~SplineElement();

	void PreCalc();

	bool ReadCommon(FILE *f);

	//Fun��es
	void FillNodes();														//Atualiza as vari�veis internas da superf�cie, para pegarem info do pilot node para uso posterior com posi��o atualizada
	void SaveConfiguration();												//Salva vetores de configura��o convergida
	void CenterPoint(Matrix* center, double* radii);						//Retorna coordenadas globais do ponto central da superf�cie a ser utilizado para c�lculos grosseiros de sua localiza��o (pinball)
	bool Check();															//Checa inconsist�ncias para evitar erros de execu��o
	void SetMinMaxRange();
	void SplinePoint(double & zeta, Matrix& point);					//Obtem ponto da spline
	void UpdateBox();

	int *nodes;		//N�s globais - conectividade
	int n_nodes;	//N�mero de n�s da superf�cie
	int **DOFs;		//Indica para a indexa��o de cada grau de liberdade, 1 ou 0, ativo ou inativo para o elemento em quest�o
	int nDOFs;		//N�mero de GL
	//int VTK_type;	//Tipo de c�lula para VTK
	int *VTK_nodes;	//Indexa��o para converter numera��o do formato giraffe para o formato da c�lula equivalente do paraview
	int **GLs;		//Ponteiro para os GL globais utilizados na superf�cie

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

