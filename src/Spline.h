#pragma once
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "Matrix.h"
#include "SplineElement.h"

class Spline
{
public:
	Spline();
	~Spline();

	bool Read(FILE *f);
	void Write(FILE *f);	

	void PreCalc();								//C�lculos realizados uma �nica vez
	void WriteVTK_XML_SplineMesh(FILE *f);		//Malha da spline junto ao paraview
	void WriteVTK_XML_SplineRender(FILE *f);	//Malha renderizada da spline junto ao paraview
	void CalculateSpline();						//C�lculo de pontos da spline
	void CalculateSplineTangentNormal();		//C�lculo de vetor tangente e normal a spline
	void SaveConfiguration();					//Atualiza as vari�veis internas da superf�cie, para pegarem info do pilot node para uso posterior com posi��o atualizada
	bool Check();								//Checa inconsist�ncias para evitar erros de execu��o
	
	//Vari�veis internas da spline		
	bool alloc = false;			//Controle de aloca��o da spline para destrutor
	int number;					//N�mero de refer�ncia
	double radius;				//Raio que define a superf�cie externa a spline (offset)
	int nodeset;				//Idenificador do nodeset que define os pontos de controle da spline	
	int size_nodeset;			//N�mero de n�s definidos no nodeset
	int size_sp_nodes;			//N�mero de n�s utilizados para descrever a spline apenas para visualiza��o	
	int size_sp_elements;		//N�mero de elementos (trechos) de spline

	int* nodeset_list;			//Lista de n�s do nodeset	
	double* knot;				//Knot vector calculado de acordo com a quantidade de n�s do nodeset
	double* sp0;				//Coeficientes para montar a primeira itera��o da spline (p=0)
	double* sp1;				//Coeficientes para montar a segunda itera��o da spline (p=1)
	double* sp1_dd;				//Coeficientes para montar a segunda itera��o da spline (p=1) para segunda derivada
	double** sp2;				//Coeficientes para montar a terceira itera��o da spline (p=2)
	double** sp2_d;				//Coeficientes para montar a primeria derivada da spline quadrada
	double** sp2_dd;			//Coeficientes para montar a segunda derivada da spline quadrada	
	int* sp_elements_list;		//Lista de elementos (trechos) de spline
	
	Matrix** x_sp_Ai;			//Pontos da spline
	Matrix** x_sp_d;			//Pontos da primeira derivada
	Matrix** x_sp_dd;			//Pontos da segunda derivada
	Matrix** x_sp_tangent;		//Vetor de tangentes para plotagem
	Matrix** x_sp_normal;		//Vetor de normais para plortagem

	Matrix** x_Ai;				//Pontos de controle da spline

	SplineElement** sp_element;	//Lista de elementos de spline

};

