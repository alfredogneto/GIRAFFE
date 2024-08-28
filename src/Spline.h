#pragma once
#include <stdio.h>

class Matrix;
class SplineElement;

class Spline
{
public:
	Spline();
	~Spline();

	bool Read(FILE *f);
	void Write(FILE *f);	

	void PreCalc();								//Calculos realizados uma unica vez
	void WriteVTK_XML_SplineMesh(FILE *f);		//Malha da spline junto ao paraview
	void WriteVTK_XML_SplineRender(FILE *f);	//Malha renderizada da spline junto ao paraview
	void CalculateSpline();						//Calculo de pontos da spline
	void CalculateSplineTangentNormal();		//Calculo de vetor tangente e normal a spline
	void SaveConfiguration();					//Atualiza as variaveis internas da superficie, para pegarem info do pilot node para uso posterior com posição atualizada
	bool Check();								//Checa inconsistências para evitar erros de execução
	
	//Variaveis internas da spline		
	bool alloc = false;			//Controle de alocação da spline para destrutor
	int number;					//Numero de referência
	double radius;				//Raio que define a superficie externa a spline (offset)
	int nodeset;				//Idenificador do nodeset que define os pontos de controle da spline	
	int size_nodeset;			//Numero de nós definidos no nodeset
	int size_sp_nodes;			//Numero de nós utilizados para descrever a spline apenas para visualização	
	int size_sp_elements;		//Numero de elementos (trechos) de spline

	int* nodeset_list;			//Lista de nós do nodeset	
	double* knot;				//Knot vector calculado de acordo com a quantidade de nós do nodeset
	double* sp0;				//Coeficientes para montar a primeira iteração da spline (p=0)
	double* sp1;				//Coeficientes para montar a segunda iteração da spline (p=1)
	double* sp1_dd;				//Coeficientes para montar a segunda iteração da spline (p=1) para segunda derivada
	double** sp2;				//Coeficientes para montar a terceira iteração da spline (p=2)
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

