#pragma once
#include <stdio.h>

class Matrix;
class BoundingVolume;

class Boundary
{
public:
	Boundary();
	virtual ~Boundary();

	int node;				//N� de refer�ncia
	int number;				//ID da boundary
	int cs;					//ID do sistema de coordenadas
	int material;			//ID do material do contorno (para atribui��o de leis de interface)

	void Alloc();
	void Free();
	void UpdateVariables();

	//Fun��es espec�ficas para cada tipo de part�cula
	virtual bool Read(FILE *f) = 0;
	virtual void Write(FILE *f) = 0;
	virtual void WriteModifyingParameters(FILE *f, int e_material, int e_node, int e_number, int e_cs) = 0;
	virtual bool Check() = 0;								//Checa inconsist�ncias na part�cula para evitar erros de execu��o

	virtual void PreCalc() = 0;								//Pr�-c�lculo de vari�veis que � feito uma �nica vez no in�cio
	virtual void UpdateBoundingVolumes() = 0;				//Updates bounding volumes
	virtual void SaveLagrange() = 0;						//Salva vari�veis
	virtual void WriteVTK_XMLBase(FILE *f) = 0;
	virtual void WriteVTK_XMLRender(FILE *f) = 0;

	//Vari�veis - bounding volumes
	BoundingVolume* bv;												//Bounding volume
	int n_sub_bv;													//Number of sub bounding volumes
	BoundingVolume** sub_bv;										//Sub bounding volumes

	Matrix* I3;
	Matrix* Q0;														//Transforma��o de coordenadas
	//Vari�veis para calcular estado atual (nas funcoes de contato)
	Matrix* Qip;
	Matrix* x0ip;
	//Vari�veis para calcular estado anterior (nas funcoes de contato)
	Matrix* Qi;
	Matrix* x0i;

	//Explicit
	virtual void InitialEvaluations() = 0;
};

