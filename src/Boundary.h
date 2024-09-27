#pragma once
#include <stdio.h>

class Matrix;
class BoundingVolume;

class Boundary
{
public:
	Boundary();
	virtual ~Boundary();

	int node;				//Nó de referência
	int number;				//ID da boundary
	int cs;					//ID do sistema de coordenadas
	int material;			//ID do material do contorno (para atribuição de leis de interface)

	void Alloc();
	void Free();
	void UpdateVariables();

	//Funções especificas para cada tipo de particula
	virtual bool Read(FILE *f) = 0;
	virtual void Write(FILE *f) = 0;
	virtual void WriteModifyingParameters(FILE *f, int e_material, int e_node, int e_number, int e_cs) = 0;
	virtual bool Check() = 0;								//Checa inconsistências na particula para evitar erros de execução

	virtual void PreCalc() = 0;								//Pre-calculo de variaveis que e feito uma unica vez no inicio
	virtual void UpdateBoundingVolumes() = 0;				//Updates bounding volumes
	virtual void SaveLagrange() = 0;						//Salva variaveis
	virtual void WriteVTK_XMLBase(FILE *f) = 0;
	virtual void WriteVTK_XMLRender(FILE *f) = 0;

	//Variaveis - bounding volumes
	BoundingVolume* bv;												//Bounding volume
	int n_sub_bv;													//Number of sub bounding volumes
	BoundingVolume** sub_bv;										//Sub bounding volumes

	Matrix* I3;
	Matrix* Q0;														//Transformação de coordenadas
	//Variaveis para calcular estado atual (nas funcoes de contato)
	Matrix* Qip;
	Matrix* x0ip;
	//Variaveis para calcular estado anterior (nas funcoes de contato)
	Matrix* Qi;
	Matrix* x0i;

	//Explicit
	virtual void InitialEvaluations() = 0;
};

