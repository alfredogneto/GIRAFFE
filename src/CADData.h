#pragma once
#include <stdio.h>

class Matrix;

class CADData
{
public:
	CADData() {}
	virtual ~CADData() {}
	int number;
	char file[200];																//Nome do arquivo CAD
	virtual bool Read(FILE *f) = 0;												//Leitura do arquivo de entrada
	virtual void Write(FILE *f) = 0;											//Escrita do arquivo de saida
	virtual void PreCalc() = 0;													//PreCalc
	virtual void WriteVTK_XMLRender(FILE *f, Matrix& pos, Matrix& rot, int number, int matID) = 0;		//Plots CAD Data

	//Geometric evaluation functions
	virtual void EvaluateVolume() = 0;
	virtual void EvaluateCentroid() = 0;
	virtual void EvaluateInertiaTensor() = 0;
	virtual void EvaluateRadius() = 0;

	Matrix* J_O;																//Tensor de inercia em relação a origem do CS local (CAD)
	Matrix* G;																	//Posição do baricentro em relação a origem do CS local (CAD)
	double volume;																//volume

	float radius;																//radius encompassing the whole surface (w/r to the origin)
};	