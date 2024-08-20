#pragma once
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "Matrix.h"
#define PI 3.1415926535897932384626433832795
#include <vector>
class RigidBodyData
{
public:
	RigidBodyData();
	~RigidBodyData();

	int number;										//ID
	int CADData_ID;									//ID do CAD
	bool CAD_entered;
	double mass;									//massa
	Matrix* J_G;									//Tensor de inércia
	Matrix* G;										//Posição do baricentro
	
	bool Read(FILE *f);
	void Write(FILE *f);

	void PreCalc();
	void WriteVTK_XMLRender(FILE *f, int pole_node, int cs);//Plota corpo rígido - formato XML VTK - recebe o número do nó que é o pólo e o sistema de coordenadas de referência (que é atualizado de acordo com rotações sofridas pelo nó)
};

