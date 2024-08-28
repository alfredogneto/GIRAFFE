#pragma once
#include <stdio.h>

class Matrix;

class RigidBodyData
{
public:
	RigidBodyData();
	~RigidBodyData();

	int number;										//ID
	int CADData_ID;									//ID do CAD
	bool CAD_entered;
	double mass;									//massa
	Matrix* J_G;									//Tensor de inercia
	Matrix* G;										//Posi��o do baricentro
	
	bool Read(FILE *f);
	void Write(FILE *f);

	void PreCalc();
	void WriteVTK_XMLRender(FILE *f, int pole_node, int cs);//Plota corpo rigido - formato XML VTK - recebe o numero do n� que e o p�lo e o sistema de coordenadas de refer�ncia (que e atualizado de acordo com rota��es sofridas pelo n�)
};

