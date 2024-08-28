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
	Matrix* J_G;									//Tensor de in�rcia
	Matrix* G;										//Posi��o do baricentro
	
	bool Read(FILE *f);
	void Write(FILE *f);

	void PreCalc();
	void WriteVTK_XMLRender(FILE *f, int pole_node, int cs);//Plota corpo r�gido - formato XML VTK - recebe o n�mero do n� que � o p�lo e o sistema de coordenadas de refer�ncia (que � atualizado de acordo com rota��es sofridas pelo n�)
};

