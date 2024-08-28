#include "RigidBodyData.h"


#include "Matrix.h"
#include "CoordinateSystem.h"
#include "Node.h"
#include "CADData.h"
#include"Database.h"
//Variáveis globais
extern
Database db;
#define PI 3.1415926535897932384626433832795

RigidBodyData::RigidBodyData()
{
	number = 0;
	mass = 0.0;
	J_G = new Matrix(3,3);
	G = new Matrix(3);
	CADData_ID = 0;
	CAD_entered = false;
}

RigidBodyData::~RigidBodyData()
{
	delete J_G;
	delete G;
}

bool RigidBodyData::Read(FILE *f)
{
	char s[1000];
	
	fscanf(f, "%s", s);
	//Verifica a palavra chave "RB"
	if (!strcmp(s, "RBData"))
	{
		fscanf(f, "%s", s);
		number = atoi(s);
	}
	else
		return false;
	//Massa
	fscanf(f, "%s", s);
	if (!strcmp(s, "Mass"))
	{
		fscanf(f, "%s", s);
		mass = atof(s);
	}
	else
		return false;
	//Tensor de inércia
	fscanf(f, "%s", s);
	if (!strcmp(s, "J11"))
	{
		fscanf(f, "%s", s);
		(*J_G)(0, 0) = atof(s);
	}
	else
		return false;
	fscanf(f, "%s", s);
	if (!strcmp(s, "J22"))
	{
		fscanf(f, "%s", s);
		(*J_G)(1, 1) = atof(s);
	}
	else
		return false;
	fscanf(f, "%s", s);
	if (!strcmp(s, "J33"))
	{
		fscanf(f, "%s", s);
		(*J_G)(2, 2) = atof(s);
	}
	else
		return false;
	fscanf(f, "%s", s);
	if (!strcmp(s, "J12"))
	{
		fscanf(f, "%s", s);
		(*J_G)(0, 1) = -atof(s);
		(*J_G)(1, 0) = -atof(s);
	}
	else
		return false;
	fscanf(f, "%s", s);
	if (!strcmp(s, "J13"))
	{
		fscanf(f, "%s", s);
		(*J_G)(0, 2) = -atof(s);
		(*J_G)(2, 0) = -atof(s);
	}
	else
		return false;
	fscanf(f, "%s", s);
	if (!strcmp(s, "J23"))
	{
		fscanf(f, "%s", s);
		(*J_G)(1, 2) = -atof(s);
		(*J_G)(2, 1) = -atof(s);
	}
	else
		return false;
	//Baricentro
	fscanf(f, "%s", s);
	if (!strcmp(s, "Barycenter"))
	{
		fscanf(f, "%s", s);
		(*G)(0, 0) = atof(s);
		fscanf(f, "%s", s);
		(*G)(1, 0) = atof(s);
		fscanf(f, "%s", s);
		(*G)(2, 0) = atof(s);
	}
	else
		return false;
	
	//Salva a posição (stream)
	fpos_t pos;
	fgetpos(f, &pos);
	fscanf(f, "%s", s);
	//CADData
	if (!strcmp(s, "CADData"))
	{
		CAD_entered = true;
		fscanf(f, "%s", s);
		CADData_ID = atoi(s);
	}
	else
		fsetpos(f, &pos);
	
	return true;
}

void RigidBodyData::Write(FILE *f)
{
	fprintf(f, "RBData\t%d\n",number);
	fprintf(f, "Mass\t%.6f\n", mass);
	fprintf(f, "J11\t%.6f\tJ22\t%.6f\tJ33\t%.6f\tJ12\t%.6f\tJ13\t%.6f\tJ23\t%.6f\n", (*J_G)(0, 0), (*J_G)(1, 1), (*J_G)(2, 2), (*J_G)(0, 1), (*J_G)(0, 2), (*J_G)(1, 2));
	fprintf(f, "Barycenter\t%.6f\t%.6f\t%.6f\n", (*G)(0, 0), (*G)(1, 0), (*G)(2, 0));
	if (CAD_entered == true)
		fprintf(f, "CADData\t%d\n", CADData_ID);
}

void RigidBodyData::PreCalc()
{
	
}

//Plota corpo rígido - formato XML VTK - recebe o número do nó que é o pólo e o sistema de coordenadas de referência (que é atualizado de acordo com rotações sofridas pelo nó)
void RigidBodyData::WriteVTK_XMLRender(FILE *f, int pole_node, int cs)
{
	if (CAD_entered == true)
	{
		//Transformação de coordenadas - local para global
		Matrix Q(3, 3);
		Q(0, 0) = (*db.CS[cs - 1]->E1)(0, 0);
		Q(1, 0) = (*db.CS[cs - 1]->E1)(1, 0);
		Q(2, 0) = (*db.CS[cs - 1]->E1)(2, 0);
		Q(0, 1) = (*db.CS[cs - 1]->E2)(0, 0);
		Q(1, 1) = (*db.CS[cs - 1]->E2)(1, 0);
		Q(2, 1) = (*db.CS[cs - 1]->E2)(2, 0);
		Q(0, 2) = (*db.CS[cs - 1]->E3)(0, 0);
		Q(1, 2) = (*db.CS[cs - 1]->E3)(1, 0);
		Q(2, 2) = (*db.CS[cs - 1]->E3)(2, 0);
		Matrix Q_alpha(3, 3);
		Matrix alpha(3);
		Matrix I3(3, 3);
		I3(0, 0) = 1.0;
		I3(1, 1) = 1.0;
		I3(2, 2) = 1.0;
		//Cálculo do Q_alpha - rotação do slave node
		alpha(0, 0) = db.nodes[pole_node - 1]->copy_coordinates[3];
		alpha(1, 0) = db.nodes[pole_node - 1]->copy_coordinates[4];
		alpha(2, 0) = db.nodes[pole_node - 1]->copy_coordinates[5];
		double alpha_escalar = norm(alpha);
		double g = 4.0 / (4.0 + alpha_escalar * alpha_escalar);
		Q_alpha = I3 + g * (skew(alpha) + 0.5*skew(alpha)*skew(alpha));
		Q = Q_alpha * Q;	//Aplicar o operador Q para escrever os pontos no sistema de coordenadas global - rotação de corpo rígido - entrar com vetor de pontos no sistema local
		Matrix xO(3);
		xO(0, 0) = db.nodes[pole_node - 1]->copy_coordinates[0];
		xO(1, 0) = db.nodes[pole_node - 1]->copy_coordinates[1];
		xO(2, 0) = db.nodes[pole_node - 1]->copy_coordinates[2];

		db.cad_data[CADData_ID - 1]->WriteVTK_XMLRender(f, xO, Q, number);
	}
}