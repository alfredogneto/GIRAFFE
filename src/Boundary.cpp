#include "Boundary.h"

#include "Node.h"
#include "Matrix.h"
#include "BoundingVolume.h"
#include "Database.h"
//Variáveis globais
extern
Database db;


Boundary::Boundary()
{
}


Boundary::~Boundary()
{
}

void Boundary::Alloc()
{
	n_sub_bv = 0;
	sub_bv = NULL;

	node = 0;
	number = 0;
	cs = 0;
	material = 0;

	I3 = new Matrix(3, 3);
	(*I3)(0, 0) = 1.0;
	(*I3)(1, 1) = 1.0;
	(*I3)(2, 2) = 1.0;

	Q0 = new Matrix(3, 3);
	Qip = new Matrix(3, 3);
	x0ip = new Matrix(3);

	Qi = new Matrix(3, 3);
	x0i = new Matrix(3);
}

void Boundary::Free()
{
	
	delete Q0;
	delete Qip;
	delete x0ip;
	delete Qi;
	delete x0i;
}

void Boundary::UpdateVariables()
{
	//Rotation
	Matrix alpha(3);
	alpha(0, 0) = db.nodes[node - 1]->displacements[3];
	alpha(1, 0) = db.nodes[node - 1]->displacements[4];
	alpha(2, 0) = db.nodes[node - 1]->displacements[5];
	double alpha_escalar = norm(alpha);
	Matrix A = skew(alpha);
	double g = 4.0 / (4.0 + alpha_escalar * alpha_escalar);
	Matrix QdA = *I3 + g * (A + 0.5*((A)*(A)));

	if (db.nodes[node - 1]->flag_material_description)
		*Qip = (*Q0) * (*db.nodes[node - 1]->Q) * QdA;
	else
		*Qip = QdA * (*db.nodes[node - 1]->Q) * (*Q0);

	//Position
	(*x0ip)(0, 0) = db.nodes[node - 1]->copy_coordinates[0] + db.nodes[node - 1]->displacements[0];
	(*x0ip)(1, 0) = db.nodes[node - 1]->copy_coordinates[1] + db.nodes[node - 1]->displacements[1];
	(*x0ip)(2, 0) = db.nodes[node - 1]->copy_coordinates[2] + db.nodes[node - 1]->displacements[2];
}