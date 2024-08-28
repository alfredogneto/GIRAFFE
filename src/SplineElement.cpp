#include "SplineElement.h"

#include "Matrix.h"
#include "Node.h"

#include "Database.h"
extern
Database db;

SplineElement::SplineElement()
{
	nDOFs = 9;
	n_nodes = 3;
	nodes = new int[n_nodes];
	VTK_nodes = new int[n_nodes];	
	DOFs = new int *[n_nodes];
	for (int i = 0; i < n_nodes; i++)
		DOFs[i] = new int[db.number_GLs_node];
	VTK_nodes[0] = 0;
	VTK_nodes[1] = 1;
	VTK_nodes[2] = 2;

	knot_element = new double[6];

	radius = new double;

	x_Ai = new Matrix(3);
	x_Bi = new Matrix(3);
	x_Ci = new Matrix(3);

	x_Ap = new Matrix(3);
	x_Bp = new Matrix(3);
	x_Cp = new Matrix(3);

	d = new Matrix(9);
	dui = new Matrix(9);
	ddui = new Matrix(9);

	//Rotina para ativar os GLS de cada nó
	for (int i = 0; i < n_nodes; i++)
	{
		for (int j = 0; j < db.number_GLs_node; j++)
		{
			DOFs[i][j] = 0;
		}
		//Translação
		DOFs[i][0] = 1;
		DOFs[i][1] = 1;
		DOFs[i][2] = 1;
		//Rotação
		DOFs[i][3] = 0;
		DOFs[i][4] = 0;
		DOFs[i][5] = 0;
	}
	GLs = new int*[nDOFs];
	for (int i = 0; i < nDOFs; i++)
		GLs[i] = NULL;

	I3 = new Matrix(3, 3);
	(*I3)(0, 0) = 1.0;
	(*I3)(1, 1) = 1.0;
	(*I3)(2, 2) = 1.0;

	div1 = 1;
	div2 = 1;
}

SplineElement::~SplineElement()
{
	delete[] nodes;
	delete[] VTK_nodes;
	if (DOFs != NULL)
	{
		for (int i = 0; i < 3; i++)
			delete[] DOFs[i];
		delete[] DOFs;
	}
	delete[]GLs;

	delete[] knot_element;

	delete[] radius;

	delete x_Ai;
	delete x_Bi;
	delete x_Ci;

	delete x_Ap;
	delete x_Bp;
	delete x_Cp;

	delete d;
	delete dui;
	delete ddui;

	delete I3;
}

void SplineElement::SetMinMaxRange()
{
	u1_min = knot_element[2];
	u1_max = knot_element[3];

	u1_range = knot_element[3] - knot_element[2];
}

void SplineElement::SplinePoint(double & zeta, Matrix & point)
{
	double sp0[5];
	double sp1[4];
	double sp2[3];

	if (zeta != 1) {
		for (int i = 0; i < 5; i++) {
			if (knot_element[i] <= zeta && zeta < knot_element[i + 1])
				sp0[i] = 1;
			else
				sp0[i] = 0;
		}
	}

	double a, b;

	int p = 1;
	for (int i = 0; i < 4; i++) {
		if (knot_element[i + p] == knot_element[i])
			a = 0;
		else
			a = sp0[i] * (zeta - knot_element[i]) / (knot_element[i + p] - knot_element[i]);
		if (knot_element[i + p + 1] == knot_element[i + 1])
			b = 0;
		else
			b = sp0[i + 1] * (knot_element[i + p + 1] - zeta) / (knot_element[i + p + 1] - knot_element[i + 1]);
		sp1[i] = a + b;
	}

	p = 2;
	for (int i = 0; i < 3; i++) {
		if (knot_element[i + p] == knot_element[i])
			a = 0;
		else
			a = sp1[i] * (zeta - knot_element[i]) / (knot_element[i + p] - knot_element[i]);
		if (knot_element[i + p + 1] == knot_element[i + 1])
			b = 0;
		else
			b = sp1[i + 1] * (knot_element[i + p + 1] - zeta) / (knot_element[i + p + 1] - knot_element[i + 1]);
		sp2[i] = a + b;
	}

	if (zeta == 1)
	{
		sp2[0] = 0;
		sp2[1] = 0;
		sp2[2] = 1;
	}

	point(0, 0) = sp2[0] * (*x_Ap)(0, 0) + sp2[1] * (*x_Bp)(0, 0) + sp2[2] * (*x_Cp)(0, 0);
	point(1, 0) = sp2[0] * (*x_Ap)(1, 0) + sp2[1] * (*x_Bp)(1, 0) + sp2[2] * (*x_Cp)(1, 0);
	point(2, 0) = sp2[0] * (*x_Ap)(2, 0) + sp2[1] * (*x_Bp)(2, 0) + sp2[2] * (*x_Cp)(2, 0);
}

void SplineElement::UpdateBox()
{
	//Pontos da spline
	Matrix x_spA(3);
	Matrix x_spB(3);
	Matrix x_spC(3);


	SplinePoint(knot_element[3], x_spA);
	double zeta = knot_element[3] / 2 + knot_element[4] / 2;
	SplinePoint(zeta, x_spB);
	SplinePoint(knot_element[4], x_spC);

	double x_min = std::min(x_spA(0, 0), std::min(x_spB(0, 0), x_spC(0, 0))) - *radius * 4;
	double x_max = std::max(x_spA(0, 0), std::max(x_spB(0, 0), x_spC(0, 0))) + *radius * 4;
	double y_min = std::min(x_spA(1, 0), std::min(x_spB(1, 0), x_spC(1, 0))) - *radius * 4;
	double y_max = std::max(x_spA(1, 0), std::max(x_spB(1, 0), x_spC(1, 0))) + *radius * 4;
	double z_min = std::min(x_spA(2, 0), std::min(x_spB(2, 0), x_spC(2, 0))) - *radius * 4;
	double z_max = std::max(x_spA(2, 0), std::max(x_spB(2, 0), x_spC(2, 0))) + *radius * 4;

	//Bounding box
	double x[8], y[8], z[8];
	Matrix vec(3);
	//Vertex A-AA
	vec(0, 0) = x_min;
	vec(1, 0) = y_min;
	vec(2, 0) = z_min;
	x[0] = vec(0, 0);
	y[0] = vec(1, 0);
	z[0] = vec(2, 0);
	//Vertex B-AA
	vec(0, 0) = x_min;
	vec(1, 0) = y_max;
	vec(2, 0) = z_min;
	x[1] = vec(0, 0);
	y[1] = vec(1, 0);
	z[1] = vec(2, 0);
	//Vertex C-AA
	vec(0, 0) = x_min;
	vec(1, 0) = y_min;
	vec(2, 0) = z_max;
	x[2] = vec(0, 0);
	y[2] = vec(1, 0);
	z[2] = vec(2, 0);
	//Vertex D-AA
	vec(0, 0) = x_min;
	vec(1, 0) = y_max;
	vec(2, 0) = z_max;
	x[3] = vec(0, 0);
	y[3] = vec(1, 0);
	z[3] = vec(2, 0);

	//Vertex A-BA
	vec(0, 0) = x_max;
	vec(1, 0) = y_min;
	vec(2, 0) = z_min;
	x[4] = vec(0, 0);
	y[4] = vec(1, 0);
	z[4] = vec(2, 0);
	//Vertex B-BA
	vec(0, 0) = x_max;
	vec(1, 0) = y_max;
	vec(2, 0) = z_min;
	x[5] = vec(0, 0);
	y[5] = vec(1, 0);
	z[5] = vec(2, 0);
	//Vertex C-BA
	vec(0, 0) = x_max;
	vec(1, 0) = y_min;
	vec(2, 0) = z_max;
	x[6] = vec(0, 0);
	y[6] = vec(1, 0);
	z[6] = vec(2, 0);
	//Vertex D-BA
	vec(0, 0) = x_max;
	vec(1, 0) = y_max;
	vec(2, 0) = z_max;
	x[7] = vec(0, 0);
	y[7] = vec(1, 0);
	z[7] = vec(2, 0);

	//Setando os vértices
	box.SetVertices(x, y, z);
}

void SplineElement::PreCalc()
{
	//Apontando para posição que indica valor dos GLs globais
	for (int i = 0; i < 3; i++)
	{
		GLs[i] = &db.nodes[nodes[0] - 1]->GLs[i];
		GLs[i + 3] = &db.nodes[nodes[1] - 1]->GLs[i];
		GLs[i + 6] = &db.nodes[nodes[2] - 1]->GLs[i];
	}

	////Transformação de coordenadas
	//Matrix e1g(3);
	//e1g(0, 0) = 1.0;
	//Matrix e2g(3);
	//e2g(1, 0) = 1.0;
	//Matrix e3g(3);
	//e3g(2, 0) = 1.0;
	//Matrix e1l = *db.CS[csA - 1]->E1;
	//Matrix e2l = *db.CS[csA - 1]->E2;
	//Matrix e3l = *db.CS[csA - 1]->E3;

	////Salva a matriz de transformação de coordenadas (para orientar o plano da ST de acordo com a orientação de referência da superfície)
	//(*Q0A)(0, 0) = dot(e1g, e1l);
	//(*Q0A)(0, 1) = dot(e1g, e2l);
	//(*Q0A)(0, 2) = dot(e1g, e3l);

	//(*Q0A)(1, 0) = dot(e2g, e1l);
	//(*Q0A)(1, 1) = dot(e2g, e2l);
	//(*Q0A)(1, 2) = dot(e2g, e3l);

	//(*Q0A)(2, 0) = dot(e3g, e1l);
	//(*Q0A)(2, 1) = dot(e3g, e2l);
	//(*Q0A)(2, 2) = dot(e3g, e3l);

	//e1l = *db.CS[csB - 1]->E1;
	//e2l = *db.CS[csB - 1]->E2;
	//e3l = *db.CS[csB - 1]->E3;

	////Salva a matriz de transformação de coordenadas (para orientar o plano da ST de acordo com a orientação de referência da superfície)
	//(*Q0B)(0, 0) = dot(e1g, e1l);
	//(*Q0B)(0, 1) = dot(e1g, e2l);
	//(*Q0B)(0, 2) = dot(e1g, e3l);

	//(*Q0B)(1, 0) = dot(e2g, e1l);
	//(*Q0B)(1, 1) = dot(e2g, e2l);
	//(*Q0B)(1, 2) = dot(e2g, e3l);

	//(*Q0B)(2, 0) = dot(e3g, e1l);
	//(*Q0B)(2, 1) = dot(e3g, e2l);
	//(*Q0B)(2, 2) = dot(e3g, e3l);

	SaveConfiguration();

}

void SplineElement::FillNodes()
{
	Matrix uA(3);
	Matrix uB(3);
	Matrix uC(3);

	for (int i = 0; i < 3; i++)
	{
		(*d)(i + 0, 0) = db.nodes[nodes[0] - 1]->displacements[i];
		(*d)(i + 3, 0) = db.nodes[nodes[1] - 1]->displacements[i];
		(*d)(i + 6, 0) = db.nodes[nodes[2] - 1]->displacements[i];

		uA(i, 0) = db.nodes[nodes[0] - 1]->displacements[i];
		uB(i, 0) = db.nodes[nodes[1] - 1]->displacements[i];
		uC(i, 0) = db.nodes[nodes[2] - 1]->displacements[i];
	}
	//x_Ap
	*x_Ap = *x_Ai + uA;
	//x_Bp
	*x_Bp = *x_Bi + uB;
	//x_Cp
	*x_Cp = *x_Ci + uC;
}

void SplineElement::CenterPoint(Matrix* center, double* radii)
{
	double sp0[5];
	double sp1[4];
	double sp2[3];

	double zeta = (knot_element[2] + knot_element[3]) / 2;

	if (zeta != 1)
	{
		for (int i = 0; i < 5; i++) {
			if (knot_element[i] <= zeta && zeta < knot_element[i + 1])
				sp0[i] = 1;
			else
				sp0[i] = 0;
		}
	}

	double a, b;

	int p = 1;
	for (int i = 0; i < 4; i++) {
		if (knot_element[i + p] == knot_element[i])
			a = 0;
		else
			a = sp0[i] * (zeta - knot_element[i]) / (knot_element[i + p] - knot_element[i]);
		if (knot_element[i + p + 1] == knot_element[i + 1])
			b = 0;
		else
			b = sp0[i + 1] * (knot_element[i + p + 1] - zeta) / (knot_element[i + p + 1] - knot_element[i + 1]);
		sp1[i] = a + b;
	}

	p = 2;
	for (int i = 0; i < 3; i++) {
		if (knot_element[i + p] == knot_element[i])
			a = 0;
		else
			a = sp1[i] * (zeta - knot_element[i]) / (knot_element[i + p] - knot_element[i]);
		if (knot_element[i + p + 1] == knot_element[i + 1])
			b = 0;
		else
			b = sp1[i + 1] * (knot_element[i + p + 1] - zeta) / (knot_element[i + p + 1] - knot_element[i + 1]);
		sp2[i] = a + b;
	}
	if (zeta == 1)
	{
		sp2[0] = 0;
		sp2[1] = 0;
		sp2[2] = 1;
	}

	(*center)(0, 0) = sp2[0] * (*x_Ai)(0, 0) + sp2[1] * (*x_Bi)(0, 0) + sp2[2] * (*x_Ci)(0, 0);
	(*center)(1, 0) = sp2[0] * (*x_Ai)(1, 0) + sp2[1] * (*x_Bi)(1, 0) + sp2[2] * (*x_Ci)(1, 0);
	(*center)(2, 0) = sp2[0] * (*x_Ai)(2, 0) + sp2[1] * (*x_Bi)(2, 0) + sp2[2] * (*x_Ci)(2, 0);

	*radii = *radius;

}

void SplineElement::SaveConfiguration()
{
	d->clear();
	for (int i = 0; i < 3; i++)
	{
		(*dui)(i + 0, 0) = db.nodes[nodes[0] - 1]->copy_vel[i];
		(*dui)(i + 3, 0) = db.nodes[nodes[1] - 1]->copy_vel[i];
		(*dui)(i + 6, 0) = db.nodes[nodes[2] - 1]->copy_vel[i];

		(*ddui)(i + 0, 0) = db.nodes[nodes[0] - 1]->copy_accel[i];
		(*ddui)(i + 3, 0) = db.nodes[nodes[1] - 1]->copy_accel[i];
		(*ddui)(i + 6, 0) = db.nodes[nodes[2] - 1]->copy_accel[i];

		(*x_Ai)(i, 0) = db.nodes[nodes[0] - 1]->copy_coordinates[i];
		(*x_Bi)(i, 0) = db.nodes[nodes[1] - 1]->copy_coordinates[i];
		(*x_Ci)(i, 0) = db.nodes[nodes[2] - 1]->copy_coordinates[i];
	}

	*x_Ap = *x_Ai;
	*x_Bp = *x_Bi;
	*x_Cp = *x_Ci;
}

bool SplineElement::Check()
{
	for (int i = 0; i < n_nodes; i++)
	{
		if (nodes[i] > db.number_nodes)
			return false;
	}
	return true;
}

bool SplineElement::ReadCommon(FILE *f)
{
	char s[1000];

	//Salva a posição (stream)
	fpos_t pos;
	fgetpos(f, &pos);
	fscanf(f, "%s", s);
	fsetpos(f, &pos);

	return true;
}