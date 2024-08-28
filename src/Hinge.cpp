#include "Hinge.h"

#include "InitialCondition.h"
#include "Node.h"
#include "CoordinateSystem.h"
#include "Dynamic.h"

#include"Database.h"
//Variaveis globais
extern
Database db;

#define PI 3.1415926535897932384626433832795

Hinge::Hinge()
{
	n_GL = 5;						//Dois graus de liberdade (esse vinculo possui 5 multiplicadores de lagrange)
	active_lambda = new int[n_GL];
	lambda = new double[n_GL];
	copy_lambda = new double[n_GL];
	GLs = new int[n_GL];
	node_A = 0;
	node_B = 0;
	cs = 0;

	for (int i = 0; i < n_GL; i++)
	{
		active_lambda[i] = 0;
		GLs[i] = 0;
		//Chute inicial para os lambdas: valores nulos
		lambda[i] = 0.0;
		copy_lambda[i] = 0.0;
	}

	I3 = Matrix(3, 3);
	I3(0, 0) = 1.0;
	I3(1, 1) = 1.0;
	I3(2, 2) = 1.0;
	r1 = Matrix(3);

	alphaA = Matrix(3);
	alphaB = Matrix(3);
	alphaiA = Matrix(3);
	alphaiB = Matrix(3);
	A = Matrix(3, 3);
	QA = Matrix(3, 3);
	QB = Matrix(3, 3);
	ei3A = Matrix(3);
	ei1B = Matrix(3);
	ei2B = Matrix(3);
	ei1A = Matrix(3);
	ei2A = Matrix(3);

	//As contribuições estão divididas em duas partes:
	//Parte 1 - exatamente a mesma contribuição do Same Displacement
	//Parte 2 - restrições das rotações
	//Matriz de rigidez tangente e vetor residuo
	stiffness1 = new Matrix(9, 9);
	residual1 = new Matrix(9, 1);
	//Matriz de rigidez tangente e vetor residuo
	stiffness2 = new double*[8];
	for (int i = 0; i < 8; i++)
		stiffness2[i] = new double[8];
	residual2 = new double[8];
	temp_lambda = new double[2];

	//Matriz de rigidez tangente e vetor residuo - contribuição da mola (spring)
	stiffness_spring = new double*[6];
	for (int i = 0; i < 6; i++)
		stiffness_spring[i] = new double[6];
	residual_damper = new double[6];
	residual_spring = new double[6];
	//Amortecedor de torção
	stiffness_damper = new double*[6];
	for (int i = 0; i < 6; i++)
		stiffness_damper[i] = new double[6];
	for (int i = 0; i < 6; i++)
	{
		residual_damper[i] = 0.0;
		residual_spring[i] = 0.0;
		for (int j = 0; j < 6; j++)
		{
			stiffness_damper[i][j] = 0.0;
			stiffness_spring[i][j] = 0.0;
		}
	}
	for (int i = 0; i < 8; i++)
	{
		residual2[i] = 0.0;
		for (int j = 0; j < 8; j++)
			stiffness2[i][j] = 0.0;
	}

	thetai = 0.0;
	thetad = 0.0;
	stiffc = 0.0;
	dampc1 = 0.0;
	dampc2 = 0.0;

	omegaiA = Matrix(3);
	omegaiB = Matrix(3);
	domegaiA = Matrix(3);
	domegaiB = Matrix(3);
}

//Zera matrizes e vetores
void Hinge::ClearContributions()
{
	for (int i = 0; i < 6; i++)
	{
		residual_damper[i] = 0.0;
		residual_spring[i] = 0.0;
		for (int j = 0; j < 6; j++)
		{
			stiffness_damper[i][j] = 0.0;
			stiffness_spring[i][j] = 0.0;
		}
	}
	for (int i = 0; i < 8; i++)
	{
		residual2[i] = 0.0;
		for (int j = 0; j < 8; j++)
			stiffness2[i][j] = 0.0;
	}
}

Hinge::~Hinge()
{
	delete[]active_lambda;
	delete[]GLs;
	delete[]lambda;
	delete[]copy_lambda;
	delete stiffness1;
	delete residual1;
	
	for (int i = 0; i < 8; i++)
		delete []stiffness2[i];
	delete []stiffness2;
	delete []residual2;
	delete []temp_lambda;

	for (int i = 0; i < 6; i++)
		delete[]stiffness_spring[i];
	delete[]stiffness_spring;
	delete[]residual_spring;

	for (int i = 0; i < 6; i++)
		delete[]stiffness_damper[i];
	delete[]stiffness_damper;
	delete[]residual_damper;
}

//Leitura
bool Hinge::Read(FILE *f)
{
	char s[1000];
	fscanf(f, "%s", s);
	number = atoi(s);

	fscanf(f, "%s", s);
	if (!strcmp(s, "Nodes"))
	{
		fscanf(f, "%s", s);
		node_A = atoi(s);

		fscanf(f, "%s", s);
		node_B = atoi(s);
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "CS"))
	{
		fscanf(f, "%s", s);
		cs = atoi(s);
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "LinearStiffness"))
	{
		fscanf(f, "%s", s);
		stiffc = atof(s);
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "LinearDamping"))
	{
		fscanf(f, "%s", s);
		dampc1 = atof(s);
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "QuadraticDamping"))
	{
		fscanf(f, "%s", s);
		dampc2 = atof(s);
	}
	else
		return false;
	//Salva a posição (stream)
	fpos_t pos;
	fgetpos(f, &pos);
	fscanf(f, "%s", s);
	if (!strcmp(s, "BoolTable"))
		bool_table.Read(f);
	else
	{
		fsetpos(f, &pos);
		bool_table.SetDefault(true);
	}
		

	return true;
}

//Gravação
void Hinge::Write(FILE *f)
{
	fprintf(f, "HingeJoint\t%d\tNodes\t%d\t%d\tCS\t%d\tLinearStiffness\t%.6e\tLinearDamping\t%.6e\tQuadraticDamping\t%.6e\n",
		number,
		node_A,
		node_B,
		cs,
		stiffc,
		dampc1,
		dampc2);
}

//Escreve no monitor do SpecialConstraint
void Hinge::WriteMonitor(FILE *f, bool first_record, double time)
{

}

void Hinge::WriteVTK_XMLRender(FILE *f)
{
	//if (data.post_files->WriteSpecialConstraints_flag == true)
	//{
	//	int i, offsets;			            /*numero de pontos*/
	//	int tamanho, n;						/*tamanho do vetor - numero de repartições do cilindro*/
	//	double xa, ya, za, xb, yb, zb;       /*pontos extremidade do vetor da intercecção*/
	//	double x1, y1, z1, x2, y2, z2;       /*pontos de cada uma das superficies*/
	//	double *x, *y, *z;					/*pontos*/
	//	double d, R, r;                      /*Comprimento da intercecção e raio*/
	//	double alfa, beta, teta;             /*angulos rotação nos eixos x e y - mudança de coordenadas*/
	//	double *xg, *yg, *zg;				/*pontos no sistema global*/
	//	double dA, dB, tA1, tA2, tB1, tB2;
	//	double xa1, ya1, za1, xa2, ya2, za2, xb1, yb1, zb1, xb2, yb2, zb2;
	//	double xv1, yv1, zv1, xv2, yv2, zv2, xv3, yv3, zv3, xv4, yv4, zv4;           /*vertices placas*/
	//	double xv5, yv5, zv5, xv6, yv6, zv6, xv7, yv7, zv7, xv8, yv8, zv8;
	//	double xv9, yv9, zv9, xv10, yv10, zv10, xv11, yv11, zv11, xv12, yv12, zv12;
	//	double xv13, yv13, zv13, xv14, yv14, zv14, xv15, yv15, zv15, xv16, yv16, zv16;
	//	double a, b, c;


	//	for (int i = 0; i < 3; i++)
	//	{
	//		alphaiA(i, 0) = data.nodes[node_A - 1]->copy_coordinates[i + 3];	//vetor rotação acumulada (do inicio) do nó A
	//		alphaiB(i, 0) = data.nodes[node_B - 1]->copy_coordinates[i + 3];	//vetor rotação acumulada (do inicio) do nó B
	//	}

	//	alpha_escalar_i = norm(alphaiA);
	//	A = skew(alphaiA);
	//	g = 4.0 / (4.0 + alpha_escalar_i*alpha_escalar_i);
	//	QA = I3 + g*(A + 0.5*(A*A));
	//	ei1A = QA*(*data.CS[cs - 1]->E1);	//Eixo e1 no inicio do incremento
	//	ei2A = QA*(*data.CS[cs - 1]->E2);	//Eixo e2 no inicio do incremento
	//	ei3A = QA*(*data.CS[cs - 1]->E3);	//Eixo e2 no inicio do incremento

	//	alpha_escalar_i = norm(alphaiB);
	//	A = skew(alphaiB);
	//	g = 4.0 / (4.0 + alpha_escalar_i*alpha_escalar_i);
	//	QB = I3 + g*(A + 0.5*(A*A));
	//	ei1B = QB*(*data.CS[cs - 1]->E1);	//Eixo e1 no inicio do incremento
	//	ei2B = QB*(*data.CS[cs - 1]->E2);	//Eixo e2 no inicio do incremento
	//	double len = EvaluateBoundingBoxDiag() / data.number_nodes;
	//	xa = data.nodes[node_A - 1]->copy_coordinates[0] - ei3A(0, 0)*len;
	//	ya = data.nodes[node_A - 1]->copy_coordinates[1] - ei3A(1, 0)*len;
	//	za = data.nodes[node_A - 1]->copy_coordinates[2] - ei3A(2, 0)*len;
	//	xb = data.nodes[node_B - 1]->copy_coordinates[0] + ei3A(0, 0)*len;
	//	yb = data.nodes[node_B - 1]->copy_coordinates[1] + ei3A(1, 0)*len;
	//	zb = data.nodes[node_B - 1]->copy_coordinates[2] + ei3A(2, 0)*len;

	//	x1 = (ei1A(0, 0) + ei2A(0, 0))*len;
	//	y1 = (ei1A(1, 0) + ei2A(1, 0))*len;
	//	z1 = (ei1A(2, 0) + ei2A(2, 0))*len;
	//	x2 = (ei1B(0, 0) + ei2B(0, 0))*len;
	//	y2 = (ei1B(1, 0) + ei2B(1, 0))*len;
	//	z2 = (ei1B(2, 0) + ei2B(2, 0))*len;

	//	n = 40;
	//	tamanho = 2 * n;
	//	x = new double[tamanho];
	//	y = new double[tamanho];
	//	z = new double[tamanho];
	//	xg = new double[tamanho];
	//	yg = new double[tamanho];
	//	zg = new double[tamanho];

	//	d = sqrt((xb - xa)*(xb - xa) + (yb - ya)*(yb - ya) + (zb - za)*(zb - za));

	//	R = 0.1*d;
	//	r = 0.1*R;

	//	teta = 0;
	//	for (i = 0; i < n; i++)				/*n repartições*/
	//	{
	//		x[i] = R*(float)cos(teta);
	//		teta = teta + (3.14159265 / (n / 2));			/*9º em radianos*/
	//	}
	//	teta = 0;
	//	for (i = n; i < 2 * n; i++)				/*n repartições*/
	//	{
	//		x[i] = R*(float)cos(teta);
	//		teta = teta + (3.14159265 / (n / 2));			/*9º em radianos*/
	//	}

	//	teta = 0;
	//	for (i = 0; i < n; i++)
	//	{
	//		y[i] = R*(float)sin(teta);
	//		teta = teta + (3.14159265 / (n / 2)); /*9º em radianos*/
	//	}
	//	teta = 0;
	//	for (i = n; i < 2 * n; i++)
	//	{
	//		y[i] = R*(float)sin(teta);
	//		teta = teta + (3.14159265 / (n / 2)); /*9º em radianos*/
	//	}

	//	for (i = 0; i < n; i++)
	//	{
	//		z[i] = 0;
	//	}
	//	for (i = n; i < 2 * n; i++)
	//	{
	//		z[i] = d;
	//	}

	//	/*Mudança de coordenadas*/
	//	alfa = atan2((ya - yb), (zb - za));
	//	beta = asin((xb - xa) / d);
	//	for (i = 0; i < 2 * n; i++)
	//	{
	//		xg[i] = xa + x[i] * cos(beta) + z[i] * sin(beta);
	//		yg[i] = ya + x[i] * sin(alfa)*sin(beta) + y[i] * cos(alfa) - z[i] * sin(alfa)*cos(beta);
	//		zg[i] = za - x[i] * cos(alfa)*sin(beta) + y[i] * sin(alfa) + z[i] * cos(alfa)*cos(beta);
	//	}


	//	/*Plano ortogonal a AB em A*/ /*Geometria Analitica, Paulo Boulos, Ex resolvido 17-8*/
	//	dA = -((xb - xa)*xa + (yb - ya)*ya + (zb - za)*za);

	//	tA1 = (-dA - ((xb - xa)*x1 + (yb - ya)*y1 + (zb - za)*z1)) / ((xb - xa)*(xb - xa) + (yb - ya)*(yb - ya) + (zb - za)*(zb - za));
	//	xa1 = tA1*(xb - xa) + x1;
	//	ya1 = tA1*(yb - ya) + y1;
	//	za1 = tA1*(zb - za) + z1;

	//	tA2 = (-dA - ((xb - xa)*x2 + (yb - ya)*y2 + (zb - za)*z2)) / ((xb - xa)*(xb - xa) + (yb - ya)*(yb - ya) + (zb - za)*(zb - za));
	//	xa2 = tA2*(xb - xa) + x2;
	//	ya2 = tA2*(yb - ya) + y2;
	//	za2 = tA2*(zb - za) + z2;

	//	/*Plano ortogonal a* AB em B*/
	//	dB = -((xb - xa)*xb + (yb - ya)*yb + (zb - za)*zb);

	//	tB1 = (-dB - ((xb - xa)*x1 + (yb - ya)*y1 + (zb - za)*z1)) / ((xb - xa)*(xb - xa) + (yb - ya)*(yb - ya) + (zb - za)*(zb - za));
	//	xb1 = tB1*(xb - xa) + x1;
	//	yb1 = tB1*(yb - ya) + y1;
	//	zb1 = tB1*(zb - za) + z1;

	//	tB2 = (-dB - ((xb - xa)*x2 + (yb - ya)*y2 + (zb - za)*z2)) / ((xb - xa)*(xb - xa) + (yb - ya)*(yb - ya) + (zb - za)*(zb - za));
	//	xb2 = tB2*(xb - xa) + x2;
	//	yb2 = tB2*(yb - ya) + y2;
	//	zb2 = tB2*(zb - za) + z2;

	//	/*Placa - vertices A B B1 A1*/
	//	a = ((z1 - za)*(yb - ya) - (zb - za)*(y1 - ya));
	//	b = ((zb - za)*(x1 - xa) - (xb - xa)*(z1 - za));
	//	c = ((xb - xa)*(y1 - ya) - (yb - ya)*(x1 - xa));

	//	/*Placa - vertices A B*/

	//	xv1 = xa + r*a / sqrt(a*a + b*b + c*c);
	//	yv1 = ya + r*b / sqrt(a*a + b*b + c*c);
	//	zv1 = za + r*c / sqrt(a*a + b*b + c*c);

	//	xv2 = xa - r*a / sqrt(a*a + b*b + c*c);
	//	yv2 = ya - r*b / sqrt(a*a + b*b + c*c);
	//	zv2 = za - r*c / sqrt(a*a + b*b + c*c);

	//	xv3 = xb + r*a / sqrt(a*a + b*b + c*c);
	//	yv3 = yb + r*b / sqrt(a*a + b*b + c*c);
	//	zv3 = zb + r*c / sqrt(a*a + b*b + c*c);

	//	xv4 = xb - r*a / sqrt(a*a + b*b + c*c);
	//	yv4 = yb - r*b / sqrt(a*a + b*b + c*c);
	//	zv4 = zb - r*c / sqrt(a*a + b*b + c*c);

	//	/*Placa - vertices B1 A1*/

	//	xv5 = xa1 + r*a / sqrt(a*a + b*b + c*c);
	//	yv5 = ya1 + r*b / sqrt(a*a + b*b + c*c);
	//	zv5 = za1 + r*c / sqrt(a*a + b*b + c*c);

	//	xv6 = xa1 - r*a / sqrt(a*a + b*b + c*c);
	//	yv6 = ya1 - r*b / sqrt(a*a + b*b + c*c);
	//	zv6 = za1 - r*c / sqrt(a*a + b*b + c*c);

	//	xv7 = xb1 + r*a / sqrt(a*a + b*b + c*c);
	//	yv7 = yb1 + r*b / sqrt(a*a + b*b + c*c);
	//	zv7 = zb1 + r*c / sqrt(a*a + b*b + c*c);

	//	xv8 = xb1 - r*a / sqrt(a*a + b*b + c*c);
	//	yv8 = yb1 - r*b / sqrt(a*a + b*b + c*c);
	//	zv8 = zb1 - r*c / sqrt(a*a + b*b + c*c);

	//	/*Placa - vertices A B B2 A2*/
	//	/*Placa - vertices A B*/
	//	a = ((z2 - za)*(yb - ya) - (zb - za)*(y2 - ya));
	//	b = ((zb - za)*(x2 - xa) - (xb - xa)*(z2 - za));
	//	c = ((xb - xa)*(y2 - ya) - (yb - ya)*(x2 - xa));

	//	xv13 = xa + r*a / sqrt(a*a + b*b + c*c);
	//	yv13 = ya + r*b / sqrt(a*a + b*b + c*c);
	//	zv13 = za + r*c / sqrt(a*a + b*b + c*c);

	//	xv14 = xa - r*a / sqrt(a*a + b*b + c*c);
	//	yv14 = ya - r*b / sqrt(a*a + b*b + c*c);
	//	zv14 = za - r*c / sqrt(a*a + b*b + c*c);

	//	xv15 = xb + r*a / sqrt(a*a + b*b + c*c);
	//	yv15 = yb + r*b / sqrt(a*a + b*b + c*c);
	//	zv15 = zb + r*c / sqrt(a*a + b*b + c*c);

	//	xv16 = xb - r*a / sqrt(a*a + b*b + c*c);
	//	yv16 = yb - r*b / sqrt(a*a + b*b + c*c);
	//	zv16 = zb - r*c / sqrt(a*a + b*b + c*c);

	//	/*Placa - vertices B2 A2*/
	//	xv9 = xa2 + r*a / sqrt(a*a + b*b + c*c);
	//	yv9 = ya2 + r*b / sqrt(a*a + b*b + c*c);
	//	zv9 = za2 + r*c / sqrt(a*a + b*b + c*c);

	//	xv10 = xa2 - r*a / sqrt(a*a + b*b + c*c);
	//	yv10 = ya2 - r*b / sqrt(a*a + b*b + c*c);
	//	zv10 = za2 - r*c / sqrt(a*a + b*b + c*c);

	//	xv11 = xb2 + r*a / sqrt(a*a + b*b + c*c);
	//	yv11 = yb2 + r*b / sqrt(a*a + b*b + c*c);
	//	zv11 = zb2 + r*c / sqrt(a*a + b*b + c*c);

	//	xv12 = xb2 - r*a / sqrt(a*a + b*b + c*c);
	//	yv12 = yb2 - r*b / sqrt(a*a + b*b + c*c);
	//	zv12 = zb2 - r*c / sqrt(a*a + b*b + c*c);

	//	/*Piece 1 - Cilindro da intercecão*/
	//	fprintf(f, "     <Piece NumberOfPoints = \"%d\" NumberOfCells = \"%d\">\n", 2 * n, n + 2);
	//	fprintf(f, "         <Points>\n");
	//	fprintf(f, "             <DataArray type = \"Float32\" NumberOfComponents = \"3\" Format = \"ascii\">\n");
	//	for (i = 0; i < 2 * n; i++)
	//	{
	//		fprintf(f, "                 %.6f %.6f %.6f\n", xg[i], yg[i], zg[i]);
	//	}
	//	fprintf(f, "             </DataArray>\n");
	//	fprintf(f, "         </Points>\n");
	//	fprintf(f, "         <Cells>\n");
	//	fprintf(f, "             <DataArray type = \"Int32\" Name = \"connectivity\" Format = \"ascii\">\n");
	//	for (i = 0; i < n - 1; i++)
	//	{
	//		fprintf(f, "                  %d %d %d %d\n", i, i + 1, i + 1 + n, i + n);
	//	}
	//	fprintf(f, "                  %d %d %d %d\n", n - 1, 0, n, 2 * n - 1);
	//	fprintf(f, "                  ");
	//	for (i = 0; i < n; i++)
	//	{
	//		fprintf(f, "%d ", i);
	//	}
	//	fprintf(f, "                  \n");
	//	fprintf(f, "                  ");
	//	for (i = n; i < 2 * n; i++)
	//	{
	//		fprintf(f, "%d ", i);
	//	}
	//	fprintf(f, "                  \n");
	//	fprintf(f, "             </DataArray>\n");
	//	fprintf(f, "             <DataArray type = \"Int32\" Name = \"types\" Format = \"ascii\">\n");
	//	for (i = 0; i < n + 2; i++)
	//	{
	//		fprintf(f, "                 7\n");
	//	}
	//	fprintf(f, "             </DataArray>\n");
	//	fprintf(f, "             <DataArray type = \"Int32\" Name = \"offsets\" Format = \"ascii\">\n");
	//	offsets = 0;
	//	for (i = 0; i < n; i++)
	//	{
	//		offsets = offsets + 4;
	//		fprintf(f, "                 %d\n", offsets);
	//	}
	//	fprintf(f, "                 %d\n", offsets + n);
	//	fprintf(f, "                 %d\n", offsets + 2 * n);
	//	fprintf(f, "             </DataArray>\n");
	//	fprintf(f, "         </Cells>\n");
	//	fprintf(f, "     </Piece>\n");

	//	/*Piece 2 - Placa 1*/
	//	fprintf(f, "     <Piece NumberOfPoints = \"%d\" NumberOfCells = \"%d\">\n", 8, 5);
	//	fprintf(f, "         <Points>\n");
	//	fprintf(f, "             <DataArray type = \"Float32\" NumberOfComponents = \"3\" Format = \"ascii\">\n");
	//	fprintf(f, "                 %.6f %.6f %.6f\n", xv1, yv1, zv1);
	//	fprintf(f, "                 %.6f %.6f %.6f\n", xv2, yv2, zv2);
	//	fprintf(f, "                 %.6f %.6f %.6f\n", xv3, yv3, zv3);
	//	fprintf(f, "                 %.6f %.6f %.6f\n", xv4, yv4, zv4);
	//	fprintf(f, "                 %.6f %.6f %.6f\n", xv5, yv5, zv5);
	//	fprintf(f, "                 %.6f %.6f %.6f\n", xv6, yv6, zv6);
	//	fprintf(f, "                 %.6f %.6f %.6f\n", xv7, yv7, zv7);
	//	fprintf(f, "                 %.6f %.6f %.6f\n", xv8, yv8, zv8);
	//	fprintf(f, "             </DataArray>\n");
	//	fprintf(f, "         </Points>\n");
	//	fprintf(f, "         <Cells>\n");
	//	fprintf(f, "             <DataArray type = \"Int32\" Name = \"connectivity\" Format = \"ascii\">\n");
	//	fprintf(f, "                  %d %d %d %d\n", 0, 1, 5, 4);
	//	fprintf(f, "                  %d %d %d %d\n", 2, 3, 7, 6);
	//	fprintf(f, "                  %d %d %d %d\n", 4, 5, 7, 6);
	//	fprintf(f, "                  %d %d %d %d\n", 0, 4, 6, 2);
	//	fprintf(f, "                  %d %d %d %d\n", 1, 5, 7, 3);
	//	fprintf(f, "             </DataArray>\n");
	//	fprintf(f, "             <DataArray type = \"Int32\" Name = \"types\" Format = \"ascii\">\n");
	//	for (i = 0; i < 5; i++)
	//	{
	//		fprintf(f, "                 7\n");
	//	}
	//	fprintf(f, "             </DataArray>\n");
	//	fprintf(f, "             <DataArray type = \"Int32\" Name = \"offsets\" Format = \"ascii\">\n");
	//	offsets = 0;
	//	for (i = 0; i < 5; i++)
	//	{
	//		offsets = offsets + 4;
	//		fprintf(f, "                 %d\n", offsets);
	//	}
	//	fprintf(f, "             </DataArray>\n");
	//	fprintf(f, "         </Cells>\n");
	//	fprintf(f, "     </Piece>\n");

	//	/*Piece 3 - Placa 2*/
	//	fprintf(f, "     <Piece NumberOfPoints = \"%d\" NumberOfCells = \"%d\">\n", 8, 5);
	//	fprintf(f, "         <Points>\n");
	//	fprintf(f, "             <DataArray type = \"Float32\" NumberOfComponents = \"3\" Format = \"ascii\">\n");
	//	fprintf(f, "                 %.6f %.6f %.6f\n", xv13, yv13, zv13);
	//	fprintf(f, "                 %.6f %.6f %.6f\n", xv14, yv14, zv14);
	//	fprintf(f, "                 %.6f %.6f %.6f\n", xv15, yv15, zv15);
	//	fprintf(f, "                 %.6f %.6f %.6f\n", xv16, yv16, zv16);
	//	fprintf(f, "                 %.6f %.6f %.6f\n", xv9, yv9, zv9);
	//	fprintf(f, "                 %.6f %.6f %.6f\n", xv10, yv10, zv10);
	//	fprintf(f, "                 %.6f %.6f %.6f\n", xv11, yv11, zv11);
	//	fprintf(f, "                 %.6f %.6f %.6f\n", xv12, yv12, zv12);
	//	fprintf(f, "             </DataArray>\n");
	//	fprintf(f, "         </Points>\n");
	//	fprintf(f, "         <Cells>\n");
	//	fprintf(f, "             <DataArray type = \"Int32\" Name = \"connectivity\" Format = \"ascii\">\n");
	//	fprintf(f, "                  %d %d %d %d\n", 0, 1, 5, 4);
	//	fprintf(f, "                  %d %d %d %d\n", 2, 3, 7, 6);
	//	fprintf(f, "                  %d %d %d %d\n", 4, 5, 7, 6);
	//	fprintf(f, "                  %d %d %d %d\n", 0, 4, 6, 2);
	//	fprintf(f, "                  %d %d %d %d\n", 1, 5, 7, 3);
	//	fprintf(f, "             </DataArray>\n");
	//	fprintf(f, "             <DataArray type = \"Int32\" Name = \"types\" Format = \"ascii\">\n");
	//	for (i = 0; i < 5; i++)
	//	{
	//		fprintf(f, "                 7\n");
	//	}
	//	fprintf(f, "             </DataArray>\n");
	//	fprintf(f, "             <DataArray type = \"Int32\" Name = \"offsets\" Format = \"ascii\">\n");
	//	offsets = 0;
	//	for (i = 0; i < 5; i++)
	//	{
	//		offsets = offsets + 4;
	//		fprintf(f, "                 %d\n", offsets);
	//	}
	//	fprintf(f, "             </DataArray>\n");
	//	fprintf(f, "         </Cells>\n");
	//	fprintf(f, "     </Piece>\n");

	//	delete[]x;
	//	delete[]y;
	//	delete[]z;
	//	delete[]xg;
	//	delete[]yg;
	//	delete[]zg;
	//}
}

//Checa inconsistências no SC para evitar erros de execução
bool Hinge::Check()
{
	if (node_A > db.number_nodes)
		return false;
	if (node_B > db.number_nodes)
		return false;
	if (cs > db.number_CS)
		return false;

	//Checagem das condições iniciais
	int temp_node = 0;
	for (int i = 0; i < db.number_IC; i++)
	{
		temp_node = db.IC[i]->node;
		if (node_B == temp_node)
		{
			db.myprintf("Warning in Special Constraint %d.\nInitial Condition %d was prescribed to node %d (slave), leading to ignoring some of its components.\n", number, db.IC[i]->number, db.IC[i]->node);
		}
	}
	return true;
}

//Montagem dos residuos e rigidez tangente
void Hinge::Mount()
{
	ClearContributions();
	//Nesse vinculo a rigidez tangente da parte 1 não se modifica nunca. e montada no PreCalc()
	//Montagem do residuo da parte 1:
	if (active_lambda[0] == 1 && active_lambda[1] == 1 && active_lambda[2] == 1)
	{
		for (int i = 0; i < 3; i++)
			r1(i, 0) = (db.nodes[node_A - 1]->copy_coordinates[i] - db.nodes[node_A - 1]->ref_coordinates[i] + db.nodes[node_A - 1]->displacements[i]) -
			(db.nodes[node_B - 1]->copy_coordinates[i] - db.nodes[node_B - 1]->ref_coordinates[i] + db.nodes[node_B - 1]->displacements[i]);

		//Atualização do vetor de residuos
		(*residual1)(0, 0) = lambda[0];
		(*residual1)(1, 0) = lambda[1];
		(*residual1)(2, 0) = lambda[2];

		(*residual1)(3, 0) = -lambda[0];
		(*residual1)(4, 0) = -lambda[1];
		(*residual1)(5, 0) = -lambda[2];

		(*residual1)(6, 0) = r1(0, 0);
		(*residual1)(7, 0) = r1(1, 0);
		(*residual1)(8, 0) = r1(2, 0);
	}
	
	//Montagem da rigidez tangente e residuo da parte 2:
	if (active_lambda[3] == 1 && active_lambda[4] == 1)
	{
		for (int i = 0; i < 3; i++)
		{
			alphaA(i, 0) = db.nodes[node_A - 1]->displacements[i + 3];		//vetor rotação (atual) do nó A
			alphaB(i, 0) = db.nodes[node_B - 1]->displacements[i + 3];		//vetor rotação (atual) do nó B
			alphaiA(i, 0) = db.nodes[node_A - 1]->copy_coordinates[i + 3];	//vetor rotação acumulada (do inicio) do nó A
			alphaiB(i, 0) = db.nodes[node_B - 1]->copy_coordinates[i + 3];	//vetor rotação acumulada (do inicio) do nó B
		}
		alpha_escalar_i = norm(alphaiA);
		A = skew(alphaiA);
		g = 4.0 / (4.0 + alpha_escalar_i*alpha_escalar_i);
		QA = I3 + g*(A + 0.5*(A*A));
		ei3A = QA*(*db.CS[cs - 1]->E3);	//Eixo e3 no inicio do incremento

		alpha_escalar_i = norm(alphaiB);
		A = skew(alphaiB);
		g = 4.0 / (4.0 + alpha_escalar_i*alpha_escalar_i);
		QB = I3 + g*(A + 0.5*(A*A));
		ei1B = QB*(*db.CS[cs - 1]->E1);	//Eixo e1 no inicio do incremento
		ei2B = QB*(*db.CS[cs - 1]->E2);	//Eixo e2 no inicio do incremento

		//printf("cos theta = %lf\n", dot(ei1A, ei1B));

		temp_lambda[0] = lambda[3];
		temp_lambda[1] = lambda[4];
		EvaluateHingeContribution2(temp_v, residual2, stiffness2, alphaA.getMatrix(), alphaB.getMatrix(), ei3A.getMatrix(), ei1B.getMatrix(), ei2B.getMatrix(), temp_lambda);

		//Contribuição do mola de torção (spring)
		ei1A = QA*(*db.CS[cs - 1]->E1);	//Eixo e1 no inicio do incremento
		ei2A = QA*(*db.CS[cs - 1]->E2);	//Eixo e2 no inicio do incremento
		EvaluateTorsionSpring(temp_v, residual_spring, stiffness_spring, alphaA.getMatrix(), alphaB.getMatrix(), ei1A.getMatrix(), ei1B.getMatrix(), ei3A.getMatrix(), &stiffc, &thetai, &thetad);

		for (int i = 0; i < 3; i++)
		{
			omegaiA(i, 0) = db.nodes[node_A - 1]->copy_vel[i + 3];
			omegaiB(i, 0) = db.nodes[node_B - 1]->copy_vel[i + 3];
			domegaiA(i, 0) = db.nodes[node_A - 1]->copy_accel[i + 3];
			domegaiB(i, 0) = db.nodes[node_B - 1]->copy_accel[i + 3];
		}
		if (typeid(*db.solution[db.current_solution_number - 1]) == typeid(Dynamic))
		{
			Dynamic* ptr_sol = static_cast<Dynamic*>(db.solution[db.current_solution_number - 1]);
			EvaluateTorsionDamping(temp_v, residual_damper, stiffness_damper, alphaA.getMatrix(), alphaB.getMatrix(), omegaiA.getMatrix(), omegaiB.getMatrix(), domegaiA.getMatrix(), domegaiB.getMatrix(), ei3A.getMatrix(), &dampc1, &dampc2, &ptr_sol->a4, &ptr_sol->a5, &ptr_sol->a6);
		}
	}
}

//Preenche a contribuição do elemento nas matrizes globais
void Hinge::MountGlobal()
{
	//Variaveis temporarias para salvar a indexação global dos graus de liberdade a serem setados na matriz de rigidez global
	int GL_global_1 = 0;
	int GL_global_2 = 0;
	double anterior = 0;
	
	//PARTE 1 - deslocamentos
	if (active_lambda[0] == 1 && active_lambda[1] == 1 && active_lambda[2] == 1)
	{
		for (int i = 0; i < 9; i++)
		{
			//////////////MONTAGEM DO VETOR DE ESFORÇOS DESBALANCEADOS//////////////////
			//Toma do vetor de GL globais, a indexação de cada grau de liberdade global
			if (i < 3)//uA
				GL_global_1 = db.nodes[node_A - 1]->GLs[i];
			else
			{
				if (i<6)//uB
					GL_global_1 = db.nodes[node_B - 1]->GLs[i - 3];
				else
				{
					//lambda
					GL_global_1 = GLs[i - 6];
				}
			}

			//Caso o grau de liberdade seja livre:
			if (GL_global_1 > 0)
			{
				anterior = db.global_P_A(GL_global_1 - 1, 0);
				db.global_P_A(GL_global_1 - 1, 0) = anterior + (*residual1)(i, 0);
				anterior = db.global_I_A(GL_global_1 - 1, 0);
				db.global_I_A(GL_global_1 - 1, 0) = anterior + (*residual1)(i, 0);
			}
			else
			{
				anterior = db.global_P_B(-GL_global_1 - 1, 0);
				db.global_P_B(-GL_global_1 - 1, 0) = anterior + (*residual1)(i, 0);
			}
			for (int j = 0; j < 9; j++)
			{
				//////////////////////MONTAGEM DA MATRIZ DE RIGIDEZ/////////////////////////
				//Toma do vetor de GL globais, a indexação de cada grau de liberdade global
				if (j < 3)//uA
					GL_global_2 = db.nodes[node_A - 1]->GLs[j];
				else
				{
					if (j<6)//uB
						GL_global_2 = db.nodes[node_B - 1]->GLs[j - 3];
					else
					{
						//lambda
						GL_global_2 = GLs[j - 6];
					}
				}

				//Caso os graus de liberdade sejam ambos livres (Matriz Kaa)
				if (GL_global_1 > 0 && GL_global_2 > 0)
					db.global_stiffness_AA.setValue(GL_global_1 - 1, GL_global_2 - 1, (*stiffness1)(i, j));
				//Caso os graus de liberdade sejam ambos fixos (Matriz Kbb)
				if (GL_global_1 < 0 && GL_global_2 < 0)
					db.global_stiffness_BB.setValue(-GL_global_1 - 1, -GL_global_2 - 1, (*stiffness1)(i, j));
				//Caso os graus de liberdade sejam livre e fixo (Matriz Kab)
				if (GL_global_1 > 0 && GL_global_2 < 0)
					db.global_stiffness_AB.setValue(GL_global_1 - 1, -GL_global_2 - 1, (*stiffness1)(i, j));
				//Caso os graus de liberdade sejam fixo e livre (Matriz Kba)
				if (GL_global_1 < 0 && GL_global_2 > 0)
					db.global_stiffness_BA.setValue(-GL_global_1 - 1, GL_global_2 - 1, (*stiffness1)(i, j));
			}
		}
	}

	
	if (active_lambda[3] == 1 && active_lambda[4] == 1)
	{
		//PARTE 2 - rotações
		for (int i = 0; i < 8; i++)
		{
			//////////////MONTAGEM DO VETOR DE ESFORÇOS DESBALANCEADOS//////////////////
			//Toma do vetor de GL globais, a indexação de cada grau de liberdade global
			if (i < 3)//alpha_A -> i=0,1,2
				GL_global_1 = db.nodes[node_A - 1]->GLs[i + 3];
			else
			{
				if (i<6)//alpha_B -> i=3,4,5
					GL_global_1 = db.nodes[node_B - 1]->GLs[i];
				else
				{
					//lambda  -> i=6,7
					GL_global_1 = GLs[i - 3];
				}
			}

			//Caso o grau de liberdade seja livre:
			if (GL_global_1 > 0)
			{
				anterior = db.global_P_A(GL_global_1 - 1, 0);
				db.global_P_A(GL_global_1 - 1, 0) = anterior + residual2[i];
				anterior = db.global_I_A(GL_global_1 - 1, 0);
				db.global_I_A(GL_global_1 - 1, 0) = anterior + residual2[i];
			}
			else
			{
				anterior = db.global_P_B(-GL_global_1 - 1, 0);
				db.global_P_B(-GL_global_1 - 1, 0) = anterior + residual2[i];
			}
			for (int j = 0; j < 8; j++)
			{
				//////////////////////MONTAGEM DA MATRIZ DE RIGIDEZ/////////////////////////
				//Toma do vetor de GL globais, a indexação de cada grau de liberdade global
				if (j < 3)//alpha_A -> j=0,1,2
					GL_global_2 = db.nodes[node_A - 1]->GLs[j + 3];
				else
				{
					if (j<6)//alpha_B -> j=3,4,5
						GL_global_2 = db.nodes[node_B - 1]->GLs[j];
					else
					{
						//lambda  -> i=6,7
						GL_global_2 = GLs[j - 3];
					}
				}

				//Caso os graus de liberdade sejam ambos livres (Matriz Kaa)
				if (GL_global_1 > 0 && GL_global_2 > 0)
					db.global_stiffness_AA.setValue(GL_global_1 - 1, GL_global_2 - 1, stiffness2[i][j]);
				//Caso os graus de liberdade sejam ambos fixos (Matriz Kbb)
				if (GL_global_1 < 0 && GL_global_2 < 0)
					db.global_stiffness_BB.setValue(-GL_global_1 - 1, -GL_global_2 - 1, stiffness2[i][j]);
				//Caso os graus de liberdade sejam livre e fixo (Matriz Kab)
				if (GL_global_1 > 0 && GL_global_2 < 0)
					db.global_stiffness_AB.setValue(GL_global_1 - 1, -GL_global_2 - 1, stiffness2[i][j]);
				//Caso os graus de liberdade sejam fixo e livre (Matriz Kba)
				if (GL_global_1 < 0 && GL_global_2 > 0)
					db.global_stiffness_BA.setValue(-GL_global_1 - 1, GL_global_2 - 1, stiffness2[i][j]);
			}
		}

		//PARTE 3 - spring
		for (int i = 0; i < 6; i++)
		{
			//////////////MONTAGEM DO VETOR DE ESFORÇOS DESBALANCEADOS//////////////////
			//Toma do vetor de GL globais, a indexação de cada grau de liberdade global
			if (i < 3)//alpha_A -> i=0,1,2
				GL_global_1 = db.nodes[node_A - 1]->GLs[i + 3];
			else//alpha_B -> i=3,4,5
				GL_global_1 = db.nodes[node_B - 1]->GLs[i];

			//Caso o grau de liberdade seja livre:
			if (GL_global_1 > 0)
			{
				anterior = db.global_P_A(GL_global_1 - 1, 0);
				db.global_P_A(GL_global_1 - 1, 0) = anterior + residual_spring[i];
				anterior = db.global_I_A(GL_global_1 - 1, 0);
				db.global_I_A(GL_global_1 - 1, 0) = anterior + residual_spring[i];
			}
			else
			{
				anterior = db.global_P_B(-GL_global_1 - 1, 0);
				db.global_P_B(-GL_global_1 - 1, 0) = anterior + residual_spring[i];
			}
			for (int j = 0; j < 6; j++)
			{
				//////////////////////MONTAGEM DA MATRIZ DE RIGIDEZ/////////////////////////
				//Toma do vetor de GL globais, a indexação de cada grau de liberdade global
				if (j < 3)//alpha_A -> j=0,1,2
					GL_global_2 = db.nodes[node_A - 1]->GLs[j + 3];
				else//alpha_B -> j=3,4,5
					GL_global_2 = db.nodes[node_B - 1]->GLs[j];

				//Caso os graus de liberdade sejam ambos livres (Matriz Kaa)
				if (GL_global_1 > 0 && GL_global_2 > 0)
					db.global_stiffness_AA.setValue(GL_global_1 - 1, GL_global_2 - 1, stiffness_spring[i][j]);
				//Caso os graus de liberdade sejam ambos fixos (Matriz Kbb)
				if (GL_global_1 < 0 && GL_global_2 < 0)
					db.global_stiffness_BB.setValue(-GL_global_1 - 1, -GL_global_2 - 1, stiffness_spring[i][j]);
				//Caso os graus de liberdade sejam livre e fixo (Matriz Kab)
				if (GL_global_1 > 0 && GL_global_2 < 0)
					db.global_stiffness_AB.setValue(GL_global_1 - 1, -GL_global_2 - 1, stiffness_spring[i][j]);
				//Caso os graus de liberdade sejam fixo e livre (Matriz Kba)
				if (GL_global_1 < 0 && GL_global_2 > 0)
					db.global_stiffness_BA.setValue(-GL_global_1 - 1, GL_global_2 - 1, stiffness_spring[i][j]);
			}

		}
		/*PrintPtr(stiffness_spring, 6, 6);
		PrintPtr(stiffness_damper, 6, 6);*/
		//PARTE 4 - damping
		for (int i = 0; i < 6; i++)
		{
			//////////////MONTAGEM DO VETOR DE ESFORÇOS DESBALANCEADOS//////////////////
			//Toma do vetor de GL globais, a indexação de cada grau de liberdade global
			if (i < 3)//alpha_A -> i=0,1,2
				GL_global_1 = db.nodes[node_A - 1]->GLs[i + 3];
			else//alpha_B -> i=3,4,5
				GL_global_1 = db.nodes[node_B - 1]->GLs[i];

			//Caso o grau de liberdade seja livre:
			if (GL_global_1 > 0)
			{
				anterior = db.global_P_A(GL_global_1 - 1, 0);
				db.global_P_A(GL_global_1 - 1, 0) = anterior + residual_damper[i];
				anterior = db.global_I_A(GL_global_1 - 1, 0);
				db.global_I_A(GL_global_1 - 1, 0) = anterior + residual_damper[i];
			}
			else
			{
				anterior = db.global_P_B(-GL_global_1 - 1, 0);
				db.global_P_B(-GL_global_1 - 1, 0) = anterior + residual_damper[i];
			}
			for (int j = 0; j < 6; j++)
			{
				//////////////////////MONTAGEM DA MATRIZ DE RIGIDEZ/////////////////////////
				//Toma do vetor de GL globais, a indexação de cada grau de liberdade global
				if (j < 3)//alpha_A -> j=0,1,2
					GL_global_2 = db.nodes[node_A - 1]->GLs[j + 3];
				else//alpha_B -> j=3,4,5
					GL_global_2 = db.nodes[node_B - 1]->GLs[j];

				//Caso os graus de liberdade sejam ambos livres (Matriz Kaa)
				if (GL_global_1 > 0 && GL_global_2 > 0)
					db.global_stiffness_AA.setValue(GL_global_1 - 1, GL_global_2 - 1, stiffness_damper[i][j]);
				//Caso os graus de liberdade sejam ambos fixos (Matriz Kbb)
				if (GL_global_1 < 0 && GL_global_2 < 0)
					db.global_stiffness_BB.setValue(-GL_global_1 - 1, -GL_global_2 - 1, stiffness_damper[i][j]);
				//Caso os graus de liberdade sejam livre e fixo (Matriz Kab)
				if (GL_global_1 > 0 && GL_global_2 < 0)
					db.global_stiffness_AB.setValue(GL_global_1 - 1, -GL_global_2 - 1, stiffness_damper[i][j]);
				//Caso os graus de liberdade sejam fixo e livre (Matriz Kba)
				if (GL_global_1 < 0 && GL_global_2 > 0)
					db.global_stiffness_BA.setValue(-GL_global_1 - 1, GL_global_2 - 1, stiffness_damper[i][j]);
			}

		}
	}
}

void Hinge::ComputeInitialGuessDisplacements()
{
	//Se o vinculo estiver ativo
	if (bool_table.GetAt(db.current_solution_number - 1) == true)
	{
		//Percorre GL e seta condições cinematicas
		for (int i = 0; i < 3; i++)
		{
			if (active_lambda[i] == 1)
			{
				db.nodes[node_B - 1]->displacements[i] = db.nodes[node_A - 1]->displacements[i];
			}
		}
		//Projeção da velocidade angular inicial no sistema local - para impor somente velocidade nos GL restritos (não impõe mesma velocidade angular na direção do eixo da articulação)
		Matrix omega_global_A(3);
		omega_global_A(0, 0) = db.nodes[node_A - 1]->displacements[3];
		omega_global_A(1, 0) = db.nodes[node_A - 1]->displacements[4];
		omega_global_A(2, 0) = db.nodes[node_A - 1]->displacements[5];
		Matrix omega_global_B(3);
		omega_global_B(0, 0) = db.nodes[node_B - 1]->displacements[3];
		omega_global_B(1, 0) = db.nodes[node_B - 1]->displacements[4];
		omega_global_B(2, 0) = db.nodes[node_B - 1]->displacements[5];
		Matrix omega_local_A = (*db.CS[cs - 1]->Q)*omega_global_A;
		Matrix omega_local_B = (*db.CS[cs - 1]->Q)*omega_global_B;
		//Imposições de mesma condição inicial no sistema local - somente direções 1 e 2 locais
		omega_local_B(0, 0) = omega_local_A(0, 0);
		omega_local_B(1, 0) = omega_local_A(1, 0);
		//De volta para o sistema global - nó B
		omega_global_B = transp(*db.CS[cs - 1]->Q)*omega_local_B;
		if (active_lambda[3] == 1 && active_lambda[4] == 1)
		{
			db.nodes[node_B - 1]->displacements[3] = omega_global_B(0, 0);
			db.nodes[node_B - 1]->displacements[4] = omega_global_B(1, 0);
			db.nodes[node_B - 1]->displacements[5] = omega_global_B(2, 0);
		}
	}
}

//Computa efeito das condições iniciais nos nós da restrição
void Hinge::ComputeVelAccel()
{
	//Se o vinculo estiver ativo
	if (bool_table.GetAt(db.current_solution_number - 1) == true)
	{
		//Percorre GL e seta condições cinematicas
		for (int i = 0; i < 3; i++)
		{
			if (active_lambda[i] == 1)
			{
				db.nodes[node_B - 1]->vel[i] = db.nodes[node_A - 1]->vel[i];
				db.nodes[node_B - 1]->accel[i] = db.nodes[node_A - 1]->accel[i];
			}
		}
		//Projeção da velocidade angular inicial no sistema local - para impor somente velocidade nos GL restritos (não impõe mesma velocidade angular na direção do eixo da articulação)
		Matrix omega_global_A(3);
		omega_global_A(0, 0) = db.nodes[node_A - 1]->vel[3];
		omega_global_A(1, 0) = db.nodes[node_A - 1]->vel[4];
		omega_global_A(2, 0) = db.nodes[node_A - 1]->vel[5];
		Matrix omega_global_B(3);
		omega_global_B(0, 0) = db.nodes[node_B - 1]->vel[3];
		omega_global_B(1, 0) = db.nodes[node_B - 1]->vel[4];
		omega_global_B(2, 0) = db.nodes[node_B - 1]->vel[5];
		Matrix omega_local_A = (*db.CS[cs - 1]->Q)*omega_global_A;
		Matrix omega_local_B = (*db.CS[cs - 1]->Q)*omega_global_B;
		//Imposições de mesma condição inicial no sistema local - somente direções 1 e 2 locais
		omega_local_B(0, 0) = omega_local_A(0, 0);
		omega_local_B(1, 0) = omega_local_A(1, 0);
		//De volta para o sistema global - nó B
		omega_global_B = transp(*db.CS[cs - 1]->Q)*omega_local_B;
		if (active_lambda[3] == 1 && active_lambda[4] == 1)
		{
			db.nodes[node_B - 1]->vel[3] = omega_global_B(0, 0);
			db.nodes[node_B - 1]->vel[4] = omega_global_B(1, 0);
			db.nodes[node_B - 1]->vel[5] = omega_global_B(2, 0);
		}
	}
	
}

//Pre-calculo de variaveis que e feito uma unica vez no inicio
void Hinge::PreCalc()
{
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			(*stiffness1)(i, j + 6) = I3(i, j);
			(*stiffness1)(i + 3, j + 6) = -I3(i, j);
			(*stiffness1)(i + 6, j) = I3(i, j);
			(*stiffness1)(i + 6, j + 3) = -I3(i, j);
		}
	}

	//Zerando contribuição dos residuos e rigidez 2
	for (int i = 0; i < 8; i++)
	{
		residual2[i] = 0.0;
		for (int j = 0; j < 8; j++)
		{
			stiffness2[i][j] = 0.0;
		}
	}

	//Zerando contribuição dos residuos e rigidez spring
	for (int i = 0; i < 6; i++)
	{
		residual_spring[i] = 0.0;
		for (int j = 0; j < 6; j++)
		{
			stiffness_spring[i][j] = 0.0;
		}
	}
}

//Salvando variaveis da configuração convergida
void Hinge::SaveLagrange()
{
	for (int i = 0; i < n_GL; i++)
		copy_lambda[i] = lambda[i];
	//PrintPtr(lambda, 5);
	thetai = thetai + thetad;//Incremento da rotação para calculo do momento da mola de torção
}

//Checa quais multiplicadores de lagrange serão ativados,de acordo com a ativação dos GLs dos nós dos quais a special constraint participa
void Hinge::ActivateDOFs()
{
	//Ativa GLs de translação e rotação dos nós A e B
	for (int i = 0; i < 6; i++)
		db.nodes[node_A - 1]->active_GL[i] = 1;
	for (int i = 0; i < 6; i++)
		db.nodes[node_B - 1]->active_GL[i] = 1;
	if (bool_table.GetAt(db.current_solution_number - 1))
	{
		//Ativa GLs de multiplicadores de Lagrange
		active_lambda[0] = 1;
		active_lambda[1] = 1;
		active_lambda[2] = 1;
		active_lambda[3] = 1;
		active_lambda[4] = 1;
	}
	else
	{
		active_lambda[0] = 0;
		active_lambda[1] = 0;
		active_lambda[2] = 0;
		active_lambda[3] = 0;
		active_lambda[4] = 0;
	}
}
//Calcula contribuições do residuo e operador tangente - gerado no AceGen
void Hinge::EvaluateHingeContribution(double *v, double *residual
	, double **stiffness, double *alphaA, double *alphaB, double *ei3A
	, double *ei1B, double *ei2B, double *lambda)
{
	int i01; int i02;
	v[668] = ei2B[2] * lambda[1];
	v[667] = ei1B[2] * lambda[0];
	v[663] = ei2B[1] * lambda[1];
	v[662] = ei1B[1] * lambda[0];
	v[658] = ei2B[0] * lambda[1];
	v[657] = ei1B[0] * lambda[0];
	v[639] = 0.5e0*ei3A[2];
	v[638] = 0.5e0*ei3A[1];
	v[637] = 0.5e0*ei3A[0];
	v[635] = Power(alphaB[2], 2);
	v[634] = 0.5e0*alphaB[2];
	v[633] = 2e0*alphaB[2];
	v[632] = Power(alphaB[1], 2);
	v[631] = 0.5e0*alphaB[1];
	v[630] = 2e0*alphaB[1];
	v[629] = Power(alphaB[0], 2);
	v[628] = 2e0*alphaB[0];
	v[627] = 0.5e0*alphaB[0];
	v[625] = Power(alphaA[2], 2);
	v[624] = 0.5e0*alphaA[2];
	v[623] = 2e0*alphaA[2];
	v[622] = Power(alphaA[1], 2);
	v[621] = 0.5e0*alphaA[1];
	v[620] = 2e0*alphaA[1];
	v[619] = Power(alphaA[0], 2);
	v[618] = 2e0*alphaA[0];
	v[617] = 0.5e0*alphaA[0];
	v[105] = alphaA[1] * v[617];
	v[169] = -v[619] - v[622];
	v[176] = -alphaA[2] + v[105];
	v[172] = alphaA[2] + v[105];
	v[112] = alphaA[2] * v[621];
	v[174] = -alphaA[0] + v[112];
	v[650] = ei3A[0] * v[172] + ei3A[2] * v[174];
	v[171] = alphaA[0] + v[112];
	v[110] = alphaA[2] * v[617];
	v[177] = alphaA[1] + v[110];
	v[647] = ei3A[1] * v[176] + ei3A[2] * v[177];
	v[170] = -alphaA[1] + v[110];
	v[652] = ei3A[0] * v[170] + ei3A[1] * v[171];
	v[183] = 4e0 + v[619] + v[622] + v[625];
	v[221] = 1e0 / Power(v[183], 2);
	v[626] = -4e0*v[221];
	v[224] = v[623] * v[626];
	v[223] = v[620] * v[626];
	v[222] = v[618] * v[626];
	v[175] = -v[622] - v[625];
	v[173] = -v[619] - v[625];
	v[130] = alphaB[1] * v[627];
	v[164] = -v[629] - v[632];
	v[685] = 0.5e0*v[164];
	v[159] = -alphaB[2] + v[130];
	v[156] = alphaB[2] + v[130];
	v[137] = alphaB[2] * v[631];
	v[157] = -alphaB[0] + v[137];
	v[154] = alphaB[0] + v[137];
	v[135] = alphaB[2] * v[627];
	v[160] = alphaB[1] + v[135];
	v[153] = -alphaB[1] + v[135];
	v[197] = 4e0 + v[629] + v[632] + v[635];
	v[225] = 1e0 / Power(v[197], 2);
	v[636] = -4e0*v[225];
	v[228] = v[633] * v[636];
	v[227] = v[630] * v[636];
	v[226] = v[628] * v[636];
	v[163] = -v[629] - v[635];
	v[680] = 0.5e0*v[163];
	v[162] = -v[632] - v[635];
	v[675] = 0.5e0*v[162];
	v[457] = ei3A[0] * v[617];
	v[456] = ei3A[1] * v[621];
	v[461] = v[175] * v[637] + v[647];
	v[460] = v[173] * v[638] + v[650];
	v[459] = v[169] * v[639] + v[652];
	v[454] = ei3A[2] * v[624];
	v[99] = 4e0 / v[183];
	v[649] = 0.5e0*v[99];
	v[540] = -(v[637] * v[99]);
	v[524] = -(v[638] * v[99]);
	v[523] = -(v[639] * v[99]);
	v[124] = 4e0 / v[197];
	v[640] = -0.5e0*v[124];
	v[238] = v[630] * v[640];
	v[236] = v[628] * v[640];
	v[234] = v[633] * v[640];
	v[127] = 1e0 - v[162] * v[640];
	v[366] = ei2B[0] * v[127] + v[124] * (ei2B[1] * v[159] + ei2B[2] * v[160]);
	v[641] = v[366] * v[99];
	v[573] = v[366] * v[540];
	v[512] = ei3A[1] * v[641];
	v[488] = ei3A[2] * v[641];
	v[365] = ei1B[0] * v[127] + v[124] * (ei1B[1] * v[159] + ei1B[2] * v[160]);
	v[642] = v[365] * v[99];
	v[570] = v[365] * v[540];
	v[509] = ei3A[1] * v[642];
	v[485] = ei3A[2] * v[642];
	v[133] = 1e0 - v[163] * v[640];
	v[357] = ei2B[1] * v[133] + v[124] * (ei2B[0] * v[156] + ei2B[2] * v[157]);
	v[643] = v[357] * v[99];
	v[572] = v[357] * v[524];
	v[513] = ei3A[0] * v[643];
	v[514] = v[512] + v[513];
	v[450] = ei3A[2] * v[643];
	v[356] = ei1B[1] * v[133] + v[124] * (ei1B[0] * v[156] + ei1B[2] * v[157]);
	v[644] = v[356] * v[99];
	v[569] = v[356] * v[524];
	v[510] = ei3A[0] * v[644];
	v[511] = v[509] + v[510];
	v[447] = ei3A[2] * v[644];
	v[139] = 1e0 - v[164] * v[640];
	v[350] = ei2B[2] * v[139] + v[124] * (ei2B[0] * v[153] + ei2B[1] * v[154]);
	v[645] = v[350] * v[99];
	v[549] = v[350] * v[523];
	v[489] = ei3A[0] * v[645];
	v[490] = v[488] + v[489];
	v[550] = (v[350] * v[459] + v[357] * v[460] + v[366] * v[461])*v[626];
	v[700] = v[550] + v[572];
	v[451] = ei3A[1] * v[645];
	v[452] = v[450] + v[451];
	v[349] = ei1B[2] * v[139] + v[124] * (ei1B[0] * v[153] + ei1B[1] * v[154]);
	v[646] = v[349] * v[99];
	v[546] = v[349] * v[523];
	v[486] = ei3A[0] * v[646];
	v[487] = v[485] + v[486];
	v[547] = (v[349] * v[459] + v[356] * v[460] + v[365] * v[461])*v[626];
	v[699] = v[547] + v[569];
	v[448] = ei3A[1] * v[646];
	v[449] = v[447] + v[448];
	v[140] = ei3A[0] * (1e0 + v[175] * v[649]) + v[647] * v[99];
	v[659] = 0.5e0*v[140];
	v[648] = v[124] * v[140];
	v[596] = v[140] * v[640];
	v[611] = ei2B[0] * v[596];
	v[608] = ei1B[0] * v[596];
	v[402] = ei2B[2] * v[648];
	v[399] = ei1B[2] * v[648];
	v[374] = v[140] * v[627];
	v[343] = ei2B[1] * v[648];
	v[340] = ei1B[1] * v[648];
	v[141] = ei3A[1] * (1e0 + v[173] * v[649]) + v[650] * v[99];
	v[664] = 0.5e0*v[141];
	v[654] = v[141] * v[157] + v[140] * v[160];
	v[651] = v[124] * v[141];
	v[583] = v[141] * v[640];
	v[612] = ei2B[1] * v[583];
	v[609] = ei1B[1] * v[583];
	v[426] = ei2B[2] * v[651];
	v[423] = ei1B[2] * v[651];
	v[373] = v[141] * v[631];
	v[344] = ei2B[0] * v[651];
	v[345] = v[343] + v[344];
	v[341] = ei1B[0] * v[651];
	v[342] = v[340] + v[341];
	v[142] = ei3A[2] * (1e0 + v[169] * v[649]) + v[652] * v[99];
	v[669] = 0.5e0*v[142];
	v[656] = v[142] * v[153] + v[141] * v[156];
	v[655] = v[142] * v[154] + v[140] * v[159];
	v[653] = v[124] * v[142];
	v[584] = v[142] * v[640];
	v[597] = ei2B[2] * v[584];
	v[593] = ei1B[2] * v[584];
	v[427] = ei2B[1] * v[653];
	v[428] = v[426] + v[427];
	v[424] = ei1B[1] * v[653];
	v[425] = v[423] + v[424];
	v[403] = ei2B[0] * v[653];
	v[404] = v[402] + v[403];
	v[400] = ei1B[0] * v[653];
	v[401] = v[399] + v[400];
	v[378] = v[654] + v[164] * v[669];
	v[377] = v[655] + v[163] * v[664];
	v[376] = v[656] + v[162] * v[659];
	v[598] = (ei2B[0] * v[376] + ei2B[1] * v[377] + ei2B[2] * v[378])*v[636];
	v[704] = v[598] + v[612];
	v[594] = (ei1B[0] * v[376] + ei1B[1] * v[377] + ei1B[2] * v[378])*v[636];
	v[703] = v[594] + v[609];
	v[371] = v[142] * v[634];
	v[146] = v[139] * v[142] + v[124] * v[654];
	v[145] = v[133] * v[141] + v[124] * v[655];
	v[144] = v[127] * v[140] + v[124] * v[656];
	v[149] = v[657] + v[658];
	v[661] = v[141] * v[149];
	v[660] = v[142] * v[149];
	v[589] = -(v[149] * v[659]);
	v[602] = v[228] * v[589];
	v[397] = v[228] * v[660];
	v[394] = v[227] * v[660];
	v[354] = -(v[124] * v[149]);
	v[352] = -(v[354] * v[627]);
	v[338] = v[228] * v[661];
	v[206] = v[124] * v[589];
	v[204] = v[124] * v[661];
	v[199] = v[149] * v[653];
	v[150] = v[662] + v[663];
	v[684] = v[149] * v[153] + v[150] * v[154];
	v[666] = v[140] * v[150];
	v[702] = v[631] * (v[661] + v[666]);
	v[665] = v[142] * v[150];
	v[576] = -(v[150] * v[664]);
	v[603] = v[228] * v[576];
	v[421] = v[228] * v[665];
	v[418] = v[227] * v[665];
	v[363] = v[124] * v[150];
	v[358] = v[363] * v[631];
	v[337] = v[228] * v[666];
	v[339] = v[337] + v[338];
	v[587] = v[226] * v[702];
	v[207] = v[124] * v[576];
	v[205] = v[150] * v[648];
	v[195] = v[150] * v[653];
	v[151] = v[667] + v[668];
	v[679] = v[149] * v[156] + v[151] * v[157];
	v[674] = v[150] * v[159] + v[151] * v[160];
	v[673] = v[124] * v[151];
	v[671] = v[140] * v[151];
	v[670] = v[141] * v[151];
	v[577] = -(v[151] * v[669]);
	v[420] = v[228] * v[670];
	v[417] = v[227] * v[670];
	v[606] = (v[417] + v[418])*v[634];
	v[396] = v[228] * v[671];
	v[393] = v[227] * v[671];
	v[605] = v[226] * v[634] * (v[660] + v[671]);
	v[604] = v[636] * (v[151] * (v[373] + v[374]) + v[149] * (v[141] + v[142] * v[627] - v[633] * v[659]) + v[150] * (-v[140]
		+ v[142] * v[631] - v[633] * v[664]));
	v[364] = -v[363] + v[627] * v[673] + v[228] * v[674] + v[149] * (v[234] + v[228] * v[675]);
	v[672] = v[364] * v[99];
	v[567] = v[364] * v[540];
	v[506] = ei3A[1] * v[672];
	v[482] = ei3A[2] * v[672];
	v[362] = v[363] * v[627] + v[673] + v[227] * v[674] + v[149] * (v[238] + v[227] * v[675]);
	v[676] = v[362] * v[99];
	v[564] = v[362] * v[540];
	v[503] = ei3A[1] * v[676];
	v[479] = ei3A[2] * v[676];
	v[359] = v[634] * v[673];
	v[360] = v[358] + v[359] + v[226] * (v[674] + v[149] * v[675]);
	v[677] = v[360] * v[99];
	v[561] = v[360] * v[540];
	v[500] = ei3A[1] * v[677];
	v[476] = ei3A[2] * v[677];
	v[355] = -v[354] + v[631] * v[673] + v[228] * v[679] + v[150] * (v[234] + v[228] * v[680]);
	v[678] = v[355] * v[99];
	v[566] = v[355] * v[524];
	v[507] = ei3A[0] * v[678];
	v[508] = v[506] + v[507];
	v[444] = ei3A[2] * v[678];
	v[353] = v[352] + v[359] + v[227] * (v[679] + v[150] * v[680]);
	v[681] = v[353] * v[99];
	v[563] = v[353] * v[524];
	v[504] = ei3A[0] * v[681];
	v[505] = v[503] + v[504];
	v[441] = ei3A[2] * v[681];
	v[351] = -(v[354] * v[631]) - v[673] + v[226] * v[679] + v[150] * (v[236] + v[226] * v[680]);
	v[682] = v[351] * v[99];
	v[560] = v[351] * v[524];
	v[501] = ei3A[0] * v[682];
	v[502] = v[500] + v[501];
	v[438] = ei3A[2] * v[682];
	v[348] = v[352] + v[358] + v[228] * (v[684] + v[151] * v[685]);
	v[683] = v[348] * v[99];
	v[543] = v[348] * v[523];
	v[483] = ei3A[0] * v[683];
	v[484] = v[482] + v[483];
	v[544] = (v[348] * v[459] + v[355] * v[460] + v[364] * v[461])*v[626];
	v[698] = v[544] + v[566];
	v[445] = ei3A[1] * v[683];
	v[446] = v[444] + v[445];
	v[347] = v[354] + v[363] * v[634] + v[227] * v[684] + v[151] * (v[238] + v[227] * v[685]);
	v[686] = v[347] * v[99];
	v[539] = v[347] * v[523];
	v[480] = ei3A[0] * v[686];
	v[481] = v[479] + v[480];
	v[541] = (v[347] * v[459] + v[353] * v[460] + v[362] * v[461])*v[626];
	v[697] = v[541] + v[563];
	v[442] = ei3A[1] * v[686];
	v[443] = v[441] + v[442];
	v[346] = v[363] - v[354] * v[634] + v[226] * v[684] + v[151] * (v[236] + v[226] * v[685]);
	v[687] = v[346] * v[99];
	v[536] = v[346] * v[523];
	v[477] = ei3A[0] * v[687];
	v[478] = v[476] + v[477];
	v[537] = (v[346] * v[459] + v[351] * v[460] + v[360] * v[461])*v[626];
	v[696] = v[537] + v[560];
	v[439] = ei3A[1] * v[687];
	v[440] = v[438] + v[439];
	v[201] = v[124] * v[577];
	v[200] = v[124] * v[671];
	v[196] = v[124] * v[670];
	v[152] = v[204] + v[205];
	v[155] = v[139] * v[151] + v[124] * v[684];
	v[689] = ei3A[1] * v[155];
	v[688] = ei3A[0] * v[155];
	v[516] = -(v[155] * v[639]);
	v[474] = v[224] * v[688];
	v[471] = v[223] * v[688];
	v[436] = v[224] * v[689];
	v[433] = v[223] * v[689];
	v[187] = v[516] * v[99];
	v[185] = v[688] * v[99];
	v[181] = v[689] * v[99];
	v[158] = v[133] * v[150] + v[124] * v[679];
	v[694] = ei3A[2] * v[158];
	v[691] = v[158] * v[99];
	v[690] = v[158] * v[224];
	v[517] = -(v[158] * v[638]);
	v[555] = v[224] * v[517];
	v[498] = ei3A[0] * v[690];
	v[435] = ei3A[2] * v[690];
	v[432] = v[223] * v[694];
	v[553] = (v[432] + v[433])*v[624];
	v[192] = v[517] * v[99];
	v[190] = ei3A[0] * v[691];
	v[182] = ei3A[2] * v[691];
	v[161] = v[127] * v[149] + v[124] * v[674];
	v[693] = ei3A[2] * v[161];
	v[692] = ei3A[1] * v[161];
	v[695] = v[621] * (ei3A[0] * v[158] + v[692]);
	v[531] = -(v[161] * v[637]);
	v[556] = v[224] * v[531];
	v[497] = v[224] * v[692];
	v[499] = v[497] + v[498];
	v[532] = v[222] * v[695];
	v[473] = v[224] * v[693];
	v[470] = v[223] * v[693];
	v[558] = v[222] * v[624] * (v[688] + v[693]);
	v[557] = v[626] * (v[155] * (v[456] + v[457]) + v[161] * (-ei3A[1] + ei3A[2] * v[617] - v[623] * v[637]) + v[158] *
		(ei3A[0] + ei3A[2] * v[621] - v[623] * v[638]));
	v[193] = v[531] * v[99];
	v[191] = v[692] * v[99];
	v[186] = v[693] * v[99];
	v[165] = v[149] * v[376] + v[150] * v[377] + v[151] * v[378];
	v[579] = (8e0*v[165]) / Power(v[197], 3);
	v[601] = v[579] * v[633];
	v[705] = v[601] + v[604];
	v[591] = v[228] * v[577] + v[705];
	v[588] = v[227] * v[577] + v[579] * v[630] + v[636] * (v[150] * (v[371] + v[374]) + v[149] * (-v[142] + v[141] * v[627]
		- v[630] * v[659]) + v[151] * (v[140] + v[141] * v[634] - v[630] * v[669]));
	v[202] = v[165] * v[636];
	v[600] = v[202] + v[206] + v[207];
	v[586] = v[201] - v[207] + v[600];
	v[575] = -v[206] + v[207] + v[586];
	v[166] = v[199] + v[200];
	v[167] = v[195] + v[196];
	v[168] = v[181] + v[182];
	v[178] = v[155] * v[459] + v[158] * v[460] + v[161] * v[461];
	v[519] = (8e0*v[178]) / Power(v[183], 3);
	v[554] = v[519] * v[623];
	v[701] = v[554] + v[557];
	v[534] = v[224] * v[516] + v[701];
	v[530] = v[223] * v[516] + v[519] * v[620] + v[626] * (v[158] * (v[454] + v[457]) + v[161] * (ei3A[2] + ei3A[1] * v[617]
		- v[620] * v[637]) + v[155] * (-ei3A[0] + ei3A[1] * v[624] - v[620] * v[639]));
	v[188] = v[178] * v[626];
	v[552] = v[188] + v[192] + v[193];
	v[529] = v[187] - v[192] + v[552];
	v[515] = v[192] - v[193] + v[529];
	v[179] = v[185] + v[186];
	v[180] = v[190] + v[191];
	residual[0] = v[181] - v[182] + v[515] * v[618] + v[180] * v[621] + v[179] * v[624];
	residual[1] = -v[185] + v[186] + v[180] * v[617] + v[529] * v[620] + v[168] * v[624];
	residual[2] = v[190] - v[191] + v[179] * v[617] + v[168] * v[621] + v[552] * v[623];
	residual[3] = v[195] - v[196] + v[575] * v[628] + v[152] * v[631] + v[166] * v[634];
	residual[4] = -v[199] + v[200] + v[152] * v[627] + v[586] * v[630] + v[167] * v[634];
	residual[5] = v[204] - v[205] + v[166] * v[627] + v[167] * v[631] + v[600] * v[633];
	residual[6] = ei1B[0] * v[144] + ei1B[1] * v[145] + ei1B[2] * v[146];
	residual[7] = ei2B[0] * v[144] + ei2B[1] * v[145] + ei2B[2] * v[146];
	stiffness[0][0] = 2e0*v[515] + v[532] + v[558] + v[618] * (v[222] * (v[516] + v[517]) + v[519] * v[618] + v[626] *
		(v[161] * (v[454] + v[456]) + v[158] * (-ei3A[2] + ei3A[0] * v[621] - v[618] * v[638]) + v[155] * (ei3A[1]
		+ ei3A[0] * v[624] - v[618] * v[639]))) + v[222] * (v[689] - v[694]);
	stiffness[0][1] = 0.5e0*v[180] - v[432] + v[433] + (v[223] * v[517] + v[530])*v[618] + (v[470] + v[471])*v[624]
		+ v[223] * v[695];
	stiffness[0][2] = 0.5e0*v[179] - v[435] + v[436] + (v[534] + v[555])*v[618] + v[499] * v[621] + (v[473] + v[474]
		)*v[624];
	stiffness[0][3] = -v[438] + v[439] + v[502] * v[621] + v[478] * v[624] + v[618] * (v[536] + v[696]);
	stiffness[0][4] = -v[441] + v[442] + v[505] * v[621] + v[481] * v[624] + v[618] * (v[539] + v[697]);
	stiffness[0][5] = -v[444] + v[445] + v[508] * v[621] + v[484] * v[624] + v[618] * (v[543] + v[698]);
	stiffness[0][6] = -v[447] + v[448] + v[511] * v[621] + v[487] * v[624] + v[618] * (v[546] + v[699]);
	stiffness[0][7] = -v[450] + v[451] + v[514] * v[621] + v[490] * v[624] + v[618] * (v[549] + v[700]);
	stiffness[1][1] = v[470] - v[471] + 2e0*v[529] + v[532] + v[553] + (v[530] + v[223] * v[531])*v[620];
	stiffness[1][2] = 0.5e0*v[168] + v[473] - v[474] + v[499] * v[617] + (v[534] + v[556])*v[620] + (v[435] + v[436]
		)*v[624];
	stiffness[1][3] = v[476] - v[477] + v[502] * v[617] + (v[536] + v[537] + v[561])*v[620] + v[440] * v[624];
	stiffness[1][4] = v[479] - v[480] + v[505] * v[617] + (v[539] + v[541] + v[564])*v[620] + v[443] * v[624];
	stiffness[1][5] = v[482] - v[483] + v[508] * v[617] + (v[543] + v[544] + v[567])*v[620] + v[446] * v[624];
	stiffness[1][6] = v[485] - v[486] + v[511] * v[617] + (v[546] + v[547] + v[570])*v[620] + v[449] * v[624];
	stiffness[1][7] = v[488] - v[489] + v[514] * v[617] + (v[549] + v[550] + v[573])*v[620] + v[452] * v[624];
	stiffness[2][2] = -v[497] + v[498] + 2e0*v[552] + v[553] + v[558] + v[623] * (v[555] + v[556] + v[701]);
	stiffness[2][3] = -v[500] + v[501] + v[478] * v[617] + v[440] * v[621] + v[623] * (v[561] + v[696]);
	stiffness[2][4] = -v[503] + v[504] + v[481] * v[617] + v[443] * v[621] + v[623] * (v[564] + v[697]);
	stiffness[2][5] = -v[506] + v[507] + v[484] * v[617] + v[446] * v[621] + v[623] * (v[567] + v[698]);
	stiffness[2][6] = -v[509] + v[510] + v[487] * v[617] + v[449] * v[621] + v[623] * (v[570] + v[699]);
	stiffness[2][7] = -v[512] + v[513] + v[490] * v[617] + v[452] * v[621] + v[623] * (v[573] + v[700]);
	stiffness[3][3] = 2e0*v[575] + v[587] + v[605] + v[628] * (v[226] * (v[576] + v[577]) + v[579] * v[628] + v[636] *
		(v[149] * (v[371] + v[373]) + v[150] * (v[142] + v[140] * v[631] - v[628] * v[664]) + v[151] * (-v[141] + v[140] * v[634]
		- v[628] * v[669]))) + v[226] * (v[665] - v[670]);
	stiffness[3][4] = 0.5e0*v[152] - v[417] + v[418] + (v[227] * v[576] + v[588])*v[628] + (v[393] + v[394])*v[634]
		+ v[227] * v[702];
	stiffness[3][5] = 0.5e0*v[166] - v[420] + v[421] + (v[591] + v[603])*v[628] + v[339] * v[631] + (v[396] + v[397]
		)*v[634];
	stiffness[3][6] = -v[423] + v[424] + v[342] * v[631] + v[401] * v[634] + v[628] * (v[593] + v[703]);
	stiffness[3][7] = -v[426] + v[427] + v[345] * v[631] + v[404] * v[634] + v[628] * (v[597] + v[704]);
	stiffness[4][4] = v[393] - v[394] + 2e0*v[586] + v[587] + v[606] + (v[588] + v[227] * v[589])*v[630];
	stiffness[4][5] = 0.5e0*v[167] + v[396] - v[397] + v[339] * v[627] + (v[591] + v[602])*v[630] + (v[420] + v[421]
		)*v[634];
	stiffness[4][6] = v[399] - v[400] + v[342] * v[627] + (v[593] + v[594] + v[608])*v[630] + v[425] * v[634];
	stiffness[4][7] = v[402] - v[403] + v[345] * v[627] + (v[597] + v[598] + v[611])*v[630] + v[428] * v[634];
	stiffness[5][5] = -v[337] + v[338] + 2e0*v[600] + v[605] + v[606] + v[633] * (v[602] + v[603] + v[705]);
	stiffness[5][6] = -v[340] + v[341] + v[401] * v[627] + v[425] * v[631] + v[633] * (v[608] + v[703]);
	stiffness[5][7] = -v[343] + v[344] + v[404] * v[627] + v[428] * v[631] + v[633] * (v[611] + v[704]);
	stiffness[6][6] = 0e0;
	stiffness[6][7] = 0e0;
	stiffness[7][7] = 0e0;
	for (i01 = 1; i01<8; i01++){
		for (i02 = 0; i02<i01; i02++){
			stiffness[i01][i02] = stiffness[i02][i01];
		}
	};
};

//Calcula contribuições do residuo e operador tangente - gerado no AceGen (sem usar SMSD)
void Hinge::EvaluateHingeContribution2(double *v, double *residual
	, double **stiffness, double *alphaA, double *alphaB, double *ei3A
	, double *ei1B, double *ei2B, double *lambda)
{
	int i01; int i02;
	v[360] = Power(alphaB[2], 2);
	v[359] = 0.5e0*alphaB[0] * alphaB[2];
	v[358] = 0.5e0*alphaB[1];
	v[357] = Power(alphaB[1], 2);
	v[366] = v[357] + v[360];
	v[356] = alphaB[0] * v[358];
	v[355] = Power(alphaB[0], 2);
	v[354] = Power(alphaA[2], 2);
	v[353] = 0.5e0*alphaA[0] * alphaA[2];
	v[352] = 0.5e0*alphaA[1];
	v[351] = Power(alphaA[1], 2);
	v[361] = v[351] + v[354];
	v[350] = alphaA[0] * v[352];
	v[349] = Power(alphaA[0], 2);
	v[121] = alphaA[2] * v[352];
	v[140] = alphaB[2] * v[358];
	v[108] = 4e0 / (4e0 + v[349] + v[361]);
	v[365] = ei3A[0] * v[108];
	v[364] = ei3A[1] * v[108];
	v[363] = -0.5e0*v[108];
	v[362] = ei3A[2] * v[108];
	v[150] = (-alphaA[1] - v[353])*v[362] - ei3A[0] * (1e0 + v[361] * v[363]) + (alphaA[2] - v[350])*v[364];
	v[167] = -(v[108] * v[150]);
	v[149] = (-alphaA[0] + v[121])*v[362] + ei3A[1] * (1e0 + (v[349] + v[354])*v[363]) + (alphaA[2] + v[350])*v[365];
	v[165] = -(v[108] * v[149]);
	v[147] = -(ei3A[2] * (1e0 + (v[349] + v[351])*v[363])) + (-alphaA[0] - v[121])*v[364] + (alphaA[1] - v[353]
		)*v[365];
	v[153] = -(v[108] * v[147]);
	v[124] = alphaA[2] * v[363];
	v[156] = -(v[124] * v[147]);
	v[125] = -(alphaA[1] * v[363]);
	v[169] = -(v[125] * v[149]);
	v[126] = alphaA[0] * v[363];
	v[170] = -(v[126] * v[150]);
	v[127] = 4e0 / (4e0 + v[355] + v[366]);
	v[367] = -0.5e0*v[127];
	v[318] = v[127] * v[150];
	v[316] = v[127] * v[149];
	v[305] = v[127] * v[147];
	v[130] = 1e0 + v[366] * v[367];
	v[131] = v[127] * (-alphaB[2] + v[356]);
	v[132] = v[127] * (alphaB[1] + v[359]);
	v[161] = ei2B[0] * v[130] + ei2B[1] * v[131] + ei2B[2] * v[132];
	v[257] = -(v[108] * v[161]);
	v[254] = -(v[126] * v[161]);
	v[183] = v[127] * v[161];
	v[155] = ei1B[0] * v[130] + ei1B[1] * v[131] + ei1B[2] * v[132];
	v[221] = -(v[108] * v[155]);
	v[218] = -(v[126] * v[155]);
	v[179] = v[127] * v[155];
	v[134] = v[127] * (alphaB[2] + v[356]);
	v[136] = 1e0 + (v[355] + v[360])*v[367];
	v[137] = v[127] * (-alphaB[0] + v[140]);
	v[162] = ei2B[0] * v[134] + ei2B[1] * v[136] + ei2B[2] * v[137];
	v[259] = v[108] * v[162];
	v[255] = v[125] * v[162];
	v[186] = -(v[127] * v[162]);
	v[158] = ei1B[0] * v[134] + ei1B[1] * v[136] + ei1B[2] * v[137];
	v[223] = v[108] * v[158];
	v[219] = v[125] * v[158];
	v[182] = -(v[127] * v[158]);
	v[139] = v[127] * (-alphaB[1] + v[359]);
	v[141] = v[127] * (alphaB[0] + v[140]);
	v[142] = 1e0 + (v[355] + v[357])*v[367];
	v[163] = ei2B[0] * v[139] + ei2B[1] * v[141] + ei2B[2] * v[142];
	v[244] = -(v[124] * v[163]);
	v[241] = -(v[108] * v[163]);
	v[177] = v[127] * v[163];
	v[160] = ei1B[0] * v[139] + ei1B[1] * v[141] + ei1B[2] * v[142];
	v[208] = -(v[124] * v[160]);
	v[205] = -(v[108] * v[160]);
	v[175] = v[127] * v[160];
	v[143] = alphaB[2] * v[367];
	v[309] = v[143] * v[147];
	v[176] = v[143] * v[163];
	v[174] = v[143] * v[160];
	v[144] = -(alphaB[1] * v[367]);
	v[321] = v[144] * v[149];
	v[184] = -(v[144] * v[162]);
	v[325] = -(v[149] * (v[144] * v[161] + v[177])) + v[150] * (v[176] + v[184]) + v[147] * (-(v[143] * v[161]) + v[186]);
	v[180] = -(v[144] * v[158]);
	v[295] = -(v[149] * (v[144] * v[155] + v[175])) + v[150] * (v[174] + v[180]) + v[147] * (-(v[143] * v[155]) + v[182]);
	v[145] = alphaB[0] * v[367];
	v[322] = v[145] * v[150];
	v[185] = v[145] * v[161];
	v[338] = -(v[149] * (v[144] * v[163] - v[183])) + v[147] * (v[184] + v[185]) - v[150] * (v[145] * v[163] + v[186]);
	v[332] = -(v[150] * (v[145] * v[162] + v[177])) + v[147] * (-(v[143] * v[162]) + v[183]) - v[149] * (v[176] + v[185]);
	v[181] = v[145] * v[155];
	v[320] = -(v[149] * (v[144] * v[160] - v[179])) + v[147] * (v[180] + v[181]) - v[150] * (v[145] * v[160] + v[182]);
	v[308] = -(v[150] * (v[145] * v[158] + v[175])) + v[147] * (-(v[143] * v[158]) + v[179]) - v[149] * (v[174] + v[181]);
	v[146] = v[156] + v[169];
	v[148] = -(v[125] * v[150]) + v[153];
	v[274] = -(v[148] * v[161]) + v[146] * v[162];
	v[271] = -(v[148] * v[155]) + v[146] * v[158];
	v[151] = v[124] * v[150] + v[165];
	v[273] = -(v[151] * v[162]) + v[148] * v[163];
	v[272] = v[151] * v[161] - v[146] * v[163];
	v[270] = -(v[151] * v[158]) + v[148] * v[160];
	v[269] = v[151] * v[155] - v[146] * v[160];
	v[228] = -(v[146] * v[161]) - v[148] * v[162] - v[151] * v[163];
	v[192] = -(v[146] * v[155]) - v[148] * v[158] - v[151] * v[160];
	v[154] = -(v[126] * v[149]) - v[153];
	v[157] = v[156] + v[170];
	v[283] = -(v[157] * v[161]) + v[154] * v[162];
	v[280] = -(v[155] * v[157]) + v[154] * v[158];
	v[159] = -(v[124] * v[149]) + v[167];
	v[282] = -(v[159] * v[162]) + v[157] * v[163];
	v[281] = v[159] * v[161] - v[154] * v[163];
	v[279] = -(v[158] * v[159]) + v[157] * v[160];
	v[278] = v[155] * v[159] - v[154] * v[160];
	v[243] = -(v[154] * v[161]) - v[157] * v[162] - v[159] * v[163];
	v[207] = -(v[154] * v[155]) - v[157] * v[158] - v[159] * v[160];
	v[166] = v[126] * v[147] - v[165];
	v[168] = -(v[125] * v[147]) - v[167];
	v[292] = v[162] * v[166] - v[161] * v[168];
	v[289] = v[158] * v[166] - v[155] * v[168];
	v[171] = v[169] + v[170];
	v[291] = v[163] * v[168] - v[162] * v[171];
	v[290] = -(v[163] * v[166]) + v[161] * v[171];
	v[288] = v[160] * v[168] - v[158] * v[171];
	v[287] = -(v[160] * v[166]) + v[155] * v[171];
	v[253] = -(v[161] * v[166]) - v[162] * v[168] - v[163] * v[171];
	v[217] = -(v[155] * v[166]) - v[158] * v[168] - v[160] * v[171];
	v[190] = v[124] * v[155] + v[223];
	v[191] = -(v[125] * v[155]) + v[205];
	v[195] = -(v[149] * v[190]) - v[147] * v[191];
	v[193] = v[208] + v[219];
	v[197] = -(v[150] * v[190]) + v[147] * v[193];
	v[196] = v[150] * v[191] + v[149] * v[193];
	v[198] = -(v[149] * v[155]);
	v[199] = -(v[150] * v[158]);
	v[201] = -(v[147] * v[155]);
	v[202] = v[150] * v[160];
	v[204] = v[124] * v[158] + v[221];
	v[206] = v[126] * v[158] - v[205];
	v[211] = -(v[150] * v[204]) + v[147] * v[206];
	v[209] = v[208] + v[218];
	v[213] = -(v[149] * v[204]) - v[147] * v[209];
	v[212] = v[149] * v[206] + v[150] * v[209];
	v[214] = v[147] * v[158];
	v[215] = v[149] * v[160];
	v[220] = v[218] + v[219];
	v[222] = -(v[125] * v[160]) - v[221];
	v[224] = v[126] * v[160] - v[223];
	v[226] = v[124] * v[161] + v[259];
	v[227] = -(v[125] * v[161]) + v[241];
	v[231] = -(v[149] * v[226]) - v[147] * v[227];
	v[229] = v[244] + v[255];
	v[233] = -(v[150] * v[226]) + v[147] * v[229];
	v[232] = v[150] * v[227] + v[149] * v[229];
	v[234] = -(v[149] * v[161]);
	v[235] = -(v[150] * v[162]);
	v[237] = -(v[147] * v[161]);
	v[238] = v[150] * v[163];
	v[240] = v[124] * v[162] + v[257];
	v[242] = v[126] * v[162] - v[241];
	v[247] = -(v[150] * v[240]) + v[147] * v[242];
	v[245] = v[244] + v[254];
	v[249] = -(v[149] * v[240]) - v[147] * v[245];
	v[248] = v[149] * v[242] + v[150] * v[245];
	v[250] = v[147] * v[162];
	v[251] = v[149] * v[163];
	v[256] = v[254] + v[255];
	v[258] = -(v[125] * v[163]) - v[257];
	v[260] = v[126] * v[163] - v[259];
	v[296] = v[309] + v[321];
	v[297] = -(v[143] * v[150]) + v[316];
	v[327] = -(v[163] * v[296]) + v[161] * v[297];
	v[300] = -(v[160] * v[296]) + v[155] * v[297];
	v[298] = v[144] * v[150] + v[305];
	v[329] = -(v[162] * v[297]) + v[163] * v[298];
	v[328] = v[162] * v[296] - v[161] * v[298];
	v[302] = -(v[158] * v[297]) + v[160] * v[298];
	v[301] = v[158] * v[296] - v[155] * v[298];
	v[306] = v[145] * v[149] - v[305];
	v[307] = v[143] * v[149] + v[318];
	v[334] = -(v[163] * v[306]) + v[161] * v[307];
	v[312] = -(v[160] * v[306]) + v[155] * v[307];
	v[310] = v[309] + v[322];
	v[336] = -(v[162] * v[307]) + v[163] * v[310];
	v[335] = v[162] * v[306] - v[161] * v[310];
	v[314] = -(v[158] * v[307]) + v[160] * v[310];
	v[313] = v[158] * v[306] - v[155] * v[310];
	v[317] = -(v[145] * v[147]) - v[316];
	v[319] = v[144] * v[147] - v[318];
	v[323] = v[321] + v[322];
	residual[0] = lambda[0] * v[192] + lambda[1] * v[228];
	residual[1] = lambda[0] * v[207] + lambda[1] * v[243];
	residual[2] = lambda[0] * v[217] + lambda[1] * v[253];
	residual[3] = lambda[0] * v[295] + lambda[1] * v[325];
	residual[4] = lambda[0] * v[308] + lambda[1] * v[332];
	residual[5] = lambda[0] * v[320] + lambda[1] * v[338];
	residual[6] = -(v[150] * v[155]) + v[149] * v[158] - v[147] * v[160];
	residual[7] = -(v[150] * v[161]) + v[149] * v[162] - v[147] * v[163];
	stiffness[0][0] = lambda[0] * (v[126] * v[192] + v[108] * v[195] - v[125] * v[196] - v[124] * v[197]) + lambda[1] *
		(v[126] * v[228] + v[108] * v[231] - v[125] * v[232] - v[124] * v[233]);
	stiffness[0][1] = lambda[0] * (v[124] * v[195] - v[126] * v[196] + v[108] * (v[197] - 0.5e0*(alphaA[1] * v[192]
		+ v[198] + v[199]))) + lambda[1] * (v[124] * v[231] - v[126] * v[232] + v[108] * (v[233] - 0.5e0*(alphaA[1] * v[228]
		+ v[234] + v[235])));
	stiffness[0][2] = lambda[0] * (v[125] * v[195] + v[126] * v[197] + v[108] * (v[196] - 0.5e0*(alphaA[2] * v[192]
		- v[201] - v[202]))) + lambda[1] * (v[125] * v[231] + v[126] * v[233] + v[108] * (v[232] - 0.5e0*(alphaA[2] * v[228]
		- v[237] - v[238])));
	stiffness[0][3] = lambda[0] * (-(v[143] * v[269]) + v[127] * v[270] - v[144] * v[271]) + lambda[1] * (-(v[143] * v[272]
		) + v[127] * v[273] - v[144] * v[274]);
	stiffness[0][4] = lambda[0] * (v[127] * v[269] + v[143] * v[270] - v[145] * v[271]) + lambda[1] * (v[127] * v[272]
		+ v[143] * v[273] - v[145] * v[274]);
	stiffness[0][5] = lambda[0] * (v[145] * v[269] + v[144] * v[270] + v[127] * v[271]) + lambda[1] * (v[145] * v[272]
		+ v[144] * v[273] + v[127] * v[274]);
	stiffness[0][6] = v[192];
	stiffness[0][7] = v[228];
	stiffness[1][1] = lambda[0] * (-(v[125] * v[207]) + v[108] * v[211] - v[126] * v[212] + v[124] * v[213]) + lambda[1] * (-
		(v[125] * v[243]) + v[108] * v[247] - v[126] * v[248] + v[124] * v[249]);
	stiffness[1][2] = lambda[0] * (v[126] * v[211] + v[125] * v[213] + v[108] * (v[212] - 0.5e0*(alphaA[2] * v[207]
		+ v[214] + v[215]))) + lambda[1] * (v[126] * v[247] + v[125] * v[249] + v[108] * (v[248] - 0.5e0*(alphaA[2] * v[243]
		+ v[250] + v[251])));
	stiffness[1][3] = lambda[0] * (-(v[143] * v[278]) + v[127] * v[279] - v[144] * v[280]) + lambda[1] * (-(v[143] * v[281]
		) + v[127] * v[282] - v[144] * v[283]);
	stiffness[1][4] = lambda[0] * (v[127] * v[278] + v[143] * v[279] - v[145] * v[280]) + lambda[1] * (v[127] * v[281]
		+ v[143] * v[282] - v[145] * v[283]);
	stiffness[1][5] = lambda[0] * (v[145] * v[278] + v[144] * v[279] + v[127] * v[280]) + lambda[1] * (v[145] * v[281]
		+ v[144] * v[282] + v[127] * v[283]);
	stiffness[1][6] = v[207];
	stiffness[1][7] = v[243];
	stiffness[2][2] = lambda[0] * (v[124] * v[217] + v[125] * (-(v[149] * v[220]) - v[147] * v[222]) + v[126] * (-
		(v[150] * v[220]) + v[147] * v[224]) + v[108] * (v[150] * v[222] + v[149] * v[224])) + lambda[1] * (v[124] * v[253]
		+ v[125] * (-(v[149] * v[256]) - v[147] * v[258]) + v[126] * (-(v[150] * v[256]) + v[147] * v[260]) + v[108] *
		(v[150] * v[258] + v[149] * v[260]));
	stiffness[2][3] = lambda[0] * (-(v[143] * v[287]) + v[127] * v[288] - v[144] * v[289]) + lambda[1] * (-(v[143] * v[290]
		) + v[127] * v[291] - v[144] * v[292]);
	stiffness[2][4] = lambda[0] * (v[127] * v[287] + v[143] * v[288] - v[145] * v[289]) + lambda[1] * (v[127] * v[290]
		+ v[143] * v[291] - v[145] * v[292]);
	stiffness[2][5] = lambda[0] * (v[145] * v[287] + v[144] * v[288] + v[127] * v[289]) + lambda[1] * (v[145] * v[290]
		+ v[144] * v[291] + v[127] * v[292]);
	stiffness[2][6] = v[217];
	stiffness[2][7] = v[253];
	stiffness[3][3] = lambda[0] * (v[145] * v[295] - v[143] * v[300] - v[144] * v[301] + v[127] * v[302]) + lambda[1] *
		(v[145] * v[325] - v[143] * v[327] - v[144] * v[328] + v[127] * v[329]);
	stiffness[3][4] = lambda[0] * (v[127] * (-0.5e0*(-v[198] - v[199] + alphaB[1] * v[295]) + v[300]) - v[145] * v[301]
		+ v[143] * v[302]) + lambda[1] * (v[127] * (-0.5e0*(-v[234] - v[235] + alphaB[1] * v[325]) + v[327]) - v[145] * v[328]
		+ v[143] * v[329]);
	stiffness[3][5] = lambda[0] * (v[145] * v[300] + v[127] * (-0.5e0*(v[201] + v[202] + alphaB[2] * v[295]) + v[301])
		+ v[144] * v[302]) + lambda[1] * (v[145] * v[327] + v[127] * (-0.5e0*(v[237] + v[238] + alphaB[2] * v[325]) + v[328])
		+ v[144] * v[329]);
	stiffness[3][6] = v[295];
	stiffness[3][7] = v[325];
	stiffness[4][4] = lambda[0] * (-(v[144] * v[308]) + v[127] * v[312] - v[145] * v[313] + v[143] * v[314]) + lambda[1] * (-
		(v[144] * v[332]) + v[127] * v[334] - v[145] * v[335] + v[143] * v[336]);
	stiffness[4][5] = lambda[0] * (v[145] * v[312] + v[127] * (-0.5e0*(-v[214] - v[215] + alphaB[2] * v[308]) + v[313])
		+ v[144] * v[314]) + lambda[1] * (v[145] * v[334] + v[127] * (-0.5e0*(-v[250] - v[251] + alphaB[2] * v[332]) + v[335])
		+ v[144] * v[336]);
	stiffness[4][6] = v[308];
	stiffness[4][7] = v[332];
	stiffness[5][5] = lambda[0] * (v[127] * (v[158] * v[317] - v[155] * v[319]) + v[143] * v[320] + v[145] * (-
		(v[160] * v[317]) + v[155] * v[323]) + v[144] * (v[160] * v[319] - v[158] * v[323])) + lambda[1] * (v[127] *
		(v[162] * v[317] - v[161] * v[319]) + v[145] * (-(v[163] * v[317]) + v[161] * v[323]) + v[144] * (v[163] * v[319]
		- v[162] * v[323]) + v[143] * v[338]);
	stiffness[5][6] = v[320];
	stiffness[5][7] = v[338];
	stiffness[6][6] = 0e0;
	stiffness[6][7] = 0e0;
	stiffness[7][7] = 0e0;
	for (i01 = 1; i01<8; i01++){
		for (i02 = 0; i02<i01; i02++){
			stiffness[i01][i02] = stiffness[i02][i01];
		}
	};
};

//Calcula contribuições da mola de torção (spring)
void Hinge::EvaluateTorsionSpring(double *v, double *residual
	, double **stiffness, double *alphaA, double *alphaB, double *ei1A, double *ei1B, double *ei3A
	, double(*stiffc), double(*thetai)
	, double(*thetad))
{
	int i01; int i02;
	v[994] = 1e0*(*stiffc);
	v[969] = Power(ei1B[2], 2);
	v[966] = Power(ei1B[1], 2);
	v[965] = Power(ei1B[0], 2);
	v[959] = ei1B[0] * ei1B[2];
	v[955] = ei1B[1] * ei1B[2];
	v[954] = ei1B[0] * ei1B[1];
	v[952] = Power(alphaB[2], 2);
	v[951] = 0.5e0*alphaB[2];
	v[950] = 2e0*alphaB[2];
	v[949] = Power(alphaB[1], 2);
	v[948] = 0.5e0*alphaB[1];
	v[947] = 2e0*alphaB[1];
	v[946] = Power(alphaB[0], 2);
	v[945] = 2e0*alphaB[0];
	v[944] = 0.5e0*alphaB[0];
	v[942] = Power(alphaA[2], 2);
	v[941] = 0.5e0*alphaA[2];
	v[940] = 2e0*alphaA[2];
	v[939] = Power(alphaA[1], 2);
	v[938] = 0.5e0*alphaA[1];
	v[937] = 2e0*alphaA[1];
	v[936] = Power(alphaA[0], 2);
	v[935] = 2e0*alphaA[0];
	v[934] = 0.5e0*alphaA[0];
	v[76] = alphaA[1] * v[934];
	v[206] = -v[936] - v[939];
	v[1017] = 0.5e0*v[206];
	v[212] = alphaA[2] + v[76];
	v[210] = -alphaA[2] + v[76];
	v[83] = alphaA[2] * v[938];
	v[213] = -alphaA[0] + v[83];
	v[207] = alphaA[0] + v[83];
	v[81] = alphaA[2] * v[934];
	v[209] = alphaA[1] + v[81];
	v[208] = -alphaA[1] + v[81];
	v[218] = 4e0 + v[936] + v[939] + v[942];
	v[256] = 1e0 / Power(v[218], 2);
	v[943] = -4e0*v[256];
	v[259] = v[940] * v[943];
	v[960] = 0.5e0*v[259];
	v[328] = v[206] * v[960];
	v[258] = v[937] * v[943];
	v[956] = 0.5e0*v[258];
	v[257] = v[935] * v[943];
	v[958] = 0.5e0*v[257];
	v[214] = -v[936] - v[942];
	v[1019] = 0.5e0*v[214];
	v[296] = v[214] * v[956];
	v[211] = -v[939] - v[942];
	v[1018] = 0.5e0*v[211];
	v[260] = v[211] * v[958];
	v[101] = alphaB[1] * v[944];
	v[193] = -v[946] - v[949];
	v[1021] = 0.5e0*v[193];
	v[201] = alphaB[2] + v[101];
	v[197] = -alphaB[2] + v[101];
	v[108] = alphaB[2] * v[948];
	v[199] = -alphaB[0] + v[108];
	v[194] = alphaB[0] + v[108];
	v[106] = alphaB[2] * v[944];
	v[196] = alphaB[1] + v[106];
	v[195] = -alphaB[1] + v[106];
	v[232] = 4e0 + v[946] + v[949] + v[952];
	v[344] = 1e0 / Power(v[232], 2);
	v[953] = -4e0*v[344];
	v[347] = v[950] * v[953];
	v[964] = 0.5e0*v[347];
	v[392] = v[193] * v[964];
	v[999] = ei1B[2] * v[392];
	v[346] = v[947] * v[953];
	v[961] = 0.5e0*v[346];
	v[345] = v[945] * v[953];
	v[963] = 0.5e0*v[345];
	v[200] = -v[946] - v[952];
	v[1026] = 0.5e0*v[200];
	v[363] = v[200] * v[961];
	v[198] = -v[949] - v[952];
	v[1025] = 0.5e0*v[198];
	v[348] = v[198] * v[963];
	v[1015] = -(ei1B[0] * v[348]);
	v[451] = v[363] * v[954];
	v[400] = v[296] * v[954];
	v[454] = v[363] * v[955];
	v[421] = v[296] * v[955];
	v[70] = 4e0 / v[218];
	v[957] = -0.5e0*v[70];
	v[326] = v[937] * v[957];
	v[327] = v[326] + v[206] * v[956];
	v[324] = v[935] * v[957];
	v[325] = v[324] + v[206] * v[958];
	v[315] = v[207] * v[257] + v[70];
	v[307] = v[208] * v[258] - v[70];
	v[300] = v[70] * v[941];
	v[316] = v[207] * v[258] + v[300];
	v[336] = ei3A[0] * v[307] + ei3A[1] * v[316] + ei3A[2] * v[327];
	v[306] = v[208] * v[257] + v[300];
	v[335] = ei3A[0] * v[306] + ei3A[1] * v[315] + ei3A[2] * v[325];
	v[301] = v[213] * v[258] + v[300];
	v[412] = v[301] * v[955];
	v[403] = v[301] * v[959];
	v[299] = v[213] * v[257] - v[70];
	v[411] = v[299] * v[955];
	v[402] = v[299] * v[959];
	v[297] = v[940] * v[957];
	v[298] = v[297] + v[214] * v[960];
	v[422] = v[298] * v[955];
	v[401] = v[298] * v[954];
	v[295] = v[324] + v[214] * v[958];
	v[420] = v[295] * v[955];
	v[399] = v[295] * v[954];
	v[294] = v[212] * v[259] + v[70];
	v[419] = v[294] * v[959];
	v[410] = v[294] * v[954];
	v[291] = v[70] * v[938];
	v[317] = v[207] * v[259] + v[291];
	v[302] = v[213] * v[259] + v[291];
	v[987] = v[419] + v[422] + v[302] * v[969];
	v[413] = v[302] * v[955];
	v[986] = v[410] + v[413] + v[298] * v[966];
	v[404] = v[302] * v[959];
	v[985] = v[401] + v[404] + v[294] * v[965];
	v[305] = ei3A[0] * v[294] + ei3A[1] * v[298] + ei3A[2] * v[302];
	v[292] = v[212] * v[257] + v[291];
	v[991] = v[399] + v[402] + v[292] * v[965];
	v[417] = v[292] * v[959];
	v[993] = v[417] + v[420] + v[299] * v[969];
	v[408] = v[292] * v[954];
	v[992] = v[408] + v[411] + v[295] * v[966];
	v[303] = ei3A[0] * v[292] + ei3A[1] * v[295] + ei3A[2] * v[299];
	v[280] = v[70] * v[934];
	v[308] = v[208] * v[259] + v[280];
	v[337] = ei3A[0] * v[308] + ei3A[1] * v[317] + ei3A[2] * v[328];
	v[293] = v[212] * v[258] + v[280];
	v[988] = v[400] + v[403] + v[293] * v[965];
	v[418] = v[293] * v[959];
	v[990] = v[418] + v[421] + v[301] * v[969];
	v[409] = v[293] * v[954];
	v[989] = v[409] + v[412] + v[296] * v[966];
	v[304] = ei3A[0] * v[293] + ei3A[1] * v[296] + ei3A[2] * v[301];
	v[281] = v[209] * v[259] + v[280];
	v[279] = v[209] * v[258] + v[70];
	v[278] = v[209] * v[257] + v[300];
	v[271] = v[210] * v[259] - v[70];
	v[270] = v[210] * v[258] + v[280];
	v[269] = v[210] * v[257] + v[291];
	v[288] = ei3A[0] * v[260] + ei3A[1] * v[269] + ei3A[2] * v[278];
	v[262] = v[297] + v[211] * v[960];
	v[290] = ei3A[0] * v[262] + ei3A[1] * v[271] + ei3A[2] * v[281];
	v[261] = v[326] + v[211] * v[956];
	v[289] = ei3A[0] * v[261] + ei3A[1] * v[270] + ei3A[2] * v[279];
	v[73] = 1e0 - v[211] * v[957];
	v[138] = v[73] * v[954];
	v[131] = v[73] * v[959];
	v[74] = v[210] * v[70];
	v[141] = v[74] * v[954];
	v[133] = v[74] * v[955];
	v[75] = v[209] * v[70];
	v[982] = v[131] + v[133] + v[75] * v[969];
	v[142] = v[75] * v[959];
	v[978] = v[141] + v[142] + v[73] * v[965];
	v[139] = v[75] * v[955];
	v[980] = v[138] + v[139] + v[74] * v[966];
	v[134] = ei3A[0] * v[73] + ei3A[1] * v[74] + ei3A[2] * v[75];
	v[77] = v[212] * v[70];
	v[79] = 1e0 - v[214] * v[957];
	v[80] = v[213] * v[70];
	v[129] = ei3A[0] * v[77] + ei3A[1] * v[79] + ei3A[2] * v[80];
	v[82] = v[208] * v[70];
	v[149] = -(v[82] * v[954]);
	v[146] = -(v[82] * v[959]);
	v[84] = v[207] * v[70];
	v[153] = -(v[84] * v[954]);
	v[147] = -(v[84] * v[955]);
	v[85] = 1e0 - v[206] * v[957];
	v[983] = v[146] + v[147] - v[85] * v[969];
	v[154] = -(v[85] * v[959]);
	v[979] = v[153] + v[154] - v[82] * v[965];
	v[150] = -(v[85] * v[955]);
	v[981] = v[149] + v[150] - v[84] * v[966];
	v[144] = ei3A[0] * v[82] + ei3A[1] * v[84] + ei3A[2] * v[85];
	v[518] = -(v[144] * v[73]) + v[134] * v[82];
	v[176] = -(v[144] * v[74]) + v[134] * v[84];
	v[174] = -(v[144] * v[75]) + v[134] * v[85];
	v[95] = 4e0 / v[232];
	v[962] = -0.5e0*v[95];
	v[390] = v[947] * v[962];
	v[391] = v[390] + v[193] * v[961];
	v[1000] = ei1B[2] * v[391];
	v[388] = v[945] * v[962];
	v[389] = v[388] + v[193] * v[963];
	v[1001] = ei1B[2] * v[389];
	v[379] = v[194] * v[345] + v[95];
	v[1007] = ei1B[1] * v[379];
	v[371] = v[195] * v[346] - v[95];
	v[1004] = ei1B[0] * v[371];
	v[367] = v[95] * v[951];
	v[380] = v[194] * v[346] + v[367];
	v[1005] = ei1B[1] * v[380];
	v[370] = v[195] * v[345] + v[367];
	v[1006] = ei1B[0] * v[370];
	v[376] = v[144] * v[348] - v[134] * v[370];
	v[368] = v[199] * v[346] + v[367];
	v[466] = v[368] * v[955];
	v[460] = v[368] * v[959];
	v[366] = v[199] * v[345] - v[95];
	v[465] = v[366] * v[955];
	v[459] = v[366] * v[959];
	v[364] = v[950] * v[962];
	v[365] = v[364] + v[200] * v[964];
	v[455] = v[365] * v[955];
	v[452] = v[365] * v[954];
	v[362] = v[388] + v[200] * v[963];
	v[453] = v[362] * v[955];
	v[450] = v[362] * v[954];
	v[361] = v[201] * v[347] + v[95];
	v[449] = v[361] * v[959];
	v[446] = v[361] * v[954];
	v[358] = v[948] * v[95];
	v[381] = v[194] * v[347] + v[358];
	v[1003] = ei1B[1] * v[381];
	v[369] = v[199] * v[347] + v[358];
	v[975] = -v[449] - v[455] - v[369] * v[969];
	v[467] = v[369] * v[955];
	v[461] = v[369] * v[959];
	v[359] = v[201] * v[345] + v[358];
	v[447] = v[359] * v[959];
	v[977] = -v[447] - v[453] - v[366] * v[969];
	v[444] = v[359] * v[954];
	v[356] = v[944] * v[95];
	v[372] = v[195] * v[347] + v[356];
	v[1002] = ei1B[0] * v[372];
	v[360] = v[201] * v[346] + v[356];
	v[448] = v[360] * v[959];
	v[976] = -v[448] - v[454] - v[368] * v[969];
	v[445] = v[360] * v[954];
	v[357] = v[196] * v[347] + v[356];
	v[1008] = ei1B[2] * v[357];
	v[398] = v[144] * v[357] - v[134] * v[392];
	v[355] = v[196] * v[346] + v[95];
	v[1009] = ei1B[2] * v[355];
	v[397] = v[144] * v[355] - v[134] * v[391];
	v[354] = v[196] * v[345] + v[367];
	v[1010] = ei1B[2] * v[354];
	v[396] = v[144] * v[354] - v[134] * v[389];
	v[353] = v[197] * v[347] - v[95];
	v[1012] = ei1B[1] * v[353];
	v[387] = v[144] * v[353] - v[134] * v[381];
	v[352] = v[197] * v[346] + v[356];
	v[1014] = ei1B[1] * v[352];
	v[386] = v[144] * v[352] - v[134] * v[380];
	v[351] = v[197] * v[345] + v[358];
	v[1016] = ei1B[1] * v[351];
	v[385] = v[144] * v[351] - v[134] * v[379];
	v[350] = v[364] + v[198] * v[964];
	v[1011] = -(ei1B[0] * v[350]);
	v[378] = v[144] * v[350] - v[134] * v[372];
	v[349] = v[390] + v[198] * v[961];
	v[1013] = -(ei1B[0] * v[349]);
	v[377] = v[144] * v[349] - v[134] * v[371];
	v[98] = 1e0 - v[198] * v[962];
	v[972] = -(ei1B[0] * v[98]);
	v[99] = v[197] * v[95];
	v[973] = ei1B[1] * v[99];
	v[100] = v[196] * v[95];
	v[968] = ei1B[2] * v[100];
	v[102] = v[201] * v[95];
	v[104] = 1e0 - v[200] * v[962];
	v[105] = v[199] * v[95];
	v[107] = v[195] * v[95];
	v[970] = ei1B[0] * v[107];
	v[375] = -(v[107] * v[290]) + v[337] * v[98];
	v[374] = -(v[107] * v[289]) + v[336] * v[98];
	v[373] = -(v[107] * v[288]) + v[335] * v[98];
	v[188] = -(v[107] * v[134]) + v[144] * v[98];
	v[109] = v[194] * v[95];
	v[971] = ei1B[1] * v[109];
	v[384] = -(v[109] * v[290]) + v[337] * v[99];
	v[383] = -(v[109] * v[289]) + v[336] * v[99];
	v[382] = -(v[109] * v[288]) + v[335] * v[99];
	v[185] = -(v[109] * v[134]) + v[144] * v[99];
	v[110] = 1e0 - v[193] * v[962];
	v[967] = ei1B[2] * v[110];
	v[395] = -(v[110] * v[290]) + v[100] * v[337];
	v[394] = -(v[110] * v[289]) + v[100] * v[336];
	v[393] = -(v[110] * v[288]) + v[100] * v[335];
	v[181] = -(v[110] * v[134]) + v[100] * v[144];
	v[653] = v[518] * v[965];
	v[619] = v[107] * v[965] + ei1B[0] * (v[967] + v[971]);
	v[595] = -(ei1B[0] * (v[968] + v[973])) - v[965] * v[98];
	v[464] = v[452] + v[461] + v[361] * v[965];
	v[463] = v[451] + v[460] + v[360] * v[965];
	v[462] = v[450] + v[459] + v[359] * v[965];
	v[112] = v[79] * v[954];
	v[113] = v[80] * v[959];
	v[155] = -v[112] - v[113] - v[77] * v[965];
	v[572] = -(v[144] * v[155]) + v[129] * v[979];
	v[548] = v[134] * v[155] + v[129] * v[978];
	v[114] = v[77] * v[954];
	v[640] = v[176] * v[966];
	v[611] = v[109] * v[966] + ei1B[1] * (v[967] + v[970]);
	v[587] = ei1B[1] * (-v[968] + v[972]) - v[966] * v[99];
	v[470] = v[446] + v[467] + v[365] * v[966];
	v[469] = v[445] + v[466] + v[363] * v[966];
	v[468] = v[444] + v[465] + v[362] * v[966];
	v[116] = v[80] * v[955];
	v[151] = -v[114] - v[116] - v[79] * v[966];
	v[564] = -(v[144] * v[151]) + v[129] * v[981];
	v[540] = v[134] * v[151] + v[129] * v[980];
	v[117] = v[77] * v[959];
	v[118] = v[79] * v[955];
	v[974] = v[117] + v[118] + v[80] * v[969];
	v[627] = v[174] * v[969];
	v[603] = v[110] * v[969] + ei1B[2] * (v[970] + v[971]);
	v[579] = -(v[100] * v[969]) + ei1B[2] * (v[972] - v[973]);
	v[482] = -(v[155] * v[350]) - v[151] * v[353] - v[464] * v[73] - v[470] * v[74] + v[357] * v[974] + v[75] * v[975];
	v[481] = -(v[155] * v[349]) - v[151] * v[352] - v[463] * v[73] - v[469] * v[74] + v[355] * v[974] + v[75] * v[976];
	v[480] = -(v[155] * v[348]) - v[151] * v[351] - v[462] * v[73] - v[468] * v[74] + v[354] * v[974] + v[75] * v[977];
	v[476] = v[155] * v[372] + v[151] * v[381] + v[464] * v[82] + v[470] * v[84] - v[392] * v[974] - v[85] * v[975];
	v[475] = v[155] * v[371] + v[151] * v[380] + v[463] * v[82] + v[469] * v[84] - v[391] * v[974] - v[85] * v[976];
	v[474] = v[155] * v[370] + v[151] * v[379] + v[462] * v[82] + v[468] * v[84] - v[389] * v[974] - v[85] * v[977];
	v[440] = v[372] * v[978] + v[350] * v[979] + v[381] * v[980] + v[353] * v[981] + v[392] * v[982] + v[357] * v[983];
	v[488] = v[129] * v[440] + v[134] * v[476] + v[144] * v[482];
	v[437] = v[371] * v[978] + v[349] * v[979] + v[380] * v[980] + v[352] * v[981] + v[391] * v[982] + v[355] * v[983];
	v[487] = v[129] * v[437] + v[134] * v[475] + v[144] * v[481];
	v[434] = v[370] * v[978] + v[348] * v[979] + v[379] * v[980] + v[351] * v[981] + v[389] * v[982] + v[354] * v[983];
	v[486] = v[129] * v[434] + v[134] * v[474] + v[144] * v[480];
	v[431] = (v[110] * v[281] - v[100] * v[328])*v[969] + v[965] * (v[107] * v[262] - v[308] * v[98]) + v[954] *
		(v[109] * v[262] + v[107] * v[271] - v[317] * v[98] - v[308] * v[99]) + v[966] * (v[109] * v[271] - v[317] * v[99])
		+ ei1B[2] * (ei1B[0] * (v[110] * v[262] + v[107] * v[281] - v[100] * v[308] - v[328] * v[98]) + ei1B[1] * (v[110] * v[271]
		+ v[109] * v[281] - v[100] * v[317] - v[328] * v[99]));
	v[428] = (v[110] * v[279] - v[100] * v[327])*v[969] + v[965] * (v[107] * v[261] - v[307] * v[98]) + v[954] *
		(v[109] * v[261] + v[107] * v[270] - v[316] * v[98] - v[307] * v[99]) + v[966] * (v[109] * v[270] - v[316] * v[99])
		+ ei1B[2] * (ei1B[0] * (v[110] * v[261] + v[107] * v[279] - v[100] * v[307] - v[327] * v[98]) + ei1B[1] * (v[110] * v[270]
		+ v[109] * v[279] - v[100] * v[316] - v[327] * v[99]));
	v[425] = (v[110] * v[278] - v[100] * v[325])*v[969] + v[965] * (v[107] * v[260] - v[306] * v[98]) + v[954] *
		(v[109] * v[260] + v[107] * v[269] - v[315] * v[98] - v[306] * v[99]) + v[966] * (v[109] * v[269] - v[315] * v[99])
		+ ei1B[2] * (ei1B[0] * (v[110] * v[260] + v[107] * v[278] - v[100] * v[306] - v[325] * v[98]) + ei1B[1] * (v[110] * v[269]
		+ v[109] * v[278] - v[100] * v[315] - v[325] * v[99]));
	v[178] = v[107] * v[978] + v[979] * v[98] + v[109] * v[980] + v[110] * v[982] + v[100] * v[983] + v[981] * v[99];
	v[692] = ei3A[1] * v[178] + v[185] * v[966];
	v[679] = ei3A[2] * v[178] + v[181] * v[969];
	v[666] = ei3A[0] * v[178] + v[188] * v[965];
	v[556] = v[144] * v[974] + v[129] * v[983];
	v[532] = -(v[134] * v[974]) + v[129] * v[982];
	v[120] = v[102] * v[954];
	v[121] = v[102] * v[959];
	v[122] = v[104] * v[954];
	v[123] = v[104] * v[955];
	v[161] = v[121] + v[123] + v[105] * v[969];
	v[124] = v[105] * v[959];
	v[166] = v[122] + v[124] + v[102] * v[965];
	v[125] = v[105] * v[955];
	v[984] = -v[120] - v[125] - v[104] * v[966];
	v[479] = -(v[166] * v[262]) - v[161] * v[281] + v[271] * v[984] + v[98] * v[985] + v[100] * v[987] + v[986] * v[99];
	v[478] = -(v[166] * v[261]) - v[161] * v[279] + v[270] * v[984] + v[98] * v[988] + v[989] * v[99] + v[100] * v[990];
	v[477] = -(v[166] * v[260]) - v[161] * v[278] + v[269] * v[984] + v[98] * v[991] + v[99] * v[992] + v[100] * v[993];
	v[473] = v[166] * v[308] + v[161] * v[328] - v[317] * v[984] - v[107] * v[985] - v[109] * v[986] - v[110] * v[987];
	v[472] = v[166] * v[307] + v[161] * v[327] - v[316] * v[984] - v[107] * v[988] - v[109] * v[989] - v[110] * v[990];
	v[471] = v[166] * v[306] + v[161] * v[325] - v[315] * v[984] - v[107] * v[991] - v[109] * v[992] - v[110] * v[993];
	v[162] = v[109] * v[151] + v[107] * v[155] + v[166] * v[82] + v[161] * v[85] - v[110] * v[974] - v[84] * v[984];
	v[620] = ei3A[0] * v[162] - v[144] * v[166] + v[129] * v[619];
	v[612] = ei3A[1] * v[162] + v[129] * v[611] + v[144] * v[984];
	v[604] = -(v[144] * v[161]) + ei3A[2] * v[162] + v[129] * v[603];
	v[157] = -(v[166] * v[73]) - v[161] * v[75] + v[100] * v[974] - v[155] * v[98] + v[74] * v[984] - v[151] * v[99];
	v[596] = ei3A[0] * v[157] + v[134] * v[166] + v[129] * v[595];
	v[588] = ei3A[1] * v[157] + v[129] * v[587] - v[134] * v[984];
	v[580] = ei3A[2] * v[157] + v[134] * v[161] + v[129] * v[579];
	v[485] = v[162] * v[290] + v[178] * v[305] + v[157] * v[337] + v[129] * v[431] + v[134] * v[473] + v[144] * v[479];
	v[484] = v[162] * v[289] + v[178] * v[304] + v[157] * v[336] + v[129] * v[428] + v[134] * v[472] + v[144] * v[478];
	v[483] = v[162] * v[288] + v[178] * v[303] + v[157] * v[335] + v[129] * v[425] + v[134] * v[471] + v[144] * v[477];
	v[135] = v[144] * v[157] + v[134] * v[162] + v[129] * v[178];
	v[497] = 1e0 - (v[135] * v[135]);
	v[489] = 1e0 / sqrt(v[497]);
	v[500] = v[489] * v[994];
	v[126] = asin(v[135]);
	v[496] = ((*thetai) + v[126])*v[994];
	v[499] = (v[135] * v[496]) / Power(v[497], 0.15e1);
	v[996] = v[499] + v[489] * v[500];
	v[505] = v[488] * v[996];
	v[995] = ei1B[0] * v[505];
	v[578] = v[505] * v[572];
	v[926] = v[578] * v[962];
	v[570] = v[505] * v[564];
	v[921] = v[570] * v[95];
	v[562] = v[505] * v[556];
	v[915] = v[562] * v[95];
	v[554] = v[505] * v[548];
	v[914] = v[554] * v[95];
	v[546] = v[505] * v[540];
	v[902] = v[546] * v[95];
	v[538] = v[505] * v[532];
	v[530] = v[505] * v[955];
	v[524] = v[518] * v[995];
	v[652] = ei1B[1] * v[524] + v[174] * v[530] + v[505] * v[640];
	v[927] = v[652] * v[962];
	v[639] = ei1B[2] * v[524] + v[176] * v[530] + v[505] * v[627];
	v[903] = v[639] * v[95];
	v[517] = ei1B[2] * v[995];
	v[511] = ei1B[1] * v[995];
	v[665] = v[176] * v[511] + v[174] * v[517] + v[505] * v[653];
	v[922] = v[665] * v[95];
	v[504] = v[487] * v[996];
	v[997] = ei1B[0] * v[504];
	v[577] = v[504] * v[572];
	v[569] = v[504] * v[564];
	v[561] = v[504] * v[556];
	v[909] = v[561] * v[95];
	v[553] = v[504] * v[548];
	v[908] = v[553] * v[95];
	v[545] = v[504] * v[540];
	v[896] = v[545] * v[95];
	v[537] = v[504] * v[532];
	v[529] = v[504] * v[955];
	v[523] = v[518] * v[997];
	v[650] = ei1B[1] * v[523] + v[174] * v[529] + v[504] * v[640];
	v[637] = ei1B[2] * v[523] + v[176] * v[529] + v[504] * v[627];
	v[897] = v[637] * v[95];
	v[516] = ei1B[2] * v[997];
	v[510] = ei1B[1] * v[997];
	v[663] = v[176] * v[510] + v[174] * v[516] + v[504] * v[653];
	v[503] = v[486] * v[996];
	v[998] = ei1B[0] * v[503];
	v[568] = v[503] * v[564];
	v[560] = v[503] * v[556];
	v[552] = v[503] * v[548];
	v[544] = v[503] * v[540];
	v[536] = v[503] * v[532];
	v[528] = v[503] * v[955];
	v[522] = v[518] * v[998];
	v[648] = ei1B[1] * v[522] + v[174] * v[528] + v[503] * v[640];
	v[635] = ei1B[2] * v[522] + v[176] * v[528] + v[503] * v[627];
	v[515] = ei1B[2] * v[998];
	v[509] = ei1B[1] * v[998];
	v[661] = v[176] * v[509] + v[174] * v[515] + v[503] * v[653];
	v[502] = v[485] * v[996];
	v[527] = v[502] * v[955];
	v[514] = v[502] * v[959];
	v[508] = v[502] * v[954];
	v[501] = v[484] * v[996];
	v[526] = v[501] * v[955];
	v[513] = v[501] * v[959];
	v[507] = v[501] * v[954];
	v[498] = v[483] * v[996];
	v[525] = v[498] * v[955];
	v[512] = v[498] * v[959];
	v[506] = v[498] * v[954];
	v[137] = v[489] * v[496];
	v[626] = v[505] * v[620] + v[137] * (-(v[144] * v[464]) + ei3A[0] * v[476] + v[129] * (v[372] * v[965] + ei1B[0] * (v[1003]
		+ v[999])));
	v[885] = v[626] * v[957];
	v[625] = v[504] * v[620] + v[137] * (-(v[144] * v[463]) + ei3A[0] * v[475] + v[129] * (ei1B[0] * (v[1000] + v[1005])
		+ v[371] * v[965]));
	v[882] = v[625] * v[957];
	v[624] = v[503] * v[620] + v[137] * (-(v[144] * v[462]) + ei3A[0] * v[474] + v[129] * (ei1B[0] * (v[1001] + v[1007])
		+ v[370] * v[965]));
	v[879] = v[624] * v[957];
	v[623] = v[137] * (-(v[166] * v[337]) + ei3A[0] * v[473] + v[305] * v[619]) + v[502] * v[620];
	v[875] = v[623] * v[957];
	v[622] = v[137] * (-(v[166] * v[336]) + ei3A[0] * v[472] + v[304] * v[619]) + v[501] * v[620];
	v[618] = v[505] * v[612] + v[137] * (-(v[144] * v[470]) + ei3A[1] * v[476] + v[129] * (v[381] * v[966] + ei1B[1] * (v[1002]
		+ v[999])));
	v[817] = v[618] * v[70];
	v[617] = v[504] * v[612] + v[137] * (-(v[144] * v[469]) + ei3A[1] * v[475] + v[129] * (ei1B[1] * (v[1000] + v[1004])
		+ v[380] * v[966]));
	v[814] = v[617] * v[70];
	v[616] = v[503] * v[612] + v[137] * (-(v[144] * v[468]) + ei3A[1] * v[474] + v[129] * (ei1B[1] * (v[1001] + v[1006])
		+ v[379] * v[966]));
	v[811] = v[616] * v[70];
	v[615] = v[502] * v[612] + v[137] * (ei3A[1] * v[473] + v[305] * v[611] + v[337] * v[984]);
	v[870] = v[615] * v[70];
	v[614] = v[501] * v[612] + v[137] * (ei3A[1] * v[472] + v[304] * v[611] + v[336] * v[984]);
	v[613] = v[498] * v[612] + v[137] * (ei3A[1] * v[471] + v[303] * v[611] + v[335] * v[984]);
	v[610] = v[505] * v[604] + v[137] * (ei3A[2] * v[476] + v[129] * (ei1B[2] * (v[1002] + v[1003]) + v[392] * v[969])
		+ v[144] * v[975]);
	v[798] = v[610] * v[70];
	v[609] = v[504] * v[604] + v[137] * (ei3A[2] * v[475] + v[129] * (ei1B[2] * (v[1004] + v[1005]) + v[391] * v[969])
		+ v[144] * v[976]);
	v[795] = v[609] * v[70];
	v[608] = v[503] * v[604] + v[137] * (ei3A[2] * v[474] + v[129] * (ei1B[2] * (v[1006] + v[1007]) + v[389] * v[969])
		+ v[144] * v[977]);
	v[792] = v[608] * v[70];
	v[607] = v[137] * (-(v[161] * v[337]) + ei3A[2] * v[473] + v[305] * v[603]) + v[502] * v[604];
	v[855] = v[607] * v[70];
	v[606] = v[137] * (-(v[161] * v[336]) + ei3A[2] * v[472] + v[304] * v[603]) + v[501] * v[604];
	v[849] = v[606] * v[70];
	v[605] = v[137] * (-(v[161] * v[335]) + ei3A[2] * v[471] + v[303] * v[603]) + v[498] * v[604];
	v[602] = v[505] * v[596] + v[137] * (v[134] * v[464] + ei3A[0] * v[482] - v[129] * (ei1B[0] * (v[1008] + v[1012])
		+ v[350] * v[965]));
	v[799] = v[602] * v[70];
	v[800] = v[798] + v[799];
	v[601] = v[504] * v[596] + v[137] * (v[134] * v[463] + ei3A[0] * v[481] - v[129] * (ei1B[0] * (v[1009] + v[1014])
		+ v[349] * v[965]));
	v[796] = v[601] * v[70];
	v[797] = v[795] + v[796];
	v[600] = v[503] * v[596] + v[137] * (v[134] * v[462] + ei3A[0] * v[480] - v[129] * (ei1B[0] * (v[1010] + v[1016])
		+ v[348] * v[965]));
	v[793] = v[600] * v[70];
	v[794] = v[792] + v[793];
	v[599] = v[137] * (v[166] * v[290] + ei3A[0] * v[479] + v[305] * v[595]) + v[502] * v[596];
	v[854] = v[599] * v[70];
	v[598] = v[137] * (v[166] * v[289] + ei3A[0] * v[478] + v[304] * v[595]) + v[501] * v[596];
	v[848] = v[598] * v[70];
	v[597] = v[137] * (v[166] * v[288] + ei3A[0] * v[477] + v[303] * v[595]) + v[498] * v[596];
	v[594] = v[505] * v[588] + v[137] * (v[134] * v[470] + ei3A[1] * v[482] + v[129] * (ei1B[1] * (-v[1008] + v[1011])
		- v[353] * v[966]));
	v[781] = v[594] * v[70];
	v[593] = v[504] * v[588] + v[137] * (v[134] * v[469] + ei3A[1] * v[481] + v[129] * (ei1B[1] * (-v[1009] + v[1013])
		- v[352] * v[966]));
	v[778] = v[593] * v[70];
	v[592] = v[503] * v[588] + v[137] * (v[134] * v[468] + ei3A[1] * v[480] + v[129] * (ei1B[1] * (-v[1010] + v[1015])
		- v[351] * v[966]));
	v[775] = v[592] * v[70];
	v[591] = v[502] * v[588] + v[137] * (ei3A[1] * v[479] + v[305] * v[587] - v[290] * v[984]);
	v[839] = v[591] * v[70];
	v[590] = v[501] * v[588] + v[137] * (ei3A[1] * v[478] + v[304] * v[587] - v[289] * v[984]);
	v[833] = v[590] * v[70];
	v[589] = v[498] * v[588] + v[137] * (ei3A[1] * v[477] + v[303] * v[587] - v[288] * v[984]);
	v[586] = v[505] * v[580] + v[137] * (ei3A[2] * v[482] + v[129] * (ei1B[2] * (v[1011] - v[1012]) - v[357] * v[969])
		- v[134] * v[975]);
	v[864] = v[586] * v[957];
	v[585] = v[504] * v[580] + v[137] * (ei3A[2] * v[481] + v[129] * (ei1B[2] * (v[1013] - v[1014]) - v[355] * v[969])
		- v[134] * v[976]);
	v[861] = v[585] * v[957];
	v[584] = v[503] * v[580] + v[137] * (ei3A[2] * v[480] + v[129] * (ei1B[2] * (v[1015] - v[1016]) - v[354] * v[969])
		- v[134] * v[977]);
	v[858] = v[584] * v[957];
	v[583] = v[137] * (v[161] * v[290] + ei3A[2] * v[479] + v[305] * v[579]) + v[502] * v[580];
	v[582] = v[137] * (v[161] * v[289] + ei3A[2] * v[478] + v[304] * v[579]) + v[501] * v[580];
	v[581] = v[137] * (v[161] * v[288] + ei3A[2] * v[477] + v[303] * v[579]) + v[498] * v[580];
	v[184] = v[137] * v[954];
	v[180] = v[137] * v[959];
	v[678] = v[184] * v[387] + v[180] * v[398] + v[185] * v[511] + v[181] * v[517] + v[505] * v[666] + v[137] * (ei3A[0] * v[440]
		+ v[378] * v[965]);
	v[816] = v[678] * v[70];
	v[818] = v[816] + v[817];
	v[676] = v[184] * v[386] + v[180] * v[397] + v[185] * v[510] + v[181] * v[516] + v[504] * v[666] + v[137] * (ei3A[0] * v[437]
		+ v[377] * v[965]);
	v[813] = v[676] * v[70];
	v[815] = v[813] + v[814];
	v[674] = v[184] * v[385] + v[180] * v[396] + v[185] * v[509] + v[181] * v[515] + v[503] * v[666] + v[137] * (ei3A[0] * v[434]
		+ v[376] * v[965]);
	v[810] = v[674] * v[70];
	v[812] = v[810] + v[811];
	v[672] = v[184] * v[384] + v[180] * v[395] + v[185] * v[508] + v[181] * v[514] + v[502] * v[666] + v[137] * (ei3A[0] * v[431]
		+ v[375] * v[965]);
	v[871] = v[672] * v[70];
	v[670] = v[184] * v[383] + v[180] * v[394] + v[185] * v[507] + v[181] * v[513] + v[501] * v[666] + v[137] * (ei3A[0] * v[428]
		+ v[374] * v[965]);
	v[668] = v[184] * v[382] + v[180] * v[393] + v[185] * v[506] + v[181] * v[512] + v[498] * v[666] + v[137] * (ei3A[0] * v[425]
		+ v[373] * v[965]);
	v[172] = ei1B[0] * v[137] * v[518];
	v[170] = v[137] * v[955];
	v[704] = v[184] * v[378] + v[170] * v[398] + v[188] * v[511] + v[181] * v[530] + v[505] * v[692] + v[137] * (ei3A[1] * v[440]
		+ v[387] * v[966]);
	v[886] = v[704] * v[957];
	v[702] = v[184] * v[377] + v[170] * v[397] + v[188] * v[510] + v[181] * v[529] + v[504] * v[692] + v[137] * (ei3A[1] * v[437]
		+ v[386] * v[966]);
	v[883] = v[702] * v[957];
	v[700] = v[184] * v[376] + v[170] * v[396] + v[188] * v[509] + v[181] * v[528] + v[503] * v[692] + v[137] * (ei3A[1] * v[434]
		+ v[385] * v[966]);
	v[880] = v[700] * v[957];
	v[698] = v[184] * v[375] + v[170] * v[395] + v[188] * v[508] + v[181] * v[527] + v[502] * v[692] + v[137] * (ei3A[1] * v[431]
		+ v[384] * v[966]);
	v[876] = v[698] * v[957];
	v[696] = v[184] * v[374] + v[170] * v[394] + v[188] * v[507] + v[181] * v[526] + v[501] * v[692] + v[137] * (ei3A[1] * v[428]
		+ v[383] * v[966]);
	v[694] = v[184] * v[373] + v[170] * v[393] + v[188] * v[506] + v[181] * v[525] + v[498] * v[692] + v[137] * (ei3A[1] * v[425]
		+ v[382] * v[966]);
	v[691] = v[180] * v[378] + v[170] * v[387] + v[188] * v[517] + v[185] * v[530] + v[505] * v[679] + v[137] * (ei3A[2] * v[440]
		+ v[398] * v[969]);
	v[865] = (v[1017] * v[586] + v[207] * v[594] + v[208] * v[602] + v[209] * v[610] + v[210] * v[618] + v[1018] * v[626]
		+ v[212] * v[678] + v[213] * v[691] + v[1019] * v[704])*v[943];
	v[1043] = v[865] + v[886];
	v[780] = v[691] * v[70];
	v[782] = v[780] + v[781];
	v[689] = v[180] * v[377] + v[170] * v[386] + v[188] * v[516] + v[185] * v[529] + v[504] * v[679] + v[137] * (ei3A[2] * v[437]
		+ v[397] * v[969]);
	v[862] = (v[1017] * v[585] + v[207] * v[593] + v[208] * v[601] + v[209] * v[609] + v[210] * v[617] + v[1018] * v[625]
		+ v[212] * v[676] + v[213] * v[689] + v[1019] * v[702])*v[943];
	v[1042] = v[862] + v[883];
	v[777] = v[689] * v[70];
	v[779] = v[777] + v[778];
	v[687] = v[180] * v[376] + v[170] * v[385] + v[188] * v[515] + v[185] * v[528] + v[503] * v[679] + v[137] * (ei3A[2] * v[434]
		+ v[396] * v[969]);
	v[859] = (v[1017] * v[584] + v[207] * v[592] + v[208] * v[600] + v[209] * v[608] + v[210] * v[616] + v[1018] * v[624]
		+ v[212] * v[674] + v[213] * v[687] + v[1019] * v[700])*v[943];
	v[1041] = v[859] + v[880];
	v[774] = v[687] * v[70];
	v[776] = v[774] + v[775];
	v[685] = v[180] * v[375] + v[170] * v[384] + v[188] * v[514] + v[185] * v[527] + v[502] * v[679] + v[137] * (ei3A[2] * v[431]
		+ v[395] * v[969]);
	v[840] = v[685] * v[70];
	v[683] = v[180] * v[374] + v[170] * v[383] + v[188] * v[513] + v[185] * v[526] + v[501] * v[679] + v[137] * (ei3A[2] * v[428]
		+ v[394] * v[969]);
	v[834] = v[683] * v[70];
	v[681] = v[180] * v[373] + v[170] * v[382] + v[188] * v[512] + v[185] * v[525] + v[498] * v[679] + v[137] * (ei3A[2] * v[425]
		+ v[393] * v[969]);
	v[136] = v[137] * v[532];
	v[1024] = -0.5e0*v[136];
	v[236] = v[136] * v[962];
	v[140] = v[137] * v[540];
	v[900] = v[140] * v[347];
	v[1052] = v[900] + v[902];
	v[894] = v[140] * v[346];
	v[1050] = v[894] + v[896];
	v[230] = v[140] * v[95];
	v[143] = v[137] * v[548];
	v[912] = v[143] * v[347];
	v[906] = v[143] * v[346];
	v[234] = v[143] * v[95];
	v[148] = v[137] * v[556];
	v[1047] = v[143] + v[148];
	v[913] = v[148] * v[347];
	v[1055] = v[913] + v[915];
	v[907] = v[148] * v[346];
	v[1051] = v[907] + v[909];
	v[740] = v[1055] + v[912] + v[914];
	v[235] = v[148] * v[95];
	v[152] = v[137] * v[564];
	v[919] = v[152] * v[347];
	v[239] = v[152] * v[95];
	v[156] = v[137] * v[572];
	v[1022] = 0.5e0*v[156];
	v[924] = -(v[156] * v[964]);
	v[1056] = v[924] + v[926];
	v[241] = v[156] * v[962];
	v[158] = v[137] * v[580];
	v[1030] = -0.5e0*v[158];
	v[222] = v[158] * v[957];
	v[159] = v[137] * v[588];
	v[837] = v[159] * v[259];
	v[1038] = v[837] + v[839];
	v[831] = v[159] * v[258];
	v[1036] = v[831] + v[833];
	v[216] = v[159] * v[70];
	v[160] = v[137] * v[596];
	v[852] = v[160] * v[259];
	v[846] = v[160] * v[258];
	v[220] = v[160] * v[70];
	v[163] = v[137] * v[604];
	v[1033] = v[160] + v[163];
	v[853] = v[163] * v[259];
	v[1044] = v[853] + v[855];
	v[847] = v[163] * v[258];
	v[1037] = v[847] + v[849];
	v[791] = v[1044] + v[852] + v[854];
	v[221] = v[163] * v[70];
	v[165] = v[137] * v[612];
	v[868] = v[165] * v[259];
	v[225] = v[165] * v[70];
	v[167] = v[137] * v[620];
	v[1029] = 0.5e0*v[167];
	v[873] = -(v[167] * v[960]);
	v[1045] = v[873] + v[875];
	v[227] = v[167] * v[957];
	v[169] = ei1B[2] * v[172] + v[170] * v[176] + v[137] * v[627];
	v[1048] = v[140] - v[169];
	v[1028] = v[140] + v[169];
	v[901] = v[169] * v[347];
	v[895] = v[169] * v[346];
	v[722] = v[1052] + v[901] + v[903];
	v[231] = v[169] * v[95];
	v[173] = ei1B[1] * v[172] + v[170] * v[174] + v[137] * v[640];
	v[1023] = 0.5e0*v[173];
	v[1049] = -v[1023] + v[1024];
	v[925] = -(v[173] * v[964]);
	v[1054] = v[925] + v[927];
	v[242] = v[173] * v[962];
	v[177] = v[174] * v[180] + v[176] * v[184] + v[137] * v[653];
	v[1027] = v[152] + v[177];
	v[920] = v[177] * v[347];
	v[1057] = v[920] + v[922];
	v[928] = (-v[152] + v[177] + v[1021] * v[538] + v[194] * v[546] + v[195] * v[554] + v[196] * v[562] + v[197] * v[570]
		+ v[1025] * v[578] + v[199] * v[639] + v[1026] * v[652] + v[201] * v[665] + v[1047] * v[944] + v[1028] * v[948] + (-v[1022]
		- v[1023])*v[950])*v[953];
	v[758] = v[1057] + v[919] + v[921];
	v[755] = v[1027] * v[346] + (v[569] + v[663])*v[95];
	v[240] = v[177] * v[95];
	v[186] = v[180] * v[181] + v[184] * v[185] + v[137] * v[666];
	v[1031] = v[165] + v[186];
	v[869] = v[186] * v[259];
	v[1046] = v[869] + v[871];
	v[809] = v[1046] + v[868] + v[870];
	v[806] = v[1031] * v[258] + (v[614] + v[670])*v[70];
	v[226] = v[186] * v[70];
	v[187] = v[170] * v[185] + v[180] * v[188] + v[137] * v[679];
	v[1034] = v[159] - v[187];
	v[1032] = v[159] + v[187];
	v[838] = v[187] * v[259];
	v[832] = v[187] * v[258];
	v[773] = v[1038] + v[838] + v[840];
	v[217] = v[187] * v[70];
	v[189] = v[170] * v[181] + v[184] * v[188] + v[137] * v[692];
	v[1020] = -0.5e0*v[189];
	v[1035] = v[1020] + v[1030];
	v[874] = -(v[189] * v[960]);
	v[1040] = v[874] + v[876];
	v[877] = (-v[165] + v[186] + v[1017] * v[583] + v[207] * v[591] + v[208] * v[599] + v[209] * v[607] + v[210] * v[615]
		+ v[1018] * v[623] + v[212] * v[672] + v[213] * v[685] + v[1019] * v[698] + v[1033] * v[934] + v[1032] * v[938] + (v[1020]
		- v[1029])*v[940])*v[943];
	v[228] = v[1020] * v[70];
	v[190] = v[230] + v[231];
	v[191] = v[234] + v[235];
	v[192] = v[239] + v[240];
	v[202] = v[1021] * v[136] + v[140] * v[194] + v[143] * v[195] + v[148] * v[196] + v[152] * v[197] + v[1022] * v[198]
		+ v[169] * v[199] + v[1023] * v[200] + v[177] * v[201];
	v[898] = (8e0*v[202]) / Power(v[232], 3);
	v[923] = v[898] * v[950];
	v[1053] = v[923] + v[928];
	v[916] = v[1053] + v[538] * v[962] - v[136] * v[964];
	v[910] = v[1024] * v[346] + v[898] * v[947] + (-v[143] + v[148] + v[1021] * v[537] + v[194] * v[545] + v[195] * v[553]
		+ v[196] * v[561] + v[197] * v[569] + v[1025] * v[577] + v[199] * v[637] + v[1026] * v[650] + v[201] * v[663]
		+ v[1027] * v[944] + (-v[1022] + v[1024])*v[947] + v[1028] * v[951])*v[953] + v[537] * v[962];
	v[237] = v[202] * v[953];
	v[918] = v[237] + v[241] + v[242];
	v[905] = v[236] - v[242] + v[918];
	v[888] = -v[241] + v[242] + v[905];
	v[203] = v[216] + v[217];
	v[204] = v[220] + v[221];
	v[205] = v[225] + v[226];
	v[215] = v[1017] * v[158] + v[159] * v[207] + v[160] * v[208] + v[163] * v[209] + v[165] * v[210] + v[1029] * v[211]
		+ v[186] * v[212] + v[187] * v[213] - v[1020] * v[214];
	v[835] = (8e0*v[215]) / Power(v[218], 3);
	v[872] = v[835] * v[940];
	v[1039] = v[872] + v[877];
	v[856] = v[1039] + v[583] * v[957] - v[158] * v[960];
	v[850] = v[1030] * v[258] + v[835] * v[937] + (-v[160] + v[163] + v[1017] * v[582] + v[207] * v[590] + v[208] * v[598]
		+ v[209] * v[606] + v[210] * v[614] + v[1018] * v[622] + v[212] * v[670] + v[213] * v[683] + v[1019] * v[696]
		+ v[1031] * v[934] + (-v[1029] + v[1030])*v[937] + v[1032] * v[941])*v[943] + v[582] * v[957];
	v[223] = v[215] * v[943];
	v[867] = v[223] + v[227] + v[228];
	v[845] = v[222] - v[228] + v[867];
	v[825] = -v[227] + v[228] + v[845];
	residual[0] = v[216] - v[217] + v[825] * v[935] + v[205] * v[938] + v[204] * v[941];
	residual[1] = -v[220] + v[221] + v[205] * v[934] + v[845] * v[937] + v[203] * v[941];
	residual[2] = -v[225] + v[226] + v[204] * v[934] + v[203] * v[938] + v[867] * v[940];
	residual[3] = v[230] - v[231] + v[888] * v[945] + v[192] * v[948] + v[191] * v[951];
	residual[4] = -v[234] + v[235] + v[192] * v[944] + v[905] * v[947] + v[190] * v[951];
	residual[5] = -v[239] + v[240] + v[191] * v[944] + v[190] * v[948] + v[918] * v[950];
	stiffness[0][0] = v[1034] * v[257] + (v[589] - v[681])*v[70] + 2e0*v[825] + (v[1031] * v[257] + (v[613] + v[668]
		)*v[70])*v[938] + (v[1033] * v[257] + (v[597] + v[605])*v[70])*v[941] + v[935] * (v[1035] * v[257] - 0.5e0*(v[581]
		+ v[694])*v[70] + v[835] * v[935] + (v[1034] + v[1017] * v[581] + v[207] * v[589] + v[208] * v[597] + v[209] * v[605]
		+ v[210] * v[613] + v[1018] * (v[137] * (-(v[166] * v[335]) + ei3A[0] * v[471] + v[303] * v[619]) + v[498] * v[620])
		+ v[212] * v[668] + v[213] * v[681] + v[1019] * v[694] + v[1035] * v[935] + v[1031] * v[938] + v[1033] * v[941])*v[943]);
	stiffness[0][1] = v[1036] + 0.5e0*v[205] - v[832] - v[834] + v[806] * v[938] + (v[1037] + v[846] + v[848])*v[941]
		+ v[935] * (v[1020] * v[258] + v[850] + v[696] * v[957]);
	stiffness[0][2] = v[1038] + 0.5e0*v[204] - v[838] - v[840] + (v[1040] + v[856])*v[935] + v[809] * v[938]
		+ v[791] * v[941];
	stiffness[0][3] = -v[774] + v[775] + (v[1041] + v[858])*v[935] + v[812] * v[938] + v[794] * v[941];
	stiffness[0][4] = -v[777] + v[778] + (v[1042] + v[861])*v[935] + v[815] * v[938] + v[797] * v[941];
	stiffness[0][5] = -v[780] + v[781] + (v[1043] + v[864])*v[935] + v[818] * v[938] + v[800] * v[941];
	stiffness[1][1] = v[1037] + 2e0*v[845] - v[846] - v[848] + v[806] * v[934] + (v[1036] + v[832] + v[834])*v[941]
		+ v[937] * (-(v[1029] * v[258]) + v[850] + v[622] * v[957]);
	stiffness[1][2] = v[1044] + 0.5e0*v[203] - v[852] - v[854] + v[809] * v[934] + (v[1045] + v[856])*v[937]
		+ v[773] * v[941];
	stiffness[1][3] = v[792] - v[793] + v[812] * v[934] + (v[858] + v[859] + v[879])*v[937] + v[776] * v[941];
	stiffness[1][4] = v[795] - v[796] + v[815] * v[934] + (v[861] + v[862] + v[882])*v[937] + v[779] * v[941];
	stiffness[1][5] = v[798] - v[799] + v[818] * v[934] + (v[864] + v[865] + v[885])*v[937] + v[782] * v[941];
	stiffness[2][2] = v[1046] + 2e0*v[867] - v[868] - v[870] + v[791] * v[934] + v[773] * v[938] + (v[1039] + v[1040]
		+ v[1045])*v[940];
	stiffness[2][3] = v[810] - v[811] + v[794] * v[934] + v[776] * v[938] + (v[1041] + v[879])*v[940];
	stiffness[2][4] = v[813] - v[814] + v[797] * v[934] + v[779] * v[938] + (v[1042] + v[882])*v[940];
	stiffness[2][5] = v[816] - v[817] + v[800] * v[934] + v[782] * v[938] + (v[1043] + v[885])*v[940];
	stiffness[3][3] = v[1048] * v[345] + 2e0*v[888] + (v[544] - v[635])*v[95] + v[948] * (v[1027] * v[345] + (v[568]
		+ v[661])*v[95]) + (v[1047] * v[345] + (v[552] + v[560])*v[95])*v[951] + v[945] * (v[1049] * v[345] + v[898] * v[945]
		- 0.5e0*(v[536] + v[648])*v[95] + (v[1048] + v[1021] * v[536] + v[194] * v[544] + v[195] * v[552] + v[196] * v[560]
		+ v[197] * v[568] + v[1025] * v[503] * v[572] + v[199] * v[635] + v[1026] * v[648] + v[201] * v[661] + v[1049] * v[945]
		+ v[1027] * v[948] + v[1047] * v[951])*v[953]);
	stiffness[3][4] = v[1050] + 0.5e0*v[192] - v[895] - v[897] + v[755] * v[948] + (v[1051] + v[906] + v[908])*v[951]
		+ v[945] * (-(v[1023] * v[346]) + v[910] + v[650] * v[962]);
	stiffness[3][5] = v[1052] + 0.5e0*v[191] - v[901] - v[903] + (v[1054] + v[916])*v[945] + v[758] * v[948]
		+ v[740] * v[951];
	stiffness[4][4] = v[1051] + 2e0*v[905] - v[906] - v[908] + v[755] * v[944] + (v[1050] + v[895] + v[897])*v[951]
		+ v[947] * (-(v[1022] * v[346]) + v[910] + v[577] * v[962]);
	stiffness[4][5] = v[1055] + 0.5e0*v[190] - v[912] - v[914] + v[758] * v[944] + (v[1056] + v[916])*v[947]
		+ v[722] * v[951];
	stiffness[5][5] = v[1057] + 2e0*v[918] - v[919] - v[921] + v[740] * v[944] + v[722] * v[948] + (v[1053] + v[1054]
		+ v[1056])*v[950];
	for (i01 = 1; i01<6; i01++){
		for (i02 = 0; i02<i01; i02++){
			stiffness[i01][i02] = stiffness[i02][i01];
		}
	};
	(*thetad) = v[126];
};


//Calcula contribuições do amortecedor de torção
void Hinge::EvaluateTorsionDamping(double *v, double *residual
	, double **stiffness, double *alphaA, double *alphaB, double *omegaiA
	, double *omegaiB, double *domegaiA, double *domegaiB, double *ei3A, double
	(*dampc1), double(*dampc2), double(*alpha4), double(*alpha5), double(*alpha6))
{
	v[392] = (*alpha4)*alphaA[0] + (*alpha6)*domegaiA[0] + (*alpha5)*omegaiA[0];
	v[391] = (*alpha4)*alphaA[1] + (*alpha6)*domegaiA[1] + (*alpha5)*omegaiA[1];
	v[390] = (*alpha4)*alphaA[2] + (*alpha6)*domegaiA[2] + (*alpha5)*omegaiA[2];
	v[389] = (*alpha4)*alphaB[0] + (*alpha6)*domegaiB[0] + (*alpha5)*omegaiB[0];
	v[388] = (*alpha4)*alphaB[1] + (*alpha6)*domegaiB[1] + (*alpha5)*omegaiB[1];
	v[387] = (*alpha4)*alphaB[2] + (*alpha6)*domegaiB[2] + (*alpha5)*omegaiB[2];
	v[386] = 0.5e0*alphaB[2];
	v[384] = Power(alphaB[2], 2);
	v[383] = alphaB[0] * v[386];
	v[382] = 0.5e0*alphaB[1];
	v[381] = 2e0*alphaB[2];
	v[380] = Power(alphaB[1], 2);
	v[379] = alphaB[0] * v[382];
	v[378] = 2e0*alphaB[1];
	v[377] = Power(alphaB[0], 2);
	v[376] = 2e0*alphaB[0];
	v[375] = 0.5e0*alphaA[2];
	v[373] = Power(alphaA[2], 2);
	v[372] = alphaA[0] * v[375];
	v[371] = 0.5e0*alphaA[1];
	v[370] = 2e0*alphaA[2];
	v[369] = Power(alphaA[1], 2);
	v[368] = alphaA[0] * v[371];
	v[367] = 2e0*alphaA[1];
	v[366] = Power(alphaA[0], 2);
	v[365] = 2e0*alphaA[0];
	v[214] = -v[366] - v[369];
	v[189] = alphaA[2] + v[368];
	v[179] = -alphaA[2] + v[368];
	v[100] = alphaA[2] * v[371];
	v[209] = alphaA[0] + v[100];
	v[198] = -alphaA[0] + v[100];
	v[205] = -alphaA[1] + v[372];
	v[185] = alphaA[1] + v[372];
	v[193] = -v[366] - v[373];
	v[174] = -v[369] - v[373];
	v[169] = 4e0 + v[366] + v[369] + v[373];
	v[171] = 1e0 / Power(v[169], 2);
	v[374] = -4e0*v[171];
	v[173] = v[370] * v[374];
	v[395] = 0.5e0*v[173];
	v[218] = v[214] * v[395];
	v[172] = v[367] * v[374];
	v[394] = 0.5e0*v[172];
	v[223] = -(v[172] * v[375]);
	v[195] = v[193] * v[394];
	v[170] = v[365] * v[374];
	v[393] = 0.5e0*v[170];
	v[225] = v[170] * v[371];
	v[222] = -(v[170] * v[375]);
	v[175] = v[174] * v[393];
	v[270] = -v[377] - v[380];
	v[260] = alphaB[2] + v[379];
	v[245] = -alphaB[2] + v[379];
	v[128] = alphaB[2] * v[382];
	v[272] = alphaB[0] + v[128];
	v[256] = -alphaB[0] + v[128];
	v[274] = -alphaB[1] + v[383];
	v[242] = alphaB[1] + v[383];
	v[258] = -v[377] - v[384];
	v[247] = -v[380] - v[384];
	v[229] = 4e0 + v[377] + v[380] + v[384];
	v[231] = 1e0 / Power(v[229], 2);
	v[385] = -4e0*v[231];
	v[233] = v[381] * v[385];
	v[397] = 0.5e0*v[233];
	v[232] = v[378] * v[385];
	v[398] = 0.5e0*v[232];
	v[235] = -(v[232] * v[386]);
	v[230] = v[376] * v[385];
	v[399] = 0.5e0*v[230];
	v[237] = v[230] * v[382];
	v[234] = -(v[230] * v[386]);
	v[87] = 4e0 / v[169];
	v[226] = -0.5e0*v[87];
	v[228] = v[226] - alphaA[0] * v[393];
	v[227] = -v[226] + alphaA[1] * v[394];
	v[224] = v[226] - alphaA[2] * v[395];
	v[216] = v[226] * v[367];
	v[217] = v[216] + v[214] * v[394];
	v[213] = v[226] * v[365];
	v[215] = v[213] + v[214] * v[393];
	v[210] = v[170] * v[209] + v[87];
	v[207] = v[172] * v[205] - v[87];
	v[199] = v[170] * v[198] - v[87];
	v[196] = v[226] * v[370];
	v[197] = v[196] + v[193] * v[395];
	v[194] = v[213] + v[193] * v[393];
	v[192] = v[173] * v[189] + v[87];
	v[187] = v[172] * v[185] + v[87];
	v[184] = alphaA[2] * v[226];
	v[211] = -v[184] + v[172] * v[209];
	v[220] = ei3A[0] * v[207] + ei3A[1] * v[211] + ei3A[2] * v[217];
	v[206] = -v[184] + v[170] * v[205];
	v[219] = ei3A[0] * v[206] + ei3A[1] * v[210] + ei3A[2] * v[215];
	v[200] = -v[184] + v[172] * v[198];
	v[186] = -v[184] + v[170] * v[185];
	v[183] = v[173] * v[179] - v[87];
	v[181] = alphaA[0] * v[226];
	v[208] = -v[181] + v[173] * v[205];
	v[191] = -v[181] + v[172] * v[189];
	v[203] = ei3A[0] * v[191] + ei3A[1] * v[195] + ei3A[2] * v[200];
	v[188] = -v[181] + v[173] * v[185];
	v[182] = v[172] * v[179] - v[181];
	v[178] = -(alphaA[1] * v[226]);
	v[212] = v[178] + v[173] * v[209];
	v[221] = ei3A[0] * v[208] + ei3A[1] * v[212] + ei3A[2] * v[218];
	v[201] = v[178] + v[173] * v[198];
	v[204] = ei3A[0] * v[192] + ei3A[1] * v[197] + ei3A[2] * v[201];
	v[190] = v[178] + v[170] * v[189];
	v[202] = ei3A[0] * v[190] + ei3A[1] * v[194] + ei3A[2] * v[199];
	v[180] = v[178] + v[170] * v[179];
	v[283] = ei3A[0] * v[175] + ei3A[1] * v[180] + ei3A[2] * v[186];
	v[177] = v[196] + v[174] * v[395];
	v[285] = ei3A[0] * v[177] + ei3A[1] * v[183] + ei3A[2] * v[188];
	v[176] = v[216] + v[174] * v[394];
	v[284] = ei3A[0] * v[176] + ei3A[1] * v[182] + ei3A[2] * v[187];
	v[90] = 1e0 - v[174] * v[226];
	v[91] = v[179] * v[87];
	v[92] = v[185] * v[87];
	v[94] = v[189] * v[87];
	v[96] = 1e0 - v[193] * v[226];
	v[97] = v[198] * v[87];
	v[149] = ei3A[0] * v[94] + ei3A[1] * v[96] + ei3A[2] * v[97];
	v[407] = v[149] * v[172];
	v[99] = v[205] * v[87];
	v[101] = v[209] * v[87];
	v[102] = 1e0 - v[214] * v[226];
	v[150] = ei3A[1] * v[101] + ei3A[2] * v[102] + ei3A[0] * v[99];
	v[409] = v[150] * v[173];
	v[115] = 4e0 / v[229];
	v[396] = -0.5e0*v[115];
	v[276] = v[378] * v[396];
	v[269] = v[376] * v[396];
	v[266] = v[381] * v[396];
	v[250] = alphaB[0] * v[396];
	v[244] = -(alphaB[1] * v[396]);
	v[241] = alphaB[2] * v[396];
	v[240] = v[396] - alphaB[0] * v[399];
	v[239] = v[232] * v[382] - v[396];
	v[236] = -(v[233] * v[386]) + v[396];
	v[118] = 1e0 - v[247] * v[396];
	v[119] = v[115] * v[245];
	v[120] = v[115] * v[242];
	v[122] = v[115] * v[260];
	v[124] = 1e0 - v[258] * v[396];
	v[125] = v[115] * v[256];
	v[127] = v[115] * v[274];
	v[129] = v[115] * v[272];
	v[130] = 1e0 - v[270] * v[396];
	v[291] = -(v[120] * v[387]) - v[119] * v[388] - v[118] * v[389] + v[392] * v[90] + v[391] * v[91] + v[390] * v[92];
	v[289] = -(v[125] * v[387]) - v[124] * v[388] - v[122] * v[389] + v[392] * v[94] + v[391] * v[96] + v[390] * v[97];
	v[290] = -(v[130] * v[387]) - v[129] * v[388] - v[127] * v[389] + v[102] * v[390] + v[101] * v[391] + v[392] * v[99];
	v[146] = ei3A[0] * v[90] + ei3A[1] * v[91] + ei3A[2] * v[92];
	v[412] = -(v[115] * v[150]) - v[146] * v[244] - v[149] * v[250];
	v[411] = -(v[115] * v[149]) - v[146] * v[241] + v[150] * v[250];
	v[410] = -(v[115] * v[146]) + v[149] * v[241] + v[150] * v[244];
	v[408] = v[146] * v[178] + v[149] * v[181] + v[150] * v[87];
	v[406] = -(v[150] * v[181]) + v[146] * v[184] + v[149] * v[87];
	v[405] = -(v[150] * v[178]) - v[149] * v[184] + v[146] * v[87];
	v[404] = v[146] * v[170];
	v[312] = -(v[150] * ((*alpha4)*v[130] + (v[244] + v[233] * v[272])*v[388] + (-v[250] + v[233] * v[274])*v[389]
		+ v[270] * v[387] * v[397])) - v[146] * ((*alpha4)*v[120] + (v[233] * v[242] - v[250])*v[387] + (-v[115]
		+ v[233] * v[245])*v[388] + v[389] * (v[266] + v[247] * v[397])) - v[149] * ((*alpha4)*v[125] + (v[244]
		+ v[233] * v[256])*v[387] + (v[115] + v[233] * v[260])*v[389] + v[388] * (v[266] + v[258] * v[397]));
	v[308] = -(v[149] * ((*alpha4)*v[124] + (-v[241] + v[232] * v[256])*v[387] + (-v[250] + v[232] * v[260])*v[389]
		+ v[258] * v[388] * v[398])) - v[146] * ((*alpha4)*v[119] + (v[115] + v[232] * v[242])*v[387] + (v[232] * v[245]
		- v[250])*v[388] + v[389] * (v[276] + v[247] * v[398])) - v[150] * ((*alpha4)*v[129] + (-v[241] + v[232] * v[272]
		)*v[388] + (-v[115] + v[232] * v[274])*v[389] + v[387] * (v[276] + v[270] * v[398]));
	v[304] = -(v[146] * ((*alpha4)*v[118] + (-v[241] + v[230] * v[242])*v[387] + (v[244] + v[230] * v[245])*v[388]
		+ v[247] * v[389] * v[399])) - v[149] * ((*alpha4)*v[122] + (-v[115] + v[230] * v[256])*v[387] + (v[244]
		+ v[230] * v[260])*v[389] + v[388] * (v[269] + v[258] * v[399])) - v[150] * ((*alpha4)*v[127] + (v[115]
		+ v[230] * v[272])*v[388] + (-v[241] + v[230] * v[274])*v[389] + v[387] * (v[269] + v[270] * v[399]));
	v[300] = v[204] * v[289] + v[221] * v[290] + v[285] * v[291] + v[150] * ((*alpha4)*v[102] + v[218] * v[390]
		+ v[212] * v[391] + v[208] * v[392]) + v[146] * (v[188] * v[390] + v[183] * v[391] + v[177] * v[392] + (*alpha4)*v[92])
		+ v[149] * (v[201] * v[390] + v[197] * v[391] + v[192] * v[392] + (*alpha4)*v[97]);
	v[296] = v[203] * v[289] + v[220] * v[290] + v[284] * v[291] + v[150] * ((*alpha4)*v[101] + v[217] * v[390]
		+ v[211] * v[391] + v[207] * v[392]) + v[146] * (v[187] * v[390] + v[182] * v[391] + v[176] * v[392] + (*alpha4)*v[91])
		+ v[149] * (v[200] * v[390] + v[195] * v[391] + v[191] * v[392] + (*alpha4)*v[96]);
	v[292] = v[202] * v[289] + v[219] * v[290] + v[283] * v[291] + v[146] * (v[186] * v[390] + v[180] * v[391] + v[175] * v[392] +
		(*alpha4)*v[90]) + v[149] * (v[199] * v[390] + v[194] * v[391] + v[190] * v[392] + (*alpha4)*v[94]) + v[150] *
		(v[215] * v[390] + v[210] * v[391] + v[206] * v[392] + (*alpha4)*v[99]);
	v[147] = v[149] * v[289] + v[150] * v[290] + v[146] * v[291];
	v[315] = (*dampc2)*v[147] * _copysign(1.e0, v[147]);
	v[313] = fabs(v[147]);
	v[400] = (*dampc1) + (*dampc2)*v[313] + v[315];
	v[320] = v[312] * v[400];
	v[319] = v[308] * v[400];
	v[318] = v[304] * v[400];
	v[317] = v[300] * v[400];
	v[316] = v[296] * v[400];
	v[314] = v[292] * v[400];
	v[148] = v[147] * (-v[315] + v[400]);
	v[403] = v[148] * v[150];
	v[402] = v[148] * v[149];
	v[401] = v[146] * v[148];
	v[361] = -(v[235] * v[401]);
	v[360] = v[234] * v[402];
	v[352] = v[237] * v[403];
	v[337] = v[223] * v[401];
	v[336] = -(v[222] * v[402]);
	v[328] = -(v[225] * v[403]);
	residual[0] = v[148] * v[405];
	residual[1] = v[148] * v[406];
	residual[2] = v[148] * v[408];
	residual[3] = v[148] * v[410];
	residual[4] = v[148] * v[411];
	residual[5] = v[148] * v[412];
	stiffness[0][0] = v[328] + v[336] + v[148] * v[404] + v[314] * v[405];
	stiffness[0][1] = v[316] * v[405] + v[148] * (v[146] * v[172] - v[184] * v[203] - v[178] * v[220] - v[149] * v[223]
		- v[150] * v[227] + v[284] * v[87]);
	stiffness[0][2] = v[317] * v[405] + v[148] * (v[146] * v[173] - v[184] * v[204] - v[178] * v[221] + v[150] * v[223]
		- v[149] * v[224] + v[285] * v[87]);
	stiffness[0][3] = v[318] * v[405];
	stiffness[0][4] = v[319] * v[405];
	stiffness[0][5] = v[320] * v[405];
	stiffness[1][0] = v[314] * v[406] + v[148] * (v[149] * v[170] - v[181] * v[219] + v[146] * v[222] - v[150] * v[228]
		+ v[184] * v[283] + v[202] * v[87]);
	stiffness[1][1] = -v[328] + v[337] + v[316] * v[406] + v[148] * v[407];
	stiffness[1][2] = v[317] * v[406] + v[148] * (v[149] * v[173] - v[181] * v[221] - v[150] * v[222] + v[146] * v[224]
		+ v[184] * v[285] + v[204] * v[87]);
	stiffness[1][3] = v[318] * v[406];
	stiffness[1][4] = v[319] * v[406];
	stiffness[1][5] = v[320] * v[406];
	stiffness[2][0] = v[314] * v[408] + v[148] * (v[150] * v[170] + v[181] * v[202] + v[146] * v[225] + v[149] * v[228]
		+ v[178] * v[283] + v[219] * v[87]);
	stiffness[2][1] = v[316] * v[408] + v[148] * (v[150] * v[172] + v[181] * v[203] - v[149] * v[225] + v[146] * v[227]
		+ v[178] * v[284] + v[220] * v[87]);
	stiffness[2][2] = -v[336] - v[337] + v[317] * v[408] + v[148] * v[409];
	stiffness[2][3] = v[318] * v[408];
	stiffness[2][4] = v[319] * v[408];
	stiffness[2][5] = v[320] * v[408];
	stiffness[3][0] = v[148] * (v[202] * v[241] + v[219] * v[244] - v[115] * v[283]) + v[314] * v[410];
	stiffness[3][1] = v[148] * (v[203] * v[241] + v[220] * v[244] - v[115] * v[284]) + v[316] * v[410];
	stiffness[3][2] = v[148] * (v[204] * v[241] + v[221] * v[244] - v[115] * v[285]) + v[317] * v[410];
	stiffness[3][3] = v[352] + v[360] - v[230] * v[401] + v[318] * v[410];
	stiffness[3][4] = v[148] * (-(v[146] * v[232]) + v[149] * v[235] + v[150] * v[239]) + v[319] * v[410];
	stiffness[3][5] = v[148] * (-(v[146] * v[233]) - v[150] * v[235] + v[149] * v[236]) + v[320] * v[410];
	stiffness[4][0] = v[148] * (-(v[115] * v[202]) + v[219] * v[250] - v[241] * v[283]) + v[314] * v[411];
	stiffness[4][1] = v[148] * (-(v[115] * v[203]) + v[220] * v[250] - v[241] * v[284]) + v[316] * v[411];
	stiffness[4][2] = v[148] * (-(v[115] * v[204]) + v[221] * v[250] - v[241] * v[285]) + v[317] * v[411];
	stiffness[4][3] = v[148] * (-(v[149] * v[230]) - v[146] * v[234] + v[150] * v[240]) + v[318] * v[411];
	stiffness[4][4] = -v[352] + v[361] - v[232] * v[402] + v[319] * v[411];
	stiffness[4][5] = v[148] * (-(v[149] * v[233]) + v[150] * v[234] - v[146] * v[236]) + v[320] * v[411];
	stiffness[5][0] = v[148] * (-(v[115] * v[219]) - v[202] * v[250] - v[244] * v[283]) + v[314] * v[412];
	stiffness[5][1] = v[148] * (-(v[115] * v[220]) - v[203] * v[250] - v[244] * v[284]) + v[316] * v[412];
	stiffness[5][2] = v[148] * (-(v[115] * v[221]) - v[204] * v[250] - v[244] * v[285]) + v[317] * v[412];
	stiffness[5][3] = v[148] * (-(v[150] * v[230]) - v[146] * v[237] - v[149] * v[240]) + v[318] * v[412];
	stiffness[5][4] = v[148] * (-(v[150] * v[232]) + v[149] * v[237] - v[146] * v[239]) + v[319] * v[412];
	stiffness[5][5] = -v[360] - v[361] - v[233] * v[403] + v[320] * v[412];
};




