#include "Plane.h"
#include "Database.h"
#include "Matrix.h"
//Variáveis globais
extern
Database db;

Plane::Plane()
{
	N = new Matrix(3);
	T1 = new Matrix(3);
	T2 = new Matrix(3);
	P = new Matrix(3);
}


Plane::~Plane()
{
	delete N;
	delete T1;
	delete T2;
	delete P;
}

void Plane::WriteVTK_XMLRender(FILE *f)
{
	int n_points = 4;	//quatro pontos dos cantos da supefície
	int n_cells = 1;	//quadrilátero
	Matrix point1(3);
	Matrix point2(3);
	Matrix point3(3);
	Matrix point4(3);
	//Opens Piece
	fprintf(f, "\t\t<Piece NumberOfPoints = \"%d\" NumberOfCells = \"%d\">\n", n_points, n_cells);
	//Opens Points
	fprintf(f, "\t\t\t<Points>\n");
	//Opens DataArray
	fprintf(f, "\t\t\t\t<DataArray type = \"Float32\" NumberOfComponents = \"3\" Format = \"ascii\">\n");

	//Posição de cada ponto P no plano xy (referência)
	point1 = *P + len*(*T1) + len*(*T2);
	point2 = *P - len*(*T1) + len*(*T2);
	point3 = *P - len*(*T1) - len*(*T2);
	point4 = *P + len*(*T1) - len*(*T2);
	
	fprintf(f, "\t\t\t\t\t%f\t%f\t%f\t\n", point1(0, 0), point1(1, 0), point1(2, 0));
	fprintf(f, "\t\t\t\t\t%f\t%f\t%f\t\n", point2(0, 0), point2(1, 0), point2(2, 0));
	fprintf(f, "\t\t\t\t\t%f\t%f\t%f\t\n", point3(0, 0), point3(1, 0), point3(2, 0));
	fprintf(f, "\t\t\t\t\t%f\t%f\t%f\t\n", point4(0, 0), point4(1, 0), point4(2, 0));
	//Closes DataArray
	fprintf(f, "\t\t\t\t</DataArray>\n");
	//Closes Points
	fprintf(f, "\t\t\t</Points>\n");

	//Opens Cells
	fprintf(f, "\t\t\t<Cells>\n");
	//Opens DataArray
	fprintf(f, "\t\t\t\t<DataArray type = \"Int32\" Name = \"connectivity\" Format = \"ascii\">\n");
	//Linhas
	fprintf(f, "\t\t\t\t\t%d\t%d\t%d\t%d\n", 0,1,2,3);
	//Closes DataArray
	fprintf(f, "\t\t\t\t</DataArray>\n");
	//Opens DataArray
	fprintf(f, "\t\t\t\t<DataArray type = \"Int32\" Name = \"types\" Format = \"ascii\">\n");
	fprintf(f, "\t\t\t\t\t%d\n", 9);
	//Closes DataArray
	fprintf(f, "\t\t\t\t</DataArray>\n");
	//Opens DataArray
	fprintf(f, "\t\t\t\t<DataArray type = \"Int32\" Name = \"offsets\" Format = \"ascii\">\n");
	fprintf(f, "\t\t\t\t\t%d\n", 4);
	//Closes DataArray
	fprintf(f, "\t\t\t\t</DataArray>\n");
	//Closes Cells
	fprintf(f, "\t\t\t</Cells>\n");
	//Closes Piece
	fprintf(f, "\t\t</Piece>\n");
}

bool Plane::Read(FILE *f)
{
	char s[1000];
	fscanf(f, "%s", s);
	number = atoi(s);

	fscanf(f, "%s", s);
	if (!strcmp(s, "N"))
	{
		fscanf(f, "%s", s);
		(*N)(0, 0) = atof(s);
		fscanf(f, "%s", s);
		(*N)(1, 0) = atof(s);
		fscanf(f, "%s", s);
		(*N)(2, 0) = atof(s);
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "P"))
	{
		fscanf(f, "%s", s);
		(*P)(0, 0) = atof(s);
		fscanf(f, "%s", s);
		(*P)(1, 0) = atof(s);
		fscanf(f, "%s", s);
		(*P)(2, 0) = atof(s);
	}
	else
		return false;
	return true;
}

void Plane::Write(FILE *f)
{
	fprintf(f, "Plane\t%d\tN\t%.6e\t%.6e\t%.6e\tP\t%.6e\t%.6e\t%.6e\n",
		number,
		(*N)(0, 0),
		(*N)(1, 0),
		(*N)(2, 0),
		(*P)(0, 0),
		(*P)(1, 0),
		(*P)(2, 0));
}

//Calculates the distance between the Analytical Surface and a given point in space
double Plane::Distance(Matrix &O)
{
	if (O.getLines() != 3)
	{
		printf("Error - point has to be dimension 3\n");
		return 0;
	}
	else	//Returns distance
		return (dot(O - *P, *N));
}
//Calculates and returns the normal direction of the closest point on analytical suface to point O in space
void Plane::N_O(Matrix *O, Matrix *Normal)
{
	(*Normal) = (*N);
}
//Calculates and returns the T1 direction of the closest point on analytical suface to point O in space
void Plane::T1_O(Matrix *O, Matrix *Tang1)
{
	(*Tang1) = (*T1);
}
//Calculates and returns the T2 direction of the closest point on analytical suface to point O in space
void Plane::T2_O(Matrix *O, Matrix *Tang2)
{
	(*Tang2) = (*T2);
}
void Plane::PreCalc()
{
	//Calculando t1 e t2, obrigatoriamente ortogonais a n
	//Passo 1 - normalização da normal
	(*N) = 1.0 / (norm(*N))*(*N);
	//Passo 2 - definição de direção t1
	Matrix t1_temp(3, 1);

	//Chute - componentes do vetor t1
	t1_temp(0, 0) = 1.0;
	t1_temp(1, 0) = 0.0;
	t1_temp(2, 0) = 0.0;

	double min_norm = 1e-3;
	//Se existe alguma componente LD com n
	if (dot(*N, t1_temp) != 0.0)
	{
		//Calcula componente ortho a n e subtrai de t1_temp
		t1_temp = t1_temp - dot(t1_temp, *N)*(*N);

		//Verifica se sobrou vetor grande o bastante, maior que tol
		if (norm(t1_temp) >= min_norm)
			t1_temp = 1.0 / (norm(t1_temp))*t1_temp;
		else
		{
			//Dá um chute melhor e faz o processo
			t1_temp(0, 0) = 0.0;
			t1_temp(1, 0) = 1.0;
			t1_temp(2, 0) = 0.0;
			//Calcula componente ortho a n e subtrai de t1_temp
			t1_temp = t1_temp - dot(t1_temp, *N)*(*N);
			t1_temp = 1.0 / (norm(t1_temp))*t1_temp;
		}

	}
	(*T1) = t1_temp;
	(*T2) = cross(*N, *T1);
	len = db.EvaluateBoundingBoxDiag();
}