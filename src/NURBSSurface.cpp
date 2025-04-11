#include "NURBSSurface.h"
#include"Database.h"
#include <cmath>
#include <vector>
#include "Particle.h"
#include "NURBSParticle.h"
#include "Encoding.h"
#include <Eigen/dense>
#include "MatrixFloat.h"

#include <iostream>

//Variaveis globais
extern
Database db;

using namespace std;

NURBSSurface::NURBSSurface()
{
	G = new Matrix(3);
	J_O = new Matrix(3, 3);

	number = 0;

	//Internal and storing variables
	U_dim = 0;
	U_order = 0;
	V_dim = 0;
	V_order = 0;
	U_knot_vector = NULL;
	V_knot_vector = NULL;
	weights = NULL;
	control_points = NULL;
	Pw = NULL;
	Bin = NULL;
	//Marina
	//n_sub = NULL;
	parameter_sub = NULL;
	box_sub_coord = NULL;
	box_sub_center = NULL;
	box_sub_points = NULL;
	halfedge_lengths = NULL;
}


NURBSSurface::~NURBSSurface()
{
	delete G;
	delete J_O;

	if (U_knot_vector != NULL)
		delete[] U_knot_vector;
	if (V_knot_vector != NULL)
		delete[] V_knot_vector;
	if (weights != NULL)
	{
		for (int i = 0; i < U_dim; i++)
			delete[] weights[i];
		delete[]weights;
	}

	if (control_points != NULL)
	{
		for (int i = 0; i < U_dim; i++)
			delete[] control_points[i];
		delete[] control_points;
	}
	if (Pw != NULL)
	{
		for (int i = 0; i < U_dim; i++)
			delete[] Pw[i];
		delete[] Pw;
	}

	if (Bin != NULL)
	{
		for (int i = 0; i < U_order + V_order + 4; i++)
			delete[] Bin[i];
		delete[] Bin;
	}

	//Marina
	//delete n_sub;
	if (parameter_sub != NULL)
	{
		for (int i = 0; i < subdivisions[0]; i++)
		{
			for (int j = 0; j < subdivisions[1]; j++)
				delete[]parameter_sub[i][j];
			delete[] parameter_sub[i];
		}
		delete[]parameter_sub;
	}
	if (box_sub_coord != NULL)
	{
		for (int i = 0; i < subdivisions[0]; i++)
		{
			for (int j = 0; j < subdivisions[1]; j++)
				delete[]box_sub_coord[i][j];
			delete[] box_sub_coord[i];
		}
		delete[]box_sub_coord;
	}
	if (box_sub_center != NULL)
	{
		for (int i = 0; i < subdivisions[0]; i++)
		{
			delete[] box_sub_center[i];
		}
		delete[]box_sub_center;
	}
	if (box_sub_points != NULL)
	{
		for (int i = 0; i < subdivisions[0]; i++)
		{
			for (int j = 0; j < subdivisions[1]; j++)
				delete[] box_sub_points[i][j];
			delete[] box_sub_points[i];
		}
		delete[] box_sub_points;
	}
	if (halfedge_lengths != NULL)
	{
		for (int i = 0; i < subdivisions[0]; i++)
		{
			delete[] halfedge_lengths[i];
		}
		delete[]halfedge_lengths;
	}


	if (n_sub_U_knot_vector != NULL)
		delete[] n_sub_U_knot_vector;
	if (n_sub_V_knot_vector != NULL)
		delete[] n_sub_V_knot_vector;
	if (n_sub_weights != NULL)
	{
		for (int i = 0; i < U_dim; i++)
			delete[] n_sub_weights[i];
		delete[]n_sub_weights;
	}
	if (n_sub_control_points != NULL)
	{
		for (int i = 0; i < U_dim; i++)
			delete[] n_sub_control_points[i];
		delete[] n_sub_control_points;
	}
	if (n_sub_Pw != NULL)
	{
		for (int i = 0; i < n_sub_U_dim; i++)
			delete[] n_sub_Pw[i];
		delete[] n_sub_Pw;
	}
	if (n_sub_Bin != NULL)
	{
		for (int i = 0; i < n_sub_U_order + n_sub_V_order + 4; i++)
			delete[] n_sub_Bin[i];
		delete[] n_sub_Bin;
	}
}

//Leitura do arquivo de entrada
bool NURBSSurface::Read(FILE *f)
{
	char s[1000];
	bool CADread;

	fscanf(f, "%s", s);
	number = atoi(s);
	fscanf(f, "%s", s);
	strcpy(file, s);
	//Leitura do arquivo CAD
	CADread = ReadCADFile();
	if (CADread == false)
	{
		printf("Error reading CADData number %d\n", number);
		return false;
	}

	return true;
}

//Escrita do arquivo de saida
void NURBSSurface::Write(FILE *f)
{
	fprintf(f, "NURBSSurface\t%d\t%s\n",
		number,
		file
	);
}

void NURBSSurface::PreCalc()
{
	//Aloca Pw
	Pw = new Matrix*[U_dim];
	for (int i = 0; i < U_dim; i++)
	{
		Pw[i] = new Matrix[V_dim];
		for (int j = 0; j < V_dim; j++)
			Pw[i][j] = Matrix(4);
	}
	//Preenche Pw
	for (int i = 0; i < U_dim; i++)
	{
		for (int j = 0; j < V_dim; j++)
		{
			(Pw[i][j])(0, 0) = weights[i][j] * (control_points[i][j])(0, 0);
			(Pw[i][j])(1, 0) = weights[i][j] * (control_points[i][j])(1, 0);
			(Pw[i][j])(2, 0) = weights[i][j] * (control_points[i][j])(2, 0);
			(Pw[i][j])(3, 0) = weights[i][j];
		}
	}
	Bin = new double*[U_order + V_order + 4];
	for (int i = 0; i < U_order + V_order + 4; i++)
		Bin[i] = new double[U_order + V_order + 4];
	for (int i = 0; i < U_order + V_order + 4; i++)
		for (int j = 0; j < U_order + V_order + 4; j++)
			Bin[i][j] = 0.0;

	//Binomial coefficients
	Bin[0][0] = 1.0;
	for (int i = 1; i < U_order + V_order + 4; i++)
	{
		for (int j = 0; j <= i; j++)
		{
			if (j == 0)
				Bin[i][j] = 1.0;
			else
				Bin[i][j] = Bin[i - 1][j - 1] + Bin[i - 1][j];
		}
	}
	//PrintPtr(Bin, U_order + V_order + 1, U_order + V_order + 1);

	Subdivision();

	ControlPointsSub();

	/* Nao interessa se pensarmos so em NURBSMultipatchSurfcae
	EvaluateRadius();

	if (ReadNurbsData()) {

	}
	else{

		EvaluateVolume();
		EvaluateCentroid();
		EvaluateInertiaTensor();
	}*/
}

//Leitura do arquivo de CAD
bool  NURBSSurface::ReadCADFile()
{
	//Faz leitura do arquivo 'file'
	//Atribui valor para variaveis de NURBS e aloca matrizes de NURBS

	FILE *f1 = NULL;
	char name_file[500];
	strcpy(name_file, db.folder_name);
	strcat(name_file, "CAD/");
	strcat(name_file, file);
	f1 = fopen(name_file, "r");
	if (f1 == NULL)
		return false;//Erro de leitura do arquivo de CAD

	//Leitura do aquivo
	char s[200];

	//Marina
	//Leitura ID
	fscanf(f1, "%s", s);
	if (!strcmp(s, "ID"))
	{
		fscanf(f1, "%s", s);
		id_nurbs = atoi(s);
	}
	else
		return false;
	//Leitura Subdivisions
	fscanf(f1, "%s", s);
	if (!strcmp(s, "Subdivisions"))
	{
		fscanf(f1, "%s", s);
		subdivisions[0] = atoi(s);
		fscanf(f1, "%s", s);
		subdivisions[1] = atoi(s);
	}
	else
		return false;
	//Leitura Connectivity
	fscanf(f1, "%s", s);
	if (!strcmp(s, "Connectivity"))
	{
		fscanf(f1, "%s", s);
		connectivity[0] = atoi(s) - 1;
		fscanf(f1, "%s", s);
		connectivity[1] = atoi(s) - 1;
		fscanf(f1, "%s", s);
		connectivity[2] = atoi(s) - 1;
		fscanf(f1, "%s", s);
		connectivity[3] = atoi(s) - 1;
	}
	else
		return false;
	//Leitura UDim
	fscanf(f1, "%s", s);
	if (!strcmp(s, "UDim"))
	{
		fscanf(f1, "%s", s);
		U_dim = atoi(s);
	}
	else
		return false;
	//Leitura UOrder
	fscanf(f1, "%s", s);
	if (!strcmp(s, "UOrder"))
	{
		fscanf(f1, "%s", s);
		U_order = atoi(s);
	}
	else
		return false;
	//Aloca U knot vector
	U_knot_vector = new double[U_dim + U_order + 1];
	//Leitura U knot vector
	fscanf(f1, "%s", s);
	if (!strcmp(s, "UKnotVector"))
	{
		for (int i = 0; i < (U_dim + U_order + 1); i++)
		{
			fscanf(f1, "%s", s);
			U_knot_vector[i] = atof(s);
		}
	}
	else
		return false;

	//Leitura VDim
	fscanf(f1, "%s", s);
	if (!strcmp(s, "VDim"))
	{
		fscanf(f1, "%s", s);
		V_dim = atoi(s);
	}
	else
		return false;
	//Leitura VOrder
	fscanf(f1, "%s", s);
	if (!strcmp(s, "VOrder"))
	{
		fscanf(f1, "%s", s);
		V_order = atoi(s);
	}
	else
		return false;
	//Aloca V knot vector
	V_knot_vector = new double[V_dim + V_order + 1];
	//Leitura V knot vector
	fscanf(f1, "%s", s);
	if (!strcmp(s, "VKnotVector"))
	{
		for (int i = 0; i < (V_dim + V_order + 1); i++)
		{
			fscanf(f1, "%s", s);
			V_knot_vector[i] = atof(s);
		}
	}
	else
		return false;

	/*
	Vector weights and control points are read assuming they are in a sequence:
	v1 -> u1,...nn, then v2 -> u1,...,un, then ..., vm -> u1,...,un
	with that, if i varies along u and j varies along v, the position related to a given i,j in the vector weights is given at position [i][j]
	The same applies for control points.
	*/
	//Aloca weights
	weights = new double*[U_dim];
	for (int i = 0; i < U_dim; i++)
		weights[i] = new double[V_dim];
	//Aloca control points
	control_points = new Matrix*[U_dim];
	for (int i = 0; i < U_dim; i++)
	{
		control_points[i] = new Matrix[V_dim];
		for (int j = 0; j < V_dim; j++)
			control_points[i][j] = Matrix(3);
	}

	//Leitura weights
	fscanf(f1, "%s", s);
	if (!strcmp(s, "Weights"))
	{
		for (int j = 0; j < V_dim; j++)
		{
			for (int i = 0; i < U_dim; i++)
			{
				fscanf(f1, "%s", s);
				weights[i][j] = atof(s);
			}
		}
	}
	else
		return false;
	//Leitura control points
	fscanf(f1, "%s", s);
	if (!strcmp(s, "ControlPoints"))
	{
		for (int j = 0; j < V_dim; j++)
		{
			for (int i = 0; i < U_dim; i++)
			{
				fscanf(f1, "%s", s);
				control_points[i][j](0, 0) = atof(s);
				fscanf(f1, "%s", s);
				control_points[i][j](1, 0) = atof(s);
				fscanf(f1, "%s", s);
				control_points[i][j](2, 0) = atof(s);
			}
		}
	}
	else
		return false;
	fclose(f1);
	return true;
}

int NURBSSurface::FindSpan(int &n, int &p, double &u, double* U)
{
	/* Determine the knot span index */
	/* Input: n,p,u,U */
	/* Return: the knot span index */
	/*Obs: n is equal to the approximation number of control points -1, because the summation starts on zero.*/

	//Prevention of problems with  span:
	if (u > U[n + 1])
		return n;
	if (u < U[0])
		return p;

	if (u == U[n + 1])
		return(n); /* Special case */
	int low = p;
	int high = n + 1;
	/* Do binary search */
	int mid = (low + high) / 2;
	while (u < U[mid] || u >= U[mid + 1])
	{
		if (u < U[mid])
			high = mid;
		else
			low = mid;
		mid = (low + high) / 2;
	}
	return mid;
}

void NURBSSurface::BasisFunctions(int &i, double &u, int &p, double* U, double* N)
{
	/* Compute the nonvanishing basis functions */
	/* Input: i: just evaluated knot span index
			  u: the value of the coordinate for the evaluation of the basis function
			  p: degree of the polynomial
			  U: complete knot vector
	*/
	/* Output: N - array with dimension p+1 - functions N_(i-p),p,...,N_(i),p*/
	/* OBS: the entire knot vector has to be input. The function makes use of knot span to decide where to evaluate basis function index i)*/
	N[0] = 1.0;
	double* left;
	double* right;
	left = new double[p + 1];
	right = new double[p + 1];
	double temp;
	double saved;

	for (int j = 1; j <= p; j++)
	{
		left[j] = u - U[i + 1 - j];
		right[j] = U[i + j] - u;
		saved = 0.0;
		for (int r = 0; r < j; r++)
		{
			temp = N[r] / (right[r + 1] + left[j - r]);
			N[r] = saved + right[r + 1] * temp;
			saved = left[j - r] * temp;
		}
		N[j] = saved;
	}
	delete[] left;
	delete[] right;
}

void NURBSSurface::DersBasisFunctions(int &i, double &u, int &p, int &n, double* U, double** ders)
{
	/* Compute nonzero basis functions and their derivatives*/
	/*
	Input: i,u,p,n,U
	i: knot span index
	u: knot value (where to evaluate the basis functions and their derivatives
	p: degree of b-splines
	n: order of the derivative required
	U: knot vector
	Output: ders
	ders[k][j], with k'th derivative of the j'th non-null basis function
	*/
	//Allocating auxiliary variables:
	double** ndu;
	double** a;
	ndu = new double*[p + 1];
	for (int i = 0; i < p + 1; i++)
		ndu[i] = new double[p + 1];
	a = new double*[2];
	for (int i = 0; i < 2; i++)
		a[i] = new double[p + 1];
	double* left;
	double* right;
	left = new double[p + 1];
	right = new double[p + 1];
	double temp;
	double saved;
	double d;
	int rk;
	int pk;
	int j1;
	int j2;

	ndu[0][0] = 1.0;
	for (int j = 1; j <= p; j++)
	{
		left[j] = u - U[i + 1 - j];
		right[j] = U[i + j] - u;
		saved = 0.0;
		for (int r = 0; r < j; r++)
		{
			ndu[j][r] = right[r + 1] + left[j - r];
			temp = ndu[r][j - 1] / ndu[j][r];
			ndu[r][j] = saved + right[r + 1] * temp;
			saved = left[j - r] * temp;
		}
		ndu[j][j] = saved;
	}
	/* Load the basis functions */
	for (int j = 0; j <= p; j++)
		ders[0][j] = ndu[j][p];
	/* This section computes the derivatives (Eq. [2.9]) */
	for (int r = 0; r <= p; r++) /* Loop over function index */
	{
		int s1 = 0;
		int s2 = 1; /* Alternate rows in array a */
		a[0][0] = 1.0;
		/* Loop to compute kth derivative */
		for (int k = 1; k <= n; k++)
		{
			d = 0.0;
			rk = r - k;
			pk = p - k;

			if (r >= k)
			{
				a[s2][0] = a[s1][0] / ndu[pk + 1][rk];
				d = a[s2][0] * ndu[rk][pk];
			}
			if (rk >= -1)
				j1 = 1;
			else
				j1 = -rk;
			if (r - 1 <= pk)
				j2 = k - 1;
			else
				j2 = p - r;
			for (int j = j1; j <= j2; j++)
			{
				a[s2][j] = (a[s1][j] - a[s1][j - 1]) / ndu[pk + 1][rk + j];
				d += a[s2][j] * ndu[rk + j][pk];
			}
			if (r <= pk)
			{
				a[s2][k] = -a[s1][k - 1] / ndu[pk + 1][r];
				d += a[s2][k] * ndu[r][pk];
			}
			ders[k][r] = d;
			int j = s1;
			s1 = s2;
			s2 = j; /* Switch rows */
		}
	}
	/* Multiply through by the correct factors */
	/* (Eq. [2.9]) */
	int r = p;
	for (int k = 1; k <= n; k++)
	{
		for (int j = 0; j <= p; j++)
			ders[k][j] *= r;
		r *= (p - k);
	}

	//De-allocating auxiliary variables
	for (int i = 0; i < p + 1; i++)
		delete[]ndu[i];
	delete[]ndu;
	for (int i = 0; i < 2; i++)
		delete[]a[i];
	delete[]a;
	delete[] left;
	delete[] right;
}

void NURBSSurface::NURBSPoint(double &uc, double &vc, Matrix &point)
{
	/*Input: u,v - knot vector coordinates
	  Output: point - vector indicating the location of the point associated with (u,v)
	 */
	int n;
	//knot span index u
	n = U_dim - 1;
	int span_u = FindSpan(n, U_order, uc, U_knot_vector);
	//knot span index v
	n = V_dim - 1;
	int span_v = FindSpan(n, V_order, vc, V_knot_vector);
	//Non-vanishing functions - u
	double* N_u = new double[U_order + 1];
	BasisFunctions(span_u, uc, U_order, U_knot_vector, N_u);
	//Non-vanishing functions - v
	double* N_v = new double[V_order + 1];
	BasisFunctions(span_v, vc, V_order, V_knot_vector, N_v);
	//Evaluating the surface point
	Matrix* temp;
	temp = new Matrix[V_order + 1];
	for (int l = 0; l <= V_order; l++)
	{
		temp[l] = Matrix(4);
		for (int k = 0; k <= U_order; k++)
			temp[l] = temp[l] + N_u[k] * Pw[span_u - U_order + k][span_v - V_order + l];
	}
	Matrix Sw(4);
	for (int l = 0; l <= V_order; l++)
		Sw = Sw + N_v[l] * temp[l];
	//Computing the surface point from Sw:
	point(0, 0) = Sw(0, 0) / Sw(3, 0);
	point(1, 0) = Sw(1, 0) / Sw(3, 0);
	point(2, 0) = Sw(2, 0) / Sw(3, 0);

	delete[] temp;
	delete[] N_u;
	delete[] N_v;
}

void NURBSSurface::NURBSDerivatives(double &uc, double &vc, Matrix** &Skl, int &d)
{
	/*Input: u,v - knot vector coordinates
			 d - order of maximum derivative 0<=(k+l)<=d
	  Output: Collection of vectors Skl[k][l] with derivatives k'th in u and l'th in v, evaluated at the location of the point associated with (u,v)
	 */
	int n;
	//knot span index u
	n = U_dim - 1;
	int span_u = FindSpan(n, U_order, uc, U_knot_vector);
	//knot span index v
	n = V_dim - 1;
	int span_v = FindSpan(n, V_order, vc, V_knot_vector);
	//Evaluating orders of derivatives
	int du = min(U_order, d);
	int dv = min(V_order, d);
	//Non-vanishing functions and derivatives - u: Nu
	double** Nu = new double*[du + 1];
	for (int i = 0; i < du + 1; i++)
		Nu[i] = new double[U_order + 1];
	DersBasisFunctions(span_u, uc, U_order, du, U_knot_vector, Nu);
	//Non-vanishing functions and derivatives - v: Nv
	double** Nv = new double*[dv + 1];
	for (int i = 0; i < dv + 1; i++)
		Nv[i] = new double[V_order + 1];
	DersBasisFunctions(span_v, vc, V_order, dv, V_knot_vector, Nv);

	Matrix* temp;
	temp = new Matrix[V_order + 1];

	Matrix** Skl1;
	Skl1 = new Matrix*[d + 1];
	for (int i = 0; i < d + 1; i++)
		Skl1[i] = new Matrix[d + 1];
	for (int i = 0; i < d + 1; i++)
		for (int j = 0; j < d + 1; j++)
			Skl1[i][j] = Matrix(4);

	//Algorithm 3.6 (pag 111) - evaluating at Skl1 - 4D (Sw)
	for (int k = 0; k <= du; k++)
	{
		{
			for (int s = 0; s <= V_order; s++)
			{
				temp[s] = Matrix(4);
				for (int r = 0; r <= U_order; r++)
					temp[s] = temp[s] + Nu[k][r] * Pw[span_u - U_order + r][span_v - V_order + s];
			}
			int dd = min(d - k, dv);
			for (int l = 0; l <= dd; l++)
			{
				Skl1[k][l] = Matrix(4);
				for (int s = 0; s <= V_order; s++)
					Skl1[k][l] = Skl1[k][l] + Nv[l][s] * temp[s];
			}
		}
	}
	//Algorithm 4.4 (pag 137) - evaluating Skl (S) from Skl1 (Sw)
	Matrix v(3);

	//Cleans Skl
	for (int i = 0; i < d + 1; i++)
	{
		for (int j = 0; j < d + 1; j++)
		{
			Skl[i][j](0, 0) = 0.0;
			Skl[i][j](1, 0) = 0.0;
			Skl[i][j](2, 0) = 0.0;
		}
	}

	for (int k = 0; k <= d; k++)
	{
		for (int l = 0; l <= d - k; l++)
		{
			v(0, 0) = Skl1[k][l](0, 0);
			v(1, 0) = Skl1[k][l](1, 0);
			v(2, 0) = Skl1[k][l](2, 0);
			for (int j = 1; j <= l; j++)
				v = v - Bin[l][j] * Skl1[0][j](3, 0) * Skl[k][l - j];
			for (int i = 1; i <= k; i++)
			{
				v = v - Bin[k][i] * Skl1[i][0](3, 0) * Skl[k - i][l];
				Matrix v2(3);
				for (int j = 1; j <= l; j++)
					v2 = v2 + Bin[l][j] * Skl1[i][j](3, 0) * Skl[k - i][l - j];
				v = v - Bin[k][i] * v2;
			}
			Skl[k][l] = 1.0 / Skl1[0][0](3, 0) * v;
		}
	}
	//Cleaning
	for (int i = 0; i < du+1; i++)
		delete[] Nu[i];
	delete[] Nu;
	for (int i = 0; i < dv+1; i++)
		delete[] Nv[i];
	delete[] Nv;
	delete[] temp;
	for (int i = 0; i < d+1; i++)
		delete[] Skl1[i];
	delete[] Skl1;
}

void NURBSSurface::WriteVTK_XMLRender(FILE *f, Matrix& pos, Matrix& rot, int number)
{
	//vetores para escrita no formato binario - usando a funcao 'enconde'
	std::vector<float> float_vector;
	std::vector<int> int_vector;

	int intervals_U = U_dim - U_order;
	int intervals_V = V_dim - V_order;
	int nU = U_order;
	int nV = V_order;
	int n_div_u = intervals_U * nU;
	int n_div_v = intervals_V * nV;

	int n_points = (2 * n_div_u + 1)*(2 * n_div_v + 1);
	int n_cells = n_div_u * n_div_v;
	
	double u, v;
	double du, dv;
	//Opens Piece
	fprintf(f, "\t\t<Piece NumberOfPoints = \"%d\" NumberOfCells = \"%d\">\n", n_points, n_cells);
	//Opens Points
	fprintf(f, "\t\t\t<Points>\n");
	//Opens DataArray
	fprintf(f, "\t\t\t\t<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"binary\">\n");
	float_vector.clear();
	
	Matrix local(3);
	for (int i = U_order; i < U_dim; i++)
	{
		du = (U_knot_vector[i + 1] - U_knot_vector[i]) / (2 * nU);
		for (int ii = 0; ii < 2 * nU; ii++)
		{
			u = U_knot_vector[i] + ii * du;
			for (int j = V_order; j < V_dim; j++)
			{
				dv = (V_knot_vector[j + 1] - V_knot_vector[j]) / (2 * nV);
				for (int jj = 0; jj < 2 * nV; jj++)
				{
					v = V_knot_vector[j] + jj * dv;
					NURBSPoint(u, v, local);
					local = rot * local;
					float_vector.push_back((float)(pos(0, 0) + (local)(0, 0)));
					float_vector.push_back((float)(pos(1, 0) + (local)(1, 0)));
					float_vector.push_back((float)(pos(2, 0) + (local)(2, 0)));
				}
			}
			v = V_knot_vector[V_dim];
			NURBSPoint(u, v, local);
			local = rot * local;
			float_vector.push_back((float)(pos(0, 0) + (local)(0, 0)));
			float_vector.push_back((float)(pos(1, 0) + (local)(1, 0)));
			float_vector.push_back((float)(pos(2, 0) + (local)(2, 0)));
		}
	}
	u = U_knot_vector[U_dim];
	for (int j = V_order; j < V_dim; j++)
	{
		dv = (V_knot_vector[j + 1] - V_knot_vector[j]) / (2 * nV);
		for (int jj = 0; jj < 2 * nV; jj++)
		{
			v = V_knot_vector[j] + jj * dv;
			NURBSPoint(u, v, local);
			local = rot * local;
			float_vector.push_back((float)(pos(0, 0) + (local)(0, 0)));
			float_vector.push_back((float)(pos(1, 0) + (local)(1, 0)));
			float_vector.push_back((float)(pos(2, 0) + (local)(2, 0)));
		}
	}
	v = V_knot_vector[V_dim];
	NURBSPoint(u, v, local);
	local = rot * local;
	float_vector.push_back((float)(pos(0, 0) + (local)(0, 0)));
	float_vector.push_back((float)(pos(1, 0) + (local)(1, 0)));
	float_vector.push_back((float)(pos(2, 0) + (local)(2, 0)));


	fprintf(f, encodeData<float>(float_vector).c_str());
	fprintf(f, "\n");
	//Closes DataArray
	fprintf(f, "\t\t\t\t</DataArray>\n");
	//Closes Points
	fprintf(f, "\t\t\t</Points>\n");

	//Opens Cells
	fprintf(f, "\t\t\t<Cells>\n");
	//Opens DataArray
	fprintf(f, "\t\t\t\t<DataArray type=\"Int32\" Name=\"connectivity\" format=\"binary\">\n");
	//Linhas
	int nj = (2 * n_div_v + 1);
	int_vector.clear();
	for (int i = 0; i < (2 * n_div_u - 1); i = i + 2)
		for (int j = 0; j < (2 * n_div_v - 1); j = j + 2)
		{
			int_vector.push_back(j + 0 + (i + 0)*nj);
			int_vector.push_back(j + 0 + (i + 2)*nj);
			int_vector.push_back(j + 2 + (i + 2)*nj);
			int_vector.push_back(j + 2 + (i + 0)*nj);
			int_vector.push_back(j + 0 + (i + 1)*nj);
			int_vector.push_back(j + 1 + (i + 2)*nj);
			int_vector.push_back(j + 2 + (i + 1)*nj);
			int_vector.push_back(j + 1 + (i + 0)*nj);
		}
	fprintf(f, encodeData(int_vector).c_str());
	fprintf(f, "\n");
	//Closes DataArray
	fprintf(f, "\t\t\t\t</DataArray>\n");
	//Opens DataArray
	fprintf(f, "\t\t\t\t<DataArray type=\"Int32\" Name=\"types\" format=\"binary\">\n");
	int_vector.clear();
	for (int i = 0; i < n_cells; i++)
		int_vector.push_back(23);
	fprintf(f, encodeData(int_vector).c_str());
	fprintf(f, "\n");
	//Closes DataArray
	fprintf(f, "\t\t\t\t</DataArray>\n");
	//Opens DataArray
	fprintf(f, "\t\t\t\t<DataArray type=\"Int32\" Name=\"offsets\" format=\"binary\">\n");
	int_vector.clear();
	for (int i = 0; i < n_cells; i++)
		int_vector.push_back((i + 1) * 8);
	fprintf(f, encodeData(int_vector).c_str());
	fprintf(f, "\n");
	//Closes DataArray
	fprintf(f, "\t\t\t\t</DataArray>\n");
	//Closes Cells
	fprintf(f, "\t\t\t</Cells>\n");

	//Opens CellData
	fprintf(f, "\t\t\t<CellData FieldData=\"SurfaceData\">\n");
	int_vector.clear();
	//Opens DataArray
	fprintf(f, "\t\t\t\t<DataArray Name=\"SurfaceProperties\" type=\"Int32\" NumberOfComponents=\"%d\" format=\"binary\">\n", 1);
	for (int cell = 0; cell < n_cells; cell++)
	{
		int_vector.push_back(number);		//Surface ID
	}
	fprintf(f, encodeData(int_vector).c_str());
	fprintf(f, "\n");
	//Closes DataArray
	fprintf(f, "\t\t\t\t</DataArray>\n");
	//Closes CellData
	fprintf(f, "\t\t\t</CellData>\n");

	//Closes Piece
	fprintf(f, "\t\t</Piece>\n");
}

void NURBSSurface::WriteVTK_XMLRender(FILE *f, Matrix& pos, Matrix& rot, int number, int number2)
{
	//vetores para escrita no formato binario - usando a funcao 'enconde'
	std::vector<float> float_vector;
	std::vector<int> int_vector;

	int intervals_U = U_dim - U_order;
	int intervals_V = V_dim - V_order;
	int nU = U_order;
	int nV = V_order;
	int n_div_u = intervals_U * nU;
	int n_div_v = intervals_V * nV;

	int n_points = (2 * n_div_u + 1)*(2 * n_div_v + 1);
	int n_cells = n_div_u * n_div_v;

	double u, v;
	double du, dv;
	//Opens Piece
	fprintf(f, "\t\t<Piece NumberOfPoints = \"%d\" NumberOfCells = \"%d\">\n", n_points, n_cells);
	//Opens Points
	fprintf(f, "\t\t\t<Points>\n");
	//Opens DataArray
	fprintf(f, "\t\t\t\t<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"binary\">\n");
	float_vector.clear();

	Matrix local(3);
	for (int i = U_order; i < U_dim; i++)
	{
		du = (U_knot_vector[i + 1] - U_knot_vector[i]) / (2 * nU);
		for (int ii = 0; ii < 2 * nU; ii++)
		{
			u = U_knot_vector[i] + ii * du;
			for (int j = V_order; j < V_dim; j++)
			{
				dv = (V_knot_vector[j + 1] - V_knot_vector[j]) / (2 * nV);
				for (int jj = 0; jj < 2 * nV; jj++)
				{
					v = V_knot_vector[j] + jj * dv;
					NURBSPoint(u, v, local);
					local = rot * local;
					float_vector.push_back((float)(pos(0, 0) + (local)(0, 0)));
					float_vector.push_back((float)(pos(1, 0) + (local)(1, 0)));
					float_vector.push_back((float)(pos(2, 0) + (local)(2, 0)));
				}
			}
			v = V_knot_vector[V_dim];
			NURBSPoint(u, v, local);
			local = rot * local;
			float_vector.push_back((float)(pos(0, 0) + (local)(0, 0)));
			float_vector.push_back((float)(pos(1, 0) + (local)(1, 0)));
			float_vector.push_back((float)(pos(2, 0) + (local)(2, 0)));
		}
	}
	u = U_knot_vector[U_dim];
	for (int j = V_order; j < V_dim; j++)
	{
		dv = (V_knot_vector[j + 1] - V_knot_vector[j]) / (2 * nV);
		for (int jj = 0; jj < 2 * nV; jj++)
		{
			v = V_knot_vector[j] + jj * dv;
			NURBSPoint(u, v, local);
			local = rot * local;
			float_vector.push_back((float)(pos(0, 0) + (local)(0, 0)));
			float_vector.push_back((float)(pos(1, 0) + (local)(1, 0)));
			float_vector.push_back((float)(pos(2, 0) + (local)(2, 0)));
		}
	}
	v = V_knot_vector[V_dim];
	NURBSPoint(u, v, local);
	local = rot * local;
	float_vector.push_back((float)(pos(0, 0) + (local)(0, 0)));
	float_vector.push_back((float)(pos(1, 0) + (local)(1, 0)));
	float_vector.push_back((float)(pos(2, 0) + (local)(2, 0)));


	fprintf(f, encodeData<float>(float_vector).c_str());
	fprintf(f, "\n");
	//Closes DataArray
	fprintf(f, "\t\t\t\t</DataArray>\n");
	//Closes Points
	fprintf(f, "\t\t\t</Points>\n");

	//Opens Cells
	fprintf(f, "\t\t\t<Cells>\n");
	//Opens DataArray
	fprintf(f, "\t\t\t\t<DataArray type=\"Int32\" Name=\"connectivity\" format=\"binary\">\n");
	//Linhas
	int nj = (2 * n_div_v + 1);
	int_vector.clear();
	for (int i = 0; i < (2 * n_div_u - 1); i = i + 2)
		for (int j = 0; j < (2 * n_div_v - 1); j = j + 2)
		{
			int_vector.push_back(j + 0 + (i + 0)*nj);
			int_vector.push_back(j + 0 + (i + 2)*nj);
			int_vector.push_back(j + 2 + (i + 2)*nj);
			int_vector.push_back(j + 2 + (i + 0)*nj);
			int_vector.push_back(j + 0 + (i + 1)*nj);
			int_vector.push_back(j + 1 + (i + 2)*nj);
			int_vector.push_back(j + 2 + (i + 1)*nj);
			int_vector.push_back(j + 1 + (i + 0)*nj);
		}
	fprintf(f, encodeData(int_vector).c_str());
	fprintf(f, "\n");
	//Closes DataArray
	fprintf(f, "\t\t\t\t</DataArray>\n");
	//Opens DataArray
	fprintf(f, "\t\t\t\t<DataArray type=\"Int32\" Name=\"types\" format=\"binary\">\n");
	int_vector.clear();
	for (int i = 0; i < n_cells; i++)
		int_vector.push_back(23);
	fprintf(f, encodeData(int_vector).c_str());
	fprintf(f, "\n");
	//Closes DataArray
	fprintf(f, "\t\t\t\t</DataArray>\n");
	//Opens DataArray
	fprintf(f, "\t\t\t\t<DataArray type=\"Int32\" Name=\"offsets\" format=\"binary\">\n");
	int_vector.clear();
	for (int i = 0; i < n_cells; i++)
		int_vector.push_back((i + 1) * 8);
	fprintf(f, encodeData(int_vector).c_str());
	fprintf(f, "\n");
	//Closes DataArray
	fprintf(f, "\t\t\t\t</DataArray>\n");
	//Closes Cells
	fprintf(f, "\t\t\t</Cells>\n");

	//Opens CellData
	fprintf(f, "\t\t\t<CellData FieldData=\"SurfaceData\">\n");
	int_vector.clear();
	//Opens DataArray
	fprintf(f, "\t\t\t\t<DataArray Name=\"SurfaceProperties\" type=\"Int32\" NumberOfComponents=\"%d\" format=\"binary\">\n", 1);
	for (int cell = 0; cell < n_cells; cell++)
	{
		int_vector.push_back(number2);		//CADData ID
	}
	fprintf(f, encodeData(int_vector).c_str());
	fprintf(f, "\n");
	//Closes DataArray
	fprintf(f, "\t\t\t\t</DataArray>\n");
	int_vector.clear();
	//Opens DataArray
	fprintf(f, "\t\t\t\t<DataArray Name=\"CAD\" type=\"Int32\" NumberOfComponents=\"%d\" format=\"binary\">\n", 1);
	for (int cell = 0; cell < n_cells; cell++)
	{
		int_vector.push_back(number);		//CADData ID
	}
	fprintf(f, encodeData(int_vector).c_str());
	fprintf(f, "\n");
	//Closes DataArray
	fprintf(f, "\t\t\t\t</DataArray>\n");
	//Closes CellData
	fprintf(f, "\t\t\t</CellData>\n");

	//Closes Piece
	fprintf(f, "\t\t</Piece>\n");
}

void NURBSSurface::SurfaceKnotIns(int np, int p, double * UP, int rp, int mp, int q, double * VP, int sp, Matrix ** Pw, int dir, double uv, int k, int s, int r, int & nq, double * UQ, int & mq, double * VQ, Matrix ** Qw)
{
	/* A5.3 */
	/* Surface knot insertion */
	/* Input: np, p, UP, rp, mp, q, VP, sp, Pw, dir, uv, k, s, r */
	/* Output: nq, UQ, rq, mq, VQ, sq, Qw */

	if (dir == 1) { //u-direction
		// Load new knot vector 
		for (int i = 0; i <= k; i++) UQ[i] = UP[i];
		for (int i = 1; i <= r; i++) UQ[k + i] = uv;
		for (int i = k + 1; i <= rp; i++) UQ[i + r] = UP[i];
		// Copy v-vector into VQ 
		for (int i = 0; i < sp; i++)
		{
			VQ[i] = VP[i];
		}
		// Save the alphas
		int L = 0;
		double **alpha;
		alpha = new double*[p - 1 - s + 1];
		for (int i = 0; i < (p - 1 - s + 1); i++)
		{
			alpha[i] = new double[r];
		}
		for (int j = 1; j <= r; j++)
		{
			L = k - p + j;
			for (int i = 0; i <= p - j - s; i++)
				alpha[i][j] = (uv - UP[L + i]) / (UP[i + k + 1] - UP[L + i]);
		}
		// For each row do
		double **Rw;
		Rw = new double*[4];
		for (int i = 0; i < 4; i++)
		{
			Rw[i] = new double[p - s + 1];
		}
		for (int row = 0; row <= (mp - 1); row++)
		{
			// Save unaltered control points
			for (int i = 0; i <= k - p; i++)
			{
				Qw[i][row](0, 0) = Pw[i][row](0, 0);
				Qw[i][row](1, 0) = Pw[i][row](1, 0);
				Qw[i][row](2, 0) = Pw[i][row](2, 0);
				Qw[i][row](3, 0) = Pw[i][row](3, 0);
			}
			for (int i = k - s; i <= (np - 1); i++)
			{
				Qw[i + r][row](0, 0) = Pw[i][row](0, 0);
				Qw[i + r][row](1, 0) = Pw[i][row](1, 0);
				Qw[i + r][row](2, 0) = Pw[i][row](2, 0);
				Qw[i + r][row](3, 0) = Pw[i][row](3, 0);
			}
			// Load auxiliary control points
			for (int i = 0; i <= p - s; i++)
			{
				Rw[0][i] = Pw[k - p + i][row](0, 0);
				Rw[1][i] = Pw[k - p + i][row](1, 0);
				Rw[2][i] = Pw[k - p + i][row](2, 0);
				Rw[3][i] = Pw[k - p + i][row](3, 0);
			}
			// Insert the knot r times
			for (int j = 1; j <= r; j++)
			{
				L = k - p + j;
				for (int i = 0; i <= p - j - s; i++)
				{
					Rw[0][i] = alpha[i][j] * Rw[0][i + 1] + (1.0 - alpha[i][j])*Rw[0][i];
					Rw[1][i] = alpha[i][j] * Rw[1][i + 1] + (1.0 - alpha[i][j])*Rw[1][i];
					Rw[2][i] = alpha[i][j] * Rw[2][i + 1] + (1.0 - alpha[i][j])*Rw[2][i];
					Rw[3][i] = alpha[i][j] * Rw[3][i + 1] + (1.0 - alpha[i][j])*Rw[3][i];
				}
				Qw[L][row](0, 0) = Rw[0][0];
				Qw[L][row](1, 0) = Rw[1][0];
				Qw[L][row](2, 0) = Rw[2][0];
				Qw[L][row](3, 0) = Rw[3][0];
				Qw[k + r - j - s][row](0, 0) = Rw[0][p - j - s];
				Qw[k + r - j - s][row](1, 0) = Rw[1][p - j - s];
				Qw[k + r - j - s][row](2, 0) = Rw[2][p - j - s];
				Qw[k + r - j - s][row](3, 0) = Rw[3][p - j - s];
			}
			// Load the remaining control points
			for (int i = L + 1; i < k - s; i++) {
				Qw[i][row](0, 0) = Rw[0][i - L];
				Qw[i][row](1, 0) = Rw[1][i - L];
				Qw[i][row](2, 0) = Rw[2][i - L];
				Qw[i][row](3, 0) = Rw[3][i - L];
			}
		}
		nq = np + r;
	}
	else if (dir == 2) { //v-direction
		// Load new knot vector 
		for (int i = 0; i <= k; i++) VQ[i] = VP[i];
		for (int i = 1; i <= r; i++) VQ[k + i] = uv;
		for (int i = k + 1; i <= sp; i++) VQ[i + r] = VP[i];
		// Copy u-vector into UQ 
		for (int i = 0; i < rp; i++)
		{
			UQ[i] = UP[i];
		}
		// Save the alphas
		int L = 0;
		double **alpha;
		alpha = new double*[q - 1 - s + 1];
		for (int i = 0; i < (q - 1 - s + 1); i++)
		{
			alpha[i] = new double[r];
		}
		for (int j = 1; j <= r; j++)
		{
			L = k - q + j;
			for (int i = 0; i <= q - j - s; i++)
				alpha[i][j] = (uv - VP[L + i]) / (VP[i + k + 1] - VP[L + i]);
		}
		// For each row do
		double **Rw;
		Rw = new double*[4];
		for (int i = 0; i < 4; i++)
		{
			Rw[i] = new double[q - s + 1];
		}
		for (int row = 0; row <= (np - 1); row++)
		{
			// Save unaltered control points
			for (int i = 0; i <= k - q; i++)
			{
				Qw[row][i](0, 0) = Pw[row][i](0, 0);
				Qw[row][i](1, 0) = Pw[row][i](1, 0);
				Qw[row][i](2, 0) = Pw[row][i](2, 0);
				Qw[row][i](3, 0) = Pw[row][i](3, 0);
			}
			for (int i = k - s; i <= (mp - 1); i++)
			{
				Qw[row][i + r](0, 0) = Pw[row][i](0, 0);
				Qw[row][i + r](1, 0) = Pw[row][i](1, 0);
				Qw[row][i + r](2, 0) = Pw[row][i](2, 0);
				Qw[row][i + r](3, 0) = Pw[row][i](3, 0);
			}
			// Load auxiliary control points
			for (int i = 0; i <= q - s; i++)
			{
				Rw[0][i] = Pw[row][k - q + i](0, 0);
				Rw[1][i] = Pw[row][k - q + i](1, 0);
				Rw[2][i] = Pw[row][k - q + i](2, 0);
				Rw[3][i] = Pw[row][k - q + i](3, 0);
			}
			// Insert the knot r times
			for (int j = 1; j <= r; j++)
			{
				L = k - q + j;
				for (int i = 0; i <= q - j - s; i++)
				{
					Rw[0][i] = alpha[i][j] * Rw[0][i + 1] + (1.0 - alpha[i][j])*Rw[0][i];
					Rw[1][i] = alpha[i][j] * Rw[1][i + 1] + (1.0 - alpha[i][j])*Rw[1][i];
					Rw[2][i] = alpha[i][j] * Rw[2][i + 1] + (1.0 - alpha[i][j])*Rw[2][i];
					Rw[3][i] = alpha[i][j] * Rw[3][i + 1] + (1.0 - alpha[i][j])*Rw[3][i];
				}
				Qw[row][L](0, 0) = Rw[0][0];
				Qw[row][L](1, 0) = Rw[1][0];
				Qw[row][L](2, 0) = Rw[2][0];
				Qw[row][L](3, 0) = Rw[3][0];
				Qw[row][k + r - j - s](0, 0) = Rw[0][q - j - s];
				Qw[row][k + r - j - s](1, 0) = Rw[1][q - j - s];
				Qw[row][k + r - j - s](2, 0) = Rw[2][q - j - s];
				Qw[row][k + r - j - s](3, 0) = Rw[3][q - j - s];
			}
			// Load the remaining control points
			for (int i = L + 1; i < k - s; i++) {
				Qw[row][i](0, 0) = Rw[0][i - L];
				Qw[row][i](1, 0) = Rw[1][i - L];
				Qw[row][i](2, 0) = Rw[2][i - L];
				Qw[row][i](3, 0) = Rw[3][i - L];
			}
		}
		mq = mp + r;
	}
}

void NURBSSurface::PointInversion(Matrix p, double &uf, double &vf, Matrix p_f, bool &conv, double &tol, double &err, int uv)
{

	double uvf = 0.0;
	if (uv == 0) {
		uvf = uf;
	}
	else {
		uvf = vf;
	}

	//Numero de nos da NURBS
	int r = (U_dim - 1) + U_order + 1;
	int s = (V_dim - 1) + V_order + 1;

	//Derivadas da superficie no ponto inicial
	int d = 2;
	Matrix **Skl2;
	Skl2 = new Matrix*[d + 1];
	for (int i = 0; i < d + 1; i++)
	{
		Skl2[i] = new Matrix[d + 1];
		for (int j = 0; j < d + 1; j++)
			Skl2[i][j] = Matrix(3);
	}
	NURBSDerivatives(uf, vf, Skl2, d);

	double delta = 0.0;

	double delta_u = 0.0; //Passo em u
	double delta_v = 0.0; //Passo em v

	double step_n = sqrt(delta_u * delta_u + delta_v * delta_v); //Norma do passo

	double delta_u_b = 0.0;
	double delta_v_b = 0.0;

	double delta_u_u = 0.0;
	double delta_v_u = 0.0;

	double uf_test = 0.0;

	double vf_test = 0.0;

	//Gradiente no ponto inicial
	double g[2];
	g[0] = -(Skl2[1][0](0, 0) * (Skl2[0][0](0, 0) - p(0, 0)) + Skl2[1][0](1, 0) * (Skl2[0][0](1, 0) - p(1, 0)) + Skl2[1][0](2, 0) * (Skl2[0][0](2, 0) - p(2, 0)));
	g[1] = -(Skl2[0][1](0, 0) * (Skl2[0][0](0, 0) - p(0, 0)) + Skl2[0][1](1, 0) * (Skl2[0][0](1, 0) - p(1, 0)) + Skl2[0][1](2, 0) * (Skl2[0][0](2, 0) - p(2, 0)));

	double g_abs = sqrt(g[0] * g[0] + g[1] * g[1]);
	err = g_abs;
	//if (err == 0) {
		//cout << "wrong" << endl;
	//}

	//Hessiana no ponto inicial
	Matrix H(2, 2);
	H(0, 0) = Skl2[2][0](0, 0) * Skl2[0][0](0, 0) + Skl2[1][0](0, 0) * Skl2[1][0](0, 0) - Skl2[2][0](0, 0) * p(0, 0)
		+ Skl2[2][0](1, 0) * Skl2[0][0](1, 0) + Skl2[1][0](1, 0) * Skl2[1][0](1, 0) - Skl2[2][0](1, 0) * p(1, 0)
		+ Skl2[2][0](2, 0) * Skl2[0][0](2, 0) + Skl2[1][0](2, 0) * Skl2[1][0](2, 0) - Skl2[2][0](2, 0) * p(2, 0);

	H(0, 1) = Skl2[1][1](0, 0) * Skl2[0][0](0, 0) + Skl2[1][0](0, 0) * Skl2[0][1](0, 0) - Skl2[1][1](0, 0) * p(0, 0)
		+ Skl2[1][1](1, 0) * Skl2[0][0](1, 0) + Skl2[1][0](1, 0) * Skl2[0][1](1, 0) - Skl2[1][1](1, 0) * p(1, 0)
		+ Skl2[1][1](2, 0) * Skl2[0][0](2, 0) + Skl2[1][0](2, 0) * Skl2[0][1](2, 0) - Skl2[1][1](2, 0) * p(2, 0);

	H(1, 0) = H(0, 1);

	H(1, 1) = Skl2[0][2](0, 0) * Skl2[0][0](0, 0) + Skl2[0][1](0, 0) * Skl2[0][1](0, 0) - Skl2[0][2](0, 0) * p(0, 0)
		+ Skl2[0][2](1, 0) * Skl2[0][0](1, 0) + Skl2[0][1](1, 0) * Skl2[0][1](1, 0) - Skl2[0][2](1, 0) *p(1, 0)
		+ Skl2[0][2](2, 0) * Skl2[0][0](2, 0) + Skl2[0][1](2, 0) * Skl2[0][1](2, 0) - Skl2[0][2](2, 0) * p(2, 0);

	//Autovalores da Hessiana no ponto inicial
	Matrix P(2, 2);
	Matrix D(2, 2);
	double tol_eig = 1e-14;
	fulleigen1(H, P, D, tol_eig);
	double max_eig = D(0, 0);
	double min_eig = D(1, 1);
	int index = 1;
	if (D(1, 1) > max_eig) {
		min_eig = D(0, 0);
		max_eig = D(1, 1);
		index = 0;
	}

	double error = pow(10, -12) * abs(max_eig); //Tolerancia para o erro no metodo de Newton

	tol = error;

	//Variaveis - TR

	double delta_max = sqrt((U_knot_vector[r - U_order] - U_knot_vector[U_order]) * (U_knot_vector[r - U_order] - U_knot_vector[U_order]) / 100 + (V_knot_vector[s - V_order] - V_knot_vector[V_order]) * (V_knot_vector[s - V_order] - V_knot_vector[V_order]) / 100);

	if (delta == 0) {
		delta = delta_max / 10;
	}

	double eta = 0.25;

	double rt = 0.0;

	double obj1 = 0.0;
	double obj2 = 0.0;

	double quad1 = 0.0;
	double quad2 = 0.0;

	double cp_u = 0.0;
	double cp_v = 0.0;

	double tau = 0.0;
	double gbg = 0.0;

	//Minimizacao em subespaco bidimensional - variaveis

	double A[2][2];
	A[0][0] = 0.0;
	A[0][1] = 0.0;
	A[1][0] = 0.0;
	A[1][1] = 0.0;

	double M[2][2];
	M[0][0] = 0.0;
	M[0][1] = 0.0;
	M[1][0] = 0.0;
	M[1][1] = 0.0;

	double w[2];
	w[0] = 0.0;
	w[1] = 0.0;

	double w_n = sqrt(w[0] * w[0] + w[1] * w[1]);

	double Uo[2];
	Uo[0] = 0.0;
	Uo[1] = 0.0;

	double G[2][2];
	G[0][0] = 0.0;
	G[0][1] = 0.0;
	G[1][0] = 0.0;
	G[1][1] = 0.0;

	double c1 = 0.0;
	double c2 = 0.0;
	double c3 = 0.0;
	double c4 = 0.0;
	double c5 = 0.0;

	double l = 0.0;
	double o = 0.0;

	double lf = 0.0;
	double of = 0.0;

	double quadr1 = 0.0;

	double quadr2 = 0.0;

	double alpha = 0.0;
	double gamma = 0.0;

	double qo[2];
	qo[0] = 0.0;
	qo[1] = 0.0;

	bool first = true;

	//Polinomio 2o grau

	double a_pol = 0.0;
	double b_pol = 0.0;
	double c_pol = 0.0;
	double x1_pol = 0.0;
	double x2_pol = 0.0;

	int it = 1; //Iteracao

	//Newton com TR
	while (g_abs > error || step_n > error)
	{

		//Se a matriz hessiana e positivo definida
		if (min_eig > tol_eig)
		{
			delta_v_b = (g[1] - g[0] * H(1, 0) / H(0, 0)) / (H(1, 1) - H(0, 1) * H(1, 0) / H(0, 0));
			delta_u_b = (g[0] - H(0, 1) * delta_v_b) / H(0, 0);
			step_n = sqrt(delta_u_b * delta_u_b + delta_v_b * delta_v_b);

			//O passo esta dentro da regiao de confianca?
			//Sim
			if (step_n <= delta) {
				delta_u = delta_u_b;
				delta_v = delta_v_b;
				step_n = sqrt(delta_u * delta_u + delta_v * delta_v);
			}
			//Nao
			else {
				//Dogleg method

				gbg = g[0] * (H(0, 0) * g[0] + H(0, 1) * g[1]) + g[1] * (H(1, 0) * g[0] + H(1, 1) * g[1]);

				delta_u_u = pow(g_abs, 2) / gbg * g[0];

				delta_v_u = pow(g_abs, 2) / gbg * g[1];

				a_pol = (pow((delta_u_b - delta_u_u), 2) + pow((delta_v_b - delta_v_u), 2));
				if (a_pol < pow(10, -15)) {
					a_pol = 0;
				}
				b_pol = (2 * delta_u_u * (delta_u_b - delta_u_u) - 2 * pow((delta_u_b - delta_u_u), 2) + 2 * delta_v_u * (delta_v_b - delta_v_u) - 2 * pow((delta_v_b - delta_v_u), 2));
				c_pol = (pow((delta_u_u), 2) - 2 * delta_u_u * (delta_u_b - delta_u_u) + pow((delta_u_b - delta_u_u), 2) + pow((delta_v_u), 2) - 2 * delta_v_u * (delta_v_b - delta_v_u) + pow((delta_v_b - delta_v_u), 2) - delta * delta);
				x1_pol = (-b_pol + sqrt(b_pol * b_pol - 4 * a_pol * c_pol)) / (2 * a_pol);
				x2_pol = (-b_pol - sqrt(b_pol * b_pol - 4 * a_pol * c_pol)) / (2 * a_pol);

				//Raiz quadrada de numero negativo: o passo pertence ao primeiro segmento de linha
				if (isnan(x1_pol) || x1_pol < 1.0 || a_pol == 0)
				{
					tau = sqrt(delta * delta / (delta_u_u * delta_u_u + delta_v_u * delta_v_u));
					delta_u = tau * delta_u_u;
					delta_v = tau * delta_v_u;
					step_n = delta;
				}
				//Caso contrario, o passo pertence ao segundo segmento de linha
				else {
					if (x1_pol >= 1.0 && x1_pol <= 2.0) {
						tau = x1_pol;
					}
					else {
						tau = x2_pol;
					}
					delta_u = delta_u_u + (tau - 1) * (delta_u_b - delta_u_u);
					delta_v = delta_v_u + (tau - 1) * (delta_v_b - delta_v_u);
					step_n = delta;
				}
			}
		}
		//Se a matriz hessiana e singular
		else if (abs(min_eig) < tol_eig || isnan(min_eig))
		{
			//Cauchy point
			cp_u = delta / g_abs * g[0];
			cp_v = delta / g_abs * g[1];

			gbg = g[0] * (H(0, 0) * g[0] + H(0, 1) * g[1]) + g[1] * (H(1, 0) * g[0] + H(1, 1) * g[1]);

			if (gbg <= 0) {
				tau = 1.0;
			}
			else {
				tau = g_abs * g_abs * g_abs / (delta * gbg);
				if (tau > 1.0) {
					tau = 1.0;
				}
			}

			delta_u = tau * cp_u;
			delta_v = tau * cp_v;
			step_n = sqrt(delta_u * delta_u + delta_v * delta_v);
		}

		//Se a matriz hessiana e indefinida
		else {
			alpha = 1.5 * abs(min_eig);

			Matrix B(2, 2);
			B(0, 0) = H(0, 0) + alpha;
			B(0, 1) = H(0, 1);
			B(1, 0) = H(1, 0);
			B(1, 1) = H(1, 1) + alpha;

			fulleigen1(B, P, D, tol_eig);
			max_eig = D(0, 0);
			min_eig = D(1, 1);
			index = 1;
			if (D(1, 1) > max_eig) {
				min_eig = D(0, 0);
				max_eig = D(1, 1);
				index = 0;
			}

			//Se B e positivo definida
			if (min_eig > tol_eig || abs(min_eig) < tol_eig) {
				delta_v_b = (g[1] - g[0] * B(1, 0) / B(0, 0)) / (B(1, 1) - B(0, 1) * B(1, 0) / B(0, 0));
				delta_u_b = (g[0] - B(0, 1) * delta_v_b) / B(0, 0);
				step_n = sqrt(delta_u_b * delta_u_b + delta_v_b * delta_v_b);

				//O passo esta dentro da regiao de confianca?
				//Sim
				if (step_n <= delta) {

					qo[0] = P(0, index);
					qo[1] = P(0, index);

					gamma = -(delta_u_b * qo[0] + delta_v_b * qo[1]) + sqrt(pow((delta_u_b * qo[0] + delta_v_b * qo[1]), 2) + delta * delta - step_n * step_n);

					delta_u = delta_u_b + gamma * qo[0];
					delta_v = delta_v_b + gamma * qo[1];
					step_n = sqrt(delta_u * delta_u + delta_v * delta_v);
				}
				//Nao
				else {
					//Subproblema do espaco bidimensional

					// -g
					A[0][0] = g[0];
					A[1][0] = g[1];

					// -step_DBG_NEWton
					A[0][1] = delta_u_b;
					A[1][1] = delta_v_b;

					M[0][0] = -A[0][0] / g_abs;
					M[1][0] = -A[1][0] / g_abs;

					w[0] = delta_u_b - (delta_u_b * (-g[0]) + delta_v_b * (-g[1])) / (g_abs * g_abs) * (-g[0]);
					w[1] = delta_v_b - (delta_u_b * (-g[0]) + delta_v_b * (-g[1])) / (g_abs * g_abs) * (-g[1]);

					w_n = sqrt(w[0] * w[0] + w[1] * w[1]);

					M[0][1] = w[0] / w_n;
					M[1][1] = w[1] / w_n;

					Uo[0] = M[0][0] * (-g[0]) + M[1][0] * (-g[1]);
					Uo[1] = M[0][1] * (-g[0]) + M[1][1] * (-g[1]);

					G[0][0] = (M[0][0] * H(0, 0) + M[1][0] * H(1, 0)) * M[0][0] + (M[0][0] * H(0, 1) + M[1][0] * H(1, 1)) * M[1][0];
					G[0][1] = (M[0][0] * H(0, 0) + M[1][0] * H(1, 0)) * M[0][1] + (M[0][0] * H(0, 1) + M[1][0] * H(1, 1)) * M[1][1];
					G[1][0] = (M[0][1] * H(0, 0) + M[1][1] * H(1, 0)) * M[0][0] + (M[0][1] * H(0, 1) + M[1][1] * H(1, 1)) * M[1][0];
					G[1][1] = (M[0][1] * H(0, 0) + M[1][1] * H(1, 0)) * M[0][1] + (M[0][1] * H(0, 1) + M[1][1] * H(1, 1)) * M[1][1];

					//Polinomio de 4o grau
					c1 = delta * (-G[0][1] * delta + Uo[0]);
					c2 = 2 * delta * (-(-G[0][0] + G[1][1]) * delta + Uo[1]);
					c3 = 6 * delta * delta * G[0][1];
					c4 = 2 * delta * ((-G[0][0] + G[1][1]) * delta + Uo[1]);
					c5 = -G[0][1] * delta * delta - Uo[0] * delta;

					//https://github.com/sasamil/Quartic

					complex<double>* solutions = solve_quartic(c2 / c1, c3 / c1, c4 / c1, c5 / c1);

					first = true;

					l = 0.0;
					o = 0.0;

					lf = 0.0;
					of = 0.0;

					for (int i = 0; i < 4; i++)
					{
						if (solutions[i].imag() != 0.0) {

						}
						else {
							l = 2 * delta * solutions[i].real() / (1 + solutions[i].real() * solutions[i].real());
							o = delta * (1 - solutions[i].real() * solutions[i].real()) / (1 + solutions[i].real() * solutions[i].real());
							if (i == 0 || first == true) {
								lf = l;
								of = o;
								//objk = DBG_NEW ObjectiveFunction(n1, p01, Q1, Q1_i, n2, p02, Q2, Q2_i, rho, phi);
								quadr1 = 0.5 * ((Skl2[0][0](0, 0) - p(0, 0))*(Skl2[0][0](0, 0) - p(0, 0)) + (Skl2[0][0](1, 0) - p(1, 0))*(Skl2[0][0](1, 0) - p(1, 0))) + (Skl2[0][0](2, 0) - p(2, 0))*(Skl2[0][0](2, 0) - p(2, 0)) + Uo[0] * l + Uo[1] * o + 0.5 * (l * (l * G[0][0] + o * G[1][0]) + o * (l * G[0][1] + o * G[1][1]));
								first = false;
							}
							else {
								//objk = DBG_NEW ObjectiveFunction(n1, p01, Q1, Q1_i, n2, p02, Q2, Q2_i, rho, phi);
								quadr2 = 0.5 * ((Skl2[0][0](0, 0) - p(0, 0))*(Skl2[0][0](0, 0) - p(0, 0)) + (Skl2[0][0](1, 0) - p(1, 0))*(Skl2[0][0](1, 0) - p(1, 0))) + (Skl2[0][0](2, 0) - p(2, 0))*(Skl2[0][0](2, 0) - p(2, 0)) + Uo[0] * l + Uo[1] * o + 0.5 * (l * (l * G[0][0] + o * G[1][0]) + o * (l * G[0][1] + o * G[1][1]));
								if (quadr2 < quadr1) {
									quadr1 = quadr2;
									lf = l;
									of = o;
								}
							}
						}
					}

					delta_u = M[0][0] * lf + M[0][1] * of;
					delta_v = M[1][0] * lf + M[1][1] * of;
					step_n = delta;

					delete[]solutions;
				}
			}
			else {
				cout << "ERROR" << endl;
			}
		}

		obj1 = 0.5 * ((Skl2[0][0](0, 0) - p(0, 0))*(Skl2[0][0](0, 0) - p(0, 0)) + (Skl2[0][0](1, 0) - p(1, 0))*(Skl2[0][0](1, 0) - p(1, 0)) + (Skl2[0][0](2, 0) - p(2, 0))*(Skl2[0][0](2, 0) - p(2, 0)));
		quad1 = obj1;
		quad2 = obj1 - g[0] * delta_u - g[1] * delta_v + 0.5 * (delta_u * (H(0, 0) * delta_u + H(0, 1) * delta_v) + delta_v * (H(1, 0) * delta_u + H(1, 1) * delta_v));
		uf_test = uf + delta_u;
		vf_test = vf + delta_v;
		NURBSDerivatives(uf_test, vf_test, Skl2, d);
		obj2 = 0.5 * ((Skl2[0][0](0, 0) - p(0, 0))*(Skl2[0][0](0, 0) - p(0, 0)) + (Skl2[0][0](1, 0) - p(1, 0))*(Skl2[0][0](1, 0) - p(1, 0)) + (Skl2[0][0](2, 0) - p(2, 0))*(Skl2[0][0](2, 0) - p(2, 0)));

		//Verificacao da funcao modelo quadratica
		rt = (obj1 - obj2) / (quad1 - quad2);

		if (abs(quad1 - quad2) < pow(10, -12)) {
			uf = uf + delta_u;
			vf = vf + delta_v;
			it++;
			break;
		}

		//if ((abs(quad1 - quad2) < pow(10, -12) || abs(obj1 - obj2) < pow(10, -12))) {
			//rt = 1;
		//}
		//else {
			if (rt < 0.25 || isnan(rt)) {
				delta = 0.25 * delta;
   
		
									  
					  
							
					   
	 
			}
			else {
				if (rt > 0.75 && step_n == delta) {
					delta = 2 * delta;
					if (delta > delta_max) {
						delta = delta_max;
					}
				}
				else {
				}
			}
   
		//}

		if (rt <= eta) {
			it++;
			if (it > 100) {
				break;
			}
		}
		else {
			uf = uf + delta_u;
			vf = vf + delta_v;

			//Derivadas da superficie no novo ponto
			NURBSDerivatives(uf, vf, Skl2, d);

			//Gradiente no novo ponto
			g[0] = -(Skl2[1][0](0, 0) * (Skl2[0][0](0, 0) - p(0, 0)) + Skl2[1][0](1, 0) * (Skl2[0][0](1, 0) - p(1, 0)) + Skl2[1][0](2, 0) * (Skl2[0][0](2, 0) - p(2, 0)));
			g[1] = -(Skl2[0][1](0, 0) * (Skl2[0][0](0, 0) - p(0, 0)) + Skl2[0][1](1, 0) * (Skl2[0][0](1, 0) - p(1, 0)) + Skl2[0][1](2, 0) * (Skl2[0][0](2, 0) - p(2, 0)));

			g_abs = sqrt(g[0] * g[0] + g[1] * g[1]);
			err = g_abs;

			//Hessiana no no novo ponto
			H(0, 0) = Skl2[2][0](0, 0) * Skl2[0][0](0, 0) + Skl2[1][0](0, 0) * Skl2[1][0](0, 0) - Skl2[2][0](0, 0) * p(0, 0)
				+ Skl2[2][0](1, 0) * Skl2[0][0](1, 0) + Skl2[1][0](1, 0) * Skl2[1][0](1, 0) - Skl2[2][0](1, 0) * p(1, 0)
				+ Skl2[2][0](2, 0) * Skl2[0][0](2, 0) + Skl2[1][0](2, 0) * Skl2[1][0](2, 0) - Skl2[2][0](2, 0) * p(2, 0);

			H(0, 1) = Skl2[1][1](0, 0) * Skl2[0][0](0, 0) + Skl2[1][0](0, 0) * Skl2[0][1](0, 0) - Skl2[1][1](0, 0) * p(0, 0)
				+ Skl2[1][1](1, 0) * Skl2[0][0](1, 0) + Skl2[1][0](1, 0) * Skl2[0][1](1, 0) - Skl2[1][1](1, 0) * p(1, 0)
				+ Skl2[1][1](2, 0) * Skl2[0][0](2, 0) + Skl2[1][0](2, 0) * Skl2[0][1](2, 0) - Skl2[1][1](2, 0) * p(2, 0);

			H(1, 0) = H(0, 1);

			H(1, 1) = Skl2[0][2](0, 0) * Skl2[0][0](0, 0) + Skl2[0][1](0, 0) * Skl2[0][1](0, 0) - Skl2[0][2](0, 0) * p(0, 0)
				+ Skl2[0][2](1, 0) * Skl2[0][0](1, 0) + Skl2[0][1](1, 0) * Skl2[0][1](1, 0) - Skl2[0][2](1, 0) *p(1, 0)
				+ Skl2[0][2](2, 0) * Skl2[0][0](2, 0) + Skl2[0][1](2, 0) * Skl2[0][1](2, 0) - Skl2[0][2](2, 0) * p(2, 0);


			it++;
			if (it > 100) {
				break;
			}
		}
	}

	//Se nao convergiu com numero de iteracões dentro do limite

	if (it > 100) {
		conv = false;
		uf = 0.0;
		vf = 0.0;
		p_f(0, 0) = 0.0;
		p_f(1, 0) = 0.0;
		p_f(2, 0) = 0.0;
	}
	else {
		if (uf < U_knot_vector[U_order]) {
			uf = U_knot_vector[U_order];
		}
		if (uf > U_knot_vector[r - U_order]) {
			uf = U_knot_vector[r - U_order];
		}
		if (vf < V_knot_vector[V_order]) {
			vf = V_knot_vector[V_order];
		}
		if (vf > V_knot_vector[s - V_order]) {
			vf = V_knot_vector[s - V_order];
		}
		NURBSDerivatives(uf, vf, Skl2, d);
		p_f(0, 0) = Skl2[0][0](0, 0);
		p_f(1, 0) = Skl2[0][0](1, 0);
		p_f(2, 0) = Skl2[0][0](2, 0);

		// Antigo, mas que talvez seja util

		//Se convergiu com numero de iteracões dentro do limite, mas parametros estao fora do range
		/*if (uf < U_knot_vector[U_order] || uf > U_knot_vector[r - U_order] || vf < V_knot_vector[V_order] || vf > V_knot_vector[s - V_order])
		{
			if ((uf < U_knot_vector[U_order] || uf > U_knot_vector[r - U_order]) && (vf > V_knot_vector[V_order] && vf < V_knot_vector[s - V_order])) {
				if (uv == 0 /*&& ((abs((uf - uvf)/uvf) < pow(10,-12)) || ((abs((uf - uvf))) < pow(10, -12) && uvf == 0))*//*) {
					uf = uvf;
					NURBSDerivatives(uf, vf, Skl2, d);
					p_f(0, 0) = Skl2[0][0](0, 0);
					p_f(1, 0) = Skl2[0][0](1, 0);
					p_f(2, 0) = Skl2[0][0](2, 0);
				}
				else {
					conv = false;
					uf = 0.0;
					vf = 0.0;
					p_f(0, 0) = 0.0;
					p_f(1, 0) = 0.0;
					p_f(2, 0) = 0.0;
				}
			}
			else if ((vf < V_knot_vector[V_order] || vf > V_knot_vector[s - V_order]) && (uf > U_knot_vector[U_order] && uf < U_knot_vector[r - U_order])) {
				if (uv == 1/* && ((abs((vf - uvf) / uvf) < pow(10, -12)) || ((abs((vf - uvf))) < pow(10, -12) && uvf == 0))*//*) {
					vf = uvf;
					NURBSDerivatives(uf, vf, Skl2, d);
					p_f(0, 0) = Skl2[0][0](0, 0);
					p_f(1, 0) = Skl2[0][0](1, 0);
					p_f(2, 0) = Skl2[0][0](2, 0);
				}
				else {
					NURBSDerivatives(uf, vf, Skl2, d);
					p_f(0, 0) = Skl2[0][0](0, 0);
					p_f(1, 0) = Skl2[0][0](1, 0);
					p_f(2, 0) = Skl2[0][0](2, 0);
				}
			}
			else {
				conv = false;
				uf = 0.0;
				vf = 0.0;
				p_f(0, 0) = 0.0;
				p_f(1, 0) = 0.0;
				p_f(2, 0) = 0.0;
			}
		}
		//Se convergiu
		else {
			NURBSDerivatives(uf, vf, Skl2, d);
			p_f(0, 0) = Skl2[0][0](0, 0);
			p_f(1, 0) = Skl2[0][0](1, 0);
			p_f(2, 0) = Skl2[0][0](2, 0);
		}*/
	}

	for (int i = 0; i < d + 1; i++)
	{
		delete[]Skl2[i];
	}
	delete[]Skl2;
}

void NURBSSurface::NURBSPointProjection(double p[3], double & uf, double & vf, double * p_f, bool & conv, double & tol, double & err, bool &conv_total, Matrix * P_0)
{
	int d_n = 2;
	Matrix** data;
	data = new Matrix*[d_n + 1];
	double der_zeta_zeta;
	double der_zeta_theta;
	double der_theta_theta;
	for (int i = 0; i < d_n + 1; i++)
	{
		data[i] = new Matrix[d_n + 1];
		for (int j = 0; j < d_n + 1; j++)
			data[i][j] = Matrix(3);
	}
	der_zeta_zeta = sqrt(data[2][0](0, 0) * data[2][0](0, 0) + data[2][0](1, 0) * data[2][0](1, 0) + data[2][0](2, 0) * data[2][0](2, 0));
	der_zeta_theta = sqrt(data[1][1](0, 0) * data[1][1](0, 0) + data[1][1](1, 0) * data[1][1](1, 0) + data[1][1](2, 0) * data[1][1](2, 0));
	der_theta_theta = sqrt(data[0][2](0, 0) * data[0][2](0, 0) + data[0][2](1, 0) * data[0][2](1, 0) + data[0][2](2, 0) * data[0][2](2, 0));

	bool plane = false;

	if (der_zeta_zeta < 1e-10 /*&& der_zeta_theta < tol_small_1*/ && der_theta_theta < 1e-10) {
		plane = true;
	}

	for (int i = 0; i < d_n + 1; i++)
	{
		delete[]data[i];
	}
	delete[]data;
	/*-----------------------------------------------------------------------------------------------------------------*/

	/*-----------------------------------------------------------------------------------------------------------------*/
	// Convergencia
	conv = true;
	conv_total = true; // O processo de otimizacao chegou no final ou os parametros sairam do range durante o processo?
	/*-----------------------------------------------------------------------------------------------------------------*/

	/*-----------------------------------------------------------------------------------------------------------------*/
	// Consideracões iniciais
	// Degeneracao
	int s_free = P_0->getColumns();
	// Numero de nos da NURBS
	int r = (U_dim - 1) + U_order + 1;
	int s = (V_dim - 1) + V_order + 1;

	// Derivadas da superficie no ponto inicial
	int d = 2;
	Matrix **Skl2; // Matriz para o calculo dos pontos e variaveis NURBS
	Skl2 = new Matrix*[d + 1];
	for (int i = 0; i < d + 1; i++)
	{
		Skl2[i] = new Matrix[d + 1];
		for (int j = 0; j < d + 1; j++)
			Skl2[i][j] = Matrix(3);
	}

	// Degeneracao
	// Degeneracao dos dois parametros (u e v) da superficie NURBS: solucao conhecida
	if (s_free == 0) {
		NURBSDerivatives(uf, vf, Skl2, d);
		p_f[0] = Skl2[0][0](0, 0);
		p_f[1] = Skl2[0][0](1, 0);
		p_f[2] = Skl2[0][0](2, 0);
	}
	// Sem degeneracao ou degeneracao de apenas um dos parametros da superficie NURBS (u ou v): resolver projecao
	else {
		// Variaveis iniciais
		int it = 1; // Iteracões

		Matrix xk(2, 1); // Vetor de parametros da superficie NURBS
		xk(0, 0) = uf;
		xk(1, 0) = vf;
		Matrix pk(s_free, 1); // Passo nas variaveis u e v (ou so u ou so v, se houver degeneracao)
		double step_n = norm(pk); // Norma do passo
		Matrix pb(s_free, 1); // Passo nas variaveis u e v pelo metodo de Newton (ou so u ou so v, se houver degeneracao)
		Matrix pu(s_free, 1); // Passo nas variaveis u e v pelo metodo dogleg (ou so u ou so v, se houver degeneracao)

		NURBSDerivatives(uf, vf, Skl2, d);

		// Gradiente no ponto inicial
		Matrix grad(2, 1); // "-" Gradiente da funcao suporte
		grad(0, 0) = -(Skl2[1][0](0, 0) * (Skl2[0][0](0, 0) - p[0]) + Skl2[1][0](1, 0) * (Skl2[0][0](1, 0) - p[1]) + Skl2[1][0](2, 0) * (Skl2[0][0](2, 0) - p[2]));
		grad(1, 0) = -(Skl2[0][1](0, 0) * (Skl2[0][0](0, 0) - p[0]) + Skl2[0][1](1, 0) * (Skl2[0][0](1, 0) - p[1]) + Skl2[0][1](2, 0) * (Skl2[0][0](2, 0) - p[2]));
		Matrix grad_deg(s_free, 1); // "-" Gradiente degenerado da funcao suporte
		grad_deg = transp(*P_0) * grad;
		double g_abs = norm(grad_deg); // Norma do gradiente degenerado
		Matrix grad_deg_aux(s_free, 1); // Copia do gradiente degenerado
		grad_deg_aux = grad_deg;
		err = g_abs;

		// Hessiana no ponto inicial
		Matrix H(2, 2);
		H(0, 0) = Skl2[2][0](0, 0) * Skl2[0][0](0, 0) + Skl2[1][0](0, 0) * Skl2[1][0](0, 0) - Skl2[2][0](0, 0) * p[0]
			+ Skl2[2][0](1, 0) * Skl2[0][0](1, 0) + Skl2[1][0](1, 0) * Skl2[1][0](1, 0) - Skl2[2][0](1, 0) * p[1]
			+ Skl2[2][0](2, 0) * Skl2[0][0](2, 0) + Skl2[1][0](2, 0) * Skl2[1][0](2, 0) - Skl2[2][0](2, 0) * p[2];
		H(0, 1) = Skl2[1][1](0, 0) * Skl2[0][0](0, 0) + Skl2[1][0](0, 0) * Skl2[0][1](0, 0) - Skl2[1][1](0, 0) * p[0]
			+ Skl2[1][1](1, 0) * Skl2[0][0](1, 0) + Skl2[1][0](1, 0) * Skl2[0][1](1, 0) - Skl2[1][1](1, 0) * p[1]
			+ Skl2[1][1](2, 0) * Skl2[0][0](2, 0) + Skl2[1][0](2, 0) * Skl2[0][1](2, 0) - Skl2[1][1](2, 0) * p[2];
		H(1, 0) = H(0, 1);
		H(1, 1) = Skl2[0][2](0, 0) * Skl2[0][0](0, 0) + Skl2[0][1](0, 0) * Skl2[0][1](0, 0) - Skl2[0][2](0, 0) * p[0]
			+ Skl2[0][2](1, 0) * Skl2[0][0](1, 0) + Skl2[0][1](1, 0) * Skl2[0][1](1, 0) - Skl2[0][2](1, 0) * p[1]
			+ Skl2[0][2](2, 0) * Skl2[0][0](2, 0) + Skl2[0][1](2, 0) * Skl2[0][1](2, 0) - Skl2[0][2](2, 0) * p[2];
		Matrix H_deg(s_free, s_free); // Matriz hessiana degenerada da funcao suporte
		H_deg = transp(*P_0) * H * (*P_0);
		Matrix H_deg_aux(s_free, s_free); // Copia da matriz hessiana degenerada da funcao suporte
		H_deg_aux = H_deg;

		Matrix P_deg(s_free, s_free);
		Matrix D_deg(s_free, s_free);
		double tol_eig = 1e-14;
		fulleigen1(H_deg, P_deg, D_deg, tol_eig); // Autovalores da matriz hessiana degenerada
		double max_eig = D_deg(0, 0);
		double min_eig = D_deg(0, 0);
		int index = 0;
		for (int i = 1; i < s_free; i++)
		{
			if (D_deg(i, i) > max_eig) {
				max_eig = D_deg(i, i);
			}
			if (D_deg(i, i) < min_eig) {
				min_eig = D_deg(i, i);
				index = 1;
			}
		}
		double error = pow(10, -12) * abs(max_eig);
		tol = error;

		// Variaveis - TR
		double delta = 0.0;
		double delta_max = sqrt((U_knot_vector[r - U_order] - U_knot_vector[U_order]) * (U_knot_vector[r - U_order] - U_knot_vector[U_order]) / 100 + (V_knot_vector[s - V_order] - V_knot_vector[V_order]) * (V_knot_vector[s - V_order] - V_knot_vector[V_order]) / 100);

		if (delta == 0) {
			delta = delta_max / 10;
		}

		double eta = 0.25;
		double rt = 0.0;

		double obj1 = 0.0;
		double obj2 = 0.0;

		double quad1 = 0.0;
		double quad2 = 0.0;

		Matrix cp(s_free, 1); // Cauchy
		Matrix xk_test(s_free, 1); // Variaveis u e v (ou so u ou so v, se houver degeneracao)

		double tau = 0.0;
		double gbg = 0.0;

		// Polinomio 2o grau
		double a_pol = 0.0;
		double b_pol = 0.0;
		double c_pol = 0.0;
		double x1_pol = 0.0;
		double x2_pol = 0.0;
		double somatorio = 0.0;

		// Minimizacao em subespaco bidimensional - variaveis
		Matrix I(s_free, s_free);
		for (int i = 0; i < s_free; i++)
		{
			for (int j = 0; j < s_free; j++)
			{
				if (i == j) {
					I(i, j) = 1.0;
				}
			}
		}
		Matrix B(s_free, s_free);
		Matrix B_aux(s_free, s_free);
		Matrix A(2, 2);
		Matrix M(2, 2);
		Matrix w(2, 1);
		double w_n = norm(w);
		Matrix Uo(2, 1);
		Matrix G(2, 2);
		Matrix qo(2, 1);

		double c1 = 0.0;
		double c2 = 0.0;
		double c3 = 0.0;
		double c4 = 0.0;
		double c5 = 0.0;

		double l = 0.0;
		double o = 0.0;

		double lf = 0.0;
		double of = 0.0;

		double quadr1 = 0.0;

		double quadr2 = 0.0;

		double alpha = 0.0;
		double gamma = 0.0;

		bool first = true;

		//Newton com TR
		int *flag = new int(0);
		while (g_abs > error || step_n > error)
		{
			// Se a matriz hessiana e positivo definida
			if (min_eig > tol_eig)
			{
				pb = fullsystem(H_deg, grad_deg, flag);
				step_n = norm(pb);
				H_deg = H_deg_aux;
				grad_deg = grad_deg_aux;

				// O passo esta dentro da regiao de confianca?
				// Sim
				if (step_n <= delta) {
					pk = pb;
					step_n = norm(pk);
				}
				// Nao
				else {
					// Dogleg method
					gbg = norm(transp(grad_deg) * H_deg * grad_deg);
					pu = pow(g_abs, 2) / gbg * grad_deg;

					a_pol = norm(pb - pu) * norm(pb - pu);
					if (a_pol < pow(10, -14)) {
						a_pol = 0;
					}
					somatorio = 0;
					for (int i = 0; i < s_free; i++)
					{
						somatorio = somatorio + (pb(i, 0) - pu(i, 0)) * pu(i, 0);
					}
					b_pol = -2 * norm(pb - pu) * norm(pb - pu) + 2 * somatorio;
					c_pol = -delta * delta + norm(pb - pu) * norm(pb - pu) + norm(pu) * norm(pu) - 2 * somatorio;
					x1_pol = (-b_pol + sqrt(b_pol * b_pol - 4 * a_pol * c_pol)) / (2 * a_pol);
					x2_pol = (-b_pol - sqrt(b_pol * b_pol - 4 * a_pol * c_pol)) / (2 * a_pol);

					// Raiz quadrada de numero negativo: o passo pertence ao primeiro segmento de linha
					if (isnan(x1_pol) || x1_pol < 1.0 || a_pol == 0)
					{
						tau = sqrt(delta * delta / (norm(pu) * norm(pu)));
						pk = tau * pu;
						step_n = delta;
					}
					// Caso contrario, o passo pertence ao segundo segmento de linha
					else {
						if (x1_pol >= 1.0 && x1_pol <= 2.0) {
							tau = x1_pol;
						}
						else {
							tau = x2_pol;
						}
						pk = pu + (tau - 1) * (pb - pu);
						step_n = delta;
					}
				}
			}
			// Se a matriz hessiana e singular
			else if (abs(min_eig) < tol_eig || isnan(min_eig))
			{
				// Cauchy point
				cp = delta / g_abs * grad_deg;

				gbg = norm(transp(grad_deg) * H_deg * grad_deg);

				if (gbg <= 0) {
					tau = 1.0;
				}
				else {
					tau = g_abs * g_abs * g_abs / (delta * gbg);
					if (tau > 1.0) {
						tau = 1.0;
					}
				}

				pk = tau * cp;
				step_n = norm(pk);
			}

			// Se a matriz hessiana e indefinida
			else {
				alpha = 1.5 * abs(min_eig);
				B = H_deg + alpha * I;
				B_aux = B;

				fulleigen1(B, P_deg, D_deg, tol_eig);
				max_eig = D_deg(0, 0);
				min_eig = D_deg(0, 0);
				index = 0;
				for (int i = 1; i < s_free; i++)
				{
					if (D_deg(i, i) > max_eig) {
						max_eig = D_deg(i, i);
					}
					if (D_deg(i, i) < min_eig) {
						min_eig = D_deg(i, i);
						index = 1;
					}
				}
				if (s_free != 2) {
					pk = fullsystem(B, grad_deg, flag);
					step_n = norm(pk);
					B = B_aux;
					grad_deg = grad_deg_aux;
				}
				else {
					// Se B e positivo definida
					if (min_eig > tol_eig || abs(min_eig) < tol_eig) {
						pb = fullsystem(B, grad_deg, flag);
						step_n = norm(pb);
						B = B_aux;
						grad_deg = grad_deg_aux;

						// O passo esta dentro da regiao de confianca?
						// Sim
						if (step_n <= delta) {
							qo(0, 0) = P_deg(0, index);
							qo(1, 0) = P_deg(1, index);

							gamma = -1.0 * norm(transp(pb) * qo) + sqrt(norm(transp(pb) * qo) * norm(transp(pb) * qo) + delta * delta - norm(pb) * norm(pb));

							pk = pb + gamma * qo;
							step_n = norm(pk);
						}
						// Nao
						else {
							// Subproblema do espaco bidimensional

							// -g
							A(0, 0) = grad_deg(0, 0);
							A(1, 0) = grad_deg(1, 0);

							// -step_newton
							A(0, 1) = pb(0, 0);
							A(1, 1) = pb(1, 0);

							M(0, 0) = -A(0, 0) / g_abs;
							M(1, 0) = -A(1, 0) / g_abs;

							w = pb + norm(transp(pb) * (-1.0 * grad_deg)) / (g_abs * g_abs) * grad_deg;

							w_n = norm(w);

							M(0, 1) = w(0, 0) / w_n;
							M(1, 1) = w(1, 0) / w_n;

							Uo = transp(M) * (-1.0 * grad_deg);

							G = transp(M) * H_deg * M;

							// Polinomio de 4o grau
							c1 = delta * (-G(0, 1) * delta + Uo(0, 0));
							c2 = 2 * delta * (-(-G(0, 0) + G(1, 1)) * delta + Uo(1, 0));
							c3 = 6 * delta * delta * G(0, 1);
							c4 = 2 * delta * ((-G(0, 0) + G(1, 1)) * delta + Uo(1, 0));
							c5 = -G(0, 1) * delta * delta - Uo(0, 0) * delta;

							//https://github.com/sasamil/Quartic

							complex<double>* solutions = solve_quartic(c2 / c1, c3 / c1, c4 / c1, c5 / c1);

							first = true;

							l = 0.0;
							o = 0.0;

							lf = 0.0;
							of = 0.0;

							for (int i = 0; i < 4; i++)
							{
								if (solutions[i].imag() != 0.0) {

								}
								else {
									l = 2 * delta * solutions[i].real() / (1 + solutions[i].real() * solutions[i].real());
									o = delta * (1 - solutions[i].real() * solutions[i].real()) / (1 + solutions[i].real() * solutions[i].real());
									if (i == 0 || first == true) {
										lf = l;
										of = o;
										quadr1 = 0.5 * ((Skl2[0][0](0, 0) - p[0])*(Skl2[0][0](0, 0) - p[0]) + (Skl2[0][0](1, 0) - p[1])*(Skl2[0][0](1, 0) - p[1]) + (Skl2[0][0](2, 0) - p[2])*(Skl2[0][0](2, 0) - p[2])) + Uo(0, 0) * l + Uo(1, 0) * o + 0.5 * (l * (l * G(0, 0) + o * G(1, 0)) + o * (l * G(0, 1) + o * G(1, 1)));
										first = false;
									}
									else {
										quadr2 = 0.5 * ((Skl2[0][0](0, 0) - p[0])*(Skl2[0][0](0, 0) - p[0]) + (Skl2[0][0](1, 0) - p[1])*(Skl2[0][0](1, 0) - p[1]) + (Skl2[0][0](2, 0) - p[2])*(Skl2[0][0](2, 0) - p[2])) + Uo(0, 0) * l + Uo(1, 0) * o + 0.5 * (l * (l * G(0, 0) + o * G(1, 0)) + o * (l * G(0, 1) + o * G(1, 1)));
										if (quadr2 < quadr1) {
											quadr1 = quadr2;
											lf = l;
											of = o;
										}
									}
								}
							}

							pk(0, 0) = M(0, 0) * lf + M(0, 1) * of;
							pk(1, 0) = M(1, 0) * lf + M(1, 1) * of;
							step_n = delta;

							delete[]solutions;
						}
					}
					else {
						cout << "ERROR" << endl;
					}
				}
			}

			obj1 = 0.5 * ((Skl2[0][0](0, 0) - p[0])*(Skl2[0][0](0, 0) - p[0]) + (Skl2[0][0](1, 0) - p[1])*(Skl2[0][0](1, 0) - p[1]) + (Skl2[0][0](2, 0) - p[2])*(Skl2[0][0](2, 0) - p[2]));
			quad1 = obj1;
			quad2 = obj1 - norm(transp(grad_deg) * pk) + 0.5 * norm(transp(pk) * H_deg * pk);
			xk_test = xk + (*P_0) * pk;
			NURBSDerivatives(xk_test(0, 0), xk_test(1, 0), Skl2, d);
			obj2 = 0.5 * ((Skl2[0][0](0, 0) - p[0])*(Skl2[0][0](0, 0) - p[0]) + (Skl2[0][0](1, 0) - p[1])*(Skl2[0][0](1, 0) - p[1]) + (Skl2[0][0](2, 0) - p[2])*(Skl2[0][0](2, 0) - p[2]));

			//Verificacao da funcao modelo quadratica
			rt = (obj1 - obj2) / (quad1 - quad2);

			if (abs(quad1 - quad2) < pow(10, -12)) {
				xk = xk + (*P_0) * pk;
				uf = xk(0, 0);
				vf = xk(1, 0);
				if (plane == false && (uf < U_knot_vector[U_order] || uf > U_knot_vector[U_dim] || vf < V_knot_vector[V_order] || vf > V_knot_vector[V_dim])) {
					conv_total = false;
				}
				it++;
				break;
			}

			//if ((abs(quad1 - quad2) < pow(10, -12) || abs(obj1 - obj2) < pow(10, -12))) {
				//rt = 1;
			//}
			//else {
				if (rt < 0.25 || isnan(rt)) {
					delta = 0.25 * delta;
	
		 
									   
					   
							 
						
	  
				}
				else {
					if (rt > 0.75 && step_n == delta) {
						delta = 2 * delta;
						if (delta > delta_max) {
							delta = delta_max;
						}
					}
					else {
					}
				}
	
			//}

			if (rt <= eta) {
				it++;
				if (it > 100) {
					break;
				}
			}
			else {
				xk = xk + (*P_0) * pk;
				uf = xk(0, 0);
				vf = xk(1, 0);

				if (plane == false && (uf < U_knot_vector[U_order] || uf > U_knot_vector[U_dim] || vf < V_knot_vector[V_order] || vf > V_knot_vector[V_dim])) {
					conv_total = false;
					it++;
					break;
				}
				else {
					// Derivadas da superficie no novo ponto
					NURBSDerivatives(uf, vf, Skl2, d);

					// Gradiente no novo ponto
					grad(0, 0) = -(Skl2[1][0](0, 0) * (Skl2[0][0](0, 0) - p[0]) + Skl2[1][0](1, 0) * (Skl2[0][0](1, 0) - p[1]) + Skl2[1][0](2, 0) * (Skl2[0][0](2, 0) - p[2]));
					grad(1, 0) = -(Skl2[0][1](0, 0) * (Skl2[0][0](0, 0) - p[0]) + Skl2[0][1](1, 0) * (Skl2[0][0](1, 0) - p[1]) + Skl2[0][1](2, 0) * (Skl2[0][0](2, 0) - p[2]));

					grad_deg = transp(*P_0) * grad;
					g_abs = norm(grad_deg);
					grad_deg_aux = grad_deg;
					err = g_abs;

					//Hessiana no no novo ponto
					H(0, 0) = Skl2[2][0](0, 0) * Skl2[0][0](0, 0) + Skl2[1][0](0, 0) * Skl2[1][0](0, 0) - Skl2[2][0](0, 0) * p[0]
						+ Skl2[2][0](1, 0) * Skl2[0][0](1, 0) + Skl2[1][0](1, 0) * Skl2[1][0](1, 0) - Skl2[2][0](1, 0) * p[1]
						+ Skl2[2][0](2, 0) * Skl2[0][0](2, 0) + Skl2[1][0](2, 0) * Skl2[1][0](2, 0) - Skl2[2][0](2, 0) * p[2];

					H(0, 1) = Skl2[1][1](0, 0) * Skl2[0][0](0, 0) + Skl2[1][0](0, 0) * Skl2[0][1](0, 0) - Skl2[1][1](0, 0) * p[0]
						+ Skl2[1][1](1, 0) * Skl2[0][0](1, 0) + Skl2[1][0](1, 0) * Skl2[0][1](1, 0) - Skl2[1][1](1, 0) * p[1]
						+ Skl2[1][1](2, 0) * Skl2[0][0](2, 0) + Skl2[1][0](2, 0) * Skl2[0][1](2, 0) - Skl2[1][1](2, 0) * p[2];

					H(1, 0) = H(0, 1);

					H(1, 1) = Skl2[0][2](0, 0) * Skl2[0][0](0, 0) + Skl2[0][1](0, 0) * Skl2[0][1](0, 0) - Skl2[0][2](0, 0) * p[0]
						+ Skl2[0][2](1, 0) * Skl2[0][0](1, 0) + Skl2[0][1](1, 0) * Skl2[0][1](1, 0) - Skl2[0][2](1, 0) * p[1]
						+ Skl2[0][2](2, 0) * Skl2[0][0](2, 0) + Skl2[0][1](2, 0) * Skl2[0][1](2, 0) - Skl2[0][2](2, 0) * p[2];
					H_deg = transp(*P_0) * H * (*P_0);
					H_deg_aux = H_deg;

					fulleigen1(H_deg, P_deg, D_deg, tol_eig);
					max_eig = D_deg(0, 0);
					min_eig = D_deg(0, 0);
					index = 0;
					for (int i = 1; i < s_free; i++)
					{
						if (D_deg(i, i) > max_eig) {
							max_eig = D_deg(i, i);
						}
						if (D_deg(i, i) < min_eig) {
							min_eig = D_deg(i, i);
							index = 1;
						}
					}

					it++;
					if (it > 100) {
						break;
					}
				}
			}
		}

		if (plane && (uf < U_knot_vector[U_order] || uf > U_knot_vector[U_dim] || vf < V_knot_vector[V_order] || vf > V_knot_vector[V_dim])) {
			conv_total = false;
		}

		//Se nao convergiu com numero de iteracões dentro do limite

		if (it > 100 || isnan(uf) || isnan(vf)) {
			if (isnan(uf) || isnan(vf)) {
				conv = false;
				uf = 0.0;
				vf = 0.0;
				p_f[0] = 0.0;
				p_f[1] = 0.0;
				p_f[2] = 0.0;
			}
			else {
				conv = false;
				NURBSDerivatives(uf, vf, Skl2, d);
				p_f[0] = Skl2[0][0](0, 0);
				p_f[1] = Skl2[0][0](1, 0);
				p_f[2] = Skl2[0][0](2, 0);
			}
		}
		else {
			NURBSDerivatives(uf, vf, Skl2, d);
			p_f[0] = Skl2[0][0](0, 0);
			p_f[1] = Skl2[0][0](1, 0);
			p_f[2] = Skl2[0][0](2, 0);
			if (conv_total == false)
				conv = false;
		}
		delete flag;
	}

	for (int i = 0; i < d + 1; i++)
	{
		delete[]Skl2[i];
	}
	delete[]Skl2;
}

void NURBSSurface::NURBSSupportFunction(double* v_dir, double &uf, double &vf, double* p_f, double &obj, bool &conv, bool &conv_total, double &delta, Matrix* P_0)
{
	/*-----------------------------------------------------------------------------------------------------------------*/
	// Convergencia
	conv = true;
	conv_total = true; // O processo de otimizacao chegou no final ou os parametros sairam do range durante o processo?
	/*-----------------------------------------------------------------------------------------------------------------*/

	/*-----------------------------------------------------------------------------------------------------------------*/
	// Consideracões iniciais
	int s_free = P_0->getColumns();
	int r = (U_dim - 1) + U_order + 1; // Nos em u
	int s = (V_dim - 1) + V_order + 1; // Nos em v
	int d = U_order + V_order;
	Matrix **Skl2; // Matriz para o calculo dos pontos e variaveis NURBS
	Skl2 = new Matrix*[d + 1];
	for (int i = 0; i < d + 1; i++)
	{
		Skl2[i] = new Matrix[d + 1];
		for (int j = 0; j < d + 1; j++)
			Skl2[i][j] = Matrix(3);
	}
	/*-----------------------------------------------------------------------------------------------------------------*/

	/*-----------------------------------------------------------------------------------------------------------------*/
	// Degeneracao
	// Degeneracao dos dois parametros (u e v) da superficie NURBS: solucao conhecida
	if (s_free == 0) {
		NURBSDerivatives(uf, vf, Skl2, d);
		p_f[0] = Skl2[0][0](0, 0);
		p_f[1] = Skl2[0][0](1, 0);
		p_f[2] = Skl2[0][0](2, 0);
		obj = (Skl2[0][0](0, 0) * v_dir[0] + Skl2[0][0](1, 0) * v_dir[1] + Skl2[0][0](2, 0) * v_dir[2]);
	}
	// Sem degeneracao ou degeneracao de apenas um dos parametros da superficie NURBS (u ou v): resolver processo de otimizacao (minimizacao) da funcao suporte
	else {
		// Variaveis iniciais
		int it = 1; // Iteracões
		Matrix xk(2, 1); // Vetor de parametros da superficie NURBS
		xk(0, 0) = uf;
		xk(1, 0) = vf;
		Matrix pk(s_free, 1); // Passo nas variaveis u e v (ou so u ou so v, se houver degeneracao)
		double step_n = norm(pk); // Norma do passo
		Matrix pb(s_free, 1); // Passo nas variaveis u e v pelo metodo de Newton (ou so u ou so v, se houver degeneracao)
		Matrix pu(s_free, 1); // Passo nas variaveis u e v pelo metodo dogleg (ou so u ou so v, se houver degeneracao)

		NURBSDerivatives(uf, vf, Skl2, d);

		Matrix grad(2, 1); // "-" Gradiente da funcao suporte
		grad(0, 0) = (Skl2[1][0](0, 0) * v_dir[0] + Skl2[1][0](1, 0) * v_dir[1] + Skl2[1][0](2, 0) * v_dir[2]);
		grad(1, 0) = (Skl2[0][1](0, 0) * v_dir[0] + Skl2[0][1](1, 0) * v_dir[1] + Skl2[0][1](2, 0) * v_dir[2]);
		Matrix grad_deg(s_free, 1); // "-" Gradiente degenerado da funcao suporte
		grad_deg = transp(*P_0) * grad;
		double g_abs = norm(grad_deg); // Norma do gradiente degenerado
		Matrix grad_deg_aux(s_free, 1); // Copia do gradiente degenerado
		grad_deg_aux = grad_deg;

		Matrix H(2, 2); // Matriz hessiana da funcao suporte
		H(0, 0) = -(Skl2[2][0](0, 0) * v_dir[0] + Skl2[2][0](1, 0) * v_dir[1] + Skl2[2][0](2, 0) * v_dir[2]);
		H(0, 1) = -(Skl2[1][1](0, 0) * v_dir[0] + Skl2[1][1](1, 0) * v_dir[1] + Skl2[1][1](2, 0) * v_dir[2]);
		H(1, 0) = H(0, 1);
		H(1, 1) = -(Skl2[0][2](0, 0) * v_dir[0] + Skl2[0][2](1, 0) * v_dir[1] + Skl2[0][2](2, 0) * v_dir[2]);
		Matrix H_deg(s_free, s_free); // Matriz hessiana degenerada da funcao suporte
		H_deg = transp(*P_0) * H * (*P_0);
		Matrix H_deg_aux(s_free, s_free); // Copia da matriz hessiana degenerada da funcao suporte
		H_deg_aux = H_deg;

		Matrix P_deg(s_free, s_free);
		Matrix D_deg(s_free, s_free);
		double tol_eig = 1e-14;
		fulleigen1(H_deg, P_deg, D_deg, tol_eig); // Autovalores da matriz hessiana degenerada
		double max_eig = D_deg(0, 0);
		double min_eig = D_deg(0, 0);
		int index = 0;
		for (int i = 1; i < s_free; i++)
		{
			if (D_deg(i, i) > max_eig) {
				max_eig = D_deg(i, i);
			}
			if (D_deg(i, i) < min_eig) {
				min_eig = D_deg(i, i);
				index = 1;
			}
		}
		double error = pow(10, -12) * abs(max_eig);

		// Variaveis - TR
		double delta_max = sqrt((U_knot_vector[r - U_order] - U_knot_vector[U_order]) * (U_knot_vector[r - U_order] - U_knot_vector[U_order]) / 100 + (V_knot_vector[s - V_order] - V_knot_vector[V_order]) * (V_knot_vector[s - V_order] - V_knot_vector[V_order]) / 100);

		if (delta == 0) {
			delta = delta_max / 10;
		}

		double eta = 0.25;
		double rt = 0.0;

		double obj1 = 0.0;
		double obj2 = 0.0;

		double quad1 = 0.0;
		double quad2 = 0.0;

		Matrix cp(s_free, 1); // Cauchy
		Matrix xk_test(s_free, 1); // Variaveis u e v (ou so u ou so v, se houver degeneracao)

		double tau = 0.0;
		double gbg = 0.0;

		// Polinomio 2o grau
		double a_pol = 0.0;
		double b_pol = 0.0;
		double c_pol = 0.0;
		double x1_pol = 0.0;
		double x2_pol = 0.0;
		double somatorio = 0.0;

		//Minimizacao em subespaco bidimensional - variaveis
		Matrix I(s_free, s_free);
		for (int i = 0; i < s_free; i++)
		{
			for (int j = 0; j < s_free; j++)
			{
				if (i == j) {
					I(i, j) = 1.0;
				}
			}
		}
		Matrix B(s_free, s_free);
		Matrix B_aux(s_free, s_free);
		Matrix A(2, 2);
		Matrix M(2, 2);
		Matrix w(2, 1);
		double w_n = norm(w);
		Matrix Uo(2, 1);
		Matrix G(2, 2);
		Matrix qo(2, 1);

		double c1 = 0.0;
		double c2 = 0.0;
		double c3 = 0.0;
		double c4 = 0.0;
		double c5 = 0.0;

		double l = 0.0;
		double o = 0.0;
		double lf = 0.0;
		double of = 0.0;

		double quadr1 = 0.0;
		double quadr2 = 0.0;

		double alpha = 0.0;
		double gamma = 0.0;

		bool first = true;

		// Newton com TR
		int *flag = new int(0);
		while (g_abs > error || step_n > error)
		{
			//Se a matriz hessiana e positivo definida
			if (min_eig > tol_eig)
			{
				pb = fullsystem(H_deg, grad_deg, flag);
				step_n = norm(pb);
				H_deg = H_deg_aux;
				grad_deg = grad_deg_aux;

				//O passo esta dentro da regiao de confianca?
				//Sim
				if (step_n <= delta) {
					pk = pb;
					step_n = norm(pk);
				}
				//Nao
				else {
					//Dogleg method
					gbg = norm(transp(grad_deg) * H_deg * grad_deg);
					pu = pow(g_abs, 2) / gbg * grad_deg;

					a_pol = norm(pb - pu) * norm(pb - pu);
					if (a_pol < pow(10, -14)) {
						a_pol = 0;
					}
					somatorio = 0;
					for (int i = 0; i < s_free; i++)
					{
						somatorio = somatorio + (pb(i, 0) - pu(i, 0)) * pu(i, 0);
					}
					b_pol = -2 * norm(pb - pu) * norm(pb - pu) + 2 * somatorio;
					c_pol = -delta * delta + norm(pb - pu) * norm(pb - pu) + norm(pu) * norm(pu) - 2 * somatorio;
					x1_pol = (-b_pol + sqrt(b_pol * b_pol - 4 * a_pol * c_pol)) / (2 * a_pol);
					x2_pol = (-b_pol - sqrt(b_pol * b_pol - 4 * a_pol * c_pol)) / (2 * a_pol);

					//Raiz quadrada de numero negativo: o passo pertence ao primeiro segmento de linha
					if (isnan(x1_pol) || x1_pol < 1.0 || a_pol == 0)
					{
						tau = sqrt(delta * delta / (norm(pu) * norm(pu)));
						pk = tau * pu;
						step_n = delta;
					}
					//Caso contrario, o passo pertence ao segundo segmento de linha
					else {
						if (x1_pol >= 1.0 && x1_pol <= 2.0) {
							tau = x1_pol;
						}
						else {
							tau = x2_pol;
						}
						pk = pu + (tau - 1) * (pb - pu);
						step_n = delta;
					}
				}
			}
			//Se a matriz hessiana e singular
			else if (abs(min_eig) < tol_eig || isnan(min_eig))
			{
				//Cauchy point
				cp = delta / g_abs * grad_deg;

				gbg = norm(transp(grad_deg) * H_deg * grad_deg);

				if (gbg <= 0) {
					tau = 1.0;
				}
				else {
					tau = g_abs * g_abs * g_abs / (delta * gbg);
					if (tau > 1.0) {
						tau = 1.0;
					}
				}

				pk = tau * cp;
				step_n = norm(pk);
			}
			//Se a matriz hessiana e indefinida
			else {
				alpha = 1.5 * abs(min_eig);
				B = H_deg + alpha * I;
				B_aux = B;

				fulleigen1(B, P_deg, D_deg, tol_eig);
				max_eig = D_deg(0, 0);
				min_eig = D_deg(0, 0);
				index = 0;
				for (int i = 1; i < s_free; i++)
				{
					if (D_deg(i, i) > max_eig) {
						max_eig = D_deg(i, i);
					}
					if (D_deg(i, i) < min_eig) {
						min_eig = D_deg(i, i);
						index = 1;
					}
				}

				if (s_free != 2) {
					pk = fullsystem(B, grad_deg, flag);
					step_n = norm(pk);
					B = B_aux;
					grad_deg = grad_deg_aux;
				}
				else {
					//Se B e positivo definida
					if (min_eig > tol_eig || abs(min_eig) < tol_eig) {
						pb = fullsystem(B, grad_deg, flag);
						step_n = norm(pb);
						B = B_aux;
						grad_deg = grad_deg_aux;

						//O passo esta dentro da regiao de confianca?
						//Sim
						if (step_n <= delta) {

							qo(0, 0) = P_deg(0, index);
							qo(1, 0) = P_deg(1, index);

							gamma = -1.0 * norm(transp(pb) * qo) + sqrt(norm(transp(pb) * qo) * norm(transp(pb) * qo) + delta * delta - norm(pb) * norm(pb));

							pk = pb + gamma * qo;
							step_n = norm(pk);
						}
						//Nao
						else {
							//Subproblema do espaco bidimensional

							// -g
							A(0, 0) = grad_deg(0, 0);
							A(1, 0) = grad_deg(1, 0);

							// -step_newton
							A(0, 1) = pb(0, 0);
							A(1, 1) = pb(1, 0);

							M(0, 0) = -A(0, 0) / g_abs;
							M(1, 0) = -A(1, 0) / g_abs;

							w = pb + norm(transp(pb) * (-1.0 * grad_deg)) / (g_abs * g_abs) * grad_deg;

							w_n = norm(w);

							M(0, 1) = w(0, 0) / w_n;
							M(1, 1) = w(1, 0) / w_n;

							Uo = transp(M) * (-1.0 * grad_deg);

							G = transp(M) * H_deg * M;

							//Polinomio de 4o grau
							c1 = delta * (-G(0, 1) * delta + Uo(0, 0));
							c2 = 2 * delta * (-(-G(0, 0) + G(1, 1)) * delta + Uo(1, 0));
							c3 = 6 * delta * delta * G(0, 1);
							c4 = 2 * delta * ((-G(0, 0) + G(1, 1)) * delta + Uo(1, 0));
							c5 = -G(0, 1) * delta * delta - Uo(0, 0) * delta;

							//https://github.com/sasamil/Quartic

							complex<double>* solutions = solve_quartic(c2 / c1, c3 / c1, c4 / c1, c5 / c1);

							first = true;

							l = 0.0;
							o = 0.0;

							lf = 0.0;
							of = 0.0;

							for (int i = 0; i < 4; i++)
							{
								if (solutions[i].imag() != 0.0) {

								}
								else {
									l = 2 * delta * solutions[i].real() / (1 + solutions[i].real() * solutions[i].real());
									o = delta * (1 - solutions[i].real() * solutions[i].real()) / (1 + solutions[i].real() * solutions[i].real());
									if (i == 0 || first == true) {
										lf = l;
										of = o;
										quadr1 = -(Skl2[0][0](0, 0) * v_dir[0] + Skl2[0][0](1, 0) * v_dir[1] + Skl2[0][0](2, 0) * v_dir[2]) + Uo(0, 0) * l + Uo(1, 0) * o + 0.5 * (l * (l * G(0, 0) + o * G(1, 0)) + o * (l * G(0, 1) + o * G(1, 1)));
										first = false;
									}
									else {
										quadr2 = -(Skl2[0][0](0, 0) * v_dir[0] + Skl2[0][0](1, 0) * v_dir[1] + Skl2[0][0](2, 0) * v_dir[2]) + Uo(0, 0) * l + Uo(1, 0) * o + 0.5 * (l * (l * G(0, 0) + o * G(1, 0)) + o * (l * G(0, 1) + o * G(1, 1)));
										if (quadr2 < quadr1) {
											quadr1 = quadr2;
											lf = l;
											of = o;
										}
									}
								}
							}

							pk(0, 0) = M(0, 0) * lf + M(0, 1) * of;
							pk(1, 0) = M(1, 0) * lf + M(1, 1) * of;
							step_n = delta;

							delete[]solutions;
						}
					}
					else {
						cout << "ERROR" << endl;
					}
				}
			}

			obj1 = -(Skl2[0][0](0, 0) * v_dir[0] + Skl2[0][0](1, 0) * v_dir[1] + Skl2[0][0](2, 0) * v_dir[2]);
			quad1 = obj1;
			quad2 = obj1 - norm(transp(grad_deg) * pk) + 0.5 * norm(transp(pk) * H_deg * pk);
			xk_test = xk + (*P_0) * pk;
			NURBSDerivatives(xk_test(0, 0), xk_test(1, 0), Skl2, d);
			obj2 = -(Skl2[0][0](0, 0)  * v_dir[0] + Skl2[0][0](1, 0) * v_dir[1] + Skl2[0][0](2, 0) * v_dir[2]);

			//Verificacao da funcao modelo quadratica
			rt = (obj1 - obj2) / (quad1 - quad2);

			if (abs(quad1 - quad2) < pow(10, -12)) {
				xk = xk + (*P_0) * pk;
				uf = xk(0, 0);
				vf = xk(1, 0);
				if (uf < U_knot_vector[U_order] || uf > U_knot_vector[U_dim] || vf < V_knot_vector[V_order] || vf > V_knot_vector[V_dim]) {
					conv_total = false;
				}
				it++;
				break;
			}

			//if ((abs(quad1 - quad2) < pow(10, -12) || abs(obj1 - obj2) < pow(10, -12))) {
				//rt = 1;
			//}
			//else {
				if (rt < 0.25 || isnan(rt)/*|| (objk->obj - objk_s->obj) < 0.0*/) {
					delta = 0.25 * delta;
	
		 
									   
					   
							 
						
	  
				}
				else {
					if (rt > 0.75 && step_n == delta) {
						delta = 2 * delta;
						if (delta > delta_max) {
							delta = delta_max;
						}
					}
					else {
					}
				}
	
			//}

			if (rt <= eta /*|| (objk->obj - objk_s->obj) < 0.0*/) {
				it++;
				if (it > 100) {
					break;
				}
			}
			else {
				xk = xk + (*P_0) * pk;
				uf = xk(0, 0);
				vf = xk(1, 0);

				if (uf < U_knot_vector[U_order] || uf > U_knot_vector[U_dim] || vf < V_knot_vector[V_order] || vf > V_knot_vector[V_dim]) {
					conv_total = false;
					it++;
					break;
				}
				else {
					NURBSDerivatives(uf, vf, Skl2, d);

					grad(0, 0) = (Skl2[1][0](0, 0) * v_dir[0] + Skl2[1][0](1, 0) * v_dir[1] + Skl2[1][0](2, 0) * v_dir[2]);
					grad(1, 0) = (Skl2[0][1](0, 0) * v_dir[0] + Skl2[0][1](1, 0) * v_dir[1] + Skl2[0][1](2, 0) * v_dir[2]);
					grad_deg = transp(*P_0) * grad;
					g_abs = norm(grad_deg);
					grad_deg_aux = grad_deg;

					H(0, 0) = -(Skl2[2][0](0, 0) * v_dir[0] + Skl2[2][0](1, 0) * v_dir[1] + Skl2[2][0](2, 0) * v_dir[2]);
					H(0, 1) = -(Skl2[1][1](0, 0) * v_dir[0] + Skl2[1][1](1, 0) * v_dir[1] + Skl2[1][1](2, 0) * v_dir[2]);
					H(1, 0) = H(0, 1);
					H(1, 1) = -(Skl2[0][2](0, 0) * v_dir[0] + Skl2[0][2](1, 0) * v_dir[1] + Skl2[0][2](2, 0) * v_dir[2]);
					H_deg = transp(*P_0) * H * (*P_0);
					H_deg_aux = H_deg;

					fulleigen1(H_deg, P_deg, D_deg, tol_eig);
					max_eig = D_deg(0, 0);
					min_eig = D_deg(0, 0);
					index = 0;
					for (int i = 1; i < s_free; i++)
					{
						if (D_deg(i, i) > max_eig) {
							max_eig = D_deg(i, i);
						}
						if (D_deg(i, i) < min_eig) {
							min_eig = D_deg(i, i);
							index = 1;
						}
					}

					it++;
					if (it > 100) {
						break;
					}
				}
			}
		}

		if (it > 100 || isnan(uf) || isnan(vf)) {
			if (isnan(uf) || isnan(vf)) {
				conv = false;
				uf = 0.0;
				vf = 0.0;
				p_f[0] = 0.0;
				p_f[1] = 0.0;
				p_f[2] = 0.0;
			}
			else {
				conv = false;
				NURBSDerivatives(uf, vf, Skl2, d);
				p_f[0] = Skl2[0][0](0, 0);
				p_f[1] = Skl2[0][0](1, 0);
				p_f[2] = Skl2[0][0](2, 0);
				obj = (Skl2[0][0](0, 0) * v_dir[0] + Skl2[0][0](1, 0) * v_dir[1] + Skl2[0][0](2, 0) * v_dir[2]);
			}
		}
		else {
			NURBSDerivatives(uf, vf, Skl2, d);
			p_f[0] = Skl2[0][0](0, 0);
			p_f[1] = Skl2[0][0](1, 0);
			p_f[2] = Skl2[0][0](2, 0);
			obj = (Skl2[0][0](0, 0) * v_dir[0] + Skl2[0][0](1, 0) * v_dir[1] + Skl2[0][0](2, 0) * v_dir[2]);
		}

		delete flag;
	}

	for (int i = 0; i < d + 1; i++)
	{
		delete[]Skl2[i];
	}
	delete[]Skl2;
}

//Normal exterior à superficie na posicao escolhida
void NURBSSurface::NormalExt(Matrix * Qp, double * zeta, double * theta, Matrix * n)
{
	Matrix** data;
	int d = 1;
	data = new Matrix*[d + 1];
	for (int i = 0; i < d + 1; i++)
	{
		data[i] = new Matrix[d + 1];
		for (int j = 0; j < d + 1; j++)
			data[i][j] = Matrix(3);
	}
	//material point (zeta,theta)
	NURBSDerivatives(*zeta, *theta, data, d);
	Matrix t1 = 1.0 / norm(data[1][0])*data[1][0];
	Matrix t2 = 1.0 / norm(data[0][1])*data[0][1];
	*n = 1.0 / norm(cross(t1, t2))*cross(t1, t2);
	//if (invert_normal == true)
		//*n = -1.0* (*n);
	*n = (*Qp) * (*n);
	for (int i = 0; i < d + 1; i++)
		delete[] data[i];
	delete[]data;
}

//Obtem ponto da superficie
void NURBSSurface::SurfacePoint(Matrix * Qp, Matrix * xp, double & zeta, double & theta, Matrix & point)
{
	Matrix p(3);
	NURBSPoint(zeta, theta, p);
	point = (*xp) + (*Qp) * p;
}

void NURBSSurface::EvaluateVolume()
{
	//TODO
	volume = 0.0;
}
void NURBSSurface::EvaluateCentroid()
{
	//TODO
	//Centroid in local coordinate system
	(*G)(0, 0) = 0.0;
	(*G)(1, 0) = 0.0;
	(*G)(2, 0) = 0.0;
}
void NURBSSurface::EvaluateInertiaTensor()
{
	//TODO
	//Inertia tensor with respect to the origin in local coordinate system
					//      | Jxx -Jxy -Jxz|     //
					//      |-Jxy  Jyy -Jyz|     //
					//      |-Jxz -Jyz  Jzz|     //
	(*J_O)(0, 0) = 0.0;
	(*J_O)(1, 0) = 0.0;
	(*J_O)(2, 0) = 0.0;
	(*J_O)(0, 1) = 0.0;
	(*J_O)(1, 1) = 0.0;
	(*J_O)(2, 1) = 0.0;
	(*J_O)(0, 2) = 0.0;
	(*J_O)(1, 2) = 0.0;
	(*J_O)(2, 2) = 0.0;

}
void NURBSSurface::EvaluateRadius()
{
	//Marina
	//Evaluates the particle radius (w/r to the origin)
	//Loop on control points
	radius = -1.0;
	double temp_len = 0.0;

	for (int j = 0; j < U_dim; j++)
	{
		for (int k = 0; k < V_dim; k++)
		{
			temp_len = sqrt(pow((control_points[j][k](0, 0) - (*G)(0, 0)), 2) + pow((control_points[j][k](1, 0) - (*G)(1, 0)), 2) + +pow((control_points[j][k](2, 0) - (*G)(2, 0)), 2));
			if (temp_len > radius) {
				radius = temp_len;
			}
			else {

			}
		}
	}
}

//Marina
bool NURBSSurface::ReadNurbsData()
{
	FILE *f1 = NULL;
	char name_file_1[500];
	strcpy(name_file_1, db.folder_name);
	strcat(name_file_1, "CAD/");
	strcat(name_file_1, file);
	//Correction on file termination - eliminating the ".txt" and replacing for ".snrb"
	char *temp;
	temp = strrchr(name_file_1, '.');   //Get the pointer to the last occurrence to the character '.'
	*temp = '\0';  //Replace token with null char
	strcat(name_file_1, ".snrb"); //Simplified/summarized file

	//Reading the file
	f1 = fopen(name_file_1, "r");
	if (f1 == NULL)
	{
		char name_file_2[500];
		strcpy(name_file_2, db.folder_name);
		strcat(name_file_2, "CAD/");
		strcat(name_file_2, file);
		//Correction on file termination - eliminating the ".txt" and replacing for ".rnrb"
		char *temp;
		temp = strrchr(name_file_2, '.');   //Get the pointer to the last occurrence to the character '.'
		*temp = '\0';  //Replace token with null char
		strcat(name_file_2, ".rnrb"); //Rhino file

		f1 = fopen(name_file_2, "r");

		if (f1 == NULL)
		{
			//Not available - return false
			return false;
		}
		else {
			//TO DO
			//Rhino file
			char ss[1000];

			fscanf(f1, "%s", ss);
			if (!strcmp(ss, "Volume"))
			{
				fscanf(f1, "%s", ss);
				fscanf(f1, "%s", ss);
				volume = atof(ss);
				fscanf(f1, "%s", ss);
				fscanf(f1, "%s", ss);
				fscanf(f1, "%s", ss);
				fscanf(f1, "%s", ss);
				fscanf(f1, "%s", ss);
				fscanf(f1, "%s", ss);
				if (!strcmp(ss, "Centroid"))
				{
					fscanf(f1, "%s", ss);
					fscanf(f1, "%[^,]s", ss);
					(*G)(0, 0) = atof(ss);
					fscanf(f1, "%c", ss);
					fscanf(f1, "%[^,]s", ss);
					(*G)(1, 0) = atof(ss);
					fscanf(f1, "%c", ss);
					fscanf(f1, "%s", ss);
					(*G)(2, 0) = atof(ss);
					fscanf(f1, "%s", ss);
					fscanf(f1, "%s", ss);
					fgets(ss, 200, f1);
					fgets(ss, 200, f1);
					fgets(ss, 200, f1);
					fgets(ss, 200, f1);
					fgets(ss, 200, f1);
					fgets(ss, 200, f1);
					fgets(ss, 200, f1);
					fgets(ss, 200, f1);
					fgets(ss, 200, f1);
					fgets(ss, 200, f1);
					fgets(ss, 200, f1);
					fscanf(f1, "%s", ss);
					if (!strcmp(ss, "xy:"))
					{
						fscanf(f1, "%s", ss);
						(*J_O)(0, 1) = -atof(ss);
						(*J_O)(1, 0) = (*J_O)(0, 1);
						fscanf(f1, "%s", ss);
						fscanf(f1, "%s", ss);
						fscanf(f1, "%s", ss);
						if (!strcmp(ss, "yz:"))
						{
							fscanf(f1, "%s", ss);
							(*J_O)(2, 1) = -atof(ss);
							(*J_O)(1, 2) = (*J_O)(2, 1);
							fscanf(f1, "%s", ss);
							fscanf(f1, "%s", ss);
							fscanf(f1, "%s", ss);
							if (!strcmp(ss, "zx:"))
							{
								fscanf(f1, "%s", ss);
								(*J_O)(0, 2) = -atof(ss);
								(*J_O)(2, 0) = (*J_O)(0, 2);
								fscanf(f1, "%s", ss);
								fscanf(f1, "%s", ss);
								fscanf(f1, "%s", ss);
								fgets(ss, 200, f1);
								fgets(ss, 200, f1);
								fgets(ss, 200, f1);
								fgets(ss, 200, f1);
								fgets(ss, 200, f1);
								fgets(ss, 200, f1);
								fgets(ss, 200, f1);
								fgets(ss, 200, f1);
								fgets(ss, 200, f1);
								fscanf(f1, "%s", ss);
								if (!strcmp(ss, "Ix:"))
								{
									fscanf(f1, "%s", ss);
									(*J_O)(0, 0) = atof(ss);
									fscanf(f1, "%s", ss);
									fscanf(f1, "%s", ss);
									fscanf(f1, "%s", ss);
									if (!strcmp(ss, "Iy:"))
									{
										fscanf(f1, "%s", ss);
										(*J_O)(1, 1) = atof(ss);
										fscanf(f1, "%s", ss);
										fscanf(f1, "%s", ss);
										fscanf(f1, "%s", ss);
										if (!strcmp(ss, "Iz:"))
										{
											fscanf(f1, "%s", ss);
											(*J_O)(2, 2) = atof(ss);
										}
										else {
											db.myprintf("\nError reading mesh file for CAD ID %d\n\n", number);
											fclose(f1);
											return false;
										}
									}
									else {
										db.myprintf("\nError reading mesh file for CAD ID %d\n\n", number);
										fclose(f1);
										return false;
									}
								}
								else {
									db.myprintf("\nError reading mesh file for CAD ID %d\n\n", number);
									fclose(f1);
									return false;
								}
							}
							else {
								db.myprintf("\nError reading mesh file for CAD ID %d\n\n", number);
								fclose(f1);
								return false;
							}
						}
						else {
							db.myprintf("\nError reading mesh file for CAD ID %d\n\n", number);
							fclose(f1);
							return false;
						}
					}
					else {
						db.myprintf("\nError reading mesh file for CAD ID %d\n\n", number);
						fclose(f1);
						return false;
					}
				}
				else {
					db.myprintf("\nError reading mesh file for CAD ID %d\n\n", number);
					fclose(f1);
					return false;
				}
			}
			else {
				db.myprintf("\nError reading mesh file for CAD ID %d\n\n", number);
				fclose(f1);
				return false;
			}
		}
		fclose(f1);
		return true;
	}
	else
	{
		//Simplified/summarized file
		char ss[1000];

		fscanf(f1, "%s", ss);
		if (!strcmp(ss, "VOLUME"))
		{
			fscanf(f1, "%s", ss);
			volume = (float)atof(ss);
			fscanf(f1, "%s", ss);
			if (!strcmp(ss, "CENTROID"))
			{
				fscanf(f1, "%s", ss);
				(*G)(0, 0) = (float)atof(ss);
				fscanf(f1, "%s", ss);
				(*G)(1, 0) = (float)atof(ss);
				fscanf(f1, "%s", ss);
				(*G)(2, 0) = (float)atof(ss);
				fscanf(f1, "%s", ss);
				fscanf(f1, "%s", ss);
				fscanf(f1, "%s", ss);
				if (!strcmp(ss, "INERTIA"))
				{
					fscanf(f1, "%s", ss);
					if (!strcmp(ss, "Ix"))
					{
						fscanf(f1, "%s", ss);
						(*J_O)(0, 0) = (float)atof(ss);
						fscanf(f1, "%s", ss);
						if (!strcmp(ss, "Iy"))
						{
							fscanf(f1, "%s", ss);
							(*J_O)(1, 1) = (float)atof(ss);
							fscanf(f1, "%s", ss);
							if (!strcmp(ss, "Iz"))
							{
								fscanf(f1, "%s", ss);
								(*J_O)(2, 2) = (float)atof(ss);
								fscanf(f1, "%s", ss);
								if (!strcmp(ss, "Ixy"))
								{
									fscanf(f1, "%s", ss);
									(*J_O)(0, 1) = -(float)atof(ss);
									(*J_O)(1, 0) = (*J_O)(0, 1);
									fscanf(f1, "%s", ss);
									if (!strcmp(ss, "Iyz"))
									{
										fscanf(f1, "%s", ss);
										(*J_O)(2, 1) = -(float)atof(ss);
										(*J_O)(1, 2) = (*J_O)(2, 1);
										fscanf(f1, "%s", ss);
										if (!strcmp(ss, "Izx"))
										{
											fscanf(f1, "%s", ss);
											(*J_O)(0, 2) = -(float)atof(ss);
											(*J_O)(2, 0) = (*J_O)(0, 2);
										}
										else {
											db.myprintf("\nError reading mesh file for CAD ID %d\n\n", number);
											fclose(f1);
											return false;
										}
									}
									else {
										db.myprintf("\nError reading mesh file for CAD ID %d\n\n", number);
										fclose(f1);
										return false;
									}
								}
								else {
									db.myprintf("\nError reading mesh file for CAD ID %d\n\n", number);
									fclose(f1);
									return false;
								}
							}
							else {
								db.myprintf("\nError reading mesh file for CAD ID %d\n\n", number);
								fclose(f1);
								return false;
							}
						}
						else {
							db.myprintf("\nError reading mesh file for CAD ID %d\n\n", number);
							fclose(f1);
							return false;
						}
					}
					else {
						db.myprintf("\nError reading mesh file for CAD ID %d\n\n", number);
						fclose(f1);
						return false;
					}
				}
				else {
					db.myprintf("\nError reading mesh file for CAD ID %d\n\n", number);
					fclose(f1);
					return false;
				}
			}
			else {
				db.myprintf("\nError reading mesh file for CAD ID %d\n\n", number);
				fclose(f1);
				return false;
			}
		}
		else
		{
			db.myprintf("\nError reading mesh file for CAD ID %d\n\n", number);
			fclose(f1);
			return false;
		}
		fclose(f1);
		return true;
	}

	//Marina
	delete f1;
	delete temp;
}

//Marina
void NURBSSurface::Subdivision()
{
	//Marina
	//n_sub = new NURBSSurface();
	//Parameter interval for each subdivision
	parameter_sub = new double**[subdivisions[0]];
	for (int i = 0; i < subdivisions[0]; i++)
	{
		parameter_sub[i] = new double*[subdivisions[1]];
	}
	for (int i = 0; i < subdivisions[0]; i++)
	{
		for (int j = 0; j < subdivisions[1]; j++)
		{
			parameter_sub[i][j] = new double[6];
		}

	}
	//Subdivisions - u
	double *u_b;
	double *u_multi;
	double *u_insert;
	if (subdivisions[0] == 1) {
		u_b = new double[subdivisions[0]];
		u_multi = new double[subdivisions[0]];
		u_insert = new double[subdivisions[0]];
	}
	else {
		u_b = new double[subdivisions[0] - 1];
		u_multi = new double[subdivisions[0] - 1];
		u_insert = new double[subdivisions[0] - 1];
		for (int i = 0; i < (subdivisions[0] - 1); i++)
		{
			u_multi[i] = 0;
			u_b[i] = U_knot_vector[U_order] + (i + 1) * (U_knot_vector[U_order] + U_knot_vector[U_order + (U_dim - 1) + 1 - U_order]) / subdivisions[0];
			for (int j = U_order; j <= (U_order + (U_dim - 1) + 1 - U_order); j++)
			{
				if (u_b[i] > (U_knot_vector[j] - pow(10, -16)) && u_b[i] < (U_knot_vector[j] + pow(10, -16))) {
					u_multi[i] = u_multi[i] + 1;
					u_b[i] = U_knot_vector[j];
				}
			}
			u_insert[i] = U_order - u_multi[i];
		}
	}
	//Subdivisions - v
	double *v_b;
	double *v_multi;
	double *v_insert;
	if (subdivisions[1] == 1) {
		v_b = new double[subdivisions[1]];
		v_multi = new double[subdivisions[1]];
		v_insert = new double[subdivisions[1]];
	}
	else {
		v_b = new double[subdivisions[1] - 1];
		v_multi = new double[subdivisions[1] - 1];
		v_insert = new double[subdivisions[1] - 1];
		for (int i = 0; i < (subdivisions[1] - 1); i++)
		{
			v_multi[i] = 0;
			v_b[i] = V_knot_vector[V_order] + (i + 1) * (V_knot_vector[V_order] + V_knot_vector[V_order + (V_dim - 1) + 1 - V_order]) / subdivisions[1];
			for (int j = V_order; j <= (V_order + (V_dim - 1) + 1 - V_order); j++)
			{
				if (v_b[i] > (V_knot_vector[j] - pow(10, -15)) && v_b[i] < (V_knot_vector[j] + pow(10, -15))) {
					v_multi[i] = v_multi[i] + 1;
					v_b[i] = V_knot_vector[j];
				}
			}
			v_insert[i] = V_order - v_multi[i];
		}
	}
	//parameter_sub
	if (subdivisions[0] == 1 && subdivisions[1] == 1) {
		parameter_sub[0][0][0] = U_knot_vector[U_order];
		parameter_sub[0][0][1] = U_knot_vector[U_dim];
		parameter_sub[0][0][2] = U_dim - 1;
		parameter_sub[0][0][3] = V_knot_vector[V_order];
		parameter_sub[0][0][4] = V_knot_vector[V_dim];
		parameter_sub[0][0][5] = V_dim - 1;
	}
	else if (subdivisions[0] == 1) {
		for (int j = 0; j < subdivisions[1]; j++)
		{
			parameter_sub[0][j][0] = U_knot_vector[U_order];
			parameter_sub[0][j][1] = U_knot_vector[U_dim];
			parameter_sub[0][j][2] = U_dim - 1;
			if (j == 0) {
				parameter_sub[0][j][3] = V_knot_vector[V_order];
				parameter_sub[0][j][4] = v_b[j];
			}
			else if (j == (subdivisions[1] - 1)) {
				parameter_sub[0][j][4] = V_knot_vector[V_dim];
				parameter_sub[0][j][3] = v_b[j - 1];
			}
			else {
				parameter_sub[0][j][3] = v_b[j - 1];
				parameter_sub[0][j][4] = v_b[j];
			}
		}
	}
	else if (subdivisions[1] == 1) {
		for (int i = 0; i < subdivisions[0]; i++)
		{
			parameter_sub[i][0][3] = V_knot_vector[V_order];
			parameter_sub[i][0][4] = V_knot_vector[V_dim];
			parameter_sub[i][0][5] = V_dim - 1;
			if (i == 0) {
				parameter_sub[i][0][0] = U_knot_vector[U_order];
				parameter_sub[i][0][1] = u_b[i];
			}
			else if (i == (subdivisions[0] - 1)) {
				parameter_sub[i][0][1] = U_knot_vector[U_dim];
				parameter_sub[i][0][0] = u_b[i - 1];
			}
			else {
				parameter_sub[i][0][0] = u_b[i - 1];
				parameter_sub[i][0][1] = u_b[i];
			}
		}
	}
	else {
		for (int i = 0; i < subdivisions[0]; i++)
		{
			for (int j = 0; j < subdivisions[1]; j++)
			{
				if (i == 0) {
					parameter_sub[i][j][0] = U_knot_vector[U_order];
					parameter_sub[i][j][1] = u_b[i];
				}
				else if (i == (subdivisions[0] - 1)) {
					parameter_sub[i][j][1] = U_knot_vector[U_dim];
					parameter_sub[i][j][0] = u_b[i - 1];
				}
				else {
					parameter_sub[i][j][0] = u_b[i - 1];
					parameter_sub[i][j][1] = u_b[i];
				}
				if (j == 0) {
					parameter_sub[i][j][3] = V_knot_vector[V_order];
					parameter_sub[i][j][4] = v_b[j];
				}
				else if (j == (subdivisions[1] - 1)) {
					parameter_sub[i][j][4] = V_knot_vector[V_dim];
					parameter_sub[i][j][3] = v_b[j - 1];
				}
				else {
					parameter_sub[i][j][3] = v_b[j - 1];
					parameter_sub[i][j][4] = v_b[j];
				}
			}
		}
	}

	double *U_knot_vector_sub_temp;
	double *V_knot_vector_sub_temp;
	Matrix **Pw_sub_temp;

	if (subdivisions[0] == 1) {
		n_sub_U_dim = U_dim;
		n_sub_U_order = U_order;
		n_sub_U_knot_vector = new double[n_sub_U_dim + n_sub_U_order + 1];
		U_knot_vector_sub_temp = new double[n_sub_U_dim + n_sub_U_order + 1];
		for (int i = 0; i < (U_dim + U_order + 1); i++)
		{
			n_sub_U_knot_vector[i] = U_knot_vector[i];
			U_knot_vector_sub_temp[i] = U_knot_vector[i];
		}
	}
	else {
		int ins_u = 0;
		for (int i = 0; i < (subdivisions[0] - 1); i++) {
			ins_u = ins_u + u_insert[i];
		}
		n_sub_U_dim = U_dim + ins_u;
		n_sub_U_order = U_order;
		n_sub_U_knot_vector = new double[n_sub_U_dim + n_sub_U_order + 1];
		U_knot_vector_sub_temp = new double[n_sub_U_dim + n_sub_U_order + 1];
		for (int i = 0; i < (U_dim + U_order + 1); i++)
		{
			n_sub_U_knot_vector[i] = U_knot_vector[i];
			U_knot_vector_sub_temp[i] = U_knot_vector[i];
		}

	}

	if (subdivisions[1] == 1) {
		n_sub_V_dim = V_dim;
		n_sub_V_order = V_order;
		n_sub_V_knot_vector = new double[n_sub_V_dim + n_sub_V_order + 1];
		V_knot_vector_sub_temp = new double[n_sub_V_dim + n_sub_V_order + 1];
		for (int i = 0; i < (V_dim + V_order + 1); i++)
		{
			n_sub_V_knot_vector[i] = V_knot_vector[i];
			V_knot_vector_sub_temp[i] = V_knot_vector[i];
		}
	}
	else {
		int ins_v = 0;
		for (int i = 0; i < (subdivisions[1] - 1); i++) {
			ins_v = ins_v + v_insert[i];
		}
		n_sub_V_dim = V_dim + ins_v;
		n_sub_V_order = V_order;
		n_sub_V_knot_vector = new double[n_sub_V_dim + n_sub_V_order + 1];
		V_knot_vector_sub_temp = new double[n_sub_V_dim + n_sub_V_order + 1];
		for (int i = 0; i < (V_dim + V_order + 1); i++)
		{
			n_sub_V_knot_vector[i] = V_knot_vector[i];
			V_knot_vector_sub_temp[i] = V_knot_vector[i];
		}

	}

	n_sub_weights = new double*[n_sub_U_dim];
	for (int i = 0; i < n_sub_U_dim; i++)
		n_sub_weights[i] = new double[n_sub_V_dim];

	n_sub_control_points = new Matrix*[n_sub_U_dim];
	for (int i = 0; i < n_sub_U_dim; i++)
	{
		n_sub_control_points[i] = new Matrix[n_sub_V_dim];
		for (int j = 0; j < n_sub_V_dim; j++)
			n_sub_control_points[i][j] = Matrix(3);
	}

	n_sub_Pw = new Matrix*[n_sub_U_dim];
	for (int i = 0; i < n_sub_U_dim; i++)
	{
		n_sub_Pw[i] = new Matrix[n_sub_V_dim];
		for (int j = 0; j < n_sub_V_dim; j++)
			n_sub_Pw[i][j] = Matrix(4);
	}

	Pw_sub_temp = new Matrix*[n_sub_U_dim];
	for (int i = 0; i < n_sub_U_dim; i++)
	{
		Pw_sub_temp[i] = new Matrix[n_sub_V_dim];
		for (int j = 0; j < n_sub_V_dim; j++)
			Pw_sub_temp[i][j] = Matrix(4);
	}

	for (int k = 0; k < U_dim; k++)
	{
		for (int l = 0; l < V_dim; l++)
		{
			n_sub_Pw[k][l](0, 0) = Pw[k][l](0, 0);
			n_sub_Pw[k][l](1, 0) = Pw[k][l](1, 0);
			n_sub_Pw[k][l](2, 0) = Pw[k][l](2, 0);
			n_sub_Pw[k][l](3, 0) = Pw[k][l](3, 0);
		}
	}

	n_sub_Bin = new double*[n_sub_U_order + n_sub_V_order + 4];
	for (int i = 0; i < n_sub_U_order + n_sub_V_order + 4; i++)
		n_sub_Bin[i] = new double[n_sub_U_order + n_sub_V_order + 4];
	for (int i = 0; i < n_sub_U_order + n_sub_V_order + 4; i++)
		for (int j = 0; j < n_sub_U_order + n_sub_V_order + 4; j++)
			n_sub_Bin[i][j] = 0.0;

	n_sub_Bin[0][0] = 1.0;
	for (int i = 1; i < n_sub_U_order + n_sub_V_order + 4; i++)
	{
		for (int j = 0; j <= i; j++)
		{
			if (j == 0)
				n_sub_Bin[i][j] = 1.0;
			else
				n_sub_Bin[i][j] = n_sub_Bin[i - 1][j - 1] + n_sub_Bin[i - 1][j];
		}
	}

	int span = 0;
	int aux = 0;
	if (subdivisions[0] == 1 && subdivisions[1] == 1) {
		for (int i = 0; i < U_dim; i++) {
			for (int j = 0; j < V_dim; j++) {
				n_sub_weights[i][j] = weights[i][j];
				n_sub_control_points[i][j](0, 0) = control_points[i][j](0, 0);
				n_sub_control_points[i][j](1, 0) = control_points[i][j](1, 0);
				n_sub_control_points[i][j](2, 0) = control_points[i][j](2, 0);
				n_sub_Pw[i][j](0, 0) = Pw[i][j](0, 0);
				n_sub_Pw[i][j](1, 0) = Pw[i][j](1, 0);
				n_sub_Pw[i][j](2, 0) = Pw[i][j](2, 0);
				n_sub_Pw[i][j](3, 0) = Pw[i][j](3, 0);
			}
		}
	}
	else {
		if (subdivisions[1] == 1) {
			n_sub_U_dim = U_dim;
			for (int i = 0; i < (subdivisions[0] - 1); i++)
			{
				aux = n_sub_U_dim - 1;
				span = FindSpan(aux, n_sub_U_order, u_b[i], n_sub_U_knot_vector);
				parameter_sub[i][0][2] = span - u_multi[i];
				if (i == (subdivisions[0] - 2)) {
					parameter_sub[subdivisions[0] - 1][0][2] = span - u_multi[i] + U_order;
				}
				SurfaceKnotIns(n_sub_U_dim, n_sub_U_order, n_sub_U_knot_vector, ((n_sub_U_dim - 1) + n_sub_U_order + 1), n_sub_V_dim, n_sub_V_order, n_sub_V_knot_vector, ((n_sub_V_dim - 1) + n_sub_V_order + 1), n_sub_Pw, 1, u_b[i], span, u_multi[i], u_insert[i], n_sub_U_dim, U_knot_vector_sub_temp, n_sub_V_dim, V_knot_vector_sub_temp, Pw_sub_temp);
				for (int j = 0; j < (n_sub_U_dim + n_sub_U_order + 1); j++)
				{
					n_sub_U_knot_vector[j] = U_knot_vector_sub_temp[j];
				}
				for (int k = 0; k < n_sub_U_dim; k++)
				{
					for (int l = 0; l < n_sub_V_dim; l++)
					{
						n_sub_Pw[k][l](0, 0) = Pw_sub_temp[k][l](0, 0);
						n_sub_Pw[k][l](1, 0) = Pw_sub_temp[k][l](1, 0);
						n_sub_Pw[k][l](2, 0) = Pw_sub_temp[k][l](2, 0);
						n_sub_Pw[k][l](3, 0) = Pw_sub_temp[k][l](3, 0);
					}
				}
			}
		}
		else if (subdivisions[0] == 1) {
			n_sub_V_dim = V_dim;
			for (int i = 0; i < (subdivisions[1] - 1); i++)
			{
				aux = n_sub_V_dim - 1;
				span = FindSpan(aux, n_sub_V_order, v_b[i], n_sub_V_knot_vector);
				parameter_sub[0][i][5] = span - v_multi[i];
				if (i == (subdivisions[1] - 2)) {
					parameter_sub[0][subdivisions[1] - 1][5] = span - v_multi[i] + V_order;
				}
				SurfaceKnotIns(n_sub_U_dim, n_sub_U_order, n_sub_U_knot_vector, ((n_sub_U_dim - 1) + n_sub_U_order + 1), n_sub_V_dim, n_sub_V_order, n_sub_V_knot_vector, ((n_sub_V_dim - 1) + n_sub_V_order + 1), n_sub_Pw, 2, v_b[i], span, v_multi[i], v_insert[i], n_sub_U_dim, U_knot_vector_sub_temp, n_sub_V_dim, V_knot_vector_sub_temp, Pw_sub_temp);
				for (int j = 0; j < (n_sub_V_dim + n_sub_V_order + 1); j++)
				{
					n_sub_V_knot_vector[j] = V_knot_vector_sub_temp[j];
				}
				for (int k = 0; k < n_sub_U_dim; k++)
				{
					for (int l = 0; l < n_sub_V_dim; l++)
					{
						n_sub_Pw[k][l](0, 0) = Pw_sub_temp[k][l](0, 0);
						n_sub_Pw[k][l](1, 0) = Pw_sub_temp[k][l](1, 0);
						n_sub_Pw[k][l](2, 0) = Pw_sub_temp[k][l](2, 0);
						n_sub_Pw[k][l](3, 0) = Pw_sub_temp[k][l](3, 0);
					}
				}
			}
		}
		else {
			n_sub_U_dim = U_dim;
			n_sub_V_dim = V_dim;
			for (int i = 0; i < (subdivisions[0] - 1); i++)
			{
				aux = n_sub_U_dim - 1;
				span = FindSpan(aux, n_sub_U_order, u_b[i], n_sub_U_knot_vector);
				for (int j = 0; j < (subdivisions[1]); j++)
				{
					parameter_sub[i][j][2] = span - u_multi[i];
					if (i == (subdivisions[0] - 2)) {
						parameter_sub[subdivisions[0] - 1][j][2] = span - u_multi[i] + U_order;
					}
				}
				SurfaceKnotIns(n_sub_U_dim, n_sub_U_order, n_sub_U_knot_vector, ((n_sub_U_dim - 1) + n_sub_U_order + 1), n_sub_V_dim, n_sub_V_order, n_sub_V_knot_vector, ((n_sub_V_dim - 1) + n_sub_V_order + 1), n_sub_Pw, 1, u_b[i], span, u_multi[i], u_insert[i], n_sub_U_dim, U_knot_vector_sub_temp, n_sub_V_dim, V_knot_vector_sub_temp, Pw_sub_temp);
				for (int j = 0; j < (n_sub_U_dim + n_sub_U_order + 1); j++)
				{
					n_sub_U_knot_vector[j] = U_knot_vector_sub_temp[j];
				}
				for (int k = 0; k < n_sub_U_dim; k++)
				{
					for (int l = 0; l < n_sub_V_dim; l++)
					{
						n_sub_Pw[k][l](0, 0) = Pw_sub_temp[k][l](0, 0);
						n_sub_Pw[k][l](1, 0) = Pw_sub_temp[k][l](1, 0);
						n_sub_Pw[k][l](2, 0) = Pw_sub_temp[k][l](2, 0);
						n_sub_Pw[k][l](3, 0) = Pw_sub_temp[k][l](3, 0);
					}
				}
			}
			for (int i = 0; i < (subdivisions[1] - 1); i++)
			{
				aux = n_sub_V_dim - 1;
				span = FindSpan(aux, n_sub_V_order, v_b[i], n_sub_V_knot_vector);
				for (int j = 0; j < (subdivisions[0]); j++)
				{
					parameter_sub[j][i][5] = span - v_multi[i];
					if (i == (subdivisions[1] - 2)) {
						parameter_sub[j][subdivisions[1] - 1][5] = span - v_multi[i] + V_order;
					}
				}
				SurfaceKnotIns(n_sub_U_dim, n_sub_U_order, n_sub_U_knot_vector, ((n_sub_U_dim - 1) + n_sub_U_order + 1), n_sub_V_dim, n_sub_V_order, n_sub_V_knot_vector, ((n_sub_V_dim - 1) + n_sub_V_order + 1), n_sub_Pw, 2, v_b[i], span, v_multi[i], v_insert[i], n_sub_U_dim, U_knot_vector_sub_temp, n_sub_V_dim, V_knot_vector_sub_temp, Pw_sub_temp);
				for (int j = 0; j < (n_sub_V_dim + n_sub_V_order + 1); j++)
				{
					n_sub_V_knot_vector[j] = V_knot_vector_sub_temp[j];
				}
				for (int k = 0; k < n_sub_U_dim; k++)
				{
					for (int l = 0; l < n_sub_V_dim; l++)
					{
						n_sub_Pw[k][l](0, 0) = Pw_sub_temp[k][l](0, 0);
						n_sub_Pw[k][l](1, 0) = Pw_sub_temp[k][l](1, 0);
						n_sub_Pw[k][l](2, 0) = Pw_sub_temp[k][l](2, 0);
						n_sub_Pw[k][l](3, 0) = Pw_sub_temp[k][l](3, 0);
					}
				}
			}
		}
		for (int k = 0; k < n_sub_U_dim; k++)
		{
			for (int l = 0; l < n_sub_V_dim; l++)
			{
				n_sub_weights[k][l] = n_sub_Pw[k][l](3, 0);
				n_sub_control_points[k][l](0, 0) = n_sub_Pw[k][l](0, 0) / n_sub_Pw[k][l](3, 0);
				n_sub_control_points[k][l](1, 0) = n_sub_Pw[k][l](1, 0) / n_sub_Pw[k][l](3, 0);
				n_sub_control_points[k][l](2, 0) = n_sub_Pw[k][l](2, 0) / n_sub_Pw[k][l](3, 0);
			}
		}
	}

	/*for (int k = 0; k < n_sub_U_dim; k++)
	{
		for (int l = 0; l < n_sub_V_dim; l++)
		{
			cout << n_sub_Pw[k][l](0, 0) << " " << n_sub_Pw[k][l](1, 0) << " " << n_sub_Pw[k][l](2, 0) << " " << n_sub_Pw[k][l](3, 0) << endl;
		}
	}

	cout << endl;*/

	if (u_b != NULL)
		delete[] u_b;
	if (u_multi != NULL)
		delete[] u_multi;
	if (u_insert != NULL)
		delete[] u_insert;
	if (v_b != NULL)
		delete[] v_b;
	if (v_multi != NULL)
		delete[] v_multi;
	if (v_insert != NULL)
		delete[] v_insert;
	if (U_knot_vector_sub_temp != NULL)
		delete[] U_knot_vector_sub_temp;
	if (V_knot_vector_sub_temp != NULL)
		delete[] V_knot_vector_sub_temp;
	if (Pw_sub_temp != NULL)
	{
		for (int i = 0; i < n_sub_U_dim; i++)
			delete[] Pw_sub_temp[i];
		delete[] Pw_sub_temp;
	}
}

//Marina
void NURBSSurface::ControlPointsSub()
{
	box_sub_coord = new double**[subdivisions[0]];
	box_sub_center = new MatrixFloat*[subdivisions[0]];
	box_sub_points = new MatrixFloat **[subdivisions[0]];
	halfedge_lengths = new MatrixFloat *[subdivisions[0]];
	for (int i = 0; i < subdivisions[0]; i++)
	{
		box_sub_coord[i] = new double*[subdivisions[1]];
		box_sub_center[i] = new MatrixFloat[subdivisions[1]];
		box_sub_points[i] = new MatrixFloat*[subdivisions[1]];
		halfedge_lengths[i] = new MatrixFloat[subdivisions[1]];
	}

	for (int i = 0; i < subdivisions[0]; i++)
	{
		for (int j = 0; j < subdivisions[1]; j++) {
			box_sub_coord[i][j] = new double[6];
			box_sub_center[i][j] = MatrixFloat(3);
			box_sub_points[i][j] = new MatrixFloat[8];
			halfedge_lengths[i][j] = MatrixFloat(3);
		}
	}

	for (int i = 0; i < subdivisions[0]; i++)
	{
		for (int j = 0; j < subdivisions[1]; j++) {
			for (int k = 0; k < 8; k++)
			{
				box_sub_points[i][j][k] = MatrixFloat(3);
			}
		}
	}

	float aux_x_min;
	float aux_x_max;;

	float aux_y_min;
	float aux_y_max;

	float aux_z_min;
	float aux_z_max;

	for (int i = 0; i < subdivisions[0]; i++)
	{
		for (int j = 0; j < subdivisions[1]; j++) {
			aux_x_min = n_sub_control_points[(int)parameter_sub[i][j][2]][(int)parameter_sub[i][j][5]](0, 0);
			aux_x_max = n_sub_control_points[(int)parameter_sub[i][j][2]][(int)parameter_sub[i][j][5]](0, 0);

			aux_y_min = n_sub_control_points[(int)parameter_sub[i][j][2]][(int)parameter_sub[i][j][5]](1, 0);
			aux_y_max = n_sub_control_points[(int)parameter_sub[i][j][2]][(int)parameter_sub[i][j][5]](1, 0);

			aux_z_min = n_sub_control_points[(int)parameter_sub[i][j][2]][(int)parameter_sub[i][j][5]](2, 0);
			aux_z_max = n_sub_control_points[(int)parameter_sub[i][j][2]][(int)parameter_sub[i][j][5]](2, 0);

			if (i == 0) {
				for (int k = 0; k <= (int)parameter_sub[i][j][2]; k++)
				{
					if (j == 0) {
						for (int l = 0; l <= (int)parameter_sub[i][j][5]; l++)
						{
							if (n_sub_control_points[k][l](0, 0) < aux_x_min) {
								aux_x_min = n_sub_control_points[k][l](0, 0);
							}
							if (n_sub_control_points[k][l](0, 0) > aux_x_max) {
								aux_x_max = n_sub_control_points[k][l](0, 0);
							}
							if (n_sub_control_points[k][l](1, 0) < aux_y_min) {
								aux_y_min = n_sub_control_points[k][l](1, 0);
							}
							if (n_sub_control_points[k][l](1, 0) > aux_y_max) {
								aux_y_max = n_sub_control_points[k][l](1, 0);
							}
							if (n_sub_control_points[k][l](2, 0) < aux_z_min) {
								aux_z_min = n_sub_control_points[k][l](2, 0);
							}
							if (n_sub_control_points[k][l](2, 0) > aux_z_max) {
								aux_z_max = n_sub_control_points[k][l](2, 0);
							}
						}
					}
					else {
						for (int l = (int)parameter_sub[i][j - 1][5]; l < (int)parameter_sub[i][j][5]; l++)
						{
							if (n_sub_control_points[k][l](0, 0) < aux_x_min) {
								aux_x_min = n_sub_control_points[k][l](0, 0);
							}
							if (n_sub_control_points[k][l](0, 0) > aux_x_max) {
								aux_x_max = n_sub_control_points[k][l](0, 0);
							}
							if (n_sub_control_points[k][l](1, 0) < aux_y_min) {
								aux_y_min = n_sub_control_points[k][l](1, 0);
							}
							if (n_sub_control_points[k][l](1, 0) > aux_y_max) {
								aux_y_max = n_sub_control_points[k][l](1, 0);
							}
							if (n_sub_control_points[k][l](2, 0) < aux_z_min) {
								aux_z_min = n_sub_control_points[k][l](2, 0);
							}
							if (n_sub_control_points[k][l](2, 0) > aux_z_max) {
								aux_z_max = n_sub_control_points[k][l](2, 0);
							}
						}
					}

				}
			}
			else {
				for (int k = (int)parameter_sub[i - 1][j][2]; k < (int)parameter_sub[i][j][2]; k++)
				{
					if (j == 0) {
						for (int l = 0; l < (int)parameter_sub[i][j][5]; l++)
						{
							if (n_sub_control_points[k][l](0, 0) < aux_x_min) {
								aux_x_min = n_sub_control_points[k][l](0, 0);
							}
							if (n_sub_control_points[k][l](0, 0) > aux_x_max) {
								aux_x_max = n_sub_control_points[k][l](0, 0);
							}
							if (n_sub_control_points[k][l](1, 0) < aux_y_min) {
								aux_y_min = n_sub_control_points[k][l](1, 0);
							}
							if (n_sub_control_points[k][l](1, 0) > aux_y_max) {
								aux_y_max = n_sub_control_points[k][l](1, 0);
							}
							if (n_sub_control_points[k][l](2, 0) < aux_z_min) {
								aux_z_min = n_sub_control_points[k][l](2, 0);
							}
							if (n_sub_control_points[k][l](2, 0) > aux_z_max) {
								aux_z_max = n_sub_control_points[k][l](2, 0);
							}
						}
					}
					else {
						for (int l = (int)parameter_sub[i][j - 1][5]; l < (int)parameter_sub[i][j][5]; l++)
						{
							if (n_sub_control_points[k][l](0, 0) < aux_x_min) {
								aux_x_min = n_sub_control_points[k][l](0, 0);
							}
							if (n_sub_control_points[k][l](0, 0) > aux_x_max) {
								aux_x_max = n_sub_control_points[k][l](0, 0);
							}
							if (n_sub_control_points[k][l](1, 0) < aux_y_min) {
								aux_y_min = n_sub_control_points[k][l](1, 0);
							}
							if (n_sub_control_points[k][l](1, 0) > aux_y_max) {
								aux_y_max = n_sub_control_points[k][l](1, 0);
							}
							if (n_sub_control_points[k][l](2, 0) < aux_z_min) {
								aux_z_min = n_sub_control_points[k][l](2, 0);
							}
							if (n_sub_control_points[k][l](2, 0) > aux_z_max) {
								aux_z_max = n_sub_control_points[k][l](2, 0);
							}
						}
					}
				}
			}

			double diagonal = sqrt((aux_x_max - aux_x_min) * (aux_x_max - aux_x_min) + (aux_y_max - aux_y_min) * (aux_y_max - aux_y_min) + (aux_z_max - aux_z_min) * (aux_z_max - aux_z_min));
			if (aux_x_min == aux_x_max) {
				aux_x_min = aux_x_min - 0.001 * diagonal;
				aux_x_max = aux_x_max + 0.001 * diagonal;
			}
			box_sub_coord[i][j][0] = aux_x_min;
			box_sub_coord[i][j][1] = aux_x_max;
			if (aux_y_min == aux_y_max) {
				aux_y_min = aux_y_min - 0.001 * diagonal;
				aux_y_max = aux_y_max + 0.001 * diagonal;
			}
			box_sub_coord[i][j][2] = aux_y_min;
			box_sub_coord[i][j][3] = aux_y_max;
			if (aux_z_min == aux_z_max) {
				aux_z_min = aux_z_min - 0.001 * diagonal;
				aux_z_max = aux_z_max + 0.001 * diagonal;
			}
			box_sub_coord[i][j][4] = aux_z_min;
			box_sub_coord[i][j][5] = aux_z_max;

			box_sub_center[i][j](0, 0) = (aux_x_min + aux_x_max) / 2;
			box_sub_center[i][j](1, 0) = (aux_y_min + aux_y_max) / 2;
			box_sub_center[i][j](2, 0) = (aux_z_min + aux_z_max) / 2;

			box_sub_points[i][j][0](0, 0) = box_sub_coord[i][j][0];
			box_sub_points[i][j][0](1, 0) = box_sub_coord[i][j][2];
			box_sub_points[i][j][0](2, 0) = box_sub_coord[i][j][4];

			box_sub_points[i][j][1](0, 0) = box_sub_coord[i][j][1];
			box_sub_points[i][j][1](1, 0) = box_sub_coord[i][j][2];
			box_sub_points[i][j][1](2, 0) = box_sub_coord[i][j][4];

			box_sub_points[i][j][2](0, 0) = box_sub_coord[i][j][1];
			box_sub_points[i][j][2](1, 0) = box_sub_coord[i][j][3];
			box_sub_points[i][j][2](2, 0) = box_sub_coord[i][j][4];

			box_sub_points[i][j][3](0, 0) = box_sub_coord[i][j][0];
			box_sub_points[i][j][3](1, 0) = box_sub_coord[i][j][3];
			box_sub_points[i][j][3](2, 0) = box_sub_coord[i][j][4];

			box_sub_points[i][j][4](0, 0) = box_sub_coord[i][j][0];
			box_sub_points[i][j][4](1, 0) = box_sub_coord[i][j][2];
			box_sub_points[i][j][4](2, 0) = box_sub_coord[i][j][5];

			box_sub_points[i][j][5](0, 0) = box_sub_coord[i][j][1];
			box_sub_points[i][j][5](1, 0) = box_sub_coord[i][j][2];
			box_sub_points[i][j][5](2, 0) = box_sub_coord[i][j][5];

			box_sub_points[i][j][6](0, 0) = box_sub_coord[i][j][1];
			box_sub_points[i][j][6](1, 0) = box_sub_coord[i][j][3];
			box_sub_points[i][j][6](2, 0) = box_sub_coord[i][j][5];

			box_sub_points[i][j][7](0, 0) = box_sub_coord[i][j][0];
			box_sub_points[i][j][7](1, 0) = box_sub_coord[i][j][3];
			box_sub_points[i][j][7](2, 0) = box_sub_coord[i][j][5];

			halfedge_lengths[i][j](0, 0) = norm(box_sub_points[i][j][1] - box_sub_points[i][j][0]) / 2;
			halfedge_lengths[i][j](1, 0) = norm(box_sub_points[i][j][2] - box_sub_points[i][j][1]) / 2;
			halfedge_lengths[i][j](2, 0) = norm(box_sub_points[i][j][4] - box_sub_points[i][j][0]) / 2;

		}
	}
}
