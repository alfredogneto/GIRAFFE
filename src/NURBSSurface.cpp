#include "NURBSSurface.h"
#include <math.h>

#include "Encoding.h"

#include"Database.h"
//Variáveis globais
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

//Escrita do arquivo de saída
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

	EvaluateRadius();
	EvaluateVolume();
	EvaluateCentroid();
	EvaluateInertiaTensor();
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
	for (int i = 0; i < du; i++)
		delete[] Nu[i];
	delete[] Nu;
	for (int i = 0; i < dv; i++)
		delete[] Nv[i];
	delete[] Nv;
	delete[] temp;
	for (int i = 0; i < d; i++)
		delete[] Skl1[i];
	delete[] Skl1;
}

void NURBSSurface::WriteVTK_XMLRender(FILE *f, Matrix& pos, Matrix& rot, int number)
{
	//vetores para escrita no formato binário - usando a função 'enconde'
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
	//TODO
	radius = 0.0;
}
