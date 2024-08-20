#include "SecSuperEllipse.h"


SecSuperEllipse::SecSuperEllipse()
{
	sec_details = new SolidSection();

	AC = Matrix(2);
	aero_length = 0;
}


SecSuperEllipse::~SecSuperEllipse()
{
	delete sec_details;
}

bool SecSuperEllipse::Read(FILE *f)
{
	char s[1000];
	fscanf(f, "%s", s);
	number = atoi(s);
	fscanf(f, "%s", s);
	if (!strcmp(s, "A"))
	{
		fscanf(f, "%s", s);
		a = atof(s);
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "B"))
	{
		fscanf(f, "%s", s);
		b = atof(s);
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "N"))
	{
		fscanf(f, "%s", s);
		n = atof(s);
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "AMeshFDM"))
	{
		fscanf(f, "%s", s);
		div_mesh_a = atoi(s);
	}
	else
		return false;

	//Salva a posição (stream)
	fpos_t pos;
	fgetpos(f, &pos);
	fscanf(f, "%s", s);
	if (!strcmp(s, "AD"))
	{
		fscanf(f, "%s", s);
		aerodynamicdataID = atoi(s);

		fscanf(f, "%s", s);
		if (!strcmp(s, "AC"))
		{
			fscanf(f, "%s", s);
			AC(0, 0) = atof(s);
			fscanf(f, "%s", s);
			AC(1, 0) = atof(s);
		}
		else
			return false;

		fscanf(f, "%s", s);
		if (!strcmp(s, "AeroLength"))
		{
			fscanf(f, "%s", s);
			aero_length = atof(s);
		}
		else
			return false;
	}
	else
		fsetpos(f, &pos);	//volta à posição anterior


	return true;
}

void SecSuperEllipse::Write(FILE *f)
{
	fprintf(f, "SuperEllipse\t%d\tA\t%.6e\tB\t%.6e\tN\t%.6e\tAMeshFDM\t%d\tAD\t%d\n", number, a, b, n,div_mesh_a, aerodynamicdataID);
}

void SecSuperEllipse::PreCalc()
{
	SolidSection* ptr_sd = static_cast<SolidSection*>(sec_details);
	double eps = 2.0 / n;
	A = 2.0*a*b*eps*Beta(eps / 2.0, (eps+2.0)/2.0);
	I11 = 0.5*a*b*b*b*eps*Beta(3.0*eps / 2.0, eps / 2.0);
	I22 = 0.5*b*a*a*a*eps*Beta(3.0*eps / 2.0, eps / 2.0);
	I12 = 0.0;
	I33 = I11 + I22;

	int n_circ = 24;
	ptr_sd->Alloc(n_circ);
	double theta = 0;
	double phi;
	for (int index = 0; index < n_circ; index++)
	{
		theta = (index * 2 * PI) / n_circ;
		phi = 1.0 / (pow((pow(abs(sin(theta)), n) + pow(abs(cos(theta)), n)),1.0/n));
		ptr_sd->points[index][0] = phi*a*cos(theta);
		ptr_sd->points[index][1] = phi*b*sin(theta);
	}

	It = MDF_SaintVenantSE();
	printf("\nSaint-Venant torsion stiffness for section %d: %.6e\n", this->number, It);
}
double SecSuperEllipse::Beta(double p1, double p2)
{
	double product = 1.0;
	double previous = 2.0;
	double error = 1e-12;
	int k = 1;
	while (abs(product-previous)>error)
	{
		previous = product;
		product = product*(1.0 + (p1 + p2) / k) / ((1 + p1 / k)*(1 + p2 / k));
		k++;
	}
		
	return ((p1 + p2) / (p1*p2))*product;
}

double SecSuperEllipse::MDF_SaintVenantSE()
{
	//Malha do MDF
	int nx = div_mesh_a;
	int ny = (int)(div_mesh_a*b/a);

	//Matriz com informações condensadas de todos os pontos do grid
	double *** info_p;
	info_p = new double **[nx + 1];
	for (int i = 0; i < nx+1; i++)
		info_p[i] = new double *[ny + 1];
	for (int i = 0; i < nx+1; i++)
	for (int j = 0; j < ny+1; j++)
		info_p[i][j] = new double[6];

	//Índices utilizados
	//0 - não utilizado
	//1 - x
	//2 - y
	//3 - phi(x, y)
	//4 - status(0) - inativo - (1) ativo
	//5 - número da linha(equação)

	//%%%%%%%%%%%%%%%%%%%%%%%%%%GERAÇÃO DE MALHA%%%%%%%%%%%%%%%%%%%%%%%
	//percorre por linhas
	for (int j = 0; j <= ny;j++)
	{
		double y_cur = j*b / ny;
		//percorre por colunas
		for (int i = 0; i <= nx; i++)
		{
			info_p[i][j][1] = (a*i / nx);		//coordenada x
			info_p[i][j][2] = y_cur;			//coordenada y
			double xmax = a*pow(1.0 - pow(y_cur / b,n),1/n);
			if ((a*i / nx) <= xmax)
			{
				//ponto ativo
				info_p[i][j][4] = 1;
			}
			else
			{
				//ponto inativo
				info_p[i][j][4] = 0;
			}
		}
	}
	int n_GL = 0;
	for (int j = 0; j <= ny; j++)
	{
		for (int i = 0; i <= nx; i++)
		{
			if (info_p[i][j][4] == 1)
			{
				info_p[i][j][5] = n_GL;
				n_GL = n_GL + 1;
			}
		}
	}
	//Banda da matriz A
	int larger = nx;
	if (ny > nx)
		larger = ny;
	SparseMatrix A(n_GL,n_GL,larger*3);
	
	Matrix B(n_GL);
	//%%%%%%%%%%%%%%%%%%%%%%%%%%MONTAGEM DE MATRIZES%%%%%%%%%%%%%%%%%%%%%%%
	double hy = b / ny;
	double hx = a / nx;
	double coefx = 1.0 / (hx*hx);
	double coefy = 1.0 / (hy*hy);
	//percorre por linhas
	for (int j = 0; j <= ny; j++)
	{
		//percorre por colunas
		for (int i = 0; i <= nx; i++)
		{
			//Se o ponto se encontra no domínio
			if (info_p[i][j][4] == 1)
			{
				A.setValue((long)(info_p[i][j][5]), (long)(info_p[i][j][5]), -2.0 * coefx - 2.0 * coefy);
				B((long)(info_p[i][j][5]), 0) = -2;
				//ponto à esquerda - simetria - eixo y
				if (i == 0)
				{
					if (i + 1 <= nx)
					{
						if (info_p[i + 1][j][4] == 1)
							A.setValue((long)(info_p[i][j][5]), (long)(info_p[i + 1][j][5]), 2.0 * coefx);
					}
				}
				//ponto no interior do domínio ou no contorno
				else
				{
					if (i + 1 <= nx)
					{
						if (info_p[i + 1][j][4] == 1)
						{
							A.setValue((long)(info_p[i][j][5]), (long)(info_p[i + 1][j][5]), 1 * coefx);
						}
					}
					if (info_p[i - 1][j][4] == 1)
					{
						A.setValue((long)(info_p[i][j][5]), (long)(info_p[i - 1][j][5]), 1 * coefx);
					}
				}

				//ponto em y = 0 - simetria - eixo x
				if (j == 0)
				{
					if (j + 1 <= ny)
					{
						if (info_p[i][j + 1][4] == 1)
						{
							A.setValue((long)(info_p[i][j][5]), (long)(info_p[i][j + 1][5]), 2 * coefy);
						}
					}
				}
				//ponto no interior do domínio ou no contorno
				else
				{
					if (j + 1 <= ny)
					{
						if (info_p[i][j + 1][4] == 1)
						{
							A.setValue((long)(info_p[i][j][5]), (long)(info_p[i][j + 1][5]), 1 * coefy);
						}
					}
					if (info_p[i][j - 1][4] == 1)
					{
						A.setValue((long)(info_p[i][j][5]), (long)(info_p[i][j - 1][5]), 1 * coefy);
					}
				}
			}
		}
	}
	
	int flag;
	Matrix phi = sparsesystem(A, B, &flag,1,0);

	//composição do phi em info_p
	//percorre por linhas
	for (int j = 0; j <= ny; j++)
	{
		//percorre por colunas
		for (int i = 0; i <= nx; i++)
		{
			if (info_p[i][j][4] == 1)
				info_p[i][j][3] = phi((long)(info_p[i][j][5]), 0);
		}

	}

	//Integração
	double temp_int = 0;
	for (int j = 0; j <= ny; j++)
	{
		//percorre por colunas
		for (int i = 0; i <= nx; i++)
		{
			// Se o ponto se encontra no domínio
			if (info_p[i][j][4] == 1)
			{
				if ((i == 0 && j == 0) || (i == 0 && j == ny) || (i == nx && j == 0) || (i == nx && j == ny)) // ponto nos cantos
				{
					temp_int = temp_int + info_p[i][j][3] * (hx*hy) / 4;
				}
				else
				{
					if (i == 0 || j == 0 || j == ny || i == nx) // ponto nas bordas
					{
						temp_int = temp_int + info_p[i][j][3] * (hx*hy) / 2;
					}
					else
					{
						temp_int = temp_int + info_p[i][j][3] * (hx*hy);
					}
				}
			}
		}
	}

	for (int i = 0; i < nx+1; i++)
	for (int j = 0; j < ny+1; j++)
		delete[] info_p[i][j];
	for (int i = 0; i < nx + 1; i++)
		delete[] info_p[i];
	delete[] info_p;

	return (2 * temp_int * 4);
}
