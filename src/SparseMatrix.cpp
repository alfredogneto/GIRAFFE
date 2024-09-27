#include "SparseMatrix.h"
#include <iostream>
#include <mkl.h>

#include <arpack.hpp>
#include <debug_c.hpp>  // debug arpack.
#include <stat_c.hpp>   // arpack statistics.

#include "Matrix.h"
#include "SolverOptions.h"
#include "Database.h"
//Variaveis globais
extern
Database db;

SparseMatrix::SparseMatrix()
{
	mounted = false;
	rows = 1;
	cols = 1;
	non_null_estimative = 1;
	m_matrix.resize(1, 1);				//tamanho da matriz
	tripletList.reserve(1);				//pre-aloca o triplet list
}
SparseMatrix::SparseMatrix(int e_rows,int e_cols, int e_non_null_estimative)
{
	rows = e_rows;
	cols = e_cols;
	non_null_estimative = e_non_null_estimative;
	m_matrix.resize(e_rows, e_cols);					//tamanho da matriz
	m_matrix.reserve(e_non_null_estimative);
	tripletList.reserve(2*e_non_null_estimative);		//pre-aloca o triplet list
	mounted = false;
}
SparseMatrix::SparseMatrix(SparseMatrix &copied)
{
	this->rows = copied.rows;
	this->cols = copied.cols;
	this->m_matrix = copied.m_matrix;
	this->non_null_estimative = copied.non_null_estimative;
	this->mounted = copied.mounted;
	this->tripletList = copied.tripletList;
}

SparseMatrix::~SparseMatrix()
{
	//Não ha variaveis a serem desalocadas manualmente (todas são automaticas)
}

//Zera coeficientes da matriz esparsa (alocados na tripletList)
void SparseMatrix::Clear()
{
	tripletList.erase(tripletList.begin(), tripletList.end());
	mounted = false;
}

//Seta na posição (i,j) da matriz quadrada o valor v
void SparseMatrix::setValue(int i, int j, double v)
{
	if (i + 1 > m_matrix.rows() || j + 1 > m_matrix.cols())
		printf("Error assigning sparse matrix value\n");
	else
		tripletList.push_back(T(i, j, v));
	if (mounted)//se houve inserção após montagem, seta montagem como falso
		mounted = false;
}
void SparseMatrix::Mount()
{
	m_matrix.setFromTriplets(tripletList.begin(), tripletList.end());		//preparação da matriz - montagem a partir do triplet list
	mounted = true;
}

//Resolve o sistema linear da forma Ax=b
Matrix sparsesystem(SparseMatrix &A, Matrix &b, int *info_fail,int processors,int solver_type)	
{
	if (A.mounted == false)
		A.Mount();
	
	//Parametros - PARDISO
	int			mtype = 11;						//Tipo de solução -- Real unsymmetric matrix
	int			nrhs = 1;						//Numero de colunas do lado direito da equacao
	void*		pt[64];							//Ponteiro de memoria interno do PARDISO
	int			iparm[64];						//Parametros de controle do PARDISO
	int			maxfct = 1;						//Numero maximo de fatorizacoes numericas salvas na memória
	int			mnum = 1;						//Qual fatorizacao utilizar
	int			phase;							//Fase do processo de solução
	int			msglvl = 0;						//Imprime informações estatisticas
	int			solver;							//0 - Esparso 
												//1 - Iterativo
	if (solver_type != 0 && solver_type != 1)
	{
		printf("Bad solver choice. Sparse Solver is adopted\n");
		solver = 0;
	}
	else
	{
		solver = solver_type;
	}
	//Variaveis auxiliares
	double		ddum;										//Variavel dummy
	int			idum;										//Variavel dummy

	double*		a = A.m_matrix.valuePtr();					//Valores da matriz
	int*		ia = A.m_matrix.outerIndexPtr();			//indice dos valores "outer" da matriz
	int*		ja = A.m_matrix.innerIndexPtr();			//indice do primeiro valor da linha da matriz

	int			n = (int)A.m_matrix.rows();						//Numero de linhas da matriz
	int			nnz = ia[n];								//Numero de valores na matriz
	Matrix		x(n);										//Matriz de retorno - solução
	if (n == 0)
		return n;
	///////////////////////////////////////////////////////////////////////////
	//     Inicializa parametros de controle                                 //
	///////////////////////////////////////////////////////////////////////////
	pardisoinit(pt, &mtype, iparm);
	//Setando algumas variaveis especificas de controle
	iparm[1] = 3;
	if (solver == 1)//solve iterativo
		iparm[3] = 31;
	///////////////////////////////////////////////////////////////////////////
	//     Adiciona 1 ao indice - 1-index (FORTRAN)                          //
	///////////////////////////////////////////////////////////////////////////
	for (size_t i = 0, iLen = n + 1; i < iLen; i++) 
	{
		ia[i] += 1;
	}

	for (size_t i = 0, iLen = nnz; i < iLen; i++) 
	{
		ja[i] += 1;
	}
	///////////////////////////////////////////////////////////////////////////
	//     Fase de analise, fatorização numerica, solução e refinamento      //
	//     iterativo                                                         //
	///////////////////////////////////////////////////////////////////////////
	phase = 13;

	pardiso(pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, b.getMatrix(), x.getMatrix(), info_fail);

	if (*info_fail != 0) 
	{
		std::cout << "ERROR during analysis, numerical factorization, solver, iterative refinement: " << *info_fail << std::endl;
	}
	///////////////////////////////////////////////////////////////////////////
	//     Subtrai 1 do indice - 0-index (C++)                               //
	///////////////////////////////////////////////////////////////////////////
	for (size_t i = 0, iLen = n + 1; i < iLen; i++) 
	{
		ia[i] -= 1;
	}

	for (size_t i = 0, iLen = nnz; i < iLen; i++) 
	{
		ja[i] -= 1;
	}
	///////////////////////////////////////////////////////////////////////////
	//     Fase de limpeza de memória                                        //
	///////////////////////////////////////////////////////////////////////////
	phase = -1;
	int temp_fail = *info_fail;
	pardiso(pt, &maxfct, &mnum, &mtype, &phase, &n, &ddum, ia, ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, info_fail);
	if (temp_fail != 0 || *info_fail != 0)
		*info_fail = 1; //falha

	return x;
}

//Operador Multiplicacao de matrizes (matrix2 deve necessariamente ser um vetor)
Matrix operator * (SparseMatrix &matrix1, Matrix &matrix2)
{
	if (matrix1.mounted == false)
		matrix1.Mount();
	//Verificação da possibilidade de multiplicação
	if (matrix1.m_matrix.cols() != matrix2.getLines())
	{
		printf("Impossible to multiply matrices. Dimensions are not compatible\n");
		return NULL;
	}
	else
	{
		Matrix return_m((long)matrix1.m_matrix.rows(), 1);

		double*		a = matrix1.m_matrix.valuePtr();					//Valores da matriz
		int*		ia = matrix1.m_matrix.outerIndexPtr();			//indice dos valores das colunas da matriz
		int*		ja = matrix1.m_matrix.innerIndexPtr();			//indice do primeiro valor da linha da matriz
		//Multiplicação
		for (int i = 0; i < matrix1.m_matrix.rows(); i++)
		{
			for (int pos = ia[i]; pos < ia[i + 1]; pos++)
				return_m(i, 0) += a[pos] * matrix2(ja[pos],0);
		}
		return return_m;
	}
}

//Operador Multiplicacao de matriz esparsa por escalar
SparseMatrix operator * (double sigma, SparseMatrix &matrix1)
{
	SparseMatrix return_m;
	if (matrix1.mounted == false)
		matrix1.Mount();
	return_m.rows = matrix1.rows;
	return_m.cols = matrix1.cols;
	return_m.m_matrix = sigma*matrix1.m_matrix;
	return_m.non_null_estimative = matrix1.non_null_estimative;
	return_m.mounted = matrix1.mounted;
	return_m.tripletList = matrix1.tripletList;
	return return_m;
}

//Operador Soma de matrizes esparsas
SparseMatrix operator + (SparseMatrix& matrix1, SparseMatrix& matrix2)
{
	SparseMatrix return_m;
	if (matrix1.mounted == false)
		matrix1.Mount();
	if (matrix2.mounted == false)
		matrix2.Mount();
	//Verificação da possibilidade de soma
	if (matrix1.m_matrix.rows() != matrix2.m_matrix.rows() || matrix1.m_matrix.cols() != matrix2.m_matrix.cols())
	{
		printf("Impossible to sum matrices. Dimensions are not compatible\n");
		return return_m;
	}
	else
	{
		return_m.rows = matrix1.rows;
		return_m.cols = matrix1.cols;
		return_m.mounted = matrix1.mounted;
		return_m.non_null_estimative = matrix1.non_null_estimative;
		return_m.m_matrix = matrix1.m_matrix + matrix2.m_matrix;
		return return_m;
	}
}

//Função de escrita da matriz em arquivo de texto
void SparseMatrix::WriteMatrix(char* name)
{
	FILE *f = fopen(name, "w");
	//Loop para percorrer as linhas
	fprintf(f, "rows\t%d\n",(int)m_matrix.rows());
	fprintf(f, "columns\t%d\n",(int)m_matrix.cols());
	fprintf(f, "row\t");
	fprintf(f, "column\t");
	fprintf(f, "value\n");
	for (int k = 0; k<m_matrix.outerSize(); ++k)
	for (Eigen::SparseMatrix<double, 1, int>::InnerIterator it(m_matrix, k); it; ++it)
	{
		fprintf(f, "%d\t", (int) it.row() + 1);
		fprintf(f, "%d\t", (int) it.col() + 1);
		fprintf(f, "%.20e\n", it.value());
	}
	fclose(f);
}

//Operador de atribuição
SparseMatrix &SparseMatrix::operator = (SparseMatrix const &matrix1)
{
	rows = matrix1.rows;
	cols = matrix1.cols;
	m_matrix = matrix1.m_matrix;
	non_null_estimative = matrix1.non_null_estimative;
	mounted = matrix1.mounted;
	tripletList = matrix1.tripletList;
	//Retorna esta matriz
	return *this;
}

//Calcula freq. naturais usando código ARPACK - Shifted inverse mode
Matrix sparseeigen(SparseMatrix &K, SparseMatrix &M, Matrix &z, int n_e, bool eigenvectors, int* ret, double tolerance)
{
	*ret = 0;
	int ido = 0;				//Parametro de controle
	char bmat = 'G';			//Autovalor generalizado
	int n = K.rows;				//Ordem do problema
	char which[2];
	which[0] = 'L';				//Larger module parts of the eigenvalues (por conta do shift inverse,depois serão transformados para os menores)
	which[1] = 'M';
	int nev = n_e;				//Numero de autovalores requeridos
	double tol = tolerance;		//Tolerancia - valor nulo para ativar o calculo interno da precisão de maquina, via Lapack
	Matrix resid(n);			//Residuo
	int ncv = 2 * n_e + 1;		//Dimensão da base utilizada para o calculo aproximado
	Matrix v(n, ncv);
	int ldv = n;
	int* iparam;
	iparam = new int[11];
	iparam[0] = 1;				//Shift mode
	iparam[2] = 300;			//Max iterations
	iparam[6] = 3;				//Mode
	int* ipntr;
	ipntr = new int[14];		//Ponteiro para saidas
	Matrix workd(3 * n);
	int lworkl = 3 * ncv*ncv + 18 * ncv;
	Matrix workl(lworkl);
	int info = 0;				//Indica para o ARPACK realizar um chute inicial aleatório

	double sigma = 0;			//Shift value - manter zero, para capturar os menores autovalores do problema (transformação espectral de interesse)
	double sigmar = sigma;
	double sigmai = 0;

	Matrix return_m(n_e, 2);

	//Verificações de tamanhos de matrizes
	if (ncv > n)
	{
		*ret = -1;
		printf("The requested number of eigenvalues is too large. Decrease it.\n");
		return return_m;
	}
		

	SparseMatrix C = K + ((-sigma)*M);
	
	Matrix temp1(n);
	Matrix Awd(n);
	///LOOP chamando a função dnaupd do ARPACK////
	while (true)
	{
		arpack::internal::dnaupd_c(&ido, &bmat, n, which, nev, tol, resid.getMatrix(),ncv, v.getMatrix(), ldv, iparam, ipntr, workd.getMatrix(), workl.getMatrix(), lworkl, &info);
		//printf("Iteration \t %lf\n",norm(resid));//Plota norma do residuo
		if (ido == -1)
		{
			//printf("IDO = -1\n");
			for (int i = 0; i<n; i++)
				temp1(i, 0) = workd(ipntr[0] - 1 + i, 0);
			Awd = M*temp1;
			int info_fail = 0;
			Awd = sparsesystem(C, Awd, &info_fail, db.solver_options->processors, db.solver_options->solver);	//Resolve o sistema linear
			for (int i = 0; i<n; i++)
				workd(ipntr[1] - 1 + i, 0) = Awd(i, 0);
		}
		else if (ido == 1)
		{
			//printf("IDO = 1\n");
			for (int i = 0; i<n; i++)
				temp1(i, 0) = workd(ipntr[2] - 1 + i, 0);
			int info_fail = 0;
			temp1 = sparsesystem(C, temp1, &info_fail, 1, 0);	//Resolve o sistema linear
			for (int i = 0; i<n; i++)
				workd(ipntr[1] - 1 + i, 0) = temp1(i, 0);
		}
		else if (ido == 2)
		{
			//printf("IDO = 2\n");
			for (int i = 0; i<n; i++)
				temp1(i, 0) = workd(ipntr[0] - 1 + i, 0);
			Awd = M*temp1;
			for (int i = 0; i<n; i++)
				workd(ipntr[1] - 1 + i, 0) = Awd(i, 0);
		}
		else if (info < 0)
		{
			printf("Error with dnaupd in ARPACK method, info = %d.\n", info);
			*ret = -1;
			return return_m;
		}
		else//Convergência ocorreu - saida do loop
			break;
	}//end of while

	//Pós - processamento do ARPACK	
	//Os autovetores serão armazenados em z (parametro de entrada da função), se requeridos
	int rvec = int (eigenvectors);		//Computar autovetores ou não
	char Howmny = 'A';
	int *select;
	select = new int[ncv];
	Matrix dr(nev + 1);		//parte real dos autovalores
	Matrix di(nev + 1);		//parte imaginaria dos autovalores
	Matrix res_calc(nev + 1);//para salvar o residuo calculado dos autovetores
	int ldz = n;
	Matrix Workev(3 * ncv);
	arpack::internal::dneupd_c(rvec, &Howmny, select, dr.getMatrix(), di.getMatrix(), z.getMatrix(),ldz, sigmar, sigmai, Workev.getMatrix(), &bmat, n, which, nev,
			tol, resid.getMatrix(), ncv, v.getMatrix(), ldv, iparam, ipntr, workd.getMatrix(), workl.getMatrix(), lworkl, &info);
	
	if (info != 0)
	{
		printf("Error with dneupd in ARPACK method, info = %d.\n", info);
		*ret = -1;
	}
	else
	{
		for (int i = 0; i < n_e; i++)
		{
			return_m(i, 0) = dr(i, 0);
			return_m(i, 1) = di(i, 0);
		}
		//Autovetores
		if (eigenvectors == true)
		{
			Matrix xr(n);
			Matrix xi(n);
			//Percorre cada um dos autovetores calculados para calculo do erro cometido
			bool first = true;
			for (int j = 0; j < iparam[4]; j++)
			{
				//Autovalor real
				if (di(j, 0) == 0.0)
				{
					for (int cp = 0; cp < n; cp++)
						xr(cp, 0) = z(cp, j);
					Awd = K*xr - dr(j, 0)*(M*xr);
					res_calc(j, 0) = norm(Awd)/abs(dr(j,0));
				}
				else//Autovalor complexo
				{
					if (first)//se for o primeiro dos autovetores (o próximo sera complexo conjugado)
					{
						for (int cp = 0; cp < n; cp++)
						{
							xr(cp, 0) = z(cp, j);		//parte real
							xi(cp, 0) = z(cp, j + 1);	//parte imaginaria
						}
						//parte real
						Awd = K*xr - dr(j, 0)*(M*xr) + di(j, 0)*(M*xi);
						res_calc(j, 0) = norm(Awd);
						//parte imaginaria
						Awd = K*xi - dr(j, 0)*(M*xi) - di(j, 0)*(M*xr);
						double normi = norm(Awd);
						res_calc(j, 0) = sqrt(res_calc(j, 0)*res_calc(j, 0) + normi*normi);
						res_calc(j, 0) = res_calc(j, 0) / sqrt(dr(j,0)*dr(j,0)+di(j,0)*di(j,0));
						res_calc(j + 1, 0) = res_calc(j, 0);
						first = false;
					}
					else
						first = true;
				}
			}
		}
	}

	//Check final de execução
	double res = norm(res_calc);
	if (res > tolerance)
		*ret = -2;
	
	//Autovalores e residuos calculados
	db.myprintf("\nARPACK output:\n");
	for (int i = 0; i < n_e; i++)
		db.myprintf("%.6f + %.6f i -> residual %.2e\n", dr(i, 0), di(i, 0), res_calc(i,0));
	db.myprintf("\nConverged Eigenvalues: %d\n", iparam[4]);

	delete[] iparam;
	delete[] ipntr;
	delete[] select;	
	return return_m;

	//Valores de retorno
	//-1 - ocorreu erro durante a execução do dneupd
	//-2 - chegou a pós-processar, mas não houve convergência com a tolerancia especificada
	// 0 - execução correta
}
