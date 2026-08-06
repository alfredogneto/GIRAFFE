#include "SparseMatrix.h"
#include <iostream>
#ifdef GIRAFFE_USE_MKL
#include <mkl.h>
#else
#include <Eigen/SparseLU>
#endif

#include <arpack.hpp>
#include <debug_c.hpp>  // debug arpack.
#include <stat_c.hpp>   // arpack statistics.

#include "Matrix.h"
#include "SolverOptions.h"
#include "Database.h"
//Global variables
extern
Database db;

SparseMatrix::SparseMatrix()
{
	mounted = false;
	rows = 1;
	cols = 1;
	non_null_estimative = 1;
	m_matrix.resize(1, 1);				//size of the matrix
	tripletList.reserve(1);				//pre-allocates the triplet list
}
SparseMatrix::SparseMatrix(int e_rows,int e_cols, int e_non_null_estimative)
{
	rows = e_rows;
	cols = e_cols;
	non_null_estimative = e_non_null_estimative;
	m_matrix.resize(e_rows, e_cols);					//size of the matrix
	m_matrix.reserve(e_non_null_estimative);
	tripletList.reserve(2*e_non_null_estimative);		//pre-allocates the triplet list
	mounted = false;
}
SparseMatrix::SparseMatrix(const SparseMatrix &copied)
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
	//N�o there is variables the to be deallocated manualmente (all s�o automaticas)
}

//Zeroes coefficients of the sparse matrix (allocated in the tripletList)
void SparseMatrix::Clear()
{
	tripletList.erase(tripletList.begin(), tripletList.end());
	mounted = false;
}

//Sets in the posi��o (i,j) of the square matrix the value v
void SparseMatrix::setValue(int i, int j, double v)
{
	if (i + 1 > m_matrix.rows() || j + 1 > m_matrix.cols())
		printf("Error assigning sparse matrix value\n");
	else
		tripletList.push_back(T(i, j, v));
	if (mounted)//if there was inser��o ap�s assembly, sets assembly as false
		mounted = false;
}
void SparseMatrix::Mount() const
{
	m_matrix.setFromTriplets(tripletList.begin(), tripletList.end());		//prepara��o of the matrix - assembly the partir of the triplet list
	mounted = true;
}

//Solves the linear system of the way Ax=b
Matrix sparsesystem(SparseMatrix &A, Matrix &b, int *info_fail,int processors,int solver_type)
{
	if (A.mounted == false)
		A.Mount();

	int			n = (int)A.m_matrix.rows();						//Number of rows of the matrix
	Matrix		x(n);										//Matrix of return - solu��o
	if (n == 0)
		return n;

#ifdef GIRAFFE_USE_MKL
	//Parameters - PARDISO
	int			mtype = 11;						//Type of solu��o -- Real unsymmetric matrix
	int			nrhs = 1;						//Number of columns of the lado direito of the equacao
	void*		pt[64];							//Pointer of memory internal of the PARDISO
	int			iparm[64];						//Parameters of control of the PARDISO
	int			maxfct = 1;						//Maximum number of fatorizacoes numericas saved in the mem�ria
	int			mnum = 1;						//Qual fatorizacao use
	int			phase;							//Fase of the processo of solu��o
	int			msglvl = 0;						//Prints informa��es estatisticas
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
	//Auxiliary variables
	double		ddum;										//Variable dummy
	int			idum;										//Variable dummy

	double*		a = A.m_matrix.valuePtr();					//Values of the matrix
	int*		ia = A.m_matrix.outerIndexPtr();			//indice of the values "outer" of the matrix
	int*		ja = A.m_matrix.innerIndexPtr();			//indice of the first value of the linha of the matrix

	int			nnz = ia[n];								//Number of values in the matrix
	///////////////////////////////////////////////////////////////////////////
	//     Initializes parameters of control                                 //
	///////////////////////////////////////////////////////////////////////////
	pardisoinit(pt, &mtype, iparm);
	//Setting some variables especificas of control
	iparm[1] = 3;
	if (solver == 1)//solve iterativo
		iparm[3] = 31;
	///////////////////////////////////////////////////////////////////////////
	//     Adds 1 to the indice - 1-index (FORTRAN)                          //
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
	//     Fase of analysis, fatoriza��o numerica, solu��o and refinamento      //
	//     iterativo                                                         //
	///////////////////////////////////////////////////////////////////////////
	phase = 13;

	pardiso(pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, b.getMatrix(), x.getMatrix(), info_fail);

	if (*info_fail != 0)
	{
		std::cout << "ERROR during analysis, numerical factorization, solver, iterative refinement: " << *info_fail << std::endl;
	}
	///////////////////////////////////////////////////////////////////////////
	//     Subtrai 1 of the indice - 0-index (C++)                               //
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
	//     Fase of limpeza of mem�ria                                        //
	///////////////////////////////////////////////////////////////////////////
	phase = -1;
	int temp_fail = *info_fail;
	pardiso(pt, &maxfct, &mnum, &mtype, &phase, &n, &ddum, ia, ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, info_fail);
	if (temp_fail != 0 || *info_fail != 0)
		*info_fail = 1; //falha
#else
	//Non-MKL platforms (macOS/Linux): no PARDISO available, so solve the sparse
	//system with Eigen's SparseLU instead. Note: this always does a direct
	//solve - the PARDISO "Iterative" solver_type option has no equivalent here.
	(void)processors;
	(void)solver_type;
	Eigen::SparseMatrix<double, Eigen::ColMajor, int> A_colmajor = A.m_matrix;
	Eigen::SparseLU<Eigen::SparseMatrix<double, Eigen::ColMajor, int>> solver;
	solver.compute(A_colmajor);
	if (solver.info() != Eigen::Success)
	{
		std::cout << "ERROR during sparse LU factorization (Eigen SparseLU)." << std::endl;
		*info_fail = 1;
		return x;
	}
	Eigen::Map<Eigen::VectorXd> b_vec(b.getMatrix(), n);
	Eigen::VectorXd x_vec = solver.solve(b_vec);
	if (solver.info() != Eigen::Success)
	{
		std::cout << "ERROR during sparse solve (Eigen SparseLU)." << std::endl;
		*info_fail = 1;
		return x;
	}
	for (int i = 0; i < n; i++)
		x(i, 0) = x_vec(i);
	*info_fail = 0;
#endif

	return x;
}

//Operator Multiplication of matrices (matrix2 must necessariamente be a vector)
Matrix operator * (const SparseMatrix &matrix1, const Matrix &matrix2)
{
	if (matrix1.mounted == false)
		matrix1.Mount();
	//Verifica��o of the possibility of multiplica��o
	if (matrix1.m_matrix.cols() != matrix2.getLines())
	{
		printf("Impossible to multiply matrices. Dimensions are not compatible\n");
		return NULL;
	}
	else
	{
		Matrix return_m((long)matrix1.m_matrix.rows(), 1);

		double*		a = matrix1.m_matrix.valuePtr();					//Values of the matrix
		int*		ia = matrix1.m_matrix.outerIndexPtr();			//indice of the values of the columns of the matrix
		int*		ja = matrix1.m_matrix.innerIndexPtr();			//indice of the first value of the linha of the matrix
		//Multiplica��o
		for (int i = 0; i < matrix1.m_matrix.rows(); i++)
		{
			for (int pos = ia[i]; pos < ia[i + 1]; pos++)
				return_m(i, 0) += a[pos] * matrix2(ja[pos],0);
		}
		return return_m;
	}
}

//Operator Multiplication of sparse matrix by scalar
SparseMatrix operator * (double sigma, const SparseMatrix &matrix1)
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

//Operator Sum of matrices esparsas
SparseMatrix operator + (const SparseMatrix& matrix1, const SparseMatrix& matrix2)
{
	SparseMatrix return_m;
	if (matrix1.mounted == false)
		matrix1.Mount();
	if (matrix2.mounted == false)
		matrix2.Mount();
	//Verifica��o of the possibility of sum
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

//Fun��o of matrix writing in file of text
void SparseMatrix::WriteMatrix(char* name)
{
	FILE *f = fopen(name, "w");
	//Loop for iterate over the rows
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

//Operator of atribui��o
SparseMatrix &SparseMatrix::operator = (SparseMatrix const &matrix1)
{
	rows = matrix1.rows;
	cols = matrix1.cols;
	m_matrix = matrix1.m_matrix;
	non_null_estimative = matrix1.non_null_estimative;
	mounted = matrix1.mounted;
	tripletList = matrix1.tripletList;
	//Returns this matrix
	return *this;
}

//Calculates freq. naturais using c�digo ARPACK - Shifted inverse mode
Matrix sparseeigen(SparseMatrix &K, SparseMatrix &M, Matrix &z, int n_e, bool eigenvectors, int* ret, double tolerance)
{
	*ret = 0;
	int ido = 0;				//Parameter of control
	char bmat = 'G';			//Eigenvalue generalizado
	int n = K.rows;				//Order of the problem
	char which[2];
	which[0] = 'L';				//Larger module parts of the eigenvalues (by conta of the shift inverse,depois ser�o transformados for the menores)
	which[1] = 'M';
	int nev = n_e;				//Number of eigenvalues requeridos
	double tol = tolerance;		//Tolerance - value null for activate the calculation internal of the precis�o of maquina, via Lapack
	Matrix resid(n);			//Residual
	int ncv = 2 * n_e + 1;		//Dimens�o of the base used for the calculation aproximado
	Matrix v(n, ncv);
	int ldv = n;
	int* iparam;
	iparam = new int[11];
	iparam[0] = 1;				//Shift mode
	iparam[2] = 300;			//Max iterations
	iparam[6] = 3;				//Mode
	int* ipntr;
	ipntr = new int[14];		//Pointer for saidas
	Matrix workd(3 * n);
	int lworkl = 3 * ncv*ncv + 18 * ncv;
	Matrix workl(lworkl);
	int info = 0;				//Indicates for the ARPACK perform a initial guess aleat�rio

	double sigma = 0;			//Shift value - keep zero, for capturar the menores eigenvalues of the problem (transforma��o espectral of interest)
	double sigmar = sigma;
	double sigmai = 0;

	Matrix return_m(n_e, 2);

	//Verifica��es of tamanhos of matrices
	if (ncv > n)
	{
		*ret = -1;
		printf("The requested number of eigenvalues is too large. Decrease it.\n");
		return return_m;
	}
		

	SparseMatrix C = K + ((-sigma)*M);

#ifndef GIRAFFE_USE_MKL
	// C never changes across ARPACK's shift-invert iterations - only the
	// right-hand side does - so factorize it once here instead of letting
	// every ido==-1/ido==1 solve below re-factorize it from scratch via
	// sparsesystem(). With PARDISO (Windows/MKL) a from-scratch factorize is
	// fast enough not to matter; with Eigen::SparseLU it's the dominant cost
	// of this whole function (one call could take tens of seconds on a model
	// with a few thousand DOFs and ~300 ARPACK iterations).
	if (C.mounted == false)
		C.Mount();
	Eigen::SparseMatrix<double, Eigen::ColMajor, int> C_colmajor = C.m_matrix;
	Eigen::SparseLU<Eigen::SparseMatrix<double, Eigen::ColMajor, int>> C_solver;
	C_solver.compute(C_colmajor);
	if (C_solver.info() != Eigen::Success)
	{
		printf("Error during sparse LU factorization (Eigen SparseLU) in sparseeigen.\n");
		*ret = -1;
		return return_m;
	}
#endif

	Matrix temp1(n);
	Matrix Awd(n);
	///LOOP chamando a fun��o dnaupd do ARPACK////
	while (true)
	{
		arpack::internal::dnaupd_c(&ido, &bmat, n, which, nev, tol, resid.getMatrix(),ncv, v.getMatrix(), ldv, iparam, ipntr, workd.getMatrix(), workl.getMatrix(), lworkl, &info);
		//printf("Iteration \t %lf\n",norm(resid));//Plots norm of the residual
		if (ido == -1)
		{
			//printf("IDO = -1\n");
			for (int i = 0; i<n; i++)
				temp1(i, 0) = workd(ipntr[0] - 1 + i, 0);
			Awd = M*temp1;
#ifdef GIRAFFE_USE_MKL
			int info_fail = 0;
			Awd = sparsesystem(C, Awd, &info_fail, db.solver_options->processors, db.solver_options->solver);	//Solves the linear system
#else
			Eigen::Map<Eigen::VectorXd> Awd_vec(Awd.getMatrix(), n);
			Awd_vec = C_solver.solve(Awd_vec);
#endif
			for (int i = 0; i<n; i++)
				workd(ipntr[1] - 1 + i, 0) = Awd(i, 0);
		}
		else if (ido == 1)
		{
			//printf("IDO = 1\n");
			for (int i = 0; i<n; i++)
				temp1(i, 0) = workd(ipntr[2] - 1 + i, 0);
#ifdef GIRAFFE_USE_MKL
			int info_fail = 0;
			temp1 = sparsesystem(C, temp1, &info_fail, 1, 0);	//Solves the linear system
#else
			Eigen::Map<Eigen::VectorXd> temp1_vec(temp1.getMatrix(), n);
			temp1_vec = C_solver.solve(temp1_vec);
#endif
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
		else//Converg�ncia ocorreu - output of the loop
			break;
	}//end of while

	//P�s - processamento do ARPACK	
	//The eigenvectors ser�o armazenados in z (parameter of input of the fun��o), if requeridos
	int rvec = int (eigenvectors);		//Compute eigenvectors or n�o
	char Howmny = 'A';
	int *select;
	select = new int[ncv];
	Matrix dr(nev + 1);		//part real of the eigenvalues
	Matrix di(nev + 1);		//part imaginaria of the eigenvalues
	Matrix res_calc(nev + 1);//for save the residual calculated of the eigenvectors
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
		//Eigenvectors
		if (eigenvectors == true)
		{
			Matrix xr(n);
			Matrix xi(n);
			//Iterates over each a of the eigenvectors calculated for calculation of the error cometido
			bool first = true;
			for (int j = 0; j < iparam[4]; j++)
			{
				//Eigenvalue real
				if (di(j, 0) == 0.0)
				{
					for (int cp = 0; cp < n; cp++)
						xr(cp, 0) = z(cp, j);
					Awd = K*xr - dr(j, 0)*(M*xr);
					res_calc(j, 0) = norm(Awd)/abs(dr(j,0));
				}
				else//Eigenvalue complexo
				{
					if (first)//if for the first of the eigenvectors (the pr�ximo will be complexo conjugado)
					{
						for (int cp = 0; cp < n; cp++)
						{
							xr(cp, 0) = z(cp, j);		//part real
							xi(cp, 0) = z(cp, j + 1);	//part imaginaria
						}
						//part real
						Awd = K*xr - dr(j, 0)*(M*xr) + di(j, 0)*(M*xi);
						res_calc(j, 0) = norm(Awd);
						//part imaginaria
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

	//Check final of execu��o
	double res = norm(res_calc);
	if (res > tolerance)
		*ret = -2;
	
	//Eigenvalues and residuals calculated
	db.myprintf("\nARPACK output:\n");
	for (int i = 0; i < n_e; i++)
		db.myprintf("%.6f + %.6f i -> residual %.2e\n", dr(i, 0), di(i, 0), res_calc(i,0));
	db.myprintf("\nConverged Eigenvalues: %d\n", iparam[4]);

	delete[] iparam;
	delete[] ipntr;
	delete[] select;	
	return return_m;

	//Values of return
	//-1 - ocorreu error durante the execu��o of the dneupd
	//-2 - chegou the p�s-processar, but n�o there was converg�ncia with the tolerance specified
	// 0 - execu��o correct
}
