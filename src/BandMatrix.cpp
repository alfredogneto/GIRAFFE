#include "BandMatrix.h"
#ifdef GIRAFFE_USE_MKL
#include <mkl.h>
#else
// Non-MKL LAPACK (OpenBLAS/Accelerate/reference LAPACK) exports the classic
// Fortran-mangled names (trailing underscore) instead of MKL's bare C aliases.
extern "C" {
	void dgbtrf_(int *m, int *n, int *kl, int *ku, double *ab, int *ldab, int *ipiv, int *info);
	void dgbtrs_(char *trans, int *n, int *kl, int *ku, int *nrhs, double *ab, int *ldab, int *ipiv, double *b, int *ldb, int *info);
}
static inline void dgbtrf(int *m, int *n, int *kl, int *ku, double *ab, int *ldab, int *ipiv, int *info)
{
	dgbtrf_(m, n, kl, ku, ab, ldab, ipiv, info);
}
static inline void dgbtrs(char *trans, int *n, int *kl, int *ku, int *nrhs, double *ab, int *ldab, int *ipiv, double *b, int *ldb, int *info)
{
	dgbtrs_(trans, n, kl, ku, nrhs, ab, ldab, ipiv, b, ldb, info);
}
#endif

//Standard constructor
BandMatrix::BandMatrix(void)
{
	m_lines_deleted = true;
	m_alloced_lines = 0;
	m_columns = 1;
	m_band_s = 0;
	m_band_i = 0;
	//Initializes the matrix as null and tries allocates-it
	m_matrix = NULL;
	if(!alloc())
		printf("Not possible to alloc matrix\n");
}
//Parametric constructor
BandMatrix::BandMatrix(long band_s,long band_i,long columns)
{
	m_lines_deleted = true;
	m_alloced_lines = 0;
	m_columns = columns;
	m_band_s = band_s;
	m_band_i = band_i;
	//Initializes the matrix as null and tries allocates-it
	m_matrix = NULL;
	if(!alloc())
		printf("Not possible to alloc matrix\n");
}
//Constructor of copy
BandMatrix::BandMatrix(BandMatrix &copied)
{
	//Checks dimensions of the matrix - if necessary, does re-alocation
	//if (copied.m_alloced_lines != m_alloced_lines)
//	{
		flush();
		m_alloced_lines = 0;
		m_columns = copied.m_columns;
		m_band_s = copied.m_band_s;
		m_band_i = copied.m_band_i;
		//Initializes the matrix as null and tries allocates-it
		m_matrix = NULL;
		if(!alloc())
			printf("Not possible to alloc matrix\n");
//	}
	//Copies the values in the new matrix
	for(long i=0; i < copied.m_alloced_lines; i++)
		m_matrix[i] = copied.m_matrix[i];	
}
BandMatrix::~BandMatrix(void)
{
	flush();
}
//Prints the matrix in the tela
void BandMatrix::print()
{
	printf("\n");
	for(long i=0; i < (2*m_band_i+m_band_s+1); i++)
	{
		printf("|");
		for(long j=0; j < m_columns; j++)
		{
			if(m_matrix[i+j*(2*m_band_i+m_band_s+1)] >= 0.0)
				printf(" %.2f\t",m_matrix[i+j*(2*m_band_i+m_band_s+1)]);
			else
				printf(" %.1f\t",m_matrix[i+j*(2*m_band_i+m_band_s+1)]);
		}
		printf("| \n");
	}
	printf("\n");
}
//Prints the matrix in a file of text, whose name this in the char s
void BandMatrix::fprint(char* s)
{
	FILE *file1 = fopen(s,"w");
	
	for(long i=0; i < (2*m_band_i+m_band_s+1); i++)
	{
		for(long j=0; j < m_columns; j++)
		{
			if(m_matrix[i+j*(2*m_band_i+m_band_s+1)] >= 0.0)
				fprintf(file1," %e\t",m_matrix[i+j*(2*m_band_i+m_band_s+1)]);
			else
				fprintf(file1," %e\t",m_matrix[i+j*(2*m_band_i+m_band_s+1)]);
		}
		fprintf(file1,"\n");
	}
	
	fclose(file1);
}

//Prints the matrix in a file of text, whose name this in the char s (in formato quadrada)
void BandMatrix::fprintSquare(char* s)
{
	FILE *file1 = fopen(s,"w");
	for(long i=0; i < m_columns; i++)
	{
		for(long j=0; j < m_columns; j++)
		{
			//If this in the regi�o of values n�o null
			if(i-j <= m_band_i && j-i <= m_band_s)
				fprintf(file1," %.20e\t",this->getValue(i,j));
			else
				fprintf(file1," 0\t");
		}
		fprintf(file1,"\n");
	}
	fclose(file1);
}


//Allocates the matrix
bool BandMatrix::alloc()
{
	flush();
	//Tries allocate the matrix
	if(!(m_matrix = new double[m_columns*(2*m_band_i+m_band_s+1)]))
			//Falha
			return false;
	//Case consiga
	else
	{
		m_alloced_lines = m_columns*(2*m_band_i+m_band_s+1);
		m_lines_deleted = false;
		
		for(long i=0; i < m_alloced_lines; i++)
				m_matrix[i] = 0.0;
		return true;
	}	
}
//Frees the mem�ria occupied through the matrix
bool BandMatrix::flush()
{
	//Iterates over the rows
	//If there is rows for deallocate
	if (m_lines_deleted == false)
	{
		delete[]m_matrix;
		//Desalocada with sucesso
		m_lines_deleted = true;	//for avoid desaloca��o novamente
		m_matrix = NULL;
		m_alloced_lines = 0;
		//m_columns=0;
		//m_band_i=0;
		//m_band_s=0;
	}
	return true;
}
//Zeroes the matrix, keeping the dimens�es current
void BandMatrix::clear()
{
	for(int i=0; i < m_alloced_lines; i++)
			m_matrix[i] = 0.0;
}
//Returns the band lower of the matrix
int BandMatrix::getBandi() const
{
	return m_band_i;
}
//Returns the band upper of the matrix
int BandMatrix::getBands() const
{
	return m_band_s;
}
//Returns the number of rows of the matrix
int BandMatrix::getColumns() const
{
	return m_columns;
}
//Sets the band of the matrix
void BandMatrix::setBandi(int bandi)
{
	m_band_i = bandi;
}
//Sets the band of the matrix
void BandMatrix::setBands(int bands)
{
	m_band_s = bands;
}
//Sets the number of columns of the matrix
void BandMatrix::setColumns(int columns)
{	
	m_columns = columns;
}
//Return of the value in the posi��o specified - input i and j in the posi��o of the square matrix original
double BandMatrix::getValue(long i, long j) const
{
	long column,line;
	bool flag_write = true;
	//Tests for ver if the posi��o this inside of the band allocated
	if (i>j)//part lower of the diagonal
	{
		if (i-j > m_band_i)
		{
			printf("Invalid position at band matrix\n");
			flag_write = false;
		}
	}
	else
	{
		if (j-i > m_band_s)
		{
			printf("Invalid position at band matrix\n");
			flag_write = false;
		}
	}
	
	if ( flag_write == false )
	{
		return 0;
	}
	else
	{
		line = m_band_i + m_band_s + i - j;
		column = j;
		
		return m_matrix[line + column*(2*m_band_i+m_band_s+1)];
	}
}
//Sets in the posi��o (i,j) of the square matrix the value v
void BandMatrix::setValue(long i, long j, double v)
{
	long column,line;

	bool flag_write = true;
	//Tests for ver if the posi��o this inside of the band allocated
	if (i>j)//part lower of the diagonal
	{
		if (i-j > m_band_i)
		{
			printf("Invalid position at band matrix\n");
			flag_write = false;
		}
	}
	else
	{
		if (j-i > m_band_s)
		{
			printf("Invalid position at band matrix\n");
			flag_write = false;
		}
	}
	
	if ( flag_write == true )
	{
		line = m_band_i + m_band_s + i - j;
		column = j;
		m_matrix[line + column*(2*m_band_i+m_band_s+1)] = v;
	}
}
//Returns the endere�o of a matrix
double* BandMatrix::getMatrix() const
{
	return m_matrix;
}
//Solves the linear system of the way Ax=b (being A a band matrix)
Matrix bandsystem(BandMatrix &A, Matrix &b,int *info_fail)
{
	int n = A.getColumns();
	int kl = A.getBandi();
	int ku = A.getBands();
	char trans = 'N';
	int nrhs = 1;
	int ldab = 2*A.getBandi() + A.getBands() + 1;
	int info = 1;
	int *ipiv;
	ipiv = new int[b.getLines()];
	//Does LU factorization
#ifdef GIRAFFE_USE_MKL
	MKL_Set_Num_Threads(1);
	//mkl_set_dynamic(1);
#endif
	dgbtrf(&n,&n,&kl,&ku,A.m_matrix,&ldab,ipiv,&info);
	if (info != 0)
	{
		printf("Error in LU factorization - bandsystem function\n");
		*info_fail = 1;
	}
	//Solves the linear system
	dgbtrs(&trans,&n,&kl,&ku,&nrhs,A.m_matrix,&ldab,ipiv,b.getMatrix(),&n,&info);
	if (info != 0)
	{
		printf("Error in linear system - bandsystem function\n");
		*info_fail = 1;
	}
	
	delete []ipiv;
	return b;
}
//Operator Multiplication of matrices (matrix2 must necessariamente be a vector)
Matrix operator * (const BandMatrix &matrix1, const Matrix &matrix2)
{
	//Verifica��o of the possibility of multiplica��o
	if (matrix1.getColumns() != matrix2.getLines())
	{
		printf("Impossible to multiply matrices. Dimensions are uncompatible\n");
		return NULL;
	}
	else
	{
		Matrix return_m(matrix1.getColumns(),1);
		
		//Multiplica��o
		int c_index = 0;
		for (int i = 0; i < matrix1.getColumns();i++)
		{
			for (c_index = i-matrix1.getBandi(); c_index <= i+matrix1.getBands();c_index++)
			{
				if (c_index >= 0 && c_index < matrix1.getColumns())
					return_m(i,0) +=matrix1.getValue(i,c_index)*matrix2(c_index,0);
			}
		}
		return return_m;
	}
}
//Returns the square matrix (FULL)
Matrix BandMatrix::Full()
{

	Matrix return_m(m_columns,m_columns);
	
	for(long i=0; i < m_columns; i++)
	{
		for(long j=0; j < m_columns; j++)
		{
			//If this in the regi�o of values n�o null
			if(i-j <= m_band_i && j-i <= m_band_s)
				return_m(i,j) = this->getValue(i,j);
			else
				return_m(i,j) = 0;
		}
	}
	
	return return_m;
}
//Operator Multiplication of band matrix by scalar
BandMatrix operator * (double sigma, const BandMatrix &matrix1)
{
	BandMatrix return_m(matrix1.getBands(),matrix1.getBandi(),matrix1.getColumns());
	
	//Multiplica��o
	for (int i=0;i<return_m.m_alloced_lines;i++)
		return_m.getMatrix()[i] = matrix1.getMatrix()[i]*sigma;

	return return_m;
}
//Operator sum of matrices of band
BandMatrix operator + (const BandMatrix &matrix1, const BandMatrix &matrix2)
{
	BandMatrix return_m(matrix1.getBands(),matrix1.getBandi(),matrix1.getColumns());

	//Verifica��o of the possibility of multiplica��o
	if (matrix1.getColumns() != matrix2.getColumns() || matrix1.getBandi() != matrix2.getBandi() || matrix1.getBands() != matrix2.getBands())
	{
		printf("Impossible to sum matrices. Dimensions are uncompatible\n");
	}
	else
	{
		//SUM
		for (int i=0;i<return_m.m_alloced_lines;i++)
			return_m.getMatrix()[i] = matrix1.getMatrix()[i] + matrix2.getMatrix()[i];
	}
	return return_m;
}

long BandMatrix::getAllocedLines() const
{
	return m_alloced_lines;
}
