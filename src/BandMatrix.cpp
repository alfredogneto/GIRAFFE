#include "BandMatrix.h"
#include <mkl.h>

//Construtor Padrão
BandMatrix::BandMatrix(void)
{
	m_lines_deleted = true;
	m_alloced_lines = 0;
	m_columns = 1;
	m_band_s = 0;
	m_band_i = 0;
	//Inicializa a matriz como null e tenta aloca-la
	m_matrix = NULL;
	if(!alloc())
		printf("Not possible to alloc matrix\n");
}
//Construtor Parametrico
BandMatrix::BandMatrix(long band_s,long band_i,long columns)
{
	m_lines_deleted = true;
	m_alloced_lines = 0;
	m_columns = columns;
	m_band_s = band_s;
	m_band_i = band_i;
	//Inicializa a matriz como null e tenta aloca-la
	m_matrix = NULL;
	if(!alloc())
		printf("Not possible to alloc matrix\n");
}
//Construtor de cópia
BandMatrix::BandMatrix(BandMatrix &copied)
{
	//Verifica dimensões da matriz - se necessario, faz re-alocação
	//if (copied.m_alloced_lines != m_alloced_lines)
//	{
		flush();
		m_alloced_lines = 0;
		m_columns = copied.m_columns;
		m_band_s = copied.m_band_s;
		m_band_i = copied.m_band_i;
		//Inicializa a matriz como null e tenta aloca-la
		m_matrix = NULL;
		if(!alloc())
			printf("Not possible to alloc matrix\n");
//	}
	//Copia os valores na nova matriz
	for(long i=0; i < copied.m_alloced_lines; i++)
		m_matrix[i] = copied.m_matrix[i];	
}
BandMatrix::~BandMatrix(void)
{
	flush();
}
//Imprime a matriz na tela
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
//Imprime a matriz em um arquivo de texto, cujo nome esta no char s
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

//Imprime a matriz em um arquivo de texto, cujo nome esta no char s (em formato quadrada)
void BandMatrix::fprintSquare(char* s)
{
	FILE *file1 = fopen(s,"w");
	for(long i=0; i < m_columns; i++)
	{
		for(long j=0; j < m_columns; j++)
		{
			//Se esta na região de valores não nulos
			if(i-j <= m_band_i && j-i <= m_band_s)
				fprintf(file1," %.20e\t",this->getValue(i,j));
			else
				fprintf(file1," 0\t");
		}
		fprintf(file1,"\n");
	}
	fclose(file1);
}


//Aloca a matriz
bool BandMatrix::alloc()
{
	flush();
	//Tenta alocar a matriz
	if(!(m_matrix = new double[m_columns*(2*m_band_i+m_band_s+1)]))
			//Falha
			return false;
	//Caso consiga
	else
	{
		m_alloced_lines = m_columns*(2*m_band_i+m_band_s+1);
		m_lines_deleted = false;
		
		for(long i=0; i < m_alloced_lines; i++)
				m_matrix[i] = 0.0;
		return true;
	}	
}
//Libera a memória ocupada pela matriz
bool BandMatrix::flush()
{
	//Percorre as linhas
	//Se ha linhas para desalocar
	if (m_lines_deleted == false)
	{
		delete[]m_matrix;
		//Desalocada com sucesso
		m_lines_deleted = true;	//para evitar desalocação novamente
		m_matrix = NULL;
		m_alloced_lines = 0;
		//m_columns=0;
		//m_band_i=0;
		//m_band_s=0;
	}
	return true;
}
//Zera a matriz, mantando as dimensões atuais
void BandMatrix::clear()
{
	for(int i=0; i < m_alloced_lines; i++)
			m_matrix[i] = 0.0;
}
//Retorna a banda inferior da matriz
int BandMatrix::getBandi()
{
	return m_band_i;
}
//Retorna a banda superior da matriz
int BandMatrix::getBands()
{
	return m_band_s;
}
//Retorna o numero de linhas da matriz
int BandMatrix::getColumns()
{
	return m_columns;
}
//Seta a banda da matriz
void BandMatrix::setBandi(int bandi)
{
	m_band_i = bandi;
}
//Seta a banda da matriz
void BandMatrix::setBands(int bands)
{
	m_band_s = bands;
}
//Seta o numero de colunas da matriz
void BandMatrix::setColumns(int columns)
{	
	m_columns = columns;
}
//Retorno do valor na posição especificada - entrada i e j na posição da matriz quadrada original
double BandMatrix::getValue(long i, long j)
{
	long column,line;
	bool flag_write = true;
	//Testa para ver se a posição esta dentro da banda alocada
	if (i>j)//parte inferior da diagonal
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
//Seta na posição (i,j) da matriz quadrada o valor v
void BandMatrix::setValue(long i, long j, double v)
{
	long column,line;

	bool flag_write = true;
	//Testa para ver se a posição esta dentro da banda alocada
	if (i>j)//parte inferior da diagonal
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
//Retorna o endereço de uma matriz
double* BandMatrix::getMatrix()
{
	return m_matrix;
}
//Resolve o sistema linear da forma Ax=b (sendo A uma matriz de banda)
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
	//Faz fatoracao LU
	MKL_Set_Num_Threads(1);
	//mkl_set_dynamic(1);
	dgbtrf(&n,&n,&kl,&ku,A.m_matrix,&ldab,ipiv,&info);
	if (info != 0)
	{
		printf("Error in LU factorization - bandsystem function\n");
		*info_fail = 1;
	}
	//Resolve o sistema linear
	dgbtrs(&trans,&n,&kl,&ku,&nrhs,A.m_matrix,&ldab,ipiv,b.getMatrix(),&n,&info);
	if (info != 0)
	{
		printf("Error in linear system - bandsystem function\n");
		*info_fail = 1;
	}
	
	delete []ipiv;
	return b;
}
//Operador Multiplicacao de matrizes (matrix2 deve necessariamente ser um vetor)
Matrix operator * (BandMatrix &matrix1, Matrix &matrix2)
{
	//Verificação da possibilidade de multiplicação
	if (matrix1.getColumns() != matrix2.getLines())
	{
		printf("Impossible to multiply matrices. Dimensions are uncompatible\n");
		return NULL;
	}
	else
	{
		Matrix return_m(matrix1.getColumns(),1);
		
		//Multiplicação
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
//Retorna a matriz quadrada (FULL)
Matrix BandMatrix::Full()
{

	Matrix return_m(m_columns,m_columns);
	
	for(long i=0; i < m_columns; i++)
	{
		for(long j=0; j < m_columns; j++)
		{
			//Se esta na região de valores não nulos
			if(i-j <= m_band_i && j-i <= m_band_s)
				return_m(i,j) = this->getValue(i,j);
			else
				return_m(i,j) = 0;
		}
	}
	
	return return_m;
}
//Operador Multiplicacao de matriz de banda por escalar
BandMatrix operator * (double sigma, BandMatrix &matrix1)
{
	BandMatrix return_m(matrix1.getBands(),matrix1.getBandi(),matrix1.getColumns());
	
	//Multiplicação
	for (int i=0;i<return_m.m_alloced_lines;i++)
		return_m.getMatrix()[i] = matrix1.getMatrix()[i]*sigma;

	return return_m;
}
//Operador soma de matrizes de banda
BandMatrix operator + (BandMatrix &matrix1, BandMatrix &matrix2)
{
	BandMatrix return_m(matrix1.getBands(),matrix1.getBandi(),matrix1.getColumns());

	//Verificação da possibilidade de multiplicação
	if (matrix1.getColumns() != matrix2.getColumns() || matrix1.getBandi() != matrix2.getBandi() || matrix1.getBands() != matrix2.getBands())
	{
		printf("Impossible to sum matrices. Dimensions are uncompatible\n");
	}
	else
	{
		//SOMA
		for (int i=0;i<return_m.m_alloced_lines;i++)
			return_m.getMatrix()[i] = matrix1.getMatrix()[i] + matrix2.getMatrix()[i];
	}
	return return_m;
}

long BandMatrix::getAllocedLines()
{
	return m_alloced_lines;
}
