#include "MatrixFloat.h"
#include "mkl.h"
#include <float.h>    
#define PI 3.1415926535897932384626433832795

//Construtor Padr�o
MatrixFloat::MatrixFloat(void)
{
	m_lines_deleted = true;
	m_alloced_lines = 0;
	m_lines = 1;
	m_columns = 1;
	//Inicializa a matriz como null e tenta aloc�-la
	m_matrix = NULL;
	if(!alloc())
		printf("Nao foi possivel alocar matriz! \n");
}
//Construtor de matriz quadrada
MatrixFloat::MatrixFloat(long lines)
{
	m_lines_deleted = true;
	m_alloced_lines = 0;
	m_lines = lines;
	m_columns = 1;
	//Inicializa a matriz como null e tenta aloc�-la
	m_matrix = NULL;
	if(!alloc())
		printf("Nao foi possivel alocar matriz! \n");
}
//Construtor Param�trico
MatrixFloat::MatrixFloat(long lines, long columns)
{
	m_lines_deleted = true;
	m_alloced_lines = 0;
	m_lines = lines;
	m_columns = columns;
	//Inicializa a matriz como null e tenta aloc�-la
	m_matrix = NULL;
	if(!alloc())
		printf("Nao foi possivel alocar matriz! \n");
}
//Construtor de c�pia
MatrixFloat::MatrixFloat(const MatrixFloat &copied)
{
	//Verifica dimens�es da matriz - se necess�rio, faz re-aloca��o
	//if (copied.m_alloced_lines != m_alloced_lines)
	//{
		m_alloced_lines = 0;
		m_lines = copied.m_lines;
		m_columns = copied.m_columns;
		m_lines_deleted = true;
		//Inicializa a matriz como null e tenta aloc�-la
		m_matrix = NULL;
		if(!alloc())
			printf("Nao foi possivel alocar matriz! \n");
	//}
	//Copia os valores na nova matriz
	for(long i=0; i < copied.m_alloced_lines; i++)
		m_matrix[i] = copied.m_matrix[i];	
}
//Destrutor Padr�o
MatrixFloat::~MatrixFloat(void)
{
	flush();
}

//Retorna o n�mero de linhas da matriz
long MatrixFloat::getLines()
{
	return this->m_lines;
}
//Retorna o n�mero de colunas da matriz 
long MatrixFloat::getColumns()
{
	return this->m_columns;
}
//Define o n�mero de linhas da matriz
void MatrixFloat::setLines(long value)
{
	m_lines = value;
}
//Define o n�mero de colunas da matriz
void MatrixFloat::setColumns(long value)
{
	m_columns = value;
}

//Imprime a matriz na tela
void MatrixFloat::print()
{
	printf("\n");
	for(long i=0; i < this->getLines(); i++)
	{
		printf("|");
		for(long j=0; j < this->getColumns(); j++)
		{
			if(m_matrix[i+j*m_lines] >= 0.0)
				printf(" %.4e ",m_matrix[i+j*m_lines]);
			else
				printf(" %.4e ",m_matrix[i+j*m_lines]);
		}
		printf("|\n");
	}
	printf("\n");
}
//Imprime a matriz em um arquivo de texto, cujo nome est� no char s
void MatrixFloat::fprint(char* s)
{
	FILE *file1 = fopen(s,"w");
	fprintf(file1,"\n");
	for(long i=0; i < this->getLines(); i++)
	{
		for(long j=0; j < this->getColumns(); j++)
		{
			if(m_matrix[i+j*m_lines] >= 0.0)
				fprintf(file1," %.14e\t",m_matrix[i+j*m_lines]);
			else
				fprintf(file1," %.14e\t",m_matrix[i+j*m_lines]);
		}
		fprintf(file1,"\n");
	}
	fprintf(file1,"\n");
	fclose(file1);
}
//Aloca a matriz
bool MatrixFloat::alloc()
{
	flush();
	//Tenta alocar a matriz
	if(!(m_matrix = new float[m_lines*m_columns]))
			//Falha
			return false;
	//Caso consiga
	else
	{
		m_alloced_lines = m_lines*m_columns;
		this->m_lines_deleted = false;
		
		for(long i=0; i < m_alloced_lines; i++)
				m_matrix[i] = 0.0;
		return true;
	}	
}
//Libera a mem�ria ocupada pela matriz
bool MatrixFloat::flush()
{
	//Percorre as linhas
	//Se h� linhas para desalocar
	if (this->m_lines_deleted == false)
	{
		if (m_matrix)
			delete[]m_matrix;
		//Desalocada com sucesso
		m_lines_deleted=true;	//para evitar desaloca��o novamente
		m_matrix = NULL;
		m_alloced_lines = 0;
	}
	return true;
}
//Zera a matriz, mantando as dimens�es atuais
void MatrixFloat::clear()
{
	for(int i=0; i < m_alloced_lines; i++)
			m_matrix[i] = 0.0;
}

//Operador Soma
MatrixFloat operator + (MatrixFloat &matrix1, MatrixFloat &matrix2)
{
	//Verifica se as dimens�es das matrizes s�o compat�veis
	if((matrix1.getLines() != matrix2.getLines()) || (matrix1.getColumns() != matrix2.getColumns()))
	{
		//Mensagem de erro - Dimens�es Incompat�veis
		printf("Matrizes devem possuir a mesma dimensao! \n");
		//Retorna vazio
		return NULL;
	}
	//Caso as dimens�es sejam compat�veis
	else
	{
		//Cria uma matriz de retorno
		MatrixFloat return_matrix(matrix1.getLines(),matrix1.getColumns());
		
		//Subtrai elemento a elemento da matriz
		for(long i=0; i < matrix1.m_alloced_lines; i++)
			return_matrix.m_matrix[i] = matrix1.m_matrix[i] + matrix2.m_matrix[i];
		
		//Retorna a matriz subtra��o
		return return_matrix;
	}
}
//Operador Soma
MatrixFloat operator + (Matrix &matrix1, MatrixFloat &matrix2)			
{
	//Verifica se as dimens�es das matrizes s�o compat�veis
	if ((matrix1.getLines() != matrix2.getLines()) || (matrix1.getColumns() != matrix2.getColumns()))
	{
		//Mensagem de erro - Dimens�es Incompat�veis
		printf("Matrizes devem possuir a mesma dimensao! \n");
		//Retorna vazio
		return NULL;
	}
	//Caso as dimens�es sejam compat�veis
	else
	{
		//Cria uma matriz de retorno
		MatrixFloat return_matrix(matrix1.getLines(), matrix1.getColumns());

		//Subtrai elemento a elemento da matriz
		for (long i = 0; i < matrix1.m_alloced_lines; i++)
			return_matrix.m_matrix[i] = (float)matrix1.m_matrix[i] + matrix2.m_matrix[i];

		//Retorna a matriz subtra��o
		return return_matrix;
	}
}

//Operador Subtra��o
MatrixFloat operator - (MatrixFloat &matrix1, MatrixFloat &matrix2)
{
	//Verifica se as dimens�es das matrizes s�o compat�veis
	if((matrix1.getLines() != matrix2.getLines()) || (matrix1.getColumns() != matrix2.getColumns()))
	{
		//Mensagem de erro - Dimens�es Incompat�veis
		printf("Matrizes devem possuir a mesma dimensao! \n");
		//Retorna vazio
		return NULL;
	}
	//Caso as dimens�es sejam compat�veis
	else
	{
		//Cria uma matriz de retorno
		MatrixFloat return_matrix(matrix1.getLines(),matrix1.getColumns());
		
		//Subtrai elemento a elemento da matriz
		for(long i=0; i < matrix1.m_alloced_lines; i++)
			return_matrix.m_matrix[i] = matrix1.m_matrix[i] - matrix2.m_matrix[i];
		
		//Retorna a matriz subtra��o
		return return_matrix;
	}
}
//Operador Multiplicacao de matrizes
MatrixFloat operator * (MatrixFloat &matrix1, MatrixFloat &matrix2)
{
	//Verifica��o da possibilidade de multiplica��o
	if (matrix1.m_columns != matrix2.m_lines)
	{
		//Produto escalar
		if (matrix1.m_lines == matrix2.m_lines)
		{
			MatrixFloat ret(1,1);
			for (int i = 0; i < matrix1.m_lines; i++)
				ret(0,0) += matrix1(i, 0)*matrix2(i, 0);
			return ret;
		}
		else
		{
			printf("Nao e possivel multiplicar as matrizes. Dimensoes incompativeis!");
			return NULL;
		}
	}
	else
	{
		
		MatrixFloat return_m(matrix1.m_lines,matrix2.m_columns);
		//Se a segunda matriz for um vetor
		if (matrix2.m_columns == 1)
			cblas_sgemv(CblasColMajor,CblasNoTrans,matrix1.m_lines,matrix1.m_columns,1.0,matrix1.m_matrix,matrix1.m_lines,matrix2.m_matrix,1,0.0,return_m.m_matrix,1);
		//Se n�o
		else
			cblas_sgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,matrix1.m_lines,matrix2.m_columns,matrix2.m_lines,1.0,matrix1.m_matrix,matrix1.m_lines,matrix2.m_matrix,matrix2.m_lines,1,return_m.m_matrix,return_m.m_lines);
		//return_m.print();
		return return_m;
	}
}

//Operador Multiplicacao de matrizes
MatrixFloat operator * (Matrix &matrix1, MatrixFloat &matrix2)
{
	//Verifica��o da possibilidade de multiplica��o
	if (matrix1.m_columns != matrix2.m_lines)
	{
		//Produto escalar
		if (matrix1.m_lines == matrix2.m_lines)
		{
			MatrixFloat ret(1, 1);
			for (int i = 0; i < matrix1.m_alloced_lines; i++)
				ret(0, 0) += ((float)matrix1.m_matrix[i]) * matrix2.m_matrix[i];
			return ret;
		}
		else
		{
			printf("Nao e possivel multiplicar as matrizes. Dimensoes incompativeis!");
			return NULL;
		}
	}
	else
	{
		float* tempmatrix = new float[matrix1.m_alloced_lines];
		for (int i = 0; i < matrix1.m_alloced_lines; i++)
			tempmatrix[i] = (float)matrix1.m_matrix[i];

		MatrixFloat return_m(matrix1.m_lines, matrix2.m_columns);
		//Se a segunda matriz for um vetor
		if (matrix2.m_columns == 1)
			cblas_sgemv(CblasColMajor, CblasNoTrans, matrix1.m_lines, matrix1.m_columns, 1.0, tempmatrix, matrix1.m_lines, matrix2.m_matrix, 1, 0.0, return_m.m_matrix, 1);
		//Se n�o
		else
			cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, matrix1.m_lines, matrix2.m_columns, matrix2.m_lines, 1.0, tempmatrix, matrix1.m_lines, matrix2.m_matrix, matrix2.m_lines, 1, return_m.m_matrix, return_m.m_lines);
		
		delete[]tempmatrix;
		//return_m.print();
		return return_m;
	}
}

//Operador Multiplicacao por escalar
MatrixFloat operator * (float escalar, MatrixFloat &matrix1)
{
	//Cria uma matriz de retorno
	MatrixFloat return_matrix(matrix1.getLines(),matrix1.getColumns());

	//Faz o produto de cada elemento da matriz pelo escalar
	for(long i=0; i < return_matrix.m_alloced_lines; i++)
		return_matrix.m_matrix[i] = matrix1.m_matrix[i]*escalar;
	
	//Retorna a matriz produto por escalar
	return return_matrix;

	
}
//Operador Multiplicacao por escalar
MatrixFloat operator * (MatrixFloat &matrix1, float escalar)
{
	//Retorna a fun��o anterior, j� que esta opera��o � comutativa
	return escalar*matrix1;
}
//Verifica��o de igualdade
bool operator == (MatrixFloat &matrix1, MatrixFloat &matrix2)
{
	//Verifica a igualdade de dimens�es
	if((matrix1.getLines() != matrix2.getLines()) || (matrix1.getColumns() != matrix2.getColumns()))
		return false;
	else
	{
		//Varre a matriz e verifica a igualdade de elementos
		for(long i=0; i < matrix1.m_alloced_lines; i++)
				if(matrix1.m_matrix[i] != matrix2.m_matrix[i]) return false;
		//Matrizes iguais
		return true;
	}
}
//Verifica��o de inegualdade
bool operator != (MatrixFloat &matrix1, MatrixFloat &matrix2)
{
	return !(matrix1 == matrix2);
}
//Operador de Atribui��o	
MatrixFloat &MatrixFloat::operator = (MatrixFloat const &matrix1)
{
	//Verifica dimens�es da matriz - se necess�rio, faz re-aloca��o
	if (matrix1.m_alloced_lines != m_alloced_lines)
	{
		this->flush();
		m_alloced_lines = 0;
		m_lines = matrix1.m_lines;
		m_columns = matrix1.m_columns;
		m_lines_deleted = matrix1.m_lines_deleted;
		//Inicializa a matriz como null e tenta aloc�-la
		m_matrix = NULL;
		if(!alloc())
			printf("Nao foi possivel alocar matriz! \n");
	}
	//Copia os valores na nova matriz
	for(long i=0; i < matrix1.m_alloced_lines; i++)
		m_matrix[i] = matrix1.m_matrix[i];

	//Retorna esta matriz
	return *this;
}

//Operador de Atribui��o 2
MatrixFloat &MatrixFloat::operator = (Matrix const &matrix1)
{
	//Verifica dimens�es da matriz - se necess�rio, faz re-aloca��o
	if (matrix1.m_alloced_lines != m_alloced_lines)
	{
		this->flush();
		m_alloced_lines = 0;
		m_lines = matrix1.m_lines;
		m_columns = matrix1.m_columns;
		m_lines_deleted = matrix1.m_lines_deleted;
		//Inicializa a matriz como null e tenta aloc�-la
		m_matrix = NULL;
		if (!alloc())
			printf("Nao foi possivel alocar matriz! \n");
	}
	//Copia os valores na nova matriz
	for (long i = 0; i < matrix1.m_alloced_lines; i++)
		m_matrix[i] = (float)matrix1.m_matrix[i];

	//Retorna esta matriz
	return *this;
}

//Retorno do valor na posi��o especificada
float &MatrixFloat::operator() (long line, long column)
{
	//Verifica se a posi��o � v�lida
	if(line > this->m_lines-1 || column > this->m_columns-1 || line < 0 || column < 0)
	{
		printf("Not valid position accessed in matrix! (%d,%d)\n",line,column);
		float* ret = new float[1];
		ret[0] = 0;
		return ret[0];
	}
	else
		//Retorna o valor na posi��o desejada
		return this->m_matrix[line + column*m_lines];
}
//Operador produto escalar entre dois vetores
float dot(MatrixFloat &matrix1, MatrixFloat &matrix2)
{
	//Verifica��o da possibilidade do produto
	if (matrix1.m_lines != matrix2.m_lines)
	{
		printf("Nao e possivel calcular o produto escalar. Dimensoes incompativeis!");
		return 0;
	}
	else
	{
		float return_value = 0.0;
		for (int i = 0; i < matrix1.m_alloced_lines; i++)
			return_value+=matrix1.m_matrix[i]*matrix2.m_matrix[i];
		return return_value;
	}
}

//Operador produto escalar entre dois vetores
double dot(MatrixFloat &matrix1, Matrix &matrix2)
{
	//Verifica��o da possibilidade do produto
	if (matrix1.m_lines != matrix2.m_lines)
	{
		printf("Nao e possivel calcular o produto escalar. Dimensoes incompativeis!");
		return 0;
	}
	else
	{
		double return_value = 0.0;
		for (int i = 0; i < matrix1.m_alloced_lines; i++)
			return_value += matrix1.m_matrix[i] * matrix2.m_matrix[i];
		return return_value;
	}
}

//Operador produto vetorial entre dois vetores
MatrixFloat cross(MatrixFloat &matrix1, MatrixFloat &matrix2)
{
	//Verifica��o da possibilidade do produto
	if (matrix1.m_columns != 1 || matrix2.m_columns != 1 || matrix1.m_lines != 3 || matrix2.m_lines != 3)
	{
		printf("Nao e possivel calcular o produto vetorial. Dimensoes incompativeis!");
		return 0;
	}
	else
	{
		MatrixFloat return_m(3);
		return_m(0,0) = matrix1(1,0)*matrix2(2,0) - matrix1(2,0)*matrix2(1,0);
		return_m(1,0) = matrix1(2,0)*matrix2(0,0) - matrix1(0,0)*matrix2(2,0);
		return_m(2,0) = matrix1(0,0)*matrix2(1,0) - matrix1(1,0)*matrix2(0,0);
		return return_m;
	}
}


//Operador produto tensorial entre dois vetores
MatrixFloat dyadic(MatrixFloat &matrix1, MatrixFloat &matrix2)
{
	//Verifica��o da possibilidade do produto
	if (matrix1.m_columns != 1 || matrix2.m_columns != 1 || matrix1.m_lines != matrix2.m_lines)
	{
		printf("Nao e possivel calcular o produto tensorial. Dimensoes incompativeis!");
		return 0;
	}
	else
	{
		int order = matrix1.m_lines;
		MatrixFloat return_m(order, order);
		for (int i = 0; i < order; i++)
		for (int j = 0; j < order; j++)
				return_m(i,j)=matrix1(i,0)*matrix2(j,0);
		return return_m;
	}
}
//Operador skew de um vetor
MatrixFloat skew(MatrixFloat &matrix1)
{
	//Verifica��o da possibilidade do produto
	if (matrix1.m_columns != 1 || matrix1.m_lines != 3)
	{
		printf("Nao e possivel calcular o produto escalar. Dimensoes incompativeis!");
		return 0;
	}
	else
	{
		MatrixFloat return_m(3,3);
		return_m(0,1)=-matrix1(2,0);
		return_m(0,2)=+matrix1(1,0);
		return_m(1,2)=-matrix1(0,0);

		return_m(1,0)=+matrix1(2,0);
		return_m(2,0)=-matrix1(1,0);
		return_m(2,1)=+matrix1(0,0);
		return return_m;
	}
}
//Operador axial de um vetor
MatrixFloat axial(MatrixFloat &matrix1)
{
	//Verifica��o da possibilidade do produto
	if (matrix1.m_columns != 3 || matrix1.m_lines != 3)
	{
		printf("Nao e possivel calcular o produto escalar. Dimensoes incompativeis!");
		return 0;
	}
	else
	{
		MatrixFloat return_m(3);
		return_m(0,0)=-matrix1(1,2);
		return_m(1,0)=+matrix1(0,2);
		return_m(2,0)=-matrix1(0,1);
		return return_m;
	}
}

//Retorna o endere�o de uma matriz
float* MatrixFloat::getMatrix()
{
	return m_matrix;
}

//Retorna a norma de um vetor
float norm(MatrixFloat &matrix1)
{
	if (matrix1.getColumns() != 1)
		printf("Dimensao nao consistente para calculo da norma");
	else
	{
		if (matrix1.getLines() != 3)
		{
			if (matrix1.getLines() == 2)
			{
				float return_value = sqrt(matrix1(0, 0)*matrix1(0, 0) +
					matrix1(1, 0)*matrix1(1, 0));
				return return_value;
			}
			if (matrix1.getLines() == 4)
			{
				float return_value = sqrt(matrix1(0, 0)*matrix1(0, 0) +
					matrix1(1, 0)*matrix1(1, 0) +
					matrix1(2, 0)*matrix1(2, 0) + 
					matrix1(3, 0)*matrix1(3, 0) );
				return return_value;
			}
			if (matrix1.getLines() == 6)
			{
				float return_value = sqrt(matrix1(0, 0)*matrix1(0, 0) +
					matrix1(1, 0)*matrix1(1, 0) +
					matrix1(2, 0)*matrix1(2, 0) +
					matrix1(3, 0)*matrix1(3, 0) + 
					matrix1(4, 0)*matrix1(4, 0) +
					matrix1(5, 0)*matrix1(5, 0));
				return return_value;
			}
			//printf("Norma infinito\n");
			float max = 0;
			for (int i=0; i< matrix1.getLines(); i++)
			{
				//Detec��o de NaN
				if (matrix1(i, 0) == matrix1(i, 0))
				{
					//Detec��o de infinito
					if (matrix1(i, 0) >= FLT_MAX || matrix1(i, 0) <= -FLT_MAX)
						return FLT_MAX;
					else
					{
						if (max < abs(matrix1(i, 0)))
						{
							max = abs(matrix1(i, 0));
							//printf("GL %d\n",i + 1);
						}
							
					}
					
				}
				else
					return FLT_MAX;//valor muito alto, pois detectou NaN
				
			}
			return max;
		}
		else
		{
			float return_value = sqrt( matrix1(0,0)*matrix1(0,0) +
										matrix1(1,0)*matrix1(1,0) +
										matrix1(2,0)*matrix1(2,0) );
			return return_value;
		}
	}
	return 0;
}

//Retorna a transposta de uma matriz
MatrixFloat transp(MatrixFloat &matrix1)
{
	MatrixFloat answer(matrix1.getColumns(),matrix1.getLines());
	for (int j=0;j<matrix1.getColumns();j++)
	{
		for (int i=0;i<matrix1.getLines();i++)
		{	
			answer(j,i)=matrix1(i,j);
		}
	}
	return answer;
}
//Zera a matriz
void zeros(MatrixFloat* matrix1)
{
	for (int j=0;j<matrix1->getColumns();j++)
	{
		for (int i=0;i<matrix1->getLines();i++)
		{	
			(*matrix1)(i,j) = 0.0;
		}
	}
}
