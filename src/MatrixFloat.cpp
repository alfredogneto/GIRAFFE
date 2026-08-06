#include "MatrixFloat.h"
#ifdef GIRAFFE_USE_MKL
#include <mkl.h>
#else
#include <cblas.h>
#endif
#include <float.h>
#define PI 3.1415926535897932384626433832795

//Constructor Padr�o
MatrixFloat::MatrixFloat(void)
{
	m_lines_deleted = true;
	m_alloced_lines = 0;
	m_lines = 1;
	m_columns = 1;
	//Initializes the matrix as null and tries allocates-it
	m_matrix = NULL;
	if(!alloc())
		printf("Nao foi possivel alocar matriz! \n");
}
//Constructor of square matrix
MatrixFloat::MatrixFloat(long lines)
{
	m_lines_deleted = true;
	m_alloced_lines = 0;
	m_lines = lines;
	m_columns = 1;
	//Initializes the matrix as null and tries allocates-it
	m_matrix = NULL;
	if(!alloc())
		printf("Nao foi possivel alocar matriz! \n");
}
//Parametric constructor
MatrixFloat::MatrixFloat(long lines, long columns)
{
	m_lines_deleted = true;
	m_alloced_lines = 0;
	m_lines = lines;
	m_columns = columns;
	//Initializes the matrix as null and tries allocates-it
	m_matrix = NULL;
	if(!alloc())
		printf("Nao foi possivel alocar matriz! \n");
}
//Constructor of c�pia
MatrixFloat::MatrixFloat(const MatrixFloat &copied)
{
	//Checks dimens�es of the matrix - if necessary, does re-aloca��o
	//if (copied.m_alloced_lines != m_alloced_lines)
	//{
		m_alloced_lines = 0;
		m_lines = copied.m_lines;
		m_columns = copied.m_columns;
		m_lines_deleted = true;
		//Initializes the matrix as null and tries allocates-it
		m_matrix = NULL;
		if(!alloc())
			printf("Nao foi possivel alocar matriz! \n");
	//}
	//Copies the values in the new matrix
	for(long i=0; i < copied.m_alloced_lines; i++)
		m_matrix[i] = copied.m_matrix[i];	
}
//Destructor Padr�o
MatrixFloat::~MatrixFloat(void)
{
	flush();
}

//Returns the number of rows of the matrix
long MatrixFloat::getLines() const
{
	return this->m_lines;
}
//Returns the number of columns of the matrix 
long MatrixFloat::getColumns() const
{
	return this->m_columns;
}
//Defines the number of rows of the matrix
void MatrixFloat::setLines(long value)
{
	m_lines = value;
}
//Defines the number of columns of the matrix
void MatrixFloat::setColumns(long value)
{
	m_columns = value;
}

//Prints the matrix in the tela
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
//Prints the matrix in a file of text, whose name this in the char s
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
//Allocates the matrix
bool MatrixFloat::alloc()
{
	flush();
	//Tries allocate the matrix
	if(!(m_matrix = new float[m_lines*m_columns]))
			//Falha
			return false;
	//Case consiga
	else
	{
		m_alloced_lines = m_lines*m_columns;
		this->m_lines_deleted = false;
		
		for(long i=0; i < m_alloced_lines; i++)
				m_matrix[i] = 0.0;
		return true;
	}	
}
//Frees the mem�ria occupied through the matrix
bool MatrixFloat::flush()
{
	//Iterates over the rows
	//If there is rows for deallocate
	if (this->m_lines_deleted == false)
	{
		if (m_matrix)
			delete[]m_matrix;
		//Desalocada with sucesso
		m_lines_deleted=true;	//for avoid desaloca��o novamente
		m_matrix = NULL;
		m_alloced_lines = 0;
	}
	return true;
}
//Zeroes the matrix, keeping the dimens�es current
void MatrixFloat::clear()
{
	for(int i=0; i < m_alloced_lines; i++)
			m_matrix[i] = 0.0;
}

//Operator Sum
MatrixFloat operator + (const MatrixFloat &matrix1, const MatrixFloat &matrix2)
{
	//Checks if the dimens�es of the matrices s�o compatible
	if((matrix1.getLines() != matrix2.getLines()) || (matrix1.getColumns() != matrix2.getColumns()))
	{
		//Mensagem of error - Dimens�es Incompatible
		printf("Matrizes devem possuir a mesma dimensao! \n");
		//Returns vazio
		return NULL;
	}
	//Case the dimens�es be compatible
	else
	{
		//Creates a matrix of return
		MatrixFloat return_matrix(matrix1.getLines(),matrix1.getColumns());
		
		//Subtrai element the element of the matrix
		for(long i=0; i < matrix1.m_alloced_lines; i++)
			return_matrix.m_matrix[i] = matrix1.m_matrix[i] + matrix2.m_matrix[i];
		
		//Returns the matrix subtra��o
		return return_matrix;
	}
}
//Operator Sum
MatrixFloat operator + (const Matrix &matrix1, const MatrixFloat &matrix2)			
{
	//Checks if the dimens�es of the matrices s�o compatible
	if ((matrix1.getLines() != matrix2.getLines()) || (matrix1.getColumns() != matrix2.getColumns()))
	{
		//Mensagem of error - Dimens�es Incompatible
		printf("Matrizes devem possuir a mesma dimensao! \n");
		//Returns vazio
		return NULL;
	}
	//Case the dimens�es be compatible
	else
	{
		//Creates a matrix of return
		MatrixFloat return_matrix(matrix1.getLines(), matrix1.getColumns());

		//Subtrai element the element of the matrix
		for (long i = 0; i < matrix1.m_alloced_lines; i++)
			return_matrix.m_matrix[i] = (float)matrix1.m_matrix[i] + matrix2.m_matrix[i];

		//Returns the matrix subtra��o
		return return_matrix;
	}
}

//Operator Subtra��o
MatrixFloat operator - (const MatrixFloat &matrix1, const MatrixFloat &matrix2)
{
	//Checks if the dimens�es of the matrices s�o compatible
	if((matrix1.getLines() != matrix2.getLines()) || (matrix1.getColumns() != matrix2.getColumns()))
	{
		//Mensagem of error - Dimens�es Incompatible
		printf("Matrizes devem possuir a mesma dimensao! \n");
		//Returns vazio
		return NULL;
	}
	//Case the dimens�es be compatible
	else
	{
		//Creates a matrix of return
		MatrixFloat return_matrix(matrix1.getLines(),matrix1.getColumns());
		
		//Subtrai element the element of the matrix
		for(long i=0; i < matrix1.m_alloced_lines; i++)
			return_matrix.m_matrix[i] = matrix1.m_matrix[i] - matrix2.m_matrix[i];
		
		//Returns the matrix subtra��o
		return return_matrix;
	}
}
//Operator Multiplication of matrices
MatrixFloat operator * (const MatrixFloat &matrix1, const MatrixFloat &matrix2)
{
	//Verifica��o of the possibility of multiplica��o
	if (matrix1.m_columns != matrix2.m_lines)
	{
		//Dot product
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
		//If the second matrix for a vector
		if (matrix2.m_columns == 1)
			cblas_sgemv(CblasColMajor,CblasNoTrans,matrix1.m_lines,matrix1.m_columns,1.0,matrix1.m_matrix,matrix1.m_lines,matrix2.m_matrix,1,0.0,return_m.m_matrix,1);
		//Se n�o
		else
			cblas_sgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,matrix1.m_lines,matrix2.m_columns,matrix2.m_lines,1.0,matrix1.m_matrix,matrix1.m_lines,matrix2.m_matrix,matrix2.m_lines,1,return_m.m_matrix,return_m.m_lines);
		//return_m.print();
		return return_m;
	}
}

//Operator Multiplication of matrices
MatrixFloat operator * (const Matrix &matrix1, const MatrixFloat &matrix2)
{
	//Verifica��o of the possibility of multiplica��o
	if (matrix1.m_columns != matrix2.m_lines)
	{
		//Dot product
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
		//If the second matrix for a vector
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

//Operator Multiplication by scalar
MatrixFloat operator * (float escalar, const MatrixFloat &matrix1)
{
	//Creates a matrix of return
	MatrixFloat return_matrix(matrix1.getLines(),matrix1.getColumns());

	//Does the product of each element of the matrix through the scalar
	for(long i=0; i < return_matrix.m_alloced_lines; i++)
		return_matrix.m_matrix[i] = matrix1.m_matrix[i]*escalar;
	
	//Returns the matrix product by scalar
	return return_matrix;

	
}
//Operator Multiplication by scalar
MatrixFloat operator * (const MatrixFloat &matrix1, float escalar)
{
	//Returns the fun��o previous, already that this opera��o and comutativa
	return escalar*matrix1;
}
//Verifica��o of equality
bool operator == (const MatrixFloat &matrix1, const MatrixFloat &matrix2)
{
	//Checks the equality of dimens�es
	if((matrix1.getLines() != matrix2.getLines()) || (matrix1.getColumns() != matrix2.getColumns()))
		return false;
	else
	{
		//Varre the matrix and checks the equality of elements
		for(long i=0; i < matrix1.m_alloced_lines; i++)
				if(matrix1.m_matrix[i] != matrix2.m_matrix[i]) return false;
		//Matrices iguais
		return true;
	}
}
//Verifica��o of inequality
bool operator != (const MatrixFloat &matrix1, const MatrixFloat &matrix2)
{
	return !(matrix1 == matrix2);
}
//Operator of Atribui��o	
MatrixFloat &MatrixFloat::operator = (MatrixFloat const &matrix1)
{
	//Checks dimens�es of the matrix - if necessary, does re-aloca��o
	if (matrix1.m_alloced_lines != m_alloced_lines)
	{
		this->flush();
		m_alloced_lines = 0;
		m_lines = matrix1.m_lines;
		m_columns = matrix1.m_columns;
		m_lines_deleted = matrix1.m_lines_deleted;
		//Initializes the matrix as null and tries allocates-it
		m_matrix = NULL;
		if(!alloc())
			printf("Nao foi possivel alocar matriz! \n");
	}
	//Copies the values in the new matrix
	for(long i=0; i < matrix1.m_alloced_lines; i++)
		m_matrix[i] = matrix1.m_matrix[i];

	//Returns this matrix
	return *this;
}

//Operator of Atribui��o 2
MatrixFloat &MatrixFloat::operator = (Matrix const &matrix1)
{
	//Checks dimens�es of the matrix - if necessary, does re-aloca��o
	if (matrix1.m_alloced_lines != m_alloced_lines)
	{
		this->flush();
		m_alloced_lines = 0;
		m_lines = matrix1.m_lines;
		m_columns = matrix1.m_columns;
		m_lines_deleted = matrix1.m_lines_deleted;
		//Initializes the matrix as null and tries allocates-it
		m_matrix = NULL;
		if (!alloc())
			printf("Nao foi possivel alocar matriz! \n");
	}
	//Copies the values in the new matrix
	for (long i = 0; i < matrix1.m_alloced_lines; i++)
		m_matrix[i] = (float)matrix1.m_matrix[i];

	//Returns this matrix
	return *this;
}

//Return of the value in the posi��o specified
float &MatrixFloat::operator() (long line, long column)
{
	//Checks if the posi��o and valida
	if(line > this->m_lines-1 || column > this->m_columns-1 || line < 0 || column < 0)
	{
		printf("Not valid position accessed in matrix! (%d,%d)\n",line,column);
		float* ret = new float[1];
		ret[0] = 0;
		return ret[0];
	}
	else
		//Returns the value in the posi��o desired
		return this->m_matrix[line + column*m_lines];
}
//Return of the value in the posi��o specified (version const)
float MatrixFloat::operator() (long line, long column) const
{
	if(line > this->m_lines-1 || column > this->m_columns-1 || line < 0 || column < 0)
	{
		printf("Not valid position accessed in matrix! (%d,%d)\n",line,column);
		return 0;
	}
	else
		return this->m_matrix[line + column*m_lines];
}
//Operator dot product between two vectors
float dot(const MatrixFloat &matrix1, const MatrixFloat &matrix2)
{
	//Verifica��o of the possibility of the product
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

//Operator dot product between two vectors
double dot(const MatrixFloat &matrix1, const Matrix &matrix2)
{
	//Verifica��o of the possibility of the product
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

//Operator cross product between two vectors
MatrixFloat cross(const MatrixFloat &matrix1, const MatrixFloat &matrix2)
{
	//Verifica��o of the possibility of the product
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


//Operator tensor (dyadic) product between two vectors
MatrixFloat dyadic(const MatrixFloat &matrix1, const MatrixFloat &matrix2)
{
	//Verifica��o of the possibility of the product
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
//Operator skew of a vector
MatrixFloat skew(const MatrixFloat &matrix1)
{
	//Verifica��o of the possibility of the product
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
//Operator axial of a vector
MatrixFloat axial(const MatrixFloat &matrix1)
{
	//Verifica��o of the possibility of the product
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

//Returns the endere�o of a matrix
float* MatrixFloat::getMatrix() const
{
	return m_matrix;
}

//Returns the norm of a vector
float norm(const MatrixFloat &matrix1)
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
			//printf("Norm infinito\n");
			float max = 0;
			for (int i=0; i< matrix1.getLines(); i++)
			{
				//Detec��o of NaN
				if (matrix1(i, 0) == matrix1(i, 0))
				{
					//Detec��o of infinito
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
					return FLT_MAX;//value very alto, since detectou NaN
				
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

//Returns the transposta of a matrix
MatrixFloat transp(const MatrixFloat &matrix1)
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
//Zeroes the matrix
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
