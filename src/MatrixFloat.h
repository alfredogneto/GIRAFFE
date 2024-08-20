#pragma once
#include "mkl.h"
#include "mkl_blas.h"
#include "mkl_types.h"
#include "mkl_lapack.h"
#include "mkl_dss.h"
#include "stdio.h"
#include <cmath>
#include "arpack.h"
#include "Matrix.h"

class MatrixFloat
{
public:
	//Constructors and Destrutor
	MatrixFloat(void);													//Construtor Padr�o
	MatrixFloat(long lines);											//Construtor de matriz coluna
	MatrixFloat(long lines, long columns);								//Construtor Param�trico
	MatrixFloat(const MatrixFloat &copied);									//Construtor de c�pia
	~MatrixFloat(void);													//Destrutor Padr�o

	//Gets and Sets
	long getLines();												//Retorna o n�mero de linhas da matriz
	long getColumns();												//Retorna o n�mero de colunas da matriz 

	void setLines(long value);										//Define o n�mero de linhas da matriz 
	void setColumns(long value);									//Define o n�mero de colunas da matriz 
	float* getMatrix();												//Retorna o endere�o de uma matriz
	
	//General Functions
	void print();													//Imprime a matriz no console
	void fprint(char* s);											//Imprime a matriz em um arquivo de texto, cujo nome est� no char s
	bool alloc();													//Aloca a matriz
	bool flush();													//Libera a mem�ria ocupada pela matriz
	void clear();													//Zera a matriz, mantando as dimens�es atuais

	
	float &operator() (long line, long column);						//Retorno do valor na posi��o especificada
	MatrixFloat &operator = (MatrixFloat const &matrix1);			//Operador de Atribui��o
	MatrixFloat &operator = (Matrix const &matrix1);				//Operador de Atribui��o 2

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	float*  m_matrix;												//Matriz unidimensional
	long     m_lines;												//N�mero de linhas
	long     m_columns;												//N�mero de colunas
	long	 m_alloced_lines;										//N�mero de linhas atualmente alocadas
	bool	 m_lines_deleted;										//Flag que indica se as linhas foram desalocadas				
};

MatrixFloat operator + (MatrixFloat &matrix1, MatrixFloat &matrix2);				//Operador Soma
MatrixFloat operator + (Matrix &matrix1, MatrixFloat &matrix2);						//Operador Soma
MatrixFloat operator - (MatrixFloat &matrix1, MatrixFloat &matrix2);				//Operador Subtra��o
MatrixFloat operator * (MatrixFloat &matrix1, MatrixFloat &matrix2);				//Operador Multiplicacao de matrizes
MatrixFloat operator * (Matrix &matrix1, MatrixFloat &matrix2);				//Operador Multiplicacao de matrizes
MatrixFloat operator * (float escalar, MatrixFloat &matrix1);				//Operador Multiplicacao por escalar
MatrixFloat operator * (MatrixFloat &matrix1, float escalar);				//Operador Multiplicacao por escalar
bool operator == (MatrixFloat &matrix1, MatrixFloat &matrix2);				//Verifica��o de igualdade
bool operator != (MatrixFloat &matrix1, MatrixFloat &matrix2);				//Verifica��o de inegualdade
float dot(MatrixFloat &matrix1, MatrixFloat &matrix2);						//Produto escalar entre dois vetores
double dot(MatrixFloat &matrix1, Matrix &matrix2);							//Produto escalar entre dois vetores
MatrixFloat cross(MatrixFloat &matrix1, MatrixFloat &matrix2);						//Operador produto vetorial entre duas matrizes
MatrixFloat dyadic(MatrixFloat &matrix1, MatrixFloat &matrix2);					//Operador produto vetorial entre duas matrizes
MatrixFloat skew(MatrixFloat &matrix1);										//Operador skew de um vetor
MatrixFloat axial(MatrixFloat &matrix1);										//Operador axial de um vetor

float norm(MatrixFloat &matrix1);										//Retorna a norma de um vetor
MatrixFloat transp(MatrixFloat &matrix1);								//Retorna a transposta de uma matriz
void zeros(MatrixFloat* matrix1);										//Zera a matriz
