#pragma once
#include "Matrix.h"

class MatrixFloat
{
public:
	//Constructors and Destructor
	MatrixFloat(void);													//Constructor Padr�o
	MatrixFloat(long lines);											//Constructor of column matrix
	MatrixFloat(long lines, long columns);								//Parametric constructor
	MatrixFloat(const MatrixFloat &copied);									//Constructor of c�pia
	~MatrixFloat(void);													//Destructor Padr�o

	//Gets and Sets
	long getLines() const;											//Returns the number of rows of the matrix
	long getColumns() const;										//Returns the number of columns of the matrix

	void setLines(long value);										//Defines the number of rows of the matrix
	void setColumns(long value);									//Defines the number of columns of the matrix
	float* getMatrix() const;										//Returns the endere�o of a matrix

	//General Functions
	void print();													//Prints the matrix in the console
	void fprint(char* s);											//Prints the matrix in a file of text, whose name this in the char s
	bool alloc();													//Allocates the matrix
	bool flush();													//Frees the mem�ria occupied through the matrix
	void clear();													//Zeroes the matrix, keeping the dimens�es current


	float &operator() (long line, long column);						//Return of the value in the posi��o specified
	float operator() (long line, long column) const;					//Return of the value in the posi��o specified (version const)
	MatrixFloat &operator = (MatrixFloat const &matrix1);			//Operator of Atribui��o
	MatrixFloat &operator = (Matrix const &matrix1);				//Operator of Atribui��o 2

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	float*  m_matrix;												//Matrix one-dimensional
	long     m_lines;												//Number of rows
	long     m_columns;												//Number of columns
	long	 m_alloced_lines;										//Number of rows currently allocated
	bool	 m_lines_deleted;										//Flag that indicates if the rows were deallocated				
};

MatrixFloat operator + (const MatrixFloat &matrix1, const MatrixFloat &matrix2);				//Operator Sum
MatrixFloat operator + (const Matrix &matrix1, const MatrixFloat &matrix2);						//Operator Sum
MatrixFloat operator - (const MatrixFloat &matrix1, const MatrixFloat &matrix2);				//Operator Subtra��o
MatrixFloat operator * (const MatrixFloat &matrix1, const MatrixFloat &matrix2);				//Operator Multiplication of matrices
MatrixFloat operator * (const Matrix &matrix1, const MatrixFloat &matrix2);				//Operator Multiplication of matrices
MatrixFloat operator * (float escalar, const MatrixFloat &matrix1);				//Operator Multiplication by scalar
MatrixFloat operator * (const MatrixFloat &matrix1, float escalar);				//Operator Multiplication by scalar
bool operator == (const MatrixFloat &matrix1, const MatrixFloat &matrix2);				//Verifica��o of equality
bool operator != (const MatrixFloat &matrix1, const MatrixFloat &matrix2);				//Verifica��o of inequality
float dot(const MatrixFloat &matrix1, const MatrixFloat &matrix2);						//Dot product between two vectors
double dot(const MatrixFloat &matrix1, const Matrix &matrix2);							//Dot product between two vectors
MatrixFloat cross(const MatrixFloat &matrix1, const MatrixFloat &matrix2);						//Operator cross product between two matrices
MatrixFloat dyadic(const MatrixFloat &matrix1, const MatrixFloat &matrix2);					//Operator cross product between two matrices
MatrixFloat skew(const MatrixFloat &matrix1);										//Operator skew of a vector
MatrixFloat axial(const MatrixFloat &matrix1);										//Operator axial of a vector

float norm(const MatrixFloat &matrix1);										//Returns the norm of a vector
MatrixFloat transp(const MatrixFloat &matrix1);								//Returns the transposta of a matrix
void zeros(MatrixFloat* matrix1);										//Zeroes the matrix
