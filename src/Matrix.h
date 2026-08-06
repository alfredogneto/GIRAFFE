#pragma once
#include <stdio.h>
#include <cmath>

class Matrix
{
public:
	//Constructors and Destructor
	Matrix(void);													//Constructor Padr�o
	Matrix(long lines);												//Constructor of column matrix
	Matrix(long lines, long columns);								//Parametric constructor
	Matrix(const Matrix &copied);									//Constructor of c�pia
	~Matrix(void);													//Destructor Padr�o

	//Gets and Sets
	long getLines() const;											//Returns the number of rows of the matrix
	long getColumns() const;										//Returns the number of columns of the matrix

	void setLines(long value);										//Defines the number of rows of the matrix
	void setColumns(long value);									//Defines the number of columns of the matrix
	double* getMatrix() const;										//Returns the endere�o of a matrix

	//General Functions
	void print();													//Prints the matrix in the console
	void fprint(char* s);											//Prints the matrix in a file of text, whose name this in the char s
	bool alloc();													//Allocates the matrix
	bool flush();													//Frees the mem�ria occupied through the matrix
	void clear();													//Zeroes the matrix, keeping the dimens�es current

	//Operators
	void MatrixToPtr(double** ptr, int order) const;				//Saves pointer double** in ptr
	void PtrToMatrix(double** ptr, int order);						//Saves in the matrix the content of the double**
	void PtrToMatrix(double** ptr, int lines, int columns);			//Saves in the matrix the content of the double**
	double &operator() (long line, long column);					//Return of the value in the posi��o specified
	double operator() (long line, long column) const;				//Return of the value in the posi��o specified (version const)
	Matrix &operator = (Matrix const &matrix1);						//Operator of Atribui��o

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	double*  m_matrix;												//Matrix one-dimensional
	long     m_lines;												//Number of rows
	long     m_columns;												//Number of columns
	long	 m_alloced_lines;										//Number of rows currently allocated
	bool	 m_lines_deleted;										//Flag that indicates if the rows were deallocated				
};

Matrix operator + (const Matrix &matrix1, const Matrix &matrix2);				//Operator Sum
Matrix operator - (const Matrix &matrix1, const Matrix &matrix2);				//Operator Subtra��o
Matrix operator * (const Matrix &matrix1, const Matrix &matrix2);				//Operator Multiplication of matrices
Matrix operator * (double escalar, const Matrix &matrix1);				//Operator Multiplication by scalar
Matrix operator * (const Matrix &matrix1, double escalar);				//Operator Multiplication by scalar
bool operator == (const Matrix &matrix1, const Matrix &matrix2);				//Verifica��o of equality
bool operator != (const Matrix &matrix1, const Matrix &matrix2);				//Verifica��o of inequality
double dot(const Matrix &matrix1, const Matrix &matrix2);						//Dot product between two vectors
Matrix cross(const Matrix &matrix1, const Matrix &matrix2);						//Operator cross product between two matrices
Matrix dyadic(const Matrix &matrix1, const Matrix &matrix2);					//Operator cross product between two matrices
Matrix skew(const Matrix &matrix1);										//Operator skew of a vector
Matrix axial(const Matrix &matrix1);										//Operator axial of a vector
Matrix fullsystem(Matrix &A, Matrix b, int *flag_error);			//Solves the linear system of the way Ax=b (mutates A in-place via LAPACK; b is taken by value since it's used the scratch space and returned)
double norm(const Matrix &matrix1);										//Returns the norm of a vector
double norm4(const Matrix &matrix1);										//Returns the norm of a vector considerando only the 4 primeiros degrees of freedom
Matrix transp(const Matrix &matrix1);										//Returns the transposta of a matrix
void zeros(Matrix* matrix1);										//Zeroes the matrix
Matrix invert2x2(const Matrix &matrix);									//Inverte a matrix 2x2
Matrix invert3x3(const Matrix &matrix);									//Inverte a matrix 3x3
Matrix invert4x4(const Matrix &matrix);									//Inverte a matrix 4x4
Matrix invert5x5(const Matrix &matrix);									//Inverte a matrix 5x5
Matrix invert6x6(const Matrix &matrix);									//Inverte a matrix 6x6
Matrix invert(const Matrix &matrix);										//Inverte a matrix of 2x2 the 6x6 escolhendo automaticamente the fun��o correct

int fulleigen1(Matrix &A, Matrix &P, Matrix &D, double abstol);		//Calculates the eigenvalues and eigenvectors of a matrix simetrica
int fulleigen2(Matrix &A, Matrix &P, Matrix &D);					//Calculates the eigenvalues and eigenvectors of a matrix simetrica
double mineigen(Matrix &A, Matrix &P, Matrix &D, double abstol);	//Calculates the smallest eigenvalue of a matrix simetrica
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Matrix V(Matrix x, Matrix t, double alpha_escalar);					//Fun��o for the calculation of the operator V (for assembly of the geometric stiffness matrix)
Matrix d_V(Matrix x, Matrix d_x, Matrix t, double alpha_escalar);	//Fun��o for the calculation of the operator d_V (for assembly of the geometric stiffness matrix)

double ArcReduction(double arc);									//Calculates redu��o the first turn [-pi,pi]
double ArcReduction2p(double arc);									//Calculates redu��o the first turn [0,2*pi]
//Fun��es for communicate with the Mathematica
Matrix List(double a, double b, double c);							//Returns a vector with these components the,b,c
double Power(double a, double b);
double Power(Matrix a, double b);
double Sin(double a);
double Cos(double a);
Matrix Dot(const Matrix &matrix1, const Matrix &matrix2);						//Product of matrices
double operator + (double a, const Matrix &matrix2);						//Operator Sum of matrix with a element and a vector