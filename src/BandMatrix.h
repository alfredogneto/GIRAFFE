#pragma once
#include"Matrix.h"


class BandMatrix 
{
public:
	BandMatrix(void);
	BandMatrix(long band_s,long band_i,long columns);				//Parametric constructor
	BandMatrix(BandMatrix &copied);									//Constructor of c�pia
	~BandMatrix(void);

	//General Functions
	void print();													//Prints the matrix in the console
	void fprint(char* s);											//Prints the matrix in a file of text, whose name this in the char s
	void fprintSquare(char* s);										//Prints the matrix in a file of text, whose name this in the char s (quadrada)
	bool alloc();													//Allocates the matrix
	bool flush();													//Frees the mem�ria occupied through the matrix
	void clear();													//Zeroes the matrix, keeping the dimens�es current
	Matrix Full();													//Returns the square matrix (FULL)

	friend Matrix bandsystem(BandMatrix &A, Matrix &b,int *info_fail);	//Solves the linear system of the way Ax=b (mutates A and b in-place via LAPACK)
	friend Matrix operator * (const BandMatrix &matrix1, const Matrix &matrix2);		//Operator Multiplication of matrices (matrix2 must necessariamente be a vector)
	friend BandMatrix operator * (double sigma, const BandMatrix &matrix1);		//Operator Multiplication of band matrix by scalar
	friend BandMatrix operator + (const BandMatrix &matrix1, const BandMatrix &matrix2);//Operator sum of matrices of band

	double getValue(long i, long j) const;							//Return of the value in the posi��o specified - input i and j
																	//in the posi��o of the square matrix original
	void setValue(long i, long j, double v);						//Sets in the posi��o (i,j) of the square matrix the value v

	int getBandi() const;											//Returns the band lower of the matrix
	int getBands() const;											//Returns the band upper of the matrix
	int getColumns() const;											//Returns the number of columns of the matrix
	void setBandi(int bandi);										//Sets the band lower of the matrix
	void setBands(int bands);										//Sets the band upper of the matrix
	void setColumns(int lines);										//Sets the number of columns of the matrix
	double* getMatrix() const;										//Returns the endere�o of a matrix

	long getAllocedLines() const;											//Returns the rows allocated

protected:
	double*  m_matrix;												//Matrix one-dimensional
	long     m_columns;												//Number of rows
	long     m_band_s;												//Band upper of the matrix
	long     m_band_i;												//Band lower of the matrix
	long	 m_alloced_lines;										//Number of rows currently allocated (since m_matrix and one-dimensional)
	bool	 m_lines_deleted;										//Flag that indicates if the rows were deallocated
};
