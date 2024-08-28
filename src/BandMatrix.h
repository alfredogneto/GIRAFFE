#pragma once
#include"Matrix.h"


class BandMatrix 
{
public:
	BandMatrix(void);
	BandMatrix(long band_s,long band_i,long columns);				//Construtor Parametrico
	BandMatrix(BandMatrix &copied);									//Construtor de cópia
	~BandMatrix(void);

	//General Functions
	void print();													//Imprime a matriz no console
	void fprint(char* s);											//Imprime a matriz em um arquivo de texto, cujo nome esta no char s
	void fprintSquare(char* s);										//Imprime a matriz em um arquivo de texto, cujo nome esta no char s (quadrada)
	bool alloc();													//Aloca a matriz
	bool flush();													//Libera a memória ocupada pela matriz
	void clear();													//Zera a matriz, mantando as dimensões atuais
	Matrix Full();													//Retorna a matriz quadrada (FULL)

	friend Matrix bandsystem(BandMatrix &A, Matrix &b,int *info_fail);	//Resolve o sistema linear da forma Ax=b
	friend Matrix operator * (BandMatrix &matrix1, Matrix &matrix2);		//Operador Multiplicacao de matrizes (matrix2 deve necessariamente ser um vetor)	
	friend BandMatrix operator * (double sigma, BandMatrix &matrix1);		//Operador Multiplicacao de matriz de banda por escalar
	friend BandMatrix operator + (BandMatrix &matrix1, BandMatrix &matrix2);//Operador soma de matrizes de banda

	double getValue(long i, long j);								//Retorno do valor na posição especificada - entrada i e j
																	//na posição da matriz quadrada original
	void setValue(long i, long j, double v);						//Seta na posição (i,j) da matriz quadrada o valor v

	int getBandi();													//Retorna a banda inferior da matriz
	int getBands();													//Retorna a banda superior da matriz
	int getColumns();												//Retorna o numero de colunas da matriz
	void setBandi(int bandi);										//Seta a banda inferior da matriz
	void setBands(int bands);										//Seta a banda superior da matriz
	void setColumns(int lines);										//Seta o numero de colunas da matriz
	double* getMatrix();											//Retorna o endereço de uma matriz

	long getAllocedLines();											//Retorna as linhas alocadas

protected:
	double*  m_matrix;												//Matriz unidimensional
	long     m_columns;												//Numero de linhas
	long     m_band_s;												//Banda superior da matriz
	long     m_band_i;												//Banda inferior da matriz
	long	 m_alloced_lines;										//Numero de linhas atualmente alocadas (pois m_matrix e unidimensional)
	bool	 m_lines_deleted;										//Flag que indica se as linhas foram desalocadas
};
