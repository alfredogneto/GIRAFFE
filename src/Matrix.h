#pragma once
#include <stdio.h>
#include <cmath>
#include "arpack.h"
#include <complex>
#define QUARTIC_H_INCLUDED

#include "mkl.h"
#include "mkl_blas.h"
#include "mkl_types.h"
#include "mkl_lapack.h"
#include "mkl_dss.h"
#define PI 3.1415926535897932384626433832795

class Matrix
{
public:
	//Constructors and Destrutor
	Matrix(void);													//Construtor Padrão
	Matrix(long lines);												//Construtor de matriz coluna
	Matrix(long lines, long columns);								//Construtor Parametrico
	Matrix(const Matrix &copied);									//Construtor de cópia
	~Matrix(void);													//Destrutor Padrão

	//Gets and Sets
	long getLines();												//Retorna o numero de linhas da matriz
	long getColumns();												//Retorna o numero de colunas da matriz 

	void setLines(long value);										//Define o numero de linhas da matriz 
	void setColumns(long value);									//Define o numero de colunas da matriz 
	double* getMatrix();											//Retorna o endereço de uma matriz
	
	//General Functions
	void print();													//Imprime a matriz no console
	void fprint(char* s);											//Imprime a matriz em um arquivo de texto, cujo nome esta no char s
	bool alloc();													//Aloca a matriz
	bool flush();													//Libera a memória ocupada pela matriz
	void clear();													//Zera a matriz, mantando as dimensões atuais

	//Operators
	void MatrixToPtr(double** ptr, int order);						//Salva ponteiro double** em ptr
	void PtrToMatrix(double** ptr, int order);						//Salva na matrix o conteudo do double**
	void PtrToMatrix(double** ptr, int lines, int columns);			//Salva na matrix o conteudo do double**
	double &operator() (long line, long column);					//Retorno do valor na posição especificada
	Matrix &operator = (Matrix const &matrix1);						//Operador de Atribuição

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	double*  m_matrix;												//Matriz unidimensional
	long     m_lines;												//Numero de linhas
	long     m_columns;												//Numero de colunas
	long	 m_alloced_lines;										//Numero de linhas atualmente alocadas
	bool	 m_lines_deleted;										//Flag que indica se as linhas foram desalocadas				
};

Matrix operator + (Matrix &matrix1, Matrix &matrix2);				//Operador Soma
Matrix operator - (Matrix &matrix1, Matrix &matrix2);				//Operador Subtração
Matrix operator * (Matrix &matrix1, Matrix &matrix2);				//Operador Multiplicacao de matrizes
Matrix operator * (double escalar, Matrix &matrix1);				//Operador Multiplicacao por escalar
Matrix operator * (Matrix &matrix1, double escalar);				//Operador Multiplicacao por escalar
bool operator == (Matrix &matrix1, Matrix &matrix2);				//Verificação de igualdade
bool operator != (Matrix &matrix1, Matrix &matrix2);				//Verificação de inegualdade
double dot(Matrix &matrix1, Matrix &matrix2);						//Produto escalar entre dois vetores
Matrix cross(Matrix &matrix1, Matrix &matrix2);						//Operador produto vetorial entre duas matrizes
Matrix dyadic(Matrix &matrix1, Matrix &matrix2);					//Operador produto vetorial entre duas matrizes
Matrix skew(Matrix &matrix1);										//Operador skew de um vetor
Matrix axial(Matrix &matrix1);										//Operador axial de um vetor
Matrix fullsystem(Matrix &A, Matrix &b, int *flag_error);			//Resolve o sistema linear da forma Ax=b
double norm(Matrix &matrix1);										//Retorna a norma de um vetor
double norm4(Matrix &matrix1);										//Retorna a norma de um vetor considerando somente os 4 primeiros graus de liberdade
Matrix transp(Matrix &matrix1);										//Retorna a transposta de uma matriz
void zeros(Matrix* matrix1);										//Zera a matriz
Matrix invert2x2(Matrix &matrix);									//Inverte uma matriz 2x2
Matrix invert3x3(Matrix &matrix);									//Inverte uma matriz 3x3
Matrix invert4x4(Matrix &matrix);									//Inverte uma matriz 4x4
Matrix invert5x5(Matrix &matrix);									//Inverte uma matriz 5x5
Matrix invert6x6(Matrix &matrix);									//Inverte uma matriz 6x6	
Matrix invert(Matrix &matrix);										//Inverte uma matriz de 2x2 a 6x6 escolhendo automaticamente a função correta

int fulleigen1(Matrix &A, Matrix &P, Matrix &D, double abstol);		//Calcula os autovalores e autovetores de uma matriz simetrica
int fulleigen2(Matrix &A, Matrix &P, Matrix &D);					//Calcula os autovalores e autovetores de uma matriz simetrica
double mineigen(Matrix &A, Matrix &P, Matrix &D, double abstol);	//Calcula o menor autovalor de uma matriz simetrica
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Matrix V(Matrix x, Matrix t, double alpha_escalar);					//Função para o calculo do operador V (para montagem da matriz de rigidez geometrica)
Matrix d_V(Matrix x, Matrix d_x, Matrix t, double alpha_escalar);	//Função para o calculo do operador d_V (para montagem da matriz de rigidez geometrica)

double ArcReduction(double arc);									//Calcula redução a primeira volta [-pi,pi]
double ArcReduction2p(double arc);									//Calcula redução a primeira volta [0,2*pi]
//Funções para conversar com o Mathematica
Matrix List(double a, double b, double c);							//Retorna um vetor com esses componentes a,b,c
double Power(double a, double b);
double Power(Matrix a, double b);
double Sin(double a);
double Cos(double a);
Matrix Dot(Matrix &matrix1, Matrix &matrix2);						//Produto de matrizes
double operator + (double a, Matrix &matrix2);						//Operador Soma de matriz com um elemento e um vetor


//Marina
//const double PI = 3.141592653589793238463L;
const double M_2PI = 2 * PI;
const double eps = 1e-12;
typedef std::complex<double> DComplex;
//Polinômio de 3o grau
//---------------------------------------------------------------------------
// x - array of size 3
// In case 3 real roots: => x[0], x[1], x[2], return 3
//         2 real roots: x[0], x[1],          return 2
//         1 real root : x[0], x[1] ± i*x[2], return 1
unsigned int solveP3(double* x, double a, double b, double c);

//Polinômio 4o grau
//---------------------------------------------------------------------------
// solve quartic equation x^4 + a*x^3 + b*x^2 + c*x + d
// Attention - this function returns dynamically allocated array. It has to be released afterwards.
DComplex* solve_quartic(double a, double b, double c, double d);