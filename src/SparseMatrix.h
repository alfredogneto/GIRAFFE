#pragma once
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mkl.h"
#include "Eigen\SparseCore"
#include"Matrix.h"
//#include <unsupported/Eigen/SparseExtra>

class SparseMatrix
{
public:
	SparseMatrix();												//Construtor padrão
	SparseMatrix(int rows, int cols, int non_null_estimative);	//Construtor paramétrico
	SparseMatrix(SparseMatrix &copied);							//Construtor de cópia
	~SparseMatrix();
	void Clear();												//Zera coeficientes da matriz
	void setValue(int i, int j, double v);						//Seta na posição (i,j) da matriz quadrada o valor v
	void Mount();												//Seta matriz a partir da lista de triplets

	//Matrix esparsa criada via Eigen
	Eigen::SparseMatrix<double,1,int> m_matrix;					//matriz esparsa da biblioteca Eigen - RowMajor
	typedef Eigen::Triplet<double> T;							//cria typedef para a estrutura de triplets do Eigen
	std::vector<T> tripletList;									//cria vector de triplets (tripletList)
	bool mounted;												//Indica que a matriz está montada - conteúdo da tripletList foi transferido para a matriz, de fato
	int rows;													//Número de linhas
	int cols;													//Número de colunas
	int non_null_estimative;									//Estimativa de coeficientes não nulos

	void WriteMatrix(char* name);																				//Função de escrita da matriz em arquivo de texto
	friend Matrix sparsesystem(SparseMatrix &A, Matrix &b, int *info_fail, int processors, int solver_type);	//Resolve o sistema linear da forma Ax=b
	friend Matrix operator * (SparseMatrix& matrix1, Matrix& matrix2);											//Operador Multiplicacao de matrizes (matrix2 deve necessariamente ser um vetor)
	friend SparseMatrix operator + (SparseMatrix& matrix1, SparseMatrix& matrix2);								//Operador Soma de matrizes esparsas
	friend SparseMatrix operator * (double sigma, SparseMatrix &matrix1);										//Operador Multiplicacao de matriz esparsa por escalar
	SparseMatrix &operator = (SparseMatrix const &matrix1);														//Operador de Atribuição
	friend Matrix sparseeigen(SparseMatrix &K, SparseMatrix &M, Matrix &z, int n_e, bool eigenvectors, int* ret, double tolerance);			//Função que calcula autovalores de matriz esparsa - ARPACK - retorna em ret 0 se a execução foi correta
};

