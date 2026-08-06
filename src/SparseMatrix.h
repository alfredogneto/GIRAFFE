#pragma once
#include <Eigen/SparseCore>

class Matrix;
//#include <unsupported/Eigen/SparseExtra>

class SparseMatrix
{
public:
	SparseMatrix();												//Constructor padr�o
	SparseMatrix(int rows, int cols, int non_null_estimative);	//Parametric constructor
	SparseMatrix(const SparseMatrix &copied);							//Constructor of c�pia
	~SparseMatrix();
	void Clear();												//Zeroes coefficients of the matrix
	void setValue(int i, int j, double v);						//Sets in the posi��o (i,j) of the square matrix the value v
	void Mount() const;											//Sets matrix the partir of the list of triplets (lazy init - ver mutable below)

	//Matrix esparsa created via Eigen
	//m_matrix/mounted are "mutable" since Mount() does a initialization preguicosa (lazy) the partir
	//of the tripletList, the that needs funcionar same when the object and acessado via const&
	mutable Eigen::SparseMatrix<double,1,int> m_matrix;					//sparse matrix of the biblioteca Eigen - RowMajor
	typedef Eigen::Triplet<double> T;							//creates typedef for the structure of triplets of the Eigen
	std::vector<T> tripletList;									//creates vector of triplets (tripletList)
	mutable bool mounted;												//Indicates that the matrix this montada - content of the tripletList was transferido for the matrix, of fato
	int rows;													//Number of rows
	int cols;													//Number of columns
	int non_null_estimative;									//Estimativa of coefficients n�o null

	void WriteMatrix(char* name);																				//Fun��o of matrix writing in file of text
	friend Matrix sparsesystem(SparseMatrix &A, Matrix &b, int *info_fail, int processors, int solver_type);	//Solves the linear system of the way Ax=b
	friend Matrix operator * (const SparseMatrix& matrix1, const Matrix& matrix2);											//Operator Multiplication of matrices (matrix2 must necessariamente be a vector)
	friend SparseMatrix operator + (const SparseMatrix& matrix1, const SparseMatrix& matrix2);								//Operator Sum of matrices esparsas
	friend SparseMatrix operator * (double sigma, const SparseMatrix &matrix1);										//Operator Multiplication of sparse matrix by scalar
	SparseMatrix &operator = (SparseMatrix const &matrix1);														//Operator of Atribui��o
	friend Matrix sparseeigen(SparseMatrix &K, SparseMatrix &M, Matrix &z, int n_e, bool eigenvectors, int* ret, double tolerance);			//Fun��o that calculates eigenvalues of sparse matrix - ARPACK - returns in ret 0 if the execu��o was correct
};

