#pragma once

class Matrix;

class NSContactData
{
public:
	NSContactData();
	~NSContactData();
	int n_solutions;
	double** convective;
	double** copy_convective;
	int* return_value;
	bool* repeated;

	double* g_n;								//gap normal
	double* copy_g_n;							//gap normal
	Matrix** g_t;								//gap tangencial - atual
	Matrix** copy_g_t;							//gap tangencial - c�pia
	Matrix** G_p;
	Matrix** t1_p;
	Matrix** t2_p;
	Matrix** n_p;
	Matrix** G_i;
	Matrix** G_ip;
	Matrix** t1_i;
	Matrix** t2_i;
	Matrix** n_i;
	bool alloced;								//Booleana que indica se esta ou n�o alocado - controle de aloca��o
	
	void CheckRepeated(double tol_coordinate_value);	//Checa repeti��o de raizes e salva a informa��o na matriz 'repeated'
	void Plot();
	void Alloc();										//Aloca matrizes
	void Free();										//Desaloca matrizes
};	
