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
	Matrix** g_t;								//gap tangential - current
	Matrix** copy_g_t;							//gap tangential - copies
	Matrix** G_p;
	Matrix** t1_p;
	Matrix** t2_p;
	Matrix** n_p;
	Matrix** G_i;
	Matrix** G_ip;
	Matrix** t1_i;
	Matrix** t2_i;
	Matrix** n_i;
	bool alloced;								//Boolean that indicates if this or not allocated - control of allocation
	
	void CheckRepeated(double tol_coordinate_value);	//Checks repetição of raizes and saves the information in the matrix 'repeated'
	void Plot();
	void Alloc();										//Allocates matrices
	void Free();										//Deallocates matrices
};	
