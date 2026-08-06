#pragma once

class Matrix;

class LagrangeSaveLRAS
{
public:
	LagrangeSaveLRAS(int e_elements);
	~LagrangeSaveLRAS();

	int **contact_status;								//Indicates true if there is contact, and false if not there is. For control of the point of start of contact
	Matrix*** x0prev;									//Point previous of contact
	double **gt1s;										//Deslizamentos acumulados in the directions 1 and 2
	double **gt2s;										//Deslizamentos acumulados in the directions 1 and 2
	double **gn;										//gn
	double **Fx;										//Force total of contact in x
	double **Fy;										//Force total of contact in y
	double **Fz;										//Force total of contact in z

	int n_elements;
};

