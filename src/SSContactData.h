class Matrix;
class SurfaceData;

#pragma once
class SSContactData
{
public:
	SSContactData();
	~SSContactData();
	int n_solutions;
	double** convective;
	double** copy_convective;
	double** copy_deg_coordinates;				//degenerated coordinates copy
	double** initial_guess;						//initial guess to be used in SSMP
	
	bool** deg_control;							//control of degeneration
	Matrix** P;									//degeneration basis
	Matrix** P_0;								//Degenerative operator
	void MountDegenerativeOperator();			//Using info from deg_control and P, establishes P_0 by selecting appropriate columns

	int* characterization_index;
	int* copy_characterization_index;
	int* return_value;
	int* copy_return_value;
	bool* copy_convergedLCP;
	bool* convergedLCP;
	bool* repeated;
	double*** invHessian;						//Inverse of the Hessian matrix - determined during the FindMinimumParameters routine
	
	bool* copy_degenerated;						//true or false for indicate if the contact was or not degenerated in the step previous
	bool* degenerated;							//true or false for indicate if the contact and or not degenerated
	
	double* g_n;								//gap normal
	double* copy_g_n;							//gap normal
	Matrix** g_t;								//gap tangential - current
	Matrix** copy_g_t;							//gap tangential - copies
	Matrix** g;									//gap cross - current
	Matrix** copy_g;							//gap cross - copies
	Matrix ** n;								//contact normal direction
	Matrix ** copy_n;							//contact normal direction - copies
	bool* stick;								//boolean - indicates stick (true) or slide (false)
	bool* copy_stick;							//boolean - indicates stick (true) or slide (false) - copies
	
	SurfaceData* surf1;							//Data - surface 1
	SurfaceData* surf2;							//Data - surface 2

	bool alloced;								//Boolean that indicates if this or not allocated - control of allocation
	void CheckRepeated(double tol_coordinate_value);	//Checks repetição of raizes and saves the information in the matrix 'repeated'
	void Plot();
	void PlotSmallReport();
	void Alloc();										//Allocates matrices
	void Free();										//Deallocates matrices	
};	

