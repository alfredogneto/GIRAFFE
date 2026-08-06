#pragma once
#include "Contact.h"
class Matrix;

////////Contact of the type LRLR - "Line region to line region"///////////
/////Desenvolvido for trabalhar in set with elements Pipe_1/////
class LRLR :
	public Contact
{
public:
	bool Check();				//Checks inconsistencies in the element for avoid errors of execution
	LRLR();
	~LRLR();
	int n_LR1;
	int n_LR2;
	double mu;
	double ept;
	double epn;
	double pinball;
	double sum_Fx;												//Variables that armazenam somas vetoriais of the internal forces (monitor)
	double sum_Fy;
	double sum_Fz;
	bool Read(FILE *f);
	void Write(FILE *f);
	void WriteResults(FILE *f);									//Writes file of results
	void WriteMonitor(FILE *f, bool first_record, double time);	//Writes in the monitor of the contact
	void WriteVTK_XMLRender(FILE *f);
	void WriteVTK_XMLForces(FILE *f);
	void SaveLagrange();										//Saves variables for description lagrangiana updated
	void Mount();
	void MountGlobal();											//Fills the contribution of the contact in the matrices global
	void Band(int* band_fixed, int* band_free);					//Calculates the band generated in the matrix global through the contact
	void PreCalc();												//Pre-calculation of variables that and done a single time in the start
	void PinballCheck();										//Checks proximity between each beam from LR to LR
	void BeginStepCheck();											//check initial of the contact  - start of each increment
	void Alloc(int e_elements1, int e_elements2);				//Allocates in the memory the variables that dependem of the number of elements
	bool HaveErrors();											//Returns 0 - not there is cruzamento 1 - there is cruzamento
	void MountDyn();												//Assemblies - Newmark
	//Functions internal
	void FillNodes(int e_element1, int e_element2);				//Fills the variables of the nodes with values atualizados
	int FindMinimumParameters(int i, int j);					//Determines csi1 and csi2 of the point of minimum distance
	void EvaluateParameters(int i, int j);						//Calculates parameters as a function of csi_1 and csi_2
	void CalculateLengthsAndTangents();							//Calculates the length of the elements //Calculates the vectors t1 and t2 of the elements
	void PlotContactStatus(int i, int j);						//Plots caracteristicas of the contact, independente of ter converged or not
	//Variables internal:
	int n_elements1;											//Number of elements of the line region corresponding to the contact
	int n_elements2;											//Number of elements of the line region corresponding to the contact
	int temp_element1;
	int temp_element2;
	int temp_node;
	double error;
	double tol_NR;												//Tolerência NR
	int max_it;													//Maximum number of iterations
	double tol_ortho;											//Tolerance of error the orthogonality
	double dot_test;
	double gt1, gt2;
	double l1, l2;
	double gte1, gte2;
	double Fat_max;
	double t1t2;
	double tt1;
	double tt2;
	double m1;
	double m2;
	double Fat_try;
	double Fat;
	double Fat1;
	double Fat2;
	double delta_lambda, delta_lambda_1, delta_lambda_2;
	double d;
	double N_1_1, N_3_1, dN_1_1, dN_3_1, ddN_1_1, ddN_3_1;		//Variables relacionadas the functions of way of the element 1
	double N_1_2, N_3_2, dN_1_2, dN_3_2, ddN_1_2, ddN_3_2;		//Variables relacionadas the functions of way of the element 2
	double D_;
	double dot11;
	double dot22;
	double dot12;
	Matrix* mean_position1;										//Position of the media aritmetica of the nodes of each element - calculated for perform the pinball search
	Matrix* mean_position2;										//Position of the media aritmetica of the nodes of each element - calculated for perform the pinball search
	Matrix* gte;
	Matrix* node_1_1;											//Coordinates of the node 1 - element 1
	Matrix* node_3_1;											//Coordinates of the node 3 - element 1
	Matrix* node_1_2;											//Coordinates of the node 1 - element 2
	Matrix* node_3_2;											//Coordinates of the node 3 - element 2
	Matrix* N1;
	Matrix* dN1;
	Matrix* ddN1;
	Matrix* N2;
	Matrix* dN2;
	Matrix* ddN2;
	Matrix* x1;													//Positions nodal in a single vector:
	Matrix* x2;													//Positions nodal in a single vector:
	Matrix* z1;
	Matrix* z2;
	Matrix* t1ext;
	Matrix* t2ext;
	Matrix* I361;
	Matrix* I362;
	Matrix* Ntio;
	Matrix* G1;
	Matrix* G1b;
	Matrix* G2;
	Matrix* G2b;
	Matrix* G3;
	Matrix* G3b;
	Matrix* aux1;
	Matrix* aux2;
	Matrix* G4;
	Matrix* G4b;
	Matrix* aux3;
	Matrix* G5;
	Matrix* aux4;
	Matrix* aux5;
	Matrix* aux6;
	Matrix* aux7;
	Matrix* Q1;
	Matrix* Q2;
	Matrix* Q3;
	Matrix* Q4;
	Matrix* G6;
	Matrix* G6b;
	Matrix* G7;
	Matrix* G7b;
	Matrix* I3;													//Identidade of order 3
	Matrix* n;
	Matrix* dz1;		
	Matrix* dz2;		
	Matrix* ddz1;	
	Matrix* ddz2;	
	Matrix* N_ext;
	Matrix* dN_ext;
	Matrix* d1;
	Matrix* d2;
	Matrix* E1;
	Matrix* E2;
	Matrix* ST;													//Contributions of the part of the atrito to the operator tangent
	Matrix* STb;												//Contributions of the part of the atrito to the operator tangent
	Matrix* R;													//Residual
	Matrix* A;													//Matrix Jacobiana - that tb. and the matrix A, once converged the NR
	Matrix* B;
	Matrix* C;
	Matrix* D;
	Matrix* E;
	Matrix* F;
	Matrix* G;
	Matrix* t1;													//Directions tangenciais of the elements 1 and 2
	Matrix* t2;
	Matrix* delta_csi;
	//Variables dependentes of the number of elements
	Matrix ***last_dif_pos;
	Matrix ***c_stiffness;										//Stiffness matrix
	Matrix ***c_loading;										//Vector of internal forces external
	Matrix*** z2z1;
	bool** first_mount;
	bool** flag_cross;
	bool** sticking;
	double** last_g_n;
	double** csi_1_0;
	double** csi_2_0;
	double** copy_csi_1_converged;
	double** copy_csi_2_converged;
	double** gt1s;
	double** gt2s;
	double** g_t1s_temp;
	double** g_t2s_temp;
	int** return_value;
	int** last_return_value;
	double** g_n;
	double* L1;													//Length
	double* r1;													//Radius
	double* L2;													//Length
	double* r2;													//Radius
	double** csi_1;
	double** csi_2;
	double** factor_1;
	double** factor_2;
	//Variables for check of the function Band
	int temp_band_free;
	int temp_band_fixed;
	int lowest_free_global_DOF;
	int highest_free_global_DOF;
	int lowest_fixed_global_DOF;
	int highest_fixed_global_DOF;
	bool free_marked;
	bool fixed_marked;
};

