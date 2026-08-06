#pragma once
#include "Solution.h"
#include "Matrix.h"

class Modal:
	public Solution
{
public:
	Modal();
	~Modal();

	bool Read(FILE *f);											//Reading of input file
	void Write(FILE *f);										//Recording of output file
	bool Solve();												//Solves solution routine
	void CreateOutputFolder();									//Creates output folder
	void WriteResults(double time_value);						//Writes results
	void WriteMatrices();										//Writes the matrices in the folder modal
	void WriteModes();											//Writes files with the modos of vibrar
	void WriteVTK_XML();										//Writes modos of vibrar of all the elements

	//Functions to evaluate mode displacements
	void ComputeModalDisplacement(int node, int DOF, int mode, double* Re, double* Im);	//Calculates the value of the displacement modal in the node and degree of freedom in question
	double ComputeModeNorm(int mode);													//Computes the norm of the eigenvector
	
	//Input variables
	bool export_matrices;
	int number_modes;
	double tolerance;
	bool compute_eigenvectors;
	int number_frames;

	//Internal variables
	char copy_name[1000];										//Address of the directory for save files
	Matrix eig;													//Matrix with eigenvalues complexos
	Matrix eigenvectors;										//Matrix with eigenvectors complexos
	double mode_factor;											//Factor multiplicativo for the modo of vibrar
	bool concomitant_solution;									//Flag that indicates if and analysis concomitante
	
	//Special functions for matrices
	void ZerosLocalMatrices();				//Zeroes matrices of the elements (for use in the analysis modal)
	void MountMassModal();					//Assembles mass matrix for analysis modal
	void MountDampingModal();				//Assembles damping matrix for analysis modal
	void MountDynModal();					//Assemblies for analysis Modal

};

