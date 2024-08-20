#pragma once
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "Matrix.h"
#include <math.h>

class SuperNode
{
public:
	SuperNode();
	~SuperNode();
	
	double super_node_coordinates[3];	//Input super node coordinates (reference configuration)

	double* ref_coordinates;	//Reference coordinates
	double* copy_coordinates;	//Last converged coordinates
	double* displacements;		//Displacements on DOFs
	double* vel;				//Velocidades on DOFs
	double* copy_vel;			//Velocidades on DOFs (copy of last converged)
	double* accel;				//Accelerations on DOFs 
	double* copy_accel;			//Accelerations on DOFs (copy of last converged)

	double* ref_temperatures;	//Reference temperatures
	double* copy_temperatures;	//Last converged temperatures
	double* temperatures;		//Temperatures

	bool* constraints;			//true:		constrained
								//false:	unconstrained
	int* DOFs;					//Number of global DOFs corresponding to all super node DOFs
	
	int n_DOFs;						//Number of DOFs (total)
	int n_displacement_DOFs;		//Number of DOFs of displacements
	int n_temperature_DOFs;			//Number of DOFs of temperatures
	//..other natures of DOFs

	int ID;						//SuperNode ID number

	bool alloced;
	void Alloc(int e_displacement_DOFs, int e_temperature_DOFs);		//Allocates the memory for the super node
	void Free();				//Frees the memory for the super node
	bool Read(FILE *f);			//Read function
	void Write(FILE *f);		//Write function
	void SaveConfiguration();	//Saves a converged configuration
};

