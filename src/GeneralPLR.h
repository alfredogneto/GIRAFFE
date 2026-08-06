#pragma once
#include "Contact.h"
#include <stdio.h>

#include "Matrix.h"

class GeneralPLR :
	public Contact
{
public:
	GeneralPLR();
	~GeneralPLR();

	bool Check();				//Checks inconsistencies in the element for avoid errors of execution
	void Mount();													//Assembles contacts
	void MountGlobal();												//Fills the contribution of the contact in the matrices global
	bool Read(FILE *f);
	void Write(FILE *f);
	void WriteResults(FILE *f);										//Writes file of results
	void WriteMonitor(FILE *f, bool first_record, double time);		//Writes in the monitor of the contact
	void PreCalc();													//Pre-calculation of variables that and done a single time in the start
	void PinballCheck();											//Checks proximity between each beam from LR to AS
	void BeginStepCheck();											//check initial of the contact  - start of each increment
	void SaveLagrange();											//Saves variables for description lagrangiana updated
	bool HaveErrors();												//Returns 1 - there is some error, same that tenha converged. Ex: penetração excessiva
	void Band(int* band_fixed, int* band_free);						//Calculates the band generated in the matrix global through the contact
	void MountDyn();												//Assemblies - Newmark
	void Alloc();													//Allocates structure of matrices for endereçar vectors and matrices devidos to the contact
	void AllocSpecific(int i, int j);								//Allocates matrices for contact specific
	void FreeSpecific(int i, int j);								//Allocates matrices for contact specific
	void WriteVTK_XMLRender(FILE *f);
	void WriteVTK_XMLForces(FILE *f);
	//Variables internal
	double mu, epn, ept, pinball;
	double c;														//Coefficient of dissipação of energy
	Matrix z1, z2;													//Positions centrais of the particles in the PinballSearch
	int n_particles;												//Number of particles
	int n_elements;													//Number of elements inside of the line region
	int n_LR;														//Number (indice) of the line region associated to the contact
	Matrix ****c_stiffness;											//Stiffness matrix
	Matrix ****c_damping;											//Damping matrix
	Matrix ****c_loading;											//Vector of internal forces external
	bool** alloc_control;											//Control of allocation dynamics
	Matrix I3;
};

