#pragma once
#include <stdio.h>
#include "BoolTable.h"

class Contact
{
public:
	Contact() {}
	virtual ~Contact() {}

	char* type_name;											//Name of the type of the element
	int number;													//Number 
	bool** activate;											//Based on pinball check, contains 1 or 0 info to activate/not the contact
	bool* typeOK1;												//Contains 1 or 0, based on verification if the type of the underlying element is compatible with the kind of contact
	bool* typeOK2;												//Contains 1 or 0, based on verification if the type of the underlying element is compatible with the kind of contact
	BoolTable bool_table;										//Bool table that active or desativa in the seq. of solutions
	virtual void Mount() = 0;									//Assembles contacts
	virtual void MountGlobal() = 0;								//Fills the contribution of the contact in the matrices global
	virtual void Band(int* band_fixed, int* band_free) = 0;		//Calculates the band generated in the matrix global through the contact

	virtual bool Read(FILE *f) = 0;
	virtual void Write(FILE *f) = 0;
	virtual void WriteResults(FILE *f) = 0;										//Writes file of results
	virtual void WriteMonitor(FILE *f, bool first_record, double time) = 0;		//Writes in the monitor of the contact
	virtual void WriteVTK_XMLRender(FILE *f) = 0;
	virtual void WriteVTK_XMLForces(FILE *f) = 0;
	virtual void PreCalc() = 0;													//Pre-calculation of variables that and done a single time in the start
	virtual void PinballCheck() = 0;											//Checks proximity between each beam from LR to AS
	virtual void BeginStepCheck() = 0;											//check initial of the contact  - start of each increment
	virtual void SaveLagrange() = 0;											//Saves variables for description lagrangiana updated
	virtual bool HaveErrors() = 0;												//Returns 1 - there is some error, same that tenha converged. Ex: penetração excessiva
	virtual void MountDyn() = 0;												//Assemblies - Newmark
	virtual bool Check() = 0;													//Checks inconsistencies in the element for avoid errors of execution
};

