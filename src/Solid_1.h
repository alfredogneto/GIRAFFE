#pragma once
#include "Element.h"

class Solid_1 :
	public Element
{
public:
	Solid_1();
	~Solid_1();
	bool Check();				//Checks inconsistencies in the element for avoid errors of execution
	bool Read(FILE *f);
	void Write(FILE *f);
	void WriteResults(FILE *f);//Writes file of results
	void WriteMonitor(FILE *f, bool first_record, double time);	//Writes in the monitor of the element
	void WriteVTK_XMLBase(std::vector<float> *float_vector);
	void WriteVTK_XMLRender(FILE *f);
	void Mount();			//Assembles elements
	void MountElementLoads();			//Assembles loadings associated to the element
	void MountMass();					//Assembles the mass matrix
	void MountMassModal();				//Assembles the mass matrix for performing of the analysis modal
	void MountDampingModal();			//Assembles the damping matrix for performing of the analysis modal
	void MountDamping(bool update_rayleigh);				//Assembles the damping matrix
	void MountDyn();					//Assemblies - Newmark
	void MountDynModal();				//Assemblies for analysis modal - insertion of the mass matrix and damping in the stiffness matrix for subsequent assembly global
	void TransformMatrix();	//Assembles transformation matrix of coordinates
	void MountGlobal();		//Fills the contribution of the element in the matrices global
	void Zeros();			//Zeroes matrices local of the element

	void SaveLagrange();	//Saves variables in the gauss points useful for description lagrangiana updated
	void PreCalc();		//Pre-calculation of variables that and done a single time in the start
};

