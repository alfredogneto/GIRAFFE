#pragma once
#include "Element.h"

class TwoNodeConnector_1 :
	public Element
{
public:
	TwoNodeConnector_1();
	~TwoNodeConnector_1();

	bool Check();												//Checks inconsistencies in the element for avoid errors of execution
	bool Read(FILE *f);
	void Write(FILE *f);
	void WriteResults(FILE *f);									//Writes file of results
	void WriteMonitor(FILE *f, bool first_record, double time);	//Writes in the monitor of the element
	void WriteVTK_XMLBase(std::vector<float> *float_vector);
	void WriteVTK_XMLRender(FILE *f);
	void Mount();												//Assembles element
	void MountElementLoads();									//Assembles loadings associated to the element
	void MountMass();											//Assembles the mass matrix
	void MountMassModal();										//Assembles the mass matrix for performing of the analysis modal
	void MountDampingModal();									//Assembles the damping matrix for performing of the analysis modal
	void MountDamping(bool update_rayleigh);					//Assembles the damping matrix
	void MountDyn();											//Assemblies - Newmark
	void MountDynModal();										//Assemblies for analysis modal - insertion of the mass matrix and damping in the stiffness matrix for subsequent assembly global
	void TransformMatrix();										//Assembles transformation matrix of coordinates
	void MountGlobal();											//Fills the contribution of the element in the matrices global
	void Zeros();												//Zeroes some matrices used in the calculations

	void SaveLagrange();										//Saves variables in the gauss points useful for description lagrangiana updated
	void PreCalc();												//Pre-calculation of variables that and done a single time in the start

	//Variables of the element
	Matrix* K;													//Stiffness matrix 12x12
	Matrix* C;													//Damping matrix 12x12
	Matrix* Q;													//Transformation of coordinates 12x12

	Matrix* K_m;												//Stiffness matrix 6x6
	Matrix* C_m;												//Damping matrix 6x6
	Matrix* Q_m;												//Transformation of coordinates 6x6

	Matrix* f;													//Force for compor global
	Matrix* fd;													//Force - part of the damping
	Matrix* Kt;													//Operator tangent for compor global
	Matrix* Ct;													//Operator tangent - part of the damping
};