#pragma once
#include "Element.h"

class RigidBody_1 :
	public Element
{
public:
	RigidBody_1();
	~RigidBody_1();
	bool Check();													//Checks inconsistencies in the element for avoid errors of execution
	bool Read(FILE *f);												//Reading of the input file
	void Write(FILE *f);											//File writing of output
	void WriteResults(FILE *f);										//Writes file of results
	void WriteMonitor(FILE *f, bool first_record, double time);		//Writes in the monitor of the element
	void WriteVTK_XMLBase(std::vector<float> *float_vector);		//Writing VTK XML - used in the file base
	void WriteVTK_XMLRender(FILE *f);								//Writing VTK XML - used in the file render
	void Mount();													//Assembles contributions of the element for stiffness and residual
	void MountElementLoads();										//Assembles loadings associated to the element
	void MountMass();												//Assembles the mass matrix
	void MountMassModal();											//Assembles the mass matrix for performing of the analysis modal
	void MountDampingModal();										//Assembles the damping matrix for performing of the analysis modal
	void MountDamping(bool update_rayleigh);						//Assembles the damping matrix
	void MountDyn();												//Assemblies - Newmark
	void MountDynModal();											//Assemblies for analysis modal - insertion of the mass matrix and damping in the stiffness matrix for subsequent assembly global
	void TransformMatrix();											//Assembles transformation matrix of coordinates
	void MountGlobal();												//Fills the contribution of the element in the matrices global
	void Zeros();													//Zeroes some matrices used in the calculations
	void SaveLagrange();											//Saves data for configuration converged
	void PreCalc();													//Pre-calculation of variables that and done a single time in the start
	void InertialContributions();				                    //Calculates contributions of inertia of the Rigid Body
	void MountFieldLoading();											//Calculates the contribution of the internal forces of field (acceleration of the gravidade)

	//Variables of the element
	int RB_data_ID;													//ID of the parameters of the Rigid Body
	double v[1000];													//Variable necessary for the AceGen			
	double** Jr;													//Tensor of inertia - formato double**
	double* br;														//Vector br - formato double*
	double** DdT;													//Operator tangent
	double* dT;														//Residual
	double** Ddfield;												//Loading of field linearizado
	double* dfield;													//Residual loading of field
	
	double*	T1;														//Kinetic energy translacional
	double*	T2;														//Kinetic energy rotacional
	double*	T3;														//Kinetic energy acoplamento
	double*	T;														//Kinetic energy
	double* L;														//Momentum Linear 
	double* magL;													//Momentum Linear 
	double* HG;													    //Momentum Angular around Barycenter
	double* magHG;												    //Momentum Angular around Barycenter

	
	//Variables cinematicas - declaração
	double alphai[3];
	double ud[3];
	double alphad[3];
	double dui[3];
	double ddui[3];
	double omegai[3];
	double domegai[3];
		
};

