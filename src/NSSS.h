#pragma once
#include "Contact.h"

class Matrix;
class NSContactData;

////////Contact of the type NSSS - "Node Set to Surface Set"///////////
class NSSS:
	public Contact
{
public:
	NSSS();
	~NSSS();
	int n_NS;		//Node set ID
	int n_SS;		//Surface set ID
	double mu;		//Coef. of atrito
	double ept;		//penalidade tangential
	double epn;		//penalidade normal
	double cn;		//coef. of dissipação normal
	double ct;		//coef. of dissipação tangential
	double pinball;	//pinball radius
	double radius;

	bool Check();				//Checks inconsistencies in the element for avoid errors of execution
	bool Read(FILE *f);
	void Write(FILE *f);
	void WriteResults(FILE *f);									//Writes file of results
	void WriteMonitor(FILE *f, bool first_record, double time);	//Writes in the monitor of the contact
	void WriteVTK_XMLRender(FILE *f);							//Plots contact surfaces
	void WriteVTK_XMLForces(FILE *f);							//Plots contact forces
	void SaveLagrange();										//Saves variables for description lagrangiana updated
	bool HaveErrors();											//Returns 1 - there is some error, same that tenha converged. Ex: penetração excessiva
	void MountDyn();											//Assemblies - Newmark
	void Mount();
	void MountGlobal();								//Fills the contribution of the contact in the matrices global
	void Band(int* band_fixed, int* band_free);		//Calculates the band generated in the matrix global through the contact
	void PreCalc();									//Pre-calculation of variables that and done a single time in the start
	void PinballCheck();							//Checks proximity between each beam from LR to AS
	void Alloc();									//Allocates in the memory
	void BeginStepCheck();											//check initial of the contact  - start of each increment
	
	void AllocSpecific(int node_index, int surface_index, int sol_index);				//Allocation specifies for the contact node_index, surface_index
	void FreeSpecific (int node_index, int surface_index, int sol_index);				//Liberação specifies for the contact node_index, surface_index
	void EvaluateTangentialGap(int node_index, int surf_index, int sol_index);			//Calculates matrix with gap tangential

	void ReportContact(int node_index, int surf_index, int sol_index);					//Prints in the tela information of the contact pair
	int number_nodes;								//Number of nodes of the NodeSet
	int number_surfaces;							//Number of surfaces of the SurfaceSet

	bool*** alloc_control;							//Control of allocation dynamics
	
	//Variables internal
	int number_pointwise_interactions;			//Maximum number of interações pontuais by pair node-surface
	NSContactData*** cd;						//Data of the solutions of contact
	Matrix*** xS_p;
	Matrix*** QS;
	Matrix*** alphaS;
	Matrix*** uS;
	Matrix**** i_loading;						//Vector of internal forces internal
	double***** contact_stiffness;				//Stiffness matrix
	
	int * DOFs_surfaces;						//Number of degrees of freedom that each a of the surfaces in the Surface Set has - influência in the vector of residuals and in the regidez tangent			

	Matrix* I3;

	
};

