#pragma once
#include "Contact.h"
//Types of pairs of surface:

class SplineElement_SplineElement;
class SPContactData;
class SplineElementPair;
class Matrix;

////////Contact of the type SPSP - "Spline to Spline"///////////
class SPSP :
	public Contact
{
public:
	SPSP();
	~SPSP();
	int n_SP1;		//Spline ID 1
	int n_SP2;		//Spline ID 2
	double mus;		//Coef. of atrito estatico
	double mud;		//Coef. of atrito dinamico
	double ept;		//penalidade tangential
	double epn;		//penalidade normal
	double epn_n;	//potência of the penalidade normal
	double ct;		//coef. of dissipação tangential
	double cn;		//coef. of dissipação normal
	double pinball;	//pinball radius
	bool write_report;
	bool write_report_diverged;
	
	bool Check();												//Checks inconsistencies in the element for avoid errors of execution
	bool Read(FILE *f);											//Reading of the input file
	void Write(FILE *f);										//Recording of the output file
	void WriteResults(FILE *f);									//Writes file of results
	void WriteMonitor(FILE *f, bool first_record, double time);	//Writes in the monitor of the contact
	void WriteVTK_XMLRender(FILE *f);							//Plots contact surfaces
	void WriteVTK_XMLForces(FILE *f);							//Plots contact forces
	void SaveLagrange();										//Saves variables for description lagrangiana updated
	bool HaveErrors();											//Returns 1 - there is some error, same that tenha converged. Ex: penetração excessiva
	void MountDyn();											//Assemblies referring to the part dependente of velocities - dynamics
	void Mount();												//Assembly - contributions local of contact
	void MountGlobal();								//Fills the contribution of the contact in the matrices global
	void Band(int* band_fixed, int* band_free);		//Calculates the band generated in the matrix global through the contact
	void PreCalc();									//Pre-calculation of variables that and done a single time in the start
	void PinballCheck();							//Checks proximity between Surface to Surface
	void Alloc();									//Allocates in the memory
	void SetPairs();								//Assigns pairs of surfaces, according to sweep done between the surface sets 1 and 2										
	void BeginStepCheck();											//check initial of the contact  - start of each increment

	void AllocSpecific(int surface1_index, int surface2_index, int sol_index);				//Allocation specifies for the contact surface1, surface2
	void FreeSpecific(int surface1_index, int surface2_index, int sol_index);				//Liberação specifies for the contact surface1, surface2
	void EvaluateNormalGap(int surf1_index, int surf2_index, int sol_index);				//Calculates matrix with gap normal between surface1, surface2
	void ReportContact(int surface1_index, int surface2_index, int sol_index);				//Prints in the tela information of the contact pair surface1, surface2
	
	
	//Variables internal
	int number_surfaces1;														//Number of surfaces of the SurfaceSet 1
	int number_surfaces2;														//Number of surfaces of the SurfaceSet 2
	bool*** alloc_control;														//Control of allocation dynamics
	SPContactData*** cd;														//Data of the solutions of contact
	int number_pointwise_interactions;			//Maximum number of interações pontuais by pair surface-surface
	Matrix**** i_loading;						//Vector of internal forces internal
	double***** contact_stiffness;				//Stiffness matrix
	int * DOFs_surfaces1;						//Number of degrees of freedom that each a of the surfaces in the Surface Set has - influência in the vector of residuals and in the regidez tangent
	int * DOFs_surfaces2;						//Number of degrees of freedom that each a of the surfaces in the Surface Set has - influência in the vector of residuals and in the regidez tangent			
	Matrix* I3;									//Identity matrix of order 3
	SplineElementPair*** surf_pair;				//Pairs of surfaces for atribuições of functions especificas for each pair

	//Saving contact forces (for post-processing)
	Matrix**** fn;
	Matrix**** ft;

	//Tolerance for contact detection purposes (to avoid zero gaps at the begining of the simulation, for example)
	double pen_tol = 1e-14;
	bool specialLCP;
};

