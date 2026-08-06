#pragma once

#include "BoolTable.h"

class SpecialConstraint
{
public:
	SpecialConstraint() {}
	virtual ~SpecialConstraint() {}

	virtual bool Read(FILE *f) = 0;		//Reading
	virtual void Write(FILE *f) = 0;	//Recording
	virtual void Mount() = 0;			//Assembly of the residuals and stiffness tangent
	virtual void MountGlobal() = 0;		//Fills the contribution of the element in the matrices global
	virtual void PreCalc() = 0;			//Pre-calculation of variables that and done a single time in the start
	virtual void SaveLagrange() = 0;	//Saving variables of the configuration converged
	virtual void ActivateDOFs() = 0;	//Performs activation of the GLs of the nodes of the which the special constraint participates
	virtual bool Check() = 0;			//Checks inconsistencies in the SC for avoid errors of execution
	virtual void ComputeVelAccel() = 0;	//Computes effect of the initial conditions in the nodes of the constraint
	virtual void WriteVTK_XMLRender(FILE *f) = 0;
	virtual void ComputeInitialGuessDisplacements() = 0;
	virtual void WriteMonitor(FILE *f, bool first_record, double time) = 0;	//Writes in the monitor of the SpecialConstraint
	int number;							//ID of the Special Constraint

	int* active_lambda;					//Control of the activation or not of the Lagrange multipliers, according to the activation or not of the GLs of the nodes associated to the constraint					
	//1 - GL active
	//0 - GL inactive
	int* GLs;							//Numbers of the GLs of the Lagrange multipliers
	int n_GL;							//Number of Lagrange multipliers
	double* lambda;						//Lagrange multipliers of the constraint
	double* copy_lambda;				//Lagrange multipliers of the constraint (last configuration converged)
	BoolTable bool_table;				//Bool table that active or desativa in the seq. of solutions
};

