#pragma once
#include "SpecialConstraint.h"
#include "BoolTable.h"
#include "Matrix.h"

class NodalConstraintDOF:
	public SpecialConstraint
{
public:
	NodalConstraintDOF();
	~NodalConstraintDOF();

	bool Read(FILE* f);		//Reading
	void Write(FILE* f);	//Recording
	void Mount();			//Assembly of the residuals and stiffness tangent
	void MountGlobal();		//Fills the contribution of the element in the matrices global
	void PreCalc();			//Pre-calculation of variables that and done a single time in the start
	void SaveLagrange();	//Saving variables of the configuration converged
	void ActivateDOFs();	//Checks which Lagrange multipliers will be ativados,according to the activation of the GLs of the nodes of the which the special constraint participates
	bool Check();			//Checks inconsistencies in the SC for avoid errors of execution
	void ComputeVelAccel();		//Computes effect of the initial conditions in the nodes of the constraint
	void ComputeInitialGuessDisplacements();
	void WriteVTK_XMLRender(FILE* f);
	void WriteMonitor(FILE *f, bool first_record, double time);	//Writes in the monitor of the SpecialConstraint
	//Variables
	int node_A;
	int node_B;

	//BoolTables for each kind of DOF
	BoolTable UX_table;
	BoolTable UY_table;
	BoolTable UZ_table;
	BoolTable ROTX_table;
	BoolTable ROTY_table;
	BoolTable ROTZ_table;

	Matrix I6;
	//Stiffness matrix tangent and vector residual
	Matrix* stiffness;
	Matrix* residual;
	Matrix r;

};

