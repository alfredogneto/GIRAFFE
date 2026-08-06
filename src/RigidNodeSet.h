#pragma once
#include "SpecialConstraint.h"

#include "Matrix.h"
class RigidNodeSet :
	public SpecialConstraint
{
public:
	RigidNodeSet();
	~RigidNodeSet();

	bool Read(FILE *f);		//Reading
	void Write(FILE *f);	//Recording
	void Mount();			//Assembly of the residuals and stiffness tangent
	void MountGlobal();		//Fills the contribution of the element in the matrices global
	void PreCalc();			//Pre-calculation of variables that and done a single time in the start
	void SaveLagrange();	//Saving variables of the configuration converged
	void ActivateDOFs();	//Checks which Lagrange multipliers will be ativados,according to the activation of the GLs of the nodes of the which the special constraint participates
	void Alloc();			//Allocates matrices
	void EvaluateM2M3(double* v, double** M2, double** M3, double** Qp, double* alphap, double* di, double* lamb);
	bool Check();			//Checks inconsistencies in the SC for avoid errors of execution
	void ComputeVelAccel();		//Computes effect of the initial conditions in the nodes of the constraint
	void ComputeInitialGuessDisplacements();
	void WriteVTK_XMLRender(FILE *f);
	void WriteMonitor(FILE *f, bool first_record, double time);	//Writes in the monitor of the SpecialConstraint
	//Variables
	int pilot_node;			//Node piloto of the rigid body
	int node_set_ID;		//Number of the node set vinculado to the node piloto
	int n_nodes_set;		//Number of nodes of the set vinculado to the node piloto
	//Flag - activation of GL of rotation
	bool flag_rotation;
	
	Matrix I3;
	//Stiffness matrix tangent and vector residual
	Matrix** stiffness1;
	Matrix** residual1;
	Matrix** stiffness2;
	Matrix** residual2;
	Matrix** distancei;

	double v[600];
};

