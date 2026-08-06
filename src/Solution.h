#pragma once
#include <stdio.h>

using namespace std;

class Solution
{
public:
	//Tracking solution - index and start/end times
	int solution_number;					//Solution number
	double start_time;						//Start time
	double end_time;						//End time

	Solution();
	virtual ~Solution();
	virtual bool Read(FILE *f) = 0;			//Reads input file
	virtual void Write(FILE *f) = 0;		//Writes output file
	virtual bool Solve() = 0;				//Solves solution routine

	void SetGlobalDOFs();					//Performs the numbering of the global degrees of freedom
	void DOFsActive();						//Checks and points the degrees of freedom active and inativos according to the elements, special constraints and constraints
	void MountLocal();						//Assembly of the elements and particles (information local)
	void MountElementLoads();				//Assembly of the loadings of field in elements (information local)
	void MountContacts();					//Assembly of the contacts (information local)
	void MountSpecialConstraints();			//Assembly of the special constraints (information local)
	void MountGlobal();						//Espalhamento of the information local in the matrices/vectors global
	void MountDisplacements();				//Inclui info of displacements impostos in the matrices/vectors global
	void MountLoads();						//Inclui info of loadings impostos in the matrices/vectors global
	void UpdateDisps();						//Updates the displacements nodal and Lagrange multipliers (Newton-Raphson)
	void Zeros();							//Zeroes variables iterativas for initial guess null
	void ZerosVelAccel();					//Zeroes velocities and acelerações for performing of analysis estatica (seguida of the dynamics)
	void Clear();							//Limpa the matrices global
	void SaveConfiguration();				//Saves configuration converged
	void RestoreConfiguration();			//Restaura the last configuration that converged
	void SetGlobalSize();					//Sets the size of the stiffness matrix global
	void PinballCheck();					//Varre contacts and checks pinball
	void BeginStepCheck();					//check initial of start of step for contacts
	void MountMass();						//Assembles mass matrix and internal forces inerciais
	void MountDamping(bool update_rayleigh);//Assembles damping matrix
	void MountDyn();						//Composicao of stiffness and internal forces for analysis dynamics (information local)
	void ComputeInitialConditions(bool zero_ICs);		//Computes the conditions inciais nodal impostas
	bool HaveErrors();						//Checks errors that impedem avança of the analysis, same in case of convergence
	void MountSparse();						//Assembles matrices esparsas
};

