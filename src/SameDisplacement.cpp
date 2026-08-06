#include "SameDisplacement.h"

#include "Node.h"
#include "InitialCondition.h"
#include "Database.h"
//Global variables
extern
Database db;

SameDisplacement::SameDisplacement()
{
	
	n_GL = 3;						//Três degrees of freedom (this constraint has 3 Lagrange multipliers)
	active_lambda = new int[n_GL];
	lambda = new double[n_GL];
	copy_lambda = new double[n_GL];
	GLs = new int[n_GL];
	node_A = 0;
	node_B = 0;

	for (int i = 0; i < n_GL;i++)
	{
		active_lambda[i] = 0;
		GLs[i] = 0;
		//Initial guess for the lambdas: values null
		lambda[i] = 0.0;
		copy_lambda[i] = 0.0;
	}
	
	
	
	I3 = Matrix(3, 3);
	I3(0, 0) = 1.0;
	I3(1, 1) = 1.0;
	I3(2, 2) = 1.0;
	r = Matrix(3);
	//Stiffness matrix tangent and vector residual
	stiffness = new Matrix(9,9);
	residual = new Matrix(9,1);

}

SameDisplacement::~SameDisplacement()
{
	delete []active_lambda;
	delete []GLs;
	delete []lambda;
	delete []copy_lambda;
	delete stiffness;
	delete residual;
}

//Reading
bool SameDisplacement::Read(FILE *f)
{
	char s[1000];
	fscanf(f, "%s", s);
	number = atoi(s);

	fscanf(f, "%s", s);
	if (!strcmp(s, "Nodes"))
	{
		fscanf(f, "%s", s);
		node_A = atoi(s);

		fscanf(f, "%s", s);
		node_B = atoi(s);
	}
	else
		return false;
	//Saves the position (stream)
	fpos_t pos;
	fgetpos(f, &pos);
	fscanf(f, "%s", s);
	if (!strcmp(s, "BoolTable"))
		bool_table.Read(f);
	else
	{
		fsetpos(f, &pos);
		bool_table.SetDefault(true);
	}

	return true;
}

//Recording
void SameDisplacement::Write(FILE *f)
{
	fprintf(f, "SameDisplacement\t%d\tNodes\t%d\t%d\n",
		number,
		node_A,
		node_B);
}

//Writes in the monitor of the SpecialConstraint
void SameDisplacement::WriteMonitor(FILE *f, bool first_record, double time)
{

}

void SameDisplacement::WriteVTK_XMLRender(FILE *f)
{

}

//Checks inconsistencies in the SC for avoid errors of execution
bool SameDisplacement::Check()
{
	if (node_A > db.number_nodes)
		return false;
	if (node_B > db.number_nodes)
		return false;
	//Check of the initial conditions
	int temp_node = 0;
	for (int i = 0; i < db.number_IC; i++)
	{
		temp_node = db.IC[i]->node;
		if (node_B == temp_node)
		{
			db.myprintf("Warning in Special Constraint %d.\nInitial Condition %d was prescribed to node %d (slave), leading to ignoring some of its components.\n", number, db.IC[i]->number, db.IC[i]->node);
		}
	}
	return true;
}

//Assembly of the residuals and stiffness tangent
void SameDisplacement::Mount()
{
	//In this constraint the stiffness tangent not if modifies never. and montada in the PreCalc()
	//If the constraint is active, performs the assembly
	if (active_lambda[0] == 1 && active_lambda[1] == 1 && active_lambda[2] == 1)
	{
		for (int i = 0; i < 3; i++)
			r(i, 0) = db.nodes[node_A - 1]->displacements[i] - db.nodes[node_B - 1]->displacements[i];

		//Update of the vector of residuals
		(*residual)(0, 0) = lambda[0];
		(*residual)(1, 0) = lambda[1];
		(*residual)(2, 0) = lambda[2];

		(*residual)(3, 0) = -lambda[0];
		(*residual)(4, 0) = -lambda[1];
		(*residual)(5, 0) = -lambda[2];

		(*residual)(6, 0) = r(0, 0);
		(*residual)(7, 0) = r(1, 0);
		(*residual)(8, 0) = r(2, 0);
	}
}

//Fills the contribution of the element in the matrices global
void SameDisplacement::MountGlobal()
{
	//Temporary variables for save the indexing global of the degrees of freedom the to be set in the stiffness matrix global
	int GL_global_1 = 0;
	int GL_global_2 = 0;
	double anterior = 0;
	//If the constraint is active, performs the assembly
	if (active_lambda[0] == 1 && active_lambda[1] == 1 && active_lambda[2] == 1)
	{
		for (int i = 0; i < 9; i++)
		{
			//////////////ASSEMBLY OF THE VECTOR OF INTERNAL FORCES UNBALANCED//////////////////
			//Takes of the global-DOF vector, the indexing of each global degree of freedom
			if (i < 3)//uA
				GL_global_1 = db.nodes[node_A - 1]->GLs[i];
			else
			{
				if (i<6)//uB
					GL_global_1 = db.nodes[node_B - 1]->GLs[i - 3];
				else
				{
					//lambda
					GL_global_1 = GLs[i - 6];
				}
			}

			//Case the degree of freedom be free:
			if (GL_global_1 > 0)
			{
				anterior = db.global_P_A(GL_global_1 - 1, 0);
				db.global_P_A(GL_global_1 - 1, 0) = anterior + (*residual)(i, 0);
				anterior = db.global_I_A(GL_global_1 - 1, 0);
				db.global_I_A(GL_global_1 - 1, 0) = anterior + (*residual)(i, 0);
			}
			else
			{
				anterior = db.global_P_B(-GL_global_1 - 1, 0);
				db.global_P_B(-GL_global_1 - 1, 0) = anterior + (*residual)(i, 0);
			}
			for (int j = 0; j < 9; j++)
			{
				//////////////////////ASSEMBLY OF THE STIFFNESS MATRIX/////////////////////////
				//Takes of the global-DOF vector, the indexing of each global degree of freedom
				if (j < 3)//uA
					GL_global_2 = db.nodes[node_A - 1]->GLs[j];
				else
				{
					if (j<6)//uB
						GL_global_2 = db.nodes[node_B - 1]->GLs[j - 3];
					else
					{
						//lambda
						GL_global_2 = GLs[j - 6];
					}
				}

				//Case the degrees of freedom be both free (Matrix Kaa)
				if (GL_global_1 > 0 && GL_global_2 > 0)
					db.global_stiffness_AA.setValue(GL_global_1 - 1, GL_global_2 - 1, (*stiffness)(i, j));
				//Case the degrees of freedom be both fixed (Matrix Kbb)
				if (GL_global_1 < 0 && GL_global_2 < 0)
					db.global_stiffness_BB.setValue(-GL_global_1 - 1, -GL_global_2 - 1, (*stiffness)(i, j));
				//Case the degrees of freedom be free and fixed (Matrix Kab)
				if (GL_global_1 > 0 && GL_global_2 < 0)
					db.global_stiffness_AB.setValue(GL_global_1 - 1, -GL_global_2 - 1, (*stiffness)(i, j));
				//Case the degrees of freedom be fixed and free (Matrix Kba)
				if (GL_global_1 < 0 && GL_global_2 > 0)
					db.global_stiffness_BA.setValue(-GL_global_1 - 1, GL_global_2 - 1, (*stiffness)(i, j));
			}
		}
	}
}

void SameDisplacement::ComputeInitialGuessDisplacements()
{
	//If for the step of creation of the constraint initializes initial conditions
	if (bool_table.GetAt(db.current_solution_number - 1) == true)
	{
		//Iterates over GL and sets initial conditions
		for (int i = 0; i < 3; i++)
		{
			if (active_lambda[i] == 1)
			{
				db.nodes[node_B - 1]->displacements[i] = db.nodes[node_A - 1]->displacements[i];
			}

		}
	}
}

//Computes effect of the initial conditions in the nodes of the constraint
void SameDisplacement::ComputeVelAccel()
{
	//If for the step of creation of the constraint initializes initial conditions
	if (bool_table.GetAt(db.current_solution_number - 1))
	{
		//Iterates over GL and sets initial conditions
		for (int i = 0; i < 3; i++)
		{
			if (active_lambda[i] == 1)
			{
				db.nodes[node_B - 1]->vel[i] = db.nodes[node_A - 1]->vel[i];
				db.nodes[node_B - 1]->accel[i] = db.nodes[node_A - 1]->accel[i];
			}

		}
	}

	
	
}

//Pre-calculation of variables that and done a single time in the start
void SameDisplacement::PreCalc()
{
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			(*stiffness)(i, j + 6) = I3(i, j);
			(*stiffness)(i + 3, j + 6) = -I3(i, j);
			(*stiffness)(i + 6, j) = I3(i, j);
			(*stiffness)(i + 6, j + 3) = -I3(i, j);
		}
	}
}

//Saving variables of the configuration converged
void SameDisplacement::SaveLagrange()
{
	//PrintPtr(lambda, 3);
	for (int i = 0; i < n_GL; i++)
		copy_lambda[i] = lambda[i];
}

//Checks which Lagrange multipliers will be ativados,according to the activation of the GLs of the nodes of the which the special constraint participates
void SameDisplacement::ActivateDOFs()
{
	//Active GLs of translation of the nodes A and B
	for (int i = 0; i < 3; i++)
		db.nodes[node_A - 1]->active_GL[i] = 1;
	for (int i = 0; i < 3; i++)
		db.nodes[node_B - 1]->active_GL[i] = 1;

	if (bool_table.GetAt(db.current_solution_number - 1))
	{
		
		//Active GLs of Lagrange multipliers
		active_lambda[0] = 1;
		active_lambda[1] = 1;
		active_lambda[2] = 1;
	}
	else
	{
		active_lambda[0] = 0;
		active_lambda[1] = 0;
		active_lambda[2] = 0;
	}
}
