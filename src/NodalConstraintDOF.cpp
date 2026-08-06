#include "NodalConstraintDOF.h"

#include "Node.h"
#include "InitialCondition.h"


#include "IO.h"

NodalConstraintDOF::NodalConstraintDOF()
{
	number = 0;
	//setting default values for constraints (false) - in case of not reading
	UX_table.SetDefault(false);
	UY_table.SetDefault(false);
	UZ_table.SetDefault(false);
	ROTX_table.SetDefault(false);
	ROTY_table.SetDefault(false);
	ROTZ_table.SetDefault(false);

	n_GL = 6;						//Três degrees of freedom (this constraint has 3 Lagrange multipliers)
	active_lambda = new int[n_GL];
	lambda = new double[n_GL];
	copy_lambda = new double[n_GL];
	GLs = new int[n_GL];
	node_A = 0;
	node_B = 0;

	for (int i = 0; i < n_GL; i++)
	{
		active_lambda[i] = 0;
		GLs[i] = 0;
		//Initial guess for the lambdas: values null
		lambda[i] = 0.0;
		copy_lambda[i] = 0.0;
	}

	I6 = Matrix(6, 6);
	I6(0, 0) = 1.0;
	I6(1, 1) = 1.0;
	I6(2, 2) = 1.0;
	I6(3, 3) = 1.0;
	I6(4, 4) = 1.0;
	I6(5, 5) = 1.0;
	r = Matrix(6);
	//Stiffness matrix tangent and vector residual
	stiffness = new Matrix(18, 18);
	residual = new Matrix(18, 1);

}

NodalConstraintDOF::~NodalConstraintDOF()
{
	delete[]active_lambda;
	delete[]GLs;
	delete[]lambda;
	delete[]copy_lambda;
	delete stiffness;
	delete residual;
}

bool NodalConstraintDOF::Read(FILE* f)
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

	bool DOF_OK = false;
	bool flag_continue = true;
	fpos_t pos;
	while (flag_continue == true)
	{
		DOF_OK = false;
		TryComment(f);
		fgetpos(f, &pos);
		if (fscanf(f, "%s", s) == EOF)
			return true;
		if (!strcmp(s, "UX"))
		{
			fscanf(f, "%s", s);
			if (!strcmp(s, "BoolTable"))
				UX_table.Read(f);
			else
				return false;
			DOF_OK = true;
		}
		if (!strcmp(s, "UY"))
		{
			fscanf(f, "%s", s);
			if (!strcmp(s, "BoolTable"))
				UY_table.Read(f);
			else
				return false;
			DOF_OK = true;
		}
		if (!strcmp(s, "UZ"))
		{
			fscanf(f, "%s", s);
			if (!strcmp(s, "BoolTable"))
				UZ_table.Read(f);
			else
				return false;
			DOF_OK = true;
		}
		if (!strcmp(s, "ROTX"))
		{
			fscanf(f, "%s", s);
			if (!strcmp(s, "BoolTable"))
				ROTX_table.Read(f);
			else
				return false;
			DOF_OK = true;
		}
		if (!strcmp(s, "ROTY"))
		{
			fscanf(f, "%s", s);
			if (!strcmp(s, "BoolTable"))
				ROTY_table.Read(f);
			else
				return false;
			DOF_OK = true;
		}
		if (!strcmp(s, "ROTZ"))
		{
			fscanf(f, "%s", s);
			if (!strcmp(s, "BoolTable"))
				ROTZ_table.Read(f);
			else
				return false;
			DOF_OK = true;
		}
		if (DOF_OK == true)
			flag_continue = true;
		else
		{
			fsetpos(f, &pos);
			flag_continue = false;
		}
	}
	return true;
}

void NodalConstraintDOF::Write(FILE* f)
{
	fprintf(f, "NodalConstraintDOF\t%d\tNodes\t%d\t%d\n", number, node_A, node_B);
	fprintf(f, "UX\t");
	UX_table.Write(f);
	fprintf(f, "UY\t");
	UY_table.Write(f);
	fprintf(f, "UZ\t");
	UZ_table.Write(f);
	fprintf(f, "ROTX\t");
	ROTX_table.Write(f);
	fprintf(f, "ROTY\t");
	ROTY_table.Write(f);
	fprintf(f, "ROTZ\t");
	ROTZ_table.Write(f);
}

//Writes in the monitor of the SpecialConstraint
void NodalConstraintDOF::WriteMonitor(FILE *f, bool first_record, double time)
{

}

//Checks which Lagrange multipliers will be ativados,according to the activation of the GLs of the nodes of the which the special constraint participates
void NodalConstraintDOF::ActivateDOFs()
{
	//Active GLs of the nodes A and B and the Lagrange multipliers
	if (UX_table.GetAt(db.current_solution_number - 1) == 1) {
		db.nodes[node_A - 1]->active_GL[0] = 1;
		db.nodes[node_B - 1]->active_GL[0] = 1;
		active_lambda[0] = 1;
	}
	else
	{
		active_lambda[0] = 0;
	}

	if (UY_table.GetAt(db.current_solution_number - 1) == 1) {
		db.nodes[node_A - 1]->active_GL[1] = 1;
		db.nodes[node_B - 1]->active_GL[1] = 1;
		active_lambda[1] = 1;
	}
	else
	{
		active_lambda[1] = 0;
	}

	if (UZ_table.GetAt(db.current_solution_number - 1) == 1) {
		db.nodes[node_A - 1]->active_GL[2] = 1;
		db.nodes[node_B - 1]->active_GL[2] = 1;
		active_lambda[2] = 1;
	}
	else
	{
		active_lambda[2] = 0;
	}

	if (ROTX_table.GetAt(db.current_solution_number - 1) == 1) {
		db.nodes[node_A - 1]->active_GL[3] = 1;
		db.nodes[node_B - 1]->active_GL[3] = 1;
		active_lambda[3] = 1;
	}
	else
	{
		active_lambda[3] = 0;
	}

	if (ROTY_table.GetAt(db.current_solution_number - 1) == 1) {
		db.nodes[node_A - 1]->active_GL[4] = 1;
		db.nodes[node_B - 1]->active_GL[4] = 1;
		active_lambda[4] = 1;
	}
	else
	{
		active_lambda[4] = 0;
	}

	if (ROTZ_table.GetAt(db.current_solution_number - 1) == 1) {
		db.nodes[node_A - 1]->active_GL[5] = 1;
		db.nodes[node_B - 1]->active_GL[5] = 1;
		active_lambda[5] = 1;
	}
	else
	{
		active_lambda[5] = 0;
	}
}



bool NodalConstraintDOF::Check()
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

void NodalConstraintDOF::ComputeVelAccel()
{
	//If for the step of creation of the constraint initializes initial conditions
	if (UX_table.GetAt(db.current_solution_number - 1))
	{
		//Iterates over GL and sets initial conditions
		if (active_lambda[0] == 1)
		{
			db.nodes[node_B - 1]->vel[0] = db.nodes[node_A - 1]->vel[0];
			db.nodes[node_B - 1]->accel[0] = db.nodes[node_A - 1]->accel[0];
		}
	}
	if (UY_table.GetAt(db.current_solution_number - 1))
	{
		//Iterates over GL and sets initial conditions
		if (active_lambda[1] == 1)
		{
			db.nodes[node_B - 1]->vel[1] = db.nodes[node_A - 1]->vel[1];
			db.nodes[node_B - 1]->accel[1] = db.nodes[node_A - 1]->accel[1];
		}
	}
	if (UZ_table.GetAt(db.current_solution_number - 1))
	{
		//Iterates over GL and sets initial conditions
		if (active_lambda[2] == 1)
		{
			db.nodes[node_B - 1]->vel[2] = db.nodes[node_A - 1]->vel[2];
			db.nodes[node_B - 1]->accel[2] = db.nodes[node_A - 1]->accel[2];
		}
	}
	if (ROTX_table.GetAt(db.current_solution_number - 1))
	{
		//Iterates over GL and sets initial conditions
		if (active_lambda[3] == 1)
		{
			db.nodes[node_B - 1]->vel[3] = db.nodes[node_A - 1]->vel[3];
			db.nodes[node_B - 1]->accel[3] = db.nodes[node_A - 1]->accel[3];
		}
	}
	if (ROTY_table.GetAt(db.current_solution_number - 1))
	{
		//Iterates over GL and sets initial conditions
		if (active_lambda[4] == 1)
		{
			db.nodes[node_B - 1]->vel[4] = db.nodes[node_A - 1]->vel[4];
			db.nodes[node_B - 1]->accel[4] = db.nodes[node_A - 1]->accel[4];
		}
	}
	if (ROTZ_table.GetAt(db.current_solution_number - 1))
	{
		//Iterates over GL and sets initial conditions
		if (active_lambda[5] == 1)
		{
			db.nodes[node_B - 1]->vel[5] = db.nodes[node_A - 1]->vel[5];
			db.nodes[node_B - 1]->accel[5] = db.nodes[node_A - 1]->accel[5];
		}
	}
}

void NodalConstraintDOF::ComputeInitialGuessDisplacements()
{
	//If for the step of creation of the constraint initializes initial conditions
	if (UX_table.GetAt(db.current_solution_number - 1) == true)
	{
		//Iterates over GL and sets initial conditions
		if (active_lambda[0] == 1)
		{
			db.nodes[node_B - 1]->displacements[0] = db.nodes[node_A - 1]->displacements[0];
		}
	}
	if (UY_table.GetAt(db.current_solution_number - 1) == true)
	{
		//Iterates over GL and sets initial conditions
		if (active_lambda[1] == 1)
		{
			db.nodes[node_B - 1]->displacements[1] = db.nodes[node_A - 1]->displacements[1];
		}
	}
	if (UZ_table.GetAt(db.current_solution_number - 1) == true)
	{
		//Iterates over GL and sets initial conditions
		if (active_lambda[2] == 1)
		{
			db.nodes[node_B - 1]->displacements[2] = db.nodes[node_A - 1]->displacements[2];
		}
	}
	if (ROTX_table.GetAt(db.current_solution_number - 1) == true)
	{
		//Iterates over GL and sets initial conditions
		if (active_lambda[3] == 1)
		{
			db.nodes[node_B - 1]->displacements[3] = db.nodes[node_A - 1]->displacements[3];
		}
	}
	if (ROTY_table.GetAt(db.current_solution_number - 1) == true)
	{
		//Iterates over GL and sets initial conditions
		if (active_lambda[4] == 1)
		{
			db.nodes[node_B - 1]->displacements[4] = db.nodes[node_A - 1]->displacements[4];
		}
	}
	if (ROTZ_table.GetAt(db.current_solution_number - 1) == true)
	{
		//Iterates over GL and sets initial conditions
		if (active_lambda[5] == 1)
		{
			db.nodes[node_B - 1]->displacements[5] = db.nodes[node_A - 1]->displacements[5];
		}
	}
}

void NodalConstraintDOF::WriteVTK_XMLRender(FILE* f)
{
}

//Saving variables of the configuration converged
void NodalConstraintDOF::SaveLagrange()
{
	//PrintPtr(lambda, 3);
	for (int i = 0; i < n_GL; i++)
		copy_lambda[i] = lambda[i];
}

//Assembly of the residuals and stiffness tangent
void NodalConstraintDOF::Mount()
{
	//In this constraint the stiffness tangent not if modifies never. and montada in the PreCalc()
	//If the constraint is active, performs the assembly
	for (int i = 0; i < n_GL; i++) {
		if (active_lambda[i] == 1)
		{
			r(i, 0) = db.nodes[node_A - 1]->displacements[i] - db.nodes[node_B - 1]->displacements[i];
			//Update of the vector of residuals
			(*residual)(i, 0) = lambda[i];
			(*residual)(6 + i, 0) = -lambda[i];
			(*residual)(12 + i, 0) = r(i, 0);
		}
	}
}

//Fills the contribution of the element in the matrices global
void NodalConstraintDOF::MountGlobal()
{
	//Temporary variables for save the indexing global of the degrees of freedom the to be set in the stiffness matrix global
	int GL_global_1 = 0;
	int GL_global_2 = 0;
	double anterior = 0;
	//If the constraint is active, performs the assembly
	if (active_lambda[0] == 1 || active_lambda[1] == 1 || active_lambda[2] == 1 || active_lambda[3] == 1 || active_lambda[4] == 1 || active_lambda[5] == 1)
	{
		for (int i = 0; i < 18; i++)
		{
			if (active_lambda[i % 6] == 1)
			{
				//////////////ASSEMBLY OF THE VECTOR OF INTERNAL FORCES UNBALANCED//////////////////
				//Takes of the global-DOF vector, the indexing of each global degree of freedom
				if (i < 6)//uA
				{
					GL_global_1 = db.nodes[node_A - 1]->GLs[i];
				}
				else
				{
					if (i < 12)//uB
					{
						GL_global_1 = db.nodes[node_B - 1]->GLs[i - 6];
					}
					else
					{
						//lambda
						GL_global_1 = GLs[i - 12];
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

				for (int j = 0; j < 18; j++)
				{
					if (active_lambda[j % 6] == 1)
					{
						//////////////////////ASSEMBLY OF THE STIFFNESS MATRIX/////////////////////////
						//Takes of the global-DOF vector, the indexing of each global degree of freedom
						if (j < 6)//uA
						{
							GL_global_2 = db.nodes[node_A - 1]->GLs[j];
						}
						else
						{
							if (j < 12)//uB
							{
								GL_global_2 = db.nodes[node_B - 1]->GLs[j - 6];
							}
							else
							{
								//lambda
								GL_global_2 = GLs[j - 12];
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
	}
}



//Pre-calculation of variables that and done a single time in the start
void NodalConstraintDOF::PreCalc()
{
	for (int i = 0; i < 6; i++)
	{
		for (int j = 0; j < 6; j++)
		{
			(*stiffness)(i, j + 12) = I6(i, j);
			(*stiffness)(i + 6, j + 12) = -I6(i, j);
			(*stiffness)(i + 12, j) = I6(i, j);
			(*stiffness)(i + 12, j + 6) = -I6(i, j);
		}
	}
}
