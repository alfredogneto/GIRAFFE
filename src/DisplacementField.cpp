#include "DisplacementField.h"

#include "Node.h"
#include "CoordinateSystem.h"
#include "Solution.h"
#include "IO.h"

DisplacementField::DisplacementField()
{
	number = 0;
	cs = 0;
	n_nodes = 0;
	nodes = NULL;
	displacements = NULL;
	n_nodes = 0;;
	n_values = 6;
	alloced = false;

	Q = Matrix(3, 3);
	I = Matrix(3, 3);
	I(0, 0) = 1.0;
	I(1, 1) = 1.0;
	I(2, 2) = 1.0;

	sprintf(angular_parameters, "EulerVector");
}

DisplacementField::~DisplacementField()
{
	Free();
}

bool  DisplacementField::Read(FILE *f)
{
	char s[1000];
	fscanf(f, "%s", s);
	number = atoi(s);
	fscanf(f, "%s", s);
	if (!strcmp(s, "NNodes"))
	{
		fscanf(f, "%s", s);
		n_nodes = atoi(s);
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "CS"))
	{
		fscanf(f, "%s", s);
		cs = atoi(s);
	}
	else
		return false;

	//Salva a posição (stream)
	fpos_t pos;
	fgetpos(f, &pos);
	fscanf(f, "%s", s);
	if (!strcmp(s, "AngularParameters"))
	{
		fscanf(f, "%s", s);
		sprintf(angular_parameters, s);
	}
	else
	{
		fsetpos(f, &pos);
		sprintf(angular_parameters, "EulerAngles");
	}

	fscanf(f, "%s", s);
	if (!strcmp(s, "SolutionStep"))
	{
		fscanf(f, "%s", s);
		solution_step = atoi(s);
	}
	else
		return false;

	


	TryComment(f);
	Alloc();
	//Reads displacement data
	for (int i = 0; i < n_nodes; i++)
	{
		//Node
		fscanf(f, "%s", s);
		nodes[i] = atoi(s);
		//Displacement
		for (int j = 0; j < n_values; j++)
		{
			fscanf(f, "%s", s);
			displacements[i][j] = atof(s);
		}
	}
	return true;
}

//Alocação dinâmica de variáveis
void DisplacementField::Alloc()
{
	Free();
	//Alocando
	nodes = new int[n_nodes];
	displacements = new double*[n_nodes];
	for (int i = 0; i < n_nodes; i++)
		displacements[i] = new double[n_values];
	//Zerando valores
	for (int i = 0; i < n_nodes; i++)
		for (int j = 0; j < n_values; j++)
			displacements[i][j] = 0.0;
	alloced = true;
}

//Desalocação dinâmica de variáveis
void DisplacementField::Free()
{
	if (alloced == true)
	{
		delete[]nodes;
		for (int i = 0; i < n_nodes; i++)
			delete[]displacements[i];
		delete[]displacements;

		nodes = NULL;
		displacements = NULL;
	}
}

//Checking inconsistencies
bool DisplacementField::Check()
{
	if (cs > db.number_CS)
		return false;
	for (int i = 0; i < n_nodes; i++)
	{
		if (nodes[i] > db.number_nodes)
			return false;
	}
	if (solution_step > db.number_solutions)
		return false;
	return true;
}

void  DisplacementField::Write(FILE *f)
{
	fprintf(f, "DisplacementField\t%d\tNNodes\t%d\tCS\t%d\tAngularParameters\t%s\tSolutionStep\t%d\n", number, n_nodes, cs, angular_parameters,solution_step);
	//Write displacement data
	for (int i = 0; i < n_nodes; i++)
	{
		fprintf(f, "%d\t", nodes[i]);
		for (int j = 0; j < n_values; j++)
			fprintf(f, "%.6e\t", displacements[i][j]);
		fprintf(f, "\n");
	}
}

//Writes VTK XML data for post-processing
void DisplacementField::WriteVTK_XML(FILE *f)
{

}

//Pre-calculus
void DisplacementField::PreCalc()
{
	
}

void DisplacementField::Mount()
{
	if (db.current_solution_number == solution_step)
	{
		int temp_node = 0;
		int GL;
		double factor = (db.last_converged_time + db.current_time_step - db.solution[solution_step - 1]->start_time) / (db.solution[solution_step - 1]->end_time - db.solution[solution_step - 1]->start_time);
		//Varredura no node set pelos nós definidos, para imposição dos deslocamentos
		for (int index_node = 0; index_node < n_nodes; index_node++)
		{
			temp_node = nodes[index_node];
			//Translações
			Matrix disp(3, 1);
			for (int j = 0; j < 3; j++)
				disp(j, 0) = displacements[index_node][j] * factor;
			//Coordinate transformation (to global)
			disp = transp(*db.CS[cs - 1]->Q)*disp;
			for (int j = 0; j < 3; j++)
			{
				if (db.nodes[temp_node - 1]->constraints[j] == 1 && db.nodes[temp_node - 1]->active_GL[j] == 1)
				{
					GL = -db.nodes[temp_node - 1]->GLs[j] - 1;
					db.global_X_B(GL, 0) = disp(j, 0);
					db.nodes[temp_node - 1]->displacements[j] = db.global_X_B(GL, 0);
					db.nodes[temp_node - 1]->vel[j] = db.global_X_B(GL, 0) / db.current_time_step;	//Velocidade prescrita
					db.nodes[temp_node - 1]->accel[j] = (db.nodes[temp_node - 1]->vel[j] - db.nodes[temp_node - 1]->copy_vel[j]) / db.current_time_step;	//Aceleração prescrita
				}
			}
			
			//Rotações
			Matrix euler(3, 1);
			Matrix rodrigues(3, 1);
			MountEulerVector(&euler, &rodrigues, temp_node, index_node, factor);

			for (int j = 3; j < 6; j++)
			{
				//Se o grau de liberdade for fixo
				if (db.nodes[temp_node - 1]->constraints[j] == 1 && db.nodes[temp_node - 1]->active_GL[j] == 1)
				{
					GL = -db.nodes[temp_node - 1]->GLs[j] - 1;
					db.global_X_B(GL, 0) = rodrigues(j - 3, 0);
					db.nodes[temp_node - 1]->displacements[j] = db.global_X_B(GL, 0);
					db.nodes[temp_node - 1]->vel[j] = euler(j - 3, 0) / db.current_time_step;	//Velocidade prescrita
					db.nodes[temp_node - 1]->accel[j] = (db.nodes[temp_node - 1]->vel[j] - db.nodes[temp_node - 1]->copy_vel[j]) / db.current_time_step;	//Aceleração prescrita
				}
			}
		}
	}
}


void DisplacementField::MountEulerVector(Matrix* euler, Matrix* rodrigues, int temp_node, int index_node, double factor)
{
	//Rotações
	if (!strcmp(angular_parameters, "EulerVector"))
	{
		Matrix euler_ip(3, 1);
		if (db.nodes[temp_node - 1]->constraints[3] == 1)
			euler_ip(0, 0) = displacements[index_node][3] * factor;
		if (db.nodes[temp_node - 1]->constraints[4] == 1)
			euler_ip(1, 0) = displacements[index_node][4] * factor;
		if (db.nodes[temp_node - 1]->constraints[5] == 1)
			euler_ip(2, 0) = displacements[index_node][5] * factor;
		Matrix direction(3, 1);
		double euler_escalar = norm(euler_ip);
		if (euler_escalar != 0.0)
			direction = (1.0 / euler_escalar)*euler_ip;
		double rod_escalar = 2.0*tan(euler_escalar / 2.0);
		Matrix rodrigues_ip(3, 1);
		rodrigues_ip = rod_escalar * direction;

		////////////////////////////////////////////////////////////////////////////////////////
		Matrix A = skew(rodrigues_ip);
		double g = 4.0 / (4.0 + rod_escalar * rod_escalar);
		Q = I + g * (A + 0.5*A*A);
		////////////////////////////////////////////////////////////////////////////////////////
	}
	//Rotações
	if (!strcmp(angular_parameters, "EulerAngles"))
	{
		double phi = 0.0;
		double theta = 0.0;
		double psi = 0.0;
		if (db.nodes[temp_node - 1]->constraints[3] == 1)
			phi = displacements[index_node][3] * factor;
		if (db.nodes[temp_node - 1]->constraints[4] == 1)
			theta = displacements[index_node][4] * factor;
		if (db.nodes[temp_node - 1]->constraints[5] == 1)
			psi = displacements[index_node][5] * factor;

		Q(0, 0) = cos(theta)*cos(psi);
		Q(0, 1) = sin(phi)*sin(theta)*cos(psi) - cos(phi)*sin(psi);
		Q(0, 2) = cos(phi)*sin(theta)*cos(psi) + sin(phi)*sin(psi);

		Q(1, 0) = cos(theta)*sin(psi);
		Q(1, 1) = sin(phi)*sin(theta)*sin(psi) + cos(phi)*cos(psi);
		Q(1, 2) = cos(phi)*sin(theta)*sin(psi) - sin(phi)*cos(psi);

		Q(2, 0) = -sin(theta);
		Q(2, 1) = sin(phi)*cos(theta);
		Q(2, 2) = cos(phi)*cos(theta);
	}
	//Coordinate transformation (to global)
	Q = transp(*db.CS[cs - 1]->Q) * Q * (*db.CS[cs - 1]->Q);
	//Previous rotation vector (from node)
	Matrix Qi = *(db.nodes[temp_node - 1]->Q);
	Matrix Qd = Q * transp(Qi);
	double trQd = Qd(0, 0) + Qd(1, 1) + Qd(2, 2);
	Matrix skewQd = 0.5 * (Qd - transp(Qd));

	//Extraction of euler
	double theta = 2.0 * asin(sqrt((3.0 - trQd) / 4.0));
	if (theta != 0.0)
		(*euler) = (theta / sin(theta)) * axial(skewQd);
	else
	{
		(*euler)(0, 0) = 0.0;
		(*euler)(1, 0) = 0.0;
		(*euler)(2, 0) = 0.0;
	}

	//Extraction of rodrigues
	double alpha = 2.0 * sqrt((3.0 - trQd) / (1.0 + trQd));
	(*rodrigues) = ((4.0 + alpha * alpha) / 4.0) * axial(skewQd);
}

void DisplacementField::EvaluateExplicit(double t)
{

}