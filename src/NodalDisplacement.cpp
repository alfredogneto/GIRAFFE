#include "NodalDisplacement.h"
#include "IO.h"

NodalDisplacement::NodalDisplacement()
{
	number = 0;
	table = NULL;
	mcode = NULL;
	n_times = 0;
	n_values = 6;

	sprintf(angular_parameters, "EulerAngles");

	//Single node
	single_node = 0;

	Q = Matrix(3, 3);
	I = Matrix(3, 3);
	I(0, 0) = 1.0;
	I(1, 1) = 1.0;
	I(2, 2) = 1.0;
}

NodalDisplacement::~NodalDisplacement()
{
	if (table != NULL)
		delete table;
	if (mcode != NULL)
		delete mcode;
}

bool  NodalDisplacement::Read(FILE *f)
{
	char s[1000];
	fscanf(f, "%s", s);
	number = atoi(s);
	fscanf(f, "%s", s);
	if (!strcmp(s, "NodeSet"))
	{
		fscanf(f, "%s", s);
		node_set = atoi(s);
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
	//////////////////////////////////////////
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
		sprintf(angular_parameters, "EulerAngles");//default
	}
	///////////////////////////////////////////
	//Salva a posição (stream)
	fgetpos(f, &pos);
	fscanf(f, "%s", s);
	if (!strcmp(s, "BoolTable"))
		bool_table.Read(f);
	else
	{
		fsetpos(f, &pos);
		bool_table.SetDefault(true);
	}
	///////////////////////////////////////////
	fscanf(f, "%s", s);
	if (!strcmp(s, "NTimes"))
	{
		fscanf(f, "%s", s);
		n_times = atoi(s);
		//Allocating table
		table = new Table(n_times, n_values);
		TryComment(f);
		//Reads table data
		if (!table->Read(f))
			return false;
			
	}
	else
	{
		if (!strcmp(s, "MathCode"))
		{
			//Allocating MathCode
			mcode = new MathCode(n_values);
			//Reads MathCode
			if (!mcode->Read(f))
				return false;
		}
		else
			return false;
	}
	
	return true;
}

//Checking inconsistencies
bool NodalDisplacement::Check()
{
	if (cs > db.number_CS)
		return false;
	if (node_set > db.number_node_sets)
		return false;
	return true;
}

void  NodalDisplacement::Write(FILE *f)
{
	fprintf(f, "NodalDisplacement\t%d\tNodeSet\t%d\tCS\t%d\tAngularParameters\t%s\t", number, node_set, cs, angular_parameters);
	bool_table.Write(f);
	if (table != NULL)
	{
		fprintf(f, "NTimes\t%d\n", n_times);
		table->Write(f);
	}
		
	if (mcode != NULL)
	{
		fprintf(f, "MathCode\n");
		mcode->Write(f);
	}
		
}

//Writes VTK XML data for post-processing
void NodalDisplacement::WriteVTK_XML(FILE *f)
{
	
}

//Pre-calculus
void NodalDisplacement::PreCalc()
{

}

void NodalDisplacement::Mount()
{
	if (bool_table.GetAt(db.current_solution_number - 1))
	{
		int temp_node = 0;
		int GL;
		double time_value = db.last_converged_time + db.current_time_step;
		double dt = db.current_time_step;
		//Varredura no node set pelos nós definidos, para imposição dos deslocamentos
		for (int index_node = 0; index_node < db.node_sets[node_set - 1]->n_nodes; index_node++)
		{
			temp_node = db.node_sets[node_set - 1]->node_list[index_node];
			//Translações
			Matrix disp(3, 1);
			for (int j = 0; j < 3; j++)
				disp(j, 0) = GetValueAt(time_value, j) - GetValueAt(time_value - dt, j);

			//Coordinate transformation (to global)
			disp = transp(*db.CS[cs - 1]->Q)*disp;
			for (int j = 0; j < 3; j++)
			{
				if (db.nodes[temp_node - 1]->constraints[j] == 1 && db.nodes[temp_node - 1]->active_GL[j] == 1)
				{
					GL = -db.nodes[temp_node - 1]->GLs[j] - 1;
					db.global_X_B(GL, 0) = disp(j, 0);
					db.nodes[temp_node - 1]->displacements[j] = db.global_X_B(GL, 0);
					db.nodes[temp_node - 1]->vel[j] = db.global_X_B(GL, 0) / dt;	//Velocidade prescrita
					db.nodes[temp_node - 1]->accel[j] = (db.nodes[temp_node - 1]->vel[j] - db.nodes[temp_node - 1]->copy_vel[j]) / dt;	//Aceleração prescrita
				}
			}

			//Rotações
			Matrix euler(3, 1);
			Matrix rodrigues(3, 1);
			MountEulerVector(&euler, &rodrigues, temp_node, time_value);

			for (int j = 3; j < 6; j++)
			{
				//Se o grau de liberdade for fixo
				if (db.nodes[temp_node - 1]->constraints[j] == 1 && db.nodes[temp_node - 1]->active_GL[j] == 1)
				{
					GL = -db.nodes[temp_node - 1]->GLs[j] - 1;
					db.global_X_B(GL, 0) = rodrigues(j - 3, 0);
					db.nodes[temp_node - 1]->displacements[j] = db.global_X_B(GL, 0);
					db.nodes[temp_node - 1]->vel[j] = euler(j - 3, 0) / dt;	//Velocidade prescrita
					db.nodes[temp_node - 1]->accel[j] = (db.nodes[temp_node - 1]->vel[j] - db.nodes[temp_node - 1]->copy_vel[j]) / dt;	//Aceleração prescrita
				}
			}
		}
	}
}

void NodalDisplacement::EvaluateExplicit(double t)
{
	if (bool_table.GetAt(db.current_solution_number - 1))
	{
		int temp_node = 0;
		
		
		double t0 = db.last_converged_time;
		if (t == t0)
			t0 = db.last_converged_time - db.current_time_step;
		//Varredura no node set pelos nós definidos, para imposição dos deslocamentos
		for (int index_node = 0; index_node < db.node_sets[node_set - 1]->n_nodes; index_node++)
		{
			temp_node = db.node_sets[node_set - 1]->node_list[index_node];
			//Translações
			Matrix disp(3, 1);
			for (int j = 0; j < 3; j++)
				disp(j, 0) = GetValueAt(t, j) - GetValueAt(t0, j);

			//Coordinate transformation 
			disp = transp( (*db.CS[cs - 1]->Q)) *  disp;
			for (int j = 0; j < 3; j++)
			{
				if (db.nodes[temp_node - 1]->constraints[j] == 1 && db.nodes[temp_node - 1]->active_GL[j] == 1)
				{
					//Variables to be used during time-integration
					(*db.nodes[temp_node - 1]->u)(j, 0) = disp(j, 0);
					double v0 = db.nodes[temp_node - 1]->copy_vel[j];
					(*db.nodes[temp_node - 1]->du)(j, 0) = (*db.nodes[temp_node - 1]->u)(j, 0) / (t - t0);
					(*db.nodes[temp_node - 1]->ddu)(j, 0) = ((*db.nodes[temp_node - 1]->du)(j, 0) - v0) / (t - t0);
					
					//Variables to be used after time-integration
					db.nodes[temp_node - 1]->displacements[j] = (*db.nodes[temp_node - 1]->u)(j, 0);
					db.nodes[temp_node - 1]->vel[j] = (*db.nodes[temp_node - 1]->du)(j, 0);
					db.nodes[temp_node - 1]->accel[j] = (*db.nodes[temp_node - 1]->ddu)(j, 0);
				}
					
			}

			//Rotações
			Matrix euler(3, 1);
			Matrix rodrigues(3, 1);
			MountEulerVector(&euler, &rodrigues, temp_node, t);

			for (int j = 0; j < 3; j++)
			{
				//Se o grau de liberdade for fixo
				if (db.nodes[temp_node - 1]->constraints[j + 3] == 1 && db.nodes[temp_node - 1]->active_GL[j + 3] == 1)
				{
					(*db.nodes[temp_node - 1]->alpha)(j, 0) = rodrigues(j, 0);
					double omega0 = db.nodes[temp_node - 1]->copy_vel[j + 3];
					(*db.nodes[temp_node - 1]->omega)(j, 0) = (*db.nodes[temp_node - 1]->alpha)(j, 0) / (t - t0);
					(*db.nodes[temp_node - 1]->domega)(j, 0) = ((*db.nodes[temp_node - 1]->omega)(j, 0) - omega0) / (t - t0);

					//Variables to be used after time-integration
					db.nodes[temp_node - 1]->displacements[j + 3] = (*db.nodes[temp_node - 1]->alpha)(j, 0);
					db.nodes[temp_node - 1]->vel[j + 3] = (*db.nodes[temp_node - 1]->omega)(j, 0);
					db.nodes[temp_node - 1]->accel[j + 3] = (*db.nodes[temp_node - 1]->domega)(j, 0);
				}
					
			}
		}
	}
}

void NodalDisplacement::MountSingleNodeDisplacement()
{
	if (bool_table.GetAt(db.current_solution_number - 1))
	{
		int temp_node = single_node;
		int GL;
		double time_value = db.last_converged_time + db.current_time_step;
		double dt = db.current_time_step;

		//Translações
		Matrix disp(3, 1);
		for (int j = 0; j < 3; j++)
			disp(j, 0) = GetValueAt(time_value, j) - GetValueAt(time_value - dt, j);
		for (int j = 0; j < 3; j++)
		{
			if (db.nodes[temp_node - 1]->constraints[j] == 1 && db.nodes[temp_node - 1]->active_GL[j] == 1)
			{
				GL = -db.nodes[temp_node - 1]->GLs[j] - 1;
				db.global_X_B(GL, 0) = disp(j, 0);
				db.nodes[temp_node - 1]->displacements[j] = db.global_X_B(GL, 0);
				db.nodes[temp_node - 1]->vel[j] = db.global_X_B(GL, 0) / dt;	//Velocidade prescrita
				db.nodes[temp_node - 1]->accel[j] = (db.nodes[temp_node - 1]->vel[j] - db.nodes[temp_node - 1]->copy_vel[j]) / dt;	//Aceleração prescrita
			}
		}
		//Rotações
		Matrix euler(3, 1);
		Matrix rodrigues(3, 1);
		MountEulerVector(&euler, &rodrigues, temp_node, time_value);

		for (int j = 3; j < 6; j++)
		{
			//Se o grau de liberdade for fixo
			if (db.nodes[temp_node - 1]->constraints[j] == 1 && db.nodes[temp_node - 1]->active_GL[j] == 1)
			{
				GL = -db.nodes[temp_node - 1]->GLs[j] - 1;
				db.global_X_B(GL, 0) = rodrigues(j - 3, 0);
				db.nodes[temp_node - 1]->displacements[j] = db.global_X_B(GL, 0);
				db.nodes[temp_node - 1]->vel[j] = euler(j - 3, 0) / dt;	//Velocidade prescrita
				db.nodes[temp_node - 1]->accel[j] = (db.nodes[temp_node - 1]->vel[j] - db.nodes[temp_node - 1]->copy_vel[j]) / dt;	//Aceleração prescrita
			}
		}
	}
}
double NodalDisplacement::GetValueAt(double t, int position)
{
	if (table != NULL)
		return table->GetValueAt(t, position);
	if (mcode != NULL)
		return mcode->GetValueAt(t, position);
	return 0.0;
}

void NodalDisplacement::MountEulerVector(Matrix* euler, Matrix* rodrigues, int temp_node, double time)
{
	//Rotações
	if (!strcmp(angular_parameters, "EulerVector"))
	{
		Matrix euler_ip(3, 1);
		if (db.nodes[temp_node - 1]->constraints[3] == 1)
			euler_ip(0, 0) = GetValueAt(time, 3);
		if (db.nodes[temp_node - 1]->constraints[4] == 1)
			euler_ip(1, 0) = GetValueAt(time, 4);
		if (db.nodes[temp_node - 1]->constraints[5] == 1)
			euler_ip(2, 0) = GetValueAt(time, 5);
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
			phi = GetValueAt(time, 3);
		if (db.nodes[temp_node - 1]->constraints[4] == 1)
			theta = GetValueAt(time, 4);
		if (db.nodes[temp_node - 1]->constraints[5] == 1)
			psi = GetValueAt(time, 5);

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
	//Q.print();
	//Coordinate transformation (to global)
	if (single_node == 0)
		Q = transp(*db.CS[cs - 1]->Q) * Q * (*db.CS[cs - 1]->Q);
	//Q.print();
	//Previous rotation vector (from node)
	Matrix Qi = *(db.nodes[temp_node - 1]->Q);
	//Material or global descriptions
	Matrix Qd;
	if (db.nodes[temp_node-1]->flag_material_description == false)
		Qd = Q * transp(Qi);
	else
		Qd = transp(Qi) * Q;
	double trQd = Qd(0, 0) + Qd(1, 1) + Qd(2, 2);
	if (trQd > 3.0)
		trQd = 3.0;
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