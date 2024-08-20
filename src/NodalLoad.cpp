#include "NodalLoad.h"
#include "IO.h"

NodalLoad::NodalLoad()
{
	number = 0;
	table = NULL;
	mcode = NULL;
	n_times = 0;
	n_values = 6;

	n_nodes_f = new int[3];			//número de nós para divisão de forças - 3 componentes
	n_nodes_m = new int[3];			//número de nós para divisão de momentos - 3 componentes
	mult_f = new double[3];			//multiplicador para os esforços de força
	mult_m = new double[3];			//multiplicador para os esforços de momento
}

NodalLoad::~NodalLoad()
{
	if (table != NULL)
		delete table;
	if (mcode != NULL)
		delete mcode;

	delete[]n_nodes_f;
	delete[]n_nodes_m;
	delete[]mult_f;
	delete[]mult_m;
}
//Reads input file
bool NodalLoad::Read(FILE *f)
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
bool NodalLoad::Check()
{
	if (cs > db.number_CS)
		return false;
	if (node_set > db.number_node_sets)
		return false;
	return true;
}

//Writes output file
void NodalLoad::Write(FILE *f)
{
	fprintf(f, "NodalLoad\t%d\tNodeSet\t%d\tCS\t%d\t", number, node_set, cs);
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
void NodalLoad::WriteVTK_XML(FILE *f)
{
	//vetores para escrita no formato binário - usando a função 'enconde'
	std::vector<float> float_vector;
	std::vector<int> int_vector;
	int count = db.node_sets[node_set - 1]->n_nodes;
	//Opens Piece
	fprintf(f, "\t\t<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", count, count);
	//Opens Points
	fprintf(f, "\t\t\t<Points>\n");
	//Opens DataArray
	fprintf(f, "\t\t\t\t<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"binary\">\n");
	float_vector.clear();
	//Preenchendo as coordenadas dos pontos
	int node;
	for (int index = 0; index < db.node_sets[node_set - 1]->n_nodes; index++)
	{
		//Node number
		node = db.node_sets[node_set - 1]->node_list[index];
		float_vector.push_back((float)(db.nodes[node - 1]->copy_coordinates[0]));
		float_vector.push_back((float)(db.nodes[node - 1]->copy_coordinates[1]));
		float_vector.push_back((float)(db.nodes[node - 1]->copy_coordinates[2]));
	}
	fprintf(f, encodeData<float>(float_vector).c_str());
	fprintf(f, "\n");
	//Closes DataArray
	fprintf(f, "\t\t\t\t</DataArray>\n");
	//Closes Points
	fprintf(f, "\t\t\t</Points>\n");
	//Opens Cells
	fprintf(f, "\t\t\t<Cells>\n");
	//Opens DataArray
	fprintf(f, "\t\t\t\t<DataArray type=\"Int32\" Name=\"connectivity\" format=\"binary\">\n");
	int_vector.clear();
	for (int cell = 0; cell < count; cell++)
		int_vector.push_back(cell);
	fprintf(f, encodeData(int_vector).c_str());
	fprintf(f, "\n");
	//Closes DataArray
	fprintf(f, "\t\t\t\t</DataArray>\n");
	//Opens DataArray
	int_vector.clear();
	fprintf(f, "\t\t\t\t<DataArray type=\"Int32\" Name=\"types\" format=\"binary\">\n");
	for (int cell = 0; cell < count; cell++)
		int_vector.push_back(1);
	fprintf(f, encodeData(int_vector).c_str());
	fprintf(f, "\n");
	//Closes DataArray
	fprintf(f, "\t\t\t\t</DataArray>\n");
	//Opens DataArray
	fprintf(f, "\t\t\t\t<DataArray type=\"Int32\" Name=\"offsets\" format=\"binary\">\n");
	int_vector.clear();
	for (int cell = 0; cell < count; cell++)
		int_vector.push_back(cell + 1);
	fprintf(f, encodeData(int_vector).c_str());
	fprintf(f, "\n");
	//Closes DataArray
	fprintf(f, "\t\t\t\t</DataArray>\n");
	//Closes Cells
	fprintf(f, "\t\t\t</Cells>\n");

	//Opens PointData
	fprintf(f, "\t\t\t<PointData Vectors = \"Loads\">\n");

	float_vector.clear();
	//Opens DataArray
	fprintf(f, "\t\t\t\t<DataArray Name = \"Force\" type = \"Float32\" NumberOfComponents= \"3\" format=\"binary\">\n");
	Matrix force(3);
	for (int index = 0; index < db.node_sets[node_set - 1]->n_nodes; index++)
	{
		//Node number
		node = db.node_sets[node_set - 1]->node_list[index];
		//Evaluating input data at current time
		for (int i = 0; i < 3; i++)
			force(i, 0) = mult_f[i] * GetValueAt(db.last_converged_time + db.current_time_step, i);
		//Coordinate transformation (to global)
		force = transp(*db.CS[cs - 1]->Q)*force;
		float_vector.push_back((float)(force(0, 0)));
		float_vector.push_back((float)(force(1, 0)));
		float_vector.push_back((float)(force(2, 0)));
	}
	fprintf(f, encodeData(float_vector).c_str());
	fprintf(f, "\n");
	//Closes DataArray
	fprintf(f, "\t\t\t\t</DataArray>\n");
	float_vector.clear();
	//Opens DataArray
	fprintf(f, "\t\t\t\t<DataArray Name = \"Moment\" type = \"Float32\" NumberOfComponents= \"3\" format=\"binary\">\n");
	Matrix moment(3);
	for (int index = 0; index < db.node_sets[node_set - 1]->n_nodes; index++)
	{
		//Node number
		node = db.node_sets[node_set - 1]->node_list[index];
		//Evaluating input data at current time
		for (int i = 0; i < 3; i++)
			moment(i, 0) = mult_m[i] * GetValueAt(db.last_converged_time + db.current_time_step, i + 3);
		//Coordinate transformation (to global)
		moment = transp(*db.CS[cs - 1]->Q)*moment;
		float_vector.push_back((float)(moment(0, 0)));
		float_vector.push_back((float)(moment(1, 0)));
		float_vector.push_back((float)(moment(2, 0)));
	}
	fprintf(f, encodeData(float_vector).c_str());
	fprintf(f, "\n");
	//Closes DataArray
	fprintf(f, "\t\t\t\t</DataArray>\n");

	//Closes PointData
	fprintf(f, "\t\t\t</PointData>\n");

	//Closes Piece
	fprintf(f, "\t\t</Piece>\n");
}
//Pre-calculus
void NodalLoad::PreCalc()
{
	
}

//Atualiza dados necessários e que sejam dependentes de DOFs ativos/inativos - chamado no início de cada solution step
void NodalLoad::UpdateforSolutionStep()
{
	int node;
	int n_nodes = db.node_sets[node_set - 1]->n_nodes;
	for (int j = 0; j < 3; j++)	//zerando os contadores
	{
		n_nodes_f[j] = 0;
		n_nodes_m[j] = 0;
		mult_f[j] = 0.0;
		mult_m[j] = 0.0;
	}
	//Percorre os nós do NodeSet
	for (int i = 0; i < n_nodes; i++)
	{
		node = db.node_sets[node_set - 1]->node_list[i];
		for (int j = 0; j < 3; j++)
		{
			if (db.nodes[node - 1]->active_GL[j] != 0)	//se o GL em questão é ativo
				n_nodes_f[j]++;
			if (db.nodes[node - 1]->active_GL[j + 3] != 0)	//se o GL em questão é ativo
				n_nodes_m[j]++;
		}
	}
	//Calculando os multiplicadores
	for (int j = 0; j < 3; j++)
	{
		mult_f[j] = 1.0 / n_nodes_f[j];
		mult_m[j] = 1.0 / n_nodes_m[j];
	}
}


void NodalLoad::EvaluateExplicit(double t)
{
	int node;
	Matrix f(3);
	Matrix m(3);
	//For all nodes of the node set - Mounts contributions
	for (int index = 0; index < db.node_sets[node_set - 1]->n_nodes; index++)
	{
		//Node number
		node = db.node_sets[node_set - 1]->node_list[index];
		//Evaluating input data at current time
		for (int i = 0; i < 3; i++)
		{
			f(i, 0) = mult_f[i] * GetValueAt(t, i);
			m(i, 0) = mult_m[i] * GetValueAt(t, i + 3);
		}
		//Coordinate transformation (to material (node))
		f = transp(*db.nodes[node - 1]->Q0)*transp(*db.CS[cs - 1]->Q)*f;
		m = transp(*db.nodes[node - 1]->Q0)*transp(*db.CS[cs - 1]->Q)*m;
		//Global contributions
		int GL_lin;
		//Forces
		for (int lin = 0; lin < 3; lin++)
		{
			GL_lin = db.nodes[node - 1]->GLs[lin];
			if (db.nodes[node - 1]->active_GL[lin] == 1)
			{
				if (GL_lin > 0)	//Grau de liberdade livre e ativo
				{
					db.global_P_A(GL_lin - 1, 0) += 1.0*f(lin, 0);
					db.global_I_A(GL_lin - 1, 0) += 1.0*f(lin, 0);
				}	
				else
					db.global_P_B(-GL_lin - 1, 0) += 1.0*f(lin, 0);
			}
		}
		//Moments
		for (int lin = 0; lin < 3; lin++)
		{
			GL_lin = db.nodes[node - 1]->GLs[lin + 3];
			if (db.nodes[node - 1]->active_GL[lin + 3] == 1)
			{
				if (GL_lin > 0)	//Grau de liberdade livre e ativo
				{
					db.global_P_A(GL_lin - 1, 0) += 1.0*m(lin, 0);
					db.global_I_A(GL_lin - 1, 0) += 1.0*m(lin, 0);
				}
				else
					db.global_P_B(-GL_lin - 1, 0) += 1.0*m(lin, 0);
			}
		}
	}
}

void NodalLoad::Mount()
{
	int node;
	Matrix f(3);
	Matrix m(3);
	Matrix I(3, 3);
	I(0, 0) = 1.0; I(1, 1) = 1.0; I(2, 2) = 1.0;			//Identidade de ordem 3
	//For all nodes of the node set - Mounts contributions
	for (int index = 0; index < db.node_sets[node_set - 1]->n_nodes; index++)
	{
		//Node number
		node = db.node_sets[node_set - 1]->node_list[index];
		//Evaluating input data at current time
		for (int i = 0; i < 3; i++)
		{
			f(i, 0) = mult_f[i] * GetValueAt(db.last_converged_time + db.current_time_step, i);
			m(i, 0) = mult_m[i] * GetValueAt(db.last_converged_time + db.current_time_step, i + 3);
		}
		//Coordinate transformation (to global)
		f = transp(*db.CS[cs - 1]->Q)*f;
		m = transp(*db.CS[cs - 1]->Q)*m;
		//Calcula e insere o momento concentrado nos nós, de acordo com a rotação pré-existente nesses nós (procedimento de pseudo-momentos)
		//Rotação do nó em questão:
		Matrix alpha_delta(3, 1);
		alpha_delta(0, 0) = db.nodes[node - 1]->displacements[3];
		alpha_delta(1, 0) = db.nodes[node - 1]->displacements[4];
		alpha_delta(2, 0) = db.nodes[node - 1]->displacements[5];
		//Calculando o operador Xi
		double alpha_escalar = norm(alpha_delta);				//Valor escalar do parametro alpha
		Matrix A = skew(alpha_delta);							//Matriz A
		double g = 4.0 / (4.0 + alpha_escalar*alpha_escalar);	//função g(alpha) - em algumas ref. tb. chamado de h(alpha)
		Matrix Xi = g*(I + 0.5*A);
		//Operador Xi:
		m = transp(Xi)*m;
		//Contribuição para a matriz de rigidez (somente devido ao momento)
		Matrix stiff_moment = V(alpha_delta, m, alpha_escalar);
		//Global contributions
		int GL_lin;
		int GL_col;
		//Forces
		for (int lin = 0; lin < 3; lin++)
		{
			GL_lin = db.nodes[node - 1]->GLs[lin];
			if (db.nodes[node - 1]->active_GL[lin] == 1)
			{
				if (GL_lin > 0)	//Grau de liberdade livre e ativo
					db.global_P_A(GL_lin - 1, 0) += -1.0*f(lin, 0);
				else
					db.global_P_B(-GL_lin - 1, 0) += -1.0*f(lin, 0);
			}
		}
		//Moments
		for (int lin = 0; lin < 3; lin++)
		{
			GL_lin = db.nodes[node - 1]->GLs[lin + 3];
			if (db.nodes[node - 1]->active_GL[lin + 3] == 1)
			{
				if (GL_lin > 0)	//Grau de liberdade livre e ativo
					db.global_P_A(GL_lin - 1, 0) += -1.0*m(lin, 0);
				else
					db.global_P_B(-GL_lin - 1, 0) += -1.0*m(lin, 0);
			}
			for (int col = 0; col < 3; col++)
			{
				GL_col = db.nodes[node - 1]->GLs[col + 3];
				if (db.nodes[node - 1]->active_GL[col + 3] == 1)
				{
					if (GL_lin > 0 && GL_col > 0)	//Kaa
						db.global_stiffness_AA.setValue(GL_lin - 1, GL_col - 1, -1.0*stiff_moment(lin, col));
					if (GL_lin < 0 && GL_col > 0)	//Kba
						db.global_stiffness_BA.setValue(-GL_lin - 1, GL_col - 1, -1.0*stiff_moment(lin, col));
					if (GL_lin > 0 && GL_col < 0)	//Kab
						db.global_stiffness_AB.setValue(GL_lin - 1, -GL_col - 1, -1.0*stiff_moment(lin, col));
					if (GL_lin < 0 && GL_col < 0)	//Kbb
						db.global_stiffness_BB.setValue(-GL_lin - 1, -GL_col - 1, -1.0*stiff_moment(lin, col));
				}
			}
		}
	}
}
