#include "NodalFollowerLoad.h"

#include "Table.h"
#include "MathCode.h"
#include "NodeSet.h"
#include "Node.h"
#include "CoordinateSystem.h"
#include "Encoding.h"
#include "IO.h"
#include"Database.h"
//Variáveis globais
extern
Database db;

NodalFollowerLoad::NodalFollowerLoad()
{
	//Auxiliary variables
	Q = Matrix(3, 3);
	Xi = Matrix(3, 3);
	I = Matrix(3, 3);
	A = Matrix(3, 3);
	alpha = Matrix(3);
	g = 0;
	alpha_escalar = 0;
	q = NULL;
	dqdd = Matrix(6, 6);
	f = Matrix(3);
	m = Matrix(3);

	number = 0;
	table = NULL;
	mcode = NULL;
	n_times = 0;
	n_values = 6;

	I(0, 0) = 1.0;
	I(1, 1) = 1.0;
	I(2, 2) = 1.0;

	n_nodes_f = new int[3];			//número de nós para divisão de forças - 3 componentes
	n_nodes_m = new int[3];			//número de nós para divisão de momentos - 3 componentes
	mult_f = new double[3];			//multiplicador para os esforços de força
	mult_m = new double[3];			//multiplicador para os esforços de momento
}

NodalFollowerLoad::~NodalFollowerLoad()
{
	if (q != NULL)
	{
		for (int i = 0; i < n_nodes_copy; i++)
			delete q[i];
		delete[] q;
	}
	if (table != NULL)
		delete table;
	if (mcode != NULL)
		delete mcode;

	delete[]n_nodes_f;
	delete[]n_nodes_m;
	delete[]mult_f;
	delete[]mult_m;
}

bool NodalFollowerLoad::Read(FILE *f)
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
bool NodalFollowerLoad::Check()
{
	if (cs > db.number_CS)
		return false;
	if (node_set > db.number_node_sets)
		return false;
	return true;
}

void NodalFollowerLoad::Write(FILE *f)
{
	fprintf(f, "NodalFollowerLoad\t%d\tNodeSet\t%d\tCS\t%d\t", number, node_set, cs);
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

//Realiza pré-cálculos do esforço seguidor
void NodalFollowerLoad::PreCalc()
{
	n_nodes_copy = db.node_sets[node_set - 1]->n_nodes;
	//Alocando q
	q = new Matrix*[n_nodes_copy];
	for (int i = 0; i < db.node_sets[node_set - 1]->n_nodes; i++)
		q[i] = new Matrix(6);
}

//Atualiza dados necessários e que sejam dependentes de DOFs ativos/inativos - chamado no início de cada solution step
void NodalFollowerLoad::UpdateforSolutionStep()
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

void NodalFollowerLoad::EvaluateExplicit(double t)
{
	int node;
	Matrix Qi;
	//For all nodes of the node set - Mounts contributions
	for (int index = 0; index < db.node_sets[node_set - 1]->n_nodes; index++)
	{
		//Node number
		node = db.node_sets[node_set - 1]->node_list[index];
		
		//Coordinate transformation inclusion
		Qi = (*db.nodes[node - 1]->Q) * transp(*db.CS[cs - 1]->Q);
		//Rotations 
		alpha_escalar = norm(*db.nodes[node-1]->alpha);
		A = skew(alpha);
		g = 4.0 / (4.0 + alpha_escalar * alpha_escalar);
		//Cálculo da matriz de rotação Q
		Q = I + g * (A + 0.5*(A*A));
		
		//Evaluating input data at current time
		for (int i = 0; i < 3; i++)
		{
			f(i, 0) = mult_f[i] * GetValueAt(t, i);
			m(i, 0) = mult_m[i] * GetValueAt(t, i + 3);
		}
		//Forces and moments on current (trial) configuration
		Matrix fip = Q*Qi*f;
		Matrix mip = Q*Qi*m;
		
		//Local contributions
		for (int i = 0; i < 3; i++)
		{
			(*q[index])(i, 0) = fip(i, 0);
			(*q[index])(i + 3, 0) = mip(i, 0);
		}
		//Global contributions
		int GL_lin;
		//int GL_col;
		for (int lin = 0; lin < 6; lin++)
		{
			GL_lin = db.nodes[node - 1]->GLs[lin];
			if (db.nodes[node - 1]->active_GL[lin] == 1)
			{
				if (GL_lin > 0)	//Grau de liberdade livre e ativo
				{
					db.global_P_A(GL_lin - 1, 0) += 1.0*(*q[index])(lin, 0);
					db.global_I_A(GL_lin - 1, 0) += 1.0*(*q[index])(lin, 0);
				}

				else
					db.global_P_B(-GL_lin - 1, 0) += 1.0*(*q[index])(lin, 0);
			}
		}
	}
}


//Calcula o esforço e o operador tangente
void NodalFollowerLoad::Mount()
{
	int node;
	Matrix Qi;
	//For all nodes of the node set - Mounts contributions
	for (int index = 0; index < db.node_sets[node_set - 1]->n_nodes; index++)
	{
		//Node number
		node = db.node_sets[node_set - 1]->node_list[index];
		//Rotations - last converged configuration
		alpha(0, 0) = db.nodes[node - 1]->copy_coordinates[3];
		alpha(1, 0) = db.nodes[node - 1]->copy_coordinates[4];
		alpha(2, 0) = db.nodes[node - 1]->copy_coordinates[5];
		alpha_escalar = norm(alpha);
		A = skew(alpha);
		g = 4.0 / (4.0 + alpha_escalar*alpha_escalar);
		Qi = I + g*(A + 0.5*(A*A));
		//Coordinate transformation inclusion
		Qi = Qi*transp(*db.CS[cs - 1]->Q);
		//Rotations - last incremental rotation (still not converged)
		alpha(0, 0) = db.nodes[node - 1]->displacements[3];
		alpha(1, 0) = db.nodes[node - 1]->displacements[4];
		alpha(2, 0) = db.nodes[node - 1]->displacements[5];
		alpha_escalar = norm(alpha);
		A = skew(alpha);
		g = 4.0 / (4.0 + alpha_escalar*alpha_escalar);
		//Cálculo da matriz de rotação Q
		Q = I + g*(A + 0.5*(A*A));
		//Calculando o operador Xi	
		Xi = g*(I + 0.5*A);
		//Evaluating input data at current time
		for (int i = 0; i < 3; i++)
		{
			f(i, 0) = mult_f[i] * GetValueAt(db.last_converged_time+db.current_time_step, i);
			m(i, 0) = mult_m[i] * GetValueAt(db.last_converged_time + db.current_time_step, i + 3);
		}
		//Forces and moments on current (trial) configuration
		Matrix fip = Q*Qi*f;
		Matrix mip = Xi*Qi*m;
		Matrix K12 = -1.0*skew(fip)*Xi;
		Matrix K22 = -0.5*g*(skew(Qi*m) + Xi*(dyadic(Qi*m, alpha)));
		//Local contributions
		for (int i = 0; i < 3; i++)
		{
			(*q[index])(i, 0) = fip(i, 0);
			(*q[index])(i + 3, 0) = mip(i, 0);
			for (int j = 0; j < 3; j++)
			{
				dqdd(i, j + 3) = K12(i, j);
				dqdd(i + 3, j + 3) = K22(i, j);
			}
		}
		//Global contributions
		int GL_lin;
		int GL_col;
		for (int lin = 0; lin < 6; lin++)
		{
			GL_lin = db.nodes[node - 1]->GLs[lin];
			if (db.nodes[node - 1]->active_GL[lin] == 1)
			{
				if (GL_lin > 0)	//Grau de liberdade livre e ativo
					db.global_P_A(GL_lin - 1, 0) += -1.0*(*q[index])(lin, 0);
				else
					db.global_P_B(-GL_lin - 1, 0) += -1.0*(*q[index])(lin, 0);
			}
			for (int col = 0; col < 6; col++)
			{
				GL_col = db.nodes[node - 1]->GLs[col];
				if (db.nodes[node - 1]->active_GL[col] == 1)
				{
					if (GL_lin > 0 && GL_col > 0)	//Kaa
						db.global_stiffness_AA.setValue(GL_lin - 1, GL_col - 1, -1.0*dqdd(lin, col));
					if (GL_lin < 0 && GL_col > 0)	//Kba
						db.global_stiffness_BA.setValue(-GL_lin - 1, GL_col - 1, -1.0*dqdd(lin, col));
					if (GL_lin > 0 && GL_col < 0)	//Kab
						db.global_stiffness_AB.setValue(GL_lin - 1, -GL_col - 1, -1.0*dqdd(lin, col));
					if (GL_lin < 0 && GL_col < 0)	//Kbb
						db.global_stiffness_BB.setValue(-GL_lin - 1, -GL_col - 1, -1.0*dqdd(lin, col));
				}
			}
		}
	}
}

void NodalFollowerLoad::WriteVTK_XML(FILE *f)
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
	
	for (int index = 0; index < db.node_sets[node_set - 1]->n_nodes; index++)
	{
		float_vector.push_back((float)((*q[index])(0, 0)));
		float_vector.push_back((float)((*q[index])(1, 0)));
		float_vector.push_back((float)((*q[index])(2, 0)));
	}
	fprintf(f, encodeData(float_vector).c_str());
	fprintf(f, "\n");
	//Closes DataArray
	fprintf(f, "\t\t\t\t</DataArray>\n");
	float_vector.clear();
	//Opens DataArray
	fprintf(f, "\t\t\t\t<DataArray Name = \"Moment\" type = \"Float32\" NumberOfComponents= \"3\" format=\"binary\">\n");
	for (int index = 0; index < db.node_sets[node_set - 1]->n_nodes; index++)
	{
		float_vector.push_back((float)((*q[index])(3, 0)));
		float_vector.push_back((float)((*q[index])(4, 0)));
		float_vector.push_back((float)((*q[index])(5, 0)));
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