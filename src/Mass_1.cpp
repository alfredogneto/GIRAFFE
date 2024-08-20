#include "Mass_1.h"
#include"Database.h"
//Variáveis globais
extern
Database db;

Mass_1::Mass_1()
{
	strain_energy = 0.0;
	kinetic_energy = 0.0;
	potential_gravitational_energy = 0.0;

	nDOFs = 3;
	material = 0;
	section = 0;
	n_nodes = 1;
	number = 0;
	nodes = new int[n_nodes];
	VTK_nodes = new int[n_nodes];
	DOFs = new int *[n_nodes];
	for (int i = 0; i < n_nodes;i++)
		DOFs[i] = new int[db.number_GLs_node];
	type_name = new char[20];//Nome do tipo do elemento
	sprintf(type_name, "Mass_1");
	VTK_type = 1;
	VTK_nodes[0] = 0;
	
	//Rotina para ativar os GLS de cada nó do elemento
	for (int i = 0; i < n_nodes; i++)
	{
		for (int j = 0; j < db.number_GLs_node; j++)
		{
			DOFs[i][j] = 0;
		}
		DOFs[i][0] = 1;
		DOFs[i][1] = 1;
		DOFs[i][2] = 1;
	}

	//Zerando coeficientes do elemento
	m = 0.0;
	
	c_mass = Matrix(3, 3);									//Matriz de massa
	c_loading = Matrix(3, 1);								//Vetor de esforços
	I3 = Matrix(3, 3);
	I3(0, 0) = 1.0;
	I3(1, 1) = 1.0;
	I3(2, 2) = 1.0;
	r = 0;
}

Mass_1::~Mass_1()
{
	delete[] nodes;
	delete[] VTK_nodes;
	if (DOFs != NULL)
	{
		for (int i = 0; i < n_nodes; i++)
			delete[] DOFs[i];
		delete[] DOFs;
	}
	delete[]type_name;
}

bool Mass_1::Read(FILE *f)
{
	char s[1000];
	fscanf(f, "%s", s);
	number = atoi(s);

	fscanf(f, "%s", s);
	if (!strcmp(s, "Mass"))
	{
		fscanf(f, "%s", s);
		m = atof(s);
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "Node"))
	{
		for (int n = 0; n < n_nodes; n++)
		{
			fscanf(f, "%s", s);
			nodes[n] = atoi(s);
		}
	}
	else
		return false;
	return true;
}

//Checa inconsistências no elemento para evitar erros de execução
bool Mass_1::Check()
{
	for (int i = 0; i < n_nodes; i++)
	{
		if (nodes[i] > db.number_nodes)
			return false;
	}
	return true;
}

void Mass_1::Write(FILE *f)
{
	fprintf(f, "Mass_1\t%d\tMass\t%.6e\tNode\t%d\n",number,m,nodes[0]);
}
//Escreve arquivo de resultados
void Mass_1::WriteResults(FILE *f)
{
	//DOES NOTHING
}

//Escreve no monitor do elemento//Escreve no monitor do elemento
void Mass_1::WriteMonitor(FILE *f, bool first_record, double time)
{
	//DOES NOTHING
}

void Mass_1::WriteVTK_XMLBase(std::vector<float> *float_vector)
{
	//Imprime os resultados do elemento
	int res_element = 0;
	//Imprime valores nulos para resultados que não fazem sentido para esse tipo de elemento
	for (int i = res_element; i < db.n_element_results; i++)
		float_vector->push_back(0.0);
}

void Mass_1::WriteVTK_XMLRender(FILE *f)
{
	//vetores para escrita no formato binário - usando a função 'enconde'
	std::vector<float> float_vector;
	std::vector<int> int_vector;

	double p0x, p0y, p0z, px, py, pz;
	//p0x, p0y,p0z - centro da esfera
	p0x = db.nodes[nodes[0] - 1]->copy_coordinates[0];
	p0y = db.nodes[nodes[0] - 1]->copy_coordinates[1];
	p0z = db.nodes[nodes[0] - 1]->copy_coordinates[2];

	int n = 10;		//discretização da renderização
	
	int sq = n*n;
	int points = (n + 1)*(n + 1);

	
	fprintf(f, "<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n<Points>\n<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"binary\">\n", points, 2 * sq);
	double theta = 0, phi = 0;
	float_vector.clear();
	for (int i = 0; i <= n; i = i + 1)
	{
		theta = i * 2 * PI / n;
		px = p0x + r*sin(phi)*cos(theta);
		py = p0y + r*sin(phi)*sin(theta); 
		pz = p0z + r*cos(phi);
		float_vector.push_back((float)(px));
		float_vector.push_back((float)(py));
		float_vector.push_back((float)(pz));
		for (int j = 1; j <= n; j++)
		{
			phi = j*PI / n;
			px = p0x + r*sin(phi)*cos(theta);
			py = p0y + r*sin(phi)*sin(theta);
			pz = p0z + r*cos(phi);
			float_vector.push_back((float)(px));
			float_vector.push_back((float)(py));
			float_vector.push_back((float)(pz));
		}
		phi = 0;
	}
	fprintf(f, encodeData<float>(float_vector).c_str());
	fprintf(f, "\n");
	fprintf(f, "\n</DataArray>\n</Points>\n<Cells><DataArray type=\"Int32\" Name=\"connectivity\" format=\"binary\">\n");
	int_vector.clear();
	for (int l = 1; l <= n; l++)
	{
		for (int k = (l - 1)*(n + 1) + 1; k<(n + 1)*l; k++)
		{
			int_vector.push_back(k - 1);
			int_vector.push_back(k);
			int_vector.push_back(k + n);
			int_vector.push_back(k);
			int_vector.push_back(k+n);
			int_vector.push_back(k + n + 1);

			//fprintf(f, "%d %d %d\n", k - 1, k, k + n);
			//fprintf(f, "%d %d %d\n", k, k + n, k + n + 1);
		}
	}
	fprintf(f, encodeData(int_vector).c_str());
	fprintf(f, "\n");
	int_vector.clear();
	fprintf(f, "</DataArray>\n<DataArray type=\"Int32\" Name=\"offsets\" format=\"binary\">");

	for (int w = 0; w<12 * sq; w++)
	{
		int_vector.push_back(3 + (3 * w));
		//fprintf(f, "%d ", 3 + (3 * w));
	}
	fprintf(f, encodeData(int_vector).c_str());
	fprintf(f, "\n");
	int_vector.clear();
	fprintf(f, "</DataArray>\n<DataArray type=\"Int32\" Name=\"types\" format=\"binary\">");
	for (int w = 0; w<6 * sq; w++)
	{
		int_vector.push_back(5);
		//fprintf(f, "5 ");
	}
	fprintf(f, encodeData(int_vector).c_str());
	fprintf(f, "\n");
	fprintf(f, "</DataArray>\n</Cells>\n");
	
	/////////////////////////////////////////////////////////////////////
	//Opens CellData
	fprintf(f, "\t\t\t<CellData FieldData=\"ElementData\">\n");
	float_vector.clear();
	//Opens DataArray
	int n_cells = 2 * sq;
	fprintf(f, "\t\t\t\t<DataArray Name=\"ElementResults\" type=\"Float32\" NumberOfComponents=\"%d\" format=\"binary\">\n", db.n_element_results);
	for (int cell = 0; cell < n_cells; cell++)
	{
		//Imprime os resultados do elemento
		int res_element = 0;
		//Imprime valores nulos para resultados que não fazem sentido para esse tipo de elemento
		for (int i = res_element; i < db.n_element_results; i++)
			float_vector.push_back(0.0);
	}
	fprintf(f, encodeData(float_vector).c_str());
	fprintf(f, "\n");
	//Closes DataArray
	fprintf(f, "\t\t\t\t</DataArray>\n");
	int_vector.clear();
	//Opens DataArray
	fprintf(f, "\t\t\t\t<DataArray Name=\"ElementProperties\" type=\"Int32\" NumberOfComponents=\"%d\" format=\"binary\">\n", 4);
	for (int cell = 0; cell < n_cells; cell++)
	{
		int_vector.push_back(4);		//Element ID
		int_vector.push_back(0);
		int_vector.push_back(0);
		int_vector.push_back(0);
	}
	fprintf(f, encodeData(int_vector).c_str());
	fprintf(f, "\n");
	//Closes DataArray
	fprintf(f, "\t\t\t\t</DataArray>\n");
	//Closes CellData
	fprintf(f, "\t\t\t</CellData>\n");
	/////////////////////////////////////////////////////////////////////
	fprintf(f,"</Piece>\n");
}

//Monta carregamentos associados ao elemento
void Mass_1::MountElementLoads()
{
	zeros(&c_loading);
	//////////////////Efeito do peso próprio////////////////////////////
	if (db.environment != NULL)
	{
		//Se existe campo gravitacional
		if (db.environment->g_exist == true)
		{
			load_multiplier = 1.0;
			l_factor = db.environment->bool_g.GetLinearFactorAtCurrentTime();
			mult = l_factor*load_multiplier*m;
			c_loading(0, 0) = -mult * db.environment->G(0, 0);
			c_loading(1, 0) = -mult * db.environment->G(1, 0);
			c_loading(2, 0) = -mult * db.environment->G(2, 0);
		}
	}
}

//Monta elementos
void Mass_1::Mount()
{
	//DOES NOTHING
}
//Monta matriz de transformação de coordenadas
void Mass_1::TransformMatrix()
{
	//DOES NOTHING
}
//Preenche a contribuição do elemento nas matrizes globais
void Mass_1::MountGlobal()
{
	int GL_global_1 = 0;
	int GL_global_2 = 0;
	double anterior = 0;
	for (int i = 0; i < 3; i++)
	{
		GL_global_1 = db.nodes[nodes[0] - 1]->GLs[i];
		//Caso o grau de liberdade seja livre:
		if (GL_global_1 > 0)
		{
			anterior = db.global_P_A(GL_global_1 - 1, 0);
			db.global_P_A(GL_global_1 - 1, 0) = anterior + c_loading(i, 0);
		}
		else
		{
			anterior = db.global_P_B(-GL_global_1 - 1, 0);
			db.global_P_B(-GL_global_1 - 1, 0) = anterior + c_loading(i, 0);
		}
		for (int j = 0; j < 3; j++)
		{
			GL_global_2 = db.nodes[nodes[0] - 1]->GLs[j];
			//Caso os graus de liberdade sejam ambos livres (Matriz Kaa)
			if (GL_global_1 > 0 && GL_global_2 > 0)
				db.global_stiffness_AA.setValue(GL_global_1 - 1, GL_global_2 - 1, c_mass(i, j));
			//Caso os graus de liberdade sejam ambos fixos (Matriz Kbb)
			if (GL_global_1 < 0 && GL_global_2 < 0)
				db.global_stiffness_BB.setValue(-GL_global_1 - 1, -GL_global_2 - 1, c_mass(i, j));
			//Caso os graus de liberdade sejam livre e fixo (Matriz Kab)
			if (GL_global_1 > 0 && GL_global_2 < 0)
				db.global_stiffness_AB.setValue(GL_global_1 - 1, -GL_global_2 - 1, c_mass(i, j));
			//Caso os graus de liberdade sejam fixo e livre (Matriz Kba)
			if (GL_global_1 < 0 && GL_global_2 > 0)
				db.global_stiffness_BA.setValue(-GL_global_1 - 1, GL_global_2 - 1, c_mass(i, j));
		}
	}
}
//Salva variáveis nos pontos de Gauss úteis para descrição lagrangiana atualizada
void Mass_1::SaveLagrange()
{
	//DOES NOTHING
}
//Pré-cálculo de variáveis que é feito uma única vez no início
void Mass_1::PreCalc()
{
	//Tenta tomar como referência algum material existente no modelo
	double rhomax = 0.0;
	for (int i = 0; i < db.number_materials; i++)
	{
		if (rhomax < db.materials[i]->rho)
			rhomax = db.materials[i]->rho;
	}
	if (rhomax != 0.0)
		r = pow(3 * m / (4 * PI*rhomax), 0.3333333333333333);
	else
		r = 0.2*db.EvaluateBoundingBoxDiag() / db.number_nodes;		//r = raio
}

//Monta a matriz de massa
void Mass_1::MountMass()
{
	c_mass = m*I3;
}
//Monta a matriz de massa
void Mass_1::MountMassModal()
{
	c_mass = m*I3;
}

//Monta a matriz de amortecimento para realização da análise modal
void Mass_1::MountDampingModal()
{
	Zeros();
}

//Monta a matriz de amortecimento
void Mass_1::MountDamping(bool update_rayleigh)
{
	//DOES NOTHING
}

//Montagens - Newmark
void Mass_1::MountDyn()
{
	Matrix accel(3);
	for (int ind = 0; ind < 3; ind++)
		accel(ind, 0) = db.nodes[nodes[0] - 1]->accel[ind];
	//Modificações da dinâmica nos esforços - presença das forças de inércia
	c_loading = c_loading + c_mass*accel;
	//Modificações da dinâmica na matriz de rigidez:
	Dynamic* ptr_sol = static_cast<Dynamic*>(db.solution[db.current_solution_number - 1]);
	c_mass = ptr_sol->a1*c_mass;
}
//Zera matrizes locais do elemento
void Mass_1::Zeros()
{
	c_mass.clear();
	c_loading.clear();
	kinetic_energy = 0.0;
	strain_energy = 0.0;
	potential_gravitational_energy = 0.0;
}
//Montagens para análise modal - inserção da matriz de massa e amortecimento na matriz de rigidez para posterior montagem global
void Mass_1::MountDynModal()
{
	//DOES NOTHING
}