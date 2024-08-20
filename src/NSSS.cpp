#include "NSSS.h"
#define PI 3.1415926535897932384626433832795
#include"Database.h"
//Variáveis globais
extern
Database db;

NSSS::NSSS()
{
	type_name = new char[20];		//Nome do tipo do contato
	sprintf(type_name, "NSSS");
	number = 0;
	n_SS = 0;
	n_NS = 0;
	mu = 0.0;
	ept = 0.0;
	epn = 0.0;
	cn = 0.0;
	ct = 0.0;
	pinball = 0.0;
	radius = 0.0;
	number_pointwise_interactions = 5;	
	

	//Abaixo as variáveis que dependem do número de elementos para serem alocadas - Alocação feita na função Alloc, quando chamada durante o PreCalc
	typeOK1 = NULL;
	typeOK2 = NULL;
	activate = NULL;
	i_loading = NULL;
	contact_stiffness = NULL;
	alloc_control = NULL;
	DOFs_surfaces = NULL;
	
	cd = NULL;
	xS_p = NULL;
	QS = NULL;
	alphaS = NULL;
	uS = NULL;

	I3 = new Matrix(3, 3);
	(*I3)(0, 0) = 1.0;
	(*I3)(1, 1) = 1.0;
	(*I3)(2, 2) = 1.0;

	number_surfaces = 0;
	number_nodes = 0;
}

void NSSS::Alloc()
{
	activate = new bool*[number_nodes];
	for (int i = 0; i < number_nodes; i++)
		activate[i] = new bool[number_surfaces];
	
	alloc_control = new bool**[number_nodes];
	for (int i = 0; i < number_nodes; i++)
	{
		alloc_control[i] = new bool*[number_surfaces];
		for (int j = 0; j < number_surfaces; j++)
			alloc_control[i][j] = new bool[number_pointwise_interactions];
	}

	cd = new NSContactData**[number_nodes];
	for (int i = 0; i < number_nodes; i++)
	{
		cd[i] = new NSContactData*[number_surfaces];
		for (int j = 0; j < number_surfaces; j++)
		{
			cd[i][j] = new NSContactData();
			cd[i][j]->n_solutions = number_pointwise_interactions;
			//Alocação de cada um dos cd[i][j] será feita de acordo com a proximidade de contato - função Pinball
		}
	}
	
	DOFs_surfaces = new int[number_surfaces];

	//Alocação parcial de matrizes
	i_loading = new Matrix***[number_nodes];
	for (int i = 0; i < number_nodes; i++)
	{
		i_loading[i] = new Matrix**[number_surfaces];
		for (int j = 0; j < number_surfaces; j++)
			i_loading[i][j] = new Matrix*[number_pointwise_interactions];
	}
		
	contact_stiffness = new double****[number_nodes];
	for (int i = 0; i < number_nodes; i++)
	{
		contact_stiffness[i] = new double***[number_surfaces];
		for (int j = 0; j < number_surfaces; j++)
			contact_stiffness[i][j] = new double**[number_pointwise_interactions];
	}
		
	xS_p = new Matrix**[number_nodes];
	for (int i = 0; i < number_nodes; i++)
	{
		xS_p[i] = new Matrix*[number_surfaces];
		for (int j = 0; j < number_surfaces; j++)
			xS_p[i][j] = new Matrix(3);
	}

	QS = new Matrix**[number_nodes];
	for (int i = 0; i < number_nodes; i++)
	{
		QS[i] = new Matrix*[number_surfaces];
		for (int j = 0; j < number_surfaces; j++)
			QS[i][j] = new Matrix(3, 3);
	}

	alphaS = new Matrix**[number_nodes];
	for (int i = 0; i < number_nodes; i++)
	{
		alphaS[i] = new Matrix*[number_surfaces];
		for (int j = 0; j < number_surfaces; j++)
			alphaS[i][j] = new Matrix(3);
	}

	uS = new Matrix**[number_nodes];
	for (int i = 0; i < number_nodes; i++)
	{
		uS[i] = new Matrix*[number_surfaces];
		for (int j = 0; j < number_surfaces; j++)
			uS[i][j] = new Matrix(3);
	}

	//Atribuição de valores iniciais
	for (int i = 0; i < number_nodes; i++)
	{
		for (int j = 0; j < number_surfaces; j++)
		{
			activate[i][j] = false;
			for (int k = 0; k < number_pointwise_interactions; k++)
				alloc_control[i][j][k] = false;	//Não há alocação no início
		}
	}
	for (int j = 0; j < number_surfaces; j++)
		DOFs_surfaces[j] = 0;
}

//Alocação específica para o contato node_index, surface_index
void NSSS::AllocSpecific(int node_index, int surface_index, int sol_index)
{
	//Caso não esteja ainda alocado, realiza a alocação
	if (alloc_control[node_index][surface_index][sol_index] == false)
	{
		i_loading[node_index][surface_index][sol_index] = new Matrix(6 + DOFs_surfaces[surface_index]);
		contact_stiffness[node_index][surface_index][sol_index] = new double*[6 + DOFs_surfaces[surface_index]];
		for (int i = 0; i < (6 + DOFs_surfaces[surface_index]); i++)
			contact_stiffness[node_index][surface_index][sol_index][i] = new double[6 + DOFs_surfaces[surface_index]];
		for (int ni = 0; ni < (6 + DOFs_surfaces[surface_index]); ni++)
		for (int nj = 0; nj < (6 + DOFs_surfaces[surface_index]); nj++)
			contact_stiffness[node_index][surface_index][sol_index][ni][nj] = 0.0;
		alloc_control[node_index][surface_index][sol_index] = true;
	}
}

//Liberação específica para o contato node_index, surface_index
void NSSS::FreeSpecific(int node_index, int surface_index, int sol_index)
{
	//Caso esteja alocado, realiza a desalocação
	if (alloc_control[node_index][surface_index][sol_index] == true)
	{
		delete i_loading[node_index][surface_index][sol_index];
		for (int i = 0; i < (6 + DOFs_surfaces[surface_index]); i++)
			delete[] contact_stiffness[node_index][surface_index][sol_index][i];
		delete[] contact_stiffness[node_index][surface_index][sol_index];
		alloc_control[node_index][surface_index][sol_index] = false;
	}
}

NSSS::~NSSS()
{
	delete[] type_name;
	delete I3;

	if (activate != NULL)
	{
		for (int i = 0; i < number_nodes; i++)
			delete [] activate[i];
		delete[]activate;
	}
	for (int i = 0; i < number_nodes;i++)
	for (int j = 0; j < number_surfaces;j++)
	for (int k = 0; k < number_pointwise_interactions;k++)
		FreeSpecific(i,j,k);

	if (i_loading != NULL)
	{
		for (int i = 0; i < number_nodes; i++)
		{
			for (int j = 0; j < number_surfaces; j++)
				delete[]i_loading[i][j];
			delete[] i_loading[i];
		}
		delete[]i_loading;
	}
	if (contact_stiffness != NULL)
	{
		for (int i = 0; i < number_nodes; i++)
		{
			for (int j = 0; j < number_surfaces; j++)
				delete[] contact_stiffness[i][j];
			delete[] contact_stiffness[i];
		}
		delete[]contact_stiffness;
	}
	if (alloc_control != NULL)
	{
		for (int i = 0; i < number_nodes; i++)
		{
			for (int j = 0; j < number_surfaces; j++)
				delete[]alloc_control[i][j];
			delete[] alloc_control[i];
		}
		delete[]alloc_control;
	}

	if (cd != NULL)
	{
		for (int i = 0; i < number_nodes; i++)
		{
			for (int j = 0; j < number_surfaces; j++)
				delete cd[i][j];
			delete[]cd[i];
		}
			
		delete[]cd;
	}

	if (DOFs_surfaces != NULL)
		delete[]DOFs_surfaces;

	if (xS_p != NULL)
	{
		for (int i = 0; i < number_nodes; i++)
		{
			for (int j = 0; j < number_surfaces; j++)
				delete xS_p[i][j];
			delete[] xS_p[i];
		}
		delete[]xS_p;
	}
	if (QS != NULL)
	{
		for (int i = 0; i < number_nodes; i++)
		{
			for (int j = 0; j < number_surfaces; j++)
				delete QS[i][j];
			delete[] QS[i];
		}
		delete[]QS;
	}
	if (alphaS != NULL)
	{
		for (int i = 0; i < number_nodes; i++)
		{
			for (int j = 0; j < number_surfaces; j++)
				delete alphaS[i][j];
			delete[] alphaS[i];
		}
		delete[]alphaS;
	}
	if (uS != NULL)
	{
		for (int i = 0; i < number_nodes; i++)
		{
			for (int j = 0; j < number_surfaces; j++)
				delete uS[i][j];
			delete[] uS[i];
		}
		delete[]uS;
	}
}

void NSSS::WriteVTK_XMLRender(FILE *f)
{
	//Plota superfícies
	int temp_surf = 0;
	for (int j = 0; j < number_surfaces; j++)
	{
		temp_surf = db.surface_sets[n_SS - 1]->surf_list[j];
		db.surfaces[temp_surf - 1]->WriteVTK_XMLRender(f);
	}
		
}

void NSSS::WriteVTK_XMLForces(FILE *f)
{
	//Plotagem de forças de contato
	if (db.post_files->WriteContactForces_flag == true)
	{
		//vetores para escrita no formato binário - usando a função 'enconde'
		std::vector<float> float_vector;
		std::vector<int> int_vector;
		int count = 0;
		for (int surf_index = 0; surf_index < number_surfaces; surf_index++)
		for (int node_index = 0; node_index < number_nodes; node_index++)
		for (int sol_index = 0; sol_index < number_pointwise_interactions; sol_index++)
			if (activate[node_index][surf_index] == true && cd[node_index][surf_index]->repeated[sol_index] == false && (cd[node_index][surf_index]->return_value[sol_index] == 0) && cd[node_index][surf_index]->g_n[sol_index] <= 0.0)
			count++;
		//Opens Piece
		fprintf(f, "\t\t<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", count, count);
		//Opens Points
		fprintf(f, "\t\t\t<Points>\n");
		//Opens DataArray
		fprintf(f, "\t\t\t\t<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"binary\">\n");
		float_vector.clear();
		//Preenchendo as coordenadas dos pontos
		for (int surf_index = 0; surf_index < number_surfaces; surf_index++)
		for (int node_index = 0; node_index < number_nodes; node_index++)
		for (int sol_index = 0; sol_index < number_pointwise_interactions; sol_index++)
		if (activate[node_index][surf_index] == true && cd[node_index][surf_index]->repeated[sol_index] == false && (cd[node_index][surf_index]->return_value[sol_index] == 0) && cd[node_index][surf_index] ->g_n[sol_index] <= 0.0)
		{
			float_vector.push_back((float)(*cd[node_index][surf_index]->G_p[sol_index])(0, 0));
			float_vector.push_back((float)(*cd[node_index][surf_index]->G_p[sol_index])(1, 0));
			float_vector.push_back((float)(*cd[node_index][surf_index]->G_p[sol_index])(2, 0));
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
		fprintf(f, "\t\t\t<PointData Vectors = \"ContactForces\">\n");
		Matrix normal(3);
		float_vector.clear();
		//Opens DataArray
		fprintf(f, "\t\t\t\t<DataArray Name = \"Normal\" type = \"Float32\" NumberOfComponents= \"3\" format=\"binary\">\n");
		for (int surf_index = 0; surf_index < number_surfaces; surf_index++)
		for (int node_index = 0; node_index < number_nodes; node_index++)
		for (int sol_index = 0; sol_index < number_pointwise_interactions; sol_index++)
		if (activate[node_index][surf_index] == true && cd[node_index][surf_index]->repeated[sol_index] == false && (cd[node_index][surf_index]->return_value[sol_index] == 0) && cd[node_index][surf_index]->g_n[sol_index] <= 0.0)
		{
			normal = -epn*cd[node_index][surf_index]->copy_g_n[sol_index] * (*cd[node_index][surf_index]->n_p[sol_index]);
			float_vector.push_back((float)normal(0, 0));
			float_vector.push_back((float)normal(1, 0));
			float_vector.push_back((float)normal(2, 0));
		}
		fprintf(f, encodeData(float_vector).c_str());
		fprintf(f, "\n");
		//Closes DataArray
		fprintf(f, "\t\t\t\t</DataArray>\n");
		Matrix friction(3);
		float_vector.clear();
		//Opens DataArray
		fprintf(f, "\t\t\t\t<DataArray Name = \"Friction\" type = \"Float32\" NumberOfComponents= \"3\" format=\"binary\">\n");
		for (int surf_index = 0; surf_index < number_surfaces; surf_index++)
		for (int node_index = 0; node_index < number_nodes; node_index++)
		for (int sol_index = 0; sol_index < number_pointwise_interactions; sol_index++)
		if (activate[node_index][surf_index] == true && cd[node_index][surf_index]->repeated[sol_index] == false && (cd[node_index][surf_index]->return_value[sol_index] == 0) && cd[node_index][surf_index]->g_n[sol_index] <= 0.0)
		{
			friction = -ept*(*cd[node_index][surf_index]->copy_g_t[sol_index]);
			float_vector.push_back((float)friction(0, 0));
			float_vector.push_back((float)friction(1, 0));
			float_vector.push_back((float)friction(2, 0));
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
}

bool NSSS::Read(FILE *f)
{
	char s[1000];
	fscanf(f, "%s", s);
	number = atoi(s);
	fscanf(f, "%s", s);
	if (!strcmp(s, "NodeSet"))
	{
		fscanf(f, "%s", s);
		n_NS = atoi(s);
	}
	else
		return false;
	fscanf(f, "%s", s);
	if (!strcmp(s, "SurfaceSet"))
	{
		fscanf(f, "%s", s);
		n_SS = atoi(s);
	}
	else
		return false;
	fscanf(f, "%s", s);
	if (!strcmp(s, "MU"))
	{
		fscanf(f, "%s", s);
		mu = atof(s);
	}
	else
		return false;
	fscanf(f, "%s", s);
	if (!strcmp(s, "EPN"))
	{
		fscanf(f, "%s", s);
		epn = atof(s);
	}
	else
		return false;
	fscanf(f, "%s", s);
	if (!strcmp(s, "CN"))
	{
		fscanf(f, "%s", s);
		cn = atof(s);
	}
	else
		return false;
	fscanf(f, "%s", s);
	if (!strcmp(s, "EPT"))
	{
		fscanf(f, "%s", s);
		ept = atof(s);
	}
	else
		return false;
	fscanf(f, "%s", s);
	if (!strcmp(s, "CT"))
	{
		fscanf(f, "%s", s);
		ct = atof(s);
	}
	else
		return false;
	fscanf(f, "%s", s);
	if (!strcmp(s, "Pinball"))
	{
		fscanf(f, "%s", s);
		pinball = atof(s);
	}
	else
		return false;
	fscanf(f, "%s", s);
	if (!strcmp(s, "Radius"))
	{
		fscanf(f, "%s", s);
		radius = atof(s);
	}
	else
		return false;
	fscanf(f, "%s", s);
	if (!strcmp(s, "MaxPointwiseInt"))
	{
		fscanf(f, "%s", s);
		number_pointwise_interactions = atoi(s);
	}
	else
		return false;
	//Salva a posição (stream)
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

void NSSS::Write(FILE *f)
{
	fprintf(f, "NSSS\t%d\tNodeSet\t%d\tSurfaceSet\t%d\tMU\t%.6e\tEPN\t%.6e\tCN\t%.6e\tEPT\t%.6e\tCT\t%.6e\tPinball\t%.6e\tRadius\t%.6e\tMaxPointwiseInt\t%d\t",
		number, n_NS, n_SS, mu, epn, cn, ept, ct, pinball, radius,number_pointwise_interactions);
	bool_table.Write(f);
}

//Checa inconsistências no elemento para evitar erros de execução
bool NSSS::Check()
{
	if (n_SS > db.number_surface_sets || n_NS > db.number_node_sets)
	{
		printf("Error in Contact %d. Invalid Surface Set or Node Set assigned. Check the input file.\n", number);
		return false;
	}
	return true;
}

//Escreve arquivo de resultados
void NSSS::WriteResults(FILE *f)
{
	//TODO
}

//Escreve no monitor do contato
void NSSS::WriteMonitor(FILE *f, bool first_record, double time)
{
	//Cabeçalho
	if (first_record == true)
	{
		fprintf(f, "\t");
		for (int i = 0; i < number_nodes; i++)
			fprintf(f, "Node %d\t\t\t\t\t\t\t\t\t", db.node_sets[n_NS - 1]->node_list[i]);
		fprintf(f, "\n");
		fprintf(f, "TIME");
		for (int i = 0; i < number_nodes; i++)
			fprintf(f, "\tFnx\tFny\tFnz\tFn\tFatx\tFaty\tFatz\tFat\tNForces");
		fprintf(f, "\n");
	}
	fprintf(f, "%.6e\t", time);	//TIME
	for (int i = 0; i < number_nodes; i++)
	{
		int count = 0;
		Matrix Fn(3);
		Matrix Fat(3);
		for (int j = 0; j < number_surfaces; j++)
		{

			if (activate[i][j] == 1)
			{
				for (int k = 0; k < number_pointwise_interactions; k++)
				{
					if (cd[i][j]->return_value[k] == 0 && cd[i][j]->g_n[k] < 0.0 && cd[i][j]->repeated[k] == false)
					{
						count++;
						Fn = Fn - epn*cd[i][j]->copy_g_n[k] * (*cd[i][j]->n_p[k]);
						Fat = Fat - ept*(*cd[i][j]->copy_g_t[k]);
					}
				}
				
			}
		}
		fprintf(f, "%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%d\t", Fn(0, 0), Fn(1, 0), Fn(2, 0), norm(Fn), Fat(0, 0), Fat(1, 0), Fat(2, 0), norm(Fat),count);
	}
		
	fprintf(f, "\n");
}

void NSSS::Mount()
{
#pragma omp parallel 
	{
		//Loop in nodes to perform contact Mount for each one
#pragma omp for
		for (int surf_index = 0; surf_index < number_surfaces; surf_index++)
		{
			int tempsurf = db.surface_sets[n_SS - 1]->surf_list[surf_index];
			for (int node_index = 0; node_index < number_nodes; node_index++)
			{
				if (activate[node_index][surf_index] == true)
				{
					int tempnode = db.node_sets[n_NS - 1]->node_list[node_index];
					(*xS_p[node_index][surf_index])(0, 0) = db.nodes[tempnode - 1]->copy_coordinates[0] + db.nodes[tempnode - 1]->displacements[0];
					(*xS_p[node_index][surf_index])(1, 0) = db.nodes[tempnode - 1]->copy_coordinates[1] + db.nodes[tempnode - 1]->displacements[1];
					(*xS_p[node_index][surf_index])(2, 0) = db.nodes[tempnode - 1]->copy_coordinates[2] + db.nodes[tempnode - 1]->displacements[2];
					//Seta na primeira posição do 'convective' os últimos valores convergidos
					cd[node_index][surf_index]->convective[0][0] = cd[node_index][surf_index]->copy_convective[0][0];
					cd[node_index][surf_index]->convective[0][1] = cd[node_index][surf_index]->copy_convective[0][1];
					//Chute inicial para a determinação da projeção ortogonal - seta somente na primeira posição do convective
					if (cd[node_index][surf_index]->copy_g_n[0] == 1.0)
						db.surfaces[tempsurf - 1]->InitialGuess(xS_p[node_index][surf_index], cd[node_index][surf_index]->convective, cd[node_index][surf_index]->n_solutions);
					//Solução do problema de mínima distância
					db.surfaces[tempsurf - 1]->FindMinimimumParameters(xS_p[node_index][surf_index], cd[node_index][surf_index]);
					//Se a projeção cair na faixa do domínio de interesse ou nas proximidades
					//Varredura nos possíveis pontos de interação pointwise
					for (int sol_index = 0; sol_index < number_pointwise_interactions; sol_index++)
					{
						//Se não for uma raiz repetida (já computada anteriormente)
						if (cd[node_index][surf_index]->repeated[sol_index] == false)
						{
							//printf("not repeated node %d surf %d sol_index %d\n", tempnode, tempsurf, sol_index);
							if (cd[node_index][surf_index]->return_value[sol_index] == 0 || cd[node_index][surf_index]->return_value[sol_index] == 3)
							{
								//Cálculo da função gap (escalar)
								db.surfaces[tempsurf - 1]->Gamma_and_Triad(cd[node_index][surf_index]->G_p[sol_index], cd[node_index][surf_index]->t1_p[sol_index],
									cd[node_index][surf_index]->t2_p[sol_index], cd[node_index][surf_index]->n_p[sol_index], cd[node_index][surf_index]->G_i[sol_index],
									cd[node_index][surf_index]->t1_i[sol_index], cd[node_index][surf_index]->t2_i[sol_index], cd[node_index][surf_index]->n_i[sol_index],
									cd[node_index][surf_index]->G_ip[sol_index], &cd[node_index][surf_index]->copy_convective[sol_index][0],
									&cd[node_index][surf_index]->copy_convective[sol_index][1], &cd[node_index][surf_index]->convective[sol_index][0], &cd[node_index][surf_index]->convective[sol_index][1]);
								cd[node_index][surf_index]->g_n[sol_index] = dot(*cd[node_index][surf_index]->n_p[sol_index], *xS_p[node_index][surf_index] - *cd[node_index][surf_index]->G_p[sol_index]) - radius;
							}
							else
								cd[node_index][surf_index]->g_n[sol_index] = 1.0;//valor para indicar que não há contato
							//se houver penetração e o range de coordenadas convectivas for válido monta a contribuição do contato na forma fraca e operador tangente
							if (cd[node_index][surf_index]->g_n[sol_index] <= 0.0 && (cd[node_index][surf_index]->return_value[sol_index] == 0 || cd[node_index][surf_index]->return_value[sol_index] == 3))
							{
								//ReportContact(node_index, surf_index, sol_index);
								//Aloca espaço para salvar a matriz de rigidez e vetor de esforços
								AllocSpecific(node_index, surf_index,sol_index);
								//Se o contato já tiver sido estabelecido - avalia gap tangencial. Caso contrário, o gap tangencial é nulo.
								if (cd[node_index][surf_index]->copy_g_n[sol_index] < 0.0)
									EvaluateTangentialGap(node_index, surf_index, sol_index);
								else
								{
									//Seta cópias de convectivas iguais às atuais, para que o gap tangencial seja nulo
									cd[node_index][surf_index]->copy_convective[sol_index][0] = cd[node_index][surf_index]->convective[sol_index][0];
									cd[node_index][surf_index]->copy_convective[sol_index][1] = cd[node_index][surf_index]->convective[sol_index][1];
								}
								//Verificação de stick/slide
								double FatMax = mu*epn*abs(cd[node_index][surf_index]->g_n[sol_index]);
								double FatTry = ept*norm(*cd[node_index][surf_index]->g_t[sol_index]);
								//Monta esforços de contato e rigidez - efeito da penetração e atrito
								//if (FatTry < FatMax)//sticking
								if (FatTry < FatMax && cd[node_index][surf_index]->copy_g_n[sol_index] < 0.0)//sticking e já havia contato prévio (mesmo de elemento vizinho)
								{
									//Função contact stick
									db.surfaces[tempsurf - 1]->ContactSphereSurfaceSticking(i_loading[node_index][surf_index][sol_index]->getMatrix(), contact_stiffness[node_index][surf_index][sol_index], cd[node_index][surf_index]->convective[sol_index][0], cd[node_index][surf_index]->convective[sol_index][1], cd[node_index][surf_index]->copy_convective[sol_index][0], cd[node_index][surf_index]->copy_convective[sol_index][1], cd[node_index][surf_index]->copy_g_t[sol_index]->getMatrix(), tempnode, &epn, &ept, &cn, &ct, &mu, &radius);
								}
								else
								{
									//if (mu == 0.0)
									if (mu == 0.0 || cd[node_index][surf_index]->copy_g_n[sol_index] >= 0.0)//atrito nulo ou primeiro contato
									{
										double temp_ept = 0.0;
										double temp_ct = 0.0;
										//Função contact stick - chamada com ept e ct nulos
										db.surfaces[tempsurf - 1]->ContactSphereSurfaceSticking(i_loading[node_index][surf_index][sol_index]->getMatrix(), contact_stiffness[node_index][surf_index][sol_index], cd[node_index][surf_index]->convective[sol_index][0], cd[node_index][surf_index]->convective[sol_index][1], cd[node_index][surf_index]->copy_convective[sol_index][0], cd[node_index][surf_index]->copy_convective[sol_index][1], cd[node_index][surf_index]->copy_g_t[sol_index]->getMatrix(), tempnode, &epn, &temp_ept, &cn, &temp_ct, &mu, &radius);
										zeros(cd[node_index][surf_index]->g_t[sol_index]);
									}
									else
									{
										//Função contact slide
										db.surfaces[tempsurf - 1]->ContactSphereSurfaceSliding(i_loading[node_index][surf_index][sol_index]->getMatrix(), contact_stiffness[node_index][surf_index][sol_index], cd[node_index][surf_index]->convective[sol_index][0], cd[node_index][surf_index]->convective[sol_index][1], cd[node_index][surf_index]->copy_convective[sol_index][0], cd[node_index][surf_index]->copy_convective[sol_index][1], cd[node_index][surf_index]->copy_g_t[sol_index]->getMatrix(), tempnode, &epn, &ept, &cn, &ct, &mu, &radius);
										//Atualização do gap tangencial acumulado (elástico) - tirando o efeito do deslizamento ocorrido nesse passo
										double delta_lambda = (FatTry - FatMax) / ept;
										Matrix gTslide = (delta_lambda / norm(*cd[node_index][surf_index]->g_t[sol_index]))*(*cd[node_index][surf_index]->g_t[sol_index]);	//deslizamento
										*cd[node_index][surf_index]->g_t[sol_index] = *cd[node_index][surf_index]->g_t[sol_index] - gTslide;								//atualização do gap - tirando efeito de deslizamentos
									}
								}
							}
							else
								FreeSpecific(node_index, surf_index,sol_index);
						}
					}
				}
			}
		}
	}
}

//Montagens - Newmark
void NSSS::MountDyn()
{
	//NOT USED (dynamic effects added inside the function Mount())
}
//Preenche a contribuição do contato nas matrizes globais
void NSSS::MountGlobal()
{
	//Variáveis temporárias para salvar a indexação global dos graus de liberdade a serem setados na matriz de rigidez global
	int GL_global_1 = 0;
	int GL_global_2 = 0;
	double anterior = 0;
	int temp_node;
	int tempsurf;
	for (int surf_index = 0; surf_index < number_surfaces; surf_index++)
	{
		tempsurf = db.surface_sets[n_SS - 1]->surf_list[surf_index];
		for (int node_index = 0; node_index < number_nodes; node_index++)
		{
			if (activate[node_index][surf_index] == 1)
			{
				temp_node = db.node_sets[n_NS - 1]->node_list[node_index];	//Nó 1 - referente ao node set
				for (int sol_index = 0; sol_index < number_pointwise_interactions; sol_index++)
				{
					//Se não for uma raiz repetida (já computada anteriormente)
					if (cd[node_index][surf_index]->repeated[sol_index] == false)
					{
						//ReportContact(node_index, surf_index, sol_index);
						if (cd[node_index][surf_index]->return_value[sol_index] == 0 && cd[node_index][surf_index]->g_n[sol_index] < 0.0)
						{
							//ReportContact(node_index, surf_index, sol_index);
							for (int i = 0; i < (6 + DOFs_surfaces[surf_index]); i++)
							{
								//Nó da esfera
								if (i < 6)
									GL_global_1 = db.nodes[temp_node - 1]->GLs[i];
								else
									GL_global_1 = *db.surfaces[tempsurf - 1]->GLs[i - 6];
								//Caso o grau de liberdade seja livre:
								if (GL_global_1 > 0)
								{
									anterior = db.global_P_A(GL_global_1 - 1, 0);
									db.global_P_A(GL_global_1 - 1, 0) = anterior + (*i_loading[node_index][surf_index][sol_index])(i, 0);
									anterior = db.global_I_A(GL_global_1 - 1, 0);
									db.global_I_A(GL_global_1 - 1, 0) = anterior + (*i_loading[node_index][surf_index][sol_index])(i, 0);
								}
								else
								{
									if (GL_global_1 != 0)//se o GL é ativo
									{
										anterior = db.global_P_B(-GL_global_1 - 1, 0);
										db.global_P_B(-GL_global_1 - 1, 0) = anterior + (*i_loading[node_index][surf_index][sol_index])(i, 0);
									}

								}
								for (int j = 0; j < (6 + DOFs_surfaces[surf_index]); j++)
								{
									//Nó da esfera
									if (j < 6)
										GL_global_2 = db.nodes[temp_node - 1]->GLs[j];
									else
										GL_global_2 = *db.surfaces[tempsurf - 1]->GLs[j - 6];
									//Caso os graus de liberdade sejam ambos livres (Matriz Kaa)
									if (GL_global_1 > 0 && GL_global_2 > 0)
										db.global_stiffness_AA.setValue(GL_global_1 - 1, GL_global_2 - 1, contact_stiffness[node_index][surf_index][sol_index][i][j]);
									//Caso os graus de liberdade sejam ambos fixos (Matriz Kbb)
									if (GL_global_1 < 0 && GL_global_2 < 0)
										db.global_stiffness_BB.setValue(-GL_global_1 - 1, -GL_global_2 - 1, contact_stiffness[node_index][surf_index][sol_index][i][j]);
									//Caso os graus de liberdade sejam livre e fixo (Matriz Kab)
									if (GL_global_1 > 0 && GL_global_2 < 0)
										db.global_stiffness_AB.setValue(GL_global_1 - 1, -GL_global_2 - 1, contact_stiffness[node_index][surf_index][sol_index][i][j]);
									//Caso os graus de liberdade sejam fixo e livre (Matriz Kba)
									if (GL_global_1 < 0 && GL_global_2 > 0)
										db.global_stiffness_BA.setValue(-GL_global_1 - 1, GL_global_2 - 1, contact_stiffness[node_index][surf_index][sol_index][i][j]);
								}
							}
						}
					}
				}
			}
		}
	}
}

//Calcula a banda na matriz global devido a esse contato
void NSSS::Band(int* band_fixed, int* band_free)
{
	//NOT IMPLEMENTED
}

//Pré-cálculo de variáveis que é feito uma única vez no início
void NSSS::PreCalc()
{
	number_nodes = db.node_sets[n_NS - 1]->n_nodes;
	number_surfaces = db.surface_sets[n_SS - 1]->n_surf;
	Alloc();	//Alocação dinâmica
	int temp_surf = 0;
	for (int j = 0; j < number_surfaces; j++)
	{
		temp_surf = db.surface_sets[n_SS - 1]->surf_list[j];
		DOFs_surfaces[j] = db.surfaces[temp_surf - 1]->nDOFs;
	}
}

//checagem inicial do contato  - início de cada incremento
void NSSS::BeginStepCheck()
{

}

//Checks proximity between each node and surface
void NSSS::PinballCheck()
{
#pragma omp parallel
	{
	//Searching	
#pragma omp for
		for (int surf_index = 0; surf_index < number_surfaces; surf_index++)
		{
			Matrix distance(3);
			Matrix center(3);
			int tempsurf = db.surface_sets[n_SS - 1]->surf_list[surf_index];
			db.surfaces[tempsurf - 1]->CenterPoint(&center);
			for (int node_index = 0; node_index < number_nodes; node_index++)
			{
				int tempnode = db.node_sets[n_NS - 1]->node_list[node_index];
				//Check for proximity - calcular "distance"
				distance(0, 0) = center(0, 0) - db.nodes[tempnode - 1]->copy_coordinates[0];
				distance(1, 0) = center(1, 0) - db.nodes[tempnode - 1]->copy_coordinates[1];
				distance(2, 0) = center(2, 0) - db.nodes[tempnode - 1]->copy_coordinates[2];
				if (norm(distance) <= pinball)
				{
					activate[node_index][surf_index] = true;	//Near to contact
					cd[node_index][surf_index]->Alloc();		//Aloca - caso não haja pré-alocação, atribui valores iniciais a todas as variáveis
				}
					
				else
				{
					activate[node_index][surf_index] = false;				//Far to contact
					for (int sol_index = 0; sol_index < number_pointwise_interactions;sol_index++)
						FreeSpecific(node_index, surf_index,sol_index);		//Libera memória RAM
					cd[node_index][surf_index]->Free();						//Desaloca
				}
			}
		}
	}
}
//Salva variáveis para descrição lagrangiana atualizada
void NSSS::SaveLagrange()
{
	int temp_surf = 0;
	for (int surf_index = 0; surf_index < number_surfaces; surf_index++)
	{
		for (int node_index = 0; node_index < number_nodes; node_index++)
		{
			if (activate[node_index][surf_index] == true)
			{
				for (int sol_index = 0; sol_index < number_pointwise_interactions; sol_index++)
				{
					if (cd[node_index][surf_index]->g_n[sol_index] >= 0.0 || cd[node_index][surf_index]->copy_g_n[sol_index] >= 0.0 )
					{
						zeros(cd[node_index][surf_index]->g_t[sol_index]);
						zeros(cd[node_index][surf_index]->copy_g_t[sol_index]);
					}
					else
						*cd[node_index][surf_index]->copy_g_t[sol_index] = *cd[node_index][surf_index]->g_t[sol_index];
					cd[node_index][surf_index]->copy_g_n[sol_index] = cd[node_index][surf_index]->g_n[sol_index];
					cd[node_index][surf_index]->copy_convective[sol_index][0] = cd[node_index][surf_index]->convective[sol_index][0];
					cd[node_index][surf_index]->copy_convective[sol_index][1] = cd[node_index][surf_index]->convective[sol_index][1];
				}
			}
		}
	}
}

//Pode ser inserido algum critério impeditivo de convergência aqui
bool NSSS::HaveErrors()
{
	for (int surf_index = 0; surf_index < number_surfaces; surf_index++)
	{
		for (int node_index = 0; node_index < number_nodes; node_index++)
		{
			if (activate[node_index][surf_index] == true)
			{
				for (int sol = 0; sol < number_pointwise_interactions; sol++)
				{
					if (cd[node_index][surf_index]->repeated[sol] == false)
					{
						//Se algum dos pares ativos apresentou divergência do método de otimização
						if (cd[node_index][surf_index]->return_value[sol] == 1)
						{
							//data.myprintf("Node %d\n", data.node_sets[n_NS - 1]->node_list[node_index]);
							return true;
						}
					}
				}
			}
		}
	}
	return false;
}

//Calcula matriz com gap tangencial
void NSSS::EvaluateTangentialGap(int node_index, int surf_index, int sol_index)
{
	Matrix gtdelta(3);
	Matrix QM(3, 3);
	int temp_node = db.node_sets[n_NS - 1]->node_list[node_index];
	//Cálculo do QS - rotação do slave node
	(*alphaS[node_index][surf_index])(0, 0) = db.nodes[temp_node - 1]->displacements[3];
	(*alphaS[node_index][surf_index])(1, 0) = db.nodes[temp_node - 1]->displacements[4];
	(*alphaS[node_index][surf_index])(2, 0) = db.nodes[temp_node - 1]->displacements[5];
	double alpha_escalar = norm(*alphaS[node_index][surf_index]);
	double g = 4.0 / (4.0 + alpha_escalar*alpha_escalar);
	*QS[node_index][surf_index] = *I3 + g*(skew(*alphaS[node_index][surf_index]) + 0.5*(skew(*alphaS[node_index][surf_index])*skew(*alphaS[node_index][surf_index])));
	(*uS[node_index][surf_index])(0, 0) = db.nodes[temp_node - 1]->displacements[0];
	(*uS[node_index][surf_index])(1, 0) = db.nodes[temp_node - 1]->displacements[1];
	(*uS[node_index][surf_index])(2, 0) = db.nodes[temp_node - 1]->displacements[2];
	gtdelta = *uS[node_index][surf_index] + radius*(*I3 - *QS[node_index][surf_index])*(*cd[node_index][surf_index]->n_p[sol_index]) + *cd[node_index][surf_index]->G_i[sol_index] - *cd[node_index][surf_index]->G_ip[sol_index];
	gtdelta = gtdelta - dot(gtdelta, *cd[node_index][surf_index]->n_p[sol_index])*(*cd[node_index][surf_index]->n_p[sol_index]);
	//Rotação da normal de contato para convecção do gap tangencial acumulado dos passos anteriores
	Matrix to_i = -1.0*cross(*cd[node_index][surf_index]->t1_i[sol_index], *cd[node_index][surf_index]->n_i[sol_index]);
	Matrix to_p = -1.0*cross(*cd[node_index][surf_index]->t1_p[sol_index], *cd[node_index][surf_index]->n_p[sol_index]);
	QM(0, 0) = dot(*cd[node_index][surf_index]->t1_p[sol_index], *cd[node_index][surf_index]->t1_i[sol_index]);
	QM(0, 1) = dot(to_p, *cd[node_index][surf_index]->t1_i[sol_index]);
	QM(0, 2) = dot(*cd[node_index][surf_index]->n_p[sol_index], *cd[node_index][surf_index]->t1_i[sol_index]);
	QM(1, 0) = dot(*cd[node_index][surf_index]->t1_p[sol_index], to_i);
	QM(1, 1) = dot(to_p, to_i);
	QM(1, 2) = dot(*cd[node_index][surf_index]->n_p[sol_index], to_i);
	QM(2, 0) = dot(*cd[node_index][surf_index]->t1_p[sol_index], *cd[node_index][surf_index]->n_i[sol_index]);
	QM(2, 1) = dot(to_p, *cd[node_index][surf_index]->n_i[sol_index]);
	QM(2, 2) = dot(*cd[node_index][surf_index]->n_p[sol_index], *cd[node_index][surf_index]->n_i[sol_index]);
	Matrix Q(3, 3);
	Q(0, 0) = (*cd[node_index][surf_index]->t1_i[sol_index])(0, 0);
	Q(1, 0) = (*cd[node_index][surf_index]->t1_i[sol_index])(1, 0);
	Q(2, 0) = (*cd[node_index][surf_index]->t1_i[sol_index])(2, 0);
	Q(0, 1) = (to_i)(0, 0);
	Q(1, 1) = (to_i)(1, 0);
	Q(2, 1) = (to_i)(2, 0);
	Q(0, 2) = (*cd[node_index][surf_index]->n_i[sol_index])(0, 0);
	Q(1, 2) = (*cd[node_index][surf_index]->n_i[sol_index])(1, 0);
	Q(2, 2) = (*cd[node_index][surf_index]->n_i[sol_index])(2, 0);
	*cd[node_index][surf_index]->g_t[sol_index] = gtdelta + (Q*QM*transp(Q))*(*cd[node_index][surf_index]->copy_g_t[sol_index]);
}

//Imprime na tela informações do par de contato
void NSSS::ReportContact(int node_index, int surf_index, int sol_index)
{
	db.myprintf("Contact number %d\n", number);
	int tempnode = db.node_sets[n_NS - 1]->node_list[node_index];
	int tempsurf = db.surface_sets[n_SS - 1]->surf_list[surf_index];
	db.myprintf("Node %d Surface %d Solution %d\n", tempnode, tempsurf, sol_index + 1);
	db.myprintf("Return Value: %d\n", cd[node_index][surf_index]->return_value[sol_index]);
	db.myprintf("Convective coordinates: %.6f\t%.6f\n", cd[node_index][surf_index]->convective[sol_index][0], cd[node_index][surf_index]->convective[sol_index][1]);
	db.myprintf("Copy Convective coordinates: %.6f\t%.6f\n", cd[node_index][surf_index]->copy_convective[sol_index][0], cd[node_index][surf_index]->copy_convective[sol_index][1]);
	db.myprintf("Normal Gap: %.6f\n", cd[node_index][surf_index]->g_n[sol_index]);
	db.myprintf("Copy Normal Gap: %.6f\n", cd[node_index][surf_index]->copy_g_n[sol_index]);
	db.myprintf("Tangential Gap: %.6f\t%.6f\t%.6f\n", (*cd[node_index][surf_index]->g_t[sol_index])(0, 0), (*cd[node_index][surf_index]->g_t[sol_index])(1, 0), (*cd[node_index][surf_index]->g_t[sol_index])(2, 0));
	db.myprintf("Copy Tangential Gap: %.6f\t%.6f\t%.6f\n", (*cd[node_index][surf_index]->copy_g_t[sol_index])(0, 0), (*cd[node_index][surf_index]->copy_g_t[sol_index])(1, 0), (*cd[node_index][surf_index]->copy_g_t[sol_index])(2, 0));
	db.myprintf("Friction force: %.6f\n", norm(*cd[node_index][surf_index]->copy_g_t[sol_index])*ept);
}