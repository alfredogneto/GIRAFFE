#include "SSSS.h"
#define PI 3.1415926535897932384626433832795
#include"Database.h"
//Variáveis globais
extern
Database db;

SSSS::SSSS()
{
	type_name = new char[20];		//Nome do tipo do contato
	sprintf(type_name, "SSSS");
	number = 0;
	n_SS1 = 0;
	n_SS2 = 0;
	mus = 0.0;
	mud = 0.0;
	ept = 0.0;
	epn = 0.0;
	epn0 = 1.0;
	ct = 0.0;
	cn = 0.0;
	pinball = 0.0;
	number_pointwise_interactions = 1;	
	write_report = false;
	
	//Abaixo as variáveis que dependem do número de elementos para serem alocadas - Alocação feita na função Alloc, quando chamada durante o PreCalc
	typeOK1 = NULL;
	typeOK2 = NULL;
	activate = NULL;
	i_loading = NULL;
	contact_stiffness = NULL;
	alloc_control = NULL;
	DOFs_surfaces1 = NULL;
	DOFs_surfaces2 = NULL;
	cd = NULL;
	surf_pair = NULL;
	fn = NULL;
	ft = NULL;
	
	I3 = new Matrix(3, 3);
	(*I3)(0, 0) = 1.0;
	(*I3)(1, 1) = 1.0;
	(*I3)(2, 2) = 1.0;

	number_surfaces1 = 0;
	number_surfaces2 = 0;

	specialLCP = false;
}

void SSSS::Alloc()
{
	activate = new bool*[number_surfaces1];
	for (int i = 0; i < number_surfaces1; i++)
		activate[i] = new bool[number_surfaces2];
	
	alloc_control = new bool**[number_surfaces1];
	for (int i = 0; i < number_surfaces1; i++)
	{
		alloc_control[i] = new bool*[number_surfaces2];
		for (int j = 0; j < number_surfaces2; j++)
			alloc_control[i][j] = new bool[number_pointwise_interactions];
	}

	cd = new SSContactData**[number_surfaces1];
	for (int i = 0; i < number_surfaces1; i++)
	{
		cd[i] = new SSContactData*[number_surfaces2];
		for (int j = 0; j < number_surfaces2; j++)
		{
			cd[i][j] = new SSContactData();
			cd[i][j]->n_solutions = number_pointwise_interactions;
			//Alocação de cada um dos cd[i][j] será feita de acordo com a proximidade de contato - função Pinball
		}
	}

	DOFs_surfaces1 = new int[number_surfaces1];
	DOFs_surfaces2 = new int[number_surfaces2];

	//Alocação parcial de matrizes
	i_loading = new Matrix***[number_surfaces1];
	for (int i = 0; i < number_surfaces1; i++)
	{
		i_loading[i] = new Matrix**[number_surfaces2];
		for (int j = 0; j < number_surfaces2; j++)
			i_loading[i][j] = new Matrix*[number_pointwise_interactions];
	}
		
	contact_stiffness = new double****[number_surfaces1];
	for (int i = 0; i < number_surfaces1; i++)
	{
		contact_stiffness[i] = new double***[number_surfaces2];
		for (int j = 0; j < number_surfaces2; j++)
			contact_stiffness[i][j] = new double**[number_pointwise_interactions];
	}
	
	//Atribuição de valores iniciais
	for (int i = 0; i < number_surfaces1; i++)
	{
		for (int j = 0; j < number_surfaces2; j++)
		{
			activate[i][j] = false;
			for (int k = 0; k < number_pointwise_interactions; k++)
				alloc_control[i][j][k] = false;	//Não há alocação no início
		}
	}
	for (int j = 0; j < number_surfaces1; j++)
		DOFs_surfaces1[j] = 0;
	for (int j = 0; j < number_surfaces2; j++)
		DOFs_surfaces2[j] = 0;

	surf_pair = new SurfacePair**[number_surfaces1];

	for (int i = 0; i < number_surfaces1; i++)
	{
		surf_pair[i] = new SurfacePair*[number_surfaces2];
		for (int j = 0; j < number_surfaces2; j++)
			surf_pair[i][j] = NULL;
	}

	//Alocação de matrizes
	fn = new Matrix***[number_surfaces1];
	ft = new Matrix***[number_surfaces1];
	
	for (int i = 0; i < number_surfaces1; i++)
	{
		fn[i] = new Matrix**[number_surfaces2];
		ft[i] = new Matrix**[number_surfaces2];
	
		for (int j = 0; j < number_surfaces2; j++)
		{
			fn[i][j] = new Matrix*[number_pointwise_interactions];
			ft[i][j] = new Matrix*[number_pointwise_interactions];
			

			for (int k = 0; k < number_pointwise_interactions; k++)
			{
				fn[i][j][k] = new Matrix(3);
				ft[i][j][k] = new Matrix(3);
				
			}
		}	
	}
		
}

//Alocação específica para o contato surface1_index, surface2_index
void SSSS::AllocSpecific(int surface1_index, int surface2_index, int sol_index)
{
	//Caso não esteja ainda alocado, realiza a alocação
	if (alloc_control[surface1_index][surface2_index][sol_index] == false)
	{
		i_loading[surface1_index][surface2_index][sol_index] = new Matrix(DOFs_surfaces1[surface1_index] + DOFs_surfaces2[surface2_index]);
		contact_stiffness[surface1_index][surface2_index][sol_index] = new double*[DOFs_surfaces1[surface1_index] + DOFs_surfaces2[surface2_index]];
		for (int i = 0; i < (DOFs_surfaces1[surface1_index] + DOFs_surfaces2[surface2_index]); i++)
			contact_stiffness[surface1_index][surface2_index][sol_index][i] = new double[DOFs_surfaces1[surface1_index] + DOFs_surfaces2[surface2_index]];
		for (int ni = 0; ni < (DOFs_surfaces1[surface1_index] + DOFs_surfaces2[surface2_index]); ni++)
		for (int nj = 0; nj < (DOFs_surfaces1[surface1_index] + DOFs_surfaces2[surface2_index]); nj++)
			contact_stiffness[surface1_index][surface2_index][sol_index][ni][nj] = 0.0;
		
		alloc_control[surface1_index][surface2_index][sol_index] = true;
	}
}

//Liberação específica para o surface1_index, surface2_index
void SSSS::FreeSpecific(int surface1_index, int surface2_index, int sol_index)
{
	//Caso esteja alocado, realiza a desalocação
	if (alloc_control[surface1_index][surface2_index][sol_index] == true)
	{
		delete i_loading[surface1_index][surface2_index][sol_index];
		for (int i = 0; i < (DOFs_surfaces1[surface1_index] + DOFs_surfaces2[surface2_index]); i++)
			delete[] contact_stiffness[surface1_index][surface2_index][sol_index][i];
		delete[] contact_stiffness[surface1_index][surface2_index][sol_index];
		alloc_control[surface1_index][surface2_index][sol_index] = false;
	}
}

SSSS::~SSSS()
{
	delete[] type_name;
	delete I3;

	if (activate != NULL)
	{
		for (int i = 0; i < number_surfaces1; i++)
			delete [] activate[i];
		delete[]activate;
	}
	for (int i = 0; i < number_surfaces1;i++)
	for (int j = 0; j < number_surfaces2;j++)
	for (int k = 0; k < number_pointwise_interactions;k++)
		FreeSpecific(i,j,k);

	if (i_loading != NULL)
	{
		for (int i = 0; i < number_surfaces1; i++)
		{
			for (int j = 0; j < number_surfaces2; j++)
				delete[]i_loading[i][j];
			delete[] i_loading[i];
		}
		delete[]i_loading;
	}
	if (contact_stiffness != NULL)
	{
		for (int i = 0; i < number_surfaces1; i++)
		{
			for (int j = 0; j < number_surfaces2; j++)
				delete[] contact_stiffness[i][j];
			delete[] contact_stiffness[i];
		}
		delete[]contact_stiffness;
	}
	if (alloc_control != NULL)
	{
		for (int i = 0; i < number_surfaces1; i++)
		{
			for (int j = 0; j < number_surfaces2; j++)
				delete[]alloc_control[i][j];
			delete[] alloc_control[i];
		}
		delete[]alloc_control;
	}
	if (cd != NULL)
	{
		for (int i = 0; i < number_surfaces1; i++)
		{
			for (int j = 0; j < number_surfaces2; j++)
				delete cd[i][j];
			delete[]cd[i];
		}

		delete[]cd;
	}
	if (DOFs_surfaces1 != NULL)
		delete[]DOFs_surfaces1;
	if (DOFs_surfaces2 != NULL)
		delete[]DOFs_surfaces2;

	if (surf_pair != NULL)
	{
		for (int i = 0; i < number_surfaces1; i++)
		{
			for (int j = 0; j < number_surfaces2; j++)
			{
				if (surf_pair[i][j] != NULL)
					delete surf_pair[i][j];
			}
			delete[] surf_pair[i];
		}
		delete surf_pair;
	}

	if (fn != NULL)
	{
		for (int i = 0; i < number_surfaces1; i++)
		{
			for (int j = 0; j < number_surfaces2; j++)
			{
				for (int k = 0; k < number_pointwise_interactions; k++)
				{
					delete fn[i][j][k];
				}
				delete[]fn[i][j];
			}
			delete[] fn[i];
		}
		delete[]fn;
	}
	if (ft != NULL)
	{
		for (int i = 0; i < number_surfaces1; i++)
		{
			for (int j = 0; j < number_surfaces2; j++)
			{
				for (int k = 0; k < number_pointwise_interactions; k++)
				{
					delete ft[i][j][k];
				}
				delete[]ft[i][j];
			}
			delete[] ft[i];
		}
		delete[]ft;
	}
	
}

void SSSS::WriteVTK_XMLRender(FILE *f)
{
	//Plota superfícies
	int temp_surf1 = 0;
	for (int j = 0; j < number_surfaces1; j++)
	{
		temp_surf1 = db.surface_sets[n_SS1 - 1]->surf_list[j];
		db.surfaces[temp_surf1 - 1]->WriteVTK_XMLRender(f);
	}
	int temp_surf2 = 0;
	for (int j = 0; j < number_surfaces2; j++)
	{
		temp_surf2 = db.surface_sets[n_SS2 - 1]->surf_list[j];
		db.surfaces[temp_surf2 - 1]->WriteVTK_XMLRender(f);
	}
	
}

void SSSS::WriteVTK_XMLForces(FILE *f)
{
	//Plotagem de forças de contato
	if (db.post_files->WriteContactForces_flag == true)
	{
		//vetores para escrita no formato binário - usando a função 'enconde'
		std::vector<float> float_vector;
		std::vector<int> int_vector;
		//Verificação dos pontos a serem plotados
		int count = 0;
		for (int surf1_index = 0; surf1_index < number_surfaces1; surf1_index++)
		for (int surf2_index = 0; surf2_index < number_surfaces2; surf2_index++)
		for (int sol_index = 0; sol_index < number_pointwise_interactions; sol_index++)
			if (activate[surf1_index][surf2_index] == true && surf_pair[surf1_index][surf2_index] != NULL && cd[surf1_index][surf2_index]->repeated[sol_index] == false && (cd[surf1_index][surf2_index]->return_value[sol_index] == 0 && cd[surf1_index][surf2_index]->g_n[sol_index] < 0.0))
				count=count+2;
		//Opens Piece
		fprintf(f, "\t\t<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", count, count);
		//Opens Points
		fprintf(f, "\t\t\t<Points>\n");
		//Opens DataArray
		fprintf(f, "\t\t\t\t<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"binary\">\n");
		float_vector.clear();
		//Preenchendo as coordenadas dos pontos
		for (int surf1_index = 0; surf1_index < number_surfaces1; surf1_index++)
		for (int surf2_index = 0; surf2_index < number_surfaces2; surf2_index++)
		for (int sol_index = 0; sol_index < number_pointwise_interactions; sol_index++)
		if (activate[surf1_index][surf2_index] == true && surf_pair[surf1_index][surf2_index] != NULL && cd[surf1_index][surf2_index]->repeated[sol_index] == false && (cd[surf1_index][surf2_index]->return_value[sol_index] == 0 && cd[surf1_index][surf2_index]->g_n[sol_index] < 0.0))
		{
			float_vector.push_back((float)(*cd[surf1_index][surf2_index]->surf1->G_p[sol_index])(0, 0));
			float_vector.push_back((float)(*cd[surf1_index][surf2_index]->surf1->G_p[sol_index])(1, 0));
			float_vector.push_back((float)(*cd[surf1_index][surf2_index]->surf1->G_p[sol_index])(2, 0));

			float_vector.push_back((float)(*cd[surf1_index][surf2_index]->surf2->G_p[sol_index])(0, 0));
			float_vector.push_back((float)(*cd[surf1_index][surf2_index]->surf2->G_p[sol_index])(1, 0));
			float_vector.push_back((float)(*cd[surf1_index][surf2_index]->surf2->G_p[sol_index])(2, 0));
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
		fprintf(f, "\t\t\t\t<DataArray Name = \"Normal\" type = \"Float32\" NumberOfComponents= \"3\" format=\"binary\">\n");
		for (int surf1_index = 0; surf1_index < number_surfaces1; surf1_index++)
		for (int surf2_index = 0; surf2_index < number_surfaces2; surf2_index++)
		for (int sol_index = 0; sol_index < number_pointwise_interactions; sol_index++)
		if (activate[surf1_index][surf2_index] == true && surf_pair[surf1_index][surf2_index] != NULL && cd[surf1_index][surf2_index]->repeated[sol_index] == false && (cd[surf1_index][surf2_index]->return_value[sol_index] == 0 && cd[surf1_index][surf2_index]->g_n[sol_index] < 0.0))
		{
			float_vector.push_back(-(float)(*fn[surf1_index][surf2_index][sol_index])(0, 0));
			float_vector.push_back(-(float)(*fn[surf1_index][surf2_index][sol_index])(1, 0));
			float_vector.push_back(-(float)(*fn[surf1_index][surf2_index][sol_index])(2, 0));

			float_vector.push_back(+(float)(*fn[surf1_index][surf2_index][sol_index])(0, 0));
			float_vector.push_back(+(float)(*fn[surf1_index][surf2_index][sol_index])(1, 0));
			float_vector.push_back(+(float)(*fn[surf1_index][surf2_index][sol_index])(2, 0));
		}
		fprintf(f, encodeData(float_vector).c_str());
		fprintf(f, "\n");
		//Closes DataArray
		fprintf(f, "\t\t\t\t</DataArray>\n");

		float_vector.clear();
		fprintf(f, "\t\t\t\t<DataArray Name = \"Friction\" type = \"Float32\" NumberOfComponents= \"3\" format=\"binary\">\n");
		for (int surf1_index = 0; surf1_index < number_surfaces1; surf1_index++)
			for (int surf2_index = 0; surf2_index < number_surfaces2; surf2_index++)
				for (int sol_index = 0; sol_index < number_pointwise_interactions; sol_index++)
					if (activate[surf1_index][surf2_index] == true && surf_pair[surf1_index][surf2_index] != NULL && cd[surf1_index][surf2_index]->repeated[sol_index] == false && (cd[surf1_index][surf2_index]->return_value[sol_index] == 0 && cd[surf1_index][surf2_index]->g_n[sol_index] < 0.0))
					{
						float_vector.push_back(-(float)(*ft[surf1_index][surf2_index][sol_index])(0, 0));
						float_vector.push_back(-(float)(*ft[surf1_index][surf2_index][sol_index])(1, 0));
						float_vector.push_back(-(float)(*ft[surf1_index][surf2_index][sol_index])(2, 0));

						float_vector.push_back(+(float)(*ft[surf1_index][surf2_index][sol_index])(0, 0));
						float_vector.push_back(+(float)(*ft[surf1_index][surf2_index][sol_index])(1, 0));
						float_vector.push_back(+(float)(*ft[surf1_index][surf2_index][sol_index])(2, 0));
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

bool SSSS::Read(FILE *f)
{
	char s[1000];
	fscanf(f, "%s", s);
	number = atoi(s);
	fscanf(f, "%s", s);
	if (!strcmp(s, "SurfaceSet1"))
	{
		fscanf(f, "%s", s);
		n_SS1 = atoi(s);
	}
	else
		return false;
	fscanf(f, "%s", s);
	if (!strcmp(s, "SurfaceSet2"))
	{
		fscanf(f, "%s", s);
		n_SS2 = atoi(s);
	}
	else
		return false;

	//Salva a posição (stream)
	fpos_t pos;
	fgetpos(f, &pos);
	fscanf(f, "%s", s);
	if (!strcmp(s, "MU"))
	{
		fscanf(f, "%s", s);
		mus = atof(s);
		mud = atof(s);
	}
	else
	{
		fsetpos(f, &pos);
		fscanf(f, "%s", s);
		if (!strcmp(s, "MUS"))
		{
			fscanf(f, "%s", s);
			mus = atof(s);
		}
		else
			return false;
		fscanf(f, "%s", s);
		if (!strcmp(s, "MUD"))
		{
			fscanf(f, "%s", s);
			mud = atof(s);
		}
		else
			return false;
	}

	fscanf(f, "%s", s);
	if (!strcmp(s, "EPN"))
	{
		fscanf(f, "%s", s);
		epn = atof(s);
	}
	else
		return false;
	//Salva a posição (stream)
	fgetpos(f, &pos);
	fscanf(f, "%s", s);
	if (!strcmp(s, "EPN0"))
	{
		fscanf(f, "%s", s);
		epn0 = atof(s);
	}
	else
	{
		fsetpos(f, &pos);
		epn0 = 1.0;
	}
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

	//Salva a posição (stream)
	fgetpos(f, &pos);
	fscanf(f, "%s", s);
	if (!strcmp(s, "MaxPointwiseInt"))
	{
		fscanf(f, "%s", s);
		number_pointwise_interactions = atoi(s);
	}
	else
	{
		fsetpos(f, &pos);
		number_pointwise_interactions = 1;
	}

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

	//Salva a posição (stream)
	fgetpos(f, &pos);
	fscanf(f, "%s", s);
	if (!strcmp(s, "SpecialLCP"))
		specialLCP = true;
	else
	{
		fsetpos(f, &pos);
		specialLCP = false;
	}

	//Salva a posição (stream)
	fgetpos(f, &pos);
	fscanf(f, "%s", s);
	if (!strcmp(s, "WriteReport"))
		write_report = true;
	else
	{
		fsetpos(f, &pos);
		write_report = false;
	}

	//Salva a posição (stream)
	fgetpos(f, &pos);
	fscanf(f, "%s", s);
	if (!strcmp(s, "WriteReportDiverged"))
		write_report_diverged = true;
	else
	{
		fsetpos(f, &pos);
		write_report_diverged = false;
	}
		
	return true;
}

void SSSS::Write(FILE *f)
{
	fprintf(f, "SSSS\t%d\tSurfaceSet1\t%d\tSurfaceSet2\t%d\tMUS\t%.6e\tMUD\t%.6e\tEPN\t%.6e\tCN\t%.6e\tEPT\t%.6e\tCT\t%.6e\tPinball\t%.6e\tMaxPointwiseInt\t%d\t",
		number, n_SS1, n_SS2, mus, mud, epn, cn, ept, ct, pinball, number_pointwise_interactions);
	bool_table.Write(f);
}

//Checa inconsistências no elemento para evitar erros de execução
bool SSSS::Check()
{
	if (n_SS1 > db.number_surface_sets || n_SS2 > db.number_surface_sets)
	{
		printf("Error in Contact %d. Invalid Surface Set assigned. Check the input file.\n", number);
		return false;
	}
	return true;
}

//Escreve arquivo de resultados
void SSSS::WriteResults(FILE *f)
{
	//TODO
}

//Escreve no monitor do contato
void SSSS::WriteMonitor(FILE *f, bool first_record, double time)
{
	if (db.monitor->contact_special_output)
	{
		if (first_record == true)
		{
			fprintf(f, "Active contacts output only\nTIME\tSurf1ID\tTheta1\tSurf2ID\tTheta2\tFn\n");
		}
		for (int i = 0; i < number_surfaces1; i++)
		{
			int count = 0;
			Matrix Fn(3);
			Matrix Fat(3);
			for (int j = 0; j < number_surfaces2; j++)
			{
				 
				if (activate[i][j] == 1)
				{
					for (int k = 0; k < number_pointwise_interactions; k++)
					{
						if (cd[i][j]->return_value[k] == 0 && cd[i][j]->g_n[k] < 0.0 && cd[i][j]->repeated[k] == false)
							fprintf(f, "%.6e\t%d\t%.12e\t%d\t%.12e\t%.6e\n", time, db.surface_sets[n_SS1 - 1]->surf_list[i], cd[i][j]->convective[k][1], db.surface_sets[n_SS2 - 1]->surf_list[j], cd[i][j]->convective[k][3], norm(*fn[i][j][k]));
					}

				}
			}
		}
	}
	else
	{
		//Cabeçalho
		if (first_record == true)
		{
			fprintf(f, "\t");
			for (int i = 0; i < number_surfaces1; i++)
				fprintf(f, "Surface %d\t\t\t\t\t\t\t\t\t", db.surface_sets[n_SS1 - 1]->surf_list[i]);
			for (int i = 0; i < number_surfaces2; i++)
				fprintf(f, "Surface %d\t\t\t\t\t\t\t\t\t", db.surface_sets[n_SS2 - 1]->surf_list[i]);
			fprintf(f, "\n");
			fprintf(f, "TIME");
			for (int i = 0; i < number_surfaces1; i++)
				fprintf(f, "\tFnx\tFny\tFnz\tFn\tFatx\tFaty\tFatz\tFat\tNForces");
			for (int i = 0; i < number_surfaces2; i++)
				fprintf(f, "\tFnx\tFny\tFnz\tFn\tFatx\tFaty\tFatz\tFat\tNForces");
			fprintf(f, "\n");
		}
		fprintf(f, "%.6e\t", time);	//TIME
		for (int i = 0; i < number_surfaces1; i++)
		{
			int count = 0;
			Matrix Fn(3);
			Matrix Fat(3);
			for (int j = 0; j < number_surfaces2; j++)
			{

				if (activate[i][j] == 1)
				{
					for (int k = 0; k < number_pointwise_interactions; k++)
					{
						if (cd[i][j]->return_value[k] == 0 && cd[i][j]->g_n[k] < 0.0 && cd[i][j]->repeated[k] == false)
						{
							count++;
							Fn = Fn - *fn[i][j][k];
							Fat = Fat - *ft[i][j][k];
						}
					}

				}
			}
			fprintf(f, "%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%d\t", Fn(0, 0), Fn(1, 0), Fn(2, 0), norm(Fn), Fat(0, 0), Fat(1, 0), Fat(2, 0), norm(Fat), count);
		}
		for (int i = 0; i < number_surfaces2; i++)
		{
			int count = 0;
			Matrix Fn(3);
			Matrix Fat(3);

			for (int j = 0; j < number_surfaces1; j++)
			{

				if (activate[j][i] == 1)
				{
					for (int k = 0; k < number_pointwise_interactions; k++)
					{
						if (cd[j][i]->return_value[k] == 0 && cd[j][i]->g_n[k] < 0.0 && cd[j][i]->repeated[k] == false)
						{
							count++;
							Fn = -1.0*Fn + *fn[j][i][k];
							Fat = -1.0*Fat + *ft[j][i][k];
						}
					}

				}
			}
			fprintf(f, "%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%d\t", Fn(0, 0), Fn(1, 0), Fn(2, 0), norm(Fn), Fat(0, 0), Fat(1, 0), Fat(2, 0), norm(Fat), count);
		}

		fprintf(f, "\n");
	}
}

void SSSS::Mount()
{
	//Loop in nodes to perform contact Mount for each one
#pragma omp parallel for
	for (int surf1_index = 0; surf1_index < number_surfaces1; surf1_index++)
	{
		int const tempsurf1 = db.surface_sets[n_SS1 - 1]->surf_list[surf1_index];
		for (int surf2_index = 0; surf2_index < number_surfaces2; surf2_index++)
		{
			int const tempsurf2 = db.surface_sets[n_SS2 - 1]->surf_list[surf2_index];
			if (activate[surf1_index][surf2_index] == true && surf_pair[surf1_index][surf2_index] != NULL)
			{	
				//Solução do problema de mínima distância
				surf_pair[surf1_index][surf2_index]->SolveLCP(cd[surf1_index][surf2_index]);
				//Varredura nos possíveis pontos de interação pointwise
				for (int sol_index = 0; sol_index < number_pointwise_interactions; sol_index++)
				{
					//Se não for uma raiz repetida (já computada anteriormente)
					if (cd[surf1_index][surf2_index]->repeated[sol_index] == false)
					{
						if (cd[surf1_index][surf2_index]->return_value[sol_index] == 0 || cd[surf1_index][surf2_index]->return_value[sol_index] == 4)
							EvaluateNormalGap(surf1_index, surf2_index, sol_index);
						else
							cd[surf1_index][surf2_index]->g_n[sol_index] = 1.0;	//valor para indicar que não há contato
						/*if (cd[surf1_index][surf2_index]->return_value[sol_index] == 0)
							ReportContact(surf1_index, surf2_index, sol_index);*/
						//Se houver penetração e o range de coordenadas convectivas for válido monta a contribuição do contato na forma fraca e operador tangente
						if (cd[surf1_index][surf2_index]->g_n[sol_index] < -pen_tol && (cd[surf1_index][surf2_index]->return_value[sol_index] == 0 || cd[surf1_index][surf2_index]->return_value[sol_index] == 4))
						{
							//Aloca espaço para salvar a matriz de rigidez e vetor de esforços
							AllocSpecific(surf1_index, surf2_index, sol_index);
							bool previouscontact = false;
							if (cd[surf1_index][surf2_index]->copy_g_n[sol_index] < -pen_tol)
								previouscontact = true;
							//Avalia contribuições de contato
							surf_pair[surf1_index][surf2_index]->ContactSS(&cd[surf1_index][surf2_index]->copy_stick[sol_index], &cd[surf1_index][surf2_index]->stick[sol_index], &previouscontact,
								i_loading[surf1_index][surf2_index][sol_index]->getMatrix(), contact_stiffness[surf1_index][surf2_index][sol_index], 
								cd[surf1_index][surf2_index]->invHessian[sol_index], cd[surf1_index][surf2_index]->convective[sol_index], 
								cd[surf1_index][surf2_index]->copy_convective[sol_index], cd[surf1_index][surf2_index]->copy_g_t[sol_index]->getMatrix(), 
								cd[surf1_index][surf2_index]->g_t[sol_index]->getMatrix(),
								&epn,&epn0, &ept, &cn, &ct, &mus, &mud,
								fn[surf1_index][surf2_index][sol_index]->getMatrix(), ft[surf1_index][surf2_index][sol_index]->getMatrix());
						}
						else
							FreeSpecific(surf1_index, surf2_index, sol_index);
					}
				}
			}
		}
	}
}

//Montagens - Newmark
void SSSS::MountDyn()
{
	//NOT USED (dynamic effects added inside the function Mount())
}
//Preenche a contribuição do contato nas matrizes globais
void SSSS::MountGlobal()
{
	//Variáveis temporárias para salvar a indexação global dos graus de liberdade a serem setados na matriz de rigidez global
	int GL_global_1 = 0;
	int GL_global_2 = 0;
	double anterior = 0;
	int tempsurf1, tempsurf2;
	for (int surf1_index = 0; surf1_index < number_surfaces1; surf1_index++)
	{
		//Superfície 1
		tempsurf1 = db.surface_sets[n_SS1 - 1]->surf_list[surf1_index];
		for (int surf2_index = 0; surf2_index < number_surfaces2; surf2_index++)
		{
			if (activate[surf1_index][surf2_index] == true && surf_pair[surf1_index][surf2_index] != NULL)
			{
				//Superfície 2
				tempsurf2 = db.surface_sets[n_SS2 - 1]->surf_list[surf2_index];
				for (int sol_index = 0; sol_index < number_pointwise_interactions; sol_index++)
				{
					//Se não for uma raiz repetida (já computada anteriormente)
					if (cd[surf1_index][surf2_index]->repeated[sol_index] == false)
					{
						/*if (cd[surf1_index][surf2_index]->return_value[sol_index] == 0 || cd[surf1_index][surf2_index]->return_value[sol_index] == 4)
							ReportContact(surf1_index, surf2_index, sol_index);*/
						if (cd[surf1_index][surf2_index]->return_value[sol_index] == 0 && cd[surf1_index][surf2_index]->g_n[sol_index] < -pen_tol)
						{
							//ReportContact(surf1_index, surf2_index, sol_index);
							//i_loading[surf1_index][surf2_index][sol_index]->print();
							//PrintPtr(contact_stiffness[surf1_index][surf2_index][sol_index], 12, 12);
							for (int i = 0; i < (DOFs_surfaces1[surf1_index] + DOFs_surfaces2[surf2_index]); i++)
							{
								
								if (i < DOFs_surfaces1[surf1_index])
									GL_global_1 = *db.surfaces[tempsurf1 - 1]->GLs[i];									//GLs da superfície 1
								else
									GL_global_1 = *db.surfaces[tempsurf2 - 1]->GLs[i - DOFs_surfaces1[surf1_index]];		//GLs da superfície 2
								//Caso o grau de liberdade seja livre:
								if (GL_global_1 > 0)
								{
									anterior = db.global_P_A(GL_global_1 - 1, 0);
									db.global_P_A(GL_global_1 - 1, 0) = anterior + (*i_loading[surf1_index][surf2_index][sol_index])(i, 0);
									anterior = db.global_I_A(GL_global_1 - 1, 0);
									db.global_I_A(GL_global_1 - 1, 0) = anterior + (*i_loading[surf1_index][surf2_index][sol_index])(i, 0);
								}
								else
								{
									if (GL_global_1 != 0)//se o GL é ativo
									{
										anterior = db.global_P_B(-GL_global_1 - 1, 0);
										db.global_P_B(-GL_global_1 - 1, 0) = anterior + (*i_loading[surf1_index][surf2_index][sol_index])(i, 0);
									}

								}
								for (int j = 0; j < (DOFs_surfaces1[surf1_index] + DOFs_surfaces2[surf2_index]); j++)
								{
									if (j < DOFs_surfaces1[surf1_index])
										GL_global_2 = *db.surfaces[tempsurf1 - 1]->GLs[j];									//GLs da superfície 1
									else
										GL_global_2 = *db.surfaces[tempsurf2 - 1]->GLs[j - DOFs_surfaces1[surf1_index]];		//GLs da superfície 2

									//Caso os graus de liberdade sejam ambos livres (Matriz Kaa)
									if (GL_global_1 > 0 && GL_global_2 > 0)
										db.global_stiffness_AA.setValue(GL_global_1 - 1, GL_global_2 - 1, contact_stiffness[surf1_index][surf2_index][sol_index][i][j]);
									//Caso os graus de liberdade sejam ambos fixos (Matriz Kbb)
									if (GL_global_1 < 0 && GL_global_2 < 0)
										db.global_stiffness_BB.setValue(-GL_global_1 - 1, -GL_global_2 - 1, contact_stiffness[surf1_index][surf2_index][sol_index][i][j]);
									//Caso os graus de liberdade sejam livre e fixo (Matriz Kab)
									if (GL_global_1 > 0 && GL_global_2 < 0)
										db.global_stiffness_AB.setValue(GL_global_1 - 1, -GL_global_2 - 1, contact_stiffness[surf1_index][surf2_index][sol_index][i][j]);
									//Caso os graus de liberdade sejam fixo e livre (Matriz Kba)
									if (GL_global_1 < 0 && GL_global_2 > 0)
										db.global_stiffness_BA.setValue(-GL_global_1 - 1, GL_global_2 - 1, contact_stiffness[surf1_index][surf2_index][sol_index][i][j]);
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
void SSSS::Band(int* band_fixed, int* band_free)
{
	//NOT IMPLEMENTED
}

//Pré-cálculo de variáveis que é feito uma única vez no início
void SSSS::PreCalc()
{
	number_surfaces1 = db.surface_sets[n_SS1 - 1]->n_surf;
	number_surfaces2 = db.surface_sets[n_SS2 - 1]->n_surf;
	int temp_surf = 0;
	
	//Neste ponto já foi chamado o PreCalc() das Surfaces. Com isso já há informação de degeneração calculada em cada uma delas.

	//Setando os número de interações pontuais por par de superfície - de acordo com o número de degenerações presentes
	int n1, n2;
	for (int i = 0; i < number_surfaces1; i++)
	{
		temp_surf = db.surface_sets[n_SS1 - 1]->surf_list[i];
		n1 = db.surfaces[temp_surf - 1]->div1*db.surfaces[temp_surf - 1]->div2;
		for (int j = 0; j < number_surfaces2; j++)
		{
			temp_surf = db.surface_sets[n_SS2 - 1]->surf_list[j];
			n2 = db.surfaces[temp_surf - 1]->div1*db.surfaces[temp_surf - 1]->div2;
			number_pointwise_interactions = n1 * n2;
		}
	}
	
	Alloc();	//Alocação dinâmica
	
	for (int j = 0; j < number_surfaces1; j++)
	{
		temp_surf = db.surface_sets[n_SS1 - 1]->surf_list[j];
		DOFs_surfaces1[j] = db.surfaces[temp_surf - 1]->nDOFs;
	}
	for (int j = 0; j < number_surfaces2; j++)
	{
		temp_surf = db.surface_sets[n_SS2 - 1]->surf_list[j];
		DOFs_surfaces2[j] = db.surfaces[temp_surf - 1]->nDOFs;
	}
	//Atribui e aloca pares de superfícies
	SetPairs();
}

//checagem inicial do contato  - início de cada incremento
void SSSS::BeginStepCheck()
{
	for (int surf1_index = 0; surf1_index < number_surfaces1; surf1_index++)
	for (int surf2_index = 0; surf2_index < number_surfaces2; surf2_index++)
	{
		//Se o contato for ativo - realiza estudo de degeneração do início do step
		if (activate[surf1_index][surf2_index] == true)
			surf_pair[surf1_index][surf2_index]->BeginStepCheck(cd[surf1_index][surf2_index]);
	}
		
}

//Checks proximity between each node and surface
void SSSS::PinballCheck()
{
#pragma omp parallel
	{
	//Searching	
#pragma omp for
		for (int surf1_index = 0; surf1_index < number_surfaces1; surf1_index++)
		{
			Matrix center1(3);
			Matrix center2(3);
			int tempsurf1 = db.surface_sets[n_SS1 - 1]->surf_list[surf1_index];
			db.surfaces[tempsurf1 - 1]->CenterPoint(&center1);
			for (int surf2_index = 0; surf2_index < number_surfaces2; surf2_index++)
			{
				int tempsurf2 = db.surface_sets[n_SS2 - 1]->surf_list[surf2_index];
				//Atualiza informações com os últimos deslocamentos ocorridos - prepara parametrização da superfície para montagem do contato
				db.surfaces[tempsurf2 - 1]->CenterPoint(&center2);
				//Verifica se é o mesmo ID de superfície
				if (tempsurf1 == tempsurf2)
					activate[surf1_index][surf2_index] = false;					//DESATIVA
				else
				{
					if (norm(center1 - center2) <= pinball)
					{
						activate[surf1_index][surf2_index] = true;					//Near to contact
						cd[surf1_index][surf2_index]->Alloc();						//Aloca - caso não haja pré-alocação, atribui valores iniciais a todas as variáveis
						surf_pair[surf1_index][surf2_index]->Alloc(cd[surf1_index][surf2_index]);
					}
					else
					{
						activate[surf1_index][surf2_index] = false;					//Far to contact
						for (int sol_index = 0; sol_index < number_pointwise_interactions; sol_index++)
							FreeSpecific(surf1_index, surf2_index, sol_index);		//Libera memória RAM
						cd[surf1_index][surf2_index]->Free();						//Desaloca
						surf_pair[surf1_index][surf2_index]->Free();
					}
				}
				
			}
		}
	}
}
//Salva variáveis para descrição lagrangiana atualizada
void SSSS::SaveLagrange()
{
	for (int surf1_index = 0; surf1_index < number_surfaces1; surf1_index++)
	{
		for (int surf2_index = 0; surf2_index < number_surfaces2; surf2_index++)
		{
			if (activate[surf1_index][surf2_index] == true && surf_pair[surf1_index][surf2_index] != NULL)
			{
				for (int sol_index = 0; sol_index < number_pointwise_interactions; sol_index++)
				{
					//Zera o gap tangencial se:
					//1 - gn for positivo (não há contato)
					//2 - copy_gn for positivo (é a primeira ocorrência de contato - não há como acumular gap tangencial pois o contato acaba de começar)
					//3 - return value é 2 - contato não estabelecido - e não é strong candidate, mas pode ocorrer em próximos instantes
					if (cd[surf1_index][surf2_index]->g_n[sol_index] >= pen_tol || cd[surf1_index][surf2_index]->copy_g_n[sol_index] >= pen_tol || cd[surf1_index][surf2_index]->return_value[sol_index] == 2)
					{
						zeros(cd[surf1_index][surf2_index]->g_t[sol_index]);
						zeros(cd[surf1_index][surf2_index]->copy_g_t[sol_index]);
					}
					else
						*cd[surf1_index][surf2_index]->copy_g_t[sol_index] = *cd[surf1_index][surf2_index]->g_t[sol_index];
					
					cd[surf1_index][surf2_index]->copy_g_n[sol_index] = cd[surf1_index][surf2_index]->g_n[sol_index];
					*cd[surf1_index][surf2_index]->copy_g[sol_index] = *cd[surf1_index][surf2_index]->g[sol_index];
					*cd[surf1_index][surf2_index]->copy_n[sol_index] = *cd[surf1_index][surf2_index]->n[sol_index];
					cd[surf1_index][surf2_index]->copy_return_value[sol_index] = cd[surf1_index][surf2_index]->return_value[sol_index];
					cd[surf1_index][surf2_index]->copy_convective[sol_index][0] = cd[surf1_index][surf2_index]->convective[sol_index][0];
					cd[surf1_index][surf2_index]->copy_convective[sol_index][1] = cd[surf1_index][surf2_index]->convective[sol_index][1];
					cd[surf1_index][surf2_index]->copy_convective[sol_index][2] = cd[surf1_index][surf2_index]->convective[sol_index][2];
					cd[surf1_index][surf2_index]->copy_convective[sol_index][3] = cd[surf1_index][surf2_index]->convective[sol_index][3];
					cd[surf1_index][surf2_index]->copy_degenerated[sol_index] = cd[surf1_index][surf2_index]->degenerated[sol_index];
					cd[surf1_index][surf2_index]->copy_stick[sol_index] = cd[surf1_index][surf2_index]->stick[sol_index];
				}
			}
		}
	}
}

//Pode ser inserido algum critério impeditivo de convergência aqui
bool SSSS::HaveErrors()
{
	for (int i = 0; i < number_surfaces1; i++)
	{
		for (int j = 0; j < number_surfaces2; j++)
		{
			if (activate[i][j] == true && surf_pair[i][j] != NULL)
			{
				//ReportContact(i, j, 0);
				//A função abaixo é encarregada de checar de acordo com a conveniência do par de superfícies as condições para cut-back
				if (surf_pair[i][j]->EndStepCheck(cd[i][j]) == true)
					return true;
			}
		}
	}
	return false;
}

//Calcula matriz com gap normal entre surface1, surface2
void SSSS::EvaluateNormalGap(int surf1_index, int surf2_index, int sol_index)
{
	int const tempsurf1 = db.surface_sets[n_SS1 - 1]->surf_list[surf1_index];
	int const tempsurf2 = db.surface_sets[n_SS2 - 1]->surf_list[surf2_index];
	Matrix temp_conv(4);
	for (int co = 0; co < 4; co++)
		temp_conv(co, 0) = cd[surf1_index][surf2_index]->convective[sol_index][co];
	//Cálculo da função gap (escalar)
	Matrix GammaA(3);
	Matrix GammaB(3);
	db.surfaces[tempsurf1 - 1]->SurfacePoint(temp_conv(0, 0), temp_conv(1, 0), GammaA);
	db.surfaces[tempsurf2 - 1]->SurfacePoint(temp_conv(2, 0), temp_conv(3, 0), GammaB);
	//Copies to surface data
	*cd[surf1_index][surf2_index]->surf1->G_p[sol_index] = GammaA;
	*cd[surf1_index][surf2_index]->surf2->G_p[sol_index] = GammaB;
	//Gap vetorial
	*cd[surf1_index][surf2_index]->g[sol_index] = GammaA - GammaB;
	//Normal do contato
	if (norm(*cd[surf1_index][surf2_index]->g[sol_index]) != 0.0)
		*cd[surf1_index][surf2_index]->n[sol_index] = (1.0 / norm(*cd[surf1_index][surf2_index]->g[sol_index]))*(*cd[surf1_index][surf2_index]->g[sol_index]);
	else
		zeros(cd[surf1_index][surf2_index]->n[sol_index]);

	//Gap escalar - tentativa de calculo com base na normal anterior (se a mesma existir)
	if (cd[surf1_index][surf2_index]->copy_return_value[sol_index] != 2)
	{
		bool previouscontact = false;
		bool contact = false;
		if (cd[surf1_index][surf2_index]->copy_g_n[sol_index] <= 0)
			previouscontact = true;
		if (dot(*cd[surf1_index][surf2_index]->n[sol_index], *cd[surf1_index][surf2_index]->copy_n[sol_index]) >= 0)
			contact = previouscontact;
		else
		{
			if (previouscontact)
				contact = false;
			else
				contact = true;
		}
		if (contact)
			cd[surf1_index][surf2_index]->g_n[sol_index] = -1.0*norm(*cd[surf1_index][surf2_index]->g[sol_index]);
		else
			cd[surf1_index][surf2_index]->g_n[sol_index] = +1.0*norm(*cd[surf1_index][surf2_index]->g[sol_index]);
	}
	//Gap escalar - calculo com base nas normais das superficies
	else
	{
		//Verificação de normais 
		bool passed_normal_test = false;
		Matrix n1(3);
		Matrix n2(3);
		db.surfaces[tempsurf1 - 1]->NormalExt(&temp_conv(0, 0), &temp_conv(1, 0), &n1);
		db.surfaces[tempsurf2 - 1]->NormalExt(&temp_conv(2, 0), &temp_conv(3, 0), &n2);
		if (dot(n1, n2) < 0)
			passed_normal_test = true;
		if (surf_pair[surf1_index][surf2_index]->Gap(temp_conv, false, n1, n2) <= 0.0  && passed_normal_test)
			cd[surf1_index][surf2_index]->g_n[sol_index] = -1.0*norm(*cd[surf1_index][surf2_index]->g[sol_index]);
		else
			cd[surf1_index][surf2_index]->g_n[sol_index] = +1.0*norm(*cd[surf1_index][surf2_index]->g[sol_index]);
	}
}

//Calcula matriz com gap tangencial
void SSSS::EvaluateTangentialGap(int surf1_index, int surf2_index, int sol_index)
{
	int const tempsurf1 = db.surface_sets[n_SS1 - 1]->surf_list[surf1_index];
	int const tempsurf2 = db.surface_sets[n_SS2 - 1]->surf_list[surf2_index];
	//Previous convective coordinates
	Matrix copy_conv(4);
	for (int co = 0; co < 4; co++)
		copy_conv(co, 0) = cd[surf1_index][surf2_index]->copy_convective[sol_index][co];
	Matrix nA(3);
	Matrix nB(3);
	Matrix GammaA(3);
	Matrix GammaB(3);
	db.surfaces[tempsurf1 - 1]->SurfacePoint(copy_conv(0, 0), copy_conv(1, 0), GammaA);
	db.surfaces[tempsurf2 - 1]->SurfacePoint(copy_conv(2, 0), copy_conv(3, 0), GammaB);
	Matrix nip = *cd[surf1_index][surf2_index]->n[sol_index];
	Matrix ni = *cd[surf1_index][surf2_index]->copy_n[sol_index];
	Matrix gtdelta(3);
	gtdelta = (*I3 - dyadic(nip, nip))*(GammaA - GammaB);
	
	//Rotação da normal de contato para rotacao do gap tangencial acumulado dos passos anteriores
	Matrix Qdelta(3, 3);
	Matrix esin = cross(ni, nip);
	double norm_esin = norm(esin);
	double TOL = 1e-8;
	if (norm_esin > TOL)
	{
		Matrix e = (1.0 / norm_esin)*esin;
		double theta_scalar = asin(norm_esin);
		double alpha_scalar = 2.0*tan(theta_scalar / 2);
		Matrix alpha = alpha_scalar * e;
		Matrix A = skew(alpha);
		double g = 4.0 / (4.0 + alpha_scalar * alpha_scalar);
		Qdelta = *I3 + g * (A + 0.5*A*A);			//Tensor de rotação
	}
	else
		Qdelta = *I3;
	
	*cd[surf1_index][surf2_index]->g_t[sol_index] = gtdelta + Qdelta *(*cd[surf1_index][surf2_index]->copy_g_t[sol_index]);
}

//Imprime na tela informações do par de contato
void SSSS::ReportContact(int surf1_index, int surf2_index, int sol_index)
{
	db.myprintf("Contact number %d\n", number);
	int tempsurf1 = db.surface_sets[n_SS1 - 1]->surf_list[surf1_index];
	int tempsurf2 = db.surface_sets[n_SS2 - 1]->surf_list[surf2_index];
	db.myprintf("Surface %d Surface %d Solution %d\n", tempsurf1, tempsurf2, sol_index);
	db.myprintf("Return Value: %d\n", cd[surf1_index][surf2_index]->return_value[sol_index]);
	db.myprintf("Copy Return Value: %d\n", cd[surf1_index][surf2_index]->copy_return_value[sol_index]);
	db.myprintf("Convective coordinates: %.6e\t%.6e\t%.6e\t%.6e\n", cd[surf1_index][surf2_index]->convective[sol_index][0], cd[surf1_index][surf2_index]->convective[sol_index][1], cd[surf1_index][surf2_index]->convective[sol_index][2], cd[surf1_index][surf2_index]->convective[sol_index][3]);
	db.myprintf("Copy Convective coordinates: %.6e\t%.6e\t%.6e\t%.6e\n", cd[surf1_index][surf2_index]->copy_convective[sol_index][0], cd[surf1_index][surf2_index]->copy_convective[sol_index][1], cd[surf1_index][surf2_index]->copy_convective[sol_index][2], cd[surf1_index][surf2_index]->copy_convective[sol_index][3]);
	db.myprintf("Normal Gap: %.6f\n", cd[surf1_index][surf2_index]->g_n[sol_index]);
	
	db.myprintf("Copy Normal Gap: %.6e\n", cd[surf1_index][surf2_index]->copy_g_n[sol_index]);
	db.myprintf("Tangential Gap: %.6e\t%.6e\t%.6e\n", (*cd[surf1_index][surf2_index]->g_t[sol_index])(0, 0), (*cd[surf1_index][surf2_index]->g_t[sol_index])(1, 0), (*cd[surf1_index][surf2_index]->g_t[sol_index])(2, 0));
	db.myprintf("Vetor Gap: %.6e\t%.6e\t%.6e\n", (*cd[surf1_index][surf2_index]->g[sol_index])(0, 0), (*cd[surf1_index][surf2_index]->g[sol_index])(1, 0), (*cd[surf1_index][surf2_index]->g[sol_index])(2, 0));
	db.myprintf("Copy Tangential Gap: %.6e\t%.6e\t%.6e\n", (*cd[surf1_index][surf2_index]->copy_g_t[sol_index])(0, 0), (*cd[surf1_index][surf2_index]->copy_g_t[sol_index])(1, 0), (*cd[surf1_index][surf2_index]->copy_g_t[sol_index])(2, 0));
	db.myprintf("Friction force: %.6e\n", norm(*cd[surf1_index][surf2_index]->copy_g_t[sol_index])*ept);
}

//Atribui pares de superfícies, de acordo com varredura feita entre os surface sets 1 e 2
void SSSS::SetPairs()
{ 
	int temp_ID_1, temp_ID_2;
	bool control = false;
	for (int i = 0; i < number_surfaces1; i++)
	{
		for (int j = 0; j < number_surfaces2; j++)
		{
			control = false;
			temp_ID_1 = db.surface_sets[n_SS1 - 1]->surf_list[i];
			temp_ID_2 = db.surface_sets[n_SS2 - 1]->surf_list[j];
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			/////////////////////////////////////////////////////////////FlexibleSECylinder_1_FlexibleSECylinder_1//////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			if (typeid(*db.surfaces[temp_ID_1 - 1]) == typeid(FlexibleSECylinder_1) && typeid(*db.surfaces[temp_ID_2 - 1]) == typeid(FlexibleSECylinder_1))
			{
				surf_pair[i][j] = new FlexibleSECylinder_1_FlexibleSECylinder_1();
				surf_pair[i][j]->inverted = false;
				surf_pair[i][j]->surf1_ID = temp_ID_1;
				surf_pair[i][j]->surf2_ID = temp_ID_2;
				
				control = true;
			}
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			/////////////////////////////////////////////////////////////FlexibleArcExtrusion_1_RigidArcRevolution_1//////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			if (typeid(*db.surfaces[temp_ID_1 - 1]) == typeid(FlexibleArcExtrusion_1) && typeid(*db.surfaces[temp_ID_2 - 1]) == typeid(RigidArcRevolution_1))
			{
				surf_pair[i][j] = new FlexibleArcExtrusion_1_RigidArcRevolution_1();
				surf_pair[i][j]->inverted = false;
				surf_pair[i][j]->surf1_ID = temp_ID_1;
				surf_pair[i][j]->surf2_ID = temp_ID_2;

				control = true;
			}
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			/////////////////////////////////////////////////////////////RigidNURBS_1_RigidNURBS_1//////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			if (typeid(*db.surfaces[temp_ID_1 - 1]) == typeid(RigidNURBS_1) && typeid(*db.surfaces[temp_ID_2 - 1]) == typeid(RigidNURBS_1))
			{
				surf_pair[i][j] = new RigidNURBS_1_RigidNURBS_1();
				surf_pair[i][j]->inverted = false;
				surf_pair[i][j]->surf1_ID = temp_ID_1;
				surf_pair[i][j]->surf2_ID = temp_ID_2;

				control = true;
			}
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			/////////////////////////////////////////////////////////////FlexibleSECylinder_1_FlexibleTriangularSurface_2///////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			if (typeid(*db.surfaces[temp_ID_1 - 1]) == typeid(FlexibleSECylinder_1) && typeid(*db.surfaces[temp_ID_2 - 1]) == typeid(FlexibleTriangularSurface_2))
			{
				surf_pair[i][j] = new FlexibleSECylinder_1_FlexibleTriangularSurface_2();
				surf_pair[i][j]->inverted = false;
				surf_pair[i][j]->surf1_ID = temp_ID_1;
				surf_pair[i][j]->surf2_ID = temp_ID_2;

				control = true;
			}
			/*if (typeid(*db.surfaces[temp_ID_1 - 1]) == typeid(RigidArcRevolution_1) && typeid(*db.surfaces[temp_ID_2 - 1]) == typeid(FlexibleArcExtrusion_1))
			{
				surf_pair[i][j] = new FlexibleArcExtrusion_1_RigidArcRevolution_1();
				surf_pair[i][j]->inverted = true;
				surf_pair[i][j]->surf1_ID = temp_ID_2;
				surf_pair[i][j]->surf2_ID = temp_ID_1;

				control = true;
			}*/

			if (control == false)
			{
				db.myprintf("Warning! No Surface Pair exists for surfaces %d and %d in contact SSSS %d. This pair is ignored during analysis.", temp_ID_1, temp_ID_2, number);
			}
			//Acrescentar aqui outros possíveis pares de superfícies
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			surf_pair[i][j]->write_report = write_report;
			surf_pair[i][j]->write_report_diverged = write_report_diverged;
			surf_pair[i][j]->specialLCP = specialLCP;
			surf_pair[i][j]->n_pointwise = number_pointwise_interactions;
			surf_pair[i][j]->PreCalc();
		}
	}
}