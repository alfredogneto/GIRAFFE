#include "SolidSection.h"

#include"Database.h"

//Variáveis globais
extern
Database db;


SolidSection::SolidSection()
{
	n_points = 0;
	points = NULL;
	number = 0;				//ID
	axis_position[0] = 0.0;	//Ref position
	axis_position[1] = 0.0;	//Ref position
	sprintf(section_type, "SolidSection");
}
//Construtor de cópia
SolidSection::SolidSection(SolidSection &copied)
{
	n_points = 0;
	points = NULL;
	number = 0;				//ID
	axis_position[0] = 0.0;	//Ref position
	axis_position[1] = 0.0;	//Ref position
	sprintf(section_type, "SolidSection");

	Alloc(copied.n_points);
	//Copia os valores dos pontos
	for (long i = 0; i < copied.n_points; i++)
	{
		points[i][0] = copied.points[i][0];
		points[i][1] = copied.points[i][1];
	}
	strcpy(section_type, copied.section_type);
	axis_position[0] = copied.axis_position[0];
	axis_position[1] = copied.axis_position[1];
}

void SolidSection::Alloc(int e_points)
{
	//Desaloca matrizes
	if (points != NULL)
	{
		for (int i = 0; i < n_points; i++)
			delete points[i];
		delete[]points;
	}
	n_points = e_points;
	//Aloca matrizes
	points = new double*[n_points];
	for (int i = 0; i < n_points; i++)
		points[i] = new double[2];
}

SolidSection::~SolidSection()
{
	//Desaloca matrizes
	if (points != NULL)
	{
		for (int i = 0; i < n_points; i++)
			delete points[i];
		delete[]points;
	}
}
bool SolidSection::Read(FILE *f)
{
	char s[1000];
	fscanf(f, "%s", s);
	number = atoi(s);
	
	fscanf(f, "%s", s);
	if (!strcmp(s, "AxisPosition"))
	{
		fscanf(f, "%s", s);
		axis_position[0] = atof(s);
		fscanf(f, "%s", s);
		axis_position[1] = atof(s);
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "NPoints"))
	{
		fscanf(f, "%s", s);
		n_points = atoi(s);
	}
	else
		return false;

	Alloc(n_points);
	//Leitura dos pontos
	double x, y;
	int index;
	for (int p = 0; p < n_points; p++)
	{
		fscanf(f, "%s", s);
		if (!strcmp(s, "Point"))
		{
			fscanf(f, "%s", s);
			index = atoi(s);
			fscanf(f, "%s", s);
			x = atof(s);
			fscanf(f, "%s", s);
			y = atof(s);
			points[p][0] = x;
			points[p][1] = y;
		}
		else
			return false;
	}
	return true;
}
void SolidSection::Write(FILE *f)
{
	fprintf(f, section_type);
	fprintf(f, "\t%d\t", number);
	fprintf(f, "\tAxisPosition\t%.6f\t%.6f\t", axis_position[0], axis_position[1]);
	fprintf(f, "NPoints\t%d\n", n_points);
	for (int p = 0; p < n_points; p++)
	{
		fprintf(f, "Point\t%d\t%.6f\t%.6f\n", p + 1, points[p][0], points[p][1]);
	}
}

void SolidSection::WriteVTK_XMLRender(FILE *f, Beam_1* elem)
{
	//vetores para escrita no formato binário - usando a função 'enconde'
	std::vector<float> float_vector;
	std::vector<int> int_vector;
	
	Matrix Q(3, 3);
	Matrix alpha_1(3);
	double alpha_escalar;
	Matrix A;
	Matrix vec_P(3);
	//Número de pontos a serem gerados
	int n_p = elem->n_nodes*n_points;
	//Número de células a serem geradas - laterais + duas "tampas"
	int n_cells = n_points*(elem->n_nodes - 1) + 2;
	//Opens Piece
	fprintf(f, "\t\t<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", n_p, n_cells);
	//Opens Points
	fprintf(f, "\t\t\t<Points>\n");
	//Opens DataArray
	fprintf(f, "\t\t\t\t<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"binary\">\n");
	float_vector.clear();
	//Preenchendo as coordenadas dos pontos
	for (int i = 0; i < elem->n_nodes; i++)//percorrendo os 3 nós do elemento
	{
		//Para cada nó do elemento, calcula o tensor rotação
		alpha_1(0, 0) = db.nodes[elem->nodes[i] - 1]->copy_coordinates[3];
		alpha_1(1, 0) = db.nodes[elem->nodes[i] - 1]->copy_coordinates[4];
		alpha_1(2, 0) = db.nodes[elem->nodes[i] - 1]->copy_coordinates[5];
		alpha_escalar = norm(alpha_1);
		A = skew(alpha_1);
		double g = 4.0 / (4.0 + alpha_escalar*alpha_escalar);
		Q = *elem->I3 + g*(A + 0.5*A*A);
		Q = Q*transp(*elem->transform3);//Matriz de transformação para trazer o vetor da ST do plano xy para o plano atual em que ela se encontra

		//Percorre os nós que descrevem o perímetro da ST
		for (int point = 0; point < n_points; point++)
		{
			//Posição de cada ponto P no plano xy (referência)
			vec_P(0, 0) = points[point][0] - db.sections[elem->section - 1]->sec_details->axis_position[0];
			vec_P(1, 0) = points[point][1] - db.sections[elem->section - 1]->sec_details->axis_position[1];
			vec_P(2, 0) = 0.0;
			vec_P = Q*vec_P;//Operando rotacionando para o sistema da barra
			for (int c = 0; c < 3; c++)
				vec_P(c, 0) += db.nodes[elem->nodes[i] - 1]->ref_coordinates[c] + db.post_files->mag_factor*(db.nodes[elem->nodes[i] - 1]->copy_coordinates[c] - db.nodes[elem->nodes[i] - 1]->ref_coordinates[c]);//Translação - soma a posição do centro da barra (do ponto em questão)
			float_vector.push_back((float)(vec_P(0, 0)));
			float_vector.push_back((float)(vec_P(1, 0)));
			float_vector.push_back((float)(vec_P(2, 0)));
		}
	}
	fprintf(f, encodeData<float>(float_vector).c_str());
	fprintf(f, "\n");
	//Closes DataArray
	fprintf(f, "\t\t\t\t</DataArray>\n");
	//Closes Points
	fprintf(f, "\t\t\t</Points>\n");

	int nodes[8];
	//Opens Cells
	fprintf(f, "\t\t\t<Cells>\n");
	//Opens DataArray
	fprintf(f, "\t\t\t\t<DataArray type=\"Int32\" Name=\"connectivity\" format=\"binary\">\n");
	int_vector.clear();
	for (int index = 0; index < n_points; index++)
	{
		//Células com os pontos do perímetro
		if (index != n_points-1)
		{
			nodes[0] = index;
			nodes[1] = index + n_points;
			nodes[2] = index + n_points + 1;
			nodes[3] = index + 1;

			nodes[4] = index + n_points;
			nodes[5] = index + 2*n_points;
			nodes[6] = index + 2*n_points + 1;
			nodes[7] = index + n_points + 1;
		}
		//última célula do perímetro (relacionada com a numeração da primeira - para fechar a figura)
		else
		{
			nodes[0] = index;
			nodes[1] = index + n_points;
			nodes[2] = n_points;
			nodes[3] = 0;

			nodes[4] = index + n_points;
			nodes[5] = index + 2 * n_points;
			nodes[6] = 2 * n_points;
			nodes[7] = n_points;
		}
		
		int_vector.push_back(nodes[0]);
		int_vector.push_back(nodes[1]);
		int_vector.push_back(nodes[2]);
		int_vector.push_back(nodes[3]);
		int_vector.push_back(nodes[4]);
		int_vector.push_back(nodes[5]);
		int_vector.push_back(nodes[6]);
		int_vector.push_back(nodes[7]);
	}
	//Células com os pontos das tampas
	for (int index = 0; index < n_points; index++)
		int_vector.push_back(index);
	int_vector.push_back(0);
	for (int index = 0; index < n_points; index++)
		int_vector.push_back(index + 2*n_points);
	int_vector.push_back(2 * n_points);

	fprintf(f, encodeData(int_vector).c_str());
	fprintf(f, "\n");
	//Closes DataArray
	fprintf(f, "\t\t\t\t</DataArray>\n");
	//Opens DataArray
	int_vector.clear();
	fprintf(f, "\t\t\t\t<DataArray type=\"Int32\" Name=\"types\" format=\"binary\">\n");
	for (int cell = 0; cell < n_cells-2; cell++)
	{
		int_vector.push_back(9);
		//fprintf(f, "\t\t\t\t\t%d\n", 13);
	}
	//tampas
	int_vector.push_back(7);
	int_vector.push_back(7);

	fprintf(f, encodeData(int_vector).c_str());
	fprintf(f, "\n");
	//Closes DataArray
	fprintf(f, "\t\t\t\t</DataArray>\n");
	//Opens DataArray
	fprintf(f, "\t\t\t\t<DataArray type=\"Int32\" Name=\"offsets\" format=\"binary\">\n");
	int cur_off = 0;
	int_vector.clear();
	for (int cell = 0; cell < n_cells-2; cell++)
	{
		cur_off += 4;
		int_vector.push_back(cur_off);
	}
	//tampas
	cur_off += n_points+1;
	int_vector.push_back(cur_off);
	cur_off += n_points+1;
	int_vector.push_back(cur_off);
	
	fprintf(f, encodeData(int_vector).c_str());
	fprintf(f, "\n");
	//Closes DataArray
	fprintf(f, "\t\t\t\t</DataArray>\n");
	//Closes Cells
	fprintf(f, "\t\t\t</Cells>\n");

	//Opens CellData
	fprintf(f, "\t\t\t<CellData FieldData=\"ElementData\">\n");
	float_vector.clear();
	//Opens DataArray
	fprintf(f, "\t\t\t\t<DataArray Name=\"ElementResults\" type=\"Float32\" NumberOfComponents=\"%d\" format=\"binary\">\n", db.n_element_results);
	for (int cell = 0; cell < n_cells; cell++)
	{
		//Imprime os resultados do elemento
		int res_element = 6;
		float_vector.push_back((float)((*elem->sigma_r[0])(0, 0) + (*elem->sigma_r[1])(0, 0)) / 2);
		float_vector.push_back((float)((*elem->sigma_r[0])(1, 0) + (*elem->sigma_r[1])(1, 0)) / 2);
		float_vector.push_back((float)((*elem->sigma_r[0])(2, 0) + (*elem->sigma_r[1])(2, 0)) / 2);
		float_vector.push_back((float)((*elem->sigma_r[0])(3, 0) + (*elem->sigma_r[1])(3, 0)) / 2);
		float_vector.push_back((float)((*elem->sigma_r[0])(4, 0) + (*elem->sigma_r[1])(4, 0)) / 2);
		float_vector.push_back((float)((*elem->sigma_r[0])(5, 0) + (*elem->sigma_r[1])(5, 0)) / 2);
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
		int_vector.push_back(1);				//Element ID - viga
		int_vector.push_back(elem->material);
		int_vector.push_back(elem->section);
		int_vector.push_back(elem->cs);
	}
	fprintf(f, encodeData(int_vector).c_str());
	fprintf(f, "\n");
	//Closes DataArray
	fprintf(f, "\t\t\t\t</DataArray>\n");
	//Closes CellData
	fprintf(f, "\t\t\t</CellData>\n");
	//Closes Piece
	fprintf(f, "\t\t</Piece>\n");
}