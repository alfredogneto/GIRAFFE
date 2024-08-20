#include "Spline.h"
#include "Database.h"
//Variáveis globais
extern
Database db;

Spline::Spline()
{
	number = 0;
	radius = 0;
	nodeset = 0;
	size_nodeset = 0;
	size_sp_nodes = 0;
	size_sp_elements = 0;
	sp_element = NULL;
}

Spline::~Spline()
{
	if (alloc == true) 
	{
	for (int i = 0; i < 3; i++)
	{
		delete x_sp_Ai[i];
		delete x_sp_d[i];
		delete x_sp_dd[i];
		delete x_sp_tangent[i];
		delete x_sp_normal[i];
		delete x_Ai[i];
	}

	delete[] x_sp_Ai;
	delete[] x_sp_d;
	delete[] x_sp_dd;
	delete[] x_sp_tangent;
	delete[] x_sp_normal;
	delete[] x_Ai;

	delete[] knot;
	delete[] sp0;
	delete[] sp1;
	delete[] sp2;
	delete[] sp2_d;
	delete[] sp1_dd;
	delete[] sp2_dd;
	/*if (sp2 != NULL) {
		for (int i = 0; i < size_sp_nodes; i++)
			delete sp2[i];
		delete[] sp2;
	}*/
	if (sp_element != NULL)
	{
		for (int i = 0; i < size_sp_elements; i++)
			delete sp_element[i];
		delete[] sp_element;
	}
	delete[] sp_elements_list;
}
}

bool Spline::Read(FILE *f)
{
	char s[1000];

	fscanf(f, "%s", s);
	if (!strcmp(s, "Spline"))
	{
		fscanf(f, "%s", s);
		number = atoi(s);
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "R"))
	{
		fscanf(f, "%s", s);
		radius = atof(s);
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "NodeSet"))
	{
		fscanf(f, "%s", s);
		nodeset = atoi(s);
	}
	else
		return false;
	return true;
}

void Spline::Write(FILE *f)
{
	fprintf(f, "Spline\t%d\tR\t%.6e\tNodeSet\t%d\n",
		number,
		radius,
		nodeset
	);
	//Escrever o knot vector
	//for (int i = 0; i < (6 + db.node_sets[nodeset - 1]->n_nodes - 3); i++)
	//{
	//	fprintf(f, "%.6e\n", knot[i]);
	//}
}

void Spline::PreCalc()
{
	//Número de nós do nodeset
	size_nodeset = db.node_sets[nodeset - 1]->n_nodes;
	//Lista de nós do nodeset
	nodeset_list = db.node_sets[nodeset - 1]->node_list;
	//Número de elementos de spline (segmentos/trechos)
	size_sp_elements = db.node_sets[nodeset - 1]->n_nodes - 2;

	//Vetor de elementos de spline
	sp_element = new SplineElement*[size_sp_elements];
	//ID de cada elemento de spline
	sp_elements_list = new int[size_sp_elements];
	for (int i = 0; i < size_sp_elements; i++) {
		sp_element[i] = new SplineElement();
		sp_elements_list[i] = i;
	}

	//Criando o knot vector entre 0 e 1 de acordo com o número de pontos do nodeset
	knot = new double[size_nodeset + 3];
	knot[0] = 0;
	knot[1] = 0;
	knot[2] = 0;
	for (int i = 0; i < (size_nodeset - 3); i++)
	{
		knot[3 + i] = (i + 1) / (1.0*size_nodeset - 2);
	}
	knot[size_nodeset] = 1;
	knot[size_nodeset + 1] = 1;
	knot[size_nodeset + 2] = 1;

	//Definindo as variáveis de cada elemento da spline
	for (int i = 0; i < size_sp_elements; i++) {

		*sp_element[i]->radius = radius;

		sp_element[i]->knot_element[0] = knot[i];
		sp_element[i]->knot_element[1] = knot[i + 1];
		sp_element[i]->knot_element[2] = knot[i + 2];
		sp_element[i]->knot_element[3] = knot[i + 3];
		sp_element[i]->knot_element[4] = knot[i + 4];
		sp_element[i]->knot_element[5] = knot[i + 5];

		sp_element[i]->nodes[0] = nodeset_list[i];
		sp_element[i]->nodes[1] = nodeset_list[i + 1];
		sp_element[i]->nodes[2] = nodeset_list[i + 2];

		(*sp_element[i]->x_Ai)(0, 0) = db.nodes[nodeset_list[i] - 1]->copy_coordinates[0];
		(*sp_element[i]->x_Ai)(1, 0) = db.nodes[nodeset_list[i] - 1]->copy_coordinates[1];
		(*sp_element[i]->x_Ai)(2, 0) = db.nodes[nodeset_list[i] - 1]->copy_coordinates[2];

		(*sp_element[i]->x_Bi)(0, 0) = db.nodes[nodeset_list[i + 1] - 1]->copy_coordinates[0];
		(*sp_element[i]->x_Bi)(1, 0) = db.nodes[nodeset_list[i + 1] - 1]->copy_coordinates[1];
		(*sp_element[i]->x_Bi)(2, 0) = db.nodes[nodeset_list[i + 1] - 1]->copy_coordinates[2];

		(*sp_element[i]->x_Ci)(0, 0) = db.nodes[nodeset_list[i + 2] - 1]->copy_coordinates[0];
		(*sp_element[i]->x_Ci)(1, 0) = db.nodes[nodeset_list[i + 2] - 1]->copy_coordinates[1];
		(*sp_element[i]->x_Ci)(2, 0) = db.nodes[nodeset_list[i + 2] - 1]->copy_coordinates[2];
	}

	//Quantidade total de pontos para plotagem da spline
	size_sp_nodes = 2 * size_nodeset - 3;

	x_Ai = new Matrix*[size_nodeset];
	//x_Ap = new Matrix*[size_nodeset];
	//u_A = new Matrix*[size_nodeset];

	for (int i = 0; i < size_nodeset; i++)
	{
		x_Ai[i] = new Matrix(3);
		//x_Ap[i] = new Matrix(3);
		//u_A[i] = new Matrix(3);
	}

	for (int i = 0; i < size_nodeset; i++)
	{
		(*x_Ai[i])(0, 0) = db.nodes[nodeset_list[i] - 1]->copy_coordinates[0];
		(*x_Ai[i])(1, 0) = db.nodes[nodeset_list[i] - 1]->copy_coordinates[1];
		(*x_Ai[i])(2, 0) = db.nodes[nodeset_list[i] - 1]->copy_coordinates[2];
	}

	x_sp_Ai = new Matrix*[size_sp_nodes];
	for (int i = 0; i < size_sp_nodes; i++)
		x_sp_Ai[i] = new Matrix(3);


	x_sp_d = new Matrix*[size_sp_elements * 3];
	x_sp_dd = new Matrix*[size_sp_elements * 3];
	for (int i = 0; i < size_sp_elements * 3; i++) {
		x_sp_d[i] = new Matrix(3);
		x_sp_dd[i] = new Matrix(3);
	}

	x_sp_tangent = new Matrix*[size_sp_elements * 3];
	x_sp_normal = new Matrix*[size_sp_elements * 3];
	for (int i = 0; i < size_sp_elements * 3; i++) {
		x_sp_tangent[i] = new Matrix(3);
		x_sp_normal[i] = new Matrix(3);
	}

	//Calculando os pesos para todos os trechos da spline
	sp0 = new double[size_nodeset + 2];

	sp1 = new double[size_nodeset + 1];
	sp1_dd = new double[size_nodeset + 1];	//Segunda derivada

	sp2 = new double*[size_sp_nodes];
	sp2_d = new double*[size_sp_elements * 3];		//Primeira derivada
	sp2_dd = new double*[size_sp_elements * 3];		//Segunda derivada

	for (int i = 0; i < size_sp_nodes; i++)
		sp2[i] = new double[size_nodeset];

	for (int i = 0; i < size_sp_elements * 3; i++)
	{
		sp2_d[i] = new double[size_nodeset];
		sp2_dd[i] = new double[size_nodeset];
	}

	//Calculo dos coeficientes (pesos) para os pontos da spline
	for (int j = 0; j < size_sp_nodes; j++)
	{
		double zeta = j / (1.0*size_sp_nodes - 1);
		double tol = 0.0000001;

		if (j == 1)
			zeta = zeta - tol;

		if (zeta != 1) {
			for (int i = 0; i < (size_nodeset + 2); i++) {
				if (knot[i] <= zeta && zeta < knot[i + 1])
					sp0[i] = 1;
				else
					sp0[i] = 0;
			}
		}

		double a, b;

		int p = 1;
		for (int i = 0; i < size_nodeset + 1; i++) {
			if (knot[i + p] == knot[i]) {
				a = 0;
			}
			else {
				a = sp0[i] * (zeta - knot[i]) / (knot[i + p] - knot[i]);
			}
			if (knot[i + p + 1] == knot[i + 1]) {
				b = 0;
			}
			else {
				b = sp0[i + 1] * (knot[i + p + 1] - zeta) / (knot[i + p + 1] - knot[i + 1]);
			}
			sp1[i] = a + b;
		}

		p = 2;
		for (int i = 0; i < size_nodeset + 1; i++) {
			if (knot[i + p] == knot[i]) {
				a = 0;
			}
			else
			{
				a = sp1[i] * (zeta - knot[i]) / (knot[i + p] - knot[i]);
			}

			if (knot[i + p + 1] == knot[i + 1])
			{
				b = 0;
			}
			else
			{
				b = sp1[i + 1] * (knot[i + p + 1] - zeta) / (knot[i + p + 1] - knot[i + 1]);
			}
			sp2[j][i] = a + b;
		}
	}

	//Calculo dos coeficientes (pesos) para os pontos da derivada e da segunda derivada da spline
	for (int j = 0; j < size_sp_elements * 3; j++)
	{
		double zeta_step = 1 / (1.0*size_sp_nodes - 1);
		int zeta_aux = j / 3;
		double zeta;
		double tol = 0.0000001;

		zeta = zeta_step * (j - zeta_aux);

		if ((j + 3) % 3 == 0)
			zeta = zeta + tol;

		if ((j + 1) % 3 == 0)
			zeta = zeta - tol;

		if (zeta != 1) {
			for (int i = 0; i < (size_nodeset + 2); i++) {
				if (knot[i] <= zeta && zeta < knot[i + 1])
					sp0[i] = 1;
				else
					sp0[i] = 0;
			}
		}

		double a, b;
		double a_d, b_d; //primeira derivada
		double a_dd, b_dd; //segunda derivada

		int p = 1;
		for (int i = 0; i < size_nodeset + 1; i++) {
			if (knot[i + p] == knot[i]) {
				a = 0;
				a_dd = 0;
			}
			else {
				a = sp0[i] * (zeta - knot[i]) / (knot[i + p] - knot[i]);
				a_dd = sp0[i] * (p) / (knot[i + p] - knot[i]);
			}
			if (knot[i + p + 1] == knot[i + 1]) {
				b = 0;
				b_dd = 0;
			}
			else {
				b = sp0[i + 1] * (knot[i + p + 1] - zeta) / (knot[i + p + 1] - knot[i + 1]);
				b_dd = sp0[i + 1] * (p) / (knot[i + p + 1] - knot[i + 1]);
			}
			sp1[i] = a + b;
			sp1_dd[i] = a_dd - b_dd;
		}

		p = 2;
		for (int i = 0; i < size_nodeset + 1; i++) {
			if (knot[i + p] == knot[i]) {
				a_d = 0;
				a_dd = 0;
			}
			else
			{
				a_d = sp1[i] * (p) / (knot[i + p] - knot[i]);
				a_dd = sp1_dd[i] * (p) / (knot[i + p] - knot[i]);
			}

			if (knot[i + p + 1] == knot[i + 1])
			{
				b_d = 0;
				b_dd = 0;
			}
			else
			{
				b_d = sp1[i + 1] * (p) / (knot[i + p + 1] - knot[i + 1]);
				b_dd = sp1_dd[i + 1] * (p) / (knot[i + p + 1] - knot[i + 1]);
			}
			sp2_d[j][i] = a_d - b_d;
			sp2_dd[j][i] = a_dd - b_dd;
		}
	}
	////Limpando o vetor de Splines
	//(*x_sp_Ai[j])(0, 0) = 0;
	//(*x_sp_Ai[j])(1, 0) = 0;
	//(*x_sp_Ai[j])(2, 0) = 0;

	//for (int k = 0; k < (size_nodeset); k++)
	//{
	//	(*x_sp_Ai[j])(0, 0) += n2[k] * (*x_Ai[k])(0, 0);
	//	(*x_sp_Ai[j])(1, 0) += n2[k] * (*x_Ai[k])(1, 0);
	//	(*x_sp_Ai[j])(2, 0) += n2[k] * (*x_Ai[k])(2, 0);
	//}

	alloc = true;
}

void Spline::CalculateSpline() {

	for (int j = 0; j < size_sp_nodes; j++)
	{
		//Limpando o vetor de Splines
		(*x_sp_Ai[j])(0, 0) = 0;
		(*x_sp_Ai[j])(1, 0) = 0;
		(*x_sp_Ai[j])(2, 0) = 0;

		for (int k = 0; k < (size_nodeset); k++)
		{
			(*x_sp_Ai[j])(0, 0) += sp2[j][k] * (*x_Ai[k])(0, 0);
			(*x_sp_Ai[j])(1, 0) += sp2[j][k] * (*x_Ai[k])(1, 0);
			(*x_sp_Ai[j])(2, 0) += sp2[j][k] * (*x_Ai[k])(2, 0);
		}
	}
}

void Spline::CalculateSplineTangentNormal()
{
	//Preenchendo os pontos com as derivadas
	for (int j = 0; j < size_sp_elements * 3; j++)
	{
		//Limpando a primeira derivada
		(*x_sp_d[j])(0, 0) = 0;
		(*x_sp_d[j])(1, 0) = 0;
		(*x_sp_d[j])(2, 0) = 0;

		//Limpando a segunda derivada
		(*x_sp_dd[j])(0, 0) = 0;
		(*x_sp_dd[j])(1, 0) = 0;
		(*x_sp_dd[j])(2, 0) = 0;

		for (int k = 0; k < size_nodeset; k++)
		{
			//Primeira derivada
			(*x_sp_d[j])(0, 0) += sp2_d[j][k] * (*x_Ai[k])(0, 0);
			(*x_sp_d[j])(1, 0) += sp2_d[j][k] * (*x_Ai[k])(1, 0);
			(*x_sp_d[j])(2, 0) += sp2_d[j][k] * (*x_Ai[k])(2, 0);

			//Segunda derivada
			(*x_sp_dd[j])(0, 0) += sp2_dd[j][k] * (*x_Ai[k])(0, 0);
			(*x_sp_dd[j])(1, 0) += sp2_dd[j][k] * (*x_Ai[k])(1, 0);
			(*x_sp_dd[j])(2, 0) += sp2_dd[j][k] * (*x_Ai[k])(2, 0);
		}
	}
	for (int j = 0; j < size_sp_elements; j++)
	{
		//Testando colinearidade entre pontos
		double tol = 0.00000001;
		double collinear_test = norm(cross((*x_Ai[j + 1] - *x_Ai[j]), (*x_Ai[j + 2] - *x_Ai[j])));
		if (collinear_test <tol && collinear_test >-tol) {
			//Calculando tangente e normal para pontos colineares
			//Tangente
			*x_sp_tangent[j * 3] = (*x_Ai[j + 1] - *x_Ai[j]) * (1 / norm((*x_Ai[j + 1] - *x_Ai[j])));
			*x_sp_tangent[j * 3 + 1] = *x_sp_tangent[j * 3];
			*x_sp_tangent[j * 3 + 2] = *x_sp_tangent[j * 3];

			//Normal
			Matrix unit(3); //Ponto arbitrário
			unit(0, 0) = 1;
			unit(1, 0) = 1;
			unit(2, 0) = 1;
			if (norm(cross(*x_sp_tangent[j * 3], *x_sp_tangent[j * 3] + unit)) == 0) //Testando alinhamento dos vetores
				unit(1, 0) = 1 + tol;
			*x_sp_normal[j * 3] = cross(*x_sp_tangent[j * 3], *x_sp_tangent[j * 3] + unit);
			*x_sp_normal[j * 3] = *x_sp_normal[j * 3] * (1 / norm(*x_sp_normal[j * 3]));
			*x_sp_normal[j * 3 + 1] = *x_sp_normal[j * 3];
			*x_sp_normal[j * 3 + 2] = *x_sp_normal[j * 3];
		}
		else
		{
			//Calculando tangente e normal pelo Triedro de Frenet
			//Tangente
			*x_sp_tangent[j * 3] = *x_sp_d[j * 3] * (1 / norm(*x_sp_d[j * 3]));
			*x_sp_tangent[j * 3 + 1] = *x_sp_d[j * 3 + 1] * (1 / norm(*x_sp_d[j * 3 + 1]));
			*x_sp_tangent[j * 3 + 2] = *x_sp_d[j * 3 + 2] * (1 / norm(*x_sp_d[j * 3 + 2]));

			////Normal
			*x_sp_normal[j * 3] = (cross(*x_sp_d[j * 3], cross(*x_sp_dd[j * 3], *x_sp_d[j * 3])))*(1 / (norm(*x_sp_d[j * 3])*norm(cross(*x_sp_dd[j * 3], *x_sp_d[j * 3]))));
			*x_sp_normal[j * 3 + 1] = (cross(*x_sp_d[j * 3 + 1], cross(*x_sp_dd[j * 3 + 1], *x_sp_d[j * 3 + 1])))*(1 / (norm(*x_sp_d[j * 3 + 1])*norm(cross(*x_sp_dd[j * 3 + 1], *x_sp_d[j * 3 + 1]))));
			*x_sp_normal[j * 3 + 2] = (cross(*x_sp_d[j * 3 + 2], cross(*x_sp_dd[j * 3 + 2], *x_sp_d[j * 3 + 2])))*(1 / (norm(*x_sp_d[j * 3 + 2])*norm(cross(*x_sp_dd[j * 3 + 2], *x_sp_d[j * 3 + 2]))));
		}
	}
}

//void Spline::FillNodes()
//{
//	//for (int i = 0; i < size_nodeset; i++)
//	//{
//	//	(*u_A[i])(0, 0) = db.nodes[nodeset_list[i] - 1]->displacements[0];
//	//	(*u_A[i])(1, 0) = db.nodes[nodeset_list[i] - 1]->displacements[1];
//	//	(*u_A[i])(2, 0) = db.nodes[nodeset_list[i] - 1]->displacements[2];
//	//}
//	//for (int i = 0; i < size_nodeset; i++)
//	//{
//	//	*x_Ap[i] = *x_Ai[i] + *u_A[i];
//	//}
//}

void Spline::SaveConfiguration()
{
	for (int i = 0; i < size_nodeset; i++)
	{
		(*x_Ai[i])(0, 0) = db.nodes[nodeset_list[i] - 1]->copy_coordinates[0];
		(*x_Ai[i])(1, 0) = db.nodes[nodeset_list[i] - 1]->copy_coordinates[1];
		(*x_Ai[i])(2, 0) = db.nodes[nodeset_list[i] - 1]->copy_coordinates[2];

	}
	//for (int i = 0; i < size_sp_elements; i++) {

	//	sp_element[i]->SaveConfiguration();
	//}
}

void Spline::WriteVTK_XML_SplineMesh(FILE *f) {
	if (db.post_files->WriteSpline_flag == true)
	{
		CalculateSpline();
		//vetores para escrita no formato binário - usando a função 'enconde'
		std::vector<float> float_vector;
		std::vector<int> int_vector;
		Matrix* vec_P;
		vec_P = new Matrix(3);

		//Número de pontos a serem gerados
		int n_points = size_sp_nodes;
		//Número de células a serem geradas
		int n_cells = (size_sp_nodes - 1) / 2;
		//Opens Piece
		fprintf(f, "\t\t<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", n_points, n_cells);
		//Opens Points
		fprintf(f, "\t\t\t<Points>\n");
		//Opens DataArray
		fprintf(f, "\t\t\t\t<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"binary\">\n");
		float_vector.clear();
		//Preenchendo as coordenadas dos pontos
		for (int i = 0; i < size_sp_nodes; i++)//percorrendo os nós
		{
			//Posição de cada ponto 
			(*vec_P)(0, 0) = (*x_sp_Ai[i])(0, 0);
			(*vec_P)(1, 0) = (*x_sp_Ai[i])(1, 0);
			(*vec_P)(2, 0) = (*x_sp_Ai[i])(2, 0);
			//(*vec_P).fprint("vec.txt");
			float_vector.push_back((float)((*vec_P)(0, 0)));
			float_vector.push_back((float)((*vec_P)(1, 0)));
			float_vector.push_back((float)((*vec_P)(2, 0)));
		}
		fprintf(f, encodeData<float>(float_vector).c_str());
		fprintf(f, "\n");
		//Closes DataArray
		fprintf(f, "\t\t\t\t</DataArray>\n");
		//Closes Points
		fprintf(f, "\t\t\t</Points>\n");

		int nodes[3];
		//Opens Cells
		fprintf(f, "\t\t\t<Cells>\n");
		//Opens DataArray
		fprintf(f, "\t\t\t\t<DataArray type=\"Int32\" Name=\"connectivity\" format=\"binary\">\n");
		int_vector.clear();
		for (int cell = 0; cell < n_cells; cell++)
		{
			nodes[0] = cell * 2;
			nodes[1] = cell * 2 + 2;
			nodes[2] = cell * 2 + 1;

			int_vector.push_back(nodes[0]);
			int_vector.push_back(nodes[1]);
			int_vector.push_back(nodes[2]);
		}

		fprintf(f, encodeData(int_vector).c_str());
		fprintf(f, "\n");
		//Closes DataArray
		fprintf(f, "\t\t\t\t</DataArray>\n");
		//Opens DataArray
		int_vector.clear();
		fprintf(f, "\t\t\t\t<DataArray type=\"Int32\" Name=\"types\" format=\"binary\">\n");
		for (int cell = 0; cell < n_cells; cell++)
			int_vector.push_back(21);
		fprintf(f, encodeData(int_vector).c_str());
		fprintf(f, "\n");
		//Closes DataArray
		fprintf(f, "\t\t\t\t</DataArray>\n");
		//Opens DataArray
		fprintf(f, "\t\t\t\t<DataArray type=\"Int32\" Name=\"offsets\" format=\"binary\">\n");
		int cur_off = 0;
		int_vector.clear();
		for (int cell = 0; cell < n_cells; cell++)
		{
			cur_off += 3;
			int_vector.push_back(cur_off);
		}
		fprintf(f, encodeData(int_vector).c_str());
		fprintf(f, "\n");
		//Closes DataArray
		fprintf(f, "\t\t\t\t</DataArray>\n");
		//Closes Cells
		fprintf(f, "\t\t\t</Cells>\n");
		//Closes Piece
		fprintf(f, "\t\t</Piece>\n");
	}
}

void Spline::WriteVTK_XML_SplineRender(FILE * f)
{
	if (db.post_files->WriteRenderSpline_flag == true && db.splines_exist == true)
	{
		CalculateSpline();
		CalculateSplineTangentNormal();
		int n_circ = 16;
		double theta = 0;

		//vetores para escrita no formato binário - usando a função 'enconde'
		std::vector<float> float_vector;
		std::vector<int> int_vector;

		Matrix vec_x(3);
		Matrix vec_P(3);
		Matrix vec_P_radius(3);
		//Número de pontos a serem gerados
		int n_points = size_sp_elements * 40;
		//Número de células a serem geradas
		int n_cells = size_sp_elements * 8;
		//Opens Piece
		fprintf(f, "\t\t<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", n_points, n_cells);
		//Opens Points
		fprintf(f, "\t\t\t<Points>\n");
		//Opens DataArray
		fprintf(f, "\t\t\t\t<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"binary\">\n");
		float_vector.clear();
		//Preenchendo as coordenadas dos pontos
		for (int i = 0; i < size_sp_elements * 3; i++)//percorrendo os cada ponto da spline
		{
			theta = 0;
			vec_P_radius = *x_sp_normal[i] * radius;
			int i_aux = i / 3;
			vec_x = *x_sp_Ai[i - i_aux];
			if ((i + 2) % 3 == 0) {
				//Pontos intermediários do trecho de spline
				for (int point = 0; point < n_circ; point = point + 2)//Percorre os nós que descrevem o perímetro da ST
				{
					//Atualização do angulo theta
					theta = 2 * (2 * PI) / n_circ;
					//Rotação da normal a spline pela fórmula de rodrigues em torno da tangente
					vec_P_radius = vec_P_radius * cos(theta) + cross(*x_sp_tangent[i], vec_P_radius)*sin(theta) + *x_sp_tangent[i] * dot(*x_sp_tangent[i], vec_P_radius)*(1 - cos(theta));

					//Posição de cada ponto P da spline junto ao ponto da superfície
					vec_P(0, 0) = vec_x(0, 0) + vec_P_radius(0, 0);
					vec_P(1, 0) = vec_x(1, 0) + vec_P_radius(1, 0);
					vec_P(2, 0) = vec_x(2, 0) + vec_P_radius(2, 0);

					float_vector.push_back((float)(vec_P(0, 0)));
					float_vector.push_back((float)(vec_P(1, 0)));
					float_vector.push_back((float)(vec_P(2, 0)));
				}
			}
			else
			{
				//Ponto inicial ou final do trecho de spline
				for (int point = 0; point < n_circ; point++)//Percorre os nós que descrevem o perímetro da ST
				{
					//Atualização do angulo theta
					theta = (2 * PI) / n_circ;
					//Rotação da normal a spline pela fórmula de rodrigues em torno da tangente
					vec_P_radius = vec_P_radius * cos(theta) + cross(*x_sp_tangent[i], vec_P_radius)*sin(theta) + *x_sp_tangent[i] * dot(*x_sp_tangent[i], vec_P_radius)*(1 - cos(theta));

					//Posição de cada ponto P da spline junto ao ponto da superfície
					vec_P(0, 0) = vec_x(0, 0) + vec_P_radius(0, 0);
					vec_P(1, 0) = vec_x(1, 0) + vec_P_radius(1, 0);
					vec_P(2, 0) = vec_x(2, 0) + vec_P_radius(2, 0);

					float_vector.push_back((float)(vec_P(0, 0)));
					float_vector.push_back((float)(vec_P(1, 0)));
					float_vector.push_back((float)(vec_P(2, 0)));
				}
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
		for (int cell = 0; cell < size_sp_elements; cell++)
		{
			for (int cell_aux = 0; cell_aux < 8; cell_aux++) {
				if (cell_aux == 7) {
					nodes[0] = cell * 40 + cell_aux * 2;
					nodes[4] = cell * 40 + cell_aux * 2 + 1;
					nodes[1] = cell * 40 + 0;
					nodes[7] = cell * 40 + 16 + cell_aux;
					nodes[5] = cell * 40 + 16;
					nodes[3] = cell * 40 + 24 + cell_aux * 2;
					nodes[6] = cell * 40 + 24 + cell_aux * 2 + 1;
					nodes[2] = cell * 40 + 24;
				}
				else {
					nodes[0] = cell * 40 + cell_aux * 2;
					nodes[4] = cell * 40 + cell_aux * 2 + 1;
					nodes[1] = cell * 40 + cell_aux * 2 + 2;
					nodes[7] = cell * 40 + 16 + cell_aux;
					nodes[5] = cell * 40 + 16 + cell_aux + 1;
					nodes[3] = cell * 40 + 24 + cell_aux * 2;
					nodes[6] = cell * 40 + 24 + cell_aux * 2 + 1;
					nodes[2] = cell * 40 + 24 + cell_aux * 2 + 2;
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
		}

		fprintf(f, encodeData(int_vector).c_str());
		fprintf(f, "\n");
		//Closes DataArray
		fprintf(f, "\t\t\t\t</DataArray>\n");
		//Opens DataArray
		int_vector.clear();
		fprintf(f, "\t\t\t\t<DataArray type=\"Int32\" Name=\"types\" format=\"binary\">\n");
		for (int cell = 0; cell < n_cells; cell++)
			int_vector.push_back(23);
		fprintf(f, encodeData(int_vector).c_str());
		fprintf(f, "\n");
		//Closes DataArray
		fprintf(f, "\t\t\t\t</DataArray>\n");
		//Opens DataArray
		fprintf(f, "\t\t\t\t<DataArray type=\"Int32\" Name=\"offsets\" format=\"binary\">\n");
		int cur_off = 0;
		int_vector.clear();
		for (int cell = 0; cell < n_cells; cell++)
		{
			cur_off += 8;
			int_vector.push_back(cur_off);
		}
		fprintf(f, encodeData(int_vector).c_str());
		fprintf(f, "\n");
		//Closes DataArray
		fprintf(f, "\t\t\t\t</DataArray>\n");
		//Closes Cells
		fprintf(f, "\t\t\t</Cells>\n");
		//Closes Piece
		fprintf(f, "\t\t</Piece>\n");
		//for (int i = 0; i < n_circ; i++)
		//	delete[] points[i];
		//delete[]points;
	}
}

bool Spline::Check()
{
	//if ((db.node_sets[nodeset - 1]->n_nodes - 1) % 2 != 0 || (db.node_sets[nodeset - 1]->n_nodes - 1) < 3 || nodeset > db.number_node_sets)
	//{
	//	db.myprintf("Check number of nodes in nodeset %d.\n", this->nodeset);
	//	return false;
	//}
	return true;
}

