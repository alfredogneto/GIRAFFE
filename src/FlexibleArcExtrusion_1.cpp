#include "FlexibleArcExtrusion_1.h"

#include "Dynamic.h"
#include "ArcCirc.h"
#include "PostFiles.h"
#include "Node.h"
#include "CoordinateSystem.h"
#include "Encoding.h"

#include"Database.h"
//Variáveis globais
extern
Database db;

#define PI 3.1415926535897932384626433832795

FlexibleArcExtrusion_1::FlexibleArcExtrusion_1()
{
	nDOFs = 12;
	n_nodes = 2;
	number = 0;
	cs = 0;
	nodes = new int[n_nodes];
	VTK_nodes = new int[n_nodes];
	DOFs = new int *[n_nodes];
	for (int i = 0; i < n_nodes; i++)
		DOFs[i] = new int[db.number_GLs_node];
	VTK_nodes[0] = 0;
	VTK_nodes[1] = 1;
	//Rotina para ativar os GLS de cada nó
	for (int i = 0; i < n_nodes; i++)
	{
		for (int j = 0; j < db.number_GLs_node; j++)
		{
			DOFs[i][j] = 0;
		}
		//Translação
		DOFs[i][0] = 1;
		DOFs[i][1] = 1;
		DOFs[i][2] = 1;
		//Rotação
		DOFs[i][3] = 1;
		DOFs[i][4] = 1;
		DOFs[i][5] = 1;
	}
	GLs = new int*[nDOFs];
	for (int i = 0; i < nDOFs; i++)
		GLs[i] = NULL;


	x_AAi = new Matrix(3);
	x_BAi = new Matrix(3);
	Q_AAi = new Matrix(3, 3);
	Q_BAi = new Matrix(3, 3);
	Q_AAic = new Matrix(3, 3);
	Q_BAic = new Matrix(3, 3);
	Q0A = new Matrix(3, 3);
	d_A = new Matrix(12);
	dui_A = new Matrix(12);
	ddui_A = new Matrix(12);
	alpha_AA = new Matrix(3);
	alpha_BA = new Matrix(3);
	aQ_AAi = new double*[3];
	for (int i = 0; i < 3; i++)
		aQ_AAi[i] = new double[3];
	aQ_BAi = new double*[3];
	for (int i = 0; i < 3; i++)
		aQ_BAi[i] = new double[3];
	I3 = new Matrix(3, 3);
	(*I3)(0, 0) = 1.0;
	(*I3)(1, 1) = 1.0;
	(*I3)(2, 2) = 1.0;

	x_AAp = new Matrix(3);
	x_BAp = new Matrix(3);
	Q_AAp = new Matrix(3, 3);
	Q_BAp = new Matrix(3, 3);

	InitializeDegeneration();
}

void FlexibleArcExtrusion_1::SetMinMaxRange()
{
	u1_min = -1.0;
	u1_max = +1.0;
	u2_min = db.arcs[arc_ID - 1]->theta_i;
	u2_max = db.arcs[arc_ID - 1]->theta_f;
	
	u1_range = 2.0;

	if (u2_max >= u2_min)
		u2_range = abs(u2_max - u2_min);
	else
		u2_range = 2 * PI - abs(u2_max - u2_min);
}

FlexibleArcExtrusion_1::~FlexibleArcExtrusion_1()
{
	delete[] nodes;
	delete[] VTK_nodes;
	if (DOFs != NULL)
	{
		for (int i = 0; i < n_nodes; i++)
			delete[] DOFs[i];
		delete[] DOFs;
	}
	delete[]GLs;


	delete x_AAi;
	delete x_BAi;
	delete Q_AAi;
	delete Q_BAi;
	delete Q_AAic;
	delete Q_BAic;
	delete d_A;
	delete dui_A;
	delete ddui_A;
	delete alpha_AA;
	delete alpha_BA;

	for (int i = 0; i < 3; i++)
		delete[]aQ_AAi[i];
	delete[]aQ_AAi;
	for (int i = 0; i < 3; i++)
		delete[]aQ_BAi[i];
	delete[]aQ_BAi;
	delete I3;
	delete Q0A;

	delete x_AAp;
	delete x_BAp;
	delete Q_AAp;
	delete Q_BAp;
	FreeDegeneration();
}

//Obtem ponto da superficie
void FlexibleArcExtrusion_1::SurfacePoint(double& zeta, double& theta, Matrix& point)
{
	
	Matrix ar(3);
	ar(0, 0) = (*radius) *cos(theta) + (*c_point)(0, 0);
	ar(1, 0) = (*radius) *sin(theta) + (*c_point)(1, 0);
	ar(2, 0) = 0.0;
	//Gamma na posicao atual
	point = (1 - zeta) / 2 * (*x_AAp + (*Q_AAp) * ar) + (1 + zeta) / 2 * (*x_BAp + (*Q_BAp) * ar);
}

//Normal exterior à superfície na posição escolhida
void FlexibleArcExtrusion_1::NormalExt(double* zeta, double* theta, Matrix* n)
{
	double *d = d_A->getMatrix();		//ponteiro para o vetor d
	double* xAAi = x_AAi->getMatrix();
	double* xBAi = x_BAi->getMatrix();
	double** QAAi = aQ_AAi;
	double** QBAi = aQ_BAi;


	double *rad = radius;
	double *cpoint = c_point->getMatrix();

	double *next = n->getMatrix();		//ponteiro para o vetor np

	bool* normalint = &flag_normal_int;	//ponteiro para booleana flag_normal_int
	double v[1000];
	
	v[171] = 0.5e0*(1e0 + (*zeta));
	v[170] = 0.5e0 - 0.5e0*(*zeta);
	v[165] = (*rad)*sin((*theta));
	v[164] = (*rad)*cos((*theta));
	v[163] = Power(d[11], 2);
	v[162] = 0.5e0*d[11] * d[9];
	v[161] = 0.5e0*d[10];
	v[160] = Power(d[10], 2);
	v[168] = v[160] + v[163];
	v[159] = d[9] * v[161];
	v[158] = Power(d[9], 2);
	v[157] = Power(d[5], 2);
	v[156] = 0.5e0*d[3] * d[5];
	v[155] = 0.5e0*d[4];
	v[154] = Power(d[4], 2);
	v[166] = v[154] + v[157];
	v[153] = d[3] * v[155];
	v[152] = Power(d[3], 2);
	v[68] = d[5] * v[155];
	v[84] = d[11] * v[161];
	v[55] = 4e0 / (4e0 + v[152] + v[166]);
	v[167] = -0.5e0*v[55];
	v[58] = 1e0 + v[166] * v[167];
	v[59] = (-d[5] + v[153])*v[55];
	v[60] = (d[4] + v[156])*v[55];
	v[62] = (d[5] + v[153])*v[55];
	v[64] = 1e0 + (v[152] + v[157])*v[167];
	v[65] = v[55] * (-d[3] + v[68]);
	v[67] = (-d[4] + v[156])*v[55];
	v[69] = v[55] * (d[3] + v[68]);
	v[70] = 1e0 + (v[152] + v[154])*v[167];
	v[71] = 4e0 / (4e0 + v[158] + v[168]);
	v[169] = -0.5e0*v[71];
	v[74] = 1e0 + v[168] * v[169];
	v[75] = (-d[11] + v[159])*v[71];
	v[76] = (d[10] + v[162])*v[71];
	v[78] = (d[11] + v[159])*v[71];
	v[80] = 1e0 + (v[158] + v[163])*v[169];
	v[81] = v[71] * (-d[9] + v[84]);
	v[83] = (-d[10] + v[162])*v[71];
	v[85] = v[71] * (d[9] + v[84]);
	v[86] = 1e0 + (v[158] + v[160])*v[169];
	v[87] = QAAi[0][0] * v[58] + QAAi[1][0] * v[59] + QAAi[2][0] * v[60];
	v[88] = QAAi[0][1] * v[58] + QAAi[1][1] * v[59] + QAAi[2][1] * v[60];
	v[90] = QAAi[0][0] * v[62] + QAAi[1][0] * v[64] + QAAi[2][0] * v[65];
	v[91] = QAAi[0][1] * v[62] + QAAi[1][1] * v[64] + QAAi[2][1] * v[65];
	v[93] = QAAi[0][0] * v[67] + QAAi[1][0] * v[69] + QAAi[2][0] * v[70];
	v[94] = QAAi[0][1] * v[67] + QAAi[1][1] * v[69] + QAAi[2][1] * v[70];
	v[96] = QBAi[0][0] * v[74] + QBAi[1][0] * v[75] + QBAi[2][0] * v[76];
	v[97] = QBAi[0][1] * v[74] + QBAi[1][1] * v[75] + QBAi[2][1] * v[76];
	v[99] = QBAi[0][0] * v[78] + QBAi[1][0] * v[80] + QBAi[2][0] * v[81];
	v[100] = QBAi[0][1] * v[78] + QBAi[1][1] * v[80] + QBAi[2][1] * v[81];
	v[102] = QBAi[0][0] * v[83] + QBAi[1][0] * v[85] + QBAi[2][0] * v[86];
	v[103] = QBAi[0][1] * v[83] + QBAi[1][1] * v[85] + QBAi[2][1] * v[86];
	v[113] = cpoint[0] + v[164];
	v[114] = cpoint[1] + v[165];
	v[120] = v[170] * (-(v[165] * v[87]) + v[164] * v[88]) + v[171] * (-(v[165] * v[96]) + v[164] * v[97]);
	v[121] = v[170] * (-(v[165] * v[90]) + v[164] * v[91]) + v[171] * (v[100] * v[164] - v[165] * v[99]);
	v[122] = (v[103] * v[164] - v[102] * v[165])*v[171] + v[170] * (-(v[165] * v[93]) + v[164] * v[94]);
	v[126] = 1e0 / sqrt(Power(v[120], 2) + Power(v[121], 2) + Power(v[122], 2));
	v[125] = v[120] * v[126];
	v[127] = v[121] * v[126];
	v[128] = v[122] * v[126];
	v[131] = -0.5e0*(d[0] + v[113] * v[87] + v[114] * v[88] + xAAi[0]) + 0.5e0*(d[6] + v[113] * v[96] + v[114] * v[97]
		+ xBAi[0]);
	v[134] = -0.5e0*(d[1] + v[113] * v[90] + v[114] * v[91] + xAAi[1]) + 0.5e0*(d[7] + v[100] * v[114] + v[113] * v[99]
		+ xBAi[1]);
	v[137] = -0.5e0*(d[2] + v[113] * v[93] + v[114] * v[94] + xAAi[2]) + 0.5e0*(d[8] + v[102] * v[113] + v[103] * v[114]
		+ xBAi[2]);
	v[139] = 1e0 / sqrt(Power(v[131], 2) + Power(v[134], 2) + Power(v[137], 2));
	v[138] = v[131] * v[139];
	v[140] = v[134] * v[139];
	v[149] = -(v[127] * v[138]) + v[125] * v[140];
	v[141] = v[137] * v[139];
	v[147] = v[128] * v[138] - v[125] * v[141];
	if ((*normalint)) {
		v[143] = -1e0;
	}
	else {
		v[143] = 1e0;
	};
	v[144] = -(v[128] * v[140]) + v[127] * v[141];
	v[146] = 1e0 / sqrt(Power(v[144], 2) + Power(v[147], 2) + Power(v[149], 2));
	v[172] = v[143] * v[146];
	next[0] = v[144] * v[172];
	next[1] = v[147] * v[172];
	next[2] = v[149] * v[172];
}

void FlexibleArcExtrusion_1::WriteVTK_XMLRender(FILE *f)
{
	if (db.post_files->WriteFlexibleContactSurfaces_flag == true)
	{
		int n_circ = 24;

		//Número de pontos a serem gerados
		int n_points = n_circ + 1;
		double** points;
		points = new double*[n_points];
		for (int i = 0; i < n_points; i++)
			points[i] = new double[2];
		int index = 0;
		double theta = 0;

		//
		double theta_imod;
		double theta_fmod;
		//

		//Se os ângulos inicial e final esiverem do segundo para o terceiro quadrante
		if ((((db.arcs[arc_ID - 1]->theta_i) >= PI / 2 && (db.arcs[arc_ID - 1]->theta_i) <= PI)) && (((db.arcs[arc_ID - 1]->theta_f) < -PI / 2 && (db.arcs[arc_ID - 1]->theta_f) >= -PI)))
		{
			if ((db.arcs[arc_ID - 1]->theta_i) < 0)
				theta_imod = 2 * PI + (db.arcs[arc_ID - 1]->theta_i);
			else
				theta_imod = (db.arcs[arc_ID - 1]->theta_i);
			if ((db.arcs[arc_ID - 1]->theta_f) < 0)
				theta_fmod = 2 * PI + (db.arcs[arc_ID - 1]->theta_f);
			else
				theta_fmod = (db.arcs[arc_ID - 1]->theta_f);

			theta = theta_imod;

			for (int n = 0; n <= n_circ; n++)
			{
				points[index][0] = (db.arcs[arc_ID - 1]->radius)*cos(theta) + (db.arcs[arc_ID - 1]->c_point)(0, 0);
				points[index][1] = (db.arcs[arc_ID - 1]->radius)*sin(theta) + (db.arcs[arc_ID - 1]->c_point)(1, 0);
				index += 1;
				theta += (theta_fmod - theta_imod) / (n_circ - 1);
			}
		}
		else
		{
			for (int n = 0; n <= n_circ; n++)
			{
				theta = (db.arcs[arc_ID - 1]->theta_i) + (index * ((db.arcs[arc_ID - 1]->theta_f) - (db.arcs[arc_ID - 1]->theta_i))) / (n_circ - 1);
				points[index][0] = (db.arcs[arc_ID - 1]->radius)*cos(theta) + (db.arcs[arc_ID - 1]->c_point)(0, 0);
				points[index][1] = (db.arcs[arc_ID - 1]->radius)*sin(theta) + (db.arcs[arc_ID - 1]->c_point)(1, 0);
				index += 1;
			}
		}
		//vetores para escrita no formato binário - usando a função 'enconde'
		std::vector<float> float_vector;
		std::vector<int> int_vector;

		Matrix vec_P(3);
		//Número de pontos a serem gerados
		n_points = n_circ * 2;
		//Número de células a serem geradas
		int n_cells = n_circ;
		//Opens Piece
		fprintf(f, "\t\t<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", n_points, n_cells);
		//Opens Points
		fprintf(f, "\t\t\t<Points>\n");
		//Opens DataArray
		fprintf(f, "\t\t\t\t<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"binary\">\n");
		float_vector.clear();
		//Preenchendo as coordenadas dos pontos
		for (int i = 0; i < 2; i++)//percorrendo os 2 nós das extremidades
		{
			for (int point = 0; point < n_circ; point++)//Percorre os nós que descrevem o perímetro da ST
			{
				//Posição de cada ponto P no plano xy (referência)
				vec_P(0, 0) = points[point][0];
				vec_P(1, 0) = points[point][1];
				vec_P(2, 0) = 0.0;
				if (i == 0)
					vec_P = (*Q_AAi)*vec_P;//Operando rotacionando para o sistema da barra
				else
					vec_P = (*Q_BAi)*vec_P;//Operando rotacionando para o sistema da barra
				for (int c = 0; c < 3; c++)
					vec_P(c, 0) += db.nodes[nodes[i] - 1]->ref_coordinates[c] + db.post_files->mag_factor*(db.nodes[nodes[i] - 1]->copy_coordinates[c] - db.nodes[nodes[i] - 1]->ref_coordinates[c]);//Translação - soma a posição do centro da barra (do ponto em questão)
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

		int nodes[4];
		//Opens Cells
		fprintf(f, "\t\t\t<Cells>\n");
		//Opens DataArray
		fprintf(f, "\t\t\t\t<DataArray type=\"Int32\" Name=\"connectivity\" format=\"binary\">\n");
		int_vector.clear();
		for (int cell = 0; cell < n_cells; cell++)
		{
			if (cell != n_cells - 1)
			{
				nodes[0] = cell;
				nodes[1] = cell + n_cells;
				nodes[2] = cell + n_cells + 1;
				nodes[3] = cell + 1;
			}
			//else
			//{
			//	nodes[0] = cell;
			//	nodes[1] = 2 * n_cells - 1;
			//	nodes[2] = n_cells;
			//	nodes[3] = 0;
			//}
			int_vector.push_back(nodes[0]);
			int_vector.push_back(nodes[1]);
			int_vector.push_back(nodes[2]);
			int_vector.push_back(nodes[3]);
		}

		fprintf(f, encodeData(int_vector).c_str());
		fprintf(f, "\n");
		//Closes DataArray
		fprintf(f, "\t\t\t\t</DataArray>\n");
		//Opens DataArray
		int_vector.clear();
		fprintf(f, "\t\t\t\t<DataArray type=\"Int32\" Name=\"types\" format=\"binary\">\n");
		for (int cell = 0; cell < n_cells; cell++)
			int_vector.push_back(9);
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
			cur_off += 4;
			int_vector.push_back(cur_off);
		}
		fprintf(f, encodeData(int_vector).c_str());
		fprintf(f, "\n");
		//Closes DataArray
		fprintf(f, "\t\t\t\t</DataArray>\n");
		//Closes Cells
		fprintf(f, "\t\t\t</Cells>\n");
		
		
		//Opens CellData
		fprintf(f, "\t\t\t<CellData FieldData=\"SurfaceData\">\n");
		int_vector.clear();
		//Opens DataArray
		fprintf(f, "\t\t\t\t<DataArray Name=\"SurfaceProperties\" type=\"Int32\" NumberOfComponents=\"%d\" format=\"binary\">\n", 1);
		for (int cell = 0; cell < n_cells; cell++)
		{
			int_vector.push_back(number);		//Surface ID
		}
		fprintf(f, encodeData(int_vector).c_str());
		fprintf(f, "\n");
		//Closes DataArray
		fprintf(f, "\t\t\t\t</DataArray>\n");
		//Closes CellData
		fprintf(f, "\t\t\t</CellData>\n");


		//Closes Piece
		fprintf(f, "\t\t</Piece>\n");
		for (int i = 0; i < n_circ; i++)
			delete[] points[i];
		delete[]points;
	}
}

bool FlexibleArcExtrusion_1::Read(FILE *f)
{
	char s[1000];
	fscanf(f, "%s", s);
	number = atoi(s);

	fscanf(f, "%s", s);
	if (!strcmp(s, "Arc"))
	{
		fscanf(f, "%s", s);
		arc_ID = atoi(s);
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
	if (!strcmp(s, "Nodes"))
	{
		fscanf(f, "%s", s);
		nodes[0] = atoi(s);
		fscanf(f, "%s", s);
		nodes[1] = atoi(s);
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "Concave"))
		flag_normal_int = true;
	else
	{
		if (!strcmp(s, "Convex"))
			flag_normal_int = false;
		else
			return false;
	}

	ReadCommon(f);

	return true;
}

void FlexibleArcExtrusion_1::Write(FILE *f)
{
	char s[20];
	if (flag_normal_int == true)
		sprintf(s, "Concave");
	else
		sprintf(s, "Convex");
	fprintf(f, "FlexibleArcExtrusion_1\t%d\tArc\t%d\tCS\t%d\tNodes\t%d\t%d\t%s\n",
		number,
		arc_ID,
		cs,
		nodes[0],
		nodes[1],
		s);
}

//Checa inconsistências para evitar erros de execução
bool FlexibleArcExtrusion_1::Check()
{
	//Check nodes
	for (int i = 0; i < n_nodes; i++)
	{
		if (nodes[i] > db.number_nodes)
			return false;
	}
	//Check CS
	if (cs > db.number_CS)
		return false;
	//Check arc
	if (arc_ID > db.number_arcs)
		return false;
	return true;
}

//Realiza chute inicial para as variáveis zeta e theta
void FlexibleArcExtrusion_1::InitialGuess(Matrix* xS, double** convective, int n_solutions)
{
	//TODO
}

void FlexibleArcExtrusion_1::PreCalc()
{
	//Atribuindo valores do arco - ponteiros
	radius = &db.arcs[arc_ID - 1]->radius;
	c_point = &db.arcs[arc_ID - 1]->c_point;
	i_point = &db.arcs[arc_ID - 1]->i_point;
	f_point = &db.arcs[arc_ID - 1]->f_point;
	theta_i = &db.arcs[arc_ID - 1]->theta_i;
	theta_f = &db.arcs[arc_ID - 1]->theta_f;

	//Apontando para posição que indica valor dos GLs globais
	for (int i = 0; i < 6; i++)
	{
		GLs[i] = &db.nodes[nodes[0] - 1]->GLs[i];
		GLs[i + 6] = &db.nodes[nodes[1] - 1]->GLs[i];
	}
	
	//Transformação de coordenadas
	Matrix e1g(3);
	e1g(0, 0) = 1.0;
	Matrix e2g(3);
	e2g(1, 0) = 1.0;
	Matrix e3g(3);
	e3g(2, 0) = 1.0;
	Matrix e1l = *db.CS[cs - 1]->E1;
	Matrix e2l = *db.CS[cs - 1]->E2;
	Matrix e3l = *db.CS[cs - 1]->E3;

	//Salva a matriz de transformação de coordenadas (para orientar o plano da ST de acordo com a orientação de referência da superfície)
	(*Q0A)(0, 0) = dot(e1g, e1l);
	(*Q0A)(0, 1) = dot(e1g, e2l);
	(*Q0A)(0, 2) = dot(e1g, e3l);

	(*Q0A)(1, 0) = dot(e2g, e1l);
	(*Q0A)(1, 1) = dot(e2g, e2l);
	(*Q0A)(1, 2) = dot(e2g, e3l);

	(*Q0A)(2, 0) = dot(e3g, e1l);
	(*Q0A)(2, 1) = dot(e3g, e2l);
	(*Q0A)(2, 2) = dot(e3g, e3l);

	SaveConfiguration();

	DegenerationPreCalc();

}

//Retorna as coordenadas da superfície para um par (zeta,theta) - configuração anterior convergida
void FlexibleArcExtrusion_1::Gamma_and_Triad(Matrix* G_p, Matrix* t1_p, Matrix* t2_p, Matrix* n_p, Matrix* G_i, Matrix* t1_i, Matrix* t2_i, Matrix* n_i, Matrix* G_ip, double* zetai, double* thi, double* zetap, double* thp)
{

	double *d = d_A->getMatrix();		//ponteiro para o vetor d
	double* xAAi = x_AAi->getMatrix();
	double* xBAi = x_BAi->getMatrix();
	double** QAAi = aQ_AAi;
	double** QBAi = aQ_BAi;
	

	double *rad = radius;
	double *cpoint = c_point->getMatrix();
	

		
	double *Gi = G_i->getMatrix();		//ponteiro para o vetor Gi
	double *t1i = t1_i->getMatrix();	//ponteiro para o vetor t1i 
	double *t2i = t2_i->getMatrix();	//ponteiro para o vetor t2i
	double *ni = n_i->getMatrix();		//ponteiro para o vetor ni
	

	double *Gp = G_p->getMatrix();		//ponteiro para o vetor Gp 
	double *t1p = t1_p->getMatrix();	//ponteiro para o vetor t1p 
	double *t2p = t2_p->getMatrix();	//ponteiro para o vetor t2p
	double *np = n_p->getMatrix();		//ponteiro para o vetor np


	double *Gip = G_ip->getMatrix();	//ponteiro para o vetor Gip

	bool* normalint = &flag_normal_int;	//ponteiro para booleana flag_normal_int
	double v[1000];
	
	//AceGen
	//int b184, b225;
	v[274] = 0.5e0*(1e0 + (*zetap));
	v[273] = 0.5e0 - 0.5e0*(*zetap);
	v[272] = 0.5e0*(1e0 + (*zetai));
	v[271] = 0.5e0 - 0.5e0*(*zetai);
	v[270] = d[8] + xBAi[2];
	v[269] = d[7] + xBAi[1];
	v[268] = d[6] + xBAi[0];
	v[267] = d[2] + xAAi[2];
	v[266] = d[1] + xAAi[1];
	v[265] = d[0] + xAAi[0];
	v[260] = (*rad)*sin((*thi));
	v[259] = (*rad)*cos((*thi));
	v[258] = (*rad)*sin((*thp));
	v[257] = (*rad)*cos((*thp));
	v[256] = Power(d[11], 2);
	v[255] = 0.5e0*d[11] * d[9];
	v[254] = 0.5e0*d[10];
	v[253] = Power(d[10], 2);
	v[263] = v[253] + v[256];
	v[252] = d[9] * v[254];
	v[251] = Power(d[9], 2);
	v[250] = Power(d[5], 2);
	v[249] = 0.5e0*d[3] * d[5];
	v[248] = 0.5e0*d[4];
	v[247] = Power(d[4], 2);
	v[261] = v[247] + v[250];
	v[246] = d[3] * v[248];
	v[245] = Power(d[3], 2);
	v[94] = d[5] * v[248];
	v[110] = d[11] * v[254];
	v[81] = 4e0 / (4e0 + v[245] + v[261]);
	v[262] = -0.5e0*v[81];
	v[84] = 1e0 + v[261] * v[262];
	v[85] = (-d[5] + v[246])*v[81];
	v[86] = (d[4] + v[249])*v[81];
	v[88] = (d[5] + v[246])*v[81];
	v[90] = 1e0 + (v[245] + v[250])*v[262];
	v[91] = v[81] * (-d[3] + v[94]);
	v[93] = (-d[4] + v[249])*v[81];
	v[95] = v[81] * (d[3] + v[94]);
	v[96] = 1e0 + (v[245] + v[247])*v[262];
	v[97] = 4e0 / (4e0 + v[251] + v[263]);
	v[264] = -0.5e0*v[97];
	v[100] = 1e0 + v[263] * v[264];
	v[101] = (-d[11] + v[252])*v[97];
	v[102] = (d[10] + v[255])*v[97];
	v[104] = (d[11] + v[252])*v[97];
	v[106] = 1e0 + (v[251] + v[256])*v[264];
	v[107] = (-d[9] + v[110])*v[97];
	v[109] = (-d[10] + v[255])*v[97];
	v[111] = (d[9] + v[110])*v[97];
	v[112] = 1e0 + (v[251] + v[253])*v[264];
	v[113] = QAAi[0][0] * v[84] + QAAi[1][0] * v[85] + QAAi[2][0] * v[86];
	v[114] = QAAi[0][1] * v[84] + QAAi[1][1] * v[85] + QAAi[2][1] * v[86];
	v[116] = QAAi[0][0] * v[88] + QAAi[1][0] * v[90] + QAAi[2][0] * v[91];
	v[117] = QAAi[0][1] * v[88] + QAAi[1][1] * v[90] + QAAi[2][1] * v[91];
	v[119] = QAAi[0][0] * v[93] + QAAi[1][0] * v[95] + QAAi[2][0] * v[96];
	v[120] = QAAi[0][1] * v[93] + QAAi[1][1] * v[95] + QAAi[2][1] * v[96];
	v[122] = QBAi[0][0] * v[100] + QBAi[1][0] * v[101] + QBAi[2][0] * v[102];
	v[123] = QBAi[0][1] * v[100] + QBAi[1][1] * v[101] + QBAi[2][1] * v[102];
	v[125] = QBAi[0][0] * v[104] + QBAi[1][0] * v[106] + QBAi[2][0] * v[107];
	v[126] = QBAi[0][1] * v[104] + QBAi[1][1] * v[106] + QBAi[2][1] * v[107];
	v[128] = QBAi[0][0] * v[109] + QBAi[1][0] * v[111] + QBAi[2][0] * v[112];
	v[129] = QBAi[0][1] * v[109] + QBAi[1][1] * v[111] + QBAi[2][1] * v[112];
	v[141] = cpoint[0] + v[259];
	v[142] = cpoint[1] + v[260];
	v[161] = QBAi[2][0] * v[141] + QBAi[2][1] * v[142] + xBAi[2];
	v[160] = QAAi[2][0] * v[141] + QAAi[2][1] * v[142] + xAAi[2];
	v[158] = QBAi[1][0] * v[141] + QBAi[1][1] * v[142] + xBAi[1];
	v[157] = QAAi[1][0] * v[141] + QAAi[1][1] * v[142] + xAAi[1];
	v[155] = QBAi[0][0] * v[141] + QBAi[0][1] * v[142] + xBAi[0];
	v[154] = QAAi[0][0] * v[141] + QAAi[0][1] * v[142] + xAAi[0];
	v[143] = cpoint[0] + v[257];
	v[144] = cpoint[1] + v[258];
	v[202] = v[128] * v[143] + v[129] * v[144] + v[270];
	v[201] = v[119] * v[143] + v[120] * v[144] + v[267];
	v[199] = v[125] * v[143] + v[126] * v[144] + v[269];
	v[198] = v[116] * v[143] + v[117] * v[144] + v[266];
	v[196] = v[122] * v[143] + v[123] * v[144] + v[268];
	v[195] = v[113] * v[143] + v[114] * v[144] + v[265];
	v[156] = -0.5e0*v[154] + 0.5e0*v[155];
	v[159] = -0.5e0*v[157] + 0.5e0*v[158];
	v[162] = -0.5e0*v[160] + 0.5e0*v[161];
	v[164] = 1e0 / sqrt(Power(v[156], 2) + Power(v[159], 2) + Power(v[162], 2));
	v[163] = v[156] * v[164];
	v[165] = v[159] * v[164];
	v[166] = v[162] * v[164];
	v[175] = (QAAi[0][1] * v[259] - QAAi[0][0] * v[260])*v[271] + (QBAi[0][1] * v[259] - QBAi[0][0] * v[260])*v[272];
	v[176] = (QAAi[1][1] * v[259] - QAAi[1][0] * v[260])*v[271] + (QBAi[1][1] * v[259] - QBAi[1][0] * v[260])*v[272];
	v[177] = (QAAi[2][1] * v[259] - QAAi[2][0] * v[260])*v[271] + (QBAi[2][1] * v[259] - QBAi[2][0] * v[260])*v[272];
	v[181] = 1e0 / sqrt(Power(v[175], 2) + Power(v[176], 2) + Power(v[177], 2));
	v[180] = v[175] * v[181];
	v[182] = v[176] * v[181];
	v[275] = -(v[165] * v[180]) + v[163] * v[182];
	v[183] = v[177] * v[181];
	v[278] = -(v[166] * v[182]) + v[165] * v[183];
	v[277] = v[166] * v[180] - v[163] * v[183];
	v[276] = 1e0 / sqrt(Power(v[275], 2) + Power(v[277], 2) + Power(v[278], 2));
	v[194] = v[275] * v[276];
	v[193] = v[276] * v[277];
	v[192] = v[276] * v[278];
	if ((*normalint) == 1){
		v[186] = v[192];
		v[189] = v[193];
		v[191] = v[194];
	}
	else {
		v[186] = -v[192];
		v[189] = -v[193];
		v[191] = -v[194];
	};
	v[197] = -0.5e0*v[195] + 0.5e0*v[196];
	v[200] = -0.5e0*v[198] + 0.5e0*v[199];
	v[203] = -0.5e0*v[201] + 0.5e0*v[202];
	v[205] = 1e0 / sqrt(Power(v[197], 2) + Power(v[200], 2) + Power(v[203], 2));
	v[204] = v[197] * v[205];
	v[206] = v[200] * v[205];
	v[207] = v[203] * v[205];
	v[216] = (v[114] * v[257] - v[113] * v[258])*v[273] + (v[123] * v[257] - v[122] * v[258])*v[274];
	v[217] = (v[117] * v[257] - v[116] * v[258])*v[273] + (v[126] * v[257] - v[125] * v[258])*v[274];
	v[218] = (v[120] * v[257] - v[119] * v[258])*v[273] + (v[129] * v[257] - v[128] * v[258])*v[274];
	v[222] = 1e0 / sqrt(Power(v[216], 2) + Power(v[217], 2) + Power(v[218], 2));
	v[221] = v[216] * v[222];
	v[223] = v[217] * v[222];
	v[279] = -(v[206] * v[221]) + v[204] * v[223];
	v[224] = v[218] * v[222];
	v[282] = -(v[207] * v[223]) + v[206] * v[224];
	v[281] = v[207] * v[221] - v[204] * v[224];
	v[280] = 1e0 / sqrt(Power(v[279], 2) + Power(v[281], 2) + Power(v[282], 2));
	v[235] = v[279] * v[280];
	v[234] = v[280] * v[281];
	v[233] = v[280] * v[282];
	if ((*normalint) == 1){
		v[227] = v[233];
		v[230] = v[234];
		v[232] = v[235];
	}
	else {
		v[227] = -v[233];
		v[230] = -v[234];
		v[232] = -v[235];
	};
	Gip[0] = (v[113] * v[141] + v[114] * v[142] + v[265])*v[271] + (v[122] * v[141] + v[123] * v[142] + v[268])*v[272];
	Gip[1] = (v[116] * v[141] + v[117] * v[142] + v[266])*v[271] + (v[125] * v[141] + v[126] * v[142] + v[269])*v[272];
	Gip[2] = (v[119] * v[141] + v[120] * v[142] + v[267])*v[271] + (v[128] * v[141] + v[129] * v[142] + v[270])*v[272];
	t1i[0] = v[163];
	t1i[1] = v[165];
	t1i[2] = v[166];
	t2i[0] = v[180];
	t2i[1] = v[182];
	t2i[2] = v[183];
	ni[0] = v[186];
	ni[1] = v[189];
	ni[2] = v[191];
	Gp[0] = v[195] * v[273] + v[196] * v[274];
	Gp[1] = v[198] * v[273] + v[199] * v[274];
	Gp[2] = v[201] * v[273] + v[202] * v[274];
	t1p[0] = v[204];
	t1p[1] = v[206];
	t1p[2] = v[207];
	t2p[0] = v[221];
	t2p[1] = v[223];
	t2p[2] = v[224];
	np[0] = v[227];
	np[1] = v[230];
	np[2] = v[232];
	Gi[0] = v[154] * v[271] + v[155] * v[272];
	Gi[1] = v[157] * v[271] + v[158] * v[272];
	Gi[2] = v[160] * v[271] + v[161] * v[272];
	
}

//Dado o ponto xS, calcula as coordenadas (zeta,theta) referentes à mínima distância
void FlexibleArcExtrusion_1::FindMinimimumParameters(Matrix* xS, NSContactData* cd)
{

}

//Atualiza as variáveis internas da superfície, para pegarem info do pilot node para uso posterior com posição atualizada
void FlexibleArcExtrusion_1::FillNodes()
{
	Matrix alphaAA(3);
	Matrix alphaBA(3);
	Matrix uA(3);
	Matrix uB(3);

	for (int i = 0; i < 3; i++)
	{
		(*d_A)(i + 0, 0) = db.nodes[nodes[0] - 1]->displacements[i];
		(*d_A)(i + 3, 0) = db.nodes[nodes[0] - 1]->displacements[i + 3];
		(*d_A)(i + 6, 0) = db.nodes[nodes[1] - 1]->displacements[i];
		(*d_A)(i + 9, 0) = db.nodes[nodes[1] - 1]->displacements[i + 3];

		uA(i, 0) = db.nodes[nodes[0] - 1]->displacements[i];
		uB(i, 0) = db.nodes[nodes[1] - 1]->displacements[i];
		alphaAA(i, 0) = db.nodes[nodes[0] - 1]->displacements[i + 3];
		alphaBA(i, 0) = db.nodes[nodes[1] - 1]->displacements[i + 3];
	}

	//Q_AAp
	double alpha = norm(alphaAA);							//Valor escalar do parametro alpha
	Matrix A = skew(alphaAA);								//Matriz A
	double g = 4.0 / (4.0 + alpha * alpha);					//função g(alpha) - em algumas ref. tb. chamado de h(alpha)
	Matrix Qd = *I3 + g * (A + 0.5*(A*A));					//Tensor de rotação
	*Q_AAp = Qd * (*Q_AAi);
	//Q_BAp
	alpha = norm(alphaBA);							//Valor escalar do parametro alpha
	A = skew(alphaBA);								//Matriz A
	g = 4.0 / (4.0 + alpha * alpha);				//função g(alpha) - em algumas ref. tb. chamado de h(alpha)
	Qd = *I3 + g * (A + 0.5*(A*A));					//Tensor de rotação
	*Q_BAp = Qd * (*Q_BAi);
	//x_AAp
	*x_AAp = *x_AAi + uA;
	//x_BAp
	*x_BAp = *x_BAi + uB;
}

//Retorna coordenadas globais do ponto central da superfície a ser utilizado para cálculos grosseiros de sua localização (pinball)
void FlexibleArcExtrusion_1::CenterPoint(Matrix* center)
{
	*center = 0.5*(*x_AAi + *x_BAi);
}

//Salva vetores de configuração convergida
void FlexibleArcExtrusion_1::SaveConfiguration()
{
	d_A->clear();
	for (int i = 0; i < 3; i++)
	{
		(*x_AAi)(i, 0) = db.nodes[nodes[0] - 1]->copy_coordinates[i];
		(*x_BAi)(i, 0) = db.nodes[nodes[1] - 1]->copy_coordinates[i];
		(*alpha_AA)(i, 0) = db.nodes[nodes[0] - 1]->copy_coordinates[i + 3];
		(*alpha_BA)(i, 0) = db.nodes[nodes[1] - 1]->copy_coordinates[i + 3];

		(*dui_A)(i + 0, 0) = db.nodes[nodes[0] - 1]->copy_vel[i];
		(*dui_A)(i + 3, 0) = db.nodes[nodes[0] - 1]->copy_vel[i + 3];
		(*dui_A)(i + 6, 0) = db.nodes[nodes[1] - 1]->copy_vel[i];
		(*dui_A)(i + 9, 0) = db.nodes[nodes[1] - 1]->copy_vel[i + 3];

		(*ddui_A)(i + 0, 0) = db.nodes[nodes[0] - 1]->copy_accel[i];
		(*ddui_A)(i + 3, 0) = db.nodes[nodes[0] - 1]->copy_accel[i + 3];
		(*ddui_A)(i + 6, 0) = db.nodes[nodes[1] - 1]->copy_accel[i];
		(*ddui_A)(i + 9, 0) = db.nodes[nodes[1] - 1]->copy_accel[i + 3];
	}
	
	//Q_AAi
	double alpha = norm(*alpha_AA);							//Valor escalar do parametro alpha
	Matrix A = skew(*alpha_AA);								//Matriz A
	double g = 4.0 / (4.0 + alpha*alpha);					//função g(alpha) - em algumas ref. tb. chamado de h(alpha)
	*Q_AAi = *I3 + g*(A + 0.5*(A*A));						//Tensor de rotação
	*Q_AAic = *Q_AAi;										//Cópia - antes da transformação embutida
	*Q_AAi = (*Q_AAi)*(*Q0A);								//Para já realizar a transformação da parametrização do arco para o sistema global
	Q_AAi->MatrixToPtr(aQ_AAi, 3);
	//Q_BAi
	alpha = norm(*alpha_BA);								//Valor escalar do parametro alpha
	A = skew(*alpha_BA);									//Matriz A
	g = 4.0 / (4.0 + alpha*alpha);							//função g(alpha) - em algumas ref. tb. chamado de h(alpha)
	*Q_BAi = *I3 + g*(A + 0.5*(A*A));						//Tensor de rotação
	*Q_BAic = *Q_BAi;										//Cópia - antes da transformação embutida
	*Q_BAi = (*Q_BAi)*(*Q0A);								//Para já realizar a transformação da parametrização do arco para o sistema global
	Q_BAi->MatrixToPtr(aQ_BAi, 3);
}

void FlexibleArcExtrusion_1::UpdateBox()
{
	//Objetivo dessa funÁ„o È obter o minimo box que englobe o arco

	//Atribuindo valores do arco
	radius = &db.arcs[arc_ID - 1]->radius;
	c_point = &db.arcs[arc_ID - 1]->c_point;
	i_point = &db.arcs[arc_ID - 1]->i_point;
	f_point = &db.arcs[arc_ID - 1]->f_point;
	theta_i = &db.arcs[arc_ID - 1]->theta_i;
	theta_f = &db.arcs[arc_ID - 1]->theta_f;

	//Movendo para origem
	Matrix c_point0(2);
	c_point0(0, 0) = 0.0;
	c_point0(1, 0) = 0.0;
	Matrix i_point0 = *i_point - *c_point;
	Matrix f_point0 = *f_point - *c_point;

	//Localizando pontos que cruzam os quadrantes (m·ximo 4) - isso ser· feito considerando que o centro de curvatura est· na origem
	Matrix p1(2); //no eixo x positivo
	Matrix p2(2); //no eixo y positivo
	Matrix p3(2); //no eixo x negativo
	Matrix p4(2); //no eixo y negativo

	//Iniciando os pontos com o mesmo valor do initial point, pois assim mesmo que o ponto pn n„o exista, ele n„o ir· interferir na localizaÁ„o
	//de m·ximo e mÌnimos, pois o ponto i sempre ir· existir...
	p1 = i_point0;
	p2 = i_point0;
	p3 = i_point0;
	p4 = i_point0;

	//Criando thetas que v„o de 0 a 2PI, para viabilizar as verificaÁıes da exisitÍncia dos pontos p1-4
	double theta_imod;
	double theta_fmod;

	//Caso seja um retorno do terceiro ou quarto quadrante (negativo na funÁ„o atan2), somar 2PI para obter o ‚ngulo completo. Caso contr·tio, usar o ‚ngulo original.
	if (*theta_i < 0)// && *theta_i >= -PI)
		theta_imod = 2 * PI + *theta_i;
	else
		theta_imod = *theta_i;
	if (*theta_f < 0)// && *theta_f >= -PI)
		theta_fmod = 2 * PI + *theta_f;
	else
		theta_fmod = *theta_f;


	//Ponto inicial no quarto quadrante e final no primeiro
	if ((sin(theta_imod) < 0 && cos(theta_imod) > 0) && (sin(theta_fmod) > 0 && cos(theta_fmod) > 0))
	{
		p1(0, 0) = *radius;
		p1(1, 0) = 0.0;
	}
	//Ponto inicial no quarto quadrante e final no segundo
	if ((sin(theta_imod) < 0 && cos(theta_imod) > 0) && (sin(theta_fmod) > 0 && cos(theta_fmod) < 0))
	{
		p1(0, 0) = *radius;
		p1(1, 0) = 0.0;
		p2(0, 0) = 0.0;
		p2(1, 0) = *radius;
	}
	//Ponto inicial no quarto quadrante e final no terceiro
	if ((sin(theta_imod) < 0 && cos(theta_imod) > 0) && (sin(theta_fmod) < 0 && cos(theta_fmod) < 0))
	{
		p1(0, 0) = *radius;
		p1(1, 0) = 0.0;
		p2(0, 0) = 0.0;
		p2(1, 0) = *radius;
		p3(0, 0) = -(*radius);
		p3(1, 0) = 0.0;
	}
	//Ponto inicial no quarto quadrante e final no quarto
	if ((sin(theta_imod) < 0 && cos(theta_imod) > 0) && (sin(theta_fmod) < 0 && cos(theta_fmod) > 0) && (theta_fmod < theta_imod))
	{
		p1(0, 0) = *radius;
		p1(1, 0) = 0.0;
		p2(0, 0) = 0.0;
		p2(1, 0) = *radius;
		p3(0, 0) = -(*radius);
		p3(1, 0) = 0.0;
		p4(0, 0) = 0.0;
		p4(1, 0) = -(*radius);
	}

	//

	//Ponto inicial no primeiro quadrante e final no primeiro
	if ((sin(theta_imod) > 0 && cos(theta_imod) > 0) && (sin(theta_fmod) > 0 && cos(theta_fmod) > 0) && (theta_fmod < theta_imod))
	{
		p1(0, 0) = *radius;
		p1(1, 0) = 0.0;
		p2(0, 0) = 0.0;
		p2(1, 0) = *radius;
		p3(0, 0) = -(*radius);
		p3(1, 0) = 0.0;
		p4(0, 0) = 0.0;
		p4(1, 0) = -(*radius);
	}
	//Ponto inicial no primeiro quadrante e final no segundo
	if ((sin(theta_imod) > 0 && cos(theta_imod) > 0) && (sin(theta_fmod) > 0 && cos(theta_fmod) < 0))
	{
		p2(0, 0) = 0.0;
		p2(1, 0) = *radius;
	}
	//Ponto inicial no primeiro quadrante e final no terceiro
	if ((sin(theta_imod) > 0 && cos(theta_imod) > 0) && (sin(theta_fmod) < 0 && cos(theta_fmod) < 0))
	{
		p2(0, 0) = 0.0;
		p2(1, 0) = *radius;
		p3(0, 0) = -(*radius);
		p3(1, 0) = 0.0;
	}
	//Ponto inicial no primeiro quadrante e final no quarto
	if ((sin(theta_imod) > 0 && cos(theta_imod) > 0) && (sin(theta_fmod) < 0 && cos(theta_fmod) > 0))
	{
		p2(0, 0) = 0.0;
		p2(1, 0) = *radius;
		p3(0, 0) = -(*radius);
		p3(1, 0) = 0.0;
		p4(0, 0) = 0.0;
		p4(1, 0) = -(*radius);
	}

	//

	//Ponto inicial no segundo quadrante e final no primeiro
	if ((sin(theta_imod) > 0 && cos(theta_imod) < 0) && (sin(theta_fmod) > 0 && cos(theta_fmod) > 0))
	{
		p1(0, 0) = *radius;
		p1(1, 0) = 0.0;
		p3(0, 0) = -(*radius);
		p3(1, 0) = 0.0;
		p4(0, 0) = 0.0;
		p4(1, 0) = -(*radius);
	}
	//Ponto inicial no segundo quadrante e final no segundo
	if ((sin(theta_imod) > 0 && cos(theta_imod) < 0) && (sin(theta_fmod) > 0 && cos(theta_fmod) < 0) && (theta_fmod < theta_imod))
	{
		p1(0, 0) = *radius;
		p1(1, 0) = 0.0;
		p2(0, 0) = 0.0;
		p2(1, 0) = *radius;
		p3(0, 0) = -(*radius);
		p3(1, 0) = 0.0;
		p4(0, 0) = 0.0;
		p4(1, 0) = -(*radius);
	}
	//Ponto inicial no segundo quadrante e final no terceiro
	if ((sin(theta_imod) > 0 && cos(theta_imod) < 0) && (sin(theta_fmod) < 0 && cos(theta_fmod) < 0))
	{
		p3(0, 0) = -(*radius);
		p3(1, 0) = 0.0;
	}
	//Ponto inicial no segundo quadrante e final no quarto
	if ((sin(theta_imod) > 0 && cos(theta_imod) < 0) && (sin(theta_fmod) < 0 && cos(theta_fmod) > 0))
	{
		p3(0, 0) = -(*radius);
		p3(1, 0) = 0.0;
		p4(0, 0) = 0.0;
		p4(1, 0) = -(*radius);
	}

	//

	//Ponto inicial no terceiro quadrante e final no primeiro
	if ((sin(theta_imod) < 0 && cos(theta_imod) < 0) && (sin(theta_fmod) > 0 && cos(theta_fmod) > 0))
	{
		p1(0, 0) = *radius;
		p1(1, 0) = 0.0;
		p4(0, 0) = 0.0;
		p4(1, 0) = -(*radius);
	}
	//Ponto inicial no terceiro quadrante e final no segundo
	if ((sin(theta_imod) < 0 && cos(theta_imod) < 0) && (sin(theta_fmod) > 0 && cos(theta_fmod) < 0))
	{
		p1(0, 0) = *radius;
		p1(1, 0) = 0.0;
		p2(0, 0) = 0.0;
		p2(1, 0) = *radius;
		p4(0, 0) = 0.0;
		p4(1, 0) = -(*radius);
	}
	//Ponto inicial no terceiro quadrante e final no terceiro
	if ((sin(theta_imod) < 0 && cos(theta_imod) < 0) && (sin(theta_fmod) < 0 && cos(theta_fmod) < 0) && (theta_fmod < theta_imod))
	{
		p1(0, 0) = *radius;
		p1(1, 0) = 0.0;
		p2(0, 0) = 0.0;
		p2(1, 0) = *radius;
		p3(0, 0) = -(*radius);
		p3(1, 0) = 0.0;
		p4(0, 0) = 0.0;
		p4(1, 0) = -(*radius);
	}
	//Ponto inicial no terceiro quadrante e final no quarto
	if ((sin(theta_imod) < 0 && cos(theta_imod) < 0) && (sin(theta_fmod) < 0 && cos(theta_fmod) > 0))
	{
		p4(0, 0) = 0.0;
		p4(1, 0) = -(*radius);
	}


	//Organizando os pontos em uma tabela para comparaÁ„o (centro, inicial, final, p1, p2, p3 e p4)
	Matrix points(7, 2);
	points(0, 0) = i_point0(0, 0); //Mudei de ideia, melhor n„o incluir o center point na comparaÁ„o, pois a principio parece desnecess·rio
	points(0, 1) = i_point0(1, 0); //No entanto, se identificar que È melhor fazer o box incluindo o centro È sÛ trocar i_point0 por c_point0 nessas duas linhas

	points(1, 0) = i_point0(0, 0);
	points(1, 1) = i_point0(1, 0);

	points(2, 0) = f_point0(0, 0);
	points(2, 1) = f_point0(1, 0);

	points(3, 0) = p1(0, 0);
	points(3, 1) = p1(1, 0);

	points(4, 0) = p2(0, 0);
	points(4, 1) = p2(1, 0);

	points(5, 0) = p3(0, 0);
	points(5, 1) = p3(1, 0);

	points(6, 0) = p4(0, 0);
	points(6, 1) = p4(1, 0);

	//Encontrar o m·ximo x e y e o mÌnimo x e y
	Matrix maxcorner0(2);
	Matrix mincorner0(2);

	maxcorner0 = i_point0;
	mincorner0 = i_point0;

	for (int i = 0; i < 7; i++)
	{
		if (points(i, 0) > maxcorner0(0, 0))
			maxcorner0(0, 0) = points(i, 0);
		if (points(i, 1) > maxcorner0(1, 0))
			maxcorner0(1, 0) = points(i, 1);
		if (points(i, 0) < mincorner0(0, 0))
			mincorner0(0, 0) = points(i, 0);
		if (points(i, 1) < mincorner0(1, 0))
			mincorner0(1, 0) = points(i, 1);
	}

	//Definindo os quatro vertices do retangulo e contabilizando a posiÁ„o do centro de curvatura
	Matrix A0(2);
	Matrix B0(2);
	Matrix C0(2);
	Matrix D0(2);

	A0(0, 0) = maxcorner0(0, 0);
	A0(1, 0) = maxcorner0(1, 0);

	B0(0, 0) = mincorner0(0, 0);
	B0(1, 0) = maxcorner0(1, 0);

	C0(0, 0) = mincorner0(0, 0);
	C0(1, 0) = mincorner0(1, 0);

	D0(0, 0) = maxcorner0(0, 0);
	D0(1, 0) = mincorner0(1, 0);

	Matrix A = A0 + *c_point;
	Matrix B = B0 + *c_point;
	Matrix C = C0 + *c_point;
	Matrix D = D0 + *c_point;

	//Criando vetores com os corners do box para serem interpolados entre os nÛs da viga
	Matrix vA(3);
	vA(0, 0) = A(0, 0);
	vA(1, 0) = A(1, 0);
	vA(2, 0) = 0.0;

	Matrix vB(3);
	vB(0, 0) = B(0, 0);
	vB(1, 0) = B(1, 0);
	vB(2, 0) = 0.0;

	Matrix vC(3);
	vC(0, 0) = C(0, 0);
	vC(1, 0) = C(1, 0);
	vC(2, 0) = 0.0;

	Matrix vD(3);
	vD(0, 0) = D(0, 0);
	vD(1, 0) = D(1, 0);
	vD(2, 0) = 0.0;

	//InterpolaÁ„o
	Matrix xAAi(3);
	xAAi = *x_AAi;
	Matrix xBAi(3);
	xBAi = *x_BAi;
	Matrix QAAi(3, 3);
	QAAi = *Q_AAi;
	Matrix QBAi(3, 3);
	QBAi = *Q_BAi;

	double x[8], y[8], z[8];

	//NÛ A superfÌcie A
	Matrix v1(3);
	v1 = QAAi * vA + xAAi;
	x[0] = v1(0, 0);
	y[0] = v1(1, 0);
	z[0] = v1(2, 0);

	Matrix v2(3);
	v2 = QAAi * vB + xAAi;
	x[1] = v2(0, 0);
	y[1] = v2(1, 0);
	z[1] = v2(2, 0);

	Matrix v3(3);
	v3 = QAAi * vC + xAAi;
	x[2] = v3(0, 0);
	y[2] = v3(1, 0);
	z[2] = v3(2, 0);

	Matrix v4(3);
	v4 = QAAi * vD + xAAi;
	x[3] = v4(0, 0);
	y[3] = v4(1, 0);
	z[3] = v4(2, 0);

	//NÛ B superfÌcie A
	Matrix v5(3);
	v5 = QBAi * vA + xBAi;
	x[4] = v5(0, 0);
	y[4] = v5(1, 0);
	z[4] = v5(2, 0);

	Matrix v6(3);
	v6 = QBAi * vB + xBAi;
	x[5] = v6(0, 0);
	y[5] = v6(1, 0);
	z[5] = v6(2, 0);

	Matrix v7(3);
	v7 = QBAi * vC + xBAi;
	x[6] = v7(0, 0);
	y[6] = v7(1, 0);
	z[6] = v7(2, 0);

	Matrix v8(3);
	v8 = QBAi * vD + xBAi;
	x[7] = v8(0, 0);
	y[7] = v8(1, 0);
	z[7] = v8(2, 0);

	//Setando os vÈrtices
	box.SetVertices(x, y, z);
}
//Calcula contribuições de contato entre esfera e superfície
void FlexibleArcExtrusion_1::ContactSphereSurfaceSticking(double* Rc, double** Kc, double zetap, double thetap, double zetai, double thetai, double* gti, int node, double* epsn, double* epst, double* cn, double* ct, double* mu, double* radius)
{
	//TODO
}

//Calcula contribuições de contato entre esfera e superfície
void FlexibleArcExtrusion_1::ContactSphereSurfaceSliding(double* Rc, double** Kc, double zetap, double thetap, double zetai, double thetai, double* gti, int node, double* epsn, double* epst, double* cn, double* ct, double* mu, double* radius)
{
	//TODO
}

