#include "RigidArcRevolution_1.h"

#include "ArcCirc.h"
#include "PostFiles.h"
#include "Node.h"
#include "Encoding.h"
#include "CoordinateSystem.h"

#include"Database.h"
//Variáveis globais
extern
Database db;
#define PI 3.1415926535897932384626433832795

RigidArcRevolution_1::RigidArcRevolution_1()
{
	nDOFs = 6;
	n_nodes = 1;
	number = 0;
	cs = 0;
	nodes = new int[n_nodes];
	VTK_nodes = new int[n_nodes];
	DOFs = new int *[n_nodes];
	for (int i = 0; i < n_nodes; i++)
		DOFs[i] = new int[db.number_GLs_node];
	VTK_nodes[0] = 0;
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

	rev_Ang = 2 * PI;
	x_fac = 1;
	z_fac = 1;


	x_ABi = new Matrix(3);
	Q_ABi = new Matrix(3, 3);
	Q_ABic = new Matrix(3, 3);
	Q0B = new Matrix(3, 3);
	d_B = new Matrix(6);
	dui_B = new Matrix(6);
	ddui_B = new Matrix(6);
	alpha_AB = new Matrix(3);
	aQ_ABi = new double*[3];
	for (int i = 0; i < 3; i++)
		aQ_ABi[i] = new double[3];
	I3 = new Matrix(3, 3);
	(*I3)(0, 0) = 1.0;
	(*I3)(1, 1) = 1.0;
	(*I3)(2, 2) = 1.0;

	x_ABp = new Matrix(3);
	Q_ABp = new Matrix(3, 3);

	InitializeDegeneration();
}

RigidArcRevolution_1::~RigidArcRevolution_1()
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


	delete x_ABi;
	delete Q_ABi;
	delete Q_ABic;
	delete d_B;
	delete dui_B;
	delete ddui_B;
	delete alpha_AB;
	
	for (int i = 0; i < 3; i++)
		delete[]aQ_ABi[i];
	delete[]aQ_ABi;
	delete I3;
	delete Q0B;

	delete x_ABp;
	delete Q_ABp;
	FreeDegeneration();
}

void RigidArcRevolution_1::SetMinMaxRange()
{
	u1_min = 0;
	u1_max = rev_Ang;
	u2_min = db.arcs[arc_ID - 1]->theta_i;
	u2_max = db.arcs[arc_ID - 1]->theta_f;

	u1_range = rev_Ang;

	if (u2_max >= u2_min)
		u2_range = abs(u2_max - u2_min);
	else
		u2_range = 2 * PI - abs(u2_max - u2_min);
}

//Obtem ponto da superficie
void RigidArcRevolution_1::SurfacePoint(double& zeta, double& theta, Matrix& point)
{
	Matrix ar(3);
	ar(0, 0) = ((*radius) * cos(theta) + (*c_point)(0, 0))*x_fac*cos(zeta);
	ar(1, 0) = (*radius) * sin(theta) + (*c_point)(1, 0);
	ar(2, 0) = -((*radius) * cos(theta) + (*c_point)(0, 0))*z_fac*sin(zeta);
	//Gamma na posicao atual
	point = *x_ABp + (*Q_ABp) * ar;
}

//Normal exterior à superfície na posição escolhida
void RigidArcRevolution_1::NormalExt(double* zeta, double* theta, Matrix* n)
{
	double *d = d_B->getMatrix();		//ponteiro para o vetor d
	double* xABi = x_ABi->getMatrix();
	double** QABi = aQ_ABi;

	double *rad = radius;
	double *cpoint = c_point->getMatrix();
	double *xfac = &x_fac;
	double *zfac = &z_fac;

	double *next = n->getMatrix();		//ponteiro para o vetor n

	double* phi = zeta;

	bool* normalint = &flag_normal_int;	//ponteiro para booleana flag_normal_int
	double v[1000];
	
	v[108] = (*rad)*sin((*theta));
	v[107] = (*rad)*cos((*theta));
	v[106] = Power(d[5], 2);
	v[105] = 0.5e0*d[3] * d[5];
	v[104] = 0.5e0*d[4];
	v[103] = Power(d[4], 2);
	v[109] = v[103] + v[106];
	v[102] = d[3] * v[104];
	v[101] = Power(d[3], 2);
	v[100] = sin((*phi));
	v[99] = cos((*phi));
	v[52] = d[5] * v[104];
	v[76] = -(v[108] * v[99] * (*xfac));
	v[79] = v[100] * v[108] * (*zfac);
	v[69] = cpoint[0] + v[107];
	v[85] = -(v[69] * v[99] * (*zfac));
	v[83] = -(v[100] * v[69] * (*xfac));
	v[39] = 4e0 / (4e0 + v[101] + v[109]);
	v[110] = -0.5e0*v[39];
	v[42] = 1e0 + v[109] * v[110];
	v[43] = (-d[5] + v[102])*v[39];
	v[44] = (d[4] + v[105])*v[39];
	v[46] = (d[5] + v[102])*v[39];
	v[48] = 1e0 + (v[101] + v[106])*v[110];
	v[49] = v[39] * (-d[3] + v[52]);
	v[51] = (-d[4] + v[105])*v[39];
	v[53] = v[39] * (d[3] + v[52]);
	v[54] = 1e0 + (v[101] + v[103])*v[110];
	v[55] = QABi[0][0] * v[42] + QABi[1][0] * v[43] + QABi[2][0] * v[44];
	v[57] = QABi[0][2] * v[42] + QABi[1][2] * v[43] + QABi[2][2] * v[44];
	v[58] = QABi[0][0] * v[46] + QABi[1][0] * v[48] + QABi[2][0] * v[49];
	v[60] = QABi[0][2] * v[46] + QABi[1][2] * v[48] + QABi[2][2] * v[49];
	v[61] = QABi[0][0] * v[51] + QABi[1][0] * v[53] + QABi[2][0] * v[54];
	v[63] = QABi[0][2] * v[51] + QABi[1][2] * v[53] + QABi[2][2] * v[54];
	v[80] = v[107] * (QABi[0][1] * v[42] + QABi[1][1] * v[43] + QABi[2][1] * v[44]) + v[55] * v[76] + v[57] * v[79];
	v[81] = v[107] * (QABi[0][1] * v[46] + QABi[1][1] * v[48] + QABi[2][1] * v[49]) + v[58] * v[76] + v[60] * v[79];
	v[82] = v[107] * (QABi[0][1] * v[51] + QABi[1][1] * v[53] + QABi[2][1] * v[54]) + v[61] * v[76] + v[63] * v[79];
	v[86] = v[55] * v[83] + v[57] * v[85];
	v[87] = v[58] * v[83] + v[60] * v[85];
	v[96] = v[81] * v[86] - v[80] * v[87];
	v[88] = v[61] * v[83] + v[63] * v[85];
	v[94] = -(v[82] * v[86]) + v[80] * v[88];
	if ((*normalint)) {
		v[90] = -1e0;
	}
	else {
		v[90] = 1e0;
	};
	v[91] = v[82] * v[87] - v[81] * v[88];
	v[93] = 1e0 / sqrt(Power(v[91], 2) + Power(v[94], 2) + Power(v[96], 2));
	v[111] = v[90] * v[93];
	next[0] = v[111] * v[91];
	next[1] = v[111] * v[94];
	next[2] = v[111] * v[96];
}

void RigidArcRevolution_1::WriteVTK_XMLRender(FILE *f)
{
	if (db.post_files->WriteRigidContactSurfaces_flag == true)
	{
		int n_circ1 = 24; //número de divisões ao longo da revolução
		int n_circ2 = 24; //número de divisões ao longo do arco
		//Número de pontos a serem gerados
		int n_points = (2*n_circ1 + 1) * (2*n_circ2 + 1);
		//Número de células a serem geradas
		int n_cells = n_circ1*n_circ2;
		double** points;
		points = new double*[n_points];
		for (int i = 0; i < n_points; i++)
			points[i] = new double[3];
		int index = 0;
		double phi = 0;
		double phi_f = this->rev_Ang;	//ângulo de revolução
		double theta = 0;

		//
		double theta_imod;
		double theta_fmod;
		//

		double a_ellipse = x_fac;
		double b_ellipse = z_fac;
		for (int m = 0; m < 2*n_circ1+1; m++)
		{
			if ((((*theta_i) >= PI / 2 && (*theta_i) <= PI)) && (((*theta_f) <= -PI / 2 && (*theta_f) >= -PI)) || (((*theta_f) >= PI / 2 && (*theta_f) <= PI)) && (((*theta_i) <= -PI / 2 && (*theta_i) >= -PI)) || (((*theta_i <= -PI / 2) && (*theta_i <= -PI)) && ((*theta_f <= -PI / 2) && (*theta_f <= -PI))))
			{
				if (*theta_i < 0)
					theta_imod = 2 * PI + *theta_i;
				else
					theta_imod = *theta_i;
				if (*theta_f < 0)
					theta_fmod = 2 * PI + *theta_f;
				else
					theta_fmod = *theta_f;

				theta = theta_imod;


				for (int n = 0; n < 2 * n_circ2 + 1; n++)
				{
					points[index][0] = ((*radius)*cos(theta) + (*c_point)(0, 0))*a_ellipse*cos(phi);
					points[index][1] = (*radius)*sin(theta) + (*c_point)(1, 0);
					points[index][2] = -1.0*((*radius)*cos(theta) + (*c_point)(0, 0))*b_ellipse*sin(phi);

					index += 1;
					theta += (theta_fmod - theta_imod) / (2 * n_circ2);
				}

				phi += phi_f / (2 * n_circ1);
			}
			else
			{
				theta = *theta_i;
				for (int n = 0; n < 2 * n_circ2 + 1; n++)
				{
					points[index][0] = ((*radius)*cos(theta) + (*c_point)(0, 0))*a_ellipse*cos(phi);
					points[index][1] = (*radius)*sin(theta) + (*c_point)(1, 0);
					points[index][2] = -1.0*((*radius)*cos(theta) + (*c_point)(0, 0))*b_ellipse*sin(phi);

					index += 1;
					theta += (*theta_f - *theta_i) / (2 * n_circ2);
				}

				phi += phi_f / (2 * n_circ1);
			}


			////Se os ângulos inicial e final estiverem no primeiro ou no quarto quadrante
			//if (((*theta_i < 0 && *theta_i >= -PI / 2) || (*theta_i >= 0 && *theta_i <= PI / 2)) && ((*theta_f < 0 && *theta_f >= -PI / 2) || (*theta_f >= 0 && *theta_f <= PI / 2)))
			//{
			//	theta = *theta_i;
			//	for (int n = 0; n < 2*n_circ2+1; n++)
			//	{
			//		points[index][0] = ((*radius)*cos(theta) + (*c_point)(0, 0))*a_ellipse*cos(phi);
			//		points[index][1] = (*radius)*sin(theta) + (*c_point)(1, 0);
			//		points[index][2] = -1.0*((*radius)*cos(theta) + (*c_point)(0, 0))*b_ellipse*sin(phi);

			//		index += 1;
			//		theta += (*theta_f - *theta_i) / (2*n_circ2);
			//	}

			//	phi += phi_f / (2*n_circ1);

			//}
			//
			////Se os ângulos inicial e final esiverem no segundo ou terceiro quadrante
			//if (((*theta_i < -PI / 2 && *theta_i >= -PI) || (*theta_i >= PI / 2 && *theta_i <= PI)) && ((*theta_f < -PI / 2 && *theta_f >= -PI) || (*theta_f >= PI / 2 && *theta_f <= PI)))
			//{
			//	if (*theta_i<0)
			//		theta_imod = 2 * PI + *theta_i;
			//	else
			//		theta_imod = *theta_i;
			//	if (*theta_f<0)
			//		theta_fmod = 2 * PI + *theta_f;
			//	else
			//		theta_fmod = *theta_f;

			//	theta = theta_imod;

			//	for (int n = 0; n < 2 * n_circ2 + 1; n++)
			//	{
			//		points[index][0] = ((*radius)*cos(theta) + (*c_point)(0, 0))*a_ellipse*cos(phi);
			//		points[index][1] = (*radius)*sin(theta) + (*c_point)(1, 0);
			//		points[index][2] = -1.0*((*radius)*cos(theta) + (*c_point)(0, 0))*b_ellipse*sin(phi);

			//		index += 1;
			//		theta += (theta_fmod - theta_imod) / (2 * n_circ2);
			//	}

			//	phi += phi_f / (2 * n_circ1);
			//}
						

		}
		//vetores para escrita no formato binário - usando a função 'enconde'
		std::vector<float> float_vector;
		std::vector<int> int_vector;

		Matrix vec_P(3);

		//Opens Piece
		fprintf(f, "\t\t<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", n_points, n_cells);
		//Opens Points
		fprintf(f, "\t\t\t<Points>\n");
		//Opens DataArray
		fprintf(f, "\t\t\t\t<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"binary\">\n");
		float_vector.clear();
		//Preenchendo as coordenadas dos pontos
		Matrix xO(3);
		xO(0, 0) = db.nodes[nodes[0] - 1]->copy_coordinates[0];
		xO(1, 0) = db.nodes[nodes[0] - 1]->copy_coordinates[1];
		xO(2, 0) = db.nodes[nodes[0] - 1]->copy_coordinates[2];
		for (int point = 0; point <n_points; point++)//Percorre os nós que descrevem o perímetro da ST
		{
			//Posição de cada ponto P no plano xy (referência)
			vec_P(0, 0) = points[point][0];
			vec_P(1, 0) = points[point][1];// -2 * (*c_point)(1, 0);
			vec_P(2, 0) = points[point][2];
			vec_P = xO + (*Q_ABi)*vec_P;//Operando rotacionando para o sistema da barra e translação

			float_vector.push_back((float)(vec_P(0, 0)));
			float_vector.push_back((float)(vec_P(1, 0)));
			float_vector.push_back((float)(vec_P(2, 0)));

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
		//Linhas
		int nj = (2 * n_circ2 + 1);
		int_vector.clear();
		for (int i = 0; i < (2 * n_circ1 - 1); i = i + 2)
			for (int j = 0; j < (2 * n_circ2 - 1); j = j + 2)
			{
				/*fprintf(f, "\t\t\t\t\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
				j + 0 + (i + 0)*nj,
				j + 0 + (i + 2)*nj,
				j + 2 + (i + 2)*nj,
				j + 2 + (i + 0)*nj,
				j + 0 + (i + 1)*nj,
				j + 1 + (i + 2)*nj,
				j + 2 + (i + 1)*nj,
				j + 1 + (i + 0)*nj);*/
				int_vector.push_back(j + 0 + (i + 0)*nj);
				int_vector.push_back(j + 0 + (i + 2)*nj);
				int_vector.push_back(j + 2 + (i + 2)*nj);
				int_vector.push_back(j + 2 + (i + 0)*nj);
				int_vector.push_back(j + 0 + (i + 1)*nj);
				int_vector.push_back(j + 1 + (i + 2)*nj);
				int_vector.push_back(j + 2 + (i + 1)*nj);
				int_vector.push_back(j + 1 + (i + 0)*nj);
			}
		fprintf(f, encodeData(int_vector).c_str());
		fprintf(f, "\n");
		//Closes DataArray
		fprintf(f, "\t\t\t\t</DataArray>\n");
		//Opens DataArray
		fprintf(f, "\t\t\t\t<DataArray type=\"Int32\" Name=\"types\" format=\"binary\">\n");
		int_vector.clear();
		for (int i = 0; i < n_cells; i++)
		{
			int_vector.push_back(23);
			//fprintf(f, "\t\t\t\t\t%d\n", 23);
		}
		fprintf(f, encodeData(int_vector).c_str());
		fprintf(f, "\n");
		//Closes DataArray
		fprintf(f, "\t\t\t\t</DataArray>\n");
		//Opens DataArray
		fprintf(f, "\t\t\t\t<DataArray type=\"Int32\" Name=\"offsets\" format=\"binary\">\n");
		int_vector.clear();
		for (int i = 0; i < n_cells; i++)
		{
			int_vector.push_back((i + 1) * 8);
			//fprintf(f, "\t\t\t\t\t%d\n", (i + 1) * 8);
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
		for (int i = 0; i < n_points; i++)
			delete[] points[i];
		delete[]points;
	}
}

bool RigidArcRevolution_1::Read(FILE *f)
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
	if (!strcmp(s, "Node"))
	{
		fscanf(f, "%s", s);
		nodes[0] = atoi(s);
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

	fpos_t pos1;
	fgetpos(f, &pos1);
	fscanf(f, "%s", s);
	if (!strcmp(s, "RevolutionAngle"))
	{
		fscanf(f, "%s", s);
		rev_Ang = atof(s);
	}
	else
		fsetpos(f, &pos1);

	fpos_t pos2;
	fgetpos(f, &pos2);
	fscanf(f, "%s", s);
	if (!strcmp(s, "FactorX"))
	{
		fscanf(f, "%s", s);
		x_fac = atof(s);
	}
	else
		fsetpos(f, &pos2);

	fpos_t pos3;
	fgetpos(f, &pos3);
	fscanf(f, "%s", s);
	if (!strcmp(s, "FactorZ"))
	{
		fscanf(f, "%s", s);
		z_fac = atof(s);
	}
	else
		fsetpos(f, &pos3);

	ReadCommon(f);

	return true;
}

void RigidArcRevolution_1::Write(FILE *f)
{
	char s[20];
	if (flag_normal_int == true)
		sprintf(s, "Concave");
	else
		sprintf(s, "Convex");

	fprintf(f, "RigidArcRevolution_1\t%d\tArc\t%d\tCS\t%d\tNode\t%d\t%s\tRevolutionAngle\t%.6e\tFactorX\t%.6e\tFactorZ\t%.6e\n",
		number,
		arc_ID,
		cs,
		nodes[0],
		s,
		rev_Ang,
		x_fac,
		z_fac);
}

//Checa inconsistências para evitar erros de execução
bool RigidArcRevolution_1::Check()
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
void RigidArcRevolution_1::InitialGuess(Matrix* xS, double** convective, int n_solutions)
{

}

void RigidArcRevolution_1::PreCalc()
{
	//Atribuindo valor do raio, centro de curvatura do arco e fatores de revolução (ponteiros)
	radius = &db.arcs[arc_ID - 1]->radius;
	c_point = &db.arcs[arc_ID - 1]->c_point;
	i_point = &db.arcs[arc_ID - 1]->i_point;
	f_point = &db.arcs[arc_ID - 1]->f_point;
	theta_i = &db.arcs[arc_ID - 1]->theta_i;
	theta_f = &db.arcs[arc_ID - 1]->theta_f;
	
	//Apontando para posição que indica valor dos GLs globais
	for (int i = 0; i < 6; i++)
		GLs[i] = &db.nodes[nodes[0] - 1]->GLs[i];

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
	(*Q0B)(0, 0) = dot(e1g, e1l);
	(*Q0B)(0, 1) = dot(e1g, e2l);
	(*Q0B)(0, 2) = dot(e1g, e3l);

	(*Q0B)(1, 0) = dot(e2g, e1l);
	(*Q0B)(1, 1) = dot(e2g, e2l);
	(*Q0B)(1, 2) = dot(e2g, e3l);

	(*Q0B)(2, 0) = dot(e3g, e1l);
	(*Q0B)(2, 1) = dot(e3g, e2l);
	(*Q0B)(2, 2) = dot(e3g, e3l);

	SaveConfiguration();

	DegenerationPreCalc();

}

//Retorna as coordenadas da superfície para um par (zeta,theta) - configuração anterior convergida
void RigidArcRevolution_1::Gamma_and_Triad(Matrix* G_p, Matrix* t1_p, Matrix* t2_p, Matrix* n_p, Matrix* G_i, Matrix* t1_i, Matrix* t2_i, Matrix* n_i, Matrix* G_ip, double* phii, double* thi, double* phip, double* thp)
{
	double *d = d_B->getMatrix();		//ponteiro para o vetor d
	double* xABi = x_ABi->getMatrix();
	double** QABi = aQ_ABi;

	double *rad = radius;
	double *cpoint = c_point->getMatrix();
	double *xfac = &x_fac;
	double *zfac = &z_fac;


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
	//int b136, b173;
	v[211] = d[2] + xABi[2];
	v[210] = d[1] + xABi[1];
	v[209] = d[0] + xABi[0];
	v[206] = (*rad)*sin((*thi));
	v[205] = (*rad)*cos((*thi));
	v[204] = (*rad)*sin((*thp));
	v[203] = (*rad)*cos((*thp));
	v[202] = Power(d[5], 2);
	v[201] = 0.5e0*d[3] * d[5];
	v[200] = 0.5e0*d[4];
	v[199] = Power(d[4], 2);
	v[207] = v[199] + v[202];
	v[198] = d[3] * v[200];
	v[197] = Power(d[3], 2);
	v[196] = sin((*phip));
	v[195] = cos((*phip));
	v[194] = sin((*phii));
	v[193] = cos((*phii));
	v[78] = d[5] * v[200];
	v[161] = -(v[195] * v[204] * (*xfac));
	v[124] = -(v[193] * v[206] * (*xfac));
	v[163] = v[196] * v[204] * (*zfac);
	v[126] = v[194] * v[206] * (*zfac);
	v[99] = cpoint[0] + v[203];
	v[215] = -(v[99] * (*zfac));
	v[214] = v[99] * (*xfac);
	v[150] = v[195] * v[215];
	v[148] = -(v[196] * v[214]);
	v[95] = cpoint[0] + v[205];
	v[213] = -(v[95] * (*zfac));
	v[212] = v[95] * (*xfac);
	v[113] = v[193] * v[213];
	v[111] = -(v[194] * v[212]);
	v[65] = 4e0 / (4e0 + v[197] + v[207]);
	v[208] = -0.5e0*v[65];
	v[68] = 1e0 + v[207] * v[208];
	v[69] = (-d[5] + v[198])*v[65];
	v[70] = (d[4] + v[201])*v[65];
	v[72] = (d[5] + v[198])*v[65];
	v[74] = 1e0 + (v[197] + v[202])*v[208];
	v[75] = v[65] * (-d[3] + v[78]);
	v[77] = (-d[4] + v[201])*v[65];
	v[79] = v[65] * (d[3] + v[78]);
	v[80] = 1e0 + (v[197] + v[199])*v[208];
	v[81] = QABi[0][0] * v[68] + QABi[1][0] * v[69] + QABi[2][0] * v[70];
	v[82] = QABi[0][1] * v[68] + QABi[1][1] * v[69] + QABi[2][1] * v[70];
	v[83] = QABi[0][2] * v[68] + QABi[1][2] * v[69] + QABi[2][2] * v[70];
	v[84] = QABi[0][0] * v[72] + QABi[1][0] * v[74] + QABi[2][0] * v[75];
	v[85] = QABi[0][1] * v[72] + QABi[1][1] * v[74] + QABi[2][1] * v[75];
	v[86] = QABi[0][2] * v[72] + QABi[1][2] * v[74] + QABi[2][2] * v[75];
	v[87] = QABi[0][0] * v[77] + QABi[1][0] * v[79] + QABi[2][0] * v[80];
	v[88] = QABi[0][1] * v[77] + QABi[1][1] * v[79] + QABi[2][1] * v[80];
	v[89] = QABi[0][2] * v[77] + QABi[1][2] * v[79] + QABi[2][2] * v[80];
	v[93] = v[193] * v[212];
	v[94] = cpoint[1] + v[206];
	v[96] = v[194] * v[213];
	v[97] = v[195] * v[214];
	v[98] = cpoint[1] + v[204];
	v[100] = v[196] * v[215];
	v[114] = QABi[0][0] * v[111] + QABi[0][2] * v[113];
	v[115] = QABi[1][0] * v[111] + QABi[1][2] * v[113];
	v[116] = QABi[2][0] * v[111] + QABi[2][2] * v[113];
	v[120] = 1e0 / sqrt(Power(v[114], 2) + Power(v[115], 2) + Power(v[116], 2));
	v[119] = v[114] * v[120];
	v[121] = v[115] * v[120];
	v[122] = v[116] * v[120];
	v[127] = QABi[0][0] * v[124] + QABi[0][2] * v[126] + QABi[0][1] * v[205];
	v[128] = QABi[1][0] * v[124] + QABi[1][2] * v[126] + QABi[1][1] * v[205];
	v[129] = QABi[2][0] * v[124] + QABi[2][2] * v[126] + QABi[2][1] * v[205];
	v[133] = 1e0 / sqrt(Power(v[127], 2) + Power(v[128], 2) + Power(v[129], 2));
	v[132] = v[127] * v[133];
	v[134] = v[128] * v[133];
	v[216] = -(v[121] * v[132]) + v[119] * v[134];
	v[135] = v[129] * v[133];
	v[219] = -(v[122] * v[134]) + v[121] * v[135];
	v[218] = v[122] * v[132] - v[119] * v[135];
	v[217] = 1e0 / sqrt(Power(v[216], 2) + Power(v[218], 2) + Power(v[219], 2));
	v[146] = -(v[216] * v[217]);
	v[145] = -(v[217] * v[218]);
	v[144] = -(v[217] * v[219]);
	if ((*normalint) == 1){
		v[138] = v[144];
		v[141] = v[145];
		v[143] = v[146];
	}
	else {
		v[138] = -v[144];
		v[141] = -v[145];
		v[143] = -v[146];
	};
	v[151] = v[148] * v[81] + v[150] * v[83];
	v[152] = v[148] * v[84] + v[150] * v[86];
	v[153] = v[148] * v[87] + v[150] * v[89];
	v[157] = 1e0 / sqrt(Power(v[151], 2) + Power(v[152], 2) + Power(v[153], 2));
	v[156] = v[151] * v[157];
	v[158] = v[152] * v[157];
	v[159] = v[153] * v[157];
	v[164] = v[161] * v[81] + v[203] * v[82] + v[163] * v[83];
	v[165] = v[161] * v[84] + v[203] * v[85] + v[163] * v[86];
	v[166] = v[161] * v[87] + v[203] * v[88] + v[163] * v[89];
	v[170] = 1e0 / sqrt(Power(v[164], 2) + Power(v[165], 2) + Power(v[166], 2));
	v[169] = v[164] * v[170];
	v[171] = v[165] * v[170];
	v[220] = -(v[158] * v[169]) + v[156] * v[171];
	v[172] = v[166] * v[170];
	v[223] = -(v[159] * v[171]) + v[158] * v[172];
	v[222] = v[159] * v[169] - v[156] * v[172];
	v[221] = 1e0 / sqrt(Power(v[220], 2) + Power(v[222], 2) + Power(v[223], 2));
	v[183] = -(v[220] * v[221]);
	v[182] = -(v[221] * v[222]);
	v[181] = -(v[221] * v[223]);
	if ((*normalint) == 1){
		v[175] = v[181];
		v[178] = v[182];
		v[180] = v[183];
	}
	else {
		v[175] = -v[181];
		v[178] = -v[182];
		v[180] = -v[183];
	};
	Gip[0] = v[209] + v[81] * v[93] + v[82] * v[94] + v[83] * v[96];
	Gip[1] = v[210] + v[84] * v[93] + v[85] * v[94] + v[86] * v[96];
	Gip[2] = v[211] + v[87] * v[93] + v[88] * v[94] + v[89] * v[96];
	t1i[0] = v[119];
	t1i[1] = v[121];
	t1i[2] = v[122];
	t2i[0] = v[132];
	t2i[1] = v[134];
	t2i[2] = v[135];
	ni[0] = v[138];
	ni[1] = v[141];
	ni[2] = v[143];
	Gp[0] = v[209] + v[100] * v[83] + v[81] * v[97] + v[82] * v[98];
	Gp[1] = v[210] + v[100] * v[86] + v[84] * v[97] + v[85] * v[98];
	Gp[2] = v[211] + v[100] * v[89] + v[87] * v[97] + v[88] * v[98];
	t1p[0] = v[156];
	t1p[1] = v[158];
	t1p[2] = v[159];
	t2p[0] = v[169];
	t2p[1] = v[171];
	t2p[2] = v[172];
	np[0] = v[175];
	np[1] = v[178];
	np[2] = v[180];
	Gi[0] = QABi[0][0] * v[93] + QABi[0][1] * v[94] + QABi[0][2] * v[96] + xABi[0];
	Gi[1] = QABi[1][0] * v[93] + QABi[1][1] * v[94] + QABi[1][2] * v[96] + xABi[1];
	Gi[2] = QABi[2][0] * v[93] + QABi[2][1] * v[94] + QABi[2][2] * v[96] + xABi[2];
		
}

//Dado o ponto xS, calcula as coordenadas (zeta,theta) referentes à mínima distância
void RigidArcRevolution_1::FindMinimimumParameters(Matrix* xS, NSContactData* cd)
{

}
//Atualiza as variáveis internas da superfície, para pegarem info do pilot node para uso posterior com posição atualizada
void RigidArcRevolution_1::FillNodes()
{
	Matrix alphaB(3);
	Matrix uB(3);

	for (int i = 0; i < 3; i++)
	{
		(*d_B)(i + 0, 0) = db.nodes[nodes[0] - 1]->displacements[i];
		(*d_B)(i + 3, 0) = db.nodes[nodes[0] - 1]->displacements[i + 3];
		
		uB(i, 0) = db.nodes[nodes[0] - 1]->displacements[i];
		alphaB(i, 0) = db.nodes[nodes[0] - 1]->displacements[i + 3];
	}

	//Q_ABp
	double alpha = norm(alphaB);							//Valor escalar do parametro alpha
	Matrix A = skew(alphaB);								//Matriz A
	double g = 4.0 / (4.0 + alpha * alpha);					//função g(alpha) - em algumas ref. tb. chamado de h(alpha)
	Matrix Qd = *I3 + g * (A + 0.5*(A*A));					//Tensor de rotação
	*Q_ABp = Qd * (*Q_ABi);
	//x_ABp
	*x_ABp = *x_ABi + uB;
}

//Retorna coordenadas globais do ponto central da superfície a ser utilizado para cálculos grosseiros de sua localização (pinball)
void RigidArcRevolution_1::CenterPoint(Matrix* center)
{
	//*center = *x_ABi;

	int arc2id = arc_ID;

	Matrix chord = 0.5*(db.arcs[arc2id - 1]->f_point - db.arcs[arc2id - 1]->i_point);
	Matrix midarc = (db.arcs[arc2id - 1]->i_point + chord);

	Matrix centerarc2(3);
	centerarc2(0, 0) = midarc(0, 0);
	centerarc2(1, 0) = midarc(1, 0);
	centerarc2(2, 0) = 0;

	Matrix centerarc2glob;
	centerarc2glob = *Q_ABi*centerarc2;

	*center = centerarc2glob;
}

//Salva vetores de configuração convergida
void RigidArcRevolution_1::SaveConfiguration()
{
	d_B->clear();
	for (int i = 0; i < 3; i++)
	{
		(*x_ABi)(i, 0) = db.nodes[nodes[0] - 1]->copy_coordinates[i];
		(*alpha_AB)(i, 0) = db.nodes[nodes[0] - 1]->copy_coordinates[i + 3];

		(*dui_B)(i + 0, 0) = db.nodes[nodes[0] - 1]->copy_vel[i];
		(*dui_B)(i + 3, 0) = db.nodes[nodes[0] - 1]->copy_vel[i + 3];

		(*ddui_B)(i + 0, 0) = db.nodes[nodes[0] - 1]->copy_accel[i];
		(*ddui_B)(i + 3, 0) = db.nodes[nodes[0] - 1]->copy_accel[i + 3];
	}
	
	//Q_ABi
	double alpha = norm(*alpha_AB);							//Valor escalar do parametro alpha
	Matrix A = skew(*alpha_AB);								//Matriz A
	double g = 4.0 / (4.0 + alpha*alpha);					//função g(alpha) - em algumas ref. tb. chamado de h(alpha)
	*Q_ABi = *I3 + g*(A + 0.5*(A*A));						//Tensor de rotação
	*Q_ABic = *Q_ABi;										//Cópia - antes da transformação embutida
	*Q_ABi = (*Q_ABi)*(*Q0B);								//Para já realizar a transformação da parametrização do arco para o sistema global
	Q_ABi->MatrixToPtr(aQ_ABi, 3);
	
}

void RigidArcRevolution_1::UpdateBox()
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

	//Criando vetores com os corners do box que englobem o arco j· revolucionado considerando os fatores de ovalizaÁ„o
	Matrix v1(3);
	v1(0, 0) = A(0, 0)*x_fac;
	v1(1, 0) = A(1, 0);
	v1(2, 0) = *radius*z_fac;

	Matrix v2(3);
	v2(0, 0) = D(0, 0)*x_fac;
	v2(1, 0) = D(1, 0);
	v2(2, 0) = *radius*z_fac;

	Matrix v3(3);
	v3(0, 0) = D(0, 0)*x_fac;
	v3(1, 0) = D(1, 0);
	v3(2, 0) = -(*radius)*z_fac;

	Matrix v4(3);
	v4(0, 0) = A(0, 0)*x_fac;
	v4(1, 0) = A(1, 0);
	v4(2, 0) = -(*radius)*z_fac;

	Matrix v5(3);
	v5(0, 0) = -A(0, 0)*x_fac;
	v5(1, 0) = A(1, 0);
	v5(2, 0) = *radius*z_fac;

	Matrix v6(3);
	v6(0, 0) = -D(0, 0)*x_fac;
	v6(1, 0) = D(1, 0);
	v6(2, 0) = *radius*z_fac;

	Matrix v7(3);
	v7(0, 0) = -D(0, 0)*x_fac;
	v7(1, 0) = D(1, 0);
	v7(2, 0) = -(*radius)*z_fac;

	Matrix v8(3);
	v8(0, 0) = -A(0, 0)*x_fac;
	v8(1, 0) = A(1, 0);
	v8(2, 0) = -(*radius)*z_fac;

	//AssociaÁ„o com o nÛ
	Matrix xABi(3);
	xABi = *x_ABi;
	Matrix QABi(3, 3);
	QABi = *Q_ABi;

	double x[8], y[8], z[8];

	v1 = QABi * v1 + xABi;
	x[0] = v1(0, 0);
	y[0] = v1(1, 0);
	z[0] = v1(2, 0);

	v2 = QABi * v2 + xABi;
	x[1] = v2(0, 0);
	y[1] = v2(1, 0);
	z[1] = v2(2, 0);

	v3 = QABi * v3 + xABi;
	x[2] = v3(0, 0);
	y[2] = v3(1, 0);
	z[2] = v3(2, 0);

	v4 = QABi * v4 + xABi;
	x[3] = v4(0, 0);
	y[3] = v4(1, 0);
	z[3] = v4(2, 0);

	v5 = QABi * v5 + xABi;
	x[4] = v5(0, 0);
	y[4] = v5(1, 0);
	z[4] = v5(2, 0);

	v6 = QABi * v6 + xABi;
	x[5] = v6(0, 0);
	y[5] = v6(1, 0);
	z[5] = v6(2, 0);

	v7 = QABi * v7 + xABi;
	x[6] = v7(0, 0);
	y[6] = v7(1, 0);
	z[6] = v7(2, 0);

	v8 = QABi * v8 + xABi;
	x[7] = v8(0, 0);
	y[7] = v8(1, 0);
	z[7] = v8(2, 0);

	//Setando os vÈrtices
	box.SetVertices(x, y, z);
}

//Calcula contribuições de contato entre esfera e superfície
void RigidArcRevolution_1::ContactSphereSurfaceSticking(double* Rc, double** Kc, double zetap, double thetap, double zetai, double thetai, double* gti, int node, double* epsn, double* epst, double* cn, double* ct, double* mu, double* radius)
{

}

//Calcula contribuições de contato entre esfera e superfície
void RigidArcRevolution_1::ContactSphereSurfaceSliding(double* Rc, double** Kc, double zetap, double thetap, double zetai, double thetai, double* gti, int node, double* epsn, double* epst, double* cn, double* ct, double* mu, double* radius)
{

}