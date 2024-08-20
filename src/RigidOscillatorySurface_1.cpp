#include "RigidOscillatorySurface_1.h"

#include"Database.h"
//Variáveis globais
extern
Database db;

RigidOscillatorySurface_1::RigidOscillatorySurface_1()
{
	nDOFs = 6;
	pilot_node = 1;
	pilot_is_used = true;
	vNR = new Matrix(2);
	n_nodes = 1;

	xP_i = new Matrix(3);
	e1_i = new Matrix(3);
	e2_i = new Matrix(3);

	xP_p = new Matrix(3);
	e1_p = new Matrix(3);
	e2_p = new Matrix(3);

	number_CS = 0;
	A_1 = 0.0;
	A_2 = 0.0;
	A_12 = 0.0;
	lambda_1 = 0;
	lambda_2 = 0;
	phi_1 = 0;
	phi_2 = 0;

	waves_1 = 1.0;
	waves_2 = 1.0;

	I3 = new Matrix(3, 3);
	(*I3)(0, 0) = 1.0;
	(*I3)(1, 1) = 1.0;
	(*I3)(2, 2) = 1.0;
	GLs = new int*[nDOFs];
	for (int i = 0; i < nDOFs; i++)
		GLs[i] = NULL;

	InitializeDegeneration();
}

RigidOscillatorySurface_1::~RigidOscillatorySurface_1()
{
	delete vNR;
	delete xP_i;
	delete e1_i;
	delete e2_i;
	delete xP_p;
	delete e1_p;
	delete e2_p;
	delete I3;
	delete[]GLs;

	FreeDegeneration();
}

void RigidOscillatorySurface_1::SetMinMaxRange()
{
	u1_min = -0.5*waves_1*lambda_1;
	u1_max = +0.5*waves_1*lambda_1;
	u1_range = waves_1 * lambda_1;

	u2_min = -0.5*waves_2 * lambda_2;
	u2_max = +0.5*waves_2 * lambda_2;
	u2_range = waves_2 * lambda_2;
}

//Obtem ponto da superficie
void RigidOscillatorySurface_1::SurfacePoint(double& zeta, double& theta, Matrix& point)
{
	//TODO
}

//Normal exterior à superfície na posição escolhida
void RigidOscillatorySurface_1::NormalExt(double* zeta, double* theta, Matrix* n)
{

}

//Atualiza bounding box
void RigidOscillatorySurface_1::UpdateBox()
{

}

void RigidOscillatorySurface_1::WriteVTK_XMLRender(FILE *f)
{
	if (db.post_files->WriteRigidContactSurfaces_flag == true)
	{
		//vetores para escrita no formato binário - usando a função 'enconde'
		std::vector<float> float_vector;
		std::vector<int> int_vector;

		int n_div_wave_1 = 12;
		int n_div_wave_2 = 12;
		int divisions1 = (int)(waves_1*n_div_wave_1);
		int divisions2 = (int)(waves_2*n_div_wave_2);

		int n_points = (2 * divisions1 + 1)*(2 * divisions2 + 1);
		int n_cells = divisions1*divisions2;	//triângulo + pilot
		Matrix pilot(3);
		Matrix QM(3, 3);
		Matrix e3_p = cross(*e1_p, *e2_p);
		Matrix local(3);
		QM(0, 0) = (*e1_p)(0, 0);
		QM(0, 1) = (*e2_p)(0, 0);
		QM(0, 2) = (e3_p)(0, 0);

		QM(1, 0) = (*e1_p)(1, 0);
		QM(1, 1) = (*e2_p)(1, 0);
		QM(1, 2) = (e3_p)(1, 0);

		QM(2, 0) = (*e1_p)(2, 0);
		QM(2, 1) = (*e2_p)(2, 0);
		QM(2, 2) = (e3_p)(2, 0);
		double zeta, theta;
		//Opens Piece
		fprintf(f, "\t\t<Piece NumberOfPoints = \"%d\" NumberOfCells = \"%d\">\n", n_points, n_cells);
		//Opens Points
		fprintf(f, "\t\t\t<Points>\n");
		//Opens DataArray
		fprintf(f, "\t\t\t\t<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"binary\">\n");
		float_vector.clear();
		//Pilot node
		//Posição de cada ponto P no plano xy (referência)
		pilot(0, 0) = db.post_files->mag_factor*(db.nodes[pilot_node - 1]->copy_coordinates[0] - db.nodes[pilot_node - 1]->ref_coordinates[0]) + db.nodes[pilot_node - 1]->ref_coordinates[0];
		pilot(1, 0) = db.post_files->mag_factor*(db.nodes[pilot_node - 1]->copy_coordinates[1] - db.nodes[pilot_node - 1]->ref_coordinates[1]) + db.nodes[pilot_node - 1]->ref_coordinates[1];
		pilot(2, 0) = db.post_files->mag_factor*(db.nodes[pilot_node - 1]->copy_coordinates[2] - db.nodes[pilot_node - 1]->ref_coordinates[2]) + db.nodes[pilot_node - 1]->ref_coordinates[2];
		for (int i = 0; i < (2 * divisions1 + 1); i++)
		{
			for (int j = 0; j < (2 * divisions2 + 1); j++)
			{
				zeta = -0.5*waves_1*lambda_1 + waves_1*lambda_1*i / (2 * divisions1);
				theta = -0.5*waves_2*lambda_2 + waves_2*lambda_2*j / (2 * divisions2);
				local(0, 0) = zeta;
				local(1, 0) = theta;
				local(2, 0) = A_1*sin(2 * PI * zeta / lambda_1 + phi_1) + A_2*sin(2 * PI * theta / lambda_2 + phi_2) + A_12*sin(2 * PI * zeta / lambda_1 + phi_1)*sin(2 * PI * theta / lambda_2 + phi_2);
				local = QM*local;
				float_vector.push_back((float)(pilot(0, 0) + local(0, 0)));
				float_vector.push_back((float)(pilot(1, 0) + local(1, 0)));
				float_vector.push_back((float)(pilot(2, 0) + local(2, 0)));
				//fprintf(f, "\t\t\t\t\t%f\t%f\t%f\n", pilot(0, 0) + local(0, 0), pilot(1, 0) + local(1, 0), pilot(2, 0) + local(2, 0));
			}
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
		int nj = (2 * divisions2 + 1);
		int_vector.clear();
		for (int i = 0; i < (2 * divisions1 - 1); i = i + 2)
		for (int j = 0; j < (2 * divisions2 - 1); j = j + 2)
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
		//Closes Piece
		fprintf(f, "\t\t</Piece>\n");
	}
}

bool RigidOscillatorySurface_1::Read(FILE *f)
{
	char s[1000];
	fscanf(f, "%s", s);
	number = atoi(s);

	fscanf(f, "%s", s);
	if (!strcmp(s, "A1"))
	{
		fscanf(f, "%s", s);
		A_1 = atof(s);
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "A2"))
	{
		fscanf(f, "%s", s);
		A_2 = atof(s);
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "A12"))
	{
		fscanf(f, "%s", s);
		A_12 = atof(s);
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "Lambda1"))
	{
		fscanf(f, "%s", s);
		lambda_1 = atof(s);
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "Lambda2"))
	{
		fscanf(f, "%s", s);
		lambda_2 = atof(s);
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "Phi1"))
	{
		fscanf(f, "%s", s);
		phi_1 = atof(s);
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "Phi2"))
	{
		fscanf(f, "%s", s);
		phi_2 = atof(s);
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "Waves1"))
	{
		fscanf(f, "%s", s);
		waves_1 = atof(s);
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "Waves2"))
	{
		fscanf(f, "%s", s);
		waves_2 = atof(s);
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "CS"))
	{
		fscanf(f, "%s", s);
		number_CS = atoi(s);
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "PilotNode"))
	{
		fscanf(f, "%s", s);
		pilot_node = atoi(s);
	}
	else
		return false;

	ReadCommon(f);

	return true;
}

void RigidOscillatorySurface_1::Write(FILE *f)
{
	fprintf(f, "RigidOscillatorySurface_1\t%d\tA1\t%.6e\tA2\t%.6e\tA12\t%.6e\tLambda1\t%.6e\tLambda2\t%.6e\tPhi1\t%.6e\tPhi2\t%.6e\tWaves1\t%.6e\tWaves2\t%.6e\tCS\t%d\tPilotNode\t%d\n", number,A_1,A_2,A_12,lambda_1,lambda_2,phi_1,phi_2,waves_1,waves_2,number_CS,pilot_node);
}

//Checa inconsistências para evitar erros de execução
bool RigidOscillatorySurface_1::Check()
{
	if (pilot_node > db.number_nodes)
		return false;
	if (number_CS > db.number_CS)
		return false;
	return true;
}

void RigidOscillatorySurface_1::PreCalc()
{
	//Setando valores dos versores locais - de acordo com a orientação do sistema de coordenadas local atribuído à superfície
	*e1_i = *(db.CS[number_CS - 1]->E1);
	*e2_i = *(db.CS[number_CS - 1]->E2);
	(*xP_i)(0, 0) = db.nodes[pilot_node - 1]->ref_coordinates[0];
	(*xP_i)(1, 0) = db.nodes[pilot_node - 1]->ref_coordinates[1];
	(*xP_i)(2, 0) = db.nodes[pilot_node - 1]->ref_coordinates[2];
	*xP_p = *xP_i;
	*e1_p = *e1_i;
	*e2_p = *e2_i;

	//Apontando para posição que indica valor dos GLs globais
	for (int i = 0; i < 6; i++)
		GLs[i] = &db.nodes[pilot_node - 1]->GLs[i];

	DegenerationPreCalc();
}

//Retorna as coordenadas da superfície para um par (zeta,theta) - configuração anterior convergida
void RigidOscillatorySurface_1::Gamma_and_Triad(Matrix* G_p, Matrix* t1_p, Matrix* t2_p, Matrix* n_p, Matrix* G_i, Matrix* t1_i, Matrix* t2_i, Matrix* n_i, Matrix* G_ip, double* zi, double* thi, double* zp, double* thp)
{
	double *Gp = G_p->getMatrix();		//ponteiro para o vetor Gp 
	double *t1p = t1_p->getMatrix();	//ponteiro para o vetor t1p 
	double *t2p = t2_p->getMatrix();	//ponteiro para o vetor t2p
	double *np = n_p->getMatrix();		//ponteiro para o vetor np

	double *t1i = t1_i->getMatrix();	//ponteiro para o vetor t1i 
	double *t2i = t2_i->getMatrix();	//ponteiro para o vetor t2i
	double *ni = n_i->getMatrix();		//ponteiro para o vetor ni

	double *Gip = G_ip->getMatrix();	//ponteiro para o vetor Gip
	double *Gi = G_i->getMatrix();		//ponteiro para o vetor Gi

	double* e1i = e1_i->getMatrix();
	double* e2i = e2_i->getMatrix();
	double* e1p = e1_p->getMatrix();
	double* e2p = e2_p->getMatrix();

	double* xPi = xP_i->getMatrix();
	double* xPp = xP_p->getMatrix();

	double* A1 = &A_1;
	double* A2 = &A_2;
	double* A12 = &A_12;
	double* lambda1 = &lambda_1;
	double* lambda2 = &lambda_2;
	double* phi1 = &phi_1;
	double* phi2 = &phi_2;
	double v[1000];

	v[151] = -(e1p[1] * e2p[0]) + e1p[0] * e2p[1];
	v[150] = e1p[2] * e2p[0] - e1p[0] * e2p[2];
	v[149] = -(e1p[2] * e2p[1]) + e1p[1] * e2p[2];
	v[148] = -(e1i[1] * e2i[0]) + e1i[0] * e2i[1];
	v[147] = e1i[2] * e2i[0] - e1i[0] * e2i[2];
	v[146] = -(e1i[2] * e2i[1]) + e1i[1] * e2i[2];
	v[145] = 1e0 / (*lambda2);
	v[144] = 0.6283185307179586e1;
	v[143] = v[144] / (*lambda1);
	v[71] = v[144] * v[145];
	v[107] = (*phi1) + v[143] * (*zp);
	v[108] = v[143] * cos(v[107]);
	v[80] = (*phi1) + v[143] * (*zi);
	v[81] = v[143] * cos(v[80]);
	v[117] = (*phi2) + (*thp)*v[71];
	v[118] = v[71] * cos(v[117]);
	v[90] = (*phi2) + (*thi)*v[71];
	v[91] = v[71] * cos(v[90]);
	v[63] = sin(v[80]);
	v[152] = (*A2) + (*A12)*v[63];
	v[92] = v[152] * v[91];
	v[64] = sin(v[90]);
	v[82] = ((*A1) + (*A12)*v[64])*v[81];
	v[66] = (*A1)*v[63] + v[152] * v[64];
	v[70] = sin(v[107]);
	v[153] = (*A2) + (*A12)*v[70];
	v[119] = v[118] * v[153];
	v[72] = sin(v[117]);
	v[109] = v[108] * ((*A1) + (*A12)*v[72]);
	v[74] = (*A1)*v[70] + v[153] * v[72];
	v[83] = e1i[0] + v[146] * v[82];
	v[84] = e1i[1] + v[147] * v[82];
	v[85] = e1i[2] + v[148] * v[82];
	v[87] = 1e0 / sqrt(Power(v[83], 2) + Power(v[84], 2) + Power(v[85], 2));
	v[86] = v[83] * v[87];
	v[88] = v[84] * v[87];
	v[89] = v[85] * v[87];
	v[93] = e2i[0] + v[146] * v[92];
	v[94] = e2i[1] + v[147] * v[92];
	v[95] = e2i[2] + v[148] * v[92];
	v[97] = 1e0 / sqrt(Power(v[93], 2) + Power(v[94], 2) + Power(v[95], 2));
	v[96] = v[93] * v[97];
	v[98] = v[94] * v[97];
	v[105] = -(v[88] * v[96]) + v[86] * v[98];
	v[99] = v[95] * v[97];
	v[103] = v[89] * v[96] - v[86] * v[99];
	v[100] = -(v[89] * v[98]) + v[88] * v[99];
	v[102] = 1e0 / sqrt(Power(v[100], 2) + Power(v[103], 2) + Power(v[105], 2));
	v[110] = e1p[0] + v[109] * v[149];
	v[111] = e1p[1] + v[109] * v[150];
	v[112] = e1p[2] + v[109] * v[151];
	v[114] = 1e0 / sqrt(Power(v[110], 2) + Power(v[111], 2) + Power(v[112], 2));
	v[113] = v[110] * v[114];
	v[115] = v[111] * v[114];
	v[116] = v[112] * v[114];
	v[120] = e2p[0] + v[119] * v[149];
	v[121] = e2p[1] + v[119] * v[150];
	v[122] = e2p[2] + v[119] * v[151];
	v[124] = 1e0 / sqrt(Power(v[120], 2) + Power(v[121], 2) + Power(v[122], 2));
	v[123] = v[120] * v[124];
	v[125] = v[121] * v[124];
	v[132] = -(v[115] * v[123]) + v[113] * v[125];
	v[126] = v[122] * v[124];
	v[130] = v[116] * v[123] - v[113] * v[126];
	v[127] = -(v[116] * v[125]) + v[115] * v[126];
	v[129] = 1e0 / sqrt(Power(v[127], 2) + Power(v[130], 2) + Power(v[132], 2));
	Gip[0] = e2p[0] * (*thi) + v[149] * v[66] + xPp[0] + e1p[0] * (*zi);
	Gip[1] = e2p[1] * (*thi) + v[150] * v[66] + xPp[1] + e1p[1] * (*zi);
	Gip[2] = e2p[2] * (*thi) + v[151] * v[66] + xPp[2] + e1p[2] * (*zi);
	t1i[0] = v[86];
	t1i[1] = v[88];
	t1i[2] = v[89];
	t2i[0] = v[96];
	t2i[1] = v[98];
	t2i[2] = v[99];
	ni[0] = v[100] * v[102];
	ni[1] = v[102] * v[103];
	ni[2] = v[102] * v[105];
	Gp[0] = e2p[0] * (*thp) + v[149] * v[74] + xPp[0] + e1p[0] * (*zp);
	Gp[1] = e2p[1] * (*thp) + v[150] * v[74] + xPp[1] + e1p[1] * (*zp);
	Gp[2] = e2p[2] * (*thp) + v[151] * v[74] + xPp[2] + e1p[2] * (*zp);
	t1p[0] = v[113];
	t1p[1] = v[115];
	t1p[2] = v[116];
	t2p[0] = v[123];
	t2p[1] = v[125];
	t2p[2] = v[126];
	np[0] = v[127] * v[129];
	np[1] = v[129] * v[130];
	np[2] = v[129] * v[132];
	Gi[0] = e2i[0] * (*thi) + v[146] * v[66] + xPi[0] + e1i[0] * (*zi);
	Gi[1] = e2i[1] * (*thi) + v[147] * v[66] + xPi[1] + e1i[1] * (*zi);
	Gi[2] = e2i[2] * (*thi) + v[148] * v[66] + xPi[2] + e1i[2] * (*zi);
}

//Realiza chute inicial para as variáveis zeta e theta
void RigidOscillatorySurface_1::InitialGuess(Matrix* xS, double** convective, int n_solutions)
{
	Matrix QMt(3, 3);
	Matrix e3_p = cross(*e1_p, *e2_p);
	QMt(0, 0) = (*e1_p)(0, 0);
	QMt(1, 0) = (*e2_p)(0, 0);
	QMt(2, 0) = (e3_p)(0, 0);

	QMt(0, 1) = (*e1_p)(1, 0);
	QMt(1, 1) = (*e2_p)(1, 0);
	QMt(2, 1) = (e3_p)(1, 0);

	QMt(0, 2) = (*e1_p)(2, 0);
	QMt(1, 2) = (*e2_p)(2, 0);
	QMt(2, 2) = (e3_p)(2, 0);
	Matrix temp = QMt*(*xS - *xP_p);
	convective[0][0] = temp(0, 0);
	convective[0][1] = temp(1, 0);
}

//Retorna em Gamma o valor de Gamma para zeta e theta (tomando posição atual)
void RigidOscillatorySurface_1::Gamma(Matrix* Gamma, double *zeta, double *theta)
{
	double* A1 = &A_1;
	double* A2 = &A_2;
	double* A12 = &A_12;
	double* lambda1 = &lambda_1;
	double* lambda2 = &lambda_2;
	double* phi1 = &phi_1;
	double* phi2 = &phi_2;
	double* e1 = e1_p->getMatrix();
	double* e2 = e2_p->getMatrix();
	double* G = Gamma->getMatrix();
	double* xP = xP_p->getMatrix();
	double v[100];

	v[34] = 1e0 / (*lambda2);
	v[33] = 0.6283185307179586e1;
	v[32] = sin((*phi1) + (v[33] * (*zeta)) / (*lambda1));
	v[26] = sin((*phi2) + (*theta)*v[33] * v[34]);
	v[28] = (*A1)*v[32] + v[26] * ((*A2) + (*A12)*v[32]);
	G[0] = e2[0] * (*theta) + (-(e1[2] * e2[1]) + e1[1] * e2[2])*v[28] + xP[0] + e1[0] * (*zeta);
	G[1] = e2[1] * (*theta) + (e1[2] * e2[0] - e1[0] * e2[2])*v[28] + xP[1] + e1[1] * (*zeta);
	G[2] = e2[2] * (*theta) + (-(e1[1] * e2[0]) + e1[0] * e2[1])*v[28] + xP[2] + e1[2] * (*zeta);
}

//Dado o ponto xS, calcula as coordenadas (zeta,theta) referentes à mínima distância
void RigidOscillatorySurface_1::FindMinimimumParameters(Matrix* xS, NSContactData* cd)
{
	//Parametros NR
	double tol_ortho = 1e-12;
	double max_delta = 0.0;
	int max_it = 20;
	double v[1000];
	double** J;
	//Alocação J
	J = new double*[2];
	for (int i = 0; i < 2; i++)
		J[i] = new double[2];
	Matrix Jacobian(2, 2);
	Matrix residual(2);
	Matrix delta(2);
	double *r = residual.getMatrix();		//ponteiro para o vetor residual
	double *vi = vNR->getMatrix();			//ponteiro para vNR

	double* e1 = e1_p->getMatrix();
	double* e2 = e2_p->getMatrix();
	double* xP = xP_p->getMatrix();
	double* A1 = &A_1;
	double* A2 = &A_2;
	double* A12 = &A_12;
	double* lambda1 = &lambda_1;
	double* lambda2 = &lambda_2;
	double* phi1 = &phi_1;
	double* phi2 = &phi_2;
	double *pxS = xS->getMatrix();			//ponteiro para o vetor xS
	
	for (int sol_index = 0; sol_index < cd->n_solutions; sol_index++)
	{
		//Inicialização de chute inicial
		(*vNR)(0, 0) = cd->convective[sol_index][0];
		(*vNR)(1, 0) = cd->convective[sol_index][1];
		
		if (sol_index == 1 && cd->repeated[sol_index] == true)//tentativa de encontrar outras raízes
		{
			(*vNR)(0, 0) = cd->convective[0][0] + lambda_1;
			(*vNR)(1, 0) = cd->convective[0][1];
		}
		if (sol_index == 2 && cd->repeated[sol_index] == true)//tentativa de encontrar outras raízes
		{
			(*vNR)(0, 0) = cd->convective[0][0] - lambda_1;
			(*vNR)(1, 0) = cd->convective[0][1];
		}
		if (sol_index == 3 && cd->repeated[sol_index] == true)//tentativa de encontrar outras raízes
		{
			(*vNR)(0, 0) = cd->convective[0][0];
			(*vNR)(1, 0) = cd->convective[0][1] + lambda_2;
		}
		if (sol_index == 4 && cd->repeated[sol_index] == true)//tentativa de encontrar outras raízes
		{
			(*vNR)(0, 0) = cd->convective[0][0];
			(*vNR)(1, 0) = cd->convective[0][1] - lambda_2;
		}
		if (sol_index > 4 && cd->repeated[sol_index] == true)//tentativa de encontrar outras raízes
		{
			(*vNR)(0, 0) = cd->convective[0][0];
			(*vNR)(1, 0) = cd->convective[0][1];
		}

		double error = tol_ortho + 1.0;	//Forçando a entrar no loop
		int it = 1;
		int flag_error = 0;
		while (error > tol_ortho && it <= max_it)
		{
			if (error > tol_ortho)
			{
				//Cálculo do resíduo e jacobiano
				v[96] = -(e2[2] * vi[1]);
				v[95] = -(e1[2] * vi[0]);
				v[94] = -(e2[1] * vi[1]);
				v[93] = -(e1[1] * vi[0]);
				v[92] = -(e2[0] * vi[1]);
				v[91] = -(e1[0] * vi[0]);
				v[88] = -(e1[1] * e2[0]) + e1[0] * e2[1];
				v[87] = e1[2] * e2[0] - e1[0] * e2[2];
				v[86] = -(e1[2] * e2[1]) + e1[1] * e2[2];
				v[85] = 1e0 / (*lambda2);
				v[84] = 0.6283185307179586e1;
				v[83] = v[84] / (*lambda1);
				v[45] = v[84] * v[85];
				v[38] = (*phi1) + v[83] * vi[0];
				v[56] = sin(v[38]);
				v[89] = (*A2) + (*A12)*v[56];
				v[57] = -(v[56] * (v[83] * v[83]));
				v[39] = v[83] * cos(v[38]);
				v[46] = (*phi2) + v[45] * vi[1];
				v[58] = sin(v[46]);
				v[90] = (*A1) + (*A12)*v[58];
				v[59] = -((v[45] * v[45])*v[58]);
				v[47] = v[45] * cos(v[46]);
				v[60] = (*A12)*v[39] * v[47];
				v[61] = v[59] * v[89];
				v[48] = v[47] * v[89];
				v[68] = e2[2] + v[48] * v[88];
				v[66] = e2[1] + v[48] * v[87];
				v[64] = e2[0] + v[48] * v[86];
				v[62] = v[57] * v[90];
				v[40] = v[39] * v[90];
				v[67] = e1[2] + v[40] * v[88];
				v[65] = e1[1] + v[40] * v[87];
				v[63] = e1[0] + v[40] * v[86];
				v[34] = (*A1)*v[56] + v[58] * v[89];
				v[52] = pxS[0] - v[34] * v[86] + v[91] + v[92] - xP[0];
				v[53] = pxS[1] - v[34] * v[87] + v[93] + v[94] - xP[1];
				v[54] = pxS[2] - v[34] * v[88] + v[95] + v[96] - xP[2];
				v[97] = v[52] * v[86] + v[53] * v[87] + v[54] * v[88];
				v[79] = -(v[63] * v[64]) - v[65] * v[66] - v[67] * v[68] + v[60] * v[97];
				r[0] = v[52] * v[63] + v[53] * v[65] + v[54] * v[67];
				r[1] = v[52] * v[64] + v[53] * v[66] + v[54] * v[68];
				J[0][0] = -(v[63] * v[63]) - (v[65] * v[65]) - (v[67] * v[67]) + v[62] * v[97];
				J[0][1] = v[79];
				J[1][0] = v[79];
				J[1][1] = -(v[64] * v[64]) - (v[66] * v[66]) - (v[68] * v[68]) + v[61] * v[97];
				Jacobian.PtrToMatrix(J, 2);									//Conversão do Jacobiano para matriz
				delta = fullsystem(Jacobian, -1.0*residual, &flag_error);	//Resolve sistema linear
				if (flag_error == 0)										//Se conseguiu fazer o sistema linear
				{
					(*vNR) = (*vNR) + delta;								//Atualização das variáveis
					error = norm(delta);									//Norma do residuo
					it++;
				}
				else
				{
					it = max_it + 1;//Força saída - divergência
				}
			}
		}//end while NR
		
		

		//Convergiu - ainda há ações a verificar...
		if (error <= tol_ortho && flag_error == 0)
		{
			if (norm(delta) > max_delta)
				max_delta = norm(delta);
			cd->convective[sol_index][0] = (*vNR)(0, 0);
			cd->convective[sol_index][1] = (*vNR)(1, 0);
			//Se está no range local de interesse - domínio físico da superfície triangular
			if (abs((*vNR)(0, 0)) <= 0.5*waves_1*lambda_1 && abs((*vNR)(1, 0)) <= 0.5*waves_2*lambda_2)
				cd->return_value[sol_index] = 0;
			else
			{
				//Se está em região próxima, mas não no range local de interesse
				if (abs((*vNR)(0, 0)) <= 0.6*waves_1*lambda_1 && abs((*vNR)(1, 0)) <= 0.6*waves_2*lambda_2)
					cd->return_value[sol_index] = 3;
				//Se não estiver no range de interesse
				else
					cd->return_value[sol_index] = 2;
			}
			
		}
		else
			cd->return_value[sol_index] = 1;
		//Retornos da função
		//0 - Convergiu e está no range de interesse para contato
		//1 - Não houve convergência (pode ou não estar no range para contato) - retorno problemático!!
		//2 - Houve convergência, mas não está no range para contato
		//3 - Houve convergência, está fora do range para contato, mas próximo
	}
	//Clean
	for (int i = 0; i < 2; i++)
		delete[]J[i];
	delete[]J;

	//cd->CheckRepeated(1e-6);
	//cd->CheckRepeated(1e-3);
	cd->CheckRepeated(max_delta*1e3);
	//cd->Plot();
	//printf("%d %d %d %d %d\n", cd->repeated[0], cd->repeated[1], cd->repeated[2], cd->repeated[3], cd->repeated[4]);
}
//Atualiza as variáveis internas da superfície, para pegarem info do pilot node para uso posterior com posição atualizada
void RigidOscillatorySurface_1::FillNodes()
{
	//Atualização do Pilot Node
	(*xP_p)(0, 0) = (*xP_i)(0, 0) + db.nodes[pilot_node - 1]->displacements[0];
	(*xP_p)(1, 0) = (*xP_i)(1, 0) + db.nodes[pilot_node - 1]->displacements[1];
	(*xP_p)(2, 0) = (*xP_i)(2, 0) + db.nodes[pilot_node - 1]->displacements[2];

	//Atualização dos versores orientativos da superfície
	Matrix alpha_1(3);
	alpha_1(0, 0) = db.nodes[pilot_node - 1]->displacements[3];
	alpha_1(1, 0) = db.nodes[pilot_node - 1]->displacements[4];
	alpha_1(2, 0) = db.nodes[pilot_node - 1]->displacements[5];
	double alpha_escalar = norm(alpha_1);
	Matrix A = skew(alpha_1);
	double g = 4.0 / (4.0 + alpha_escalar*alpha_escalar);
	Matrix Q = *I3 + g*(A + 0.5*((A)*(A)));
	*e1_p = Q*(*e1_i);
	*e2_p = Q*(*e2_i);
}

//Retorna coordenadas globais do ponto central da superfície a ser utilizado para cálculos grosseiros de sua localização (pinball)
void RigidOscillatorySurface_1::CenterPoint(Matrix* center)
{
	*center = *xP_i;
}

//Salva vetores de configuração convergida
void RigidOscillatorySurface_1::SaveConfiguration()
{
	*xP_i = *xP_p;
	*e1_i = *e1_p;
	*e2_i = *e2_p;
}

//Calcula contribuições de contato entre esfera e superfície
void RigidOscillatorySurface_1::ContactSphereSurfaceSticking(double* Rc, double** Kc, double zetap, double thetap, double zetai, double thetai, double* gti, int node, double* epsn, double* epst, double* cn, double* ct, double* mu, double* radius)
{
	double *xPi = xP_i->getMatrix();	//ponteiro para o vetor xPi
	double* e1i = e1_i->getMatrix();
	double* e2i = e2_i->getMatrix();
	double* A1 = &A_1;
	double* A2 = &A_2;
	double* A12 = &A_12;
	double* lambda1 = &lambda_1;
	double* lambda2 = &lambda_2;
	double* phi1 = &phi_1;
	double* phi2 = &phi_2;

	double ci[2];
	ci[0] = zetai;
	ci[1] = thetai;
	double cp[2];
	cp[0] = zetap;
	cp[1] = thetap;
	double xSi[3];
	double d[12];
	for (int i = 0; i < 3; i++)
	{
		/*(*1 - 6: sphere*)
		(*7 - 12: surface*)*/
		xSi[i] = db.nodes[node - 1]->copy_coordinates[i];
		d[i] = db.nodes[node - 1]->displacements[i];
		d[i + 3] = db.nodes[node - 1]->displacements[i + 3];
		d[i + 6] = db.nodes[pilot_node - 1]->displacements[i];
		d[i + 9] = db.nodes[pilot_node - 1]->displacements[i + 3];
	}
	
	double* a4;
	double* a5;
	double* a6;
	double value = 0.0;
	if (typeid(*db.solution[db.current_solution_number - 1]) == typeid(Dynamic))
	{
		Dynamic* ptr_sol = static_cast<Dynamic*>(db.solution[db.current_solution_number - 1]);
		a4 = &ptr_sol->a4;
		a5 = &ptr_sol->a5;
		a6 = &ptr_sol->a6;
	}
	else
	{
		a4 = &value;
		a5 = &value;
		a6 = &value;
	}
	double omegaiS[3];
	double domegaiS[3];
	double duiS[3];
	double dduiS[3];
	double omegaiP[3];
	double domegaiP[3];
	double duiP[3];
	double dduiP[3];
	for (int i = 0; i < 3; i++)
	{
		duiS[i] = db.nodes[node - 1]->copy_vel[i];
		omegaiS[i] = db.nodes[node - 1]->copy_vel[i + 3];
		dduiS[i] = db.nodes[node - 1]->copy_accel[i];
		domegaiS[i] = db.nodes[node - 1]->copy_accel[i + 3];

		duiP[i] = db.nodes[pilot_node - 1]->copy_vel[i];
		omegaiP[i] = db.nodes[pilot_node - 1]->copy_vel[i + 3];
		dduiP[i] = db.nodes[pilot_node - 1]->copy_accel[i];
		domegaiP[i] = db.nodes[pilot_node - 1]->copy_accel[i + 3];
	}
	double v[10000];
	//AceGen
	v[4794] = (*a4)*d[2] + (*a6)*dduiS[2] + (*a5)*duiS[2];
	v[4793] = (*a4)*d[1] + (*a6)*dduiS[1] + (*a5)*duiS[1];
	v[4792] = (*a4)*d[0] + (*a6)*dduiS[0] + (*a5)*duiS[0];
	v[4787] = (*a4)*d[8] + (*a6)*dduiP[2] + (*a5)*duiP[2];
	v[4797] = -v[4787] + v[4794];
	v[4786] = (*a4)*d[7] + (*a6)*dduiP[1] + (*a5)*duiP[1];
	v[4796] = -v[4786] + v[4793];
	v[4785] = (*a4)*d[6] + (*a6)*dduiP[0] + (*a5)*duiP[0];
	v[4795] = -v[4785] + v[4792];
	v[4536] = -(e1i[1] * e2i[0]) + e1i[0] * e2i[1];
	v[4535] = e1i[2] * e2i[0] - e1i[0] * e2i[2];
	v[4534] = -(e1i[2] * e2i[1]) + e1i[1] * e2i[2];
	v[4520] = d[2] + xSi[2];
	v[4519] = d[1] + xSi[1];
	v[4518] = d[0] + xSi[0];
	v[4517] = d[8] + xPi[2];
	v[4516] = d[7] + xPi[1];
	v[4515] = d[6] + xPi[0];
	v[4511] = 1e0 / (*lambda2);
	v[4510] = 0.6283185307179586e1;
	v[4509] = v[4510] / (*lambda1);
	v[4508] = (*a4)*d[11] + (*a6)*domegaiP[2] + (*a5)*omegaiP[2];
	v[4507] = (*a4)*d[10] + (*a6)*domegaiP[1] + (*a5)*omegaiP[1];
	v[4506] = (*a4)*d[9] + (*a6)*domegaiP[0] + (*a5)*omegaiP[0];
	v[4505] = (*a4)*d[5] + (*a6)*domegaiS[2] + (*a5)*omegaiS[2];
	v[4504] = (*a4)*d[4] + (*a6)*domegaiS[1] + (*a5)*omegaiS[1];
	v[4503] = (*a4)*d[3] + (*a6)*domegaiS[0] + (*a5)*omegaiS[0];
	v[4490] = Power(d[11], 2);
	v[4489] = 0.5e0*d[11];
	v[4488] = 2e0*d[11];
	v[4487] = Power(d[10], 2);
	v[4486] = 0.5e0*d[10];
	v[4485] = 2e0*d[10];
	v[4484] = Power(d[9], 2);
	v[4483] = 2e0*d[9];
	v[4482] = 0.5e0*d[9];
	v[4475] = Power(d[5], 2);
	v[4474] = 0.5e0*d[5];
	v[4473] = 2e0*d[5];
	v[4472] = Power(d[4], 2);
	v[4471] = 0.5e0*d[4];
	v[4470] = 2e0*d[4];
	v[4469] = Power(d[3], 2);
	v[4468] = 2e0*d[3];
	v[4467] = 0.5e0*d[3];
	v[268] = d[4] * v[4467];
	v[1274] = -v[4469] - v[4472];
	v[4545] = 0.5e0*v[1274];
	v[1258] = d[5] + v[268];
	v[1239] = -d[5] + v[268];
	v[275] = d[5] * v[4471];
	v[1276] = d[3] + v[275];
	v[1256] = -d[3] + v[275];
	v[273] = d[5] * v[4467];
	v[1278] = -d[4] + v[273];
	v[1236] = d[4] + v[273];
	v[1254] = -v[4469] - v[4475];
	v[4543] = 0.5e0*v[1254];
	v[1233] = -v[4472] - v[4475];
	v[4542] = 0.5e0*v[1233];
	v[1192] = 4e0 + v[4469] + v[4472] + v[4475];
	v[2735] = 1e0 / Power(v[1192], 3);
	v[4476] = -2e0*v[2735];
	v[2736] = v[4470] * v[4476];
	v[4478] = -4e0*v[2736];
	v[2739] = v[4473] * v[4478];
	v[2734] = v[4468] * v[4476];
	v[4477] = -4e0*v[2734];
	v[2741] = v[4470] * v[4477];
	v[2738] = v[4473] * v[4477];
	v[1194] = 1e0 / Power(v[1192], 2);
	v[4480] = -4e0*v[1194];
	v[2742] = -8e0*v[1194];
	v[2744] = v[2742] + v[4468] * v[4477];
	v[2743] = v[2742] + v[4470] * v[4478];
	v[2740] = v[2742] + 8e0*v[2735] * (v[4473] * v[4473]);
	v[1196] = v[4473] * v[4480];
	v[4479] = -0.5e0*v[1196];
	v[2843] = 2e0*v[1196];
	v[2817] = v[4473] * v[4479];
	v[2780] = -(v[1196] * v[4474]);
	v[2760] = v[1196] * v[4471];
	v[2757] = -(v[1196] * v[4467]);
	v[2753] = v[4468] * v[4479];
	v[2751] = v[4470] * v[4479];
	v[1195] = v[4470] * v[4480];
	v[4544] = v[1195] + v[2757];
	v[4481] = -0.5e0*v[1195];
	v[2853] = -2e0*v[1195];
	v[2823] = v[4468] * v[4481];
	v[2820] = v[1195] * v[4471];
	v[2799] = v[4470] * v[4481];
	v[2756] = -(v[1195] * v[4467]);
	v[4547] = v[1196] + v[2756];
	v[1261] = -(v[1254] * v[4481]);
	v[1193] = v[4468] * v[4480];
	v[4546] = v[1193] - v[2760];
	v[4527] = -0.5e0*v[1193];
	v[3132] = -(v[1193] * v[4467]);
	v[2831] = 2e0*v[1193];
	v[2804] = v[4468] * v[4527];
	v[249] = d[10] * v[4482];
	v[468] = -v[4484] - v[4487];
	v[4492] = 0.5e0*v[468];
	v[482] = -d[11] + v[249];
	v[478] = d[11] + v[249];
	v[256] = d[11] * v[4486];
	v[473] = -d[9] + v[256];
	v[471] = d[9] + v[256];
	v[254] = d[11] * v[4482];
	v[479] = d[10] + v[254];
	v[472] = -d[10] + v[254];
	v[489] = 4e0 + v[4484] + v[4487] + v[4490];
	v[829] = 1e0 / Power(v[489], 3);
	v[4491] = -2e0*v[829];
	v[831] = v[4488] * v[4491];
	v[2427] = -4e0*v[831];
	v[830] = v[4485] * v[4491];
	v[2426] = -4e0*v[830];
	v[833] = v[2426] * v[4488];
	v[828] = v[4483] * v[4491];
	v[2425] = -4e0*v[828];
	v[838] = v[2425] * v[4485];
	v[832] = v[2425] * v[4488];
	v[562] = 1e0 / Power(v[489], 2);
	v[4493] = -4e0*v[562];
	v[839] = -8e0*v[562];
	v[841] = v[2425] * v[4483] + v[839];
	v[840] = v[2426] * v[4485] + v[839];
	v[834] = v[2427] * v[4488] + v[839];
	v[837] = v[4492] * v[834];
	v[565] = v[4488] * v[4493];
	v[3108] = -(v[4489] * v[565]);
	v[2352] = -0.5e0*v[565];
	v[907] = 2e0*v[565];
	v[908] = v[482] * v[834] - v[907];
	v[878] = v[478] * v[834] + v[907];
	v[872] = v[2352] * v[4488];
	v[849] = v[2352] * v[4485];
	v[844] = v[2352] * v[4483];
	v[836] = v[4492] * v[833] + v[849];
	v[835] = v[4492] * v[832] + v[844];
	v[598] = -(v[2352] * v[468]);
	v[564] = v[4485] * v[4493];
	v[2351] = -0.5e0*v[564];
	v[900] = v[4486] * v[564];
	v[880] = -2e0*v[564];
	v[881] = v[479] * v[840] - v[880];
	v[875] = v[4489] * v[564];
	v[926] = v[473] * v[834] - v[849];
	v[923] = v[471] * v[834] - v[849];
	v[899] = v[473] * v[840] - v[849];
	v[885] = v[471] * v[840] - v[849];
	v[862] = v[472] * v[840] + v[880];
	v[853] = v[2351] * v[4485];
	v[563] = v[4483] * v[4493];
	v[4498] = -v[563] + v[875];
	v[4494] = v[563] + v[875];
	v[2350] = -0.5e0*v[563];
	v[905] = v[4498] + v[482] * v[832];
	v[903] = v[4482] * v[563];
	v[879] = v[4494] + v[479] * v[838];
	v[876] = v[4494] + v[478] * v[832];
	v[867] = v[4489] * v[563];
	v[4500] = v[564] + v[867];
	v[4495] = -v[564] + v[867];
	v[919] = v[479] * v[834] - v[844];
	v[910] = v[472] * v[834] - v[844];
	v[906] = v[4495] + v[482] * v[833];
	v[902] = v[479] * v[841] - v[844];
	v[893] = v[472] * v[841] - v[844];
	v[877] = v[4500] + v[478] * v[833];
	v[868] = v[4495] + v[473] * v[838];
	v[865] = 2e0*v[563];
	v[866] = v[473] * v[841] - v[865];
	v[863] = v[4486] * v[563];
	v[4499] = v[565] + v[863];
	v[4497] = -v[565] + v[863];
	v[4496] = 2e0*v[863];
	v[934] = v[4496] + v[482] * v[841];
	v[931] = v[4496] + v[478] * v[841];
	v[921] = v[4496] + v[482] * v[840];
	v[912] = v[4496] + v[478] * v[840];
	v[882] = v[4499] + v[479] * v[833];
	v[869] = v[4497] + v[473] * v[832];
	v[864] = v[4497] + v[472] * v[833];
	v[861] = v[4498] + v[472] * v[838];
	v[860] = v[4499] + v[471] * v[832];
	v[859] = v[4500] + v[471] * v[838];
	v[858] = v[471] * v[841] + v[865];
	v[856] = v[2350] * v[4483];
	v[852] = -2e0*v[4496] + v[4492] * v[838];
	v[485] = -v[4487] - v[4490];
	v[4501] = 0.5e0*v[485];
	v[942] = v[4501] * v[833] + 2e0*v[849];
	v[845] = v[4501] * v[832] + v[844];
	v[843] = -v[4496] + v[4501] * v[838];
	v[842] = v[4501] * v[841];
	v[566] = v[4501] * v[563];
	v[476] = -v[4484] - v[4490];
	v[4502] = 0.5e0*v[476];
	v[871] = v[4502] * v[832] + 2e0*v[844];
	v[850] = v[4502] * v[833] + v[849];
	v[848] = v[4502] * v[840];
	v[847] = -v[4496] + v[4502] * v[838];
	v[582] = v[4502] * v[564];
	v[4076] = v[4503] / 2e0;
	v[4073] = -v[4504] / 2e0;
	v[4052] = v[4505] / 2e0;
	v[3457] = v[4506] / 2e0;
	v[3454] = -v[4507] / 2e0;
	v[3432] = v[4508] / 2e0;
	v[301] = v[4510] * v[4511];
	v[337] = (*phi1) + cp[0] * v[4509];
	v[4538] = sin(v[337]);
	v[4539] = (*A2) + (*A12)*v[4538];
	v[4512] = cos(v[337]);
	v[1473] = -(Power(v[4509], 3)*v[4512]);
	v[717] = -((v[4509] * v[4509])*v[4538]);
	v[338] = v[4509] * v[4512];
	v[310] = (*phi1) + ci[0] * v[4509];
	v[311] = v[4509] * cos(v[310]);
	v[347] = (*phi2) + cp[1] * v[301];
	v[4540] = sin(v[347]);
	v[4541] = (*A1) + (*A12)*v[4540];
	v[4513] = cos(v[347]);
	v[1475] = -(Power(v[301], 3)*v[4513]);
	v[718] = -((v[301] * v[301])*v[4540]);
	v[1477] = (*A12)*v[338] * v[718];
	v[348] = v[301] * v[4513];
	v[4514] = (*A12)*v[348];
	v[1476] = v[4514] * v[717];
	v[376] = v[338] * v[4514];
	v[320] = (*phi2) + ci[1] * v[301];
	v[321] = v[301] * cos(v[320]);
	v[946] = e1i[2] * v[882] + e1i[1] * v[906] + e1i[0] * v[942];
	v[938] = e1i[0] * v[842] + e1i[2] * v[902] + e1i[1] * v[934];
	v[929] = e1i[2] * v[869] + e1i[1] * v[871] + e1i[0] * v[876];
	v[925] = e1i[2] * v[837] + e1i[0] * v[910] + e1i[1] * v[923];
	v[917] = e1i[1] * v[848] + e1i[2] * v[899] + e1i[0] * v[912];
	v[890] = e1i[2] * v[852] + e1i[1] * v[859] + e1i[0] * v[861];
	v[944] = e2i[2] * v[882] + e2i[1] * v[906] + e2i[0] * v[942];
	v[935] = e2i[0] * v[842] + e2i[2] * v[902] + e2i[1] * v[934];
	v[927] = e2i[2] * v[869] + e2i[1] * v[871] + e2i[0] * v[876];
	v[924] = e2i[2] * v[837] + e2i[0] * v[910] + e2i[1] * v[923];
	v[914] = e2i[1] * v[848] + e2i[2] * v[899] + e2i[0] * v[912];
	v[887] = e2i[2] * v[852] + e2i[1] * v[859] + e2i[0] * v[861];
	v[243] = 4e0 / v[489];
	v[3112] = (v[243] * v[243]);
	v[4790] = -(v[4499] / v[3112]);
	v[4788] = -(v[4495] / v[3112]);
	v[883] = -0.5e0*v[243];
	v[909] = v[883] - v[903];
	v[4521] = v[900] - v[909];
	v[920] = v[4521] + v[482] * v[838];
	v[939] = e1i[0] * v[843] + e1i[2] * v[879] + e1i[1] * v[920];
	v[936] = e2i[0] * v[843] + e2i[2] * v[879] + e2i[1] * v[920];
	v[911] = v[4521] + v[478] * v[838];
	v[916] = e1i[1] * v[847] + e1i[2] * v[868] + e1i[0] * v[911];
	v[913] = e2i[1] * v[847] + e2i[2] * v[868] + e2i[0] * v[911];
	v[884] = v[3108] + v[883];
	v[4523] = -v[884] + v[900];
	v[4522] = -v[884] + v[903];
	v[904] = v[4522] + v[479] * v[832];
	v[940] = e1i[0] * v[845] + e1i[2] * v[904] + e1i[1] * v[905];
	v[937] = e2i[0] * v[845] + e2i[2] * v[904] + e2i[1] * v[905];
	v[901] = v[4523] + v[473] * v[833];
	v[918] = e1i[1] * v[850] + e1i[0] * v[877] + e1i[2] * v[901];
	v[915] = e2i[1] * v[850] + e2i[0] * v[877] + e2i[2] * v[901];
	v[894] = v[4522] + v[472] * v[832];
	v[898] = e1i[2] * v[835] + e1i[1] * v[860] + e1i[0] * v[894];
	v[896] = e2i[2] * v[835] + e2i[1] * v[860] + e2i[0] * v[894];
	v[886] = v[4523] + v[471] * v[833];
	v[892] = e1i[2] * v[836] + e1i[0] * v[864] + e1i[1] * v[886];
	v[889] = e2i[2] * v[836] + e2i[0] * v[864] + e2i[1] * v[886];
	v[870] = -v[243] + v[872];
	v[4524] = v[870] + v[872];
	v[943] = v[4524] + v[4501] * v[834];
	v[947] = e1i[1] * v[908] + e1i[2] * v[919] + e1i[0] * v[943];
	v[945] = e2i[1] * v[908] + e2i[2] * v[919] + e2i[0] * v[943];
	v[873] = v[4524] + v[4502] * v[834];
	v[930] = e1i[1] * v[873] + e1i[0] * v[878] + e1i[2] * v[926];
	v[928] = e2i[1] * v[873] + e2i[0] * v[878] + e2i[2] * v[926];
	v[855] = -v[243] + v[856];
	v[4525] = v[855] + v[856];
	v[874] = v[4525] + v[4502] * v[841];
	v[933] = e1i[2] * v[866] + e1i[1] * v[874] + e1i[0] * v[931];
	v[932] = e2i[2] * v[866] + e2i[1] * v[874] + e2i[0] * v[931];
	v[857] = v[4525] + v[4492] * v[841];
	v[897] = e1i[2] * v[857] + e1i[1] * v[858] + e1i[0] * v[893];
	v[895] = e2i[2] * v[857] + e2i[1] * v[858] + e2i[0] * v[893];
	v[851] = -v[243] + v[853];
	v[4526] = v[851] + v[853];
	v[948] = v[4526] + v[4501] * v[840];
	v[950] = e1i[2] * v[881] + e1i[1] * v[921] + e1i[0] * v[948];
	v[949] = e2i[2] * v[881] + e2i[1] * v[921] + e2i[0] * v[948];
	v[854] = v[4526] + v[4492] * v[840];
	v[891] = e1i[2] * v[854] + e1i[0] * v[862] + e1i[1] * v[885];
	v[888] = e2i[2] * v[854] + e2i[0] * v[862] + e2i[1] * v[885];
	v[596] = v[4485] * v[883];
	v[597] = -(v[2351] * v[468]) + v[596];
	v[594] = v[4483] * v[883];
	v[595] = -(v[2350] * v[468]) + v[594];
	v[591] = v[243] + v[471] * v[563];
	v[589] = -v[243] + v[472] * v[564];
	v[585] = -v[243] + v[473] * v[563];
	v[583] = v[4488] * v[883];
	v[584] = v[4502] * v[565] + v[583];
	v[581] = v[4502] * v[563] + v[594];
	v[580] = v[243] + v[478] * v[565];
	v[576] = v[243] + v[479] * v[564];
	v[574] = -(v[243] * v[4489]);
	v[592] = v[471] * v[564] - v[574];
	v[615] = e2i[0] * v[589] + e2i[1] * v[592] + e2i[2] * v[597];
	v[606] = e1i[0] * v[589] + e1i[1] * v[592] + e1i[2] * v[597];
	v[588] = v[472] * v[563] - v[574];
	v[614] = e2i[0] * v[588] + e2i[1] * v[591] + e2i[2] * v[595];
	v[605] = e1i[0] * v[588] + e1i[1] * v[591] + e1i[2] * v[595];
	v[586] = v[473] * v[564] - v[574];
	v[575] = v[479] * v[563] - v[574];
	v[573] = -v[243] + v[482] * v[565];
	v[571] = -(v[243] * v[4482]);
	v[590] = v[472] * v[565] - v[571];
	v[579] = v[478] * v[564] - v[571];
	v[612] = e2i[0] * v[579] + e2i[1] * v[582] + e2i[2] * v[586];
	v[4530] = 2e0*v[612];
	v[603] = e1i[0] * v[579] + e1i[1] * v[582] + e1i[2] * v[586];
	v[4531] = 2e0*v[603];
	v[577] = v[479] * v[565] - v[571];
	v[572] = v[482] * v[564] - v[571];
	v[569] = v[243] * v[4486];
	v[593] = v[471] * v[565] + v[569];
	v[616] = e2i[0] * v[590] + e2i[1] * v[593] + e2i[2] * v[598];
	v[607] = e1i[0] * v[590] + e1i[1] * v[593] + e1i[2] * v[598];
	v[587] = v[473] * v[565] + v[569];
	v[613] = e2i[0] * v[580] + e2i[1] * v[584] + e2i[2] * v[587];
	v[4532] = 2e0*v[613];
	v[604] = e1i[0] * v[580] + e1i[1] * v[584] + e1i[2] * v[587];
	v[4533] = 2e0*v[604];
	v[578] = v[478] * v[563] + v[569];
	v[611] = e2i[0] * v[578] + e2i[1] * v[581] + e2i[2] * v[585];
	v[4528] = 2e0*v[611];
	v[602] = e1i[0] * v[578] + e1i[1] * v[581] + e1i[2] * v[585];
	v[4529] = 2e0*v[602];
	v[570] = v[482] * v[563] + v[569];
	v[608] = e2i[0] * v[566] + e2i[1] * v[570] + e2i[2] * v[575];
	v[599] = e1i[0] * v[566] + e1i[1] * v[570] + e1i[2] * v[575];
	v[568] = v[4501] * v[565] + v[583];
	v[610] = e2i[0] * v[568] + e2i[1] * v[573] + e2i[2] * v[577];
	v[601] = e1i[0] * v[568] + e1i[1] * v[573] + e1i[2] * v[577];
	v[567] = v[4501] * v[564] + v[596];
	v[609] = e2i[0] * v[567] + e2i[1] * v[572] + e2i[2] * v[576];
	v[600] = e1i[0] * v[567] + e1i[1] * v[572] + e1i[2] * v[576];
	v[246] = 1e0 - v[485] * v[883];
	v[247] = v[243] * v[482];
	v[248] = v[243] * v[479];
	v[250] = v[243] * v[478];
	v[252] = 1e0 - v[476] * v[883];
	v[253] = v[243] * v[473];
	v[255] = v[243] * v[472];
	v[257] = v[243] * v[471];
	v[258] = 1e0 - v[468] * v[883];
	v[262] = 4e0 / v[1192];
	v[3136] = (v[262] * v[262]);
	v[4798] = v[4544] / v[3136];
	v[2758] = -0.5e0*v[262];
	v[2759] = -v[2758] + v[2820];
	v[2755] = v[2758] + v[3132];
	v[4548] = v[2755] + v[2780];
	v[2754] = -v[262] + v[2817];
	v[2752] = -v[262] + v[2804];
	v[2750] = -v[262] + v[2799];
	v[1281] = v[2758] * v[4470];
	v[1277] = v[1193] * v[1276] + v[262];
	v[1273] = v[2758] * v[4468];
	v[1275] = v[1273] - v[1274] * v[4527];
	v[1265] = v[2758] * v[4473];
	v[1257] = v[1193] * v[1256] - v[262];
	v[1255] = v[1273] - v[1254] * v[4527];
	v[1244] = -(v[262] * v[4467]);
	v[1245] = v[1195] * v[1239] - v[1244];
	v[1238] = v[262] * v[4471];
	v[1240] = v[1238] + v[1193] * v[1239];
	v[1235] = -(v[262] * v[4474]);
	v[1283] = -v[1235] + v[1195] * v[1276];
	v[1237] = -v[1235] + v[1193] * v[1236];
	v[1168] = -(v[1274] * v[2758]);
	v[1166] = -(v[1254] * v[2758]);
	v[1164] = -(v[1233] * v[2758]);
	v[266] = v[1239] * v[262];
	v[267] = v[1236] * v[262];
	v[269] = v[1258] * v[262];
	v[271] = 1e0 + v[1166];
	v[272] = v[1256] * v[262];
	v[274] = v[1278] * v[262];
	v[276] = v[1276] * v[262];
	v[277] = 1e0 + v[1168];
	v[281] = e1i[0] * v[246] + e1i[1] * v[247] + e1i[2] * v[248];
	v[282] = e1i[0] * v[250] + e1i[1] * v[252] + e1i[2] * v[253];
	v[283] = e1i[0] * v[255] + e1i[1] * v[257] + e1i[2] * v[258];
	v[284] = e2i[0] * v[246] + e2i[1] * v[247] + e2i[2] * v[248];
	v[285] = e2i[0] * v[250] + e2i[1] * v[252] + e2i[2] * v[253];
	v[956] = v[4528] * v[599] - v[4529] * v[608] + v[281] * v[932] - v[284] * v[933] - v[282] * v[935] + v[285] * v[938];
	v[955] = v[4530] * v[600] - v[4531] * v[609] + v[281] * v[914] - v[284] * v[917] - v[282] * v[949] + v[285] * v[950];
	v[954] = -(v[603] * v[608]) - v[602] * v[609] + v[600] * v[611] + v[599] * v[612] + v[281] * v[913] - v[284] * v[916]
		- v[282] * v[936] + v[285] * v[939];
	v[953] = v[4532] * v[601] - v[4533] * v[610] + v[281] * v[928] - v[284] * v[930] - v[282] * v[945] + v[285] * v[947];
	v[952] = -(v[604] * v[609]) - v[603] * v[610] + v[601] * v[612] + v[600] * v[613] + v[281] * v[915] - v[284] * v[918]
		- v[282] * v[944] + v[285] * v[946];
	v[951] = -(v[604] * v[608]) - v[602] * v[610] + v[601] * v[611] + v[599] * v[613] + v[281] * v[927] - v[284] * v[929]
		- v[282] * v[937] + v[285] * v[940];
	v[625] = v[285] * v[601] - v[284] * v[604] - v[282] * v[610] + v[281] * v[613];
	v[721] = v[376] * v[625];
	v[624] = v[285] * v[600] - v[284] * v[603] - v[282] * v[609] + v[281] * v[612];
	v[725] = v[376] * v[624];
	v[623] = v[285] * v[599] - v[284] * v[602] - v[282] * v[608] + v[281] * v[611];
	v[729] = v[376] * v[623];
	v[286] = e2i[0] * v[255] + e2i[1] * v[257] + e2i[2] * v[258];
	v[968] = -(v[4528] * v[605]) + v[4529] * v[614] + v[282] * v[895] - v[285] * v[897] - v[283] * v[932] + v[286] * v[933];
	v[967] = -(v[4530] * v[606]) + v[4531] * v[615] + v[282] * v[888] - v[285] * v[891] - v[283] * v[914] + v[286] * v[917];
	v[966] = -(v[606] * v[611]) - v[605] * v[612] + v[603] * v[614] + v[602] * v[615] + v[282] * v[887] - v[285] * v[890]
		- v[283] * v[913] + v[286] * v[916];
	v[965] = -(v[4532] * v[607]) + v[4533] * v[616] + v[282] * v[924] - v[285] * v[925] - v[283] * v[928] + v[286] * v[930];
	v[964] = -(v[607] * v[612]) - v[606] * v[613] + v[604] * v[615] + v[603] * v[616] + v[282] * v[889] - v[285] * v[892]
		- v[283] * v[915] + v[286] * v[918];
	v[963] = -(v[607] * v[611]) - v[605] * v[613] + v[604] * v[614] + v[602] * v[616] + v[282] * v[896] - v[285] * v[898]
		- v[283] * v[927] + v[286] * v[929];
	v[962] = 2e0*v[605] * v[608] - 2e0*v[599] * v[614] - v[281] * v[895] + v[284] * v[897] + v[283] * v[935] - v[286] * v[938];
	v[961] = 2e0*v[606] * v[609] - 2e0*v[600] * v[615] - v[281] * v[888] + v[284] * v[891] + v[283] * v[949] - v[286] * v[950];
	v[960] = v[606] * v[608] + v[605] * v[609] - v[600] * v[614] - v[599] * v[615] - v[281] * v[887] + v[284] * v[890]
		+ v[283] * v[936] - v[286] * v[939];
	v[959] = 2e0*v[607] * v[610] - 2e0*v[601] * v[616] - v[281] * v[924] + v[284] * v[925] + v[283] * v[945] - v[286] * v[947];
	v[958] = v[607] * v[609] + v[606] * v[610] - v[601] * v[615] - v[600] * v[616] - v[281] * v[889] + v[284] * v[892]
		+ v[283] * v[944] - v[286] * v[946];
	v[957] = v[607] * v[608] + v[605] * v[610] - v[601] * v[614] - v[599] * v[616] - v[281] * v[896] + v[284] * v[898]
		+ v[283] * v[937] - v[286] * v[940];
	v[622] = -(v[286] * v[601]) + v[284] * v[607] + v[283] * v[610] - v[281] * v[616];
	v[733] = v[376] * v[622];
	v[621] = -(v[286] * v[600]) + v[284] * v[606] + v[283] * v[609] - v[281] * v[615];
	v[737] = v[376] * v[621];
	v[620] = -(v[286] * v[599]) + v[284] * v[605] + v[283] * v[608] - v[281] * v[614];
	v[741] = v[376] * v[620];
	v[619] = v[286] * v[604] - v[285] * v[607] - v[283] * v[613] + v[282] * v[616];
	v[745] = v[376] * v[619];
	v[618] = v[286] * v[603] - v[285] * v[606] - v[283] * v[612] + v[282] * v[615];
	v[749] = v[376] * v[618];
	v[617] = v[286] * v[602] - v[285] * v[605] - v[283] * v[611] + v[282] * v[614];
	v[753] = v[376] * v[617];
	v[287] = -(v[283] * v[285]) + v[282] * v[286];
	v[1479] = v[1477] * v[287];
	v[1478] = v[1476] * v[287];
	v[386] = v[287] * v[376];
	v[288] = v[283] * v[284] - v[281] * v[286];
	v[1481] = v[1477] * v[288];
	v[1480] = v[1476] * v[288];
	v[388] = v[288] * v[376];
	v[289] = -(v[282] * v[284]) + v[281] * v[285];
	v[1483] = v[1477] * v[289];
	v[1482] = v[1476] * v[289];
	v[390] = v[289] * v[376];
	v[1527] = 2e0*((v[386] * v[386]) + (v[388] * v[388]) + (v[390] * v[390]));
	v[293] = sin(v[310]);
	v[4537] = (*A2) + (*A12)*v[293];
	v[322] = v[321] * v[4537];
	v[294] = sin(v[320]);
	v[312] = ((*A1) + (*A12)*v[294])*v[311];
	v[296] = (*A1)*v[293] + v[294] * v[4537];
	v[2779] = ci[1] * v[935] + ci[0] * v[938] + v[296] * v[968];
	v[2778] = ci[1] * v[949] + ci[0] * v[950] + v[296] * v[967];
	v[2777] = ci[1] * v[936] + ci[0] * v[939] + v[296] * v[966];
	v[2776] = ci[1] * v[945] + ci[0] * v[947] + v[296] * v[965];
	v[2775] = ci[1] * v[944] + ci[0] * v[946] + v[296] * v[964];
	v[2774] = ci[1] * v[937] + ci[0] * v[940] + v[296] * v[963];
	v[2773] = ci[1] * v[932] + ci[0] * v[933] + v[296] * v[962];
	v[2772] = ci[1] * v[914] + ci[0] * v[917] + v[296] * v[961];
	v[2771] = ci[1] * v[913] + ci[0] * v[916] + v[296] * v[960];
	v[2770] = ci[1] * v[928] + ci[0] * v[930] + v[296] * v[959];
	v[2769] = ci[1] * v[915] + ci[0] * v[918] + v[296] * v[958];
	v[2768] = ci[1] * v[927] + ci[0] * v[929] + v[296] * v[957];
	v[2767] = ci[1] * v[895] + ci[0] * v[897] + v[296] * v[956];
	v[2766] = ci[1] * v[888] + ci[0] * v[891] + v[296] * v[955];
	v[2765] = ci[1] * v[887] + ci[0] * v[890] + v[296] * v[954];
	v[2764] = ci[1] * v[924] + ci[0] * v[925] + v[296] * v[953];
	v[2763] = ci[1] * v[889] + ci[0] * v[892] + v[296] * v[952];
	v[2762] = ci[1] * v[896] + ci[0] * v[898] + v[296] * v[951];
	v[1292] = ci[0] * v[607] + ci[1] * v[616] + v[296] * v[625];
	v[1291] = ci[0] * v[606] + ci[1] * v[615] + v[296] * v[624];
	v[1290] = ci[0] * v[605] + ci[1] * v[614] + v[296] * v[623];
	v[1272] = ci[0] * v[604] + ci[1] * v[613] + v[296] * v[622];
	v[1271] = ci[0] * v[603] + ci[1] * v[612] + v[296] * v[621];
	v[1270] = ci[0] * v[602] + ci[1] * v[611] + v[296] * v[620];
	v[1253] = ci[0] * v[601] + ci[1] * v[610] + v[296] * v[619];
	v[1252] = ci[0] * v[600] + ci[1] * v[609] + v[296] * v[618];
	v[1251] = ci[0] * v[599] + ci[1] * v[608] + v[296] * v[617];
	v[1484] = v[1475] * v[4539];
	v[1487] = v[1484] * v[287];
	v[1486] = v[1484] * v[288];
	v[1485] = v[1484] * v[289];
	v[377] = v[4539] * v[718];
	v[398] = v[289] * v[377];
	v[397] = v[288] * v[377];
	v[396] = v[287] * v[377];
	v[349] = v[348] * v[4539];
	v[652] = v[616] + v[349] * v[625];
	v[651] = v[615] + v[349] * v[624];
	v[650] = v[614] + v[349] * v[623];
	v[649] = v[613] + v[349] * v[622];
	v[648] = v[612] + v[349] * v[621];
	v[647] = v[611] + v[349] * v[620];
	v[646] = v[610] + v[349] * v[619];
	v[645] = v[609] + v[349] * v[618];
	v[644] = v[608] + v[349] * v[617];
	v[384] = v[286] + v[289] * v[349];
	v[4572] = 2e0*v[384];
	v[382] = v[285] + v[288] * v[349];
	v[4571] = 2e0*v[382];
	v[380] = v[284] + v[287] * v[349];
	v[4570] = 2e0*v[380];
	v[1488] = v[1473] * v[4541];
	v[1491] = v[1488] * v[287];
	v[1490] = v[1488] * v[288];
	v[1489] = v[1488] * v[289];
	v[378] = v[4541] * v[717];
	v[389] = v[289] * v[378];
	v[387] = v[288] * v[378];
	v[385] = v[287] * v[378];
	v[339] = v[338] * v[4541];
	v[634] = v[607] + v[339] * v[625];
	v[633] = v[606] + v[339] * v[624];
	v[632] = v[605] + v[339] * v[623];
	v[631] = v[604] + v[339] * v[622];
	v[630] = v[603] + v[339] * v[621];
	v[629] = v[602] + v[339] * v[620];
	v[628] = v[601] + v[339] * v[619];
	v[627] = v[600] + v[339] * v[618];
	v[626] = v[599] + v[339] * v[617];
	v[383] = v[283] + v[289] * v[339];
	v[4566] = 2e0*v[383];
	v[381] = v[282] + v[288] * v[339];
	v[4567] = 2e0*v[381];
	v[379] = v[281] + v[287] * v[339];
	v[4568] = 2e0*v[379];
	v[304] = (*A1)*v[4538] + v[4539] * v[4540];
	v[688] = cp[0] * v[607] + cp[1] * v[616] + v[304] * v[625];
	v[686] = cp[0] * v[606] + cp[1] * v[615] + v[304] * v[624];
	v[684] = cp[0] * v[605] + cp[1] * v[614] + v[304] * v[623];
	v[682] = cp[0] * v[604] + cp[1] * v[613] + v[304] * v[622];
	v[680] = cp[0] * v[603] + cp[1] * v[612] + v[304] * v[621];
	v[678] = cp[0] * v[602] + cp[1] * v[611] + v[304] * v[620];
	v[676] = cp[0] * v[601] + cp[1] * v[610] + v[304] * v[619];
	v[674] = cp[0] * v[600] + cp[1] * v[609] + v[304] * v[618];
	v[672] = cp[0] * v[599] + cp[1] * v[608] + v[304] * v[617];
	v[303] = cp[0] * v[281] + cp[1] * v[284] + v[287] * v[304] + v[4515];
	v[368] = -v[303] + v[4518];
	v[305] = cp[0] * v[282] + cp[1] * v[285] + v[288] * v[304] + v[4516];
	v[369] = -v[305] + v[4519];
	v[306] = cp[0] * v[283] + cp[1] * v[286] + v[289] * v[304] + v[4517];
	v[370] = -v[306] + v[4520];
	v[313] = e1i[0] + v[312] * v[4534];
	v[314] = e1i[1] + v[312] * v[4535];
	v[315] = e1i[2] + v[312] * v[4536];
	v[317] = 1e0 / sqrt(Power(v[313], 2) + Power(v[314], 2) + Power(v[315], 2));
	v[316] = v[313] * v[317];
	v[318] = v[314] * v[317];
	v[319] = v[315] * v[317];
	v[323] = e2i[0] + v[322] * v[4534];
	v[324] = e2i[1] + v[322] * v[4535];
	v[325] = e2i[2] + v[322] * v[4536];
	v[327] = 1e0 / sqrt(Power(v[323], 2) + Power(v[324], 2) + Power(v[325], 2));
	v[326] = v[323] * v[327];
	v[328] = v[324] * v[327];
	v[335] = -(v[318] * v[326]) + v[316] * v[328];
	v[329] = v[325] * v[327];
	v[333] = v[319] * v[326] - v[316] * v[329];
	v[330] = -(v[319] * v[328]) + v[318] * v[329];
	v[332] = 1e0 / sqrt(Power(v[330], 2) + Power(v[333], 2) + Power(v[335], 2));
	v[331] = v[330] * v[332];
	v[334] = v[332] * v[333];
	v[336] = v[332] * v[335];
	v[2860] = (*radius)*(-((v[1239] * v[2744] - v[2823])*v[334]) - (v[1236] * v[2744] - v[2753])*v[336]
		- v[2744] * v[331] * v[4542]);
	v[2856] = (*radius)*(-((v[1239] * v[2743] - v[2823])*v[334]) - (v[1236] * v[2743] - v[2853])*v[336] - v[331] *
		(v[2750] + v[2799] + v[2743] * v[4542]));
	v[2851] = (*radius)*(-((v[1239] * v[2741] - v[2755] + v[2820])*v[334]) - (v[1193] + v[1236] * v[2741] + v[2760]
		)*v[336] - v[331] * (v[2823] + v[2741] * v[4542]));
	v[2847] = (*radius)*(-((v[1239] * v[2740] - v[2843])*v[334]) - (v[1236] * v[2740] - v[2753])*v[336] - v[331] *
		(v[2754] + v[2817] + v[2740] * v[4542]));
	v[2842] = (*radius)*(-((v[1196] + v[1236] * v[2739] - v[2756])*v[336]) + v[331] * (-2e0*v[2751] - v[2739] * v[4542]
		) + v[334] * (-(v[1239] * v[2739]) + v[4544]));
	v[2838] = (*radius)*(v[331] * (-v[2753] - v[2738] * v[4542]) + v[334] * (-(v[1239] * v[2738]) + v[4546]) + v[336] * (-
		(v[1236] * v[2738]) + v[4548]));
	v[2834] = (*radius)*(-((v[1258] * v[2744] - v[2823])*v[331]) - (v[1256] * v[2744] - v[2831])*v[336] - v[334] *
		(v[2752] + v[2804] + v[2744] * v[4543]));
	v[2829] = (*radius)*(-((v[1258] * v[2743] - v[2823])*v[331]) - (v[1256] * v[2743] - v[2751])*v[336]
		- v[2743] * v[334] * v[4543]);
	v[2825] = (*radius)*(-((v[1258] * v[2741] - v[2755] + v[2820])*v[331]) - v[334] * (v[2823] + v[2741] * v[4543])
		+ v[336] * (-(v[1256] * v[2741]) + v[4544]));
	v[2819] = (*radius)*(-((v[1258] * v[2740] + v[2843])*v[331]) - (v[1256] * v[2740] - v[2751])*v[336] - v[334] *
		(v[2754] + v[2817] + v[2740] * v[4543]));
	v[2814] = (*radius)*(-((v[1195] + v[1258] * v[2739] - v[2757])*v[331]) - (v[1256] * v[2739] + v[2759] - v[2780]
		)*v[336] - v[334] * (v[2751] + v[2739] * v[4543]));
	v[2810] = (*radius)*(-((v[1193] + v[1258] * v[2738] + v[2760])*v[331]) + v[334] * (-2e0*v[2753] - v[2738] * v[4543]
		) + v[336] * (-(v[1256] * v[2738]) + v[4547]));
	v[2806] = (*radius)*(-((v[1278] * v[2744] - v[2753])*v[331]) - (v[1276] * v[2744] + v[2831])*v[334] - v[336] *
		(v[2752] + v[2804] + v[2744] * v[4545]));
	v[2801] = (*radius)*(-((v[1278] * v[2743] + v[2853])*v[331]) - (v[1276] * v[2743] - v[2751])*v[334] - v[336] *
		(v[2750] + v[2799] + v[2743] * v[4545]));
	v[2796] = (*radius)*(-((v[1195] + v[1276] * v[2741] - v[2757])*v[334]) + v[336] * (-2e0*v[2823] - v[2741] * v[4545]
		) + v[331] * (-(v[1278] * v[2741]) + v[4546]));
	v[2792] = (*radius)*(-((v[1278] * v[2740] - v[2753])*v[331]) - (v[1276] * v[2740] - v[2751])*v[334]
		- v[2740] * v[336] * v[4545]);
	v[2788] = (*radius)*(-((v[1276] * v[2739] + v[2759] - v[2780])*v[334]) - v[336] * (v[2751] + v[2739] * v[4545])
		+ v[331] * (-(v[1278] * v[2739]) + v[4547]));
	v[2784] = (*radius)*(-((v[1196] + v[1276] * v[2738] - v[2756])*v[334]) + v[336] * (-v[2753] - v[2738] * v[4545])
		+ v[331] * (-(v[1278] * v[2738]) + v[4548]));
	v[1289] = (*radius)*((v[1244] - v[1196] * v[1278])*v[331] - (v[1238] + v[1196] * v[1276])*v[334]
		+ v[1274] * v[336] * v[4479]);
	v[1285] = (*radius)*((-(v[1195] * v[1278]) + v[262])*v[331] - v[1283] * v[334] - v[336] * (v[1281]
		- v[1274] * v[4481]));
	v[1280] = -((*radius)*((-v[1235] + v[1193] * v[1278])*v[331] + v[1277] * v[334] + v[1275] * v[336]));
	v[1269] = (*radius)*(-((v[1196] * v[1258] + v[262])*v[331]) - (v[1238] + v[1196] * v[1256])*v[336] - v[334] *
		(v[1265] - v[1254] * v[4479]));
	v[1264] = (*radius)*((v[1244] - v[1195] * v[1258])*v[331] - v[1261] * v[334] + (v[1235] - v[1195] * v[1256])*v[336]
		);
	v[1260] = -((*radius)*((v[1238] + v[1193] * v[1258])*v[331] + v[1255] * v[334] + v[1257] * v[336]));
	v[1250] = (*radius)*((-(v[1196] * v[1239]) + v[262])*v[334] - (v[1196] * v[1236] - v[1244])*v[336] - v[331] *
		(v[1265] - v[1233] * v[4479]));
	v[1246] = (*radius)*(-(v[1245] * v[334]) - (v[1195] * v[1236] + v[262])*v[336] - v[331] * (v[1281]
		- v[1233] * v[4481]));
	v[1241] = (*radius)*(-(v[1240] * v[334]) - v[1237] * v[336] + v[1233] * v[331] * v[4527]);
	v[526] = 2e0*(v[379] * v[386] + v[381] * v[388] + v[383] * v[390]);
	v[525] = 2e0*(v[379] * v[385] + v[381] * v[387] + v[383] * v[389]);
	v[392] = (v[379] * v[379]) + (v[381] * v[381]) + (v[383] * v[383]);
	v[4550] = 1e0 / Power(v[392], 2);
	v[4195] = 1e0 / sqrt(v[392]);
	v[4556] = v[378] * v[4195];
	v[4551] = -(v[392] * v[4195]);
	v[4549] = v[4195] / 2e0;
	v[4197] = -v[4195] / (4e0*v[392]);
	v[4621] = 4e0*v[4197];
	v[4553] = 4e0*v[4197];
	v[4552] = -(v[392] * v[4197]);
	v[3550] = v[4549] * v[526];
	v[3548] = v[4549] * v[525];
	v[1500] = v[4550] * ((v[1491] * v[379] + v[1490] * v[381] + v[1489] * v[383] + (v[385] * v[385]) + (v[387] * v[387]) +
		(v[389] * v[389]))*v[4551] + v[3548] * v[525] + v[4552] * (v[525] * v[525]));
	v[1499] = v[4550] * (v[3550] * v[526] - v[392] * (v[4549] * (v[1527] + v[1483] * v[4566] + v[1481] * v[4567]
		+ v[1479] * v[4568]) + v[4197] * (v[526] * v[526])));
	v[1497] = v[4550] * ((v[1478] * v[379] + v[1480] * v[381] + v[1482] * v[383] + v[385] * v[386] + v[387] * v[388]
		+ v[389] * v[390])*v[4551] + v[525] * (v[3550] + v[4552] * v[526]));
	v[528] = -(v[3550] / v[392]);
	v[4555] = 2e0*v[528];
	v[527] = -(v[3548] / v[392]);
	v[4554] = 2e0*v[527];
	v[2025] = v[4553] * (v[368] * v[628] + v[369] * v[631] + v[370] * v[634] - v[379] * v[676] - v[381] * v[682]
		- v[383] * v[688]);
	v[2023] = v[4553] * (v[368] * v[627] + v[369] * v[630] + v[370] * v[633] - v[379] * v[674] - v[381] * v[680]
		- v[383] * v[686]);
	v[2021] = v[4621] * (v[368] * v[626] + v[369] * v[629] + v[370] * v[632] - v[379] * v[672] - v[381] * v[678]
		- v[383] * v[684]);
	v[1522] = v[1500] * v[383] + v[1489] * v[4195] + v[389] * v[4554];
	v[1520] = v[1500] * v[381] + v[1490] * v[4195] + v[387] * v[4554];
	v[1518] = v[1500] * v[379] + v[1491] * v[4195] + v[385] * v[4554];
	v[1524] = v[1518] * v[316] + v[1520] * v[318] + v[1522] * v[319];
	v[1523] = v[1518] * v[331] + v[1520] * v[334] + v[1522] * v[336];
	v[1512] = v[1499] * v[383] + v[1483] * v[4195] + v[390] * v[4555];
	v[1510] = v[1497] * v[383] + v[1482] * v[4195] + v[389] * v[4555];
	v[1508] = v[1499] * v[381] + v[1481] * v[4195] + v[388] * v[4555];
	v[1506] = v[1497] * v[381] + v[1480] * v[4195] + v[387] * v[4555];
	v[1504] = v[1499] * v[379] + v[1479] * v[4195] + v[386] * v[4555];
	v[1516] = v[1504] * v[316] + v[1508] * v[318] + v[1512] * v[319];
	v[1514] = v[1504] * v[331] + v[1508] * v[334] + v[1512] * v[336];
	v[1502] = v[1497] * v[379] + v[1478] * v[4195] + v[385] * v[4555];
	v[1515] = v[1502] * v[316] + v[1506] * v[318] + v[1510] * v[319];
	v[1513] = v[1502] * v[331] + v[1506] * v[334] + v[1510] * v[336];
	v[1004] = v[4195] * (v[938] + v[339] * v[968]);
	v[1002] = v[4195] * (v[950] + v[339] * v[967]);
	v[1000] = v[4195] * (v[939] + v[339] * v[966]);
	v[998] = v[4195] * (v[947] + v[339] * v[965]);
	v[996] = v[4195] * (v[946] + v[339] * v[964]);
	v[994] = v[4195] * (v[940] + v[339] * v[963]);
	v[992] = v[4195] * (v[933] + v[339] * v[962]);
	v[990] = v[4195] * (v[917] + v[339] * v[961]);
	v[988] = v[4195] * (v[916] + v[339] * v[960]);
	v[986] = v[4195] * (v[930] + v[339] * v[959]);
	v[984] = v[4195] * (v[918] + v[339] * v[958]);
	v[982] = v[4195] * (v[929] + v[339] * v[957]);
	v[980] = v[4195] * (v[897] + v[339] * v[956]);
	v[2872] = v[1004] * v[316] + v[319] * v[980] + v[318] * v[992];
	v[2871] = v[1004] * v[331] + v[336] * v[980] + v[334] * v[992];
	v[978] = v[4195] * (v[891] + v[339] * v[955]);
	v[2870] = v[1002] * v[316] + v[319] * v[978] + v[318] * v[990];
	v[2868] = v[1002] * v[331] + v[336] * v[978] + v[334] * v[990];
	v[976] = v[4195] * (v[890] + v[339] * v[954]);
	v[2869] = v[1000] * v[316] + v[319] * v[976] + v[318] * v[988];
	v[2867] = v[1000] * v[331] + v[336] * v[976] + v[334] * v[988];
	v[974] = v[4195] * (v[925] + v[339] * v[953]);
	v[2866] = v[319] * v[974] + v[318] * v[986] + v[316] * v[998];
	v[2863] = v[336] * v[974] + v[334] * v[986] + v[331] * v[998];
	v[972] = v[4195] * (v[892] + v[339] * v[952]);
	v[2865] = v[319] * v[972] + v[318] * v[984] + v[316] * v[996];
	v[2862] = v[336] * v[972] + v[334] * v[984] + v[331] * v[996];
	v[970] = v[4195] * (v[898] + v[339] * v[951]);
	v[2864] = v[319] * v[970] + v[318] * v[982] + v[316] * v[994];
	v[2861] = v[336] * v[970] + v[334] * v[982] + v[331] * v[994];
	v[754] = v[528] * v[626] + v[4195] * v[753];
	v[752] = v[4556] * v[617] + v[527] * v[626];
	v[750] = v[528] * v[627] + v[4195] * v[749];
	v[748] = v[4556] * v[618] + v[527] * v[627];
	v[746] = v[528] * v[628] + v[4195] * v[745];
	v[744] = v[4556] * v[619] + v[527] * v[628];
	v[742] = v[528] * v[629] + v[4195] * v[741];
	v[740] = v[4556] * v[620] + v[527] * v[629];
	v[738] = v[528] * v[630] + v[4195] * v[737];
	v[736] = v[4556] * v[621] + v[527] * v[630];
	v[734] = v[528] * v[631] + v[4195] * v[733];
	v[732] = v[4556] * v[622] + v[527] * v[631];
	v[730] = v[528] * v[632] + v[4195] * v[729];
	v[1918] = v[319] * v[730] + v[318] * v[742] + v[316] * v[754];
	v[1915] = v[336] * v[730] + v[334] * v[742] + v[331] * v[754];
	v[728] = v[4556] * v[623] + v[527] * v[632];
	v[1924] = v[319] * v[728] + v[318] * v[740] + v[316] * v[752];
	v[1921] = v[336] * v[728] + v[334] * v[740] + v[331] * v[752];
	v[726] = v[528] * v[633] + v[4195] * v[725];
	v[1919] = v[319] * v[726] + v[318] * v[738] + v[316] * v[750];
	v[1916] = v[336] * v[726] + v[334] * v[738] + v[331] * v[750];
	v[724] = v[4556] * v[624] + v[527] * v[633];
	v[1925] = v[319] * v[724] + v[318] * v[736] + v[316] * v[748];
	v[1922] = v[336] * v[724] + v[334] * v[736] + v[331] * v[748];
	v[722] = v[528] * v[634] + v[4195] * v[721];
	v[1920] = v[319] * v[722] + v[318] * v[734] + v[316] * v[746];
	v[1917] = v[336] * v[722] + v[334] * v[734] + v[331] * v[746];
	v[720] = v[4556] * v[625] + v[527] * v[634];
	v[1926] = v[319] * v[720] + v[318] * v[732] + v[316] * v[744];
	v[1923] = v[336] * v[720] + v[334] * v[732] + v[331] * v[744];
	v[643] = v[4195] * v[634];
	v[4615] = 2e0*v[643];
	v[642] = v[4195] * v[633];
	v[4614] = 2e0*v[642];
	v[641] = v[4195] * v[632];
	v[4613] = 2e0*v[641];
	v[640] = v[4195] * v[631];
	v[4611] = 2e0*v[640];
	v[639] = v[4195] * v[630];
	v[4610] = 2e0*v[639];
	v[638] = v[4195] * v[629];
	v[4609] = 2e0*v[638];
	v[637] = v[4195] * v[628];
	v[4602] = 2e0*v[637];
	v[1226] = v[331] * v[637] + v[334] * v[640] + v[336] * v[643];
	v[1208] = v[316] * v[637] + v[318] * v[640] + v[319] * v[643];
	v[636] = v[4195] * v[627];
	v[4604] = 2e0*v[636];
	v[1225] = v[331] * v[636] + v[334] * v[639] + v[336] * v[642];
	v[1207] = v[316] * v[636] + v[318] * v[639] + v[319] * v[642];
	v[635] = v[4195] * v[626];
	v[4606] = 2e0*v[635];
	v[1224] = v[331] * v[635] + v[334] * v[638] + v[336] * v[641];
	v[1206] = v[316] * v[635] + v[318] * v[638] + v[319] * v[641];
	v[409] = v[386] * v[4195] + v[379] * v[528];
	v[408] = v[388] * v[4195] + v[381] * v[528];
	v[4619] = 2e0*v[408];
	v[407] = v[390] * v[4195] + v[383] * v[528];
	v[1369] = v[336] * v[407] + v[334] * v[408] + v[331] * v[409];
	v[1357] = v[319] * v[407] + v[318] * v[408] + v[316] * v[409];
	v[406] = v[385] * v[4195] + v[379] * v[527];
	v[405] = v[387] * v[4195] + v[381] * v[527];
	v[4617] = 2e0*v[405];
	v[404] = v[389] * v[4195] + v[383] * v[527];
	v[1368] = v[336] * v[404] + v[334] * v[405] + v[331] * v[406];
	v[1356] = v[319] * v[404] + v[318] * v[405] + v[316] * v[406];
	v[343] = v[379] * v[4195];
	v[4563] = -(v[343] * v[626]);
	v[4560] = -(v[343] * v[627]);
	v[4557] = -(v[343] * v[628]);
	v[345] = v[381] * v[4195];
	v[4564] = -(v[345] * v[629]);
	v[4561] = -(v[345] * v[630]);
	v[4558] = -(v[345] * v[631]);
	v[346] = v[383] * v[4195];
	v[4565] = -(v[346] * v[632]);
	v[4562] = -(v[346] * v[633]);
	v[4559] = -(v[346] * v[634]);
	v[2004] = -(v[380] * v[637]) - v[382] * v[640] - v[384] * v[643] - v[343] * v[646] - v[345] * v[649] - v[346] * v[652]
		- v[409] * v[676] - v[408] * v[682] - v[407] * v[688] + v[370] * v[722] + v[369] * v[734] + v[368] * v[746];
	v[2003] = -(v[380] * v[636]) - v[382] * v[639] - v[384] * v[642] - v[343] * v[645] - v[345] * v[648] - v[346] * v[651]
		- v[409] * v[674] - v[408] * v[680] - v[407] * v[686] + v[370] * v[726] + v[369] * v[738] + v[368] * v[750];
	v[2002] = -(v[380] * v[635]) - v[382] * v[638] - v[384] * v[641] - v[343] * v[644] - v[345] * v[647] - v[346] * v[650]
		- v[409] * v[672] - v[408] * v[678] - v[407] * v[684] + v[370] * v[730] + v[369] * v[742] + v[368] * v[754];
	v[1998] = 2e0*v[4557] + 2e0*v[4558] + 2e0*v[4559] - v[406] * v[676] - v[405] * v[682] - v[404] * v[688] + v[370] * v[720]
		+ v[369] * v[732] + v[368] * v[744];
	v[1997] = 2e0*v[4560] + 2e0*v[4561] + 2e0*v[4562] - v[406] * v[674] - v[405] * v[680] - v[404] * v[686] + v[370] * v[724]
		+ v[369] * v[736] + v[368] * v[748];
	v[1996] = 2e0*v[4563] + 2e0*v[4564] + 2e0*v[4565] - v[406] * v[672] - v[405] * v[678] - v[404] * v[684] + v[370] * v[728]
		+ v[369] * v[740] + v[368] * v[752];
	v[1611] = v[1504] * v[368] + v[1508] * v[369] + v[1512] * v[370] - v[343] * v[396] - v[345] * v[397] - v[346] * v[398]
		- v[409] * v[4570] - v[408] * v[4571] - v[407] * v[4572];
	v[1608] = v[1502] * v[368] + v[1506] * v[369] + v[1510] * v[370] - v[343] * v[386] - v[345] * v[388] - v[346] * v[390]
		- v[384] * v[404] - v[382] * v[405] - v[380] * v[406] - v[383] * v[407] - v[381] * v[408] - v[379] * v[409];
	v[1607] = v[1518] * v[368] + v[1520] * v[369] + v[1522] * v[370] - v[343] * v[385] - v[345] * v[387] - v[346] * v[389]
		- v[404] * v[4566] - v[405] * v[4567] - v[406] * v[4568];
	v[530] = 2e0*(v[380] * v[396] + v[382] * v[397] + v[384] * v[398]);
	v[529] = 2e0*(v[380] * v[386] + v[382] * v[388] + v[384] * v[390]);
	v[400] = (v[380] * v[380]) + (v[382] * v[382]) + (v[384] * v[384]);
	v[4573] = 1e0 / Power(v[400], 2);
	v[4205] = 1e0 / sqrt(v[400]);
	v[4579] = v[377] * v[4205];
	v[4574] = -(v[400] * v[4205]);
	v[4569] = v[4205] / 2e0;
	v[4207] = -v[4205] / (4e0*v[400]);
	v[4620] = 4e0*v[4207];
	v[4576] = 4e0*v[4207];
	v[4575] = -(v[400] * v[4207]);
	v[3580] = v[4569] * v[530];
	v[3578] = v[4569] * v[529];
	v[1534] = v[4573] * (v[3578] * v[529] - v[400] * (v[4569] * (v[1527] + v[1478] * v[4570] + v[1480] * v[4571]
		+ v[1482] * v[4572]) + v[4207] * (v[529] * v[529])));
	v[1533] = v[4573] * ((v[1487] * v[380] + v[1486] * v[382] + v[1485] * v[384] + (v[396] * v[396]) + (v[397] * v[397]) +
		(v[398] * v[398]))*v[4574] + v[3580] * v[530] + v[4575] * (v[530] * v[530]));
	v[1531] = v[4573] * ((v[1479] * v[380] + v[1481] * v[382] + v[1483] * v[384] + v[386] * v[396] + v[388] * v[397]
		+ v[390] * v[398])*v[4574] + v[529] * (v[3580] + v[4575] * v[530]));
	v[532] = -(v[3580] / v[400]);
	v[4578] = 2e0*v[532];
	v[531] = -(v[3578] / v[400]);
	v[4577] = 2e0*v[531];
	v[2016] = v[4576] * (v[368] * v[646] + v[369] * v[649] + v[370] * v[652] - v[380] * v[676] - v[382] * v[682]
		- v[384] * v[688]);
	v[2014] = v[4576] * (v[368] * v[645] + v[369] * v[648] + v[370] * v[651] - v[380] * v[674] - v[382] * v[680]
		- v[384] * v[686]);
	v[2012] = v[4620] * (v[368] * v[644] + v[369] * v[647] + v[370] * v[650] - v[380] * v[672] - v[382] * v[678]
		- v[384] * v[684]);
	v[1552] = v[1534] * v[384] + v[1482] * v[4205] + v[390] * v[4577];
	v[1550] = v[1534] * v[382] + v[1480] * v[4205] + v[388] * v[4577];
	v[1548] = v[1534] * v[380] + v[1478] * v[4205] + v[386] * v[4577];
	v[1546] = v[1533] * v[384] + v[1485] * v[4205] + v[398] * v[4578];
	v[1544] = v[1531] * v[384] + v[1483] * v[4205] + v[390] * v[4578];
	v[1542] = v[1533] * v[382] + v[1486] * v[4205] + v[397] * v[4578];
	v[1540] = v[1531] * v[382] + v[1481] * v[4205] + v[388] * v[4578];
	v[1538] = v[1533] * v[380] + v[1487] * v[4205] + v[396] * v[4578];
	v[1536] = v[1531] * v[380] + v[1479] * v[4205] + v[386] * v[4578];
	v[1040] = v[4205] * (v[935] + v[349] * v[968]);
	v[1038] = v[4205] * (v[949] + v[349] * v[967]);
	v[1036] = v[4205] * (v[936] + v[349] * v[966]);
	v[1034] = v[4205] * (v[945] + v[349] * v[965]);
	v[1032] = v[4205] * (v[944] + v[349] * v[964]);
	v[1030] = v[4205] * (v[937] + v[349] * v[963]);
	v[1028] = v[4205] * (v[932] + v[349] * v[962]);
	v[1026] = v[4205] * (v[914] + v[349] * v[961]);
	v[1024] = v[4205] * (v[913] + v[349] * v[960]);
	v[1022] = v[4205] * (v[928] + v[349] * v[959]);
	v[1020] = v[4205] * (v[915] + v[349] * v[958]);
	v[1018] = v[4205] * (v[927] + v[349] * v[957]);
	v[1016] = v[4205] * (v[895] + v[349] * v[956]);
	v[1014] = v[4205] * (v[888] + v[349] * v[955]);
	v[1012] = v[4205] * (v[887] + v[349] * v[954]);
	v[1010] = v[4205] * (v[924] + v[349] * v[953]);
	v[1008] = v[4205] * (v[889] + v[349] * v[952]);
	v[1006] = v[4205] * (v[896] + v[349] * v[951]);
	v[781] = v[4579] * v[617] + v[532] * v[644];
	v[779] = v[531] * v[644] + v[4205] * v[753];
	v[778] = v[4579] * v[618] + v[532] * v[645];
	v[776] = v[531] * v[645] + v[4205] * v[749];
	v[775] = v[4579] * v[619] + v[532] * v[646];
	v[773] = v[531] * v[646] + v[4205] * v[745];
	v[772] = v[4579] * v[620] + v[532] * v[647];
	v[770] = v[531] * v[647] + v[4205] * v[741];
	v[769] = v[4579] * v[621] + v[532] * v[648];
	v[767] = v[531] * v[648] + v[4205] * v[737];
	v[766] = v[4579] * v[622] + v[532] * v[649];
	v[764] = v[531] * v[649] + v[4205] * v[733];
	v[763] = v[4579] * v[623] + v[532] * v[650];
	v[761] = v[531] * v[650] + v[4205] * v[729];
	v[760] = v[4579] * v[624] + v[532] * v[651];
	v[758] = v[531] * v[651] + v[4205] * v[725];
	v[757] = v[4579] * v[625] + v[532] * v[652];
	v[755] = v[531] * v[652] + v[4205] * v[721];
	v[664] = v[4205] * v[652];
	v[663] = v[4205] * v[651];
	v[662] = v[4205] * v[650];
	v[658] = v[4205] * v[649];
	v[657] = v[4205] * v[648];
	v[656] = v[4205] * v[647];
	v[655] = v[4205] * v[646];
	v[4601] = 2e0*v[655];
	v[654] = v[4205] * v[645];
	v[4603] = 2e0*v[654];
	v[653] = v[4205] * v[644];
	v[4605] = 2e0*v[653];
	v[415] = v[396] * v[4205] + v[380] * v[532];
	v[4590] = -2e0*v[415];
	v[414] = v[397] * v[4205] + v[382] * v[532];
	v[4589] = 2e0*v[414];
	v[413] = v[398] * v[4205] + v[384] * v[532];
	v[4594] = 2e0*v[413];
	v[412] = v[386] * v[4205] + v[380] * v[531];
	v[4592] = -2e0*v[412];
	v[411] = v[388] * v[4205] + v[382] * v[531];
	v[4591] = 2e0*v[411];
	v[410] = v[390] * v[4205] + v[384] * v[531];
	v[4593] = 2e0*v[410];
	v[353] = v[380] * v[4205];
	v[4586] = -(v[353] * v[644]);
	v[4583] = -(v[353] * v[645]);
	v[4580] = -(v[353] * v[646]);
	v[355] = v[382] * v[4205];
	v[4587] = -(v[355] * v[647]);
	v[4584] = -(v[355] * v[648]);
	v[4581] = -(v[355] * v[649]);
	v[1555] = v[1550] * v[343] - v[1548] * v[345] - v[1520] * v[353] + v[1518] * v[355] + v[406] * v[4591] + v[405] * v[4592];
	v[1554] = v[1542] * v[343] - v[1538] * v[345] - v[1508] * v[353] + v[1504] * v[355] + v[409] * v[4589] + v[408] * v[4590];
	v[1553] = v[1540] * v[343] - v[1536] * v[345] - v[1506] * v[353] + v[1502] * v[355] + v[409] * v[411] - v[408] * v[412]
		+ v[406] * v[414] - v[405] * v[415];
	v[661] = v[355] * v[637] - v[353] * v[640] - v[345] * v[655] + v[343] * v[658];
	v[660] = v[355] * v[636] - v[353] * v[639] - v[345] * v[654] + v[343] * v[657];
	v[659] = v[355] * v[635] - v[353] * v[638] - v[345] * v[653] + v[343] * v[656];
	v[534] = -(v[353] * v[408]) + v[355] * v[409] + v[343] * v[414] - v[345] * v[415];
	v[533] = -(v[353] * v[405]) + v[355] * v[406] + v[343] * v[411] - v[345] * v[412];
	v[362] = -(v[345] * v[353]) + v[343] * v[355];
	v[356] = v[384] * v[4205];
	v[4588] = -(v[356] * v[650]);
	v[4585] = -(v[356] * v[651]);
	v[4582] = -(v[356] * v[652]);
	v[2007] = 2e0*v[4580] + 2e0*v[4581] + 2e0*v[4582] - v[415] * v[676] - v[414] * v[682] - v[413] * v[688] + v[370] * v[757]
		+ v[369] * v[766] + v[368] * v[775];
	v[2006] = 2e0*v[4583] + 2e0*v[4584] + 2e0*v[4585] - v[415] * v[674] - v[414] * v[680] - v[413] * v[686] + v[370] * v[760]
		+ v[369] * v[769] + v[368] * v[778];
	v[2005] = 2e0*v[4586] + 2e0*v[4587] + 2e0*v[4588] - v[415] * v[672] - v[414] * v[678] - v[413] * v[684] + v[370] * v[763]
		+ v[369] * v[772] + v[368] * v[781];
	v[2001] = -(v[353] * v[628]) - v[355] * v[631] - v[356] * v[634] - v[379] * v[655] - v[381] * v[658] - v[383] * v[664]
		- v[412] * v[676] - v[411] * v[682] - v[410] * v[688] + v[370] * v[755] + v[369] * v[764] + v[368] * v[773];
	v[2000] = -(v[353] * v[627]) - v[355] * v[630] - v[356] * v[633] - v[379] * v[654] - v[381] * v[657] - v[383] * v[663]
		- v[412] * v[674] - v[411] * v[680] - v[410] * v[686] + v[370] * v[758] + v[369] * v[767] + v[368] * v[776];
	v[1999] = -(v[353] * v[626]) - v[355] * v[629] - v[356] * v[632] - v[379] * v[653] - v[381] * v[656] - v[383] * v[662]
		- v[412] * v[672] - v[411] * v[678] - v[410] * v[684] + v[370] * v[761] + v[369] * v[770] + v[368] * v[779];
	v[1612] = v[1538] * v[368] + v[1542] * v[369] + v[1546] * v[370] - v[353] * v[396] - v[355] * v[397] - v[356] * v[398]
		- v[413] * v[4572] - v[382] * v[4589] + v[380] * v[4590];
	v[1610] = v[1536] * v[368] + v[1540] * v[369] + v[1544] * v[370] - v[353] * v[386] - v[355] * v[388] - v[356] * v[390]
		- v[384] * v[410] - v[382] * v[411] - v[380] * v[412] - v[383] * v[413] - v[381] * v[414] - v[379] * v[415];
	v[1609] = v[1548] * v[368] + v[1550] * v[369] + v[1552] * v[370] - v[353] * v[385] - v[355] * v[387] - v[356] * v[389]
		- v[410] * v[4566] - v[381] * v[4591] + v[379] * v[4592];
	v[1561] = -(v[1552] * v[343]) + v[1548] * v[346] + v[1522] * v[353] - v[1518] * v[356] - v[404] * v[4592]
		- v[406] * v[4593];
	v[1560] = -(v[1546] * v[343]) + v[1538] * v[346] + v[1512] * v[353] - v[1504] * v[356] - v[407] * v[4590]
		- v[409] * v[4594];
	v[1559] = -(v[1544] * v[343]) + v[1536] * v[346] + v[1510] * v[353] - v[1502] * v[356] - v[409] * v[410] + v[407] * v[412]
		- v[406] * v[413] + v[404] * v[415];
	v[1558] = v[1552] * v[345] - v[1550] * v[346] - v[1522] * v[355] + v[1520] * v[356] - v[404] * v[4591] + v[405] * v[4593];
	v[1557] = v[1546] * v[345] - v[1542] * v[346] - v[1512] * v[355] + v[1508] * v[356] - v[407] * v[4589] + v[408] * v[4594];
	v[1556] = v[1544] * v[345] - v[1540] * v[346] - v[1510] * v[355] + v[1506] * v[356] + v[408] * v[410] - v[407] * v[411]
		+ v[405] * v[413] - v[404] * v[414];
	v[670] = v[356] * v[640] - v[355] * v[643] - v[346] * v[658] + v[345] * v[664];
	v[669] = v[356] * v[639] - v[355] * v[642] - v[346] * v[657] + v[345] * v[663];
	v[668] = v[356] * v[638] - v[355] * v[641] - v[346] * v[656] + v[345] * v[662];
	v[667] = -(v[356] * v[637]) + v[353] * v[643] + v[346] * v[655] - v[343] * v[664];
	v[666] = -(v[356] * v[636]) + v[353] * v[642] + v[346] * v[654] - v[343] * v[663];
	v[665] = -(v[356] * v[635]) + v[353] * v[641] + v[346] * v[653] - v[343] * v[662];
	v[538] = -(v[355] * v[407]) + v[356] * v[408] + v[345] * v[413] - v[346] * v[414];
	v[537] = -(v[355] * v[404]) + v[356] * v[405] + v[345] * v[410] - v[346] * v[411];
	v[536] = v[353] * v[407] - v[356] * v[409] - v[343] * v[413] + v[346] * v[415];
	v[535] = v[353] * v[404] - v[356] * v[406] - v[343] * v[410] + v[346] * v[412];
	v[360] = v[346] * v[353] - v[343] * v[356];
	v[357] = -(v[346] * v[355]) + v[345] * v[356];
	v[783] = 2e0*(v[362] * v[534] + v[360] * v[536] + v[357] * v[538]);
	v[782] = 2e0*(v[362] * v[533] + v[360] * v[535] + v[357] * v[537]);
	v[540] = (v[357] * v[357]) + (v[360] * v[360]) + (v[362] * v[362]);
	v[4596] = 1e0 / Power(v[540], 2);
	v[4215] = 1e0 / sqrt(v[540]);
	v[4597] = -(v[4215] * v[540]);
	v[4595] = v[4215] / 2e0;
	v[4217] = -v[4215] / (4e0*v[540]);
	v[4598] = -(v[4217] * v[540]);
	v[3587] = v[4595] * v[783];
	v[3584] = v[4595] * v[782];
	v[1570] = v[4596] * (v[4597] * (v[1558] * v[357] + v[1561] * v[360] + v[1555] * v[362] + (v[533] * v[533]) +
		(v[535] * v[535]) + (v[537] * v[537])) + v[3584] * v[782] + v[4598] * (v[782] * v[782]));
	v[1569] = v[4596] * (v[4597] * (v[1557] * v[357] + v[1560] * v[360] + v[1554] * v[362] + (v[534] * v[534]) +
		(v[536] * v[536]) + (v[538] * v[538])) + v[3587] * v[783] + v[4598] * (v[783] * v[783]));
	v[1567] = v[4596] * (v[4597] * (v[1556] * v[357] + v[1559] * v[360] + v[1553] * v[362] + v[533] * v[534] + v[535] * v[536]
		+ v[537] * v[538]) + v[782] * (v[3587] + v[4598] * v[783]));
	v[785] = -(v[3587] / v[540]);
	v[4600] = 2e0*v[785];
	v[784] = -(v[3584] / v[540]);
	v[4599] = 2e0*v[784];
	v[1589] = v[1570] * v[357] + v[1558] * v[4215] + v[4599] * v[537];
	v[3597] = -((*radius)*v[1589]) - v[385];
	v[1583] = v[1569] * v[357] + v[1557] * v[4215] + v[4600] * v[538];
	v[3596] = -((*radius)*v[1583]) - v[396];
	v[1581] = v[1567] * v[357] + v[1556] * v[4215] + v[538] * v[784] + v[537] * v[785];
	v[3595] = -((*radius)*v[1581]) - v[386];
	v[1580] = v[1570] * v[360] + v[1561] * v[4215] + v[4599] * v[535];
	v[3594] = -((*radius)*v[1580]) - v[387];
	v[1578] = v[1569] * v[360] + v[1560] * v[4215] + v[4600] * v[536];
	v[3593] = -((*radius)*v[1578]) - v[397];
	v[1576] = v[1567] * v[360] + v[1559] * v[4215] + v[536] * v[784] + v[535] * v[785];
	v[3592] = -((*radius)*v[1576]) - v[388];
	v[1575] = v[1570] * v[362] + v[1555] * v[4215] + v[4599] * v[533];
	v[3591] = -((*radius)*v[1575]) - v[389];
	v[1591] = v[1589] * v[316] + v[1580] * v[318] + v[1575] * v[319];
	v[1590] = v[1589] * v[331] + v[1580] * v[334] + v[1575] * v[336];
	v[1573] = v[1569] * v[362] + v[1554] * v[4215] + v[4600] * v[534];
	v[3590] = -((*radius)*v[1573]) - v[398];
	v[1587] = v[1583] * v[316] + v[1578] * v[318] + v[1573] * v[319];
	v[1585] = v[1583] * v[331] + v[1578] * v[334] + v[1573] * v[336];
	v[1571] = v[1567] * v[362] + v[1553] * v[4215] + v[534] * v[784] + v[533] * v[785];
	v[3589] = -((*radius)*v[1571]) - v[390];
	v[1586] = v[1581] * v[316] + v[1576] * v[318] + v[1571] * v[319];
	v[1584] = v[1581] * v[331] + v[1576] * v[334] + v[1571] * v[336];
	v[1103] = v[4215] * (v[1022] * v[343] - v[1034] * v[345] - v[4601] * v[640] + v[4602] * v[658] - v[353] * v[986]
		+ v[355] * v[998]);
	v[3369] = -((*radius)*v[1103]) - cp[1] * v[924] - cp[0] * v[925] - v[304] * v[953];
	v[1101] = v[4215] * (-(v[1010] * v[343]) + v[1034] * v[346] + v[4601] * v[643] - v[4602] * v[664] + v[353] * v[974]
		- v[356] * v[998]);
	v[3376] = -((*radius)*v[1101]) - cp[1] * v[928] - cp[0] * v[930] - v[304] * v[959];
	v[1099] = v[4215] * (v[1010] * v[345] - v[1022] * v[346] - v[4615] * v[658] + v[4611] * v[664] - v[355] * v[974]
		+ v[356] * v[986]);
	v[3383] = -((*radius)*v[1099]) - cp[1] * v[945] - cp[0] * v[947] - v[304] * v[965];
	v[2878] = v[1099] * v[316] + v[1101] * v[318] + v[1103] * v[319];
	v[2875] = v[1099] * v[331] + v[1101] * v[334] + v[1103] * v[336];
	v[1097] = v[4215] * (v[1020] * v[343] - v[1032] * v[345] - v[640] * v[654] - v[639] * v[655] + v[637] * v[657]
		+ v[636] * v[658] - v[353] * v[984] + v[355] * v[996]);
	v[3346] = -((*radius)*v[1097]) - cp[1] * v[889] - cp[0] * v[892] - v[304] * v[952];
	v[1095] = v[4215] * (v[1026] * v[343] - v[1038] * v[345] + v[1002] * v[355] - v[4603] * v[639] + v[4604] * v[657]
		- v[353] * v[990]);
	v[3344] = -((*radius)*v[1095]) - cp[1] * v[888] - cp[0] * v[891] - v[304] * v[955];
	v[1093] = v[4215] * (-(v[1008] * v[343]) + v[1032] * v[346] + v[643] * v[654] + v[642] * v[655] - v[637] * v[663]
		- v[636] * v[664] + v[353] * v[972] - v[356] * v[996]);
	v[3354] = -((*radius)*v[1093]) - cp[1] * v[915] - cp[0] * v[918] - v[304] * v[958];
	v[1091] = v[4215] * (-(v[1014] * v[343]) + v[1038] * v[346] - v[1002] * v[356] + v[4603] * v[642] - v[4604] * v[663]
		+ v[353] * v[978]);
	v[3352] = -((*radius)*v[1091]) - cp[1] * v[914] - cp[0] * v[917] - v[304] * v[961];
	v[1089] = v[4215] * (v[1008] * v[345] - v[1020] * v[346] - v[643] * v[657] - v[642] * v[658] + v[640] * v[663]
		+ v[639] * v[664] - v[355] * v[972] + v[356] * v[984]);
	v[3362] = -((*radius)*v[1089]) - cp[1] * v[944] - cp[0] * v[946] - v[304] * v[964];
	v[2877] = v[1089] * v[316] + v[1093] * v[318] + v[1097] * v[319];
	v[2874] = v[1089] * v[331] + v[1093] * v[334] + v[1097] * v[336];
	v[1087] = v[4215] * (v[1014] * v[345] - v[1026] * v[346] - v[4614] * v[657] + v[4610] * v[663] - v[355] * v[978]
		+ v[356] * v[990]);
	v[3360] = -((*radius)*v[1087]) - cp[1] * v[949] - cp[0] * v[950] - v[304] * v[967];
	v[2882] = v[1087] * v[316] + v[1091] * v[318] + v[1095] * v[319];
	v[2880] = v[1087] * v[331] + v[1091] * v[334] + v[1095] * v[336];
	v[1085] = v[4215] * (v[1018] * v[343] - v[1030] * v[345] - v[640] * v[653] - v[638] * v[655] + v[637] * v[656]
		+ v[635] * v[658] - v[353] * v[982] + v[355] * v[994]);
	v[3320] = -((*radius)*v[1085]) - cp[1] * v[896] - cp[0] * v[898] - v[304] * v[951];
	v[1083] = v[4215] * (v[1024] * v[343] - v[1036] * v[345] + v[1000] * v[355] - v[639] * v[653] - v[638] * v[654]
		+ v[636] * v[656] + v[635] * v[657] - v[353] * v[988]);
	v[3318] = -((*radius)*v[1083]) - cp[1] * v[887] - cp[0] * v[890] - v[304] * v[954];
	v[1081] = v[4215] * (v[1028] * v[343] - v[1040] * v[345] + v[1004] * v[355] - v[4605] * v[638] + v[4606] * v[656]
		- v[353] * v[992]);
	v[3316] = -((*radius)*v[1081]) - cp[1] * v[895] - cp[0] * v[897] - v[304] * v[956];
	v[1079] = v[4215] * (-(v[1006] * v[343]) + v[1030] * v[346] + v[643] * v[653] + v[641] * v[655] - v[637] * v[662]
		- v[635] * v[664] + v[353] * v[970] - v[356] * v[994]);
	v[3329] = -((*radius)*v[1079]) - cp[1] * v[927] - cp[0] * v[929] - v[304] * v[957];
	v[1077] = v[4215] * (-(v[1012] * v[343]) + v[1036] * v[346] - v[1000] * v[356] + v[642] * v[653] + v[641] * v[654]
		- v[636] * v[662] - v[635] * v[663] + v[353] * v[976]);
	v[3327] = -((*radius)*v[1077]) - cp[1] * v[913] - cp[0] * v[916] - v[304] * v[960];
	v[1075] = v[4215] * (-(v[1016] * v[343]) + v[1040] * v[346] - v[1004] * v[356] + v[4605] * v[641] - v[4606] * v[662]
		+ v[353] * v[980]);
	v[3325] = -((*radius)*v[1075]) - cp[1] * v[932] - cp[0] * v[933] - v[304] * v[962];
	v[1073] = v[4215] * (v[1006] * v[345] - v[1018] * v[346] - v[643] * v[656] - v[641] * v[658] + v[640] * v[662]
		+ v[638] * v[664] - v[355] * v[970] + v[356] * v[982]);
	v[3338] = -((*radius)*v[1073]) - cp[1] * v[937] - cp[0] * v[940] - v[304] * v[963];
	v[2876] = v[1073] * v[316] + v[1079] * v[318] + v[1085] * v[319];
	v[2873] = v[1073] * v[331] + v[1079] * v[334] + v[1085] * v[336];
	v[1071] = v[4215] * (v[1012] * v[345] - v[1024] * v[346] - v[642] * v[656] - v[641] * v[657] + v[639] * v[662]
		+ v[638] * v[663] - v[355] * v[976] + v[356] * v[988]);
	v[3336] = -((*radius)*v[1071]) - cp[1] * v[936] - cp[0] * v[939] - v[304] * v[966];
	v[2881] = v[1071] * v[316] + v[1077] * v[318] + v[1083] * v[319];
	v[2879] = v[1071] * v[331] + v[1077] * v[334] + v[1083] * v[336];
	v[1069] = v[4215] * (v[1016] * v[345] - v[1028] * v[346] - v[4613] * v[656] + v[4609] * v[662] - v[355] * v[980]
		+ v[356] * v[992]);
	v[3334] = -((*radius)*v[1069]) - cp[1] * v[935] - cp[0] * v[938] - v[304] * v[968];
	v[2884] = v[1069] * v[316] + v[1075] * v[318] + v[1081] * v[319];
	v[2883] = v[1069] * v[331] + v[1075] * v[334] + v[1081] * v[336];
	v[827] = v[4215] * (v[414] * v[637] - v[415] * v[640] - v[408] * v[655] + v[409] * v[658] - v[353] * v[734] + v[355] * v[746]
		+ v[343] * v[766] - v[345] * v[775]) + v[661] * v[785];
	v[3207] = -v[652] - (*radius)*v[827];
	v[826] = v[4215] * (v[411] * v[637] - v[412] * v[640] - v[405] * v[655] + v[406] * v[658] - v[353] * v[732] + v[355] * v[744]
		+ v[343] * v[764] - v[345] * v[773]) + v[661] * v[784];
	v[3210] = -v[634] - (*radius)*v[826];
	v[825] = v[4215] * (-(v[413] * v[637]) + v[415] * v[643] + v[407] * v[655] - v[409] * v[664] + v[353] * v[722]
		- v[356] * v[746] - v[343] * v[757] + v[346] * v[775]) + v[667] * v[785];
	v[3213] = -v[649] - (*radius)*v[825];
	v[824] = v[4215] * (-(v[410] * v[637]) + v[412] * v[643] + v[404] * v[655] - v[406] * v[664] + v[353] * v[720]
		- v[356] * v[744] - v[343] * v[755] + v[346] * v[773]) + v[667] * v[784];
	v[3216] = -v[631] - (*radius)*v[824];
	v[823] = v[4215] * (v[413] * v[640] - v[414] * v[643] - v[407] * v[658] + v[408] * v[664] - v[355] * v[722] + v[356] * v[734]
		+ v[345] * v[757] - v[346] * v[766]) + v[670] * v[785];
	v[3219] = -v[646] - (*radius)*v[823];
	v[1959] = v[316] * v[823] + v[318] * v[825] + v[319] * v[827];
	v[1956] = v[331] * v[823] + v[334] * v[825] + v[336] * v[827];
	v[822] = v[4215] * (v[410] * v[640] - v[411] * v[643] - v[404] * v[658] + v[405] * v[664] - v[355] * v[720] + v[356] * v[732]
		+ v[345] * v[755] - v[346] * v[764]) + v[670] * v[784];
	v[3222] = -v[628] - (*radius)*v[822];
	v[1965] = v[316] * v[822] + v[318] * v[824] + v[319] * v[826];
	v[1962] = v[331] * v[822] + v[334] * v[824] + v[336] * v[826];
	v[821] = v[4215] * (v[414] * v[636] - v[415] * v[639] - v[408] * v[654] + v[409] * v[657] - v[353] * v[738] + v[355] * v[750]
		+ v[343] * v[769] - v[345] * v[778]) + v[660] * v[785];
	v[3206] = -v[651] - (*radius)*v[821];
	v[820] = v[4215] * (v[411] * v[636] - v[412] * v[639] - v[405] * v[654] + v[406] * v[657] - v[353] * v[736] + v[355] * v[748]
		+ v[343] * v[767] - v[345] * v[776]) + v[660] * v[784];
	v[3209] = -v[633] - (*radius)*v[820];
	v[819] = v[4215] * (-(v[413] * v[636]) + v[415] * v[642] + v[407] * v[654] - v[409] * v[663] + v[353] * v[726]
		- v[356] * v[750] - v[343] * v[760] + v[346] * v[778]) + v[666] * v[785];
	v[3212] = -v[648] - (*radius)*v[819];
	v[818] = v[4215] * (-(v[410] * v[636]) + v[412] * v[642] + v[404] * v[654] - v[406] * v[663] + v[353] * v[724]
		- v[356] * v[748] - v[343] * v[758] + v[346] * v[776]) + v[666] * v[784];
	v[3215] = -v[630] - (*radius)*v[818];
	v[817] = v[4215] * (v[413] * v[639] - v[414] * v[642] - v[407] * v[657] + v[408] * v[663] - v[355] * v[726] + v[356] * v[738]
		+ v[345] * v[760] - v[346] * v[769]) + v[669] * v[785];
	v[3218] = -v[645] - (*radius)*v[817];
	v[1958] = v[316] * v[817] + v[318] * v[819] + v[319] * v[821];
	v[1955] = v[331] * v[817] + v[334] * v[819] + v[336] * v[821];
	v[816] = v[4215] * (v[410] * v[639] - v[411] * v[642] - v[404] * v[657] + v[405] * v[663] - v[355] * v[724] + v[356] * v[736]
		+ v[345] * v[758] - v[346] * v[767]) + v[669] * v[784];
	v[3221] = -v[627] - (*radius)*v[816];
	v[1964] = v[316] * v[816] + v[318] * v[818] + v[319] * v[820];
	v[1961] = v[331] * v[816] + v[334] * v[818] + v[336] * v[820];
	v[815] = v[4215] * (v[414] * v[635] - v[415] * v[638] - v[408] * v[653] + v[409] * v[656] - v[353] * v[742] + v[355] * v[754]
		+ v[343] * v[772] - v[345] * v[781]) + v[659] * v[785];
	v[3205] = -v[650] - (*radius)*v[815];
	v[814] = v[4215] * (v[411] * v[635] - v[412] * v[638] - v[405] * v[653] + v[406] * v[656] - v[353] * v[740] + v[355] * v[752]
		+ v[343] * v[770] - v[345] * v[779]) + v[659] * v[784];
	v[3208] = -v[632] - (*radius)*v[814];
	v[813] = v[4215] * (-(v[413] * v[635]) + v[415] * v[641] + v[407] * v[653] - v[409] * v[662] + v[353] * v[730]
		- v[356] * v[754] - v[343] * v[763] + v[346] * v[781]) + v[665] * v[785];
	v[3211] = -v[647] - (*radius)*v[813];
	v[812] = v[4215] * (-(v[410] * v[635]) + v[412] * v[641] + v[404] * v[653] - v[406] * v[662] + v[353] * v[728]
		- v[356] * v[752] - v[343] * v[761] + v[346] * v[779]) + v[665] * v[784];
	v[3214] = -v[629] - (*radius)*v[812];
	v[811] = v[4215] * (v[413] * v[638] - v[414] * v[641] - v[407] * v[656] + v[408] * v[662] - v[355] * v[730] + v[356] * v[742]
		+ v[345] * v[763] - v[346] * v[772]) + v[668] * v[785];
	v[3217] = -v[644] - (*radius)*v[811];
	v[1957] = v[316] * v[811] + v[318] * v[813] + v[319] * v[815];
	v[1954] = v[331] * v[811] + v[334] * v[813] + v[336] * v[815];
	v[810] = v[4215] * (v[410] * v[638] - v[411] * v[641] - v[404] * v[656] + v[405] * v[662] - v[355] * v[728] + v[356] * v[740]
		+ v[345] * v[761] - v[346] * v[770]) + v[668] * v[784];
	v[3220] = -v[626] - (*radius)*v[810];
	v[1963] = v[316] * v[810] + v[318] * v[812] + v[319] * v[814];
	v[1960] = v[331] * v[810] + v[334] * v[812] + v[336] * v[814];
	v[687] = v[4215] * v[661];
	v[4757] = -2e0*v[687];
	v[1067] = -((*radius)*v[687]) - v[688];
	v[685] = v[4215] * v[660];
	v[4756] = -2e0*v[685];
	v[1066] = -((*radius)*v[685]) - v[686];
	v[683] = v[4215] * v[659];
	v[4755] = -2e0*v[683];
	v[1065] = -((*radius)*v[683]) - v[684];
	v[681] = v[4215] * v[667];
	v[4760] = -2e0*v[681];
	v[1064] = -((*radius)*v[681]) - v[682];
	v[679] = v[4215] * v[666];
	v[4759] = -2e0*v[679];
	v[1063] = -((*radius)*v[679]) - v[680];
	v[677] = v[4215] * v[665];
	v[4758] = -2e0*v[677];
	v[1062] = -((*radius)*v[677]) - v[678];
	v[675] = v[4215] * v[670];
	v[4763] = -2e0*v[675];
	v[2905] = v[1241] * v[675] + v[1260] * v[681] + v[1280] * v[687];
	v[2901] = v[1246] * v[675] + v[1264] * v[681] + v[1285] * v[687];
	v[2896] = v[1250] * v[675] + v[1269] * v[681] + v[1289] * v[687];
	v[1232] = v[331] * v[675] + v[334] * v[681] + v[336] * v[687];
	v[1214] = v[316] * v[675] + v[318] * v[681] + v[319] * v[687];
	v[1061] = -((*radius)*v[675]) - v[676];
	v[673] = v[4215] * v[669];
	v[4762] = -2e0*v[673];
	v[2904] = v[1241] * v[673] + v[1260] * v[679] + v[1280] * v[685];
	v[2900] = v[1246] * v[673] + v[1264] * v[679] + v[1285] * v[685];
	v[2895] = v[1250] * v[673] + v[1269] * v[679] + v[1289] * v[685];
	v[1231] = v[331] * v[673] + v[334] * v[679] + v[336] * v[685];
	v[1213] = v[316] * v[673] + v[318] * v[679] + v[319] * v[685];
	v[1060] = -((*radius)*v[673]) - v[674];
	v[671] = v[4215] * v[668];
	v[4761] = -2e0*v[671];
	v[2903] = v[1241] * v[671] + v[1260] * v[677] + v[1280] * v[683];
	v[2899] = v[1246] * v[671] + v[1264] * v[677] + v[1285] * v[683];
	v[2894] = v[1250] * v[671] + v[1269] * v[677] + v[1289] * v[683];
	v[1230] = v[331] * v[671] + v[334] * v[677] + v[336] * v[683];
	v[1212] = v[316] * v[671] + v[318] * v[677] + v[319] * v[683];
	v[1059] = -((*radius)*v[671]) - v[672];
	v[549] = v[4215] * v[534] + v[362] * v[785];
	v[4764] = 2e0*v[549];
	v[809] = -v[384] - (*radius)*v[549];
	v[548] = v[4215] * v[533] + v[362] * v[784];
	v[4766] = 2e0*v[548];
	v[808] = -v[383] - (*radius)*v[548];
	v[547] = v[4215] * v[536] + v[360] * v[785];
	v[4618] = 2e0*v[547];
	v[807] = -v[382] - (*radius)*v[547];
	v[546] = v[4215] * v[535] + v[360] * v[784];
	v[4616] = 2e0*v[546];
	v[806] = -v[381] - (*radius)*v[546];
	v[545] = v[4215] * v[538] + v[357] * v[785];
	v[4765] = -2e0*v[545];
	v[2676] = v[1250] * v[545] + v[1269] * v[547] + v[1289] * v[549];
	v[2675] = v[1246] * v[545] + v[1264] * v[547] + v[1285] * v[549];
	v[2674] = v[1241] * v[545] + v[1260] * v[547] + v[1280] * v[549];
	v[1373] = v[331] * v[545] + v[334] * v[547] + v[336] * v[549];
	v[1361] = v[316] * v[545] + v[318] * v[547] + v[319] * v[549];
	v[805] = -v[380] - (*radius)*v[545];
	v[544] = v[4215] * v[537] + v[357] * v[784];
	v[4767] = -2e0*v[544];
	v[2682] = v[1250] * v[544] + v[1269] * v[546] + v[1289] * v[548];
	v[2681] = v[1246] * v[544] + v[1264] * v[546] + v[1285] * v[548];
	v[2680] = v[1241] * v[544] + v[1260] * v[546] + v[1280] * v[548];
	v[1372] = v[331] * v[544] + v[334] * v[546] + v[336] * v[548];
	v[1360] = v[316] * v[544] + v[318] * v[546] + v[319] * v[548];
	v[804] = -v[379] - (*radius)*v[544];
	v[358] = v[357] * v[4215];
	v[4607] = -2e0*v[358];
	v[3010] = v[4607] * v[675];
	v[3009] = v[4607] * v[673];
	v[3008] = v[4607] * v[671];
	v[2710] = v[4607] * v[545];
	v[2686] = v[4607] * v[544];
	v[361] = v[360] * v[4215];
	v[4608] = -2e0*v[361];
	v[3019] = v[4608] * v[681];
	v[3018] = v[4608] * v[679];
	v[3017] = v[4608] * v[677];
	v[3013] = -(v[361] * v[675]) - v[358] * v[681];
	v[3012] = -(v[361] * v[673]) - v[358] * v[679];
	v[3011] = -(v[361] * v[671]) - v[358] * v[677];
	v[2890] = -(v[1075] * v[343]) + v[1069] * v[345] - v[1004] * v[361] + v[4609] * v[671] - v[4606] * v[677]
		+ v[358] * v[992];
	v[2889] = -(v[1091] * v[343]) + v[1087] * v[345] - v[1002] * v[361] + v[4610] * v[673] - v[4604] * v[679]
		+ v[358] * v[990];
	v[2888] = -(v[1077] * v[343]) + v[1071] * v[345] - v[1000] * v[361] + v[639] * v[671] + v[638] * v[673] - v[636] * v[677]
		- v[635] * v[679] + v[358] * v[988];
	v[2887] = -(v[1101] * v[343]) + v[1099] * v[345] + v[4611] * v[675] - v[4602] * v[681] + v[358] * v[986] - v[361] * v[998];
	v[2886] = -(v[1093] * v[343]) + v[1089] * v[345] + v[640] * v[673] + v[639] * v[675] - v[637] * v[679] - v[636] * v[681]
		+ v[358] * v[984] - v[361] * v[996];
	v[2885] = -(v[1079] * v[343]) + v[1073] * v[345] + v[640] * v[671] + v[638] * v[675] - v[637] * v[677] - v[635] * v[681]
		+ v[358] * v[982] - v[361] * v[994];
	v[2719] = v[4608] * v[547];
	v[2711] = -(v[361] * v[545]) - v[358] * v[547];
	v[2695] = v[4608] * v[546];
	v[2687] = -(v[361] * v[544]) - v[358] * v[546];
	v[1971] = -(v[546] * v[637]) + v[544] * v[640] + v[405] * v[675] - v[406] * v[681] + v[358] * v[732] - v[361] * v[744]
		+ v[345] * v[822] - v[343] * v[824];
	v[1970] = -(v[546] * v[636]) + v[544] * v[639] + v[405] * v[673] - v[406] * v[679] + v[358] * v[736] - v[361] * v[748]
		+ v[345] * v[816] - v[343] * v[818];
	v[1969] = -(v[546] * v[635]) + v[544] * v[638] + v[405] * v[671] - v[406] * v[677] + v[358] * v[740] - v[361] * v[752]
		+ v[345] * v[810] - v[343] * v[812];
	v[1968] = -(v[547] * v[637]) + v[545] * v[640] + v[408] * v[675] - v[409] * v[681] + v[358] * v[734] - v[361] * v[746]
		+ v[345] * v[823] - v[343] * v[825];
	v[1967] = -(v[547] * v[636]) + v[545] * v[639] + v[408] * v[673] - v[409] * v[679] + v[358] * v[738] - v[361] * v[750]
		+ v[345] * v[817] - v[343] * v[819];
	v[1966] = -(v[547] * v[635]) + v[545] * v[638] + v[408] * v[671] - v[409] * v[677] + v[358] * v[742] - v[361] * v[754]
		+ v[345] * v[811] - v[343] * v[813];
	v[1594] = -(v[1580] * v[343]) + v[1589] * v[345] + v[1520] * v[358] - v[1518] * v[361] - v[406] * v[4616]
		+ v[4617] * v[544];
	v[1593] = -(v[1578] * v[343]) + v[1583] * v[345] + v[1508] * v[358] - v[1504] * v[361] - v[409] * v[4618]
		+ v[4619] * v[545];
	v[1592] = -(v[1576] * v[343]) + v[1581] * v[345] + v[1506] * v[358] - v[1502] * v[361] + v[408] * v[544] + v[405] * v[545]
		- v[409] * v[546] - v[406] * v[547];
	v[1355] = v[358] * v[408] - v[361] * v[409] + v[345] * v[545] - v[343] * v[547];
	v[1354] = v[358] * v[405] - v[361] * v[406] + v[345] * v[544] - v[343] * v[546];
	v[1205] = -(v[361] * v[637]) + v[358] * v[640] + v[345] * v[675] - v[343] * v[681];
	v[1204] = -(v[361] * v[636]) + v[358] * v[639] + v[345] * v[673] - v[343] * v[679];
	v[1203] = -(v[361] * v[635]) + v[358] * v[638] + v[345] * v[671] - v[343] * v[677];
	v[363] = v[362] * v[4215];
	v[4612] = -2e0*v[363];
	v[3025] = v[4612] * v[687];
	v[3024] = v[4612] * v[685];
	v[3023] = v[4612] * v[683];
	v[3022] = -(v[363] * v[681]) - v[361] * v[687];
	v[3021] = -(v[363] * v[679]) - v[361] * v[685];
	v[3020] = -(v[363] * v[677]) - v[361] * v[683];
	v[3016] = -(v[363] * v[675]) - v[358] * v[687];
	v[3015] = -(v[363] * v[673]) - v[358] * v[685];
	v[3014] = -(v[363] * v[671]) - v[358] * v[683];
	v[2927] = -(v[1081] * v[345]) + v[1075] * v[346] + v[4613] * v[677] - v[4609] * v[683] + v[361] * v[980] - v[363] * v[992];
	v[2922] = -(v[1095] * v[345]) + v[1091] * v[346] + v[4614] * v[679] - v[4610] * v[685] + v[361] * v[978] - v[363] * v[990];
	v[2921] = -(v[1083] * v[345]) + v[1077] * v[346] + v[642] * v[677] + v[641] * v[679] - v[639] * v[683] - v[638] * v[685]
		+ v[361] * v[976] - v[363] * v[988];
	v[2914] = -(v[1103] * v[345]) + v[1101] * v[346] + v[4615] * v[681] - v[4611] * v[687] + v[361] * v[974] - v[363] * v[986];
	v[2913] = -(v[1097] * v[345]) + v[1093] * v[346] + v[643] * v[679] + v[642] * v[681] - v[640] * v[685] - v[639] * v[687]
		+ v[361] * v[972] - v[363] * v[984];
	v[2912] = -(v[1085] * v[345]) + v[1079] * v[346] + v[643] * v[677] + v[641] * v[681] - v[640] * v[683] - v[638] * v[687]
		+ v[361] * v[970] - v[363] * v[982];
	v[2911] = v[1081] * v[343] - v[1069] * v[346] + v[1004] * v[363] - v[4613] * v[671] + v[4606] * v[683] - v[358] * v[980];
	v[2929] = v[2927] * v[316] + v[2911] * v[318] + v[2890] * v[319];
	v[2928] = v[2927] * v[331] + v[2911] * v[334] + v[2890] * v[336];
	v[2910] = v[1095] * v[343] - v[1087] * v[346] + v[1002] * v[363] - v[4614] * v[673] + v[4604] * v[685] - v[358] * v[978];
	v[2926] = v[2922] * v[316] + v[2910] * v[318] + v[2889] * v[319];
	v[2924] = v[2922] * v[331] + v[2910] * v[334] + v[2889] * v[336];
	v[2909] = v[1083] * v[343] - v[1071] * v[346] + v[1000] * v[363] - v[642] * v[671] - v[641] * v[673] + v[636] * v[683]
		+ v[635] * v[685] - v[358] * v[976];
	v[2925] = v[2921] * v[316] + v[2909] * v[318] + v[2888] * v[319];
	v[2923] = v[2921] * v[331] + v[2909] * v[334] + v[2888] * v[336];
	v[2908] = v[1103] * v[343] - v[1099] * v[346] - v[4615] * v[675] + v[4602] * v[687] - v[358] * v[974] + v[363] * v[998];
	v[2920] = v[2914] * v[316] + v[2908] * v[318] + v[2887] * v[319];
	v[2917] = v[2914] * v[331] + v[2908] * v[334] + v[2887] * v[336];
	v[2907] = v[1097] * v[343] - v[1089] * v[346] - v[643] * v[673] - v[642] * v[675] + v[637] * v[685] + v[636] * v[687]
		- v[358] * v[972] + v[363] * v[996];
	v[2919] = v[2913] * v[316] + v[2907] * v[318] + v[2886] * v[319];
	v[2916] = v[2913] * v[331] + v[2907] * v[334] + v[2886] * v[336];
	v[2906] = v[1085] * v[343] - v[1073] * v[346] - v[643] * v[671] - v[641] * v[675] + v[637] * v[683] + v[635] * v[687]
		- v[358] * v[970] + v[363] * v[994];
	v[2918] = v[2912] * v[316] + v[2906] * v[318] + v[2885] * v[319];
	v[2915] = v[2912] * v[331] + v[2906] * v[334] + v[2885] * v[336];
	v[2902] = v[2860] * v[358] + v[2834] * v[361] + v[2806] * v[363];
	v[3727] = v[2806] - v[2902] * v[363];
	v[3721] = v[2834] - v[2902] * v[361];
	v[3715] = v[2860] - v[2902] * v[358];
	v[2898] = v[2856] * v[358] + v[2829] * v[361] + v[2801] * v[363];
	v[3743] = v[2801] - v[2898] * v[363];
	v[3738] = v[2829] - v[2898] * v[361];
	v[3733] = v[2856] - v[2898] * v[358];
	v[2897] = v[2851] * v[358] + v[2825] * v[361] + v[2796] * v[363];
	v[3728] = v[2796] - v[2897] * v[363];
	v[3722] = v[2825] - v[2897] * v[361];
	v[3716] = v[2851] - v[2897] * v[358];
	v[2893] = v[2847] * v[358] + v[2819] * v[361] + v[2792] * v[363];
	v[3756] = v[2792] - v[2893] * v[363];
	v[3752] = v[2819] - v[2893] * v[361];
	v[3748] = v[2847] - v[2893] * v[358];
	v[2892] = v[2842] * v[358] + v[2814] * v[361] + v[2788] * v[363];
	v[3744] = v[2788] - v[2892] * v[363];
	v[3739] = v[2814] - v[2892] * v[361];
	v[3734] = v[2842] - v[2892] * v[358];
	v[2891] = v[2838] * v[358] + v[2810] * v[361] + v[2784] * v[363];
	v[3729] = v[2784] - v[2891] * v[363];
	v[3723] = v[2810] - v[2891] * v[361];
	v[3717] = v[2838] - v[2891] * v[358];
	v[2727] = v[4612] * v[549];
	v[2720] = -(v[363] * v[547]) - v[361] * v[549];
	v[2712] = -(v[363] * v[545]) - v[358] * v[549];
	v[2703] = v[4612] * v[548];
	v[2696] = -(v[363] * v[546]) - v[361] * v[548];
	v[2688] = -(v[363] * v[544]) - v[358] * v[548];
	v[1989] = -(v[548] * v[640]) + v[546] * v[643] + v[404] * v[681] - v[405] * v[687] + v[361] * v[720] - v[363] * v[732]
		+ v[346] * v[824] - v[345] * v[826];
	v[1988] = -(v[548] * v[639]) + v[546] * v[642] + v[404] * v[679] - v[405] * v[685] + v[361] * v[724] - v[363] * v[736]
		+ v[346] * v[818] - v[345] * v[820];
	v[1987] = -(v[548] * v[638]) + v[546] * v[641] + v[404] * v[677] - v[405] * v[683] + v[361] * v[728] - v[363] * v[740]
		+ v[346] * v[812] - v[345] * v[814];
	v[1980] = -(v[549] * v[640]) + v[547] * v[643] + v[407] * v[681] - v[408] * v[687] + v[361] * v[722] - v[363] * v[734]
		+ v[346] * v[825] - v[345] * v[827];
	v[1979] = -(v[549] * v[639]) + v[547] * v[642] + v[407] * v[679] - v[408] * v[685] + v[361] * v[726] - v[363] * v[738]
		+ v[346] * v[819] - v[345] * v[821];
	v[1978] = -(v[549] * v[638]) + v[547] * v[641] + v[407] * v[677] - v[408] * v[683] + v[361] * v[730] - v[363] * v[742]
		+ v[346] * v[813] - v[345] * v[815];
	v[1977] = v[548] * v[637] - v[544] * v[643] - v[404] * v[675] + v[406] * v[687] - v[358] * v[720] + v[363] * v[744]
		- v[346] * v[822] + v[343] * v[826];
	v[1995] = v[1989] * v[316] + v[1977] * v[318] + v[1971] * v[319];
	v[1992] = v[1989] * v[331] + v[1977] * v[334] + v[1971] * v[336];
	v[1976] = v[548] * v[636] - v[544] * v[642] - v[404] * v[673] + v[406] * v[685] - v[358] * v[724] + v[363] * v[748]
		- v[346] * v[816] + v[343] * v[820];
	v[1994] = v[1988] * v[316] + v[1976] * v[318] + v[1970] * v[319];
	v[1991] = v[1988] * v[331] + v[1976] * v[334] + v[1970] * v[336];
	v[1975] = v[548] * v[635] - v[544] * v[641] - v[404] * v[671] + v[406] * v[683] - v[358] * v[728] + v[363] * v[752]
		- v[346] * v[810] + v[343] * v[814];
	v[1993] = v[1987] * v[316] + v[1975] * v[318] + v[1969] * v[319];
	v[1990] = v[1987] * v[331] + v[1975] * v[334] + v[1969] * v[336];
	v[1974] = v[549] * v[637] - v[545] * v[643] - v[407] * v[675] + v[409] * v[687] - v[358] * v[722] + v[363] * v[746]
		- v[346] * v[823] + v[343] * v[827];
	v[1986] = v[1980] * v[316] + v[1974] * v[318] + v[1968] * v[319];
	v[1983] = v[1980] * v[331] + v[1974] * v[334] + v[1968] * v[336];
	v[1973] = v[549] * v[636] - v[545] * v[642] - v[407] * v[673] + v[409] * v[685] - v[358] * v[726] + v[363] * v[750]
		- v[346] * v[817] + v[343] * v[821];
	v[1985] = v[1979] * v[316] + v[1973] * v[318] + v[1967] * v[319];
	v[1982] = v[1979] * v[331] + v[1973] * v[334] + v[1967] * v[336];
	v[1972] = v[549] * v[635] - v[545] * v[641] - v[407] * v[671] + v[409] * v[683] - v[358] * v[730] + v[363] * v[754]
		- v[346] * v[811] + v[343] * v[815];
	v[1984] = v[1978] * v[316] + v[1972] * v[318] + v[1966] * v[319];
	v[1981] = v[1978] * v[331] + v[1972] * v[334] + v[1966] * v[336];
	v[1604] = -(v[1575] * v[345]) + v[1580] * v[346] + v[1522] * v[361] - v[1520] * v[363] + v[404] * v[4616]
		- v[4617] * v[548];
	v[1599] = -(v[1573] * v[345]) + v[1578] * v[346] + v[1512] * v[361] - v[1508] * v[363] + v[407] * v[4618]
		- v[4619] * v[549];
	v[1598] = -(v[1571] * v[345]) + v[1576] * v[346] + v[1510] * v[361] - v[1506] * v[363] + v[407] * v[546] + v[404] * v[547]
		- v[408] * v[548] - v[405] * v[549];
	v[1597] = v[1575] * v[343] - v[1589] * v[346] - v[1522] * v[358] + v[1518] * v[363] + v[406] * v[4766] + v[404] * v[4767];
	v[1606] = v[1604] * v[316] + v[1597] * v[318] + v[1594] * v[319];
	v[1605] = v[1604] * v[331] + v[1597] * v[334] + v[1594] * v[336];
	v[1596] = v[1573] * v[343] - v[1583] * v[346] - v[1512] * v[358] + v[1504] * v[363] + v[409] * v[4764] + v[407] * v[4765];
	v[1603] = v[1599] * v[316] + v[1596] * v[318] + v[1593] * v[319];
	v[1601] = v[1599] * v[331] + v[1596] * v[334] + v[1593] * v[336];
	v[1595] = v[1571] * v[343] - v[1581] * v[346] - v[1510] * v[358] + v[1502] * v[363] - v[407] * v[544] - v[404] * v[545]
		+ v[409] * v[548] + v[406] * v[549];
	v[1602] = v[1598] * v[316] + v[1595] * v[318] + v[1592] * v[319];
	v[1600] = v[1598] * v[331] + v[1595] * v[334] + v[1592] * v[336];
	v[1353] = -(v[358] * v[407]) + v[363] * v[409] - v[346] * v[545] + v[343] * v[549];
	v[1352] = -(v[358] * v[404]) + v[363] * v[406] - v[346] * v[544] + v[343] * v[548];
	v[1351] = v[361] * v[407] - v[363] * v[408] + v[346] * v[547] - v[345] * v[549];
	v[1371] = v[1351] * v[331] + v[1353] * v[334] + v[1355] * v[336];
	v[1359] = v[1351] * v[316] + v[1353] * v[318] + v[1355] * v[319];
	v[1350] = v[361] * v[404] - v[363] * v[405] + v[346] * v[546] - v[345] * v[548];
	v[1370] = v[1350] * v[331] + v[1352] * v[334] + v[1354] * v[336];
	v[1358] = v[1350] * v[316] + v[1352] * v[318] + v[1354] * v[319];
	v[1295] = v[1250] * v[358] + v[1269] * v[361] + v[1289] * v[363];
	v[3759] = -(v[2896] * v[363]) - v[1295] * v[687];
	v[3758] = -(v[2895] * v[363]) - v[1295] * v[685];
	v[3757] = -(v[2894] * v[363]) - v[1295] * v[683];
	v[3755] = -(v[2896] * v[361]) - v[1295] * v[681];
	v[3754] = -(v[2895] * v[361]) - v[1295] * v[679];
	v[3753] = -(v[2894] * v[361]) - v[1295] * v[677];
	v[3751] = -(v[2896] * v[358]) - v[1295] * v[675];
	v[3750] = -(v[2895] * v[358]) - v[1295] * v[673];
	v[3749] = -(v[2894] * v[358]) - v[1295] * v[671];
	v[2730] = -(v[2676] * v[363]) - v[1295] * v[549];
	v[2723] = -(v[2676] * v[361]) - v[1295] * v[547];
	v[2715] = -(v[2676] * v[358]) - v[1295] * v[545];
	v[2706] = -(v[2682] * v[363]) - v[1295] * v[548];
	v[2699] = -(v[2682] * v[361]) - v[1295] * v[546];
	v[2691] = -(v[2682] * v[358]) - v[1295] * v[544];
	v[1294] = v[1246] * v[358] + v[1264] * v[361] + v[1285] * v[363];
	v[3747] = -(v[2901] * v[363]) - v[1294] * v[687];
	v[3746] = -(v[2900] * v[363]) - v[1294] * v[685];
	v[3745] = -(v[2899] * v[363]) - v[1294] * v[683];
	v[3742] = -(v[2901] * v[361]) - v[1294] * v[681];
	v[3741] = -(v[2900] * v[361]) - v[1294] * v[679];
	v[3740] = -(v[2899] * v[361]) - v[1294] * v[677];
	v[3737] = -(v[2901] * v[358]) - v[1294] * v[675];
	v[3736] = -(v[2900] * v[358]) - v[1294] * v[673];
	v[3735] = -(v[2899] * v[358]) - v[1294] * v[671];
	v[2729] = -(v[2675] * v[363]) - v[1294] * v[549];
	v[2722] = -(v[2675] * v[361]) - v[1294] * v[547];
	v[2714] = -(v[2675] * v[358]) - v[1294] * v[545];
	v[2705] = -(v[2681] * v[363]) - v[1294] * v[548];
	v[2698] = -(v[2681] * v[361]) - v[1294] * v[546];
	v[2690] = -(v[2681] * v[358]) - v[1294] * v[544];
	v[1293] = v[1241] * v[358] + v[1260] * v[361] + v[1280] * v[363];
	v[3732] = -(v[2905] * v[363]) - v[1293] * v[687];
	v[3731] = -(v[2904] * v[363]) - v[1293] * v[685];
	v[3730] = -(v[2903] * v[363]) - v[1293] * v[683];
	v[3726] = -(v[2905] * v[361]) - v[1293] * v[681];
	v[3725] = -(v[2904] * v[361]) - v[1293] * v[679];
	v[3724] = -(v[2903] * v[361]) - v[1293] * v[677];
	v[3720] = -(v[2905] * v[358]) - v[1293] * v[675];
	v[3719] = -(v[2904] * v[358]) - v[1293] * v[673];
	v[3718] = -(v[2903] * v[358]) - v[1293] * v[671];
	v[2728] = -(v[2674] * v[363]) - v[1293] * v[549];
	v[2721] = -(v[2674] * v[361]) - v[1293] * v[547];
	v[2713] = -(v[2674] * v[358]) - v[1293] * v[545];
	v[2704] = -(v[2680] * v[363]) - v[1293] * v[548];
	v[2697] = -(v[2680] * v[361]) - v[1293] * v[546];
	v[2689] = -(v[2680] * v[358]) - v[1293] * v[544];
	v[1202] = v[363] * v[637] - v[358] * v[643] - v[346] * v[675] + v[343] * v[687];
	v[1201] = v[363] * v[636] - v[358] * v[642] - v[346] * v[673] + v[343] * v[685];
	v[1200] = v[363] * v[635] - v[358] * v[641] - v[346] * v[671] + v[343] * v[683];
	v[1199] = -(v[363] * v[640]) + v[361] * v[643] + v[346] * v[681] - v[345] * v[687];
	v[1229] = v[1199] * v[331] + v[1202] * v[334] + v[1205] * v[336];
	v[1211] = v[1199] * v[316] + v[1202] * v[318] + v[1205] * v[319];
	v[1198] = -(v[363] * v[639]) + v[361] * v[642] + v[346] * v[679] - v[345] * v[685];
	v[1228] = v[1198] * v[331] + v[1201] * v[334] + v[1204] * v[336];
	v[1210] = v[1198] * v[316] + v[1201] * v[318] + v[1204] * v[319];
	v[1197] = -(v[363] * v[638]) + v[361] * v[641] + v[346] * v[677] - v[345] * v[683];
	v[1227] = v[1197] * v[331] + v[1200] * v[334] + v[1203] * v[336];
	v[1209] = v[1197] * v[316] + v[1200] * v[318] + v[1203] * v[319];
	v[364] = -((*radius)*v[358]) + v[368];
	v[365] = -((*radius)*v[361]) + v[369];
	v[366] = -((*radius)*v[363]) + v[370];
	v[416] = -(v[343] * v[379]) - v[345] * v[381] - v[346] * v[383] + v[370] * v[404] + v[369] * v[405] + v[368] * v[406];
	v[4734] = v[353] * v[416];
	v[4731] = v[415] * v[416];
	v[4730] = v[355] * v[416];
	v[4727] = v[414] * v[416];
	v[4724] = v[356] * v[416];
	v[4721] = v[413] * v[416];
	v[417] = -(v[353] * v[379]) - v[355] * v[381] - v[356] * v[383] + v[370] * v[410] + v[369] * v[411] + v[368] * v[412];
	v[4733] = v[343] * v[417];
	v[4732] = -(v[409] * v[417]);
	v[4729] = v[345] * v[417];
	v[4728] = -(v[408] * v[417]);
	v[4723] = v[346] * v[417];
	v[4722] = -(v[407] * v[417]);
	v[418] = -(v[343] * v[380]) - v[345] * v[382] - v[346] * v[384] + v[370] * v[407] + v[369] * v[408] + v[368] * v[409];
	v[4750] = v[412] * v[418];
	v[4748] = v[411] * v[418];
	v[4746] = v[410] * v[418];
	v[419] = -(v[353] * v[380]) - v[355] * v[382] - v[356] * v[384] + v[370] * v[413] + v[369] * v[414] + v[368] * v[415];
	v[4754] = v[356] * v[418] - v[346] * v[419];
	v[4753] = v[355] * v[418] - v[345] * v[419];
	v[4752] = v[353] * v[418] - v[343] * v[419];
	v[4751] = v[406] * v[419];
	v[4749] = v[405] * v[419];
	v[4747] = v[404] * v[419];
	v[1784] = -(v[417] * v[418]) + v[416] * v[419];
	v[4743] = v[1999] / v[1784];
	v[4742] = -(v[1996] / v[1784]);
	v[4741] = v[2000] / v[1784];
	v[4740] = -(v[1997] / v[1784]);
	v[4737] = v[2001] / v[1784];
	v[4736] = -(v[1998] / v[1784]);
	v[4726] = v[1608] / v[1784];
	v[4725] = -(v[1610] / v[1784]);
	v[4720] = -(v[2005] / v[1784]);
	v[4719] = v[2002] / v[1784];
	v[4716] = -(v[2006] / v[1784]);
	v[4715] = v[2003] / v[1784];
	v[4713] = -(v[2007] / v[1784]);
	v[4712] = v[2004] / v[1784];
	v[4680] = -(v[419] / v[1784]);
	v[4679] = v[418] / v[1784];
	v[4638] = v[417] / v[1784];
	v[4637] = -(v[416] / v[1784]);
	v[1786] = 1e0 / Power(v[1784], 2);
	v[2481] = -(v[1786] * (v[2007] * v[416] - v[2004] * v[417] - v[2001] * v[418] + v[1998] * v[419]));
	v[4690] = -(v[2481] * v[419]) + v[4713];
	v[4689] = v[2481] * v[418] + v[4712];
	v[4652] = v[2481] * v[417] + v[4737];
	v[4651] = -(v[2481] * v[416]) + v[4736];
	v[4636] = v[1784] * v[2481];
	v[2479] = -(v[1786] * (v[2006] * v[416] - v[2003] * v[417] - v[2000] * v[418] + v[1997] * v[419]));
	v[4684] = v[2479] * v[418];
	v[4692] = v[4684] + v[4715];
	v[4683] = v[2479] * v[419];
	v[4693] = -v[4683] + v[4716];
	v[4644] = v[2479] * v[416];
	v[4654] = -v[4644] + v[4740];
	v[4643] = v[2479] * v[417];
	v[4655] = v[4643] + v[4741];
	v[2477] = -(v[1786] * (v[2005] * v[416] - v[2002] * v[417] - v[1999] * v[418] + v[1996] * v[419]));
	v[4695] = v[2477] * v[418];
	v[4694] = v[2477] * v[419];
	v[4686] = -v[4694] + v[4720];
	v[4685] = v[4695] + v[4719];
	v[4657] = v[2477] * v[416];
	v[4656] = v[2477] * v[417];
	v[4646] = v[4656] + v[4743];
	v[4645] = -v[4657] + v[4742];
	v[2475] = -(v[1786] * (v[4721] + v[4722] - v[4746] + v[4747]));
	v[4664] = v[1784] * v[2475];
	v[4704] = v[413] + v[419] * v[4664];
	v[4703] = v[407] + v[418] * v[4664];
	v[4671] = v[410] + v[417] * v[4664];
	v[4670] = v[404] + v[416] * v[4664];
	v[2473] = -(v[1786] * (v[4727] + v[4728] - v[4748] + v[4749]));
	v[4665] = v[1784] * v[2473];
	v[4706] = v[414] + v[419] * v[4665];
	v[4705] = v[408] + v[418] * v[4665];
	v[4673] = v[411] + v[417] * v[4665];
	v[4672] = v[405] + v[416] * v[4665];
	v[2471] = -(v[1786] * (v[4731] + v[4732] - v[4750] + v[4751]));
	v[4666] = v[1784] * v[2471];
	v[4708] = v[415] + v[419] * v[4666];
	v[4707] = v[409] + v[418] * v[4666];
	v[4675] = v[412] + v[417] * v[4666];
	v[4674] = v[406] + v[416] * v[4666];
	v[1787] = -(v[1786] * (v[1612] * v[416] - v[1611] * v[417] - v[1610] * v[418] + v[1608] * v[419]));
	v[4735] = v[1784] * v[1787];
	v[4745] = v[418] * v[4735];
	v[4744] = v[419] * v[4735];
	v[4739] = v[1612] + v[4744];
	v[4738] = v[1611] + v[4745];
	v[1785] = -(v[1786] * (v[1610] * v[416] - v[1608] * v[417] - v[1609] * v[418] + v[1607] * v[419]));
	v[4714] = v[1784] * v[1785];
	v[4718] = v[1609] + v[417] * v[4714];
	v[4717] = v[1607] + v[416] * v[4714];
	v[424] = (v[368] * v[380] + v[369] * v[382] + v[370] * v[384])*v[4620];
	v[2040] = v[2016] * v[384] + v[424] * v[652] - v[4205] * v[688];
	v[2085] = v[2040] * v[349] - v[304] * v[664];
	v[2039] = v[2014] * v[384] + v[424] * v[651] - v[4205] * v[686];
	v[2084] = v[2039] * v[349] - v[304] * v[663];
	v[2038] = v[2012] * v[384] + v[424] * v[650] - v[4205] * v[684];
	v[2083] = v[2038] * v[349] - v[304] * v[662];
	v[2036] = v[2016] * v[382] + v[424] * v[649] - v[4205] * v[682];
	v[2077] = v[2036] * v[349] - v[304] * v[658];
	v[2035] = v[2014] * v[382] + v[424] * v[648] - v[4205] * v[680];
	v[2076] = v[2035] * v[349] - v[304] * v[657];
	v[2034] = v[2012] * v[382] + v[424] * v[647] - v[4205] * v[678];
	v[2075] = v[2034] * v[349] - v[304] * v[656];
	v[2031] = v[2016] * v[380] + v[424] * v[646] - v[4205] * v[676];
	v[2067] = v[2031] * v[349] - v[304] * v[655];
	v[2030] = v[2014] * v[380] + v[424] * v[645] - v[4205] * v[674];
	v[2066] = v[2030] * v[349] - v[304] * v[654];
	v[2029] = v[2012] * v[380] + v[424] * v[644] - v[4205] * v[672];
	v[2065] = v[2029] * v[349] - v[304] * v[653];
	v[429] = (v[368] * v[379] + v[369] * v[381] + v[370] * v[383])*v[4621];
	v[2055] = v[2025] * v[383] + v[429] * v[634] - v[4195] * v[688];
	v[2081] = v[2055] * v[339] - v[304] * v[643];
	v[2054] = v[2023] * v[383] + v[429] * v[633] - v[4195] * v[686];
	v[2080] = v[2054] * v[339] - v[304] * v[642];
	v[2053] = v[2021] * v[383] + v[429] * v[632] - v[4195] * v[684];
	v[2079] = v[2053] * v[339] - v[304] * v[641];
	v[2051] = v[2025] * v[381] + v[429] * v[631] - v[4195] * v[682];
	v[2072] = v[2051] * v[339] - v[304] * v[640];
	v[2050] = v[2023] * v[381] + v[429] * v[630] - v[4195] * v[680];
	v[2071] = v[2050] * v[339] - v[304] * v[639];
	v[2049] = v[2021] * v[381] + v[429] * v[629] - v[4195] * v[678];
	v[2070] = v[2049] * v[339] - v[304] * v[638];
	v[2046] = v[2025] * v[379] + v[429] * v[628] - v[4195] * v[676];
	v[2061] = v[2046] * v[339] - v[304] * v[637];
	v[2045] = v[2023] * v[379] + v[429] * v[627] - v[4195] * v[674];
	v[2060] = v[2045] * v[339] - v[304] * v[636];
	v[2044] = v[2021] * v[379] + v[429] * v[626] - v[4195] * v[672];
	v[2059] = v[2044] * v[339] - v[304] * v[635];
	v[423] = v[368] * v[4205] + v[380] * v[424];
	v[425] = v[369] * v[4205] + v[382] * v[424];
	v[426] = v[370] * v[4205] + v[384] * v[424];
	v[428] = v[368] * v[4195] + v[379] * v[429];
	v[430] = v[369] * v[4195] + v[381] * v[429];
	v[431] = v[370] * v[4195] + v[383] * v[429];
	v[432] = -(v[304] * v[343]) + v[339] * v[428];
	v[433] = -(v[304] * v[353]) + v[349] * v[423];
	v[434] = -(v[304] * v[345]) + v[339] * v[430];
	v[2146] = v[2055] + v[2072] * v[284] - v[2061] * v[285] + v[434] * v[610] - v[432] * v[613] - cp[0] * v[643];
	v[2145] = v[2054] + v[2071] * v[284] - v[2060] * v[285] + v[434] * v[609] - v[432] * v[612] - cp[0] * v[642];
	v[2144] = v[2053] + v[2070] * v[284] - v[2059] * v[285] + v[434] * v[608] - v[432] * v[611] - cp[0] * v[641];
	v[2116] = -(v[2072] * v[281]) + v[2061] * v[282] - v[434] * v[601] + v[432] * v[604] - cp[1] * v[643];
	v[2211] = e2i[0] * v[2116] + e1i[0] * v[2146];
	v[2187] = e2i[1] * v[2116] + e1i[1] * v[2146];
	v[2175] = e2i[2] * v[2116] + e1i[2] * v[2146];
	v[2115] = -(v[2071] * v[281]) + v[2060] * v[282] - v[434] * v[600] + v[432] * v[603] - cp[1] * v[642];
	v[2210] = e2i[0] * v[2115] + e1i[0] * v[2145];
	v[2186] = e2i[1] * v[2115] + e1i[1] * v[2145];
	v[2174] = e2i[2] * v[2115] + e1i[2] * v[2145];
	v[2114] = -(v[2070] * v[281]) + v[2059] * v[282] - v[434] * v[599] + v[432] * v[602] - cp[1] * v[641];
	v[2209] = e2i[0] * v[2114] + e1i[0] * v[2144];
	v[2185] = e2i[1] * v[2114] + e1i[1] * v[2144];
	v[2173] = e2i[2] * v[2114] + e1i[2] * v[2144];
	v[435] = -(v[304] * v[355]) + v[349] * v[425];
	v[2153] = v[2077] * v[284] - v[2067] * v[285] + v[435] * v[610] - v[433] * v[613] - cp[0] * v[664];
	v[2152] = v[2076] * v[284] - v[2066] * v[285] + v[435] * v[609] - v[433] * v[612] - cp[0] * v[663];
	v[2151] = v[2075] * v[284] - v[2065] * v[285] + v[435] * v[608] - v[433] * v[611] - cp[0] * v[662];
	v[2123] = v[2040] - v[2077] * v[281] + v[2067] * v[282] - v[435] * v[601] + v[433] * v[604] - cp[1] * v[664];
	v[2223] = e2i[0] * v[2123] + e1i[0] * v[2153];
	v[2199] = e2i[1] * v[2123] + e1i[1] * v[2153];
	v[2181] = e2i[2] * v[2123] + e1i[2] * v[2153];
	v[2122] = v[2039] - v[2076] * v[281] + v[2066] * v[282] - v[435] * v[600] + v[433] * v[603] - cp[1] * v[663];
	v[2222] = e2i[0] * v[2122] + e1i[0] * v[2152];
	v[2198] = e2i[1] * v[2122] + e1i[1] * v[2152];
	v[2180] = e2i[2] * v[2122] + e1i[2] * v[2152];
	v[2121] = v[2038] - v[2075] * v[281] + v[2065] * v[282] - v[435] * v[599] + v[433] * v[602] - cp[1] * v[662];
	v[2221] = e2i[0] * v[2121] + e1i[0] * v[2151];
	v[2197] = e2i[1] * v[2121] + e1i[1] * v[2151];
	v[2179] = e2i[2] * v[2121] + e1i[2] * v[2151];
	v[436] = -(v[304] * v[346]) + v[339] * v[431];
	v[2161] = v[2051] - v[2081] * v[284] + v[2061] * v[286] - v[436] * v[610] + v[432] * v[616] - cp[0] * v[640];
	v[2160] = v[2050] - v[2080] * v[284] + v[2060] * v[286] - v[436] * v[609] + v[432] * v[615] - cp[0] * v[639];
	v[2159] = v[2049] - v[2079] * v[284] + v[2059] * v[286] - v[436] * v[608] + v[432] * v[614] - cp[0] * v[638];
	v[2131] = v[2081] * v[281] - v[2061] * v[283] + v[436] * v[601] - v[432] * v[607] - cp[1] * v[640];
	v[2271] = e2i[0] * v[2131] + e1i[0] * v[2161];
	v[2259] = e2i[1] * v[2131] + e1i[1] * v[2161];
	v[2235] = e2i[2] * v[2131] + e1i[2] * v[2161];
	v[2130] = v[2080] * v[281] - v[2060] * v[283] + v[436] * v[600] - v[432] * v[606] - cp[1] * v[639];
	v[2270] = e2i[0] * v[2130] + e1i[0] * v[2160];
	v[2258] = e2i[1] * v[2130] + e1i[1] * v[2160];
	v[2234] = e2i[2] * v[2130] + e1i[2] * v[2160];
	v[2129] = v[2079] * v[281] - v[2059] * v[283] + v[436] * v[599] - v[432] * v[605] - cp[1] * v[638];
	v[2269] = e2i[0] * v[2129] + e1i[0] * v[2159];
	v[2257] = e2i[1] * v[2129] + e1i[1] * v[2159];
	v[2233] = e2i[2] * v[2129] + e1i[2] * v[2159];
	v[2103] = v[2046] + v[2081] * v[285] - v[2072] * v[286] + v[436] * v[613] - v[434] * v[616] - cp[0] * v[637];
	v[2102] = v[2045] + v[2080] * v[285] - v[2071] * v[286] + v[436] * v[612] - v[434] * v[615] - cp[0] * v[636];
	v[2101] = v[2044] + v[2079] * v[285] - v[2070] * v[286] + v[436] * v[611] - v[434] * v[614] - cp[0] * v[635];
	v[2091] = -(v[2081] * v[282]) + v[2072] * v[283] - v[436] * v[604] + v[434] * v[607] - cp[1] * v[637];
	v[2343] = e2i[0] * v[2091] + e1i[0] * v[2103];
	v[2319] = e2i[1] * v[2091] + e1i[1] * v[2103];
	v[2295] = e2i[2] * v[2091] + e1i[2] * v[2103];
	v[2090] = -(v[2080] * v[282]) + v[2071] * v[283] - v[436] * v[603] + v[434] * v[606] - cp[1] * v[636];
	v[2342] = e2i[0] * v[2090] + e1i[0] * v[2102];
	v[2318] = e2i[1] * v[2090] + e1i[1] * v[2102];
	v[2294] = e2i[2] * v[2090] + e1i[2] * v[2102];
	v[2089] = -(v[2079] * v[282]) + v[2070] * v[283] - v[436] * v[602] + v[434] * v[605] - cp[1] * v[635];
	v[2317] = e2i[1] * v[2089] + e1i[1] * v[2101];
	v[2293] = e2i[2] * v[2089] + e1i[2] * v[2101];
	v[437] = -(v[304] * v[356]) + v[349] * v[426];
	v[2169] = -(v[2085] * v[284]) + v[2067] * v[286] - v[437] * v[610] + v[433] * v[616] - cp[0] * v[658];
	v[2168] = -(v[2084] * v[284]) + v[2066] * v[286] - v[437] * v[609] + v[433] * v[615] - cp[0] * v[657];
	v[2167] = -(v[2083] * v[284]) + v[2065] * v[286] - v[437] * v[608] + v[433] * v[614] - cp[0] * v[656];
	v[2139] = v[2036] + v[2085] * v[281] - v[2067] * v[283] + v[437] * v[601] - v[433] * v[607] - cp[1] * v[658];
	v[2283] = e2i[0] * v[2139] + e1i[0] * v[2169];
	v[2265] = e2i[1] * v[2139] + e1i[1] * v[2169];
	v[2247] = e2i[2] * v[2139] + e1i[2] * v[2169];
	v[2138] = v[2035] + v[2084] * v[281] - v[2066] * v[283] + v[437] * v[600] - v[433] * v[606] - cp[1] * v[657];
	v[2282] = e2i[0] * v[2138] + e1i[0] * v[2168];
	v[2264] = e2i[1] * v[2138] + e1i[1] * v[2168];
	v[2246] = e2i[2] * v[2138] + e1i[2] * v[2168];
	v[2137] = v[2034] + v[2083] * v[281] - v[2065] * v[283] + v[437] * v[599] - v[433] * v[605] - cp[1] * v[656];
	v[2281] = e2i[0] * v[2137] + e1i[0] * v[2167];
	v[2263] = e2i[1] * v[2137] + e1i[1] * v[2167];
	v[2245] = e2i[2] * v[2137] + e1i[2] * v[2167];
	v[2109] = v[2085] * v[285] - v[2077] * v[286] + v[437] * v[613] - v[435] * v[616] - cp[0] * v[655];
	v[2108] = v[2084] * v[285] - v[2076] * v[286] + v[437] * v[612] - v[435] * v[615] - cp[0] * v[654];
	v[2107] = v[2083] * v[285] - v[2075] * v[286] + v[437] * v[611] - v[435] * v[614] - cp[0] * v[653];
	v[2097] = v[2031] - v[2085] * v[282] + v[2077] * v[283] - v[437] * v[604] + v[435] * v[607] - cp[1] * v[655];
	v[2349] = e2i[0] * v[2097] + e1i[0] * v[2109];
	v[2331] = e2i[1] * v[2097] + e1i[1] * v[2109];
	v[2307] = e2i[2] * v[2097] + e1i[2] * v[2109];
	v[2096] = v[2030] - v[2084] * v[282] + v[2076] * v[283] - v[437] * v[603] + v[435] * v[606] - cp[1] * v[654];
	v[2348] = e2i[0] * v[2096] + e1i[0] * v[2108];
	v[2330] = e2i[1] * v[2096] + e1i[1] * v[2108];
	v[2306] = e2i[2] * v[2096] + e1i[2] * v[2108];
	v[2095] = v[2029] - v[2083] * v[282] + v[2075] * v[283] - v[437] * v[602] + v[435] * v[605] - cp[1] * v[653];
	v[2329] = e2i[1] * v[2095] + e1i[1] * v[2107];
	v[2305] = e2i[2] * v[2095] + e1i[2] * v[2107];
	v[438] = -(cp[1] * v[343]) + v[283] * v[434] - v[282] * v[436];
	v[439] = -(cp[1] * v[353]) + v[423] + v[283] * v[435] - v[282] * v[437];
	v[440] = -(cp[0] * v[343]) + v[428] - v[286] * v[434] + v[285] * v[436];
	v[441] = -(cp[0] * v[353]) - v[286] * v[435] + v[285] * v[437];
	v[442] = -(cp[1] * v[346]) + v[282] * v[432] - v[281] * v[434];
	v[443] = -(cp[1] * v[356]) + v[426] + v[282] * v[433] - v[281] * v[435];
	v[444] = -(cp[1] * v[345]) - v[283] * v[432] + v[281] * v[436];
	v[445] = -(cp[1] * v[355]) + v[425] - v[283] * v[433] + v[281] * v[437];
	v[446] = -(cp[0] * v[346]) + v[431] - v[285] * v[432] + v[284] * v[434];
	v[447] = -(cp[0] * v[356]) - v[285] * v[433] + v[284] * v[435];
	v[448] = -(cp[0] * v[345]) + v[430] + v[286] * v[432] - v[284] * v[436];
	v[449] = -(cp[0] * v[355]) + v[286] * v[433] - v[284] * v[437];
	v[450] = e2i[2] * v[442] + e1i[2] * v[446];
	v[451] = e2i[2] * v[443] + e1i[2] * v[447];
	v[452] = e2i[1] * v[442] + e1i[1] * v[446];
	v[2193] = v[2187] * v[243] + v[452] * v[565];
	v[2192] = v[2186] * v[243] + v[452] * v[564];
	v[497] = v[243] * v[452];
	v[453] = e2i[1] * v[443] + e1i[1] * v[447];
	v[2205] = v[2199] * v[243] + v[453] * v[565];
	v[2204] = v[2198] * v[243] + v[453] * v[564];
	v[510] = v[243] * v[453];
	v[454] = e2i[0] * v[442] + e1i[0] * v[446];
	v[2217] = v[2211] * v[243] + v[454] * v[565];
	v[2216] = v[2210] * v[243] + v[454] * v[564];
	v[500] = v[243] * v[454];
	v[455] = e2i[0] * v[443] + e1i[0] * v[447];
	v[2229] = v[2223] * v[243] + v[455] * v[565];
	v[2228] = v[2222] * v[243] + v[455] * v[564];
	v[513] = v[243] * v[455];
	v[456] = e2i[2] * v[444] + e1i[2] * v[448];
	v[4634] = v[452] - v[456];
	v[4623] = v[452] + v[456];
	v[2241] = v[2235] * v[243] + v[456] * v[565];
	v[2394] = v[2193] + v[2241];
	v[2240] = v[2234] * v[243] + v[456] * v[564];
	v[498] = v[243] * v[456];
	v[457] = e2i[2] * v[445] + e1i[2] * v[449];
	v[4630] = v[453] - v[457];
	v[4625] = v[453] + v[457];
	v[2253] = v[2247] * v[243] + v[457] * v[565];
	v[2400] = v[2205] + v[2253];
	v[2252] = v[2246] * v[243] + v[457] * v[564];
	v[511] = v[243] * v[457];
	v[458] = e2i[1] * v[444] + e1i[1] * v[448];
	v[4635] = v[450] + v[458];
	v[459] = e2i[1] * v[445] + e1i[1] * v[449];
	v[4631] = v[451] + v[459];
	v[460] = e2i[0] * v[444] + e1i[0] * v[448];
	v[2277] = v[2271] * v[243] + v[460] * v[565];
	v[505] = v[243] * v[460];
	v[461] = e2i[0] * v[445] + e1i[0] * v[449];
	v[2289] = v[2283] * v[243] + v[461] * v[565];
	v[518] = v[243] * v[461];
	v[462] = e2i[2] * v[438] + e1i[2] * v[440];
	v[4633] = v[454] + v[462];
	v[2301] = v[2295] * v[243] + v[462] * v[565];
	v[2406] = v[2217] + v[2301];
	v[2300] = v[2294] * v[243] + v[462] * v[564];
	v[501] = v[243] * v[462];
	v[463] = e2i[2] * v[439] + e1i[2] * v[441];
	v[4629] = v[455] + v[463];
	v[2313] = v[2307] * v[243] + v[463] * v[565];
	v[2412] = v[2229] + v[2313];
	v[2312] = v[2306] * v[243] + v[463] * v[564];
	v[514] = v[243] * v[463];
	v[464] = e2i[1] * v[438] + e1i[1] * v[440];
	v[4622] = v[460] + v[464];
	v[2325] = v[2319] * v[243] + v[464] * v[565];
	v[2418] = v[2277] + v[2325];
	v[2417] = (v[2270] + v[2318])*v[243] + v[4622] * v[564];
	v[506] = v[243] * v[464];
	v[465] = e2i[1] * v[439] + e1i[1] * v[441];
	v[4624] = v[461] + v[465];
	v[2337] = v[2331] * v[243] + v[465] * v[565];
	v[2424] = v[2289] + v[2337];
	v[2423] = (v[2282] + v[2330])*v[243] + v[4624] * v[564];
	v[519] = v[243] * v[465];
	v[466] = e2i[0] * v[438] + e1i[0] * v[440];
	v[467] = e2i[0] * v[439] + e1i[0] * v[441];
	v[2388] = v[2352] * v[450] + v[2175] * v[883];
	v[2387] = v[2351] * v[450] + v[2174] * v[883];
	v[2382] = v[2352] * v[458] + v[2259] * v[883];
	v[2376] = v[2352] * v[466] + v[2343] * v[883];
	v[2370] = v[2352] * v[451] + v[2181] * v[883];
	v[2369] = v[2351] * v[451] + v[2180] * v[883];
	v[2364] = v[2352] * v[459] + v[2265] * v[883];
	v[2358] = v[2352] * v[467] + v[2349] * v[883];
	v[521] = v[467] * v[883];
	v[520] = v[459] * v[883];
	v[515] = v[451] * v[883];
	v[508] = v[466] * v[883];
	v[507] = v[458] * v[883];
	v[502] = v[450] * v[883];
	v[474] = v[497] + v[498];
	v[475] = v[510] + v[511];
	v[480] = v[500] + v[501];
	v[481] = v[513] + v[514];
	v[483] = v[505] + v[506];
	v[484] = v[518] + v[519];
	v[487] = v[4492] * v[450] + v[4502] * v[458] + v[4501] * v[466] + v[452] * v[471] + v[454] * v[472] + v[456] * v[473]
		+ v[460] * v[478] + v[462] * v[479] + v[464] * v[482];
	v[488] = v[4492] * v[451] + v[4502] * v[459] + v[4501] * v[467] + v[453] * v[471] + v[455] * v[472] + v[457] * v[473]
		+ v[461] * v[478] + v[463] * v[479] + v[465] * v[482];
	v[2451] = v[4493] * (v[2175] * v[4492] + v[2343] * v[4501] + v[2259] * v[4502] + v[460] + v[4486] * v[4623]
		+ v[4482] * v[4633] - v[464] + d[11] * (-v[458] - v[466]) + v[2187] * v[471] + v[2211] * v[472] + v[2235] * v[473]
		+ v[2271] * v[478] + v[2295] * v[479] + v[2319] * v[482]) + v[2427] * v[487];
	v[4627] = v[2376] + v[2451];
	v[2449] = v[4493] * (v[2174] * v[4492] + v[2342] * v[4501] + v[2258] * v[4502] - v[454] + v[462] + v[4482] * v[4622]
		+ v[4489] * v[4623] + d[10] * (-v[450] - v[466]) + v[2186] * v[471] + v[2210] * v[472] + v[2234] * v[473] + v[2270] * v[478]
		+ v[2294] * v[479] + v[2318] * v[482]) + v[2426] * v[487];
	v[4632] = v[2387] + v[2449];
	v[2439] = v[4493] * (v[2181] * v[4492] + v[2349] * v[4501] + v[2265] * v[4502] + v[461] + v[4486] * v[4625]
		+ v[4482] * v[4629] - v[465] + d[11] * (-v[459] - v[467]) + v[2199] * v[471] + v[2223] * v[472] + v[2247] * v[473]
		+ v[2283] * v[478] + v[2307] * v[479] + v[2331] * v[482]) + v[2427] * v[488];
	v[4626] = v[2358] + v[2439];
	v[2437] = v[4493] * (v[2180] * v[4492] + v[2348] * v[4501] + v[2264] * v[4502] - v[455] + v[4482] * v[4624]
		+ v[4489] * v[4625] + v[463] + d[10] * (-v[451] - v[467]) + v[2198] * v[471] + v[2222] * v[472] + v[2246] * v[473]
		+ v[2282] * v[478] + v[2306] * v[479] + v[2330] * v[482]) + v[2426] * v[488];
	v[4628] = v[2369] + v[2437];
	v[516] = v[4493] * v[488];
	v[2468] = v[516] + v[520] + v[521];
	v[2465] = v[2468] + v[515] - v[520];
	v[2461] = v[2465] + v[520] - v[521];
	v[503] = v[4493] * v[487];
	v[2459] = v[503] + v[507] + v[508];
	v[2456] = v[2459] + v[502] - v[507];
	v[2452] = v[2456] + v[507] - v[508];
	v[2469] = v[2289] - v[2337] + 2e0*v[2468] + v[2412] * v[4482] + v[2400] * v[4486] + v[4488] * (v[2364] + v[4626]);
	v[2460] = v[2277] - v[2325] + 2e0*v[2459] + v[2406] * v[4482] + v[2394] * v[4486] + v[4488] * (v[2382] + v[4627]);
	v[2467] = -v[2229] + v[2313] + v[2424] * v[4482] + v[2400] * v[4489] + v[4485] * (v[2370] + v[4626]) + 0.5e0*v[475];
	v[2466] = -v[2228] + v[2312] + 2e0*v[2465] + v[2423] * v[4482] + (v[2204] + v[2252])*v[4489] + v[4485] * (v[4628]
		+ v[2351] * v[467] + v[2348] * v[883]);
	v[2458] = -v[2217] + v[2301] + v[2418] * v[4482] + v[2394] * v[4489] + v[4485] * (v[2388] + v[4627]) + 0.5e0*v[474];
	v[2457] = -v[2216] + v[2300] + 2e0*v[2456] + v[2417] * v[4482] + (v[2192] + v[2240])*v[4489] + v[4485] * (v[4632]
		+ v[2351] * v[466] + v[2342] * v[883]);
	v[2464] = v[2205] - v[2253] + (v[2364] + v[2370] + v[2439])*v[4483] + v[2424] * v[4486] + v[2412] * v[4489]
		+ 0.5e0*v[481];
	v[2463] = v[2204] - v[2252] + v[2423] * v[4486] + (v[2228] + v[2312])*v[4489] + 0.5e0*v[484] + v[4483] *
		(v[2351] * v[459] + v[4628] + v[2264] * v[883]);
	v[2462] = (v[2197] - v[2245])*v[243] + 2e0*v[2461] + v[4630] * v[563] + v[4486] * ((v[2281] + v[2329])*v[243]
		+ v[4624] * v[563]) + v[4489] * ((v[2221] + v[2305])*v[243] + v[4629] * v[563]) + v[4483] * (v[2350] * v[4631]
			+ v[4493] * (v[2179] * v[4492] + (e2i[0] * v[2095] + e1i[0] * v[2107])*v[4501] + v[2263] * v[4502] + v[4486] * v[4624]
				+ v[4489] * v[4629] + v[4630] - d[9] * v[4631] + v[2197] * v[471] + v[2221] * v[472] + v[2245] * v[473] + v[2281] * v[478]
				+ v[2305] * v[479] + v[2329] * v[482]) + v[2425] * v[488] + (v[2179] + v[2263])*v[883]);
	v[2455] = v[2193] - v[2241] + (v[2382] + v[2388] + v[2451])*v[4483] + v[2418] * v[4486] + v[2406] * v[4489]
		+ 0.5e0*v[480];
	v[2454] = v[2192] - v[2240] + v[2417] * v[4486] + (v[2216] + v[2300])*v[4489] + 0.5e0*v[483] + v[4483] *
		(v[2351] * v[458] + v[4632] + v[2258] * v[883]);
	v[2453] = (v[2185] - v[2233])*v[243] + 2e0*v[2452] + v[4634] * v[563] + v[4486] * ((v[2269] + v[2317])*v[243]
		+ v[4622] * v[563]) + v[4489] * ((v[2209] + v[2293])*v[243] + v[4633] * v[563]) + v[4483] * (v[2350] * v[4635]
			+ v[4493] * (v[2173] * v[4492] + (e2i[0] * v[2089] + e1i[0] * v[2101])*v[4501] + v[2257] * v[4502] + v[4486] * v[4622]
				+ v[4489] * v[4633] + v[4634] - d[9] * v[4635] + v[2185] * v[471] + v[2209] * v[472] + v[2233] * v[473] + v[2269] * v[478]
				+ v[2293] * v[479] + v[2317] * v[482]) + v[2425] * v[487] + (v[2173] + v[2257])*v[883]);
	v[499] = v[2452] * v[4483] + v[4489] * v[480] + v[4486] * v[483] + v[497] - v[498];
	v[4658] = v[499] / v[1784];
	v[504] = v[2456] * v[4485] + v[4489] * v[474] + v[4482] * v[483] - v[500] + v[501];
	v[4647] = v[504] / v[1784];
	v[509] = v[2459] * v[4488] + v[4486] * v[474] + v[4482] * v[480] + v[505] - v[506];
	v[4639] = v[509] / v[1784];
	v[512] = v[2461] * v[4483] + v[4489] * v[481] + v[4486] * v[484] + v[510] - v[511];
	v[4696] = -(v[419] * v[499]) + v[418] * v[512];
	v[4660] = v[417] * v[499] - v[416] * v[512];
	v[4659] = -(v[512] / v[1784]);
	v[517] = v[2465] * v[4485] + v[4489] * v[475] + v[4482] * v[484] - v[513] + v[514];
	v[4687] = -(v[419] * v[504]) + v[418] * v[517];
	v[4649] = v[417] * v[504] - v[416] * v[517];
	v[4648] = -(v[517] / v[1784]);
	v[522] = v[2468] * v[4488] + v[4486] * v[475] + v[4482] * v[481] + v[518] - v[519];
	v[4681] = -(v[419] * v[509]) + v[418] * v[522];
	v[4641] = v[417] * v[509] - v[416] * v[522];
	v[4640] = -(v[522] / v[1784]);
	v[1144] = v[1059] * v[1061] + v[1062] * v[1064] + v[1065] * v[1067] + v[3338] * v[364] + v[3329] * v[365]
		+ v[3320] * v[366];
	v[1135] = v[1059] * v[1060] + v[1062] * v[1063] + v[1065] * v[1066] + v[3336] * v[364] + v[3327] * v[365]
		+ v[3318] * v[366];
	v[1124] = v[3217] * v[364] + v[3211] * v[365] + v[3205] * v[366] + v[1059] * v[805] + v[1062] * v[807] + v[1065] * v[809];
	v[1123] = v[3220] * v[364] + v[3214] * v[365] + v[3208] * v[366] + v[1059] * v[804] + v[1062] * v[806] + v[1065] * v[808];
	v[1146] = v[1060] * v[1061] + v[1063] * v[1064] + v[1066] * v[1067] + v[3362] * v[364] + v[3354] * v[365]
		+ v[3346] * v[366];
	v[1132] = v[3218] * v[364] + v[3212] * v[365] + v[3206] * v[366] + v[1060] * v[805] + v[1063] * v[807] + v[1066] * v[809];
	v[1131] = v[3221] * v[364] + v[3215] * v[365] + v[3209] * v[366] + v[1060] * v[804] + v[1063] * v[806] + v[1066] * v[808];
	v[1141] = v[3219] * v[364] + v[3213] * v[365] + v[3207] * v[366] + v[1061] * v[805] + v[1064] * v[807] + v[1067] * v[809];
	v[1140] = v[3222] * v[364] + v[3216] * v[365] + v[3210] * v[366] + v[1061] * v[804] + v[1064] * v[806] + v[1067] * v[808];
	v[2601] = (-(v[2469] * v[416]) + v[2460] * v[417] + (v[2001] + v[417] * v[4636])*v[509] - (v[1998] + v[416] * v[4636]
		)*v[522]) / v[1784];
	v[2599] = v[2467] * v[4637];
	v[2598] = v[2458] * v[4638];
	v[4642] = v[2598] + v[2599];
	v[2600] = v[4642] + v[4655] * v[509] + v[4654] * v[522];
	v[2596] = v[2464] * v[4637];
	v[2595] = v[2455] * v[4638];
	v[4650] = v[2595] + v[2596];
	v[2597] = v[4650] + v[4646] * v[509] + v[4645] * v[522];
	v[2593] = v[4637] * v[664];
	v[2592] = v[4638] * v[643];
	v[4661] = v[2592] + v[2593];
	v[2594] = v[410] * v[4639] + v[404] * v[4640] + v[2475] * v[4641] + v[4661];
	v[2590] = v[4637] * v[658];
	v[2589] = v[4638] * v[640];
	v[4667] = v[2589] + v[2590];
	v[2591] = v[411] * v[4639] + v[405] * v[4640] + v[2473] * v[4641] + v[4667];
	v[2587] = v[4637] * v[655];
	v[2586] = v[4638] * v[637];
	v[4676] = v[2586] + v[2587];
	v[2588] = v[412] * v[4639] + v[406] * v[4640] + v[2471] * v[4641] + v[4676];
	v[2585] = v[4642] + v[4652] * v[504] + v[4651] * v[517];
	v[2584] = (-(v[2466] * v[416]) + v[2457] * v[417] + (v[2000] + v[1784] * v[4643])*v[504] - (v[1997] + v[1784] * v[4644]
		)*v[517]) / v[1784];
	v[2582] = v[2463] * v[4637];
	v[2581] = v[2454] * v[4638];
	v[4653] = v[2581] + v[2582];
	v[2583] = v[4653] + v[4646] * v[504] + v[4645] * v[517];
	v[2579] = v[4637] * v[663];
	v[2578] = v[4638] * v[642];
	v[4662] = v[2578] + v[2579];
	v[2580] = v[410] * v[4647] + v[404] * v[4648] + v[2475] * v[4649] + v[4662];
	v[2576] = v[4637] * v[657];
	v[2575] = v[4638] * v[639];
	v[4668] = v[2575] + v[2576];
	v[2577] = v[411] * v[4647] + v[405] * v[4648] + v[2473] * v[4649] + v[4668];
	v[2573] = v[4637] * v[654];
	v[2572] = v[4638] * v[636];
	v[4677] = v[2572] + v[2573];
	v[2574] = v[412] * v[4647] + v[406] * v[4648] + v[2471] * v[4649] + v[4677];
	v[2571] = v[4650] + v[4652] * v[499] + v[4651] * v[512];
	v[2570] = v[4653] + v[4655] * v[499] + v[4654] * v[512];
	v[2569] = (-(v[2462] * v[416]) + v[2453] * v[417] + (v[1999] + v[1784] * v[4656])*v[499] - (v[1996] + v[1784] * v[4657]
		)*v[512]) / v[1784];
	v[2567] = v[4637] * v[662];
	v[2566] = v[4638] * v[641];
	v[4663] = v[2566] + v[2567];
	v[2568] = v[410] * v[4658] + v[404] * v[4659] + v[2475] * v[4660] + v[4663];
	v[2564] = v[4637] * v[656];
	v[2563] = v[4638] * v[638];
	v[4669] = v[2563] + v[2564];
	v[2565] = v[411] * v[4658] + v[405] * v[4659] + v[2473] * v[4660] + v[4669];
	v[2561] = v[4637] * v[653];
	v[2560] = v[4638] * v[635];
	v[4678] = v[2560] + v[2561];
	v[2562] = v[412] * v[4658] + v[406] * v[4659] + v[2471] * v[4660] + v[4678];
	v[2559] = v[356] * v[4651] + v[346] * v[4652] + v[4661];
	v[2558] = v[356] * v[4654] + v[346] * v[4655] + v[4662];
	v[2557] = v[356] * v[4645] + v[346] * v[4646] + v[4663];
	v[2556] = (-(v[356] * v[4670]) + v[346] * v[4671]) / v[1784];
	v[2555] = (-(v[356] * v[4672]) + v[346] * v[4673]) / v[1784];
	v[2554] = (-(v[356] * v[4674]) + v[346] * v[4675]) / v[1784];
	v[2553] = v[355] * v[4651] + v[345] * v[4652] + v[4667];
	v[2552] = v[355] * v[4654] + v[345] * v[4655] + v[4668];
	v[2551] = v[355] * v[4645] + v[345] * v[4646] + v[4669];
	v[2550] = (-(v[355] * v[4670]) + v[345] * v[4671]) / v[1784];
	v[2549] = (-(v[355] * v[4672]) + v[345] * v[4673]) / v[1784];
	v[2548] = (-(v[355] * v[4674]) + v[345] * v[4675]) / v[1784];
	v[2547] = v[353] * v[4651] + v[343] * v[4652] + v[4676];
	v[2546] = v[353] * v[4654] + v[343] * v[4655] + v[4677];
	v[2545] = v[353] * v[4645] + v[343] * v[4646] + v[4678];
	v[2544] = (-(v[353] * v[4670]) + v[343] * v[4671]) / v[1784];
	v[2543] = (-(v[353] * v[4672]) + v[343] * v[4673]) / v[1784];
	v[2542] = (-(v[353] * v[4674]) + v[343] * v[4675]) / v[1784];
	v[2541] = (v[2469] * v[418] - v[2460] * v[419] - (v[2007] + v[419] * v[4636])*v[509] + (v[2004] + v[418] * v[4636]
		)*v[522]) / v[1784];
	v[2539] = v[2467] * v[4679];
	v[2538] = v[2458] * v[4680];
	v[4682] = v[2538] + v[2539];
	v[2540] = v[4682] + v[4693] * v[509] + v[4692] * v[522];
	v[2536] = v[2464] * v[4679];
	v[2535] = v[2455] * v[4680];
	v[4688] = v[2535] + v[2536];
	v[2537] = v[4688] + v[4686] * v[509] + v[4685] * v[522];
	v[2533] = v[4679] * v[664];
	v[2532] = v[4680] * v[643];
	v[4697] = v[2532] + v[2533];
	v[2534] = -(v[413] * v[4639]) - v[407] * v[4640] + v[2475] * v[4681] + v[4697];
	v[3380] = v[2534] * v[804] + v[2594] * v[805];
	v[3373] = v[2534] * v[806] + v[2594] * v[807];
	v[3366] = v[2534] * v[808] + v[2594] * v[809];
	v[2530] = v[4679] * v[658];
	v[2529] = v[4680] * v[640];
	v[4700] = v[2529] + v[2530];
	v[2531] = -(v[414] * v[4639]) - v[408] * v[4640] + v[2473] * v[4681] + v[4700];
	v[3379] = v[2531] * v[804] + v[2591] * v[805];
	v[3372] = v[2531] * v[806] + v[2591] * v[807];
	v[3365] = v[2531] * v[808] + v[2591] * v[809];
	v[2527] = v[4679] * v[655];
	v[2526] = v[4680] * v[637];
	v[4709] = v[2526] + v[2527];
	v[2528] = -(v[415] * v[4639]) - v[409] * v[4640] + v[2471] * v[4681] + v[4709];
	v[3378] = v[2528] * v[804] + v[2588] * v[805];
	v[3371] = v[2528] * v[806] + v[2588] * v[807];
	v[3364] = v[2528] * v[808] + v[2588] * v[809];
	v[2525] = v[4682] + v[4690] * v[504] + v[4689] * v[517];
	v[2524] = (v[2466] * v[418] - v[2457] * v[419] - (v[2006] + v[1784] * v[4683])*v[504] + (v[2003] + v[1784] * v[4684]
		)*v[517]) / v[1784];
	v[2522] = v[2463] * v[4679];
	v[2521] = v[2454] * v[4680];
	v[4691] = v[2521] + v[2522];
	v[2523] = v[4691] + v[4686] * v[504] + v[4685] * v[517];
	v[2519] = v[4679] * v[663];
	v[2518] = v[4680] * v[642];
	v[4698] = v[2518] + v[2519];
	v[2520] = -(v[413] * v[4647]) - v[407] * v[4648] + v[2475] * v[4687] + v[4698];
	v[3358] = v[2520] * v[804] + v[2580] * v[805];
	v[3350] = v[2520] * v[806] + v[2580] * v[807];
	v[3342] = v[2520] * v[808] + v[2580] * v[809];
	v[2516] = v[4679] * v[657];
	v[2515] = v[4680] * v[639];
	v[4701] = v[2515] + v[2516];
	v[2517] = -(v[414] * v[4647]) - v[408] * v[4648] + v[2473] * v[4687] + v[4701];
	v[3357] = v[2517] * v[804] + v[2577] * v[805];
	v[3349] = v[2517] * v[806] + v[2577] * v[807];
	v[3341] = v[2517] * v[808] + v[2577] * v[809];
	v[2513] = v[4679] * v[654];
	v[2512] = v[4680] * v[636];
	v[4710] = v[2512] + v[2513];
	v[2514] = -(v[415] * v[4647]) - v[409] * v[4648] + v[2471] * v[4687] + v[4710];
	v[3356] = v[2514] * v[804] + v[2574] * v[805];
	v[3348] = v[2514] * v[806] + v[2574] * v[807];
	v[3340] = v[2514] * v[808] + v[2574] * v[809];
	v[2511] = v[4688] + v[4690] * v[499] + v[4689] * v[512];
	v[2510] = v[4691] + v[4693] * v[499] + v[4692] * v[512];
	v[2509] = (v[2462] * v[418] - v[2453] * v[419] - (v[2005] + v[1784] * v[4694])*v[499] + (v[2002] + v[1784] * v[4695]
		)*v[512]) / v[1784];
	v[2507] = v[4679] * v[662];
	v[2506] = v[4680] * v[641];
	v[4699] = v[2506] + v[2507];
	v[2508] = -(v[413] * v[4658]) - v[407] * v[4659] + v[2475] * v[4696] + v[4699];
	v[3333] = v[2508] * v[804] + v[2568] * v[805];
	v[3324] = v[2508] * v[806] + v[2568] * v[807];
	v[3315] = v[2508] * v[808] + v[2568] * v[809];
	v[2504] = v[4679] * v[656];
	v[2503] = v[4680] * v[638];
	v[4702] = v[2503] + v[2504];
	v[2505] = -(v[414] * v[4658]) - v[408] * v[4659] + v[2473] * v[4696] + v[4702];
	v[3332] = v[2505] * v[804] + v[2565] * v[805];
	v[3323] = v[2505] * v[806] + v[2565] * v[807];
	v[3314] = v[2505] * v[808] + v[2565] * v[809];
	v[2501] = v[4679] * v[653];
	v[2500] = v[4680] * v[635];
	v[4711] = v[2500] + v[2501];
	v[2502] = -(v[415] * v[4658]) - v[409] * v[4659] + v[2471] * v[4696] + v[4711];
	v[3331] = v[2502] * v[804] + v[2562] * v[805];
	v[3322] = v[2502] * v[806] + v[2562] * v[807];
	v[3313] = v[2502] * v[808] + v[2562] * v[809];
	v[2499] = v[356] * v[4689] + v[346] * v[4690] + v[4697];
	v[2498] = v[356] * v[4692] + v[346] * v[4693] + v[4698];
	v[2497] = v[356] * v[4685] + v[346] * v[4686] + v[4699];
	v[2496] = (v[356] * v[4703] - v[346] * v[4704]) / v[1784];
	v[3309] = v[2496] * v[804] + v[2556] * v[805];
	v[3303] = v[2496] * v[806] + v[2556] * v[807];
	v[3297] = v[2496] * v[808] + v[2556] * v[809];
	v[2495] = (v[356] * v[4705] - v[346] * v[4706]) / v[1784];
	v[3308] = v[2495] * v[804] + v[2555] * v[805];
	v[3302] = v[2495] * v[806] + v[2555] * v[807];
	v[3296] = v[2495] * v[808] + v[2555] * v[809];
	v[2494] = (v[356] * v[4707] - v[346] * v[4708]) / v[1784];
	v[3307] = v[2494] * v[804] + v[2554] * v[805];
	v[3301] = v[2494] * v[806] + v[2554] * v[807];
	v[3295] = v[2494] * v[808] + v[2554] * v[809];
	v[2493] = v[355] * v[4689] + v[345] * v[4690] + v[4700];
	v[2492] = v[355] * v[4692] + v[345] * v[4693] + v[4701];
	v[2491] = v[355] * v[4685] + v[345] * v[4686] + v[4702];
	v[2490] = (v[355] * v[4703] - v[345] * v[4704]) / v[1784];
	v[3291] = v[2490] * v[804] + v[2550] * v[805];
	v[3285] = v[2490] * v[806] + v[2550] * v[807];
	v[3279] = v[2490] * v[808] + v[2550] * v[809];
	v[2489] = (v[355] * v[4705] - v[345] * v[4706]) / v[1784];
	v[3290] = v[2489] * v[804] + v[2549] * v[805];
	v[3284] = v[2489] * v[806] + v[2549] * v[807];
	v[3278] = v[2489] * v[808] + v[2549] * v[809];
	v[2488] = (v[355] * v[4707] - v[345] * v[4708]) / v[1784];
	v[3289] = v[2488] * v[804] + v[2548] * v[805];
	v[3283] = v[2488] * v[806] + v[2548] * v[807];
	v[3277] = v[2488] * v[808] + v[2548] * v[809];
	v[2487] = v[353] * v[4689] + v[343] * v[4690] + v[4709];
	v[2486] = v[353] * v[4692] + v[343] * v[4693] + v[4710];
	v[2485] = v[353] * v[4685] + v[343] * v[4686] + v[4711];
	v[2484] = (v[353] * v[4703] - v[343] * v[4704]) / v[1784];
	v[3273] = v[2484] * v[804] + v[2544] * v[805];
	v[3267] = v[2484] * v[806] + v[2544] * v[807];
	v[3261] = v[2484] * v[808] + v[2544] * v[809];
	v[2483] = (v[353] * v[4705] - v[343] * v[4706]) / v[1784];
	v[3272] = v[2483] * v[804] + v[2543] * v[805];
	v[3266] = v[2483] * v[806] + v[2543] * v[807];
	v[3260] = v[2483] * v[808] + v[2543] * v[809];
	v[2482] = (v[353] * v[4707] - v[343] * v[4708]) / v[1784];
	v[3271] = v[2482] * v[804] + v[2542] * v[805];
	v[3265] = v[2482] * v[806] + v[2542] * v[807];
	v[3259] = v[2482] * v[808] + v[2542] * v[809];
	v[1822] = -(v[1610] * v[4639]);
	v[1821] = -(v[1608] * v[4640]);
	v[1823] = -v[1821] - v[1822] + v[1787] * v[4641] + v[417] * v[4712] + v[416] * v[4713];
	v[1820] = (-(v[2001] * v[416]) + v[1998] * v[417] + v[4718] * v[509] - v[4717] * v[522]) / v[1784];
	v[1818] = -(v[1610] * v[4647]);
	v[1817] = -(v[1608] * v[4648]);
	v[1819] = -v[1817] - v[1818] + v[1787] * v[4649] + v[417] * v[4715] + v[416] * v[4716];
	v[1816] = (-(v[2000] * v[416]) + v[1997] * v[417] + v[4718] * v[504] - v[4717] * v[517]) / v[1784];
	v[1814] = -(v[1610] * v[4658]);
	v[1813] = -(v[1608] * v[4659]);
	v[1815] = -v[1813] - v[1814] + v[1787] * v[4660] + v[417] * v[4719] + v[416] * v[4720];
	v[1812] = (-(v[1999] * v[416]) + v[1996] * v[417] + v[4718] * v[499] - v[4717] * v[512]) / v[1784];
	v[1810] = v[346] * v[4725];
	v[1809] = v[356] * v[4726];
	v[1811] = (-v[4721] - v[4722] - v[1784] * (v[1809] + v[1810] - v[1787] * v[4723] + v[1787] * v[4724])) / v[1784];
	v[1808] = (v[1609] * v[346] - v[1607] * v[356] - v[410] * v[416] + v[404] * v[417] + v[4714] * v[4723] - v[4714] * v[4724])
		/ v[1784];
	v[1806] = v[345] * v[4725];
	v[1805] = v[355] * v[4726];
	v[1807] = (-v[4727] - v[4728] - v[1784] * (v[1805] + v[1806] - v[1787] * v[4729] + v[1787] * v[4730])) / v[1784];
	v[1804] = (v[1609] * v[345] - v[1607] * v[355] - v[411] * v[416] + v[405] * v[417] + v[4714] * v[4729] - v[4714] * v[4730])
		/ v[1784];
	v[1802] = v[343] * v[4725];
	v[1801] = v[353] * v[4726];
	v[1803] = (-v[4731] - v[4732] - v[1784] * (v[1801] + v[1802] - v[1787] * v[4733] + v[1787] * v[4734])) / v[1784];
	v[1800] = (v[1609] * v[343] - v[1607] * v[353] - v[412] * v[416] + v[406] * v[417] + v[4714] * v[4733] - v[4714] * v[4734])
		/ v[1784];
	v[1799] = (v[2007] * v[418] - v[2004] * v[419] - v[4739] * v[509] + v[4738] * v[522]) / v[1784];
	v[1798] = v[1821] + v[1822] + v[419] * (v[4736] - v[1785] * v[509]) + v[418] * (v[4737] + v[1785] * v[522]);
	v[1797] = (v[2006] * v[418] - v[2003] * v[419] - v[4739] * v[504] + v[4738] * v[517]) / v[1784];
	v[1796] = v[1817] + v[1818] + v[419] * (v[4740] - v[1785] * v[504]) + v[418] * (v[4741] + v[1785] * v[517]);
	v[1795] = (v[2005] * v[418] - v[2002] * v[419] - v[4739] * v[499] + v[4738] * v[512]) / v[1784];
	v[1794] = v[1813] + v[1814] + v[419] * (v[4742] - v[1785] * v[499]) + v[418] * (v[4743] + v[1785] * v[512]);
	v[1793] = (-(v[1612] * v[346]) + v[1611] * v[356] + v[413] * v[418] - v[407] * v[419] - v[346] * v[4744] + v[356] * v[4745]
		) / v[1784];
	v[1792] = v[1809] + v[1810] + v[4746] / v[1784] - v[4747] / v[1784] + v[1785] * v[4754];
	v[1791] = (-(v[1612] * v[345]) + v[1611] * v[355] + v[414] * v[418] - v[408] * v[419] - v[345] * v[4744] + v[355] * v[4745]
		) / v[1784];
	v[1790] = v[1805] + v[1806] + v[4748] / v[1784] - v[4749] / v[1784] + v[1785] * v[4753];
	v[1789] = (-(v[1612] * v[343]) + v[1611] * v[353] + v[415] * v[418] - v[409] * v[419] - v[343] * v[4744] + v[353] * v[4745]
		) / v[1784];
	v[1788] = v[1801] + v[1802] + v[4750] / v[1784] - v[4751] / v[1784] + v[1785] * v[4752];
	v[699] = v[4752] / v[1784];
	v[3036] = v[2696] * v[699];
	v[3034] = v[2688] * v[699];
	v[3030] = v[2687] * v[699];
	v[700] = v[4753] / v[1784];
	v[3978] = v[2688] * v[700];
	v[3891] = v[2687] * v[700];
	v[3049] = v[2696] * v[700];
	v[701] = v[4754] / v[1784];
	v[3965] = v[2696] * v[701];
	v[3962] = v[2688] * v[701];
	v[3877] = v[2687] * v[701];
	v[702] = v[4696] / v[1784];
	v[703] = v[4687] / v[1784];
	v[704] = v[4681] / v[1784];
	v[705] = (v[4733] - v[4734]) / v[1784];
	v[3997] = v[2706] * v[699] + v[2730] * v[705];
	v[3996] = v[2705] * v[699] + v[2729] * v[705];
	v[3995] = v[2704] * v[699] + v[2728] * v[705];
	v[3910] = v[2699] * v[699] + v[2723] * v[705];
	v[3909] = v[2698] * v[699] + v[2722] * v[705];
	v[3908] = v[2697] * v[699] + v[2721] * v[705];
	v[3828] = v[2691] * v[699] + v[2715] * v[705];
	v[3827] = v[2690] * v[699] + v[2714] * v[705];
	v[3826] = v[2689] * v[699] + v[2713] * v[705];
	v[3617] = v[3595] * v[699] + v[3596] * v[705] + v[1789] * v[804] + v[1803] * v[805];
	v[3616] = v[3597] * v[699] + v[3595] * v[705] + v[1788] * v[804] + v[1800] * v[805];
	v[3615] = v[3592] * v[699] + v[3593] * v[705] + v[1789] * v[806] + v[1803] * v[807];
	v[3614] = v[3594] * v[699] + v[3592] * v[705] + v[1788] * v[806] + v[1800] * v[807];
	v[3613] = v[3589] * v[699] + v[3590] * v[705] + v[1789] * v[808] + v[1803] * v[809];
	v[3612] = v[3591] * v[699] + v[3589] * v[705] + v[1788] * v[808] + v[1800] * v[809];
	v[3276] = v[3222] * v[699] + v[3219] * v[705] + v[2487] * v[804] + v[2547] * v[805];
	v[3275] = v[3221] * v[699] + v[3218] * v[705] + v[2486] * v[804] + v[2546] * v[805];
	v[3274] = v[3220] * v[699] + v[3217] * v[705] + v[2485] * v[804] + v[2545] * v[805];
	v[3270] = v[3216] * v[699] + v[3213] * v[705] + v[2487] * v[806] + v[2547] * v[807];
	v[3269] = v[3215] * v[699] + v[3212] * v[705] + v[2486] * v[806] + v[2546] * v[807];
	v[3268] = v[3214] * v[699] + v[3211] * v[705] + v[2485] * v[806] + v[2545] * v[807];
	v[3264] = v[3210] * v[699] + v[3207] * v[705] + v[2487] * v[808] + v[2547] * v[809];
	v[3263] = v[3209] * v[699] + v[3206] * v[705] + v[2486] * v[808] + v[2546] * v[809];
	v[3262] = v[3208] * v[699] + v[3205] * v[705] + v[2485] * v[808] + v[2545] * v[809];
	v[3167] = v[699] * v[808] + v[705] * v[809];
	v[3160] = v[699] * v[806] + v[705] * v[807];
	v[3153] = 1e0 + v[699] * v[804] + v[705] * v[805];
	v[3037] = v[2720] * v[705];
	v[4771] = v[3036] + v[3037];
	v[3035] = v[2712] * v[705];
	v[4772] = v[3034] + v[3035];
	v[3031] = v[2711] * v[705];
	v[4768] = v[3030] + v[3031];
	v[706] = (v[4729] - v[4730]) / v[1784];
	v[3985] = v[2706] * v[700] + v[2730] * v[706];
	v[3984] = v[2705] * v[700] + v[2729] * v[706];
	v[3983] = v[2704] * v[700] + v[2728] * v[706];
	v[3979] = v[2712] * v[706];
	v[4774] = v[3978] + v[3979];
	v[3898] = v[2699] * v[700] + v[2723] * v[706];
	v[3897] = v[2698] * v[700] + v[2722] * v[706];
	v[3896] = v[2697] * v[700] + v[2721] * v[706];
	v[3892] = v[2711] * v[706];
	v[4769] = v[3891] + v[3892];
	v[3816] = v[2691] * v[700] + v[2715] * v[706];
	v[3815] = v[2690] * v[700] + v[2714] * v[706];
	v[3814] = v[2689] * v[700] + v[2713] * v[706];
	v[3623] = v[3595] * v[700] + v[3596] * v[706] + v[1791] * v[804] + v[1807] * v[805];
	v[3622] = v[3597] * v[700] + v[3595] * v[706] + v[1790] * v[804] + v[1804] * v[805];
	v[3621] = v[3592] * v[700] + v[3593] * v[706] + v[1791] * v[806] + v[1807] * v[807];
	v[3620] = v[3594] * v[700] + v[3592] * v[706] + v[1790] * v[806] + v[1804] * v[807];
	v[3619] = v[3589] * v[700] + v[3590] * v[706] + v[1791] * v[808] + v[1807] * v[809];
	v[3618] = v[3591] * v[700] + v[3589] * v[706] + v[1790] * v[808] + v[1804] * v[809];
	v[3294] = v[3222] * v[700] + v[3219] * v[706] + v[2493] * v[804] + v[2553] * v[805];
	v[3293] = v[3221] * v[700] + v[3218] * v[706] + v[2492] * v[804] + v[2552] * v[805];
	v[3292] = v[3220] * v[700] + v[3217] * v[706] + v[2491] * v[804] + v[2551] * v[805];
	v[3288] = v[3216] * v[700] + v[3213] * v[706] + v[2493] * v[806] + v[2553] * v[807];
	v[3287] = v[3215] * v[700] + v[3212] * v[706] + v[2492] * v[806] + v[2552] * v[807];
	v[3286] = v[3214] * v[700] + v[3211] * v[706] + v[2491] * v[806] + v[2551] * v[807];
	v[3282] = v[3210] * v[700] + v[3207] * v[706] + v[2493] * v[808] + v[2553] * v[809];
	v[3281] = v[3209] * v[700] + v[3206] * v[706] + v[2492] * v[808] + v[2552] * v[809];
	v[3280] = v[3208] * v[700] + v[3205] * v[706] + v[2491] * v[808] + v[2551] * v[809];
	v[3168] = v[700] * v[808] + v[706] * v[809];
	v[3161] = 1e0 + v[700] * v[806] + v[706] * v[807];
	v[3154] = v[700] * v[804] + v[706] * v[805];
	v[3050] = v[2720] * v[706];
	v[4773] = v[3049] + v[3050];
	v[707] = (v[4723] - v[4724]) / v[1784];
	v[3971] = v[2706] * v[701] + v[2730] * v[707];
	v[3970] = v[2705] * v[701] + v[2729] * v[707];
	v[3969] = v[2704] * v[701] + v[2728] * v[707];
	v[3966] = v[2720] * v[707];
	v[4775] = v[3965] + v[3966];
	v[3963] = v[2712] * v[707];
	v[4776] = v[3962] + v[3963];
	v[3884] = v[2699] * v[701] + v[2723] * v[707];
	v[3883] = v[2698] * v[701] + v[2722] * v[707];
	v[3882] = v[2697] * v[701] + v[2721] * v[707];
	v[3878] = v[2711] * v[707];
	v[4770] = v[3877] + v[3878];
	v[3804] = v[2691] * v[701] + v[2715] * v[707];
	v[3803] = v[2690] * v[701] + v[2714] * v[707];
	v[3802] = v[2689] * v[701] + v[2713] * v[707];
	v[3629] = v[3595] * v[701] + v[3596] * v[707] + v[1793] * v[804] + v[1811] * v[805];
	v[3628] = v[3597] * v[701] + v[3595] * v[707] + v[1792] * v[804] + v[1808] * v[805];
	v[3627] = v[3592] * v[701] + v[3593] * v[707] + v[1793] * v[806] + v[1811] * v[807];
	v[3626] = v[3594] * v[701] + v[3592] * v[707] + v[1792] * v[806] + v[1808] * v[807];
	v[3625] = v[3589] * v[701] + v[3590] * v[707] + v[1793] * v[808] + v[1811] * v[809];
	v[3624] = v[3591] * v[701] + v[3589] * v[707] + v[1792] * v[808] + v[1808] * v[809];
	v[3312] = v[3222] * v[701] + v[3219] * v[707] + v[2499] * v[804] + v[2559] * v[805];
	v[3311] = v[3221] * v[701] + v[3218] * v[707] + v[2498] * v[804] + v[2558] * v[805];
	v[3310] = v[3220] * v[701] + v[3217] * v[707] + v[2497] * v[804] + v[2557] * v[805];
	v[3306] = v[3216] * v[701] + v[3213] * v[707] + v[2499] * v[806] + v[2559] * v[807];
	v[3305] = v[3215] * v[701] + v[3212] * v[707] + v[2498] * v[806] + v[2558] * v[807];
	v[3304] = v[3214] * v[701] + v[3211] * v[707] + v[2497] * v[806] + v[2557] * v[807];
	v[3300] = v[3210] * v[701] + v[3207] * v[707] + v[2499] * v[808] + v[2559] * v[809];
	v[3299] = v[3209] * v[701] + v[3206] * v[707] + v[2498] * v[808] + v[2558] * v[809];
	v[3298] = v[3208] * v[701] + v[3205] * v[707] + v[2497] * v[808] + v[2557] * v[809];
	v[3169] = 1e0 + v[701] * v[808] + v[707] * v[809];
	v[3162] = v[701] * v[806] + v[707] * v[807];
	v[3155] = v[701] * v[804] + v[707] * v[805];
	v[708] = v[4660] / v[1784];
	v[3954] = v[3757] + v[2706] * v[702] + v[2730] * v[708];
	v[3953] = v[3745] + v[2705] * v[702] + v[2729] * v[708];
	v[3952] = v[3730] + v[2704] * v[702] + v[2728] * v[708];
	v[3949] = v[3020] + v[2696] * v[702] + v[2720] * v[708];
	v[3947] = v[3014] + v[2688] * v[702] + v[2712] * v[708];
	v[3869] = v[3753] + v[2699] * v[702] + v[2723] * v[708];
	v[3868] = v[3740] + v[2698] * v[702] + v[2722] * v[708];
	v[3867] = v[3724] + v[2697] * v[702] + v[2721] * v[708];
	v[3863] = v[3011] + v[2687] * v[702] + v[2711] * v[708];
	v[3791] = v[3749] + v[2691] * v[702] + v[2715] * v[708];
	v[3790] = v[3735] + v[2690] * v[702] + v[2714] * v[708];
	v[3789] = v[3718] + v[2689] * v[702] + v[2713] * v[708];
	v[3635] = v[3217] + v[3595] * v[702] + v[3596] * v[708] + v[1795] * v[804] + v[1815] * v[805];
	v[3634] = v[3220] + v[3597] * v[702] + v[3595] * v[708] + v[1794] * v[804] + v[1812] * v[805];
	v[3633] = v[3211] + v[3592] * v[702] + v[3593] * v[708] + v[1795] * v[806] + v[1815] * v[807];
	v[3632] = v[3214] + v[3594] * v[702] + v[3592] * v[708] + v[1794] * v[806] + v[1812] * v[807];
	v[3631] = v[3205] + v[3589] * v[702] + v[3590] * v[708] + v[1795] * v[808] + v[1815] * v[809];
	v[3630] = v[3208] + v[3591] * v[702] + v[3589] * v[708] + v[1794] * v[808] + v[1812] * v[809];
	v[3339] = v[3338] + v[3222] * v[702] + v[3219] * v[708] + v[2511] * v[804] + v[2571] * v[805];
	v[3337] = v[3336] + v[3221] * v[702] + v[3218] * v[708] + v[2510] * v[804] + v[2570] * v[805];
	v[3335] = v[3334] + v[3220] * v[702] + v[3217] * v[708] + v[2509] * v[804] + v[2569] * v[805];
	v[3330] = v[3329] + v[3216] * v[702] + v[3213] * v[708] + v[2511] * v[806] + v[2571] * v[807];
	v[3328] = v[3327] + v[3215] * v[702] + v[3212] * v[708] + v[2510] * v[806] + v[2570] * v[807];
	v[3326] = v[3325] + v[3214] * v[702] + v[3211] * v[708] + v[2509] * v[806] + v[2569] * v[807];
	v[3321] = v[3320] + v[3210] * v[702] + v[3207] * v[708] + v[2511] * v[808] + v[2571] * v[809];
	v[3319] = v[3318] + v[3209] * v[702] + v[3206] * v[708] + v[2510] * v[808] + v[2570] * v[809];
	v[3317] = v[3316] + v[3208] * v[702] + v[3205] * v[708] + v[2509] * v[808] + v[2569] * v[809];
	v[3170] = v[1065] + v[702] * v[808] + v[708] * v[809];
	v[3163] = v[1062] + v[702] * v[806] + v[708] * v[807];
	v[3156] = v[1059] + v[702] * v[804] + v[708] * v[805];
	v[709] = v[4649] / v[1784];
	v[3939] = v[3758] + v[2706] * v[703] + v[2730] * v[709];
	v[3938] = v[3746] + v[2705] * v[703] + v[2729] * v[709];
	v[3937] = v[3731] + v[2704] * v[703] + v[2728] * v[709];
	v[3934] = v[3021] + v[2696] * v[703] + v[2720] * v[709];
	v[3932] = v[3015] + v[2688] * v[703] + v[2712] * v[709];
	v[3855] = v[3754] + v[2699] * v[703] + v[2723] * v[709];
	v[3854] = v[3741] + v[2698] * v[703] + v[2722] * v[709];
	v[3853] = v[3725] + v[2697] * v[703] + v[2721] * v[709];
	v[3849] = v[3012] + v[2687] * v[703] + v[2711] * v[709];
	v[3778] = v[3750] + v[2691] * v[703] + v[2715] * v[709];
	v[3777] = v[3736] + v[2690] * v[703] + v[2714] * v[709];
	v[3776] = v[3719] + v[2689] * v[703] + v[2713] * v[709];
	v[3641] = v[3218] + v[3595] * v[703] + v[3596] * v[709] + v[1797] * v[804] + v[1819] * v[805];
	v[3640] = v[3221] + v[3597] * v[703] + v[3595] * v[709] + v[1796] * v[804] + v[1816] * v[805];
	v[3639] = v[3212] + v[3592] * v[703] + v[3593] * v[709] + v[1797] * v[806] + v[1819] * v[807];
	v[3638] = v[3215] + v[3594] * v[703] + v[3592] * v[709] + v[1796] * v[806] + v[1816] * v[807];
	v[3637] = v[3206] + v[3589] * v[703] + v[3590] * v[709] + v[1797] * v[808] + v[1819] * v[809];
	v[3636] = v[3209] + v[3591] * v[703] + v[3589] * v[709] + v[1796] * v[808] + v[1816] * v[809];
	v[3363] = v[3362] + v[3222] * v[703] + v[3219] * v[709] + v[2525] * v[804] + v[2585] * v[805];
	v[3361] = v[3360] + v[3221] * v[703] + v[3218] * v[709] + v[2524] * v[804] + v[2584] * v[805];
	v[3359] = v[3336] + v[3220] * v[703] + v[3217] * v[709] + v[2523] * v[804] + v[2583] * v[805];
	v[3355] = v[3354] + v[3216] * v[703] + v[3213] * v[709] + v[2525] * v[806] + v[2585] * v[807];
	v[3353] = v[3352] + v[3215] * v[703] + v[3212] * v[709] + v[2524] * v[806] + v[2584] * v[807];
	v[3351] = v[3327] + v[3214] * v[703] + v[3211] * v[709] + v[2523] * v[806] + v[2583] * v[807];
	v[3347] = v[3346] + v[3210] * v[703] + v[3207] * v[709] + v[2525] * v[808] + v[2585] * v[809];
	v[3345] = v[3344] + v[3209] * v[703] + v[3206] * v[709] + v[2524] * v[808] + v[2584] * v[809];
	v[3343] = v[3318] + v[3208] * v[703] + v[3205] * v[709] + v[2523] * v[808] + v[2583] * v[809];
	v[3171] = v[1066] + v[703] * v[808] + v[709] * v[809];
	v[3164] = v[1063] + v[703] * v[806] + v[709] * v[807];
	v[3157] = v[1060] + v[703] * v[804] + v[709] * v[805];
	v[710] = v[4641] / v[1784];
	v[3924] = v[3759] + v[2706] * v[704] + v[2730] * v[710];
	v[3923] = v[3747] + v[2705] * v[704] + v[2729] * v[710];
	v[3922] = v[3732] + v[2704] * v[704] + v[2728] * v[710];
	v[3919] = v[3022] + v[2696] * v[704] + v[2720] * v[710];
	v[3917] = v[3016] + v[2688] * v[704] + v[2712] * v[710];
	v[3841] = v[3755] + v[2699] * v[704] + v[2723] * v[710];
	v[3840] = v[3742] + v[2698] * v[704] + v[2722] * v[710];
	v[3839] = v[3726] + v[2697] * v[704] + v[2721] * v[710];
	v[3835] = v[3013] + v[2687] * v[704] + v[2711] * v[710];
	v[3765] = v[3751] + v[2691] * v[704] + v[2715] * v[710];
	v[3764] = v[3737] + v[2690] * v[704] + v[2714] * v[710];
	v[3763] = v[3720] + v[2689] * v[704] + v[2713] * v[710];
	v[3647] = v[3219] + v[3595] * v[704] + v[3596] * v[710] + v[1799] * v[804] + v[1823] * v[805];
	v[3646] = v[3222] + v[3597] * v[704] + v[3595] * v[710] + v[1798] * v[804] + v[1820] * v[805];
	v[3645] = v[3213] + v[3592] * v[704] + v[3593] * v[710] + v[1799] * v[806] + v[1823] * v[807];
	v[3644] = v[3216] + v[3594] * v[704] + v[3592] * v[710] + v[1798] * v[806] + v[1820] * v[807];
	v[3643] = v[3207] + v[3589] * v[704] + v[3590] * v[710] + v[1799] * v[808] + v[1823] * v[809];
	v[3642] = v[3210] + v[3591] * v[704] + v[3589] * v[710] + v[1798] * v[808] + v[1820] * v[809];
	v[3384] = v[3383] + v[3222] * v[704] + v[3219] * v[710] + v[2541] * v[804] + v[2601] * v[805];
	v[3382] = v[3362] + v[3221] * v[704] + v[3218] * v[710] + v[2540] * v[804] + v[2600] * v[805];
	v[3381] = v[3338] + v[3220] * v[704] + v[3217] * v[710] + v[2537] * v[804] + v[2597] * v[805];
	v[3377] = v[3376] + v[3216] * v[704] + v[3213] * v[710] + v[2541] * v[806] + v[2601] * v[807];
	v[3375] = v[3354] + v[3215] * v[704] + v[3212] * v[710] + v[2540] * v[806] + v[2600] * v[807];
	v[3374] = v[3329] + v[3214] * v[704] + v[3211] * v[710] + v[2537] * v[806] + v[2597] * v[807];
	v[3370] = v[3369] + v[3210] * v[704] + v[3207] * v[710] + v[2541] * v[808] + v[2601] * v[809];
	v[3368] = v[3346] + v[3209] * v[704] + v[3206] * v[710] + v[2540] * v[808] + v[2600] * v[809];
	v[3367] = v[3320] + v[3208] * v[704] + v[3205] * v[710] + v[2537] * v[808] + v[2597] * v[809];
	v[3172] = v[1067] + v[704] * v[808] + v[710] * v[809];
	v[3165] = v[1064] + v[704] * v[806] + v[710] * v[807];
	v[3158] = v[1061] + v[704] * v[804] + v[710] * v[805];
	v[1149] = v[319] * v[334] - v[318] * v[336];
	v[1150] = -(v[319] * v[331]) + v[316] * v[336];
	v[1151] = v[318] * v[331] - v[316] * v[334];
	v[2998] = v[1004] * v[1149] + v[1151] * v[980] + v[1150] * v[992];
	v[3001] = v[1149] * v[2998] + v[2872] * v[316] + v[2871] * v[331];
	v[3000] = v[1150] * v[2998] + v[2872] * v[318] + v[2871] * v[334];
	v[2999] = v[1151] * v[2998] + v[2872] * v[319] + v[2871] * v[336];
	v[2991] = v[1002] * v[1149] + v[1151] * v[978] + v[1150] * v[990];
	v[2997] = v[1149] * v[2991] + v[2870] * v[316] + v[2868] * v[331];
	v[2995] = v[1150] * v[2991] + v[2870] * v[318] + v[2868] * v[334];
	v[2993] = v[1151] * v[2991] + v[2870] * v[319] + v[2868] * v[336];
	v[2990] = v[1000] * v[1149] + v[1151] * v[976] + v[1150] * v[988];
	v[2996] = v[1149] * v[2990] + v[2869] * v[316] + v[2867] * v[331];
	v[2994] = v[1150] * v[2990] + v[2869] * v[318] + v[2867] * v[334];
	v[2992] = v[1151] * v[2990] + v[2869] * v[319] + v[2867] * v[336];
	v[2980] = v[1151] * v[974] + v[1150] * v[986] + v[1149] * v[998];
	v[2989] = v[1149] * v[2980] + v[2866] * v[316] + v[2863] * v[331];
	v[2986] = v[1150] * v[2980] + v[2866] * v[318] + v[2863] * v[334];
	v[2983] = v[1151] * v[2980] + v[2866] * v[319] + v[2863] * v[336];
	v[2979] = v[1151] * v[972] + v[1150] * v[984] + v[1149] * v[996];
	v[2988] = v[1149] * v[2979] + v[2865] * v[316] + v[2862] * v[331];
	v[2985] = v[1150] * v[2979] + v[2865] * v[318] + v[2862] * v[334];
	v[2982] = v[1151] * v[2979] + v[2865] * v[319] + v[2862] * v[336];
	v[2978] = v[1151] * v[970] + v[1150] * v[982] + v[1149] * v[994];
	v[2987] = v[1149] * v[2978] + v[2864] * v[316] + v[2861] * v[331];
	v[2984] = v[1150] * v[2978] + v[2864] * v[318] + v[2861] * v[334];
	v[2981] = v[1151] * v[2978] + v[2864] * v[319] + v[2861] * v[336];
	v[2974] = v[1151] * v[2890] + v[1150] * v[2911] + v[1149] * v[2927];
	v[2977] = v[1149] * v[2974] + v[2929] * v[316] + v[2928] * v[331];
	v[2976] = v[1150] * v[2974] + v[2929] * v[318] + v[2928] * v[334];
	v[2975] = v[1151] * v[2974] + v[2929] * v[319] + v[2928] * v[336];
	v[2967] = v[1151] * v[2889] + v[1150] * v[2910] + v[1149] * v[2922];
	v[2973] = v[1149] * v[2967] + v[2926] * v[316] + v[2924] * v[331];
	v[2971] = v[1150] * v[2967] + v[2926] * v[318] + v[2924] * v[334];
	v[2969] = v[1151] * v[2967] + v[2926] * v[319] + v[2924] * v[336];
	v[2966] = v[1151] * v[2888] + v[1150] * v[2909] + v[1149] * v[2921];
	v[2972] = v[1149] * v[2966] + v[2925] * v[316] + v[2923] * v[331];
	v[2970] = v[1150] * v[2966] + v[2925] * v[318] + v[2923] * v[334];
	v[2968] = v[1151] * v[2966] + v[2925] * v[319] + v[2923] * v[336];
	v[2956] = v[1151] * v[2887] + v[1150] * v[2908] + v[1149] * v[2914];
	v[2965] = v[1149] * v[2956] + v[2920] * v[316] + v[2917] * v[331];
	v[2962] = v[1150] * v[2956] + v[2920] * v[318] + v[2917] * v[334];
	v[2959] = v[1151] * v[2956] + v[2920] * v[319] + v[2917] * v[336];
	v[2955] = v[1151] * v[2886] + v[1150] * v[2907] + v[1149] * v[2913];
	v[2964] = v[1149] * v[2955] + v[2919] * v[316] + v[2916] * v[331];
	v[2961] = v[1150] * v[2955] + v[2919] * v[318] + v[2916] * v[334];
	v[2958] = v[1151] * v[2955] + v[2919] * v[319] + v[2916] * v[336];
	v[2954] = v[1151] * v[2885] + v[1150] * v[2906] + v[1149] * v[2912];
	v[2963] = v[1149] * v[2954] + v[2918] * v[316] + v[2915] * v[331];
	v[2960] = v[1150] * v[2954] + v[2918] * v[318] + v[2915] * v[334];
	v[2957] = v[1151] * v[2954] + v[2918] * v[319] + v[2915] * v[336];
	v[2950] = v[1069] * v[1149] + v[1075] * v[1150] + v[1081] * v[1151];
	v[2953] = v[1149] * v[2950] + v[2884] * v[316] + v[2883] * v[331];
	v[2952] = v[1150] * v[2950] + v[2884] * v[318] + v[2883] * v[334];
	v[2951] = v[1151] * v[2950] + v[2884] * v[319] + v[2883] * v[336];
	v[2943] = v[1087] * v[1149] + v[1091] * v[1150] + v[1095] * v[1151];
	v[2949] = v[1149] * v[2943] + v[2882] * v[316] + v[2880] * v[331];
	v[2947] = v[1150] * v[2943] + v[2882] * v[318] + v[2880] * v[334];
	v[2945] = v[1151] * v[2943] + v[2882] * v[319] + v[2880] * v[336];
	v[2942] = v[1071] * v[1149] + v[1077] * v[1150] + v[1083] * v[1151];
	v[2948] = v[1149] * v[2942] + v[2881] * v[316] + v[2879] * v[331];
	v[2946] = v[1150] * v[2942] + v[2881] * v[318] + v[2879] * v[334];
	v[2944] = v[1151] * v[2942] + v[2881] * v[319] + v[2879] * v[336];
	v[2932] = v[1099] * v[1149] + v[1101] * v[1150] + v[1103] * v[1151];
	v[2941] = v[1149] * v[2932] + v[2878] * v[316] + v[2875] * v[331];
	v[2938] = v[1150] * v[2932] + v[2878] * v[318] + v[2875] * v[334];
	v[2935] = v[1151] * v[2932] + v[2878] * v[319] + v[2875] * v[336];
	v[2931] = v[1089] * v[1149] + v[1093] * v[1150] + v[1097] * v[1151];
	v[2940] = v[1149] * v[2931] + v[2877] * v[316] + v[2874] * v[331];
	v[2937] = v[1150] * v[2931] + v[2877] * v[318] + v[2874] * v[334];
	v[2934] = v[1151] * v[2931] + v[2877] * v[319] + v[2874] * v[336];
	v[2930] = v[1073] * v[1149] + v[1079] * v[1150] + v[1085] * v[1151];
	v[2939] = v[1149] * v[2930] + v[2876] * v[316] + v[2873] * v[331];
	v[2936] = v[1150] * v[2930] + v[2876] * v[318] + v[2873] * v[334];
	v[2933] = v[1151] * v[2930] + v[2876] * v[319] + v[2873] * v[336];
	v[2664] = v[1151] * v[720] + v[1150] * v[732] + v[1149] * v[744];
	v[2673] = v[1149] * v[2664] + v[1926] * v[316] + v[1923] * v[331];
	v[2670] = v[1150] * v[2664] + v[1926] * v[318] + v[1923] * v[334];
	v[2667] = v[1151] * v[2664] + v[1926] * v[319] + v[1923] * v[336];
	v[2663] = v[1151] * v[724] + v[1150] * v[736] + v[1149] * v[748];
	v[2672] = v[1149] * v[2663] + v[1925] * v[316] + v[1922] * v[331];
	v[2669] = v[1150] * v[2663] + v[1925] * v[318] + v[1922] * v[334];
	v[2666] = v[1151] * v[2663] + v[1925] * v[319] + v[1922] * v[336];
	v[2662] = v[1151] * v[728] + v[1150] * v[740] + v[1149] * v[752];
	v[2671] = v[1149] * v[2662] + v[1924] * v[316] + v[1921] * v[331];
	v[2668] = v[1150] * v[2662] + v[1924] * v[318] + v[1921] * v[334];
	v[2665] = v[1151] * v[2662] + v[1924] * v[319] + v[1921] * v[336];
	v[2652] = v[1151] * v[722] + v[1150] * v[734] + v[1149] * v[746];
	v[2661] = v[1149] * v[2652] + v[1920] * v[316] + v[1917] * v[331];
	v[2658] = v[1150] * v[2652] + v[1920] * v[318] + v[1917] * v[334];
	v[2655] = v[1151] * v[2652] + v[1920] * v[319] + v[1917] * v[336];
	v[2651] = v[1151] * v[726] + v[1150] * v[738] + v[1149] * v[750];
	v[2660] = v[1149] * v[2651] + v[1919] * v[316] + v[1916] * v[331];
	v[2657] = v[1150] * v[2651] + v[1919] * v[318] + v[1916] * v[334];
	v[2654] = v[1151] * v[2651] + v[1919] * v[319] + v[1916] * v[336];
	v[2650] = v[1151] * v[730] + v[1150] * v[742] + v[1149] * v[754];
	v[2659] = v[1149] * v[2650] + v[1918] * v[316] + v[1915] * v[331];
	v[2656] = v[1150] * v[2650] + v[1918] * v[318] + v[1915] * v[334];
	v[2653] = v[1151] * v[2650] + v[1918] * v[319] + v[1915] * v[336];
	v[2640] = v[1151] * v[1971] + v[1150] * v[1977] + v[1149] * v[1989];
	v[2649] = v[1149] * v[2640] + v[1995] * v[316] + v[1992] * v[331];
	v[2646] = v[1150] * v[2640] + v[1995] * v[318] + v[1992] * v[334];
	v[2643] = v[1151] * v[2640] + v[1995] * v[319] + v[1992] * v[336];
	v[2639] = v[1151] * v[1970] + v[1150] * v[1976] + v[1149] * v[1988];
	v[2648] = v[1149] * v[2639] + v[1994] * v[316] + v[1991] * v[331];
	v[2645] = v[1150] * v[2639] + v[1994] * v[318] + v[1991] * v[334];
	v[2642] = v[1151] * v[2639] + v[1994] * v[319] + v[1991] * v[336];
	v[2638] = v[1151] * v[1969] + v[1150] * v[1975] + v[1149] * v[1987];
	v[2647] = v[1149] * v[2638] + v[1993] * v[316] + v[1990] * v[331];
	v[2644] = v[1150] * v[2638] + v[1993] * v[318] + v[1990] * v[334];
	v[2641] = v[1151] * v[2638] + v[1993] * v[319] + v[1990] * v[336];
	v[2628] = v[1151] * v[1968] + v[1150] * v[1974] + v[1149] * v[1980];
	v[2637] = v[1149] * v[2628] + v[1986] * v[316] + v[1983] * v[331];
	v[2634] = v[1150] * v[2628] + v[1986] * v[318] + v[1983] * v[334];
	v[2631] = v[1151] * v[2628] + v[1986] * v[319] + v[1983] * v[336];
	v[2627] = v[1151] * v[1967] + v[1150] * v[1973] + v[1149] * v[1979];
	v[2636] = v[1149] * v[2627] + v[1985] * v[316] + v[1982] * v[331];
	v[2633] = v[1150] * v[2627] + v[1985] * v[318] + v[1982] * v[334];
	v[2630] = v[1151] * v[2627] + v[1985] * v[319] + v[1982] * v[336];
	v[2626] = v[1151] * v[1966] + v[1150] * v[1972] + v[1149] * v[1978];
	v[2635] = v[1149] * v[2626] + v[1984] * v[316] + v[1981] * v[331];
	v[2632] = v[1150] * v[2626] + v[1984] * v[318] + v[1981] * v[334];
	v[2629] = v[1151] * v[2626] + v[1984] * v[319] + v[1981] * v[336];
	v[2616] = v[1149] * v[822] + v[1150] * v[824] + v[1151] * v[826];
	v[2625] = v[1149] * v[2616] + v[1965] * v[316] + v[1962] * v[331];
	v[2622] = v[1150] * v[2616] + v[1965] * v[318] + v[1962] * v[334];
	v[2619] = v[1151] * v[2616] + v[1965] * v[319] + v[1962] * v[336];
	v[2615] = v[1149] * v[816] + v[1150] * v[818] + v[1151] * v[820];
	v[2624] = v[1149] * v[2615] + v[1964] * v[316] + v[1961] * v[331];
	v[2621] = v[1150] * v[2615] + v[1964] * v[318] + v[1961] * v[334];
	v[2618] = v[1151] * v[2615] + v[1964] * v[319] + v[1961] * v[336];
	v[2614] = v[1149] * v[810] + v[1150] * v[812] + v[1151] * v[814];
	v[2623] = v[1149] * v[2614] + v[1963] * v[316] + v[1960] * v[331];
	v[2620] = v[1150] * v[2614] + v[1963] * v[318] + v[1960] * v[334];
	v[2617] = v[1151] * v[2614] + v[1963] * v[319] + v[1960] * v[336];
	v[2604] = v[1149] * v[823] + v[1150] * v[825] + v[1151] * v[827];
	v[2613] = v[1149] * v[2604] + v[1959] * v[316] + v[1956] * v[331];
	v[2610] = v[1150] * v[2604] + v[1959] * v[318] + v[1956] * v[334];
	v[2607] = v[1151] * v[2604] + v[1959] * v[319] + v[1956] * v[336];
	v[2603] = v[1149] * v[817] + v[1150] * v[819] + v[1151] * v[821];
	v[2612] = v[1149] * v[2603] + v[1958] * v[316] + v[1955] * v[331];
	v[2609] = v[1150] * v[2603] + v[1958] * v[318] + v[1955] * v[334];
	v[2606] = v[1151] * v[2603] + v[1958] * v[319] + v[1955] * v[336];
	v[2602] = v[1149] * v[811] + v[1150] * v[813] + v[1151] * v[815];
	v[2611] = v[1149] * v[2602] + v[1957] * v[316] + v[1954] * v[331];
	v[2608] = v[1150] * v[2602] + v[1957] * v[318] + v[1954] * v[334];
	v[2605] = v[1151] * v[2602] + v[1957] * v[319] + v[1954] * v[336];
	v[1856] = v[1149] * v[1518] + v[1150] * v[1520] + v[1151] * v[1522];
	v[1859] = v[1149] * v[1856] + v[1524] * v[316] + v[1523] * v[331];
	v[1858] = v[1150] * v[1856] + v[1524] * v[318] + v[1523] * v[334];
	v[1857] = v[1151] * v[1856] + v[1524] * v[319] + v[1523] * v[336];
	v[1849] = v[1149] * v[1504] + v[1150] * v[1508] + v[1151] * v[1512];
	v[1855] = v[1149] * v[1849] + v[1516] * v[316] + v[1514] * v[331];
	v[1853] = v[1150] * v[1849] + v[1516] * v[318] + v[1514] * v[334];
	v[1851] = v[1151] * v[1849] + v[1516] * v[319] + v[1514] * v[336];
	v[1848] = v[1149] * v[1502] + v[1150] * v[1506] + v[1151] * v[1510];
	v[1854] = v[1149] * v[1848] + v[1515] * v[316] + v[1513] * v[331];
	v[1852] = v[1150] * v[1848] + v[1515] * v[318] + v[1513] * v[334];
	v[1850] = v[1151] * v[1848] + v[1515] * v[319] + v[1513] * v[336];
	v[1844] = v[1151] * v[1594] + v[1150] * v[1597] + v[1149] * v[1604];
	v[1847] = v[1149] * v[1844] + v[1606] * v[316] + v[1605] * v[331];
	v[1846] = v[1150] * v[1844] + v[1606] * v[318] + v[1605] * v[334];
	v[1845] = v[1151] * v[1844] + v[1606] * v[319] + v[1605] * v[336];
	v[1837] = v[1151] * v[1593] + v[1150] * v[1596] + v[1149] * v[1599];
	v[1843] = v[1149] * v[1837] + v[1603] * v[316] + v[1601] * v[331];
	v[1841] = v[1150] * v[1837] + v[1603] * v[318] + v[1601] * v[334];
	v[1839] = v[1151] * v[1837] + v[1603] * v[319] + v[1601] * v[336];
	v[1836] = v[1151] * v[1592] + v[1150] * v[1595] + v[1149] * v[1598];
	v[1842] = v[1149] * v[1836] + v[1602] * v[316] + v[1600] * v[331];
	v[1840] = v[1150] * v[1836] + v[1602] * v[318] + v[1600] * v[334];
	v[1838] = v[1151] * v[1836] + v[1602] * v[319] + v[1600] * v[336];
	v[1832] = v[1151] * v[1575] + v[1150] * v[1580] + v[1149] * v[1589];
	v[1835] = v[1149] * v[1832] + v[1591] * v[316] + v[1590] * v[331];
	v[1834] = v[1150] * v[1832] + v[1591] * v[318] + v[1590] * v[334];
	v[1833] = v[1151] * v[1832] + v[1591] * v[319] + v[1590] * v[336];
	v[1825] = v[1151] * v[1573] + v[1150] * v[1578] + v[1149] * v[1583];
	v[1831] = v[1149] * v[1825] + v[1587] * v[316] + v[1585] * v[331];
	v[1829] = v[1150] * v[1825] + v[1587] * v[318] + v[1585] * v[334];
	v[1827] = v[1151] * v[1825] + v[1587] * v[319] + v[1585] * v[336];
	v[1824] = v[1151] * v[1571] + v[1150] * v[1576] + v[1149] * v[1581];
	v[1830] = v[1149] * v[1824] + v[1586] * v[316] + v[1584] * v[331];
	v[1828] = v[1150] * v[1824] + v[1586] * v[318] + v[1584] * v[334];
	v[1826] = v[1151] * v[1824] + v[1586] * v[319] + v[1584] * v[336];
	v[1367] = v[1149] * v[545] + v[1150] * v[547] + v[1151] * v[549];
	v[1393] = v[1151] * v[1367] + v[1361] * v[319] + v[1373] * v[336];
	v[1387] = v[1150] * v[1367] + v[1361] * v[318] + v[1373] * v[334];
	v[1381] = v[1149] * v[1367] + v[1361] * v[316] + v[1373] * v[331];
	v[1366] = v[1149] * v[544] + v[1150] * v[546] + v[1151] * v[548];
	v[1392] = v[1151] * v[1366] + v[1360] * v[319] + v[1372] * v[336];
	v[1386] = v[1150] * v[1366] + v[1360] * v[318] + v[1372] * v[334];
	v[1380] = v[1149] * v[1366] + v[1360] * v[316] + v[1372] * v[331];
	v[1365] = v[1149] * v[1351] + v[1150] * v[1353] + v[1151] * v[1355];
	v[1391] = v[1151] * v[1365] + v[1359] * v[319] + v[1371] * v[336];
	v[1385] = v[1150] * v[1365] + v[1359] * v[318] + v[1371] * v[334];
	v[1379] = v[1149] * v[1365] + v[1359] * v[316] + v[1371] * v[331];
	v[1364] = v[1149] * v[1350] + v[1150] * v[1352] + v[1151] * v[1354];
	v[1390] = v[1151] * v[1364] + v[1358] * v[319] + v[1370] * v[336];
	v[1384] = v[1150] * v[1364] + v[1358] * v[318] + v[1370] * v[334];
	v[1378] = v[1149] * v[1364] + v[1358] * v[316] + v[1370] * v[331];
	v[1363] = v[1151] * v[407] + v[1150] * v[408] + v[1149] * v[409];
	v[1389] = v[1151] * v[1363] + v[1357] * v[319] + v[1369] * v[336];
	v[1383] = v[1150] * v[1363] + v[1357] * v[318] + v[1369] * v[334];
	v[1377] = v[1149] * v[1363] + v[1357] * v[316] + v[1369] * v[331];
	v[1362] = v[1151] * v[404] + v[1150] * v[405] + v[1149] * v[406];
	v[1388] = v[1151] * v[1362] + v[1356] * v[319] + v[1368] * v[336];
	v[1382] = v[1150] * v[1362] + v[1356] * v[318] + v[1368] * v[334];
	v[1376] = v[1149] * v[1362] + v[1356] * v[316] + v[1368] * v[331];
	v[1223] = v[1149] * v[675] + v[1150] * v[681] + v[1151] * v[687];
	v[1325] = v[1151] * v[1223] + v[1214] * v[319] + v[1232] * v[336];
	v[1316] = v[1150] * v[1223] + v[1214] * v[318] + v[1232] * v[334];
	v[1307] = v[1149] * v[1223] + v[1214] * v[316] + v[1232] * v[331];
	v[1222] = v[1149] * v[673] + v[1150] * v[679] + v[1151] * v[685];
	v[1324] = v[1151] * v[1222] + v[1213] * v[319] + v[1231] * v[336];
	v[1315] = v[1150] * v[1222] + v[1213] * v[318] + v[1231] * v[334];
	v[1306] = v[1149] * v[1222] + v[1213] * v[316] + v[1231] * v[331];
	v[1221] = v[1149] * v[671] + v[1150] * v[677] + v[1151] * v[683];
	v[1323] = v[1151] * v[1221] + v[1212] * v[319] + v[1230] * v[336];
	v[1314] = v[1150] * v[1221] + v[1212] * v[318] + v[1230] * v[334];
	v[1305] = v[1149] * v[1221] + v[1212] * v[316] + v[1230] * v[331];
	v[1220] = v[1149] * v[1199] + v[1150] * v[1202] + v[1151] * v[1205];
	v[1322] = v[1151] * v[1220] + v[1211] * v[319] + v[1229] * v[336];
	v[1313] = v[1150] * v[1220] + v[1211] * v[318] + v[1229] * v[334];
	v[1304] = v[1149] * v[1220] + v[1211] * v[316] + v[1229] * v[331];
	v[1219] = v[1149] * v[1198] + v[1150] * v[1201] + v[1151] * v[1204];
	v[1321] = v[1151] * v[1219] + v[1210] * v[319] + v[1228] * v[336];
	v[1312] = v[1150] * v[1219] + v[1210] * v[318] + v[1228] * v[334];
	v[1303] = v[1149] * v[1219] + v[1210] * v[316] + v[1228] * v[331];
	v[1218] = v[1149] * v[1197] + v[1150] * v[1200] + v[1151] * v[1203];
	v[1320] = v[1151] * v[1218] + v[1209] * v[319] + v[1227] * v[336];
	v[1311] = v[1150] * v[1218] + v[1209] * v[318] + v[1227] * v[334];
	v[1302] = v[1149] * v[1218] + v[1209] * v[316] + v[1227] * v[331];
	v[1217] = v[1149] * v[637] + v[1150] * v[640] + v[1151] * v[643];
	v[1319] = v[1151] * v[1217] + v[1208] * v[319] + v[1226] * v[336];
	v[1310] = v[1150] * v[1217] + v[1208] * v[318] + v[1226] * v[334];
	v[1301] = v[1149] * v[1217] + v[1208] * v[316] + v[1226] * v[331];
	v[1216] = v[1149] * v[636] + v[1150] * v[639] + v[1151] * v[642];
	v[1318] = v[1151] * v[1216] + v[1207] * v[319] + v[1225] * v[336];
	v[1309] = v[1150] * v[1216] + v[1207] * v[318] + v[1225] * v[334];
	v[1300] = v[1149] * v[1216] + v[1207] * v[316] + v[1225] * v[331];
	v[1215] = v[1149] * v[635] + v[1150] * v[638] + v[1151] * v[641];
	v[1317] = v[1151] * v[1215] + v[1206] * v[319] + v[1224] * v[336];
	v[1308] = v[1150] * v[1215] + v[1206] * v[318] + v[1224] * v[334];
	v[1299] = v[1149] * v[1215] + v[1206] * v[316] + v[1224] * v[331];
	v[1152] = v[346] * v[361] - v[345] * v[363];
	v[1153] = -(v[346] * v[358]) + v[343] * v[363];
	v[1154] = v[345] * v[358] - v[343] * v[361];
	v[1155] = v[316] * v[343] + v[318] * v[345] + v[319] * v[346];
	v[1156] = v[1152] * v[316] + v[1153] * v[318] + v[1154] * v[319];
	v[1157] = v[316] * v[358] + v[318] * v[361] + v[319] * v[363];
	v[1158] = v[1149] * v[343] + v[1150] * v[345] + v[1151] * v[346];
	v[1159] = v[1149] * v[1152] + v[1150] * v[1153] + v[1151] * v[1154];
	v[1160] = v[1149] * v[358] + v[1150] * v[361] + v[1151] * v[363];
	v[1161] = v[331] * v[343] + v[334] * v[345] + v[336] * v[346];
	v[1162] = v[1152] * v[331] + v[1153] * v[334] + v[1154] * v[336];
	v[1163] = v[331] * v[358] + v[334] * v[361] + v[336] * v[363];
	v[1165] = d[0] + ci[0] * (e1i[0] - v[281]) + ci[1] * (e2i[0] - v[284]) + (*radius)*(-(v[1164] * v[331]) - v[266] * v[334]
		- v[267] * v[336]) - v[4515] + v[296] * (-v[287] + v[4534]) + xPi[0];
	v[1167] = d[1] + ci[0] * (e1i[1] - v[282]) + ci[1] * (e2i[1] - v[285]) + (*radius)*(-(v[269] * v[331]) - v[1166] * v[334]
		- v[272] * v[336]) - v[4516] + v[296] * (-v[288] + v[4535]) + xPi[1];
	v[1169] = d[2] + ci[0] * (e1i[2] - v[283]) + ci[1] * (e2i[2] - v[286]) + (*radius)*(-(v[274] * v[331]) - v[276] * v[334]
		- v[1168] * v[336]) - v[4517] + v[296] * (-v[289] + v[4536]) + xPi[2];
	v[3007] = v[1069] * v[1165] + v[1075] * v[1167] + v[1081] * v[1169] - v[2779] * v[358] - v[2773] * v[361]
		- v[2767] * v[363] + v[1290] * v[4755] + v[1270] * v[4758] + v[1251] * v[4761];
	v[3006] = v[1087] * v[1165] + v[1091] * v[1167] + v[1095] * v[1169] - v[2778] * v[358] - v[2772] * v[361]
		- v[2766] * v[363] + v[1291] * v[4756] + v[1271] * v[4759] + v[1252] * v[4762];
	v[3005] = v[1071] * v[1165] + v[1077] * v[1167] + v[1083] * v[1169] - v[2777] * v[358] - v[2771] * v[361]
		- v[2765] * v[363] - v[1252] * v[671] - v[1251] * v[673] - v[1271] * v[677] - v[1270] * v[679] - v[1291] * v[683]
		- v[1290] * v[685];
	v[3004] = v[1099] * v[1165] + v[1101] * v[1167] + v[1103] * v[1169] - v[2776] * v[358] - v[2770] * v[361]
		- v[2764] * v[363] + v[1292] * v[4757] + v[1272] * v[4760] + v[1253] * v[4763];
	v[3003] = v[1089] * v[1165] + v[1093] * v[1167] + v[1097] * v[1169] - v[2775] * v[358] - v[2769] * v[361]
		- v[2763] * v[363] - v[1253] * v[673] - v[1252] * v[675] - v[1272] * v[679] - v[1271] * v[681] - v[1292] * v[685]
		- v[1291] * v[687];
	v[3002] = v[1073] * v[1165] + v[1079] * v[1167] + v[1085] * v[1169] - v[2774] * v[358] - v[2768] * v[361]
		- v[2762] * v[363] - v[1253] * v[671] - v[1251] * v[675] - v[1272] * v[677] - v[1270] * v[681] - v[1292] * v[683]
		- v[1290] * v[687];
	v[2685] = -(v[1253] * v[544]) - v[1272] * v[546] - v[1292] * v[548] + v[1165] * v[822] + v[1167] * v[824]
		+ v[1169] * v[826];
	v[2684] = -(v[1252] * v[544]) - v[1271] * v[546] - v[1291] * v[548] + v[1165] * v[816] + v[1167] * v[818]
		+ v[1169] * v[820];
	v[2683] = -(v[1251] * v[544]) - v[1270] * v[546] - v[1290] * v[548] + v[1165] * v[810] + v[1167] * v[812]
		+ v[1169] * v[814];
	v[2679] = -(v[1253] * v[545]) - v[1272] * v[547] - v[1292] * v[549] + v[1165] * v[823] + v[1167] * v[825]
		+ v[1169] * v[827];
	v[2678] = -(v[1252] * v[545]) - v[1271] * v[547] - v[1291] * v[549] + v[1165] * v[817] + v[1167] * v[819]
		+ v[1169] * v[821];
	v[2677] = -(v[1251] * v[545]) - v[1270] * v[547] - v[1290] * v[549] + v[1165] * v[811] + v[1167] * v[813]
		+ v[1169] * v[815];
	v[1862] = v[1169] * v[1575] + v[1167] * v[1580] + v[1165] * v[1589];
	v[1861] = v[1169] * v[1573] + v[1167] * v[1578] + v[1165] * v[1583];
	v[1860] = v[1169] * v[1571] + v[1167] * v[1576] + v[1165] * v[1581];
	v[1375] = v[1165] * v[545] + v[1167] * v[547] + v[1169] * v[549];
	v[1374] = v[1165] * v[544] + v[1167] * v[546] + v[1169] * v[548];
	v[1298] = -(v[1253] * v[358]) - v[1272] * v[361] - v[1292] * v[363] + v[1165] * v[675] + v[1167] * v[681]
		+ v[1169] * v[687];
	v[1297] = -(v[1252] * v[358]) - v[1271] * v[361] - v[1291] * v[363] + v[1165] * v[673] + v[1167] * v[679]
		+ v[1169] * v[685];
	v[1296] = -(v[1251] * v[358]) - v[1270] * v[361] - v[1290] * v[363] + v[1165] * v[671] + v[1167] * v[677]
		+ v[1169] * v[683];
	v[1177] = v[1165] * v[358] + v[1167] * v[361] + v[1169] * v[363];
	v[3955] = -(v[1081] * v[1177]) - v[2767] + gti[0] * (v[1149] * v[2975] + v[2999] * v[316] + v[2951] * v[331]) + gti[1] *
		(v[1150] * v[2975] + v[2999] * v[318] + v[2951] * v[334]) + gti[2] * (v[1151] * v[2975] + v[2999] * v[319]
			+ v[2951] * v[336]) - v[3007] * v[363] + v[1296] * v[4755];
	v[3942] = -(v[1095] * v[1177]) - v[2766] + gti[0] * (v[1149] * v[2969] + v[2993] * v[316] + v[2945] * v[331]) + gti[1] *
		(v[1150] * v[2969] + v[2993] * v[318] + v[2945] * v[334]) + gti[2] * (v[1151] * v[2969] + v[2993] * v[319]
			+ v[2945] * v[336]) - v[3006] * v[363] + v[1297] * v[4756];
	v[3940] = -(v[1083] * v[1177]) - v[2765] + gti[0] * (v[1149] * v[2968] + v[2992] * v[316] + v[2944] * v[331]) + gti[1] *
		(v[1150] * v[2968] + v[2992] * v[318] + v[2944] * v[334]) + gti[2] * (v[1151] * v[2968] + v[2992] * v[319]
			+ v[2944] * v[336]) - v[3005] * v[363] - v[1297] * v[683] - v[1296] * v[685];
	v[3929] = -(v[1103] * v[1177]) - v[2764] + gti[0] * (v[1149] * v[2959] + v[2983] * v[316] + v[2935] * v[331]) + gti[1] *
		(v[1150] * v[2959] + v[2983] * v[318] + v[2935] * v[334]) + gti[2] * (v[1151] * v[2959] + v[2983] * v[319]
			+ v[2935] * v[336]) - v[3004] * v[363] + v[1298] * v[4757];
	v[3927] = -(v[1097] * v[1177]) - v[2763] + gti[0] * (v[1149] * v[2958] + v[2982] * v[316] + v[2934] * v[331]) + gti[1] *
		(v[1150] * v[2958] + v[2982] * v[318] + v[2934] * v[334]) + gti[2] * (v[1151] * v[2958] + v[2982] * v[319]
			+ v[2934] * v[336]) - v[3003] * v[363] - v[1298] * v[685] - v[1297] * v[687];
	v[3925] = -(v[1085] * v[1177]) - v[2762] + gti[0] * (v[1149] * v[2957] + v[2981] * v[316] + v[2933] * v[331]) + gti[1] *
		(v[1150] * v[2957] + v[2981] * v[318] + v[2933] * v[334]) + gti[2] * (v[1151] * v[2957] + v[2981] * v[319]
			+ v[2933] * v[336]) - v[3002] * v[363] - v[1298] * v[683] - v[1296] * v[687];
	v[3870] = -(v[1075] * v[1177]) - v[2773] + gti[0] * (v[1149] * v[2976] + v[3000] * v[316] + v[2952] * v[331]) + gti[1] *
		(v[1150] * v[2976] + v[3000] * v[318] + v[2952] * v[334]) + gti[2] * (v[1151] * v[2976] + v[3000] * v[319]
			+ v[2952] * v[336]) - v[3007] * v[361] + v[1296] * v[4758];
	v[3858] = -(v[1091] * v[1177]) - v[2772] + gti[0] * (v[1149] * v[2971] + v[2995] * v[316] + v[2947] * v[331]) + gti[1] *
		(v[1150] * v[2971] + v[2995] * v[318] + v[2947] * v[334]) + gti[2] * (v[1151] * v[2971] + v[2995] * v[319]
			+ v[2947] * v[336]) - v[3006] * v[361] + v[1297] * v[4759];
	v[3856] = -(v[1077] * v[1177]) - v[2771] + gti[0] * (v[1149] * v[2970] + v[2994] * v[316] + v[2946] * v[331]) + gti[1] *
		(v[1150] * v[2970] + v[2994] * v[318] + v[2946] * v[334]) + gti[2] * (v[1151] * v[2970] + v[2994] * v[319]
			+ v[2946] * v[336]) - v[3005] * v[361] - v[1297] * v[677] - v[1296] * v[679];
	v[3846] = -(v[1101] * v[1177]) - v[2770] + gti[0] * (v[1149] * v[2962] + v[2986] * v[316] + v[2938] * v[331]) + gti[1] *
		(v[1150] * v[2962] + v[2986] * v[318] + v[2938] * v[334]) + gti[2] * (v[1151] * v[2962] + v[2986] * v[319]
			+ v[2938] * v[336]) - v[3004] * v[361] + v[1298] * v[4760];
	v[3844] = -(v[1093] * v[1177]) - v[2769] + gti[0] * (v[1149] * v[2961] + v[2985] * v[316] + v[2937] * v[331]) + gti[1] *
		(v[1150] * v[2961] + v[2985] * v[318] + v[2937] * v[334]) + gti[2] * (v[1151] * v[2961] + v[2985] * v[319]
			+ v[2937] * v[336]) - v[3003] * v[361] - v[1298] * v[679] - v[1297] * v[681];
	v[3842] = -(v[1079] * v[1177]) - v[2768] + gti[0] * (v[1149] * v[2960] + v[2984] * v[316] + v[2936] * v[331]) + gti[1] *
		(v[1150] * v[2960] + v[2984] * v[318] + v[2936] * v[334]) + gti[2] * (v[1151] * v[2960] + v[2984] * v[319]
			+ v[2936] * v[336]) - v[3002] * v[361] - v[1298] * v[677] - v[1296] * v[681];
	v[3792] = -(v[1069] * v[1177]) - v[2779] + gti[0] * (v[1149] * v[2977] + v[3001] * v[316] + v[2953] * v[331]) + gti[1] *
		(v[1150] * v[2977] + v[3001] * v[318] + v[2953] * v[334]) + gti[2] * (v[1151] * v[2977] + v[3001] * v[319]
			+ v[2953] * v[336]) - v[3007] * v[358] + v[1296] * v[4761];
	v[3781] = -(v[1087] * v[1177]) - v[2778] + gti[0] * (v[1149] * v[2973] + v[2997] * v[316] + v[2949] * v[331]) + gti[1] *
		(v[1150] * v[2973] + v[2997] * v[318] + v[2949] * v[334]) + gti[2] * (v[1151] * v[2973] + v[2997] * v[319]
			+ v[2949] * v[336]) - v[3006] * v[358] + v[1297] * v[4762];
	v[3779] = -(v[1071] * v[1177]) - v[2777] + gti[0] * (v[1149] * v[2972] + v[2996] * v[316] + v[2948] * v[331]) + gti[1] *
		(v[1150] * v[2972] + v[2996] * v[318] + v[2948] * v[334]) + gti[2] * (v[1151] * v[2972] + v[2996] * v[319]
			+ v[2948] * v[336]) - v[3005] * v[358] - v[1297] * v[671] - v[1296] * v[673];
	v[3770] = -(v[1099] * v[1177]) - v[2776] + gti[0] * (v[1149] * v[2965] + v[2989] * v[316] + v[2941] * v[331]) + gti[1] *
		(v[1150] * v[2965] + v[2989] * v[318] + v[2941] * v[334]) + gti[2] * (v[1151] * v[2965] + v[2989] * v[319]
			+ v[2941] * v[336]) - v[3004] * v[358] + v[1298] * v[4763];
	v[3768] = -(v[1089] * v[1177]) - v[2775] + gti[0] * (v[1149] * v[2964] + v[2988] * v[316] + v[2940] * v[331]) + gti[1] *
		(v[1150] * v[2964] + v[2988] * v[318] + v[2940] * v[334]) + gti[2] * (v[1151] * v[2964] + v[2988] * v[319]
			+ v[2940] * v[336]) - v[3003] * v[358] - v[1298] * v[673] - v[1297] * v[675];
	v[3766] = -(v[1073] * v[1177]) - v[2774] + gti[0] * (v[1149] * v[2963] + v[2987] * v[316] + v[2939] * v[331]) + gti[1] *
		(v[1150] * v[2963] + v[2987] * v[318] + v[2939] * v[334]) + gti[2] * (v[1151] * v[2963] + v[2987] * v[319]
			+ v[2939] * v[336]) - v[3002] * v[358] - v[1298] * v[671] - v[1296] * v[675];
	v[2733] = gti[0] * (v[1149] * v[2631] + v[2655] * v[316] + v[2607] * v[331]) + gti[1] * (v[1150] * v[2631]
		+ v[2655] * v[318] + v[2607] * v[334]) + gti[2] * (v[1151] * v[2631] + v[2655] * v[319] + v[2607] * v[336])
		- v[2679] * v[363] - v[1298] * v[549] - v[1375] * v[687] - v[1177] * v[827];
	v[2732] = gti[0] * (v[1149] * v[2630] + v[2654] * v[316] + v[2606] * v[331]) + gti[1] * (v[1150] * v[2630]
		+ v[2654] * v[318] + v[2606] * v[334]) + gti[2] * (v[1151] * v[2630] + v[2654] * v[319] + v[2606] * v[336])
		- v[2678] * v[363] - v[1297] * v[549] - v[1375] * v[685] - v[1177] * v[821];
	v[2731] = gti[0] * (v[1149] * v[2629] + v[2653] * v[316] + v[2605] * v[331]) + gti[1] * (v[1150] * v[2629]
		+ v[2653] * v[318] + v[2605] * v[334]) + gti[2] * (v[1151] * v[2629] + v[2653] * v[319] + v[2605] * v[336])
		- v[2677] * v[363] - v[1296] * v[549] - v[1375] * v[683] - v[1177] * v[815];
	v[2726] = gti[0] * (v[1149] * v[2634] + v[2658] * v[316] + v[2610] * v[331]) + gti[1] * (v[1150] * v[2634]
		+ v[2658] * v[318] + v[2610] * v[334]) + gti[2] * (v[1151] * v[2634] + v[2658] * v[319] + v[2610] * v[336])
		- v[2679] * v[361] - v[1298] * v[547] - v[1375] * v[681] - v[1177] * v[825];
	v[2725] = gti[0] * (v[1149] * v[2633] + v[2657] * v[316] + v[2609] * v[331]) + gti[1] * (v[1150] * v[2633]
		+ v[2657] * v[318] + v[2609] * v[334]) + gti[2] * (v[1151] * v[2633] + v[2657] * v[319] + v[2609] * v[336])
		- v[2678] * v[361] - v[1297] * v[547] - v[1375] * v[679] - v[1177] * v[819];
	v[2724] = gti[0] * (v[1149] * v[2632] + v[2656] * v[316] + v[2608] * v[331]) + gti[1] * (v[1150] * v[2632]
		+ v[2656] * v[318] + v[2608] * v[334]) + gti[2] * (v[1151] * v[2632] + v[2656] * v[319] + v[2608] * v[336])
		- v[2677] * v[361] - v[1296] * v[547] - v[1375] * v[677] - v[1177] * v[813];
	v[2718] = gti[0] * (v[1149] * v[2637] + v[2661] * v[316] + v[2613] * v[331]) + gti[1] * (v[1150] * v[2637]
		+ v[2661] * v[318] + v[2613] * v[334]) + gti[2] * (v[1151] * v[2637] + v[2661] * v[319] + v[2613] * v[336])
		- v[2679] * v[358] - v[1298] * v[545] - v[1375] * v[675] - v[1177] * v[823];
	v[2717] = gti[0] * (v[1149] * v[2636] + v[2660] * v[316] + v[2612] * v[331]) + gti[1] * (v[1150] * v[2636]
		+ v[2660] * v[318] + v[2612] * v[334]) + gti[2] * (v[1151] * v[2636] + v[2660] * v[319] + v[2612] * v[336])
		- v[2678] * v[358] - v[1297] * v[545] - v[1375] * v[673] - v[1177] * v[817];
	v[2716] = gti[0] * (v[1149] * v[2635] + v[2659] * v[316] + v[2611] * v[331]) + gti[1] * (v[1150] * v[2635]
		+ v[2659] * v[318] + v[2611] * v[334]) + gti[2] * (v[1151] * v[2635] + v[2659] * v[319] + v[2611] * v[336])
		- v[2677] * v[358] - v[1296] * v[545] - v[1375] * v[671] - v[1177] * v[811];
	v[2709] = gti[0] * (v[1149] * v[2643] + v[2667] * v[316] + v[2619] * v[331]) + gti[1] * (v[1150] * v[2643]
		+ v[2667] * v[318] + v[2619] * v[334]) + gti[2] * (v[1151] * v[2643] + v[2667] * v[319] + v[2619] * v[336])
		- v[2685] * v[363] - v[1298] * v[548] - v[1374] * v[687] - v[1177] * v[826];
	v[2708] = gti[0] * (v[1149] * v[2642] + v[2666] * v[316] + v[2618] * v[331]) + gti[1] * (v[1150] * v[2642]
		+ v[2666] * v[318] + v[2618] * v[334]) + gti[2] * (v[1151] * v[2642] + v[2666] * v[319] + v[2618] * v[336])
		- v[2684] * v[363] - v[1297] * v[548] - v[1374] * v[685] - v[1177] * v[820];
	v[2707] = gti[0] * (v[1149] * v[2641] + v[2665] * v[316] + v[2617] * v[331]) + gti[1] * (v[1150] * v[2641]
		+ v[2665] * v[318] + v[2617] * v[334]) + gti[2] * (v[1151] * v[2641] + v[2665] * v[319] + v[2617] * v[336])
		- v[2683] * v[363] - v[1296] * v[548] - v[1374] * v[683] - v[1177] * v[814];
	v[2702] = gti[0] * (v[1149] * v[2646] + v[2670] * v[316] + v[2622] * v[331]) + gti[1] * (v[1150] * v[2646]
		+ v[2670] * v[318] + v[2622] * v[334]) + gti[2] * (v[1151] * v[2646] + v[2670] * v[319] + v[2622] * v[336])
		- v[2685] * v[361] - v[1298] * v[546] - v[1374] * v[681] - v[1177] * v[824];
	v[2701] = gti[0] * (v[1149] * v[2645] + v[2669] * v[316] + v[2621] * v[331]) + gti[1] * (v[1150] * v[2645]
		+ v[2669] * v[318] + v[2621] * v[334]) + gti[2] * (v[1151] * v[2645] + v[2669] * v[319] + v[2621] * v[336])
		- v[2684] * v[361] - v[1297] * v[546] - v[1374] * v[679] - v[1177] * v[818];
	v[2700] = gti[0] * (v[1149] * v[2644] + v[2668] * v[316] + v[2620] * v[331]) + gti[1] * (v[1150] * v[2644]
		+ v[2668] * v[318] + v[2620] * v[334]) + gti[2] * (v[1151] * v[2644] + v[2668] * v[319] + v[2620] * v[336])
		- v[2683] * v[361] - v[1296] * v[546] - v[1374] * v[677] - v[1177] * v[812];
	v[2694] = gti[0] * (v[1149] * v[2649] + v[2673] * v[316] + v[2625] * v[331]) + gti[1] * (v[1150] * v[2649]
		+ v[2673] * v[318] + v[2625] * v[334]) + gti[2] * (v[1151] * v[2649] + v[2673] * v[319] + v[2625] * v[336])
		- v[2685] * v[358] - v[1298] * v[544] - v[1374] * v[675] - v[1177] * v[822];
	v[2693] = gti[0] * (v[1149] * v[2648] + v[2672] * v[316] + v[2624] * v[331]) + gti[1] * (v[1150] * v[2648]
		+ v[2672] * v[318] + v[2624] * v[334]) + gti[2] * (v[1151] * v[2648] + v[2672] * v[319] + v[2624] * v[336])
		- v[2684] * v[358] - v[1297] * v[544] - v[1374] * v[673] - v[1177] * v[816];
	v[2692] = gti[0] * (v[1149] * v[2647] + v[2671] * v[316] + v[2623] * v[331]) + gti[1] * (v[1150] * v[2647]
		+ v[2671] * v[318] + v[2623] * v[334]) + gti[2] * (v[1151] * v[2647] + v[2671] * v[319] + v[2623] * v[336])
		- v[2683] * v[358] - v[1296] * v[544] - v[1374] * v[671] - v[1177] * v[810];
	v[1871] = -(v[1177] * v[1573]) + gti[0] * (v[1149] * v[1839] + v[1851] * v[316] + v[1827] * v[331]) + gti[1] *
		(v[1150] * v[1839] + v[1851] * v[318] + v[1827] * v[334]) + gti[2] * (v[1151] * v[1839] + v[1851] * v[319]
			+ v[1827] * v[336]) - v[1861] * v[363] - v[1375] * v[4764];
	v[1870] = -(v[1177] * v[1578]) + gti[0] * (v[1149] * v[1841] + v[1853] * v[316] + v[1829] * v[331]) + gti[1] *
		(v[1150] * v[1841] + v[1853] * v[318] + v[1829] * v[334]) + gti[2] * (v[1151] * v[1841] + v[1853] * v[319]
			+ v[1829] * v[336]) - v[1861] * v[361] - v[1375] * v[4618];
	v[1869] = -(v[1177] * v[1583]) + gti[0] * (v[1149] * v[1843] + v[1855] * v[316] + v[1831] * v[331]) + gti[1] *
		(v[1150] * v[1843] + v[1855] * v[318] + v[1831] * v[334]) + gti[2] * (v[1151] * v[1843] + v[1855] * v[319]
			+ v[1831] * v[336]) - v[1861] * v[358] + v[1375] * v[4765];
	v[1868] = -(v[1177] * v[1571]) + gti[0] * (v[1149] * v[1838] + v[1850] * v[316] + v[1826] * v[331]) + gti[1] *
		(v[1150] * v[1838] + v[1850] * v[318] + v[1826] * v[334]) + gti[2] * (v[1151] * v[1838] + v[1850] * v[319]
			+ v[1826] * v[336]) - v[1860] * v[363] - v[1375] * v[548] - v[1374] * v[549];
	v[1867] = -(v[1177] * v[1575]) + gti[0] * (v[1149] * v[1845] + v[1857] * v[316] + v[1833] * v[331]) + gti[1] *
		(v[1150] * v[1845] + v[1857] * v[318] + v[1833] * v[334]) + gti[2] * (v[1151] * v[1845] + v[1857] * v[319]
			+ v[1833] * v[336]) - v[1862] * v[363] - v[1374] * v[4766];
	v[1866] = -(v[1177] * v[1576]) + gti[0] * (v[1149] * v[1840] + v[1852] * v[316] + v[1828] * v[331]) + gti[1] *
		(v[1150] * v[1840] + v[1852] * v[318] + v[1828] * v[334]) + gti[2] * (v[1151] * v[1840] + v[1852] * v[319]
			+ v[1828] * v[336]) - v[1860] * v[361] - v[1375] * v[546] - v[1374] * v[547];
	v[1865] = -(v[1177] * v[1580]) + gti[0] * (v[1149] * v[1846] + v[1858] * v[316] + v[1834] * v[331]) + gti[1] *
		(v[1150] * v[1846] + v[1858] * v[318] + v[1834] * v[334]) + gti[2] * (v[1151] * v[1846] + v[1858] * v[319]
			+ v[1834] * v[336]) - v[1862] * v[361] - v[1374] * v[4616];
	v[1864] = -(v[1177] * v[1581]) + gti[0] * (v[1149] * v[1842] + v[1854] * v[316] + v[1830] * v[331]) + gti[1] *
		(v[1150] * v[1842] + v[1854] * v[318] + v[1830] * v[334]) + gti[2] * (v[1151] * v[1842] + v[1854] * v[319]
			+ v[1830] * v[336]) - v[1860] * v[358] - v[1375] * v[544] - v[1374] * v[545];
	v[1863] = -(v[1177] * v[1589]) + gti[0] * (v[1149] * v[1847] + v[1859] * v[316] + v[1835] * v[331]) + gti[1] *
		(v[1150] * v[1847] + v[1859] * v[318] + v[1835] * v[334]) + gti[2] * (v[1151] * v[1847] + v[1859] * v[319]
			+ v[1835] * v[336]) - v[1862] * v[358] + v[1374] * v[4767];
	v[1170] = v[1149] * v[1158] + v[1155] * v[316] + v[1161] * v[331];
	v[1171] = v[1149] * v[1159] + v[1156] * v[316] + v[1162] * v[331];
	v[1172] = v[1149] * v[1160] + v[1157] * v[316] + v[1163] * v[331];
	v[1174] = v[1150] * v[1158] + v[1155] * v[318] + v[1161] * v[334];
	v[1175] = v[1150] * v[1159] + v[1156] * v[318] + v[1162] * v[334];
	v[1176] = v[1150] * v[1160] + v[1157] * v[318] + v[1163] * v[334];
	v[1179] = v[1151] * v[1158] + v[1155] * v[319] + v[1161] * v[336];
	v[1180] = v[1151] * v[1159] + v[1156] * v[319] + v[1162] * v[336];
	v[1181] = v[1151] * v[1160] + v[1157] * v[319] + v[1163] * v[336];
	v[1183] = (*epst)*(v[1165] + gti[0] * (v[1149] * v[1171] + v[1170] * v[316] + v[1172] * v[331]) + gti[1] *
		(v[1150] * v[1171] + v[1170] * v[318] + v[1172] * v[334]) + gti[2] * (v[1151] * v[1171] + v[1170] * v[319]
			+ v[1172] * v[336]) - v[1177] * v[358]);
	v[1184] = (*epst)*(v[1167] + gti[0] * (v[1149] * v[1175] + v[1174] * v[316] + v[1176] * v[331]) + gti[1] *
		(v[1150] * v[1175] + v[1174] * v[318] + v[1176] * v[334]) + gti[2] * (v[1151] * v[1175] + v[1174] * v[319]
			+ v[1176] * v[336]) - v[1177] * v[361]);
	v[1185] = (*epst)*(v[1169] + gti[0] * (v[1149] * v[1180] + v[1179] * v[316] + v[1181] * v[331]) + gti[1] *
		(v[1150] * v[1180] + v[1179] * v[318] + v[1181] * v[334]) + gti[2] * (v[1151] * v[1180] + v[1179] * v[319]
			+ v[1181] * v[336]) - v[1177] * v[363]);
	v[3102] = v[1183] * v[2693] + v[1184] * v[2701] + v[1185] * v[2708];
	v[3101] = v[1183] * v[2717] + v[1184] * v[2725] + v[1185] * v[2732];
	v[3095] = v[1183] * v[2692] + v[1184] * v[2700] + v[1185] * v[2707];
	v[3094] = v[1183] * v[2716] + v[1184] * v[2724] + v[1185] * v[2731];
	v[3088] = v[1183] * v[2691] + v[1184] * v[2699] + v[1185] * v[2706];
	v[3087] = v[1183] * v[2715] + v[1184] * v[2723] + v[1185] * v[2730];
	v[3082] = v[1183] * v[2690] + v[1184] * v[2698] + v[1185] * v[2705];
	v[3081] = v[1183] * v[2714] + v[1184] * v[2722] + v[1185] * v[2729];
	v[3075] = v[1183] * v[2689] + v[1184] * v[2697] + v[1185] * v[2704];
	v[3074] = v[1183] * v[2713] + v[1184] * v[2721] + v[1185] * v[2728];
	v[3066] = v[1183] * v[2712] + v[1184] * v[2720] + v[1185] * v[2727];
	v[3065] = v[1183] * v[2688] + v[1184] * v[2696] + v[1185] * v[2703];
	v[3051] = v[1183] * v[2711] + v[1184] * v[2719] + v[1185] * v[2720];
	v[3048] = v[1183] * v[2687] + v[1184] * v[2695] + v[1185] * v[2696];
	v[3032] = v[1183] * v[2710] + v[1184] * v[2711] + v[1185] * v[2712];
	v[3029] = v[1183] * v[2686] + v[1184] * v[2687] + v[1185] * v[2688];
	v[1327] = -(v[358] * v[361]);
	v[1328] = -(v[358] * v[363]);
	v[1330] = -(v[361] * v[363]);
	v[1332] = v[1241] - v[1293] * v[358];
	v[1333] = v[1260] - v[1293] * v[361];
	v[1334] = v[1280] - v[1293] * v[363];
	v[1335] = v[1246] - v[1294] * v[358];
	v[1336] = v[1264] - v[1294] * v[361];
	v[1337] = v[1285] - v[1294] * v[363];
	v[1338] = v[1250] - v[1295] * v[358];
	v[1339] = v[1269] - v[1295] * v[361];
	v[1340] = v[1289] - v[1295] * v[363];
	v[1394] = gti[0] * (v[1149] * v[1378] + v[1376] * v[316] + v[1380] * v[331]) + gti[1] * (v[1150] * v[1378]
		+ v[1376] * v[318] + v[1380] * v[334]) + gti[2] * (v[1151] * v[1378] + v[1376] * v[319] + v[1380] * v[336])
		- v[1374] * v[358] - v[1177] * v[544];
	v[1395] = gti[0] * (v[1149] * v[1384] + v[1382] * v[316] + v[1386] * v[331]) + gti[1] * (v[1150] * v[1384]
		+ v[1382] * v[318] + v[1386] * v[334]) + gti[2] * (v[1151] * v[1384] + v[1382] * v[319] + v[1386] * v[336])
		- v[1374] * v[361] - v[1177] * v[546];
	v[1396] = gti[0] * (v[1149] * v[1390] + v[1388] * v[316] + v[1392] * v[331]) + gti[1] * (v[1150] * v[1390]
		+ v[1388] * v[318] + v[1392] * v[334]) + gti[2] * (v[1151] * v[1390] + v[1388] * v[319] + v[1392] * v[336])
		- v[1374] * v[363] - v[1177] * v[548];
	v[1397] = gti[0] * (v[1149] * v[1379] + v[1377] * v[316] + v[1381] * v[331]) + gti[1] * (v[1150] * v[1379]
		+ v[1377] * v[318] + v[1381] * v[334]) + gti[2] * (v[1151] * v[1379] + v[1377] * v[319] + v[1381] * v[336])
		- v[1375] * v[358] - v[1177] * v[545];
	v[4247] = v[1394] * v[1789] + v[1397] * v[1803] + v[1864] * v[699] + v[1869] * v[705];
	v[4248] = v[2710] + v[4247];
	v[4245] = v[1394] * v[1788] + v[1397] * v[1800] + v[1863] * v[699] + v[1864] * v[705];
	v[4246] = v[2686] + v[4245];
	v[4243] = v[1394] * v[1791] + v[1397] * v[1807] + v[1864] * v[700] + v[1869] * v[706];
	v[4244] = v[2711] + v[4243];
	v[4241] = v[1394] * v[1790] + v[1397] * v[1804] + v[1863] * v[700] + v[1864] * v[706];
	v[4242] = v[2687] + v[4241];
	v[4239] = v[1394] * v[1793] + v[1397] * v[1811] + v[1864] * v[701] + v[1869] * v[707];
	v[4240] = v[2712] + v[4239];
	v[4237] = v[1394] * v[1792] + v[1397] * v[1808] + v[1863] * v[701] + v[1864] * v[707];
	v[4238] = v[2688] + v[4237];
	v[4235] = v[1394] * v[1795] + v[1397] * v[1815] + v[1864] * v[702] + v[1869] * v[708];
	v[4236] = v[2716] + v[4235];
	v[4233] = v[1394] * v[1794] + v[1397] * v[1812] + v[1863] * v[702] + v[1864] * v[708];
	v[4234] = v[2692] + v[4233];
	v[4231] = v[1394] * v[1797] + v[1397] * v[1819] + v[1864] * v[703] + v[1869] * v[709];
	v[4232] = v[2717] + v[4231];
	v[4229] = v[1394] * v[1796] + v[1397] * v[1816] + v[1863] * v[703] + v[1864] * v[709];
	v[4230] = v[2693] + v[4229];
	v[4227] = v[1394] * v[1799] + v[1397] * v[1823] + v[1864] * v[704] + v[1869] * v[710];
	v[4228] = v[2718] + v[4227];
	v[4225] = v[1394] * v[1798] + v[1397] * v[1820] + v[1863] * v[704] + v[1864] * v[710];
	v[4226] = v[2694] + v[4225];
	v[3833] = v[1394] * v[2487] + v[1397] * v[2547] + v[2694] * v[699] + v[2718] * v[705];
	v[3834] = v[3010] + v[3833];
	v[3831] = v[1394] * v[2486] + v[1397] * v[2546] + v[2693] * v[699] + v[2717] * v[705];
	v[3832] = v[3009] + v[3831];
	v[3829] = v[1394] * v[2485] + v[1397] * v[2545] + v[2692] * v[699] + v[2716] * v[705];
	v[3830] = v[3008] + v[3829];
	v[3825] = v[1394] * v[2484] + v[1397] * v[2544] + v[4772];
	v[3824] = v[1394] * v[2483] + v[1397] * v[2543] + v[4768];
	v[3823] = v[1394] * v[2482] + v[1397] * v[2542] + v[2686] * v[699] + v[2710] * v[705];
	v[3821] = v[1394] * v[2493] + v[1397] * v[2553] + v[2694] * v[700] + v[2718] * v[706];
	v[3822] = v[3013] + v[3821];
	v[3819] = v[1394] * v[2492] + v[1397] * v[2552] + v[2693] * v[700] + v[2717] * v[706];
	v[3820] = v[3012] + v[3819];
	v[3817] = v[1394] * v[2491] + v[1397] * v[2551] + v[2692] * v[700] + v[2716] * v[706];
	v[3818] = v[3011] + v[3817];
	v[3813] = v[1394] * v[2490] + v[1397] * v[2550] + v[4774];
	v[3812] = v[1394] * v[2489] + v[1397] * v[2549] + v[4769];
	v[3811] = v[1394] * v[2488] + v[1397] * v[2548] + v[2686] * v[700] + v[2710] * v[706];
	v[3809] = v[1394] * v[2499] + v[1397] * v[2559] + v[2694] * v[701] + v[2718] * v[707];
	v[3810] = v[3016] + v[3809];
	v[3807] = v[1394] * v[2498] + v[1397] * v[2558] + v[2693] * v[701] + v[2717] * v[707];
	v[3808] = v[3015] + v[3807];
	v[3805] = v[1394] * v[2497] + v[1397] * v[2557] + v[2692] * v[701] + v[2716] * v[707];
	v[3806] = v[3014] + v[3805];
	v[3801] = v[1394] * v[2496] + v[1397] * v[2556] + v[4776];
	v[3800] = v[1394] * v[2495] + v[1397] * v[2555] + v[4770];
	v[3799] = v[1394] * v[2494] + v[1397] * v[2554] + v[2686] * v[701] + v[2710] * v[707];
	v[3797] = v[1394] * v[2511] + v[1397] * v[2571] + v[2694] * v[702] + v[2718] * v[708];
	v[3798] = v[3766] + v[3797];
	v[3795] = v[1394] * v[2510] + v[1397] * v[2570] + v[2693] * v[702] + v[2717] * v[708];
	v[3796] = v[3779] + v[3795];
	v[3793] = v[1394] * v[2509] + v[1397] * v[2569] + v[2692] * v[702] + v[2716] * v[708];
	v[3794] = v[3792] + v[3793];
	v[3788] = v[1394] * v[2508] + v[1397] * v[2568] + v[3947];
	v[3787] = v[1394] * v[2505] + v[1397] * v[2565] + v[3863];
	v[3786] = v[1394] * v[2502] + v[1397] * v[2562] + v[3008] + v[2686] * v[702] + v[2710] * v[708];
	v[3784] = v[1394] * v[2525] + v[1397] * v[2585] + v[2694] * v[703] + v[2718] * v[709];
	v[3785] = v[3768] + v[3784];
	v[3782] = v[1394] * v[2524] + v[1397] * v[2584] + v[2693] * v[703] + v[2717] * v[709];
	v[3783] = v[3781] + v[3782];
	v[3780] = v[1394] * v[2523] + v[1397] * v[2583] + v[3779] + v[2692] * v[703] + v[2716] * v[709];
	v[3775] = v[1394] * v[2520] + v[1397] * v[2580] + v[3932];
	v[3774] = v[1394] * v[2517] + v[1397] * v[2577] + v[3849];
	v[3773] = v[1394] * v[2514] + v[1397] * v[2574] + v[3009] + v[2686] * v[703] + v[2710] * v[709];
	v[3771] = v[1394] * v[2541] + v[1397] * v[2601] + v[2694] * v[704] + v[2718] * v[710];
	v[3772] = v[3770] + v[3771];
	v[3769] = v[1394] * v[2540] + v[1397] * v[2600] + v[3768] + v[2693] * v[704] + v[2717] * v[710];
	v[3767] = v[1394] * v[2537] + v[1397] * v[2597] + v[3766] + v[2692] * v[704] + v[2716] * v[710];
	v[3762] = v[1394] * v[2534] + v[1397] * v[2594] + v[3917];
	v[3761] = v[1394] * v[2531] + v[1397] * v[2591] + v[3835];
	v[3760] = v[1394] * v[2528] + v[1397] * v[2588] + v[3010] + v[2686] * v[704] + v[2710] * v[710];
	v[1432] = -v[1253] + gti[0] * (v[1149] * v[1304] + v[1301] * v[316] + v[1307] * v[331]) + gti[1] * (v[1150] * v[1304]
		+ v[1301] * v[318] + v[1307] * v[334]) + gti[2] * (v[1151] * v[1304] + v[1301] * v[319] + v[1307] * v[336])
		- v[1298] * v[358] - v[1177] * v[675] + v[1394] * v[704] + v[1397] * v[710];
	v[1428] = -v[1252] + gti[0] * (v[1149] * v[1303] + v[1300] * v[316] + v[1306] * v[331]) + gti[1] * (v[1150] * v[1303]
		+ v[1300] * v[318] + v[1306] * v[334]) + gti[2] * (v[1151] * v[1303] + v[1300] * v[319] + v[1306] * v[336])
		- v[1297] * v[358] - v[1177] * v[673] + v[1394] * v[703] + v[1397] * v[709];
	v[1424] = -v[1251] + gti[0] * (v[1149] * v[1302] + v[1299] * v[316] + v[1305] * v[331]) + gti[1] * (v[1150] * v[1302]
		+ v[1299] * v[318] + v[1305] * v[334]) + gti[2] * (v[1151] * v[1302] + v[1299] * v[319] + v[1305] * v[336])
		- v[1296] * v[358] - v[1177] * v[671] + v[1394] * v[702] + v[1397] * v[708];
	v[1417] = v[1328] + v[1394] * v[701] + v[1397] * v[707];
	v[1413] = v[1327] + v[1394] * v[700] + v[1397] * v[706];
	v[1409] = 1e0 - (v[358] * v[358]) + v[1394] * v[699] + v[1397] * v[705];
	v[1398] = gti[0] * (v[1149] * v[1385] + v[1383] * v[316] + v[1387] * v[331]) + gti[1] * (v[1150] * v[1385]
		+ v[1383] * v[318] + v[1387] * v[334]) + gti[2] * (v[1151] * v[1385] + v[1383] * v[319] + v[1387] * v[336])
		- v[1375] * v[361] - v[1177] * v[547];
	v[4271] = v[1395] * v[1789] + v[1398] * v[1803] + v[1866] * v[699] + v[1870] * v[705];
	v[4272] = v[2711] + v[4271];
	v[4269] = v[1395] * v[1788] + v[1398] * v[1800] + v[1865] * v[699] + v[1866] * v[705];
	v[4270] = v[2687] + v[4269];
	v[4267] = v[1395] * v[1791] + v[1398] * v[1807] + v[1866] * v[700] + v[1870] * v[706];
	v[4268] = v[2719] + v[4267];
	v[4265] = v[1395] * v[1790] + v[1398] * v[1804] + v[1865] * v[700] + v[1866] * v[706];
	v[4266] = v[2695] + v[4265];
	v[4263] = v[1395] * v[1793] + v[1398] * v[1811] + v[1866] * v[701] + v[1870] * v[707];
	v[4264] = v[2720] + v[4263];
	v[4261] = v[1395] * v[1792] + v[1398] * v[1808] + v[1865] * v[701] + v[1866] * v[707];
	v[4262] = v[2696] + v[4261];
	v[4259] = v[1395] * v[1795] + v[1398] * v[1815] + v[1866] * v[702] + v[1870] * v[708];
	v[4260] = v[2724] + v[4259];
	v[4257] = v[1395] * v[1794] + v[1398] * v[1812] + v[1865] * v[702] + v[1866] * v[708];
	v[4258] = v[2700] + v[4257];
	v[4255] = v[1395] * v[1797] + v[1398] * v[1819] + v[1866] * v[703] + v[1870] * v[709];
	v[4256] = v[2725] + v[4255];
	v[4253] = v[1395] * v[1796] + v[1398] * v[1816] + v[1865] * v[703] + v[1866] * v[709];
	v[4254] = v[2701] + v[4253];
	v[4251] = v[1395] * v[1799] + v[1398] * v[1823] + v[1866] * v[704] + v[1870] * v[710];
	v[4252] = v[2726] + v[4251];
	v[4249] = v[1395] * v[1798] + v[1398] * v[1820] + v[1865] * v[704] + v[1866] * v[710];
	v[4250] = v[2702] + v[4249];
	v[3915] = v[1395] * v[2487] + v[1398] * v[2547] + v[2702] * v[699] + v[2726] * v[705];
	v[3916] = v[3013] + v[3915];
	v[3913] = v[1395] * v[2486] + v[1398] * v[2546] + v[2701] * v[699] + v[2725] * v[705];
	v[3914] = v[3012] + v[3913];
	v[3911] = v[1395] * v[2485] + v[1398] * v[2545] + v[2700] * v[699] + v[2724] * v[705];
	v[3912] = v[3011] + v[3911];
	v[3907] = v[1395] * v[2484] + v[1398] * v[2544] + v[4771];
	v[3906] = v[1395] * v[2483] + v[1398] * v[2543] + v[2695] * v[699] + v[2719] * v[705];
	v[3905] = v[1395] * v[2482] + v[1398] * v[2542] + v[4768];
	v[3903] = v[1395] * v[2493] + v[1398] * v[2553] + v[2702] * v[700] + v[2726] * v[706];
	v[3904] = v[3019] + v[3903];
	v[3901] = v[1395] * v[2492] + v[1398] * v[2552] + v[2701] * v[700] + v[2725] * v[706];
	v[3902] = v[3018] + v[3901];
	v[3899] = v[1395] * v[2491] + v[1398] * v[2551] + v[2700] * v[700] + v[2724] * v[706];
	v[3900] = v[3017] + v[3899];
	v[3895] = v[1395] * v[2490] + v[1398] * v[2550] + v[4773];
	v[3894] = v[1395] * v[2489] + v[1398] * v[2549] + v[2695] * v[700] + v[2719] * v[706];
	v[3893] = v[1395] * v[2488] + v[1398] * v[2548] + v[4769];
	v[3889] = v[1395] * v[2499] + v[1398] * v[2559] + v[2702] * v[701] + v[2726] * v[707];
	v[3890] = v[3022] + v[3889];
	v[3887] = v[1395] * v[2498] + v[1398] * v[2558] + v[2701] * v[701] + v[2725] * v[707];
	v[3888] = v[3021] + v[3887];
	v[3885] = v[1395] * v[2497] + v[1398] * v[2557] + v[2700] * v[701] + v[2724] * v[707];
	v[3886] = v[3020] + v[3885];
	v[3881] = v[1395] * v[2496] + v[1398] * v[2556] + v[4775];
	v[3880] = v[1395] * v[2495] + v[1398] * v[2555] + v[2695] * v[701] + v[2719] * v[707];
	v[3879] = v[1395] * v[2494] + v[1398] * v[2554] + v[4770];
	v[3875] = v[1395] * v[2511] + v[1398] * v[2571] + v[2702] * v[702] + v[2726] * v[708];
	v[3876] = v[3842] + v[3875];
	v[3873] = v[1395] * v[2510] + v[1398] * v[2570] + v[2701] * v[702] + v[2725] * v[708];
	v[3874] = v[3856] + v[3873];
	v[3871] = v[1395] * v[2509] + v[1398] * v[2569] + v[2700] * v[702] + v[2724] * v[708];
	v[3872] = v[3870] + v[3871];
	v[3866] = v[1395] * v[2508] + v[1398] * v[2568] + v[3949];
	v[3865] = v[1395] * v[2505] + v[1398] * v[2565] + v[3017] + v[2695] * v[702] + v[2719] * v[708];
	v[3864] = v[1395] * v[2502] + v[1398] * v[2562] + v[3863];
	v[3861] = v[1395] * v[2525] + v[1398] * v[2585] + v[2702] * v[703] + v[2726] * v[709];
	v[3862] = v[3844] + v[3861];
	v[3859] = v[1395] * v[2524] + v[1398] * v[2584] + v[2701] * v[703] + v[2725] * v[709];
	v[3860] = v[3858] + v[3859];
	v[3857] = v[1395] * v[2523] + v[1398] * v[2583] + v[3856] + v[2700] * v[703] + v[2724] * v[709];
	v[3852] = v[1395] * v[2520] + v[1398] * v[2580] + v[3934];
	v[3851] = v[1395] * v[2517] + v[1398] * v[2577] + v[3018] + v[2695] * v[703] + v[2719] * v[709];
	v[3850] = v[1395] * v[2514] + v[1398] * v[2574] + v[3849];
	v[3847] = v[1395] * v[2541] + v[1398] * v[2601] + v[2702] * v[704] + v[2726] * v[710];
	v[3848] = v[3846] + v[3847];
	v[3845] = v[1395] * v[2540] + v[1398] * v[2600] + v[3844] + v[2701] * v[704] + v[2725] * v[710];
	v[3843] = v[1395] * v[2537] + v[1398] * v[2597] + v[3842] + v[2700] * v[704] + v[2724] * v[710];
	v[3838] = v[1395] * v[2534] + v[1398] * v[2594] + v[3919];
	v[3837] = v[1395] * v[2531] + v[1398] * v[2591] + v[3019] + v[2695] * v[704] + v[2719] * v[710];
	v[3836] = v[1395] * v[2528] + v[1398] * v[2588] + v[3835];
	v[1433] = -v[1272] + gti[0] * (v[1149] * v[1313] + v[1310] * v[316] + v[1316] * v[331]) + gti[1] * (v[1150] * v[1313]
		+ v[1310] * v[318] + v[1316] * v[334]) + gti[2] * (v[1151] * v[1313] + v[1310] * v[319] + v[1316] * v[336])
		- v[1298] * v[361] - v[1177] * v[681] + v[1395] * v[704] + v[1398] * v[710];
	v[1429] = -v[1271] + gti[0] * (v[1149] * v[1312] + v[1309] * v[316] + v[1315] * v[331]) + gti[1] * (v[1150] * v[1312]
		+ v[1309] * v[318] + v[1315] * v[334]) + gti[2] * (v[1151] * v[1312] + v[1309] * v[319] + v[1315] * v[336])
		- v[1297] * v[361] - v[1177] * v[679] + v[1395] * v[703] + v[1398] * v[709];
	v[1425] = -v[1270] + gti[0] * (v[1149] * v[1311] + v[1308] * v[316] + v[1314] * v[331]) + gti[1] * (v[1150] * v[1311]
		+ v[1308] * v[318] + v[1314] * v[334]) + gti[2] * (v[1151] * v[1311] + v[1308] * v[319] + v[1314] * v[336])
		- v[1296] * v[361] - v[1177] * v[677] + v[1395] * v[702] + v[1398] * v[708];
	v[1418] = v[1330] + v[1395] * v[701] + v[1398] * v[707];
	v[1414] = 1e0 - (v[361] * v[361]) + v[1395] * v[700] + v[1398] * v[706];
	v[1410] = v[1327] + v[1395] * v[699] + v[1398] * v[705];
	v[1399] = gti[0] * (v[1149] * v[1391] + v[1389] * v[316] + v[1393] * v[331]) + gti[1] * (v[1150] * v[1391]
		+ v[1389] * v[318] + v[1393] * v[334]) + gti[2] * (v[1151] * v[1391] + v[1389] * v[319] + v[1393] * v[336])
		- v[1375] * v[363] - v[1177] * v[549];
	v[4295] = v[1396] * v[1789] + v[1399] * v[1803] + v[1868] * v[699] + v[1871] * v[705];
	v[4296] = v[2712] + v[4295];
	v[4293] = v[1396] * v[1788] + v[1399] * v[1800] + v[1867] * v[699] + v[1868] * v[705];
	v[4294] = v[2688] + v[4293];
	v[4291] = v[1396] * v[1791] + v[1399] * v[1807] + v[1868] * v[700] + v[1871] * v[706];
	v[4292] = v[2720] + v[4291];
	v[4289] = v[1396] * v[1790] + v[1399] * v[1804] + v[1867] * v[700] + v[1868] * v[706];
	v[4290] = v[2696] + v[4289];
	v[4287] = v[1396] * v[1793] + v[1399] * v[1811] + v[1868] * v[701] + v[1871] * v[707];
	v[4288] = v[2727] + v[4287];
	v[4285] = v[1396] * v[1792] + v[1399] * v[1808] + v[1867] * v[701] + v[1868] * v[707];
	v[4286] = v[2703] + v[4285];
	v[4283] = v[1396] * v[1795] + v[1399] * v[1815] + v[1868] * v[702] + v[1871] * v[708];
	v[4284] = v[2731] + v[4283];
	v[4281] = v[1396] * v[1794] + v[1399] * v[1812] + v[1867] * v[702] + v[1868] * v[708];
	v[4282] = v[2707] + v[4281];
	v[4279] = v[1396] * v[1797] + v[1399] * v[1819] + v[1868] * v[703] + v[1871] * v[709];
	v[4280] = v[2732] + v[4279];
	v[4277] = v[1396] * v[1796] + v[1399] * v[1816] + v[1867] * v[703] + v[1868] * v[709];
	v[4278] = v[2708] + v[4277];
	v[4275] = v[1396] * v[1799] + v[1399] * v[1823] + v[1868] * v[704] + v[1871] * v[710];
	v[4276] = v[2733] + v[4275];
	v[4273] = v[1396] * v[1798] + v[1399] * v[1820] + v[1867] * v[704] + v[1868] * v[710];
	v[4274] = v[2709] + v[4273];
	v[4002] = v[1396] * v[2487] + v[1399] * v[2547] + v[2709] * v[699] + v[2733] * v[705];
	v[4003] = v[3016] + v[4002];
	v[4000] = v[1396] * v[2486] + v[1399] * v[2546] + v[2708] * v[699] + v[2732] * v[705];
	v[4001] = v[3015] + v[4000];
	v[3998] = v[1396] * v[2485] + v[1399] * v[2545] + v[2707] * v[699] + v[2731] * v[705];
	v[3999] = v[3014] + v[3998];
	v[3994] = v[1396] * v[2484] + v[1399] * v[2544] + v[2703] * v[699] + v[2727] * v[705];
	v[3993] = v[1396] * v[2483] + v[1399] * v[2543] + v[4771];
	v[3992] = v[1396] * v[2482] + v[1399] * v[2542] + v[4772];
	v[3990] = v[1396] * v[2493] + v[1399] * v[2553] + v[2709] * v[700] + v[2733] * v[706];
	v[3991] = v[3022] + v[3990];
	v[3988] = v[1396] * v[2492] + v[1399] * v[2552] + v[2708] * v[700] + v[2732] * v[706];
	v[3989] = v[3021] + v[3988];
	v[3986] = v[1396] * v[2491] + v[1399] * v[2551] + v[2707] * v[700] + v[2731] * v[706];
	v[3987] = v[3020] + v[3986];
	v[3982] = v[1396] * v[2490] + v[1399] * v[2550] + v[2703] * v[700] + v[2727] * v[706];
	v[3981] = v[1396] * v[2489] + v[1399] * v[2549] + v[4773];
	v[3980] = v[1396] * v[2488] + v[1399] * v[2548] + v[4774];
	v[3976] = v[1396] * v[2499] + v[1399] * v[2559] + v[2709] * v[701] + v[2733] * v[707];
	v[3977] = v[3025] + v[3976];
	v[3974] = v[1396] * v[2498] + v[1399] * v[2558] + v[2708] * v[701] + v[2732] * v[707];
	v[3975] = v[3024] + v[3974];
	v[3972] = v[1396] * v[2497] + v[1399] * v[2557] + v[2707] * v[701] + v[2731] * v[707];
	v[3973] = v[3023] + v[3972];
	v[3968] = v[1396] * v[2496] + v[1399] * v[2556] + v[2703] * v[701] + v[2727] * v[707];
	v[3967] = v[1396] * v[2495] + v[1399] * v[2555] + v[4775];
	v[3964] = v[1396] * v[2494] + v[1399] * v[2554] + v[4776];
	v[3960] = v[1396] * v[2511] + v[1399] * v[2571] + v[2709] * v[702] + v[2733] * v[708];
	v[3961] = v[3925] + v[3960];
	v[3958] = v[1396] * v[2510] + v[1399] * v[2570] + v[2708] * v[702] + v[2732] * v[708];
	v[3959] = v[3940] + v[3958];
	v[3956] = v[1396] * v[2509] + v[1399] * v[2569] + v[2707] * v[702] + v[2731] * v[708];
	v[3957] = v[3955] + v[3956];
	v[3951] = v[1396] * v[2508] + v[1399] * v[2568] + v[3023] + v[2703] * v[702] + v[2727] * v[708];
	v[3950] = v[1396] * v[2505] + v[1399] * v[2565] + v[3949];
	v[3948] = v[1396] * v[2502] + v[1399] * v[2562] + v[3947];
	v[3945] = v[1396] * v[2525] + v[1399] * v[2585] + v[2709] * v[703] + v[2733] * v[709];
	v[3946] = v[3927] + v[3945];
	v[3943] = v[1396] * v[2524] + v[1399] * v[2584] + v[2708] * v[703] + v[2732] * v[709];
	v[3944] = v[3942] + v[3943];
	v[3941] = v[1396] * v[2523] + v[1399] * v[2583] + v[3940] + v[2707] * v[703] + v[2731] * v[709];
	v[3936] = v[1396] * v[2520] + v[1399] * v[2580] + v[3024] + v[2703] * v[703] + v[2727] * v[709];
	v[3935] = v[1396] * v[2517] + v[1399] * v[2577] + v[3934];
	v[3933] = v[1396] * v[2514] + v[1399] * v[2574] + v[3932];
	v[3930] = v[1396] * v[2541] + v[1399] * v[2601] + v[2709] * v[704] + v[2733] * v[710];
	v[3931] = v[3929] + v[3930];
	v[3928] = v[1396] * v[2540] + v[1399] * v[2600] + v[3927] + v[2708] * v[704] + v[2732] * v[710];
	v[3926] = v[1396] * v[2537] + v[1399] * v[2597] + v[3925] + v[2707] * v[704] + v[2731] * v[710];
	v[3921] = v[1396] * v[2534] + v[1399] * v[2594] + v[3025] + v[2703] * v[704] + v[2727] * v[710];
	v[3920] = v[1396] * v[2531] + v[1399] * v[2591] + v[3919];
	v[3918] = v[1396] * v[2528] + v[1399] * v[2588] + v[3917];
	v[3100] = v[1183] * v[4231] + v[1184] * v[4255] + v[1185] * v[4279];
	v[4850] = v[3100] + v[3101];
	v[3099] = v[1183] * v[4229] + v[1184] * v[4253] + v[1185] * v[4277];
	v[4849] = v[3099] + v[3102];
	v[3093] = v[1183] * v[4235] + v[1184] * v[4259] + v[1185] * v[4283];
	v[4783] = v[3093] + v[3094];
	v[3092] = v[1183] * v[4233] + v[1184] * v[4257] + v[1185] * v[4281];
	v[4784] = v[3092] + v[3095];
	v[3064] = v[1183] * v[4239] + v[1184] * v[4263] + v[1185] * v[4287];
	v[4782] = v[3064] + v[3066];
	v[3063] = v[1183] * v[4237] + v[1184] * v[4261] + v[1185] * v[4285];
	v[4781] = v[3063] + v[3065];
	v[3047] = v[1183] * v[4243] + v[1184] * v[4267] + v[1185] * v[4291];
	v[4780] = v[3047] + v[3051];
	v[3046] = v[1183] * v[4241] + v[1184] * v[4265] + v[1185] * v[4289];
	v[4779] = v[3046] + v[3048];
	v[3028] = v[1183] * v[4247] + v[1184] * v[4271] + v[1185] * v[4295];
	v[4778] = v[3028] + v[3032];
	v[3027] = v[1183] * v[4245] + v[1184] * v[4269] + v[1185] * v[4293];
	v[4777] = v[3027] + v[3029];
	v[1434] = -v[1292] + gti[0] * (v[1149] * v[1322] + v[1319] * v[316] + v[1325] * v[331]) + gti[1] * (v[1150] * v[1322]
		+ v[1319] * v[318] + v[1325] * v[334]) + gti[2] * (v[1151] * v[1322] + v[1319] * v[319] + v[1325] * v[336])
		- v[1298] * v[363] - v[1177] * v[687] + v[1396] * v[704] + v[1399] * v[710];
	v[1430] = -v[1291] + gti[0] * (v[1149] * v[1321] + v[1318] * v[316] + v[1324] * v[331]) + gti[1] * (v[1150] * v[1321]
		+ v[1318] * v[318] + v[1324] * v[334]) + gti[2] * (v[1151] * v[1321] + v[1318] * v[319] + v[1324] * v[336])
		- v[1297] * v[363] - v[1177] * v[685] + v[1396] * v[703] + v[1399] * v[709];
	v[1426] = -v[1290] + gti[0] * (v[1149] * v[1320] + v[1317] * v[316] + v[1323] * v[331]) + gti[1] * (v[1150] * v[1320]
		+ v[1317] * v[318] + v[1323] * v[334]) + gti[2] * (v[1151] * v[1320] + v[1317] * v[319] + v[1323] * v[336])
		- v[1296] * v[363] - v[1177] * v[683] + v[1396] * v[702] + v[1399] * v[708];
	v[1419] = 1e0 - (v[363] * v[363]) + v[1396] * v[701] + v[1399] * v[707];
	v[1415] = v[1330] + v[1396] * v[700] + v[1399] * v[706];
	v[1411] = v[1328] + v[1396] * v[699] + v[1399] * v[705];
	v[1416] = (*epst)*(v[1409] * v[1413] + v[1410] * v[1414] + v[1411] * v[1415]);
	v[1420] = (*epst)*(v[1409] * v[1417] + v[1410] * v[1418] + v[1411] * v[1419]);
	v[1421] = (*epst)*(v[1332] * v[1409] + v[1333] * v[1410] + v[1334] * v[1411]);
	v[1422] = (*epst)*(v[1335] * v[1409] + v[1336] * v[1410] + v[1337] * v[1411]);
	v[1423] = (*epst)*(v[1338] * v[1409] + v[1339] * v[1410] + v[1340] * v[1411]);
	v[1427] = (*epst)*(v[1409] * v[1424] + v[1410] * v[1425] + v[1411] * v[1426]);
	v[1431] = (*epst)*(v[1409] * v[1428] + v[1410] * v[1429] + v[1411] * v[1430]);
	v[1435] = (*epst)*(v[1409] * v[1432] + v[1410] * v[1433] + v[1411] * v[1434]);
	v[1437] = (*epst)*(v[1413] * v[1417] + v[1414] * v[1418] + v[1415] * v[1419]);
	v[1438] = (*epst)*(v[1332] * v[1413] + v[1333] * v[1414] + v[1334] * v[1415]);
	v[1439] = (*epst)*(v[1335] * v[1413] + v[1336] * v[1414] + v[1337] * v[1415]);
	v[1440] = (*epst)*(v[1338] * v[1413] + v[1339] * v[1414] + v[1340] * v[1415]);
	v[1441] = (*epst)*(v[1413] * v[1424] + v[1414] * v[1425] + v[1415] * v[1426]);
	v[1442] = (*epst)*(v[1413] * v[1428] + v[1414] * v[1429] + v[1415] * v[1430]);
	v[1443] = (*epst)*(v[1413] * v[1432] + v[1414] * v[1433] + v[1415] * v[1434]);
	v[1445] = (*epst)*(v[1332] * v[1417] + v[1333] * v[1418] + v[1334] * v[1419]);
	v[1446] = (*epst)*(v[1335] * v[1417] + v[1336] * v[1418] + v[1337] * v[1419]);
	v[1447] = (*epst)*(v[1338] * v[1417] + v[1339] * v[1418] + v[1340] * v[1419]);
	v[1448] = (*epst)*(v[1417] * v[1424] + v[1418] * v[1425] + v[1419] * v[1426]);
	v[1449] = (*epst)*(v[1417] * v[1428] + v[1418] * v[1429] + v[1419] * v[1430]);
	v[1450] = (*epst)*(v[1417] * v[1432] + v[1418] * v[1433] + v[1419] * v[1434]);
	v[1452] = (*epst)*(v[1332] * v[1335] + v[1333] * v[1336] + v[1334] * v[1337]);
	v[1453] = (*epst)*(v[1332] * v[1338] + v[1333] * v[1339] + v[1334] * v[1340]);
	v[1454] = (*epst)*(v[1332] * v[1424] + v[1333] * v[1425] + v[1334] * v[1426]);
	v[1455] = (*epst)*(v[1332] * v[1428] + v[1333] * v[1429] + v[1334] * v[1430]);
	v[1456] = (*epst)*(v[1332] * v[1432] + v[1333] * v[1433] + v[1334] * v[1434]);
	v[1458] = (*epst)*(v[1335] * v[1338] + v[1336] * v[1339] + v[1337] * v[1340]);
	v[1459] = (*epst)*(v[1335] * v[1424] + v[1336] * v[1425] + v[1337] * v[1426]);
	v[1460] = (*epst)*(v[1335] * v[1428] + v[1336] * v[1429] + v[1337] * v[1430]);
	v[1461] = (*epst)*(v[1335] * v[1432] + v[1336] * v[1433] + v[1337] * v[1434]);
	v[1463] = (*epst)*(v[1338] * v[1424] + v[1339] * v[1425] + v[1340] * v[1426]);
	v[1464] = (*epst)*(v[1338] * v[1428] + v[1339] * v[1429] + v[1340] * v[1430]);
	v[1465] = (*epst)*(v[1338] * v[1432] + v[1339] * v[1433] + v[1340] * v[1434]);
	v[1467] = (*epst)*(v[1424] * v[1428] + v[1425] * v[1429] + v[1426] * v[1430]);
	v[1468] = (*epst)*(v[1424] * v[1432] + v[1425] * v[1433] + v[1426] * v[1434]);
	v[1470] = (*epst)*(v[1428] * v[1432] + v[1429] * v[1433] + v[1430] * v[1434]);
	v[3033] = v[1183] * v[3824] + v[1184] * v[3906] + v[1185] * v[3993] + v[4777] * v[700] + v[4778] * v[706];
	v[4804] = v[1416] + v[3033];
	v[3038] = v[1183] * v[3825] + v[1184] * v[3907] + v[1185] * v[3994] + v[4777] * v[701] + v[4778] * v[707];
	v[4807] = v[1420] + v[3038];
	v[3039] = v[1183] * v[3826] + v[1184] * v[3908] + v[1185] * v[3995];
	v[4811] = v[1421] + v[3039];
	v[3040] = v[1183] * v[3827] + v[1184] * v[3909] + v[1185] * v[3996];
	v[4814] = v[1422] + v[3040];
	v[3041] = v[1183] * v[3828] + v[1184] * v[3910] + v[1185] * v[3997];
	v[4817] = v[1423] + v[3041];
	v[3042] = v[1183] * v[3830] + v[1184] * v[3912] + v[1185] * v[3999] + v[4777] * v[702] + v[4778] * v[708];
	v[4820] = v[1427] + v[3042];
	v[3043] = v[1183] * v[3832] + v[1184] * v[3914] + v[1185] * v[4001] + v[4777] * v[703] + v[4778] * v[709];
	v[4825] = v[1431] + v[3043];
	v[3044] = v[1183] * v[3834] + v[1184] * v[3916] + v[1185] * v[4003] + v[4777] * v[704] + v[4778] * v[710];
	v[4830] = v[1435] + v[3044];
	v[3052] = v[1183] * v[3813] + v[1184] * v[3895] + v[1185] * v[3982] + v[4779] * v[701] + v[4780] * v[707];
	v[4808] = v[1437] + v[3052];
	v[3053] = v[1183] * v[3814] + v[1184] * v[3896] + v[1185] * v[3983];
	v[4812] = v[1438] + v[3053];
	v[3054] = v[1183] * v[3815] + v[1184] * v[3897] + v[1185] * v[3984];
	v[4815] = v[1439] + v[3054];
	v[3055] = v[1183] * v[3816] + v[1184] * v[3898] + v[1185] * v[3985];
	v[4818] = v[1440] + v[3055];
	v[3056] = v[1183] * v[3818] + v[1184] * v[3900] + v[1185] * v[3987] + v[4779] * v[702] + v[4780] * v[708];
	v[4821] = v[1441] + v[3056];
	v[3057] = v[1183] * v[3820] + v[1184] * v[3902] + v[1185] * v[3989] + v[4779] * v[703] + v[4780] * v[709];
	v[4826] = v[1442] + v[3057];
	v[3058] = v[1183] * v[3822] + v[1184] * v[3904] + v[1185] * v[3991] + v[4779] * v[704] + v[4780] * v[710];
	v[4831] = v[1443] + v[3058];
	v[3060] = v[1183] * v[3802] + v[1184] * v[3882] + v[1185] * v[3969];
	v[4813] = v[1445] + v[3060];
	v[3061] = v[1183] * v[3803] + v[1184] * v[3883] + v[1185] * v[3970];
	v[4816] = v[1446] + v[3061];
	v[3062] = v[1183] * v[3804] + v[1184] * v[3884] + v[1185] * v[3971];
	v[4819] = v[1447] + v[3062];
	v[3067] = v[1183] * v[3806] + v[1184] * v[3886] + v[1185] * v[3973] + v[4781] * v[702] + v[4782] * v[708];
	v[4824] = v[1448] + v[3067];
	v[3068] = v[1183] * v[3808] + v[1184] * v[3888] + v[1185] * v[3975] + v[4781] * v[703] + v[4782] * v[709];
	v[4829] = v[1449] + v[3068];
	v[3069] = v[1183] * v[3810] + v[1184] * v[3890] + v[1185] * v[3977] + v[4781] * v[704] + v[4782] * v[710];
	v[4834] = v[1450] + v[3069];
	v[3071] = v[1183] * v[3716] + v[1184] * v[3722] + v[1185] * v[3728];
	v[3072] = v[1183] * v[3717] + v[1184] * v[3723] + v[1185] * v[3729];
	v[3073] = v[1183] * v[3718] + v[1184] * v[3724] + v[1185] * v[3730] + v[3075] * v[702] + v[3074] * v[708];
	v[4837] = v[1454] + v[3073];
	v[3076] = v[1183] * v[3719] + v[1184] * v[3725] + v[1185] * v[3731] + v[3075] * v[703] + v[3074] * v[709];
	v[4838] = v[1455] + v[3076];
	v[3077] = v[1183] * v[3720] + v[1184] * v[3726] + v[1185] * v[3732] + v[3075] * v[704] + v[3074] * v[710];
	v[4839] = v[1456] + v[3077];
	v[3079] = v[1183] * v[3734] + v[1184] * v[3739] + v[1185] * v[3744];
	v[3080] = v[1183] * v[3735] + v[1184] * v[3740] + v[1185] * v[3745] + v[3082] * v[702] + v[3081] * v[708];
	v[4841] = v[1459] + v[3080];
	v[3083] = v[1183] * v[3736] + v[1184] * v[3741] + v[1185] * v[3746] + v[3082] * v[703] + v[3081] * v[709];
	v[4842] = v[1460] + v[3083];
	v[3084] = v[1183] * v[3737] + v[1184] * v[3742] + v[1185] * v[3747] + v[3082] * v[704] + v[3081] * v[710];
	v[4843] = v[1461] + v[3084];
	v[3086] = v[1183] * v[3749] + v[1184] * v[3753] + v[1185] * v[3757] + v[3088] * v[702] + v[3087] * v[708];
	v[4844] = v[1463] + v[3086];
	v[3089] = v[1183] * v[3750] + v[1184] * v[3754] + v[1185] * v[3758] + v[3088] * v[703] + v[3087] * v[709];
	v[4845] = v[1464] + v[3089];
	v[3090] = v[1183] * v[3751] + v[1184] * v[3755] + v[1185] * v[3759] + v[3088] * v[704] + v[3087] * v[710];
	v[4846] = v[1465] + v[3090];
	v[3096] = v[1183] * v[3796] + v[1184] * v[3874] + v[1185] * v[3959] + v[4784] * v[703] + v[4783] * v[709];
	v[4847] = v[1467] + v[3096];
	v[3097] = v[1183] * v[3798] + v[1184] * v[3876] + v[1185] * v[3961] + v[4784] * v[704] + v[4783] * v[710];
	v[4848] = v[1468] + v[3097];
	v[3103] = v[1183] * v[3785] + v[1184] * v[3862] + v[1185] * v[3946] + v[4849] * v[704] + v[4850] * v[710];
	v[4851] = v[1470] + v[3103];
	v[3117] = v[257] * v[4788];
	v[3124] = v[253] * v[4790];
	v[3113] = v[3112] - v[903];
	v[4789] = v[3113] / v[3112];
	v[3453] = v[2352] / v[3112];
	v[3435] = v[3453] * v[4508];
	v[3431] = v[2351] / v[3112];
	v[3461] = v[3431] * v[4507];
	v[3413] = v[2350] / v[3112];
	v[3460] = v[3413] * v[4506];
	v[3412] = v[3124] + v[258] * v[4788] + v[248] * v[4789];
	v[3414] = (*a4)*v[3412] + v[3454] + v[3413] * v[4508];
	v[3409] = v[3117] + v[247] * v[4789] + v[252] * v[4790];
	v[3411] = (*a4)*v[3409] + v[3432] + v[3413] * v[4507];
	v[3408] = v[3435] + v[3461] + v[4482] * v[4506] + (*a4)*v[4789];
	v[3114] = v[3409] * v[4507] + v[3412] * v[4508] + v[4506] * v[4789];
	v[3116] = v[3112] - v[900];
	v[4791] = v[3116] / v[3112];
	v[3433] = -v[3432] + v[3431] * v[4506] + (*a4)*v[4790];
	v[3127] = -((v[248] * v[4497]) / v[3112]);
	v[3437] = v[3127] - (v[258] * v[4494]) / v[3112] + v[253] * v[4791];
	v[3438] = (*a4)*v[3437] + v[3457] + v[3431] * v[4508];
	v[3436] = v[3435] + v[3460] + v[4486] * v[4507] + (*a4)*v[4791];
	v[3122] = v[3437] * v[4508] + v[4506] * v[4790] + v[4507] * v[4791];
	v[3123] = v[3108] + v[3112];
	v[3455] = -v[3454] + v[3453] * v[4506] + (*a4)*v[4788];
	v[3459] = v[3123] / v[3112];
	v[3462] = (*a4)*v[3459] + v[3460] + v[3461] + v[4489] * v[4508];
	v[3456] = v[3127] + v[257] * v[3459] - (v[252] * v[4498]) / v[3112];
	v[3458] = (*a4)*v[3456] - v[3457] + v[3453] * v[4507];
	v[3128] = v[3456] * v[4507] + v[3459] * v[4508] + v[4506] * v[4788];
	v[4102] = (*a4)*v[1419] + v[3128] * v[3921] + v[3122] * v[3936] + v[3114] * v[3951] + v[3994] * v[4795]
		+ v[3982] * v[4796] + v[3968] * v[4797];
	v[4101] = (*a4)*v[1415] + v[3128] * v[3920] + v[3122] * v[3935] + v[3114] * v[3950] + v[3993] * v[4795]
		+ v[3981] * v[4796] + v[3967] * v[4797];
	v[4100] = (*a4)*v[1411] + v[3128] * v[3918] + v[3122] * v[3933] + v[3114] * v[3948] + v[3992] * v[4795]
		+ v[3980] * v[4796] + v[3964] * v[4797];
	v[4093] = (*a4)*v[1418] + v[3128] * v[3838] + v[3122] * v[3852] + v[3114] * v[3866] + v[3907] * v[4795]
		+ v[3895] * v[4796] + v[3881] * v[4797];
	v[4092] = (*a4)*v[1414] + v[3128] * v[3837] + v[3122] * v[3851] + v[3114] * v[3865] + v[3906] * v[4795]
		+ v[3894] * v[4796] + v[3880] * v[4797];
	v[4091] = (*a4)*v[1410] + v[3128] * v[3836] + v[3122] * v[3850] + v[3114] * v[3864] + v[3905] * v[4795]
		+ v[3893] * v[4796] + v[3879] * v[4797];
	v[4084] = (*a4)*v[1417] + v[3128] * v[3762] + v[3122] * v[3775] + v[3114] * v[3788] + v[3825] * v[4795]
		+ v[3813] * v[4796] + v[3801] * v[4797];
	v[4083] = (*a4)*v[1413] + v[3128] * v[3761] + v[3122] * v[3774] + v[3114] * v[3787] + v[3824] * v[4795]
		+ v[3812] * v[4796] + v[3800] * v[4797];
	v[4082] = (*a4)*v[1409] + v[3128] * v[3760] + v[3122] * v[3773] + v[3114] * v[3786] + v[3823] * v[4795]
		+ v[3811] * v[4796] + v[3799] * v[4797];
	v[3653] = v[3114] * v[3631] + v[3122] * v[3637] + v[3128] * v[3643] + v[3613] * v[4795] + v[3619] * v[4796]
		+ v[3625] * v[4797];
	v[3652] = v[3114] * v[3630] + v[3122] * v[3636] + v[3128] * v[3642] + v[3612] * v[4795] + v[3618] * v[4796]
		+ v[3624] * v[4797];
	v[3651] = v[3114] * v[3633] + v[3122] * v[3639] + v[3128] * v[3645] + v[3615] * v[4795] + v[3621] * v[4796]
		+ v[3627] * v[4797];
	v[3650] = v[3114] * v[3632] + v[3122] * v[3638] + v[3128] * v[3644] + v[3614] * v[4795] + v[3620] * v[4796]
		+ v[3626] * v[4797];
	v[3649] = v[3114] * v[3635] + v[3122] * v[3641] + v[3128] * v[3647] + v[3617] * v[4795] + v[3623] * v[4796]
		+ v[3629] * v[4797];
	v[3648] = v[3114] * v[3634] + v[3122] * v[3640] + v[3128] * v[3646] + v[3616] * v[4795] + v[3622] * v[4796]
		+ v[3628] * v[4797];
	v[3480] = v[3114] * v[3321] + v[3122] * v[3347] + v[3128] * v[3370] + v[3170] * v[3414] + v[3171] * v[3438]
		+ v[3172] * v[3462] + v[3264] * v[4795] + v[3282] * v[4796] + v[3300] * v[4797];
	v[3479] = v[3114] * v[3319] + v[3122] * v[3345] + v[3128] * v[3368] + v[3170] * v[3411] + v[3171] * v[3436]
		+ v[3172] * v[3458] + v[3263] * v[4795] + v[3281] * v[4796] + v[3299] * v[4797];
	v[3478] = v[3114] * v[3317] + v[3122] * v[3343] + v[3128] * v[3367] + v[3170] * v[3408] + v[3171] * v[3433]
		+ v[3172] * v[3455] + v[3262] * v[4795] + v[3280] * v[4796] + v[3298] * v[4797];
	v[3477] = (*a4)*v[3169] + v[3114] * v[3315] + v[3122] * v[3342] + v[3128] * v[3366] + v[3261] * v[4795]
		+ v[3279] * v[4796] + v[3297] * v[4797];
	v[3476] = (*a4)*v[3168] + v[3114] * v[3314] + v[3122] * v[3341] + v[3128] * v[3365] + v[3260] * v[4795]
		+ v[3278] * v[4796] + v[3296] * v[4797];
	v[3475] = (*a4)*v[3167] + v[3114] * v[3313] + v[3122] * v[3340] + v[3128] * v[3364] + v[3259] * v[4795]
		+ v[3277] * v[4796] + v[3295] * v[4797];
	v[3474] = v[3114] * v[3330] + v[3122] * v[3355] + v[3128] * v[3377] + v[3163] * v[3414] + v[3164] * v[3438]
		+ v[3165] * v[3462] + v[3270] * v[4795] + v[3288] * v[4796] + v[3306] * v[4797];
	v[3473] = v[3114] * v[3328] + v[3122] * v[3353] + v[3128] * v[3375] + v[3163] * v[3411] + v[3164] * v[3436]
		+ v[3165] * v[3458] + v[3269] * v[4795] + v[3287] * v[4796] + v[3305] * v[4797];
	v[3472] = v[3114] * v[3326] + v[3122] * v[3351] + v[3128] * v[3374] + v[3163] * v[3408] + v[3164] * v[3433]
		+ v[3165] * v[3455] + v[3268] * v[4795] + v[3286] * v[4796] + v[3304] * v[4797];
	v[3471] = (*a4)*v[3162] + v[3114] * v[3324] + v[3122] * v[3350] + v[3128] * v[3373] + v[3267] * v[4795]
		+ v[3285] * v[4796] + v[3303] * v[4797];
	v[3470] = (*a4)*v[3161] + v[3114] * v[3323] + v[3122] * v[3349] + v[3128] * v[3372] + v[3266] * v[4795]
		+ v[3284] * v[4796] + v[3302] * v[4797];
	v[3469] = (*a4)*v[3160] + v[3114] * v[3322] + v[3122] * v[3348] + v[3128] * v[3371] + v[3265] * v[4795]
		+ v[3283] * v[4796] + v[3301] * v[4797];
	v[3468] = v[3114] * v[3339] + v[3122] * v[3363] + v[3128] * v[3384] + v[3156] * v[3414] + v[3157] * v[3438]
		+ v[3158] * v[3462] + v[3276] * v[4795] + v[3294] * v[4796] + v[3312] * v[4797];
	v[3467] = v[3114] * v[3337] + v[3122] * v[3361] + v[3128] * v[3382] + v[3156] * v[3411] + v[3157] * v[3436]
		+ v[3158] * v[3458] + v[3275] * v[4795] + v[3293] * v[4796] + v[3311] * v[4797];
	v[3466] = v[3114] * v[3335] + v[3122] * v[3359] + v[3128] * v[3381] + v[3156] * v[3408] + v[3157] * v[3433]
		+ v[3158] * v[3455] + v[3274] * v[4795] + v[3292] * v[4796] + v[3310] * v[4797];
	v[3465] = (*a4)*v[3155] + v[3114] * v[3333] + v[3122] * v[3358] + v[3128] * v[3380] + v[3273] * v[4795]
		+ v[3291] * v[4796] + v[3309] * v[4797];
	v[3464] = (*a4)*v[3154] + v[3114] * v[3332] + v[3122] * v[3357] + v[3128] * v[3379] + v[3272] * v[4795]
		+ v[3290] * v[4796] + v[3308] * v[4797];
	v[3463] = (*a4)*v[3153] + v[3114] * v[3331] + v[3122] * v[3356] + v[3128] * v[3378] + v[3271] * v[4795]
		+ v[3289] * v[4796] + v[3307] * v[4797];
	v[3141] = v[276] * v[4798];
	v[3135] = -v[1196] + v[2756];
	v[4799] = v[3135] / v[3136];
	v[3148] = v[272] * v[4799];
	v[3137] = v[3132] + v[3136];
	v[4800] = v[3137] / v[3136];
	v[4072] = v[4479] / v[3136];
	v[4055] = v[4072] * v[4505];
	v[4051] = v[4481] / v[3136];
	v[4080] = v[4051] * v[4504];
	v[4032] = v[4527] / v[3136];
	v[4079] = v[4032] * v[4503];
	v[4031] = v[3148] + v[277] * v[4798] + v[267] * v[4800];
	v[4033] = (*a4)*v[4031] + v[4073] + v[4032] * v[4505];
	v[4028] = v[3141] + v[271] * v[4799] + v[266] * v[4800];
	v[4030] = (*a4)*v[4028] + v[4052] + v[4032] * v[4504];
	v[4027] = v[4055] + v[4080] + v[4467] * v[4503] + (*a4)*v[4800];
	v[3138] = v[4028] * v[4504] + v[4031] * v[4505] + v[4503] * v[4800];
	v[3140] = -v[2820] + v[3136];
	v[4801] = v[3140] / v[3136];
	v[4053] = -v[4052] + v[4051] * v[4503] + (*a4)*v[4799];
	v[3151] = (v[267] * v[4547]) / v[3136];
	v[3144] = -v[1193] - v[2760];
	v[4057] = (v[277] * v[3144]) / v[3136] + v[3151] + v[272] * v[4801];
	v[4058] = (*a4)*v[4057] + v[4076] + v[4051] * v[4505];
	v[4056] = v[4055] + v[4079] + v[4471] * v[4504] + (*a4)*v[4801];
	v[3146] = v[4057] * v[4505] + v[4503] * v[4799] + v[4504] * v[4801];
	v[3147] = v[2780] + v[3136];
	v[4074] = -v[4073] + v[4072] * v[4503] + (*a4)*v[4798];
	v[4078] = v[3147] / v[3136];
	v[4081] = (*a4)*v[4078] + v[4079] + v[4080] + v[4474] * v[4505];
	v[4075] = v[3151] + v[276] * v[4078] + (v[271] * v[4546]) / v[3136];
	v[4077] = (*a4)*v[4075] - v[4076] + v[4072] * v[4504];
	v[3152] = v[4075] * v[4504] + v[4078] * v[4505] + v[4503] * v[4798];
	v[4302] = v[2728] * v[3138] + v[2729] * v[3146] + v[2730] * v[3152] + v[3128] * v[4276] + v[3122] * v[4280]
		+ v[3114] * v[4284] + v[4296] * v[4795] + v[4292] * v[4796] + v[4288] * v[4797];
	v[4301] = v[2704] * v[3138] + v[2705] * v[3146] + v[2706] * v[3152] + v[3128] * v[4274] + v[3122] * v[4278]
		+ v[3114] * v[4282] + v[4294] * v[4795] + v[4290] * v[4796] + v[4286] * v[4797];
	v[4300] = v[2721] * v[3138] + v[2722] * v[3146] + v[2723] * v[3152] + v[3128] * v[4252] + v[3122] * v[4256]
		+ v[3114] * v[4260] + v[4272] * v[4795] + v[4268] * v[4796] + v[4264] * v[4797];
	v[4299] = v[2697] * v[3138] + v[2698] * v[3146] + v[2699] * v[3152] + v[3128] * v[4250] + v[3122] * v[4254]
		+ v[3114] * v[4258] + v[4270] * v[4795] + v[4266] * v[4796] + v[4262] * v[4797];
	v[4298] = v[2713] * v[3138] + v[2714] * v[3146] + v[2715] * v[3152] + v[3128] * v[4228] + v[3122] * v[4232]
		+ v[3114] * v[4236] + v[4248] * v[4795] + v[4244] * v[4796] + v[4240] * v[4797];
	v[4297] = v[2689] * v[3138] + v[2690] * v[3146] + v[2691] * v[3152] + v[3128] * v[4226] + v[3122] * v[4230]
		+ v[3114] * v[4234] + v[4246] * v[4795] + v[4242] * v[4796] + v[4238] * v[4797];
	v[4108] = v[1426] * v[3414] + v[1430] * v[3438] + v[1434] * v[3462] + v[3138] * v[3732] + v[3146] * v[3747]
		+ v[3152] * v[3759] + v[3128] * v[3931] + v[3122] * v[3946] + v[3114] * v[3961] + v[4003] * v[4795] + v[3991] * v[4796]
		+ v[3977] * v[4797];
	v[4107] = v[1426] * v[3411] + v[1430] * v[3436] + v[1434] * v[3458] + v[3138] * v[3731] + v[3146] * v[3746]
		+ v[3152] * v[3758] + v[3128] * v[3928] + v[3122] * v[3944] + v[3114] * v[3959] + v[4001] * v[4795] + v[3989] * v[4796]
		+ v[3975] * v[4797];
	v[4106] = v[1426] * v[3408] + v[1430] * v[3433] + v[1434] * v[3455] + v[3138] * v[3730] + v[3146] * v[3745]
		+ v[3152] * v[3757] + v[3128] * v[3926] + v[3122] * v[3941] + v[3114] * v[3957] + v[3999] * v[4795] + v[3987] * v[4796]
		+ v[3973] * v[4797];
	v[4105] = v[3138] * v[3729] + v[3146] * v[3744] + v[3152] * v[3756] + v[3128] * v[3924] + v[3122] * v[3939]
		+ v[3114] * v[3954] + v[1334] * v[4033] + v[1337] * v[4058] + v[1340] * v[4081] + v[3997] * v[4795] + v[3985] * v[4796]
		+ v[3971] * v[4797];
	v[4104] = v[3138] * v[3728] + v[3146] * v[3743] + v[3152] * v[3744] + v[3128] * v[3923] + v[3122] * v[3938]
		+ v[3114] * v[3953] + v[1334] * v[4030] + v[1337] * v[4056] + v[1340] * v[4077] + v[3996] * v[4795] + v[3984] * v[4796]
		+ v[3970] * v[4797];
	v[4103] = v[3138] * v[3727] + v[3146] * v[3728] + v[3152] * v[3729] + v[3128] * v[3922] + v[3122] * v[3937]
		+ v[3114] * v[3952] + v[1334] * v[4027] + v[1337] * v[4053] + v[1340] * v[4074] + v[3995] * v[4795] + v[3983] * v[4796]
		+ v[3969] * v[4797];
	v[4099] = v[1425] * v[3414] + v[1429] * v[3438] + v[1433] * v[3462] + v[3138] * v[3726] + v[3146] * v[3742]
		+ v[3152] * v[3755] + v[3128] * v[3848] + v[3122] * v[3862] + v[3114] * v[3876] + v[3916] * v[4795] + v[3904] * v[4796]
		+ v[3890] * v[4797];
	v[4098] = v[1425] * v[3411] + v[1429] * v[3436] + v[1433] * v[3458] + v[3138] * v[3725] + v[3146] * v[3741]
		+ v[3152] * v[3754] + v[3128] * v[3845] + v[3122] * v[3860] + v[3114] * v[3874] + v[3914] * v[4795] + v[3902] * v[4796]
		+ v[3888] * v[4797];
	v[4097] = v[1425] * v[3408] + v[1429] * v[3433] + v[1433] * v[3455] + v[3138] * v[3724] + v[3146] * v[3740]
		+ v[3152] * v[3753] + v[3128] * v[3843] + v[3122] * v[3857] + v[3114] * v[3872] + v[3912] * v[4795] + v[3900] * v[4796]
		+ v[3886] * v[4797];
	v[4096] = v[3138] * v[3723] + v[3146] * v[3739] + v[3152] * v[3752] + v[3128] * v[3841] + v[3122] * v[3855]
		+ v[3114] * v[3869] + v[1333] * v[4033] + v[1336] * v[4058] + v[1339] * v[4081] + v[3910] * v[4795] + v[3898] * v[4796]
		+ v[3884] * v[4797];
	v[4095] = v[3138] * v[3722] + v[3146] * v[3738] + v[3152] * v[3739] + v[3128] * v[3840] + v[3122] * v[3854]
		+ v[3114] * v[3868] + v[1333] * v[4030] + v[1336] * v[4056] + v[1339] * v[4077] + v[3909] * v[4795] + v[3897] * v[4796]
		+ v[3883] * v[4797];
	v[4094] = v[3138] * v[3721] + v[3146] * v[3722] + v[3152] * v[3723] + v[3128] * v[3839] + v[3122] * v[3853]
		+ v[3114] * v[3867] + v[1333] * v[4027] + v[1336] * v[4053] + v[1339] * v[4074] + v[3908] * v[4795] + v[3896] * v[4796]
		+ v[3882] * v[4797];
	v[4090] = v[1424] * v[3414] + v[1428] * v[3438] + v[1432] * v[3462] + v[3138] * v[3720] + v[3146] * v[3737]
		+ v[3152] * v[3751] + v[3128] * v[3772] + v[3122] * v[3785] + v[3114] * v[3798] + v[3834] * v[4795] + v[3822] * v[4796]
		+ v[3810] * v[4797];
	v[4089] = v[1424] * v[3411] + v[1428] * v[3436] + v[1432] * v[3458] + v[3138] * v[3719] + v[3146] * v[3736]
		+ v[3152] * v[3750] + v[3128] * v[3769] + v[3122] * v[3783] + v[3114] * v[3796] + v[3832] * v[4795] + v[3820] * v[4796]
		+ v[3808] * v[4797];
	v[4088] = v[1424] * v[3408] + v[1428] * v[3433] + v[1432] * v[3455] + v[3138] * v[3718] + v[3146] * v[3735]
		+ v[3152] * v[3749] + v[3128] * v[3767] + v[3122] * v[3780] + v[3114] * v[3794] + v[3830] * v[4795] + v[3818] * v[4796]
		+ v[3806] * v[4797];
	v[4087] = v[3138] * v[3717] + v[3146] * v[3734] + v[3152] * v[3748] + v[3128] * v[3765] + v[3122] * v[3778]
		+ v[3114] * v[3791] + v[1332] * v[4033] + v[1335] * v[4058] + v[1338] * v[4081] + v[3828] * v[4795] + v[3816] * v[4796]
		+ v[3804] * v[4797];
	v[4086] = v[3138] * v[3716] + v[3146] * v[3733] + v[3152] * v[3734] + v[3128] * v[3764] + v[3122] * v[3777]
		+ v[3114] * v[3790] + v[1332] * v[4030] + v[1335] * v[4056] + v[1338] * v[4077] + v[3827] * v[4795] + v[3815] * v[4796]
		+ v[3803] * v[4797];
	v[4085] = v[3138] * v[3715] + v[3146] * v[3716] + v[3152] * v[3717] + v[3128] * v[3763] + v[3122] * v[3776]
		+ v[3114] * v[3789] + v[1332] * v[4027] + v[1335] * v[4053] + v[1338] * v[4074] + v[3826] * v[4795] + v[3814] * v[4796]
		+ v[3802] * v[4797];
	v[3159] = v[3114] * v[3156] + v[3122] * v[3157] + v[3128] * v[3158] + v[3153] * v[4795] + v[3154] * v[4796]
		+ v[3155] * v[4797];
	v[3166] = v[3114] * v[3163] + v[3122] * v[3164] + v[3128] * v[3165] + v[3160] * v[4795] + v[3161] * v[4796]
		+ v[3162] * v[4797];
	v[3173] = v[3114] * v[3170] + v[3122] * v[3171] + v[3128] * v[3172] + v[3167] * v[4795] + v[3168] * v[4796]
		+ v[3169] * v[4797];
	v[3659] = (*cn)*(v[3173] * v[3625] + v[3166] * v[3627] + v[3159] * v[3629] + v[3155] * v[3649] + v[3162] * v[3651]
		+ v[3169] * v[3653]);
	v[3658] = (*cn)*(v[3173] * v[3624] + v[3166] * v[3626] + v[3159] * v[3628] + v[3155] * v[3648] + v[3162] * v[3650]
		+ v[3169] * v[3652]);
	v[3657] = (*cn)*(v[3173] * v[3619] + v[3166] * v[3621] + v[3159] * v[3623] + v[3154] * v[3649] + v[3161] * v[3651]
		+ v[3168] * v[3653]);
	v[3656] = (*cn)*(v[3173] * v[3618] + v[3166] * v[3620] + v[3159] * v[3622] + v[3154] * v[3648] + v[3161] * v[3650]
		+ v[3168] * v[3652]);
	v[3655] = (*cn)*(v[3173] * v[3613] + v[3166] * v[3615] + v[3159] * v[3617] + v[3153] * v[3649] + v[3160] * v[3651]
		+ v[3167] * v[3653]);
	v[3654] = (*cn)*(v[3173] * v[3612] + v[3166] * v[3614] + v[3159] * v[3616] + v[3153] * v[3648] + v[3160] * v[3650]
		+ v[3167] * v[3652]);
	v[3174] = v[1424] * v[3114] + v[1428] * v[3122] + v[1432] * v[3128] + v[1332] * v[3138] + v[1335] * v[3146]
		+ v[1338] * v[3152] + v[1409] * v[4795] + v[1413] * v[4796] + v[1417] * v[4797];
	v[3175] = v[1425] * v[3114] + v[1429] * v[3122] + v[1433] * v[3128] + v[1333] * v[3138] + v[1336] * v[3146]
		+ v[1339] * v[3152] + v[1410] * v[4795] + v[1414] * v[4796] + v[1418] * v[4797];
	v[3176] = v[1426] * v[3114] + v[1430] * v[3122] + v[1434] * v[3128] + v[1334] * v[3138] + v[1337] * v[3146]
		+ v[1340] * v[3152] + v[1411] * v[4795] + v[1415] * v[4796] + v[1419] * v[4797];
	v[4308] = (*ct)*(v[3174] * v[4240] + v[3175] * v[4264] + v[3176] * v[4288] + v[1417] * v[4298] + v[1418] * v[4300]
		+ v[1419] * v[4302]);
	v[4810] = v[3659] + v[4308];
	v[4307] = (*ct)*(v[3174] * v[4238] + v[3175] * v[4262] + v[3176] * v[4286] + v[1417] * v[4297] + v[1418] * v[4299]
		+ v[1419] * v[4301]);
	v[4809] = v[3658] + v[4307];
	v[4306] = (*ct)*(v[3174] * v[4244] + v[3175] * v[4268] + v[3176] * v[4292] + v[1413] * v[4298] + v[1414] * v[4300]
		+ v[1415] * v[4302]);
	v[4806] = v[3657] + v[4306];
	v[4305] = (*ct)*(v[3174] * v[4242] + v[3175] * v[4266] + v[3176] * v[4290] + v[1413] * v[4297] + v[1414] * v[4299]
		+ v[1415] * v[4301]);
	v[4805] = v[3656] + v[4305];
	v[4304] = (*ct)*(v[3174] * v[4248] + v[3175] * v[4272] + v[3176] * v[4296] + v[1409] * v[4298] + v[1410] * v[4300]
		+ v[1411] * v[4302]);
	v[4803] = v[3655] + v[4304];
	v[4303] = (*ct)*(v[3174] * v[4246] + v[3175] * v[4270] + v[3176] * v[4294] + v[1409] * v[4297] + v[1410] * v[4299]
		+ v[1411] * v[4301]);
	v[4802] = v[3654] + v[4303];
	v[4169] = (*ct)*(v[3174] * v[3734] + v[3175] * v[3739] + v[3176] * v[3744]);
	v[4840] = v[1458] + v[3079] + v[4169];
	v[4167] = (*ct)*(v[3174] * v[3717] + v[3175] * v[3723] + v[3176] * v[3729]);
	v[4836] = v[1453] + v[3072] + v[4167];
	v[4160] = (*ct)*(v[3174] * v[3716] + v[3175] * v[3722] + v[3176] * v[3728]);
	v[4835] = v[1452] + v[3071] + v[4160];
	v[3660] = (*cn)*(v[3173] * v[3630] + v[3166] * v[3632] + v[3159] * v[3634] + v[3156] * v[3648] + v[3163] * v[3650]
		+ v[3170] * v[3652]);
	v[3661] = (*cn)*(v[3173] * v[3636] + v[3166] * v[3638] + v[3159] * v[3640] + v[3157] * v[3648] + v[3164] * v[3650]
		+ v[3171] * v[3652]);
	v[3662] = (*cn)*(v[3173] * v[3642] + v[3166] * v[3644] + v[3159] * v[3646] + v[3158] * v[3648] + v[3165] * v[3650]
		+ v[3172] * v[3652]);
	v[3663] = (*cn)*(v[3173] * v[3631] + v[3166] * v[3633] + v[3159] * v[3635] + v[3156] * v[3649] + v[3163] * v[3651]
		+ v[3170] * v[3653]);
	v[3664] = (*cn)*(v[3173] * v[3637] + v[3166] * v[3639] + v[3159] * v[3641] + v[3157] * v[3649] + v[3164] * v[3651]
		+ v[3171] * v[3653]);
	v[3665] = (*cn)*(v[3173] * v[3643] + v[3166] * v[3645] + v[3159] * v[3647] + v[3158] * v[3649] + v[3165] * v[3651]
		+ v[3172] * v[3653]);
	v[4309] = (*ct)*(v[2689] * v[3174] + v[2697] * v[3175] + v[2704] * v[3176] + v[1332] * v[4297] + v[1333] * v[4299]
		+ v[1334] * v[4301]);
	v[4310] = (*ct)*(v[2690] * v[3174] + v[2698] * v[3175] + v[2705] * v[3176] + v[1335] * v[4297] + v[1336] * v[4299]
		+ v[1337] * v[4301]);
	v[4311] = (*ct)*(v[2691] * v[3174] + v[2699] * v[3175] + v[2706] * v[3176] + v[1338] * v[4297] + v[1339] * v[4299]
		+ v[1340] * v[4301]);
	v[4312] = (*ct)*(v[3174] * v[4234] + v[3175] * v[4258] + v[3176] * v[4282] + v[1424] * v[4297] + v[1425] * v[4299]
		+ v[1426] * v[4301]);
	v[4822] = v[3660] + v[4312];
	v[4313] = (*ct)*(v[3174] * v[4230] + v[3175] * v[4254] + v[3176] * v[4278] + v[1428] * v[4297] + v[1429] * v[4299]
		+ v[1430] * v[4301]);
	v[4827] = v[3661] + v[4313];
	v[4314] = (*ct)*(v[3174] * v[4226] + v[3175] * v[4250] + v[3176] * v[4274] + v[1432] * v[4297] + v[1433] * v[4299]
		+ v[1434] * v[4301]);
	v[4832] = v[3662] + v[4314];
	v[4315] = (*ct)*(v[2713] * v[3174] + v[2721] * v[3175] + v[2728] * v[3176] + v[1332] * v[4298] + v[1333] * v[4300]
		+ v[1334] * v[4302]);
	v[4316] = (*ct)*(v[2714] * v[3174] + v[2722] * v[3175] + v[2729] * v[3176] + v[1335] * v[4298] + v[1336] * v[4300]
		+ v[1337] * v[4302]);
	v[4317] = (*ct)*(v[2715] * v[3174] + v[2723] * v[3175] + v[2730] * v[3176] + v[1338] * v[4298] + v[1339] * v[4300]
		+ v[1340] * v[4302]);
	v[4318] = (*ct)*(v[3174] * v[4236] + v[3175] * v[4260] + v[3176] * v[4284] + v[1424] * v[4298] + v[1425] * v[4300]
		+ v[1426] * v[4302]);
	v[4823] = v[3663] + v[4318];
	v[4319] = (*ct)*(v[3174] * v[4232] + v[3175] * v[4256] + v[3176] * v[4280] + v[1428] * v[4298] + v[1429] * v[4300]
		+ v[1430] * v[4302]);
	v[4828] = v[3664] + v[4319];
	v[4320] = (*ct)*(v[3174] * v[4228] + v[3175] * v[4252] + v[3176] * v[4276] + v[1432] * v[4298] + v[1433] * v[4300]
		+ v[1434] * v[4302]);
	v[4833] = v[3665] + v[4320];
	v[4375] = v[1183] * v[1409] + v[1184] * v[1410] + v[1185] * v[1411] + (*cn)*(v[3153] * v[3159] + v[3160] * v[3166]
		+ v[3167] * v[3173]) + (*ct)*(v[1409] * v[3174] + v[1410] * v[3175] + v[1411] * v[3176]) + (*epsn)*v[364];
	v[4376] = v[1183] * v[1413] + v[1184] * v[1414] + v[1185] * v[1415] + (*cn)*(v[3154] * v[3159] + v[3161] * v[3166]
		+ v[3168] * v[3173]) + (*ct)*(v[1413] * v[3174] + v[1414] * v[3175] + v[1415] * v[3176]) + (*epsn)*v[365];
	v[4377] = v[1183] * v[1417] + v[1184] * v[1418] + v[1185] * v[1419] + (*cn)*(v[3155] * v[3159] + v[3162] * v[3166]
		+ v[3169] * v[3173]) + (*ct)*(v[1417] * v[3174] + v[1418] * v[3175] + v[1419] * v[3176]) + (*epsn)*v[366];
	v[4384] = (*epst)*((v[1409] * v[1409]) + (v[1410] * v[1410]) + (v[1411] * v[1411])) + (*epsn)*v[3153] + (*cn)*
		(v[3173] * v[3259] + v[3166] * v[3265] + v[3159] * v[3271] + v[3153] * v[3463] + v[3160] * v[3469] + v[3167] * v[3475])
		+ v[1183] * v[3823] + v[1184] * v[3905] + v[1185] * v[3992] + (*ct)*(v[3174] * v[3823] + v[3175] * v[3905]
			+ v[3176] * v[3992] + v[1409] * v[4082] + v[1410] * v[4091] + v[1411] * v[4100]) + (v[4777] + v[4802])*v[699] + (v[4778]
				+ v[4803])*v[705];
	v[4385] = (*epsn)*v[3154] + (*cn)*(v[3173] * v[3260] + v[3166] * v[3266] + v[3159] * v[3272] + v[3153] * v[3464]
		+ v[3160] * v[3470] + v[3167] * v[3476]) + (*ct)*(v[3174] * v[3824] + v[3175] * v[3906] + v[3176] * v[3993]
			+ v[1409] * v[4083] + v[1410] * v[4092] + v[1411] * v[4101]) + v[4804] + v[4802] * v[700] + v[4803] * v[706];
	v[4386] = (*epsn)*v[3155] + (*cn)*(v[3173] * v[3261] + v[3166] * v[3267] + v[3159] * v[3273] + v[3153] * v[3465]
		+ v[3160] * v[3471] + v[3167] * v[3477]) + (*ct)*(v[3174] * v[3825] + v[3175] * v[3907] + v[3176] * v[3994]
			+ v[1409] * v[4084] + v[1410] * v[4093] + v[1411] * v[4102]) + v[4807] + v[4802] * v[701] + v[4803] * v[707];
	v[4387] = (*ct)*(v[3174] * v[3826] + v[3175] * v[3908] + v[3176] * v[3995] + v[1409] * v[4085] + v[1410] * v[4094]
		+ v[1411] * v[4103]) + v[4811];
	v[4388] = (*ct)*(v[3174] * v[3827] + v[3175] * v[3909] + v[3176] * v[3996] + v[1409] * v[4086] + v[1410] * v[4095]
		+ v[1411] * v[4104]) + v[4814];
	v[4389] = (*ct)*(v[3174] * v[3828] + v[3175] * v[3910] + v[3176] * v[3997] + v[1409] * v[4087] + v[1410] * v[4096]
		+ v[1411] * v[4105]) + v[4817];
	v[4390] = (*epsn)*v[3156] + (*cn)*(v[3173] * v[3262] + v[3166] * v[3268] + v[3159] * v[3274] + v[3153] * v[3466]
		+ v[3160] * v[3472] + v[3167] * v[3478]) + (*ct)*(v[3174] * v[3830] + v[3175] * v[3912] + v[3176] * v[3999]
			+ v[1409] * v[4088] + v[1410] * v[4097] + v[1411] * v[4106]) + v[4820] + v[4802] * v[702] + v[4803] * v[708];
	v[4391] = (*epsn)*v[3157] + (*cn)*(v[3173] * v[3263] + v[3166] * v[3269] + v[3159] * v[3275] + v[3153] * v[3467]
		+ v[3160] * v[3473] + v[3167] * v[3479]) + (*ct)*(v[3174] * v[3832] + v[3175] * v[3914] + v[3176] * v[4001]
			+ v[1409] * v[4089] + v[1410] * v[4098] + v[1411] * v[4107]) + v[4825] + v[4802] * v[703] + v[4803] * v[709];
	v[4392] = (*epsn)*v[3158] + (*cn)*(v[3173] * v[3264] + v[3166] * v[3270] + v[3159] * v[3276] + v[3153] * v[3468]
		+ v[3160] * v[3474] + v[3167] * v[3480]) + (*ct)*(v[3174] * v[3834] + v[3175] * v[3916] + v[3176] * v[4003]
			+ v[1409] * v[4090] + v[1410] * v[4099] + v[1411] * v[4108]) + v[4830] + v[4802] * v[704] + v[4803] * v[710];
	v[4393] = (*epsn)*v[3160] + (*cn)*(v[3173] * v[3277] + v[3166] * v[3283] + v[3159] * v[3289] + v[3154] * v[3463]
		+ v[3161] * v[3469] + v[3168] * v[3475]) + (*ct)*(v[3174] * v[3811] + v[3175] * v[3893] + v[3176] * v[3980]
			+ v[1413] * v[4082] + v[1414] * v[4091] + v[1415] * v[4100]) + v[4804] + v[4805] * v[699] + v[4806] * v[705];
	v[4394] = (*epst)*((v[1413] * v[1413]) + (v[1414] * v[1414]) + (v[1415] * v[1415])) + (*epsn)*v[3161] + (*cn)*
		(v[3173] * v[3278] + v[3166] * v[3284] + v[3159] * v[3290] + v[3154] * v[3464] + v[3161] * v[3470] + v[3168] * v[3476])
		+ v[1183] * v[3812] + v[1184] * v[3894] + v[1185] * v[3981] + (*ct)*(v[3174] * v[3812] + v[3175] * v[3894]
			+ v[3176] * v[3981] + v[1413] * v[4083] + v[1414] * v[4092] + v[1415] * v[4101]) + (v[4779] + v[4805])*v[700] + (v[4780]
				+ v[4806])*v[706];
	v[4395] = (*epsn)*v[3162] + (*cn)*(v[3173] * v[3279] + v[3166] * v[3285] + v[3159] * v[3291] + v[3154] * v[3465]
		+ v[3161] * v[3471] + v[3168] * v[3477]) + (*ct)*(v[3174] * v[3813] + v[3175] * v[3895] + v[3176] * v[3982]
			+ v[1413] * v[4084] + v[1414] * v[4093] + v[1415] * v[4102]) + v[4808] + v[4805] * v[701] + v[4806] * v[707];
	v[4396] = (*ct)*(v[3174] * v[3814] + v[3175] * v[3896] + v[3176] * v[3983] + v[1413] * v[4085] + v[1414] * v[4094]
		+ v[1415] * v[4103]) + v[4812];
	v[4397] = (*ct)*(v[3174] * v[3815] + v[3175] * v[3897] + v[3176] * v[3984] + v[1413] * v[4086] + v[1414] * v[4095]
		+ v[1415] * v[4104]) + v[4815];
	v[4398] = (*ct)*(v[3174] * v[3816] + v[3175] * v[3898] + v[3176] * v[3985] + v[1413] * v[4087] + v[1414] * v[4096]
		+ v[1415] * v[4105]) + v[4818];
	v[4399] = (*epsn)*v[3163] + (*cn)*(v[3173] * v[3280] + v[3166] * v[3286] + v[3159] * v[3292] + v[3154] * v[3466]
		+ v[3161] * v[3472] + v[3168] * v[3478]) + (*ct)*(v[3174] * v[3818] + v[3175] * v[3900] + v[3176] * v[3987]
			+ v[1413] * v[4088] + v[1414] * v[4097] + v[1415] * v[4106]) + v[4821] + v[4805] * v[702] + v[4806] * v[708];
	v[4400] = (*epsn)*v[3164] + (*cn)*(v[3173] * v[3281] + v[3166] * v[3287] + v[3159] * v[3293] + v[3154] * v[3467]
		+ v[3161] * v[3473] + v[3168] * v[3479]) + (*ct)*(v[3174] * v[3820] + v[3175] * v[3902] + v[3176] * v[3989]
			+ v[1413] * v[4089] + v[1414] * v[4098] + v[1415] * v[4107]) + v[4826] + v[4805] * v[703] + v[4806] * v[709];
	v[4401] = (*epsn)*v[3165] + (*cn)*(v[3173] * v[3282] + v[3166] * v[3288] + v[3159] * v[3294] + v[3154] * v[3468]
		+ v[3161] * v[3474] + v[3168] * v[3480]) + (*ct)*(v[3174] * v[3822] + v[3175] * v[3904] + v[3176] * v[3991]
			+ v[1413] * v[4090] + v[1414] * v[4099] + v[1415] * v[4108]) + v[4831] + v[4805] * v[704] + v[4806] * v[710];
	v[4402] = (*epsn)*v[3167] + (*cn)*(v[3173] * v[3295] + v[3166] * v[3301] + v[3159] * v[3307] + v[3155] * v[3463]
		+ v[3162] * v[3469] + v[3169] * v[3475]) + (*ct)*(v[3174] * v[3799] + v[3175] * v[3879] + v[3176] * v[3964]
			+ v[1417] * v[4082] + v[1418] * v[4091] + v[1419] * v[4100]) + v[4807] + v[4809] * v[699] + v[4810] * v[705];
	v[4403] = (*epsn)*v[3168] + (*cn)*(v[3173] * v[3296] + v[3166] * v[3302] + v[3159] * v[3308] + v[3155] * v[3464]
		+ v[3162] * v[3470] + v[3169] * v[3476]) + (*ct)*(v[3174] * v[3800] + v[3175] * v[3880] + v[3176] * v[3967]
			+ v[1417] * v[4083] + v[1418] * v[4092] + v[1419] * v[4101]) + v[4808] + v[4809] * v[700] + v[4810] * v[706];
	v[4404] = (*epst)*((v[1417] * v[1417]) + (v[1418] * v[1418]) + (v[1419] * v[1419])) + (*epsn)*v[3169] + (*cn)*
		(v[3173] * v[3297] + v[3166] * v[3303] + v[3159] * v[3309] + v[3155] * v[3465] + v[3162] * v[3471] + v[3169] * v[3477])
		+ v[1183] * v[3801] + v[1184] * v[3881] + v[1185] * v[3968] + (*ct)*(v[3174] * v[3801] + v[3175] * v[3881]
			+ v[3176] * v[3968] + v[1417] * v[4084] + v[1418] * v[4093] + v[1419] * v[4102]) + (v[4781] + v[4809])*v[701] + (v[4782]
				+ v[4810])*v[707];
	v[4405] = (*ct)*(v[3174] * v[3802] + v[3175] * v[3882] + v[3176] * v[3969] + v[1417] * v[4085] + v[1418] * v[4094]
		+ v[1419] * v[4103]) + v[4813];
	v[4406] = (*ct)*(v[3174] * v[3803] + v[3175] * v[3883] + v[3176] * v[3970] + v[1417] * v[4086] + v[1418] * v[4095]
		+ v[1419] * v[4104]) + v[4816];
	v[4407] = (*ct)*(v[3174] * v[3804] + v[3175] * v[3884] + v[3176] * v[3971] + v[1417] * v[4087] + v[1418] * v[4096]
		+ v[1419] * v[4105]) + v[4819];
	v[4408] = (*epsn)*v[3170] + (*cn)*(v[3173] * v[3298] + v[3166] * v[3304] + v[3159] * v[3310] + v[3155] * v[3466]
		+ v[3162] * v[3472] + v[3169] * v[3478]) + (*ct)*(v[3174] * v[3806] + v[3175] * v[3886] + v[3176] * v[3973]
			+ v[1417] * v[4088] + v[1418] * v[4097] + v[1419] * v[4106]) + v[4824] + v[4809] * v[702] + v[4810] * v[708];
	v[4409] = (*epsn)*v[3171] + (*cn)*(v[3173] * v[3299] + v[3166] * v[3305] + v[3159] * v[3311] + v[3155] * v[3467]
		+ v[3162] * v[3473] + v[3169] * v[3479]) + (*ct)*(v[3174] * v[3808] + v[3175] * v[3888] + v[3176] * v[3975]
			+ v[1417] * v[4089] + v[1418] * v[4098] + v[1419] * v[4107]) + v[4829] + v[4809] * v[703] + v[4810] * v[709];
	v[4410] = (*epsn)*v[3172] + (*cn)*(v[3173] * v[3300] + v[3166] * v[3306] + v[3159] * v[3312] + v[3155] * v[3468]
		+ v[3162] * v[3474] + v[3169] * v[3480]) + (*ct)*(v[3174] * v[3810] + v[3175] * v[3890] + v[3176] * v[3977]
			+ v[1417] * v[4090] + v[1418] * v[4099] + v[1419] * v[4108]) + v[4834] + v[4809] * v[704] + v[4810] * v[710];
	v[4411] = (*ct)*(v[1332] * v[4082] + v[1333] * v[4091] + v[1334] * v[4100]) + v[4811] + v[4309] * v[699]
		+ v[4315] * v[705];
	v[4412] = (*ct)*(v[1332] * v[4083] + v[1333] * v[4092] + v[1334] * v[4101]) + v[4812] + v[4309] * v[700]
		+ v[4315] * v[706];
	v[4413] = (*ct)*(v[1332] * v[4084] + v[1333] * v[4093] + v[1334] * v[4102]) + v[4813] + v[4309] * v[701]
		+ v[4315] * v[707];
	v[4420] = (*ct)*(v[1335] * v[4082] + v[1336] * v[4091] + v[1337] * v[4100]) + v[4814] + v[4310] * v[699]
		+ v[4316] * v[705];
	v[4421] = (*ct)*(v[1335] * v[4083] + v[1336] * v[4092] + v[1337] * v[4101]) + v[4815] + v[4310] * v[700]
		+ v[4316] * v[706];
	v[4422] = (*ct)*(v[1335] * v[4084] + v[1336] * v[4093] + v[1337] * v[4102]) + v[4816] + v[4310] * v[701]
		+ v[4316] * v[707];
	v[4429] = (*ct)*(v[1338] * v[4082] + v[1339] * v[4091] + v[1340] * v[4100]) + v[4817] + v[4311] * v[699]
		+ v[4317] * v[705];
	v[4430] = (*ct)*(v[1338] * v[4083] + v[1339] * v[4092] + v[1340] * v[4101]) + v[4818] + v[4311] * v[700]
		+ v[4317] * v[706];
	v[4431] = (*ct)*(v[1338] * v[4084] + v[1339] * v[4093] + v[1340] * v[4102]) + v[4819] + v[4311] * v[701]
		+ v[4317] * v[707];
	v[4438] = (*cn)*(v[3173] * v[3313] + v[3166] * v[3322] + v[3159] * v[3331] + v[3156] * v[3463] + v[3163] * v[3469]
		+ v[3170] * v[3475]) + (*ct)*(v[3174] * v[3786] + v[3175] * v[3864] + v[3176] * v[3948] + v[1424] * v[4082]
			+ v[1425] * v[4091] + v[1426] * v[4100]) + v[4820] + v[4822] * v[699] + v[4823] * v[705] + (*epsn)*(v[1059]
				+ v[1123] * v[699] + v[1124] * v[705]);
	v[4439] = (*cn)*(v[3173] * v[3314] + v[3166] * v[3323] + v[3159] * v[3332] + v[3156] * v[3464] + v[3163] * v[3470]
		+ v[3170] * v[3476]) + (*ct)*(v[3174] * v[3787] + v[3175] * v[3865] + v[3176] * v[3950] + v[1424] * v[4083]
			+ v[1425] * v[4092] + v[1426] * v[4101]) + v[4821] + v[4822] * v[700] + v[4823] * v[706] + (*epsn)*(v[1062]
				+ v[1123] * v[700] + v[1124] * v[706]);
	v[4440] = (*cn)*(v[3173] * v[3315] + v[3166] * v[3324] + v[3159] * v[3333] + v[3156] * v[3465] + v[3163] * v[3471]
		+ v[3170] * v[3477]) + (*ct)*(v[3174] * v[3788] + v[3175] * v[3866] + v[3176] * v[3951] + v[1424] * v[4084]
			+ v[1425] * v[4093] + v[1426] * v[4102]) + v[4824] + v[4822] * v[701] + v[4823] * v[707] + (*epsn)*(v[1065]
				+ v[1123] * v[701] + v[1124] * v[707]);
	v[4447] = (*cn)*(v[3173] * v[3340] + v[3166] * v[3348] + v[3159] * v[3356] + v[3157] * v[3463] + v[3164] * v[3469]
		+ v[3171] * v[3475]) + (*ct)*(v[3174] * v[3773] + v[3175] * v[3850] + v[3176] * v[3933] + v[1428] * v[4082]
			+ v[1429] * v[4091] + v[1430] * v[4100]) + v[4825] + v[4827] * v[699] + v[4828] * v[705] + (*epsn)*(v[1060]
				+ v[1131] * v[699] + v[1132] * v[705]);
	v[4448] = (*cn)*(v[3173] * v[3341] + v[3166] * v[3349] + v[3159] * v[3357] + v[3157] * v[3464] + v[3164] * v[3470]
		+ v[3171] * v[3476]) + (*ct)*(v[3174] * v[3774] + v[3175] * v[3851] + v[3176] * v[3935] + v[1428] * v[4083]
			+ v[1429] * v[4092] + v[1430] * v[4101]) + v[4826] + v[4827] * v[700] + v[4828] * v[706] + (*epsn)*(v[1063]
				+ v[1131] * v[700] + v[1132] * v[706]);
	v[4449] = (*cn)*(v[3173] * v[3342] + v[3166] * v[3350] + v[3159] * v[3358] + v[3157] * v[3465] + v[3164] * v[3471]
		+ v[3171] * v[3477]) + (*ct)*(v[3174] * v[3775] + v[3175] * v[3852] + v[3176] * v[3936] + v[1428] * v[4084]
			+ v[1429] * v[4093] + v[1430] * v[4102]) + v[4829] + v[4827] * v[701] + v[4828] * v[707] + (*epsn)*(v[1066]
				+ v[1131] * v[701] + v[1132] * v[707]);
	v[4456] = (*cn)*(v[3173] * v[3364] + v[3166] * v[3371] + v[3159] * v[3378] + v[3158] * v[3463] + v[3165] * v[3469]
		+ v[3172] * v[3475]) + (*ct)*(v[3174] * v[3760] + v[3175] * v[3836] + v[3176] * v[3918] + v[1432] * v[4082]
			+ v[1433] * v[4091] + v[1434] * v[4100]) + v[4830] + v[4832] * v[699] + v[4833] * v[705] + (*epsn)*(v[1061]
				+ v[1140] * v[699] + v[1141] * v[705]);
	v[4457] = (*cn)*(v[3173] * v[3365] + v[3166] * v[3372] + v[3159] * v[3379] + v[3158] * v[3464] + v[3165] * v[3470]
		+ v[3172] * v[3476]) + (*ct)*(v[3174] * v[3761] + v[3175] * v[3837] + v[3176] * v[3920] + v[1432] * v[4083]
			+ v[1433] * v[4092] + v[1434] * v[4101]) + v[4831] + v[4832] * v[700] + v[4833] * v[706] + (*epsn)*(v[1064]
				+ v[1140] * v[700] + v[1141] * v[706]);
	v[4458] = (*cn)*(v[3173] * v[3366] + v[3166] * v[3373] + v[3159] * v[3380] + v[3158] * v[3465] + v[3165] * v[3471]
		+ v[3172] * v[3477]) + (*ct)*(v[3174] * v[3762] + v[3175] * v[3838] + v[3176] * v[3921] + v[1432] * v[4084]
			+ v[1433] * v[4093] + v[1434] * v[4102]) + v[4834] + v[4832] * v[701] + v[4833] * v[707] + (*epsn)*(v[1067]
				+ v[1140] * v[701] + v[1141] * v[707]);
	Rc[0] = v[4375];
	Rc[1] = v[4376];
	Rc[2] = v[4377];
	Rc[3] = v[1183] * v[1332] + v[1184] * v[1333] + v[1185] * v[1334] + (*ct)*(v[1332] * v[3174] + v[1333] * v[3175]
		+ v[1334] * v[3176]);
	Rc[4] = v[1183] * v[1335] + v[1184] * v[1336] + v[1185] * v[1337] + (*ct)*(v[1335] * v[3174] + v[1336] * v[3175]
		+ v[1337] * v[3176]);
	Rc[5] = v[1183] * v[1338] + v[1184] * v[1339] + v[1185] * v[1340] + (*ct)*(v[1338] * v[3174] + v[1339] * v[3175]
		+ v[1340] * v[3176]);
	Rc[6] = -v[4375];
	Rc[7] = -v[4376];
	Rc[8] = -v[4377];
	Rc[9] = v[1183] * v[1424] + v[1184] * v[1425] + v[1185] * v[1426] + (*cn)*(v[3156] * v[3159] + v[3163] * v[3166]
		+ v[3170] * v[3173]) + (*ct)*(v[1424] * v[3174] + v[1425] * v[3175] + v[1426] * v[3176]) + (*epsn)*(v[1059] * v[364]
			+ v[1062] * v[365] + v[1065] * v[366]);
	Rc[10] = v[1183] * v[1428] + v[1184] * v[1429] + v[1185] * v[1430] + (*cn)*(v[3157] * v[3159] + v[3164] * v[3166]
		+ v[3171] * v[3173]) + (*ct)*(v[1428] * v[3174] + v[1429] * v[3175] + v[1430] * v[3176]) + (*epsn)*(v[1060] * v[364]
			+ v[1063] * v[365] + v[1066] * v[366]);
	Rc[11] = v[1183] * v[1432] + v[1184] * v[1433] + v[1185] * v[1434] + (*cn)*(v[3158] * v[3159] + v[3165] * v[3166]
		+ v[3172] * v[3173]) + (*ct)*(v[1432] * v[3174] + v[1433] * v[3175] + v[1434] * v[3176]) + (*epsn)*(v[1061] * v[364]
			+ v[1064] * v[365] + v[1067] * v[366]);
	Kc[0][0] = v[4384];
	Kc[0][1] = v[4385];
	Kc[0][2] = v[4386];
	Kc[0][3] = v[4387];
	Kc[0][4] = v[4388];
	Kc[0][5] = v[4389];
	Kc[0][6] = -v[4384];
	Kc[0][7] = -v[4385];
	Kc[0][8] = -v[4386];
	Kc[0][9] = v[4390];
	Kc[0][10] = v[4391];
	Kc[0][11] = v[4392];
	Kc[1][0] = v[4393];
	Kc[1][1] = v[4394];
	Kc[1][2] = v[4395];
	Kc[1][3] = v[4396];
	Kc[1][4] = v[4397];
	Kc[1][5] = v[4398];
	Kc[1][6] = -v[4393];
	Kc[1][7] = -v[4394];
	Kc[1][8] = -v[4395];
	Kc[1][9] = v[4399];
	Kc[1][10] = v[4400];
	Kc[1][11] = v[4401];
	Kc[2][0] = v[4402];
	Kc[2][1] = v[4403];
	Kc[2][2] = v[4404];
	Kc[2][3] = v[4405];
	Kc[2][4] = v[4406];
	Kc[2][5] = v[4407];
	Kc[2][6] = -v[4402];
	Kc[2][7] = -v[4403];
	Kc[2][8] = -v[4404];
	Kc[2][9] = v[4408];
	Kc[2][10] = v[4409];
	Kc[2][11] = v[4410];
	Kc[3][0] = v[4411];
	Kc[3][1] = v[4412];
	Kc[3][2] = v[4413];
	Kc[3][3] = (*epst)*((v[1332] * v[1332]) + (v[1333] * v[1333]) + (v[1334] * v[1334])) + v[1183] * v[3715]
		+ v[1184] * v[3721] + v[1185] * v[3727] + (*ct)*(v[3174] * v[3715] + v[3175] * v[3721] + v[3176] * v[3727]
			+ v[1332] * v[4085] + v[1333] * v[4094] + v[1334] * v[4103]);
	Kc[3][4] = (*ct)*(v[1332] * v[4086] + v[1333] * v[4095] + v[1334] * v[4104]) + v[4835];
	Kc[3][5] = (*ct)*(v[1332] * v[4087] + v[1333] * v[4096] + v[1334] * v[4105]) + v[4836];
	Kc[3][6] = -v[4411];
	Kc[3][7] = -v[4412];
	Kc[3][8] = -v[4413];
	Kc[3][9] = (*ct)*(v[3174] * v[3718] + v[3175] * v[3724] + v[3176] * v[3730] + v[1332] * v[4088] + v[1333] * v[4097]
		+ v[1334] * v[4106]) + v[4837] + v[4309] * v[702] + v[4315] * v[708];
	Kc[3][10] = (*ct)*(v[3174] * v[3719] + v[3175] * v[3725] + v[3176] * v[3731] + v[1332] * v[4089] + v[1333] * v[4098]
		+ v[1334] * v[4107]) + v[4838] + v[4309] * v[703] + v[4315] * v[709];
	Kc[3][11] = (*ct)*(v[3174] * v[3720] + v[3175] * v[3726] + v[3176] * v[3732] + v[1332] * v[4090] + v[1333] * v[4099]
		+ v[1334] * v[4108]) + v[4839] + v[4309] * v[704] + v[4315] * v[710];
	Kc[4][0] = v[4420];
	Kc[4][1] = v[4421];
	Kc[4][2] = v[4422];
	Kc[4][3] = (*ct)*(v[1335] * v[4085] + v[1336] * v[4094] + v[1337] * v[4103]) + v[4835];
	Kc[4][4] = (*epst)*((v[1335] * v[1335]) + (v[1336] * v[1336]) + (v[1337] * v[1337])) + v[1183] * v[3733]
		+ v[1184] * v[3738] + v[1185] * v[3743] + (*ct)*(v[3174] * v[3733] + v[3175] * v[3738] + v[3176] * v[3743]
			+ v[1335] * v[4086] + v[1336] * v[4095] + v[1337] * v[4104]);
	Kc[4][5] = (*ct)*(v[1335] * v[4087] + v[1336] * v[4096] + v[1337] * v[4105]) + v[4840];
	Kc[4][6] = -v[4420];
	Kc[4][7] = -v[4421];
	Kc[4][8] = -v[4422];
	Kc[4][9] = (*ct)*(v[3174] * v[3735] + v[3175] * v[3740] + v[3176] * v[3745] + v[1335] * v[4088] + v[1336] * v[4097]
		+ v[1337] * v[4106]) + v[4841] + v[4310] * v[702] + v[4316] * v[708];
	Kc[4][10] = (*ct)*(v[3174] * v[3736] + v[3175] * v[3741] + v[3176] * v[3746] + v[1335] * v[4089] + v[1336] * v[4098]
		+ v[1337] * v[4107]) + v[4842] + v[4310] * v[703] + v[4316] * v[709];
	Kc[4][11] = (*ct)*(v[3174] * v[3737] + v[3175] * v[3742] + v[3176] * v[3747] + v[1335] * v[4090] + v[1336] * v[4099]
		+ v[1337] * v[4108]) + v[4843] + v[4310] * v[704] + v[4316] * v[710];
	Kc[5][0] = v[4429];
	Kc[5][1] = v[4430];
	Kc[5][2] = v[4431];
	Kc[5][3] = (*ct)*(v[1338] * v[4085] + v[1339] * v[4094] + v[1340] * v[4103]) + v[4836];
	Kc[5][4] = (*ct)*(v[1338] * v[4086] + v[1339] * v[4095] + v[1340] * v[4104]) + v[4840];
	Kc[5][5] = (*epst)*((v[1338] * v[1338]) + (v[1339] * v[1339]) + (v[1340] * v[1340])) + v[1183] * v[3748]
		+ v[1184] * v[3752] + v[1185] * v[3756] + (*ct)*(v[3174] * v[3748] + v[3175] * v[3752] + v[3176] * v[3756]
			+ v[1338] * v[4087] + v[1339] * v[4096] + v[1340] * v[4105]);
	Kc[5][6] = -v[4429];
	Kc[5][7] = -v[4430];
	Kc[5][8] = -v[4431];
	Kc[5][9] = (*ct)*(v[3174] * v[3749] + v[3175] * v[3753] + v[3176] * v[3757] + v[1338] * v[4088] + v[1339] * v[4097]
		+ v[1340] * v[4106]) + v[4844] + v[4311] * v[702] + v[4317] * v[708];
	Kc[5][10] = (*ct)*(v[3174] * v[3750] + v[3175] * v[3754] + v[3176] * v[3758] + v[1338] * v[4089] + v[1339] * v[4098]
		+ v[1340] * v[4107]) + v[4845] + v[4311] * v[703] + v[4317] * v[709];
	Kc[5][11] = (*ct)*(v[3174] * v[3751] + v[3175] * v[3755] + v[3176] * v[3759] + v[1338] * v[4090] + v[1339] * v[4099]
		+ v[1340] * v[4108]) + v[4846] + v[4311] * v[704] + v[4317] * v[710];
	Kc[6][0] = -v[4384];
	Kc[6][1] = -v[4385];
	Kc[6][2] = -v[4386];
	Kc[6][3] = -v[4387];
	Kc[6][4] = -v[4388];
	Kc[6][5] = -v[4389];
	Kc[6][6] = v[4384];
	Kc[6][7] = v[4385];
	Kc[6][8] = v[4386];
	Kc[6][9] = -v[4390];
	Kc[6][10] = -v[4391];
	Kc[6][11] = -v[4392];
	Kc[7][0] = -v[4393];
	Kc[7][1] = -v[4394];
	Kc[7][2] = -v[4395];
	Kc[7][3] = -v[4396];
	Kc[7][4] = -v[4397];
	Kc[7][5] = -v[4398];
	Kc[7][6] = v[4393];
	Kc[7][7] = v[4394];
	Kc[7][8] = v[4395];
	Kc[7][9] = -v[4399];
	Kc[7][10] = -v[4400];
	Kc[7][11] = -v[4401];
	Kc[8][0] = -v[4402];
	Kc[8][1] = -v[4403];
	Kc[8][2] = -v[4404];
	Kc[8][3] = -v[4405];
	Kc[8][4] = -v[4406];
	Kc[8][5] = -v[4407];
	Kc[8][6] = v[4402];
	Kc[8][7] = v[4403];
	Kc[8][8] = v[4404];
	Kc[8][9] = -v[4408];
	Kc[8][10] = -v[4409];
	Kc[8][11] = -v[4410];
	Kc[9][0] = v[4438];
	Kc[9][1] = v[4439];
	Kc[9][2] = v[4440];
	Kc[9][3] = (*ct)*(v[3174] * v[3789] + v[3175] * v[3867] + v[3176] * v[3952] + v[1424] * v[4085] + v[1425] * v[4094]
		+ v[1426] * v[4103]) + v[4837];
	Kc[9][4] = (*ct)*(v[3174] * v[3790] + v[3175] * v[3868] + v[3176] * v[3953] + v[1424] * v[4086] + v[1425] * v[4095]
		+ v[1426] * v[4104]) + v[4841];
	Kc[9][5] = (*ct)*(v[3174] * v[3791] + v[3175] * v[3869] + v[3176] * v[3954] + v[1424] * v[4087] + v[1425] * v[4096]
		+ v[1426] * v[4105]) + v[4844];
	Kc[9][6] = -v[4438];
	Kc[9][7] = -v[4439];
	Kc[9][8] = -v[4440];
	Kc[9][9] = (*epst)*((v[1424] * v[1424]) + (v[1425] * v[1425]) + (v[1426] * v[1426])) + (*cn)*(v[3173] * v[3317]
		+ v[3166] * v[3326] + v[3159] * v[3335] + v[3156] * v[3466] + v[3163] * v[3472] + v[3170] * v[3478]) + v[1183] * v[3794]
		+ v[1184] * v[3872] + v[1185] * v[3957] + (*ct)*(v[3174] * v[3794] + v[3175] * v[3872] + v[3176] * v[3957]
			+ v[1424] * v[4088] + v[1425] * v[4097] + v[1426] * v[4106]) + (v[4784] + v[4822])*v[702] + (v[4783] + v[4823])*v[708]
		+ (*epsn)*((v[1059] * v[1059]) + (v[1062] * v[1062]) + (v[1065] * v[1065]) + v[3334] * v[364] + v[3325] * v[365]
			+ v[3316] * v[366] + v[1123] * v[702] + v[1124] * v[708]);
	Kc[9][10] = (*cn)*(v[3173] * v[3319] + v[3166] * v[3328] + v[3159] * v[3337] + v[3156] * v[3467] + v[3163] * v[3473]
		+ v[3170] * v[3479]) + (*ct)*(v[3174] * v[3796] + v[3175] * v[3874] + v[3176] * v[3959] + v[1424] * v[4089]
			+ v[1425] * v[4098] + v[1426] * v[4107]) + v[4847] + v[4822] * v[703] + v[4823] * v[709] + (*epsn)*(v[1135]
				+ v[1123] * v[703] + v[1124] * v[709]);
	Kc[9][11] = (*cn)*(v[3173] * v[3321] + v[3166] * v[3330] + v[3159] * v[3339] + v[3156] * v[3468] + v[3163] * v[3474]
		+ v[3170] * v[3480]) + (*ct)*(v[3174] * v[3798] + v[3175] * v[3876] + v[3176] * v[3961] + v[1424] * v[4090]
			+ v[1425] * v[4099] + v[1426] * v[4108]) + v[4848] + v[4822] * v[704] + v[4823] * v[710] + (*epsn)*(v[1144]
				+ v[1123] * v[704] + v[1124] * v[710]);
	Kc[10][0] = v[4447];
	Kc[10][1] = v[4448];
	Kc[10][2] = v[4449];
	Kc[10][3] = (*ct)*(v[3174] * v[3776] + v[3175] * v[3853] + v[3176] * v[3937] + v[1428] * v[4085] + v[1429] * v[4094]
		+ v[1430] * v[4103]) + v[4838];
	Kc[10][4] = (*ct)*(v[3174] * v[3777] + v[3175] * v[3854] + v[3176] * v[3938] + v[1428] * v[4086] + v[1429] * v[4095]
		+ v[1430] * v[4104]) + v[4842];
	Kc[10][5] = (*ct)*(v[3174] * v[3778] + v[3175] * v[3855] + v[3176] * v[3939] + v[1428] * v[4087] + v[1429] * v[4096]
		+ v[1430] * v[4105]) + v[4845];
	Kc[10][6] = -v[4447];
	Kc[10][7] = -v[4448];
	Kc[10][8] = -v[4449];
	Kc[10][9] = (*cn)*(v[3173] * v[3343] + v[3166] * v[3351] + v[3159] * v[3359] + v[3157] * v[3466] + v[3164] * v[3472]
		+ v[3171] * v[3478]) + (*ct)*(v[3174] * v[3780] + v[3175] * v[3857] + v[3176] * v[3941] + v[1428] * v[4088]
			+ v[1429] * v[4097] + v[1430] * v[4106]) + v[4847] + v[4827] * v[702] + v[4828] * v[708] + (*epsn)*(v[1135]
				+ v[1131] * v[702] + v[1132] * v[708]);
	Kc[10][10] = (*epst)*((v[1428] * v[1428]) + (v[1429] * v[1429]) + (v[1430] * v[1430])) + (*cn)*(v[3173] * v[3345]
		+ v[3166] * v[3353] + v[3159] * v[3361] + v[3157] * v[3467] + v[3164] * v[3473] + v[3171] * v[3479]) + v[1183] * v[3783]
		+ v[1184] * v[3860] + v[1185] * v[3944] + (*ct)*(v[3174] * v[3783] + v[3175] * v[3860] + v[3176] * v[3944]
			+ v[1428] * v[4089] + v[1429] * v[4098] + v[1430] * v[4107]) + (v[4827] + v[4849])*v[703] + (v[4828] + v[4850])*v[709]
		+ (*epsn)*((v[1060] * v[1060]) + (v[1063] * v[1063]) + (v[1066] * v[1066]) + v[3360] * v[364] + v[3352] * v[365]
			+ v[3344] * v[366] + v[1131] * v[703] + v[1132] * v[709]);
	Kc[10][11] = (*cn)*(v[3173] * v[3347] + v[3166] * v[3355] + v[3159] * v[3363] + v[3157] * v[3468] + v[3164] * v[3474]
		+ v[3171] * v[3480]) + (*ct)*(v[3174] * v[3785] + v[3175] * v[3862] + v[3176] * v[3946] + v[1428] * v[4090]
			+ v[1429] * v[4099] + v[1430] * v[4108]) + v[4851] + v[4827] * v[704] + v[4828] * v[710] + (*epsn)*(v[1146]
				+ v[1131] * v[704] + v[1132] * v[710]);
	Kc[11][0] = v[4456];
	Kc[11][1] = v[4457];
	Kc[11][2] = v[4458];
	Kc[11][3] = (*ct)*(v[3174] * v[3763] + v[3175] * v[3839] + v[3176] * v[3922] + v[1432] * v[4085] + v[1433] * v[4094]
		+ v[1434] * v[4103]) + v[4839];
	Kc[11][4] = (*ct)*(v[3174] * v[3764] + v[3175] * v[3840] + v[3176] * v[3923] + v[1432] * v[4086] + v[1433] * v[4095]
		+ v[1434] * v[4104]) + v[4843];
	Kc[11][5] = (*ct)*(v[3174] * v[3765] + v[3175] * v[3841] + v[3176] * v[3924] + v[1432] * v[4087] + v[1433] * v[4096]
		+ v[1434] * v[4105]) + v[4846];
	Kc[11][6] = -v[4456];
	Kc[11][7] = -v[4457];
	Kc[11][8] = -v[4458];
	Kc[11][9] = (*cn)*(v[3173] * v[3367] + v[3166] * v[3374] + v[3159] * v[3381] + v[3158] * v[3466] + v[3165] * v[3472]
		+ v[3172] * v[3478]) + (*ct)*(v[3174] * v[3767] + v[3175] * v[3843] + v[3176] * v[3926] + v[1432] * v[4088]
			+ v[1433] * v[4097] + v[1434] * v[4106]) + v[4848] + v[4832] * v[702] + v[4833] * v[708] + (*epsn)*(v[1144]
				+ v[1140] * v[702] + v[1141] * v[708]);
	Kc[11][10] = (*cn)*(v[3173] * v[3368] + v[3166] * v[3375] + v[3159] * v[3382] + v[3158] * v[3467] + v[3165] * v[3473]
		+ v[3172] * v[3479]) + (*ct)*(v[3174] * v[3769] + v[3175] * v[3845] + v[3176] * v[3928] + v[1432] * v[4089]
			+ v[1433] * v[4098] + v[1434] * v[4107]) + v[4851] + v[4832] * v[703] + v[4833] * v[709] + (*epsn)*(v[1146]
				+ v[1140] * v[703] + v[1141] * v[709]);
	Kc[11][11] = (*epst)*((v[1432] * v[1432]) + (v[1433] * v[1433]) + (v[1434] * v[1434])) + (*cn)*(v[3173] * v[3370]
		+ v[3166] * v[3377] + v[3159] * v[3384] + v[3158] * v[3468] + v[3165] * v[3474] + v[3172] * v[3480]) + v[1183] * v[3772]
		+ v[1184] * v[3848] + v[1185] * v[3931] + (*ct)*(v[3174] * v[3772] + v[3175] * v[3848] + v[3176] * v[3931]
			+ v[1432] * v[4090] + v[1433] * v[4099] + v[1434] * v[4108]) + (v[1183] * v[4226] + v[1184] * v[4250] + v[1185] * v[4274]
				+ v[4832])*v[704] + (v[1183] * v[4228] + v[1184] * v[4252] + v[1185] * v[4276] + v[4833])*v[710] + (*epsn)*(
				(v[1061] * v[1061]) + (v[1064] * v[1064]) + (v[1067] * v[1067]) + v[3383] * v[364] + v[3376] * v[365] + v[3369] * v[366]
					+ v[1140] * v[704] + v[1141] * v[710]);
}

//Calcula contribuições de contato entre esfera e superfície
void RigidOscillatorySurface_1::ContactSphereSurfaceSliding(double* Rc, double** Kc, double zetap, double thetap, double zetai, double thetai, double* gti, int node, double* epsn, double* epst, double* cn, double* ct, double* mu, double* radius)
{
	double *xPi = xP_i->getMatrix();	//ponteiro para o vetor xPi
	double* e1i = e1_i->getMatrix();
	double* e2i = e2_i->getMatrix();
	double* A1 = &A_1;
	double* A2 = &A_2;
	double* A12 = &A_12;
	double* lambda1 = &lambda_1;
	double* lambda2 = &lambda_2;
	double* phi1 = &phi_1;
	double* phi2 = &phi_2;

	double ci[2];
	ci[0] = zetai;
	ci[1] = thetai;
	double cp[2];
	cp[0] = zetap;
	cp[1] = thetap;
	double xSi[3];
	double d[12];
	for (int i = 0; i < 3; i++)
	{
		/*(*1 - 6: sphere*)
		(*7 - 12: surface*)*/
		xSi[i] = db.nodes[node - 1]->copy_coordinates[i];
		d[i] = db.nodes[node - 1]->displacements[i];
		d[i + 3] = db.nodes[node - 1]->displacements[i + 3];
		d[i + 6] = db.nodes[pilot_node - 1]->displacements[i];
		d[i + 9] = db.nodes[pilot_node - 1]->displacements[i + 3];
	}

	double* a4;
	double* a5;
	double* a6;
	double value = 0.0;
	if (typeid(*db.solution[db.current_solution_number - 1]) == typeid(Dynamic))
	{
		Dynamic* ptr_sol = static_cast<Dynamic*>(db.solution[db.current_solution_number - 1]);
		a4 = &ptr_sol->a4;
		a5 = &ptr_sol->a5;
		a6 = &ptr_sol->a6;
	}
	else
	{
		a4 = &value;
		a5 = &value;
		a6 = &value;
	}
	double omegaiS[3];
	double domegaiS[3];
	double duiS[3];
	double dduiS[3];
	double omegaiP[3];
	double domegaiP[3];
	double duiP[3];
	double dduiP[3];
	for (int i = 0; i < 3; i++)
	{
		duiS[i] = db.nodes[node - 1]->copy_vel[i];
		omegaiS[i] = db.nodes[node - 1]->copy_vel[i + 3];
		dduiS[i] = db.nodes[node - 1]->copy_accel[i];
		domegaiS[i] = db.nodes[node - 1]->copy_accel[i + 3];

		duiP[i] = db.nodes[pilot_node - 1]->copy_vel[i];
		omegaiP[i] = db.nodes[pilot_node - 1]->copy_vel[i + 3];
		dduiP[i] = db.nodes[pilot_node - 1]->copy_accel[i];
		domegaiP[i] = db.nodes[pilot_node - 1]->copy_accel[i + 3];
	}
	double v[10000];
	//AceGen
	v[4358] = (*a4)*d[2] + (*a6)*dduiS[2] + (*a5)*duiS[2];
	v[4357] = (*a4)*d[1] + (*a6)*dduiS[1] + (*a5)*duiS[1];
	v[4356] = (*a4)*d[0] + (*a6)*dduiS[0] + (*a5)*duiS[0];
	v[4351] = (*a4)*d[8] + (*a6)*dduiP[2] + (*a5)*duiP[2];
	v[4361] = -v[4351] + v[4358];
	v[4350] = (*a4)*d[7] + (*a6)*dduiP[1] + (*a5)*duiP[1];
	v[4360] = -v[4350] + v[4357];
	v[4349] = (*a4)*d[6] + (*a6)*dduiP[0] + (*a5)*duiP[0];
	v[4359] = -v[4349] + v[4356];
	v[4297] = (*epsn)*(*mu);
	v[4076] = -(e1i[1] * e2i[0]) + e1i[0] * e2i[1];
	v[4075] = e1i[2] * e2i[0] - e1i[0] * e2i[2];
	v[4074] = -(e1i[2] * e2i[1]) + e1i[1] * e2i[2];
	v[4061] = d[2] + xSi[2];
	v[4060] = d[1] + xSi[1];
	v[4059] = d[0] + xSi[0];
	v[4058] = d[8] + xPi[2];
	v[4057] = d[7] + xPi[1];
	v[4056] = d[6] + xPi[0];
	v[4052] = 1e0 / (*lambda2);
	v[4051] = 0.6283185307179586e1;
	v[4050] = v[4051] / (*lambda1);
	v[4049] = (*a4)*d[11] + (*a6)*domegaiP[2] + (*a5)*omegaiP[2];
	v[4048] = (*a4)*d[10] + (*a6)*domegaiP[1] + (*a5)*omegaiP[1];
	v[4047] = (*a4)*d[9] + (*a6)*domegaiP[0] + (*a5)*omegaiP[0];
	v[4034] = Power(d[11], 2);
	v[4033] = 0.5e0*d[11];
	v[4032] = 2e0*d[11];
	v[4031] = Power(d[10], 2);
	v[4030] = 0.5e0*d[10];
	v[4029] = 2e0*d[10];
	v[4028] = Power(d[9], 2);
	v[4027] = 2e0*d[9];
	v[4026] = 0.5e0*d[9];
	v[4019] = Power(d[5], 2);
	v[4018] = 0.5e0*d[5];
	v[4017] = 2e0*d[5];
	v[4016] = Power(d[4], 2);
	v[4015] = 0.5e0*d[4];
	v[4014] = 2e0*d[4];
	v[4013] = Power(d[3], 2);
	v[4012] = 2e0*d[3];
	v[4011] = 0.5e0*d[3];
	v[268] = d[4] * v[4011];
	v[1329] = -v[4013] - v[4016];
	v[4085] = 0.5e0*v[1329];
	v[1313] = d[5] + v[268];
	v[1294] = -d[5] + v[268];
	v[275] = d[5] * v[4015];
	v[1331] = d[3] + v[275];
	v[1311] = -d[3] + v[275];
	v[273] = d[5] * v[4011];
	v[1333] = -d[4] + v[273];
	v[1291] = d[4] + v[273];
	v[1309] = -v[4013] - v[4019];
	v[4083] = 0.5e0*v[1309];
	v[1288] = -v[4016] - v[4019];
	v[4082] = 0.5e0*v[1288];
	v[1247] = 4e0 + v[4013] + v[4016] + v[4019];
	v[2988] = 1e0 / Power(v[1247], 3);
	v[4020] = -2e0*v[2988];
	v[2989] = v[4014] * v[4020];
	v[4022] = -4e0*v[2989];
	v[2992] = v[4017] * v[4022];
	v[2987] = v[4012] * v[4020];
	v[4021] = -4e0*v[2987];
	v[2994] = v[4014] * v[4021];
	v[2991] = v[4017] * v[4021];
	v[1249] = 1e0 / Power(v[1247], 2);
	v[4024] = -4e0*v[1249];
	v[2995] = -8e0*v[1249];
	v[2997] = v[2995] + v[4012] * v[4021];
	v[2996] = v[2995] + v[4014] * v[4022];
	v[2993] = v[2995] + 8e0*v[2988] * (v[4017] * v[4017]);
	v[1251] = v[4017] * v[4024];
	v[4023] = -0.5e0*v[1251];
	v[3096] = 2e0*v[1251];
	v[3070] = v[4017] * v[4023];
	v[3033] = -(v[1251] * v[4018]);
	v[3013] = v[1251] * v[4015];
	v[3010] = -(v[1251] * v[4011]);
	v[3006] = v[4012] * v[4023];
	v[3004] = v[4014] * v[4023];
	v[1250] = v[4014] * v[4024];
	v[4084] = v[1250] + v[3010];
	v[4025] = -0.5e0*v[1250];
	v[3106] = -2e0*v[1250];
	v[3076] = v[4012] * v[4025];
	v[3073] = v[1250] * v[4015];
	v[3052] = v[4014] * v[4025];
	v[3009] = -(v[1250] * v[4011]);
	v[4087] = v[1251] + v[3009];
	v[1248] = v[4012] * v[4024];
	v[4089] = 0.5e0*v[1248];
	v[4086] = v[1248] - v[3013];
	v[3084] = 2e0*v[1248];
	v[3057] = -(v[4012] * v[4089]);
	v[249] = d[10] * v[4026];
	v[468] = -v[4028] - v[4031];
	v[4036] = 0.5e0*v[468];
	v[482] = -d[11] + v[249];
	v[478] = d[11] + v[249];
	v[256] = d[11] * v[4030];
	v[473] = -d[9] + v[256];
	v[471] = d[9] + v[256];
	v[254] = d[11] * v[4026];
	v[479] = d[10] + v[254];
	v[472] = -d[10] + v[254];
	v[489] = 4e0 + v[4028] + v[4031] + v[4034];
	v[829] = 1e0 / Power(v[489], 3);
	v[4035] = -2e0*v[829];
	v[831] = v[4032] * v[4035];
	v[2680] = -4e0*v[831];
	v[830] = v[4029] * v[4035];
	v[2679] = -4e0*v[830];
	v[833] = v[2679] * v[4032];
	v[828] = v[4027] * v[4035];
	v[2678] = -4e0*v[828];
	v[838] = v[2678] * v[4029];
	v[832] = v[2678] * v[4032];
	v[562] = 1e0 / Power(v[489], 2);
	v[4037] = -4e0*v[562];
	v[839] = -8e0*v[562];
	v[841] = v[2678] * v[4027] + v[839];
	v[840] = v[2679] * v[4029] + v[839];
	v[834] = v[2680] * v[4032] + v[839];
	v[837] = v[4036] * v[834];
	v[565] = v[4032] * v[4037];
	v[3361] = -(v[4033] * v[565]);
	v[2605] = -0.5e0*v[565];
	v[907] = 2e0*v[565];
	v[908] = v[482] * v[834] - v[907];
	v[878] = v[478] * v[834] + v[907];
	v[872] = v[2605] * v[4032];
	v[849] = v[2605] * v[4029];
	v[844] = v[2605] * v[4027];
	v[836] = v[4036] * v[833] + v[849];
	v[835] = v[4036] * v[832] + v[844];
	v[598] = -(v[2605] * v[468]);
	v[564] = v[4029] * v[4037];
	v[2604] = -0.5e0*v[564];
	v[900] = v[4030] * v[564];
	v[880] = -2e0*v[564];
	v[881] = v[479] * v[840] - v[880];
	v[875] = v[4033] * v[564];
	v[926] = v[473] * v[834] - v[849];
	v[923] = v[471] * v[834] - v[849];
	v[899] = v[473] * v[840] - v[849];
	v[885] = v[471] * v[840] - v[849];
	v[862] = v[472] * v[840] + v[880];
	v[853] = v[2604] * v[4029];
	v[563] = v[4027] * v[4037];
	v[4042] = -v[563] + v[875];
	v[4038] = v[563] + v[875];
	v[2603] = -0.5e0*v[563];
	v[905] = v[4042] + v[482] * v[832];
	v[903] = v[4026] * v[563];
	v[879] = v[4038] + v[479] * v[838];
	v[876] = v[4038] + v[478] * v[832];
	v[867] = v[4033] * v[563];
	v[4044] = v[564] + v[867];
	v[4039] = -v[564] + v[867];
	v[919] = v[479] * v[834] - v[844];
	v[910] = v[472] * v[834] - v[844];
	v[906] = v[4039] + v[482] * v[833];
	v[902] = v[479] * v[841] - v[844];
	v[893] = v[472] * v[841] - v[844];
	v[877] = v[4044] + v[478] * v[833];
	v[868] = v[4039] + v[473] * v[838];
	v[865] = 2e0*v[563];
	v[866] = v[473] * v[841] - v[865];
	v[863] = v[4030] * v[563];
	v[4043] = v[565] + v[863];
	v[4041] = -v[565] + v[863];
	v[4040] = 2e0*v[863];
	v[934] = v[4040] + v[482] * v[841];
	v[931] = v[4040] + v[478] * v[841];
	v[921] = v[4040] + v[482] * v[840];
	v[912] = v[4040] + v[478] * v[840];
	v[882] = v[4043] + v[479] * v[833];
	v[869] = v[4041] + v[473] * v[832];
	v[864] = v[4041] + v[472] * v[833];
	v[861] = v[4042] + v[472] * v[838];
	v[860] = v[4043] + v[471] * v[832];
	v[859] = v[4044] + v[471] * v[838];
	v[858] = v[471] * v[841] + v[865];
	v[856] = v[2603] * v[4027];
	v[852] = -2e0*v[4040] + v[4036] * v[838];
	v[485] = -v[4031] - v[4034];
	v[4045] = 0.5e0*v[485];
	v[942] = v[4045] * v[833] + 2e0*v[849];
	v[845] = v[4045] * v[832] + v[844];
	v[843] = -v[4040] + v[4045] * v[838];
	v[842] = v[4045] * v[841];
	v[566] = v[4045] * v[563];
	v[476] = -v[4028] - v[4034];
	v[4046] = 0.5e0*v[476];
	v[871] = v[4046] * v[832] + 2e0*v[844];
	v[850] = v[4046] * v[833] + v[849];
	v[848] = v[4046] * v[840];
	v[847] = -v[4040] + v[4046] * v[838];
	v[582] = v[4046] * v[564];
	v[3680] = v[4047] / 2e0;
	v[3677] = -v[4048] / 2e0;
	v[3655] = v[4049] / 2e0;
	v[301] = v[4051] * v[4052];
	v[337] = (*phi1) + cp[0] * v[4050];
	v[4078] = sin(v[337]);
	v[4079] = (*A2) + (*A12)*v[4078];
	v[4053] = cos(v[337]);
	v[1726] = -(Power(v[4050], 3)*v[4053]);
	v[717] = -((v[4050] * v[4050])*v[4078]);
	v[338] = v[4050] * v[4053];
	v[310] = (*phi1) + ci[0] * v[4050];
	v[311] = v[4050] * cos(v[310]);
	v[347] = (*phi2) + cp[1] * v[301];
	v[4080] = sin(v[347]);
	v[4081] = (*A1) + (*A12)*v[4080];
	v[4054] = cos(v[347]);
	v[1728] = -(Power(v[301], 3)*v[4054]);
	v[718] = -((v[301] * v[301])*v[4080]);
	v[1730] = (*A12)*v[338] * v[718];
	v[348] = v[301] * v[4054];
	v[4055] = (*A12)*v[348];
	v[1729] = v[4055] * v[717];
	v[376] = v[338] * v[4055];
	v[320] = (*phi2) + ci[1] * v[301];
	v[321] = v[301] * cos(v[320]);
	v[946] = e1i[2] * v[882] + e1i[1] * v[906] + e1i[0] * v[942];
	v[938] = e1i[0] * v[842] + e1i[2] * v[902] + e1i[1] * v[934];
	v[929] = e1i[2] * v[869] + e1i[1] * v[871] + e1i[0] * v[876];
	v[925] = e1i[2] * v[837] + e1i[0] * v[910] + e1i[1] * v[923];
	v[917] = e1i[1] * v[848] + e1i[2] * v[899] + e1i[0] * v[912];
	v[890] = e1i[2] * v[852] + e1i[1] * v[859] + e1i[0] * v[861];
	v[944] = e2i[2] * v[882] + e2i[1] * v[906] + e2i[0] * v[942];
	v[935] = e2i[0] * v[842] + e2i[2] * v[902] + e2i[1] * v[934];
	v[927] = e2i[2] * v[869] + e2i[1] * v[871] + e2i[0] * v[876];
	v[924] = e2i[2] * v[837] + e2i[0] * v[910] + e2i[1] * v[923];
	v[914] = e2i[1] * v[848] + e2i[2] * v[899] + e2i[0] * v[912];
	v[887] = e2i[2] * v[852] + e2i[1] * v[859] + e2i[0] * v[861];
	v[243] = 4e0 / v[489];
	v[3365] = (v[243] * v[243]);
	v[4354] = -(v[4043] / v[3365]);
	v[4352] = -(v[4039] / v[3365]);
	v[883] = -0.5e0*v[243];
	v[909] = v[883] - v[903];
	v[4062] = v[900] - v[909];
	v[920] = v[4062] + v[482] * v[838];
	v[939] = e1i[0] * v[843] + e1i[2] * v[879] + e1i[1] * v[920];
	v[936] = e2i[0] * v[843] + e2i[2] * v[879] + e2i[1] * v[920];
	v[911] = v[4062] + v[478] * v[838];
	v[916] = e1i[1] * v[847] + e1i[2] * v[868] + e1i[0] * v[911];
	v[913] = e2i[1] * v[847] + e2i[2] * v[868] + e2i[0] * v[911];
	v[884] = v[3361] + v[883];
	v[4064] = -v[884] + v[900];
	v[4063] = -v[884] + v[903];
	v[904] = v[4063] + v[479] * v[832];
	v[940] = e1i[0] * v[845] + e1i[2] * v[904] + e1i[1] * v[905];
	v[937] = e2i[0] * v[845] + e2i[2] * v[904] + e2i[1] * v[905];
	v[901] = v[4064] + v[473] * v[833];
	v[918] = e1i[1] * v[850] + e1i[0] * v[877] + e1i[2] * v[901];
	v[915] = e2i[1] * v[850] + e2i[0] * v[877] + e2i[2] * v[901];
	v[894] = v[4063] + v[472] * v[832];
	v[898] = e1i[2] * v[835] + e1i[1] * v[860] + e1i[0] * v[894];
	v[896] = e2i[2] * v[835] + e2i[1] * v[860] + e2i[0] * v[894];
	v[886] = v[4064] + v[471] * v[833];
	v[892] = e1i[2] * v[836] + e1i[0] * v[864] + e1i[1] * v[886];
	v[889] = e2i[2] * v[836] + e2i[0] * v[864] + e2i[1] * v[886];
	v[870] = -v[243] + v[872];
	v[4065] = v[870] + v[872];
	v[943] = v[4065] + v[4045] * v[834];
	v[947] = e1i[1] * v[908] + e1i[2] * v[919] + e1i[0] * v[943];
	v[945] = e2i[1] * v[908] + e2i[2] * v[919] + e2i[0] * v[943];
	v[873] = v[4065] + v[4046] * v[834];
	v[930] = e1i[1] * v[873] + e1i[0] * v[878] + e1i[2] * v[926];
	v[928] = e2i[1] * v[873] + e2i[0] * v[878] + e2i[2] * v[926];
	v[855] = -v[243] + v[856];
	v[4066] = v[855] + v[856];
	v[874] = v[4066] + v[4046] * v[841];
	v[933] = e1i[2] * v[866] + e1i[1] * v[874] + e1i[0] * v[931];
	v[932] = e2i[2] * v[866] + e2i[1] * v[874] + e2i[0] * v[931];
	v[857] = v[4066] + v[4036] * v[841];
	v[897] = e1i[2] * v[857] + e1i[1] * v[858] + e1i[0] * v[893];
	v[895] = e2i[2] * v[857] + e2i[1] * v[858] + e2i[0] * v[893];
	v[851] = -v[243] + v[853];
	v[4067] = v[851] + v[853];
	v[948] = v[4067] + v[4045] * v[840];
	v[950] = e1i[2] * v[881] + e1i[1] * v[921] + e1i[0] * v[948];
	v[949] = e2i[2] * v[881] + e2i[1] * v[921] + e2i[0] * v[948];
	v[854] = v[4067] + v[4036] * v[840];
	v[891] = e1i[2] * v[854] + e1i[0] * v[862] + e1i[1] * v[885];
	v[888] = e2i[2] * v[854] + e2i[0] * v[862] + e2i[1] * v[885];
	v[596] = v[4029] * v[883];
	v[597] = -(v[2604] * v[468]) + v[596];
	v[594] = v[4027] * v[883];
	v[595] = -(v[2603] * v[468]) + v[594];
	v[591] = v[243] + v[471] * v[563];
	v[589] = -v[243] + v[472] * v[564];
	v[585] = -v[243] + v[473] * v[563];
	v[583] = v[4032] * v[883];
	v[584] = v[4046] * v[565] + v[583];
	v[581] = v[4046] * v[563] + v[594];
	v[580] = v[243] + v[478] * v[565];
	v[576] = v[243] + v[479] * v[564];
	v[574] = -(v[243] * v[4033]);
	v[592] = v[471] * v[564] - v[574];
	v[615] = e2i[0] * v[589] + e2i[1] * v[592] + e2i[2] * v[597];
	v[606] = e1i[0] * v[589] + e1i[1] * v[592] + e1i[2] * v[597];
	v[588] = v[472] * v[563] - v[574];
	v[614] = e2i[0] * v[588] + e2i[1] * v[591] + e2i[2] * v[595];
	v[605] = e1i[0] * v[588] + e1i[1] * v[591] + e1i[2] * v[595];
	v[586] = v[473] * v[564] - v[574];
	v[575] = v[479] * v[563] - v[574];
	v[573] = -v[243] + v[482] * v[565];
	v[571] = -(v[243] * v[4026]);
	v[590] = v[472] * v[565] - v[571];
	v[579] = v[478] * v[564] - v[571];
	v[612] = e2i[0] * v[579] + e2i[1] * v[582] + e2i[2] * v[586];
	v[4070] = 2e0*v[612];
	v[603] = e1i[0] * v[579] + e1i[1] * v[582] + e1i[2] * v[586];
	v[4071] = 2e0*v[603];
	v[577] = v[479] * v[565] - v[571];
	v[572] = v[482] * v[564] - v[571];
	v[569] = v[243] * v[4030];
	v[593] = v[471] * v[565] + v[569];
	v[616] = e2i[0] * v[590] + e2i[1] * v[593] + e2i[2] * v[598];
	v[607] = e1i[0] * v[590] + e1i[1] * v[593] + e1i[2] * v[598];
	v[587] = v[473] * v[565] + v[569];
	v[613] = e2i[0] * v[580] + e2i[1] * v[584] + e2i[2] * v[587];
	v[4072] = 2e0*v[613];
	v[604] = e1i[0] * v[580] + e1i[1] * v[584] + e1i[2] * v[587];
	v[4073] = 2e0*v[604];
	v[578] = v[478] * v[563] + v[569];
	v[611] = e2i[0] * v[578] + e2i[1] * v[581] + e2i[2] * v[585];
	v[4068] = 2e0*v[611];
	v[602] = e1i[0] * v[578] + e1i[1] * v[581] + e1i[2] * v[585];
	v[4069] = 2e0*v[602];
	v[570] = v[482] * v[563] + v[569];
	v[608] = e2i[0] * v[566] + e2i[1] * v[570] + e2i[2] * v[575];
	v[599] = e1i[0] * v[566] + e1i[1] * v[570] + e1i[2] * v[575];
	v[568] = v[4045] * v[565] + v[583];
	v[610] = e2i[0] * v[568] + e2i[1] * v[573] + e2i[2] * v[577];
	v[601] = e1i[0] * v[568] + e1i[1] * v[573] + e1i[2] * v[577];
	v[567] = v[4045] * v[564] + v[596];
	v[609] = e2i[0] * v[567] + e2i[1] * v[572] + e2i[2] * v[576];
	v[600] = e1i[0] * v[567] + e1i[1] * v[572] + e1i[2] * v[576];
	v[246] = 1e0 - v[485] * v[883];
	v[247] = v[243] * v[482];
	v[248] = v[243] * v[479];
	v[250] = v[243] * v[478];
	v[252] = 1e0 - v[476] * v[883];
	v[253] = v[243] * v[473];
	v[255] = v[243] * v[472];
	v[257] = v[243] * v[471];
	v[258] = 1e0 - v[468] * v[883];
	v[262] = 4e0 / v[1247];
	v[4288] = (*radius)*v[262];
	v[3011] = -0.5e0*v[262];
	v[3012] = -v[3011] + v[3073];
	v[3008] = v[3011] - v[1248] * v[4011];
	v[4088] = v[3008] + v[3033];
	v[3007] = -v[262] + v[3070];
	v[3005] = -v[262] + v[3057];
	v[3003] = -v[262] + v[3052];
	v[1336] = v[3011] * v[4014];
	v[1328] = v[3011] * v[4012];
	v[1320] = v[3011] * v[4017];
	v[1299] = -(v[262] * v[4011]);
	v[1293] = v[262] * v[4015];
	v[1290] = -(v[262] * v[4018]);
	v[281] = e1i[0] * v[246] + e1i[1] * v[247] + e1i[2] * v[248];
	v[282] = e1i[0] * v[250] + e1i[1] * v[252] + e1i[2] * v[253];
	v[283] = e1i[0] * v[255] + e1i[1] * v[257] + e1i[2] * v[258];
	v[284] = e2i[0] * v[246] + e2i[1] * v[247] + e2i[2] * v[248];
	v[285] = e2i[0] * v[250] + e2i[1] * v[252] + e2i[2] * v[253];
	v[956] = v[4068] * v[599] - v[4069] * v[608] + v[281] * v[932] - v[284] * v[933] - v[282] * v[935] + v[285] * v[938];
	v[955] = v[4070] * v[600] - v[4071] * v[609] + v[281] * v[914] - v[284] * v[917] - v[282] * v[949] + v[285] * v[950];
	v[954] = -(v[603] * v[608]) - v[602] * v[609] + v[600] * v[611] + v[599] * v[612] + v[281] * v[913] - v[284] * v[916]
		- v[282] * v[936] + v[285] * v[939];
	v[953] = v[4072] * v[601] - v[4073] * v[610] + v[281] * v[928] - v[284] * v[930] - v[282] * v[945] + v[285] * v[947];
	v[952] = -(v[604] * v[609]) - v[603] * v[610] + v[601] * v[612] + v[600] * v[613] + v[281] * v[915] - v[284] * v[918]
		- v[282] * v[944] + v[285] * v[946];
	v[951] = -(v[604] * v[608]) - v[602] * v[610] + v[601] * v[611] + v[599] * v[613] + v[281] * v[927] - v[284] * v[929]
		- v[282] * v[937] + v[285] * v[940];
	v[625] = v[285] * v[601] - v[284] * v[604] - v[282] * v[610] + v[281] * v[613];
	v[721] = v[376] * v[625];
	v[624] = v[285] * v[600] - v[284] * v[603] - v[282] * v[609] + v[281] * v[612];
	v[725] = v[376] * v[624];
	v[623] = v[285] * v[599] - v[284] * v[602] - v[282] * v[608] + v[281] * v[611];
	v[729] = v[376] * v[623];
	v[286] = e2i[0] * v[255] + e2i[1] * v[257] + e2i[2] * v[258];
	v[968] = -(v[4068] * v[605]) + v[4069] * v[614] + v[282] * v[895] - v[285] * v[897] - v[283] * v[932] + v[286] * v[933];
	v[967] = -(v[4070] * v[606]) + v[4071] * v[615] + v[282] * v[888] - v[285] * v[891] - v[283] * v[914] + v[286] * v[917];
	v[966] = -(v[606] * v[611]) - v[605] * v[612] + v[603] * v[614] + v[602] * v[615] + v[282] * v[887] - v[285] * v[890]
		- v[283] * v[913] + v[286] * v[916];
	v[965] = -(v[4072] * v[607]) + v[4073] * v[616] + v[282] * v[924] - v[285] * v[925] - v[283] * v[928] + v[286] * v[930];
	v[964] = -(v[607] * v[612]) - v[606] * v[613] + v[604] * v[615] + v[603] * v[616] + v[282] * v[889] - v[285] * v[892]
		- v[283] * v[915] + v[286] * v[918];
	v[963] = -(v[607] * v[611]) - v[605] * v[613] + v[604] * v[614] + v[602] * v[616] + v[282] * v[896] - v[285] * v[898]
		- v[283] * v[927] + v[286] * v[929];
	v[962] = 2e0*v[605] * v[608] - 2e0*v[599] * v[614] - v[281] * v[895] + v[284] * v[897] + v[283] * v[935] - v[286] * v[938];
	v[961] = 2e0*v[606] * v[609] - 2e0*v[600] * v[615] - v[281] * v[888] + v[284] * v[891] + v[283] * v[949] - v[286] * v[950];
	v[960] = v[606] * v[608] + v[605] * v[609] - v[600] * v[614] - v[599] * v[615] - v[281] * v[887] + v[284] * v[890]
		+ v[283] * v[936] - v[286] * v[939];
	v[959] = 2e0*v[607] * v[610] - 2e0*v[601] * v[616] - v[281] * v[924] + v[284] * v[925] + v[283] * v[945] - v[286] * v[947];
	v[958] = v[607] * v[609] + v[606] * v[610] - v[601] * v[615] - v[600] * v[616] - v[281] * v[889] + v[284] * v[892]
		+ v[283] * v[944] - v[286] * v[946];
	v[957] = v[607] * v[608] + v[605] * v[610] - v[601] * v[614] - v[599] * v[616] - v[281] * v[896] + v[284] * v[898]
		+ v[283] * v[937] - v[286] * v[940];
	v[622] = -(v[286] * v[601]) + v[284] * v[607] + v[283] * v[610] - v[281] * v[616];
	v[733] = v[376] * v[622];
	v[621] = -(v[286] * v[600]) + v[284] * v[606] + v[283] * v[609] - v[281] * v[615];
	v[737] = v[376] * v[621];
	v[620] = -(v[286] * v[599]) + v[284] * v[605] + v[283] * v[608] - v[281] * v[614];
	v[741] = v[376] * v[620];
	v[619] = v[286] * v[604] - v[285] * v[607] - v[283] * v[613] + v[282] * v[616];
	v[745] = v[376] * v[619];
	v[618] = v[286] * v[603] - v[285] * v[606] - v[283] * v[612] + v[282] * v[615];
	v[749] = v[376] * v[618];
	v[617] = v[286] * v[602] - v[285] * v[605] - v[283] * v[611] + v[282] * v[614];
	v[753] = v[376] * v[617];
	v[287] = -(v[283] * v[285]) + v[282] * v[286];
	v[1732] = v[1730] * v[287];
	v[1731] = v[1729] * v[287];
	v[386] = v[287] * v[376];
	v[288] = v[283] * v[284] - v[281] * v[286];
	v[1734] = v[1730] * v[288];
	v[1733] = v[1729] * v[288];
	v[388] = v[288] * v[376];
	v[289] = -(v[282] * v[284]) + v[281] * v[285];
	v[1736] = v[1730] * v[289];
	v[1735] = v[1729] * v[289];
	v[390] = v[289] * v[376];
	v[1780] = 2e0*((v[386] * v[386]) + (v[388] * v[388]) + (v[390] * v[390]));
	v[293] = sin(v[310]);
	v[4077] = (*A2) + (*A12)*v[293];
	v[322] = v[321] * v[4077];
	v[294] = sin(v[320]);
	v[312] = ((*A1) + (*A12)*v[294])*v[311];
	v[296] = (*A1)*v[293] + v[294] * v[4077];
	v[3032] = ci[1] * v[935] + ci[0] * v[938] + v[296] * v[968];
	v[3031] = ci[1] * v[949] + ci[0] * v[950] + v[296] * v[967];
	v[3030] = ci[1] * v[936] + ci[0] * v[939] + v[296] * v[966];
	v[3029] = ci[1] * v[945] + ci[0] * v[947] + v[296] * v[965];
	v[3028] = ci[1] * v[944] + ci[0] * v[946] + v[296] * v[964];
	v[3027] = ci[1] * v[937] + ci[0] * v[940] + v[296] * v[963];
	v[3026] = ci[1] * v[932] + ci[0] * v[933] + v[296] * v[962];
	v[3025] = ci[1] * v[914] + ci[0] * v[917] + v[296] * v[961];
	v[3024] = ci[1] * v[913] + ci[0] * v[916] + v[296] * v[960];
	v[3023] = ci[1] * v[928] + ci[0] * v[930] + v[296] * v[959];
	v[3022] = ci[1] * v[915] + ci[0] * v[918] + v[296] * v[958];
	v[3021] = ci[1] * v[927] + ci[0] * v[929] + v[296] * v[957];
	v[3020] = ci[1] * v[895] + ci[0] * v[897] + v[296] * v[956];
	v[3019] = ci[1] * v[888] + ci[0] * v[891] + v[296] * v[955];
	v[3018] = ci[1] * v[887] + ci[0] * v[890] + v[296] * v[954];
	v[3017] = ci[1] * v[924] + ci[0] * v[925] + v[296] * v[953];
	v[3016] = ci[1] * v[889] + ci[0] * v[892] + v[296] * v[952];
	v[3015] = ci[1] * v[896] + ci[0] * v[898] + v[296] * v[951];
	v[1347] = ci[0] * v[607] + ci[1] * v[616] + v[296] * v[625];
	v[1346] = ci[0] * v[606] + ci[1] * v[615] + v[296] * v[624];
	v[1345] = ci[0] * v[605] + ci[1] * v[614] + v[296] * v[623];
	v[1327] = ci[0] * v[604] + ci[1] * v[613] + v[296] * v[622];
	v[1326] = ci[0] * v[603] + ci[1] * v[612] + v[296] * v[621];
	v[1325] = ci[0] * v[602] + ci[1] * v[611] + v[296] * v[620];
	v[1308] = ci[0] * v[601] + ci[1] * v[610] + v[296] * v[619];
	v[1307] = ci[0] * v[600] + ci[1] * v[609] + v[296] * v[618];
	v[1306] = ci[0] * v[599] + ci[1] * v[608] + v[296] * v[617];
	v[1737] = v[1728] * v[4079];
	v[1740] = v[1737] * v[287];
	v[1739] = v[1737] * v[288];
	v[1738] = v[1737] * v[289];
	v[377] = v[4079] * v[718];
	v[398] = v[289] * v[377];
	v[397] = v[288] * v[377];
	v[396] = v[287] * v[377];
	v[349] = v[348] * v[4079];
	v[652] = v[616] + v[349] * v[625];
	v[651] = v[615] + v[349] * v[624];
	v[650] = v[614] + v[349] * v[623];
	v[649] = v[613] + v[349] * v[622];
	v[648] = v[612] + v[349] * v[621];
	v[647] = v[611] + v[349] * v[620];
	v[646] = v[610] + v[349] * v[619];
	v[645] = v[609] + v[349] * v[618];
	v[644] = v[608] + v[349] * v[617];
	v[384] = v[286] + v[289] * v[349];
	v[4109] = 2e0*v[384];
	v[382] = v[285] + v[288] * v[349];
	v[4108] = 2e0*v[382];
	v[380] = v[284] + v[287] * v[349];
	v[4107] = 2e0*v[380];
	v[1741] = v[1726] * v[4081];
	v[1744] = v[1741] * v[287];
	v[1743] = v[1741] * v[288];
	v[1742] = v[1741] * v[289];
	v[378] = v[4081] * v[717];
	v[389] = v[289] * v[378];
	v[387] = v[288] * v[378];
	v[385] = v[287] * v[378];
	v[339] = v[338] * v[4081];
	v[634] = v[607] + v[339] * v[625];
	v[633] = v[606] + v[339] * v[624];
	v[632] = v[605] + v[339] * v[623];
	v[631] = v[604] + v[339] * v[622];
	v[630] = v[603] + v[339] * v[621];
	v[629] = v[602] + v[339] * v[620];
	v[628] = v[601] + v[339] * v[619];
	v[627] = v[600] + v[339] * v[618];
	v[626] = v[599] + v[339] * v[617];
	v[383] = v[283] + v[289] * v[339];
	v[4104] = 2e0*v[383];
	v[381] = v[282] + v[288] * v[339];
	v[4105] = 2e0*v[381];
	v[379] = v[281] + v[287] * v[339];
	v[4106] = 2e0*v[379];
	v[304] = (*A1)*v[4078] + v[4079] * v[4080];
	v[688] = cp[0] * v[607] + cp[1] * v[616] + v[304] * v[625];
	v[686] = cp[0] * v[606] + cp[1] * v[615] + v[304] * v[624];
	v[684] = cp[0] * v[605] + cp[1] * v[614] + v[304] * v[623];
	v[682] = cp[0] * v[604] + cp[1] * v[613] + v[304] * v[622];
	v[680] = cp[0] * v[603] + cp[1] * v[612] + v[304] * v[621];
	v[678] = cp[0] * v[602] + cp[1] * v[611] + v[304] * v[620];
	v[676] = cp[0] * v[601] + cp[1] * v[610] + v[304] * v[619];
	v[674] = cp[0] * v[600] + cp[1] * v[609] + v[304] * v[618];
	v[672] = cp[0] * v[599] + cp[1] * v[608] + v[304] * v[617];
	v[303] = cp[0] * v[281] + cp[1] * v[284] + v[287] * v[304] + v[4056];
	v[368] = -v[303] + v[4059];
	v[305] = cp[0] * v[282] + cp[1] * v[285] + v[288] * v[304] + v[4057];
	v[369] = -v[305] + v[4060];
	v[306] = cp[0] * v[283] + cp[1] * v[286] + v[289] * v[304] + v[4058];
	v[370] = -v[306] + v[4061];
	v[313] = e1i[0] + v[312] * v[4074];
	v[314] = e1i[1] + v[312] * v[4075];
	v[315] = e1i[2] + v[312] * v[4076];
	v[317] = 1e0 / sqrt(Power(v[313], 2) + Power(v[314], 2) + Power(v[315], 2));
	v[316] = v[313] * v[317];
	v[318] = v[314] * v[317];
	v[319] = v[315] * v[317];
	v[323] = e2i[0] + v[322] * v[4074];
	v[324] = e2i[1] + v[322] * v[4075];
	v[325] = e2i[2] + v[322] * v[4076];
	v[327] = 1e0 / sqrt(Power(v[323], 2) + Power(v[324], 2) + Power(v[325], 2));
	v[326] = v[323] * v[327];
	v[328] = v[324] * v[327];
	v[335] = -(v[318] * v[326]) + v[316] * v[328];
	v[329] = v[325] * v[327];
	v[333] = v[319] * v[326] - v[316] * v[329];
	v[330] = -(v[319] * v[328]) + v[318] * v[329];
	v[332] = 1e0 / sqrt(Power(v[330], 2) + Power(v[333], 2) + Power(v[335], 2));
	v[331] = v[330] * v[332];
	v[4287] = -(v[1288] * v[331]);
	v[334] = v[332] * v[333];
	v[4289] = v[1309] * v[334];
	v[336] = v[332] * v[335];
	v[4290] = v[1329] * v[336];
	v[3113] = (*radius)*(-((v[1294] * v[2997] - v[3076])*v[334]) - (v[1291] * v[2997] - v[3006])*v[336]
		- v[2997] * v[331] * v[4082]);
	v[3109] = (*radius)*(-((v[1294] * v[2996] - v[3076])*v[334]) - (v[1291] * v[2996] - v[3106])*v[336] - v[331] *
		(v[3003] + v[3052] + v[2996] * v[4082]));
	v[3104] = (*radius)*(-((v[1294] * v[2994] - v[3008] + v[3073])*v[334]) - (v[1248] + v[1291] * v[2994] + v[3013]
		)*v[336] - v[331] * (v[3076] + v[2994] * v[4082]));
	v[3100] = (*radius)*(-((v[1294] * v[2993] - v[3096])*v[334]) - (v[1291] * v[2993] - v[3006])*v[336] - v[331] *
		(v[3007] + v[3070] + v[2993] * v[4082]));
	v[3095] = (*radius)*(-((v[1251] + v[1291] * v[2992] - v[3009])*v[336]) + v[331] * (-2e0*v[3004] - v[2992] * v[4082]
		) + v[334] * (-(v[1294] * v[2992]) + v[4084]));
	v[3091] = (*radius)*(v[331] * (-v[3006] - v[2991] * v[4082]) + v[334] * (-(v[1294] * v[2991]) + v[4086]) + v[336] * (-
		(v[1291] * v[2991]) + v[4088]));
	v[3087] = (*radius)*(-((v[1313] * v[2997] - v[3076])*v[331]) - (v[1311] * v[2997] - v[3084])*v[336] - v[334] *
		(v[3005] + v[3057] + v[2997] * v[4083]));
	v[3082] = (*radius)*(-((v[1313] * v[2996] - v[3076])*v[331]) - (v[1311] * v[2996] - v[3004])*v[336]
		- v[2996] * v[334] * v[4083]);
	v[3078] = (*radius)*(-((v[1313] * v[2994] - v[3008] + v[3073])*v[331]) - v[334] * (v[3076] + v[2994] * v[4083])
		+ v[336] * (-(v[1311] * v[2994]) + v[4084]));
	v[3072] = (*radius)*(-((v[1313] * v[2993] + v[3096])*v[331]) - (v[1311] * v[2993] - v[3004])*v[336] - v[334] *
		(v[3007] + v[3070] + v[2993] * v[4083]));
	v[3067] = (*radius)*(-((v[1250] + v[1313] * v[2992] - v[3010])*v[331]) - (v[1311] * v[2992] + v[3012] - v[3033]
		)*v[336] - v[334] * (v[3004] + v[2992] * v[4083]));
	v[3063] = (*radius)*(-((v[1248] + v[1313] * v[2991] + v[3013])*v[331]) + v[334] * (-2e0*v[3006] - v[2991] * v[4083]
		) + v[336] * (-(v[1311] * v[2991]) + v[4087]));
	v[3059] = (*radius)*(-((v[1333] * v[2997] - v[3006])*v[331]) - (v[1331] * v[2997] + v[3084])*v[334] - v[336] *
		(v[3005] + v[3057] + v[2997] * v[4085]));
	v[3054] = (*radius)*(-((v[1333] * v[2996] + v[3106])*v[331]) - (v[1331] * v[2996] - v[3004])*v[334] - v[336] *
		(v[3003] + v[3052] + v[2996] * v[4085]));
	v[3049] = (*radius)*(-((v[1250] + v[1331] * v[2994] - v[3010])*v[334]) + v[336] * (-2e0*v[3076] - v[2994] * v[4085]
		) + v[331] * (-(v[1333] * v[2994]) + v[4086]));
	v[3045] = (*radius)*(-((v[1333] * v[2993] - v[3006])*v[331]) - (v[1331] * v[2993] - v[3004])*v[334]
		- v[2993] * v[336] * v[4085]);
	v[3041] = (*radius)*(-((v[1331] * v[2992] + v[3012] - v[3033])*v[334]) - v[336] * (v[3004] + v[2992] * v[4085])
		+ v[331] * (-(v[1333] * v[2992]) + v[4087]));
	v[3037] = (*radius)*(-((v[1251] + v[1331] * v[2991] - v[3009])*v[334]) + v[336] * (-v[3006] - v[2991] * v[4085])
		+ v[331] * (-(v[1333] * v[2991]) + v[4088]));
	v[1344] = (*radius)*((v[1299] - v[1251] * v[1333])*v[331] - (v[1293] + v[1251] * v[1331])*v[334]
		+ v[4023] * v[4290]);
	v[1340] = (*radius)*((-(v[1250] * v[1333]) + v[262])*v[331] + (v[1290] - v[1250] * v[1331])*v[334] - v[336] *
		(v[1336] - v[1329] * v[4025]));
	v[1335] = (*radius)*((v[1290] - v[1248] * v[1333])*v[331] - (v[1248] * v[1331] + v[262])*v[334] - 1e0*v[336] *
		(v[1328] + v[1329] * v[4089]));
	v[1324] = (*radius)*(-((v[1251] * v[1313] + v[262])*v[331]) - (v[1293] + v[1251] * v[1311])*v[336] - v[334] *
		(v[1320] - v[1309] * v[4023]));
	v[1319] = (*radius)*((v[1299] - v[1250] * v[1313])*v[331] + (v[1290] - v[1250] * v[1311])*v[336]
		+ v[4025] * v[4289]);
	v[1315] = (*radius)*(-((v[1293] + v[1248] * v[1313])*v[331]) + (-(v[1248] * v[1311]) + v[262])*v[336] - v[334] *
		(v[1328] + v[1309] * v[4089]));
	v[1305] = (*radius)*((-(v[1251] * v[1294]) + v[262])*v[334] - (v[1251] * v[1291] - v[1299])*v[336] - v[331] *
		(v[1320] - v[1288] * v[4023]));
	v[1301] = (*radius)*(-((v[1250] * v[1294] - v[1299])*v[334]) - (v[1250] * v[1291] + v[262])*v[336] - v[331] *
		(v[1336] - v[1288] * v[4025]));
	v[1296] = (*radius)*(-((v[1293] + v[1248] * v[1294])*v[334]) + (v[1290] - v[1248] * v[1291])*v[336]
		+ v[4089] * v[4287]);
	v[526] = 2e0*(v[379] * v[386] + v[381] * v[388] + v[383] * v[390]);
	v[525] = 2e0*(v[379] * v[385] + v[381] * v[387] + v[383] * v[389]);
	v[392] = (v[379] * v[379]) + (v[381] * v[381]) + (v[383] * v[383]);
	v[1748] = 1e0 / sqrt(v[392]);
	v[4094] = v[1748] * v[378];
	v[1751] = (0.75e0*v[1748]) / Power(v[392], 2);
	v[1749] = -v[1748] / (2e0*v[392]);
	v[4153] = 2e0*v[1749];
	v[4091] = 2e0*v[1749];
	v[4090] = 2e0*v[1749];
	v[1753] = (v[1744] * v[379] + v[1743] * v[381] + v[1742] * v[383] + (v[385] * v[385]) + (v[387] * v[387]) +
		(v[389] * v[389]))*v[4090] + v[1751] * (v[525] * v[525]);
	v[1752] = v[1749] * (v[1780] + v[1736] * v[4104] + v[1734] * v[4105] + v[1732] * v[4106]) + v[1751] * (v[526] * v[526]);
	v[1750] = (v[1731] * v[379] + v[1733] * v[381] + v[1735] * v[383] + v[385] * v[386] + v[387] * v[388] + v[389] * v[390]
		)*v[4090] + v[1751] * v[525] * v[526];
	v[528] = v[1749] * v[526];
	v[4093] = 2e0*v[528];
	v[527] = v[1749] * v[525];
	v[4092] = 2e0*v[527];
	v[2278] = v[4091] * (v[368] * v[628] + v[369] * v[631] + v[370] * v[634] - v[379] * v[676] - v[381] * v[682]
		- v[383] * v[688]);
	v[2276] = v[4091] * (v[368] * v[627] + v[369] * v[630] + v[370] * v[633] - v[379] * v[674] - v[381] * v[680]
		- v[383] * v[686]);
	v[2274] = v[4153] * (v[368] * v[626] + v[369] * v[629] + v[370] * v[632] - v[379] * v[672] - v[381] * v[678]
		- v[383] * v[684]);
	v[1775] = v[1742] * v[1748] + v[1753] * v[383] + v[389] * v[4092];
	v[1773] = v[1743] * v[1748] + v[1753] * v[381] + v[387] * v[4092];
	v[1771] = v[1744] * v[1748] + v[1753] * v[379] + v[385] * v[4092];
	v[1777] = v[1771] * v[316] + v[1773] * v[318] + v[1775] * v[319];
	v[1776] = v[1771] * v[331] + v[1773] * v[334] + v[1775] * v[336];
	v[1765] = v[1736] * v[1748] + v[1752] * v[383] + v[390] * v[4093];
	v[1763] = v[1735] * v[1748] + v[1750] * v[383] + v[389] * v[4093];
	v[1761] = v[1734] * v[1748] + v[1752] * v[381] + v[388] * v[4093];
	v[1759] = v[1733] * v[1748] + v[1750] * v[381] + v[387] * v[4093];
	v[1757] = v[1732] * v[1748] + v[1752] * v[379] + v[386] * v[4093];
	v[1769] = v[1757] * v[316] + v[1761] * v[318] + v[1765] * v[319];
	v[1767] = v[1757] * v[331] + v[1761] * v[334] + v[1765] * v[336];
	v[1755] = v[1731] * v[1748] + v[1750] * v[379] + v[385] * v[4093];
	v[1768] = v[1755] * v[316] + v[1759] * v[318] + v[1763] * v[319];
	v[1766] = v[1755] * v[331] + v[1759] * v[334] + v[1763] * v[336];
	v[1004] = v[1748] * (v[938] + v[339] * v[968]);
	v[1002] = v[1748] * (v[950] + v[339] * v[967]);
	v[1000] = v[1748] * (v[939] + v[339] * v[966]);
	v[998] = v[1748] * (v[947] + v[339] * v[965]);
	v[996] = v[1748] * (v[946] + v[339] * v[964]);
	v[994] = v[1748] * (v[940] + v[339] * v[963]);
	v[992] = v[1748] * (v[933] + v[339] * v[962]);
	v[990] = v[1748] * (v[917] + v[339] * v[961]);
	v[988] = v[1748] * (v[916] + v[339] * v[960]);
	v[986] = v[1748] * (v[930] + v[339] * v[959]);
	v[984] = v[1748] * (v[918] + v[339] * v[958]);
	v[982] = v[1748] * (v[929] + v[339] * v[957]);
	v[980] = v[1748] * (v[897] + v[339] * v[956]);
	v[3125] = v[1004] * v[316] + v[319] * v[980] + v[318] * v[992];
	v[3124] = v[1004] * v[331] + v[336] * v[980] + v[334] * v[992];
	v[978] = v[1748] * (v[891] + v[339] * v[955]);
	v[3123] = v[1002] * v[316] + v[319] * v[978] + v[318] * v[990];
	v[3121] = v[1002] * v[331] + v[336] * v[978] + v[334] * v[990];
	v[976] = v[1748] * (v[890] + v[339] * v[954]);
	v[3122] = v[1000] * v[316] + v[319] * v[976] + v[318] * v[988];
	v[3120] = v[1000] * v[331] + v[336] * v[976] + v[334] * v[988];
	v[974] = v[1748] * (v[925] + v[339] * v[953]);
	v[3119] = v[319] * v[974] + v[318] * v[986] + v[316] * v[998];
	v[3116] = v[336] * v[974] + v[334] * v[986] + v[331] * v[998];
	v[972] = v[1748] * (v[892] + v[339] * v[952]);
	v[3118] = v[319] * v[972] + v[318] * v[984] + v[316] * v[996];
	v[3115] = v[336] * v[972] + v[334] * v[984] + v[331] * v[996];
	v[970] = v[1748] * (v[898] + v[339] * v[951]);
	v[3117] = v[319] * v[970] + v[318] * v[982] + v[316] * v[994];
	v[3114] = v[336] * v[970] + v[334] * v[982] + v[331] * v[994];
	v[754] = v[528] * v[626] + v[1748] * v[753];
	v[752] = v[4094] * v[617] + v[527] * v[626];
	v[750] = v[528] * v[627] + v[1748] * v[749];
	v[748] = v[4094] * v[618] + v[527] * v[627];
	v[746] = v[528] * v[628] + v[1748] * v[745];
	v[744] = v[4094] * v[619] + v[527] * v[628];
	v[742] = v[528] * v[629] + v[1748] * v[741];
	v[740] = v[4094] * v[620] + v[527] * v[629];
	v[738] = v[528] * v[630] + v[1748] * v[737];
	v[736] = v[4094] * v[621] + v[527] * v[630];
	v[734] = v[528] * v[631] + v[1748] * v[733];
	v[732] = v[4094] * v[622] + v[527] * v[631];
	v[730] = v[528] * v[632] + v[1748] * v[729];
	v[2171] = v[319] * v[730] + v[318] * v[742] + v[316] * v[754];
	v[2168] = v[336] * v[730] + v[334] * v[742] + v[331] * v[754];
	v[728] = v[4094] * v[623] + v[527] * v[632];
	v[2177] = v[319] * v[728] + v[318] * v[740] + v[316] * v[752];
	v[2174] = v[336] * v[728] + v[334] * v[740] + v[331] * v[752];
	v[726] = v[528] * v[633] + v[1748] * v[725];
	v[2172] = v[319] * v[726] + v[318] * v[738] + v[316] * v[750];
	v[2169] = v[336] * v[726] + v[334] * v[738] + v[331] * v[750];
	v[724] = v[4094] * v[624] + v[527] * v[633];
	v[2178] = v[319] * v[724] + v[318] * v[736] + v[316] * v[748];
	v[2175] = v[336] * v[724] + v[334] * v[736] + v[331] * v[748];
	v[722] = v[528] * v[634] + v[1748] * v[721];
	v[2173] = v[319] * v[722] + v[318] * v[734] + v[316] * v[746];
	v[2170] = v[336] * v[722] + v[334] * v[734] + v[331] * v[746];
	v[720] = v[4094] * v[625] + v[527] * v[634];
	v[2179] = v[319] * v[720] + v[318] * v[732] + v[316] * v[744];
	v[2176] = v[336] * v[720] + v[334] * v[732] + v[331] * v[744];
	v[643] = v[1748] * v[634];
	v[4146] = 2e0*v[643];
	v[642] = v[1748] * v[633];
	v[4145] = 2e0*v[642];
	v[641] = v[1748] * v[632];
	v[4144] = 2e0*v[641];
	v[640] = v[1748] * v[631];
	v[4142] = 2e0*v[640];
	v[639] = v[1748] * v[630];
	v[4141] = 2e0*v[639];
	v[638] = v[1748] * v[629];
	v[4140] = 2e0*v[638];
	v[637] = v[1748] * v[628];
	v[4134] = 2e0*v[637];
	v[1281] = v[331] * v[637] + v[334] * v[640] + v[336] * v[643];
	v[1263] = v[316] * v[637] + v[318] * v[640] + v[319] * v[643];
	v[636] = v[1748] * v[627];
	v[4136] = 2e0*v[636];
	v[1280] = v[331] * v[636] + v[334] * v[639] + v[336] * v[642];
	v[1262] = v[316] * v[636] + v[318] * v[639] + v[319] * v[642];
	v[635] = v[1748] * v[626];
	v[4138] = 2e0*v[635];
	v[1279] = v[331] * v[635] + v[334] * v[638] + v[336] * v[641];
	v[1261] = v[316] * v[635] + v[318] * v[638] + v[319] * v[641];
	v[409] = v[1748] * v[386] + v[379] * v[528];
	v[408] = v[1748] * v[388] + v[381] * v[528];
	v[4151] = 2e0*v[408];
	v[407] = v[1748] * v[390] + v[383] * v[528];
	v[1210] = v[336] * v[407] + v[334] * v[408] + v[331] * v[409];
	v[1198] = v[319] * v[407] + v[318] * v[408] + v[316] * v[409];
	v[406] = v[1748] * v[385] + v[379] * v[527];
	v[405] = v[1748] * v[387] + v[381] * v[527];
	v[4149] = 2e0*v[405];
	v[404] = v[1748] * v[389] + v[383] * v[527];
	v[1209] = v[336] * v[404] + v[334] * v[405] + v[331] * v[406];
	v[1197] = v[319] * v[404] + v[318] * v[405] + v[316] * v[406];
	v[343] = v[1748] * v[379];
	v[4101] = -(v[343] * v[626]);
	v[4098] = -(v[343] * v[627]);
	v[4095] = -(v[343] * v[628]);
	v[345] = v[1748] * v[381];
	v[4102] = -(v[345] * v[629]);
	v[4099] = -(v[345] * v[630]);
	v[4096] = -(v[345] * v[631]);
	v[346] = v[1748] * v[383];
	v[4103] = -(v[346] * v[632]);
	v[4100] = -(v[346] * v[633]);
	v[4097] = -(v[346] * v[634]);
	v[2257] = -(v[380] * v[637]) - v[382] * v[640] - v[384] * v[643] - v[343] * v[646] - v[345] * v[649] - v[346] * v[652]
		- v[409] * v[676] - v[408] * v[682] - v[407] * v[688] + v[370] * v[722] + v[369] * v[734] + v[368] * v[746];
	v[2256] = -(v[380] * v[636]) - v[382] * v[639] - v[384] * v[642] - v[343] * v[645] - v[345] * v[648] - v[346] * v[651]
		- v[409] * v[674] - v[408] * v[680] - v[407] * v[686] + v[370] * v[726] + v[369] * v[738] + v[368] * v[750];
	v[2255] = -(v[380] * v[635]) - v[382] * v[638] - v[384] * v[641] - v[343] * v[644] - v[345] * v[647] - v[346] * v[650]
		- v[409] * v[672] - v[408] * v[678] - v[407] * v[684] + v[370] * v[730] + v[369] * v[742] + v[368] * v[754];
	v[2251] = 2e0*v[4095] + 2e0*v[4096] + 2e0*v[4097] - v[406] * v[676] - v[405] * v[682] - v[404] * v[688] + v[370] * v[720]
		+ v[369] * v[732] + v[368] * v[744];
	v[2250] = 2e0*v[4098] + 2e0*v[4099] + 2e0*v[4100] - v[406] * v[674] - v[405] * v[680] - v[404] * v[686] + v[370] * v[724]
		+ v[369] * v[736] + v[368] * v[748];
	v[2249] = 2e0*v[4101] + 2e0*v[4102] + 2e0*v[4103] - v[406] * v[672] - v[405] * v[678] - v[404] * v[684] + v[370] * v[728]
		+ v[369] * v[740] + v[368] * v[752];
	v[1864] = v[1757] * v[368] + v[1761] * v[369] + v[1765] * v[370] - v[343] * v[396] - v[345] * v[397] - v[346] * v[398]
		- v[409] * v[4107] - v[408] * v[4108] - v[407] * v[4109];
	v[1861] = v[1755] * v[368] + v[1759] * v[369] + v[1763] * v[370] - v[343] * v[386] - v[345] * v[388] - v[346] * v[390]
		- v[384] * v[404] - v[382] * v[405] - v[380] * v[406] - v[383] * v[407] - v[381] * v[408] - v[379] * v[409];
	v[1860] = v[1771] * v[368] + v[1773] * v[369] + v[1775] * v[370] - v[343] * v[385] - v[345] * v[387] - v[346] * v[389]
		- v[404] * v[4104] - v[405] * v[4105] - v[406] * v[4106];
	v[530] = 2e0*(v[380] * v[396] + v[382] * v[397] + v[384] * v[398]);
	v[529] = 2e0*(v[380] * v[386] + v[382] * v[388] + v[384] * v[390]);
	v[400] = (v[380] * v[380]) + (v[382] * v[382]) + (v[384] * v[384]);
	v[1782] = 1e0 / sqrt(v[400]);
	v[4114] = v[1782] * v[377];
	v[1785] = (0.75e0*v[1782]) / Power(v[400], 2);
	v[1783] = -v[1782] / (2e0*v[400]);
	v[4152] = 2e0*v[1783];
	v[4111] = 2e0*v[1783];
	v[4110] = 2e0*v[1783];
	v[1787] = v[1783] * (v[1780] + v[1731] * v[4107] + v[1733] * v[4108] + v[1735] * v[4109]) + v[1785] * (v[529] * v[529]);
	v[1786] = (v[1740] * v[380] + v[1739] * v[382] + v[1738] * v[384] + (v[396] * v[396]) + (v[397] * v[397]) +
		(v[398] * v[398]))*v[4110] + v[1785] * (v[530] * v[530]);
	v[1784] = (v[1732] * v[380] + v[1734] * v[382] + v[1736] * v[384] + v[386] * v[396] + v[388] * v[397] + v[390] * v[398]
		)*v[4110] + v[1785] * v[529] * v[530];
	v[532] = v[1783] * v[530];
	v[4113] = 2e0*v[532];
	v[531] = v[1783] * v[529];
	v[4112] = 2e0*v[531];
	v[2269] = v[4111] * (v[368] * v[646] + v[369] * v[649] + v[370] * v[652] - v[380] * v[676] - v[382] * v[682]
		- v[384] * v[688]);
	v[2267] = v[4111] * (v[368] * v[645] + v[369] * v[648] + v[370] * v[651] - v[380] * v[674] - v[382] * v[680]
		- v[384] * v[686]);
	v[2265] = v[4152] * (v[368] * v[644] + v[369] * v[647] + v[370] * v[650] - v[380] * v[672] - v[382] * v[678]
		- v[384] * v[684]);
	v[1805] = v[1735] * v[1782] + v[1787] * v[384] + v[390] * v[4112];
	v[1803] = v[1733] * v[1782] + v[1787] * v[382] + v[388] * v[4112];
	v[1801] = v[1731] * v[1782] + v[1787] * v[380] + v[386] * v[4112];
	v[1799] = v[1738] * v[1782] + v[1786] * v[384] + v[398] * v[4113];
	v[1797] = v[1736] * v[1782] + v[1784] * v[384] + v[390] * v[4113];
	v[1795] = v[1739] * v[1782] + v[1786] * v[382] + v[397] * v[4113];
	v[1793] = v[1734] * v[1782] + v[1784] * v[382] + v[388] * v[4113];
	v[1791] = v[1740] * v[1782] + v[1786] * v[380] + v[396] * v[4113];
	v[1789] = v[1732] * v[1782] + v[1784] * v[380] + v[386] * v[4113];
	v[1040] = v[1782] * (v[935] + v[349] * v[968]);
	v[1038] = v[1782] * (v[949] + v[349] * v[967]);
	v[1036] = v[1782] * (v[936] + v[349] * v[966]);
	v[1034] = v[1782] * (v[945] + v[349] * v[965]);
	v[1032] = v[1782] * (v[944] + v[349] * v[964]);
	v[1030] = v[1782] * (v[937] + v[349] * v[963]);
	v[1028] = v[1782] * (v[932] + v[349] * v[962]);
	v[1026] = v[1782] * (v[914] + v[349] * v[961]);
	v[1024] = v[1782] * (v[913] + v[349] * v[960]);
	v[1022] = v[1782] * (v[928] + v[349] * v[959]);
	v[1020] = v[1782] * (v[915] + v[349] * v[958]);
	v[1018] = v[1782] * (v[927] + v[349] * v[957]);
	v[1016] = v[1782] * (v[895] + v[349] * v[956]);
	v[1014] = v[1782] * (v[888] + v[349] * v[955]);
	v[1012] = v[1782] * (v[887] + v[349] * v[954]);
	v[1010] = v[1782] * (v[924] + v[349] * v[953]);
	v[1008] = v[1782] * (v[889] + v[349] * v[952]);
	v[1006] = v[1782] * (v[896] + v[349] * v[951]);
	v[781] = v[4114] * v[617] + v[532] * v[644];
	v[779] = v[531] * v[644] + v[1782] * v[753];
	v[778] = v[4114] * v[618] + v[532] * v[645];
	v[776] = v[531] * v[645] + v[1782] * v[749];
	v[775] = v[4114] * v[619] + v[532] * v[646];
	v[773] = v[531] * v[646] + v[1782] * v[745];
	v[772] = v[4114] * v[620] + v[532] * v[647];
	v[770] = v[531] * v[647] + v[1782] * v[741];
	v[769] = v[4114] * v[621] + v[532] * v[648];
	v[767] = v[531] * v[648] + v[1782] * v[737];
	v[766] = v[4114] * v[622] + v[532] * v[649];
	v[764] = v[531] * v[649] + v[1782] * v[733];
	v[763] = v[4114] * v[623] + v[532] * v[650];
	v[761] = v[531] * v[650] + v[1782] * v[729];
	v[760] = v[4114] * v[624] + v[532] * v[651];
	v[758] = v[531] * v[651] + v[1782] * v[725];
	v[757] = v[4114] * v[625] + v[532] * v[652];
	v[755] = v[531] * v[652] + v[1782] * v[721];
	v[664] = v[1782] * v[652];
	v[663] = v[1782] * v[651];
	v[662] = v[1782] * v[650];
	v[658] = v[1782] * v[649];
	v[657] = v[1782] * v[648];
	v[656] = v[1782] * v[647];
	v[655] = v[1782] * v[646];
	v[4133] = 2e0*v[655];
	v[654] = v[1782] * v[645];
	v[4135] = 2e0*v[654];
	v[653] = v[1782] * v[644];
	v[4137] = 2e0*v[653];
	v[415] = v[1782] * v[396] + v[380] * v[532];
	v[4125] = -2e0*v[415];
	v[414] = v[1782] * v[397] + v[382] * v[532];
	v[4124] = 2e0*v[414];
	v[413] = v[1782] * v[398] + v[384] * v[532];
	v[4129] = 2e0*v[413];
	v[412] = v[1782] * v[386] + v[380] * v[531];
	v[4127] = -2e0*v[412];
	v[411] = v[1782] * v[388] + v[382] * v[531];
	v[4126] = 2e0*v[411];
	v[410] = v[1782] * v[390] + v[384] * v[531];
	v[4128] = 2e0*v[410];
	v[353] = v[1782] * v[380];
	v[4121] = -(v[353] * v[644]);
	v[4118] = -(v[353] * v[645]);
	v[4115] = -(v[353] * v[646]);
	v[355] = v[1782] * v[382];
	v[4122] = -(v[355] * v[647]);
	v[4119] = -(v[355] * v[648]);
	v[4116] = -(v[355] * v[649]);
	v[1808] = v[1803] * v[343] - v[1801] * v[345] - v[1773] * v[353] + v[1771] * v[355] + v[406] * v[4126] + v[405] * v[4127];
	v[1807] = v[1795] * v[343] - v[1791] * v[345] - v[1761] * v[353] + v[1757] * v[355] + v[409] * v[4124] + v[408] * v[4125];
	v[1806] = v[1793] * v[343] - v[1789] * v[345] - v[1759] * v[353] + v[1755] * v[355] + v[409] * v[411] - v[408] * v[412]
		+ v[406] * v[414] - v[405] * v[415];
	v[661] = v[355] * v[637] - v[353] * v[640] - v[345] * v[655] + v[343] * v[658];
	v[660] = v[355] * v[636] - v[353] * v[639] - v[345] * v[654] + v[343] * v[657];
	v[659] = v[355] * v[635] - v[353] * v[638] - v[345] * v[653] + v[343] * v[656];
	v[534] = -(v[353] * v[408]) + v[355] * v[409] + v[343] * v[414] - v[345] * v[415];
	v[533] = -(v[353] * v[405]) + v[355] * v[406] + v[343] * v[411] - v[345] * v[412];
	v[362] = -(v[345] * v[353]) + v[343] * v[355];
	v[356] = v[1782] * v[384];
	v[4123] = -(v[356] * v[650]);
	v[4120] = -(v[356] * v[651]);
	v[4117] = -(v[356] * v[652]);
	v[2260] = 2e0*v[4115] + 2e0*v[4116] + 2e0*v[4117] - v[415] * v[676] - v[414] * v[682] - v[413] * v[688] + v[370] * v[757]
		+ v[369] * v[766] + v[368] * v[775];
	v[2259] = 2e0*v[4118] + 2e0*v[4119] + 2e0*v[4120] - v[415] * v[674] - v[414] * v[680] - v[413] * v[686] + v[370] * v[760]
		+ v[369] * v[769] + v[368] * v[778];
	v[2258] = 2e0*v[4121] + 2e0*v[4122] + 2e0*v[4123] - v[415] * v[672] - v[414] * v[678] - v[413] * v[684] + v[370] * v[763]
		+ v[369] * v[772] + v[368] * v[781];
	v[2254] = -(v[353] * v[628]) - v[355] * v[631] - v[356] * v[634] - v[379] * v[655] - v[381] * v[658] - v[383] * v[664]
		- v[412] * v[676] - v[411] * v[682] - v[410] * v[688] + v[370] * v[755] + v[369] * v[764] + v[368] * v[773];
	v[2253] = -(v[353] * v[627]) - v[355] * v[630] - v[356] * v[633] - v[379] * v[654] - v[381] * v[657] - v[383] * v[663]
		- v[412] * v[674] - v[411] * v[680] - v[410] * v[686] + v[370] * v[758] + v[369] * v[767] + v[368] * v[776];
	v[2252] = -(v[353] * v[626]) - v[355] * v[629] - v[356] * v[632] - v[379] * v[653] - v[381] * v[656] - v[383] * v[662]
		- v[412] * v[672] - v[411] * v[678] - v[410] * v[684] + v[370] * v[761] + v[369] * v[770] + v[368] * v[779];
	v[1865] = v[1791] * v[368] + v[1795] * v[369] + v[1799] * v[370] - v[353] * v[396] - v[355] * v[397] - v[356] * v[398]
		- v[382] * v[4124] + v[380] * v[4125] - v[4109] * v[413];
	v[1863] = v[1789] * v[368] + v[1793] * v[369] + v[1797] * v[370] - v[353] * v[386] - v[355] * v[388] - v[356] * v[390]
		- v[384] * v[410] - v[382] * v[411] - v[380] * v[412] - v[383] * v[413] - v[381] * v[414] - v[379] * v[415];
	v[1862] = v[1801] * v[368] + v[1803] * v[369] + v[1805] * v[370] - v[353] * v[385] - v[355] * v[387] - v[356] * v[389]
		- v[410] * v[4104] - v[381] * v[4126] + v[379] * v[4127];
	v[1814] = -(v[1805] * v[343]) + v[1801] * v[346] + v[1775] * v[353] - v[1771] * v[356] - v[404] * v[4127]
		- v[406] * v[4128];
	v[1813] = -(v[1799] * v[343]) + v[1791] * v[346] + v[1765] * v[353] - v[1757] * v[356] - v[407] * v[4125]
		- v[409] * v[4129];
	v[1812] = -(v[1797] * v[343]) + v[1789] * v[346] + v[1763] * v[353] - v[1755] * v[356] - v[409] * v[410] + v[407] * v[412]
		- v[406] * v[413] + v[404] * v[415];
	v[1811] = v[1805] * v[345] - v[1803] * v[346] - v[1775] * v[355] + v[1773] * v[356] - v[404] * v[4126] + v[405] * v[4128];
	v[1810] = v[1799] * v[345] - v[1795] * v[346] - v[1765] * v[355] + v[1761] * v[356] - v[407] * v[4124] + v[408] * v[4129];
	v[1809] = v[1797] * v[345] - v[1793] * v[346] - v[1763] * v[355] + v[1759] * v[356] + v[408] * v[410] - v[407] * v[411]
		+ v[405] * v[413] - v[404] * v[414];
	v[670] = v[356] * v[640] - v[355] * v[643] - v[346] * v[658] + v[345] * v[664];
	v[669] = v[356] * v[639] - v[355] * v[642] - v[346] * v[657] + v[345] * v[663];
	v[668] = v[356] * v[638] - v[355] * v[641] - v[346] * v[656] + v[345] * v[662];
	v[667] = -(v[356] * v[637]) + v[353] * v[643] + v[346] * v[655] - v[343] * v[664];
	v[666] = -(v[356] * v[636]) + v[353] * v[642] + v[346] * v[654] - v[343] * v[663];
	v[665] = -(v[356] * v[635]) + v[353] * v[641] + v[346] * v[653] - v[343] * v[662];
	v[538] = -(v[355] * v[407]) + v[356] * v[408] + v[345] * v[413] - v[346] * v[414];
	v[537] = -(v[355] * v[404]) + v[356] * v[405] + v[345] * v[410] - v[346] * v[411];
	v[536] = v[353] * v[407] - v[356] * v[409] - v[343] * v[413] + v[346] * v[415];
	v[535] = v[353] * v[404] - v[356] * v[406] - v[343] * v[410] + v[346] * v[412];
	v[360] = v[346] * v[353] - v[343] * v[356];
	v[357] = -(v[346] * v[355]) + v[345] * v[356];
	v[783] = 2e0*(v[362] * v[534] + v[360] * v[536] + v[357] * v[538]);
	v[782] = 2e0*(v[362] * v[533] + v[360] * v[535] + v[357] * v[537]);
	v[540] = (v[357] * v[357]) + (v[360] * v[360]) + (v[362] * v[362]);
	v[1818] = 1e0 / sqrt(v[540]);
	v[1821] = (0.75e0*v[1818]) / Power(v[540], 2);
	v[1819] = -v[1818] / (2e0*v[540]);
	v[4130] = 2e0*v[1819];
	v[1823] = v[4130] * (v[1811] * v[357] + v[1814] * v[360] + v[1808] * v[362] + (v[533] * v[533]) + (v[535] * v[535]) +
		(v[537] * v[537])) + v[1821] * (v[782] * v[782]);
	v[1822] = v[4130] * (v[1810] * v[357] + v[1813] * v[360] + v[1807] * v[362] + (v[534] * v[534]) + (v[536] * v[536]) +
		(v[538] * v[538])) + v[1821] * (v[783] * v[783]);
	v[1820] = v[4130] * (v[1809] * v[357] + v[1812] * v[360] + v[1806] * v[362] + v[533] * v[534] + v[535] * v[536]
		+ v[537] * v[538]) + v[1821] * v[782] * v[783];
	v[785] = v[1819] * v[783];
	v[4132] = 2e0*v[785];
	v[784] = v[1819] * v[782];
	v[4131] = 2e0*v[784];
	v[1842] = v[1811] * v[1818] + v[1823] * v[357] + v[4131] * v[537];
	v[3820] = -((*radius)*v[1842]) - v[385];
	v[1836] = v[1810] * v[1818] + v[1822] * v[357] + v[4132] * v[538];
	v[3819] = -((*radius)*v[1836]) - v[396];
	v[1834] = v[1809] * v[1818] + v[1820] * v[357] + v[538] * v[784] + v[537] * v[785];
	v[3818] = -((*radius)*v[1834]) - v[386];
	v[1833] = v[1814] * v[1818] + v[1823] * v[360] + v[4131] * v[535];
	v[3817] = -((*radius)*v[1833]) - v[387];
	v[1831] = v[1813] * v[1818] + v[1822] * v[360] + v[4132] * v[536];
	v[3816] = -((*radius)*v[1831]) - v[397];
	v[1829] = v[1812] * v[1818] + v[1820] * v[360] + v[536] * v[784] + v[535] * v[785];
	v[3815] = -((*radius)*v[1829]) - v[388];
	v[1828] = v[1808] * v[1818] + v[1823] * v[362] + v[4131] * v[533];
	v[3814] = -((*radius)*v[1828]) - v[389];
	v[1844] = v[1842] * v[316] + v[1833] * v[318] + v[1828] * v[319];
	v[1843] = v[1842] * v[331] + v[1833] * v[334] + v[1828] * v[336];
	v[1826] = v[1807] * v[1818] + v[1822] * v[362] + v[4132] * v[534];
	v[3813] = -((*radius)*v[1826]) - v[398];
	v[1840] = v[1836] * v[316] + v[1831] * v[318] + v[1826] * v[319];
	v[1838] = v[1836] * v[331] + v[1831] * v[334] + v[1826] * v[336];
	v[1824] = v[1806] * v[1818] + v[1820] * v[362] + v[534] * v[784] + v[533] * v[785];
	v[3812] = -((*radius)*v[1824]) - v[390];
	v[1839] = v[1834] * v[316] + v[1829] * v[318] + v[1824] * v[319];
	v[1837] = v[1834] * v[331] + v[1829] * v[334] + v[1824] * v[336];
	v[1103] = v[1818] * (v[1022] * v[343] - v[1034] * v[345] - v[4133] * v[640] + v[4134] * v[658] - v[353] * v[986]
		+ v[355] * v[998]);
	v[3592] = -((*radius)*v[1103]) - cp[1] * v[924] - cp[0] * v[925] - v[304] * v[953];
	v[1101] = v[1818] * (-(v[1010] * v[343]) + v[1034] * v[346] + v[4133] * v[643] - v[4134] * v[664] + v[353] * v[974]
		- v[356] * v[998]);
	v[3599] = -((*radius)*v[1101]) - cp[1] * v[928] - cp[0] * v[930] - v[304] * v[959];
	v[1099] = v[1818] * (v[1010] * v[345] - v[1022] * v[346] - v[4146] * v[658] + v[4142] * v[664] - v[355] * v[974]
		+ v[356] * v[986]);
	v[3606] = -((*radius)*v[1099]) - cp[1] * v[945] - cp[0] * v[947] - v[304] * v[965];
	v[3131] = v[1099] * v[316] + v[1101] * v[318] + v[1103] * v[319];
	v[3128] = v[1099] * v[331] + v[1101] * v[334] + v[1103] * v[336];
	v[1097] = v[1818] * (v[1020] * v[343] - v[1032] * v[345] - v[640] * v[654] - v[639] * v[655] + v[637] * v[657]
		+ v[636] * v[658] - v[353] * v[984] + v[355] * v[996]);
	v[3569] = -((*radius)*v[1097]) - cp[1] * v[889] - cp[0] * v[892] - v[304] * v[952];
	v[1095] = v[1818] * (v[1026] * v[343] - v[1038] * v[345] + v[1002] * v[355] - v[4135] * v[639] + v[4136] * v[657]
		- v[353] * v[990]);
	v[3567] = -((*radius)*v[1095]) - cp[1] * v[888] - cp[0] * v[891] - v[304] * v[955];
	v[1093] = v[1818] * (-(v[1008] * v[343]) + v[1032] * v[346] + v[643] * v[654] + v[642] * v[655] - v[637] * v[663]
		- v[636] * v[664] + v[353] * v[972] - v[356] * v[996]);
	v[3577] = -((*radius)*v[1093]) - cp[1] * v[915] - cp[0] * v[918] - v[304] * v[958];
	v[1091] = v[1818] * (-(v[1014] * v[343]) + v[1038] * v[346] - v[1002] * v[356] + v[4135] * v[642] - v[4136] * v[663]
		+ v[353] * v[978]);
	v[3575] = -((*radius)*v[1091]) - cp[1] * v[914] - cp[0] * v[917] - v[304] * v[961];
	v[1089] = v[1818] * (v[1008] * v[345] - v[1020] * v[346] - v[643] * v[657] - v[642] * v[658] + v[640] * v[663]
		+ v[639] * v[664] - v[355] * v[972] + v[356] * v[984]);
	v[3585] = -((*radius)*v[1089]) - cp[1] * v[944] - cp[0] * v[946] - v[304] * v[964];
	v[3130] = v[1089] * v[316] + v[1093] * v[318] + v[1097] * v[319];
	v[3127] = v[1089] * v[331] + v[1093] * v[334] + v[1097] * v[336];
	v[1087] = v[1818] * (v[1014] * v[345] - v[1026] * v[346] - v[4145] * v[657] + v[4141] * v[663] - v[355] * v[978]
		+ v[356] * v[990]);
	v[3583] = -((*radius)*v[1087]) - cp[1] * v[949] - cp[0] * v[950] - v[304] * v[967];
	v[3135] = v[1087] * v[316] + v[1091] * v[318] + v[1095] * v[319];
	v[3133] = v[1087] * v[331] + v[1091] * v[334] + v[1095] * v[336];
	v[1085] = v[1818] * (v[1018] * v[343] - v[1030] * v[345] - v[640] * v[653] - v[638] * v[655] + v[637] * v[656]
		+ v[635] * v[658] - v[353] * v[982] + v[355] * v[994]);
	v[3543] = -((*radius)*v[1085]) - cp[1] * v[896] - cp[0] * v[898] - v[304] * v[951];
	v[1083] = v[1818] * (v[1024] * v[343] - v[1036] * v[345] + v[1000] * v[355] - v[639] * v[653] - v[638] * v[654]
		+ v[636] * v[656] + v[635] * v[657] - v[353] * v[988]);
	v[3541] = -((*radius)*v[1083]) - cp[1] * v[887] - cp[0] * v[890] - v[304] * v[954];
	v[1081] = v[1818] * (v[1028] * v[343] - v[1040] * v[345] + v[1004] * v[355] - v[4137] * v[638] + v[4138] * v[656]
		- v[353] * v[992]);
	v[3539] = -((*radius)*v[1081]) - cp[1] * v[895] - cp[0] * v[897] - v[304] * v[956];
	v[1079] = v[1818] * (-(v[1006] * v[343]) + v[1030] * v[346] + v[643] * v[653] + v[641] * v[655] - v[637] * v[662]
		- v[635] * v[664] + v[353] * v[970] - v[356] * v[994]);
	v[3552] = -((*radius)*v[1079]) - cp[1] * v[927] - cp[0] * v[929] - v[304] * v[957];
	v[1077] = v[1818] * (-(v[1012] * v[343]) + v[1036] * v[346] - v[1000] * v[356] + v[642] * v[653] + v[641] * v[654]
		- v[636] * v[662] - v[635] * v[663] + v[353] * v[976]);
	v[3550] = -((*radius)*v[1077]) - cp[1] * v[913] - cp[0] * v[916] - v[304] * v[960];
	v[1075] = v[1818] * (-(v[1016] * v[343]) + v[1040] * v[346] - v[1004] * v[356] + v[4137] * v[641] - v[4138] * v[662]
		+ v[353] * v[980]);
	v[3548] = -((*radius)*v[1075]) - cp[1] * v[932] - cp[0] * v[933] - v[304] * v[962];
	v[1073] = v[1818] * (v[1006] * v[345] - v[1018] * v[346] - v[643] * v[656] - v[641] * v[658] + v[640] * v[662]
		+ v[638] * v[664] - v[355] * v[970] + v[356] * v[982]);
	v[3561] = -((*radius)*v[1073]) - cp[1] * v[937] - cp[0] * v[940] - v[304] * v[963];
	v[3129] = v[1073] * v[316] + v[1079] * v[318] + v[1085] * v[319];
	v[3126] = v[1073] * v[331] + v[1079] * v[334] + v[1085] * v[336];
	v[1071] = v[1818] * (v[1012] * v[345] - v[1024] * v[346] - v[642] * v[656] - v[641] * v[657] + v[639] * v[662]
		+ v[638] * v[663] - v[355] * v[976] + v[356] * v[988]);
	v[3559] = -((*radius)*v[1071]) - cp[1] * v[936] - cp[0] * v[939] - v[304] * v[966];
	v[3134] = v[1071] * v[316] + v[1077] * v[318] + v[1083] * v[319];
	v[3132] = v[1071] * v[331] + v[1077] * v[334] + v[1083] * v[336];
	v[1069] = v[1818] * (v[1016] * v[345] - v[1028] * v[346] - v[4144] * v[656] + v[4140] * v[662] - v[355] * v[980]
		+ v[356] * v[992]);
	v[3557] = -((*radius)*v[1069]) - cp[1] * v[935] - cp[0] * v[938] - v[304] * v[968];
	v[3137] = v[1069] * v[316] + v[1075] * v[318] + v[1081] * v[319];
	v[3136] = v[1069] * v[331] + v[1075] * v[334] + v[1081] * v[336];
	v[827] = v[1818] * (v[414] * v[637] - v[415] * v[640] - v[408] * v[655] + v[409] * v[658] - v[353] * v[734] + v[355] * v[746]
		+ v[343] * v[766] - v[345] * v[775]) + v[661] * v[785];
	v[3430] = -v[652] - (*radius)*v[827];
	v[826] = v[1818] * (v[411] * v[637] - v[412] * v[640] - v[405] * v[655] + v[406] * v[658] - v[353] * v[732] + v[355] * v[744]
		+ v[343] * v[764] - v[345] * v[773]) + v[661] * v[784];
	v[3433] = -v[634] - (*radius)*v[826];
	v[825] = v[1818] * (-(v[413] * v[637]) + v[415] * v[643] + v[407] * v[655] - v[409] * v[664] + v[353] * v[722]
		- v[356] * v[746] - v[343] * v[757] + v[346] * v[775]) + v[667] * v[785];
	v[3436] = -v[649] - (*radius)*v[825];
	v[824] = v[1818] * (-(v[410] * v[637]) + v[412] * v[643] + v[404] * v[655] - v[406] * v[664] + v[353] * v[720]
		- v[356] * v[744] - v[343] * v[755] + v[346] * v[773]) + v[667] * v[784];
	v[3439] = -v[631] - (*radius)*v[824];
	v[823] = v[1818] * (v[413] * v[640] - v[414] * v[643] - v[407] * v[658] + v[408] * v[664] - v[355] * v[722] + v[356] * v[734]
		+ v[345] * v[757] - v[346] * v[766]) + v[670] * v[785];
	v[3442] = -v[646] - (*radius)*v[823];
	v[2212] = v[316] * v[823] + v[318] * v[825] + v[319] * v[827];
	v[2209] = v[331] * v[823] + v[334] * v[825] + v[336] * v[827];
	v[822] = v[1818] * (v[410] * v[640] - v[411] * v[643] - v[404] * v[658] + v[405] * v[664] - v[355] * v[720] + v[356] * v[732]
		+ v[345] * v[755] - v[346] * v[764]) + v[670] * v[784];
	v[3445] = -v[628] - (*radius)*v[822];
	v[2218] = v[316] * v[822] + v[318] * v[824] + v[319] * v[826];
	v[2215] = v[331] * v[822] + v[334] * v[824] + v[336] * v[826];
	v[821] = v[1818] * (v[414] * v[636] - v[415] * v[639] - v[408] * v[654] + v[409] * v[657] - v[353] * v[738] + v[355] * v[750]
		+ v[343] * v[769] - v[345] * v[778]) + v[660] * v[785];
	v[3429] = -v[651] - (*radius)*v[821];
	v[820] = v[1818] * (v[411] * v[636] - v[412] * v[639] - v[405] * v[654] + v[406] * v[657] - v[353] * v[736] + v[355] * v[748]
		+ v[343] * v[767] - v[345] * v[776]) + v[660] * v[784];
	v[3432] = -v[633] - (*radius)*v[820];
	v[819] = v[1818] * (-(v[413] * v[636]) + v[415] * v[642] + v[407] * v[654] - v[409] * v[663] + v[353] * v[726]
		- v[356] * v[750] - v[343] * v[760] + v[346] * v[778]) + v[666] * v[785];
	v[3435] = -v[648] - (*radius)*v[819];
	v[818] = v[1818] * (-(v[410] * v[636]) + v[412] * v[642] + v[404] * v[654] - v[406] * v[663] + v[353] * v[724]
		- v[356] * v[748] - v[343] * v[758] + v[346] * v[776]) + v[666] * v[784];
	v[3438] = -v[630] - (*radius)*v[818];
	v[817] = v[1818] * (v[413] * v[639] - v[414] * v[642] - v[407] * v[657] + v[408] * v[663] - v[355] * v[726] + v[356] * v[738]
		+ v[345] * v[760] - v[346] * v[769]) + v[669] * v[785];
	v[3441] = -v[645] - (*radius)*v[817];
	v[2211] = v[316] * v[817] + v[318] * v[819] + v[319] * v[821];
	v[2208] = v[331] * v[817] + v[334] * v[819] + v[336] * v[821];
	v[816] = v[1818] * (v[410] * v[639] - v[411] * v[642] - v[404] * v[657] + v[405] * v[663] - v[355] * v[724] + v[356] * v[736]
		+ v[345] * v[758] - v[346] * v[767]) + v[669] * v[784];
	v[3444] = -v[627] - (*radius)*v[816];
	v[2217] = v[316] * v[816] + v[318] * v[818] + v[319] * v[820];
	v[2214] = v[331] * v[816] + v[334] * v[818] + v[336] * v[820];
	v[815] = v[1818] * (v[414] * v[635] - v[415] * v[638] - v[408] * v[653] + v[409] * v[656] - v[353] * v[742] + v[355] * v[754]
		+ v[343] * v[772] - v[345] * v[781]) + v[659] * v[785];
	v[3428] = -v[650] - (*radius)*v[815];
	v[814] = v[1818] * (v[411] * v[635] - v[412] * v[638] - v[405] * v[653] + v[406] * v[656] - v[353] * v[740] + v[355] * v[752]
		+ v[343] * v[770] - v[345] * v[779]) + v[659] * v[784];
	v[3431] = -v[632] - (*radius)*v[814];
	v[813] = v[1818] * (-(v[413] * v[635]) + v[415] * v[641] + v[407] * v[653] - v[409] * v[662] + v[353] * v[730]
		- v[356] * v[754] - v[343] * v[763] + v[346] * v[781]) + v[665] * v[785];
	v[3434] = -v[647] - (*radius)*v[813];
	v[812] = v[1818] * (-(v[410] * v[635]) + v[412] * v[641] + v[404] * v[653] - v[406] * v[662] + v[353] * v[728]
		- v[356] * v[752] - v[343] * v[761] + v[346] * v[779]) + v[665] * v[784];
	v[3437] = -v[629] - (*radius)*v[812];
	v[811] = v[1818] * (v[413] * v[638] - v[414] * v[641] - v[407] * v[656] + v[408] * v[662] - v[355] * v[730] + v[356] * v[742]
		+ v[345] * v[763] - v[346] * v[772]) + v[668] * v[785];
	v[3440] = -v[644] - (*radius)*v[811];
	v[2210] = v[316] * v[811] + v[318] * v[813] + v[319] * v[815];
	v[2207] = v[331] * v[811] + v[334] * v[813] + v[336] * v[815];
	v[810] = v[1818] * (v[410] * v[638] - v[411] * v[641] - v[404] * v[656] + v[405] * v[662] - v[355] * v[728] + v[356] * v[740]
		+ v[345] * v[761] - v[346] * v[770]) + v[668] * v[784];
	v[3443] = -v[626] - (*radius)*v[810];
	v[2216] = v[316] * v[810] + v[318] * v[812] + v[319] * v[814];
	v[2213] = v[331] * v[810] + v[334] * v[812] + v[336] * v[814];
	v[687] = v[1818] * v[661];
	v[4346] = -2e0*v[687];
	v[1067] = -((*radius)*v[687]) - v[688];
	v[685] = v[1818] * v[660];
	v[4343] = -2e0*v[685];
	v[1066] = -((*radius)*v[685]) - v[686];
	v[683] = v[1818] * v[659];
	v[4342] = -2e0*v[683];
	v[1065] = -((*radius)*v[683]) - v[684];
	v[681] = v[1818] * v[667];
	v[4341] = -2e0*v[681];
	v[1064] = -((*radius)*v[681]) - v[682];
	v[679] = v[1818] * v[666];
	v[4340] = -2e0*v[679];
	v[1063] = -((*radius)*v[679]) - v[680];
	v[677] = v[1818] * v[665];
	v[4337] = -2e0*v[677];
	v[1062] = -((*radius)*v[677]) - v[678];
	v[675] = v[1818] * v[670];
	v[4336] = -2e0*v[675];
	v[3176] = v[1296] * v[675] + v[1315] * v[681] + v[1335] * v[687];
	v[3172] = v[1301] * v[675] + v[1319] * v[681] + v[1340] * v[687];
	v[3167] = v[1305] * v[675] + v[1324] * v[681] + v[1344] * v[687];
	v[1287] = v[331] * v[675] + v[334] * v[681] + v[336] * v[687];
	v[1269] = v[316] * v[675] + v[318] * v[681] + v[319] * v[687];
	v[1061] = -((*radius)*v[675]) - v[676];
	v[673] = v[1818] * v[669];
	v[4335] = -2e0*v[673];
	v[3175] = v[1296] * v[673] + v[1315] * v[679] + v[1335] * v[685];
	v[3171] = v[1301] * v[673] + v[1319] * v[679] + v[1340] * v[685];
	v[3166] = v[1305] * v[673] + v[1324] * v[679] + v[1344] * v[685];
	v[1286] = v[331] * v[673] + v[334] * v[679] + v[336] * v[685];
	v[1268] = v[316] * v[673] + v[318] * v[679] + v[319] * v[685];
	v[1060] = -((*radius)*v[673]) - v[674];
	v[671] = v[1818] * v[668];
	v[4334] = -2e0*v[671];
	v[3174] = v[1296] * v[671] + v[1315] * v[677] + v[1335] * v[683];
	v[3170] = v[1301] * v[671] + v[1319] * v[677] + v[1340] * v[683];
	v[3165] = v[1305] * v[671] + v[1324] * v[677] + v[1344] * v[683];
	v[1285] = v[331] * v[671] + v[334] * v[677] + v[336] * v[683];
	v[1267] = v[316] * v[671] + v[318] * v[677] + v[319] * v[683];
	v[1059] = -((*radius)*v[671]) - v[672];
	v[549] = v[1818] * v[534] + v[362] * v[785];
	v[4294] = 2e0*v[549];
	v[809] = -v[384] - (*radius)*v[549];
	v[548] = v[1818] * v[533] + v[362] * v[784];
	v[4293] = 2e0*v[548];
	v[808] = -v[383] - (*radius)*v[548];
	v[547] = v[1818] * v[536] + v[360] * v[785];
	v[4150] = 2e0*v[547];
	v[807] = -v[382] - (*radius)*v[547];
	v[546] = v[1818] * v[535] + v[360] * v[784];
	v[4148] = 2e0*v[546];
	v[806] = -v[381] - (*radius)*v[546];
	v[545] = v[1818] * v[538] + v[357] * v[785];
	v[4292] = -2e0*v[545];
	v[2929] = v[1305] * v[545] + v[1324] * v[547] + v[1344] * v[549];
	v[2928] = v[1301] * v[545] + v[1319] * v[547] + v[1340] * v[549];
	v[2927] = v[1296] * v[545] + v[1315] * v[547] + v[1335] * v[549];
	v[1214] = v[331] * v[545] + v[334] * v[547] + v[336] * v[549];
	v[1202] = v[316] * v[545] + v[318] * v[547] + v[319] * v[549];
	v[805] = -v[380] - (*radius)*v[545];
	v[544] = v[1818] * v[537] + v[357] * v[784];
	v[4291] = -2e0*v[544];
	v[2935] = v[1305] * v[544] + v[1324] * v[546] + v[1344] * v[548];
	v[2934] = v[1301] * v[544] + v[1319] * v[546] + v[1340] * v[548];
	v[2933] = v[1296] * v[544] + v[1315] * v[546] + v[1335] * v[548];
	v[1213] = v[331] * v[544] + v[334] * v[546] + v[336] * v[548];
	v[1201] = v[316] * v[544] + v[318] * v[546] + v[319] * v[548];
	v[804] = -v[379] - (*radius)*v[544];
	v[358] = v[1818] * v[357];
	v[4311] = v[358] * v[4297];
	v[4139] = -2e0*v[358];
	v[2980] = v[4139] * v[544];
	v[2973] = v[4139] * v[545];
	v[361] = v[1818] * v[360];
	v[4309] = v[361] * v[4297];
	v[4143] = -2e0*v[361];
	v[3152] = -(v[1075] * v[343]) + v[1069] * v[345] - v[1004] * v[361] + v[4140] * v[671] - v[4138] * v[677]
		+ v[358] * v[992];
	v[3151] = -(v[1091] * v[343]) + v[1087] * v[345] - v[1002] * v[361] + v[4141] * v[673] - v[4136] * v[679]
		+ v[358] * v[990];
	v[3150] = -(v[1077] * v[343]) + v[1071] * v[345] - v[1000] * v[361] + v[639] * v[671] + v[638] * v[673] - v[636] * v[677]
		- v[635] * v[679] + v[358] * v[988];
	v[3149] = -(v[1101] * v[343]) + v[1099] * v[345] + v[4142] * v[675] - v[4134] * v[681] + v[358] * v[986] - v[361] * v[998];
	v[3148] = -(v[1093] * v[343]) + v[1089] * v[345] + v[640] * v[673] + v[639] * v[675] - v[637] * v[679] - v[636] * v[681]
		+ v[358] * v[984] - v[361] * v[996];
	v[3147] = -(v[1079] * v[343]) + v[1073] * v[345] + v[640] * v[671] + v[638] * v[675] - v[637] * v[677] - v[635] * v[681]
		+ v[358] * v[982] - v[361] * v[994];
	v[3146] = -(v[361] * v[675]) - v[358] * v[681];
	v[3145] = -(v[361] * v[673]) - v[358] * v[679];
	v[3144] = -(v[361] * v[671]) - v[358] * v[677];
	v[2966] = v[4143] * v[546];
	v[2965] = -(v[361] * v[544]) - v[358] * v[546];
	v[2958] = v[4143] * v[547];
	v[2957] = -(v[361] * v[545]) - v[358] * v[547];
	v[2224] = -(v[546] * v[637]) + v[544] * v[640] + v[405] * v[675] - v[406] * v[681] + v[358] * v[732] - v[361] * v[744]
		+ v[345] * v[822] - v[343] * v[824];
	v[2223] = -(v[546] * v[636]) + v[544] * v[639] + v[405] * v[673] - v[406] * v[679] + v[358] * v[736] - v[361] * v[748]
		+ v[345] * v[816] - v[343] * v[818];
	v[2222] = -(v[546] * v[635]) + v[544] * v[638] + v[405] * v[671] - v[406] * v[677] + v[358] * v[740] - v[361] * v[752]
		+ v[345] * v[810] - v[343] * v[812];
	v[2221] = -(v[547] * v[637]) + v[545] * v[640] + v[408] * v[675] - v[409] * v[681] + v[358] * v[734] - v[361] * v[746]
		+ v[345] * v[823] - v[343] * v[825];
	v[2220] = -(v[547] * v[636]) + v[545] * v[639] + v[408] * v[673] - v[409] * v[679] + v[358] * v[738] - v[361] * v[750]
		+ v[345] * v[817] - v[343] * v[819];
	v[2219] = -(v[547] * v[635]) + v[545] * v[638] + v[408] * v[671] - v[409] * v[677] + v[358] * v[742] - v[361] * v[754]
		+ v[345] * v[811] - v[343] * v[813];
	v[1847] = -(v[1833] * v[343]) + v[1842] * v[345] + v[1773] * v[358] - v[1771] * v[361] - v[406] * v[4148]
		+ v[4149] * v[544];
	v[1846] = -(v[1831] * v[343]) + v[1836] * v[345] + v[1761] * v[358] - v[1757] * v[361] - v[409] * v[4150]
		+ v[4151] * v[545];
	v[1845] = -(v[1829] * v[343]) + v[1834] * v[345] + v[1759] * v[358] - v[1755] * v[361] + v[408] * v[544] + v[405] * v[545]
		- v[409] * v[546] - v[406] * v[547];
	v[1446] = -(v[358] * v[361]);
	v[1260] = -(v[361] * v[637]) + v[358] * v[640] + v[345] * v[675] - v[343] * v[681];
	v[1259] = -(v[361] * v[636]) + v[358] * v[639] + v[345] * v[673] - v[343] * v[679];
	v[1258] = -(v[361] * v[635]) + v[358] * v[638] + v[345] * v[671] - v[343] * v[677];
	v[1196] = v[358] * v[408] - v[361] * v[409] + v[345] * v[545] - v[343] * v[547];
	v[1195] = v[358] * v[405] - v[361] * v[406] + v[345] * v[544] - v[343] * v[546];
	v[363] = v[1818] * v[362];
	v[4307] = v[363] * v[4297];
	v[4147] = -2e0*v[363];
	v[3198] = -(v[1081] * v[345]) + v[1075] * v[346] + v[4144] * v[677] - v[4140] * v[683] + v[361] * v[980] - v[363] * v[992];
	v[3193] = -(v[1095] * v[345]) + v[1091] * v[346] + v[4145] * v[679] - v[4141] * v[685] + v[361] * v[978] - v[363] * v[990];
	v[3192] = -(v[1083] * v[345]) + v[1077] * v[346] + v[642] * v[677] + v[641] * v[679] - v[639] * v[683] - v[638] * v[685]
		+ v[361] * v[976] - v[363] * v[988];
	v[3185] = -(v[1103] * v[345]) + v[1101] * v[346] + v[4146] * v[681] - v[4142] * v[687] + v[361] * v[974] - v[363] * v[986];
	v[3184] = -(v[1097] * v[345]) + v[1093] * v[346] + v[643] * v[679] + v[642] * v[681] - v[640] * v[685] - v[639] * v[687]
		+ v[361] * v[972] - v[363] * v[984];
	v[3183] = -(v[1085] * v[345]) + v[1079] * v[346] + v[643] * v[677] + v[641] * v[681] - v[640] * v[683] - v[638] * v[687]
		+ v[361] * v[970] - v[363] * v[982];
	v[3182] = v[1081] * v[343] - v[1069] * v[346] + v[1004] * v[363] - v[4144] * v[671] + v[4138] * v[683] - v[358] * v[980];
	v[3200] = v[318] * v[3182] + v[3152] * v[319] + v[316] * v[3198];
	v[3199] = v[3198] * v[331] + v[3182] * v[334] + v[3152] * v[336];
	v[3181] = v[1095] * v[343] - v[1087] * v[346] + v[1002] * v[363] - v[4145] * v[673] + v[4136] * v[685] - v[358] * v[978];
	v[3197] = v[318] * v[3181] + v[3151] * v[319] + v[316] * v[3193];
	v[3195] = v[3193] * v[331] + v[3181] * v[334] + v[3151] * v[336];
	v[3180] = v[1083] * v[343] - v[1071] * v[346] + v[1000] * v[363] - v[642] * v[671] - v[641] * v[673] + v[636] * v[683]
		+ v[635] * v[685] - v[358] * v[976];
	v[3196] = v[318] * v[3180] + v[3150] * v[319] + v[316] * v[3192];
	v[3194] = v[3192] * v[331] + v[3180] * v[334] + v[3150] * v[336];
	v[3179] = v[1103] * v[343] - v[1099] * v[346] - v[4146] * v[675] + v[4134] * v[687] - v[358] * v[974] + v[363] * v[998];
	v[3191] = v[3179] * v[318] + v[316] * v[3185] + v[3149] * v[319];
	v[3188] = v[3185] * v[331] + v[3179] * v[334] + v[3149] * v[336];
	v[3178] = v[1097] * v[343] - v[1089] * v[346] - v[643] * v[673] - v[642] * v[675] + v[637] * v[685] + v[636] * v[687]
		- v[358] * v[972] + v[363] * v[996];
	v[3190] = v[3178] * v[318] + v[316] * v[3184] + v[3148] * v[319];
	v[3187] = v[3184] * v[331] + v[3178] * v[334] + v[3148] * v[336];
	v[3177] = v[1085] * v[343] - v[1073] * v[346] - v[643] * v[671] - v[641] * v[675] + v[637] * v[683] + v[635] * v[687]
		- v[358] * v[970] + v[363] * v[994];
	v[3189] = v[3177] * v[318] + v[316] * v[3183] + v[3147] * v[319];
	v[3186] = v[3183] * v[331] + v[3177] * v[334] + v[3147] * v[336];
	v[3173] = v[3113] * v[358] + v[3087] * v[361] + v[3059] * v[363];
	v[3169] = v[3109] * v[358] + v[3082] * v[361] + v[3054] * v[363];
	v[3168] = v[3104] * v[358] + v[3078] * v[361] + v[3049] * v[363];
	v[3164] = v[3100] * v[358] + v[3072] * v[361] + v[3045] * v[363];
	v[3163] = v[3095] * v[358] + v[3067] * v[361] + v[3041] * v[363];
	v[3162] = v[3091] * v[358] + v[3063] * v[361] + v[3037] * v[363];
	v[3161] = -(v[363] * v[675]) - v[358] * v[687];
	v[3160] = -(v[363] * v[673]) - v[358] * v[685];
	v[3159] = -(v[363] * v[671]) - v[358] * v[683];
	v[3158] = -(v[363] * v[681]) - v[361] * v[687];
	v[3157] = -(v[363] * v[679]) - v[361] * v[685];
	v[3156] = -(v[363] * v[677]) - v[361] * v[683];
	v[2950] = v[4147] * v[548];
	v[2949] = -(v[363] * v[546]) - v[361] * v[548];
	v[2948] = -(v[363] * v[544]) - v[358] * v[548];
	v[2941] = v[4147] * v[549];
	v[2940] = -(v[363] * v[547]) - v[361] * v[549];
	v[2939] = -(v[363] * v[545]) - v[358] * v[549];
	v[2242] = -(v[548] * v[640]) + v[546] * v[643] + v[404] * v[681] - v[405] * v[687] + v[361] * v[720] - v[363] * v[732]
		+ v[346] * v[824] - v[345] * v[826];
	v[2241] = -(v[548] * v[639]) + v[546] * v[642] + v[404] * v[679] - v[405] * v[685] + v[361] * v[724] - v[363] * v[736]
		+ v[346] * v[818] - v[345] * v[820];
	v[2240] = -(v[548] * v[638]) + v[546] * v[641] + v[404] * v[677] - v[405] * v[683] + v[361] * v[728] - v[363] * v[740]
		+ v[346] * v[812] - v[345] * v[814];
	v[2233] = -(v[549] * v[640]) + v[547] * v[643] + v[407] * v[681] - v[408] * v[687] + v[361] * v[722] - v[363] * v[734]
		+ v[346] * v[825] - v[345] * v[827];
	v[2232] = -(v[549] * v[639]) + v[547] * v[642] + v[407] * v[679] - v[408] * v[685] + v[361] * v[726] - v[363] * v[738]
		+ v[346] * v[819] - v[345] * v[821];
	v[2231] = -(v[549] * v[638]) + v[547] * v[641] + v[407] * v[677] - v[408] * v[683] + v[361] * v[730] - v[363] * v[742]
		+ v[346] * v[813] - v[345] * v[815];
	v[2230] = v[548] * v[637] - v[544] * v[643] - v[404] * v[675] + v[406] * v[687] - v[358] * v[720] + v[363] * v[744]
		- v[346] * v[822] + v[343] * v[826];
	v[2248] = v[2242] * v[316] + v[2230] * v[318] + v[2224] * v[319];
	v[2245] = v[2242] * v[331] + v[2230] * v[334] + v[2224] * v[336];
	v[2229] = v[548] * v[636] - v[544] * v[642] - v[404] * v[673] + v[406] * v[685] - v[358] * v[724] + v[363] * v[748]
		- v[346] * v[816] + v[343] * v[820];
	v[2247] = v[2241] * v[316] + v[2229] * v[318] + v[2223] * v[319];
	v[2244] = v[2241] * v[331] + v[2229] * v[334] + v[2223] * v[336];
	v[2228] = v[548] * v[635] - v[544] * v[641] - v[404] * v[671] + v[406] * v[683] - v[358] * v[728] + v[363] * v[752]
		- v[346] * v[810] + v[343] * v[814];
	v[2246] = v[2240] * v[316] + v[2228] * v[318] + v[2222] * v[319];
	v[2243] = v[2240] * v[331] + v[2228] * v[334] + v[2222] * v[336];
	v[2227] = v[549] * v[637] - v[545] * v[643] - v[407] * v[675] + v[409] * v[687] - v[358] * v[722] + v[363] * v[746]
		- v[346] * v[823] + v[343] * v[827];
	v[2239] = v[2233] * v[316] + v[2227] * v[318] + v[2221] * v[319];
	v[2236] = v[2233] * v[331] + v[2227] * v[334] + v[2221] * v[336];
	v[2226] = v[549] * v[636] - v[545] * v[642] - v[407] * v[673] + v[409] * v[685] - v[358] * v[726] + v[363] * v[750]
		- v[346] * v[817] + v[343] * v[821];
	v[2238] = v[2232] * v[316] + v[2226] * v[318] + v[2220] * v[319];
	v[2235] = v[2232] * v[331] + v[2226] * v[334] + v[2220] * v[336];
	v[2225] = v[549] * v[635] - v[545] * v[641] - v[407] * v[671] + v[409] * v[683] - v[358] * v[730] + v[363] * v[754]
		- v[346] * v[811] + v[343] * v[815];
	v[2237] = v[2231] * v[316] + v[2225] * v[318] + v[2219] * v[319];
	v[2234] = v[2231] * v[331] + v[2225] * v[334] + v[2219] * v[336];
	v[1857] = -(v[1828] * v[345]) + v[1833] * v[346] + v[1775] * v[361] - v[1773] * v[363] + v[404] * v[4148]
		- v[4149] * v[548];
	v[1852] = -(v[1826] * v[345]) + v[1831] * v[346] + v[1765] * v[361] - v[1761] * v[363] + v[407] * v[4150]
		- v[4151] * v[549];
	v[1851] = -(v[1824] * v[345]) + v[1829] * v[346] + v[1763] * v[361] - v[1759] * v[363] + v[407] * v[546] + v[404] * v[547]
		- v[408] * v[548] - v[405] * v[549];
	v[1850] = v[1828] * v[343] - v[1842] * v[346] - v[1775] * v[358] + v[1771] * v[363] + v[404] * v[4291] + v[406] * v[4293];
	v[1859] = v[1857] * v[316] + v[1850] * v[318] + v[1847] * v[319];
	v[1858] = v[1857] * v[331] + v[1850] * v[334] + v[1847] * v[336];
	v[1849] = v[1826] * v[343] - v[1836] * v[346] - v[1765] * v[358] + v[1757] * v[363] + v[407] * v[4292] + v[409] * v[4294];
	v[1856] = v[1852] * v[316] + v[1849] * v[318] + v[1846] * v[319];
	v[1854] = v[1852] * v[331] + v[1849] * v[334] + v[1846] * v[336];
	v[1848] = v[1824] * v[343] - v[1834] * v[346] - v[1763] * v[358] + v[1755] * v[363] - v[407] * v[544] - v[404] * v[545]
		+ v[409] * v[548] + v[406] * v[549];
	v[1855] = v[1851] * v[316] + v[1848] * v[318] + v[1845] * v[319];
	v[1853] = v[1851] * v[331] + v[1848] * v[334] + v[1845] * v[336];
	v[1455] = -(v[361] * v[363]);
	v[1447] = -(v[358] * v[363]);
	v[1350] = v[1305] * v[358] + v[1324] * v[361] + v[1344] * v[363];
	v[2983] = -(v[2935] * v[358]) - v[1350] * v[544];
	v[2976] = -(v[2929] * v[358]) - v[1350] * v[545];
	v[2969] = -(v[2935] * v[361]) - v[1350] * v[546];
	v[2961] = -(v[2929] * v[361]) - v[1350] * v[547];
	v[2953] = -(v[2935] * v[363]) - v[1350] * v[548];
	v[2944] = -(v[2929] * v[363]) - v[1350] * v[549];
	v[1465] = v[1344] - v[1350] * v[363];
	v[1458] = v[1324] - v[1350] * v[361];
	v[1450] = v[1305] - v[1350] * v[358];
	v[1349] = v[1301] * v[358] + v[1319] * v[361] + v[1340] * v[363];
	v[2982] = -(v[2934] * v[358]) - v[1349] * v[544];
	v[2975] = -(v[2928] * v[358]) - v[1349] * v[545];
	v[2968] = -(v[2934] * v[361]) - v[1349] * v[546];
	v[2960] = -(v[2928] * v[361]) - v[1349] * v[547];
	v[2952] = -(v[2934] * v[363]) - v[1349] * v[548];
	v[2943] = -(v[2928] * v[363]) - v[1349] * v[549];
	v[1464] = v[1340] - v[1349] * v[363];
	v[1457] = v[1319] - v[1349] * v[361];
	v[1449] = v[1301] - v[1349] * v[358];
	v[1348] = v[1296] * v[358] + v[1315] * v[361] + v[1335] * v[363];
	v[2981] = -(v[2933] * v[358]) - v[1348] * v[544];
	v[2974] = -(v[2927] * v[358]) - v[1348] * v[545];
	v[2967] = -(v[2933] * v[361]) - v[1348] * v[546];
	v[2959] = -(v[2927] * v[361]) - v[1348] * v[547];
	v[2951] = -(v[2933] * v[363]) - v[1348] * v[548];
	v[2942] = -(v[2927] * v[363]) - v[1348] * v[549];
	v[1463] = v[1335] - v[1348] * v[363];
	v[1456] = v[1315] - v[1348] * v[361];
	v[1448] = v[1296] - v[1348] * v[358];
	v[1257] = v[363] * v[637] - v[358] * v[643] - v[346] * v[675] + v[343] * v[687];
	v[1256] = v[363] * v[636] - v[358] * v[642] - v[346] * v[673] + v[343] * v[685];
	v[1255] = v[363] * v[635] - v[358] * v[641] - v[346] * v[671] + v[343] * v[683];
	v[1254] = -(v[363] * v[640]) + v[361] * v[643] + v[346] * v[681] - v[345] * v[687];
	v[1284] = v[1254] * v[331] + v[1257] * v[334] + v[1260] * v[336];
	v[1266] = v[1254] * v[316] + v[1257] * v[318] + v[1260] * v[319];
	v[1253] = -(v[363] * v[639]) + v[361] * v[642] + v[346] * v[679] - v[345] * v[685];
	v[1283] = v[1253] * v[331] + v[1256] * v[334] + v[1259] * v[336];
	v[1265] = v[1253] * v[316] + v[1256] * v[318] + v[1259] * v[319];
	v[1252] = -(v[363] * v[638]) + v[361] * v[641] + v[346] * v[677] - v[345] * v[683];
	v[1282] = v[1252] * v[331] + v[1255] * v[334] + v[1258] * v[336];
	v[1264] = v[1252] * v[316] + v[1255] * v[318] + v[1258] * v[319];
	v[1194] = -(v[358] * v[407]) + v[363] * v[409] - v[346] * v[545] + v[343] * v[549];
	v[1193] = -(v[358] * v[404]) + v[363] * v[406] - v[346] * v[544] + v[343] * v[548];
	v[1192] = v[361] * v[407] - v[363] * v[408] + v[346] * v[547] - v[345] * v[549];
	v[1212] = v[1192] * v[331] + v[1194] * v[334] + v[1196] * v[336];
	v[1200] = v[1192] * v[316] + v[1194] * v[318] + v[1196] * v[319];
	v[1191] = v[361] * v[404] - v[363] * v[405] + v[346] * v[546] - v[345] * v[548];
	v[1211] = v[1191] * v[331] + v[1193] * v[334] + v[1195] * v[336];
	v[1199] = v[1191] * v[316] + v[1193] * v[318] + v[1195] * v[319];
	v[364] = -((*radius)*v[358]) + v[368];
	v[365] = -((*radius)*v[361]) + v[369];
	v[366] = -((*radius)*v[363]) + v[370];
	v[1188] = v[4297] * sqrt(Power(v[364], 2) + Power(v[365], 2) + Power(v[366], 2));
	v[416] = -(v[343] * v[379]) - v[345] * v[381] - v[346] * v[383] + v[370] * v[404] + v[369] * v[405] + v[368] * v[406];
	v[4266] = v[353] * v[416];
	v[4263] = v[415] * v[416];
	v[4262] = v[355] * v[416];
	v[4259] = v[414] * v[416];
	v[4256] = v[356] * v[416];
	v[4253] = v[413] * v[416];
	v[417] = -(v[353] * v[379]) - v[355] * v[381] - v[356] * v[383] + v[370] * v[410] + v[369] * v[411] + v[368] * v[412];
	v[4265] = v[343] * v[417];
	v[4264] = -(v[409] * v[417]);
	v[4261] = v[345] * v[417];
	v[4260] = -(v[408] * v[417]);
	v[4255] = v[346] * v[417];
	v[4254] = -(v[407] * v[417]);
	v[418] = -(v[343] * v[380]) - v[345] * v[382] - v[346] * v[384] + v[370] * v[407] + v[369] * v[408] + v[368] * v[409];
	v[4282] = v[412] * v[418];
	v[4280] = v[411] * v[418];
	v[4278] = v[410] * v[418];
	v[419] = -(v[353] * v[380]) - v[355] * v[382] - v[356] * v[384] + v[370] * v[413] + v[369] * v[414] + v[368] * v[415];
	v[4286] = v[356] * v[418] - v[346] * v[419];
	v[4285] = v[355] * v[418] - v[345] * v[419];
	v[4284] = v[353] * v[418] - v[343] * v[419];
	v[4283] = v[406] * v[419];
	v[4281] = v[405] * v[419];
	v[4279] = v[404] * v[419];
	v[2037] = -(v[417] * v[418]) + v[416] * v[419];
	v[4275] = v[2252] / v[2037];
	v[4274] = -(v[2249] / v[2037]);
	v[4273] = v[2253] / v[2037];
	v[4272] = -(v[2250] / v[2037]);
	v[4269] = v[2254] / v[2037];
	v[4268] = -(v[2251] / v[2037]);
	v[4258] = v[1861] / v[2037];
	v[4257] = -(v[1863] / v[2037]);
	v[4252] = -(v[2258] / v[2037]);
	v[4251] = v[2255] / v[2037];
	v[4248] = -(v[2259] / v[2037]);
	v[4247] = v[2256] / v[2037];
	v[4245] = -(v[2260] / v[2037]);
	v[4244] = v[2257] / v[2037];
	v[4212] = -(v[419] / v[2037]);
	v[4211] = v[418] / v[2037];
	v[4170] = v[417] / v[2037];
	v[4169] = -(v[416] / v[2037]);
	v[2039] = 1e0 / Power(v[2037], 2);
	v[2734] = -(v[2039] * (v[2260] * v[416] - v[2257] * v[417] - v[2254] * v[418] + v[2251] * v[419]));
	v[4222] = -(v[2734] * v[419]) + v[4245];
	v[4221] = v[2734] * v[418] + v[4244];
	v[4184] = v[2734] * v[417] + v[4269];
	v[4183] = -(v[2734] * v[416]) + v[4268];
	v[4168] = v[2037] * v[2734];
	v[2732] = -(v[2039] * (v[2259] * v[416] - v[2256] * v[417] - v[2253] * v[418] + v[2250] * v[419]));
	v[4216] = v[2732] * v[418];
	v[4224] = v[4216] + v[4247];
	v[4215] = v[2732] * v[419];
	v[4225] = -v[4215] + v[4248];
	v[4176] = v[2732] * v[416];
	v[4186] = -v[4176] + v[4272];
	v[4175] = v[2732] * v[417];
	v[4187] = v[4175] + v[4273];
	v[2730] = -(v[2039] * (v[2258] * v[416] - v[2255] * v[417] - v[2252] * v[418] + v[2249] * v[419]));
	v[4227] = v[2730] * v[418];
	v[4226] = v[2730] * v[419];
	v[4218] = -v[4226] + v[4252];
	v[4217] = v[4227] + v[4251];
	v[4189] = v[2730] * v[416];
	v[4188] = v[2730] * v[417];
	v[4178] = v[4188] + v[4275];
	v[4177] = -v[4189] + v[4274];
	v[2728] = -(v[2039] * (v[4253] + v[4254] - v[4278] + v[4279]));
	v[4196] = v[2037] * v[2728];
	v[4236] = v[413] + v[419] * v[4196];
	v[4235] = v[407] + v[418] * v[4196];
	v[4203] = v[410] + v[417] * v[4196];
	v[4202] = v[404] + v[416] * v[4196];
	v[2726] = -(v[2039] * (v[4259] + v[4260] - v[4280] + v[4281]));
	v[4197] = v[2037] * v[2726];
	v[4238] = v[414] + v[419] * v[4197];
	v[4237] = v[408] + v[418] * v[4197];
	v[4205] = v[411] + v[417] * v[4197];
	v[4204] = v[405] + v[416] * v[4197];
	v[2724] = -(v[2039] * (v[4263] + v[4264] - v[4282] + v[4283]));
	v[4198] = v[2037] * v[2724];
	v[4240] = v[415] + v[419] * v[4198];
	v[4239] = v[409] + v[418] * v[4198];
	v[4207] = v[412] + v[417] * v[4198];
	v[4206] = v[406] + v[416] * v[4198];
	v[2040] = -(v[2039] * (v[1865] * v[416] - v[1864] * v[417] - v[1863] * v[418] + v[1861] * v[419]));
	v[4267] = v[2037] * v[2040];
	v[4277] = v[418] * v[4267];
	v[4276] = v[419] * v[4267];
	v[4271] = v[1865] + v[4276];
	v[4270] = v[1864] + v[4277];
	v[2038] = -(v[2039] * (v[1863] * v[416] - v[1861] * v[417] - v[1862] * v[418] + v[1860] * v[419]));
	v[4246] = v[2037] * v[2038];
	v[4250] = v[1862] + v[417] * v[4246];
	v[4249] = v[1860] + v[416] * v[4246];
	v[424] = (v[368] * v[380] + v[369] * v[382] + v[370] * v[384])*v[4152];
	v[2293] = v[2269] * v[384] + v[424] * v[652] - v[1782] * v[688];
	v[2338] = v[2293] * v[349] - v[304] * v[664];
	v[2292] = v[2267] * v[384] + v[424] * v[651] - v[1782] * v[686];
	v[2337] = v[2292] * v[349] - v[304] * v[663];
	v[2291] = v[2265] * v[384] + v[424] * v[650] - v[1782] * v[684];
	v[2336] = v[2291] * v[349] - v[304] * v[662];
	v[2289] = v[2269] * v[382] + v[424] * v[649] - v[1782] * v[682];
	v[2330] = v[2289] * v[349] - v[304] * v[658];
	v[2288] = v[2267] * v[382] + v[424] * v[648] - v[1782] * v[680];
	v[2329] = v[2288] * v[349] - v[304] * v[657];
	v[2287] = v[2265] * v[382] + v[424] * v[647] - v[1782] * v[678];
	v[2328] = v[2287] * v[349] - v[304] * v[656];
	v[2284] = v[2269] * v[380] + v[424] * v[646] - v[1782] * v[676];
	v[2320] = v[2284] * v[349] - v[304] * v[655];
	v[2283] = v[2267] * v[380] + v[424] * v[645] - v[1782] * v[674];
	v[2319] = v[2283] * v[349] - v[304] * v[654];
	v[2282] = v[2265] * v[380] + v[424] * v[644] - v[1782] * v[672];
	v[2318] = v[2282] * v[349] - v[304] * v[653];
	v[429] = (v[368] * v[379] + v[369] * v[381] + v[370] * v[383])*v[4153];
	v[2308] = v[2278] * v[383] + v[429] * v[634] - v[1748] * v[688];
	v[2334] = v[2308] * v[339] - v[304] * v[643];
	v[2307] = v[2276] * v[383] + v[429] * v[633] - v[1748] * v[686];
	v[2333] = v[2307] * v[339] - v[304] * v[642];
	v[2306] = v[2274] * v[383] + v[429] * v[632] - v[1748] * v[684];
	v[2332] = v[2306] * v[339] - v[304] * v[641];
	v[2304] = v[2278] * v[381] + v[429] * v[631] - v[1748] * v[682];
	v[2325] = v[2304] * v[339] - v[304] * v[640];
	v[2303] = v[2276] * v[381] + v[429] * v[630] - v[1748] * v[680];
	v[2324] = v[2303] * v[339] - v[304] * v[639];
	v[2302] = v[2274] * v[381] + v[429] * v[629] - v[1748] * v[678];
	v[2323] = v[2302] * v[339] - v[304] * v[638];
	v[2299] = v[2278] * v[379] + v[429] * v[628] - v[1748] * v[676];
	v[2314] = v[2299] * v[339] - v[304] * v[637];
	v[2298] = v[2276] * v[379] + v[429] * v[627] - v[1748] * v[674];
	v[2313] = v[2298] * v[339] - v[304] * v[636];
	v[2297] = v[2274] * v[379] + v[429] * v[626] - v[1748] * v[672];
	v[2312] = v[2297] * v[339] - v[304] * v[635];
	v[423] = v[1782] * v[368] + v[380] * v[424];
	v[425] = v[1782] * v[369] + v[382] * v[424];
	v[426] = v[1782] * v[370] + v[384] * v[424];
	v[428] = v[1748] * v[368] + v[379] * v[429];
	v[430] = v[1748] * v[369] + v[381] * v[429];
	v[431] = v[1748] * v[370] + v[383] * v[429];
	v[432] = -(v[304] * v[343]) + v[339] * v[428];
	v[433] = -(v[304] * v[353]) + v[349] * v[423];
	v[434] = -(v[304] * v[345]) + v[339] * v[430];
	v[2399] = v[2308] + v[2325] * v[284] - v[2314] * v[285] + v[434] * v[610] - v[432] * v[613] - cp[0] * v[643];
	v[2398] = v[2307] + v[2324] * v[284] - v[2313] * v[285] + v[434] * v[609] - v[432] * v[612] - cp[0] * v[642];
	v[2397] = v[2306] + v[2323] * v[284] - v[2312] * v[285] + v[434] * v[608] - v[432] * v[611] - cp[0] * v[641];
	v[2369] = -(v[2325] * v[281]) + v[2314] * v[282] - v[434] * v[601] + v[432] * v[604] - cp[1] * v[643];
	v[2464] = e2i[0] * v[2369] + e1i[0] * v[2399];
	v[2440] = e2i[1] * v[2369] + e1i[1] * v[2399];
	v[2428] = e2i[2] * v[2369] + e1i[2] * v[2399];
	v[2368] = -(v[2324] * v[281]) + v[2313] * v[282] - v[434] * v[600] + v[432] * v[603] - cp[1] * v[642];
	v[2463] = e2i[0] * v[2368] + e1i[0] * v[2398];
	v[2439] = e2i[1] * v[2368] + e1i[1] * v[2398];
	v[2427] = e2i[2] * v[2368] + e1i[2] * v[2398];
	v[2367] = -(v[2323] * v[281]) + v[2312] * v[282] - v[434] * v[599] + v[432] * v[602] - cp[1] * v[641];
	v[2462] = e2i[0] * v[2367] + e1i[0] * v[2397];
	v[2438] = e2i[1] * v[2367] + e1i[1] * v[2397];
	v[2426] = e2i[2] * v[2367] + e1i[2] * v[2397];
	v[435] = -(v[304] * v[355]) + v[349] * v[425];
	v[2406] = v[2330] * v[284] - v[2320] * v[285] + v[435] * v[610] - v[433] * v[613] - cp[0] * v[664];
	v[2405] = v[2329] * v[284] - v[2319] * v[285] + v[435] * v[609] - v[433] * v[612] - cp[0] * v[663];
	v[2404] = v[2328] * v[284] - v[2318] * v[285] + v[435] * v[608] - v[433] * v[611] - cp[0] * v[662];
	v[2376] = v[2293] - v[2330] * v[281] + v[2320] * v[282] - v[435] * v[601] + v[433] * v[604] - cp[1] * v[664];
	v[2476] = e2i[0] * v[2376] + e1i[0] * v[2406];
	v[2452] = e2i[1] * v[2376] + e1i[1] * v[2406];
	v[2434] = e2i[2] * v[2376] + e1i[2] * v[2406];
	v[2375] = v[2292] - v[2329] * v[281] + v[2319] * v[282] - v[435] * v[600] + v[433] * v[603] - cp[1] * v[663];
	v[2475] = e2i[0] * v[2375] + e1i[0] * v[2405];
	v[2451] = e2i[1] * v[2375] + e1i[1] * v[2405];
	v[2433] = e2i[2] * v[2375] + e1i[2] * v[2405];
	v[2374] = v[2291] - v[2328] * v[281] + v[2318] * v[282] - v[435] * v[599] + v[433] * v[602] - cp[1] * v[662];
	v[2474] = e2i[0] * v[2374] + e1i[0] * v[2404];
	v[2450] = e2i[1] * v[2374] + e1i[1] * v[2404];
	v[2432] = e2i[2] * v[2374] + e1i[2] * v[2404];
	v[436] = -(v[304] * v[346]) + v[339] * v[431];
	v[2414] = v[2304] - v[2334] * v[284] + v[2314] * v[286] - v[436] * v[610] + v[432] * v[616] - cp[0] * v[640];
	v[2413] = v[2303] - v[2333] * v[284] + v[2313] * v[286] - v[436] * v[609] + v[432] * v[615] - cp[0] * v[639];
	v[2412] = v[2302] - v[2332] * v[284] + v[2312] * v[286] - v[436] * v[608] + v[432] * v[614] - cp[0] * v[638];
	v[2384] = v[2334] * v[281] - v[2314] * v[283] + v[436] * v[601] - v[432] * v[607] - cp[1] * v[640];
	v[2524] = e2i[0] * v[2384] + e1i[0] * v[2414];
	v[2512] = e2i[1] * v[2384] + e1i[1] * v[2414];
	v[2488] = e2i[2] * v[2384] + e1i[2] * v[2414];
	v[2383] = v[2333] * v[281] - v[2313] * v[283] + v[436] * v[600] - v[432] * v[606] - cp[1] * v[639];
	v[2523] = e2i[0] * v[2383] + e1i[0] * v[2413];
	v[2511] = e2i[1] * v[2383] + e1i[1] * v[2413];
	v[2487] = e2i[2] * v[2383] + e1i[2] * v[2413];
	v[2382] = v[2332] * v[281] - v[2312] * v[283] + v[436] * v[599] - v[432] * v[605] - cp[1] * v[638];
	v[2522] = e2i[0] * v[2382] + e1i[0] * v[2412];
	v[2510] = e2i[1] * v[2382] + e1i[1] * v[2412];
	v[2486] = e2i[2] * v[2382] + e1i[2] * v[2412];
	v[2356] = v[2299] + v[2334] * v[285] - v[2325] * v[286] + v[436] * v[613] - v[434] * v[616] - cp[0] * v[637];
	v[2355] = v[2298] + v[2333] * v[285] - v[2324] * v[286] + v[436] * v[612] - v[434] * v[615] - cp[0] * v[636];
	v[2354] = v[2297] + v[2332] * v[285] - v[2323] * v[286] + v[436] * v[611] - v[434] * v[614] - cp[0] * v[635];
	v[2344] = -(v[2334] * v[282]) + v[2325] * v[283] - v[436] * v[604] + v[434] * v[607] - cp[1] * v[637];
	v[2596] = e2i[0] * v[2344] + e1i[0] * v[2356];
	v[2572] = e2i[1] * v[2344] + e1i[1] * v[2356];
	v[2548] = e2i[2] * v[2344] + e1i[2] * v[2356];
	v[2343] = -(v[2333] * v[282]) + v[2324] * v[283] - v[436] * v[603] + v[434] * v[606] - cp[1] * v[636];
	v[2595] = e2i[0] * v[2343] + e1i[0] * v[2355];
	v[2571] = e2i[1] * v[2343] + e1i[1] * v[2355];
	v[2547] = e2i[2] * v[2343] + e1i[2] * v[2355];
	v[2342] = -(v[2332] * v[282]) + v[2323] * v[283] - v[436] * v[602] + v[434] * v[605] - cp[1] * v[635];
	v[2570] = e2i[1] * v[2342] + e1i[1] * v[2354];
	v[2546] = e2i[2] * v[2342] + e1i[2] * v[2354];
	v[437] = -(v[304] * v[356]) + v[349] * v[426];
	v[2422] = -(v[2338] * v[284]) + v[2320] * v[286] - v[437] * v[610] + v[433] * v[616] - cp[0] * v[658];
	v[2421] = -(v[2337] * v[284]) + v[2319] * v[286] - v[437] * v[609] + v[433] * v[615] - cp[0] * v[657];
	v[2420] = -(v[2336] * v[284]) + v[2318] * v[286] - v[437] * v[608] + v[433] * v[614] - cp[0] * v[656];
	v[2392] = v[2289] + v[2338] * v[281] - v[2320] * v[283] + v[437] * v[601] - v[433] * v[607] - cp[1] * v[658];
	v[2536] = e2i[0] * v[2392] + e1i[0] * v[2422];
	v[2518] = e2i[1] * v[2392] + e1i[1] * v[2422];
	v[2500] = e2i[2] * v[2392] + e1i[2] * v[2422];
	v[2391] = v[2288] + v[2337] * v[281] - v[2319] * v[283] + v[437] * v[600] - v[433] * v[606] - cp[1] * v[657];
	v[2535] = e2i[0] * v[2391] + e1i[0] * v[2421];
	v[2517] = e2i[1] * v[2391] + e1i[1] * v[2421];
	v[2499] = e2i[2] * v[2391] + e1i[2] * v[2421];
	v[2390] = v[2287] + v[2336] * v[281] - v[2318] * v[283] + v[437] * v[599] - v[433] * v[605] - cp[1] * v[656];
	v[2534] = e2i[0] * v[2390] + e1i[0] * v[2420];
	v[2516] = e2i[1] * v[2390] + e1i[1] * v[2420];
	v[2498] = e2i[2] * v[2390] + e1i[2] * v[2420];
	v[2362] = v[2338] * v[285] - v[2330] * v[286] + v[437] * v[613] - v[435] * v[616] - cp[0] * v[655];
	v[2361] = v[2337] * v[285] - v[2329] * v[286] + v[437] * v[612] - v[435] * v[615] - cp[0] * v[654];
	v[2360] = v[2336] * v[285] - v[2328] * v[286] + v[437] * v[611] - v[435] * v[614] - cp[0] * v[653];
	v[2350] = v[2284] - v[2338] * v[282] + v[2330] * v[283] - v[437] * v[604] + v[435] * v[607] - cp[1] * v[655];
	v[2602] = e2i[0] * v[2350] + e1i[0] * v[2362];
	v[2584] = e2i[1] * v[2350] + e1i[1] * v[2362];
	v[2560] = e2i[2] * v[2350] + e1i[2] * v[2362];
	v[2349] = v[2283] - v[2337] * v[282] + v[2329] * v[283] - v[437] * v[603] + v[435] * v[606] - cp[1] * v[654];
	v[2601] = e2i[0] * v[2349] + e1i[0] * v[2361];
	v[2583] = e2i[1] * v[2349] + e1i[1] * v[2361];
	v[2559] = e2i[2] * v[2349] + e1i[2] * v[2361];
	v[2348] = v[2282] - v[2336] * v[282] + v[2328] * v[283] - v[437] * v[602] + v[435] * v[605] - cp[1] * v[653];
	v[2582] = e2i[1] * v[2348] + e1i[1] * v[2360];
	v[2558] = e2i[2] * v[2348] + e1i[2] * v[2360];
	v[438] = -(cp[1] * v[343]) + v[283] * v[434] - v[282] * v[436];
	v[439] = -(cp[1] * v[353]) + v[423] + v[283] * v[435] - v[282] * v[437];
	v[440] = -(cp[0] * v[343]) + v[428] - v[286] * v[434] + v[285] * v[436];
	v[441] = -(cp[0] * v[353]) - v[286] * v[435] + v[285] * v[437];
	v[442] = -(cp[1] * v[346]) + v[282] * v[432] - v[281] * v[434];
	v[443] = -(cp[1] * v[356]) + v[426] + v[282] * v[433] - v[281] * v[435];
	v[444] = -(cp[1] * v[345]) - v[283] * v[432] + v[281] * v[436];
	v[445] = -(cp[1] * v[355]) + v[425] - v[283] * v[433] + v[281] * v[437];
	v[446] = -(cp[0] * v[346]) + v[431] - v[285] * v[432] + v[284] * v[434];
	v[447] = -(cp[0] * v[356]) - v[285] * v[433] + v[284] * v[435];
	v[448] = -(cp[0] * v[345]) + v[430] + v[286] * v[432] - v[284] * v[436];
	v[449] = -(cp[0] * v[355]) + v[286] * v[433] - v[284] * v[437];
	v[450] = e2i[2] * v[442] + e1i[2] * v[446];
	v[451] = e2i[2] * v[443] + e1i[2] * v[447];
	v[452] = e2i[1] * v[442] + e1i[1] * v[446];
	v[2446] = v[243] * v[2440] + v[452] * v[565];
	v[2445] = v[243] * v[2439] + v[452] * v[564];
	v[497] = v[243] * v[452];
	v[453] = e2i[1] * v[443] + e1i[1] * v[447];
	v[2458] = v[243] * v[2452] + v[453] * v[565];
	v[2457] = v[243] * v[2451] + v[453] * v[564];
	v[510] = v[243] * v[453];
	v[454] = e2i[0] * v[442] + e1i[0] * v[446];
	v[2470] = v[243] * v[2464] + v[454] * v[565];
	v[2469] = v[243] * v[2463] + v[454] * v[564];
	v[500] = v[243] * v[454];
	v[455] = e2i[0] * v[443] + e1i[0] * v[447];
	v[2482] = v[243] * v[2476] + v[455] * v[565];
	v[2481] = v[243] * v[2475] + v[455] * v[564];
	v[513] = v[243] * v[455];
	v[456] = e2i[2] * v[444] + e1i[2] * v[448];
	v[4166] = v[452] - v[456];
	v[4155] = v[452] + v[456];
	v[2494] = v[243] * v[2488] + v[456] * v[565];
	v[2647] = v[2446] + v[2494];
	v[2493] = v[243] * v[2487] + v[456] * v[564];
	v[498] = v[243] * v[456];
	v[457] = e2i[2] * v[445] + e1i[2] * v[449];
	v[4162] = v[453] - v[457];
	v[4157] = v[453] + v[457];
	v[2506] = v[243] * v[2500] + v[457] * v[565];
	v[2653] = v[2458] + v[2506];
	v[2505] = v[243] * v[2499] + v[457] * v[564];
	v[511] = v[243] * v[457];
	v[458] = e2i[1] * v[444] + e1i[1] * v[448];
	v[4167] = v[450] + v[458];
	v[459] = e2i[1] * v[445] + e1i[1] * v[449];
	v[4163] = v[451] + v[459];
	v[460] = e2i[0] * v[444] + e1i[0] * v[448];
	v[2530] = v[243] * v[2524] + v[460] * v[565];
	v[505] = v[243] * v[460];
	v[461] = e2i[0] * v[445] + e1i[0] * v[449];
	v[2542] = v[243] * v[2536] + v[461] * v[565];
	v[518] = v[243] * v[461];
	v[462] = e2i[2] * v[438] + e1i[2] * v[440];
	v[4165] = v[454] + v[462];
	v[2554] = v[243] * v[2548] + v[462] * v[565];
	v[2659] = v[2470] + v[2554];
	v[2553] = v[243] * v[2547] + v[462] * v[564];
	v[501] = v[243] * v[462];
	v[463] = e2i[2] * v[439] + e1i[2] * v[441];
	v[4161] = v[455] + v[463];
	v[2566] = v[243] * v[2560] + v[463] * v[565];
	v[2665] = v[2482] + v[2566];
	v[2565] = v[243] * v[2559] + v[463] * v[564];
	v[514] = v[243] * v[463];
	v[464] = e2i[1] * v[438] + e1i[1] * v[440];
	v[4154] = v[460] + v[464];
	v[2578] = v[243] * v[2572] + v[464] * v[565];
	v[2671] = v[2530] + v[2578];
	v[2670] = v[243] * (v[2523] + v[2571]) + v[4154] * v[564];
	v[506] = v[243] * v[464];
	v[465] = e2i[1] * v[439] + e1i[1] * v[441];
	v[4156] = v[461] + v[465];
	v[2590] = v[243] * v[2584] + v[465] * v[565];
	v[2677] = v[2542] + v[2590];
	v[2676] = v[243] * (v[2535] + v[2583]) + v[4156] * v[564];
	v[519] = v[243] * v[465];
	v[466] = e2i[0] * v[438] + e1i[0] * v[440];
	v[467] = e2i[0] * v[439] + e1i[0] * v[441];
	v[2641] = v[2605] * v[450] + v[2428] * v[883];
	v[2640] = v[2604] * v[450] + v[2427] * v[883];
	v[2635] = v[2605] * v[458] + v[2512] * v[883];
	v[2629] = v[2605] * v[466] + v[2596] * v[883];
	v[2623] = v[2605] * v[451] + v[2434] * v[883];
	v[2622] = v[2604] * v[451] + v[2433] * v[883];
	v[2617] = v[2605] * v[459] + v[2518] * v[883];
	v[2611] = v[2605] * v[467] + v[2602] * v[883];
	v[521] = v[467] * v[883];
	v[520] = v[459] * v[883];
	v[515] = v[451] * v[883];
	v[508] = v[466] * v[883];
	v[507] = v[458] * v[883];
	v[502] = v[450] * v[883];
	v[474] = v[497] + v[498];
	v[475] = v[510] + v[511];
	v[480] = v[500] + v[501];
	v[481] = v[513] + v[514];
	v[483] = v[505] + v[506];
	v[484] = v[518] + v[519];
	v[487] = v[4036] * v[450] + v[4046] * v[458] + v[4045] * v[466] + v[452] * v[471] + v[454] * v[472] + v[456] * v[473]
		+ v[460] * v[478] + v[462] * v[479] + v[464] * v[482];
	v[488] = v[4036] * v[451] + v[4046] * v[459] + v[4045] * v[467] + v[453] * v[471] + v[455] * v[472] + v[457] * v[473]
		+ v[461] * v[478] + v[463] * v[479] + v[465] * v[482];
	v[2704] = v[4037] * (v[2428] * v[4036] + v[2596] * v[4045] + v[2512] * v[4046] + v[4030] * v[4155] + v[4026] * v[4165]
		+ v[460] - v[464] + d[11] * (-v[458] - v[466]) + v[2440] * v[471] + v[2464] * v[472] + v[2488] * v[473] + v[2524] * v[478]
		+ v[2548] * v[479] + v[2572] * v[482]) + v[2680] * v[487];
	v[4159] = v[2629] + v[2704];
	v[2702] = v[4037] * (v[2427] * v[4036] + v[2595] * v[4045] + v[2511] * v[4046] + v[4026] * v[4154] + v[4033] * v[4155]
		- v[454] + v[462] + d[10] * (-v[450] - v[466]) + v[2439] * v[471] + v[2463] * v[472] + v[2487] * v[473] + v[2523] * v[478]
		+ v[2547] * v[479] + v[2571] * v[482]) + v[2679] * v[487];
	v[4164] = v[2640] + v[2702];
	v[2692] = v[4037] * (v[2434] * v[4036] + v[2602] * v[4045] + v[2518] * v[4046] + v[4030] * v[4157] + v[4026] * v[4161]
		+ v[461] - v[465] + d[11] * (-v[459] - v[467]) + v[2452] * v[471] + v[2476] * v[472] + v[2500] * v[473] + v[2536] * v[478]
		+ v[2560] * v[479] + v[2584] * v[482]) + v[2680] * v[488];
	v[4158] = v[2611] + v[2692];
	v[2690] = v[4037] * (v[2433] * v[4036] + v[2601] * v[4045] + v[2517] * v[4046] + v[4026] * v[4156] + v[4033] * v[4157]
		- v[455] + v[463] + d[10] * (-v[451] - v[467]) + v[2451] * v[471] + v[2475] * v[472] + v[2499] * v[473] + v[2535] * v[478]
		+ v[2559] * v[479] + v[2583] * v[482]) + v[2679] * v[488];
	v[4160] = v[2622] + v[2690];
	v[516] = v[4037] * v[488];
	v[2721] = v[516] + v[520] + v[521];
	v[2718] = v[2721] + v[515] - v[520];
	v[2714] = v[2718] + v[520] - v[521];
	v[503] = v[4037] * v[487];
	v[2712] = v[503] + v[507] + v[508];
	v[2709] = v[2712] + v[502] - v[507];
	v[2705] = v[2709] + v[507] - v[508];
	v[2722] = v[2542] - v[2590] + 2e0*v[2721] + v[2665] * v[4026] + v[2653] * v[4030] + v[4032] * (v[2617] + v[4158]);
	v[2713] = v[2530] - v[2578] + 2e0*v[2712] + v[2659] * v[4026] + v[2647] * v[4030] + v[4032] * (v[2635] + v[4159]);
	v[2720] = -v[2482] + v[2566] + v[2677] * v[4026] + v[2653] * v[4033] + v[4029] * (v[2623] + v[4158]) + 0.5e0*v[475];
	v[2719] = -v[2481] + v[2565] + 2e0*v[2718] + v[2676] * v[4026] + (v[2457] + v[2505])*v[4033] + v[4029] * (v[4160]
		+ v[2604] * v[467] + v[2601] * v[883]);
	v[2711] = -v[2470] + v[2554] + v[2671] * v[4026] + v[2647] * v[4033] + v[4029] * (v[2641] + v[4159]) + 0.5e0*v[474];
	v[2710] = -v[2469] + v[2553] + 2e0*v[2709] + v[2670] * v[4026] + (v[2445] + v[2493])*v[4033] + v[4029] * (v[4164]
		+ v[2604] * v[466] + v[2595] * v[883]);
	v[2717] = v[2458] - v[2506] + (v[2617] + v[2623] + v[2692])*v[4027] + v[2677] * v[4030] + v[2665] * v[4033]
		+ 0.5e0*v[481];
	v[2716] = v[2457] - v[2505] + v[2676] * v[4030] + (v[2481] + v[2565])*v[4033] + 0.5e0*v[484] + v[4027] * (v[4160]
		+ v[2604] * v[459] + v[2517] * v[883]);
	v[2715] = v[243] * (v[2450] - v[2498]) + 2e0*v[2714] + v[4162] * v[563] + v[4030] * (v[243] * (v[2534] + v[2582])
		+ v[4156] * v[563]) + v[4033] * (v[243] * (v[2474] + v[2558]) + v[4161] * v[563]) + v[4027] * (v[2603] * v[4163]
			+ v[4037] * (v[2432] * v[4036] + (e2i[0] * v[2348] + e1i[0] * v[2360])*v[4045] + v[2516] * v[4046] + v[4030] * v[4156]
				+ v[4033] * v[4161] + v[4162] - d[9] * v[4163] + v[2450] * v[471] + v[2474] * v[472] + v[2498] * v[473] + v[2534] * v[478]
				+ v[2558] * v[479] + v[2582] * v[482]) + v[2678] * v[488] + (v[2432] + v[2516])*v[883]);
	v[2708] = v[2446] - v[2494] + (v[2635] + v[2641] + v[2704])*v[4027] + v[2671] * v[4030] + v[2659] * v[4033]
		+ 0.5e0*v[480];
	v[2707] = v[2445] - v[2493] + v[2670] * v[4030] + (v[2469] + v[2553])*v[4033] + 0.5e0*v[483] + v[4027] * (v[4164]
		+ v[2604] * v[458] + v[2511] * v[883]);
	v[2706] = v[243] * (v[2438] - v[2486]) + 2e0*v[2705] + v[4166] * v[563] + v[4030] * (v[243] * (v[2522] + v[2570])
		+ v[4154] * v[563]) + v[4033] * (v[243] * (v[2462] + v[2546]) + v[4165] * v[563]) + v[4027] * (v[2603] * v[4167]
			+ v[4037] * (v[2426] * v[4036] + (e2i[0] * v[2342] + e1i[0] * v[2354])*v[4045] + v[2510] * v[4046] + v[4030] * v[4154]
				+ v[4033] * v[4165] + v[4166] - d[9] * v[4167] + v[2438] * v[471] + v[2462] * v[472] + v[2486] * v[473] + v[2522] * v[478]
				+ v[2546] * v[479] + v[2570] * v[482]) + v[2678] * v[487] + (v[2426] + v[2510])*v[883]);
	v[499] = v[2705] * v[4027] + v[4033] * v[480] + v[4030] * v[483] + v[497] - v[498];
	v[4190] = v[499] / v[2037];
	v[504] = v[2709] * v[4029] + v[4033] * v[474] + v[4026] * v[483] - v[500] + v[501];
	v[4179] = v[504] / v[2037];
	v[509] = v[2712] * v[4032] + v[4030] * v[474] + v[4026] * v[480] + v[505] - v[506];
	v[4171] = v[509] / v[2037];
	v[512] = v[2714] * v[4027] + v[4033] * v[481] + v[4030] * v[484] + v[510] - v[511];
	v[4228] = -(v[419] * v[499]) + v[418] * v[512];
	v[4192] = v[417] * v[499] - v[416] * v[512];
	v[4191] = -(v[512] / v[2037]);
	v[517] = v[2718] * v[4029] + v[4033] * v[475] + v[4026] * v[484] - v[513] + v[514];
	v[4219] = -(v[419] * v[504]) + v[418] * v[517];
	v[4181] = v[417] * v[504] - v[416] * v[517];
	v[4180] = -(v[517] / v[2037]);
	v[522] = v[2721] * v[4032] + v[4030] * v[475] + v[4026] * v[481] + v[518] - v[519];
	v[4213] = -(v[419] * v[509]) + v[418] * v[522];
	v[4173] = v[417] * v[509] - v[416] * v[522];
	v[4172] = -(v[522] / v[2037]);
	v[1144] = v[1059] * v[1061] + v[1062] * v[1064] + v[1065] * v[1067] + v[3561] * v[364] + v[3552] * v[365]
		+ v[3543] * v[366];
	v[1135] = v[1059] * v[1060] + v[1062] * v[1063] + v[1065] * v[1066] + v[3559] * v[364] + v[3550] * v[365]
		+ v[3541] * v[366];
	v[1124] = v[3440] * v[364] + v[3434] * v[365] + v[3428] * v[366] + v[1059] * v[805] + v[1062] * v[807] + v[1065] * v[809];
	v[1123] = v[3443] * v[364] + v[3437] * v[365] + v[3431] * v[366] + v[1059] * v[804] + v[1062] * v[806] + v[1065] * v[808];
	v[1146] = v[1060] * v[1061] + v[1063] * v[1064] + v[1066] * v[1067] + v[3585] * v[364] + v[3577] * v[365]
		+ v[3569] * v[366];
	v[1132] = v[3441] * v[364] + v[3435] * v[365] + v[3429] * v[366] + v[1060] * v[805] + v[1063] * v[807] + v[1066] * v[809];
	v[1131] = v[3444] * v[364] + v[3438] * v[365] + v[3432] * v[366] + v[1060] * v[804] + v[1063] * v[806] + v[1066] * v[808];
	v[1141] = v[3442] * v[364] + v[3436] * v[365] + v[3430] * v[366] + v[1061] * v[805] + v[1064] * v[807] + v[1067] * v[809];
	v[1140] = v[3445] * v[364] + v[3439] * v[365] + v[3433] * v[366] + v[1061] * v[804] + v[1064] * v[806] + v[1067] * v[808];
	v[2854] = (-(v[2722] * v[416]) + v[2713] * v[417] + (v[2254] + v[4168] * v[417])*v[509] - (v[2251] + v[416] * v[4168]
		)*v[522]) / v[2037];
	v[2852] = v[2720] * v[4169];
	v[2851] = v[2711] * v[4170];
	v[4174] = v[2851] + v[2852];
	v[2853] = v[4174] + v[4187] * v[509] + v[4186] * v[522];
	v[2849] = v[2717] * v[4169];
	v[2848] = v[2708] * v[4170];
	v[4182] = v[2848] + v[2849];
	v[2850] = v[4182] + v[4178] * v[509] + v[4177] * v[522];
	v[2846] = v[4169] * v[664];
	v[2845] = v[4170] * v[643];
	v[4193] = v[2845] + v[2846];
	v[2847] = v[410] * v[4171] + v[404] * v[4172] + v[2728] * v[4173] + v[4193];
	v[2843] = v[4169] * v[658];
	v[2842] = v[4170] * v[640];
	v[4199] = v[2842] + v[2843];
	v[2844] = v[411] * v[4171] + v[405] * v[4172] + v[2726] * v[4173] + v[4199];
	v[2840] = v[4169] * v[655];
	v[2839] = v[4170] * v[637];
	v[4208] = v[2839] + v[2840];
	v[2841] = v[412] * v[4171] + v[406] * v[4172] + v[2724] * v[4173] + v[4208];
	v[2838] = v[4174] + v[4184] * v[504] + v[4183] * v[517];
	v[2837] = (-(v[2719] * v[416]) + v[2710] * v[417] + (v[2253] + v[2037] * v[4175])*v[504] - (v[2250] + v[2037] * v[4176]
		)*v[517]) / v[2037];
	v[2835] = v[2716] * v[4169];
	v[2834] = v[2707] * v[4170];
	v[4185] = v[2834] + v[2835];
	v[2836] = v[4185] + v[4178] * v[504] + v[4177] * v[517];
	v[2832] = v[4169] * v[663];
	v[2831] = v[4170] * v[642];
	v[4194] = v[2831] + v[2832];
	v[2833] = v[410] * v[4179] + v[404] * v[4180] + v[2728] * v[4181] + v[4194];
	v[2829] = v[4169] * v[657];
	v[2828] = v[4170] * v[639];
	v[4200] = v[2828] + v[2829];
	v[2830] = v[411] * v[4179] + v[405] * v[4180] + v[2726] * v[4181] + v[4200];
	v[2826] = v[4169] * v[654];
	v[2825] = v[4170] * v[636];
	v[4209] = v[2825] + v[2826];
	v[2827] = v[412] * v[4179] + v[406] * v[4180] + v[2724] * v[4181] + v[4209];
	v[2824] = v[4182] + v[4184] * v[499] + v[4183] * v[512];
	v[2823] = v[4185] + v[4187] * v[499] + v[4186] * v[512];
	v[2822] = (-(v[2715] * v[416]) + v[2706] * v[417] + (v[2252] + v[2037] * v[4188])*v[499] - (v[2249] + v[2037] * v[4189]
		)*v[512]) / v[2037];
	v[2820] = v[4169] * v[662];
	v[2819] = v[4170] * v[641];
	v[4195] = v[2819] + v[2820];
	v[2821] = v[410] * v[4190] + v[404] * v[4191] + v[2728] * v[4192] + v[4195];
	v[2817] = v[4169] * v[656];
	v[2816] = v[4170] * v[638];
	v[4201] = v[2816] + v[2817];
	v[2818] = v[411] * v[4190] + v[405] * v[4191] + v[2726] * v[4192] + v[4201];
	v[2814] = v[4169] * v[653];
	v[2813] = v[4170] * v[635];
	v[4210] = v[2813] + v[2814];
	v[2815] = v[412] * v[4190] + v[406] * v[4191] + v[2724] * v[4192] + v[4210];
	v[2812] = v[356] * v[4183] + v[346] * v[4184] + v[4193];
	v[2811] = v[356] * v[4186] + v[346] * v[4187] + v[4194];
	v[2810] = v[356] * v[4177] + v[346] * v[4178] + v[4195];
	v[2809] = (-(v[356] * v[4202]) + v[346] * v[4203]) / v[2037];
	v[2808] = (-(v[356] * v[4204]) + v[346] * v[4205]) / v[2037];
	v[2807] = (-(v[356] * v[4206]) + v[346] * v[4207]) / v[2037];
	v[2806] = v[355] * v[4183] + v[345] * v[4184] + v[4199];
	v[2805] = v[355] * v[4186] + v[345] * v[4187] + v[4200];
	v[2804] = v[355] * v[4177] + v[345] * v[4178] + v[4201];
	v[2803] = (-(v[355] * v[4202]) + v[345] * v[4203]) / v[2037];
	v[2802] = (-(v[355] * v[4204]) + v[345] * v[4205]) / v[2037];
	v[2801] = (-(v[355] * v[4206]) + v[345] * v[4207]) / v[2037];
	v[2800] = v[353] * v[4183] + v[343] * v[4184] + v[4208];
	v[2799] = v[353] * v[4186] + v[343] * v[4187] + v[4209];
	v[2798] = v[353] * v[4177] + v[343] * v[4178] + v[4210];
	v[2797] = (-(v[353] * v[4202]) + v[343] * v[4203]) / v[2037];
	v[2796] = (-(v[353] * v[4204]) + v[343] * v[4205]) / v[2037];
	v[2795] = (-(v[353] * v[4206]) + v[343] * v[4207]) / v[2037];
	v[2794] = (v[2722] * v[418] - v[2713] * v[419] - (v[2260] + v[4168] * v[419])*v[509] + (v[2257] + v[4168] * v[418]
		)*v[522]) / v[2037];
	v[2792] = v[2720] * v[4211];
	v[2791] = v[2711] * v[4212];
	v[4214] = v[2791] + v[2792];
	v[2793] = v[4214] + v[4225] * v[509] + v[4224] * v[522];
	v[2789] = v[2717] * v[4211];
	v[2788] = v[2708] * v[4212];
	v[4220] = v[2788] + v[2789];
	v[2790] = v[4220] + v[4218] * v[509] + v[4217] * v[522];
	v[2786] = v[4211] * v[664];
	v[2785] = v[4212] * v[643];
	v[4229] = v[2785] + v[2786];
	v[2787] = -(v[413] * v[4171]) - v[407] * v[4172] + v[2728] * v[4213] + v[4229];
	v[3603] = v[2787] * v[804] + v[2847] * v[805];
	v[3596] = v[2787] * v[806] + v[2847] * v[807];
	v[3589] = v[2787] * v[808] + v[2847] * v[809];
	v[2783] = v[4211] * v[658];
	v[2782] = v[4212] * v[640];
	v[4232] = v[2782] + v[2783];
	v[2784] = -(v[414] * v[4171]) - v[408] * v[4172] + v[2726] * v[4213] + v[4232];
	v[3602] = v[2784] * v[804] + v[2844] * v[805];
	v[3595] = v[2784] * v[806] + v[2844] * v[807];
	v[3588] = v[2784] * v[808] + v[2844] * v[809];
	v[2780] = v[4211] * v[655];
	v[2779] = v[4212] * v[637];
	v[4241] = v[2779] + v[2780];
	v[2781] = -(v[415] * v[4171]) - v[409] * v[4172] + v[2724] * v[4213] + v[4241];
	v[3601] = v[2781] * v[804] + v[2841] * v[805];
	v[3594] = v[2781] * v[806] + v[2841] * v[807];
	v[3587] = v[2781] * v[808] + v[2841] * v[809];
	v[2778] = v[4214] + v[4222] * v[504] + v[4221] * v[517];
	v[2777] = (v[2719] * v[418] - v[2710] * v[419] - (v[2259] + v[2037] * v[4215])*v[504] + (v[2256] + v[2037] * v[4216]
		)*v[517]) / v[2037];
	v[2775] = v[2716] * v[4211];
	v[2774] = v[2707] * v[4212];
	v[4223] = v[2774] + v[2775];
	v[2776] = v[4223] + v[4218] * v[504] + v[4217] * v[517];
	v[2772] = v[4211] * v[663];
	v[2771] = v[4212] * v[642];
	v[4230] = v[2771] + v[2772];
	v[2773] = -(v[413] * v[4179]) - v[407] * v[4180] + v[2728] * v[4219] + v[4230];
	v[3581] = v[2773] * v[804] + v[2833] * v[805];
	v[3573] = v[2773] * v[806] + v[2833] * v[807];
	v[3565] = v[2773] * v[808] + v[2833] * v[809];
	v[2769] = v[4211] * v[657];
	v[2768] = v[4212] * v[639];
	v[4233] = v[2768] + v[2769];
	v[2770] = -(v[414] * v[4179]) - v[408] * v[4180] + v[2726] * v[4219] + v[4233];
	v[3580] = v[2770] * v[804] + v[2830] * v[805];
	v[3572] = v[2770] * v[806] + v[2830] * v[807];
	v[3564] = v[2770] * v[808] + v[2830] * v[809];
	v[2766] = v[4211] * v[654];
	v[2765] = v[4212] * v[636];
	v[4242] = v[2765] + v[2766];
	v[2767] = -(v[415] * v[4179]) - v[409] * v[4180] + v[2724] * v[4219] + v[4242];
	v[3579] = v[2767] * v[804] + v[2827] * v[805];
	v[3571] = v[2767] * v[806] + v[2827] * v[807];
	v[3563] = v[2767] * v[808] + v[2827] * v[809];
	v[2764] = v[4220] + v[4222] * v[499] + v[4221] * v[512];
	v[2763] = v[4223] + v[4225] * v[499] + v[4224] * v[512];
	v[2762] = (v[2715] * v[418] - v[2706] * v[419] - (v[2258] + v[2037] * v[4226])*v[499] + (v[2255] + v[2037] * v[4227]
		)*v[512]) / v[2037];
	v[2760] = v[4211] * v[662];
	v[2759] = v[4212] * v[641];
	v[4231] = v[2759] + v[2760];
	v[2761] = -(v[413] * v[4190]) - v[407] * v[4191] + v[2728] * v[4228] + v[4231];
	v[3556] = v[2761] * v[804] + v[2821] * v[805];
	v[3547] = v[2761] * v[806] + v[2821] * v[807];
	v[3538] = v[2761] * v[808] + v[2821] * v[809];
	v[2757] = v[4211] * v[656];
	v[2756] = v[4212] * v[638];
	v[4234] = v[2756] + v[2757];
	v[2758] = -(v[414] * v[4190]) - v[408] * v[4191] + v[2726] * v[4228] + v[4234];
	v[3555] = v[2758] * v[804] + v[2818] * v[805];
	v[3546] = v[2758] * v[806] + v[2818] * v[807];
	v[3537] = v[2758] * v[808] + v[2818] * v[809];
	v[2754] = v[4211] * v[653];
	v[2753] = v[4212] * v[635];
	v[4243] = v[2753] + v[2754];
	v[2755] = -(v[415] * v[4190]) - v[409] * v[4191] + v[2724] * v[4228] + v[4243];
	v[3554] = v[2755] * v[804] + v[2815] * v[805];
	v[3545] = v[2755] * v[806] + v[2815] * v[807];
	v[3536] = v[2755] * v[808] + v[2815] * v[809];
	v[2752] = v[356] * v[4221] + v[346] * v[4222] + v[4229];
	v[2751] = v[356] * v[4224] + v[346] * v[4225] + v[4230];
	v[2750] = v[356] * v[4217] + v[346] * v[4218] + v[4231];
	v[2749] = (v[356] * v[4235] - v[346] * v[4236]) / v[2037];
	v[3532] = v[2749] * v[804] + v[2809] * v[805];
	v[3526] = v[2749] * v[806] + v[2809] * v[807];
	v[3520] = v[2749] * v[808] + v[2809] * v[809];
	v[2748] = (v[356] * v[4237] - v[346] * v[4238]) / v[2037];
	v[3531] = v[2748] * v[804] + v[2808] * v[805];
	v[3525] = v[2748] * v[806] + v[2808] * v[807];
	v[3519] = v[2748] * v[808] + v[2808] * v[809];
	v[2747] = (v[356] * v[4239] - v[346] * v[4240]) / v[2037];
	v[3530] = v[2747] * v[804] + v[2807] * v[805];
	v[3524] = v[2747] * v[806] + v[2807] * v[807];
	v[3518] = v[2747] * v[808] + v[2807] * v[809];
	v[2746] = v[355] * v[4221] + v[345] * v[4222] + v[4232];
	v[2745] = v[355] * v[4224] + v[345] * v[4225] + v[4233];
	v[2744] = v[355] * v[4217] + v[345] * v[4218] + v[4234];
	v[2743] = (v[355] * v[4235] - v[345] * v[4236]) / v[2037];
	v[3514] = v[2743] * v[804] + v[2803] * v[805];
	v[3508] = v[2743] * v[806] + v[2803] * v[807];
	v[3502] = v[2743] * v[808] + v[2803] * v[809];
	v[2742] = (v[355] * v[4237] - v[345] * v[4238]) / v[2037];
	v[3513] = v[2742] * v[804] + v[2802] * v[805];
	v[3507] = v[2742] * v[806] + v[2802] * v[807];
	v[3501] = v[2742] * v[808] + v[2802] * v[809];
	v[2741] = (v[355] * v[4239] - v[345] * v[4240]) / v[2037];
	v[3512] = v[2741] * v[804] + v[2801] * v[805];
	v[3506] = v[2741] * v[806] + v[2801] * v[807];
	v[3500] = v[2741] * v[808] + v[2801] * v[809];
	v[2740] = v[353] * v[4221] + v[343] * v[4222] + v[4241];
	v[2739] = v[353] * v[4224] + v[343] * v[4225] + v[4242];
	v[2738] = v[353] * v[4217] + v[343] * v[4218] + v[4243];
	v[2737] = (v[353] * v[4235] - v[343] * v[4236]) / v[2037];
	v[3496] = v[2737] * v[804] + v[2797] * v[805];
	v[3490] = v[2737] * v[806] + v[2797] * v[807];
	v[3484] = v[2737] * v[808] + v[2797] * v[809];
	v[2736] = (v[353] * v[4237] - v[343] * v[4238]) / v[2037];
	v[3495] = v[2736] * v[804] + v[2796] * v[805];
	v[3489] = v[2736] * v[806] + v[2796] * v[807];
	v[3483] = v[2736] * v[808] + v[2796] * v[809];
	v[2735] = (v[353] * v[4239] - v[343] * v[4240]) / v[2037];
	v[3494] = v[2735] * v[804] + v[2795] * v[805];
	v[3488] = v[2735] * v[806] + v[2795] * v[807];
	v[3482] = v[2735] * v[808] + v[2795] * v[809];
	v[2075] = -(v[1863] * v[4171]);
	v[2074] = -(v[1861] * v[4172]);
	v[2076] = -v[2074] - v[2075] + v[2040] * v[4173] + v[417] * v[4244] + v[416] * v[4245];
	v[2073] = (-(v[2254] * v[416]) + v[2251] * v[417] + v[4250] * v[509] - v[4249] * v[522]) / v[2037];
	v[2071] = -(v[1863] * v[4179]);
	v[2070] = -(v[1861] * v[4180]);
	v[2072] = -v[2070] - v[2071] + v[2040] * v[4181] + v[417] * v[4247] + v[416] * v[4248];
	v[2069] = (-(v[2253] * v[416]) + v[2250] * v[417] + v[4250] * v[504] - v[4249] * v[517]) / v[2037];
	v[2067] = -(v[1863] * v[4190]);
	v[2066] = -(v[1861] * v[4191]);
	v[2068] = -v[2066] - v[2067] + v[2040] * v[4192] + v[417] * v[4251] + v[416] * v[4252];
	v[2065] = (-(v[2252] * v[416]) + v[2249] * v[417] + v[4250] * v[499] - v[4249] * v[512]) / v[2037];
	v[2063] = v[346] * v[4257];
	v[2062] = v[356] * v[4258];
	v[2064] = (-v[4253] - v[4254] - v[2037] * (v[2062] + v[2063] - v[2040] * v[4255] + v[2040] * v[4256])) / v[2037];
	v[2061] = (v[1862] * v[346] - v[1860] * v[356] - v[410] * v[416] + v[404] * v[417] + v[4246] * v[4255] - v[4246] * v[4256])
		/ v[2037];
	v[2059] = v[345] * v[4257];
	v[2058] = v[355] * v[4258];
	v[2060] = (-v[4259] - v[4260] - v[2037] * (v[2058] + v[2059] - v[2040] * v[4261] + v[2040] * v[4262])) / v[2037];
	v[2057] = (v[1862] * v[345] - v[1860] * v[355] - v[411] * v[416] + v[405] * v[417] + v[4246] * v[4261] - v[4246] * v[4262])
		/ v[2037];
	v[2055] = v[343] * v[4257];
	v[2054] = v[353] * v[4258];
	v[2056] = (-v[4263] - v[4264] - v[2037] * (v[2054] + v[2055] - v[2040] * v[4265] + v[2040] * v[4266])) / v[2037];
	v[2053] = (v[1862] * v[343] - v[1860] * v[353] - v[412] * v[416] + v[406] * v[417] + v[4246] * v[4265] - v[4246] * v[4266])
		/ v[2037];
	v[2052] = (v[2260] * v[418] - v[2257] * v[419] - v[4271] * v[509] + v[4270] * v[522]) / v[2037];
	v[2051] = v[2074] + v[2075] + v[419] * (v[4268] - v[2038] * v[509]) + v[418] * (v[4269] + v[2038] * v[522]);
	v[2050] = (v[2259] * v[418] - v[2256] * v[419] - v[4271] * v[504] + v[4270] * v[517]) / v[2037];
	v[2049] = v[2070] + v[2071] + v[419] * (v[4272] - v[2038] * v[504]) + v[418] * (v[4273] + v[2038] * v[517]);
	v[2048] = (v[2258] * v[418] - v[2255] * v[419] - v[4271] * v[499] + v[4270] * v[512]) / v[2037];
	v[2047] = v[2066] + v[2067] + v[419] * (v[4274] - v[2038] * v[499]) + v[418] * (v[4275] + v[2038] * v[512]);
	v[2046] = (-(v[1865] * v[346]) + v[1864] * v[356] + v[413] * v[418] - v[407] * v[419] - v[346] * v[4276] + v[356] * v[4277]
		) / v[2037];
	v[2045] = v[2062] + v[2063] + v[4278] / v[2037] - v[4279] / v[2037] + v[2038] * v[4286];
	v[2044] = (-(v[1865] * v[345]) + v[1864] * v[355] + v[414] * v[418] - v[408] * v[419] - v[345] * v[4276] + v[355] * v[4277]
		) / v[2037];
	v[2043] = v[2058] + v[2059] + v[4280] / v[2037] - v[4281] / v[2037] + v[2038] * v[4285];
	v[2042] = (-(v[1865] * v[343]) + v[1864] * v[353] + v[415] * v[418] - v[409] * v[419] - v[343] * v[4276] + v[353] * v[4277]
		) / v[2037];
	v[2041] = v[2054] + v[2055] + v[4282] / v[2037] - v[4283] / v[2037] + v[2038] * v[4284];
	v[699] = v[4284] / v[2037];
	v[3290] = v[2949] * v[699];
	v[3288] = v[2948] * v[699];
	v[3283] = v[2965] * v[699];
	v[700] = v[4285] / v[2037];
	v[3302] = v[2949] * v[700];
	v[701] = v[4286] / v[2037];
	v[702] = v[4228] / v[2037];
	v[703] = v[4219] / v[2037];
	v[704] = v[4213] / v[2037];
	v[705] = (v[4265] - v[4266]) / v[2037];
	v[3840] = v[3818] * v[699] + v[3819] * v[705] + v[2042] * v[804] + v[2056] * v[805];
	v[3839] = v[3820] * v[699] + v[3818] * v[705] + v[2041] * v[804] + v[2053] * v[805];
	v[3838] = v[3815] * v[699] + v[3816] * v[705] + v[2042] * v[806] + v[2056] * v[807];
	v[3837] = v[3817] * v[699] + v[3815] * v[705] + v[2041] * v[806] + v[2053] * v[807];
	v[3836] = v[3812] * v[699] + v[3813] * v[705] + v[2042] * v[808] + v[2056] * v[809];
	v[3835] = v[3814] * v[699] + v[3812] * v[705] + v[2041] * v[808] + v[2053] * v[809];
	v[3499] = v[3445] * v[699] + v[3442] * v[705] + v[2740] * v[804] + v[2800] * v[805];
	v[3498] = v[3444] * v[699] + v[3441] * v[705] + v[2739] * v[804] + v[2799] * v[805];
	v[3497] = v[3443] * v[699] + v[3440] * v[705] + v[2738] * v[804] + v[2798] * v[805];
	v[3493] = v[3439] * v[699] + v[3436] * v[705] + v[2740] * v[806] + v[2800] * v[807];
	v[3492] = v[3438] * v[699] + v[3435] * v[705] + v[2739] * v[806] + v[2799] * v[807];
	v[3491] = v[3437] * v[699] + v[3434] * v[705] + v[2738] * v[806] + v[2798] * v[807];
	v[3487] = v[3433] * v[699] + v[3430] * v[705] + v[2740] * v[808] + v[2800] * v[809];
	v[3486] = v[3432] * v[699] + v[3429] * v[705] + v[2739] * v[808] + v[2799] * v[809];
	v[3485] = v[3431] * v[699] + v[3428] * v[705] + v[2738] * v[808] + v[2798] * v[809];
	v[3289] = v[2940] * v[705];
	v[4333] = v[3289] + v[3290];
	v[3287] = v[2939] * v[705];
	v[4363] = v[3287] + v[3288];
	v[3282] = v[2957] * v[705];
	v[4362] = v[3282] + v[3283];
	v[1557] = v[699] * v[808] + v[705] * v[809];
	v[1556] = v[699] * v[806] + v[705] * v[807];
	v[1555] = 1e0 + v[699] * v[804] + v[705] * v[805];
	v[706] = (v[4261] - v[4262]) / v[2037];
	v[3846] = v[3818] * v[700] + v[3819] * v[706] + v[2044] * v[804] + v[2060] * v[805];
	v[3845] = v[3820] * v[700] + v[3818] * v[706] + v[2043] * v[804] + v[2057] * v[805];
	v[3844] = v[3815] * v[700] + v[3816] * v[706] + v[2044] * v[806] + v[2060] * v[807];
	v[3843] = v[3817] * v[700] + v[3815] * v[706] + v[2043] * v[806] + v[2057] * v[807];
	v[3842] = v[3812] * v[700] + v[3813] * v[706] + v[2044] * v[808] + v[2060] * v[809];
	v[3841] = v[3814] * v[700] + v[3812] * v[706] + v[2043] * v[808] + v[2057] * v[809];
	v[3517] = v[3445] * v[700] + v[3442] * v[706] + v[2746] * v[804] + v[2806] * v[805];
	v[3516] = v[3444] * v[700] + v[3441] * v[706] + v[2745] * v[804] + v[2805] * v[805];
	v[3515] = v[3443] * v[700] + v[3440] * v[706] + v[2744] * v[804] + v[2804] * v[805];
	v[3511] = v[3439] * v[700] + v[3436] * v[706] + v[2746] * v[806] + v[2806] * v[807];
	v[3510] = v[3438] * v[700] + v[3435] * v[706] + v[2745] * v[806] + v[2805] * v[807];
	v[3509] = v[3437] * v[700] + v[3434] * v[706] + v[2744] * v[806] + v[2804] * v[807];
	v[3505] = v[3433] * v[700] + v[3430] * v[706] + v[2746] * v[808] + v[2806] * v[809];
	v[3504] = v[3432] * v[700] + v[3429] * v[706] + v[2745] * v[808] + v[2805] * v[809];
	v[3503] = v[3431] * v[700] + v[3428] * v[706] + v[2744] * v[808] + v[2804] * v[809];
	v[3301] = v[2940] * v[706];
	v[4364] = v[3301] + v[3302];
	v[1563] = v[700] * v[808] + v[706] * v[809];
	v[1562] = 1e0 + v[700] * v[806] + v[706] * v[807];
	v[1561] = v[700] * v[804] + v[706] * v[805];
	v[707] = (v[4255] - v[4256]) / v[2037];
	v[3852] = v[3818] * v[701] + v[3819] * v[707] + v[2046] * v[804] + v[2064] * v[805];
	v[3851] = v[3820] * v[701] + v[3818] * v[707] + v[2045] * v[804] + v[2061] * v[805];
	v[3850] = v[3815] * v[701] + v[3816] * v[707] + v[2046] * v[806] + v[2064] * v[807];
	v[3849] = v[3817] * v[701] + v[3815] * v[707] + v[2045] * v[806] + v[2061] * v[807];
	v[3848] = v[3812] * v[701] + v[3813] * v[707] + v[2046] * v[808] + v[2064] * v[809];
	v[3847] = v[3814] * v[701] + v[3812] * v[707] + v[2045] * v[808] + v[2061] * v[809];
	v[3535] = v[3445] * v[701] + v[3442] * v[707] + v[2752] * v[804] + v[2812] * v[805];
	v[3534] = v[3444] * v[701] + v[3441] * v[707] + v[2751] * v[804] + v[2811] * v[805];
	v[3533] = v[3443] * v[701] + v[3440] * v[707] + v[2750] * v[804] + v[2810] * v[805];
	v[3529] = v[3439] * v[701] + v[3436] * v[707] + v[2752] * v[806] + v[2812] * v[807];
	v[3528] = v[3438] * v[701] + v[3435] * v[707] + v[2751] * v[806] + v[2811] * v[807];
	v[3527] = v[3437] * v[701] + v[3434] * v[707] + v[2750] * v[806] + v[2810] * v[807];
	v[3523] = v[3433] * v[701] + v[3430] * v[707] + v[2752] * v[808] + v[2812] * v[809];
	v[3522] = v[3432] * v[701] + v[3429] * v[707] + v[2751] * v[808] + v[2811] * v[809];
	v[3521] = v[3431] * v[701] + v[3428] * v[707] + v[2750] * v[808] + v[2810] * v[809];
	v[1573] = 1e0 + v[701] * v[808] + v[707] * v[809];
	v[1572] = v[701] * v[806] + v[707] * v[807];
	v[1571] = v[701] * v[804] + v[707] * v[805];
	v[708] = v[4192] / v[2037];
	v[3858] = v[3440] + v[3818] * v[702] + v[3819] * v[708] + v[2048] * v[804] + v[2068] * v[805];
	v[3857] = v[3443] + v[3820] * v[702] + v[3818] * v[708] + v[2047] * v[804] + v[2065] * v[805];
	v[3856] = v[3434] + v[3815] * v[702] + v[3816] * v[708] + v[2048] * v[806] + v[2068] * v[807];
	v[3855] = v[3437] + v[3817] * v[702] + v[3815] * v[708] + v[2047] * v[806] + v[2065] * v[807];
	v[3854] = v[3428] + v[3812] * v[702] + v[3813] * v[708] + v[2048] * v[808] + v[2068] * v[809];
	v[3853] = v[3431] + v[3814] * v[702] + v[3812] * v[708] + v[2047] * v[808] + v[2065] * v[809];
	v[3562] = v[3561] + v[3445] * v[702] + v[3442] * v[708] + v[2764] * v[804] + v[2824] * v[805];
	v[3560] = v[3559] + v[3444] * v[702] + v[3441] * v[708] + v[2763] * v[804] + v[2823] * v[805];
	v[3558] = v[3557] + v[3443] * v[702] + v[3440] * v[708] + v[2762] * v[804] + v[2822] * v[805];
	v[3553] = v[3552] + v[3439] * v[702] + v[3436] * v[708] + v[2764] * v[806] + v[2824] * v[807];
	v[3551] = v[3550] + v[3438] * v[702] + v[3435] * v[708] + v[2763] * v[806] + v[2823] * v[807];
	v[3549] = v[3548] + v[3437] * v[702] + v[3434] * v[708] + v[2762] * v[806] + v[2822] * v[807];
	v[3544] = v[3543] + v[3433] * v[702] + v[3430] * v[708] + v[2764] * v[808] + v[2824] * v[809];
	v[3542] = v[3541] + v[3432] * v[702] + v[3429] * v[708] + v[2763] * v[808] + v[2823] * v[809];
	v[3540] = v[3539] + v[3431] * v[702] + v[3428] * v[708] + v[2762] * v[808] + v[2822] * v[809];
	v[1586] = v[1065] + v[702] * v[808] + v[708] * v[809];
	v[1585] = v[1062] + v[702] * v[806] + v[708] * v[807];
	v[1584] = v[1059] + v[702] * v[804] + v[708] * v[805];
	v[4302] = (v[1584] * v[358] + v[1585] * v[361] + v[1586] * v[363])*v[4297];
	v[709] = v[4181] / v[2037];
	v[3864] = v[3441] + v[3818] * v[703] + v[3819] * v[709] + v[2050] * v[804] + v[2072] * v[805];
	v[3863] = v[3444] + v[3820] * v[703] + v[3818] * v[709] + v[2049] * v[804] + v[2069] * v[805];
	v[3862] = v[3435] + v[3815] * v[703] + v[3816] * v[709] + v[2050] * v[806] + v[2072] * v[807];
	v[3861] = v[3438] + v[3817] * v[703] + v[3815] * v[709] + v[2049] * v[806] + v[2069] * v[807];
	v[3860] = v[3429] + v[3812] * v[703] + v[3813] * v[709] + v[2050] * v[808] + v[2072] * v[809];
	v[3859] = v[3432] + v[3814] * v[703] + v[3812] * v[709] + v[2049] * v[808] + v[2069] * v[809];
	v[3586] = v[3585] + v[3445] * v[703] + v[3442] * v[709] + v[2778] * v[804] + v[2838] * v[805];
	v[3584] = v[3583] + v[3444] * v[703] + v[3441] * v[709] + v[2777] * v[804] + v[2837] * v[805];
	v[3582] = v[3559] + v[3443] * v[703] + v[3440] * v[709] + v[2776] * v[804] + v[2836] * v[805];
	v[3578] = v[3577] + v[3439] * v[703] + v[3436] * v[709] + v[2778] * v[806] + v[2838] * v[807];
	v[3576] = v[3575] + v[3438] * v[703] + v[3435] * v[709] + v[2777] * v[806] + v[2837] * v[807];
	v[3574] = v[3550] + v[3437] * v[703] + v[3434] * v[709] + v[2776] * v[806] + v[2836] * v[807];
	v[3570] = v[3569] + v[3433] * v[703] + v[3430] * v[709] + v[2778] * v[808] + v[2838] * v[809];
	v[3568] = v[3567] + v[3432] * v[703] + v[3429] * v[709] + v[2777] * v[808] + v[2837] * v[809];
	v[3566] = v[3541] + v[3431] * v[703] + v[3428] * v[709] + v[2776] * v[808] + v[2836] * v[809];
	v[1593] = v[1066] + v[703] * v[808] + v[709] * v[809];
	v[1592] = v[1063] + v[703] * v[806] + v[709] * v[807];
	v[1591] = v[1060] + v[703] * v[804] + v[709] * v[805];
	v[4300] = (v[1591] * v[358] + v[1592] * v[361] + v[1593] * v[363])*v[4297];
	v[710] = v[4173] / v[2037];
	v[3870] = v[3442] + v[3818] * v[704] + v[3819] * v[710] + v[2052] * v[804] + v[2076] * v[805];
	v[3869] = v[3445] + v[3820] * v[704] + v[3818] * v[710] + v[2051] * v[804] + v[2073] * v[805];
	v[3868] = v[3436] + v[3815] * v[704] + v[3816] * v[710] + v[2052] * v[806] + v[2076] * v[807];
	v[3867] = v[3439] + v[3817] * v[704] + v[3815] * v[710] + v[2051] * v[806] + v[2073] * v[807];
	v[3866] = v[3430] + v[3812] * v[704] + v[3813] * v[710] + v[2052] * v[808] + v[2076] * v[809];
	v[3865] = v[3433] + v[3814] * v[704] + v[3812] * v[710] + v[2051] * v[808] + v[2073] * v[809];
	v[3607] = v[3606] + v[3445] * v[704] + v[3442] * v[710] + v[2794] * v[804] + v[2854] * v[805];
	v[3605] = v[3585] + v[3444] * v[704] + v[3441] * v[710] + v[2793] * v[804] + v[2853] * v[805];
	v[3604] = v[3561] + v[3443] * v[704] + v[3440] * v[710] + v[2790] * v[804] + v[2850] * v[805];
	v[3600] = v[3599] + v[3439] * v[704] + v[3436] * v[710] + v[2794] * v[806] + v[2854] * v[807];
	v[3598] = v[3577] + v[3438] * v[704] + v[3435] * v[710] + v[2793] * v[806] + v[2853] * v[807];
	v[3597] = v[3552] + v[3437] * v[704] + v[3434] * v[710] + v[2790] * v[806] + v[2850] * v[807];
	v[3593] = v[3592] + v[3433] * v[704] + v[3430] * v[710] + v[2794] * v[808] + v[2854] * v[809];
	v[3591] = v[3569] + v[3432] * v[704] + v[3429] * v[710] + v[2793] * v[808] + v[2853] * v[809];
	v[3590] = v[3543] + v[3431] * v[704] + v[3428] * v[710] + v[2790] * v[808] + v[2850] * v[809];
	v[1600] = v[1067] + v[704] * v[808] + v[710] * v[809];
	v[1599] = v[1064] + v[704] * v[806] + v[710] * v[807];
	v[1598] = v[1061] + v[704] * v[804] + v[710] * v[805];
	v[4298] = (v[1598] * v[358] + v[1599] * v[361] + v[1600] * v[363])*v[4297];
	v[1149] = v[319] * v[334] - v[318] * v[336];
	v[1150] = -(v[319] * v[331]) + v[316] * v[336];
	v[1151] = v[318] * v[331] - v[316] * v[334];
	v[3269] = v[1004] * v[1149] + v[1151] * v[980] + v[1150] * v[992];
	v[3272] = v[3125] * v[316] + v[1149] * v[3269] + v[3124] * v[331];
	v[3271] = v[3125] * v[318] + v[1150] * v[3269] + v[3124] * v[334];
	v[3270] = v[3125] * v[319] + v[1151] * v[3269] + v[3124] * v[336];
	v[3262] = v[1002] * v[1149] + v[1151] * v[978] + v[1150] * v[990];
	v[3268] = v[3123] * v[316] + v[1149] * v[3262] + v[3121] * v[331];
	v[3266] = v[3123] * v[318] + v[1150] * v[3262] + v[3121] * v[334];
	v[3264] = v[3123] * v[319] + v[1151] * v[3262] + v[3121] * v[336];
	v[3261] = v[1000] * v[1149] + v[1151] * v[976] + v[1150] * v[988];
	v[3267] = v[3122] * v[316] + v[1149] * v[3261] + v[3120] * v[331];
	v[3265] = v[3122] * v[318] + v[1150] * v[3261] + v[3120] * v[334];
	v[3263] = v[3122] * v[319] + v[1151] * v[3261] + v[3120] * v[336];
	v[3251] = v[1151] * v[974] + v[1150] * v[986] + v[1149] * v[998];
	v[3260] = v[3119] * v[316] + v[1149] * v[3251] + v[3116] * v[331];
	v[3257] = v[3119] * v[318] + v[1150] * v[3251] + v[3116] * v[334];
	v[3254] = v[3119] * v[319] + v[1151] * v[3251] + v[3116] * v[336];
	v[3250] = v[1151] * v[972] + v[1150] * v[984] + v[1149] * v[996];
	v[3259] = v[3118] * v[316] + v[1149] * v[3250] + v[3115] * v[331];
	v[3256] = v[3118] * v[318] + v[1150] * v[3250] + v[3115] * v[334];
	v[3253] = v[3118] * v[319] + v[1151] * v[3250] + v[3115] * v[336];
	v[3249] = v[1151] * v[970] + v[1150] * v[982] + v[1149] * v[994];
	v[3258] = v[3117] * v[316] + v[1149] * v[3249] + v[3114] * v[331];
	v[3255] = v[3117] * v[318] + v[1150] * v[3249] + v[3114] * v[334];
	v[3252] = v[3117] * v[319] + v[1151] * v[3249] + v[3114] * v[336];
	v[3245] = v[1151] * v[3152] + v[1150] * v[3182] + v[1149] * v[3198];
	v[3248] = v[316] * v[3200] + v[1149] * v[3245] + v[3199] * v[331];
	v[3247] = v[318] * v[3200] + v[1150] * v[3245] + v[3199] * v[334];
	v[3246] = v[319] * v[3200] + v[1151] * v[3245] + v[3199] * v[336];
	v[3238] = v[1151] * v[3151] + v[1150] * v[3181] + v[1149] * v[3193];
	v[3244] = v[316] * v[3197] + v[1149] * v[3238] + v[3195] * v[331];
	v[3242] = v[318] * v[3197] + v[1150] * v[3238] + v[3195] * v[334];
	v[3240] = v[319] * v[3197] + v[1151] * v[3238] + v[3195] * v[336];
	v[3237] = v[1151] * v[3150] + v[1150] * v[3180] + v[1149] * v[3192];
	v[3243] = v[316] * v[3196] + v[1149] * v[3237] + v[3194] * v[331];
	v[3241] = v[318] * v[3196] + v[1150] * v[3237] + v[3194] * v[334];
	v[3239] = v[319] * v[3196] + v[1151] * v[3237] + v[3194] * v[336];
	v[3227] = v[1151] * v[3149] + v[1150] * v[3179] + v[1149] * v[3185];
	v[3236] = v[316] * v[3191] + v[1149] * v[3227] + v[3188] * v[331];
	v[3233] = v[318] * v[3191] + v[1150] * v[3227] + v[3188] * v[334];
	v[3230] = v[319] * v[3191] + v[1151] * v[3227] + v[3188] * v[336];
	v[3226] = v[1151] * v[3148] + v[1150] * v[3178] + v[1149] * v[3184];
	v[3235] = v[316] * v[3190] + v[1149] * v[3226] + v[3187] * v[331];
	v[3232] = v[318] * v[3190] + v[1150] * v[3226] + v[3187] * v[334];
	v[3229] = v[319] * v[3190] + v[1151] * v[3226] + v[3187] * v[336];
	v[3225] = v[1151] * v[3147] + v[1150] * v[3177] + v[1149] * v[3183];
	v[3234] = v[316] * v[3189] + v[1149] * v[3225] + v[3186] * v[331];
	v[3231] = v[318] * v[3189] + v[1150] * v[3225] + v[3186] * v[334];
	v[3228] = v[3189] * v[319] + v[1151] * v[3225] + v[3186] * v[336];
	v[3221] = v[1069] * v[1149] + v[1075] * v[1150] + v[1081] * v[1151];
	v[3224] = v[3137] * v[316] + v[1149] * v[3221] + v[3136] * v[331];
	v[3223] = v[3137] * v[318] + v[1150] * v[3221] + v[3136] * v[334];
	v[3222] = v[3137] * v[319] + v[1151] * v[3221] + v[3136] * v[336];
	v[3214] = v[1087] * v[1149] + v[1091] * v[1150] + v[1095] * v[1151];
	v[3220] = v[3135] * v[316] + v[1149] * v[3214] + v[3133] * v[331];
	v[3218] = v[3135] * v[318] + v[1150] * v[3214] + v[3133] * v[334];
	v[3216] = v[3135] * v[319] + v[1151] * v[3214] + v[3133] * v[336];
	v[3213] = v[1071] * v[1149] + v[1077] * v[1150] + v[1083] * v[1151];
	v[3219] = v[3134] * v[316] + v[1149] * v[3213] + v[3132] * v[331];
	v[3217] = v[3134] * v[318] + v[1150] * v[3213] + v[3132] * v[334];
	v[3215] = v[3134] * v[319] + v[1151] * v[3213] + v[3132] * v[336];
	v[3203] = v[1099] * v[1149] + v[1101] * v[1150] + v[1103] * v[1151];
	v[3212] = v[3131] * v[316] + v[1149] * v[3203] + v[3128] * v[331];
	v[3209] = v[3131] * v[318] + v[1150] * v[3203] + v[3128] * v[334];
	v[3206] = v[3131] * v[319] + v[1151] * v[3203] + v[3128] * v[336];
	v[3202] = v[1089] * v[1149] + v[1093] * v[1150] + v[1097] * v[1151];
	v[3211] = v[3130] * v[316] + v[1149] * v[3202] + v[3127] * v[331];
	v[3208] = v[3130] * v[318] + v[1150] * v[3202] + v[3127] * v[334];
	v[3205] = v[3130] * v[319] + v[1151] * v[3202] + v[3127] * v[336];
	v[3201] = v[1073] * v[1149] + v[1079] * v[1150] + v[1085] * v[1151];
	v[3210] = v[3129] * v[316] + v[1149] * v[3201] + v[3126] * v[331];
	v[3207] = v[3129] * v[318] + v[1150] * v[3201] + v[3126] * v[334];
	v[3204] = v[3129] * v[319] + v[1151] * v[3201] + v[3126] * v[336];
	v[2917] = v[1151] * v[720] + v[1150] * v[732] + v[1149] * v[744];
	v[2926] = v[1149] * v[2917] + v[2179] * v[316] + v[2176] * v[331];
	v[2923] = v[1150] * v[2917] + v[2179] * v[318] + v[2176] * v[334];
	v[2920] = v[1151] * v[2917] + v[2179] * v[319] + v[2176] * v[336];
	v[2916] = v[1151] * v[724] + v[1150] * v[736] + v[1149] * v[748];
	v[2925] = v[1149] * v[2916] + v[2178] * v[316] + v[2175] * v[331];
	v[2922] = v[1150] * v[2916] + v[2178] * v[318] + v[2175] * v[334];
	v[2919] = v[1151] * v[2916] + v[2178] * v[319] + v[2175] * v[336];
	v[2915] = v[1151] * v[728] + v[1150] * v[740] + v[1149] * v[752];
	v[2924] = v[1149] * v[2915] + v[2177] * v[316] + v[2174] * v[331];
	v[2921] = v[1150] * v[2915] + v[2177] * v[318] + v[2174] * v[334];
	v[2918] = v[1151] * v[2915] + v[2177] * v[319] + v[2174] * v[336];
	v[2905] = v[1151] * v[722] + v[1150] * v[734] + v[1149] * v[746];
	v[2914] = v[1149] * v[2905] + v[2173] * v[316] + v[2170] * v[331];
	v[2911] = v[1150] * v[2905] + v[2173] * v[318] + v[2170] * v[334];
	v[2908] = v[1151] * v[2905] + v[2173] * v[319] + v[2170] * v[336];
	v[2904] = v[1151] * v[726] + v[1150] * v[738] + v[1149] * v[750];
	v[2913] = v[1149] * v[2904] + v[2172] * v[316] + v[2169] * v[331];
	v[2910] = v[1150] * v[2904] + v[2172] * v[318] + v[2169] * v[334];
	v[2907] = v[1151] * v[2904] + v[2172] * v[319] + v[2169] * v[336];
	v[2903] = v[1151] * v[730] + v[1150] * v[742] + v[1149] * v[754];
	v[2912] = v[1149] * v[2903] + v[2171] * v[316] + v[2168] * v[331];
	v[2909] = v[1150] * v[2903] + v[2171] * v[318] + v[2168] * v[334];
	v[2906] = v[1151] * v[2903] + v[2171] * v[319] + v[2168] * v[336];
	v[2893] = v[1151] * v[2224] + v[1150] * v[2230] + v[1149] * v[2242];
	v[2902] = v[1149] * v[2893] + v[2248] * v[316] + v[2245] * v[331];
	v[2899] = v[1150] * v[2893] + v[2248] * v[318] + v[2245] * v[334];
	v[2896] = v[1151] * v[2893] + v[2248] * v[319] + v[2245] * v[336];
	v[2892] = v[1151] * v[2223] + v[1150] * v[2229] + v[1149] * v[2241];
	v[2901] = v[1149] * v[2892] + v[2247] * v[316] + v[2244] * v[331];
	v[2898] = v[1150] * v[2892] + v[2247] * v[318] + v[2244] * v[334];
	v[2895] = v[1151] * v[2892] + v[2247] * v[319] + v[2244] * v[336];
	v[2891] = v[1151] * v[2222] + v[1150] * v[2228] + v[1149] * v[2240];
	v[2900] = v[1149] * v[2891] + v[2246] * v[316] + v[2243] * v[331];
	v[2897] = v[1150] * v[2891] + v[2246] * v[318] + v[2243] * v[334];
	v[2894] = v[1151] * v[2891] + v[2246] * v[319] + v[2243] * v[336];
	v[2881] = v[1151] * v[2221] + v[1150] * v[2227] + v[1149] * v[2233];
	v[2890] = v[1149] * v[2881] + v[2239] * v[316] + v[2236] * v[331];
	v[2887] = v[1150] * v[2881] + v[2239] * v[318] + v[2236] * v[334];
	v[2884] = v[1151] * v[2881] + v[2239] * v[319] + v[2236] * v[336];
	v[2880] = v[1151] * v[2220] + v[1150] * v[2226] + v[1149] * v[2232];
	v[2889] = v[1149] * v[2880] + v[2238] * v[316] + v[2235] * v[331];
	v[2886] = v[1150] * v[2880] + v[2238] * v[318] + v[2235] * v[334];
	v[2883] = v[1151] * v[2880] + v[2238] * v[319] + v[2235] * v[336];
	v[2879] = v[1151] * v[2219] + v[1150] * v[2225] + v[1149] * v[2231];
	v[2888] = v[1149] * v[2879] + v[2237] * v[316] + v[2234] * v[331];
	v[2885] = v[1150] * v[2879] + v[2237] * v[318] + v[2234] * v[334];
	v[2882] = v[1151] * v[2879] + v[2237] * v[319] + v[2234] * v[336];
	v[2869] = v[1149] * v[822] + v[1150] * v[824] + v[1151] * v[826];
	v[2878] = v[1149] * v[2869] + v[2218] * v[316] + v[2215] * v[331];
	v[2875] = v[1150] * v[2869] + v[2218] * v[318] + v[2215] * v[334];
	v[2872] = v[1151] * v[2869] + v[2218] * v[319] + v[2215] * v[336];
	v[2868] = v[1149] * v[816] + v[1150] * v[818] + v[1151] * v[820];
	v[2877] = v[1149] * v[2868] + v[2217] * v[316] + v[2214] * v[331];
	v[2874] = v[1150] * v[2868] + v[2217] * v[318] + v[2214] * v[334];
	v[2871] = v[1151] * v[2868] + v[2217] * v[319] + v[2214] * v[336];
	v[2867] = v[1149] * v[810] + v[1150] * v[812] + v[1151] * v[814];
	v[2876] = v[1149] * v[2867] + v[2216] * v[316] + v[2213] * v[331];
	v[2873] = v[1150] * v[2867] + v[2216] * v[318] + v[2213] * v[334];
	v[2870] = v[1151] * v[2867] + v[2216] * v[319] + v[2213] * v[336];
	v[2857] = v[1149] * v[823] + v[1150] * v[825] + v[1151] * v[827];
	v[2866] = v[1149] * v[2857] + v[2212] * v[316] + v[2209] * v[331];
	v[2863] = v[1150] * v[2857] + v[2212] * v[318] + v[2209] * v[334];
	v[2860] = v[1151] * v[2857] + v[2212] * v[319] + v[2209] * v[336];
	v[2856] = v[1149] * v[817] + v[1150] * v[819] + v[1151] * v[821];
	v[2865] = v[1149] * v[2856] + v[2211] * v[316] + v[2208] * v[331];
	v[2862] = v[1150] * v[2856] + v[2211] * v[318] + v[2208] * v[334];
	v[2859] = v[1151] * v[2856] + v[2211] * v[319] + v[2208] * v[336];
	v[2855] = v[1149] * v[811] + v[1150] * v[813] + v[1151] * v[815];
	v[2864] = v[1149] * v[2855] + v[2210] * v[316] + v[2207] * v[331];
	v[2861] = v[1150] * v[2855] + v[2210] * v[318] + v[2207] * v[334];
	v[2858] = v[1151] * v[2855] + v[2210] * v[319] + v[2207] * v[336];
	v[2109] = v[1149] * v[1771] + v[1150] * v[1773] + v[1151] * v[1775];
	v[2112] = v[1149] * v[2109] + v[1777] * v[316] + v[1776] * v[331];
	v[2111] = v[1150] * v[2109] + v[1777] * v[318] + v[1776] * v[334];
	v[2110] = v[1151] * v[2109] + v[1777] * v[319] + v[1776] * v[336];
	v[2102] = v[1149] * v[1757] + v[1150] * v[1761] + v[1151] * v[1765];
	v[2108] = v[1149] * v[2102] + v[1769] * v[316] + v[1767] * v[331];
	v[2106] = v[1150] * v[2102] + v[1769] * v[318] + v[1767] * v[334];
	v[2104] = v[1151] * v[2102] + v[1769] * v[319] + v[1767] * v[336];
	v[2101] = v[1149] * v[1755] + v[1150] * v[1759] + v[1151] * v[1763];
	v[2107] = v[1149] * v[2101] + v[1768] * v[316] + v[1766] * v[331];
	v[2105] = v[1150] * v[2101] + v[1768] * v[318] + v[1766] * v[334];
	v[2103] = v[1151] * v[2101] + v[1768] * v[319] + v[1766] * v[336];
	v[2097] = v[1151] * v[1847] + v[1150] * v[1850] + v[1149] * v[1857];
	v[2100] = v[1149] * v[2097] + v[1859] * v[316] + v[1858] * v[331];
	v[2099] = v[1150] * v[2097] + v[1859] * v[318] + v[1858] * v[334];
	v[2098] = v[1151] * v[2097] + v[1859] * v[319] + v[1858] * v[336];
	v[2090] = v[1151] * v[1846] + v[1150] * v[1849] + v[1149] * v[1852];
	v[2096] = v[1149] * v[2090] + v[1856] * v[316] + v[1854] * v[331];
	v[2094] = v[1150] * v[2090] + v[1856] * v[318] + v[1854] * v[334];
	v[2092] = v[1151] * v[2090] + v[1856] * v[319] + v[1854] * v[336];
	v[2089] = v[1151] * v[1845] + v[1150] * v[1848] + v[1149] * v[1851];
	v[2095] = v[1149] * v[2089] + v[1855] * v[316] + v[1853] * v[331];
	v[2093] = v[1150] * v[2089] + v[1855] * v[318] + v[1853] * v[334];
	v[2091] = v[1151] * v[2089] + v[1855] * v[319] + v[1853] * v[336];
	v[2085] = v[1151] * v[1828] + v[1150] * v[1833] + v[1149] * v[1842];
	v[2088] = v[1149] * v[2085] + v[1844] * v[316] + v[1843] * v[331];
	v[2087] = v[1150] * v[2085] + v[1844] * v[318] + v[1843] * v[334];
	v[2086] = v[1151] * v[2085] + v[1844] * v[319] + v[1843] * v[336];
	v[2078] = v[1151] * v[1826] + v[1150] * v[1831] + v[1149] * v[1836];
	v[2084] = v[1149] * v[2078] + v[1840] * v[316] + v[1838] * v[331];
	v[2082] = v[1150] * v[2078] + v[1840] * v[318] + v[1838] * v[334];
	v[2080] = v[1151] * v[2078] + v[1840] * v[319] + v[1838] * v[336];
	v[2077] = v[1151] * v[1824] + v[1150] * v[1829] + v[1149] * v[1834];
	v[2083] = v[1149] * v[2077] + v[1839] * v[316] + v[1837] * v[331];
	v[2081] = v[1150] * v[2077] + v[1839] * v[318] + v[1837] * v[334];
	v[2079] = v[1151] * v[2077] + v[1839] * v[319] + v[1837] * v[336];
	v[1278] = v[1149] * v[675] + v[1150] * v[681] + v[1151] * v[687];
	v[1380] = v[1151] * v[1278] + v[1269] * v[319] + v[1287] * v[336];
	v[1371] = v[1150] * v[1278] + v[1269] * v[318] + v[1287] * v[334];
	v[1362] = v[1149] * v[1278] + v[1269] * v[316] + v[1287] * v[331];
	v[1277] = v[1149] * v[673] + v[1150] * v[679] + v[1151] * v[685];
	v[1379] = v[1151] * v[1277] + v[1268] * v[319] + v[1286] * v[336];
	v[1370] = v[1150] * v[1277] + v[1268] * v[318] + v[1286] * v[334];
	v[1361] = v[1149] * v[1277] + v[1268] * v[316] + v[1286] * v[331];
	v[1276] = v[1149] * v[671] + v[1150] * v[677] + v[1151] * v[683];
	v[1378] = v[1151] * v[1276] + v[1267] * v[319] + v[1285] * v[336];
	v[1369] = v[1150] * v[1276] + v[1267] * v[318] + v[1285] * v[334];
	v[1360] = v[1149] * v[1276] + v[1267] * v[316] + v[1285] * v[331];
	v[1275] = v[1149] * v[1254] + v[1150] * v[1257] + v[1151] * v[1260];
	v[1377] = v[1151] * v[1275] + v[1266] * v[319] + v[1284] * v[336];
	v[1368] = v[1150] * v[1275] + v[1266] * v[318] + v[1284] * v[334];
	v[1359] = v[1149] * v[1275] + v[1266] * v[316] + v[1284] * v[331];
	v[1274] = v[1149] * v[1253] + v[1150] * v[1256] + v[1151] * v[1259];
	v[1376] = v[1151] * v[1274] + v[1265] * v[319] + v[1283] * v[336];
	v[1367] = v[1150] * v[1274] + v[1265] * v[318] + v[1283] * v[334];
	v[1358] = v[1149] * v[1274] + v[1265] * v[316] + v[1283] * v[331];
	v[1273] = v[1149] * v[1252] + v[1150] * v[1255] + v[1151] * v[1258];
	v[1375] = v[1151] * v[1273] + v[1264] * v[319] + v[1282] * v[336];
	v[1366] = v[1150] * v[1273] + v[1264] * v[318] + v[1282] * v[334];
	v[1357] = v[1149] * v[1273] + v[1264] * v[316] + v[1282] * v[331];
	v[1272] = v[1149] * v[637] + v[1150] * v[640] + v[1151] * v[643];
	v[1374] = v[1151] * v[1272] + v[1263] * v[319] + v[1281] * v[336];
	v[1365] = v[1150] * v[1272] + v[1263] * v[318] + v[1281] * v[334];
	v[1356] = v[1149] * v[1272] + v[1263] * v[316] + v[1281] * v[331];
	v[1271] = v[1149] * v[636] + v[1150] * v[639] + v[1151] * v[642];
	v[1373] = v[1151] * v[1271] + v[1262] * v[319] + v[1280] * v[336];
	v[1364] = v[1150] * v[1271] + v[1262] * v[318] + v[1280] * v[334];
	v[1355] = v[1149] * v[1271] + v[1262] * v[316] + v[1280] * v[331];
	v[1270] = v[1149] * v[635] + v[1150] * v[638] + v[1151] * v[641];
	v[1372] = v[1151] * v[1270] + v[1261] * v[319] + v[1279] * v[336];
	v[1363] = v[1150] * v[1270] + v[1261] * v[318] + v[1279] * v[334];
	v[1354] = v[1149] * v[1270] + v[1261] * v[316] + v[1279] * v[331];
	v[1208] = v[1149] * v[545] + v[1150] * v[547] + v[1151] * v[549];
	v[1234] = v[1151] * v[1208] + v[1202] * v[319] + v[1214] * v[336];
	v[1228] = v[1150] * v[1208] + v[1202] * v[318] + v[1214] * v[334];
	v[1222] = v[1149] * v[1208] + v[1202] * v[316] + v[1214] * v[331];
	v[1207] = v[1149] * v[544] + v[1150] * v[546] + v[1151] * v[548];
	v[1233] = v[1151] * v[1207] + v[1201] * v[319] + v[1213] * v[336];
	v[1227] = v[1150] * v[1207] + v[1201] * v[318] + v[1213] * v[334];
	v[1221] = v[1149] * v[1207] + v[1201] * v[316] + v[1213] * v[331];
	v[1206] = v[1149] * v[1192] + v[1150] * v[1194] + v[1151] * v[1196];
	v[1232] = v[1151] * v[1206] + v[1200] * v[319] + v[1212] * v[336];
	v[1226] = v[1150] * v[1206] + v[1200] * v[318] + v[1212] * v[334];
	v[1220] = v[1149] * v[1206] + v[1200] * v[316] + v[1212] * v[331];
	v[1205] = v[1149] * v[1191] + v[1150] * v[1193] + v[1151] * v[1195];
	v[1231] = v[1151] * v[1205] + v[1199] * v[319] + v[1211] * v[336];
	v[1225] = v[1150] * v[1205] + v[1199] * v[318] + v[1211] * v[334];
	v[1219] = v[1149] * v[1205] + v[1199] * v[316] + v[1211] * v[331];
	v[1204] = v[1151] * v[407] + v[1150] * v[408] + v[1149] * v[409];
	v[1230] = v[1151] * v[1204] + v[1198] * v[319] + v[1210] * v[336];
	v[1224] = v[1150] * v[1204] + v[1198] * v[318] + v[1210] * v[334];
	v[1218] = v[1149] * v[1204] + v[1198] * v[316] + v[1210] * v[331];
	v[1203] = v[1151] * v[404] + v[1150] * v[405] + v[1149] * v[406];
	v[1229] = v[1151] * v[1203] + v[1197] * v[319] + v[1209] * v[336];
	v[1223] = v[1150] * v[1203] + v[1197] * v[318] + v[1209] * v[334];
	v[1217] = v[1149] * v[1203] + v[1197] * v[316] + v[1209] * v[331];
	v[1152] = v[346] * v[361] - v[345] * v[363];
	v[1153] = -(v[346] * v[358]) + v[343] * v[363];
	v[1154] = v[345] * v[358] - v[343] * v[361];
	v[1155] = v[316] * v[343] + v[318] * v[345] + v[319] * v[346];
	v[1156] = v[1152] * v[316] + v[1153] * v[318] + v[1154] * v[319];
	v[1157] = v[316] * v[358] + v[318] * v[361] + v[319] * v[363];
	v[1158] = v[1149] * v[343] + v[1150] * v[345] + v[1151] * v[346];
	v[1159] = v[1149] * v[1152] + v[1150] * v[1153] + v[1151] * v[1154];
	v[1160] = v[1149] * v[358] + v[1150] * v[361] + v[1151] * v[363];
	v[1161] = v[331] * v[343] + v[334] * v[345] + v[336] * v[346];
	v[1162] = v[1152] * v[331] + v[1153] * v[334] + v[1154] * v[336];
	v[1163] = v[331] * v[358] + v[334] * v[361] + v[336] * v[363];
	v[1165] = d[0] + ci[0] * (e1i[0] - v[281]) + ci[1] * (e2i[0] - v[284]) - v[4056] + v[296] * (-v[287] + v[4074]) + (-
		(v[1294] * v[334]) - v[1291] * v[336] + 0.5e0*v[4287])*v[4288] + xPi[0];
	v[1167] = d[1] + ci[0] * (e1i[1] - v[282]) + ci[1] * (e2i[1] - v[285]) - v[4057] + v[296] * (-v[288] + v[4075]) + v[4288] * (
		-(v[1313] * v[331]) - v[1311] * v[336] - 0.5e0*v[4289]) + xPi[1];
	v[1169] = d[2] + ci[0] * (e1i[2] - v[283]) + ci[1] * (e2i[2] - v[286]) - v[4058] + v[296] * (-v[289] + v[4076]) + v[4288] * (
		-(v[1333] * v[331]) - v[1331] * v[334] - 0.5e0*v[4290]) + xPi[2];
	v[3278] = v[1069] * v[1165] + v[1075] * v[1167] + v[1081] * v[1169] - v[3032] * v[358] - v[3026] * v[361]
		- v[3020] * v[363] + v[1306] * v[4334] + v[1325] * v[4337] + v[1345] * v[4342];
	v[3277] = v[1087] * v[1165] + v[1091] * v[1167] + v[1095] * v[1169] - v[3031] * v[358] - v[3025] * v[361]
		- v[3019] * v[363] + v[1307] * v[4335] + v[1326] * v[4340] + v[1346] * v[4343];
	v[3276] = v[1071] * v[1165] + v[1077] * v[1167] + v[1083] * v[1169] - v[3030] * v[358] - v[3024] * v[361]
		- v[3018] * v[363] - v[1307] * v[671] - v[1306] * v[673] - v[1326] * v[677] - v[1325] * v[679] - v[1346] * v[683]
		- v[1345] * v[685];
	v[3275] = v[1099] * v[1165] + v[1101] * v[1167] + v[1103] * v[1169] - v[3029] * v[358] - v[3023] * v[361]
		- v[3017] * v[363] + v[1308] * v[4336] + v[1327] * v[4341] + v[1347] * v[4346];
	v[3274] = v[1089] * v[1165] + v[1093] * v[1167] + v[1097] * v[1169] - v[3028] * v[358] - v[3022] * v[361]
		- v[3016] * v[363] - v[1308] * v[673] - v[1307] * v[675] - v[1327] * v[679] - v[1326] * v[681] - v[1347] * v[685]
		- v[1346] * v[687];
	v[3273] = v[1073] * v[1165] + v[1079] * v[1167] + v[1085] * v[1169] - v[3027] * v[358] - v[3021] * v[361]
		- v[3015] * v[363] - v[1308] * v[671] - v[1306] * v[675] - v[1327] * v[677] - v[1325] * v[681] - v[1347] * v[683]
		- v[1345] * v[687];
	v[2938] = -(v[1308] * v[544]) - v[1327] * v[546] - v[1347] * v[548] + v[1165] * v[822] + v[1167] * v[824]
		+ v[1169] * v[826];
	v[2937] = -(v[1307] * v[544]) - v[1326] * v[546] - v[1346] * v[548] + v[1165] * v[816] + v[1167] * v[818]
		+ v[1169] * v[820];
	v[2936] = -(v[1306] * v[544]) - v[1325] * v[546] - v[1345] * v[548] + v[1165] * v[810] + v[1167] * v[812]
		+ v[1169] * v[814];
	v[2932] = -(v[1308] * v[545]) - v[1327] * v[547] - v[1347] * v[549] + v[1165] * v[823] + v[1167] * v[825]
		+ v[1169] * v[827];
	v[2931] = -(v[1307] * v[545]) - v[1326] * v[547] - v[1346] * v[549] + v[1165] * v[817] + v[1167] * v[819]
		+ v[1169] * v[821];
	v[2930] = -(v[1306] * v[545]) - v[1325] * v[547] - v[1345] * v[549] + v[1165] * v[811] + v[1167] * v[813]
		+ v[1169] * v[815];
	v[2115] = v[1169] * v[1828] + v[1167] * v[1833] + v[1165] * v[1842];
	v[2114] = v[1169] * v[1826] + v[1167] * v[1831] + v[1165] * v[1836];
	v[2113] = v[1169] * v[1824] + v[1167] * v[1829] + v[1165] * v[1834];
	v[1353] = -(v[1308] * v[358]) - v[1327] * v[361] - v[1347] * v[363] + v[1165] * v[675] + v[1167] * v[681]
		+ v[1169] * v[687];
	v[1352] = -(v[1307] * v[358]) - v[1326] * v[361] - v[1346] * v[363] + v[1165] * v[673] + v[1167] * v[679]
		+ v[1169] * v[685];
	v[1351] = -(v[1306] * v[358]) - v[1325] * v[361] - v[1345] * v[363] + v[1165] * v[671] + v[1167] * v[677]
		+ v[1169] * v[683];
	v[1216] = v[1165] * v[545] + v[1167] * v[547] + v[1169] * v[549];
	v[1215] = v[1165] * v[544] + v[1167] * v[546] + v[1169] * v[548];
	v[1177] = v[1165] * v[358] + v[1167] * v[361] + v[1169] * v[363];
	v[2986] = gti[0] * (v[1149] * v[2902] + v[2926] * v[316] + v[2878] * v[331]) + gti[1] * (v[1150] * v[2902]
		+ v[2926] * v[318] + v[2878] * v[334]) + gti[2] * (v[1151] * v[2902] + v[2926] * v[319] + v[2878] * v[336])
		- v[2938] * v[358] - v[1353] * v[544] - v[1215] * v[675] - v[1177] * v[822];
	v[2985] = gti[0] * (v[1149] * v[2901] + v[2925] * v[316] + v[2877] * v[331]) + gti[1] * (v[1150] * v[2901]
		+ v[2925] * v[318] + v[2877] * v[334]) + gti[2] * (v[1151] * v[2901] + v[2925] * v[319] + v[2877] * v[336])
		- v[2937] * v[358] - v[1352] * v[544] - v[1215] * v[673] - v[1177] * v[816];
	v[2984] = gti[0] * (v[1149] * v[2900] + v[2924] * v[316] + v[2876] * v[331]) + gti[1] * (v[1150] * v[2900]
		+ v[2924] * v[318] + v[2876] * v[334]) + gti[2] * (v[1151] * v[2900] + v[2924] * v[319] + v[2876] * v[336])
		- v[2936] * v[358] - v[1351] * v[544] - v[1215] * v[671] - v[1177] * v[810];
	v[2979] = gti[0] * (v[1149] * v[2890] + v[2914] * v[316] + v[2866] * v[331]) + gti[1] * (v[1150] * v[2890]
		+ v[2914] * v[318] + v[2866] * v[334]) + gti[2] * (v[1151] * v[2890] + v[2914] * v[319] + v[2866] * v[336])
		- v[2932] * v[358] - v[1353] * v[545] - v[1216] * v[675] - v[1177] * v[823];
	v[2978] = gti[0] * (v[1149] * v[2889] + v[2913] * v[316] + v[2865] * v[331]) + gti[1] * (v[1150] * v[2889]
		+ v[2913] * v[318] + v[2865] * v[334]) + gti[2] * (v[1151] * v[2889] + v[2913] * v[319] + v[2865] * v[336])
		- v[2931] * v[358] - v[1352] * v[545] - v[1216] * v[673] - v[1177] * v[817];
	v[2977] = gti[0] * (v[1149] * v[2888] + v[2912] * v[316] + v[2864] * v[331]) + gti[1] * (v[1150] * v[2888]
		+ v[2912] * v[318] + v[2864] * v[334]) + gti[2] * (v[1151] * v[2888] + v[2912] * v[319] + v[2864] * v[336])
		- v[2930] * v[358] - v[1351] * v[545] - v[1216] * v[671] - v[1177] * v[811];
	v[2972] = gti[0] * (v[1149] * v[2899] + v[2923] * v[316] + v[2875] * v[331]) + gti[1] * (v[1150] * v[2899]
		+ v[2923] * v[318] + v[2875] * v[334]) + gti[2] * (v[1151] * v[2899] + v[2923] * v[319] + v[2875] * v[336])
		- v[2938] * v[361] - v[1353] * v[546] - v[1215] * v[681] - v[1177] * v[824];
	v[2971] = gti[0] * (v[1149] * v[2898] + v[2922] * v[316] + v[2874] * v[331]) + gti[1] * (v[1150] * v[2898]
		+ v[2922] * v[318] + v[2874] * v[334]) + gti[2] * (v[1151] * v[2898] + v[2922] * v[319] + v[2874] * v[336])
		- v[2937] * v[361] - v[1352] * v[546] - v[1215] * v[679] - v[1177] * v[818];
	v[2970] = gti[0] * (v[1149] * v[2897] + v[2921] * v[316] + v[2873] * v[331]) + gti[1] * (v[1150] * v[2897]
		+ v[2921] * v[318] + v[2873] * v[334]) + gti[2] * (v[1151] * v[2897] + v[2921] * v[319] + v[2873] * v[336])
		- v[2936] * v[361] - v[1351] * v[546] - v[1215] * v[677] - v[1177] * v[812];
	v[2964] = gti[0] * (v[1149] * v[2887] + v[2911] * v[316] + v[2863] * v[331]) + gti[1] * (v[1150] * v[2887]
		+ v[2911] * v[318] + v[2863] * v[334]) + gti[2] * (v[1151] * v[2887] + v[2911] * v[319] + v[2863] * v[336])
		- v[2932] * v[361] - v[1353] * v[547] - v[1216] * v[681] - v[1177] * v[825];
	v[2963] = gti[0] * (v[1149] * v[2886] + v[2910] * v[316] + v[2862] * v[331]) + gti[1] * (v[1150] * v[2886]
		+ v[2910] * v[318] + v[2862] * v[334]) + gti[2] * (v[1151] * v[2886] + v[2910] * v[319] + v[2862] * v[336])
		- v[2931] * v[361] - v[1352] * v[547] - v[1216] * v[679] - v[1177] * v[819];
	v[2962] = gti[0] * (v[1149] * v[2885] + v[2909] * v[316] + v[2861] * v[331]) + gti[1] * (v[1150] * v[2885]
		+ v[2909] * v[318] + v[2861] * v[334]) + gti[2] * (v[1151] * v[2885] + v[2909] * v[319] + v[2861] * v[336])
		- v[2930] * v[361] - v[1351] * v[547] - v[1216] * v[677] - v[1177] * v[813];
	v[2956] = gti[0] * (v[1149] * v[2896] + v[2920] * v[316] + v[2872] * v[331]) + gti[1] * (v[1150] * v[2896]
		+ v[2920] * v[318] + v[2872] * v[334]) + gti[2] * (v[1151] * v[2896] + v[2920] * v[319] + v[2872] * v[336])
		- v[2938] * v[363] - v[1353] * v[548] - v[1215] * v[687] - v[1177] * v[826];
	v[2955] = gti[0] * (v[1149] * v[2895] + v[2919] * v[316] + v[2871] * v[331]) + gti[1] * (v[1150] * v[2895]
		+ v[2919] * v[318] + v[2871] * v[334]) + gti[2] * (v[1151] * v[2895] + v[2919] * v[319] + v[2871] * v[336])
		- v[2937] * v[363] - v[1352] * v[548] - v[1215] * v[685] - v[1177] * v[820];
	v[2954] = gti[0] * (v[1149] * v[2894] + v[2918] * v[316] + v[2870] * v[331]) + gti[1] * (v[1150] * v[2894]
		+ v[2918] * v[318] + v[2870] * v[334]) + gti[2] * (v[1151] * v[2894] + v[2918] * v[319] + v[2870] * v[336])
		- v[2936] * v[363] - v[1351] * v[548] - v[1215] * v[683] - v[1177] * v[814];
	v[2947] = gti[0] * (v[1149] * v[2884] + v[2908] * v[316] + v[2860] * v[331]) + gti[1] * (v[1150] * v[2884]
		+ v[2908] * v[318] + v[2860] * v[334]) + gti[2] * (v[1151] * v[2884] + v[2908] * v[319] + v[2860] * v[336])
		- v[2932] * v[363] - v[1353] * v[549] - v[1216] * v[687] - v[1177] * v[827];
	v[2946] = gti[0] * (v[1149] * v[2883] + v[2907] * v[316] + v[2859] * v[331]) + gti[1] * (v[1150] * v[2883]
		+ v[2907] * v[318] + v[2859] * v[334]) + gti[2] * (v[1151] * v[2883] + v[2907] * v[319] + v[2859] * v[336])
		- v[2931] * v[363] - v[1352] * v[549] - v[1216] * v[685] - v[1177] * v[821];
	v[2945] = gti[0] * (v[1149] * v[2882] + v[2906] * v[316] + v[2858] * v[331]) + gti[1] * (v[1150] * v[2882]
		+ v[2906] * v[318] + v[2858] * v[334]) + gti[2] * (v[1151] * v[2882] + v[2906] * v[319] + v[2858] * v[336])
		- v[2930] * v[363] - v[1351] * v[549] - v[1216] * v[683] - v[1177] * v[815];
	v[2124] = -(v[1177] * v[1842]) + gti[0] * (v[1149] * v[2100] + v[2112] * v[316] + v[2088] * v[331]) + gti[1] *
		(v[1150] * v[2100] + v[2112] * v[318] + v[2088] * v[334]) + gti[2] * (v[1151] * v[2100] + v[2112] * v[319]
			+ v[2088] * v[336]) - v[2115] * v[358] + v[1215] * v[4291];
	v[2123] = -(v[1177] * v[1836]) + gti[0] * (v[1149] * v[2096] + v[2108] * v[316] + v[2084] * v[331]) + gti[1] *
		(v[1150] * v[2096] + v[2108] * v[318] + v[2084] * v[334]) + gti[2] * (v[1151] * v[2096] + v[2108] * v[319]
			+ v[2084] * v[336]) - v[2114] * v[358] + v[1216] * v[4292];
	v[2122] = -(v[1177] * v[1834]) + gti[0] * (v[1149] * v[2095] + v[2107] * v[316] + v[2083] * v[331]) + gti[1] *
		(v[1150] * v[2095] + v[2107] * v[318] + v[2083] * v[334]) + gti[2] * (v[1151] * v[2095] + v[2107] * v[319]
			+ v[2083] * v[336]) - v[2113] * v[358] - v[1216] * v[544] - v[1215] * v[545];
	v[2121] = -(v[1177] * v[1833]) + gti[0] * (v[1149] * v[2099] + v[2111] * v[316] + v[2087] * v[331]) + gti[1] *
		(v[1150] * v[2099] + v[2111] * v[318] + v[2087] * v[334]) + gti[2] * (v[1151] * v[2099] + v[2111] * v[319]
			+ v[2087] * v[336]) - v[2115] * v[361] - v[1215] * v[4148];
	v[2120] = -(v[1177] * v[1831]) + gti[0] * (v[1149] * v[2094] + v[2106] * v[316] + v[2082] * v[331]) + gti[1] *
		(v[1150] * v[2094] + v[2106] * v[318] + v[2082] * v[334]) + gti[2] * (v[1151] * v[2094] + v[2106] * v[319]
			+ v[2082] * v[336]) - v[2114] * v[361] - v[1216] * v[4150];
	v[2119] = -(v[1177] * v[1829]) + gti[0] * (v[1149] * v[2093] + v[2105] * v[316] + v[2081] * v[331]) + gti[1] *
		(v[1150] * v[2093] + v[2105] * v[318] + v[2081] * v[334]) + gti[2] * (v[1151] * v[2093] + v[2105] * v[319]
			+ v[2081] * v[336]) - v[2113] * v[361] - v[1216] * v[546] - v[1215] * v[547];
	v[2118] = -(v[1177] * v[1828]) + gti[0] * (v[1149] * v[2098] + v[2110] * v[316] + v[2086] * v[331]) + gti[1] *
		(v[1150] * v[2098] + v[2110] * v[318] + v[2086] * v[334]) + gti[2] * (v[1151] * v[2098] + v[2110] * v[319]
			+ v[2086] * v[336]) - v[2115] * v[363] - v[1215] * v[4293];
	v[2117] = -(v[1177] * v[1826]) + gti[0] * (v[1149] * v[2092] + v[2104] * v[316] + v[2080] * v[331]) + gti[1] *
		(v[1150] * v[2092] + v[2104] * v[318] + v[2080] * v[334]) + gti[2] * (v[1151] * v[2092] + v[2104] * v[319]
			+ v[2080] * v[336]) - v[2114] * v[363] - v[1216] * v[4294];
	v[2116] = -(v[1177] * v[1824]) + gti[0] * (v[1149] * v[2091] + v[2103] * v[316] + v[2079] * v[331]) + gti[1] *
		(v[1150] * v[2091] + v[2103] * v[318] + v[2079] * v[334]) + gti[2] * (v[1151] * v[2091] + v[2103] * v[319]
			+ v[2079] * v[336]) - v[2113] * v[363] - v[1216] * v[548] - v[1215] * v[549];
	v[1414] = gti[0] * (v[1149] * v[1232] + v[1230] * v[316] + v[1234] * v[331]) + gti[1] * (v[1150] * v[1232]
		+ v[1230] * v[318] + v[1234] * v[334]) + gti[2] * (v[1151] * v[1232] + v[1230] * v[319] + v[1234] * v[336])
		- v[1216] * v[363] - v[1177] * v[549];
	v[1413] = gti[0] * (v[1149] * v[1231] + v[1229] * v[316] + v[1233] * v[331]) + gti[1] * (v[1150] * v[1231]
		+ v[1229] * v[318] + v[1233] * v[334]) + gti[2] * (v[1151] * v[1231] + v[1229] * v[319] + v[1233] * v[336])
		- v[1215] * v[363] - v[1177] * v[548];
	v[1603] = -v[1347] + gti[0] * (v[1149] * v[1377] + v[1374] * v[316] + v[1380] * v[331]) + gti[1] * (v[1150] * v[1377]
		+ v[1374] * v[318] + v[1380] * v[334]) + gti[2] * (v[1151] * v[1377] + v[1374] * v[319] + v[1380] * v[336])
		- v[1353] * v[363] - v[1177] * v[687] + v[1413] * v[704] + v[1414] * v[710];
	v[1596] = -v[1346] + gti[0] * (v[1149] * v[1376] + v[1373] * v[316] + v[1379] * v[331]) + gti[1] * (v[1150] * v[1376]
		+ v[1373] * v[318] + v[1379] * v[334]) + gti[2] * (v[1151] * v[1376] + v[1373] * v[319] + v[1379] * v[336])
		- v[1352] * v[363] - v[1177] * v[685] + v[1413] * v[703] + v[1414] * v[709];
	v[1589] = -v[1345] + gti[0] * (v[1149] * v[1375] + v[1372] * v[316] + v[1378] * v[331]) + gti[1] * (v[1150] * v[1375]
		+ v[1372] * v[318] + v[1378] * v[334]) + gti[2] * (v[1151] * v[1375] + v[1372] * v[319] + v[1378] * v[336])
		- v[1351] * v[363] - v[1177] * v[683] + v[1413] * v[702] + v[1414] * v[708];
	v[1576] = 1e0 - (v[363] * v[363]) + v[1413] * v[701] + v[1414] * v[707];
	v[1567] = v[1455] + v[1413] * v[700] + v[1414] * v[706];
	v[1559] = v[1447] + v[1413] * v[699] + v[1414] * v[705];
	v[1412] = gti[0] * (v[1149] * v[1226] + v[1224] * v[316] + v[1228] * v[331]) + gti[1] * (v[1150] * v[1226]
		+ v[1224] * v[318] + v[1228] * v[334]) + gti[2] * (v[1151] * v[1226] + v[1224] * v[319] + v[1228] * v[336])
		- v[1216] * v[361] - v[1177] * v[547];
	v[1411] = gti[0] * (v[1149] * v[1225] + v[1223] * v[316] + v[1227] * v[331]) + gti[1] * (v[1150] * v[1225]
		+ v[1223] * v[318] + v[1227] * v[334]) + gti[2] * (v[1151] * v[1225] + v[1223] * v[319] + v[1227] * v[336])
		- v[1215] * v[361] - v[1177] * v[546];
	v[1602] = -v[1327] + gti[0] * (v[1149] * v[1368] + v[1365] * v[316] + v[1371] * v[331]) + gti[1] * (v[1150] * v[1368]
		+ v[1365] * v[318] + v[1371] * v[334]) + gti[2] * (v[1151] * v[1368] + v[1365] * v[319] + v[1371] * v[336])
		- v[1353] * v[361] - v[1177] * v[681] + v[1411] * v[704] + v[1412] * v[710];
	v[1595] = -v[1326] + gti[0] * (v[1149] * v[1367] + v[1364] * v[316] + v[1370] * v[331]) + gti[1] * (v[1150] * v[1367]
		+ v[1364] * v[318] + v[1370] * v[334]) + gti[2] * (v[1151] * v[1367] + v[1364] * v[319] + v[1370] * v[336])
		- v[1352] * v[361] - v[1177] * v[679] + v[1411] * v[703] + v[1412] * v[709];
	v[1588] = -v[1325] + gti[0] * (v[1149] * v[1366] + v[1363] * v[316] + v[1369] * v[331]) + gti[1] * (v[1150] * v[1366]
		+ v[1363] * v[318] + v[1369] * v[334]) + gti[2] * (v[1151] * v[1366] + v[1363] * v[319] + v[1369] * v[336])
		- v[1351] * v[361] - v[1177] * v[677] + v[1411] * v[702] + v[1412] * v[708];
	v[1575] = v[1455] + v[1411] * v[701] + v[1412] * v[707];
	v[1565] = 1e0 - (v[361] * v[361]) + v[1411] * v[700] + v[1412] * v[706];
	v[1558] = v[1446] + v[1411] * v[699] + v[1412] * v[705];
	v[1410] = gti[0] * (v[1149] * v[1220] + v[1218] * v[316] + v[1222] * v[331]) + gti[1] * (v[1150] * v[1220]
		+ v[1218] * v[318] + v[1222] * v[334]) + gti[2] * (v[1151] * v[1220] + v[1218] * v[319] + v[1222] * v[336])
		- v[1216] * v[358] - v[1177] * v[545];
	v[1409] = gti[0] * (v[1149] * v[1219] + v[1217] * v[316] + v[1221] * v[331]) + gti[1] * (v[1150] * v[1219]
		+ v[1217] * v[318] + v[1221] * v[334]) + gti[2] * (v[1151] * v[1219] + v[1217] * v[319] + v[1221] * v[336])
		- v[1215] * v[358] - v[1177] * v[544];
	v[1601] = -v[1308] + gti[0] * (v[1149] * v[1359] + v[1356] * v[316] + v[1362] * v[331]) + gti[1] * (v[1150] * v[1359]
		+ v[1356] * v[318] + v[1362] * v[334]) + gti[2] * (v[1151] * v[1359] + v[1356] * v[319] + v[1362] * v[336])
		- v[1353] * v[358] - v[1177] * v[675] + v[1409] * v[704] + v[1410] * v[710];
	v[1594] = -v[1307] + gti[0] * (v[1149] * v[1358] + v[1355] * v[316] + v[1361] * v[331]) + gti[1] * (v[1150] * v[1358]
		+ v[1355] * v[318] + v[1361] * v[334]) + gti[2] * (v[1151] * v[1358] + v[1355] * v[319] + v[1361] * v[336])
		- v[1352] * v[358] - v[1177] * v[673] + v[1409] * v[703] + v[1410] * v[709];
	v[1587] = -v[1306] + gti[0] * (v[1149] * v[1357] + v[1354] * v[316] + v[1360] * v[331]) + gti[1] * (v[1150] * v[1357]
		+ v[1354] * v[318] + v[1360] * v[334]) + gti[2] * (v[1151] * v[1357] + v[1354] * v[319] + v[1360] * v[336])
		- v[1351] * v[358] - v[1177] * v[671] + v[1409] * v[702] + v[1410] * v[708];
	v[1574] = v[1447] + v[1409] * v[701] + v[1410] * v[707];
	v[1564] = v[1446] + v[1409] * v[700] + v[1410] * v[706];
	v[1554] = 1e0 - (v[358] * v[358]) + v[1409] * v[699] + v[1410] * v[705];
	v[1170] = v[1149] * v[1158] + v[1155] * v[316] + v[1161] * v[331];
	v[1171] = v[1149] * v[1159] + v[1156] * v[316] + v[1162] * v[331];
	v[1172] = v[1149] * v[1160] + v[1157] * v[316] + v[1163] * v[331];
	v[1173] = v[1165] + gti[0] * (v[1149] * v[1171] + v[1170] * v[316] + v[1172] * v[331]) + gti[1] * (v[1150] * v[1171]
		+ v[1170] * v[318] + v[1172] * v[334]) + gti[2] * (v[1151] * v[1171] + v[1170] * v[319] + v[1172] * v[336])
		- v[1177] * v[358];
	v[1174] = v[1150] * v[1158] + v[1155] * v[318] + v[1161] * v[334];
	v[1175] = v[1150] * v[1159] + v[1156] * v[318] + v[1162] * v[334];
	v[1176] = v[1150] * v[1160] + v[1157] * v[318] + v[1163] * v[334];
	v[1178] = v[1167] + gti[0] * (v[1149] * v[1175] + v[1174] * v[316] + v[1176] * v[331]) + gti[1] * (v[1150] * v[1175]
		+ v[1174] * v[318] + v[1176] * v[334]) + gti[2] * (v[1151] * v[1175] + v[1174] * v[319] + v[1176] * v[336])
		- v[1177] * v[361];
	v[1179] = v[1151] * v[1158] + v[1155] * v[319] + v[1161] * v[336];
	v[1180] = v[1151] * v[1159] + v[1156] * v[319] + v[1162] * v[336];
	v[1181] = v[1151] * v[1160] + v[1157] * v[319] + v[1163] * v[336];
	v[1182] = v[1169] + gti[0] * (v[1149] * v[1180] + v[1179] * v[316] + v[1181] * v[331]) + gti[1] * (v[1150] * v[1180]
		+ v[1179] * v[318] + v[1181] * v[334]) + gti[2] * (v[1151] * v[1180] + v[1179] * v[319] + v[1181] * v[336])
		- v[1177] * v[363];
	v[1415] = 1e0 / sqrt(Power(v[1173], 2) + Power(v[1178], 2) + Power(v[1182], 2));
	v[4296] = 1e0*v[1415];
	v[4295] = v[1188] * v[1415];
	v[1719] = v[1601] * v[4295];
	v[1718] = v[1602] * v[4295];
	v[1717] = v[1603] * v[4295];
	v[1707] = v[1594] * v[4295];
	v[1706] = v[1595] * v[4295];
	v[1705] = v[1596] * v[4295];
	v[1692] = v[1587] * v[4295];
	v[1691] = v[1588] * v[4295];
	v[1690] = v[1589] * v[4295];
	v[1674] = v[1449] * v[4295];
	v[1673] = v[1457] * v[4295];
	v[1672] = v[1464] * v[4295];
	v[1662] = v[1448] * v[4295];
	v[1661] = v[1456] * v[4295];
	v[1660] = v[1463] * v[4295];
	v[1650] = v[1574] * v[4295];
	v[1649] = v[1575] * v[4295];
	v[1648] = v[1576] * v[4295];
	v[1625] = v[1564] * v[4295];
	v[1623] = v[1565] * v[4295];
	v[1621] = v[1567] * v[4295];
	v[1581] = v[1554] * v[4295];
	v[1580] = v[1558] * v[4295];
	v[1579] = v[1559] * v[4295];
	v[1183] = v[1173] * v[4296];
	v[4329] = -(v[1183] * v[1554]);
	v[4327] = -(v[1183] * v[1564]);
	v[4325] = -(v[1183] * v[1574]);
	v[4323] = -(v[1183] * v[1448]);
	v[4321] = -(v[1183] * v[1449]);
	v[4319] = -(v[1183] * v[1450]);
	v[4317] = -(v[1183] * v[1587]);
	v[4315] = -(v[1183] * v[1594]);
	v[4313] = -(v[1183] * v[1601]);
	v[1569] = 1e0 - (v[1183] * v[1183]);
	v[1185] = v[1178] * v[4296];
	v[4330] = -(v[1185] * v[1558]);
	v[4328] = -(v[1185] * v[1565]);
	v[4326] = v[1185] * v[1575];
	v[4324] = v[1185] * v[1456];
	v[4322] = v[1185] * v[1457];
	v[4320] = v[1185] * v[1458];
	v[4318] = v[1185] * v[1588];
	v[4316] = v[1185] * v[1595];
	v[4314] = v[1185] * v[1602];
	v[1568] = 1e0 - (v[1185] * v[1185]);
	v[1186] = v[1182] * v[4296];
	v[4312] = -(v[1186] * v[1559]);
	v[4310] = -(v[1186] * v[1567]);
	v[4308] = -(v[1186] * v[1576]);
	v[4306] = v[1186] * v[1463];
	v[4305] = v[1186] * v[1464];
	v[4304] = v[1186] * v[1465];
	v[4303] = -(v[1186] * v[1589]);
	v[4301] = -(v[1186] * v[1596]);
	v[4299] = -(v[1186] * v[1603]);
	v[3355] = v[1186] * v[2955] + v[1185] * v[2971] + v[1183] * v[2985];
	v[3354] = v[1186] * v[2946] + v[1185] * v[2963] + v[1183] * v[2978];
	v[3353] = v[1186] * (v[1413] * v[2049] + v[1414] * v[2069] + v[2118] * v[703] + v[2116] * v[709]) + v[1185] *
		(v[1411] * v[2049] + v[1412] * v[2069] + v[2121] * v[703] + v[2119] * v[709]) + v[1183] * (v[1409] * v[2049]
			+ v[1410] * v[2069] + v[2124] * v[703] + v[2122] * v[709]);
	v[4366] = v[3353] + v[3355];
	v[3352] = v[1186] * (v[1413] * v[2050] + v[1414] * v[2072] + v[2116] * v[703] + v[2117] * v[709]) + v[1185] *
		(v[1411] * v[2050] + v[1412] * v[2072] + v[2119] * v[703] + v[2120] * v[709]) + v[1183] * (v[1409] * v[2050]
			+ v[1410] * v[2072] + v[2122] * v[703] + v[2123] * v[709]);
	v[4365] = v[3352] + v[3354];
	v[3348] = v[1186] * v[2954] + v[1185] * v[2970] + v[1183] * v[2984];
	v[3347] = v[1186] * v[2945] + v[1185] * v[2962] + v[1183] * v[2977];
	v[3346] = v[1186] * (v[1413] * v[2047] + v[1414] * v[2065] + v[2118] * v[702] + v[2116] * v[708]) + v[1185] *
		(v[1411] * v[2047] + v[1412] * v[2065] + v[2121] * v[702] + v[2119] * v[708]) + v[1183] * (v[1409] * v[2047]
			+ v[1410] * v[2065] + v[2124] * v[702] + v[2122] * v[708]);
	v[4348] = v[3346] + v[3348];
	v[3345] = v[1186] * (v[1413] * v[2048] + v[1414] * v[2068] + v[2116] * v[702] + v[2117] * v[708]) + v[1185] *
		(v[1411] * v[2048] + v[1412] * v[2068] + v[2119] * v[702] + v[2120] * v[708]) + v[1183] * (v[1409] * v[2048]
			+ v[1410] * v[2068] + v[2122] * v[702] + v[2123] * v[708]);
	v[4347] = v[3345] + v[3347];
	v[3341] = v[1186] * v[2953] + v[1185] * v[2969] + v[1183] * v[2983];
	v[3340] = v[1186] * v[2944] + v[1185] * v[2961] + v[1183] * v[2976];
	v[3335] = v[1186] * v[2952] + v[1185] * v[2968] + v[1183] * v[2982];
	v[3334] = v[1186] * v[2943] + v[1185] * v[2960] + v[1183] * v[2975];
	v[3328] = v[1186] * v[2951] + v[1185] * v[2967] + v[1183] * v[2981];
	v[3327] = v[1186] * v[2942] + v[1185] * v[2959] + v[1183] * v[2974];
	v[3319] = v[1183] * v[2948] + v[1185] * v[2949] + v[1186] * v[2950];
	v[3318] = v[1183] * v[2939] + v[1185] * v[2940] + v[1186] * v[2941];
	v[3317] = v[1186] * (v[1413] * v[2045] + v[1414] * v[2061] + v[2118] * v[701] + v[2116] * v[707]) + v[1185] *
		(v[1411] * v[2045] + v[1412] * v[2061] + v[2121] * v[701] + v[2119] * v[707]) + v[1183] * (v[1409] * v[2045]
			+ v[1410] * v[2061] + v[2124] * v[701] + v[2122] * v[707]);
	v[4345] = v[3317] + v[3319];
	v[3316] = v[1186] * (v[1413] * v[2046] + v[1414] * v[2064] + v[2116] * v[701] + v[2117] * v[707]) + v[1185] *
		(v[1411] * v[2046] + v[1412] * v[2064] + v[2119] * v[701] + v[2120] * v[707]) + v[1183] * (v[1409] * v[2046]
			+ v[1410] * v[2064] + v[2122] * v[701] + v[2123] * v[707]);
	v[4344] = v[3316] + v[3318];
	v[3304] = v[1186] * v[2949] + v[1183] * v[2965] + v[1185] * v[2966];
	v[3303] = v[1186] * v[2940] + v[1183] * v[2957] + v[1185] * v[2958];
	v[3300] = v[1186] * (v[1413] * v[2043] + v[1414] * v[2057] + v[2118] * v[700] + v[2116] * v[706]) + v[1185] *
		(v[1411] * v[2043] + v[1412] * v[2057] + v[2121] * v[700] + v[2119] * v[706]) + v[1183] * (v[1409] * v[2043]
			+ v[1410] * v[2057] + v[2124] * v[700] + v[2122] * v[706]);
	v[4339] = v[3300] + v[3304];
	v[3299] = v[1186] * (v[1413] * v[2044] + v[1414] * v[2060] + v[2116] * v[700] + v[2117] * v[706]) + v[1185] *
		(v[1411] * v[2044] + v[1412] * v[2060] + v[2119] * v[700] + v[2120] * v[706]) + v[1183] * (v[1409] * v[2044]
			+ v[1410] * v[2060] + v[2122] * v[700] + v[2123] * v[706]);
	v[4338] = v[3299] + v[3303];
	v[3285] = v[1186] * v[2948] + v[1185] * v[2965] + v[1183] * v[2980];
	v[3284] = v[1186] * v[2939] + v[1185] * v[2957] + v[1183] * v[2973];
	v[3281] = v[1186] * (v[1413] * v[2041] + v[1414] * v[2053] + v[2118] * v[699] + v[2116] * v[705]) + v[1185] *
		(v[1411] * v[2041] + v[1412] * v[2053] + v[2121] * v[699] + v[2119] * v[705]) + v[1183] * (v[1409] * v[2041]
			+ v[1410] * v[2053] + v[2124] * v[699] + v[2122] * v[705]);
	v[4332] = v[3281] + v[3285];
	v[3280] = v[1186] * (v[1413] * v[2042] + v[1414] * v[2056] + v[2116] * v[699] + v[2117] * v[705]) + v[1185] *
		(v[1411] * v[2042] + v[1412] * v[2056] + v[2119] * v[699] + v[2120] * v[705]) + v[1183] * (v[1409] * v[2042]
			+ v[1410] * v[2056] + v[2122] * v[699] + v[2123] * v[705]);
	v[4331] = v[3280] + v[3284];
	v[1641] = v[1185] * v[4298] + v[4295] * (v[1568] * v[1602] + v[1185] * (v[4299] + v[4313]));
	v[1640] = v[1183] * v[4298] + v[4295] * (v[1569] * v[1601] - v[1183] * (-v[4299] + v[4314]));
	v[1637] = v[1185] * v[4300] + v[4295] * (v[1568] * v[1595] + v[1185] * (v[4301] + v[4315]));
	v[1636] = v[1183] * v[4300] + v[4295] * (v[1569] * v[1594] - v[1183] * (-v[4301] + v[4316]));
	v[1633] = v[1185] * v[4302] + v[4295] * (v[1568] * v[1588] + v[1185] * (v[4303] + v[4317]));
	v[1632] = v[1183] * v[4302] + v[4295] * (v[1569] * v[1587] - v[1183] * (-v[4303] + v[4318]));
	v[1630] = v[1450] * v[1569] - v[1183] * (v[4304] + v[4320]);
	v[1696] = v[1630] * v[4295];
	v[1629] = v[1458] * v[1568] + v[1185] * (-v[4304] + v[4319]);
	v[1695] = v[1629] * v[4295];
	v[1626] = v[1449] * v[1569] - v[1183] * (v[4305] + v[4322]);
	v[1624] = v[1457] * v[1568] + v[1185] * (-v[4305] + v[4321]);
	v[1619] = v[1448] * v[1569] - v[1183] * (v[4306] + v[4324]);
	v[1618] = v[1456] * v[1568] + v[1185] * (-v[4306] + v[4323]);
	v[1614] = v[1185] * v[4307] + v[4295] * (v[1568] * v[1575] + v[1185] * (v[4308] + v[4325]));
	v[1613] = v[1183] * v[4307] + v[4295] * (v[1569] * v[1574] - v[1183] * (-v[4308] + v[4326]));
	v[1610] = v[1185] * v[4309] + v[4295] * (v[1565] * v[1568] + v[1185] * (v[4310] + v[4327]));
	v[1609] = v[1183] * v[4309] + v[4295] * (v[1564] * v[1569] + v[1183] * (v[4310] + v[4328]));
	v[1606] = v[1185] * v[4311] + v[4295] * (v[1558] * v[1568] + v[1185] * (v[4312] + v[4329]));
	v[1605] = v[1183] * v[4311] + v[4295] * (v[1554] * v[1569] + v[1183] * (v[4312] + v[4330]));
	v[1566] = 1e0 - (v[1186] * v[1186]);
	v[1642] = v[1186] * v[4298] + v[4295] * (v[1566] * v[1603] + v[1186] * (v[4313] - v[4314]));
	v[1638] = v[1186] * v[4300] + v[4295] * (v[1566] * v[1596] + v[1186] * (v[4315] - v[4316]));
	v[1634] = v[1186] * v[4302] + v[4295] * (v[1566] * v[1589] + v[1186] * (v[4317] - v[4318]));
	v[1628] = v[1465] * v[1566] + v[1186] * (v[4319] - v[4320]);
	v[1694] = v[1628] * v[4295];
	v[1622] = v[1464] * v[1566] + v[1186] * (v[4321] - v[4322]);
	v[1617] = v[1463] * v[1566] + v[1186] * (v[4323] - v[4324]);
	v[1615] = v[1186] * v[4307] + v[4295] * (v[1566] * v[1576] + v[1186] * (v[4325] - v[4326]));
	v[1611] = v[1186] * v[4309] + v[4295] * (v[1566] * v[1567] + v[1186] * (v[4327] + v[4328]));
	v[1607] = v[1186] * v[4311] + v[4295] * (v[1559] * v[1566] + v[1186] * (v[4329] + v[4330]));
	v[1187] = v[1183] * v[1188];
	v[1189] = v[1185] * v[1188];
	v[1190] = v[1186] * v[1188];
	v[3286] = v[1188] * (v[1186] * (v[1413] * v[2736] + v[1414] * v[2796] + v[4333]) + v[1183] * (v[1409] * v[2736]
		+ v[1410] * v[2796] + v[4362]) + v[4332] * v[700] + v[1185] * (v[1411] * v[2736] + v[1412] * v[2796] + v[2966] * v[699]
			+ v[2958] * v[705]) + v[4331] * v[706]);
	v[3291] = v[1188] * (v[1185] * (v[1411] * v[2737] + v[1412] * v[2797] + v[4333]) + v[1183] * (v[1409] * v[2737]
		+ v[1410] * v[2797] + v[4363]) + v[4332] * v[701] + v[1186] * (v[1413] * v[2737] + v[1414] * v[2797] + v[2950] * v[699]
			+ v[2941] * v[705]) + v[4331] * v[707]);
	v[3292] = v[1188] * (v[3328] * v[699] + v[3327] * v[705]);
	v[3293] = v[1188] * (v[3335] * v[699] + v[3334] * v[705]);
	v[3294] = v[1188] * (v[3341] * v[699] + v[3340] * v[705]);
	v[3295] = v[1188] * (v[4332] * v[702] + v[1186] * (v[1413] * v[2738] + v[1414] * v[2798] + v[3159] + v[2954] * v[699]
		+ v[2945] * v[705]) + v[1185] * (v[1411] * v[2738] + v[1412] * v[2798] + v[3144] + v[2970] * v[699] + v[2962] * v[705])
		+ v[1183] * (v[1409] * v[2738] + v[1410] * v[2798] + v[358] * v[4334] + v[2984] * v[699] + v[2977] * v[705])
		+ v[4331] * v[708]);
	v[3296] = v[1188] * (v[4332] * v[703] + v[1186] * (v[1413] * v[2739] + v[1414] * v[2799] + v[3160] + v[2955] * v[699]
		+ v[2946] * v[705]) + v[1185] * (v[1411] * v[2739] + v[1412] * v[2799] + v[3145] + v[2971] * v[699] + v[2963] * v[705])
		+ v[1183] * (v[1409] * v[2739] + v[1410] * v[2799] + v[358] * v[4335] + v[2985] * v[699] + v[2978] * v[705])
		+ v[4331] * v[709]);
	v[3297] = v[1188] * (v[4332] * v[704] + v[1186] * (v[1413] * v[2740] + v[1414] * v[2800] + v[3161] + v[2956] * v[699]
		+ v[2947] * v[705]) + v[1185] * (v[1411] * v[2740] + v[1412] * v[2800] + v[3146] + v[2972] * v[699] + v[2964] * v[705])
		+ v[1183] * (v[1409] * v[2740] + v[1410] * v[2800] + v[358] * v[4336] + v[2986] * v[699] + v[2979] * v[705])
		+ v[4331] * v[710]);
	v[3305] = v[1188] * (v[1185] * (v[1411] * v[2743] + v[1412] * v[2803] + v[4364]) + v[4339] * v[701] + v[1183] *
		(v[1409] * v[2743] + v[1410] * v[2803] + v[2948] * v[700] + v[2939] * v[706]) + v[1186] * (v[1413] * v[2743]
			+ v[1414] * v[2803] + v[2950] * v[700] + v[2941] * v[706]) + v[4338] * v[707]);
	v[3306] = v[1188] * (v[3328] * v[700] + v[3327] * v[706]);
	v[3307] = v[1188] * (v[3335] * v[700] + v[3334] * v[706]);
	v[3308] = v[1188] * (v[3341] * v[700] + v[3340] * v[706]);
	v[3309] = v[1188] * (v[4339] * v[702] + v[1186] * (v[1413] * v[2744] + v[1414] * v[2804] + v[3156] + v[2954] * v[700]
		+ v[2945] * v[706]) + v[1185] * (v[1411] * v[2744] + v[1412] * v[2804] + v[361] * v[4337] + v[2970] * v[700]
			+ v[2962] * v[706]) + v[1183] * (v[1409] * v[2744] + v[1410] * v[2804] + v[3144] + v[2984] * v[700] + v[2977] * v[706])
		+ v[4338] * v[708]);
	v[3310] = v[1188] * (v[4339] * v[703] + v[1186] * (v[1413] * v[2745] + v[1414] * v[2805] + v[3157] + v[2955] * v[700]
		+ v[2946] * v[706]) + v[1185] * (v[1411] * v[2745] + v[1412] * v[2805] + v[361] * v[4340] + v[2971] * v[700]
			+ v[2963] * v[706]) + v[1183] * (v[1409] * v[2745] + v[1410] * v[2805] + v[3145] + v[2985] * v[700] + v[2978] * v[706])
		+ v[4338] * v[709]);
	v[3311] = v[1188] * (v[4339] * v[704] + v[1186] * (v[1413] * v[2746] + v[1414] * v[2806] + v[3158] + v[2956] * v[700]
		+ v[2947] * v[706]) + v[1185] * (v[1411] * v[2746] + v[1412] * v[2806] + v[361] * v[4341] + v[2972] * v[700]
			+ v[2964] * v[706]) + v[1183] * (v[1409] * v[2746] + v[1410] * v[2806] + v[3146] + v[2986] * v[700] + v[2979] * v[706])
		+ v[4338] * v[710]);
	v[3313] = v[1188] * (v[3328] * v[701] + v[3327] * v[707]);
	v[3314] = v[1188] * (v[3335] * v[701] + v[3334] * v[707]);
	v[3315] = v[1188] * (v[3341] * v[701] + v[3340] * v[707]);
	v[3320] = v[1188] * (v[4345] * v[702] + v[1186] * (v[1413] * v[2750] + v[1414] * v[2810] + v[363] * v[4342]
		+ v[2954] * v[701] + v[2945] * v[707]) + v[1185] * (v[1411] * v[2750] + v[1412] * v[2810] + v[3156] + v[2970] * v[701]
			+ v[2962] * v[707]) + v[1183] * (v[1409] * v[2750] + v[1410] * v[2810] + v[3159] + v[2984] * v[701] + v[2977] * v[707])
		+ v[4344] * v[708]);
	v[3321] = v[1188] * (v[4345] * v[703] + v[1186] * (v[1413] * v[2751] + v[1414] * v[2811] + v[363] * v[4343]
		+ v[2955] * v[701] + v[2946] * v[707]) + v[1185] * (v[1411] * v[2751] + v[1412] * v[2811] + v[3157] + v[2971] * v[701]
			+ v[2963] * v[707]) + v[1183] * (v[1409] * v[2751] + v[1410] * v[2811] + v[3160] + v[2985] * v[701] + v[2978] * v[707])
		+ v[4344] * v[709]);
	v[3322] = v[1188] * (v[4345] * v[704] + v[1186] * (v[1413] * v[2752] + v[1414] * v[2812] + v[363] * v[4346]
		+ v[2956] * v[701] + v[2947] * v[707]) + v[1185] * (v[1411] * v[2752] + v[1412] * v[2812] + v[3158] + v[2972] * v[701]
			+ v[2964] * v[707]) + v[1183] * (v[1409] * v[2752] + v[1410] * v[2812] + v[3161] + v[2986] * v[701] + v[2979] * v[707])
		+ v[4344] * v[710]);
	v[3326] = v[1188] * (-(v[1183] * (v[3174] * v[358] + v[1348] * v[671])) - v[1185] * (v[3174] * v[361] + v[1348] * v[677])
		- v[1186] * (v[3174] * v[363] + v[1348] * v[683]) + v[3328] * v[702] + v[3327] * v[708]);
	v[3329] = v[1188] * (-(v[1183] * (v[3175] * v[358] + v[1348] * v[673])) - v[1185] * (v[3175] * v[361] + v[1348] * v[679])
		- v[1186] * (v[3175] * v[363] + v[1348] * v[685]) + v[3328] * v[703] + v[3327] * v[709]);
	v[3330] = v[1188] * (-(v[1183] * (v[3176] * v[358] + v[1348] * v[675])) - v[1185] * (v[3176] * v[361] + v[1348] * v[681])
		- v[1186] * (v[3176] * v[363] + v[1348] * v[687]) + v[3328] * v[704] + v[3327] * v[710]);
	v[3333] = v[1188] * (-(v[1183] * (v[3170] * v[358] + v[1349] * v[671])) - v[1185] * (v[3170] * v[361] + v[1349] * v[677])
		- v[1186] * (v[3170] * v[363] + v[1349] * v[683]) + v[3335] * v[702] + v[3334] * v[708]);
	v[3336] = v[1188] * (-(v[1183] * (v[3171] * v[358] + v[1349] * v[673])) - v[1185] * (v[3171] * v[361] + v[1349] * v[679])
		- v[1186] * (v[3171] * v[363] + v[1349] * v[685]) + v[3335] * v[703] + v[3334] * v[709]);
	v[3337] = v[1188] * (-(v[1183] * (v[3172] * v[358] + v[1349] * v[675])) - v[1185] * (v[3172] * v[361] + v[1349] * v[681])
		- v[1186] * (v[3172] * v[363] + v[1349] * v[687]) + v[3335] * v[704] + v[3334] * v[710]);
	v[3339] = v[1188] * (-(v[1183] * (v[3165] * v[358] + v[1350] * v[671])) - v[1185] * (v[3165] * v[361] + v[1350] * v[677])
		- v[1186] * (v[3165] * v[363] + v[1350] * v[683]) + v[3341] * v[702] + v[3340] * v[708]);
	v[3342] = v[1188] * (-(v[1183] * (v[3166] * v[358] + v[1350] * v[673])) - v[1185] * (v[3166] * v[361] + v[1350] * v[679])
		- v[1186] * (v[3166] * v[363] + v[1350] * v[685]) + v[3341] * v[703] + v[3340] * v[709]);
	v[3343] = v[1188] * (-(v[1183] * (v[3167] * v[358] + v[1350] * v[675])) - v[1185] * (v[3167] * v[361] + v[1350] * v[681])
		- v[1186] * (v[3167] * v[363] + v[1350] * v[687]) + v[3341] * v[704] + v[3340] * v[710]);
	v[3349] = v[1188] * (v[4348] * v[703] + v[1186] * (-(v[1083] * v[1177]) + v[1413] * v[2763] + v[1414] * v[2823] - v[3018]
		+ gti[0] * (v[1149] * v[3239] + v[316] * v[3263] + v[3215] * v[331]) + gti[1] * (v[1150] * v[3239] + v[318] * v[3263]
			+ v[3215] * v[334]) + gti[2] * (v[1151] * v[3239] + v[319] * v[3263] + v[3215] * v[336]) - v[3276] * v[363]
		- v[1352] * v[683] - v[1351] * v[685] + v[2955] * v[702] + v[2946] * v[708]) + v[1185] * (-(v[1077] * v[1177])
			+ v[1411] * v[2763] + v[1412] * v[2823] - v[3024] + gti[0] * (v[1149] * v[3241] + v[316] * v[3265] + v[3217] * v[331])
			+ gti[1] * (v[1150] * v[3241] + v[318] * v[3265] + v[3217] * v[334]) + gti[2] * (v[1151] * v[3241] + v[319] * v[3265]
				+ v[3217] * v[336]) - v[3276] * v[361] - v[1352] * v[677] - v[1351] * v[679] + v[2971] * v[702] + v[2963] * v[708])
		+ v[1183] * (-(v[1071] * v[1177]) + v[1409] * v[2763] + v[1410] * v[2823] - v[3030] + gti[0] * (v[1149] * v[3243]
			+ v[316] * v[3267] + v[3219] * v[331]) + gti[1] * (v[1150] * v[3243] + v[318] * v[3267] + v[3219] * v[334]) + gti[2] *
			(v[1151] * v[3243] + v[319] * v[3267] + v[3219] * v[336]) - v[3276] * v[358] - v[1352] * v[671] - v[1351] * v[673]
			+ v[2985] * v[702] + v[2978] * v[708]) + v[4347] * v[709]);
	v[3350] = v[1188] * (v[4348] * v[704] + v[1186] * (-(v[1085] * v[1177]) + v[1413] * v[2764] + v[1414] * v[2824] - v[3015]
		+ gti[0] * (v[1149] * v[3228] + v[316] * v[3252] + v[3204] * v[331]) + gti[1] * (v[1150] * v[3228] + v[318] * v[3252]
			+ v[3204] * v[334]) + gti[2] * (v[1151] * v[3228] + v[319] * v[3252] + v[3204] * v[336]) - v[3273] * v[363]
		- v[1353] * v[683] - v[1351] * v[687] + v[2956] * v[702] + v[2947] * v[708]) + v[1185] * (-(v[1079] * v[1177])
			+ v[1411] * v[2764] + v[1412] * v[2824] - v[3021] + gti[0] * (v[1149] * v[3231] + v[316] * v[3255] + v[3207] * v[331])
			+ gti[1] * (v[1150] * v[3231] + v[318] * v[3255] + v[3207] * v[334]) + gti[2] * (v[1151] * v[3231] + v[319] * v[3255]
				+ v[3207] * v[336]) - v[3273] * v[361] - v[1353] * v[677] - v[1351] * v[681] + v[2972] * v[702] + v[2964] * v[708])
		+ v[1183] * (-(v[1073] * v[1177]) + v[1409] * v[2764] + v[1410] * v[2824] - v[3027] + gti[0] * (v[1149] * v[3234]
			+ v[316] * v[3258] + v[3210] * v[331]) + gti[1] * (v[1150] * v[3234] + v[318] * v[3258] + v[3210] * v[334]) + gti[2] *
			(v[1151] * v[3234] + v[319] * v[3258] + v[3210] * v[336]) - v[3273] * v[358] - v[1353] * v[671] - v[1351] * v[675]
			+ v[2986] * v[702] + v[2979] * v[708]) + v[4347] * v[710]);
	v[3356] = v[1188] * (v[4366] * v[704] + v[1186] * (-(v[1097] * v[1177]) + v[1413] * v[2778] + v[1414] * v[2838] - v[3016]
		+ gti[0] * (v[1149] * v[3229] + v[316] * v[3253] + v[3205] * v[331]) + gti[1] * (v[1150] * v[3229] + v[318] * v[3253]
			+ v[3205] * v[334]) + gti[2] * (v[1151] * v[3229] + v[319] * v[3253] + v[3205] * v[336]) - v[3274] * v[363]
		- v[1353] * v[685] - v[1352] * v[687] + v[2956] * v[703] + v[2947] * v[709]) + v[1185] * (-(v[1093] * v[1177])
			+ v[1411] * v[2778] + v[1412] * v[2838] - v[3022] + gti[0] * (v[1149] * v[3232] + v[316] * v[3256] + v[3208] * v[331])
			+ gti[1] * (v[1150] * v[3232] + v[318] * v[3256] + v[3208] * v[334]) + gti[2] * (v[1151] * v[3232] + v[319] * v[3256]
				+ v[3208] * v[336]) - v[3274] * v[361] - v[1353] * v[679] - v[1352] * v[681] + v[2972] * v[703] + v[2964] * v[709])
		+ v[1183] * (-(v[1089] * v[1177]) + v[1409] * v[2778] + v[1410] * v[2838] - v[3028] + gti[0] * (v[1149] * v[3235]
			+ v[316] * v[3259] + v[3211] * v[331]) + gti[1] * (v[1150] * v[3235] + v[318] * v[3259] + v[3211] * v[334]) + gti[2] *
			(v[1151] * v[3235] + v[319] * v[3259] + v[3211] * v[336]) - v[3274] * v[358] - v[1353] * v[673] - v[1352] * v[675]
			+ v[2986] * v[703] + v[2979] * v[709]) + v[4365] * v[710]);
	v[3370] = v[257] * v[4352];
	v[3377] = v[253] * v[4354];
	v[3366] = v[3365] - v[903];
	v[4353] = v[3366] / v[3365];
	v[3676] = v[2605] / v[3365];
	v[3658] = v[3676] * v[4049];
	v[3654] = v[2604] / v[3365];
	v[3684] = v[3654] * v[4048];
	v[3636] = v[2603] / v[3365];
	v[3683] = v[3636] * v[4047];
	v[3635] = v[3377] + v[258] * v[4352] + v[248] * v[4353];
	v[3637] = (*a4)*v[3635] + v[3677] + v[3636] * v[4049];
	v[3632] = v[3370] + v[247] * v[4353] + v[252] * v[4354];
	v[3634] = (*a4)*v[3632] + v[3655] + v[3636] * v[4048];
	v[3631] = v[3658] + v[3684] + v[4026] * v[4047] + (*a4)*v[4353];
	v[3367] = v[3632] * v[4048] + v[3635] * v[4049] + v[4047] * v[4353];
	v[3369] = v[3365] - v[900];
	v[4355] = v[3369] / v[3365];
	v[3656] = -v[3655] + v[3654] * v[4047] + (*a4)*v[4354];
	v[3380] = -((v[248] * v[4041]) / v[3365]);
	v[3660] = v[3380] - (v[258] * v[4038]) / v[3365] + v[253] * v[4355];
	v[3661] = (*a4)*v[3660] + v[3680] + v[3654] * v[4049];
	v[3659] = v[3658] + v[3683] + v[4030] * v[4048] + (*a4)*v[4355];
	v[3375] = v[3660] * v[4049] + v[4047] * v[4354] + v[4048] * v[4355];
	v[3376] = v[3361] + v[3365];
	v[3678] = -v[3677] + v[3676] * v[4047] + (*a4)*v[4352];
	v[3682] = v[3376] / v[3365];
	v[3685] = (*a4)*v[3682] + v[3683] + v[3684] + v[4033] * v[4049];
	v[3679] = v[3380] + v[257] * v[3682] - (v[252] * v[4042]) / v[3365];
	v[3681] = (*a4)*v[3679] - v[3680] + v[3676] * v[4048];
	v[3381] = v[3679] * v[4048] + v[3682] * v[4049] + v[4047] * v[4352];
	v[3876] = v[3367] * v[3854] + v[3375] * v[3860] + v[3381] * v[3866] + v[3836] * v[4359] + v[3842] * v[4360]
		+ v[3848] * v[4361];
	v[3875] = v[3367] * v[3853] + v[3375] * v[3859] + v[3381] * v[3865] + v[3835] * v[4359] + v[3841] * v[4360]
		+ v[3847] * v[4361];
	v[3874] = v[3367] * v[3856] + v[3375] * v[3862] + v[3381] * v[3868] + v[3838] * v[4359] + v[3844] * v[4360]
		+ v[3850] * v[4361];
	v[3873] = v[3367] * v[3855] + v[3375] * v[3861] + v[3381] * v[3867] + v[3837] * v[4359] + v[3843] * v[4360]
		+ v[3849] * v[4361];
	v[3872] = v[3367] * v[3858] + v[3375] * v[3864] + v[3381] * v[3870] + v[3840] * v[4359] + v[3846] * v[4360]
		+ v[3852] * v[4361];
	v[3871] = v[3367] * v[3857] + v[3375] * v[3863] + v[3381] * v[3869] + v[3839] * v[4359] + v[3845] * v[4360]
		+ v[3851] * v[4361];
	v[3703] = v[3367] * v[3544] + v[3375] * v[3570] + v[3381] * v[3593] + v[1586] * v[3637] + v[1593] * v[3661]
		+ v[1600] * v[3685] + v[3487] * v[4359] + v[3505] * v[4360] + v[3523] * v[4361];
	v[3702] = v[3367] * v[3542] + v[3375] * v[3568] + v[3381] * v[3591] + v[1586] * v[3634] + v[1593] * v[3659]
		+ v[1600] * v[3681] + v[3486] * v[4359] + v[3504] * v[4360] + v[3522] * v[4361];
	v[3701] = v[3367] * v[3540] + v[3375] * v[3566] + v[3381] * v[3590] + v[1586] * v[3631] + v[1593] * v[3656]
		+ v[1600] * v[3678] + v[3485] * v[4359] + v[3503] * v[4360] + v[3521] * v[4361];
	v[3700] = (*a4)*v[1573] + v[3367] * v[3538] + v[3375] * v[3565] + v[3381] * v[3589] + v[3484] * v[4359]
		+ v[3502] * v[4360] + v[3520] * v[4361];
	v[3699] = (*a4)*v[1563] + v[3367] * v[3537] + v[3375] * v[3564] + v[3381] * v[3588] + v[3483] * v[4359]
		+ v[3501] * v[4360] + v[3519] * v[4361];
	v[3698] = (*a4)*v[1557] + v[3367] * v[3536] + v[3375] * v[3563] + v[3381] * v[3587] + v[3482] * v[4359]
		+ v[3500] * v[4360] + v[3518] * v[4361];
	v[3697] = v[3367] * v[3553] + v[3375] * v[3578] + v[3381] * v[3600] + v[1585] * v[3637] + v[1592] * v[3661]
		+ v[1599] * v[3685] + v[3493] * v[4359] + v[3511] * v[4360] + v[3529] * v[4361];
	v[3696] = v[3367] * v[3551] + v[3375] * v[3576] + v[3381] * v[3598] + v[1585] * v[3634] + v[1592] * v[3659]
		+ v[1599] * v[3681] + v[3492] * v[4359] + v[3510] * v[4360] + v[3528] * v[4361];
	v[3695] = v[3367] * v[3549] + v[3375] * v[3574] + v[3381] * v[3597] + v[1585] * v[3631] + v[1592] * v[3656]
		+ v[1599] * v[3678] + v[3491] * v[4359] + v[3509] * v[4360] + v[3527] * v[4361];
	v[3694] = (*a4)*v[1572] + v[3367] * v[3547] + v[3375] * v[3573] + v[3381] * v[3596] + v[3490] * v[4359]
		+ v[3508] * v[4360] + v[3526] * v[4361];
	v[3693] = (*a4)*v[1562] + v[3367] * v[3546] + v[3375] * v[3572] + v[3381] * v[3595] + v[3489] * v[4359]
		+ v[3507] * v[4360] + v[3525] * v[4361];
	v[3692] = (*a4)*v[1556] + v[3367] * v[3545] + v[3375] * v[3571] + v[3381] * v[3594] + v[3488] * v[4359]
		+ v[3506] * v[4360] + v[3524] * v[4361];
	v[3691] = v[3367] * v[3562] + v[3375] * v[3586] + v[3381] * v[3607] + v[1584] * v[3637] + v[1591] * v[3661]
		+ v[1598] * v[3685] + v[3499] * v[4359] + v[3517] * v[4360] + v[3535] * v[4361];
	v[3690] = v[3367] * v[3560] + v[3375] * v[3584] + v[3381] * v[3605] + v[1584] * v[3634] + v[1591] * v[3659]
		+ v[1598] * v[3681] + v[3498] * v[4359] + v[3516] * v[4360] + v[3534] * v[4361];
	v[3689] = v[3367] * v[3558] + v[3375] * v[3582] + v[3381] * v[3604] + v[1584] * v[3631] + v[1591] * v[3656]
		+ v[1598] * v[3678] + v[3497] * v[4359] + v[3515] * v[4360] + v[3533] * v[4361];
	v[3688] = (*a4)*v[1571] + v[3367] * v[3556] + v[3375] * v[3581] + v[3381] * v[3603] + v[3496] * v[4359]
		+ v[3514] * v[4360] + v[3532] * v[4361];
	v[3687] = (*a4)*v[1561] + v[3367] * v[3555] + v[3375] * v[3580] + v[3381] * v[3602] + v[3495] * v[4359]
		+ v[3513] * v[4360] + v[3531] * v[4361];
	v[3686] = (*a4)*v[1555] + v[3367] * v[3554] + v[3375] * v[3579] + v[3381] * v[3601] + v[3494] * v[4359]
		+ v[3512] * v[4360] + v[3530] * v[4361];
	v[3406] = v[1584] * v[3367] + v[1591] * v[3375] + v[1598] * v[3381] + v[1555] * v[4359] + v[1561] * v[4360]
		+ v[1571] * v[4361];
	v[3407] = v[1585] * v[3367] + v[1592] * v[3375] + v[1599] * v[3381] + v[1556] * v[4359] + v[1562] * v[4360]
		+ v[1572] * v[4361];
	v[3408] = v[1586] * v[3367] + v[1593] * v[3375] + v[1600] * v[3381] + v[1557] * v[4359] + v[1563] * v[4360]
		+ v[1573] * v[4361];
	v[3882] = (*cn)*(v[3408] * v[3848] + v[3407] * v[3850] + v[3406] * v[3852] + v[1571] * v[3872] + v[1572] * v[3874]
		+ v[1573] * v[3876]);
	v[3881] = (*cn)*(v[3408] * v[3847] + v[3407] * v[3849] + v[3406] * v[3851] + v[1571] * v[3871] + v[1572] * v[3873]
		+ v[1573] * v[3875]);
	v[3880] = (*cn)*(v[3408] * v[3842] + v[3407] * v[3844] + v[3406] * v[3846] + v[1561] * v[3872] + v[1562] * v[3874]
		+ v[1563] * v[3876]);
	v[3879] = (*cn)*(v[3408] * v[3841] + v[3407] * v[3843] + v[3406] * v[3845] + v[1561] * v[3871] + v[1562] * v[3873]
		+ v[1563] * v[3875]);
	v[3878] = (*cn)*(v[3408] * v[3836] + v[3407] * v[3838] + v[3406] * v[3840] + v[1555] * v[3872] + v[1556] * v[3874]
		+ v[1557] * v[3876]);
	v[3877] = (*cn)*(v[3408] * v[3835] + v[3407] * v[3837] + v[3406] * v[3839] + v[1555] * v[3871] + v[1556] * v[3873]
		+ v[1557] * v[3875]);
	v[3883] = (*cn)*(v[3408] * v[3853] + v[3407] * v[3855] + v[3406] * v[3857] + v[1584] * v[3871] + v[1585] * v[3873]
		+ v[1586] * v[3875]);
	v[3884] = (*cn)*(v[3408] * v[3859] + v[3407] * v[3861] + v[3406] * v[3863] + v[1591] * v[3871] + v[1592] * v[3873]
		+ v[1593] * v[3875]);
	v[3885] = (*cn)*(v[3408] * v[3865] + v[3407] * v[3867] + v[3406] * v[3869] + v[1598] * v[3871] + v[1599] * v[3873]
		+ v[1600] * v[3875]);
	v[3886] = (*cn)*(v[3408] * v[3854] + v[3407] * v[3856] + v[3406] * v[3858] + v[1584] * v[3872] + v[1585] * v[3874]
		+ v[1586] * v[3876]);
	v[3887] = (*cn)*(v[3408] * v[3860] + v[3407] * v[3862] + v[3406] * v[3864] + v[1591] * v[3872] + v[1592] * v[3874]
		+ v[1593] * v[3876]);
	v[3888] = (*cn)*(v[3408] * v[3866] + v[3407] * v[3868] + v[3406] * v[3870] + v[1598] * v[3872] + v[1599] * v[3874]
		+ v[1600] * v[3876]);
	v[3925] = v[1187] * v[1554] + v[1189] * v[1558] + v[1190] * v[1559] + (*cn)*(v[1555] * v[3406] + v[1556] * v[3407]
		+ v[1557] * v[3408]) + (*epsn)*v[364];
	v[3926] = v[1187] * v[1564] + v[1189] * v[1565] + v[1190] * v[1567] + (*cn)*(v[1561] * v[3406] + v[1562] * v[3407]
		+ v[1563] * v[3408]) + (*epsn)*v[365];
	v[3927] = v[1187] * v[1574] + v[1189] * v[1575] + v[1190] * v[1576] + (*cn)*(v[1571] * v[3406] + v[1572] * v[3407]
		+ v[1573] * v[3408]) + (*epsn)*v[366];
	v[3931] = (*epsn)*v[1555] + v[1554] * v[1605] + v[1558] * v[1606] + v[1559] * v[1607] + (*cn)*(v[3408] * v[3482]
		+ v[3407] * v[3488] + v[3406] * v[3494] + v[1555] * v[3686] + v[1556] * v[3692] + v[1557] * v[3698]) + v[3877] * v[699]
		+ v[3878] * v[705] + v[1188] * (v[1185] * (v[1411] * v[2735] + v[1412] * v[2795] + v[4362]) + v[1186] * (v[1413] * v[2735]
			+ v[1414] * v[2795] + v[4363]) + v[4332] * v[699] + v[4331] * v[705] + v[1183] * (v[1409] * v[2735] + v[1410] * v[2795]
				+ v[2980] * v[699] + v[2973] * v[705]));
	v[3932] = (*epsn)*v[1561] + v[1554] * v[1609] + v[1558] * v[1610] + v[1559] * v[1611] + v[3286] + (*cn)*
		(v[3408] * v[3483] + v[3407] * v[3489] + v[3406] * v[3495] + v[1555] * v[3687] + v[1556] * v[3693] + v[1557] * v[3699])
		+ v[3877] * v[700] + v[3878] * v[706];
	v[3933] = (*epsn)*v[1571] + v[1554] * v[1613] + v[1558] * v[1614] + v[1559] * v[1615] + v[3291] + (*cn)*
		(v[3408] * v[3484] + v[3407] * v[3490] + v[3406] * v[3496] + v[1555] * v[3688] + v[1556] * v[3694] + v[1557] * v[3700])
		+ v[3877] * v[701] + v[3878] * v[707];
	v[3934] = v[1579] * v[1617] + v[1580] * v[1618] + v[1581] * v[1619] + v[3292];
	v[3935] = v[1579] * v[1622] + v[1580] * v[1624] + v[1581] * v[1626] + v[3293];
	v[3936] = v[1579] * v[1628] + v[1580] * v[1629] + v[1581] * v[1630] + v[3294];
	v[3937] = (*epsn)*v[1584] + v[1554] * v[1632] + v[1558] * v[1633] + v[1559] * v[1634] + v[3295] + (*cn)*
		(v[3408] * v[3485] + v[3407] * v[3491] + v[3406] * v[3497] + v[1555] * v[3689] + v[1556] * v[3695] + v[1557] * v[3701])
		+ v[3877] * v[702] + v[3878] * v[708];
	v[3938] = (*epsn)*v[1591] + v[1554] * v[1636] + v[1558] * v[1637] + v[1559] * v[1638] + v[3296] + (*cn)*
		(v[3408] * v[3486] + v[3407] * v[3492] + v[3406] * v[3498] + v[1555] * v[3690] + v[1556] * v[3696] + v[1557] * v[3702])
		+ v[3877] * v[703] + v[3878] * v[709];
	v[3939] = (*epsn)*v[1598] + v[1554] * v[1640] + v[1558] * v[1641] + v[1559] * v[1642] + v[3297] + (*cn)*
		(v[3408] * v[3487] + v[3407] * v[3493] + v[3406] * v[3499] + v[1555] * v[3691] + v[1556] * v[3697] + v[1557] * v[3703])
		+ v[3877] * v[704] + v[3878] * v[710];
	v[3940] = (*epsn)*v[1556] + v[1564] * v[1605] + v[1565] * v[1606] + v[1567] * v[1607] + v[3286] + (*cn)*
		(v[3408] * v[3500] + v[3407] * v[3506] + v[3406] * v[3512] + v[1561] * v[3686] + v[1562] * v[3692] + v[1563] * v[3698])
		+ v[3879] * v[699] + v[3880] * v[705];
	v[3941] = (*epsn)*v[1562] + v[1564] * v[1609] + v[1565] * v[1610] + v[1567] * v[1611] + (*cn)*(v[3408] * v[3501]
		+ v[3407] * v[3507] + v[3406] * v[3513] + v[1561] * v[3687] + v[1562] * v[3693] + v[1563] * v[3699]) + v[3879] * v[700]
		+ v[3880] * v[706] + v[1188] * (v[1186] * (v[1413] * v[2742] + v[1414] * v[2802] + v[4364]) + v[4339] * v[700]
			+ v[4338] * v[706] + v[1183] * (v[1409] * v[2742] + v[1410] * v[2802] + v[2965] * v[700] + v[2957] * v[706]) + v[1185] *
			(v[1411] * v[2742] + v[1412] * v[2802] + v[2966] * v[700] + v[2958] * v[706]));
	v[3942] = (*epsn)*v[1572] + v[1564] * v[1613] + v[1565] * v[1614] + v[1567] * v[1615] + v[3305] + (*cn)*
		(v[3408] * v[3502] + v[3407] * v[3508] + v[3406] * v[3514] + v[1561] * v[3688] + v[1562] * v[3694] + v[1563] * v[3700])
		+ v[3879] * v[701] + v[3880] * v[707];
	v[3943] = v[1617] * v[1621] + v[1618] * v[1623] + v[1619] * v[1625] + v[3306];
	v[3944] = v[1621] * v[1622] + v[1623] * v[1624] + v[1625] * v[1626] + v[3307];
	v[3945] = v[1621] * v[1628] + v[1623] * v[1629] + v[1625] * v[1630] + v[3308];
	v[3946] = (*epsn)*v[1585] + v[1564] * v[1632] + v[1565] * v[1633] + v[1567] * v[1634] + v[3309] + (*cn)*
		(v[3408] * v[3503] + v[3407] * v[3509] + v[3406] * v[3515] + v[1561] * v[3689] + v[1562] * v[3695] + v[1563] * v[3701])
		+ v[3879] * v[702] + v[3880] * v[708];
	v[3947] = (*epsn)*v[1592] + v[1564] * v[1636] + v[1565] * v[1637] + v[1567] * v[1638] + v[3310] + (*cn)*
		(v[3408] * v[3504] + v[3407] * v[3510] + v[3406] * v[3516] + v[1561] * v[3690] + v[1562] * v[3696] + v[1563] * v[3702])
		+ v[3879] * v[703] + v[3880] * v[709];
	v[3948] = (*epsn)*v[1599] + v[1564] * v[1640] + v[1565] * v[1641] + v[1567] * v[1642] + v[3311] + (*cn)*
		(v[3408] * v[3505] + v[3407] * v[3511] + v[3406] * v[3517] + v[1561] * v[3691] + v[1562] * v[3697] + v[1563] * v[3703])
		+ v[3879] * v[704] + v[3880] * v[710];
	v[3949] = (*epsn)*v[1557] + v[1574] * v[1605] + v[1575] * v[1606] + v[1576] * v[1607] + v[3291] + (*cn)*
		(v[3408] * v[3518] + v[3407] * v[3524] + v[3406] * v[3530] + v[1571] * v[3686] + v[1572] * v[3692] + v[1573] * v[3698])
		+ v[3881] * v[699] + v[3882] * v[705];
	v[3950] = (*epsn)*v[1563] + v[1574] * v[1609] + v[1575] * v[1610] + v[1576] * v[1611] + v[3305] + (*cn)*
		(v[3408] * v[3519] + v[3407] * v[3525] + v[3406] * v[3531] + v[1571] * v[3687] + v[1572] * v[3693] + v[1573] * v[3699])
		+ v[3881] * v[700] + v[3882] * v[706];
	v[3951] = (*epsn)*v[1573] + v[1574] * v[1613] + v[1575] * v[1614] + v[1576] * v[1615] + (*cn)*(v[3408] * v[3520]
		+ v[3407] * v[3526] + v[3406] * v[3532] + v[1571] * v[3688] + v[1572] * v[3694] + v[1573] * v[3700]) + v[3881] * v[701]
		+ v[3882] * v[707] + v[1188] * (v[4345] * v[701] + v[4344] * v[707] + v[1183] * (v[1409] * v[2749] + v[1410] * v[2809]
			+ v[2948] * v[701] + v[2939] * v[707]) + v[1185] * (v[1411] * v[2749] + v[1412] * v[2809] + v[2949] * v[701]
				+ v[2940] * v[707]) + v[1186] * (v[1413] * v[2749] + v[1414] * v[2809] + v[2950] * v[701] + v[2941] * v[707]));
	v[3952] = v[1617] * v[1648] + v[1618] * v[1649] + v[1619] * v[1650] + v[3313];
	v[3953] = v[1622] * v[1648] + v[1624] * v[1649] + v[1626] * v[1650] + v[3314];
	v[3954] = v[1628] * v[1648] + v[1629] * v[1649] + v[1630] * v[1650] + v[3315];
	v[3955] = (*epsn)*v[1586] + v[1574] * v[1632] + v[1575] * v[1633] + v[1576] * v[1634] + v[3320] + (*cn)*
		(v[3408] * v[3521] + v[3407] * v[3527] + v[3406] * v[3533] + v[1571] * v[3689] + v[1572] * v[3695] + v[1573] * v[3701])
		+ v[3881] * v[702] + v[3882] * v[708];
	v[3956] = (*epsn)*v[1593] + v[1574] * v[1636] + v[1575] * v[1637] + v[1576] * v[1638] + v[3321] + (*cn)*
		(v[3408] * v[3522] + v[3407] * v[3528] + v[3406] * v[3534] + v[1571] * v[3690] + v[1572] * v[3696] + v[1573] * v[3702])
		+ v[3881] * v[703] + v[3882] * v[709];
	v[3957] = (*epsn)*v[1600] + v[1574] * v[1640] + v[1575] * v[1641] + v[1576] * v[1642] + v[3322] + (*cn)*
		(v[3408] * v[3523] + v[3407] * v[3529] + v[3406] * v[3535] + v[1571] * v[3691] + v[1572] * v[3697] + v[1573] * v[3703])
		+ v[3881] * v[704] + v[3882] * v[710];
	v[3958] = v[1448] * v[1605] + v[1456] * v[1606] + v[1463] * v[1607] + v[3292];
	v[3959] = v[1448] * v[1609] + v[1456] * v[1610] + v[1463] * v[1611] + v[3306];
	v[3960] = v[1448] * v[1613] + v[1456] * v[1614] + v[1463] * v[1615] + v[3313];
	v[3962] = v[1622] * v[1660] + v[1624] * v[1661] + v[1626] * v[1662] + v[1188] * (v[1183] * (v[3104] - v[3168] * v[358])
		+ v[1185] * (v[3078] - v[3168] * v[361]) + v[1186] * (v[3049] - v[3168] * v[363]));
	v[3963] = v[1628] * v[1660] + v[1629] * v[1661] + v[1630] * v[1662] + v[1188] * (v[1183] * (v[3091] - v[3162] * v[358])
		+ v[1185] * (v[3063] - v[3162] * v[361]) + v[1186] * (v[3037] - v[3162] * v[363]));
	v[3967] = v[1449] * v[1605] + v[1457] * v[1606] + v[1464] * v[1607] + v[3293];
	v[3968] = v[1449] * v[1609] + v[1457] * v[1610] + v[1464] * v[1611] + v[3307];
	v[3969] = v[1449] * v[1613] + v[1457] * v[1614] + v[1464] * v[1615] + v[3314];
	v[3971] = v[1628] * v[1672] + v[1629] * v[1673] + v[1630] * v[1674] + v[1188] * (v[1183] * (v[3095] - v[3163] * v[358])
		+ v[1185] * (v[3067] - v[3163] * v[361]) + v[1186] * (v[3041] - v[3163] * v[363]));
	v[3975] = v[1450] * v[1605] + v[1458] * v[1606] + v[1465] * v[1607] + v[3294];
	v[3976] = v[1450] * v[1609] + v[1458] * v[1610] + v[1465] * v[1611] + v[3308];
	v[3977] = v[1450] * v[1613] + v[1458] * v[1614] + v[1465] * v[1615] + v[3315];
	v[3982] = v[1587] * v[1605] + v[1588] * v[1606] + v[1589] * v[1607] + v[3295] + (*cn)*(v[3408] * v[3536]
		+ v[3407] * v[3545] + v[3406] * v[3554] + v[1584] * v[3686] + v[1585] * v[3692] + v[1586] * v[3698]) + v[3883] * v[699]
		+ v[3886] * v[705] + (*epsn)*(v[1059] + v[1123] * v[699] + v[1124] * v[705]);
	v[3983] = v[1587] * v[1609] + v[1588] * v[1610] + v[1589] * v[1611] + v[3309] + (*cn)*(v[3408] * v[3537]
		+ v[3407] * v[3546] + v[3406] * v[3555] + v[1584] * v[3687] + v[1585] * v[3693] + v[1586] * v[3699]) + v[3883] * v[700]
		+ v[3886] * v[706] + (*epsn)*(v[1062] + v[1123] * v[700] + v[1124] * v[706]);
	v[3984] = v[1587] * v[1613] + v[1588] * v[1614] + v[1589] * v[1615] + v[3320] + (*cn)*(v[3408] * v[3538]
		+ v[3407] * v[3547] + v[3406] * v[3556] + v[1584] * v[3688] + v[1585] * v[3694] + v[1586] * v[3700]) + v[3883] * v[701]
		+ v[3886] * v[707] + (*epsn)*(v[1065] + v[1123] * v[701] + v[1124] * v[707]);
	v[3991] = v[1594] * v[1605] + v[1595] * v[1606] + v[1596] * v[1607] + v[3296] + (*cn)*(v[3408] * v[3563]
		+ v[3407] * v[3571] + v[3406] * v[3579] + v[1591] * v[3686] + v[1592] * v[3692] + v[1593] * v[3698]) + v[3884] * v[699]
		+ v[3887] * v[705] + (*epsn)*(v[1060] + v[1131] * v[699] + v[1132] * v[705]);
	v[3992] = v[1594] * v[1609] + v[1595] * v[1610] + v[1596] * v[1611] + v[3310] + (*cn)*(v[3408] * v[3564]
		+ v[3407] * v[3572] + v[3406] * v[3580] + v[1591] * v[3687] + v[1592] * v[3693] + v[1593] * v[3699]) + v[3884] * v[700]
		+ v[3887] * v[706] + (*epsn)*(v[1063] + v[1131] * v[700] + v[1132] * v[706]);
	v[3993] = v[1594] * v[1613] + v[1595] * v[1614] + v[1596] * v[1615] + v[3321] + (*cn)*(v[3408] * v[3565]
		+ v[3407] * v[3573] + v[3406] * v[3581] + v[1591] * v[3688] + v[1592] * v[3694] + v[1593] * v[3700]) + v[3884] * v[701]
		+ v[3887] * v[707] + (*epsn)*(v[1066] + v[1131] * v[701] + v[1132] * v[707]);
	v[4000] = v[1601] * v[1605] + v[1602] * v[1606] + v[1603] * v[1607] + v[3297] + (*cn)*(v[3408] * v[3587]
		+ v[3407] * v[3594] + v[3406] * v[3601] + v[1598] * v[3686] + v[1599] * v[3692] + v[1600] * v[3698]) + v[3885] * v[699]
		+ v[3888] * v[705] + (*epsn)*(v[1061] + v[1140] * v[699] + v[1141] * v[705]);
	v[4001] = v[1601] * v[1609] + v[1602] * v[1610] + v[1603] * v[1611] + v[3311] + (*cn)*(v[3408] * v[3588]
		+ v[3407] * v[3595] + v[3406] * v[3602] + v[1598] * v[3687] + v[1599] * v[3693] + v[1600] * v[3699]) + v[3885] * v[700]
		+ v[3888] * v[706] + (*epsn)*(v[1064] + v[1140] * v[700] + v[1141] * v[706]);
	v[4002] = v[1601] * v[1613] + v[1602] * v[1614] + v[1603] * v[1615] + v[3322] + (*cn)*(v[3408] * v[3589]
		+ v[3407] * v[3596] + v[3406] * v[3603] + v[1598] * v[3688] + v[1599] * v[3694] + v[1600] * v[3700]) + v[3885] * v[701]
		+ v[3888] * v[707] + (*epsn)*(v[1067] + v[1140] * v[701] + v[1141] * v[707]);
	Rc[0] = v[3925];
	Rc[1] = v[3926];
	Rc[2] = v[3927];
	Rc[3] = v[1187] * v[1448] + v[1189] * v[1456] + v[1190] * v[1463];
	Rc[4] = v[1187] * v[1449] + v[1189] * v[1457] + v[1190] * v[1464];
	Rc[5] = v[1187] * v[1450] + v[1189] * v[1458] + v[1190] * v[1465];
	Rc[6] = -v[3925];
	Rc[7] = -v[3926];
	Rc[8] = -v[3927];
	Rc[9] = v[1187] * v[1587] + v[1189] * v[1588] + v[1190] * v[1589] + (*cn)*(v[1584] * v[3406] + v[1585] * v[3407]
		+ v[1586] * v[3408]) + (*epsn)*(v[1059] * v[364] + v[1062] * v[365] + v[1065] * v[366]);
	Rc[10] = v[1187] * v[1594] + v[1189] * v[1595] + v[1190] * v[1596] + (*cn)*(v[1591] * v[3406] + v[1592] * v[3407]
		+ v[1593] * v[3408]) + (*epsn)*(v[1060] * v[364] + v[1063] * v[365] + v[1066] * v[366]);
	Rc[11] = v[1187] * v[1601] + v[1189] * v[1602] + v[1190] * v[1603] + (*cn)*(v[1598] * v[3406] + v[1599] * v[3407]
		+ v[1600] * v[3408]) + (*epsn)*(v[1061] * v[364] + v[1064] * v[365] + v[1067] * v[366]);
	Kc[0][0] = v[3931];
	Kc[0][1] = v[3932];
	Kc[0][2] = v[3933];
	Kc[0][3] = v[3934];
	Kc[0][4] = v[3935];
	Kc[0][5] = v[3936];
	Kc[0][6] = -v[3931];
	Kc[0][7] = -v[3932];
	Kc[0][8] = -v[3933];
	Kc[0][9] = v[3937];
	Kc[0][10] = v[3938];
	Kc[0][11] = v[3939];
	Kc[1][0] = v[3940];
	Kc[1][1] = v[3941];
	Kc[1][2] = v[3942];
	Kc[1][3] = v[3943];
	Kc[1][4] = v[3944];
	Kc[1][5] = v[3945];
	Kc[1][6] = -v[3940];
	Kc[1][7] = -v[3941];
	Kc[1][8] = -v[3942];
	Kc[1][9] = v[3946];
	Kc[1][10] = v[3947];
	Kc[1][11] = v[3948];
	Kc[2][0] = v[3949];
	Kc[2][1] = v[3950];
	Kc[2][2] = v[3951];
	Kc[2][3] = v[3952];
	Kc[2][4] = v[3953];
	Kc[2][5] = v[3954];
	Kc[2][6] = -v[3949];
	Kc[2][7] = -v[3950];
	Kc[2][8] = -v[3951];
	Kc[2][9] = v[3955];
	Kc[2][10] = v[3956];
	Kc[2][11] = v[3957];
	Kc[3][0] = v[3958];
	Kc[3][1] = v[3959];
	Kc[3][2] = v[3960];
	Kc[3][3] = v[1617] * v[1660] + v[1618] * v[1661] + v[1619] * v[1662] + v[1188] * (v[1183] * (v[3113] - v[3173] * v[358])
		+ v[1185] * (v[3087] - v[3173] * v[361]) + v[1186] * (v[3059] - v[3173] * v[363]));
	Kc[3][4] = v[3962];
	Kc[3][5] = v[3963];
	Kc[3][6] = -v[3958];
	Kc[3][7] = -v[3959];
	Kc[3][8] = -v[3960];
	Kc[3][9] = v[1448] * v[1632] + v[1456] * v[1633] + v[1463] * v[1634] + v[3326];
	Kc[3][10] = v[1448] * v[1636] + v[1456] * v[1637] + v[1463] * v[1638] + v[3329];
	Kc[3][11] = v[1448] * v[1640] + v[1456] * v[1641] + v[1463] * v[1642] + v[3330];
	Kc[4][0] = v[3967];
	Kc[4][1] = v[3968];
	Kc[4][2] = v[3969];
	Kc[4][3] = v[3962];
	Kc[4][4] = v[1622] * v[1672] + v[1624] * v[1673] + v[1626] * v[1674] + v[1188] * (v[1183] * (v[3109] - v[3169] * v[358])
		+ v[1185] * (v[3082] - v[3169] * v[361]) + v[1186] * (v[3054] - v[3169] * v[363]));
	Kc[4][5] = v[3971];
	Kc[4][6] = -v[3967];
	Kc[4][7] = -v[3968];
	Kc[4][8] = -v[3969];
	Kc[4][9] = v[1449] * v[1632] + v[1457] * v[1633] + v[1464] * v[1634] + v[3333];
	Kc[4][10] = v[1449] * v[1636] + v[1457] * v[1637] + v[1464] * v[1638] + v[3336];
	Kc[4][11] = v[1449] * v[1640] + v[1457] * v[1641] + v[1464] * v[1642] + v[3337];
	Kc[5][0] = v[3975];
	Kc[5][1] = v[3976];
	Kc[5][2] = v[3977];
	Kc[5][3] = v[3963];
	Kc[5][4] = v[3971];
	Kc[5][5] = v[1465] * v[1694] + v[1458] * v[1695] + v[1450] * v[1696] + v[1188] * (v[1183] * (v[3100] - v[3164] * v[358])
		+ v[1185] * (v[3072] - v[3164] * v[361]) + v[1186] * (v[3045] - v[3164] * v[363]));
	Kc[5][6] = -v[3975];
	Kc[5][7] = -v[3976];
	Kc[5][8] = -v[3977];
	Kc[5][9] = v[1450] * v[1632] + v[1458] * v[1633] + v[1465] * v[1634] + v[3339];
	Kc[5][10] = v[1450] * v[1636] + v[1458] * v[1637] + v[1465] * v[1638] + v[3342];
	Kc[5][11] = v[1450] * v[1640] + v[1458] * v[1641] + v[1465] * v[1642] + v[3343];
	Kc[6][0] = -v[3931];
	Kc[6][1] = -v[3932];
	Kc[6][2] = -v[3933];
	Kc[6][3] = -v[3934];
	Kc[6][4] = -v[3935];
	Kc[6][5] = -v[3936];
	Kc[6][6] = v[3931];
	Kc[6][7] = v[3932];
	Kc[6][8] = v[3933];
	Kc[6][9] = -v[3937];
	Kc[6][10] = -v[3938];
	Kc[6][11] = -v[3939];
	Kc[7][0] = -v[3940];
	Kc[7][1] = -v[3941];
	Kc[7][2] = -v[3942];
	Kc[7][3] = -v[3943];
	Kc[7][4] = -v[3944];
	Kc[7][5] = -v[3945];
	Kc[7][6] = v[3940];
	Kc[7][7] = v[3941];
	Kc[7][8] = v[3942];
	Kc[7][9] = -v[3946];
	Kc[7][10] = -v[3947];
	Kc[7][11] = -v[3948];
	Kc[8][0] = -v[3949];
	Kc[8][1] = -v[3950];
	Kc[8][2] = -v[3951];
	Kc[8][3] = -v[3952];
	Kc[8][4] = -v[3953];
	Kc[8][5] = -v[3954];
	Kc[8][6] = v[3949];
	Kc[8][7] = v[3950];
	Kc[8][8] = v[3951];
	Kc[8][9] = -v[3955];
	Kc[8][10] = -v[3956];
	Kc[8][11] = -v[3957];
	Kc[9][0] = v[3982];
	Kc[9][1] = v[3983];
	Kc[9][2] = v[3984];
	Kc[9][3] = v[1617] * v[1690] + v[1618] * v[1691] + v[1619] * v[1692] + v[3326];
	Kc[9][4] = v[1622] * v[1690] + v[1624] * v[1691] + v[1626] * v[1692] + v[3333];
	Kc[9][5] = v[1589] * v[1694] + v[1588] * v[1695] + v[1587] * v[1696] + v[3339];
	Kc[9][6] = -v[3982];
	Kc[9][7] = -v[3983];
	Kc[9][8] = -v[3984];
	Kc[9][9] = v[1587] * v[1632] + v[1588] * v[1633] + v[1589] * v[1634] + (*cn)*(v[3408] * v[3540] + v[3407] * v[3549]
		+ v[3406] * v[3558] + v[1584] * v[3689] + v[1585] * v[3695] + v[1586] * v[3701]) + v[3883] * v[702] + v[3886] * v[708] +
		(*epsn)*((v[1059] * v[1059]) + (v[1062] * v[1062]) + (v[1065] * v[1065]) + v[3557] * v[364] + v[3548] * v[365]
			+ v[3539] * v[366] + v[1123] * v[702] + v[1124] * v[708]) + v[1188] * (v[4348] * v[702] + v[4347] * v[708] + v[1186] * (-
			(v[1081] * v[1177]) + v[1413] * v[2762] + v[1414] * v[2822] - v[3020] + gti[0] * (v[1149] * v[3246] + v[316] * v[3270]
				+ v[3222] * v[331]) + gti[1] * (v[1150] * v[3246] + v[318] * v[3270] + v[3222] * v[334]) + gti[2] * (v[1151] * v[3246]
					+ v[319] * v[3270] + v[3222] * v[336]) - v[3278] * v[363] + v[1351] * v[4342] + v[2954] * v[702] + v[2945] * v[708])
				+ v[1185] * (-(v[1075] * v[1177]) + v[1411] * v[2762] + v[1412] * v[2822] - v[3026] + gti[0] * (v[1149] * v[3247]
					+ v[316] * v[3271] + v[3223] * v[331]) + gti[1] * (v[1150] * v[3247] + v[318] * v[3271] + v[3223] * v[334]) + gti[2] *
					(v[1151] * v[3247] + v[319] * v[3271] + v[3223] * v[336]) - v[3278] * v[361] + v[1351] * v[4337] + v[2970] * v[702]
					+ v[2962] * v[708]) + v[1183] * (-(v[1069] * v[1177]) + v[1409] * v[2762] + v[1410] * v[2822] - v[3032] + gti[0] *
					(v[1149] * v[3248] + v[316] * v[3272] + v[3224] * v[331]) + gti[1] * (v[1150] * v[3248] + v[318] * v[3272]
						+ v[3224] * v[334]) + gti[2] * (v[1151] * v[3248] + v[319] * v[3272] + v[3224] * v[336]) - v[3278] * v[358]
						+ v[1351] * v[4334] + v[2984] * v[702] + v[2977] * v[708]));
	Kc[9][10] = v[1587] * v[1636] + v[1588] * v[1637] + v[1589] * v[1638] + v[3349] + (*cn)*(v[3408] * v[3542]
		+ v[3407] * v[3551] + v[3406] * v[3560] + v[1584] * v[3690] + v[1585] * v[3696] + v[1586] * v[3702]) + v[3883] * v[703]
		+ v[3886] * v[709] + (*epsn)*(v[1135] + v[1123] * v[703] + v[1124] * v[709]);
	Kc[9][11] = v[1587] * v[1640] + v[1588] * v[1641] + v[1589] * v[1642] + v[3350] + (*cn)*(v[3408] * v[3544]
		+ v[3407] * v[3553] + v[3406] * v[3562] + v[1584] * v[3691] + v[1585] * v[3697] + v[1586] * v[3703]) + v[3883] * v[704]
		+ v[3886] * v[710] + (*epsn)*(v[1144] + v[1123] * v[704] + v[1124] * v[710]);
	Kc[10][0] = v[3991];
	Kc[10][1] = v[3992];
	Kc[10][2] = v[3993];
	Kc[10][3] = v[1617] * v[1705] + v[1618] * v[1706] + v[1619] * v[1707] + v[3329];
	Kc[10][4] = v[1622] * v[1705] + v[1624] * v[1706] + v[1626] * v[1707] + v[3336];
	Kc[10][5] = v[1596] * v[1694] + v[1595] * v[1695] + v[1594] * v[1696] + v[3342];
	Kc[10][6] = -v[3991];
	Kc[10][7] = -v[3992];
	Kc[10][8] = -v[3993];
	Kc[10][9] = v[1594] * v[1632] + v[1595] * v[1633] + v[1596] * v[1634] + v[3349] + (*cn)*(v[3408] * v[3566]
		+ v[3407] * v[3574] + v[3406] * v[3582] + v[1591] * v[3689] + v[1592] * v[3695] + v[1593] * v[3701]) + v[3884] * v[702]
		+ v[3887] * v[708] + (*epsn)*(v[1135] + v[1131] * v[702] + v[1132] * v[708]);
	Kc[10][10] = v[1594] * v[1636] + v[1595] * v[1637] + v[1596] * v[1638] + (*cn)*(v[3408] * v[3568] + v[3407] * v[3576]
		+ v[3406] * v[3584] + v[1591] * v[3690] + v[1592] * v[3696] + v[1593] * v[3702]) + v[3884] * v[703] + v[3887] * v[709] +
		(*epsn)*((v[1060] * v[1060]) + (v[1063] * v[1063]) + (v[1066] * v[1066]) + v[3583] * v[364] + v[3575] * v[365]
			+ v[3567] * v[366] + v[1131] * v[703] + v[1132] * v[709]) + v[1188] * (v[4366] * v[703] + v[4365] * v[709] + v[1186] * (-
			(v[1095] * v[1177]) + v[1413] * v[2777] + v[1414] * v[2837] - v[3019] + gti[0] * (v[1149] * v[3240] + v[316] * v[3264]
				+ v[3216] * v[331]) + gti[1] * (v[1150] * v[3240] + v[318] * v[3264] + v[3216] * v[334]) + gti[2] * (v[1151] * v[3240]
					+ v[319] * v[3264] + v[3216] * v[336]) - v[3277] * v[363] + v[1352] * v[4343] + v[2955] * v[703] + v[2946] * v[709])
				+ v[1185] * (-(v[1091] * v[1177]) + v[1411] * v[2777] + v[1412] * v[2837] - v[3025] + gti[0] * (v[1149] * v[3242]
					+ v[316] * v[3266] + v[3218] * v[331]) + gti[1] * (v[1150] * v[3242] + v[318] * v[3266] + v[3218] * v[334]) + gti[2] *
					(v[1151] * v[3242] + v[319] * v[3266] + v[3218] * v[336]) - v[3277] * v[361] + v[1352] * v[4340] + v[2971] * v[703]
					+ v[2963] * v[709]) + v[1183] * (-(v[1087] * v[1177]) + v[1409] * v[2777] + v[1410] * v[2837] - v[3031] + gti[0] *
					(v[1149] * v[3244] + v[316] * v[3268] + v[3220] * v[331]) + gti[1] * (v[1150] * v[3244] + v[318] * v[3268]
						+ v[3220] * v[334]) + gti[2] * (v[1151] * v[3244] + v[319] * v[3268] + v[3220] * v[336]) - v[3277] * v[358]
						+ v[1352] * v[4335] + v[2985] * v[703] + v[2978] * v[709]));
	Kc[10][11] = v[1594] * v[1640] + v[1595] * v[1641] + v[1596] * v[1642] + v[3356] + (*cn)*(v[3408] * v[3570]
		+ v[3407] * v[3578] + v[3406] * v[3586] + v[1591] * v[3691] + v[1592] * v[3697] + v[1593] * v[3703]) + v[3884] * v[704]
		+ v[3887] * v[710] + (*epsn)*(v[1146] + v[1131] * v[704] + v[1132] * v[710]);
	Kc[11][0] = v[4000];
	Kc[11][1] = v[4001];
	Kc[11][2] = v[4002];
	Kc[11][3] = v[1617] * v[1717] + v[1618] * v[1718] + v[1619] * v[1719] + v[3330];
	Kc[11][4] = v[1622] * v[1717] + v[1624] * v[1718] + v[1626] * v[1719] + v[3337];
	Kc[11][5] = v[1603] * v[1694] + v[1602] * v[1695] + v[1601] * v[1696] + v[3343];
	Kc[11][6] = -v[4000];
	Kc[11][7] = -v[4001];
	Kc[11][8] = -v[4002];
	Kc[11][9] = v[1601] * v[1632] + v[1602] * v[1633] + v[1603] * v[1634] + v[3350] + (*cn)*(v[3408] * v[3590]
		+ v[3407] * v[3597] + v[3406] * v[3604] + v[1598] * v[3689] + v[1599] * v[3695] + v[1600] * v[3701]) + v[3885] * v[702]
		+ v[3888] * v[708] + (*epsn)*(v[1144] + v[1140] * v[702] + v[1141] * v[708]);
	Kc[11][10] = v[1601] * v[1636] + v[1602] * v[1637] + v[1603] * v[1638] + v[3356] + (*cn)*(v[3408] * v[3591]
		+ v[3407] * v[3598] + v[3406] * v[3605] + v[1598] * v[3690] + v[1599] * v[3696] + v[1600] * v[3702]) + v[3885] * v[703]
		+ v[3888] * v[709] + (*epsn)*(v[1146] + v[1140] * v[703] + v[1141] * v[709]);
	Kc[11][11] = v[1601] * v[1640] + v[1602] * v[1641] + v[1603] * v[1642] + (*cn)*(v[3408] * v[3593] + v[3407] * v[3600]
		+ v[3406] * v[3607] + v[1598] * v[3691] + v[1599] * v[3697] + v[1600] * v[3703]) + v[3885] * v[704] + v[3888] * v[710] +
		(*epsn)*((v[1061] * v[1061]) + (v[1064] * v[1064]) + (v[1067] * v[1067]) + v[3606] * v[364] + v[3599] * v[365]
			+ v[3592] * v[366] + v[1140] * v[704] + v[1141] * v[710]) + v[1188] * (v[1186] * (-(v[1103] * v[1177]) + v[1413] * v[2794]
				+ v[1414] * v[2854] - v[3017] + gti[0] * (v[1149] * v[3230] + v[316] * v[3254] + v[3206] * v[331]) + gti[1] *
				(v[1150] * v[3230] + v[318] * v[3254] + v[3206] * v[334]) + gti[2] * (v[1151] * v[3230] + v[319] * v[3254]
					+ v[3206] * v[336]) - v[3275] * v[363] + v[1353] * v[4346] + v[2956] * v[704] + v[2947] * v[710]) + v[1185] * (-
					(v[1101] * v[1177]) + v[1411] * v[2794] + v[1412] * v[2854] - v[3023] + gti[0] * (v[1149] * v[3233] + v[316] * v[3257]
						+ v[3209] * v[331]) + gti[1] * (v[1150] * v[3233] + v[318] * v[3257] + v[3209] * v[334]) + gti[2] * (v[1151] * v[3233]
							+ v[319] * v[3257] + v[3209] * v[336]) - v[3275] * v[361] + v[1353] * v[4341] + v[2972] * v[704] + v[2964] * v[710])
				+ v[1183] * (-(v[1099] * v[1177]) + v[1409] * v[2794] + v[1410] * v[2854] - v[3029] + gti[0] * (v[1149] * v[3236]
					+ v[316] * v[3260] + v[3212] * v[331]) + gti[1] * (v[1150] * v[3236] + v[318] * v[3260] + v[3212] * v[334]) + gti[2] *
					(v[1151] * v[3236] + v[319] * v[3260] + v[3212] * v[336]) - v[3275] * v[358] + v[1353] * v[4336] + v[2986] * v[704]
					+ v[2979] * v[710]) + v[704] * (v[1186] * (v[1413] * v[2051] + v[1414] * v[2073] + v[2956] + v[2118] * v[704]
						+ v[2116] * v[710]) + v[1185] * (v[1411] * v[2051] + v[1412] * v[2073] + v[2972] + v[2121] * v[704] + v[2119] * v[710])
						+ v[1183] * (v[1409] * v[2051] + v[1410] * v[2073] + v[2986] + v[2124] * v[704] + v[2122] * v[710])) + v[710] * (v[1186] *
						(v[1413] * v[2052] + v[1414] * v[2076] + v[2947] + v[2116] * v[704] + v[2117] * v[710]) + v[1185] * (v[1411] * v[2052]
							+ v[1412] * v[2076] + v[2964] + v[2119] * v[704] + v[2120] * v[710]) + v[1183] * (v[1409] * v[2052] + v[1410] * v[2076]
								+ v[2979] + v[2122] * v[704] + v[2123] * v[710])));

}