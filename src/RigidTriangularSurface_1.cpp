#include "RigidTriangularSurface_1.h"

#include"Database.h"
//Variáveis globais
extern
Database db;

RigidTriangularSurface_1::RigidTriangularSurface_1()
{
	nDOFs = 6;
	pilot_node = 1;
	pilot_is_used = true;
	n_nodes = 1;
	points = new int[3];
	dA_i = new Matrix(3);
	dB_i = new Matrix(3);
	dC_i = new Matrix(3);
	xP_i = new Matrix(3);

	dA_p = new Matrix(3);
	dB_p = new Matrix(3);
	dC_p = new Matrix(3);
	xP_p = new Matrix(3);

	I3 = new Matrix(3, 3);
	(*I3)(0, 0) = 1.0;
	(*I3)(1, 1) = 1.0;
	(*I3)(2, 2) = 1.0;
	vNR = new Matrix(2);
	GLs = new int*[nDOFs];
	for (int i = 0; i < nDOFs; i++)
		GLs[i] = NULL;

	InitializeDegeneration();
}

RigidTriangularSurface_1::~RigidTriangularSurface_1()
{
	delete[]points;
	delete dA_i;
	delete dB_i;
	delete dC_i;
	delete xP_i;

	delete dA_p;
	delete dB_p;
	delete dC_p;
	delete xP_p;

	delete I3;
	delete vNR;
	delete[]GLs;

	FreeDegeneration();
}

void RigidTriangularSurface_1::SetMinMaxRange()
{
	u1_min = -1.0;
	u1_max = +1.0;
	u1_range = 2.0;

	u2_min = -1.0;
	u2_max = +1.0;
	u2_range = 2.0;
}

//Obtem ponto da superficie
void RigidTriangularSurface_1::SurfacePoint(double& zeta, double& theta, Matrix& point)
{
	//TODO
}

//Normal exterior à superfície na posição escolhida
void RigidTriangularSurface_1::NormalExt(double* zeta, double* theta, Matrix* n)
{

}

//Atualiza bounding box
void RigidTriangularSurface_1::UpdateBox()
{

}

void RigidTriangularSurface_1::WriteVTK_XMLRender(FILE *f)
{
	if (db.post_files->WriteRigidContactSurfaces_flag == true)
	{
		//vetores para escrita no formato binário - usando a função 'enconde'
		std::vector<float> float_vector;
		std::vector<int> int_vector;

		int n_points = 4;	//três da superfície + pilot node
		int n_cells = 2;	//triângulo + pilot
		Matrix pilot(3);
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
		float_vector.push_back((float)(pilot(0, 0)));
		float_vector.push_back((float)(pilot(1, 0)));
		float_vector.push_back((float)(pilot(2, 0)));
		//fprintf(f, "\t\t\t\t\t%f\t%f\t%f\n", pilot(0, 0), pilot(1, 0), pilot(2, 0));
		//Pontos A,B,C
		//fprintf(f, "\t\t\t\t\t%f\t%f\t%f\n", pilot(0, 0) + (*dA_i)(0, 0), pilot(1, 0) + (*dA_i)(1, 0), pilot(2, 0) + (*dA_i)(2, 0));
		//fprintf(f, "\t\t\t\t\t%f\t%f\t%f\n", pilot(0, 0) + (*dB_i)(0, 0), pilot(1, 0) + (*dB_i)(1, 0), pilot(2, 0) + (*dB_i)(2, 0));
		//fprintf(f, "\t\t\t\t\t%f\t%f\t%f\n", pilot(0, 0) + (*dC_i)(0, 0), pilot(1, 0) + (*dC_i)(1, 0), pilot(2, 0) + (*dC_i)(2, 0));
		float_vector.push_back((float)(pilot(0, 0) + (*dA_i)(0, 0)));
		float_vector.push_back((float)(pilot(1, 0) + (*dA_i)(1, 0)));
		float_vector.push_back((float)(pilot(2, 0) + (*dA_i)(2, 0)));
		float_vector.push_back((float)(pilot(0, 0) + (*dB_i)(0, 0)));
		float_vector.push_back((float)(pilot(1, 0) + (*dB_i)(1, 0)));
		float_vector.push_back((float)(pilot(2, 0) + (*dB_i)(2, 0)));
		float_vector.push_back((float)(pilot(0, 0) + (*dC_i)(0, 0)));
		float_vector.push_back((float)(pilot(1, 0) + (*dC_i)(1, 0)));
		float_vector.push_back((float)(pilot(2, 0) + (*dC_i)(2, 0)));
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
		int_vector.clear();
		int_vector.push_back(0);
		int_vector.push_back(1);
		int_vector.push_back(2);
		int_vector.push_back(3);
		//fprintf(f, "\t\t\t\t\t%d\n", 0);
		//fprintf(f, "\t\t\t\t\t%d\t%d\t%d\n", 1,2,3);
		fprintf(f, encodeData(int_vector).c_str());
		fprintf(f, "\n");
		//Closes DataArray
		fprintf(f, "\t\t\t\t</DataArray>\n");
		//Opens DataArray
		fprintf(f, "\t\t\t\t<DataArray type=\"Int32\" Name=\"types\" format=\"binary\">\n");
		int_vector.clear();
		int_vector.push_back(1);
		int_vector.push_back(5);
		//fprintf(f, "\t\t\t\t\t%d\n", 1);
		//fprintf(f, "\t\t\t\t\t%d\n", 5);
		fprintf(f, encodeData(int_vector).c_str());
		fprintf(f, "\n");
		//Closes DataArray
		fprintf(f, "\t\t\t\t</DataArray>\n");
		//Opens DataArray
		fprintf(f, "\t\t\t\t<DataArray type=\"Int32\" Name=\"offsets\" format=\"binary\">\n");
		int_vector.clear();
		int_vector.push_back(1);
		int_vector.push_back(4);
		//fprintf(f, "\t\t\t\t\t%d\n", 1);
		//fprintf(f, "\t\t\t\t\t%d\n", 4);
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

bool RigidTriangularSurface_1::Read(FILE *f)
{
	char s[1000];
	fscanf(f, "%s", s);
	number = atoi(s);

	fscanf(f, "%s", s);
	if (!strcmp(s, "Points"))
	{
		fscanf(f, "%s", s);
		points[0] = atoi(s);
		fscanf(f, "%s", s);
		points[1] = atoi(s);
		fscanf(f, "%s", s);
		points[2] = atoi(s);
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

void RigidTriangularSurface_1::Write(FILE *f)
{
	fprintf(f, "RigidTriangularSurface_1\t%d\tPoints\t%d\t%d\t%d\tPilotNode\t%d\n", number,points[0], points[1], points[2], pilot_node);
}

//Checa inconsistências para evitar erros de execução
bool RigidTriangularSurface_1::Check()
{
	for (int i = 0; i < 3; i++)
	{
		if (points[i] > db.number_points)
			return false;		
	}
	if (pilot_node > db.number_nodes)
		return false;
	return true;
}

//Realiza chute inicial para as variáveis zeta e theta
void RigidTriangularSurface_1::InitialGuess(Matrix* xS, double** convective, int n_solutions)
{
	convective[0][0] = 0.0;
	convective[0][1] = 0.0;
}

void RigidTriangularSurface_1::PreCalc()
{
	Matrix xA = db.points[points[0] - 1]->coordinates;
	Matrix xB = db.points[points[1] - 1]->coordinates;
	Matrix xC = db.points[points[2] - 1]->coordinates;
	//Atribuindo posição inicial do pilot node
	(*xP_i)(0, 0) = db.nodes[pilot_node - 1]->ref_coordinates[0];
	(*xP_i)(1, 0) = db.nodes[pilot_node - 1]->ref_coordinates[1];
	(*xP_i)(2, 0) = db.nodes[pilot_node - 1]->ref_coordinates[2];
	//Cálculo dos vetores dAi, dBi e dCi
	*dA_i = xA - *xP_i;
	*dB_i = xB - *xP_i;
	*dC_i = xC - *xP_i;

	//Apontando para posição que indica valor dos GLs globais
	for (int i = 0; i < 6; i++)
		GLs[i] = &db.nodes[pilot_node - 1]->GLs[i];

	DegenerationPreCalc();
}

//Retorna as coordenadas da superfície para um par (zeta,theta) - configuração anterior convergida
void RigidTriangularSurface_1::Gamma_and_Triad(Matrix* G_p, Matrix* t1_p, Matrix* t2_p, Matrix* n_p, Matrix* G_i, Matrix* t1_i, Matrix* t2_i, Matrix* n_i, Matrix* G_ip, double* zi, double* thi, double* zp, double* thp)
{
	double *dAp = dA_p->getMatrix();	//ponteiro para o vetor dA
	double *dBp = dB_p->getMatrix();	//ponteiro para o vetor dB
	double *dCp = dC_p->getMatrix();	//ponteiro para o vetor dC
	double *xPp = xP_p->getMatrix();	//ponteiro para o vetor xP
	
	double *dAi = dA_i->getMatrix();	//ponteiro para o vetor dA
	double *dBi = dB_i->getMatrix();	//ponteiro para o vetor dB
	double *dCi = dC_i->getMatrix();	//ponteiro para o vetor dC
	double *xPi = xP_i->getMatrix();	//ponteiro para o vetor xP
	
	double *Gp = G_p->getMatrix();		//ponteiro para o vetor Gp 
	double *t1p = t1_p->getMatrix();	//ponteiro para o vetor t1p 
	double *t2p = t2_p->getMatrix();	//ponteiro para o vetor t2p
	double *np = n_p->getMatrix();		//ponteiro para o vetor np

	double *t1i = t1_i->getMatrix();	//ponteiro para o vetor t1i 
	double *t2i = t2_i->getMatrix();	//ponteiro para o vetor t2i
	double *ni = n_i->getMatrix();		//ponteiro para o vetor ni

	double *Gip = G_ip->getMatrix();	//ponteiro para o vetor Gip
	double *Gi = G_i->getMatrix();		//ponteiro para o vetor Gi
	double v[1000];
	//Código - AceGen
	v[169] = (1e0 + (*thp)) / 2e0;
	v[168] = (1e0 + (*zp)) / 2e0;
	v[167] = (-(*thp) - (*zp)) / 2e0;
	v[166] = (1e0 + (*thi)) / 2e0;
	v[165] = (1e0 + (*zi)) / 2e0;
	v[164] = (-(*thi) - (*zi)) / 2e0;
	v[163] = dCp[2] + xPp[2];
	v[162] = dCp[1] + xPp[1];
	v[161] = dCp[0] + xPp[0];
	v[160] = dBp[2] + xPp[2];
	v[159] = dBp[1] + xPp[1];
	v[158] = dBp[0] + xPp[0];
	v[157] = dAp[2] + xPp[2];
	v[156] = dAp[1] + xPp[1];
	v[155] = dAp[0] + xPp[0];
	v[154] = dCi[2] + xPi[2];
	v[153] = dCi[1] + xPi[1];
	v[152] = dCi[0] + xPi[0];
	v[151] = dBi[2] + xPi[2];
	v[150] = dBi[1] + xPi[1];
	v[149] = dBi[0] + xPi[0];
	v[148] = dAi[2] + xPi[2];
	v[147] = dAi[1] + xPi[1];
	v[146] = dAi[0] + xPi[0];
	v[96] = -v[146] / 2e0;
	v[98] = -v[147] / 2e0;
	v[100] = -v[148] / 2e0;
	v[120] = -v[155] / 2e0;
	v[122] = -v[156] / 2e0;
	v[124] = -v[157] / 2e0;
	v[89] = v[149] / 2e0 + v[96];
	v[90] = v[150] / 2e0 + v[98];
	v[91] = v[100] + v[151] / 2e0;
	v[93] = 1e0 / sqrt(Power(v[89], 2) + Power(v[90], 2) + Power(v[91], 2));
	v[92] = v[89] * v[93];
	v[94] = v[90] * v[93];
	v[95] = v[91] * v[93];
	v[97] = v[152] / 2e0 + v[96];
	v[99] = v[153] / 2e0 + v[98];
	v[101] = v[100] + v[154] / 2e0;
	v[103] = 1e0 / sqrt(Power(v[101], 2) + Power(v[97], 2) + Power(v[99], 2));
	v[102] = v[103] * v[97];
	v[104] = v[103] * v[99];
	v[111] = v[104] * v[92] - v[102] * v[94];
	v[105] = v[101] * v[103];
	v[109] = -(v[105] * v[92]) + v[102] * v[95];
	v[106] = v[105] * v[94] - v[104] * v[95];
	v[108] = 1e0 / sqrt(Power(v[106], 2) + Power(v[109], 2) + Power(v[111], 2));
	v[113] = v[120] + v[158] / 2e0;
	v[114] = v[122] + v[159] / 2e0;
	v[115] = v[124] + v[160] / 2e0;
	v[117] = 1e0 / sqrt(Power(v[113], 2) + Power(v[114], 2) + Power(v[115], 2));
	v[116] = v[113] * v[117];
	v[118] = v[114] * v[117];
	v[119] = v[115] * v[117];
	v[121] = v[120] + v[161] / 2e0;
	v[123] = v[122] + v[162] / 2e0;
	v[125] = v[124] + v[163] / 2e0;
	v[127] = 1e0 / sqrt(Power(v[121], 2) + Power(v[123], 2) + Power(v[125], 2));
	v[126] = v[121] * v[127];
	v[128] = v[123] * v[127];
	v[135] = -(v[118] * v[126]) + v[116] * v[128];
	v[129] = v[125] * v[127];
	v[133] = v[119] * v[126] - v[116] * v[129];
	v[130] = -(v[119] * v[128]) + v[118] * v[129];
	v[132] = 1e0 / sqrt(Power(v[130], 2) + Power(v[133], 2) + Power(v[135], 2));
	Gip[0] = v[155] * v[164] + v[158] * v[165] + v[161] * v[166];
	Gip[1] = v[156] * v[164] + v[159] * v[165] + v[162] * v[166];
	Gip[2] = v[157] * v[164] + v[160] * v[165] + v[163] * v[166];
	t1i[0] = v[92];
	t1i[1] = v[94];
	t1i[2] = v[95];
	t2i[0] = v[102];
	t2i[1] = v[104];
	t2i[2] = v[105];
	ni[0] = v[106] * v[108];
	ni[1] = v[108] * v[109];
	ni[2] = v[108] * v[111];
	Gp[0] = v[155] * v[167] + v[158] * v[168] + v[161] * v[169];
	Gp[1] = v[156] * v[167] + v[159] * v[168] + v[162] * v[169];
	Gp[2] = v[157] * v[167] + v[160] * v[168] + v[163] * v[169];
	t1p[0] = v[116];
	t1p[1] = v[118];
	t1p[2] = v[119];
	t2p[0] = v[126];
	t2p[1] = v[128];
	t2p[2] = v[129];
	np[0] = v[130] * v[132];
	np[1] = v[132] * v[133];
	np[2] = v[132] * v[135];
	Gi[0] = v[146] * v[164] + v[149] * v[165] + v[152] * v[166];
	Gi[1] = v[147] * v[164] + v[150] * v[165] + v[153] * v[166];
	Gi[2] = v[148] * v[164] + v[151] * v[165] + v[154] * v[166];
}

//Dado o ponto xS, calcula as coordenadas (zeta,theta) referentes à mínima distância
void RigidTriangularSurface_1::FindMinimimumParameters(Matrix* xS, NSContactData* cd)
{
	
	//Analitical solution
	Matrix A = *xP_p + *dA_p;
	Matrix B = *xP_p + *dB_p;
	Matrix C = *xP_p + *dC_p;
	//Normal
	Matrix n = (1.0 / norm(cross(B - A, C - A)))*cross((B - A), (C - A));
	//Projection
	Matrix xS_A = (*xS) - dot(*xS, n)*n;
	double d = 0.25*(B - A)(0, 0)*(C - A)(2, 0) - 0.25*(B - A)(2, 0)*(C - A)(0, 0);
	double dzeta = 0.5*(xS_A - 0.5*(B + C))(0, 0)*(C - A)(2, 0) - 0.5*(xS_A - 0.5*(B + C))(2, 0)*(C - A)(0, 0);
	double dtheta = 0.5*(B - A)(0, 0)*(xS_A - 0.5*(B + C))(2, 0) - 0.5*(B - A)(2, 0)*(xS_A - 0.5*(B + C))(0, 0);
	double zeta = dzeta / d;
	double theta = dtheta/d;

	////Salva nas variáveis
	//cd->convective[0][0] = zeta;
	//cd->convective[0][1] = theta;
	//
	////Se está no range local de interesse - domínio físico da superfície triangular
	//if (abs(zeta) <= 1.0 && abs(theta) <= 1.0 && theta <= -zeta)
	//	cd->return_value[0] = 0;
	//else
	//{
	//	//Se está em região próxima, mas não no range local de interesse
	//	if (abs(zeta) < 1.2 && abs(theta) < 1.2 && theta < -zeta + 0.2)
	//		cd->return_value[0] = 3;
	//	//Se não estiver no range de interesse
	//	else
	//		cd->return_value[0] = 2;
	//}

	//cd->repeated[0] = false;
	//for (int i = 1; i < cd->n_solutions; i++)
	//	cd->repeated[i] = true;
	
	//Parametros NR
	double tol_ortho = 1e-10;
	int max_it = 20;

	double** J;
	//Alocação J
	J = new double*[2];
	for (int i = 0; i < 2; i++)
		J[i] = new double[2];
	Matrix Jacobian(2, 2);
	Matrix residual(2);
	Matrix delta(2);
	double *dA = dA_p->getMatrix();		//ponteiro para o vetor dA
	double *dB = dB_p->getMatrix();		//ponteiro para o vetor dB
	double *dC = dC_p->getMatrix();		//ponteiro para o vetor dC
	double *xP = xP_p->getMatrix();		//ponteiro para o vetor xP
	double *pxS = xS->getMatrix();			//ponteiro para o vetor xS
	double *r = residual.getMatrix();		//ponteiro para o vetor residual
	double *vi = vNR->getMatrix();			//pointeiro para vNR
	
	//Inicialização de chute inicial
	(*vNR)(0, 0) = cd->convective[0][0];
	(*vNR)(1, 0) = cd->convective[0][1];
	double error = tol_ortho + 1.0;	//Forçando entrar no loop
	int it = 1;
	int flag_error;
	double v[200];
	while (error > tol_ortho && it <= max_it)
	{
		if (error > tol_ortho)
		{
			//Cálculo do resíduo e jacobiano
			v[70] = (1e0 + vi[1]) / 2e0;
			v[69] = (1e0 + vi[0]) / 2e0;
			v[68] = (-vi[0] - vi[1]) / 2e0;
			v[67] = dC[2] + xP[2];
			v[66] = dC[1] + xP[1];
			v[65] = dC[0] + xP[0];
			v[64] = dB[2] + xP[2];
			v[63] = dB[1] + xP[1];
			v[62] = dB[0] + xP[0];
			v[61] = dA[2] + xP[2];
			v[60] = dA[1] + xP[1];
			v[59] = dA[0] + xP[0];
			v[47] = -v[59] / 2e0;
			v[45] = -v[60] / 2e0;
			v[43] = -v[61] / 2e0;
			v[42] = v[47] + v[62] / 2e0;
			v[41] = v[45] + v[63] / 2e0;
			v[40] = v[43] + v[64] / 2e0;
			v[48] = v[47] + v[65] / 2e0;
			v[46] = v[45] + v[66] / 2e0;
			v[44] = v[43] + v[67] / 2e0;
			v[71] = v[40] * v[44] + v[41] * v[46] + v[42] * v[48];
			v[51] = -(v[59] * v[68]) - v[62] * v[69] - v[65] * v[70] + pxS[0];
			v[74] = -2e0*v[51];
			v[50] = -(v[60] * v[68]) - v[63] * v[69] - v[66] * v[70] + pxS[1];
			v[73] = -2e0*v[50];
			v[49] = -(v[61] * v[68]) - v[64] * v[69] - v[67] * v[70] + pxS[2];
			v[72] = -2e0*v[49];
			r[0] = 0.5e0*(v[40] * v[72] + v[41] * v[73] + v[42] * v[74]);
			r[1] = 0.5e0*(v[44] * v[72] + v[46] * v[73] + v[48] * v[74]);
			J[0][0] = 1e0*((v[40] * v[40]) + (v[41] * v[41]) + (v[42] * v[42]));
			J[0][1] = v[71];
			J[1][0] = v[71];
			J[1][1] = 1e0*((v[44] * v[44]) + (v[46] * v[46]) + (v[48] * v[48]));
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
			/*(*vNR).print();
			delta.print();
			residual.print();*/
		}
	}//end while NR

	//Clean
	for (int i = 0; i < 2; i++)
		delete[]J[i];
	delete[]J;
	
	//Convergiu - ainda há ações a verificar...
	if (error <= tol_ortho && flag_error == 0)
	{
		//Salva nas variáveis
		cd->convective[0][0] = (*vNR)(0, 0);
		cd->convective[0][1] = (*vNR)(1, 0);
		//Se está no range local de interesse - domínio físico da superfície triangular
		if (abs((*vNR)(0, 0)) <= 1.0 && abs((*vNR)(1, 0)) <= 1.0 && (*vNR)(1, 0) <= -(*vNR)(0, 0))
			cd->return_value[0] = 0;
		else
		{
			//Se está em região próxima, mas não no range local de interesse
			if (abs((*vNR)(0, 0)) < 1.2 && abs((*vNR)(1, 0)) < 1.2 && (*vNR)(1, 0) < -(*vNR)(0, 0) + 0.2)
				cd->return_value[0] = 3;
			//Se não estiver no range de interesse
			else
				cd->return_value[0] = 2;
		}
		
	}
	else
		cd->return_value[0] = 1;
	//Retornos da função
	//0 - Convergiu e está no range de interesse para contato
	//1 - Não houve convergência (pode ou não estar no range para contato) - retorno problemático!!
	//2 - Houve convergência, mas não está no range para contato
	//3 - Houve convergência, está fora do range para contato, mas próximo

	cd->repeated[0] = false;
	for (int i = 1; i < cd->n_solutions; i++)
		cd->repeated[i] = true;
}
//Atualiza as variáveis internas da superfície, para pegarem info do pilot node para uso posterior com posição atualizada
void RigidTriangularSurface_1::FillNodes()
{
	//Atualização do Pilot Node
	(*xP_p)(0, 0) = (*xP_i)(0, 0) + db.nodes[pilot_node - 1]->displacements[0];
	(*xP_p)(1, 0) = (*xP_i)(1, 0) + db.nodes[pilot_node - 1]->displacements[1];
	(*xP_p)(2, 0) = (*xP_i)(2, 0) + db.nodes[pilot_node - 1]->displacements[2];

	//Atualização dos vetores dAtemp, dBtemp e dCtemp
	Matrix alpha_1(3);
	alpha_1(0, 0) = db.nodes[pilot_node - 1]->displacements[3];
	alpha_1(1, 0) = db.nodes[pilot_node - 1]->displacements[4];
	alpha_1(2, 0) = db.nodes[pilot_node - 1]->displacements[5];
	double alpha_escalar = norm(alpha_1);
	Matrix A = skew(alpha_1);
	double g = 4.0 / (4.0 + alpha_escalar*alpha_escalar);
	Matrix Q = *I3 + g*(A + 0.5*((A)*(A)));
	*dA_p = Q*(*dA_i);
	*dB_p = Q*(*dB_i);
	*dC_p = Q*(*dC_i);
}

//Retorna coordenadas globais do ponto central da superfície a ser utilizado para cálculos grosseiros de sua localização (pinball)
void RigidTriangularSurface_1::CenterPoint(Matrix* center)
{
	*center = *xP_i + 0.3333333333333333333333*(*dA_i + *dB_i + *dC_i);
}

//Salva vetores de configuração convergida
void RigidTriangularSurface_1::SaveConfiguration()
{
	*xP_i = *xP_p;
	*dA_i = *dA_p;
	*dB_i = *dB_p;
	*dC_i = *dC_p;
}


//Calcula contribuições de contato entre esfera e superfície
void RigidTriangularSurface_1::ContactSphereSurfaceSticking(double* Rc, double** Kc, double zetap, double thetap, double zetai, double thetai, double* gti, int node, double* epsn, double* epst, double* cn, double* ct, double* mu, double* radius)
{
	double *dAi = dA_i->getMatrix();	//ponteiro para o vetor dA
	double *dBi = dB_i->getMatrix();	//ponteiro para o vetor dB
	double *dCi = dC_i->getMatrix();	//ponteiro para o vetor dC
	double *xPi = xP_i->getMatrix();	//ponteiro para o vetor xP
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
	//Código - AceGen
	v[1059] = 0.5e0*d[3];
	v[1057] = 2e0*d[3];
	v[266] = Power(d[3], 2);
	v[1060] = 2e0*d[4];
	v[1058] = 0.5e0*d[4];
	v[264] = 0.5e0*d[3] * d[4];
	v[259] = Power(d[4], 2);
	v[1068] = -v[259] - v[266];
	v[2794] = 0.5e0*v[1068];
	v[2792] = -0.5e0*v[1068];
	v[1094] = d[5] + v[264];
	v[1084] = -d[5] + v[264];
	v[1062] = 2e0*d[5];
	v[1061] = 0.5e0*d[5];
	v[271] = 0.5e0*d[4] * d[5];
	v[1106] = d[3] + v[271];
	v[1098] = -d[3] + v[271];
	v[269] = 0.5e0*d[3] * d[5];
	v[1102] = -d[4] + v[269];
	v[1090] = d[4] + v[269];
	v[260] = Power(d[5], 2);
	v[1077] = -v[259] - v[260];
	v[2796] = 0.5e0*v[1077];
	v[2795] = 0.5e0*v[1077];
	v[1073] = -v[260] - v[266];
	v[2789] = -0.5e0*v[1073];
	v[1063] = 4e0 + v[259] + v[260] + v[266];
	v[1279] = 1e0 / Power(v[1063], 3);
	v[1280] = -2e0*v[1060] * v[1279];
	v[1283] = -4e0*v[1062] * v[1280];
	v[1278] = -2e0*v[1057] * v[1279];
	v[1288] = -4e0*v[1060] * v[1278];
	v[1282] = -4e0*v[1062] * v[1278];
	v[1065] = 1e0 / Power(v[1063], 2);
	v[1289] = -8e0*v[1065];
	v[1295] = -4e0*v[1057] * v[1278] + v[1289];
	v[2481] = -v[1282] - v[1058] * v[1295];
	v[2475] = v[1288] - v[1061] * v[1295];
	v[1290] = -4e0*v[1060] * v[1280] + v[1289];
	v[1284] = 8e0*(v[1062] * v[1062])*v[1279] + v[1289];
	v[1067] = -4e0*v[1062] * v[1065];
	v[1685] = -(v[1061] * v[1067]);
	v[1351] = -0.5e0*v[1062] * v[1067];
	v[1334] = 2e0*v[1067];
	v[1299] = -0.5e0*v[1057] * v[1067];
	v[1293] = -0.5e0*v[1060] * v[1067];
	v[1066] = -4e0*v[1060] * v[1065];
	v[1354] = -0.5e0*v[1060] * v[1066];
	v[1327] = v[1058] * v[1066];
	v[1318] = -2e0*v[1066];
	v[1313] = v[1061] * v[1066];
	v[1297] = -0.5e0*v[1057] * v[1066];
	v[1064] = -4e0*v[1057] * v[1065];
	v[2791] = v[1064] - v[1313];
	v[2464] = 0.5e0*v[1064];
	v[1359] = -0.5e0*v[1057] * v[1064];
	v[1330] = v[1059] * v[1064];
	v[1310] = v[1061] * v[1064];
	v[2827] = -v[1066] - v[1310];
	v[2790] = v[1066] - v[1310];
	v[1308] = 2e0*v[1064];
	v[1306] = v[1058] * v[1064];
	v[2793] = v[1067] - v[1306];
	v[472] = 0.5e0*d[9];
	v[470] = 2e0*d[9];
	v[247] = Power(d[9], 2);
	v[473] = 2e0*d[10];
	v[471] = 0.5e0*d[10];
	v[245] = 0.5e0*d[10] * d[9];
	v[240] = Power(d[10], 2);
	v[415] = -v[240] - v[247];
	v[2786] = 0.5e0*v[415];
	v[475] = 2e0*d[11];
	v[474] = 0.5e0*d[11];
	v[429] = -d[11] + v[245];
	v[425] = d[11] + v[245];
	v[252] = 0.5e0*d[10] * d[11];
	v[420] = -d[9] + v[252];
	v[418] = d[9] + v[252];
	v[250] = 0.5e0*d[11] * d[9];
	v[426] = d[10] + v[250];
	v[419] = -d[10] + v[250];
	v[241] = Power(d[11], 2);
	v[436] = 4e0 + v[240] + v[241] + v[247];
	v[688] = 1e0 / Power(v[436], 3);
	v[690] = -2e0*v[475] * v[688];
	v[1975] = -4e0*v[690];
	v[689] = -2e0*v[473] * v[688];
	v[1974] = -4e0*v[689];
	v[692] = -4e0*v[475] * v[689];
	v[687] = -2e0*v[470] * v[688];
	v[1973] = -4e0*v[687];
	v[697] = -4e0*v[473] * v[687];
	v[691] = -4e0*v[475] * v[687];
	v[476] = 1e0 / Power(v[436], 2);
	v[698] = -8e0*v[476];
	v[700] = -4e0*v[470] * v[687] + v[698];
	v[699] = -4e0*v[473] * v[689] + v[698];
	v[693] = -4e0*v[475] * v[690] + v[698];
	v[696] = 0.5e0*v[415] * v[693];
	v[479] = -4e0*v[475] * v[476];
	v[1748] = -0.5e0*v[479];
	v[1661] = -(v[474] * v[479]);
	v[786] = 2e0*v[479];
	v[787] = v[429] * v[693] - v[786];
	v[737] = v[425] * v[693] + v[786];
	v[731] = -0.5e0*v[475] * v[479];
	v[708] = -0.5e0*v[473] * v[479];
	v[821] = v[420] * v[693] - v[708];
	v[814] = v[418] * v[693] - v[708];
	v[778] = v[420] * v[699] - v[708];
	v[744] = v[418] * v[699] - v[708];
	v[703] = -0.5e0*v[470] * v[479];
	v[810] = v[426] * v[693] - v[703];
	v[789] = v[419] * v[693] - v[703];
	v[781] = v[426] * v[700] - v[703];
	v[764] = v[419] * v[700] - v[703];
	v[695] = 0.5e0*v[415] * v[692] + v[708];
	v[694] = 0.5e0*v[415] * v[691] + v[703];
	v[512] = 0.5e0*v[415] * v[479];
	v[478] = -4e0*v[473] * v[476];
	v[1735] = 0.5e0*v[478];
	v[779] = v[471] * v[478];
	v[739] = -2e0*v[478];
	v[740] = v[426] * v[699] - v[739];
	v[734] = v[474] * v[478];
	v[721] = v[419] * v[699] + v[739];
	v[712] = -0.5e0*v[473] * v[478];
	v[711] = -1e0*v[470] * v[478] + v[2786] * v[697];
	v[477] = -4e0*v[470] * v[476];
	v[2777] = -v[477] + v[734];
	v[2774] = v[477] + v[734];
	v[1741] = 0.5e0*v[477];
	v[784] = v[2777] + v[429] * v[691];
	v[782] = v[472] * v[477];
	v[738] = v[2774] + v[426] * v[697];
	v[735] = v[2774] + v[425] * v[691];
	v[726] = v[474] * v[477];
	v[2779] = v[478] + v[726];
	v[2775] = -v[478] + v[726];
	v[785] = v[2775] + v[429] * v[692];
	v[736] = v[2779] + v[425] * v[692];
	v[727] = v[2775] + v[420] * v[697];
	v[724] = 2e0*v[477];
	v[725] = v[420] * v[700] - v[724];
	v[722] = v[471] * v[477];
	v[2778] = v[479] + v[722];
	v[2776] = -v[479] + v[722];
	v[741] = v[2778] + v[426] * v[692];
	v[728] = v[2776] + v[420] * v[691];
	v[723] = v[2776] + v[419] * v[692];
	v[720] = v[2777] + v[419] * v[697];
	v[719] = v[2778] + v[418] * v[691];
	v[718] = v[2779] + v[418] * v[697];
	v[717] = v[418] * v[700] + v[724];
	v[715] = -0.5e0*v[470] * v[477];
	v[705] = -0.5e0*v[473] * v[477];
	v[841] = v[429] * v[700] - v[705];
	v[834] = v[425] * v[700] - v[705];
	v[812] = v[429] * v[699] - v[705];
	v[791] = v[425] * v[699] - v[705];
	v[432] = -v[240] - v[241];
	v[2788] = 0.5e0*v[432];
	v[860] = v[2788] * v[692] + 2e0*v[708];
	v[704] = 0.5e0*v[432] * v[691] + v[703];
	v[702] = 0.5e0*v[432] * v[697] + v[705];
	v[701] = 0.5e0*v[432] * v[700];
	v[480] = 0.5e0*v[432] * v[477];
	v[423] = -v[241] - v[247];
	v[2787] = 0.5e0*v[423];
	v[730] = v[2787] * v[691] + 2e0*v[703];
	v[709] = 0.5e0*v[423] * v[692] + v[708];
	v[707] = 0.5e0*v[423] * v[699];
	v[706] = 0.5e0*v[423] * v[697] + v[705];
	v[496] = 0.5e0*v[423] * v[478];
	v[1692] = (*a4)*d[3] + (*a6)*domegaiS[0] + (*a5)*omegaiS[0];
	v[2531] = v[1692] / 2e0;
	v[1696] = (*a4)*d[4] + (*a6)*domegaiS[1] + (*a5)*omegaiS[1];
	v[2527] = -v[1696] / 2e0;
	v[1698] = (*a4)*d[5] + (*a6)*domegaiS[2] + (*a5)*omegaiS[2];
	v[2510] = v[1698] / 2e0;
	v[1668] = (*a4)*d[9] + (*a6)*domegaiP[0] + (*a5)*omegaiP[0];
	v[2213] = v[1668] / 2e0;
	v[1672] = (*a4)*d[10] + (*a6)*domegaiP[1] + (*a5)*omegaiP[1];
	v[2209] = -v[1672] / 2e0;
	v[1674] = (*a4)*d[11] + (*a6)*domegaiP[2] + (*a5)*omegaiP[2];
	v[2192] = v[1674] / 2e0;
	v[866] = dAi[2] * v[741] + dAi[1] * v[785] + dAi[0] * v[860];
	v[868] = -v[866] / 2e0;
	v[848] = dAi[0] * v[701] + dAi[2] * v[781] + dAi[1] * v[841];
	v[851] = -v[848] / 2e0;
	v[826] = dAi[2] * v[728] + dAi[1] * v[730] + dAi[0] * v[735];
	v[828] = -v[826] / 2e0;
	v[817] = dAi[2] * v[696] + dAi[0] * v[789] + dAi[1] * v[814];
	v[818] = -v[817] / 2e0;
	v[799] = dAi[1] * v[707] + dAi[2] * v[778] + dAi[0] * v[791];
	v[802] = -v[799] / 2e0;
	v[752] = dAi[2] * v[711] + dAi[1] * v[718] + dAi[0] * v[720];
	v[755] = -v[752] / 2e0;
	v[864] = dBi[2] * v[741] + dBi[1] * v[785] + dBi[0] * v[860];
	v[845] = dBi[0] * v[701] + dBi[2] * v[781] + dBi[1] * v[841];
	v[824] = dBi[2] * v[728] + dBi[1] * v[730] + dBi[0] * v[735];
	v[816] = dBi[2] * v[696] + dBi[0] * v[789] + dBi[1] * v[814];
	v[796] = dBi[1] * v[707] + dBi[2] * v[778] + dBi[0] * v[791];
	v[749] = dBi[2] * v[711] + dBi[1] * v[718] + dBi[0] * v[720];
	v[862] = dCi[2] * v[741] + dCi[1] * v[785] + dCi[0] * v[860];
	v[842] = dCi[0] * v[701] + dCi[2] * v[781] + dCi[1] * v[841];
	v[822] = dCi[2] * v[728] + dCi[1] * v[730] + dCi[0] * v[735];
	v[815] = dCi[2] * v[696] + dCi[0] * v[789] + dCi[1] * v[814];
	v[793] = dCi[1] * v[707] + dCi[2] * v[778] + dCi[0] * v[791];
	v[746] = dCi[2] * v[711] + dCi[1] * v[718] + dCi[0] * v[720];
	v[224] = d[6] + xPi[0];
	v[225] = d[7] + xPi[1];
	v[226] = d[8] + xPi[2];
	v[239] = 4e0 / v[436];
	v[1665] = (v[239] * v[239]);
	v[742] = -0.5e0*v[239];
	v[2780] = -v[742] + v[779] + v[782];
	v[811] = v[2780] + v[429] * v[697];
	v[849] = dAi[0] * v[702] + dAi[2] * v[738] + dAi[1] * v[811];
	v[852] = -v[849] / 2e0;
	v[846] = dBi[0] * v[702] + dBi[2] * v[738] + dBi[1] * v[811];
	v[843] = dCi[0] * v[702] + dCi[2] * v[738] + dCi[1] * v[811];
	v[790] = v[2780] + v[425] * v[697];
	v[798] = dAi[1] * v[706] + dAi[2] * v[727] + dAi[0] * v[790];
	v[801] = -v[798] / 2e0;
	v[795] = dBi[1] * v[706] + dBi[2] * v[727] + dBi[0] * v[790];
	v[792] = dCi[1] * v[706] + dCi[2] * v[727] + dCi[0] * v[790];
	v[743] = v[1661] + v[742];
	v[2782] = -v[743] + v[779];
	v[2781] = -v[743] + v[782];
	v[783] = v[2781] + v[426] * v[691];
	v[850] = dAi[0] * v[704] + dAi[2] * v[783] + dAi[1] * v[784];
	v[853] = -v[850] / 2e0;
	v[847] = dBi[0] * v[704] + dBi[2] * v[783] + dBi[1] * v[784];
	v[844] = dCi[0] * v[704] + dCi[2] * v[783] + dCi[1] * v[784];
	v[780] = v[2782] + v[420] * v[692];
	v[800] = dAi[1] * v[709] + dAi[0] * v[736] + dAi[2] * v[780];
	v[803] = -v[800] / 2e0;
	v[797] = dBi[1] * v[709] + dBi[0] * v[736] + dBi[2] * v[780];
	v[794] = dCi[1] * v[709] + dCi[0] * v[736] + dCi[2] * v[780];
	v[765] = v[2781] + v[419] * v[691];
	v[771] = dAi[2] * v[694] + dAi[1] * v[719] + dAi[0] * v[765];
	v[773] = -v[771] / 2e0;
	v[769] = dBi[2] * v[694] + dBi[1] * v[719] + dBi[0] * v[765];
	v[767] = dCi[2] * v[694] + dCi[1] * v[719] + dCi[0] * v[765];
	v[745] = v[2782] + v[418] * v[692];
	v[754] = dAi[2] * v[695] + dAi[0] * v[723] + dAi[1] * v[745];
	v[757] = -v[754] / 2e0;
	v[751] = dBi[2] * v[695] + dBi[0] * v[723] + dBi[1] * v[745];
	v[748] = dCi[2] * v[695] + dCi[0] * v[723] + dCi[1] * v[745];
	v[2783] = -v[239] + 2e0*v[731];
	v[861] = v[2783] + 0.5e0*v[432] * v[693];
	v[867] = dAi[1] * v[787] + dAi[2] * v[810] + dAi[0] * v[861];
	v[869] = -v[867] / 2e0;
	v[865] = dBi[1] * v[787] + dBi[2] * v[810] + dBi[0] * v[861];
	v[863] = dCi[1] * v[787] + dCi[2] * v[810] + dCi[0] * v[861];
	v[732] = v[2783] + 0.5e0*v[423] * v[693];
	v[827] = dAi[1] * v[732] + dAi[0] * v[737] + dAi[2] * v[821];
	v[829] = -v[827] / 2e0;
	v[825] = dBi[1] * v[732] + dBi[0] * v[737] + dBi[2] * v[821];
	v[823] = dCi[1] * v[732] + dCi[0] * v[737] + dCi[2] * v[821];
	v[2784] = -v[239] + 2e0*v[715];
	v[733] = v[2784] + 0.5e0*v[423] * v[700];
	v[837] = dAi[2] * v[725] + dAi[1] * v[733] + dAi[0] * v[834];
	v[838] = -v[837] / 2e0;
	v[836] = dBi[2] * v[725] + dBi[1] * v[733] + dBi[0] * v[834];
	v[835] = dCi[2] * v[725] + dCi[1] * v[733] + dCi[0] * v[834];
	v[716] = v[2784] + 0.5e0*v[415] * v[700];
	v[770] = dAi[2] * v[716] + dAi[1] * v[717] + dAi[0] * v[764];
	v[772] = -v[770] / 2e0;
	v[768] = dBi[2] * v[716] + dBi[1] * v[717] + dBi[0] * v[764];
	v[766] = dCi[2] * v[716] + dCi[1] * v[717] + dCi[0] * v[764];
	v[2785] = -v[239] + 2e0*v[712];
	v[874] = v[2785] + 0.5e0*v[432] * v[699];
	v[877] = dAi[2] * v[740] + dAi[1] * v[812] + dAi[0] * v[874];
	v[878] = -v[877] / 2e0;
	v[876] = dBi[2] * v[740] + dBi[1] * v[812] + dBi[0] * v[874];
	v[875] = dCi[2] * v[740] + dCi[1] * v[812] + dCi[0] * v[874];
	v[713] = v[2785] + 0.5e0*v[415] * v[699];
	v[753] = dAi[2] * v[713] + dAi[0] * v[721] + dAi[1] * v[744];
	v[756] = -v[753] / 2e0;
	v[750] = dBi[2] * v[713] + dBi[0] * v[721] + dBi[1] * v[744];
	v[747] = dCi[2] * v[713] + dCi[0] * v[721] + dCi[1] * v[744];
	v[510] = -0.5e0*v[239] * v[473];
	v[511] = v[2786] * v[478] + v[510];
	v[508] = -0.5e0*v[239] * v[470];
	v[509] = 0.5e0*v[415] * v[477] + v[508];
	v[505] = v[239] + v[418] * v[477];
	v[503] = -v[239] + v[419] * v[478];
	v[499] = -v[239] + v[420] * v[477];
	v[497] = -0.5e0*v[239] * v[475];
	v[498] = v[2787] * v[479] + v[497];
	v[495] = 0.5e0*v[423] * v[477] + v[508];
	v[494] = v[239] + v[425] * v[479];
	v[490] = v[239] + v[426] * v[478];
	v[488] = -(v[239] * v[474]);
	v[506] = v[418] * v[478] - v[488];
	v[538] = dCi[0] * v[503] + dCi[1] * v[506] + dCi[2] * v[511];
	v[529] = dBi[0] * v[503] + dBi[1] * v[506] + dBi[2] * v[511];
	v[520] = dAi[0] * v[503] + dAi[1] * v[506] + dAi[2] * v[511];
	v[547] = -v[520] / 2e0;
	v[565] = v[538] / 2e0 + v[547];
	v[556] = v[529] / 2e0 + v[547];
	v[502] = v[419] * v[477] - v[488];
	v[537] = dCi[0] * v[502] + dCi[1] * v[505] + dCi[2] * v[509];
	v[528] = dBi[0] * v[502] + dBi[1] * v[505] + dBi[2] * v[509];
	v[519] = dAi[0] * v[502] + dAi[1] * v[505] + dAi[2] * v[509];
	v[546] = -v[519] / 2e0;
	v[564] = v[537] / 2e0 + v[546];
	v[555] = v[528] / 2e0 + v[546];
	v[500] = v[420] * v[478] - v[488];
	v[489] = v[426] * v[477] - v[488];
	v[487] = -v[239] + v[429] * v[479];
	v[485] = -(v[239] * v[472]);
	v[504] = v[419] * v[479] - v[485];
	v[493] = v[425] * v[478] - v[485];
	v[535] = dCi[0] * v[493] + dCi[1] * v[496] + dCi[2] * v[500];
	v[526] = dBi[0] * v[493] + dBi[1] * v[496] + dBi[2] * v[500];
	v[517] = dAi[0] * v[493] + dAi[1] * v[496] + dAi[2] * v[500];
	v[544] = -v[517] / 2e0;
	v[562] = v[535] / 2e0 + v[544];
	v[553] = v[526] / 2e0 + v[544];
	v[491] = v[426] * v[479] - v[485];
	v[486] = v[429] * v[478] - v[485];
	v[483] = v[239] * v[471];
	v[507] = v[418] * v[479] + v[483];
	v[539] = dCi[0] * v[504] + dCi[1] * v[507] + dCi[2] * v[512];
	v[530] = dBi[0] * v[504] + dBi[1] * v[507] + dBi[2] * v[512];
	v[521] = dAi[0] * v[504] + dAi[1] * v[507] + dAi[2] * v[512];
	v[548] = -v[521] / 2e0;
	v[566] = v[539] / 2e0 + v[548];
	v[557] = v[530] / 2e0 + v[548];
	v[501] = v[420] * v[479] + v[483];
	v[536] = dCi[0] * v[494] + dCi[1] * v[498] + dCi[2] * v[501];
	v[527] = dBi[0] * v[494] + dBi[1] * v[498] + dBi[2] * v[501];
	v[518] = dAi[0] * v[494] + dAi[1] * v[498] + dAi[2] * v[501];
	v[545] = -v[518] / 2e0;
	v[563] = v[536] / 2e0 + v[545];
	v[554] = v[527] / 2e0 + v[545];
	v[492] = v[425] * v[477] + v[483];
	v[534] = dCi[0] * v[492] + dCi[1] * v[495] + dCi[2] * v[499];
	v[525] = dBi[0] * v[492] + dBi[1] * v[495] + dBi[2] * v[499];
	v[516] = dAi[0] * v[492] + dAi[1] * v[495] + dAi[2] * v[499];
	v[543] = -v[516] / 2e0;
	v[561] = v[534] / 2e0 + v[543];
	v[552] = v[525] / 2e0 + v[543];
	v[484] = v[429] * v[477] + v[483];
	v[531] = dCi[0] * v[480] + dCi[1] * v[484] + dCi[2] * v[489];
	v[522] = dBi[0] * v[480] + dBi[1] * v[484] + dBi[2] * v[489];
	v[513] = dAi[0] * v[480] + dAi[1] * v[484] + dAi[2] * v[489];
	v[540] = -v[513] / 2e0;
	v[558] = v[531] / 2e0 + v[540];
	v[549] = v[522] / 2e0 + v[540];
	v[482] = v[2788] * v[479] + v[497];
	v[533] = dCi[0] * v[482] + dCi[1] * v[487] + dCi[2] * v[491];
	v[524] = dBi[0] * v[482] + dBi[1] * v[487] + dBi[2] * v[491];
	v[515] = dAi[0] * v[482] + dAi[1] * v[487] + dAi[2] * v[491];
	v[542] = -v[515] / 2e0;
	v[560] = v[533] / 2e0 + v[542];
	v[551] = v[524] / 2e0 + v[542];
	v[481] = 0.5e0*v[432] * v[478] + v[510];
	v[532] = dCi[0] * v[481] + dCi[1] * v[486] + dCi[2] * v[490];
	v[523] = dBi[0] * v[481] + dBi[1] * v[486] + dBi[2] * v[490];
	v[514] = dAi[0] * v[481] + dAi[1] * v[486] + dAi[2] * v[490];
	v[541] = -v[514] / 2e0;
	v[559] = v[532] / 2e0 + v[541];
	v[550] = v[523] / 2e0 + v[541];
	v[242] = 1e0 + 0.5e0*v[239] * v[432];
	v[243] = v[239] * v[429];
	v[244] = v[239] * v[426];
	v[246] = v[239] * v[425];
	v[248] = 1e0 + 0.5e0*v[239] * v[423];
	v[249] = v[239] * v[420];
	v[251] = v[239] * v[419];
	v[253] = v[239] * v[418];
	v[254] = 1e0 + 0.5e0*v[239] * v[415];
	v[258] = 4e0 / v[1063];
	v[2468] = 2e0*v[1064] * v[258];
	v[2487] = -(v[1059] * v[1295]) - v[2464] + v[2468];
	v[1689] = (v[258] * v[258]);
	v[2472] = 1e0 / Power(v[1689], 2);
	v[2471] = -(v[2468] * v[2472]);
	v[1358] = v[1359] - v[258];
	v[1353] = v[1354] - v[258];
	v[1349] = v[1351] - v[258];
	v[1336] = -0.5e0*v[258];
	v[1337] = -v[1330] + v[1336];
	v[1321] = v[1336] + v[1685];
	v[1107] = v[1064] * v[1106] + v[258];
	v[1099] = v[1064] * v[1098] - v[258];
	v[1089] = -(v[1061] * v[258]);
	v[1091] = -v[1089] + v[1064] * v[1090];
	v[1086] = -(v[1059] * v[258]);
	v[1083] = v[1058] * v[258];
	v[1085] = v[1083] + v[1064] * v[1084];
	v[1081] = -0.5e0*v[1062] * v[258];
	v[1079] = -0.5e0*v[1060] * v[258];
	v[1072] = -0.5e0*v[1057] * v[258];
	v[1074] = v[1072] + 0.5e0*v[1064] * v[1073];
	v[1069] = 0.5e0*v[1064] * v[1068] + v[1072];
	v[1039] = 0.5e0*v[1068] * v[258];
	v[1037] = 0.5e0*v[1073] * v[258];
	v[1035] = 0.5e0*v[1077] * v[258];
	v[262] = v[1084] * v[258];
	v[263] = v[1090] * v[258];
	v[267] = 1e0 + v[1037];
	v[268] = v[1098] * v[258];
	v[272] = v[1106] * v[258];
	v[273] = 1e0 + v[1039];
	v[286] = dAi[0] + xPi[0];
	v[326] = -v[286] / 2e0;
	v[287] = dAi[1] + xPi[1];
	v[328] = -v[287] / 2e0;
	v[288] = dAi[2] + xPi[2];
	v[330] = -v[288] / 2e0;
	v[289] = dBi[0] + xPi[0];
	v[319] = v[289] / 2e0 + v[326];
	v[290] = dBi[1] + xPi[1];
	v[320] = v[290] / 2e0 + v[328];
	v[291] = dBi[2] + xPi[2];
	v[321] = v[291] / 2e0 + v[330];
	v[323] = 1e0 / sqrt((v[319] * v[319]) + (v[320] * v[320]) + (v[321] * v[321]));
	v[898] = v[323] * (v[845] / 2e0 + v[851]);
	v[897] = v[323] * (v[876] / 2e0 + v[878]);
	v[896] = v[323] * (v[846] / 2e0 + v[852]);
	v[895] = v[323] * (v[865] / 2e0 + v[869]);
	v[894] = v[323] * (v[864] / 2e0 + v[868]);
	v[893] = v[323] * (v[847] / 2e0 + v[853]);
	v[892] = v[323] * (v[836] / 2e0 + v[838]);
	v[891] = v[323] * (v[796] / 2e0 + v[802]);
	v[890] = v[323] * (v[795] / 2e0 + v[801]);
	v[889] = v[323] * (v[825] / 2e0 + v[829]);
	v[888] = v[323] * (v[797] / 2e0 + v[803]);
	v[887] = v[323] * (v[824] / 2e0 + v[828]);
	v[886] = v[323] * (v[768] / 2e0 + v[772]);
	v[885] = v[323] * (v[750] / 2e0 + v[756]);
	v[884] = v[323] * (v[749] / 2e0 + v[755]);
	v[883] = v[323] * (v[816] / 2e0 + v[818]);
	v[882] = v[323] * (v[751] / 2e0 + v[757]);
	v[881] = v[323] * (v[769] / 2e0 + v[773]);
	v[584] = v[323] * v[557];
	v[583] = v[323] * v[556];
	v[582] = v[323] * v[555];
	v[581] = v[323] * v[554];
	v[580] = v[323] * v[553];
	v[579] = v[323] * v[552];
	v[578] = v[323] * v[551];
	v[577] = v[323] * v[550];
	v[576] = v[323] * v[549];
	v[292] = dCi[0] + xPi[0];
	v[327] = v[292] / 2e0 + v[326];
	v[293] = dCi[1] + xPi[1];
	v[329] = v[293] / 2e0 + v[328];
	v[294] = dCi[2] + xPi[2];
	v[331] = v[294] / 2e0 + v[330];
	v[333] = 1e0 / sqrt((v[327] * v[327]) + (v[329] * v[329]) + (v[331] * v[331]));
	v[916] = v[333] * (v[842] / 2e0 + v[851]);
	v[915] = v[333] * (v[875] / 2e0 + v[878]);
	v[914] = v[333] * (v[843] / 2e0 + v[852]);
	v[913] = v[333] * (v[863] / 2e0 + v[869]);
	v[912] = v[333] * (v[862] / 2e0 + v[868]);
	v[911] = v[333] * (v[844] / 2e0 + v[853]);
	v[910] = v[333] * (v[835] / 2e0 + v[838]);
	v[909] = v[333] * (v[793] / 2e0 + v[802]);
	v[908] = v[333] * (v[792] / 2e0 + v[801]);
	v[907] = v[333] * (v[823] / 2e0 + v[829]);
	v[906] = v[333] * (v[794] / 2e0 + v[803]);
	v[905] = v[333] * (v[822] / 2e0 + v[828]);
	v[904] = v[333] * (v[766] / 2e0 + v[772]);
	v[903] = v[333] * (v[747] / 2e0 + v[756]);
	v[902] = v[333] * (v[746] / 2e0 + v[755]);
	v[901] = v[333] * (v[815] / 2e0 + v[818]);
	v[900] = v[333] * (v[748] / 2e0 + v[757]);
	v[899] = v[333] * (v[767] / 2e0 + v[773]);
	v[593] = v[333] * v[566];
	v[592] = v[333] * v[565];
	v[591] = v[333] * v[564];
	v[590] = v[333] * v[563];
	v[589] = v[333] * v[562];
	v[588] = v[333] * v[561];
	v[587] = v[333] * v[560];
	v[586] = v[333] * v[559];
	v[585] = v[333] * v[558];
	v[295] = v[224] + dAi[0] * v[242] + dAi[1] * v[243] + dAi[2] * v[244];
	v[349] = -v[295] / 2e0;
	v[296] = v[225] + dAi[0] * v[246] + dAi[1] * v[248] + dAi[2] * v[249];
	v[351] = -v[296] / 2e0;
	v[297] = v[226] + dAi[0] * v[251] + dAi[1] * v[253] + dAi[2] * v[254];
	v[353] = -v[297] / 2e0;
	v[298] = v[224] + dBi[0] * v[242] + dBi[1] * v[243] + dBi[2] * v[244];
	v[343] = v[298] / 2e0 + v[349];
	v[299] = v[225] + dBi[0] * v[246] + dBi[1] * v[248] + dBi[2] * v[249];
	v[344] = v[299] / 2e0 + v[351];
	v[300] = v[226] + dBi[0] * v[251] + dBi[1] * v[253] + dBi[2] * v[254];
	v[345] = v[300] / 2e0 + v[353];
	v[301] = v[224] + dCi[0] * v[242] + dCi[1] * v[243] + dCi[2] * v[244];
	v[350] = v[301] / 2e0 + v[349];
	v[302] = v[225] + dCi[0] * v[246] + dCi[1] * v[248] + dCi[2] * v[249];
	v[352] = v[302] / 2e0 + v[351];
	v[303] = v[226] + dCi[0] * v[251] + dCi[1] * v[253] + dCi[2] * v[254];
	v[354] = v[303] / 2e0 + v[353];
	v[304] = (-ci[0] - ci[1]) / 2e0;
	v[305] = (1e0 + ci[0]) / 2e0;
	v[306] = (1e0 + ci[1]) / 2e0;
	v[1381] = v[306] * v[842] + v[305] * v[845] + v[304] * v[848];
	v[1380] = v[306] * v[875] + v[305] * v[876] + v[304] * v[877];
	v[1379] = v[306] * v[843] + v[305] * v[846] + v[304] * v[849];
	v[1378] = v[306] * v[863] + v[305] * v[865] + v[304] * v[867];
	v[1377] = v[306] * v[862] + v[305] * v[864] + v[304] * v[866];
	v[1376] = v[306] * v[844] + v[305] * v[847] + v[304] * v[850];
	v[1375] = v[306] * v[835] + v[305] * v[836] + v[304] * v[837];
	v[1374] = v[306] * v[793] + v[305] * v[796] + v[304] * v[799];
	v[1373] = v[306] * v[792] + v[305] * v[795] + v[304] * v[798];
	v[1372] = v[306] * v[823] + v[305] * v[825] + v[304] * v[827];
	v[1371] = v[306] * v[794] + v[305] * v[797] + v[304] * v[800];
	v[1370] = v[306] * v[822] + v[305] * v[824] + v[304] * v[826];
	v[1369] = v[306] * v[766] + v[305] * v[768] + v[304] * v[770];
	v[1368] = v[306] * v[747] + v[305] * v[750] + v[304] * v[753];
	v[1367] = v[306] * v[746] + v[305] * v[749] + v[304] * v[752];
	v[1366] = v[306] * v[815] + v[305] * v[816] + v[304] * v[817];
	v[1365] = v[306] * v[748] + v[305] * v[751] + v[304] * v[754];
	v[1364] = v[306] * v[767] + v[305] * v[769] + v[304] * v[771];
	v[1118] = v[304] * v[521] + v[305] * v[530] + v[306] * v[539];
	v[1117] = v[304] * v[520] + v[305] * v[529] + v[306] * v[538];
	v[1116] = v[304] * v[519] + v[305] * v[528] + v[306] * v[537];
	v[1115] = v[304] * v[518] + v[305] * v[527] + v[306] * v[536];
	v[1114] = v[304] * v[517] + v[305] * v[526] + v[306] * v[535];
	v[1113] = v[304] * v[516] + v[305] * v[525] + v[306] * v[534];
	v[1112] = v[304] * v[515] + v[305] * v[524] + v[306] * v[533];
	v[1111] = v[304] * v[514] + v[305] * v[523] + v[306] * v[532];
	v[1110] = v[304] * v[513] + v[305] * v[522] + v[306] * v[531];
	v[307] = (-cp[0] - cp[1]) / 2e0;
	v[308] = (1e0 + cp[0]) / 2e0;
	v[1804] = -(v[308] * v[587]);
	v[1803] = -(v[308] * v[586]);
	v[1802] = -(v[308] * v[585]);
	v[1792] = -(v[308] * v[590]);
	v[1791] = -(v[308] * v[589]);
	v[1790] = -(v[308] * v[588]);
	v[1780] = -(v[308] * v[593]);
	v[1779] = -(v[308] * v[592]);
	v[1778] = -(v[308] * v[591]);
	v[309] = (1e0 + cp[1]) / 2e0;
	v[1795] = -(v[309] * v[578]);
	v[1794] = -(v[309] * v[577]);
	v[1793] = -(v[309] * v[576]);
	v[1783] = -(v[309] * v[581]);
	v[1782] = -(v[309] * v[580]);
	v[1781] = -(v[309] * v[579]);
	v[1771] = -(v[309] * v[584]);
	v[1770] = -(v[309] * v[583]);
	v[1769] = -(v[309] * v[582]);
	v[575] = v[307] * v[521] + v[308] * v[530] + v[309] * v[539];
	v[1768] = -(v[323] * v[575]) / 2e0;
	v[1807] = -v[1768] - v[307] * v[584];
	v[1777] = v[1768] - v[308] * v[584];
	v[1833] = dCi[0] * v[1771] + dBi[0] * v[1777] + dAi[0] * v[1807];
	v[1819] = dCi[1] * v[1771] + dBi[1] * v[1777] + dAi[1] * v[1807];
	v[1813] = dCi[2] * v[1771] + dBi[2] * v[1777] + dAi[2] * v[1807];
	v[1765] = -(v[333] * v[575]) / 2e0;
	v[1810] = -v[1765] - v[307] * v[593];
	v[1774] = v[1765] - v[309] * v[593];
	v[1840] = dCi[0] * v[1774] + dBi[0] * v[1780] + dAi[0] * v[1810];
	v[1826] = dCi[1] * v[1774] + dBi[1] * v[1780] + dAi[1] * v[1810];
	v[1816] = dCi[2] * v[1774] + dBi[2] * v[1780] + dAi[2] * v[1810];
	v[574] = v[307] * v[520] + v[308] * v[529] + v[309] * v[538];
	v[1767] = -(v[323] * v[574]) / 2e0;
	v[1806] = -v[1767] - v[307] * v[583];
	v[1776] = v[1767] - v[308] * v[583];
	v[1832] = dCi[0] * v[1770] + dBi[0] * v[1776] + dAi[0] * v[1806];
	v[1818] = dCi[1] * v[1770] + dBi[1] * v[1776] + dAi[1] * v[1806];
	v[1812] = dCi[2] * v[1770] + dBi[2] * v[1776] + dAi[2] * v[1806];
	v[1764] = -(v[333] * v[574]) / 2e0;
	v[1809] = -v[1764] - v[307] * v[592];
	v[1773] = v[1764] - v[309] * v[592];
	v[1839] = dCi[0] * v[1773] + dBi[0] * v[1779] + dAi[0] * v[1809];
	v[1825] = dCi[1] * v[1773] + dBi[1] * v[1779] + dAi[1] * v[1809];
	v[1815] = dCi[2] * v[1773] + dBi[2] * v[1779] + dAi[2] * v[1809];
	v[573] = v[307] * v[519] + v[308] * v[528] + v[309] * v[537];
	v[1766] = -(v[323] * v[573]) / 2e0;
	v[1805] = -v[1766] - v[307] * v[582];
	v[1775] = v[1766] - v[308] * v[582];
	v[1831] = dCi[0] * v[1769] + dBi[0] * v[1775] + dAi[0] * v[1805];
	v[1817] = dCi[1] * v[1769] + dBi[1] * v[1775] + dAi[1] * v[1805];
	v[1811] = dCi[2] * v[1769] + dBi[2] * v[1775] + dAi[2] * v[1805];
	v[1763] = -(v[333] * v[573]) / 2e0;
	v[1808] = -v[1763] - v[307] * v[591];
	v[1772] = v[1763] - v[309] * v[591];
	v[1838] = dCi[0] * v[1772] + dBi[0] * v[1778] + dAi[0] * v[1808];
	v[1824] = dCi[1] * v[1772] + dBi[1] * v[1778] + dAi[1] * v[1808];
	v[1814] = dCi[2] * v[1772] + dBi[2] * v[1778] + dAi[2] * v[1808];
	v[572] = v[307] * v[518] + v[308] * v[527] + v[309] * v[536];
	v[1762] = -(v[323] * v[572]) / 2e0;
	v[1847] = -v[1762] - v[307] * v[581];
	v[1789] = v[1762] - v[308] * v[581];
	v[1873] = dCi[0] * v[1783] + dBi[0] * v[1789] + dAi[0] * v[1847];
	v[1867] = dCi[1] * v[1783] + dBi[1] * v[1789] + dAi[1] * v[1847];
	v[1853] = dCi[2] * v[1783] + dBi[2] * v[1789] + dAi[2] * v[1847];
	v[1759] = -(v[333] * v[572]) / 2e0;
	v[1850] = -v[1759] - v[307] * v[590];
	v[1786] = v[1759] - v[309] * v[590];
	v[1879] = dCi[0] * v[1786] + dBi[0] * v[1792] + dAi[0] * v[1850];
	v[1870] = dCi[1] * v[1786] + dBi[1] * v[1792] + dAi[1] * v[1850];
	v[1860] = dCi[2] * v[1786] + dBi[2] * v[1792] + dAi[2] * v[1850];
	v[571] = v[307] * v[517] + v[308] * v[526] + v[309] * v[535];
	v[1761] = -(v[323] * v[571]) / 2e0;
	v[1846] = -v[1761] - v[307] * v[580];
	v[1788] = v[1761] - v[308] * v[580];
	v[1872] = dCi[0] * v[1782] + dBi[0] * v[1788] + dAi[0] * v[1846];
	v[1866] = dCi[1] * v[1782] + dBi[1] * v[1788] + dAi[1] * v[1846];
	v[1852] = dCi[2] * v[1782] + dBi[2] * v[1788] + dAi[2] * v[1846];
	v[1758] = -(v[333] * v[571]) / 2e0;
	v[1849] = -v[1758] - v[307] * v[589];
	v[1785] = v[1758] - v[309] * v[589];
	v[1878] = dCi[0] * v[1785] + dBi[0] * v[1791] + dAi[0] * v[1849];
	v[1869] = dCi[1] * v[1785] + dBi[1] * v[1791] + dAi[1] * v[1849];
	v[1859] = dCi[2] * v[1785] + dBi[2] * v[1791] + dAi[2] * v[1849];
	v[570] = v[307] * v[516] + v[308] * v[525] + v[309] * v[534];
	v[1760] = -(v[323] * v[570]) / 2e0;
	v[1845] = -v[1760] - v[307] * v[579];
	v[1787] = v[1760] - v[308] * v[579];
	v[1871] = dCi[0] * v[1781] + dBi[0] * v[1787] + dAi[0] * v[1845];
	v[1865] = dCi[1] * v[1781] + dBi[1] * v[1787] + dAi[1] * v[1845];
	v[1851] = dCi[2] * v[1781] + dBi[2] * v[1787] + dAi[2] * v[1845];
	v[1757] = -(v[333] * v[570]) / 2e0;
	v[1848] = -v[1757] - v[307] * v[588];
	v[1784] = v[1757] - v[309] * v[588];
	v[1877] = dCi[0] * v[1784] + dBi[0] * v[1790] + dAi[0] * v[1848];
	v[1868] = dCi[1] * v[1784] + dBi[1] * v[1790] + dAi[1] * v[1848];
	v[1858] = dCi[2] * v[1784] + dBi[2] * v[1790] + dAi[2] * v[1848];
	v[569] = v[307] * v[515] + v[308] * v[524] + v[309] * v[533];
	v[1756] = -(v[323] * v[569]) / 2e0;
	v[1885] = -v[1756] - v[307] * v[578];
	v[1801] = v[1756] - v[308] * v[578];
	v[1915] = dCi[0] * v[1795] + dBi[0] * v[1801] + dAi[0] * v[1885];
	v[1903] = dCi[1] * v[1795] + dBi[1] * v[1801] + dAi[1] * v[1885];
	v[1891] = dCi[2] * v[1795] + dBi[2] * v[1801] + dAi[2] * v[1885];
	v[1752] = -(v[333] * v[569]) / 2e0;
	v[1888] = -v[1752] - v[307] * v[587];
	v[1798] = v[1752] - v[309] * v[587];
	v[1918] = dCi[0] * v[1798] + dBi[0] * v[1804] + dAi[0] * v[1888];
	v[1909] = dCi[1] * v[1798] + dBi[1] * v[1804] + dAi[1] * v[1888];
	v[1897] = dCi[2] * v[1798] + dBi[2] * v[1804] + dAi[2] * v[1888];
	v[568] = v[307] * v[514] + v[308] * v[523] + v[309] * v[532];
	v[1755] = -(v[323] * v[568]) / 2e0;
	v[1884] = -v[1755] - v[307] * v[577];
	v[1800] = v[1755] - v[308] * v[577];
	v[1914] = dCi[0] * v[1794] + dBi[0] * v[1800] + dAi[0] * v[1884];
	v[1902] = dCi[1] * v[1794] + dBi[1] * v[1800] + dAi[1] * v[1884];
	v[1890] = dCi[2] * v[1794] + dBi[2] * v[1800] + dAi[2] * v[1884];
	v[1751] = -(v[333] * v[568]) / 2e0;
	v[1887] = -v[1751] - v[307] * v[586];
	v[1797] = v[1751] - v[309] * v[586];
	v[1917] = dCi[0] * v[1797] + dBi[0] * v[1803] + dAi[0] * v[1887];
	v[1908] = dCi[1] * v[1797] + dBi[1] * v[1803] + dAi[1] * v[1887];
	v[1896] = dCi[2] * v[1797] + dBi[2] * v[1803] + dAi[2] * v[1887];
	v[567] = v[307] * v[513] + v[308] * v[522] + v[309] * v[531];
	v[1754] = -(v[323] * v[567]) / 2e0;
	v[1883] = -v[1754] - v[307] * v[576];
	v[1799] = v[1754] - v[308] * v[576];
	v[1901] = dCi[1] * v[1793] + dBi[1] * v[1799] + dAi[1] * v[1883];
	v[1889] = dCi[2] * v[1793] + dBi[2] * v[1799] + dAi[2] * v[1883];
	v[1750] = -(v[333] * v[567]) / 2e0;
	v[1886] = -v[1750] - v[307] * v[585];
	v[1796] = v[1750] - v[309] * v[585];
	v[1907] = dCi[1] * v[1796] + dBi[1] * v[1802] + dAi[1] * v[1886];
	v[1895] = dCi[2] * v[1796] + dBi[2] * v[1802] + dAi[2] * v[1886];
	v[365] = d[0] - v[295] * v[307] - v[298] * v[308] - v[301] * v[309] + xSi[0];
	v[652] = -v[365] / 2e0;
	v[407] = (v[333] * v[365]) / 2e0;
	v[405] = (v[323] * v[365]) / 2e0;
	v[366] = d[1] - v[296] * v[307] - v[299] * v[308] - v[302] * v[309] + xSi[1];
	v[649] = -v[366] / 2e0;
	v[397] = (v[333] * v[366]) / 2e0;
	v[395] = (v[323] * v[366]) / 2e0;
	v[367] = d[2] - v[297] * v[307] - v[300] * v[308] - v[303] * v[309] + xSi[2];
	v[646] = -v[367] / 2e0;
	v[387] = (v[333] * v[367]) / 2e0;
	v[385] = (v[323] * v[367]) / 2e0;
	v[322] = v[319] * v[323];
	v[324] = v[320] * v[323];
	v[325] = v[321] * v[323];
	v[1387] = v[325] * v[886] + v[324] * v[892] + v[322] * v[898];
	v[1386] = v[325] * v[885] + v[324] * v[891] + v[322] * v[897];
	v[1385] = v[325] * v[884] + v[324] * v[890] + v[322] * v[896];
	v[1384] = v[325] * v[883] + v[324] * v[889] + v[322] * v[895];
	v[1383] = v[325] * v[882] + v[324] * v[888] + v[322] * v[894];
	v[1382] = v[325] * v[881] + v[324] * v[887] + v[322] * v[893];
	v[1130] = v[322] * v[578] + v[324] * v[581] + v[325] * v[584];
	v[1129] = v[322] * v[577] + v[324] * v[580] + v[325] * v[583];
	v[1128] = v[322] * v[576] + v[324] * v[579] + v[325] * v[582];
	v[332] = v[327] * v[333];
	v[334] = v[329] * v[333];
	v[341] = -(v[324] * v[332]) + v[322] * v[334];
	v[335] = v[331] * v[333];
	v[339] = v[325] * v[332] - v[322] * v[335];
	v[336] = -(v[325] * v[334]) + v[324] * v[335];
	v[338] = 1e0 / sqrt((v[336] * v[336]) + (v[339] * v[339]) + (v[341] * v[341]));
	v[337] = v[336] * v[338];
	v[340] = v[338] * v[339];
	v[342] = v[338] * v[341];
	v[1411] = v[342] * v[886] + v[340] * v[892] + v[337] * v[898];
	v[1410] = v[342] * v[885] + v[340] * v[891] + v[337] * v[897];
	v[1409] = v[342] * v[884] + v[340] * v[890] + v[337] * v[896];
	v[1408] = v[342] * v[883] + v[340] * v[889] + v[337] * v[895];
	v[1407] = v[342] * v[882] + v[340] * v[888] + v[337] * v[894];
	v[1406] = v[342] * v[881] + v[340] * v[887] + v[337] * v[893];
	v[1405] = (*radius)*(-0.5e0*v[1077] * v[1295] * v[337] + (-(v[1084] * v[1295]) + v[1297])*v[340] + (-
		(v[1090] * v[1295]) + v[1299])*v[342]);
	v[1404] = (*radius)*((-0.5e0*v[1077] * v[1290] - v[1353] - v[1354])*v[337] + (-(v[1084] * v[1290]) + v[1297]
		)*v[340] + (-(v[1090] * v[1290]) + v[1318])*v[342]);
	v[1403] = (*radius)*(-((0.5e0*v[1077] * v[1288] + v[1297])*v[337]) - (v[1084] * v[1288] + v[1327] - v[1337]
		)*v[340] - (v[1064] + v[1090] * v[1288] + v[1313])*v[342]);
	v[1402] = (*radius)*((-0.5e0*v[1077] * v[1284] - v[1349] - v[1351])*v[337] + (-(v[1084] * v[1284]) + v[1334]
		)*v[340] + (-(v[1090] * v[1284]) + v[1299])*v[342]);
	v[1401] = (*radius)*((-2e0*v[1293] - v[1283] * v[2795])*v[337] + (-(v[1084] * v[1283]) + v[2790])*v[340] -
		(v[1067] + v[1090] * v[1283] + v[1306])*v[342]);
	v[1400] = (*radius)*(-((0.5e0*v[1077] * v[1282] + v[1299])*v[337]) + (-(v[1084] * v[1282]) + v[2791])*v[340] -
		(v[1090] * v[1282] - v[1321] + v[1330])*v[342]);
	v[1399] = (*radius)*((-(v[1094] * v[1295]) + v[1297])*v[337] + (-v[1358] - v[1359] + v[1295] * v[2789])*v[340] + (-
		(v[1098] * v[1295]) + v[1308])*v[342]);
	v[1398] = (*radius)*((-(v[1094] * v[1290]) + v[1297])*v[337] + v[1290] * v[2789] * v[340] + (-(v[1098] * v[1290])
		+ v[1293])*v[342]);
	v[1397] = (*radius)*(-((v[1094] * v[1288] + v[1327] - v[1337])*v[337]) - (0.5e0*v[1073] * v[1288] + v[1297]
		)*v[340] + (-(v[1098] * v[1288]) + v[2790])*v[342]);
	v[1396] = (*radius)*(-((v[1094] * v[1284] + v[1334])*v[337]) - (0.5e0*v[1073] * v[1284] + v[1349] + v[1351]
		)*v[340] - (v[1098] * v[1284] - v[1293])*v[342]);
	v[1395] = (*radius)*(-((v[1066] + v[1094] * v[1283] + v[1310])*v[337]) - (0.5e0*v[1073] * v[1283] + v[1293]
		)*v[340] - (v[1098] * v[1283] - v[1321] + v[1327])*v[342]);
	v[1394] = (*radius)*(-((v[1064] + v[1094] * v[1282] + v[1313])*v[337]) + (-0.5e0*v[1073] * v[1282] - 2e0*v[1299]
		)*v[340] + (-(v[1098] * v[1282]) + v[2793])*v[342]);
	v[1393] = (*radius)*(-((v[1102] * v[1295] - v[1299])*v[337]) - (v[1106] * v[1295] + v[1308])*v[340] -
		(0.5e0*v[1068] * v[1295] + v[1358] + v[1359])*v[342]);
	v[1392] = (*radius)*(-((v[1102] * v[1290] + v[1318])*v[337]) - (v[1106] * v[1290] - v[1293])*v[340] -
		(0.5e0*v[1068] * v[1290] + v[1353] + v[1354])*v[342]);
	v[1391] = (*radius)*((-(v[1102] * v[1288]) + v[2791])*v[337] - (v[1066] + v[1106] * v[1288] + v[1310])*v[340] + (
		-2e0*v[1297] + v[1288] * v[2792])*v[342]);
	v[1390] = (*radius)*((-(v[1102] * v[1284]) + v[1299])*v[337] + (-(v[1106] * v[1284]) + v[1293])*v[340]
		+ v[1284] * v[2792] * v[342]);
	v[1389] = (*radius)*((-(v[1102] * v[1283]) + v[2793])*v[337] - (v[1106] * v[1283] - v[1321] + v[1327])*v[340] -
		(0.5e0*v[1068] * v[1283] + v[1293])*v[342]);
	v[1388] = (*radius)*(-((v[1102] * v[1282] - v[1321] + v[1330])*v[337]) - (v[1067] + v[1106] * v[1282] + v[1306]
		)*v[340] - (v[1299] + v[1282] * v[2794])*v[342]);
	v[1163] = (*radius)*((v[1086] - v[1067] * v[1102])*v[337] - (v[1083] + v[1067] * v[1106])*v[340]
		- v[1067] * v[2794] * v[342]);
	v[1162] = (*radius)*((-(v[1066] * v[1102]) + v[258])*v[337] - (-v[1089] + v[1066] * v[1106])*v[340] -
		(0.5e0*v[1066] * v[1068] + v[1079])*v[342]);
	v[1161] = -((*radius)*((-v[1089] + v[1064] * v[1102])*v[337] + v[1107] * v[340] + v[1069] * v[342]));
	v[1160] = (*radius)*(-((v[1067] * v[1094] + v[258])*v[337]) - (0.5e0*v[1067] * v[1073] + v[1081])*v[340] -
		(v[1083] + v[1067] * v[1098])*v[342]);
	v[1159] = (*radius)*((v[1086] - v[1066] * v[1094])*v[337] - 0.5e0*v[1066] * v[1073] * v[340] + (v[1089]
		- v[1066] * v[1098])*v[342]);
	v[1158] = -((*radius)*((v[1083] + v[1064] * v[1094])*v[337] + v[1074] * v[340] + v[1099] * v[342]));
	v[1157] = (*radius)*(-((v[1081] + v[1067] * v[2795])*v[337]) + (-(v[1067] * v[1084]) + v[258])*v[340] + (v[1086]
		- v[1067] * v[1090])*v[342]);
	v[1156] = (*radius)*(-((v[1079] + v[1066] * v[2796])*v[337]) - (v[1066] * v[1084] - v[1086])*v[340] -
		(v[1066] * v[1090] + v[258])*v[342]);
	v[1155] = (*radius)*(-(v[1064] * v[2796] * v[337]) - v[1085] * v[340] - v[1091] * v[342]);
	v[1148] = v[337] * v[578] + v[340] * v[581] + v[342] * v[584];
	v[1147] = v[337] * v[577] + v[340] * v[580] + v[342] * v[583];
	v[1146] = v[337] * v[576] + v[340] * v[579] + v[342] * v[582];
	v[346] = v[323] * v[343];
	v[347] = v[323] * v[344];
	v[348] = v[323] * v[345];
	v[370] = -(v[346] * v[350]) - v[347] * v[352] - v[348] * v[354];
	v[369] = -(v[343] * v[346]) - v[344] * v[347] - v[345] * v[348];
	v[355] = v[333] * v[350];
	v[356] = v[333] * v[352];
	v[945] = v[338] * (-2e0*v[579] * v[585] + 2e0*v[576] * v[588] - v[355] * v[892] + v[356] * v[898] + v[346] * v[910]
		- v[347] * v[916]);
	v[946] = -(v[309] * v[766]) - v[308] * v[768] - v[307] * v[770] - (*radius)*v[945];
	v[942] = v[338] * (-2e0*v[580] * v[586] + 2e0*v[577] * v[589] - v[355] * v[891] + v[356] * v[897] + v[346] * v[909]
		- v[347] * v[915]);
	v[944] = -(v[309] * v[747]) - v[308] * v[750] - v[307] * v[753] - (*radius)*v[942];
	v[941] = v[338] * (-(v[580] * v[585]) - v[579] * v[586] + v[577] * v[588] + v[576] * v[589] - v[355] * v[890]
		+ v[356] * v[896] + v[346] * v[908] - v[347] * v[914]);
	v[943] = -(v[309] * v[746]) - v[308] * v[749] - v[307] * v[752] - (*radius)*v[941];
	v[937] = v[338] * (-2e0*v[581] * v[587] + 2e0*v[578] * v[590] - v[355] * v[889] + v[356] * v[895] + v[346] * v[907]
		- v[347] * v[913]);
	v[940] = -(v[309] * v[815]) - v[308] * v[816] - v[307] * v[817] - (*radius)*v[937];
	v[936] = v[338] * (-(v[581] * v[586]) - v[580] * v[587] + v[578] * v[589] + v[577] * v[590] - v[355] * v[888]
		+ v[356] * v[894] + v[346] * v[906] - v[347] * v[912]);
	v[939] = -(v[309] * v[748]) - v[308] * v[751] - v[307] * v[754] - (*radius)*v[936];
	v[935] = v[338] * (-(v[581] * v[585]) - v[579] * v[587] + v[578] * v[588] + v[576] * v[590] - v[355] * v[887]
		+ v[356] * v[893] + v[346] * v[905] - v[347] * v[911]);
	v[938] = -(v[309] * v[767]) - v[308] * v[769] - v[307] * v[771] - (*radius)*v[935];
	v[602] = v[338] * (v[356] * v[578] - v[355] * v[581] - v[347] * v[587] + v[346] * v[590]);
	v[2800] = 2e0*v[602];
	v[611] = -v[575] - (*radius)*v[602];
	v[601] = v[338] * (v[356] * v[577] - v[355] * v[580] - v[347] * v[586] + v[346] * v[589]);
	v[2801] = 2e0*v[601];
	v[610] = -v[574] - (*radius)*v[601];
	v[600] = v[338] * (v[356] * v[576] - v[355] * v[579] - v[347] * v[585] + v[346] * v[588]);
	v[2802] = 2e0*v[600];
	v[609] = -v[573] - (*radius)*v[600];
	v[357] = v[333] * v[354];
	v[969] = v[338] * (-2e0*v[582] * v[588] + 2e0*v[579] * v[591] - v[356] * v[886] + v[357] * v[892] + v[347] * v[904]
		- v[348] * v[910]);
	v[970] = -(v[309] * v[842]) - v[308] * v[845] - v[307] * v[848] - (*radius)*v[969];
	v[966] = v[338] * (-2e0*v[583] * v[589] + 2e0*v[580] * v[592] - v[356] * v[885] + v[357] * v[891] + v[347] * v[903]
		- v[348] * v[909]);
	v[968] = -(v[309] * v[875]) - v[308] * v[876] - v[307] * v[877] - (*radius)*v[966];
	v[965] = v[338] * (-(v[583] * v[588]) - v[582] * v[589] + v[580] * v[591] + v[579] * v[592] - v[356] * v[884]
		+ v[357] * v[890] + v[347] * v[902] - v[348] * v[908]);
	v[967] = -(v[309] * v[843]) - v[308] * v[846] - v[307] * v[849] - (*radius)*v[965];
	v[961] = v[338] * (-2e0*v[584] * v[590] + 2e0*v[581] * v[593] - v[356] * v[883] + v[357] * v[889] + v[347] * v[901]
		- v[348] * v[907]);
	v[964] = -(v[309] * v[863]) - v[308] * v[865] - v[307] * v[867] - (*radius)*v[961];
	v[960] = v[338] * (-(v[584] * v[589]) - v[583] * v[590] + v[581] * v[592] + v[580] * v[593] - v[356] * v[882]
		+ v[357] * v[888] + v[347] * v[900] - v[348] * v[906]);
	v[963] = -(v[309] * v[862]) - v[308] * v[864] - v[307] * v[866] - (*radius)*v[960];
	v[959] = v[338] * (-(v[584] * v[588]) - v[582] * v[590] + v[581] * v[591] + v[579] * v[593] - v[356] * v[881]
		+ v[357] * v[887] + v[347] * v[899] - v[348] * v[905]);
	v[962] = -(v[309] * v[844]) - v[308] * v[847] - v[307] * v[850] - (*radius)*v[959];
	v[957] = v[338] * (2e0*v[582] * v[585] - 2e0*v[576] * v[591] + v[355] * v[886] - v[357] * v[898] - v[346] * v[904]
		+ v[348] * v[916]);
	v[1423] = v[325] * v[945] + v[324] * v[957] + v[322] * v[969];
	v[1422] = v[342] * v[945] + v[340] * v[957] + v[337] * v[969];
	v[958] = -(v[309] * v[835]) - v[308] * v[836] - v[307] * v[837] - (*radius)*v[957];
	v[954] = v[338] * (2e0*v[583] * v[586] - 2e0*v[577] * v[592] + v[355] * v[885] - v[357] * v[897] - v[346] * v[903]
		+ v[348] * v[915]);
	v[1421] = v[325] * v[942] + v[324] * v[954] + v[322] * v[966];
	v[1419] = v[342] * v[942] + v[340] * v[954] + v[337] * v[966];
	v[956] = -(v[309] * v[793]) - v[308] * v[796] - v[307] * v[799] - (*radius)*v[954];
	v[953] = v[338] * (v[583] * v[585] + v[582] * v[586] - v[577] * v[591] - v[576] * v[592] + v[355] * v[884] - v[357] * v[896]
		- v[346] * v[902] + v[348] * v[914]);
	v[1420] = v[325] * v[941] + v[324] * v[953] + v[322] * v[965];
	v[1418] = v[342] * v[941] + v[340] * v[953] + v[337] * v[965];
	v[955] = -(v[309] * v[792]) - v[308] * v[795] - v[307] * v[798] - (*radius)*v[953];
	v[949] = v[338] * (2e0*v[584] * v[587] - 2e0*v[578] * v[593] + v[355] * v[883] - v[357] * v[895] - v[346] * v[901]
		+ v[348] * v[913]);
	v[1417] = v[325] * v[937] + v[324] * v[949] + v[322] * v[961];
	v[1414] = v[342] * v[937] + v[340] * v[949] + v[337] * v[961];
	v[952] = -(v[309] * v[823]) - v[308] * v[825] - v[307] * v[827] - (*radius)*v[949];
	v[948] = v[338] * (v[584] * v[586] + v[583] * v[587] - v[578] * v[592] - v[577] * v[593] + v[355] * v[882] - v[357] * v[894]
		- v[346] * v[900] + v[348] * v[912]);
	v[1416] = v[325] * v[936] + v[324] * v[948] + v[322] * v[960];
	v[1413] = v[342] * v[936] + v[340] * v[948] + v[337] * v[960];
	v[951] = -(v[309] * v[794]) - v[308] * v[797] - v[307] * v[800] - (*radius)*v[948];
	v[947] = v[338] * (v[584] * v[585] + v[582] * v[587] - v[578] * v[591] - v[576] * v[593] + v[355] * v[881] - v[357] * v[893]
		- v[346] * v[899] + v[348] * v[911]);
	v[1415] = v[325] * v[935] + v[324] * v[947] + v[322] * v[959];
	v[1412] = v[342] * v[935] + v[340] * v[947] + v[337] * v[959];
	v[950] = -(v[309] * v[822]) - v[308] * v[824] - v[307] * v[826] - (*radius)*v[947];
	v[599] = v[338] * (-(v[357] * v[578]) + v[355] * v[584] + v[348] * v[587] - v[346] * v[593]);
	v[2797] = -2e0*v[599];
	v[608] = -v[572] - (*radius)*v[599];
	v[598] = v[338] * (-(v[357] * v[577]) + v[355] * v[583] + v[348] * v[586] - v[346] * v[592]);
	v[2798] = -2e0*v[598];
	v[607] = -v[571] - (*radius)*v[598];
	v[597] = v[338] * (-(v[357] * v[576]) + v[355] * v[582] + v[348] * v[585] - v[346] * v[591]);
	v[2799] = -2e0*v[597];
	v[606] = -v[570] - (*radius)*v[597];
	v[596] = v[338] * (v[357] * v[581] - v[356] * v[584] - v[348] * v[590] + v[347] * v[593]);
	v[1495] = v[1155] * v[596] + v[1158] * v[599] + v[1161] * v[602];
	v[1476] = v[1156] * v[596] + v[1159] * v[599] + v[1162] * v[602];
	v[1453] = v[1157] * v[596] + v[1160] * v[599] + v[1163] * v[602];
	v[1154] = v[337] * v[596] + v[340] * v[599] + v[342] * v[602];
	v[1136] = v[322] * v[596] + v[324] * v[599] + v[325] * v[602];
	v[605] = -v[569] - (*radius)*v[596];
	v[595] = v[338] * (v[357] * v[580] - v[356] * v[583] - v[348] * v[589] + v[347] * v[592]);
	v[1494] = v[1155] * v[595] + v[1158] * v[598] + v[1161] * v[601];
	v[1475] = v[1156] * v[595] + v[1159] * v[598] + v[1162] * v[601];
	v[1452] = v[1157] * v[595] + v[1160] * v[598] + v[1163] * v[601];
	v[1153] = v[337] * v[595] + v[340] * v[598] + v[342] * v[601];
	v[1135] = v[322] * v[595] + v[324] * v[598] + v[325] * v[601];
	v[604] = -v[568] - (*radius)*v[595];
	v[594] = v[338] * (v[357] * v[579] - v[356] * v[582] - v[348] * v[588] + v[347] * v[591]);
	v[1493] = v[1155] * v[594] + v[1158] * v[597] + v[1161] * v[600];
	v[1474] = v[1156] * v[594] + v[1159] * v[597] + v[1162] * v[600];
	v[1451] = v[1157] * v[594] + v[1160] * v[597] + v[1163] * v[600];
	v[1152] = v[337] * v[594] + v[340] * v[597] + v[342] * v[600];
	v[1134] = v[322] * v[594] + v[324] * v[597] + v[325] * v[600];
	v[603] = -v[567] - (*radius)*v[594];
	v[372] = -(v[350] * v[355]) - v[352] * v[356] - v[354] * v[357];
	v[371] = -(v[343] * v[355]) - v[344] * v[356] - v[345] * v[357];
	v[358] = v[338] * (-(v[348] * v[356]) + v[347] * v[357]);
	v[1426] = -2e0*v[358] * v[596];
	v[1425] = -2e0*v[358] * v[595];
	v[1424] = -2e0*v[358] * v[594];
	v[1179] = 1e0 - (v[358] * v[358]);
	v[359] = v[338] * (v[348] * v[355] - v[346] * v[357]);
	v[1438] = v[2799] * v[576] + 2e0*v[579] * v[594] + v[358] * v[892] - v[359] * v[898] - v[346] * v[957] + v[347] * v[969];
	v[1437] = v[2798] * v[577] + 2e0*v[580] * v[595] + v[358] * v[891] - v[359] * v[897] - v[346] * v[954] + v[347] * v[966];
	v[1436] = v[580] * v[594] + v[579] * v[595] - v[577] * v[597] - v[576] * v[598] + v[358] * v[890] - v[359] * v[896]
		- v[346] * v[953] + v[347] * v[965];
	v[1435] = v[2797] * v[578] + 2e0*v[581] * v[596] + v[358] * v[889] - v[359] * v[895] - v[346] * v[949] + v[347] * v[961];
	v[1434] = v[581] * v[595] + v[580] * v[596] - v[578] * v[598] - v[577] * v[599] + v[358] * v[888] - v[359] * v[894]
		- v[346] * v[948] + v[347] * v[960];
	v[1433] = v[581] * v[594] + v[579] * v[596] - v[578] * v[597] - v[576] * v[599] + v[358] * v[887] - v[359] * v[893]
		- v[346] * v[947] + v[347] * v[959];
	v[1432] = -(v[359] * v[596]) - v[358] * v[599];
	v[1431] = -(v[359] * v[595]) - v[358] * v[598];
	v[1430] = -(v[359] * v[594]) - v[358] * v[597];
	v[1429] = v[2797] * v[359];
	v[1428] = v[2798] * v[359];
	v[1427] = v[2799] * v[359];
	v[1197] = 1e0 - (v[359] * v[359]);
	v[1180] = -(v[358] * v[359]);
	v[1239] = (v[1180] * v[1180]);
	v[1127] = -(v[359] * v[578]) + v[358] * v[581] + v[347] * v[596] - v[346] * v[599];
	v[1126] = -(v[359] * v[577]) + v[358] * v[580] + v[347] * v[595] - v[346] * v[598];
	v[1125] = -(v[359] * v[576]) + v[358] * v[579] + v[347] * v[594] - v[346] * v[597];
	v[360] = v[338] * (-(v[347] * v[355]) + v[346] * v[356]);
	v[1529] = 2e0*v[582] * v[597] - 2e0*v[579] * v[600] + v[359] * v[886] - v[360] * v[892] - v[347] * v[945]
		+ v[348] * v[957];
	v[1524] = 2e0*v[583] * v[598] - 2e0*v[580] * v[601] + v[359] * v[885] - v[360] * v[891] - v[347] * v[942]
		+ v[348] * v[954];
	v[1523] = v[583] * v[597] + v[582] * v[598] - v[580] * v[600] - v[579] * v[601] + v[359] * v[884] - v[360] * v[890]
		- v[347] * v[941] + v[348] * v[953];
	v[1516] = 2e0*v[584] * v[599] - 2e0*v[581] * v[602] + v[359] * v[883] - v[360] * v[889] - v[347] * v[937]
		+ v[348] * v[949];
	v[1515] = v[584] * v[598] + v[583] * v[599] - v[581] * v[601] - v[580] * v[602] + v[359] * v[882] - v[360] * v[888]
		- v[347] * v[936] + v[348] * v[948];
	v[1514] = v[584] * v[597] + v[582] * v[599] - v[581] * v[600] - v[579] * v[602] + v[359] * v[881] - v[360] * v[887]
		- v[347] * v[935] + v[348] * v[947];
	v[1513] = v[2802] * v[576] - 2e0*v[582] * v[594] - v[358] * v[886] + v[360] * v[898] + v[346] * v[945] - v[348] * v[969];
	v[1531] = v[1529] * v[322] + v[1513] * v[324] + v[1438] * v[325];
	v[1530] = v[1529] * v[337] + v[1513] * v[340] + v[1438] * v[342];
	v[1512] = v[2801] * v[577] - 2e0*v[583] * v[595] - v[358] * v[885] + v[360] * v[897] + v[346] * v[942] - v[348] * v[966];
	v[1528] = v[1524] * v[322] + v[1512] * v[324] + v[1437] * v[325];
	v[1526] = v[1524] * v[337] + v[1512] * v[340] + v[1437] * v[342];
	v[1511] = -(v[583] * v[594]) - v[582] * v[595] + v[577] * v[600] + v[576] * v[601] - v[358] * v[884] + v[360] * v[896]
		+ v[346] * v[941] - v[348] * v[965];
	v[1527] = v[1523] * v[322] + v[1511] * v[324] + v[1436] * v[325];
	v[1525] = v[1523] * v[337] + v[1511] * v[340] + v[1436] * v[342];
	v[1510] = v[2800] * v[578] - 2e0*v[584] * v[596] - v[358] * v[883] + v[360] * v[895] + v[346] * v[937] - v[348] * v[961];
	v[1522] = v[1516] * v[322] + v[1510] * v[324] + v[1435] * v[325];
	v[1519] = v[1516] * v[337] + v[1510] * v[340] + v[1435] * v[342];
	v[1509] = -(v[584] * v[595]) - v[583] * v[596] + v[578] * v[601] + v[577] * v[602] - v[358] * v[882] + v[360] * v[894]
		+ v[346] * v[936] - v[348] * v[960];
	v[1521] = v[1515] * v[322] + v[1509] * v[324] + v[1434] * v[325];
	v[1518] = v[1515] * v[337] + v[1509] * v[340] + v[1434] * v[342];
	v[1508] = -(v[584] * v[594]) - v[582] * v[596] + v[578] * v[600] + v[576] * v[602] - v[358] * v[881] + v[360] * v[893]
		+ v[346] * v[935] - v[348] * v[959];
	v[1520] = v[1514] * v[322] + v[1508] * v[324] + v[1433] * v[325];
	v[1517] = v[1514] * v[337] + v[1508] * v[340] + v[1433] * v[342];
	v[1492] = v[1405] * v[358] + v[1399] * v[359] + v[1393] * v[360];
	v[1504] = v[1405] - v[1492] * v[358];
	v[1500] = v[1399] - v[1492] * v[359];
	v[1496] = v[1393] - v[1492] * v[360];
	v[1473] = v[1404] * v[358] + v[1398] * v[359] + v[1392] * v[360];
	v[1488] = v[1404] - v[1473] * v[358];
	v[1483] = v[1398] - v[1473] * v[359];
	v[1478] = v[1392] - v[1473] * v[360];
	v[1472] = v[1403] * v[358] + v[1397] * v[359] + v[1391] * v[360];
	v[1487] = v[1403] - v[1472] * v[358];
	v[1482] = v[1397] - v[1472] * v[359];
	v[1477] = v[1391] - v[1472] * v[360];
	v[1450] = v[1402] * v[358] + v[1396] * v[359] + v[1390] * v[360];
	v[1468] = v[1402] - v[1450] * v[358];
	v[1462] = v[1396] - v[1450] * v[359];
	v[1456] = v[1390] - v[1450] * v[360];
	v[1449] = v[1401] * v[358] + v[1395] * v[359] + v[1389] * v[360];
	v[1467] = v[1401] - v[1449] * v[358];
	v[1461] = v[1395] - v[1449] * v[359];
	v[1455] = v[1389] - v[1449] * v[360];
	v[1448] = v[1400] * v[358] + v[1394] * v[359] + v[1388] * v[360];
	v[1466] = v[1400] - v[1448] * v[358];
	v[1460] = v[1394] - v[1448] * v[359];
	v[1454] = v[1388] - v[1448] * v[360];
	v[1447] = -(v[360] * v[596]) - v[358] * v[602];
	v[1446] = -(v[360] * v[595]) - v[358] * v[601];
	v[1445] = -(v[360] * v[594]) - v[358] * v[600];
	v[1444] = -(v[360] * v[599]) - v[359] * v[602];
	v[1443] = -(v[360] * v[598]) - v[359] * v[601];
	v[1442] = -(v[360] * v[597]) - v[359] * v[600];
	v[1441] = -(v[2800] * v[360]);
	v[1440] = -(v[2801] * v[360]);
	v[1439] = -(v[2802] * v[360]);
	v[1214] = 1e0 - (v[360] * v[360]);
	v[1198] = -(v[359] * v[360]);
	v[1249] = (v[1198] * v[1198]);
	v[1181] = -(v[358] * v[360]);
	v[1248] = (v[1181] * v[1181]);
	v[1166] = v[1157] * v[358] + v[1160] * v[359] + v[1163] * v[360];
	v[1471] = -(v[1453] * v[358]) - v[1166] * v[596];
	v[1470] = -(v[1452] * v[358]) - v[1166] * v[595];
	v[1469] = -(v[1451] * v[358]) - v[1166] * v[594];
	v[1465] = -(v[1453] * v[359]) - v[1166] * v[599];
	v[1464] = -(v[1452] * v[359]) - v[1166] * v[598];
	v[1463] = -(v[1451] * v[359]) - v[1166] * v[597];
	v[1459] = -(v[1453] * v[360]) - v[1166] * v[602];
	v[1458] = -(v[1452] * v[360]) - v[1166] * v[601];
	v[1457] = -(v[1451] * v[360]) - v[1166] * v[600];
	v[1217] = v[1163] - v[1166] * v[360];
	v[1201] = v[1160] - v[1166] * v[359];
	v[1184] = v[1157] - v[1166] * v[358];
	v[1165] = v[1156] * v[358] + v[1159] * v[359] + v[1162] * v[360];
	v[1491] = -(v[1476] * v[358]) - v[1165] * v[596];
	v[1490] = -(v[1475] * v[358]) - v[1165] * v[595];
	v[1489] = -(v[1474] * v[358]) - v[1165] * v[594];
	v[1486] = -(v[1476] * v[359]) - v[1165] * v[599];
	v[1485] = -(v[1475] * v[359]) - v[1165] * v[598];
	v[1484] = -(v[1474] * v[359]) - v[1165] * v[597];
	v[1481] = -(v[1476] * v[360]) - v[1165] * v[602];
	v[1480] = -(v[1475] * v[360]) - v[1165] * v[601];
	v[1479] = -(v[1474] * v[360]) - v[1165] * v[600];
	v[1216] = v[1162] - v[1165] * v[360];
	v[1200] = v[1159] - v[1165] * v[359];
	v[1183] = v[1156] - v[1165] * v[358];
	v[1164] = v[1155] * v[358] + v[1158] * v[359] + v[1161] * v[360];
	v[1507] = -(v[1495] * v[358]) - v[1164] * v[596];
	v[1506] = -(v[1494] * v[358]) - v[1164] * v[595];
	v[1505] = -(v[1493] * v[358]) - v[1164] * v[594];
	v[1503] = -(v[1495] * v[359]) - v[1164] * v[599];
	v[1502] = -(v[1494] * v[359]) - v[1164] * v[598];
	v[1501] = -(v[1493] * v[359]) - v[1164] * v[597];
	v[1499] = -(v[1495] * v[360]) - v[1164] * v[602];
	v[1498] = -(v[1494] * v[360]) - v[1164] * v[601];
	v[1497] = -(v[1493] * v[360]) - v[1164] * v[600];
	v[1215] = v[1161] - v[1164] * v[360];
	v[1199] = v[1158] - v[1164] * v[359];
	v[1182] = v[1155] - v[1164] * v[358];
	v[1124] = v[360] * v[578] - v[358] * v[584] - v[348] * v[596] + v[346] * v[602];
	v[1123] = v[360] * v[577] - v[358] * v[583] - v[348] * v[595] + v[346] * v[601];
	v[1122] = v[360] * v[576] - v[358] * v[582] - v[348] * v[594] + v[346] * v[600];
	v[1121] = -(v[360] * v[581]) + v[359] * v[584] + v[348] * v[599] - v[347] * v[602];
	v[1151] = v[1121] * v[337] + v[1124] * v[340] + v[1127] * v[342];
	v[1133] = v[1121] * v[322] + v[1124] * v[324] + v[1127] * v[325];
	v[1120] = -(v[360] * v[580]) + v[359] * v[583] + v[348] * v[598] - v[347] * v[601];
	v[1150] = v[1120] * v[337] + v[1123] * v[340] + v[1126] * v[342];
	v[1132] = v[1120] * v[322] + v[1123] * v[324] + v[1126] * v[325];
	v[1119] = -(v[360] * v[579]) + v[359] * v[582] + v[348] * v[597] - v[347] * v[600];
	v[1149] = v[1119] * v[337] + v[1122] * v[340] + v[1125] * v[342];
	v[1131] = v[1119] * v[322] + v[1122] * v[324] + v[1125] * v[325];
	v[361] = -((*radius)*v[358]) + v[365];
	v[362] = -((*radius)*v[359]) + v[366];
	v[363] = -((*radius)*v[360]) + v[367];
	v[373] = -(v[309] * v[348]);
	v[374] = -(v[309] * v[357]) + v[387];
	v[375] = -(v[308] * v[348]) + v[385];
	v[376] = -(v[308] * v[357]);
	v[377] = -(v[309] * v[347]);
	v[378] = -(v[309] * v[356]) + v[397];
	v[379] = -(v[308] * v[347]) + v[395];
	v[380] = -(v[308] * v[356]);
	v[381] = -(v[309] * v[346]);
	v[382] = -(v[309] * v[355]) + v[407];
	v[383] = -(v[308] * v[346]) + v[405];
	v[384] = -(v[308] * v[355]);
	v[386] = -(v[307] * v[348]) - v[385];
	v[388] = -(v[307] * v[357]) - v[387];
	v[389] = dCi[2] * v[373] + dBi[2] * v[375] + dAi[2] * v[386];
	v[390] = dCi[2] * v[374] + dBi[2] * v[376] + dAi[2] * v[388];
	v[391] = dCi[1] * v[373] + dBi[1] * v[375] + dAi[1] * v[386];
	v[1823] = v[1819] * v[239] + v[391] * v[479];
	v[1822] = v[1818] * v[239] + v[391] * v[478];
	v[444] = v[239] * v[391];
	v[392] = dCi[1] * v[374] + dBi[1] * v[376] + dAi[1] * v[388];
	v[1830] = v[1826] * v[239] + v[392] * v[479];
	v[1829] = v[1825] * v[239] + v[392] * v[478];
	v[457] = v[239] * v[392];
	v[393] = dCi[0] * v[373] + dBi[0] * v[375] + dAi[0] * v[386];
	v[1837] = v[1833] * v[239] + v[393] * v[479];
	v[1836] = v[1832] * v[239] + v[393] * v[478];
	v[447] = v[239] * v[393];
	v[394] = dCi[0] * v[374] + dBi[0] * v[376] + dAi[0] * v[388];
	v[1844] = v[1840] * v[239] + v[394] * v[479];
	v[1843] = v[1839] * v[239] + v[394] * v[478];
	v[460] = v[239] * v[394];
	v[396] = -(v[307] * v[347]) - v[395];
	v[398] = -(v[307] * v[356]) - v[397];
	v[399] = dCi[2] * v[377] + dBi[2] * v[379] + dAi[2] * v[396];
	v[2815] = v[391] - v[399];
	v[2804] = v[391] + v[399];
	v[1857] = v[1853] * v[239] + v[399] * v[479];
	v[1945] = v[1823] + v[1857];
	v[1856] = v[1852] * v[239] + v[399] * v[478];
	v[445] = v[239] * v[399];
	v[400] = dCi[2] * v[378] + dBi[2] * v[380] + dAi[2] * v[398];
	v[2811] = v[392] - v[400];
	v[2806] = v[392] + v[400];
	v[1864] = v[1860] * v[239] + v[400] * v[479];
	v[1948] = v[1830] + v[1864];
	v[1863] = v[1859] * v[239] + v[400] * v[478];
	v[458] = v[239] * v[400];
	v[401] = dCi[1] * v[377] + dBi[1] * v[379] + dAi[1] * v[396];
	v[2816] = -v[389] - v[401];
	v[402] = dCi[1] * v[378] + dBi[1] * v[380] + dAi[1] * v[398];
	v[2812] = -v[390] - v[402];
	v[403] = dCi[0] * v[377] + dBi[0] * v[379] + dAi[0] * v[396];
	v[1876] = v[1873] * v[239] + v[403] * v[479];
	v[452] = v[239] * v[403];
	v[404] = dCi[0] * v[378] + dBi[0] * v[380] + dAi[0] * v[398];
	v[1882] = v[1879] * v[239] + v[404] * v[479];
	v[465] = v[239] * v[404];
	v[406] = -(v[307] * v[346]) - v[405];
	v[408] = -(v[307] * v[355]) - v[407];
	v[409] = dCi[2] * v[381] + dBi[2] * v[383] + dAi[2] * v[406];
	v[2814] = v[393] + v[409];
	v[1894] = v[1891] * v[239] + v[409] * v[479];
	v[1951] = v[1837] + v[1894];
	v[1893] = v[1890] * v[239] + v[409] * v[478];
	v[448] = v[239] * v[409];
	v[410] = dCi[2] * v[382] + dBi[2] * v[384] + dAi[2] * v[408];
	v[2810] = v[394] + v[410];
	v[1900] = v[1897] * v[239] + v[410] * v[479];
	v[1954] = v[1844] + v[1900];
	v[1899] = v[1896] * v[239] + v[410] * v[478];
	v[461] = v[239] * v[410];
	v[411] = dCi[1] * v[381] + dBi[1] * v[383] + dAi[1] * v[406];
	v[2803] = v[403] + v[411];
	v[1906] = v[1903] * v[239] + v[411] * v[479];
	v[1957] = v[1876] + v[1906];
	v[1956] = (v[1872] + v[1902])*v[239] + v[2803] * v[478];
	v[453] = v[239] * v[411];
	v[412] = dCi[1] * v[382] + dBi[1] * v[384] + dAi[1] * v[408];
	v[2805] = v[404] + v[412];
	v[1912] = v[1909] * v[239] + v[412] * v[479];
	v[1960] = v[1882] + v[1912];
	v[1959] = (v[1878] + v[1908])*v[239] + v[2805] * v[478];
	v[466] = v[239] * v[412];
	v[413] = dCi[0] * v[381] + dBi[0] * v[383] + dAi[0] * v[406];
	v[414] = dCi[0] * v[382] + dBi[0] * v[384] + dAi[0] * v[408];
	v[1942] = v[1748] * v[389] + v[1813] * v[742];
	v[1938] = v[1748] * v[401] + v[1867] * v[742];
	v[1930] = v[1748] * v[390] + v[1816] * v[742];
	v[1926] = v[1748] * v[402] + v[1870] * v[742];
	v[468] = v[414] * v[742];
	v[467] = v[402] * v[742];
	v[455] = v[413] * v[742];
	v[454] = v[401] * v[742];
	v[421] = v[444] + v[445];
	v[422] = v[457] + v[458];
	v[427] = v[447] + v[448];
	v[428] = v[460] + v[461];
	v[430] = v[452] + v[453];
	v[431] = v[465] + v[466];
	v[434] = v[2786] * v[389] + v[2787] * v[401] + v[2788] * v[413] + v[391] * v[418] + v[393] * v[419] + v[399] * v[420]
		+ v[403] * v[425] + v[409] * v[426] + v[411] * v[429];
	v[435] = v[2786] * v[390] + v[2787] * v[402] + v[2788] * v[414] + v[392] * v[418] + v[394] * v[419] + v[400] * v[420]
		+ v[404] * v[425] + v[410] * v[426] + v[412] * v[429];
	v[437] = -4e0*v[476];
	v[1987] = v[1975] * v[434] + v[437] * (v[1813] * v[2786] + v[1867] * v[2787] + v[1915] * v[2788] + v[403] - v[411] + d[11] *
		(-v[401] - v[413]) + v[1819] * v[418] + v[1833] * v[419] + v[1853] * v[420] + v[1873] * v[425] + v[1891] * v[426]
		+ v[1903] * v[429] + v[2804] * v[471] + v[2814] * v[472]);
	v[2808] = v[1987] + v[1748] * v[413] + v[1915] * v[742];
	v[2813] = -(v[1735] * v[389]) + v[1974] * v[434] + v[437] * (v[1812] * v[2786] + v[1866] * v[2787] + v[1914] * v[2788]
		- v[393] + v[409] + d[10] * (-v[389] - v[413]) + v[1818] * v[418] + v[1832] * v[419] + v[1852] * v[420] + v[1872] * v[425]
		+ v[1890] * v[426] + v[1902] * v[429] + v[2803] * v[472] + v[2804] * v[474]) + v[1812] * v[742];
	v[1981] = v[1975] * v[435] + v[437] * (v[1816] * v[2786] + v[1870] * v[2787] + v[1918] * v[2788] + v[404] - v[412] + d[11] *
		(-v[402] - v[414]) + v[1826] * v[418] + v[1840] * v[419] + v[1860] * v[420] + v[1879] * v[425] + v[1897] * v[426]
		+ v[1909] * v[429] + v[2806] * v[471] + v[2810] * v[472]);
	v[2807] = v[1981] + v[1748] * v[414] + v[1918] * v[742];
	v[2809] = -(v[1735] * v[390]) + v[1974] * v[435] + v[437] * (v[1815] * v[2786] + v[1869] * v[2787] + v[1917] * v[2788]
		- v[394] + v[410] + d[10] * (-v[390] - v[414]) + v[1825] * v[418] + v[1839] * v[419] + v[1859] * v[420] + v[1878] * v[425]
		+ v[1896] * v[426] + v[1908] * v[429] + v[2805] * v[472] + v[2806] * v[474]) + v[1815] * v[742];
	v[2004] = v[435] * v[437] + v[467] + v[468];
	v[2001] = v[2004] - v[467] + v[390] * v[742];
	v[1997] = v[2001] + v[467] - v[468];
	v[1995] = v[434] * v[437] + v[454] + v[455];
	v[1992] = v[1995] - v[454] + v[389] * v[742];
	v[1988] = v[1992] + v[454] - v[455];
	v[2005] = v[1882] - v[1912] + 2e0*v[2004] + v[1948] * v[471] + v[1954] * v[472] + (v[1926] + v[2807])*v[475];
	v[1996] = v[1876] - v[1906] + 2e0*v[1995] + v[1945] * v[471] + v[1951] * v[472] + (v[1938] + v[2808])*v[475];
	v[2003] = -v[1844] + v[1900] + 0.5e0*v[422] + v[1960] * v[472] + (v[1930] + v[2807])*v[473] + v[1948] * v[474];
	v[2002] = -v[1843] + v[1899] + 2e0*v[2001] + v[1959] * v[472] + (v[1829] + v[1863])*v[474] + v[473] * (v[2809]
		- v[1735] * v[414] + v[1917] * v[742]);
	v[1994] = -v[1837] + v[1894] + 0.5e0*v[421] + v[1957] * v[472] + (v[1942] + v[2808])*v[473] + v[1945] * v[474];
	v[1993] = -v[1836] + v[1893] + 2e0*v[1992] + v[1956] * v[472] + (v[1822] + v[1856])*v[474] + v[473] * (v[2813]
		- v[1735] * v[413] + v[1914] * v[742]);
	v[2000] = v[1830] - v[1864] + 0.5e0*v[428] + (v[1926] + v[1930] + v[1981])*v[470] + v[1960] * v[471]
		+ v[1954] * v[474];
	v[1999] = v[1829] - v[1863] + 0.5e0*v[431] + v[1959] * v[471] + (v[1843] + v[1899])*v[474] + v[470] * (v[2809]
		- v[1735] * v[402] + v[1869] * v[742]);
	v[1998] = 2e0*v[1997] + (v[1824] - v[1858])*v[239] + v[2811] * v[477] + v[471] * ((v[1877] + v[1907])*v[239]
		+ v[2805] * v[477]) + v[474] * ((v[1838] + v[1895])*v[239] + v[2810] * v[477]) + v[470] * (v[1741] * v[2812]
			+ v[1973] * v[435] + v[437] * (v[1814] * v[2786] + v[1868] * v[2787] + (dCi[0] * v[1796] + dBi[0] * v[1802]
				+ dAi[0] * v[1886])*v[2788] + v[2811] + d[9] * v[2812] + v[1824] * v[418] + v[1838] * v[419] + v[1858] * v[420]
				+ v[1877] * v[425] + v[1895] * v[426] + v[1907] * v[429] + v[2805] * v[471] + v[2810] * v[474]) + (v[1814] + v[1868]
					)*v[742]);
	v[1991] = v[1823] - v[1857] + 0.5e0*v[427] + (v[1938] + v[1942] + v[1987])*v[470] + v[1957] * v[471]
		+ v[1951] * v[474];
	v[1990] = v[1822] - v[1856] + 0.5e0*v[430] + v[1956] * v[471] + (v[1836] + v[1893])*v[474] + v[470] * (v[2813]
		- v[1735] * v[401] + v[1866] * v[742]);
	v[1989] = 2e0*v[1988] + (v[1817] - v[1851])*v[239] + v[2815] * v[477] + v[471] * ((v[1871] + v[1901])*v[239]
		+ v[2803] * v[477]) + v[474] * ((v[1831] + v[1889])*v[239] + v[2814] * v[477]) + v[470] * (v[1741] * v[2816]
			+ v[1973] * v[434] + v[437] * (v[1811] * v[2786] + v[1865] * v[2787] + (dCi[0] * v[1793] + dBi[0] * v[1799]
				+ dAi[0] * v[1883])*v[2788] + v[2815] + d[9] * v[2816] + v[1817] * v[418] + v[1831] * v[419] + v[1851] * v[420]
				+ v[1871] * v[425] + v[1889] * v[426] + v[1901] * v[429] + v[2803] * v[471] + v[2814] * v[474]) + (v[1811] + v[1865]
					)*v[742]);
	v[446] = v[444] - v[445] + v[1988] * v[470] + v[430] * v[471] + v[427] * v[474];
	v[451] = -v[447] + v[448] + v[430] * v[472] + v[1992] * v[473] + v[421] * v[474];
	v[456] = v[452] - v[453] + v[421] * v[471] + v[427] * v[472] + v[1995] * v[475];
	v[459] = v[457] - v[458] + v[1997] * v[470] + v[431] * v[471] + v[428] * v[474];
	v[464] = -v[460] + v[461] + v[431] * v[472] + v[2001] * v[473] + v[422] * v[474];
	v[469] = v[465] - v[466] + v[422] * v[471] + v[428] * v[472] + v[2004] * v[475];
	v[612] = 1e0 / (-(v[370] * v[371]) + v[369] * v[372]);
	v[2116] = (-(v[2005] * v[369]) + v[1996] * v[371])*v[612];
	v[2082] = (-(v[2003] * v[369]) + v[1994] * v[371])*v[612];
	v[2081] = (-(v[2002] * v[369]) + v[1993] * v[371])*v[612];
	v[2050] = (-(v[2000] * v[369]) + v[1991] * v[371])*v[612];
	v[2049] = (-(v[1999] * v[369]) + v[1990] * v[371])*v[612];
	v[2048] = (-(v[1998] * v[369]) + v[1989] * v[371])*v[612];
	v[2044] = (v[371] * v[584] - v[369] * v[593])*v[612];
	v[2043] = (v[371] * v[583] - v[369] * v[592])*v[612];
	v[2042] = (v[371] * v[582] - v[369] * v[591])*v[612];
	v[2035] = (v[371] * v[581] - v[369] * v[590])*v[612];
	v[2034] = (v[371] * v[580] - v[369] * v[589])*v[612];
	v[2033] = (v[371] * v[579] - v[369] * v[588])*v[612];
	v[2023] = (v[371] * v[578] - v[369] * v[587])*v[612];
	v[2022] = (v[371] * v[577] - v[369] * v[586])*v[612];
	v[2021] = (v[371] * v[576] - v[369] * v[585])*v[612];
	v[2020] = (v[2005] * v[370] - v[1996] * v[372])*v[612];
	v[2019] = (v[2003] * v[370] - v[1994] * v[372])*v[612];
	v[2150] = -(v[2019] * v[343]) - v[2082] * v[350] + v[963];
	v[2140] = -(v[2019] * v[344]) - v[2082] * v[352] + v[951];
	v[2128] = -(v[2019] * v[345]) - v[2082] * v[354] + v[939];
	v[2018] = (v[2002] * v[370] - v[1993] * v[372])*v[612];
	v[2017] = (v[2000] * v[370] - v[1991] * v[372])*v[612];
	v[2148] = -(v[2017] * v[343]) - v[2050] * v[350] + v[962];
	v[2138] = -(v[2017] * v[344]) - v[2050] * v[352] + v[950];
	v[2126] = -(v[2017] * v[345]) - v[2050] * v[354] + v[938];
	v[2016] = (v[1999] * v[370] - v[1990] * v[372])*v[612];
	v[2112] = -(v[2016] * v[343]) - v[2049] * v[350] + v[967];
	v[2103] = -(v[2016] * v[344]) - v[2049] * v[352] + v[955];
	v[2092] = -(v[2016] * v[345]) - v[2049] * v[354] + v[943];
	v[2015] = (v[1998] * v[370] - v[1989] * v[372])*v[612];
	v[2014] = (-(v[372] * v[584]) + v[370] * v[593])*v[612];
	v[2147] = -(v[2014] * v[343]) - v[2044] * v[350];
	v[2137] = -(v[2014] * v[344]) - v[2044] * v[352];
	v[2125] = -(v[2014] * v[345]) - v[2044] * v[354];
	v[2013] = (-(v[372] * v[583]) + v[370] * v[592])*v[612];
	v[2111] = -(v[2013] * v[343]) - v[2043] * v[350];
	v[2102] = -(v[2013] * v[344]) - v[2043] * v[352];
	v[2091] = -(v[2013] * v[345]) - v[2043] * v[354];
	v[2012] = (-(v[372] * v[582]) + v[370] * v[591])*v[612];
	v[2077] = -(v[2012] * v[343]) - v[2042] * v[350];
	v[2069] = -(v[2012] * v[344]) - v[2042] * v[352];
	v[2059] = -(v[2012] * v[345]) - v[2042] * v[354];
	v[2011] = (-(v[372] * v[581]) + v[370] * v[590])*v[612];
	v[2146] = -(v[2011] * v[343]) - v[2035] * v[350];
	v[2136] = -(v[2011] * v[344]) - v[2035] * v[352];
	v[2122] = -(v[2011] * v[345]) - v[2035] * v[354];
	v[2010] = (-(v[372] * v[580]) + v[370] * v[589])*v[612];
	v[2110] = -(v[2010] * v[343]) - v[2034] * v[350];
	v[2101] = -(v[2010] * v[344]) - v[2034] * v[352];
	v[2088] = -(v[2010] * v[345]) - v[2034] * v[354];
	v[2009] = (-(v[372] * v[579]) + v[370] * v[588])*v[612];
	v[2076] = -(v[2009] * v[343]) - v[2033] * v[350];
	v[2068] = -(v[2009] * v[344]) - v[2033] * v[352];
	v[2056] = -(v[2009] * v[345]) - v[2033] * v[354];
	v[2008] = (-(v[372] * v[578]) + v[370] * v[587])*v[612];
	v[2145] = -(v[2008] * v[343]) - v[2023] * v[350];
	v[2133] = -(v[2008] * v[344]) - v[2023] * v[352];
	v[2119] = -(v[2008] * v[345]) - v[2023] * v[354];
	v[2007] = (-(v[372] * v[577]) + v[370] * v[586])*v[612];
	v[2109] = -(v[2007] * v[343]) - v[2022] * v[350];
	v[2098] = -(v[2007] * v[344]) - v[2022] * v[352];
	v[2085] = -(v[2007] * v[345]) - v[2022] * v[354];
	v[2006] = (-(v[372] * v[576]) + v[370] * v[585])*v[612];
	v[2075] = -(v[2006] * v[343]) - v[2021] * v[350];
	v[2065] = -(v[2006] * v[344]) - v[2021] * v[352];
	v[2053] = -(v[2006] * v[345]) - v[2021] * v[354];
	v[613] = (v[355] * v[370] - v[346] * v[372])*v[612];
	v[614] = (v[356] * v[370] - v[347] * v[372])*v[612];
	v[615] = (v[357] * v[370] - v[348] * v[372])*v[612];
	v[616] = (-(v[372] * v[446]) + v[370] * v[459])*v[612];
	v[617] = (-(v[372] * v[451]) + v[370] * v[464])*v[612];
	v[618] = (-(v[372] * v[456]) + v[370] * v[469])*v[612];
	v[619] = (-(v[355] * v[369]) + v[346] * v[371])*v[612];
	v[2032] = v[2133] - v[554] * v[613] - v[563] * v[619];
	v[2031] = v[2098] - v[553] * v[613] - v[562] * v[619];
	v[2030] = v[2065] - v[552] * v[613] - v[561] * v[619];
	v[2029] = v[2119] - v[557] * v[613] - v[566] * v[619];
	v[2028] = v[2085] - v[556] * v[613] - v[565] * v[619];
	v[2027] = v[2053] - v[555] * v[613] - v[564] * v[619];
	v[2026] = v[2145] - v[551] * v[613] - v[560] * v[619];
	v[2025] = v[2109] - v[550] * v[613] - v[559] * v[619];
	v[2024] = v[2075] - v[549] * v[613] - v[558] * v[619];
	v[971] = 1e0 - v[343] * v[613] - v[350] * v[619];
	v[2829] = (*a4)*v[971];
	v[628] = -(v[345] * v[613]) - v[354] * v[619];
	v[2233] = (*a4)*v[628];
	v[626] = -(v[344] * v[613]) - v[352] * v[619];
	v[2226] = (*a4)*v[626];
	v[620] = (-(v[356] * v[369]) + v[347] * v[371])*v[612];
	v[2041] = v[2122] - v[557] * v[614] - v[566] * v[620];
	v[2040] = v[2088] - v[556] * v[614] - v[565] * v[620];
	v[2039] = v[2056] - v[555] * v[614] - v[564] * v[620];
	v[2038] = v[2136] - v[554] * v[614] - v[563] * v[620];
	v[2037] = v[2101] - v[553] * v[614] - v[562] * v[620];
	v[2036] = v[2068] - v[552] * v[614] - v[561] * v[620];
	v[981] = 1e0 - v[344] * v[614] - v[352] * v[620];
	v[2831] = (*a4)*v[981];
	v[629] = -(v[345] * v[614]) - v[354] * v[620];
	v[2235] = (*a4)*v[629];
	v[621] = (-(v[357] * v[369]) + v[348] * v[371])*v[612];
	v[2047] = v[2125] - v[557] * v[615] - v[566] * v[621];
	v[2046] = v[2091] - v[556] * v[615] - v[565] * v[621];
	v[2045] = v[2059] - v[555] * v[615] - v[564] * v[621];
	v[990] = 1e0 - v[345] * v[615] - v[354] * v[621];
	v[2834] = (*a4)*v[990];
	v[622] = (v[371] * v[446] - v[369] * v[459])*v[612];
	v[2080] = v[2148] - v[551] * v[616] - v[560] * v[622];
	v[2079] = v[2112] - v[550] * v[616] - v[559] * v[622];
	v[2078] = -(v[2015] * v[343]) - v[2048] * v[350] - v[549] * v[616] - v[558] * v[622] + v[970];
	v[2072] = v[2138] - v[554] * v[616] - v[563] * v[622];
	v[2071] = v[2103] - v[553] * v[616] - v[562] * v[622];
	v[2070] = -(v[2015] * v[344]) - v[2048] * v[352] - v[552] * v[616] - v[561] * v[622] + v[958];
	v[2062] = v[2126] - v[557] * v[616] - v[566] * v[622];
	v[2061] = v[2092] - v[556] * v[616] - v[565] * v[622];
	v[2060] = -(v[2015] * v[345]) - v[2048] * v[354] - v[555] * v[616] - v[564] * v[622] + v[946];
	v[992] = v[609] - v[345] * v[616] - v[354] * v[622];
	v[984] = v[606] - v[344] * v[616] - v[352] * v[622];
	v[975] = v[603] - v[343] * v[616] - v[350] * v[622];
	v[623] = (v[371] * v[451] - v[369] * v[464])*v[612];
	v[2115] = v[2150] - v[551] * v[617] - v[560] * v[623];
	v[2114] = -(v[2018] * v[343]) - v[2081] * v[350] - v[550] * v[617] - v[559] * v[623] + v[968];
	v[2113] = v[2112] - v[549] * v[617] - v[558] * v[623];
	v[2106] = v[2140] - v[554] * v[617] - v[563] * v[623];
	v[2105] = -(v[2018] * v[344]) - v[2081] * v[352] - v[553] * v[617] - v[562] * v[623] + v[956];
	v[2104] = v[2103] - v[552] * v[617] - v[561] * v[623];
	v[2095] = v[2128] - v[557] * v[617] - v[566] * v[623];
	v[2094] = -(v[2018] * v[345]) - v[2081] * v[354] - v[556] * v[617] - v[565] * v[623] + v[944];
	v[2093] = v[2092] - v[555] * v[617] - v[564] * v[623];
	v[994] = v[610] - v[345] * v[617] - v[354] * v[623];
	v[986] = v[607] - v[344] * v[617] - v[352] * v[623];
	v[977] = v[604] - v[343] * v[617] - v[350] * v[623];
	v[624] = (v[371] * v[456] - v[369] * v[469])*v[612];
	v[2152] = -(v[2020] * v[343]) - v[2116] * v[350] - v[551] * v[618] - v[560] * v[624] + v[964];
	v[2151] = v[2150] - v[550] * v[618] - v[559] * v[624];
	v[2149] = v[2148] - v[549] * v[618] - v[558] * v[624];
	v[2142] = -(v[2020] * v[344]) - v[2116] * v[352] - v[554] * v[618] - v[563] * v[624] + v[952];
	v[2141] = v[2140] - v[553] * v[618] - v[562] * v[624];
	v[2139] = v[2138] - v[552] * v[618] - v[561] * v[624];
	v[2130] = -(v[2020] * v[345]) - v[2116] * v[354] - v[557] * v[618] - v[566] * v[624] + v[940];
	v[2129] = v[2128] - v[556] * v[618] - v[565] * v[624];
	v[2127] = v[2126] - v[555] * v[618] - v[564] * v[624];
	v[996] = v[611] - v[345] * v[618] - v[354] * v[624];
	v[988] = v[608] - v[344] * v[618] - v[352] * v[624];
	v[979] = v[605] - v[343] * v[618] - v[350] * v[624];
	v[640] = v[309] * v[345];
	v[641] = v[309] * v[354] + v[646];
	v[642] = v[309] * v[344];
	v[643] = v[309] * v[352] + v[649];
	v[644] = v[309] * v[343];
	v[645] = v[309] * v[350] + v[652];
	v[647] = v[308] * v[345] + v[646];
	v[648] = v[308] * v[354];
	v[650] = v[308] * v[344] + v[649];
	v[651] = v[308] * v[352];
	v[653] = v[308] * v[343] + v[652];
	v[654] = v[308] * v[350];
	v[655] = v[307] * v[345] - v[646];
	v[656] = v[307] * v[354] - v[646];
	v[657] = dCi[2] * v[640] + dBi[2] * v[647] + dAi[2] * v[655];
	v[658] = dCi[2] * v[641] + dBi[2] * v[648] + dAi[2] * v[656];
	v[659] = dCi[1] * v[640] + dBi[1] * v[647] + dAi[1] * v[655];
	v[660] = dCi[1] * v[641] + dBi[1] * v[648] + dAi[1] * v[656];
	v[661] = dCi[0] * v[640] + dBi[0] * v[647] + dAi[0] * v[655];
	v[999] = v[239] * v[661];
	v[662] = dCi[0] * v[641] + dBi[0] * v[648] + dAi[0] * v[656];
	v[1001] = v[239] * v[662];
	v[663] = v[307] * v[344] - v[649];
	v[664] = v[307] * v[352] - v[649];
	v[665] = dCi[2] * v[642] + dBi[2] * v[650] + dAi[2] * v[663];
	v[666] = dCi[2] * v[643] + dBi[2] * v[651] + dAi[2] * v[664];
	v[667] = dCi[1] * v[642] + dBi[1] * v[650] + dAi[1] * v[663];
	v[668] = dCi[1] * v[643] + dBi[1] * v[651] + dAi[1] * v[664];
	v[669] = dCi[0] * v[642] + dBi[0] * v[650] + dAi[0] * v[663];
	v[1004] = v[239] * v[669];
	v[670] = dCi[0] * v[643] + dBi[0] * v[651] + dAi[0] * v[664];
	v[1008] = v[239] * v[670];
	v[671] = v[307] * v[343] - v[652];
	v[672] = v[307] * v[350] - v[652];
	v[673] = dCi[2] * v[644] + dBi[2] * v[653] + dAi[2] * v[671];
	v[1000] = v[239] * v[673];
	v[674] = dCi[2] * v[645] + dBi[2] * v[654] + dAi[2] * v[672];
	v[1002] = v[239] * v[674];
	v[675] = dCi[1] * v[644] + dBi[1] * v[653] + dAi[1] * v[671];
	v[1005] = v[239] * v[675];
	v[676] = dCi[1] * v[645] + dBi[1] * v[654] + dAi[1] * v[672];
	v[1009] = v[239] * v[676];
	v[677] = dCi[0] * v[644] + dBi[0] * v[653] + dAi[0] * v[671];
	v[678] = dCi[0] * v[645] + dBi[0] * v[654] + dAi[0] * v[672];
	v[679] = v[239] * (v[659] + v[665]);
	v[680] = v[239] * (v[660] + v[666]);
	v[2817] = v[437] * (v[2786] * v[657] + v[418] * v[659] + v[419] * v[661] + v[420] * v[665] + v[2787] * v[667]
		+ v[425] * v[669] + v[426] * v[673] + v[429] * v[675] + v[2788] * v[677]) + v[677] * v[742];
	v[1016] = v[1004] - v[1005] + v[471] * v[679] + v[475] * (v[2817] + v[667] * v[742]) + v[472] * (v[1000] + v[999]);
	v[1013] = v[1000] + (v[1004] + v[1005])*v[472] + v[474] * v[679] + v[473] * (v[2817] + v[657] * v[742]) - v[999];
	v[2818] = v[437] * (v[2786] * v[658] + v[418] * v[660] + v[419] * v[662] + v[420] * v[666] + v[2787] * v[668]
		+ v[425] * v[670] + v[426] * v[674] + v[429] * v[676] + v[2788] * v[678]) + v[678] * v[742];
	v[1017] = v[1008] - v[1009] + (v[1001] + v[1002])*v[472] + v[471] * v[680] + v[475] * (v[2818] + v[668] * v[742]);
	v[1014] = -v[1001] + v[1002] + (v[1008] + v[1009])*v[472] + v[474] * v[680] + v[473] * (v[2818] + v[658] * v[742]);
	v[1020] = v[325] * v[340] - v[324] * v[342];
	v[1021] = -(v[325] * v[337]) + v[322] * v[342];
	v[1022] = v[324] * v[337] - v[322] * v[340];
	v[1600] = v[1022] * v[886] + v[1021] * v[892] + v[1020] * v[898];
	v[1603] = v[1020] * v[1600] + v[1387] * v[322] + v[1411] * v[337];
	v[1602] = v[1021] * v[1600] + v[1387] * v[324] + v[1411] * v[340];
	v[1601] = v[1022] * v[1600] + v[1387] * v[325] + v[1411] * v[342];
	v[1593] = v[1022] * v[885] + v[1021] * v[891] + v[1020] * v[897];
	v[1599] = v[1020] * v[1593] + v[1386] * v[322] + v[1410] * v[337];
	v[1597] = v[1021] * v[1593] + v[1386] * v[324] + v[1410] * v[340];
	v[1595] = v[1022] * v[1593] + v[1386] * v[325] + v[1410] * v[342];
	v[1592] = v[1022] * v[884] + v[1021] * v[890] + v[1020] * v[896];
	v[1598] = v[1020] * v[1592] + v[1385] * v[322] + v[1409] * v[337];
	v[1596] = v[1021] * v[1592] + v[1385] * v[324] + v[1409] * v[340];
	v[1594] = v[1022] * v[1592] + v[1385] * v[325] + v[1409] * v[342];
	v[1582] = v[1022] * v[883] + v[1021] * v[889] + v[1020] * v[895];
	v[1591] = v[1020] * v[1582] + v[1384] * v[322] + v[1408] * v[337];
	v[1588] = v[1021] * v[1582] + v[1384] * v[324] + v[1408] * v[340];
	v[1585] = v[1022] * v[1582] + v[1384] * v[325] + v[1408] * v[342];
	v[1581] = v[1022] * v[882] + v[1021] * v[888] + v[1020] * v[894];
	v[1590] = v[1020] * v[1581] + v[1383] * v[322] + v[1407] * v[337];
	v[1587] = v[1021] * v[1581] + v[1383] * v[324] + v[1407] * v[340];
	v[1584] = v[1022] * v[1581] + v[1383] * v[325] + v[1407] * v[342];
	v[1580] = v[1022] * v[881] + v[1021] * v[887] + v[1020] * v[893];
	v[1589] = v[1020] * v[1580] + v[1382] * v[322] + v[1406] * v[337];
	v[1586] = v[1021] * v[1580] + v[1382] * v[324] + v[1406] * v[340];
	v[1583] = v[1022] * v[1580] + v[1382] * v[325] + v[1406] * v[342];
	v[1576] = v[1022] * v[1438] + v[1021] * v[1513] + v[1020] * v[1529];
	v[1579] = v[1020] * v[1576] + v[1531] * v[322] + v[1530] * v[337];
	v[1578] = v[1021] * v[1576] + v[1531] * v[324] + v[1530] * v[340];
	v[1577] = v[1022] * v[1576] + v[1531] * v[325] + v[1530] * v[342];
	v[1569] = v[1022] * v[1437] + v[1021] * v[1512] + v[1020] * v[1524];
	v[1575] = v[1020] * v[1569] + v[1528] * v[322] + v[1526] * v[337];
	v[1573] = v[1021] * v[1569] + v[1528] * v[324] + v[1526] * v[340];
	v[1571] = v[1022] * v[1569] + v[1528] * v[325] + v[1526] * v[342];
	v[1568] = v[1022] * v[1436] + v[1021] * v[1511] + v[1020] * v[1523];
	v[1574] = v[1020] * v[1568] + v[1527] * v[322] + v[1525] * v[337];
	v[1572] = v[1021] * v[1568] + v[1527] * v[324] + v[1525] * v[340];
	v[1570] = v[1022] * v[1568] + v[1527] * v[325] + v[1525] * v[342];
	v[1558] = v[1022] * v[1435] + v[1021] * v[1510] + v[1020] * v[1516];
	v[1567] = v[1020] * v[1558] + v[1522] * v[322] + v[1519] * v[337];
	v[1564] = v[1021] * v[1558] + v[1522] * v[324] + v[1519] * v[340];
	v[1561] = v[1022] * v[1558] + v[1522] * v[325] + v[1519] * v[342];
	v[1557] = v[1022] * v[1434] + v[1021] * v[1509] + v[1020] * v[1515];
	v[1566] = v[1020] * v[1557] + v[1521] * v[322] + v[1518] * v[337];
	v[1563] = v[1021] * v[1557] + v[1521] * v[324] + v[1518] * v[340];
	v[1560] = v[1022] * v[1557] + v[1521] * v[325] + v[1518] * v[342];
	v[1556] = v[1022] * v[1433] + v[1021] * v[1508] + v[1020] * v[1514];
	v[1565] = v[1020] * v[1556] + v[1520] * v[322] + v[1517] * v[337];
	v[1562] = v[1021] * v[1556] + v[1520] * v[324] + v[1517] * v[340];
	v[1559] = v[1022] * v[1556] + v[1520] * v[325] + v[1517] * v[342];
	v[1552] = v[1022] * v[945] + v[1021] * v[957] + v[1020] * v[969];
	v[1555] = v[1020] * v[1552] + v[1423] * v[322] + v[1422] * v[337];
	v[1554] = v[1021] * v[1552] + v[1423] * v[324] + v[1422] * v[340];
	v[1553] = v[1022] * v[1552] + v[1423] * v[325] + v[1422] * v[342];
	v[1545] = v[1022] * v[942] + v[1021] * v[954] + v[1020] * v[966];
	v[1551] = v[1020] * v[1545] + v[1421] * v[322] + v[1419] * v[337];
	v[1549] = v[1021] * v[1545] + v[1421] * v[324] + v[1419] * v[340];
	v[1547] = v[1022] * v[1545] + v[1421] * v[325] + v[1419] * v[342];
	v[1544] = v[1022] * v[941] + v[1021] * v[953] + v[1020] * v[965];
	v[1550] = v[1020] * v[1544] + v[1420] * v[322] + v[1418] * v[337];
	v[1548] = v[1021] * v[1544] + v[1420] * v[324] + v[1418] * v[340];
	v[1546] = v[1022] * v[1544] + v[1420] * v[325] + v[1418] * v[342];
	v[1534] = v[1022] * v[937] + v[1021] * v[949] + v[1020] * v[961];
	v[1543] = v[1020] * v[1534] + v[1417] * v[322] + v[1414] * v[337];
	v[1540] = v[1021] * v[1534] + v[1417] * v[324] + v[1414] * v[340];
	v[1537] = v[1022] * v[1534] + v[1417] * v[325] + v[1414] * v[342];
	v[1533] = v[1022] * v[936] + v[1021] * v[948] + v[1020] * v[960];
	v[1542] = v[1020] * v[1533] + v[1416] * v[322] + v[1413] * v[337];
	v[1539] = v[1021] * v[1533] + v[1416] * v[324] + v[1413] * v[340];
	v[1536] = v[1022] * v[1533] + v[1416] * v[325] + v[1413] * v[342];
	v[1532] = v[1022] * v[935] + v[1021] * v[947] + v[1020] * v[959];
	v[1541] = v[1020] * v[1532] + v[1415] * v[322] + v[1412] * v[337];
	v[1538] = v[1021] * v[1532] + v[1415] * v[324] + v[1412] * v[340];
	v[1535] = v[1022] * v[1532] + v[1415] * v[325] + v[1412] * v[342];
	v[1145] = v[1020] * v[596] + v[1021] * v[599] + v[1022] * v[602];
	v[1213] = v[1022] * v[1145] + v[1136] * v[325] + v[1154] * v[342];
	v[1196] = v[1021] * v[1145] + v[1136] * v[324] + v[1154] * v[340];
	v[1178] = v[1020] * v[1145] + v[1136] * v[322] + v[1154] * v[337];
	v[1144] = v[1020] * v[595] + v[1021] * v[598] + v[1022] * v[601];
	v[1212] = v[1022] * v[1144] + v[1135] * v[325] + v[1153] * v[342];
	v[1195] = v[1021] * v[1144] + v[1135] * v[324] + v[1153] * v[340];
	v[1177] = v[1020] * v[1144] + v[1135] * v[322] + v[1153] * v[337];
	v[1143] = v[1020] * v[594] + v[1021] * v[597] + v[1022] * v[600];
	v[1211] = v[1022] * v[1143] + v[1134] * v[325] + v[1152] * v[342];
	v[1194] = v[1021] * v[1143] + v[1134] * v[324] + v[1152] * v[340];
	v[1176] = v[1020] * v[1143] + v[1134] * v[322] + v[1152] * v[337];
	v[1142] = v[1020] * v[1121] + v[1021] * v[1124] + v[1022] * v[1127];
	v[1210] = v[1022] * v[1142] + v[1133] * v[325] + v[1151] * v[342];
	v[1193] = v[1021] * v[1142] + v[1133] * v[324] + v[1151] * v[340];
	v[1175] = v[1020] * v[1142] + v[1133] * v[322] + v[1151] * v[337];
	v[1141] = v[1020] * v[1120] + v[1021] * v[1123] + v[1022] * v[1126];
	v[1209] = v[1022] * v[1141] + v[1132] * v[325] + v[1150] * v[342];
	v[1192] = v[1021] * v[1141] + v[1132] * v[324] + v[1150] * v[340];
	v[1174] = v[1020] * v[1141] + v[1132] * v[322] + v[1150] * v[337];
	v[1140] = v[1020] * v[1119] + v[1021] * v[1122] + v[1022] * v[1125];
	v[1208] = v[1022] * v[1140] + v[1131] * v[325] + v[1149] * v[342];
	v[1191] = v[1021] * v[1140] + v[1131] * v[324] + v[1149] * v[340];
	v[1173] = v[1020] * v[1140] + v[1131] * v[322] + v[1149] * v[337];
	v[1139] = v[1020] * v[578] + v[1021] * v[581] + v[1022] * v[584];
	v[1207] = v[1022] * v[1139] + v[1130] * v[325] + v[1148] * v[342];
	v[1190] = v[1021] * v[1139] + v[1130] * v[324] + v[1148] * v[340];
	v[1172] = v[1020] * v[1139] + v[1130] * v[322] + v[1148] * v[337];
	v[1138] = v[1020] * v[577] + v[1021] * v[580] + v[1022] * v[583];
	v[1206] = v[1022] * v[1138] + v[1129] * v[325] + v[1147] * v[342];
	v[1189] = v[1021] * v[1138] + v[1129] * v[324] + v[1147] * v[340];
	v[1171] = v[1020] * v[1138] + v[1129] * v[322] + v[1147] * v[337];
	v[1137] = v[1020] * v[576] + v[1021] * v[579] + v[1022] * v[582];
	v[1205] = v[1022] * v[1137] + v[1128] * v[325] + v[1146] * v[342];
	v[1188] = v[1021] * v[1137] + v[1128] * v[324] + v[1146] * v[340];
	v[1170] = v[1020] * v[1137] + v[1128] * v[322] + v[1146] * v[337];
	v[1023] = v[348] * v[359] - v[347] * v[360];
	v[1024] = -(v[348] * v[358]) + v[346] * v[360];
	v[1025] = v[347] * v[358] - v[346] * v[359];
	v[1026] = v[322] * v[346] + v[324] * v[347] + v[325] * v[348];
	v[1027] = v[1023] * v[322] + v[1024] * v[324] + v[1025] * v[325];
	v[1028] = v[322] * v[358] + v[324] * v[359] + v[325] * v[360];
	v[1029] = v[1020] * v[346] + v[1021] * v[347] + v[1022] * v[348];
	v[1030] = v[1020] * v[1023] + v[1021] * v[1024] + v[1022] * v[1025];
	v[1031] = v[1020] * v[358] + v[1021] * v[359] + v[1022] * v[360];
	v[1032] = v[337] * v[346] + v[340] * v[347] + v[342] * v[348];
	v[1033] = v[1023] * v[337] + v[1024] * v[340] + v[1025] * v[342];
	v[1034] = v[337] * v[358] + v[340] * v[359] + v[342] * v[360];
	v[1036] = d[0] + (v[286] - v[295])*v[304] + (v[289] - v[298])*v[305] + (v[292] - v[301])*v[306] + (*radius)*(-
		(v[1035] * v[337]) - v[262] * v[340] - v[263] * v[342]);
	v[1038] = d[1] + (v[287] - v[296])*v[304] + (v[290] - v[299])*v[305] + (v[293] - v[302])*v[306] + (*radius)*(-
		(v[1094] * v[258] * v[337]) - v[1037] * v[340] - v[268] * v[342]);
	v[1040] = d[2] + (v[288] - v[297])*v[304] + (v[291] - v[300])*v[305] + (v[294] - v[303])*v[306] + (*radius)*(-
		(v[1102] * v[258] * v[337]) - v[272] * v[340] - v[1039] * v[342]);
	v[1609] = -(v[1381] * v[358]) - v[1375] * v[359] - v[1369] * v[360] - 2e0*v[1110] * v[594] - 2e0*v[1113] * v[597]
		- 2e0*v[1116] * v[600] + v[1040] * v[945] + v[1038] * v[957] + v[1036] * v[969];
	v[1608] = -(v[1380] * v[358]) - v[1374] * v[359] - v[1368] * v[360] - 2e0*v[1111] * v[595] - 2e0*v[1114] * v[598]
		- 2e0*v[1117] * v[601] + v[1040] * v[942] + v[1038] * v[954] + v[1036] * v[966];
	v[1607] = -(v[1379] * v[358]) - v[1373] * v[359] - v[1367] * v[360] - v[1111] * v[594] - v[1110] * v[595]
		- v[1114] * v[597] - v[1113] * v[598] - v[1117] * v[600] - v[1116] * v[601] + v[1040] * v[941] + v[1038] * v[953]
		+ v[1036] * v[965];
	v[1606] = -(v[1378] * v[358]) - v[1372] * v[359] - v[1366] * v[360] - 2e0*v[1112] * v[596] - 2e0*v[1115] * v[599]
		- 2e0*v[1118] * v[602] + v[1040] * v[937] + v[1038] * v[949] + v[1036] * v[961];
	v[1605] = -(v[1377] * v[358]) - v[1371] * v[359] - v[1365] * v[360] - v[1112] * v[595] - v[1111] * v[596]
		- v[1115] * v[598] - v[1114] * v[599] - v[1118] * v[601] - v[1117] * v[602] + v[1040] * v[936] + v[1038] * v[948]
		+ v[1036] * v[960];
	v[1604] = -(v[1376] * v[358]) - v[1370] * v[359] - v[1364] * v[360] - v[1112] * v[594] - v[1110] * v[596]
		- v[1115] * v[597] - v[1113] * v[599] - v[1118] * v[600] - v[1116] * v[602] + v[1040] * v[935] + v[1038] * v[947]
		+ v[1036] * v[959];
	v[1169] = -(v[1112] * v[358]) - v[1115] * v[359] - v[1118] * v[360] + v[1036] * v[596] + v[1038] * v[599]
		+ v[1040] * v[602];
	v[1168] = -(v[1111] * v[358]) - v[1114] * v[359] - v[1117] * v[360] + v[1036] * v[595] + v[1038] * v[598]
		+ v[1040] * v[601];
	v[1167] = -(v[1110] * v[358]) - v[1113] * v[359] - v[1116] * v[360] + v[1036] * v[594] + v[1038] * v[597]
		+ v[1040] * v[600];
	v[1048] = v[1036] * v[358] + v[1038] * v[359] + v[1040] * v[360];
	v[1627] = -v[1381] + gti[0] * (v[1020] * v[1579] + v[1603] * v[322] + v[1555] * v[337]) + gti[1] * (v[1021] * v[1579]
		+ v[1603] * v[324] + v[1555] * v[340]) + gti[2] * (v[1022] * v[1579] + v[1603] * v[325] + v[1555] * v[342])
		- v[1609] * v[358] - 2e0*v[1167] * v[594] - v[1048] * v[969];
	v[1626] = -v[1380] + gti[0] * (v[1020] * v[1575] + v[1599] * v[322] + v[1551] * v[337]) + gti[1] * (v[1021] * v[1575]
		+ v[1599] * v[324] + v[1551] * v[340]) + gti[2] * (v[1022] * v[1575] + v[1599] * v[325] + v[1551] * v[342])
		- v[1608] * v[358] - 2e0*v[1168] * v[595] - v[1048] * v[966];
	v[1625] = -v[1379] + gti[0] * (v[1020] * v[1574] + v[1598] * v[322] + v[1550] * v[337]) + gti[1] * (v[1021] * v[1574]
		+ v[1598] * v[324] + v[1550] * v[340]) + gti[2] * (v[1022] * v[1574] + v[1598] * v[325] + v[1550] * v[342])
		- v[1607] * v[358] - v[1168] * v[594] - v[1167] * v[595] - v[1048] * v[965];
	v[1624] = -v[1378] + gti[0] * (v[1020] * v[1567] + v[1591] * v[322] + v[1543] * v[337]) + gti[1] * (v[1021] * v[1567]
		+ v[1591] * v[324] + v[1543] * v[340]) + gti[2] * (v[1022] * v[1567] + v[1591] * v[325] + v[1543] * v[342])
		- v[1606] * v[358] - 2e0*v[1169] * v[596] - v[1048] * v[961];
	v[1623] = -v[1377] + gti[0] * (v[1020] * v[1566] + v[1590] * v[322] + v[1542] * v[337]) + gti[1] * (v[1021] * v[1566]
		+ v[1590] * v[324] + v[1542] * v[340]) + gti[2] * (v[1022] * v[1566] + v[1590] * v[325] + v[1542] * v[342])
		- v[1605] * v[358] - v[1169] * v[595] - v[1168] * v[596] - v[1048] * v[960];
	v[1622] = -v[1376] + gti[0] * (v[1020] * v[1565] + v[1589] * v[322] + v[1541] * v[337]) + gti[1] * (v[1021] * v[1565]
		+ v[1589] * v[324] + v[1541] * v[340]) + gti[2] * (v[1022] * v[1565] + v[1589] * v[325] + v[1541] * v[342])
		- v[1604] * v[358] - v[1169] * v[594] - v[1167] * v[596] - v[1048] * v[959];
	v[1621] = -v[1375] + gti[0] * (v[1020] * v[1578] + v[1602] * v[322] + v[1554] * v[337]) + gti[1] * (v[1021] * v[1578]
		+ v[1602] * v[324] + v[1554] * v[340]) + gti[2] * (v[1022] * v[1578] + v[1602] * v[325] + v[1554] * v[342])
		- v[1609] * v[359] - 2e0*v[1167] * v[597] - v[1048] * v[957];
	v[1620] = -v[1374] + gti[0] * (v[1020] * v[1573] + v[1597] * v[322] + v[1549] * v[337]) + gti[1] * (v[1021] * v[1573]
		+ v[1597] * v[324] + v[1549] * v[340]) + gti[2] * (v[1022] * v[1573] + v[1597] * v[325] + v[1549] * v[342])
		- v[1608] * v[359] - 2e0*v[1168] * v[598] - v[1048] * v[954];
	v[1619] = -v[1373] + gti[0] * (v[1020] * v[1572] + v[1596] * v[322] + v[1548] * v[337]) + gti[1] * (v[1021] * v[1572]
		+ v[1596] * v[324] + v[1548] * v[340]) + gti[2] * (v[1022] * v[1572] + v[1596] * v[325] + v[1548] * v[342])
		- v[1607] * v[359] - v[1168] * v[597] - v[1167] * v[598] - v[1048] * v[953];
	v[1618] = -v[1372] + gti[0] * (v[1020] * v[1564] + v[1588] * v[322] + v[1540] * v[337]) + gti[1] * (v[1021] * v[1564]
		+ v[1588] * v[324] + v[1540] * v[340]) + gti[2] * (v[1022] * v[1564] + v[1588] * v[325] + v[1540] * v[342])
		- v[1606] * v[359] - 2e0*v[1169] * v[599] - v[1048] * v[949];
	v[1617] = -v[1371] + gti[0] * (v[1020] * v[1563] + v[1587] * v[322] + v[1539] * v[337]) + gti[1] * (v[1021] * v[1563]
		+ v[1587] * v[324] + v[1539] * v[340]) + gti[2] * (v[1022] * v[1563] + v[1587] * v[325] + v[1539] * v[342])
		- v[1605] * v[359] - v[1169] * v[598] - v[1168] * v[599] - v[1048] * v[948];
	v[1616] = -v[1370] + gti[0] * (v[1020] * v[1562] + v[1586] * v[322] + v[1538] * v[337]) + gti[1] * (v[1021] * v[1562]
		+ v[1586] * v[324] + v[1538] * v[340]) + gti[2] * (v[1022] * v[1562] + v[1586] * v[325] + v[1538] * v[342])
		- v[1604] * v[359] - v[1169] * v[597] - v[1167] * v[599] - v[1048] * v[947];
	v[1615] = -v[1369] + gti[0] * (v[1020] * v[1577] + v[1601] * v[322] + v[1553] * v[337]) + gti[1] * (v[1021] * v[1577]
		+ v[1601] * v[324] + v[1553] * v[340]) + gti[2] * (v[1022] * v[1577] + v[1601] * v[325] + v[1553] * v[342])
		- v[1609] * v[360] - 2e0*v[1167] * v[600] - v[1048] * v[945];
	v[1614] = -v[1368] + gti[0] * (v[1020] * v[1571] + v[1595] * v[322] + v[1547] * v[337]) + gti[1] * (v[1021] * v[1571]
		+ v[1595] * v[324] + v[1547] * v[340]) + gti[2] * (v[1022] * v[1571] + v[1595] * v[325] + v[1547] * v[342])
		- v[1608] * v[360] - 2e0*v[1168] * v[601] - v[1048] * v[942];
	v[1613] = -v[1367] + gti[0] * (v[1020] * v[1570] + v[1594] * v[322] + v[1546] * v[337]) + gti[1] * (v[1021] * v[1570]
		+ v[1594] * v[324] + v[1546] * v[340]) + gti[2] * (v[1022] * v[1570] + v[1594] * v[325] + v[1546] * v[342])
		- v[1607] * v[360] - v[1168] * v[600] - v[1167] * v[601] - v[1048] * v[941];
	v[1612] = -v[1366] + gti[0] * (v[1020] * v[1561] + v[1585] * v[322] + v[1537] * v[337]) + gti[1] * (v[1021] * v[1561]
		+ v[1585] * v[324] + v[1537] * v[340]) + gti[2] * (v[1022] * v[1561] + v[1585] * v[325] + v[1537] * v[342])
		- v[1606] * v[360] - 2e0*v[1169] * v[602] - v[1048] * v[937];
	v[1611] = -v[1365] + gti[0] * (v[1020] * v[1560] + v[1584] * v[322] + v[1536] * v[337]) + gti[1] * (v[1021] * v[1560]
		+ v[1584] * v[324] + v[1536] * v[340]) + gti[2] * (v[1022] * v[1560] + v[1584] * v[325] + v[1536] * v[342])
		- v[1605] * v[360] - v[1169] * v[601] - v[1168] * v[602] - v[1048] * v[936];
	v[1610] = -v[1364] + gti[0] * (v[1020] * v[1559] + v[1583] * v[322] + v[1535] * v[337]) + gti[1] * (v[1021] * v[1559]
		+ v[1583] * v[324] + v[1535] * v[340]) + gti[2] * (v[1022] * v[1559] + v[1583] * v[325] + v[1535] * v[342])
		- v[1604] * v[360] - v[1169] * v[600] - v[1167] * v[602] - v[1048] * v[935];
	v[1220] = -v[1118] + gti[0] * (v[1020] * v[1210] + v[1207] * v[322] + v[1213] * v[337]) + gti[1] * (v[1021] * v[1210]
		+ v[1207] * v[324] + v[1213] * v[340]) + gti[2] * (v[1022] * v[1210] + v[1207] * v[325] + v[1213] * v[342])
		- v[1169] * v[360] - v[1048] * v[602];
	v[1219] = -v[1117] + gti[0] * (v[1020] * v[1209] + v[1206] * v[322] + v[1212] * v[337]) + gti[1] * (v[1021] * v[1209]
		+ v[1206] * v[324] + v[1212] * v[340]) + gti[2] * (v[1022] * v[1209] + v[1206] * v[325] + v[1212] * v[342])
		- v[1168] * v[360] - v[1048] * v[601];
	v[1218] = -v[1116] + gti[0] * (v[1020] * v[1208] + v[1205] * v[322] + v[1211] * v[337]) + gti[1] * (v[1021] * v[1208]
		+ v[1205] * v[324] + v[1211] * v[340]) + gti[2] * (v[1022] * v[1208] + v[1205] * v[325] + v[1211] * v[342])
		- v[1167] * v[360] - v[1048] * v[600];
	v[1204] = -v[1115] + gti[0] * (v[1020] * v[1193] + v[1190] * v[322] + v[1196] * v[337]) + gti[1] * (v[1021] * v[1193]
		+ v[1190] * v[324] + v[1196] * v[340]) + gti[2] * (v[1022] * v[1193] + v[1190] * v[325] + v[1196] * v[342])
		- v[1169] * v[359] - v[1048] * v[599];
	v[1203] = -v[1114] + gti[0] * (v[1020] * v[1192] + v[1189] * v[322] + v[1195] * v[337]) + gti[1] * (v[1021] * v[1192]
		+ v[1189] * v[324] + v[1195] * v[340]) + gti[2] * (v[1022] * v[1192] + v[1189] * v[325] + v[1195] * v[342])
		- v[1168] * v[359] - v[1048] * v[598];
	v[1202] = -v[1113] + gti[0] * (v[1020] * v[1191] + v[1188] * v[322] + v[1194] * v[337]) + gti[1] * (v[1021] * v[1191]
		+ v[1188] * v[324] + v[1194] * v[340]) + gti[2] * (v[1022] * v[1191] + v[1188] * v[325] + v[1194] * v[342])
		- v[1167] * v[359] - v[1048] * v[597];
	v[1187] = -v[1112] + gti[0] * (v[1020] * v[1175] + v[1172] * v[322] + v[1178] * v[337]) + gti[1] * (v[1021] * v[1175]
		+ v[1172] * v[324] + v[1178] * v[340]) + gti[2] * (v[1022] * v[1175] + v[1172] * v[325] + v[1178] * v[342])
		- v[1169] * v[358] - v[1048] * v[596];
	v[1186] = -v[1111] + gti[0] * (v[1020] * v[1174] + v[1171] * v[322] + v[1177] * v[337]) + gti[1] * (v[1021] * v[1174]
		+ v[1171] * v[324] + v[1177] * v[340]) + gti[2] * (v[1022] * v[1174] + v[1171] * v[325] + v[1177] * v[342])
		- v[1168] * v[358] - v[1048] * v[595];
	v[1185] = -v[1110] + gti[0] * (v[1020] * v[1173] + v[1170] * v[322] + v[1176] * v[337]) + gti[1] * (v[1021] * v[1173]
		+ v[1170] * v[324] + v[1176] * v[340]) + gti[2] * (v[1022] * v[1173] + v[1170] * v[325] + v[1176] * v[342])
		- v[1167] * v[358] - v[1048] * v[594];
	v[1041] = v[1020] * v[1029] + v[1026] * v[322] + v[1032] * v[337];
	v[1042] = v[1020] * v[1030] + v[1027] * v[322] + v[1033] * v[337];
	v[1043] = v[1020] * v[1031] + v[1028] * v[322] + v[1034] * v[337];
	v[1045] = v[1021] * v[1029] + v[1026] * v[324] + v[1032] * v[340];
	v[1046] = v[1021] * v[1030] + v[1027] * v[324] + v[1033] * v[340];
	v[1047] = v[1021] * v[1031] + v[1028] * v[324] + v[1034] * v[340];
	v[1050] = v[1022] * v[1029] + v[1026] * v[325] + v[1032] * v[342];
	v[1051] = v[1022] * v[1030] + v[1027] * v[325] + v[1033] * v[342];
	v[1052] = v[1022] * v[1031] + v[1028] * v[325] + v[1034] * v[342];
	v[1054] = (*epst)*(v[1036] + gti[0] * (v[1020] * v[1042] + v[1041] * v[322] + v[1043] * v[337]) + gti[1] *
		(v[1021] * v[1042] + v[1041] * v[324] + v[1043] * v[340]) + gti[2] * (v[1022] * v[1042] + v[1041] * v[325]
			+ v[1043] * v[342]) - v[1048] * v[358]);
	v[1055] = (*epst)*(v[1038] + gti[0] * (v[1020] * v[1046] + v[1045] * v[322] + v[1047] * v[337]) + gti[1] *
		(v[1021] * v[1046] + v[1045] * v[324] + v[1047] * v[340]) + gti[2] * (v[1022] * v[1046] + v[1045] * v[325]
			+ v[1047] * v[342]) - v[1048] * v[359]);
	v[1056] = (*epst)*(v[1040] + gti[0] * (v[1020] * v[1051] + v[1050] * v[322] + v[1052] * v[337]) + gti[1] *
		(v[1021] * v[1051] + v[1050] * v[324] + v[1052] * v[340]) + gti[2] * (v[1022] * v[1051] + v[1050] * v[325]
			+ v[1052] * v[342]) - v[1048] * v[360]);
	v[1233] = (*epst)*(v[1179] * v[1182] + v[1180] * v[1199] + v[1181] * v[1215]);
	v[1234] = (*epst)*(v[1179] * v[1183] + v[1180] * v[1200] + v[1181] * v[1216]);
	v[1235] = (*epst)*(v[1179] * v[1184] + v[1180] * v[1201] + v[1181] * v[1217]);
	v[1242] = (*epst)*(v[1180] * v[1182] + v[1197] * v[1199] + v[1198] * v[1215]);
	v[1243] = (*epst)*(v[1180] * v[1183] + v[1197] * v[1200] + v[1198] * v[1216]);
	v[1244] = (*epst)*(v[1180] * v[1184] + v[1197] * v[1201] + v[1198] * v[1217]);
	v[1251] = (*epst)*(v[1181] * v[1182] + v[1198] * v[1199] + v[1214] * v[1215]);
	v[1252] = (*epst)*(v[1181] * v[1183] + v[1198] * v[1200] + v[1214] * v[1216]);
	v[1253] = (*epst)*(v[1181] * v[1184] + v[1198] * v[1201] + v[1214] * v[1217]);
	v[2820] = -(v[2775] / v[1665]);
	v[2819] = -(v[2778] / v[1665]);
	v[2211] = v[1748] / v[1665];
	v[2196] = v[1674] * v[2211];
	v[2193] = -(v[1735] / v[1665]);
	v[2218] = v[1672] * v[2193];
	v[2177] = -(v[1741] / v[1665]);
	v[2217] = v[1668] * v[2177];
	v[2176] = -(v[2779] / v[1665]);
	v[2178] = (*a4)*v[2176] + v[1674] * v[2177] + v[2209];
	v[2174] = -(v[2776] / v[1665]);
	v[2175] = (*a4)*v[2174] + v[1672] * v[2177] + v[2192];
	v[2172] = (v[1665] - v[782]) / v[1665];
	v[2173] = (*a4)*v[2172] + v[2196] + v[2218] + v[1668] * v[472];
	v[1667] = v[1668] * v[2172] + v[1672] * v[2174] + v[1674] * v[2176];
	v[2194] = -v[2192] + v[1668] * v[2193] + (*a4)*v[2819];
	v[2198] = -(v[2777] / v[1665]);
	v[2199] = v[1674] * v[2193] + (*a4)*v[2198] + v[2213];
	v[2195] = (v[1665] - v[779]) / v[1665];
	v[2197] = (*a4)*v[2195] + v[2196] + v[2217] + v[1672] * v[471];
	v[1675] = v[1672] * v[2195] + v[1674] * v[2198] + v[1668] * v[2819];
	v[2212] = -v[2209] + v[1668] * v[2211] + (*a4)*v[2820];
	v[2216] = (v[1661] + v[1665]) / v[1665];
	v[2219] = (*a4)*v[2216] + v[2217] + v[2218] + v[1674] * v[474];
	v[2214] = -(v[2774] / v[1665]);
	v[2215] = v[1672] * v[2211] - v[2213] + (*a4)*v[2214];
	v[1681] = v[1672] * v[2214] + v[1674] * v[2216] + v[1668] * v[2820];
	v[2555] = (*a4)*v[1214] + v[1439] * v[1667] + v[1440] * v[1675] + v[1441] * v[1681];
	v[2548] = (*a4)*v[1198] + v[1442] * v[1667] + v[1443] * v[1675] + v[1444] * v[1681];
	v[2584] = (*ct)*v[1198] * v[2548];
	v[2547] = (*a4)*v[1197] + v[1427] * v[1667] + v[1428] * v[1675] + v[1429] * v[1681];
	v[2540] = (*a4)*v[1181] + v[1445] * v[1667] + v[1446] * v[1675] + v[1447] * v[1681];
	v[2583] = (*ct)*v[1181] * v[2540];
	v[2539] = (*a4)*v[1180] + v[1430] * v[1667] + v[1431] * v[1675] + v[1432] * v[1681];
	v[2572] = (*ct)*v[1180] * v[2539];
	v[2538] = (*a4)*v[1179] + v[1424] * v[1667] + v[1425] * v[1675] + v[1426] * v[1681];
	v[2401] = -(v[1667] * v[564]) - v[1675] * v[565] - v[1681] * v[566];
	v[2400] = -(v[1667] * v[555]) - v[1675] * v[556] - v[1681] * v[557];
	v[2399] = -(v[1667] * v[561]) - v[1675] * v[562] - v[1681] * v[563];
	v[2398] = -(v[1667] * v[552]) - v[1675] * v[553] - v[1681] * v[554];
	v[2397] = -(v[1667] * v[558]) - v[1675] * v[559] - v[1681] * v[560];
	v[2407] = (*cn)*(v[2397] * v[628] + v[2399] * v[629] + v[2401] * v[990]);
	v[2405] = (*cn)*(v[2397] * v[626] + v[2401] * v[629] + v[2399] * v[981]);
	v[2403] = (*cn)*(v[2399] * v[626] + v[2401] * v[628] + v[2397] * v[971]);
	v[2396] = -(v[1667] * v[549]) - v[1675] * v[550] - v[1681] * v[551];
	v[2406] = (*cn)*(v[2396] * v[628] + v[2398] * v[629] + v[2400] * v[990]);
	v[2404] = (*cn)*(v[2396] * v[626] + v[2400] * v[629] + v[2398] * v[981]);
	v[2402] = (*cn)*(v[2398] * v[626] + v[2400] * v[628] + v[2396] * v[971]);
	v[2237] = v[1667] * v[2059] + v[1675] * v[2091] + v[1681] * v[2125] + v[2834];
	v[2236] = v[1667] * v[2056] + v[1675] * v[2088] + v[1681] * v[2122] + v[2235];
	v[2234] = v[1667] * v[2053] + v[1675] * v[2085] + v[1681] * v[2119] + v[2233];
	v[2229] = v[1667] * v[2069] + v[1675] * v[2102] + v[1681] * v[2137] + v[2235];
	v[2228] = v[1667] * v[2068] + v[1675] * v[2101] + v[1681] * v[2136] + v[2831];
	v[2227] = v[1667] * v[2065] + v[1675] * v[2098] + v[1681] * v[2133] + v[2226];
	v[2222] = v[1667] * v[2077] + v[1675] * v[2111] + v[1681] * v[2147] + v[2233];
	v[2833] = (*epst)*(v[1180] * v[1181] + v[1198] * (v[1197] + v[1214])) + (*cn)*v[2235] + (*epsn)*v[629];
	v[2832] = (*epst)*(v[1180] * v[1198] + v[1181] * (v[1179] + v[1214])) + (*cn)*v[2233] + (*epsn)*v[628];
	v[2221] = v[1667] * v[2076] + v[1675] * v[2110] + v[1681] * v[2146] + v[2226];
	v[2830] = (*epst)*(v[1180] * (v[1179] + v[1197]) + v[1181] * v[1198]) + (*cn)*v[2226] + (*epsn)*v[626];
	v[2220] = v[1667] * v[2075] + v[1675] * v[2109] + v[1681] * v[2145] + v[2829];
	v[2821] = (*a4)*d[0] - (*a4)*d[6] - (*a6)*dduiP[0] + (*a6)*dduiS[0] - (*a5)*duiP[0] + (*a5)*duiS[0];
	v[2822] = (*a4)*d[1] - (*a4)*d[7] - (*a6)*dduiP[1] + (*a6)*dduiS[1] - (*a5)*duiP[1] + (*a5)*duiS[1];
	v[2823] = (*a4)*d[2] - (*a4)*d[8] - (*a6)*dduiP[2] + (*a6)*dduiS[2] - (*a5)*duiP[2] + (*a5)*duiS[2];
	v[2240] = v[1667] * v[2062] + v[1675] * v[2095] + v[1681] * v[2130] + v[2029] * v[2821] + v[2041] * v[2822]
		+ v[2047] * v[2823] + v[2178] * v[992] + v[2199] * v[994] + v[2219] * v[996];
	v[2239] = v[1667] * v[2061] + v[1675] * v[2094] + v[1681] * v[2129] + v[2028] * v[2821] + v[2040] * v[2822]
		+ v[2046] * v[2823] + v[2175] * v[992] + v[2197] * v[994] + v[2215] * v[996];
	v[2238] = v[1667] * v[2060] + v[1675] * v[2093] + v[1681] * v[2127] + v[2027] * v[2821] + v[2039] * v[2822]
		+ v[2045] * v[2823] + v[2173] * v[992] + v[2194] * v[994] + v[2212] * v[996];
	v[2232] = v[1667] * v[2072] + v[1675] * v[2106] + v[1681] * v[2142] + v[2032] * v[2821] + v[2038] * v[2822]
		+ v[2041] * v[2823] + v[2178] * v[984] + v[2199] * v[986] + v[2219] * v[988];
	v[2231] = v[1667] * v[2071] + v[1675] * v[2105] + v[1681] * v[2141] + v[2031] * v[2821] + v[2037] * v[2822]
		+ v[2040] * v[2823] + v[2175] * v[984] + v[2197] * v[986] + v[2215] * v[988];
	v[2230] = v[1667] * v[2070] + v[1675] * v[2104] + v[1681] * v[2139] + v[2030] * v[2821] + v[2036] * v[2822]
		+ v[2039] * v[2823] + v[2173] * v[984] + v[2194] * v[986] + v[2212] * v[988];
	v[2225] = v[1667] * v[2080] + v[1675] * v[2115] + v[1681] * v[2152] + v[2026] * v[2821] + v[2032] * v[2822]
		+ v[2029] * v[2823] + v[2178] * v[975] + v[2199] * v[977] + v[2219] * v[979];
	v[2224] = v[1667] * v[2079] + v[1675] * v[2114] + v[1681] * v[2151] + v[2025] * v[2821] + v[2031] * v[2822]
		+ v[2028] * v[2823] + v[2175] * v[975] + v[2197] * v[977] + v[2215] * v[979];
	v[2223] = v[1667] * v[2078] + v[1675] * v[2113] + v[1681] * v[2149] + v[2024] * v[2821] + v[2030] * v[2822]
		+ v[2027] * v[2823] + v[2173] * v[975] + v[2194] * v[977] + v[2212] * v[979];
	v[2828] = v[2790] / v[1689];
	v[2824] = -(v[1072] * v[2790]);
	v[1688] = -v[1067] - v[1306];
	v[2826] = v[1688] / v[1689];
	v[1690] = -v[1330] + v[1689];
	v[2825] = -(v[1072] * v[1690]);
	v[2529] = (v[1099] * v[1688] + v[1091] * v[1690] + v[2487] * v[263] + v[2481] * v[268] + v[2475] * v[273]
		+ v[1069] * v[2790] + v[273] * v[2824] + v[263] * v[2825] + (v[1689] * v[1689])*v[2471] * v[268] * v[2826]) / v[1689];
	v[2514] = v[1698] * v[2529];
	v[2511] = (v[1074] * v[1688] + v[1085] * v[1690] + v[2487] * v[262] + v[1688] * v[1689] * v[2471] * v[267]
		+ v[2481] * v[267] + v[2475] * v[272] + v[1107] * v[2790] + v[272] * v[2824] + v[262] * v[2825]) / v[1689];
	v[2536] = v[1696] * v[2511];
	v[2495] = -(v[2464] / v[1689]);
	v[2535] = v[1692] * v[2495];
	v[2494] = v[2827] / v[1689];
	v[2496] = (*a4)*v[2494] + v[1698] * v[2495] + v[2527];
	v[2492] = v[2793] / v[1689];
	v[2493] = (*a4)*v[2492] + v[1696] * v[2495] + v[2510];
	v[2490] = v[1690] / v[1689];
	v[2491] = v[1059] * v[1692] + (*a4)*v[2490] + v[2514] + v[2536];
	v[1691] = v[1692] * v[2490] + v[1696] * v[2492] + v[1698] * v[2494];
	v[2512] = -v[2510] + v[1692] * v[2511] + (*a4)*v[2826];
	v[2516] = v[2791] / v[1689];
	v[2517] = v[1698] * v[2511] + (*a4)*v[2516] + v[2531];
	v[2513] = (-v[1327] + v[1689]) / v[1689];
	v[2515] = v[1058] * v[1696] + (*a4)*v[2513] + v[2514] + v[2535];
	v[1699] = v[1696] * v[2513] + v[1698] * v[2516] + v[1692] * v[2826];
	v[2530] = -v[2527] + v[1692] * v[2529] + (*a4)*v[2828];
	v[2534] = (v[1685] + v[1689]) / v[1689];
	v[2537] = v[1061] * v[1698] + (*a4)*v[2534] + v[2535] + v[2536];
	v[2532] = (-v[1064] - v[1313]) / v[1689];
	v[2533] = v[1696] * v[2529] - v[2531] + (*a4)*v[2532];
	v[1705] = v[1696] * v[2532] + v[1698] * v[2534] + v[1692] * v[2828];
	v[2561] = v[1610] * v[1667] + v[1611] * v[1675] + v[1612] * v[1681] + v[1499] * v[1691] + v[1481] * v[1699]
		+ v[1459] * v[1705] + v[1218] * v[2178] + v[1219] * v[2199] + v[1220] * v[2219] + v[1447] * v[2821] + v[1444] * v[2822]
		+ v[1441] * v[2823];
	v[2560] = v[1613] * v[1667] + v[1614] * v[1675] + v[1611] * v[1681] + v[1498] * v[1691] + v[1480] * v[1699]
		+ v[1458] * v[1705] + v[1218] * v[2175] + v[1219] * v[2197] + v[1220] * v[2215] + v[1446] * v[2821] + v[1443] * v[2822]
		+ v[1440] * v[2823];
	v[2559] = v[1615] * v[1667] + v[1613] * v[1675] + v[1610] * v[1681] + v[1497] * v[1691] + v[1479] * v[1699]
		+ v[1457] * v[1705] + v[1218] * v[2173] + v[1219] * v[2194] + v[1220] * v[2212] + v[1445] * v[2821] + v[1442] * v[2822]
		+ v[1439] * v[2823];
	v[2558] = v[1457] * v[1667] + v[1458] * v[1675] + v[1459] * v[1681] + v[1454] * v[1691] + v[1455] * v[1699]
		+ v[1456] * v[1705] + v[1215] * v[2496] + v[1216] * v[2517] + v[1217] * v[2537];
	v[2557] = v[1479] * v[1667] + v[1480] * v[1675] + v[1481] * v[1681] + v[1477] * v[1691] + v[1478] * v[1699]
		+ v[1455] * v[1705] + v[1215] * v[2493] + v[1216] * v[2515] + v[1217] * v[2533];
	v[2556] = v[1497] * v[1667] + v[1498] * v[1675] + v[1499] * v[1681] + v[1496] * v[1691] + v[1477] * v[1699]
		+ v[1454] * v[1705] + v[1215] * v[2491] + v[1216] * v[2512] + v[1217] * v[2530];
	v[2554] = v[1616] * v[1667] + v[1617] * v[1675] + v[1618] * v[1681] + v[1503] * v[1691] + v[1486] * v[1699]
		+ v[1465] * v[1705] + v[1202] * v[2178] + v[1203] * v[2199] + v[1204] * v[2219] + v[1432] * v[2821] + v[1429] * v[2822]
		+ v[1444] * v[2823];
	v[2553] = v[1619] * v[1667] + v[1620] * v[1675] + v[1617] * v[1681] + v[1502] * v[1691] + v[1485] * v[1699]
		+ v[1464] * v[1705] + v[1202] * v[2175] + v[1203] * v[2197] + v[1204] * v[2215] + v[1431] * v[2821] + v[1428] * v[2822]
		+ v[1443] * v[2823];
	v[2552] = v[1621] * v[1667] + v[1619] * v[1675] + v[1616] * v[1681] + v[1501] * v[1691] + v[1484] * v[1699]
		+ v[1463] * v[1705] + v[1202] * v[2173] + v[1203] * v[2194] + v[1204] * v[2212] + v[1430] * v[2821] + v[1427] * v[2822]
		+ v[1442] * v[2823];
	v[2551] = v[1463] * v[1667] + v[1464] * v[1675] + v[1465] * v[1681] + v[1460] * v[1691] + v[1461] * v[1699]
		+ v[1462] * v[1705] + v[1199] * v[2496] + v[1200] * v[2517] + v[1201] * v[2537];
	v[2550] = v[1484] * v[1667] + v[1485] * v[1675] + v[1486] * v[1681] + v[1482] * v[1691] + v[1483] * v[1699]
		+ v[1461] * v[1705] + v[1199] * v[2493] + v[1200] * v[2515] + v[1201] * v[2533];
	v[2549] = v[1501] * v[1667] + v[1502] * v[1675] + v[1503] * v[1681] + v[1500] * v[1691] + v[1482] * v[1699]
		+ v[1460] * v[1705] + v[1199] * v[2491] + v[1200] * v[2512] + v[1201] * v[2530];
	v[2546] = v[1622] * v[1667] + v[1623] * v[1675] + v[1624] * v[1681] + v[1507] * v[1691] + v[1491] * v[1699]
		+ v[1471] * v[1705] + v[1185] * v[2178] + v[1186] * v[2199] + v[1187] * v[2219] + v[1426] * v[2821] + v[1432] * v[2822]
		+ v[1447] * v[2823];
	v[2545] = v[1625] * v[1667] + v[1626] * v[1675] + v[1623] * v[1681] + v[1506] * v[1691] + v[1490] * v[1699]
		+ v[1470] * v[1705] + v[1185] * v[2175] + v[1186] * v[2197] + v[1187] * v[2215] + v[1425] * v[2821] + v[1431] * v[2822]
		+ v[1446] * v[2823];
	v[2544] = v[1627] * v[1667] + v[1625] * v[1675] + v[1622] * v[1681] + v[1505] * v[1691] + v[1489] * v[1699]
		+ v[1469] * v[1705] + v[1185] * v[2173] + v[1186] * v[2194] + v[1187] * v[2212] + v[1424] * v[2821] + v[1430] * v[2822]
		+ v[1445] * v[2823];
	v[2543] = v[1469] * v[1667] + v[1470] * v[1675] + v[1471] * v[1681] + v[1466] * v[1691] + v[1467] * v[1699]
		+ v[1468] * v[1705] + v[1182] * v[2496] + v[1183] * v[2517] + v[1184] * v[2537];
	v[2542] = v[1489] * v[1667] + v[1490] * v[1675] + v[1491] * v[1681] + v[1487] * v[1691] + v[1488] * v[1699]
		+ v[1467] * v[1705] + v[1182] * v[2493] + v[1183] * v[2515] + v[1184] * v[2533];
	v[2541] = v[1505] * v[1667] + v[1506] * v[1675] + v[1507] * v[1681] + v[1504] * v[1691] + v[1487] * v[1699]
		+ v[1466] * v[1705] + v[1182] * v[2491] + v[1183] * v[2512] + v[1184] * v[2530];
	v[1706] = v[2822] * v[626] + v[2823] * v[628] + v[2821] * v[971] + v[1667] * v[975] + v[1675] * v[977] + v[1681] * v[979];
	v[1707] = v[2821] * v[626] + v[2823] * v[629] + v[2822] * v[981] + v[1667] * v[984] + v[1675] * v[986] + v[1681] * v[988];
	v[1708] = v[2821] * v[628] + v[2822] * v[629] + v[2823] * v[990] + v[1667] * v[992] + v[1675] * v[994] + v[1681] * v[996];
	v[2413] = (*cn)*(-(v[1706] * v[560]) - v[1707] * v[563] - v[1708] * v[566] + v[2397] * v[979] + v[2399] * v[988]
		+ v[2401] * v[996]);
	v[2412] = (*cn)*(-(v[1706] * v[551]) - v[1707] * v[554] - v[1708] * v[557] + v[2396] * v[979] + v[2398] * v[988]
		+ v[2400] * v[996]);
	v[2411] = (*cn)*(-(v[1706] * v[559]) - v[1707] * v[562] - v[1708] * v[565] + v[2397] * v[977] + v[2399] * v[986]
		+ v[2401] * v[994]);
	v[2410] = (*cn)*(-(v[1706] * v[550]) - v[1707] * v[553] - v[1708] * v[556] + v[2396] * v[977] + v[2398] * v[986]
		+ v[2400] * v[994]);
	v[2409] = (*cn)*(-(v[1706] * v[558]) - v[1707] * v[561] - v[1708] * v[564] + v[2397] * v[975] + v[2399] * v[984]
		+ v[2401] * v[992]);
	v[2408] = (*cn)*(-(v[1706] * v[549]) - v[1707] * v[552] - v[1708] * v[555] + v[2396] * v[975] + v[2398] * v[984]
		+ v[2400] * v[992]);
	v[1709] = v[1185] * v[1667] + v[1186] * v[1675] + v[1187] * v[1681] + v[1182] * v[1691] + v[1183] * v[1699]
		+ v[1184] * v[1705] + v[1179] * v[2821] + v[1180] * v[2822] + v[1181] * v[2823];
	v[1710] = v[1202] * v[1667] + v[1203] * v[1675] + v[1204] * v[1681] + v[1199] * v[1691] + v[1200] * v[1699]
		+ v[1201] * v[1705] + v[1180] * v[2821] + v[1197] * v[2822] + v[1198] * v[2823];
	v[1711] = v[1218] * v[1667] + v[1219] * v[1675] + v[1220] * v[1681] + v[1215] * v[1691] + v[1216] * v[1699]
		+ v[1217] * v[1705] + v[1181] * v[2821] + v[1198] * v[2822] + v[1214] * v[2823];
	v[2858] = (*epst)*(v[1186] * v[1187] + v[1203] * v[1204] + v[1219] * v[1220]) + v[1056] * v[1611] + v[1055] * v[1617]
		+ v[1054] * v[1623] + (*ct)*(v[1623] * v[1709] + v[1617] * v[1710] + v[1611] * v[1711]) + (*epsn)*(v[604] * v[605]
			+ v[607] * v[608] + v[610] * v[611] + v[1016] * v[617] + v[1017] * v[623] + v[363] * v[939] + v[362] * v[951]
			+ v[361] * v[963]);
	v[2857] = (*epst)*(v[1185] * v[1187] + v[1202] * v[1204] + v[1218] * v[1220]) + v[1056] * v[1610] + v[1055] * v[1616]
		+ v[1054] * v[1622] + (*ct)*(v[1622] * v[1709] + v[1616] * v[1710] + v[1610] * v[1711]) + (*epsn)*(v[603] * v[605]
			+ v[606] * v[608] + v[609] * v[611] + v[1016] * v[616] + v[1017] * v[622] + v[363] * v[938] + v[362] * v[950]
			+ v[361] * v[962]);
	v[2855] = (*epst)*(v[1184] * v[1187] + v[1201] * v[1204] + v[1217] * v[1220]) + v[1056] * v[1459] + v[1055] * v[1465]
		+ v[1054] * v[1471] + (*ct)*(v[1471] * v[1709] + v[1465] * v[1710] + v[1459] * v[1711]);
	v[2852] = (*epst)*(v[1183] * v[1187] + v[1200] * v[1204] + v[1216] * v[1220]) + v[1056] * v[1481] + v[1055] * v[1486]
		+ v[1054] * v[1491] + (*ct)*(v[1491] * v[1709] + v[1486] * v[1710] + v[1481] * v[1711]);
	v[2848] = (*epst)*(v[1182] * v[1187] + v[1199] * v[1204] + v[1215] * v[1220]) + v[1056] * v[1499] + v[1055] * v[1503]
		+ v[1054] * v[1507] + (*ct)*(v[1507] * v[1709] + v[1503] * v[1710] + v[1499] * v[1711]);
	v[2843] = (*epst)*(v[1181] * v[1187] + v[1198] * v[1204] + v[1214] * v[1220]) + v[1056] * v[1441] + v[1055] * v[1444]
		+ v[1054] * v[1447] + (*ct)*(v[1447] * v[1709] + v[1444] * v[1710] + v[1441] * v[1711]) + (*epsn)*v[996];
	v[2842] = (*epst)*(v[1180] * v[1187] + v[1197] * v[1204] + v[1198] * v[1220]) + v[1055] * v[1429] + v[1054] * v[1432]
		+ v[1056] * v[1444] + (*ct)*(v[1432] * v[1709] + v[1429] * v[1710] + v[1444] * v[1711]) + (*epsn)*v[988];
	v[2841] = (*epst)*(v[1179] * v[1187] + v[1180] * v[1204] + v[1181] * v[1220]) + v[1054] * v[1426] + v[1055] * v[1432]
		+ v[1056] * v[1447] + (*ct)*(v[1426] * v[1709] + v[1432] * v[1710] + v[1447] * v[1711]) + (*epsn)*v[979];
	v[2856] = (*epst)*(v[1185] * v[1186] + v[1202] * v[1203] + v[1218] * v[1219]) + v[1056] * v[1613] + v[1055] * v[1619]
		+ v[1054] * v[1625] + (*ct)*(v[1625] * v[1709] + v[1619] * v[1710] + v[1613] * v[1711]) + (*epsn)*(v[603] * v[604]
			+ v[606] * v[607] + v[609] * v[610] + v[1013] * v[616] + v[1014] * v[622] + v[363] * v[943] + v[362] * v[955]
			+ v[361] * v[967]);
	v[2854] = (*epst)*(v[1184] * v[1186] + v[1201] * v[1203] + v[1217] * v[1219]) + v[1056] * v[1458] + v[1055] * v[1464]
		+ v[1054] * v[1470] + (*ct)*(v[1470] * v[1709] + v[1464] * v[1710] + v[1458] * v[1711]);
	v[2851] = (*epst)*(v[1183] * v[1186] + v[1200] * v[1203] + v[1216] * v[1219]) + v[1056] * v[1480] + v[1055] * v[1485]
		+ v[1054] * v[1490] + (*ct)*(v[1490] * v[1709] + v[1485] * v[1710] + v[1480] * v[1711]);
	v[2847] = (*epst)*(v[1182] * v[1186] + v[1199] * v[1203] + v[1215] * v[1219]) + v[1056] * v[1498] + v[1055] * v[1502]
		+ v[1054] * v[1506] + (*ct)*(v[1506] * v[1709] + v[1502] * v[1710] + v[1498] * v[1711]);
	v[2840] = (*epst)*(v[1181] * v[1186] + v[1198] * v[1203] + v[1214] * v[1219]) + v[1056] * v[1440] + v[1055] * v[1443]
		+ v[1054] * v[1446] + (*ct)*(v[1446] * v[1709] + v[1443] * v[1710] + v[1440] * v[1711]) + (*epsn)*v[994];
	v[2839] = (*epst)*(v[1180] * v[1186] + v[1197] * v[1203] + v[1198] * v[1219]) + v[1055] * v[1428] + v[1054] * v[1431]
		+ v[1056] * v[1443] + (*ct)*(v[1431] * v[1709] + v[1428] * v[1710] + v[1443] * v[1711]) + (*epsn)*v[986];
	v[2838] = (*epst)*(v[1179] * v[1186] + v[1180] * v[1203] + v[1181] * v[1219]) + v[1054] * v[1425] + v[1055] * v[1431]
		+ v[1056] * v[1446] + (*ct)*(v[1425] * v[1709] + v[1431] * v[1710] + v[1446] * v[1711]) + (*epsn)*v[977];
	v[2853] = (*epst)*(v[1184] * v[1185] + v[1201] * v[1202] + v[1217] * v[1218]) + v[1056] * v[1457] + v[1055] * v[1463]
		+ v[1054] * v[1469] + (*ct)*(v[1469] * v[1709] + v[1463] * v[1710] + v[1457] * v[1711]);
	v[2850] = (*epst)*(v[1183] * v[1185] + v[1200] * v[1202] + v[1216] * v[1218]) + v[1056] * v[1479] + v[1055] * v[1484]
		+ v[1054] * v[1489] + (*ct)*(v[1489] * v[1709] + v[1484] * v[1710] + v[1479] * v[1711]);
	v[2846] = (*epst)*(v[1182] * v[1185] + v[1199] * v[1202] + v[1215] * v[1218]) + v[1056] * v[1497] + v[1055] * v[1501]
		+ v[1054] * v[1505] + (*ct)*(v[1505] * v[1709] + v[1501] * v[1710] + v[1497] * v[1711]);
	v[2837] = (*epst)*(v[1181] * v[1185] + v[1198] * v[1202] + v[1214] * v[1218]) + v[1056] * v[1439] + v[1055] * v[1442]
		+ v[1054] * v[1445] + (*ct)*(v[1445] * v[1709] + v[1442] * v[1710] + v[1439] * v[1711]) + (*epsn)*v[992];
	v[2836] = (*epst)*(v[1180] * v[1185] + v[1197] * v[1202] + v[1198] * v[1218]) + v[1055] * v[1427] + v[1054] * v[1430]
		+ v[1056] * v[1442] + (*ct)*(v[1430] * v[1709] + v[1427] * v[1710] + v[1442] * v[1711]) + (*epsn)*v[984];
	v[2835] = (*epst)*(v[1179] * v[1185] + v[1180] * v[1202] + v[1181] * v[1218]) + v[1054] * v[1424] + v[1055] * v[1430]
		+ v[1056] * v[1445] + (*ct)*(v[1424] * v[1709] + v[1430] * v[1710] + v[1445] * v[1711]) + (*epsn)*v[975];
	v[2849] = (*epst)*(v[1183] * v[1184] + v[1200] * v[1201] + v[1216] * v[1217]) + v[1056] * v[1455] + v[1055] * v[1461]
		+ v[1054] * v[1467] + (*ct)*(v[1467] * v[1709] + v[1461] * v[1710] + v[1455] * v[1711]);
	v[2845] = (*epst)*(v[1182] * v[1184] + v[1199] * v[1201] + v[1215] * v[1217]) + v[1056] * v[1454] + v[1055] * v[1460]
		+ v[1054] * v[1466] + (*ct)*(v[1466] * v[1709] + v[1460] * v[1710] + v[1454] * v[1711]);
	v[2844] = (*epst)*(v[1182] * v[1183] + v[1199] * v[1200] + v[1215] * v[1216]) + v[1056] * v[1477] + v[1055] * v[1482]
		+ v[1054] * v[1487] + (*ct)*(v[1487] * v[1709] + v[1482] * v[1710] + v[1477] * v[1711]);
	v[2670] = v[1054] * v[1179] + v[1055] * v[1180] + v[1056] * v[1181] + (*ct)*(v[1179] * v[1709] + v[1180] * v[1710]
		+ v[1181] * v[1711]) + (*cn)*(v[1707] * v[626] + v[1708] * v[628] + v[1706] * v[971]) + (*epsn)*(v[362] * v[626]
			+ v[363] * v[628] + v[361] * v[971]);
	v[2671] = v[1054] * v[1180] + v[1055] * v[1197] + v[1056] * v[1198] + (*ct)*(v[1180] * v[1709] + v[1197] * v[1710]
		+ v[1198] * v[1711]) + (*cn)*(v[1706] * v[626] + v[1708] * v[629] + v[1707] * v[981]) + (*epsn)*(v[361] * v[626]
			+ v[363] * v[629] + v[362] * v[981]);
	v[2672] = v[1054] * v[1181] + v[1055] * v[1198] + v[1056] * v[1214] + (*ct)*(v[1181] * v[1709] + v[1198] * v[1710]
		+ v[1214] * v[1711]) + (*cn)*(v[1706] * v[628] + v[1707] * v[629] + v[1708] * v[990]) + (*epsn)*(v[361] * v[628]
			+ v[362] * v[629] + v[363] * v[990]);
	v[2679] = (*epst)*((v[1179] * v[1179]) + v[1239] + v[1248]) + (*ct)*v[1179] * v[2538] + v[2572] + v[2583] + (*cn
		)*v[2829] + v[2402] * v[613] + v[2403] * v[619] + (*epsn)*v[971];
	v[2680] = (*ct)*(v[1179] * v[2539] + v[1180] * v[2547] + v[1181] * v[2548]) + v[2830] + v[2402] * v[614]
		+ v[2403] * v[620];
	v[2681] = (*ct)*(v[1179] * v[2540] + v[1180] * v[2548] + v[1181] * v[2555]) + v[2832] + v[2402] * v[615]
		+ v[2403] * v[621];
	v[2682] = v[1233] + (*ct)*(v[1179] * v[2541] + v[1180] * v[2549] + v[1181] * v[2556]);
	v[2683] = v[1234] + (*ct)*(v[1179] * v[2542] + v[1180] * v[2550] + v[1181] * v[2557]);
	v[2684] = v[1235] + (*ct)*(v[1179] * v[2543] + v[1180] * v[2551] + v[1181] * v[2558]);
	v[2685] = (*ct)*(v[1179] * v[2544] + v[1180] * v[2552] + v[1181] * v[2559]) + v[2835] + v[2402] * v[616]
		+ v[2403] * v[622] + (*cn)*(v[1706] * v[2024] + v[1708] * v[2027] + v[1707] * v[2030] + v[2230] * v[626]
			+ v[2238] * v[628] + v[2223] * v[971]);
	v[2686] = (*ct)*(v[1179] * v[2545] + v[1180] * v[2553] + v[1181] * v[2560]) + v[2838] + v[2402] * v[617]
		+ v[2403] * v[623] + (*cn)*(v[1706] * v[2025] + v[1708] * v[2028] + v[1707] * v[2031] + v[2231] * v[626]
			+ v[2239] * v[628] + v[2224] * v[971]);
	v[2687] = (*ct)*(v[1179] * v[2546] + v[1180] * v[2554] + v[1181] * v[2561]) + v[2841] + v[2402] * v[618]
		+ v[2403] * v[624] + (*cn)*(v[1706] * v[2026] + v[1708] * v[2029] + v[1707] * v[2032] + v[2232] * v[626]
			+ v[2240] * v[628] + v[2225] * v[971]);
	v[2688] = (*ct)*(v[1180] * v[2538] + v[1197] * v[2539] + v[1198] * v[2540]) + v[2830] + v[2404] * v[613]
		+ v[2405] * v[619];
	v[2689] = (*epst)*((v[1197] * v[1197]) + v[1239] + v[1249]) + (*ct)*v[1197] * v[2547] + v[2572] + v[2584] + (*cn
		)*v[2831] + v[2404] * v[614] + v[2405] * v[620] + (*epsn)*v[981];
	v[2690] = (*ct)*(v[1180] * v[2540] + v[1197] * v[2548] + v[1198] * v[2555]) + v[2833] + v[2404] * v[615]
		+ v[2405] * v[621];
	v[2691] = v[1242] + (*ct)*(v[1180] * v[2541] + v[1197] * v[2549] + v[1198] * v[2556]);
	v[2692] = v[1243] + (*ct)*(v[1180] * v[2542] + v[1197] * v[2550] + v[1198] * v[2557]);
	v[2693] = v[1244] + (*ct)*(v[1180] * v[2543] + v[1197] * v[2551] + v[1198] * v[2558]);
	v[2694] = (*ct)*(v[1180] * v[2544] + v[1197] * v[2552] + v[1198] * v[2559]) + v[2836] + v[2404] * v[616]
		+ v[2405] * v[622] + (*cn)*(v[1706] * v[2030] + v[1707] * v[2036] + v[1708] * v[2039] + v[2223] * v[626]
			+ v[2238] * v[629] + v[2230] * v[981]);
	v[2695] = (*ct)*(v[1180] * v[2545] + v[1197] * v[2553] + v[1198] * v[2560]) + v[2839] + v[2404] * v[617]
		+ v[2405] * v[623] + (*cn)*(v[1706] * v[2031] + v[1707] * v[2037] + v[1708] * v[2040] + v[2224] * v[626]
			+ v[2239] * v[629] + v[2231] * v[981]);
	v[2696] = (*ct)*(v[1180] * v[2546] + v[1197] * v[2554] + v[1198] * v[2561]) + v[2842] + v[2404] * v[618]
		+ v[2405] * v[624] + (*cn)*(v[1706] * v[2032] + v[1707] * v[2038] + v[1708] * v[2041] + v[2225] * v[626]
			+ v[2240] * v[629] + v[2232] * v[981]);
	v[2697] = (*ct)*(v[1181] * v[2538] + v[1198] * v[2539] + v[1214] * v[2540]) + v[2832] + v[2406] * v[613]
		+ v[2407] * v[619];
	v[2698] = (*ct)*(v[1181] * v[2539] + v[1198] * v[2547] + v[1214] * v[2548]) + v[2833] + v[2406] * v[614]
		+ v[2407] * v[620];
	v[2699] = (*epst)*((v[1214] * v[1214]) + v[1248] + v[1249]) + (*ct)*v[1214] * v[2555] + v[2583] + v[2584] + (*cn
		)*v[2834] + v[2406] * v[615] + v[2407] * v[621] + (*epsn)*v[990];
	v[2700] = v[1251] + (*ct)*(v[1181] * v[2541] + v[1198] * v[2549] + v[1214] * v[2556]);
	v[2701] = v[1252] + (*ct)*(v[1181] * v[2542] + v[1198] * v[2550] + v[1214] * v[2557]);
	v[2702] = v[1253] + (*ct)*(v[1181] * v[2543] + v[1198] * v[2551] + v[1214] * v[2558]);
	v[2703] = (*ct)*(v[1181] * v[2544] + v[1198] * v[2552] + v[1214] * v[2559]) + v[2837] + v[2406] * v[616]
		+ v[2407] * v[622] + (*cn)*(v[1706] * v[2027] + v[1707] * v[2039] + v[1708] * v[2045] + v[2223] * v[628]
			+ v[2230] * v[629] + v[2238] * v[990]);
	v[2704] = (*ct)*(v[1181] * v[2545] + v[1198] * v[2553] + v[1214] * v[2560]) + v[2840] + v[2406] * v[617]
		+ v[2407] * v[623] + (*cn)*(v[1706] * v[2028] + v[1707] * v[2040] + v[1708] * v[2046] + v[2224] * v[628]
			+ v[2231] * v[629] + v[2239] * v[990]);
	v[2705] = (*ct)*(v[1181] * v[2546] + v[1198] * v[2554] + v[1214] * v[2561]) + v[2843] + v[2406] * v[618]
		+ v[2407] * v[624] + (*cn)*(v[1706] * v[2029] + v[1707] * v[2041] + v[1708] * v[2047] + v[2225] * v[628]
			+ v[2232] * v[629] + v[2240] * v[990]);
	v[2706] = v[1233] + (*ct)*(v[1182] * v[2538] + v[1199] * v[2539] + v[1215] * v[2540]);
	v[2707] = v[1242] + (*ct)*(v[1182] * v[2539] + v[1199] * v[2547] + v[1215] * v[2548]);
	v[2708] = v[1251] + (*ct)*(v[1182] * v[2540] + v[1199] * v[2548] + v[1215] * v[2555]);
	v[2715] = v[1234] + (*ct)*(v[1183] * v[2538] + v[1200] * v[2539] + v[1216] * v[2540]);
	v[2716] = v[1243] + (*ct)*(v[1183] * v[2539] + v[1200] * v[2547] + v[1216] * v[2548]);
	v[2717] = v[1252] + (*ct)*(v[1183] * v[2540] + v[1200] * v[2548] + v[1216] * v[2555]);
	v[2724] = v[1235] + (*ct)*(v[1184] * v[2538] + v[1201] * v[2539] + v[1217] * v[2540]);
	v[2725] = v[1244] + (*ct)*(v[1184] * v[2539] + v[1201] * v[2547] + v[1217] * v[2548]);
	v[2726] = v[1253] + (*ct)*(v[1184] * v[2540] + v[1201] * v[2548] + v[1217] * v[2555]);
	v[2734] = (*ct)*(v[1185] * v[2538] + v[1202] * v[2539] + v[1218] * v[2540]) + v[2835] + v[2408] * v[613]
		+ v[2409] * v[619] + (*cn)*(v[1708] * v[2053] + v[1707] * v[2065] + v[1706] * v[2075] + v[2220] * v[975]
			+ v[2227] * v[984] + v[2234] * v[992]);
	v[2736] = (*ct)*(v[1185] * v[2539] + v[1202] * v[2547] + v[1218] * v[2548]) + v[2836] + v[2408] * v[614]
		+ v[2409] * v[620] + (*cn)*(v[1708] * v[2056] + v[1707] * v[2068] + v[1706] * v[2076] + v[2221] * v[975]
			+ v[2228] * v[984] + v[2236] * v[992]);
	v[2738] = (*ct)*(v[1185] * v[2540] + v[1202] * v[2548] + v[1218] * v[2555]) + v[2837] + v[2408] * v[615]
		+ v[2409] * v[621] + (*cn)*(v[1708] * v[2059] + v[1707] * v[2069] + v[1706] * v[2077] + v[2222] * v[975]
			+ v[2229] * v[984] + v[2237] * v[992]);
	v[2746] = (*ct)*(v[1186] * v[2538] + v[1203] * v[2539] + v[1219] * v[2540]) + v[2838] + v[2410] * v[613]
		+ v[2411] * v[619] + (*cn)*(v[1708] * v[2085] + v[1707] * v[2098] + v[1706] * v[2109] + v[2220] * v[977]
			+ v[2227] * v[986] + v[2234] * v[994]);
	v[2748] = (*ct)*(v[1186] * v[2539] + v[1203] * v[2547] + v[1219] * v[2548]) + v[2839] + v[2410] * v[614]
		+ v[2411] * v[620] + (*cn)*(v[1708] * v[2088] + v[1707] * v[2101] + v[1706] * v[2110] + v[2221] * v[977]
			+ v[2228] * v[986] + v[2236] * v[994]);
	v[2750] = (*ct)*(v[1186] * v[2540] + v[1203] * v[2548] + v[1219] * v[2555]) + v[2840] + v[2410] * v[615]
		+ v[2411] * v[621] + (*cn)*(v[1708] * v[2091] + v[1707] * v[2102] + v[1706] * v[2111] + v[2222] * v[977]
			+ v[2229] * v[986] + v[2237] * v[994]);
	v[2759] = (*ct)*(v[1187] * v[2538] + v[1204] * v[2539] + v[1220] * v[2540]) + v[2841] + v[2412] * v[613]
		+ v[2413] * v[619] + (*cn)*(v[1708] * v[2119] + v[1707] * v[2133] + v[1706] * v[2145] + v[2220] * v[979]
			+ v[2227] * v[988] + v[2234] * v[996]);
	v[2761] = (*ct)*(v[1187] * v[2539] + v[1204] * v[2547] + v[1220] * v[2548]) + v[2842] + v[2412] * v[614]
		+ v[2413] * v[620] + (*cn)*(v[1708] * v[2122] + v[1707] * v[2136] + v[1706] * v[2146] + v[2221] * v[979]
			+ v[2228] * v[988] + v[2236] * v[996]);
	v[2763] = (*ct)*(v[1187] * v[2540] + v[1204] * v[2548] + v[1220] * v[2555]) + v[2843] + v[2412] * v[615]
		+ v[2413] * v[621] + (*cn)*(v[1708] * v[2125] + v[1707] * v[2137] + v[1706] * v[2147] + v[2222] * v[979]
			+ v[2229] * v[988] + v[2237] * v[996]);
	Rc[0] = v[2670];
	Rc[1] = v[2671];
	Rc[2] = v[2672];
	Rc[3] = v[1054] * v[1182] + v[1055] * v[1199] + v[1056] * v[1215] + (*ct)*(v[1182] * v[1709] + v[1199] * v[1710]
		+ v[1215] * v[1711]);
	Rc[4] = v[1054] * v[1183] + v[1055] * v[1200] + v[1056] * v[1216] + (*ct)*(v[1183] * v[1709] + v[1200] * v[1710]
		+ v[1216] * v[1711]);
	Rc[5] = v[1054] * v[1184] + v[1055] * v[1201] + v[1056] * v[1217] + (*ct)*(v[1184] * v[1709] + v[1201] * v[1710]
		+ v[1217] * v[1711]);
	Rc[6] = -v[2670];
	Rc[7] = -v[2671];
	Rc[8] = -v[2672];
	Rc[9] = v[1054] * v[1185] + v[1055] * v[1202] + v[1056] * v[1218] + (*ct)*(v[1185] * v[1709] + v[1202] * v[1710]
		+ v[1218] * v[1711]) + (*cn)*(v[1706] * v[975] + v[1707] * v[984] + v[1708] * v[992]) + (*epsn)*(v[361] * v[975]
			+ v[362] * v[984] + v[363] * v[992]);
	Rc[10] = v[1054] * v[1186] + v[1055] * v[1203] + v[1056] * v[1219] + (*ct)*(v[1186] * v[1709] + v[1203] * v[1710]
		+ v[1219] * v[1711]) + (*cn)*(v[1706] * v[977] + v[1707] * v[986] + v[1708] * v[994]) + (*epsn)*(v[361] * v[977]
			+ v[362] * v[986] + v[363] * v[994]);
	Rc[11] = v[1054] * v[1187] + v[1055] * v[1204] + v[1056] * v[1220] + (*ct)*(v[1187] * v[1709] + v[1204] * v[1710]
		+ v[1220] * v[1711]) + (*cn)*(v[1706] * v[979] + v[1707] * v[988] + v[1708] * v[996]) + (*epsn)*(v[361] * v[979]
			+ v[362] * v[988] + v[363] * v[996]);
	Kc[0][0] = v[2679];
	Kc[0][1] = v[2680];
	Kc[0][2] = v[2681];
	Kc[0][3] = v[2682];
	Kc[0][4] = v[2683];
	Kc[0][5] = v[2684];
	Kc[0][6] = -v[2679];
	Kc[0][7] = -v[2680];
	Kc[0][8] = -v[2681];
	Kc[0][9] = v[2685];
	Kc[0][10] = v[2686];
	Kc[0][11] = v[2687];
	Kc[1][0] = v[2688];
	Kc[1][1] = v[2689];
	Kc[1][2] = v[2690];
	Kc[1][3] = v[2691];
	Kc[1][4] = v[2692];
	Kc[1][5] = v[2693];
	Kc[1][6] = -v[2688];
	Kc[1][7] = -v[2689];
	Kc[1][8] = -v[2690];
	Kc[1][9] = v[2694];
	Kc[1][10] = v[2695];
	Kc[1][11] = v[2696];
	Kc[2][0] = v[2697];
	Kc[2][1] = v[2698];
	Kc[2][2] = v[2699];
	Kc[2][3] = v[2700];
	Kc[2][4] = v[2701];
	Kc[2][5] = v[2702];
	Kc[2][6] = -v[2697];
	Kc[2][7] = -v[2698];
	Kc[2][8] = -v[2699];
	Kc[2][9] = v[2703];
	Kc[2][10] = v[2704];
	Kc[2][11] = v[2705];
	Kc[3][0] = v[2706];
	Kc[3][1] = v[2707];
	Kc[3][2] = v[2708];
	Kc[3][3] = (*epst)*((v[1182] * v[1182]) + (v[1199] * v[1199]) + (v[1215] * v[1215])) + v[1056] * v[1496]
		+ v[1055] * v[1500] + v[1054] * v[1504] + (*ct)*(v[1504] * v[1709] + v[1500] * v[1710] + v[1496] * v[1711]
			+ v[1182] * v[2541] + v[1199] * v[2549] + v[1215] * v[2556]);
	Kc[3][4] = (*ct)*(v[1182] * v[2542] + v[1199] * v[2550] + v[1215] * v[2557]) + v[2844];
	Kc[3][5] = (*ct)*(v[1182] * v[2543] + v[1199] * v[2551] + v[1215] * v[2558]) + v[2845];
	Kc[3][6] = -v[2706];
	Kc[3][7] = -v[2707];
	Kc[3][8] = -v[2708];
	Kc[3][9] = (*ct)*(v[1182] * v[2544] + v[1199] * v[2552] + v[1215] * v[2559]) + v[2846];
	Kc[3][10] = (*ct)*(v[1182] * v[2545] + v[1199] * v[2553] + v[1215] * v[2560]) + v[2847];
	Kc[3][11] = (*ct)*(v[1182] * v[2546] + v[1199] * v[2554] + v[1215] * v[2561]) + v[2848];
	Kc[4][0] = v[2715];
	Kc[4][1] = v[2716];
	Kc[4][2] = v[2717];
	Kc[4][3] = (*ct)*(v[1183] * v[2541] + v[1200] * v[2549] + v[1216] * v[2556]) + v[2844];
	Kc[4][4] = (*epst)*((v[1183] * v[1183]) + (v[1200] * v[1200]) + (v[1216] * v[1216])) + v[1056] * v[1478]
		+ v[1055] * v[1483] + v[1054] * v[1488] + (*ct)*(v[1488] * v[1709] + v[1483] * v[1710] + v[1478] * v[1711]
			+ v[1183] * v[2542] + v[1200] * v[2550] + v[1216] * v[2557]);
	Kc[4][5] = (*ct)*(v[1183] * v[2543] + v[1200] * v[2551] + v[1216] * v[2558]) + v[2849];
	Kc[4][6] = -v[2715];
	Kc[4][7] = -v[2716];
	Kc[4][8] = -v[2717];
	Kc[4][9] = (*ct)*(v[1183] * v[2544] + v[1200] * v[2552] + v[1216] * v[2559]) + v[2850];
	Kc[4][10] = (*ct)*(v[1183] * v[2545] + v[1200] * v[2553] + v[1216] * v[2560]) + v[2851];
	Kc[4][11] = (*ct)*(v[1183] * v[2546] + v[1200] * v[2554] + v[1216] * v[2561]) + v[2852];
	Kc[5][0] = v[2724];
	Kc[5][1] = v[2725];
	Kc[5][2] = v[2726];
	Kc[5][3] = (*ct)*(v[1184] * v[2541] + v[1201] * v[2549] + v[1217] * v[2556]) + v[2845];
	Kc[5][4] = (*ct)*(v[1184] * v[2542] + v[1201] * v[2550] + v[1217] * v[2557]) + v[2849];
	Kc[5][5] = (*epst)*((v[1184] * v[1184]) + (v[1201] * v[1201]) + (v[1217] * v[1217])) + v[1056] * v[1456]
		+ v[1055] * v[1462] + v[1054] * v[1468] + (*ct)*(v[1468] * v[1709] + v[1462] * v[1710] + v[1456] * v[1711]
			+ v[1184] * v[2543] + v[1201] * v[2551] + v[1217] * v[2558]);
	Kc[5][6] = -v[2724];
	Kc[5][7] = -v[2725];
	Kc[5][8] = -v[2726];
	Kc[5][9] = (*ct)*(v[1184] * v[2544] + v[1201] * v[2552] + v[1217] * v[2559]) + v[2853];
	Kc[5][10] = (*ct)*(v[1184] * v[2545] + v[1201] * v[2553] + v[1217] * v[2560]) + v[2854];
	Kc[5][11] = (*ct)*(v[1184] * v[2546] + v[1201] * v[2554] + v[1217] * v[2561]) + v[2855];
	Kc[6][0] = -v[2679];
	Kc[6][1] = -v[2680];
	Kc[6][2] = -v[2681];
	Kc[6][3] = -v[2682];
	Kc[6][4] = -v[2683];
	Kc[6][5] = -v[2684];
	Kc[6][6] = v[2679];
	Kc[6][7] = v[2680];
	Kc[6][8] = v[2681];
	Kc[6][9] = -v[2685];
	Kc[6][10] = -v[2686];
	Kc[6][11] = -v[2687];
	Kc[7][0] = -v[2688];
	Kc[7][1] = -v[2689];
	Kc[7][2] = -v[2690];
	Kc[7][3] = -v[2691];
	Kc[7][4] = -v[2692];
	Kc[7][5] = -v[2693];
	Kc[7][6] = v[2688];
	Kc[7][7] = v[2689];
	Kc[7][8] = v[2690];
	Kc[7][9] = -v[2694];
	Kc[7][10] = -v[2695];
	Kc[7][11] = -v[2696];
	Kc[8][0] = -v[2697];
	Kc[8][1] = -v[2698];
	Kc[8][2] = -v[2699];
	Kc[8][3] = -v[2700];
	Kc[8][4] = -v[2701];
	Kc[8][5] = -v[2702];
	Kc[8][6] = v[2697];
	Kc[8][7] = v[2698];
	Kc[8][8] = v[2699];
	Kc[8][9] = -v[2703];
	Kc[8][10] = -v[2704];
	Kc[8][11] = -v[2705];
	Kc[9][0] = v[2734];
	Kc[9][1] = v[2736];
	Kc[9][2] = v[2738];
	Kc[9][3] = (*ct)*(v[1185] * v[2541] + v[1202] * v[2549] + v[1218] * v[2556]) + v[2846];
	Kc[9][4] = (*ct)*(v[1185] * v[2542] + v[1202] * v[2550] + v[1218] * v[2557]) + v[2850];
	Kc[9][5] = (*ct)*(v[1185] * v[2543] + v[1202] * v[2551] + v[1218] * v[2558]) + v[2853];
	Kc[9][6] = -v[2734];
	Kc[9][7] = -v[2736];
	Kc[9][8] = -v[2738];
	Kc[9][9] = (*epst)*((v[1185] * v[1185]) + (v[1202] * v[1202]) + (v[1218] * v[1218])) + v[1056] * v[1615]
		+ v[1055] * v[1621] + v[1054] * v[1627] + (*ct)*(v[1627] * v[1709] + v[1621] * v[1710] + v[1615] * v[1711]
			+ v[1185] * v[2544] + v[1202] * v[2552] + v[1218] * v[2559]) + v[2408] * v[616] + v[2409] * v[622] + (*epsn)*(
			(v[603] * v[603]) + (v[606] * v[606]) + (v[609] * v[609]) + (-(v[361] * v[549]) - v[362] * v[552] - v[363] * v[555]
				- v[343] * v[603] - v[344] * v[606] - v[345] * v[609])*v[616] + (-(v[361] * v[558]) - v[362] * v[561] - v[363] * v[564]
					- v[350] * v[603] - v[352] * v[606] - v[354] * v[609])*v[622] + v[363] * v[946] + v[362] * v[958] + v[361] * v[970]) + (*cn
						)*(v[1708] * v[2060] + v[1707] * v[2070] + v[1706] * v[2078] + v[2223] * v[975] + v[2230] * v[984] + v[2238] * v[992]);
	Kc[9][10] = (*ct)*(v[1185] * v[2545] + v[1202] * v[2553] + v[1218] * v[2560]) + v[2856] + v[2408] * v[617]
		+ v[2409] * v[623] + (*cn)*(v[1708] * v[2061] + v[1707] * v[2071] + v[1706] * v[2079] + v[2224] * v[975]
			+ v[2231] * v[984] + v[2239] * v[992]);
	Kc[9][11] = (*ct)*(v[1185] * v[2546] + v[1202] * v[2554] + v[1218] * v[2561]) + v[2857] + v[2408] * v[618]
		+ v[2409] * v[624] + (*cn)*(v[1708] * v[2062] + v[1707] * v[2072] + v[1706] * v[2080] + v[2225] * v[975]
			+ v[2232] * v[984] + v[2240] * v[992]);
	Kc[10][0] = v[2746];
	Kc[10][1] = v[2748];
	Kc[10][2] = v[2750];
	Kc[10][3] = (*ct)*(v[1186] * v[2541] + v[1203] * v[2549] + v[1219] * v[2556]) + v[2847];
	Kc[10][4] = (*ct)*(v[1186] * v[2542] + v[1203] * v[2550] + v[1219] * v[2557]) + v[2851];
	Kc[10][5] = (*ct)*(v[1186] * v[2543] + v[1203] * v[2551] + v[1219] * v[2558]) + v[2854];
	Kc[10][6] = -v[2746];
	Kc[10][7] = -v[2748];
	Kc[10][8] = -v[2750];
	Kc[10][9] = (*ct)*(v[1186] * v[2544] + v[1203] * v[2552] + v[1219] * v[2559]) + v[2856] + v[2410] * v[616]
		+ v[2411] * v[622] + (*cn)*(v[1708] * v[2093] + v[1707] * v[2104] + v[1706] * v[2113] + v[2223] * v[977]
			+ v[2230] * v[986] + v[2238] * v[994]);
	Kc[10][10] = (*epst)*((v[1186] * v[1186]) + (v[1203] * v[1203]) + (v[1219] * v[1219])) + v[1056] * v[1614]
		+ v[1055] * v[1620] + v[1054] * v[1626] + (*ct)*(v[1626] * v[1709] + v[1620] * v[1710] + v[1614] * v[1711]
			+ v[1186] * v[2545] + v[1203] * v[2553] + v[1219] * v[2560]) + v[2410] * v[617] + v[2411] * v[623] + (*epsn)*(
			(v[604] * v[604]) + (v[607] * v[607]) + (v[610] * v[610]) + v[1013] * v[617] + v[1014] * v[623] + v[363] * v[944]
				+ v[362] * v[956] + v[361] * v[968]) + (*cn)*(v[1708] * v[2094] + v[1707] * v[2105] + v[1706] * v[2114] + v[2224] * v[977]
					+ v[2231] * v[986] + v[2239] * v[994]);
	Kc[10][11] = (*ct)*(v[1186] * v[2546] + v[1203] * v[2554] + v[1219] * v[2561]) + v[2858] + v[2410] * v[618]
		+ v[2411] * v[624] + (*cn)*(v[1708] * v[2095] + v[1707] * v[2106] + v[1706] * v[2115] + v[2225] * v[977]
			+ v[2232] * v[986] + v[2240] * v[994]);
	Kc[11][0] = v[2759];
	Kc[11][1] = v[2761];
	Kc[11][2] = v[2763];
	Kc[11][3] = (*ct)*(v[1187] * v[2541] + v[1204] * v[2549] + v[1220] * v[2556]) + v[2848];
	Kc[11][4] = (*ct)*(v[1187] * v[2542] + v[1204] * v[2550] + v[1220] * v[2557]) + v[2852];
	Kc[11][5] = (*ct)*(v[1187] * v[2543] + v[1204] * v[2551] + v[1220] * v[2558]) + v[2855];
	Kc[11][6] = -v[2759];
	Kc[11][7] = -v[2761];
	Kc[11][8] = -v[2763];
	Kc[11][9] = (*ct)*(v[1187] * v[2544] + v[1204] * v[2552] + v[1220] * v[2559]) + v[2857] + v[2412] * v[616]
		+ v[2413] * v[622] + (*cn)*(v[1708] * v[2127] + v[1707] * v[2139] + v[1706] * v[2149] + v[2223] * v[979]
			+ v[2230] * v[988] + v[2238] * v[996]);
	Kc[11][10] = (*ct)*(v[1187] * v[2545] + v[1204] * v[2553] + v[1220] * v[2560]) + v[2858] + v[2412] * v[617]
		+ v[2413] * v[623] + (*cn)*(v[1708] * v[2129] + v[1707] * v[2141] + v[1706] * v[2151] + v[2224] * v[979]
			+ v[2231] * v[988] + v[2239] * v[996]);
	Kc[11][11] = (*epst)*((v[1187] * v[1187]) + (v[1204] * v[1204]) + (v[1220] * v[1220])) + v[1056] * v[1612]
		+ v[1055] * v[1618] + v[1054] * v[1624] + (*ct)*(v[1624] * v[1709] + v[1618] * v[1710] + v[1612] * v[1711]
			+ v[1187] * v[2546] + v[1204] * v[2554] + v[1220] * v[2561]) + v[2412] * v[618] + v[2413] * v[624] + (*epsn)*(
			(v[605] * v[605]) + (v[608] * v[608]) + (v[611] * v[611]) + v[1016] * v[618] + v[1017] * v[624] + v[363] * v[940]
				+ v[362] * v[952] + v[361] * v[964]) + (*cn)*(v[1708] * v[2130] + v[1707] * v[2142] + v[1706] * v[2152] + v[2225] * v[979]
					+ v[2232] * v[988] + v[2240] * v[996]);
}

//Calcula contribuições de contato entre esfera e superfície
void RigidTriangularSurface_1::ContactSphereSurfaceSliding(double* Rc, double** Kc, double zetap, double thetap, double zetai, double thetai, double* gti, int node, double* epsn, double* epst, double* cn, double* ct, double* mu, double* radius)
{
	double *dAi = dA_i->getMatrix();	//ponteiro para o vetor dA
	double *dBi = dB_i->getMatrix();	//ponteiro para o vetor dB
	double *dCi = dC_i->getMatrix();	//ponteiro para o vetor dC
	double *xPi = xP_i->getMatrix();	//ponteiro para o vetor xP
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
	v[2736] = (*epsn)*(*mu);
	v[1064] = 0.5e0*d[3];
	v[1062] = 2e0*d[3];
	v[266] = Power(d[3], 2);
	v[1065] = 2e0*d[4];
	v[1063] = 0.5e0*d[4];
	v[264] = 0.5e0*d[3] * d[4];
	v[259] = Power(d[4], 2);
	v[1073] = -v[259] - v[266];
	v[2707] = 0.5e0*v[1073];
	v[2705] = -0.5e0*v[1073];
	v[1099] = d[5] + v[264];
	v[1089] = -d[5] + v[264];
	v[1067] = 2e0*d[5];
	v[1066] = 0.5e0*d[5];
	v[271] = 0.5e0*d[4] * d[5];
	v[1111] = d[3] + v[271];
	v[1103] = -d[3] + v[271];
	v[269] = 0.5e0*d[3] * d[5];
	v[1107] = -d[4] + v[269];
	v[1095] = d[4] + v[269];
	v[260] = Power(d[5], 2);
	v[1082] = -v[259] - v[260];
	v[2711] = 0.5e0*v[1082];
	v[2710] = 0.5e0*v[1082];
	v[1078] = -v[260] - v[266];
	v[2709] = 0.5e0*v[1078];
	v[2708] = 0.5e0*v[1078];
	v[2702] = -0.5e0*v[1078];
	v[1068] = 4e0 + v[259] + v[260] + v[266];
	v[1460] = 1e0 / Power(v[1068], 3);
	v[1461] = -2e0*v[1065] * v[1460];
	v[1464] = -4e0*v[1067] * v[1461];
	v[1459] = -2e0*v[1062] * v[1460];
	v[1469] = -4e0*v[1065] * v[1459];
	v[1463] = -4e0*v[1067] * v[1459];
	v[1070] = 1e0 / Power(v[1068], 2);
	v[1470] = -8e0*v[1070];
	v[1476] = -4e0*v[1062] * v[1459] + v[1470];
	v[1471] = -4e0*v[1065] * v[1461] + v[1470];
	v[1465] = 8e0*(v[1067] * v[1067])*v[1460] + v[1470];
	v[1072] = -4e0*v[1067] * v[1070];
	v[1532] = -0.5e0*v[1067] * v[1072];
	v[1515] = 2e0*v[1072];
	v[1480] = -0.5e0*v[1062] * v[1072];
	v[1474] = -0.5e0*v[1065] * v[1072];
	v[1071] = -4e0*v[1065] * v[1070];
	v[1535] = -0.5e0*v[1065] * v[1071];
	v[1508] = v[1063] * v[1071];
	v[1499] = -2e0*v[1071];
	v[1494] = v[1066] * v[1071];
	v[1478] = -0.5e0*v[1062] * v[1071];
	v[1069] = -4e0*v[1062] * v[1070];
	v[2704] = v[1069] - v[1494];
	v[1540] = -0.5e0*v[1062] * v[1069];
	v[1511] = v[1064] * v[1069];
	v[1491] = v[1066] * v[1069];
	v[2703] = v[1071] - v[1491];
	v[1489] = 2e0*v[1069];
	v[1487] = v[1063] * v[1069];
	v[2706] = v[1072] - v[1487];
	v[472] = 0.5e0*d[9];
	v[470] = 2e0*d[9];
	v[247] = Power(d[9], 2);
	v[473] = 2e0*d[10];
	v[471] = 0.5e0*d[10];
	v[245] = 0.5e0*d[10] * d[9];
	v[240] = Power(d[10], 2);
	v[415] = -v[240] - v[247];
	v[2699] = 0.5e0*v[415];
	v[475] = 2e0*d[11];
	v[474] = 0.5e0*d[11];
	v[429] = -d[11] + v[245];
	v[425] = d[11] + v[245];
	v[252] = 0.5e0*d[10] * d[11];
	v[420] = -d[9] + v[252];
	v[418] = d[9] + v[252];
	v[250] = 0.5e0*d[11] * d[9];
	v[426] = d[10] + v[250];
	v[419] = -d[10] + v[250];
	v[241] = Power(d[11], 2);
	v[436] = 4e0 + v[240] + v[241] + v[247];
	v[688] = 1e0 / Power(v[436], 3);
	v[690] = -2e0*v[475] * v[688];
	v[2144] = -4e0*v[690];
	v[689] = -2e0*v[473] * v[688];
	v[2143] = -4e0*v[689];
	v[692] = -4e0*v[475] * v[689];
	v[687] = -2e0*v[470] * v[688];
	v[2142] = -4e0*v[687];
	v[697] = -4e0*v[473] * v[687];
	v[691] = -4e0*v[475] * v[687];
	v[476] = 1e0 / Power(v[436], 2);
	v[698] = -8e0*v[476];
	v[700] = -4e0*v[470] * v[687] + v[698];
	v[699] = -4e0*v[473] * v[689] + v[698];
	v[693] = -4e0*v[475] * v[690] + v[698];
	v[696] = 0.5e0*v[415] * v[693];
	v[479] = -4e0*v[475] * v[476];
	v[1917] = -0.5e0*v[479];
	v[1842] = -(v[474] * v[479]);
	v[786] = 2e0*v[479];
	v[787] = v[429] * v[693] - v[786];
	v[737] = v[425] * v[693] + v[786];
	v[731] = -0.5e0*v[475] * v[479];
	v[708] = -0.5e0*v[473] * v[479];
	v[821] = v[420] * v[693] - v[708];
	v[814] = v[418] * v[693] - v[708];
	v[778] = v[420] * v[699] - v[708];
	v[744] = v[418] * v[699] - v[708];
	v[703] = -0.5e0*v[470] * v[479];
	v[810] = v[426] * v[693] - v[703];
	v[789] = v[419] * v[693] - v[703];
	v[781] = v[426] * v[700] - v[703];
	v[764] = v[419] * v[700] - v[703];
	v[695] = 0.5e0*v[415] * v[692] + v[708];
	v[694] = 0.5e0*v[415] * v[691] + v[703];
	v[512] = 0.5e0*v[415] * v[479];
	v[478] = -4e0*v[473] * v[476];
	v[1904] = 0.5e0*v[478];
	v[779] = v[471] * v[478];
	v[739] = -2e0*v[478];
	v[740] = v[426] * v[699] - v[739];
	v[734] = v[474] * v[478];
	v[721] = v[419] * v[699] + v[739];
	v[712] = -0.5e0*v[473] * v[478];
	v[711] = -1e0*v[470] * v[478] + v[2699] * v[697];
	v[477] = -4e0*v[470] * v[476];
	v[2690] = -v[477] + v[734];
	v[2687] = v[477] + v[734];
	v[1910] = 0.5e0*v[477];
	v[784] = v[2690] + v[429] * v[691];
	v[782] = v[472] * v[477];
	v[738] = v[2687] + v[426] * v[697];
	v[735] = v[2687] + v[425] * v[691];
	v[726] = v[474] * v[477];
	v[2692] = v[478] + v[726];
	v[2688] = -v[478] + v[726];
	v[785] = v[2688] + v[429] * v[692];
	v[736] = v[2692] + v[425] * v[692];
	v[727] = v[2688] + v[420] * v[697];
	v[724] = 2e0*v[477];
	v[725] = v[420] * v[700] - v[724];
	v[722] = v[471] * v[477];
	v[2691] = v[479] + v[722];
	v[2689] = -v[479] + v[722];
	v[741] = v[2691] + v[426] * v[692];
	v[728] = v[2689] + v[420] * v[691];
	v[723] = v[2689] + v[419] * v[692];
	v[720] = v[2690] + v[419] * v[697];
	v[719] = v[2691] + v[418] * v[691];
	v[718] = v[2692] + v[418] * v[697];
	v[717] = v[418] * v[700] + v[724];
	v[715] = -0.5e0*v[470] * v[477];
	v[705] = -0.5e0*v[473] * v[477];
	v[841] = v[429] * v[700] - v[705];
	v[834] = v[425] * v[700] - v[705];
	v[812] = v[429] * v[699] - v[705];
	v[791] = v[425] * v[699] - v[705];
	v[432] = -v[240] - v[241];
	v[2701] = 0.5e0*v[432];
	v[860] = v[2701] * v[692] + 2e0*v[708];
	v[704] = 0.5e0*v[432] * v[691] + v[703];
	v[702] = 0.5e0*v[432] * v[697] + v[705];
	v[701] = 0.5e0*v[432] * v[700];
	v[480] = 0.5e0*v[432] * v[477];
	v[423] = -v[241] - v[247];
	v[2700] = 0.5e0*v[423];
	v[730] = v[2700] * v[691] + 2e0*v[703];
	v[709] = 0.5e0*v[423] * v[692] + v[708];
	v[707] = 0.5e0*v[423] * v[699];
	v[706] = 0.5e0*v[423] * v[697] + v[705];
	v[496] = 0.5e0*v[423] * v[478];
	v[1849] = (*a4)*d[9] + (*a6)*domegaiP[0] + (*a5)*omegaiP[0];
	v[2382] = v[1849] / 2e0;
	v[1853] = (*a4)*d[10] + (*a6)*domegaiP[1] + (*a5)*omegaiP[1];
	v[2378] = -v[1853] / 2e0;
	v[1855] = (*a4)*d[11] + (*a6)*domegaiP[2] + (*a5)*omegaiP[2];
	v[2361] = v[1855] / 2e0;
	v[866] = dAi[2] * v[741] + dAi[1] * v[785] + dAi[0] * v[860];
	v[868] = -v[866] / 2e0;
	v[848] = dAi[0] * v[701] + dAi[2] * v[781] + dAi[1] * v[841];
	v[851] = -v[848] / 2e0;
	v[826] = dAi[2] * v[728] + dAi[1] * v[730] + dAi[0] * v[735];
	v[828] = -v[826] / 2e0;
	v[817] = dAi[2] * v[696] + dAi[0] * v[789] + dAi[1] * v[814];
	v[818] = -v[817] / 2e0;
	v[799] = dAi[1] * v[707] + dAi[2] * v[778] + dAi[0] * v[791];
	v[802] = -v[799] / 2e0;
	v[752] = dAi[2] * v[711] + dAi[1] * v[718] + dAi[0] * v[720];
	v[755] = -v[752] / 2e0;
	v[864] = dBi[2] * v[741] + dBi[1] * v[785] + dBi[0] * v[860];
	v[845] = dBi[0] * v[701] + dBi[2] * v[781] + dBi[1] * v[841];
	v[824] = dBi[2] * v[728] + dBi[1] * v[730] + dBi[0] * v[735];
	v[816] = dBi[2] * v[696] + dBi[0] * v[789] + dBi[1] * v[814];
	v[796] = dBi[1] * v[707] + dBi[2] * v[778] + dBi[0] * v[791];
	v[749] = dBi[2] * v[711] + dBi[1] * v[718] + dBi[0] * v[720];
	v[862] = dCi[2] * v[741] + dCi[1] * v[785] + dCi[0] * v[860];
	v[842] = dCi[0] * v[701] + dCi[2] * v[781] + dCi[1] * v[841];
	v[822] = dCi[2] * v[728] + dCi[1] * v[730] + dCi[0] * v[735];
	v[815] = dCi[2] * v[696] + dCi[0] * v[789] + dCi[1] * v[814];
	v[793] = dCi[1] * v[707] + dCi[2] * v[778] + dCi[0] * v[791];
	v[746] = dCi[2] * v[711] + dCi[1] * v[718] + dCi[0] * v[720];
	v[224] = d[6] + xPi[0];
	v[225] = d[7] + xPi[1];
	v[226] = d[8] + xPi[2];
	v[239] = 4e0 / v[436];
	v[1846] = (v[239] * v[239]);
	v[742] = -0.5e0*v[239];
	v[2693] = -v[742] + v[779] + v[782];
	v[811] = v[2693] + v[429] * v[697];
	v[849] = dAi[0] * v[702] + dAi[2] * v[738] + dAi[1] * v[811];
	v[852] = -v[849] / 2e0;
	v[846] = dBi[0] * v[702] + dBi[2] * v[738] + dBi[1] * v[811];
	v[843] = dCi[0] * v[702] + dCi[2] * v[738] + dCi[1] * v[811];
	v[790] = v[2693] + v[425] * v[697];
	v[798] = dAi[1] * v[706] + dAi[2] * v[727] + dAi[0] * v[790];
	v[801] = -v[798] / 2e0;
	v[795] = dBi[1] * v[706] + dBi[2] * v[727] + dBi[0] * v[790];
	v[792] = dCi[1] * v[706] + dCi[2] * v[727] + dCi[0] * v[790];
	v[743] = v[1842] + v[742];
	v[2695] = -v[743] + v[779];
	v[2694] = -v[743] + v[782];
	v[783] = v[2694] + v[426] * v[691];
	v[850] = dAi[0] * v[704] + dAi[2] * v[783] + dAi[1] * v[784];
	v[853] = -v[850] / 2e0;
	v[847] = dBi[0] * v[704] + dBi[2] * v[783] + dBi[1] * v[784];
	v[844] = dCi[0] * v[704] + dCi[2] * v[783] + dCi[1] * v[784];
	v[780] = v[2695] + v[420] * v[692];
	v[800] = dAi[1] * v[709] + dAi[0] * v[736] + dAi[2] * v[780];
	v[803] = -v[800] / 2e0;
	v[797] = dBi[1] * v[709] + dBi[0] * v[736] + dBi[2] * v[780];
	v[794] = dCi[1] * v[709] + dCi[0] * v[736] + dCi[2] * v[780];
	v[765] = v[2694] + v[419] * v[691];
	v[771] = dAi[2] * v[694] + dAi[1] * v[719] + dAi[0] * v[765];
	v[773] = -v[771] / 2e0;
	v[769] = dBi[2] * v[694] + dBi[1] * v[719] + dBi[0] * v[765];
	v[767] = dCi[2] * v[694] + dCi[1] * v[719] + dCi[0] * v[765];
	v[745] = v[2695] + v[418] * v[692];
	v[754] = dAi[2] * v[695] + dAi[0] * v[723] + dAi[1] * v[745];
	v[757] = -v[754] / 2e0;
	v[751] = dBi[2] * v[695] + dBi[0] * v[723] + dBi[1] * v[745];
	v[748] = dCi[2] * v[695] + dCi[0] * v[723] + dCi[1] * v[745];
	v[2696] = -v[239] + 2e0*v[731];
	v[861] = v[2696] + 0.5e0*v[432] * v[693];
	v[867] = dAi[1] * v[787] + dAi[2] * v[810] + dAi[0] * v[861];
	v[869] = -v[867] / 2e0;
	v[865] = dBi[1] * v[787] + dBi[2] * v[810] + dBi[0] * v[861];
	v[863] = dCi[1] * v[787] + dCi[2] * v[810] + dCi[0] * v[861];
	v[732] = v[2696] + 0.5e0*v[423] * v[693];
	v[827] = dAi[1] * v[732] + dAi[0] * v[737] + dAi[2] * v[821];
	v[829] = -v[827] / 2e0;
	v[825] = dBi[1] * v[732] + dBi[0] * v[737] + dBi[2] * v[821];
	v[823] = dCi[1] * v[732] + dCi[0] * v[737] + dCi[2] * v[821];
	v[2697] = -v[239] + 2e0*v[715];
	v[733] = v[2697] + 0.5e0*v[423] * v[700];
	v[837] = dAi[2] * v[725] + dAi[1] * v[733] + dAi[0] * v[834];
	v[838] = -v[837] / 2e0;
	v[836] = dBi[2] * v[725] + dBi[1] * v[733] + dBi[0] * v[834];
	v[835] = dCi[2] * v[725] + dCi[1] * v[733] + dCi[0] * v[834];
	v[716] = v[2697] + 0.5e0*v[415] * v[700];
	v[770] = dAi[2] * v[716] + dAi[1] * v[717] + dAi[0] * v[764];
	v[772] = -v[770] / 2e0;
	v[768] = dBi[2] * v[716] + dBi[1] * v[717] + dBi[0] * v[764];
	v[766] = dCi[2] * v[716] + dCi[1] * v[717] + dCi[0] * v[764];
	v[2698] = -v[239] + 2e0*v[712];
	v[874] = v[2698] + 0.5e0*v[432] * v[699];
	v[877] = dAi[2] * v[740] + dAi[1] * v[812] + dAi[0] * v[874];
	v[878] = -v[877] / 2e0;
	v[876] = dBi[2] * v[740] + dBi[1] * v[812] + dBi[0] * v[874];
	v[875] = dCi[2] * v[740] + dCi[1] * v[812] + dCi[0] * v[874];
	v[713] = v[2698] + 0.5e0*v[415] * v[699];
	v[753] = dAi[2] * v[713] + dAi[0] * v[721] + dAi[1] * v[744];
	v[756] = -v[753] / 2e0;
	v[750] = dBi[2] * v[713] + dBi[0] * v[721] + dBi[1] * v[744];
	v[747] = dCi[2] * v[713] + dCi[0] * v[721] + dCi[1] * v[744];
	v[510] = -0.5e0*v[239] * v[473];
	v[511] = v[2699] * v[478] + v[510];
	v[508] = -0.5e0*v[239] * v[470];
	v[509] = 0.5e0*v[415] * v[477] + v[508];
	v[505] = v[239] + v[418] * v[477];
	v[503] = -v[239] + v[419] * v[478];
	v[499] = -v[239] + v[420] * v[477];
	v[497] = -0.5e0*v[239] * v[475];
	v[498] = v[2700] * v[479] + v[497];
	v[495] = 0.5e0*v[423] * v[477] + v[508];
	v[494] = v[239] + v[425] * v[479];
	v[490] = v[239] + v[426] * v[478];
	v[488] = -(v[239] * v[474]);
	v[506] = v[418] * v[478] - v[488];
	v[538] = dCi[0] * v[503] + dCi[1] * v[506] + dCi[2] * v[511];
	v[529] = dBi[0] * v[503] + dBi[1] * v[506] + dBi[2] * v[511];
	v[520] = dAi[0] * v[503] + dAi[1] * v[506] + dAi[2] * v[511];
	v[547] = -v[520] / 2e0;
	v[565] = v[538] / 2e0 + v[547];
	v[556] = v[529] / 2e0 + v[547];
	v[502] = v[419] * v[477] - v[488];
	v[537] = dCi[0] * v[502] + dCi[1] * v[505] + dCi[2] * v[509];
	v[528] = dBi[0] * v[502] + dBi[1] * v[505] + dBi[2] * v[509];
	v[519] = dAi[0] * v[502] + dAi[1] * v[505] + dAi[2] * v[509];
	v[546] = -v[519] / 2e0;
	v[564] = v[537] / 2e0 + v[546];
	v[555] = v[528] / 2e0 + v[546];
	v[500] = v[420] * v[478] - v[488];
	v[489] = v[426] * v[477] - v[488];
	v[487] = -v[239] + v[429] * v[479];
	v[485] = -(v[239] * v[472]);
	v[504] = v[419] * v[479] - v[485];
	v[493] = v[425] * v[478] - v[485];
	v[535] = dCi[0] * v[493] + dCi[1] * v[496] + dCi[2] * v[500];
	v[526] = dBi[0] * v[493] + dBi[1] * v[496] + dBi[2] * v[500];
	v[517] = dAi[0] * v[493] + dAi[1] * v[496] + dAi[2] * v[500];
	v[544] = -v[517] / 2e0;
	v[562] = v[535] / 2e0 + v[544];
	v[553] = v[526] / 2e0 + v[544];
	v[491] = v[426] * v[479] - v[485];
	v[486] = v[429] * v[478] - v[485];
	v[483] = v[239] * v[471];
	v[507] = v[418] * v[479] + v[483];
	v[539] = dCi[0] * v[504] + dCi[1] * v[507] + dCi[2] * v[512];
	v[530] = dBi[0] * v[504] + dBi[1] * v[507] + dBi[2] * v[512];
	v[521] = dAi[0] * v[504] + dAi[1] * v[507] + dAi[2] * v[512];
	v[548] = -v[521] / 2e0;
	v[566] = v[539] / 2e0 + v[548];
	v[557] = v[530] / 2e0 + v[548];
	v[501] = v[420] * v[479] + v[483];
	v[536] = dCi[0] * v[494] + dCi[1] * v[498] + dCi[2] * v[501];
	v[527] = dBi[0] * v[494] + dBi[1] * v[498] + dBi[2] * v[501];
	v[518] = dAi[0] * v[494] + dAi[1] * v[498] + dAi[2] * v[501];
	v[545] = -v[518] / 2e0;
	v[563] = v[536] / 2e0 + v[545];
	v[554] = v[527] / 2e0 + v[545];
	v[492] = v[425] * v[477] + v[483];
	v[534] = dCi[0] * v[492] + dCi[1] * v[495] + dCi[2] * v[499];
	v[525] = dBi[0] * v[492] + dBi[1] * v[495] + dBi[2] * v[499];
	v[516] = dAi[0] * v[492] + dAi[1] * v[495] + dAi[2] * v[499];
	v[543] = -v[516] / 2e0;
	v[561] = v[534] / 2e0 + v[543];
	v[552] = v[525] / 2e0 + v[543];
	v[484] = v[429] * v[477] + v[483];
	v[531] = dCi[0] * v[480] + dCi[1] * v[484] + dCi[2] * v[489];
	v[522] = dBi[0] * v[480] + dBi[1] * v[484] + dBi[2] * v[489];
	v[513] = dAi[0] * v[480] + dAi[1] * v[484] + dAi[2] * v[489];
	v[540] = -v[513] / 2e0;
	v[558] = v[531] / 2e0 + v[540];
	v[549] = v[522] / 2e0 + v[540];
	v[482] = v[2701] * v[479] + v[497];
	v[533] = dCi[0] * v[482] + dCi[1] * v[487] + dCi[2] * v[491];
	v[524] = dBi[0] * v[482] + dBi[1] * v[487] + dBi[2] * v[491];
	v[515] = dAi[0] * v[482] + dAi[1] * v[487] + dAi[2] * v[491];
	v[542] = -v[515] / 2e0;
	v[560] = v[533] / 2e0 + v[542];
	v[551] = v[524] / 2e0 + v[542];
	v[481] = 0.5e0*v[432] * v[478] + v[510];
	v[532] = dCi[0] * v[481] + dCi[1] * v[486] + dCi[2] * v[490];
	v[523] = dBi[0] * v[481] + dBi[1] * v[486] + dBi[2] * v[490];
	v[514] = dAi[0] * v[481] + dAi[1] * v[486] + dAi[2] * v[490];
	v[541] = -v[514] / 2e0;
	v[559] = v[532] / 2e0 + v[541];
	v[550] = v[523] / 2e0 + v[541];
	v[242] = 1e0 + 0.5e0*v[239] * v[432];
	v[243] = v[239] * v[429];
	v[244] = v[239] * v[426];
	v[246] = v[239] * v[425];
	v[248] = 1e0 + 0.5e0*v[239] * v[423];
	v[249] = v[239] * v[420];
	v[251] = v[239] * v[419];
	v[253] = v[239] * v[418];
	v[254] = 1e0 + 0.5e0*v[239] * v[415];
	v[258] = 4e0 / v[1068];
	v[2730] = (*radius)*v[258];
	v[1539] = v[1540] - v[258];
	v[1534] = v[1535] - v[258];
	v[1530] = v[1532] - v[258];
	v[1517] = -0.5e0*v[258];
	v[1518] = -v[1511] + v[1517];
	v[1502] = -(v[1066] * v[1072]) + v[1517];
	v[1094] = -(v[1066] * v[258]);
	v[1091] = -(v[1064] * v[258]);
	v[1088] = v[1063] * v[258];
	v[1086] = -0.5e0*v[1067] * v[258];
	v[1084] = -0.5e0*v[1065] * v[258];
	v[1077] = -0.5e0*v[1062] * v[258];
	v[286] = dAi[0] + xPi[0];
	v[326] = -v[286] / 2e0;
	v[287] = dAi[1] + xPi[1];
	v[328] = -v[287] / 2e0;
	v[288] = dAi[2] + xPi[2];
	v[330] = -v[288] / 2e0;
	v[289] = dBi[0] + xPi[0];
	v[319] = v[289] / 2e0 + v[326];
	v[290] = dBi[1] + xPi[1];
	v[320] = v[290] / 2e0 + v[328];
	v[291] = dBi[2] + xPi[2];
	v[321] = v[291] / 2e0 + v[330];
	v[323] = 1e0 / sqrt((v[319] * v[319]) + (v[320] * v[320]) + (v[321] * v[321]));
	v[898] = v[323] * (v[845] / 2e0 + v[851]);
	v[897] = v[323] * (v[876] / 2e0 + v[878]);
	v[896] = v[323] * (v[846] / 2e0 + v[852]);
	v[895] = v[323] * (v[865] / 2e0 + v[869]);
	v[894] = v[323] * (v[864] / 2e0 + v[868]);
	v[893] = v[323] * (v[847] / 2e0 + v[853]);
	v[892] = v[323] * (v[836] / 2e0 + v[838]);
	v[891] = v[323] * (v[796] / 2e0 + v[802]);
	v[890] = v[323] * (v[795] / 2e0 + v[801]);
	v[889] = v[323] * (v[825] / 2e0 + v[829]);
	v[888] = v[323] * (v[797] / 2e0 + v[803]);
	v[887] = v[323] * (v[824] / 2e0 + v[828]);
	v[886] = v[323] * (v[768] / 2e0 + v[772]);
	v[885] = v[323] * (v[750] / 2e0 + v[756]);
	v[884] = v[323] * (v[749] / 2e0 + v[755]);
	v[883] = v[323] * (v[816] / 2e0 + v[818]);
	v[882] = v[323] * (v[751] / 2e0 + v[757]);
	v[881] = v[323] * (v[769] / 2e0 + v[773]);
	v[584] = v[323] * v[557];
	v[583] = v[323] * v[556];
	v[582] = v[323] * v[555];
	v[581] = v[323] * v[554];
	v[580] = v[323] * v[553];
	v[579] = v[323] * v[552];
	v[578] = v[323] * v[551];
	v[577] = v[323] * v[550];
	v[576] = v[323] * v[549];
	v[292] = dCi[0] + xPi[0];
	v[327] = v[292] / 2e0 + v[326];
	v[293] = dCi[1] + xPi[1];
	v[329] = v[293] / 2e0 + v[328];
	v[294] = dCi[2] + xPi[2];
	v[331] = v[294] / 2e0 + v[330];
	v[333] = 1e0 / sqrt((v[327] * v[327]) + (v[329] * v[329]) + (v[331] * v[331]));
	v[916] = v[333] * (v[842] / 2e0 + v[851]);
	v[915] = v[333] * (v[875] / 2e0 + v[878]);
	v[914] = v[333] * (v[843] / 2e0 + v[852]);
	v[913] = v[333] * (v[863] / 2e0 + v[869]);
	v[912] = v[333] * (v[862] / 2e0 + v[868]);
	v[911] = v[333] * (v[844] / 2e0 + v[853]);
	v[910] = v[333] * (v[835] / 2e0 + v[838]);
	v[909] = v[333] * (v[793] / 2e0 + v[802]);
	v[908] = v[333] * (v[792] / 2e0 + v[801]);
	v[907] = v[333] * (v[823] / 2e0 + v[829]);
	v[906] = v[333] * (v[794] / 2e0 + v[803]);
	v[905] = v[333] * (v[822] / 2e0 + v[828]);
	v[904] = v[333] * (v[766] / 2e0 + v[772]);
	v[903] = v[333] * (v[747] / 2e0 + v[756]);
	v[902] = v[333] * (v[746] / 2e0 + v[755]);
	v[901] = v[333] * (v[815] / 2e0 + v[818]);
	v[900] = v[333] * (v[748] / 2e0 + v[757]);
	v[899] = v[333] * (v[767] / 2e0 + v[773]);
	v[593] = v[333] * v[566];
	v[592] = v[333] * v[565];
	v[591] = v[333] * v[564];
	v[590] = v[333] * v[563];
	v[589] = v[333] * v[562];
	v[588] = v[333] * v[561];
	v[587] = v[333] * v[560];
	v[586] = v[333] * v[559];
	v[585] = v[333] * v[558];
	v[295] = v[224] + dAi[0] * v[242] + dAi[1] * v[243] + dAi[2] * v[244];
	v[349] = -v[295] / 2e0;
	v[296] = v[225] + dAi[0] * v[246] + dAi[1] * v[248] + dAi[2] * v[249];
	v[351] = -v[296] / 2e0;
	v[297] = v[226] + dAi[0] * v[251] + dAi[1] * v[253] + dAi[2] * v[254];
	v[353] = -v[297] / 2e0;
	v[298] = v[224] + dBi[0] * v[242] + dBi[1] * v[243] + dBi[2] * v[244];
	v[343] = v[298] / 2e0 + v[349];
	v[299] = v[225] + dBi[0] * v[246] + dBi[1] * v[248] + dBi[2] * v[249];
	v[344] = v[299] / 2e0 + v[351];
	v[300] = v[226] + dBi[0] * v[251] + dBi[1] * v[253] + dBi[2] * v[254];
	v[345] = v[300] / 2e0 + v[353];
	v[301] = v[224] + dCi[0] * v[242] + dCi[1] * v[243] + dCi[2] * v[244];
	v[350] = v[301] / 2e0 + v[349];
	v[302] = v[225] + dCi[0] * v[246] + dCi[1] * v[248] + dCi[2] * v[249];
	v[352] = v[302] / 2e0 + v[351];
	v[303] = v[226] + dCi[0] * v[251] + dCi[1] * v[253] + dCi[2] * v[254];
	v[354] = v[303] / 2e0 + v[353];
	v[304] = (-ci[0] - ci[1]) / 2e0;
	v[305] = (1e0 + ci[0]) / 2e0;
	v[306] = (1e0 + ci[1]) / 2e0;
	v[1562] = v[306] * v[842] + v[305] * v[845] + v[304] * v[848];
	v[1561] = v[306] * v[875] + v[305] * v[876] + v[304] * v[877];
	v[1560] = v[306] * v[843] + v[305] * v[846] + v[304] * v[849];
	v[1559] = v[306] * v[863] + v[305] * v[865] + v[304] * v[867];
	v[1558] = v[306] * v[862] + v[305] * v[864] + v[304] * v[866];
	v[1557] = v[306] * v[844] + v[305] * v[847] + v[304] * v[850];
	v[1556] = v[306] * v[835] + v[305] * v[836] + v[304] * v[837];
	v[1555] = v[306] * v[793] + v[305] * v[796] + v[304] * v[799];
	v[1554] = v[306] * v[792] + v[305] * v[795] + v[304] * v[798];
	v[1553] = v[306] * v[823] + v[305] * v[825] + v[304] * v[827];
	v[1552] = v[306] * v[794] + v[305] * v[797] + v[304] * v[800];
	v[1551] = v[306] * v[822] + v[305] * v[824] + v[304] * v[826];
	v[1550] = v[306] * v[766] + v[305] * v[768] + v[304] * v[770];
	v[1549] = v[306] * v[747] + v[305] * v[750] + v[304] * v[753];
	v[1548] = v[306] * v[746] + v[305] * v[749] + v[304] * v[752];
	v[1547] = v[306] * v[815] + v[305] * v[816] + v[304] * v[817];
	v[1546] = v[306] * v[748] + v[305] * v[751] + v[304] * v[754];
	v[1545] = v[306] * v[767] + v[305] * v[769] + v[304] * v[771];
	v[1123] = v[304] * v[521] + v[305] * v[530] + v[306] * v[539];
	v[1122] = v[304] * v[520] + v[305] * v[529] + v[306] * v[538];
	v[1121] = v[304] * v[519] + v[305] * v[528] + v[306] * v[537];
	v[1120] = v[304] * v[518] + v[305] * v[527] + v[306] * v[536];
	v[1119] = v[304] * v[517] + v[305] * v[526] + v[306] * v[535];
	v[1118] = v[304] * v[516] + v[305] * v[525] + v[306] * v[534];
	v[1117] = v[304] * v[515] + v[305] * v[524] + v[306] * v[533];
	v[1116] = v[304] * v[514] + v[305] * v[523] + v[306] * v[532];
	v[1115] = v[304] * v[513] + v[305] * v[522] + v[306] * v[531];
	v[307] = (-cp[0] - cp[1]) / 2e0;
	v[308] = (1e0 + cp[0]) / 2e0;
	v[1973] = -(v[308] * v[587]);
	v[1972] = -(v[308] * v[586]);
	v[1971] = -(v[308] * v[585]);
	v[1961] = -(v[308] * v[590]);
	v[1960] = -(v[308] * v[589]);
	v[1959] = -(v[308] * v[588]);
	v[1949] = -(v[308] * v[593]);
	v[1948] = -(v[308] * v[592]);
	v[1947] = -(v[308] * v[591]);
	v[309] = (1e0 + cp[1]) / 2e0;
	v[1964] = -(v[309] * v[578]);
	v[1963] = -(v[309] * v[577]);
	v[1962] = -(v[309] * v[576]);
	v[1952] = -(v[309] * v[581]);
	v[1951] = -(v[309] * v[580]);
	v[1950] = -(v[309] * v[579]);
	v[1940] = -(v[309] * v[584]);
	v[1939] = -(v[309] * v[583]);
	v[1938] = -(v[309] * v[582]);
	v[575] = v[307] * v[521] + v[308] * v[530] + v[309] * v[539];
	v[1937] = -(v[323] * v[575]) / 2e0;
	v[1976] = -v[1937] - v[307] * v[584];
	v[1946] = v[1937] - v[308] * v[584];
	v[2002] = dCi[0] * v[1940] + dBi[0] * v[1946] + dAi[0] * v[1976];
	v[1988] = dCi[1] * v[1940] + dBi[1] * v[1946] + dAi[1] * v[1976];
	v[1982] = dCi[2] * v[1940] + dBi[2] * v[1946] + dAi[2] * v[1976];
	v[1934] = -(v[333] * v[575]) / 2e0;
	v[1979] = -v[1934] - v[307] * v[593];
	v[1943] = v[1934] - v[309] * v[593];
	v[2009] = dCi[0] * v[1943] + dBi[0] * v[1949] + dAi[0] * v[1979];
	v[1995] = dCi[1] * v[1943] + dBi[1] * v[1949] + dAi[1] * v[1979];
	v[1985] = dCi[2] * v[1943] + dBi[2] * v[1949] + dAi[2] * v[1979];
	v[574] = v[307] * v[520] + v[308] * v[529] + v[309] * v[538];
	v[1936] = -(v[323] * v[574]) / 2e0;
	v[1975] = -v[1936] - v[307] * v[583];
	v[1945] = v[1936] - v[308] * v[583];
	v[2001] = dCi[0] * v[1939] + dBi[0] * v[1945] + dAi[0] * v[1975];
	v[1987] = dCi[1] * v[1939] + dBi[1] * v[1945] + dAi[1] * v[1975];
	v[1981] = dCi[2] * v[1939] + dBi[2] * v[1945] + dAi[2] * v[1975];
	v[1933] = -(v[333] * v[574]) / 2e0;
	v[1978] = -v[1933] - v[307] * v[592];
	v[1942] = v[1933] - v[309] * v[592];
	v[2008] = dCi[0] * v[1942] + dBi[0] * v[1948] + dAi[0] * v[1978];
	v[1994] = dCi[1] * v[1942] + dBi[1] * v[1948] + dAi[1] * v[1978];
	v[1984] = dCi[2] * v[1942] + dBi[2] * v[1948] + dAi[2] * v[1978];
	v[573] = v[307] * v[519] + v[308] * v[528] + v[309] * v[537];
	v[1935] = -(v[323] * v[573]) / 2e0;
	v[1974] = -v[1935] - v[307] * v[582];
	v[1944] = v[1935] - v[308] * v[582];
	v[2000] = dCi[0] * v[1938] + dBi[0] * v[1944] + dAi[0] * v[1974];
	v[1986] = dCi[1] * v[1938] + dBi[1] * v[1944] + dAi[1] * v[1974];
	v[1980] = dCi[2] * v[1938] + dBi[2] * v[1944] + dAi[2] * v[1974];
	v[1932] = -(v[333] * v[573]) / 2e0;
	v[1977] = -v[1932] - v[307] * v[591];
	v[1941] = v[1932] - v[309] * v[591];
	v[2007] = dCi[0] * v[1941] + dBi[0] * v[1947] + dAi[0] * v[1977];
	v[1993] = dCi[1] * v[1941] + dBi[1] * v[1947] + dAi[1] * v[1977];
	v[1983] = dCi[2] * v[1941] + dBi[2] * v[1947] + dAi[2] * v[1977];
	v[572] = v[307] * v[518] + v[308] * v[527] + v[309] * v[536];
	v[1931] = -(v[323] * v[572]) / 2e0;
	v[2016] = -v[1931] - v[307] * v[581];
	v[1958] = v[1931] - v[308] * v[581];
	v[2042] = dCi[0] * v[1952] + dBi[0] * v[1958] + dAi[0] * v[2016];
	v[2036] = dCi[1] * v[1952] + dBi[1] * v[1958] + dAi[1] * v[2016];
	v[2022] = dCi[2] * v[1952] + dBi[2] * v[1958] + dAi[2] * v[2016];
	v[1928] = -(v[333] * v[572]) / 2e0;
	v[2019] = -v[1928] - v[307] * v[590];
	v[1955] = v[1928] - v[309] * v[590];
	v[2048] = dCi[0] * v[1955] + dBi[0] * v[1961] + dAi[0] * v[2019];
	v[2039] = dCi[1] * v[1955] + dBi[1] * v[1961] + dAi[1] * v[2019];
	v[2029] = dCi[2] * v[1955] + dBi[2] * v[1961] + dAi[2] * v[2019];
	v[571] = v[307] * v[517] + v[308] * v[526] + v[309] * v[535];
	v[1930] = -(v[323] * v[571]) / 2e0;
	v[2015] = -v[1930] - v[307] * v[580];
	v[1957] = v[1930] - v[308] * v[580];
	v[2041] = dCi[0] * v[1951] + dBi[0] * v[1957] + dAi[0] * v[2015];
	v[2035] = dCi[1] * v[1951] + dBi[1] * v[1957] + dAi[1] * v[2015];
	v[2021] = dCi[2] * v[1951] + dBi[2] * v[1957] + dAi[2] * v[2015];
	v[1927] = -(v[333] * v[571]) / 2e0;
	v[2018] = -v[1927] - v[307] * v[589];
	v[1954] = v[1927] - v[309] * v[589];
	v[2047] = dCi[0] * v[1954] + dBi[0] * v[1960] + dAi[0] * v[2018];
	v[2038] = dCi[1] * v[1954] + dBi[1] * v[1960] + dAi[1] * v[2018];
	v[2028] = dCi[2] * v[1954] + dBi[2] * v[1960] + dAi[2] * v[2018];
	v[570] = v[307] * v[516] + v[308] * v[525] + v[309] * v[534];
	v[1929] = -(v[323] * v[570]) / 2e0;
	v[2014] = -v[1929] - v[307] * v[579];
	v[1956] = v[1929] - v[308] * v[579];
	v[2040] = dCi[0] * v[1950] + dBi[0] * v[1956] + dAi[0] * v[2014];
	v[2034] = dCi[1] * v[1950] + dBi[1] * v[1956] + dAi[1] * v[2014];
	v[2020] = dCi[2] * v[1950] + dBi[2] * v[1956] + dAi[2] * v[2014];
	v[1926] = -(v[333] * v[570]) / 2e0;
	v[2017] = -v[1926] - v[307] * v[588];
	v[1953] = v[1926] - v[309] * v[588];
	v[2046] = dCi[0] * v[1953] + dBi[0] * v[1959] + dAi[0] * v[2017];
	v[2037] = dCi[1] * v[1953] + dBi[1] * v[1959] + dAi[1] * v[2017];
	v[2027] = dCi[2] * v[1953] + dBi[2] * v[1959] + dAi[2] * v[2017];
	v[569] = v[307] * v[515] + v[308] * v[524] + v[309] * v[533];
	v[1925] = -(v[323] * v[569]) / 2e0;
	v[2054] = -v[1925] - v[307] * v[578];
	v[1970] = v[1925] - v[308] * v[578];
	v[2084] = dCi[0] * v[1964] + dBi[0] * v[1970] + dAi[0] * v[2054];
	v[2072] = dCi[1] * v[1964] + dBi[1] * v[1970] + dAi[1] * v[2054];
	v[2060] = dCi[2] * v[1964] + dBi[2] * v[1970] + dAi[2] * v[2054];
	v[1921] = -(v[333] * v[569]) / 2e0;
	v[2057] = -v[1921] - v[307] * v[587];
	v[1967] = v[1921] - v[309] * v[587];
	v[2087] = dCi[0] * v[1967] + dBi[0] * v[1973] + dAi[0] * v[2057];
	v[2078] = dCi[1] * v[1967] + dBi[1] * v[1973] + dAi[1] * v[2057];
	v[2066] = dCi[2] * v[1967] + dBi[2] * v[1973] + dAi[2] * v[2057];
	v[568] = v[307] * v[514] + v[308] * v[523] + v[309] * v[532];
	v[1924] = -(v[323] * v[568]) / 2e0;
	v[2053] = -v[1924] - v[307] * v[577];
	v[1969] = v[1924] - v[308] * v[577];
	v[2083] = dCi[0] * v[1963] + dBi[0] * v[1969] + dAi[0] * v[2053];
	v[2071] = dCi[1] * v[1963] + dBi[1] * v[1969] + dAi[1] * v[2053];
	v[2059] = dCi[2] * v[1963] + dBi[2] * v[1969] + dAi[2] * v[2053];
	v[1920] = -(v[333] * v[568]) / 2e0;
	v[2056] = -v[1920] - v[307] * v[586];
	v[1966] = v[1920] - v[309] * v[586];
	v[2086] = dCi[0] * v[1966] + dBi[0] * v[1972] + dAi[0] * v[2056];
	v[2077] = dCi[1] * v[1966] + dBi[1] * v[1972] + dAi[1] * v[2056];
	v[2065] = dCi[2] * v[1966] + dBi[2] * v[1972] + dAi[2] * v[2056];
	v[567] = v[307] * v[513] + v[308] * v[522] + v[309] * v[531];
	v[1923] = -(v[323] * v[567]) / 2e0;
	v[2052] = -v[1923] - v[307] * v[576];
	v[1968] = v[1923] - v[308] * v[576];
	v[2070] = dCi[1] * v[1962] + dBi[1] * v[1968] + dAi[1] * v[2052];
	v[2058] = dCi[2] * v[1962] + dBi[2] * v[1968] + dAi[2] * v[2052];
	v[1919] = -(v[333] * v[567]) / 2e0;
	v[2055] = -v[1919] - v[307] * v[585];
	v[1965] = v[1919] - v[309] * v[585];
	v[2076] = dCi[1] * v[1965] + dBi[1] * v[1971] + dAi[1] * v[2055];
	v[2064] = dCi[2] * v[1965] + dBi[2] * v[1971] + dAi[2] * v[2055];
	v[365] = d[0] - v[295] * v[307] - v[298] * v[308] - v[301] * v[309] + xSi[0];
	v[652] = -v[365] / 2e0;
	v[407] = (v[333] * v[365]) / 2e0;
	v[405] = (v[323] * v[365]) / 2e0;
	v[366] = d[1] - v[296] * v[307] - v[299] * v[308] - v[302] * v[309] + xSi[1];
	v[649] = -v[366] / 2e0;
	v[397] = (v[333] * v[366]) / 2e0;
	v[395] = (v[323] * v[366]) / 2e0;
	v[367] = d[2] - v[297] * v[307] - v[300] * v[308] - v[303] * v[309] + xSi[2];
	v[646] = -v[367] / 2e0;
	v[387] = (v[333] * v[367]) / 2e0;
	v[385] = (v[323] * v[367]) / 2e0;
	v[322] = v[319] * v[323];
	v[324] = v[320] * v[323];
	v[325] = v[321] * v[323];
	v[1568] = v[325] * v[886] + v[324] * v[892] + v[322] * v[898];
	v[1567] = v[325] * v[885] + v[324] * v[891] + v[322] * v[897];
	v[1566] = v[325] * v[884] + v[324] * v[890] + v[322] * v[896];
	v[1565] = v[325] * v[883] + v[324] * v[889] + v[322] * v[895];
	v[1564] = v[325] * v[882] + v[324] * v[888] + v[322] * v[894];
	v[1563] = v[325] * v[881] + v[324] * v[887] + v[322] * v[893];
	v[1135] = v[322] * v[578] + v[324] * v[581] + v[325] * v[584];
	v[1134] = v[322] * v[577] + v[324] * v[580] + v[325] * v[583];
	v[1133] = v[322] * v[576] + v[324] * v[579] + v[325] * v[582];
	v[332] = v[327] * v[333];
	v[334] = v[329] * v[333];
	v[341] = -(v[324] * v[332]) + v[322] * v[334];
	v[335] = v[331] * v[333];
	v[339] = v[325] * v[332] - v[322] * v[335];
	v[336] = -(v[325] * v[334]) + v[324] * v[335];
	v[338] = 1e0 / sqrt((v[336] * v[336]) + (v[339] * v[339]) + (v[341] * v[341]));
	v[337] = v[336] * v[338];
	v[2729] = -(v[2711] * v[337]);
	v[340] = v[338] * v[339];
	v[2732] = -(v[2708] * v[340]);
	v[342] = v[338] * v[341];
	v[2734] = -(v[2707] * v[342]);
	v[1592] = v[342] * v[886] + v[340] * v[892] + v[337] * v[898];
	v[1591] = v[342] * v[885] + v[340] * v[891] + v[337] * v[897];
	v[1590] = v[342] * v[884] + v[340] * v[890] + v[337] * v[896];
	v[1589] = v[342] * v[883] + v[340] * v[889] + v[337] * v[895];
	v[1588] = v[342] * v[882] + v[340] * v[888] + v[337] * v[894];
	v[1587] = v[342] * v[881] + v[340] * v[887] + v[337] * v[893];
	v[1586] = (*radius)*(-0.5e0*v[1082] * v[1476] * v[337] + (-(v[1089] * v[1476]) + v[1478])*v[340] + (-
		(v[1095] * v[1476]) + v[1480])*v[342]);
	v[1585] = (*radius)*((-0.5e0*v[1082] * v[1471] - v[1534] - v[1535])*v[337] + (-(v[1089] * v[1471]) + v[1478]
		)*v[340] + (-(v[1095] * v[1471]) + v[1499])*v[342]);
	v[1584] = (*radius)*(-((0.5e0*v[1082] * v[1469] + v[1478])*v[337]) - (v[1089] * v[1469] + v[1508] - v[1518]
		)*v[340] - (v[1069] + v[1095] * v[1469] + v[1494])*v[342]);
	v[1583] = (*radius)*((-0.5e0*v[1082] * v[1465] - v[1530] - v[1532])*v[337] + (-(v[1089] * v[1465]) + v[1515]
		)*v[340] + (-(v[1095] * v[1465]) + v[1480])*v[342]);
	v[1582] = (*radius)*((-2e0*v[1474] - v[1464] * v[2710])*v[337] + (-(v[1089] * v[1464]) + v[2703])*v[340] -
		(v[1072] + v[1095] * v[1464] + v[1487])*v[342]);
	v[1581] = (*radius)*(-((0.5e0*v[1082] * v[1463] + v[1480])*v[337]) + (-(v[1089] * v[1463]) + v[2704])*v[340] -
		(v[1095] * v[1463] - v[1502] + v[1511])*v[342]);
	v[1580] = (*radius)*((-(v[1099] * v[1476]) + v[1478])*v[337] + (-v[1539] - v[1540] + v[1476] * v[2702])*v[340] + (-
		(v[1103] * v[1476]) + v[1489])*v[342]);
	v[1579] = (*radius)*((-(v[1099] * v[1471]) + v[1478])*v[337] + v[1471] * v[2702] * v[340] + (-(v[1103] * v[1471])
		+ v[1474])*v[342]);
	v[1578] = (*radius)*(-((v[1099] * v[1469] + v[1508] - v[1518])*v[337]) - (0.5e0*v[1078] * v[1469] + v[1478]
		)*v[340] + (-(v[1103] * v[1469]) + v[2703])*v[342]);
	v[1577] = (*radius)*(-((v[1099] * v[1465] + v[1515])*v[337]) - (0.5e0*v[1078] * v[1465] + v[1530] + v[1532]
		)*v[340] - (v[1103] * v[1465] - v[1474])*v[342]);
	v[1576] = (*radius)*(-((v[1071] + v[1099] * v[1464] + v[1491])*v[337]) - (0.5e0*v[1078] * v[1464] + v[1474]
		)*v[340] - (v[1103] * v[1464] - v[1502] + v[1508])*v[342]);
	v[1575] = (*radius)*(-((v[1069] + v[1099] * v[1463] + v[1494])*v[337]) + (-2e0*v[1480] - v[1463] * v[2709]
		)*v[340] + (-(v[1103] * v[1463]) + v[2706])*v[342]);
	v[1574] = (*radius)*(-((v[1107] * v[1476] - v[1480])*v[337]) - (v[1111] * v[1476] + v[1489])*v[340] -
		(0.5e0*v[1073] * v[1476] + v[1539] + v[1540])*v[342]);
	v[1573] = (*radius)*(-((v[1107] * v[1471] + v[1499])*v[337]) - (v[1111] * v[1471] - v[1474])*v[340] -
		(0.5e0*v[1073] * v[1471] + v[1534] + v[1535])*v[342]);
	v[1572] = (*radius)*((-(v[1107] * v[1469]) + v[2704])*v[337] - (v[1071] + v[1111] * v[1469] + v[1491])*v[340] + (
		-2e0*v[1478] + v[1469] * v[2705])*v[342]);
	v[1571] = (*radius)*((-(v[1107] * v[1465]) + v[1480])*v[337] + (-(v[1111] * v[1465]) + v[1474])*v[340]
		+ v[1465] * v[2705] * v[342]);
	v[1570] = (*radius)*((-(v[1107] * v[1464]) + v[2706])*v[337] - (v[1111] * v[1464] - v[1502] + v[1508])*v[340] -
		(0.5e0*v[1073] * v[1464] + v[1474])*v[342]);
	v[1569] = (*radius)*(-((v[1107] * v[1463] - v[1502] + v[1511])*v[337]) - (v[1072] + v[1111] * v[1463] + v[1487]
		)*v[340] - (v[1480] + v[1463] * v[2707])*v[342]);
	v[1168] = (*radius)*(v[1072] * v[2734] + (v[1091] - v[1072] * v[1107])*v[337] - (v[1088] + v[1072] * v[1111]
		)*v[340]);
	v[1167] = (*radius)*((-(v[1071] * v[1107]) + v[258])*v[337] + (v[1094] - v[1071] * v[1111])*v[340] -
		(0.5e0*v[1071] * v[1073] + v[1084])*v[342]);
	v[1166] = (*radius)*((v[1094] - v[1069] * v[1107])*v[337] - (v[1069] * v[1111] + v[258])*v[340] -
		(0.5e0*v[1069] * v[1073] + v[1077])*v[342]);
	v[1165] = (*radius)*(-((v[1072] * v[1099] + v[258])*v[337]) - (v[1086] + v[1072] * v[2708])*v[340] - (v[1088]
		+ v[1072] * v[1103])*v[342]);
	v[1164] = (*radius)*(v[1071] * v[2732] + (v[1091] - v[1071] * v[1099])*v[337] + (v[1094] - v[1071] * v[1103]
		)*v[342]);
	v[1163] = (*radius)*(-((v[1088] + v[1069] * v[1099])*v[337]) - 1e0*(v[1077] + v[1069] * v[2709])*v[340] + (-
		(v[1069] * v[1103]) + v[258])*v[342]);
	v[1162] = (*radius)*(-((v[1086] + v[1072] * v[2710])*v[337]) + (-(v[1072] * v[1089]) + v[258])*v[340] + (v[1091]
		- v[1072] * v[1095])*v[342]);
	v[1161] = (*radius)*(-((v[1084] + v[1071] * v[2711])*v[337]) - (v[1071] * v[1089] - v[1091])*v[340] -
		(v[1071] * v[1095] + v[258])*v[342]);
	v[1160] = (*radius)*(v[1069] * v[2729] - (v[1088] + v[1069] * v[1089])*v[340] + (v[1094] - v[1069] * v[1095]
		)*v[342]);
	v[1153] = v[337] * v[578] + v[340] * v[581] + v[342] * v[584];
	v[1152] = v[337] * v[577] + v[340] * v[580] + v[342] * v[583];
	v[1151] = v[337] * v[576] + v[340] * v[579] + v[342] * v[582];
	v[346] = v[323] * v[343];
	v[347] = v[323] * v[344];
	v[348] = v[323] * v[345];
	v[370] = -(v[346] * v[350]) - v[347] * v[352] - v[348] * v[354];
	v[369] = -(v[343] * v[346]) - v[344] * v[347] - v[345] * v[348];
	v[355] = v[333] * v[350];
	v[356] = v[333] * v[352];
	v[945] = v[338] * (-2e0*v[579] * v[585] + 2e0*v[576] * v[588] - v[355] * v[892] + v[356] * v[898] + v[346] * v[910]
		- v[347] * v[916]);
	v[946] = -(v[309] * v[766]) - v[308] * v[768] - v[307] * v[770] - (*radius)*v[945];
	v[942] = v[338] * (-2e0*v[580] * v[586] + 2e0*v[577] * v[589] - v[355] * v[891] + v[356] * v[897] + v[346] * v[909]
		- v[347] * v[915]);
	v[944] = -(v[309] * v[747]) - v[308] * v[750] - v[307] * v[753] - (*radius)*v[942];
	v[941] = v[338] * (-(v[580] * v[585]) - v[579] * v[586] + v[577] * v[588] + v[576] * v[589] - v[355] * v[890]
		+ v[356] * v[896] + v[346] * v[908] - v[347] * v[914]);
	v[943] = -(v[309] * v[746]) - v[308] * v[749] - v[307] * v[752] - (*radius)*v[941];
	v[937] = v[338] * (-2e0*v[581] * v[587] + 2e0*v[578] * v[590] - v[355] * v[889] + v[356] * v[895] + v[346] * v[907]
		- v[347] * v[913]);
	v[940] = -(v[309] * v[815]) - v[308] * v[816] - v[307] * v[817] - (*radius)*v[937];
	v[936] = v[338] * (-(v[581] * v[586]) - v[580] * v[587] + v[578] * v[589] + v[577] * v[590] - v[355] * v[888]
		+ v[356] * v[894] + v[346] * v[906] - v[347] * v[912]);
	v[939] = -(v[309] * v[748]) - v[308] * v[751] - v[307] * v[754] - (*radius)*v[936];
	v[935] = v[338] * (-(v[581] * v[585]) - v[579] * v[587] + v[578] * v[588] + v[576] * v[590] - v[355] * v[887]
		+ v[356] * v[893] + v[346] * v[905] - v[347] * v[911]);
	v[938] = -(v[309] * v[767]) - v[308] * v[769] - v[307] * v[771] - (*radius)*v[935];
	v[602] = v[338] * (v[356] * v[578] - v[355] * v[581] - v[347] * v[587] + v[346] * v[590]);
	v[611] = -v[575] - (*radius)*v[602];
	v[601] = v[338] * (v[356] * v[577] - v[355] * v[580] - v[347] * v[586] + v[346] * v[589]);
	v[610] = -v[574] - (*radius)*v[601];
	v[600] = v[338] * (v[356] * v[576] - v[355] * v[579] - v[347] * v[585] + v[346] * v[588]);
	v[2772] = -2e0*v[600];
	v[609] = -v[573] - (*radius)*v[600];
	v[357] = v[333] * v[354];
	v[969] = v[338] * (-2e0*v[582] * v[588] + 2e0*v[579] * v[591] - v[356] * v[886] + v[357] * v[892] + v[347] * v[904]
		- v[348] * v[910]);
	v[970] = -(v[309] * v[842]) - v[308] * v[845] - v[307] * v[848] - (*radius)*v[969];
	v[966] = v[338] * (-2e0*v[583] * v[589] + 2e0*v[580] * v[592] - v[356] * v[885] + v[357] * v[891] + v[347] * v[903]
		- v[348] * v[909]);
	v[968] = -(v[309] * v[875]) - v[308] * v[876] - v[307] * v[877] - (*radius)*v[966];
	v[965] = v[338] * (-(v[583] * v[588]) - v[582] * v[589] + v[580] * v[591] + v[579] * v[592] - v[356] * v[884]
		+ v[357] * v[890] + v[347] * v[902] - v[348] * v[908]);
	v[967] = -(v[309] * v[843]) - v[308] * v[846] - v[307] * v[849] - (*radius)*v[965];
	v[961] = v[338] * (-2e0*v[584] * v[590] + 2e0*v[581] * v[593] - v[356] * v[883] + v[357] * v[889] + v[347] * v[901]
		- v[348] * v[907]);
	v[964] = -(v[309] * v[863]) - v[308] * v[865] - v[307] * v[867] - (*radius)*v[961];
	v[960] = v[338] * (-(v[584] * v[589]) - v[583] * v[590] + v[581] * v[592] + v[580] * v[593] - v[356] * v[882]
		+ v[357] * v[888] + v[347] * v[900] - v[348] * v[906]);
	v[963] = -(v[309] * v[862]) - v[308] * v[864] - v[307] * v[866] - (*radius)*v[960];
	v[959] = v[338] * (-(v[584] * v[588]) - v[582] * v[590] + v[581] * v[591] + v[579] * v[593] - v[356] * v[881]
		+ v[357] * v[887] + v[347] * v[899] - v[348] * v[905]);
	v[962] = -(v[309] * v[844]) - v[308] * v[847] - v[307] * v[850] - (*radius)*v[959];
	v[957] = v[338] * (2e0*v[582] * v[585] - 2e0*v[576] * v[591] + v[355] * v[886] - v[357] * v[898] - v[346] * v[904]
		+ v[348] * v[916]);
	v[1604] = v[325] * v[945] + v[324] * v[957] + v[322] * v[969];
	v[1603] = v[342] * v[945] + v[340] * v[957] + v[337] * v[969];
	v[958] = -(v[309] * v[835]) - v[308] * v[836] - v[307] * v[837] - (*radius)*v[957];
	v[954] = v[338] * (2e0*v[583] * v[586] - 2e0*v[577] * v[592] + v[355] * v[885] - v[357] * v[897] - v[346] * v[903]
		+ v[348] * v[915]);
	v[1602] = v[325] * v[942] + v[324] * v[954] + v[322] * v[966];
	v[1600] = v[342] * v[942] + v[340] * v[954] + v[337] * v[966];
	v[956] = -(v[309] * v[793]) - v[308] * v[796] - v[307] * v[799] - (*radius)*v[954];
	v[953] = v[338] * (v[583] * v[585] + v[582] * v[586] - v[577] * v[591] - v[576] * v[592] + v[355] * v[884] - v[357] * v[896]
		- v[346] * v[902] + v[348] * v[914]);
	v[1601] = v[325] * v[941] + v[324] * v[953] + v[322] * v[965];
	v[1599] = v[342] * v[941] + v[340] * v[953] + v[337] * v[965];
	v[955] = -(v[309] * v[792]) - v[308] * v[795] - v[307] * v[798] - (*radius)*v[953];
	v[949] = v[338] * (2e0*v[584] * v[587] - 2e0*v[578] * v[593] + v[355] * v[883] - v[357] * v[895] - v[346] * v[901]
		+ v[348] * v[913]);
	v[1598] = v[325] * v[937] + v[324] * v[949] + v[322] * v[961];
	v[1595] = v[342] * v[937] + v[340] * v[949] + v[337] * v[961];
	v[952] = -(v[309] * v[823]) - v[308] * v[825] - v[307] * v[827] - (*radius)*v[949];
	v[948] = v[338] * (v[584] * v[586] + v[583] * v[587] - v[578] * v[592] - v[577] * v[593] + v[355] * v[882] - v[357] * v[894]
		- v[346] * v[900] + v[348] * v[912]);
	v[1597] = v[325] * v[936] + v[324] * v[948] + v[322] * v[960];
	v[1594] = v[342] * v[936] + v[340] * v[948] + v[337] * v[960];
	v[951] = -(v[309] * v[794]) - v[308] * v[797] - v[307] * v[800] - (*radius)*v[948];
	v[947] = v[338] * (v[584] * v[585] + v[582] * v[587] - v[578] * v[591] - v[576] * v[593] + v[355] * v[881] - v[357] * v[893]
		- v[346] * v[899] + v[348] * v[911]);
	v[1596] = v[325] * v[935] + v[324] * v[947] + v[322] * v[959];
	v[1593] = v[342] * v[935] + v[340] * v[947] + v[337] * v[959];
	v[950] = -(v[309] * v[822]) - v[308] * v[824] - v[307] * v[826] - (*radius)*v[947];
	v[599] = v[338] * (-(v[357] * v[578]) + v[355] * v[584] + v[348] * v[587] - v[346] * v[593]);
	v[608] = -v[572] - (*radius)*v[599];
	v[598] = v[338] * (-(v[357] * v[577]) + v[355] * v[583] + v[348] * v[586] - v[346] * v[592]);
	v[607] = -v[571] - (*radius)*v[598];
	v[597] = v[338] * (-(v[357] * v[576]) + v[355] * v[582] + v[348] * v[585] - v[346] * v[591]);
	v[2770] = -2e0*v[597];
	v[606] = -v[570] - (*radius)*v[597];
	v[596] = v[338] * (v[357] * v[581] - v[356] * v[584] - v[348] * v[590] + v[347] * v[593]);
	v[1676] = v[1160] * v[596] + v[1163] * v[599] + v[1166] * v[602];
	v[1657] = v[1161] * v[596] + v[1164] * v[599] + v[1167] * v[602];
	v[1634] = v[1162] * v[596] + v[1165] * v[599] + v[1168] * v[602];
	v[1159] = v[337] * v[596] + v[340] * v[599] + v[342] * v[602];
	v[1141] = v[322] * v[596] + v[324] * v[599] + v[325] * v[602];
	v[605] = -v[569] - (*radius)*v[596];
	v[595] = v[338] * (v[357] * v[580] - v[356] * v[583] - v[348] * v[589] + v[347] * v[592]);
	v[1675] = v[1160] * v[595] + v[1163] * v[598] + v[1166] * v[601];
	v[1656] = v[1161] * v[595] + v[1164] * v[598] + v[1167] * v[601];
	v[1633] = v[1162] * v[595] + v[1165] * v[598] + v[1168] * v[601];
	v[1158] = v[337] * v[595] + v[340] * v[598] + v[342] * v[601];
	v[1140] = v[322] * v[595] + v[324] * v[598] + v[325] * v[601];
	v[604] = -v[568] - (*radius)*v[595];
	v[594] = v[338] * (v[357] * v[579] - v[356] * v[582] - v[348] * v[588] + v[347] * v[591]);
	v[2768] = -2e0*v[594];
	v[1674] = v[1160] * v[594] + v[1163] * v[597] + v[1166] * v[600];
	v[1655] = v[1161] * v[594] + v[1164] * v[597] + v[1167] * v[600];
	v[1632] = v[1162] * v[594] + v[1165] * v[597] + v[1168] * v[600];
	v[1157] = v[337] * v[594] + v[340] * v[597] + v[342] * v[600];
	v[1139] = v[322] * v[594] + v[324] * v[597] + v[325] * v[600];
	v[603] = -v[567] - (*radius)*v[594];
	v[372] = -(v[350] * v[355]) - v[352] * v[356] - v[354] * v[357];
	v[371] = -(v[343] * v[355]) - v[344] * v[356] - v[345] * v[357];
	v[358] = v[338] * (-(v[348] * v[356]) + v[347] * v[357]);
	v[1184] = 1e0 - (v[358] * v[358]);
	v[359] = v[338] * (v[348] * v[355] - v[346] * v[357]);
	v[1619] = 2e0*v[579] * v[594] - 2e0*v[576] * v[597] + v[358] * v[892] - v[359] * v[898] - v[346] * v[957]
		+ v[347] * v[969];
	v[1618] = 2e0*v[580] * v[595] - 2e0*v[577] * v[598] + v[358] * v[891] - v[359] * v[897] - v[346] * v[954]
		+ v[347] * v[966];
	v[1617] = v[580] * v[594] + v[579] * v[595] - v[577] * v[597] - v[576] * v[598] + v[358] * v[890] - v[359] * v[896]
		- v[346] * v[953] + v[347] * v[965];
	v[1616] = 2e0*v[581] * v[596] - 2e0*v[578] * v[599] + v[358] * v[889] - v[359] * v[895] - v[346] * v[949]
		+ v[347] * v[961];
	v[1615] = v[581] * v[595] + v[580] * v[596] - v[578] * v[598] - v[577] * v[599] + v[358] * v[888] - v[359] * v[894]
		- v[346] * v[948] + v[347] * v[960];
	v[1614] = v[581] * v[594] + v[579] * v[596] - v[578] * v[597] - v[576] * v[599] + v[358] * v[887] - v[359] * v[893]
		- v[346] * v[947] + v[347] * v[959];
	v[1613] = -(v[359] * v[596]) - v[358] * v[599];
	v[1612] = -(v[359] * v[595]) - v[358] * v[598];
	v[1611] = -(v[359] * v[594]) - v[358] * v[597];
	v[1202] = 1e0 - (v[359] * v[359]);
	v[1185] = -(v[358] * v[359]);
	v[1132] = -(v[359] * v[578]) + v[358] * v[581] + v[347] * v[596] - v[346] * v[599];
	v[1131] = -(v[359] * v[577]) + v[358] * v[580] + v[347] * v[595] - v[346] * v[598];
	v[1130] = -(v[359] * v[576]) + v[358] * v[579] + v[347] * v[594] - v[346] * v[597];
	v[360] = v[338] * (-(v[347] * v[355]) + v[346] * v[356]);
	v[2746] = v[2736] * v[360];
	v[1710] = 2e0*v[582] * v[597] - 2e0*v[579] * v[600] + v[359] * v[886] - v[360] * v[892] - v[347] * v[945]
		+ v[348] * v[957];
	v[1705] = 2e0*v[583] * v[598] - 2e0*v[580] * v[601] + v[359] * v[885] - v[360] * v[891] - v[347] * v[942]
		+ v[348] * v[954];
	v[1704] = v[583] * v[597] + v[582] * v[598] - v[580] * v[600] - v[579] * v[601] + v[359] * v[884] - v[360] * v[890]
		- v[347] * v[941] + v[348] * v[953];
	v[1697] = 2e0*v[584] * v[599] - 2e0*v[581] * v[602] + v[359] * v[883] - v[360] * v[889] - v[347] * v[937]
		+ v[348] * v[949];
	v[1696] = v[584] * v[598] + v[583] * v[599] - v[581] * v[601] - v[580] * v[602] + v[359] * v[882] - v[360] * v[888]
		- v[347] * v[936] + v[348] * v[948];
	v[1695] = v[584] * v[597] + v[582] * v[599] - v[581] * v[600] - v[579] * v[602] + v[359] * v[881] - v[360] * v[887]
		- v[347] * v[935] + v[348] * v[947];
	v[1694] = -2e0*v[582] * v[594] + 2e0*v[576] * v[600] - v[358] * v[886] + v[360] * v[898] + v[346] * v[945]
		- v[348] * v[969];
	v[1712] = v[1710] * v[322] + v[1694] * v[324] + v[1619] * v[325];
	v[1711] = v[1710] * v[337] + v[1694] * v[340] + v[1619] * v[342];
	v[1693] = -2e0*v[583] * v[595] + 2e0*v[577] * v[601] - v[358] * v[885] + v[360] * v[897] + v[346] * v[942]
		- v[348] * v[966];
	v[1709] = v[1705] * v[322] + v[1693] * v[324] + v[1618] * v[325];
	v[1707] = v[1705] * v[337] + v[1693] * v[340] + v[1618] * v[342];
	v[1692] = -(v[583] * v[594]) - v[582] * v[595] + v[577] * v[600] + v[576] * v[601] - v[358] * v[884] + v[360] * v[896]
		+ v[346] * v[941] - v[348] * v[965];
	v[1708] = v[1704] * v[322] + v[1692] * v[324] + v[1617] * v[325];
	v[1706] = v[1704] * v[337] + v[1692] * v[340] + v[1617] * v[342];
	v[1691] = -2e0*v[584] * v[596] + 2e0*v[578] * v[602] - v[358] * v[883] + v[360] * v[895] + v[346] * v[937]
		- v[348] * v[961];
	v[1703] = v[1697] * v[322] + v[1691] * v[324] + v[1616] * v[325];
	v[1700] = v[1697] * v[337] + v[1691] * v[340] + v[1616] * v[342];
	v[1690] = -(v[584] * v[595]) - v[583] * v[596] + v[578] * v[601] + v[577] * v[602] - v[358] * v[882] + v[360] * v[894]
		+ v[346] * v[936] - v[348] * v[960];
	v[1702] = v[1696] * v[322] + v[1690] * v[324] + v[1615] * v[325];
	v[1699] = v[1696] * v[337] + v[1690] * v[340] + v[1615] * v[342];
	v[1689] = -(v[584] * v[594]) - v[582] * v[596] + v[578] * v[600] + v[576] * v[602] - v[358] * v[881] + v[360] * v[893]
		+ v[346] * v[935] - v[348] * v[959];
	v[1701] = v[1695] * v[322] + v[1689] * v[324] + v[1614] * v[325];
	v[1698] = v[1695] * v[337] + v[1689] * v[340] + v[1614] * v[342];
	v[1673] = v[1586] * v[358] + v[1580] * v[359] + v[1574] * v[360];
	v[1654] = v[1585] * v[358] + v[1579] * v[359] + v[1573] * v[360];
	v[1653] = v[1584] * v[358] + v[1578] * v[359] + v[1572] * v[360];
	v[1631] = v[1583] * v[358] + v[1577] * v[359] + v[1571] * v[360];
	v[1630] = v[1582] * v[358] + v[1576] * v[359] + v[1570] * v[360];
	v[1629] = v[1581] * v[358] + v[1575] * v[359] + v[1569] * v[360];
	v[1628] = -(v[360] * v[596]) - v[358] * v[602];
	v[1627] = -(v[360] * v[595]) - v[358] * v[601];
	v[1626] = -(v[360] * v[594]) - v[358] * v[600];
	v[1625] = -(v[360] * v[599]) - v[359] * v[602];
	v[1624] = -(v[360] * v[598]) - v[359] * v[601];
	v[1623] = -(v[360] * v[597]) - v[359] * v[600];
	v[1219] = 1e0 - (v[360] * v[360]);
	v[1203] = -(v[359] * v[360]);
	v[1186] = -(v[358] * v[360]);
	v[1171] = v[1162] * v[358] + v[1165] * v[359] + v[1168] * v[360];
	v[1222] = v[1168] - v[1171] * v[360];
	v[1206] = v[1165] - v[1171] * v[359];
	v[1189] = v[1162] - v[1171] * v[358];
	v[1170] = v[1161] * v[358] + v[1164] * v[359] + v[1167] * v[360];
	v[1221] = v[1167] - v[1170] * v[360];
	v[1205] = v[1164] - v[1170] * v[359];
	v[1188] = v[1161] - v[1170] * v[358];
	v[1169] = v[1160] * v[358] + v[1163] * v[359] + v[1166] * v[360];
	v[1220] = v[1166] - v[1169] * v[360];
	v[1204] = v[1163] - v[1169] * v[359];
	v[1187] = v[1160] - v[1169] * v[358];
	v[1129] = v[360] * v[578] - v[358] * v[584] - v[348] * v[596] + v[346] * v[602];
	v[1128] = v[360] * v[577] - v[358] * v[583] - v[348] * v[595] + v[346] * v[601];
	v[1127] = v[360] * v[576] - v[358] * v[582] - v[348] * v[594] + v[346] * v[600];
	v[1126] = -(v[360] * v[581]) + v[359] * v[584] + v[348] * v[599] - v[347] * v[602];
	v[1156] = v[1126] * v[337] + v[1129] * v[340] + v[1132] * v[342];
	v[1138] = v[1126] * v[322] + v[1129] * v[324] + v[1132] * v[325];
	v[1125] = -(v[360] * v[580]) + v[359] * v[583] + v[348] * v[598] - v[347] * v[601];
	v[1155] = v[1125] * v[337] + v[1128] * v[340] + v[1131] * v[342];
	v[1137] = v[1125] * v[322] + v[1128] * v[324] + v[1131] * v[325];
	v[1124] = -(v[360] * v[579]) + v[359] * v[582] + v[348] * v[597] - v[347] * v[600];
	v[1154] = v[1124] * v[337] + v[1127] * v[340] + v[1130] * v[342];
	v[1136] = v[1124] * v[322] + v[1127] * v[324] + v[1130] * v[325];
	v[361] = -((*radius)*v[358]) + v[365];
	v[362] = -((*radius)*v[359]) + v[366];
	v[363] = -((*radius)*v[360]) + v[367];
	v[1059] = v[2736] * sqrt((v[361] * v[361]) + (v[362] * v[362]) + (v[363] * v[363]));
	v[373] = -(v[309] * v[348]);
	v[374] = -(v[309] * v[357]) + v[387];
	v[375] = -(v[308] * v[348]) + v[385];
	v[376] = -(v[308] * v[357]);
	v[377] = -(v[309] * v[347]);
	v[378] = -(v[309] * v[356]) + v[397];
	v[379] = -(v[308] * v[347]) + v[395];
	v[380] = -(v[308] * v[356]);
	v[381] = -(v[309] * v[346]);
	v[382] = -(v[309] * v[355]) + v[407];
	v[383] = -(v[308] * v[346]) + v[405];
	v[384] = -(v[308] * v[355]);
	v[386] = -(v[307] * v[348]) - v[385];
	v[388] = -(v[307] * v[357]) - v[387];
	v[389] = dCi[2] * v[373] + dBi[2] * v[375] + dAi[2] * v[386];
	v[390] = dCi[2] * v[374] + dBi[2] * v[376] + dAi[2] * v[388];
	v[391] = dCi[1] * v[373] + dBi[1] * v[375] + dAi[1] * v[386];
	v[1992] = v[1988] * v[239] + v[391] * v[479];
	v[1991] = v[1987] * v[239] + v[391] * v[478];
	v[444] = v[239] * v[391];
	v[392] = dCi[1] * v[374] + dBi[1] * v[376] + dAi[1] * v[388];
	v[1999] = v[1995] * v[239] + v[392] * v[479];
	v[1998] = v[1994] * v[239] + v[392] * v[478];
	v[457] = v[239] * v[392];
	v[393] = dCi[0] * v[373] + dBi[0] * v[375] + dAi[0] * v[386];
	v[2006] = v[2002] * v[239] + v[393] * v[479];
	v[2005] = v[2001] * v[239] + v[393] * v[478];
	v[447] = v[239] * v[393];
	v[394] = dCi[0] * v[374] + dBi[0] * v[376] + dAi[0] * v[388];
	v[2013] = v[2009] * v[239] + v[394] * v[479];
	v[2012] = v[2008] * v[239] + v[394] * v[478];
	v[460] = v[239] * v[394];
	v[396] = -(v[307] * v[347]) - v[395];
	v[398] = -(v[307] * v[356]) - v[397];
	v[399] = dCi[2] * v[377] + dBi[2] * v[379] + dAi[2] * v[396];
	v[2724] = v[391] - v[399];
	v[2713] = v[391] + v[399];
	v[2026] = v[2022] * v[239] + v[399] * v[479];
	v[2114] = v[1992] + v[2026];
	v[2025] = v[2021] * v[239] + v[399] * v[478];
	v[445] = v[239] * v[399];
	v[400] = dCi[2] * v[378] + dBi[2] * v[380] + dAi[2] * v[398];
	v[2720] = v[392] - v[400];
	v[2715] = v[392] + v[400];
	v[2033] = v[2029] * v[239] + v[400] * v[479];
	v[2117] = v[1999] + v[2033];
	v[2032] = v[2028] * v[239] + v[400] * v[478];
	v[458] = v[239] * v[400];
	v[401] = dCi[1] * v[377] + dBi[1] * v[379] + dAi[1] * v[396];
	v[2725] = -v[389] - v[401];
	v[402] = dCi[1] * v[378] + dBi[1] * v[380] + dAi[1] * v[398];
	v[2721] = -v[390] - v[402];
	v[403] = dCi[0] * v[377] + dBi[0] * v[379] + dAi[0] * v[396];
	v[2045] = v[2042] * v[239] + v[403] * v[479];
	v[452] = v[239] * v[403];
	v[404] = dCi[0] * v[378] + dBi[0] * v[380] + dAi[0] * v[398];
	v[2051] = v[2048] * v[239] + v[404] * v[479];
	v[465] = v[239] * v[404];
	v[406] = -(v[307] * v[346]) - v[405];
	v[408] = -(v[307] * v[355]) - v[407];
	v[409] = dCi[2] * v[381] + dBi[2] * v[383] + dAi[2] * v[406];
	v[2723] = v[393] + v[409];
	v[2063] = v[2060] * v[239] + v[409] * v[479];
	v[2120] = v[2006] + v[2063];
	v[2062] = v[2059] * v[239] + v[409] * v[478];
	v[448] = v[239] * v[409];
	v[410] = dCi[2] * v[382] + dBi[2] * v[384] + dAi[2] * v[408];
	v[2719] = v[394] + v[410];
	v[2069] = v[2066] * v[239] + v[410] * v[479];
	v[2123] = v[2013] + v[2069];
	v[2068] = v[2065] * v[239] + v[410] * v[478];
	v[461] = v[239] * v[410];
	v[411] = dCi[1] * v[381] + dBi[1] * v[383] + dAi[1] * v[406];
	v[2712] = v[403] + v[411];
	v[2075] = v[2072] * v[239] + v[411] * v[479];
	v[2126] = v[2045] + v[2075];
	v[2125] = (v[2041] + v[2071])*v[239] + v[2712] * v[478];
	v[453] = v[239] * v[411];
	v[412] = dCi[1] * v[382] + dBi[1] * v[384] + dAi[1] * v[408];
	v[2714] = v[404] + v[412];
	v[2081] = v[2078] * v[239] + v[412] * v[479];
	v[2129] = v[2051] + v[2081];
	v[2128] = (v[2047] + v[2077])*v[239] + v[2714] * v[478];
	v[466] = v[239] * v[412];
	v[413] = dCi[0] * v[381] + dBi[0] * v[383] + dAi[0] * v[406];
	v[414] = dCi[0] * v[382] + dBi[0] * v[384] + dAi[0] * v[408];
	v[2111] = v[1917] * v[389] + v[1982] * v[742];
	v[2107] = v[1917] * v[401] + v[2036] * v[742];
	v[2099] = v[1917] * v[390] + v[1985] * v[742];
	v[2095] = v[1917] * v[402] + v[2039] * v[742];
	v[468] = v[414] * v[742];
	v[467] = v[402] * v[742];
	v[455] = v[413] * v[742];
	v[454] = v[401] * v[742];
	v[421] = v[444] + v[445];
	v[422] = v[457] + v[458];
	v[427] = v[447] + v[448];
	v[428] = v[460] + v[461];
	v[430] = v[452] + v[453];
	v[431] = v[465] + v[466];
	v[434] = v[2699] * v[389] + v[2700] * v[401] + v[2701] * v[413] + v[391] * v[418] + v[393] * v[419] + v[399] * v[420]
		+ v[403] * v[425] + v[409] * v[426] + v[411] * v[429];
	v[435] = v[2699] * v[390] + v[2700] * v[402] + v[2701] * v[414] + v[392] * v[418] + v[394] * v[419] + v[400] * v[420]
		+ v[404] * v[425] + v[410] * v[426] + v[412] * v[429];
	v[437] = -4e0*v[476];
	v[2156] = v[2144] * v[434] + v[437] * (v[1982] * v[2699] + v[2036] * v[2700] + v[2084] * v[2701] + v[403] - v[411] + d[11] *
		(-v[401] - v[413]) + v[1988] * v[418] + v[2002] * v[419] + v[2022] * v[420] + v[2042] * v[425] + v[2060] * v[426]
		+ v[2072] * v[429] + v[2713] * v[471] + v[2723] * v[472]);
	v[2717] = v[2156] + v[1917] * v[413] + v[2084] * v[742];
	v[2722] = -(v[1904] * v[389]) + v[2143] * v[434] + v[437] * (v[1981] * v[2699] + v[2035] * v[2700] + v[2083] * v[2701]
		- v[393] + v[409] + d[10] * (-v[389] - v[413]) + v[1987] * v[418] + v[2001] * v[419] + v[2021] * v[420] + v[2041] * v[425]
		+ v[2059] * v[426] + v[2071] * v[429] + v[2712] * v[472] + v[2713] * v[474]) + v[1981] * v[742];
	v[2150] = v[2144] * v[435] + v[437] * (v[1985] * v[2699] + v[2039] * v[2700] + v[2087] * v[2701] + v[404] - v[412] + d[11] *
		(-v[402] - v[414]) + v[1995] * v[418] + v[2009] * v[419] + v[2029] * v[420] + v[2048] * v[425] + v[2066] * v[426]
		+ v[2078] * v[429] + v[2715] * v[471] + v[2719] * v[472]);
	v[2716] = v[2150] + v[1917] * v[414] + v[2087] * v[742];
	v[2718] = -(v[1904] * v[390]) + v[2143] * v[435] + v[437] * (v[1984] * v[2699] + v[2038] * v[2700] + v[2086] * v[2701]
		- v[394] + v[410] + d[10] * (-v[390] - v[414]) + v[1994] * v[418] + v[2008] * v[419] + v[2028] * v[420] + v[2047] * v[425]
		+ v[2065] * v[426] + v[2077] * v[429] + v[2714] * v[472] + v[2715] * v[474]) + v[1984] * v[742];
	v[2173] = v[435] * v[437] + v[467] + v[468];
	v[2170] = v[2173] - v[467] + v[390] * v[742];
	v[2166] = v[2170] + v[467] - v[468];
	v[2164] = v[434] * v[437] + v[454] + v[455];
	v[2161] = v[2164] - v[454] + v[389] * v[742];
	v[2157] = v[2161] + v[454] - v[455];
	v[2174] = v[2051] - v[2081] + 2e0*v[2173] + v[2117] * v[471] + v[2123] * v[472] + (v[2095] + v[2716])*v[475];
	v[2165] = v[2045] - v[2075] + 2e0*v[2164] + v[2114] * v[471] + v[2120] * v[472] + (v[2107] + v[2717])*v[475];
	v[2172] = -v[2013] + v[2069] + 0.5e0*v[422] + v[2129] * v[472] + (v[2099] + v[2716])*v[473] + v[2117] * v[474];
	v[2171] = -v[2012] + v[2068] + 2e0*v[2170] + v[2128] * v[472] + (v[1998] + v[2032])*v[474] + v[473] * (v[2718]
		- v[1904] * v[414] + v[2086] * v[742]);
	v[2163] = -v[2006] + v[2063] + 0.5e0*v[421] + v[2126] * v[472] + (v[2111] + v[2717])*v[473] + v[2114] * v[474];
	v[2162] = -v[2005] + v[2062] + 2e0*v[2161] + v[2125] * v[472] + (v[1991] + v[2025])*v[474] + v[473] * (v[2722]
		- v[1904] * v[413] + v[2083] * v[742]);
	v[2169] = v[1999] - v[2033] + 0.5e0*v[428] + (v[2095] + v[2099] + v[2150])*v[470] + v[2129] * v[471]
		+ v[2123] * v[474];
	v[2168] = v[1998] - v[2032] + 0.5e0*v[431] + v[2128] * v[471] + (v[2012] + v[2068])*v[474] + v[470] * (v[2718]
		- v[1904] * v[402] + v[2038] * v[742]);
	v[2167] = 2e0*v[2166] + (v[1993] - v[2027])*v[239] + v[2720] * v[477] + v[471] * ((v[2046] + v[2076])*v[239]
		+ v[2714] * v[477]) + v[474] * ((v[2007] + v[2064])*v[239] + v[2719] * v[477]) + v[470] * (v[1910] * v[2721]
			+ v[2142] * v[435] + v[437] * (v[1983] * v[2699] + v[2037] * v[2700] + (dCi[0] * v[1965] + dBi[0] * v[1971]
				+ dAi[0] * v[2055])*v[2701] + v[2720] + d[9] * v[2721] + v[1993] * v[418] + v[2007] * v[419] + v[2027] * v[420]
				+ v[2046] * v[425] + v[2064] * v[426] + v[2076] * v[429] + v[2714] * v[471] + v[2719] * v[474]) + (v[1983] + v[2037]
					)*v[742]);
	v[2160] = v[1992] - v[2026] + 0.5e0*v[427] + (v[2107] + v[2111] + v[2156])*v[470] + v[2126] * v[471]
		+ v[2120] * v[474];
	v[2159] = v[1991] - v[2025] + 0.5e0*v[430] + v[2125] * v[471] + (v[2005] + v[2062])*v[474] + v[470] * (v[2722]
		- v[1904] * v[401] + v[2035] * v[742]);
	v[2158] = 2e0*v[2157] + (v[1986] - v[2020])*v[239] + v[2724] * v[477] + v[471] * ((v[2040] + v[2070])*v[239]
		+ v[2712] * v[477]) + v[474] * ((v[2000] + v[2058])*v[239] + v[2723] * v[477]) + v[470] * (v[1910] * v[2725]
			+ v[2142] * v[434] + v[437] * (v[1980] * v[2699] + v[2034] * v[2700] + (dCi[0] * v[1962] + dBi[0] * v[1968]
				+ dAi[0] * v[2052])*v[2701] + v[2724] + d[9] * v[2725] + v[1986] * v[418] + v[2000] * v[419] + v[2020] * v[420]
				+ v[2040] * v[425] + v[2058] * v[426] + v[2070] * v[429] + v[2712] * v[471] + v[2723] * v[474]) + (v[1980] + v[2034]
					)*v[742]);
	v[446] = v[444] - v[445] + v[2157] * v[470] + v[430] * v[471] + v[427] * v[474];
	v[451] = -v[447] + v[448] + v[430] * v[472] + v[2161] * v[473] + v[421] * v[474];
	v[456] = v[452] - v[453] + v[421] * v[471] + v[427] * v[472] + v[2164] * v[475];
	v[459] = v[457] - v[458] + v[2166] * v[470] + v[431] * v[471] + v[428] * v[474];
	v[464] = -v[460] + v[461] + v[431] * v[472] + v[2170] * v[473] + v[422] * v[474];
	v[469] = v[465] - v[466] + v[422] * v[471] + v[428] * v[472] + v[2173] * v[475];
	v[612] = 1e0 / (-(v[370] * v[371]) + v[369] * v[372]);
	v[2285] = (-(v[2174] * v[369]) + v[2165] * v[371])*v[612];
	v[2251] = (-(v[2172] * v[369]) + v[2163] * v[371])*v[612];
	v[2250] = (-(v[2171] * v[369]) + v[2162] * v[371])*v[612];
	v[2219] = (-(v[2169] * v[369]) + v[2160] * v[371])*v[612];
	v[2218] = (-(v[2168] * v[369]) + v[2159] * v[371])*v[612];
	v[2217] = (-(v[2167] * v[369]) + v[2158] * v[371])*v[612];
	v[2213] = (v[371] * v[584] - v[369] * v[593])*v[612];
	v[2212] = (v[371] * v[583] - v[369] * v[592])*v[612];
	v[2211] = (v[371] * v[582] - v[369] * v[591])*v[612];
	v[2204] = (v[371] * v[581] - v[369] * v[590])*v[612];
	v[2203] = (v[371] * v[580] - v[369] * v[589])*v[612];
	v[2202] = (v[371] * v[579] - v[369] * v[588])*v[612];
	v[2192] = (v[371] * v[578] - v[369] * v[587])*v[612];
	v[2191] = (v[371] * v[577] - v[369] * v[586])*v[612];
	v[2190] = (v[371] * v[576] - v[369] * v[585])*v[612];
	v[2189] = (v[2174] * v[370] - v[2165] * v[372])*v[612];
	v[2188] = (v[2172] * v[370] - v[2163] * v[372])*v[612];
	v[2319] = -(v[2188] * v[343]) - v[2251] * v[350] + v[963];
	v[2309] = -(v[2188] * v[344]) - v[2251] * v[352] + v[951];
	v[2297] = -(v[2188] * v[345]) - v[2251] * v[354] + v[939];
	v[2187] = (v[2171] * v[370] - v[2162] * v[372])*v[612];
	v[2186] = (v[2169] * v[370] - v[2160] * v[372])*v[612];
	v[2317] = -(v[2186] * v[343]) - v[2219] * v[350] + v[962];
	v[2307] = -(v[2186] * v[344]) - v[2219] * v[352] + v[950];
	v[2295] = -(v[2186] * v[345]) - v[2219] * v[354] + v[938];
	v[2185] = (v[2168] * v[370] - v[2159] * v[372])*v[612];
	v[2281] = -(v[2185] * v[343]) - v[2218] * v[350] + v[967];
	v[2272] = -(v[2185] * v[344]) - v[2218] * v[352] + v[955];
	v[2261] = -(v[2185] * v[345]) - v[2218] * v[354] + v[943];
	v[2184] = (v[2167] * v[370] - v[2158] * v[372])*v[612];
	v[2183] = (-(v[372] * v[584]) + v[370] * v[593])*v[612];
	v[2316] = -(v[2183] * v[343]) - v[2213] * v[350];
	v[2306] = -(v[2183] * v[344]) - v[2213] * v[352];
	v[2294] = -(v[2183] * v[345]) - v[2213] * v[354];
	v[2182] = (-(v[372] * v[583]) + v[370] * v[592])*v[612];
	v[2280] = -(v[2182] * v[343]) - v[2212] * v[350];
	v[2271] = -(v[2182] * v[344]) - v[2212] * v[352];
	v[2260] = -(v[2182] * v[345]) - v[2212] * v[354];
	v[2181] = (-(v[372] * v[582]) + v[370] * v[591])*v[612];
	v[2246] = -(v[2181] * v[343]) - v[2211] * v[350];
	v[2238] = -(v[2181] * v[344]) - v[2211] * v[352];
	v[2228] = -(v[2181] * v[345]) - v[2211] * v[354];
	v[2180] = (-(v[372] * v[581]) + v[370] * v[590])*v[612];
	v[2315] = -(v[2180] * v[343]) - v[2204] * v[350];
	v[2305] = -(v[2180] * v[344]) - v[2204] * v[352];
	v[2291] = -(v[2180] * v[345]) - v[2204] * v[354];
	v[2179] = (-(v[372] * v[580]) + v[370] * v[589])*v[612];
	v[2279] = -(v[2179] * v[343]) - v[2203] * v[350];
	v[2270] = -(v[2179] * v[344]) - v[2203] * v[352];
	v[2257] = -(v[2179] * v[345]) - v[2203] * v[354];
	v[2178] = (-(v[372] * v[579]) + v[370] * v[588])*v[612];
	v[2245] = -(v[2178] * v[343]) - v[2202] * v[350];
	v[2237] = -(v[2178] * v[344]) - v[2202] * v[352];
	v[2225] = -(v[2178] * v[345]) - v[2202] * v[354];
	v[2177] = (-(v[372] * v[578]) + v[370] * v[587])*v[612];
	v[2314] = -(v[2177] * v[343]) - v[2192] * v[350];
	v[2302] = -(v[2177] * v[344]) - v[2192] * v[352];
	v[2288] = -(v[2177] * v[345]) - v[2192] * v[354];
	v[2176] = (-(v[372] * v[577]) + v[370] * v[586])*v[612];
	v[2278] = -(v[2176] * v[343]) - v[2191] * v[350];
	v[2267] = -(v[2176] * v[344]) - v[2191] * v[352];
	v[2254] = -(v[2176] * v[345]) - v[2191] * v[354];
	v[2175] = (-(v[372] * v[576]) + v[370] * v[585])*v[612];
	v[2244] = -(v[2175] * v[343]) - v[2190] * v[350];
	v[2234] = -(v[2175] * v[344]) - v[2190] * v[352];
	v[2222] = -(v[2175] * v[345]) - v[2190] * v[354];
	v[613] = (v[355] * v[370] - v[346] * v[372])*v[612];
	v[614] = (v[356] * v[370] - v[347] * v[372])*v[612];
	v[615] = (v[357] * v[370] - v[348] * v[372])*v[612];
	v[616] = (-(v[372] * v[446]) + v[370] * v[459])*v[612];
	v[617] = (-(v[372] * v[451]) + v[370] * v[464])*v[612];
	v[618] = (-(v[372] * v[456]) + v[370] * v[469])*v[612];
	v[619] = (-(v[355] * v[369]) + v[346] * v[371])*v[612];
	v[2201] = v[2302] - v[554] * v[613] - v[563] * v[619];
	v[2200] = v[2267] - v[553] * v[613] - v[562] * v[619];
	v[2199] = v[2234] - v[552] * v[613] - v[561] * v[619];
	v[2198] = v[2288] - v[557] * v[613] - v[566] * v[619];
	v[2197] = v[2254] - v[556] * v[613] - v[565] * v[619];
	v[2196] = v[2222] - v[555] * v[613] - v[564] * v[619];
	v[2195] = v[2314] - v[551] * v[613] - v[560] * v[619];
	v[2194] = v[2278] - v[550] * v[613] - v[559] * v[619];
	v[2193] = v[2244] - v[549] * v[613] - v[558] * v[619];
	v[971] = 1e0 - v[343] * v[613] - v[350] * v[619];
	v[2779] = (*a4)*v[971];
	v[628] = -(v[345] * v[613]) - v[354] * v[619];
	v[2402] = (*a4)*v[628];
	v[626] = -(v[344] * v[613]) - v[352] * v[619];
	v[2395] = (*a4)*v[626];
	v[620] = (-(v[356] * v[369]) + v[347] * v[371])*v[612];
	v[2210] = v[2291] - v[557] * v[614] - v[566] * v[620];
	v[2209] = v[2257] - v[556] * v[614] - v[565] * v[620];
	v[2208] = v[2225] - v[555] * v[614] - v[564] * v[620];
	v[2207] = v[2305] - v[554] * v[614] - v[563] * v[620];
	v[2206] = v[2270] - v[553] * v[614] - v[562] * v[620];
	v[2205] = v[2237] - v[552] * v[614] - v[561] * v[620];
	v[981] = 1e0 - v[344] * v[614] - v[352] * v[620];
	v[2781] = (*a4)*v[981];
	v[629] = -(v[345] * v[614]) - v[354] * v[620];
	v[2404] = (*a4)*v[629];
	v[621] = (-(v[357] * v[369]) + v[348] * v[371])*v[612];
	v[2216] = v[2294] - v[557] * v[615] - v[566] * v[621];
	v[2215] = v[2260] - v[556] * v[615] - v[565] * v[621];
	v[2214] = v[2228] - v[555] * v[615] - v[564] * v[621];
	v[990] = 1e0 - v[345] * v[615] - v[354] * v[621];
	v[2784] = (*a4)*v[990];
	v[622] = (v[371] * v[446] - v[369] * v[459])*v[612];
	v[2249] = v[2317] - v[551] * v[616] - v[560] * v[622];
	v[2248] = v[2281] - v[550] * v[616] - v[559] * v[622];
	v[2247] = -(v[2184] * v[343]) - v[2217] * v[350] - v[549] * v[616] - v[558] * v[622] + v[970];
	v[2241] = v[2307] - v[554] * v[616] - v[563] * v[622];
	v[2240] = v[2272] - v[553] * v[616] - v[562] * v[622];
	v[2239] = -(v[2184] * v[344]) - v[2217] * v[352] - v[552] * v[616] - v[561] * v[622] + v[958];
	v[2231] = v[2295] - v[557] * v[616] - v[566] * v[622];
	v[2230] = v[2261] - v[556] * v[616] - v[565] * v[622];
	v[2229] = -(v[2184] * v[345]) - v[2217] * v[354] - v[555] * v[616] - v[564] * v[622] + v[946];
	v[992] = v[609] - v[345] * v[616] - v[354] * v[622];
	v[984] = v[606] - v[344] * v[616] - v[352] * v[622];
	v[975] = v[603] - v[343] * v[616] - v[350] * v[622];
	v[2741] = v[2736] * (v[358] * v[975] + v[359] * v[984] + v[360] * v[992]);
	v[623] = (v[371] * v[451] - v[369] * v[464])*v[612];
	v[2284] = v[2319] - v[551] * v[617] - v[560] * v[623];
	v[2283] = -(v[2187] * v[343]) - v[2250] * v[350] - v[550] * v[617] - v[559] * v[623] + v[968];
	v[2282] = v[2281] - v[549] * v[617] - v[558] * v[623];
	v[2275] = v[2309] - v[554] * v[617] - v[563] * v[623];
	v[2274] = -(v[2187] * v[344]) - v[2250] * v[352] - v[553] * v[617] - v[562] * v[623] + v[956];
	v[2273] = v[2272] - v[552] * v[617] - v[561] * v[623];
	v[2264] = v[2297] - v[557] * v[617] - v[566] * v[623];
	v[2263] = -(v[2187] * v[345]) - v[2250] * v[354] - v[556] * v[617] - v[565] * v[623] + v[944];
	v[2262] = v[2261] - v[555] * v[617] - v[564] * v[623];
	v[994] = v[610] - v[345] * v[617] - v[354] * v[623];
	v[986] = v[607] - v[344] * v[617] - v[352] * v[623];
	v[977] = v[604] - v[343] * v[617] - v[350] * v[623];
	v[2739] = v[2736] * (v[358] * v[977] + v[359] * v[986] + v[360] * v[994]);
	v[624] = (v[371] * v[456] - v[369] * v[469])*v[612];
	v[2321] = -(v[2189] * v[343]) - v[2285] * v[350] - v[551] * v[618] - v[560] * v[624] + v[964];
	v[2320] = v[2319] - v[550] * v[618] - v[559] * v[624];
	v[2318] = v[2317] - v[549] * v[618] - v[558] * v[624];
	v[2311] = -(v[2189] * v[344]) - v[2285] * v[352] - v[554] * v[618] - v[563] * v[624] + v[952];
	v[2310] = v[2309] - v[553] * v[618] - v[562] * v[624];
	v[2308] = v[2307] - v[552] * v[618] - v[561] * v[624];
	v[2299] = -(v[2189] * v[345]) - v[2285] * v[354] - v[557] * v[618] - v[566] * v[624] + v[940];
	v[2298] = v[2297] - v[556] * v[618] - v[565] * v[624];
	v[2296] = v[2295] - v[555] * v[618] - v[564] * v[624];
	v[996] = v[611] - v[345] * v[618] - v[354] * v[624];
	v[988] = v[608] - v[344] * v[618] - v[352] * v[624];
	v[979] = v[605] - v[343] * v[618] - v[350] * v[624];
	v[2737] = v[2736] * (v[358] * v[979] + v[359] * v[988] + v[360] * v[996]);
	v[640] = v[309] * v[345];
	v[641] = v[309] * v[354] + v[646];
	v[642] = v[309] * v[344];
	v[643] = v[309] * v[352] + v[649];
	v[644] = v[309] * v[343];
	v[645] = v[309] * v[350] + v[652];
	v[647] = v[308] * v[345] + v[646];
	v[648] = v[308] * v[354];
	v[650] = v[308] * v[344] + v[649];
	v[651] = v[308] * v[352];
	v[653] = v[308] * v[343] + v[652];
	v[654] = v[308] * v[350];
	v[655] = v[307] * v[345] - v[646];
	v[656] = v[307] * v[354] - v[646];
	v[657] = dCi[2] * v[640] + dBi[2] * v[647] + dAi[2] * v[655];
	v[658] = dCi[2] * v[641] + dBi[2] * v[648] + dAi[2] * v[656];
	v[659] = dCi[1] * v[640] + dBi[1] * v[647] + dAi[1] * v[655];
	v[660] = dCi[1] * v[641] + dBi[1] * v[648] + dAi[1] * v[656];
	v[661] = dCi[0] * v[640] + dBi[0] * v[647] + dAi[0] * v[655];
	v[999] = v[239] * v[661];
	v[662] = dCi[0] * v[641] + dBi[0] * v[648] + dAi[0] * v[656];
	v[1001] = v[239] * v[662];
	v[663] = v[307] * v[344] - v[649];
	v[664] = v[307] * v[352] - v[649];
	v[665] = dCi[2] * v[642] + dBi[2] * v[650] + dAi[2] * v[663];
	v[666] = dCi[2] * v[643] + dBi[2] * v[651] + dAi[2] * v[664];
	v[667] = dCi[1] * v[642] + dBi[1] * v[650] + dAi[1] * v[663];
	v[668] = dCi[1] * v[643] + dBi[1] * v[651] + dAi[1] * v[664];
	v[669] = dCi[0] * v[642] + dBi[0] * v[650] + dAi[0] * v[663];
	v[1004] = v[239] * v[669];
	v[670] = dCi[0] * v[643] + dBi[0] * v[651] + dAi[0] * v[664];
	v[1008] = v[239] * v[670];
	v[671] = v[307] * v[343] - v[652];
	v[672] = v[307] * v[350] - v[652];
	v[673] = dCi[2] * v[644] + dBi[2] * v[653] + dAi[2] * v[671];
	v[1000] = v[239] * v[673];
	v[674] = dCi[2] * v[645] + dBi[2] * v[654] + dAi[2] * v[672];
	v[1002] = v[239] * v[674];
	v[675] = dCi[1] * v[644] + dBi[1] * v[653] + dAi[1] * v[671];
	v[1005] = v[239] * v[675];
	v[676] = dCi[1] * v[645] + dBi[1] * v[654] + dAi[1] * v[672];
	v[1009] = v[239] * v[676];
	v[677] = dCi[0] * v[644] + dBi[0] * v[653] + dAi[0] * v[671];
	v[678] = dCi[0] * v[645] + dBi[0] * v[654] + dAi[0] * v[672];
	v[679] = v[239] * (v[659] + v[665]);
	v[680] = v[239] * (v[660] + v[666]);
	v[2726] = v[437] * (v[2699] * v[657] + v[418] * v[659] + v[419] * v[661] + v[420] * v[665] + v[2700] * v[667]
		+ v[425] * v[669] + v[426] * v[673] + v[429] * v[675] + v[2701] * v[677]) + v[677] * v[742];
	v[1016] = v[1004] - v[1005] + v[471] * v[679] + v[475] * (v[2726] + v[667] * v[742]) + v[472] * (v[1000] + v[999]);
	v[1013] = v[1000] + (v[1004] + v[1005])*v[472] + v[474] * v[679] + v[473] * (v[2726] + v[657] * v[742]) - v[999];
	v[2727] = v[437] * (v[2699] * v[658] + v[418] * v[660] + v[419] * v[662] + v[420] * v[666] + v[2700] * v[668]
		+ v[425] * v[670] + v[426] * v[674] + v[429] * v[676] + v[2701] * v[678]) + v[678] * v[742];
	v[1017] = v[1008] - v[1009] + (v[1001] + v[1002])*v[472] + v[471] * v[680] + v[475] * (v[2727] + v[668] * v[742]);
	v[1014] = -v[1001] + v[1002] + (v[1008] + v[1009])*v[472] + v[474] * v[680] + v[473] * (v[2727] + v[658] * v[742]);
	v[1020] = v[325] * v[340] - v[324] * v[342];
	v[1021] = -(v[325] * v[337]) + v[322] * v[342];
	v[1022] = v[324] * v[337] - v[322] * v[340];
	v[1781] = v[1022] * v[886] + v[1021] * v[892] + v[1020] * v[898];
	v[1784] = v[1020] * v[1781] + v[1568] * v[322] + v[1592] * v[337];
	v[1783] = v[1021] * v[1781] + v[1568] * v[324] + v[1592] * v[340];
	v[1782] = v[1022] * v[1781] + v[1568] * v[325] + v[1592] * v[342];
	v[1774] = v[1022] * v[885] + v[1021] * v[891] + v[1020] * v[897];
	v[1780] = v[1020] * v[1774] + v[1567] * v[322] + v[1591] * v[337];
	v[1778] = v[1021] * v[1774] + v[1567] * v[324] + v[1591] * v[340];
	v[1776] = v[1022] * v[1774] + v[1567] * v[325] + v[1591] * v[342];
	v[1773] = v[1022] * v[884] + v[1021] * v[890] + v[1020] * v[896];
	v[1779] = v[1020] * v[1773] + v[1566] * v[322] + v[1590] * v[337];
	v[1777] = v[1021] * v[1773] + v[1566] * v[324] + v[1590] * v[340];
	v[1775] = v[1022] * v[1773] + v[1566] * v[325] + v[1590] * v[342];
	v[1763] = v[1022] * v[883] + v[1021] * v[889] + v[1020] * v[895];
	v[1772] = v[1020] * v[1763] + v[1565] * v[322] + v[1589] * v[337];
	v[1769] = v[1021] * v[1763] + v[1565] * v[324] + v[1589] * v[340];
	v[1766] = v[1022] * v[1763] + v[1565] * v[325] + v[1589] * v[342];
	v[1762] = v[1022] * v[882] + v[1021] * v[888] + v[1020] * v[894];
	v[1771] = v[1020] * v[1762] + v[1564] * v[322] + v[1588] * v[337];
	v[1768] = v[1021] * v[1762] + v[1564] * v[324] + v[1588] * v[340];
	v[1765] = v[1022] * v[1762] + v[1564] * v[325] + v[1588] * v[342];
	v[1761] = v[1022] * v[881] + v[1021] * v[887] + v[1020] * v[893];
	v[1770] = v[1020] * v[1761] + v[1563] * v[322] + v[1587] * v[337];
	v[1767] = v[1021] * v[1761] + v[1563] * v[324] + v[1587] * v[340];
	v[1764] = v[1022] * v[1761] + v[1563] * v[325] + v[1587] * v[342];
	v[1757] = v[1022] * v[1619] + v[1021] * v[1694] + v[1020] * v[1710];
	v[1760] = v[1020] * v[1757] + v[1712] * v[322] + v[1711] * v[337];
	v[1759] = v[1021] * v[1757] + v[1712] * v[324] + v[1711] * v[340];
	v[1758] = v[1022] * v[1757] + v[1712] * v[325] + v[1711] * v[342];
	v[1750] = v[1022] * v[1618] + v[1021] * v[1693] + v[1020] * v[1705];
	v[1756] = v[1020] * v[1750] + v[1709] * v[322] + v[1707] * v[337];
	v[1754] = v[1021] * v[1750] + v[1709] * v[324] + v[1707] * v[340];
	v[1752] = v[1022] * v[1750] + v[1709] * v[325] + v[1707] * v[342];
	v[1749] = v[1022] * v[1617] + v[1021] * v[1692] + v[1020] * v[1704];
	v[1755] = v[1020] * v[1749] + v[1708] * v[322] + v[1706] * v[337];
	v[1753] = v[1021] * v[1749] + v[1708] * v[324] + v[1706] * v[340];
	v[1751] = v[1022] * v[1749] + v[1708] * v[325] + v[1706] * v[342];
	v[1739] = v[1022] * v[1616] + v[1021] * v[1691] + v[1020] * v[1697];
	v[1748] = v[1020] * v[1739] + v[1703] * v[322] + v[1700] * v[337];
	v[1745] = v[1021] * v[1739] + v[1703] * v[324] + v[1700] * v[340];
	v[1742] = v[1022] * v[1739] + v[1703] * v[325] + v[1700] * v[342];
	v[1738] = v[1022] * v[1615] + v[1021] * v[1690] + v[1020] * v[1696];
	v[1747] = v[1020] * v[1738] + v[1702] * v[322] + v[1699] * v[337];
	v[1744] = v[1021] * v[1738] + v[1702] * v[324] + v[1699] * v[340];
	v[1741] = v[1022] * v[1738] + v[1702] * v[325] + v[1699] * v[342];
	v[1737] = v[1022] * v[1614] + v[1021] * v[1689] + v[1020] * v[1695];
	v[1746] = v[1020] * v[1737] + v[1701] * v[322] + v[1698] * v[337];
	v[1743] = v[1021] * v[1737] + v[1701] * v[324] + v[1698] * v[340];
	v[1740] = v[1022] * v[1737] + v[1701] * v[325] + v[1698] * v[342];
	v[1733] = v[1022] * v[945] + v[1021] * v[957] + v[1020] * v[969];
	v[1736] = v[1020] * v[1733] + v[1604] * v[322] + v[1603] * v[337];
	v[1735] = v[1021] * v[1733] + v[1604] * v[324] + v[1603] * v[340];
	v[1734] = v[1022] * v[1733] + v[1604] * v[325] + v[1603] * v[342];
	v[1726] = v[1022] * v[942] + v[1021] * v[954] + v[1020] * v[966];
	v[1732] = v[1020] * v[1726] + v[1602] * v[322] + v[1600] * v[337];
	v[1730] = v[1021] * v[1726] + v[1602] * v[324] + v[1600] * v[340];
	v[1728] = v[1022] * v[1726] + v[1602] * v[325] + v[1600] * v[342];
	v[1725] = v[1022] * v[941] + v[1021] * v[953] + v[1020] * v[965];
	v[1731] = v[1020] * v[1725] + v[1601] * v[322] + v[1599] * v[337];
	v[1729] = v[1021] * v[1725] + v[1601] * v[324] + v[1599] * v[340];
	v[1727] = v[1022] * v[1725] + v[1601] * v[325] + v[1599] * v[342];
	v[1715] = v[1022] * v[937] + v[1021] * v[949] + v[1020] * v[961];
	v[1724] = v[1020] * v[1715] + v[1598] * v[322] + v[1595] * v[337];
	v[1721] = v[1021] * v[1715] + v[1598] * v[324] + v[1595] * v[340];
	v[1718] = v[1022] * v[1715] + v[1598] * v[325] + v[1595] * v[342];
	v[1714] = v[1022] * v[936] + v[1021] * v[948] + v[1020] * v[960];
	v[1723] = v[1020] * v[1714] + v[1597] * v[322] + v[1594] * v[337];
	v[1720] = v[1021] * v[1714] + v[1597] * v[324] + v[1594] * v[340];
	v[1717] = v[1022] * v[1714] + v[1597] * v[325] + v[1594] * v[342];
	v[1713] = v[1022] * v[935] + v[1021] * v[947] + v[1020] * v[959];
	v[1722] = v[1020] * v[1713] + v[1596] * v[322] + v[1593] * v[337];
	v[1719] = v[1021] * v[1713] + v[1596] * v[324] + v[1593] * v[340];
	v[1716] = v[1022] * v[1713] + v[1596] * v[325] + v[1593] * v[342];
	v[1150] = v[1020] * v[596] + v[1021] * v[599] + v[1022] * v[602];
	v[1218] = v[1022] * v[1150] + v[1141] * v[325] + v[1159] * v[342];
	v[1201] = v[1021] * v[1150] + v[1141] * v[324] + v[1159] * v[340];
	v[1183] = v[1020] * v[1150] + v[1141] * v[322] + v[1159] * v[337];
	v[1149] = v[1020] * v[595] + v[1021] * v[598] + v[1022] * v[601];
	v[1217] = v[1022] * v[1149] + v[1140] * v[325] + v[1158] * v[342];
	v[1200] = v[1021] * v[1149] + v[1140] * v[324] + v[1158] * v[340];
	v[1182] = v[1020] * v[1149] + v[1140] * v[322] + v[1158] * v[337];
	v[1148] = v[1020] * v[594] + v[1021] * v[597] + v[1022] * v[600];
	v[1216] = v[1022] * v[1148] + v[1139] * v[325] + v[1157] * v[342];
	v[1199] = v[1021] * v[1148] + v[1139] * v[324] + v[1157] * v[340];
	v[1181] = v[1020] * v[1148] + v[1139] * v[322] + v[1157] * v[337];
	v[1147] = v[1020] * v[1126] + v[1021] * v[1129] + v[1022] * v[1132];
	v[1215] = v[1022] * v[1147] + v[1138] * v[325] + v[1156] * v[342];
	v[1198] = v[1021] * v[1147] + v[1138] * v[324] + v[1156] * v[340];
	v[1180] = v[1020] * v[1147] + v[1138] * v[322] + v[1156] * v[337];
	v[1146] = v[1020] * v[1125] + v[1021] * v[1128] + v[1022] * v[1131];
	v[1214] = v[1022] * v[1146] + v[1137] * v[325] + v[1155] * v[342];
	v[1197] = v[1021] * v[1146] + v[1137] * v[324] + v[1155] * v[340];
	v[1179] = v[1020] * v[1146] + v[1137] * v[322] + v[1155] * v[337];
	v[1145] = v[1020] * v[1124] + v[1021] * v[1127] + v[1022] * v[1130];
	v[1213] = v[1022] * v[1145] + v[1136] * v[325] + v[1154] * v[342];
	v[1196] = v[1021] * v[1145] + v[1136] * v[324] + v[1154] * v[340];
	v[1178] = v[1020] * v[1145] + v[1136] * v[322] + v[1154] * v[337];
	v[1144] = v[1020] * v[578] + v[1021] * v[581] + v[1022] * v[584];
	v[1212] = v[1022] * v[1144] + v[1135] * v[325] + v[1153] * v[342];
	v[1195] = v[1021] * v[1144] + v[1135] * v[324] + v[1153] * v[340];
	v[1177] = v[1020] * v[1144] + v[1135] * v[322] + v[1153] * v[337];
	v[1143] = v[1020] * v[577] + v[1021] * v[580] + v[1022] * v[583];
	v[1211] = v[1022] * v[1143] + v[1134] * v[325] + v[1152] * v[342];
	v[1194] = v[1021] * v[1143] + v[1134] * v[324] + v[1152] * v[340];
	v[1176] = v[1020] * v[1143] + v[1134] * v[322] + v[1152] * v[337];
	v[1142] = v[1020] * v[576] + v[1021] * v[579] + v[1022] * v[582];
	v[1210] = v[1022] * v[1142] + v[1133] * v[325] + v[1151] * v[342];
	v[1193] = v[1021] * v[1142] + v[1133] * v[324] + v[1151] * v[340];
	v[1175] = v[1020] * v[1142] + v[1133] * v[322] + v[1151] * v[337];
	v[1023] = v[348] * v[359] - v[347] * v[360];
	v[1024] = -(v[348] * v[358]) + v[346] * v[360];
	v[1025] = v[347] * v[358] - v[346] * v[359];
	v[1026] = v[322] * v[346] + v[324] * v[347] + v[325] * v[348];
	v[1027] = v[1023] * v[322] + v[1024] * v[324] + v[1025] * v[325];
	v[1028] = v[322] * v[358] + v[324] * v[359] + v[325] * v[360];
	v[1029] = v[1020] * v[346] + v[1021] * v[347] + v[1022] * v[348];
	v[1030] = v[1020] * v[1023] + v[1021] * v[1024] + v[1022] * v[1025];
	v[1031] = v[1020] * v[358] + v[1021] * v[359] + v[1022] * v[360];
	v[1032] = v[337] * v[346] + v[340] * v[347] + v[342] * v[348];
	v[1033] = v[1023] * v[337] + v[1024] * v[340] + v[1025] * v[342];
	v[1034] = v[337] * v[358] + v[340] * v[359] + v[342] * v[360];
	v[1036] = d[0] + (v[286] - v[295])*v[304] + (v[289] - v[298])*v[305] + (v[292] - v[301])*v[306] + v[2730] * (v[2729]
		- v[1089] * v[340] - v[1095] * v[342]);
	v[1038] = d[1] + (v[287] - v[296])*v[304] + (v[290] - v[299])*v[305] + (v[293] - v[302])*v[306] + v[2730] * (v[2732]
		- v[1099] * v[337] - v[1103] * v[342]);
	v[1040] = d[2] + (v[288] - v[297])*v[304] + (v[291] - v[300])*v[305] + (v[294] - v[303])*v[306] + v[2730] * (v[2734]
		- v[1107] * v[337] - v[1111] * v[340]);
	v[1790] = v[1115] * v[2768] + v[1118] * v[2770] + v[1121] * v[2772] - v[1562] * v[358] - v[1556] * v[359]
		- v[1550] * v[360] + v[1040] * v[945] + v[1038] * v[957] + v[1036] * v[969];
	v[1789] = -(v[1561] * v[358]) - v[1555] * v[359] - v[1549] * v[360] - 2e0*v[1116] * v[595] - 2e0*v[1119] * v[598]
		- 2e0*v[1122] * v[601] + v[1040] * v[942] + v[1038] * v[954] + v[1036] * v[966];
	v[1788] = -(v[1560] * v[358]) - v[1554] * v[359] - v[1548] * v[360] - v[1116] * v[594] - v[1115] * v[595]
		- v[1119] * v[597] - v[1118] * v[598] - v[1122] * v[600] - v[1121] * v[601] + v[1040] * v[941] + v[1038] * v[953]
		+ v[1036] * v[965];
	v[1787] = -(v[1559] * v[358]) - v[1553] * v[359] - v[1547] * v[360] - 2e0*v[1117] * v[596] - 2e0*v[1120] * v[599]
		- 2e0*v[1123] * v[602] + v[1040] * v[937] + v[1038] * v[949] + v[1036] * v[961];
	v[1786] = -(v[1558] * v[358]) - v[1552] * v[359] - v[1546] * v[360] - v[1117] * v[595] - v[1116] * v[596]
		- v[1120] * v[598] - v[1119] * v[599] - v[1123] * v[601] - v[1122] * v[602] + v[1040] * v[936] + v[1038] * v[948]
		+ v[1036] * v[960];
	v[1785] = -(v[1557] * v[358]) - v[1551] * v[359] - v[1545] * v[360] - v[1117] * v[594] - v[1115] * v[596]
		- v[1120] * v[597] - v[1118] * v[599] - v[1123] * v[600] - v[1121] * v[602] + v[1040] * v[935] + v[1038] * v[947]
		+ v[1036] * v[959];
	v[1174] = -(v[1117] * v[358]) - v[1120] * v[359] - v[1123] * v[360] + v[1036] * v[596] + v[1038] * v[599]
		+ v[1040] * v[602];
	v[2799] = -2e0*v[1174];
	v[1173] = -(v[1116] * v[358]) - v[1119] * v[359] - v[1122] * v[360] + v[1036] * v[595] + v[1038] * v[598]
		+ v[1040] * v[601];
	v[2797] = -2e0*v[1173];
	v[1172] = -(v[1115] * v[358]) - v[1118] * v[359] - v[1121] * v[360] + v[1036] * v[594] + v[1038] * v[597]
		+ v[1040] * v[600];
	v[2794] = -2e0*v[1172];
	v[1048] = v[1036] * v[358] + v[1038] * v[359] + v[1040] * v[360];
	v[1225] = -v[1123] + gti[0] * (v[1020] * v[1215] + v[1212] * v[322] + v[1218] * v[337]) + gti[1] * (v[1021] * v[1215]
		+ v[1212] * v[324] + v[1218] * v[340]) + gti[2] * (v[1022] * v[1215] + v[1212] * v[325] + v[1218] * v[342])
		- v[1174] * v[360] - v[1048] * v[602];
	v[1224] = -v[1122] + gti[0] * (v[1020] * v[1214] + v[1211] * v[322] + v[1217] * v[337]) + gti[1] * (v[1021] * v[1214]
		+ v[1211] * v[324] + v[1217] * v[340]) + gti[2] * (v[1022] * v[1214] + v[1211] * v[325] + v[1217] * v[342])
		- v[1173] * v[360] - v[1048] * v[601];
	v[1223] = -v[1121] + gti[0] * (v[1020] * v[1213] + v[1210] * v[322] + v[1216] * v[337]) + gti[1] * (v[1021] * v[1213]
		+ v[1210] * v[324] + v[1216] * v[340]) + gti[2] * (v[1022] * v[1213] + v[1210] * v[325] + v[1216] * v[342])
		- v[1172] * v[360] - v[1048] * v[600];
	v[1209] = -v[1120] + gti[0] * (v[1020] * v[1198] + v[1195] * v[322] + v[1201] * v[337]) + gti[1] * (v[1021] * v[1198]
		+ v[1195] * v[324] + v[1201] * v[340]) + gti[2] * (v[1022] * v[1198] + v[1195] * v[325] + v[1201] * v[342])
		- v[1174] * v[359] - v[1048] * v[599];
	v[1208] = -v[1119] + gti[0] * (v[1020] * v[1197] + v[1194] * v[322] + v[1200] * v[337]) + gti[1] * (v[1021] * v[1197]
		+ v[1194] * v[324] + v[1200] * v[340]) + gti[2] * (v[1022] * v[1197] + v[1194] * v[325] + v[1200] * v[342])
		- v[1173] * v[359] - v[1048] * v[598];
	v[1207] = -v[1118] + gti[0] * (v[1020] * v[1196] + v[1193] * v[322] + v[1199] * v[337]) + gti[1] * (v[1021] * v[1196]
		+ v[1193] * v[324] + v[1199] * v[340]) + gti[2] * (v[1022] * v[1196] + v[1193] * v[325] + v[1199] * v[342])
		- v[1172] * v[359] - v[1048] * v[597];
	v[1192] = -v[1117] + gti[0] * (v[1020] * v[1180] + v[1177] * v[322] + v[1183] * v[337]) + gti[1] * (v[1021] * v[1180]
		+ v[1177] * v[324] + v[1183] * v[340]) + gti[2] * (v[1022] * v[1180] + v[1177] * v[325] + v[1183] * v[342])
		- v[1174] * v[358] - v[1048] * v[596];
	v[1191] = -v[1116] + gti[0] * (v[1020] * v[1179] + v[1176] * v[322] + v[1182] * v[337]) + gti[1] * (v[1021] * v[1179]
		+ v[1176] * v[324] + v[1182] * v[340]) + gti[2] * (v[1022] * v[1179] + v[1176] * v[325] + v[1182] * v[342])
		- v[1173] * v[358] - v[1048] * v[595];
	v[1190] = -v[1115] + gti[0] * (v[1020] * v[1178] + v[1175] * v[322] + v[1181] * v[337]) + gti[1] * (v[1021] * v[1178]
		+ v[1175] * v[324] + v[1181] * v[340]) + gti[2] * (v[1022] * v[1178] + v[1175] * v[325] + v[1181] * v[342])
		- v[1172] * v[358] - v[1048] * v[594];
	v[1041] = v[1020] * v[1029] + v[1026] * v[322] + v[1032] * v[337];
	v[1042] = v[1020] * v[1030] + v[1027] * v[322] + v[1033] * v[337];
	v[1043] = v[1020] * v[1031] + v[1028] * v[322] + v[1034] * v[337];
	v[1044] = v[1036] + gti[0] * (v[1020] * v[1042] + v[1041] * v[322] + v[1043] * v[337]) + gti[1] * (v[1021] * v[1042]
		+ v[1041] * v[324] + v[1043] * v[340]) + gti[2] * (v[1022] * v[1042] + v[1041] * v[325] + v[1043] * v[342])
		- v[1048] * v[358];
	v[1045] = v[1021] * v[1029] + v[1026] * v[324] + v[1032] * v[340];
	v[1046] = v[1021] * v[1030] + v[1027] * v[324] + v[1033] * v[340];
	v[1047] = v[1021] * v[1031] + v[1028] * v[324] + v[1034] * v[340];
	v[1049] = v[1038] + gti[0] * (v[1020] * v[1046] + v[1045] * v[322] + v[1047] * v[337]) + gti[1] * (v[1021] * v[1046]
		+ v[1045] * v[324] + v[1047] * v[340]) + gti[2] * (v[1022] * v[1046] + v[1045] * v[325] + v[1047] * v[342])
		- v[1048] * v[359];
	v[1050] = v[1022] * v[1029] + v[1026] * v[325] + v[1032] * v[342];
	v[1051] = v[1022] * v[1030] + v[1027] * v[325] + v[1033] * v[342];
	v[1052] = v[1022] * v[1031] + v[1028] * v[325] + v[1034] * v[342];
	v[1053] = v[1040] + gti[0] * (v[1020] * v[1051] + v[1050] * v[322] + v[1052] * v[337]) + gti[1] * (v[1021] * v[1051]
		+ v[1050] * v[324] + v[1052] * v[340]) + gti[2] * (v[1022] * v[1051] + v[1050] * v[325] + v[1052] * v[342])
		- v[1048] * v[360];
	v[1249] = 1e0 / sqrt((v[1044] * v[1044]) + (v[1049] * v[1049]) + (v[1053] * v[1053]));
	v[2735] = v[1059] * v[1249];
	v[1453] = v[1225] * v[2735];
	v[1452] = v[1209] * v[2735];
	v[1451] = v[1192] * v[2735];
	v[1441] = v[1224] * v[2735];
	v[1440] = v[1208] * v[2735];
	v[1439] = v[1191] * v[2735];
	v[1429] = v[1223] * v[2735];
	v[1428] = v[1207] * v[2735];
	v[1427] = v[1190] * v[2735];
	v[1411] = v[1221] * v[2735];
	v[1410] = v[1205] * v[2735];
	v[1409] = v[1188] * v[2735];
	v[1399] = v[1220] * v[2735];
	v[1398] = v[1204] * v[2735];
	v[1397] = v[1187] * v[2735];
	v[1387] = v[1219] * v[2735];
	v[1364] = v[1203] * v[2735];
	v[1362] = v[1202] * v[2735];
	v[1339] = v[1186] * v[2735];
	v[1338] = v[1185] * v[2735];
	v[1337] = v[1184] * v[2735];
	v[1054] = 1e0*v[1044] * v[1249];
	v[2769] = v[1054] * v[358];
	v[2767] = -(v[1054] * v[1184]);
	v[2764] = -(v[1054] * v[1185]);
	v[2762] = -(v[1054] * v[1187]);
	v[2760] = -(v[1054] * v[1188]);
	v[2758] = -(v[1054] * v[1189]);
	v[2756] = -(v[1054] * v[1190]);
	v[2754] = -(v[1054] * v[1191]);
	v[2752] = -(v[1054] * v[1192]);
	v[2751] = v[1054] * v[2736];
	v[1328] = 1e0 - (v[1054] * v[1054]);
	v[1056] = 1e0*v[1049] * v[1249];
	v[2771] = v[1056] * v[359];
	v[2765] = v[1056] * v[1202];
	v[2763] = v[1056] * v[1204];
	v[2761] = v[1056] * v[1205];
	v[2759] = v[1056] * v[1206];
	v[2757] = v[1056] * v[1207];
	v[2755] = v[1056] * v[1208];
	v[2753] = v[1056] * v[1209];
	v[2749] = v[1056] * v[2736];
	v[1330] = 1e0 - (v[1056] * v[1056]);
	v[1329] = v[1056] * v[2764];
	v[1057] = 1e0*v[1053] * v[1249];
	v[2773] = v[1057] * v[360];
	v[2766] = v[1057] * v[2736];
	v[2750] = -(v[1057] * v[1186]);
	v[2748] = v[1057] * v[1203];
	v[2747] = -(v[1057] * v[1219]);
	v[2745] = -(v[1057] * v[1220]);
	v[2744] = -(v[1057] * v[1221]);
	v[2743] = -(v[1057] * v[1222]);
	v[2742] = -(v[1057] * v[1223]);
	v[2740] = -(v[1057] * v[1224]);
	v[2738] = -(v[1057] * v[1225]);
	v[1380] = v[1056] * v[2737] + v[2735] * (v[1209] * v[1330] + v[1056] * (v[2738] + v[2752]));
	v[1379] = v[1054] * v[2737] - v[2735] * (-(v[1192] * v[1328]) + v[1054] * (-v[2738] + v[2753]));
	v[1376] = v[1056] * v[2739] + v[2735] * (v[1208] * v[1330] + v[1056] * (v[2740] + v[2754]));
	v[1375] = v[1054] * v[2739] - v[2735] * (-(v[1191] * v[1328]) + v[1054] * (-v[2740] + v[2755]));
	v[1372] = v[1056] * v[2741] + v[2735] * (v[1207] * v[1330] + v[1056] * (v[2742] + v[2756]));
	v[1371] = v[1054] * v[2741] - v[2735] * (-(v[1190] * v[1328]) + v[1054] * (-v[2742] + v[2757]));
	v[1368] = v[1206] * v[1330] + v[1056] * (v[2743] + v[2758]);
	v[1367] = v[1189] * v[1328] - v[1054] * (-v[2743] + v[2759]);
	v[1363] = v[1205] * v[1330] + v[1056] * (v[2744] + v[2760]);
	v[1361] = v[1188] * v[1328] - v[1054] * (-v[2744] + v[2761]);
	v[1358] = v[1204] * v[1330] + v[1056] * (v[2745] + v[2762]);
	v[1357] = v[1187] * v[1328] - v[1054] * (-v[2745] + v[2763]);
	v[1354] = v[1056] * v[2746] + v[2735] * (v[1203] * v[1330] + v[1056] * (-(v[1054] * v[1186]) + v[2747]));
	v[1353] = v[1054] * v[2746] - v[2735] * (-(v[1186] * v[1328]) + v[1054] * (v[1056] * v[1203] - v[2747]));
	v[1349] = -(v[2735] * (-(v[1185] * v[1328]) + v[1054] * (v[2748] + v[2765]))) + v[2751] * v[359];
	v[1346] = v[2735] * (v[1185] * v[1330] + v[1056] * (v[2750] + v[2767])) + v[2749] * v[358];
	v[1334] = -(v[1056] * v[2748]);
	v[1350] = (v[1329] + v[1202] * v[1330] + v[1334])*v[2735] + v[2749] * v[359];
	v[1333] = v[1054] * v[2750];
	v[1345] = (v[1184] * v[1328] + v[1329] + v[1333])*v[2735] + v[2751] * v[358];
	v[1331] = 1e0 - (v[1057] * v[1057]);
	v[1381] = v[1057] * v[2737] + v[2735] * (v[1225] * v[1331] + v[1057] * (v[2752] - v[2753]));
	v[1377] = v[1057] * v[2739] + v[2735] * (v[1224] * v[1331] + v[1057] * (v[2754] - v[2755]));
	v[1373] = v[1057] * v[2741] + v[2735] * (v[1223] * v[1331] + v[1057] * (v[2756] - v[2757]));
	v[1369] = v[1222] * v[1331] + v[1057] * (v[2758] - v[2759]);
	v[1365] = v[1221] * v[1331] + v[1057] * (v[2760] - v[2761]);
	v[1359] = v[1220] * v[1331] + v[1057] * (v[2762] - v[2763]);
	v[1355] = (v[1219] * v[1331] + v[1333] + v[1334])*v[2735] + v[1057] * v[2746];
	v[1351] = v[2735] * (v[1203] * v[1331] + v[1057] * (v[2764] - v[2765])) + v[2766] * v[359];
	v[1347] = v[2735] * (v[1186] * v[1331] + v[1057] * (-(v[1056] * v[1185]) + v[2767])) + v[2766] * v[358];
	v[1058] = v[1054] * v[1059];
	v[1060] = v[1056] * v[1059];
	v[1061] = v[1057] * v[1059];
	v[1336] = v[1337] * v[1357] + v[1338] * v[1358] + v[1339] * v[1359];
	v[1340] = v[1337] * v[1361] + v[1338] * v[1363] + v[1339] * v[1365];
	v[1341] = v[1337] * v[1367] + v[1338] * v[1368] + v[1339] * v[1369];
	v[1360] = v[1338] * v[1357] + v[1358] * v[1362] + v[1359] * v[1364];
	v[1366] = v[1338] * v[1361] + v[1362] * v[1363] + v[1364] * v[1365];
	v[1370] = v[1338] * v[1367] + v[1362] * v[1368] + v[1364] * v[1369];
	v[1386] = v[1339] * v[1357] + v[1358] * v[1364] + v[1359] * v[1387];
	v[1388] = v[1339] * v[1361] + v[1363] * v[1364] + v[1365] * v[1387];
	v[1389] = v[1339] * v[1367] + v[1364] * v[1368] + v[1369] * v[1387];
	v[1393] = v[1187] * v[1345] + v[1204] * v[1346] + v[1220] * v[1347];
	v[1394] = v[1187] * v[1349] + v[1204] * v[1350] + v[1220] * v[1351];
	v[1395] = v[1187] * v[1353] + v[1204] * v[1354] + v[1220] * v[1355];
	v[1405] = v[1188] * v[1345] + v[1205] * v[1346] + v[1221] * v[1347];
	v[1406] = v[1188] * v[1349] + v[1205] * v[1350] + v[1221] * v[1351];
	v[1407] = v[1188] * v[1353] + v[1205] * v[1354] + v[1221] * v[1355];
	v[1416] = v[1189] * v[1345] + v[1206] * v[1346] + v[1222] * v[1347];
	v[1417] = v[1189] * v[1349] + v[1206] * v[1350] + v[1222] * v[1351];
	v[1418] = v[1189] * v[1353] + v[1206] * v[1354] + v[1222] * v[1355];
	v[2785] = v[1059] * (v[1056] * v[1611] + v[1057] * v[1626] + v[2768] * v[2769]) + (*epsn)*v[975];
	v[2788] = v[1059] * (v[1056] * v[1612] + v[1057] * v[1627] - 2e0*v[2769] * v[595]) + (*epsn)*v[977];
	v[2791] = v[1059] * (v[1056] * v[1613] + v[1057] * v[1628] - 2e0*v[2769] * v[596]) + (*epsn)*v[979];
	v[2786] = v[1059] * (v[1054] * v[1611] + v[1057] * v[1623] + v[2770] * v[2771]) + (*epsn)*v[984];
	v[2789] = v[1059] * (v[1054] * v[1612] + v[1057] * v[1624] - 2e0*v[2771] * v[598]) + (*epsn)*v[986];
	v[2792] = v[1059] * (v[1054] * v[1613] + v[1057] * v[1625] - 2e0*v[2771] * v[599]) + (*epsn)*v[988];
	v[2787] = v[1059] * (v[1056] * v[1623] + v[1054] * v[1626] + v[2772] * v[2773]) + (*epsn)*v[992];
	v[2790] = v[1059] * (v[1056] * v[1624] + v[1054] * v[1627] - 2e0*v[2773] * v[601]) + (*epsn)*v[994];
	v[2793] = v[1059] * (v[1056] * v[1625] + v[1054] * v[1628] - 2e0*v[2773] * v[602]) + (*epsn)*v[996];
	v[1821] = v[1059] * (-(v[1054] * (v[1674] * v[358] + v[1169] * v[594])) - v[1056] * (v[1674] * v[359] + v[1169] * v[597])
		- v[1057] * (v[1674] * v[360] + v[1169] * v[600]));
	v[1822] = v[1059] * (-(v[1054] * (v[1675] * v[358] + v[1169] * v[595])) - v[1056] * (v[1675] * v[359] + v[1169] * v[598])
		- v[1057] * (v[1675] * v[360] + v[1169] * v[601]));
	v[1823] = v[1059] * (-(v[1054] * (v[1676] * v[358] + v[1169] * v[596])) - v[1056] * (v[1676] * v[359] + v[1169] * v[599])
		- v[1057] * (v[1676] * v[360] + v[1169] * v[602]));
	v[1826] = v[1059] * (-(v[1054] * (v[1655] * v[358] + v[1170] * v[594])) - v[1056] * (v[1655] * v[359] + v[1170] * v[597])
		- v[1057] * (v[1655] * v[360] + v[1170] * v[600]));
	v[1827] = v[1059] * (-(v[1054] * (v[1656] * v[358] + v[1170] * v[595])) - v[1056] * (v[1656] * v[359] + v[1170] * v[598])
		- v[1057] * (v[1656] * v[360] + v[1170] * v[601]));
	v[1828] = v[1059] * (-(v[1054] * (v[1657] * v[358] + v[1170] * v[596])) - v[1056] * (v[1657] * v[359] + v[1170] * v[599])
		- v[1057] * (v[1657] * v[360] + v[1170] * v[602]));
	v[1830] = v[1059] * (-(v[1054] * (v[1632] * v[358] + v[1171] * v[594])) - v[1056] * (v[1632] * v[359] + v[1171] * v[597])
		- v[1057] * (v[1632] * v[360] + v[1171] * v[600]));
	v[1831] = v[1059] * (-(v[1054] * (v[1633] * v[358] + v[1171] * v[595])) - v[1056] * (v[1633] * v[359] + v[1171] * v[598])
		- v[1057] * (v[1633] * v[360] + v[1171] * v[601]));
	v[1832] = v[1059] * (-(v[1054] * (v[1634] * v[358] + v[1171] * v[596])) - v[1056] * (v[1634] * v[359] + v[1171] * v[599])
		- v[1057] * (v[1634] * v[360] + v[1171] * v[602]));
	v[2795] = v[1059] * (v[1057] * (-v[1548] + gti[0] * (v[1020] * v[1751] + v[1775] * v[322] + v[1727] * v[337]) + gti[1] *
		(v[1021] * v[1751] + v[1775] * v[324] + v[1727] * v[340]) + gti[2] * (v[1022] * v[1751] + v[1775] * v[325]
			+ v[1727] * v[342]) - v[1788] * v[360] - v[1173] * v[600] - v[1172] * v[601] - v[1048] * v[941]) + v[1056] * (-v[1554]
				+ gti[0] * (v[1020] * v[1753] + v[1777] * v[322] + v[1729] * v[337]) + gti[1] * (v[1021] * v[1753] + v[1777] * v[324]
					+ v[1729] * v[340]) + gti[2] * (v[1022] * v[1753] + v[1777] * v[325] + v[1729] * v[342]) - v[1788] * v[359]
				- v[1173] * v[597] - v[1172] * v[598] - v[1048] * v[953]) + v[1054] * (-v[1560] + gti[0] * (v[1020] * v[1755]
					+ v[1779] * v[322] + v[1731] * v[337]) + gti[1] * (v[1021] * v[1755] + v[1779] * v[324] + v[1731] * v[340]) + gti[2] *
					(v[1022] * v[1755] + v[1779] * v[325] + v[1731] * v[342]) - v[1788] * v[358] - v[1173] * v[594] - v[1172] * v[595]
					- v[1048] * v[965])) + (*epsn)*(v[603] * v[604] + v[606] * v[607] + v[609] * v[610] + v[1013] * v[616] + v[1014] * v[622]
						+ v[363] * v[943] + v[362] * v[955] + v[361] * v[967]);
	v[2796] = v[1059] * (v[1057] * (-v[1545] + gti[0] * (v[1020] * v[1740] + v[1764] * v[322] + v[1716] * v[337]) + gti[1] *
		(v[1021] * v[1740] + v[1764] * v[324] + v[1716] * v[340]) + gti[2] * (v[1022] * v[1740] + v[1764] * v[325]
			+ v[1716] * v[342]) - v[1785] * v[360] - v[1174] * v[600] - v[1172] * v[602] - v[1048] * v[935]) + v[1056] * (-v[1551]
				+ gti[0] * (v[1020] * v[1743] + v[1767] * v[322] + v[1719] * v[337]) + gti[1] * (v[1021] * v[1743] + v[1767] * v[324]
					+ v[1719] * v[340]) + gti[2] * (v[1022] * v[1743] + v[1767] * v[325] + v[1719] * v[342]) - v[1785] * v[359]
				- v[1174] * v[597] - v[1172] * v[599] - v[1048] * v[947]) + v[1054] * (-v[1557] + gti[0] * (v[1020] * v[1746]
					+ v[1770] * v[322] + v[1722] * v[337]) + gti[1] * (v[1021] * v[1746] + v[1770] * v[324] + v[1722] * v[340]) + gti[2] *
					(v[1022] * v[1746] + v[1770] * v[325] + v[1722] * v[342]) - v[1785] * v[358] - v[1174] * v[594] - v[1172] * v[596]
					- v[1048] * v[959])) + (*epsn)*(v[603] * v[605] + v[606] * v[608] + v[609] * v[611] + v[1016] * v[616] + v[1017] * v[622]
						+ v[363] * v[938] + v[362] * v[950] + v[361] * v[962]);
	v[2798] = v[1059] * (v[1057] * (-v[1546] + gti[0] * (v[1020] * v[1741] + v[1765] * v[322] + v[1717] * v[337]) + gti[1] *
		(v[1021] * v[1741] + v[1765] * v[324] + v[1717] * v[340]) + gti[2] * (v[1022] * v[1741] + v[1765] * v[325]
			+ v[1717] * v[342]) - v[1786] * v[360] - v[1174] * v[601] - v[1173] * v[602] - v[1048] * v[936]) + v[1056] * (-v[1552]
				+ gti[0] * (v[1020] * v[1744] + v[1768] * v[322] + v[1720] * v[337]) + gti[1] * (v[1021] * v[1744] + v[1768] * v[324]
					+ v[1720] * v[340]) + gti[2] * (v[1022] * v[1744] + v[1768] * v[325] + v[1720] * v[342]) - v[1786] * v[359]
				- v[1174] * v[598] - v[1173] * v[599] - v[1048] * v[948]) + v[1054] * (-v[1558] + gti[0] * (v[1020] * v[1747]
					+ v[1771] * v[322] + v[1723] * v[337]) + gti[1] * (v[1021] * v[1747] + v[1771] * v[324] + v[1723] * v[340]) + gti[2] *
					(v[1022] * v[1747] + v[1771] * v[325] + v[1723] * v[342]) - v[1786] * v[358] - v[1174] * v[595] - v[1173] * v[596]
					- v[1048] * v[960])) + (*epsn)*(v[604] * v[605] + v[607] * v[608] + v[610] * v[611] + v[1016] * v[617] + v[1017] * v[623]
						+ v[363] * v[939] + v[362] * v[951] + v[361] * v[963]);
	v[2775] = -(v[2688] / v[1846]);
	v[2774] = -(v[2691] / v[1846]);
	v[2380] = v[1917] / v[1846];
	v[2365] = v[1855] * v[2380];
	v[2362] = -(v[1904] / v[1846]);
	v[2387] = v[1853] * v[2362];
	v[2346] = -(v[1910] / v[1846]);
	v[2386] = v[1849] * v[2346];
	v[2345] = -(v[2692] / v[1846]);
	v[2347] = (*a4)*v[2345] + v[1855] * v[2346] + v[2378];
	v[2343] = -(v[2689] / v[1846]);
	v[2344] = (*a4)*v[2343] + v[1853] * v[2346] + v[2361];
	v[2341] = (v[1846] - v[782]) / v[1846];
	v[2342] = (*a4)*v[2341] + v[2365] + v[2387] + v[1849] * v[472];
	v[1848] = v[1849] * v[2341] + v[1853] * v[2343] + v[1855] * v[2345];
	v[2363] = -v[2361] + v[1849] * v[2362] + (*a4)*v[2774];
	v[2367] = -(v[2690] / v[1846]);
	v[2368] = v[1855] * v[2362] + (*a4)*v[2367] + v[2382];
	v[2364] = (v[1846] - v[779]) / v[1846];
	v[2366] = (*a4)*v[2364] + v[2365] + v[2386] + v[1853] * v[471];
	v[1856] = v[1853] * v[2364] + v[1855] * v[2367] + v[1849] * v[2774];
	v[2381] = -v[2378] + v[1849] * v[2380] + (*a4)*v[2775];
	v[2385] = (v[1842] + v[1846]) / v[1846];
	v[2388] = (*a4)*v[2385] + v[2386] + v[2387] + v[1855] * v[474];
	v[2383] = -(v[2687] / v[1846]);
	v[2384] = v[1853] * v[2380] - v[2382] + (*a4)*v[2383];
	v[1862] = v[1853] * v[2383] + v[1855] * v[2385] + v[1849] * v[2775];
	v[2570] = -(v[1848] * v[564]) - v[1856] * v[565] - v[1862] * v[566];
	v[2569] = -(v[1848] * v[555]) - v[1856] * v[556] - v[1862] * v[557];
	v[2568] = -(v[1848] * v[561]) - v[1856] * v[562] - v[1862] * v[563];
	v[2567] = -(v[1848] * v[552]) - v[1856] * v[553] - v[1862] * v[554];
	v[2566] = -(v[1848] * v[558]) - v[1856] * v[559] - v[1862] * v[560];
	v[2576] = (*cn)*(v[2566] * v[628] + v[2568] * v[629] + v[2570] * v[990]);
	v[2574] = (*cn)*(v[2566] * v[626] + v[2570] * v[629] + v[2568] * v[981]);
	v[2572] = (*cn)*(v[2568] * v[626] + v[2570] * v[628] + v[2566] * v[971]);
	v[2565] = -(v[1848] * v[549]) - v[1856] * v[550] - v[1862] * v[551];
	v[2575] = (*cn)*(v[2565] * v[628] + v[2567] * v[629] + v[2569] * v[990]);
	v[2573] = (*cn)*(v[2565] * v[626] + v[2569] * v[629] + v[2567] * v[981]);
	v[2571] = (*cn)*(v[2567] * v[626] + v[2569] * v[628] + v[2565] * v[971]);
	v[2406] = v[1848] * v[2228] + v[1856] * v[2260] + v[1862] * v[2294] + v[2784];
	v[2405] = v[1848] * v[2225] + v[1856] * v[2257] + v[1862] * v[2291] + v[2404];
	v[2403] = v[1848] * v[2222] + v[1856] * v[2254] + v[1862] * v[2288] + v[2402];
	v[2398] = v[1848] * v[2238] + v[1856] * v[2271] + v[1862] * v[2306] + v[2404];
	v[2397] = v[1848] * v[2237] + v[1856] * v[2270] + v[1862] * v[2305] + v[2781];
	v[2396] = v[1848] * v[2234] + v[1856] * v[2267] + v[1862] * v[2302] + v[2395];
	v[2391] = v[1848] * v[2246] + v[1856] * v[2280] + v[1862] * v[2316] + v[2402];
	v[2783] = (*cn)*v[2404] + (*epsn)*v[629];
	v[2782] = (*cn)*v[2402] + (*epsn)*v[628];
	v[2390] = v[1848] * v[2245] + v[1856] * v[2279] + v[1862] * v[2315] + v[2395];
	v[2780] = (*cn)*v[2395] + (*epsn)*v[626];
	v[2389] = v[1848] * v[2244] + v[1856] * v[2278] + v[1862] * v[2314] + v[2779];
	v[2776] = (*a4)*d[0] - (*a4)*d[6] - (*a6)*dduiP[0] + (*a6)*dduiS[0] - (*a5)*duiP[0] + (*a5)*duiS[0];
	v[2777] = (*a4)*d[1] - (*a4)*d[7] - (*a6)*dduiP[1] + (*a6)*dduiS[1] - (*a5)*duiP[1] + (*a5)*duiS[1];
	v[2778] = (*a4)*d[2] - (*a4)*d[8] - (*a6)*dduiP[2] + (*a6)*dduiS[2] - (*a5)*duiP[2] + (*a5)*duiS[2];
	v[2409] = v[1848] * v[2231] + v[1856] * v[2264] + v[1862] * v[2299] + v[2198] * v[2776] + v[2210] * v[2777]
		+ v[2216] * v[2778] + v[2347] * v[992] + v[2368] * v[994] + v[2388] * v[996];
	v[2408] = v[1848] * v[2230] + v[1856] * v[2263] + v[1862] * v[2298] + v[2197] * v[2776] + v[2209] * v[2777]
		+ v[2215] * v[2778] + v[2344] * v[992] + v[2366] * v[994] + v[2384] * v[996];
	v[2407] = v[1848] * v[2229] + v[1856] * v[2262] + v[1862] * v[2296] + v[2196] * v[2776] + v[2208] * v[2777]
		+ v[2214] * v[2778] + v[2342] * v[992] + v[2363] * v[994] + v[2381] * v[996];
	v[2401] = v[1848] * v[2241] + v[1856] * v[2275] + v[1862] * v[2311] + v[2201] * v[2776] + v[2207] * v[2777]
		+ v[2210] * v[2778] + v[2347] * v[984] + v[2368] * v[986] + v[2388] * v[988];
	v[2400] = v[1848] * v[2240] + v[1856] * v[2274] + v[1862] * v[2310] + v[2200] * v[2776] + v[2206] * v[2777]
		+ v[2209] * v[2778] + v[2344] * v[984] + v[2366] * v[986] + v[2384] * v[988];
	v[2399] = v[1848] * v[2239] + v[1856] * v[2273] + v[1862] * v[2308] + v[2199] * v[2776] + v[2205] * v[2777]
		+ v[2208] * v[2778] + v[2342] * v[984] + v[2363] * v[986] + v[2381] * v[988];
	v[2394] = v[1848] * v[2249] + v[1856] * v[2284] + v[1862] * v[2321] + v[2195] * v[2776] + v[2201] * v[2777]
		+ v[2198] * v[2778] + v[2347] * v[975] + v[2368] * v[977] + v[2388] * v[979];
	v[2393] = v[1848] * v[2248] + v[1856] * v[2283] + v[1862] * v[2320] + v[2194] * v[2776] + v[2200] * v[2777]
		+ v[2197] * v[2778] + v[2344] * v[975] + v[2366] * v[977] + v[2384] * v[979];
	v[2392] = v[1848] * v[2247] + v[1856] * v[2282] + v[1862] * v[2318] + v[2193] * v[2776] + v[2199] * v[2777]
		+ v[2196] * v[2778] + v[2342] * v[975] + v[2363] * v[977] + v[2381] * v[979];
	v[1887] = v[2777] * v[626] + v[2778] * v[628] + v[2776] * v[971] + v[1848] * v[975] + v[1856] * v[977] + v[1862] * v[979];
	v[1888] = v[2776] * v[626] + v[2778] * v[629] + v[2777] * v[981] + v[1848] * v[984] + v[1856] * v[986] + v[1862] * v[988];
	v[1889] = v[2776] * v[628] + v[2777] * v[629] + v[2778] * v[990] + v[1848] * v[992] + v[1856] * v[994] + v[1862] * v[996];
	v[2582] = (*cn)*(-(v[1887] * v[560]) - v[1888] * v[563] - v[1889] * v[566] + v[2566] * v[979] + v[2568] * v[988]
		+ v[2570] * v[996]);
	v[2581] = (*cn)*(-(v[1887] * v[551]) - v[1888] * v[554] - v[1889] * v[557] + v[2565] * v[979] + v[2567] * v[988]
		+ v[2569] * v[996]);
	v[2580] = (*cn)*(-(v[1887] * v[559]) - v[1888] * v[562] - v[1889] * v[565] + v[2566] * v[977] + v[2568] * v[986]
		+ v[2570] * v[994]);
	v[2579] = (*cn)*(-(v[1887] * v[550]) - v[1888] * v[553] - v[1889] * v[556] + v[2565] * v[977] + v[2567] * v[986]
		+ v[2569] * v[994]);
	v[2578] = (*cn)*(-(v[1887] * v[558]) - v[1888] * v[561] - v[1889] * v[564] + v[2566] * v[975] + v[2568] * v[984]
		+ v[2570] * v[992]);
	v[2577] = (*cn)*(-(v[1887] * v[549]) - v[1888] * v[552] - v[1889] * v[555] + v[2565] * v[975] + v[2567] * v[984]
		+ v[2569] * v[992]);
	v[2619] = v[1058] * v[1184] + v[1060] * v[1185] + v[1061] * v[1186] + (*cn)*(v[1888] * v[626] + v[1889] * v[628]
		+ v[1887] * v[971]) + (*epsn)*(v[362] * v[626] + v[363] * v[628] + v[361] * v[971]);
	v[2620] = v[1058] * v[1185] + v[1060] * v[1202] + v[1061] * v[1203] + (*cn)*(v[1887] * v[626] + v[1889] * v[629]
		+ v[1888] * v[981]) + (*epsn)*(v[361] * v[626] + v[363] * v[629] + v[362] * v[981]);
	v[2621] = v[1058] * v[1186] + v[1060] * v[1203] + v[1061] * v[1219] + (*cn)*(v[1887] * v[628] + v[1888] * v[629]
		+ v[1889] * v[990]) + (*epsn)*(v[361] * v[628] + v[362] * v[629] + v[363] * v[990]);
	v[2625] = v[1184] * v[1345] + v[1185] * v[1346] + v[1186] * v[1347] + (*cn)*v[2779] + v[2571] * v[613] + v[2572] * v[619]
		+ (*epsn)*v[971];
	v[2626] = v[1184] * v[1349] + v[1185] * v[1350] + v[1186] * v[1351] + v[2780] + v[2571] * v[614] + v[2572] * v[620];
	v[2627] = v[1184] * v[1353] + v[1185] * v[1354] + v[1186] * v[1355] + v[2782] + v[2571] * v[615] + v[2572] * v[621];
	v[2628] = v[1184] * v[1371] + v[1185] * v[1372] + v[1186] * v[1373] + v[2785] + v[2571] * v[616] + v[2572] * v[622] + (*cn
		)*(v[1887] * v[2193] + v[1889] * v[2196] + v[1888] * v[2199] + v[2399] * v[626] + v[2407] * v[628] + v[2392] * v[971]);
	v[2629] = v[1184] * v[1375] + v[1185] * v[1376] + v[1186] * v[1377] + v[2788] + v[2571] * v[617] + v[2572] * v[623] + (*cn
		)*(v[1887] * v[2194] + v[1889] * v[2197] + v[1888] * v[2200] + v[2400] * v[626] + v[2408] * v[628] + v[2393] * v[971]);
	v[2630] = v[1184] * v[1379] + v[1185] * v[1380] + v[1186] * v[1381] + v[2791] + v[2571] * v[618] + v[2572] * v[624] + (*cn
		)*(v[1887] * v[2195] + v[1889] * v[2198] + v[1888] * v[2201] + v[2401] * v[626] + v[2409] * v[628] + v[2394] * v[971]);
	v[2631] = v[1185] * v[1345] + v[1202] * v[1346] + v[1203] * v[1347] + v[2780] + v[2573] * v[613] + v[2574] * v[619];
	v[2632] = v[1185] * v[1349] + v[1202] * v[1350] + v[1203] * v[1351] + (*cn)*v[2781] + v[2573] * v[614] + v[2574] * v[620]
		+ (*epsn)*v[981];
	v[2633] = v[1185] * v[1353] + v[1202] * v[1354] + v[1203] * v[1355] + v[2783] + v[2573] * v[615] + v[2574] * v[621];
	v[2634] = v[1185] * v[1371] + v[1202] * v[1372] + v[1203] * v[1373] + v[2786] + v[2573] * v[616] + v[2574] * v[622] + (*cn
		)*(v[1887] * v[2199] + v[1888] * v[2205] + v[1889] * v[2208] + v[2392] * v[626] + v[2407] * v[629] + v[2399] * v[981]);
	v[2635] = v[1185] * v[1375] + v[1202] * v[1376] + v[1203] * v[1377] + v[2789] + v[2573] * v[617] + v[2574] * v[623] + (*cn
		)*(v[1887] * v[2200] + v[1888] * v[2206] + v[1889] * v[2209] + v[2393] * v[626] + v[2408] * v[629] + v[2400] * v[981]);
	v[2636] = v[1185] * v[1379] + v[1202] * v[1380] + v[1203] * v[1381] + v[2792] + v[2573] * v[618] + v[2574] * v[624] + (*cn
		)*(v[1887] * v[2201] + v[1888] * v[2207] + v[1889] * v[2210] + v[2394] * v[626] + v[2409] * v[629] + v[2401] * v[981]);
	v[2637] = v[1186] * v[1345] + v[1203] * v[1346] + v[1219] * v[1347] + v[2782] + v[2575] * v[613] + v[2576] * v[619];
	v[2638] = v[1186] * v[1349] + v[1203] * v[1350] + v[1219] * v[1351] + v[2783] + v[2575] * v[614] + v[2576] * v[620];
	v[2639] = v[1186] * v[1353] + v[1203] * v[1354] + v[1219] * v[1355] + (*cn)*v[2784] + v[2575] * v[615] + v[2576] * v[621]
		+ (*epsn)*v[990];
	v[2640] = v[1186] * v[1371] + v[1203] * v[1372] + v[1219] * v[1373] + v[2787] + v[2575] * v[616] + v[2576] * v[622] + (*cn
		)*(v[1887] * v[2196] + v[1888] * v[2208] + v[1889] * v[2214] + v[2392] * v[628] + v[2399] * v[629] + v[2407] * v[990]);
	v[2641] = v[1186] * v[1375] + v[1203] * v[1376] + v[1219] * v[1377] + v[2790] + v[2575] * v[617] + v[2576] * v[623] + (*cn
		)*(v[1887] * v[2197] + v[1888] * v[2209] + v[1889] * v[2215] + v[2393] * v[628] + v[2400] * v[629] + v[2408] * v[990]);
	v[2642] = v[1186] * v[1379] + v[1203] * v[1380] + v[1219] * v[1381] + v[2793] + v[2575] * v[618] + v[2576] * v[624] + (*cn
		)*(v[1887] * v[2198] + v[1888] * v[2210] + v[1889] * v[2216] + v[2394] * v[628] + v[2401] * v[629] + v[2409] * v[990]);
	v[2644] = v[1361] * v[1397] + v[1363] * v[1398] + v[1365] * v[1399] + v[1059] * (v[1054] * (v[1584] - v[1653] * v[358])
		+ v[1056] * (v[1578] - v[1653] * v[359]) + v[1057] * (v[1572] - v[1653] * v[360]));
	v[2645] = v[1367] * v[1397] + v[1368] * v[1398] + v[1369] * v[1399] + v[1059] * (v[1054] * (v[1581] - v[1629] * v[358])
		+ v[1056] * (v[1575] - v[1629] * v[359]) + v[1057] * (v[1569] - v[1629] * v[360]));
	v[2650] = v[1367] * v[1409] + v[1368] * v[1410] + v[1369] * v[1411] + v[1059] * (v[1054] * (v[1582] - v[1630] * v[358])
		+ v[1056] * (v[1576] - v[1630] * v[359]) + v[1057] * (v[1570] - v[1630] * v[360]));
	v[2658] = v[1190] * v[1345] + v[1207] * v[1346] + v[1223] * v[1347] + v[2785] + v[2577] * v[613] + v[2578] * v[619] + (*cn
		)*(v[1889] * v[2222] + v[1888] * v[2234] + v[1887] * v[2244] + v[2389] * v[975] + v[2396] * v[984] + v[2403] * v[992]);
	v[2659] = v[1190] * v[1349] + v[1207] * v[1350] + v[1223] * v[1351] + v[2786] + v[2577] * v[614] + v[2578] * v[620] + (*cn
		)*(v[1889] * v[2225] + v[1888] * v[2237] + v[1887] * v[2245] + v[2390] * v[975] + v[2397] * v[984] + v[2405] * v[992]);
	v[2660] = v[1190] * v[1353] + v[1207] * v[1354] + v[1223] * v[1355] + v[2787] + v[2577] * v[615] + v[2578] * v[621] + (*cn
		)*(v[1889] * v[2228] + v[1888] * v[2238] + v[1887] * v[2246] + v[2391] * v[975] + v[2398] * v[984] + v[2406] * v[992]);
	v[2667] = v[1191] * v[1345] + v[1208] * v[1346] + v[1224] * v[1347] + v[2788] + v[2579] * v[613] + v[2580] * v[619] + (*cn
		)*(v[1889] * v[2254] + v[1888] * v[2267] + v[1887] * v[2278] + v[2389] * v[977] + v[2396] * v[986] + v[2403] * v[994]);
	v[2668] = v[1191] * v[1349] + v[1208] * v[1350] + v[1224] * v[1351] + v[2789] + v[2579] * v[614] + v[2580] * v[620] + (*cn
		)*(v[1889] * v[2257] + v[1888] * v[2270] + v[1887] * v[2279] + v[2390] * v[977] + v[2397] * v[986] + v[2405] * v[994]);
	v[2669] = v[1191] * v[1353] + v[1208] * v[1354] + v[1224] * v[1355] + v[2790] + v[2579] * v[615] + v[2580] * v[621] + (*cn
		)*(v[1889] * v[2260] + v[1888] * v[2271] + v[1887] * v[2280] + v[2391] * v[977] + v[2398] * v[986] + v[2406] * v[994]);
	v[2676] = v[1192] * v[1345] + v[1209] * v[1346] + v[1225] * v[1347] + v[2791] + v[2581] * v[613] + v[2582] * v[619] + (*cn
		)*(v[1889] * v[2288] + v[1888] * v[2302] + v[1887] * v[2314] + v[2389] * v[979] + v[2396] * v[988] + v[2403] * v[996]);
	v[2677] = v[1192] * v[1349] + v[1209] * v[1350] + v[1225] * v[1351] + v[2792] + v[2581] * v[614] + v[2582] * v[620] + (*cn
		)*(v[1889] * v[2291] + v[1888] * v[2305] + v[1887] * v[2315] + v[2390] * v[979] + v[2397] * v[988] + v[2405] * v[996]);
	v[2678] = v[1192] * v[1353] + v[1209] * v[1354] + v[1225] * v[1355] + v[2793] + v[2581] * v[615] + v[2582] * v[621] + (*cn
		)*(v[1889] * v[2294] + v[1888] * v[2306] + v[1887] * v[2316] + v[2391] * v[979] + v[2398] * v[988] + v[2406] * v[996]);
	Rc[0] = v[2619];
	Rc[1] = v[2620];
	Rc[2] = v[2621];
	Rc[3] = v[1058] * v[1187] + v[1060] * v[1204] + v[1061] * v[1220];
	Rc[4] = v[1058] * v[1188] + v[1060] * v[1205] + v[1061] * v[1221];
	Rc[5] = v[1058] * v[1189] + v[1060] * v[1206] + v[1061] * v[1222];
	Rc[6] = -v[2619];
	Rc[7] = -v[2620];
	Rc[8] = -v[2621];
	Rc[9] = v[1058] * v[1190] + v[1060] * v[1207] + v[1061] * v[1223] + (*cn)*(v[1887] * v[975] + v[1888] * v[984]
		+ v[1889] * v[992]) + (*epsn)*(v[361] * v[975] + v[362] * v[984] + v[363] * v[992]);
	Rc[10] = v[1058] * v[1191] + v[1060] * v[1208] + v[1061] * v[1224] + (*cn)*(v[1887] * v[977] + v[1888] * v[986]
		+ v[1889] * v[994]) + (*epsn)*(v[361] * v[977] + v[362] * v[986] + v[363] * v[994]);
	Rc[11] = v[1058] * v[1192] + v[1060] * v[1209] + v[1061] * v[1225] + (*cn)*(v[1887] * v[979] + v[1888] * v[988]
		+ v[1889] * v[996]) + (*epsn)*(v[361] * v[979] + v[362] * v[988] + v[363] * v[996]);
	Kc[0][0] = v[2625];
	Kc[0][1] = v[2626];
	Kc[0][2] = v[2627];
	Kc[0][3] = v[1336];
	Kc[0][4] = v[1340];
	Kc[0][5] = v[1341];
	Kc[0][6] = -v[2625];
	Kc[0][7] = -v[2626];
	Kc[0][8] = -v[2627];
	Kc[0][9] = v[2628];
	Kc[0][10] = v[2629];
	Kc[0][11] = v[2630];
	Kc[1][0] = v[2631];
	Kc[1][1] = v[2632];
	Kc[1][2] = v[2633];
	Kc[1][3] = v[1360];
	Kc[1][4] = v[1366];
	Kc[1][5] = v[1370];
	Kc[1][6] = -v[2631];
	Kc[1][7] = -v[2632];
	Kc[1][8] = -v[2633];
	Kc[1][9] = v[2634];
	Kc[1][10] = v[2635];
	Kc[1][11] = v[2636];
	Kc[2][0] = v[2637];
	Kc[2][1] = v[2638];
	Kc[2][2] = v[2639];
	Kc[2][3] = v[1386];
	Kc[2][4] = v[1388];
	Kc[2][5] = v[1389];
	Kc[2][6] = -v[2637];
	Kc[2][7] = -v[2638];
	Kc[2][8] = -v[2639];
	Kc[2][9] = v[2640];
	Kc[2][10] = v[2641];
	Kc[2][11] = v[2642];
	Kc[3][0] = v[1393];
	Kc[3][1] = v[1394];
	Kc[3][2] = v[1395];
	Kc[3][3] = v[1357] * v[1397] + v[1358] * v[1398] + v[1359] * v[1399] + v[1059] * (v[1054] * (v[1586] - v[1673] * v[358])
		+ v[1056] * (v[1580] - v[1673] * v[359]) + v[1057] * (v[1574] - v[1673] * v[360]));
	Kc[3][4] = v[2644];
	Kc[3][5] = v[2645];
	Kc[3][6] = -v[1393];
	Kc[3][7] = -v[1394];
	Kc[3][8] = -v[1395];
	Kc[3][9] = v[1187] * v[1371] + v[1204] * v[1372] + v[1220] * v[1373] + v[1821];
	Kc[3][10] = v[1187] * v[1375] + v[1204] * v[1376] + v[1220] * v[1377] + v[1822];
	Kc[3][11] = v[1187] * v[1379] + v[1204] * v[1380] + v[1220] * v[1381] + v[1823];
	Kc[4][0] = v[1405];
	Kc[4][1] = v[1406];
	Kc[4][2] = v[1407];
	Kc[4][3] = v[2644];
	Kc[4][4] = v[1361] * v[1409] + v[1363] * v[1410] + v[1365] * v[1411] + v[1059] * (v[1054] * (v[1585] - v[1654] * v[358])
		+ v[1056] * (v[1579] - v[1654] * v[359]) + v[1057] * (v[1573] - v[1654] * v[360]));
	Kc[4][5] = v[2650];
	Kc[4][6] = -v[1405];
	Kc[4][7] = -v[1406];
	Kc[4][8] = -v[1407];
	Kc[4][9] = v[1188] * v[1371] + v[1205] * v[1372] + v[1221] * v[1373] + v[1826];
	Kc[4][10] = v[1188] * v[1375] + v[1205] * v[1376] + v[1221] * v[1377] + v[1827];
	Kc[4][11] = v[1188] * v[1379] + v[1205] * v[1380] + v[1221] * v[1381] + v[1828];
	Kc[5][0] = v[1416];
	Kc[5][1] = v[1417];
	Kc[5][2] = v[1418];
	Kc[5][3] = v[2645];
	Kc[5][4] = v[2650];
	Kc[5][5] = v[1059] * (v[1249] * (v[1189] * v[1367] + v[1206] * v[1368] + v[1222] * v[1369]) + v[1054] * (v[1583]
		- v[1631] * v[358]) + v[1056] * (v[1577] - v[1631] * v[359]) + v[1057] * (v[1571] - v[1631] * v[360]));
	Kc[5][6] = -v[1416];
	Kc[5][7] = -v[1417];
	Kc[5][8] = -v[1418];
	Kc[5][9] = v[1189] * v[1371] + v[1206] * v[1372] + v[1222] * v[1373] + v[1830];
	Kc[5][10] = v[1189] * v[1375] + v[1206] * v[1376] + v[1222] * v[1377] + v[1831];
	Kc[5][11] = v[1189] * v[1379] + v[1206] * v[1380] + v[1222] * v[1381] + v[1832];
	Kc[6][0] = -v[2625];
	Kc[6][1] = -v[2626];
	Kc[6][2] = -v[2627];
	Kc[6][3] = -v[1336];
	Kc[6][4] = -v[1340];
	Kc[6][5] = -v[1341];
	Kc[6][6] = v[2625];
	Kc[6][7] = v[2626];
	Kc[6][8] = v[2627];
	Kc[6][9] = -v[2628];
	Kc[6][10] = -v[2629];
	Kc[6][11] = -v[2630];
	Kc[7][0] = -v[2631];
	Kc[7][1] = -v[2632];
	Kc[7][2] = -v[2633];
	Kc[7][3] = -v[1360];
	Kc[7][4] = -v[1366];
	Kc[7][5] = -v[1370];
	Kc[7][6] = v[2631];
	Kc[7][7] = v[2632];
	Kc[7][8] = v[2633];
	Kc[7][9] = -v[2634];
	Kc[7][10] = -v[2635];
	Kc[7][11] = -v[2636];
	Kc[8][0] = -v[2637];
	Kc[8][1] = -v[2638];
	Kc[8][2] = -v[2639];
	Kc[8][3] = -v[1386];
	Kc[8][4] = -v[1388];
	Kc[8][5] = -v[1389];
	Kc[8][6] = v[2637];
	Kc[8][7] = v[2638];
	Kc[8][8] = v[2639];
	Kc[8][9] = -v[2640];
	Kc[8][10] = -v[2641];
	Kc[8][11] = -v[2642];
	Kc[9][0] = v[2658];
	Kc[9][1] = v[2659];
	Kc[9][2] = v[2660];
	Kc[9][3] = v[1357] * v[1427] + v[1358] * v[1428] + v[1359] * v[1429] + v[1821];
	Kc[9][4] = v[1361] * v[1427] + v[1363] * v[1428] + v[1365] * v[1429] + v[1826];
	Kc[9][5] = v[1367] * v[1427] + v[1368] * v[1428] + v[1369] * v[1429] + v[1830];
	Kc[9][6] = -v[2658];
	Kc[9][7] = -v[2659];
	Kc[9][8] = -v[2660];
	Kc[9][9] = v[1190] * v[1371] + v[1207] * v[1372] + v[1223] * v[1373] + v[2577] * v[616] + v[2578] * v[622] + v[1059] *
		(v[1057] * (-v[1550] + gti[0] * (v[1020] * v[1758] + v[1782] * v[322] + v[1734] * v[337]) + gti[1] * (v[1021] * v[1758]
			+ v[1782] * v[324] + v[1734] * v[340]) + gti[2] * (v[1022] * v[1758] + v[1782] * v[325] + v[1734] * v[342])
			- v[1790] * v[360] + v[2794] * v[600] - v[1048] * v[945]) + v[1056] * (-v[1556] + gti[0] * (v[1020] * v[1759]
				+ v[1783] * v[322] + v[1735] * v[337]) + gti[1] * (v[1021] * v[1759] + v[1783] * v[324] + v[1735] * v[340]) + gti[2] *
				(v[1022] * v[1759] + v[1783] * v[325] + v[1735] * v[342]) - v[1790] * v[359] - 2e0*v[1172] * v[597] - v[1048] * v[957])
			+ v[1054] * (-v[1562] + gti[0] * (v[1020] * v[1760] + v[1784] * v[322] + v[1736] * v[337]) + gti[1] * (v[1021] * v[1760]
				+ v[1784] * v[324] + v[1736] * v[340]) + gti[2] * (v[1022] * v[1760] + v[1784] * v[325] + v[1736] * v[342])
				- v[1790] * v[358] + v[2794] * v[594] - v[1048] * v[969])) + (*epsn)*((v[603] * v[603]) + (v[606] * v[606]) +
				(v[609] * v[609]) + (-(v[361] * v[549]) - v[362] * v[552] - v[363] * v[555] - v[343] * v[603] - v[344] * v[606]
					- v[345] * v[609])*v[616] + (-(v[361] * v[558]) - v[362] * v[561] - v[363] * v[564] - v[350] * v[603] - v[352] * v[606]
						- v[354] * v[609])*v[622] + v[363] * v[946] + v[362] * v[958] + v[361] * v[970]) + (*cn)*(v[1889] * v[2229]
							+ v[1888] * v[2239] + v[1887] * v[2247] + v[2392] * v[975] + v[2399] * v[984] + v[2407] * v[992]);
	Kc[9][10] = v[1190] * v[1375] + v[1207] * v[1376] + v[1223] * v[1377] + v[2795] + v[2577] * v[617] + v[2578] * v[623] +
		(*cn)*(v[1889] * v[2230] + v[1888] * v[2240] + v[1887] * v[2248] + v[2393] * v[975] + v[2400] * v[984] + v[2408] * v[992]
			);
	Kc[9][11] = v[1190] * v[1379] + v[1207] * v[1380] + v[1223] * v[1381] + v[2796] + v[2577] * v[618] + v[2578] * v[624] +
		(*cn)*(v[1889] * v[2231] + v[1888] * v[2241] + v[1887] * v[2249] + v[2394] * v[975] + v[2401] * v[984] + v[2409] * v[992]
			);
	Kc[10][0] = v[2667];
	Kc[10][1] = v[2668];
	Kc[10][2] = v[2669];
	Kc[10][3] = v[1357] * v[1439] + v[1358] * v[1440] + v[1359] * v[1441] + v[1822];
	Kc[10][4] = v[1361] * v[1439] + v[1363] * v[1440] + v[1365] * v[1441] + v[1827];
	Kc[10][5] = v[1367] * v[1439] + v[1368] * v[1440] + v[1369] * v[1441] + v[1831];
	Kc[10][6] = -v[2667];
	Kc[10][7] = -v[2668];
	Kc[10][8] = -v[2669];
	Kc[10][9] = v[1191] * v[1371] + v[1208] * v[1372] + v[1224] * v[1373] + v[2795] + v[2579] * v[616] + v[2580] * v[622] +
		(*cn)*(v[1889] * v[2262] + v[1888] * v[2273] + v[1887] * v[2282] + v[2392] * v[977] + v[2399] * v[986] + v[2407] * v[994]
			);
	Kc[10][10] = v[1191] * v[1375] + v[1208] * v[1376] + v[1224] * v[1377] + v[2579] * v[617] + v[2580] * v[623] + v[1059] *
		(v[1057] * (-v[1549] + gti[0] * (v[1020] * v[1752] + v[1776] * v[322] + v[1728] * v[337]) + gti[1] * (v[1021] * v[1752]
			+ v[1776] * v[324] + v[1728] * v[340]) + gti[2] * (v[1022] * v[1752] + v[1776] * v[325] + v[1728] * v[342])
			- v[1789] * v[360] + v[2797] * v[601] - v[1048] * v[942]) + v[1056] * (-v[1555] + gti[0] * (v[1020] * v[1754]
				+ v[1778] * v[322] + v[1730] * v[337]) + gti[1] * (v[1021] * v[1754] + v[1778] * v[324] + v[1730] * v[340]) + gti[2] *
				(v[1022] * v[1754] + v[1778] * v[325] + v[1730] * v[342]) - v[1789] * v[359] - 2e0*v[1173] * v[598] - v[1048] * v[954])
			+ v[1054] * (-v[1561] + gti[0] * (v[1020] * v[1756] + v[1780] * v[322] + v[1732] * v[337]) + gti[1] * (v[1021] * v[1756]
				+ v[1780] * v[324] + v[1732] * v[340]) + gti[2] * (v[1022] * v[1756] + v[1780] * v[325] + v[1732] * v[342])
				- v[1789] * v[358] + v[2797] * v[595] - v[1048] * v[966])) + (*epsn)*((v[604] * v[604]) + (v[607] * v[607]) +
				(v[610] * v[610]) + v[1013] * v[617] + v[1014] * v[623] + v[363] * v[944] + v[362] * v[956] + v[361] * v[968]) + (*cn)*
					(v[1889] * v[2263] + v[1888] * v[2274] + v[1887] * v[2283] + v[2393] * v[977] + v[2400] * v[986] + v[2408] * v[994]);
	Kc[10][11] = v[1191] * v[1379] + v[1208] * v[1380] + v[1224] * v[1381] + v[2798] + v[2579] * v[618] + v[2580] * v[624] +
		(*cn)*(v[1889] * v[2264] + v[1888] * v[2275] + v[1887] * v[2284] + v[2394] * v[977] + v[2401] * v[986] + v[2409] * v[994]
			);
	Kc[11][0] = v[2676];
	Kc[11][1] = v[2677];
	Kc[11][2] = v[2678];
	Kc[11][3] = v[1357] * v[1451] + v[1358] * v[1452] + v[1359] * v[1453] + v[1823];
	Kc[11][4] = v[1361] * v[1451] + v[1363] * v[1452] + v[1365] * v[1453] + v[1828];
	Kc[11][5] = v[1367] * v[1451] + v[1368] * v[1452] + v[1369] * v[1453] + v[1832];
	Kc[11][6] = -v[2676];
	Kc[11][7] = -v[2677];
	Kc[11][8] = -v[2678];
	Kc[11][9] = v[1192] * v[1371] + v[1209] * v[1372] + v[1225] * v[1373] + v[2796] + v[2581] * v[616] + v[2582] * v[622] +
		(*cn)*(v[1889] * v[2296] + v[1888] * v[2308] + v[1887] * v[2318] + v[2392] * v[979] + v[2399] * v[988] + v[2407] * v[996]
			);
	Kc[11][10] = v[1192] * v[1375] + v[1209] * v[1376] + v[1225] * v[1377] + v[2798] + v[2581] * v[617] + v[2582] * v[623] +
		(*cn)*(v[1889] * v[2298] + v[1888] * v[2310] + v[1887] * v[2320] + v[2393] * v[979] + v[2400] * v[988] + v[2408] * v[996]
			);
	Kc[11][11] = v[1192] * v[1379] + v[1209] * v[1380] + v[1225] * v[1381] + v[2581] * v[618] + v[2582] * v[624] + v[1059] *
		(v[1057] * (-v[1547] + gti[0] * (v[1020] * v[1742] + v[1766] * v[322] + v[1718] * v[337]) + gti[1] * (v[1021] * v[1742]
			+ v[1766] * v[324] + v[1718] * v[340]) + gti[2] * (v[1022] * v[1742] + v[1766] * v[325] + v[1718] * v[342])
			- v[1787] * v[360] + v[2799] * v[602] - v[1048] * v[937]) + v[1056] * (-v[1553] + gti[0] * (v[1020] * v[1745]
				+ v[1769] * v[322] + v[1721] * v[337]) + gti[1] * (v[1021] * v[1745] + v[1769] * v[324] + v[1721] * v[340]) + gti[2] *
				(v[1022] * v[1745] + v[1769] * v[325] + v[1721] * v[342]) - v[1787] * v[359] - 2e0*v[1174] * v[599] - v[1048] * v[949])
			+ v[1054] * (-v[1559] + gti[0] * (v[1020] * v[1748] + v[1772] * v[322] + v[1724] * v[337]) + gti[1] * (v[1021] * v[1748]
				+ v[1772] * v[324] + v[1724] * v[340]) + gti[2] * (v[1022] * v[1748] + v[1772] * v[325] + v[1724] * v[342])
				- v[1787] * v[358] + v[2799] * v[596] - v[1048] * v[961])) + (*epsn)*((v[605] * v[605]) + (v[608] * v[608]) +
				(v[611] * v[611]) + v[1016] * v[618] + v[1017] * v[624] + v[363] * v[940] + v[362] * v[952] + v[361] * v[964]) + (*cn)*
					(v[1889] * v[2299] + v[1888] * v[2311] + v[1887] * v[2321] + v[2394] * v[979] + v[2401] * v[988] + v[2409] * v[996]);
}

