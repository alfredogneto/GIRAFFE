#include "Pipe_1.h"

#include "LagrangeSave.h"
#include "CoordinateSystem.h"
#include "Node.h"
#include "PipeSection.h"
#include "SolidSection.h"
#include "PostFiles.h"
#include "Encoding.h"
#include "Environment.h"
#include "Dynamic.h"
#include "Load.h"


#include"Database.h"
//Variáveis globais
extern
Database db;
#define PI 3.1415926535897932384626433832795

Pipe_1::Pipe_1()
{
	strain_energy = 0.0;
	kinetic_energy = 0.0;
	potential_gravitational_energy = 0.0;

	VTK_type = 21;
	nDOFs = 18;
	material = 0;
	section = 0;
	n_nodes = 3;
	number = 0;
	nodes = new int[n_nodes];
	VTK_nodes = new int[n_nodes];
	DOFs = new int *[n_nodes];
	for (int i = 0; i < n_nodes; i++)
		DOFs[i] = new int[db.number_GLs_node];
	type_name = new char[20];//Nome do tipo do elemento
	sprintf(type_name, "Pipe_1");
	
	VTK_nodes[0] = 0;
	VTK_nodes[1] = 2;
	VTK_nodes[2] = 1;

	rho_adt = 0;
	rho_adn = 0;

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
		DOFs[i][3] = 1;
		DOFs[i][4] = 1;
		DOFs[i][5] = 1;
	}
	//Variáveis internas do elemento
	alpha1 = 0;
	length = 0;
	jacobian = 0;
	alpha_escalar_delta = 0;
	g = 0;
	//Cada variável carregará a informação para seu ponto de Gauss
	N1 = new double[2];
	N2 = new double[2];
	N3 = new double[2];
	dN1 = new double[2];
	dN2 = new double[2];
	dN3 = new double[2];
	csi = new double[2];
	deltaN = new Matrix*[2];
	N = new Matrix*[2];
	alpha_delta = new Matrix*[2];
	u_delta = new Matrix*[2];
	d_alpha_delta = new Matrix*[2];
	d_u_delta = new Matrix*[2];
	A_delta = new Matrix*[2];
	Q_delta = new Matrix*[2];
	Xi_delta = new Matrix*[2];
	d_A_delta = new Matrix*[2];
	d_Xi_delta = new Matrix*[2];
	d_z = new Matrix*[2];
	d_Z = new Matrix*[2];
	eta_r = new Matrix*[2];
	kappa_r = new Matrix*[2];
	epsilon_r = new Matrix*[2];
	sigma_r = new Matrix*[2];
	n_r = new Matrix*[2];
	m_r = new Matrix*[2];
	n = new Matrix*[2];
	m = new Matrix*[2];
	//Loop nos pontos de Gauss
	for (int i = 0; i < 2; i++)
	{
		deltaN[i] = new Matrix(9, 18);
		N[i] = new Matrix(6, 18);
		alpha_delta[i] = new Matrix(3);
		u_delta[i] = new Matrix(3);
		d_alpha_delta[i] = new Matrix(3);
		d_u_delta[i] = new Matrix(3);
		A_delta[i] = new Matrix(3, 3);
		Q_delta[i] = new Matrix(3, 3);
		Xi_delta[i] = new Matrix(3, 3);
		d_A_delta[i] = new Matrix(3, 3);
		d_Xi_delta[i] = new Matrix(3, 3);
		d_z[i] = new Matrix(3);
		d_Z[i] = new Matrix(3, 3);
		eta_r[i] = new Matrix(3);
		kappa_r[i] = new Matrix(3);
		epsilon_r[i] = new Matrix(6);
		sigma_r[i] = new Matrix(6);
		n_r[i] = new Matrix(3);
		m_r[i] = new Matrix(3);
		n[i] = new Matrix(3);
		m[i] = new Matrix(3);
	}
	D = new Matrix(6, 6);
	I3 = new Matrix(3, 3);
	e3r = new Matrix(3);
	B1 = new Matrix(6, 6);
	Qtransp = new Matrix(3, 3);
	B2 = new Matrix(6, 9);
	B2temp = new Matrix(3, 3);
	stiffness = new Matrix(18, 18);
	mass = new Matrix(18, 18);
	damping = new Matrix(18, 18);
	mass_modal = new Matrix(18, 18);
	damping_modal = new Matrix(18, 18);
	rayleigh_damping = new Matrix(18, 18);
	inertial_loading = new Matrix(18);
	i_loading = new Matrix(18);
	e_loading = new Matrix(18);
	P_loading = new Matrix(18);
	damping_loading = new Matrix(18);
	constitutive_stiffness = new Matrix(18, 18);
	loading_stiffness = new Matrix(18, 18);
	geometric_stiffness = new Matrix(18, 18);
	transform = new Matrix(18, 18);
	transform3 = new Matrix(3, 3);
	V_alpha_dz_n = new Matrix(3, 3);
	V_alpha_m = new Matrix(3, 3);
	d_V_dalpha_apha_m = new Matrix(3, 3);
	G_d_u_alpha = new Matrix(3, 3);
	G_d_u_alpha_transp = new Matrix(3, 3);
	G_alpha_alpha = new Matrix(3, 3);
	G_alpha_d_alpha = new Matrix(3, 3);
	G_alpha_d_alpha_transp = new Matrix(3, 3);
	B = new Matrix(6, 9);
	G = new Matrix(9, 9);
	lag_save = new LagrangeSave();
	//Variáveis internas para o esforço de correnteza marítima
	e3ip = new Matrix*[2];
	zi = new Matrix*[2];
	vel = new Matrix*[2];
	velr = new Matrix*[2];
	element_vel = new Matrix*[2];
	ut = new Matrix*[2];
	un = new Matrix*[2];
	d_e3_d_alpha = new Matrix*[2];
	Lt = new Matrix*[2];
	Ln = new Matrix*[2];
	L_u_alpha = new Matrix*[2];
	f_current = new Matrix*[2];
	L = new Matrix*[2];
	t_e = new Matrix*[2];
	n_e = new Matrix*[2];
	vtr = new Matrix*[2];
	vnr = new Matrix*[2];
	Mdt = new Matrix*[2];
	Mdn = new Matrix*[2];
	Md2 = new Matrix*[2];
	Un_ = 0;
	un_ = 0;
	Ut_ = 0;
	ut_ = 0;
	C1t = 0;
	C1n = 0;
	//Loop nos pontos de Gauss
	for (int i = 0; i < 2; i++)
	{
		e3ip[i] = new Matrix(3,1);
		zi[i] = new Matrix(3, 1);
		vel[i] = new Matrix(3, 1);
		velr[i] = new Matrix(3, 1);
		element_vel[i] = new Matrix(3, 1);
		ut[i] = new Matrix(3, 1);
		un[i] = new Matrix(3, 1);
		d_e3_d_alpha[i] = new Matrix(3, 3);
		Lt[i] = new Matrix(3, 3);
		Ln[i] = new Matrix(3, 3);
		L_u_alpha[i] = new Matrix(3, 3);
		f_current[i] = new Matrix(3, 1);
		L[i] = new Matrix(6, 6);
		t_e[i] = new Matrix(3, 1);
		n_e[i] = new Matrix(3, 1);
		vtr[i] = new Matrix(3, 1);
		vnr[i] = new Matrix(3, 1);
		Mdt[i] = new Matrix(3, 3);
		Mdn[i] = new Matrix(3, 3);
		Md2[i] = new Matrix(6, 6);
	}
	e3rg = new Matrix(3,1);
	signal_t = 0;
	Cdt = 0;
	Cdn = 0;
	Aext = 0;
	rho_f = 0;
	//Variáveis internas para o esforço pipe load
	kip = new Matrix*[2];
	temp_f = new Matrix*[2];
	temp_m = new Matrix*[2];
	temp_l = new Matrix*[2];
	Kip = new Matrix*[2];
	E3ip = new Matrix*[2];
	UpsilonN = new Matrix*[2];
	//Loop nos pontos de Gauss
	for (int i = 0; i < 2; i++)
	{
		kip[i] = new Matrix(3,1);
		temp_f[i] = new Matrix(3, 1);
		temp_m[i] = new Matrix(3, 1);
		temp_l[i] = new Matrix(6, 1);
		Kip[i] = new Matrix(3, 3);
		E3ip[i] = new Matrix(3, 3);
		UpsilonN[i] = new Matrix(12, 18);
	}
	O1 = new Matrix(3, 3);
	K1ua = new Matrix(3, 3);
	K1aa = new Matrix(3, 3);
	K2ua = new Matrix(3, 3);
	K2au = new Matrix(3, 3);
	Kext = new Matrix(6, 12);
	p0i = 0;
	p0e = 0;
	rhoi = 0;
	rhoe = 0;
	Aint = 0;
	//Variáveis internas para uso na dinâmica
	alpha_dot = new Matrix*[2];
	Xi_dot = new Matrix*[2];
	Mip = new Matrix*[2];
	Jip = new Matrix*[2];
	M = new Matrix*[2];
	Md1 = new Matrix*[2];
	//Loop nos pontos de Gauss
	for (int i = 0; i < 2; i++)
	{
		alpha_dot[i] = new Matrix(3, 1);
		Xi_dot[i] = new Matrix(3, 1);
		Mip[i] = new Matrix(3, 3);
		Jip[i] = new Matrix(3, 3);
		M[i] = new Matrix(6, 6);
		Md1[i] = new Matrix(6, 6);
	}
	Mr = new Matrix(3, 3);
	Jr = new Matrix(3, 3);

	//Ponteiros double** - conversão de matriz
	pJr = new double*[3];
	for (int i = 0; i < 3; i++)
		pJr[i] = new double[3];
	pMr = new double*[3];
	for (int i = 0; i < 3; i++)
		pMr[i] = new double[3];
	pDdT = new double*[6];
	for (int i = 0; i < 6; i++)
		pDdT[i] = new double[6];
	br = new Matrix(3);
	DdT = new Matrix(6,6);
	dT = new Matrix(6);

	morison_loading = new Matrix(18);
	pL = new double*[6];
	for (int i = 0; i < 6; i++)
		pL[i] = new double[6];
	

}

Pipe_1::~Pipe_1()
{
	delete[] nodes;
	delete[] VTK_nodes;
	if (DOFs != NULL)
	{
		for (int i = 0; i < n_nodes; i++)
			delete[] DOFs[i];
		delete[] DOFs;
	}
	delete[] N1;
	delete[] N2;
	delete[] N3;
	delete[] dN1;
	delete[] dN2;
	delete[] dN3;
	delete[] csi;

	delete[] type_name;
	//Loop nos pontos de Gauss
	for (int i = 0; i < 2; i++)
	{
		delete deltaN[i];
		delete N[i];
		delete alpha_delta[i];
		delete u_delta[i];
		delete d_alpha_delta[i];
		delete d_u_delta[i];
		delete A_delta[i];
		delete Q_delta[i];
		delete Xi_delta[i];
		delete d_A_delta[i];
		delete d_Xi_delta[i];
		delete d_z[i];
		delete d_Z[i];
		delete eta_r[i];
		delete kappa_r[i];
		delete epsilon_r[i];
		delete sigma_r[i];
		delete n_r[i];
		delete m_r[i];
		delete n[i];
		delete m[i];
	}
	delete[] deltaN;
	delete[] N;
	delete[] alpha_delta;
	delete[] u_delta;
	delete[] d_alpha_delta;
	delete[] d_u_delta;
	delete[] A_delta;
	delete[] Q_delta;
	delete[] Xi_delta;
	delete[] d_A_delta;
	delete[] d_Xi_delta;
	delete[] d_z;
	delete[] d_Z;
	delete[] eta_r;
	delete[] kappa_r;
	delete[] epsilon_r;
	delete[] sigma_r;
	delete[] n_r;
	delete[] m_r;
	delete[] n;
	delete[] m;
	delete D;
	delete I3;
	delete e3r;
	delete B1;
	delete Qtransp;
	delete B2;
	delete B2temp;
	delete stiffness;
	delete mass;
	delete damping;
	delete mass_modal;
	delete damping_modal;
	delete rayleigh_damping;
	delete e_loading;
	delete i_loading;
	delete inertial_loading;
	delete P_loading;
	delete damping_loading;
	delete constitutive_stiffness;
	delete loading_stiffness;
	delete geometric_stiffness;
	delete transform;
	delete transform3;
	delete V_alpha_dz_n;
	delete V_alpha_m;
	delete d_V_dalpha_apha_m;
	delete G_d_u_alpha;
	delete G_d_u_alpha_transp;
	delete G_alpha_alpha;
	delete G_alpha_d_alpha;
	delete G_alpha_d_alpha_transp;
	delete B;
	delete G;
	delete lag_save;
	//Variáveis internas para o esforço de correnteza marítima
	//Loop nos pontos de Gauss
	for (int i = 0; i < 2; i++)
	{
		delete e3ip[i];
		delete zi[i];
		delete vel[i];
		delete velr[i];
		delete element_vel[i];
		delete ut[i];
		delete un[i];
		delete d_e3_d_alpha[i];
		delete Lt[i];
		delete Ln[i];
		delete L_u_alpha[i];
		delete f_current[i];
		delete L[i];
		delete t_e[i];
		delete n_e[i];
		delete vtr[i];
		delete vnr[i];
		delete Mdt[i];
		delete Mdn[i];
		delete Md2[i];
	}
	delete[] e3ip;
	delete[] zi;
	delete[] vel;
	delete[] velr;
	delete[] element_vel;
	delete[] ut;
	delete[] un;
	delete[] d_e3_d_alpha;
	delete[] Lt;
	delete[] Ln;
	delete[] L_u_alpha;
	delete[] f_current;
	delete[] L;
	delete[] t_e;
	delete[] n_e;
	delete[] vtr;
	delete[] vnr;
	delete[] Mdt;
	delete[] Mdn;
	delete[] Md2;
	delete e3rg;
	//Variáveis internas para o esforço pipe load
	//Loop nos pontos de Gauss
	for (int i = 0; i < 2; i++)
	{
		delete kip[i];
		delete temp_f[i];
		delete temp_m[i];
		delete temp_l[i];
		delete Kip[i];
		delete E3ip[i];
		delete UpsilonN[i];
	}
	delete[] kip;
	delete[] temp_f;
	delete[] temp_m;
	delete[] temp_l;
	delete[] Kip;
	delete[] E3ip;
	delete[] UpsilonN;
	delete O1;
	delete K1ua;
	delete K1aa;
	delete K2ua;
	delete K2au;
	delete Kext;
	//Variáveis internas para uso na dinâmica
	//Loop nos pontos de Gauss
	for (int i = 0; i < 2; i++)
	{
		delete alpha_dot[i];
		delete Xi_dot[i];


		delete Mip[i];
		delete Jip[i];
		delete M[i];
		delete Md1[i];
	}
	delete[] alpha_dot;
	delete[] Xi_dot;
	delete Mr;
	delete Jr;
	delete[] Mip;
	delete[] Jip;
	delete[] M;
	delete[] Md1;

	delete br;
	delete DdT;
	delete dT;

	//Ponteiros de double**
	for (int i = 0; i < 3; i++)
		delete[]pMr[i];
	delete[] pMr;
	for (int i = 0; i < 3; i++)
		delete[]pJr[i];
	delete[] pJr;
	for (int i = 0; i < 6; i++)
		delete[]pDdT[i];
	delete[] pDdT;

	delete morison_loading;
	for (int i = 0; i < 6; i++)
		delete[]pL[i];
	delete[] pL;
	
}

//Checa inconsistências no elemento para evitar erros de execução
bool Pipe_1::Check()
{
	for (int i = 0; i < n_nodes; i++)
	{
		if (nodes[i] > db.number_nodes)
			return false;
	}
	if (cs > db.number_CS)
		return false;
	if (section > db.number_pipe_sections)
		return false;
	//Checagem do alinhamento dos eixos do CS com direção dos nós
	Matrix e1r = *db.CS[cs - 1]->E1;
	Matrix e2r = *db.CS[cs - 1]->E2;
	Matrix e3r = *db.CS[cs - 1]->E3;
	double tolortho = 1e-8;
	Matrix e3_check(3);
	//Versor e3_check (configuração de referência)
	e3_check(0, 0) = db.nodes[nodes[2] - 1]->ref_coordinates[0] - db.nodes[nodes[0] - 1]->ref_coordinates[0];
	e3_check(1, 0) = db.nodes[nodes[2] - 1]->ref_coordinates[1] - db.nodes[nodes[0] - 1]->ref_coordinates[1];
	e3_check(2, 0) = db.nodes[nodes[2] - 1]->ref_coordinates[2] - db.nodes[nodes[0] - 1]->ref_coordinates[2];
	e3_check = (1.0 / (norm(e3_check)))*(e3_check);
	if (abs(dot(e3_check, e1r)) > tolortho || abs(dot(e3_check, e2r)) > tolortho)
	{
		db.myprintf("Check CS alignment and nodes definitions in element %d.\n", this->number);
		return false;
	}
	if (dot(e3_check, e3r) < 0.0)
	{
		db.myprintf("Check consistency between the sequence of nodes defined for element %d and direction E3 of the chosen CS (number %d).\n", this->number, cs);
		return false;
	}
	return true;
}

bool Pipe_1::Read(FILE *f)
{
	char s[1000];
	fscanf(f, "%s", s);
	number = atoi(s);

	fscanf(f, "%s", s);
	if (!strcmp(s, "PipeSec"))
	{
		fscanf(f, "%s", s);
		section = atoi(s);
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
		fscanf(f, "%s", s);
		nodes[2] = atoi(s);
	}
	else
		return false;
	return true;
}

void Pipe_1::Write(FILE *f)
{
	fprintf(f, "Pipe_1\t%d\tPipeSec\t%d\tCS\t%d\tNodes\t%d\t%d\t%d\n",
		number,
		section,
		cs,
		nodes[0],
		nodes[1],
		nodes[2]);
}

//Escreve arquivo de resultados
void Pipe_1::WriteResults(FILE *f)
{
	fprintf(f, "Pipe_1\t%d\tStress_r\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\tStrain_r\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\tRef_Coor\t%.6e\t%.6e\t%.6e\tCoor\t%.6e\t%.6e\t%.6e\n",
		number,
		((*sigma_r[0])(0, 0) + (*sigma_r[1])(0, 0)) / 2, ((*sigma_r[0])(1, 0) + (*sigma_r[1])(1, 0)) / 2, ((*sigma_r[0])(2, 0) + (*sigma_r[1])(2, 0)) / 2,
		((*sigma_r[0])(3, 0) + (*sigma_r[1])(3, 0)) / 2, ((*sigma_r[0])(4, 0) + (*sigma_r[1])(4, 0)) / 2, ((*sigma_r[0])(5, 0) + (*sigma_r[1])(5, 0)) / 2,
		((*epsilon_r[0])(0, 0) + (*epsilon_r[1])(0, 0)) / 2, ((*epsilon_r[0])(1, 0) + (*epsilon_r[1])(1, 0)) / 2, ((*epsilon_r[0])(2, 0) + (*epsilon_r[1])(2, 0)) / 2,
		((*epsilon_r[0])(3, 0) + (*epsilon_r[1])(3, 0)) / 2, ((*epsilon_r[0])(4, 0) + (*epsilon_r[1])(4, 0)) / 2, ((*epsilon_r[0])(5, 0) + (*epsilon_r[1])(5, 0)) / 2,
		db.nodes[nodes[1] - 1]->ref_coordinates[0], db.nodes[nodes[1] - 1]->ref_coordinates[1], db.nodes[nodes[1] - 1]->ref_coordinates[2],
		db.nodes[nodes[1] - 1]->copy_coordinates[0] + db.nodes[nodes[1] - 1]->displacements[0], 
		db.nodes[nodes[1] - 1]->copy_coordinates[1] + db.nodes[nodes[1] - 1]->displacements[1], 
		db.nodes[nodes[1] - 1]->copy_coordinates[2] + db.nodes[nodes[1] - 1]->displacements[2]
		);
}

void Pipe_1::WriteVTK_XMLBase(std::vector<float> *float_vector)
{
	//Imprime os resultados do elemento
	int res_element = 6;
	float_vector->push_back((float)((*sigma_r[0])(0, 0) + (*sigma_r[1])(0, 0)) / 2);
	float_vector->push_back((float)((*sigma_r[0])(1, 0) + (*sigma_r[1])(1, 0)) / 2);
	float_vector->push_back((float)((*sigma_r[0])(2, 0) + (*sigma_r[1])(2, 0)) / 2);
	float_vector->push_back((float)((*sigma_r[0])(3, 0) + (*sigma_r[1])(3, 0)) / 2);
	float_vector->push_back((float)((*sigma_r[0])(4, 0) + (*sigma_r[1])(4, 0)) / 2);
	float_vector->push_back((float)((*sigma_r[0])(5, 0) + (*sigma_r[1])(5, 0)) / 2);
	//Imprime valores nulos para resultados que não fazem sentido para esse tipo de elemento
	for (int i = res_element; i < db.n_element_results; i++)
		float_vector->push_back(0.0);
}

void Pipe_1::WriteVTK_XMLRender(FILE *f)
{
	//vetores para escrita no formato binário - usando a função 'enconde'
	std::vector<float> float_vector;
	std::vector<int> int_vector;

	Matrix Q(3, 3);
	Matrix alpha_1(3);
	double alpha_escalar;
	Matrix A;
	Matrix vec_P(3);
	double factor = db.pipe_sections[section - 1]->Di / db.pipe_sections[section - 1]->De;
	//Número de pontos a serem gerados
	int n_points = n_nodes*db.pipe_sections[section - 1]->sec_details->n_points * 2;
	//Número de células a serem geradas
	int n_cells = db.pipe_sections[section - 1]->sec_details->n_points*(n_nodes - 1);
	//Opens Piece
	fprintf(f, "\t\t<Piece NumberOfPoints = \"%d\" NumberOfCells = \"%d\">\n", n_points, n_cells);
		//Opens Points
		fprintf(f, "\t\t\t<Points>\n");
			//Opens DataArray
			fprintf(f, "\t\t\t\t<DataArray type = \"Float32\" NumberOfComponents = \"3\" format=\"binary\">\n");
			float_vector.clear();
			//Preenchendo as coordenadas dos pontos
			for (int i = 0; i < n_nodes; i++)//percorrendo os 3 nós do elemento
			{
				//Para cada nó do elemento, calcula o tensor rotação
				alpha_1(0, 0) = db.nodes[nodes[i] - 1]->copy_coordinates[3];
				alpha_1(1, 0) = db.nodes[nodes[i] - 1]->copy_coordinates[4];
				alpha_1(2, 0) = db.nodes[nodes[i] - 1]->copy_coordinates[5];
				alpha_escalar = norm(alpha_1);
				A = skew(alpha_1);
				g = 4.0 / (4.0 + alpha_escalar*alpha_escalar);
				Q = *I3 + g*(A + 0.5*A*A);
				Q = Q*transp(*transform3);//Matriz de transformação para trazer o vetor da ST do plano xy para o plano atual em que ela se encontra

				for (int point = 0; point < db.pipe_sections[section - 1]->sec_details->n_points; point++)//Percorre os nós que descrevem o perímetro da ST - diâmetro externo
				{
					//Posição de cada ponto P no plano xy (referência)
					vec_P(0, 0) = db.pipe_sections[section - 1]->sec_details->points[point][0];
					vec_P(1, 0) = db.pipe_sections[section - 1]->sec_details->points[point][1];
					vec_P(2, 0) = 0.0;
					vec_P = Q*vec_P;//Operando rotacionando para o sistema da barra
					for (int c = 0; c < 3; c++)
						vec_P(c, 0) += db.nodes[nodes[i] - 1]->ref_coordinates[c] + db.post_files->mag_factor*(db.nodes[nodes[i] - 1]->copy_coordinates[c] - db.nodes[nodes[i] - 1]->ref_coordinates[c]);//Translação - soma a posição do centro da barra (do ponto em questão)
					float_vector.push_back((float)(vec_P(0, 0)));
					float_vector.push_back((float)(vec_P(1, 0)));
					float_vector.push_back((float)(vec_P(2, 0)));
				}
				for (int point = 0; point < db.pipe_sections[section - 1]->sec_details->n_points; point++)//Percorre os nós que descrevem o perímetro da ST - diâmetro interno
				{
					//Posição de cada ponto P no plano xy (referência)
					vec_P(0, 0) = db.pipe_sections[section - 1]->sec_details->points[point][0] * factor;
					vec_P(1, 0) = db.pipe_sections[section - 1]->sec_details->points[point][1] * factor;
					vec_P(2, 0) = 0.0;
					vec_P = Q*vec_P;//Operando rotacionando para o sistema da barra
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

		int nodes[8];
		int offset = db.pipe_sections[section - 1]->sec_details->n_points * 2;
		//Opens Cells
		fprintf(f, "\t\t\t<Cells>\n");
		//Opens DataArray
		fprintf(f, "\t\t\t\t<DataArray type = \"Int32\" Name = \"connectivity\" format =\"binary\">\n");
		int_vector.clear();
		for (int cell = 0; cell < db.pipe_sections[section - 1]->sec_details->n_points; cell++)
		{
			//Se não for a última célula ao longo do perímetro
			if (cell != db.pipe_sections[section - 1]->sec_details->n_points-1)
			{
				nodes[0] = cell;
				nodes[1] = cell + 1;
				nodes[2] = cell + db.pipe_sections[section - 1]->sec_details->n_points + 1;
				nodes[3] = cell + db.pipe_sections[section - 1]->sec_details->n_points;
				nodes[4] = cell + offset;
				nodes[5] = cell + 1 + offset;
				nodes[6] = cell + db.pipe_sections[section - 1]->sec_details->n_points + 1 + offset;
				nodes[7] = cell + db.pipe_sections[section - 1]->sec_details->n_points + offset;
			}
			else
			{
				nodes[0] = cell;
				nodes[1] = 0;
				nodes[2] = db.pipe_sections[section - 1]->sec_details->n_points;
				nodes[3] = cell + db.pipe_sections[section - 1]->sec_details->n_points;
				nodes[4] = cell + offset;
				nodes[5] = 0 + offset;
				nodes[6] = db.pipe_sections[section - 1]->sec_details->n_points + offset;
				nodes[7] = cell + db.pipe_sections[section - 1]->sec_details->n_points + offset;
			}
			int_vector.push_back(nodes[0]);
			int_vector.push_back(nodes[1]);
			int_vector.push_back(nodes[2]);
			int_vector.push_back(nodes[3]);
			int_vector.push_back(nodes[4]);
			int_vector.push_back(nodes[5]);
			int_vector.push_back(nodes[6]);
			int_vector.push_back(nodes[7]);
			//fprintf(f, "\t\t\t\t\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", nodes[0], nodes[1], nodes[2], nodes[3], nodes[4], nodes[5], nodes[6], nodes[7]);

			//Incrementa todos os pontos de offset para fazer a segunda metade dos polígonos
			for (int i = 0; i < 8; i++)
				nodes[i] += offset;
			//fprintf(f, "\t\t\t\t\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", nodes[0], nodes[1], nodes[2], nodes[3], nodes[4], nodes[5], nodes[6], nodes[7]);
			int_vector.push_back(nodes[0]);
			int_vector.push_back(nodes[1]);
			int_vector.push_back(nodes[2]);
			int_vector.push_back(nodes[3]);
			int_vector.push_back(nodes[4]);
			int_vector.push_back(nodes[5]);
			int_vector.push_back(nodes[6]);
			int_vector.push_back(nodes[7]);
		}
		fprintf(f, encodeData(int_vector).c_str());
		fprintf(f, "\n");
		//Closes DataArray
		fprintf(f, "\t\t\t\t</DataArray>\n");
		//Opens DataArray
		fprintf(f, "\t\t\t\t<DataArray type = \"Int32\" Name = \"types\" format=\"binary\">\n");
		int_vector.clear();
		for (int cell = 0; cell < n_cells; cell++)
			int_vector.push_back(12);
			//fprintf(f, "\t\t\t\t\t%d\n", 12);
		fprintf(f, encodeData(int_vector).c_str());
		fprintf(f, "\n");
		//Closes DataArray
		fprintf(f, "\t\t\t\t</DataArray>\n");
		//Opens DataArray
		fprintf(f, "\t\t\t\t<DataArray type = \"Int32\" Name = \"offsets\" format=\"binary\">\n");
		int_vector.clear();
		int cur_off = 0;
		for (int cell = 0; cell < n_cells; cell++)
		{
			cur_off += 8;
			int_vector.push_back(cur_off);
			//fprintf(f, "\t\t\t\t\t%d\n", cur_off);
		}
		fprintf(f, encodeData(int_vector).c_str());
		fprintf(f, "\n");
		//Closes DataArray
		fprintf(f, "\t\t\t\t</DataArray>\n");
		//Closes Cells
		fprintf(f, "\t\t\t</Cells>\n");

		//Opens CellData
		fprintf(f, "\t\t\t<CellData FieldData = \"ElementData\">\n");
		float_vector.clear();
		//Opens DataArray
		fprintf(f, "\t\t\t\t<DataArray Name = \"ElementResults\" type = \"Float32\" NumberOfComponents=\"%d\" format=\"binary\">\n", db.n_element_results);
		for (int cell = 0; cell < n_cells; cell++)
		{
			//Imprime os resultados do elemento
			int res_element = 6;
			float_vector.push_back((float)((*sigma_r[0])(0, 0) + (*sigma_r[1])(0, 0)) / 2);
			float_vector.push_back((float)((*sigma_r[0])(1, 0) + (*sigma_r[1])(1, 0)) / 2);
			float_vector.push_back((float)((*sigma_r[0])(2, 0) + (*sigma_r[1])(2, 0)) / 2);
			float_vector.push_back((float)((*sigma_r[0])(3, 0) + (*sigma_r[1])(3, 0)) / 2);
			float_vector.push_back((float)((*sigma_r[0])(4, 0) + (*sigma_r[1])(4, 0)) / 2);
			float_vector.push_back((float)((*sigma_r[0])(5, 0) + (*sigma_r[1])(5, 0)) / 2);
			/*fprintf(f, "\t\t\t\t\t%f\t%f\t%f\t%f\t%f\t%f\t",
				((*sigma_r[0])(0, 0) + (*sigma_r[1])(0, 0)) / 2, ((*sigma_r[0])(1, 0) + (*sigma_r[1])(1, 0)) / 2, ((*sigma_r[0])(2, 0) + (*sigma_r[1])(2, 0)) / 2,
				((*sigma_r[0])(3, 0) + (*sigma_r[1])(3, 0)) / 2, ((*sigma_r[0])(4, 0) + (*sigma_r[1])(4, 0)) / 2, ((*sigma_r[0])(5, 0) + (*sigma_r[1])(5, 0)) / 2
				);*/
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
			int_vector.push_back(2);		//Element ID
			int_vector.push_back(0);
			int_vector.push_back(section);
			int_vector.push_back(cs);
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

//Escreve no monitor do elemento//Escreve no monitor do elemento
void Pipe_1::WriteMonitor(FILE *f, bool first_record, double time)
{
	//Cabeçalho
	if (first_record == true)
		fprintf(f, "TIME\tnx1r\tny1r\tnz1r\tmx1r\tmy1r\tmz1r\tnx2r\tny2r\tnz2r\tmx2r\tmy2r\tmz2r\n");
	//Informações a serem salvas
	fprintf(f, "%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\n",
		time,
		(*sigma_r[0])(0, 0), (*sigma_r[0])(1, 0), (*sigma_r[0])(2, 0), (*sigma_r[0])(3, 0), (*sigma_r[0])(4, 0), (*sigma_r[0])(5, 0),
		(*sigma_r[1])(0, 0), (*sigma_r[1])(1, 0), (*sigma_r[1])(2, 0), (*sigma_r[1])(3, 0), (*sigma_r[1])(4, 0), (*sigma_r[1])(5, 0)
		);
}

//Monta elementos
void Pipe_1::Mount()
{
	Zeros();	//Zera matrizes
	//Loop nos pontos de Gauss
	for (int gauss = 0; gauss < 2; gauss++)
	{
		//alpha_delta calculado no ponto de Gauss							
		(*alpha_delta[gauss])(0, 0) =	(db.nodes[nodes[0] - 1]->displacements[3])* N1[gauss] +
										(db.nodes[nodes[1] - 1]->displacements[3])* N2[gauss] +
										(db.nodes[nodes[2] - 1]->displacements[3])* N3[gauss] ;
		(*alpha_delta[gauss])(1, 0) =	(db.nodes[nodes[0] - 1]->displacements[4])* N1[gauss] +
										(db.nodes[nodes[1] - 1]->displacements[4])* N2[gauss] +
										(db.nodes[nodes[2] - 1]->displacements[4])* N3[gauss] ;
		(*alpha_delta[gauss])(2, 0) =	(db.nodes[nodes[0] - 1]->displacements[5])* N1[gauss] +
										(db.nodes[nodes[1] - 1]->displacements[5])* N2[gauss] +
										(db.nodes[nodes[2] - 1]->displacements[5])* N3[gauss] ;
		//d_alpha_delta calculado no ponto de Gauss							
		(*d_alpha_delta[gauss])(0, 0) = (db.nodes[nodes[0] - 1]->displacements[3])* dN1[gauss] +
										(db.nodes[nodes[1] - 1]->displacements[3])* dN2[gauss] +
										(db.nodes[nodes[2] - 1]->displacements[3])* dN3[gauss] ;
		(*d_alpha_delta[gauss])(1, 0) = (db.nodes[nodes[0] - 1]->displacements[4])* dN1[gauss] +
										(db.nodes[nodes[1] - 1]->displacements[4])* dN2[gauss] +
										(db.nodes[nodes[2] - 1]->displacements[4])* dN3[gauss] ;
		(*d_alpha_delta[gauss])(2, 0) = (db.nodes[nodes[0] - 1]->displacements[5])* dN1[gauss] +
										(db.nodes[nodes[1] - 1]->displacements[5])* dN2[gauss] +
										(db.nodes[nodes[2] - 1]->displacements[5])* dN3[gauss] ;
		//u_delta calculado no ponto de Gauss	
		(*u_delta[gauss])(0, 0) =	(db.nodes[nodes[0] - 1]->displacements[0])* N1[gauss] +
									(db.nodes[nodes[1] - 1]->displacements[0])* N2[gauss] +
									(db.nodes[nodes[2] - 1]->displacements[0])* N3[gauss] ;
		(*u_delta[gauss])(1, 0) =	(db.nodes[nodes[0] - 1]->displacements[1])* N1[gauss] +
									(db.nodes[nodes[1] - 1]->displacements[1])* N2[gauss] +
									(db.nodes[nodes[2] - 1]->displacements[1])* N3[gauss] ;
		(*u_delta[gauss])(2, 0) =	(db.nodes[nodes[0] - 1]->displacements[2])* N1[gauss] +
									(db.nodes[nodes[1] - 1]->displacements[2])* N2[gauss] +
									(db.nodes[nodes[2] - 1]->displacements[2])* N3[gauss] ;
		//d_u_delta calculado no ponto de Gauss	
		(*d_u_delta[gauss])(0, 0) =	(db.nodes[nodes[0] - 1]->displacements[0])* dN1[gauss] +
									(db.nodes[nodes[1] - 1]->displacements[0])* dN2[gauss] +
									(db.nodes[nodes[2] - 1]->displacements[0])* dN3[gauss] ;
		(*d_u_delta[gauss])(1, 0) = (db.nodes[nodes[0] - 1]->displacements[1])* dN1[gauss] +
									(db.nodes[nodes[1] - 1]->displacements[1])* dN2[gauss] +
									(db.nodes[nodes[2] - 1]->displacements[1])* dN3[gauss] ;
		(*d_u_delta[gauss])(2, 0) = (db.nodes[nodes[0] - 1]->displacements[2])* dN1[gauss] +
									(db.nodes[nodes[1] - 1]->displacements[2])* dN2[gauss] +
									(db.nodes[nodes[2] - 1]->displacements[2])* dN3[gauss] ;
		//Transformação para o sistema local do elemento
		*alpha_delta[gauss] = (*transform3)*(*alpha_delta[gauss]);
		*d_alpha_delta[gauss] = (*transform3)*(*d_alpha_delta[gauss]);
		*u_delta[gauss] = (*transform3)*(*u_delta[gauss]);
		*d_u_delta[gauss] = (*transform3)*(*d_u_delta[gauss]);
		//Cálculo de tensores de rotação e outros
		alpha_escalar_delta = norm(*alpha_delta[gauss]);																	//Valor escalar do parametro alpha
		*A_delta[gauss] = skew(*alpha_delta[gauss]);																		//Matriz A_delta
		g = 4.0 / (4.0 + alpha_escalar_delta*alpha_escalar_delta);															//função g(alpha) - em algumas ref. tb. chamado de h(alpha)
		*Q_delta[gauss] = *I3 + g*(*A_delta[gauss] + 0.5*((*A_delta[gauss])*(*A_delta[gauss])));							//Tensor de rotação
		*Xi_delta[gauss] = g*(*I3 + 0.5*(*A_delta[gauss]));
		*d_A_delta[gauss] = skew(*d_alpha_delta[gauss]);																	//Derivada do vetor d_alpha
		*d_Xi_delta[gauss] = -0.5*g*((dot(*alpha_delta[gauss], *d_alpha_delta[gauss]))*(*Xi_delta[gauss]) - (*d_A_delta[gauss]));		//Derivada do tensor Xi
		*d_z[gauss] = *d_u_delta[gauss] + (*lag_save->dz_i[gauss]);															//Vetor derivada de z (calculada em i+1)
		*d_Z[gauss] = skew(*d_z[gauss]);																					//Tensor anti-simétrico, cuja axial é z'
		*Qtransp = transp((*Q_delta[gauss])*(*lag_save->Q_i[gauss]));	//Q transposta (i+1)
		for (int i = 0; i<3; i++)
			for (int j = 0; j<3; j++)
				(*B1)(i, j) = (*Qtransp)(i, j);
		for (int i = 0; i<3; i++)
			for (int j = 0; j<3; j++)
				(*B1)(i + 3, j + 3) = (*Qtransp)(i, j);			//Matriz B1 preenchida
		*B2temp = (*d_Z[gauss])*(*Xi_delta[gauss]);
		for (int i = 0; i<3; i++)
			for (int j = 0; j<3; j++)
				(*B2)(i, j + 6) = (*B2temp)(i, j);
		for (int i = 0; i<3; i++)
			for (int j = 0; j<3; j++)
				(*B2)(i + 3, j + 3) = (*Xi_delta[gauss])(i, j);
		for (int i = 0; i<3; i++)
			for (int j = 0; j<3; j++)
				(*B2)(i + 3, j + 6) = (*d_Xi_delta[gauss])(i, j);
		(*B) = (*B1)*(*B2);
		//Matriz de rigidez constitutiva
		(*constitutive_stiffness) = (transp(*deltaN[gauss])) * ((((transp(*B)) * (*D)) * (*B)) * (*deltaN[gauss]));
		//Deformações retro-rotacionadas
		*eta_r[gauss] = transp((*Q_delta[gauss]) * (*lag_save->Q_i[gauss]))*(*d_z[gauss]) - *e3r;
		*kappa_r[gauss] = transp(*lag_save->Q_i[gauss])*transp(*Xi_delta[gauss])*(*d_alpha_delta[gauss]) + (*lag_save->kappa_i_ref[gauss]);
		//Vetor deformação generalizada retro-rotacionada
		(*epsilon_r[gauss])(0, 0) = (*eta_r[gauss])(0, 0);
		(*epsilon_r[gauss])(1, 0) = (*eta_r[gauss])(1, 0);
		(*epsilon_r[gauss])(2, 0) = (*eta_r[gauss])(2, 0);
		(*epsilon_r[gauss])(3, 0) = (*kappa_r[gauss])(0, 0);
		(*epsilon_r[gauss])(4, 0) = (*kappa_r[gauss])(1, 0);
		(*epsilon_r[gauss])(5, 0) = (*kappa_r[gauss])(2, 0);
		//Forças e momentos internos retro-rotacionados
		(*sigma_r[gauss]) = (*D)*(*epsilon_r[gauss]);
		(*n_r[gauss])(0, 0) = (*sigma_r[gauss])(0, 0);
		(*n_r[gauss])(1, 0) = (*sigma_r[gauss])(1, 0);
		(*n_r[gauss])(2, 0) = (*sigma_r[gauss])(2, 0);
		(*m_r[gauss])(0, 0) = (*sigma_r[gauss])(3, 0);
		(*m_r[gauss])(1, 0) = (*sigma_r[gauss])(4, 0);
		(*m_r[gauss])(2, 0) = (*sigma_r[gauss])(5, 0);
		//Forças e momentos internos não retro-rotacionados
		(*n[gauss]) = ((*Q_delta[gauss])* (*lag_save->Q_i[gauss]))*(*n_r[gauss]);
		(*m[gauss]) = ((*Q_delta[gauss])* (*lag_save->Q_i[gauss]))*(*m_r[gauss]);
		//Cálculo da matriz G
		*V_alpha_dz_n = V(*alpha_delta[gauss], (*d_Z[gauss])* (*n[gauss]), alpha_escalar_delta);
		*V_alpha_m = V(*alpha_delta[gauss], (*m[gauss]), alpha_escalar_delta);
		*d_V_dalpha_apha_m = d_V(*alpha_delta[gauss], *d_alpha_delta[gauss], *m[gauss], alpha_escalar_delta);
		*G_d_u_alpha = -1.0* skew(*n[gauss]) * (*Xi_delta[gauss]);
		*G_d_u_alpha_transp = 1.0*transp(*Xi_delta[gauss])*skew(*n[gauss]);
		*G_alpha_alpha = 1.0*(transp(*Xi_delta[gauss])) * ((*d_Z[gauss])* (skew(*n[gauss])))*(*Xi_delta[gauss]) - 1.0* (*V_alpha_dz_n) + (*d_V_dalpha_apha_m) - (transp(*d_Xi_delta[gauss]))*((skew(*m[gauss]))*(*Xi_delta[gauss]));
		*G_alpha_d_alpha = *V_alpha_m - 1.0*(transp(*Xi_delta[gauss]))*((skew(*m[gauss]))*(*Xi_delta[gauss]));
		*G_alpha_d_alpha_transp = *V_alpha_m;
		//Cálculo da matriz G:
		for (int i = 0; i<3; i++)
			for (int j = 0; j<3; j++)
				(*G)(i, j + 6) = (*G_d_u_alpha)(i, j);
		for (int i = 0; i<3; i++)
			for (int j = 0; j<3; j++)
				(*G)(i + 6, j) = (*G_d_u_alpha_transp)(i, j);
		for (int i = 0; i<3; i++)
			for (int j = 0; j<3; j++)
				(*G)(i + 6, j + 6) = (*G_alpha_alpha)(i, j);
		for (int i = 0; i<3; i++)
			for (int j = 0; j<3; j++)
				(*G)(i + 3, j + 6) = (*G_alpha_d_alpha)(i, j);
		for (int i = 0; i<3; i++)
			for (int j = 0; j<3; j++)
				(*G)(i + 6, j + 3) = (*G_alpha_d_alpha_transp)(i, j);
		//Matriz de rigidez geométrica
		(*geometric_stiffness) = (transp(*deltaN[gauss]))*((*G)*(*deltaN[gauss]));
		//Atualiza a matriz de rigidez do elemento
		(*stiffness) = (*stiffness) + (alpha1*jacobian)*((*constitutive_stiffness) + (*geometric_stiffness));
		//Esforços internos calculados
		(*i_loading) = (*i_loading) + (alpha1*jacobian)*(((transp(*deltaN[gauss])) * (transp(*B))) * (*sigma_r[gauss]));
	}//end of gauss loop
	//Transformações de coordenadas, para o sistema global - matriz de rigidez
	(*stiffness) = (transp(*transform)*(*stiffness))*(*transform);
	(*i_loading) = transp(*transform)*(*i_loading);
}
//Monta carregamentos associados ao elemento
void Pipe_1::MountElementLoads()
{
	MountFieldLoading();		//Inclui atualizações no vetor e_loading
	MountSeaCurrentLoading();	//Inclui atualizações no vetor e_loading e matriz loading_stiffness e stiffness
	(*P_loading) = (*i_loading) - (*e_loading);		//Vetor esforço desbalanceado
}
//Monta matriz de transformação de coordenadas
void Pipe_1::TransformMatrix()
{
	Matrix e1(3);
	e1(0, 0) = 1.0;
	Matrix e2(3);
	e2(1, 0) = 1.0;
	Matrix e3(3);
	e3(2, 0) = 1.0;
	Matrix e1r = *db.CS[cs - 1]->E1;
	Matrix e2r = *db.CS[cs - 1]->E2;
	Matrix e3r = *db.CS[cs - 1]->E3;
	//Preenche a matriz de transformação de coordenadas
	for (int i = 0; i<18; i = i + 3)
	{
		//Preenche também a matriz transform3
		if (i == 0)
		{
			(*transform3)(0, 0) = dot(e1r, e1);
			(*transform3)(0, 1) = dot(e1r, e2);
			(*transform3)(0, 2) = dot(e1r, e3);

			(*transform3)(1, 0) = dot(e2r, e1);
			(*transform3)(1, 1) = dot(e2r, e2);
			(*transform3)(1, 2) = dot(e2r, e3);

			(*transform3)(2, 0) = dot(e3r, e1);
			(*transform3)(2, 1) = dot(e3r, e2);
			(*transform3)(2, 2) = dot(e3r, e3);
		}

		(*transform)(0 + i, 0 + i) = dot(e1r, e1);
		(*transform)(0 + i, 1 + i) = dot(e1r, e2);
		(*transform)(0 + i, 2 + i) = dot(e1r, e3);

		(*transform)(1 + i, 0 + i) = dot(e2r, e1);
		(*transform)(1 + i, 1 + i) = dot(e2r, e2);
		(*transform)(1 + i, 2 + i) = dot(e2r, e3);

		(*transform)(2 + i, 0 + i) = dot(e3r, e1);
		(*transform)(2 + i, 1 + i) = dot(e3r, e2);
		(*transform)(2 + i, 2 + i) = dot(e3r, e3);
	}
}

//Preenche a contribuição do elemento nas matrizes globais
void Pipe_1::MountGlobal()
{
	//Variáveis temporárias para salvar a indexação global dos graus de liberdade a serem setados na matriz de rigidez global
	int GL_global_1 = 0;
	int GL_global_2 = 0;
	double anterior = 0;
	for (int i = 0; i<18; i++)
	{
		//////////////MONTAGEM DO VETOR DE ESFORÇOS DESBALANCEADOS//////////////////
		//Toma do vetor de GL globais, a indexação de cada grau de liberdade global
		if (i<6)
			GL_global_1 = db.nodes[nodes[0] - 1]->GLs[i];
		if (i > 5 && i < 12)
			GL_global_1 = db.nodes[nodes[1] - 1]->GLs[i - 6];
		if (i > 11)
			GL_global_1 = db.nodes[nodes[2] - 1]->GLs[i - 12];

		//Caso o grau de liberdade seja livre:
		if (GL_global_1 > 0)
		{
			anterior = db.global_P_A(GL_global_1 - 1, 0);
			db.global_P_A(GL_global_1 - 1, 0) = anterior + (*P_loading)(i, 0);
			anterior = db.global_I_A(GL_global_1 - 1, 0);
			db.global_I_A(GL_global_1 - 1, 0) = anterior + (*P_loading)(i, 0);
		}
		else
		{
			anterior = db.global_P_B(-GL_global_1 - 1, 0);
			db.global_P_B(-GL_global_1 - 1, 0) = anterior + (*P_loading)(i, 0);
		}
		for (int j = 0; j<18; j++)
		{
			//////////////////////MONTAGEM DA MATRIZ DE RIGIDEZ/////////////////////////
			//Toma do vetor de GL globais, a indexação de cada grau de liberdade global
			if (j<6)
				GL_global_2 = db.nodes[nodes[0] - 1]->GLs[j];
			if (j > 5 && j < 12)
				GL_global_2 = db.nodes[nodes[1] - 1]->GLs[j - 6];
			if (j > 11)
				GL_global_2 = db.nodes[nodes[2] - 1]->GLs[j - 12];

			//Caso os graus de liberdade sejam ambos livres (Matriz Kaa)
			if (GL_global_1 > 0 && GL_global_2 > 0)
			{
				db.global_stiffness_AA.setValue(GL_global_1 - 1, GL_global_2 - 1, (*stiffness)(i, j));
			}
			//Caso os graus de liberdade sejam ambos fixos (Matriz Kbb)
			if (GL_global_1 < 0 && GL_global_2 < 0)
			{
				db.global_stiffness_BB.setValue(-GL_global_1 - 1, -GL_global_2 - 1, (*stiffness)(i, j));
			}
			//Caso os graus de liberdade sejam livre e fixo (Matriz Kab)
			if (GL_global_1 > 0 && GL_global_2 < 0)
			{
				db.global_stiffness_AB.setValue(GL_global_1 - 1, -GL_global_2 - 1, (*stiffness)(i, j));
			}
			//Caso os graus de liberdade sejam fixo e livre (Matriz Kba)
			if (GL_global_1 < 0 && GL_global_2 > 0)
			{
				db.global_stiffness_BA.setValue(-GL_global_1 - 1, GL_global_2 - 1, (*stiffness)(i, j));
			}
		}
	}
}
//Salva variáveis nos pontos de Gauss úteis para descrição lagrangiana atualizada
void Pipe_1::SaveLagrange()
{
	//Loop nos pontos de Gauss
	for (int gauss = 0; gauss < 2; gauss++)
	{
		(*lag_save->Q_i[gauss]) = (*Q_delta[gauss])*(*lag_save->Q_i[gauss]);										//Tensor de rotação
		(*lag_save->u_i[gauss]) = (*lag_save->u_i[gauss]) + *u_delta[gauss];										//Deslocamento
		(*lag_save->kappa_i_ref[gauss]) = *kappa_r[gauss];															//Rotação por unidade de comprimento 
		(*lag_save->alpha_i[gauss]) = 4.0 / (4.0 - dot(*alpha_delta[gauss], *lag_save->alpha_i[gauss]))* (*alpha_delta[gauss] +
			(*lag_save->alpha_i[gauss]) + 0.5*cross(*alpha_delta[gauss], *lag_save->alpha_i[gauss]));				//Vetor rotação de Rodrigues
		(*lag_save->dz_i[gauss]) = *d_z[gauss];
	}
}
//Pré-cálculo de variáveis que é feito uma única vez no início
void Pipe_1::PreCalc()
{
	//Operador constitutivo
	(*D)(0, 0) = db.pipe_sections[section-1]->GA;
	(*D)(1, 1) = db.pipe_sections[section - 1]->GA;
	(*D)(2, 2) = db.pipe_sections[section - 1]->EA;
	(*D)(3, 3) = db.pipe_sections[section - 1]->EI;
	(*D)(4, 4) = db.pipe_sections[section - 1]->EI;
	(*D)(5, 5) = db.pipe_sections[section - 1]->GJ;
	//Identidade de ordem 3
	(*I3)(0, 0) = 1.0;
	(*I3)(1, 1) = 1.0;
	(*I3)(2, 2) = 1.0;
	(*B2)(0, 0) = 1.0;
	(*B2)(1, 1) = 1.0;
	(*B2)(2, 2) = 1.0;
	//Monta as matrizes de rotação
	TransformMatrix();
	//Versor e3r (configuração de referência) - escrito no sistema local do elemento
	(*e3r)(0, 0) = db.nodes[nodes[2] - 1]->ref_coordinates[0] - db.nodes[nodes[0] - 1]->ref_coordinates[0];
	(*e3r)(1, 0) = db.nodes[nodes[2] - 1]->ref_coordinates[1] - db.nodes[nodes[0] - 1]->ref_coordinates[1];
	(*e3r)(2, 0) = db.nodes[nodes[2] - 1]->ref_coordinates[2] - db.nodes[nodes[0] - 1]->ref_coordinates[2];
	*e3r = (1.0 / (norm(*e3r)))*(*e3r);							//Normalização
	*e3rg = *e3r;												//Cópia do e3r (preservado no sistema global)
	*e3r = (*transform3)*(*e3r);								//e3r escrito no sistema local do elemento
	//Massa por unidade de comprimento na configuração de referência
	(*Mr)(0, 0) = db.pipe_sections[section - 1]->Rho;
	(*Mr)(1, 1) = db.pipe_sections[section - 1]->Rho;
	(*Mr)(2, 2) = db.pipe_sections[section - 1]->Rho;
	Mr->MatrixToPtr(pMr, 3);//salvando ptr para uso na dinâmica
	//Inércia por unidade de comprimento na configuração de referência
	(*Jr)(0, 0) = (db.pipe_sections[section - 1]->Rho*((db.pipe_sections[section - 1]->De / 2.0)*(db.pipe_sections[section - 1]->De / 2.0) - (db.pipe_sections[section - 1]->Di / 2.0)*(db.pipe_sections[section - 1]->Di / 2.0)) / 4.0);
	(*Jr)(1, 1) = (db.pipe_sections[section - 1]->Rho*((db.pipe_sections[section - 1]->De / 2.0)*(db.pipe_sections[section - 1]->De / 2.0) - (db.pipe_sections[section - 1]->Di / 2.0)*(db.pipe_sections[section - 1]->Di / 2.0)) / 4.0);
	(*Jr)(2, 2) = (db.pipe_sections[section - 1]->Rho*((db.pipe_sections[section - 1]->De / 2.0)*(db.pipe_sections[section - 1]->De / 2.0) - (db.pipe_sections[section - 1]->Di / 2.0)*(db.pipe_sections[section - 1]->Di / 2.0)) / 2.0);
	Jr->MatrixToPtr(pJr, 3);//salvando ptr para uso na dinâmica
	//Termos de transporte - alterar em situações em que o baricentro não é o polo utilizado no cálculo do momento de inércia
	(*br)(0, 0) = 0.0;
	(*br)(1, 0) = 0.0;
	(*br)(2, 0) = 0.0;
	length = CalculateLength();	//comprimento
	jacobian = length / 2.0;	//jacobiano
	alpha1 = 1.0;				//Peso do método de quadratura Gaussiana
	if (db.environment_exist == true)
		rho_f = db.environment->rho_fluid;
	Aext = PI*De()*De() / 4.0;	//área externa
	Aint = PI*(db.pipe_sections[section - 1]->Di)*(db.pipe_sections[section - 1]->Di) / 4.0;	//área interna
	Cdt = db.pipe_sections[section - 1]->CDt;	//coef. arrasto
	Cdn = db.pipe_sections[section - 1]->CDn;	//coef. arrasto
	C1t = 0.5*De()*rho_f*Cdt;
	C1n = 0.5*De()*rho_f*Cdn;
	if (db.environment_exist == true)
	{
		rho_adt = db.environment->rho_fluid*PI*(db.pipe_sections[section - 1]->De*db.pipe_sections[section - 1]->De)*db.pipe_sections[section - 1]->CAt / 4.0;
		rho_adn = db.environment->rho_fluid*PI*(db.pipe_sections[section - 1]->De*db.pipe_sections[section - 1]->De)*db.pipe_sections[section - 1]->CAn / 4.0;
	}
	else
	{
		rho_adt = 0.0;
		rho_adn = 0.0;
	}
	//Percorre pontos de Gauss para salvar algumas matrizes/vetores de interesse nos cálculos do elemento
	for (int gauss = 0; gauss < 2; gauss++)
	{
		//Ponto localizado nas coordenadas naturais  0.577350269189626 (2 pontos de Gauss)
		if (gauss == 0)
			csi[gauss] = -0.577350269189626;
		if (gauss == 1)
			csi[gauss] = +0.577350269189626;

		N1[gauss] = 0.5*csi[gauss] * (csi[gauss] - 1.0);
		N2[gauss] = 1.0 - csi[gauss] * csi[gauss];
		N3[gauss] = 0.5*csi[gauss] * (1.0 + csi[gauss]);
		dN1[gauss] = (1.0 / jacobian)*(csi[gauss] - 0.5);
		dN2[gauss] = (1.0 / jacobian)*(-2.0*csi[gauss]);
		dN3[gauss] = (1.0 / jacobian)*(0.5 + csi[gauss]);

		(*deltaN[gauss])(0, 0) = dN1[gauss];
		(*deltaN[gauss])(1, 1) = dN1[gauss];
		(*deltaN[gauss])(2, 2) = dN1[gauss];
		(*deltaN[gauss])(3, 3) = dN1[gauss];
		(*deltaN[gauss])(4, 4) = dN1[gauss];
		(*deltaN[gauss])(5, 5) = dN1[gauss];
		(*deltaN[gauss])(6, 3) = N1[gauss];
		(*deltaN[gauss])(7, 4) = N1[gauss];
		(*deltaN[gauss])(8, 5) = N1[gauss];

		(*deltaN[gauss])(0, 6) = dN2[gauss];
		(*deltaN[gauss])(1, 7) = dN2[gauss];
		(*deltaN[gauss])(2, 8) = dN2[gauss];
		(*deltaN[gauss])(3, 9) = dN2[gauss];
		(*deltaN[gauss])(4, 10) = dN2[gauss];
		(*deltaN[gauss])(5, 11) = dN2[gauss];
		(*deltaN[gauss])(6, 9) = N2[gauss];
		(*deltaN[gauss])(7, 10) = N2[gauss];
		(*deltaN[gauss])(8, 11) = N2[gauss];

		(*deltaN[gauss])(0, 12) = dN3[gauss];
		(*deltaN[gauss])(1, 13) = dN3[gauss];
		(*deltaN[gauss])(2, 14) = dN3[gauss];
		(*deltaN[gauss])(3, 15) = dN3[gauss];
		(*deltaN[gauss])(4, 16) = dN3[gauss];
		(*deltaN[gauss])(5, 17) = dN3[gauss];
		(*deltaN[gauss])(6, 15) = N3[gauss];
		(*deltaN[gauss])(7, 16) = N3[gauss];
		(*deltaN[gauss])(8, 17) = N3[gauss];

		(*N[gauss])(0, 0) = N1[gauss];
		(*N[gauss])(1, 1) = N1[gauss];
		(*N[gauss])(2, 2) = N1[gauss];
		(*N[gauss])(3, 3) = N1[gauss];
		(*N[gauss])(4, 4) = N1[gauss];
		(*N[gauss])(5, 5) = N1[gauss];

		(*N[gauss])(0, 6) = N2[gauss];
		(*N[gauss])(1, 7) = N2[gauss];
		(*N[gauss])(2, 8) = N2[gauss];
		(*N[gauss])(3, 9) = N2[gauss];
		(*N[gauss])(4, 10) = N2[gauss];
		(*N[gauss])(5, 11) = N2[gauss];

		(*N[gauss])(0, 12) = N3[gauss];
		(*N[gauss])(1, 13) = N3[gauss];
		(*N[gauss])(2, 14) = N3[gauss];
		(*N[gauss])(3, 15) = N3[gauss];
		(*N[gauss])(4, 16) = N3[gauss];
		(*N[gauss])(5, 17) = N3[gauss];

		(*UpsilonN[gauss])(0, 0) = N1[gauss];
		(*UpsilonN[gauss])(1, 1) = N1[gauss];
		(*UpsilonN[gauss])(2, 2) = N1[gauss];
		(*UpsilonN[gauss])(3, 3) = N1[gauss];
		(*UpsilonN[gauss])(4, 4) = N1[gauss];
		(*UpsilonN[gauss])(5, 5) = N1[gauss];

		(*UpsilonN[gauss])(0, 6) = N2[gauss];
		(*UpsilonN[gauss])(1, 7) = N2[gauss];
		(*UpsilonN[gauss])(2, 8) = N2[gauss];
		(*UpsilonN[gauss])(3, 9) = N2[gauss];
		(*UpsilonN[gauss])(4, 10) = N2[gauss];
		(*UpsilonN[gauss])(5, 11) = N2[gauss];

		(*UpsilonN[gauss])(0, 12) = N3[gauss];
		(*UpsilonN[gauss])(1, 13) = N3[gauss];
		(*UpsilonN[gauss])(2, 14) = N3[gauss];
		(*UpsilonN[gauss])(3, 15) = N3[gauss];
		(*UpsilonN[gauss])(4, 16) = N3[gauss];
		(*UpsilonN[gauss])(5, 17) = N3[gauss];

		(*UpsilonN[gauss])(6, 0) = dN1[gauss];
		(*UpsilonN[gauss])(7, 1) = dN1[gauss];
		(*UpsilonN[gauss])(8, 2) = dN1[gauss];
		(*UpsilonN[gauss])(9, 3) = dN1[gauss];
		(*UpsilonN[gauss])(10, 4) = dN1[gauss];
		(*UpsilonN[gauss])(11, 5) = dN1[gauss];

		(*UpsilonN[gauss])(6, 6) = dN2[gauss];
		(*UpsilonN[gauss])(7, 7) = dN2[gauss];
		(*UpsilonN[gauss])(8, 8) = dN2[gauss];
		(*UpsilonN[gauss])(9, 9) = dN2[gauss];
		(*UpsilonN[gauss])(10, 10) = dN2[gauss];
		(*UpsilonN[gauss])(11, 11) = dN2[gauss];

		(*UpsilonN[gauss])(6, 12) = dN3[gauss];
		(*UpsilonN[gauss])(7, 13) = dN3[gauss];
		(*UpsilonN[gauss])(8, 14) = dN3[gauss];
		(*UpsilonN[gauss])(9, 15) = dN3[gauss];
		(*UpsilonN[gauss])(10, 16) = dN3[gauss];
		(*UpsilonN[gauss])(11, 17) = dN3[gauss];
	}//end of gauss
}

//Monta carregamentos de campo
void Pipe_1::MountFieldLoading()
{
	if (db.environment_exist == true)
	{
		if (db.environment->g_exist == true)
		{
			load_multiplier = 1.0;
			l_factor = db.environment->bool_g.GetLinearFactorAtCurrentTime();

			//Loop nos pontos de Gauss
			for (int gauss = 0; gauss < 2; gauss++)
			{
				if (db.environment->ocean_data_exist == true)
				{
					//Posição do ponto de Gauss
					(*zi[gauss])(0, 0) = (db.nodes[nodes[0] - 1]->copy_coordinates[0])* N1[gauss] +
						(db.nodes[nodes[1] - 1]->copy_coordinates[0])* N2[gauss] +
						(db.nodes[nodes[2] - 1]->copy_coordinates[0])* N3[gauss] +
						(db.nodes[nodes[0] - 1]->displacements[0])* N1[gauss] +
						(db.nodes[nodes[1] - 1]->displacements[0])* N2[gauss] +
						(db.nodes[nodes[2] - 1]->displacements[0])* N3[gauss];
					(*zi[gauss])(1, 0) = (db.nodes[nodes[0] - 1]->copy_coordinates[1])* N1[gauss] +
						(db.nodes[nodes[1] - 1]->copy_coordinates[1])* N2[gauss] +
						(db.nodes[nodes[2] - 1]->copy_coordinates[1])* N3[gauss] +
						(db.nodes[nodes[0] - 1]->displacements[1])* N1[gauss] +
						(db.nodes[nodes[1] - 1]->displacements[1])* N2[gauss] +
						(db.nodes[nodes[2] - 1]->displacements[1])* N3[gauss];
					(*zi[gauss])(2, 0) = (db.nodes[nodes[0] - 1]->copy_coordinates[2])* N1[gauss] +
						(db.nodes[nodes[1] - 1]->copy_coordinates[2])* N2[gauss] +
						(db.nodes[nodes[2] - 1]->copy_coordinates[2])* N3[gauss] +
						(db.nodes[nodes[0] - 1]->displacements[2])* N1[gauss] +
						(db.nodes[nodes[1] - 1]->displacements[2])* N2[gauss] +
						(db.nodes[nodes[2] - 1]->displacements[2])* N3[gauss];

					//Verifica profundidade e assinala o valor correto para o rho
					depth = dot((*zi[gauss] - db.environment->surface_position), db.environment->G)*(1.0 / norm(db.environment->G));
					if (depth > 0)
						mult = l_factor*load_multiplier*alpha1*jacobian*(db.pipe_sections[section - 1]->Rho - db.environment->rho_fluid*Aext);
					else
						mult = l_factor*load_multiplier*alpha1*jacobian*(db.pipe_sections[section - 1]->Rho);
				}
				else
					mult = l_factor*load_multiplier*alpha1*jacobian*(db.pipe_sections[section - 1]->Rho);

				if (mult != 0.0)
				{
					(*e_loading)(0, 0) += mult * N1[gauss] * db.environment->G(0, 0);
					(*e_loading)(1, 0) += mult * N1[gauss] * db.environment->G(1, 0);
					(*e_loading)(2, 0) += mult * N1[gauss] * db.environment->G(2, 0);

					(*e_loading)(6, 0) += mult * N2[gauss] * db.environment->G(0, 0);
					(*e_loading)(7, 0) += mult * N2[gauss] * db.environment->G(1, 0);
					(*e_loading)(8, 0) += mult * N2[gauss] * db.environment->G(2, 0);

					(*e_loading)(12, 0) += mult * N3[gauss] * db.environment->G(0, 0);
					(*e_loading)(13, 0) += mult * N3[gauss] * db.environment->G(1, 0);
					(*e_loading)(14, 0) += mult * N3[gauss] * db.environment->G(2, 0);
				}
			}//end of gauss loop
		}
	}
}

//Monta a parte estática dos esforços de Morison - correnteza marítima
void Pipe_1::MountSeaCurrentLoading()
{
	if (db.environment_exist == true)
	{
		if (db.environment->ocean_data_exist == true)
		{
			load_multiplier = 1.0;
			l_factor = db.environment->bool_current.GetLinearFactorAtCurrentTime();

			mult = l_factor*load_multiplier*alpha1*jacobian;
			if (mult != 0.0)
			{
				Matrix du_i(3);
				Matrix ddu_i(3);
				Matrix force(6);
				double temp_a4, temp_a5, temp_a6;
				if (typeid(*db.solution[db.current_solution_number - 1]) == typeid(Dynamic))
				{
					Dynamic* ptr_sol = static_cast<Dynamic*>(db.solution[db.current_solution_number - 1]);

					temp_a4 = ptr_sol->a4;
					temp_a5 = ptr_sol->a5;
					temp_a6 = ptr_sol->a6;
				}
				else
				{
					temp_a4 = 0;
					temp_a5 = 0;
					temp_a6 = 0;
				}
				//Loop nos pontos de Gauss
				for (int gauss = 0; gauss < 2; gauss++)
				{
					//Posição do ponto de Gauss
					(*zi[gauss])(0, 0) = (db.nodes[nodes[0] - 1]->copy_coordinates[0])* N1[gauss] +
						(db.nodes[nodes[1] - 1]->copy_coordinates[0])* N2[gauss] +
						(db.nodes[nodes[2] - 1]->copy_coordinates[0])* N3[gauss] +
						(db.nodes[nodes[0] - 1]->displacements[0])* N1[gauss] +
						(db.nodes[nodes[1] - 1]->displacements[0])* N2[gauss] +
						(db.nodes[nodes[2] - 1]->displacements[0])* N3[gauss];
					(*zi[gauss])(1, 0) = (db.nodes[nodes[0] - 1]->copy_coordinates[1])* N1[gauss] +
						(db.nodes[nodes[1] - 1]->copy_coordinates[1])* N2[gauss] +
						(db.nodes[nodes[2] - 1]->copy_coordinates[1])* N3[gauss] +
						(db.nodes[nodes[0] - 1]->displacements[1])* N1[gauss] +
						(db.nodes[nodes[1] - 1]->displacements[1])* N2[gauss] +
						(db.nodes[nodes[2] - 1]->displacements[1])* N3[gauss];
					(*zi[gauss])(2, 0) = (db.nodes[nodes[0] - 1]->copy_coordinates[2])* N1[gauss] +
						(db.nodes[nodes[1] - 1]->copy_coordinates[2])* N2[gauss] +
						(db.nodes[nodes[2] - 1]->copy_coordinates[2])* N3[gauss] +
						(db.nodes[nodes[0] - 1]->displacements[2])* N1[gauss] +
						(db.nodes[nodes[1] - 1]->displacements[2])* N2[gauss] +
						(db.nodes[nodes[2] - 1]->displacements[2])* N3[gauss];
					//Velocidade da correnteza no ponto de Gauss (cota z)
					//Verifica profundidade
					depth = dot((*zi[gauss] - db.environment->surface_position), db.environment->G)*(1.0 / norm(db.environment->G));
					if (depth > 0)
					{
						*vel[gauss] = db.environment->VelocityAt(depth);
						(du_i)(0, 0) = (db.nodes[nodes[0] - 1]->copy_vel[0])* N1[gauss] +
							(db.nodes[nodes[1] - 1]->copy_vel[0])* N2[gauss] +
							(db.nodes[nodes[2] - 1]->copy_vel[0])* N3[gauss];
						(du_i)(1, 0) = (db.nodes[nodes[0] - 1]->copy_vel[1])* N1[gauss] +
							(db.nodes[nodes[1] - 1]->copy_vel[1])* N2[gauss] +
							(db.nodes[nodes[2] - 1]->copy_vel[1])* N3[gauss];
						(du_i)(2, 0) = (db.nodes[nodes[0] - 1]->copy_vel[2])* N1[gauss] +
							(db.nodes[nodes[1] - 1]->copy_vel[2])* N2[gauss] +
							(db.nodes[nodes[2] - 1]->copy_vel[2])* N3[gauss];

						(ddu_i)(0, 0) = (db.nodes[nodes[0] - 1]->copy_accel[0])* N1[gauss] +
							(db.nodes[nodes[1] - 1]->copy_accel[0])* N2[gauss] +
							(db.nodes[nodes[2] - 1]->copy_accel[0])* N3[gauss];
						(ddu_i)(1, 0) = (db.nodes[nodes[0] - 1]->copy_accel[1])* N1[gauss] +
							(db.nodes[nodes[1] - 1]->copy_accel[1])* N2[gauss] +
							(db.nodes[nodes[2] - 1]->copy_accel[1])* N3[gauss];
						(ddu_i)(2, 0) = (db.nodes[nodes[0] - 1]->copy_accel[2])* N1[gauss] +
							(db.nodes[nodes[1] - 1]->copy_accel[2])* N2[gauss] +
							(db.nodes[nodes[2] - 1]->copy_accel[2])* N3[gauss];

						(*vel[gauss]) = (*transform3)*(*vel[gauss]);
						du_i = (*transform3)*du_i;
						ddu_i = (*transform3)*ddu_i;

						EvaluateMorisonContributions(temp_v, &temp_a4, &temp_a5, &temp_a6, &C1t, &C1n, vel[gauss]->getMatrix(), e3r->getMatrix(), lag_save->alpha_i[gauss]->getMatrix(),
							alpha_delta[gauss]->getMatrix(), u_delta[gauss]->getMatrix(), du_i.getMatrix(), ddu_i.getMatrix(), force.getMatrix(), pL);
						L[gauss]->PtrToMatrix(pL, 6);

						(*loading_stiffness) = (*loading_stiffness) + mult*(transp(*N[gauss]))*(*L[gauss])*(*N[gauss]);
						(*morison_loading) = (*morison_loading) + mult*transp(*N[gauss])*force;
					}
				
				}//end of gauss loop
				(*stiffness) = (*stiffness) - transp(*transform)*(*loading_stiffness)*(*transform);	//Acréscimo na matriz de rigidez
				(*e_loading) = (*e_loading) + transp(*transform)*(*morison_loading);				//Forças de Morison
			}
		}
	}
}

//Monta carregamentos de pressão interna/externa no tubo
void Pipe_1::MountPipeSpecialLoads(int l_number)
{
	//Função para montar carregamentos de pressão em tubos. Atende simulações estática e dinâmica
	//////////////////Efeito da Pressão////////////////////////////
	p0i = db.loads[l_number - 1]->GetValueAt(db.last_converged_time + db.current_time_step, 0);
	p0e = db.loads[l_number - 1]->GetValueAt(db.last_converged_time + db.current_time_step, 1);
	rhoi = db.loads[l_number - 1]->GetValueAt(db.last_converged_time + db.current_time_step, 2);
	rhoe = db.loads[l_number - 1]->GetValueAt(db.last_converged_time + db.current_time_step, 3);

	//Multiplicador de carregamentos - caso estático
	mult = alpha1*jacobian;

	//Loop nos pontos de Gauss
	for (int gauss = 0; gauss < 2; gauss++)
	{
		*kip[gauss] = ((*Q_delta[gauss])* (*lag_save->Q_i[gauss]))*(*kappa_r[gauss]);	//Rotação específica no ponto de gauss no instante i+1 (Lag. Atualizada)
		*e3ip[gauss] = ((*Q_delta[gauss])* (*lag_save->Q_i[gauss]))*(*e3r);				//Orientação do elemento no ponto de gauss no instante i+1 (Lag. Atualizada)
		//Forças - sistema local (elemento)
		*temp_f[gauss] = -p0i*Aint*cross(*kip[gauss], *e3ip[gauss]);
		//Momentos - sistema local (elemento)
		*temp_m[gauss] = -p0i*Aint*transp(*Xi_delta[gauss])*cross(*d_z[gauss], *e3ip[gauss]);
		//Vetor de esforços externos generalizados - sistema local (elemento)
		(*temp_l[gauss])(0, 0) = (*temp_f[gauss])(0, 0);
		(*temp_l[gauss])(1, 0) = (*temp_f[gauss])(1, 0);
		(*temp_l[gauss])(2, 0) = (*temp_f[gauss])(2, 0);
		(*temp_l[gauss])(3, 0) = (*temp_m[gauss])(0, 0);
		(*temp_l[gauss])(4, 0) = (*temp_m[gauss])(1, 0);
		(*temp_l[gauss])(5, 0) = (*temp_m[gauss])(2, 0);
		*E3ip[gauss] = skew(*e3ip[gauss]);
		*Kip[gauss] = skew(*kip[gauss]);
		*O1 = -0.5*g*(transp(*Xi_delta[gauss])*(cross(*e3ip[gauss], *d_z[gauss]))*transp(*alpha_delta[gauss]) - skew(cross(*e3ip[gauss], *d_z[gauss])));
		*K1ua = (*Kip[gauss]) * (*E3ip[gauss])* (*Xi_delta[gauss]) + (*E3ip[gauss])*(*d_Xi_delta[gauss]) - 1.0*(*E3ip[gauss])*skew(*kip[gauss])*(*Xi_delta[gauss]);
		*K1aa = *O1 + transp((*Xi_delta[gauss]))*(*d_Z[gauss])*(*E3ip[gauss])*(*Xi_delta[gauss]);
		*K2ua = (*E3ip[gauss])*(*Xi_delta[gauss]);
		*K2au = transp(*Xi_delta[gauss])*(*E3ip[gauss]);
		for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			(*Kext)(i, j + 3) = (*K1ua)(i, j);
		for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			(*Kext)(i + 3, j + 3) = (*K1aa)(i, j);
		for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			(*Kext)(i, j + 9) = (*K2ua)(i, j);
		for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			(*Kext)(i + 3, j + 6) = (*K2au)(i, j);
		//Carregamentos no elemento
		(*P_loading) = (*P_loading) - mult* transp(*transform)*(transp(*N[gauss]))*(*temp_l[gauss]);	//Esforço desbalanceado
		(*e_loading) = (*e_loading) + mult* transp(*transform)*(transp(*N[gauss]))*(*temp_l[gauss]);	//Esforço externo
		//Contribuição para a matriz de rigidez
		(*stiffness) = (*stiffness) - mult*p0i*Aint* transp(*transform) * ((transp(*N[gauss]))*(*Kext)*(*UpsilonN[gauss])) * (*transform);
	}
}
//Comprimento do elemento
double Pipe_1::CalculateLength()
{
	double x1 = db.nodes[nodes[0] - 1]->ref_coordinates[0];
	double y1 = db.nodes[nodes[0] - 1]->ref_coordinates[1];
	double z1 = db.nodes[nodes[0] - 1]->ref_coordinates[2];
	double x2 = db.nodes[nodes[2] - 1]->ref_coordinates[0];
	double y2 = db.nodes[nodes[2] - 1]->ref_coordinates[1];
	double z2 = db.nodes[nodes[2] - 1]->ref_coordinates[2];
	return sqrt((x2 - x1)*(x2 - x1) + (y2 - y1)*(y2 - y1) + (z2 - z1)*(z2 - z1));
}
void Pipe_1::Zeros()
{
	zeros(stiffness);
	zeros(loading_stiffness);
	zeros(mass);
	zeros(damping);
	zeros(damping_loading);
	zeros(i_loading);
	zeros(inertial_loading);
	zeros(morison_loading);
	zeros(e_loading);
	zeros(P_loading);
	kinetic_energy = 0.0;
	strain_energy = 0.0;
	potential_gravitational_energy = 0.0;
}

//Retorna o diâmetro externo
double Pipe_1::De()
{
	return db.pipe_sections[section - 1]->De;
}

//Monta a matriz de massa
void Pipe_1::MountMassModal()
{
	zeros(mass_modal);
	zeros(damping_modal);
	//Loop nos pontos de Gauss
	for (int gauss = 0; gauss < 2; gauss++)
	{
		if (db.environment_exist == true && db.environment->ocean_data_exist == true && db.environment->g_exist)
		{
			//Avalia se o tubo está imerso ou não para modificar a massa adicional
			//Posição do ponto de Gauss
			(*zi[gauss])(0, 0) = (db.nodes[nodes[0] - 1]->copy_coordinates[0])* N1[gauss] +
				(db.nodes[nodes[1] - 1]->copy_coordinates[0])* N2[gauss] +
				(db.nodes[nodes[2] - 1]->copy_coordinates[0])* N3[gauss] +
				(db.nodes[nodes[0] - 1]->displacements[0])* N1[gauss] +
				(db.nodes[nodes[1] - 1]->displacements[0])* N2[gauss] +
				(db.nodes[nodes[2] - 1]->displacements[0])* N3[gauss];
			(*zi[gauss])(1, 0) = (db.nodes[nodes[0] - 1]->copy_coordinates[1])* N1[gauss] +
				(db.nodes[nodes[1] - 1]->copy_coordinates[1])* N2[gauss] +
				(db.nodes[nodes[2] - 1]->copy_coordinates[1])* N3[gauss] +
				(db.nodes[nodes[0] - 1]->displacements[1])* N1[gauss] +
				(db.nodes[nodes[1] - 1]->displacements[1])* N2[gauss] +
				(db.nodes[nodes[2] - 1]->displacements[1])* N3[gauss];
			(*zi[gauss])(2, 0) = (db.nodes[nodes[0] - 1]->copy_coordinates[2])* N1[gauss] +
				(db.nodes[nodes[1] - 1]->copy_coordinates[2])* N2[gauss] +
				(db.nodes[nodes[2] - 1]->copy_coordinates[2])* N3[gauss] +
				(db.nodes[nodes[0] - 1]->displacements[2])* N1[gauss] +
				(db.nodes[nodes[1] - 1]->displacements[2])* N2[gauss] +
				(db.nodes[nodes[2] - 1]->displacements[2])* N3[gauss];
			//Velocidade da correnteza no ponto de Gauss (cota z)
			//Verifica profundidade
			depth = dot((*zi[gauss] - db.environment->surface_position), db.environment->G)*(1.0 / norm(db.environment->G));
			if (depth > 0)
			{
				//Massa por unidade de comprimento na configuração de referência
				(*Mr)(0, 0) = db.pipe_sections[section - 1]->Rho + rho_adn;
				(*Mr)(1, 1) = db.pipe_sections[section - 1]->Rho + rho_adn;
				(*Mr)(2, 2) = db.pipe_sections[section - 1]->Rho + rho_adt;
			}
			else
			{
				//Massa por unidade de comprimento na configuração de referência
				(*Mr)(0, 0) = db.pipe_sections[section - 1]->Rho;
				(*Mr)(1, 1) = db.pipe_sections[section - 1]->Rho;
				(*Mr)(2, 2) = db.pipe_sections[section - 1]->Rho;
			}
			
		}
		else
		{
			//Massa por unidade de comprimento na configuração de referência
			(*Mr)(0, 0) = db.pipe_sections[section - 1]->Rho;
			(*Mr)(1, 1) = db.pipe_sections[section - 1]->Rho;
			(*Mr)(2, 2) = db.pipe_sections[section - 1]->Rho;
		}
		Mr->MatrixToPtr(pMr, 3);//salvando ptr para uso na dinâmica
		//Cálculo do integrando no ponto de Gauss
		EvaluateMassModal(temp_v, lag_save->alpha_i[gauss]->getMatrix(), pJr, pMr, br->getMatrix(), pDdT);
		//Transformando o operador tangente em Matrix
		DdT->PtrToMatrix(pDdT, 6);
		(*mass_modal) = (*mass_modal) + (alpha1*jacobian)*transp(*N[gauss])*(*DdT)*(*N[gauss]);
	}//end of gauss
	//Conversão para o sistema global
	(*mass_modal) = (transp(*transform)*(*mass_modal))*(*transform);
}

//Monta a matriz de amortecimento para realização da análise modal
void Pipe_1::MountDampingModal()
{
	zeros(mass_modal);
	zeros(damping_modal);
	//TODO
}

void Pipe_1::EvaluateMassModal(double* v, double* alphai, double** Jr, double** Mr, double* br, double** matrixm)
{
	v[140] = Power(alphai[2], 2);
	v[139] = 0.5e0*alphai[0] * alphai[2];
	v[138] = 0.5e0*alphai[1];
	v[137] = Power(alphai[1], 2);
	v[141] = v[137] + v[140];
	v[136] = alphai[0] * v[138];
	v[135] = Power(alphai[0], 2);
	v[83] = alphai[2] * v[138];
	v[70] = 4e0 / (4e0 + v[135] + v[141]);
	v[142] = -0.5e0*v[70];
	v[73] = 1e0 + v[141] * v[142];
	v[74] = (-alphai[2] + v[136])*v[70];
	v[75] = (alphai[1] + v[139])*v[70];
	v[107] = Jr[0][2] * v[73] + Jr[1][2] * v[74] + Jr[2][2] * v[75];
	v[106] = Jr[0][1] * v[73] + Jr[1][1] * v[74] + Jr[2][1] * v[75];
	v[105] = Jr[0][0] * v[73] + Jr[1][0] * v[74] + Jr[2][0] * v[75];
	v[89] = Mr[0][2] * v[73] + Mr[1][2] * v[74] + Mr[2][2] * v[75];
	v[88] = Mr[0][1] * v[73] + Mr[1][1] * v[74] + Mr[2][1] * v[75];
	v[87] = Mr[0][0] * v[73] + Mr[1][0] * v[74] + Mr[2][0] * v[75];
	v[77] = (alphai[2] + v[136])*v[70];
	v[79] = 1e0 + (v[135] + v[140])*v[142];
	v[80] = v[70] * (-alphai[0] + v[83]);
	v[113] = Jr[0][2] * v[77] + Jr[1][2] * v[79] + Jr[2][2] * v[80];
	v[112] = Jr[0][1] * v[77] + Jr[1][1] * v[79] + Jr[2][1] * v[80];
	v[111] = Jr[0][0] * v[77] + Jr[1][0] * v[79] + Jr[2][0] * v[80];
	v[95] = Mr[0][2] * v[77] + Mr[1][2] * v[79] + Mr[2][2] * v[80];
	v[94] = Mr[0][1] * v[77] + Mr[1][1] * v[79] + Mr[2][1] * v[80];
	v[93] = Mr[0][0] * v[77] + Mr[1][0] * v[79] + Mr[2][0] * v[80];
	v[82] = (-alphai[1] + v[139])*v[70];
	v[84] = v[70] * (alphai[0] + v[83]);
	v[85] = 1e0 + (v[135] + v[137])*v[142];
	v[119] = Jr[0][2] * v[82] + Jr[1][2] * v[84] + Jr[2][2] * v[85];
	v[118] = Jr[0][1] * v[82] + Jr[1][1] * v[84] + Jr[2][1] * v[85];
	v[117] = Jr[0][0] * v[82] + Jr[1][0] * v[84] + Jr[2][0] * v[85];
	v[101] = Mr[0][2] * v[82] + Mr[1][2] * v[84] + Mr[2][2] * v[85];
	v[100] = Mr[0][1] * v[82] + Mr[1][1] * v[84] + Mr[2][1] * v[85];
	v[99] = Mr[0][0] * v[82] + Mr[1][0] * v[84] + Mr[2][0] * v[85];
	v[86] = v[73] * v[87] + v[74] * v[88] + v[75] * v[89];
	v[90] = v[77] * v[87] + v[79] * v[88] + v[80] * v[89];
	v[91] = v[82] * v[87] + v[84] * v[88] + v[85] * v[89];
	v[92] = v[73] * v[93] + v[74] * v[94] + v[75] * v[95];
	v[96] = v[77] * v[93] + v[79] * v[94] + v[80] * v[95];
	v[97] = v[82] * v[93] + v[84] * v[94] + v[85] * v[95];
	v[98] = v[100] * v[74] + v[101] * v[75] + v[73] * v[99];
	v[102] = v[100] * v[79] + v[101] * v[80] + v[77] * v[99];
	v[103] = v[100] * v[84] + v[101] * v[85] + v[82] * v[99];
	v[122] = br[0] * v[73] + br[1] * v[74] + br[2] * v[75];
	v[123] = br[0] * v[77] + br[1] * v[79] + br[2] * v[80];
	v[124] = br[0] * v[82] + br[1] * v[84] + br[2] * v[85];
	v[125] = v[124] * v[90] - v[123] * v[91];
	v[126] = -(v[124] * v[86]) + v[122] * v[91];
	v[127] = v[123] * v[86] - v[122] * v[90];
	v[128] = v[124] * v[96] - v[123] * v[97];
	v[129] = -(v[124] * v[92]) + v[122] * v[97];
	v[130] = v[123] * v[92] - v[122] * v[96];
	v[131] = -(v[103] * v[123]) + v[102] * v[124];
	v[132] = v[103] * v[122] - v[124] * v[98];
	v[133] = -(v[102] * v[122]) + v[123] * v[98];
	matrixm[0][0] = v[86];
	matrixm[0][1] = v[90];
	matrixm[0][2] = v[91];
	matrixm[0][3] = v[125];
	matrixm[0][4] = v[126];
	matrixm[0][5] = v[127];
	matrixm[1][0] = v[92];
	matrixm[1][1] = v[96];
	matrixm[1][2] = v[97];
	matrixm[1][3] = v[128];
	matrixm[1][4] = v[129];
	matrixm[1][5] = v[130];
	matrixm[2][0] = v[98];
	matrixm[2][1] = v[102];
	matrixm[2][2] = v[103];
	matrixm[2][3] = v[131];
	matrixm[2][4] = v[132];
	matrixm[2][5] = v[133];
	matrixm[3][0] = -v[125];
	matrixm[3][1] = -v[126];
	matrixm[3][2] = -v[127];
	matrixm[3][3] = v[105] * v[73] + v[106] * v[74] + v[107] * v[75];
	matrixm[3][4] = v[105] * v[77] + v[106] * v[79] + v[107] * v[80];
	matrixm[3][5] = v[105] * v[82] + v[106] * v[84] + v[107] * v[85];
	matrixm[4][0] = -v[128];
	matrixm[4][1] = -v[129];
	matrixm[4][2] = -v[130];
	matrixm[4][3] = v[111] * v[73] + v[112] * v[74] + v[113] * v[75];
	matrixm[4][4] = v[111] * v[77] + v[112] * v[79] + v[113] * v[80];
	matrixm[4][5] = v[111] * v[82] + v[112] * v[84] + v[113] * v[85];
	matrixm[5][0] = -v[131];
	matrixm[5][1] = -v[132];
	matrixm[5][2] = -v[133];
	matrixm[5][3] = v[117] * v[73] + v[118] * v[74] + v[119] * v[75];
	matrixm[5][4] = v[117] * v[77] + v[118] * v[79] + v[119] * v[80];
	matrixm[5][5] = v[117] * v[82] + v[118] * v[84] + v[119] * v[85];
};

//Monta a matriz de massa
void Pipe_1::MountMass()
{
	Matrix omega_i(3);
	Matrix domega_i(3);
	Matrix du_i(3);
	Matrix ddu_i(3);
	//Loop nos pontos de Gauss
	for (int gauss = 0; gauss < 2; gauss++)
	{
		//Calculando as variáveis necessárias para obter o valor do integrando							
		(omega_i)(0, 0) =	(db.nodes[nodes[0] - 1]->copy_vel[3])* N1[gauss] +
							(db.nodes[nodes[1] - 1]->copy_vel[3])* N2[gauss] +
							(db.nodes[nodes[2] - 1]->copy_vel[3])* N3[gauss];
		(omega_i)(1, 0) =	(db.nodes[nodes[0] - 1]->copy_vel[4])* N1[gauss] +
							(db.nodes[nodes[1] - 1]->copy_vel[4])* N2[gauss] +
							(db.nodes[nodes[2] - 1]->copy_vel[4])* N3[gauss];
		(omega_i)(2, 0) =	(db.nodes[nodes[0] - 1]->copy_vel[5])* N1[gauss] +
							(db.nodes[nodes[1] - 1]->copy_vel[5])* N2[gauss] +
							(db.nodes[nodes[2] - 1]->copy_vel[5])* N3[gauss];

		(domega_i)(0, 0) =	(db.nodes[nodes[0] - 1]->copy_accel[3])* N1[gauss] +
							(db.nodes[nodes[1] - 1]->copy_accel[3])* N2[gauss] +
							(db.nodes[nodes[2] - 1]->copy_accel[3])* N3[gauss];
		(domega_i)(1, 0) =	(db.nodes[nodes[0] - 1]->copy_accel[4])* N1[gauss] +
							(db.nodes[nodes[1] - 1]->copy_accel[4])* N2[gauss] +
							(db.nodes[nodes[2] - 1]->copy_accel[4])* N3[gauss];
		(domega_i)(2, 0) =	(db.nodes[nodes[0] - 1]->copy_accel[5])* N1[gauss] +
							(db.nodes[nodes[1] - 1]->copy_accel[5])* N2[gauss] +
							(db.nodes[nodes[2] - 1]->copy_accel[5])* N3[gauss];

		(du_i)(0, 0) =	(db.nodes[nodes[0] - 1]->copy_vel[0])* N1[gauss] +
						(db.nodes[nodes[1] - 1]->copy_vel[0])* N2[gauss] +
						(db.nodes[nodes[2] - 1]->copy_vel[0])* N3[gauss];
		(du_i)(1, 0) =	(db.nodes[nodes[0] - 1]->copy_vel[1])* N1[gauss] +
						(db.nodes[nodes[1] - 1]->copy_vel[1])* N2[gauss] +
						(db.nodes[nodes[2] - 1]->copy_vel[1])* N3[gauss];
		(du_i)(2, 0) =	(db.nodes[nodes[0] - 1]->copy_vel[2])* N1[gauss] +
						(db.nodes[nodes[1] - 1]->copy_vel[2])* N2[gauss] +
						(db.nodes[nodes[2] - 1]->copy_vel[2])* N3[gauss];

		(ddu_i)(0, 0) = (db.nodes[nodes[0] - 1]->copy_accel[0])* N1[gauss] +
						(db.nodes[nodes[1] - 1]->copy_accel[0])* N2[gauss] +
						(db.nodes[nodes[2] - 1]->copy_accel[0])* N3[gauss];
		(ddu_i)(1, 0) = (db.nodes[nodes[0] - 1]->copy_accel[1])* N1[gauss] +
						(db.nodes[nodes[1] - 1]->copy_accel[1])* N2[gauss] +
						(db.nodes[nodes[2] - 1]->copy_accel[1])* N3[gauss];
		(ddu_i)(2, 0) = (db.nodes[nodes[0] - 1]->copy_accel[2])* N1[gauss] +
						(db.nodes[nodes[1] - 1]->copy_accel[2])* N2[gauss] +
						(db.nodes[nodes[2] - 1]->copy_accel[2])* N3[gauss];
		
		//Para o sistema do elemento
		omega_i = (*transform3)*omega_i;
		domega_i = (*transform3)*domega_i;
		du_i = (*transform3)*du_i;
		ddu_i = (*transform3)*ddu_i;
		if (db.environment_exist == true && db.environment->ocean_data_exist == true && db.environment->g_exist)
		{
			//Avalia se o tubo está imerso ou não para modificar a massa adicional
			//Posição do ponto de Gauss
			(*zi[gauss])(0, 0) = (db.nodes[nodes[0] - 1]->copy_coordinates[0])* N1[gauss] +
				(db.nodes[nodes[1] - 1]->copy_coordinates[0])* N2[gauss] +
				(db.nodes[nodes[2] - 1]->copy_coordinates[0])* N3[gauss] +
				(db.nodes[nodes[0] - 1]->displacements[0])* N1[gauss] +
				(db.nodes[nodes[1] - 1]->displacements[0])* N2[gauss] +
				(db.nodes[nodes[2] - 1]->displacements[0])* N3[gauss];
			(*zi[gauss])(1, 0) = (db.nodes[nodes[0] - 1]->copy_coordinates[1])* N1[gauss] +
				(db.nodes[nodes[1] - 1]->copy_coordinates[1])* N2[gauss] +
				(db.nodes[nodes[2] - 1]->copy_coordinates[1])* N3[gauss] +
				(db.nodes[nodes[0] - 1]->displacements[1])* N1[gauss] +
				(db.nodes[nodes[1] - 1]->displacements[1])* N2[gauss] +
				(db.nodes[nodes[2] - 1]->displacements[1])* N3[gauss];
			(*zi[gauss])(2, 0) = (db.nodes[nodes[0] - 1]->copy_coordinates[2])* N1[gauss] +
				(db.nodes[nodes[1] - 1]->copy_coordinates[2])* N2[gauss] +
				(db.nodes[nodes[2] - 1]->copy_coordinates[2])* N3[gauss] +
				(db.nodes[nodes[0] - 1]->displacements[2])* N1[gauss] +
				(db.nodes[nodes[1] - 1]->displacements[2])* N2[gauss] +
				(db.nodes[nodes[2] - 1]->displacements[2])* N3[gauss];
			//Velocidade da correnteza no ponto de Gauss (cota z)
			//Verifica profundidade
			depth = dot((*zi[gauss] - db.environment->surface_position), db.environment->G)*(1.0 / norm(db.environment->G));
			if (depth > 0)
			{
				//Massa por unidade de comprimento na configuração de referência
				(*Mr)(0, 0) = db.pipe_sections[section - 1]->Rho + rho_adn;
				(*Mr)(1, 1) = db.pipe_sections[section - 1]->Rho + rho_adn;
				(*Mr)(2, 2) = db.pipe_sections[section - 1]->Rho + rho_adt;
			}
			else
			{
				//Massa por unidade de comprimento na configuração de referência
				(*Mr)(0, 0) = db.pipe_sections[section - 1]->Rho;
				(*Mr)(1, 1) = db.pipe_sections[section - 1]->Rho;
				(*Mr)(2, 2) = db.pipe_sections[section - 1]->Rho;
			}
		}
		else
		{
			//Massa por unidade de comprimento na configuração de referência
			(*Mr)(0, 0) = db.pipe_sections[section - 1]->Rho;
			(*Mr)(1, 1) = db.pipe_sections[section - 1]->Rho;
			(*Mr)(2, 2) = db.pipe_sections[section - 1]->Rho;
		}
		Dynamic* ptr_sol = static_cast<Dynamic*>(db.solution[db.current_solution_number - 1]);
		Mr->MatrixToPtr(pMr, 3);//salvando ptr para uso na dinâmica
		//Cálculo do integrando no ponto de Gauss
		EvaluateInertialContributions(temp_v, &ptr_sol->a1, &ptr_sol->a2, &ptr_sol->a3, &ptr_sol->a4,
			&ptr_sol->a5, &ptr_sol->a6, lag_save->alpha_i[gauss]->getMatrix(), alpha_delta[gauss]->getMatrix(),
			lag_save->u_i[gauss]->getMatrix(), u_delta[gauss]->getMatrix(), omega_i.getMatrix(), domega_i.getMatrix(), du_i.getMatrix(), ddu_i.getMatrix(), pJr, pMr, br->getMatrix(), dT->getMatrix(), pDdT);
		//Transformando o operador tangente em Matrix
		DdT->PtrToMatrix(pDdT, 6);
		(*inertial_loading) = (*inertial_loading) + (alpha1*jacobian)*transp(*N[gauss])*(*dT);
		(*mass) = (*mass) + (alpha1*jacobian)*transp(*N[gauss])*(*DdT)*(*N[gauss]);
	}//end of gauss
	//Conversão para o sistema global
	(*inertial_loading) = transp(*transform)*(*inertial_loading);
	(*mass) = (transp(*transform)*(*mass))*(*transform);
}

//Monta a matriz de amortecimento
void Pipe_1::MountDamping(bool update_rayleigh)
{
	if (typeid(*db.solution[db.current_solution_number - 1]) == typeid(Dynamic))
	{
		Dynamic* ptr_sol = static_cast<Dynamic*>(db.solution[db.current_solution_number - 1]);
		//Amortecimento de Rayleigh
		if (update_rayleigh == true)
		{
			MountMassModal();
			(*rayleigh_damping) = ptr_sol->alpha*(*mass_modal) + ptr_sol->beta*(*stiffness);
		}
		(*damping) = (*damping) + ptr_sol->a4*(*rayleigh_damping);	//Atualizando a matriz de amortecimento - inclusão do efeito de Rayleigh

		Matrix v_ipp(18);
		//Cálculo dos esforços de amortecimento - utilização de informação das velocidades nos GL do elemento
		for (int node = 0; node < 3; node++)
		{
			for (int GL = 0; GL<6; GL++)
				v_ipp(node * 6 + GL, 0) = db.nodes[nodes[node] - 1]->vel[GL];
		}
		(*damping_loading) = (*damping_loading) + (*rayleigh_damping)*v_ipp;
	}
}

//Composição da forma fraca e operador tangente - contribuições inerciais e amortecimento
void Pipe_1::MountDyn()
{
	(*P_loading) = (*P_loading) + (*inertial_loading) + (*damping_loading);
	(*i_loading) = (*i_loading) + (*inertial_loading) + (*damping_loading);
	(*stiffness) = (*stiffness) + (*mass) + (*damping);
}
//Montagens para análise modal - inserção da matriz de massa e amortecimento na matriz de rigidez para posterior montagem global
void Pipe_1::MountDynModal()
{
	(*stiffness) = (*mass_modal) + (*damping_modal);
}

//Calcula as contribuições inerciais para a análise dinâmica - inclui todas as contribuições para a forma fraca e para o operador tangente
void Pipe_1::EvaluateInertialContributions(double* v, double(*a1)
	, double(*a2), double(*a3), double(*a4), double(*a5), double(*a6)
	, double* alphai, double* alphad, double* ui, double* ud, double* omegai
	, double* domegai, double* dui, double* ddui, double** Jr, double** Mr
	, double* br, double* dT, double** DdT)
{
	v[604] = -((*a3)*ddui[2]) - (*a2)*dui[2] + (*a1)*ud[2];
	v[603] = -((*a3)*ddui[1]) - (*a2)*dui[1] + (*a1)*ud[1];
	v[602] = -((*a3)*ddui[0]) - (*a2)*dui[0] + (*a1)*ud[0];
	v[596] = (*a4)*alphad[2] + (*a6)*domegai[2] + (*a5)*omegai[2];
	v[595] = (*a1)*alphad[2] - (*a3)*domegai[2] - (*a2)*omegai[2];
	v[594] = (*a4)*alphad[1] + (*a6)*domegai[1] + (*a5)*omegai[1];
	v[593] = (*a1)*alphad[1] - (*a3)*domegai[1] - (*a2)*omegai[1];
	v[592] = (*a4)*alphad[0] + (*a6)*domegai[0] + (*a5)*omegai[0];
	v[591] = (*a1)*alphad[0] - (*a3)*domegai[0] - (*a2)*omegai[0];
	v[589] = Power(alphad[2], 2);
	v[588] = 0.5e0*alphad[2];
	v[587] = 2e0*alphad[2];
	v[586] = Power(alphad[1], 2);
	v[585] = 0.5e0*alphad[1];
	v[584] = 2e0*alphad[1];
	v[583] = Power(alphad[0], 2);
	v[582] = 2e0*alphad[0];
	v[581] = 0.5e0*alphad[0];
	v[580] = Power(alphai[2], 2);
	v[579] = 0.5e0*alphai[0] * alphai[2];
	v[578] = 0.5e0*alphai[1];
	v[577] = Power(alphai[1], 2);
	v[600] = v[577] + v[580];
	v[576] = alphai[0] * v[578];
	v[575] = Power(alphai[0], 2);
	v[132] = alphai[2] * v[578];
	v[109] = alphad[1] * v[581];
	v[288] = -v[583] - v[586];
	v[266] = alphad[2] + v[109];
	v[256] = -alphad[2] + v[109];
	v[116] = alphad[2] * v[585];
	v[283] = alphad[0] + v[116];
	v[275] = -alphad[0] + v[116];
	v[114] = alphad[2] * v[581];
	v[279] = -alphad[1] + v[114];
	v[262] = alphad[1] + v[114];
	v[270] = -v[583] - v[589];
	v[251] = -v[586] - v[589];
	v[246] = 4e0 + v[583] + v[586] + v[589];
	v[248] = 1e0 / Power(v[246], 2);
	v[590] = -4e0*v[248];
	v[250] = v[587] * v[590];
	v[599] = 0.5e0*v[250];
	v[292] = v[288] * v[599];
	v[249] = v[584] * v[590];
	v[597] = 0.5e0*v[249];
	v[477] = -(v[249] * v[588]);
	v[272] = v[270] * v[597];
	v[247] = v[582] * v[590];
	v[598] = 0.5e0*v[247];
	v[479] = v[247] * v[585];
	v[476] = -(v[247] * v[588]);
	v[252] = v[251] * v[598];
	v[103] = 4e0 / v[246];
	v[480] = -0.5e0*v[103];
	v[482] = v[480] - v[247] * v[581];
	v[481] = -v[480] + v[249] * v[585];
	v[478] = v[480] - v[250] * v[588];
	v[290] = v[480] * v[584];
	v[291] = v[290] + v[288] * v[597];
	v[287] = v[480] * v[582];
	v[289] = v[287] + v[288] * v[598];
	v[284] = v[103] + v[247] * v[283];
	v[281] = -v[103] + v[249] * v[279];
	v[276] = -v[103] + v[247] * v[275];
	v[273] = v[480] * v[587];
	v[274] = v[273] + v[270] * v[599];
	v[271] = v[287] + v[270] * v[598];
	v[269] = v[103] + v[250] * v[266];
	v[264] = v[103] + v[249] * v[262];
	v[261] = -(v[103] * v[588]);
	v[285] = -v[261] + v[249] * v[283];
	v[280] = -v[261] + v[247] * v[279];
	v[277] = -v[261] + v[249] * v[275];
	v[263] = -v[261] + v[247] * v[262];
	v[260] = -v[103] + v[250] * v[256];
	v[258] = -(v[103] * v[581]);
	v[282] = -v[258] + v[250] * v[279];
	v[268] = -v[258] + v[249] * v[266];
	v[265] = -v[258] + v[250] * v[262];
	v[259] = v[249] * v[256] - v[258];
	v[255] = v[103] * v[585];
	v[286] = v[255] + v[250] * v[283];
	v[278] = v[255] + v[250] * v[275];
	v[267] = v[255] + v[247] * v[266];
	v[257] = v[255] + v[247] * v[256];
	v[254] = v[273] + v[251] * v[599];
	v[253] = v[290] + v[251] * v[597];
	v[106] = 1e0 - v[251] * v[480];
	v[365] = (*a1)*v[106] + v[252] * v[591] + v[257] * v[593] + v[263] * v[595];
	v[356] = (*a4)*v[106] + v[252] * v[592] + v[257] * v[594] + v[263] * v[596];
	v[107] = v[103] * v[256];
	v[366] = (*a1)*v[107] + v[253] * v[591] + v[259] * v[593] + v[264] * v[595];
	v[357] = (*a4)*v[107] + v[253] * v[592] + v[259] * v[594] + v[264] * v[596];
	v[108] = v[103] * v[262];
	v[367] = (*a1)*v[108] + v[254] * v[591] + v[260] * v[593] + v[265] * v[595];
	v[358] = (*a4)*v[108] + v[254] * v[592] + v[260] * v[594] + v[265] * v[596];
	v[110] = v[103] * v[266];
	v[368] = (*a1)*v[110] + v[267] * v[591] + v[271] * v[593] + v[276] * v[595];
	v[359] = (*a4)*v[110] + v[267] * v[592] + v[271] * v[594] + v[276] * v[596];
	v[112] = 1e0 - v[270] * v[480];
	v[369] = (*a1)*v[112] + v[268] * v[591] + v[272] * v[593] + v[277] * v[595];
	v[360] = (*a4)*v[112] + v[268] * v[592] + v[272] * v[594] + v[277] * v[596];
	v[113] = v[103] * v[275];
	v[370] = (*a1)*v[113] + v[269] * v[591] + v[274] * v[593] + v[278] * v[595];
	v[361] = (*a4)*v[113] + v[269] * v[592] + v[274] * v[594] + v[278] * v[596];
	v[115] = v[103] * v[279];
	v[371] = (*a1)*v[115] + v[280] * v[591] + v[284] * v[593] + v[289] * v[595];
	v[362] = (*a4)*v[115] + v[280] * v[592] + v[284] * v[594] + v[289] * v[596];
	v[117] = v[103] * v[283];
	v[372] = (*a1)*v[117] + v[281] * v[591] + v[285] * v[593] + v[291] * v[595];
	v[363] = (*a4)*v[117] + v[281] * v[592] + v[285] * v[594] + v[291] * v[596];
	v[118] = 1e0 - v[288] * v[480];
	v[373] = (*a1)*v[118] + v[282] * v[591] + v[286] * v[593] + v[292] * v[595];
	v[364] = (*a4)*v[118] + v[282] * v[592] + v[286] * v[594] + v[292] * v[596];
	v[119] = 4e0 / (4e0 + v[575] + v[600]);
	v[601] = -0.5e0*v[119];
	v[122] = 1e0 + v[600] * v[601];
	v[123] = v[119] * (-alphai[2] + v[576]);
	v[124] = v[119] * (alphai[1] + v[579]);
	v[126] = v[119] * (alphai[2] + v[576]);
	v[128] = 1e0 + (v[575] + v[580])*v[601];
	v[129] = v[119] * (-alphai[0] + v[132]);
	v[131] = v[119] * (-alphai[1] + v[579]);
	v[331] = v[122] * v[282] + v[126] * v[286] + v[131] * v[292];
	v[330] = v[122] * v[281] + v[126] * v[285] + v[131] * v[291];
	v[329] = v[122] * v[280] + v[126] * v[284] + v[131] * v[289];
	v[313] = v[122] * v[269] + v[126] * v[274] + v[131] * v[278];
	v[312] = v[122] * v[268] + v[126] * v[272] + v[131] * v[277];
	v[311] = v[122] * v[267] + v[126] * v[271] + v[131] * v[276];
	v[295] = v[122] * v[254] + v[126] * v[260] + v[131] * v[265];
	v[294] = v[122] * v[253] + v[126] * v[259] + v[131] * v[264];
	v[293] = v[122] * v[252] + v[126] * v[257] + v[131] * v[263];
	v[133] = v[119] * (alphai[0] + v[132]);
	v[334] = v[123] * v[282] + v[128] * v[286] + v[133] * v[292];
	v[333] = v[123] * v[281] + v[128] * v[285] + v[133] * v[291];
	v[332] = v[123] * v[280] + v[128] * v[284] + v[133] * v[289];
	v[316] = v[123] * v[269] + v[128] * v[274] + v[133] * v[278];
	v[315] = v[123] * v[268] + v[128] * v[272] + v[133] * v[277];
	v[314] = v[123] * v[267] + v[128] * v[271] + v[133] * v[276];
	v[298] = v[123] * v[254] + v[128] * v[260] + v[133] * v[265];
	v[297] = v[123] * v[253] + v[128] * v[259] + v[133] * v[264];
	v[296] = v[123] * v[252] + v[128] * v[257] + v[133] * v[263];
	v[134] = 1e0 + (v[575] + v[577])*v[601];
	v[337] = v[124] * v[282] + v[129] * v[286] + v[134] * v[292];
	v[475] = Jr[0][0] * v[331] + Jr[1][0] * v[334] + Jr[2][0] * v[337];
	v[472] = Jr[0][1] * v[331] + Jr[1][1] * v[334] + Jr[2][1] * v[337];
	v[469] = Jr[0][2] * v[331] + Jr[1][2] * v[334] + Jr[2][2] * v[337];
	v[355] = br[0] * v[331] + br[1] * v[334] + br[2] * v[337];
	v[346] = Mr[0][0] * v[331] + Mr[1][0] * v[334] + Mr[2][0] * v[337];
	v[343] = Mr[0][1] * v[331] + Mr[1][1] * v[334] + Mr[2][1] * v[337];
	v[340] = Mr[0][2] * v[331] + Mr[1][2] * v[334] + Mr[2][2] * v[337];
	v[336] = v[124] * v[281] + v[129] * v[285] + v[134] * v[291];
	v[474] = Jr[0][0] * v[330] + Jr[1][0] * v[333] + Jr[2][0] * v[336];
	v[471] = Jr[0][1] * v[330] + Jr[1][1] * v[333] + Jr[2][1] * v[336];
	v[468] = Jr[0][2] * v[330] + Jr[1][2] * v[333] + Jr[2][2] * v[336];
	v[354] = br[0] * v[330] + br[1] * v[333] + br[2] * v[336];
	v[345] = Mr[0][0] * v[330] + Mr[1][0] * v[333] + Mr[2][0] * v[336];
	v[342] = Mr[0][1] * v[330] + Mr[1][1] * v[333] + Mr[2][1] * v[336];
	v[339] = Mr[0][2] * v[330] + Mr[1][2] * v[333] + Mr[2][2] * v[336];
	v[335] = v[124] * v[280] + v[129] * v[284] + v[134] * v[289];
	v[473] = Jr[0][0] * v[329] + Jr[1][0] * v[332] + Jr[2][0] * v[335];
	v[470] = Jr[0][1] * v[329] + Jr[1][1] * v[332] + Jr[2][1] * v[335];
	v[467] = Jr[0][2] * v[329] + Jr[1][2] * v[332] + Jr[2][2] * v[335];
	v[353] = br[0] * v[329] + br[1] * v[332] + br[2] * v[335];
	v[344] = Mr[0][0] * v[329] + Mr[1][0] * v[332] + Mr[2][0] * v[335];
	v[341] = Mr[0][1] * v[329] + Mr[1][1] * v[332] + Mr[2][1] * v[335];
	v[338] = Mr[0][2] * v[329] + Mr[1][2] * v[332] + Mr[2][2] * v[335];
	v[319] = v[124] * v[269] + v[129] * v[274] + v[134] * v[278];
	v[466] = Jr[0][0] * v[313] + Jr[1][0] * v[316] + Jr[2][0] * v[319];
	v[463] = Jr[0][1] * v[313] + Jr[1][1] * v[316] + Jr[2][1] * v[319];
	v[460] = Jr[0][2] * v[313] + Jr[1][2] * v[316] + Jr[2][2] * v[319];
	v[352] = br[0] * v[313] + br[1] * v[316] + br[2] * v[319];
	v[328] = Mr[0][0] * v[313] + Mr[1][0] * v[316] + Mr[2][0] * v[319];
	v[325] = Mr[0][1] * v[313] + Mr[1][1] * v[316] + Mr[2][1] * v[319];
	v[322] = Mr[0][2] * v[313] + Mr[1][2] * v[316] + Mr[2][2] * v[319];
	v[318] = v[124] * v[268] + v[129] * v[272] + v[134] * v[277];
	v[465] = Jr[0][0] * v[312] + Jr[1][0] * v[315] + Jr[2][0] * v[318];
	v[462] = Jr[0][1] * v[312] + Jr[1][1] * v[315] + Jr[2][1] * v[318];
	v[459] = Jr[0][2] * v[312] + Jr[1][2] * v[315] + Jr[2][2] * v[318];
	v[351] = br[0] * v[312] + br[1] * v[315] + br[2] * v[318];
	v[327] = Mr[0][0] * v[312] + Mr[1][0] * v[315] + Mr[2][0] * v[318];
	v[324] = Mr[0][1] * v[312] + Mr[1][1] * v[315] + Mr[2][1] * v[318];
	v[321] = Mr[0][2] * v[312] + Mr[1][2] * v[315] + Mr[2][2] * v[318];
	v[317] = v[124] * v[267] + v[129] * v[271] + v[134] * v[276];
	v[464] = Jr[0][0] * v[311] + Jr[1][0] * v[314] + Jr[2][0] * v[317];
	v[461] = Jr[0][1] * v[311] + Jr[1][1] * v[314] + Jr[2][1] * v[317];
	v[458] = Jr[0][2] * v[311] + Jr[1][2] * v[314] + Jr[2][2] * v[317];
	v[350] = br[0] * v[311] + br[1] * v[314] + br[2] * v[317];
	v[326] = Mr[0][0] * v[311] + Mr[1][0] * v[314] + Mr[2][0] * v[317];
	v[323] = Mr[0][1] * v[311] + Mr[1][1] * v[314] + Mr[2][1] * v[317];
	v[320] = Mr[0][2] * v[311] + Mr[1][2] * v[314] + Mr[2][2] * v[317];
	v[301] = v[124] * v[254] + v[129] * v[260] + v[134] * v[265];
	v[457] = Jr[0][0] * v[295] + Jr[1][0] * v[298] + Jr[2][0] * v[301];
	v[454] = Jr[0][1] * v[295] + Jr[1][1] * v[298] + Jr[2][1] * v[301];
	v[451] = Jr[0][2] * v[295] + Jr[1][2] * v[298] + Jr[2][2] * v[301];
	v[349] = br[0] * v[295] + br[1] * v[298] + br[2] * v[301];
	v[310] = Mr[0][0] * v[295] + Mr[1][0] * v[298] + Mr[2][0] * v[301];
	v[307] = Mr[0][1] * v[295] + Mr[1][1] * v[298] + Mr[2][1] * v[301];
	v[304] = Mr[0][2] * v[295] + Mr[1][2] * v[298] + Mr[2][2] * v[301];
	v[300] = v[124] * v[253] + v[129] * v[259] + v[134] * v[264];
	v[456] = Jr[0][0] * v[294] + Jr[1][0] * v[297] + Jr[2][0] * v[300];
	v[453] = Jr[0][1] * v[294] + Jr[1][1] * v[297] + Jr[2][1] * v[300];
	v[450] = Jr[0][2] * v[294] + Jr[1][2] * v[297] + Jr[2][2] * v[300];
	v[348] = br[0] * v[294] + br[1] * v[297] + br[2] * v[300];
	v[309] = Mr[0][0] * v[294] + Mr[1][0] * v[297] + Mr[2][0] * v[300];
	v[306] = Mr[0][1] * v[294] + Mr[1][1] * v[297] + Mr[2][1] * v[300];
	v[303] = Mr[0][2] * v[294] + Mr[1][2] * v[297] + Mr[2][2] * v[300];
	v[299] = v[124] * v[252] + v[129] * v[257] + v[134] * v[263];
	v[455] = Jr[0][0] * v[293] + Jr[1][0] * v[296] + Jr[2][0] * v[299];
	v[452] = Jr[0][1] * v[293] + Jr[1][1] * v[296] + Jr[2][1] * v[299];
	v[449] = Jr[0][2] * v[293] + Jr[1][2] * v[296] + Jr[2][2] * v[299];
	v[347] = br[0] * v[293] + br[1] * v[296] + br[2] * v[299];
	v[308] = Mr[0][0] * v[293] + Mr[1][0] * v[296] + Mr[2][0] * v[299];
	v[305] = Mr[0][1] * v[293] + Mr[1][1] * v[296] + Mr[2][1] * v[299];
	v[302] = Mr[0][2] * v[293] + Mr[1][2] * v[296] + Mr[2][2] * v[299];
	v[135] = v[106] * v[122] + v[107] * v[126] + v[108] * v[131];
	v[136] = v[106] * v[123] + v[107] * v[128] + v[108] * v[133];
	v[137] = v[106] * v[124] + v[107] * v[129] + v[108] * v[134];
	v[168] = Jr[0][2] * v[135] + Jr[1][2] * v[136] + Jr[2][2] * v[137];
	v[167] = Jr[0][1] * v[135] + Jr[1][1] * v[136] + Jr[2][1] * v[137];
	v[166] = Jr[0][0] * v[135] + Jr[1][0] * v[136] + Jr[2][0] * v[137];
	v[485] = v[166] * v[295] + v[167] * v[298] + v[168] * v[301] + v[137] * v[451] + v[136] * v[454] + v[135] * v[457];
	v[484] = v[166] * v[294] + v[167] * v[297] + v[168] * v[300] + v[137] * v[450] + v[136] * v[453] + v[135] * v[456];
	v[483] = v[166] * v[293] + v[167] * v[296] + v[168] * v[299] + v[137] * v[449] + v[136] * v[452] + v[135] * v[455];
	v[150] = Mr[0][2] * v[135] + Mr[1][2] * v[136] + Mr[2][2] * v[137];
	v[149] = Mr[0][1] * v[135] + Mr[1][1] * v[136] + Mr[2][1] * v[137];
	v[148] = Mr[0][0] * v[135] + Mr[1][0] * v[136] + Mr[2][0] * v[137];
	v[400] = v[148] * v[295] + v[149] * v[298] + v[150] * v[301] + v[137] * v[304] + v[136] * v[307] + v[135] * v[310];
	v[397] = v[148] * v[294] + v[149] * v[297] + v[150] * v[300] + v[137] * v[303] + v[136] * v[306] + v[135] * v[309];
	v[394] = v[148] * v[293] + v[149] * v[296] + v[150] * v[299] + v[137] * v[302] + v[136] * v[305] + v[135] * v[308];
	v[138] = v[110] * v[122] + v[112] * v[126] + v[113] * v[131];
	v[139] = v[110] * v[123] + v[112] * v[128] + v[113] * v[133];
	v[140] = v[110] * v[124] + v[112] * v[129] + v[113] * v[134];
	v[488] = v[166] * v[313] + v[167] * v[316] + v[168] * v[319] + v[140] * v[451] + v[139] * v[454] + v[138] * v[457];
	v[487] = v[166] * v[312] + v[167] * v[315] + v[168] * v[318] + v[140] * v[450] + v[139] * v[453] + v[138] * v[456];
	v[486] = v[166] * v[311] + v[167] * v[314] + v[168] * v[317] + v[140] * v[449] + v[139] * v[452] + v[138] * v[455];
	v[399] = v[140] * v[304] + v[139] * v[307] + v[138] * v[310] + v[148] * v[313] + v[149] * v[316] + v[150] * v[319];
	v[396] = v[140] * v[303] + v[139] * v[306] + v[138] * v[309] + v[148] * v[312] + v[149] * v[315] + v[150] * v[318];
	v[393] = v[140] * v[302] + v[139] * v[305] + v[138] * v[308] + v[148] * v[311] + v[149] * v[314] + v[150] * v[317];
	v[174] = Jr[0][2] * v[138] + Jr[1][2] * v[139] + Jr[2][2] * v[140];
	v[173] = Jr[0][1] * v[138] + Jr[1][1] * v[139] + Jr[2][1] * v[140];
	v[172] = Jr[0][0] * v[138] + Jr[1][0] * v[139] + Jr[2][0] * v[140];
	v[497] = v[172] * v[313] + v[173] * v[316] + v[174] * v[319] + v[140] * v[460] + v[139] * v[463] + v[138] * v[466];
	v[496] = v[172] * v[312] + v[173] * v[315] + v[174] * v[318] + v[140] * v[459] + v[139] * v[462] + v[138] * v[465];
	v[495] = v[172] * v[311] + v[173] * v[314] + v[174] * v[317] + v[140] * v[458] + v[139] * v[461] + v[138] * v[464];
	v[494] = v[172] * v[295] + v[173] * v[298] + v[174] * v[301] + v[137] * v[460] + v[136] * v[463] + v[135] * v[466];
	v[493] = v[172] * v[294] + v[173] * v[297] + v[174] * v[300] + v[137] * v[459] + v[136] * v[462] + v[135] * v[465];
	v[492] = v[172] * v[293] + v[173] * v[296] + v[174] * v[299] + v[137] * v[458] + v[136] * v[461] + v[135] * v[464];
	v[156] = Mr[0][2] * v[138] + Mr[1][2] * v[139] + Mr[2][2] * v[140];
	v[155] = Mr[0][1] * v[138] + Mr[1][1] * v[139] + Mr[2][1] * v[140];
	v[154] = Mr[0][0] * v[138] + Mr[1][0] * v[139] + Mr[2][0] * v[140];
	v[409] = v[154] * v[295] + v[155] * v[298] + v[156] * v[301] + v[137] * v[322] + v[136] * v[325] + v[135] * v[328];
	v[408] = v[154] * v[313] + v[155] * v[316] + v[156] * v[319] + v[140] * v[322] + v[139] * v[325] + v[138] * v[328];
	v[406] = v[154] * v[294] + v[155] * v[297] + v[156] * v[300] + v[137] * v[321] + v[136] * v[324] + v[135] * v[327];
	v[405] = v[154] * v[312] + v[155] * v[315] + v[156] * v[318] + v[140] * v[321] + v[139] * v[324] + v[138] * v[327];
	v[403] = v[154] * v[293] + v[155] * v[296] + v[156] * v[299] + v[137] * v[320] + v[136] * v[323] + v[135] * v[326];
	v[402] = v[154] * v[311] + v[155] * v[314] + v[156] * v[317] + v[140] * v[320] + v[139] * v[323] + v[138] * v[326];
	v[141] = v[115] * v[122] + v[117] * v[126] + v[118] * v[131];
	v[142] = v[115] * v[123] + v[117] * v[128] + v[118] * v[133];
	v[143] = v[115] * v[124] + v[117] * v[129] + v[118] * v[134];
	v[500] = v[172] * v[331] + v[173] * v[334] + v[174] * v[337] + v[143] * v[460] + v[142] * v[463] + v[141] * v[466];
	v[499] = v[172] * v[330] + v[173] * v[333] + v[174] * v[336] + v[143] * v[459] + v[142] * v[462] + v[141] * v[465];
	v[498] = v[172] * v[329] + v[173] * v[332] + v[174] * v[335] + v[143] * v[458] + v[142] * v[461] + v[141] * v[464];
	v[491] = v[166] * v[331] + v[167] * v[334] + v[168] * v[337] + v[143] * v[451] + v[142] * v[454] + v[141] * v[457];
	v[490] = v[166] * v[330] + v[167] * v[333] + v[168] * v[336] + v[143] * v[450] + v[142] * v[453] + v[141] * v[456];
	v[489] = v[166] * v[329] + v[167] * v[332] + v[168] * v[335] + v[143] * v[449] + v[142] * v[452] + v[141] * v[455];
	v[407] = v[143] * v[322] + v[142] * v[325] + v[141] * v[328] + v[154] * v[331] + v[155] * v[334] + v[156] * v[337];
	v[404] = v[143] * v[321] + v[142] * v[324] + v[141] * v[327] + v[154] * v[330] + v[155] * v[333] + v[156] * v[336];
	v[401] = v[143] * v[320] + v[142] * v[323] + v[141] * v[326] + v[154] * v[329] + v[155] * v[332] + v[156] * v[335];
	v[398] = v[143] * v[304] + v[142] * v[307] + v[141] * v[310] + v[148] * v[331] + v[149] * v[334] + v[150] * v[337];
	v[395] = v[143] * v[303] + v[142] * v[306] + v[141] * v[309] + v[148] * v[330] + v[149] * v[333] + v[150] * v[336];
	v[392] = v[143] * v[302] + v[142] * v[305] + v[141] * v[308] + v[148] * v[329] + v[149] * v[332] + v[150] * v[335];
	v[180] = Jr[0][2] * v[141] + Jr[1][2] * v[142] + Jr[2][2] * v[143];
	v[179] = Jr[0][1] * v[141] + Jr[1][1] * v[142] + Jr[2][1] * v[143];
	v[178] = Jr[0][0] * v[141] + Jr[1][0] * v[142] + Jr[2][0] * v[143];
	v[509] = v[178] * v[331] + v[179] * v[334] + v[180] * v[337] + v[143] * v[469] + v[142] * v[472] + v[141] * v[475];
	v[508] = v[178] * v[330] + v[179] * v[333] + v[180] * v[336] + v[143] * v[468] + v[142] * v[471] + v[141] * v[474];
	v[507] = v[178] * v[329] + v[179] * v[332] + v[180] * v[335] + v[143] * v[467] + v[142] * v[470] + v[141] * v[473];
	v[506] = v[178] * v[313] + v[179] * v[316] + v[180] * v[319] + v[140] * v[469] + v[139] * v[472] + v[138] * v[475];
	v[505] = v[178] * v[312] + v[179] * v[315] + v[180] * v[318] + v[140] * v[468] + v[139] * v[471] + v[138] * v[474];
	v[504] = v[178] * v[311] + v[179] * v[314] + v[180] * v[317] + v[140] * v[467] + v[139] * v[470] + v[138] * v[473];
	v[503] = v[178] * v[295] + v[179] * v[298] + v[180] * v[301] + v[137] * v[469] + v[136] * v[472] + v[135] * v[475];
	v[502] = v[178] * v[294] + v[179] * v[297] + v[180] * v[300] + v[137] * v[468] + v[136] * v[471] + v[135] * v[474];
	v[501] = v[178] * v[293] + v[179] * v[296] + v[180] * v[299] + v[137] * v[467] + v[136] * v[470] + v[135] * v[473];
	v[162] = Mr[0][2] * v[141] + Mr[1][2] * v[142] + Mr[2][2] * v[143];
	v[161] = Mr[0][1] * v[141] + Mr[1][1] * v[142] + Mr[2][1] * v[143];
	v[160] = Mr[0][0] * v[141] + Mr[1][0] * v[142] + Mr[2][0] * v[143];
	v[418] = v[160] * v[295] + v[161] * v[298] + v[162] * v[301] + v[137] * v[340] + v[136] * v[343] + v[135] * v[346];
	v[417] = v[160] * v[313] + v[161] * v[316] + v[162] * v[319] + v[140] * v[340] + v[139] * v[343] + v[138] * v[346];
	v[416] = v[160] * v[331] + v[161] * v[334] + v[162] * v[337] + v[143] * v[340] + v[142] * v[343] + v[141] * v[346];
	v[415] = v[160] * v[294] + v[161] * v[297] + v[162] * v[300] + v[137] * v[339] + v[136] * v[342] + v[135] * v[345];
	v[414] = v[160] * v[312] + v[161] * v[315] + v[162] * v[318] + v[140] * v[339] + v[139] * v[342] + v[138] * v[345];
	v[413] = v[160] * v[330] + v[161] * v[333] + v[162] * v[336] + v[143] * v[339] + v[142] * v[342] + v[141] * v[345];
	v[412] = v[160] * v[293] + v[161] * v[296] + v[162] * v[299] + v[137] * v[338] + v[136] * v[341] + v[135] * v[344];
	v[411] = v[160] * v[311] + v[161] * v[314] + v[162] * v[317] + v[140] * v[338] + v[139] * v[341] + v[138] * v[344];
	v[410] = v[160] * v[329] + v[161] * v[332] + v[162] * v[335] + v[143] * v[338] + v[142] * v[341] + v[141] * v[344];
	v[147] = v[135] * v[148] + v[136] * v[149] + v[137] * v[150];
	v[151] = v[138] * v[148] + v[139] * v[149] + v[140] * v[150];
	v[152] = v[141] * v[148] + v[142] * v[149] + v[143] * v[150];
	v[153] = v[135] * v[154] + v[136] * v[155] + v[137] * v[156];
	v[157] = v[138] * v[154] + v[139] * v[155] + v[140] * v[156];
	v[158] = v[141] * v[154] + v[142] * v[155] + v[143] * v[156];
	v[159] = v[135] * v[160] + v[136] * v[161] + v[137] * v[162];
	v[163] = v[138] * v[160] + v[139] * v[161] + v[140] * v[162];
	v[164] = v[141] * v[160] + v[142] * v[161] + v[143] * v[162];
	v[165] = v[135] * v[166] + v[136] * v[167] + v[137] * v[168];
	v[169] = v[138] * v[166] + v[139] * v[167] + v[140] * v[168];
	v[170] = v[141] * v[166] + v[142] * v[167] + v[143] * v[168];
	v[171] = v[135] * v[172] + v[136] * v[173] + v[137] * v[174];
	v[175] = v[138] * v[172] + v[139] * v[173] + v[140] * v[174];
	v[176] = v[141] * v[172] + v[142] * v[173] + v[143] * v[174];
	v[177] = v[135] * v[178] + v[136] * v[179] + v[137] * v[180];
	v[181] = v[138] * v[178] + v[139] * v[179] + v[140] * v[180];
	v[182] = v[141] * v[178] + v[142] * v[179] + v[143] * v[180];
	v[183] = br[0] * v[135] + br[1] * v[136] + br[2] * v[137];
	v[429] = (*a1)*v[183];
	v[184] = br[0] * v[138] + br[1] * v[139] + br[2] * v[140];
	v[428] = -((*a1)*v[184]);
	v[439] = -(v[147] * v[428]) - v[151] * v[429];
	v[436] = -(v[153] * v[428]) - v[157] * v[429];
	v[433] = -(v[159] * v[428]) - v[163] * v[429];
	v[185] = br[0] * v[141] + br[1] * v[142] + br[2] * v[143];
	v[430] = (*a1)*v[185];
	v[438] = v[152] * v[429] - v[147] * v[430];
	v[437] = v[152] * v[428] + v[151] * v[430];
	v[435] = v[158] * v[429] - v[153] * v[430];
	v[434] = v[158] * v[428] + v[157] * v[430];
	v[432] = v[164] * v[429] - v[159] * v[430];
	v[431] = v[164] * v[428] + v[163] * v[430];
	v[521] = -(v[352] * v[602]) + v[349] * v[603];
	v[520] = -(v[351] * v[602]) + v[348] * v[603];
	v[519] = -(v[350] * v[602]) + v[347] * v[603];
	v[527] = -(v[355] * v[603]) + v[352] * v[604];
	v[526] = -(v[354] * v[603]) + v[351] * v[604];
	v[525] = -(v[353] * v[603]) + v[350] * v[604];
	v[524] = v[355] * v[602] - v[349] * v[604];
	v[523] = v[354] * v[602] - v[348] * v[604];
	v[522] = v[353] * v[602] - v[347] * v[604];
	v[192] = v[106] * v[592] + v[107] * v[594] + v[108] * v[596];
	v[196] = v[110] * v[592] + v[112] * v[594] + v[113] * v[596];
	v[385] = -(v[196] * v[349]) + v[192] * v[352] + v[184] * v[358] - v[183] * v[361];
	v[384] = -(v[196] * v[348]) + v[192] * v[351] + v[184] * v[357] - v[183] * v[360];
	v[383] = -(v[196] * v[347]) + v[192] * v[350] + v[184] * v[356] - v[183] * v[359];
	v[197] = v[115] * v[592] + v[117] * v[594] + v[118] * v[596];
	v[548] = v[177] * v[358] + v[181] * v[361] + v[182] * v[364] + v[192] * v[503] + v[196] * v[506] + v[197] * v[509];
	v[547] = v[177] * v[357] + v[181] * v[360] + v[182] * v[363] + v[192] * v[502] + v[196] * v[505] + v[197] * v[508];
	v[546] = v[177] * v[356] + v[181] * v[359] + v[182] * v[362] + v[192] * v[501] + v[196] * v[504] + v[197] * v[507];
	v[542] = v[171] * v[358] + v[175] * v[361] + v[176] * v[364] + v[192] * v[494] + v[196] * v[497] + v[197] * v[500];
	v[541] = v[171] * v[357] + v[175] * v[360] + v[176] * v[363] + v[192] * v[493] + v[196] * v[496] + v[197] * v[499];
	v[540] = v[171] * v[356] + v[175] * v[359] + v[176] * v[362] + v[192] * v[492] + v[196] * v[495] + v[197] * v[498];
	v[539] = v[165] * v[358] + v[169] * v[361] + v[170] * v[364] + v[192] * v[485] + v[196] * v[488] + v[197] * v[491];
	v[538] = v[165] * v[357] + v[169] * v[360] + v[170] * v[363] + v[192] * v[484] + v[196] * v[487] + v[197] * v[490];
	v[537] = v[165] * v[356] + v[169] * v[359] + v[170] * v[362] + v[192] * v[483] + v[196] * v[486] + v[197] * v[489];
	v[379] = -(v[197] * v[352]) + v[196] * v[355] + v[185] * v[361] - v[184] * v[364];
	v[378] = -(v[197] * v[351]) + v[196] * v[354] + v[185] * v[360] - v[184] * v[363];
	v[377] = -(v[197] * v[350]) + v[196] * v[353] + v[185] * v[359] - v[184] * v[362];
	v[376] = v[197] * v[349] - v[192] * v[355] - v[185] * v[358] + v[183] * v[364];
	v[375] = v[197] * v[348] - v[192] * v[354] - v[185] * v[357] + v[183] * v[363];
	v[374] = v[197] * v[347] - v[192] * v[353] - v[185] * v[356] + v[183] * v[362];
	v[198] = v[106] * v[591] + v[107] * v[593] + v[108] * v[595];
	v[202] = v[110] * v[591] + v[112] * v[593] + v[113] * v[595];
	v[203] = v[115] * v[591] + v[117] * v[593] + v[118] * v[595];
	v[518] = v[165] * v[367] + v[169] * v[370] + v[170] * v[373] + v[198] * v[485] + v[202] * v[488] + v[203] * v[491];
	v[517] = v[165] * v[366] + v[169] * v[369] + v[170] * v[372] + v[198] * v[484] + v[202] * v[487] + v[203] * v[490];
	v[516] = v[165] * v[365] + v[169] * v[368] + v[170] * v[371] + v[198] * v[483] + v[202] * v[486] + v[203] * v[489];
	v[515] = v[171] * v[367] + v[175] * v[370] + v[176] * v[373] + v[198] * v[494] + v[202] * v[497] + v[203] * v[500];
	v[514] = v[171] * v[366] + v[175] * v[369] + v[176] * v[372] + v[198] * v[493] + v[202] * v[496] + v[203] * v[499];
	v[513] = v[171] * v[365] + v[175] * v[368] + v[176] * v[371] + v[198] * v[492] + v[202] * v[495] + v[203] * v[498];
	v[512] = v[177] * v[367] + v[181] * v[370] + v[182] * v[373] + v[198] * v[503] + v[202] * v[506] + v[203] * v[509];
	v[511] = v[177] * v[366] + v[181] * v[369] + v[182] * v[372] + v[198] * v[502] + v[202] * v[505] + v[203] * v[508];
	v[510] = v[177] * v[365] + v[181] * v[368] + v[182] * v[371] + v[198] * v[501] + v[202] * v[504] + v[203] * v[507];
	v[228] = v[177] * v[198] + v[181] * v[202] + v[182] * v[203];
	v[561] = -(v[228] * v[479]);
	v[227] = v[171] * v[198] + v[175] * v[202] + v[176] * v[203];
	v[568] = -(v[227] * v[476]);
	v[226] = v[165] * v[198] + v[169] * v[202] + v[170] * v[203];
	v[571] = v[226] * v[477];
	v[204] = -(v[185] * v[192]) + v[183] * v[197];
	v[205] = v[185] * v[196] - v[184] * v[197];
	v[382] = -(v[202] * v[349]) + v[198] * v[352] + v[204] * v[358] - v[205] * v[361] + v[184] * v[367] - v[183] * v[370]
		+ v[192] * v[376] - v[196] * v[379];
	v[381] = -(v[202] * v[348]) + v[198] * v[351] + v[204] * v[357] - v[205] * v[360] + v[184] * v[366] - v[183] * v[369]
		+ v[192] * v[375] - v[196] * v[378];
	v[380] = -(v[202] * v[347]) + v[198] * v[350] + v[204] * v[356] - v[205] * v[359] + v[184] * v[365] - v[183] * v[368]
		+ v[192] * v[374] - v[196] * v[377];
	v[208] = v[184] * v[198] - v[183] * v[202] + v[192] * v[204] - v[196] * v[205] + v[604];
	v[206] = v[184] * v[192] - v[183] * v[196];
	v[391] = v[203] * v[349] - v[198] * v[355] - v[206] * v[358] + v[205] * v[364] - v[185] * v[367] + v[183] * v[373]
		+ v[197] * v[379] - v[192] * v[385];
	v[390] = v[203] * v[348] - v[198] * v[354] - v[206] * v[357] + v[205] * v[363] - v[185] * v[366] + v[183] * v[372]
		+ v[197] * v[378] - v[192] * v[384];
	v[389] = v[203] * v[347] - v[198] * v[353] - v[206] * v[356] + v[205] * v[362] - v[185] * v[365] + v[183] * v[371]
		+ v[197] * v[377] - v[192] * v[383];
	v[388] = -(v[203] * v[352]) + v[202] * v[355] + v[206] * v[361] - v[204] * v[364] + v[185] * v[370] - v[184] * v[373]
		- v[197] * v[376] + v[196] * v[385];
	v[387] = -(v[203] * v[351]) + v[202] * v[354] + v[206] * v[360] - v[204] * v[363] + v[185] * v[369] - v[184] * v[372]
		- v[197] * v[375] + v[196] * v[384];
	v[386] = -(v[203] * v[350]) + v[202] * v[353] + v[206] * v[359] - v[204] * v[362] + v[185] * v[368] - v[184] * v[371]
		- v[197] * v[374] + v[196] * v[383];
	v[210] = v[185] * v[202] - v[184] * v[203] - v[197] * v[204] + v[196] * v[206] + v[602];
	v[209] = -(v[185] * v[198]) + v[183] * v[203] + v[197] * v[205] - v[192] * v[206] + v[603];
	v[213] = -(v[184] * v[602]) + v[183] * v[603];
	v[214] = v[185] * v[602] - v[183] * v[604];
	v[215] = -(v[185] * v[603]) + v[184] * v[604];
	v[536] = v[213] * v[398] + v[214] * v[399] + v[215] * v[400] + v[152] * v[521] + v[151] * v[524] + v[147] * v[527];
	v[535] = v[213] * v[395] + v[214] * v[396] + v[215] * v[397] + v[152] * v[520] + v[151] * v[523] + v[147] * v[526];
	v[534] = v[213] * v[392] + v[214] * v[393] + v[215] * v[394] + v[152] * v[519] + v[151] * v[522] + v[147] * v[525];
	v[533] = v[213] * v[407] + v[214] * v[408] + v[215] * v[409] + v[158] * v[521] + v[157] * v[524] + v[153] * v[527];
	v[532] = v[213] * v[404] + v[214] * v[405] + v[215] * v[406] + v[158] * v[520] + v[157] * v[523] + v[153] * v[526];
	v[531] = v[213] * v[401] + v[214] * v[402] + v[215] * v[403] + v[158] * v[519] + v[157] * v[522] + v[153] * v[525];
	v[530] = v[213] * v[416] + v[214] * v[417] + v[215] * v[418] + v[164] * v[521] + v[163] * v[524] + v[159] * v[527];
	v[529] = v[213] * v[413] + v[214] * v[414] + v[215] * v[415] + v[164] * v[520] + v[163] * v[523] + v[159] * v[526];
	v[528] = v[213] * v[410] + v[214] * v[411] + v[215] * v[412] + v[164] * v[519] + v[163] * v[522] + v[159] * v[525];
	v[222] = v[164] * v[213] + v[163] * v[214] + v[159] * v[215];
	v[559] = -(v[222] * v[479]);
	v[221] = v[158] * v[213] + v[157] * v[214] + v[153] * v[215];
	v[566] = -(v[221] * v[476]);
	v[220] = v[152] * v[213] + v[151] * v[214] + v[147] * v[215];
	v[569] = v[220] * v[477];
	v[216] = v[165] * v[192] + v[169] * v[196] + v[170] * v[197];
	v[217] = v[171] * v[192] + v[175] * v[196] + v[176] * v[197];
	v[545] = v[217] * v[358] - v[216] * v[361] - v[196] * v[539] + v[192] * v[542];
	v[614] = v[512] + v[530] + v[545];
	v[544] = v[217] * v[357] - v[216] * v[360] - v[196] * v[538] + v[192] * v[541];
	v[611] = v[511] + v[529] + v[544];
	v[543] = v[217] * v[356] - v[216] * v[359] - v[196] * v[537] + v[192] * v[540];
	v[608] = v[510] + v[528] + v[543];
	v[223] = -(v[196] * v[216]) + v[192] * v[217];
	v[606] = v[222] + v[223] + v[228];
	v[560] = -(v[223] * v[479]);
	v[218] = v[177] * v[192] + v[181] * v[196] + v[182] * v[197];
	v[554] = -(v[218] * v[358]) + v[216] * v[364] + v[197] * v[539] - v[192] * v[548];
	v[616] = v[515] + v[533] + v[554];
	v[553] = -(v[218] * v[357]) + v[216] * v[363] + v[197] * v[538] - v[192] * v[547];
	v[613] = v[514] + v[532] + v[553];
	v[552] = -(v[218] * v[356]) + v[216] * v[362] + v[197] * v[537] - v[192] * v[546];
	v[610] = v[513] + v[531] + v[552];
	v[551] = v[218] * v[361] - v[217] * v[364] - v[197] * v[542] + v[196] * v[548];
	v[615] = v[518] + v[536] + v[551];
	v[550] = v[218] * v[360] - v[217] * v[363] - v[197] * v[541] + v[196] * v[547];
	v[612] = v[517] + v[535] + v[550];
	v[549] = v[218] * v[359] - v[217] * v[362] - v[197] * v[540] + v[196] * v[546];
	v[609] = v[516] + v[534] + v[549];
	v[225] = -(v[197] * v[217]) + v[196] * v[218];
	v[605] = v[220] + v[225] + v[226];
	v[570] = v[225] * v[477];
	v[224] = v[197] * v[216] - v[192] * v[218];
	v[607] = v[221] + v[224] + v[227];
	v[567] = -(v[224] * v[476]);
	dT[0] = v[152] * v[208] + v[151] * v[209] + v[147] * v[210];
	dT[1] = v[158] * v[208] + v[157] * v[209] + v[153] * v[210];
	dT[2] = v[164] * v[208] + v[163] * v[209] + v[159] * v[210];
	dT[3] = v[103] * v[605] - v[255] * v[606] - v[261] * v[607];
	dT[4] = v[261] * v[605] - v[258] * v[606] + v[103] * v[607];
	dT[5] = v[255] * v[605] + v[103] * v[606] + v[258] * v[607];
	DdT[0][0] = (*a1)*v[147];
	DdT[0][1] = (*a1)*v[151];
	DdT[0][2] = (*a1)*v[152];
	DdT[0][3] = v[152] * v[380] + v[147] * v[386] + v[151] * v[389] + v[208] * v[392] + v[209] * v[393] + v[210] * v[394];
	DdT[0][4] = v[152] * v[381] + v[147] * v[387] + v[151] * v[390] + v[208] * v[395] + v[209] * v[396] + v[210] * v[397];
	DdT[0][5] = v[152] * v[382] + v[147] * v[388] + v[151] * v[391] + v[208] * v[398] + v[209] * v[399] + v[210] * v[400];
	DdT[1][0] = (*a1)*v[153];
	DdT[1][1] = (*a1)*v[157];
	DdT[1][2] = (*a1)*v[158];
	DdT[1][3] = v[158] * v[380] + v[153] * v[386] + v[157] * v[389] + v[208] * v[401] + v[209] * v[402] + v[210] * v[403];
	DdT[1][4] = v[158] * v[381] + v[153] * v[387] + v[157] * v[390] + v[208] * v[404] + v[209] * v[405] + v[210] * v[406];
	DdT[1][5] = v[158] * v[382] + v[153] * v[388] + v[157] * v[391] + v[208] * v[407] + v[209] * v[408] + v[210] * v[409];
	DdT[2][0] = (*a1)*v[159];
	DdT[2][1] = (*a1)*v[163];
	DdT[2][2] = (*a1)*v[164];
	DdT[2][3] = v[164] * v[380] + v[159] * v[386] + v[163] * v[389] + v[208] * v[410] + v[209] * v[411] + v[210] * v[412];
	DdT[2][4] = v[164] * v[381] + v[159] * v[387] + v[163] * v[390] + v[208] * v[413] + v[209] * v[414] + v[210] * v[415];
	DdT[2][5] = v[164] * v[382] + v[159] * v[388] + v[163] * v[391] + v[208] * v[416] + v[209] * v[417] + v[210] * v[418];
	DdT[3][0] = -(v[255] * v[431]) - v[261] * v[434] + v[103] * v[437];
	DdT[3][1] = -(v[255] * v[432]) - v[261] * v[435] + v[103] * v[438];
	DdT[3][2] = -(v[255] * v[433]) - v[261] * v[436] + v[103] * v[439];
	DdT[3][3] = v[559] + v[560] + v[561] + v[566] + v[567] + v[568] + v[247] * v[605] - v[255] * v[608] + v[103] * v[609]
		- v[261] * v[610];
	DdT[3][4] = v[249] * v[605] - v[481] * v[606] - v[477] * v[607] - v[255] * v[611] + v[103] * v[612] - v[261] * v[613];
	DdT[3][5] = v[250] * v[605] + v[477] * v[606] - v[478] * v[607] - v[255] * v[614] + v[103] * v[615] - v[261] * v[616];
	DdT[4][0] = -(v[258] * v[431]) + v[103] * v[434] + v[261] * v[437];
	DdT[4][1] = -(v[258] * v[432]) + v[103] * v[435] + v[261] * v[438];
	DdT[4][2] = -(v[258] * v[433]) + v[103] * v[436] + v[261] * v[439];
	DdT[4][3] = v[476] * v[605] - v[482] * v[606] + v[247] * v[607] - v[258] * v[608] + v[261] * v[609] + v[103] * v[610];
	DdT[4][4] = -v[559] - v[560] - v[561] + v[569] + v[570] + v[571] + v[249] * v[607] - v[258] * v[611] + v[261] * v[612]
		+ v[103] * v[613];
	DdT[4][5] = v[478] * v[605] - v[476] * v[606] + v[250] * v[607] - v[258] * v[614] + v[261] * v[615] + v[103] * v[616];
	DdT[5][0] = v[103] * v[431] + v[258] * v[434] + v[255] * v[437];
	DdT[5][1] = v[103] * v[432] + v[258] * v[435] + v[255] * v[438];
	DdT[5][2] = v[103] * v[433] + v[258] * v[436] + v[255] * v[439];
	DdT[5][3] = v[479] * v[605] + v[247] * v[606] + v[482] * v[607] + v[103] * v[608] + v[255] * v[609] + v[258] * v[610];
	DdT[5][4] = v[481] * v[605] + v[249] * v[606] - v[479] * v[607] + v[103] * v[611] + v[255] * v[612] + v[258] * v[613];
	DdT[5][5] = -v[566] - v[567] - v[568] - v[569] - v[570] - v[571] + v[250] * v[606] + v[103] * v[614] + v[255] * v[615]
		+ v[258] * v[616];
};

void Pipe_1::EvaluateMorisonContributions(double* v, double(*a4), double
	(*a5), double(*a6), double(*C1t), double(*C1n), double* dU, double* e3r
	, double* alphai, double* alphad, double* ud, double* dui, double* ddui
	, double* force, double** stiffness)
{
	//int b135, b140;
	v[362] = -((*a6)*ddui[2]) + dU[2] - (*a5)*dui[2] - (*a4)*ud[2];
	v[361] = -((*a6)*ddui[1]) + dU[1] - (*a5)*dui[1] - (*a4)*ud[1];
	v[360] = -((*a6)*ddui[0]) + dU[0] - (*a5)*dui[0] - (*a4)*ud[0];
	v[352] = Power(alphad[2], 2);
	v[351] = 0.5e0*alphad[0] * alphad[2];
	v[350] = 0.5e0*alphad[1];
	v[349] = 2e0*alphad[2];
	v[348] = Power(alphad[1], 2);
	v[347] = alphad[0] * v[350];
	v[346] = 2e0*alphad[1];
	v[345] = Power(alphad[0], 2);
	v[344] = 2e0*alphad[0];
	v[343] = Power(alphai[2], 2);
	v[342] = 0.5e0*alphai[0] * alphai[2];
	v[341] = 0.5e0*alphai[1];
	v[340] = Power(alphai[1], 2);
	v[358] = v[340] + v[343];
	v[339] = alphai[0] * v[341];
	v[338] = Power(alphai[0], 2);
	v[107] = alphai[2] * v[341];
	v[248] = -v[345] - v[348];
	v[225] = alphad[2] + v[347];
	v[215] = -alphad[2] + v[347];
	v[91] = alphad[2] * v[350];
	v[243] = alphad[0] + v[91];
	v[234] = -alphad[0] + v[91];
	v[239] = -alphad[1] + v[351];
	v[219] = alphad[1] + v[351];
	v[229] = -v[345] - v[352];
	v[211] = -v[348] - v[352];
	v[206] = 4e0 + v[345] + v[348] + v[352];
	v[208] = 1e0 / Power(v[206], 2);
	v[353] = -4e0*v[208];
	v[210] = v[349] * v[353];
	v[357] = 0.5e0*v[210];
	v[252] = v[248] * v[357];
	v[209] = v[346] * v[353];
	v[354] = 0.5e0*v[209];
	v[231] = v[229] * v[354];
	v[207] = v[344] * v[353];
	v[356] = 0.5e0*v[207];
	v[212] = v[211] * v[356];
	v[78] = 4e0 / v[206];
	v[355] = -0.5e0*v[78];
	v[250] = v[346] * v[355];
	v[251] = v[250] + v[248] * v[354];
	v[247] = v[344] * v[355];
	v[249] = v[247] + v[248] * v[356];
	v[244] = v[207] * v[243] + v[78];
	v[241] = v[209] * v[239] - v[78];
	v[236] = -(alphad[2] * v[355]);
	v[245] = v[236] + v[209] * v[243];
	v[240] = v[236] + v[207] * v[239];
	v[237] = v[209] * v[234] + v[236];
	v[235] = v[207] * v[234] - v[78];
	v[232] = v[349] * v[355];
	v[233] = v[232] + v[229] * v[357];
	v[230] = v[247] + v[229] * v[356];
	v[228] = v[210] * v[225] + v[78];
	v[224] = -(alphad[1] * v[355]);
	v[246] = v[224] + v[210] * v[243];
	v[238] = v[224] + v[210] * v[234];
	v[226] = v[224] + v[207] * v[225];
	v[222] = -(alphad[0] * v[355]);
	v[242] = v[222] + v[210] * v[239];
	v[227] = v[222] + v[209] * v[225];
	v[223] = v[210] * v[219] + v[222];
	v[221] = v[209] * v[219] + v[78];
	v[220] = v[207] * v[219] + v[236];
	v[218] = v[210] * v[215] - v[78];
	v[217] = v[209] * v[215] + v[222];
	v[216] = v[207] * v[215] + v[224];
	v[214] = v[232] + v[211] * v[357];
	v[213] = v[250] + v[211] * v[354];
	v[81] = 1e0 - v[211] * v[355];
	v[82] = v[215] * v[78];
	v[83] = v[219] * v[78];
	v[85] = v[225] * v[78];
	v[87] = 1e0 - v[229] * v[355];
	v[88] = v[234] * v[78];
	v[90] = v[239] * v[78];
	v[92] = v[243] * v[78];
	v[93] = 1e0 - v[248] * v[355];
	v[94] = 4e0 / (4e0 + v[338] + v[358]);
	v[359] = -0.5e0*v[94];
	v[97] = 1e0 + v[358] * v[359];
	v[98] = (-alphai[2] + v[339])*v[94];
	v[99] = (alphai[1] + v[342])*v[94];
	v[101] = (alphai[2] + v[339])*v[94];
	v[103] = 1e0 + (v[338] + v[343])*v[359];
	v[104] = (-alphai[0] + v[107])*v[94];
	v[106] = (-alphai[1] + v[342])*v[94];
	v[108] = (alphai[0] + v[107])*v[94];
	v[109] = 1e0 + (v[338] + v[340])*v[359];
	v[288] = e3r[0] * (v[101] * v[246] + v[106] * v[252] + v[242] * v[97]) + e3r[1] * (v[103] * v[246] + v[108] * v[252]
		+ v[242] * v[98]) + e3r[2] * (v[104] * v[246] + v[109] * v[252] + v[242] * v[99]);
	v[284] = e3r[0] * (v[101] * v[245] + v[106] * v[251] + v[241] * v[97]) + e3r[1] * (v[103] * v[245] + v[108] * v[251]
		+ v[241] * v[98]) + e3r[2] * (v[104] * v[245] + v[109] * v[251] + v[241] * v[99]);
	v[280] = e3r[0] * (v[101] * v[244] + v[106] * v[249] + v[240] * v[97]) + e3r[1] * (v[103] * v[244] + v[108] * v[249]
		+ v[240] * v[98]) + e3r[2] * (v[104] * v[244] + v[109] * v[249] + v[240] * v[99]);
	v[276] = e3r[0] * (v[101] * v[233] + v[106] * v[238] + v[228] * v[97]) + e3r[1] * (v[103] * v[233] + v[108] * v[238]
		+ v[228] * v[98]) + e3r[2] * (v[104] * v[233] + v[109] * v[238] + v[228] * v[99]);
	v[272] = e3r[0] * (v[101] * v[231] + v[106] * v[237] + v[227] * v[97]) + e3r[1] * (v[103] * v[231] + v[108] * v[237]
		+ v[227] * v[98]) + e3r[2] * (v[104] * v[231] + v[109] * v[237] + v[227] * v[99]);
	v[268] = e3r[0] * (v[101] * v[230] + v[106] * v[235] + v[226] * v[97]) + e3r[1] * (v[103] * v[230] + v[108] * v[235]
		+ v[226] * v[98]) + e3r[2] * (v[104] * v[230] + v[109] * v[235] + v[226] * v[99]);
	v[264] = e3r[0] * (v[101] * v[218] + v[106] * v[223] + v[214] * v[97]) + e3r[1] * (v[103] * v[218] + v[108] * v[223]
		+ v[214] * v[98]) + e3r[2] * (v[104] * v[218] + v[109] * v[223] + v[214] * v[99]);
	v[260] = e3r[0] * (v[101] * v[217] + v[106] * v[221] + v[213] * v[97]) + e3r[1] * (v[103] * v[217] + v[108] * v[221]
		+ v[213] * v[98]) + e3r[2] * (v[104] * v[217] + v[109] * v[221] + v[213] * v[99]);
	v[256] = e3r[0] * (v[101] * v[216] + v[106] * v[220] + v[212] * v[97]) + e3r[1] * (v[103] * v[216] + v[108] * v[220]
		+ v[212] * v[98]) + e3r[2] * (v[104] * v[216] + v[109] * v[220] + v[212] * v[99]);
	v[122] = e3r[0] * (v[101] * v[82] + v[106] * v[83] + v[81] * v[97]) + e3r[1] * (v[103] * v[82] + v[108] * v[83] + v[81] * v[98]
		) + e3r[2] * (v[104] * v[82] + v[109] * v[83] + v[81] * v[99]);
	v[151] = -((*a4)*(v[122] * v[122]));
	v[157] = -(*a4) - v[151];
	v[123] = e3r[0] * (v[101] * v[87] + v[106] * v[88] + v[85] * v[97]) + e3r[1] * (v[103] * v[87] + v[108] * v[88] + v[85] * v[98]
		) + e3r[2] * (v[104] * v[87] + v[109] * v[88] + v[85] * v[99]);
	v[149] = -((*a4)*v[123]);
	v[154] = v[123] * v[149];
	v[158] = -(*a4) - v[154];
	v[152] = v[122] * v[149];
	v[124] = e3r[0] * (v[101] * v[92] + v[106] * v[93] + v[90] * v[97]) + e3r[1] * (v[103] * v[92] + v[108] * v[93] + v[90] * v[98]
		) + e3r[2] * (v[104] * v[92] + v[109] * v[93] + v[90] * v[99]);
	v[150] = -((*a4)*v[124]);
	v[156] = v[124] * v[150];
	v[159] = -(*a4) - v[156];
	v[155] = v[123] * v[150];
	v[153] = v[122] * v[150];
	v[291] = v[264] * v[360] + v[276] * v[361] + v[288] * v[362];
	v[290] = v[260] * v[360] + v[272] * v[361] + v[284] * v[362];
	v[289] = v[256] * v[360] + v[268] * v[361] + v[280] * v[362];
	v[129] = v[122] * v[360] + v[123] * v[361] + v[124] * v[362];
	v[300] = v[129] * v[288] + v[124] * v[291];
	v[299] = v[129] * v[284] + v[124] * v[290];
	v[298] = v[129] * v[280] + v[124] * v[289];
	v[297] = v[129] * v[276] + v[123] * v[291];
	v[296] = v[129] * v[272] + v[123] * v[290];
	v[295] = v[129] * v[268] + v[123] * v[289];
	v[294] = v[129] * v[264] + v[122] * v[291];
	v[293] = v[129] * v[260] + v[122] * v[290];
	v[292] = v[129] * v[256] + v[122] * v[289];
	v[128] = v[122] * v[129];
	v[130] = v[123] * v[129];
	v[131] = v[124] * v[129];
	v[132] = -v[128] + v[360];
	v[133] = -v[130] + v[361];
	v[134] = -v[131] + v[362];
	if (sqrt(Power(v[128], 2) + Power(v[130], 2) + Power(v[131], 2)) == 0e0){
		v[160] = 0e0;
		v[161] = 0e0;
		v[162] = 0e0;
		v[301] = 0e0;
		v[302] = 0e0;
		v[303] = 0e0;
		v[136] = 0e0;
		v[163] = 0e0;
		v[164] = 0e0;
		v[165] = 0e0;
		v[304] = 0e0;
		v[305] = 0e0;
		v[306] = 0e0;
		v[137] = 0e0;
		v[166] = 0e0;
		v[167] = 0e0;
		v[168] = 0e0;
		v[307] = 0e0;
		v[308] = 0e0;
		v[309] = 0e0;
		v[138] = 0e0;
	}
	else {
		v[363] = sqrt(Power(v[128], 2) + Power(v[130], 2) + Power(v[131], 2));
		v[364] = (*C1t)*v[363];
		v[310] = 1e0 / v[363];
		v[313] = (v[128] * v[294] + v[130] * v[297] + v[131] * v[300])*v[310];
		v[312] = (v[128] * v[293] + v[130] * v[296] + v[131] * v[299])*v[310];
		v[311] = (v[128] * v[292] + v[130] * v[295] + v[131] * v[298])*v[310];
		v[173] = (v[128] * v[153] + v[130] * v[155] + v[131] * v[156])*v[310];
		v[172] = (v[128] * v[152] + v[130] * v[154] + v[131] * v[155])*v[310];
		v[160] = (*C1t)*(v[128] * (v[128] * v[151] + v[130] * v[152] + v[131] * v[153])*v[310] + v[151] * v[363]);
		v[161] = (*C1t)*(v[128] * v[172] + v[152] * v[363]);
		v[162] = (*C1t)*(v[128] * v[173] + v[153] * v[363]);
		v[301] = (*C1t)*(v[128] * v[311] + v[292] * v[363]);
		v[302] = (*C1t)*(v[128] * v[312] + v[293] * v[363]);
		v[303] = (*C1t)*(v[128] * v[313] + v[294] * v[363]);
		v[136] = v[128] * v[364];
		v[163] = v[161];
		v[164] = (*C1t)*(v[130] * v[172] + v[154] * v[363]);
		v[165] = (*C1t)*(v[130] * v[173] + v[155] * v[363]);
		v[304] = (*C1t)*(v[130] * v[311] + v[295] * v[363]);
		v[305] = (*C1t)*(v[130] * v[312] + v[296] * v[363]);
		v[306] = (*C1t)*(v[130] * v[313] + v[297] * v[363]);
		v[137] = v[130] * v[364];
		v[166] = v[162];
		v[167] = v[165];
		v[168] = (*C1t)*(v[131] * v[173] + v[156] * v[363]);
		v[307] = (*C1t)*(v[131] * v[311] + v[298] * v[363]);
		v[308] = (*C1t)*(v[131] * v[312] + v[299] * v[363]);
		v[309] = (*C1t)*(v[131] * v[313] + v[300] * v[363]);
		v[138] = v[131] * v[364];
	};
	if (sqrt(Power(v[132], 2) + Power(v[133], 2) + Power(v[134], 2)) == 0e0){
		v[174] = 0e0;
		v[175] = 0e0;
		v[176] = 0e0;
		v[314] = 0e0;
		v[315] = 0e0;
		v[316] = 0e0;
		v[141] = 0e0;
		v[177] = 0e0;
		v[178] = 0e0;
		v[179] = 0e0;
		v[317] = 0e0;
		v[318] = 0e0;
		v[319] = 0e0;
		v[142] = 0e0;
		v[180] = 0e0;
		v[181] = 0e0;
		v[182] = 0e0;
		v[320] = 0e0;
		v[321] = 0e0;
		v[322] = 0e0;
		v[143] = 0e0;
	}
	else {
		v[369] = (*C1n)*v[134];
		v[368] = (*C1n)*v[133];
		v[367] = (*C1n)*v[132];
		v[365] = sqrt(Power(v[132], 2) + Power(v[133], 2) + Power(v[134], 2));
		v[366] = -((*C1n)*v[365]);
		v[323] = 1e0 / v[365];
		v[326] = -((v[132] * v[294] + v[133] * v[297] + v[134] * v[300])*v[323]);
		v[325] = -((v[132] * v[293] + v[133] * v[296] + v[134] * v[299])*v[323]);
		v[324] = -((v[132] * v[292] + v[133] * v[295] + v[134] * v[298])*v[323]);
		v[187] = -((v[132] * v[153] + v[133] * v[155] - v[134] * v[159])*v[323]);
		v[186] = -((v[132] * v[152] + v[134] * v[155] - v[133] * v[158])*v[323]);
		v[184] = -((v[133] * v[152] + v[134] * v[153] - v[132] * v[157])*v[323]);
		v[190] = v[155] * v[366];
		v[189] = v[153] * v[366];
		v[188] = v[152] * v[366];
		v[174] = (*C1n)*(v[132] * v[184] + v[157] * v[365]);
		v[175] = v[188] + v[186] * v[367];
		v[176] = v[189] + v[187] * v[367];
		v[314] = (*C1n)*(v[132] * v[324] - v[292] * v[365]);
		v[315] = (*C1n)*(v[132] * v[325] - v[293] * v[365]);
		v[316] = (*C1n)*(v[132] * v[326] - v[294] * v[365]);
		v[141] = -(v[132] * v[366]);
		v[177] = v[188] + v[184] * v[368];
		v[178] = (*C1n)*(v[133] * v[186] + v[158] * v[365]);
		v[179] = v[190] + v[187] * v[368];
		v[317] = (*C1n)*(v[133] * v[324] - v[295] * v[365]);
		v[318] = (*C1n)*(v[133] * v[325] - v[296] * v[365]);
		v[319] = (*C1n)*(v[133] * v[326] - v[297] * v[365]);
		v[142] = -(v[133] * v[366]);
		v[180] = v[189] + v[184] * v[369];
		v[181] = v[190] + v[186] * v[369];
		v[182] = (*C1n)*(v[134] * v[187] + v[159] * v[365]);
		v[320] = (*C1n)*(v[134] * v[324] - v[298] * v[365]);
		v[321] = (*C1n)*(v[134] * v[325] - v[299] * v[365]);
		v[322] = (*C1n)*(v[134] * v[326] - v[300] * v[365]);
		v[143] = -(v[134] * v[366]);
	};
	force[0] = v[136] + v[141];
	force[1] = v[137] + v[142];
	force[2] = v[138] + v[143];
	force[3] = 0e0;
	force[4] = 0e0;
	force[5] = 0e0;
	stiffness[0][0] = v[160] + v[174];
	stiffness[0][1] = v[161] + v[175];
	stiffness[0][2] = v[162] + v[176];
	stiffness[0][3] = v[301] + v[314];
	stiffness[0][4] = v[302] + v[315];
	stiffness[0][5] = v[303] + v[316];
	stiffness[1][0] = v[163] + v[177];
	stiffness[1][1] = v[164] + v[178];
	stiffness[1][2] = v[165] + v[179];
	stiffness[1][3] = v[304] + v[317];
	stiffness[1][4] = v[305] + v[318];
	stiffness[1][5] = v[306] + v[319];
	stiffness[2][0] = v[166] + v[180];
	stiffness[2][1] = v[167] + v[181];
	stiffness[2][2] = v[168] + v[182];
	stiffness[2][3] = v[307] + v[320];
	stiffness[2][4] = v[308] + v[321];
	stiffness[2][5] = v[309] + v[322];
	stiffness[3][0] = 0e0;
	stiffness[3][1] = 0e0;
	stiffness[3][2] = 0e0;
	stiffness[3][3] = 0e0;
	stiffness[3][4] = 0e0;
	stiffness[3][5] = 0e0;
	stiffness[4][0] = 0e0;
	stiffness[4][1] = 0e0;
	stiffness[4][2] = 0e0;
	stiffness[4][3] = 0e0;
	stiffness[4][4] = 0e0;
	stiffness[4][5] = 0e0;
	stiffness[5][0] = 0e0;
	stiffness[5][1] = 0e0;
	stiffness[5][2] = 0e0;
	stiffness[5][3] = 0e0;
	stiffness[5][4] = 0e0;
	stiffness[5][5] = 0e0;
};


