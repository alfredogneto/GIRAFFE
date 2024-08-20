#include "Shell_1.h"

#include"Database.h"
//Variáveis globais
extern
Database db;

Shell_1::Shell_1()
{
	strain_energy = 0.0;
	kinetic_energy = 0.0;
	potential_gravitational_energy = 0.0;

	tempkin = new double;

	VTK_type = 22;
	nDOFs = 27;
	//CORNERS são os nós 1,2,3 e MIDS são os nós 4,5,6
	material = 0;
	section = 0;
	n_nodes = 6;
	number = 0;
	nodes = new int[n_nodes];
	VTK_nodes = new int[n_nodes];
	DOFs = new int *[n_nodes];
	for (int i = 0; i < n_nodes; i++)
		DOFs[i] = new int[db.number_GLs_node];

	VTK_nodes[0] = 0;
	VTK_nodes[1] = 1;
	VTK_nodes[2] = 2;
	VTK_nodes[3] = 3;
	VTK_nodes[4] = 4;
	VTK_nodes[5] = 5;
	//Rotina para ativar os GLS de cada nó do elemento
	for (int i = 0; i < n_nodes; i++)
	{
		for (int j = 0; j < db.number_GLs_node; j++)
		{
			DOFs[i][j] = 0;
		}
		//Só translação
		DOFs[i][0] = 1;
		DOFs[i][1] = 1;
		DOFs[i][2] = 1;
		//Rotação
		if (i > 2)
		{
			DOFs[i][3] = 1;
			DOFs[i][4] = 1;
			DOFs[i][5] = 1;
		}
	}
	type_name = new char[20];//Nome do tipo do elemento
	sprintf(type_name, "Shell_1");

	//Alocações de memória
	//Cada variável carregará a informação para seu ponto de Gauss (3 pontos de Gauss)
	N1u = new double[3];
	N2u = new double[3];
	N3u = new double[3];
	N4u = new double[3];
	N5u = new double[3];
	N6u = new double[3];
	N4a = new double[3];
	N5a = new double[3];
	N6a = new double[3];

	N1u_x1 = new double[3];
	N2u_x1 = new double[3];
	N3u_x1 = new double[3];
	N4u_x1 = new double[3];
	N5u_x1 = new double[3];
	N6u_x1 = new double[3];
	N4a_x1 = new double[3];
	N5a_x1 = new double[3];
	N6a_x1 = new double[3];

	N1u_x2 = new double[3];
	N2u_x2 = new double[3];
	N3u_x2 = new double[3];
	N4u_x2 = new double[3];
	N5u_x2 = new double[3];
	N6u_x2 = new double[3];
	N4a_x2 = new double[3];
	N5a_x2 = new double[3];
	N6a_x2 = new double[3];

	N = new Matrix*[3];
	deltaN = new Matrix*[3];
	alpha_delta = new Matrix*[3];
	alpha_delta_x1 = new Matrix*[3];
	alpha_delta_x2 = new Matrix*[3];
	u_delta = new Matrix*[3];
	u_delta_x1 = new Matrix*[3];
	u_delta_x2 = new Matrix*[3];

	Q_i = new Matrix*[3];
	Q_delta = new Matrix*[3];
	Xi_delta = new Matrix*[3];
	u_i = new Matrix*[3];
	kappa_r1_i = new Matrix*[3];
	kappa_r2_i = new Matrix*[3];
	kappa_r1_delta = new Matrix*[3];
	kappa_r2_delta = new Matrix*[3];
	alpha_i = new Matrix*[3];
	z_x1_i = new Matrix*[3];
	z_x2_i = new Matrix*[3];

	eta_r1 = new Matrix*[3];			//deformação em er1
	eta_r2 = new Matrix*[3];			//deformação em er2
	kappa_r1 = new Matrix*[3];			//rot/comp em er1
	kappa_r2 = new Matrix*[3];			//rot/comp em er2
	n_r1 = new Matrix*[3];				//força em er1
	n_r2 = new Matrix*[3];				//força em er2
	m_r1 = new Matrix*[3];				//momento em er1
	m_r2 = new Matrix*[3];				//momento em er2

	eta_r1_global = new Matrix*[3];			//deformação em er1
	eta_r2_global = new Matrix*[3];			//deformação em er2
	kappa_r1_global = new Matrix*[3];			//rot/comp em er1
	kappa_r2_global = new Matrix*[3];			//rot/comp em er2
	n_r1_global = new Matrix*[3];				//força em er1
	n_r2_global = new Matrix*[3];				//força em er2
	m_r1_global = new Matrix*[3];				//momento em er1
	m_r2_global = new Matrix*[3];				//momento em er2

	//Loop nos pontos de Gauss
	for (int i = 0; i < 3; i++)
	{
		N[i] = new Matrix(6, 27);
		deltaN[i] = new Matrix(15, 27);
		alpha_delta[i] = new Matrix(3);
		alpha_delta_x1[i] = new Matrix(3);
		alpha_delta_x2[i] = new Matrix(3);
		u_delta[i] = new Matrix(3);
		u_delta_x1[i] = new Matrix(3);
		u_delta_x2[i] = new Matrix(3);

		Q_i[i] = new Matrix(3, 3);
		Xi_delta[i] = new Matrix(3, 3);
		Q_delta[i] = new Matrix(3, 3);
		u_i[i] = new Matrix(3);
		kappa_r1_i[i] = new Matrix(3);
		kappa_r2_i[i] = new Matrix(3);
		kappa_r1_delta[i] = new Matrix(3);
		kappa_r2_delta[i] = new Matrix(3);
		alpha_i[i] = new Matrix(3);
		z_x1_i[i] = new Matrix(3);
		z_x2_i[i] = new Matrix(3);
		
		eta_r1[i] = new Matrix(3);				//deformação em er1
		eta_r2[i] = new Matrix(3);				//deformação em er2
		kappa_r1[i] = new Matrix(3);			//rot/comp em er1
		kappa_r2[i] = new Matrix(3);			//rot/comp em er2
		n_r1[i] = new Matrix(3);				//força em er1
		n_r2[i] = new Matrix(3);				//força em er2
		m_r1[i] = new Matrix(3);				//momento em er1
		m_r2[i] = new Matrix(3);				//momento em er2

		eta_r1_global[i] = new Matrix(3);				//deformação em er1
		eta_r2_global[i] = new Matrix(3);				//deformação em er2
		kappa_r1_global[i] = new Matrix(3);			//rot/comp em er1
		kappa_r2_global[i] = new Matrix(3);			//rot/comp em er2
		n_r1_global[i] = new Matrix(3);				//força em er1
		n_r2_global[i] = new Matrix(3);				//força em er2
		m_r1_global[i] = new Matrix(3);				//momento em er1
		m_r2_global[i] = new Matrix(3);				//momento em er2
	}

	//Variáveis para integração com 6 pontos de Gauss (precisão de polinômios de ordem 4)
	N1u4 = new double[6];
	N2u4 = new double[6];
	N3u4 = new double[6];
	N4u4 = new double[6];
	N5u4 = new double[6];
	N6u4 = new double[6];
	N4a4 = new double[6];
	N5a4 = new double[6];
	N6a4 = new double[6];

	N1u4_x1 = new double[6];
	N2u4_x1 = new double[6];
	N3u4_x1 = new double[6];
	N4u4_x1 = new double[6];
	N5u4_x1 = new double[6];
	N6u4_x1 = new double[6];

	N1u4_x2 = new double[6];
	N2u4_x2 = new double[6];
	N3u4_x2 = new double[6];
	N4u4_x2 = new double[6];
	N5u4_x2 = new double[6];
	N6u4_x2 = new double[6];

	w4 = new double[6];
	N4 = new Matrix*[6];
	N4_u_x1 = new Matrix*[6];
	N4_u_x2 = new Matrix*[6];
	alpha_i4 = new Matrix*[6];
	//Loop nos pontos de Gauss
	for (int i = 0; i < 6; i++)
	{
		N4[i] = new Matrix(6, 27);
		N4_u_x1[i] = new Matrix(3, 18);
		N4_u_x2[i] = new Matrix(3, 18);
		alpha_i4[i] = new Matrix(3);
	}

	stiffness = new Matrix(27, 27);
	mass = new Matrix(27, 27);
	mass_modal = new Matrix(27, 27);
	damping = new Matrix(27, 27);
	damping_modal = new Matrix(27, 27);
	rayleigh_damping = new Matrix(27, 27);
	damping_loading = new Matrix(27);
	i_loading = new Matrix(27);
	inertial_loading = new Matrix(27);
	P_loading = new Matrix(27);
	e_loading = new Matrix(27);
	e1r = new Matrix(3);
	e2r = new Matrix(3);
	e3r = new Matrix(3);
	e1rlocal = new Matrix(3);
	e2rlocal = new Matrix(3);
	e3rlocal = new Matrix(3);
	transform = new Matrix(27,27);
	transform3 = new Matrix(3, 3);

	I3 = Matrix(3, 3);
	I3(0, 0) = 1.0;
	I3(1, 1) = 1.0;
	I3(2, 2) = 1.0;

	n_1 = Matrix(3, 1);
	n_2 = Matrix(3, 1);
	m_1 = Matrix(3, 1);
	m_2 = Matrix(3, 1);

	e_1 = Matrix(3, 1);
	e_2 = Matrix(3, 1);
	k_1 = Matrix(3, 1);
	k_2 = Matrix(3, 1);

	//Variáveis internas para uso na dinâmica
	alpha_dot = new Matrix*[3];
	Xi_dot = new Matrix*[3];
	Mip = new Matrix*[3];
	Jip = new Matrix*[3];
	M = new Matrix*[3];
	Md1 = new Matrix*[3];
	//Loop nos pontos de Gauss
	for (int i = 0; i < 3; i++)
	{
		alpha_dot[i] = new Matrix(3, 1);
		Xi_dot[i] = new Matrix(3, 1);
		Mip[i] = new Matrix(3, 3);
		Jip[i] = new Matrix(3, 3);
		M[i] = new Matrix(6, 6);
		Md1[i] = new Matrix(6, 6);
	}

	v_ipp = Matrix(27);//Estimativa da velocidade no instante posterior

	//Ponteiros double** - conversão de matriz
	pDdT = new double*[6];
	for (int i = 0; i < 6; i++)
		pDdT[i] = new double[6];
	DdT = new Matrix(6, 6);
	dT = new Matrix(6);

	//Composite material
	cs = 0;
	thick_comp = 0;
	rho_comp = 0;
	D_comp = new Matrix(12, 12);
}

Shell_1::~Shell_1()
{
	delete[] nodes;
	delete[] VTK_nodes;
	delete[] type_name;
	if (DOFs != NULL)
	{
		for (int i = 0; i < n_nodes; i++)
			delete[] DOFs[i];
		delete[] DOFs;
	}

	delete[]N1u;
	delete[]N2u;
	delete[]N3u;
	delete[]N4u;
	delete[]N5u;
	delete[]N6u;
	delete[]N4a;
	delete[]N5a;
	delete[]N6a;

	delete[]N1u_x1;
	delete[]N2u_x1;
	delete[]N3u_x1;
	delete[]N4u_x1;
	delete[]N5u_x1;
	delete[]N6u_x1;
	delete[]N4a_x1;
	delete[]N5a_x1;
	delete[]N6a_x1;

	delete[]N1u_x2;
	delete[]N2u_x2;
	delete[]N3u_x2;
	delete[]N4u_x2;
	delete[]N5u_x2;
	delete[]N6u_x2;
	delete[]N4a_x2;
	delete[]N5a_x2;
	delete[]N6a_x2;

	////////////////////////
	delete[]N1u4;
	delete[]N2u4;
	delete[]N3u4;
	delete[]N4u4;
	delete[]N5u4;
	delete[]N6u4;
	delete[]N4a4;
	delete[]N5a4;
	delete[]N6a4;

	delete[]N1u4_x1;
	delete[]N2u4_x1;
	delete[]N3u4_x1;
	delete[]N4u4_x1;
	delete[]N5u4_x1;
	delete[]N6u4_x1;

	delete[]N1u4_x2;
	delete[]N2u4_x2;
	delete[]N3u4_x2;
	delete[]N4u4_x2;
	delete[]N5u4_x2;
	delete[]N6u4_x2;

	delete[]w4;
	for (int i = 0; i < 6; i++)
	{
		delete N4[i];
		delete N4_u_x1[i];
		delete N4_u_x2[i];
		delete alpha_i4[i];
	}
	delete[] N4;
	delete[] N4_u_x1;
	delete[] N4_u_x2;
	delete[] alpha_i4;
	////////////////////////////

	for (int i = 0; i < 3; i++)
	{
		delete N[i];
		delete deltaN[i];
		delete u_delta[i];
		delete u_delta_x1[i];
		delete u_delta_x2[i];
		delete alpha_delta[i];
		delete alpha_delta_x1[i];
		delete alpha_delta_x2[i];
		delete Q_i[i];
		delete Q_delta[i];
		delete Xi_delta[i];
		delete u_i[i];
		delete kappa_r1_i[i];
		delete kappa_r2_i[i];
		delete kappa_r1_delta[i];
		delete kappa_r2_delta[i];
		delete alpha_i[i];
		delete z_x1_i[i];
		delete z_x2_i[i];
		
		delete eta_r1[i];			//deformação em er1
		delete eta_r2[i];			//deformação em er2
		delete kappa_r1[i];			//rot/comp em er1
		delete kappa_r2[i];			//rot/comp em er2
		delete n_r1[i];				//força em er1
		delete n_r2[i];				//força em er2
		delete m_r1[i];				//momento em er1
		delete m_r2[i];				//momento em er2

		delete eta_r1_global[i];			//deformação em er1
		delete eta_r2_global[i];			//deformação em er2
		delete kappa_r1_global[i];			//rot/comp em er1
		delete kappa_r2_global[i];			//rot/comp em er2
		delete n_r1_global[i];				//força em er1
		delete n_r2_global[i];				//força em er2
		delete m_r1_global[i];				//momento em er1
		delete m_r2_global[i];				//momento em er2
	}
	delete[] N;
	delete[] deltaN;
	delete[] u_delta;
	delete[] u_delta_x1;
	delete[] u_delta_x2;
	delete[] alpha_delta;
	delete[] alpha_delta_x1;
	delete[] alpha_delta_x2;
	delete[] Q_i;
	delete[] Q_delta;
	delete[] Xi_delta;
	delete[] u_i;
	delete[] kappa_r1_i;
	delete[] kappa_r2_i;
	delete[] kappa_r1_delta;
	delete[] kappa_r2_delta;
	delete[] alpha_i;
	delete[] z_x1_i;
	delete[] z_x2_i;
	
	delete[] eta_r1;			//deformação em er1
	delete[] eta_r2;			//deformação em er2
	delete[] kappa_r1;			//rot/comp em er1
	delete[] kappa_r2;			//rot/comp em er2
	delete[] n_r1;				//força em er1
	delete[] n_r2;				//força em er2
	delete[] m_r1;				//momento em er1
	delete[] m_r2;				//momento em er2

	delete[] eta_r1_global;				//deformação em er1
	delete[] eta_r2_global;				//deformação em er2
	delete[] kappa_r1_global;			//rot/comp em er1
	delete[] kappa_r2_global;			//rot/comp em er2
	delete[] n_r1_global;				//força em er1
	delete[] n_r2_global;				//força em er2
	delete[] m_r1_global;				//momento em er1
	delete[] m_r2_global;				//momento em er2

	delete stiffness;
	delete mass;
	delete mass_modal;
	delete damping;
	delete damping_modal;
	delete rayleigh_damping;
	delete damping_loading;
	delete i_loading;
	delete inertial_loading;
	delete P_loading;
	delete e_loading;
	delete e1r;
	delete e2r;
	delete e3r;
	delete e1rlocal;
	delete e2rlocal;
	delete e3rlocal;
	delete transform;
	delete transform3;

	delete tempkin;

	//Variáveis internas para uso na dinâmica
	//Loop nos pontos de Gauss
	for (int i = 0; i < 3; i++)
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
	delete[] Mip;
	delete[] Jip;
	delete[] M;
	delete[] Md1;

	delete DdT;
	delete dT;
	//Ponteiros de double**
	for (int i = 0; i < 6; i++)
		delete[]pDdT[i];
	delete[] pDdT;

	delete D_comp;
}

//Checa inconsistências no elemento para evitar erros de execução
bool Shell_1::Check()
{
	for (int i = 0; i < n_nodes; i++)
	{
		if (nodes[i] > db.number_nodes)
		{
			//printf("NODES %d %d",data.number_nodes,nodes[i]);
			return false;
		}
			
	}
	if (cs > db.number_CS)
	{
		printf("CS %d %d", db.number_CS, cs);
		return false;
	}
		
	if (section > db.number_shell_sections)
	{
		//printf("SEC %d %d", data.number_shell_sections, section);
		return false;
	}
		
	if (material > db.number_materials)
	{
		//printf("MAT %d %d", data.number_materials, material);
		return false;
	}
		
	return true;
}

bool Shell_1::Read(FILE *f)
{
	char s[1000];
	fscanf(f, "%s", s);
	number = atoi(s);

	fscanf(f, "%s", s);
	if (!strcmp(s, "Mat"))
	{
		fscanf(f, "%s", s);
		material = atoi(s);
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "Sec"))
	{
		fscanf(f, "%s", s);
		section = atoi(s);
	}
	else
		return false;

	//Composite
	fpos_t pos;
	fgetpos(f, &pos);
	fscanf(f, "%s", s);
	if (!strcmp(s, "CS"))
	{
		fscanf(f, "%s", s);
		cs = atoi(s);
	}
	else
		fsetpos(f, &pos);

	fscanf(f, "%s", s);
	if (!strcmp(s, "Nodes"))
	{
		fscanf(f, "%s", s);
		nodes[0] = atoi(s);
		fscanf(f, "%s", s);
		nodes[1] = atoi(s);
		fscanf(f, "%s", s);
		nodes[2] = atoi(s);
		fscanf(f, "%s", s);
		nodes[3] = atoi(s);
		fscanf(f, "%s", s);
		nodes[4] = atoi(s);
		fscanf(f, "%s", s);
		nodes[5] = atoi(s);
	}
	else
		return false;
	return true;
}

void Shell_1::Write(FILE *f)
{
	if (cs != 0) {
		fprintf(f, "Shell_1\t%d\tMat\t%d\tSec\t%d\tCS\t%d\tNodes\t%d\t%d\t%d\t%d\t%d\t%d\n",
			number,
			material,
			section,
			cs,
			nodes[0],
			nodes[1],
			nodes[2],
			nodes[3],
			nodes[4],
			nodes[5]);
	}
	else {
		fprintf(f, "Shell_1\t%d\tMat\t%d\tSec\t%d\tNodes\t%d\t%d\t%d\t%d\t%d\t%d\n",
			number,
			material,
			section,
			nodes[0],
			nodes[1],
			nodes[2],
			nodes[3],
			nodes[4],
			nodes[5]);
	}
}

//Escreve arquivo de resultados
void Shell_1::WriteResults(FILE *f)
{
	//Zerando os vetores de resultados
	zeros(&n_1);
	zeros(&n_2);
	zeros(&m_1);
	zeros(&m_2);
	zeros(&e_1);
	zeros(&e_2);
	zeros(&k_1);
	zeros(&k_2);
	//Percorre os pontos de Gauss e calcula um único vetor de resultados, por elemento (média)
	for (int gauss = 0; gauss < 3; gauss++)
	{
		n_1 = n_1 + (1.0 / 3.0)*(*n_r1[gauss]);
		n_2 = n_2 + (1.0 / 3.0)*(*n_r2[gauss]);
		m_1 = m_1 + (1.0 / 3.0)*(*m_r1[gauss]);
		m_2 = m_2 + (1.0 / 3.0)*(*m_r2[gauss]);

		e_1 = e_1 + (1.0 / 3.0)*(*eta_r1[gauss]);
		e_2 = e_2 + (1.0 / 3.0)*(*eta_r2[gauss]);
		k_1 = k_1 + (1.0 / 3.0)*(*kappa_r1[gauss]);
		k_2 = k_2 + (1.0 / 3.0)*(*kappa_r2[gauss]);
	}

	fprintf(f, "Shell_1\t%d\tn1\t%.6e\t%.6e\t%.6e\tm1\t%.6e\t%.6e\t%.6e\tn2\t%.6e\t%.6e\t%.6e\tm2\t%.6e\t%.6e\t%.6e\teta1\t%.6e\t%.6e\t%.6e\tkappa1\t%.6e\t%.6e\t%.6e\teta2\t%.6e\t%.6e\t%.6e\tkappa2\t%.6e\t%.6e\t%.6e\tRef_Coor\t%.6e\t%.6e\t%.6e\tCoor\t%.6e\t%.6e\t%.6e\n",
		number,
		n_1(0, 0), n_1(1, 0),n_1(2, 0),
		m_1(0, 0), m_1(1, 0), m_1(2, 0),
		n_2(0, 0), n_2(1, 0), n_2(2, 0),
		m_2(0, 0), m_2(1, 0), m_2(2, 0),
		e_1(0, 0), e_1(1, 0), e_1(2, 0),
		k_1(0, 0), k_1(1, 0), k_1(2, 0),
		e_2(0, 0), e_2(1, 0), e_2(2, 0),
		k_2(0, 0), k_2(1, 0), k_2(2, 0),
		(db.nodes[nodes[0] - 1]->ref_coordinates[0] + db.nodes[nodes[1] - 1]->ref_coordinates[0] + db.nodes[nodes[2] - 1]->ref_coordinates[0]) / 3,
		(db.nodes[nodes[0] - 1]->ref_coordinates[1] + db.nodes[nodes[1] - 1]->ref_coordinates[1] + db.nodes[nodes[2] - 1]->ref_coordinates[1]) / 3,
		(db.nodes[nodes[0] - 1]->ref_coordinates[2] + db.nodes[nodes[1] - 1]->ref_coordinates[2] + db.nodes[nodes[2] - 1]->ref_coordinates[2]) / 3,
		(db.nodes[nodes[0] - 1]->copy_coordinates[0] + db.nodes[nodes[1] - 1]->copy_coordinates[0] + db.nodes[nodes[2] - 1]->copy_coordinates[0]) / 3 + (db.nodes[nodes[0] - 1]->displacements[0] + db.nodes[nodes[1] - 1]->displacements[0] + db.nodes[nodes[2] - 1]->displacements[0]) / 3,
		(db.nodes[nodes[0] - 1]->copy_coordinates[1] + db.nodes[nodes[1] - 1]->copy_coordinates[1] + db.nodes[nodes[2] - 1]->copy_coordinates[1]) / 3 + (db.nodes[nodes[0] - 1]->displacements[1] + db.nodes[nodes[1] - 1]->displacements[1] + db.nodes[nodes[2] - 1]->displacements[1]) / 3,
		(db.nodes[nodes[0] - 1]->copy_coordinates[2] + db.nodes[nodes[1] - 1]->copy_coordinates[2] + db.nodes[nodes[2] - 1]->copy_coordinates[2]) / 3 + (db.nodes[nodes[0] - 1]->displacements[2] + db.nodes[nodes[1] - 1]->displacements[2] + db.nodes[nodes[2] - 1]->displacements[2]) / 3
		);
}

void Shell_1::WriteVTK_XMLBase(std::vector<float> *float_vector)
{
	//Zerando os vetores de resultados
	zeros(&n_1);
	zeros(&n_2);
	zeros(&m_1);
	zeros(&m_2);
	//Percorre os pontos de Gauss e calcula um único vetor de resultados, por elemento (média)
	for (int gauss = 0; gauss < 3; gauss++)
	{
		n_1 = n_1 + (1.0 / 3.0)*(*n_r1[gauss]);
		n_2 = n_2 + (1.0 / 3.0)*(*n_r2[gauss]);
		m_1 = m_1 + (1.0 / 3.0)*(*m_r1[gauss]);
		m_2 = m_2 + (1.0 / 3.0)*(*m_r2[gauss]);
	}
	//Imprime os resultados do elemento
	int res_element = 12;
	float_vector->push_back((float)(n_1(0, 0)));
	float_vector->push_back((float)(n_1(1, 0)));
	float_vector->push_back((float)(n_1(2, 0)));

	float_vector->push_back((float)(m_1(0, 0)));
	float_vector->push_back((float)(m_1(1, 0)));
	float_vector->push_back((float)(m_1(2, 0)));

	float_vector->push_back((float)(n_2(0, 0)));
	float_vector->push_back((float)(n_2(1, 0)));
	float_vector->push_back((float)(n_2(2, 0)));

	float_vector->push_back((float)(m_2(0, 0)));
	float_vector->push_back((float)(m_2(1, 0)));
	float_vector->push_back((float)(m_2(2, 0)));
	/*fprintf(f, "\t\t\t\t\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t", 
		n_1(0, 0), n_1(1, 0), n_1(2, 0), m_1(0, 0), m_1(1, 0), m_1(2, 0), n_2(0, 0), n_2(1, 0), n_2(2, 0), m_2(0, 0), m_2(1, 0), m_2(2, 0));*/
	//Imprime valores nulos para resultados que não fazem sentido para esse tipo de elemento
	for (int i = res_element; i < db.n_element_results; i++)
		float_vector->push_back(0.0);
		//fprintf(f, "%f\t", 0.0);
	//fprintf(f, "\n");
}
void Shell_1::WriteVTK_XMLRender(FILE *f)
{
	//vetores para escrita no formato binário - usando a função 'enconde'
	std::vector<float> float_vector;
	std::vector<int> int_vector;

	double p1x, p1y, p1z, p2x, p2y, p2z, p3x, p3y, p3z, p0x, p0y, p0z, p4x, p4y, p4z, p5x, p5y, p5z, p6x, p6y, p6z, p7x, p7y, p7z, p8x, p8y, p8z, p9x, p9y, p9z, p10x, p10y, p10z, a, m, n;
	//Atribuições
	p0x = db.nodes[nodes[0] - 1]->ref_coordinates[0] + db.post_files->mag_factor*(db.nodes[nodes[0] - 1]->copy_coordinates[0] - db.nodes[nodes[0] - 1]->ref_coordinates[0]);
	p0y = db.nodes[nodes[0] - 1]->ref_coordinates[1] + db.post_files->mag_factor*(db.nodes[nodes[0] - 1]->copy_coordinates[1] - db.nodes[nodes[0] - 1]->ref_coordinates[1]);
	p0z = db.nodes[nodes[0] - 1]->ref_coordinates[2] + db.post_files->mag_factor*(db.nodes[nodes[0] - 1]->copy_coordinates[2] - db.nodes[nodes[0] - 1]->ref_coordinates[2]);
	p1x = db.nodes[nodes[1] - 1]->ref_coordinates[0] + db.post_files->mag_factor*(db.nodes[nodes[1] - 1]->copy_coordinates[0] - db.nodes[nodes[1] - 1]->ref_coordinates[0]);
	p1y = db.nodes[nodes[1] - 1]->ref_coordinates[1] + db.post_files->mag_factor*(db.nodes[nodes[1] - 1]->copy_coordinates[1] - db.nodes[nodes[1] - 1]->ref_coordinates[1]);
	p1z = db.nodes[nodes[1] - 1]->ref_coordinates[2] + db.post_files->mag_factor*(db.nodes[nodes[1] - 1]->copy_coordinates[2] - db.nodes[nodes[1] - 1]->ref_coordinates[2]);
	p2x = db.nodes[nodes[2] - 1]->ref_coordinates[0] + db.post_files->mag_factor*(db.nodes[nodes[2] - 1]->copy_coordinates[0] - db.nodes[nodes[2] - 1]->ref_coordinates[0]);
	p2y = db.nodes[nodes[2] - 1]->ref_coordinates[1] + db.post_files->mag_factor*(db.nodes[nodes[2] - 1]->copy_coordinates[1] - db.nodes[nodes[2] - 1]->ref_coordinates[1]);
	p2z = db.nodes[nodes[2] - 1]->ref_coordinates[2] + db.post_files->mag_factor*(db.nodes[nodes[2] - 1]->copy_coordinates[2] - db.nodes[nodes[2] - 1]->ref_coordinates[2]);
	a = db.shell_sections[section - 1]->thickness;
	if (composite_shell == true) {
		a = thick_comp;
	}
	float_vector.clear();
	fprintf(f, "<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n<Points>\n<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"binary\">\n", 6, 1);
	a = a / 2;
	p3x = ((p1y - p0y)*(p2z - p0z)) - ((p2y - p0y)*(p1z - p0z));  //produto vetorial para achar o vetor normal p3
	p3y = ((p1z - p0z)*(p2x - p0x)) - ((p1x - p0x)*(p2z - p0z));
	p3z = ((p1x - p0x)*(p2y - p0y)) - ((p1y - p0y)*(p2x - p0x));

	m = pow(p3x, 2) + pow(p3y, 2) + pow(p3z, 2);
	n = pow(m, 0.5);  //módulo do vetor p3

	p4x = p3x / n;
	p4y = p3y / n;
	p4z = p3z / n;  //normalizo p3

	p5x = p0x + a*p4x;  //encontro os novos pontos, acima ou abaixo dos pontos dados, paralelos a p3
	p5y = p0y + a*p4y;
	p5z = p0z + a*p4z;

	p6x = p1x + a*p4x;
	p6y = p1y + a*p4y;
	p6z = p1z + a*p4z;

	p7x = p2x + a*p4x;
	p7y = p2y + a*p4y;
	p7z = p2z + a*p4z;

	p8x = p0x - a*p4x;
	p8y = p0y - a*p4y;
	p8z = p0z - a*p4z;

	p9x = p1x - a*p4x;
	p9y = p1y - a*p4y;
	p9z = p1z - a*p4z;

	p10x = p2x - a*p4x;
	p10y = p2y - a*p4y;
	p10z = p2z - a*p4z;
	float_vector.push_back((float)(p5x));
	float_vector.push_back((float)(p5y));
	float_vector.push_back((float)(p5z));
	float_vector.push_back((float)(p6x));
	float_vector.push_back((float)(p6y));
	float_vector.push_back((float)(p6z));
	float_vector.push_back((float)(p7x));
	float_vector.push_back((float)(p7y));
	float_vector.push_back((float)(p7z));
	float_vector.push_back((float)(p8x));
	float_vector.push_back((float)(p8y));
	float_vector.push_back((float)(p8z));
	float_vector.push_back((float)(p9x));
	float_vector.push_back((float)(p9y));
	float_vector.push_back((float)(p9z));
	float_vector.push_back((float)(p10x));
	float_vector.push_back((float)(p10y));
	float_vector.push_back((float)(p10z));
	//fprintf(f, "%f %f %f\n%f %f %f\n%f %f %f\n%f %f %f\n%f %f %f\n%f %f %f\n", p5x, p5y, p5z, p6x, p6y, p6z, p7x, p7y, p7z, p8x, p8y, p8z, p9x, p9y, p9z, p10x, p10y, p10z);
	fprintf(f, encodeData(float_vector).c_str());
	fprintf(f, "\n");
	int_vector.clear();
	fprintf(f, "</DataArray>\n</Points>\n<Cells>\n<DataArray type=\"Int32\" Name=\"connectivity\" format=\"binary\">\n");
	int_vector.push_back(0);
	int_vector.push_back(1);
	int_vector.push_back(2);
	int_vector.push_back(3);
	int_vector.push_back(4);
	int_vector.push_back(5);
	fprintf(f, encodeData(int_vector).c_str());
	fprintf(f, "\n");
	int_vector.clear();
	fprintf(f, "</DataArray>\n<DataArray type=\"Int32\" Name=\"offsets\" format=\"binary\">\n");
	int_vector.push_back(6);
	fprintf(f, encodeData(int_vector).c_str());
	fprintf(f, "\n");
	int_vector.clear();
	fprintf(f, "</DataArray>\n<DataArray type=\"Int32\" Name=\"types\" format=\"binary\">\n");
	int_vector.push_back(13);
	fprintf(f, encodeData(int_vector).c_str());
	fprintf(f, "\n");
	fprintf(f, "</DataArray>\n</Cells>\n");
	//Opens CellData
	fprintf(f, "<CellData FieldData = \"ElementData\">\n");
	//Opens DataArray
	fprintf(f, "<DataArray Name = \"ElementResults\" type = \"Float32\" NumberOfComponents=\"%d\" format=\"binary\">\n", db.n_element_results);
	//Imprime os resultados do elemento
	int res_element = 12;
	float_vector.clear();
	float_vector.push_back((float)(n_1(0, 0)));
	float_vector.push_back((float)(n_1(1, 0)));
	float_vector.push_back((float)(n_1(2, 0)));

	float_vector.push_back((float)(m_1(0, 0)));
	float_vector.push_back((float)(m_1(1, 0)));
	float_vector.push_back((float)(m_1(2, 0)));

	float_vector.push_back((float)(n_2(0, 0)));
	float_vector.push_back((float)(n_2(1, 0)));
	float_vector.push_back((float)(n_2(2, 0)));

	float_vector.push_back((float)(m_2(0, 0)));
	float_vector.push_back((float)(m_2(1, 0)));
	float_vector.push_back((float)(m_2(2, 0)));

	/*fprintf(f, "\t\t\t\t\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t",
		n_1(0, 0), n_1(1, 0), n_1(2, 0), m_1(0, 0), m_1(1, 0), m_1(2, 0), n_2(0, 0), n_2(1, 0), n_2(2, 0), m_2(0, 0), m_2(1, 0), m_2(2, 0));*/
	//Imprime valores nulos para resultados que não fazem sentido para esse tipo de elemento
	for (int i = res_element; i < db.n_element_results; i++)
		float_vector.push_back((float)(0.0));
	fprintf(f, encodeData(float_vector).c_str());
	fprintf(f, "\n");
	//Closes DataArray
	fprintf(f, "</DataArray>\n");

	int_vector.clear();
	//Opens DataArray
	fprintf(f, "<DataArray Name=\"ElementProperties\" type=\"Int32\" NumberOfComponents=\"%d\" format=\"binary\">\n", 4);
	int_vector.push_back(3);		//Element ID
	int_vector.push_back(material);
	int_vector.push_back(section);
	int_vector.push_back(0);
	fprintf(f, encodeData(int_vector).c_str());
	fprintf(f, "\n");
	//Closes DataArray
	fprintf(f, "</DataArray>\n");


	//Closes CellData
	fprintf(f, "</CellData>\n");
	fprintf(f, "</Piece>\n");
}

//Escreve no monitor do elemento//Escreve no monitor do elemento
void Shell_1::WriteMonitor(FILE *f, bool first_record, double time)
{
	//Zerando os vetores de resultados
	zeros(&n_1);
	zeros(&n_2);
	zeros(&m_1);
	zeros(&m_2);
	zeros(&e_1);
	zeros(&e_2);
	zeros(&k_1);
	zeros(&k_2);
	//Percorre os pontos de Gauss e calcula um único vetor de resultados, por elemento (média)
	for (int gauss = 0; gauss < 3; gauss++)
	{
		n_1 = n_1 + 1.0 / 3.0*(*n_r1[gauss]);
		n_2 = n_2 + 1.0 / 3.0*(*n_r2[gauss]);
		m_1 = m_1 + 1.0 / 3.0*(*m_r1[gauss]);
		m_2 = m_2 + 1.0 / 3.0*(*m_r2[gauss]);

		/*e_1 = e_1 + 1.0 / 3.0*(*eta_r1[gauss]);
		e_2 = e_2 + 1.0 / 3.0*(*eta_r2[gauss]);
		k_1 = k_1 + 1.0 / 3.0*(*kappa_r1[gauss]);
		k_2 = k_2 + 1.0 / 3.0*(*kappa_r2[gauss]);*/
	}

	//Cabeçalho
	if (first_record == true)
		fprintf(f, "TIME\tn1x\tn1y\tn1z\tn2x\tn2y\tn2z\tm1x\tm1y\tm1z\tm2x\tm2y\tm2z\n");
	//Informações a serem salvas
	fprintf(f, "%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\n",
		time,
		n_1(0, 0), n_1(1, 0), n_1(2, 0),
		n_2(0, 0), n_2(1, 0), n_2(2, 0),
		m_1(0, 0), m_1(1, 0), m_1(2, 0),
		m_2(0, 0), m_2(1, 0), m_2(2, 0)
		);
}

//Monta elemento
void Shell_1::Mount()
{
	Zeros();
	//Loop nos pontos de Gauss
	for (int gauss = 0; gauss < 3; gauss++)
	{
		//u_delta calculado no ponto de Gauss							
		(*u_delta[gauss])(0, 0)		=		(db.nodes[nodes[0] - 1]->displacements[0])* N1u[gauss] +
											(db.nodes[nodes[1] - 1]->displacements[0])* N2u[gauss] +
											(db.nodes[nodes[2] - 1]->displacements[0])* N3u[gauss] + 
											(db.nodes[nodes[3] - 1]->displacements[0])* N4u[gauss] + 
											(db.nodes[nodes[4] - 1]->displacements[0])* N5u[gauss] + 
											(db.nodes[nodes[5] - 1]->displacements[0])* N6u[gauss] ;
		(*u_delta[gauss])(1, 0)		=		(db.nodes[nodes[0] - 1]->displacements[1])* N1u[gauss] +
											(db.nodes[nodes[1] - 1]->displacements[1])* N2u[gauss] +
											(db.nodes[nodes[2] - 1]->displacements[1])* N3u[gauss] +
											(db.nodes[nodes[3] - 1]->displacements[1])* N4u[gauss] +
											(db.nodes[nodes[4] - 1]->displacements[1])* N5u[gauss] +
											(db.nodes[nodes[5] - 1]->displacements[1])* N6u[gauss] ;
		(*u_delta[gauss])(2, 0)		=		(db.nodes[nodes[0] - 1]->displacements[2])* N1u[gauss] +
											(db.nodes[nodes[1] - 1]->displacements[2])* N2u[gauss] +
											(db.nodes[nodes[2] - 1]->displacements[2])* N3u[gauss] +
											(db.nodes[nodes[3] - 1]->displacements[2])* N4u[gauss] +
											(db.nodes[nodes[4] - 1]->displacements[2])* N5u[gauss] +
											(db.nodes[nodes[5] - 1]->displacements[2])* N6u[gauss] ;
		//u_delta_x1 calculado no ponto de Gauss							
		(*u_delta_x1[gauss])(0, 0)	=		(db.nodes[nodes[0] - 1]->displacements[0])* N1u_x1[gauss] +
											(db.nodes[nodes[1] - 1]->displacements[0])* N2u_x1[gauss] +
											(db.nodes[nodes[2] - 1]->displacements[0])* N3u_x1[gauss] +
											(db.nodes[nodes[3] - 1]->displacements[0])* N4u_x1[gauss] +
											(db.nodes[nodes[4] - 1]->displacements[0])* N5u_x1[gauss] +
											(db.nodes[nodes[5] - 1]->displacements[0])* N6u_x1[gauss] ;
		(*u_delta_x1[gauss])(1, 0)	=		(db.nodes[nodes[0] - 1]->displacements[1])* N1u_x1[gauss] +
											(db.nodes[nodes[1] - 1]->displacements[1])* N2u_x1[gauss] +
											(db.nodes[nodes[2] - 1]->displacements[1])* N3u_x1[gauss] +
											(db.nodes[nodes[3] - 1]->displacements[1])* N4u_x1[gauss] +
											(db.nodes[nodes[4] - 1]->displacements[1])* N5u_x1[gauss] +
											(db.nodes[nodes[5] - 1]->displacements[1])* N6u_x1[gauss] ;
		(*u_delta_x1[gauss])(2, 0)	=		(db.nodes[nodes[0] - 1]->displacements[2])* N1u_x1[gauss] +
											(db.nodes[nodes[1] - 1]->displacements[2])* N2u_x1[gauss] +
											(db.nodes[nodes[2] - 1]->displacements[2])* N3u_x1[gauss] +
											(db.nodes[nodes[3] - 1]->displacements[2])* N4u_x1[gauss] +
											(db.nodes[nodes[4] - 1]->displacements[2])* N5u_x1[gauss] +
											(db.nodes[nodes[5] - 1]->displacements[2])* N6u_x1[gauss] ;
		//u_delta_x2 calculado no ponto de Gauss							
		(*u_delta_x2[gauss])(0, 0)	=		(db.nodes[nodes[0] - 1]->displacements[0])* N1u_x2[gauss] +
											(db.nodes[nodes[1] - 1]->displacements[0])* N2u_x2[gauss] +
											(db.nodes[nodes[2] - 1]->displacements[0])* N3u_x2[gauss] +
											(db.nodes[nodes[3] - 1]->displacements[0])* N4u_x2[gauss] +
											(db.nodes[nodes[4] - 1]->displacements[0])* N5u_x2[gauss] +
											(db.nodes[nodes[5] - 1]->displacements[0])* N6u_x2[gauss] ;
		(*u_delta_x2[gauss])(1, 0)	=		(db.nodes[nodes[0] - 1]->displacements[1])* N1u_x2[gauss] +
											(db.nodes[nodes[1] - 1]->displacements[1])* N2u_x2[gauss] +
											(db.nodes[nodes[2] - 1]->displacements[1])* N3u_x2[gauss] +
											(db.nodes[nodes[3] - 1]->displacements[1])* N4u_x2[gauss] +
											(db.nodes[nodes[4] - 1]->displacements[1])* N5u_x2[gauss] +
											(db.nodes[nodes[5] - 1]->displacements[1])* N6u_x2[gauss] ;
		(*u_delta_x2[gauss])(2, 0)	=		(db.nodes[nodes[0] - 1]->displacements[2])* N1u_x2[gauss] +
											(db.nodes[nodes[1] - 1]->displacements[2])* N2u_x2[gauss] +
											(db.nodes[nodes[2] - 1]->displacements[2])* N3u_x2[gauss] +
											(db.nodes[nodes[3] - 1]->displacements[2])* N4u_x2[gauss] +
											(db.nodes[nodes[4] - 1]->displacements[2])* N5u_x2[gauss] +
											(db.nodes[nodes[5] - 1]->displacements[2])* N6u_x2[gauss] ;
		//alpha_delta calculado no ponto de Gauss							
		(*alpha_delta[gauss])(0, 0) =		(db.nodes[nodes[3] - 1]->displacements[3])* N4a[gauss] +
											(db.nodes[nodes[4] - 1]->displacements[3])* N5a[gauss] +
											(db.nodes[nodes[5] - 1]->displacements[3])* N6a[gauss] ;
		(*alpha_delta[gauss])(1, 0) =		(db.nodes[nodes[3] - 1]->displacements[4])* N4a[gauss] +
											(db.nodes[nodes[4] - 1]->displacements[4])* N5a[gauss] +
											(db.nodes[nodes[5] - 1]->displacements[4])* N6a[gauss] ;
		(*alpha_delta[gauss])(2, 0) =		(db.nodes[nodes[3] - 1]->displacements[5])* N4a[gauss] +
											(db.nodes[nodes[4] - 1]->displacements[5])* N5a[gauss] +
											(db.nodes[nodes[5] - 1]->displacements[5])* N6a[gauss] ;
		//alpha_delta_x1 calculado no ponto de Gauss							
		(*alpha_delta_x1[gauss])(0, 0) =	(db.nodes[nodes[3] - 1]->displacements[3])* N4a_x1[gauss] +
											(db.nodes[nodes[4] - 1]->displacements[3])* N5a_x1[gauss] +
											(db.nodes[nodes[5] - 1]->displacements[3])* N6a_x1[gauss] ;
		(*alpha_delta_x1[gauss])(1, 0) =	(db.nodes[nodes[3] - 1]->displacements[4])* N4a_x1[gauss] +
											(db.nodes[nodes[4] - 1]->displacements[4])* N5a_x1[gauss] +
											(db.nodes[nodes[5] - 1]->displacements[4])* N6a_x1[gauss] ;
		(*alpha_delta_x1[gauss])(2, 0) =	(db.nodes[nodes[3] - 1]->displacements[5])* N4a_x1[gauss] +
											(db.nodes[nodes[4] - 1]->displacements[5])* N5a_x1[gauss] +
											(db.nodes[nodes[5] - 1]->displacements[5])* N6a_x1[gauss] ;
		//alpha_delta_x2 calculado no ponto de Gauss							
		(*alpha_delta_x2[gauss])(0, 0) =	(db.nodes[nodes[3] - 1]->displacements[3])* N4a_x2[gauss] +
											(db.nodes[nodes[4] - 1]->displacements[3])* N5a_x2[gauss] +
											(db.nodes[nodes[5] - 1]->displacements[3])* N6a_x2[gauss] ;
		(*alpha_delta_x2[gauss])(1, 0) =	(db.nodes[nodes[3] - 1]->displacements[4])* N4a_x2[gauss] +
											(db.nodes[nodes[4] - 1]->displacements[4])* N5a_x2[gauss] +
											(db.nodes[nodes[5] - 1]->displacements[4])* N6a_x2[gauss] ;
		(*alpha_delta_x2[gauss])(2, 0) =	(db.nodes[nodes[3] - 1]->displacements[5])* N4a_x2[gauss] +
											(db.nodes[nodes[4] - 1]->displacements[5])* N5a_x2[gauss] +
											(db.nodes[nodes[5] - 1]->displacements[5])* N6a_x2[gauss] ;
		
		*u_delta[gauss] = (*transform3)*(*u_delta[gauss]);					//TO LOCAL
		*u_delta_x1[gauss] = (*transform3)*(*u_delta_x1[gauss]);			//TO LOCAL
		*u_delta_x2[gauss] = (*transform3)*(*u_delta_x2[gauss]);			//TO LOCAL
		*alpha_delta[gauss] = (*transform3)*(*alpha_delta[gauss]);			//TO LOCAL
		*alpha_delta_x1[gauss] = (*transform3)*(*alpha_delta_x1[gauss]);	//TO LOCAL
		*alpha_delta_x2[gauss] = (*transform3)*(*alpha_delta_x2[gauss]);	//TO LOCAL
		
		//Cálculo de tensores de rotação e outros
		double alpha_escalar_delta = norm(*alpha_delta[gauss]);																//Valor escalar do parametro alpha
		Matrix A_delta = skew(*alpha_delta[gauss]);																			//Matriz A_delta
		double g = 4.0 / (4.0 + alpha_escalar_delta*alpha_escalar_delta);													//função g(alpha) - em algumas ref. tb. chamado de h(alpha)
		*Q_delta[gauss] = I3 + g*(A_delta + 0.5*((A_delta)*(A_delta)));														//Tensor Q_delta
		*Xi_delta[gauss] = g*(I3 + 0.5*A_delta);																			//Tensor Xi
		Matrix A_delta_x1 = skew(*alpha_delta_x1[gauss]);																	//Derivada do vetor d_alpha em x1
		Matrix A_delta_x2 = skew(*alpha_delta_x2[gauss]);																	//Derivada do vetor d_alpha em x2
		Matrix Xi_delta_x1 = -0.5*g*(dot(*alpha_delta[gauss], *alpha_delta_x1[gauss])*(*Xi_delta[gauss]) - A_delta_x1);		//Derivada do tensor Xi em x1
		Matrix Xi_delta_x2 = -0.5*g*(dot(*alpha_delta[gauss], *alpha_delta_x2[gauss])*(*Xi_delta[gauss]) - A_delta_x2);		//Derivada do tensor Xi em x2
		Matrix z_x1 = *u_delta_x1[gauss] + (*z_x1_i[gauss]);																//Vetor derivada de z em x1 (calculada em i+1)
		Matrix z_x2 = *u_delta_x2[gauss] + (*z_x2_i[gauss]);																//Vetor derivada de z em x2 (calculada em i+1)
		Matrix Z_x1 = skew(z_x1);																							//Tensor anti-simétrico, cuja axial é z_x1
		Matrix Z_x2 = skew(z_x2);																							//Tensor anti-simétrico, cuja axial é z_x2
		Matrix Qtransp = transp((*Q_delta[gauss])*(*Q_i[gauss]));															//Q transposta (i+1)

		//Deformações retro-rotacionadas
		*eta_r1[gauss] = Qtransp*z_x1 - *e1rlocal;
		*eta_r2[gauss] = Qtransp*z_x2 - *e2rlocal;
		*kappa_r1[gauss] = transp(*Q_i[gauss])*transp(*Xi_delta[gauss])*(*alpha_delta_x1[gauss]) + *kappa_r1_i[gauss];
		*kappa_r2[gauss] = transp(*Q_i[gauss])*transp(*Xi_delta[gauss])*(*alpha_delta_x2[gauss]) + *kappa_r2_i[gauss];
		
		//Variáveis a serem calculadas na integração na espessura:
		Matrix gamma_r1(3);			//Deformações em er1
		Matrix gamma_r2(3);			//Deformações em er2
		zeros(n_r1[gauss]);				//força em er1
		zeros(n_r2[gauss]);				//força em er2
		zeros(m_r1[gauss]);				//momento em er1
		zeros(m_r2[gauss]);				//momento em er2
		Matrix C_11r(3, 3);
		Matrix C_12r(3, 3);
		Matrix C_21r(3, 3);
		Matrix C_22r(3, 3);

		Matrix d_n1_d_eta1(3, 3);
		Matrix d_n1_d_eta2(3, 3);
		Matrix d_n1_d_kappa1(3, 3);
		Matrix d_n1_d_kappa2(3, 3);

		Matrix d_n2_d_eta1(3, 3);
		Matrix d_n2_d_eta2(3, 3);
		Matrix d_n2_d_kappa1(3, 3);
		Matrix d_n2_d_kappa2(3, 3);

		Matrix d_m1_d_eta1(3, 3);
		Matrix d_m1_d_eta2(3, 3);
		Matrix d_m1_d_kappa1(3, 3);
		Matrix d_m1_d_kappa2(3, 3);

		Matrix d_m2_d_eta1(3, 3);
		Matrix d_m2_d_eta2(3, 3);
		Matrix d_m2_d_kappa1(3, 3);
		Matrix d_m2_d_kappa2(3, 3);
		double g11, g12, g13, g21, g22, g23;

		double psi_integ_thick = 0.0;
		if (composite_shell == false)
		{
			//Integração ao longo da espessura
			double csi, alpha2, jacobian, j_bar, v_, dv_,j_F,I1,psi;
			double thick = db.shell_sections[section - 1]->thickness;
			jacobian = thick / 2.0;
			Matrix tau_r1(3);
			Matrix tau_r2(3);
			for (int index = 0; index < 3; index++)
			{
				if (index == 0)
				{
					csi = -0.77459666924148337703585307995648;
					alpha2 = 0.55555555555555555555555555555556;
				}
				if (index == 1)
				{
					csi = 0.0;
					alpha2 = 0.88888888888888888888888888888889;
				}
				if (index == 2)
				{
					csi = +0.77459666924148337703585307995648;
					alpha2 = 0.55555555555555555555555555555556;
				}
				double zeta = thick * csi / 2.0;
				Matrix E3r = skew(*e3rlocal);
				gamma_r1 = *eta_r1[gauss] + zeta * cross(*kappa_r1[gauss], *e3rlocal);
				gamma_r2 = *eta_r2[gauss] + zeta * cross(*kappa_r2[gauss], *e3rlocal);
				g11 = dot(gamma_r1, *e1rlocal);
				g12 = dot(gamma_r1, *e2rlocal);
				g13 = dot(gamma_r1, *e3rlocal);
				g21 = dot(gamma_r2, *e1rlocal);
				g22 = dot(gamma_r2, *e2rlocal);
				g23 = dot(gamma_r2, *e3rlocal);


				j_bar = (1.0 + g11)*(1.0 + g22) - g12 * g21;
				v_ = (lambda*(j_bar*j_bar*j_bar - 1.0) + 2.0 * mu*(j_bar - 1.0)) / (lambda*j_bar*j_bar*j_bar + 2.0 * mu*j_bar);
				dv_ = ((lambda + 2.0 * mu)*(3.0 * lambda*j_bar*j_bar + 2.0 * mu)) / (j_bar*j_bar*(lambda*j_bar*j_bar + 2.0 * mu)*(lambda*j_bar*j_bar + 2.0 * mu));

				tau_r1(0, 0) = mu * v_*(1.0 + g22) + mu * (g11 - g22);
				tau_r1(1, 0) = mu * v_*(-g21) + mu * (g12 + g21);
				tau_r1(2, 0) = 0 + mu * g13;

				tau_r2(0, 0) = mu * v_*(-g12) + mu * (g12 + g21);
				tau_r2(1, 0) = mu * v_*(1.0 + g11) + mu * (g22 - g11);
				tau_r2(2, 0) = 0 + mu * g23;

				C_11r(0, 0) = mu * ((1.0 + g22)*(1.0 + g22)*dv_ + 1.0);
				C_11r(0, 1) = -mu * (1.0 + g22)*g21*dv_;
				C_11r(1, 0) = -mu * (1.0 + g22)*g21*dv_;
				C_11r(1, 1) = mu * (g21*g21*dv_ + 1.0);
				C_11r(2, 2) = mu;

				C_22r(0, 0) = mu * (g12*g12*dv_ + 1.0);
				C_22r(0, 1) = -mu * ((1.0 + g11)*g12*dv_);
				C_22r(1, 0) = -mu * ((1.0 + g11)*g12*dv_);
				C_22r(1, 1) = mu * ((1.0 + g11)*(1.0 + g11)*dv_ + 1.0);
				C_22r(2, 2) = mu;

				C_12r(0, 0) = -mu * ((1.0 + g22)*g12*dv_);
				C_12r(0, 1) = mu * (v_ - 1.0 + (1.0 + g11)*(1.0 + g22)*dv_);
				C_12r(1, 0) = mu * (1.0 - v_ + g12 * g21*dv_);
				C_12r(1, 1) = -mu * ((1.0 + g11)*g21*dv_);

				C_21r = transp(C_12r);

				//Integrando
				d_n1_d_eta1 = d_n1_d_eta1 + alpha2 * jacobian*C_11r;
				d_n1_d_eta2 = d_n1_d_eta2 + alpha2 * jacobian*C_12r;
				d_n1_d_kappa1 = d_n1_d_kappa1 - alpha2 * jacobian*zeta*(C_11r*E3r);
				d_n1_d_kappa2 = d_n1_d_kappa2 - alpha2 * jacobian*zeta*(C_12r*E3r);

				d_n2_d_eta1 = d_n2_d_eta1 + alpha2 * jacobian*C_21r;
				d_n2_d_eta2 = d_n2_d_eta2 + alpha2 * jacobian*C_22r;
				d_n2_d_kappa1 = d_n2_d_kappa1 - alpha2 * jacobian*zeta*(C_21r*E3r);
				d_n2_d_kappa2 = d_n2_d_kappa2 - alpha2 * jacobian*zeta*(C_22r*E3r);

				d_m1_d_eta1 = d_m1_d_eta1 + alpha2 * jacobian*zeta*(E3r*C_11r);
				d_m1_d_eta2 = d_m1_d_eta2 + alpha2 * jacobian*zeta*(E3r*C_12r);
				d_m1_d_kappa1 = d_m1_d_kappa1 - alpha2 * jacobian*zeta*zeta*(E3r*C_11r*E3r);
				d_m1_d_kappa2 = d_m1_d_kappa2 - alpha2 * jacobian*zeta*zeta*(E3r*C_12r*E3r);

				d_m2_d_eta1 = d_m2_d_eta1 + alpha2 * jacobian*zeta*(E3r*C_21r);
				d_m2_d_eta2 = d_m2_d_eta2 + alpha2 * jacobian*zeta*(E3r*C_22r);
				d_m2_d_kappa1 = d_m2_d_kappa1 - alpha2 * jacobian*zeta*zeta*(E3r*C_21r*E3r);
				d_m2_d_kappa2 = d_m2_d_kappa2 - alpha2 * jacobian*zeta*zeta*(E3r*C_22r*E3r);

				//Integrando:
				*n_r1[gauss] = *n_r1[gauss] + alpha2 * jacobian*tau_r1;
				*n_r2[gauss] = *n_r2[gauss] + alpha2 * jacobian*tau_r2;
				*m_r1[gauss] = *m_r1[gauss] + alpha2 * jacobian*zeta*(cross(*e3rlocal, tau_r1));
				*m_r2[gauss] = *m_r2[gauss] + alpha2 * jacobian*zeta*(cross(*e3rlocal, tau_r2));

				//Gamma33 - eq.75 paper Campello 2003
				double g33 = sqrt((lambda + 2.0 * mu) / (lambda * j_bar * j_bar + 2.0 * mu)) - 1.0;
				//Strain invariants 
				//J 
				j_F = j_bar * (1 + g33);
				//I1
				I1 = (1.0 + g11) * (1.0 + g11) + g12 * g12 + g13 * g13 + g21 * g21 + (1.0 + g22) * (1.0 + g22) + g23 * g23 + (1.0 + g33) * (1.0 + g33);
				//Specific strain energy (Psi)
				psi = 0.5 * lambda * (0.5 * (j_F * j_F - 1.0) - log(j_F)) + 0.5 * mu * (I1 - 3.0 - 2.0 * log(j_F));
				//Integrando na espessura
				psi_integ_thick += alpha2 * jacobian * psi;
			}
		}

		Matrix D(12, 12);			//Matriz constitutiva	
		Matrix sigma_r(12);			//Tensões generalizadas
		Matrix eta_kappa_r(12);		//Vetor das deformações

		if (composite_shell == false)
		{
			for (int i = 0; i < 3; i++)
			{
				for (int j = 0; j < 3; j++)
				{
					D(i + 0, j + 0) = d_n1_d_eta1(i, j);
					D(i + 0, j + 3) = d_n1_d_kappa1(i, j);
					D(i + 0, j + 6) = d_n1_d_eta2(i, j);
					D(i + 0, j + 9) = d_n1_d_kappa2(i, j);

					D(i + 3, j + 0) = d_m1_d_eta1(i, j);
					D(i + 3, j + 3) = d_m1_d_kappa1(i, j);
					D(i + 3, j + 6) = d_m1_d_eta2(i, j);
					D(i + 3, j + 9) = d_m1_d_kappa2(i, j);

					D(i + 6, j + 0) = d_n2_d_eta1(i, j);
					D(i + 6, j + 3) = d_n2_d_kappa1(i, j);
					D(i + 6, j + 6) = d_n2_d_eta2(i, j);
					D(i + 6, j + 9) = d_n2_d_kappa2(i, j);

					D(i + 9, j + 0) = d_m2_d_eta1(i, j);
					D(i + 9, j + 3) = d_m2_d_kappa1(i, j);
					D(i + 9, j + 6) = d_m2_d_eta2(i, j);
					D(i + 9, j + 9) = d_m2_d_kappa2(i, j);
				}
			}
			//D.fprint("D_Homogeneous.txt");
		}
		else
		{
			D = *D_comp;
			for (int i = 0; i < 3; i++)
			{
				eta_kappa_r(i + 0, 0) = (*eta_r1[gauss])(i, 0);
				eta_kappa_r(i + 3, 0) = (*kappa_r1[gauss])(i, 0);
				eta_kappa_r(i + 6, 0) = (*eta_r2[gauss])(i, 0);
				eta_kappa_r(i + 9, 0) = (*kappa_r2[gauss])(i, 0);
			}
			//(*D_comp).fprint("D_Composite.txt");
		}

		if (composite_shell == false)
		{
			//Evitando singularidade - inclusão de rigidez de drilling da casca
			D(5, 5) = stiff_drill;
			D(11, 11) = stiff_drill;
			(*m_r1[gauss])(2, 0) = stiff_drill * (*kappa_r1[gauss])(2, 0);
			(*m_r2[gauss])(2, 0) = stiff_drill * (*kappa_r2[gauss])(2, 0);

			for (int i = 0; i < 3; i++)
			{
				sigma_r(i + 0, 0) = (*n_r1[gauss])(i, 0);
				sigma_r(i + 3, 0) = (*m_r1[gauss])(i, 0);
				sigma_r(i + 6, 0) = (*n_r2[gauss])(i, 0);
				sigma_r(i + 9, 0) = (*m_r2[gauss])(i, 0);
			}
		}
		else
		{
			//Evitando singularidade - inclusão de rigidez de drilling da casca
			D(5, 5) = stiff_drill_comp1;
			D(11, 11) = stiff_drill_comp2;
			(*m_r1[gauss])(2, 0) = stiff_drill_comp1 * (*kappa_r1[gauss])(2, 0);
			(*m_r2[gauss])(2, 0) = stiff_drill_comp2 * (*kappa_r2[gauss])(2, 0);

			sigma_r = D * eta_kappa_r;
		}


		Matrix Psi(12, 15);
		Matrix tM(3, 3);
		//Preenchimento de Psi
		for (int i = 0; i<3; i++)
		for (int j = 0; j<3; j++)
			Psi(i, j) = Qtransp(i, j);
		for (int i = 0; i<3; i++)
		for (int j = 0; j<3; j++)
			Psi(i+6, j+6) = Qtransp(i, j);
		tM = Qtransp*(*Xi_delta[gauss]);
		for (int i = 0; i<3; i++)
		for (int j = 0; j<3; j++)
			Psi(i+3, j+3) = tM(i,j);
		for (int i = 0; i<3; i++)
		for (int j = 0; j<3; j++)
			Psi(i + 9, j + 9) = tM(i, j);
		tM = Qtransp*Z_x1*(*Xi_delta[gauss]);
		for (int i = 0; i<3; i++)
		for (int j = 0; j<3; j++)
			Psi(i , j + 12) = tM(i, j);
		tM = Qtransp*Xi_delta_x1;
		for (int i = 0; i<3; i++)
		for (int j = 0; j<3; j++)
			Psi(i+3, j + 12) = tM(i, j);
		tM = Qtransp*Z_x2*(*Xi_delta[gauss]);
		for (int i = 0; i<3; i++)
		for (int j = 0; j<3; j++)
			Psi(i + 6, j + 12) = tM(i, j);
		tM = Qtransp*Xi_delta_x2;
		for (int i = 0; i<3; i++)
		for (int j = 0; j<3; j++)
			Psi(i + 9, j + 12) = tM(i, j);

		//Matrix de rigidez - constribuição constitutiva
		(*stiffness) = (*stiffness) + alpha1*(transp(*deltaN[gauss])) * ((((transp(Psi)) * (D)) * (Psi)) * (*deltaN[gauss]));
		//deltaN[gauss]->fprint("deltaN.txt");
		
		//Forças e momentos internos não retro-rotacionados
		Matrix n_1 = ((*Q_delta[gauss])*(*Q_i[gauss])) * (*n_r1[gauss]);
		Matrix n_2 = ((*Q_delta[gauss])*(*Q_i[gauss])) * (*n_r2[gauss]);
		Matrix m_1 = ((*Q_delta[gauss])*(*Q_i[gauss])) * (*m_r1[gauss]);
		Matrix m_2 = ((*Q_delta[gauss])*(*Q_i[gauss])) * (*m_r2[gauss]);
		
		//Variáveis auxiliares para montagem de G
		Matrix V_alpha_Z_x1_n1 = V(*alpha_delta[gauss], Z_x1*n_1, alpha_escalar_delta);
		Matrix V_alpha_Z_x2_n2 = V(*alpha_delta[gauss], Z_x2*n_2, alpha_escalar_delta);
		Matrix V_alpha_m1 = V(*alpha_delta[gauss], m_1, alpha_escalar_delta);
		Matrix V_alpha_m2 = V(*alpha_delta[gauss], m_2, alpha_escalar_delta);
		Matrix dV_1_dalpha_m1 = d_V(*alpha_delta[gauss], *alpha_delta_x1[gauss], m_1, alpha_escalar_delta);
		Matrix dV_2_dalpha_m2 = d_V(*alpha_delta[gauss], *alpha_delta_x2[gauss], m_2, alpha_escalar_delta);
		
		Matrix G_du_alpha_1 = -1.0* skew(n_1) * (*Xi_delta[gauss]);
		Matrix G_alpha_du_1 = transp(G_du_alpha_1);
		Matrix G_du_alpha_2 = -1.0* skew(n_2) * (*Xi_delta[gauss]);
		Matrix G_alpha_du_2 = transp(G_du_alpha_2);
		
		Matrix G_alpha_dalpha_1 = (V_alpha_m1);
		Matrix G_dalpha_alpha_1 = transp(G_alpha_dalpha_1);
		
		Matrix G_alpha_dalpha_2 = (V_alpha_m2);
		Matrix G_dalpha_alpha_2 = transp(G_alpha_dalpha_2);
		
		Matrix G_alpha_alpha_1 = +1.0*(transp(*Xi_delta[gauss])) * (Z_x1* skew(n_1)) * (*Xi_delta[gauss]) - 1.0*V_alpha_Z_x1_n1 + dV_1_dalpha_m1 - transp(Xi_delta_x1) * (skew(m_1) * (*Xi_delta[gauss]));
		Matrix G_alpha_alpha_2 = +1.0*(transp(*Xi_delta[gauss])) * (Z_x2* skew(n_2)) * (*Xi_delta[gauss]) - 1.0*V_alpha_Z_x2_n2 + dV_2_dalpha_m2 - transp(Xi_delta_x2) * (skew(m_2) * (*Xi_delta[gauss]));
		
		//Cálculo da matriz G
		Matrix G(15, 15);
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				G(i + 0, j + 12) = G_du_alpha_1(i, j);
				G(i + 3, j + 12) = G_dalpha_alpha_1(i, j);
				G(i + 6, j + 12) = G_du_alpha_2(i, j);
				G(i + 9, j + 12) = G_dalpha_alpha_2(i, j);
				G(i + 12, j + 0) = G_alpha_du_1(i, j);
				G(i + 12, j + 3) = G_alpha_dalpha_1(i, j);
				G(i + 12, j + 6) = G_alpha_du_2(i, j);
				G(i + 12, j + 9) = G_alpha_dalpha_2(i, j);
				G(i + 12, j + 12) = G_alpha_alpha_1(i, j) + G_alpha_alpha_2(i, j);
			}
		}
		//Matrix de rigidez - constribuição geométrica
		(*stiffness) = (*stiffness) + alpha1*(transp(*deltaN[gauss]))*((G)*(*deltaN[gauss]));
		//Vetor de esforços internos
		(*i_loading) = (*i_loading) + alpha1*(transp(*deltaN[gauss]) * transp(Psi) * sigma_r);
		//Energia de deformação
		strain_energy += (alpha1 * psi_integ_thick);
		
	}
	//Transformações de coordenadas, para o sistema global - matriz de rigidez e esforços internos
	(*stiffness) = (transp(*transform)*(*stiffness))*(*transform);
	(*i_loading) = transp(*transform)*(*i_loading);
}

//Monta carregamentos de campo (ex: peso próprio)
void Shell_1::MountFieldLoads()
{
	if (db.environment != NULL)
	{
		//Se existe campo gravitacional
		if (db.environment->g_exist == true)
		{
			load_multiplier = 1.0;
			l_factor = db.environment->bool_g.GetLinearFactorAtCurrentTime();

			//Fator multiplicativo para o carregamento (contém o fator, rho e a espessura das casca) = fator * peso por unidade de área
			mult = l_factor*load_multiplier*(db.materials[material - 1]->rho)*db.shell_sections[section-1]->thickness;
			Matrix q_weight(6);
			q_weight(0, 0) = mult*db.environment->G(0, 0);
			q_weight(1, 0) = mult*db.environment->G(1, 0);
			q_weight(2, 0) = mult*db.environment->G(2, 0);
			//Loop nos pontos de Gauss
			for (int gauss = 0; gauss < 6; gauss++)
				(*e_loading) = (*e_loading) + w4[gauss] * transp(*N4[gauss])*q_weight;
		}

		//Se existe campo gravitacional
		if (db.environment->g_exist == true)
		{
			load_multiplier = 1.0;
			l_factor = db.environment->bool_g.GetLinearFactorAtCurrentTime();

			//Fator multiplicativo para o carregamento (contém o fator, rho e a espessura das casca) = fator * peso por unidade de área
			if (composite_shell == false)	
				mult = l_factor * load_multiplier*(db.materials[material - 1]->rho)*db.shell_sections[section - 1]->thickness;
			else
				mult = l_factor * load_multiplier*(rho_comp)*thick_comp;
			
			Matrix q_weight(6);
			q_weight(0, 0) = mult * db.environment->G(0, 0);
			q_weight(1, 0) = mult * db.environment->G(1, 0);
			q_weight(2, 0) = mult * db.environment->G(2, 0);
			//Loop nos pontos de Gauss
			for (int gauss = 0; gauss < 6; gauss++)
				(*e_loading) = (*e_loading) + w4[gauss] * transp(*N4[gauss])*q_weight;
		}


		
	}

}

//Monta carregamentos associados ao elemento
void Shell_1::MountElementLoads()
{
	MountFieldLoads();
	//Chamar aqui a função para montar esforços externos consistentes da casca (ex. peso próprio)
	(*P_loading) = (*i_loading) - (*e_loading);		//Vetor esforço desbalanceado
}

//Monta carregamentos de pressão na casca
void Shell_1::MountShellSpecialLoads(int l_number)
{
	//////////////////Efeito da Pressão////////////////////////////
	ShellLoad* ptr = static_cast<ShellLoad*>(db.loads[l_number - 1]);
	double pressure = ptr->GetValueAt(db.last_converged_time + db.current_time_step, 0);
	bool area_update = ptr->area_update;
	mult = 1.0;
	
	//Esforço externo generalizado por unidade de área de referência
	Matrix q_pres(6);
	//Direções tangentes
	Matrix t1_p(3);
	Matrix t2_p(3);
	//Direção normal
	Matrix n_p(3);
	//Graus de liberdade
	Matrix pu_delta(18);
	Matrix pu_i(18);
	Matrix pu_ip(18);
	//Preenchimento com graus de liberdade						
	for (int i = 0; i < 6; i++)
	{
		pu_delta(3*i + 0, 0) = db.nodes[nodes[i] - 1]->displacements[0];
		pu_delta(3*i + 1, 0) = db.nodes[nodes[i] - 1]->displacements[1];
		pu_delta(3*i + 2, 0) = db.nodes[nodes[i] - 1]->displacements[2];

		pu_i(3 * i + 0, 0) = db.nodes[nodes[i] - 1]->copy_coordinates[0] - db.nodes[nodes[i] - 1]->ref_coordinates[0];
		pu_i(3 * i + 1, 0) = db.nodes[nodes[i] - 1]->copy_coordinates[1] - db.nodes[nodes[i] - 1]->ref_coordinates[1];
		pu_i(3 * i + 2, 0) = db.nodes[nodes[i] - 1]->copy_coordinates[2] - db.nodes[nodes[i] - 1]->ref_coordinates[2];
	}
	pu_ip = pu_i + pu_delta;
	//Loop nos pontos de Gauss
	for (int gauss = 0; gauss < 6; gauss++)
	{
		//Direções tangentes
		t1_p = *e1r + (*N4_u_x1[gauss])*pu_ip;
		t2_p = *e2r + (*N4_u_x2[gauss])*pu_ip;
		//Direção normal - já normalizada
		n_p = (1.0 / norm(cross(t1_p, t2_p)))*cross(t1_p, t2_p);
		//Carregamento
		q_pres(0, 0) = -mult*pressure*n_p(0, 0);
		q_pres(1, 0) = -mult*pressure*n_p(1, 0);
		q_pres(2, 0) = -mult*pressure*n_p(2, 0);

		//Caso 1 - area update false
		if (area_update == false)
		{
			(*P_loading) = (*P_loading) - w4[gauss] * transp(*N4[gauss])*q_pres;	//esforço desbalanceado
			(*e_loading) = (*e_loading) + w4[gauss] * transp(*N4[gauss])*q_pres;	//carregamento externo
			//Contribuição para a matriz de rigidez
			Matrix Kpaa(3, 18);
			Matrix Kp(6, 27);
			Kpaa = (1.0 / norm(cross(t1_p, t2_p)))*(I3 - dyadic(n_p, n_p))*(skew(t1_p)*(*N4_u_x2[gauss]) - skew(t2_p)*(*N4_u_x1[gauss]));
			for (int i = 0; i < 3; i++)
			for (int j = 0; j < 18; j++)
				Kp(i, j) = Kpaa(i, j);
			(*stiffness) = (*stiffness) + w4[gauss] * mult*pressure*transp(*N4[gauss])*Kp;
		}
		//Caso 2 - area update true
		else
		{
			double jac = norm(cross(t1_p, t2_p));
			(*P_loading) = (*P_loading) - w4[gauss] * jac * transp(*N4[gauss])*q_pres;	//esforço desbalanceado
			(*e_loading) = (*e_loading) + w4[gauss] * jac * transp(*N4[gauss])*q_pres;	//carregamento externo
			//Contribuição para a matriz de rigidez
			Matrix Kpaa(3, 18);
			Matrix Kp(6, 27);
			Kpaa = (skew(t1_p)*(*N4_u_x2[gauss]) - skew(t2_p)*(*N4_u_x1[gauss]));
			for (int i = 0; i < 3; i++)
			for (int j = 0; j < 18; j++)
				Kp(i, j) = Kpaa(i, j);
			(*stiffness) = (*stiffness) + w4[gauss] * mult*pressure*transp(*N4[gauss])*Kp;
		}
		
	}
}

//Monta matriz de transformação de coordenadas
void Shell_1::TransformMatrix()
{
	Matrix e1(3);
	e1(0, 0) = 1.0;
	Matrix e2(3);
	e2(1, 0) = 1.0;
	Matrix e3(3);
	e3(2, 0) = 1.0;
	//Preenche a matriz de transformação de coordenadas
	for (int i = 0; i<27; i = i + 3)
	{
		//Preenche também a matriz transform3
		if (i == 0)
		{
			(*transform3)(0, 0) = dot(*e1r, e1);
			(*transform3)(0, 1) = dot(*e1r, e2);
			(*transform3)(0, 2) = dot(*e1r, e3);

			(*transform3)(1, 0) = dot(*e2r, e1);
			(*transform3)(1, 1) = dot(*e2r, e2);
			(*transform3)(1, 2) = dot(*e2r, e3);

			(*transform3)(2, 0) = dot(*e3r, e1);
			(*transform3)(2, 1) = dot(*e3r, e2);
			(*transform3)(2, 2) = dot(*e3r, e3);
		}

		(*transform)(0 + i, 0 + i) = dot(*e1r, e1);
		(*transform)(0 + i, 1 + i) = dot(*e1r, e2);
		(*transform)(0 + i, 2 + i) = dot(*e1r, e3);

		(*transform)(1 + i, 0 + i) = dot(*e2r, e1);
		(*transform)(1 + i, 1 + i) = dot(*e2r, e2);
		(*transform)(1 + i, 2 + i) = dot(*e2r, e3);

		(*transform)(2 + i, 0 + i) = dot(*e3r, e1);
		(*transform)(2 + i, 1 + i) = dot(*e3r, e2);
		(*transform)(2 + i, 2 + i) = dot(*e3r, e3);
	}
	//transform->fprint("d_transf.txt");
	
}
//Preenche a contribuição do elemento nas matrizes globais
void Shell_1::MountGlobal()
{
	//Variáveis temporárias para salvar a indexação global dos graus de liberdade
	int GL_global_1 = 0;
	int GL_global_2 = 0;
	double anterior = 0;
	for (int i = 0; i < 27; i++)
	{
		//////////////MONTAGEM DO VETOR DE ESFORÇOS DESBALANCEADOS//////////////////
		//Toma do vetor de GL globais, a indexação de cada grau de liberdade global
		if (i < 3)//node 1 - displacements
			GL_global_1 = db.nodes[nodes[0] - 1]->GLs[i];
		else	
		{
			if (i < 6)//node 2 - displacements
				GL_global_1 = db.nodes[nodes[1] - 1]->GLs[i - 3];
			else
			{
				if (i < 9)//node 3 - displacements
					GL_global_1 = db.nodes[nodes[2] - 1]->GLs[i - 6];
				else
				{
					if (i < 12)//node 4 - displacements
						GL_global_1 = db.nodes[nodes[3] - 1]->GLs[i - 9];
					else
					{
						if (i < 15)//node 5 - displacements
							GL_global_1 = db.nodes[nodes[4] - 1]->GLs[i - 12];
						else
						{
							if (i < 18)//node 6 - displacements
								GL_global_1 = db.nodes[nodes[5] - 1]->GLs[i - 15];
							else
							{
								if (i < 21)//node 4 - rotations
									GL_global_1 = db.nodes[nodes[3] - 1]->GLs[i - 18 + 3];
								else
								{
									if (i < 24)//node 5 - rotations
										GL_global_1 = db.nodes[nodes[4] - 1]->GLs[i - 21 + 3];
									else      
									{
										//node 6 - rotations
										GL_global_1 = db.nodes[nodes[5] - 1]->GLs[i - 24 + 3];
									}
								}
							}
						}
					}
				}
			}
		}

		//Caso o grau de liberdade seja livre:
		if (GL_global_1 > 0)
		{
			anterior = db.global_P_A(GL_global_1 - 1, 0);
			db.global_P_A(GL_global_1 - 1, 0) = anterior + (*P_loading)(i,0);
			anterior = db.global_I_A(GL_global_1 - 1, 0);
			db.global_I_A(GL_global_1 - 1, 0) = anterior + (*P_loading)(i, 0);
		}
		else
		{
			anterior = db.global_P_B(-GL_global_1 - 1, 0);
			db.global_P_B(-GL_global_1 - 1, 0) = anterior + (*P_loading)(i, 0);
		}
		for (int j = 0; j < 27; j++)
		{
			//////////////////////MONTAGEM DA MATRIZ DE RIGIDEZ/////////////////////////
			//Toma do vetor de GL globais, a indexação de cada grau de liberdade global
			if (j < 3)//node 1 - displacements
				GL_global_2 = db.nodes[nodes[0] - 1]->GLs[j];
			else	
			{
				if (j < 6)//node 2 - displacements
					GL_global_2 = db.nodes[nodes[1] - 1]->GLs[j - 3];
				else
				{
					if (j < 9)//node 3 - displacements
						GL_global_2 = db.nodes[nodes[2] - 1]->GLs[j - 6];
					else
					{
						if (j < 12)//node 4 - displacements
							GL_global_2 = db.nodes[nodes[3] - 1]->GLs[j - 9];
						else
						{
							if (j < 15)//node 5 - displacements
								GL_global_2 = db.nodes[nodes[4] - 1]->GLs[j - 12];
							else
							{
								if (j < 18)//node 6 - displacements
									GL_global_2 = db.nodes[nodes[5] - 1]->GLs[j - 15];
								else
								{
									if (j < 21)//node 4 - rotations
										GL_global_2 = db.nodes[nodes[3] - 1]->GLs[j - 18 + 3];
									else
									{
										if (j < 24)//node 5 - rotations
											GL_global_2 = db.nodes[nodes[4] - 1]->GLs[j - 21 + 3];
										else
										{
											//node 6 - rotations
											GL_global_2 = db.nodes[nodes[5] - 1]->GLs[j - 24 + 3];
										}
									}
								}
							}
						}
					}
				}
			}

			//Caso os graus de liberdade sejam ambos livres (Matriz Kaa)
			if (GL_global_1 > 0 && GL_global_2 > 0)
			{
				db.global_stiffness_AA.setValue(GL_global_1 - 1, GL_global_2 - 1, (*stiffness)(i,j));
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
		}//segundo loop
	}//primeiro loop
}
//Salva variáveis nos pontos de Gauss úteis para descrição lagrangiana atualizada
void Shell_1::SaveLagrange()
{
	//Loop nos pontos de Gauss
	for (int gauss = 0; gauss < 3; gauss++)
	{
		*kappa_r1_i[gauss] = transp(*Q_i[gauss])*transp(*Xi_delta[gauss])*(*alpha_delta_x1[gauss]) + *kappa_r1_i[gauss];
		*kappa_r2_i[gauss] = transp(*Q_i[gauss])*transp(*Xi_delta[gauss])*(*alpha_delta_x2[gauss]) + *kappa_r2_i[gauss];
		*Q_i[gauss] = (*Q_delta[gauss])*(*Q_i[gauss]);														//Tensor de rotação
		*u_i[gauss] = (*u_i[gauss]) + *u_delta[gauss];														//Deslocamento
		(*alpha_i[gauss]) = 4.0 / (4.0 - dot(*alpha_delta[gauss], *alpha_i[gauss]))* (*alpha_delta[gauss] +
			(*alpha_i[gauss]) + 0.5*cross(*alpha_delta[gauss], *alpha_i[gauss]));							//Vetor rotação de Rodrigues
		(*z_x1_i[gauss]) = *u_delta_x1[gauss] + *z_x1_i[gauss];
		(*z_x2_i[gauss]) = *u_delta_x2[gauss] + *z_x2_i[gauss];
	}
}
//Pré-cálculo de variáveis que é feito uma única vez no início
void Shell_1::PreCalc()
{
	Matrix xp(3);		//Posição do ponto de integração
	Matrix x1(3);		//Posição do nó 1
	Matrix x2(3);		//Posição do nó 2
	Matrix x3(3);		//Posição do nó 3
	Matrix x4(3);		//Posição do nó 4
	Matrix x5(3);		//Posição do nó 5
	Matrix x6(3);		//Posição do nó 6
	
	//Preenchendo as coordenadas - sistema global
	for (int coord = 0; coord < 3; coord++)
	{
		x1(coord, 0) = db.nodes[nodes[0] - 1]->ref_coordinates[coord];
		x2(coord, 0) = db.nodes[nodes[1] - 1]->ref_coordinates[coord];
		x3(coord, 0) = db.nodes[nodes[2] - 1]->ref_coordinates[coord];
		x4(coord, 0) = db.nodes[nodes[3] - 1]->ref_coordinates[coord];
		x5(coord, 0) = db.nodes[nodes[4] - 1]->ref_coordinates[coord];
		x6(coord, 0) = db.nodes[nodes[5] - 1]->ref_coordinates[coord];
	}
	
	double L1, L2, L3, A1, A2, A3, A;
	A = 0.5 * norm(cross(x2 - x1, x3 - x1));
	area = A;

	//Composite Shell
	composite_shell = false;
	if (typeid(*db.shell_sections[section - 1]) == typeid(ShellSectionComposite))
	{
		composite_shell = true;
		//Material - cálculo das propriedades do material compósito			
		ShellSectionComposite* shell_section = static_cast<ShellSectionComposite*>(db.shell_sections[section - 1]);
		int n_lam = shell_section->n_laminas; //Número de lâminas do compósito
		//Matriz com os dados de cálculo da matriz ABD (identificador, espessura e ângulo)
		Lamina *db_lam = shell_section->db_laminas;

		//Elementos da Matriz ABD de um laminado compósito
		double Axx = 0;
		double Ayy = 0;
		double Axy = 0;
		double Axs = 0;
		double Ays = 0;
		double Ass = 0;

		double Bxx = 0;
		double Byy = 0;
		double Bxy = 0;
		double Bxs = 0;
		double Bys = 0;
		double Bss = 0;

		double Dxx = 0;
		double Dyy = 0;
		double Dxy = 0;
		double Dxs = 0;
		double Dys = 0;
		double Dss = 0;

		double Axzxz = 0;
		double Ayzyz = 0;
		double Axzyz = 0;

		//Espessura total do compósito
		for (int i = 0; i < n_lam; i++)
		{
			thick_comp += db_lam[i].thickness;
		}

		double hi;
		double hf;
		double E1_comp = 0;
		double E2_comp = 0;

		for (int i = 0; i < n_lam; i++)
		{
			Orthotropic* ortho = static_cast<Orthotropic*>(db.materials[db_lam[i].id - 1]);

			double E1 = ortho->E1;
			double E2 = ortho->E2;
			double G12 = ortho->G12;
			double G23 = ortho->G23;
			double nu12 = ortho->nu12;
			double nu21 = nu12 * (E2 / E1);
			double rho = ortho->rho;

			//Rho Equivalente
			rho_comp += rho * db_lam[i].thickness / thick_comp;

			//Menor "E" medio
			if (i == 0)
			{
				E1_comp = E1;
				E2_comp = E2;
			}
			else
			{
				if (E1 < E1_comp)
					E1_comp = E1;
				if (E2 < E2_comp)
					E2_comp = E2;
			}

			//Cálculo da Matriz Q
			double Q11 = E1 / (1 - nu12 * nu21);
			double Q22 = E2 / (1 - nu12 * nu21);
			double Q12 = (E1*nu21) / (1 - nu12 * nu21);
			double Q66 = G12;
			double Q55 = G12;
			double Q44 = G23;

			//Cálculo da Matriz Q Rotacionada
			double m = cos(db_lam[i].angle * PI / 180);
			double n = sin(db_lam[i].angle * PI / 180);

			double m2 = pow(m, 2);
			double m3 = pow(m, 3);
			double m4 = pow(m, 4);

			double n2 = pow(n, 2);
			double n3 = pow(n, 3);
			double n4 = pow(n, 4);

			double Qxx = m4 * Q11 + n4 * Q22 + 2 * m2*n2*Q12 + 4 * m2*n2*Q66;
			double Qyy = n4 * Q11 + m4 * Q22 + 2 * m2*n2*Q12 + 4 * m2*n2*Q66;
			double Qxy = m2 * n2*Q11 + m2 * n2*Q22 + (m4 + n4)*Q12 - 4 * m2*n2*Q66;

			double Qxs = m3 * n*Q11 - m * n3*Q22 - m * n*(m2 - n2)*Q12 - 2 * m*n*(m2 - n2)*Q66;
			double Qys = m * n3*Q11 - m3 * n*Q22 + m * n*(m2 - n2)*Q12 + 2 * m*n*(m2 - n2)*Q66;
			double Qss = m2 * n2*Q11 + m2 * n2*Q22 - 2 * m2*n2*Q12 + (m2 - n2)*(m2 - n2)*Q66;

			double Qxzxz = m2 * Q55 + n2 * Q44;
			double Qyzyz = n2 * Q55 + m2 * Q44;
			double Qxzyz = m * n*(Q55 - Q44);

			//Coordenada z para cada lâmina


			if (i == 0) {
				hi = -thick_comp / 2;
			}
			else {
				hi = hf;
			}

			hf = hi + db_lam[i].thickness;

			//Coeficientes da matriz ABD
			double CA = hf - hi;
			double CB = (pow(hf, 2) - pow(hi, 2)) / 2;
			double CD = (pow(hf, 3) - pow(hi, 3)) / 3;

			//Cáclulo da Matriz ABD
			Axx += Qxx * CA;
			Ayy += Qyy * CA;
			Axy += Qxy * CA;
			Axs += Qxs * CA;
			Ays += Qys * CA;
			Ass += Qss * CA;

			Bxx += Qxx * CB;
			Byy += Qyy * CB;
			Bxy += Qxy * CB;
			Bxs += Qxs * CB;
			Bys += Qys * CB;
			Bss += Qss * CB;

			Dxx += Qxx * CD;
			Dyy += Qyy * CD;
			Dxy += Qxy * CD;
			Dxs += Qxs * CD;
			Dys += Qys * CD;
			Dss += Qss * CD;

			Axzxz += Qxzxz * CA;
			Ayzyz += Qyzyz * CA;
			Axzyz += Qxzyz * CA;
		}

		//Matriz Constitutiva Compósito
		int ni = 0; //contador auxiliar das linhas
		int nj = 0; //contador auxiliar das colunas


		//dn1_deta1
		(*D_comp)(0 + ni, 0 + nj) = Axx;
		(*D_comp)(1 + ni, 0 + nj) = Axs;
		(*D_comp)(0 + ni, 1 + nj) = Axs;
		(*D_comp)(1 + ni, 1 + nj) = Ass;
		(*D_comp)(2 + ni, 2 + nj) = Axzxz;

		nj = 3;
		//dn1_dkappa1
		(*D_comp)(0 + ni, 0 + nj) = -Bxs;
		(*D_comp)(1 + ni, 0 + nj) = -Bss;
		(*D_comp)(0 + ni, 1 + nj) = Bxx;
		(*D_comp)(1 + ni, 1 + nj) = Bxs;

		nj = 6;
		//dn1_deta2
		(*D_comp)(0 + ni, 0 + nj) = Axs;
		(*D_comp)(1 + ni, 0 + nj) = Ass;
		(*D_comp)(0 + ni, 1 + nj) = Axy;
		(*D_comp)(1 + ni, 1 + nj) = Ays;
		(*D_comp)(2 + ni, 2 + nj) = Axzyz;

		nj = 9;
		//dn1_dkappa2
		(*D_comp)(0 + ni, 0 + nj) = -Bxy;
		(*D_comp)(1 + ni, 0 + nj) = -Bys;
		(*D_comp)(0 + ni, 1 + nj) = Bxs;
		(*D_comp)(1 + ni, 1 + nj) = Bss;

		ni = 3;
		nj = 0;
		//dm1_deta1
		(*D_comp)(0 + ni, 0 + nj) = -Bxs;
		(*D_comp)(1 + ni, 0 + nj) = Bxx;
		(*D_comp)(0 + ni, 1 + nj) = -Bss;
		(*D_comp)(1 + ni, 1 + nj) = Bxs;

		nj = 3;
		//dm1_dkappa1
		(*D_comp)(0 + ni, 0 + nj) = Dss;
		(*D_comp)(1 + ni, 0 + nj) = -Dxs;
		(*D_comp)(0 + ni, 1 + nj) = -Dxs;
		(*D_comp)(1 + ni, 1 + nj) = Dxx;

		nj = 6;
		//dm1_deta2
		(*D_comp)(0 + ni, 0 + nj) = -Bss;
		(*D_comp)(1 + ni, 0 + nj) = Bxs;
		(*D_comp)(0 + ni, 1 + nj) = -Bys;
		(*D_comp)(1 + ni, 1 + nj) = Bxy;

		nj = 9;
		//dm1_dkappa2
		(*D_comp)(0 + ni, 0 + nj) = Dys;
		(*D_comp)(1 + ni, 0 + nj) = -Dxy;
		(*D_comp)(0 + ni, 1 + nj) = -Dss;
		(*D_comp)(1 + ni, 1 + nj) = Dxs;

		ni = 6;
		nj = 0;
		//dn2_deta1
		(*D_comp)(0 + ni, 0 + nj) = Axs;
		(*D_comp)(1 + ni, 0 + nj) = Axy;
		(*D_comp)(0 + ni, 1 + nj) = Ass;
		(*D_comp)(1 + ni, 1 + nj) = Ays;
		(*D_comp)(2 + ni, 2 + nj) = Axzyz;

		nj = 3;
		//dn2_dkappa1
		(*D_comp)(0 + ni, 0 + nj) = -Bss;
		(*D_comp)(1 + ni, 0 + nj) = -Bys;
		(*D_comp)(0 + ni, 1 + nj) = Bxs;
		(*D_comp)(1 + ni, 1 + nj) = Bxy;

		nj = 6;
		//dn2_deta2
		(*D_comp)(0 + ni, 0 + nj) = Ass;
		(*D_comp)(1 + ni, 0 + nj) = Ays;
		(*D_comp)(0 + ni, 1 + nj) = Ays;
		(*D_comp)(1 + ni, 1 + nj) = Ayy;
		(*D_comp)(2 + ni, 2 + nj) = Ayzyz;

		nj = 9;
		//dn2_dkappa2
		(*D_comp)(0 + ni, 0 + nj) = -Bys;
		(*D_comp)(1 + ni, 0 + nj) = -Byy;
		(*D_comp)(0 + ni, 1 + nj) = Bss;
		(*D_comp)(1 + ni, 1 + nj) = Bys;

		ni = 9;
		nj = 0;
		//dm2_deta1
		(*D_comp)(0 + ni, 0 + nj) = -Bxy;
		(*D_comp)(1 + ni, 0 + nj) = Bxs;
		(*D_comp)(0 + ni, 1 + nj) = -Bys;
		(*D_comp)(1 + ni, 1 + nj) = Bss;

		nj = 3;
		//dm2_dkappa1
		(*D_comp)(0 + ni, 0 + nj) = Dys;
		(*D_comp)(1 + ni, 0 + nj) = -Dss;
		(*D_comp)(0 + ni, 1 + nj) = -Dxy;
		(*D_comp)(1 + ni, 1 + nj) = Dxs;

		nj = 6;
		//dm2_deta2
		(*D_comp)(0 + ni, 0 + nj) = -Bys;
		(*D_comp)(1 + ni, 0 + nj) = Bss;
		(*D_comp)(0 + ni, 1 + nj) = -Byy;
		(*D_comp)(1 + ni, 1 + nj) = Bys;

		nj = 9;
		//dm2_dkappa2
		(*D_comp)(0 + ni, 0 + nj) = Dyy;
		(*D_comp)(1 + ni, 0 + nj) = -Dys;
		(*D_comp)(0 + ni, 1 + nj) = -Dys;
		(*D_comp)(1 + ni, 1 + nj) = Dss;

		////Material - cálculo das constantes para material compósito
		//stiff_drill = E_comp * thick_comp*thick_comp*thick_comp;	//Rigidez ao drilling
		stiff_drill_comp1 = E1_comp * thick_comp*thick_comp*thick_comp;
		stiff_drill_comp2 = E2_comp * thick_comp*thick_comp*thick_comp;
		coef1 = thick_comp * rho_comp;
		coef2 = (1.0 / 12.0)*thick_comp*thick_comp*thick_comp*rho_comp;
		coef3 = rho_comp * thick_comp*area / (3 * PI);		//rotação de drilling, por unidade de área
	}


	if (composite_shell == false)
	{
		//Criação de um sistema local - e1r alinhado com a direção x2 - x1 e e2r ortogonal
		//*e1r = (1.0 / norm(x2 - x1))*(x2 - x1);
		//*e3r = (1.0 / norm(cross(x2 - x1, x3 - x1)))*cross(x2 - x1, x3 - x1);
		//*e2r = cross(*e3r,*e1r);

		//Criação de um sistema local - primeiro é definido o e3r (normal à casca)
		*e3r = (1.0 / norm(cross(x2 - x1, x3 - x1)))*cross(x2 - x1, x3 - x1);
		//O vetor e1r é a projeção do x global no plano da casca. Se for projeção nula, será a projeção do y global no plano da casca
		Matrix eg(3);
		double eps = 1e-4;
		eg(0, 0) = 1.0;	//X
		if (abs(dot(eg, *e3r)) >= 1.0 - eps)//Se houver alinhamento entre e3r e o eixo x, então define o eixo global como sendo o y
		{
			eg(0, 0) = 0.0;
			eg(1, 0) = 1.0;
		}
		*e1r = eg - dot(eg, *e3r)*(*e3r);
		*e1r = (1.0 / norm(*e1r))*(*e1r);//normatização
		*e2r = cross(*e3r, *e1r);
	}
	else
	{
		*e1r = *db.CS[cs - 1]->E1;
		*e2r = *db.CS[cs - 1]->E2;
		*e3r = *db.CS[cs - 1]->E3;
	}

	(*e1rlocal)(0, 0) = 1.0;
	(*e2rlocal)(1, 0) = 1.0;
	(*e3rlocal)(2, 0) = 1.0;

	//Monta as matrizes de transformação de coordenadas
	TransformMatrix();

	if (composite_shell == false)
	{
		//Material - cálculo das constantes de Lamé
		Hooke* hooke = static_cast<Hooke*>(db.materials[material - 1]);
		double E = hooke->E;
		double nu = hooke->nu;
		mu = E / (2.0 * (1 + nu));
		lambda = 2.0*mu*nu / (1 - 2.0*nu);
		stiff_drill = E * db.shell_sections[section - 1]->thickness*db.shell_sections[section - 1]->thickness*db.shell_sections[section - 1]->thickness;	//Rigidez ao drilling
		coef1 = db.shell_sections[section - 1]->thickness*hooke->rho;
		coef2 = (1.0 / 12.0)*db.shell_sections[section - 1]->thickness*db.shell_sections[section - 1]->thickness*db.shell_sections[section - 1]->thickness*hooke->rho;
		coef3 = hooke->rho*db.shell_sections[section - 1]->thickness*area / (3 * PI);		//rotação de drilling, por unidade de área
	}

	
	//Percorre pontos de Gauss para salvar algumas matrizes/vetores de interesse nos cálculos do elemento
	for (int gauss = 0; gauss < 3; gauss++)
	{
		//Ponto localizado sobre o nó 4
		if (gauss == 0)
			xp = x4;
		//Ponto localizado sobre o nó 5
		if (gauss == 1)
			xp = x5;
		//Ponto localizado sobre o nó 6
		if (gauss == 2)
			xp = x6;
		//Cálculo das funções de forma e de suas derivadas
		A1 = 0.5 * norm(cross(x2 - xp, x3 - xp));
		A2 = 0.5 * norm(cross(x3 - xp, x1 - xp));
		A3 = 0.5 * norm(cross(x1 - xp, x2 - xp));
		L1 = A1 / A;
		L2 = A2 / A;
		L3 = A3 / A;
		double b1, b2, b3, c1, c2, c3;
		b1 = dot(x2 - x3, *e2r);
		b2 = dot(x3 - x1, *e2r);
		b3 = dot(x1 - x2, *e2r);
		c1 = dot(x3 - x2, *e1r);
		c2 = dot(x1 - x3, *e1r);
		c3 = dot(x2 - x1, *e1r);
		double L1_x1, L2_x1, L3_x1, L1_x2, L2_x2, L3_x2;
		L1_x1 = 0.5 * b1 / A;
		L2_x1 = 0.5 * b2 / A;
		L3_x1 = 0.5 * b3 / A;
		L1_x2 = 0.5 * c1 / A;
		L2_x2 = 0.5 * c2 / A;
		L3_x2 = 0.5 * c3 / A;
		//Atribuindo os valores das funções de forma e derivadas
		N1u[gauss] = (2 * L1 - 1)*L1;
		N2u[gauss] = (2 * L2 - 1)*L2;
		N3u[gauss] = (2 * L3 - 1)*L3;
		N4u[gauss] = 4 * L1 * L2;
		N5u[gauss] = 4 * L2 * L3;
		N6u[gauss] = 4 * L3 * L1;
		N4a[gauss] = 1 - 2 * L3;
		N5a[gauss] = 1 - 2 * L1;
		N6a[gauss] = 1 - 2 * L2;
		N1u_x1[gauss] = 4 * L1_x1*L1 - L1_x1;
		N2u_x1[gauss] = 4 * L2_x1*L2 - L2_x1;
		N3u_x1[gauss] = 4 * L3_x1*L3 - L3_x1;
		N4u_x1[gauss] = 4 * L1_x1*L2 + 4 * L1*L2_x1;
		N5u_x1[gauss] = 4 * L2_x1*L3 + 4 * L2*L3_x1;
		N6u_x1[gauss] = 4 * L3_x1*L1 + 4 * L3*L1_x1;
		N4a_x1[gauss] = -2 * L3_x1;
		N5a_x1[gauss] = -2 * L1_x1;
		N6a_x1[gauss] = -2 * L2_x1;
		N1u_x2[gauss] = 4 * L1_x2*L1 - L1_x2;
		N2u_x2[gauss] = 4 * L2_x2*L2 - L2_x2;
		N3u_x2[gauss] = 4 * L3_x2*L3 - L3_x2;
		N4u_x2[gauss] = 4 * L1_x2*L2 + 4 * L1*L2_x2;
		N5u_x2[gauss] = 4 * L2_x2*L3 + 4 * L2*L3_x2;
		N6u_x2[gauss] = 4 * L3_x2*L1 + 4 * L3*L1_x2;
		N4a_x2[gauss] = -2 * L3_x2;
		N5a_x2[gauss] = -2 * L1_x2;
		N6a_x2[gauss] = -2 * L2_x2;
		//Preenchendo a matriz N (funções de forma):
		(*N[gauss])(0, 0 + 0) = N1u[gauss];
		(*N[gauss])(1, 1 + 0) = N1u[gauss];
		(*N[gauss])(2, 2 + 0) = N1u[gauss];
		(*N[gauss])(0, 0 + 3) = N2u[gauss];
		(*N[gauss])(1, 1 + 3) = N2u[gauss];
		(*N[gauss])(2, 2 + 3) = N2u[gauss];
		(*N[gauss])(0, 0 + 6) = N3u[gauss];
		(*N[gauss])(1, 1 + 6) = N3u[gauss];
		(*N[gauss])(2, 2 + 6) = N3u[gauss];
		(*N[gauss])(0, 0 + 9) = N4u[gauss];
		(*N[gauss])(1, 1 + 9) = N4u[gauss];
		(*N[gauss])(2, 2 + 9) = N4u[gauss];
		(*N[gauss])(0, 0 + 12) = N5u[gauss];
		(*N[gauss])(1, 1 + 12) = N5u[gauss];
		(*N[gauss])(2, 2 + 12) = N5u[gauss];
		(*N[gauss])(0, 0 + 15) = N6u[gauss];
		(*N[gauss])(1, 1 + 15) = N6u[gauss];
		(*N[gauss])(2, 2 + 15) = N6u[gauss];
		(*N[gauss])(3, 0 + 18) = N4a[gauss];
		(*N[gauss])(4, 1 + 18) = N4a[gauss];
		(*N[gauss])(5, 2 + 18) = N4a[gauss];
		(*N[gauss])(3, 0 + 21) = N5a[gauss];
		(*N[gauss])(4, 1 + 21) = N5a[gauss];
		(*N[gauss])(5, 2 + 21) = N5a[gauss];
		(*N[gauss])(3, 0 + 24) = N6a[gauss];
		(*N[gauss])(4, 1 + 24) = N6a[gauss];
		(*N[gauss])(5, 2 + 24) = N6a[gauss];
		//Prenchendo a matriz deltaN (com funções de forma e suas derivadas)
		(*deltaN[gauss])(0, 0 + 0) = N1u_x1[gauss];
		(*deltaN[gauss])(1, 1 + 0) = N1u_x1[gauss];
		(*deltaN[gauss])(2, 2 + 0) = N1u_x1[gauss];
		(*deltaN[gauss])(0, 0 + 3) = N2u_x1[gauss];
		(*deltaN[gauss])(1, 1 + 3) = N2u_x1[gauss];
		(*deltaN[gauss])(2, 2 + 3) = N2u_x1[gauss];
		(*deltaN[gauss])(0, 0 + 6) = N3u_x1[gauss];
		(*deltaN[gauss])(1, 1 + 6) = N3u_x1[gauss];
		(*deltaN[gauss])(2, 2 + 6) = N3u_x1[gauss];
		(*deltaN[gauss])(0, 0 + 9) = N4u_x1[gauss];
		(*deltaN[gauss])(1, 1 + 9) = N4u_x1[gauss];
		(*deltaN[gauss])(2, 2 + 9) = N4u_x1[gauss];
		(*deltaN[gauss])(0, 0 +12) = N5u_x1[gauss];
		(*deltaN[gauss])(1, 1 +12) = N5u_x1[gauss];
		(*deltaN[gauss])(2, 2 +12) = N5u_x1[gauss];
		(*deltaN[gauss])(0, 0 +15) = N6u_x1[gauss];
		(*deltaN[gauss])(1, 1 +15) = N6u_x1[gauss];
		(*deltaN[gauss])(2, 2 +15) = N6u_x1[gauss];
		(*deltaN[gauss])(3, 0 +18) = N4a_x1[gauss];
		(*deltaN[gauss])(4, 1 +18) = N4a_x1[gauss];
		(*deltaN[gauss])(5, 2 +18) = N4a_x1[gauss];
		(*deltaN[gauss])(3, 0 +21) = N5a_x1[gauss];
		(*deltaN[gauss])(4, 1 +21) = N5a_x1[gauss];
		(*deltaN[gauss])(5, 2 +21) = N5a_x1[gauss];
		(*deltaN[gauss])(3, 0 +24) = N6a_x1[gauss];
		(*deltaN[gauss])(4, 1 +24) = N6a_x1[gauss];
		(*deltaN[gauss])(5, 2 +24) = N6a_x1[gauss];
		(*deltaN[gauss])(6, 0 + 0) = N1u_x2[gauss];
		(*deltaN[gauss])(7, 1 + 0) = N1u_x2[gauss];
		(*deltaN[gauss])(8, 2 + 0) = N1u_x2[gauss];
		(*deltaN[gauss])(6, 0 + 3) = N2u_x2[gauss];
		(*deltaN[gauss])(7, 1 + 3) = N2u_x2[gauss];
		(*deltaN[gauss])(8, 2 + 3) = N2u_x2[gauss];
		(*deltaN[gauss])(6, 0 + 6) = N3u_x2[gauss];
		(*deltaN[gauss])(7, 1 + 6) = N3u_x2[gauss];
		(*deltaN[gauss])(8, 2 + 6) = N3u_x2[gauss];
		(*deltaN[gauss])(6, 0 + 9) = N4u_x2[gauss];
		(*deltaN[gauss])(7, 1 + 9) = N4u_x2[gauss];
		(*deltaN[gauss])(8, 2 + 9) = N4u_x2[gauss];
		(*deltaN[gauss])(6, 0 + 12) = N5u_x2[gauss];
		(*deltaN[gauss])(7, 1 + 12) = N5u_x2[gauss];
		(*deltaN[gauss])(8, 2 + 12) = N5u_x2[gauss];
		(*deltaN[gauss])(6, 0 + 15) = N6u_x2[gauss];
		(*deltaN[gauss])(7, 1 + 15) = N6u_x2[gauss];
		(*deltaN[gauss])(8, 2 + 15) = N6u_x2[gauss];
		(*deltaN[gauss])(9, 0 + 18) = N4a_x2[gauss];
		(*deltaN[gauss])(10, 1 + 18) = N4a_x2[gauss];
		(*deltaN[gauss])(11, 2 + 18) = N4a_x2[gauss];
		(*deltaN[gauss])(9, 0 + 21) = N5a_x2[gauss];
		(*deltaN[gauss])(10, 1 + 21) = N5a_x2[gauss];
		(*deltaN[gauss])(11, 2 + 21) = N5a_x2[gauss];
		(*deltaN[gauss])(9, 0 + 24) = N6a_x2[gauss];
		(*deltaN[gauss])(10, 1 + 24) = N6a_x2[gauss];
		(*deltaN[gauss])(11, 2 + 24) = N6a_x2[gauss];
		(*deltaN[gauss])(12, 0 + 18) = N4a[gauss];
		(*deltaN[gauss])(13, 1 + 18) = N4a[gauss];
		(*deltaN[gauss])(14, 2 + 18) = N4a[gauss];
		(*deltaN[gauss])(12, 0 + 21) = N5a[gauss];
		(*deltaN[gauss])(13, 1 + 21) = N5a[gauss];
		(*deltaN[gauss])(14, 2 + 21) = N5a[gauss];
		(*deltaN[gauss])(12, 0 + 24) = N6a[gauss];
		(*deltaN[gauss])(13, 1 + 24) = N6a[gauss];
		(*deltaN[gauss])(14, 2 + 24) = N6a[gauss];
	}

	//Percorre pontos de Gauss para salvar algumas matrizes/vetores de interesse nos cálculos do elemento
	double c1, c2, c3;
	for (int gauss = 0; gauss < 6; gauss++)
	{
		//Pontos de Gauss localizados com base no artigo de Cowper (1973)
		//Ponto 1
		if (gauss == 0)
		{
			c1 = 0.816847572980459;
			c2 = 0.091576213509771;
			c3 = 0.091576213509771;
			w4[gauss] = area * 0.109951743655322;
		}
		//Ponto 2
		if (gauss == 1)
		{
			c1 = 0.091576213509771;
			c2 = 0.816847572980459;
			c3 = 0.091576213509771;
			w4[gauss] = area * 0.109951743655322;
		}
		//Ponto 3
		if (gauss == 2)
		{
			c1 = 0.091576213509771;
			c2 = 0.091576213509771;
			c3 = 0.816847572980459;
			w4[gauss] = area * 0.109951743655322;
		}
		//Ponto 4
		if (gauss == 3)
		{
			c1 = 0.108103018168070;
			c2 = 0.445948490915965;
			c3 = 0.445948490915965;
			w4[gauss] = area * 0.223381589678011;
		}
		//Ponto 5
		if (gauss == 4)
		{
			c1 = 0.445948490915965;
			c2 = 0.108103018168070;
			c3 = 0.445948490915965;
			w4[gauss] = area * 0.223381589678011;
		}
		//Ponto 6
		if (gauss == 5)
		{
			c1 = 0.445948490915965;
			c2 = 0.445948490915965;
			c3 = 0.108103018168070;
			w4[gauss] = area * 0.223381589678011;
		}
		//xp = c1*x1 + c2*x2 + c3*x3;
		//Cálculo das funções de forma
		//A1 = 0.5 * norm(cross(x2 - xp, x3 - xp));
		//A2 = 0.5 * norm(cross(x3 - xp, x1 - xp));
		//A3 = 0.5 * norm(cross(x1 - xp, x2 - xp));
		//L1 = A1 / A;
		//L2 = A2 / A;
		//L3 = A3 / A;
		L1 = c1;
		L2 = c2;
		L3 = c3;

		double b1, b2, b3, c1, c2, c3;
		b1 = dot(x2 - x3, *e2r);
		b2 = dot(x3 - x1, *e2r);
		b3 = dot(x1 - x2, *e2r);
		c1 = dot(x3 - x2, *e1r);
		c2 = dot(x1 - x3, *e1r);
		c3 = dot(x2 - x1, *e1r);
		double L1_x1, L2_x1, L3_x1, L1_x2, L2_x2, L3_x2;
		L1_x1 = 0.5 * b1 / A;
		L2_x1 = 0.5 * b2 / A;
		L3_x1 = 0.5 * b3 / A;
		L1_x2 = 0.5 * c1 / A;
		L2_x2 = 0.5 * c2 / A;
		L3_x2 = 0.5 * c3 / A;
		//Atribuindo os valores das funções de forma e derivadas
		N1u4[gauss] = (2 * L1 - 1)*L1;
		N2u4[gauss] = (2 * L2 - 1)*L2;
		N3u4[gauss] = (2 * L3 - 1)*L3;
		N4u4[gauss] = 4 * L1 * L2;
		N5u4[gauss] = 4 * L2 * L3;
		N6u4[gauss] = 4 * L3 * L1;
		N4a4[gauss] = 1 - 2 * L3;
		N5a4[gauss] = 1 - 2 * L1;
		N6a4[gauss] = 1 - 2 * L2;

		N1u4_x1[gauss] = 4 * L1_x1*L1 - L1_x1;
		N2u4_x1[gauss] = 4 * L2_x1*L2 - L2_x1;
		N3u4_x1[gauss] = 4 * L3_x1*L3 - L3_x1;
		N4u4_x1[gauss] = 4 * L1_x1*L2 + 4 * L1*L2_x1;
		N5u4_x1[gauss] = 4 * L2_x1*L3 + 4 * L2*L3_x1;
		N6u4_x1[gauss] = 4 * L3_x1*L1 + 4 * L3*L1_x1;
		
		N1u4_x2[gauss] = 4 * L1_x2*L1 - L1_x2;
		N2u4_x2[gauss] = 4 * L2_x2*L2 - L2_x2;
		N3u4_x2[gauss] = 4 * L3_x2*L3 - L3_x2;
		N4u4_x2[gauss] = 4 * L1_x2*L2 + 4 * L1*L2_x2;
		N5u4_x2[gauss] = 4 * L2_x2*L3 + 4 * L2*L3_x2;
		N6u4_x2[gauss] = 4 * L3_x2*L1 + 4 * L3*L1_x2;
		
		//Preenchendo a matriz N4 (funções de forma):
		(*N4[gauss])(0, 0 + 0) = N1u4[gauss];
		(*N4[gauss])(1, 1 + 0) = N1u4[gauss];
		(*N4[gauss])(2, 2 + 0) = N1u4[gauss];
		(*N4[gauss])(0, 0 + 3) = N2u4[gauss];
		(*N4[gauss])(1, 1 + 3) = N2u4[gauss];
		(*N4[gauss])(2, 2 + 3) = N2u4[gauss];
		(*N4[gauss])(0, 0 + 6) = N3u4[gauss];
		(*N4[gauss])(1, 1 + 6) = N3u4[gauss];
		(*N4[gauss])(2, 2 + 6) = N3u4[gauss];
		(*N4[gauss])(0, 0 + 9) = N4u4[gauss];
		(*N4[gauss])(1, 1 + 9) = N4u4[gauss];
		(*N4[gauss])(2, 2 + 9) = N4u4[gauss];
		(*N4[gauss])(0, 0 + 12) = N5u4[gauss];
		(*N4[gauss])(1, 1 + 12) = N5u4[gauss];
		(*N4[gauss])(2, 2 + 12) = N5u4[gauss];
		(*N4[gauss])(0, 0 + 15) = N6u4[gauss];
		(*N4[gauss])(1, 1 + 15) = N6u4[gauss];
		(*N4[gauss])(2, 2 + 15) = N6u4[gauss];
		(*N4[gauss])(3, 0 + 18) = N4a4[gauss];
		(*N4[gauss])(4, 1 + 18) = N4a4[gauss];
		(*N4[gauss])(5, 2 + 18) = N4a4[gauss];
		(*N4[gauss])(3, 0 + 21) = N5a4[gauss];
		(*N4[gauss])(4, 1 + 21) = N5a4[gauss];
		(*N4[gauss])(5, 2 + 21) = N5a4[gauss];
		(*N4[gauss])(3, 0 + 24) = N6a4[gauss];
		(*N4[gauss])(4, 1 + 24) = N6a4[gauss];
		(*N4[gauss])(5, 2 + 24) = N6a4[gauss];

		//Preenchendo a matriz N4u_x1 (derivadas das funções de forma):
		(*N4_u_x1[gauss])(0, 0 + 0) = N1u4_x1[gauss];
		(*N4_u_x1[gauss])(1, 1 + 0) = N1u4_x1[gauss];
		(*N4_u_x1[gauss])(2, 2 + 0) = N1u4_x1[gauss];
		(*N4_u_x1[gauss])(0, 0 + 3) = N2u4_x1[gauss];
		(*N4_u_x1[gauss])(1, 1 + 3) = N2u4_x1[gauss];
		(*N4_u_x1[gauss])(2, 2 + 3) = N2u4_x1[gauss];
		(*N4_u_x1[gauss])(0, 0 + 6) = N3u4_x1[gauss];
		(*N4_u_x1[gauss])(1, 1 + 6) = N3u4_x1[gauss];
		(*N4_u_x1[gauss])(2, 2 + 6) = N3u4_x1[gauss];
		(*N4_u_x1[gauss])(0, 0 + 9) = N4u4_x1[gauss];
		(*N4_u_x1[gauss])(1, 1 + 9) = N4u4_x1[gauss];
		(*N4_u_x1[gauss])(2, 2 + 9) = N4u4_x1[gauss];
		(*N4_u_x1[gauss])(0, 0 + 12) = N5u4_x1[gauss];
		(*N4_u_x1[gauss])(1, 1 + 12) = N5u4_x1[gauss];
		(*N4_u_x1[gauss])(2, 2 + 12) = N5u4_x1[gauss];
		(*N4_u_x1[gauss])(0, 0 + 15) = N6u4_x1[gauss];
		(*N4_u_x1[gauss])(1, 1 + 15) = N6u4_x1[gauss];
		(*N4_u_x1[gauss])(2, 2 + 15) = N6u4_x1[gauss];
		//Preenchendo a matriz N4u_x2 (derivadas das funções de forma):
		(*N4_u_x2[gauss])(0, 0 + 0) = N1u4_x2[gauss];
		(*N4_u_x2[gauss])(1, 1 + 0) = N1u4_x2[gauss];
		(*N4_u_x2[gauss])(2, 2 + 0) = N1u4_x2[gauss];
		(*N4_u_x2[gauss])(0, 0 + 3) = N2u4_x2[gauss];
		(*N4_u_x2[gauss])(1, 1 + 3) = N2u4_x2[gauss];
		(*N4_u_x2[gauss])(2, 2 + 3) = N2u4_x2[gauss];
		(*N4_u_x2[gauss])(0, 0 + 6) = N3u4_x2[gauss];
		(*N4_u_x2[gauss])(1, 1 + 6) = N3u4_x2[gauss];
		(*N4_u_x2[gauss])(2, 2 + 6) = N3u4_x2[gauss];
		(*N4_u_x2[gauss])(0, 0 + 9) = N4u4_x2[gauss];
		(*N4_u_x2[gauss])(1, 1 + 9) = N4u4_x2[gauss];
		(*N4_u_x2[gauss])(2, 2 + 9) = N4u4_x2[gauss];
		(*N4_u_x2[gauss])(0, 0 + 12) = N5u4_x2[gauss];
		(*N4_u_x2[gauss])(1, 1 + 12) = N5u4_x2[gauss];
		(*N4_u_x2[gauss])(2, 2 + 12) = N5u4_x2[gauss];
		(*N4_u_x2[gauss])(0, 0 + 15) = N6u4_x2[gauss];
		(*N4_u_x2[gauss])(1, 1 + 15) = N6u4_x2[gauss];
		(*N4_u_x2[gauss])(2, 2 + 15) = N6u4_x2[gauss];
	}

	//Inicialização de variáveis acumuladas nos pontos de Gauss (Lagrangiano atualizado)
	for (int gauss = 0; gauss < 3; gauss++)
	{
		*z_x1_i[gauss] = *e1rlocal;
		*z_x2_i[gauss] = *e2rlocal;
		*Q_i[gauss] = I3;
	}

	alpha1 = area / 3.0; //peso dos pontos de integração com 3 pontos
}

//Monta a matriz de massa para análise modal
void Shell_1::MountMassModal()
{
	zeros(mass_modal);
	zeros(damping_modal);
	//Loop nos pontos de Gauss
	for (int gauss = 0; gauss < 6; gauss++)
	{
		//alpha_i4 calculado no ponto de Gauss							
		(*alpha_i4[gauss])(0, 0) =	(db.nodes[nodes[3] - 1]->copy_coordinates[3])* N4a4[gauss] +
									(db.nodes[nodes[4] - 1]->copy_coordinates[3])* N5a4[gauss] +
									(db.nodes[nodes[5] - 1]->copy_coordinates[3])* N6a4[gauss];
		(*alpha_i4[gauss])(1, 0) =	(db.nodes[nodes[3] - 1]->copy_coordinates[4])* N4a4[gauss] +
									(db.nodes[nodes[4] - 1]->copy_coordinates[4])* N5a4[gauss] +
									(db.nodes[nodes[5] - 1]->copy_coordinates[4])* N6a4[gauss];
		(*alpha_i4[gauss])(2, 0) =	(db.nodes[nodes[3] - 1]->copy_coordinates[5])* N4a4[gauss] +
									(db.nodes[nodes[4] - 1]->copy_coordinates[5])* N5a4[gauss] +
									(db.nodes[nodes[5] - 1]->copy_coordinates[5])* N6a4[gauss];
		*alpha_i4[gauss] = (*transform3)*(*alpha_i4[gauss]);					//TO LOCAL
		//Cálculo do integrando no ponto de Gauss
		EvaluateMassModal(temp_v, alpha_i4[gauss]->getMatrix(), &coef1, &coef2, &coef3, e3rlocal->getMatrix(), pDdT);
		//Transformando o operador tangente em Matrix
		DdT->PtrToMatrix(pDdT, 6);
		//DdT->print();
		(*mass_modal) = (*mass_modal) + w4[gauss]*transp(*N4[gauss])*(*DdT)*(*N4[gauss]);
	}//end of gauss
	//Conversão para o sistema global
	(*mass_modal) = (transp(*transform)*(*mass_modal))*(*transform);
}

//Monta a matriz de amortecimento para realização da análise modal
void Shell_1::MountDampingModal()
{
	zeros(mass_modal);
	zeros(damping_modal);
	//TODO
}

//Monta a matriz de massa
void Shell_1::MountMass()
{
	Matrix omega_i(3);
	Matrix domega_i(3);
	Matrix du_i(3);
	Matrix ddu_i(3);
	//Loop nos pontos de Gauss
	for (int gauss = 0; gauss < 3; gauss++)
	{
		//Calculando as variáveis necessárias para obter o valor do integrando							
		(omega_i)(0, 0) = (db.nodes[nodes[3] - 1]->copy_vel[3])* N4a[gauss] +
			(db.nodes[nodes[4] - 1]->copy_vel[3])* N5a[gauss] +
			(db.nodes[nodes[5] - 1]->copy_vel[3])* N6a[gauss];
		(omega_i)(1, 0) = (db.nodes[nodes[3] - 1]->copy_vel[4])* N4a[gauss] +
			(db.nodes[nodes[4] - 1]->copy_vel[4])* N5a[gauss] +
			(db.nodes[nodes[5] - 1]->copy_vel[4])* N6a[gauss];
		(omega_i)(2, 0) = (db.nodes[nodes[3] - 1]->copy_vel[5])* N4a[gauss] +
			(db.nodes[nodes[4] - 1]->copy_vel[5])* N5a[gauss] +
			(db.nodes[nodes[5] - 1]->copy_vel[5])* N6a[gauss];

		(domega_i)(0, 0) = (db.nodes[nodes[3] - 1]->copy_accel[3])* N4a[gauss] +
			(db.nodes[nodes[4] - 1]->copy_accel[3])* N5a[gauss] +
			(db.nodes[nodes[5] - 1]->copy_accel[3])* N6a[gauss];
		(domega_i)(1, 0) = (db.nodes[nodes[3] - 1]->copy_accel[4])* N4a[gauss] +
			(db.nodes[nodes[4] - 1]->copy_accel[4])* N5a[gauss] +
			(db.nodes[nodes[5] - 1]->copy_accel[4])* N6a[gauss];
		(domega_i)(2, 0) = (db.nodes[nodes[3] - 1]->copy_accel[5])* N4a[gauss] +
			(db.nodes[nodes[4] - 1]->copy_accel[5])* N5a[gauss] +
			(db.nodes[nodes[5] - 1]->copy_accel[5])* N6a[gauss];

		(du_i)(0, 0) = (db.nodes[nodes[0] - 1]->copy_vel[0])* N1u[gauss] +
			(db.nodes[nodes[1] - 1]->copy_vel[0])* N2u[gauss] +
			(db.nodes[nodes[2] - 1]->copy_vel[0])* N3u[gauss] +
			(db.nodes[nodes[3] - 1]->copy_vel[0])* N4u[gauss] +
			(db.nodes[nodes[4] - 1]->copy_vel[0])* N5u[gauss] +
			(db.nodes[nodes[5] - 1]->copy_vel[0])* N6u[gauss];
		(du_i)(1, 0) = (db.nodes[nodes[0] - 1]->copy_vel[1])* N1u[gauss] +
			(db.nodes[nodes[1] - 1]->copy_vel[1])* N2u[gauss] +
			(db.nodes[nodes[2] - 1]->copy_vel[1])* N3u[gauss] +
			(db.nodes[nodes[3] - 1]->copy_vel[1])* N4u[gauss] +
			(db.nodes[nodes[4] - 1]->copy_vel[1])* N5u[gauss] +
			(db.nodes[nodes[5] - 1]->copy_vel[1])* N6u[gauss];
		(du_i)(2, 0) = (db.nodes[nodes[0] - 1]->copy_vel[2])* N1u[gauss] +
			(db.nodes[nodes[1] - 1]->copy_vel[2])* N2u[gauss] +
			(db.nodes[nodes[2] - 1]->copy_vel[2])* N3u[gauss] +
			(db.nodes[nodes[3] - 1]->copy_vel[2])* N4u[gauss] +
			(db.nodes[nodes[4] - 1]->copy_vel[2])* N5u[gauss] +
			(db.nodes[nodes[5] - 1]->copy_vel[2])* N6u[gauss];

		(ddu_i)(0, 0) = (db.nodes[nodes[0] - 1]->copy_accel[0])* N1u[gauss] +
			(db.nodes[nodes[1] - 1]->copy_accel[0])* N2u[gauss] +
			(db.nodes[nodes[2] - 1]->copy_accel[0])* N3u[gauss] +
			(db.nodes[nodes[3] - 1]->copy_accel[0])* N4u[gauss] +
			(db.nodes[nodes[4] - 1]->copy_accel[0])* N5u[gauss] +
			(db.nodes[nodes[5] - 1]->copy_accel[0])* N6u[gauss];
		(ddu_i)(1, 0) = (db.nodes[nodes[0] - 1]->copy_accel[1])* N1u[gauss] +
			(db.nodes[nodes[1] - 1]->copy_accel[1])* N2u[gauss] +
			(db.nodes[nodes[2] - 1]->copy_accel[1])* N3u[gauss] +
			(db.nodes[nodes[3] - 1]->copy_accel[1])* N4u[gauss] +
			(db.nodes[nodes[4] - 1]->copy_accel[1])* N5u[gauss] +
			(db.nodes[nodes[5] - 1]->copy_accel[1])* N6u[gauss];
		(ddu_i)(2, 0) = (db.nodes[nodes[0] - 1]->copy_accel[2])* N1u[gauss] +
			(db.nodes[nodes[1] - 1]->copy_accel[2])* N2u[gauss] +
			(db.nodes[nodes[2] - 1]->copy_accel[2])* N3u[gauss] +
			(db.nodes[nodes[3] - 1]->copy_accel[2])* N4u[gauss] +
			(db.nodes[nodes[4] - 1]->copy_accel[2])* N5u[gauss] +
			(db.nodes[nodes[5] - 1]->copy_accel[2])* N6u[gauss];

		//Para o sistema do elemento
		omega_i = (*transform3)*omega_i;
		domega_i = (*transform3)*domega_i;
		du_i = (*transform3)*du_i;
		ddu_i = (*transform3)*ddu_i;

		if (typeid(*db.solution[db.current_solution_number - 1]) == typeid(Dynamic))
		{
			Dynamic* ptr_sol = static_cast<Dynamic*>(db.solution[db.current_solution_number - 1]);
			//Cálculo do integrando no ponto de Gauss
			EvaluateInertialContributions(temp_v, &ptr_sol->a1, &ptr_sol->a2, &ptr_sol->a3, &ptr_sol->a4,
				&ptr_sol->a5, &ptr_sol->a6, alpha_i[gauss]->getMatrix(), alpha_delta[gauss]->getMatrix(),
				u_i[gauss]->getMatrix(), u_delta[gauss]->getMatrix(), omega_i.getMatrix(), domega_i.getMatrix(), du_i.getMatrix(), ddu_i.getMatrix(), e3rlocal->getMatrix(), &coef1, &coef2, dT->getMatrix(), pDdT);
			//Transformando o operador tangente em Matrix
			DdT->PtrToMatrix(pDdT, 6);
			(*inertial_loading) = (*inertial_loading) + alpha1*transp(*N[gauss])*(*dT);
			(*mass) = (*mass) + alpha1*transp(*N[gauss])*(*DdT)*(*N[gauss]);
			kinetic_energy += (*tempkin) * (alpha1);
		}
		
	}//end of gauss
	//Conversão para o sistema global
	(*inertial_loading) = transp(*transform)*(*inertial_loading);
	(*mass) = (transp(*transform)*(*mass))*(*transform);
}

//Monta a matriz de amortecimento
void Shell_1::MountDamping(bool update_rayleigh)
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
		(*damping) = (*damping) + ptr_sol->a4*(*rayleigh_damping);											//Atualizando a matriz de amortecimento - inclusão do efeito de Rayleigh

		//Cálculo dos esforços de amortecimento - utilização de informação das velocidades nos GL do elemento
		for (int node = 0; node<6; node++)
		{
			//Percorre os 3 GL por nó
			for (int GL = 0; GL<3; GL++)
				v_ipp(node * 3 + GL, 0) = db.nodes[nodes[node] - 1]->vel[GL];
		}
		//Rotações no sistema local
		for (int node = 3; node<6; node++)
		{
			//Percorre os 3 GL por nó
			for (int GL = 3; GL<6; GL++)
			{
				v_ipp(18 + (node - 3) * 3 + GL - 3, 0) = db.nodes[nodes[node] - 1]->vel[GL];
			}
		}
		(*damping_loading) = (*damping_loading) + (*rayleigh_damping)*v_ipp;
	}
	
}

//Montagens - Newmark
void Shell_1::MountDyn()
{
	(*P_loading) = (*P_loading) + (*inertial_loading) + (*damping_loading);
	(*i_loading) = (*i_loading) + (*inertial_loading) + (*damping_loading);
	(*stiffness) = (*stiffness) + (*mass) + (*damping);
}

//Montagens para análise modal - inserção da matriz de massa e amortecimento na matriz de rigidez para posterior montagem global
void Shell_1::MountDynModal()
{
	(*stiffness) = (*mass_modal) + (*damping_modal);
}

//Zera matrizes
void Shell_1::Zeros()
{
	zeros(stiffness);
	zeros(mass);
	zeros(damping);
	zeros(damping_loading);
	zeros(i_loading);
	zeros(inertial_loading);
	zeros(e_loading);
	zeros(P_loading);
	zeros(mass_modal);
	zeros(damping_modal);
	kinetic_energy = 0.0;
	strain_energy = 0.0;
	potential_gravitational_energy = 0.0;
}

//Calcula as contribuições inerciais para a análise dinâmica - inclui todas as contribuições para a forma fraca e para o operador tangente
void Shell_1::EvaluateInertialContributions(double* v, double(*a1)
	, double(*a2), double(*a3), double(*a4), double(*a5), double(*a6)
	, double* alphai, double* alphad, double* ui, double* ud, double* omegai
	, double* domegai, double* dui, double* ddui, double* e3r, double(*coef1)
	, double(*coef2), double* dT, double** DdT)
{
	v[2] = (*a1);
	v[3] = (*a2);
	v[4] = (*a3);
	v[5] = (*a4);
	v[6] = (*a5);
	v[7] = (*a6);
	v[8] = alphai[0];
	v[112] = (v[8] * v[8]);
	v[9] = alphai[1];
	v[343] = 0.5e0 * v[9];
	v[110] = v[343] * v[8];
	v[105] = (v[9] * v[9]);
	v[10] = alphai[2];
	v[117] = v[10] * v[343];
	v[115] = 0.5e0 * v[10] * v[8];
	v[106] = (v[10] * v[10]);
	v[353] = v[105] + v[106];
	v[11] = alphad[0];
	v[177] = 0.5e0 * v[11];
	v[175] = 2e0 * v[11];
	v[96] = (v[11] * v[11]);
	v[12] = alphad[1];
	v[178] = 2e0 * v[12];
	v[176] = 0.5e0 * v[12];
	v[346] = 1e0 * v[176];
	v[94] = v[11] * v[346];
	v[89] = (v[12] * v[12]);
	v[223] = -v[89] - v[96];
	v[13] = alphad[2];
	v[201] = v[13] + v[94];
	v[191] = -v[13] + v[94];
	v[180] = 2e0 * v[13];
	v[179] = 0.5e0 * v[13];
	v[344] = 1e0 * v[179];
	v[101] = v[12] * v[344];
	v[218] = v[101] + v[11];
	v[210] = v[101] - v[11];
	v[99] = v[11] * v[344];
	v[214] = -v[12] + v[99];
	v[197] = v[12] + v[99];
	v[90] = (v[13] * v[13]);
	v[205] = -v[90] - v[96];
	v[186] = -v[89] - v[90];
	v[352] = 1e0 * v[186];
	v[181] = 4e0 + v[89] + v[90] + v[96];
	v[345] = -4e0 / (v[181] * v[181]);
	v[185] = v[180] * v[345];
	v[351] = 0.5e0 * v[185];
	v[227] = v[223] * v[351];
	v[184] = v[178] * v[345];
	v[348] = 0.5e0 * v[184];
	v[256] = -1e0 * v[184] * v[344];
	v[207] = v[205] * v[348];
	v[182] = v[175] * v[345];
	v[350] = 0.5e0 * v[182];
	v[347] = -1e0 * v[182];
	v[258] = v[182] * v[346];
	v[255] = v[344] * v[347];
	v[187] = v[186] * v[350];
	v[17] = ud[0];
	v[18] = ud[1];
	v[19] = ud[2];
	v[20] = omegai[0];
	v[21] = omegai[1];
	v[22] = omegai[2];
	v[23] = domegai[0];
	v[148] = v[11] * v[2] - v[20] * v[3] - v[23] * v[4];
	v[142] = v[11] * v[5] + v[20] * v[6] + v[23] * v[7];
	v[24] = domegai[1];
	v[149] = v[12] * v[2] - v[21] * v[3] - v[24] * v[4];
	v[143] = v[12] * v[5] + v[21] * v[6] + v[24] * v[7];
	v[25] = domegai[2];
	v[150] = v[13] * v[2] - v[22] * v[3] - v[25] * v[4];
	v[144] = v[13] * v[5] + v[22] * v[6] + v[25] * v[7];
	v[26] = dui[0];
	v[27] = dui[1];
	v[28] = dui[2];
	v[29] = ddui[0];
	v[30] = ddui[1];
	v[31] = ddui[2];
	v[32] = e3r[0];
	v[33] = e3r[1];
	v[34] = e3r[2];
	v[35] = (*coef1);
	v[174] = v[2] * v[35];
	v[36] = (*coef2);
	v[88] = 4e0 / v[181];
	v[259] = -0.5e0 * v[88];
	v[349] = 1e0 * v[259];
	v[261] = v[259] + v[177] * v[347];
	v[260] = -v[259] + 1e0 * v[184] * v[346];
	v[257] = v[259] - 1e0 * v[185] * v[344];
	v[225] = v[178] * v[349];
	v[226] = v[225] + v[223] * v[348];
	v[222] = v[175] * v[349];
	v[224] = v[222] + v[223] * v[350];
	v[219] = v[182] * v[218] + v[88];
	v[216] = v[184] * v[214] - v[88];
	v[211] = v[182] * v[210] - v[88];
	v[208] = 1e0 * v[180] * v[349];
	v[209] = v[208] + v[205] * v[351];
	v[206] = v[222] + 1e0 * v[205] * v[350];
	v[204] = v[185] * v[201] + v[88];
	v[199] = v[184] * v[197] + v[88];
	v[196] = -(v[179] * v[88]);
	v[220] = -v[196] + v[184] * v[218];
	v[215] = -v[196] + v[182] * v[214];
	v[212] = -v[196] + v[184] * v[210];
	v[198] = -v[196] + v[182] * v[197];
	v[195] = v[185] * v[191] - v[88];
	v[193] = -(v[177] * v[88]);
	v[217] = -v[193] + v[185] * v[214];
	v[203] = -v[193] + v[184] * v[201];
	v[200] = -v[193] + v[185] * v[197];
	v[194] = v[184] * v[191] - v[193];
	v[190] = v[176] * v[88];
	v[221] = v[190] + v[185] * v[218];
	v[213] = v[190] + v[185] * v[210];
	v[202] = v[190] + v[182] * v[201];
	v[192] = v[190] + v[182] * v[191];
	v[189] = v[208] + v[351] * v[352];
	v[188] = v[225] + v[348] * v[352];
	v[91] = 1e0 - 1e0 * v[259] * v[352];
	v[289] = v[148] * v[187] + v[149] * v[192] + v[150] * v[198] + v[2] * v[91];
	v[271] = v[142] * v[187] + v[143] * v[192] + v[144] * v[198] + v[5] * v[91];
	v[92] = v[191] * v[88];
	v[290] = v[148] * v[188] + v[149] * v[194] + v[150] * v[199] + v[2] * v[92];
	v[272] = v[142] * v[188] + v[143] * v[194] + v[144] * v[199] + v[5] * v[92];
	v[93] = v[197] * v[88];
	v[291] = v[148] * v[189] + v[149] * v[195] + v[150] * v[200] + v[2] * v[93];
	v[273] = v[142] * v[189] + v[143] * v[195] + v[144] * v[200] + v[5] * v[93];
	v[95] = v[201] * v[88];
	v[292] = v[148] * v[202] + v[149] * v[206] + v[150] * v[211] + v[2] * v[95];
	v[274] = v[142] * v[202] + v[143] * v[206] + v[144] * v[211] + v[5] * v[95];
	v[97] = 1e0 - 1e0 * v[205] * v[349];
	v[293] = v[148] * v[203] + v[149] * v[207] + v[150] * v[212] + v[2] * v[97];
	v[275] = v[142] * v[203] + v[143] * v[207] + v[144] * v[212] + v[5] * v[97];
	v[98] = v[210] * v[88];
	v[294] = v[148] * v[204] + v[149] * v[209] + v[150] * v[213] + v[2] * v[98];
	v[276] = v[142] * v[204] + v[143] * v[209] + v[144] * v[213] + v[5] * v[98];
	v[100] = v[214] * v[88];
	v[295] = v[100] * v[2] + v[148] * v[215] + v[149] * v[219] + v[150] * v[224];
	v[280] = v[142] * v[215] + v[143] * v[219] + v[144] * v[224] + v[100] * v[5];
	v[102] = v[218] * v[88];
	v[296] = v[102] * v[2] + v[148] * v[216] + v[149] * v[220] + v[150] * v[226];
	v[281] = v[142] * v[216] + v[143] * v[220] + v[144] * v[226] + v[102] * v[5];
	v[103] = 1e0 - 1e0 * v[223] * v[349];
	v[297] = v[103] * v[2] + v[148] * v[217] + v[149] * v[221] + v[150] * v[227];
	v[282] = v[142] * v[217] + v[143] * v[221] + v[144] * v[227] + v[103] * v[5];
	v[104] = 4e0 / (4e0 + v[112] + v[353]);
	v[354] = -0.5e0 * v[104];
	v[107] = 1e0 + v[353] * v[354];
	v[108] = v[104] * (-v[10] + v[110]);
	v[109] = v[104] * (v[115] + v[9]);
	v[111] = v[104] * (v[10] + v[110]);
	v[113] = 1e0 + (v[106] + v[112]) * v[354];
	v[114] = v[104] * (v[117] - v[8]);
	v[116] = v[104] * (v[115] - v[9]);
	v[118] = v[104] * (v[117] + v[8]);
	v[119] = 1e0 + 1e0 * (v[105] + v[112]) * v[354];
	v[270] = (v[107] * v[217] + v[111] * v[221] + v[116] * v[227]) * v[32] + (v[108] * v[217] + v[113] * v[221] + v[118] * v[227]
		) * v[33] + (v[109] * v[217] + v[114] * v[221] + v[119] * v[227]) * v[34];
	v[269] = (v[107] * v[216] + v[111] * v[220] + v[116] * v[226]) * v[32] + (v[108] * v[216] + v[113] * v[220] + v[118] * v[226]
		) * v[33] + (v[109] * v[216] + v[114] * v[220] + v[119] * v[226]) * v[34];
	v[268] = (v[107] * v[215] + v[111] * v[219] + v[116] * v[224]) * v[32] + (v[108] * v[215] + v[113] * v[219] + v[118] * v[224]
		) * v[33] + (v[109] * v[215] + v[114] * v[219] + v[119] * v[224]) * v[34];
	v[267] = (v[107] * v[204] + v[111] * v[209] + v[116] * v[213]) * v[32] + (v[108] * v[204] + v[113] * v[209] + v[118] * v[213]
		) * v[33] + (v[109] * v[204] + v[114] * v[209] + v[119] * v[213]) * v[34];
	v[266] = (v[107] * v[203] + v[111] * v[207] + v[116] * v[212]) * v[32] + (v[108] * v[203] + v[113] * v[207] + v[118] * v[212]
		) * v[33] + (v[109] * v[203] + v[114] * v[207] + v[119] * v[212]) * v[34];
	v[265] = (v[107] * v[202] + v[111] * v[206] + v[116] * v[211]) * v[32] + (v[108] * v[202] + v[113] * v[206] + v[118] * v[211]
		) * v[33] + (v[109] * v[202] + v[114] * v[206] + v[119] * v[211]) * v[34];
	v[264] = (v[107] * v[189] + v[111] * v[195] + v[116] * v[200]) * v[32] + (v[108] * v[189] + v[113] * v[195] + v[118] * v[200]
		) * v[33] + (v[109] * v[189] + v[114] * v[195] + v[119] * v[200]) * v[34];
	v[263] = (v[107] * v[188] + v[111] * v[194] + v[116] * v[199]) * v[32] + (v[108] * v[188] + v[113] * v[194] + v[118] * v[199]
		) * v[33] + (v[109] * v[188] + v[114] * v[194] + v[119] * v[199]) * v[34];
	v[262] = (v[107] * v[187] + v[111] * v[192] + v[116] * v[198]) * v[32] + (v[108] * v[187] + v[113] * v[192] + v[118] * v[198]
		) * v[33] + (v[109] * v[187] + v[114] * v[192] + v[119] * v[198]) * v[34];
	v[132] = v[32] * (v[107] * v[91] + v[111] * v[92] + v[116] * v[93]) + v[33] * (v[108] * v[91] + v[113] * v[92] + v[118] * v[93]
		) + v[34] * (v[109] * v[91] + v[114] * v[92] + v[119] * v[93]);
	v[375] = v[132] * v[289];
	v[372] = v[132] * v[290];
	v[369] = v[132] * v[291];
	v[133] = v[32] * (v[107] * v[95] + v[111] * v[97] + v[116] * v[98]) + v[33] * (v[108] * v[95] + v[113] * v[97] + v[118] * v[98]
		) + v[34] * (v[109] * v[95] + v[114] * v[97] + v[119] * v[98]);
	v[367] = v[133] * v[292];
	v[364] = v[133] * v[293];
	v[360] = v[133] * v[294];
	v[134] = (v[100] * v[107] + v[102] * v[111] + v[103] * v[116]) * v[32] + (v[100] * v[108] + v[102] * v[113] + v[103] * v[118]
		) * v[33] + (v[100] * v[109] + v[102] * v[114] + v[103] * v[119]) * v[34];
	v[376] = v[134] * v[295];
	v[373] = v[134] * v[296];
	v[370] = v[134] * v[297];
	v[141] = v[142] * v[91] + v[143] * v[92] + v[144] * v[93];
	v[145] = v[142] * v[95] + v[143] * v[97] + v[144] * v[98];
	v[279] = v[145] * v[264] - v[141] * v[267] - v[133] * v[273] + v[132] * v[276];
	v[278] = v[145] * v[263] - v[141] * v[266] - v[133] * v[272] + v[132] * v[275];
	v[277] = v[145] * v[262] - v[141] * v[265] - v[133] * v[271] + v[132] * v[274];
	v[157] = -(v[133] * v[141]) + v[132] * v[145];
	v[146] = v[100] * v[142] + v[102] * v[143] + v[103] * v[144];
	v[288] = -(v[146] * v[264]) + v[141] * v[270] + v[134] * v[273] - v[132] * v[282];
	v[287] = -(v[146] * v[263]) + v[141] * v[269] + v[134] * v[272] - v[132] * v[281];
	v[286] = -(v[146] * v[262]) + v[141] * v[268] + v[134] * v[271] - v[132] * v[280];
	v[285] = v[146] * v[267] - v[145] * v[270] - v[134] * v[276] + v[133] * v[282];
	v[284] = v[146] * v[266] - v[145] * v[269] - v[134] * v[275] + v[133] * v[281];
	v[283] = v[146] * v[265] - v[145] * v[268] - v[134] * v[274] + v[133] * v[280];
	v[162] = -(v[134] * v[145]) + v[133] * v[146];
	v[303] = -(v[157] * v[273]) - v[141] * v[279] + v[162] * v[282] + v[146] * v[285];
	v[302] = -(v[157] * v[272]) - v[141] * v[278] + v[162] * v[281] + v[146] * v[284];
	v[301] = -(v[157] * v[271]) - v[141] * v[277] + v[162] * v[280] + v[146] * v[283];
	v[161] = v[134] * v[141] - v[132] * v[146];
	v[316] = v[161] * v[273] - v[162] * v[276] - v[145] * v[285] + v[141] * v[288];
	v[315] = v[161] * v[272] - v[162] * v[275] - v[145] * v[284] + v[141] * v[287];
	v[314] = v[161] * v[271] - v[162] * v[274] - v[145] * v[283] + v[141] * v[286];
	v[300] = v[157] * v[276] + v[145] * v[279] - v[161] * v[282] - v[146] * v[288];
	v[299] = v[157] * v[275] + v[145] * v[278] - v[161] * v[281] - v[146] * v[287];
	v[298] = v[157] * v[274] + v[145] * v[277] - v[161] * v[280] - v[146] * v[286];
	v[147] = v[148] * v[91] + v[149] * v[92] + v[150] * v[93];
	v[365] = v[147] * v[262];
	v[362] = v[147] * v[263];
	v[358] = v[147] * v[264];
	v[355] = -(v[132] * v[147]);
	v[151] = v[148] * v[95] + v[149] * v[97] + v[150] * v[98];
	v[366] = v[151] * v[265];
	v[363] = v[151] * v[266];
	v[359] = v[151] * v[267];
	v[356] = -(v[133] * v[151]);
	v[152] = v[100] * v[148] + v[102] * v[149] + v[103] * v[150];
	v[374] = v[152] * v[268];
	v[371] = v[152] * v[269];
	v[368] = v[152] * v[270];
	v[361] = -(v[134] * v[152]);
	v[357] = 2e0 * v[152];
	v[156] = v[145] * v[157] - v[146] * v[161];
	v[158] = -(v[141] * v[157]) + v[146] * v[162];
	v[159] = (v[132] * v[132]);
	v[160] = (v[133] * v[133]);
	v[310] = v[159] + v[160];
	v[313] = -(v[158] * v[264]) + v[156] * v[267] + v[133] * v[300] - v[132] * v[303] + v[297] * v[310] + v[270] * v[355]
		+ v[270] * v[356] - v[134] * (v[270] * v[357] + v[358] + v[359] + v[360] + v[369]);
	v[312] = -(v[158] * v[263]) + v[156] * v[266] + v[133] * v[299] - v[132] * v[302] + v[296] * v[310] + v[269] * v[355]
		+ v[269] * v[356] - v[134] * (v[269] * v[357] + v[362] + v[363] + v[364] + v[372]);
	v[311] = -(v[158] * v[262]) + v[156] * v[265] + v[133] * v[298] - v[132] * v[301] + v[295] * v[310] + v[268] * v[355]
		+ v[268] * v[356] - v[134] * (v[268] * v[357] + v[365] + v[366] + v[367] + v[375]);
	v[166] = v[133] * (-(v[134] * v[151]) + v[156]) - v[132] * (v[134] * v[147] + v[158]) + v[152] * v[310];
	v[332] = -(v[166] * v[258]);
	v[163] = v[141] * v[161] - v[145] * v[162];
	v[164] = (v[134] * v[134]);
	v[324] = v[160] + v[164];
	v[327] = -(v[163] * v[267]) + v[158] * v[270] + v[134] * v[303] - v[133] * (v[151] * v[264] + v[316]) + v[291] * v[324]
		+ v[264] * v[361] - v[132] * (2e0 * v[358] + v[359] + v[360] + v[368] + v[370]);
	v[326] = -(v[163] * v[266]) + v[158] * v[269] + v[134] * v[302] - v[133] * (v[151] * v[263] + v[315]) + v[290] * v[324]
		+ v[263] * v[361] - v[132] * (2e0 * v[362] + v[363] + v[364] + v[371] + v[373]);
	v[325] = -(v[163] * v[265]) + v[158] * v[268] + v[134] * v[301] - v[133] * (v[151] * v[262] + v[314]) + v[289] * v[324]
		+ v[262] * v[361] - v[132] * (2e0 * v[365] + v[366] + v[367] + v[374] + v[376]);
	v[320] = v[159] + v[164];
	v[323] = v[163] * v[264] - v[156] * v[270] - v[134] * v[300] + v[132] * v[316] + v[294] * v[320] + v[267] * v[355]
		+ v[267] * v[361] - v[133] * (v[358] + 2e0 * v[359] + v[368] + v[369] + v[370]);
	v[322] = v[163] * v[263] - v[156] * v[269] - v[134] * v[299] + v[132] * v[315] + v[293] * v[320] + v[266] * v[355]
		+ v[266] * v[361] - v[133] * (v[362] + 2e0 * v[363] + v[371] + v[372] + v[373]);
	v[321] = v[163] * v[262] - v[156] * v[268] - v[134] * v[298] + v[132] * v[314] + v[292] * v[320] + v[265] * v[355]
		+ v[265] * v[361] - v[133] * (v[365] + 2e0 * v[366] + v[374] + v[375] + v[376]);
	v[168] = v[134] * (-(v[133] * v[152]) - v[156]) + v[132] * (-(v[133] * v[147]) + v[163]) + v[151] * v[320];
	v[337] = -(v[168] * v[255]);
	v[167] = v[134] * (-(v[132] * v[152]) + v[158]) + v[133] * (-(v[132] * v[151]) - v[163]) + v[147] * v[324];
	v[338] = v[167] * v[256];
	dT[0] = v[35] * (v[17] * v[2] - v[26] * v[3] - v[29] * v[4]);
	dT[1] = v[35] * (v[18] * v[2] - v[27] * v[3] - v[30] * v[4]);
	dT[2] = v[35] * (v[19] * v[2] - v[28] * v[3] - v[31] * v[4]);
	dT[3] = v[36] * (-(v[166] * v[190]) - v[168] * v[196] + v[167] * v[88]);
	dT[4] = v[36] * (-(v[166] * v[193]) + v[167] * v[196] + v[168] * v[88]);
	dT[5] = v[36] * (v[167] * v[190] + v[168] * v[193] + v[166] * v[88]);
	DdT[0][0] = v[174];
	DdT[0][1] = 0e0;
	DdT[0][2] = 0e0;
	DdT[0][3] = 0e0;
	DdT[0][4] = 0e0;
	DdT[0][5] = 0e0;
	DdT[1][0] = 0e0;
	DdT[1][1] = v[174];
	DdT[1][2] = 0e0;
	DdT[1][3] = 0e0;
	DdT[1][4] = 0e0;
	DdT[1][5] = 0e0;
	DdT[2][0] = 0e0;
	DdT[2][1] = 0e0;
	DdT[2][2] = v[174];
	DdT[2][3] = 0e0;
	DdT[2][4] = 0e0;
	DdT[2][5] = 0e0;
	DdT[3][0] = 0e0;
	DdT[3][1] = 0e0;
	DdT[3][2] = 0e0;
	DdT[3][3] = v[36] * (v[167] * v[182] - v[190] * v[311] - v[196] * v[321] + v[332] + v[337] + v[325] * v[88]);
	DdT[3][4] = v[36] * (v[167] * v[184] - v[168] * v[256] - v[166] * v[260] - v[190] * v[312] - v[196] * v[322] + v[326] * v[88]
		);
	DdT[3][5] = v[36] * (v[167] * v[185] + v[166] * v[256] - v[168] * v[257] - v[190] * v[313] - v[196] * v[323] + v[327] * v[88]
		);
	DdT[4][0] = 0e0;
	DdT[4][1] = 0e0;
	DdT[4][2] = 0e0;
	DdT[4][3] = v[36] * (v[168] * v[182] + v[167] * v[255] - v[166] * v[261] - v[193] * v[311] + v[196] * v[325] + v[321] * v[88]
		);
	DdT[4][4] = v[36] * (v[168] * v[184] - v[193] * v[312] + v[196] * v[326] - v[332] + v[338] + v[322] * v[88]);
	DdT[4][5] = v[36] * (v[168] * v[185] - v[166] * v[255] + v[167] * v[257] - v[193] * v[313] + v[196] * v[327] + v[323] * v[88]
		);
	DdT[5][0] = 0e0;
	DdT[5][1] = 0e0;
	DdT[5][2] = 0e0;
	DdT[5][3] = v[36] * (v[166] * v[182] + v[167] * v[258] + v[168] * v[261] + v[193] * v[321] + v[190] * v[325] + v[311] * v[88]
		);
	DdT[5][4] = v[36] * (v[166] * v[184] - v[168] * v[258] + v[167] * v[260] + v[193] * v[322] + v[190] * v[326] + v[312] * v[88]
		);
	DdT[5][5] = v[36] * (v[166] * v[185] + v[193] * v[323] + v[190] * v[327] - v[337] - v[338] + v[313] * v[88]);
	(*tempkin) = 0.5e0 * (((v[157] * v[157]) + (v[161] * v[161]) + (v[162] * v[162])) * v[36] + v[35] * (Power(v[17] * v[5]
		+ v[26] * v[6] + v[29] * v[7], 2) + Power(v[18] * v[5] + v[27] * v[6] + v[30] * v[7], 2) + Power(v[19] * v[5] + v[28] * v[6]
			+ v[31] * v[7], 2)));
	/*
	v[388] = (*a1)*(*coef1);
	v[360] = (*a4)*alphad[2] + (*a6)*domegai[2] + (*a5)*omegai[2];
	v[359] = (*a1)*alphad[2] - (*a3)*domegai[2] - (*a2)*omegai[2];
	v[358] = (*a4)*alphad[1] + (*a6)*domegai[1] + (*a5)*omegai[1];
	v[357] = (*a1)*alphad[1] - (*a3)*domegai[1] - (*a2)*omegai[1];
	v[356] = (*a4)*alphad[0] + (*a6)*domegai[0] + (*a5)*omegai[0];
	v[355] = (*a1)*alphad[0] - (*a3)*domegai[0] - (*a2)*omegai[0];
	v[354] = 0.5e0*alphad[2];
	v[352] = Power(alphad[2], 2);
	v[351] = alphad[0] * v[354];
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
	v[364] = v[340] + v[343];
	v[339] = alphai[0] * v[341];
	v[338] = Power(alphai[0], 2);
	v[116] = alphai[2] * v[341];
	v[219] = -v[345] - v[348];
	v[197] = alphad[2] + v[347];
	v[187] = -alphad[2] + v[347];
	v[100] = alphad[2] * v[350];
	v[214] = alphad[0] + v[100];
	v[206] = -alphad[0] + v[100];
	v[210] = -alphad[1] + v[351];
	v[193] = alphad[1] + v[351];
	v[201] = -v[345] - v[352];
	v[182] = -v[348] - v[352];
	v[177] = 4e0 + v[345] + v[348] + v[352];
	v[179] = 1e0 / Power(v[177], 2);
	v[353] = -4e0*v[179];
	v[181] = v[349] * v[353];
	v[363] = 0.5e0*v[181];
	v[223] = v[219] * v[363];
	v[180] = v[346] * v[353];
	v[362] = 0.5e0*v[180];
	v[225] = -(v[180] * v[354]);
	v[203] = v[201] * v[362];
	v[178] = v[344] * v[353];
	v[361] = 0.5e0*v[178];
	v[227] = v[178] * v[350];
	v[224] = -(v[178] * v[354]);
	v[183] = v[182] * v[361];
	v[87] = 4e0 / v[177];
	v[228] = -0.5e0*v[87];
	v[230] = v[228] - alphad[0] * v[361];
	v[229] = -v[228] + alphad[1] * v[362];
	v[226] = v[228] - alphad[2] * v[363];
	v[221] = v[228] * v[346];
	v[222] = v[221] + v[219] * v[362];
	v[218] = v[228] * v[344];
	v[220] = v[218] + v[219] * v[361];
	v[215] = v[178] * v[214] + v[87];
	v[212] = v[180] * v[210] - v[87];
	v[207] = v[178] * v[206] - v[87];
	v[204] = v[228] * v[349];
	v[205] = v[204] + v[201] * v[363];
	v[202] = v[218] + v[201] * v[361];
	v[200] = v[181] * v[197] + v[87];
	v[195] = v[180] * v[193] + v[87];
	v[192] = alphad[2] * v[228];
	v[216] = -v[192] + v[180] * v[214];
	v[211] = -v[192] + v[178] * v[210];
	v[208] = -v[192] + v[180] * v[206];
	v[194] = -v[192] + v[178] * v[193];
	v[191] = v[181] * v[187] - v[87];
	v[189] = alphad[0] * v[228];
	v[213] = -v[189] + v[181] * v[210];
	v[199] = -v[189] + v[180] * v[197];
	v[196] = -v[189] + v[181] * v[193];
	v[190] = v[180] * v[187] - v[189];
	v[186] = -(alphad[1] * v[228]);
	v[217] = v[186] + v[181] * v[214];
	v[209] = v[186] + v[181] * v[206];
	v[198] = v[186] + v[178] * v[197];
	v[188] = v[186] + v[178] * v[187];
	v[185] = v[204] + v[182] * v[363];
	v[184] = v[221] + v[182] * v[362];
	v[90] = 1e0 - v[182] * v[228];
	v[285] = v[183] * v[355] + v[188] * v[357] + v[194] * v[359] + (*a1)*v[90];
	v[267] = v[183] * v[356] + v[188] * v[358] + v[194] * v[360] + (*a4)*v[90];
	v[91] = v[187] * v[87];
	v[286] = v[184] * v[355] + v[190] * v[357] + v[195] * v[359] + (*a1)*v[91];
	v[268] = v[184] * v[356] + v[190] * v[358] + v[195] * v[360] + (*a4)*v[91];
	v[92] = v[193] * v[87];
	v[287] = v[185] * v[355] + v[191] * v[357] + v[196] * v[359] + (*a1)*v[92];
	v[269] = v[185] * v[356] + v[191] * v[358] + v[196] * v[360] + (*a4)*v[92];
	v[94] = v[197] * v[87];
	v[288] = v[198] * v[355] + v[202] * v[357] + v[207] * v[359] + (*a1)*v[94];
	v[270] = v[198] * v[356] + v[202] * v[358] + v[207] * v[360] + (*a4)*v[94];
	v[96] = 1e0 - v[201] * v[228];
	v[289] = v[199] * v[355] + v[203] * v[357] + v[208] * v[359] + (*a1)*v[96];
	v[271] = v[199] * v[356] + v[203] * v[358] + v[208] * v[360] + (*a4)*v[96];
	v[97] = v[206] * v[87];
	v[290] = v[200] * v[355] + v[205] * v[357] + v[209] * v[359] + (*a1)*v[97];
	v[272] = v[200] * v[356] + v[205] * v[358] + v[209] * v[360] + (*a4)*v[97];
	v[99] = v[210] * v[87];
	v[291] = v[211] * v[355] + v[215] * v[357] + v[220] * v[359] + (*a1)*v[99];
	v[276] = v[211] * v[356] + v[215] * v[358] + v[220] * v[360] + (*a4)*v[99];
	v[101] = v[214] * v[87];
	v[292] = (*a1)*v[101] + v[212] * v[355] + v[216] * v[357] + v[222] * v[359];
	v[277] = (*a4)*v[101] + v[212] * v[356] + v[216] * v[358] + v[222] * v[360];
	v[102] = 1e0 - v[219] * v[228];
	v[293] = (*a1)*v[102] + v[213] * v[355] + v[217] * v[357] + v[223] * v[359];
	v[278] = (*a4)*v[102] + v[213] * v[356] + v[217] * v[358] + v[223] * v[360];
	v[103] = 4e0 / (4e0 + v[338] + v[364]);
	v[365] = -0.5e0*v[103];
	v[106] = 1e0 + v[364] * v[365];
	v[107] = v[103] * (-alphai[2] + v[339]);
	v[108] = v[103] * (alphai[1] + v[342]);
	v[110] = v[103] * (alphai[2] + v[339]);
	v[112] = 1e0 + (v[338] + v[343])*v[365];
	v[113] = v[103] * (-alphai[0] + v[116]);
	v[115] = v[103] * (-alphai[1] + v[342]);
	v[117] = v[103] * (alphai[0] + v[116]);
	v[118] = 1e0 + (v[338] + v[340])*v[365];
	v[266] = e3r[0] * (v[106] * v[213] + v[110] * v[217] + v[115] * v[223]) + e3r[1] * (v[107] * v[213] + v[112] * v[217]
		+ v[117] * v[223]) + e3r[2] * (v[108] * v[213] + v[113] * v[217] + v[118] * v[223]);
	v[262] = e3r[0] * (v[106] * v[212] + v[110] * v[216] + v[115] * v[222]) + e3r[1] * (v[107] * v[212] + v[112] * v[216]
		+ v[117] * v[222]) + e3r[2] * (v[108] * v[212] + v[113] * v[216] + v[118] * v[222]);
	v[258] = e3r[0] * (v[106] * v[211] + v[110] * v[215] + v[115] * v[220]) + e3r[1] * (v[107] * v[211] + v[112] * v[215]
		+ v[117] * v[220]) + e3r[2] * (v[108] * v[211] + v[113] * v[215] + v[118] * v[220]);
	v[254] = e3r[0] * (v[106] * v[200] + v[110] * v[205] + v[115] * v[209]) + e3r[1] * (v[107] * v[200] + v[112] * v[205]
		+ v[117] * v[209]) + e3r[2] * (v[108] * v[200] + v[113] * v[205] + v[118] * v[209]);
	v[250] = e3r[0] * (v[106] * v[199] + v[110] * v[203] + v[115] * v[208]) + e3r[1] * (v[107] * v[199] + v[112] * v[203]
		+ v[117] * v[208]) + e3r[2] * (v[108] * v[199] + v[113] * v[203] + v[118] * v[208]);
	v[246] = e3r[0] * (v[106] * v[198] + v[110] * v[202] + v[115] * v[207]) + e3r[1] * (v[107] * v[198] + v[112] * v[202]
		+ v[117] * v[207]) + e3r[2] * (v[108] * v[198] + v[113] * v[202] + v[118] * v[207]);
	v[242] = e3r[0] * (v[106] * v[185] + v[110] * v[191] + v[115] * v[196]) + e3r[1] * (v[107] * v[185] + v[112] * v[191]
		+ v[117] * v[196]) + e3r[2] * (v[108] * v[185] + v[113] * v[191] + v[118] * v[196]);
	v[238] = e3r[0] * (v[106] * v[184] + v[110] * v[190] + v[115] * v[195]) + e3r[1] * (v[107] * v[184] + v[112] * v[190]
		+ v[117] * v[195]) + e3r[2] * (v[108] * v[184] + v[113] * v[190] + v[118] * v[195]);
	v[234] = e3r[0] * (v[106] * v[183] + v[110] * v[188] + v[115] * v[194]) + e3r[1] * (v[107] * v[183] + v[112] * v[188]
		+ v[117] * v[194]) + e3r[2] * (v[108] * v[183] + v[113] * v[188] + v[118] * v[194]);
	v[131] = e3r[0] * (v[106] * v[90] + v[110] * v[91] + v[115] * v[92]) + e3r[1] * (v[107] * v[90] + v[112] * v[91]
		+ v[117] * v[92]) + e3r[2] * (v[108] * v[90] + v[113] * v[91] + v[118] * v[92]);
	v[386] = v[131] * v[285];
	v[383] = v[131] * v[286];
	v[380] = v[131] * v[287];
	v[132] = e3r[0] * (v[106] * v[94] + v[110] * v[96] + v[115] * v[97]) + e3r[1] * (v[107] * v[94] + v[112] * v[96]
		+ v[117] * v[97]) + e3r[2] * (v[108] * v[94] + v[113] * v[96] + v[118] * v[97]);
	v[378] = v[132] * v[288];
	v[375] = v[132] * v[289];
	v[371] = v[132] * v[290];
	v[133] = e3r[0] * (v[101] * v[110] + v[102] * v[115] + v[106] * v[99]) + e3r[1] * (v[101] * v[112] + v[102] * v[117]
		+ v[107] * v[99]) + e3r[2] * (v[101] * v[113] + v[102] * v[118] + v[108] * v[99]);
	v[387] = v[133] * v[291];
	v[384] = v[133] * v[292];
	v[381] = v[133] * v[293];
	v[140] = v[356] * v[90] + v[358] * v[91] + v[360] * v[92];
	v[144] = v[356] * v[94] + v[358] * v[96] + v[360] * v[97];
	v[275] = v[144] * v[242] - v[140] * v[254] - v[132] * v[269] + v[131] * v[272];
	v[274] = v[144] * v[238] - v[140] * v[250] - v[132] * v[268] + v[131] * v[271];
	v[273] = v[144] * v[234] - v[140] * v[246] - v[132] * v[267] + v[131] * v[270];
	v[156] = -(v[132] * v[140]) + v[131] * v[144];
	v[145] = v[101] * v[358] + v[102] * v[360] + v[356] * v[99];
	v[284] = -(v[145] * v[242]) + v[140] * v[266] + v[133] * v[269] - v[131] * v[278];
	v[283] = -(v[145] * v[238]) + v[140] * v[262] + v[133] * v[268] - v[131] * v[277];
	v[282] = -(v[145] * v[234]) + v[140] * v[258] + v[133] * v[267] - v[131] * v[276];
	v[281] = v[145] * v[254] - v[144] * v[266] - v[133] * v[272] + v[132] * v[278];
	v[280] = v[145] * v[250] - v[144] * v[262] - v[133] * v[271] + v[132] * v[277];
	v[279] = v[145] * v[246] - v[144] * v[258] - v[133] * v[270] + v[132] * v[276];
	v[161] = -(v[133] * v[144]) + v[132] * v[145];
	v[299] = -(v[156] * v[269]) - v[140] * v[275] + v[161] * v[278] + v[145] * v[281];
	v[298] = -(v[156] * v[268]) - v[140] * v[274] + v[161] * v[277] + v[145] * v[280];
	v[297] = -(v[156] * v[267]) - v[140] * v[273] + v[161] * v[276] + v[145] * v[279];
	v[160] = v[133] * v[140] - v[131] * v[145];
	v[312] = v[160] * v[269] - v[161] * v[272] - v[144] * v[281] + v[140] * v[284];
	v[311] = v[160] * v[268] - v[161] * v[271] - v[144] * v[280] + v[140] * v[283];
	v[310] = v[160] * v[267] - v[161] * v[270] - v[144] * v[279] + v[140] * v[282];
	v[296] = v[156] * v[272] + v[144] * v[275] - v[160] * v[278] - v[145] * v[284];
	v[295] = v[156] * v[271] + v[144] * v[274] - v[160] * v[277] - v[145] * v[283];
	v[294] = v[156] * v[270] + v[144] * v[273] - v[160] * v[276] - v[145] * v[282];
	v[146] = v[355] * v[90] + v[357] * v[91] + v[359] * v[92];
	v[376] = v[146] * v[234];
	v[373] = v[146] * v[238];
	v[369] = v[146] * v[242];
	v[366] = -(v[131] * v[146]);
	v[150] = v[355] * v[94] + v[357] * v[96] + v[359] * v[97];
	v[377] = v[150] * v[246];
	v[374] = v[150] * v[250];
	v[370] = v[150] * v[254];
	v[367] = -(v[132] * v[150]);
	v[151] = v[101] * v[357] + v[102] * v[359] + v[355] * v[99];
	v[385] = v[151] * v[258];
	v[382] = v[151] * v[262];
	v[379] = v[151] * v[266];
	v[372] = -(v[133] * v[151]);
	v[368] = 2e0*v[151];
	v[155] = v[144] * v[156] - v[145] * v[160];
	v[157] = -(v[140] * v[156]) + v[145] * v[161];
	v[158] = (v[131] * v[131]);
	v[159] = (v[132] * v[132]);
	v[306] = v[158] + v[159];
	v[309] = -(v[157] * v[242]) + v[155] * v[254] + v[132] * v[296] - v[131] * v[299] + v[293] * v[306] + v[266] * v[366]
		+ v[266] * v[367] - v[133] * (v[266] * v[368] + v[369] + v[370] + v[371] + v[380]);
	v[308] = -(v[157] * v[238]) + v[155] * v[250] + v[132] * v[295] - v[131] * v[298] + v[292] * v[306] + v[262] * v[366]
		+ v[262] * v[367] - v[133] * (v[262] * v[368] + v[373] + v[374] + v[375] + v[383]);
	v[307] = -(v[157] * v[234]) + v[155] * v[246] + v[132] * v[294] - v[131] * v[297] + v[291] * v[306] + v[258] * v[366]
		+ v[258] * v[367] - v[133] * (v[258] * v[368] + v[376] + v[377] + v[378] + v[386]);
	v[165] = v[132] * (-(v[133] * v[150]) + v[155]) - v[131] * (v[133] * v[146] + v[157]) + v[151] * v[306];
	v[328] = -(v[165] * v[227]);
	v[162] = v[140] * v[160] - v[144] * v[161];
	v[163] = (v[133] * v[133]);
	v[320] = v[159] + v[163];
	v[323] = -(v[162] * v[254]) + v[157] * v[266] + v[133] * v[299] - v[132] * (v[150] * v[242] + v[312]) + v[287] * v[320]
		+ v[242] * v[372] - v[131] * (2e0*v[369] + v[370] + v[371] + v[379] + v[381]);
	v[322] = -(v[162] * v[250]) + v[157] * v[262] + v[133] * v[298] - v[132] * (v[150] * v[238] + v[311]) + v[286] * v[320]
		+ v[238] * v[372] - v[131] * (2e0*v[373] + v[374] + v[375] + v[382] + v[384]);
	v[321] = -(v[162] * v[246]) + v[157] * v[258] + v[133] * v[297] - v[132] * (v[150] * v[234] + v[310]) + v[285] * v[320]
		+ v[234] * v[372] - v[131] * (2e0*v[376] + v[377] + v[378] + v[385] + v[387]);
	v[316] = v[158] + v[163];
	v[319] = v[162] * v[242] - v[155] * v[266] - v[133] * v[296] + v[131] * v[312] + v[290] * v[316] + v[254] * v[366]
		+ v[254] * v[372] - v[132] * (v[369] + 2e0*v[370] + v[379] + v[380] + v[381]);
	v[318] = v[162] * v[238] - v[155] * v[262] - v[133] * v[295] + v[131] * v[311] + v[289] * v[316] + v[250] * v[366]
		+ v[250] * v[372] - v[132] * (v[373] + 2e0*v[374] + v[382] + v[383] + v[384]);
	v[317] = v[162] * v[234] - v[155] * v[258] - v[133] * v[294] + v[131] * v[310] + v[288] * v[316] + v[246] * v[366]
		+ v[246] * v[372] - v[132] * (v[376] + 2e0*v[377] + v[385] + v[386] + v[387]);
	v[167] = v[133] * (-(v[132] * v[151]) - v[155]) + v[131] * (-(v[132] * v[146]) + v[162]) + v[150] * v[316];
	v[333] = -(v[167] * v[224]);
	v[166] = v[133] * (-(v[131] * v[151]) + v[157]) + v[132] * (-(v[131] * v[150]) - v[162]) + v[146] * v[320];
	v[334] = v[166] * v[225];
	dT[0] = (*coef1)*(-((*a3)*ddui[0]) - (*a2)*dui[0] + (*a1)*ud[0]);
	dT[1] = (*coef1)*(-((*a3)*ddui[1]) - (*a2)*dui[1] + (*a1)*ud[1]);
	dT[2] = (*coef1)*(-((*a3)*ddui[2]) - (*a2)*dui[2] + (*a1)*ud[2]);
	dT[3] = (*coef2)*(-(v[165] * v[186]) - v[167] * v[192] + v[166] * v[87]);
	dT[4] = (*coef2)*(-(v[165] * v[189]) + v[166] * v[192] + v[167] * v[87]);
	dT[5] = (*coef2)*(v[166] * v[186] + v[167] * v[189] + v[165] * v[87]);
	DdT[0][0] = v[388];
	DdT[0][1] = 0e0;
	DdT[0][2] = 0e0;
	DdT[0][3] = 0e0;
	DdT[0][4] = 0e0;
	DdT[0][5] = 0e0;
	DdT[1][0] = 0e0;
	DdT[1][1] = v[388];
	DdT[1][2] = 0e0;
	DdT[1][3] = 0e0;
	DdT[1][4] = 0e0;
	DdT[1][5] = 0e0;
	DdT[2][0] = 0e0;
	DdT[2][1] = 0e0;
	DdT[2][2] = v[388];
	DdT[2][3] = 0e0;
	DdT[2][4] = 0e0;
	DdT[2][5] = 0e0;
	DdT[3][0] = 0e0;
	DdT[3][1] = 0e0;
	DdT[3][2] = 0e0;
	DdT[3][3] = (*coef2)*(v[166] * v[178] - v[186] * v[307] - v[192] * v[317] + v[328] + v[333] + v[321] * v[87]);
	DdT[3][4] = (*coef2)*(v[166] * v[180] - v[167] * v[225] - v[165] * v[229] - v[186] * v[308] - v[192] * v[318]
		+ v[322] * v[87]);
	DdT[3][5] = (*coef2)*(v[166] * v[181] + v[165] * v[225] - v[167] * v[226] - v[186] * v[309] - v[192] * v[319]
		+ v[323] * v[87]);
	DdT[4][0] = 0e0;
	DdT[4][1] = 0e0;
	DdT[4][2] = 0e0;
	DdT[4][3] = (*coef2)*(v[167] * v[178] + v[166] * v[224] - v[165] * v[230] - v[189] * v[307] + v[192] * v[321]
		+ v[317] * v[87]);
	DdT[4][4] = (*coef2)*(v[167] * v[180] - v[189] * v[308] + v[192] * v[322] - v[328] + v[334] + v[318] * v[87]);
	DdT[4][5] = (*coef2)*(v[167] * v[181] - v[165] * v[224] + v[166] * v[226] - v[189] * v[309] + v[192] * v[323]
		+ v[319] * v[87]);
	DdT[5][0] = 0e0;
	DdT[5][1] = 0e0;
	DdT[5][2] = 0e0;
	DdT[5][3] = (*coef2)*(v[165] * v[178] + v[166] * v[227] + v[167] * v[230] + v[189] * v[317] + v[186] * v[321]
		+ v[307] * v[87]);
	DdT[5][4] = (*coef2)*(v[165] * v[180] - v[167] * v[227] + v[166] * v[229] + v[189] * v[318] + v[186] * v[322]
		+ v[308] * v[87]);
	DdT[5][5] = (*coef2)*(v[165] * v[181] + v[189] * v[319] + v[186] * v[323] - v[333] - v[334] + v[309] * v[87]);
	*/
};
void Shell_1::EvaluateMassModal(double* v, double* alphai, double* coef1, double* coef2, double* coef3, double* e3r, double** matrixm)
{
	int i01; int i02;
	v[97] = -(*coef2) + (*coef3);
	v[90] = Power(alphai[2], 2);
	v[89] = 0.5e0*alphai[0] * alphai[2];
	v[88] = 0.5e0*alphai[1];
	v[87] = Power(alphai[1], 2);
	v[91] = v[87] + v[90];
	v[86] = alphai[0] * v[88];
	v[85] = Power(alphai[0], 2);
	v[68] = alphai[2] * v[88];
	v[55] = 4e0 / (4e0 + v[85] + v[91]);
	v[95] = e3r[0] * v[55];
	v[94] = e3r[1] * v[55];
	v[93] = -0.5e0*v[55];
	v[92] = e3r[2] * v[55];
	v[71] = (alphai[1] + v[89])*v[92] + e3r[0] * (1e0 + v[91] * v[93]) + (-alphai[2] + v[86])*v[94];
	v[78] = (v[71] * v[71]);
	v[72] = (-alphai[0] + v[68])*v[92] + e3r[1] * (1e0 + (v[85] + v[90])*v[93]) + (alphai[2] + v[86])*v[95];
	v[96] = v[72] * v[97];
	v[77] = (v[72] * v[72]);
	v[73] = e3r[2] * (1e0 + (v[85] + v[87])*v[93]) + (alphai[0] + v[68])*v[94] + (-alphai[1] + v[89])*v[95];
	v[79] = (v[73] * v[73]);
	matrixm[0][0] = (*coef1);
	matrixm[0][1] = 0e0;
	matrixm[0][2] = 0e0;
	matrixm[0][3] = 0e0;
	matrixm[0][4] = 0e0;
	matrixm[0][5] = 0e0;
	matrixm[1][1] = (*coef1);
	matrixm[1][2] = 0e0;
	matrixm[1][3] = 0e0;
	matrixm[1][4] = 0e0;
	matrixm[1][5] = 0e0;
	matrixm[2][2] = (*coef1);
	matrixm[2][3] = 0e0;
	matrixm[2][4] = 0e0;
	matrixm[2][5] = 0e0;
	matrixm[3][3] = (*coef3)*v[78] + (*coef2)*(v[77] + v[79]);
	matrixm[3][4] = v[71] * v[96];
	matrixm[3][5] = v[71] * v[73] * v[97];
	matrixm[4][4] = (*coef3)*v[77] + (*coef2)*(v[78] + v[79]);
	matrixm[4][5] = v[73] * v[96];
	matrixm[5][5] = (*coef2)*(v[77] + v[78]) + (*coef3)*v[79];
	for (i01 = 1; i01<6; i01++){
		for (i02 = 0; i02<i01; i02++){
			matrixm[i01][i02] = matrixm[i02][i01];
		}
	};
};

