#include "Beam_1.h"
#include"Database.h"

//Variáveis globais
extern
Database db;

Beam_1::Beam_1(void)
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
	VTK_nodes[0] = 0;
	VTK_nodes[1] = 2;
	VTK_nodes[2] = 1;

	DOFs = new int *[n_nodes];
	for (int i = 0; i < n_nodes; i++)
		DOFs[i] = new int[db.number_GLs_node];
	type_name = new char[20];//Nome do tipo do elemento
	sprintf(type_name, "Beam_1");
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
	for (int i = 0; i < 2;i++) 
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
	B1 = new Matrix(6,6);
	Qtransp = new Matrix(3,3);
	B2 = new Matrix(6,9);
	B2temp = new Matrix(3, 3);
	stiffness = new Matrix(18, 18);
	mass = new Matrix(18, 18);
	mass_modal = new Matrix(18, 18);
	damping = new Matrix(18, 18);
	damping_modal = new Matrix(18, 18);
	rayleigh_damping = new Matrix(18, 18);
	i_loading = new Matrix(18);
	inertial_loading = new Matrix(18);
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
	e3rg = new Matrix(3, 1);

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
	DdT = new Matrix(6, 6);
	dT = new Matrix(6);

	tempkin = new double;

	T0 = 0.0;//Pre tension
}

Beam_1::~Beam_1(void)
{
	delete [] nodes;
	delete[] VTK_nodes;
	if (DOFs != NULL)
	{
		for (int i = 0; i < n_nodes; i++)
			delete[] DOFs[i];
		delete[] DOFs;
	}
	delete [] N1;
	delete [] N2;
	delete [] N3;
	delete [] dN1;
	delete [] dN2;
	delete [] dN3;
	delete [] csi;
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
	delete e3rg;

	delete br;
	delete DdT;
	delete dT;

	delete tempkin;

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
}

//Checa inconsistências no elemento para evitar erros de execução
bool Beam_1::Check()
{
	for (int i = 0; i < n_nodes; i++)
	{
		if (nodes[i] > db.number_nodes)
			return false;
	}	
	if (cs > db.number_CS)
		return false;
	if (section > db.number_sections)
		return false;
	if (material > db.number_materials)
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

bool Beam_1::Read(FILE *f)
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

	//Salva a posição (stream)
	fpos_t pos;
	fgetpos(f, &pos);
	fscanf(f, "%s", s);
	if (!strcmp(s, "PreTension"))
	{
		fscanf(f, "%s", s);
		T0 = atof(s);
	}
	else
		fsetpos(f, &pos);

	return true;
}

void Beam_1::Write(FILE *f)
{
	fprintf(f, "Beam_1\t%d\tMat\t%d\tSec\t%d\tCS\t%d\tNodes\t%d\t%d\t%dPreTension\t%.6e\n",
		number,
		material,
		section,
		cs,
		nodes[0],
		nodes[1],
		nodes[2],
		T0);
}
//Escreve arquivo de resultados
void Beam_1::WriteResults(FILE *f)
{
	fprintf(f, "Beam_1\t%d\tStress_r\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\tStrain_r\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\tRef_Coor\t%.6e\t%.6e\t%.6e\tCoor\t%.6e\t%.6e\t%.6e\n",
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

void Beam_1::WriteVTK_XMLBase(std::vector<float> *float_vector)
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
void Beam_1::WriteVTK_XMLRender(FILE *f)
{
	db.sections[section - 1]->sec_details->WriteVTK_XMLRender(f, this);
}

//Escreve no monitor do elemento//Escreve no monitor do elemento
void Beam_1::WriteMonitor(FILE *f, bool first_record, double time)
{
	Matrix e1_0 = transp(*transform3)*(*lag_save->Q_i[0])*(*transform3)*(*db.CS[cs - 1]->E1);
	Matrix e2_0 = transp(*transform3)*(*lag_save->Q_i[0])*(*transform3)*(*db.CS[cs - 1]->E2);
	Matrix e3_0 = transp(*transform3)*(*lag_save->Q_i[0])*(*transform3)*(*db.CS[cs - 1]->E3);
	Matrix e1_1 = transp(*transform3)*(*lag_save->Q_i[1])*(*transform3)*(*db.CS[cs - 1]->E1);
	Matrix e2_1 = transp(*transform3)*(*lag_save->Q_i[1])*(*transform3)*(*db.CS[cs - 1]->E2);
	Matrix e3_1 = transp(*transform3)*(*lag_save->Q_i[1])*(*transform3)*(*db.CS[cs - 1]->E3);
	//Cabeçalho
	if (first_record == true)
		fprintf(f, "TIME\tnx1r\tny1r\tnz1r\tmx1r\tmy1r\tmz1r\tnx2r\tny2r\tnz2r\tmx2r\tmy2r\tmz2r\tex1\t\t\tey1\t\t\tez1\t\t\tex2\t\t\tey2\t\t\tez2\t\t\t\n");
	//Informações a serem salvas
	fprintf(f, "%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t\n",
		time,
		(*sigma_r[0])(0, 0), (*sigma_r[0])(1, 0), (*sigma_r[0])(2, 0), (*sigma_r[0])(3, 0), (*sigma_r[0])(4, 0), (*sigma_r[0])(5, 0),
		(*sigma_r[1])(0, 0), (*sigma_r[1])(1, 0), (*sigma_r[1])(2, 0), (*sigma_r[1])(3, 0), (*sigma_r[1])(4, 0), (*sigma_r[1])(5, 0),
		e1_0(0, 0), e1_0(1, 0), e1_0(2, 0), e2_0(0, 0), e2_0(1, 0), e2_0(2, 0), e3_0(0, 0), e3_0(1, 0), e3_0(2, 0),
		e1_1(0, 0), e1_1(1, 0), e1_1(2, 0), e2_1(0, 0), e2_1(1, 0), e2_1(2, 0), e3_1(0, 0), e3_1(1, 0), e3_1(2, 0)
		);
}

//Pré-cálculo de variáveis que é feito uma única vez no início
void Beam_1::PreCalc()
{
	bool special_section = false;
	//Operador constitutivo - cast feito para acessar membros da classe filha
	//Se a seção transversal for "UserDefined"
	if (typeid(*db.sections[section - 1]) == typeid(SecUserDefined))
	{
		SecUserDefined* ptr_section = static_cast<SecUserDefined*>(db.sections[section - 1]);	//ptr_section is a pointer to the User Defined Section
		(*D)(0, 0) = ptr_section->GA;
		(*D)(1, 1) = ptr_section->GA;
		(*D)(2, 2) = ptr_section->EA;
		(*D)(3, 3) = ptr_section->EI11;
		(*D)(4, 4) = ptr_section->EI22;
		(*D)(5, 5) = ptr_section->GJT;
		(*D)(0, 5) = ptr_section->GS1S - ptr_section->GS1;
		(*D)(5, 0) = ptr_section->GS1S - ptr_section->GS1;
		(*D)(1, 5) = ptr_section->GS2S - ptr_section->GS2;
		(*D)(5, 1) = ptr_section->GS2S - ptr_section->GS2;
		(*D)(2, 3) = ptr_section->ES1;
		(*D)(3, 2) = ptr_section->ES1;
		(*D)(2, 4) = ptr_section->ES2;
		(*D)(4, 2) = ptr_section->ES2;
		(*D)(3, 4) = ptr_section->EI12;
		(*D)(4, 3) = ptr_section->EI12;
		
		//Massa por unidade de comprimento na configuração de referência
		(*Mr)(0, 0) = ptr_section->Rho;
		(*Mr)(1, 1) = ptr_section->Rho;
		(*Mr)(2, 2) = ptr_section->Rho;
		Mr->MatrixToPtr(pMr, 3);//salvando ptr para uso na dinâmica
		//Inércia por unidade de comprimento na configuração de referência
		(*Jr)(0, 0) = ptr_section->J11;
		(*Jr)(1, 1) = ptr_section->J22;
		(*Jr)(2, 2) = ptr_section->J11 + ptr_section->J22;
		(*Jr)(0, 1) = ptr_section->J12;
		(*Jr)(1, 0) = ptr_section->J12;
		Jr->MatrixToPtr(pJr, 3);//salvando ptr para uso na dinâmica
		//Termos de transporte - alterar em situações em que o baricentro não é o polo utilizado no cálculo do momento de inércia
		(*br)(0, 0) = ptr_section->BC(0, 0);
		(*br)(1, 0) = ptr_section->BC(1, 0);
		(*br)(2, 0) = 0.0;
		special_section = true;
	}
	//Se a seção transversal for "SecHelicalFiber"
	if (typeid(*db.sections[section - 1]) == typeid(SecHelicalFiber))
	{
		SecHelicalFiber* ptr_section = static_cast<SecHelicalFiber*>(db.sections[section - 1]);	//ptr_section is a pointer to the SecHelicalFiber Section
		ptr_section->ComputeStiffnessMass(*D, *Mr, *Jr, db.materials[material - 1]->number);
		/*D->print();
		Mr->print();
		Jr->print();*/
		Mr->MatrixToPtr(pMr, 3);//salvando ptr para uso na dinâmica
		Jr->MatrixToPtr(pJr, 3);//salvando ptr para uso na dinâmica
		//Termos de transporte
		(*br)(0, 0) = 0.0;
		(*br)(1, 0) = 0.0;
		(*br)(2, 0) = 0.0;
		special_section = true;
	}
	if (special_section == false)
	{
		Hooke* hooke = static_cast<Hooke*>(db.materials[material - 1]);
		double E = hooke->E;
		double nu = hooke->nu;
		double G = E / (2 * (1 + nu));
		//double sf = 10 * (1 + nu) / (12 + 11 * nu);
		double sf = 1.0;
		double A = db.sections[section - 1]->A;
		double I1 = db.sections[section - 1]->I11;
		double I2 = db.sections[section - 1]->I22;
		double I12 = db.sections[section - 1]->I12;
		double It = db.sections[section - 1]->It;
		(*D)(0, 0) = sf*G*A;
		(*D)(1, 1) = sf*G*A;
		(*D)(2, 2) = E*A;
		(*D)(3, 3) = E*I1;
		(*D)(4, 4) = E*I2;
		(*D)(3, 4) = E*I12;
		(*D)(4, 3) = E*I12;
		(*D)(5, 5) = G*It;
		//Massa por unidade de comprimento na configuração de referência
		(*Mr)(0, 0) = db.materials[material - 1]->rho*db.sections[section - 1]->A;
		(*Mr)(1, 1) = db.materials[material - 1]->rho*db.sections[section - 1]->A;
		(*Mr)(2, 2) = db.materials[material - 1]->rho*db.sections[section - 1]->A;
		Mr->MatrixToPtr(pMr, 3);//salvando ptr para uso na dinâmica
		//Inércia por unidade de comprimento na configuração de referência
		(*Jr)(0, 0) = hooke->rho*db.sections[section - 1]->I11;
		(*Jr)(1, 1) = hooke->rho*db.sections[section - 1]->I22;
		(*Jr)(2, 2) = hooke->rho*db.sections[section - 1]->I33;
		(*Jr)(0, 1) = hooke->rho*db.sections[section - 1]->I12;
		(*Jr)(1, 0) = hooke->rho*db.sections[section - 1]->I12;
		Jr->MatrixToPtr(pJr, 3);//salvando ptr para uso na dinâmica
		//Termos de transporte - alterar em situações em que o baricentro não é o polo utilizado no cálculo do momento de inércia
		(*br)(0, 0) = 0.0;
		(*br)(1, 0) = 0.0;
		(*br)(2, 0) = 0.0;
	}
	
	//Identidade de ordem 3
	(*I3)(0, 0) = 1.0;
	(*I3)(1, 1) = 1.0;
	(*I3)(2, 2) = 1.0;
	(*B2)(0, 0) = 1.0;
	(*B2)(1, 1) = 1.0;
	(*B2)(2, 2) = 1.0;
	//Monta as matrizes de rotação
	TransformMatrix();
	//Versor e3r (configuração de referência)
	(*e3r)(0, 0) = db.nodes[nodes[2] - 1]->ref_coordinates[0] - db.nodes[nodes[0] - 1]->ref_coordinates[0];
	(*e3r)(1, 0) = db.nodes[nodes[2] - 1]->ref_coordinates[1] - db.nodes[nodes[0] - 1]->ref_coordinates[1];
	(*e3r)(2, 0) = db.nodes[nodes[2] - 1]->ref_coordinates[2] - db.nodes[nodes[0] - 1]->ref_coordinates[2];	
	*e3r = (1.0 / (norm(*e3r)))*(*e3r);							//Normalização unitária
	*e3rg = *e3r;												//Cópia do e3r (preservado no sistema global)
	*e3r = (*transform3)*(*e3r);

	//Pre-strain (due to pre-tension T0)
	double du0 = T0/(*D)(2, 2);
	(*lag_save->dz_i[0])(2, 0) = 1.0 + du0;
	(*lag_save->dz_i[1])(2, 0) = 1.0 + du0;
	//Reference (undeformed) length
	length = CalculateLength()/(1.0+du0);

	jacobian = length / 2.0;
	alpha1 = 1.0; //Peso do método de quadratura Gaussiana
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
	}//end of gauss
}

//Monta elementos
void Beam_1::Mount()
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
		(*d_alpha_delta[gauss])(0, 0) =		(db.nodes[nodes[0] - 1]->displacements[3])* dN1[gauss] +
											(db.nodes[nodes[1] - 1]->displacements[3])* dN2[gauss] +
											(db.nodes[nodes[2] - 1]->displacements[3])* dN3[gauss] ;
		(*d_alpha_delta[gauss])(1, 0) =		(db.nodes[nodes[0] - 1]->displacements[4])* dN1[gauss] +
											(db.nodes[nodes[1] - 1]->displacements[4])* dN2[gauss] +
											(db.nodes[nodes[2] - 1]->displacements[4])* dN3[gauss] ;
		(*d_alpha_delta[gauss])(2, 0) =		(db.nodes[nodes[0] - 1]->displacements[5])* dN1[gauss] +
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
		(*d_u_delta[gauss])(0, 0) = (db.nodes[nodes[0] - 1]->displacements[0])* dN1[gauss] +
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
		//Energia de deformação
		strain_energy += 0.5 * (alpha1*jacobian) * (transp(*sigma_r[gauss]) * (*epsilon_r[gauss]))(0,0);
	}//end of gauss loop
	//Transformações de coordenadas, para o sistema global - matriz de rigidez
	(*stiffness) = (transp(*transform)*(*stiffness))*(*transform);
	(*i_loading) = transp(*transform)*(*i_loading);
}

//Monta carregamentos associados ao elemento
void Beam_1::MountElementLoads()
{
	MountFieldLoading();
	(*P_loading) = (*i_loading) - (*e_loading);			//Vetor esforço desbalanceado
	(*stiffness) = (*stiffness) - (*loading_stiffness);	//Rigidez de carregamentos
}

//Monta carregamentos de campo
void Beam_1::MountFieldLoading()
{
	if (db.environment != NULL)
	{
		//Se existe campo gravitacional
		if (db.environment->g_exist == true)
		{
			load_multiplier = 1.0;
			l_factor = db.environment->bool_g.GetLinearFactorAtCurrentTime();
			
			//Seção UserDefined - calcula o peso próprio com o rho por unidade de comprimento que vem da seção (não usa propriedade do material)
			if (typeid(*db.sections[section - 1]) == typeid(SecUserDefined))
			{
				SecUserDefined* ptr_section = static_cast<SecUserDefined*>(db.sections[section - 1]);	//ptr_section is a pointer to the User Defined Section
				mult = l_factor*load_multiplier*alpha1*jacobian*ptr_section->Rho;
			}
			else
				mult = l_factor*load_multiplier*alpha1*jacobian*(db.materials[material - 1]->rho)*db.sections[section - 1]->A;
			//////////////////Efeito do peso próprio////////////////////////////
			//Loop nos pontos de Gauss
			for (int gauss = 0; gauss < 2; gauss++)
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

				//Posição do ponto de Gauss
				Matrix pos_gauss(3);
				pos_gauss(0, 0) = (db.nodes[nodes[0] - 1]->copy_coordinates[0])* N1[gauss] +
					(db.nodes[nodes[1] - 1]->copy_coordinates[0])* N2[gauss] +
					(db.nodes[nodes[2] - 1]->copy_coordinates[0])* N3[gauss] +
					(db.nodes[nodes[0] - 1]->displacements[0])* N1[gauss] +
					(db.nodes[nodes[1] - 1]->displacements[0])* N2[gauss] +
					(db.nodes[nodes[2] - 1]->displacements[0])* N3[gauss];
				pos_gauss(1, 0) = (db.nodes[nodes[0] - 1]->copy_coordinates[1])* N1[gauss] +
					(db.nodes[nodes[1] - 1]->copy_coordinates[1])* N2[gauss] +
					(db.nodes[nodes[2] - 1]->copy_coordinates[1])* N3[gauss] +
					(db.nodes[nodes[0] - 1]->displacements[1])* N1[gauss] +
					(db.nodes[nodes[1] - 1]->displacements[1])* N2[gauss] +
					(db.nodes[nodes[2] - 1]->displacements[1])* N3[gauss];
				pos_gauss(2, 0) = (db.nodes[nodes[0] - 1]->copy_coordinates[2])* N1[gauss] +
					(db.nodes[nodes[1] - 1]->copy_coordinates[2])* N2[gauss] +
					(db.nodes[nodes[2] - 1]->copy_coordinates[2])* N3[gauss] +
					(db.nodes[nodes[0] - 1]->displacements[2])* N1[gauss] +
					(db.nodes[nodes[1] - 1]->displacements[2])* N2[gauss] +
					(db.nodes[nodes[2] - 1]->displacements[2])* N3[gauss];
				//Velocidade da correnteza no ponto de Gauss (cota z)
				//Verifica profundidade
				double hg = -1.0 * dot(pos_gauss , db.environment->G);
				potential_gravitational_energy += mult * hg;
			}//end of gauss loop
			//Verifica se é necessário inserir o carregamento do momento induzido pelo peso (eixo fora do baricentro)
			if (typeid(*db.sections[section - 1]) == typeid(SecUserDefined) && norm(*br) != 0.0)
			{
				//Definindo vetores v, alphai, alphad, brglobal, gamglobal, dwe, ddwe:
				double* v = temp_v;
				double alphad[3];
				double alphai[3];
				SecUserDefined* ptr_section = static_cast<SecUserDefined*>(db.sections[section - 1]);	//ptr_section is a pointer to the User Defined Section
				Matrix mbrglobal = ptr_section->BC(0, 0)*(*db.CS[cs - 1]->E1) + ptr_section->BC(1, 0)*(*db.CS[cs - 1]->E2);
				double brglobal[3];
				brglobal[0] = mbrglobal(0, 0);
				brglobal[1] = mbrglobal(1, 0);
				brglobal[2] = mbrglobal(2, 0);
				double gamglobal[3];
				gamglobal[0] = mult*db.environment->G(0, 0); //já tem o jacobiano, peso de integração e fatores de carregamento embutidos
				gamglobal[1] = mult*db.environment->G(1, 0); //já tem o jacobiano, peso de integração e fatores de carregamento embutidos
				gamglobal[2] = mult*db.environment->G(2, 0); //já tem o jacobiano, peso de integração e fatores de carregamento embutidos
				
				double dwe[3];
				double ddwe[3][3];
				
				//Loop nos pontos de Gauss
				for (int gauss = 0; gauss < 2; gauss++)
				{
					for (int i = 0; i < 3; i++)
					{
						//alphad calculado no ponto de Gauss							
						alphad[i] = (db.nodes[nodes[0] - 1]->displacements[3 + i])* N1[gauss] +
							(db.nodes[nodes[1] - 1]->displacements[3 + i])* N2[gauss] +
							(db.nodes[nodes[2] - 1]->displacements[3 + i])* N3[gauss];
						//alphai calculado no ponto de Gauss							
						alphai[i] = (db.nodes[nodes[0] - 1]->copy_coordinates[3 + i])* N1[gauss] +
							(db.nodes[nodes[1] - 1]->copy_coordinates[3 + i])* N2[gauss] +
							(db.nodes[nodes[2] - 1]->copy_coordinates[3 + i])* N3[gauss];
					}
					//AceGen
					#pragma region acegen
					v[219] = 0.5e0*alphad[2];
					v[217] = Power(alphad[2], 2);
					v[216] = alphad[0] * v[219];
					v[215] = 0.5e0*alphad[1];
					v[214] = 2e0*alphad[2];
					v[213] = Power(alphad[1], 2);
					v[212] = alphad[0] * v[215];
					v[211] = 2e0*alphad[1];
					v[210] = Power(alphad[0], 2);
					v[209] = 2e0*alphad[0];
					v[208] = Power(alphai[2], 2);
					v[207] = 0.5e0*alphai[0] * alphai[2];
					v[206] = 0.5e0*alphai[1];
					v[205] = Power(alphai[1], 2);
					v[223] = v[205] + v[208];
					v[204] = alphai[0] * v[206];
					v[203] = Power(alphai[0], 2);
					v[63] = alphai[2] * v[206];
					v[135] = -v[210] - v[213];
					v[113] = alphad[2] + v[212];
					v[103] = -alphad[2] + v[212];
					v[47] = alphad[2] * v[215];
					v[130] = alphad[0] + v[47];
					v[122] = -alphad[0] + v[47];
					v[126] = -alphad[1] + v[216];
					v[109] = alphad[1] + v[216];
					v[117] = -v[210] - v[217];
					v[98] = -v[213] - v[217];
					v[93] = 4e0 + v[210] + v[213] + v[217];
					v[95] = 1e0 / Power(v[93], 2);
					v[218] = -4e0*v[95];
					v[97] = v[214] * v[218];
					v[222] = 0.5e0*v[97];
					v[139] = v[135] * v[222];
					v[96] = v[211] * v[218];
					v[221] = 0.5e0*v[96];
					v[119] = v[117] * v[221];
					v[94] = v[209] * v[218];
					v[220] = 0.5e0*v[94];
					v[143] = v[215] * v[94];
					v[140] = -(v[219] * v[94]);
					v[99] = v[220] * v[98];
					v[34] = 4e0 / v[93];
					v[144] = -0.5e0*v[34];
					v[146] = v[144] - alphad[0] * v[220];
					v[137] = v[144] * v[211];
					v[138] = v[137] + v[135] * v[221];
					v[134] = v[144] * v[209];
					v[136] = v[134] + v[135] * v[220];
					v[131] = v[34] + v[130] * v[94];
					v[128] = -v[34] + v[126] * v[96];
					v[123] = -v[34] + v[122] * v[94];
					v[120] = v[144] * v[214];
					v[121] = v[120] + v[117] * v[222];
					v[118] = v[134] + v[117] * v[220];
					v[116] = v[34] + v[113] * v[97];
					v[111] = v[34] + v[109] * v[96];
					v[108] = alphad[2] * v[144];
					v[132] = -v[108] + v[130] * v[96];
					v[127] = -v[108] + v[126] * v[94];
					v[124] = -v[108] + v[122] * v[96];
					v[110] = -v[108] + v[109] * v[94];
					v[107] = -v[34] + v[103] * v[97];
					v[105] = alphad[0] * v[144];
					v[129] = -v[105] + v[126] * v[97];
					v[115] = -v[105] + v[113] * v[96];
					v[112] = -v[105] + v[109] * v[97];
					v[106] = -v[105] + v[103] * v[96];
					v[102] = -(alphad[1] * v[144]);
					v[133] = v[102] + v[130] * v[97];
					v[125] = v[102] + v[122] * v[97];
					v[114] = v[102] + v[113] * v[94];
					v[104] = v[102] + v[103] * v[94];
					v[101] = v[120] + v[222] * v[98];
					v[100] = v[137] + v[221] * v[98];
					v[37] = 1e0 - v[144] * v[98];
					v[38] = v[103] * v[34];
					v[39] = v[109] * v[34];
					v[41] = v[113] * v[34];
					v[43] = 1e0 - v[117] * v[144];
					v[44] = v[122] * v[34];
					v[46] = v[126] * v[34];
					v[48] = v[130] * v[34];
					v[49] = 1e0 - v[135] * v[144];
					v[50] = 4e0 / (4e0 + v[203] + v[223]);
					v[224] = -0.5e0*v[50];
					v[53] = 1e0 + v[223] * v[224];
					v[54] = (-alphai[2] + v[204])*v[50];
					v[55] = (alphai[1] + v[207])*v[50];
					v[57] = (alphai[2] + v[204])*v[50];
					v[59] = 1e0 + (v[203] + v[208])*v[224];
					v[60] = v[50] * (-alphai[0] + v[63]);
					v[62] = (-alphai[1] + v[207])*v[50];
					v[64] = v[50] * (alphai[0] + v[63]);
					v[65] = 1e0 + (v[203] + v[205])*v[224];
					v[185] = brglobal[0] * (v[129] * v[53] + v[133] * v[57] + v[139] * v[62]) + brglobal[1] * (v[129] * v[54] + v[133] * v[59]
						+ v[139] * v[64]) + brglobal[2] * (v[129] * v[55] + v[133] * v[60] + v[139] * v[65]);
					v[181] = brglobal[0] * (v[128] * v[53] + v[132] * v[57] + v[138] * v[62]) + brglobal[1] * (v[128] * v[54] + v[132] * v[59]
						+ v[138] * v[64]) + brglobal[2] * (v[128] * v[55] + v[132] * v[60] + v[138] * v[65]);
					v[177] = brglobal[0] * (v[127] * v[53] + v[131] * v[57] + v[136] * v[62]) + brglobal[1] * (v[127] * v[54] + v[131] * v[59]
						+ v[136] * v[64]) + brglobal[2] * (v[127] * v[55] + v[131] * v[60] + v[136] * v[65]);
					v[170] = brglobal[0] * (v[116] * v[53] + v[121] * v[57] + v[125] * v[62]) + brglobal[1] * (v[116] * v[54] + v[121] * v[59]
						+ v[125] * v[64]) + brglobal[2] * (v[116] * v[55] + v[121] * v[60] + v[125] * v[65]);
					v[166] = brglobal[0] * (v[115] * v[53] + v[119] * v[57] + v[124] * v[62]) + brglobal[1] * (v[115] * v[54] + v[119] * v[59]
						+ v[124] * v[64]) + brglobal[2] * (v[115] * v[55] + v[119] * v[60] + v[124] * v[65]);
					v[187] = gamglobal[2] * v[166] - gamglobal[1] * v[181];
					v[162] = brglobal[0] * (v[114] * v[53] + v[118] * v[57] + v[123] * v[62]) + brglobal[1] * (v[114] * v[54] + v[118] * v[59]
						+ v[123] * v[64]) + brglobal[2] * (v[114] * v[55] + v[118] * v[60] + v[123] * v[65]);
					v[186] = gamglobal[2] * v[162] - gamglobal[1] * v[177];
					v[158] = brglobal[0] * (v[101] * v[53] + v[107] * v[57] + v[112] * v[62]) + brglobal[1] * (v[101] * v[54] + v[107] * v[59]
						+ v[112] * v[64]) + brglobal[2] * (v[101] * v[55] + v[107] * v[60] + v[112] * v[65]);
					v[154] = brglobal[0] * (v[100] * v[53] + v[106] * v[57] + v[111] * v[62]) + brglobal[1] * (v[100] * v[54] + v[106] * v[59]
						+ v[111] * v[64]) + brglobal[2] * (v[100] * v[55] + v[106] * v[60] + v[111] * v[65]);
					v[190] = -(gamglobal[2] * v[154]) + gamglobal[0] * v[181];
					v[172] = gamglobal[1] * v[154] - gamglobal[0] * v[166];
					v[150] = brglobal[0] * (v[104] * v[57] + v[110] * v[62] + v[53] * v[99]) + brglobal[1] * (v[104] * v[59] + v[110] * v[64]
						+ v[54] * v[99]) + brglobal[2] * (v[104] * v[60] + v[110] * v[65] + v[55] * v[99]);
					v[189] = -(gamglobal[2] * v[150]) + gamglobal[0] * v[177];
					v[171] = gamglobal[1] * v[150] - gamglobal[0] * v[162];
					v[78] = brglobal[0] * (v[37] * v[53] + v[38] * v[57] + v[39] * v[62]) + brglobal[1] * (v[37] * v[54] + v[38] * v[59]
						+ v[39] * v[64]) + brglobal[2] * (v[37] * v[55] + v[38] * v[60] + v[39] * v[65]);
					v[79] = brglobal[0] * (v[41] * v[53] + v[43] * v[57] + v[44] * v[62]) + brglobal[1] * (v[41] * v[54] + v[43] * v[59]
						+ v[44] * v[64]) + brglobal[2] * (v[41] * v[55] + v[43] * v[60] + v[44] * v[65]);
					v[82] = gamglobal[1] * v[78] - gamglobal[0] * v[79];
					v[195] = -(v[143] * v[82]);
					v[80] = brglobal[0] * (v[46] * v[53] + v[48] * v[57] + v[49] * v[62]) + brglobal[1] * (v[46] * v[54] + v[48] * v[59]
						+ v[49] * v[64]) + brglobal[2] * (v[46] * v[55] + v[48] * v[60] + v[49] * v[65]);
					v[84] = gamglobal[2] * v[79] - gamglobal[1] * v[80];
					v[199] = -(alphad[2] * v[221] * v[84]);
					v[83] = -(gamglobal[2] * v[78]) + gamglobal[0] * v[80];
					v[198] = -(v[140] * v[83]);
					v[193] = -(v[105] * v[171]) + v[108] * v[186] + v[189] * v[34] - v[146] * v[82] + v[140] * v[84] + v[83] * v[94];
					v[194] = v[102] * v[186] + v[105] * v[189] + v[171] * v[34] + v[146] * v[83] + v[143] * v[84] + v[82] * v[94];
					v[197] = v[102] * v[187] + v[105] * v[190] + v[172] * v[34] - v[143] * v[83] + (-v[144] + alphad[1] * v[221])*v[84]
						+ v[82] * v[96];
					dwe[0] = -(v[102] * v[82]) - v[108] * v[83] + v[34] * v[84];
					dwe[1] = -(v[105] * v[82]) + v[34] * v[83] + v[108] * v[84];
					dwe[2] = v[34] * v[82] + v[105] * v[83] + v[102] * v[84];
					ddwe[0][0] = -(v[102] * v[171]) - v[108] * v[189] + v[195] + v[198] + v[186] * v[34] + v[84] * v[94];
					ddwe[0][1] = v[193];
					ddwe[0][2] = v[194];
					ddwe[1][0] = v[193];
					ddwe[1][1] = -(v[105] * v[172]) + v[108] * v[187] - v[195] + v[199] + v[190] * v[34] + v[83] * v[96];
					ddwe[1][2] = v[197];
					ddwe[2][0] = v[194];
					ddwe[2][1] = v[197];
					ddwe[2][2] = v[105] * (-(gamglobal[2] * v[158]) + gamglobal[0] * v[185]) + v[102] * (gamglobal[2] * v[170]
						- gamglobal[1] * v[185]) - v[198] - v[199] + (gamglobal[1] * v[158] - gamglobal[0] * v[170])*v[34] + v[82] * v[97];
#pragma endregion				
					
					//Contribuição para carregamento e rigidez
					Matrix dwe6(6);		//Matriz auxiliar
					Matrix ddwe6(6, 6);	//Matriz auxiliar
					for (int i = 0; i < 3; i++)
					{
						dwe6(3 + i, 0) = dwe[i];
						for (int j = 0; j < 3; j++)
							ddwe6(3 + i, 3 + j) = ddwe[i][j];
					}

					(*e_loading) = (*e_loading) + transp(*N[gauss])*dwe6;
					(*loading_stiffness) = (*loading_stiffness) + transp(*N[gauss])*ddwe6*(*N[gauss]);
				}
			}
		}
		//Se existe carregamento de vento
		if (db.environment->wind_data_exist == true)
		{
			//Se há dados de curvas aerodinâmicas associadas à seção transversal em questão - calcula esforços aerodinâmicos
			if (db.sections[section - 1]->aerodynamicdataID != 0)
			{
				load_multiplier = 1.0;
				l_factor = 1.0;

				//Loop nos pontos de Gauss
 				for (int gauss = 0; gauss < 2; gauss++)
				{
					//Variáveis para rotações
					//incremental - Delta
					double alpha_e;
					Matrix A_d(3,3);
					Matrix Q_d(3,3);
					Matrix Xi_d(3,3);
					//"i" - até a configuração anterior convergida
					double alpha_e_i;
					Matrix A_i(3, 3);
					Matrix Q_i(3, 3);
					//Variáveis para facilitar acesso:
					double rho_air = db.environment->rho_air;					//rho ar
					double chord = db.sections[section - 1]->aero_length;		//corda
					//Cálculo da posição do ponto de Gauss para avaliação do vetor velocidade nesse ponto
					Matrix x_Gauss(3);
					//Velocidade do ponto de Gauss - tomada no início do time step
					Matrix u_gauss(3);
					//Rotação incremental (para o cálculo do pseudo-momento aerodinâmico)
					Matrix alpha_d(3);
					//Rotação total (para o cálculo do tensor rotação para atualizar orientações da ST)
					Matrix alpha_i(3);
					//Avaliação de posição do ponto de Gauss na configuração "i"
					for (int i = 0; i < 3; i++)
					{
						//x_Gauss calculado no ponto de Gauss
						x_Gauss(i, 0) = (db.nodes[nodes[0] - 1]->copy_coordinates[i]) * N1[gauss]
							+ (db.nodes[nodes[1] - 1]->copy_coordinates[i]) * N2[gauss]
							+ (db.nodes[nodes[2] - 1]->copy_coordinates[i]) * N3[gauss];
						//u_Gauss calculado no ponto de Gauss
						u_gauss(i, 0) = (db.nodes[nodes[0] - 1]->copy_vel[i]) * N1[gauss]
							+ (db.nodes[nodes[1] - 1]->copy_vel[i]) * N2[gauss]
							+ (db.nodes[nodes[2] - 1]->copy_vel[i]) * N3[gauss];
						//alpha_d calculado no ponto de Gauss							
						alpha_d(i, 0) = (db.nodes[nodes[0] - 1]->displacements[3+i])* N1[gauss] +
							(db.nodes[nodes[1] - 1]->displacements[3+i])* N2[gauss] +
							(db.nodes[nodes[2] - 1]->displacements[3+i])* N3[gauss];
						//alpha_i calculado no ponto de Gauss							
						alpha_i(i, 0) = (db.nodes[nodes[0] - 1]->copy_coordinates[3 + i])* N1[gauss] +
							(db.nodes[nodes[1] - 1]->copy_coordinates[3 + i])* N2[gauss] +
							(db.nodes[nodes[2] - 1]->copy_coordinates[3 + i])* N3[gauss];
					}
					//Tensor rotação "i"
					alpha_e_i = norm(alpha_i);							//Valor escalar do parametro alpha
					A_i = skew(alpha_i);								//Matriz A_i
					g = 4.0 / (4.0 + alpha_e_i*alpha_e_i);				//função g(alpha)
					Q_i = *I3 + g*(A_i + 0.5*A_i*A_i);					//Tensor de rotação

					//Velocidade ao longe do vento - tomada no início do time step
					Matrix u_inf(3);
					u_inf = db.environment->WindVelocityAt(x_Gauss, db.last_converged_time);
					//Velocidade relativa
					Matrix u_rel(3);
					//Atual direção da linha de corda - direção e1 atualizada - tomada no início do time step
					Matrix e1 = Q_i*(*db.CS[cs - 1]->E1);
					//Atual direção ortogonal à da linha de corda - direção e2 atualizada - tomada no início do time step
					Matrix e2 = Q_i*(*db.CS[cs - 1]->E2);
					//Atual direção do eixo da barra - direção e3 atualizada - tomada no início do time step
					Matrix e3 = Q_i*(*db.CS[cs - 1]->E3);
					//Forças e momentos - vetoriais
					Matrix drag(3);
					Matrix lift(3);
					Matrix moment(3);
					double alpha;//ângulo de ataque
					double CL, CD, CM;//coeficientes aerodinâmicos
					
					///////////////////////////////////Se não é definido o BEM/////////////////////////////////////////
					//***********************************************************************************************//
					///////////////////////////////////////////////////////////////////////////////////////////////////
					if (db.bem_exist == false)
					{
						//Velocidade relativa - incidente no ponto de Gauss
						u_rel = u_inf - u_gauss;
						//Correção da velocidade relativa - retirando componente na direção da barra
						u_rel = u_rel - dot(u_rel, e3)*e3;
						//Cálculo do ângulo de ataque - calculado utilizando-se a orientação atual da ST no início do time step, direção do vento ambiente incidente e velocidade
						double alpha = 0.0;
						if (norm(u_rel) != 0.0)
						{
							double sin_alpha = dot(cross(e1, u_rel), e3) / norm(u_rel);
							//Verificações de erros de arredondamento para evitar erros na função asin (argumento fora do range esperado)
							if (sin_alpha > 1.0)
								sin_alpha = 1.0;
							if (sin_alpha < -1.0)
								sin_alpha = -1.0;
							alpha = asin(sin_alpha);			//em radianos
						}	
						double alpha_deg = 180 * alpha / PI;								//em graus
						//ID para acesso às curvas aerodinâmicas do perfil
						int aero_ID = db.sections[section - 1]->aerodynamicdataID;
						//Coeficientes aerodinâmicos de acordo com o ângulo de ataque
						CL = db.aerodynamic_data[aero_ID - 1]->CL->GetValueAt(alpha_deg, 0);	//Interpolação na curva aerodinâmica CL
						CD = db.aerodynamic_data[aero_ID - 1]->CD->GetValueAt(alpha_deg, 0);	//Interpolação na curva aerodinâmica CD
						CM = db.aerodynamic_data[aero_ID - 1]->CM->GetValueAt(alpha_deg, 0);	//Interpolação na curva aerodinâmica CM
					}	
					///////////////////////////////////Se é definido o BEM/////////////////////////////////////////////
					//***********************************************************************************************//
					///////////////////////////////////////////////////////////////////////////////////////////////////
					else
					{
						//Correção do u_inf - eliminando skew - tomando somente componente na direção axial do rotor
						//Obs: plano do rotor é o yz do sistema de coordenadas local definido para o rotor no BEM
						u_inf = dot(u_inf, *db.CS[db.bem->CS_rotor - 1]->E1)*(*db.CS[db.bem->CS_rotor - 1]->E1);
						//Posição do rotor
						Matrix x_Rotor(3);
						x_Rotor(0, 0) = db.nodes[db.bem->node_rotor - 1]->copy_coordinates[0];
						x_Rotor(1, 0) = db.nodes[db.bem->node_rotor - 1]->copy_coordinates[1];
						x_Rotor(2, 0) = db.nodes[db.bem->node_rotor - 1]->copy_coordinates[2];
						//Raio - projeção no plano do rotor
						Matrix r_proj = (x_Gauss - x_Rotor) - dot(x_Gauss - x_Rotor, *db.CS[db.bem->CS_rotor - 1]->E1)*(*db.CS[db.bem->CS_rotor - 1]->E1);
						//Raio do ponto de Gauss
						double radius_Gauss = norm(r_proj);
						//Direção radial
						Matrix versor_rad = (1.0 / radius_Gauss)*(r_proj);
						//Direção circunferencial
						Matrix versor_cir = cross(versor_rad, *db.CS[db.bem->CS_rotor - 1]->E1);
						//lambda_r
						double lambda_r;
						if (norm(u_inf) != 0.0)
							lambda_r = abs(dot(u_gauss, versor_cir)) / norm(u_inf);
						else
							lambda_r = 0.0;
						//sigma_l
						double sigma_l = db.bem->B*chord / (2 * PI*radius_Gauss);
						//beta - inclinação da linha de corda em relação ao plano do rotor
						double cos_gamma = dot(e1, *db.CS[db.bem->CS_rotor - 1]->E1);
						//Verificações de erros de arredondamento para evitar erros na função acos (argumento fora do range esperado)
						if (cos_gamma > 1.0)
							cos_gamma = 1.0;
						if (cos_gamma < -1.0)
							cos_gamma = -1.0;
						double gamma = abs(acos(cos_gamma));
						double beta = PI/2 - gamma;
						//Procedimento iterativo para cálculo de a e a'
						bool conv = false;
						double a, a_l, CT, F_tip, F_hub, F, phi;

						//Estimativas iniciais dos fatores a_l e a:
						//Cálculo do a_l
						a_l = 0.0;
						//Cálculo do a
						a = 0.25*(2.0 + PI*lambda_r*sigma_l - sqrt(4.0 - 4.0*PI*lambda_r*sigma_l + lambda_r*sigma_l*(8.0*beta + PI*sigma_l)));
						double prev_a, prev_a_l;	//cópias para critério de parada
						while (conv == false)
						{
							//Velocidade circunferencial - já está invertendo o sinal!
							Matrix u_circ = dot(-1.0*u_gauss, versor_cir)*(1+a_l)*versor_cir;
							//Velocidade relativa total
							u_rel = u_circ + u_inf*(1 - a);
							//Cálculo do ângulo de ataque - calculado utilizando-se a orientação atual da ST no início do time step, direção do vento ambiente incidente e velocidade
							alpha = 0.0;
							if (norm(u_rel) != 0.0)
							{
								double sin_alpha = dot(cross(e1, u_rel), e3) / norm(u_rel);
								//Verificações de erros de arredondamento para evitar erros na função asin (argumento fora do range esperado)
								if (sin_alpha > 1.0)
									sin_alpha = 1.0;
								if (sin_alpha < -1.0)
									sin_alpha = -1.0;
								alpha = asin(sin_alpha);			//em radianos
							}
							double alpha_deg = 180 * alpha / PI;								//em graus
							//ID para acesso às curvas aerodinâmicas do perfil
							int aero_ID = db.sections[section - 1]->aerodynamicdataID;
							//Coeficientes aerodinâmicos de acordo com o ângulo de ataque
							CL = db.aerodynamic_data[aero_ID - 1]->CL->GetValueAt(alpha_deg, 0);	//Interpolação na curva aerodinâmica CL
							CD = db.aerodynamic_data[aero_ID - 1]->CD->GetValueAt(alpha_deg, 0);	//Interpolação na curva aerodinâmica CD
							CM = db.aerodynamic_data[aero_ID - 1]->CM->GetValueAt(alpha_deg, 0);	//Interpolação na curva aerodinâmica CM
							//Cálculo do phi - ângulo de entrada do escoamento relativo no plano do rotor (já corrigido pelo a e a')
							phi = atan2(norm(u_inf)*(1 - a), (norm(u_circ)*(1 + a_l)));
							//CT
							if (sin(phi) != 0.0)
								CT = sigma_l*(1 - a)*(1 - a)*(CL*cos(phi) + CD*sin(phi)) / (sin(phi)*sin(phi));
							else
								CT = 0.0;
							//Ftip
							F_tip = (2.0 / PI)*acos(exp(-db.bem->B*(db.bem->R - radius_Gauss) / (2.0 * radius_Gauss*abs(sin(phi)))));
							F_hub = (2.0 / PI)*acos(exp(-db.bem->B*(radius_Gauss - db.bem->Rhub) / (2.0 * db.bem->Rhub*abs(sin(phi)))));
							F = F_hub*F_tip;
							//F = 1.0;
							//Salva cópias das anteriores
							prev_a = a;
							prev_a_l = a_l;
							//Glauert
							if (CT > 0.96*F)
								a = (18.0 * F - 20.0 - 3.0 * sqrt(CT*(50.0 - 36.0 * F) + 12.0 * F*(3.0 * F - 4.0))) / (36.0 * F - 50.0);
							//Standard BEM
							else
								a = 1.0 / (1.0 + 4.0*F*sin(phi)*sin(phi) / (sigma_l*(CL*cos(phi) + CD*sin(phi))));
							//Cálculo do a_l
							a_l = 1.0 / (-1.0 + 4.0*F*sin(phi)*cos(phi) / (sigma_l*(CL*sin(phi) - CD*cos(phi))));
							//Critério de parada
							if (abs((prev_a - a) < db.bem->tol_bem && abs((prev_a_l - a_l) < db.bem->tol_bem)))
								conv = true;
						}
					}
					
					//Cálculo dos esforços aerodinâmicos
					double FL = 0.5*rho_air*norm(u_rel)*norm(u_rel)*chord*CL;			//Força de sustentação
					double FD = 0.5*rho_air*norm(u_rel)*norm(u_rel)*chord*CD;			//Força de arrasto
					double MA = 0.5*rho_air*norm(u_rel)*norm(u_rel)*chord*chord*CM;		//Momento aerodinâmico
					//Decomposição dos esforços nas direções para compor os carregamentos
					drag = FD*cos(alpha)*e1 + FD*sin(alpha)*e2;
					lift = -FL*sin(alpha)*e1 + FL*cos(alpha)*e2;
					//AC_O é o vetor que vai de O (eixo da barra) até a posição do centro aerodinâmico - é utilizado para o transporte das forças de arrasto e sustentação para o eixo e calcular o binário de transporte
					Matrix AC_O = db.sections[section - 1]->AC(0, 0)*e1 + db.sections[section - 1]->AC(1, 0)*e2;
					moment = -MA*e3 + cross(AC_O, drag + lift);

					//Esforços nodais equivalentes - forças aerodinâmicas
					(*e_loading)(0, 0) += l_factor*load_multiplier*alpha1*jacobian* N1[gauss] * (lift(0, 0) + drag(0, 0));
					(*e_loading)(1, 0) += l_factor*load_multiplier*alpha1*jacobian* N1[gauss] * (lift(1, 0) + drag(1, 0));
					(*e_loading)(2, 0) += l_factor*load_multiplier*alpha1*jacobian* N1[gauss] * (lift(2, 0) + drag(2, 0));

					(*e_loading)(6, 0) += l_factor*load_multiplier*alpha1*jacobian* N2[gauss] * (lift(0, 0) + drag(0, 0));
					(*e_loading)(7, 0) += l_factor*load_multiplier*alpha1*jacobian* N2[gauss] * (lift(1, 0) + drag(1, 0));
					(*e_loading)(8, 0) += l_factor*load_multiplier*alpha1*jacobian* N2[gauss] * (lift(2, 0) + drag(2, 0));

					(*e_loading)(12, 0) += l_factor*load_multiplier*alpha1*jacobian* N3[gauss] * (lift(0, 0) + drag(0, 0));
					(*e_loading)(13, 0) += l_factor*load_multiplier*alpha1*jacobian* N3[gauss] * (lift(1, 0) + drag(1, 0));
					(*e_loading)(14, 0) += l_factor*load_multiplier*alpha1*jacobian* N3[gauss] * (lift(2, 0) + drag(2, 0));

					//Esforços nodais equivalentes - momentos aerodinâmicos
					alpha_e = norm(alpha_d);							//Valor escalar do parametro alpha
					A_d = skew(alpha_d);								//Matriz A_d
					g = 4.0 / (4.0 + alpha_e*alpha_e);					//função g(alpha)
					Q_d = *I3 + g*(A_d + 0.5*A_d*A_d);					//Tensor de rotação
					Xi_d = g*(*I3 + 0.5*A_d);
					Matrix pseudo_moment = transp(Xi_d)*moment;

					//Reporting aerodynamic data
					//data.myprintf("Element %d  Phi  %.3f  CL  %.3f  CD  %.3f  CM  %.3f  FL  %.3e  FD  %.3e  MA  %.3e\n",number,alpha_deg,CL,CD,CM,FL,FD,MA);
					//data.myprintf("Uinf  %.3f  %.3f  %.3f   Urel  %.3f  %.3f  %.3f \n", u_inf(0, 0), u_inf(1, 0), u_inf(2, 0), u_rel(0, 0), u_rel(1, 0), u_rel(2, 0));
					//data.myprintf("e1  %.3f  %.3f  %.3f  e2  %.3f  %.3f  %.3f   e3  %.3f  %.3f  %.3f \n", e1(0, 0), e1(1, 0), e1(2, 0), e2(0, 0), e2(1, 0), e2(2, 0), e3(0, 0), e3(1, 0), e3(2, 0));
					
					//Esforços nodais equivalentes - momentos aerodinâmicos
					(*e_loading)(3, 0) += l_factor*load_multiplier*alpha1*jacobian* N1[gauss] * pseudo_moment(0, 0);
					(*e_loading)(4, 0) += l_factor*load_multiplier*alpha1*jacobian* N1[gauss] * pseudo_moment(1, 0);
					(*e_loading)(5, 0) += l_factor*load_multiplier*alpha1*jacobian* N1[gauss] * pseudo_moment(2, 0);

					(*e_loading)(9, 0)  += l_factor*load_multiplier*alpha1*jacobian* N2[gauss] * pseudo_moment(0, 0);
					(*e_loading)(10, 0) += l_factor*load_multiplier*alpha1*jacobian* N2[gauss] * pseudo_moment(1, 0);
					(*e_loading)(11, 0) += l_factor*load_multiplier*alpha1*jacobian* N2[gauss] * pseudo_moment(2, 0);

					(*e_loading)(15, 0) += l_factor*load_multiplier*alpha1*jacobian* N3[gauss] * pseudo_moment(0, 0);
					(*e_loading)(16, 0) += l_factor*load_multiplier*alpha1*jacobian* N3[gauss] * pseudo_moment(1, 0);
					(*e_loading)(17, 0) += l_factor*load_multiplier*alpha1*jacobian* N3[gauss] * pseudo_moment(2, 0);

					//Contribuição do operador tangente - somente por conta do pseudo-momento
					Matrix stiff_moment = V(alpha_d, moment, alpha_e);	//contribuição não nula - rigidez devido aos momentos
					Matrix stiff_load(6, 6);
					//Compondo operador stiff_load (somente contribuições devido à rigidez de momentos)
					for (int i = 0; i < 3; i++)
					for (int j = 0; j < 3; j++)
						stiff_load(i + 3, j + 3) = stiff_moment(i, j);
					(*loading_stiffness) = (*loading_stiffness) + (l_factor*load_multiplier*alpha1*jacobian)*transp(*N[gauss])*stiff_load*(*N[gauss]);
				}//Gauss
			}//if wind_data_exist
		}
		//Se existe carregamentos hidrodinâmicos
		if (db.environment->ocean_data_exist == true)
		{
			//Inserir aqui os carregamentos hidrodinâmicos
		}
	}
}

//Monta matriz de transformação de coordenadas
void Beam_1::TransformMatrix()
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
void Beam_1::MountGlobal()
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
			GL_global_1 = db.nodes[nodes[1] - 1]->GLs[i-6];
		if (i > 11)
			GL_global_1 = db.nodes[nodes[2] - 1]->GLs[i-12];

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
void Beam_1::SaveLagrange()
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

//Comprimento do elemento
double Beam_1::CalculateLength()
{
	double x1 = db.nodes[nodes[0] - 1]->ref_coordinates[0];
	double y1 = db.nodes[nodes[0] - 1]->ref_coordinates[1];
	double z1 = db.nodes[nodes[0] - 1]->ref_coordinates[2];
	double x2 = db.nodes[nodes[2] - 1]->ref_coordinates[0];
	double y2 = db.nodes[nodes[2] - 1]->ref_coordinates[1];
	double z2 = db.nodes[nodes[2] - 1]->ref_coordinates[2];
	return sqrt((x2 - x1)*(x2 - x1) + (y2 - y1)*(y2 - y1) + (z2 - z1)*(z2 - z1));
}
void Beam_1::Zeros()
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
	zeros(loading_stiffness);
	kinetic_energy = 0.0;
	strain_energy = 0.0;
	potential_gravitational_energy = 0.0;
}

//Monta a matriz de massa
void Beam_1::MountMassModal()
{
	zeros(mass_modal);
	zeros(damping_modal);
	//Loop nos pontos de Gauss
	for (int gauss = 0; gauss < 2; gauss++)
	{
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
void Beam_1::MountDampingModal()
{
	zeros(mass_modal);
	zeros(damping_modal);
	//TODO
}

//Monta a matriz de massa
void Beam_1::MountMass()
{
	if (typeid(*db.solution[db.current_solution_number - 1]) == typeid(Dynamic))
	{
		Dynamic* ptr_sol = static_cast<Dynamic*>(db.solution[db.current_solution_number - 1]);

		Matrix omega_i(3);
		Matrix domega_i(3);
		Matrix du_i(3);
		Matrix ddu_i(3);
		//Loop nos pontos de Gauss
		for (int gauss = 0; gauss < 2; gauss++)
		{
			//Calculando as variáveis necessárias para obter o valor do integrando							
			(omega_i)(0, 0) = (db.nodes[nodes[0] - 1]->copy_vel[3])* N1[gauss] +
				(db.nodes[nodes[1] - 1]->copy_vel[3])* N2[gauss] +
				(db.nodes[nodes[2] - 1]->copy_vel[3])* N3[gauss];
			(omega_i)(1, 0) = (db.nodes[nodes[0] - 1]->copy_vel[4])* N1[gauss] +
				(db.nodes[nodes[1] - 1]->copy_vel[4])* N2[gauss] +
				(db.nodes[nodes[2] - 1]->copy_vel[4])* N3[gauss];
			(omega_i)(2, 0) = (db.nodes[nodes[0] - 1]->copy_vel[5])* N1[gauss] +
				(db.nodes[nodes[1] - 1]->copy_vel[5])* N2[gauss] +
				(db.nodes[nodes[2] - 1]->copy_vel[5])* N3[gauss];

			(domega_i)(0, 0) = (db.nodes[nodes[0] - 1]->copy_accel[3])* N1[gauss] +
				(db.nodes[nodes[1] - 1]->copy_accel[3])* N2[gauss] +
				(db.nodes[nodes[2] - 1]->copy_accel[3])* N3[gauss];
			(domega_i)(1, 0) = (db.nodes[nodes[0] - 1]->copy_accel[4])* N1[gauss] +
				(db.nodes[nodes[1] - 1]->copy_accel[4])* N2[gauss] +
				(db.nodes[nodes[2] - 1]->copy_accel[4])* N3[gauss];
			(domega_i)(2, 0) = (db.nodes[nodes[0] - 1]->copy_accel[5])* N1[gauss] +
				(db.nodes[nodes[1] - 1]->copy_accel[5])* N2[gauss] +
				(db.nodes[nodes[2] - 1]->copy_accel[5])* N3[gauss];

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

			//Para o sistema do elemento
			omega_i = (*transform3)*omega_i;
			domega_i = (*transform3)*domega_i;
			du_i = (*transform3)*du_i;
			ddu_i = (*transform3)*ddu_i;
			
			//Cálculo do integrando no ponto de Gauss
			EvaluateInertialContributions(temp_v, &ptr_sol->a1, &ptr_sol->a2, &ptr_sol->a3, &ptr_sol->a4,
				&ptr_sol->a5, &ptr_sol->a6, lag_save->alpha_i[gauss]->getMatrix(), alpha_delta[gauss]->getMatrix(),
				lag_save->u_i[gauss]->getMatrix(), u_delta[gauss]->getMatrix(), omega_i.getMatrix(), domega_i.getMatrix(), du_i.getMatrix(), ddu_i.getMatrix(), pJr, pMr, br->getMatrix(), dT->getMatrix(), pDdT);
			//Transformando o operador tangente em Matrix
			DdT->PtrToMatrix(pDdT, 6);
			(*inertial_loading) = (*inertial_loading) + (alpha1*jacobian)*transp(*N[gauss])*(*dT);
			(*mass) = (*mass) + (alpha1*jacobian)*transp(*N[gauss])*(*DdT)*(*N[gauss]);
			kinetic_energy += (*tempkin) * (alpha1*jacobian);
		}//end of gauss
		//Conversão para o sistema global
		(*inertial_loading) = transp(*transform)*(*inertial_loading);
		(*mass) = (transp(*transform)*(*mass))*(*transform);
	}
}

//Monta a matriz de amortecimento
void Beam_1::MountDamping(bool update_rayleigh)
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
void Beam_1::MountDyn()
{
	(*P_loading) = (*P_loading) + (*inertial_loading) + (*damping_loading);
	(*i_loading) = (*i_loading) + (*inertial_loading) + (*damping_loading);
	(*stiffness) = (*stiffness) + (*mass) + (*damping);
}

//Montagens para análise modal - inserção da matriz de massa e amortecimento na matriz de rigidez para posterior montagem global
void Beam_1::MountDynModal()
{
	(*stiffness) = (*mass_modal) + (*damping_modal);
}

void Beam_1::EvaluateMassModal(double* v, double* alphai, double** Jr, double** Mr, double* br, double** matrixm)
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
	matrixm[0][3] = -v[125];
	matrixm[0][4] = -v[126];
	matrixm[0][5] = -v[127];
	matrixm[1][0] = v[92];
	matrixm[1][1] = v[96];
	matrixm[1][2] = v[97];
	matrixm[1][3] = -v[128];
	matrixm[1][4] = -v[129];
	matrixm[1][5] = -v[130];
	matrixm[2][0] = v[98];
	matrixm[2][1] = v[102];
	matrixm[2][2] = v[103];
	matrixm[2][3] = -v[131];
	matrixm[2][4] = -v[132];
	matrixm[2][5] = -v[133];
	matrixm[3][0] = v[125];
	matrixm[3][1] = v[126];
	matrixm[3][2] = v[127];
	matrixm[3][3] = v[105] * v[73] + v[106] * v[74] + v[107] * v[75];
	matrixm[3][4] = v[105] * v[77] + v[106] * v[79] + v[107] * v[80];
	matrixm[3][5] = v[105] * v[82] + v[106] * v[84] + v[107] * v[85];
	matrixm[4][0] = v[128];
	matrixm[4][1] = v[129];
	matrixm[4][2] = v[130];
	matrixm[4][3] = v[111] * v[73] + v[112] * v[74] + v[113] * v[75];
	matrixm[4][4] = v[111] * v[77] + v[112] * v[79] + v[113] * v[80];
	matrixm[4][5] = v[111] * v[82] + v[112] * v[84] + v[113] * v[85];
	matrixm[5][0] = v[131];
	matrixm[5][1] = v[132];
	matrixm[5][2] = v[133];
	matrixm[5][3] = v[117] * v[73] + v[118] * v[74] + v[119] * v[75];
	matrixm[5][4] = v[117] * v[77] + v[118] * v[79] + v[119] * v[80];
	matrixm[5][5] = v[117] * v[82] + v[118] * v[84] + v[119] * v[85];
};


//Calcula as contribuições inerciais para a análise dinâmica - inclui todas as contribuições para a forma fraca e para o operador tangente
void Beam_1::EvaluateInertialContributions(double* v, double(*a1)
	, double(*a2), double(*a3), double(*a4), double(*a5), double(*a6)
	, double* alphai, double* alphad, double* ui, double* ud, double* omegai
	, double* domegai, double* dui, double* ddui, double** Jr, double** Mr
	, double* br, double* dT, double** DdT)
{
	v[2] = (*a1);
	v[3] = (*a2);
	v[4] = (*a3);
	v[5] = (*a4);
	v[6] = (*a5);
	v[7] = (*a6);
	v[8] = alphai[0];
	v[128] = (v[8] * v[8]);
	v[9] = alphai[1];
	v[605] = 0.5e0*v[9];
	v[126] = v[605] * v[8];
	v[121] = (v[9] * v[9]);
	v[10] = alphai[2];
	v[133] = v[10] * v[605];
	v[131] = 0.5e0*v[10] * v[8];
	v[122] = (v[10] * v[10]);
	v[613] = v[121] + v[122];
	v[11] = alphad[0];
	v[247] = 0.5e0*v[11];
	v[245] = 2e0*v[11];
	v[112] = (v[11] * v[11]);
	v[12] = alphad[1];
	v[248] = 2e0*v[12];
	v[246] = 0.5e0*v[12];
	v[110] = 1e0*v[11] * v[246];
	v[105] = (v[12] * v[12]);
	v[293] = -v[105] - v[112];
	v[13] = alphad[2];
	v[271] = v[110] + v[13];
	v[261] = v[110] - v[13];
	v[250] = 2e0*v[13];
	v[249] = 0.5e0*v[13];
	v[606] = 1e0*v[249];
	v[117] = v[12] * v[606];
	v[288] = v[11] + v[117];
	v[280] = -v[11] + v[117];
	v[115] = v[11] * v[606];
	v[284] = v[115] - v[12];
	v[267] = v[115] + v[12];
	v[106] = (v[13] * v[13]);
	v[275] = -v[106] - v[112];
	v[256] = -v[105] - v[106];
	v[612] = 1e0*v[256];
	v[251] = 4e0 + v[105] + v[106] + v[112];
	v[607] = -4e0 / (v[251] * v[251]);
	v[255] = v[250] * v[607];
	v[611] = 0.5e0*v[255];
	v[297] = v[293] * v[611];
	v[254] = v[248] * v[607];
	v[608] = 0.5e0*v[254];
	v[459] = -(v[247] * v[254]);
	v[455] = -(v[249] * v[254]);
	v[277] = v[275] * v[608];
	v[252] = v[245] * v[607];
	v[610] = 0.5e0*v[252];
	v[454] = -(v[249] * v[252]);
	v[257] = v[256] * v[610];
	v[17] = ud[0];
	v[18] = ud[1];
	v[19] = ud[2];
	v[20] = omegai[0];
	v[21] = omegai[1];
	v[22] = omegai[2];
	v[23] = domegai[0];
	v[200] = v[11] * v[2] - v[20] * v[3] - v[23] * v[4];
	v[194] = v[11] * v[5] + v[20] * v[6] + v[23] * v[7];
	v[24] = domegai[1];
	v[201] = v[12] * v[2] - v[21] * v[3] - v[24] * v[4];
	v[195] = v[12] * v[5] + v[21] * v[6] + v[24] * v[7];
	v[25] = domegai[2];
	v[202] = v[13] * v[2] - v[22] * v[3] - v[25] * v[4];
	v[196] = v[13] * v[5] + v[22] * v[6] + v[25] * v[7];
	v[26] = dui[0];
	v[27] = dui[1];
	v[28] = dui[2];
	v[29] = ddui[0];
	v[30] = ddui[1];
	v[31] = ddui[2];
	v[32] = Jr[0][0];
	v[33] = Jr[0][1];
	v[34] = Jr[0][2];
	v[35] = Jr[1][0];
	v[36] = Jr[1][1];
	v[37] = Jr[1][2];
	v[38] = Jr[2][0];
	v[39] = Jr[2][1];
	v[40] = Jr[2][2];
	v[41] = Mr[0][0];
	v[42] = Mr[0][1];
	v[43] = Mr[0][2];
	v[44] = Mr[1][0];
	v[45] = Mr[1][1];
	v[46] = Mr[1][2];
	v[47] = Mr[2][0];
	v[48] = Mr[2][1];
	v[49] = Mr[2][2];
	v[50] = br[0];
	v[51] = br[1];
	v[52] = br[2];
	v[104] = 4e0 / v[251];
	v[457] = -0.5e0*v[104];
	v[609] = 1e0*v[457];
	v[460] = v[246] * v[254] - v[457];
	v[458] = -(v[247] * v[252]) + v[457];
	v[456] = -(v[249] * v[255]) + v[457];
	v[295] = v[248] * v[609];
	v[296] = v[295] + v[293] * v[608];
	v[292] = v[245] * v[609];
	v[294] = v[292] + v[293] * v[610];
	v[289] = v[104] + v[252] * v[288];
	v[286] = -v[104] + v[254] * v[284];
	v[281] = -v[104] + v[252] * v[280];
	v[278] = 1e0*v[250] * v[609];
	v[279] = v[278] + v[275] * v[611];
	v[276] = v[292] + 1e0*v[275] * v[610];
	v[274] = v[104] + v[255] * v[271];
	v[269] = v[104] + v[254] * v[267];
	v[266] = -(v[104] * v[249]);
	v[290] = -v[266] + v[254] * v[288];
	v[285] = -v[266] + v[252] * v[284];
	v[282] = -v[266] + v[254] * v[280];
	v[268] = -v[266] + v[252] * v[267];
	v[265] = -v[104] + v[255] * v[261];
	v[263] = -(v[104] * v[247]);
	v[287] = -v[263] + v[255] * v[284];
	v[273] = -v[263] + v[254] * v[271];
	v[270] = -v[263] + v[255] * v[267];
	v[264] = v[254] * v[261] - v[263];
	v[260] = v[104] * v[246];
	v[291] = v[260] + v[255] * v[288];
	v[283] = v[260] + v[255] * v[280];
	v[272] = v[260] + v[252] * v[271];
	v[262] = v[260] + v[252] * v[261];
	v[259] = v[278] + v[611] * v[612];
	v[258] = v[295] + v[608] * v[612];
	v[107] = 1e0 - 1e0*v[457] * v[612];
	v[397] = v[107] * v[2] + v[200] * v[257] + v[201] * v[262] + v[202] * v[268];
	v[388] = v[194] * v[257] + v[195] * v[262] + v[196] * v[268] + v[107] * v[5];
	v[108] = v[104] * v[261];
	v[398] = v[108] * v[2] + v[200] * v[258] + v[201] * v[264] + v[202] * v[269];
	v[389] = v[194] * v[258] + v[195] * v[264] + v[196] * v[269] + v[108] * v[5];
	v[109] = v[104] * v[267];
	v[399] = v[109] * v[2] + v[200] * v[259] + v[201] * v[265] + v[202] * v[270];
	v[390] = v[194] * v[259] + v[195] * v[265] + v[196] * v[270] + v[109] * v[5];
	v[111] = v[104] * v[271];
	v[400] = v[111] * v[2] + v[200] * v[272] + v[201] * v[276] + v[202] * v[281];
	v[391] = v[194] * v[272] + v[195] * v[276] + v[196] * v[281] + v[111] * v[5];
	v[113] = 1e0 - 1e0*v[275] * v[609];
	v[401] = v[113] * v[2] + v[200] * v[273] + v[201] * v[277] + v[202] * v[282];
	v[392] = v[194] * v[273] + v[195] * v[277] + v[196] * v[282] + v[113] * v[5];
	v[114] = v[104] * v[280];
	v[402] = v[114] * v[2] + v[200] * v[274] + v[201] * v[279] + v[202] * v[283];
	v[393] = v[194] * v[274] + v[195] * v[279] + v[196] * v[283] + v[114] * v[5];
	v[116] = v[104] * v[284];
	v[403] = v[116] * v[2] + v[200] * v[285] + v[201] * v[289] + v[202] * v[294];
	v[394] = v[194] * v[285] + v[195] * v[289] + v[196] * v[294] + v[116] * v[5];
	v[118] = v[104] * v[288];
	v[404] = v[118] * v[2] + v[200] * v[286] + v[201] * v[290] + v[202] * v[296];
	v[395] = v[194] * v[286] + v[195] * v[290] + v[196] * v[296] + v[118] * v[5];
	v[119] = 1e0 - 1e0*v[293] * v[609];
	v[405] = v[119] * v[2] + v[200] * v[287] + v[201] * v[291] + v[202] * v[297];
	v[396] = v[194] * v[287] + v[195] * v[291] + v[196] * v[297] + v[119] * v[5];
	v[120] = 4e0 / (4e0 + v[128] + v[613]);
	v[614] = -0.5e0*v[120];
	v[123] = 1e0 + v[613] * v[614];
	v[124] = v[120] * (-v[10] + v[126]);
	v[125] = v[120] * (v[131] + v[9]);
	v[127] = v[120] * (v[10] + v[126]);
	v[129] = 1e0 + (v[122] + v[128])*v[614];
	v[130] = v[120] * (v[133] - v[8]);
	v[132] = v[120] * (v[131] - v[9]);
	v[336] = v[123] * v[287] + v[127] * v[291] + v[132] * v[297];
	v[335] = v[123] * v[286] + v[127] * v[290] + v[132] * v[296];
	v[334] = v[123] * v[285] + v[127] * v[289] + v[132] * v[294];
	v[318] = v[123] * v[274] + v[127] * v[279] + v[132] * v[283];
	v[317] = v[123] * v[273] + v[127] * v[277] + v[132] * v[282];
	v[316] = v[123] * v[272] + v[127] * v[276] + v[132] * v[281];
	v[300] = v[123] * v[259] + v[127] * v[265] + v[132] * v[270];
	v[299] = v[123] * v[258] + v[127] * v[264] + v[132] * v[269];
	v[298] = v[123] * v[257] + v[127] * v[262] + v[132] * v[268];
	v[134] = v[120] * (v[133] + v[8]);
	v[339] = v[124] * v[287] + v[129] * v[291] + v[134] * v[297];
	v[338] = v[124] * v[286] + v[129] * v[290] + v[134] * v[296];
	v[337] = v[124] * v[285] + v[129] * v[289] + v[134] * v[294];
	v[321] = v[124] * v[274] + v[129] * v[279] + v[134] * v[283];
	v[320] = v[124] * v[273] + v[129] * v[277] + v[134] * v[282];
	v[319] = v[124] * v[272] + v[129] * v[276] + v[134] * v[281];
	v[303] = v[124] * v[259] + v[129] * v[265] + v[134] * v[270];
	v[302] = v[124] * v[258] + v[129] * v[264] + v[134] * v[269];
	v[301] = v[124] * v[257] + v[129] * v[262] + v[134] * v[268];
	v[135] = 1e0 + 1e0*(v[121] + v[128])*v[614];
	v[342] = v[125] * v[287] + v[130] * v[291] + v[135] * v[297];
	v[487] = v[32] * v[336] + v[339] * v[35] + v[342] * v[38];
	v[484] = v[33] * v[336] + v[339] * v[36] + v[342] * v[39];
	v[481] = v[336] * v[34] + v[339] * v[37] + v[342] * v[40];
	v[387] = v[336] * v[50] + v[339] * v[51] + v[342] * v[52];
	v[351] = v[336] * v[41] + v[339] * v[44] + v[342] * v[47];
	v[348] = v[336] * v[42] + v[339] * v[45] + v[342] * v[48];
	v[345] = v[336] * v[43] + v[339] * v[46] + v[342] * v[49];
	v[341] = v[125] * v[286] + v[130] * v[290] + v[135] * v[296];
	v[486] = v[32] * v[335] + v[338] * v[35] + v[341] * v[38];
	v[483] = v[33] * v[335] + v[338] * v[36] + v[341] * v[39];
	v[480] = v[335] * v[34] + v[338] * v[37] + v[341] * v[40];
	v[386] = v[335] * v[50] + v[338] * v[51] + v[341] * v[52];
	v[350] = v[335] * v[41] + v[338] * v[44] + v[341] * v[47];
	v[347] = v[335] * v[42] + v[338] * v[45] + v[341] * v[48];
	v[344] = v[335] * v[43] + v[338] * v[46] + v[341] * v[49];
	v[340] = v[125] * v[285] + v[130] * v[289] + v[135] * v[294];
	v[485] = v[32] * v[334] + v[337] * v[35] + v[340] * v[38];
	v[482] = v[33] * v[334] + v[337] * v[36] + v[340] * v[39];
	v[479] = v[334] * v[34] + v[337] * v[37] + v[340] * v[40];
	v[385] = v[334] * v[50] + v[337] * v[51] + v[340] * v[52];
	v[349] = v[334] * v[41] + v[337] * v[44] + v[340] * v[47];
	v[346] = v[334] * v[42] + v[337] * v[45] + v[340] * v[48];
	v[343] = v[334] * v[43] + v[337] * v[46] + v[340] * v[49];
	v[324] = v[125] * v[274] + v[130] * v[279] + v[135] * v[283];
	v[478] = v[318] * v[32] + v[321] * v[35] + v[324] * v[38];
	v[475] = v[318] * v[33] + v[321] * v[36] + v[324] * v[39];
	v[472] = v[318] * v[34] + v[321] * v[37] + v[324] * v[40];
	v[384] = v[318] * v[50] + v[321] * v[51] + v[324] * v[52];
	v[333] = v[318] * v[41] + v[321] * v[44] + v[324] * v[47];
	v[330] = v[318] * v[42] + v[321] * v[45] + v[324] * v[48];
	v[327] = v[318] * v[43] + v[321] * v[46] + v[324] * v[49];
	v[323] = v[125] * v[273] + v[130] * v[277] + v[135] * v[282];
	v[477] = v[317] * v[32] + v[320] * v[35] + v[323] * v[38];
	v[474] = v[317] * v[33] + v[320] * v[36] + v[323] * v[39];
	v[471] = v[317] * v[34] + v[320] * v[37] + v[323] * v[40];
	v[383] = v[317] * v[50] + v[320] * v[51] + v[323] * v[52];
	v[332] = v[317] * v[41] + v[320] * v[44] + v[323] * v[47];
	v[329] = v[317] * v[42] + v[320] * v[45] + v[323] * v[48];
	v[326] = v[317] * v[43] + v[320] * v[46] + v[323] * v[49];
	v[322] = v[125] * v[272] + v[130] * v[276] + v[135] * v[281];
	v[476] = v[316] * v[32] + v[319] * v[35] + v[322] * v[38];
	v[473] = v[316] * v[33] + v[319] * v[36] + v[322] * v[39];
	v[470] = v[316] * v[34] + v[319] * v[37] + v[322] * v[40];
	v[382] = v[316] * v[50] + v[319] * v[51] + v[322] * v[52];
	v[331] = v[316] * v[41] + v[319] * v[44] + v[322] * v[47];
	v[328] = v[316] * v[42] + v[319] * v[45] + v[322] * v[48];
	v[325] = v[316] * v[43] + v[319] * v[46] + v[322] * v[49];
	v[306] = v[125] * v[259] + v[130] * v[265] + v[135] * v[270];
	v[469] = v[300] * v[32] + v[303] * v[35] + v[306] * v[38];
	v[466] = v[300] * v[33] + v[303] * v[36] + v[306] * v[39];
	v[463] = v[300] * v[34] + v[303] * v[37] + v[306] * v[40];
	v[381] = v[300] * v[50] + v[303] * v[51] + v[306] * v[52];
	v[315] = v[300] * v[41] + v[303] * v[44] + v[306] * v[47];
	v[312] = v[300] * v[42] + v[303] * v[45] + v[306] * v[48];
	v[309] = v[300] * v[43] + v[303] * v[46] + v[306] * v[49];
	v[305] = v[125] * v[258] + v[130] * v[264] + v[135] * v[269];
	v[468] = v[299] * v[32] + v[302] * v[35] + v[305] * v[38];
	v[465] = v[299] * v[33] + v[302] * v[36] + v[305] * v[39];
	v[462] = v[299] * v[34] + v[302] * v[37] + v[305] * v[40];
	v[380] = v[299] * v[50] + v[302] * v[51] + v[305] * v[52];
	v[314] = v[299] * v[41] + v[302] * v[44] + v[305] * v[47];
	v[311] = v[299] * v[42] + v[302] * v[45] + v[305] * v[48];
	v[308] = v[299] * v[43] + v[302] * v[46] + v[305] * v[49];
	v[304] = v[125] * v[257] + v[130] * v[262] + v[135] * v[268];
	v[467] = v[298] * v[32] + v[301] * v[35] + v[304] * v[38];
	v[464] = v[298] * v[33] + v[301] * v[36] + v[304] * v[39];
	v[461] = v[298] * v[34] + v[301] * v[37] + v[304] * v[40];
	v[379] = v[298] * v[50] + v[301] * v[51] + v[304] * v[52];
	v[313] = v[298] * v[41] + v[301] * v[44] + v[304] * v[47];
	v[310] = v[298] * v[42] + v[301] * v[45] + v[304] * v[48];
	v[307] = v[298] * v[43] + v[301] * v[46] + v[304] * v[49];
	v[136] = v[107] * v[123] + v[108] * v[127] + v[109] * v[132];
	v[137] = v[107] * v[124] + v[108] * v[129] + v[109] * v[134];
	v[138] = v[107] * v[125] + v[108] * v[130] + v[109] * v[135];
	v[169] = v[136] * v[34] + v[137] * v[37] + v[138] * v[40];
	v[168] = v[136] * v[33] + v[137] * v[36] + v[138] * v[39];
	v[167] = v[136] * v[32] + v[137] * v[35] + v[138] * v[38];
	v[526] = v[167] * v[300] + v[168] * v[303] + v[169] * v[306] + v[138] * v[463] + v[137] * v[466] + v[136] * v[469];
	v[525] = v[167] * v[299] + v[168] * v[302] + v[169] * v[305] + v[138] * v[462] + v[137] * v[465] + v[136] * v[468];
	v[524] = v[167] * v[298] + v[168] * v[301] + v[169] * v[304] + v[138] * v[461] + v[137] * v[464] + v[136] * v[467];
	v[151] = v[136] * v[43] + v[137] * v[46] + v[138] * v[49];
	v[150] = v[136] * v[42] + v[137] * v[45] + v[138] * v[48];
	v[149] = v[136] * v[41] + v[137] * v[44] + v[138] * v[47];
	v[354] = v[149] * v[300] + v[150] * v[303] + v[151] * v[306] + v[138] * v[309] + v[137] * v[312] + v[136] * v[315];
	v[353] = v[149] * v[299] + v[150] * v[302] + v[151] * v[305] + v[138] * v[308] + v[137] * v[311] + v[136] * v[314];
	v[352] = v[149] * v[298] + v[150] * v[301] + v[151] * v[304] + v[138] * v[307] + v[137] * v[310] + v[136] * v[313];
	v[139] = v[111] * v[123] + v[113] * v[127] + v[114] * v[132];
	v[140] = v[111] * v[124] + v[113] * v[129] + v[114] * v[134];
	v[141] = v[111] * v[125] + v[113] * v[130] + v[114] * v[135];
	v[529] = v[167] * v[318] + v[168] * v[321] + v[169] * v[324] + v[141] * v[463] + v[140] * v[466] + v[139] * v[469];
	v[528] = v[167] * v[317] + v[168] * v[320] + v[169] * v[323] + v[141] * v[462] + v[140] * v[465] + v[139] * v[468];
	v[527] = v[167] * v[316] + v[168] * v[319] + v[169] * v[322] + v[141] * v[461] + v[140] * v[464] + v[139] * v[467];
	v[357] = v[141] * v[309] + v[140] * v[312] + v[139] * v[315] + v[149] * v[318] + v[150] * v[321] + v[151] * v[324];
	v[356] = v[141] * v[308] + v[140] * v[311] + v[139] * v[314] + v[149] * v[317] + v[150] * v[320] + v[151] * v[323];
	v[355] = v[141] * v[307] + v[140] * v[310] + v[139] * v[313] + v[149] * v[316] + v[150] * v[319] + v[151] * v[322];
	v[175] = v[139] * v[34] + v[140] * v[37] + v[141] * v[40];
	v[174] = v[139] * v[33] + v[140] * v[36] + v[141] * v[39];
	v[173] = v[139] * v[32] + v[140] * v[35] + v[141] * v[38];
	v[538] = v[173] * v[318] + v[174] * v[321] + v[175] * v[324] + v[141] * v[472] + v[140] * v[475] + v[139] * v[478];
	v[537] = v[173] * v[317] + v[174] * v[320] + v[175] * v[323] + v[141] * v[471] + v[140] * v[474] + v[139] * v[477];
	v[536] = v[173] * v[316] + v[174] * v[319] + v[175] * v[322] + v[141] * v[470] + v[140] * v[473] + v[139] * v[476];
	v[535] = v[173] * v[300] + v[174] * v[303] + v[175] * v[306] + v[138] * v[472] + v[137] * v[475] + v[136] * v[478];
	v[534] = v[173] * v[299] + v[174] * v[302] + v[175] * v[305] + v[138] * v[471] + v[137] * v[474] + v[136] * v[477];
	v[533] = v[173] * v[298] + v[174] * v[301] + v[175] * v[304] + v[138] * v[470] + v[137] * v[473] + v[136] * v[476];
	v[157] = v[139] * v[43] + v[140] * v[46] + v[141] * v[49];
	v[156] = v[139] * v[42] + v[140] * v[45] + v[141] * v[48];
	v[155] = v[139] * v[41] + v[140] * v[44] + v[141] * v[47];
	v[366] = v[155] * v[318] + v[156] * v[321] + v[157] * v[324] + v[141] * v[327] + v[140] * v[330] + v[139] * v[333];
	v[365] = v[155] * v[317] + v[156] * v[320] + v[157] * v[323] + v[141] * v[326] + v[140] * v[329] + v[139] * v[332];
	v[364] = v[155] * v[316] + v[156] * v[319] + v[157] * v[322] + v[141] * v[325] + v[140] * v[328] + v[139] * v[331];
	v[363] = v[155] * v[300] + v[156] * v[303] + v[157] * v[306] + v[138] * v[327] + v[137] * v[330] + v[136] * v[333];
	v[362] = v[155] * v[299] + v[156] * v[302] + v[157] * v[305] + v[138] * v[326] + v[137] * v[329] + v[136] * v[332];
	v[361] = v[155] * v[298] + v[156] * v[301] + v[157] * v[304] + v[138] * v[325] + v[137] * v[328] + v[136] * v[331];
	v[142] = v[116] * v[123] + v[118] * v[127] + v[119] * v[132];
	v[143] = v[116] * v[124] + v[118] * v[129] + v[119] * v[134];
	v[144] = v[116] * v[125] + v[118] * v[130] + v[119] * v[135];
	v[541] = v[173] * v[336] + v[174] * v[339] + v[175] * v[342] + v[144] * v[472] + v[143] * v[475] + v[142] * v[478];
	v[540] = v[173] * v[335] + v[174] * v[338] + v[175] * v[341] + v[144] * v[471] + v[143] * v[474] + v[142] * v[477];
	v[539] = v[173] * v[334] + v[174] * v[337] + v[175] * v[340] + v[144] * v[470] + v[143] * v[473] + v[142] * v[476];
	v[532] = v[167] * v[336] + v[168] * v[339] + v[169] * v[342] + v[144] * v[463] + v[143] * v[466] + v[142] * v[469];
	v[531] = v[167] * v[335] + v[168] * v[338] + v[169] * v[341] + v[144] * v[462] + v[143] * v[465] + v[142] * v[468];
	v[530] = v[167] * v[334] + v[168] * v[337] + v[169] * v[340] + v[144] * v[461] + v[143] * v[464] + v[142] * v[467];
	v[369] = v[144] * v[327] + v[143] * v[330] + v[142] * v[333] + v[155] * v[336] + v[156] * v[339] + v[157] * v[342];
	v[368] = v[144] * v[326] + v[143] * v[329] + v[142] * v[332] + v[155] * v[335] + v[156] * v[338] + v[157] * v[341];
	v[367] = v[144] * v[325] + v[143] * v[328] + v[142] * v[331] + v[155] * v[334] + v[156] * v[337] + v[157] * v[340];
	v[360] = v[144] * v[309] + v[143] * v[312] + v[142] * v[315] + v[149] * v[336] + v[150] * v[339] + v[151] * v[342];
	v[359] = v[144] * v[308] + v[143] * v[311] + v[142] * v[314] + v[149] * v[335] + v[150] * v[338] + v[151] * v[341];
	v[358] = v[144] * v[307] + v[143] * v[310] + v[142] * v[313] + v[149] * v[334] + v[150] * v[337] + v[151] * v[340];
	v[181] = v[142] * v[34] + v[143] * v[37] + v[144] * v[40];
	v[180] = v[142] * v[33] + v[143] * v[36] + v[144] * v[39];
	v[179] = v[142] * v[32] + v[143] * v[35] + v[144] * v[38];
	v[550] = v[179] * v[336] + v[180] * v[339] + v[181] * v[342] + v[144] * v[481] + v[143] * v[484] + v[142] * v[487];
	v[549] = v[179] * v[335] + v[180] * v[338] + v[181] * v[341] + v[144] * v[480] + v[143] * v[483] + v[142] * v[486];
	v[548] = v[179] * v[334] + v[180] * v[337] + v[181] * v[340] + v[144] * v[479] + v[143] * v[482] + v[142] * v[485];
	v[547] = v[179] * v[318] + v[180] * v[321] + v[181] * v[324] + v[141] * v[481] + v[140] * v[484] + v[139] * v[487];
	v[546] = v[179] * v[317] + v[180] * v[320] + v[181] * v[323] + v[141] * v[480] + v[140] * v[483] + v[139] * v[486];
	v[545] = v[179] * v[316] + v[180] * v[319] + v[181] * v[322] + v[141] * v[479] + v[140] * v[482] + v[139] * v[485];
	v[544] = v[179] * v[300] + v[180] * v[303] + v[181] * v[306] + v[138] * v[481] + v[137] * v[484] + v[136] * v[487];
	v[543] = v[179] * v[299] + v[180] * v[302] + v[181] * v[305] + v[138] * v[480] + v[137] * v[483] + v[136] * v[486];
	v[542] = v[179] * v[298] + v[180] * v[301] + v[181] * v[304] + v[138] * v[479] + v[137] * v[482] + v[136] * v[485];
	v[163] = v[142] * v[43] + v[143] * v[46] + v[144] * v[49];
	v[162] = v[142] * v[42] + v[143] * v[45] + v[144] * v[48];
	v[161] = v[142] * v[41] + v[143] * v[44] + v[144] * v[47];
	v[378] = v[161] * v[336] + v[162] * v[339] + v[163] * v[342] + v[144] * v[345] + v[143] * v[348] + v[142] * v[351];
	v[377] = v[161] * v[335] + v[162] * v[338] + v[163] * v[341] + v[144] * v[344] + v[143] * v[347] + v[142] * v[350];
	v[376] = v[161] * v[334] + v[162] * v[337] + v[163] * v[340] + v[144] * v[343] + v[143] * v[346] + v[142] * v[349];
	v[375] = v[161] * v[318] + v[162] * v[321] + v[163] * v[324] + v[141] * v[345] + v[140] * v[348] + v[139] * v[351];
	v[374] = v[161] * v[317] + v[162] * v[320] + v[163] * v[323] + v[141] * v[344] + v[140] * v[347] + v[139] * v[350];
	v[373] = v[161] * v[316] + v[162] * v[319] + v[163] * v[322] + v[141] * v[343] + v[140] * v[346] + v[139] * v[349];
	v[372] = v[161] * v[300] + v[162] * v[303] + v[163] * v[306] + v[138] * v[345] + v[137] * v[348] + v[136] * v[351];
	v[371] = v[161] * v[299] + v[162] * v[302] + v[163] * v[305] + v[138] * v[344] + v[137] * v[347] + v[136] * v[350];
	v[370] = v[161] * v[298] + v[162] * v[301] + v[163] * v[304] + v[138] * v[343] + v[137] * v[346] + v[136] * v[349];
	v[148] = v[136] * v[149] + v[137] * v[150] + v[138] * v[151];
	v[492] = -(v[148] * v[455]);
	v[152] = v[139] * v[149] + v[140] * v[150] + v[141] * v[151];
	v[504] = -(v[152] * v[455]);
	v[153] = v[142] * v[149] + v[143] * v[150] + v[144] * v[151];
	v[516] = -(v[153] * v[455]);
	v[154] = v[136] * v[155] + v[137] * v[156] + v[138] * v[157];
	v[495] = v[154] * v[454];
	v[158] = v[139] * v[155] + v[140] * v[156] + v[141] * v[157];
	v[507] = v[158] * v[454];
	v[159] = v[142] * v[155] + v[143] * v[156] + v[144] * v[157];
	v[519] = v[159] * v[454];
	v[160] = v[136] * v[161] + v[137] * v[162] + v[138] * v[163];
	v[496] = -(v[160] * v[459]);
	v[451] = v[104] * v[160] + v[148] * v[260] + v[154] * v[263];
	v[445] = v[104] * v[154] - v[160] * v[263] + v[148] * v[266];
	v[439] = v[104] * v[148] - v[160] * v[260] - v[154] * v[266];
	v[164] = v[139] * v[161] + v[140] * v[162] + v[141] * v[163];
	v[508] = -(v[164] * v[459]);
	v[449] = v[104] * v[164] + v[152] * v[260] + v[158] * v[263];
	v[443] = v[104] * v[158] - v[164] * v[263] + v[152] * v[266];
	v[437] = v[104] * v[152] - v[164] * v[260] - v[158] * v[266];
	v[165] = v[142] * v[161] + v[143] * v[162] + v[144] * v[163];
	v[520] = -(v[165] * v[459]);
	v[448] = v[104] * v[165] + v[153] * v[260] + v[159] * v[263];
	v[442] = v[104] * v[159] - v[165] * v[263] + v[153] * v[266];
	v[436] = v[104] * v[153] - v[165] * v[260] - v[159] * v[266];
	v[166] = v[136] * v[167] + v[137] * v[168] + v[138] * v[169];
	v[170] = v[139] * v[167] + v[140] * v[168] + v[141] * v[169];
	v[171] = v[142] * v[167] + v[143] * v[168] + v[144] * v[169];
	v[172] = v[136] * v[173] + v[137] * v[174] + v[138] * v[175];
	v[176] = v[139] * v[173] + v[140] * v[174] + v[141] * v[175];
	v[177] = v[142] * v[173] + v[143] * v[174] + v[144] * v[175];
	v[178] = v[136] * v[179] + v[137] * v[180] + v[138] * v[181];
	v[182] = v[139] * v[179] + v[140] * v[180] + v[141] * v[181];
	v[183] = v[142] * v[179] + v[143] * v[180] + v[144] * v[181];
	v[184] = v[136] * v[50] + v[137] * v[51] + v[138] * v[52];
	v[434] = v[184] * v[2];
	v[185] = v[139] * v[50] + v[140] * v[51] + v[141] * v[52];
	v[433] = -(v[185] * v[2]);
	v[186] = v[142] * v[50] + v[143] * v[51] + v[144] * v[52];
	v[435] = -(v[186] * v[2]);
	v[187] = v[17] * v[5] + v[26] * v[6] + v[29] * v[7];
	v[188] = v[18] * v[5] + v[27] * v[6] + v[30] * v[7];
	v[189] = v[19] * v[5] + v[28] * v[6] + v[31] * v[7];
	v[233] = v[148] * v[187] + v[152] * v[188] + v[153] * v[189];
	v[232] = v[154] * v[187] + v[158] * v[188] + v[159] * v[189];
	v[231] = v[160] * v[187] + v[164] * v[188] + v[165] * v[189];
	v[190] = v[17] * v[2] - v[26] * v[3] - v[29] * v[4];
	v[191] = v[18] * v[2] - v[27] * v[3] - v[30] * v[4];
	v[553] = v[191] * v[381] - v[190] * v[384];
	v[552] = v[191] * v[380] - v[190] * v[383];
	v[551] = v[191] * v[379] - v[190] * v[382];
	v[218] = -(v[185] * v[190]) + v[184] * v[191];
	v[192] = v[19] * v[2] - v[28] * v[3] - v[31] * v[4];
	v[559] = -(v[192] * v[381]) + v[190] * v[387];
	v[558] = -(v[192] * v[380]) + v[190] * v[386];
	v[557] = -(v[192] * v[379]) + v[190] * v[385];
	v[556] = v[192] * v[384] - v[191] * v[387];
	v[555] = v[192] * v[383] - v[191] * v[386];
	v[554] = v[192] * v[382] - v[191] * v[385];
	v[220] = -(v[186] * v[191]) + v[185] * v[192];
	v[219] = v[186] * v[190] - v[184] * v[192];
	v[193] = v[107] * v[194] + v[108] * v[195] + v[109] * v[196];
	v[197] = v[111] * v[194] + v[113] * v[195] + v[114] * v[196];
	v[417] = -(v[197] * v[381]) + v[193] * v[384] + v[185] * v[390] - v[184] * v[393];
	v[416] = -(v[197] * v[380]) + v[193] * v[383] + v[185] * v[389] - v[184] * v[392];
	v[415] = -(v[197] * v[379]) + v[193] * v[382] + v[185] * v[388] - v[184] * v[391];
	v[198] = v[116] * v[194] + v[118] * v[195] + v[119] * v[196];
	v[580] = v[178] * v[390] + v[182] * v[393] + v[183] * v[396] + v[193] * v[544] + v[197] * v[547] + v[198] * v[550];
	v[579] = v[178] * v[389] + v[182] * v[392] + v[183] * v[395] + v[193] * v[543] + v[197] * v[546] + v[198] * v[549];
	v[578] = v[178] * v[388] + v[182] * v[391] + v[183] * v[394] + v[193] * v[542] + v[197] * v[545] + v[198] * v[548];
	v[574] = v[172] * v[390] + v[176] * v[393] + v[177] * v[396] + v[193] * v[535] + v[197] * v[538] + v[198] * v[541];
	v[573] = v[172] * v[389] + v[176] * v[392] + v[177] * v[395] + v[193] * v[534] + v[197] * v[537] + v[198] * v[540];
	v[572] = v[172] * v[388] + v[176] * v[391] + v[177] * v[394] + v[193] * v[533] + v[197] * v[536] + v[198] * v[539];
	v[571] = v[166] * v[390] + v[170] * v[393] + v[171] * v[396] + v[193] * v[526] + v[197] * v[529] + v[198] * v[532];
	v[570] = v[166] * v[389] + v[170] * v[392] + v[171] * v[395] + v[193] * v[525] + v[197] * v[528] + v[198] * v[531];
	v[569] = v[166] * v[388] + v[170] * v[391] + v[171] * v[394] + v[193] * v[524] + v[197] * v[527] + v[198] * v[530];
	v[411] = -(v[198] * v[384]) + v[197] * v[387] + v[186] * v[393] - v[185] * v[396];
	v[410] = -(v[198] * v[383]) + v[197] * v[386] + v[186] * v[392] - v[185] * v[395];
	v[409] = -(v[198] * v[382]) + v[197] * v[385] + v[186] * v[391] - v[185] * v[394];
	v[408] = v[198] * v[381] - v[193] * v[387] - v[186] * v[390] + v[184] * v[396];
	v[407] = v[198] * v[380] - v[193] * v[386] - v[186] * v[389] + v[184] * v[395];
	v[406] = v[198] * v[379] - v[193] * v[385] - v[186] * v[388] + v[184] * v[394];
	v[199] = v[107] * v[200] + v[108] * v[201] + v[109] * v[202];
	v[203] = v[111] * v[200] + v[113] * v[201] + v[114] * v[202];
	v[204] = v[116] * v[200] + v[118] * v[201] + v[119] * v[202];
	v[226] = v[178] * v[199] + v[182] * v[203] + v[183] * v[204];
	v[592] = v[226] * v[459];
	v[225] = v[172] * v[199] + v[176] * v[203] + v[177] * v[204];
	v[598] = -(v[225] * v[454]);
	v[224] = v[166] * v[199] + v[170] * v[203] + v[171] * v[204];
	v[600] = v[224] * v[455];
	v[205] = -(v[186] * v[193]) + v[184] * v[198];
	v[206] = v[186] * v[197] - v[185] * v[198];
	v[414] = -(v[203] * v[381]) + v[199] * v[384] + v[205] * v[390] - v[206] * v[393] + v[185] * v[399] - v[184] * v[402]
		+ v[193] * v[408] - v[197] * v[411];
	v[413] = -(v[203] * v[380]) + v[199] * v[383] + v[205] * v[389] - v[206] * v[392] + v[185] * v[398] - v[184] * v[401]
		+ v[193] * v[407] - v[197] * v[410];
	v[412] = -(v[203] * v[379]) + v[199] * v[382] + v[205] * v[388] - v[206] * v[391] + v[185] * v[397] - v[184] * v[400]
		+ v[193] * v[406] - v[197] * v[409];
	v[209] = v[192] + v[185] * v[199] - v[184] * v[203] + v[193] * v[205] - v[197] * v[206];
	v[207] = v[185] * v[193] - v[184] * v[197];
	v[423] = v[204] * v[381] - v[199] * v[387] - v[207] * v[390] + v[206] * v[396] - v[186] * v[399] + v[184] * v[405]
		+ v[198] * v[411] - v[193] * v[417];
	v[422] = v[204] * v[380] - v[199] * v[386] - v[207] * v[389] + v[206] * v[395] - v[186] * v[398] + v[184] * v[404]
		+ v[198] * v[410] - v[193] * v[416];
	v[421] = v[204] * v[379] - v[199] * v[385] - v[207] * v[388] + v[206] * v[394] - v[186] * v[397] + v[184] * v[403]
		+ v[198] * v[409] - v[193] * v[415];
	v[420] = -(v[204] * v[384]) + v[203] * v[387] + v[207] * v[393] - v[205] * v[396] + v[186] * v[402] - v[185] * v[405]
		- v[198] * v[408] + v[197] * v[417];
	v[419] = -(v[204] * v[383]) + v[203] * v[386] + v[207] * v[392] - v[205] * v[395] + v[186] * v[401] - v[185] * v[404]
		- v[198] * v[407] + v[197] * v[416];
	v[418] = -(v[204] * v[382]) + v[203] * v[385] + v[207] * v[391] - v[205] * v[394] + v[186] * v[400] - v[185] * v[403]
		- v[198] * v[406] + v[197] * v[415];
	v[211] = v[190] + v[186] * v[203] - v[185] * v[204] - v[198] * v[205] + v[197] * v[207];
	v[210] = v[191] - v[186] * v[199] + v[184] * v[204] + v[198] * v[206] - v[193] * v[207];
	v[214] = v[166] * v[193] + v[170] * v[197] + v[171] * v[198];
	v[215] = v[172] * v[193] + v[176] * v[197] + v[177] * v[198];
	v[624] = v[215] * v[390] - v[214] * v[393] + v[178] * v[399] + v[182] * v[402] + v[183] * v[405] + v[199] * v[544]
		+ v[203] * v[547] + v[204] * v[550] - v[197] * v[571] + v[193] * v[574];
	v[621] = v[215] * v[389] - v[214] * v[392] + v[178] * v[398] + v[182] * v[401] + v[183] * v[404] + v[199] * v[543]
		+ v[203] * v[546] + v[204] * v[549] - v[197] * v[570] + v[193] * v[573];
	v[618] = v[215] * v[388] - v[214] * v[391] + v[178] * v[397] + v[182] * v[400] + v[183] * v[403] + v[199] * v[542]
		+ v[203] * v[545] + v[204] * v[548] - v[197] * v[569] + v[193] * v[572];
	v[221] = -(v[197] * v[214]) + v[193] * v[215];
	v[616] = v[221] + v[226];
	v[591] = v[221] * v[459];
	v[216] = v[178] * v[193] + v[182] * v[197] + v[183] * v[198];
	v[626] = -(v[216] * v[390]) + v[214] * v[396] + v[172] * v[399] + v[176] * v[402] + v[177] * v[405] + v[199] * v[535]
		+ v[203] * v[538] + v[204] * v[541] + v[198] * v[571] - v[193] * v[580];
	v[623] = -(v[216] * v[389]) + v[214] * v[395] + v[172] * v[398] + v[176] * v[401] + v[177] * v[404] + v[199] * v[534]
		+ v[203] * v[537] + v[204] * v[540] + v[198] * v[570] - v[193] * v[579];
	v[620] = -(v[216] * v[388]) + v[214] * v[394] + v[172] * v[397] + v[176] * v[400] + v[177] * v[403] + v[199] * v[533]
		+ v[203] * v[536] + v[204] * v[539] + v[198] * v[569] - v[193] * v[578];
	v[625] = v[216] * v[393] - v[215] * v[396] + v[166] * v[399] + v[170] * v[402] + v[171] * v[405] + v[199] * v[526]
		+ v[203] * v[529] + v[204] * v[532] - v[198] * v[574] + v[197] * v[580];
	v[622] = v[216] * v[392] - v[215] * v[395] + v[166] * v[398] + v[170] * v[401] + v[171] * v[404] + v[199] * v[525]
		+ v[203] * v[528] + v[204] * v[531] - v[198] * v[573] + v[197] * v[579];
	v[619] = v[216] * v[391] - v[215] * v[394] + v[166] * v[397] + v[170] * v[400] + v[171] * v[403] + v[199] * v[524]
		+ v[203] * v[527] + v[204] * v[530] - v[198] * v[572] + v[197] * v[578];
	v[223] = -(v[198] * v[215]) + v[197] * v[216];
	v[615] = v[223] + v[224];
	v[599] = v[223] * v[455];
	v[222] = v[198] * v[214] - v[193] * v[216];
	v[617] = v[222] + v[225];
	v[597] = -(v[222] * v[454]);
	dT[0] = v[153] * v[209] + v[152] * v[210] + v[148] * v[211];
	dT[1] = v[159] * v[209] + v[158] * v[210] + v[154] * v[211];
	dT[2] = v[165] * v[209] + v[164] * v[210] + v[160] * v[211];
	dT[3] = v[218] * v[436] + v[219] * v[437] + v[220] * v[439] + v[104] * v[615] - v[260] * v[616] - v[266] * v[617];
	dT[4] = v[218] * v[442] + v[219] * v[443] + v[220] * v[445] + v[266] * v[615] - v[263] * v[616] + v[104] * v[617];
	dT[5] = v[218] * v[448] + v[219] * v[449] + v[220] * v[451] + v[260] * v[615] + v[104] * v[616] + v[263] * v[617];
	DdT[0][0] = v[148] * v[2];
	DdT[0][1] = v[152] * v[2];
	DdT[0][2] = v[153] * v[2];
	DdT[0][3] = v[211] * v[352] + v[210] * v[355] + v[209] * v[358] + v[153] * v[412] + v[148] * v[418] + v[152] * v[421];
	DdT[0][4] = v[211] * v[353] + v[210] * v[356] + v[209] * v[359] + v[153] * v[413] + v[148] * v[419] + v[152] * v[422];
	DdT[0][5] = v[211] * v[354] + v[210] * v[357] + v[209] * v[360] + v[153] * v[414] + v[148] * v[420] + v[152] * v[423];
	DdT[1][0] = v[154] * v[2];
	DdT[1][1] = v[158] * v[2];
	DdT[1][2] = v[159] * v[2];
	DdT[1][3] = v[211] * v[361] + v[210] * v[364] + v[209] * v[367] + v[159] * v[412] + v[154] * v[418] + v[158] * v[421];
	DdT[1][4] = v[211] * v[362] + v[210] * v[365] + v[209] * v[368] + v[159] * v[413] + v[154] * v[419] + v[158] * v[422];
	DdT[1][5] = v[211] * v[363] + v[210] * v[366] + v[209] * v[369] + v[159] * v[414] + v[154] * v[420] + v[158] * v[423];
	DdT[2][0] = v[160] * v[2];
	DdT[2][1] = v[164] * v[2];
	DdT[2][2] = v[165] * v[2];
	DdT[2][3] = v[211] * v[370] + v[210] * v[373] + v[209] * v[376] + v[165] * v[412] + v[160] * v[418] + v[164] * v[421];
	DdT[2][4] = v[211] * v[371] + v[210] * v[374] + v[209] * v[377] + v[165] * v[413] + v[160] * v[419] + v[164] * v[422];
	DdT[2][5] = v[211] * v[372] + v[210] * v[375] + v[209] * v[378] + v[165] * v[414] + v[160] * v[420] + v[164] * v[423];
	DdT[3][0] = v[433] * v[436] - v[435] * v[437];
	DdT[3][1] = v[434] * v[436] + v[435] * v[439];
	DdT[3][2] = -(v[434] * v[437]) - v[433] * v[439];
	DdT[3][3] = v[220] * (v[148] * v[252] + v[104] * v[352] - v[266] * v[361] - v[260] * v[370] - v[495] - v[496]) + v[219] *
		(v[152] * v[252] + v[104] * v[355] - v[266] * v[364] - v[260] * v[373] - v[507] - v[508]) + v[218] * (v[153] * v[252]
			+ v[104] * v[358] - v[266] * v[367] - v[260] * v[376] - v[519] - v[520]) + v[436] * v[551] + v[439] * v[554] + v[437] * v[557]
		+ v[591] + v[592] + v[597] + v[598] + v[252] * v[615] - v[260] * v[618] + v[104] * v[619] - v[266] * v[620];
	DdT[3][4] = v[220] * (v[148] * v[254] + v[104] * v[353] - v[266] * v[362] - v[260] * v[371] - v[154] * v[455]
		- v[160] * v[460]) + v[219] * (v[152] * v[254] + v[104] * v[356] - v[266] * v[365] - v[260] * v[374] - v[158] * v[455]
			- v[164] * v[460]) + v[218] * (v[153] * v[254] + v[104] * v[359] - v[266] * v[368] - v[260] * v[377] - v[159] * v[455]
				- v[165] * v[460]) + v[436] * v[552] + v[439] * v[555] + v[437] * v[558] + v[254] * v[615] - v[460] * v[616] - v[455] * v[617]
		- v[260] * v[621] + v[104] * v[622] - v[266] * v[623];
	DdT[3][5] = v[220] * (v[148] * v[255] + v[104] * v[354] - v[266] * v[363] - v[260] * v[372] + v[160] * v[455]
		- v[154] * v[456]) + v[219] * (v[152] * v[255] + v[104] * v[357] - v[266] * v[366] - v[260] * v[375] + v[164] * v[455]
			- v[158] * v[456]) + v[218] * (v[153] * v[255] + v[104] * v[360] - v[266] * v[369] - v[260] * v[378] + v[165] * v[455]
				- v[159] * v[456]) + v[436] * v[553] + v[439] * v[556] + v[437] * v[559] + v[255] * v[615] + v[455] * v[616] - v[456] * v[617]
		- v[260] * v[624] + v[104] * v[625] - v[266] * v[626];
	DdT[4][0] = v[433] * v[442] - v[435] * v[443];
	DdT[4][1] = v[434] * v[442] + v[435] * v[445];
	DdT[4][2] = -(v[434] * v[443]) - v[433] * v[445];
	DdT[4][3] = v[220] * (v[154] * v[252] + v[266] * v[352] + v[104] * v[361] - v[263] * v[370] + v[148] * v[454]
		- v[160] * v[458]) + v[219] * (v[158] * v[252] + v[266] * v[355] + v[104] * v[364] - v[263] * v[373] + v[152] * v[454]
			- v[164] * v[458]) + v[218] * (v[159] * v[252] + v[266] * v[358] + v[104] * v[367] - v[263] * v[376] + v[153] * v[454]
				- v[165] * v[458]) + v[442] * v[551] + v[445] * v[554] + v[443] * v[557] + v[454] * v[615] - v[458] * v[616] + v[252] * v[617]
		- v[263] * v[618] + v[266] * v[619] + v[104] * v[620];
	DdT[4][4] = v[220] * (v[154] * v[254] + v[266] * v[353] + v[104] * v[362] - v[263] * v[371] - v[492] + v[496]) + v[219] *
		(v[158] * v[254] + v[266] * v[356] + v[104] * v[365] - v[263] * v[374] - v[504] + v[508]) + v[218] * (v[159] * v[254]
			+ v[266] * v[359] + v[104] * v[368] - v[263] * v[377] - v[516] + v[520]) + v[442] * v[552] + v[445] * v[555] + v[443] * v[558]
		- v[591] - v[592] + v[599] + v[600] + v[254] * v[617] - v[263] * v[621] + v[266] * v[622] + v[104] * v[623];
	DdT[4][5] = v[220] * (v[154] * v[255] + v[266] * v[354] + v[104] * v[363] - v[263] * v[372] - v[160] * v[454]
		+ v[148] * v[456]) + v[219] * (v[158] * v[255] + v[266] * v[357] + v[104] * v[366] - v[263] * v[375] - v[164] * v[454]
			+ v[152] * v[456]) + v[218] * (v[159] * v[255] + v[266] * v[360] + v[104] * v[369] - v[263] * v[378] - v[165] * v[454]
				+ v[153] * v[456]) + v[442] * v[553] + v[445] * v[556] + v[443] * v[559] + v[456] * v[615] - v[454] * v[616] + v[255] * v[617]
		- v[263] * v[624] + v[266] * v[625] + v[104] * v[626];
	DdT[5][0] = v[433] * v[448] - v[435] * v[449];
	DdT[5][1] = v[434] * v[448] + v[435] * v[451];
	DdT[5][2] = -(v[434] * v[449]) - v[433] * v[451];
	DdT[5][3] = v[220] * (v[160] * v[252] + v[260] * v[352] + v[263] * v[361] + v[104] * v[370] + v[154] * v[458]
		- v[148] * v[459]) + v[219] * (v[164] * v[252] + v[260] * v[355] + v[263] * v[364] + v[104] * v[373] + v[158] * v[458]
			- v[152] * v[459]) + v[218] * (v[165] * v[252] + v[260] * v[358] + v[263] * v[367] + v[104] * v[376] + v[159] * v[458]
				- v[153] * v[459]) + v[448] * v[551] + v[451] * v[554] + v[449] * v[557] - v[459] * v[615] + v[252] * v[616] + v[458] * v[617]
		+ v[104] * v[618] + v[260] * v[619] + v[263] * v[620];
	DdT[5][4] = v[220] * (v[160] * v[254] + v[260] * v[353] + v[263] * v[362] + v[104] * v[371] + v[154] * v[459]
		+ v[148] * v[460]) + v[219] * (v[164] * v[254] + v[260] * v[356] + v[263] * v[365] + v[104] * v[374] + v[158] * v[459]
			+ v[152] * v[460]) + v[218] * (v[165] * v[254] + v[260] * v[359] + v[263] * v[368] + v[104] * v[377] + v[159] * v[459]
				+ v[153] * v[460]) + v[448] * v[552] + v[451] * v[555] + v[449] * v[558] + v[460] * v[615] + v[254] * v[616] + v[459] * v[617]
		+ v[104] * v[621] + v[260] * v[622] + v[263] * v[623];
	DdT[5][5] = v[220] * (v[160] * v[255] + v[260] * v[354] + v[263] * v[363] + v[104] * v[372] + v[492] + v[495]) + v[219] *
		(v[164] * v[255] + v[260] * v[357] + v[263] * v[366] + v[104] * v[375] + v[504] + v[507]) + v[218] * (v[165] * v[255]
			+ v[260] * v[360] + v[263] * v[369] + v[104] * v[378] + v[516] + v[519]) + v[448] * v[553] + v[451] * v[556] + v[449] * v[559]
		- v[597] - v[598] - v[599] - v[600] + v[255] * v[616] + v[104] * v[624] + v[260] * v[625] + v[263] * v[626];
	(*tempkin) = 0.5e0*(v[193] * v[214] + v[197] * v[215] + v[198] * v[216]) + v[207] * v[231] + v[205] * v[232]
		+ v[206] * v[233] + 0.5e0*(v[189] * v[231] + v[188] * v[232] + v[187] * v[233]);
};






