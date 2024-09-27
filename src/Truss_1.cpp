#include "Truss_1.h"


#include "Environment.h"
#include "Node.h"
#include "Dynamic.h"
#include "ElasticPlasticIsoHardening.h"
#include "Hooke.h"
#include "Section.h"
#include "PipeSection.h"
#include "Database.h"
//Variaveis globais
extern
Database db;
#define PI 3.1415926535897932384626433832795

Truss_1::Truss_1()
{
	strain_energy = 0.0;
	kinetic_energy = 0.0;
	potential_gravitational_energy = 0.0;

	VTK_type = 3;
	nDOFs = 6;
	material = 0;
	section = 0;
	n_nodes = 2;
	number = 0;
	nodes = new int[n_nodes];
	VTK_nodes = new int[n_nodes];
	VTK_nodes[0] = 0;
	VTK_nodes[1] = 1;

	C1t = 0.0;
	C1n = 0.0;

	T0 = 0.0;

	rho_adt = 0.0;
	rho_adn = 0.0;

	DOFs = new int *[n_nodes];
	for (int i = 0; i < n_nodes; i++)
		DOFs[i] = new int[db.number_GLs_node];

	type_name = new char[20];//Nome do tipo do elemento
	sprintf(type_name, "Truss_1");

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
	}

	//Variaveis do elemento
	T = 0.0;										
	L = 0.0;										
	l = 0.0;										
	tau = 0.0;										
	A = 0.0;										
	a = 0.0;										
	Vol = 0.0;
	Ahydro = 0.0;
	
	E = 0.0;
	H = 0.0; 
	nu = 0.0; 
	rho = 0.0;		
	tau_y_0 = 0.0;

	c_internal_loads = Matrix(6);
	c_inertial_loads = Matrix(6);
	c_damping_loads = Matrix(6);
	c_external_loads = Matrix(6);

	c_stiffness_matrix = Matrix(6,6);
	c_external_loads_stiffness = Matrix(6, 6);
	c_mass_matrix = Matrix(6,6);
	c_damping_matrix = Matrix(6, 6);

	c_mass_modal = Matrix(6, 6);
	c_rayleigh_damping = Matrix(6, 6);
	
	normal = Matrix(3);
	I3 = Matrix(3, 3);
	I3(0, 0) = 1.0;
	I3(1, 1) = 1.0;
	I3(2, 2) = 1.0;

	epsb = 0.0;
	epsp = 0.0;
	l_p = 0.0;

	epsb_i = 0.0;
	epsp_i = 0.0;
	l_p_i = 0.0;

	flag_plastic = false;
	flag_hydro = false;
	pipe_sec = 0;

	//Cada variavel carregara a informação para seu ponto de Gauss
	int ngauss = 2;
	N1 = new double[ngauss];
	N2 = new double[ngauss];
	N = new Matrix*[2];
	//Loop nos pontos de Gauss
	for (int i = 0; i < ngauss; i++)
		N[i] = new Matrix(3, 6);
}

Truss_1::~Truss_1()
{
	delete[] nodes;
	delete[] VTK_nodes;
	if (DOFs != NULL)
	{
		for (int i = 0; i < n_nodes; i++)
			delete[] DOFs[i];
		delete[] DOFs;
	}
	delete[]type_name;

	delete[] N1;
	delete[] N2;
	//Loop nos pontos de Gauss
	for (int i = 0; i < 2; i++)
		delete N[i];
	delete[] N;
}

//Checa inconsistências no elemento para evitar erros de execução
bool Truss_1::Check()
{
	for (int i = 0; i < n_nodes; i++)
	{
		if (nodes[i] > db.number_nodes)
			return false;
	}
	if (material > db.number_materials)
		return false;
	if (section > db.number_sections)
		return false;
	return true;
}

bool Truss_1::Read(FILE *f)
{
	char s[1000];
	fscanf(f, "%s", s);
	number = atoi(s);

	fscanf(f, "%s", s);
	if (!strcmp(s, "Mat"))
	{
		fscanf(f, "%s", s);
		material = atoi(s);

		flag_hydro = false;			//não vai computar esforços hidrodinamicos (não ha coef. disponiveis)

		fscanf(f, "%s", s);
		if (!strcmp(s, "Sec"))
		{
			fscanf(f, "%s", s);
			section = atoi(s);
		}
		else
			return false;
	}
	else
	{
		if (!strcmp(s, "PipeSec"))
		{
			fscanf(f, "%s", s);
			pipe_sec = atoi(s);

			flag_hydro = true;		//vai computar esforços hidrodinamicos (ha coef. disponiveis)
		}
		else
			return false;
	}
		
	fscanf(f, "%s", s);
	if (!strcmp(s, "Nodes"))
	{
		for (int n = 0; n < n_nodes; n++)
		{
			fscanf(f, "%s", s);
			nodes[n] = atoi(s);
		}
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

void Truss_1::Write(FILE *f)
{
	if (flag_hydro == false)
		fprintf(f, "Truss_1\t%d\tMat\t%d\tSec\t%d\tNodes\t%d\t%d\tPreTension\t%.6e\n",
			number,
			material,
			section,
			nodes[0],
			nodes[1],
			T0);
	else
		fprintf(f, "Truss_1\t%d\tPipeSec\t%d\tNodes\t%d\t%d\tPreTension\t%.6e\n",
		number,
		pipe_sec,
		nodes[0],
		nodes[1],
			T0);
}
//Escreve arquivo de resultados
void Truss_1::WriteResults(FILE *f)
{
	fprintf(f, "Truss_1\t%d\tNormalForce\t%.14e\tDeformedArea\t%.14e\tDeformedLength\t%.14e\tDeformedLengthPlastic\t%.14e\tKirschhoffStress\t%.14e\n",number,T,a,l,l_p,tau);
}

//Escreve no monitor do elemento
void Truss_1::WriteMonitor(FILE *f, bool first_record, double time)
{
	//Cabeçalho
	if (first_record == true)
		fprintf(f, "TIME\tTension\tDeformedArea\tDeformedLength\tDeformedLengthPlastic\tKirschhoffStress\t\n");
	//Informações a serem salvas
	fprintf(f, "%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\n",
		time, T,a,l,l_p,tau);
}

void Truss_1::WriteVTK_XMLBase(std::vector<float> *float_vector)
{
	//Imprime os resultados do elemento
	int res_element = 5;
	float_vector->push_back((float)T);
	float_vector->push_back((float)a);
	float_vector->push_back((float)l);
	float_vector->push_back((float)l_p);
	float_vector->push_back((float)tau);
	//Imprime valores nulos para resultados que não fazem sentido para esse tipo de elemento
	for (int i = res_element; i < db.n_element_results; i++)
		float_vector->push_back(0.0);
}
void Truss_1::WriteVTK_XMLRender(FILE *f)
{
	//TODO
	//element index 7
}

//Monta carregamentos associados ao elemento
void Truss_1::MountElementLoads()
{
	////////////////////////////////////////////////////////////////Peso próprio (efetivo)////////////////////////////////////////////////////////
	if (db.environment_exist == true)
	{
		if (db.environment->g_exist == true)
		{
			load_multiplier = 1.0;
			l_factor = db.environment->bool_g.GetLinearFactorAtCurrentTime();

			//Loop nos pontos de Gauss
			Matrix pos_gauss(3);
			double depth;
			double alpha1 = 1.0;
			double jacobian = L / 2;
			for (int gauss = 0; gauss < 2; gauss++)
			{
				if (db.environment->ocean_data_exist == true)
				{
					//Posição do ponto de Gauss
					for (int i = 0; i < 3; i++)
					{
						pos_gauss(i, 0) = (db.nodes[nodes[0] - 1]->copy_coordinates[i])* N1[gauss] + (db.nodes[nodes[1] - 1]->copy_coordinates[i])* N2[gauss] +
							(db.nodes[nodes[0] - 1]->displacements[i])* N1[gauss] + (db.nodes[nodes[1] - 1]->displacements[i])* N2[gauss];
					}
					
					//Verifica profundidade e assinala o valor correto para o rho
					depth = dot((pos_gauss - db.environment->surface_position), db.environment->G)*(1.0 / norm(db.environment->G));
					if (depth > 0)
						mult = l_factor*load_multiplier*alpha1*jacobian*(rho_len - db.environment->rho_fluid*Ahydro);
					else
						mult = l_factor*load_multiplier*alpha1*jacobian*(rho_len);
				}
				else
					mult = l_factor*load_multiplier*alpha1*jacobian*(rho_len);

				if (mult != 0.0)
				{
					c_external_loads(0, 0) += mult * N1[gauss] * db.environment->G(0, 0);
					c_external_loads(1, 0) += mult * N1[gauss] * db.environment->G(1, 0);
					c_external_loads(2, 0) += mult * N1[gauss] * db.environment->G(2, 0);

					c_external_loads(3, 0) += mult * N2[gauss] * db.environment->G(0, 0);
					c_external_loads(4, 0) += mult * N2[gauss] * db.environment->G(1, 0);
					c_external_loads(5, 0) += mult * N2[gauss] * db.environment->G(2, 0);
				}
			}//end of gauss loop
		}
	}
	MountSeaCurrentLoading();
	c_internal_loads = c_internal_loads - c_external_loads;					//Vetor esforço desbalanceado
	c_stiffness_matrix = c_stiffness_matrix - c_external_loads_stiffness;	//Matriz de rigidez tangente
}

//Monta carregamento de correnteza/amortecimento hidro (Morison)
void Truss_1::MountSeaCurrentLoading()
{
	////////////////////////////////////////////////////////////////Sea current load////////////////////////////////////////////////////////
	if (db.environment_exist == true)
	{
		if (db.environment->ocean_data_exist == true)
		{
			load_multiplier = 1.0;
			l_factor = db.environment->bool_current.GetLinearFactorAtCurrentTime();

			//Loop nos pontos de Gauss
			Matrix pos_gauss(3);
			double depth;
			double alpha1 = 1.0;
			double jacobian = L / 2;
			mult = l_factor*load_multiplier*alpha1*jacobian;
			if (mult != 0.0)
			{
				Matrix du_i(3);
				Matrix ddu_i(3);
				Matrix u_d(3);
				Matrix force(3);
				Matrix xA(3);
				Matrix xB(3);
				Matrix uA(3);
				Matrix uB(3);
				Matrix m_stiffness1(3, 3);
				Matrix m_stiffness2(3, 6);
				double** stiff1;
				double** stiff2;
				stiff1 = new double*[3];
				for (int i = 0; i < 3; i++)
					stiff1[i] = new double[3];
				stiff2 = new double*[3];
				for (int i = 0; i < 3; i++)
					stiff2[i] = new double[6];

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
				//Calculando posições e deslocamentos das extremidades
				for (int i = 0; i < 3; i++)
				{
					xA(i, 0) = db.nodes[nodes[0] - 1]->copy_coordinates[i];
					xB(i, 0) = db.nodes[nodes[1] - 1]->copy_coordinates[i];
					uA(i, 0) = db.nodes[nodes[0] - 1]->displacements[i];
					uB(i, 0) = db.nodes[nodes[1] - 1]->displacements[i];
				}
				//Loop nos pontos de Gauss
				for (int gauss = 0; gauss < 2; gauss++)
				{
					//Posição do ponto de Gauss
					for (int i = 0; i < 3; i++)
						pos_gauss(i, 0) = (db.nodes[nodes[0] - 1]->copy_coordinates[i])* N1[gauss] + (db.nodes[nodes[1] - 1]->copy_coordinates[i])* N2[gauss] +
							(db.nodes[nodes[0] - 1]->displacements[i])* N1[gauss] + (db.nodes[nodes[1] - 1]->displacements[i])* N2[gauss];
					
					depth = dot((pos_gauss - db.environment->surface_position), db.environment->G)*(1.0 / norm(db.environment->G));
					Matrix Uinf(3);
					if (depth > 0)
					{
						//velocidade ao longe (correnteza)
						Uinf = db.environment->VelocityAt(depth);
						for (int i = 0; i < 3; i++)
						{
							u_d(i, 0) = (db.nodes[nodes[0] - 1]->displacements[i])* N1[gauss] + (db.nodes[nodes[1] - 1]->displacements[i])* N2[gauss];
							du_i(i, 0) = (db.nodes[nodes[0] - 1]->copy_vel[i])* N1[gauss] + (db.nodes[nodes[1] - 1]->copy_vel[i])* N2[gauss];
							ddu_i(i, 0) = (db.nodes[nodes[0] - 1]->copy_accel[i])* N1[gauss] + (db.nodes[nodes[1] - 1]->copy_accel[i])* N2[gauss];
						}
						EvaluateMorisonContributions(temp_v, &temp_a4, &temp_a5, &temp_a6, &C1t, &C1n, Uinf.getMatrix(), xA.getMatrix(), xB.getMatrix(), uA.getMatrix(), uB.getMatrix(), u_d.getMatrix(),
							du_i.getMatrix(), ddu_i.getMatrix(), force.getMatrix(), stiff1, stiff2);
						m_stiffness1.PtrToMatrix(stiff1, 3);
						m_stiffness2.PtrToMatrix(stiff2, 3, 6);
						
						c_external_loads_stiffness = c_external_loads_stiffness + mult*(transp(*N[gauss]))*m_stiffness1*(*N[gauss]) + mult*(transp(*N[gauss]))*m_stiffness2;
						c_external_loads = c_external_loads + mult*transp(*N[gauss])*force;
					}

				}//end of gauss loop
				//Cleaning
				for (int i = 0; i < 3; i++)
					delete[]stiff1[i];
				delete[]stiff1;
				for (int i = 0; i < 3; i++)
					delete[]stiff2[i];
				delete[]stiff2;
			}
		}
	}
}
//Monta elementos
void Truss_1::Mount()
{
	Zeros();	//Zera matrizes
	//Contributions due to internal loads
	//Position of nodes
	Matrix xa(3);
	Matrix xb(3);
	for (int i = 0; i < 3; i++)
	{
		xa(i, 0) = db.nodes[nodes[0] - 1]->copy_coordinates[i] + db.nodes[nodes[0] - 1]->displacements[i];
		xb(i, 0) = db.nodes[nodes[1] - 1]->copy_coordinates[i] + db.nodes[nodes[1] - 1]->displacements[i];
	}
	//Deformed length
	l = norm(xb - xa);
	//Normal direction
	normal = (1.0/l)*(xb - xa);

	double k = 0.0;//axial stiffness
	//Plasticity
	if (flag_plastic == true)
	{
		//Elastic stretch (try)
		double lambda_try = l / l_p_i;
		double eps_try = log(lambda_try);
		double tau_try = E*eps_try;
		//Evaluation of the yield function
		double f = YieldingFunction(tau_try, epsb_i);
		double E_algo;//Algorithmic Young Modulus
		//Verification - plasticity - return mapping algorithm
		if (f <= 0)
		{
			//No plasticity
			tau = tau_try;
			epsp = epsp_i;
			epsb = epsb_i;
			E_algo = E;
			l_p = l_p_i;
		}
		else
		{
			//Plasticity
			double delta_gamma = f / (E + H);
			epsb = epsb_i + delta_gamma;
			epsp = epsp_i + sign(tau_try)*delta_gamma;
			tau = tau_try - E*sign(tau_try)*delta_gamma;
			E_algo = E*H / (E + H);
			l_p = L*exp(epsp);
		}
		//Jacobian
		double J = pow(l / l_p, 1 - 2 * nu);
		//Cross sectional area
		a = J*Vol / l;
		//Normal force on the truss
		T = tau*a / J;

		k = (Vol / (l*l))*(E_algo - 2 * tau);
	}
	else
	{
		//Elastic stretch
		double lambda = l / L;
		double eps = log(lambda);
		tau = E*eps;
		//Jacobian
		double J = pow(l / L, 1 - 2 * nu);
		//Cross sectional area (no changes)
		a = J*Vol / l;
		//Normal force on the truss
		T = tau*a / J;
		k = (Vol / (l*l))*(E - 2 * tau);
		//No plasticity
		epsp = 0;
		epsb = 0;
		l_p = L;
	}
	Matrix K = k*dyadic(normal, normal) + (T / l)*I3;
	double eps = 1e-6*k;
	if (abs(T) < eps)
		K = k*dyadic(normal, normal) + (eps / l)*I3;
	//Evaluating residual and stiffness contributions
	for (int i = 0; i < 3; i++)
	{
		c_internal_loads(i, 0) = -T*normal(i, 0);
		c_internal_loads(i+3, 0) = +T*normal(i, 0);
		for (int j = 0; j < 3; j++)
		{
			c_stiffness_matrix(i, j) = K(i, j);
			c_stiffness_matrix(i+3, j+3) = K(i, j);
			c_stiffness_matrix(i+3, j) = -K(i, j);
			c_stiffness_matrix(i, j+3) = -K(i, j);
		}
	}
}
//Monta matriz de transformação de coordenadas
void Truss_1::TransformMatrix()
{
	//DOES NOTHING
}
//Zera matrizes locais do elemento
void Truss_1::Zeros()
{
	zeros(&c_external_loads);
	zeros(&c_external_loads_stiffness);
	kinetic_energy = 0.0;
	strain_energy = 0.0;
	potential_gravitational_energy = 0.0;
}
//Preenche a contribuição do elemento nas matrizes globais
void Truss_1::MountGlobal()
{
	int GL_global_1 = 0;
	int GL_global_2 = 0;
	double anterior = 0;
	for (int i = 0; i < 6; i++)
	{
		//Nó 1
		if (i<3)
			GL_global_1 = db.nodes[nodes[0] - 1]->GLs[i];
		//Nó 2
		else
			GL_global_1 = db.nodes[nodes[1] - 1]->GLs[i - 3];

		//Caso o grau de liberdade seja livre:
		if (GL_global_1 > 0)
		{
			anterior = db.global_P_A(GL_global_1 - 1, 0);
			db.global_P_A(GL_global_1 - 1, 0) = anterior + c_internal_loads(i, 0);
			anterior = db.global_I_A(GL_global_1 - 1, 0);
			db.global_I_A(GL_global_1 - 1, 0) = anterior + c_internal_loads(i, 0);
		}
		else
		{
			anterior = db.global_P_B(-GL_global_1 - 1, 0);
			db.global_P_B(-GL_global_1 - 1, 0) = anterior + c_internal_loads(i, 0);
		}
		for (int j = 0; j < 6; j++)
		{
			//Nó 1
			if (j<3)
				GL_global_2 = db.nodes[nodes[0] - 1]->GLs[j];
			//Nó 2
			else
				GL_global_2 = db.nodes[nodes[1] - 1]->GLs[j - 3];

			//Caso os graus de liberdade sejam ambos livres (Matriz Kaa)
			if (GL_global_1 > 0 && GL_global_2 > 0)
				db.global_stiffness_AA.setValue(GL_global_1 - 1, GL_global_2 - 1, c_stiffness_matrix(i, j));
			//Caso os graus de liberdade sejam ambos fixos (Matriz Kbb)
			if (GL_global_1 < 0 && GL_global_2 < 0)
				db.global_stiffness_BB.setValue(-GL_global_1 - 1, -GL_global_2 - 1, c_stiffness_matrix(i, j));
			//Caso os graus de liberdade sejam livre e fixo (Matriz Kab)
			if (GL_global_1 > 0 && GL_global_2 < 0)
				db.global_stiffness_AB.setValue(GL_global_1 - 1, -GL_global_2 - 1, c_stiffness_matrix(i, j));
			//Caso os graus de liberdade sejam fixo e livre (Matriz Kba)
			if (GL_global_1 < 0 && GL_global_2 > 0)
				db.global_stiffness_BA.setValue(-GL_global_1 - 1, GL_global_2 - 1, c_stiffness_matrix(i, j));
		}
	}
}
//Salva variaveis nos pontos de Gauss uteis para descrição lagrangiana atualizada
void Truss_1::SaveLagrange()
{
	//Saving copies of plasticity history variables
	epsb_i = epsb;
	epsp_i = epsp;
	l_p_i = l_p;
}
//Pre-calculo de variaveis e feito uma unica vez no inicio
void Truss_1::PreCalc()
{
	//Position of nodes
	Matrix xa(3);
	Matrix xb(3);
	for (int i = 0; i < 3; i++)
	{
		xa(i, 0) = db.nodes[nodes[0] - 1]->ref_coordinates[i];
		xb(i, 0) = db.nodes[nodes[1] - 1]->ref_coordinates[i];
	}
	//Length
	l = norm(xb - xa);

	if (flag_hydro == false)
	{
		//ElasticPlasticIsoHardening
		if (typeid(*db.materials[material - 1]) == typeid(ElasticPlasticIsoHardening))
		{
			flag_plastic = true;
			ElasticPlasticIsoHardening* mat = static_cast<ElasticPlasticIsoHardening*>(db.materials[material - 1]);
			E = mat->E;
			nu = mat->nu;
			rho = mat->rho;
			H = mat->H;
			tau_y_0 = mat->sigma_y_0;		//assumption - yielding strength
		}
		//Hooke
		if (typeid(*db.materials[material - 1]) == typeid(Hooke))
		{
			flag_plastic = false;
			Hooke* mat = static_cast<Hooke*>(db.materials[material - 1]);
			E = mat->E;
			nu = mat->nu;
			rho = mat->rho;
		}
		//Cross sectional area
		A = db.sections[section - 1]->A;
		Ahydro = A;
		rho_len = rho*A;
		rho_adt = 0.0;
		rho_adn = 0.0;
	}
	else
	{
		E = db.pipe_sections[pipe_sec - 1]->EA;
		A = 1.0;
		Ahydro = PI*db.pipe_sections[pipe_sec - 1]->De*db.pipe_sections[pipe_sec - 1]->De / 4.0;	//area externa
		rho_len = db.pipe_sections[pipe_sec - 1]->Rho;
		double rho_f = 0.0;
		if (db.environment_exist == true)
			rho_f = db.environment->rho_fluid;
		C1t = 0.5*db.pipe_sections[pipe_sec - 1]->De*rho_f*db.pipe_sections[pipe_sec - 1]->CDt;
		C1n = 0.5*db.pipe_sections[pipe_sec - 1]->De*rho_f*db.pipe_sections[pipe_sec - 1]->CDn;
		if (db.environment_exist == true)
		{
			rho_adt = db.environment->rho_fluid*PI*(db.pipe_sections[pipe_sec - 1]->De*db.pipe_sections[pipe_sec - 1]->De)*db.pipe_sections[pipe_sec - 1]->CAt / 4.0;
			rho_adn = db.environment->rho_fluid*PI*(db.pipe_sections[pipe_sec - 1]->De*db.pipe_sections[pipe_sec - 1]->De)*db.pipe_sections[pipe_sec - 1]->CAn / 4.0;
		}
		else
		{
			rho_adt = 0.0;
			rho_adn = 0.0;
		}
	}
	//Undeformed length
	if (T0 == 0.0)
		L = l;
	else
	{
		//PreTension
		double coef = T0 / (E*A);
		L = l / exp(coef);
	}
	//Undeformed volume
	Vol = A * L;
	
	//Zeros on plasticity history variables
	epsb = 0.0;
	epsp = 0.0;
	
	epsb_i = 0.0;
	epsp_i = 0.0;

	l_p = L;
	l_p_i = L;

	//Percorre pontos de Gauss para salvar algumas matrizes/vetores de interesse nos calculos do elemento
	int ngauss = 2;
	double csi = 0.0;
	for (int gauss = 0; gauss < ngauss; gauss++)
	{
		//Ponto localizado nas coordenadas naturais  0.577350269189626 (2 pontos de Gauss)
		if (gauss == 0)
			csi = -0.577350269189626;
		if (gauss == 1)
			csi = +0.577350269189626;
		N1[gauss] = 0.5*(1.0-csi);
		N2[gauss] = 0.5*(1.0+csi);
		(*N[gauss])(0, 0) = N1[gauss];
		(*N[gauss])(1, 1) = N1[gauss];
		(*N[gauss])(2, 2) = N1[gauss];
		(*N[gauss])(0, 3) = N2[gauss];
		(*N[gauss])(1, 4) = N2[gauss];
		(*N[gauss])(2, 5) = N2[gauss];
	}//end of gauss
}

//Monta a matriz de massa
void Truss_1::MountMass()
{
	if (typeid(*db.solution[db.current_solution_number - 1]) == typeid(Dynamic))
	{
		Dynamic* ptr_sol = static_cast<Dynamic*>(db.solution[db.current_solution_number - 1]);
		double coef = 0.0;
		if (flag_hydro == false)
			coef = rho*A;
		else
			coef = rho_len;
		Matrix M1 = (L / 3.0)*I3;
		Matrix M2 = (L / 6.0)*I3;
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				c_mass_matrix(i, j) = coef*M1(i, j);
				c_mass_matrix(i + 3, j + 3) = coef*M1(i, j);
				c_mass_matrix(i + 3, j) = coef*M2(i, j);
				c_mass_matrix(i, j + 3) = coef*M2(i, j);
			}
		}
		Matrix accel(6);
		for (int i = 0; i < 3; i++)
		{
			accel(i, 0) = db.nodes[nodes[0] - 1]->accel[i];
			accel(i + 3, 0) = db.nodes[nodes[1] - 1]->accel[i];
		}


		//Residual contribution - takes current accel
		c_inertial_loads = c_mass_matrix*accel;
		//To tangent operator
		c_mass_matrix = (ptr_sol->a1)*c_mass_matrix;
		MountAddedMass();
	}
}

//Monta a matriz de massa
void Truss_1::MountMassModal()
{
	double coef = 0.0;
	if (flag_hydro == false)
		coef = rho*A;
	else
		coef = rho_len;
	Matrix M1 = (L / 3.0)*I3;
	Matrix M2 = (L / 6.0)*I3;
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			c_mass_modal(i, j) = coef*M1(i, j);
			c_mass_modal(i + 3, j + 3) = coef*M1(i, j);
			c_mass_modal(i + 3, j) = coef*M2(i, j);
			c_mass_modal(i, j + 3) = coef*M2(i, j);
		}
	}
	MountAddedMassModal();
}

//Monta a matriz de amortecimento para realização da analise modal
void Truss_1::MountDampingModal()
{
	//TODO
}

//Monta a matriz de amortecimento
void Truss_1::MountDamping(bool update_rayleigh)
{
	if (typeid(*db.solution[db.current_solution_number - 1]) == typeid(Dynamic))
	{
		Dynamic* ptr_sol = static_cast<Dynamic*>(db.solution[db.current_solution_number - 1]);

		//Amortecimento de Rayleigh
		if (update_rayleigh == true)
		{
			MountMassModal();
			c_rayleigh_damping = ptr_sol->alpha*c_mass_modal + ptr_sol->beta*c_stiffness_matrix;
		}
		c_damping_matrix = ptr_sol->a4*c_rayleigh_damping;	//Matriz de amortecimento - inclusão do efeito de Rayleigh
		Matrix vipp(6);
		for (int i = 0; i < 3; i++)
		{
			vipp(i, 0) = db.nodes[nodes[0] - 1]->vel[i];
			vipp(i + 3, 0) = db.nodes[nodes[1] - 1]->vel[i];
		}
		c_damping_loads = c_rayleigh_damping * vipp;
	}
}

//Montagens - Newmark
void Truss_1::MountDyn()
{
	c_internal_loads = c_internal_loads + c_inertial_loads + c_damping_loads;
	c_stiffness_matrix = c_stiffness_matrix + c_mass_matrix + c_damping_matrix;
}

//Montagens para analise modal - inserção da matriz de massa e amortecimento na matriz de rigidez para posterior montagem global
void Truss_1::MountDynModal()
{
	c_stiffness_matrix = c_mass_modal;
}

//Returns the value of the Yieding function for given tau
double Truss_1::YieldingFunction(double tau, double epsb)
{
	double f = abs(tau) - (tau_y_0 + H*epsb);
	return f;
}
//Returns the sign of a number
double Truss_1::sign(double number)
{
	if (number >= 0.0)
		return +1.0;
	else
		return -1.0;
}

//Monta a matriz de massa adicional
void Truss_1::MountAddedMass()
{
	if (db.environment_exist == true && db.environment->ocean_data_exist == true && db.environment->g_exist)
	{
		//Kinematic variables
		Matrix du_i(3);
		Matrix ddu_i(3);
		Matrix u_d(3);
		Matrix pos_gauss(3);
		//Nodal variables
		Matrix xA(3);
		Matrix xB(3);
		Matrix uA(3);
		Matrix uB(3);
		//force and stiffness variables (to be evaluated by AceGen routine)
		Matrix force(3);
		Matrix m_stiffness1(3, 3);
		Matrix m_stiffness2(3, 6);
		//force and stiffness pointers
		double** stiff1;
		double** stiff2;
		stiff1 = new double*[3];
		for (int i = 0; i < 3; i++)
			stiff1[i] = new double[3];
		stiff2 = new double*[3];
		for (int i = 0; i < 3; i++)
			stiff2[i] = new double[6];
		double temp_a1, temp_a2, temp_a3;
		if (typeid(*db.solution[db.current_solution_number - 1]) == typeid(Dynamic))
		{
			Dynamic* ptr_sol = static_cast<Dynamic*>(db.solution[db.current_solution_number - 1]);
			temp_a1 = ptr_sol->a1;
			temp_a2 = ptr_sol->a2;
			temp_a3 = ptr_sol->a3;
		}
		else
		{
			temp_a1 = 0;
			temp_a2 = 0;
			temp_a3 = 0;
		}
		//Additional auxiliary variables
		double depth;
		double alpha1 = 1.0;
		double jacobian = L / 2;
		mult = alpha1*jacobian;
		//Nodal variables evaluation
		for (int i = 0; i < 3; i++)
		{
			xA(i, 0) = db.nodes[nodes[0] - 1]->copy_coordinates[i];
			xB(i, 0) = db.nodes[nodes[1] - 1]->copy_coordinates[i];
			uA(i, 0) = db.nodes[nodes[0] - 1]->displacements[i];
			uB(i, 0) = db.nodes[nodes[1] - 1]->displacements[i];
		}
		//Loop nos pontos de Gauss
		for (int gauss = 0; gauss < 2; gauss++)
		{
			//Posição do ponto de Gauss
			for (int i = 0; i < 3; i++)
				pos_gauss(i, 0) = (db.nodes[nodes[0] - 1]->copy_coordinates[i])* N1[gauss] + (db.nodes[nodes[1] - 1]->copy_coordinates[i])* N2[gauss] +
				(db.nodes[nodes[0] - 1]->displacements[i])* N1[gauss] + (db.nodes[nodes[1] - 1]->displacements[i])* N2[gauss];

			depth = dot((pos_gauss - db.environment->surface_position), db.environment->G)*(1.0 / norm(db.environment->G));
			if (depth > 0)
			{
				for (int i = 0; i < 3; i++)
				{
					u_d(i, 0) = (db.nodes[nodes[0] - 1]->displacements[i])* N1[gauss] + (db.nodes[nodes[1] - 1]->displacements[i])* N2[gauss];
					du_i(i, 0) = (db.nodes[nodes[0] - 1]->copy_vel[i])* N1[gauss] + (db.nodes[nodes[1] - 1]->copy_vel[i])* N2[gauss];
					ddu_i(i, 0) = (db.nodes[nodes[0] - 1]->copy_accel[i])* N1[gauss] + (db.nodes[nodes[1] - 1]->copy_accel[i])* N2[gauss];
				}
				EvaluateAddedMassContributions(temp_v, &temp_a1, &temp_a2, &temp_a3, &rho_adt, &rho_adn, xA.getMatrix(), xB.getMatrix(), uA.getMatrix(), uB.getMatrix(), u_d.getMatrix(),
					du_i.getMatrix(), ddu_i.getMatrix(), force.getMatrix(), stiff1, stiff2);
				m_stiffness1.PtrToMatrix(stiff1, 3);
				m_stiffness2.PtrToMatrix(stiff2, 3, 6);

				c_mass_matrix = c_mass_matrix + mult*(transp(*N[gauss]))*m_stiffness1*(*N[gauss]) + mult*(transp(*N[gauss]))*m_stiffness2;
				c_inertial_loads = c_inertial_loads + mult*transp(*N[gauss])*force;
			}
		}
		//Cleaning
		for (int i = 0; i < 3; i++)
			delete[]stiff1[i];
		delete[]stiff1;
		for (int i = 0; i < 3; i++)
			delete[]stiff2[i];
		delete[]stiff2;
	}	
}

//Monta a matriz de massa adicional para analise modal
void Truss_1::MountAddedMassModal()
{
	if (db.environment_exist == true && db.environment->ocean_data_exist == true && db.environment->g_exist)
	{
		//Nodal variables
		Matrix xa(3);
		Matrix xb(3);
		Matrix pos_gauss(3);
		//Additional auxiliary variables
		double depth;
		double alpha1 = 1.0;
		double jacobian = L / 2;
		mult = alpha1*jacobian;
		//Position of nodes
		for (int i = 0; i < 3; i++)
		{
			xa(i, 0) = db.nodes[nodes[0] - 1]->copy_coordinates[i] + db.nodes[nodes[0] - 1]->displacements[i];
			xb(i, 0) = db.nodes[nodes[1] - 1]->copy_coordinates[i] + db.nodes[nodes[1] - 1]->displacements[i];
		}
		//Deformed length
		l = norm(xb - xa);
		//Tangent direction
		Matrix t = (1.0 / l)*(xb - xa);
		//Loop nos pontos de Gauss
		for (int gauss = 0; gauss < 2; gauss++)
		{
			//Posição do ponto de Gauss
			for (int i = 0; i < 3; i++)
				pos_gauss(i, 0) = (db.nodes[nodes[0] - 1]->copy_coordinates[i])* N1[gauss] + (db.nodes[nodes[1] - 1]->copy_coordinates[i])* N2[gauss] +
				(db.nodes[nodes[0] - 1]->displacements[i])* N1[gauss] + (db.nodes[nodes[1] - 1]->displacements[i])* N2[gauss];

			depth = dot((pos_gauss - db.environment->surface_position), db.environment->G)*(1.0 / norm(db.environment->G));
			if (depth > 0)
				c_mass_modal = c_mass_modal + mult*(transp(*N[gauss]))*(rho_adt*dyadic(t, t) + rho_adn*(I3 - dyadic(t, t)))*(*N[gauss]);
		}
	}
}

void Truss_1::EvaluateMorisonContributions(double* v, double(*a4), double
	(*a5), double(*a6), double(*C1t), double(*C1n), double* dUinf, double* xA, double* xB
	, double* uA, double* uB, double* ud, double* dui, double* ddui
	, double* force, double** stiffness1, double** stiffness2)
{
	//int b80, b85;
	v[244] = -((*a6)*ddui[2]) - (*a5)*dui[2] + dUinf[2] - (*a4)*ud[2];
	v[243] = -((*a6)*ddui[1]) - (*a5)*dui[1] + dUinf[1] - (*a4)*ud[1];
	v[242] = -((*a6)*ddui[0]) - (*a5)*dui[0] + dUinf[0] - (*a4)*ud[0];
	v[237] = uA[2] - uB[2] + xA[2] - xB[2];
	v[236] = uA[1] - uB[1] + xA[1] - xB[1];
	v[235] = uA[0] - uB[0] + xA[0] - xB[0];
	v[239] = sqrt(Power(v[235], 2) + Power(v[236], 2) + Power(v[237], 2));
	v[238] = 1e0 / v[239];
	v[147] = v[237] * v[238];
	v[146] = v[236] * v[238];
	v[145] = v[235] * v[238];
	v[149] = 1e0 / Power(v[239], 2);
	v[241] = -(v[146] * v[149]);
	v[240] = -(v[147] * v[149]);
	v[154] = 1e0 / v[239] + v[237] * v[240];
	v[153] = v[236] * v[240];
	v[152] = 1e0 / v[239] + v[236] * v[241];
	v[151] = v[235] * v[240];
	v[150] = v[235] * v[241];
	v[148] = -(v[145] * v[149] * v[235]) + 1e0 / v[239];
	v[96] = -((*a4)*(v[145] * v[145]));
	v[102] = -(*a4) - v[96];
	v[94] = -((*a4)*v[146]);
	v[99] = v[146] * v[94];
	v[103] = -(*a4) - v[99];
	v[97] = v[145] * v[94];
	v[95] = -((*a4)*v[147]);
	v[101] = v[147] * v[95];
	v[104] = -(*a4) - v[101];
	v[100] = v[146] * v[95];
	v[98] = v[145] * v[95];
	v[157] = v[151] * v[242] + v[153] * v[243] + v[154] * v[244];
	v[156] = v[150] * v[242] + v[152] * v[243] + v[153] * v[244];
	v[155] = v[148] * v[242] + v[150] * v[243] + v[151] * v[244];
	v[74] = v[145] * v[242] + v[146] * v[243] + v[147] * v[244];
	v[169] = v[147] * v[157] + v[154] * v[74];
	v[167] = v[153] * v[74];
	v[168] = v[147] * v[156] + v[167];
	v[165] = v[151] * v[74];
	v[166] = v[147] * v[155] + v[165];
	v[164] = v[146] * v[157] + v[167];
	v[163] = v[146] * v[156] + v[152] * v[74];
	v[161] = v[150] * v[74];
	v[162] = v[146] * v[155] + v[161];
	v[160] = v[145] * v[157] + v[165];
	v[159] = v[145] * v[156] + v[161];
	v[158] = v[145] * v[155] + v[148] * v[74];
	v[73] = v[145] * v[74];
	v[75] = v[146] * v[74];
	v[76] = v[147] * v[74];
	v[77] = v[242] - v[73];
	v[78] = v[243] - v[75];
	v[79] = v[244] - v[76];
	if (sqrt(Power(v[73], 2) + Power(v[75], 2) + Power(v[76], 2)) == 0e0){
		v[105] = 0e0;
		v[106] = 0e0;
		v[107] = 0e0;
		v[170] = 0e0;
		v[171] = 0e0;
		v[172] = 0e0;
		v[173] = 0e0;
		v[174] = 0e0;
		v[175] = 0e0;
		v[81] = 0e0;
		v[108] = 0e0;
		v[109] = 0e0;
		v[110] = 0e0;
		v[176] = 0e0;
		v[177] = 0e0;
		v[178] = 0e0;
		v[179] = 0e0;
		v[180] = 0e0;
		v[181] = 0e0;
		v[82] = 0e0;
		v[111] = 0e0;
		v[112] = 0e0;
		v[113] = 0e0;
		v[182] = 0e0;
		v[183] = 0e0;
		v[184] = 0e0;
		v[185] = 0e0;
		v[186] = 0e0;
		v[187] = 0e0;
		v[83] = 0e0;
	}
	else {
		v[245] = sqrt(Power(v[73], 2) + Power(v[75], 2) + Power(v[76], 2));
		v[246] = (*C1t)*v[245];
		v[188] = 1e0 / v[245];
		v[191] = v[188] * (v[160] * v[73] + v[164] * v[75] + v[169] * v[76]);
		v[190] = v[188] * (v[159] * v[73] + v[163] * v[75] + v[168] * v[76]);
		v[189] = v[188] * (v[158] * v[73] + v[162] * v[75] + v[166] * v[76]);
		v[118] = v[188] * (v[100] * v[75] + v[101] * v[76] + v[73] * v[98]);
		v[117] = v[188] * (v[100] * v[76] + v[73] * v[97] + v[75] * v[99]);
		v[105] = (*C1t)*(v[245] * v[96] + v[188] * v[73] * (v[73] * v[96] + v[75] * v[97] + v[76] * v[98]));
		v[106] = (*C1t)*(v[117] * v[73] + v[245] * v[97]);
		v[107] = (*C1t)*(v[118] * v[73] + v[245] * v[98]);
		v[170] = (*C1t)*(v[158] * v[245] + v[189] * v[73]);
		v[171] = (*C1t)*(v[159] * v[245] + v[190] * v[73]);
		v[172] = (*C1t)*(v[160] * v[245] + v[191] * v[73]);
		v[173] = -v[170];
		v[174] = -v[171];
		v[175] = -v[172];
		v[81] = v[246] * v[73];
		v[108] = v[106];
		v[109] = (*C1t)*(v[117] * v[75] + v[245] * v[99]);
		v[110] = (*C1t)*(v[100] * v[245] + v[118] * v[75]);
		v[176] = (*C1t)*(v[162] * v[245] + v[189] * v[75]);
		v[177] = (*C1t)*(v[163] * v[245] + v[190] * v[75]);
		v[178] = (*C1t)*(v[164] * v[245] + v[191] * v[75]);
		v[179] = -v[176];
		v[180] = -v[177];
		v[181] = -v[178];
		v[82] = v[246] * v[75];
		v[111] = v[107];
		v[112] = v[110];
		v[113] = (*C1t)*(v[101] * v[245] + v[118] * v[76]);
		v[182] = (*C1t)*(v[166] * v[245] + v[189] * v[76]);
		v[183] = (*C1t)*(v[168] * v[245] + v[190] * v[76]);
		v[184] = (*C1t)*(v[169] * v[245] + v[191] * v[76]);
		v[185] = -v[182];
		v[186] = -v[183];
		v[187] = -v[184];
		v[83] = v[246] * v[76];
	};
	if (sqrt(Power(v[77], 2) + Power(v[78], 2) + Power(v[79], 2)) == 0e0){
		v[119] = 0e0;
		v[120] = 0e0;
		v[121] = 0e0;
		v[192] = 0e0;
		v[193] = 0e0;
		v[194] = 0e0;
		v[195] = 0e0;
		v[196] = 0e0;
		v[197] = 0e0;
		v[86] = 0e0;
		v[122] = 0e0;
		v[123] = 0e0;
		v[124] = 0e0;
		v[198] = 0e0;
		v[199] = 0e0;
		v[200] = 0e0;
		v[201] = 0e0;
		v[202] = 0e0;
		v[203] = 0e0;
		v[87] = 0e0;
		v[125] = 0e0;
		v[126] = 0e0;
		v[127] = 0e0;
		v[204] = 0e0;
		v[205] = 0e0;
		v[206] = 0e0;
		v[207] = 0e0;
		v[208] = 0e0;
		v[209] = 0e0;
		v[88] = 0e0;
	}
	else {
		v[251] = (*C1n)*v[79];
		v[250] = (*C1n)*v[78];
		v[249] = (*C1n)*v[77];
		v[247] = sqrt(Power(v[77], 2) + Power(v[78], 2) + Power(v[79], 2));
		v[248] = -((*C1n)*v[247]);
		v[210] = 1e0 / v[247];
		v[213] = -(v[210] * (v[160] * v[77] + v[164] * v[78] + v[169] * v[79]));
		v[212] = -(v[210] * (v[159] * v[77] + v[163] * v[78] + v[168] * v[79]));
		v[211] = -(v[210] * (v[158] * v[77] + v[162] * v[78] + v[166] * v[79]));
		v[132] = -(v[210] * (v[100] * v[78] - v[104] * v[79] + v[77] * v[98]));
		v[131] = -(v[210] * (-(v[103] * v[78]) + v[100] * v[79] + v[77] * v[97]));
		v[129] = -(v[210] * (-(v[102] * v[77]) + v[78] * v[97] + v[79] * v[98]));
		v[135] = v[100] * v[248];
		v[134] = v[248] * v[98];
		v[133] = v[248] * v[97];
		v[119] = (*C1n)*(v[102] * v[247] + v[129] * v[77]);
		v[120] = v[133] + v[131] * v[249];
		v[121] = v[134] + v[132] * v[249];
		v[192] = (*C1n)*(-(v[158] * v[247]) + v[211] * v[77]);
		v[193] = (*C1n)*(-(v[159] * v[247]) + v[212] * v[77]);
		v[194] = (*C1n)*(-(v[160] * v[247]) + v[213] * v[77]);
		v[195] = -v[192];
		v[196] = -v[193];
		v[197] = -v[194];
		v[86] = -(v[248] * v[77]);
		v[122] = v[133] + v[129] * v[250];
		v[123] = (*C1n)*(v[103] * v[247] + v[131] * v[78]);
		v[124] = v[135] + v[132] * v[250];
		v[198] = (*C1n)*(-(v[162] * v[247]) + v[211] * v[78]);
		v[199] = (*C1n)*(-(v[163] * v[247]) + v[212] * v[78]);
		v[200] = (*C1n)*(-(v[164] * v[247]) + v[213] * v[78]);
		v[201] = -v[198];
		v[202] = -v[199];
		v[203] = -v[200];
		v[87] = -(v[248] * v[78]);
		v[125] = v[134] + v[129] * v[251];
		v[126] = v[135] + v[131] * v[251];
		v[127] = (*C1n)*(v[104] * v[247] + v[132] * v[79]);
		v[204] = (*C1n)*(-(v[166] * v[247]) + v[211] * v[79]);
		v[205] = (*C1n)*(-(v[168] * v[247]) + v[212] * v[79]);
		v[206] = (*C1n)*(-(v[169] * v[247]) + v[213] * v[79]);
		v[207] = -v[204];
		v[208] = -v[205];
		v[209] = -v[206];
		v[88] = -(v[248] * v[79]);
	};
	force[0] = v[81] + v[86];
	force[1] = v[82] + v[87];
	force[2] = v[83] + v[88];
	stiffness1[0][0] = v[105] + v[119];
	stiffness1[0][1] = v[106] + v[120];
	stiffness1[0][2] = v[107] + v[121];
	stiffness1[1][0] = v[108] + v[122];
	stiffness1[1][1] = v[109] + v[123];
	stiffness1[1][2] = v[110] + v[124];
	stiffness1[2][0] = v[111] + v[125];
	stiffness1[2][1] = v[112] + v[126];
	stiffness1[2][2] = v[113] + v[127];
	stiffness2[0][0] = v[170] + v[192];
	stiffness2[0][1] = v[171] + v[193];
	stiffness2[0][2] = v[172] + v[194];
	stiffness2[0][3] = v[173] + v[195];
	stiffness2[0][4] = v[174] + v[196];
	stiffness2[0][5] = v[175] + v[197];
	stiffness2[1][0] = v[176] + v[198];
	stiffness2[1][1] = v[177] + v[199];
	stiffness2[1][2] = v[178] + v[200];
	stiffness2[1][3] = v[179] + v[201];
	stiffness2[1][4] = v[180] + v[202];
	stiffness2[1][5] = v[181] + v[203];
	stiffness2[2][0] = v[182] + v[204];
	stiffness2[2][1] = v[183] + v[205];
	stiffness2[2][2] = v[184] + v[206];
	stiffness2[2][3] = v[185] + v[207];
	stiffness2[2][4] = v[186] + v[208];
	stiffness2[2][5] = v[187] + v[209];
}

void Truss_1::EvaluateAddedMassContributions(double* v, double(*a1), double(*a2), double(*a3)
	, double(*rhoadt), double(*rhoadn), double* xA, double* xB
	, double* uA, double* uB, double* ud, double* dui, double* ddui
	, double* force, double** stiffness1, double** stiffness2)
{
	v[149] = -(*rhoadn) + (*rhoadt);
	v[143] = -((*a3)*ddui[2]) - (*a2)*dui[2] + (*a1)*ud[2];
	v[154] = v[143] * v[149];
	v[142] = -((*a3)*ddui[1]) - (*a2)*dui[1] + (*a1)*ud[1];
	v[153] = v[142] * v[149];
	v[141] = -((*a3)*ddui[0]) - (*a2)*dui[0] + (*a1)*ud[0];
	v[152] = v[141] * v[149];
	v[139] = uA[2] - uB[2] + xA[2] - xB[2];
	v[138] = uA[1] - uB[1] + xA[1] - xB[1];
	v[137] = uA[0] - uB[0] + xA[0] - xB[0];
	v[144] = sqrt(Power(v[137], 2) + Power(v[138], 2) + Power(v[139], 2));
	v[140] = 1e0 / v[144];
	v[87] = v[139] * v[140];
	v[150] = 2e0*v[87];
	v[86] = v[138] * v[140];
	v[151] = v[149] * v[86];
	v[148] = 2e0*v[86];
	v[85] = v[137] * v[140];
	v[147] = 2e0*v[85];
	v[89] = 1e0 / Power(v[144], 2);
	v[146] = -(v[86] * v[89]);
	v[145] = -(v[87] * v[89]);
	v[99] = 1e0 / v[144] + v[139] * v[145];
	v[93] = v[138] * v[145];
	v[92] = 1e0 / v[144] + v[138] * v[146];
	v[91] = v[137] * v[145];
	v[90] = v[137] * v[146];
	v[88] = 1e0 / v[144] - v[137] * v[85] * v[89];
	v[106] = v[147] * v[91];
	v[105] = v[147] * v[90];
	v[104] = v[147] * v[88];
	v[112] = v[148] * v[93];
	v[111] = v[148] * v[92];
	v[110] = v[148] * v[90];
	v[98] = v[148] * v[149] * v[91];
	v[131] = v[141] * v[98];
	v[129] = v[142] * v[98];
	v[125] = v[143] * v[98];
	v[95] = v[149] * (v[86] * v[90] + v[85] * v[92]);
	v[94] = v[149] * (v[86] * v[88] + v[85] * v[90]);
	v[69] = v[151] * v[85];
	v[118] = v[150] * v[99];
	v[117] = v[150] * v[93];
	v[116] = v[150] * v[91];
	v[103] = v[149] * (v[87] * v[91] + v[85] * v[99]);
	v[102] = v[149] * (v[87] * v[88] + v[85] * v[91]);
	v[101] = v[149] * (v[87] * v[93] + v[86] * v[99]);
	v[100] = v[149] * (v[87] * v[92] + v[86] * v[93]);
	v[73] = v[151] * v[87];
	v[72] = v[149] * v[85] * v[87];
	v[67] = (v[85] * v[85]);
	v[76] = -((*rhoadn)*(-1e0 + v[67])) + (*rhoadt)*v[67];
	v[70] = (v[86] * v[86]);
	v[80] = -((*rhoadn)*(-1e0 + v[70])) + (*rhoadt)*v[70];
	v[74] = (v[87] * v[87]);
	v[83] = -((*rhoadn)*(-1e0 + v[74])) + (*rhoadt)*v[74];
	v[78] = (*a1)*v[69];
	v[79] = (*a1)*v[72];
	v[82] = (*a1)*v[73];
	v[122] = v[102] * v[143] + v[104] * v[152] + v[142] * v[94];
	v[123] = v[125] + v[110] * v[153] + v[141] * v[94];
	v[124] = v[129] + v[102] * v[141] + v[116] * v[154];
	v[126] = v[125] + v[105] * v[152] + v[142] * v[95];
	v[127] = v[100] * v[143] + v[111] * v[153] + v[141] * v[95];
	v[128] = v[131] + v[100] * v[142] + v[117] * v[154];
	v[130] = v[129] + v[103] * v[143] + v[106] * v[152];
	v[132] = v[131] + v[101] * v[143] + v[112] * v[153];
	v[133] = v[103] * v[141] + v[101] * v[142] + v[118] * v[154];
	force[0] = v[142] * v[69] + v[143] * v[72] + v[141] * v[76];
	force[1] = v[141] * v[69] + v[143] * v[73] + v[142] * v[80];
	force[2] = v[141] * v[72] + v[142] * v[73] + v[143] * v[83];
	stiffness1[0][0] = (*a1)*v[76];
	stiffness1[0][1] = v[78];
	stiffness1[0][2] = v[79];
	stiffness1[1][0] = v[78];
	stiffness1[1][1] = (*a1)*v[80];
	stiffness1[1][2] = v[82];
	stiffness1[2][0] = v[79];
	stiffness1[2][1] = v[82];
	stiffness1[2][2] = (*a1)*v[83];
	stiffness2[0][0] = v[122];
	stiffness2[0][1] = v[126];
	stiffness2[0][2] = v[130];
	stiffness2[0][3] = -v[122];
	stiffness2[0][4] = -v[126];
	stiffness2[0][5] = -v[130];
	stiffness2[1][0] = v[123];
	stiffness2[1][1] = v[127];
	stiffness2[1][2] = v[132];
	stiffness2[1][3] = -v[123];
	stiffness2[1][4] = -v[127];
	stiffness2[1][5] = -v[132];
	stiffness2[2][0] = v[124];
	stiffness2[2][1] = v[128];
	stiffness2[2][2] = v[133];
	stiffness2[2][3] = -v[124];
	stiffness2[2][4] = -v[128];
	stiffness2[2][5] = -v[133];
}