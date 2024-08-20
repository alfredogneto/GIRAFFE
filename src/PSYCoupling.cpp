#include "PSYCoupling.h"
#include"Database.h"
//Variáveis globais
extern
Database db;

PSYCoupling::PSYCoupling()
{
	couple_by_file = true;

	psy_displacements = NULL;
	number_psy_displacements = 0;
	psy_constraints = NULL;
	number_psy_constraints = 0;
}


PSYCoupling::~PSYCoupling()
{
	if (psy_displacements != NULL)
	{
		for (int i = 0; i < number_psy_displacements; i++)
			delete psy_displacements[i];
		delete[] psy_displacements;
	}

	if (psy_constraints != NULL)
	{
		for (int i = 0; i < number_psy_constraints; i++)
			delete psy_constraints[i];
		delete[] psy_constraints;
	}
}

//Leitura no arquivo de entrada do Giraffe
bool PSYCoupling::Read(FILE *f)
{
	char s[1000];
	fscanf(f, "%s", s);
	if (!strcmp(s, "File"))
		couple_by_file = true;
	else
	{
		if (!strcmp(s, "Binary"))
			couple_by_file = false;
		else
			return false;
	}
	fscanf(f, "%s", s);
	if (!strcmp(s, "BoolTable"))
		PSY_bool.Read(f);
	else
		return false;

	return true;
}

//Saída no arquivo de entrada do Giraffe
void PSYCoupling::Write(FILE *f)
{
	fprintf(f, "PSYCoupling\t");
	if (couple_by_file == true)
		fprintf(f, "File\t");
	else
		fprintf(f, "Binary\t");
	PSY_bool.Write(f);
}

//PreCalc
void PSYCoupling::PreCalc()
{
	number_psy_displacements = db.number_particles;
	number_psy_constraints = db.number_particles;
	//Alocação do vetor
	psy_displacements = new Displacement*[number_psy_displacements];
	psy_constraints = new Constraint*[number_psy_constraints];
	for (int i = 0; i < number_psy_displacements; i++)
	{
		psy_displacements[i] = new NodalDisplacement();		//Alocação de cada objeto
		psy_constraints[i] = new NodalConstraint();			//Alocação de cada objeto
		NodalDisplacement* ptr_disp = static_cast<NodalDisplacement*>(psy_displacements[i]);
		NodalConstraint* ptr_constr = static_cast<NodalConstraint*>(psy_constraints[i]);

		ptr_disp->single_node = db.particles[i]->node;
		ptr_disp->table = new Table(2, 6);				//Alocação da tabela dentro de cada objeto (o destrutor do objeto se encarrega de deletar no fim da execução)

		ptr_constr->single_node = db.particles[i]->node;
		//All particles are constrained - during all solution steps
		ptr_constr->UX_table.SetDefault(1);
		ptr_constr->UY_table.SetDefault(1);
		ptr_constr->UZ_table.SetDefault(1);
		ptr_constr->ROTX_table.SetDefault(1);
		ptr_constr->ROTY_table.SetDefault(1);
		ptr_constr->ROTZ_table.SetDefault(1);
	}
}

void PSYCoupling::Couple()
{
	if (couple_by_file == true)
		CoupleByFile();
	else
		CoupleByBinary();

	for (int i = 0; i < number_psy_displacements; i++)
	{
		NodalDisplacement* ptr_disp = static_cast<NodalDisplacement*>(psy_displacements[i]);
		ptr_disp->MountSingleNodeDisplacement();
	}
}

void PSYCoupling::SetConstraints()
{
	for (int i = 0; i < number_psy_displacements; i++)
	{
		NodalConstraint* ptr_constr = static_cast<NodalConstraint*>(psy_constraints[i]);
		ptr_constr->MountSingleNodeConstraint();
	}
}

void PSYCoupling::CoupleByFile()
{
	WritePSYFile();
	cout << "Running PSY - coupling by file..." << "\n";
	_spawnl(P_WAIT, "PSY.exe", "1", NULL);
	ReadPSYFile();
}

void PSYCoupling::CoupleByBinary()
{
	//Criar cópias em variáveis-espelho para passar como referência à dll do PSY - se basear na função WritePSYFile()
	cout << "Running PSY - coupling by binary..." << "\n";
	//Chamar função do PSY - com parâmetros de entrada como sendo a estrutura toda psy_info
	//Ler variáveis-espelho escritas pelo PSY e salvar no Giraffe - se basear na função ReadPSYFile()
}

//Lê o arquivo do PSY
void PSYCoupling::ReadPSYFile()
{
	double time_value = db.last_converged_time + db.current_time_step;
	double dt = db.current_time_step;

	//Leitura e preenchimento da tabela de deslocamentos prescritos
	FILE *f = NULL;
	f = fopen("gir_psy.rst", "r");
	char s[1000];
	fscanf(f, "%s", s);
	while (strcmp(s, "vz"))
		fscanf(f, "%s", s);

	int node = 0;
	Matrix coordinates(3);
	Matrix velocity(3);
	//Leitura dos dados das partículas
	for (int i = 0; i < db.number_particles; i++)
	{
		//Número
		fscanf(f, "%s", s);
		node = db.particles[atoi(s) - 1]->node;
		fscanf(f, "%s", s);
		coordinates(0, 0) = atof(s);
		fscanf(f, "%s", s);
		coordinates(1, 0) = atof(s);
		fscanf(f, "%s", s);
		coordinates(2, 0) = atof(s);
		fscanf(f, "%s", s);
		velocity(0, 0) = atof(s);
		fscanf(f, "%s", s);
		velocity(1, 0) = atof(s);
		fscanf(f, "%s", s);
		velocity(2, 0) = atof(s);

		NodalDisplacement* ptr_disp = static_cast<NodalDisplacement*>(psy_displacements[i]);
		//Deslocamentos prescritos - preenchimento da tabela
		ptr_disp->table->SetTime(0, time_value - dt);
		ptr_disp->table->SetTime(1, time_value);
		ptr_disp->table->SetValue(0, 0, 0);
		ptr_disp->table->SetValue(0, 1, 0);
		ptr_disp->table->SetValue(0, 2, 0);
		ptr_disp->table->SetValue(1, 0, coordinates(0, 0) - db.nodes[node - 1]->copy_coordinates[0]); //deslocamento prescrito em x
		ptr_disp->table->SetValue(1, 1, coordinates(1, 0) - db.nodes[node - 1]->copy_coordinates[1]); //deslocamento prescrito em y
		ptr_disp->table->SetValue(1, 2, coordinates(2, 0) - db.nodes[node - 1]->copy_coordinates[2]); //deslocamento prescrito em z
		//no caso de querer impor rotacoes, colocar aqui.
	}
	fclose(f);
}

//Escreve o arquivo do PSY
void PSYCoupling::WritePSYFile()
{
	double time_value = db.last_converged_time + db.current_time_step;
	double dt = db.current_time_step;

	//Escrita do arquivo PSY
	FILE *f_PSY;
	f_PSY = fopen("gir_psy.inp", "w");
	//Cabeçalho do arquivo PSY
	fprintf(f_PSY, "$system_name\n");
	fprintf(f_PSY, "\t'File generated by Giraffe'\n");
	fprintf(f_PSY, "$no_particles\n");
	fprintf(f_PSY, "\t%d\n", db.number_particles);
	fprintf(f_PSY, "$rotational_dofs !on/off\n");
	fprintf(f_PSY, "\toff\n");
	fprintf(f_PSY, "$particle_attributes_and_radius  !number, kind, radius, material_set_number, pressurf_number\n");
	//Impressão das partículas
	for (int i = 0; i < db.number_particles; i++)
	{
		//Partícula esférica
		//Ponteiro para o Sphere
		Sphere* ptr = static_cast<Sphere*>(db.particles[i]);
		fprintf(f_PSY, "\t%d\tjet_sphere\t%.6e\t%d\t%d\n", ptr->number, ptr->radius, ptr->material, 0);
	}
	//Impressão dos materiais
	fprintf(f_PSY, "$no_material_sets\n");
	fprintf(f_PSY, "\t%d\n", db.number_materials);
	fprintf(f_PSY, "$material_sets  !set_number, mass_dens, charge_dens, E, ni, damping_rate, e, static_friction, dynamic_friction\n");
	for (int i = 0; i < db.number_materials; i++)
	{
		//Ponteiro para o Material
		Hooke* ptr = static_cast<Hooke*>(db.materials[i]);
		double damping_rate = 9.000E-01;
		double e = 1;
		double static_friction = 0.1;
		double dynamic_friction = 0.1;
		fprintf(f_PSY, "\t%d\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\n", ptr->number, ptr->rho, 0.0, ptr->E, ptr->nu, damping_rate, e, static_friction, dynamic_friction);
	}
	//Impressão das coordenadas das partículas e velocidades iniciais
	fprintf(f_PSY, "$particle_coordinates_and_initial_velocities  ! number, x, y, z, vx, vy, vz, wx, wy, wz\n");
	for (int i = 0; i < db.number_particles; i++)
	{
		fprintf(f_PSY, "\t%d\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\n",
			db.particles[i]->number,
			db.nodes[db.particles[i]->node - 1]->copy_coordinates[0],
			db.nodes[db.particles[i]->node - 1]->copy_coordinates[1],
			db.nodes[db.particles[i]->node - 1]->copy_coordinates[2],
			db.nodes[db.particles[i]->node - 1]->copy_vel[0],
			db.nodes[db.particles[i]->node - 1]->copy_vel[1],
			db.nodes[db.particles[i]->node - 1]->copy_vel[2],
			0.0, 0.0, 0.0);
	}
	fprintf(f_PSY, "$no_constrained_particles\n");
	fprintf(f_PSY, "\t%d\n", 0);
	fprintf(f_PSY, "$no_particles_with_initial_given_forces_and_moments\n");
	fprintf(f_PSY, "\t%d\n", db.number_particles);
	//Impressão das forças nas partículas (reações vinculares da simulação do Giraffe)
	fprintf(f_PSY, "$particle_initial_given_forces_and_moments\n");
	int node = 0;
	Matrix forces(3);
	for (int i = 0; i < db.number_particles; i++)
	{
		node = db.particles[i]->node;
		//Esforços
		for (int j = 0; j < 3; j++)
		{
			if (db.nodes[node - 1]->GLs[j] < 0 && db.nodes[node - 1]->active_GL[j] == 1)	//Se o grau de liberdade for fixo e ativo
				forces(j, 0) = db.global_P_B(-db.nodes[node - 1]->GLs[j] - 1, 0);
			else
			{
				printf("Error writing PSY file. DOF supposed to be constrained was set free.\n");
				forces(j, 0) = 0.0;
			}

		}
		fprintf(f_PSY, "\t%d\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\n", db.particles[i]->number, -forces(0, 0), -forces(1, 0), -forces(2, 0), 0.0, 0.0, 0.0);
	}
	//Impressão dos rigid walls
	fprintf(f_PSY, "$no_rigid_walls\n");
	fprintf(f_PSY, "%d\n", 5);
	//fprintf(f_PSY, "%d\n", 1);
	fprintf(f_PSY, "$rigid_wall_properties  !w_num, charge, damp_rate, e, mis, mid, gap_thresh, vx, vy, vz\n");
	//fprintf(f_PSY, "\t1    0.000E+00   1.000E-01   1.000E+00   2.000E-01   2.000E-01   1.000E-11     0.000E+00   0.000E+00   0.000E+00\n");
	//fprintf(f_PSY, "\t1    0.000E+00   1.000E-01   1.000E+00   2.000E-01   2.000E-01   1.000E-11     0.000E+00   0.000E+00   0.000E+00\n");
	fprintf(f_PSY, "\t1    0.000E+00   1.000E-01   1.000E+00   2.000E-01   2.000E-01   1.000E-11     0.000E+00   0.000E+00   0.000E+00\n");
	fprintf(f_PSY, "\t2    0.000E+00   1.000E-01   1.000E+00   2.000E-01   2.000E-01   1.000E-11     0.000E+00   0.000E+00   0.000E+00\n");
	fprintf(f_PSY, "\t3    0.000E+00   1.000E-01   1.000E+00   2.000E-01   2.000E-01   1.000E-11     0.000E+00   0.000E+00   0.000E+00\n");
	fprintf(f_PSY, "\t4    0.000E+00   1.000E-01   1.000E+00   2.000E-01   2.000E-01   1.000E-11     0.000E+00   0.000E+00   0.000E+00\n");
	fprintf(f_PSY, "\t5    0.000E+00   1.000E-01   1.000E+00   2.000E-01   2.000E-01   1.000E-11     0.000E+00   0.000E+00   0.000E+00\n");

	fprintf(f_PSY, "$wall_points_coordinates !wall_number, P1xyz, P2xyz, P3xyz (note: e1=P2-P1;e2=P3-P1;n=e2xe1)\n");
	//fprintf(f_PSY, "\t1\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\n", 7.0E-1, -7.0E-1, 0.0, 7.0E-1, 7.0E-1, 0.0, -7.0E-1, -7.0E-1, 0.0);
	//fprintf(f_PSY, "\t1\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\n", 7.0E-1, -7.0E-1, 0.06, -7.0E-1, -7.0E-1, 0.06, 7.0E-1, 7.0E-1, 0.06);
	fprintf(f_PSY, "\t1\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\n", 7.0E-1, -7.0E-1, -1.0, 7.0E-1, 7.0E-1, -1.0, -7.0E-1, -7.0E-1, -1.0);
	fprintf(f_PSY, "\t2\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\n", 7.0E-1, -7.0E-1, -1.0, 7.0E-1, -7.0E-1, 2.0, 7.0E-1, 7.0E-1, -1.0);
	fprintf(f_PSY, "\t3\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\n", 7.0E-1, 7.0E-1, -1.0, 7.0E-1, 7.0E-1, 2.0, -7.0E-1, 7.0E-1, -1.0);
	fprintf(f_PSY, "\t4\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\n", -7.0E-1, 7.0E-1, -1.0, -7.0E-1, 7.0E-1, 2.0, -7.0E-1, -7.0E-1, -1.0);
	fprintf(f_PSY, "\t5\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\n", -7.0E-1, -7.0E-1, -1.0, -7.0E-1, -7.0E-1, 2.0, 7.0E-1, -7.0E-1, -1.0);

	fprintf(f_PSY, "$periodic_boundary_conditions  !on/off\n");
	fprintf(f_PSY, "\toff\n");
	fprintf(f_PSY, "$environment_forces_data  !gravity_accel_vector, global_damping_coeff\n");
	fprintf(f_PSY, "\t%.6e\t%.6e\t%.6e\t%.6e\n", 0.0, 0.0, 0.0,0.0);
	fprintf(f_PSY, "$near_field_forces  !on/off\n");
	fprintf(f_PSY, "\toff\n");
	fprintf(f_PSY, "$no_pressure_surfaces   !allowed only for a few sphere kinds\n");
	fprintf(f_PSY, "\t0\n");
	fprintf(f_PSY, "$no_springs   !allowed only for a few sphere kinds\n");
	fprintf(f_PSY, "\t0\n");
	fprintf(f_PSY, "$solution_control_variables_1  !problem_type, solver_type, contact_algorithm_type\n");
	fprintf(f_PSY, "\ttransient_dynamics    euler_explicit_solver_for_coupling_with_giraffe    hertzian_continuous_slide\n");
	fprintf(f_PSY, "$solution_control_variables_2   ! tolR, tolV, desired_no_iterations, max_no_iterations\n");
	fprintf(f_PSY, "\t%.6e\t%.6e\t%d\t%d\n", 1.0E-03, 1.0E-03, 6, 200);
	fprintf(f_PSY, "$no_steps\n");
	fprintf(f_PSY, "\t1\n");
	fprintf(f_PSY, "$step_control_variables   !step, initial_dt, max_dt, final_time, final_factor, adaptivity, collission_duraration_parameter\n");
	fprintf(f_PSY, "\t1\t%.6e\t%.6e\t%.6e\t1\toff\t%.6e\n", dt, dt, dt, 1e-2);
	fprintf(f_PSY, "$output_control_variables  !results_file_format, dt_for_results_print, print_system_properties, print_rotational_dofs\n");
	fprintf(f_PSY, "\tpsy_format\t%.6e\tno\tno\n", dt);
	fprintf(f_PSY, "$end\n");
	fclose(f_PSY);
}