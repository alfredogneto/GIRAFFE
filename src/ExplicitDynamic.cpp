#include "ExplicitDynamic.h"
#include <iostream>
using namespace std;
#include <chrono>
using namespace std::chrono;

#include "Matrix.h"
#include "ConcomitantSolution.h"
#include "PostFiles.h"
#include "Monitor.h"
#include "GeneralContactSearch.h"
#include "ConfigurationSave.h"
#include "Node.h"
#include "Particle.h"
#include "Load.h"
#include "Displacement.h"
#include "Boundary.h"

#include"Database.h"
//Variaveis globais
extern
Database db;

ExplicitDynamic::ExplicitDynamic()
{
	file_index = 1;								//Numero do arquivo para salvar resultados
	zero_IC_flag = false;
	solution_number = 0;

	start_time = 0.0;
	end_time = 0;
	i_time_step = 0;
	max_time_step = 0;
	min_time_step = 0;
	sample = 0;

	//Damping
	alpha = 0;
	beta = 0;
	update = 0;	//booleana

	I3 = new Matrix(3, 3);
	(*I3)(0, 0) = 1.0;
	(*I3)(1, 1) = 1.0;
	(*I3)(2, 2) = 1.0;

	method = new char[20];
	sprintf(method, "Euler");
}

ExplicitDynamic::~ExplicitDynamic()
{
	delete I3;
	delete[] method;
}

//Reads input file
bool ExplicitDynamic::Read(FILE *f)
{
	char s[1000];
	fscanf(f, "%s", s);
	solution_number = atoi(s);

	fscanf(f, "%s", s);
	if (!strcmp(s, "EndTime"))
	{
		fscanf(f, "%s", s);
		end_time = atof(s);
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "TimeStep"))
	{
		fscanf(f, "%s", s);
		i_time_step = atof(s);
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "MaxTimeStep"))
	{
		fscanf(f, "%s", s);
		max_time_step = atof(s);
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "MinTimeStep"))
	{
		fscanf(f, "%s", s);
		min_time_step = atof(s);
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "Method"))
	{
		bool method_OK = false;
		fscanf(f, "%s", s);
		if (!strcmp(s, "Euler"))
		{
			sprintf(method, "Euler");
			method_OK = true;
		}
		if (!strcmp(s, "RungeKutta4"))
		{
			sprintf(method, "RungeKutta4");
			method_OK = true;
		}
		if (method_OK == false)
			return false;
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "Sample"))
	{
		fscanf(f, "%s", s);
		sample = atoi(s);
	}
	else
		return false;

	

	//Salva a posição (stream)
	fpos_t pos;
	fgetpos(f, &pos);
	fscanf(f, "%s", s);
	if (!strcmp(s, "ZeroIC"))
		zero_IC_flag = true;
	else
		fsetpos(f, &pos);

	return true;
}

//Writes output file
void ExplicitDynamic::Write(FILE *f)
{
	fprintf(f, "ExplicitDynamic\t%d\nEndTime\t%.6e\nTimeStep\t%.6e\nMaxTimeStep\t%.6e\nMinTimeStep\t%.6e\nMethod\t%s\nSample\t%d\n\n",
		solution_number, end_time, i_time_step, max_time_step, min_time_step, method, sample);
	if (zero_IC_flag)
		fprintf(f, "ZeroIC\n");
}

//Escreve resultados
void ExplicitDynamic::WriteResults()
{
	//Atualiza arquivos de pós-processamento
	db.post_files->UpdateSinglePartPostFiles(solution_number, time, file_index);
	db.post_files->WriteConfigurationResults(solution_number, time, file_index);
	file_index++;
}

void ExplicitDynamic::SetGlobalSizeExplicit()
{
	int temp_size_A = 0;
	int temp_size_B = 0;
	int temp_DOF_free = 0;
	int temp_DOF_fixed = 0;
	int temp_node = 0;
	int temp_element = 0;
	
	//FORÇAS
	db.global_P_A.setLines(db.n_GL_free);
	db.global_P_A.setColumns(1);
	db.global_P_A.alloc();

	db.global_I_A.setLines(db.n_GL_free);
	db.global_I_A.setColumns(1);
	db.global_I_A.alloc();

	db.global_P_B.setLines(db.n_GL_fixed);
	db.global_P_B.setColumns(1);
	db.global_P_B.alloc();

	//DESLOCAMENTOS PRESCRITOS
	db.global_X_B.setLines(db.n_GL_fixed);
	db.global_X_B.setColumns(1);
	db.global_X_B.alloc();
}


//Solves solution routine
bool ExplicitDynamic::Solve()
{
	bool aborted = false;						//flag para abortar simulação
	db.last_converged_time = start_time;
	time = db.last_converged_time;
	//Atualiza monitor - somente se for o primeira solution
	if (db.monitor_exist == true && solution_number == 1)
		db.monitor->UpdateMonitor(db.last_converged_time);
	//Atualiza analise concomitante
	if (db.concomitant_solution_exist == true)
		db.concomitant_solution->UpdateConcomitantSolution(db.last_converged_time);
	WriteResults();								//salvando resultados
	high_resolution_clock::time_point t1 = high_resolution_clock::now();//Inicia a marcação de tempo de execução
	//Plotagem do dia da simulação
	system_clock::time_point today = system_clock::now();
	time_t tt;
	tt = system_clock::to_time_t(today);
	db.myprintf("\nSolution step %d started at %s\n", solution_number, ctime(&tt));
	DOFsActive();								//Para cada nó ativa DOFs - tambem opera sobre multiplicadores de Lagrange de SpecialConstraints
	SetGlobalDOFs();							//Numeração de graus de liberdade
	SetGlobalSizeExplicit();					//Calcula o tamanho de contribuições globais - com base nos GLs livres e fixos

	//Time variables
	time_step = i_time_step;
	InitialEvaluations();
	ComputeInitialConditions(zero_IC_flag);	//Computa as condições iniciais nodais impostas
	bool first = true;						//Flag para saber quando computar Rayleigh Damping

	if (db.gcs_exist)
		if (db.gcs->bool_table.GetAt(db.current_solution_number - 1))
			db.gcs->SolutionStepInitialCheck();					//General contact initial check

	aborted = false;	//Flag que indica se o caso deve ser abortado
	
	int steps_to_increase = 3;
	double inc_factor = 1.5;

	
	int steps_computed = 0;
	int steps_control = 0;
	////////////////////////////////////////////////////////////////////////////////
	while (time < end_time  && aborted == false)
	{
		//Setando variaveis no database
		db.last_converged_time = time;
		time += time_step;	//Incremento do time step, de acordo com a progressão da solução
		if (time > end_time)//Se ja estiver alem do ultimo instante de interesse	
		{
			time = end_time;
			time_step = end_time - db.last_converged_time;
		}
		db.current_time_step = time_step;
		db.myprintf("\nTime: %.9f\t\tTime step: %.9f\n\n", time, time_step);
		Zeros();	//zera deslocamentos e rotações incrementais nos nós
		
		if (!strcmp(method, "Euler"))
			Euler();
		if (!strcmp(method, "RungeKutta4"))
			RungeKutta4();
		
		//Teste para avaliar que o time step esta muito pequeno - ABORT SIMULATION
		if (time_step < min_time_step && time != end_time)
		{
			db.myprintf("\nAborting simulation. Step size is too small!\n\n");
			aborted = true;
		}
		FinalUpdateContactsExplicit(time);

		/////////////////////////////////SALVANDO CONFIGURAÇÃO/////////////////////////////////
		if (HaveErrors() != true)
		{
			steps_computed++;
			steps_control++;
			SaveConfiguration();	//Salva configuração convergida
			
			//Analise da facilidade de convergência:
			if (steps_computed >= steps_to_increase)
			{
				time_step = time_step * inc_factor;
				if (time_step > max_time_step)
					time_step = max_time_step;
				steps_computed = 0;	//Zera o convergence counter
			}

			//Escrita em arquivos de resultados
			if (aborted == false)	//Salva o arquivo, de acordo com amostragem
			{
				if ((steps_control%sample == 0 || time == end_time))
					WriteResults();								//salvando resultados
				//Atualiza monitor
				if (db.monitor_exist == true && (steps_control%db.monitor->sample == 0 || time == end_time))
					db.monitor->UpdateMonitor(time);
				//Atualiza analise concomitante
				if (db.concomitant_solution_exist == true && (steps_control%db.concomitant_solution->sample == 0 || time == end_time))
					db.concomitant_solution->UpdateConcomitantSolution(time);
				//Atualiza configuration save
				if (db.config_save_exist == true && (steps_control%db.config_save->sample == 0 || time == end_time))
					db.config_save->ExportConfiguration(time);
			}
		}
		/////////////////////////////////RETOMANDO CONFIGURAÇÃO/////////////////////////////////
		else
		{
			RestoreConfiguration();	//Restaura a ultima configuração que convergiu
			//Modificação do loading factor increment:
			time -= time_step;
			time_step = time_step / 2.0;	 //Bissecção
			steps_computed = 0;
		}
	}//Fim
	high_resolution_clock::time_point t2 = high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
	db.myprintf("\nSolution step %d time:\t   %lf sec.\n", solution_number, duration / 1e6);
	printf("\nSolution step %d time:\t   %lf sec.\n", solution_number, duration / 1e6);
	//Plotagem do dia da simulação
	today = system_clock::now();
	tt = system_clock::to_time_t(today);
	db.myprintf("\nSolution step %d finished at %s\n", solution_number, ctime(&tt));

	if (aborted == true)
		return false;//Erro
	else
		return true;//Não erro
}

void ExplicitDynamic::Euler()
{
	//Evaluates prescribed displacements
	MountDisplacementsExplicit(time - time_step);
	//Evaluates inertial and distributed loads contributions on particles and bodies
	MountExplicit();
	//Evaluates applied loads contributions on global vector of loads
	MountLoadsExplicit(time - time_step);
	MontSpecialConstraintsExplicit(time - time_step);
	MountContactsExplicit(time - time_step);
	MountGlobalExplicit();

	//Evaluating accelerations - contributions stemming from particles
	//Evaluating accelerations - contributions stemming from bodies
	//Loop on bodies (TO DO)
	//...
#pragma omp parallel
	{
		//updates variables associated with each particle (centroid and rotation tensor)
#pragma omp for schedule(dynamic)
		for (int i = 0; i < db.number_particles; i++)
			db.particles[i]->EvaluateAccelerations();
	}
#pragma omp parallel
	{
		//Time-integration
#pragma omp for schedule(dynamic)
		for (int i = 0; i < db.number_nodes; i++)
		{
			//Rotation variables
			Matrix Qp;
			//Auxiliary variables
			Qp = (*db.nodes[i]->Q0) *(*db.nodes[i]->Q);
			//Displacements
			(*db.nodes[i]->du) = (*db.nodes[i]->du) + time_step * Qp * (*db.nodes[i]->ddu);
			(*db.nodes[i]->u) = time_step * (*db.nodes[i]->du);
			//Rotations
			(*db.nodes[i]->omega) = (*db.nodes[i]->omega) + time_step * (*db.nodes[i]->domega);
			(*db.nodes[i]->alpha) = time_step * (*db.nodes[i]->omega);

		}
	}

#pragma omp parallel
	{
		//Final updates - nodes
#pragma omp for schedule(dynamic)
		for (int i = 0; i < db.number_nodes; i++)
		{
			//Rotation variables
			Matrix Qp;
			//Auxiliary variables
			Qp = (*db.nodes[i]->Q0) *(*db.nodes[i]->Q);
			//To global (acceleration)
			(*db.nodes[i]->ddu) = Qp * (*db.nodes[i]->ddu);
			//Updating
			for (int k = 0; k < 3; k++)
			{
				if (db.nodes[i]->constraints[k] == 0 && db.nodes[i]->active_GL[k] == 1)
				{
					db.nodes[i]->displacements[k] = (*db.nodes[i]->u)(k, 0);
					db.nodes[i]->vel[k] = (*db.nodes[i]->du)(k, 0);
					db.nodes[i]->accel[k] = (*db.nodes[i]->ddu)(k, 0);

				}
				else
				{
					(*db.nodes[i]->u)(k, 0) = db.nodes[i]->displacements[k];
					(*db.nodes[i]->du)(k, 0) = db.nodes[i]->vel[k];
					(*db.nodes[i]->ddu)(k, 0) = db.nodes[i]->accel[k];
				}
			}
			//Updating
			for (int k = 0; k < 3; k++)
			{
				if (db.nodes[i]->constraints[k + 3] == 0 && db.nodes[i]->active_GL[k + 3] == 1)
				{
					db.nodes[i]->displacements[k + 3] = (*db.nodes[i]->alpha)(k, 0);
					db.nodes[i]->vel[k + 3] = (*db.nodes[i]->omega)(k, 0);
					db.nodes[i]->accel[k + 3] = (*db.nodes[i]->domega)(k, 0);
				}
				else
				{
					(*db.nodes[i]->alpha)(k, 0) = db.nodes[i]->displacements[k + 3];
					(*db.nodes[i]->omega)(k, 0) = db.nodes[i]->vel[k + 3];
					(*db.nodes[i]->domega)(k, 0) = db.nodes[i]->accel[k + 3];
				}
			}
		}
	}
}

void ExplicitDynamic::RungeKutta4()
{
	//Rotation variables
	Matrix* rk1;
	Matrix* rK1;
	Matrix* rk2;
	Matrix* rK2;
	Matrix* rk3;
	Matrix* rK3;
	Matrix* rk4;
	Matrix* rK4;
	//Displacement variables
	Matrix* dk1;
	Matrix* dK1;
	Matrix* dk2;
	Matrix* dK2;
	Matrix* dk3;
	Matrix* dK3;
	Matrix* dk4;
	Matrix* dK4;
	//Rotation variables
	rk1 = new Matrix[db.number_nodes];
	rK1 = new Matrix[db.number_nodes];
	rk2 = new Matrix[db.number_nodes];
	rK2 = new Matrix[db.number_nodes];
	rk3 = new Matrix[db.number_nodes];
	rK3 = new Matrix[db.number_nodes];
	rk4 = new Matrix[db.number_nodes];
	rK4 = new Matrix[db.number_nodes];
	//Displacement variables
	dk1 = new Matrix[db.number_nodes];
	dK1 = new Matrix[db.number_nodes];
	dk2 = new Matrix[db.number_nodes];
	dK2 = new Matrix[db.number_nodes];
	dk3 = new Matrix[db.number_nodes];
	dK3 = new Matrix[db.number_nodes];
	dk4 = new Matrix[db.number_nodes];
	dK4 = new Matrix[db.number_nodes];
	for (int i = 0; i < db.number_nodes; i++)
	{
		rk1[i] = Matrix(3);
		rK1[i] = Matrix(3);
		rk2[i] = Matrix(3);
		rK2[i] = Matrix(3);
		rk3[i] = Matrix(3);
		rK3[i] = Matrix(3);
		rk4[i] = Matrix(3);
		rK4[i] = Matrix(3);

		dk1[i] = Matrix(3);
		dK1[i] = Matrix(3);
		dk2[i] = Matrix(3);
		dK2[i] = Matrix(3);
		dk3[i] = Matrix(3);
		dK3[i] = Matrix(3);
		dk4[i] = Matrix(3);
		dK4[i] = Matrix(3);
	}
	
	//Loop 1
	//Evaluates prescribed displacements
	MountDisplacementsExplicit(time - time_step);
	//Evaluates inertial and distributed loads contributions on particles and bodies
	MountExplicit();
	//Evaluates applied loads contributions on global vector of loads
	MountLoadsExplicit(time - time_step);
	MontSpecialConstraintsExplicit(time - time_step);
	MountContactsExplicit(time - time_step);
	MountGlobalExplicit();
	//db.global_P_A.print();
#pragma omp parallel
	{
		//Evaluating accelerations - contributions stemming from particles
#pragma omp for schedule(dynamic)
		for (int i = 0; i < db.number_particles; i++)
			db.particles[i]->EvaluateAccelerations();
	}
	//Evaluating accelerations - contributions stemming from bodies
	//Loop on bodies (TO DO)
	//...

#pragma omp parallel
	{
#pragma omp for schedule(dynamic)
		for (int i = 0; i < db.number_nodes; i++)
		{
			//Rotation variables
			Matrix A;
			double g;
			Matrix Qdelta;
			Matrix Qp;
			//Auxiliary variables
			A = skew(*db.nodes[i]->alpha);
			g = 4.0 / (4.0 + norm(*db.nodes[i]->alpha)*norm(*db.nodes[i]->alpha));
			Qdelta = *I3 + g * (A + 0.5*(A)*(A));
			Qp = (*db.nodes[i]->Q0) * (*db.nodes[i]->Q) * Qdelta;
			//Displacements
			dk1[i] = time_step * Qp * (*db.nodes[i]->ddu);
			dK1[i] = time_step * (*db.nodes[i]->copy_du);
			*db.nodes[i]->u = 0.5* dK1[i];
			*db.nodes[i]->du = *db.nodes[i]->copy_du + 0.5 * dk1[i];
			//Rotations
			rk1[i] = time_step * (*db.nodes[i]->domega);
			rK1[i] = time_step * (*db.nodes[i]->copy_omega);
			*db.nodes[i]->alpha = 0.5* rK1[i];
			*db.nodes[i]->omega = *db.nodes[i]->copy_omega + 0.5 * rk1[i];

			//UPDATING 
			for (int k = 0; k < 3; k++)
			{
				if (db.nodes[i]->constraints[k] == 0 && db.nodes[i]->active_GL[k] == 1)
					db.nodes[i]->displacements[k] = (*db.nodes[i]->u)(k, 0);
				else
					(*db.nodes[i]->u)(k, 0) = db.nodes[i]->displacements[k];
				if (db.nodes[i]->constraints[k + 3] == 0 && db.nodes[i]->active_GL[k + 3] == 1)
					db.nodes[i]->displacements[k + 3] = (*db.nodes[i]->alpha)(k, 0);
				else
					(*db.nodes[i]->alpha)(k, 0) = db.nodes[i]->displacements[k + 3];
			}
		}
	}
	//Loop 2
	//Evaluates prescribed displacements
	MountDisplacementsExplicit(time - time_step/2);
	//Evaluates inertial and distributed loads contributions on particles and bodies
	MountExplicit();
	//Evaluates applied loads contributions on global vector of loads
	MountLoadsExplicit(time - time_step/2);
	MontSpecialConstraintsExplicit(time - time_step/2);
	MountContactsExplicit(time - time_step/2);
	MountGlobalExplicit();
	//db.global_P_A.print();
#pragma omp parallel
	{
		//Evaluating accelerations - contributions stemming from particles
#pragma omp for schedule(dynamic)
		for (int i = 0; i < db.number_particles; i++)
			db.particles[i]->EvaluateAccelerations();
	}
	//Evaluating accelerations - contributions stemming from bodies
	//Loop on bodies (TO DO)
	//...

#pragma omp parallel
	{
#pragma omp for schedule(dynamic)
		for (int i = 0; i < db.number_nodes; i++)
		{
			//Rotation variables
			Matrix A;
			double g;
			Matrix Qdelta;
			Matrix Qp;
			//Auxiliary variables
			A = skew(*db.nodes[i]->alpha);
			g = 4.0 / (4.0 + norm(*db.nodes[i]->alpha)*norm(*db.nodes[i]->alpha));
			Qdelta = *I3 + g * (A + 0.5*(A)*(A));
			Qp = (*db.nodes[i]->Q0) * (*db.nodes[i]->Q) * Qdelta;
			//Displacements
			dk2[i] = time_step * Qp * (*db.nodes[i]->ddu);
			dK2[i] = time_step * (*db.nodes[i]->copy_du + 0.5 * dk1[i]);
			*db.nodes[i]->u = 0.5* dK2[i];
			*db.nodes[i]->du = *db.nodes[i]->copy_du + 0.5 * dk2[i];
			//Rotations
			rk2[i] = time_step * (*db.nodes[i]->domega);
			rK2[i] = time_step * db.nodes[i]->InvXi(0.5 * rK1[i]) * (*db.nodes[i]->copy_omega + 0.5 * rk1[i]);
			*db.nodes[i]->alpha = 0.5* rK2[i];
			*db.nodes[i]->omega = *db.nodes[i]->copy_omega + 0.5 * rk2[i];

			//UPDATING 
			for (int k = 0; k < 3; k++)
			{
				if (db.nodes[i]->constraints[k] == 0 && db.nodes[i]->active_GL[k] == 1)
					db.nodes[i]->displacements[k] = (*db.nodes[i]->u)(k, 0);
				else
					(*db.nodes[i]->u)(k, 0) = db.nodes[i]->displacements[k];
				if (db.nodes[i]->constraints[k + 3] == 0 && db.nodes[i]->active_GL[k + 3] == 1)
					db.nodes[i]->displacements[k + 3] = (*db.nodes[i]->alpha)(k, 0);
				else
					(*db.nodes[i]->alpha)(k, 0) = db.nodes[i]->displacements[k + 3];
			}
		}
	}
	//Loop 3
	//Evaluates prescribed displacements
	MountDisplacementsExplicit(time - time_step/2);
	//Evaluates inertial and distributed loads contributions on particles and bodies
	MountExplicit();
	//Evaluates applied loads contributions on global vector of loads
	MountLoadsExplicit(time - time_step/2);
	MontSpecialConstraintsExplicit(time - time_step/2);
	MountContactsExplicit(time - time_step/2);
	MountGlobalExplicit();
	//db.global_P_A.print();
	//Evaluating accelerations - contributions stemming from particles
#pragma omp parallel
	{
#pragma omp for schedule(dynamic)
		for (int i = 0; i < db.number_particles; i++)
			db.particles[i]->EvaluateAccelerations();
	}
	//Evaluating accelerations - contributions stemming from bodies
	//Loop on bodies (TO DO)
	//...
#pragma omp parallel
	{
#pragma omp for schedule(dynamic)
		for (int i = 0; i < db.number_nodes; i++)
		{
			//Rotation variables
			Matrix A;
			double g;
			Matrix Qdelta;
			Matrix Qp;
			//Auxiliary variables
			A = skew(*db.nodes[i]->alpha);
			g = 4.0 / (4.0 + norm(*db.nodes[i]->alpha)*norm(*db.nodes[i]->alpha));
			Qdelta = *I3 + g * (A + 0.5*(A)*(A));
			Qp = (*db.nodes[i]->Q0) * (*db.nodes[i]->Q) * Qdelta;
			//Displacements
			dk3[i] = time_step * Qp * (*db.nodes[i]->ddu);
			dK3[i] = time_step * (*db.nodes[i]->copy_du + 0.5 * dk2[i]);
			*db.nodes[i]->u = dK3[i];
			*db.nodes[i]->du = *db.nodes[i]->copy_du + dk3[i];
			//Rotations
			rk3[i] = time_step * (*db.nodes[i]->domega);
			rK3[i] = time_step * db.nodes[i]->InvXi(0.5 * rK2[i]) * (*db.nodes[i]->copy_omega + 0.5 * rk2[i]);
			*db.nodes[i]->alpha = rK3[i];
			*db.nodes[i]->omega = *db.nodes[i]->copy_omega + rk3[i];

			//UPDATING 
			for (int k = 0; k < 3; k++)
			{
				if (db.nodes[i]->constraints[k] == 0 && db.nodes[i]->active_GL[k] == 1)
					db.nodes[i]->displacements[k] = (*db.nodes[i]->u)(k, 0);
				else
					(*db.nodes[i]->u)(k, 0) = db.nodes[i]->displacements[k];
				if (db.nodes[i]->constraints[k + 3] == 0 && db.nodes[i]->active_GL[k + 3] == 1)
					db.nodes[i]->displacements[k + 3] = (*db.nodes[i]->alpha)(k, 0);
				else
					(*db.nodes[i]->alpha)(k, 0) = db.nodes[i]->displacements[k + 3];
			}
		}
	}
	//Loop 4
	//Evaluates prescribed displacements
	MountDisplacementsExplicit(time);
	//Evaluates inertial and distributed loads contributions on particles and bodies
	MountExplicit();
	//Evaluates applied loads contributions on global vector of loads
	MountLoadsExplicit(time);
	MontSpecialConstraintsExplicit(time);
	MountContactsExplicit(time);
	MountGlobalExplicit();
	//db.global_P_A.print();
#pragma omp parallel
	{
		//Evaluating accelerations - contributions stemming from particles
#pragma omp for schedule(dynamic)
		for (int i = 0; i < db.number_particles; i++)
			db.particles[i]->EvaluateAccelerations();
	}
	//Evaluating accelerations - contributions stemming from bodies
	//Loop on bodies (TO DO)
	//...
#pragma omp parallel
	{
#pragma omp for schedule(dynamic)
		for (int i = 0; i < db.number_nodes; i++)
		{
			//Rotation variables
			Matrix A;
			double g;
			Matrix Qdelta;
			Matrix Qp;
			//Auxiliary variables
			A = skew(*db.nodes[i]->alpha);
			g = 4.0 / (4.0 + norm(*db.nodes[i]->alpha)*norm(*db.nodes[i]->alpha));
			Qdelta = *I3 + g * (A + 0.5*(A)*(A));
			Qp = (*db.nodes[i]->Q0) * (*db.nodes[i]->Q) * Qdelta;
			//Displacements
			dk4[i] = time_step * Qp * (*db.nodes[i]->ddu);
			dK4[i] = time_step * (*db.nodes[i]->copy_du + dk3[i]);
			//Rotations
			rk4[i] = time_step * (*db.nodes[i]->domega);
			rK4[i] = time_step * db.nodes[i]->InvXi(rK3[i]) * (*db.nodes[i]->copy_omega + rk3[i]);
		}
	}

#pragma omp parallel
	{
#pragma omp for schedule(dynamic)
		for (int i = 0; i < db.number_nodes; i++)
		{
			//Final evaluations
			*db.nodes[i]->u = 1.0 / 6 * (dK1[i] + 2.0 * dK2[i] + 2.0 * dK3[i] + dK4[i]);
			*db.nodes[i]->du = *db.nodes[i]->copy_du + 1.0 / 6 * (dk1[i] + 2.0 * dk2[i] + 2.0 * dk3[i] + dk4[i]);
			*db.nodes[i]->alpha = 1.0 / 6 * (rK1[i] + 2.0 * rK2[i] + 2.0 * rK3[i] + rK4[i]);
			*db.nodes[i]->omega = *db.nodes[i]->copy_omega + 1.0 / 6 * (rk1[i] + 2.0 * rk2[i] + 2.0 * rk3[i] + rk4[i]);
		}
	}

#pragma omp parallel
	{
		//Final updates - nodes
#pragma omp for schedule(dynamic)
		for (int i = 0; i < db.number_nodes; i++)
		{
			//Rotation variables
			Matrix A;
			double g;
			Matrix Qdelta;
			Matrix Qp;
			//Auxiliary variables
			A = skew(*db.nodes[i]->alpha);
			g = 4.0 / (4.0 + norm(*db.nodes[i]->alpha)*norm(*db.nodes[i]->alpha));
			Qdelta = *I3 + g * (A + 0.5*(A)*(A));
			Qp = (*db.nodes[i]->Q0) * (*db.nodes[i]->Q) * Qdelta;
			//To global (acceleration)
			(*db.nodes[i]->ddu) = Qp * (*db.nodes[i]->ddu);
			//Updating
			for (int k = 0; k < 3; k++)
			{
				if (db.nodes[i]->constraints[k] == 0 && db.nodes[i]->active_GL[k] == 1)
				{
					db.nodes[i]->displacements[k] = (*db.nodes[i]->u)(k, 0);
					db.nodes[i]->vel[k] = (*db.nodes[i]->du)(k, 0);
					db.nodes[i]->accel[k] = (*db.nodes[i]->ddu)(k, 0);

				}
				else
				{
					(*db.nodes[i]->u)(k, 0) = db.nodes[i]->displacements[k];
					(*db.nodes[i]->du)(k, 0) = db.nodes[i]->vel[k];
					(*db.nodes[i]->ddu)(k, 0) = db.nodes[i]->accel[k];
				}
			}
			//Updating
			for (int k = 0; k < 3; k++)
			{
				if (db.nodes[i]->constraints[k + 3] == 0 && db.nodes[i]->active_GL[k + 3] == 1)
				{
					db.nodes[i]->displacements[k + 3] = (*db.nodes[i]->alpha)(k, 0);
					db.nodes[i]->vel[k + 3] = (*db.nodes[i]->omega)(k, 0);
					db.nodes[i]->accel[k + 3] = (*db.nodes[i]->domega)(k, 0);
				}
				else
				{
					(*db.nodes[i]->alpha)(k, 0) = db.nodes[i]->displacements[k + 3];
					(*db.nodes[i]->omega)(k, 0) = db.nodes[i]->vel[k + 3];
					(*db.nodes[i]->domega)(k, 0) = db.nodes[i]->accel[k + 3];
				}
			}
		}
	}
	

	//Rotation variables
	delete[] rk1;
	delete[] rK1;
	delete[] rk2;
	delete[] rK2;
	delete[] rk3;
	delete[] rK3;
	delete[] rk4;
	delete[] rK4;
	//Displacement variables
	delete[] dk1;
	delete[] dK1;
	delete[] dk2;
	delete[] dK2;
	delete[] dk3;
	delete[] dK3;
	delete[] dk4;
	delete[] dK4;
}



void ExplicitDynamic::MountExplicit()
{
	Clear();	//Clean global contributions	
	high_resolution_clock::time_point t_last = high_resolution_clock::now();

#pragma omp parallel
	{
#pragma omp for
		for (int i = 0; i < db.number_particles; i++)
			db.particles[i]->EvaluateExplicit();
	}

	//Loop on bodies (TO DO)
	//...

	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(high_resolution_clock::now() - t_last).count();
	if (db.plot_times == true)
		cout << "MountLocal duration:\t" << duration / 1e6 << " sec." << "\n";
}

void ExplicitDynamic::MountLoadsExplicit(double t)
{
	for (int i = 0; i < db.number_loads; i++)
		db.loads[i]->EvaluateExplicit(t);
}

void ExplicitDynamic::MountDisplacementsExplicit(double t)
{
	for (int i = 0; i < db.number_displacements; i++)
		db.displacements[i]->EvaluateExplicit(t);
}

void ExplicitDynamic::MountContactsExplicit(double t)
{
	//Paralelização dentro da função SolveContacts
	if (db.gcs_exist)
		if (db.gcs->bool_table.GetAt(db.current_solution_number - 1))
			db.gcs->MountContactsExplicit(t);
}
void ExplicitDynamic::FinalUpdateContactsExplicit(double t)
{
	//Paralelização dentro da função do gcs
	if (db.gcs_exist)
		if (db.gcs->bool_table.GetAt(db.current_solution_number - 1))
			db.gcs->FinalUpdateContactsExplicit(t);
}

void ExplicitDynamic::MontSpecialConstraintsExplicit(double t)
{

}

void ExplicitDynamic::InitialEvaluations()
{
	//Contributions stemming from particles
	for (int i = 0; i < db.number_particles; i++)
		db.particles[i]->InitialEvaluations();

	//Contributions stemming from particles
	for (int i = 0; i < db.number_boundaries; i++)
		db.boundaries[i]->InitialEvaluations();

	//Loop on bodies (TO DO)
	//...
}
void ExplicitDynamic::MountGlobalExplicit()
{
	//GeneralContactSearch
	if (db.gcs_exist)
		if (db.gcs->bool_table.GetAt(db.current_solution_number - 1))
			db.gcs->MountContactsGlobalExplicit();

}
