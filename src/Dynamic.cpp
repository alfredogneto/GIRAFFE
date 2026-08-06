#include "Dynamic.h"
#include <chrono>
#include <iostream>

#include "ConcomitantSolution.h"
#include "ConvergenceCriteria.h"
#include "Monitor.h"
#include "GeneralContactSearch.h"
#include "SolverOptions.h"
#include "ConfigurationSave.h"
#include "PostFiles.h"
#include "Node.h"
#include "SuperNode.h"
#include "SpecialConstraint.h"

#include "Database.h"
//Global variables
extern
Database db;

using namespace std;
using namespace std::chrono;


Dynamic::Dynamic()
{
	solution_number = 0;

	start_time = 0.0;
	end_time = 0;
	i_time_step = 0;
	max_time_step = 0;
	min_time_step = 0;
	max_it = 0;
	min_it = 0;
	conv_increase = 0;
	inc_factor = 0;
	sample = 0;

	//Damping
	alpha = 0;
	beta = 0;
	update = 0;	//boolean

	a1 = 0.0;
	a2 = 0.0;
	a3 = 0.0;
	a4 = 0.0;
	a5 = 0.0;
	a6 = 0.0;

	file_index = 1;								//Number of the file for save results

	contact_impact_control = false;
	n_steps_impact = 0;

	zero_IC_flag = false;
}


Dynamic::~Dynamic()
{
}

//Reads input file
bool Dynamic::Read(FILE *f)
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
	if (!strcmp(s, "MaxIt"))
	{
		fscanf(f, "%s", s);
		max_it = atoi(s);
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "MinIt"))
	{
		fscanf(f, "%s", s);
		min_it = atoi(s);
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "ConvIncrease"))
	{
		fscanf(f, "%s", s);
		conv_increase = atoi(s);
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "IncFactor"))
	{
		fscanf(f, "%s", s);
		inc_factor = atof(s);
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

	fscanf(f, "%s", s);
	if (!strcmp(s, "RayleighDamping"))
	{
		fscanf(f, "%s", s);
		if (!strcmp(s, "Alpha"))
		{
			fscanf(f, "%s", s);
			alpha = atof(s);
		}
		else
			return false;

		fscanf(f, "%s", s);
		if (!strcmp(s, "Beta"))
		{
			fscanf(f, "%s", s);
			beta = atof(s);
		}
		else
			return false;

		fscanf(f, "%s", s);
		if (!strcmp(s, "Update"))
		{
			fscanf(f, "%s", s);
			update = atoi(s);
		}
		else
			return false;
	}
	else
		return false;
	fscanf(f, "%s", s);
	if (!strcmp(s, "NewmarkCoefficients"))
	{
		fscanf(f, "%s", s);
		if (!strcmp(s, "Beta"))
		{
			fscanf(f, "%s", s);
			beta_new = atof(s);
		}
		else
			return false;

		fscanf(f, "%s", s);
		if (!strcmp(s, "Gamma"))
		{
			fscanf(f, "%s", s);
			gamma_new = atof(s);
		}
		else
			return false;
	}
	else
		return false;

	//Saves the position (stream)
	fpos_t pos;
	fgetpos(f, &pos);
	fscanf(f, "%s", s);
	if (!strcmp(s, "ZeroIC"))
		zero_IC_flag = true;
	else
		fsetpos(f, &pos);



	//Saves the position (stream)
	//fpos_t pos;
	fgetpos(f, &pos);
	fscanf(f, "%s", s);
	if (!strcmp(s, "ContactImpactControl"))
	{
		contact_impact_control = true;
		fscanf(f, "%s", s);
		if (!strcmp(s, "NumberStepsImpact"))
		{
			fscanf(f, "%s", s);
			n_steps_impact = atoi(s);
		}
		else
			return false;
	}
	else
		fsetpos(f, &pos);

	return true;
	
}

//Writes output file
void Dynamic::Write(FILE *f)
{
	fprintf(f, "Dynamic\t%d\nEndTime\t%.6e\nTimeStep\t%.6e\nMaxTimeStep\t%.6e\nMinTimeStep\t%.6e\nMaxIt\t%d\nMinIt\t%d\nConvIncrease\t%d\nIncFactor\t%.6e\nSample\t%d\nRayleighDamping\tAlpha\t%.6e\tBeta\t%.6e\tUpdate\t%d\nNewmarkCoefficients\tBeta\t%.6e\tGamma\t%.6e\n",
		solution_number, end_time, i_time_step, max_time_step, min_time_step, max_it, min_it, conv_increase, inc_factor, sample, alpha, beta, update, beta_new, gamma_new);
	if (zero_IC_flag)
		fprintf(f, "ZeroIC\n");
	if (contact_impact_control)
		fprintf(f, "ContactImpactControl\tNumberStepsImpact\t%d\n", n_steps_impact);
}

//Solves solution routine
bool Dynamic::Solve()
{
	bool aborted = false;						//flag for abortar simulation
	int convergence_counter = 0;				//contador of number of times seguidas that there was convergence, for possibilitar the load increment
	int converged_number = 0;					//Contador absoluto of the number of configurações convergidas
	int counter_iterations;						//Contador of the number of iterations
	db.conv_criteria->n_conv_evaluated = 0;		//Zeroing histórico of the contador of convergence criterion
	db.last_converged_time = start_time;
	//Updates monitor - only if for the first solution
	if (db.monitor_exist == true && solution_number == 1)
		db.monitor->UpdateMonitor(db.last_converged_time);
	//Updates analysis concomitante
	if (db.concomitant_solution_exist == true)
		db.concomitant_solution->UpdateConcomitantSolution(db.last_converged_time);
	WriteResults(start_time);								//saving results
	high_resolution_clock::time_point t1 = high_resolution_clock::now();//Inicia the marcação of time of execution
	//Plotting of the dia of the simulation
	system_clock::time_point today = system_clock::now();
	time_t tt;
	tt = system_clock::to_time_t(today);
	db.myprintf("\nSolution step %d started at %s\n", solution_number, ctime(&tt));
	DOFsActive();								//For each node active DOFs - tambem opera on Lagrange multipliers of SpecialConstraints
	SetGlobalDOFs();							//Numbering of degrees of freedom
	SetGlobalSize();							//Calculates the size of the matrix global - with base in the GLs free and fixed
	
	//Time variables
	double time_step = i_time_step;
	double time = db.last_converged_time;	//last converged time
	ComputeInitialConditions(zero_IC_flag);				//Computes the initial conditions nodal impostas
	bool first = true;						//Flag for saber when compute Rayleigh Damping
	
	if (db.gcs_exist)
		if (db.gcs->bool_table.GetAt(db.current_solution_number - 1))
			db.gcs->SolutionStepInitialCheck();					//General contact initial check

	aborted = false;	//Flag that indicates if the case must be abortado
	db.conv_criteria->diverged = false;	//Flag that indicates that not divergiu (false)
	////////////////////////////////////////////////////////////////////////////////
	while ((time < end_time || db.conv_criteria->diverged == true) && aborted == false)
	{
		//Setting variables in the database
		db.last_converged_time = time;
		time += time_step;	//Increment of the time step, according to the progressão of the solution
		if (time > 0.99999999*end_time)//If already is very next of the last instant of interest	
		{
			time = end_time;
			time_step = end_time - db.last_converged_time;
		}
		db.current_time_step = time_step;
		//Calculates the coefficients of the Newmark method
		CalculateNewmarkCoeff(time_step);
		db.myprintf("\nTime: %.9f\t\tTime step: %.9f\n\n", time, time_step);
		//Loop of iterations of the Newton-Raphson
		counter_iterations = 1;
		db.current_iteration_number = counter_iterations;
		Zeros();
		//Imposição of displacements
		MountDisplacements();
		PinballCheck();			//Varre contacts - check of pinball
		BeginStepCheck();
		UpdateDyn();		//Performs the first integration to the longo of the time
		bool res_converged = false;
		bool GL_converged = false;
		db.conv_criteria->diverged = false;	//Flag that indicates that not divergiu (false)
		//Newton Raphson
		while ((res_converged == false || GL_converged == false) && (db.conv_criteria->diverged == false) && (aborted == false))
		{
			db.current_iteration_number = counter_iterations;
			Clear();										//Clean global matrices	
			MountSpecialConstraints();						//Assembly of the special constraints
			MountContacts();								//Assembly of the contacts
			MountLocal();									//Assembly of the elements and particles (information local)
			MountElementLoads();							//Assembly of loadings of field in elements
			MountLoads();				
			MountMass();									//Assembly of the mass matrix of the elements and particles
			if (first == true || update == 1)
			{
				MountDamping(true);							//Assembly of the damping matrix of the elements
				first = false;
			}
			else
				MountDamping(false);						//Assembly of the damping matrix of the elements
			MountDyn();										
			MountGlobal();									//Espalhamento of the information local in the matrices/vectors global
			db.global_P_A = -1.0*db.global_P_A;			//Inverte the sign of the internal force unbalanced
			MountSparse();
			//Only in the first iteration inserts the displacements impostos
			high_resolution_clock::time_point t_last = high_resolution_clock::now();
			if (counter_iterations == 1)
				db.global_P_A = db.global_P_A - 1.0*(db.global_stiffness_AB*db.global_X_B);
			auto duration = std::chrono::duration_cast<std::chrono::microseconds>(high_resolution_clock::now() - t_last).count();
			if (db.plot_times == true)
				cout << "ImposeDisplacements duration:\t" << duration / 1e6 << " sec." << "\n\n";
			if (counter_iterations == 1)
				db.conv_criteria->EstablishResidualCriteria();	//Estabelece criteria of parada, with base in the vector of internal forces unbalanced initial of the increment
			res_converged = db.conv_criteria->CheckResidualConvergence();
			int info_fail = 0;
			//If there was convergence of the residual calculated with the ultimos incrementos avaliados, not and necessary evaluate novamente
			if (!(res_converged && GL_converged))
			{
				db.myprintf("It.: %d\n", counter_iterations);
				t_last = high_resolution_clock::now();
				db.global_P_A = sparsesystem(db.global_stiffness_AA, db.global_P_A, &info_fail, db.solver_options->processors, db.solver_options->solver);	//Solves the linear system
				duration = std::chrono::duration_cast<std::chrono::microseconds>(high_resolution_clock::now() - t_last).count();
				if (db.plot_times == true)
					cout << "SparseSystem duration:\t" << duration / 1e6 << " sec." << "\n\n";
				res_converged = false;//Sets fomo false, for obrigar the enter in the while more once for recalcular residual with inc of displacements current
				UpdateDisps();														//Updates the displacements nodal
				UpdateDyn();														//Atualizações - Newmark
				GL_converged = db.conv_criteria->CheckGLConvergence();	//With the GL already atualizados evaluates the convergence criterion
			}
			//Divergence: some of the criteria not converged and atingiu the #max_it or there was falha in the iteration (ex: NaN)
			if ((counter_iterations == max_it && !(res_converged && GL_converged)) || info_fail == 1)
				db.conv_criteria->diverged = true;
			//Test for evaluate the divergence
			if (db.conv_criteria->diverged == true)
				db.conv_criteria->PlotDivergenceReport();
			//Test for evaluate that the time step this very small - ABORT SIMULATION
			if (time_step < min_time_step && time != end_time)
			{
				db.myprintf("\nAborting simulation. Step size is too small!\n\n");
				aborted = true;
			}
			counter_iterations++;
		}//End of iterations
		//Check of errors(impeditivos of convergence)
		if (db.conv_criteria->diverged == false)
			if (HaveErrors() == true)
				db.conv_criteria->diverged = true;
		/////////////////////////////////SAVING CONFIGURATION CONVERGED/////////////////////////////////
		if (db.conv_criteria->diverged == false)
		{
			converged_number++;
			SaveConfiguration();	//Saves configuration converged
			convergence_counter++;	//Increments the number of convergências seguidas
			
			//Writing in files of results
			if (aborted == false)	//Saves the file, according to amostragem
			{
				if ((converged_number%sample == 0 || time == end_time))
				{
					WriteResults(time);								//saving results
				}
				//Updates monitor
				if (db.monitor_exist == true && (converged_number%db.monitor->sample == 0 || time == end_time))
					db.monitor->UpdateMonitor(time);
				//Updates analysis concomitante
				if (db.concomitant_solution_exist == true && (converged_number%db.concomitant_solution->sample == 0 || time == end_time))
					db.concomitant_solution->UpdateConcomitantSolution(time);
				//Updates configuration save
				if (db.config_save_exist == true && (converged_number%db.config_save->sample == 0 || time == end_time))
					db.config_save->ExportConfiguration(time);
			}
			
			//Analysis of the facilidade of convergence:
			if (((counter_iterations) <= min_it || convergence_counter >= conv_increase))
			{
				time_step = time_step * inc_factor;
				if (time_step > max_time_step)
					time_step = max_time_step;
				convergence_counter = 0;	//Zeroes the convergence counter
			}
			//Contact impact control
			if (contact_impact_control)
			{
				if (db.gcs_exist)
				{
					//Analyzing limitation on time steps due to contact-impact well-resolution
					double time_step_cimpact = db.gcs->TimeStepControl();
					if (time_step_cimpact < time_step)
					{
						if (time_step_cimpact > 10.0*min_time_step)
						{
							time_step = time_step_cimpact;
							db.myprintf("\nTime step determined by Contact-impact control: %.9f\n", time_step);
						}
						else
						{
							time_step = 10.0*min_time_step;
							db.myprintf("\nTime step close to minimum. Check results carefully: %.9f\n", time_step);
						}
						
					}
						
				}
			}
		}
		/////////////////////////////////RETOMANDO CONFIGURATION CONVERGED/////////////////////////////////
		else
		{
			RestoreConfiguration();	//Restaura the last configuration that converged
			//Modificação of the loading factor increment:
			time -= time_step;
			time_step = time_step / 2.0;	 //Bissecção
			db.conv_criteria->diverged = false;
			convergence_counter = 0;		//Zeroes the convergence counter
		}
		high_resolution_clock::time_point t2 = high_resolution_clock::now();
		auto duration = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
		db.myprintf("\nElapsed time:\t   %lf sec.\n", duration / 1e6);
	}//End of sub steps
	high_resolution_clock::time_point t2 = high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
	db.myprintf("\nSolution step %d time:\t   %lf sec.\n", solution_number, duration / 1e6);
	//Plotting of the dia of the simulation
	today = system_clock::now();
	tt = system_clock::to_time_t(today);
	db.myprintf("\nSolution step %d finished at %s\n", solution_number, ctime(&tt));

	if (aborted == true)
		return false;//Error
	else
		return true;//Not error
}

//Writes results
void Dynamic::WriteResults(double time_value)
{
	//Updates files of pós-processamento
	db.post_files->UpdateSinglePartPostFiles(solution_number, time_value, file_index);
	db.post_files->WriteConfigurationResults(solution_number, time_value, file_index);
	file_index++;
}

//Atualizações - Newmark
void Dynamic::UpdateDyn()
{
	Matrix alpha_delta(3);
	Matrix A_delta;
	double alpha_escalar_delta;
	double g;
	Matrix I(3, 3);		//Identidade of order 3
	I(0, 0) = 1.0;
	I(1, 1) = 1.0;
	I(2, 2) = 1.0;
	Matrix Q_delta(3, 3);
	//Auxiliary variables
	Matrix vel_aux(3, 1);
	Matrix ace_aux(3, 1);
	for (int i = 0; i < db.number_nodes; i++)
	{
		//Translações
		for (int j = 0; j<3; j++)
		{
			/////////////VELOCITIES/////////////////
			//If the GL for free:
			if (db.nodes[i]->GLs[j] > 0)
				db.nodes[i]->vel[j] = db.nodes[i]->displacements[j] * a4 +
				db.nodes[i]->copy_vel[j] * a5 +
				db.nodes[i]->copy_accel[j] * a6;
			/////////////ACELERAÇÕES/////////////////
			//If the GL for free:
			if (db.nodes[i]->GLs[j] > 0)
				db.nodes[i]->accel[j] = db.nodes[i]->displacements[j] * a1 -
				db.nodes[i]->copy_vel[j] * a2 -
				db.nodes[i]->copy_accel[j] * a3;
		}
		//Rotations
		//Calculations referring to to the tensor rotation Q delta - for perform the Newmark in the same espaço tangent (paper of the Ibrahimbegovic, 2000)
		alpha_delta(0, 0) = db.nodes[i]->displacements[3];
		alpha_delta(1, 0) = db.nodes[i]->displacements[4];
		alpha_delta(2, 0) = db.nodes[i]->displacements[5];
		alpha_escalar_delta = norm(alpha_delta);						//Value scalar of the parameter alpha
		A_delta = skew(alpha_delta);									//Matrix A
		g = 4.0 / (4.0 + alpha_escalar_delta*alpha_escalar_delta);		//function g(alpha) - in some ref. tb. called of h(alpha)
		Q_delta = I + g*(A_delta + 0.5*(A_delta*A_delta));				//Tensor of rotation

		for (int j = 3; j<6; j++)
		{
			/////////////VELOCITIES/////////////////
			//If the GL for free:
			if (db.nodes[i]->GLs[j] > 0)
				vel_aux(j - 3, 0) = db.nodes[i]->displacements[j] * a4 +
				db.nodes[i]->copy_vel[j] * a5 +
				db.nodes[i]->copy_accel[j] * a6;
			/////////////ACELERAÇÕES/////////////////
			//If the GL for free:
			if (db.nodes[i]->GLs[j] > 0)
				ace_aux(j - 3, 0) = db.nodes[i]->displacements[j] * a1 -
				db.nodes[i]->copy_vel[j] * a2 -
				db.nodes[i]->copy_accel[j] * a3;
		}

		vel_aux = Q_delta*vel_aux;
		ace_aux = Q_delta*ace_aux;

		for (int j = 3; j<6; j++)
		{
			/////////////VELOCITIES/////////////////
			//If the GL for free:
			if (db.nodes[i]->GLs[j] > 0)
				db.nodes[i]->vel[j] = vel_aux(j - 3, 0);
			/////////////ACELERAÇÕES/////////////////
			//If the GL for free:
			if (db.nodes[i]->GLs[j] > 0)
				db.nodes[i]->accel[j] = ace_aux(j - 3, 0);
		}
	}

	//Super nodes
	for (int i = 0; i < db.number_super_nodes; i++)
	{
		for (int gl = 0; gl < db.super_nodes[i]->n_displacement_DOFs; gl++)
		{
			//Free DOF
			if (db.super_nodes[i]->DOFs[gl] > 0)
			{
				/////////////VELOCITIES/////////////////
				db.super_nodes[i]->vel[gl] = db.super_nodes[i]->displacements[gl] * a4 +
					db.super_nodes[i]->copy_vel[gl] * a5 +
					db.super_nodes[i]->copy_accel[gl] * a6;
					
				/////////////ACELERAÇÕES/////////////////
				db.super_nodes[i]->accel[gl] = db.super_nodes[i]->displacements[gl] * a1 -
					db.super_nodes[i]->copy_vel[gl] * a2 -
					db.super_nodes[i]->copy_accel[gl] * a3;
			}
		}
		
	}
	//TO INVESTIGATE!!!
	//Computes velocity and acceleration due the imposição of constraints
	for (int i = 0; i < db.number_special_constraints; i++)
		db.special_constraints[i]->ComputeVelAccel();

}

void Dynamic::CalculateNewmarkCoeff(double time_step)
{
	a1 = 1.0 / (time_step*time_step*beta_new);
	a2 = 1.0 / (time_step*beta_new);
	a3 = 1.0 / (2.0*beta_new) - 1.0;
	a4 = gamma_new / (time_step*beta_new);
	a5 = 1.0 - gamma_new / beta_new;
	a6 = time_step*(1.0 - gamma_new / (2.0*beta_new));
}
