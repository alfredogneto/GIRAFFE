#include "Static.h"
#include <chrono>
using namespace std::chrono;
#include <iostream>

#include "ConvergenceCriteria.h"
#include "Monitor.h"
#include "ConcomitantSolution.h"
#include "GeneralContactSearch.h"
#include "SolverOptions.h"
#include "ConfigurationSave.h"
#include "PostFiles.h"
#include"Database.h"
//Global variables
extern
Database db;


Static::Static()
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

	file_index = 1;								//Number of the file for save results
}


Static::~Static()
{
}

//Reads input file
bool Static::Read(FILE *f)
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
	return true;
}

//Writes output file
void Static::Write(FILE *f)
{
	fprintf(f, "Static\t%d\nEndTime\t%.6e\nTimeStep\t%.6e\nMaxTimeStep\t%.6e\nMinTimeStep\t%.6e\nMaxIt\t%d\nMinIt\t%d\nConvIncrease\t%d\nIncFactor\t%.6e\nSample\t%d\n",
		solution_number, end_time, i_time_step, max_time_step, min_time_step, max_it, min_it, conv_increase, inc_factor, sample);
}

//Solves solution routine
bool Static::Solve()
{
	bool aborted = false;						//flag for abortar simulation
	int convergence_counter = 0;				//contador of number of times seguidas that there was convergence, for possibilitar the load increment
	int converged_number = 0;					//Contador absoluto of the number of configurações convergidas
	int counter_iterations;						//Contador of the number of iterations
	db.conv_criteria->n_conv_evaluated = 0;		//Zeroing histórico of the contador of convergence criterion
	db.last_converged_time = start_time;		//setting initial (and previously converged time)
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
	ZerosVelAccel();
	//Time variables
	double time_step = i_time_step;
	double time = db.last_converged_time;	//last converged time
	aborted = false;	//Flag that indicates if the case must be abortado
	db.conv_criteria->diverged = false;	//Flag that indicates that divergiu (false)

	if (db.gcs_exist)
		if (db.gcs->bool_table.GetAt(db.current_solution_number - 1))
			db.gcs->SolutionStepInitialCheck();					//General contact initial check

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
		db.myprintf("\nTime: %.9f\t\tTime step: %.9f\n\n", time, time_step);
		//Loop of iterations of the Newton-Raphson
		counter_iterations = 1;
		db.current_iteration_number = counter_iterations;
		Zeros();
		//Imposição of displacements
		MountDisplacements();	//Inclui info of displacements impostos in the matrices/vectors global
		PinballCheck();			//Varre contacts - check of pinball
		BeginStepCheck();
		bool res_converged = false;
		bool GL_converged = false;
		db.conv_criteria->diverged = false;	//Flag that indicates that not divergiu (false)
		//Newton Raphson
		while ((res_converged == false || GL_converged == false) && (db.conv_criteria->diverged == false) && (aborted == false))
		{
			db.current_iteration_number = counter_iterations;
			Clear();
			MountSpecialConstraints();						//Assembly of the special constraints
			MountContacts();								//Assembly of the contacts
			MountLocal();									//Assembly of the elements and particles (information local)
			MountElementLoads();							//Assembly of loadings of field in elements
			MountLoads();
			MountGlobal();									//Espalhamento of the information local in the matrices/vectors global
			db.global_P_A = -1.0*db.global_P_A;				//Inverte the sign of the internal force unbalanced
			//Only in the first iteration inserts the displacements impostos
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
			//Analysis of the facilidade of convergence:
			if (((counter_iterations) <= min_it || convergence_counter >= conv_increase))
			{
				time_step = time_step*inc_factor;
				if (time_step > max_time_step)
					time_step = max_time_step;
				convergence_counter = 0;	//Zeroes the convergence counter
			}
			//Writing in files of results
			if (aborted == false)	//Saves the file, according to amostragem
			{
				if ((converged_number%sample == 0 || time == end_time))
					WriteResults(time);								//saving results
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
	db.myprintf("\nSolution step %d time:\t   %lf sec.\n",solution_number, duration / 1e6);
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
void Static::WriteResults(double time_value)
{
	//Updates files of pós-processamento
	db.post_files->UpdateSinglePartPostFiles(solution_number, time_value, file_index);
	db.post_files->WriteConfigurationResults(solution_number, time_value, file_index);
	file_index++;
}

