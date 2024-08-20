#include "Dynamic.h"
#include"Database.h"
//Vari�veis globais
extern
Database db;

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
	update = 0;	//booleana

	a1 = 0.0;
	a2 = 0.0;
	a3 = 0.0;
	a4 = 0.0;
	a5 = 0.0;
	a6 = 0.0;

	file_index = 1;								//N�mero do arquivo para salvar resultados

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

	//Salva a posi��o (stream)
	fpos_t pos;
	fgetpos(f, &pos);
	fscanf(f, "%s", s);
	if (!strcmp(s, "ZeroIC"))
		zero_IC_flag = true;
	else
		fsetpos(f, &pos);



	//Salva a posi��o (stream)
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
	bool aborted = false;						//flag para abortar simula��o
	int convergence_counter = 0;				//contador de n�mero de vezes seguidas que houve converg�ncia, para possibilitar o incremento de carga
	int converged_number = 0;					//Contador absoluto do n�mero de configura��es convergidas
	int counter_iterations;						//Contador do numero de itera�oes
	db.conv_criteria->n_conv_evaluated = 0;		//Zerando hist�rico do contador de crit�rio de converg�ncia
	db.last_converged_time = start_time;
	//Atualiza monitor - somente se for o primeira solution
	if (db.monitor_exist == true && solution_number == 1)
		db.monitor->UpdateMonitor(db.last_converged_time);
	//Atualiza an�lise concomitante
	if (db.concomitant_solution_exist == true)
		db.concomitant_solution->UpdateConcomitantSolution(db.last_converged_time);
	WriteResults(start_time);								//salvando resultados
	high_resolution_clock::time_point t1 = high_resolution_clock::now();//Inicia a marca��o de tempo de execu��o
	//Plotagem do dia da simula��o
	system_clock::time_point today = system_clock::now();
	time_t tt;
	tt = system_clock::to_time_t(today);
	db.myprintf("\nSolution step %d started at %s\n", solution_number, ctime(&tt));
	DOFsActive();								//Para cada n� ativa DOFs - tamb�m opera sobre multiplicadores de Lagrange de SpecialConstraints
	SetGlobalDOFs();							//Numera��o de graus de liberdade
	SetGlobalSize();							//Calcula o tamanho da matriz global - com base nos GLs livres e fixos
	
	//Time variables
	double time_step = i_time_step;
	double time = db.last_converged_time;	//last converged time
	ComputeInitialConditions(zero_IC_flag);				//Computa as condi��es iniciais nodais impostas
	bool first = true;						//Flag para saber quando computar Rayleigh Damping
	
	if (db.gcs_exist)
		if (db.gcs->bool_table.GetAt(db.current_solution_number - 1))
			db.gcs->SolutionStepInitialCheck();					//General contact initial check

	aborted = false;	//Flag que indica se o caso deve ser abortado
	db.conv_criteria->diverged = false;	//Flag que indica que n�o divergiu (false)
	////////////////////////////////////////////////////////////////////////////////
	while ((time < end_time || db.conv_criteria->diverged == true) && aborted == false)
	{
		//Setando vari�veis no database
		db.last_converged_time = time;
		time += time_step;	//Incremento do time step, de acordo com a progress�o da solu��o
		if (time > 0.99999999*end_time)//Se j� estiver muito pr�ximo do �ltimo instante de interesse	
		{
			time = end_time;
			time_step = end_time - db.last_converged_time;
		}
		db.current_time_step = time_step;
		//Calcula os coeficientes do M�todo de Newmark
		CalculateNewmarkCoeff(time_step);
		db.myprintf("\nTime: %.9f\t\tTime step: %.9f\n\n", time, time_step);
		//Loop de itera��es do Newton-Raphson
		counter_iterations = 1;
		db.current_iteration_number = counter_iterations;
		Zeros();
		//Imposi��o de deslocamentos
		MountDisplacements();
		PinballCheck();			//Varre contatos - verifica��o de pinball
		BeginStepCheck();
		UpdateDyn();		//Realiza a primeira integra��o ao longo do tempo
		bool res_converged = false;
		bool GL_converged = false;
		db.conv_criteria->diverged = false;	//Flag que indica que n�o divergiu (false)
		//Newton Raphson
		while ((res_converged == false || GL_converged == false) && (db.conv_criteria->diverged == false) && (aborted == false))
		{
			db.current_iteration_number = counter_iterations;
			Clear();										//Clean global matrices	
			MountSpecialConstraints();						//Montagem das special constraints
			MountContacts();								//Montagem dos contatos
			MountLocal();									//Montagem dos elementos e part�culas (informa��es locais)
			MountElementLoads();							//Montagem de carregamentos de campo em elementos
			MountLoads();				
			MountMass();									//Montagem da matriz de massa dos elementos e part�culas
			if (first == true || update == 1)
			{
				MountDamping(true);							//Montagem da matriz de amortecimento dos elementos
				first = false;
			}
			else
				MountDamping(false);						//Montagem da matriz de amortecimento dos elementos
			MountDyn();										
			MountGlobal();									//Espalhamento das informa��es locais nas matrizes/vetores globais
			db.global_P_A = -1.0*db.global_P_A;			//Inverte o sinal do esfor�o desbalanceado
			MountSparse();
			//Apenas na primeira itera��o insere os deslocamentos impostos
			high_resolution_clock::time_point t_last = high_resolution_clock::now();
			if (counter_iterations == 1)
				db.global_P_A = db.global_P_A - 1.0*(db.global_stiffness_AB*db.global_X_B);
			auto duration = std::chrono::duration_cast<std::chrono::microseconds>(high_resolution_clock::now() - t_last).count();
			if (db.plot_times == true)
				cout << "ImposeDisplacements duration:\t" << duration / 1e6 << " sec." << "\n\n";
			if (counter_iterations == 1)
				db.conv_criteria->EstablishResidualCriteria();	//Estabelece crit�rios de parada, com base no vetor de esfor�os desbalanceados iniciais do incremento
			res_converged = db.conv_criteria->CheckResidualConvergence();
			int info_fail = 0;
			//Se houve converg�ncia do res�duo calculado com os �ltimos incrementos avaliados, n�o � necess�rio avaliar novamente
			if (!(res_converged && GL_converged))
			{
				db.myprintf("It.: %d\n", counter_iterations);
				t_last = high_resolution_clock::now();
				db.global_P_A = sparsesystem(db.global_stiffness_AA, db.global_P_A, &info_fail, db.solver_options->processors, db.solver_options->solver);	//Resolve o sistema linear
				duration = std::chrono::duration_cast<std::chrono::microseconds>(high_resolution_clock::now() - t_last).count();
				if (db.plot_times == true)
					cout << "SparseSystem duration:\t" << duration / 1e6 << " sec." << "\n\n";
				res_converged = false;//Seta fomo false, para obrigar a entrar no while mais uma vez para recalcular res�duo com inc de deslocamentos atual
				UpdateDisps();														//Atualiza os deslocamentos nodais
				UpdateDyn();														//Atualiza��es - Newmark
				GL_converged = db.conv_criteria->CheckGLConvergence();	//Com os GL j� atualizados avalia o crit�rio de converg�ncia
			}
			//Diverg�ncia: algum dos crit�rios n�o convergiu e atingiu a #max_it ou houve falha na itera��o (ex: NaN)
			if ((counter_iterations == max_it && !(res_converged && GL_converged)) || info_fail == 1)
				db.conv_criteria->diverged = true;
			//Teste para avaliar a diverg�ncia
			if (db.conv_criteria->diverged == true)
				db.conv_criteria->PlotDivergenceReport();
			//Teste para avaliar que o time step est� muito pequeno - ABORT SIMULATION
			if (time_step < min_time_step && time != end_time)
			{
				db.myprintf("\nAborting simulation. Step size is too small!\n\n");
				aborted = true;
			}
			counter_iterations++;
		}//Fim de itera��es
		//Verifica��o de erros(impeditivos de converg�ncia)
		if (db.conv_criteria->diverged == false)
			if (HaveErrors() == true)
				db.conv_criteria->diverged = true;
		/////////////////////////////////SALVANDO CONFIGURA��O CONVERGIDA/////////////////////////////////
		if (db.conv_criteria->diverged == false)
		{
			converged_number++;
			SaveConfiguration();	//Salva configura��o convergida
			convergence_counter++;	//Incrementa o n�mero de converg�ncias seguidas
			
			//Escrita em arquivos de resultados
			if (aborted == false)	//Salva o arquivo, de acordo com amostragem
			{
				if ((converged_number%sample == 0 || time == end_time))
				{
					WriteResults(time);								//salvando resultados
				}
				//Atualiza monitor
				if (db.monitor_exist == true && (converged_number%db.monitor->sample == 0 || time == end_time))
					db.monitor->UpdateMonitor(time);
				//Atualiza an�lise concomitante
				if (db.concomitant_solution_exist == true && (converged_number%db.concomitant_solution->sample == 0 || time == end_time))
					db.concomitant_solution->UpdateConcomitantSolution(time);
				//Atualiza configuration save
				if (db.config_save_exist == true && (converged_number%db.config_save->sample == 0 || time == end_time))
					db.config_save->ExportConfiguration(time);
			}
			
			//An�lise da facilidade de converg�ncia:
			if (((counter_iterations) <= min_it || convergence_counter >= conv_increase))
			{
				time_step = time_step * inc_factor;
				if (time_step > max_time_step)
					time_step = max_time_step;
				convergence_counter = 0;	//Zera o convergence counter
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
		/////////////////////////////////RETOMANDO CONFIGURA��O CONVERGIDA/////////////////////////////////
		else
		{
			RestoreConfiguration();	//Restaura a �ltima configura��o que convergiu
			//Modifica��o do loading factor increment:
			time -= time_step;
			time_step = time_step / 2.0;	 //Bissec��o
			db.conv_criteria->diverged = false;
			convergence_counter = 0;		//Zera o convergence counter
		}
		high_resolution_clock::time_point t2 = high_resolution_clock::now();
		auto duration = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
		db.myprintf("\nElapsed time:\t   %lf sec.\n", duration / 1e6);
	}//Fim de sub steps
	high_resolution_clock::time_point t2 = high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
	db.myprintf("\nSolution step %d time:\t   %lf sec.\n", solution_number, duration / 1e6);
	//Plotagem do dia da simula��o
	today = system_clock::now();
	tt = system_clock::to_time_t(today);
	db.myprintf("\nSolution step %d finished at %s\n", solution_number, ctime(&tt));

	if (aborted == true)
		return false;//Erro
	else
		return true;//N�o erro
}

//Escreve resultados
void Dynamic::WriteResults(double time_value)
{
	//Atualiza arquivos de p�s-processamento
	db.post_files->UpdateSinglePartPostFiles(solution_number, time_value, file_index);
	db.post_files->WriteConfigurationResults(solution_number, time_value, file_index);
	file_index++;
}

//Atualiza��es - Newmark
void Dynamic::UpdateDyn()
{
	Matrix alpha_delta(3);
	Matrix A_delta;
	double alpha_escalar_delta;
	double g;
	Matrix I(3, 3);		//Identidade de ordem 3
	I(0, 0) = 1.0;
	I(1, 1) = 1.0;
	I(2, 2) = 1.0;
	Matrix Q_delta(3, 3);
	//Vari�veis auxiliares
	Matrix vel_aux(3, 1);
	Matrix ace_aux(3, 1);
	for (int i = 0; i < db.number_nodes; i++)
	{
		//Transla��es
		for (int j = 0; j<3; j++)
		{
			/////////////VELOCIDADES/////////////////
			//Se o GL for livre:
			if (db.nodes[i]->GLs[j] > 0)
				db.nodes[i]->vel[j] = db.nodes[i]->displacements[j] * a4 +
				db.nodes[i]->copy_vel[j] * a5 +
				db.nodes[i]->copy_accel[j] * a6;
			/////////////ACELERA��ES/////////////////
			//Se o GL for livre:
			if (db.nodes[i]->GLs[j] > 0)
				db.nodes[i]->accel[j] = db.nodes[i]->displacements[j] * a1 -
				db.nodes[i]->copy_vel[j] * a2 -
				db.nodes[i]->copy_accel[j] * a3;
		}
		//Rota��es
		//C�lculos referentes ao tensor rota��o Q delta - para realizar o Newmark no mesmo espa�o tangente (paper do Ibrahimbegovic, 2000)
		alpha_delta(0, 0) = db.nodes[i]->displacements[3];
		alpha_delta(1, 0) = db.nodes[i]->displacements[4];
		alpha_delta(2, 0) = db.nodes[i]->displacements[5];
		alpha_escalar_delta = norm(alpha_delta);						//Valor escalar do parametro alpha
		A_delta = skew(alpha_delta);									//Matriz A
		g = 4.0 / (4.0 + alpha_escalar_delta*alpha_escalar_delta);		//fun��o g(alpha) - em algumas ref. tb. chamado de h(alpha)
		Q_delta = I + g*(A_delta + 0.5*(A_delta*A_delta));				//Tensor de rota��o

		for (int j = 3; j<6; j++)
		{
			/////////////VELOCIDADES/////////////////
			//Se o GL for livre:
			if (db.nodes[i]->GLs[j] > 0)
				vel_aux(j - 3, 0) = db.nodes[i]->displacements[j] * a4 +
				db.nodes[i]->copy_vel[j] * a5 +
				db.nodes[i]->copy_accel[j] * a6;
			/////////////ACELERA��ES/////////////////
			//Se o GL for livre:
			if (db.nodes[i]->GLs[j] > 0)
				ace_aux(j - 3, 0) = db.nodes[i]->displacements[j] * a1 -
				db.nodes[i]->copy_vel[j] * a2 -
				db.nodes[i]->copy_accel[j] * a3;
		}

		vel_aux = Q_delta*vel_aux;
		ace_aux = Q_delta*ace_aux;

		for (int j = 3; j<6; j++)
		{
			/////////////VELOCIDADES/////////////////
			//Se o GL for livre:
			if (db.nodes[i]->GLs[j] > 0)
				db.nodes[i]->vel[j] = vel_aux(j - 3, 0);
			/////////////ACELERA��ES/////////////////
			//Se o GL for livre:
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
				/////////////VELOCIDADES/////////////////
				db.super_nodes[i]->vel[gl] = db.super_nodes[i]->displacements[gl] * a4 +
					db.super_nodes[i]->copy_vel[gl] * a5 +
					db.super_nodes[i]->copy_accel[gl] * a6;
					
				/////////////ACELERA��ES/////////////////
				db.super_nodes[i]->accel[gl] = db.super_nodes[i]->displacements[gl] * a1 -
					db.super_nodes[i]->copy_vel[gl] * a2 -
					db.super_nodes[i]->copy_accel[gl] * a3;
			}
		}
		
	}
	//TO INVESTIGATE!!!
	//Computa velocidade e aceleracao devido � imposi��o de restri��es
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
