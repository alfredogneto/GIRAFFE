#include "Static.h"
#include"Database.h"
//Variáveis globais
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

	file_index = 1;								//Número do arquivo para salvar resultados
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
	bool aborted = false;						//flag para abortar simulação
	int convergence_counter = 0;				//contador de número de vezes seguidas que houve convergência, para possibilitar o incremento de carga
	int converged_number = 0;					//Contador absoluto do número de configurações convergidas
	int counter_iterations;						//Contador do numero de iteraçoes
	db.conv_criteria->n_conv_evaluated = 0;		//Zerando histórico do contador de critério de convergência
	db.last_converged_time = start_time;		//setting initial (and previously converged time)
	//Atualiza monitor - somente se for o primeira solution
	if (db.monitor_exist == true && solution_number == 1)
		db.monitor->UpdateMonitor(db.last_converged_time);
	//Atualiza análise concomitante
	if (db.concomitant_solution_exist == true)
		db.concomitant_solution->UpdateConcomitantSolution(db.last_converged_time);
	WriteResults(start_time);								//salvando resultados
	high_resolution_clock::time_point t1 = high_resolution_clock::now();//Inicia a marcação de tempo de execução
	//Plotagem do dia da simulação
	system_clock::time_point today = system_clock::now();
	time_t tt;
	tt = system_clock::to_time_t(today);
	db.myprintf("\nSolution step %d started at %s\n", solution_number, ctime(&tt));
	DOFsActive();								//Para cada nó ativa DOFs - também opera sobre multiplicadores de Lagrange de SpecialConstraints
	SetGlobalDOFs();							//Numeração de graus de liberdade
	SetGlobalSize();							//Calcula o tamanho da matriz global - com base nos GLs livres e fixos
	ZerosVelAccel();
	//Time variables
	double time_step = i_time_step;
	double time = db.last_converged_time;	//last converged time
	aborted = false;	//Flag que indica se o caso deve ser abortado
	db.conv_criteria->diverged = false;	//Flag que indica que divergiu (false)

	if (db.gcs_exist)
		if (db.gcs->bool_table.GetAt(db.current_solution_number - 1))
			db.gcs->SolutionStepInitialCheck();					//General contact initial check

	////////////////////////////////////////////////////////////////////////////////
	while ((time < end_time || db.conv_criteria->diverged == true) && aborted == false)
	{
		//Setando variáveis no database
		db.last_converged_time = time;
		time += time_step;	//Incremento do time step, de acordo com a progressão da solução
		if (time > 0.99999999*end_time)//Se já estiver muito próximo do último instante de interesse	
		{
			time = end_time;
			time_step = end_time - db.last_converged_time;
		}
		db.current_time_step = time_step;
		db.myprintf("\nTime: %.9f\t\tTime step: %.9f\n\n", time, time_step);
		//Loop de iterações do Newton-Raphson
		counter_iterations = 1;
		db.current_iteration_number = counter_iterations;
		Zeros();
		//Imposição de deslocamentos
		MountDisplacements();	//Inclui info de deslocamentos impostos nas matrizes/vetores globais
		PinballCheck();			//Varre contatos - verificação de pinball
		BeginStepCheck();
		bool res_converged = false;
		bool GL_converged = false;
		db.conv_criteria->diverged = false;	//Flag que indica que não divergiu (false)
		//Newton Raphson
		while ((res_converged == false || GL_converged == false) && (db.conv_criteria->diverged == false) && (aborted == false))
		{
			db.current_iteration_number = counter_iterations;
			Clear();
			MountSpecialConstraints();						//Montagem das special constraints
			MountContacts();								//Montagem dos contatos
			MountLocal();									//Montagem dos elementos e partículas (informações locais)
			MountElementLoads();							//Montagem de carregamentos de campo em elementos
			MountLoads();
			MountGlobal();									//Espalhamento das informações locais nas matrizes/vetores globais
			db.global_P_A = -1.0*db.global_P_A;				//Inverte o sinal do esforço desbalanceado
			//Apenas na primeira iteração insere os deslocamentos impostos
			MountSparse();
			//Apenas na primeira iteração insere os deslocamentos impostos
			high_resolution_clock::time_point t_last = high_resolution_clock::now();
			if (counter_iterations == 1)
				db.global_P_A = db.global_P_A - 1.0*(db.global_stiffness_AB*db.global_X_B);
			auto duration = std::chrono::duration_cast<std::chrono::microseconds>(high_resolution_clock::now() - t_last).count();
			if (db.plot_times == true)
				cout << "ImposeDisplacements duration:\t" << duration / 1e6 << " sec." << "\n\n";
			if (counter_iterations == 1)
				db.conv_criteria->EstablishResidualCriteria();	//Estabelece critérios de parada, com base no vetor de esforços desbalanceados iniciais do incremento
			res_converged = db.conv_criteria->CheckResidualConvergence();
			int info_fail = 0;
			//Se houve convergência do resíduo calculado com os últimos incrementos avaliados, não é necessário avaliar novamente
			if (!(res_converged && GL_converged))
			{
				db.myprintf("It.: %d\n", counter_iterations);
				t_last = high_resolution_clock::now();
				db.global_P_A = sparsesystem(db.global_stiffness_AA, db.global_P_A, &info_fail, db.solver_options->processors, db.solver_options->solver);	//Resolve o sistema linear
				duration = std::chrono::duration_cast<std::chrono::microseconds>(high_resolution_clock::now() - t_last).count();
				if (db.plot_times == true)
					cout << "SparseSystem duration:\t" << duration / 1e6 << " sec." << "\n\n";
				res_converged = false;//Seta fomo false, para obrigar a entrar no while mais uma vez para recalcular resíduo com inc de deslocamentos atual
				UpdateDisps();														//Atualiza os deslocamentos nodais
				GL_converged = db.conv_criteria->CheckGLConvergence();	//Com os GL já atualizados avalia o critério de convergência
			}
			//Divergência: algum dos critérios não convergiu e atingiu a #max_it ou houve falha na iteração (ex: NaN)
			if ((counter_iterations == max_it && !(res_converged && GL_converged)) || info_fail == 1)
				db.conv_criteria->diverged = true;
			//Teste para avaliar a divergência
			if (db.conv_criteria->diverged == true)
				db.conv_criteria->PlotDivergenceReport();
			//Teste para avaliar que o time step está muito pequeno - ABORT SIMULATION
			if (time_step < min_time_step && time != end_time)
			{
				db.myprintf("\nAborting simulation. Step size is too small!\n\n");
				aborted = true;
			}
			counter_iterations++;
		}//Fim de iterações
		//Verificação de erros(impeditivos de convergência)
		if (db.conv_criteria->diverged == false)
			if (HaveErrors() == true)
				db.conv_criteria->diverged = true;
		/////////////////////////////////SALVANDO CONFIGURAÇÃO CONVERGIDA/////////////////////////////////
		if (db.conv_criteria->diverged == false)
		{
			converged_number++;
			SaveConfiguration();	//Salva configuração convergida
			convergence_counter++;	//Incrementa o número de convergências seguidas
			//Análise da facilidade de convergência:
			if (((counter_iterations) <= min_it || convergence_counter >= conv_increase))
			{
				time_step = time_step*inc_factor;
				if (time_step > max_time_step)
					time_step = max_time_step;
				convergence_counter = 0;	//Zera o convergence counter
			}
			//Escrita em arquivos de resultados
			if (aborted == false)	//Salva o arquivo, de acordo com amostragem
			{
				if ((converged_number%sample == 0 || time == end_time))
					WriteResults(time);								//salvando resultados
				//Atualiza monitor
				if (db.monitor_exist == true && (converged_number%db.monitor->sample == 0 || time == end_time))
					db.monitor->UpdateMonitor(time);
				//Atualiza análise concomitante
				if (db.concomitant_solution_exist == true && (converged_number%db.concomitant_solution->sample == 0 || time == end_time))
					db.concomitant_solution->UpdateConcomitantSolution(time);
				//Atualiza configuration save
				if (db.config_save_exist == true && (converged_number%db.config_save->sample == 0 || time == end_time))
					db.config_save->ExportConfiguration(time);
			}
		}
		/////////////////////////////////RETOMANDO CONFIGURAÇÃO CONVERGIDA/////////////////////////////////
		else
		{
			RestoreConfiguration();	//Restaura a última configuração que convergiu
			//Modificação do loading factor increment:
			time -= time_step;
			time_step = time_step / 2.0;	 //Bissecção
			db.conv_criteria->diverged = false;
			convergence_counter = 0;		//Zera o convergence counter
		}
		high_resolution_clock::time_point t2 = high_resolution_clock::now();
		auto duration = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
		db.myprintf("\nElapsed time:\t   %lf sec.\n", duration / 1e6);
	}//Fim de sub steps
	high_resolution_clock::time_point t2 = high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
	db.myprintf("\nSolution step %d time:\t   %lf sec.\n",solution_number, duration / 1e6);
	//Plotagem do dia da simulação
	today = system_clock::now();
	tt = system_clock::to_time_t(today);
	db.myprintf("\nSolution step %d finished at %s\n", solution_number, ctime(&tt));

	if (aborted == true)
		return false;//Erro
	else
		return true;//Não erro
}

//Escreve resultados
void Static::WriteResults(double time_value)
{
	//Atualiza arquivos de pós-processamento
	db.post_files->UpdateSinglePartPostFiles(solution_number, time_value, file_index);
	db.post_files->WriteConfigurationResults(solution_number, time_value, file_index);
	file_index++;
}

