#include "Solution.h"
#include"Database.h"

//Variáveis globais
extern
Database db;


Solution::Solution()
{
}

Solution::~Solution()
{
}

//Realiza a numeração dos graus de liberdade globais
void Solution::SetGlobalDOFs()
{
	db.flag_nGL_changed = false;
	
	int GL_fixed = 0;
	int GL_free = 0;

	for (int i = 0; i < db.number_nodes; i++)
	{
		db.nodes[i]->n_GL_free = 0;
		db.nodes[i]->n_GL_fixed = 0;
	}
	//Varredura de nós
	for (int i = 0; i < db.number_nodes; i++)
	{
		for (int j = 0; j < db.number_GLs_node; j++)
		{
			//se o GL for livre e ativo
			if (db.nodes[i]->constraints[j] == 0 && db.nodes[i]->active_GL[j] == 1)
			{
				GL_free++;//global
				db.nodes[i]->GLs[j] = GL_free;//numeração do GL global
				db.nodes[i]->n_GL_free++;//incrementa o número de gls livres do nó
			}
			//se o GL for fixo e ativo
			if (db.nodes[i]->constraints[j] == 1 && db.nodes[i]->active_GL[j] == 1)
			{
				GL_fixed--;//global
				db.nodes[i]->GLs[j] = GL_fixed;//numeração do GL global
				db.nodes[i]->n_GL_fixed++;//incrementa o número de gls fixos do nó
			}
		}
	}

	//Varredura de super nodes
	for (int i = 0; i < db.number_super_nodes; i++)
	{
		//If the super node is alloced (some entity has activated it)
		if (db.super_nodes[i]->alloced == true)
		{
			for (int j = 0; j < db.super_nodes[i]->n_DOFs; j++)
			{
				//se o GL for livre e ativo
				if (db.super_nodes[i]->constraints[j] == false)
				{
					GL_free++;//global
					db.super_nodes[i]->DOFs[j] = GL_free;//numeração do GL global
				}
				else
				{
					GL_fixed--;//global
					db.super_nodes[i]->DOFs[j] = GL_fixed;//numeração do GL global
				}
			}
		}
	}

	//Varredura de special constraints - multiplicadores de lagrange (GLs livres)
	for (int i = 0; i < db.number_special_constraints; i++)
	{
		for (int j = 0; j < db.special_constraints[i]->n_GL; j++)
		{
			if (db.special_constraints[i]->active_lambda[j] == 1)
			{
				GL_free++;//global
				db.special_constraints[i]->GLs[j] = GL_free;//numeração do GL global
			}
			else
				db.special_constraints[i]->GLs[j] = 0;
		}
	}

	if (db.n_GL_free != GL_free)
	{
		db.n_GL_free = GL_free;
		db.flag_nGL_changed = true;
	}
	db.n_GL_fixed = -GL_fixed;
}

//Checa e aponta os graus de liberdade ativos e inativos de acordo com os elementos conectados
void Solution::DOFsActive()
{
	//Desativa todos os graus de liberdade
	for (int i = 0; i < db.number_nodes; i++)
	{
		for (int k = 0; k < db.number_GLs_node; k++)
		{
			db.nodes[i]->GLs[k] = 0;
			db.nodes[i]->active_GL[k] = 0;
			db.nodes[i]->constraints[k] = 0;
		}
			
	}
	//Desativa restrições dos super nodes
	for (int i = 0; i < db.number_super_nodes; i++)
	{
		for (int k = 0; k < db.super_nodes[i]->n_DOFs; k++)
			db.super_nodes[i]->constraints[k] = false;
	}
	int temp_node = 0;
	//Percorre elementos para ativar DOFs de acordo com a necessidade, para cada nó envolvido
	for (int i = 0; i < db.number_elements; i++)
	{
		//Percorre os nós do elemento "i"
		for (int j = 0; j < db.elements[i]->n_nodes; j++)
		{
			temp_node = db.elements[i]->nodes[j];
			for (int k = 0; k < db.number_GLs_node; k++)	//Ativa DOFs
			{
				if (db.elements[i]->DOFs[j][k] == 1)
					db.nodes[temp_node - 1]->active_GL[k] = 1;
			}
		}
	}
	//Percorre partículas para ativar DOFs de acordo com a necessidade, para cada nó envolvido
	for (int i = 0; i < db.number_particles; i++)
	{
		temp_node = db.particles[i]->node;
		if (temp_node != 0)
		{
			for (int k = 0; k < db.number_GLs_node; k++)	//Ativa DOFs
			{
				if (db.particles[i]->DOFs[k] == 1)
					db.nodes[temp_node - 1]->active_GL[k] = 1;
			}
		}
		//No caso de partículas associadas a super nodes, não há necessidade de ativar ou desativar. 
		//No futuro pode ser modificado, introduzindo BoolTable na partícula e ativando-a de acordo com o step, se for o caso
	}
	//Percorre superfícies
	for (int i = 0; i < db.number_surfaces; i++)
	{
		if (db.surfaces[i]->pilot_is_used == true)
		{
			temp_node = db.surfaces[i]->pilot_node;
			for (int k = 0; k < 6; k++)	//Ativa DOFs
				db.nodes[temp_node - 1]->active_GL[k] = 1;
		}

	}

	//Percorre body contact boundaries
	for (int b = 0; b < db.number_body_geometries; b++)
	{
		//Percorre geometries dentro do body contact boundary
		for (int i = 0; i < db.body_geometries[b]->n_items; i++)
		{
			int temp_geometry = db.body_geometries[b]->list_items[i];
			//Percorre nodes do geometry
			for (int j = 0; j < db.geometries[temp_geometry - 1]->n_nodes; j++)
			{
				temp_node = db.geometries[temp_geometry - 1]->nodes[j];
				for (int k = 0; k < db.number_GLs_node; k++)	//Ativa DOFs
				{
					if (db.geometries[temp_geometry - 1]->DOFs[j][k] == 1)
						db.nodes[temp_node - 1]->active_GL[k] = 1;
				}
			}
		}
	}
	
	//Percorre boundaries
	for (int i = 0; i < db.number_boundaries; i++)
	{
		temp_node = db.boundaries[i]->node;
		for (int k = 0; k < 6; k++)	//Ativa DOFs
			db.nodes[temp_node - 1]->active_GL[k] = 1;
	}
	//Percorre special constraints
	for (int i = 0; i < db.number_special_constraints; i++)
		db.special_constraints[i]->ActivateDOFs();
	
	//Percorre Constraints para setar graus de liberdade que serao prescritos na atual solution
	for (int i = 0; i < db.number_constraints; i++)
		db.constraints[i]->Mount();

	//PSY
	if (db.psy_coupling_exist)
		db.psy_coupling->SetConstraints();

	//Percorre carregamentos para setar dependência de graus de liberdade ativos 
	for (int i = 0; i < db.number_loads; i++)
		db.loads[i]->UpdateforSolutionStep();
} 

//Montagem dos elementos (informações locais)
void Solution::MountLocal()
{
	high_resolution_clock::time_point t_last = high_resolution_clock::now();

#pragma omp parallel
	{
#pragma omp for
		for (int i = 0; i < db.number_elements; i++)
			db.elements[i]->Mount();
	}

#pragma omp parallel
	{
#pragma omp for
		for (int i = 0; i < db.number_particles; i++)
			db.particles[i]->Mount();
	}

	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(high_resolution_clock::now() - t_last).count();
	if (db.plot_times == true)
		cout << "MountLocal duration:\t" << duration / 1e6 << " sec." << "\n";
}

//Montagem dos elementos (informações locais)
void Solution::MountElementLoads()
{
	high_resolution_clock::time_point t_last = high_resolution_clock::now();

#pragma omp parallel
	{
#pragma omp for
		for (int i = 0; i < db.number_elements; i++)
			db.elements[i]->MountElementLoads();
	}

	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(high_resolution_clock::now() - t_last).count();
	if (db.plot_times == true)
		cout << "MountLocalLoads duration:\t" << duration / 1e6 << " sec." << "\n";
}

//Montagem das special constraints (informações locais)
void Solution::MountSpecialConstraints()
{
	high_resolution_clock::time_point t_last = high_resolution_clock::now();

#pragma omp parallel
	{
#pragma omp for
		for (int i = 0; i < db.number_special_constraints; i++)
		{
			if (db.special_constraints[i]->bool_table.GetAt(db.current_solution_number - 1))
				db.special_constraints[i]->Mount();
		}
			
	}

	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(high_resolution_clock::now() - t_last).count();
	if (db.plot_times == true)
		cout << "MountSpecialConstraints duration:\t" << duration / 1e6 << " sec." << "\n";
}

//Montagem dos contatos (informações locais)
void Solution::MountContacts()
{
	high_resolution_clock::time_point t_last = high_resolution_clock::now();
	//Montagem das superfícies de contato - FillNodes()
#pragma omp parallel
	{
#pragma omp for
		for (int i = 0; i < db.number_surfaces; i++)
			db.surfaces[i]->FillNodes();
		//Splines
		for (int i = 0; i < db.number_splines; i++)
			for (int j = 0; j < db.splines[i]->size_sp_elements; j++) {
				db.splines[i]->sp_element[j]->FillNodes();
			}
	}
	//Paralelização feita dentro da função Mount
	for (int i = 0; i < db.number_contacts; i++)
	{
		if (db.contacts[i]->bool_table.GetAt(db.current_solution_number - 1))
			db.contacts[i]->Mount();
	}

	//Paralelização dentro da função SolveContacts
	if (db.gcs_exist)
		if (db.gcs->bool_table.GetAt(db.current_solution_number - 1))
			db.gcs->MountContacts();

	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(high_resolution_clock::now() - t_last).count();
	if (db.plot_times == true)
		cout << "MountContacts duration:\t" << duration / 1e6 << " sec." << "\n";
}

//Espalhamento das informações locais nas matrizes/vetores globais
void Solution::MountGlobal()
{
	high_resolution_clock::time_point t_last = high_resolution_clock::now();

	//Montagem da matriz global, passando por cada um dos elementos existentes na malha
	for (int i = 0; i < db.number_elements; i++)
		db.elements[i]->MountGlobal();
	//Montagem da matriz global, passando por cada uma das partículas existentes
	for (int i = 0; i < db.number_particles; i++)
		db.particles[i]->MountGlobal();
	//Para cada par de contato monta contribuições na matriz global
	for (int contact_index = 0; contact_index < db.number_contacts; contact_index++)
		if (db.contacts[contact_index]->bool_table.GetAt(db.current_solution_number - 1))
			db.contacts[contact_index]->MountGlobal();
	//Para cada special constraint monta contribuições na matriz global
	for (int i = 0; i < db.number_special_constraints; i++)	
		if (db.special_constraints[i]->bool_table.GetAt(db.current_solution_number - 1))
			db.special_constraints[i]->MountGlobal();
	 
	//GeneralContactSearch
	if (db.gcs_exist)
		if (db.gcs->bool_table.GetAt(db.current_solution_number - 1))
			db.gcs->MountContactsGlobal();

	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(high_resolution_clock::now() - t_last).count();
	if (db.plot_times == true)
		cout << "MountGlobal duration:\t" << duration / 1e6 << " sec." << "\n";
}

//Inclui info de Displacements
void Solution::MountDisplacements()
{
	for (int i = 0; i < db.number_displacements; i++)
		db.displacements[i]->Mount();

	//Computa deslocamentos prescritos (chute inicial) em GL livres, mas afetados pelas special constraints
	/*for (int i = 0; i < db.number_special_constraints; i++)
		db.special_constraints[i]->ComputeInitialGuessDisplacements();*/

	//PSY
	if (db.psy_coupling_exist)
	{
		if (db.psy_coupling->PSY_bool.GetAt(db.current_solution_number - 1))
		{
			//Forces influencing particles comes already from the last converged time-step - already filling the vector of reaction loads db.global_P_B
			//MountLoads();
			/*{
				for (int i = 0; i < db.number_particles; i++)
				{
					db.particles[i]->Mount();
					db.particles[i]->MountGlobal();
				}
			}*/
			db.psy_coupling->Couple();
		}
			
	}
		
}

//Inclui info de Loads
void Solution::MountLoads()
{
	for (int i = 0; i < db.number_loads; i++)
		db.loads[i]->Mount();
}

//Atualiza os deslocamentos nodais livres e multiplicadores de lagrange
void Solution::UpdateDisps()
{
	//Percorre todos os nós, buscando os graus de liberdade para sua atualização - Newton-Raphson
	for (int i = 0; i < db.number_nodes; i++)
	{
		for (int j = 0; j < db.number_GLs_node; j++)
		{
			//Se o grau de liberdade for livre e ativo
			if (db.nodes[i]->constraints[j] == 0 && db.nodes[i]->active_GL[j] == 1)
				db.nodes[i]->displacements[j] += db.global_P_A(db.nodes[i]->GLs[j] - 1, 0);
		}

	}
	//Percorre todos os super nodes, buscando os graus de liberdade para sua atualização - Newton-Raphson
	for (int i = 0; i < db.number_super_nodes; i++)
	{
		for (int j = 0; j < db.super_nodes[i]->n_displacement_DOFs; j++)
		{
			//Se o grau de liberdade for livre e ativo
			if (db.super_nodes[i]->constraints[j] == 0)
				db.super_nodes[i]->displacements[j] += db.global_P_A(db.super_nodes[i]->DOFs[j] - 1, 0);
		}

	}
	//Percorre os special constraints para atualização dos graus de liberdade referentes aos multiplicadores de lagrange
	for (int i = 0; i < db.number_special_constraints; i++)
	{
		for (int j = 0; j < db.special_constraints[i]->n_GL; j++)
		{
			if (db.special_constraints[i]->active_lambda[j] == 1)
				db.special_constraints[i]->lambda[j] += db.global_P_A(db.special_constraints[i]->GLs[j] - 1, 0);
		}
	}
}

//Salva configuração convergida
void Solution::SaveConfiguration()
{
#pragma omp parallel
	{
#pragma omp for
		for (int i = 0; i < db.number_nodes; i++)
		{
			db.nodes[i]->SaveConfiguration();
		}
	}

#pragma omp parallel
	{
#pragma omp for
		for (int i = 0; i < db.number_super_nodes; i++)
		{
			db.super_nodes[i]->SaveConfiguration();
		}
	}

#pragma omp parallel
	{
#pragma omp for
		//Salvando variáveis de interesse de cada elemento, referentes à descrição lagrangiana atualizada - valores nos pontos de Gauss
		for (int i = 0; i < db.number_elements; i++)
		{
			db.elements[i]->SaveLagrange();
		}
	}

#pragma omp parallel
	{
#pragma omp for
		//Salvando variáveis de interesse de cada contato, referentes à descrição lagrangiana atualizada
		for (int i = 0; i < db.number_contacts; i++)
		{
			if (db.contacts[i]->bool_table.GetAt(db.current_solution_number - 1))
				db.contacts[i]->SaveLagrange();
		}
	}

#pragma omp parallel
	{
#pragma omp for
		//Salvando variáveis de interesse de cada elemento, referentes à descrição lagrangiana atualizada - valores nos pontos de Gauss
		for (int i = 0; i < db.number_surfaces; i++)
		{
			db.surfaces[i]->SaveConfiguration();
		}
		//Spline
		for (int i = 0; i < db.number_splines; i++) {
			db.splines[i]->SaveConfiguration();
			for (int j = 0; j < db.splines[i]->size_sp_elements; j++) {
				db.splines[i]->sp_element[j]->SaveConfiguration();
			}
		}
	}

#pragma omp parallel
	{
#pragma omp for
		//Salvando variáveis de interesse de cada contato, referentes à descrição lagrangiana atualizada
		for (int i = 0; i < db.number_particles; i++)
		{
			db.particles[i]->SaveLagrange();
		}
	}

#pragma omp parallel
	{
#pragma omp for
		//Salvando variáveis de interesse de cada contato, referentes à descrição lagrangiana atualizada
		for (int i = 0; i < db.number_boundaries; i++)
		{
			db.boundaries[i]->SaveLagrange();
		}
	}

#pragma omp parallel
	{
#pragma omp for
		//Salvando variáveis de interesse de cada body_contact_boundary
		for (int i = 0; i < db.number_body_geometries; i++)
		{
			db.body_geometries[i]->SaveLagrange();
		}
	}

#pragma omp parallel
	{
#pragma omp for
		//Salvando variáveis de interesse de cada special constraint
		for (int i = 0; i < db.number_special_constraints; i++)
		{
			if (db.special_constraints[i]->bool_table.GetAt(db.current_solution_number - 1))
				db.special_constraints[i]->SaveLagrange();
		}
	}

	if (db.gcs_exist)
		if (db.gcs->bool_table.GetAt(db.current_solution_number - 1))
			db.gcs->SaveConfiguration();

}

//Restaura a última configuração que convergiu
void Solution::RestoreConfiguration()
{
	//Zera todos os deslocamentos dos nós, para que na próxima iteração, o chute inicial seja esse (Lag. Atualizado)
	for (int i = 0; i < db.number_nodes; i++)
	{
		for (int j = 0; j < 6; j++)
			db.nodes[i]->displacements[j] = 0.0;
		for (int j = 0; j < db.number_GLs_node; j++)
		{
			db.nodes[i]->vel[j] = db.nodes[i]->copy_vel[j];
			db.nodes[i]->accel[j] = db.nodes[i]->copy_accel[j];
		}
		//Explicit
		for (int k = 0; k < 3; k++)
		{
			(*db.nodes[i]->u)(k, 0) = 0.0;
			(*db.nodes[i]->du)(k, 0) = (*db.nodes[i]->copy_du)(k, 0);
			(*db.nodes[i]->ddu)(k, 0) = db.nodes[i]->copy_accel[k];
			(*db.nodes[i]->alpha)(k, 0) = 0.0;
			(*db.nodes[i]->omega)(k, 0) = (*db.nodes[i]->copy_omega)(k, 0);
			(*db.nodes[i]->domega)(k, 0) = db.nodes[i]->copy_accel[k + 3];
		}
		
	}

	//Multiplicadores de lagrange das special constraints - restaurando o último valor convergido
	for (int i = 0; i < db.number_special_constraints; i++)
	{
		for (int j = 0; j < db.special_constraints[i]->n_GL; j++)
			db.special_constraints[i]->lambda[j] = db.special_constraints[i]->copy_lambda[j];
	}
#pragma omp parallel
	{
#pragma omp for
		for (int i = 0; i < db.number_surfaces; i++)
			db.surfaces[i]->FillNodes();
		//Spline
		for (int i = 0; i < db.number_splines; i++)
			for (int j = 0; j < db.splines[i]->size_sp_elements; j++) {
				db.splines[i]->sp_element[j]->FillNodes();
			}
	}
}

//Seta o tamanho estimado da matriz de rigidez
void Solution::SetGlobalSize()
{
	int temp_size_AA = 0;
	int temp_size_BB = 0;
	int temp_size_AB = 0;
	int temp_DOF_free = 0;
	int temp_DOF_fixed = 0;
	int temp_node = 0;
	int temp_element = 0;
	//Estimativa do tamanho da matriz de rigidez
	//Varredura de elementos para preencher nDOFs livres e fixos
	for (int el = 0; el < db.number_elements; el++)
	{
		temp_DOF_free = 0;
		temp_DOF_fixed = 0;
		for (int nd = 0; nd < db.elements[el]->n_nodes; nd++)
		{
			temp_node = db.elements[el]->nodes[nd];
			for (int iGL = 0; iGL < db.number_GLs_node; iGL++)
			{
				if (db.nodes[temp_node - 1]->GLs[iGL] > 0)//free
					temp_DOF_free++;
				if (db.nodes[temp_node - 1]->GLs[iGL] < 0)//fixed
					temp_DOF_fixed++;
				//0 - inativo
			}

		}
		temp_size_AA += temp_DOF_free*temp_DOF_free;
		temp_size_BB += temp_DOF_fixed*temp_DOF_fixed;
		temp_size_AB += temp_DOF_free*temp_DOF_fixed;
	}

	db.size_AA = temp_size_AA;
	db.size_BB = temp_size_BB;
	db.size_AB = temp_size_AB;

	//RIGIDEZ
	db.global_stiffness_AA = SparseMatrix(db.n_GL_free, db.n_GL_free, temp_size_AA);
	db.global_stiffness_BB = SparseMatrix(db.n_GL_fixed, db.n_GL_fixed, temp_size_BB);
	db.global_stiffness_AB = SparseMatrix(db.n_GL_free, db.n_GL_fixed, temp_size_AB);
	db.global_stiffness_BA = SparseMatrix(db.n_GL_fixed, db.n_GL_free, temp_size_AB);

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

	//Vetores para avaliação de critérios de parada
	/*db.global_ABS_P_A.setLines(db.n_GL_free);
	db.global_ABS_P_A.setColumns(1);
	db.global_ABS_P_A.alloc();

	db.global_ABS_P_B.setLines(db.n_GL_fixed);
	db.global_ABS_P_B.setColumns(1);
	db.global_ABS_P_B.alloc();

	db.global_COUNT_ABS_P_A.setLines(db.n_GL_free);
	db.global_COUNT_ABS_P_A.setColumns(1);
	db.global_COUNT_ABS_P_A.alloc();

	db.global_COUNT_ABS_P_B.setLines(db.n_GL_fixed);
	db.global_COUNT_ABS_P_B.setColumns(1);
	db.global_COUNT_ABS_P_B.alloc();*/
}

void Solution::PinballCheck()
{
	if (db.contacts_exist == true)
	{
		high_resolution_clock::time_point t_last = high_resolution_clock::now();
		//Pinball check para contatos - paralelização feita dentro de cada função do Pinball
		for (int cc = 0; cc < db.number_contacts; cc++)
			db.contacts[cc]->PinballCheck();

		auto duration = std::chrono::duration_cast<std::chrono::microseconds>(high_resolution_clock::now() - t_last).count();
		if (db.plot_times == true)
			cout << "PinballSearch duration:\t" << duration / 1e6 << " sec." << "\n";
	}
}

//checagem inicial de início de step para contatos
void Solution::BeginStepCheck()
{
	if (db.contacts_exist == true)
	{
		//Updating bounding boxes of surfaces
		for (int s = 0; s < db.number_surfaces; s++)
			db.surfaces[s]->UpdateBox();
		//Splines
		for (int i = 0; i < db.number_splines; i++)
			for (int j = 0; j < db.splines[i]->size_sp_elements; j++) {
				db.splines[i]->sp_element[j]->UpdateBox();
			}

		high_resolution_clock::time_point t_last = high_resolution_clock::now();
#pragma omp parallel
		{
#pragma omp for
			for (int cc = 0; cc < db.number_contacts; cc++)
			{
				if (db.contacts[cc]->bool_table.GetAt(db.current_solution_number - 1))
					db.contacts[cc]->BeginStepCheck();
			}
		}

//		//Montagem das superfícies de contato - FillNodes()
//#pragma omp parallel
//		{
//#pragma omp for
//			for (int i = 0; i < db.number_surfaces; i++)
//				db.surfaces[i]->FillNodes();
//		}

		auto duration = std::chrono::duration_cast<std::chrono::microseconds>(high_resolution_clock::now() - t_last).count();
		if (db.plot_times == true)
			cout << "PinballSearch duration:\t" << duration / 1e6 << " sec." << "\n";
	}
}

//Monta matriz de massa
void Solution::MountMass()
{
	high_resolution_clock::time_point t_last = high_resolution_clock::now();

#pragma omp parallel
	{
#pragma omp for
		for (int i = 0; i < db.number_elements; i++)
			db.elements[i]->MountMass();
	}

	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(high_resolution_clock::now() - t_last).count();
	if (db.plot_times == true)
		cout << "MountMass duration:\t" << duration / 1e6 << " sec." << "\n";
}

//Monta matriz de amortecimento
void Solution::MountDamping(bool update_rayleigh)
{
	high_resolution_clock::time_point t_last = high_resolution_clock::now();

#pragma omp parallel
	{
#pragma omp for
		for (int i = 0; i < db.number_elements; i++)
			db.elements[i]->MountDamping(update_rayleigh);
	}

	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(high_resolution_clock::now() - t_last).count();
	if (db.plot_times == true)
		cout << "MountDamping duration:\t" << duration / 1e6 << " sec." << "\n";
}

//Montagens - dinamica
void Solution::MountDyn()
{
	high_resolution_clock::time_point t_last = high_resolution_clock::now();

#pragma omp parallel
	{
#pragma omp for
		for (int i = 0; i < db.number_elements; i++)
			db.elements[i]->MountDyn();
	}
	
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(high_resolution_clock::now() - t_last).count();
	if (db.plot_times == true)
		cout << "MountDyn duration:\t" << duration / 1e6 << " sec." << "\n";
}

//Verifica erros que impedem avança da análise, mesmo em caso de convergência
bool Solution::HaveErrors()
{
	//Updating bounding boxes of surfaces
	for (int s = 0; s < db.number_surfaces; s++)
		db.surfaces[s]->UpdateBox();

	for (int i = 0; i < db.number_contacts; i++)
	{
		if (db.contacts[i]->bool_table.GetAt(db.current_solution_number - 1))
		{
			if (db.contacts[i]->HaveErrors() == true)
			{
				db.myprintf("The contact pair number %d has reported errors!\n", db.contacts[i]->number);
				return true;
			}
		}
	}

	if (db.gcs_exist)
	{
		//General Contact Check
		if (db.gcs->bool_table.GetAt(db.current_solution_number - 1))
		{
			if (db.gcs->HaveErrors() == true)
			{
				db.myprintf("General Contact Search has reported errors!\n");
				return true;
			}
		}
	}
	
	
	return false;
}

//Zera variáveis iterativas para chute inicial nulo
void Solution::Zeros()
{
	for (int i = 0; i < db.number_nodes; i++)
	{
		//Zera deslocamentos e rotações para o início da próxima iteração
		for (int j = 0; j < 6; j++)
			db.nodes[i]->displacements[j] = 0.0;
		//Explicit
		db.nodes[i]->u->clear();
		db.nodes[i]->alpha->clear();
	}
	for (int i = 0; i < db.number_super_nodes; i++)
	{
		for (int d = 0; d < db.super_nodes[i]->n_displacement_DOFs; d++)
			db.super_nodes[i]->displacements[d] = 0.0;
	}
}

//Zera velocidades e acelerações para realização de análise estática (seguida da dinâmica)
void Solution::ZerosVelAccel()
{
	for (int i = 0; i < db.number_nodes; i++)
	{
		//Zera deslocamentos e rotações para o início da próxima iteração
		for (int j = 0; j < 6; j++)
		{
			db.nodes[i]->copy_vel[j] = 0.0;
			db.nodes[i]->copy_accel[j] = 0.0;
			db.nodes[i]->vel[j] = 0.0;
			db.nodes[i]->accel[j] = 0.0;
		}
	}
}

//Limpa as matrizes esparsas
void Solution::Clear()
{
	db.global_stiffness_AA.Clear();
	db.global_stiffness_BB.Clear();
	db.global_stiffness_AB.Clear();
	db.global_stiffness_BA.Clear();
	db.global_P_A.clear();
	db.global_P_B.clear();
	db.global_I_A.clear();
	db.global_X_B.clear();

	/*db.global_ABS_P_A.clear();
	db.global_ABS_P_B.clear();
	db.global_COUNT_ABS_P_A.clear();
	db.global_COUNT_ABS_P_B.clear();*/
}

//Monta matrizes esparsas
void Solution::MountSparse()
{
	high_resolution_clock::time_point t_last = high_resolution_clock::now();

	db.global_stiffness_AA.Mount();
	db.global_stiffness_BB.Mount();
	db.global_stiffness_AB.Mount();
	db.global_stiffness_BA.Mount();

	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(high_resolution_clock::now() - t_last).count();
	if (db.plot_times == true)
		cout << "MountSparse duration:\t" << duration / 1e6 << " sec." << "\n";
}

//Computa as condições iniciais nodais impostas
void Solution::ComputeInitialConditions(bool zero_IC)
{
	Zeros();
	if (zero_IC)
	{
		for (int i = 0; i < db.number_nodes; i++)
			db.nodes[i]->ZeroIC();
	}
	for (int i = 0; i < db.number_IC; i++)
		db.IC[i]->ComputeInitialCondition();
	//Computa velocidade e aceleracao devido à imposição de restrições (salva na variavel vel e accel)
	for (int i = 0; i < db.number_special_constraints; i++)
	{
		db.special_constraints[i]->ComputeVelAccel();
	}
	for (int i = 0; i < db.number_nodes; i++)
	{
		db.nodes[i]->SaveConfiguration();
	}
}