#include "Modal.h"
#include <direct.h>
#include <chrono>
#include <iostream>

#include "Node.h"
#include "PostFiles.h"
#include "Element.h"
#include "Particle.h"
#include "VEMPolyhedron.h"

#include "IO.h"
#include"Database.h"
//Variaveis globais
extern
Database db;
#define PI 3.1415926535897932384626433832795
using namespace std::chrono;
using namespace std;

Modal::Modal()
{
	solution_number = 0;
	mode_factor = 1.0;
	number_modes = 0;
	export_matrices = false;
	compute_eigenvectors = false;
	number_frames = 0;

	end_time = 0;
	start_time = 0;

	//Flag que indica se e analise concomitante - default: false
	concomitant_solution = false;									
}

Modal::~Modal()
{
	
}

bool Modal::Read(FILE *f)
{
	char s[1000];
	fscanf(f, "%s", s);
	solution_number = atoi(s);

	fscanf(f, "%s", s);
	if (!strcmp(s, "ExportMatrices"))
	{
		fscanf(f, "%s", s);
		if (atoi(s) != 0 && atoi(s) != 1)//Se não for 0 ou 1
			return false;
		if (atoi(s) == 1)
			export_matrices = true;
		else
			export_matrices = false;
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "NumberModes"))
	{
		fscanf(f, "%s", s);
		number_modes = atoi(s);
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "Tolerance"))
	{
		fscanf(f, "%s", s);
		tolerance = atof(s);
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "ComputeEigenvectors"))
	{
		fscanf(f, "%s", s);
		if (atoi(s) != 0 && atoi(s) != 1)//Se não for 0 ou 1
			return false;
		if (atoi(s) == 1)
			compute_eigenvectors = true;
		else
			compute_eigenvectors = false;
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "NumberFrames"))
	{
		fscanf(f, "%s", s);
		number_frames = atoi(s);
	}
	else
		return false;

	//Alocação do vetor que salvara os autovalores
	eig = Matrix(number_modes, 2);

	return true;
}

void Modal::Write(FILE *f)
{
	fprintf(f, "Modal\t%d\nExportMatrices\t%d\nNumberModes\t%d\nTolerance\t%.6e\nComputeEigenvectors\t%d\nNumberFrames\t%d\n",
		solution_number, export_matrices, number_modes, tolerance, compute_eigenvectors, number_frames);
}

//Creates output folder
void Modal::CreateOutputFolder()
{
	strcpy(copy_name, db.folder_name);	
	strcat(copy_name, "post/");
	_mkdir(copy_name);				//criando diretório post
	//Escrevendo o nome da solução
	char analysis[100];
	sprintf(analysis, "solution_%d", solution_number);
	strcat(copy_name, analysis);
	strcat(copy_name, "/");
	_mkdir(copy_name);			//cria a pasta com o tipo de analise em questão
	strcat(copy_name, analysis);
}

//Realiza analise modal
bool Modal::Solve()
{
	//Plotagem do dia da simulação
	system_clock::time_point today = system_clock::now();
	time_t tt;
	tt = system_clock::to_time_t(today);
	if (concomitant_solution == false)
	{
		CreateOutputFolder();
		db.myprintf("\nSolution step %d started at %s\n", solution_number, ctime(&tt));
	}
	else
		db.myprintf("\nConcomitant solution started at %s\n", ctime(&tt));
	high_resolution_clock::time_point t1 = high_resolution_clock::now();//Inicia a marcação de tempo de execução
	DOFsActive();														//Para cada nó ativa DOFs - tambem opera sobre multiplicadores de Lagrange de SpecialConstraints
	SetGlobalDOFs();													//Numeração de graus de liberdade
	SetGlobalSize();													//Calcula o tamanho da matriz global - com base nos GLs livres e fixos
	Zeros();															//Zera displacements dos nós
	//////////////////////////////////////////Montagem da matriz de massa para analise modal//////////////////////////////////////////////////////
	Clear();															//Limpeza das matrizes globais
	ZerosLocalMatrices();												//Limpeza das matrizes dos elementos (locais)
	MountMassModal();													//Montagem da matriz de massa
	MountDynModal();													//Copia informações da matriz de massa para a rigidez, para espalhamento correto no passo posterior
	MountGlobal();														//Espalhamento das informações locais nas matrizes/vetores globais
	MountSparse();														//Montagem da matriz esparsa
	db.global_mass_AA = db.global_stiffness_AA;						//Cria cópia da matriz de massa (a mesma havia sido salva na variavel referente a matriz de rigidez para se utilizar das mesmas funções de montagem de matrizes globais)
	//////////////////////////////////////////Montagem da matriz de rigidez para analise modal////////////////////////////////////////////////////
	Clear();															//Limpeza das matrizes globais
	ZerosLocalMatrices();												//Limpeza das matrizes dos elementos (locais)
	MountLocal();														//Montagem das matrizes de rigidez (locais)
	MountElementLoads();												//Montagem das matrizes de rigidez de carregamentos nos elementos (locais)
	MountSpecialConstraints();											//Montagem das special constraints
	MountContacts();													//Montagem dos contatos
	MountLoads();														//Montagem de carregamentos gerais
	MountGlobal();														//Espalhamento das informações locais nas matrizes/vetores globais
	MountSparse();														//Montagem da matriz esparsa
	//Nesse ponto possuimos as duas matrizes relevantes (massa e rigidez), salvas para analise modal. A seguir e feita a rotina de extração de autovalores.
	//Se houver requisição, exporta as matrizes para a pasta da solution
	if (export_matrices == true)
		WriteMatrices();
	int ret;
	if (compute_eigenvectors == true)
		eigenvectors = Matrix(db.global_stiffness_AA.rows, number_modes + 1);	//saida dos autovetores - alocação
	eig = sparseeigen(db.global_stiffness_AA, db.global_mass_AA, eigenvectors, number_modes, compute_eigenvectors, &ret, tolerance);
	if (ret != 0)
	{
		if (tolerance != 0)
			db.myprintf("\nEingenvectors extraction was not fulfilled with the desired tolerance requirement, check results carefully!\n");
		//return false;	//error
	}
	if (concomitant_solution == false)
		WriteResults(end_time);

	high_resolution_clock::time_point t2 = high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
	today = system_clock::now();
	tt = system_clock::to_time_t(today);
	if (concomitant_solution == false)
	{
		db.myprintf("\nSolution step %d time:\t   %lf sec.\n", solution_number, duration / 1e6);
		//Plotagem do dia da simulação
		db.myprintf("\nSolution step %d finished at %s\n", solution_number, ctime(&tt));
	}
	else
		db.myprintf("\nConcomitant solution finished at %s\n", ctime(&tt));

	return true;
}

//Escreve resultados
void Modal::WriteResults(double time_value)
{
	char name[1000];
	FILE* f;

	strcpy(name, copy_name);
	strcat(name, "_eigenvalues.txt");
	f = fopen(name, "w");
	//Escrevendo autovalores
	fprintf(f, "Eigenvalues - solution number %d\n", solution_number);
	fprintf(f, "Mode\tRe\tIm\n");
	for (int mode = 0; mode < number_modes; mode++)
	{
		fprintf(f, "%d\t%.6f\t%.6f\n", mode + 1, eig(mode, 0), eig(mode, 1));
	}
	fclose(f);
	
	if (compute_eigenvectors == true)
		WriteVTK_XML();
}

void Modal::WriteMatrices()
{
	char name[1000];
	FILE* f;
	//DOF table
	strcpy(name, copy_name);
	strcat(name, "_DOF_table.txt");
	f = fopen(name, "w");
	WriteDOFTable(f);
	fclose(f);

	//Stiffness matrix
	strcpy(name, copy_name);
	strcat(name, "_stiffness_matrix.txt");
	db.global_stiffness_AA.WriteMatrix(name);
	
	//Stiffness matrix
	strcpy(name, copy_name);
	strcat(name, "_mass_matrix.txt");
	db.global_mass_AA.WriteMatrix(name);
}

//Escreve arquivos com os modos de vibrar
void Modal::WriteModes()
{
	char file_name[400];
	double* re_mode;
	re_mode = new double[db.number_GLs_node];
	double* im_mode;
	im_mode = new double[db.number_GLs_node];
	//Percorrendo os modos para escrita
	bool first = true;	//para casos de autovalor complexo
	for (int mode = 0; mode < number_modes; mode++)
	{
		//Nome do arquivo
		sprintf(file_name, "mode_%d_solution_%d.txt", mode + 1, solution_number);
		FILE *f = fopen(file_name, "w");
		for (int node = 0; node < db.number_nodes; node++)
		{
			//Coordenadas do nó em questão
			fprintf(f, "Node\t%d\tX\t%.6e\tY\t%.6e\tZ\t%.6e\t", node+1,
				db.nodes[node]->copy_coordinates[0],
				db.nodes[node]->copy_coordinates[1], 
				db.nodes[node]->copy_coordinates[2]);
		    //Deslocamentos do nó em questão
			if (eig(mode, 1) == 0.0)//autovalor real
			{
				for (int k = 0; k < db.number_GLs_node; k++)
				{
					if (db.nodes[node]->GLs[k] > 0)
						re_mode[k] = eigenvectors(db.nodes[node]->GLs[k]-1, mode);
					else
						re_mode[k] = 0.0;
					im_mode[k] = 0.0;
				}
			}
			else//autovalor complexo
			{
				if (first == true)
				{
					for (int k = 0; k < db.number_GLs_node; k++)
					{
						if (db.nodes[node]->GLs[k] > 0)
						{
							re_mode[k] = eigenvectors(db.nodes[node]->GLs[k]-1, mode);
							im_mode[k] = eigenvectors(db.nodes[node]->GLs[k]-1, mode + 1);
						}
						else
						{
							re_mode[k] = 0.0;
							im_mode[k] = 0.0;
						}
						
					}
				}
				else
				{
					for (int k = 0; k < db.number_GLs_node; k++)
					{
						if (db.nodes[node]->GLs[k] > 0)
						{
							re_mode[k] = eigenvectors(db.nodes[node]->GLs[k]-1, mode - 1);
							im_mode[k] = -eigenvectors(db.nodes[node]->GLs[k]-1, mode);//complexo conjugado do anterior
						}
						else
						{
							re_mode[k] = 0.0;
							im_mode[k] = 0.0;
						}
					}
				}
			}
			fprintf(f, "Re_UX\t%.6e\tIm_UX\t%.6e\tRe_UY\t%.6e\tIm_UY\t%.6e\tRe_UZ\t%.6e\tIm_UZ\t%.6e\t", re_mode[0], im_mode[0], re_mode[1], im_mode[1], re_mode[2], im_mode[2]);
			if (db.number_GLs_node >= 6)
				fprintf(f, "Re_RX\t%.6e\tIm_RX\t%.6e\tRe_RY\t%.6e\tIm_RY\t%.6e\tRe_RZ\t%.6e\tIm_RZ\t%.6e\t", re_mode[3], im_mode[3], re_mode[4], im_mode[4], re_mode[5], im_mode[5]);
			fprintf(f, "\n");
		}
		//Lógica para complexos conjugados
		if (first == true && eig(mode, 1) != 0.0)
			first = false;
		else
			first = true;
		fclose(f);

	}
	delete[]re_mode;
	delete[]im_mode;
}

//Calcula o valor do deslocamento modal no nó e grau de liberdade em questão
void Modal::ComputeModalDisplacement(int node, int DOF, int mode, double* Re, double* Im)
{
	//Deslocamentos do nó em questão
	if (eig(mode, 1) == 0.0)//autovalor real
	{
		if (db.nodes[node]->GLs[DOF] > 0)
			*Re = eigenvectors(db.nodes[node]->GLs[DOF] - 1, mode);
		else
			*Re = 0.0;
		*Im = 0.0;
	}
	else//autovalor complexo
	{
		bool first = true;
		if (mode > 0)
		if (eig(mode - 1, 1) + eig(mode, 1) == 0.0)
			first = false;
		if (first == true)
		{
			if (db.nodes[node]->GLs[DOF] > 0)
			{
				*Re = eigenvectors(db.nodes[node]->GLs[DOF] - 1, mode);
				*Im = eigenvectors(db.nodes[node]->GLs[DOF] - 1, mode + 1);
			}
			else
			{
				*Re = 0.0;
				*Im = 0.0;
			}
		}
		else
		{
			if (db.nodes[node]->GLs[DOF] > 0)
			{
				*Re = eigenvectors(db.nodes[node]->GLs[DOF] - 1, mode - 1);
				*Im = -eigenvectors(db.nodes[node]->GLs[DOF] - 1, mode);//complexo conjugado do anterior
			}
			else
			{
				*Re = 0.0;
				*Im = 0.0;
			}
		}
	}
}

//Escreve modos de vibrar de todos os elementos
void Modal::WriteVTK_XML()
{
	//Creates a copy of nodal coordinates
	double** copy_copy_coordinates;
	copy_copy_coordinates = new double*[db.number_nodes];
	for (int node = 0; node < db.number_nodes; node++)
		copy_copy_coordinates[node] = new double[6];
	double** copy_ref_coordinates;
	copy_ref_coordinates = new double*[db.number_nodes];
	for (int node = 0; node < db.number_nodes; node++)
		copy_ref_coordinates[node] = new double[6];
	for (int node = 0; node < db.number_nodes; node++)
	{
		for (int i = 0; i < 6; i++)
		{
			copy_copy_coordinates[node][i] = db.nodes[node]->copy_coordinates[i];
			copy_ref_coordinates[node][i] = db.nodes[node]->ref_coordinates[i];
		}
	}

	double bbox = db.EvaluateBoundingBoxDiag() / 8;
	for (int mode = 0; mode < number_modes; mode++)
	{
		//Calcula o mode factor utilizando uma relação com a diagonal do bounding box do modelo
		mode_factor = bbox / ComputeModeNorm(mode);
		for (int c_factor = 0; c_factor < number_frames; c_factor++)
		{
			//Updating nodes
			double Re, Im;
			for (int node = 0; node < db.number_nodes; node++)
			{
				for (int i = 0; i < 6; i++)
				{
					ComputeModalDisplacement(node, i, mode, &Re, &Im);
					//Displacements of the mode
					db.nodes[node]->copy_coordinates[i] = copy_copy_coordinates[node][i] + mode_factor*Re*cos(c_factor * 2 * PI / number_frames);
					//Reference coordinates of the mode
					db.nodes[node]->ref_coordinates[i] = copy_copy_coordinates[node][i];
				}
			}
			db.post_files->UpdateMultiplePartPostFiles(solution_number, c_factor + 1, c_factor + 1, mode + 1);
		}
	}
	//Returning original coordinates to nodes
	for (int node = 0; node < db.number_nodes; node++)
	{
		for (int i = 0; i < 6; i++)
		{
			db.nodes[node]->copy_coordinates[i] = copy_copy_coordinates[node][i];
			db.nodes[node]->ref_coordinates[i] = copy_ref_coordinates[node][i];
		}
	}
	//Cleaning
	for (int node = 0; node < db.number_nodes; node++)
	{
		delete[]copy_copy_coordinates[node];
		delete[]copy_ref_coordinates[node];
	}
	delete[]copy_copy_coordinates;
	delete[]copy_ref_coordinates;
}

//Computa a norma do autovetor
double Modal::ComputeModeNorm(int mode)
{
	double ret_norm = 0;
	double temp = 0;
	//Deslocamentos do nó em questão
	if (eig(mode, 1) == 0.0)//autovalor real
	{
		for (int DOF = 0; DOF<db.global_stiffness_AA.rows; DOF++)
			if (abs(eigenvectors(DOF, mode))>ret_norm)
				ret_norm = abs(eigenvectors(DOF, mode));
	}
	else
	{
		bool first = true;
		if (mode > 0)
		if (eig(mode - 1, 1) + eig(mode, 1) == 0.0)
			first = false;
		if (first == true)
		{
			for (int DOF = 0; DOF<db.global_stiffness_AA.rows; DOF++)
			if (abs(eigenvectors(DOF, mode))>ret_norm)
				ret_norm = abs(eigenvectors(DOF, mode));
		}
		else
		{
			for (int DOF = 0; DOF<db.global_stiffness_AA.rows; DOF++)
			if (abs(eigenvectors(DOF, mode - 1))>ret_norm)
				ret_norm = abs(eigenvectors(DOF, mode - 1));
		}
	}

	return ret_norm;
}

//Monta matriz de massa para analise modal
void Modal::MountMassModal()
{
	high_resolution_clock::time_point t_last = high_resolution_clock::now();

#pragma omp parallel
	{
#pragma omp for
		for (int i = 0; i < db.number_elements; i++)
			db.elements[i]->MountMassModal();
	}

	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(high_resolution_clock::now() - t_last).count();
	if (db.plot_times == true)
		cout << "MountMassModal duration:\t" << duration / 1e6 << " sec." << "\n";
}

//Monta matriz de amortecimento para analise modal
void Modal::MountDampingModal()
{
	high_resolution_clock::time_point t_last = high_resolution_clock::now();

#pragma omp parallel
	{
#pragma omp for
		for (int i = 0; i < db.number_elements; i++)
			db.elements[i]->MountDampingModal();
	}

	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(high_resolution_clock::now() - t_last).count();
	if (db.plot_times == true)
		cout << "MountDampingModal duration:\t" << duration / 1e6 << " sec." << "\n";
}

//Montagens - Modal
void Modal::MountDynModal()
{
	high_resolution_clock::time_point t_last = high_resolution_clock::now();

#pragma omp parallel
	{
#pragma omp for
		for (int i = 0; i < db.number_elements; i++)
			db.elements[i]->MountDynModal();
	}
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(high_resolution_clock::now() - t_last).count();
	if (db.plot_times == true)
		cout << "MountDynModal duration:\t" << duration / 1e6 << " sec." << "\n";

	//Test - VEM modal analysis
	for (int i = 0; i < db.number_particles; i++)
	{
		if (typeid(*db.particles[i]) == typeid(VEMPolyhedron))
		{
			VEMPolyhedron* ptr = static_cast<VEMPolyhedron*>(db.particles[i]);
			ptr->Mount();
			for (int ii=0;ii<ptr->n_vertices*3;ii++)
				for(int jj=0;jj< ptr->n_vertices * 3;jj++)
					ptr->local_stiffness[ii][jj] = ptr->local_mass[ii][jj];
		}
	}
}

//Zera matrizes dos elementos (para uso na analise modal)
void Modal::ZerosLocalMatrices()
{
#pragma omp parallel
	{
#pragma omp for
		for (int i = 0; i < db.number_elements; i++)
			db.elements[i]->Zeros();
	}
}