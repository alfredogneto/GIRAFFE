#include "ConfigurationSave.h"

#include <direct.h>

#include "Node.h"
#include "Boundary.h"
#include "Particle.h"
#include "Polyhedron.h"

#include"Database.h"
#include"IO.h"
//Variaveis globais
extern Database db;
extern IO io;

ConfigurationSave::ConfigurationSave()
{
	acumulated_config_number = 1;
}

ConfigurationSave::~ConfigurationSave()
{
}

bool ConfigurationSave::Check()
{
	return true;
}

bool ConfigurationSave::Read(FILE *f)
{
	char s[1000];
	//Verifica a palavra chave "Sample"
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
 
void ConfigurationSave::Write(FILE *f)
{
	fprintf(f, "ConfigurationSave\tSample\t%d\n",sample);
}

void ConfigurationSave::ExportConfiguration(double time)
{
	char name[1000];
	char time_char[100];
	struct stat info;
	char config_name[100];
	sprintf(config_name, "configuration_%d_time_", acumulated_config_number);
	sprintf(time_char, "%lf", time);	//Converte o time para char
	strcat(config_name, time_char);
	
	strcpy(name, db.folder_name);	//pasta do job
	strcat(name, "config/");		//diretório de configurações salvas
	if (stat(name, &info) != 0)		//checa existência do diretório packs
		_mkdir(name);				//criando diretório packs
	strcat(name, config_name);		//nome do pack
	strcat(name, ".inp");			//criando arquivo
	FILE* f;
	f = fopen(name, "w");

	fprintf(f, "//Configuration saved automatically by GIRAFFE");
	//This function saves all the nodes in updated location and:
	//Particles with updated CS - giving current orientation
	//Boundaries with updated CS - giving current orientation

	//Nodes
	if (db.nodes_exist == true)
	{
		fprintf(f, "\nNodes\t%d\n", db.number_nodes);
		for (int i = 0; i < db.number_nodes; i++)
			db.nodes[i]->WriteUpdated(f);			//Escrita dos dados 
	}
	//Geração de boundaries
	if (db.boundaries_exist == true)
	{
		fprintf(f, "\nBoundaries\t%d\n", db.number_boundaries);
		for (int i = 0; i < db.number_boundaries; i++)
		{

			db.boundaries[i]->WriteModifyingParameters(f, db.boundaries[i]->material, db.boundaries[i]->node, db.boundaries[i]->number, db.boundaries[i]->number + 1);					//Escrita dos dados 
		}
	}
	//Geração de particulas
	if (db.particles_exist == true)
	{
		fprintf(f, "\nParticles\t%d\n", db.number_particles);
		for (int i = 0; i < db.number_particles; i++)
		{
			db.particles[i]->WriteModifyingParameters(f, db.particles[i]->material, db.particles[i]->node, db.particles[i]->number, db.particles[i]->number+1 + db.number_boundaries);	//Escrita dos dados 
		}
	}
	//Geração de Coordinate Systems
	fprintf(f, "\nCoordinateSystems\t%d\n", db.number_particles + db.number_boundaries + 1);
	//1 - global
	fprintf(f, "CS\t%d\tE1\t%.12e\t%.12e\t%.12e\tE3\t%.12e\t%.12e\t%.12e\n", 1, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0);
	//2 - boundaries
	for (int i = 0; i < db.number_boundaries; i++)
	{
		Matrix e1(3);
		e1(0, 0) = 1.0;
		Matrix e3(3);
		e3(2, 0) = 1.0;
		e1 = (*db.boundaries[i]->Qip)*e1;
		e3 = (*db.boundaries[i]->Qip)*e3;
		fprintf(f, "CS\t%d\tE1\t%.12e\t%.12e\t%.12e\tE3\t%.12e\t%.12e\t%.12e\n", i + 2, e1(0, 0), e1(1, 0), e1(2, 0), e3(0, 0), e3(1, 0), e3(2, 0));
	}
	//3 - particles
	for (int i = 0; i < db.number_particles; i++)
	{
		if (typeid(*db.particles[i]) == typeid(Polyhedron))
		{
			Polyhedron* ptr = static_cast<Polyhedron*>(db.particles[i]);
			Matrix e1(3);
			e1(0, 0) = 1.0;
			Matrix e3(3);
			e3(2, 0) = 1.0;
			if (db.particles[i]->node != 0)
			{
				e1 = (*ptr->Qip)*e1;
				e3 = (*ptr->Qip)*e3;
			}
			fprintf(f, "CS\t%d\tE1\t%.12e\t%.12e\t%.12e\tE3\t%.12e\t%.12e\t%.12e\n", i + 2 + db.number_boundaries, e1(0, 0), e1(1, 0), e1(2, 0), e3(0, 0), e3(1, 0), e3(2, 0));
		}
	}
	//Condições iniciais
	if (db.nodes_exist == true)
	{
		fprintf(f, "\nInitialConditions\t%d\n", db.number_nodes);
		for (int i = 0; i < db.number_nodes; i++)
		{
			fprintf(f, "InitialCondition\t%d\t", i + 1);
			db.nodes[i]->WriteInitialConditions(f);			//Escrita dos dados 
			fprintf(f, "SolutionStep\t%d\n", 2);
		}
	}

	//Demais itens
	io.WriteConvergenceCriteria(f);
	if (db.materials_exist)
		io.WriteMaterials(f);
	if (db.cad_data_exist)
		io.WriteCADData(f);
	if (db.node_sets_exist)
		io.WriteNodeSets(f);
	if (db.solution_exist)
		io.WriteSolutions(f);
	io.WritePostFiles(f);
	if (db.solver_options_exist)
		io.WriteSolverOptions(f);
	if (db.gcs_exist)
		io.WriteGeneralContactSearch(f);
	if (db.config_save_exist)
		io.WriteConfigurationSave(f);
	if (db.monitor_exist)
		io.WriteMonitor(f);
	if (db.contactinterfaces_exist)
		io.WriteContactInterfaces(f);
	if (db.constraints_exist)
		io.WriteConstraints(f);
	if (db.environment_exist)
		io.WriteEnvironment(f);
	

	fclose(f);
	acumulated_config_number++;
}
