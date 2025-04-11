#include "Monitor.h"
#include <direct.h>

#include "Node.h"
#include "Element.h"
#include "Contact.h"
#include "NodeSet.h"
#include "Particle.h"
#include "ContactParticleParticle.h"
#include "SurfacePairGeneralContact.h"
#include "SpecialConstraint.h"
#include "GeneralContactSearch.h"
#include "ConvergenceCriteria.h"

#include"Database.h"
#include "IO.h"
//Variaveis globais
extern
Database db;


Monitor::Monitor()
{
	nodes.clear();
	elements.clear();
	contacts.clear();
	node_sets.clear();
	global.clear();
	particles.clear();
	special_constraints.clear();

	f_nodes = NULL;
	f_elements = NULL;
	f_contacts = NULL;
	f_node_sets = NULL;
	f_global = NULL;
	f_particles = NULL;
	f_special_constraints = NULL;
	
	first_record = true;
	sample = 1;

	alloced_files = false;

	contact_special_output = false;
	print_times = true;

	monitor_nodes_exists = false;
	monitor_elements_exists = false;
	monitor_contacts_exists = false;
	monitor_node_sets_exists = false;
	monitor_global_exists = false;
	monitor_particles_exists = false;
	monitor_special_constraints_exists = false;
}

Monitor::~Monitor()
{
	FlushFiles();

	nodes.clear();
	elements.clear();
	contacts.clear();
	node_sets.clear();
	global.clear();
	particles.clear();
	special_constraints.clear();
}

void Monitor::AllocFiles()
{
	FlushFiles();
	f_nodes = new FILE*[(int)nodes.size()];
	f_elements = new FILE*[(int)elements.size()];
	f_contacts = new FILE*[(int)contacts.size()];
	f_node_sets = new FILE*[(int)node_sets.size()];
	f_particles = new FILE*[(int)particles.size()];
	f_global = new FILE*[1];
	f_special_constraints = new FILE*[(int)special_constraints.size()];
	alloced_files = true;
}

void Monitor::FlushFiles()
{
	if (alloced_files == true)
	{
		if (f_nodes != NULL)
			delete[]f_nodes;
		if (f_elements != NULL)
			delete[]f_elements;
		if (f_contacts != NULL)
			delete[]f_contacts;
		if (f_node_sets != NULL)
			delete[]f_node_sets;
		if (f_particles != NULL)
			delete[]f_particles;
		if (f_global != NULL)
			delete[]f_global;
		if (f_special_constraints != NULL)
			delete[]f_special_constraints;
	
		f_nodes = NULL;
		f_elements = NULL;
		f_contacts = NULL;
		f_node_sets = NULL;
		f_particles = NULL;
		f_global = NULL;
		f_special_constraints = NULL;
		
		alloced_files = false;
	}
}

bool Monitor::Read(FILE *f)
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

	bool option_OK = false;
	bool flag_continue = true;
	fpos_t pos;
	while (flag_continue == true)
	{
		option_OK = false;
		TryComment(f);
		fgetpos(f, &pos);
		if (fscanf(f, "%s", s) == EOF)
			return true;
		if (!strcmp(s, "MonitorNodes"))
		{
			if (!ReadIntTable(f,&nodes))
				return false;
			option_OK = true;
			monitor_nodes_exists = true;
		}
		if (!strcmp(s, "MonitorElements"))
		{
			if (!ReadIntTable(f, &elements))
				return false;
			option_OK = true;
			monitor_elements_exists = true;
		}
		if (!strcmp(s, "MonitorContacts"))
		{
			if (!ReadIntTable(f, &contacts))
				return false;
			option_OK = true;
			monitor_contacts_exists = true;
		}
		if (!strcmp(s, "MonitorContactsSpecial"))
		{
			contact_special_output = true;
			if (!ReadIntTable(f, &contacts))
				return false;
			option_OK = true;
			monitor_contacts_exists = true;
		}
		if (!strcmp(s, "MonitorNodeSets"))
		{
			if (!ReadIntTable(f, &node_sets))
				return false;
			option_OK = true;
			monitor_node_sets_exists = true;
		}
		if (!strcmp(s, "MonitorParticles"))
		{
			if (!ReadIntTable(f, &particles))
				return false;
			option_OK = true;
			monitor_particles_exists = true;
		}
		if (!strcmp(s, "MonitorGlobal"))
		{
			option_OK = true;
			monitor_global_exists = true;
		}

		if (!strcmp(s, "MonitorSpecialConstraints"))
		{
			if (!ReadIntTable(f, &special_constraints))
				return false;
			option_OK = true;
			monitor_special_constraints_exists = true;
		}

		if (option_OK == true)
			flag_continue = true;
		else
		{
			fsetpos(f, &pos);
			flag_continue = false;
		}
	}
	return true;
}

bool Monitor::ReadIntTable(FILE *f, vector<int> *data_table)
{
	char s[1000];
	fpos_t pos;
	fgetpos(f, &pos);
	bool flag_not_digit = false;
	fscanf(f, "%s", s);
	while (flag_not_digit == false)
	{
		data_table->push_back((atoi(s)));
		fgetpos(f, &pos);
		if (fscanf(f, "%s", s) == EOF)
			return true;
		if (!isdigit(s[0]))
			flag_not_digit = true;
	}
	fsetpos(f, &pos);
	return true;
}

void Monitor::Write(FILE *f)
{
	fprintf(f, "Monitor\tSample\t%d\n",sample);
	if (monitor_nodes_exists)
	{
		fprintf(f, "MonitorNodes\t");
		for (int i = 0; i < nodes.size(); i++)
			fprintf(f, "%d\t", nodes[i]);
		fprintf(f, "\n");
	}
	if (monitor_elements_exists)
	{
		fprintf(f, "MonitorElements\t");
		for (int i = 0; i < elements.size(); i++)
			fprintf(f, "%d\t", elements[i]);
		fprintf(f, "\n");
	}
	if (monitor_contacts_exists)
	{
		fprintf(f, "MonitorContacts\t");
		for (int i = 0; i < contacts.size(); i++)
			fprintf(f, "%d\t", contacts[i]);
		fprintf(f, "\n");
	}
	if (monitor_node_sets_exists)
	{
		fprintf(f, "MonitorNodeSets\t");
		for (int i = 0; i < node_sets.size(); i++)
			fprintf(f, "%d\t", node_sets[i]);
		fprintf(f, "\n");
	}
	if (monitor_particles_exists)
	{
		fprintf(f, "MonitorParticles\t");
		for (int i = 0; i < particles.size(); i++)
			fprintf(f, "%d\t", particles[i]);
		fprintf(f, "\n");
	}
	if (monitor_global_exists)
	{
		fprintf(f, "MonitorGlobal\t");
		fprintf(f, "\n");
	}

	if (monitor_special_constraints_exists)
	{
		fprintf(f, "MonitorSpecialConstraints\t");
		for (int i = 0; i < special_constraints.size(); i++)
			fprintf(f, "%d\t", special_constraints[i]);
		fprintf(f, "\n");
	}
	
}
void Monitor::StartMonitor()
{
	//Abre os arquivos dos nós
	char name[200];
	for (int i = 0; i < nodes.size(); i++)
	{
		int node = nodes[i];
		strcpy(name, db.folder_name);
		strcat(name, "monitors/");
		_mkdir(name);
		char number[20];
		_itoa(node, number, 10);	//Converte o numero inteiro para char, a fim de escrever o arquivo de resultados com o nome desejado
		strcat(name, "monitor_node_");
		strcat(name, number);
		strcat(name, ".txt");
		f_nodes[i] = fopen(name, "w");
		//fprintf(f_nodes[i], "NODE\t%d\n",node);
	}
	
	//Abre os arquivos dos elementos
	for (int i = 0; i < elements.size(); i++)
	{
		int element = elements[i];
		strcpy(name, db.folder_name);
		strcat(name, "monitors/");
		_mkdir(name);
		char number[20];
		_itoa(element, number, 10);	//Converte o numero inteiro para char, a fim de escrever o arquivo de resultados com o nome desejado
		strcat(name, "monitor_element_");
		strcat(name, number);
		strcat(name, ".txt");
		f_elements[i] = fopen(name, "w");
		//fprintf(f_elements[i], "ELEMENT\t%d\t%s\n", element,db.elements[element-1]->type_name);
	}

	//Abre os arquivos dos contatos
	for (int i = 0; i < contacts.size(); i++)
	{
		int contact = contacts[i];
		strcpy(name, db.folder_name);
		strcat(name, "monitors/");
		_mkdir(name);
		char number[20];
		_itoa(contact, number, 10);	//Converte o numero inteiro para char, a fim de escrever o arquivo de resultados com o nome desejado
		strcat(name, "monitor_contact_");
		strcat(name, number);
		strcat(name, ".txt");
		f_contacts[i] = fopen(name, "w");
		//fprintf(f_contacts[i], "CONTACT\t%d\t%s\n", contact, db.contacts[contact-1]->type_name);
	}

	//Abre os arquivos dos node sets
	for (int i = 0; i < node_sets.size(); i++)
	{
		int set = node_sets[i];
		strcpy(name, db.folder_name);
		strcat(name, "monitors/");
		_mkdir(name);
		char number[20];
		_itoa(set, number, 10);	//Converte o numero inteiro para char, a fim de escrever o arquivo de resultados com o nome desejado
		strcat(name, "monitor_nodeset_");
		strcat(name, number);
		strcat(name, ".txt");
		f_node_sets[i] = fopen(name, "w");
		//fprintf(f_node_sets[i], "NODESET\t%d\n", set);
	}

	//Abre os arquivos das particulas
	for (int i = 0; i < particles.size(); i++)
	{
		int particle = particles[i];
		strcpy(name, db.folder_name);
		strcat(name, "monitors/");
		_mkdir(name);
		char number[20];
		_itoa(particle, number, 10);	//Converte o numero inteiro para char, a fim de escrever o arquivo de resultados com o nome desejado
		strcat(name, "monitor_particle_");
		strcat(name, number);
		strcat(name, ".txt");
		f_particles[i] = fopen(name, "w");
		//fprintf(f_elements[i], "PARTICLE\t%d\t%s\n", element,db.particles[particle-1]->type_name);
	}

	//Abre o arquivo do monitor global
	strcpy(name, db.folder_name);
	strcat(name, "monitors/");
	_mkdir(name);
	strcat(name, "monitor_global");
	strcat(name, ".txt");
	f_global[0] = fopen(name, "w");

	//Abre os arquivos dos special constraints
	for (int i = 0; i < special_constraints.size(); i++)
	{
		int sc = special_constraints[i];
		strcpy(name, db.folder_name);
		strcat(name, "monitors/");
		_mkdir(name);
		char number[20];
		_itoa(sc, number, 10);	//Converte o numero inteiro para char, a fim de escrever o arquivo de resultados com o nome desejado
		strcat(name, "monitor_special_constraint_");
		strcat(name, number);
		strcat(name, ".txt");
		f_special_constraints[i] = fopen(name, "w");
		//fprintf(f_elements[i], "SPECIALCONSTRAINT\t%d\t%s\n", sc,db.special_constraints[sc-1]->type_name);
	}
}

void Monitor::UpdateMonitor(double time)
{
	//Atualizações a serem realizadas - salvas nos arquivos do monitor
	for (int i = 0; i < nodes.size(); i++)
		db.nodes[nodes[i] - 1]->WriteMonitor(f_nodes[i], first_record,time);
	for (int i = 0; i < elements.size(); i++)
		db.elements[elements[i] - 1]->WriteMonitor(f_elements[i], first_record, time);
	for (int i = 0; i < contacts.size(); i++)
		db.contacts[contacts[i] - 1]->WriteMonitor(f_contacts[i], first_record, time);
	for (int i = 0; i < node_sets.size(); i++)
		db.node_sets[node_sets[i] - 1]->WriteMonitor(f_node_sets[i], first_record, time);
	for (int i = 0; i < particles.size(); i++)
		db.particles[particles[i] - 1]->WriteMonitor(f_particles[i], first_record, time);
	UpdateGlobalMonitor(time);
	for (int i = 0; i < special_constraints.size(); i++)
		db.special_constraints[special_constraints[i] - 1]->WriteMonitor(f_special_constraints[i], first_record, time);
	first_record = false;
}

void Monitor::EndMonitor()
{
	for (int i = 0; i < nodes.size(); i++)
		fclose(f_nodes[i]);
	for (int i = 0; i < elements.size(); i++)
		fclose(f_elements[i]);
	for (int i = 0; i < contacts.size(); i++)
		fclose(f_contacts[i]);
	for (int i = 0; i < node_sets.size(); i++)
		fclose(f_node_sets[i]);
	for (int i = 0; i < special_constraints.size(); i++)
		fclose(f_special_constraints[i]);
	for (int i = 0; i < particles.size(); i++)
		fclose(f_particles[i]);
	fclose(f_global[0]);
}

void Monitor::UpdateGlobalMonitor(double time)
{
	//Cabeçalho
	if (first_record == true)
	{
		fprintf(f_global[0], "TIME\t");
		fprintf(f_global[0], "KIN_EN\t");
		fprintf(f_global[0], "STRAIN_EN\t");
		fprintf(f_global[0], "GRAV_EN\t");
		fprintf(f_global[0], "MECH_EN\t");
		/*fprintf(f_global[0], "LIN_MOM_X\t");
		fprintf(f_global[0], "LIN_MOM_Y\t");
		fprintf(f_global[0], "LIN_MOM_Z\t");
		fprintf(f_global[0], "ANG_MOM_X\t");
		fprintf(f_global[0], "ANG_MOM_Y\t");
		fprintf(f_global[0], "ANG_MOM_Z\t");*/
		if (db.gcs_exist)
		{
			fprintf(f_global[0], "FTX\t"); //Marina
			fprintf(f_global[0], "FTY\t"); //Marina
			fprintf(f_global[0], "FTZ\t"); //Marina
			fprintf(f_global[0], "MONITORED_PP\t");
			fprintf(f_global[0], "ACTIVE_PP\t");
			//Marina
			fprintf(f_global[0], "ACTIVE_PP_DEG\t");
			fprintf(f_global[0], "MONITORED_PB\t");
			fprintf(f_global[0], "ACTIVE_PB\t");
			fprintf(f_global[0], "MONITORED_BOBO\t");
			fprintf(f_global[0], "ACTIVE_BOBO\t");
			fprintf(f_global[0], "MONITORED_PBO\t");
			fprintf(f_global[0], "ACTIVE_PBO\t");

			if (print_times)
			{
				fprintf(f_global[0], "DURATION_VERLET\t");
				fprintf(f_global[0], "DURATION_LCELLS\t");
				fprintf(f_global[0], "DURATION_COLDETECTION\t");
				fprintf(f_global[0], "NUMBER_COLDETECTION\t");
				fprintf(f_global[0], "DURATION_MOUNTCONTACTS\t");
			}
		}
		
		fprintf(f_global[0], "TIMESTEP\t");
		fprintf(f_global[0], "F_CRIT\t");
		fprintf(f_global[0], "M_CRIT\t");
		fprintf(f_global[0], "DISP_CRIT\t");
		fprintf(f_global[0], "ROT_CRIT\t");
		fprintf(f_global[0], "\n");
	}
	//Informações a serem salvas
	fprintf(f_global[0], "%.6e\t",time);

	//Energia cinetica, momento linear e momento angular das particulas
	double kin = 0;
	double strain = 0;
	double grav = 0;
	/*Matrix linear_momentum(3);
	Matrix angular_momentum(3);*/
	if (db.particles_exist)
	{
		for (int i = 0; i < db.number_particles; i++)
		{
			kin += db.particles[i]->kinetic_energy;
			strain += db.particles[i]->strain_energy;
			grav += db.particles[i]->potential_g_energy;
			/*linear_momentum(0, 0) += db.particles[i]->linear_momentum[0];
			linear_momentum(1, 0) += db.particles[i]->linear_momentum[1];
			linear_momentum(2, 0) += db.particles[i]->linear_momentum[2];
			angular_momentum(0, 0) += db.particles[i]->angular_momentum_origin[0];
			angular_momentum(1, 0) += db.particles[i]->angular_momentum_origin[1];
			angular_momentum(2, 0) += db.particles[i]->angular_momentum_origin[2];*/
		}
			
	}
	if (db.elements_exist)
	{
		for (int i = 0; i < db.number_elements; i++)
		{
			kin += db.elements[i]->kinetic_energy;
			strain += db.elements[i]->strain_energy;
			grav += db.elements[i]->potential_gravitational_energy;
		}
	}
	fprintf(f_global[0], "%.6e\t", kin);
	fprintf(f_global[0], "%.6e\t", strain);
	fprintf(f_global[0], "%.6e\t", grav);
	fprintf(f_global[0], "%.6e\t", kin + strain + grav);
	/*fprintf(f_global[0], "%.6e\t%.6e\t%.6e\t", linear_momentum(0, 0), linear_momentum(1, 0), linear_momentum(2, 0));
	fprintf(f_global[0], "%.6e\t%.6e\t%.6e\t", angular_momentum(0, 0), angular_momentum(1, 0), angular_momentum(2, 0));*/

	if (db.gcs_exist)
	{
		//Marina
		double ft[3];
		ft[0] = 0;
		ft[1] = 0;
		ft[2] = 0;
		for (int i = 0; i < db.gcs->contactPP_list[0].size(); i++) {
			if (db.gcs->contactPP_list[0][i]->contact_pairs[0]->eligible) {
				ft[0] = db.gcs->contactPP_list[0][i]->contact_pairs[0]->ft[0];
				ft[1] = db.gcs->contactPP_list[0][i]->contact_pairs[0]->ft[1];
				ft[2] = db.gcs->contactPP_list[0][i]->contact_pairs[0]->ft[2];
			}
		}
		//Marina
		fprintf(f_global[0], "%.6e\t%.6e\t%.6e\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t", ft[0], ft[1], ft[2], db.gcs->n_monitoring_PP, db.gcs->n_active_PP, db.gcs->n_active_PP_deg, db.gcs->n_monitoring_PB, db.gcs->n_active_PB,
			db.gcs->n_monitoring_BOBO, db.gcs->n_active_BOBO, db.gcs->n_monitoring_PBO, db.gcs->n_active_PBO);
		if (print_times)
		{
			fprintf(f_global[0], "%d\t", db.time_step_marina); //Marina
			fprintf(f_global[0], "%d\t", (int)db.gcs->duration_verlet);
			fprintf(f_global[0], "%d\t", (int)db.gcs->duration_linkedcells);
			fprintf(f_global[0], "%d\t", (int)db.gcs->duration_collision_detection);
			fprintf(f_global[0], "%d\t", (int)db.gcs->global_n_collisiondetection);
			fprintf(f_global[0], "%d\t", (int)db.gcs->duration_mount_contact);
		}
	}
		
	fprintf(f_global[0], "%.6e\t", db.current_time_step);
	fprintf(f_global[0], "%.6e\t", db.conv_criteria->force_criterion);
	fprintf(f_global[0], "%.6e\t", db.conv_criteria->moment_criterion);
	fprintf(f_global[0], "%.6e\t", db.conv_criteria->disp_criterion);
	fprintf(f_global[0], "%.6e\t", db.conv_criteria->rot_criterion);

	

	//Quebra de linha final
	fprintf(f_global[0], "\n");

	db.myprintf("\nSystem kinetic energy:\t%6e\n", kin);
	db.myprintf("System mechanical energy:\t%6e\n", kin + strain + grav);
}