#pragma once
#include <stdio.h>
#include <vector>

using namespace std;

class Monitor
{
public:
	Monitor();
	~Monitor();
	void AllocFiles();
	void FlushFiles();
	bool alloced_files;

	vector<int> nodes;
	vector<int> elements;
	vector<int> contacts;
	vector<int> node_sets;
	vector<int> particles;
	vector<int> global;
	vector<int> special_constraints;

	FILE **f_nodes;
	FILE **f_elements;
	FILE **f_contacts;
	FILE **f_node_sets;
	FILE **f_global;
	FILE **f_particles;
	FILE **f_special_constraints;

	bool monitor_nodes_exists;
	bool monitor_elements_exists;
	bool monitor_contacts_exists;
	bool monitor_node_sets_exists;
	bool monitor_global_exists;
	bool monitor_particles_exists;
	bool monitor_special_constraints_exists;

	int sample;

	bool first_record;	//Flag que indica que é a primeira vez que a função UpdateMonitor é chamada (para gravar cabeçalhos nos arquivos)
	bool contact_special_output;
	bool print_times;

	bool Read(FILE *f);
	bool ReadIntTable(FILE *f, vector<int> *data_table);
	void Write(FILE *f);

	void StartMonitor();
	void UpdateMonitor(double time);
	void EndMonitor();

	void UpdateGlobalMonitor(double time);
};

