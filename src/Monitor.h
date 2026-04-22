#pragma once
#include <stdio.h>
#include <vector>
#include "MonitorNodesUserDefined.h"
#include "MonitorElementsUserDefined.h"

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
	bool monitor_userdef_exists;
	bool monitor_elemuserdef_exists;

	vector<UserDefMonitorEntry>     userdef_entries;
	vector<UserDefElemMonitorEntry> userdef_elem_entries;

	int sample;

	bool first_record;	//Flag que indica que e a primeira vez que a fun��o UpdateMonitor e chamada (para gravar cabe�alhos nos arquivos)
	bool contact_special_output;
	bool print_times;

	bool Read(FILE *f);
	bool ReadIntTable(FILE *f, vector<int> *data_table);
	bool ReadUserDefEntry(FILE *f);
	bool ReadUserDefParams(FILE *f, UserDefMonitorParams& p);
	bool ReadElemUserDefEntry(FILE *f);
	bool ReadElemUserDefParams(FILE *f, UserDefElemMonitorParams& p);
	static bool ReadElemUserDefParams_B(const char* s, UserDefElemMonitorParams& p);
	void Write(FILE *f);

	void StartMonitor();
	void UpdateMonitor(double time);
	void EndMonitor();

	void UpdateGlobalMonitor(double time);
};

