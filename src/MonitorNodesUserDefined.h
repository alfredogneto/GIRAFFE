#pragma once
#include <vector>
#include <stdio.h>
using namespace std;

struct UserDefMonitorParams {
	bool X   = false;
	bool Y   = false;
	bool Z   = false;
	bool dX  = false;
	bool dY  = false;
	bool dZ  = false;
	bool ddX = false;
	bool ddY = false;
	bool ddZ = false;
	bool FX  = false;
	bool FY  = false;
	bool FZ  = false;
	bool MX  = false;
	bool MY  = false;
	bool MZ  = false;
};

struct UserDefMonitorEntry {
	bool          use_node_set = false;	// true = ids s�o NodeSet IDs; false = node IDs
	vector<int>   ids;					// IDs como lidos do input (n�s ou NodeSets)
	UserDefMonitorParams params;
	double        time_start = 0.0;		// tempo a partir do qual os dados s�o gravados (opcional)
	vector<int>              resolved_node_ids;	// IDs de n�s resolvidos - preenchido em StartMonitor()
	vector<UserDefMonitorParams> node_params;	// params por n�, paralelo a resolved_node_ids
	vector<double>           node_time_starts;	// time_start por n�, paralelo a resolved_node_ids
	vector<bool>             node_first_record;	// flag de primeiro registro por n�
	vector<FILE*>            files;				// um FILE* por n� resolvido
};
