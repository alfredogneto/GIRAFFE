#include<stdio.h>
#include<math.h>
#include"IO.h"
#include "Errors.h"
#include "string.h"

//#define _CRTDBG_MAP_ALLOC
//#include <stdlib.h>
//#include <crtdbg.h>
//#ifdef _DEBUG
//#define DBG_NEW new ( _NORMAL_BLOCK , __FILE__ , __LINE__ )
//// Replace _NORMAL_BLOCK with _CLIENT_BLOCK if you want the
//// allocations to be of _CLIENT_BLOCK type
//#else
//#define DBG_NEW new
//#endif

//Variáveis globais
Database db;
IO io;									//Criação de objeto IO para entrada e saída de dados

int main(int argc, char* argv[])
{
	_setmaxstdio(2000);
	//_CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF);
	Errors errors;							//Criação de objeto Errors para checar inconsistencia de dados
	bool readOK = io.ReadFile(argc,argv);	//Lê arquivo de entrada
	bool checkOK = errors.CheckErrors();	//Checa erros
	if (readOK == true && checkOK == true)
	{
		db.myprintf("GIRAFFE simulation output report. Version %s.\nFile name: %s\n\n", db.version, db.file_name);
		db.PreCalc();						//Realiza pré-cálculos
		//Execution time
		high_resolution_clock::time_point t1 = high_resolution_clock::now();
		//Start Monitor
		if (db.monitor_exist == true)
			db.monitor->StartMonitor();
		//Start Concomitant Solution
		if (db.concomitant_solution_exist == true)
			db.concomitant_solution->StartConcomitantSolution();
		//Solution
		if (db.solution_exist == true)
		{
			db.post_files->StartPostFiles(db.number_solutions + 1);	//solução combinada - start
			//Percorre soluções sequenciais
			bool converged = true;
			int i = 0;
			while (i < db.number_solutions && converged == true)
			{
				db.post_files->StartPostFiles(i + 1);
				db.current_solution_number = i + 1;
				converged = db.solution[i]->Solve();
				db.post_files->EndPostFiles(i + 1);
				i++;
			}
			db.post_files->EndPostFiles(db.number_solutions + 1);	//solução combinada - end
		}
		//End Monitor
		if (db.monitor_exist == true)
			db.monitor->EndMonitor();
		//End Concomitant Solution
		if (db.concomitant_solution_exist == true)
			db.concomitant_solution->EndConcomitantSolution();
		//Execution time
		high_resolution_clock::time_point t2 = high_resolution_clock::now();
		auto duration = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
		db.myprintf("\nTotal solution time:\t   %lf sec.\n", duration / 1e6);
		//Escreve arquivo de saída
		io.WriteFile();					
	}
	cout << "\nGiraffe execution has finished.\n"; 
	if (argc <= 1)
		system("pause");
	return 0;
	
}