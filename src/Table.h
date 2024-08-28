#pragma once
#include <stdio.h>

class Table
{
public:
	Table();
	Table(int e_times, int e_values);
	~Table();

	
	int n_times;							//Guarda o numero de instantes salvos (linhas)
	int n_values;							//Numero de valores da tabela (colunas)
	
	void SetTime(int time_index, double time);
	void SetValue(int time_index, int value_index, double value);

	double GetTime(int time_index);
	double GetValue(int time_index, int value_index);
	double GetValueAt(double time_value, int value_index);	//Retorna o valor no instante requerido (interpolação linear)

	void flush();

	bool Read(FILE *f);
	void Write(FILE *f);

protected:
	double** values;						//Guarda os valores em função do tempo
	double* times;							//Guarda os instantes associados aos valores
};