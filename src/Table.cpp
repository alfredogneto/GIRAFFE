#include "Table.h"
#include <string>

Table::Table()
{
	values = NULL;
	times = NULL;
}

Table::Table(int e_times, int e_values)
{
	values = NULL;
	times = NULL;
	flush();
	n_times = e_times;
	n_values = e_values;

	times = new double[n_times];

	//Alocação do vetor de values
	values = new double *[n_times];
	for (int i = 0; i < n_times; i++)
	{
		values[i] = new double[n_values];
	}

	//Zera todas as componentes dos vetores recém-alocados
	for (int i = 0; i < n_times; i++)
	{
		times[i] = 0;
		for (int j = 0; j < n_values; j++)
				values[i][j] = 0.0;
	}
}

void Table::flush()
{
	if (times != NULL)
		delete [] times;
	if (values != NULL)
	{
		for (int i = 0; i < n_times; i++)
			delete [] values[i];
		delete [] values;
	}
}

Table::~Table()
{
	flush();
}

void Table::SetTime(int time_index, double time)
{
	times[time_index] = time;
}

void Table::SetValue(int time_index, int value_index, double value)
{
	values[time_index][value_index] = value;
}

double Table::GetTime(int time_index)
{
	return times[time_index];
}
double Table::GetValue(int time_index, int value_index)
{
	return values[time_index][value_index];
}

//Retorna o valor no instante requerido (interpolação linear)
double Table::GetValueAt(double time_value, int value_index)
{
	//Realiza a interpolação linear
	double return_value = 0;

	if (n_times == 0)
		return return_value;
	//Se for maior que o último time, retorna o último
	if (time_value >= times[n_times - 1])
		return values[n_times - 1][value_index];
	//Se for menor que o primeiro time, retorna o primeiro
	if (time_value <= times[0])
		return values[0][value_index];

	//Interpolação linear
	bool scape = false;
	int index = 0;
	while (scape == false)
	{
		if (time_value >= times[index])
			index++;
		else
			scape = true;
	}
	double factor = (time_value - times[index - 1]) / (times[index] - times[index - 1]);

	return_value = values[index - 1][value_index] + factor * (values[index][value_index] - values[index - 1][value_index]);

	return return_value;
}

bool  Table::Read(FILE *f)
{
	char s[1000];
	//Leitura da tabela
	for (int i = 0; i < n_times; i++)
	{
		fscanf(f, "%s", s);
		if (!strcmp(s, "Time"))//palavra chave opcional ("Time")
			fscanf(f, "%s", s);
		times[i] = atof(s);
		for (int j = 0; j < n_values; j++)
		{
			fscanf(f, "%s", s);
			values[i][j] = atof(s);
		}
	}
	return true;
}

void  Table::Write(FILE *f)
{
	//Escrita da tabela
	for (int i = 0; i < n_times; i++)
	{
		fprintf(f, "%.6e", times[i]);
		for (int j = 0; j < n_values; j++)
			fprintf(f, "\t%.6e", values[i][j]);
		fprintf(f, "\n");
	}
}