#include "BoolTable.h"
#include <stdio.h>

#include"Database.h"

//Variáveis globais
extern
Database db;

BoolTable::BoolTable()
{
	SetDefault(true);
}
BoolTable::BoolTable(BoolTable &copied)
{
	data_table.clear();
	data_table.resize(copied.data_table.size());
	for (int i = 0; i < copied.data_table.size(); i++)
		data_table[i] = copied.data_table[i];
}

//Reads input file
bool BoolTable::Read(FILE *f)
{
	data_table.clear();
	char s[1000];
	fpos_t pos;
	fgetpos(f, &pos);
	bool flag_not_digit = false;
	fscanf(f, "%s", s);
	while (flag_not_digit == false)
	{
		if (atoi(s) == 1)
			data_table.push_back(true);
		else
			data_table.push_back(false);
		fgetpos(f, &pos);
		if (fscanf(f, "%s", s) == EOF)
			return true;
		if (!isdigit(s[0]))
			flag_not_digit = true;
	}
	fsetpos(f, &pos);
	return true;
}

//Atribui valor 1
void BoolTable::SetDefault(bool def)
{
	data_table.clear();
	data_table.push_back(def);
}

//Writes output file
void BoolTable::Write(FILE *f)
{
	fprintf(f, "BoolTable ");
	for (int i = 0; i < data_table.size(); i++)
		fprintf(f, "%d ",(int)data_table[i]);
	fprintf(f, "\n");
}

BoolTable::~BoolTable()
{
	data_table.clear();
}

//Returns bool data value at 'number' index (zero-based). If number increases the size of data_table, returns the last available value
bool BoolTable::GetAt(int number)
{
	if (number >= (int) data_table.size())//larger than last position - return last available
		return data_table[data_table.size() - 1];
	else
	{
		if (number >= 0)
			return data_table[number];
		else
			return false;//smaller than the first position - return false
	}
	
}

//Returns linear interpolation factor for current time value
double BoolTable::GetLinearFactorAtCurrentTime()
{
	int sol = db.current_solution_number;
	double l_factor;
	// current solution TRUE
	if (GetAt(sol - 1))
	{
		// previous solution TRUE: keep
		if (GetAt(sol - 2) == true)
			l_factor = 1.0;
		// previous solution FALSE: ramp up
		else
			l_factor = (db.last_converged_time + db.current_time_step - db.solution[sol - 1]->start_time) / (db.solution[sol - 1]->end_time - db.solution[sol - 1]->start_time);
	}
	//current solution FALSE
	else
	{
		// previous solution TRUE: ramp down
		if (GetAt(sol - 2) == true)
			l_factor = 1.0 - (db.last_converged_time + db.current_time_step - db.solution[sol - 1]->start_time) / (db.solution[sol - 1]->end_time - db.solution[sol - 1]->start_time);
		// previous solution FALSE: keep
		else
			l_factor = 0.0;
	}
	return l_factor;
}

//Seta a tabela de dados booleanos
void BoolTable::SetBoolTable(vector<bool> data)
{
	//Limpa atual data_table
	data_table.clear();
	//Re-alocação da data_table
	data_table.resize(data.size());
	for (int i = 0; i < data.size(); i++)
		data_table[i] = data[i];
}

//Operador de Atribuição	
BoolTable &BoolTable::operator = (BoolTable const &bool1)
{
	data_table.clear();
	data_table.resize(bool1.data_table.size());
	for (int i = 0; i < bool1.data_table.size(); i++)
		data_table[i] = bool1.data_table[i];

	return *this;
}

