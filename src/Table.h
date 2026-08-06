#pragma once
#include <stdio.h>

class Table
{
public:
	Table();
	Table(int e_times, int e_values);
	~Table();

	
	int n_times;							//Stores the number of instantes saved (rows)
	int n_values;							//Number of values of the table (columns)
	
	void SetTime(int time_index, double time);
	void SetValue(int time_index, int value_index, double value);

	double GetTime(int time_index);
	double GetValue(int time_index, int value_index);
	double GetValueAt(double time_value, int value_index);	//Returns the value in the instant required (interpolation linear)

	void flush();

	bool Read(FILE *f);
	void Write(FILE *f);

protected:
	double** values;						//Stores the values in function of the time
	double* times;							//Stores the instantes associated to the values
};