#pragma once
#include <vector>
#include <string>
#include <ctype.h>
using namespace std;

class BoolTable
{
public:
	BoolTable();
	BoolTable(BoolTable &copied);
	BoolTable &operator = (BoolTable const &bool1);		//Operador de Atribuição
	~BoolTable();

	bool Read(FILE *f);			//Reads input file
	void Write(FILE *f);		//Writes output file

	//Useful functions:
	bool GetAt(int number);						//Returns bool data value at 'number' index (zero based). If number increases the size of data_table, returns the last available value
	double GetLinearFactorAtCurrentTime();		//Returns linear interpolation factor for current time value
	void SetDefault(bool def);					//Atribui valor def
	void SetBoolTable(vector<bool> data);		//Seta a tabela de dados booleanos
private:
	vector<bool> data_table;
};

