#pragma once
#include <stdio.h>
#include <string>
#include <exprtk.hpp>

class MathCode
{
public:
	MathCode();
	MathCode(int e_n_expressions);
	~MathCode();
	double GetValueAt(double time_value, int value_index);

	typedef exprtk::symbol_table<double> symbol_table_t;
	typedef exprtk::expression<double>     expression_t;
	typedef exprtk::parser<double>             parser_t;

	symbol_table_t symbol_table;
	parser_t       parser;
	expression_t*   expressions;
	double t;
	int n_expressions;
	std::string* expressions_read;

	bool Read(FILE* f);					//Leitura do MathCode
	void Write(FILE *f);				//Writes output file
};

