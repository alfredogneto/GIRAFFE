#include "MathCode.h"
#include "IO.h"

#include"Database.h"
//Variáveis globais
extern
Database db;

MathCode::MathCode()
{
	n_expressions = 0;
	expressions = NULL;
	expressions_read = NULL;
	t = 0;
}

MathCode::MathCode(int e_n_expressions)
{
	n_expressions = e_n_expressions;
	expressions = new expression_t[e_n_expressions];
	expressions_read = new string[e_n_expressions];
	symbol_table.add_variable("t", t);
	symbol_table.add_constants();
	for (int i = 0; i < n_expressions; i++)
	{
		expressions[i].register_symbol_table(symbol_table);
		expressions_read[i].clear();
	}
}

//Retorna o valor no instante requerido
double MathCode::GetValueAt(double time_value, int value_index)
{
	t = time_value;
	return expressions[value_index].value();
}

MathCode::~MathCode()
{
	if (expressions != NULL)
		delete[]expressions;
	if (expressions_read != NULL)
	{
		for (int i = 0; i < (int)expressions_read->size(); i++)
		{
			expressions_read[i].clear();
		}
		delete[]expressions_read;
	}
	
}

//Leitura do arquivo de MathCode
bool  MathCode::Read(FILE* f)
{
	char s[1000];
	char expression_code[10000];
	int limit_words = 10000;
	//Trying to find initial comments on the input file
	TryComment(f);
	//Reading codes. "n_expressions" pieces of code are expected to be found in the sequence
	for (int i = 0; i < n_expressions; i++)
	{
		TryComment(f);
		fscanf(f, "%s", s);
		if (!strcmp(s, "Begin"))
		{
			int count = 0;
			while (strcmp(s, "End") != 0 && count < limit_words)
			{
				//First reading
				if (!strcmp(s, "Begin"))
				{
					fscanf(f, "%s", s);
					if (strcmp(s, "End"))
						strcpy(expression_code, s);
				}
				//Remaining readings
				else
				{
					fscanf(f, "%s", s);
					if (strcmp(s, "End"))
						strcat(expression_code, s);
				}
				count++;
			}
		}
		else
		{
			db.myprintf("Error reading MathCode!\n");
			return false;
		}
		expressions_read[i].insert(0, expression_code);
		parser.compile(expression_code, expressions[i]);
	}
	return true;
}

//Writes output file
void MathCode::Write(FILE *f)
{
	for (int i = 0; i < n_expressions; i++)
	{
		fprintf(f, "Begin\t");
		fprintf(f, expressions_read[i].c_str());
		fprintf(f, "\tEnd\n");
	}
	fprintf(f, "\n");
}
