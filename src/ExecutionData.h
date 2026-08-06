#pragma once
#include <stdio.h>

class ExecutionData
{
public:
	ExecutionData();
	~ExecutionData();

	bool Read(FILE *f);
	void Write(FILE *f);
	bool print_console;
	bool print_file;
	bool print_contact_report;
};

