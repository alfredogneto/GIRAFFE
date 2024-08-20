#pragma once
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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

