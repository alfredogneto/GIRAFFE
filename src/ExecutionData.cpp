#include "ExecutionData.h"
#include <string>

ExecutionData::ExecutionData()
{
	print_console = true;
	print_file = true;
	print_contact_report = false;
}

ExecutionData::~ExecutionData()
{
}

bool ExecutionData::Read(FILE *f)
{
	char s[1000];
	fscanf(f, "%s", s);
	if (!strcmp(s, "PrintConsoleReport"))
	{
		fscanf(f, "%s", s);
		print_console = (bool)atoi(s);
	}
	else
		return false;
	fscanf(f, "%s", s);
	if (!strcmp(s, "PrintFileReport"))
	{
		fscanf(f, "%s", s);
		print_file = (bool)atoi(s);
	}
	else
		return false;
	fscanf(f, "%s", s);
	if (!strcmp(s, "PrintContactReport"))
	{
		fscanf(f, "%s", s);
		print_contact_report = (bool)atoi(s);
	}
	else
		return false;

	return true;
}
void ExecutionData::Write(FILE *f)
{
	fprintf(f, "ExecutionData\n");
	fprintf(f, "PrintConsoleReport\t%d\n", (int)print_console);
	fprintf(f, "PrintFileReport\t%d\n", (int)print_file);
	fprintf(f, "PrintContactReport\t%d\n\n", (int)print_contact_report);
}
