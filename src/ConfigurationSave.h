#pragma once
#include <stdio.h>

class ConfigurationSave
{
public:
	ConfigurationSave();
	~ConfigurationSave();
	bool Check();						
	bool Read(FILE *f);
	void Write(FILE *f);
	void ExportConfiguration(double time);

	int acumulated_config_number;
	int sample;
};