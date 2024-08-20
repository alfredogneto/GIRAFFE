#pragma once
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "Matrix.h"
#include <vector>
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