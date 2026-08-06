#pragma once
#include <stdio.h>


class ContactInterface
{
public:
	ContactInterface() {}
	virtual ~ContactInterface() {}

	virtual bool Check() = 0;				//Checks inconsistencies for avoid errors of execution
	virtual bool Read(FILE *f) = 0;
	virtual void Write(FILE *f) = 0;
	virtual void PreCalc() = 0;				//Pre-calculation of variables that and done a single time in the start


	int number;			//ID
	int material_1;		//Material 1 ID
	int material_2;		//Material 2 ID
};

