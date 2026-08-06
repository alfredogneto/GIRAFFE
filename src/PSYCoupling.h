#pragma once
#include <stdio.h>
#include "BoolTable.h"

class Constraint;
class Displacement;

class PSYCoupling
{
public:
	PSYCoupling();
	~PSYCoupling();

	bool Read(FILE *f);		//Reading in the input file of the Giraffe
	void Write(FILE *f);	//Output in the input file of the Giraffe

	void PreCalc();		
	void ReadPSYFile();
	void WritePSYFile();
	void CoupleByFile();
	void CoupleByBinary();
	void Couple();
	void SetConstraints();

	
	bool couple_by_file;						//Variable that indicates coupling by file
	BoolTable PSY_bool;
	int number_psy_displacements;
	Displacement** psy_displacements;			//Vector of displacements
	int number_psy_constraints;
	Constraint** psy_constraints;				//Vector of constraints
};

