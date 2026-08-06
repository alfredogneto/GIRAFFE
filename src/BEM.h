#pragma once
#include <stdio.h>

class BEM
{
public:
	BEM();
	~BEM();

	bool Read(FILE *f);					//Reading
	void Write(FILE *f);				//Recording
	bool Check();						//Checking inconsistencies
	//variables internal
	int B;		//Number of pas
	double R;		//Radius of the rotor
	double Rhub;	//Radius of the hub
	int CS_rotor;	//Coordinate system for descrever the plane of the rotor (plane xz)
	int node_rotor;	//Node associated the position of the rotor

	double tol_bem;	//tolerance of the BEM
};

