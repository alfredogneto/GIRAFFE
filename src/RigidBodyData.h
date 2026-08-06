#pragma once
#include <stdio.h>

class Matrix;

class RigidBodyData
{
public:
	RigidBodyData();
	~RigidBodyData();

	int number;										//ID
	int CADData_ID;									//ID do CAD
	bool CAD_entered;
	char CAD_file[500];								//Filename from a directly-embedded "CADFile" mesh reference (visualization only - not yet wired to a renderer)
	bool CAD_file_entered;
	double mass;									//mass
	Matrix* J_G;									//Tensor of inertia
	Matrix* G;										//Position of the baricentro

	bool Read(FILE *f);
	void Write(FILE *f);

	void PreCalc();
	void WriteVTK_XMLRender(FILE *f, int pole_node, int cs);//Plots rigid body - formato XML VTK - receives the number of the node that and the pole and the coordinate system of reference (that and updated according to rotations sofridas through the node)
};

