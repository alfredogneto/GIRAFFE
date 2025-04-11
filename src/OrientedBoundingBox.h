 #pragma once
#include "BoundingVolume.h"

class MatrixFloat;

//Marina
class OrientedBoundingBox :
	public BoundingVolume
{
public:
	OrientedBoundingBox();
	~OrientedBoundingBox();
	void SaveConfiguration();
	void Report();

	MatrixFloat* x0;
	MatrixFloat* x1;
	MatrixFloat* x2;
	MatrixFloat* x3;
	MatrixFloat* x4;
	MatrixFloat* x5;
	MatrixFloat* x6;
	MatrixFloat* x7;

	MatrixFloat* prev_x0;
	MatrixFloat* prev_x1;
	MatrixFloat* prev_x2;
	MatrixFloat* prev_x3;
	MatrixFloat* prev_x4;
	MatrixFloat* prev_x5;
	MatrixFloat* prev_x6;
	MatrixFloat* prev_x7;
	
	// Orientation of the 3D boxes

	MatrixFloat* orient;
	MatrixFloat* prev_orient;

	// Halfedge lengths

	MatrixFloat* half_dis;
	MatrixFloat* prev_half_dis;

	// Center of the 3D box

	MatrixFloat* center;
	MatrixFloat* prev_center;


};

