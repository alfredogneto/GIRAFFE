#pragma once
#include "BoundingVolume.h"

class BoundingTriangularBox :
	public BoundingVolume
{
public:
	class BoundingTriangularBox();
	~BoundingTriangularBox();
	void SaveConfiguration();
	void Report();

	float thickness;
	float ref_thickness;

	MatrixFloat* x0;
	MatrixFloat* x1;
	MatrixFloat* x2;

	MatrixFloat* prev_x0;
	MatrixFloat* prev_x1;
	MatrixFloat* prev_x2;
};