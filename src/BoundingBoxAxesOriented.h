#pragma once
#include "BoundingVolume.h"

class MatrixFloat;

class BoundingBoxAxesOriented :
	public BoundingVolume
{
public:
	class BoundingBoxAxesOriented();
	~BoundingBoxAxesOriented();
	void SaveConfiguration();
	void Report();

	MatrixFloat* center;
	MatrixFloat* x_local;
	MatrixFloat* y_local;
	MatrixFloat* z_local;
	MatrixFloat* halfwidths;

	MatrixFloat* prev_center;
	MatrixFloat* prev_x_local;
	MatrixFloat* prev_y_local;
	MatrixFloat* prev_z_local;
	MatrixFloat* prev_halfwidths;
};

