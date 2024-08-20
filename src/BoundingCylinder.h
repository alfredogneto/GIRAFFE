#pragma once
#include "BoundingVolume.h"

class BoundingCylinder :
	public BoundingVolume
{
public:
	BoundingCylinder();
	~BoundingCylinder();
	void SaveConfiguration();
	void Report();

	float radius;
	float ref_radius;

	MatrixFloat* xt;
	MatrixFloat* xb;
	MatrixFloat* prev_xt;
	MatrixFloat* prev_xb;
};

