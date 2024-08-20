#pragma once
#include "BoundingVolume.h"

class BoundingSphere :
	public BoundingVolume
{
public:
	BoundingSphere();
	~BoundingSphere();
	void SaveConfiguration();
	void Report();

	float radius;
	float ref_radius;


	MatrixFloat* center;
	MatrixFloat* prev_center;
	
};

