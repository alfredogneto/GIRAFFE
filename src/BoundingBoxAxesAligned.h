#pragma once
#include "BoundingVolume.h"

class BoundingBoxAxesAligned :
	public BoundingVolume
{
public:
	BoundingBoxAxesAligned();
	~BoundingBoxAxesAligned();
	void SaveConfiguration();
	void Report();

	float x_min, x_max, y_min, y_max, z_min, z_max;
	float prev_x_min, prev_x_max, prev_y_min, prev_y_max, prev_z_min, prev_z_max;
	float x_min_inf, x_max_inf, y_min_inf, y_max_inf, z_min_inf, z_max_inf;
};

