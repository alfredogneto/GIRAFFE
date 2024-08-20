#include "BoundingBoxAxesAligned.h"

#include "Database.h"

//Variáveis globais
extern
Database db;

BoundingBoxAxesAligned::BoundingBoxAxesAligned()
{
	bv_factor = 0.0f;
	inc_len_factor = 0.0f;

	first_set = true;

	factor_kin2 = 1.2f*1.2f;

	x_center[0] = 0.0f;
	x_center[1] = 0.0f;
	x_center[2] = 0.0f;
	size = 0.0f;

	x_min = 0.0f;
	x_max = 0.0f;
	y_min = 0.0f;
	y_max = 0.0f;
	z_min = 0.0f;
	z_max = 0.0f;

	x_min_inf = 0.0f;
	x_max_inf = 0.0f;
	y_min_inf = 0.0f;
	y_max_inf = 0.0f;
	z_min_inf = 0.0f;
	z_max_inf = 0.0f;

	prev_x_min = 0.0f;
	prev_x_max = 0.0f;
	prev_y_min = 0.0f;
	prev_y_max = 0.0f;
	prev_z_min = 0.0f;
	prev_z_max = 0.0f;

}

BoundingBoxAxesAligned::~BoundingBoxAxesAligned()
{
}

void BoundingBoxAxesAligned::SaveConfiguration()
{
	//If first set
	if (first_set)
	{
		prev_x_min = x_min;
		prev_x_max = x_max;
		prev_y_min = y_min;
		prev_y_max = y_max;
		prev_z_min = z_min;
		prev_z_max = z_max;
		
		x_center[0] = 0.5f * (x_min + x_max);
		x_center[1] = 0.5f * (y_min + y_max);
		x_center[2] = 0.5f * (z_min + z_max);
		
		first_set = false;
	}

	float prevc[3];
	prevc[0] = x_center[0];
	prevc[1] = x_center[1];
	prevc[2] = x_center[2];

	x_center[0] = 0.5f * (x_min + x_max);
	x_center[1] = 0.5f * (y_min + y_max);
	x_center[2] = 0.5f * (z_min + z_max);

	x_max_inf = x_center[0] + 0.5f*(x_max - x_min) * inc_len_factor;
	y_max_inf = x_center[1] + 0.5f*(y_max - y_min) * inc_len_factor;
	z_max_inf = x_center[2] + 0.5f*(z_max - z_min) * inc_len_factor;

	x_min_inf = x_center[0] - 0.5f*(x_max - x_min) * inc_len_factor;
	y_min_inf = x_center[1] - 0.5f*(y_max - y_min) * inc_len_factor;
	z_min_inf = x_center[2] - 0.5f*(z_max - z_min) * inc_len_factor;


	last_center_change[0] = x_center[0] - prevc[0];
	last_center_change[1] = x_center[1] - prevc[1];
	last_center_change[2] = x_center[2] - prevc[2];

	prev_x_min = x_min;
	prev_x_max = x_max;
	prev_y_min = y_min;
	prev_y_max = y_max;
	prev_z_min = z_min;
	prev_z_max = z_max;
}
void BoundingBoxAxesAligned::Report()
{
	db.myprintf("BoundingBoxAxesAligned\n");
}
