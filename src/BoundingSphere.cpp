#include "BoundingSphere.h"

#include "Database.h"

//Variáveis globais
extern
Database db;

//extern
//FILE *fdebug;

BoundingSphere::BoundingSphere()
{

	radius = 0.0;
	ref_radius = 0.0;

	bv_factor = 0.0f;
	inc_len_factor = 0.0f;

	center = new MatrixFloat(3);
	prev_center = new MatrixFloat(3);

	first_set = true;

	factor_kin2 = 1.2f*1.2f;

	x_center[0] = 0.0f;
	x_center[1] = 0.0f;
	x_center[2] = 0.0f;
	size = 0.0f;
}

BoundingSphere::~BoundingSphere()
{
	delete center;
	delete prev_center;
}
void BoundingSphere::SaveConfiguration()
{
	//If first set, prior to update radius, updates the center
	if (first_set)
	{
		*prev_center = *center;
		x_center[0] = (*center)(0, 0);
		x_center[1] = (*center)(1, 0);
		x_center[2] = (*center)(2, 0);
		ref_radius = radius;
		first_set = false;
	}
	////Updating radius - kinematics
	//MatrixFloat dif = *center - *prev_center;
	//float d2 = factor_kin2 * dot(dif,dif);
	//if (d2 > ref_radius*ref_radius)
	//	radius = sqrt(d2);
	//else
	//	radius = ref_radius;

	float prevc[3];
	prevc[0] = x_center[0];
	prevc[1] = x_center[1];
	prevc[2] = x_center[2];

	x_center[0] = (*center)(0, 0);
	x_center[1] = (*center)(1, 0);
	x_center[2] = (*center)(2, 0);

	last_center_change[0] = x_center[0] - prevc[0];
	last_center_change[1] = x_center[1] - prevc[1];
	last_center_change[2] = x_center[2] - prevc[2];

	*prev_center = *center;
}
void BoundingSphere::Report()
{
	db.myprintf("BoundingSphere\n");
	db.myprintf("Size\t%.6e\n",size);
	db.myprintf("Center\t%.6e\t%.6e\t%.6e\n", x_center[0], x_center[1], x_center[2]);
}
