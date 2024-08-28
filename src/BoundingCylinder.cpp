#include "BoundingCylinder.h"

#include "MatrixFloat.h"
#include "Database.h"
//Variáveis globais
extern
Database db;

//extern
//FILE *fdebug;

BoundingCylinder::BoundingCylinder()
{
	radius = 0.0;
	ref_radius = 0.0;

	bv_factor = 0.0f;
	inc_len_factor = 0.0f;

	xt = new MatrixFloat(3);
	xb = new MatrixFloat(3);

	prev_xt = new MatrixFloat(3);
	prev_xb = new MatrixFloat(3);

	first_set = true;

	factor_kin2 = 1.2f*1.2f;

	x_center[0] = 0.0f;
	x_center[1] = 0.0f;
	x_center[2] = 0.0f;
	size = 0.0f;
}

BoundingCylinder::~BoundingCylinder()
{
	delete xt;
	delete xb;

	delete prev_xt;
	delete prev_xb;
}
void BoundingCylinder::SaveConfiguration()
{
	//If first set, prior to update radius, updates the center
	if (first_set)
	{
		*prev_xt = *xt;
		*prev_xb = *xb;
		x_center[0] = 0.5f * ((*xt)(0, 0) + (*xb)(0, 0));
		x_center[1] = 0.5f * ((*xt)(1, 0) + (*xb)(1, 0));
		x_center[2] = 0.5f * ((*xt)(2, 0) + (*xb)(2, 0));
		ref_radius = radius;
		first_set = false;
	}
	////Updating radius - kinematics
	//MatrixFloat avg = 0.5f * (*prev_xt + *prev_xb - *xt - *xb);
	//float d2 = factor_kin2 * dot(avg, avg);
	//if (d2 > ref_radius*ref_radius)
	//	radius = sqrt(d2);
	//else
	//	radius = ref_radius;

	////Updating bounding factor - kinematics
	//bv_factor = radius / ref_radius - 1.0f;

	float prevc[3];
	prevc[0] = x_center[0];
	prevc[1] = x_center[1];
	prevc[2] = x_center[2];

	x_center[0] = 0.5f * ((*xt)(0, 0) + (*xb)(0, 0));
	x_center[1] = 0.5f * ((*xt)(1, 0) + (*xb)(1, 0));
	x_center[2] = 0.5f * ((*xt)(2, 0) + (*xb)(2, 0));

	last_center_change[0] = x_center[0] - prevc[0];
	last_center_change[1] = x_center[1] - prevc[1];
	last_center_change[2] = x_center[2] - prevc[2];

	*prev_xt = *xt;
	*prev_xb = *xb;
}
void BoundingCylinder::Report()
{
	/*float var1 = 0.5f *(bv_factor + inc_len_factor);
	MatrixFloat d1 = *xt - *xb;
	MatrixFloat m1 = 0.5f*(*xt + *xb);
	MatrixFloat c1 = *xb + d1 * (-var1);
	float size = sqrt(dot(c1 - m1, c1 - m1));
	size = sqrt(size*size + radius * radius);*/
	//Printing data to debug file
	//fprintf(fdebug, "BoundingCylinder\t%.6e\n", radius);

	//db.myprintf("\nBounding Cylinder\tradius\t%.3f\tref_radius\t%.3f\tbv_factor\t%.3f\tinc_len_factor\t%.3f\n",
	//	radius, ref_radius, bv_factor, inc_len_factor);
	db.myprintf("BoundingCylinder\n");
	db.myprintf("Size\t%.6e\n", size);
	db.myprintf("Center\t%.6e\t%.6e\t%.6e\n", x_center[0], x_center[1], x_center[2]);
}
