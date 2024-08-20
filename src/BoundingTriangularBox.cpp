#include "BoundingTriangularBox.h"

#include "Database.h"

//Variáveis globais
extern
Database db;

//extern
//FILE *fdebug;

BoundingTriangularBox::BoundingTriangularBox()
{
	thickness = 0.0f;
	ref_thickness = 0.0f;

	bv_factor = 0.0f;
	inc_len_factor = 0.0f;

	x0 = new MatrixFloat(3);
	x1 = new MatrixFloat(3);
	x2 = new MatrixFloat(3);

	prev_x0 = new MatrixFloat(3);
	prev_x1 = new MatrixFloat(3);
	prev_x2 = new MatrixFloat(3);

	first_set = true;

	factor_kin2 = 1.2f*1.2f;

	x_center[0] = 0.0f;
	x_center[1] = 0.0f;
	x_center[2] = 0.0f;
	size = 0.0f;
}

BoundingTriangularBox::~BoundingTriangularBox()
{
	delete x0;
	delete x1;
	delete x2;

	delete prev_x0; 
	delete prev_x1;
	delete prev_x2;
}
void BoundingTriangularBox::SaveConfiguration()
{
	//If first set, prior to update radius, updates the center
	if (first_set)
	{
		*prev_x0 = *x0;
		*prev_x1 = *x1;
		*prev_x2 = *x2;
		x_center[0] = 0.333333f * ((*x0)(0, 0) + (*x1)(0, 0) + (*x2)(0, 0));
		x_center[1] = 0.333333f * ((*x0)(1, 0) + (*x1)(1, 0) + (*x2)(1, 0));
		x_center[2] = 0.333333f * ((*x0)(2, 0) + (*x1)(2, 0) + (*x2)(2, 0));
		ref_thickness = thickness;
		first_set = false;
	}
	//Updating radius - kinematics
	//MatrixFloat avg = 0.333333f * (*prev_x0 + *prev_x1 + *prev_x2 - *x0 - *x1 - *x2);
	//float d2 = factor_kin2 * dot(avg, avg);
	//if (d2 > ref_thickness*ref_thickness)
	//	thickness = sqrt(d2);
	//else
	//	thickness = ref_thickness;
	////Updating bounding factor - kinematics
	//bv_factor = thickness / ref_thickness - 1.0f;

	float prevc[3];
	prevc[0] = x_center[0];
	prevc[1] = x_center[1];
	prevc[2] = x_center[2];

	x_center[0] = 0.333333f * ((*x0)(0, 0) + (*x1)(0, 0) + (*x2)(0, 0));
	x_center[1] = 0.333333f * ((*x0)(1, 0) + (*x1)(1, 0) + (*x2)(1, 0));
	x_center[2] = 0.333333f * ((*x0)(2, 0) + (*x1)(2, 0) + (*x2)(2, 0));

	last_center_change[0] = x_center[0] - prevc[0];
	last_center_change[1] = x_center[1] - prevc[1];
	last_center_change[2] = x_center[2] - prevc[2];

	*prev_x0 = *x0;
	*prev_x1 = *x1;
	*prev_x2 = *x2;

}
void BoundingTriangularBox::Report()
{
	//MatrixFloat centroid = (0.333333f)*(*x0 + *x1 + *x2);

	//MatrixFloat d0 = (*x0 - centroid);
	//MatrixFloat d1 = (*x1 - centroid);
	//MatrixFloat d2 = (*x2 - centroid);

	////Inflated corners of the triangle
	//MatrixFloat x0a = *x0 + (inc_len_factor + bv_factor)*d0;
	//MatrixFloat x1a = *x1 + (inc_len_factor + bv_factor)*d1;
	//MatrixFloat x2a = *x2 + (inc_len_factor + bv_factor)*d2;

	//float size = sqrt(dot(x0a - centroid, x0a - centroid));
	//if (size < sqrt(dot(x1a - centroid, x1a - centroid)))
	//	size = sqrt(dot(x1a - centroid, x1a - centroid));
	//if (size < sqrt(dot(x2a - centroid, x2a - centroid)))
	//	size = sqrt(dot(x2a - centroid, x2a - centroid));
	//size = sqrt(size*size + thickness * thickness);
	//Printing data to debug file
	//fprintf(fdebug, "BoundingTriangularBox\t%.6e\n", size);
	db.myprintf("BoundingTriangularBox\n");
	db.myprintf("Size\t%.6e\n", size);
	db.myprintf("Center\t%.6e\t%.6e\t%.6e\n", x_center[0], x_center[1], x_center[2]);
}
