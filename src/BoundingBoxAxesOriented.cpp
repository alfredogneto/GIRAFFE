#include "BoundingBoxAxesOriented.h"

#include "MatrixFloat.h"
#include "Database.h"

//Variaveis globais
extern
Database db;
//extern
//FILE *fdebug;

BoundingBoxAxesOriented::BoundingBoxAxesOriented()
{
	bv_factor = 0.0f;
	inc_len_factor = 0.0f;

	center = new MatrixFloat(3);
	x_local = new MatrixFloat(3);
	y_local = new MatrixFloat(3);
	z_local = new MatrixFloat(3);
	halfwidths = new MatrixFloat(3);

	prev_center = new MatrixFloat(3);
	prev_x_local = new MatrixFloat(3);
	prev_y_local = new MatrixFloat(3);
	prev_z_local = new MatrixFloat(3);
	prev_halfwidths = new MatrixFloat(3);
	
	first_set = true;

	factor_kin2 = 1.2f*1.2f;

	x_center[0] = 0.0f;
	x_center[1] = 0.0f;
	x_center[2] = 0.0f;
	size = 0.0f;
}

BoundingBoxAxesOriented::~BoundingBoxAxesOriented()
{
	delete center;
	delete x_local;
	delete y_local;
	delete z_local;
	delete halfwidths;

	delete prev_center;
	delete prev_x_local;
	delete prev_y_local;
	delete prev_z_local;
	delete prev_halfwidths;
}
void BoundingBoxAxesOriented::SaveConfiguration()
{
	*prev_center = *center;
	*prev_x_local = *x_local;
	*prev_y_local = *y_local;
	*prev_z_local = *z_local;
	*prev_halfwidths = *halfwidths;
}
void BoundingBoxAxesOriented::Report()
{
	db.myprintf("BoundingBoxAxesOriented\n");
	db.myprintf("center\n");
	center->print();
	db.myprintf("halfwidths\n");
	halfwidths->print();
	db.myprintf("x_local\n");
	x_local->print();
	db.myprintf("y_local\n");
	y_local->print();
	db.myprintf("z_local\n");
	z_local->print();
}
