#include "CollisionDetection.h"

#include <float.h>
#include <typeinfo>

#include "BoundingVolume.h"
#include "BoundingSphere.h"
#include "BoundingCylinder.h"
#include "BoundingTriangularBox.h"
#include "BoundingBoxAxesOriented.h"
#include "BoundingBoxAxesAligned.h"
#include "MatrixFloat.h"

//extern
//FILE *fdebug;

bool CollisionDetection(BoundingVolume* bv_1, BoundingVolume* bv_2)
{
	int ref_ID = 0;
	bool ret_value = false;

	if (typeid(*bv_1) == typeid(BoundingSphere) && typeid(*bv_2) == typeid(BoundingSphere))
		ref_ID = 1;
	if (typeid(*bv_1) == typeid(BoundingCylinder) && typeid(*bv_2) == typeid(BoundingCylinder))
		ref_ID = 2;
	if (typeid(*bv_1) == typeid(BoundingTriangularBox) && typeid(*bv_2) == typeid(BoundingTriangularBox))
		ref_ID = 3;
	if (typeid(*bv_1) == typeid(BoundingBoxAxesOriented) && typeid(*bv_2) == typeid(BoundingBoxAxesOriented))
		ref_ID = 4;

	if (typeid(*bv_1) == typeid(BoundingSphere) && typeid(*bv_2) == typeid(BoundingCylinder))
		ref_ID = 5;
	if (typeid(*bv_1) == typeid(BoundingSphere) && typeid(*bv_2) == typeid(BoundingTriangularBox))
		ref_ID = 6;
	if (typeid(*bv_1) == typeid(BoundingSphere) && typeid(*bv_2) == typeid(BoundingBoxAxesOriented))
		ref_ID = 7;
	if (typeid(*bv_1) == typeid(BoundingCylinder) && typeid(*bv_2) == typeid(BoundingTriangularBox))
		ref_ID = 8;
	if (typeid(*bv_1) == typeid(BoundingCylinder) && typeid(*bv_2) == typeid(BoundingBoxAxesOriented))
		ref_ID = 9;
	if (typeid(*bv_1) == typeid(BoundingTriangularBox) && typeid(*bv_2) == typeid(BoundingBoxAxesOriented))
		ref_ID = 10;

	if (typeid(*bv_1) == typeid(BoundingCylinder) && typeid(*bv_2) == typeid(BoundingSphere))
		ref_ID = 11;
	if (typeid(*bv_1) == typeid(BoundingTriangularBox) && typeid(*bv_2) == typeid(BoundingSphere))
		ref_ID = 12;
	if (typeid(*bv_1) == typeid(BoundingBoxAxesOriented) && typeid(*bv_2) == typeid(BoundingSphere))
		ref_ID = 13;
	if (typeid(*bv_1) == typeid(BoundingTriangularBox) && typeid(*bv_2) == typeid(BoundingCylinder))
		ref_ID = 14;
	if (typeid(*bv_1) == typeid(BoundingBoxAxesOriented) && typeid(*bv_2) == typeid(BoundingCylinder))
		ref_ID = 15;
	if (typeid(*bv_1) == typeid(BoundingBoxAxesOriented) && typeid(*bv_2) == typeid(BoundingTriangularBox))
		ref_ID = 16;

	if (typeid(*bv_1) == typeid(BoundingBoxAxesAligned) && typeid(*bv_2) == typeid(BoundingBoxAxesAligned))
		ref_ID = 17;
	if (typeid(*bv_1) == typeid(BoundingSphere) && typeid(*bv_2) == typeid(BoundingBoxAxesAligned))
		ref_ID = 18;
	if (typeid(*bv_1) == typeid(BoundingBoxAxesAligned) && typeid(*bv_2) == typeid(BoundingSphere))
		ref_ID = 19;

	switch (ref_ID) 
	{
		case 0:
			ret_value = false;
			break;
		case 1:
		{
			BoundingSphere* ptr_bv_1 = static_cast<BoundingSphere*>(bv_1);
			BoundingSphere* ptr_bv_2 = static_cast<BoundingSphere*>(bv_2);
			ret_value = CollisionDetection(ptr_bv_1, ptr_bv_2);
			break; 
		}
		case 2:
		{
			BoundingCylinder* ptr_bv_1 = static_cast<BoundingCylinder*>(bv_1);
			BoundingCylinder* ptr_bv_2 = static_cast<BoundingCylinder*>(bv_2);
			ret_value = CollisionDetection(ptr_bv_1, ptr_bv_2);
			break;
		}
		case 3:
		{
			BoundingTriangularBox* ptr_bv_1 = static_cast<BoundingTriangularBox*>(bv_1);
			BoundingTriangularBox* ptr_bv_2 = static_cast<BoundingTriangularBox*>(bv_2);
			ret_value = CollisionDetection(ptr_bv_2, ptr_bv_1);
			break;
		}
		case 4:
		{
			BoundingBoxAxesOriented* ptr_bv_1 = static_cast<BoundingBoxAxesOriented*>(bv_1);
			BoundingBoxAxesOriented* ptr_bv_2 = static_cast<BoundingBoxAxesOriented*>(bv_2);
			ret_value = CollisionDetection(ptr_bv_1, ptr_bv_2);
			break;
		}
		case 5:
		{
			BoundingSphere* ptr_bv_1 = static_cast<BoundingSphere*>(bv_1);
			BoundingCylinder* ptr_bv_2 = static_cast<BoundingCylinder*>(bv_2);
			ret_value = CollisionDetection(ptr_bv_1, ptr_bv_2);
			break;
		}
		case 6:
		{
			BoundingSphere* ptr_bv_1 = static_cast<BoundingSphere*>(bv_1);
			BoundingTriangularBox* ptr_bv_2 = static_cast<BoundingTriangularBox*>(bv_2);
			ret_value = CollisionDetection(ptr_bv_1, ptr_bv_2);
			break;
		}
		case 7:
		{
			BoundingSphere* ptr_bv_1 = static_cast<BoundingSphere*>(bv_1);
			BoundingBoxAxesOriented* ptr_bv_2 = static_cast<BoundingBoxAxesOriented*>(bv_2);
			ret_value = CollisionDetection(ptr_bv_1, ptr_bv_2);
			break;
		}
		case 8:
		{
			BoundingCylinder* ptr_bv_1 = static_cast<BoundingCylinder*>(bv_1);
			BoundingTriangularBox* ptr_bv_2 = static_cast<BoundingTriangularBox*>(bv_2);
			ret_value = CollisionDetection(ptr_bv_1, ptr_bv_2);
			break;
		}
		case 9:
		{
			BoundingCylinder* ptr_bv_1 = static_cast<BoundingCylinder*>(bv_1);
			BoundingBoxAxesOriented* ptr_bv_2 = static_cast<BoundingBoxAxesOriented*>(bv_2);
			ret_value = CollisionDetection(ptr_bv_1, ptr_bv_2);
			break;
		}
		case 10:
		{
			BoundingTriangularBox* ptr_bv_1 = static_cast<BoundingTriangularBox*>(bv_1);
			BoundingBoxAxesOriented* ptr_bv_2 = static_cast<BoundingBoxAxesOriented*>(bv_2);
			ret_value = CollisionDetection(ptr_bv_1, ptr_bv_2);
			break;
		}
		case 11:
		{
			BoundingCylinder* ptr_bv_1 = static_cast<BoundingCylinder*>(bv_1);
			BoundingSphere* ptr_bv_2 = static_cast<BoundingSphere*>(bv_2);
			ret_value = CollisionDetection(ptr_bv_2, ptr_bv_1);
			break;
		}
		case 12:
		{
			BoundingTriangularBox* ptr_bv_1 = static_cast<BoundingTriangularBox*>(bv_1);
			BoundingSphere* ptr_bv_2 = static_cast<BoundingSphere*>(bv_2);
			ret_value = CollisionDetection(ptr_bv_2, ptr_bv_1);
			break;
		}
		case 13:
		{
			BoundingBoxAxesOriented* ptr_bv_1 = static_cast<BoundingBoxAxesOriented*>(bv_1);
			BoundingSphere* ptr_bv_2 = static_cast<BoundingSphere*>(bv_2);
			ret_value = CollisionDetection(ptr_bv_2, ptr_bv_1);
			break;
		}
		case 14:
		{
			BoundingTriangularBox* ptr_bv_1 = static_cast<BoundingTriangularBox*>(bv_1);
			BoundingCylinder* ptr_bv_2 = static_cast<BoundingCylinder*>(bv_2);
			ret_value = CollisionDetection(ptr_bv_2, ptr_bv_1);
			break;
		}
		case 15:
		{
			BoundingBoxAxesOriented* ptr_bv_1 = static_cast<BoundingBoxAxesOriented*>(bv_1);
			BoundingCylinder* ptr_bv_2 = static_cast<BoundingCylinder*>(bv_2);
			ret_value = CollisionDetection(ptr_bv_2, ptr_bv_1);
			break;
		}
		case 16:
		{
			BoundingBoxAxesOriented* ptr_bv_1 = static_cast<BoundingBoxAxesOriented*>(bv_1);
			BoundingTriangularBox* ptr_bv_2 = static_cast<BoundingTriangularBox*>(bv_2);
			ret_value = CollisionDetection(ptr_bv_2, ptr_bv_1);
			break;
		}
		case 17:
		{
			BoundingBoxAxesAligned* ptr_bv_1 = static_cast<BoundingBoxAxesAligned*>(bv_1);
			BoundingBoxAxesAligned* ptr_bv_2 = static_cast<BoundingBoxAxesAligned*>(bv_2);
			ret_value = CollisionDetection(ptr_bv_1, ptr_bv_2);
			break;
		}
		case 18:
		{
			BoundingSphere* ptr_bv_1 = static_cast<BoundingSphere*>(bv_1);
			BoundingBoxAxesAligned* ptr_bv_2 = static_cast<BoundingBoxAxesAligned*>(bv_2);
			ret_value = CollisionDetection(ptr_bv_1, ptr_bv_2);
			break;
		}
		case 19:
		{
			BoundingBoxAxesAligned* ptr_bv_1 = static_cast<BoundingBoxAxesAligned*>(bv_1);
			BoundingSphere* ptr_bv_2 = static_cast<BoundingSphere*>(bv_2);
			ret_value = CollisionDetection(ptr_bv_2, ptr_bv_1);
			break;
		}
		
	}
	return ret_value;
}

bool CollisionDetection(BoundingSphere* bv_1, BoundingSphere* bv_2)
{
	//Distance between centers to square
	MatrixFloat d1 = *bv_1->center - *bv_2->center;
	float dist_cc_square =	dot(d1, d1);
	
	//Sum of radii to square
	float rr_square = (bv_1->radius + bv_2->radius) * (bv_1->radius + bv_2->radius);
	
	//Comparison and return values
	//False: not collision
	//True: collision
	if (dist_cc_square > rr_square)
	{
		//fprintf(fdebug, "Sphere-Sphere\t%.6e\t0\n", sqrt(rr_square));
		
		//Additional test - inversion of dot product considering current and previous center positions
		MatrixFloat d2 = *bv_1->prev_center - *bv_2->prev_center;
		if (dot(d1, d2) >= 0)
			return false;
		else
			return true;
	}
	else
	{
		//fprintf(fdebug, "Sphere-Sphere\t%.6e\t1\n", sqrt(rr_square));
		return true;
	}
		
}

// Clamp n to lie within the range [min, max]
float Clamp(float n, float min, float max) 
{
	if (n < min) return min;
	if (n > max) return max;
	return n;
}

bool CollisionDetection(BoundingCylinder* bv_1, BoundingCylinder* bv_2)
{
	//Minimum distance problem between two segments connecting the extreme points of the cylinders
	float EPSILON = FLT_EPSILON;

	bool ret_value;

	float var1 = 0.5f *(bv_1->bv_factor + bv_1->inc_len_factor);
	float var2 = 0.5f *(bv_2->bv_factor + bv_2->inc_len_factor);

	MatrixFloat d1 = *bv_1->xt - *bv_1->xb;
	MatrixFloat d2 = *bv_2->xt - *bv_2->xb;

	MatrixFloat r = *bv_1->xb - *bv_2->xb;

	float a = dot(d1, d1);
	float b = dot(d1, d2);
	float e = dot(d2, d2);
	float f = dot(d2, r);

	MatrixFloat var3;

	//Convective coordinates
	float s, t;
	//Distance to square
	float dist_square;
	float radsum = bv_1->radius + bv_2->radius;
	float rr_square = (radsum) * (radsum);


	//ERICSON page 149
	// Check if either or both segments degenerate into points
	if (a <= EPSILON && e <= EPSILON)
	{
		// Both segments degenerate into points
		s = 0.0f;
		t = 0.0f;
		var3 = *bv_1->xb - *bv_2->xb;
		dist_square = dot(var3, var3);
	}
	else
	{
		if (a <= EPSILON) 
		{
			// First segment degenerates into a point
			s = 0.0f;
			t = f / e; // s = 0 => t = (b*s + f) / e = f / e
			t = Clamp(t, 0.0f-var2, 1.0f+var2);
		}
		else 
		{
			float c = dot(d1, r);
			if (e <= EPSILON) 
			{
				// Second segment degenerates into a point
				t = 0.0f;
				s = -c / a;	// t = 0 => s = (b*t - c) / a = -c / a
				s = Clamp(s, 0.0f - var1, 1.0f + var1);
			}
			else 
			{
				// The general nondegenerate case starts here
				float b = dot(d1, d2);
				float denom = a * e - b * b; // Always nonnegative
				// If segments not parallel, compute closest point on L1 to L2 and
				// clamp to segment S1. Else pick arbitrary s (here 0)
				if (denom != 0.0f) 
				{
					s = Clamp((b*f - c * e) / denom, 0.0f - var1, 1.0f + var1);
				}
				else s = 0.0f;
				// Compute point on L2 closest to S1(s) using
				// t = Dot((P1 + D1*s) - P2,D2) / Dot(D2,D2) = (b*s + f) / e
				t = (b*s + f) / e;
				// If t in [0,1] done. Else clamp t, recompute s for the new value
				// of t using s = Dot((P2 + D2*t) - P1,D1) / Dot(D1,D1)= (t*b - c) / a
				// and clamp s to [0, 1]
				if (t < (0.0f - var2)) {
					t = 0.0f - var2;
					s = Clamp(-c / a, 0.0f - var1, 1.0f + var1);
				}
				else if (t > (1.0f+var2)) {
					t = 1.0f + var2;
					s = Clamp((b - c) / a, 0.0f - var1, 1.0f + var1);
				}
			}
		}
		MatrixFloat c1 = *bv_1->xb + d1 * s;
		MatrixFloat c2 = *bv_2->xb + d2 * t;
		var3 = c1 - c2;
		dist_square = dot(var3, var3);
	}
	if (dist_square > rr_square)
		ret_value = false;
	else
		ret_value = true;

	//fprintf(fdebug, "Cylinder-Cylinder\t%.6e\t%d\n", sqrt(dist_square), (int)ret_value);
	
	if (ret_value == false)
	{
		//Additional test - inversion of dot product considering current and previous center positions
		MatrixFloat c1c2 = 0.5f*(*bv_1->xt + *bv_1->xb) - 0.5f*(*bv_2->xt + *bv_2->xb);
		MatrixFloat prev_c1prev_c2 = 0.5f*(*bv_1->prev_xt + *bv_1->prev_xb) - 0.5f*(*bv_2->prev_xt + *bv_2->prev_xb);
		
		if (dot(c1c2, prev_c1prev_c2) >= 0)
			return false;
		else
			return true;
	}
	else
		return true;
}

bool CollisionDetection(BoundingTriangularBox* bv_1, BoundingTriangularBox* bv_2)
{
	return false;
}

bool CollisionDetection(BoundingBoxAxesOriented* bv_1, BoundingBoxAxesOriented* bv_2)
{
	//Ericson page 103

	return false;
}

bool CollisionDetection(BoundingSphere* bv_1, BoundingCylinder* bv_2)
{
	//Minimum distance problem between a point and a segment
	bool ret_value;
	MatrixFloat d = *bv_2->xt - *bv_2->xb;
	MatrixFloat xb = *bv_2->xb;
	MatrixFloat r = *bv_1->center - xb;

	float lambda = Clamp((1.0f / dot(d, d)) * dot(r, d), 0.0f - 0.5f * (bv_2->inc_len_factor + bv_2->bv_factor), 1.0f + 0.5f * (bv_2->inc_len_factor + bv_2->bv_factor));
	
	MatrixFloat xbarra = xb + lambda * d;
	//Distance between points
	MatrixFloat dd = *bv_1->center - xbarra;
	float dist_square = dot(dd, dd);

	//Sum of radii to square
	float var2 = bv_1->radius + bv_2->radius;
	float rr_square = (var2) * (var2);

	if (dist_square > rr_square)
		ret_value = false;
	else
		ret_value = true;

	//fprintf(fdebug, "Sphere-Cylinder\t%.6e\t%d\n", sqrt(dist_square), (int)ret_value);

	if (ret_value == false)
	{
		//Additional test - inversion of dot product considering current and previous center positions
		MatrixFloat c2 = 0.5f*(*bv_2->xt + *bv_2->xb);
		MatrixFloat prev_c2 = 0.5f*(*bv_2->prev_xt + *bv_2->prev_xb);

		if (dot(*bv_1->center - c2, *bv_1->prev_center - prev_c2) >= 0.0f)
			return false;
		else
			return true;
	}
	else
		return true;
}

bool CollisionDetection(BoundingSphere* bv_1, BoundingTriangularBox* bv_2)
{
	bool ret_value;

	MatrixFloat centroid = (0.333333f)*(*bv_2->x0 + *bv_2->x1 + *bv_2->x2);

	MatrixFloat d0 = (*bv_2->x0 - centroid);
	MatrixFloat d1 = (*bv_2->x1 - centroid);
	MatrixFloat d2 = (*bv_2->x2 - centroid);

	//Inflated corners of the triangle
	MatrixFloat x0 = *bv_2->x0 + (bv_2->inc_len_factor + bv_2->bv_factor)*d0;
	MatrixFloat x1 = *bv_2->x1 + (bv_2->inc_len_factor + bv_2->bv_factor)*d1;
	MatrixFloat x2 = *bv_2->x2 + (bv_2->inc_len_factor + bv_2->bv_factor)*d2;

	MatrixFloat p1 = x1 - x0;
	MatrixFloat p2 = x2 - x0;

	MatrixFloat var1 = *bv_1->center - x0;

	float a = dot(p1, p1);
	float b = dot(p2, p2);
	float c = dot(p1, p2);
	float d = a * b - c * c;
	float e = dot(var1, p1);
	float f = dot(var1, p2);

	float v1d = (1.0f / d);
	float lambda_1 = v1d *(e*b - c * f);
	float lambda_2 = v1d *(a*f - c * e);

	if (lambda_1 >= 0.0f && lambda_2 >= 0.0f && (lambda_1 + lambda_2) <= 1.0f)
	{
		MatrixFloat xbarra = x0 + lambda_1 * p1 + lambda_2 * p2;

		//Distance between points
		MatrixFloat dd = *bv_1->center - xbarra;
		float dist_square = dot(dd, dd);

		//Sum of radii to square
		float var2 = bv_1->radius + 1.0f * bv_2->thickness;
		float rr_square = (var2) * (var2);

		if (dist_square > rr_square)
			ret_value = false;
		else
			ret_value = true;
		//fprintf(fdebug, "Sphere-TriangularBox\t%.6e\t%d\n", sqrt(dist_square), (int)ret_value);
	}
	else
		ret_value = false;
	//fprintf(fdebug, "Sphere-TriangularBox\t%.6e\t%d\n", 0.0, (int)ret_value);
	
	if (ret_value == false)
	{
		//Additional test - inversion of dot product considering current and previous center positions
		
		MatrixFloat prev_centroid = (0.333333f)*(*bv_2->prev_x0 + *bv_2->prev_x1 + *bv_2->prev_x2);

		if (dot(*bv_1->center - centroid, *bv_1->prev_center - prev_centroid) >= 0)
			return false;
		else
			return true;
	}
	else
		return true;
}

bool CollisionDetection(BoundingSphere* bv_1, BoundingBoxAxesOriented* bv_2)
{
	return false;
}

bool CollisionDetection(BoundingSphere* bv_1, BoundingBoxAxesAligned* bv_2)
{
	return true;
}

bool CollisionDetection(BoundingCylinder* bv_1, BoundingTriangularBox* bv_2)
{
	MatrixFloat x1x0 = *bv_2->x1 - *bv_2->x0;
	MatrixFloat x2x1 = *bv_2->x2 - *bv_2->x1;
	MatrixFloat x0x2 = *bv_2->x0 - *bv_2->x2;
	MatrixFloat t0 = 1.0f / norm(x1x0)*(x1x0);
	MatrixFloat t1 = 1.0f / norm(x2x1)*(x2x1);
	MatrixFloat t2 = 1.0f / norm(x0x2)*(x0x2);
	MatrixFloat t0t1 = cross(t0, t1);
	MatrixFloat nt = 1.0f / norm(t0t1)*(t0t1);
	MatrixFloat nb = -1.0f * nt;
	MatrixFloat a0 = cross(t0, nt);
	MatrixFloat a1 = cross(t1, nt);
	MatrixFloat a2 = cross(t2, nt);
	MatrixFloat Pt = 0.333333f * (*bv_2->x0 + *bv_2->x1 + *bv_2->x2) + 0.5f * bv_2->thickness*nt;
	MatrixFloat Pb = Pt - bv_2->thickness*nt;
	MatrixFloat P0 = 0.5f * (*bv_2->x0 + *bv_2->x1);
	MatrixFloat P1 = 0.5f * (*bv_2->x1 + *bv_2->x2);
	MatrixFloat P2 = 0.5f * (*bv_2->x0 + *bv_2->x2);

	float m1, m2;
	
	//Test 1:
	m1 = dot(*bv_1->xb - Pt, nt);
	m2 = dot(*bv_1->xt - Pt, nt);
	if (m1 > bv_1->radius && m2 > bv_1->radius) return false;
	//Test 2:
	m1 = dot(*bv_1->xb - Pb, nb);
	m2 = dot(*bv_1->xt - Pb, nb);
	if (m1 > bv_1->radius && m2 > bv_1->radius) return false;
	//Test 3:
	m1 = dot(*bv_1->xb - P0, a0);
	m2 = dot(*bv_1->xt - P0, a0);
	if (m1 > bv_1->radius && m2 > bv_1->radius) return false;
	//Test 4:
	m1 = dot(*bv_1->xb - P1, a1);
	m2 = dot(*bv_1->xt - P1, a1);
	if (m1 > bv_1->radius && m2 > bv_1->radius) return false;
	//Test 5:
	m1 = dot(*bv_1->xb - P2, a2);
	m2 = dot(*bv_1->xt - P2, a2);
	if (m1 > bv_1->radius && m2 > bv_1->radius) return false;

	return true;
}

bool CollisionDetection(BoundingCylinder* bv_1, BoundingBoxAxesOriented* bv_2)
{
	float m1, m2;
	MatrixFloat Pi, ni;

	//Test 1:
	Pi = *bv_2->center - (*bv_2->halfwidths)(0, 0)*(*bv_2->x_local);
	ni = -1.0f * (*bv_2->x_local);
	m1 = dot(*bv_1->xb - Pi, ni);
	m2 = dot(*bv_1->xt - Pi, ni);
	if (m1 > bv_1->radius && m2 > bv_1->radius) return false;
	//Test 2:
	Pi = *bv_2->center + (*bv_2->halfwidths)(0, 0)*(*bv_2->x_local);
	ni = +1.0f * (*bv_2->x_local);
	m1 = dot(*bv_1->xb - Pi, ni);
	m2 = dot(*bv_1->xt - Pi, ni);
	if (m1 > bv_1->radius && m2 > bv_1->radius) return false;

	//Test 3:
	Pi = *bv_2->center - (*bv_2->halfwidths)(1, 0)*(*bv_2->y_local);
	ni = -1.0f * (*bv_2->y_local);
	m1 = dot(*bv_1->xb - Pi, ni);
	m2 = dot(*bv_1->xt - Pi, ni);
	if (m1 > bv_1->radius && m2 > bv_1->radius) return false;
	//Test 4:
	Pi = *bv_2->center + (*bv_2->halfwidths)(1, 0)*(*bv_2->y_local);
	ni = +1.0f * (*bv_2->y_local);
	m1 = dot(*bv_1->xb - Pi, ni);
	m2 = dot(*bv_1->xt - Pi, ni);
	if (m1 > bv_1->radius && m2 > bv_1->radius) return false;

	//Test 5:
	Pi = *bv_2->center - (*bv_2->halfwidths)(2, 0)*(*bv_2->z_local);
	ni = -1.0f * (*bv_2->z_local);
	m1 = dot(*bv_1->xb - Pi, ni);
	m2 = dot(*bv_1->xt - Pi, ni);
	if (m1 > bv_1->radius && m2 > bv_1->radius) return false;
	//Test 6:
	Pi = *bv_2->center + (*bv_2->halfwidths)(2, 0)*(*bv_2->z_local);
	ni = +1.0f * (*bv_2->z_local);
	m1 = dot(*bv_1->xb - Pi, ni);
	m2 = dot(*bv_1->xt - Pi, ni);
	if (m1 > bv_1->radius && m2 > bv_1->radius) return false;
	
	return true; // boxes presents overlap with cylinder
}

bool CollisionDetection(BoundingTriangularBox* bv_1, BoundingBoxAxesOriented* bv_2)
{
	return false;
}

//Detects if there is intersection between a triangular face and a straight line
bool RayTest(MatrixFloat& P0, MatrixFloat& P1, MatrixFloat& P2, MatrixFloat& S, MatrixFloat& T)
{
	//Returns:	true (intersection found)
	//			false (no intersection found)

	//Zeros on baricentric coordinates
	float mu0, mu1, mu2;
	
	//Algorithm to detect intersection between a ray and a triangular region 
	//(Bergen:  Book - Colision detection in interative 3D environments, page 85)
	MatrixFloat d1 = (P1 - P0);
	MatrixFloat d2 = (P2 - P0);
	MatrixFloat r = (T - S);
	MatrixFloat d1_d2 = cross(d1, d2);

	float den = dot(r, d1_d2);

	if (abs(den) <= FLT_EPSILON)
		return false;	//we assume here that paralelism means "no intersection"
	else
	{
		float coef = -1.0f / den;
		MatrixFloat b = (S - P0);
		float lambda = coef * dot(b, d1_d2);
		if (lambda >= 0.0f && lambda <= 1.0f)
		{
			MatrixFloat b_r = cross(b, r);
			mu1 = coef * dot(d2, b_r);
			mu2 = -1.0f * coef * dot(d1, b_r);
			mu0 = 1.0f - mu1 - mu2;
			if (mu1 >= 0.0f && mu2 >= 0.0f && (mu1 + mu2) <= 1.0f)
				return true;
			else
				return false;
		}
		else
			return false;
	}
}

bool CollisionDetection(BoundingBoxAxesAligned* bv_1, BoundingBoxAxesAligned* bv_2)
{
	if (bv_1->x_max_inf < bv_2->x_min_inf) return false;
	if (bv_2->x_max_inf < bv_1->x_min_inf) return false;

	if (bv_1->y_max_inf < bv_2->y_min_inf) return false;
	if (bv_2->y_max_inf < bv_1->y_min_inf) return false;

	if (bv_1->z_max_inf < bv_2->z_min_inf) return false;
	if (bv_2->z_max_inf < bv_1->z_min_inf) return false;

	return true; // boxes present overlap
}
