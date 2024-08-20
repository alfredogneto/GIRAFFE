//Collision Detection Functions

#include <typeinfo>
#include <float.h> 

#include "BoundingVolume.h"

#include "BoundingSphere.h"
#include "BoundingCylinder.h"
#include "BoundingTriangularBox.h"
#include "BoundingBoxAxesOriented.h"
#include "BoundingBoxAxesAligned.h"

bool CollisionDetection(BoundingVolume* bv_1, BoundingVolume* bv_2);

bool CollisionDetection(BoundingSphere* bv_1, BoundingSphere* bv_2);
bool CollisionDetection(BoundingCylinder* bv_1, BoundingCylinder* bv_2);
bool CollisionDetection(BoundingTriangularBox* bv_1, BoundingTriangularBox* bv_2);
bool CollisionDetection(BoundingBoxAxesOriented* bv_1, BoundingBoxAxesOriented* bv_2);
bool CollisionDetection(BoundingBoxAxesAligned* bv_1, BoundingBoxAxesAligned* bv_2);


bool CollisionDetection(BoundingSphere* bv_1, BoundingCylinder* bv_2);
bool CollisionDetection(BoundingSphere* bv_1, BoundingTriangularBox* bv_2);
bool CollisionDetection(BoundingSphere* bv_1, BoundingBoxAxesOriented* bv_2);
bool CollisionDetection(BoundingSphere* bv_1, BoundingBoxAxesAligned* bv_2);
bool CollisionDetection(BoundingCylinder* bv_1, BoundingTriangularBox* bv_2);
bool CollisionDetection(BoundingCylinder* bv_1, BoundingBoxAxesOriented* bv_2);
bool CollisionDetection(BoundingTriangularBox* bv_1, BoundingBoxAxesOriented* bv_2);

//Auxiliary functions
bool RayTest(MatrixFloat& P0, MatrixFloat& P1, MatrixFloat& P2, MatrixFloat& S, MatrixFloat& T);
float Clamp(float n, float min, float max);