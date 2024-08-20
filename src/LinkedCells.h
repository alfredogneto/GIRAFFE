#pragma once
#include "BoundingVolume.h"

class LinkedCells
{
public:
	LinkedCells();
	~LinkedCells();
	

	float Lx, Ly, Lz;				//Dimensions of the region of interest (to create cells)
	int nx, ny, nz;					//Number of divisions in each dimension
	float deltax, deltay, deltaz;	//Steps of divisions
	float max_x;					//max x
	float min_x;					//min x
	float max_y;					//max y
	float min_y;					//min y
	float max_z;					//max z
	float min_z;					//min z

	BoundingVolume**** pointers;	//Pointers to each cell BV's (start of the linked list)
	int*** n;						//Number of BV's in each cell

	//Allocation of memory and initial data
	void PreCalc(float e_max_size);
	void EmptyCells();															//Clears information (pointers and number of cells)
	bool InsertBoundingVolume(BoundingVolume* ext_ptr);							//Adds a pointer to the root of the list of the cell i,j,k
	BoundingVolume* GetCellListRoot(int i, int j, int k);						//Returns the list root of the cell i,j,k
	int GetNObjects(int i, int j, int k);										//Returns the number of objects in the cell i,j,k
	void ReportCells();															//Prints useful information
};

