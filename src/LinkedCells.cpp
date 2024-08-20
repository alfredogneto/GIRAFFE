#include "LinkedCells.h"

#include"Database.h"

//Variáveis globais
extern
Database db;

LinkedCells::LinkedCells()
{
	nx = 1;
	ny = 1;
	nz = 1;

	pointers = NULL;
	n = NULL;

	max_x = 0.0f;
	min_x = 0.0f;
	max_y = 0.0f;
	min_y = 0.0f;
	max_z = 0.0f;
	min_z = 0.0f;

	deltax = 0.0f;
	deltay = 0.0f;
	deltaz = 0.0f;

	Lx = 0.0f;
	Ly = 0.0f;
	Lz = 0.0f;
}

LinkedCells::~LinkedCells()
{
	//Clear memory
	if (pointers != NULL)
	{
		for (int i = 0; i < nx; i++)
		{
			for (int j = 0; j < ny; j++)
			{
				delete[] pointers[i][j];
			}
			delete[] pointers[i];
		}
		delete[] pointers;
	}

	if (n != NULL)
	{
		for (int i = 0; i < nx; i++)
		{
			for (int j = 0; j < ny; j++)
				delete[] n[i][j];
			delete[] n[i];
		}
		delete[] n;
	}
}

//Allocation of memory and initial data
void LinkedCells::PreCalc(float e_max_size)
{
	max_x = db.gcs->max_x;
	min_x = db.gcs->min_x;
	max_y = db.gcs->max_y;
	min_y = db.gcs->min_y;
	max_z = db.gcs->max_z;
	min_z = db.gcs->min_z;

	Lx = max_x - min_x;
	Ly = max_y - min_y;
	Lz = max_z - min_z;

	if (e_max_size != 0.0f)
	{
		nx = (int)(Lx / (2.0f*e_max_size));
		ny = (int)(Ly / (2.0f*e_max_size));
		nz = (int)(Lz / (2.0f*e_max_size));
	}
	else
	{
		nx = 1;
		ny = 1;
		nz = 1;
	}

	/*if (nx == 0)
		nx = 1;
	if (ny == 0)
		ny = 1;
	if (nz == 0)
		nz = 1;*/

	deltax = Lx / nx;
	deltay = Ly / ny;
	deltaz = Lz / nz;

	//Allocating pointer vectors
	pointers = new BoundingVolume***[nx];
	for (int i = 0; i < nx; i++)
	{
		pointers[i] = new BoundingVolume**[ny];
		for (int j = 0; j < ny; j++)
		{
			pointers[i][j] = new BoundingVolume*[nz];
			for (int k = 0; k < nz; k++)
			{
				pointers[i][j][k] = NULL;
			}
		}
	}

	n = new int**[nx];
	for (int i = 0; i < nx; i++)
	{
		n[i] = new int*[ny];
		for (int j = 0; j < ny; j++)
		{
			n[i][j] = new int[nz];
			for (int k = 0; k < nz; k++)
			{
				n[i][j][k] = 0;
			}
		}
	}
}

//Clears information
void LinkedCells::EmptyCells()
{
	for (int i = 0; i < nx; i++)
	{
		for (int j = 0; j < ny; j++)
		{
			for (int k = 0; k < nz; k++)
			{
				n[i][j][k] = 0;
			}
		}
	}
	for (int i = 0; i < nx; i++)
	{
		for (int j = 0; j < ny; j++)
		{
			for (int k = 0; k < nz; k++)
			{
				pointers[i][j][k] = NULL;
			}
		}
	}
}

//Adds a pointer to the root of the list of the cell i,j,k
bool LinkedCells::InsertBoundingVolume(BoundingVolume* ext_ptr)
{
	//Evaluating cell position:
	int i = (int)((ext_ptr->x_center[0] - min_x) / deltax);
	int j = (int)((ext_ptr->x_center[1] - min_y) / deltay); 
	int k = (int)((ext_ptr->x_center[2] - min_z) / deltaz);
	if (i >= 0 && i < nx && j >= 0 && j < ny && k >= 0 && k < nz)
	{
		ext_ptr->i_loc = i;
		ext_ptr->j_loc = j;
		ext_ptr->k_loc = k;
		if (pointers[i][j][k] != NULL)
		{
			//Sets as previous a null pointer
			ext_ptr->previous = NULL;
			//Sets as the next of the new pointer, the old root list
			ext_ptr->next = pointers[i][j][k];
			//Sets previous of pointers[i][j][k] the new pointer
			pointers[i][j][k]->previous = ext_ptr;
			//Sets the new root list
			pointers[i][j][k] = ext_ptr;
		}
		else
		{ 
			ext_ptr->previous = NULL;
			ext_ptr->next = NULL;
			//Sets the new root list
			pointers[i][j][k] = ext_ptr;
		}
		//Increments the number of BoudingVolumes in the list
		n[i][j][k]++;
		return true;
	}
	else
		return false;	//out of the cell region
}

//Returns the list root of the cell i,j,k
BoundingVolume* LinkedCells::GetCellListRoot(int i, int j, int k)
{
	if (i >= 0 && i < nx && j >= 0 && j < ny && k >= 0 && k < nz)
		return pointers[i][j][k];
	else
		return NULL;	//out of range
}

//Returns the number of objects in the cell i,j,k
int LinkedCells::GetNObjects(int i, int j, int k)
{
	if (i >= 0 && i < nx && j >= 0 && j < ny && k >= 0 && k < nz)
		return n[i][j][k];
	else
		return 0;	//out of range
}

void LinkedCells::ReportCells()
{
	db.myprintf("LinkedCellsReport");
	db.myprintf("nx = %d, ny = %d, nz = %d\n",nx, ny, nz);
	for (int i = 0; i < nx; i++)
	{
		for (int j = 0; j < ny; j++)
		{
			for (int k = 0; k < nz; k++)
			{
				db.myprintf("\ti = %d, j = %d, k = %d\n", i, j, k);
				db.myprintf("\t\t\t\t\tEntities: %d\n",n[i][j][k]);
			}
		}
	}
}