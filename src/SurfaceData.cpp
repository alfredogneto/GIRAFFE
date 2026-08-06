#include "SurfaceData.h"

#include "Matrix.h"

SurfaceData::SurfaceData(int n_solutions)
{
	n_sol = n_solutions;
	G_p = new Matrix*[n_solutions];
	for (int i = 0; i < n_solutions; i++)
		G_p[i] = new Matrix(3);

	t1_p = new Matrix*[n_solutions];
	for (int i = 0; i < n_solutions; i++)
		t1_p[i] = new Matrix(3);

	t2_p = new Matrix*[n_solutions];
	for (int i = 0; i < n_solutions; i++)
		t2_p[i] = new Matrix(3);

	n_p = new Matrix*[n_solutions];
	for (int i = 0; i < n_solutions; i++)
		n_p[i] = new Matrix(3);

	G_i = new Matrix*[n_solutions];
	for (int i = 0; i < n_solutions; i++)
		G_i[i] = new Matrix(3);

	G_ip = new Matrix*[n_solutions];
	for (int i = 0; i < n_solutions; i++)
		G_ip[i] = new Matrix(3);

	t1_i = new Matrix*[n_solutions];
	for (int i = 0; i < n_solutions; i++)
		t1_i[i] = new Matrix(3);

	t2_i = new Matrix*[n_solutions];
	for (int i = 0; i < n_solutions; i++)
		t2_i[i] = new Matrix(3);

	n_i = new Matrix*[n_solutions];
	for (int i = 0; i < n_solutions; i++)
		n_i[i] = new Matrix(3);
}

SurfaceData::~SurfaceData()
{
	if (G_p != NULL)
	{
		for (int i = 0; i < n_sol; i++)
			delete G_p[i];
		delete[]G_p;
	}
	if (G_i != NULL)
	{
		for (int i = 0; i < n_sol; i++)
			delete G_i[i];
		delete[]G_i;
	}
	if (G_ip != NULL)
	{
		for (int i = 0; i < n_sol; i++)
			delete G_ip[i];
		delete[]G_ip;
	}
	if (t1_p != NULL)
	{
		for (int i = 0; i < n_sol; i++)
			delete t1_p[i];
		delete[]t1_p;
	}
	if (t2_p != NULL)
	{
		for (int i = 0; i < n_sol; i++)
			delete t2_p[i];
		delete[]t2_p;
	}
	if (n_p != NULL)
	{
		for (int i = 0; i < n_sol; i++)
			delete n_p[i];
		delete[]n_p;
	}
	if (t1_i != NULL)
	{
		for (int i = 0; i < n_sol; i++)
			delete t1_i[i];
		delete[]t1_i;
	}
	if (t2_i != NULL)
	{
		for (int i = 0; i < n_sol; i++)
			delete t2_i[i];
		delete[]t2_i;
	}
	if (n_i != NULL)
	{
		for (int i = 0; i < n_sol; i++)
			delete n_i[i];
		delete[]n_i;
	}
}