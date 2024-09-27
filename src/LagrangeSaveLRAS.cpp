#include "LagrangeSaveLRAS.h"

#include "Matrix.h"

LagrangeSaveLRAS::LagrangeSaveLRAS(int e_elements)
{
	n_elements = e_elements;
	contact_status = new int*[n_elements];
	for (int i = 0; i < n_elements; i++)
		contact_status[i] = new int[3];
	gt1s = new double*[n_elements];
	for (int i = 0; i < n_elements; i++)
		gt1s[i] = new double[3];
	gt2s = new double*[n_elements];
	for (int i = 0; i < n_elements; i++)
		gt2s[i] = new double[3];
	x0prev = new Matrix**[n_elements];
	for (int i = 0; i < n_elements; i++)
	{
		x0prev[i] = new Matrix*[3];
		for (int j = 0; j < 3; j++)
			x0prev[i][j] = new Matrix(3);
	}
	gn = new double*[n_elements];
	for (int i = 0; i < n_elements; i++)
		gn[i] = new double[3];
	Fx = new double*[n_elements];
	for (int i = 0; i < n_elements; i++)
		Fx[i] = new double[3];
	Fy = new double*[n_elements];
	for (int i = 0; i < n_elements; i++)
		Fy[i] = new double[3];
	Fz = new double*[n_elements];
	for (int i = 0; i < n_elements; i++)
		Fz[i] = new double[3];

	for (int i = 0; i < n_elements; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			contact_status[i][j] = 0;
			gt1s[i][j] = 0;
			gt2s[i][j] = 0;
			gn[i][j] = 0;
			Fx[i][j] = 0;
			Fy[i][j] = 0;
			Fz[i][j] = 0;
		}
	}
}

LagrangeSaveLRAS::~LagrangeSaveLRAS()
{

	for (int i = 0; i < n_elements; i++)
		delete contact_status[i];
	delete[] contact_status;
	for (int i = 0; i < n_elements; i++)
		delete gt1s[i];
	delete[] gt1s;
	for (int i = 0; i < n_elements; i++)
		delete gt2s[i];
	delete[] gt2s;
	for (int i = 0; i < n_elements; i++)
	{
		for (int j = 0; j < 3; j++)
			delete x0prev[i][j];
		delete[] x0prev[i];
	}
	delete[] x0prev;
	for (int i = 0; i < n_elements; i++)
		delete gn[i];
	delete[] gn;
	for (int i = 0; i < n_elements; i++)
		delete Fx[i];
	delete[] Fx;
	for (int i = 0; i < n_elements; i++)
		delete Fy[i];
	delete[] Fy;
	for (int i = 0; i < n_elements; i++)
		delete Fz[i];
	delete[] Fz;
	
}
