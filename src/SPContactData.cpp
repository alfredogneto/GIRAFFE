#include "SPContactData.h"


#include "Matrix.h"
#include "SplineElementData.h"
#include "Database.h"
//Variáveis globais
extern
Database db;


////////////////////////////////////////////////////////////////////
#define _CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include <crtdbg.h>
#ifdef _DEBUG
#define DBG_NEW new ( _NORMAL_BLOCK , __FILE__ , __LINE__ )
// Replace _NORMAL_BLOCK with _CLIENT_BLOCK if you want the
// allocations to be of _CLIENT_BLOCK type
#else
#define DBG_NEW new
#endif
////////////////////////////////////////////////////////////////////

SPContactData::SPContactData()
{
	n_solutions = 1;
	alloced = false;
	convective = NULL;
	copy_convective = NULL;
	copy_deg_coordinates = NULL;
	initial_guess = NULL;
	return_value = NULL;
	copy_return_value = NULL;
	repeated = NULL;
	g_n = NULL;
	copy_g_n = NULL;
	g_t = NULL;
	copy_g_t = NULL;
	g = NULL;
	n = NULL;
	copy_g = NULL;
	copy_n = NULL;
	surf1 = NULL;
	surf2 = NULL;
	invHessian = NULL;
	P_0 = NULL;
	P = NULL;
	copy_degenerated = NULL;
	degenerated = NULL;
	deg_control = NULL;
	stick = NULL;
	copy_stick = NULL;

}

SPContactData::~SPContactData()
{
	Free();
}

void SPContactData::MountDegenerativeOperator()
{
	for (int i = 0; i < n_solutions; i++)
	{
		int counter = 0;
		for (int j = 0; j < 2; j++)
		{
			if (deg_control[i][j] == false)
				counter++;
		}
		P_0[i]->setColumns(counter);
		P_0[i]->setLines(2);
		P_0[i]->alloc();
		counter = 0;
		for (int j = 0; j < 2; j++)
		{
			if (deg_control[i][j] == false)
			{
				(*P_0[i])(0, counter) = (*P[i])(0, j);
				(*P_0[i])(1, counter) = (*P[i])(1, j);
				counter++;
			}
		}
	}
}

//Aloca matrizes
void SPContactData::Alloc()
{
	if (alloced == false)
	{
		convective = DBG_NEW double*[n_solutions];
		for (int i = 0; i < n_solutions; i++)
			convective[i] = new double[2];
		copy_convective = new double*[n_solutions];
		for (int i = 0; i < n_solutions; i++)
			copy_convective[i] = new double[2];
		copy_deg_coordinates = new double*[n_solutions];
		for (int i = 0; i < n_solutions; i++)
			copy_deg_coordinates[i] = new double[2];
		initial_guess = new double*[n_solutions];
		for (int i = 0; i < n_solutions; i++)
			initial_guess[i] = new double[2];
		deg_control = new bool*[n_solutions];
		for (int i = 0; i < n_solutions; i++)
			deg_control[i] = new bool[2];
		return_value = new int[n_solutions];
		copy_return_value = new int[n_solutions];
		repeated = new bool[n_solutions];
		copy_degenerated = new bool[n_solutions];
		degenerated = new bool[n_solutions];

		g_n = new double[n_solutions];

		copy_g_n = new double[n_solutions];

		g_t = new Matrix*[n_solutions];
		for (int i = 0; i < n_solutions; i++)
			g_t[i] = new Matrix(3);

		copy_g_t = new Matrix*[n_solutions];
		for (int i = 0; i < n_solutions; i++)
			copy_g_t[i] = new Matrix(3);

		g = new Matrix*[n_solutions];
		for (int i = 0; i < n_solutions; i++)
			g[i] = new Matrix(3);

		copy_g = new Matrix*[n_solutions];
		for (int i = 0; i < n_solutions; i++)
			copy_g[i] = new Matrix(3);

		n = new Matrix*[n_solutions];
		for (int i = 0; i < n_solutions; i++)
			n[i] = new Matrix(3);

		copy_n = new Matrix*[n_solutions];
		for (int i = 0; i < n_solutions; i++)
			copy_n[i] = new Matrix(3);

		surf1 = new SplineElementData(n_solutions);
		surf2 = new SplineElementData(n_solutions);

		invHessian = new double**[n_solutions];
		for (int i = 0; i < n_solutions; i++)
		{
			invHessian[i] = new double*[2];
			for (int j = 0; j < 2; j++)
				invHessian[i][j] = new double[2];
		}

		P_0 = new Matrix*[n_solutions];
		for (int i = 0; i < n_solutions; i++)
			P_0[i] = new Matrix(2, 1);
		P = new Matrix*[n_solutions];
		for (int i = 0; i < n_solutions; i++)
			P[i] = new Matrix(2, 2);

		stick = new bool[n_solutions];
		copy_stick = new bool[n_solutions];

		for (int i = 0; i < n_solutions; i++)
		{

			copy_convective[i][0] = 0.0;
			copy_convective[i][1] = 0.0;
			degenerated[i] = false;
			copy_degenerated[i] = false;
			deg_control[i][0] = false;
			deg_control[i][1] = false;
			convective[i][0] = 0.0;
			convective[i][1] = 0.0;
			copy_deg_coordinates[i][0] = 0.0;
			copy_deg_coordinates[i][1] = 0.0;
			initial_guess[i][0] = 0.0;
			initial_guess[i][1] = 0.0;
			g_n[i] = 0.0;
			copy_g_n[i] = 1.0;	//Assumindo que não há contato anterior
			return_value[i] = 2;
			copy_return_value[i] = 2;
			repeated[i] = true;
			stick[i] = true;
			copy_stick[i] = true;


			for (int ni = 0; ni < 2; ni++)
				for (int nj = 0; nj < 2; nj++)
					invHessian[i][ni][nj] = 0.0;

		}
		alloced = true;
	}
}

//Desaloca matrizes
void SPContactData::Free()
{
	if (alloced == true)
	{
		if (convective != NULL)
		{
			for (int i = 0; i < n_solutions; i++)
				delete[] convective[i];
			delete[]convective;
		}
		if (copy_convective != NULL)
		{
			for (int i = 0; i < n_solutions; i++)
				delete[] copy_convective[i];
			delete[]copy_convective;
		}
		if (copy_deg_coordinates != NULL)
		{
			for (int i = 0; i < n_solutions; i++)
				delete[] copy_deg_coordinates[i];
			delete[]copy_deg_coordinates;
		}
		if (initial_guess != NULL)
		{
			for (int i = 0; i < n_solutions; i++)
				delete[] initial_guess[i];
			delete[]initial_guess;
		}
		if (deg_control != NULL)
		{
			for (int i = 0; i < n_solutions; i++)
				delete[] deg_control[i];
			delete[]deg_control;
		}
		if (return_value != NULL)
			delete[]return_value;
		if (copy_return_value != NULL)
			delete[]copy_return_value;
		if (repeated != NULL)
			delete[]repeated;
		if (copy_degenerated != NULL)
			delete[]copy_degenerated;
		if (degenerated != NULL)
			delete[]degenerated;
		if (g_n != NULL)
			delete[]g_n;
		if (copy_g_n != NULL)
			delete[]copy_g_n;
		if (g_t != NULL)
		{
			for (int i = 0; i < n_solutions; i++)
				delete g_t[i];
			delete[]g_t;
		}
		if (copy_g_t != NULL)
		{
			for (int i = 0; i < n_solutions; i++)
				delete copy_g_t[i];
			delete[]copy_g_t;
		}
		if (g != NULL)
		{
			for (int i = 0; i < n_solutions; i++)
				delete g[i];
			delete[]g;
		}
		if (copy_g != NULL)
		{
			for (int i = 0; i < n_solutions; i++)
				delete copy_g[i];
			delete[]copy_g;
		}
		if (n != NULL)
		{
			for (int i = 0; i < n_solutions; i++)
				delete n[i];
			delete[]n;
		}
		if (copy_n != NULL)
		{
			for (int i = 0; i < n_solutions; i++)
				delete copy_n[i];
			delete[]copy_n;
		}
		if (surf1 != NULL)
			delete surf1;
		if (surf2 != NULL)
			delete surf2;

		if (invHessian != NULL)
		{
			for (int i = 0; i < n_solutions; i++)
			{
				for (int j = 0; j < 2; j++)
				{
					delete[] invHessian[i][j];
				}
				delete[]invHessian[i];
			}
			delete[]invHessian;
		}
		alloced = false;

		if (P_0 != NULL)
		{
			for (int i = 0; i < n_solutions; i++)
				delete P_0[i];
			delete[]P_0;
		}

		if (P != NULL)
		{
			for (int i = 0; i < n_solutions; i++)
				delete P[i];
			delete[]P;
		}

		if (stick != NULL)
			delete[]stick;

		if (copy_stick != NULL)
			delete[]copy_stick;
	}
	convective = NULL;
	copy_convective = NULL;
	copy_deg_coordinates = NULL;
	initial_guess = NULL;
	deg_control = NULL;
	return_value = NULL;
	copy_return_value = NULL;
	repeated = NULL;
	g_n = NULL;
	copy_g_n = NULL;
	g_t = NULL;
	copy_g_t = NULL;
	g = NULL;
	n = NULL;
	copy_g = NULL;
	copy_n = NULL;
	surf1 = NULL;
	surf2 = NULL;
	invHessian = NULL;
	P_0 = NULL;
	P = NULL;
	copy_degenerated = NULL;
	degenerated = NULL;
	stick = NULL;
	copy_stick = NULL;
}

//Checa repetição de raízes e salva a informação na matriz 'repeated'
void SPContactData::CheckRepeated(double tol_coordinate_value)
{
	for (int i = 0; i < n_solutions; i++)
	{
		repeated[i] = false;
		for (int j = 0; j < i; j++)
		{
			if (return_value[i] != 1)//se não houve divergência
			{
				if (abs(convective[i][0] - convective[j][0]) <= tol_coordinate_value && abs(convective[i][1] - convective[j][1]) <= tol_coordinate_value)
					repeated[i] = true;
			}
			else
				repeated[i] = true;
		}
	}
}
void SPContactData::Plot()
{
	for (int i = 0; i < n_solutions; i++)
	{
		db.myprintf("convective %.6e  %.6e  %.6e  %.6e  return %d\n", convective[i][0], convective[i][1], return_value[i]);
	}
}
//Using info from deg_control and P, establishes P_0 by selecting appropriate columns
//void SPContactData::MountDegenerativeOperator()
//{
//	for (int i = 0; i < n_solutions; i++)
//	{
//		int counter = 0;
//		for (int j = 0; j < 2; j++)
//		{
//			if (deg_control[i][j] == false)
//				counter++;
//		}
//		P_0[i]->setColumns(counter);
//		P_0[i]->setLines(2);
//		P_0[i]->alloc();
//		counter = 0;
//		for (int j = 0; j < 2; j++)
//		{
//			if (deg_control[i][j] == false)
//			{
//				(*P_0[i])(0, counter) = (*P[i])(0, j);
//				(*P_0[i])(1, counter) = (*P[i])(1, j);
//				counter++;
//			}
//		}
//	}
//}

