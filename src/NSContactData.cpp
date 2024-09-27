#include "NSContactData.h"

#include "Matrix.h"

NSContactData::NSContactData()
{
	n_solutions = 1;

	convective = NULL;
	copy_convective = NULL;
	return_value = NULL;
	repeated = NULL;
	g_n = NULL;
	copy_g_n = NULL;
	g_t = NULL;
	copy_g_t = NULL;
	G_p = NULL;
	t1_p = NULL;
	t2_p = NULL;
	n_p = NULL;
	G_i = NULL;
	G_ip = NULL;
	t1_i = NULL;
	t2_i = NULL;
	n_i = NULL;
	alloced = false;
}

NSContactData::~NSContactData()
{
	Free();
}

//Aloca matrizes
void NSContactData::Alloc()
{
	if (alloced == false)
	{
		convective = new double*[n_solutions];
		for (int i = 0; i < n_solutions; i++)
			convective[i] = new double[2];
		copy_convective = new double*[n_solutions];
		for (int i = 0; i < n_solutions; i++)
			copy_convective[i] = new double[2];
		return_value = new int[n_solutions];
		repeated = new bool[n_solutions];

		g_n = new double[n_solutions];

		copy_g_n = new double[n_solutions];

		g_t = new Matrix*[n_solutions];
		for (int i = 0; i < n_solutions; i++)
			g_t[i] = new Matrix(3);

		copy_g_t = new Matrix*[n_solutions];
		for (int i = 0; i < n_solutions; i++)
			copy_g_t[i] = new Matrix(3);

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

		for (int i = 0; i < n_solutions; i++)
		{
			copy_convective[i][0] = 0.0;
			copy_convective[i][1] = 0.0;
			convective[i][0] = 0.0;
			convective[i][1] = 0.0;
			g_n[i] = 0.0;
			copy_g_n[i] = 1.0;	//Assumindo que não ha contato anterior
			return_value[i] = 1;
			repeated[i] = true;
		}
		alloced = true;
	}
}

//Desaloca matrizes
void NSContactData::Free()
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
		if (return_value != NULL)
			delete[]return_value;
		if (repeated != NULL)
			delete[]repeated;
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
		if (G_p != NULL)
		{
			for (int i = 0; i < n_solutions; i++)
				delete G_p[i];
			delete[]G_p;
		}
		if (G_i != NULL)
		{
			for (int i = 0; i < n_solutions; i++)
				delete G_i[i];
			delete[]G_i;
		}
		if (G_ip != NULL)
		{
			for (int i = 0; i < n_solutions; i++)
				delete G_ip[i];
			delete[]G_ip;
		}
		if (t1_p != NULL)
		{
			for (int i = 0; i < n_solutions; i++)
				delete t1_p[i];
			delete[]t1_p;
		}
		if (t2_p != NULL)
		{
			for (int i = 0; i < n_solutions; i++)
				delete t2_p[i];
			delete[]t2_p;
		}
		if (n_p != NULL)
		{
			for (int i = 0; i < n_solutions; i++)
				delete n_p[i];
			delete[]n_p;
		}
		if (t1_i != NULL)
		{
			for (int i = 0; i < n_solutions; i++)
				delete t1_i[i];
			delete[]t1_i;
		}
		if (t2_i != NULL)
		{
			for (int i = 0; i < n_solutions; i++)
				delete t2_i[i];
			delete[]t2_i;
		}
		if (n_i != NULL)
		{
			for (int i = 0; i < n_solutions; i++)
				delete n_i[i];
			delete[]n_i;
		}
		alloced = false;
	}

	convective = NULL;
	copy_convective = NULL;
	return_value = NULL;
	repeated = NULL;
	g_n = NULL;
	copy_g_n = NULL;
	g_t = NULL;
	copy_g_t = NULL;
	G_p = NULL;
	t1_p = NULL;
	t2_p = NULL;
	n_p = NULL;
	G_i = NULL;
	G_ip = NULL;
	t1_i = NULL;
	t2_i = NULL;
	n_i = NULL;
}

//Checa repetição de raizes e salva a informação na matriz 'repeated'
void NSContactData::CheckRepeated(double tol_coordinate_value)
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
void NSContactData::Plot()
{
	for (int i = 0; i < n_solutions; i++)
	{
		printf("convective %.6e  %.6e  %d\n", convective[i][0], convective[i][1],repeated[i]);
	}
}

