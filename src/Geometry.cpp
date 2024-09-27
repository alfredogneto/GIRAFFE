#include "Geometry.h"

#include"Database.h"

//extern
//FILE *fdebug;

//Variaveis globais
extern
Database db;

Geometry::Geometry()
{
}

Geometry::~Geometry()
{
}

void Geometry::Alloc()
{
	number = 0;
	material = 0;
	super_node = 0;
	mother_entity = 0;
	type_name = new char[20];
	nodes = new int[n_nodes];
	DOFs = new int *[n_nodes];
	for (int i = 0; i < n_nodes; i++)
		DOFs[i] = new int[db.number_GLs_node];
	GLs = new int*[nDOFs];
	for (int i = 0; i < nDOFs; i++)
		GLs[i] = NULL;

	//Degeneration
	deg_u1 = false;
	deg_u2 = false;
	deg_u1_u2 = false;
	deg_u1_value = 0.0;
	deg_u2_value = 0.0;
	deg_u1_u2_relation = 0.0;
}

void Geometry::Free()
{
	delete[] type_name;
	delete[] nodes;
	for (int i = 0; i < n_nodes; i++)
		delete[] DOFs[i];
	delete[] DOFs;
	/*for (int i = 0; i < nDOFs; i++)
		delete[] GLs;*/
	delete[]GLs;
}
