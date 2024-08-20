#include "Element.h"
#define PI 3.1415926535897932384626433832795
#include"Database.h"
//Variáveis globais
extern
Database db;

Element::Element()
{
}

Element::~Element()
{
}
//Calcula a banda gerada na matriz global pelo elemento
void Element::Band(int* band_fixed, int* band_free)
{
	temp_band_free = 0;
	temp_band_fixed = 0;

	lowest_free_global_DOF = db.number_nodes*db.number_GLs_node;
	highest_free_global_DOF = 0;
	lowest_fixed_global_DOF = db.number_nodes*db.number_GLs_node;
	highest_fixed_global_DOF = 0;
	free_marked = false;
	fixed_marked = false;

	//Percorre os nós para identificar graus de liberdade globais e determinar a banda
	for (int i = 0; i < n_nodes; i++)
	{
		for (int j = 0; j < db.number_GLs_node; j++)
		{
			//se o GL for livre e ativo
			if (db.nodes[nodes[i] - 1]->constraints[j] == 0 && db.nodes[nodes[i] - 1]->active_GL[j] == 1)
			{
				if (db.nodes[nodes[i] - 1]->GLs[j] < lowest_free_global_DOF)
					lowest_free_global_DOF = db.nodes[nodes[i] - 1]->GLs[j];
				if (db.nodes[nodes[i] - 1]->GLs[j] > highest_free_global_DOF)
					highest_free_global_DOF = db.nodes[nodes[i] - 1]->GLs[j];
				free_marked = true;
			}

			//se o GL for fixo e ativo
			if (db.nodes[nodes[i] - 1]->constraints[j] == 1 && db.nodes[nodes[i] - 1]->active_GL[j] == 1)
			{
				if (-db.nodes[nodes[i] - 1]->GLs[j] < lowest_fixed_global_DOF)
					lowest_fixed_global_DOF = -db.nodes[nodes[i] - 1]->GLs[j];
				if (-db.nodes[nodes[i] - 1]->GLs[j] > highest_fixed_global_DOF)
					highest_fixed_global_DOF = -db.nodes[nodes[i] - 1]->GLs[j];
				fixed_marked = true;
			}
		}
	}
	if (free_marked == true)
		temp_band_free = highest_free_global_DOF - lowest_free_global_DOF;
	else
		temp_band_free = 0;
	if (fixed_marked == true)
		temp_band_fixed = highest_fixed_global_DOF - lowest_fixed_global_DOF;
	else
		temp_band_fixed = 0;
	*band_fixed = temp_band_fixed;
	*band_free = temp_band_free;
}
