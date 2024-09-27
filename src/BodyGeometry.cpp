#include "BodyGeometry.h"

#include "Geometry.h"
#include "BoundingBoxAxesAligned.h"
#include "GeneralContactSearch.h"
#include"Database.h"
//Variaveis globais
extern
Database db;


BodyGeometry::BodyGeometry()
{
	n_items = 0;
	list_items = NULL;

	sequence = false;
	list = false;
	initial = 0;
	increment = 0;

	bv = new BoundingBoxAxesAligned();
	ptr_geom = NULL;
	max_offset = 0.0f;

	mass = 1.0;	//default value
}

BodyGeometry::~BodyGeometry()
{
	if (list_items != NULL)
	{
		delete[] list_items;
	}
	delete bv;

	if (ptr_geom != NULL)
	{
		delete[] ptr_geom;
	}
}

bool BodyGeometry::Read(FILE *f)
{
	char s[1000];
	fscanf(f, "%s", s);
	if (!strcmp(s, "BodyGeometry"))
	{
		fscanf(f, "%s", s);
		number = atoi(s);
	}
	else
		return false;

	fscanf(f, "%s", s);
	//Optional keyword: Mass
	if (!strcmp(s, "Mass"))
	{
		fscanf(f, "%s", s);
		mass = atof(s);
		fscanf(f, "%s", s);
	}
	

	if (!strcmp(s, "Geometries"))
	{
		fscanf(f, "%s", s);
		n_items = atoi(s);
		//Alocação do vetor
		list_items = new int[n_items];
	}
	else
		return false;
	//Duas possibilidades de leitura:
	//1 - List
	//2 - Sequence
	fscanf(f, "%s", s);
	if (!strcmp(s, "List"))
	{
		list = true;
		for (int i = 0; i < n_items; i++)
		{
			fscanf(f, "%s", s);//Leitura do numero do nó
			list_items[i] = atoi(s);
		}
	}
	else
	{
		if (!strcmp(s, "Sequence"))
		{
			sequence = true;
			fscanf(f, "%s", s);
			if (!strcmp(s, "Initial"))
			{
				fscanf(f, "%s", s);
				initial = atoi(s);
			}
			else
				return false;
			fscanf(f, "%s", s);
			if (!strcmp(s, "Increment"))
			{
				fscanf(f, "%s", s);
				increment = atoi(s);
			}
			else
				return false;
			//Geração da lista de nós
			for (int i = 0; i < n_items; i++)
			{
				list_items[i] = initial + i * increment;
			}
		}
		else
			return false;
	}
	//Se atingiu esse ponto, sinal de leitura correta de tudo: retorna true
	return true;
}

void BodyGeometry::Write(FILE *f)
{
	fprintf(f, "BodyGeometry\t%d\tGeometries\t%d\tList\t", number, n_items);
	for (int i = 0; i < n_items; i++)
		fprintf(f, "%d\t", list_items[i]);
	fprintf(f, "\n");
}

void BodyGeometry::WriteVTK_XMLRender(FILE *f)
{
	for (int i = 0; i < n_items; i++)
		db.geometries[list_items[i] - 1]->WriteVTK_XMLRender(f);
}

bool BodyGeometry::Check()
{
	for (int i = 0; i < n_items; i++)
	{
		if (list_items[i] > db.number_geometries)
			return false;
	}
	return true;
}

void BodyGeometry::PreCalc()
{
	//Pointer to each Geometry
	ptr_geom = new Geometry*[n_items];
	for (int i = 0; i < n_items; i++)
		ptr_geom[i] = db.geometries[list_items[i] - 1];

	max_offset = 0.0f;
	for (int i = 0; i < n_items; i++)
	{
		db.geometries[list_items[i] - 1]->PreCalc();
		db.geometries[list_items[i] - 1]->body_mass = &mass;	//pointing to the body mass
		db.geometries[list_items[i] - 1]->bv->associated_ID = number;
		//Associated entity:
		db.geometries[list_items[i] - 1]->bv->associated_type = 'G';
		db.geometries[list_items[i] - 1]->bv->associated_ID = number;
		db.geometries[list_items[i] - 1]->bv->associated_sub_ID = i;//zero-based
		//calcula maximo offset por parametro de tamanho, raio e/ou espessura
		//avaliado por varredura de BV's de geometries, coletando informações - salva na variavel max_offset - 
		//OBS: ja inclui o inc_len_factor
		if (db.geometries[list_items[i] - 1]->bv_offset > max_offset)
			max_offset = db.geometries[list_items[i] - 1]->bv_offset;
	}
	
	if (db.gcs_exist)
		inc_len_factor = db.gcs->inc_len_factor;

	BoundingBoxAxesAligned* ptr_bv = static_cast<BoundingBoxAxesAligned*>(bv);
	//Associated entity:
	ptr_bv->associated_type = 'G';
	ptr_bv->associated_ID = number;
	ptr_bv->associated_sub_ID = 0;
	ptr_bv->inc_len_factor = inc_len_factor;

	UpdateVariables();
	UpdateBoundingVolumes();
	SaveLagrange();

	//TO DO
	//Evaluate the body geometry mass
}

void BodyGeometry::UpdateVariables()
{
	for (int i = 0; i < n_items; i++)
		db.geometries[list_items[i] - 1]->UpdateVariables();
}

void BodyGeometry::UpdateBoundingVolumes()
{
	BoundingBoxAxesAligned* ptr_bv = static_cast<BoundingBoxAxesAligned*>(bv);
	ptr_bv->x_min = 0.0f;
	ptr_bv->x_max = 0.0f;

	ptr_bv->y_min = 0.0f;
	ptr_bv->y_max = 0.0f;

	ptr_bv->z_min = 0.0f;
	ptr_bv->z_max = 0.0f;

	//Update geometry BV's
	for (int i = 0; i < n_items; i++)
	{
		db.geometries[list_items[i] - 1]->UpdateBoundingVolumes();
		//Update the own BV - based on each geometry BV and max_offset (already set on PreCalc)
		if (ptr_bv->x_max < (db.geometries[list_items[i] - 1]->max[0] + db.geometries[list_items[i] - 1]->bv_offset))
			ptr_bv->x_max = db.geometries[list_items[i] - 1]->max[0] + db.geometries[list_items[i] - 1]->bv_offset;
		if (ptr_bv->x_min > (db.geometries[list_items[i] - 1]->min[0] - db.geometries[list_items[i] - 1]->bv_offset))
			ptr_bv->x_min = db.geometries[list_items[i] - 1]->min[0] - db.geometries[list_items[i] - 1]->bv_offset;

		if (ptr_bv->y_max < (db.geometries[list_items[i] - 1]->max[1] + db.geometries[list_items[i] - 1]->bv_offset))
			ptr_bv->y_max = db.geometries[list_items[i] - 1]->max[1] + db.geometries[list_items[i] - 1]->bv_offset;
		if (ptr_bv->y_min > (db.geometries[list_items[i] - 1]->min[1] - db.geometries[list_items[i] - 1]->bv_offset))
			ptr_bv->y_min = db.geometries[list_items[i] - 1]->min[1] - db.geometries[list_items[i] - 1]->bv_offset;

		if (ptr_bv->z_max < (db.geometries[list_items[i] - 1]->max[2] + db.geometries[list_items[i] - 1]->bv_offset))
			ptr_bv->z_max = db.geometries[list_items[i] - 1]->max[2] + db.geometries[list_items[i] - 1]->bv_offset;
		if (ptr_bv->z_min > (db.geometries[list_items[i] - 1]->min[2] - db.geometries[list_items[i] - 1]->bv_offset))
			ptr_bv->z_min = db.geometries[list_items[i] - 1]->min[2] - db.geometries[list_items[i] - 1]->bv_offset;
	}
	//Updating center and size (used in Verlet/LinkedCells schemes)
	ptr_bv->size = 0.5f * sqrt((ptr_bv->x_max - ptr_bv->x_min) * (ptr_bv->x_max - ptr_bv->x_min) + 
		(ptr_bv->y_max - ptr_bv->y_min) * (ptr_bv->y_max - ptr_bv->y_min) + 
		(ptr_bv->z_max - ptr_bv->z_min) * (ptr_bv->z_max - ptr_bv->z_min));
	ptr_bv->x_center[0] = 0.5f * (ptr_bv->x_max + ptr_bv->x_min);
	ptr_bv->x_center[1] = 0.5f * (ptr_bv->y_max + ptr_bv->y_min);
	ptr_bv->x_center[2] = 0.5f * (ptr_bv->z_max + ptr_bv->z_min);
}

void BodyGeometry::SaveLagrange()
{
	for (int i = 0; i < n_items; i++)
		db.geometries[list_items[i] - 1]->SaveLagrange();
}
