#include "Solid_1.h"
#include "Beam_1.h"
#include "Encoding.h"	

#include"Database.h"
//Vari�veis globais
extern
Database db;

Solid_1::Solid_1()
{
	strain_energy = 0.0;
	kinetic_energy = 0.0;
	potential_gravitational_energy = 0.0;

	//VTK_type = 24;
	nDOFs = 24;
	material = 0;
	section = 0;
	n_nodes = 8;
	number = 0;
	nodes = new int[n_nodes];
	VTK_nodes = new int[n_nodes];
	VTK_nodes[0] = 0;
	VTK_nodes[1] = 1;
	VTK_nodes[2] = 2;
	VTK_nodes[3] = 3;
	VTK_nodes[4] = 4;
	VTK_nodes[5] = 5;
	VTK_nodes[6] = 6;
	VTK_nodes[7] = 7;
	
	DOFs = new int *[n_nodes];
	for (int i = 0; i < n_nodes;i++)
		DOFs[i] = new int[db.number_GLs_node];

	//Rotina para ativar os GLS de cada n� do elemento
	for (int i = 0; i < n_nodes; i++)
	{
		for (int j = 0; j < db.number_GLs_node; j++)
		{
			DOFs[i][j] = 0;
		}
		DOFs[i][0] = 1;
		DOFs[i][1] = 1;
		DOFs[i][2] = 1;
	}
}

Solid_1::~Solid_1()
{
	delete[] nodes;
	delete[] VTK_nodes;
	if (DOFs != NULL)
	{
		for (int i = 0; i < n_nodes; i++)
			delete[] DOFs[i];
		delete[] DOFs;
	}
}

//Checa inconsist�ncias no elemento para evitar erros de execu��o
bool Solid_1::Check()
{
	return true;
}

bool Solid_1::Read(FILE *f)
{
	char s[1000];
	fscanf(f, "%s", s);
	number = atoi(s);

	fscanf(f, "%s", s);
	if (!strcmp(s, "Mat"))
	{
		fscanf(f, "%s", s);
		material = atoi(s);
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "CS"))
	{
		fscanf(f, "%s", s);
		cs = atoi(s);
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "Nodes"))
	{
		for (int n = 0; n < n_nodes; n++)
		{
			fscanf(f, "%s", s);
			nodes[n] = atoi(s);
		}
	}
	else
		return false;
	return true;
}

void Solid_1::Write(FILE *f)
{
	fprintf(f, "Solid_1\t%d\tMat\t%d\tCS\t%d\tNodes\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
		number,
		material,
		cs,
		nodes[0],
		nodes[1],
		nodes[2],
		nodes[3], 
		nodes[4], 
		nodes[5], 
		nodes[6], 
		nodes[7] );
}
//Escreve arquivo de resultados
void Solid_1::WriteResults(FILE *f)
{

}

void Solid_1::WriteVTK_XMLBase(std::vector<float> *float_vector)
{
	//Imprime os resultados do elemento
	int res_element = 0;
	//Imprime valores nulos para resultados que n�o fazem sentido para esse tipo de elemento
	for (int i = res_element; i < db.n_element_results; i++)
		float_vector->push_back(0.0);
}

void Solid_1::WriteVTK_XMLRender(FILE *f)
{
	//DOES NOTHING
}

//Escreve no monitor do elemento//Escreve no monitor do elemento
void Solid_1::WriteMonitor(FILE *f, bool first_record, double time)
{

}

//Monta elementos
void Solid_1::Mount()
{

}
//Monta carregamentos associados ao elemento
void Solid_1::MountElementLoads()
{

}
//Monta matriz de transforma��o de coordenadas
void Solid_1::TransformMatrix()
{

}
//Preenche a contribui��o do elemento nas matrizes globais
void Solid_1::MountGlobal()
{

}
//Salva vari�veis nos pontos de Gauss �teis para descri��o lagrangiana atualizada
void Solid_1::SaveLagrange()
{

}
//Pr�-c�lculo de vari�veis que � feito uma �nica vez no in�cio
void Solid_1::PreCalc()
{

}

//Monta a matriz de massa
void Solid_1::MountMass()
{

}

//Monta a matriz de massa
void Solid_1::MountMassModal()
{
	Zeros();
}

//Monta a matriz de amortecimento para realiza��o da an�lise modal
void Solid_1::MountDampingModal()
{
	Zeros();
}

//Monta a matriz de amortecimento
void Solid_1::MountDamping(bool update_rayleigh)
{

}

//Montagens - Newmark
void Solid_1::MountDyn()
{

}

//Montagens para an�lise modal - inser��o da matriz de massa e amortecimento na matriz de rigidez para posterior montagem global
void Solid_1::MountDynModal()
{
	
}

//Zera matrizes locais do elemento
void Solid_1::Zeros()
{
	kinetic_energy = 0.0;
	strain_energy = 0.0;
	potential_gravitational_energy = 0.0;
}