#include "NodeSet.h"

#include"Database.h"
//Variáveis globais
extern
Database db;

NodeSet::NodeSet()
{
	n_nodes = 0;
	node_list = NULL;

	sequence = false;
	list = false;
	initial = 0;
	increment = 0;
}

NodeSet::~NodeSet()
{
	if (node_list != NULL)
	{
		delete[] node_list;
	}
}

bool NodeSet::Read(FILE *f)
{
	char s[1000];
	fscanf(f, "%s", s);
	//Verifica a palavra chave "Set"
	if (!strcmp(s, "NodeSet"))
	{
		fscanf(f, "%s", s);
		number = atoi(s);
	}
	else
		return false;
	fscanf(f, "%s", s);
	//Verifica a palavra chave "Nodes"
	if (!strcmp(s, "Nodes"))
	{
		fscanf(f, "%s", s);
		n_nodes = atoi(s);
		//Alocação do vetor de nós
		node_list = new int[n_nodes];
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
		for (int i = 0; i < n_nodes; i++)
		{
			fscanf(f, "%s", s);//Leitura do número do nó
			node_list[i] = atoi(s);
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
			for (int i = 0; i < n_nodes; i++)
			{
				node_list[i] = initial+i*increment;
			}
		}
		else
			return false;
	}
	//Se atingiu esse ponto, sinal de leitura correta de tudo: retorna true
	return true;
}

void NodeSet::Write(FILE *f)
{
	fprintf(f, "NodeSet\t%d\tNodes\t%d\tList\t", number, n_nodes);
	for (int i = 0; i < n_nodes; i++)
		fprintf(f, "%d\t",node_list[i]);
	fprintf(f, "\n");
}

void NodeSet::WriteMonitor(FILE *f, bool first_record, double time)
{
	//Variáveis do NodeSet para Monitor
	Matrix coordinates(3);
	Matrix rot_euler(3);
	Matrix force(3);
	Matrix moment(3);
	Matrix temp_vec(3);
	Matrix temp_coordinates(3);
	Matrix temp_force(3);
	//Percorre todos os nós de dentro do "set", salvando alguns resultados
	//Cálculo do ângulo de rotação
	int current_node;
	Node *tempnode;
	
	for (int i = 0; i < n_nodes; i++)
	{
		current_node = node_list[i];
		tempnode = db.nodes[current_node - 1];

		coordinates(0, 0) = coordinates(0, 0) + tempnode->copy_coordinates[0];
		coordinates(1, 0) = coordinates(1, 0) + tempnode->copy_coordinates[1];
		coordinates(2, 0) = coordinates(2, 0) + tempnode->copy_coordinates[2];

		rot_euler(0, 0) = rot_euler(0, 0) + tempnode->copy_rot_euler[0];
		rot_euler(1, 0) = rot_euler(1, 0) + tempnode->copy_rot_euler[1];
		rot_euler(2, 0) = rot_euler(2, 0) + tempnode->copy_rot_euler[2];
	}
	//Posição média dos nós da lista - esse será o pólo para o momento a ser transportado no monitor
	coordinates = (1.0 / n_nodes)*coordinates;
	rot_euler = (1.0 / n_nodes)*rot_euler;

	for (int i = 0; i < n_nodes; i++)
	{
		current_node = node_list[i];
		tempnode = db.nodes[current_node - 1];

		temp_coordinates(0, 0) = tempnode->copy_coordinates[0];
		temp_coordinates(1, 0) = tempnode->copy_coordinates[1];
		temp_coordinates(2, 0) = tempnode->copy_coordinates[2];

		(*tempnode->rot_rodrigues)(0, 0) = (tempnode->displacements)[3];
		(*tempnode->rot_rodrigues)(1, 0) = (tempnode->displacements)[4];
		(*tempnode->rot_rodrigues)(2, 0) = (tempnode->displacements)[5];

		//Esforços
		for (int j = 0; j < db.number_GLs_node; j++)
		{
			if (tempnode->GLs[j] < 0 && tempnode->active_GL[j] == 1)	//Se o grau de liberdade for fixo e ativo
				(*tempnode->load)(j, 0) = db.global_P_B(-tempnode->GLs[j] - 1, 0);
			if (tempnode->GLs[j] > 0 && tempnode->active_GL[j] == 1)	//Se o grau de liberdade for livre e ativo
				(*tempnode->load)(j, 0) = db.global_I_A(+tempnode->GLs[j] - 1, 0);
		}
		//Conversão do pseudo-momento para momento (parâmetros de rotação de Rodrigues)
		for (int j = 0; j < 3; j++)
			(*tempnode->moment)(j, 0) = (*tempnode->load)(j + 3, 0);
		//Calculando o operador Xi
		*tempnode->A = skew(*tempnode->rot_rodrigues);			//Matriz A
		tempnode->g = 4.0 / (4.0 + norm(*tempnode->rot_rodrigues)*norm(*tempnode->rot_rodrigues));		//função g(alpha) - em algumas ref. tb. chamado de h(alpha)
		*tempnode->Xi = tempnode->g*(*tempnode->I + 0.5*(*tempnode->A));
		*tempnode->Xi_T_inv = invert3x3(transp(*tempnode->Xi));
		*tempnode->moment = (*tempnode->Xi_T_inv)*(*tempnode->moment);
		for (int j = 0; j < 3; j++)
			(*tempnode->load)(j + 3, 0) = (*tempnode->moment)(j, 0);

		//Atualização com contribuição da força do nó atual
		for (int j = 0; j < 3; j++)
		{
			force(j, 0) = force(j, 0) + (*tempnode->load)(j, 0);		//incremento na força total
			temp_force(j, 0) = (*tempnode->load)(j, 0);					//salvando a força só do nó atual (para transporte do momento)
		}
			
		//Atualização com contribuição do momento do nó atual
		//Transporte
		temp_vec = cross(temp_coordinates - coordinates, temp_force);	//binário de transporte
		for (int j = 0; j < 3; j++)
			moment(j, 0) = moment(j, 0) + (*tempnode->load)(j + 3, 0) + temp_vec(j,0);

	}
	
	//Cabeçalho
	if (first_record == true)
		fprintf(f, "TIME\tX\tY\tZ\tROTX\tROTY\tROTZ\tFX\tFY\tFZ\tMX\tMY\tMZ\n");
	//Informações a serem salvas
	fprintf(f, "%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\n",
		time,
		coordinates(0, 0),
		coordinates(1, 0),
		coordinates(2, 0),
		rot_euler(0, 0),
		rot_euler(1, 0),
		rot_euler(2, 0),
		force(0, 0),
		force(1, 0),
		force(2, 0),
		moment(0, 0),
		moment(1, 0),
		moment(2, 0)
	);
}
