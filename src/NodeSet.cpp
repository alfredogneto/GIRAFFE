#include "NodeSet.h"

#include "Node.h"
#include "Matrix.h"

#include"Database.h"
//Global variables
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
	//Checks the keyword "Set"
	if (!strcmp(s, "NodeSet"))
	{
		fscanf(f, "%s", s);
		number = atoi(s);
	}
	else
		return false;
	fscanf(f, "%s", s);
	//Checks the keyword "Nodes"
	if (!strcmp(s, "Nodes"))
	{
		fscanf(f, "%s", s);
		n_nodes = atoi(s);
		//Allocation of the vector of nodes
		node_list = new int[n_nodes];
	}
	else
		return false;
	//Two possibilities of reading:
	//1 - List
	//2 - Sequence
	fscanf(f, "%s", s);
	if (!strcmp(s, "List"))
	{
		list = true;
		for (int i = 0; i < n_nodes; i++)
		{
			fscanf(f, "%s", s);//Reading of the number of the node
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
			//Geração of the list of nodes
			for (int i = 0; i < n_nodes; i++)
			{
				node_list[i] = initial+i*increment;
			}
		}
		else
			return false;
	}
	//If atingiu this point, sign of reading correct of tudo: returns true
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
	//Variables of the NodeSet for Monitor
	Matrix coordinates(3);
	Matrix rot_euler(3);
	Matrix force(3);
	Matrix moment(3);
	Matrix temp_vec(3);
	Matrix temp_coordinates(3);
	Matrix temp_force(3);
	//Iterates over all the nodes of inside of the "set", saving some results
	//Calculation of the angle of rotation
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
	//Position media of the nodes of the list - this will be the pole for the moment the be transportado in the monitor
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

		//Internal forces
		for (int j = 0; j < db.number_GLs_node; j++)
		{
			if (tempnode->GLs[j] < 0 && tempnode->active_GL[j] == 1)	//If the degree of freedom for fixed and active
				(*tempnode->load)(j, 0) = db.global_P_B(-tempnode->GLs[j] - 1, 0);
			if (tempnode->GLs[j] > 0 && tempnode->active_GL[j] == 1)	//If the degree of freedom for free and active
				(*tempnode->load)(j, 0) = db.global_I_A(+tempnode->GLs[j] - 1, 0);
		}
		//Conversion of the pseudo-moment for moment (parameters of rotation of Rodrigues)
		for (int j = 0; j < 3; j++)
			(*tempnode->moment)(j, 0) = (*tempnode->load)(j + 3, 0);
		//Calculating the operator Xi
		*tempnode->A = skew(*tempnode->rot_rodrigues);			//Matrix A
		tempnode->g = 4.0 / (4.0 + norm(*tempnode->rot_rodrigues)*norm(*tempnode->rot_rodrigues));		//function g(alpha) - in some ref. tb. called of h(alpha)
		*tempnode->Xi = tempnode->g*(*tempnode->I + 0.5*(*tempnode->A));
		*tempnode->Xi_T_inv = invert3x3(transp(*tempnode->Xi));
		*tempnode->moment = (*tempnode->Xi_T_inv)*(*tempnode->moment);
		for (int j = 0; j < 3; j++)
			(*tempnode->load)(j + 3, 0) = (*tempnode->moment)(j, 0);

		//Update with contribution of the force of the node current
		for (int j = 0; j < 3; j++)
		{
			force(j, 0) = force(j, 0) + (*tempnode->load)(j, 0);		//increment in the force total
			temp_force(j, 0) = (*tempnode->load)(j, 0);					//saving the force só of the node current (for transporte of the moment)
		}
			
		//Update with contribution of the moment of the node current
		//Transporte
		temp_vec = cross(temp_coordinates - coordinates, temp_force);	//binary of transporte
		for (int j = 0; j < 3; j++)
			moment(j, 0) = moment(j, 0) + (*tempnode->load)(j + 3, 0) + temp_vec(j,0);

	}
	
	//Header
	if (first_record == true)
		fprintf(f, "TIME\tX\tY\tZ\tROTX\tROTY\tROTZ\tFX\tFY\tFZ\tMX\tMY\tMZ\n");
	//Information the to be saved
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
