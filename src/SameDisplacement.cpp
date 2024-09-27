#include "SameDisplacement.h"

#include "Node.h"
#include "InitialCondition.h"
#include "Database.h"
//Variaveis globais
extern
Database db;

SameDisplacement::SameDisplacement()
{
	
	n_GL = 3;						//Três graus de liberdade (esse vinculo possui 3 multiplicadores de lagrange)
	active_lambda = new int[n_GL];
	lambda = new double[n_GL];
	copy_lambda = new double[n_GL];
	GLs = new int[n_GL];
	node_A = 0;
	node_B = 0;

	for (int i = 0; i < n_GL;i++)
	{
		active_lambda[i] = 0;
		GLs[i] = 0;
		//Chute inicial para os lambdas: valores nulos
		lambda[i] = 0.0;
		copy_lambda[i] = 0.0;
	}
	
	
	
	I3 = Matrix(3, 3);
	I3(0, 0) = 1.0;
	I3(1, 1) = 1.0;
	I3(2, 2) = 1.0;
	r = Matrix(3);
	//Matriz de rigidez tangente e vetor residuo
	stiffness = new Matrix(9,9);
	residual = new Matrix(9,1);

}

SameDisplacement::~SameDisplacement()
{
	delete []active_lambda;
	delete []GLs;
	delete []lambda;
	delete []copy_lambda;
	delete stiffness;
	delete residual;
}

//Leitura
bool SameDisplacement::Read(FILE *f)
{
	char s[1000];
	fscanf(f, "%s", s);
	number = atoi(s);

	fscanf(f, "%s", s);
	if (!strcmp(s, "Nodes"))
	{
		fscanf(f, "%s", s);
		node_A = atoi(s);

		fscanf(f, "%s", s);
		node_B = atoi(s);
	}
	else
		return false;
	//Salva a posição (stream)
	fpos_t pos;
	fgetpos(f, &pos);
	fscanf(f, "%s", s);
	if (!strcmp(s, "BoolTable"))
		bool_table.Read(f);
	else
	{
		fsetpos(f, &pos);
		bool_table.SetDefault(true);
	}

	return true;
}

//Gravação
void SameDisplacement::Write(FILE *f)
{
	fprintf(f, "SameDisplacement\t%d\tNodes\t%d\t%d\n",
		number,
		node_A,
		node_B);
}

//Escreve no monitor do SpecialConstraint
void SameDisplacement::WriteMonitor(FILE *f, bool first_record, double time)
{

}

void SameDisplacement::WriteVTK_XMLRender(FILE *f)
{

}

//Checa inconsistências no SC para evitar erros de execução
bool SameDisplacement::Check()
{
	if (node_A > db.number_nodes)
		return false;
	if (node_B > db.number_nodes)
		return false;
	//Checagem das condições iniciais
	int temp_node = 0;
	for (int i = 0; i < db.number_IC; i++)
	{
		temp_node = db.IC[i]->node;
		if (node_B == temp_node)
		{
			db.myprintf("Warning in Special Constraint %d.\nInitial Condition %d was prescribed to node %d (slave), leading to ignoring some of its components.\n", number, db.IC[i]->number, db.IC[i]->node);
		}
	}
	return true;
}

//Montagem dos residuos e rigidez tangente
void SameDisplacement::Mount()
{
	//Nesse vinculo a rigidez tangente não se modifica nunca. e montada no PreCalc()
	//Se o vinculo estiver ativo, realiza a montagem
	if (active_lambda[0] == 1 && active_lambda[1] == 1 && active_lambda[2] == 1)
	{
		for (int i = 0; i < 3; i++)
			r(i, 0) = db.nodes[node_A - 1]->displacements[i] - db.nodes[node_B - 1]->displacements[i];

		//Atualização do vetor de residuos
		(*residual)(0, 0) = lambda[0];
		(*residual)(1, 0) = lambda[1];
		(*residual)(2, 0) = lambda[2];

		(*residual)(3, 0) = -lambda[0];
		(*residual)(4, 0) = -lambda[1];
		(*residual)(5, 0) = -lambda[2];

		(*residual)(6, 0) = r(0, 0);
		(*residual)(7, 0) = r(1, 0);
		(*residual)(8, 0) = r(2, 0);
	}
}

//Preenche a contribuição do elemento nas matrizes globais
void SameDisplacement::MountGlobal()
{
	//Variaveis temporarias para salvar a indexação global dos graus de liberdade a serem setados na matriz de rigidez global
	int GL_global_1 = 0;
	int GL_global_2 = 0;
	double anterior = 0;
	//Se o vinculo estiver ativo, realiza a montagem
	if (active_lambda[0] == 1 && active_lambda[1] == 1 && active_lambda[2] == 1)
	{
		for (int i = 0; i < 9; i++)
		{
			//////////////MONTAGEM DO VETOR DE ESFORÇOS DESBALANCEADOS//////////////////
			//Toma do vetor de GL globais, a indexação de cada grau de liberdade global
			if (i < 3)//uA
				GL_global_1 = db.nodes[node_A - 1]->GLs[i];
			else
			{
				if (i<6)//uB
					GL_global_1 = db.nodes[node_B - 1]->GLs[i - 3];
				else
				{
					//lambda
					GL_global_1 = GLs[i - 6];
				}
			}

			//Caso o grau de liberdade seja livre:
			if (GL_global_1 > 0)
			{
				anterior = db.global_P_A(GL_global_1 - 1, 0);
				db.global_P_A(GL_global_1 - 1, 0) = anterior + (*residual)(i, 0);
				anterior = db.global_I_A(GL_global_1 - 1, 0);
				db.global_I_A(GL_global_1 - 1, 0) = anterior + (*residual)(i, 0);
			}
			else
			{
				anterior = db.global_P_B(-GL_global_1 - 1, 0);
				db.global_P_B(-GL_global_1 - 1, 0) = anterior + (*residual)(i, 0);
			}
			for (int j = 0; j < 9; j++)
			{
				//////////////////////MONTAGEM DA MATRIZ DE RIGIDEZ/////////////////////////
				//Toma do vetor de GL globais, a indexação de cada grau de liberdade global
				if (j < 3)//uA
					GL_global_2 = db.nodes[node_A - 1]->GLs[j];
				else
				{
					if (j<6)//uB
						GL_global_2 = db.nodes[node_B - 1]->GLs[j - 3];
					else
					{
						//lambda
						GL_global_2 = GLs[j - 6];
					}
				}

				//Caso os graus de liberdade sejam ambos livres (Matriz Kaa)
				if (GL_global_1 > 0 && GL_global_2 > 0)
					db.global_stiffness_AA.setValue(GL_global_1 - 1, GL_global_2 - 1, (*stiffness)(i, j));
				//Caso os graus de liberdade sejam ambos fixos (Matriz Kbb)
				if (GL_global_1 < 0 && GL_global_2 < 0)
					db.global_stiffness_BB.setValue(-GL_global_1 - 1, -GL_global_2 - 1, (*stiffness)(i, j));
				//Caso os graus de liberdade sejam livre e fixo (Matriz Kab)
				if (GL_global_1 > 0 && GL_global_2 < 0)
					db.global_stiffness_AB.setValue(GL_global_1 - 1, -GL_global_2 - 1, (*stiffness)(i, j));
				//Caso os graus de liberdade sejam fixo e livre (Matriz Kba)
				if (GL_global_1 < 0 && GL_global_2 > 0)
					db.global_stiffness_BA.setValue(-GL_global_1 - 1, GL_global_2 - 1, (*stiffness)(i, j));
			}
		}
	}
}

void SameDisplacement::ComputeInitialGuessDisplacements()
{
	//Se for o step de criação do vinculo inicializa condições iniciais
	if (bool_table.GetAt(db.current_solution_number - 1) == true)
	{
		//Percorre GL e seta condições iniciais
		for (int i = 0; i < 3; i++)
		{
			if (active_lambda[i] == 1)
			{
				db.nodes[node_B - 1]->displacements[i] = db.nodes[node_A - 1]->displacements[i];
			}

		}
	}
}

//Computa efeito das condições iniciais nos nós da restrição
void SameDisplacement::ComputeVelAccel()
{
	//Se for o step de criação do vinculo inicializa condições iniciais
	if (bool_table.GetAt(db.current_solution_number - 1))
	{
		//Percorre GL e seta condições iniciais
		for (int i = 0; i < 3; i++)
		{
			if (active_lambda[i] == 1)
			{
				db.nodes[node_B - 1]->vel[i] = db.nodes[node_A - 1]->vel[i];
				db.nodes[node_B - 1]->accel[i] = db.nodes[node_A - 1]->accel[i];
			}

		}
	}

	
	
}

//Pre-calculo de variaveis que e feito uma unica vez no inicio
void SameDisplacement::PreCalc()
{
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			(*stiffness)(i, j + 6) = I3(i, j);
			(*stiffness)(i + 3, j + 6) = -I3(i, j);
			(*stiffness)(i + 6, j) = I3(i, j);
			(*stiffness)(i + 6, j + 3) = -I3(i, j);
		}
	}
}

//Salvando variaveis da configuração convergida
void SameDisplacement::SaveLagrange()
{
	//PrintPtr(lambda, 3);
	for (int i = 0; i < n_GL; i++)
		copy_lambda[i] = lambda[i];
}

//Checa quais multiplicadores de lagrange serão ativados,de acordo com a ativação dos GLs dos nós dos quais a special constraint participa
void SameDisplacement::ActivateDOFs()
{
	//Ativa GLs de translação dos nós A e B
	for (int i = 0; i < 3; i++)
		db.nodes[node_A - 1]->active_GL[i] = 1;
	for (int i = 0; i < 3; i++)
		db.nodes[node_B - 1]->active_GL[i] = 1;

	if (bool_table.GetAt(db.current_solution_number - 1))
	{
		
		//Ativa GLs de multiplicadores de Lagrange
		active_lambda[0] = 1;
		active_lambda[1] = 1;
		active_lambda[2] = 1;
	}
	else
	{
		active_lambda[0] = 0;
		active_lambda[1] = 0;
		active_lambda[2] = 0;
	}
}
