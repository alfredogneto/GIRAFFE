#include "UniversalJoint.h"

#include "Node.h"
#include "InitialCondition.h"
#include "CoordinateSystem.h"
#include"Database.h"
//Variaveis globais
extern
Database db;

UniversalJoint::UniversalJoint()
{
	n_GL = 4;						//Dois graus de liberdade (esse vinculo possui 4 multiplicadores de lagrange)
	active_lambda = new int[n_GL];
	lambda = new double[n_GL];
	copy_lambda = new double[n_GL];
	GLs = new int[n_GL];
	node_A = 0;
	node_B = 0;
	csA = 0;
	csB = 0;

	for (int i = 0; i < n_GL; i++)
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
	r1 = Matrix(3);

	alphaA = Matrix(3);
	alphaB = Matrix(3);
	alphaiA = Matrix(3);
	alphaiB = Matrix(3);
	A = Matrix(3, 3);
	QA = Matrix(3, 3);
	QB = Matrix(3, 3);
	ei1A = Matrix(3);
	ei2B = Matrix(3);
	
	//As contribuições estão divididas em duas partes:
	//Parte 1 - exatamente a mesma contribuição do Same Displacement
	//Parte 2 - restrições das rotações
	//Matriz de rigidez tangente e vetor residuo
	stiffness1 = new Matrix(9, 9);
	residual1 = new Matrix(9, 1);
	//Matriz de rigidez tangente e vetor residuo
	stiffness2 = new double*[7];
	for (int i = 0; i < 7; i++)
		stiffness2[i] = new double[7];
	residual2 = new double[7];
	temp_lambda = new double[1];
}


UniversalJoint::~UniversalJoint()
{
	delete[]active_lambda;
	delete[]GLs;
	delete[]lambda;
	delete[]copy_lambda;
	delete stiffness1;
	delete residual1;
	for (int i = 0; i < 7; i++)
		delete[]stiffness2[i];
	delete[]stiffness2;
	delete[]residual2;
	delete[]temp_lambda;
}

//Leitura
bool UniversalJoint::Read(FILE *f)
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

	fscanf(f, "%s", s);
	if (!strcmp(s, "CSA"))
	{
		fscanf(f, "%s", s);
		csA = atoi(s);
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "CSB"))
	{
		fscanf(f, "%s", s);
		csB = atoi(s);
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

	//TESTE para verificar ortogonalidade entre E1A e E2B
	if (dot(*db.CS[csA - 1]->E1, *db.CS[csB - 1]->E2) != 0.0)
	{
		db.myprintf("Error in UniversalJoint %d. E1 from CSA must be orthogonal to E2 from CSB.\n",number);
		return false;
	}

	return true;
}

//Gravação
void UniversalJoint::Write(FILE *f)
{
	fprintf(f, "UniversalJoint\t%d\tNodes\t%d\t%d\tCSA\t%d\tCSB\t%d\n",
		number,
		node_A,
		node_B,
		csA,
		csB);
}

//Escreve no monitor do SpecialConstraint
void UniversalJoint::WriteMonitor(FILE *f, bool first_record, double time)
{

}

void UniversalJoint::WriteVTK_XMLRender(FILE *f)
{

}

//Checa inconsistências no SC para evitar erros de execução
bool UniversalJoint::Check()
{	
	if (node_A > db.number_nodes)
		return false;
	if (node_B > db.number_nodes)
		return false;
	if (csA > db.number_CS)
		return false;
	if (csB > db.number_CS)
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
void UniversalJoint::Mount()
{
	//Nesse vinculo a rigidez tangente da parte 1 não se modifica nunca. e montada no PreCalc()
	//Montagem do residuo da parte 1:
	if (active_lambda[0] == 1 && active_lambda[1] == 1 && active_lambda[2] == 1)
	{
		for (int i = 0; i < 3; i++)
			r1(i, 0) = (db.nodes[node_A - 1]->copy_coordinates[i] - db.nodes[node_A - 1]->ref_coordinates[i] + db.nodes[node_A - 1]->displacements[i]) -
			(db.nodes[node_B - 1]->copy_coordinates[i] - db.nodes[node_B - 1]->ref_coordinates[i] + db.nodes[node_B - 1]->displacements[i]);

		//Atualização do vetor de residuos
		(*residual1)(0, 0) = lambda[0];
		(*residual1)(1, 0) = lambda[1];
		(*residual1)(2, 0) = lambda[2];

		(*residual1)(3, 0) = -lambda[0];
		(*residual1)(4, 0) = -lambda[1];
		(*residual1)(5, 0) = -lambda[2];

		(*residual1)(6, 0) = r1(0, 0);
		(*residual1)(7, 0) = r1(1, 0);
		(*residual1)(8, 0) = r1(2, 0);
	}

	//Montagem da rigidez tangente e residuo da parte 2:
	if (active_lambda[3] == 1)
	{
		for (int i = 0; i < 3; i++)
		{
			alphaA(i, 0) = db.nodes[node_A - 1]->displacements[i + 3];	//vetor rotação (atual) do nó A
			alphaB(i, 0) = db.nodes[node_B - 1]->displacements[i + 3];	//vetor rotação (atual) do nó B
			alphaiA(i, 0) = db.nodes[node_A - 1]->copy_coordinates[i + 3];	//vetor rotação acumulada (do inicio) do nó A
			alphaiB(i, 0) = db.nodes[node_B - 1]->copy_coordinates[i + 3];	//vetor rotação acumulada (do inicio) do nó B
		}
		
		alpha_escalar_i = norm(alphaiA);
		A = skew(alphaiA);
		g = 4.0 / (4.0 + alpha_escalar_i*alpha_escalar_i);
		QA = I3 + g*(A + 0.5*(A*A));
		ei1A = QA*(*db.CS[csA - 1]->E1);	//Eixo e1 no inicio do incremento
		
		alpha_escalar_i = norm(alphaiB);
		A = skew(alphaiB);
		g = 4.0 / (4.0 + alpha_escalar_i*alpha_escalar_i);
		QB = I3 + g*(A + 0.5*(A*A));
		ei2B = QB*(*db.CS[csB - 1]->E2);	//Eixo e2 no inicio do incremento
		
		temp_lambda[0] = lambda[3];
		EvaluateUniversalJointContribution(temp_v, residual2, stiffness2, alphaA.getMatrix(), alphaB.getMatrix(), ei1A.getMatrix(), ei2B.getMatrix(), temp_lambda);
		//PrintPtr(residual2, 7);
		//printf("Error in residual %.6e\n", residual2[6]);
		//ei1A.print();
		//ei2B.print();
		//PrintPtr(stiffness2, 7, 7);
	}
}

//Preenche a contribuição do elemento nas matrizes globais
void UniversalJoint::MountGlobal()
{
	//Variaveis temporarias para salvar a indexação global dos graus de liberdade a serem setados na matriz de rigidez global
	int GL_global_1 = 0;
	int GL_global_2 = 0;
	double anterior = 0;

	//PARTE 1 - deslocamentos
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
				db.global_P_A(GL_global_1 - 1, 0) = anterior + (*residual1)(i, 0);
				anterior = db.global_I_A(GL_global_1 - 1, 0);
				db.global_I_A(GL_global_1 - 1, 0) = anterior + (*residual1)(i, 0);
			}
			else
			{
				anterior = db.global_P_B(-GL_global_1 - 1, 0);
				db.global_P_B(-GL_global_1 - 1, 0) = anterior + (*residual1)(i, 0);
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
					db.global_stiffness_AA.setValue(GL_global_1 - 1, GL_global_2 - 1, (*stiffness1)(i, j));
				//Caso os graus de liberdade sejam ambos fixos (Matriz Kbb)
				if (GL_global_1 < 0 && GL_global_2 < 0)
					db.global_stiffness_BB.setValue(-GL_global_1 - 1, -GL_global_2 - 1, (*stiffness1)(i, j));
				//Caso os graus de liberdade sejam livre e fixo (Matriz Kab)
				if (GL_global_1 > 0 && GL_global_2 < 0)
					db.global_stiffness_AB.setValue(GL_global_1 - 1, -GL_global_2 - 1, (*stiffness1)(i, j));
				//Caso os graus de liberdade sejam fixo e livre (Matriz Kba)
				if (GL_global_1 < 0 && GL_global_2 > 0)
					db.global_stiffness_BA.setValue(-GL_global_1 - 1, GL_global_2 - 1, (*stiffness1)(i, j));
			}
		}
	}

	//PARTE 2 - rotações
	if (active_lambda[3] == 1)
	{
		for (int i = 0; i < 7; i++)
		{
			//////////////MONTAGEM DO VETOR DE ESFORÇOS DESBALANCEADOS//////////////////
			//Toma do vetor de GL globais, a indexação de cada grau de liberdade global
			if (i < 3)//alpha_A -> i=0,1,2
				GL_global_1 = db.nodes[node_A - 1]->GLs[i + 3];
			else
			{
				if (i<6)//alpha_B -> i=3,4,5
					GL_global_1 = db.nodes[node_B - 1]->GLs[i];
				else
				{
					//lambda  -> i=6
					GL_global_1 = GLs[i - 3];
				}
			}

			//Caso o grau de liberdade seja livre:
			if (GL_global_1 > 0)
			{
				anterior = db.global_P_A(GL_global_1 - 1, 0);
				db.global_P_A(GL_global_1 - 1, 0) = anterior + residual2[i];
				anterior = db.global_I_A(GL_global_1 - 1, 0);
				db.global_I_A(GL_global_1 - 1, 0) = anterior + residual2[i];
			}
			else
			{
				anterior = db.global_P_B(-GL_global_1 - 1, 0);
				db.global_P_B(-GL_global_1 - 1, 0) = anterior + residual2[i];
			}
			for (int j = 0; j < 7; j++)
			{
				//////////////////////MONTAGEM DA MATRIZ DE RIGIDEZ/////////////////////////
				//Toma do vetor de GL globais, a indexação de cada grau de liberdade global
				if (j < 3)//alpha_A -> j=0,1,2
					GL_global_2 = db.nodes[node_A - 1]->GLs[j + 3];
				else
				{
					if (j<6)//alpha_B -> j=3,4,5
						GL_global_2 = db.nodes[node_B - 1]->GLs[j];
					else
					{
						//lambda  -> i=6
						GL_global_2 = GLs[j - 3];
					}
				}

				//Caso os graus de liberdade sejam ambos livres (Matriz Kaa)
				if (GL_global_1 > 0 && GL_global_2 > 0)
					db.global_stiffness_AA.setValue(GL_global_1 - 1, GL_global_2 - 1, stiffness2[i][j]);
				//Caso os graus de liberdade sejam ambos fixos (Matriz Kbb)
				if (GL_global_1 < 0 && GL_global_2 < 0)
					db.global_stiffness_BB.setValue(-GL_global_1 - 1, -GL_global_2 - 1, stiffness2[i][j]);
				//Caso os graus de liberdade sejam livre e fixo (Matriz Kab)
				if (GL_global_1 > 0 && GL_global_2 < 0)
					db.global_stiffness_AB.setValue(GL_global_1 - 1, -GL_global_2 - 1, stiffness2[i][j]);
				//Caso os graus de liberdade sejam fixo e livre (Matriz Kba)
				if (GL_global_1 < 0 && GL_global_2 > 0)
					db.global_stiffness_BA.setValue(-GL_global_1 - 1, GL_global_2 - 1, stiffness2[i][j]);
			}
		}
	}
}

void UniversalJoint::ComputeInitialGuessDisplacements()
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

		//Projeção da velocidade angular inicial no sistema local - para impor somente velocidade nos GL restritos - impõe mesma velocidade angular na direção do eixo da articulação
		Matrix omega_global_A(3);
		omega_global_A(0, 0) = db.nodes[node_A - 1]->displacements[3];
		omega_global_A(1, 0) = db.nodes[node_A - 1]->displacements[4];
		omega_global_A(2, 0) = db.nodes[node_A - 1]->displacements[5];
		Matrix omega_global_B(3);
		omega_global_B(0, 0) = db.nodes[node_B - 1]->displacements[3];
		omega_global_B(1, 0) = db.nodes[node_B - 1]->displacements[4];
		omega_global_B(2, 0) = db.nodes[node_B - 1]->displacements[5];
		Matrix omega_local_A = (*db.CS[csA - 1]->Q)*omega_global_A;
		Matrix omega_local_B = (*db.CS[csA - 1]->Q)*omega_global_B; //Obs: csA e csB possuem mesma direção E3!
		//Imposições de mesma condição inicial no sistema local - somente direções 1 e 2 locais
		omega_local_B(2, 0) = omega_local_A(2, 0);
		//De volta para o sistema global - nó B
		omega_global_B = transp(*db.CS[csA - 1]->Q)*omega_local_B;
		if (active_lambda[3] == 1 && active_lambda[4] == 1)
		{
			db.nodes[node_B - 1]->displacements[3] = omega_global_B(0, 0);
			db.nodes[node_B - 1]->displacements[4] = omega_global_B(1, 0);
			db.nodes[node_B - 1]->displacements[5] = omega_global_B(2, 0);
		}
	}
}

//Computa efeito das condições iniciais nos nós da restrição
void UniversalJoint::ComputeVelAccel()
{
	//Se for o step de criação do vinculo inicializa condições iniciais
	if (bool_table.GetAt(db.current_solution_number - 1) == true)
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

		//Projeção da velocidade angular inicial no sistema local - para impor somente velocidade nos GL restritos - impõe mesma velocidade angular na direção do eixo da articulação
		Matrix omega_global_A(3);
		omega_global_A(0, 0) = db.nodes[node_A - 1]->vel[3];
		omega_global_A(1, 0) = db.nodes[node_A - 1]->vel[4];
		omega_global_A(2, 0) = db.nodes[node_A - 1]->vel[5];
		Matrix omega_global_B(3);
		omega_global_B(0, 0) = db.nodes[node_B - 1]->vel[3];
		omega_global_B(1, 0) = db.nodes[node_B - 1]->vel[4];
		omega_global_B(2, 0) = db.nodes[node_B - 1]->vel[5];
		Matrix omega_local_A = (*db.CS[csA - 1]->Q)*omega_global_A;
		Matrix omega_local_B = (*db.CS[csA - 1]->Q)*omega_global_B; //Obs: csA e csB possuem mesma direção E3!
		//Imposições de mesma condição inicial no sistema local - somente direções 1 e 2 locais
		omega_local_B(2, 0) = omega_local_A(2, 0);
		//De volta para o sistema global - nó B
		omega_global_B = transp(*db.CS[csA - 1]->Q)*omega_local_B;
		if (active_lambda[3] == 1 && active_lambda[4] == 1)
		{
			db.nodes[node_B - 1]->vel[3] = omega_global_B(0, 0);
			db.nodes[node_B - 1]->vel[4] = omega_global_B(1, 0);
			db.nodes[node_B - 1]->vel[5] = omega_global_B(2, 0);
		}
	}

	
}

//Pre-calculo de variaveis que e feito uma unica vez no inicio
void UniversalJoint::PreCalc()
{
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			(*stiffness1)(i, j + 6) = I3(i, j);
			(*stiffness1)(i + 3, j + 6) = -I3(i, j);
			(*stiffness1)(i + 6, j) = I3(i, j);
			(*stiffness1)(i + 6, j + 3) = -I3(i, j);
		}
	}

	//Zerando contribuição dos residuos e rigidez 2
	for (int i = 0; i < 7; i++)
	{
		residual2[i] = 0.0;
		for (int j = 0; j < 7; j++)
		{
			stiffness2[i][j] = 0.0;
		}
	}
}

//Salvando variaveis da configuração convergida
void UniversalJoint::SaveLagrange()
{
	for (int i = 0; i < n_GL; i++)
		copy_lambda[i] = lambda[i];
}

//Checa quais multiplicadores de lagrange serão ativados,de acordo com a ativação dos GLs dos nós dos quais a special constraint participa
void UniversalJoint::ActivateDOFs()
{
	//Ativa GLs de translação e rotação dos nós A e B
	for (int i = 0; i < 6; i++)
		db.nodes[node_A - 1]->active_GL[i] = 1;
	for (int i = 0; i < 6; i++)
		db.nodes[node_B - 1]->active_GL[i] = 1;

	if (bool_table.GetAt(db.current_solution_number - 1))
	{
		//Ativa GLs de multiplicadores de Lagrange
		active_lambda[0] = 1;
		active_lambda[1] = 1;
		active_lambda[2] = 1;
		active_lambda[3] = 1;
	}
	else
	{
		active_lambda[0] = 0;
		active_lambda[1] = 0;
		active_lambda[2] = 0;
		active_lambda[3] = 0;
	}
}

//Calcula contribuições do residuo e operador tangente - gerado no AceGen
void UniversalJoint::EvaluateUniversalJointContribution(double *v, double *residual
	, double **stiffness, double *alphaA, double *alphaB, double *ei1A
	, double *ei2B, double *lambda)
{
	int i01; int i02;
	v[523] = 0.5e0*ei1A[2];
	v[522] = 0.5e0*ei1A[1];
	v[521] = 0.5e0*ei1A[0];
	v[519] = Power(alphaB[2], 2);
	v[518] = 0.5e0*alphaB[2];
	v[517] = 2e0*alphaB[2];
	v[516] = Power(alphaB[1], 2);
	v[515] = 0.5e0*alphaB[1];
	v[514] = 2e0*alphaB[1];
	v[513] = Power(alphaB[0], 2);
	v[512] = 2e0*alphaB[0];
	v[511] = 0.5e0*alphaB[0];
	v[509] = Power(alphaA[2], 2);
	v[508] = 0.5e0*alphaA[2];
	v[507] = 2e0*alphaA[2];
	v[506] = Power(alphaA[1], 2);
	v[505] = 0.5e0*alphaA[1];
	v[504] = 2e0*alphaA[1];
	v[503] = Power(alphaA[0], 2);
	v[502] = 2e0*alphaA[0];
	v[501] = 0.5e0*alphaA[0];
	v[85] = alphaA[1] * v[501];
	v[142] = -v[503] - v[506];
	v[149] = -alphaA[2] + v[85];
	v[145] = alphaA[2] + v[85];
	v[92] = alphaA[2] * v[505];
	v[147] = -alphaA[0] + v[92];
	v[548] = ei1A[0] * v[145] + ei1A[2] * v[147];
	v[144] = alphaA[0] + v[92];
	v[90] = alphaA[2] * v[501];
	v[150] = alphaA[1] + v[90];
	v[544] = ei1A[1] * v[149] + ei1A[2] * v[150];
	v[143] = -alphaA[1] + v[90];
	v[551] = ei1A[0] * v[143] + ei1A[1] * v[144];
	v[156] = 4e0 + v[503] + v[506] + v[509];
	v[194] = 1e0 / Power(v[156], 2);
	v[510] = -4e0*v[194];
	v[197] = v[507] * v[510];
	v[196] = v[504] * v[510];
	v[195] = v[502] * v[510];
	v[148] = -v[506] - v[509];
	v[146] = -v[503] - v[509];
	v[110] = alphaB[1] * v[511];
	v[134] = -v[513] - v[516];
	v[538] = 0.5e0*v[134];
	v[131] = -alphaB[2] + v[110];
	v[128] = alphaB[2] + v[110];
	v[117] = alphaB[2] * v[515];
	v[129] = -alphaB[0] + v[117];
	v[532] = ei2B[0] * v[128] + ei2B[2] * v[129];
	v[126] = alphaB[0] + v[117];
	v[115] = alphaB[2] * v[511];
	v[132] = alphaB[1] + v[115];
	v[527] = ei2B[1] * v[131] + ei2B[2] * v[132];
	v[125] = -alphaB[1] + v[115];
	v[537] = ei2B[0] * v[125] + ei2B[1] * v[126];
	v[170] = 4e0 + v[513] + v[516] + v[519];
	v[198] = 1e0 / Power(v[170], 2);
	v[520] = -4e0*v[198];
	v[558] = lambda[0] * v[520];
	v[201] = v[517] * v[520];
	v[200] = v[514] * v[520];
	v[199] = v[512] * v[520];
	v[554] = lambda[0] * v[199];
	v[136] = -v[513] - v[519];
	v[533] = 0.5e0*v[136];
	v[135] = -v[516] - v[519];
	v[528] = 0.5e0*v[135];
	v[366] = ei1A[0] * v[501];
	v[365] = ei1A[1] * v[505];
	v[370] = v[148] * v[521] + v[544];
	v[369] = v[146] * v[522] + v[548];
	v[368] = v[142] * v[523] + v[551];
	v[363] = ei1A[2] * v[508];
	v[79] = 4e0 / v[156];
	v[547] = 0.5e0*v[79];
	v[441] = -(v[521] * v[79]);
	v[426] = -(v[522] * v[79]);
	v[425] = -(v[523] * v[79]);
	v[104] = 4e0 / v[170];
	v[526] = -0.5e0*v[104];
	v[525] = ei2B[2] * v[104];
	v[271] = ei2B[1] * v[104];
	v[269] = v[517] * v[526];
	v[272] = lambda[0] * (-v[271] + v[511] * v[525] + v[201] * v[527] + ei2B[0] * (v[269] + v[201] * v[528]));
	v[524] = v[272] * v[79];
	v[465] = v[272] * v[441];
	v[411] = ei1A[1] * v[524];
	v[390] = ei1A[2] * v[524];
	v[265] = v[514] * v[526];
	v[268] = lambda[0] * (v[271] * v[511] + v[525] + v[200] * v[527] + ei2B[0] * (v[265] + v[200] * v[528]));
	v[529] = v[268] * v[79];
	v[462] = v[268] * v[441];
	v[408] = ei1A[1] * v[529];
	v[387] = ei1A[2] * v[529];
	v[263] = v[518] * v[525];
	v[262] = v[271] * v[515];
	v[264] = lambda[0] * (v[262] + v[263] + v[199] * (v[527] + ei2B[0] * v[528]));
	v[530] = v[264] * v[79];
	v[459] = v[264] * v[441];
	v[405] = ei1A[1] * v[530];
	v[384] = ei1A[2] * v[530];
	v[258] = -(ei2B[0] * v[104]);
	v[259] = lambda[0] * (-v[258] + v[515] * v[525] + v[201] * v[532] + ei2B[1] * (v[269] + v[201] * v[533]));
	v[531] = v[259] * v[79];
	v[464] = v[259] * v[426];
	v[412] = ei1A[0] * v[531];
	v[413] = v[411] + v[412];
	v[356] = ei1A[2] * v[531];
	v[255] = -(v[258] * v[511]);
	v[256] = lambda[0] * (v[255] + v[263] + v[200] * (v[532] + ei2B[1] * v[533]));
	v[534] = v[256] * v[79];
	v[461] = v[256] * v[426];
	v[409] = ei1A[0] * v[534];
	v[410] = v[408] + v[409];
	v[353] = ei1A[2] * v[534];
	v[251] = v[512] * v[526];
	v[253] = lambda[0] * (-(v[258] * v[515]) - v[525] + v[199] * v[532] + ei2B[1] * (v[251] + v[199] * v[533]));
	v[535] = v[253] * v[79];
	v[458] = v[253] * v[426];
	v[406] = ei1A[0] * v[535];
	v[407] = v[405] + v[406];
	v[350] = ei1A[2] * v[535];
	v[249] = lambda[0] * (v[255] + v[262] + v[201] * (v[537] + ei2B[2] * v[538]));
	v[536] = v[249] * v[79];
	v[444] = v[249] * v[425];
	v[391] = ei1A[0] * v[536];
	v[392] = v[390] + v[391];
	v[445] = (v[249] * v[368] + v[259] * v[369] + v[272] * v[370])*v[510];
	v[572] = v[445] + v[464];
	v[357] = ei1A[1] * v[536];
	v[358] = v[356] + v[357];
	v[247] = lambda[0] * (v[258] + v[271] * v[518] + v[200] * v[537] + ei2B[2] * (v[265] + v[200] * v[538]));
	v[539] = v[247] * v[79];
	v[440] = v[247] * v[425];
	v[388] = ei1A[0] * v[539];
	v[389] = v[387] + v[388];
	v[442] = (v[247] * v[368] + v[256] * v[369] + v[268] * v[370])*v[510];
	v[571] = v[442] + v[461];
	v[354] = ei1A[1] * v[539];
	v[355] = v[353] + v[354];
	v[245] = lambda[0] * (v[271] - v[258] * v[518] + v[199] * v[537] + ei2B[2] * (v[251] + v[199] * v[538]));
	v[540] = v[245] * v[79];
	v[437] = v[245] * v[425];
	v[385] = ei1A[0] * v[540];
	v[386] = v[384] + v[385];
	v[438] = (v[245] * v[368] + v[253] * v[369] + v[264] * v[370])*v[510];
	v[570] = v[438] + v[458];
	v[351] = ei1A[1] * v[540];
	v[352] = v[350] + v[351];
	v[107] = 1e0 - v[135] * v[526];
	v[273] = ei2B[0] * v[107] + v[104] * v[527];
	v[541] = v[273] * v[79];
	v[468] = v[273] * v[441];
	v[414] = ei1A[1] * v[541];
	v[393] = ei1A[2] * v[541];
	v[113] = 1e0 - v[136] * v[526];
	v[260] = ei2B[1] * v[113] + v[104] * v[532];
	v[542] = v[260] * v[79];
	v[467] = v[260] * v[426];
	v[415] = ei1A[0] * v[542];
	v[416] = v[414] + v[415];
	v[359] = ei1A[2] * v[542];
	v[119] = 1e0 - v[134] * v[526];
	v[250] = ei2B[2] * v[119] + v[104] * v[537];
	v[543] = v[250] * v[79];
	v[447] = v[250] * v[425];
	v[394] = ei1A[0] * v[543];
	v[395] = v[393] + v[394];
	v[448] = (v[250] * v[368] + v[260] * v[369] + v[273] * v[370])*v[510];
	v[573] = v[448] + v[467];
	v[360] = ei1A[1] * v[543];
	v[361] = v[359] + v[360];
	v[120] = ei1A[0] * (1e0 + v[148] * v[547]) + v[544] * v[79];
	v[556] = 0.5e0*v[120];
	v[546] = ei2B[2] * v[120];
	v[545] = v[120] * v[201];
	v[480] = -(ei2B[0] * v[556]);
	v[492] = v[201] * v[480];
	v[338] = ei2B[1] * v[545];
	v[300] = ei2B[2] * v[545];
	v[297] = v[200] * v[546];
	v[281] = v[120] * v[511];
	v[179] = v[104] * v[480];
	v[177] = v[120] * v[271];
	v[172] = v[104] * v[546];
	v[121] = ei1A[1] * (1e0 + v[146] * v[547]) + v[548] * v[79];
	v[567] = v[121] * v[129] + v[120] * v[132];
	v[557] = 0.5e0*v[121];
	v[550] = ei2B[2] * v[121];
	v[549] = ei2B[0] * v[121];
	v[576] = v[515] * (ei2B[1] * v[120] + v[549]);
	v[471] = -(ei2B[1] * v[557]);
	v[493] = v[201] * v[471];
	v[337] = v[201] * v[549];
	v[339] = lambda[0] * (v[337] + v[338]);
	v[483] = v[554] * v[576];
	v[319] = v[201] * v[550];
	v[316] = v[200] * v[550];
	v[280] = v[121] * v[515];
	v[180] = v[104] * v[471];
	v[582] = v[179] + v[180];
	v[178] = v[104] * v[549];
	v[583] = -v[177] + v[178];
	v[340] = v[177] + v[178];
	v[168] = v[104] * v[550];
	v[122] = ei1A[2] * (1e0 + v[142] * v[547]) + v[551] * v[79];
	v[566] = v[122] * v[126] + v[120] * v[131];
	v[565] = v[122] * v[125] + v[121] * v[128];
	v[555] = 0.5e0*v[122];
	v[553] = ei2B[0] * v[122];
	v[552] = ei2B[1] * v[122];
	v[472] = -(ei2B[2] * v[555]);
	v[486] = v[201] * v[472];
	v[481] = v[200] * v[472];
	v[318] = v[201] * v[552];
	v[315] = v[200] * v[552];
	v[495] = lambda[0] * (v[315] + v[316])*v[518];
	v[299] = v[201] * v[553];
	v[296] = v[200] * v[553];
	v[494] = v[518] * (v[546] + v[553])*v[554];
	v[283] = ei2B[0] * (v[135] * v[556] + v[565]) + ei2B[1] * (v[136] * v[557] + v[566]) + ei2B[2] * (v[134] * v[555] + v[567]
		);
	v[489] = v[283] * v[520];
	v[487] = (ei2B[2] * (v[280] + v[281]) + ei2B[0] * (v[121] + v[122] * v[511] - v[517] * v[556]) + ei2B[1] * (-v[120]
		+ v[122] * v[515] - v[517] * v[557]))*v[558];
	v[278] = v[122] * v[518];
	v[482] = (ei2B[1] * (v[278] + v[281]) + ei2B[2] * (v[120] + v[121] * v[518] - v[514] * v[555]) + ei2B[0] * (-v[122]
		+ v[121] * v[511] - v[514] * v[556]))*v[558];
	v[174] = v[104] * v[472];
	v[580] = v[174] + v[179];
	v[575] = v[174] + v[180];
	v[173] = v[104] * v[553];
	v[581] = v[172] - v[173];
	v[302] = v[172] + v[173];
	v[169] = v[104] * v[552];
	v[579] = -v[168] + v[169];
	v[321] = v[168] + v[169];
	v[127] = lambda[0] * v[250];
	v[560] = ei1A[1] * v[127];
	v[559] = ei1A[0] * v[127];
	v[418] = -(v[127] * v[523]);
	v[382] = v[197] * v[559];
	v[379] = v[196] * v[559];
	v[348] = v[197] * v[560];
	v[345] = v[196] * v[560];
	v[160] = v[418] * v[79];
	v[158] = v[559] * v[79];
	v[154] = v[560] * v[79];
	v[130] = lambda[0] * v[260];
	v[568] = ei1A[2] * v[130];
	v[562] = v[130] * v[79];
	v[561] = v[130] * v[197];
	v[419] = -(v[130] * v[522]);
	v[453] = v[197] * v[419];
	v[403] = ei1A[0] * v[561];
	v[347] = ei1A[2] * v[561];
	v[344] = v[196] * v[568];
	v[451] = (v[344] + v[345])*v[508];
	v[165] = v[419] * v[79];
	v[163] = ei1A[0] * v[562];
	v[155] = ei1A[2] * v[562];
	v[133] = lambda[0] * v[273];
	v[564] = ei1A[2] * v[133];
	v[563] = ei1A[1] * v[133];
	v[569] = v[505] * (ei1A[0] * v[130] + v[563]);
	v[432] = -(v[133] * v[521]);
	v[454] = v[197] * v[432];
	v[402] = v[197] * v[563];
	v[404] = v[402] + v[403];
	v[433] = v[195] * v[569];
	v[381] = v[197] * v[564];
	v[378] = v[196] * v[564];
	v[456] = v[195] * v[508] * (v[559] + v[564]);
	v[455] = v[510] * (v[127] * (v[365] + v[366]) + v[133] * (-ei1A[1] + ei1A[2] * v[501] - v[507] * v[521]) + v[130] *
		(ei1A[0] + ei1A[2] * v[505] - v[507] * v[522]));
	v[166] = v[432] * v[79];
	v[164] = v[563] * v[79];
	v[159] = v[564] * v[79];
	v[137] = lambda[0] * v[283];
	v[474] = (8e0*v[137]) / Power(v[170], 3);
	v[485] = v[474] * v[517];
	v[578] = v[485] + v[487];
	v[479] = v[474] * v[514];
	v[577] = v[479] + v[482];
	v[175] = v[137] * v[520];
	v[491] = v[175] + lambda[0] * v[582];
	v[478] = v[175] + lambda[0] * v[580];
	v[470] = v[175] + lambda[0] * v[575];
	v[138] = lambda[0] * v[302];
	v[139] = lambda[0] * v[321];
	v[140] = lambda[0] * v[340];
	v[141] = v[154] + v[155];
	v[151] = v[127] * v[368] + v[130] * v[369] + v[133] * v[370];
	v[421] = (8e0*v[151]) / Power(v[156], 3);
	v[452] = v[421] * v[507];
	v[574] = v[452] + v[455];
	v[435] = v[197] * v[418] + v[574];
	v[431] = v[196] * v[418] + v[421] * v[504] + v[510] * (v[130] * (v[363] + v[366]) + v[133] * (ei1A[2] + ei1A[1] * v[501]
		- v[504] * v[521]) + v[127] * (-ei1A[0] + ei1A[1] * v[508] - v[504] * v[523]));
	v[161] = v[151] * v[510];
	v[450] = v[161] + v[165] + v[166];
	v[430] = v[160] - v[165] + v[450];
	v[417] = v[165] - v[166] + v[430];
	v[152] = v[158] + v[159];
	v[153] = v[163] + v[164];
	residual[0] = v[154] - v[155] + v[417] * v[502] + v[153] * v[505] + v[152] * v[508];
	residual[1] = -v[158] + v[159] + v[153] * v[501] + v[430] * v[504] + v[141] * v[508];
	residual[2] = v[163] - v[164] + v[152] * v[501] + v[141] * v[505] + v[450] * v[507];
	residual[3] = v[470] * v[512] + v[140] * v[515] + v[138] * v[518] + lambda[0] * v[579];
	residual[4] = v[140] * v[511] + v[478] * v[514] + v[139] * v[518] + lambda[0] * v[581];
	residual[5] = v[138] * v[511] + v[139] * v[515] + v[491] * v[517] + lambda[0] * v[583];
	residual[6] = ei2B[0] * (v[107] * v[120] + v[104] * v[565]) + ei2B[1] * (v[113] * v[121] + v[104] * v[566]) + ei2B[2] *
		(v[119] * v[122] + v[104] * v[567]);
	stiffness[0][0] = 2e0*v[417] + v[433] + v[456] + v[502] * (v[195] * (v[418] + v[419]) + v[421] * v[502] + v[510] *
		(v[133] * (v[363] + v[365]) + v[130] * (-ei1A[2] + ei1A[0] * v[505] - v[502] * v[522]) + v[127] * (ei1A[1]
		+ ei1A[0] * v[508] - v[502] * v[523]))) + v[195] * (v[560] - v[568]);
	stiffness[0][1] = 0.5e0*v[153] - v[344] + v[345] + (v[196] * v[419] + v[431])*v[502] + (v[378] + v[379])*v[508]
		+ v[196] * v[569];
	stiffness[0][2] = 0.5e0*v[152] - v[347] + v[348] + (v[435] + v[453])*v[502] + v[404] * v[505] + (v[381] + v[382]
		)*v[508];
	stiffness[0][3] = -v[350] + v[351] + v[407] * v[505] + v[386] * v[508] + v[502] * (v[437] + v[570]);
	stiffness[0][4] = -v[353] + v[354] + v[410] * v[505] + v[389] * v[508] + v[502] * (v[440] + v[571]);
	stiffness[0][5] = -v[356] + v[357] + v[413] * v[505] + v[392] * v[508] + v[502] * (v[444] + v[572]);
	stiffness[0][6] = -v[359] + v[360] + v[416] * v[505] + v[395] * v[508] + v[502] * (v[447] + v[573]);
	stiffness[1][1] = v[378] - v[379] + 2e0*v[430] + v[433] + v[451] + (v[431] + v[196] * v[432])*v[504];
	stiffness[1][2] = 0.5e0*v[141] + v[381] - v[382] + v[404] * v[501] + (v[435] + v[454])*v[504] + (v[347] + v[348]
		)*v[508];
	stiffness[1][3] = v[384] - v[385] + v[407] * v[501] + (v[437] + v[438] + v[459])*v[504] + v[352] * v[508];
	stiffness[1][4] = v[387] - v[388] + v[410] * v[501] + (v[440] + v[442] + v[462])*v[504] + v[355] * v[508];
	stiffness[1][5] = v[390] - v[391] + v[413] * v[501] + (v[444] + v[445] + v[465])*v[504] + v[358] * v[508];
	stiffness[1][6] = v[393] - v[394] + v[416] * v[501] + (v[447] + v[448] + v[468])*v[504] + v[361] * v[508];
	stiffness[2][2] = -v[402] + v[403] + 2e0*v[450] + v[451] + v[456] + v[507] * (v[453] + v[454] + v[574]);
	stiffness[2][3] = -v[405] + v[406] + v[386] * v[501] + v[352] * v[505] + v[507] * (v[459] + v[570]);
	stiffness[2][4] = -v[408] + v[409] + v[389] * v[501] + v[355] * v[505] + v[507] * (v[462] + v[571]);
	stiffness[2][5] = -v[411] + v[412] + v[392] * v[501] + v[358] * v[505] + v[507] * (v[465] + v[572]);
	stiffness[2][6] = -v[414] + v[415] + v[395] * v[501] + v[361] * v[505] + v[507] * (v[468] + v[573]);
	stiffness[3][3] = 2e0*v[470] + v[483] + v[494] + (-v[550] + v[552])*v[554] + v[512] * (v[474] * v[512] + lambda[0] *
		(v[199] * (v[471] + v[472]) + v[520] * (ei2B[0] * (v[278] + v[280]) + ei2B[2] * (-v[121] + v[120] * v[518]
		- v[512] * v[555]) + ei2B[1] * (v[122] + v[120] * v[515] - v[512] * v[557]))));
	stiffness[3][4] = 0.5e0*v[140] + lambda[0] * (v[315] - v[316] + (v[296] + v[297])*v[518] + v[200] * v[576]) + v[512] *
		(lambda[0] * (v[200] * v[471] + v[481]) + v[577]);
	stiffness[3][5] = 0.5e0*v[138] + v[339] * v[515] + lambda[0] * (v[318] - v[319] + (v[299] + v[300])*v[518]) + v[512] *
		(lambda[0] * (v[486] + v[493]) + v[578]);
	stiffness[3][6] = v[340] * v[515] + v[302] * v[518] + v[512] * (v[489] + v[575]) + v[579];
	stiffness[4][4] = lambda[0] * (-v[296] + v[297]) + 2e0*v[478] + v[483] + v[495] + v[514] * (lambda[0] *
		(v[200] * v[480] + v[481]) + v[577]);
	stiffness[4][5] = 0.5e0*v[139] + v[339] * v[511] + lambda[0] * (-v[299] + v[300] + (v[318] + v[319])*v[518])
		+ v[514] * (lambda[0] * (v[486] + v[492]) + v[578]);
	stiffness[4][6] = v[340] * v[511] + v[321] * v[518] + v[514] * (v[489] + v[580]) + v[581];
	stiffness[5][5] = lambda[0] * (v[337] - v[338]) + 2e0*v[491] + v[494] + v[495] + v[517] * (lambda[0] * (v[492] + v[493]
		) + v[578]);
	stiffness[5][6] = v[302] * v[511] + v[321] * v[515] + v[517] * (v[489] + v[582]) + v[583];
	stiffness[6][6] = 0e0;
	for (i01 = 1; i01<7; i01++){
		for (i02 = 0; i02<i01; i02++){
			stiffness[i01][i02] = stiffness[i02][i01];
		}
	};
};

//Calcula contribuições do residuo e operador tangente - gerado no AceGen - sem usar SMSD
void UniversalJoint::EvaluateUniversalJointContribution2(double *v, double *residual
	, double **stiffness, double *alphaA, double *alphaB, double *ei1A
	, double *ei2B, double *lambda)
{
	int i01; int i02;
	v[271] = 0.5e0*alphaA[2];
	v[270] = Power(alphaB[2], 2);
	v[269] = 0.5e0*alphaB[0] * alphaB[2];
	v[268] = 0.5e0*alphaB[1];
	v[267] = Power(alphaB[1], 2);
	v[272] = v[267] + v[270];
	v[266] = alphaB[0] * v[268];
	v[265] = Power(alphaB[0], 2);
	v[264] = Power(alphaA[2], 2);
	v[263] = alphaA[0] * v[271];
	v[262] = 0.5e0*alphaA[1];
	v[261] = Power(alphaA[1], 2);
	v[260] = alphaA[0] * v[262];
	v[259] = Power(alphaA[0], 2);
	v[275] = v[259] + v[261];
	v[101] = alphaA[2] * v[262];
	v[120] = alphaB[2] * v[268];
	v[88] = 4e0 / (4e0 + v[264] + v[275]);
	v[278] = ei1A[2] * v[88];
	v[277] = ei1A[1] * v[88];
	v[276] = ei1A[0] * v[88];
	v[274] = -0.5e0*v[88];
	v[104] = -(v[271] * v[88]);
	v[105] = v[262] * v[88];
	v[106] = alphaA[0] * v[274];
	v[107] = 4e0 / (4e0 + v[265] + v[272]);
	v[273] = -0.5e0*v[107];
	v[110] = 1e0 + v[272] * v[273];
	v[111] = v[107] * (-alphaB[2] + v[266]);
	v[112] = v[107] * (alphaB[1] + v[269]);
	v[131] = ei2B[0] * v[110] + ei2B[1] * v[111] + ei2B[2] * v[112];
	v[182] = -(v[131] * v[88]);
	v[179] = -(v[106] * v[131]);
	v[145] = v[107] * v[131];
	v[114] = v[107] * (alphaB[2] + v[266]);
	v[116] = 1e0 + (v[265] + v[270])*v[273];
	v[117] = v[107] * (-alphaB[0] + v[120]);
	v[133] = ei2B[0] * v[114] + ei2B[1] * v[116] + ei2B[2] * v[117];
	v[184] = v[133] * v[88];
	v[180] = v[105] * v[133];
	v[148] = -(v[107] * v[133]);
	v[119] = v[107] * (-alphaB[1] + v[269]);
	v[121] = v[107] * (alphaB[0] + v[120]);
	v[122] = 1e0 + (v[265] + v[267])*v[273];
	v[134] = ei2B[0] * v[119] + ei2B[1] * v[121] + ei2B[2] * v[122];
	v[169] = -(v[104] * v[134]);
	v[166] = -(v[134] * v[88]);
	v[143] = v[107] * v[134];
	v[123] = alphaB[2] * v[273];
	v[142] = v[123] * v[134];
	v[124] = -(alphaB[1] * v[273]);
	v[146] = -(v[124] * v[133]);
	v[125] = alphaB[0] * v[273];
	v[147] = v[125] * v[131];
	v[126] = ei1A[2] * (1e0 + v[274] * v[275]) + (-alphaA[1] + v[263])*v[276] + (alphaA[0] + v[101])*v[277];
	v[234] = -(v[123] * v[126]);
	v[230] = -(v[107] * v[126]);
	v[132] = v[104] * v[126];
	v[130] = v[126] * v[88];
	v[127] = ei1A[1] * (1e0 + (v[259] + v[264])*v[274]) + (alphaA[2] + v[260])*v[276] + (-alphaA[0] + v[101])*v[278];
	v[246] = v[124] * v[127];
	v[241] = v[107] * v[127];
	v[204] = -(v[106] * v[127]) - v[130];
	v[138] = -(v[105] * v[127]);
	v[195] = v[132] + v[138];
	v[136] = -(v[127] * v[88]);
	v[213] = -(v[106] * v[126]) - v[136];
	v[128] = ei1A[0] * (1e0 + (v[261] + v[264])*v[274]) + (-alphaA[2] + v[260])*v[277] + (alphaA[1] + v[263])*v[278];
	v[247] = -(v[125] * v[128]);
	v[245] = -(v[127] * (v[124] * v[134] - v[145])) - v[126] * (v[146] + v[147]) + v[128] * (v[125] * v[134] + v[148]);
	v[243] = -(v[107] * v[128]);
	v[233] = v[128] * (v[125] * v[133] + v[143]) + v[126] * (v[123] * v[133] - v[145]) - v[127] * (v[142] + v[147]);
	v[220] = -(v[127] * (v[124] * v[131] + v[143])) - v[128] * (v[142] + v[146]) - v[126] * (-(v[123] * v[131]) + v[148]);
	v[194] = v[105] * v[128] + v[130];
	v[199] = -(v[131] * v[194]) + v[133] * v[195];
	v[193] = -(v[104] * v[128]) + v[136];
	v[198] = -(v[133] * v[193]) + v[134] * v[194];
	v[197] = v[131] * v[193] - v[134] * v[195];
	v[153] = -(v[134] * v[193]) - v[133] * v[194] - v[131] * v[195];
	v[139] = v[106] * v[128];
	v[211] = v[138] + v[139];
	v[215] = v[131] * v[211] - v[134] * v[213];
	v[203] = v[132] + v[139];
	v[208] = -(v[131] * v[203]) + v[133] * v[204];
	v[137] = v[128] * v[88];
	v[212] = v[105] * v[126] - v[137];
	v[217] = -(v[131] * v[212]) + v[133] * v[213];
	v[216] = -(v[133] * v[211]) + v[134] * v[212];
	v[202] = -(v[104] * v[127]) + v[137];
	v[207] = -(v[133] * v[202]) + v[134] * v[203];
	v[206] = v[131] * v[202] - v[134] * v[204];
	v[178] = -(v[134] * v[211]) - v[133] * v[212] - v[131] * v[213];
	v[168] = -(v[134] * v[202]) - v[133] * v[203] - v[131] * v[204];
	v[151] = v[104] * v[131] + v[184];
	v[152] = -(v[105] * v[131]) + v[166];
	v[156] = -(v[127] * v[151]) + v[126] * v[152];
	v[154] = v[169] + v[180];
	v[158] = v[128] * v[151] - v[126] * v[154];
	v[157] = -(v[128] * v[152]) + v[127] * v[154];
	v[159] = -(v[127] * v[131]);
	v[160] = v[128] * v[133];
	v[162] = v[126] * v[131];
	v[163] = -(v[128] * v[134]);
	v[165] = v[104] * v[133] + v[182];
	v[167] = v[106] * v[133] - v[166];
	v[172] = v[128] * v[165] - v[126] * v[167];
	v[170] = v[169] + v[179];
	v[174] = -(v[127] * v[165]) + v[126] * v[170];
	v[173] = v[127] * v[167] - v[128] * v[170];
	v[175] = -(v[126] * v[133]);
	v[176] = v[127] * v[134];
	v[181] = v[179] + v[180];
	v[183] = -(v[105] * v[134]) - v[182];
	v[185] = v[106] * v[134] - v[184];
	v[221] = v[234] + v[246];
	v[222] = v[123] * v[128] + v[241];
	v[225] = -(v[134] * v[221]) + v[131] * v[222];
	v[223] = -(v[124] * v[128]) + v[230];
	v[227] = -(v[133] * v[222]) + v[134] * v[223];
	v[226] = v[133] * v[221] - v[131] * v[223];
	v[231] = v[125] * v[127] - v[230];
	v[232] = v[123] * v[127] + v[243];
	v[237] = -(v[134] * v[231]) + v[131] * v[232];
	v[235] = v[234] + v[247];
	v[239] = -(v[133] * v[232]) + v[134] * v[235];
	v[238] = v[133] * v[231] - v[131] * v[235];
	v[242] = v[125] * v[126] - v[241];
	v[244] = -(v[124] * v[126]) - v[243];
	v[248] = v[246] + v[247];
	residual[0] = lambda[0] * v[153];
	residual[1] = lambda[0] * v[168];
	residual[2] = lambda[0] * v[178];
	residual[3] = lambda[0] * v[220];
	residual[4] = lambda[0] * v[233];
	residual[5] = lambda[0] * v[245];
	residual[6] = ei2B[0] * (v[119] * v[126] + v[114] * v[127] + v[110] * v[128]) + ei2B[1] * (v[121] * v[126]
		+ v[116] * v[127] + v[111] * v[128]) + ei2B[2] * (v[122] * v[126] + v[117] * v[127] + v[112] * v[128]);
	stiffness[0][0] = lambda[0] * (v[106] * v[153] - v[105] * v[157] - v[104] * v[158] + v[156] * v[88]);
	stiffness[0][1] = lambda[0] * (v[104] * v[156] - v[106] * v[157] + (v[158] - 0.5e0*(alphaA[1] * v[153] + v[159]
		+ v[160]))*v[88]);
	stiffness[0][2] = lambda[0] * (v[105] * v[156] + v[106] * v[158] + (v[157] - 0.5e0*(alphaA[2] * v[153] - v[162]
		- v[163]))*v[88]);
	stiffness[0][3] = lambda[0] * (-(v[123] * v[197]) + v[107] * v[198] - v[124] * v[199]);
	stiffness[0][4] = lambda[0] * (v[107] * v[197] + v[123] * v[198] - v[125] * v[199]);
	stiffness[0][5] = lambda[0] * (v[125] * v[197] + v[124] * v[198] + v[107] * v[199]);
	stiffness[0][6] = v[153];
	stiffness[1][1] = lambda[0] * (-(v[105] * v[168]) - v[106] * v[173] + v[104] * v[174] + v[172] * v[88]);
	stiffness[1][2] = lambda[0] * (v[106] * v[172] + v[105] * v[174] + (v[173] - 0.5e0*(alphaA[2] * v[168] + v[175]
		+ v[176]))*v[88]);
	stiffness[1][3] = lambda[0] * (-(v[123] * v[206]) + v[107] * v[207] - v[124] * v[208]);
	stiffness[1][4] = lambda[0] * (v[107] * v[206] + v[123] * v[207] - v[125] * v[208]);
	stiffness[1][5] = lambda[0] * (v[125] * v[206] + v[124] * v[207] + v[107] * v[208]);
	stiffness[1][6] = v[168];
	stiffness[2][2] = lambda[0] * (v[104] * v[178] + v[105] * (-(v[127] * v[181]) + v[126] * v[183]) + v[106] *
		(v[128] * v[181] - v[126] * v[185]) + (-(v[128] * v[183]) + v[127] * v[185])*v[88]);
	stiffness[2][3] = lambda[0] * (-(v[123] * v[215]) + v[107] * v[216] - v[124] * v[217]);
	stiffness[2][4] = lambda[0] * (v[107] * v[215] + v[123] * v[216] - v[125] * v[217]);
	stiffness[2][5] = lambda[0] * (v[125] * v[215] + v[124] * v[216] + v[107] * v[217]);
	stiffness[2][6] = v[178];
	stiffness[3][3] = lambda[0] * (v[125] * v[220] - v[123] * v[225] - v[124] * v[226] + v[107] * v[227]);
	stiffness[3][4] = lambda[0] * (v[107] * (-0.5e0*(-v[159] - v[160] + alphaB[1] * v[220]) + v[225]) - v[125] * v[226]
		+ v[123] * v[227]);
	stiffness[3][5] = lambda[0] * (v[125] * v[225] + v[107] * (-0.5e0*(v[162] + v[163] + alphaB[2] * v[220]) + v[226])
		+ v[124] * v[227]);
	stiffness[3][6] = v[220];
	stiffness[4][4] = lambda[0] * (-(v[124] * v[233]) + v[107] * v[237] - v[125] * v[238] + v[123] * v[239]);
	stiffness[4][5] = lambda[0] * (v[125] * v[237] + v[107] * (-0.5e0*(-v[175] - v[176] + alphaB[2] * v[233]) + v[238])
		+ v[124] * v[239]);
	stiffness[4][6] = v[233];
	stiffness[5][5] = lambda[0] * (v[107] * (v[133] * v[242] - v[131] * v[244]) + v[123] * v[245] + v[125] * (-
		(v[134] * v[242]) + v[131] * v[248]) + v[124] * (v[134] * v[244] - v[133] * v[248]));
	stiffness[5][6] = v[245];
	stiffness[6][6] = 0e0;
	for (i01 = 1; i01<7; i01++){
		for (i02 = 0; i02<i01; i02++){
			stiffness[i01][i02] = stiffness[i02][i01];
		}
	};
};


