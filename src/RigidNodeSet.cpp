#include "RigidNodeSet.h"
#include"Database.h"

//Variáveis globais
extern
Database db;

RigidNodeSet::RigidNodeSet()
{
	I3 = Matrix(3, 3);
	I3(0, 0) = 1.0;
	I3(1, 1) = 1.0;
	I3(2, 2) = 1.0;
	residual1 = NULL;
	stiffness1 = NULL;
	residual2= NULL;
	stiffness2 = NULL;
	distancei = NULL;
	active_lambda = NULL;
	GLs = NULL;
	lambda = NULL;
	copy_lambda = NULL;

	flag_rotation = true;
}

RigidNodeSet::~RigidNodeSet()
{
	if (active_lambda != NULL)
		delete[]active_lambda;
	if (GLs != NULL)
		delete[]GLs;
	if (lambda != NULL)
		delete[]lambda;
	if (copy_lambda != NULL)
		delete[]copy_lambda;
	if (distancei != NULL)
	{
		for (int i = 0; i < this->n_nodes_set; i++)
			delete distancei[i];
		delete[] distancei;
	}
	if (residual1 != NULL)
	{
		for (int i = 0; i < this->n_nodes_set; i++)
			delete residual1[i];
		delete[] residual1;
	}
	if (stiffness1 != NULL)
	{
		for (int i = 0; i < this->n_nodes_set; i++)
			delete stiffness1[i];
		delete[] stiffness1;
	}
	if (residual2 != NULL)
	{
		for (int i = 0; i < this->n_nodes_set; i++)
			delete residual2[i];
		delete[] residual2;
	}
	if (stiffness2 != NULL)
	{
		for (int i = 0; i < this->n_nodes_set; i++)
			delete stiffness2[i];
		delete[] stiffness2;
	}
}

//Leitura
bool RigidNodeSet::Read(FILE *f)
{
	char s[1000];
	fscanf(f, "%s", s);
	number = atoi(s);	//número
	fscanf(f, "%s", s);
	if (!strcmp(s, "PilotNode"))
	{
		fscanf(f, "%s", s);
		pilot_node = atoi(s);	//pilot node
	}
	else
		return false;
	fscanf(f, "%s", s);
	if (!strcmp(s, "NodeSet"))
	{
		fscanf(f, "%s", s);
		node_set_ID = atoi(s);	//node set
	}
	else
		return false;
	//Salva a posição (stream)
	fpos_t pos;
	fgetpos(f, &pos);
	fscanf(f, "%s", s);
	if (!strcmp(s, "RotationDOF"))
	{
		fscanf(f, "%s", s);
		flag_rotation = (bool)atoi(s);	//flag rotation
	}
	else
	{
		fsetpos(f, &pos);
		flag_rotation = true;
	}

	//Salva a posição (stream)
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
void RigidNodeSet::Write(FILE *f)
{
	fprintf(f, "RigidNodeSet\t%d\tPilotNode\t%d\tNodeSet\t%d\tRotationDOF\t%d\t", number, pilot_node, node_set_ID, (bool)flag_rotation);
	bool_table.Write(f);
}

//Escreve no monitor do SpecialConstraint
void RigidNodeSet::WriteMonitor(FILE *f, bool first_record, double time)
{
	//Cabeçalho
	if (first_record == true)
		fprintf(f, "TIME\tLAGRANGE\n");
	//Informações a serem salvas
	fprintf(f, "%.6e\t", time);
	for (int i = 0; i < n_GL; i++)
	{
		fprintf(f, "%.6e\t", lambda[i]);
	}
	fprintf(f, "\n");
}

void RigidNodeSet::WriteVTK_XMLRender(FILE *f)
{
	if (db.post_files->WriteSpecialConstraints_flag == true)
	{
		
		int i, offsets;			            /*numero de pontos*/
		int tamanho = n_nodes_set+1;		/*tamanho dos vetores*/
		int m = tamanho;                    /*número de vértices*/
		float *x, *y, *z;					/*pontos*/

		x = new float[tamanho];
		y = new float[tamanho];
		z = new float[tamanho];

		//Pilot node
		x[0] = (float)db.nodes[pilot_node - 1]->copy_coordinates[0];
		y[0] = (float)db.nodes[pilot_node - 1]->copy_coordinates[1];
		z[0] = (float)db.nodes[pilot_node - 1]->copy_coordinates[2];

		//Slave nodes
		for (int ei = 0; ei < n_nodes_set; ei++)
		{
			x[ei+1] = (float)db.nodes[db.node_sets[node_set_ID - 1]->node_list[ei] - 1]->copy_coordinates[0];
			y[ei+1] = (float)db.nodes[db.node_sets[node_set_ID - 1]->node_list[ei] - 1]->copy_coordinates[1];
			z[ei+1] = (float)db.nodes[db.node_sets[node_set_ID - 1]->node_list[ei] - 1]->copy_coordinates[2];
		}

		fprintf(f, "     <Piece NumberOfPoints = \"%d\" NumberOfCells = \"%d\">\n", tamanho, (m - 1) + m);
		fprintf(f, "         <Points>\n");
		fprintf(f, "             <DataArray type = \"Float32\" NumberOfComponents = \"3\" Format = \"ascii\">\n");
		for (i = 0; i<tamanho; i++)
		{
			fprintf(f, "                 %.6f %.6f %.6f\n", x[i], y[i], z[i]);
		}
		fprintf(f, "             </DataArray>\n");
		fprintf(f, "         </Points>\n");
		fprintf(f, "         <Cells>\n");
		fprintf(f, "             <DataArray type = \"Int32\" Name = \"connectivity\" Format = \"ascii\">\n");
		for (i = 0; i<m; i++)
		{
			fprintf(f, "                  %d\n", i);
		}
		for (i = 1; i<m; i++)
		{
			fprintf(f, "                  %d %d\n", 0, i);
		}
		fprintf(f, "             </DataArray>\n");
		fprintf(f, "             <DataArray type = \"Int32\" Name = \"types\" Format = \"ascii\">\n");
		for (i = 0; i<m; i++)
		{
			fprintf(f, "                 1\n");
		}
		for (i = 1; i<m; i++)
		{
			fprintf(f, "                 3\n");
		}
		fprintf(f, "             </DataArray>\n");
		fprintf(f, "             <DataArray type = \"Int32\" Name = \"offsets\" Format = \"ascii\">\n");
		offsets = 0;
		for (i = 0; i<m; i++)
		{
			offsets = offsets + 1;
			fprintf(f, "                 %d\n", offsets);
		}
		for (i = 1; i<m; i++)
		{
			offsets = offsets + 2;
			fprintf(f, "                 %d\n", offsets);
		}
		fprintf(f, "             </DataArray>\n");
		fprintf(f, "         </Cells>\n");
		fprintf(f, "     </Piece>\n");
		
		delete[]x;
		delete[]y;
		delete[]z;
	}
}

//Checa inconsistências no SC para evitar erros de execução
bool RigidNodeSet::Check()
{
	if (pilot_node > db.number_nodes)
		return false;
	if (node_set_ID > db.number_node_sets)
		return false;

	//Checagem das condições iniciais
	int temp_node1 = 0;
	int temp_node2 = 0;
	for (int i = 0; i < db.number_IC; i++)
	{
		temp_node1 = db.IC[i]->node;
		for (int nn = 0; nn < db.node_sets[node_set_ID - 1]->n_nodes; nn++)
		{
			temp_node2 = db.node_sets[node_set_ID - 1]->node_list[nn];
			if (temp_node1 == temp_node2)
			{
				db.myprintf("Warning in Special Constraint %d.\nInitial Condition %d was prescribed to node %d (slave), leading to ignoring some of its components.\n", number, db.IC[i]->number, db.IC[i]->node);
			}
		}
		
	}
	return true;
}

//Montagem dos resíduos e rigidez tangente
void RigidNodeSet::Mount()
{
	Matrix r(3);
	Matrix distanceip(3);
	Matrix xpA(3);
	Matrix xpB(3);
	Matrix Qp(3, 3);
	Matrix M2(3, 3);
	Matrix M3(3, 3);
	double** mQp;
	double** mM2;
	double** mM3;
	mQp = new double*[3];
	mM2 = new double*[3];
	mM3 = new double*[3];
	for (int i = 0; i < 3; i++)
	{
		mQp[i] = new double[3];
		mM2[i] = new double[3];
		mM3[i] = new double[3];
	}
	Matrix temp_lamb(3);
	Matrix alpha(3);

	for (int i = 0; i < 3; i++)
	{
		xpA(i, 0) = db.nodes[pilot_node - 1]->copy_coordinates[i] + db.nodes[pilot_node - 1]->displacements[i];
		alpha(i, 0) = db.nodes[pilot_node - 1]->displacements[i+3];
	}
		
	for (int ei = 0; ei < n_nodes_set; ei++)
	{
		//Nesse vínculo a rigidez tangente não se modifica nunca. É montada no PreCalc()
		//Se o vínculo estiver ativo, realiza a montagem
		if (active_lambda[ei * 3 + 0 + n_nodes_set * 3] == 1 && active_lambda[ei * 3 + 1 + n_nodes_set * 3] == 1 && active_lambda[ei * 3 + 2 + n_nodes_set * 3] == 1)
		{
			for (int i = 0; i < 3; i++)
				xpB(i, 0) = db.nodes[db.node_sets[node_set_ID - 1]->node_list[ei] - 1]->copy_coordinates[i] + db.nodes[db.node_sets[node_set_ID - 1]->node_list[ei] - 1]->displacements[i];
			distanceip = xpB - xpA;
			
			temp_lamb(0, 0) = lambda[ei * 3 + 0 + n_nodes_set * 3];
			temp_lamb(1, 0) = lambda[ei * 3 + 1 + n_nodes_set * 3];
			temp_lamb(2, 0) = lambda[ei * 3 + 2 + n_nodes_set * 3];
			
			EvaluateM2M3(v, mM2, mM3, mQp, alpha.getMatrix(), distancei[ei]->getMatrix(), temp_lamb.getMatrix());
			M2.PtrToMatrix(mM2, 3);
			M3.PtrToMatrix(mM3, 3);
			Qp.PtrToMatrix(mQp, 3);
			//M2.print();
			//M3.print();
			r = distanceip - Qp*(*distancei[ei]);
			for (int i = 0; i < 3; i++)
			{
				(*residual1[ei])(i, 0) = -temp_lamb(i,0);
				(*residual1[ei])(i + 3, 0) = (transp(temp_lamb)*M2)(0,i);
				(*residual1[ei])(i + 6, 0) = +temp_lamb(i, 0);
				(*residual1[ei])(i + 9, 0) = r(i,0);
				for (int j = 0; j < 3; j++)
				{
					(*stiffness1[ei])(i + 3, j + 3) = -M3(i, j);
					(*stiffness1[ei])(i + 3, j + 9) = transp(M2)(i, j);
					(*stiffness1[ei])(i + 9, j + 3) = (M2)(i, j);
				}
			}
		}
		if (active_lambda[ei * 3 + 0] == 1 && active_lambda[ei * 3 + 1] == 1 && active_lambda[ei * 3 + 2] == 1)
		{
			//Parte 2 - Rotações
			for (int i = 0; i < 3; i++)
				r(i, 0) = db.nodes[pilot_node - 1]->displacements[i + 3] - db.nodes[db.node_sets[node_set_ID - 1]->node_list[ei] - 1]->displacements[i + 3];
			//Atualização do vetor de resíduos
			(*residual2[ei])(0, 0) = lambda[ei * 3 + 0];
			(*residual2[ei])(1, 0) = lambda[ei * 3 + 1];
			(*residual2[ei])(2, 0) = lambda[ei * 3 + 2];
			(*residual2[ei])(3, 0) = -lambda[ei * 3 + 0];
			(*residual2[ei])(4, 0) = -lambda[ei * 3 + 1];
			(*residual2[ei])(5, 0) = -lambda[ei * 3 + 2];
			(*residual2[ei])(6, 0) = r(0, 0);
			(*residual2[ei])(7, 0) = r(1, 0);
			(*residual2[ei])(8, 0) = r(2, 0);
		}
	}
	for (int i = 0; i < 3; i++)
	{
		delete[] mQp[i];
		delete[] mM2[i];
		delete[] mM3[i];
	}
	delete[] mQp;
	delete[] mM2;
	delete[] mM3;
	
}

//Preenche a contribuição do elemento nas matrizes globais
void RigidNodeSet::MountGlobal()
{
	//Variáveis temporárias para salvar a indexação global dos graus de liberdade a serem setados na matriz de rigidez global
	int GL_global_1 = 0;
	int GL_global_2 = 0;
	double anterior = 0;

	//Parte 1 - Deslocamentos do NodeSet e rotação do Pilot Node
	for (int ei = 0; ei < n_nodes_set; ei++)
	{
		//Se o vínculo estiver ativo, realiza a montagem
		if (active_lambda[ei * 3 + 0 + n_nodes_set * 3] == 1 && active_lambda[ei * 3 + 1 + n_nodes_set * 3] == 1 && active_lambda[ei * 3 + 2 + n_nodes_set * 3] == 1)
		{
			for (int i = 0; i < 12; i++)
			{
				//////////////MONTAGEM DO VETOR DE ESFORÇOS DESBALANCEADOS//////////////////
				//Toma do vetor de GL globais, a indexação de cada grau de liberdade global
				if (i < 3)//UA
					GL_global_1 = db.nodes[pilot_node - 1]->GLs[i];
				else
				{
					if (i < 6)//ROTA
					{
						GL_global_1 = db.nodes[pilot_node - 1]->GLs[i];
					}
					else
					{
						if (i < 9)//UB
							GL_global_1 = db.nodes[db.node_sets[node_set_ID - 1]->node_list[ei] - 1]->GLs[i - 6];
						else
						{
							//lambda
							GL_global_1 = GLs[n_nodes_set * 3 + ei * 3 + (i - 9)];
						}
					}
				}

				//Caso o grau de liberdade seja livre:
				if (GL_global_1 > 0)
				{
					anterior = db.global_P_A(GL_global_1 - 1, 0);
					db.global_P_A(GL_global_1 - 1, 0) = anterior + (*residual1[ei])(i, 0);
					anterior = db.global_I_A(GL_global_1 - 1, 0);
					db.global_I_A(GL_global_1 - 1, 0) = anterior + (*residual1[ei])(i, 0);
				}
				else
				{
					anterior = db.global_P_B(-GL_global_1 - 1, 0);
					db.global_P_B(-GL_global_1 - 1, 0) = anterior + (*residual1[ei])(i, 0);
				}
				for (int j = 0; j < 12; j++)
				{
					//////////////////////MONTAGEM DA MATRIZ DE RIGIDEZ/////////////////////////
					//Toma do vetor de GL globais, a indexação de cada grau de liberdade global
					if (j < 3)//UA
						GL_global_2 = db.nodes[pilot_node - 1]->GLs[j];
					else
					{
						if (j < 6)//ROTA
						{
							GL_global_2 = db.nodes[pilot_node - 1]->GLs[j];
						}
						else
						{
							if (j < 9)//UB
								GL_global_2 = db.nodes[db.node_sets[node_set_ID - 1]->node_list[ei] - 1]->GLs[j - 6];
							else
							{
								//lambda
								GL_global_2 = GLs[n_nodes_set * 3 + ei * 3 + (j - 9)];
							}
						}
					}

					//Caso os graus de liberdade sejam ambos livres (Matriz Kaa)
					if (GL_global_1 > 0 && GL_global_2 > 0)
						db.global_stiffness_AA.setValue(GL_global_1 - 1, GL_global_2 - 1, (*stiffness1[ei])(i, j));
					//Caso os graus de liberdade sejam ambos fixos (Matriz Kbb)
					if (GL_global_1 < 0 && GL_global_2 < 0)
						db.global_stiffness_BB.setValue(-GL_global_1 - 1, -GL_global_2 - 1, (*stiffness1[ei])(i, j));
					//Caso os graus de liberdade sejam livre e fixo (Matriz Kab)
					if (GL_global_1 > 0 && GL_global_2 < 0)
						db.global_stiffness_AB.setValue(GL_global_1 - 1, -GL_global_2 - 1, (*stiffness1[ei])(i, j));
					//Caso os graus de liberdade sejam fixo e livre (Matriz Kba)
					if (GL_global_1 < 0 && GL_global_2 > 0)
						db.global_stiffness_BA.setValue(-GL_global_1 - 1, GL_global_2 - 1, (*stiffness1[ei])(i, j));
				}
			}
		}
	}

	//Parte 2 - rotações
	for (int ei = 0; ei < n_nodes_set; ei++)
	{
		//Se o vínculo estiver ativo, realiza a montagem
		if (active_lambda[ei * 3 + 0] == 1 && active_lambda[ei * 3 + 1] == 1 && active_lambda[ei * 3 + 2] == 1)
		{
			for (int i = 0; i < 9; i++)
			{
				//////////////MONTAGEM DO VETOR DE ESFORÇOS DESBALANCEADOS//////////////////
				//Toma do vetor de GL globais, a indexação de cada grau de liberdade global
				if (i < 3)//ROTA
					GL_global_1 = db.nodes[pilot_node - 1]->GLs[i + 3];
				else
				{
					if (i < 6)//ROTB
						GL_global_1 = db.nodes[db.node_sets[node_set_ID - 1]->node_list[ei] - 1]->GLs[i];
					else
					{
						//lambda
						GL_global_1 = GLs[ei * 3 + (i - 6)];
					}
				}

				//Caso o grau de liberdade seja livre:
				if (GL_global_1 > 0)
				{
					anterior = db.global_P_A(GL_global_1 - 1, 0);
					db.global_P_A(GL_global_1 - 1, 0) = anterior + (*residual2[ei])(i, 0);
					anterior = db.global_I_A(GL_global_1 - 1, 0);
					db.global_I_A(GL_global_1 - 1, 0) = anterior + (*residual2[ei])(i, 0);
				}
				else
				{
					anterior = db.global_P_B(-GL_global_1 - 1, 0);
					db.global_P_B(-GL_global_1 - 1, 0) = anterior + (*residual2[ei])(i, 0);
				}
				for (int j = 0; j < 9; j++)
				{
					//////////////////////MONTAGEM DA MATRIZ DE RIGIDEZ/////////////////////////
					//Toma do vetor de GL globais, a indexação de cada grau de liberdade global
					if (j < 3)//ROTA
						GL_global_2 = db.nodes[pilot_node - 1]->GLs[j + 3];
					else
					{
						if (j<6)//ROTB
							GL_global_2 = db.nodes[db.node_sets[node_set_ID - 1]->node_list[ei] - 1]->GLs[j];
						else
						{
							//lambda
							GL_global_2 = GLs[ei * 3 + (j - 6)];
						}
					}

					//Caso os graus de liberdade sejam ambos livres (Matriz Kaa)
					if (GL_global_1 > 0 && GL_global_2 > 0)
						db.global_stiffness_AA.setValue(GL_global_1 - 1, GL_global_2 - 1, (*stiffness2[ei])(i, j));
					//Caso os graus de liberdade sejam ambos fixos (Matriz Kbb)
					if (GL_global_1 < 0 && GL_global_2 < 0)
						db.global_stiffness_BB.setValue(-GL_global_1 - 1, -GL_global_2 - 1, (*stiffness2[ei])(i, j));
					//Caso os graus de liberdade sejam livre e fixo (Matriz Kab)
					if (GL_global_1 > 0 && GL_global_2 < 0)
						db.global_stiffness_AB.setValue(GL_global_1 - 1, -GL_global_2 - 1, (*stiffness2[ei])(i, j));
					//Caso os graus de liberdade sejam fixo e livre (Matriz Kba)
					if (GL_global_1 < 0 && GL_global_2 > 0)
						db.global_stiffness_BA.setValue(-GL_global_1 - 1, GL_global_2 - 1, (*stiffness2[ei])(i, j));
				}
			}
		}
	}
}

//Pré-cálculo de variáveis que é feito uma única vez no início
void RigidNodeSet::PreCalc()
{
	Alloc();
	Matrix xpA(3);
	Matrix xpB(3);
	for (int i = 0; i < 3; i++)
		xpA(i, 0) = db.nodes[pilot_node - 1]->ref_coordinates[i];
	for (int ei = 0; ei < n_nodes_set; ei++)
	{
		for (int i = 0; i < 3; i++)
			xpB(i, 0) = db.nodes[db.node_sets[node_set_ID - 1]->node_list[ei] - 1]->ref_coordinates[i];
		(*distancei[ei]) = xpB - xpA;
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				(*stiffness1[ei])(i, j + 9) = -I3(i, j);
				(*stiffness1[ei])(i + 9, j) = -I3(i, j);
				(*stiffness1[ei])(i + 6, j + 9) = +I3(i, j);
				(*stiffness1[ei])(i + 9, j + 6) = +I3(i, j);

				(*stiffness2[ei])(i, j + 6) = I3(i, j);
				(*stiffness2[ei])(i + 3, j + 6) = -I3(i, j);
				(*stiffness2[ei])(i + 6, j) = I3(i, j);
				(*stiffness2[ei])(i + 6, j + 3) = -I3(i, j);
			}
		}
	}
}

void RigidNodeSet::ComputeInitialGuessDisplacements()
{
	//Se for o step de criação do vínculo inicializa condições iniciais
	if (bool_table.GetAt(db.current_solution_number - 1) == true)
	{
		//Condições iniciais - transferindo as condições iniciais do pilot node para os outros nós do rigid node set
		int  temp_node = 0;
		Matrix omega_pilot(3);
		Matrix velocity_slave(3);
		Matrix velocity_pilot(3);
		
		//Velocidade - pilot node
		velocity_pilot(0, 0) = db.nodes[pilot_node - 1]->displacements[0];
		velocity_pilot(1, 0) = db.nodes[pilot_node - 1]->displacements[1];
		velocity_pilot(2, 0) = db.nodes[pilot_node - 1]->displacements[2];
		
		//Velocidade angular - pilot node
		omega_pilot(0, 0) = db.nodes[pilot_node - 1]->displacements[3];
		omega_pilot(1, 0) = db.nodes[pilot_node - 1]->displacements[4];
		omega_pilot(2, 0) = db.nodes[pilot_node - 1]->displacements[5];
		

		for (int ei = 0; ei < n_nodes_set; ei++)
		{
			temp_node = db.node_sets[node_set_ID - 1]->node_list[ei];

			//Velocidades dos nós do node set
			velocity_slave = velocity_pilot + cross(omega_pilot, *distancei[ei]);
			
			//Copiando informações calculadas
			db.nodes[temp_node - 1]->displacements[0] = velocity_slave(0, 0);
			db.nodes[temp_node - 1]->displacements[1] = velocity_slave(1, 0);
			db.nodes[temp_node - 1]->displacements[2] = velocity_slave(2, 0);
			db.nodes[temp_node - 1]->displacements[3] = omega_pilot(0, 0);
			db.nodes[temp_node - 1]->displacements[4] = omega_pilot(1, 0);
			db.nodes[temp_node - 1]->displacements[5] = omega_pilot(2, 0);

		}
	}
}

//Computa efeito das condições iniciais nos nós da restrição
void RigidNodeSet::ComputeVelAccel()
{
	//Se for o step de criação do vínculo inicializa condições iniciais
	if (bool_table.GetAt(db.current_solution_number - 1) == true)
	{
		//Condições iniciais - transferindo as condições iniciais do pilot node para os outros nós do rigid node set
		int  temp_node = 0;
		Matrix omega_pilot(3);
		Matrix velocity_slave(3);
		Matrix velocity_pilot(3);
		Matrix acc_slave(3);
		Matrix acc_pilot(3);
		Matrix omega_acc_pilot(3);
		Matrix omega_acc_slave(3);

		//Velocidade - pilot node
		velocity_pilot(0, 0) = db.nodes[pilot_node - 1]->vel[0];
		velocity_pilot(1, 0) = db.nodes[pilot_node - 1]->vel[1];
		velocity_pilot(2, 0) = db.nodes[pilot_node - 1]->vel[2];
		//Aceleração - pilot node
		acc_pilot(0, 0) = db.nodes[pilot_node - 1]->accel[0];
		acc_pilot(1, 0) = db.nodes[pilot_node - 1]->accel[1];
		acc_pilot(2, 0) = db.nodes[pilot_node - 1]->accel[2];
		//Velocidade angular - pilot node
		omega_pilot(0, 0) = db.nodes[pilot_node - 1]->vel[3];
		omega_pilot(1, 0) = db.nodes[pilot_node - 1]->vel[4];
		omega_pilot(2, 0) = db.nodes[pilot_node - 1]->vel[5];
		//Aceleração angular - pilot node
		omega_acc_pilot(0, 0) = db.nodes[pilot_node - 1]->accel[3];
		omega_acc_pilot(1, 0) = db.nodes[pilot_node - 1]->accel[4];
		omega_acc_pilot(2, 0) = db.nodes[pilot_node - 1]->accel[5];

		Matrix xpA(3);
		Matrix xpB(3);
		Matrix distanceip(3);
		for (int i = 0; i < 3; i++)
			xpA(i, 0) = db.nodes[pilot_node - 1]->copy_coordinates[i] + db.nodes[pilot_node - 1]->displacements[i];

		for (int ei = 0; ei < n_nodes_set; ei++)
		{
			temp_node = db.node_sets[node_set_ID - 1]->node_list[ei];

			for (int i = 0; i < 3; i++)
				xpB(i, 0) = db.nodes[temp_node - 1]->copy_coordinates[i] + db.nodes[temp_node - 1]->displacements[i];
			distanceip = xpB - xpA;

			//Velocidades dos nós do node set
			velocity_slave = velocity_pilot + cross(omega_pilot, distanceip);
			//Aceleração dos nós do node set
			acc_slave = acc_pilot + cross(omega_acc_pilot, distanceip) + cross(omega_pilot, cross(omega_pilot, distanceip));

			//Copiando informações calculadas
			db.nodes[temp_node - 1]->vel[0] = velocity_slave(0, 0);
			db.nodes[temp_node - 1]->vel[1] = velocity_slave(1, 0);
			db.nodes[temp_node - 1]->vel[2] = velocity_slave(2, 0);
			db.nodes[temp_node - 1]->vel[3] = omega_pilot(0, 0);
			db.nodes[temp_node - 1]->vel[4] = omega_pilot(1, 0);
			db.nodes[temp_node - 1]->vel[5] = omega_pilot(2, 0);

			db.nodes[temp_node - 1]->accel[0] = acc_slave(0, 0);
			db.nodes[temp_node - 1]->accel[1] = acc_slave(1, 0);
			db.nodes[temp_node - 1]->accel[2] = acc_slave(2, 0);
			db.nodes[temp_node - 1]->accel[3] = omega_acc_pilot(0, 0);
			db.nodes[temp_node - 1]->accel[4] = omega_acc_pilot(1, 0);
			db.nodes[temp_node - 1]->accel[5] = omega_acc_pilot(2, 0);
		}
	}
}

//Salvando variáveis da configuração convergida
void RigidNodeSet::SaveLagrange()
{
	for (int i = 0; i < n_GL; i++)
		copy_lambda[i] = lambda[i];
	Matrix xpA(3);
	Matrix xpB(3);
	for (int i = 0; i < 3; i++)
		xpA(i, 0) = db.nodes[pilot_node - 1]->copy_coordinates[i];
	for (int ei = 0; ei < n_nodes_set; ei++)
	{
		for (int i = 0; i < 3; i++)
			xpB(i, 0) = db.nodes[db.node_sets[node_set_ID - 1]->node_list[ei] - 1]->copy_coordinates[i];
		(*distancei[ei]) = xpB - xpA;
	}
	/*for (int i = 0; i < n_GL; i++)
	{
		db.myprintf("Lambda %.6f\n", lambda[i]);
	}*/
}

//Checa quais multiplicadores de lagrange serão ativados,de acordo com a ativação dos GLs dos nós dos quais a special constraint participa
void RigidNodeSet::ActivateDOFs()
{
	//Ativa GLs do Pilot Node (tanto os de translação como os de rotação)
	for (int i = 0; i < 6; i++)
		db.nodes[pilot_node - 1]->active_GL[i] = 1;

	//Ativa GLs dos nós do node set
	for (int ei = 0; ei < n_nodes_set; ei++)
	{
		db.nodes[db.node_sets[node_set_ID - 1]->node_list[ei] - 1]->active_GL[0] = 1;
		db.nodes[db.node_sets[node_set_ID - 1]->node_list[ei] - 1]->active_GL[1] = 1;
		db.nodes[db.node_sets[node_set_ID - 1]->node_list[ei] - 1]->active_GL[2] = 1;

		if (flag_rotation)
		{
			db.nodes[db.node_sets[node_set_ID - 1]->node_list[ei] - 1]->active_GL[3] = 1;
			db.nodes[db.node_sets[node_set_ID - 1]->node_list[ei] - 1]->active_GL[4] = 1;
			db.nodes[db.node_sets[node_set_ID - 1]->node_list[ei] - 1]->active_GL[5] = 1;
		}
	}

	if (bool_table.GetAt(db.current_solution_number - 1))
	{
		//Ativa GLs dos nós do node set e multiplicadores de Lagrange
		for (int ei = 0; ei < n_nodes_set; ei++)
		{
			if (flag_rotation)
			{
				active_lambda[ei * 3 + 0] = 1;
				active_lambda[ei * 3 + 1] = 1;
				active_lambda[ei * 3 + 2] = 1;
			}

			active_lambda[ei * 3 + 0 + n_nodes_set * 3] = 1;
			active_lambda[ei * 3 + 1 + n_nodes_set * 3] = 1;
			active_lambda[ei * 3 + 2 + n_nodes_set * 3] = 1;
		}
	}
	else
	{
		//Desativa GLs dos multiplicadores de Lagrange
		for (int ei = 0; ei < n_nodes_set; ei++)
		{ 
			if (flag_rotation)
			{
				active_lambda[ei * 3 + 0] = 0;
				active_lambda[ei * 3 + 1] = 0;
				active_lambda[ei * 3 + 2] = 0;
			}

			active_lambda[ei * 3 + 0 + n_nodes_set * 3] = 0;
			active_lambda[ei * 3 + 1 + n_nodes_set * 3] = 0;
			active_lambda[ei * 3 + 2 + n_nodes_set * 3] = 0;
		}
	}
}
void RigidNodeSet::Alloc()
{
	n_nodes_set = db.node_sets[node_set_ID - 1]->n_nodes;
	n_GL = n_nodes_set * 3 + n_nodes_set * 3;						//Número de graus de liberdade - multiplicadores de Lagrange 
	//n_nodes_set * 3 - lambdas devidos às vinculações de rotações do NodeSet
	//n_nodes_set * 3 - lambdas devidos às vinculações de deslocamentos do NodeSet
	active_lambda = new int[n_GL];
	GLs = new int[n_GL];
	lambda = new double[n_GL];
	copy_lambda = new double[n_GL];

	for (int i = 0; i < n_GL; i++)
	{
		active_lambda[i] = 0;
		GLs[i] = 0;
		//Chute inicial para os lambdas: valores unitários
		lambda[i] = 1.0;//Valor não nulo para evitar divergência na primeira iteração
		copy_lambda[i] = 1.0;
	}

	distancei = new Matrix*[n_nodes_set];
	residual1 = new Matrix*[n_nodes_set];
	stiffness1 = new Matrix*[n_nodes_set];
	residual2 = new Matrix*[n_nodes_set];
	stiffness2 = new Matrix*[n_nodes_set];
	for (int ei = 0; ei < n_nodes_set; ei++)
	{
		residual1[ei] = new Matrix(12);
		stiffness1[ei] = new Matrix(12, 12);
		residual2[ei] = new Matrix(9);
		stiffness2[ei] = new Matrix(9, 9);
		distancei[ei] = new Matrix(3);
	}
}
void RigidNodeSet::EvaluateM2M3(double* v, double** M2, double** M3, double** Qp, double* alphap, double* di, double* lamb)
{
	v[120] = Power(alphap[2], 2);
	v[119] = 0.5e0*alphap[0] * alphap[2];
	v[118] = 0.5e0*alphap[1];
	v[117] = Power(alphap[1], 2);
	v[121] = v[117] + v[120];
	v[116] = alphap[0] * v[118];
	v[115] = Power(alphap[0], 2);
	v[50] = alphap[2] * v[118];
	v[37] = 4e0 / (4e0 + v[115] + v[121]);
	v[122] = -0.5e0*v[37];
	v[109] = lamb[1] * v[37];
	v[104] = -(lamb[0] * v[37]);
	v[94] = -(lamb[2] * v[37]);
	v[40] = 1e0 + v[121] * v[122];
	v[41] = (-alphap[2] + v[116])*v[37];
	v[42] = (alphap[1] + v[119])*v[37];
	v[64] = -(di[0] * v[40]) - di[1] * v[41] - di[2] * v[42];
	v[85] = lamb[2] * v[64];
	v[80] = -(lamb[1] * v[64]);
	v[69] = v[37] * v[64];
	v[44] = (alphap[2] + v[116])*v[37];
	v[46] = 1e0 + (v[115] + v[120])*v[122];
	v[47] = v[37] * (-alphap[0] + v[50]);
	v[58] = di[0] * v[44] + di[1] * v[46] + di[2] * v[47];
	v[100] = lamb[2] * v[58];
	v[81] = -(lamb[0] * v[58]);
	v[90] = v[80] + v[81];
	v[67] = v[37] * v[58];
	v[49] = (-alphap[1] + v[119])*v[37];
	v[51] = v[37] * (alphap[0] + v[50]);
	v[52] = 1e0 + (v[115] + v[117])*v[122];
	v[57] = -(di[0] * v[49]) - di[1] * v[51] - di[2] * v[52];
	v[101] = lamb[1] * v[57];
	v[89] = v[100] + v[101];
	v[86] = -(lamb[0] * v[57]);
	v[88] = v[85] + v[86];
	v[61] = v[37] * v[57];
	v[53] = alphap[2] * v[122];
	v[92] = -(lamb[2] * v[53]);
	v[63] = -(v[53] * v[57]);
	v[54] = -(alphap[1] * v[122]);
	v[106] = lamb[1] * v[54];
	v[82] = -(v[53] * v[88]) + v[37] * v[89] - v[54] * v[90];
	v[71] = -(v[54] * v[58]);
	v[55] = alphap[0] * v[122];
	v[107] = -(lamb[0] * v[55]);
	v[102] = -(v[37] * v[88]) - v[53] * v[89] + v[55] * v[90];
	v[72] = -(v[55] * v[64]);
	v[74] = v[109] + lamb[0] * v[53];
	v[75] = -(lamb[0] * v[54]) + v[94];
	v[78] = -(v[58] * v[74]) - v[57] * v[75];
	v[76] = v[106] + v[92];
	v[83] = v[64] * v[75] + v[58] * v[76];
	v[79] = -(v[64] * v[74]) + v[57] * v[76];
	v[84] = v[53] * v[78] - v[55] * v[83] + v[37] * (v[79] - 0.5e0*(alphap[1] * v[82] + v[90]));
	v[87] = v[54] * v[78] + v[55] * v[79] + v[37] * (v[83] - 0.5e0*(alphap[2] * v[82] - v[85] - v[86]));
	v[91] = v[104] + lamb[1] * v[53];
	v[93] = v[107] + v[92];
	v[97] = -(v[58] * v[91]) - v[57] * v[93];
	v[95] = lamb[1] * v[55] - v[94];
	v[99] = -(v[64] * v[91]) + v[57] * v[95];
	v[98] = v[64] * v[93] + v[58] * v[95];
	v[103] = v[54] * v[97] + v[37] * (-0.5e0*(-(alphap[2] * v[102]) + v[89]) + v[98]) + v[55] * v[99];
	v[105] = -v[104] - lamb[2] * v[54];
	v[108] = v[106] + v[107];
	v[110] = -v[109] + lamb[2] * v[55];
	M2[0][0] = v[63] + v[71];
	M2[0][1] = -(v[55] * v[58]) + v[61];
	M2[0][2] = v[55] * v[57] + v[67];
	M2[1][0] = -v[61] - v[54] * v[64];
	M2[1][1] = v[63] + v[72];
	M2[1][2] = -(v[54] * v[57]) + v[69];
	M2[2][0] = v[53] * v[64] - v[67];
	M2[2][1] = -(v[53] * v[58]) - v[69];
	M2[2][2] = v[71] + v[72];
	M3[0][0] = v[37] * v[78] - v[53] * v[79] + v[55] * v[82] - v[54] * v[83];
	M3[0][1] = v[84];
	M3[0][2] = v[87];
	M3[1][0] = v[84];
	M3[1][1] = v[102] * v[54] + v[53] * v[97] - v[55] * v[98] + v[37] * v[99];
	M3[1][2] = v[103];
	M3[2][0] = v[87];
	M3[2][1] = v[103];
	M3[2][2] = v[54] * (-(v[105] * v[57]) - v[108] * v[58]) + v[37] * (v[110] * v[58] + v[105] * v[64]) + v[55] * (v[110] * v[57]
		- v[108] * v[64]) + v[53] * (v[55] * v[88] + v[54] * v[89] + v[37] * v[90]);
	Qp[0][0] = v[40];
	Qp[0][1] = v[41];
	Qp[0][2] = v[42];
	Qp[1][0] = v[44];
	Qp[1][1] = v[46];
	Qp[1][2] = v[47];
	Qp[2][0] = v[49];
	Qp[2][1] = v[51];
	Qp[2][2] = v[52];
};