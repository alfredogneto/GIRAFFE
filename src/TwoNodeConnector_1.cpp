#include "TwoNodeConnector_1.h"

#include "Node.h"
#include "CoordinateSystem.h"
#include "Dynamic.h"
#include"Database.h"
//Variaveis globais
extern
Database db;


TwoNodeConnector_1::TwoNodeConnector_1()
{
	strain_energy = 0.0;
	kinetic_energy = 0.0;
	potential_gravitational_energy = 0.0;

	VTK_type = 3;
	nDOFs = 12;
	material = 0;
	section = 0;
	n_nodes = 2;
	number = 0;
	nodes = new int[n_nodes];
	VTK_nodes = new int[n_nodes];
	VTK_nodes[0] = 0;
	VTK_nodes[1] = 1;

	DOFs = new int *[n_nodes];
	for (int i = 0; i < n_nodes; i++)
		DOFs[i] = new int[db.number_GLs_node];

	type_name = new char[20];//Nome do tipo do elemento
	sprintf(type_name, "TwoNodeConnector_1");

	//Rotina para ativar os GLS de cada nó do elemento
	for (int i = 0; i < n_nodes; i++)
	{
		for (int j = 0; j < db.number_GLs_node; j++)
		{
			DOFs[i][j] = 0;
		}
		DOFs[i][0] = 1;
		DOFs[i][1] = 1;
		DOFs[i][2] = 1;
		DOFs[i][3] = 1;
		DOFs[i][4] = 1;
		DOFs[i][5] = 1;
	}

	//Variaveis do elemento
	K = new Matrix(12, 12);													
	C = new Matrix(12, 12);	
	Q = new Matrix(12, 12);

	K_m = new Matrix(6, 6);													
	C_m = new Matrix(6, 6);
	Q_m = new Matrix(3, 3);
	
	f = new Matrix(12);
	fd = new Matrix(12);
	Kt = new Matrix(12, 12);
	Ct = new Matrix(12, 12);
}


TwoNodeConnector_1::~TwoNodeConnector_1()
{
	delete[] nodes;
	delete[] VTK_nodes;
	if (DOFs != NULL)
	{
		for (int i = 0; i < n_nodes; i++)
			delete[] DOFs[i];
		delete[] DOFs;
	}
	delete[]type_name;

	delete K;
	delete C;
	delete Q;

	delete K_m;
	delete C_m;
	delete Q_m;

	delete f;
	delete fd;
	delete Kt;
	delete Ct;
}

bool TwoNodeConnector_1::Read(FILE *f)
{
	char s[1000];
	fscanf(f, "%s", s);
	number = atoi(s);

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
		fscanf(f, "%s", s);
		nodes[0] = atoi(s);
		fscanf(f, "%s", s);
		nodes[1] = atoi(s);
	}
	else
		return false;

	//Stiffness data
	fscanf(f, "%s", s);
	if (!strcmp(s, "StiffnessData"))
	{
		for (int i = 0; i < 6; i++)
		{
			for (int j = 0; j < 6; j++)
			{
				fscanf(f, "%s", s);
				(*K_m)(i, j) = atof(s);
			}
		}
	}
	else
		return false;

	//Damping data
	fscanf(f, "%s", s);
	if (!strcmp(s, "DampingData"))
	{
		for (int i = 0; i < 6; i++)
		{
			for (int j = 0; j < 6; j++)
			{
				fscanf(f, "%s", s);
				(*C_m)(i, j) = atof(s);
			}
		}
	}
	else
		return false;

	return true;
}

//Checa inconsistências no elemento para evitar erros de execução
bool TwoNodeConnector_1::Check()
{
	for (int i = 0; i < n_nodes; i++)
	{
		if (nodes[i] > db.number_nodes)
			return false;
	}
	if (cs > db.number_CS)
		return false;
	return true;
}

void TwoNodeConnector_1::Write(FILE *f)
{
	fprintf(f, "TwoNodeConnector_1\t%d\tCS\t%d\tNodes\t%d\t%d\n", number, cs, nodes[0], nodes[1]);
	fprintf(f, "StiffnessData\n");
	for (int i = 0; i < 6; i++)
	{
		for (int j = 0; j < 6; j++)
		{
			fprintf(f, "%.6e\t", (*K_m)(i, j));
		}
		fprintf(f, "\n");
	}
	fprintf(f, "DampingData\n");
	for (int i = 0; i < 6; i++)
	{
		for (int j = 0; j < 6; j++)
		{
			fprintf(f, "%.6e\t", (*C_m)(i, j));
		}
		fprintf(f, "\n");
	}
}
//Escreve arquivo de resultados
void TwoNodeConnector_1::WriteResults(FILE *f)
{
	//DOES NOTHING
}

//Escreve no monitor do elemento//Escreve no monitor do elemento
void TwoNodeConnector_1::WriteMonitor(FILE *f, bool first_record, double time)
{
	//DOES NOTHING
}

void TwoNodeConnector_1::WriteVTK_XMLBase(std::vector<float> *float_vector)
{
	//DOES NOTHING
}

void TwoNodeConnector_1::WriteVTK_XMLRender(FILE *f)
{
	//TODO
}

//Monta carregamentos associados ao elemento
void TwoNodeConnector_1::MountElementLoads()
{
	//DOES NOTHING
}

//Monta elementos
void TwoNodeConnector_1::Mount()
{
	//Posições atuais dos GL do elemento
	Matrix xipp(12);
	Matrix uA(3);
	Matrix uB(3);
	Matrix alphaA(3);
	Matrix alphaB(3);
	//Deslocamentos
	for (int i = 0; i < 3; i++)
	{
		uA(i, 0) = db.nodes[nodes[0] - 1]->copy_coordinates[i] + db.nodes[nodes[0] - 1]->displacements[i] - db.nodes[nodes[0] - 1]->ref_coordinates[i];
		uB(i, 0) = db.nodes[nodes[1] - 1]->copy_coordinates[i] + db.nodes[nodes[1] - 1]->displacements[i] - db.nodes[nodes[1] - 1]->ref_coordinates[i];
	}
	Matrix alpha_1(3);
	Matrix alpha_2(3);
	//Rotações - nó A
	alpha_1(0, 0) = db.nodes[nodes[0] - 1]->copy_coordinates[3];
	alpha_1(1, 0) = db.nodes[nodes[0] - 1]->copy_coordinates[4];
	alpha_1(2, 0) = db.nodes[nodes[0] - 1]->copy_coordinates[5];
	alpha_2(0, 0) = db.nodes[nodes[0] - 1]->displacements[3];
	alpha_2(1, 0) = db.nodes[nodes[0] - 1]->displacements[4];
	alpha_2(2, 0) = db.nodes[nodes[0] - 1]->displacements[5];
	//Fórmula de Rodrigues
	alphaA = 4.0 / (4.0 - dot(alpha_2, alpha_1))*(alpha_2 + alpha_1 + 0.5*cross(alpha_2, alpha_1));	//Vetor rotação de Rodrigues
	//////////////////////Salvando rotações de Euler////////////////////
	//double theta_escalar = 2.0*atan(norm(alphaA) / 2.0);
	//if (norm(alphaA) != 0.0)
	//	alphaA = (theta_escalar / norm(alphaA))*(alphaA);
	//else
	//	zeros(&alphaA);
	//Rotações - nó B
	alpha_1(0, 0) = db.nodes[nodes[1] - 1]->copy_coordinates[3];
	alpha_1(1, 0) = db.nodes[nodes[1] - 1]->copy_coordinates[4];
	alpha_1(2, 0) = db.nodes[nodes[1] - 1]->copy_coordinates[5];
	alpha_2(0, 0) = db.nodes[nodes[1] - 1]->displacements[3];
	alpha_2(1, 0) = db.nodes[nodes[1] - 1]->displacements[4];
	alpha_2(2, 0) = db.nodes[nodes[1] - 1]->displacements[5];
	//Fórmula de Rodrigues
	alphaB = 4.0 / (4.0 - dot(alpha_2, alpha_1))*(alpha_2 + alpha_1 + 0.5*cross(alpha_2, alpha_1));	//Vetor rotação de Rodrigues
	//////////////////////Salvando rotações de Euler////////////////////
	//theta_escalar = 2.0*atan(norm(alphaB) / 2.0);
	//if (norm(alphaB) != 0.0)
	//	alphaB = (theta_escalar / norm(alphaB))*(alphaB);
	//else
	//	zeros(&alphaB);
	//Compondo o vetor xipp
	for (int i = 0; i < 3; i++)
	{
		xipp(i, 0) = uA(i, 0);
		xipp(i + 3, 0) = alphaA(i, 0);
		xipp(i + 6, 0) = uB(i, 0);
		xipp(i + 9, 0) = alphaB(i, 0);
	}

	//Forças elasticas
	*f = (*K)*xipp;
	//Operador tangente (parcela da rigidez)
	*Kt = (*K);
	//f->print();
	//Kt->print();
}
//Monta matriz de transformação de coordenadas
void TwoNodeConnector_1::TransformMatrix()
{
	//DOES NOTHING
}
//Preenche a contribuição do elemento nas matrizes globais
void TwoNodeConnector_1::MountGlobal()
{
	int GL_global_1 = 0;
	int GL_global_2 = 0;
	double anterior = 0;
	for (int i = 0; i < 12; i++)
	{
		//Nó 1
		if (i < 6)
			GL_global_1 = db.nodes[nodes[0] - 1]->GLs[i];
		//Nó 2
		else
			GL_global_1 = db.nodes[nodes[1] - 1]->GLs[i - 6];

		//Caso o grau de liberdade seja livre:
		if (GL_global_1 > 0)
		{
			anterior = db.global_P_A(GL_global_1 - 1, 0);
			db.global_P_A(GL_global_1 - 1, 0) = anterior + (*f)(i, 0);
			anterior = db.global_I_A(GL_global_1 - 1, 0);
			db.global_I_A(GL_global_1 - 1, 0) = anterior + (*f)(i, 0);
		}
		else
		{
			anterior = db.global_P_B(-GL_global_1 - 1, 0);
			db.global_P_B(-GL_global_1 - 1, 0) = anterior + (*f)(i, 0);
		}
		for (int j = 0; j < 12; j++)
		{
			//Nó 1
			if (j < 6)
				GL_global_2 = db.nodes[nodes[0] - 1]->GLs[j];
			//Nó 2
			else
				GL_global_2 = db.nodes[nodes[1] - 1]->GLs[j - 6];

			//Caso os graus de liberdade sejam ambos livres (Matriz Kaa)
			if (GL_global_1 > 0 && GL_global_2 > 0)
				db.global_stiffness_AA.setValue(GL_global_1 - 1, GL_global_2 - 1, (*Kt)(i, j));
			//Caso os graus de liberdade sejam ambos fixos (Matriz Kbb)
			if (GL_global_1 < 0 && GL_global_2 < 0)
				db.global_stiffness_BB.setValue(-GL_global_1 - 1, -GL_global_2 - 1, (*Kt)(i, j));
			//Caso os graus de liberdade sejam livre e fixo (Matriz Kab)
			if (GL_global_1 > 0 && GL_global_2 < 0)
				db.global_stiffness_AB.setValue(GL_global_1 - 1, -GL_global_2 - 1, (*Kt)(i, j));
			//Caso os graus de liberdade sejam fixo e livre (Matriz Kba)
			if (GL_global_1 < 0 && GL_global_2 > 0)
				db.global_stiffness_BA.setValue(-GL_global_1 - 1, GL_global_2 - 1, (*Kt)(i, j));
		}
	}
}
//Salva variaveis nos pontos de Gauss uteis para descrição lagrangiana atualizada
void TwoNodeConnector_1::SaveLagrange()
{
	//DOES NOTHING
}
//Pre-calculo de variaveis que e feito uma unica vez no inicio
void TwoNodeConnector_1::PreCalc()
{
	//Transformação de coordenadas 3x3
	(*Q_m)(0, 0) = (*db.CS[cs - 1]->E1)(0, 0);
	(*Q_m)(0, 1) = (*db.CS[cs - 1]->E1)(1, 0);
	(*Q_m)(0, 2) = (*db.CS[cs - 1]->E1)(2, 0);
	(*Q_m)(1, 0) = (*db.CS[cs - 1]->E2)(0, 0);
	(*Q_m)(1, 1) = (*db.CS[cs - 1]->E2)(1, 0);
	(*Q_m)(1, 2) = (*db.CS[cs - 1]->E2)(2, 0);
	(*Q_m)(2, 0) = (*db.CS[cs - 1]->E3)(0, 0);
	(*Q_m)(2, 1) = (*db.CS[cs - 1]->E3)(1, 0);
	(*Q_m)(2, 2) = (*db.CS[cs - 1]->E3)(2, 0);

	//Preenchendo a transformação de coordenadas 12x12
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			(*Q)(i, j) = (*Q_m)(i, j);
			(*Q)(i + 3, j + 3) = (*Q_m)(i, j);
			(*Q)(i + 6, j + 6) = (*Q_m)(i, j);
			(*Q)(i + 9, j + 9) = (*Q_m)(i, j);
		}
	}
	//Composição da matriz de rigidez e amortecimento (antes da transformação de coordenadas)
	for (int i = 0; i < 6; i++)
	{
		for (int j = 0; j < 6; j++)
		{
			(*K)(i, j) = (*K_m)(i, j);
			(*K)(i + 6, j + 6) = (*K_m)(i, j);
			(*K)(i + 6, j) = -(*K_m)(i, j);
			(*K)(i, j + 6) = -(*K_m)(i, j);

			(*C)(i, j) = (*C_m)(i, j);
			(*C)(i + 6, j + 6) = (*C_m)(i, j);
			(*C)(i + 6, j) = -(*C_m)(i, j);
			(*C)(i, j + 6) = -(*C_m)(i, j);
		}
	}
	//K->print();
	//C->print();
	//Transformação de coordenadas
	*K = transp(*Q)*(*K)*(*Q);
	*C = transp(*Q)*(*C)*(*Q);
	//K->print();
	//C->print();
}

//Monta a matriz de massa
void TwoNodeConnector_1::MountMass()
{
	//DOES NOTHING
}
//Monta a matriz de massa
void TwoNodeConnector_1::MountMassModal()
{
	//DOES NOTHING
}

//Monta a matriz de amortecimento para realização da analise modal
void TwoNodeConnector_1::MountDampingModal()
{
	//DOES NOTHING
}

//Monta a matriz de amortecimento
void TwoNodeConnector_1::MountDamping(bool update_rayleigh)
{
	if (typeid(*db.solution[db.current_solution_number - 1]) == typeid(Dynamic))
	{
		Dynamic* ptr_sol = static_cast<Dynamic*>(db.solution[db.current_solution_number - 1]);
		//Amortecimento de Rayleigh + Matriz de amortecimento definida
		*Ct = ptr_sol->beta*(*K) + (*C);
		//Velocidade nos GL do elemento
		Matrix vipp(12);
		for (int i = 0; i < 6; i++)
		{
			vipp(i, 0) = db.nodes[nodes[0] - 1]->vel[i];
			vipp(i + 6, 0) = db.nodes[nodes[1] - 1]->vel[i];
		}
		*fd = (*Ct) * vipp;
		//Damping para Newmark para compor o operador tangente
		*Ct = ptr_sol->a4*(*Ct);
	}
}

//Montagens - Newmark
void TwoNodeConnector_1::MountDyn()
{
	*f = (*f) + (*fd);
	*Kt = (*Kt) + (*Ct);
}

//Montagens para analise modal - inserção da matriz de massa e amortecimento na matriz de rigidez para posterior montagem global
void TwoNodeConnector_1::MountDynModal()
{
	//No influences are considered (no mass, neither damping)
	zeros(Kt);
	zeros(Ct);
}

//Zera algumas matrizes utilizadas nos calculos
void TwoNodeConnector_1::Zeros()
{
	zeros(f);
	zeros(Kt);
	zeros(Ct);
	kinetic_energy = 0.0;
	strain_energy = 0.0;
	potential_gravitational_energy = 0.0;
}