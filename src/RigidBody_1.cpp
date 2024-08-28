#include "RigidBody_1.h"

#include "RigidBodyData.h"
#include "Node.h"
#include "CoordinateSystem.h"
#include "Dynamic.h"
#include "Environment.h"
#include"Database.h"
//Variáveis globais
extern
Database db;

RigidBody_1::RigidBody_1()
{
	strain_energy = 0.0;
	kinetic_energy = 0.0;
	potential_gravitational_energy = 0.0;

	nDOFs = 6;
	material = 0;
	section = 0;
	n_nodes = 1;
	number = 0;		//ID
	VTK_nodes = new int[n_nodes];
	VTK_nodes[0] = 0;
	VTK_type = 1;
	nodes = new int[n_nodes];
	DOFs = new int *[n_nodes];
	for (int i = 0; i < n_nodes;i++)
		DOFs[i] = new int[db.number_GLs_node];
	type_name = new char[20];//Nome do tipo do elemento
	sprintf(type_name, "RigidBody_1");
	//Rotina para ativar os GLS de cada nó do elemento
	for (int i = 0; i < n_nodes; i++)
	{
		for (int j = 0; j < db.number_GLs_node; j++)
			DOFs[i][j] = 0;
		DOFs[i][0] = 1;
		DOFs[i][1] = 1;
		DOFs[i][2] = 1;
		DOFs[i][3] = 1;
		DOFs[i][4] = 1;
		DOFs[i][5] = 1;
	}

	DdT = new double*[6];
	for (int i = 0; i < 6; i++)
		DdT[i] = new double[6];
	dT = new double[6];

	Jr = new double*[3];
	for (int i = 0; i < 3; i++)
		Jr[i] = new double[3];
	br = new double[3];
	
	T1 = new double;
	*T1 = 0;
	T2 = new double;
	*T2 = 0;
	T3 = new double;
	*T3 = 0;
	T = new double;
	*T = 0;

	magL = new double;
	*magL = 0;
	magHG = new double;
	*magHG = 0;

	L = new double[3];
	HG = new double[3];

	Ddfield = new double*[6];
	for (int i = 0; i < 6; i++)
		Ddfield[i] = new double[6];
	dfield = new double[6];

	for (int i = 0; i < 6; i++)
	{
		dfield[i] = 0.0;
		dT[i] = 0.0;
		for (int j = 0; j < 6; j++)
		{
			Ddfield[i][j] = 0.0;
			DdT[i][j] = 0.0;
		}
	}


	for (int i = 0; i < 3; i++)
	{
		alphai[i] = 0;
		ud[i] = 0;
		alphad[i] = 0;
		dui[i] = 0;
		ddui[i] = 0;
		omegai[i] = 0;
		domegai[i] = 0;
		L[i] = 0;
		HG[i] = 0;
	}


}

RigidBody_1::~RigidBody_1()
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
	
	for (int i = 0; i < 6; i++)
		delete[]DdT[i];
	delete[] DdT;
	delete[] dT;

	for (int i = 0; i < 3; i++)
		delete[]Jr[i];
	delete[]Jr;
	delete[]br;

	delete T1;
	delete T2;
	delete T3;
	delete T;
	delete[]L;
	delete magL;
	delete[]HG;
	delete magHG;

	for (int i = 0; i < 6; i++)
		delete[]Ddfield[i];
	delete[] Ddfield;
	delete[] dfield;
	

}

bool RigidBody_1::Read(FILE *f)
{
	char s[1000];
	fscanf(f, "%s", s);
	number = atoi(s);

	fscanf(f, "%s", s);
	if (!strcmp(s, "RigidBodyData"))
	{
		fscanf(f, "%s", s);
		RB_data_ID = atoi(s);
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
	if (!strcmp(s, "Node"))
	{
		fscanf(f, "%s", s);
		nodes[0] = atoi(s);
	}
	else
		return false;

	return true;
}

//Checa inconsistências no elemento para evitar erros de execução
bool RigidBody_1::Check()
{
	for (int i = 0; i < n_nodes; i++)
	{
		if (nodes[i] > db.number_nodes)
			return false;
	}
	if (RB_data_ID > db.number_RB_data)
		return false;
	return true;
}

void RigidBody_1::Write(FILE *f)
{
	fprintf(f, "RigidBody_1\t%d\tRigidBodyData\t%d\tCS\t%d\tNode\t%d\n", number, RB_data_ID, cs, nodes[0]);
}
//Escreve arquivo de resultados
void RigidBody_1::WriteResults(FILE *f)
{
	
}

//Escreve no monitor do elemento//Escreve no monitor do elemento
void RigidBody_1::WriteMonitor(FILE *f, bool first_record, double time)
{
	//Cabeçalho
	if (first_record == true)
		fprintf(f, "TIME\tT1\tT2\tT3\tT\tLx\tLy\tLz\tL\tHGx\tHGy\tHGz\tHG\n");

	//Informações a serem salvas
	fprintf(f, "%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\n",
		time, *T1, *T2, *T3, *T, L[0], L[1], L[2],*magL, HG[0], HG[1], HG[2], *magHG);
}

void RigidBody_1::WriteVTK_XMLBase(std::vector<float> *float_vector)
{
	//Imprime os resultados do elemento
	int res_element = 0;
	//Imprime valores nulos para resultados que não fazem sentido para esse tipo de elemento
	for (int i = res_element; i < db.n_element_results; i++)
		float_vector->push_back(0.0);
}

void RigidBody_1::WriteVTK_XMLRender(FILE *f)
{
	db.RB_data[RB_data_ID - 1]->WriteVTK_XMLRender(f,nodes[0],cs);
}

//Monta carregamentos associados ao elemento
void RigidBody_1::MountElementLoads()
{
	MountFieldLoading();
}

//Monta elementos
void RigidBody_1::Mount()
{
	Zeros();	//Zera matrizes
	//Não há contribuição para rigidez nesse elemento
}
//Monta matriz de transformação de coordenadas
void RigidBody_1::TransformMatrix()
{
	
}
//Preenche a contribuição do elemento nas matrizes globais
void RigidBody_1::MountGlobal()
{
	//Variáveis temporárias para salvar a indexação global dos graus de liberdade a serem setados na matriz de rigidez global
	int GL_global_1 = 0;
	int GL_global_2 = 0;
	double anterior = 0;
	for (int i = 0; i<6; i++)
	{
		//////////////MONTAGEM DO VETOR DE ESFORÇOS DESBALANCEADOS//////////////////
		//Toma do vetor de GL globais, a indexação de cada grau de liberdade global
		GL_global_1 = db.nodes[nodes[0] - 1]->GLs[i];

		//Caso o grau de liberdade seja livre:
		if (GL_global_1 > 0)
		{
			anterior = db.global_P_A(GL_global_1 - 1, 0);
			db.global_P_A(GL_global_1 - 1, 0) = anterior + dT[i] - dfield[i];
			anterior = db.global_I_A(GL_global_1 - 1, 0);				
			db.global_I_A(GL_global_1 - 1, 0) = anterior + dT[i] - dfield[i];
		}
		else
		{
			anterior = db.global_P_B(-GL_global_1 - 1, 0);
			db.global_P_B(-GL_global_1 - 1, 0) = anterior + dT[i] - dfield[i];
		}
		for (int j = 0; j<6; j++)
		{
			//////////////////////MONTAGEM DA MATRIZ DE RIGIDEZ/////////////////////////
			//Toma do vetor de GL globais, a indexação de cada grau de liberdade global
			GL_global_2 = db.nodes[nodes[0] - 1]->GLs[j];
			//Caso os graus de liberdade sejam ambos livres (Matriz Kaa)
			if (GL_global_1 > 0 && GL_global_2 > 0)
				db.global_stiffness_AA.setValue(GL_global_1 - 1, GL_global_2 - 1, (DdT[i][j] - Ddfield[i][j]));
			//Caso os graus de liberdade sejam ambos fixos (Matriz Kbb)
			if (GL_global_1 < 0 && GL_global_2 < 0)
				db.global_stiffness_BB.setValue(-GL_global_1 - 1, -GL_global_2 - 1, (DdT[i][j] - Ddfield[i][j]));
			//Caso os graus de liberdade sejam livre e fixo (Matriz Kab)
			if (GL_global_1 > 0 && GL_global_2 < 0)
				db.global_stiffness_AB.setValue(GL_global_1 - 1, -GL_global_2 - 1, (DdT[i][j] - Ddfield[i][j]));
			//Caso os graus de liberdade sejam fixo e livre (Matriz Kba)
			if (GL_global_1 < 0 && GL_global_2 > 0)
				db.global_stiffness_BA.setValue(-GL_global_1 - 1, GL_global_2 - 1, (DdT[i][j] - Ddfield[i][j]));
		}
	}
}
//Salva variáveis nos pontos de Gauss úteis para descrição lagrangiana atualizada
void RigidBody_1::SaveLagrange()
{
	
}
//Pré-cálculo de variáveis que é feito uma única vez no início
void RigidBody_1::PreCalc()
{
	//Matriz para transformação de coordenadas  local-global
	Matrix Q(3, 3);
	Q(0, 0) = (*db.CS[cs - 1]->E1)(0, 0);
	Q(1, 0) = (*db.CS[cs - 1]->E1)(1, 0);
	Q(2, 0) = (*db.CS[cs - 1]->E1)(2, 0);
	Q(0, 1) = (*db.CS[cs - 1]->E2)(0, 0);
	Q(1, 1) = (*db.CS[cs - 1]->E2)(1, 0);
	Q(2, 1) = (*db.CS[cs - 1]->E2)(2, 0);
	Q(0, 2) = (*db.CS[cs - 1]->E3)(0, 0);
	Q(1, 2) = (*db.CS[cs - 1]->E3)(1, 0);
	Q(2, 2) = (*db.CS[cs - 1]->E3)(2, 0);

	//Calculando valores do br. Feito apenas uma vez no início. 
	//A origem do CAD coincide com o Polo, então G-O = G
	//Como G é dado no sistema do CAD, então G precisa ser transformado
	Matrix local_b = *db.RB_data[RB_data_ID - 1]->G;
	Matrix global_b = Q*local_b;
	br[0] = global_b(0, 0);
	br[1] = global_b(1, 0);
	br[2] = global_b(2, 0);

	//Atualizando J_G para o sistema global e copiando valores para Jr. Feito apenas uma vez no início
	Matrix global_J = Q*(*db.RB_data[RB_data_ID - 1]->J_G)*transp(Q);
	global_J.MatrixToPtr(Jr, 3);

	//Aplicando Teorema de Steiner
	//Momentos de Inércia
	Jr[0][0] = Jr[0][0] + (db.RB_data[RB_data_ID - 1]->mass)*((br[1] * br[1]) + (br[2] * br[2]));
	Jr[1][1] = Jr[1][1] + (db.RB_data[RB_data_ID - 1]->mass)*((br[0] * br[0]) + (br[2] * br[2]));
	Jr[2][2] = Jr[2][2] + (db.RB_data[RB_data_ID - 1]->mass)*((br[0] * br[0]) + (br[1] * br[1]));
	//Produtos de Inércia
	Jr[0][1] = Jr[0][1] - (db.RB_data[RB_data_ID - 1]->mass)*(br[0] * br[1]);
	Jr[1][0] = Jr[0][1];
	Jr[0][2] = Jr[0][2] - (db.RB_data[RB_data_ID - 1]->mass)*(br[0] * br[2]);
	Jr[2][0] = Jr[0][2];
	Jr[1][2] = Jr[1][2] - (db.RB_data[RB_data_ID - 1]->mass)*(br[1] * br[2]);
	Jr[2][1] = Jr[1][2];
}

//Monta a matriz de massa
void RigidBody_1::MountMass()
{
	InertialContributions();
}
//Monta a matriz de massa
void RigidBody_1::MountMassModal()
{
	//Ponteiro para o valor da massa
	double* m = &db.RB_data[RB_data_ID - 1]->mass;
	//Valores de variáveis cinemáticas - atribuição
	for (int i = 0; i < 3; i++)
		alphai[i] = db.nodes[nodes[0] - 1]->copy_coordinates[i + 3];
	double Mr[3][3];
	for (int i = 0; i < 3;i++)
	for (int j = 0; j < 3; j++)
		Mr[i][j] = 0.0;
	Mr[0][0] = *m;
	Mr[1][1] = *m;
	Mr[2][2] = *m;

	//AceGen
	v[140] = Power(alphai[2], 2);
	v[139] = 0.5e0*alphai[0] * alphai[2];
	v[138] = 0.5e0*alphai[1];
	v[137] = Power(alphai[1], 2);
	v[141] = v[137] + v[140];
	v[136] = alphai[0] * v[138];
	v[135] = Power(alphai[0], 2);
	v[83] = alphai[2] * v[138];
	v[70] = 4e0 / (4e0 + v[135] + v[141]);
	v[142] = -0.5e0*v[70];
	v[73] = 1e0 + v[141] * v[142];
	v[74] = (-alphai[2] + v[136])*v[70];
	v[75] = (alphai[1] + v[139])*v[70];
	v[107] = Jr[0][2] * v[73] + Jr[1][2] * v[74] + Jr[2][2] * v[75];
	v[106] = Jr[0][1] * v[73] + Jr[1][1] * v[74] + Jr[2][1] * v[75];
	v[105] = Jr[0][0] * v[73] + Jr[1][0] * v[74] + Jr[2][0] * v[75];
	v[89] = Mr[0][2] * v[73] + Mr[1][2] * v[74] + Mr[2][2] * v[75];
	v[88] = Mr[0][1] * v[73] + Mr[1][1] * v[74] + Mr[2][1] * v[75];
	v[87] = Mr[0][0] * v[73] + Mr[1][0] * v[74] + Mr[2][0] * v[75];
	v[77] = (alphai[2] + v[136])*v[70];
	v[79] = 1e0 + (v[135] + v[140])*v[142];
	v[80] = v[70] * (-alphai[0] + v[83]);
	v[113] = Jr[0][2] * v[77] + Jr[1][2] * v[79] + Jr[2][2] * v[80];
	v[112] = Jr[0][1] * v[77] + Jr[1][1] * v[79] + Jr[2][1] * v[80];
	v[111] = Jr[0][0] * v[77] + Jr[1][0] * v[79] + Jr[2][0] * v[80];
	v[95] = Mr[0][2] * v[77] + Mr[1][2] * v[79] + Mr[2][2] * v[80];
	v[94] = Mr[0][1] * v[77] + Mr[1][1] * v[79] + Mr[2][1] * v[80];
	v[93] = Mr[0][0] * v[77] + Mr[1][0] * v[79] + Mr[2][0] * v[80];
	v[82] = (-alphai[1] + v[139])*v[70];
	v[84] = v[70] * (alphai[0] + v[83]);
	v[85] = 1e0 + (v[135] + v[137])*v[142];
	v[119] = Jr[0][2] * v[82] + Jr[1][2] * v[84] + Jr[2][2] * v[85];
	v[118] = Jr[0][1] * v[82] + Jr[1][1] * v[84] + Jr[2][1] * v[85];
	v[117] = Jr[0][0] * v[82] + Jr[1][0] * v[84] + Jr[2][0] * v[85];
	v[101] = Mr[0][2] * v[82] + Mr[1][2] * v[84] + Mr[2][2] * v[85];
	v[100] = Mr[0][1] * v[82] + Mr[1][1] * v[84] + Mr[2][1] * v[85];
	v[99] = Mr[0][0] * v[82] + Mr[1][0] * v[84] + Mr[2][0] * v[85];
	v[86] = v[73] * v[87] + v[74] * v[88] + v[75] * v[89];
	v[90] = v[77] * v[87] + v[79] * v[88] + v[80] * v[89];
	v[91] = v[82] * v[87] + v[84] * v[88] + v[85] * v[89];
	v[92] = v[73] * v[93] + v[74] * v[94] + v[75] * v[95];
	v[96] = v[77] * v[93] + v[79] * v[94] + v[80] * v[95];
	v[97] = v[82] * v[93] + v[84] * v[94] + v[85] * v[95];
	v[98] = v[100] * v[74] + v[101] * v[75] + v[73] * v[99];
	v[102] = v[100] * v[79] + v[101] * v[80] + v[77] * v[99];
	v[103] = v[100] * v[84] + v[101] * v[85] + v[82] * v[99];
	v[122] = br[0] * v[73] + br[1] * v[74] + br[2] * v[75];
	v[123] = br[0] * v[77] + br[1] * v[79] + br[2] * v[80];
	v[124] = br[0] * v[82] + br[1] * v[84] + br[2] * v[85];
	v[125] = -(v[124] * v[90]) + v[123] * v[91];
	v[126] = v[124] * v[86] - v[122] * v[91];
	v[127] = -(v[123] * v[86]) + v[122] * v[90];
	v[128] = -(v[124] * v[96]) + v[123] * v[97];
	v[129] = v[124] * v[92] - v[122] * v[97];
	v[130] = -(v[123] * v[92]) + v[122] * v[96];
	v[131] = v[103] * v[123] - v[102] * v[124];
	v[132] = -(v[103] * v[122]) + v[124] * v[98];
	v[133] = v[102] * v[122] - v[123] * v[98];
	DdT[0][0] = v[86];
	DdT[0][1] = v[90];
	DdT[0][2] = v[91];
	DdT[0][3] = v[125];
	DdT[0][4] = v[126];
	DdT[0][5] = v[127];
	DdT[1][0] = v[92];
	DdT[1][1] = v[96];
	DdT[1][2] = v[97];
	DdT[1][3] = v[128];
	DdT[1][4] = v[129];
	DdT[1][5] = v[130];
	DdT[2][0] = v[98];
	DdT[2][1] = v[102];
	DdT[2][2] = v[103];
	DdT[2][3] = v[131];
	DdT[2][4] = v[132];
	DdT[2][5] = v[133];
	DdT[3][0] = -v[125];
	DdT[3][1] = -v[126];
	DdT[3][2] = -v[127];
	DdT[3][3] = v[105] * v[73] + v[106] * v[74] + v[107] * v[75];
	DdT[3][4] = v[105] * v[77] + v[106] * v[79] + v[107] * v[80];
	DdT[3][5] = v[105] * v[82] + v[106] * v[84] + v[107] * v[85];
	DdT[4][0] = -v[128];
	DdT[4][1] = -v[129];
	DdT[4][2] = -v[130];
	DdT[4][3] = v[111] * v[73] + v[112] * v[74] + v[113] * v[75];
	DdT[4][4] = v[111] * v[77] + v[112] * v[79] + v[113] * v[80];
	DdT[4][5] = v[111] * v[82] + v[112] * v[84] + v[113] * v[85];
	DdT[5][0] = -v[131];
	DdT[5][1] = -v[132];
	DdT[5][2] = -v[133];
	DdT[5][3] = v[117] * v[73] + v[118] * v[74] + v[119] * v[75];
	DdT[5][4] = v[117] * v[77] + v[118] * v[79] + v[119] * v[80];
	DdT[5][5] = v[117] * v[82] + v[118] * v[84] + v[119] * v[85];
}

//Monta a matriz de amortecimento para realização da análise modal
void RigidBody_1::MountDampingModal()
{
	
}

//Monta a matriz de amortecimento
void RigidBody_1::MountDamping(bool update_rayleigh)
{
	
}

//Montagens - Newmark
void RigidBody_1::MountDyn()
{
	
}

//Montagens para análise modal - inserção da matriz de massa e amortecimento na matriz de rigidez para posterior montagem global
void RigidBody_1::MountDynModal()
{
	
}

//Zera algumas matrizes utilizadas nos cálculos
void RigidBody_1::Zeros()
{
	for (int i = 0; i < 6; i++)
	{
		dfield[i] = 0.0;
		dT[i] = 0.0;
		for (int j = 0; j < 6; j++)
		{
			Ddfield[i][j] = 0.0;
			DdT[i][j] = 0.0;
		}
	}
	kinetic_energy = 0.0;
	strain_energy = 0.0;
	potential_gravitational_energy = 0.0;
}

//Essa função escreve nas variáveis dT e DdT as contribuições para resíduo e para o operador tangente
void RigidBody_1::InertialContributions()
{
	
	Dynamic* ptr_sol = static_cast<Dynamic*>(db.solution[db.current_solution_number - 1]);
	//Ponteiros para parâmetros do Newmark
	double* a1 = &ptr_sol->a1;
	double* a2 = &ptr_sol->a2;
	double* a3 = &ptr_sol->a3;
	double* a4 = &ptr_sol->a4;
	double* a5 = &ptr_sol->a5;
	double* a6 = &ptr_sol->a6;
	//Ponteiro para o valor da massa
	double* m = &db.RB_data[RB_data_ID - 1]->mass;
	//Valores de variáveis cinemáticas - atribuição
	for (int i = 0; i < 3; i++)
	{
		alphai[i] = db.nodes[nodes[0] - 1]->copy_coordinates[i + 3];
		ud[i] = db.nodes[nodes[0] - 1]->displacements[i];
		alphad[i] = db.nodes[nodes[0] - 1]->displacements[i+3];
		dui[i] = db.nodes[nodes[0] - 1]->copy_vel[i];
		ddui[i] = db.nodes[nodes[0] - 1]->copy_accel[i];
		omegai[i] = db.nodes[nodes[0] - 1]->copy_vel[i + 3];
		domegai[i] = db.nodes[nodes[0] - 1]->copy_accel[i + 3];
	}
	
	//AceGen
	v[544] = (*a1)*(*m);
	v[543] = -((*a3)*ddui[2]) - (*a2)*dui[2] + (*a1)*ud[2];
	v[542] = -((*a3)*ddui[1]) - (*a2)*dui[1] + (*a1)*ud[1];
	v[541] = -((*a3)*ddui[0]) - (*a2)*dui[0] + (*a1)*ud[0];
	v[540] = (*a6)*ddui[2] + (*a5)*dui[2] + (*a4)*ud[2];
	v[539] = (*a6)*ddui[1] + (*a5)*dui[1] + (*a4)*ud[1];
	v[538] = (*a6)*ddui[0] + (*a5)*dui[0] + (*a4)*ud[0];
	v[531] = (*a4)*alphad[2] + (*a6)*domegai[2] + (*a5)*omegai[2];
	v[530] = (*a1)*alphad[2] - (*a3)*domegai[2] - (*a2)*omegai[2];
	v[529] = (*a4)*alphad[1] + (*a6)*domegai[1] + (*a5)*omegai[1];
	v[528] = (*a1)*alphad[1] - (*a3)*domegai[1] - (*a2)*omegai[1];
	v[527] = (*a4)*alphad[0] + (*a6)*domegai[0] + (*a5)*omegai[0];
	v[526] = (*a1)*alphad[0] - (*a3)*domegai[0] - (*a2)*omegai[0];
	v[524] = Power(alphad[2], 2);
	v[523] = 0.5e0*alphad[2];
	v[522] = 2e0*alphad[2];
	v[521] = Power(alphad[1], 2);
	v[520] = 0.5e0*alphad[1];
	v[519] = 2e0*alphad[1];
	v[518] = Power(alphad[0], 2);
	v[517] = 2e0*alphad[0];
	v[516] = 0.5e0*alphad[0];
	v[515] = Power(alphai[2], 2);
	v[514] = 0.5e0*alphai[0] * alphai[2];
	v[513] = 0.5e0*alphai[1];
	v[512] = Power(alphai[1], 2);
	v[535] = v[512] + v[515];
	v[511] = alphai[0] * v[513];
	v[510] = Power(alphai[0], 2);
	v[133] = alphai[2] * v[513];
	v[110] = alphad[1] * v[516];
	v[251] = -v[518] - v[521];
	v[229] = alphad[2] + v[110];
	v[219] = -alphad[2] + v[110];
	v[117] = alphad[2] * v[520];
	v[246] = alphad[0] + v[117];
	v[238] = -alphad[0] + v[117];
	v[115] = alphad[2] * v[516];
	v[242] = -alphad[1] + v[115];
	v[225] = alphad[1] + v[115];
	v[233] = -v[518] - v[524];
	v[214] = -v[521] - v[524];
	v[209] = 4e0 + v[518] + v[521] + v[524];
	v[211] = 1e0 / Power(v[209], 2);
	v[525] = -4e0*v[211];
	v[213] = v[522] * v[525];
	v[534] = 0.5e0*v[213];
	v[255] = v[251] * v[534];
	v[212] = v[519] * v[525];
	v[532] = 0.5e0*v[212];
	v[374] = -(v[212] * v[523]);
	v[235] = v[233] * v[532];
	v[210] = v[517] * v[525];
	v[533] = 0.5e0*v[210];
	v[376] = v[210] * v[520];
	v[373] = -(v[210] * v[523]);
	v[215] = v[214] * v[533];
	v[104] = 4e0 / v[209];
	v[377] = -0.5e0*v[104];
	v[379] = v[377] - v[210] * v[516];
	v[378] = -v[377] + v[212] * v[520];
	v[375] = v[377] - v[213] * v[523];
	v[253] = v[377] * v[519];
	v[254] = v[253] + v[251] * v[532];
	v[250] = v[377] * v[517];
	v[252] = v[250] + v[251] * v[533];
	v[247] = v[104] + v[210] * v[246];
	v[244] = -v[104] + v[212] * v[242];
	v[239] = -v[104] + v[210] * v[238];
	v[236] = v[377] * v[522];
	v[237] = v[236] + v[233] * v[534];
	v[234] = v[250] + v[233] * v[533];
	v[232] = v[104] + v[213] * v[229];
	v[227] = v[104] + v[212] * v[225];
	v[224] = -(v[104] * v[523]);
	v[248] = -v[224] + v[212] * v[246];
	v[243] = -v[224] + v[210] * v[242];
	v[240] = -v[224] + v[212] * v[238];
	v[226] = -v[224] + v[210] * v[225];
	v[223] = -v[104] + v[213] * v[219];
	v[221] = -(v[104] * v[516]);
	v[245] = -v[221] + v[213] * v[242];
	v[231] = -v[221] + v[212] * v[229];
	v[228] = -v[221] + v[213] * v[225];
	v[222] = v[212] * v[219] - v[221];
	v[218] = v[104] * v[520];
	v[249] = v[218] + v[213] * v[246];
	v[241] = v[218] + v[213] * v[238];
	v[230] = v[218] + v[210] * v[229];
	v[220] = v[218] + v[210] * v[219];
	v[217] = v[236] + v[214] * v[534];
	v[216] = v[253] + v[214] * v[532];
	v[107] = 1e0 - v[214] * v[377];
	v[310] = (*a1)*v[107] + v[215] * v[526] + v[220] * v[528] + v[226] * v[530];
	v[292] = (*a4)*v[107] + v[215] * v[527] + v[220] * v[529] + v[226] * v[531];
	v[108] = v[104] * v[219];
	v[311] = (*a1)*v[108] + v[216] * v[526] + v[222] * v[528] + v[227] * v[530];
	v[293] = (*a4)*v[108] + v[216] * v[527] + v[222] * v[529] + v[227] * v[531];
	v[109] = v[104] * v[225];
	v[312] = (*a1)*v[109] + v[217] * v[526] + v[223] * v[528] + v[228] * v[530];
	v[294] = (*a4)*v[109] + v[217] * v[527] + v[223] * v[529] + v[228] * v[531];
	v[111] = v[104] * v[229];
	v[313] = (*a1)*v[111] + v[230] * v[526] + v[234] * v[528] + v[239] * v[530];
	v[295] = (*a4)*v[111] + v[230] * v[527] + v[234] * v[529] + v[239] * v[531];
	v[113] = 1e0 - v[233] * v[377];
	v[314] = (*a1)*v[113] + v[231] * v[526] + v[235] * v[528] + v[240] * v[530];
	v[296] = (*a4)*v[113] + v[231] * v[527] + v[235] * v[529] + v[240] * v[531];
	v[114] = v[104] * v[238];
	v[315] = (*a1)*v[114] + v[232] * v[526] + v[237] * v[528] + v[241] * v[530];
	v[297] = (*a4)*v[114] + v[232] * v[527] + v[237] * v[529] + v[241] * v[531];
	v[116] = v[104] * v[242];
	v[316] = (*a1)*v[116] + v[243] * v[526] + v[247] * v[528] + v[252] * v[530];
	v[301] = (*a4)*v[116] + v[243] * v[527] + v[247] * v[529] + v[252] * v[531];
	v[118] = v[104] * v[246];
	v[317] = (*a1)*v[118] + v[244] * v[526] + v[248] * v[528] + v[254] * v[530];
	v[302] = (*a4)*v[118] + v[244] * v[527] + v[248] * v[529] + v[254] * v[531];
	v[119] = 1e0 - v[251] * v[377];
	v[318] = (*a1)*v[119] + v[245] * v[526] + v[249] * v[528] + v[255] * v[530];
	v[303] = (*a4)*v[119] + v[245] * v[527] + v[249] * v[529] + v[255] * v[531];
	v[120] = 4e0 / (4e0 + v[510] + v[535]);
	v[536] = -0.5e0*v[120];
	v[123] = 1e0 + v[535] * v[536];
	v[124] = v[120] * (-alphai[2] + v[511]);
	v[125] = v[120] * (alphai[1] + v[514]);
	v[127] = v[120] * (alphai[2] + v[511]);
	v[129] = 1e0 + (v[510] + v[515])*v[536];
	v[130] = v[120] * (-alphai[0] + v[133]);
	v[132] = v[120] * (-alphai[1] + v[514]);
	v[290] = v[123] * v[245] + v[127] * v[249] + v[132] * v[255];
	v[286] = v[123] * v[244] + v[127] * v[248] + v[132] * v[254];
	v[282] = v[123] * v[243] + v[127] * v[247] + v[132] * v[252];
	v[278] = v[123] * v[232] + v[127] * v[237] + v[132] * v[241];
	v[274] = v[123] * v[231] + v[127] * v[235] + v[132] * v[240];
	v[270] = v[123] * v[230] + v[127] * v[234] + v[132] * v[239];
	v[266] = v[123] * v[217] + v[127] * v[223] + v[132] * v[228];
	v[262] = v[123] * v[216] + v[127] * v[222] + v[132] * v[227];
	v[258] = v[123] * v[215] + v[127] * v[220] + v[132] * v[226];
	v[134] = v[120] * (alphai[0] + v[133]);
	v[289] = v[124] * v[245] + v[129] * v[249] + v[134] * v[255];
	v[285] = v[124] * v[244] + v[129] * v[248] + v[134] * v[254];
	v[281] = v[124] * v[243] + v[129] * v[247] + v[134] * v[252];
	v[277] = v[124] * v[232] + v[129] * v[237] + v[134] * v[241];
	v[273] = v[124] * v[231] + v[129] * v[235] + v[134] * v[240];
	v[269] = v[124] * v[230] + v[129] * v[234] + v[134] * v[239];
	v[265] = v[124] * v[217] + v[129] * v[223] + v[134] * v[228];
	v[261] = v[124] * v[216] + v[129] * v[222] + v[134] * v[227];
	v[257] = v[124] * v[215] + v[129] * v[220] + v[134] * v[226];
	v[135] = 1e0 + (v[510] + v[512])*v[536];
	v[288] = v[125] * v[245] + v[130] * v[249] + v[135] * v[255];
	v[372] = Jr[2][0] * v[288] + Jr[1][0] * v[289] + Jr[0][0] * v[290];
	v[369] = Jr[2][1] * v[288] + Jr[1][1] * v[289] + Jr[0][1] * v[290];
	v[366] = Jr[2][2] * v[288] + Jr[1][2] * v[289] + Jr[0][2] * v[290];
	v[291] = br[2] * v[288] + br[1] * v[289] + br[0] * v[290];
	v[284] = v[125] * v[244] + v[130] * v[248] + v[135] * v[254];
	v[371] = Jr[2][0] * v[284] + Jr[1][0] * v[285] + Jr[0][0] * v[286];
	v[368] = Jr[2][1] * v[284] + Jr[1][1] * v[285] + Jr[0][1] * v[286];
	v[365] = Jr[2][2] * v[284] + Jr[1][2] * v[285] + Jr[0][2] * v[286];
	v[287] = br[2] * v[284] + br[1] * v[285] + br[0] * v[286];
	v[280] = v[125] * v[243] + v[130] * v[247] + v[135] * v[252];
	v[370] = Jr[2][0] * v[280] + Jr[1][0] * v[281] + Jr[0][0] * v[282];
	v[367] = Jr[2][1] * v[280] + Jr[1][1] * v[281] + Jr[0][1] * v[282];
	v[364] = Jr[2][2] * v[280] + Jr[1][2] * v[281] + Jr[0][2] * v[282];
	v[283] = br[2] * v[280] + br[1] * v[281] + br[0] * v[282];
	v[276] = v[125] * v[232] + v[130] * v[237] + v[135] * v[241];
	v[363] = Jr[2][0] * v[276] + Jr[1][0] * v[277] + Jr[0][0] * v[278];
	v[360] = Jr[2][1] * v[276] + Jr[1][1] * v[277] + Jr[0][1] * v[278];
	v[357] = Jr[2][2] * v[276] + Jr[1][2] * v[277] + Jr[0][2] * v[278];
	v[279] = br[2] * v[276] + br[1] * v[277] + br[0] * v[278];
	v[272] = v[125] * v[231] + v[130] * v[235] + v[135] * v[240];
	v[362] = Jr[2][0] * v[272] + Jr[1][0] * v[273] + Jr[0][0] * v[274];
	v[359] = Jr[2][1] * v[272] + Jr[1][1] * v[273] + Jr[0][1] * v[274];
	v[356] = Jr[2][2] * v[272] + Jr[1][2] * v[273] + Jr[0][2] * v[274];
	v[275] = br[2] * v[272] + br[1] * v[273] + br[0] * v[274];
	v[268] = v[125] * v[230] + v[130] * v[234] + v[135] * v[239];
	v[361] = Jr[2][0] * v[268] + Jr[1][0] * v[269] + Jr[0][0] * v[270];
	v[358] = Jr[2][1] * v[268] + Jr[1][1] * v[269] + Jr[0][1] * v[270];
	v[355] = Jr[2][2] * v[268] + Jr[1][2] * v[269] + Jr[0][2] * v[270];
	v[271] = br[2] * v[268] + br[1] * v[269] + br[0] * v[270];
	v[264] = v[125] * v[217] + v[130] * v[223] + v[135] * v[228];
	v[354] = Jr[2][0] * v[264] + Jr[1][0] * v[265] + Jr[0][0] * v[266];
	v[351] = Jr[2][1] * v[264] + Jr[1][1] * v[265] + Jr[0][1] * v[266];
	v[348] = Jr[2][2] * v[264] + Jr[1][2] * v[265] + Jr[0][2] * v[266];
	v[267] = br[2] * v[264] + br[1] * v[265] + br[0] * v[266];
	v[260] = v[125] * v[216] + v[130] * v[222] + v[135] * v[227];
	v[353] = Jr[2][0] * v[260] + Jr[1][0] * v[261] + Jr[0][0] * v[262];
	v[350] = Jr[2][1] * v[260] + Jr[1][1] * v[261] + Jr[0][1] * v[262];
	v[347] = Jr[2][2] * v[260] + Jr[1][2] * v[261] + Jr[0][2] * v[262];
	v[263] = br[2] * v[260] + br[1] * v[261] + br[0] * v[262];
	v[256] = v[125] * v[215] + v[130] * v[220] + v[135] * v[226];
	v[352] = Jr[2][0] * v[256] + Jr[1][0] * v[257] + Jr[0][0] * v[258];
	v[349] = Jr[2][1] * v[256] + Jr[1][1] * v[257] + Jr[0][1] * v[258];
	v[346] = Jr[2][2] * v[256] + Jr[1][2] * v[257] + Jr[0][2] * v[258];
	v[259] = br[2] * v[256] + br[1] * v[257] + br[0] * v[258];
	v[136] = v[107] * v[123] + v[108] * v[127] + v[109] * v[132];
	v[137] = v[107] * v[124] + v[108] * v[129] + v[109] * v[134];
	v[138] = v[107] * v[125] + v[108] * v[130] + v[109] * v[135];
	v[151] = Jr[0][2] * v[136] + Jr[1][2] * v[137] + Jr[2][2] * v[138];
	v[150] = Jr[0][1] * v[136] + Jr[1][1] * v[137] + Jr[2][1] * v[138];
	v[149] = Jr[0][0] * v[136] + Jr[1][0] * v[137] + Jr[2][0] * v[138];
	v[382] = v[151] * v[264] + v[150] * v[265] + v[149] * v[266] + v[138] * v[348] + v[137] * v[351] + v[136] * v[354];
	v[381] = v[151] * v[260] + v[150] * v[261] + v[149] * v[262] + v[138] * v[347] + v[137] * v[350] + v[136] * v[353];
	v[380] = v[151] * v[256] + v[150] * v[257] + v[149] * v[258] + v[138] * v[346] + v[137] * v[349] + v[136] * v[352];
	v[139] = v[111] * v[123] + v[113] * v[127] + v[114] * v[132];
	v[140] = v[111] * v[124] + v[113] * v[129] + v[114] * v[134];
	v[141] = v[111] * v[125] + v[113] * v[130] + v[114] * v[135];
	v[385] = v[151] * v[276] + v[150] * v[277] + v[149] * v[278] + v[141] * v[348] + v[140] * v[351] + v[139] * v[354];
	v[384] = v[151] * v[272] + v[150] * v[273] + v[149] * v[274] + v[141] * v[347] + v[140] * v[350] + v[139] * v[353];
	v[383] = v[151] * v[268] + v[150] * v[269] + v[149] * v[270] + v[141] * v[346] + v[140] * v[349] + v[139] * v[352];
	v[157] = Jr[0][2] * v[139] + Jr[1][2] * v[140] + Jr[2][2] * v[141];
	v[156] = Jr[0][1] * v[139] + Jr[1][1] * v[140] + Jr[2][1] * v[141];
	v[155] = Jr[0][0] * v[139] + Jr[1][0] * v[140] + Jr[2][0] * v[141];
	v[394] = v[157] * v[276] + v[156] * v[277] + v[155] * v[278] + v[141] * v[357] + v[140] * v[360] + v[139] * v[363];
	v[393] = v[157] * v[272] + v[156] * v[273] + v[155] * v[274] + v[141] * v[356] + v[140] * v[359] + v[139] * v[362];
	v[392] = v[157] * v[268] + v[156] * v[269] + v[155] * v[270] + v[141] * v[355] + v[140] * v[358] + v[139] * v[361];
	v[391] = v[157] * v[264] + v[156] * v[265] + v[155] * v[266] + v[138] * v[357] + v[137] * v[360] + v[136] * v[363];
	v[390] = v[157] * v[260] + v[156] * v[261] + v[155] * v[262] + v[138] * v[356] + v[137] * v[359] + v[136] * v[362];
	v[389] = v[157] * v[256] + v[156] * v[257] + v[155] * v[258] + v[138] * v[355] + v[137] * v[358] + v[136] * v[361];
	v[142] = v[116] * v[123] + v[118] * v[127] + v[119] * v[132];
	v[143] = v[116] * v[124] + v[118] * v[129] + v[119] * v[134];
	v[144] = v[116] * v[125] + v[118] * v[130] + v[119] * v[135];
	v[397] = v[157] * v[288] + v[156] * v[289] + v[155] * v[290] + v[144] * v[357] + v[143] * v[360] + v[142] * v[363];
	v[396] = v[157] * v[284] + v[156] * v[285] + v[155] * v[286] + v[144] * v[356] + v[143] * v[359] + v[142] * v[362];
	v[395] = v[157] * v[280] + v[156] * v[281] + v[155] * v[282] + v[144] * v[355] + v[143] * v[358] + v[142] * v[361];
	v[388] = v[151] * v[288] + v[150] * v[289] + v[149] * v[290] + v[144] * v[348] + v[143] * v[351] + v[142] * v[354];
	v[387] = v[151] * v[284] + v[150] * v[285] + v[149] * v[286] + v[144] * v[347] + v[143] * v[350] + v[142] * v[353];
	v[386] = v[151] * v[280] + v[150] * v[281] + v[149] * v[282] + v[144] * v[346] + v[143] * v[349] + v[142] * v[352];
	v[163] = Jr[0][2] * v[142] + Jr[1][2] * v[143] + Jr[2][2] * v[144];
	v[162] = Jr[0][1] * v[142] + Jr[1][1] * v[143] + Jr[2][1] * v[144];
	v[161] = Jr[0][0] * v[142] + Jr[1][0] * v[143] + Jr[2][0] * v[144];
	v[406] = v[163] * v[288] + v[162] * v[289] + v[161] * v[290] + v[144] * v[366] + v[143] * v[369] + v[142] * v[372];
	v[405] = v[163] * v[284] + v[162] * v[285] + v[161] * v[286] + v[144] * v[365] + v[143] * v[368] + v[142] * v[371];
	v[404] = v[163] * v[280] + v[162] * v[281] + v[161] * v[282] + v[144] * v[364] + v[143] * v[367] + v[142] * v[370];
	v[403] = v[163] * v[276] + v[162] * v[277] + v[161] * v[278] + v[141] * v[366] + v[140] * v[369] + v[139] * v[372];
	v[402] = v[163] * v[272] + v[162] * v[273] + v[161] * v[274] + v[141] * v[365] + v[140] * v[368] + v[139] * v[371];
	v[401] = v[163] * v[268] + v[162] * v[269] + v[161] * v[270] + v[141] * v[364] + v[140] * v[367] + v[139] * v[370];
	v[400] = v[163] * v[264] + v[162] * v[265] + v[161] * v[266] + v[138] * v[366] + v[137] * v[369] + v[136] * v[372];
	v[399] = v[163] * v[260] + v[162] * v[261] + v[161] * v[262] + v[138] * v[365] + v[137] * v[368] + v[136] * v[371];
	v[398] = v[163] * v[256] + v[162] * v[257] + v[161] * v[258] + v[138] * v[364] + v[137] * v[367] + v[136] * v[370];
	v[148] = v[136] * v[149] + v[137] * v[150] + v[138] * v[151];
	v[462] = v[148] * v[374];
	v[152] = v[139] * v[149] + v[140] * v[150] + v[141] * v[151];
	v[464] = v[152] * v[374];
	v[153] = v[142] * v[149] + v[143] * v[150] + v[144] * v[151];
	v[466] = v[153] * v[374];
	v[154] = v[136] * v[155] + v[137] * v[156] + v[138] * v[157];
	v[461] = -(v[154] * v[373]);
	v[158] = v[139] * v[155] + v[140] * v[156] + v[141] * v[157];
	v[463] = -(v[158] * v[373]);
	v[159] = v[142] * v[155] + v[143] * v[156] + v[144] * v[157];
	v[465] = -(v[159] * v[373]);
	v[160] = v[136] * v[161] + v[137] * v[162] + v[138] * v[163];
	v[449] = -(v[160] * v[376]);
	v[442] = v[104] * v[160] + v[148] * v[218] + v[154] * v[221];
	v[438] = v[104] * v[154] - v[160] * v[221] + v[148] * v[224];
	v[434] = v[104] * v[148] - v[160] * v[218] - v[154] * v[224];
	v[164] = v[139] * v[161] + v[140] * v[162] + v[141] * v[163];
	v[450] = -(v[164] * v[376]);
	v[443] = v[104] * v[164] + v[152] * v[218] + v[158] * v[221];
	v[439] = v[104] * v[158] - v[164] * v[221] + v[152] * v[224];
	v[435] = v[104] * v[152] - v[164] * v[218] - v[158] * v[224];
	v[165] = v[142] * v[161] + v[143] * v[162] + v[144] * v[163];
	v[451] = -(v[165] * v[376]);
	v[444] = v[104] * v[165] + v[153] * v[218] + v[159] * v[221];
	v[440] = v[104] * v[159] - v[165] * v[221] + v[153] * v[224];
	v[436] = v[104] * v[153] - v[165] * v[218] - v[159] * v[224];
	v[166] = br[0] * v[136] + br[1] * v[137] + br[2] * v[138];
	v[496] = (v[166] * v[166]);
	v[329] = (*a1)*v[166];
	v[344] = -(v[221] * v[329]);
	v[341] = v[104] * v[329];
	v[167] = br[0] * v[139] + br[1] * v[140] + br[2] * v[141];
	v[537] = (*m)*v[167];
	v[497] = (v[167] * v[167]);
	v[491] = v[166] * v[537];
	v[328] = -((*a1)*v[167]);
	v[343] = -(v[218] * v[328]);
	v[339] = v[104] * v[328];
	v[168] = br[0] * v[142] + br[1] * v[143] + br[2] * v[144];
	v[495] = v[168] * v[537];
	v[494] = (*m)*v[166] * v[168];
	v[492] = (v[168] * v[168]);
	v[330] = -((*a1)*v[168]);
	v[336] = v[224] * v[330];
	v[334] = -(v[104] * v[330]);
	v[408] = -(v[275] * v[541]) + v[263] * v[542];
	v[407] = -(v[271] * v[541]) + v[259] * v[542];
	v[194] = -(v[167] * v[541]) + v[166] * v[542];
	v[452] = -(v[194] * v[376]);
	v[414] = v[287] * v[541] - v[263] * v[543];
	v[413] = v[283] * v[541] - v[259] * v[543];
	v[411] = -(v[287] * v[542]) + v[275] * v[543];
	v[410] = -(v[283] * v[542]) + v[271] * v[543];
	v[196] = -(v[168] * v[542]) + v[167] * v[543];
	v[468] = v[196] * v[374];
	v[195] = v[168] * v[541] - v[166] * v[543];
	v[467] = -(v[195] * v[373]);
	v[457] = (*m)*(v[194] * v[212] - v[195] * v[376] + v[196] * v[378] + v[104] * v[408] + v[218] * v[411] + v[221] * v[414]);
	v[455] = (*m)*(v[194] * v[210] + v[196] * v[376] + v[195] * v[379] + v[104] * v[407] + v[218] * v[410] + v[221] * v[413]);
	v[446] = (*m)*(v[195] * v[210] + v[196] * v[373] - v[194] * v[379] - v[221] * v[407] + v[224] * v[410] + v[104] * v[413]);
	v[175] = v[107] * v[527] + v[108] * v[529] + v[109] * v[531];
	v[486] = v[167] * v[175];
	v[483] = -(v[168] * v[175]);
	v[475] = -(v[160] * v[175]);
	v[473] = v[154] * v[175];
	v[179] = v[111] * v[527] + v[113] * v[529] + v[114] * v[531];
	v[487] = -(v[166] * v[179]);
	v[480] = v[168] * v[179];
	v[476] = v[164] * v[179];
	v[471] = -(v[152] * v[179]);
	v[430] = v[159] * v[175] - v[153] * v[179];
	v[429] = v[158] * v[175] + v[471];
	v[428] = -(v[148] * v[179]) + v[473];
	v[300] = -(v[179] * v[267]) + v[175] * v[279] + v[167] * v[294] - v[166] * v[297];
	v[299] = -(v[179] * v[263]) + v[175] * v[275] + v[167] * v[293] - v[166] * v[296];
	v[298] = -(v[179] * v[259]) + v[175] * v[271] + v[167] * v[292] - v[166] * v[295];
	v[188] = v[486] + v[487];
	v[180] = v[116] * v[527] + v[118] * v[529] + v[119] * v[531];
	v[484] = v[166] * v[180];
	v[481] = -(v[167] * v[180]);
	v[474] = -(v[159] * v[180]);
	v[472] = v[153] * v[180];
	v[433] = v[175] * (v[154] * v[294] - v[148] * v[297] - v[179] * v[382] + v[175] * v[391]) + v[179] * (v[158] * v[294]
		- v[152] * v[297] - v[179] * v[385] + v[175] * v[394]) + v[180] * (v[159] * v[294] - v[153] * v[297] - v[179] * v[388]
		+ v[175] * v[397]) + v[294] * v[428] + v[297] * v[429] + v[303] * v[430];
	v[432] = v[175] * (v[154] * v[293] - v[148] * v[296] - v[179] * v[381] + v[175] * v[390]) + v[179] * (v[158] * v[293]
		- v[152] * v[296] - v[179] * v[384] + v[175] * v[393]) + v[180] * (v[159] * v[293] - v[153] * v[296] - v[179] * v[387]
		+ v[175] * v[396]) + v[293] * v[428] + v[296] * v[429] + v[302] * v[430];
	v[431] = v[175] * (v[154] * v[292] - v[148] * v[295] - v[179] * v[380] + v[175] * v[389]) + v[179] * (v[158] * v[292]
		- v[152] * v[295] - v[179] * v[383] + v[175] * v[392]) + v[180] * (v[159] * v[292] - v[153] * v[295] - v[179] * v[386]
		+ v[175] * v[395]) + v[292] * v[428] + v[295] * v[429] + v[301] * v[430];
	v[424] = -(v[165] * v[175]) + v[472];
	v[423] = -(v[164] * v[175]) + v[152] * v[180];
	v[422] = v[148] * v[180] + v[475];
	v[427] = v[175] * (-(v[160] * v[294]) + v[148] * v[303] + v[180] * v[382] - v[175] * v[400]) + v[179] * (-(v[164] * v[294])
		+ v[152] * v[303] + v[180] * v[385] - v[175] * v[403]) + v[180] * (-(v[165] * v[294]) + v[153] * v[303] + v[180] * v[388]
		- v[175] * v[406]) + v[294] * v[422] + v[297] * v[423] + v[303] * v[424];
	v[426] = v[175] * (-(v[160] * v[293]) + v[148] * v[302] + v[180] * v[381] - v[175] * v[399]) + v[179] * (-(v[164] * v[293])
		+ v[152] * v[302] + v[180] * v[384] - v[175] * v[402]) + v[180] * (-(v[165] * v[293]) + v[153] * v[302] + v[180] * v[387]
		- v[175] * v[405]) + v[293] * v[422] + v[296] * v[423] + v[302] * v[424];
	v[425] = v[175] * (-(v[160] * v[292]) + v[148] * v[301] + v[180] * v[380] - v[175] * v[398]) + v[179] * (-(v[164] * v[292])
		+ v[152] * v[301] + v[180] * v[383] - v[175] * v[401]) + v[180] * (-(v[165] * v[292]) + v[153] * v[301] + v[180] * v[386]
		- v[175] * v[404]) + v[292] * v[422] + v[295] * v[423] + v[301] * v[424];
	v[418] = v[165] * v[179] + v[474];
	v[417] = -(v[158] * v[180]) + v[476];
	v[416] = v[160] * v[179] - v[154] * v[180];
	v[421] = v[175] * (v[160] * v[297] - v[154] * v[303] - v[180] * v[391] + v[179] * v[400]) + v[179] * (v[164] * v[297]
		- v[158] * v[303] - v[180] * v[394] + v[179] * v[403]) + v[180] * (v[165] * v[297] - v[159] * v[303] - v[180] * v[397]
		+ v[179] * v[406]) + v[294] * v[416] + v[297] * v[417] + v[303] * v[418];
	v[420] = v[175] * (v[160] * v[296] - v[154] * v[302] - v[180] * v[390] + v[179] * v[399]) + v[179] * (v[164] * v[296]
		- v[158] * v[302] - v[180] * v[393] + v[179] * v[402]) + v[180] * (v[165] * v[296] - v[159] * v[302] - v[180] * v[396]
		+ v[179] * v[405]) + v[293] * v[416] + v[296] * v[417] + v[302] * v[418];
	v[419] = v[175] * (v[160] * v[295] - v[154] * v[301] - v[180] * v[389] + v[179] * v[398]) + v[179] * (v[164] * v[295]
		- v[158] * v[301] - v[180] * v[392] + v[179] * v[401]) + v[180] * (v[165] * v[295] - v[159] * v[301] - v[180] * v[395]
		+ v[179] * v[404]) + v[292] * v[416] + v[295] * v[417] + v[301] * v[418];
	v[309] = v[180] * v[267] - v[175] * v[291] - v[168] * v[294] + v[166] * v[303];
	v[308] = v[180] * v[263] - v[175] * v[287] - v[168] * v[293] + v[166] * v[302];
	v[307] = v[180] * v[259] - v[175] * v[283] - v[168] * v[292] + v[166] * v[301];
	v[306] = -(v[180] * v[279]) + v[179] * v[291] + v[168] * v[297] - v[167] * v[303];
	v[305] = -(v[180] * v[275]) + v[179] * v[287] + v[168] * v[296] - v[167] * v[302];
	v[304] = -(v[180] * v[271]) + v[179] * v[283] + v[168] * v[295] - v[167] * v[301];
	v[199] = v[175] * v[416] + v[179] * v[417] + v[180] * v[418];
	v[460] = v[199] * v[374];
	v[198] = v[175] * v[422] + v[179] * v[423] + v[180] * v[424];
	v[459] = -(v[198] * v[373]);
	v[197] = v[175] * v[428] + v[179] * v[429] + v[180] * v[430];
	v[448] = -(v[197] * v[376]);
	v[191] = v[480] + v[481];
	v[190] = v[483] + v[484];
	v[181] = v[107] * v[526] + v[108] * v[528] + v[109] * v[530];
	v[185] = v[111] * v[526] + v[113] * v[528] + v[114] * v[530];
	v[186] = v[116] * v[526] + v[118] * v[528] + v[119] * v[530];
	v[470] = 0.5e0*(*m)*((v[538] * v[538]) + (v[539] * v[539]) + (v[540] * v[540]));
	v[477] = 0.5e0*(v[175] * (v[148] * v[175] - v[471] + v[472]) + v[179] * (v[158] * v[179] + v[473] - v[474]) + v[180] *
		(v[165] * v[180] - v[475] + v[476]));
	v[478] = (*m)*(v[191] * v[538] + v[190] * v[539] + v[188] * v[540]);
	v[482] = (*m)*(v[191] + v[538]);
	v[485] = (*m)*(v[190] + v[539]);
	v[488] = (*m)*(v[188] + v[540]);
	v[490] = v[179] * (v[152] + v[491]) + v[180] * (v[153] + v[494]) + v[175] * (v[148] - (*m)*(v[492] + v[497]));
	v[493] = v[175] * (v[154] + v[491]) + v[180] * (v[159] + v[495]) + v[179] * (v[158] - (*m)*(v[492] + v[496]));
	v[498] = v[175] * (v[160] + v[494]) + v[179] * (v[164] + v[495]) + v[180] * (v[165] - (*m)*(v[496] + v[497]));
	dT[0] = (*m)*(v[168] * v[185] - v[167] * v[186] + v[179] * v[188] - v[180] * v[190] + v[541]);
	dT[1] = (*m)*(-(v[168] * v[181]) + v[166] * v[186] - v[175] * v[188] + v[180] * v[191] + v[542]);
	dT[2] = (*m)*(v[167] * v[181] - v[166] * v[185] + v[175] * v[190] - v[179] * v[191] + v[543]);
	dT[3] = v[104] * v[199] - v[197] * v[218] - v[198] * v[224] + (*m)*(v[104] * v[196] - v[194] * v[218] - v[195] * v[224])
		+ v[181] * v[434] + v[185] * v[435] + v[186] * v[436];
	dT[4] = v[104] * v[198] - v[197] * v[221] + v[199] * v[224] + (*m)*(v[104] * v[195] - v[194] * v[221] + v[196] * v[224])
		+ v[181] * v[438] + v[185] * v[439] + v[186] * v[440];
	dT[5] = v[104] * v[197] + v[199] * v[218] + v[198] * v[221] + (*m)*(v[104] * v[194] + v[196] * v[218] + v[195] * v[221])
		+ v[181] * v[442] + v[185] * v[443] + v[186] * v[444];
	DdT[0][0] = v[544];
	DdT[0][1] = 0e0;
	DdT[0][2] = 0e0;
	DdT[0][3] = (*m)*(-(v[186] * v[271]) + v[185] * v[283] + v[188] * v[295] + v[179] * v[298] - v[190] * v[301]
		- v[180] * v[307] + v[168] * v[313] - v[167] * v[316]);
	DdT[0][4] = (*m)*(-(v[186] * v[275]) + v[185] * v[287] + v[188] * v[296] + v[179] * v[299] - v[190] * v[302]
		- v[180] * v[308] + v[168] * v[314] - v[167] * v[317]);
	DdT[0][5] = (*m)*(-(v[186] * v[279]) + v[185] * v[291] + v[188] * v[297] + v[179] * v[300] - v[190] * v[303]
		- v[180] * v[309] + v[168] * v[315] - v[167] * v[318]);
	DdT[1][0] = 0e0;
	DdT[1][1] = v[544];
	DdT[1][2] = 0e0;
	DdT[1][3] = (*m)*(v[186] * v[259] - v[181] * v[283] - v[188] * v[292] - v[175] * v[298] + v[191] * v[301] + v[180] * v[304]
		- v[168] * v[310] + v[166] * v[316]);
	DdT[1][4] = (*m)*(v[186] * v[263] - v[181] * v[287] - v[188] * v[293] - v[175] * v[299] + v[191] * v[302] + v[180] * v[305]
		- v[168] * v[311] + v[166] * v[317]);
	DdT[1][5] = (*m)*(v[186] * v[267] - v[181] * v[291] - v[188] * v[294] - v[175] * v[300] + v[191] * v[303] + v[180] * v[306]
		- v[168] * v[312] + v[166] * v[318]);
	DdT[2][0] = 0e0;
	DdT[2][1] = 0e0;
	DdT[2][2] = v[544];
	DdT[2][3] = (*m)*(-(v[185] * v[259]) + v[181] * v[271] + v[190] * v[292] - v[191] * v[295] - v[179] * v[304]
		+ v[175] * v[307] + v[167] * v[310] - v[166] * v[313]);
	DdT[2][4] = (*m)*(-(v[185] * v[263]) + v[181] * v[275] + v[190] * v[293] - v[191] * v[296] - v[179] * v[305]
		+ v[175] * v[308] + v[167] * v[311] - v[166] * v[314]);
	DdT[2][5] = (*m)*(-(v[185] * v[267]) + v[181] * v[279] + v[190] * v[294] - v[191] * v[297] - v[179] * v[306]
		+ v[175] * v[309] + v[167] * v[312] - v[166] * v[315]);
	DdT[3][0] = (*m)*(v[336] + v[343]);
	DdT[3][1] = (*m)*(-(v[218] * v[329]) - v[334]);
	DdT[3][2] = (*m)*(v[224] * v[329] - v[339]);
	DdT[3][3] = v[199] * v[210] + v[104] * v[419] - v[224] * v[425] - v[218] * v[431] + v[310] * v[434] + v[313] * v[435]
		+ v[316] * v[436] + v[448] + v[459] + v[181] * (v[148] * v[210] + v[104] * v[380] - v[224] * v[389] - v[218] * v[398] + v[449]
		+ v[461]) + v[185] * (v[152] * v[210] + v[104] * v[383] - v[224] * v[392] - v[218] * v[401] + v[450] + v[463]) + v[186] *
		(v[153] * v[210] + v[104] * v[386] - v[224] * v[395] - v[218] * v[404] + v[451] + v[465]) + (*m)*(v[196] * v[210]
		- v[218] * v[407] + v[104] * v[410] - v[224] * v[413] + v[452] + v[467]);
	DdT[3][4] = v[199] * v[212] - v[198] * v[374] - v[197] * v[378] + v[181] * (v[148] * v[212] - v[154] * v[374]
		- v[160] * v[378] + v[104] * v[381] - v[224] * v[390] - v[218] * v[399]) + v[185] * (v[152] * v[212] - v[158] * v[374]
		- v[164] * v[378] + v[104] * v[384] - v[224] * v[393] - v[218] * v[402]) + v[186] * (v[153] * v[212] - v[159] * v[374]
		- v[165] * v[378] + v[104] * v[387] - v[224] * v[396] - v[218] * v[405]) + v[104] * v[420] - v[224] * v[426] - v[218] * v[432]
		+ v[311] * v[434] + v[314] * v[435] + v[317] * v[436] + v[446];
	DdT[3][5] = v[199] * v[213] + v[197] * v[374] - v[198] * v[375] + v[181] * (v[148] * v[213] + v[160] * v[374]
		- v[154] * v[375] + v[104] * v[382] - v[224] * v[391] - v[218] * v[400]) + v[185] * (v[152] * v[213] + v[164] * v[374]
		- v[158] * v[375] + v[104] * v[385] - v[224] * v[394] - v[218] * v[403]) + v[186] * (v[153] * v[213] + v[165] * v[374]
		- v[159] * v[375] + v[104] * v[388] - v[224] * v[397] - v[218] * v[406]) + v[104] * v[421] - v[224] * v[427] - v[218] * v[433]
		+ v[312] * v[434] + v[315] * v[435] + v[318] * v[436] + v[455];
	DdT[4][0] = (*m)*(-(v[221] * v[328]) + v[334]);
	DdT[4][1] = (*m)*(v[336] + v[344]);
	DdT[4][2] = (*m)*(-(v[224] * v[328]) - v[341]);
	DdT[4][3] = v[198] * v[210] + v[199] * v[373] - v[197] * v[379] + v[181] * (v[154] * v[210] + v[148] * v[373]
		- v[160] * v[379] + v[224] * v[380] + v[104] * v[389] - v[221] * v[398]) + v[185] * (v[158] * v[210] + v[152] * v[373]
		- v[164] * v[379] + v[224] * v[383] + v[104] * v[392] - v[221] * v[401]) + v[186] * (v[159] * v[210] + v[153] * v[373]
		- v[165] * v[379] + v[224] * v[386] + v[104] * v[395] - v[221] * v[404]) + v[224] * v[419] + v[104] * v[425] - v[221] * v[431]
		+ v[310] * v[438] + v[313] * v[439] + v[316] * v[440] + v[446];
	DdT[4][4] = v[198] * v[212] + v[224] * v[420] + v[104] * v[426] - v[221] * v[432] + v[311] * v[438] + v[314] * v[439]
		+ v[317] * v[440] - v[448] + v[460] + v[181] * (v[154] * v[212] + v[224] * v[381] + v[104] * v[390] - v[221] * v[399] - v[449]
		+ v[462]) + v[185] * (v[158] * v[212] + v[224] * v[384] + v[104] * v[393] - v[221] * v[402] - v[450] + v[464]) + v[186] *
		(v[159] * v[212] + v[224] * v[387] + v[104] * v[396] - v[221] * v[405] - v[451] + v[466]) + (*m)*(v[195] * v[212]
		- v[221] * v[408] + v[224] * v[411] + v[104] * v[414] - v[452] + v[468]);
	DdT[4][5] = v[198] * v[213] - v[197] * v[373] + v[199] * v[375] + v[181] * (v[154] * v[213] - v[160] * v[373]
		+ v[148] * v[375] + v[224] * v[382] + v[104] * v[391] - v[221] * v[400]) + v[185] * (v[158] * v[213] - v[164] * v[373]
		+ v[152] * v[375] + v[224] * v[385] + v[104] * v[394] - v[221] * v[403]) + v[186] * (v[159] * v[213] - v[165] * v[373]
		+ v[153] * v[375] + v[224] * v[388] + v[104] * v[397] - v[221] * v[406]) + v[224] * v[421] + v[104] * v[427] - v[221] * v[433]
		+ v[312] * v[438] + v[315] * v[439] + v[318] * v[440] + v[457];
	DdT[5][0] = (*m)*(-(v[221] * v[330]) + v[339]);
	DdT[5][1] = (*m)*(v[218] * v[330] + v[341]);
	DdT[5][2] = (*m)*(v[343] + v[344]);
	DdT[5][3] = v[197] * v[210] + v[199] * v[376] + v[198] * v[379] + v[181] * (v[160] * v[210] + v[148] * v[376]
		+ v[154] * v[379] + v[218] * v[380] + v[221] * v[389] + v[104] * v[398]) + v[185] * (v[164] * v[210] + v[152] * v[376]
		+ v[158] * v[379] + v[218] * v[383] + v[221] * v[392] + v[104] * v[401]) + v[186] * (v[165] * v[210] + v[153] * v[376]
		+ v[159] * v[379] + v[218] * v[386] + v[221] * v[395] + v[104] * v[404]) + v[218] * v[419] + v[221] * v[425] + v[104] * v[431]
		+ v[310] * v[442] + v[313] * v[443] + v[316] * v[444] + v[455];
	DdT[5][4] = v[197] * v[212] - v[198] * v[376] + v[199] * v[378] + v[181] * (v[160] * v[212] - v[154] * v[376]
		+ v[148] * v[378] + v[218] * v[381] + v[221] * v[390] + v[104] * v[399]) + v[185] * (v[164] * v[212] - v[158] * v[376]
		+ v[152] * v[378] + v[218] * v[384] + v[221] * v[393] + v[104] * v[402]) + v[186] * (v[165] * v[212] - v[159] * v[376]
		+ v[153] * v[378] + v[218] * v[387] + v[221] * v[396] + v[104] * v[405]) + v[218] * v[420] + v[221] * v[426] + v[104] * v[432]
		+ v[311] * v[442] + v[314] * v[443] + v[317] * v[444] + v[457];
	DdT[5][5] = v[197] * v[213] + v[218] * v[421] + v[221] * v[427] + v[104] * v[433] + v[312] * v[442] + v[315] * v[443]
		+ v[318] * v[444] - v[459] - v[460] + v[181] * (v[160] * v[213] + v[218] * v[382] + v[221] * v[391] + v[104] * v[400] - v[461]
		- v[462]) + v[185] * (v[164] * v[213] + v[218] * v[385] + v[221] * v[394] + v[104] * v[403] - v[463] - v[464]) + v[186] *
		(v[165] * v[213] + v[218] * v[388] + v[221] * v[397] + v[104] * v[406] - v[465] - v[466]) + (*m)*(v[194] * v[213] - v[467]
		- v[468] + v[104] * (-(v[279] * v[541]) + v[267] * v[542]) + v[221] * (v[291] * v[541] - v[267] * v[543]) + v[218] * (-
		(v[291] * v[542]) + v[279] * v[543]));
	(*T1) = v[470];
	(*T2) = v[477];
	(*T3) = v[478];
	(*T) = v[470] + v[477] + v[478];
	L[0] = v[482];
	L[1] = v[485];
	L[2] = v[488];
	(*magL) = sqrt(Power(v[482], 2) + Power(v[485], 2) + Power(v[488], 2));
	HG[0] = v[490];
	HG[1] = v[493];
	HG[2] = v[498];
	(*magHG) = sqrt(Power(v[490], 2) + Power(v[493], 2) + Power(v[498], 2));

	/*db.PrintPtr(dT, 6);
	db.PrintPtr(ud, 3);
	db.PrintPtr(dui, 3);
	db.PrintPtr(ddui, 3);
	db.PrintPtr(alphad, 3);
	db.PrintPtr(omegai, 3);
	db.PrintPtr(domegai, 3);
	db.PrintPtr(Jr, 3, 3);
	db.PrintPtr(br, 3);*/
}

void RigidBody_1::MountFieldLoading()
{
	//Ponteiro para o valor da massa
	double* m = &db.RB_data[RB_data_ID - 1]->mass;
	
	//Valores de variáveis cinemáticas - atribuição
	for (int i = 0; i < 3; i++)
	{
		alphai[i] = db.nodes[nodes[0] - 1]->copy_coordinates[i + 3];
		ud[i] = db.nodes[nodes[0] - 1]->displacements[i];
		alphad[i] = db.nodes[nodes[0] - 1]->displacements[i + 3];
	}

	//Variáveis para cálculo de steps
	double l_factor;
	
	if (db.environment_exist == true)
	{
		//Se existe campo gravitacional
		if (db.environment->g_exist == true)
		{
			l_factor = db.environment->bool_g.GetLinearFactorAtCurrentTime();
			double g[3];					//Vetor gravidade
			g[0] = l_factor * db.environment->G(0, 0);
			g[1] = l_factor * db.environment->G(1, 0);
			g[2] = l_factor * db.environment->G(2, 0);

			//AceGen
			int i01; int i02;
			v[257] = 0.5e0*alphad[2];
			v[255] = Power(alphad[2], 2);
			v[254] = alphad[0] * v[257];
			v[253] = 0.5e0*alphad[1];
			v[252] = 2e0*alphad[2];
			v[251] = Power(alphad[1], 2);
			v[250] = alphad[0] * v[253];
			v[249] = 2e0*alphad[1];
			v[248] = Power(alphad[0], 2);
			v[247] = 2e0*alphad[0];
			v[246] = Power(alphai[2], 2);
			v[245] = 0.5e0*alphai[0] * alphai[2];
			v[244] = 0.5e0*alphai[1];
			v[243] = Power(alphai[1], 2);
			v[261] = v[243] + v[246];
			v[242] = alphai[0] * v[244];
			v[241] = Power(alphai[0], 2);
			v[97] = alphai[2] * v[244];
			v[172] = -v[248] - v[251];
			v[150] = alphad[2] + v[250];
			v[140] = -alphad[2] + v[250];
			v[81] = alphad[2] * v[253];
			v[167] = alphad[0] + v[81];
			v[159] = -alphad[0] + v[81];
			v[163] = -alphad[1] + v[254];
			v[146] = alphad[1] + v[254];
			v[154] = -v[248] - v[255];
			v[135] = -v[251] - v[255];
			v[130] = 4e0 + v[248] + v[251] + v[255];
			v[132] = 1e0 / Power(v[130], 2);
			v[256] = -4e0*v[132];
			v[134] = v[252] * v[256];
			v[260] = 0.5e0*v[134];
			v[176] = v[172] * v[260];
			v[133] = v[249] * v[256];
			v[259] = 0.5e0*v[133];
			v[156] = v[154] * v[259];
			v[131] = v[247] * v[256];
			v[258] = 0.5e0*v[131];
			v[180] = v[131] * v[253];
			v[177] = -(v[131] * v[257]);
			v[136] = v[135] * v[258];
			v[68] = 4e0 / v[130];
			v[181] = -0.5e0*v[68];
			v[183] = v[181] - alphad[0] * v[258];
			v[174] = v[181] * v[249];
			v[175] = v[174] + v[172] * v[259];
			v[171] = v[181] * v[247];
			v[173] = v[171] + v[172] * v[258];
			v[168] = v[131] * v[167] + v[68];
			v[165] = v[133] * v[163] - v[68];
			v[160] = v[131] * v[159] - v[68];
			v[157] = v[181] * v[252];
			v[158] = v[157] + v[154] * v[260];
			v[155] = v[171] + v[154] * v[258];
			v[153] = v[134] * v[150] + v[68];
			v[148] = v[133] * v[146] + v[68];
			v[145] = alphad[2] * v[181];
			v[169] = -v[145] + v[133] * v[167];
			v[164] = -v[145] + v[131] * v[163];
			v[161] = -v[145] + v[133] * v[159];
			v[147] = -v[145] + v[131] * v[146];
			v[144] = v[134] * v[140] - v[68];
			v[142] = alphad[0] * v[181];
			v[166] = -v[142] + v[134] * v[163];
			v[152] = -v[142] + v[133] * v[150];
			v[149] = -v[142] + v[134] * v[146];
			v[143] = v[133] * v[140] - v[142];
			v[139] = -(alphad[1] * v[181]);
			v[170] = v[139] + v[134] * v[167];
			v[162] = v[139] + v[134] * v[159];
			v[151] = v[139] + v[131] * v[150];
			v[141] = v[139] + v[131] * v[140];
			v[138] = v[157] + v[135] * v[260];
			v[137] = v[174] + v[135] * v[259];
			v[71] = 1e0 - v[135] * v[181];
			v[72] = v[140] * v[68];
			v[73] = v[146] * v[68];
			v[75] = v[150] * v[68];
			v[77] = 1e0 - v[154] * v[181];
			v[78] = v[159] * v[68];
			v[80] = v[163] * v[68];
			v[82] = v[167] * v[68];
			v[83] = 1e0 - v[172] * v[181];
			v[84] = 4e0 / (4e0 + v[241] + v[261]);
			v[262] = -0.5e0*v[84];
			v[87] = 1e0 + v[261] * v[262];
			v[88] = (-alphai[2] + v[242])*v[84];
			v[89] = (alphai[1] + v[245])*v[84];
			v[91] = (alphai[2] + v[242])*v[84];
			v[93] = 1e0 + (v[241] + v[246])*v[262];
			v[94] = v[84] * (-alphai[0] + v[97]);
			v[96] = (-alphai[1] + v[245])*v[84];
			v[98] = v[84] * (alphai[0] + v[97]);
			v[99] = 1e0 + (v[241] + v[243])*v[262];
			v[222] = br[0] * (v[166] * v[87] + v[170] * v[91] + v[176] * v[96]) + br[1] * (v[166] * v[88] + v[170] * v[93] + v[176] * v[98]
				) + br[2] * (v[166] * v[89] + v[170] * v[94] + v[176] * v[99]);
			v[218] = br[0] * (v[165] * v[87] + v[169] * v[91] + v[175] * v[96]) + br[1] * (v[165] * v[88] + v[169] * v[93] + v[175] * v[98]
				) + br[2] * (v[165] * v[89] + v[169] * v[94] + v[175] * v[99]);
			v[214] = br[0] * (v[164] * v[87] + v[168] * v[91] + v[173] * v[96]) + br[1] * (v[164] * v[88] + v[168] * v[93] + v[173] * v[98]
				) + br[2] * (v[164] * v[89] + v[168] * v[94] + v[173] * v[99]);
			v[207] = br[0] * (v[153] * v[87] + v[158] * v[91] + v[162] * v[96]) + br[1] * (v[153] * v[88] + v[158] * v[93] + v[162] * v[98]
				) + br[2] * (v[153] * v[89] + v[158] * v[94] + v[162] * v[99]);
			v[203] = br[0] * (v[152] * v[87] + v[156] * v[91] + v[161] * v[96]) + br[1] * (v[152] * v[88] + v[156] * v[93] + v[161] * v[98]
				) + br[2] * (v[152] * v[89] + v[156] * v[94] + v[161] * v[99]);
			v[224] = (*m)*(g[2] * v[203] - g[1] * v[218]);
			v[199] = br[0] * (v[151] * v[87] + v[155] * v[91] + v[160] * v[96]) + br[1] * (v[151] * v[88] + v[155] * v[93] + v[160] * v[98]
				) + br[2] * (v[151] * v[89] + v[155] * v[94] + v[160] * v[99]);
			v[223] = (*m)*(g[2] * v[199] - g[1] * v[214]);
			v[195] = br[0] * (v[138] * v[87] + v[144] * v[91] + v[149] * v[96]) + br[1] * (v[138] * v[88] + v[144] * v[93] + v[149] * v[98]
				) + br[2] * (v[138] * v[89] + v[144] * v[94] + v[149] * v[99]);
			v[191] = br[0] * (v[137] * v[87] + v[143] * v[91] + v[148] * v[96]) + br[1] * (v[137] * v[88] + v[143] * v[93] + v[148] * v[98]
				) + br[2] * (v[137] * v[89] + v[143] * v[94] + v[148] * v[99]);
			v[227] = (*m)*(-(g[2] * v[191]) + g[0] * v[218]);
			v[209] = (*m)*(g[1] * v[191] - g[0] * v[203]);
			v[187] = br[0] * (v[136] * v[87] + v[141] * v[91] + v[147] * v[96]) + br[1] * (v[136] * v[88] + v[141] * v[93] + v[147] * v[98]
				) + br[2] * (v[136] * v[89] + v[141] * v[94] + v[147] * v[99]);
			v[226] = (*m)*(-(g[2] * v[187]) + g[0] * v[214]);
			v[208] = (*m)*(g[1] * v[187] - g[0] * v[199]);
			v[112] = br[0] * (v[71] * v[87] + v[72] * v[91] + v[73] * v[96]) + br[1] * (v[71] * v[88] + v[72] * v[93] + v[73] * v[98])
				+ br[2] * (v[71] * v[89] + v[72] * v[94] + v[73] * v[99]);
			v[113] = br[0] * (v[75] * v[87] + v[77] * v[91] + v[78] * v[96]) + br[1] * (v[75] * v[88] + v[77] * v[93] + v[78] * v[98])
				+ br[2] * (v[75] * v[89] + v[77] * v[94] + v[78] * v[99]);
			v[119] = (*m)*(g[1] * v[112] - g[0] * v[113]);
			v[232] = -(v[119] * v[180]);
			v[114] = br[0] * (v[80] * v[87] + v[82] * v[91] + v[83] * v[96]) + br[1] * (v[80] * v[88] + v[82] * v[93] + v[83] * v[98])
				+ br[2] * (v[80] * v[89] + v[82] * v[94] + v[83] * v[99]);
			v[121] = (*m)*(g[2] * v[113] - g[1] * v[114]);
			v[236] = -(alphad[2] * v[121] * v[259]);
			v[120] = (*m)*(-(g[2] * v[112]) + g[0] * v[114]);
			v[235] = -(v[120] * v[177]);
			dfield[0] = g[0] * (*m);
			dfield[1] = g[1] * (*m);
			dfield[2] = g[2] * (*m);
			dfield[3] = -(v[119] * v[139]) - v[120] * v[145] + v[121] * v[68];
			dfield[4] = -(v[119] * v[142]) + v[121] * v[145] + v[120] * v[68];
			dfield[5] = v[121] * v[139] + v[120] * v[142] + v[119] * v[68];
			Ddfield[0][0] = 0e0;
			Ddfield[0][1] = 0e0;
			Ddfield[0][2] = 0e0;
			Ddfield[0][3] = 0e0;
			Ddfield[0][4] = 0e0;
			Ddfield[0][5] = 0e0;
			Ddfield[1][1] = 0e0;
			Ddfield[1][2] = 0e0;
			Ddfield[1][3] = 0e0;
			Ddfield[1][4] = 0e0;
			Ddfield[1][5] = 0e0;
			Ddfield[2][2] = 0e0;
			Ddfield[2][3] = 0e0;
			Ddfield[2][4] = 0e0;
			Ddfield[2][5] = 0e0;
			Ddfield[3][3] = v[121] * v[131] - v[139] * v[208] - v[145] * v[226] + v[232] + v[235] + v[223] * v[68];
			Ddfield[3][4] = v[120] * v[131] + v[121] * v[177] - v[119] * v[183] - v[142] * v[208] + v[145] * v[223] + v[226] * v[68];
			Ddfield[3][5] = v[119] * v[131] + v[121] * v[180] + v[120] * v[183] + v[139] * v[223] + v[142] * v[226] + v[208] * v[68];
			Ddfield[4][4] = v[120] * v[133] - v[142] * v[209] + v[145] * v[224] - v[232] + v[236] + v[227] * v[68];
			Ddfield[4][5] = v[119] * v[133] - v[120] * v[180] + v[139] * v[224] + v[142] * v[227] + v[121] * (-v[181]
				+ alphad[1] * v[259]) + v[209] * v[68];
			Ddfield[5][5] = v[119] * v[134] - v[235] - v[236] + (*m)*(v[142] * (-(g[2] * v[195]) + g[0] * v[222]) + v[139] *
				(g[2] * v[207] - g[1] * v[222]) + (g[1] * v[195] - g[0] * v[207])*v[68]);
			for (i01 = 1; i01 < 6; i01++) {
				for (i02 = 0; i02 < i01; i02++) {
					Ddfield[i01][i02] = Ddfield[i02][i01];
				}
			};
		}
	}

	//db.PrintPtr(dfield, 6);

	//Outras contribuições de esforços de campo, colocar aqui...
}