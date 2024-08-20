#include "TranslationalJoint.h"
#include"Database.h"

//Variáveis globais
extern
Database db;

TranslationalJoint::TranslationalJoint()
{
	n_GL = 2;						//Dois graus de liberdade (esse vínculo possui 2 multiplicadores de lagrange)
	active_lambda = new int[n_GL];
	lambda = new double[n_GL];
	copy_lambda = new double[n_GL];
	GLs = new int[n_GL];
	node_A = 0;
	node_B = 0;
	rot_node = 0;
	cs = 0;

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

	alphaA = Matrix(3);
	alphaiA = Matrix(3);

	uA = Matrix(3);
	uB = Matrix(3);
	
	A = Matrix(3, 3);
	QA = Matrix(3, 3);
	
	ei1A = Matrix(3);
	ei2A = Matrix(3);

	//Matriz de rigidez tangente e vetor resíduo
	stiffness = new double*[11];
	for (int i = 0; i < 11; i++)
		stiffness[i] = new double[11];
	residual = new double[11];
	
	for (int i = 0; i < 11; i++)
	{
		residual[i] = 0.0;
		for (int j = 0; j < 11; j++)
			stiffness[i][j] = 0.0;
	}
}

//Zera matrizes e vetores
void TranslationalJoint::ClearContributions()
{
	for (int i = 0; i < 11; i++)
	{
		residual[i] = 0.0;
		for (int j = 0; j < 11; j++)
			stiffness[i][j] = 0.0;
	}
}

TranslationalJoint::~TranslationalJoint()
{
	delete[]active_lambda;
	delete[]GLs;
	delete[]lambda;
	delete[]copy_lambda;
	for (int i = 0; i < 11; i++)
		delete []stiffness[i];
	delete []stiffness;
	delete []residual;
}

//Leitura
bool TranslationalJoint::Read(FILE *f)
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
	if (!strcmp(s, "RotationNode"))
	{
		fscanf(f, "%s", s);
		rot_node = atoi(s);
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
void TranslationalJoint::Write(FILE *f)
{
	fprintf(f, "TranslationalJoint\t%d\tNodes\t%d\t%d\tRotationNode\t%d\tCS\t%d\n",
		number,
		node_A,
		node_B,
		rot_node,
		cs);
}

//Escreve no monitor do SpecialConstraint
void TranslationalJoint::WriteMonitor(FILE *f, bool first_record, double time)
{

}

void TranslationalJoint::WriteVTK_XMLRender(FILE *f)
{
	
}

//Checa inconsistências no SC para evitar erros de execução
bool TranslationalJoint::Check()
{
	if (node_A > db.number_nodes)
		return false;
	if (node_B > db.number_nodes)
		return false;
	if (rot_node > db.number_nodes)
		return false;
	if (cs > db.number_CS)
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

//Montagem dos resíduos e rigidez tangente
void TranslationalJoint::Mount()
{
	ClearContributions();
	//Montagem da rigidez tangente e resíduo
	if (active_lambda[0] == 1 && active_lambda[1] == 1)
	{
		for (int i = 0; i < 3; i++)
		{
			alphaA(i, 0) = db.nodes[rot_node - 1]->displacements[i + 3];		//vetor rotação (atual) do nó de ref para rotação
			alphaiA(i, 0) = db.nodes[rot_node - 1]->copy_coordinates[i + 3];	//vetor rotação acumulada (do início) do nó de ref para rotação
			uA(i, 0) = db.nodes[node_A - 1]->displacements[i];				//vetor de deslocamentos do nó A
			uB(i, 0) = db.nodes[node_B - 1]->displacements[i];				//vetor de deslocamentos do nó B
		}
		alpha_escalar_i = norm(alphaiA);
		A = skew(alphaiA);
		g = 4.0 / (4.0 + alpha_escalar_i*alpha_escalar_i);
		QA = I3 + g*(A + 0.5*(A*A));
		ei1A = QA*(*db.CS[cs - 1]->E1);	//Eixo e1 no início do incremento
		ei2A = QA*(*db.CS[cs - 1]->E2);	//Eixo e2 no início do incremento
		EvaluateTranslationalContribution(temp_v, residual, stiffness, uA.getMatrix(), uB.getMatrix(), alphaA.getMatrix(), ei1A.getMatrix(), ei2A.getMatrix(), lambda);
	}
}

//Preenche a contribuição do elemento nas matrizes globais
void TranslationalJoint::MountGlobal()
{
	//Variáveis temporárias para salvar a indexação global dos graus de liberdade a serem setados na matriz de rigidez global
	int GL_global_1 = 0;
	int GL_global_2 = 0;
	double anterior = 0;
	
	if (active_lambda[0] == 1 && active_lambda[0] == 1)
	{
		
		for (int i = 0; i < 11; i++)
		{
			//////////////MONTAGEM DO VETOR DE ESFORÇOS DESBALANCEADOS//////////////////
			//Toma do vetor de GL globais, a indexação de cada grau de liberdade global
			if (i < 3)//u_A -> i=0,1,2
				GL_global_1 = db.nodes[node_A - 1]->GLs[i];
			else
			{
				if (i < 6)//u_B -> i=3,4,5
					GL_global_1 = db.nodes[node_B - 1]->GLs[i-3];
				else
				{
					if (i < 9)//alpha_A -> i=6,7,8
						GL_global_1 = db.nodes[rot_node - 1]->GLs[i - 3];
					else
					{
						//lambda  -> i=9,10
						GL_global_1 = GLs[i - 9];
					}
					
				}
			}

			//Caso o grau de liberdade seja livre:
			if (GL_global_1 > 0)
			{
				anterior = db.global_P_A(GL_global_1 - 1, 0);
				db.global_P_A(GL_global_1 - 1, 0) = anterior + residual[i];
				anterior = db.global_I_A(GL_global_1 - 1, 0);
				db.global_I_A(GL_global_1 - 1, 0) = anterior + residual[i];
			}
			else
			{
				anterior = db.global_P_B(-GL_global_1 - 1, 0);
				db.global_P_B(-GL_global_1 - 1, 0) = anterior + residual[i];
			}
			for (int j = 0; j < 11; j++)
			{
				//////////////////////MONTAGEM DA MATRIZ DE RIGIDEZ/////////////////////////
				//Toma do vetor de GL globais, a indexação de cada grau de liberdade global
				if (j < 3)//u_A -> j=0,1,2
					GL_global_2 = db.nodes[node_A - 1]->GLs[j];
				else
				{
					if (j < 6)//u_B -> j=3,4,5
						GL_global_2 = db.nodes[node_B - 1]->GLs[j - 3];
					else
					{
						if (j < 9)//alpha_A -> j=6,7,8
							GL_global_2 = db.nodes[rot_node - 1]->GLs[j - 3];
						else
						{
							//lambda  -> j=9,10
							GL_global_2 = GLs[j - 9];
						}

					}
				}

				//Caso os graus de liberdade sejam ambos livres (Matriz Kaa)
				if (GL_global_1 > 0 && GL_global_2 > 0)
					db.global_stiffness_AA.setValue(GL_global_1 - 1, GL_global_2 - 1, stiffness[i][j]);
				//Caso os graus de liberdade sejam ambos fixos (Matriz Kbb)
				if (GL_global_1 < 0 && GL_global_2 < 0)
					db.global_stiffness_BB.setValue(-GL_global_1 - 1, -GL_global_2 - 1, stiffness[i][j]);
				//Caso os graus de liberdade sejam livre e fixo (Matriz Kab)
				if (GL_global_1 > 0 && GL_global_2 < 0)
					db.global_stiffness_AB.setValue(GL_global_1 - 1, -GL_global_2 - 1, stiffness[i][j]);
				//Caso os graus de liberdade sejam fixo e livre (Matriz Kba)
				if (GL_global_1 < 0 && GL_global_2 > 0)
					db.global_stiffness_BA.setValue(-GL_global_1 - 1, GL_global_2 - 1, stiffness[i][j]);
			}
		}
	}
}

void TranslationalJoint::ComputeInitialGuessDisplacements()
{

}

//Computa efeito das condições iniciais nos nós da restrição
void TranslationalJoint::ComputeVelAccel()
{
	if (bool_table.GetAt(db.current_solution_number - 1) == true)
	{
		//TODO
	}
}

//Pré-cálculo de variáveis que é feito uma única vez no início
void TranslationalJoint::PreCalc()
{
	//Not applicable
}

//Salvando variáveis da configuração convergida
void TranslationalJoint::SaveLagrange()
{
	for (int i = 0; i < n_GL; i++)
		copy_lambda[i] = lambda[i];
}

//Checa quais multiplicadores de lagrange serão ativados,de acordo com a ativação dos GLs dos nós dos quais a special constraint participa
void TranslationalJoint::ActivateDOFs()
{
	//Ativa GLs de translação dos nós A e B
	for (int i = 0; i < 3; i++)
		db.nodes[node_A - 1]->active_GL[i] = 1;
	for (int i = 0; i < 3; i++)
		db.nodes[node_B - 1]->active_GL[i] = 1;
	//Ativa GLs de rotação do rotation node
	for (int i = 0; i < 3; i++)
		db.nodes[rot_node - 1]->active_GL[i + 3] = 1;

	if (bool_table.GetAt(db.current_solution_number - 1))
	{
		//Ativa GLs de multiplicadores de Lagrange
		active_lambda[0] = 1;
		active_lambda[1] = 1;
	}
	else
	{
		active_lambda[0] = 0;
		active_lambda[1] = 0;
	}
}
//Calcula contribuições do resíduo e operador tangente - gerado no AceGen
void TranslationalJoint::EvaluateTranslationalContribution(double *v, double *residual, double **stiffness, double *uA, double *uB, double *alphaA, double *ei1A , double *ei2A, double *lambda)
{
	int i01; int i02;
	v[446] = ei2A[0] * lambda[1];
	v[445] = ei1A[0] * lambda[0];
	v[449] = v[445] + v[446];
	v[444] = ei2A[1] * lambda[1];
	v[443] = ei1A[1] * lambda[0];
	v[448] = v[443] + v[444];
	v[442] = ei2A[2] * lambda[1];
	v[441] = ei1A[2] * lambda[0];
	v[447] = v[441] + v[442];
	v[440] = uA[2] - uB[2];
	v[439] = uA[1] - uB[1];
	v[438] = uA[0] - uB[0];
	v[432] = Power(alphaA[2], 2);
	v[431] = 0.5e0*alphaA[2];
	v[430] = 2e0*alphaA[2];
	v[429] = Power(alphaA[1], 2);
	v[428] = 0.5e0*alphaA[1];
	v[427] = 2e0*alphaA[1];
	v[426] = Power(alphaA[0], 2);
	v[425] = 2e0*alphaA[0];
	v[424] = 0.5e0*alphaA[0];
	v[165] = alphaA[1] * v[424];
	v[203] = -v[426] - v[429];
	v[210] = -alphaA[2] + v[165];
	v[208] = alphaA[2] + v[165];
	v[172] = alphaA[2] * v[428];
	v[206] = -alphaA[0] + v[172];
	v[204] = alphaA[0] + v[172];
	v[170] = alphaA[2] * v[424];
	v[209] = alphaA[1] + v[170];
	v[205] = -alphaA[1] + v[170];
	v[215] = 4e0 + v[426] + v[429] + v[432];
	v[235] = 1e0 / Power(v[215], 2);
	v[433] = -4e0*v[235];
	v[238] = v[430] * v[433];
	v[437] = 0.5e0*v[238];
	v[283] = v[203] * v[437];
	v[237] = v[427] * v[433];
	v[434] = 0.5e0*v[237];
	v[236] = v[425] * v[433];
	v[436] = 0.5e0*v[236];
	v[211] = -v[429] - v[432];
	v[239] = v[211] * v[436];
	v[207] = -v[426] - v[432];
	v[261] = v[207] * v[434];
	v[159] = 4e0 / v[215];
	v[435] = -0.5e0*v[159];
	v[281] = v[427] * v[435];
	v[282] = v[281] + v[203] * v[434];
	v[279] = v[425] * v[435];
	v[280] = v[279] + v[203] * v[436];
	v[276] = v[159] + v[204] * v[236];
	v[274] = -v[159] + v[205] * v[237];
	v[264] = -v[159] + v[206] * v[236];
	v[262] = v[430] * v[435];
	v[263] = v[262] + v[207] * v[437];
	v[260] = v[279] + v[207] * v[436];
	v[259] = v[159] + v[208] * v[238];
	v[249] = v[159] + v[209] * v[237];
	v[247] = -(v[159] * v[431]);
	v[277] = v[204] * v[237] - v[247];
	v[288] = ei1A[0] * v[274] + ei1A[1] * v[277] + ei1A[2] * v[282];
	v[285] = ei2A[0] * v[274] + ei2A[1] * v[277] + ei2A[2] * v[282];
	v[291] = lambda[1] * v[285] + lambda[0] * v[288];
	v[273] = v[205] * v[236] - v[247];
	v[287] = ei1A[0] * v[273] + ei1A[1] * v[276] + ei1A[2] * v[280];
	v[284] = ei2A[0] * v[273] + ei2A[1] * v[276] + ei2A[2] * v[280];
	v[290] = lambda[1] * v[284] + lambda[0] * v[287];
	v[265] = v[206] * v[237] - v[247];
	v[248] = v[209] * v[236] - v[247];
	v[246] = -v[159] + v[210] * v[238];
	v[244] = -(v[159] * v[424]);
	v[275] = v[205] * v[238] - v[244];
	v[258] = v[208] * v[237] - v[244];
	v[271] = ei1A[0] * v[258] + ei1A[1] * v[261] + ei1A[2] * v[265];
	v[268] = ei2A[0] * v[258] + ei2A[1] * v[261] + ei2A[2] * v[265];
	v[294] = lambda[1] * v[268] + lambda[0] * v[271];
	v[250] = v[209] * v[238] - v[244];
	v[245] = v[210] * v[237] - v[244];
	v[242] = v[159] * v[428];
	v[278] = v[204] * v[238] + v[242];
	v[289] = ei1A[0] * v[275] + ei1A[1] * v[278] + ei1A[2] * v[283];
	v[286] = ei2A[0] * v[275] + ei2A[1] * v[278] + ei2A[2] * v[283];
	v[292] = lambda[1] * v[286] + lambda[0] * v[289];
	v[266] = v[206] * v[238] + v[242];
	v[272] = ei1A[0] * v[259] + ei1A[1] * v[263] + ei1A[2] * v[266];
	v[269] = ei2A[0] * v[259] + ei2A[1] * v[263] + ei2A[2] * v[266];
	v[295] = lambda[1] * v[269] + lambda[0] * v[272];
	v[257] = v[208] * v[236] + v[242];
	v[270] = ei1A[0] * v[257] + ei1A[1] * v[260] + ei1A[2] * v[264];
	v[267] = ei2A[0] * v[257] + ei2A[1] * v[260] + ei2A[2] * v[264];
	v[293] = lambda[1] * v[267] + lambda[0] * v[270];
	v[243] = v[210] * v[236] + v[242];
	v[254] = ei1A[0] * v[239] + ei1A[1] * v[243] + ei1A[2] * v[248];
	v[251] = ei2A[0] * v[239] + ei2A[1] * v[243] + ei2A[2] * v[248];
	v[296] = lambda[1] * v[251] + lambda[0] * v[254];
	v[241] = v[262] + v[211] * v[437];
	v[256] = ei1A[0] * v[241] + ei1A[1] * v[246] + ei1A[2] * v[250];
	v[253] = ei2A[0] * v[241] + ei2A[1] * v[246] + ei2A[2] * v[250];
	v[298] = lambda[1] * v[253] + lambda[0] * v[256];
	v[240] = v[281] + v[211] * v[434];
	v[255] = ei1A[0] * v[240] + ei1A[1] * v[245] + ei1A[2] * v[249];
	v[252] = ei2A[0] * v[240] + ei2A[1] * v[245] + ei2A[2] * v[249];
	v[297] = lambda[1] * v[252] + lambda[0] * v[255];
	v[162] = 1e0 - v[211] * v[435];
	v[163] = v[159] * v[210];
	v[164] = v[159] * v[209];
	v[189] = ei2A[0] * v[162] + ei2A[1] * v[163] + ei2A[2] * v[164];
	v[188] = ei1A[0] * v[162] + ei1A[1] * v[163] + ei1A[2] * v[164];
	v[166] = v[159] * v[208];
	v[168] = 1e0 - v[207] * v[435];
	v[169] = v[159] * v[206];
	v[186] = ei2A[0] * v[166] + ei2A[1] * v[168] + ei2A[2] * v[169];
	v[185] = ei1A[0] * v[166] + ei1A[1] * v[168] + ei1A[2] * v[169];
	v[171] = v[159] * v[205];
	v[173] = v[159] * v[204];
	v[174] = 1e0 - v[203] * v[435];
	v[183] = ei2A[0] * v[171] + ei2A[1] * v[173] + ei2A[2] * v[174];
	v[182] = ei1A[0] * v[171] + ei1A[1] * v[173] + ei1A[2] * v[174];
	v[184] = lambda[0] * v[182] + lambda[1] * v[183];
	v[187] = lambda[0] * v[185] + lambda[1] * v[186];
	v[190] = lambda[0] * v[188] + lambda[1] * v[189];
	v[191] = v[440] * v[447];
	v[450] = 0.5e0*v[191];
	v[303] = -(v[191] * v[436]);
	v[219] = v[191] * v[435];
	v[192] = v[440] * v[448];
	v[312] = v[192] * v[236];
	v[213] = v[159] * v[192];
	v[193] = v[440] * v[449];
	v[322] = v[193] * v[237];
	v[321] = v[193] * v[236];
	v[217] = v[159] * v[193];
	v[194] = v[439] * v[447];
	v[452] = v[192] + v[194];
	v[372] = v[237] * v[452];
	v[416] = v[372] * v[431];
	v[329] = v[194] * v[236];
	v[371] = v[312] + v[329];
	v[214] = v[159] * v[194];
	v[195] = v[439] * v[448];
	v[453] = 0.5e0*v[195];
	v[337] = -(v[195] * v[436]);
	v[224] = v[195] * v[435];
	v[196] = v[439] * v[449];
	v[344] = v[196] * v[236];
	v[222] = v[159] * v[196];
	v[197] = v[438] * v[447];
	v[457] = v[193] + v[197];
	v[352] = v[197] * v[237];
	v[351] = v[197] * v[236];
	v[376] = v[321] + v[351];
	v[417] = v[376] * v[431];
	v[218] = v[159] * v[197];
	v[198] = v[438] * v[448];
	v[459] = v[196] - v[198];
	v[454] = v[196] + v[198];
	v[358] = v[198] * v[236];
	v[381] = v[344] + v[358];
	v[410] = v[381] * v[428];
	v[223] = v[159] * v[198];
	v[199] = v[438] * v[449];
	v[451] = 0.5e0*v[199];
	v[456] = -v[451] - v[453];
	v[367] = -(v[199] * v[434]);
	v[366] = -(v[199] * v[436]);
	v[225] = v[199] * v[435];
	v[200] = v[213] + v[214];
	v[201] = v[217] + v[218];
	v[202] = v[222] + v[223];
	v[212] = v[192] * v[204] + v[193] * v[205] + v[194] * v[206] + v[196] * v[208] + v[197] * v[209] + v[198] * v[210]
		+ v[203] * v[450] + v[211] * v[451] + v[207] * v[453];
	v[395] = (8e0*v[212]) / Power(v[215], 3);
	v[396] = v[395] * v[427] + v[433] * (-v[193] + v[197] + v[427] * (-v[450] - v[451]) + v[431] * v[452] + v[424] * v[454]);
	v[458] = v[367] + v[396];
	v[393] = v[395] * v[425] + v[433] * (v[192] - v[194] + v[425] * (-v[450] - v[453]) + v[428] * v[454] + v[431] * v[457]);
	v[455] = v[337] + v[393];
	v[220] = v[212] * v[433];
	v[415] = v[220] + v[224] + v[225];
	v[409] = v[219] - v[224] + v[415];
	v[403] = v[224] - v[225] + v[409];
	residual[0] = v[190];
	residual[1] = v[187];
	residual[2] = v[184];
	residual[3] = -v[190];
	residual[4] = -v[187];
	residual[5] = -v[184];
	residual[6] = v[213] - v[214] + v[403] * v[425] + v[202] * v[428] + v[201] * v[431];
	residual[7] = -v[217] + v[218] + v[202] * v[424] + v[409] * v[427] + v[200] * v[431];
	residual[8] = v[222] - v[223] + v[201] * v[424] + v[200] * v[428] + v[415] * v[430];
	residual[9] = v[188] * v[438] + v[185] * v[439] + v[182] * v[440];
	residual[10] = v[189] * v[438] + v[186] * v[439] + v[183] * v[440];
	stiffness[0][0] = 0e0;
	stiffness[0][1] = 0e0;
	stiffness[0][2] = 0e0;
	stiffness[0][3] = 0e0;
	stiffness[0][4] = 0e0;
	stiffness[0][5] = 0e0;
	stiffness[0][6] = v[296];
	stiffness[0][7] = v[297];
	stiffness[0][8] = v[298];
	stiffness[0][9] = v[188];
	stiffness[0][10] = v[189];
	stiffness[1][1] = 0e0;
	stiffness[1][2] = 0e0;
	stiffness[1][3] = 0e0;
	stiffness[1][4] = 0e0;
	stiffness[1][5] = 0e0;
	stiffness[1][6] = v[293];
	stiffness[1][7] = v[294];
	stiffness[1][8] = v[295];
	stiffness[1][9] = v[185];
	stiffness[1][10] = v[186];
	stiffness[2][2] = 0e0;
	stiffness[2][3] = 0e0;
	stiffness[2][4] = 0e0;
	stiffness[2][5] = 0e0;
	stiffness[2][6] = v[290];
	stiffness[2][7] = v[291];
	stiffness[2][8] = v[292];
	stiffness[2][9] = v[182];
	stiffness[2][10] = v[183];
	stiffness[3][3] = 0e0;
	stiffness[3][4] = 0e0;
	stiffness[3][5] = 0e0;
	stiffness[3][6] = -v[296];
	stiffness[3][7] = -v[297];
	stiffness[3][8] = -v[298];
	stiffness[3][9] = -v[188];
	stiffness[3][10] = -v[189];
	stiffness[4][4] = 0e0;
	stiffness[4][5] = 0e0;
	stiffness[4][6] = -v[293];
	stiffness[4][7] = -v[294];
	stiffness[4][8] = -v[295];
	stiffness[4][9] = -v[185];
	stiffness[4][10] = -v[186];
	stiffness[5][5] = 0e0;
	stiffness[5][6] = -v[290];
	stiffness[5][7] = -v[291];
	stiffness[5][8] = -v[292];
	stiffness[5][9] = -v[182];
	stiffness[5][10] = -v[183];
	stiffness[6][6] = v[312] - v[329] + 2e0*v[403] + v[410] + v[417] + v[425] * (v[303] + v[455]);
	stiffness[6][7] = 0.5e0*v[202] - v[321] + v[351] + v[381] * v[424] + (v[303] + v[366] + v[393])*v[427]
		+ v[371] * v[431];
	stiffness[6][8] = 0.5e0*v[201] + v[344] - v[358] + v[376] * v[424] + v[371] * v[428] + v[430] * (v[366] + v[455]);
	stiffness[6][9] = v[254] * v[438] + v[270] * v[439] + v[287] * v[440];
	stiffness[6][10] = v[251] * v[438] + v[267] * v[439] + v[284] * v[440];
	stiffness[7][7] = -v[322] + v[352] + 2e0*v[409] + v[410] + v[416] + v[427] * (-(v[237] * v[450]) + v[458]);
	stiffness[7][8] = 0.5e0*v[200] + (v[322] + v[352])*v[424] + v[372] * v[428] + v[430] * (-(v[237] * v[453]) + v[458])
		+ v[237] * v[459];
	stiffness[7][9] = v[255] * v[438] + v[271] * v[439] + v[288] * v[440];
	stiffness[7][10] = v[252] * v[438] + v[268] * v[439] + v[285] * v[440];
	stiffness[8][8] = 2e0*v[415] + v[416] + v[417] + v[238] * v[459] + v[430] * (v[395] * v[430] + v[238] * v[456] + v[433] *
		(v[428] * v[452] + v[430] * v[456] + v[424] * v[457] + v[459]));
	stiffness[8][9] = v[256] * v[438] + v[272] * v[439] + v[289] * v[440];
	stiffness[8][10] = v[253] * v[438] + v[269] * v[439] + v[286] * v[440];
	stiffness[9][9] = 0e0;
	stiffness[9][10] = 0e0;
	stiffness[10][10] = 0e0;
	for (i01 = 1; i01<11; i01++){
		for (i02 = 0; i02<i01; i02++){
			stiffness[i01][i02] = stiffness[i02][i01];
		}
	};
}




