#include "LRLR.h"
#include <typeinfo>

#include "Contact.h"
#include "Pipe_1.h"
#include "Beam_1.h"
#include "SecTube.h"
#include "LagrangeSaveLRAS.h"
#include "LineRegion.h"
#include "Node.h"
#include "Matrix.h"

#define PI 3.1415926535897932384626433832795
#include"Database.h"
//Variaveis globais
extern
Database db;

LRLR::LRLR()
{
	tol_NR = 1e-16;		//Tolerência NR
	max_it = 20;		//Numero maximo de iterações
	tol_ortho = 1e-16;	//Tolerancia de erro a ortogonalidade
	type_name = new char[20];//Nome do tipo do contato
	sprintf(type_name, "LRLR");
	number = 0;
	n_LR1 = 0;
	n_LR2 = 0;
	mu = 0;
	ept = 0;
	epn = 0;
	pinball = 0;
	sum_Fx = 0;
	sum_Fy = 0;
	sum_Fz = 0;
	//Variaveis internas
	temp_element1 = 0;
	temp_element2 = 0;
	temp_node = 0;
	mean_position1 = new Matrix(3,1);
	mean_position2 = new Matrix(3,1);
	gte = new Matrix(3,1);
	node_1_1 = new Matrix(3,1);
	node_3_1 = new Matrix(3,1);
	node_1_2 = new Matrix(3,1);
	node_3_2 = new Matrix(3,1);
	N1 = new Matrix(3, 6);
	dN1 = new Matrix(3, 6);
	ddN1 = new Matrix(3, 6);
	N2 = new Matrix(3, 6);
	dN2 = new Matrix(3, 6);
	ddN2 = new Matrix(3, 6);
	x1 = new Matrix(6, 1);
	x2 = new Matrix(6, 1);
	z1 = new Matrix(3,1);
	z2 = new Matrix(3,1);
	t1ext = new Matrix(6,1);
	t2ext = new Matrix(6,1);
	I361 = new Matrix(6,6);
	I362 = new Matrix(6,6);
	Ntio = new Matrix(3,12);
	N_ext = new Matrix(6,12);
	dN_ext = new Matrix(6,12);
	G1 = new Matrix(12,1);
	G1b = new Matrix(12,1);
	G2 = new Matrix(12,1);
	G2b = new Matrix(12,1);
	G3 = new Matrix(12,12);
	G3b = new Matrix(12,12);
	aux1 = new Matrix(6,1);
	aux2 = new Matrix(6,1);
	G4 = new Matrix(12,1);
	G4b = new Matrix(12,1);
	aux3 = new Matrix(6,1);
	G5 = new Matrix(1,12);
	aux4 = new Matrix(6,1);
	aux5 = new Matrix(6, 1);
	aux6 = new Matrix(6, 1);
	aux7 = new Matrix(6, 1);
	Q1 = new Matrix(6, 6);
	Q2 = new Matrix(6, 6);
	Q3 = new Matrix(6, 6);
	Q4 = new Matrix(6, 6);
	G6 = new Matrix(6, 6);
	G6b = new Matrix(12,12);
	G7 = new Matrix(12,12);
	G7b = new Matrix(12,12);
	I3 = new Matrix(3,3);
	n = new Matrix(3,1);
	dz1 = new Matrix(3,1);
	dz2 = new Matrix(3,1);
	ddz1 = new Matrix(3,1);
	ddz2 = new Matrix(3,1);
	d1 = new Matrix(1,12);
	d2 = new Matrix(1,12);
	E1 = new Matrix(6,12);
	E2 = new Matrix(6,12);
	ST = new Matrix(12);
	STb = new Matrix(12);
	R = new Matrix(2, 1);						
	A = new Matrix(2, 2);						
	B = new Matrix(2, 6);
	C = new Matrix(2, 6);
	D = new Matrix(2, 12);
	E = new Matrix(12, 12);
	F = new Matrix(12, 12);
	G = new Matrix(12, 12);
	t1 = new Matrix(3, 1);
	t2 = new Matrix(3, 1);
	delta_csi = new Matrix(2, 1);
	//Abaixo as variaveis que dependem do numero de elementos para serem alocadas - Alocação feita na função Alloc, quando chamada durante o PreCalc
	typeOK1 = NULL;
	typeOK2 = NULL;
	activate = NULL;
	L1 = NULL;
	r1 = NULL;
	L2 = NULL;
	r2 = NULL;
	last_dif_pos = NULL;
	c_stiffness = NULL;
	c_loading = NULL;
	last_g_n = NULL;
	csi_1_0 = NULL;
	csi_2_0 = NULL;
	copy_csi_1_converged = NULL;
	copy_csi_2_converged = NULL;
	gt1s = NULL;
	gt2s = NULL;
	g_t1s_temp = NULL;
	g_t2s_temp = NULL;
	return_value = NULL;
	last_return_value = NULL;
	g_n = NULL;
	first_mount = NULL;
	z2z1 = NULL;
	csi_1 = NULL;
	csi_2 = NULL;
	flag_cross = NULL;
	sticking = NULL;
	factor_1 = NULL;
	factor_2 = NULL;

}

LRLR::~LRLR()
{
	delete[] type_name;
	delete mean_position1;
	delete mean_position2;
	delete gte;
	delete node_1_1;
	delete node_3_1;
	delete node_1_2;
	delete node_3_2;
	delete N1;
	delete dN1;
	delete ddN1;
	delete N2;
	delete dN2;
	delete ddN2;
	delete x1;
	delete x2;
	delete z1;
	delete z2;
	delete t1ext;
	delete t2ext;
	delete I361;
	delete I362;
	delete Ntio;
	delete G1;
	delete G1b;
	delete G2;
	delete G2b;
	delete G3;
	delete G3b;
	delete aux1;
	delete aux2;
	delete G4;
	delete G4b;
	delete aux3;
	delete G5;
	delete aux4;
	delete aux5;
	delete aux6;
	delete aux7;
	delete Q1;
	delete Q2;
	delete Q3;
	delete Q4;
	delete G6;
	delete G6b;
	delete G7;
	delete G7b;
	delete I3;
	delete n;
	delete dz1;
	delete dz2 ;
	delete ddz1;
	delete ddz2;
	delete N_ext;
	delete dN_ext;
	delete d1;
	delete d2;
	delete E1;
	delete E2;
	delete ST;
	delete STb;
	delete R;
	delete A;
	delete B;
	delete C;
	delete D;
	delete E;
	delete F;
	delete G;
	delete t1;
	delete t2;
	delete delta_csi;
	if (activate != NULL)
	{
		for (int i = 0; i < n_elements1; i++)
			delete activate[i];
		delete[] activate;
	}
	if (typeOK1 != NULL)
		delete[] typeOK1;
	if (L1 != NULL)
		delete[] L1;
	if (r1 != NULL)
		delete[] r1;
	if (typeOK2 != NULL)
		delete[] typeOK2;
	if (L2 != NULL)
		delete[] L2;
	if (r2 != NULL)
		delete[] r2;
	if (first_mount != NULL)
	{
		for (int i = 0; i < n_elements1; i++)
			delete first_mount[i];
		delete[] first_mount;
	}
	if (last_g_n != NULL)
	{
		for (int i = 0; i < n_elements1; i++)
			delete last_g_n[i];
		delete[] last_g_n;
	}
	if (csi_1_0 != NULL)
	{
		for (int i = 0; i < n_elements1; i++)
			delete csi_1_0[i];
		delete[] csi_1_0;
	}
	if (csi_2_0 != NULL)
	{
		for (int i = 0; i < n_elements1; i++)
			delete csi_2_0[i];
		delete[] csi_2_0;
	}
	if (copy_csi_1_converged != NULL)
	{
		for (int i = 0; i < n_elements1; i++)
			delete copy_csi_1_converged[i];
		delete[] copy_csi_1_converged;
	}
	if (copy_csi_2_converged != NULL)
	{
		for (int i = 0; i < n_elements1; i++)
			delete copy_csi_2_converged[i];
		delete[] copy_csi_2_converged;
	}
	if (gt1s != NULL)
	{
		for (int i = 0; i < n_elements1; i++)
			delete gt1s[i];
		delete[] gt1s;
	}
	if (gt2s != NULL)
	{
		for (int i = 0; i < n_elements1; i++)
			delete gt2s[i];
		delete[] gt2s;
	}
	if (g_t1s_temp != NULL)
	{
		for (int i = 0; i < n_elements1; i++)
			delete g_t1s_temp[i];
		delete[] g_t1s_temp;
	}
	if (g_t2s_temp != NULL)
	{
		for (int i = 0; i < n_elements1; i++)
			delete g_t2s_temp[i];
		delete[] g_t2s_temp;
	}
	if (return_value != NULL)
	{
		for (int i = 0; i < n_elements1; i++)
			delete return_value[i];
		delete[] return_value;
	}
	if (last_return_value != NULL)
	{
		for (int i = 0; i < n_elements1; i++)
			delete last_return_value[i];
		delete[] last_return_value;
	}
	if (g_n != NULL)
	{
		for (int i = 0; i < n_elements1; i++)
			delete g_n[i];
		delete[] g_n;
	}
	if (csi_1 != NULL)
	{
		for (int i = 0; i < n_elements1; i++)
			delete csi_1[i];
		delete[] csi_1;
	}
	if (csi_2 != NULL)
	{
		for (int i = 0; i < n_elements1; i++)
			delete csi_2[i];
		delete[] csi_2;
	}

	if (last_dif_pos != NULL)
	{
		for (int i = 0; i < n_elements1; i++)
		{
			for (int j = 0; j < n_elements2; j++)
				delete last_dif_pos[i][j];
			delete[] last_dif_pos[i];
		}
		delete[] last_dif_pos;
	}
	if (c_stiffness != NULL)
	{
		for (int i = 0; i < n_elements1; i++)
		{
			for (int j = 0; j < n_elements2; j++)
				delete c_stiffness[i][j];
			delete[] c_stiffness[i];
		}
		delete[] c_stiffness;
	}
	if (c_loading != NULL)
	{
		for (int i = 0; i < n_elements1; i++)
		{
			for (int j = 0; j < n_elements2; j++)
				delete c_loading[i][j];
			delete[] c_loading[i];
		}
		delete[] c_loading;
	}
	if (z2z1 != NULL)
	{
		for (int i = 0; i < n_elements1; i++)
		{
			for (int j = 0; j < n_elements2; j++)
				delete z2z1[i][j];
			delete[] z2z1[i];
		}
		delete[] z2z1;
	}
	if (flag_cross != NULL)
	{
		for (int i = 0; i < n_elements1; i++)
			delete flag_cross[i];
		delete[] flag_cross;
	}
	if (sticking != NULL)
	{
		for (int i = 0; i < n_elements1; i++)
			delete sticking[i];
		delete[] sticking;
	}
	if (factor_1 != NULL)
	{
		for (int i = 0; i < n_elements1; i++)
			delete factor_1[i];
		delete[] factor_1;
	}
	if (factor_2 != NULL)
	{
		for (int i = 0; i < n_elements1; i++)
			delete factor_2[i];
		delete[] factor_2;
	}
}

bool LRLR::Read(FILE *f)
{
	char s[1000];
	fscanf(f, "%s", s);
	number = atoi(s);
	fscanf(f, "%s", s);
	if (!strcmp(s, "LR"))
	{
		fscanf(f, "%s", s);
		n_LR1 = atoi(s);
	}
	else
		return false;
	fscanf(f, "%s", s);
	if (!strcmp(s, "LR"))
	{
		fscanf(f, "%s", s);
		n_LR2 = atoi(s);
	}
	else
		return false;
	fscanf(f, "%s", s);
	if (!strcmp(s, "MU"))
	{
		fscanf(f, "%s", s);
		mu = atof(s);
	}
	else
		return false;
	fscanf(f, "%s", s);
	if (!strcmp(s, "EPN"))
	{
		fscanf(f, "%s", s);
		epn = atof(s);
	}
	else
		return false;
	fscanf(f, "%s", s);
	if (!strcmp(s, "EPT"))
	{
		fscanf(f, "%s", s);
		ept = atof(s);
	}
	else
		return false;
	fscanf(f, "%s", s);
	if (!strcmp(s, "Pinball"))
	{
		fscanf(f, "%s", s);
		pinball = atof(s);
	}
	else
		return false;
	fscanf(f, "%s", s);
	if (!strcmp(s, "TolOrtho"))
	{
		fscanf(f, "%s", s);
		tol_ortho = atof(s);
	}
	else
		return false;
	return true;
}

void LRLR::Write(FILE *f)
{
	fprintf(f, "LRLR\t%d\tLR\t%d\tLR\t%d\tMU\t%.6e\tEPN\t%.6e\tEPT\t%.6e\tPinball\t%.6e\tTolOrtho\t%.6e\n",
		number, n_LR1, n_LR2, mu, epn, ept, pinball,tol_ortho);
}

//Checa inconsistências no elemento para evitar erros de execução
bool LRLR::Check()
{
	return true;
}

void LRLR::WriteVTK_XMLRender(FILE *f)
{

}

void LRLR::WriteVTK_XMLForces(FILE *f)
{
	//TODO
}

//Escreve arquivo de resultados
void LRLR::WriteResults(FILE *f)
{
	
}

//Escreve no monitor do contato
void LRLR::WriteMonitor(FILE *f, bool first_record, double time)
{
	//Cabeçalho
	if (first_record == true)
		fprintf(f, "TIME\tFx\tFy\tFz\n");
	fprintf(f, "%.6e\t", time);	//TIME
	//Soma vetorial dos esforços em todos os nós do par de contato
	sum_Fx = 0;
	sum_Fy = 0;
	sum_Fz = 0;
	for (int i = 0; i < n_elements1; i++)
	{
		for (int j = 0; j < n_elements2; j++)
		{
			sum_Fx += (*c_loading[i][j])(0, 0) + (*c_loading[i][j])(3, 0);
			sum_Fy += (*c_loading[i][j])(1, 0) + (*c_loading[i][j])(4, 0);
			sum_Fz += (*c_loading[i][j])(2, 0) + (*c_loading[i][j])(5, 0);
		}
	}
	fprintf(f, "%.6e\t%.6e\t%.6e\n", sum_Fx, sum_Fy, sum_Fz);
}

//Salva variaveis para descrição lagrangiana atualizada
void LRLR::SaveLagrange()
{
	for (int i = 0; i < n_elements1; i++)
	{
		for (int j = 0; j < n_elements2; j++)
		{
			//Salva valores para chutes iniciais do próximo passo
			copy_csi_1_converged[i][j] = csi_1[i][j];
			copy_csi_2_converged[i][j] = csi_2[i][j];
			last_g_n[i][j] = g_n[i][j];
			last_return_value[i][j] = return_value[i][j];
			//Atualização em caso de slipping
			gt1s[i][j] = g_t1s_temp[i][j];
			gt2s[i][j] = g_t2s_temp[i][j];
			(*last_dif_pos[i][j]) = (*z2z1[i][j]);
			//Imprimindo variaveis de controle
			if (activate[i][j] == 1)
			{
				if (return_value[i][j] == 0 && g_n[i][j] < 0)
					PlotContactStatus(i, j);
			}			
		}
	}
}
//Plota caracteristicas do contato, independente de ter convergido ou não
void LRLR::PlotContactStatus(int i, int j)
{
	db.myprintf("\nLRLR Contact %d\n", number);
	temp_element1 = db.line_regions[n_LR1 - 1]->elements[i];
	temp_element2 = db.line_regions[n_LR2 - 1]->elements[j];
	db.myprintf("Element 1: %d  csi_1 %lf\tElement 2: %d  csi_2 %lf\n", temp_element1, csi_1[i][j], temp_element2, csi_2[i][j]);
	db.myprintf("gn %lf\tgt1s %lf\tgt2s %lf\tRet Value %d\tsticking %d\n", g_n[i][j], gt1s[i][j], gt2s[i][j], return_value[i][j], (int)sticking[i][j]);
	db.myprintf("csi_1_0 %lf csi_2_0 %lf\n\n", csi_1_0[i][j],csi_2_0[i][j]);
}

void LRLR::Mount()
{
	for (int i = 0; i < n_elements1; i++)
	{
		for (int j = 0; j < n_elements2; j++)
		{
			zeros(c_stiffness[i][j]);
			zeros(c_loading[i][j]);
			//Se o contato esta ativo (near) - realiza a montagem
			if (activate[i][j] == 1)
			{
				temp_element1 = db.line_regions[n_LR1 - 1]->elements[i];
				temp_element2 = db.line_regions[n_LR2 - 1]->elements[j];
				FillNodes(temp_element1, temp_element2);
				//Altera csi_1, csi_2, N1, N2, dN1, dN2, ddN1 e ddN2 para o ponto do possivel contato - necessario para construção de funções gap, etc.
				return_value[i][j] = FindMinimumParameters(i, j);
				//Retorno 0 - localizou ponto de minima distancia com NR
				//Retorno 1 - Problemas de convergência NR ou ponto esta fora dos dominios das barras com tolerancia a ortogonalidade não obedecida
				*z1 = (*N1)*(*x1);	//z1
				*z2 = (*N2)*(*x2);	//z2
				*z2z1[i][j] = *z2 - *z1;
				//Calculo do gap normal
				g_n[i][j] = norm(*z2z1[i][j]) - (r1[i] + r2[j]);
				//Caso seja a primeira montagem, realiza pre-calculo da distancia entre os pontos (para ter referência para verificar cruzamento de vigas)
				if (first_mount[i][j] == true)
				{
					*last_dif_pos[i][j] = *z2z1[i][j];
					first_mount[i][j] = false;
				}
				//Verifica se mudou status em relação a ultima configuração convergida - o resultado desse algoritmo e a variavel booleana flag_cross
				dot_test = dot(*z2z1[i][j], *last_dif_pos[i][j]);
				if (last_g_n[i][j] > 0.0)
				{
					flag_cross[i][j] = false;
				}
				else
				{
					if (dot_test > 0)		//Não houve cruzamento - produto escalar positivo
						flag_cross[i][j] = false;
					else
					{
						if (dot_test == 0)	//Ha um ponto em comum - um ponto em comum entre as vigas - aqui ha falha se as vigas girarem 90 deg em um passo (muito raro)
							flag_cross[i][j] = true;
						else				//Muda o status - houve cruzamento de vigas
							flag_cross[i][j] = true;
					}
				}
				
				//Se o flag_cross e true - modifica o calculo do gap normal
				if (flag_cross[i][j] == true)
				{
					g_n[i][j] = -(norm(*z2z1[i][j]) + (r1[i] + r2[j]));
				}

				if (return_value[i][j] == 0)	//Relações de ortogonalidade obedecidas
				{
					//Se houver contato, então realiza a contribuição
					if (g_n[i][j] <= 0)
					{
						//Calculo das direções tangenciais dos elementos e comprimentos dos elementos
						CalculateLengthsAndTangents();
						//VERIFICAÇÃO DE DOMiNIO DAS BARRAS
						//double d = mu*abs(g_n[i][j])*epn / ept;
						//d = 0.001;
						////Elemento 1
						////Ponto no interior do dominio
						//if (abs(csi_1[i][j]) < 1.0 - 2.0*d / l1)
						//{
						//	factor_1 = 1.0;
						//}
						//else
						//{
						//	//Ponto fora do dominio
						//	if (abs(csi_1[i][j]) > 1.0 + 2.0*d / l1)
						//	{
						//		factor_1 = 0.0;
						//	}
						//	else//Zona de transição
						//	{
						//		factor_1 = 0.5*(1 - 0.5*l1*abs(1 - csi_1[i][j]) / d);
						//	}
						//}
						////Elemento 2
						////Ponto no interior do dominio
						//if (abs(csi_2[i][j]) < 1.0 - 2.0*d / l2)
						//{
						//	factor_2 = 1.0;
						//}
						//else
						//{
						//	//Ponto fora do dominio
						//	if (abs(csi_2[i][j]) > 1.0 + 2.0*d / l2)
						//	{
						//		factor_2 = 0.0;
						//	}
						//	else//Zona de transição
						//	{
						//		factor_2 = 0.5*(1 - 0.5*l2*abs(1 - csi_2[i][j]) / d);
						//	}
						//}
						//if (factor_1 == 0 || factor_2 == 0)
						//{
						//	factor_1 = 0.0;
						//	factor_2 = 0.0;
						//}
						
						//Direção Normal
						*n = (1.0 / norm(*z2z1[i][j]))*(*z2z1[i][j]);
						*dz1 = (*dN1)*(*x1);	//dz1
						*dz2 = (*dN2)*(*x2);	//dz2
						*ddz1 = (*ddN1)*(*x1);	//ddz1
						*ddz2 = (*ddN2)*(*x2);	//ddz2
						//Construção de matrizes auxiliares do contato
						//A matriz [A] ja esta construida - e a própria matriz Jacobiana do algoritmo de busca - NR
						//Matriz [B]
						(*B)(0, 0) = (*dz1)(0, 0);
						(*B)(0, 1) = (*dz1)(1, 0);
						(*B)(0, 2) = (*dz1)(2, 0);
						(*B)(0, 3) = -(*dz1)(0, 0);
						(*B)(0, 4) = -(*dz1)(1, 0);
						(*B)(0, 5) = -(*dz1)(2, 0);
						(*B)(1, 0) = (*dz2)(0, 0);
						(*B)(1, 1) = (*dz2)(1, 0);
						(*B)(1, 2) = (*dz2)(2, 0);
						(*B)(1, 3) = -(*dz2)(0, 0);
						(*B)(1, 4) = -(*dz2)(1, 0);
						(*B)(1, 5) = -(*dz2)(2, 0);
						//Matriz [C]
						(*C)(0, 0) = -(*z2z1[i][j])(0, 0);
						(*C)(0, 1) = -(*z2z1[i][j])(1, 0);
						(*C)(0, 2) = -(*z2z1[i][j])(2, 0);
						(*C)(1, 3) = -(*z2z1[i][j])(0, 0);
						(*C)(1, 4) = -(*z2z1[i][j])(1, 0);
						(*C)(1, 5) = -(*z2z1[i][j])(2, 0);
						//Matriz [D]
						for (int ii = 0; ii<3; ii++)
						{
							for (int jj = 0; jj<6; jj++)
							{
								(*N_ext)(ii, jj) = (*N1)(ii, jj);
								(*N_ext)(ii + 3, jj + 6) = (*N2)(ii, jj);
								(*dN_ext)(ii, jj) = (*dN1)(ii, jj);
								(*dN_ext)(ii + 3, jj + 6) = (*dN2)(ii, jj);
							}
						}
						(*D) = invert2x2(*A)*((*B)*(*N_ext) + (*C)*(*dN_ext));
						//Vetores [d1] e [d2]
						for (int ii = 0; ii<12; ii++)
						{
							(*d1)(0, ii) = (*D)(0, ii);
							(*d2)(0, ii) = (*D)(1, ii);
						}
						//Matriz [E]
						*E1 = -1.0*transp(*dN1)*((*n)*(*d1));
						*E2 = +1.0*transp(*dN2)*((*n)*(*d2));
						for (int ii = 0; ii<6; ii++)
						{
							for (int jj = 0; jj<12; jj++)
							{
								(*E)(ii, jj) = (*E1)(ii, jj);
								(*E)(ii + 6, jj) = (*E2)(ii, jj);
							}
						}
						//Matriz [F] - no caso de elementos com interpolação linear essa contribuição sera nula
						(*F) = transp(*d2)*(transp(*ddz2)*((*n)*(*d2))) - 1.0*(transp(*d1)*(transp(*n)*((*ddz1)*(*d1))));
						//Matriz [G]
						d = norm(*z2z1[i][j]);
						for (int ii = 0; ii<3; ii++)
						{
							for (int jj = 0; jj<6; jj++)
							{
								(*Ntio)(ii, jj) = -1.0*(*N1)(ii, jj);
								(*Ntio)(ii, jj + 6) = +1.0*(*N2)(ii, jj);
							}
						}
						(*G) = (1.0 / d)*((transp(*Ntio) + transp(*d2)*transp(*dz2) - transp(*d1)*transp(*dz1))*(*I3 - (*n)*transp(*n))*(*Ntio + (*dz2)*(*d2) - (*dz1)*(*d1)));
						//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
						//Direção tangencial - implementação do atrito com base no paper Zavarise e Wriggers (2000)
						//Verifica se csi0 deve ou não ser atualizado (se o contato esta começando nesse passo ou não)
						if (last_g_n[i][j] >= 0)
						{
							csi_1_0[i][j] = csi_1[i][j];
							csi_2_0[i][j] = csi_2[i][j];
							//Zera possivel acumulo de deslocamento anterior - o contato esta reiniciando - não ha deslizamento acumulado
							gt1s[i][j] = 0.0;
							gt2s[i][j] = 0.0;
						}
						else
						{
							//Se o ponto caiu fora do dominio de pelo menos uma das barras no ultimo passo convergido (pode ser indicativo de deslizamento de elemento vizinho)
							if (last_return_value[i][j] == 3)
							{
								csi_1_0[i][j] = csi_1[i][j];
								csi_2_0[i][j] = csi_2[i][j];
								//Zera possivel acumulo de deslocamento anterior - o contato esta reiniciando - não ha deslizamento acumulado
								gt1s[i][j] = 0.0;
								gt2s[i][j] = 0.0;

								////Calculo da magnitude do deslizamento entre o ponto de fora do dominio e o ponto atual (no dominio)
								//double factor = 
								;
								//double gt1_pred = 0.5*(csi_1[i][j] - copy_csi_1_converged[i][j])*l1;
								//double gt2_pred = 0.5*(csi_2[i][j] - copy_csi_2_converged[i][j])*l2;
								//Matrix g_t_pred = gt1_pred*(*t1) + gt2_pred*(*t2);
								////Se esta delizando consideravelmente - calcula csi_1_0 e csi_2_0 ficticios
								//if (norm(g_t_pred)*ept >= factor*abs(mu*epn*g_n[i][j]))
								//{
								//	t1t2 = dot(*t1, *t2);
								//	tt1 = dot((1.0 / norm(g_t_pred))*(g_t_pred), *t1);
								//	tt2 = dot((1.0 / norm(g_t_pred))*(g_t_pred), *t2);
								//	m1 = ((tt2*t1t2 - tt1) / (t1t2*t1t2 - 1));
								//	m2 = ((tt1*t1t2 - tt2) / (t1t2*t1t2 - 1));
								//	csi_1_0[i][j] = csi_1[i][j] - m1*2.0*abs(mu*epn*g_n[i][j]) / (ept*l1);
								//	csi_2_0[i][j] = csi_2[i][j] - m2*2.0*abs(mu*epn*g_n[i][j]) / (ept*l2);
								//}
								//else
								//{
								//	csi_1_0[i][j] = 0.0;
								//	csi_2_0[i][j] = 0.0;
								//}
								////Zera possivel acumulo de deslocamento anterior - o contato esta reiniciando - não ha deslizamento acumulado
								//gt1s[i][j] = 0.0;
								//gt2s[i][j] = 0.0;
							}
						}
						
						//Calculo do gap tangencial para o elemento 1 e para o elemento 2
						gt1 = 0.5*(csi_1[i][j] - csi_1_0[i][j])*l1;
						gt2 = 0.5*(csi_2[i][j] - csi_2_0[i][j])*l2;
						//Gaps tangenciais elasticos (sticking)
						gte1 = gt1 - gt1s[i][j];
						gte2 = gt2 - gt2s[i][j];
						//Força de atrito maxima
						Fat_max = abs(mu*g_n[i][j]*epn);
						sticking[i][j] = true;
						//Gaps tangenciais vetoriais
						*gte = gte1*(*t1) + gte2*(*t2);	//elastico
						//Multiplicadores (para decomposição na base não ortho formada pelas barras)
						t1t2 = dot(*t1, *t2);
						if (norm(*gte) != 0.0)
						{
							tt1 = dot((1.0 / norm(*gte))*(*gte), *t1);
							tt2 = dot((1.0 / norm(*gte))*(*gte), *t2);
						}
						else
						{
							tt1 = 0.0;
							tt2 = 0.0;
						}
						
						m1 = ((tt2*t1t2 - tt1) / (t1t2*t1t2 - 1));
						m2 = ((tt1*t1t2 - tt2) / (t1t2*t1t2 - 1));
						Fat_try = ept*norm(*gte);
						//Magnitude da força de atrito que de fato ocorrera, após as devidas verificações a serem processadas a seguir:
						//Verificação de sticking ou slipping
						////////////////////////////Sliding//////////////////////////////
						if (Fat_try >= Fat_max)
						{
							sticking[i][j] = false;
							Fat = Fat_max;
							//Calcula o deslizamento
							if (ept != 0.0)
								delta_lambda = (Fat_try - Fat_max) / ept;	//Magnitude da atualização de deslizamento
							else
								delta_lambda = 0.0;
							//Atualiza o gap de deslizamento (mas só de fato sera computado no caso de convergência do NR)
							if (norm(*gte) != 0.0)
							{
								if (delta_lambda != 0.0)
								{
									delta_lambda_1 = delta_lambda*m1;
									delta_lambda_2 = delta_lambda*m2;
								}
								else //Frictionless (ept = 0)
								{
									delta_lambda_1 = gte1;
									delta_lambda_2 = gte2;
								}
								g_t1s_temp[i][j] = gt1s[i][j] + delta_lambda_1;
								g_t2s_temp[i][j] = gt2s[i][j] + delta_lambda_2;
							}
							else
							{
								g_t1s_temp[i][j] = gt1s[i][j];
								g_t2s_temp[i][j] = gt2s[i][j];
							}
						}
						else	//Sticking
						{
							Fat = Fat_try;
							g_t1s_temp[i][j] = gt1s[i][j];
							g_t2s_temp[i][j] = gt2s[i][j];
						}
						Fat1 = 0;
						Fat2 = 0;
						if (norm(*gte) != 0.0)
						{
							Fat1 = Fat*m1;
							Fat2 = Fat*m2;
							//printf("\nFat %lf\n", norm(Fat1*(*t1) + Fat2*(*t2)));
						}
						//Nesse momento ja se sabe se stick ou slip. Agora calculo das contribuiçoes para forma fraca e operador tangente
						for (int ii = 0; ii<3; ii++)
						{
							(*t1ext)(ii, 0) = (*t1)(ii, 0);
							(*t2ext)(ii + 3, 0) = (*t2)(ii, 0);
						}
						//Matrix auxiliar para o calculo de G3 - identidade de ordem 3 com zeros no restante - I361 e I362
						for (int ii = 0; ii<3; ii++)
						{
							(*I361)(ii, ii) = 1.0;
							(*I362)(ii + 3, ii + 3) = 1.0;
						}
						//Matriz [G1]
						*G1 = 2.0*transp(*dN_ext)*(*t1ext);
						*G1b = 2.0*transp(*dN_ext)*(*t2ext);
						//Matriz [G2]
						*G2 = 0.5*(l1*transp(*d1) + (csi_1[i][j] - csi_1_0[i][j])*(*G1));
						*G2b = 0.5*(l2*transp(*d2) + (csi_2[i][j] - csi_2_0[i][j])*(*G1b));
						//Matriz [G3]
						*G3 = transp(*dN_ext)*(4.0 / l1*(*I361 - (*t1ext)*transp(*t1ext)))*(*dN_ext);
						*G3b = transp(*dN_ext)*(4.0 / l2*(*I362 - (*t2ext)*transp(*t2ext)))*(*dN_ext);
						//Matriz [G4]
						for (int ii = 0; ii<3; ii++)
						{
							(*aux1)(ii + 0, 0) = -1.0*((*dN1)*(*x1))(ii, 0);
							(*aux1)(ii + 3, 0) = +1.0*((*dN1)*(*x1))(ii, 0);
							(*aux2)(ii + 0, 0) = ((*N2)*(*x2) - (*N1)*(*x1))(ii, 0);
						}
						*G4 = transp(*N_ext)*(*aux1) + transp(*dN_ext)*(*aux2);
						zeros(aux1);
						zeros(aux2);
						for (int ii = 0; ii<3; ii++)
						{
							(*aux1)(ii + 0, 0) = -1.0*((*dN2)*(*x2))(ii, 0);
							(*aux1)(ii + 3, 0) = +1.0*((*dN2)*(*x2))(ii, 0);
							(*aux2)(ii + 3, 0) = ((*N2)*(*x2) - (*N1)*(*x1))(ii, 0);
						}
						*G4b = transp(*N_ext)*(*aux1) + transp(*dN_ext)*(*aux2);
						//Matriz [G5]
						dot11 = dot((*dN1)*(*x1), (*dN1)*(*x1));
						dot22 = dot((*dN2)*(*x2), (*dN2)*(*x2));
						dot12 = dot((*dN2)*(*x2), (*dN1)*(*x1));
						for (int ii = 0; ii<3; ii++)
						{
							(*aux3)(ii + 0, 0) = (-dot22*(*dN1)*(*x1) + dot12*(*dN2)*(*x2))(ii, 0);
							(*aux3)(ii + 3, 0) = (-dot11*(*dN2)*(*x2) + dot12*(*dN1)*(*x1))(ii, 0);
						}
						*G5 = 2.0*transp(*aux3)*(*dN_ext);
						//Matriz [G6]
						for (int ii = 0; ii<3; ii++)
						{
							(*aux4)(ii + 0, 0) = ((*dN2)*(*x2))(ii, 0);
							(*aux5)(ii + 0, 0) = ((*dN1)*(*x1))(ii, 0);
							(*aux6)(ii + 3, 0) = ((*dN1)*(*x1))(ii, 0);
							(*aux7)(ii + 3, 0) = ((*dN2)*(*x2))(ii, 0);
						}
						for (int ii = 0; ii<3; ii++)
						{
							(*Q1)(ii + 0, ii + 0) = -1.0;
							(*Q1)(ii + 3, ii + 0) = +1.0;

							(*Q2)(ii + 0, ii + 0) = -1.0;
							(*Q2)(ii + 0, ii + 3) = +1.0;

							(*Q3)(ii + 0, ii + 3) = -1.0;
							(*Q3)(ii + 3, ii + 3) = +1.0;

							(*Q4)(ii + 3, ii + 0) = -1.0;
							(*Q4)(ii + 3, ii + 3) = +1.0;
						}
						*G6 = transp(*N_ext)*(*Q1)*(*dN_ext) + transp(*dN_ext)*(*Q2)*(*N_ext) +
							transp(*dN_ext)*((*aux4)*(*d2) - 2.0*(*aux5)*(*d1) + (*aux6)*(*d2));
						*G6b = transp(*N_ext)*(*Q3)*(*dN_ext) + transp(*dN_ext)*(*Q4)*(*N_ext) +
							transp(*dN_ext)*(2.0*(*aux7)*(*d2) - (*aux6)*(*d1) - (*aux4)*(*d1));
						//Matriz [G7]
						D_ = -1.0* (dot11 * dot22) + dot12*dot12;
						*G7 = ((dot22 / (D_*D_))*(*G4) - ((1.0*dot12) / (D_*D_))*(*G4b)) * (*G5) -
							(1.0 / D_)*(2.0*(*G4)*transp(*aux7)*(*dN_ext) - (*G4b)*transp(*aux6)*(*dN_ext) - (*G4b)*transp(*aux4)*(*dN_ext) + dot22*(*G6) - dot12*(*G6b));
						*G7b = ((dot12 / (D_*D_))*(*G4) - ((1.0*dot11) / (D_*D_))*(*G4b)) * (*G5) -
							(1.0 / D_)*((*G4)*transp(*aux6)*(*dN_ext) + (*G4)*transp(*aux4)*(*dN_ext) - 2.0*(*G4b)*transp(*aux5)*(*dN_ext) + dot12*(*G6) - dot11*(*G6b));
						if (sticking[i][j] == true)
						{
							(*ST) = +ept*((*G2)*transp(*G2) + gte1*0.5*(transp(*d1)*transp(*G1) + l1*(*G7) + (*G1)*(*d1) + (csi_1[i][j] - csi_1_0[i][j])*(*G3)));
							(*STb) = +ept*((*G2b)*transp(*G2b) + gte2*0.5*(transp(*d2)*transp(*G1b) + l2*(*G7b) + (*G1b)*(*d2) + (csi_2[i][j] - csi_2_0[i][j])*(*G3b)));
						}
						else
						{
							(*ST) = -1.0*mu*epn*(m1*(*G2)*transp(*n)*(*Ntio) + m1*g_n[i][j] * 0.5* (transp(*d1)*transp(*G1) + l1*(*G7) + (*G1)*(*d1) + (csi_1[i][j] - csi_1_0[i][j])*(*G3))/*+ g_n*G2*...*/);
							(*STb) = -1.0*mu*epn*(m2*(*G2b)*transp(*n)*(*Ntio) + m2*g_n[i][j] * 0.5* (transp(*d2)*transp(*G1b) + l2*(*G7b) + (*G1b)*(*d2) + (csi_2[i][j] - csi_2_0[i][j])*(*G3b))/*+ g_n*G2b*...*/);
						}
						//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
						//Metodo de penalidades - construção dos esforços de contato e da matriz de rigidez tangente
						(*c_loading[i][j]) = 1.0*epn*g_n[i][j] * (transp(*Ntio)*(*n)) + Fat1*(*G2) + Fat2*(*G2b);
						(*c_stiffness[i][j]) = epn*(transp(*Ntio)*(*n)*(transp(*n)*(*Ntio)) + g_n[i][j] * ((*E) + transp(*E) + (*F) + (*G))) + (*ST) + (*STb);
					}//g_n < 0
					else//Se não houver contato
					{
						//Deslizamentos acumulados nulos
						g_t1s_temp[i][j] = 0.0;
						g_t2s_temp[i][j] = 0.0;
					}
				}//return_value == 0
			}
		}
	}
}

//Montagens - Newmark
void LRLR::MountDyn()
{

}
//Preenche a contribuição do contato nas matrizes globais
void LRLR::MountGlobal()
{
	for (int ei = 0; ei < n_elements1; ei++)
	{
		for (int ej = 0; ej < n_elements2; ej++)
		{
			temp_element1 = db.line_regions[n_LR1 - 1]->elements[ei];
			temp_element2 = db.line_regions[n_LR2 - 1]->elements[ej];
			if (return_value[ei][ej] == 0 && g_n[ei][ej] < 0)
			{
				//PlotContactStatus(ei, ej);
				//Variaveis temporarias para salvar a indexação global dos graus de liberdade
				int GL_global_1 = 0;
				int GL_global_2 = 0;
				double anterior = 0;
				for (int i = 0; i < 12; i++)
				{
					//////////////MONTAGEM DO VETOR DE ESFORÇOS DESBALANCEADOS//////////////////
					//Toma do vetor de GL globais, a indexação de cada grau de liberdade global
					if (i<6)//elemento 1
					{
						if (i < 3)
							GL_global_1 = db.nodes[db.elements[temp_element1 - 1]->nodes[0]-1]->GLs[i];
						else
							GL_global_1 = db.nodes[db.elements[temp_element1 - 1]->nodes[2]-1]->GLs[i-3];
					}
					else	//elemento 2
					{
						if (i < 9)
							GL_global_1 = db.nodes[db.elements[temp_element2 - 1]->nodes[0]-1]->GLs[i - 6];
						else
							GL_global_1 = db.nodes[db.elements[temp_element2 - 1]->nodes[2]-1]->GLs[i - 9];
					}

					//Caso o grau de liberdade seja livre:
					if (GL_global_1 > 0)
					{
						anterior = db.global_P_A(GL_global_1 - 1, 0);
						db.global_P_A(GL_global_1 - 1, 0) = anterior + (*c_loading[ei][ej])(i, 0);
					}
					else
					{
						anterior = db.global_P_B(-GL_global_1 - 1, 0);
						db.global_P_B(-GL_global_1 - 1, 0) = anterior + (*c_loading[ei][ej])(i, 0);
					}
					for (int j = 0; j < 12; j++)
					{
						//////////////////////MONTAGEM DA MATRIZ DE RIGIDEZ/////////////////////////
						//Toma do vetor de GL globais, a indexação de cada grau de liberdade global
						if (j<6)//elemento 1
						{
							if (j < 3)
								GL_global_2 = db.nodes[db.elements[temp_element1 - 1]->nodes[0]-1]->GLs[j];
							else
								GL_global_2 = db.nodes[db.elements[temp_element1 - 1]->nodes[2]-1]->GLs[j - 3];
						}
						else	//elemento 2
						{
							if (j < 9)
								GL_global_2 = db.nodes[db.elements[temp_element2 - 1]->nodes[0]-1]->GLs[j - 6];
							else
								GL_global_2 = db.nodes[db.elements[temp_element2 - 1]->nodes[2]-1]->GLs[j - 9];
						}

						//Caso os graus de liberdade sejam ambos livres (Matriz Kaa)
						if (GL_global_1 > 0 && GL_global_2 > 0)
						{
							db.global_stiffness_AA.setValue(GL_global_1 - 1, GL_global_2 - 1, (*c_stiffness[ei][ej])(i, j));
						}
						//Caso os graus de liberdade sejam ambos fixos (Matriz Kbb)
						if (GL_global_1 < 0 && GL_global_2 < 0)
						{
							db.global_stiffness_BB.setValue(-GL_global_1 - 1, -GL_global_2 - 1, (*c_stiffness[ei][ej])(i, j));
						}
						//Caso os graus de liberdade sejam livre e fixo (Matriz Kab)
						if (GL_global_1 > 0 && GL_global_2 < 0)
						{
							db.global_stiffness_AB.setValue(GL_global_1 - 1, -GL_global_2 - 1, (*c_stiffness[ei][ej])(i, j));
						}
						//Caso os graus de liberdade sejam fixo e livre (Matriz Kba)
						if (GL_global_1 < 0 && GL_global_2 > 0)
						{
							db.global_stiffness_BA.setValue(-GL_global_1 - 1, GL_global_2 - 1, (*c_stiffness[ei][ej])(i, j));
						}
					}//segundo loop
				}//primeiro loop
			}
		}
	}
}

//Calcula a banda gerada na matriz global pelo contato
void LRLR::Band(int* band_fixed, int* band_free)
{
	//Verifica os valores maximos dos GLs dos elementos 1 e 2
	temp_band_free = 0;
	temp_band_fixed = 0;

	lowest_free_global_DOF = db.number_nodes*db.number_GLs_node;
	highest_free_global_DOF = 0;
	lowest_fixed_global_DOF = db.number_nodes*db.number_GLs_node;
	highest_fixed_global_DOF = 0;
	free_marked = false;
	fixed_marked = false;
	for (int ei = 0; ei < n_elements1; ei++)
	{
		for (int ej = 0; ej < n_elements2; ej++)
		{
			if (return_value[ei][ej] == 0 && g_n[ei][ej] < 0)
			{
				temp_element1 = db.elements[db.line_regions[n_LR1 - 1]->elements[ei] - 1]->number;
				temp_element2 = db.elements[db.line_regions[n_LR2 - 1]->elements[ej] - 1]->number;

				for (int j = 0; j < 3; j++)//Percorre os 3 primeiros graus de liberdade de cada nó (esse contato só afeta esses graus de liberdade, não envolve as rotações)
				{
					//Checando cada um dos nós do elemento 1 
					for (int node = 0; node < 3; node++)
					{
						temp_node = db.nodes[db.elements[temp_element1 - 1]->nodes[node]-1]->number;
						//se o GL for livre e ativo
						if (db.nodes[temp_node - 1]->constraints[j] == 0 && db.nodes[temp_node - 1]->active_GL[j] == 1)
						{
							if (db.nodes[temp_node - 1]->GLs[j] < lowest_free_global_DOF)
								lowest_free_global_DOF = db.nodes[temp_node - 1]->GLs[j];
							if (db.nodes[temp_node - 1]->GLs[j] > highest_free_global_DOF)
								highest_free_global_DOF = db.nodes[temp_node - 1]->GLs[j];
							free_marked = true;
						}

						//se o GL for fixo e ativo
						if (db.nodes[temp_node - 1]->constraints[j] == 1 && db.nodes[temp_node - 1]->active_GL[j] == 1)
						{
							if (-db.nodes[temp_node - 1]->GLs[j] < lowest_fixed_global_DOF)
								lowest_fixed_global_DOF = -db.nodes[temp_node - 1]->GLs[j];
							if (-db.nodes[temp_node - 1]->GLs[j] > highest_fixed_global_DOF)
								highest_fixed_global_DOF = -db.nodes[temp_node - 1]->GLs[j];
							fixed_marked = true;
						}
					}

					//Checando cada um dos nós do elemento 2
					for (int node = 0; node < 3; node++)
					{
						temp_node = db.nodes[db.elements[temp_element2 - 1]->nodes[node]-1]->number;
						//se o GL for livre e ativo
						if (db.nodes[temp_node - 1]->constraints[j] == 0 && db.nodes[temp_node - 1]->active_GL[j] == 1)
						{
							if (db.nodes[temp_node - 1]->GLs[j] < lowest_free_global_DOF)
								lowest_free_global_DOF = db.nodes[temp_node - 1]->GLs[j];
							if (db.nodes[temp_node - 1]->GLs[j] > highest_free_global_DOF)
								highest_free_global_DOF = db.nodes[temp_node - 1]->GLs[j];
							free_marked = true;
						}

						//se o GL for fixo e ativo
						if (db.nodes[temp_node - 1]->constraints[j] == 1 && db.nodes[temp_node - 1]->active_GL[j] == 1)
						{
							if (-db.nodes[temp_node - 1]->GLs[j] < lowest_fixed_global_DOF)
								lowest_fixed_global_DOF = -db.nodes[temp_node - 1]->GLs[j];
							if (-db.nodes[temp_node - 1]->GLs[j] > highest_fixed_global_DOF)
								highest_fixed_global_DOF = -db.nodes[temp_node - 1]->GLs[j];
							fixed_marked = true;
						}
					}
				}
			}
		}
	}
	if (free_marked == true)
		temp_band_free = highest_free_global_DOF - lowest_free_global_DOF;
	else
		temp_band_free = 0;
	if (fixed_marked == true)
		temp_band_fixed = highest_fixed_global_DOF - lowest_fixed_global_DOF;
	else
		temp_band_fixed = 0;
	*band_fixed = temp_band_fixed;
	*band_free = temp_band_free;
}

//Pre-calculo de variaveis que e feito uma unica vez no inicio
void LRLR::PreCalc()
{
	(*I3)(0, 0) = 1.0;
	(*I3)(1, 1) = 1.0;
	(*I3)(2, 2) = 1.0;

	n_elements1 = db.line_regions[n_LR1 - 1]->n_elements;		//Numero de elementos na line region 1
	n_elements2 = db.line_regions[n_LR2 - 1]->n_elements;		//Numero de elementos na line region 2
	Alloc(n_elements1, n_elements2);							//Alocação dinamica, de acordo com o numero de elementos do LR1 e LR2
	//Verification of element types assigned to LR1 - only Pipe_1 elements are considered acceptable here
	for (int i = 0; i < n_elements1; i++)
	{
		temp_element1 = db.line_regions[n_LR1 - 1]->elements[i];
		if (typeid(*db.elements[temp_element1 - 1]) == typeid(Pipe_1))
		{
			typeOK1[i] = 1;
			Pipe_1* ptr_element = static_cast<Pipe_1*>(db.elements[temp_element1 - 1]);	//ptr_element is a pointer to element Pipe_1, index i from line_region
			//Calculos iniciais
			L1[i] = ptr_element->CalculateLength();
			r1[i] = ptr_element->De() / 2.0;
		}
		else
		{
			//Another possibility - Beam_1 with tubular cross section
			if (typeid(*db.elements[temp_element1 - 1]) == typeid(Beam_1) && typeid(*db.sections[db.elements[temp_element1 - 1]->section - 1]) == typeid(SecTube))
			{
				typeOK1[i] = 1;
				Beam_1* ptr_element = static_cast<Beam_1*>(db.elements[temp_element1 - 1]);	//ptr_element is a pointer to element Beam_1, index i from line_region
				SecTube* ptr_sec = static_cast<SecTube*>(db.sections[db.elements[temp_element1 - 1]->section - 1]);
				L1[i] = ptr_element->CalculateLength();
				r1[i] = ptr_sec->De / 2.0;
			}
			else
			{
				typeOK1[i] = 0;
				db.myprintf("Warning! Element defined in contact region %d is not valid!\n", number);
			}
		}
	}
	//Verification of element types assigned to LR2 - only Pipe_1 elements are considered acceptable here
	for (int i = 0; i < n_elements2; i++)
	{
		temp_element2 = db.line_regions[n_LR2 - 1]->elements[i];
		if (typeid(*db.elements[temp_element2 - 1]) == typeid(Pipe_1))
		{
			typeOK2[i] = 1;
			Pipe_1* ptr_element = static_cast<Pipe_1*>(db.elements[temp_element2 - 1]);	//ptr_element is a pointer to element Pipe_1, index i from line_region
			//Calculos iniciais
			L2[i] = ptr_element->CalculateLength();
			r2[i] = ptr_element->De() / 2.0;
		}
		else
		{
			//Another possibility - Beam_1 with tubular cross section
			if (typeid(*db.elements[temp_element2 - 1]) == typeid(Beam_1) && typeid(*db.sections[db.elements[temp_element2 - 1]->section - 1]) == typeid(SecTube))
			{
				typeOK2[i] = 1;
				Beam_1* ptr_element = static_cast<Beam_1*>(db.elements[temp_element2 - 1]);	//ptr_element is a pointer to element Beam_1, index i from line_region
				SecTube* ptr_sec = static_cast<SecTube*>(db.sections[db.elements[temp_element2 - 1]->section - 1]);
				L2[i] = ptr_element->CalculateLength();
				r2[i] = ptr_sec->De / 2.0;
			}
			else
			{
				typeOK2[i] = 0;
				db.myprintf("Warning! Element defined in contact region %d is not valid!\n", number);
			}
		}
	}
}

//checagem inicial do contato  - inicio de cada incremento
void LRLR::BeginStepCheck()
{

}

//Checks proximity between each beam from LR to AS
void LRLR::PinballCheck()
{
	//Searching is done for each element from LR1 - proximity check is made with each element from LR2
	for (int i = 0; i < n_elements1; i++)
	{
		//If the element type is compatible with the contact formulation
		if (typeOK1[i] == 1)
		{
			temp_element1 = db.line_regions[n_LR1 - 1]->elements[i];	//Current element 1
			zeros(mean_position1);	//Assigns zero values in all indexes from "mean_position1" vector
			//For the current element, the nodes are adressed to guess its position in space (roughly speaking)
			for (int node = 0; node < db.elements[temp_element1 - 1]->n_nodes; node++)
			{
				temp_node = db.elements[temp_element1 - 1]->nodes[node];
				(*mean_position1)(0, 0) += db.nodes[temp_node - 1]->copy_coordinates[0] + db.nodes[temp_node - 1]->displacements[0];
				(*mean_position1)(1, 0) += db.nodes[temp_node - 1]->copy_coordinates[1] + db.nodes[temp_node - 1]->displacements[1];
				(*mean_position1)(2, 0) += db.nodes[temp_node - 1]->copy_coordinates[2] + db.nodes[temp_node - 1]->displacements[2];
			}
			*mean_position1 = (1.0 / db.elements[temp_element1 - 1]->n_nodes)*(*mean_position1);	//Mean position is calculated
			for (int j = 0; j < n_elements2; j++)
			{
				if (typeOK2[j] == 1)
				{
					temp_element2 = db.line_regions[n_LR2 - 1]->elements[j];	//Current element 2
					zeros(mean_position2);	//Assigns zero values in all indexes from "mean_position2" vector
					for (int node = 0; node < db.elements[temp_element2 - 1]->n_nodes; node++)
					{
						temp_node = db.elements[temp_element2 - 1]->nodes[node];
						(*mean_position2)(0, 0) += db.nodes[temp_node - 1]->copy_coordinates[0] + db.nodes[temp_node - 1]->displacements[0];
						(*mean_position2)(1, 0) += db.nodes[temp_node - 1]->copy_coordinates[1] + db.nodes[temp_node - 1]->displacements[1];
						(*mean_position2)(2, 0) += db.nodes[temp_node - 1]->copy_coordinates[2] + db.nodes[temp_node - 1]->displacements[2];
					}
					*mean_position2 = (1.0 / db.elements[temp_element2 - 1]->n_nodes)*(*mean_position2);	//Mean position is calculated

					//For the current pair [i][j], assigns 1 or 0 based on comparison between the distance of mean_position to analytical surface and the pinball tolerance
					if (norm(*mean_position2 - *mean_position1) <= pinball)
						activate[i][j] = 1;	//Near to contact
					else
						activate[i][j] = 0;	//Far to contact
				}
				else
					activate[i][j] = 0;	//No contact
			}
		}
	}
}

//Aloca na memória as variaveis que dependem do numero de elementos
void LRLR::Alloc(int e_elements1, int e_elements2)
{
	n_elements1 = e_elements1;
	n_elements2 = e_elements2;
	activate = new bool*[n_elements1];
	for (int i = 0; i < n_elements1; i++)
		activate[i] = new bool[n_elements2];
	typeOK1 = new bool[n_elements1];
	L1 = new double[n_elements1];
	r1 = new double[n_elements1];
	typeOK2 = new bool[n_elements2];
	L2 = new double[n_elements2];
	r2 = new double[n_elements2];
	first_mount = new bool*[n_elements1];
	for (int i = 0; i < n_elements1; i++)
		first_mount[i] = new bool[n_elements2];
	flag_cross = new bool*[n_elements1];
	for (int i = 0; i < n_elements1; i++)
		flag_cross[i] = new bool[n_elements2];
	sticking = new bool*[n_elements1];
	for (int i = 0; i < n_elements1; i++)
		sticking[i] = new bool[n_elements2];
	last_g_n = new double*[n_elements1];
	for (int i = 0; i < n_elements1; i++)
		last_g_n[i] = new double[n_elements2];
	csi_1_0 = new double*[n_elements1];
	for (int i = 0; i < n_elements1; i++)
		csi_1_0[i] = new double[n_elements2];
	csi_2_0 = new double*[n_elements1];
	for (int i = 0; i < n_elements1; i++)
		csi_2_0[i] = new double[n_elements2];
	copy_csi_1_converged = new double*[n_elements1];
	for (int i = 0; i < n_elements1; i++)
		copy_csi_1_converged[i] = new double[n_elements2];
	copy_csi_2_converged = new double*[n_elements1];
	for (int i = 0; i < n_elements1; i++)
		copy_csi_2_converged[i] = new double[n_elements2];
	gt1s = new double*[n_elements1];
	for (int i = 0; i < n_elements1; i++)
		gt1s[i] = new double[n_elements2];
	gt2s = new double*[n_elements1];
	for (int i = 0; i < n_elements1; i++)
		gt2s[i] = new double[n_elements2];
	g_t1s_temp = new double*[n_elements1];
	for (int i = 0; i < n_elements1; i++)
		g_t1s_temp[i] = new double[n_elements2];
	g_t2s_temp = new double*[n_elements1];
	for (int i = 0; i < n_elements1; i++)
		g_t2s_temp[i] = new double[n_elements2];
	return_value = new int*[n_elements1];
	for (int i = 0; i < n_elements1; i++)
		return_value[i] = new int[n_elements2];
	last_return_value = new int*[n_elements1];
	for (int i = 0; i < n_elements1; i++)
		last_return_value[i] = new int[n_elements2];
	g_n = new double*[n_elements1];
	for (int i = 0; i < n_elements1; i++)
		g_n[i] = new double[n_elements2];
	csi_1 = new double*[n_elements1];
	for (int i = 0; i < n_elements1; i++)
		csi_1[i] = new double[n_elements2];
	csi_2 = new double*[n_elements1];
	for (int i = 0; i < n_elements1; i++)
		csi_2[i] = new double[n_elements2];
	factor_1 = new double*[n_elements1];
	for (int i = 0; i < n_elements1; i++)
		factor_1[i] = new double[n_elements2];
	factor_2 = new double*[n_elements1];
	for (int i = 0; i < n_elements1; i++)
		factor_2[i] = new double[n_elements2];

	for (int i = 0; i < n_elements1; i++)
	{
		for (int j = 0; j < n_elements2; j++)
		{
			first_mount[i][j] = true;
			flag_cross[i][j] = false;
			sticking[i][j] = false;
			last_g_n[i][j] = 0.0;
			csi_1_0[i][j] = 0.0;
			csi_2_0[i][j] = 0.0;
			copy_csi_1_converged[i][j] = 0.0;
			copy_csi_2_converged[i][j] = 0.0;
			gt1s[i][j] = 0.0;
			gt2s[i][j] = 0.0;
			g_t1s_temp[i][j] = 0.0;
			g_t2s_temp[i][j] = 0.0;
			return_value[i][j] = 0;
			last_return_value[i][j] = 0;
			g_n[i][j] = 0.0;
			csi_1[i][j] = 0.0;
			csi_2[i][j] = 0.0;
			factor_1[i][j] = 0.0;
			factor_2[i][j] = 0.0;
		}
	}
	last_dif_pos = new Matrix**[n_elements1];
	for (int i = 0; i < n_elements1; i++)
	{
		last_dif_pos[i] = new Matrix*[n_elements2];
		for (int j = 0; j < n_elements2; j++)
			last_dif_pos[i][j] = new Matrix(3,1);
	}
	c_stiffness = new Matrix**[n_elements1];
	for (int i = 0; i < n_elements1; i++)
	{
		c_stiffness[i] = new Matrix*[n_elements2];
		for (int j = 0; j < n_elements2; j++)
			c_stiffness[i][j] = new Matrix(12, 12);
	}
	c_loading = new Matrix**[n_elements1];
	for (int i = 0; i < n_elements1; i++)
	{
		c_loading[i] = new Matrix*[n_elements2];
		for (int j = 0; j < n_elements2; j++)
			c_loading[i][j] = new Matrix(12, 1);
	}
	z2z1 = new Matrix**[n_elements1];
	for (int i = 0; i < n_elements1; i++)
	{
		z2z1[i] = new Matrix*[n_elements2];
		for (int j = 0; j < n_elements2; j++)
			z2z1[i][j] = new Matrix(3, 1);
	}
}

//Preenche as variaveis dos nós com valores atualizados
void LRLR::FillNodes(int e_element1, int e_element2)
{
	temp_node = db.elements[e_element1 - 1]->nodes[0];
	(*node_1_1)(0, 0) = db.nodes[temp_node - 1]->copy_coordinates[0] + db.nodes[temp_node - 1]->displacements[0];
	(*node_1_1)(1, 0) = db.nodes[temp_node - 1]->copy_coordinates[1] + db.nodes[temp_node - 1]->displacements[1];
	(*node_1_1)(2, 0) = db.nodes[temp_node - 1]->copy_coordinates[2] + db.nodes[temp_node - 1]->displacements[2];
	temp_node = db.elements[e_element1 - 1]->nodes[2];
	(*node_3_1)(0, 0) = db.nodes[temp_node - 1]->copy_coordinates[0] + db.nodes[temp_node - 1]->displacements[0];
	(*node_3_1)(1, 0) = db.nodes[temp_node - 1]->copy_coordinates[1] + db.nodes[temp_node - 1]->displacements[1];
	(*node_3_1)(2, 0) = db.nodes[temp_node - 1]->copy_coordinates[2] + db.nodes[temp_node - 1]->displacements[2];
	temp_node = db.elements[e_element2 - 1]->nodes[0];
	(*node_1_2)(0, 0) = db.nodes[temp_node - 1]->copy_coordinates[0] + db.nodes[temp_node - 1]->displacements[0];
	(*node_1_2)(1, 0) = db.nodes[temp_node - 1]->copy_coordinates[1] + db.nodes[temp_node - 1]->displacements[1];
	(*node_1_2)(2, 0) = db.nodes[temp_node - 1]->copy_coordinates[2] + db.nodes[temp_node - 1]->displacements[2];
	temp_node = db.elements[e_element2 - 1]->nodes[2];
	(*node_3_2)(0, 0) = db.nodes[temp_node - 1]->copy_coordinates[0] + db.nodes[temp_node - 1]->displacements[0];
	(*node_3_2)(1, 0) = db.nodes[temp_node - 1]->copy_coordinates[1] + db.nodes[temp_node - 1]->displacements[1];
	(*node_3_2)(2, 0) = db.nodes[temp_node - 1]->copy_coordinates[2] + db.nodes[temp_node - 1]->displacements[2];
	////BLOCO - Posições nodais 1/////
	(*x1)(0, 0) = (*node_1_1)(0, 0);
	(*x1)(1, 0) = (*node_1_1)(1, 0);
	(*x1)(2, 0) = (*node_1_1)(2, 0);
	(*x1)(3, 0) = (*node_3_1)(0, 0);
	(*x1)(4, 0) = (*node_3_1)(1, 0);
	(*x1)(5, 0) = (*node_3_1)(2, 0);
	////BLOCO - Posições nodais 2/////
	(*x2)(0, 0) = (*node_1_2)(0, 0);
	(*x2)(1, 0) = (*node_1_2)(1, 0);
	(*x2)(2, 0) = (*node_1_2)(2, 0);
	(*x2)(3, 0) = (*node_3_2)(0, 0);
	(*x2)(4, 0) = (*node_3_2)(1, 0);
	(*x2)(5, 0) = (*node_3_2)(2, 0);
	//////////////////////////////////
}

//Algoritmo para determinação de csi_1 e csi_2 - algoritmo de minimização de distancia
//retorno 0 - retorna valores de csi_1 e csi_2
int LRLR::FindMinimumParameters(int i, int j)
{
	//Guarda cópias de variaveis convergidas na ultima iteração - serão usadas caso o NR não convirja (por exemplo, se houver paralelismo)
	double copy_csi_1 = csi_1[i][j];
	double copy_csi_2 = csi_2[i][j];
	int flag_error = 0;
	//Inicialização de chute inicial - o ultimo valor convergido no NR
	csi_1[i][j] = copy_csi_1_converged[i][j];
	csi_2[i][j] = copy_csi_2_converged[i][j];
	error = tol_ortho + 1.0;	//Forçando entrar no loop
	int it = 1;
	if (error > tol_ortho)
	{
		//LOOP principal
		while (error > tol_ortho && it <= max_it)
		{
			EvaluateParameters(i, j);
			error = norm(*R);	//Norma infinito
			//Só faz mais uma iteração, se ainda não convergiu - caso contrario todas as variaveis ja ficam salvas com o ultimo valor convergido
			//Isso e importante pois todas essas variaveis são utilizadas para a montagem da matriz de rigidez do contato
			if (error > tol_ortho)
			{
				//Resolve sistema linear
				*delta_csi = fullsystem(*A, -1.0*(*R), &flag_error);
				if (flag_error == 0)
				{
					//Atualização de csi_1
					csi_1[i][j] += (*delta_csi)(0, 0);
					csi_2[i][j] += (*delta_csi)(1, 0);
					it++;
				}
				else//força saida do loop
				{
					it = max_it + 1;
				}
			}
		}//end while NR
	}
	//Se não houver convergência
	if (it > max_it)
	{
		csi_1[i][j] = copy_csi_1;
		csi_2[i][j] = copy_csi_2;
		EvaluateParameters(i, j);
		error = norm(*R);	//Norma infinito
		if (error > tol_ortho)
			return 2;
		else
			return 0;
	}
	//Ponto fora do dominio da barra
	if (abs(csi_1[i][j]) > 1.0 || abs(csi_2[i][j]) > 1.0)
	{
		//Se for um ponto das extremidades (tips) da line region
		if (i == (n_elements1 - 1) || i == 0 || j == (n_elements2 - 1) || j == 0)
		{
			csi_1[i][j] = 0.0;
			csi_2[i][j] = 0.0;
			EvaluateParameters(i, j);
			return 0;
		}
		return 3;
		//if (csi_1[i][j] < -1.0)
		//	csi_1[i][j] = -1.0;
		//else
		//{
		//	if (csi_1[i][j] > +1.0)
		//		csi_1[i][j] = +1.0;
		//}
		//if (csi_2[i][j] < -1.0)
		//	csi_2[i][j] = -1.0;
		//else
		//{
		//	if (csi_2[i][j] > +1.0)
		//		csi_2[i][j] = +1.0;
		//}
		//EvaluateParameters(i,j);
		//error = norm(*R);	//Norma infinito
		//flag_kind = 3;
		//if (error > tol_ortho)
		//	return 3;
		//else
		//	return 0;
	}
	//Caso chegue aqui significa que retornara o ponto do Newton-Raphson, convergido devidamente e que obedece as relações de ortogonalidade
	return 0;
}

//Calcula parametros em função de csi_1 e csi_2
void LRLR::EvaluateParameters(int i, int j)
{
	////BLOCO - ELEMENTO 1 /////
	N_1_1 = 0.5*(1.0 - csi_1[i][j]);
	N_3_1 = 0.5*(1.0 + csi_1[i][j]);
	dN_1_1 = -0.5;
	dN_3_1 = +0.5;
	ddN_1_1 = 0.0;
	ddN_3_1 = 0.0;
	(*N1)(0, 0) = N_1_1;
	(*N1)(1, 1) = N_1_1;
	(*N1)(2, 2) = N_1_1;
	(*N1)(0, 3) = N_3_1;
	(*N1)(1, 4) = N_3_1;
	(*N1)(2, 5) = N_3_1;
	(*dN1)(0, 0) = dN_1_1;
	(*dN1)(1, 1) = dN_1_1;
	(*dN1)(2, 2) = dN_1_1;
	(*dN1)(0, 3) = dN_3_1;
	(*dN1)(1, 4) = dN_3_1;
	(*dN1)(2, 5) = dN_3_1;
	(*ddN1)(0, 0) = ddN_1_1;
	(*ddN1)(1, 1) = ddN_1_1;
	(*ddN1)(2, 2) = ddN_1_1;
	(*ddN1)(0, 3) = ddN_3_1;
	(*ddN1)(1, 4) = ddN_3_1;
	(*ddN1)(2, 5) = ddN_3_1;
	////BLOCO - ELEMENTO 2 /////
	N_1_2 = 0.5*(1.0 - csi_2[i][j]);
	N_3_2 = 0.5*(1.0 + csi_2[i][j]);
	dN_1_2 = -0.5;
	dN_3_2 = +0.5;
	ddN_1_2 = 0.0;
	ddN_3_2 = 0.0;
	(*N2)(0, 0) = N_1_2;
	(*N2)(1, 1) = N_1_2;
	(*N2)(2, 2) = N_1_2;
	(*N2)(0, 3) = N_3_2;
	(*N2)(1, 4) = N_3_2;
	(*N2)(2, 5) = N_3_2;
	(*dN2)(0, 0) = dN_1_2;
	(*dN2)(1, 1) = dN_1_2;
	(*dN2)(2, 2) = dN_1_2;
	(*dN2)(0, 3) = dN_3_2;
	(*dN2)(1, 4) = dN_3_2;
	(*dN2)(2, 5) = dN_3_2;
	(*ddN2)(0, 0) = ddN_1_2;
	(*ddN2)(1, 1) = ddN_1_2;
	(*ddN2)(2, 2) = ddN_1_2;
	(*ddN2)(0, 3) = ddN_3_2;
	(*ddN2)(1, 4) = ddN_3_2;
	(*ddN2)(2, 5) = ddN_3_2;
	//Residuo
	(*R)(0, 0) = dot((*N2)*(*x2) - (*N1)*(*x1), (*dN1)*(*x1));
	(*R)(1, 0) = dot((*N2)*(*x2) - (*N1)*(*x1), (*dN2)*(*x2));
	//Monta Jacobiana (A)
	(*A)(0, 0) = -1.0*dot((*dN1)*(*x1), (*dN1)*(*x1)) + dot((*N2)*(*x2) - (*N1)*(*x1), (*ddN1)*(*x1));
	(*A)(0, 1) = +dot((*dN1)*(*x1), (*dN2)*(*x2));
	(*A)(1, 0) = -dot((*dN1)*(*x1), (*dN2)*(*x2));
	(*A)(1, 1) = dot((*dN2)*(*x2), (*dN2)*(*x2)) + dot((*N2)*(*x2) - (*N1)*(*x1), (*ddN2)*(*x2));
}

//Calcula o comprimento dos elementos
//Calcula os vetores t1 e t2 dos elementos
void LRLR::CalculateLengthsAndTangents()
{
	l1 = norm((*node_3_1) - (*node_1_1));
	l2 = norm((*node_3_2) - (*node_1_2));
	(*t1) = (1.0 / l1)*((*node_3_1) - (*node_1_1));
	(*t2) = (1.0 / l2)*((*node_3_2) - (*node_1_2));
}

//Retorna 0 - não ha cruzamento 1 - ha cruzamento
bool LRLR::HaveErrors()
{
	for (int i = 0; i < n_elements1; i++)
	{
		
		for (int j = 0; j < n_elements2; j++)
		{
			if (activate[i][j] == 1 && g_n[i][j] < 0.0 && flag_cross[i][j] == true && return_value[i][j] == 0)	//Verificação de "cross"
			{
				db.myprintf("\nLRLR Contact %d\n", number);
				temp_element1 = db.line_regions[n_LR1 - 1]->elements[i];
				temp_element2 = db.line_regions[n_LR2 - 1]->elements[j];
				db.myprintf("Crossing between element %d and %d. Cutback will be applied.\n", temp_element1, temp_element2);
				//last_dif_pos[i][j]->print();
				//z2z1[i][j]->print();
				return true;
			}
		}
	}
	return false;
}