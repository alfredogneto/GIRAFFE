#include "GeneralPLR.h"
#include <direct.h>

#include "Dynamic.h"
#include "Matrix.h"
#include "LineRegion.h"
#include "Particle.h"
#include "Sphere.h"
#include "Element.h"
#include "Pipe_1.h"
#include "Beam_1.h"
#include "SecTube.h"
#include "PipeSection.h"
#include "Node.h"


//#include <omp.h>
#define PI 3.1415926535897932384626433832795
#include"Database.h"
//Variáveis globais
extern
Database db;

GeneralPLR::GeneralPLR()
{
	type_name = new char[20];//Nome do tipo do contato
	sprintf(type_name, "GeneralPLR");
	number = 0;
	n_particles = 0;
	n_LR = 0;
	n_elements = 0;
	c = 0;
	alloc_control = NULL;
	activate = NULL;
	c_loading = NULL;
	c_stiffness = NULL;
	c_damping = NULL;

	z1 = Matrix(3, 1);
	z2 = Matrix(3, 1);

	I3 = Matrix(3, 3);
	I3(0, 0) = 1.0;
	I3(1, 1) = 1.0;
	I3(2, 2) = 1.0;
}

GeneralPLR::~GeneralPLR()
{
	delete[] type_name;

	if (activate != NULL)
	{
		for (int i = 0; i < n_particles; i++)
			delete[] activate[i];
		delete[] activate;
	}

	if (c_stiffness != NULL && c_loading != NULL && c_damping != NULL)
	{
		//Desalocando as matrizes e vetores
		for (int i = 0; i < n_particles; i++)
		{
			for (int j = 0; j < n_elements; j++)
			{
				FreeSpecific(i, j);
				delete[] c_stiffness[i][j];
				delete[] c_damping[i][j];
				delete[] c_loading[i][j];
			}
			delete[] c_stiffness[i];
			delete[] c_damping[i];
			delete[] c_loading[i];
		}
		delete[] c_stiffness;
		delete[] c_damping;
		delete[] c_loading;
	}

	if (alloc_control != NULL)
	{
		for (int i = 0; i < n_particles; i++)
			delete[] alloc_control[i];
		delete[] alloc_control;
	}
}
void GeneralPLR::WriteVTK_XMLRender(FILE *f)
{

}

void GeneralPLR::WriteVTK_XMLForces(FILE *f)
{
	//TODO
}

bool GeneralPLR::Read(FILE *f)
{
	char s[1000];
	fscanf(f, "%s", s);
	number = atoi(s);
	fscanf(f, "%s", s);
	if (!strcmp(s, "LR"))
	{
		fscanf(f, "%s", s);
		n_LR = atoi(s);
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
	if (!strcmp(s, "CN"))
	{
		fscanf(f, "%s", s);
		c = atof(s);
	}
	else
		return false;
	return true;
}
void GeneralPLR::Write(FILE *f)
{
	fprintf(f, "GeneralPLR\t%d\tLR\t%d\tMU\t%.6e\tEPN\t%.6e\tEPT\t%.6e\tPinball\t%.6e\tCN\t%.6e\n",
		number, n_LR, mu, epn, ept, pinball,c);
}
//Escreve arquivo de resultados
void GeneralPLR::WriteResults(FILE *f)
{
	//DOES NOTHING
}
//Escreve no monitor do contato
void GeneralPLR::WriteMonitor(FILE *f, bool first_record, double time)
{
	//DOES NOTHING
}
//Pré-cálculo de variáveis que é feito uma única vez no início
void GeneralPLR::PreCalc()
{
	n_particles = db.number_particles;					//Número de partículas
	n_elements = db.line_regions[n_LR - 1]->n_elements;	//Número de elementos do line region
	Alloc();												//Alocação dinâmica, de acordo com o número de contatos
}

//Checa inconsistências no elemento para evitar erros de execução
bool GeneralPLR::Check()
{
	return true;
}

//Monta contatos
void GeneralPLR::Mount()
{
	int n_p1 = 0;		//particle
	int n_element = 0;	//element
	
	double r1, r2;	//Raios das partículas
	Sphere* p1;		//Ponteiro para partículas
	int temp_node;
	Matrix z1z2(3);	//distância entre centros
	double gn;		//gap normal
	Matrix n(3);	//direção normal de contato
	Matrix non(3, 3);
	double f;
	for (int i = 0; i < n_particles; i++)
	{
		n_p1 = db.particles[i]->number;
		for (int j = 0; j < db.line_regions[n_LR-1]->n_elements; j++)
		{
			n_element = db.line_regions[n_LR - 1]->elements[j];

			if (activate[i][j] == true)
			{
				/////////////////////////////////////////PARTÍCULA ESFÉRICA em contato com Viga////////////////////////////////////////////////////////
				if (typeid(*db.particles[i]) == typeid(Sphere))
				{
					bool evaluate = false;
					if (typeid(*db.elements[n_element - 1]) == typeid(Pipe_1))
					{
						r2 = db.pipe_sections[db.elements[n_element - 1]->section - 1]->De / 2;//Raio da viga
						evaluate = true;
					}
					else
					{
						if (typeid(*db.elements[n_element - 1]) == typeid(Beam_1) && typeid(*db.sections[db.elements[n_element - 1]->section - 1]) == typeid(SecTube))
						{
							SecTube* sec = static_cast<SecTube*>(db.sections[db.elements[n_element - 1]->section - 1]);
							r2 = sec->De / 2;
							evaluate = true;
						}
					}

					//Se os tipos forem adequados, computa o contato
					if (evaluate == true)
					{
						p1 = static_cast<Sphere*>(db.particles[n_p1 - 1]);
						r1 = p1->radius;

						//Posição da partícula
						temp_node = db.particles[n_p1 - 1]->node;
						z1(0, 0) = db.nodes[temp_node - 1]->copy_coordinates[0] + db.nodes[temp_node - 1]->displacements[0];
						z1(1, 0) = db.nodes[temp_node - 1]->copy_coordinates[1] + db.nodes[temp_node - 1]->displacements[1];
						z1(2, 0) = db.nodes[temp_node - 1]->copy_coordinates[2] + db.nodes[temp_node - 1]->displacements[2];
						bool alloced = false;
						for (int k = 0; k < 3; k++)
						{
							//Posição da esfera no elemento de viga
							temp_node = db.elements[n_element - 1]->nodes[k];
							z2(0, 0) = db.nodes[temp_node - 1]->copy_coordinates[0] + db.nodes[temp_node - 1]->displacements[0];
							z2(1, 0) = db.nodes[temp_node - 1]->copy_coordinates[1] + db.nodes[temp_node - 1]->displacements[1];
							z2(2, 0) = db.nodes[temp_node - 1]->copy_coordinates[2] + db.nodes[temp_node - 1]->displacements[2];
							//Cálculo da função GAP
							z1z2 = z1 - z2;
							gn = sqrt(dot(z1z2, z1z2)) - r1 - r2;
							f = (gn / (gn + r1 + r2));
							//Se houver, de fato, contato (penetração)
							if (gn < 0.0)
							{
								//printf("Contact between p1 = %d and e = %d. Gap = %.12e\n", n_p1, n_element, gn);
								
								if (alloced == false)
								{
									AllocSpecific(i, j);
									alloced = true;
								}
								//Direção normal
								n = (1.0 / norm(z1z2))*z1z2;
								non = dyadic(n, n);
								for (int l = 0; l < 3; l++)
								{
									(*c_loading[i][j][k])(l, 0) = epn*gn*n(l, 0);
									(*c_loading[i][j][k])(l + 3, 0) = -epn*gn*n(l, 0);
									
									for (int m = 0; m < 3; m++)
									{
										(*c_stiffness[i][j][k])(l, m) = epn*(f*I3(l, m) + (1 - f)*non(l, m));
										(*c_stiffness[i][j][k])(l + 3, m + 3) = epn*(f*I3(l, m) + (1 - f)*non(l, m));
										(*c_stiffness[i][j][k])(l + 3, m) = -epn*(f*I3(l, m) + (1 - f)*non(l, m));
										(*c_stiffness[i][j][k])(l, m + 3) = -epn*(f*I3(l, m) + (1 - f)*non(l, m));
										(*c_damping[i][j][k])(l, m) = c*non(l,m);
										(*c_damping[i][j][k])(l + 3, m + 3) = c*non(l, m);
										(*c_damping[i][j][k])(l + 3, m) = -c*non(l, m);
										(*c_damping[i][j][k])(l, m + 3) = -c*non(l, m);
									}
									
								}
								//(*c_loading[i][j][k]).print();
							}
							
						}
						//Caso não tenha nenhum contato (nas 3 esferas dos elemento)
						if (alloced == false)
							FreeSpecific(i, j);
					}
				}
			}//end of if activate[i][j] == true
			
		}

	}
}

//Montagens - Newmark
void GeneralPLR::MountDyn()
{
	//Varredura dos contatos
	int n_p1 = 0;		//particle
	int n_element = 0;	//element
	//Matrix disp(6);
	Matrix vel(6);
	//Matrix accel(6);
	Matrix v_ipp(6);//Estimativa da velocidade no instante posterior
	for (int i = 0; i < n_particles; i++)
	{
		n_p1 = db.particles[i]->number;
		for (int j = 0; j < db.line_regions[n_LR - 1]->n_elements; j++)
		{
			n_element = db.line_regions[n_LR - 1]->elements[j];

			//Se houver alocação do contato, computa alterações devido à dinâmica
			if (alloc_control[i][j] == true)
			{
				//Nó da partícula
				for (int ind = 0; ind < 3; ind++)
					vel(ind, 0) = db.nodes[db.particles[n_p1 - 1]->node - 1]->vel[ind];
					
				//Percorre os três nós da viga
				for (int k = 0; k < 3; k++)
				{
					//Modificações da dinâmica na matriz de rigidez:
					Dynamic* ptr_sol = static_cast<Dynamic*>(db.solution[db.current_solution_number - 1]);
					(*c_stiffness[i][j][k]) = (*c_stiffness[i][j][k]) + ptr_sol->a4*(*c_damping[i][j][k]);
					//Modificações da dinâmica nos esforços - presença das forças de amortecimento
					//Nó da viga
					for (int ind = 0; ind < 3; ind++)
						vel(ind+3, 0) = db.nodes[db.elements[n_element - 1]->nodes[k] - 1]->vel[ind];
					for (int index = 0; index < 3; index++)
					{
							v_ipp(index, 0) = vel(index, 0);
							v_ipp(index + 3, 0) = vel(index + 3, 0);
					}
					(*c_loading[i][j][k]) = (*c_loading[i][j][k]) + (*c_damping[i][j][k])*(v_ipp);
				}
			}
		}
	}
}

//Preenche a contribuição do contato nas matrizes globais
void GeneralPLR::MountGlobal()
{
	int n_p1 = 0;	//particle
	int n_element = 0;	//element
	int GL_global_1 = 0;
	int GL_global_2 = 0;
	double anterior = 0;

	for (int ni = 0; ni < n_particles; ni++)
	{
		n_p1 = db.particles[ni]->number;

		for (int nj = 0; nj < n_elements; nj++)
		{
			n_element = db.line_regions[n_LR-1]->elements[nj];

			//Espalhamento das contribuições do contato entre n_p1 e n_element
			if (alloc_control[ni][nj] == true)//Se houver alocação
			{
				//3 esferas de cada elemento
				for (int k = 0; k < 3; k++)
				{
					for (int i = 0; i < 6; i++)
					{
						//(*c_loading[ni][nj][k]).print();
						//Partícula
						if (i<3)
							GL_global_1 = db.nodes[db.particles[n_p1 - 1]->node - 1]->GLs[i];
						//Esfera do elemento
						else
							GL_global_1 = db.nodes[db.elements[n_element - 1]->nodes[k] - 1]->GLs[i - 3];

						//Caso o grau de liberdade seja livre:
						if (GL_global_1 > 0)
						{
							anterior = db.global_P_A(GL_global_1 - 1, 0);
							db.global_P_A(GL_global_1 - 1, 0) = anterior + (*c_loading[ni][nj][k])(i, 0);
						}
						else
						{
							anterior = db.global_P_B(-GL_global_1 - 1, 0);
							db.global_P_B(-GL_global_1 - 1, 0) = anterior + (*c_loading[ni][nj][k])(i, 0);
						}
						for (int j = 0; j < 6; j++)
						{
							//Partícula 1
							if (j<3)
								GL_global_2 = db.nodes[db.particles[n_p1 - 1]->node - 1]->GLs[j];
							//Partícula 2
							else
								GL_global_2 = db.nodes[db.elements[n_element - 1]->nodes[k] - 1]->GLs[j - 3];

							//Caso os graus de liberdade sejam ambos livres (Matriz Kaa)
							if (GL_global_1 > 0 && GL_global_2 > 0)
							{
								db.global_stiffness_AA.setValue(GL_global_1 - 1, GL_global_2 - 1, (*c_stiffness[ni][nj][k])(i, j));
							}
							//Caso os graus de liberdade sejam ambos fixos (Matriz Kbb)
							if (GL_global_1 < 0 && GL_global_2 < 0)
							{
								db.global_stiffness_BB.setValue(-GL_global_1 - 1, -GL_global_2 - 1, (*c_stiffness[ni][nj][k])(i, j));
							}
							//Caso os graus de liberdade sejam livre e fixo (Matriz Kab)
							if (GL_global_1 > 0 && GL_global_2 < 0)
							{
								db.global_stiffness_AB.setValue(GL_global_1 - 1, -GL_global_2 - 1, (*c_stiffness[ni][nj][k])(i, j));
							}
							//Caso os graus de liberdade sejam fixo e livre (Matriz Kba)
							if (GL_global_1 < 0 && GL_global_2 > 0)
							{
								db.global_stiffness_BA.setValue(-GL_global_1 - 1, GL_global_2 - 1, (*c_stiffness[ni][nj][k])(i, j));
							}
						}
					}
				}
			}
		}
	}
}

//Salva variáveis para descrição lagrangiana atualizada
void GeneralPLR::SaveLagrange()
{
	//DOES NOTHING
}
//Retorna 1 - há algum erro, mesmo que tenha convergido. Ex: penetração excessiva
bool GeneralPLR::HaveErrors()
{
	//DOES NOTHING
	return false;
}

//Calcula a banda gerada na matriz global pelo contato
void GeneralPLR::Band(int* band_fixed, int* band_free)
{
	//DOES NOTHING
}

//Aloca estrutura de matrizes para endereçar vetores e matrizes devidos ao contato
void GeneralPLR::Alloc()
{
	activate = new bool*[n_particles];
	for (int i = 0; i < n_particles; i++)
		activate[i] = new bool[n_elements];
	
	alloc_control = new bool*[n_particles];
	for (int i = 0; i < n_particles; i++)
		alloc_control[i] = new bool[n_elements];
	
	//Alocação parcial da matriz de rigidez
	c_stiffness = new Matrix***[n_particles];
	for (int i = 0; i < n_particles; i++)
	{
		c_stiffness[i] = new Matrix**[n_elements];
		for (int j = 0; j < n_elements; j++)
		{
			c_stiffness[i][j] = new Matrix*[3];
		}
	}

	//Alocação parcial da matriz de amortecimento
	c_damping = new Matrix***[n_particles];
	for (int i = 0; i < n_particles; i++)
	{
		c_damping[i] = new Matrix**[n_elements];
		for (int j = 0; j < n_elements; j++)
		{
			c_damping[i][j] = new Matrix*[3];
		}
	}

	//Alocação parcial do vetor de carregamentos
	c_loading = new Matrix***[n_particles];
	for (int i = 0; i < n_particles; i++)
	{
		c_loading[i] = new Matrix**[n_elements];
		for (int j = 0; j < n_elements; j++)
		{
			c_loading[i][j] = new Matrix*[3];
		}
	}

	//Incialização de variáveis
	for (int i = 0; i < n_particles; i++)
	{
		for (int j = 0; j < n_elements; j++)
		{
			alloc_control[i][j] = false;
			activate[i][j] = false;
		}
	}
}

//Aloca matrizes para contato específico
void GeneralPLR::AllocSpecific(int i, int j)
{
	if (alloc_control[i][j] == false)
	{
		for (int k = 0; k < 3; k++)
		{
			c_stiffness[i][j][k] = new Matrix(6, 6);
			c_damping[i][j][k] = new Matrix(6, 6);
			c_loading[i][j][k] = new Matrix(6, 1);
		}
		alloc_control[i][j] = true;
	}
	else
	{
		for (int k = 0; k < 3; k++)
		{
			zeros(c_stiffness[i][j][k]);
			zeros(c_damping[i][j][k]);
			zeros(c_loading[i][j][k]);
		}
	}
}

//Desaloca matrizes para contato específico
void GeneralPLR::FreeSpecific(int i, int j)
{
	if (alloc_control[i][j] == true)
	{
		for (int k = 0; k < 3; k++)
		{
			delete c_stiffness[i][j][k];
			delete c_damping[i][j][k];
			delete c_loading[i][j][k];
		}
		alloc_control[i][j] = false;
	}
}

//checagem inicial do contato  - início de cada incremento
void GeneralPLR::BeginStepCheck()
{
	
}

//Checks proximity
void GeneralPLR::PinballCheck()
{
	int n_p1 = 0;		//particle
	int n_element = 0;	//element
	int temp_node;
	Matrix z1z2(3);	//distância entre centros
	
	//Searching is done for each pair of contact
	for (int i = 0; i < n_particles; i++)
	{
		n_p1 = db.particles[i]->number;
		for (int j = 0; j < db.line_regions[n_LR - 1]->n_elements; j++)
		{
			n_element = db.line_regions[n_LR - 1]->elements[j];
			//Posição da partícula
			temp_node = db.particles[n_p1 - 1]->node;
			z1(0, 0) = db.nodes[temp_node - 1]->copy_coordinates[0] + db.nodes[temp_node - 1]->displacements[0];
			z1(1, 0) = db.nodes[temp_node - 1]->copy_coordinates[1] + db.nodes[temp_node - 1]->displacements[1];
			z1(2, 0) = db.nodes[temp_node - 1]->copy_coordinates[2] + db.nodes[temp_node - 1]->displacements[2];
			//Posição da esfera central no elemento de viga
			temp_node = db.elements[n_element - 1]->nodes[1];
			z2(0, 0) = db.nodes[temp_node - 1]->copy_coordinates[0] + db.nodes[temp_node - 1]->displacements[0];
			z2(1, 0) = db.nodes[temp_node - 1]->copy_coordinates[1] + db.nodes[temp_node - 1]->displacements[1];
			z2(2, 0) = db.nodes[temp_node - 1]->copy_coordinates[2] + db.nodes[temp_node - 1]->displacements[2];

			if (norm(z1 - z2) <= pinball)
			{
				activate[i][j] = true;			//Near to contact
				//Ainda não aloca as matrizes - alocação será feita somente se o gap normal for negativo - isso é calculado na função Mount.
			}
			else
			{
				activate[i][j] = false;	//Far to contact
				FreeSpecific(i,j);
			}

		}

	}
}
