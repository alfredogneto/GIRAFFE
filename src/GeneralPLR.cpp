#include "GeneralPLR.h"
#include "PlatformCompat.h"

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
//Global variables
extern
Database db;

GeneralPLR::GeneralPLR()
{
	type_name = new char[20];//Name of the type of the contact
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
		//Desalocando the matrices and vectors
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
	//ALL
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
//Writes file of results
void GeneralPLR::WriteResults(FILE *f)
{
	//DOES NOTHING
}
//Writes in the monitor of the contact
void GeneralPLR::WriteMonitor(FILE *f, bool first_record, double time)
{
	//DOES NOTHING
}
//Pre-calculation of variables that and done a single time in the start
void GeneralPLR::PreCalc()
{
	n_particles = db.number_particles;					//Number of particles
	n_elements = db.line_regions[n_LR - 1]->n_elements;	//Number of elements of the line region
	Alloc();												//Allocation dynamics, according to the number of contacts
}

//Checks inconsistencies in the element for avoid errors of execution
bool GeneralPLR::Check()
{
	return true;
}

//Assembles contacts
void GeneralPLR::Mount()
{
	int n_p1 = 0;		//particle
	int n_element = 0;	//element
	
	double r1, r2;	//Raios of the particles
	Sphere* p1;		//Pointer for particles
	int temp_node;
	Matrix z1z2(3);	//distance between centros
	double gn;		//gap normal
	Matrix n(3);	//normal direction of contact
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
				/////////////////////////////////////////Particle ESFeRICA in contact with Viga////////////////////////////////////////////////////////
				if (typeid(*db.particles[i]) == typeid(Sphere))
				{
					bool evaluate = false;
					if (typeid(*db.elements[n_element - 1]) == typeid(Pipe_1))
					{
						r2 = db.pipe_sections[db.elements[n_element - 1]->section - 1]->De / 2;//Radius of the viga
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

					//If the types forem adequados, computes the contact
					if (evaluate == true)
					{
						p1 = static_cast<Sphere*>(db.particles[n_p1 - 1]);
						r1 = p1->radius;

						//Position of the particle
						temp_node = db.particles[n_p1 - 1]->node;
						z1(0, 0) = db.nodes[temp_node - 1]->copy_coordinates[0] + db.nodes[temp_node - 1]->displacements[0];
						z1(1, 0) = db.nodes[temp_node - 1]->copy_coordinates[1] + db.nodes[temp_node - 1]->displacements[1];
						z1(2, 0) = db.nodes[temp_node - 1]->copy_coordinates[2] + db.nodes[temp_node - 1]->displacements[2];
						bool alloced = false;
						for (int k = 0; k < 3; k++)
						{
							//Position of the sphere in the element of viga
							temp_node = db.elements[n_element - 1]->nodes[k];
							z2(0, 0) = db.nodes[temp_node - 1]->copy_coordinates[0] + db.nodes[temp_node - 1]->displacements[0];
							z2(1, 0) = db.nodes[temp_node - 1]->copy_coordinates[1] + db.nodes[temp_node - 1]->displacements[1];
							z2(2, 0) = db.nodes[temp_node - 1]->copy_coordinates[2] + db.nodes[temp_node - 1]->displacements[2];
							//Calculation of the function GAP
							z1z2 = z1 - z2;
							gn = sqrt(dot(z1z2, z1z2)) - r1 - r2;
							f = (gn / (gn + r1 + r2));
							//If there is, of fato, contact (penetração)
							if (gn < 0.0)
							{
								//printf("Contact between p1 = %d and e = %d. Gap = %.12e\n", n_p1, n_element, gn);
								
								if (alloced == false)
								{
									AllocSpecific(i, j);
									alloced = true;
								}
								//Normal direction
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
						//Case not tenha nenhum contact (in the 3 esferas of the element)
						if (alloced == false)
							FreeSpecific(i, j);
					}
				}
			}//end of if activate[i][j] == true
			
		}

	}
}

//Assemblies - Newmark
void GeneralPLR::MountDyn()
{
	//Sweep of the contacts
	int n_p1 = 0;		//particle
	int n_element = 0;	//element
	//Matrix disp(6);
	Matrix vel(6);
	//Matrix accel(6);
	Matrix v_ipp(6);//Estimativa of the velocity in the instant subsequent
	for (int i = 0; i < n_particles; i++)
	{
		n_p1 = db.particles[i]->number;
		for (int j = 0; j < db.line_regions[n_LR - 1]->n_elements; j++)
		{
			n_element = db.line_regions[n_LR - 1]->elements[j];

			//If there is allocation of the contact, computes alterações due the dynamics
			if (alloc_control[i][j] == true)
			{
				//Node of the particle
				for (int ind = 0; ind < 3; ind++)
					vel(ind, 0) = db.nodes[db.particles[n_p1 - 1]->node - 1]->vel[ind];
					
				//Iterates over the três nodes of the viga
				for (int k = 0; k < 3; k++)
				{
					//Modificações of the dynamics in the stiffness matrix:
					Dynamic* ptr_sol = static_cast<Dynamic*>(db.solution[db.current_solution_number - 1]);
					(*c_stiffness[i][j][k]) = (*c_stiffness[i][j][k]) + ptr_sol->a4*(*c_damping[i][j][k]);
					//Modificações of the dynamics in the internal forces - presença of the forces of damping
					//Node of the viga
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

//Fills the contribution of the contact in the matrices global
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

			//Espalhamento of the contributions of the contact between n_p1 and n_element
			if (alloc_control[ni][nj] == true)//If there is allocation
			{
				//3 esferas of each element
				for (int k = 0; k < 3; k++)
				{
					for (int i = 0; i < 6; i++)
					{
						//(*c_loading[ni][nj][k]).print();
						//Particle
						if (i<3)
							GL_global_1 = db.nodes[db.particles[n_p1 - 1]->node - 1]->GLs[i];
						//Sphere of the element
						else
							GL_global_1 = db.nodes[db.elements[n_element - 1]->nodes[k] - 1]->GLs[i - 3];

						//Case the degree of freedom be free:
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
							//Particle 1
							if (j<3)
								GL_global_2 = db.nodes[db.particles[n_p1 - 1]->node - 1]->GLs[j];
							//Particle 2
							else
								GL_global_2 = db.nodes[db.elements[n_element - 1]->nodes[k] - 1]->GLs[j - 3];

							//Case the degrees of freedom be both free (Matrix Kaa)
							if (GL_global_1 > 0 && GL_global_2 > 0)
							{
								db.global_stiffness_AA.setValue(GL_global_1 - 1, GL_global_2 - 1, (*c_stiffness[ni][nj][k])(i, j));
							}
							//Case the degrees of freedom be both fixed (Matrix Kbb)
							if (GL_global_1 < 0 && GL_global_2 < 0)
							{
								db.global_stiffness_BB.setValue(-GL_global_1 - 1, -GL_global_2 - 1, (*c_stiffness[ni][nj][k])(i, j));
							}
							//Case the degrees of freedom be free and fixed (Matrix Kab)
							if (GL_global_1 > 0 && GL_global_2 < 0)
							{
								db.global_stiffness_AB.setValue(GL_global_1 - 1, -GL_global_2 - 1, (*c_stiffness[ni][nj][k])(i, j));
							}
							//Case the degrees of freedom be fixed and free (Matrix Kba)
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

//Saves variables for description lagrangiana updated
void GeneralPLR::SaveLagrange()
{
	//DOES NOTHING
}
//Returns 1 - there is some error, same that tenha converged. Ex: penetração excessiva
bool GeneralPLR::HaveErrors()
{
	//DOES NOTHING
	return false;
}

//Calculates the band generated in the matrix global through the contact
void GeneralPLR::Band(int* band_fixed, int* band_free)
{
	//DOES NOTHING
}

//Allocates structure of matrices for endereçar vectors and matrices devidos to the contact
void GeneralPLR::Alloc()
{
	activate = new bool*[n_particles];
	for (int i = 0; i < n_particles; i++)
		activate[i] = new bool[n_elements];
	
	alloc_control = new bool*[n_particles];
	for (int i = 0; i < n_particles; i++)
		alloc_control[i] = new bool[n_elements];
	
	//Allocation partial of the stiffness matrix
	c_stiffness = new Matrix***[n_particles];
	for (int i = 0; i < n_particles; i++)
	{
		c_stiffness[i] = new Matrix**[n_elements];
		for (int j = 0; j < n_elements; j++)
		{
			c_stiffness[i][j] = new Matrix*[3];
		}
	}

	//Allocation partial of the damping matrix
	c_damping = new Matrix***[n_particles];
	for (int i = 0; i < n_particles; i++)
	{
		c_damping[i] = new Matrix**[n_elements];
		for (int j = 0; j < n_elements; j++)
		{
			c_damping[i][j] = new Matrix*[3];
		}
	}

	//Allocation partial of the vector of loadings
	c_loading = new Matrix***[n_particles];
	for (int i = 0; i < n_particles; i++)
	{
		c_loading[i] = new Matrix**[n_elements];
		for (int j = 0; j < n_elements; j++)
		{
			c_loading[i][j] = new Matrix*[3];
		}
	}

	//Incialização of variables
	for (int i = 0; i < n_particles; i++)
	{
		for (int j = 0; j < n_elements; j++)
		{
			alloc_control[i][j] = false;
			activate[i][j] = false;
		}
	}
}

//Allocates matrices for contact specific
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

//Deallocates matrices for contact specific
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

//check initial of the contact  - start of each increment
void GeneralPLR::BeginStepCheck()
{
	
}

//Checks proximity
void GeneralPLR::PinballCheck()
{
	int n_p1 = 0;		//particle
	int n_element = 0;	//element
	int temp_node;
	Matrix z1z2(3);	//distance between centros
	
	//Searching is done for each pair of contact
	for (int i = 0; i < n_particles; i++)
	{
		n_p1 = db.particles[i]->number;
		for (int j = 0; j < db.line_regions[n_LR - 1]->n_elements; j++)
		{
			n_element = db.line_regions[n_LR - 1]->elements[j];
			//Position of the particle
			temp_node = db.particles[n_p1 - 1]->node;
			z1(0, 0) = db.nodes[temp_node - 1]->copy_coordinates[0] + db.nodes[temp_node - 1]->displacements[0];
			z1(1, 0) = db.nodes[temp_node - 1]->copy_coordinates[1] + db.nodes[temp_node - 1]->displacements[1];
			z1(2, 0) = db.nodes[temp_node - 1]->copy_coordinates[2] + db.nodes[temp_node - 1]->displacements[2];
			//Position of the sphere central in the element of viga
			temp_node = db.elements[n_element - 1]->nodes[1];
			z2(0, 0) = db.nodes[temp_node - 1]->copy_coordinates[0] + db.nodes[temp_node - 1]->displacements[0];
			z2(1, 0) = db.nodes[temp_node - 1]->copy_coordinates[1] + db.nodes[temp_node - 1]->displacements[1];
			z2(2, 0) = db.nodes[temp_node - 1]->copy_coordinates[2] + db.nodes[temp_node - 1]->displacements[2];

			if (norm(z1 - z2) <= pinball)
			{
				activate[i][j] = true;			//Near to contact
				//Still not allocates the matrices - allocation will be done only if the gap normal for negative - this and calculated in the function Mount.
			}
			else
			{
				activate[i][j] = false;	//Far to contact
				FreeSpecific(i,j);
			}

		}

	}
}
