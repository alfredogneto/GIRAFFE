#include "NURBSMultipatchSurface.h"
#include <Eigen/dense>

#include "MatrixFloat.h"
#include "NURBSSurface.h"
#include "Particle.h"
#include "NURBSParticle.h"

#include"Database.h"
#include <cmath>
#include "Encoding.h"
//Variaveis globais
extern
Database db;

NURBSMultipatchSurface::NURBSMultipatchSurface()
{
	G = new Matrix(3);
	J_O = new Matrix(3, 3);

	number = 0;
	n_patches = 0;

	patches.clear();
}


NURBSMultipatchSurface::~NURBSMultipatchSurface()
{
	delete G;
	delete J_O;

	for (int i = 0; i < (int)patches.size(); i++)
	{
		delete patches[i];
	}

	patches.clear();
}

//Leitura do arquivo de entrada
bool NURBSMultipatchSurface::Read(FILE *f)
{
	char s[1000];
	bool CADread;

	fscanf(f, "%s", s);
	number = atoi(s);
	fscanf(f, "%s", s);
	strcpy(file, s);
	//Leitura do arquivo CAD
	CADread = ReadCADFile();
	if (CADread == false)
	{
		printf("Error reading CADData number %d\n", number);
		return false;
	}

	return true;
}

//Escrita do arquivo de saida
void NURBSMultipatchSurface::Write(FILE *f)
{
	fprintf(f, "NURBSMultipatchSurface\t%d\t%s\n",
		number,
		file
	);
}

void NURBSMultipatchSurface::PreCalc()
{
	for (int i = 0; i < n_patches; i++)
	{
		patches[i]->PreCalc();
	}
	
	if (ReadNurbsData()) {

	}
	else {
		EvaluateVolume();
		EvaluateCentroid();
		EvaluateInertiaTensor();
	}

	EvaluateRadius();
}

//Leitura do arquivo de CAD
bool  NURBSMultipatchSurface::ReadCADFile()
{
	//Faz leitura do arquivo 'file'
	//Atribui valor para variaveis de NURBS e aloca matrizes de NURBS

	FILE *f1 = NULL;
	char name_file[500];
	strcpy(name_file, db.folder_name);
	strcat(name_file, "CAD/");
	strcat(name_file, file);
	f1 = fopen(name_file, "r");
	if (f1 == NULL)
		return false;//Erro de leitura do arquivo de CAD

	//Leitura do aquivo
	char s[200];

	NURBSSurface *ns;

	//Leitura ID
	fscanf(f1, "%s", s);
	if (!strcmp(s, "ID"))
	{
		while (!strcmp(s, "ID"))
		{
			n_patches = n_patches + 1;
			ns = new NURBSSurface();
			fscanf(f1, "%s", s);
			ns->id_nurbs = atoi(s);
			//Leitura Subdivisions
			fscanf(f1, "%s", s);
			if (!strcmp(s, "Subdivisions"))
			{	
				fscanf(f1, "%s", s);
				ns->subdivisions[0] = atoi(s);
				fscanf(f1, "%s", s);
				ns->subdivisions[1] = atoi(s);
				//Leitura Connectivity
				fscanf(f1, "%s", s);
				if (!strcmp(s, "Connectivity"))
				{
					fscanf(f1, "%s", s);
					ns->connectivity[0] = atoi(s);
					fscanf(f1, "%s", s);
					ns->connectivity[1] = atoi(s);
					fscanf(f1, "%s", s);
					ns->connectivity[2] = atoi(s);
					fscanf(f1, "%s", s);
					ns->connectivity[3] = atoi(s);
					// Leitura UDim
					fscanf(f1, "%s", s);
					if (!strcmp(s, "UDim"))
					{
						fscanf(f1, "%s", s);
						ns->U_dim = atoi(s);
						//Leitura UOrder
						fscanf(f1, "%s", s);
						if (!strcmp(s, "UOrder"))
						{
							fscanf(f1, "%s", s);
							ns->U_order = atoi(s);
							//Aloca U knot vector
							ns->U_knot_vector = new double[ns->U_dim + ns->U_order + 1];
							//Leitura U knot vector
							fscanf(f1, "%s", s);
							if (!strcmp(s, "UKnotVector"))
							{
								for (int i = 0; i < (ns->U_dim + ns->U_order + 1); i++)
								{
									fscanf(f1, "%s", s);
									ns->U_knot_vector[i] = atof(s);
								}
								//Leitura VDim
								fscanf(f1, "%s", s);
								if (!strcmp(s, "VDim"))
								{
									fscanf(f1, "%s", s);
									ns->V_dim = atoi(s);
									//Leitura VOrder
									fscanf(f1, "%s", s);
									if (!strcmp(s, "VOrder"))
									{
										fscanf(f1, "%s", s);
										ns->V_order = atoi(s);
										//Aloca V knot vector
										ns->V_knot_vector = new double[ns->V_dim + ns->V_order + 1];
										//Leitura V knot vector
										fscanf(f1, "%s", s);
										if (!strcmp(s, "VKnotVector"))
										{
											for (int i = 0; i < (ns->V_dim + ns->V_order + 1); i++)
											{
												fscanf(f1, "%s", s);
												ns->V_knot_vector[i] = atof(s);
											}
											/*
											Vector weights and control points are read assuming they are in a sequence:
											v1 -> u1,...nn, then v2 -> u1,...,un, then ..., vm -> u1,...,un
											with that, if i varies along u and j varies along v, the position related to a given i,j in the vector weights is given at position [i][j]
											The same applies for control points.
											*/
											//Aloca weights
											ns->weights = new double*[ns->U_dim];
											for (int i = 0; i < ns->U_dim; i++)
												ns->weights[i] = new double[ns->V_dim];
											//Aloca control points
											ns->control_points = new Matrix*[ns->U_dim];
											for (int i = 0; i < ns->U_dim; i++)
											{
												ns->control_points[i] = new Matrix[ns->V_dim];
												for (int j = 0; j < ns->V_dim; j++)
													ns->control_points[i][j] = Matrix(3);
											}
											//Leitura weights
											fscanf(f1, "%s", s);
											if (!strcmp(s, "Weights"))
											{
												for (int j = 0; j < ns->V_dim; j++)
												{
													for (int i = 0; i < ns->U_dim; i++)
													{
														fscanf(f1, "%s", s);
														ns->weights[i][j] = atof(s);
													}
												}
											}
											else
												return false;
											//Leitura control points
											fscanf(f1, "%s", s);
											if (!strcmp(s, "ControlPoints"))
											{
												for (int j = 0; j < ns->V_dim; j++)
												{
													for (int i = 0; i < ns->U_dim; i++)
													{
														fscanf(f1, "%s", s);
														ns->control_points[i][j](0, 0) = atof(s);
														fscanf(f1, "%s", s);
														ns->control_points[i][j](1, 0) = atof(s);
														fscanf(f1, "%s", s);
														ns->control_points[i][j](2, 0) = atof(s);
													}
												}
											}
											else
												return false;
										}
										else
											return false;
									}
									else
										return false;
								}
								else
									return false;
							}
							else
								return false;
						}
						else
							return false;
					}
					else
						return false;
				}
				else
					return false;
			}
			else
				return false;
			fscanf(f1, "%s", s);
			patches.push_back(ns);
		}
		fclose(f1);
		return true;
	}
	else
		return false;
}


void NURBSMultipatchSurface::WriteVTK_XMLRender(FILE *f, Matrix& pos, Matrix& rot, int number)
{
	for (int i = 0; i < n_patches; i++)
	{
		patches[i]->WriteVTK_XMLRender(f, pos, rot, this->number, (number - 1)*n_patches + i + 1);
	}
}

void NURBSMultipatchSurface::EvaluateVolume()
{
	//TODO
	volume = 0.0;
}
void NURBSMultipatchSurface::EvaluateCentroid()
{
	//TODO
	//Centroid in local coordinate system
	(*G)(0, 0) = 0.0;
	(*G)(1, 0) = 0.0;
	(*G)(2, 0) = 0.0;
}
void NURBSMultipatchSurface::EvaluateInertiaTensor()
{
	//TODO
	//Inertia tensor with respect to the origin in local coordinate system
					//      | Jxx -Jxy -Jxz|     //
					//      |-Jxy  Jyy -Jyz|     //
					//      |-Jxz -Jyz  Jzz|     //
	(*J_O)(0, 0) = 0.0;
	(*J_O)(1, 0) = 0.0;
	(*J_O)(2, 0) = 0.0;
	(*J_O)(0, 1) = 0.0;
	(*J_O)(1, 1) = 0.0;
	(*J_O)(2, 1) = 0.0;
	(*J_O)(0, 2) = 0.0;
	(*J_O)(1, 2) = 0.0;
	(*J_O)(2, 2) = 0.0;

}
void NURBSMultipatchSurface::EvaluateRadius()
{
	//Marina
	//Evaluates the particle radius (w/r to the origin)
	//Loop on control points
	radius = -1.0;
	double temp_len = 0.0;
	for (int i = 0; i < n_patches; i++)
	{
		for (int j = 0; j < patches[i]->U_dim; j++)
		{
			for (int k = 0; k < patches[i]->V_dim; k++)
			{
				temp_len = sqrt(pow((patches[i]->control_points[j][k](0, 0) - (*G)(0, 0)), 2) + pow((patches[i]->control_points[j][k](1, 0) - (*G)(1, 0)), 2) + +pow((patches[i]->control_points[j][k](2, 0) - (*G)(2, 0)), 2));
				if (temp_len > radius) {
					radius = temp_len;
				}
				else {

				}
			}
		}
	}
}

//Marina
bool NURBSMultipatchSurface::ReadNurbsData()
{
	FILE *f1 = NULL;
	char name_file_1[500];
	strcpy(name_file_1, db.folder_name);
	strcat(name_file_1, "CAD/");
	strcat(name_file_1, file);
	//Correction on file termination - eliminating the ".txt" and replacing for ".snrb"
	char *temp;
	temp = strrchr(name_file_1, '.');   //Get the pointer to the last occurrence to the character '.'
	*temp = '\0';  //Replace token with null char
	strcat(name_file_1, ".snrb"); //Simplified/summarized file

	//Reading the file
	f1 = fopen(name_file_1, "r");
	if (f1 == NULL)
	{
		char name_file_2[500];
		strcpy(name_file_2, db.folder_name);
		strcat(name_file_2, "CAD/");
		strcat(name_file_2, file);
		//Correction on file termination - eliminating the ".txt" and replacing for ".rnrb"
		char *temp;
		temp = strrchr(name_file_2, '.');   //Get the pointer to the last occurrence to the character '.'
		*temp = '\0';  //Replace token with null char
		strcat(name_file_2, ".rnrb"); //Rhino file

		f1 = fopen(name_file_2, "r");

		if (f1 == NULL)
		{
			//Not available - return false
			return false;
		}
		else {
			//TO DO
			//Rhino file
			char ss[1000];

			fscanf(f1, "%s", ss);
			if (!strcmp(ss, "Volume"))
			{
				fscanf(f1, "%s", ss);
				fscanf(f1, "%s", ss);
				volume = atof(ss);
				fscanf(f1, "%s", ss);
				fscanf(f1, "%s", ss);
				fscanf(f1, "%s", ss);
				fscanf(f1, "%s", ss);
				fscanf(f1, "%s", ss);
				fscanf(f1, "%s", ss);
				if (!strcmp(ss, "Centroid"))
				{
					fscanf(f1, "%s", ss);
					fscanf(f1, "%[^,]s", ss);
					(*G)(0, 0) = atof(ss);
					fscanf(f1, "%c", ss);
					fscanf(f1, "%[^,]s", ss);
					(*G)(1, 0) = atof(ss);
					fscanf(f1, "%c", ss);
					fscanf(f1, "%s", ss);
					(*G)(2, 0) = atof(ss);
					fscanf(f1, "%s", ss);
					fscanf(f1, "%s", ss);
					fgets(ss, 200, f1);
					fgets(ss, 200, f1);
					fgets(ss, 200, f1);
					fgets(ss, 200, f1);
					fgets(ss, 200, f1);
					fgets(ss, 200, f1);
					fgets(ss, 200, f1);
					fgets(ss, 200, f1);
					fgets(ss, 200, f1);
					fgets(ss, 200, f1);
					fgets(ss, 200, f1);
					fscanf(f1, "%s", ss);
					if (!strcmp(ss, "xy:"))
					{
						fscanf(f1, "%s", ss);
						(*J_O)(0, 1) = -atof(ss);
						(*J_O)(1, 0) = (*J_O)(0, 1);
						fscanf(f1, "%s", ss);
						fscanf(f1, "%s", ss);
						fscanf(f1, "%s", ss);
						if (!strcmp(ss, "yz:"))
						{
							fscanf(f1, "%s", ss);
							(*J_O)(2, 1) = -atof(ss);
							(*J_O)(1, 2) = (*J_O)(2, 1);
							fscanf(f1, "%s", ss);
							fscanf(f1, "%s", ss);
							fscanf(f1, "%s", ss);
							if (!strcmp(ss, "zx:"))
							{
								fscanf(f1, "%s", ss);
								(*J_O)(0, 2) = -atof(ss);
								(*J_O)(2, 0) = (*J_O)(0, 2);
								fscanf(f1, "%s", ss);
								fscanf(f1, "%s", ss);
								fscanf(f1, "%s", ss);
								fgets(ss, 200, f1);
								fgets(ss, 200, f1);
								fgets(ss, 200, f1);
								fgets(ss, 200, f1);
								fgets(ss, 200, f1);
								fgets(ss, 200, f1);
								fgets(ss, 200, f1);
								fgets(ss, 200, f1);
								fgets(ss, 200, f1);
								fscanf(f1, "%s", ss);
								if (!strcmp(ss, "Ix:"))
								{
									fscanf(f1, "%s", ss);
									(*J_O)(0, 0) = atof(ss);
									fscanf(f1, "%s", ss);
									fscanf(f1, "%s", ss);
									fscanf(f1, "%s", ss);
									if (!strcmp(ss, "Iy:"))
									{
										fscanf(f1, "%s", ss);
										(*J_O)(1, 1) = atof(ss);
										fscanf(f1, "%s", ss);
										fscanf(f1, "%s", ss);
										fscanf(f1, "%s", ss);
										if (!strcmp(ss, "Iz:"))
										{
											fscanf(f1, "%s", ss);
											(*J_O)(2, 2) = atof(ss);
										}
										else {
											db.myprintf("\nError reading mesh file for CAD ID %d\n\n", number);
											fclose(f1);
											return false;
										}
									}
									else {
										db.myprintf("\nError reading mesh file for CAD ID %d\n\n", number);
										fclose(f1);
										return false;
									}
								}
								else {
									db.myprintf("\nError reading mesh file for CAD ID %d\n\n", number);
									fclose(f1);
									return false;
								}
							}
							else {
								db.myprintf("\nError reading mesh file for CAD ID %d\n\n", number);
								fclose(f1);
								return false;
							}
						}
						else {
							db.myprintf("\nError reading mesh file for CAD ID %d\n\n", number);
							fclose(f1);
							return false;
						}
					}
					else {
						db.myprintf("\nError reading mesh file for CAD ID %d\n\n", number);
						fclose(f1);
						return false;
					}
				}
				else {
					db.myprintf("\nError reading mesh file for CAD ID %d\n\n", number);
					fclose(f1);
					return false;
				}
			}
			else {
				db.myprintf("\nError reading mesh file for CAD ID %d\n\n", number);
				fclose(f1);
				return false;
			}
		}
		fclose(f1);
		return true;
	}
	else
	{
		//Simplified/summarized file
		char ss[1000];

		fscanf(f1, "%s", ss);
		if (!strcmp(ss, "VOLUME"))
		{
			fscanf(f1, "%s", ss);
			volume = (float)atof(ss);
			fscanf(f1, "%s", ss);
			if (!strcmp(ss, "CENTROID"))
			{
				fscanf(f1, "%s", ss);
				(*G)(0, 0) = (float)atof(ss);
				fscanf(f1, "%s", ss);
				(*G)(1, 0) = (float)atof(ss);
				fscanf(f1, "%s", ss);
				(*G)(2, 0) = (float)atof(ss);
				fscanf(f1, "%s", ss);
				fscanf(f1, "%s", ss);
				fscanf(f1, "%s", ss);
				if (!strcmp(ss, "INERTIA"))
				{
					fscanf(f1, "%s", ss);
					if (!strcmp(ss, "Ix"))
					{
						fscanf(f1, "%s", ss);
						(*J_O)(0, 0) = (float)atof(ss);
						fscanf(f1, "%s", ss);
						if (!strcmp(ss, "Iy"))
						{
							fscanf(f1, "%s", ss);
							(*J_O)(1, 1) = (float)atof(ss);
							fscanf(f1, "%s", ss);
							if (!strcmp(ss, "Iz"))
							{
								fscanf(f1, "%s", ss);
								(*J_O)(2, 2) = (float)atof(ss);
								fscanf(f1, "%s", ss);
								if (!strcmp(ss, "Ixy"))
								{
									fscanf(f1, "%s", ss);
									(*J_O)(0, 1) = -(float)atof(ss);
									(*J_O)(1, 0) = (*J_O)(0, 1);
									fscanf(f1, "%s", ss);
									if (!strcmp(ss, "Iyz"))
									{
										fscanf(f1, "%s", ss);
										(*J_O)(2, 1) = -(float)atof(ss);
										(*J_O)(1, 2) = (*J_O)(2, 1);
										fscanf(f1, "%s", ss);
										if (!strcmp(ss, "Izx"))
										{
											fscanf(f1, "%s", ss);
											(*J_O)(0, 2) = -(float)atof(ss);
											(*J_O)(2, 0) = (*J_O)(0, 2);
										}
										else {
											db.myprintf("\nError reading mesh file for CAD ID %d\n\n", number);
											fclose(f1);
											return false;
										}
									}
									else {
										db.myprintf("\nError reading mesh file for CAD ID %d\n\n", number);
										fclose(f1);
										return false;
									}
								}
								else {
									db.myprintf("\nError reading mesh file for CAD ID %d\n\n", number);
									fclose(f1);
									return false;
								}
							}
							else {
								db.myprintf("\nError reading mesh file for CAD ID %d\n\n", number);
								fclose(f1);
								return false;
							}
						}
						else {
							db.myprintf("\nError reading mesh file for CAD ID %d\n\n", number);
							fclose(f1);
							return false;
						}
					}
					else {
						db.myprintf("\nError reading mesh file for CAD ID %d\n\n", number);
						fclose(f1);
						return false;
					}
				}
				else {
					db.myprintf("\nError reading mesh file for CAD ID %d\n\n", number);
					fclose(f1);
					return false;
				}
			}
			else {
				db.myprintf("\nError reading mesh file for CAD ID %d\n\n", number);
				fclose(f1);
				return false;
			}
		}
		else
		{
			db.myprintf("\nError reading mesh file for CAD ID %d\n\n", number);
			fclose(f1);
			return false;
		}
		fclose(f1);
		return true;
	}

	//Marina
	delete f1;
	delete temp;
}

