#include "SurfacePair.h"
#include "PlatformCompat.h"

#include "SSContactData.h"
#include "Surface.h"
#include "FlexibleArcExtrusion_1_RigidArcRevolution_1.h"
#include"Database.h"
//Global variables
extern
Database db;

SurfacePair::SurfacePair()
{
}

SurfacePair::~SurfacePair()
{
}

//Values Default of tolerancias and other variables
void SurfacePair::DefaultValues()
{
	//Tolerance for precision of convective coordinates
	tol_convective = 1e-12;
	//Tolerance for precision of a small number (machine precision related)
	tol_small_1 = 1e-14;
	//Tolerance for precision of eigenvalues extraction
	tol_eig = 1e-14;
	//Factor to be used to modify the TR Dogleg Path when ascending directions are found
	tol_ascent = 1e-4;
	//Maximum iterations for searching minimum points
	max_it_1 = 50;
	//Maximum iterations for searching saddle points
	max_it_2 = 50;
	//Flag to write convergence report for LCP's
	//write_report = false;
	//Percentage of extension on the range of validity of convective coordinates to consider a possible solution of LCP still as "strong"
	perc = 0.5;

	convective_range = Matrix(4);
	convective_max = Matrix(4);
	convective_min = Matrix(4);

	seq_number = 0;

	alloc_control = false;
}

//Allocates memory
void SurfacePair::Alloc(SSContactData* c_data)
{
	if (alloc_control == false)
	{
		n_pointwise = c_data->n_solutions;

		cNR1 = new Matrix*[n_pointwise];
		cNR2 = new Matrix*[n_pointwise];
		cdeg = new Matrix*[n_pointwise];
		f_TR_report = new FILE*[n_pointwise];

		for (int i = 0; i < n_pointwise; i++)
		{
			cNR1[i] = new Matrix(4);
			cNR2[i] = new Matrix(4);
			cdeg[i] = new Matrix(4);
		}

		if (write_report)
		{
			for (int i = 0; i < n_pointwise; i++)
			{
				InitializeTRReport(i);
			}

		}

		alloc_control = true;
	}
	
}

//Frees memory
void SurfacePair::Free()
{
	if (alloc_control == true)
	{
		//Limpeza of memory
		for (int i = 0; i < n_pointwise; i++)
		{
			delete cNR1[i];
			delete cNR2[i];
			delete cdeg[i];
		}
		delete[]cNR1;
		delete[]cNR2;
		delete[]cdeg;
		delete[]f_TR_report;

		alloc_control = false;
	}
}

void SurfacePair::OpenTRReport(int index)
{
	//Writing the name of the contact pair
	char pair_name[100];
	sprintf(pair_name, "TR_report_surf_%d_surf_%d_%d", surf1_ID, surf2_ID, index);
	strcpy(name, db.folder_name);	//folder of the job
	strcat(name, "TR/");			//directory TR
	GIRAFFE_MKDIR(name);					//creating directory TR
	strcat(name, pair_name);		//file name
	strcat(name, ".txt");			//creating file
	f_TR_report[index] = fopen(name, "a");	//opens file
}

void SurfacePair::InitializeTRReport(int index)
{
	//Writing the name of the contact pair
	char pair_name[100];
	sprintf(pair_name, "TR_report_surf_%d_surf_%d_%d", surf1_ID, surf2_ID, index);
	strcpy(name, db.folder_name);	//folder of the job
	strcat(name, "TR/");			//directory TR
	GIRAFFE_MKDIR(name);					//creating directory TR
	strcat(name, pair_name);		//file name
	strcat(name, ".txt");			//creating file
	f_TR_report[index] = fopen(name, "w");	//opens file

	fprintf(f_TR_report[index], "///////////////////////////////////////////////////////////////////////\n");
	fprintf(f_TR_report[index], "InitializeConvectiveRange\n");
	fprintf(f_TR_report[index], "Convective 1:\t%.6e\tto\t%.6e.\tRange:\t%.6e\n", convective_min(0, 0), convective_max(0, 0), convective_range(0, 0));
	fprintf(f_TR_report[index], "Convective 2:\t%.6e\tto\t%.6e.\tRange:\t%.6e\n", convective_min(1, 0), convective_max(1, 0), convective_range(1, 0));
	fprintf(f_TR_report[index], "Convective 3:\t%.6e\tto\t%.6e.\tRange:\t%.6e\n", convective_min(2, 0), convective_max(2, 0), convective_range(2, 0));
	fprintf(f_TR_report[index], "Convective 4:\t%.6e\tto\t%.6e.\tRange:\t%.6e\n", convective_min(3, 0), convective_max(3, 0), convective_range(3, 0));
	fclose(f_TR_report[index]);

	fclose(f_TR_report[index]);			//closes file	
}

void SurfacePair::PreCalc()
{
	//Initializes and sets the range of validity for the convective coordinates
	InitializeConvectiveRange();
	
	//Setting the minimum convective range
	minimum_convective_range = 1e100;
	for (int i = 0; i < 4; i++)
	{
		if (convective_range(i, 0) < minimum_convective_range && convective_range(i, 0) != 0.0)
			minimum_convective_range = convective_range(i, 0);
	}
}

void SurfacePair::BeginStepCheck(SSContactData* c_data)
{
	if (write_report)
	{
		for (int i = 0; i < n_pointwise; i++)
		{
			OpenTRReport(i);
			fprintf(f_TR_report[i], "///////////////////////////////////////////////////////////////////////\n");
			fprintf(f_TR_report[i], "\nTime\t%.6f\tIteration\t%d\t", db.last_converged_time + db.current_time_step, db.current_iteration_number);
			fprintf(f_TR_report[i], "BeginStepCheck\n");
		}
	}
	
	///////////////////////////////Pointers and variables to make access easier//////////////////////////////////////////
	Surface* surf1 = db.surfaces[surf1_ID - 1];		//Pointer to the surface 1
	Surface* surf2 = db.surfaces[surf2_ID - 1];		//Pointer to the surface 2
	//In the start of each step and tomada the decisão acerca of strong_candidate or not (with base in the check of bounding box in torno of the surfaces)
	bool converged1 = false;
	bool converged2 = false;
	bool strong_candidate = false;
	//First check - bounding box overlap
	double inflation_factor = 2.0;
	strong_candidate = true;
	strong_candidate = BoxOverlap(surf1->box, surf2->box, inflation_factor);

	if (strong_candidate == false)
	{
		for (int ip = 0; ip < n_pointwise; ip++)
			c_data->return_value[ip] = 2;	//not and strong 
		if (write_report)
		{
			for (int i = 0; i < n_pointwise; i++)
			{
				fprintf(f_TR_report[i], "Candidate is not Strong. Out of bounding boxes overlap.\n");
				fclose(f_TR_report[i]);
			}
			
		}
		return;
	}

	//Case especial - Special LCP
	if (specialLCP)
	{
		//Saving ultimas convective coordinates convergidas in the current - evitando that cases that divergiram be usados as estimativas initial of the new increment
		for (int ip = 0; ip < n_pointwise; ip++)
		{
			c_data->degenerated[ip] = false;		//marca the principio as not degenerated for all the solutions
			c_data->repeated[ip] = false;			//Indicates that the solution [ip] always will be considerada

			//Degeneration basis - canonical basis
			(*c_data->P[ip])(0, 0) = 1.0;
			(*c_data->P[ip])(1, 1) = 1.0;
			(*c_data->P[ip])(2, 2) = 1.0;
			(*c_data->P[ip])(3, 3) = 1.0;
			//Degenerative operator
			c_data->MountDegenerativeOperator();
		}

		FlexibleArcExtrusion_1_RigidArcRevolution_1* temp = static_cast<FlexibleArcExtrusion_1_RigidArcRevolution_1*>(this);
		temp->SpecialLCP(*cNR2[0]);
		c_data->convective[0][0] = (*cNR2[0])(0, 0);
		c_data->convective[0][1] = (*cNR2[0])(1, 0);
		c_data->convective[0][2] = (*cNR2[0])(2, 0);
		c_data->convective[0][3] = (*cNR2[0])(3, 0);
		
		c_data->copy_return_value[0] = 2;
		c_data->return_value[0] = 0;

		for (int ip = 1; ip < n_pointwise; ip++)
			c_data->repeated[ip] = true;			//Indicates that the solution [ip] always will be desconsiderada

		if (write_report)
		{
			for (int ip = 0; ip < n_pointwise; ip++)
				fclose(f_TR_report[ip]);
		}
		return;
	}

	//Saving ultimas convective coordinates convergidas in the current - evitando that cases that divergiram be usados as estimativas initial of the new increment
	for (int ip = 0; ip < n_pointwise; ip++)
	{
		c_data->convective[ip][0] = c_data->copy_convective[ip][0];
		c_data->convective[ip][1] = c_data->copy_convective[ip][1];
		c_data->convective[ip][2] = c_data->copy_convective[ip][2];
		c_data->convective[ip][3] = c_data->copy_convective[ip][3];
		
		c_data->degenerated[ip] = false;		//marca the principio as not degenerated for all the solutions
		c_data->repeated[ip] = false;			//Indicates that the solution [ip] always will be considerada
	}
	
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////INPUT DATA (degeneration)//////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//Checking degeneration of surfaces involved in contact
	//Degeneration of the contact - flag
	if (surf1->degeneration[0] == true || surf1->degeneration[1] == true || surf2->degeneration[0] == true || surf2->degeneration[1] == true)
	{
		for (int ip = 0; ip < n_pointwise; ip++)
			c_data->degenerated[ip] = true;
	}
	else
	{
		for (int ip = 0; ip < n_pointwise; ip++)
			c_data->degenerated[ip] = false;
	}
	//Divisions in each surface
	int div1_1 = surf1->div1;
	int div2_1 = surf1->div2;
	int div1_2 = surf2->div1;
	int div2_2 = surf2->div2;
	int index_u1_1 = 0;
	int index_u2_1 = 0;
	int index_u1_2 = 0;
	int index_u2_2 = 0;

	for (int ip = 0; ip < n_pointwise; ip++)
	{ 
		converged1 = false;
		converged2 = false;

		seq_number = ip;//for control of the report
		//Update of the indices of the solutions of the surfaces
		if (index_u2_2 + 1 >= div2_2)
		{
			index_u2_2 = 0;
			if (index_u1_2 + 1 >= div1_2)
			{
				index_u1_2 = 0;
				if (index_u2_1 + 1 >= div2_1)
				{
					index_u2_1 = 0;
					if (index_u1_1 + 1 >= div1_1)
					{
						index_u1_1 = 0;
					}
					else
						index_u1_1++;
				}
				else
					index_u2_1++;
			}
			else
				index_u1_2++;
		}
		else
			index_u2_2++;

		//Degeneration basis - canonical basis
		(*c_data->P[ip])(0, 0) = 1.0;
		(*c_data->P[ip])(1, 1) = 1.0;
		(*c_data->P[ip])(2, 2) = 1.0;
		(*c_data->P[ip])(3, 3) = 1.0;
		//Degenerated coordinates index (in the new basis)
		c_data->deg_control[ip][0] = surf1->degeneration[0];
		c_data->deg_control[ip][1] = surf1->degeneration[1];
		c_data->deg_control[ip][2] = surf2->degeneration[0];
		c_data->deg_control[ip][3] = surf2->degeneration[1];
		//Degenerated coordinates values - in case of degeneration, just fill the desired coordinate value
		c_data->copy_deg_coordinates[ip][0] = surf1->deg_coordinates[index_u1_1][index_u2_1][0];
		c_data->copy_deg_coordinates[ip][1] = surf1->deg_coordinates[index_u1_1][index_u2_1][1];
		c_data->copy_deg_coordinates[ip][2] = surf2->deg_coordinates[index_u1_2][index_u2_2][0];
		c_data->copy_deg_coordinates[ip][3] = surf2->deg_coordinates[index_u1_2][index_u2_2][1];

		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		//Fixing convective coordinates
		Matrix temp(4);
		for (int i = 0; i < 4; i++)
			if (c_data->deg_control[ip][i] == true)
				temp(i, 0) = c_data->copy_deg_coordinates[ip][i];
		//Transforming into original basis
		temp = (*c_data->P[ip])*temp;
		for (int i = 0; i < 4; i++)
			c_data->convective[ip][i] = temp(i, 0);
		//Degenerative operator
		c_data->MountDegenerativeOperator();

		
		//CASE 1: not era strong candidate (or first calculation)
		if (c_data->copy_return_value[ip] == 2)
		{
			if (write_report)
				fprintf(f_TR_report[ip], "Performing initial guess\n");
			InitialGuess(c_data); //Performs initial guess with criterion geometrico -> writes in 'convective'
			/*c_data->convective[ip][0] = 1.0;
			c_data->convective[ip][1] = 1.0;
			c_data->convective[ip][2] = 1.0;
			c_data->convective[ip][3] = 1.0;*/
			for (int i = 0; i < 4; i++)
				(*cNR1[ip])(i, 0) = c_data->convective[ip][i];
		}
		//CASE 2: or and copy_return_value 0 or 4 (already there is solution disponivel proveniente of a step previous)
		else
		{
			(*cNR1[ip])(0, 0) = c_data->copy_convective[ip][0];
			(*cNR1[ip])(1, 0) = c_data->copy_convective[ip][1];
			(*cNR1[ip])(2, 0) = c_data->copy_convective[ip][2];
			(*cNR1[ip])(3, 0) = c_data->copy_convective[ip][3];
		}

		c_data->return_value[ip] = 0;		//e strong (default)
		int info = 0;
		
		int charact1 = 3;
		////////////////NO DEGENERATION///////////////
		if (c_data->degenerated[ip] == false)
		{
			//If the contact era previamente active or vizinho - or be, if there was solution of point of saddle
			if ((c_data->copy_return_value[ip] == 0 || c_data->copy_return_value[ip] == 4) && c_data->copy_g_n[ip] < 0.0)
			{
				if (write_report)
					fprintf(f_TR_report[ip], "Direct search for saddle\n");
				//Attempt of achar direto the saddle
				*cNR2[ip] = *cNR1[ip];
				converged2 = FindSaddleSolution(c_data, cNR2[ip], info, false);
			}
			//If not entrou in the if previous or if not converged the if previous
			if (converged2 == false)
			{
				if (write_report)
					fprintf(f_TR_report[ip], "Search for minimum\n");
				//Determination of minimum or intersection
				converged1 = FindMinimumSolution(c_data, cNR1[ip], info);
				*cNR2[ip] = *cNR1[ip];
				charact1 = CharacterizeCriticalPoint(cNR1[ip]);

				//If for intersection - search the point of saddle
				if (charact1 == 1)
					converged2 = FindSaddleSolution(c_data, cNR2[ip], info, true);
				//If not, returns the minimum
				else
					converged2 = converged1;
			}
		}
		////////////////DEGENERATION///////////////
		else
		{
			//Correction in case of degeneration
			Matrix temp_coordinates(4);
			for (int index = 0; index < 4; index++)
				temp_coordinates(index, 0) = (*cNR1[ip])(index, 0);
			//Transforming into degeneration basis
			temp_coordinates = transp(*c_data->P[ip])*temp_coordinates;
			//Fixing degenerated coordinates
			for (int i = 0; i < 4; i++)
				if (c_data->deg_control[ip][i] == true)
					temp_coordinates(i, 0) = c_data->copy_deg_coordinates[ip][i];
			//Transforming into original basis
			temp_coordinates = (*c_data->P[ip])*temp_coordinates;
			//Copying info into cNR1 vector
			for (int i = 0; i < 4; i++)
				(*cNR1[ip])(i, 0) = temp_coordinates(i, 0);

			//If the contact era previamente active (return 0 or 4)
			if ((c_data->copy_return_value[ip] == 0 || c_data->copy_return_value[ip] == 4) && c_data->copy_g_n[ip] < 0.0)
			{
				if (write_report)
					fprintf(f_TR_report[ip], "Direct search for saddle or minimum\n");
				*cNR2[ip] = *cNR1[ip];
				if (c_data->P_0[ip]->getColumns() == 3)
					converged2 = FindSaddleSolutionDegenerated(c_data, c_data->P_0[ip], cNR2[ip], false);
				else
					converged2 = FindMinimumSolutionDegenerated(c_data, c_data->P_0[ip], cNR2[ip]);
			}
			//If not entrou in the if previous or if not converged the saddle of the if previous
			if (converged2 == false)
			{
				
				//Determines minimum (or intersection)
				converged1 = FindMinimumSolutionDegenerated(c_data, c_data->P_0[ip], cNR1[ip]);
				*cNR2[ip] = *cNR1[ip];
				
			
				if (c_data->P_0[ip]->getColumns() == 3)
				{
					charact1 = CharacterizeCriticalPointDegenerated(cNR1[ip], c_data->P_0[ip]);
					//Search of the point of saddle
					if (charact1 == 1)
						converged2 = FindSaddleSolutionDegenerated(c_data, c_data->P_0[ip], cNR2[ip], true);
					else
						converged2 = converged1;
				}
				else
				{
					charact1 = 0; //strict minimum is assumed
					converged2 = converged1;
				}
					
			}
		}

		///////////////////////////////VERIFICATIONS//////////////////////////////////
		Matrix nA(3);
		Matrix nB(3);
		c_data->return_value[ip] = VerifyConvectiveRange(*cNR2[ip]);

		//Strict minimum
		if (charact1 == 0)
		{
			if (VerifyConvectiveRange(*cNR1[ip]) == 2)
			{
				c_data->return_value[ip] = 2;	//not and strong
				if (write_report)
					fprintf(f_TR_report[ip], "Candidate is not Strong. Strict minimum far from range of interest.\n");
			}
			//Negative gap
			if (Gap(*cNR1[ip], false, nA, nB) < 0.0)
			{
				c_data->return_value[ip] = 2;	//not and strong
				if (write_report)
					fprintf(f_TR_report[ip], "Candidate is is not Strong. Strict minimum with negative gap\n");
			}
		}
		//Avaliando if and candidato forte cases with intersection
		if (charact1 == 1 && converged1 == true)
		{
			surf1->NormalExt(&(*cNR1[ip])(0, 0), &(*cNR1[ip])(1, 0), &nA);
			surf2->NormalExt(&(*cNR1[ip])(2, 0), &(*cNR1[ip])(3, 0), &nB);
			if (dot(nA,nB)>= 0.0)
			{
				c_data->return_value[ip] = 2;	//not and strong
				if (write_report)
				{
					fprintf(f_TR_report[ip], "Candidate is not Strong. nA.nB is %.6f.\n", dot(nA, nB));
					fclose(f_TR_report[ip]);
				}
				return;
			}
			if (converged2 == true)
			{
				if (c_data->return_value[ip] == 2)
				{
					if (write_report)
						fprintf(f_TR_report[ip], "Candidate is not Strong. Saddle far from range of interest.\n");
				}
			}
		}

		//Saves in the variables
		c_data->convective[ip][0] = (*cNR2[ip])(0, 0);
		c_data->convective[ip][1] = (*cNR2[ip])(1, 0);
		c_data->convective[ip][2] = (*cNR2[ip])(2, 0);
		c_data->convective[ip][3] = (*cNR2[ip])(3, 0);
	}

	//Performs the degeneration automatica with base in the result of the LCP initial. Acumula `the degeneration imposta in the start.
	AutomaticDegenerationProcedure();

	if (write_report)
	{
		for (int ip = 0; ip < n_pointwise; ip++)
			fclose(f_TR_report[ip]);
	}
		
}

//Determination via optimization of the convective coordinates of the pair of surfaces
void SurfacePair::SolveLCP(SSContactData* c_data)
{
		
	///////////////////////////////Pointers and variables to make access easier//////////////////////////////////////////
	Surface* surf1 = db.surfaces[surf1_ID - 1];		//Pointer to the surface 1
	Surface* surf2 = db.surfaces[surf2_ID - 1];		//Pointer to the surface 2

	if (write_report)
	{
		for (int ip = 0; ip < n_pointwise; ip++)
		{
			OpenTRReport(ip);
			fprintf(f_TR_report[ip], "///////////////////////////////////////////////////////////////////////\n");
			fprintf(f_TR_report[ip], "\nTime\t%.6f\tIteration\t%d\t", db.last_converged_time + db.current_time_step, db.current_iteration_number);
			fprintf(f_TR_report[ip], "SolveLCP\n");
		}
		
	}

	//Case especial - Special LCP
	if (specialLCP)
	{
		//CASE 1 - not and strong candidate - return imediato:
		if (c_data->return_value[0] == 2)
		{
			if (write_report)
				fprintf(f_TR_report[0], "Candidate is not Strong\n");
		}
		else
		{
			FlexibleArcExtrusion_1_RigidArcRevolution_1* temp = static_cast<FlexibleArcExtrusion_1_RigidArcRevolution_1*>(this);
			temp->SpecialLCP(*cNR2[0]);
			c_data->convective[0][0] = (*cNR2[0])(0, 0);
			c_data->convective[0][1] = (*cNR2[0])(1, 0);
			c_data->convective[0][2] = (*cNR2[0])(2, 0);
			c_data->convective[0][3] = (*cNR2[0])(3, 0);

			//Fills return value with 0, 2 or 4
			c_data->return_value[0] = VerifyConvectiveRange(*cNR2[0]);
			//Correction (for not creating a not-strong candidate here!):
			if (c_data->return_value[0] == 2)
				c_data->return_value[0] = 4;
			if (write_report)
			{
				fprintf(f_TR_report[0], "Solution: %.12e\t%.12e\t%.12e\t%.12e\n", (*cNR2[0])(0, 0), (*cNR2[0])(1, 0), (*cNR2[0])(2, 0), (*cNR2[0])(3, 0));
				fprintf(f_TR_report[0], "SolveLCP return value is %d \n", c_data->return_value[0]);
			}
				
			EvaluateInvertedHessian(c_data);	//prepara inversa of the Hessian for envio for the routine of contact - this function already does the tratamento for not convexidade
		}
		if (write_report)
		{
			for (int i = 0; i < c_data->n_solutions; i++)
				fclose(f_TR_report[i]);
		}
		return;
	}


	//Values saved through the function (indicativos of the that ocorreu in its execution)
	//0 - Converged and this in the range of interest for contact and not there was grande change of position of the point material of contact, in the case of haver contact previous - OK
	//1 - Not there was convergence (can or not to be in the range for contact) - and and strong_candidate. Return problematico responsavel by veto of step of the simulation.
	//2 - Not and strong_candidate - OK
	//3 - change abrupta in the position of contact. Return problematico responsavel by veto of step of the simulation.
	//4 - There was convergence, but not this in the range for contact. Possivelmente element vizinho. - OK

	//Sweep of the solutions active
	for (int i = 0; i < c_data->n_solutions; i++)
	{
		//CASE 1 - not and strong candidate - return imediato:
		if (c_data->return_value[i] == 2)
		{
			if (write_report)
				fprintf(f_TR_report[i], "Candidate is not Strong\n");
		}
		else
		{
			seq_number = i;//for control of the report
			//If not for solution repetida (or be, if for active)
			if (c_data->repeated[i] == false)
			{
				//CASE 2: not there was convergence or there was convergence for a solution longe of the desired
				//Takes as initial guess the last solution correct conhecida (copies of the step converged previous of the phase 2)
				if (c_data->return_value[i] == 1 || c_data->return_value[i] == 3)
				{
					(*cNR1[i])(0, 0) = c_data->copy_convective[i][0];
					(*cNR1[i])(1, 0) = c_data->copy_convective[i][1];
					(*cNR1[i])(2, 0) = c_data->copy_convective[i][2];
					(*cNR1[i])(3, 0) = c_data->copy_convective[i][3];
					//If before not era strong candidate, the copy_convective not traz info valida. Thus, calculates the initial guess (geometrico) novamente
					if (c_data->copy_return_value[i] == 2)
					{
						InitialGuess(c_data); //Performs initial guess with criterion geometrico -> writes in 'convective'
						for (int ii = 0; ii < 4; ii++)
							(*cNR1[i])(ii, 0) = c_data->convective[i][ii];
					}
				}
				//CASE 3: there was convergence
				//Takes as initial guess the last solution converged
				if (c_data->return_value[i] == 0 || c_data->return_value[i] == 4)
				{
					(*cNR1[i])(0, 0) = c_data->convective[i][0];
					(*cNR1[i])(1, 0) = c_data->convective[i][1];
					(*cNR1[i])(2, 0) = c_data->convective[i][2];
					(*cNR1[i])(3, 0) = c_data->convective[i][3];
					//*cNR2[i] = *cNR1[i];
				}

				bool converged1 = false;
				bool converged2 = false;
				int info = 0;

				////////////////NO DEGENERATION///////////////
				if (c_data->degenerated[i] == false)
				{
					//If the contact era previamente active or vizinho
					if ((c_data->return_value[i] == 0 || c_data->return_value[i] == 4) && c_data->g_n[i] < 0.0)
					{
						if (write_report)
							fprintf(f_TR_report[i], "Direct search for saddle\n");
						//Attempt of achar direto the saddle
						*cNR2[i] = *cNR1[i];
						converged2 = FindSaddleSolution(c_data, cNR2[i], info, false);
					}
					//If not entrou in the if previous or if not converged the if previous
					if (converged2 == false)
					{
						if (write_report)
							fprintf(f_TR_report[i], "Search for minimum\n");
						//Determination of minimum or intersection
						converged1 = FindMinimumSolution(c_data, cNR1[i], info);
						*cNR2[i] = *cNR1[i];
						//If for intersection - search the point of saddle
						if (CharacterizeCriticalPoint(cNR2[i]) == 1)
						{
							converged2 = FindSaddleSolution(c_data, cNR2[i], info, false);
						}
						//If not, returns the minimum
						else
						{
							converged2 = converged1;
						}
					}
				}
				////////////////DEGENERATION///////////////
				else
				{
					//Correction in case of degeneration
					Matrix temp_coordinates(4);
					for (int index = 0; index < 4; index++)
						temp_coordinates(index, 0) = (*cNR1[i])(index, 0);
					//Transforming into degeneration basis
					temp_coordinates = transp(*c_data->P[i])*temp_coordinates;
					//Fixing degenerated coordinates
					for (int index = 0; index < 4; index++)
						if (c_data->deg_control[i][index] == true)
							temp_coordinates(index, 0) = c_data->copy_deg_coordinates[i][index];
					//Transforming into original basis
					temp_coordinates = (*c_data->P[i])*temp_coordinates;
					//Copying info into cNR1 vector
					for (int index = 0; index < 4; index++)
						(*cNR1[i])(index, 0) = temp_coordinates(index, 0);

					//If the contact era previamente active, vizinho or if not era strong (or be, if there is initial guess calculated in this point)
					if ((c_data->return_value[i] == 0 || c_data->return_value[i] == 4) && c_data->g_n[i] < 0.0)
					{
						if (write_report)
							fprintf(f_TR_report[i], "Direct search for saddle or minimum\n");
						*cNR2[i] = *cNR1[i];
						if (c_data->P_0[i]->getColumns() == 3)
							converged2 = FindSaddleSolutionDegenerated(c_data, c_data->P_0[i], cNR2[i], false);
						else
							converged2 = FindMinimumSolutionDegenerated(c_data, c_data->P_0[i], cNR2[i]);
					}
					//If not entrou in the if previous or if not converged the saddle of the if previous
					if (converged2 == false)
					{
						//Determines minimum (or intersection)
						converged1 = FindMinimumSolutionDegenerated(c_data, c_data->P_0[i], cNR1[i]);
						*cNR2[i] = *cNR1[i];
						if (c_data->P_0[i]->getColumns() == 3)
						{
							//Search of the point of saddle
							if (CharacterizeCriticalPointDegenerated(cNR2[i], c_data->P_0[i]) == 1)
								converged2 = FindSaddleSolutionDegenerated(c_data, c_data->P_0[i], cNR2[i], false);
							else
								converged2 = converged1;
						}
						else
							converged2 = converged1;
					}
				}
				if (converged2 == true)
				{
					//Next convective coordinates:
					c_data->convective[i][0] = (*cNR2[i])(0, 0);
					c_data->convective[i][1] = (*cNR2[i])(1, 0);
					c_data->convective[i][2] = (*cNR2[i])(2, 0);
					c_data->convective[i][3] = (*cNR2[i])(3, 0);

					//Fills return value with 0, 2 or 4
					c_data->return_value[i] = VerifyConvectiveRange(*cNR2[i]);
					//Correction (for not creating a not-strong candidate here!):
					if (c_data->return_value[i] == 2)
						c_data->return_value[i] = 4;

					//Test - change of solution in the propria surface
					if (c_data->copy_return_value[i] == 0)
					{
						//Previous convective coordinates:
						(*cNR1[i])(0, 0) = c_data->copy_convective[i][0];
						(*cNR1[i])(1, 0) = c_data->copy_convective[i][1];
						(*cNR1[i])(2, 0) = c_data->copy_convective[i][2];
						(*cNR1[i])(3, 0) = c_data->copy_convective[i][3];
						////Check of normal 
						Matrix cur_n1(3);
						Matrix cur_n2(3);
						surf1->NormalExt(&(*cNR2[i])(0, 0), &(*cNR2[i])(1, 0), &cur_n1);
						surf2->NormalExt(&(*cNR2[i])(2, 0), &(*cNR2[i])(3, 0), &cur_n2);
						Matrix prev_n1(3);
						Matrix prev_n2(3);
						surf1->NormalExt(&(*cNR1[i])(0, 0), &(*cNR1[i])(1, 0), &prev_n1);
						surf2->NormalExt(&(*cNR1[i])(2, 0), &(*cNR1[i])(3, 0), &prev_n2);
						if (dot(cur_n1, prev_n1) < 0 || dot(cur_n2, prev_n2) < 0)
						{
							c_data->return_value[i] = 3;

							/*db.myprintf("Code 3. n1.n1 = %.6e, n2.n2 = %.6e\n", dot(cur_n1, prev_n1), dot(cur_n2, prev_n2));
							db.myprintf("cur: %.6e,  %.6e,  %.6e , %.6e\n", (*cNR2[i])(0, 0), (*cNR2[i])(1, 0), (*cNR2[i])(2, 0), (*cNR2[i])(3, 0));
							db.myprintf("pre: %.6e,  %.6e,  %.6e , %.6e\n", (*cNR1[i])(0, 0), (*cNR1[i])(1, 0), (*cNR1[i])(2, 0), (*cNR1[i])(3, 0));
							db.myprintf("Surf1 %d Surf2 %d, charac1 %d, charact2 %d\n", surf1_ID, surf2_ID,
								CharacterizeCriticalPointDegenerated(cNR1[i], c_data->P_0[i], true),
								CharacterizeCriticalPointDegenerated(cNR2[i], c_data->P_0[i], true));*/
						}
					}
					if (write_report)
						fprintf(f_TR_report[i], "SolveLCP return value is %d \n", c_data->return_value[i]);
				}
				else
				{
					c_data->return_value[i] = 1;
					if (write_report)
						fprintf(f_TR_report[i], "LCP between surfaces %d and %d has presented problems. Code 1.\n", surf1_ID, surf2_ID);
					//db.myprintf("LCP between surfaces %d and %d has presented problems. Code 1.\n", surf1_ID, surf2_ID);
					
					//Writing report based on divergence (for monitoring issues)
					if (write_report_diverged == true && write_report == false)
					{
						write_report = true;
						for (int ii = 0; ii < n_pointwise; ii++)
						{
							InitializeTRReport(ii);
							OpenTRReport(ii);
							fprintf(f_TR_report[ii], "///////////////////////////////////////////////////////////////////////\n");
							fprintf(f_TR_report[ii], "\nTime\t%.6f\tIteration\t%d\t", db.last_converged_time + db.current_time_step, db.current_iteration_number);
							fprintf(f_TR_report[ii], "Write Report automatically launched on divergence of LCP\n");
						}
					}
					
				}

			}
		}
	}
	EvaluateInvertedHessian(c_data);	//prepara inversa of the Hessian for envio for the routine of contact - this function already does the tratamento for not convexidade
	if (write_report)
	{
		for (int i = 0; i < c_data->n_solutions; i++)
			fclose(f_TR_report[i]);
	}
		
}

bool SurfacePair::EndStepCheck(SSContactData* c_data)
{

	if (write_report)
	{
		for (int i = 0; i < c_data->n_solutions; i++)
		{
			OpenTRReport(i);
			fprintf(f_TR_report[i], "///////////////////////////////////////////////////////////////////////\n");
			fprintf(f_TR_report[i], "\nTime\t%.6f\tIteration\t%d\t", db.last_converged_time + db.current_time_step, db.current_iteration_number);
			fprintf(f_TR_report[i], "EndStepCheck\n");
		}
	}
	///////////////////////////////Pointers and variables to make access easier//////////////////////////////////////////
	Surface* surf1 = db.surfaces[surf1_ID - 1];		//Pointer to the surface 1
	Surface* surf2 = db.surfaces[surf2_ID - 1];		//Pointer to the surface 2

	//Return with problems: true
	//Return without problems: false

	//Test for evaluate solutions not convergidas
	//Sweep of the solutions active
	for (int i = 0; i < c_data->n_solutions; i++)
	{
		//If not for solution repetida (or be, if for solution active)
		if (c_data->repeated[i] == false)
		{
			//If some of the pairs active apresentou divergence of the method of optimization
			if (c_data->return_value[i] == 1)
			{
				//Test of bounding box
				/*double inflation_factor = 1.0;
				bool box = BoxOverlap(surf1->box, surf2->box, inflation_factor);
				if (box == false)
				{
					db.myprintf("LCP diverged but no bounding box overlap at the end of the step for surfaces %d and %d.\n", surf1_ID, surf2_ID);
					if (write_report)
					{
						fprintf(f_TR_report, "No bounding box overlap at the end of the step for surfaces %d and %d.\n", surf1_ID, surf2_ID);
						fclose(f_TR_report);
					}
					return false;
				}*/

				db.myprintf("LCP between surfaces %d and %d has presented problems. Code 1.\n", surf1_ID, surf2_ID);
				if (write_report)
				{
					fprintf(f_TR_report[i], "LCP between surfaces %d and %d has presented problems. Code 1.\n", surf1_ID, surf2_ID);
					for (int ip = 0; ip < c_data->n_solutions; ip++)
						fclose(f_TR_report[ip]);
				}
				return true;
			}
			//If there was muita change of range of the convective coordinates and indicativo of perda of the solution desired for c_bar
			if (c_data->return_value[i] == 3)
			{
				db.myprintf("LCP between surfaces %d and %d has presented problems. Code 3.\n", surf1_ID, surf2_ID);
				if (write_report)
				{
					fprintf(f_TR_report[i], "LCP between surfaces %d and %d has presented problems. Code 3.\n", surf1_ID, surf2_ID);
					for (int ip = 0; ip < c_data->n_solutions; ip++)
						fclose(f_TR_report[ip]);
				}
				return true;
			}
		}
	}


	if (write_report)
	{
		for (int ip = 0; ip < c_data->n_solutions; ip++)
			fclose(f_TR_report[ip]);
	}
		
	return false;
}

int SurfacePair::CharacterizeCriticalPoint(Matrix* solution)
{
	///////////////////////////////Pointers and variables to make access easier//////////////////////////////////////////
	Surface* surf1 = db.surfaces[surf1_ID - 1];		//Pointer to the surface 1
	Surface* surf2 = db.surfaces[surf2_ID - 1];		//Pointer to the surface 2

	//0 - minimum strict
	//1 - minimum not strict (intersection)
	//2 - transition minimum-saddle (just-touch)
	//3 - saddle 2 negative eigenvalues
	//4 - other
	Matrix xk(4);
	Matrix Gra(4);
	Matrix Hes(4, 4);
	Matrix P(4, 4);
	Matrix D(4, 4);
	for (int i = 0; i < 4; i++)
		xk(i, 0) = (*solution)(i, 0);
	double ob = ObjectivePhase1(xk);
	//Hessian of the objective function
	HessianPhase1(xk, Hes);
	//Calculating directions main and curvatures main of the objective function
	fulleigen1(Hes, P, D, tol_eig);

	double max_eig = -1e100;
	for (int i = 0; i < 4; i++)
	{
		if (D(i, i) > max_eig)
			max_eig = D(i, i);
	}
	double tol_intersect = max_eig * tol_convective*tol_convective;
	double tol_small = max_eig * tol_small_1;

	//minimum not strict (intersection)
	if (ob < tol_intersect)
	{
		if (write_report)
			fprintf(f_TR_report[seq_number], "Intersection found (tolerance %.6e). Eigenvalues are %.6e\t%.6e\t%.6e\t%.6e\n", tol_intersect, D(0, 0), D(1, 1), D(2, 2), D(3, 3));
		return 1;
	}
	//minimum strict
	if (D(0, 0) >= tol_small && D(1, 1) >= tol_small && D(2, 2) >= tol_small && D(3, 3) >= tol_small)
	{
		if (write_report)
			fprintf(f_TR_report[seq_number], "Strict minimum found. Eigenvalues are %.6e\t%.6e\t%.6e\t%.6e\n", D(0, 0), D(1, 1), D(2, 2), D(3, 3));
		return 0;
	}
	////transition minimum-saddle (just-touch)
	//if (D(0, 0) >= -tol_small && D(0, 0) <= +tol_small && D(1, 1) >= -tol_small && D(1, 1) <= +tol_small && D(2, 2) >= tol_small && D(3, 3) >= tol_small)
	//{
	//	if (write_report)
	//		fprintf(f_TR_report, "Just-touch found. Eigenvalues are %.6e\t%.6e\t%.6e\t%.6e\n", D(0, 0), D(1, 1), D(2, 2), D(3, 3));
	//	return 2;
	//}
	//transition minimum-saddle (just-touch)
	if (D(0, 0) <= tol_small && D(1, 1) <= tol_small && D(2, 2) >= tol_small && D(3, 3) >= tol_small)
	{
		if (write_report)
			fprintf(f_TR_report[seq_number], "Saddle found. Eigenvalues are %.6e\t%.6e\t%.6e\t%.6e\n", D(0, 0), D(1, 1), D(2, 2), D(3, 3));
		return 3;
	}
	if (write_report)
		fprintf(f_TR_report[seq_number], "Other critical point found. Eigenvalues are %.6e\t%.6e\t%.6e\t%.6e. Tolerance is %.6e.\n", D(0, 0), D(1, 1), D(2, 2), D(3, 3), tol_small);
	return 4;
}

int SurfacePair::CharacterizeCriticalPointDegenerated(Matrix* solution, Matrix* P_0, bool print)
{
	///////////////////////////////Pointers and variables to make access easier//////////////////////////////////////////
	Surface* surf1 = db.surfaces[surf1_ID - 1];		//Pointer to the surface 1
	Surface* surf2 = db.surfaces[surf2_ID - 1];		//Pointer to the surface 2

	//0 - minimum strict
	//1 - minimum not strict (intersection)
	//2 - transition minimum-saddle (just-touch)
	//3 - saddle 2 negative eigenvalues
	//4 - other
	Matrix xk(4);
	Matrix Hes(4, 4);
	Matrix dHes(3, 3);
	Matrix P(3, 3);
	Matrix D(3, 3);
	for (int i = 0; i < 4; i++)
		xk(i, 0) = (*solution)(i, 0);
	double ob = ObjectivePhase1(xk);
	//Hessian of the objective function
	HessianPhase1(xk, Hes);
	dHes = transp(*P_0)*Hes*(*P_0);
	fulleigen1(dHes, P, D, tol_eig);
	if (print)
	{
		if (dHes.getLines() == 3)
			db.myprintf("CharacterizeCriticalPointDegenerated:\nEingenvalues are %.6e  %.6e  %.6e\n", D(0, 0), D(1, 1), D(2, 2));
		if (dHes.getLines() == 2)
			db.myprintf("CharacterizeCriticalPointDegenerated:\nEingenvalues are %.6e  %.6e\n", D(0, 0), D(1, 1));
	}

	double max_eig = -1e100;
	for (int i = 0; i < 3; i++)
	{
		if (D(i, i) > max_eig)
			max_eig = D(i, i);
	}
	double tol_intersect = max_eig * tol_convective*tol_convective;
	double tol_small = max_eig * tol_small_1;

	//minimum not strict (intersection)
	if (ob < tol_intersect)
	{
		if (write_report)
			fprintf(f_TR_report[seq_number], "Intersection found. Eigenvalues are %.6e\t%.6e\t%.6e\n", D(0, 0), D(1, 1), D(2, 2));
		return 1;
	}
	//minimum strict
	if (D(0, 0) >= tol_small && D(1, 1) >= tol_small && D(2, 2) >= tol_small)
	{
		if (write_report)
			fprintf(f_TR_report[seq_number], "Strict minimum found. Eigenvalues are %.6e\t%.6e\t%.6e. Tolerance is %.6e.\n", D(0, 0), D(1, 1), D(2, 2), tol_small);
		return 0;
	}
	////transition minimum-saddle (just-touch)
	//if (D(0, 0) >= -tol_small && D(0, 0) <= +tol_small && D(1, 1) >= tol_small && D(2, 2) >= tol_small)
	//{
	//	if (write_report)
	//		fprintf(f_TR_report, "Just-touch found. Eigenvalues are %.6e\t%.6e\t%.6e\n", D(0, 0), D(1, 1), D(2, 2));
	//	return 2;
	//}
	//saddle
	if (D(0, 0) <= tol_small && D(1, 1) >= tol_small && D(2, 2) >= tol_small)
	{
		if (write_report)
			fprintf(f_TR_report[seq_number], "Saddle found. Eigenvalues are %.6e\t%.6e\t%.6e\n", D(0, 0), D(1, 1), D(2, 2));
		return 3;
	}
	if (write_report)
		fprintf(f_TR_report[seq_number], "Other critical point found. Eigenvalues are %.6e\t%.6e\t%.6e\n", D(0, 0), D(1, 1), D(2, 2));
	return 4;
}

//Optimization - determination of minimum
bool SurfacePair::FindMinimumSolution(SSContactData* c_data, Matrix* solution, int &return_info)
{
	//Data - trust region
	double Deltamax = 1e4;			//maximum radius of the trust region allowed
	double Deltak = 0.1;			//current radius of trust region
	double etha = 0.15;				//value between 0 and 0.15 - indicates that the approximation and ruim and veta the increment
	double rhok = 0.0;				//reason between the difference in the objective function and the difference of the approximation used
	double actual_reduction = 0.0;
	double predicted_reduction = 0.0;
	int max_it = max_it_1;
	double last_reduction = 0.0;
	double reduction = 0.0;
	Matrix Hes(4, 4);
	Matrix Gra(4, 1);
	Matrix pGra(4, 1);
	Matrix pk(4);
	Matrix pb(4);
	Matrix pc(4);
	Matrix xk(4);
	Matrix P(4, 4);
	Matrix D(4, 4);
	Matrix cHes(4, 4);

	//Initialization of the method - initial guess
	for (int i = 0; i < 4; i++)
		xk(i, 0) = (*solution)(i, 0);

	//Stopping criterion
	HessianPhase1(xk, Hes);
	fulleigen1(Hes, P, D, tol_eig);
	double max_eig = -1e100;
	for (int i = 0; i < 4; i++)
	{
		if (D(i, i) > max_eig)
			max_eig = D(i, i);
	}
	double tol_ortho = tol_convective * abs(max_eig);
	double tol_small = tol_small_1;

	//Criterion for identify a intersection
	double tol_intersect = max_eig * tol_convective*tol_convective;

	if (write_report)
		fprintf(f_TR_report[seq_number], "FindMinimumSolution\n");
	char c = 'I';

	int it = 1;
	//Objective
	double ob = ObjectivePhase1(xk);
	//Gradient
	GradientPhase1(xk, Gra);
	//Hessian
	HessianPhase1(xk, Hes);
	//Error - forcing first input
	double error = tol_ortho + 1.0;
	//Initial guess report
	if (write_report)
		fprintf(f_TR_report[seq_number], "%d\t%.12e\t%.12e\t%.12e\t%.12e\t%.6e\t%.6e\t%.6e\t%c\t%.6e\t%.6e\t%.6e\n", it, xk(0, 0), xk(1, 0), xk(2, 0), xk(3, 0), ob, error, Deltak, c, rhok, actual_reduction, predicted_reduction);
	/////////////////////////////////////////////////////BEGIN///////////////////////////////////////////////
	
	while ((error > tol_ortho || (norm(pk) > tol_convective && ob > tol_intersect)) && it <= max_it)
	{
		//Determination of the point of Cauchy
		double gragra = (transp(Gra)*Gra)(0, 0);
		double grahesgra = (transp(Gra)*Hes*Gra)(0, 0);
		double normgra = norm(Gra);
		Matrix pc;	//direction of the Cauchy point
		/////////////////Point of Cauchy///////////////////////
		if (grahesgra <= 0.0)
			pc = -(Deltak / normgra)*Gra;
		else
			pc = -gragra / grahesgra * Gra;
		double normpc = norm(pc);
		//Cauchy point outside the TR - use a fraction of it
		if ((normpc + tol_convective) >= Deltak)
		{
			pk = (Deltak / normpc)*pc; //Steep descent
			c = 'C';
		}
		else//Cauchy point inside the TR
		{
			//Calculating directions main and curvatures main of the objective function
			cHes = Hes;
			fulleigen1(cHes, P, D, tol_eig);
			//Writing gradient in the directions main
			pGra = transp(P)*Gra;
			//Determination of the smallest eigenvalue (min_eig)
			double min_eig = 1e100;
			for (int i = 0; i < 4; i++)
			{
				if (D(i, i) < min_eig)
					min_eig = D(i, i);
			}
			//Construction of the direction of search
			//Direction of search based in NR - modificada through the smallest eigenvalue
			zeros(&pb);
			//If the smallest eigenvalue and smaller or igual the zero (tol_small) - modifies the direction of NR for garantir direction descendente
			if (min_eig < tol_small)
			{
				for (int i = 0; i < 4; i++)
					pb(i, 0) = -pGra(i, 0) / (D(i, i) - (min_eig - abs(min_eig)*tol_ascent));
			}
			//If the smallest eigenvalue and larger that zero (tol_small) - direction of NR and chosen
			else
			{
				for (int i = 0; i < 4; i++)
					pb(i, 0) = -pGra(i, 0) / D(i, i);
			}
			//Writing direction of search in the coordinates originais

			pb = P * pb;
			double normpb = norm(pb);
			double thetak;
			//Newton point inside the TR -  use it
			if (normpb <= Deltak)
			{
				pk = pb;
				c = 'N';
			}
			else
			{
				////////////////////////////Dogleg path//////////////////////////
				double a1, b1, c1;
				a1 = norm(pb - pc)*norm(pb - pc);
				b1 = 2 * (transp(pc)*(pb - pc))(0, 0);
				c1 = normpc * normpc - Deltak * Deltak;
				thetak = (-b1 + sqrt(b1*b1 - 4 * a1*c1)) / (2 * a1);
				//Determination of the path
				pk = pc + thetak * (pb - pc);
				c = 'D';
			}
		}

		//////////////////////////UPDATING SOLUTION////////////////////////////////
		//Calculation of rhok
		actual_reduction = ObjectivePhase1(xk) - ObjectivePhase1(xk + pk);
		predicted_reduction = -(transp(Gra)*pk + 0.5*transp(pk)*Hes*pk)(0, 0);
		rhok = actual_reduction / predicted_reduction;

		if (abs(actual_reduction / ObjectivePhase1(xk)) < tol_ascent) 
			rhok = 1.0;
		
		if (abs(predicted_reduction) < tol_small || abs(actual_reduction) < tol_small)
			rhok = 1.0;
			
		if (rhok >= 0)
		{
			//low reduction
			if (rhok < 0.25)
				Deltak = 0.25*norm(pk);
			else
			{
				if (rhok > 0.75 && ((norm(pk) + tol_convective) >= Deltak && (norm(pk) - tol_convective) <= Deltak))//high reduction and testing the limits of the trust region
				{
					//augments the radius of TR
					if (2.0 * Deltak < Deltamax)
						Deltak = 2.0 * Deltak;
					else
						Deltak = Deltamax;
				}
			}
			if (rhok >= etha)
				xk = xk + pk;
		}
		else
			Deltak = Deltak / 2.0;
		//Increments iterations
		it++;
		//Objective
		ob = ObjectivePhase1(xk);
		//Gradient
		GradientPhase1(xk, Gra);
		//Hessian
		HessianPhase1(xk, Hes);
		//Error - norm of the gradient
		error = norm(Gra);
		if (write_report)
			fprintf(f_TR_report[seq_number], "%d\t%.12e\t%.12e\t%.12e\t%.12e\t%.6e\t%.6e\t%.6e\t%c\t%.6e\t%.6e\t%.6e\n", it, xk(0, 0), xk(1, 0), xk(2, 0), xk(3, 0), ob, error, Deltak, c, rhok, actual_reduction, predicted_reduction);
	}
	if (write_report)
		fprintf(f_TR_report[seq_number], "Error\t%.6e\tTolerance\t%.6e\n", error, tol_ortho);

	//Return of the function
	if (error < tol_ortho)
	{
		for (int i = 0; i < 4; i++)
			(*solution)(i, 0) = xk(i, 0);
		return true;
	}
	else
		return false;
}

//Optimization - determination of saddle
bool SurfacePair::FindSaddleSolution(SSContactData* c_data, Matrix* solution, int &return_info, bool return_gap)
{
	//Data - trust region
	double Deltamax = 10;			//maximum radius of the trust region allowed
	double Deltak = 0.01;			//current radius of trust region
	double etha = 0.15;				//value between 0 and 0.15 - indicates that the approximation and ruim and veta the increment
	double rhok = 0.0;				//reason between the difference in the objective function and the difference of the approximation used
	double actual_change = 0.0;
	double predicted_change = 0.0;
	double tol_small = tol_small_1;
	int max_it = max_it_2;			//max iterations
	//Variables internal of use geral
	Matrix Hes(4, 4);
	Matrix Gra(4, 1);
	Matrix pGra(4, 1);
	Matrix P(4, 4);
	Matrix D(4, 4);
	Matrix cHes(4, 4);
	Matrix pk(4);
	Matrix pb(4);
	Matrix xk(4);
	//Initialization of the method - initial guess - obtained of the solution previous - problem of minimum distance
	for (int i = 0; i < 4; i++)
		xk(i, 0) = (*solution)(i, 0);

	//Error of the minimum distance
	GradientPhase1(xk, Gra);
	double ob = ObjectivePhase1(xk);
	//Hessian of the objective function
	HessianPhase1(xk, Hes);
	//Calculating directions main and curvatures main of the objective function
	cHes = Hes;
	fulleigen1(cHes, P, D, tol_eig);

	//Stopping criterion
	double max_eig = -1e100;
	for (int i = 0; i < 4; i++)
	{
		if (D(i, i) > max_eig)
			max_eig = D(i, i);
	}
	double tol_ortho = tol_convective * abs(max_eig);

	bool saddle_basin = false;		//variable boolean that indicates that the bacia of attraction of the saddle was localizada
	int it = 1;
	bool conv1 = false;
	saddle_basin = false;

	if (CharacterizeCriticalPoint(&xk) == 3)
	{
		if (write_report)
			fprintf(f_TR_report[seq_number], "FindSaddleSolution initial guess is already in saddle basin\n");
		saddle_basin = true;//saddle basin
	}
	else
	{
		if (write_report)
			fprintf(f_TR_report[seq_number], "FindSaddleSolution initial guess is not in saddle basin\n");

		Surface* surf1 = db.surfaces[surf1_ID - 1];		//Pointer to the surface 1
		Surface* surf2 = db.surfaces[surf2_ID - 1];		//Pointer to the surface 2
		Matrix nA(3);
		surf1->NormalExt(&xk(0, 0), &xk(1, 0), &nA);
		Matrix nB(3);
		surf2->NormalExt(&xk(2, 0), &xk(3, 0), &nB);
		bool normal_fixed = true;

		//Degeneration local
		Matrix P_0(4, 2);
		for (int i = 0; i < 4; i++)
			P_0(i, 0) = P(i, 0);//first eigenvector
		for (int i = 0; i < 4; i++)
			P_0(i, 1) = P(i, 1);//second eigenvector
		return_info = 0;
		//Search bi-direcional
		conv1 = FindMinimumGapDegenerated(c_data, &P_0, &xk, return_info, normal_fixed, nA, nB);

		/////
		//Hessian of the objective function
		HessianPhase1(xk, Hes);
		//Calculating directions main and curvatures main of the objective function
		cHes = Hes;
		fulleigen1(cHes, P, D, tol_eig);
		//Degeneration local
		for (int i = 0; i < 4; i++)
			P_0(i, 0) = P(i, 0);//first eigenvector
		for (int i = 0; i < 4; i++)
			P_0(i, 1) = P(i, 1);//second eigenvector
		conv1 = FindMinimumGapDegenerated(c_data, &P_0, &xk, return_info, normal_fixed, nA, nB);
		/////

		/////
		//Hessian of the objective function
		HessianPhase1(xk, Hes);
		//Calculating directions main and curvatures main of the objective function
		cHes = Hes;
		fulleigen1(cHes, P, D, tol_eig);
		//Degeneration local
		for (int i = 0; i < 4; i++)
			P_0(i, 0) = P(i, 0);//first eigenvector
		for (int i = 0; i < 4; i++)
			P_0(i, 1) = P(i, 1);//second eigenvector
		conv1 = FindMinimumGapDegenerated(c_data, &P_0, &xk, return_info, normal_fixed, nA, nB);
		/////


		if (conv1 == false)
			return false;
		else
		{
			if (return_gap == true)
			{
				//Saves result in solution
				for (int i = 0; i < 4; i++)
					(*solution)(i, 0) = xk(i, 0);
				return true;
			}
		}
	}
	if (CharacterizeCriticalPoint(&xk) == 3)
	{
		if (write_report)
			fprintf(f_TR_report[seq_number], "After FindMinimumGap, initial guess is in saddle basin\n");
	}
	else
	{
		if (write_report)
			fprintf(f_TR_report[seq_number], "After FindMinimumGap, initial guess is not in saddle basin\n");
	}

	//Calculations of error and objective
	GradientPhase1(xk, Gra);
	//Error - forcing first input
	double error = tol_ortho + 1.0;
	ob = ObjectivePhase1(xk);
	//Hessian
	HessianPhase1(xk, Hes);
	int flag_error = 0;
	if (write_report)
		fprintf(f_TR_report[seq_number], "FindSaddleSolution\n");
	char c = 'I';
	if (write_report)
		fprintf(f_TR_report[seq_number], "%d\t%.12e\t%.12e\t%.12e\t%.12e\t%.6e\t%.6e\t%.6e\t%c\t%.6e\t%.6e\t%.6e\n", it, xk(0, 0), xk(1, 0), xk(2, 0), xk(3, 0), ob, error, Deltak, c, rhok, actual_change, predicted_change);
	/////////////////////////////////Step of analysis of the bacia of attraction of the point of saddle/////////////////////////////////////
	while ((error > tol_ortho || norm(pk) > tol_convective) && it <= max_it)
	{
		
		cHes = Hes;
		pGra = Gra;
		//Solves linear system
		pb = fullsystem(cHes, -1.0*pGra, &flag_error);
		//Limita TR
		pk = pb;
		//Limitando TR
		double norminc = norm(pk);
		if ((norminc + tol_convective) >= Deltak)
			pk = (Deltak / norminc)*pk;
		c = 'N';
		//////////////////////////UPDATING SOLUTION////////////////////////////////
		//Calculation of rhok
		actual_change = ObjectivePhase1(xk) - ObjectivePhase1(xk + pk);
		predicted_change = -(transp(Gra)*(pk)+0.5*transp(pk)*Hes*(pk))(0, 0);
		rhok = actual_change / predicted_change;

		if (abs(actual_change / ObjectivePhase1(xk)) < tol_ascent)
			rhok = 1.0;

		if (predicted_change < tol_small)
			rhok = 1.0;

		//update - case the model function represente bem the região. This avoids exit of the bacia of attraction of the saddle
		if (rhok >= 0.9 && rhok <= 1.1)
		{
			xk = xk + pk;
			Deltak = 2 * Deltak;
			if (Deltak > Deltamax)
				Deltak = Deltamax;
		}
		else
			Deltak = Deltak / 4;
		//increments iterations
		it++;
		//Calculations that will be used for check output of the loop - or in the next iteration, case not saia
		//Objective
		ob = ObjectivePhase1(xk);
		//Gradient
		GradientPhase1(xk, Gra);
		//Hessian
		HessianPhase1(xk, Hes);
		error = norm(Gra);
		if (write_report)
			fprintf(f_TR_report[seq_number], "%d\t%.12e\t%.12e\t%.12e\t%.12e\t%.6e\t%.6e\t%.6e\t%c\t%.6e\t%.6e\t%.6e\n", it, xk(0, 0), xk(1, 0), xk(2, 0), xk(3, 0), ob, error, Deltak, c, rhok, actual_change, predicted_change);
	}
	if (write_report)
		fprintf(f_TR_report[seq_number], "Error\t%.6e\tTolerance\t%.6e\n", error, tol_ortho);
	
	////////////////////////////////Check final - if the solution and point of saddle/////////////////////////////////////
	//Return of the function
	if (error < tol_ortho)
	{
		//Saves result in solution
		for (int i = 0; i < 4; i++)
			(*solution)(i, 0) = xk(i, 0);

		//Checks the solution type encontrada
		int character = CharacterizeCriticalPoint(&xk);
		if (character == 3)
			return true;
		else
			return false;
	}
	else
		return false;
}

//Optimization - determination of minimum
bool SurfacePair::FindMinimumSolutionDegenerated(SSContactData* c_data, Matrix* P_0, Matrix* solution)
{
	//Data - trust region
	double Deltamax = 1e4;			//maximum radius of the trust region allowed
	double Deltak = 0.1;			//current radius of trust region
	double etha = 0.15;				//value between 0 and 0.15 - indicates that the approximation and ruim and veta the increment
	double rhok = 0.0;				//reason between the difference in the objective function and the difference of the approximation used
	double actual_reduction = 0.0;
	double predicted_reduction = 0.0;
	int max_it = max_it_1;
	double last_reduction = 0.0;
	double reduction = 0.0;
	
	Matrix Hes(4, 4);
	Matrix Gra(4, 1);
	Matrix xk(4);
	int order = P_0->getColumns();
	Matrix pk(order);
	Matrix pb(order);
	Matrix pc(order);
	Matrix pGra(order, 1);
	Matrix deg_Hes(order, order);
	Matrix deg_Gra(order, 1);
	Matrix P(order, order);
	Matrix D(order, order);
	Matrix cHes(order, order);

	//Initialization of the method - initial guess - obtained of the solution previous - problem of minimum distance
	for (int i = 0; i < 4; i++)
		xk(i, 0) = (*solution)(i, 0);

	//Stopping criterion
	HessianPhase1(xk, Hes);
	deg_Hes = transp(*P_0)*Hes*(*P_0);
	fulleigen1(deg_Hes, P, D, tol_eig);
	double max_eig = -1e100;
	for (int i = 0; i < order; i++)
	{
		if (D(i, i) > max_eig)
			max_eig = D(i, i);
	}
	
	double tol_ortho = tol_convective * abs(max_eig);
	//double tol_small = tol_small_1*1e4;
	double tol_small = tol_small_1;

	if (write_report)
		fprintf(f_TR_report[seq_number], "FindMinimumSolutionDegenerated\n");
	char c = 'I';

	int it = 1;
	//Objective
	double ob = ObjectivePhase1(xk);
	//Gradient
	GradientPhase1(xk, Gra);
	//Hessian
	HessianPhase1(xk, Hes);
	//Transformações - degeneration
	deg_Gra = transp(*P_0)*Gra;
	deg_Hes = transp(*P_0)*Hes*(*P_0);
	//Error - forcing first input
	double error = tol_ortho + 1.0;
	//Initial guess report
	if (write_report)
		fprintf(f_TR_report[seq_number], "%d\t%.12e\t%.12e\t%.12e\t%.12e\t%.6e\t%.6e\t%.6e\t%c\t%.6e\t%.6e\t%.6e\n", it, xk(0, 0), xk(1, 0), xk(2, 0), xk(3, 0), ob, error, Deltak, c, rhok, actual_reduction, predicted_reduction);
	/////////////////////////////////////////////////////BEGIN///////////////////////////////////////////////
	while ((error > tol_ortho || norm((*P_0)*pk) > tol_convective) && it <= max_it)
	{
		//Determination of the point of Cauchy
		double gragra = (transp(deg_Gra)*deg_Gra)(0, 0);
		double grahesgra = (transp(deg_Gra)*deg_Hes*deg_Gra)(0, 0);
		double normgra = norm(deg_Gra);
		Matrix pc;	//direction of the Cauchy point
		/////////////////Point of Cauchy///////////////////////
		if (grahesgra <= 0.0)
			pc = -(Deltak / normgra)*deg_Gra;
		else
			pc = -gragra / grahesgra * deg_Gra;
		double normpc = norm(pc);
		//Cauchy point outside the TR - use a fraction of it
		if ((normpc + tol_convective) >= Deltak)
		{
			pk = (Deltak / normpc)*pc; //Steep descent
			c = 'C';
		}

		else//Cauchy point inside the TR
		{
			//Calculating directions main and curvatures main of the objective function
			cHes = deg_Hes;
			fulleigen1(cHes, P, D, tol_eig);
			//Writing gradient in the directions main
			pGra = transp(P)*deg_Gra;
			//Determination of the smallest eigenvalue (min_eig)
			double min_eig = 1e100;
			for (int i = 0; i < order; i++)
			{
				if (D(i, i) < min_eig)
					min_eig = D(i, i);
			}
			//Construction of the direction of search
			//Direction of search based in NR - modificada through the smallest eigenvalue
			zeros(&pb);
			//If the smallest eigenvalue and smaller or igual the zero (tol_small) - modifies the direction of NR for garantir direction descendente
			if (min_eig < tol_small)
			{
				for (int i = 0; i < order; i++)
					pb(i, 0) = -pGra(i, 0) / (D(i, i) - (min_eig - abs(min_eig)*tol_ascent));
			}
			//If the smallest eigenvalue and larger that zero (tol_small) - direction of NR and chosen
			else
			{
				for (int i = 0; i < order; i++)
					pb(i, 0) = -pGra(i, 0) / D(i, i);
			}
			//Writing direction of search in the coordinates originais
			pb = P * pb;
			double normpb = norm(pb);
			double thetak;
			//Newton point inside the TR -  use it
			if (normpb <= Deltak)
			{
				pk = pb;
				c = 'N';
			}
			else
			{
				////////////////////////////Dogleg path//////////////////////////
				double a, b, c;
				a = norm(pb - pc)*norm(pb - pc);
				b = 2 * (transp(pc)*(pb - pc))(0, 0);
				c = normpc * normpc - Deltak * Deltak;
				thetak = (-b + sqrt(b*b - 4 * a*c)) / (2 * a);
				//Determination of the path
				pk = pc + thetak * (pb - pc);
				c = 'D';
			}
		}

		//////////////////////////UPDATING SOLUTION////////////////////////////////
		//Calculation of rhok
		double actual_reduction = ObjectivePhase1(xk) - ObjectivePhase1(xk + (*P_0)*pk);
		double predicted_reduction = -(transp(deg_Gra)*pk + 0.5*transp(pk)*deg_Hes*pk)(0, 0);
		rhok = actual_reduction / predicted_reduction;

		if (abs(actual_reduction / ObjectivePhase1(xk)) < tol_ascent)
			rhok = 1.0;

		if (abs(predicted_reduction) < tol_small || abs(actual_reduction) < tol_small)
			rhok = 1.0;

		if (rhok >= 0.0)
		{
			if (rhok < 0.25)//low reduction or even augmenting the objective function
				Deltak = 0.25*norm(pk);//reduce TR
			else
			{
				if (rhok > 0.75 && ((norm(pk) + tol_convective) >= Deltak && (norm(pk) - tol_convective) <= Deltak))//high reduction and testing the limits of the trust region
				{
					//augments the radius of TR
					if (2.0 * Deltak < Deltamax)
						Deltak = 2.0 * Deltak;
					else
						Deltak = Deltamax;
				}
			}
			if (rhok >= etha)
				xk = xk + (*P_0)*pk;
		}
		else
			Deltak = Deltak / 2;
		//		//Increments iterations
		it++;
		//Objective
		ob = ObjectivePhase1(xk);
		//Gradient
		GradientPhase1(xk, Gra);
		//Hessian
		HessianPhase1(xk, Hes);
		//Transformações - degeneration
		deg_Gra = transp(*P_0)*Gra;
		deg_Hes = transp(*P_0)*Hes*(*P_0);
		//Error - norm of the gradient
		error = norm(deg_Gra);
		//Initial guess report
		if (write_report)
			fprintf(f_TR_report[seq_number], "%d\t%.12e\t%.12e\t%.12e\t%.12e\t%.6e\t%.6e\t%.6e\t%c\t%.6e\t%.6e\t%.6e\n", it, xk(0, 0), xk(1, 0), xk(2, 0), xk(3, 0), ob, error, Deltak, c, rhok, actual_reduction, predicted_reduction);
	}
	if (write_report)
		fprintf(f_TR_report[seq_number], "Error\t%.6e\tTolerance\t%.6e\n", error, tol_ortho);

	//Return of the function
	if (error < tol_ortho)
	{
		//Saves result in solution
		for (int i = 0; i < 4; i++)
			(*solution)(i, 0) = xk(i, 0);
		return true;
	}
	else
		return false;
}

//Optimization - determination of minimum
bool SurfacePair::FindMinimumGapDegenerated(SSContactData* c_data, Matrix* P_0, Matrix* solution, int &return_info, bool fixed_normals, Matrix& nA, Matrix& nB)
{
	//Data - trust region   
	double Deltak = minimum_convective_range / 1000;		//current radius of trust region
	double Deltamax = 100;									//maximum radius of the trust region allowed
	double etha = 0.15;										//value between 0 and 0.15 - indicates that the approximation and ruim and veta the increment
	double rhok = 0.0;										//reason between the difference in the objective function and the difference of the approximation used
	double last_reduction = 0.0;
	double reduction = 0.0;
	double actual_reduction = 0;
	double predicted_reduction = 0;
	double tol_small = tol_small_1;
	int max_it = max_it_2;			//Max iterations
	Matrix Hes(4, 4);
	Matrix Gra(4, 1);
	Matrix xk(4);
	int order = P_0->getColumns();
	Matrix pk(order);
	Matrix pb(order);
	Matrix pc(order);
	Matrix pGra(order, 1);
	Matrix deg_Hes(order, order);
	Matrix deg_Gra(order, 1);
	Matrix P(order, order);
	Matrix D(order, order);
	Matrix cHes(order, order);

	//	//Initialization of the method - initial guess
	for (int i = 0; i < 4; i++)
		xk(i, 0) = (*solution)(i, 0);

	//Stopping criterion
	HessianGap(xk, Hes, fixed_normals, nA, nB);
	deg_Hes = transp(*P_0)*Hes*(*P_0);
	fulleigen1(deg_Hes, P, D, tol_eig);
	double max_eig = -1e100;
	for (int i = 0; i < order; i++)
	{
		if (D(i, i) > max_eig)
			max_eig = D(i, i);
	}
	double tol_ortho = tol_convective * abs(max_eig);

	if (write_report)
		fprintf(f_TR_report[seq_number], "FindMinimumGapDegenerated\n");
	char c = 'I';

	int it = 1;
	//Objective
	double ob = Gap(xk, fixed_normals, nA, nB);
	//Gradient
	GradientGap(xk, Gra, fixed_normals, nA, nB);
	//Hessian
	HessianGap(xk, Hes, fixed_normals, nA, nB);
	deg_Gra = transp(*P_0)*Gra;
	deg_Hes = transp(*P_0)*Hes*(*P_0);
	//Error - forcing first input
	double error = tol_ortho + 1.0;
	//Initial guess report
	if (write_report)
		fprintf(f_TR_report[seq_number], "%d\t%.12e\t%.12e\t%.12e\t%.12e\t%.6e\t%.6e\t%.6e\t%c\t%.6e\t%.6e\t%.6e\n", it, xk(0, 0), xk(1, 0), xk(2, 0), xk(3, 0), ob, error, Deltak, c, rhok, actual_reduction, predicted_reduction);
	/////////////////////////////////////////////////////BEGIN///////////////////////////////////////////////
	while ((error > tol_ortho || norm((*P_0)*pk) > tol_convective) && it <= max_it)
	{
		//Determination of the point of Cauchy
		double gragra = (transp(deg_Gra)*deg_Gra)(0, 0);
		double grahesgra = (transp(deg_Gra)*deg_Hes*deg_Gra)(0, 0);
		double normgra = norm(deg_Gra);
		Matrix pc;	//direction of the Cauchy point
		/////////////////Point of Cauchy///////////////////////
		if (grahesgra <= 0.0)
			pc = -(Deltak / normgra)*deg_Gra;
		else
			pc = -gragra / grahesgra * deg_Gra;
		double normpc = norm(pc);
		//Cauchy point outside the TR - use a fraction of it
		if ((normpc + tol_convective) >= Deltak)
		{
			pk = (Deltak / normpc)*pc; //Steep descent
			c = 'C';
		}
		else//Cauchy point inside the TR
		{
			//Calculating directions main and curvatures main of the objective function
			cHes = deg_Hes;
			fulleigen1(cHes, P, D, tol_eig);
			//Writing gradient in the directions main
			pGra = transp(P)*deg_Gra;
			//Determination of the smallest eigenvalue (min_eig)
			double min_eig = 1e100;
			for (int i = 0; i < order; i++)
			{
				if (D(i, i) < min_eig)
					min_eig = D(i, i);
			}
			//Construction of the direction of search
			//Direction of search based in NR - modificada through the smallest eigenvalue
			zeros(&pb);
			//If the smallest eigenvalue and smaller or igual the zero (tol_small) - modifies the direction of NR for garantir direction descendente
			if (min_eig < tol_small)
			{
				for (int i = 0; i < order; i++)
					pb(i, 0) = -pGra(i, 0) / (D(i, i) - (min_eig - abs(min_eig)*tol_ascent));
			}
			//If the smallest eigenvalue and larger that zero (tol_small) - direction of NR and chosen
			else
			{
				for (int i = 0; i < order; i++)
					pb(i, 0) = -pGra(i, 0) / D(i, i);
			}
			//Writing direction of search in the coordinates originais
			pb = P * pb;
			double normpb = norm(pb);
			double thetak;
			//Newton point inside the TR -  use it
			if (normpb <= Deltak)
			{
				pk = pb;
				c = 'N';
			}

			else
			{
				////////////////////////////Dogleg path//////////////////////////
				double a1, b1, c1;
				a1 = norm(pb - pc)*norm(pb - pc);
				b1 = 2 * (transp(pc)*(pb - pc))(0, 0);
				c1 = normpc * normpc - Deltak * Deltak;
				thetak = (-b1 + sqrt(b1*b1 - 4 * a1*c1)) / (2 * a1);
				//Determination of the path
				pk = pc + thetak * (pb - pc);
				c = 'D';
			}
		}
		//////////////////////////UPDATING SOLUTION////////////////////////////////
		//Calculation of rhok
		actual_reduction = Gap(xk, fixed_normals, nA, nB) - Gap(xk + (*P_0)*pk, fixed_normals, nA, nB);
		predicted_reduction = -(transp(deg_Gra)*pk + 0.5*transp(pk)*deg_Hes*pk)(0, 0);
		rhok = actual_reduction / predicted_reduction;

		if (abs(predicted_reduction / Gap(xk, fixed_normals, nA, nB)) < tol_ascent)
			rhok = 1.0;

		if (abs(predicted_reduction) < tol_small || abs(actual_reduction) < tol_small)
			rhok = 1.0;

		if (rhok >= 0.0)
		{
			if (rhok < 0.25)//low reduction
				Deltak = 0.25*norm(pk);//reduce TR
			else
			{
				if (rhok > 0.75 && ((norm(pk) + tol_convective) >= Deltak && (norm(pk) - tol_convective) <= Deltak))//high reduction and testing the limits of the trust region
				{
					//augments the radius of TR
					if (2.0 * Deltak < Deltamax)
						Deltak = 2.0 * Deltak;
					else
						Deltak = Deltamax;
				}
			}
			if (rhok >= etha)
				xk = xk + (*P_0)*pk;
			//error = actual_reduction;
		}
		else
			Deltak = Deltak / 2.0;

		//Increments iterations
		it++;
		//Objective
		ob = Gap(xk, fixed_normals, nA, nB);
		//Gradient
		GradientGap(xk, Gra, fixed_normals, nA, nB);
		//Hessian
		HessianGap(xk, Hes, fixed_normals, nA, nB);
		deg_Gra = transp(*P_0)*Gra;
		deg_Hes = transp(*P_0)*Hes*(*P_0);
		//Error - norm of the gradient
		error = norm(deg_Gra);
		//Initial guess report
		if (write_report)
			fprintf(f_TR_report[seq_number], "%d\t%.12e\t%.12e\t%.12e\t%.12e\t%.6e\t%.6e\t%.6e\t%c\t%.6e\t%.6e\t%.6e\n", it, xk(0, 0), xk(1, 0), xk(2, 0), xk(3, 0), ob, error, Deltak, c, rhok, actual_reduction, predicted_reduction);
	}
	if (write_report)
		fprintf(f_TR_report[seq_number], "Error\t%.6e\tTolerance\t%.6e\n", error, tol_ortho);


	//Return of the function
	if (error < tol_ortho)
	{
		//Saves result in solution
		for (int i = 0; i < 4; i++)
			(*solution)(i, 0) = xk(i, 0);
		return true;
	}
	else
		return false;
}

//Optimization - determination of saddle
bool SurfacePair::FindSaddleSolutionDegenerated(SSContactData* c_data, Matrix* P_0, Matrix* solution, bool return_gap)
{
	//Function desenvolvida for degeneration of order 3
	//Order 4 - call function FindSaddleSolution
	//Order 2 - call function FindMinimumSolution (not and necessary determine the point of saddle)

	//Data - trust region
	double Deltamax = 10;			//maximum radius of the trust region allowed
	double Deltak = 0.01;			//current radius of trust region
	double etha = 0.15;				//value between 0 and 0.15 - indicates that the approximation and ruim and veta the increment
	double rhok = 0.0;				//reason between the difference in the objective function and the difference of the approximation used
	double actual_change = 0.0;
	double predicted_change = 0.0;

	double tol_small = tol_small_1;
	int max_it = max_it_2;			//max iterations

	//Variables internal of use geral
	Matrix Hes(4, 4);
	Matrix Gra(4, 1);
	Matrix pGra(3, 1);

	Matrix dHes(3, 3);
	Matrix dGra(3, 1);

	Matrix P(3, 3);
	Matrix D(3, 3);
	Matrix cHes(3, 3);

	Matrix Hes2(4, 4);
	Matrix Gra2(4, 1);

	Matrix pk(3);
	Matrix pb(3);
	Matrix xk(4);

	if (write_report)
		fprintf(f_TR_report[seq_number], "FindSaddleSolutionDegenerated\n");

	//Initialization of the method - initial guess - obtained of the solution previous - problem of minimum distance
	for (int i = 0; i < 4; i++)
		xk(i, 0) = (*solution)(i, 0);

	//Error of the minimum distance
	GradientPhase1(xk, Gra);
	double error = norm(Gra);
	double ob = ObjectivePhase1(xk);

	bool saddle_basin = false;		//variable boolean that indicates that the bacia of attraction of the saddle was localizada
	int it = 1;

	//Hessian of the objective function
	HessianPhase1(xk, Hes);
	dHes = transp(*P_0)*Hes*(*P_0);
	//Calculating directions main and curvatures main of the objective function
	cHes = dHes;
	fulleigen1(cHes, P, D, tol_eig);

	//Stopping criterion
	double max_eig = -1e100;
	for (int i = 0; i < 3; i++)
	{
		if (D(i, i) > max_eig)
			max_eig = D(i, i);
	}
	double tol_ortho = tol_convective * abs(max_eig);

	if (CharacterizeCriticalPointDegenerated(&xk, P_0) == 3)
	{
		if (write_report)
			fprintf(f_TR_report[seq_number], "FindSaddleSolution initial guess is already in saddle basin\n");
		saddle_basin = true;//saddle basin
	}
	else
	{
		if (write_report)
			fprintf(f_TR_report[seq_number], "FindSaddleSolution initial guess is not in saddle basin\n");


		Surface* surf1 = db.surfaces[surf1_ID - 1];		//Pointer to the surface 1
		Surface* surf2 = db.surfaces[surf2_ID - 1];		//Pointer to the surface 2

		Matrix nA(3);
		if (surf1->degeneration[0] == false && surf1->degeneration[1] == false)
			surf1->NormalExt(&xk(0, 0), &xk(1, 0), &nA);

		Matrix nB(3);
		if (surf2->degeneration[0] == false && surf2->degeneration[1] == false)
			surf2->NormalExt(&xk(2, 0), &xk(3, 0), &nB);

		bool normal_fixed = true;
		int info = 0;

		//Degeneration local
		//Direction of search formada through the eigenvector associated to the smallest eigenvalue
		Matrix P_0local(3, 1);
		for (int i = 0; i < 3; i++)
			P_0local(i, 0) = P(i, 0);
		Matrix composed_P_0 = (*P_0)*P_0local;
		bool conv1 = FindMinimumGapDegenerated(c_data, &composed_P_0, &xk, info, normal_fixed, nA, nB);

		if (conv1 == false)
			return false;
		else
		{
			if (return_gap == true)
			{
				//Saves result in solution
				for (int i = 0; i < 4; i++)
					(*solution)(i, 0) = xk(i, 0);
				return true;
			}

			///////////////
			if (surf1->degeneration[0] == false && surf1->degeneration[1] == false)
				surf1->NormalExt(&xk(0, 0), &xk(1, 0), &nA);
			if (surf2->degeneration[0] == false && surf2->degeneration[1] == false)
				surf2->NormalExt(&xk(2, 0), &xk(3, 0), &nB);

			normal_fixed = true;
			info = 0;
			//Hessian of the objective function
			HessianPhase1(xk, Hes);
			dHes = transp(*P_0)*Hes*(*P_0);
			//Calculating directions main and curvatures main of the objective function
			cHes = dHes;
			fulleigen1(cHes, P, D, tol_eig);
			//Degeneration local
			//Direction of search formada through the eigenvector associated to the smallest eigenvalue
			for (int i = 0; i < 3; i++)
				P_0local(i, 0) = P(i, 0);
			composed_P_0 = (*P_0)*P_0local;
			conv1 = FindMinimumGapDegenerated(c_data, &composed_P_0, &xk, info, normal_fixed, nA, nB);

			//////////

			///////////////
			if (surf1->degeneration[0] == false && surf1->degeneration[1] == false)
				surf1->NormalExt(&xk(0, 0), &xk(1, 0), &nA);
			if (surf2->degeneration[0] == false && surf2->degeneration[1] == false)
				surf2->NormalExt(&xk(2, 0), &xk(3, 0), &nB);

			normal_fixed = true;
			info = 0;
			//Hessian of the objective function
			HessianPhase1(xk, Hes);
			dHes = transp(*P_0)*Hes*(*P_0);
			//Calculating directions main and curvatures main of the objective function
			cHes = dHes;
			fulleigen1(cHes, P, D, tol_eig);
			//Degeneration local
			//Direction of search formada through the eigenvector associated to the smallest eigenvalue
			for (int i = 0; i < 3; i++)
				P_0local(i, 0) = P(i, 0);
			composed_P_0 = (*P_0)*P_0local;
			conv1 = FindMinimumGapDegenerated(c_data, &composed_P_0, &xk, info, normal_fixed, nA, nB);

			//////////

			///////////////
			if (surf1->degeneration[0] == false && surf1->degeneration[1] == false)
				surf1->NormalExt(&xk(0, 0), &xk(1, 0), &nA);
			if (surf2->degeneration[0] == false && surf2->degeneration[1] == false)
				surf2->NormalExt(&xk(2, 0), &xk(3, 0), &nB);

			normal_fixed = true;
			info = 0;
			//Hessian of the objective function
			HessianPhase1(xk, Hes);
			dHes = transp(*P_0)*Hes*(*P_0);
			//Calculating directions main and curvatures main of the objective function
			cHes = dHes;
			fulleigen1(cHes, P, D, tol_eig);
			//Degeneration local
			//Direction of search formada through the eigenvector associated to the smallest eigenvalue
			
			for (int i = 0; i < 3; i++)
				P_0local(i, 0) = P(i, 0);
			composed_P_0 = (*P_0)*P_0local;
			conv1 = FindMinimumGapDegenerated(c_data, &composed_P_0, &xk, info, normal_fixed, nA, nB);

			//////////


		}
	}
	if (CharacterizeCriticalPointDegenerated(&xk, P_0) == 3)
	{
		if (write_report)
			fprintf(f_TR_report[seq_number], "After FindMinimumGap, initial guess is in saddle basin\n");
	}
	else
	{
		if (write_report)
			fprintf(f_TR_report[seq_number], "After FindMinimumGap, initial guess is not in saddle basin\n");
	}

	//Calculations of error and objective
	GradientPhase1(xk, Gra);
	//Gradient degenerated
	dGra = transp(*P_0)*Gra;
	//Hessian
	HessianPhase1(xk, Hes);
	//Hessian degenerated
	dHes = transp(*P_0)*Hes*(*P_0);
	//Error - forcing first input
	error = tol_ortho + 1.0;
	ob = ObjectivePhase1(xk);
	int flag_error = 0;
	if (write_report)
		fprintf(f_TR_report[seq_number], "FindSaddleSolutionDegenerated\n");
	char c = 'I';
	if (write_report)
		fprintf(f_TR_report[seq_number], "%d\t%.12e\t%.12e\t%.12e\t%.12e\t%.6e\t%.6e\t%.6e\t%c\t%.6e\t%.6e\t%.6e\n", it, xk(0, 0), xk(1, 0), xk(2, 0), xk(3, 0), ob, error, Deltak, c, rhok, actual_change, predicted_change);
	/////////////////////////////////Step of analysis of the bacia of attraction of the point of saddle/////////////////////////////////////
	while ((error > tol_ortho || norm(pk) > tol_convective) && it <= max_it)
	{
		//Calculating directions main and curvatures main of the objective function
		cHes = dHes;
		//Solves linear system - pb
		pb = fullsystem(cHes, -1.0*dGra, &flag_error);
		//Writing direction of search in the coordinates originais
		pk = (*P_0)*pb;
		double normpb = norm(pb);
		//Limita TR
		double norminc = norm(pk);
		if ((norminc + tol_convective) >= Deltak)
			pk = (Deltak / norminc)*pk;
		c = 'N';
		//////////////////////////UPDATING SOLUTION////////////////////////////////
		//Calculation of rhok
		double actual_change = ObjectivePhase1(xk) - ObjectivePhase1(xk + pk);
		double predicted_change = -(transp(Gra)*(pk)+0.5*transp(pk)*Hes*(pk))(0, 0);
		rhok = actual_change / predicted_change;

		if (abs(actual_change / ObjectivePhase1(xk)) < tol_ascent)
			rhok = 1.0;

		if (predicted_change < tol_small)
			rhok = 1.0;

		//update - case the model function represente bem the região. This avoids exit of the bacia of attraction of the saddle
		if (rhok >= 0.9 && rhok <= 1.1)
		{
			xk = xk + pk;
			Deltak = 2 * Deltak;
			if (Deltak > Deltamax)
				Deltak = Deltamax;
		}
		else
			Deltak = Deltak / 4;
		//increments iterations
		it++;
		//Calculations that will be used for check output of the loop - or in the next iteration, case not saia
		//Objective
		ob = ObjectivePhase1(xk);
		//Gradient
		GradientPhase1(xk, Gra);
		//Gradient degenerated
		dGra = transp(*P_0)*Gra;
		//Hessian
		HessianPhase1(xk, Hes);
		//Hessian degenerated
		dHes = transp(*P_0)*Hes*(*P_0);
		//Error - norm of the gradient
		error = norm(dGra);
		if (write_report)
			fprintf(f_TR_report[seq_number], "%d\t%.12e\t%.12e\t%.12e\t%.12e\t%.6e\t%.6e\t%.6e\t%c\t%.6e\t%.6e\t%.6e\n", it, xk(0, 0), xk(1, 0), xk(2, 0), xk(3, 0), ob, error, Deltak, c, rhok, actual_change, predicted_change);
	}
	if (write_report)
		fprintf(f_TR_report[seq_number], "Error\t%.6e\tTolerance\t%.6e\n", error, tol_ortho);


	////////////////////////////////Check final - if the solution and point of saddle/////////////////////////////////////
	//Return of the function
	if (error < tol_ortho)
	{
		//Saves result in solution
		for (int i = 0; i < 4; i++)
			(*solution)(i, 0) = xk(i, 0);

		//Checks the solution type encontrada
		int character = CharacterizeCriticalPointDegenerated(&xk, P_0);
		if (character == 3)
			return true;
		else
		{
			//Tests adicionais - converged but the saddle not was encontrada
			return false;
		}
	}
	else
		return false;
}

//Calculates the inversa of the Hessian
void SurfacePair::EvaluateInvertedHessian(SSContactData* c_data)
{
	for (int sol = 0; sol < c_data->n_solutions; sol++)
	{
		if (c_data->repeated[sol] == false && (c_data->return_value[sol] == 0 || c_data->return_value[sol] == 4))
		{
			//c_data->P_0[sol]->print();
			Matrix Hes(4, 4);
			Matrix xk(4);
			for (int i = 0; i < 4; i++)
				xk(i, 0) = c_data->convective[sol][i];
			HessianPhase1(xk, Hes);
			Matrix Hes_minor = transp(*c_data->P_0[sol])*Hes*(*c_data->P_0[sol]);
			int order_minor = c_data->P_0[sol]->getColumns();
			Matrix P(order_minor, order_minor);
			Matrix D(order_minor, order_minor);
			fulleigen1(Hes_minor, P, D, tol_eig);
			//D.print();
			//Inversão of the Hessian
			for (int i = 0; i < order_minor; i++)
			{
				D(i, i) = 1.0 / D(i, i);
			}
			Matrix invHes = (*c_data->P_0[sol])*P*D*transp(P)*transp(*c_data->P_0[sol]);
			for (int i = 0; i < 4; i++)
				for (int j = 0; j < 4; j++)
					c_data->invHessian[sol][i][j] = invHes(i, j);
		}
	}
}

void SurfacePair::WriteConvectiveRange()
{
	if (write_report)
	{
		//Writing the name of the contact pair
		char pair_name[100];
		sprintf(pair_name, "TR_report_surf_%d_surf_%d_%d", surf1_ID, surf2_ID, seq_number);
		strcpy(name, db.folder_name);	//folder of the job
		strcat(name, "TR/");			//directory TR
		GIRAFFE_MKDIR(name);					//creating directory TR
		strcat(name, pair_name);		//file name
		strcat(name, ".txt");			//creating file
		f_TR_report[seq_number] = fopen(name, "a");
		fprintf(f_TR_report[seq_number], "ConvectiveRange:\nMin\t%.6e\tMax\t%.6e\tRange\t%.6e\nMin\t%.6e\tMax\t%.6e\tRange\t%.6e\nMin\t%.6e\tMax\t%.6e\tRange\t%.6e\nMin\t%.6e\tMax\t%.6e\tRange\t%.6e\n",
			convective_min(0, 0), convective_max(0, 0), convective_range(0, 0),
			convective_min(1, 0), convective_max(1, 0), convective_range(1, 0),
			convective_min(2, 0), convective_max(2, 0), convective_range(2, 0),
			convective_min(3, 0), convective_max(3, 0), convective_range(3, 0));
		fclose(f_TR_report[seq_number]);
	}
}

//Performs automatic degeneration according to eigenvalues of Hessian
void SurfacePair::AutomaticDegenerationProcedure()
{
	//ALL
}