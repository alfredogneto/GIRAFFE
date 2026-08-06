#include "SplineElementPair.h"
#include "PlatformCompat.h"

#include "SPContactData.h"
#include "Surface.h"
#include "Spline.h"
#include "SplineElement.h"
#include "Database.h"
//Global variables
extern
Database db;


SplineElementPair::SplineElementPair()
{
}

SplineElementPair::~SplineElementPair()
{
}

//Values Default of tolerancias and other variables
void SplineElementPair::DefaultValues()
{
	//Tolerance for precision of convective coordinates
	tol_convective = 1e-12;
	//Tolerance for precision of a small number (machine precision related)
	tol_small_1 = 1e-3;
	//Tolerance for precision of eigenvalues extraction
	tol_eig = 1e-14;
	//Factor to be used to modify the TR Dogleg Path when ascending directions are found
	tol_ascent = 1e-4;
	//Maximum iterations for searching minimum points
	max_it_1 = 50;
	//Flag to write convergence report for LCP's
	//write_report = false;

	convective_range = Matrix(2);
	convective_max = Matrix(2);
	convective_min = Matrix(2);

	seq_number = 0;

	alloc_control = false;
}

//Allocates memory
void SplineElementPair::Alloc(SPContactData* c_data)
{
	if (alloc_control == false)
	{
		n_pointwise = c_data->n_solutions;
		cNR1 = new Matrix * [n_pointwise];
		f_TR_report = new FILE * [n_pointwise];
		f_DEG_report = new FILE * [n_pointwise];

		for (int i = 0; i < n_pointwise; i++)
		{
			cNR1[i] = new Matrix(2);
		}

		if (write_report)
		{
			for (int i = 0; i < n_pointwise; i++)
			{
				InitializeTRReport(i);
				InitializeDEGReport(i);
			}
		}
		alloc_control = true;
	}

}

//Frees memory
void SplineElementPair::Free()
{
	if (alloc_control == true)
	{
		//Limpeza of memory
		for (int i = 0; i < n_pointwise; i++)
		{
			delete cNR1[i];
		}
		delete[]cNR1;
		delete[]f_TR_report;
		delete[]f_DEG_report;

		alloc_control = false;
	}
}

void SplineElementPair::OpenTRReport(int index)
{
	//Writing the name of the contact pair
	char pair_name[100];
	sprintf(pair_name, "TR_report_Spline%d_Spline%d__ele_%d_ele_%d__%d", spline1_ID, spline2_ID, surf1_ID, surf2_ID, index);
	strcpy(name, db.folder_name);	//folder of the job
	strcat(name, "TR/");			//directory TR
	GIRAFFE_MKDIR(name);					//creating directory TR
	strcat(name, pair_name);		//file name
	strcat(name, ".txt");			//creating file
	f_TR_report[index] = fopen(name, "a");	//opens file
}

void SplineElementPair::OpenDEGReport(int index)
{
	//Writing the name of the contact pair
	char pair_name_deg[100];
	sprintf(pair_name_deg, "DEG_report_Spline%d_Spline%d__ele_%d_ele_%d__%d", spline1_ID, spline2_ID, surf1_ID, surf2_ID, index);
	strcpy(name_deg, db.folder_name);	//folder of the job
	strcat(name_deg, "DEG/");			//directory TR
	GIRAFFE_MKDIR(name_deg);					//creating directory TR
	strcat(name_deg, pair_name_deg);	//file name
	strcat(name_deg, ".txt");			//creating file
	f_DEG_report[index] = fopen(name_deg, "a");	//opens file
}

void SplineElementPair::InitializeTRReport(int index)
{
	//Writing the name of the contact pair
	char pair_name[100];
	sprintf(pair_name, "TR_report_Spline%d_Spline%d__ele_%d_ele_%d__%d", spline1_ID, spline2_ID, surf1_ID, surf2_ID, index);
	strcpy(name, db.folder_name);	//folder of the job
	strcat(name, "TR/");			//directory TR
	GIRAFFE_MKDIR(name);					//creating directory TR
	strcat(name, pair_name);		//file name
	strcat(name, ".txt");			//creating file
	f_TR_report[index] = fopen(name, "w");	//opens file

	fprintf(f_TR_report[index], "///////////////////////////////////////////////////////////////////////\n");
	fprintf(f_TR_report[index], "InitializeConvectiveRange\n");
	fprintf(f_TR_report[index], "Convective 1:\t%.6e\tto\t%.6e\tRange:\t%.6e\n", convective_min(0, 0), convective_max(0, 0), convective_range(0, 0));
	fprintf(f_TR_report[index], "Convective 2:\t%.6e\tto\t%.6e\tRange:\t%.6e\n", convective_min(1, 0), convective_max(1, 0), convective_range(1, 0));
	fclose(f_TR_report[index]);

	fclose(f_TR_report[index]);			//closes file	
}

void SplineElementPair::InitializeDEGReport(int index)
{
	//Writing the name of the contact pair
	char pair_name[100];
	sprintf(pair_name, "DEG_report_Spline%d_Spline%d__ele_%d_ele_%d__%d", spline1_ID, spline2_ID, surf1_ID, surf2_ID, index);
	strcpy(name_deg, db.folder_name);	//folder of the job
	strcat(name_deg, "DEG/");			//directory TR
	GIRAFFE_MKDIR(name_deg);					//creating directory TR
	strcat(name_deg, pair_name);		//file name
	strcat(name_deg, ".txt");			//creating file
	f_DEG_report[index] = fopen(name_deg, "w");	//opens file

	fprintf(f_DEG_report[index], "///////////////////////////////////////////////////////////////////////\n");
	fprintf(f_DEG_report[index], "InitializeConvectiveRange\n");
	fprintf(f_DEG_report[index], "Convective 1:\t%.6e\tto\t%.6e\tRange:\t%.6e\n", convective_min(0, 0), convective_max(0, 0), convective_range(0, 0));
	fprintf(f_DEG_report[index], "Convective 2:\t%.6e\tto\t%.6e\tRange:\t%.6e\n", convective_min(1, 0), convective_max(1, 0), convective_range(1, 0));

	fclose(f_DEG_report[index]);
	fclose(f_DEG_report[index]);			//closes file	
}

void SplineElementPair::PreCalc()
{
	//Initializes and sets the range of validity for the convective coordinates
	InitializeConvectiveRange();

	//Setting the minimum convective range
	minimum_convective_range = 1e100;
	for (int i = 0; i < 2; i++)
	{
		if (convective_range(i, 0) < minimum_convective_range && convective_range(i, 0) != 0.0)
			minimum_convective_range = convective_range(i, 0);
	}
}

void SplineElementPair::BeginStepCheck(SPContactData* c_data)
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
	SplineElement* surf1 = db.splines[spline1_ID - 1]->sp_element[surf1_ID];		//Pointer to the surface 1
	SplineElement* surf2 = db.splines[spline2_ID - 1]->sp_element[surf2_ID];		//Pointer to the surface 2
	Matrix GammaA(3);
	Matrix GammaB(3);
	Matrix cNR1_deg(2);
	Matrix cNR2_deg(2);
	double g_deg1;
	double g_deg2;
	//In the start of each step and tomada the decisão acerca of strong_candidate or not (with base in the check of bounding box in torno of the surfaces)
	bool converged1 = false;
	bool strong_candidate = false;
	bool deg1_converged = false;
	bool deg2_converged = false;
	//First check - bounding box overlap
	double inflation_factor = 4.0;
	strong_candidate = true;
	//strong_candidate = BoxOverlap(surf1->box, surf2->box, inflation_factor);

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

	//Saving ultimas convective coordinates convergidas in the current - evitando that cases that divergiram be usados as estimativas initial of the new increment
	for (int ip = 0; ip < n_pointwise; ip++)
	{
		c_data->convective[ip][0] = c_data->copy_convective[ip][0];
		c_data->convective[ip][1] = c_data->copy_convective[ip][1];

		c_data->repeated[ip] = false;			//Indicates that the solution [ip] always will be considerada
	}
	for (int ip = 0; ip < n_pointwise; ip++)
	{
		converged1 = false;
		int charact1 = 3;
		seq_number = ip;//for control of the report


		//Degeneration basis - canonical basis
		(*c_data->P[ip])(0, 0) = 1.0;
		(*c_data->P[ip])(1, 1) = 1.0;
		//Degenerated coordinates index (in the new basis)
		//Initially assumed as false
		c_data->degenerated[ip] = false;
		c_data->deg_control[ip][0] = false;
		c_data->deg_control[ip][1] = false;
		//Degenerated coordinates values - in case of degeneration, just fill the desired coordinate value
		c_data->copy_deg_coordinates[ip][0] = (surf1->knot_element[2] + surf1->knot_element[3]) / 2;
		c_data->copy_deg_coordinates[ip][1] = (surf2->knot_element[2] + surf2->knot_element[3]) / 2;

		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		//CASE 1: not era strong candidate (or first calculation)
		if (c_data->copy_return_value[ip] == 2)
		{
			if (write_report)
				fprintf(f_TR_report[ip], "Performing initial guess\n");
			InitialGuess(c_data); //Performs initial guess with criterion geometrico -> writes in 'convective'
			for (int i = 0; i < 2; i++)
				(*cNR1[ip])(i, 0) = c_data->convective[ip][i];
		}
		//CASE 2: or and copy_return_value 0 or 4 (already there is solution disponivel proveniente of a step previous)
		else
		{
			(*cNR1[ip])(0, 0) = c_data->copy_convective[ip][0];
			(*cNR1[ip])(1, 0) = c_data->copy_convective[ip][1];
		}

		c_data->return_value[ip] = 0;		//e strong (default)
		int info = 0;

		////////////NO DEGENERATION///////////////
		if (converged1 == false)
		{
			if (write_report)
				fprintf(f_TR_report[ip], "Search for minimum\n");
			//Determination of minimum or intersection
			converged1 = FindMinimumSolution(c_data, cNR1[ip], info); //Checks convergence of the modelo - true/false
			c_data->return_value[ip] = VerifyConvectiveRange(*cNR1[ip]); //Checks the range of the convective coordinate (2 longe, 4 next or 0 in the range) 
			charact1 = CharacterizeCriticalPoint(cNR1[ip]); //Caracteriza the point critico (0 minimum strict or 4 other type of problem)
			if (write_report)
			{
				fprintf(f_TR_report[ip], "Return value is %d \n", c_data->return_value[ip]);
				fprintf(f_TR_report[ip], "Converged %d \n", converged1);
				fprintf(f_TR_report[ip], "Critical point characterization %d \n", charact1);
			}
		}

		////////////AUTOMATIC DEGENERATION///////////////
		if (c_data->return_value[ip] == 2 || converged1 == false || charact1 == 4)
		{
			//Degeneration with respect to spline 1

			//Degenerated coordinates index (in the new basis)
			//Initially assumed as false
			c_data->degenerated[ip] = true;
			c_data->deg_control[ip][0] = true;
			c_data->deg_control[ip][1] = false;

			//Correction in case of degeneration
			Matrix temp_coordinates(2);
			for (int index = 0; index < 2; index++)
				temp_coordinates(index, 0) = (*cNR1[ip])(index, 0);
			//Transforming into degeneration basis
			temp_coordinates = transp(*c_data->P[ip]) * temp_coordinates;
			//Fixing degenerated coordinates
			for (int i = 0; i < 2; i++)
				if (c_data->deg_control[ip][i] == true)
					temp_coordinates(i, 0) = c_data->copy_deg_coordinates[ip][i];
			//Transforming into original basis
			temp_coordinates = (*c_data->P[ip]) * temp_coordinates;
			//Copying info into cNR1 vector
			for (int i = 0; i < 2; i++)
				(*cNR1[ip])(i, 0) = temp_coordinates(i, 0);

			if (write_report)
				fprintf(f_TR_report[ip], "Direct search for degenerated minimum for coordinate 1\n");

			//Degenerative operator
			c_data->MountDegenerativeOperator();

			converged1 = FindMinimumSolutionDegenerated(c_data, c_data->P_0[ip], cNR1[ip]); //Checks convergence of the modelo - true/false
			deg1_converged = converged1;

			surf1->SplinePoint((*cNR1[ip])(0, 0), GammaA);
			surf2->SplinePoint((*cNR1[ip])(1, 0), GammaB);
			g_deg1 = norm(GammaA - GammaB);
			cNR1_deg(0, 0) = (*cNR1[ip])(0, 0);
			cNR1_deg(1, 0) = (*cNR1[ip])(1, 0);

			//Degenration with respect to spline 2			

			//Degenerated coordinates index (in the new basis)
			//Initially assumed as false
			c_data->degenerated[ip] = true;
			c_data->deg_control[ip][0] = false;
			c_data->deg_control[ip][1] = true;

			//Correction in case of degeneration
			for (int index = 0; index < 2; index++)
				temp_coordinates(index, 0) = (*cNR1[ip])(index, 0);
			//Transforming into degeneration basis
			temp_coordinates = transp(*c_data->P[ip]) * temp_coordinates;
			//Fixing degenerated coordinates
			for (int i = 0; i < 2; i++)
				if (c_data->deg_control[ip][i] == true)
					temp_coordinates(i, 0) = c_data->copy_deg_coordinates[ip][i];
			//Transforming into original basis
			temp_coordinates = (*c_data->P[ip]) * temp_coordinates;
			//Copying info into cNR1 vector
			for (int i = 0; i < 2; i++)
				(*cNR1[ip])(i, 0) = temp_coordinates(i, 0);
			if (write_report)
				fprintf(f_TR_report[ip], "Direct search for degenerated minimum for coordinate 2\n");

			//Degenerative operator
			c_data->MountDegenerativeOperator();

			converged1 = FindMinimumSolutionDegenerated(c_data, c_data->P_0[ip], cNR1[ip]);
			deg2_converged = converged1;

			surf1->SplinePoint((*cNR1[ip])(0, 0), GammaA);
			surf2->SplinePoint((*cNR1[ip])(1, 0), GammaB);
			g_deg2 = norm(GammaA - GammaB);
			cNR2_deg(0, 0) = (*cNR1[ip])(0, 0);
			cNR2_deg(1, 0) = (*cNR1[ip])(1, 0);

			//Assuming solution with lower gap
			if (deg1_converged == true && deg2_converged == true)
				if (g_deg1 <= g_deg2) {
					//Initially assumed as false
					c_data->degenerated[ip] = true;
					c_data->deg_control[ip][0] = true;
					c_data->deg_control[ip][1] = false;

					//Correction in case of degeneration
			Matrix temp_coordinates(2);
			for (int index = 0; index < 2; index++)
				temp_coordinates(index, 0) = (*cNR1[ip])(index, 0);
			//Transforming into degeneration basis
			temp_coordinates = transp(*c_data->P[ip]) * temp_coordinates;
			//Fixing degenerated coordinates
			for (int i = 0; i < 2; i++)
				if (c_data->deg_control[ip][i] == true)
					temp_coordinates(i, 0) = c_data->copy_deg_coordinates[ip][i];
			//Transforming into original basis
			temp_coordinates = (*c_data->P[ip]) * temp_coordinates;
			//Copying info into cNR1 vector
			for (int i = 0; i < 2; i++)
				(*cNR1[ip])(i, 0) = temp_coordinates(i, 0);

					if (write_report)
						fprintf(f_TR_report[ip], "Direct search for degenerated minimum for coordinate 1\n");

					//Degenerative operator
					c_data->MountDegenerativeOperator();

					(*cNR1[ip])(0, 0) = cNR1_deg(0, 0);
					(*cNR1[ip])(1, 0) = cNR1_deg(1, 0);
				}
				else {
					//Initially assumed as false
					c_data->degenerated[ip] = true;
					c_data->deg_control[ip][0] = false;
					c_data->deg_control[ip][1] = true;

					//Correction in case of degeneration			
					for (int index = 0; index < 2; index++)
						temp_coordinates(index, 0) = (*cNR1[ip])(index, 0);
					//Transforming into degeneration basis
					temp_coordinates = transp(*c_data->P[ip]) * temp_coordinates;
					//Fixing degenerated coordinates
					for (int i = 0; i < 2; i++)
						if (c_data->deg_control[ip][i] == true)
							temp_coordinates(i, 0) = c_data->copy_deg_coordinates[ip][i];
					//Transforming into original basis
					temp_coordinates = (*c_data->P[ip]) * temp_coordinates;
					//Copying info into cNR1 vector
					for (int i = 0; i < 2; i++)
						(*cNR1[ip])(i, 0) = temp_coordinates(i, 0);
			if (write_report)
				fprintf(f_TR_report[ip], "Direct search for degenerated minimum for coordinate 2\n");

			//Degenerative operator
			c_data->MountDegenerativeOperator();

					(*cNR1[ip])(0, 0) = cNR2_deg(0, 0);
					(*cNR1[ip])(1, 0) = cNR2_deg(1, 0);
				}
			c_data->return_value[ip] = VerifyConvectiveRange(*cNR1[ip]);
			charact1 = CharacterizeCriticalPointDegenerated(cNR1[ip], c_data->P_0[ip], false);
		}
		//if (c_data->return_value[ip] == 2 || converged1 == false || charact1 == 4)
		//{
		//	//Degenerated coordinates index (in the new basis)
		//	//Initially assumed as false
		//	c_data->degenerated[ip] = true;
		//	c_data->deg_control[ip][0] = false;
		//	c_data->deg_control[ip][1] = true;

		//	//Correction in case of degeneration
		//	Matrix temp_coordinates(2);
		//	for (int index = 0; index < 2; index++)
		//		temp_coordinates(index, 0) = (*cNR1[ip])(index, 0);
		//	//Transforming into degeneration basis
		//	temp_coordinates = transp(*c_data->P[ip]) * temp_coordinates;
		//	//Fixing degenerated coordinates
		//	for (int i = 0; i < 2; i++)
		//		if (c_data->deg_control[ip][i] == true)
		//			temp_coordinates(i, 0) = c_data->copy_deg_coordinates[ip][i];
		//	//Transforming into original basis
		//	temp_coordinates = (*c_data->P[ip]) * temp_coordinates;
		//	//Copying info into cNR1 vector
		//	for (int i = 0; i < 2; i++)
		//		(*cNR1[ip])(i, 0) = temp_coordinates(i, 0);
		//	if (write_report)
		//		fprintf(f_TR_report[ip], "Direct search for degenerated minimum for coordinate 2\n");

		//	//Degenerative operator
		//	c_data->MountDegenerativeOperator();

		//	converged1 = FindMinimumSolutionDegenerated(c_data, c_data->P_0[ip], cNR1[ip]);
		//	c_data->return_value[ip] = VerifyConvectiveRange(*cNR1[ip]);
		//	//charact1 = CharacterizeCriticalPoint(cNR1[ip]);
		//	charact1 = CharacterizeCriticalPointDegenerated(cNR1[ip], c_data->P_0[ip], false);
		//}


		///////////////////////////////VERIFICATIONS//////////////////////////////////

		c_data->return_value[ip] = VerifyConvectiveRange(*cNR1[ip]);

		//Saves in the variables
		c_data->convective[ip][0] = (*cNR1[ip])(0, 0);
		c_data->convective[ip][1] = (*cNR1[ip])(1, 0);

		if (write_report)
			fprintf(f_TR_report[ip], "BeginStep return value is %d \n", c_data->return_value[ip]);

		if (VerifyConvectiveRange(*cNR1[ip]) == 2)
		{
			c_data->return_value[ip] = 2;	//not and strong
			if (write_report)
				fprintf(f_TR_report[ip], "Candidate is not Strong. Strict minimum far from range of interest.\n");
		}
		//if (charact1 != 0)
		//{
		//	c_data->return_value[ip] = 2;	//not and strong
		//	if (write_report)
		//		fprintf(f_TR_report[ip], "Candidate is not Strong. Critical point is not a minimum.\n");
		//}
		//if (converged1 == false) 
		//{
		//	c_data->return_value[ip] = 2;	//not and strong
		//	if (write_report)
		//		fprintf(f_TR_report[ip], "Candidate is not Strong. Solution diverged.\n");
		//}		
	}
	if (write_report)
	{
		for (int ip = 0; ip < n_pointwise; ip++)
			fclose(f_TR_report[ip]);
	}

}

//Determination via optimization of the convective coordinates of the pair of surfaces
void SplineElementPair::SolveLCP(SPContactData* c_data)
{
	///////////////////////////////Pointers and variables to make access easier//////////////////////////////////////////
	SplineElement* surf1 = db.splines[spline1_ID - 1]->sp_element[surf1_ID];		//Pointer to the surface 1
	SplineElement* surf2 = db.splines[spline2_ID - 1]->sp_element[surf2_ID];		//Pointer to the surface 2
	Matrix GammaA(3);
	Matrix GammaB(3);
	Matrix cNR1_deg(2);
	Matrix cNR2_deg(2);
	double g_deg1;
	double g_deg2;

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

	//Values saved through the function (indicativos of the that ocorreu in its execution)
	//0 - Converged and this in the range of interest for contact and not there was grande change of position of the point material of contact, in the case of haver contact previous - OK
	//1 - Not there was convergence (can or not to be in the range for contact) - and and strong_candidate. Return problematico responsavel by veto of step of the simulation.
	//2 - Not and strong_candidate - OK
	//3 - change abrupta in the position of contact. Return problematico responsavel by veto of step of the simulation.
	//4 - There was convergence, but not this in the range for contact. Possivelmente element vizinho. - OK

	//Sweep of the solutions active
	for (int i = 0; i < c_data->n_solutions; i++)
	{
		//Degeneration basis - canonical basis
		(*c_data->P[i])(0, 0) = 1.0;
		(*c_data->P[i])(1, 1) = 1.0;

		//Degenerated coordinates values - in case of degeneration, just fill the desired coordinate value
		c_data->copy_deg_coordinates[i][0] = (surf1->knot_element[2] + surf1->knot_element[3]) / 2;
		c_data->copy_deg_coordinates[i][1] = (surf2->knot_element[2] + surf2->knot_element[3]) / 2;


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

					//If before not era strong candidate, the copy_convective not traz info valida. Thus, calculates the initial guess (geometrico) novamente
					if (c_data->copy_return_value[i] == 2)
					{
						InitialGuess(c_data); //Performs initial guess with criterion geometrico -> writes in 'convective'
						for (int ii = 0; ii < 2; ii++)
							(*cNR1[i])(ii, 0) = c_data->convective[i][ii];
					}
				}
				//CASE 3: there was convergence
				//Takes as initial guess the last solution converged
				if (c_data->return_value[i] == 0 || c_data->return_value[i] == 4)
				{
					(*cNR1[i])(0, 0) = c_data->convective[i][0];
					(*cNR1[i])(1, 0) = c_data->convective[i][1];
				}

				bool converged1 = false;
				bool deg1_converged = false;
				bool deg2_converged = false;
				int charact1 = 3;
				int info = 0;

				if (c_data->degenerated[i] == false)
				{
					//////////////NO DEGENERATION///////////////
					if (converged1 == false)
					{
						if (write_report)
							fprintf(f_TR_report[i], "Search for minimum\n");

						//Degenerated coordinates index (in the new basis)
						c_data->deg_control[i][0] = false;
						c_data->deg_control[i][1] = false;

						//Degenerative operator
						c_data->MountDegenerativeOperator();

						//Determination of minimum or intersection
						converged1 = FindMinimumSolution(c_data, cNR1[i], info);
					}
				}
				else
				{
					////////////AUTOMATIC DEGENERATION///////////////
					if (converged1 == false && c_data->degenerated[i] == true)
					{
						//Degenerated coordinates index (in the new basis)
						////Initially assumed as false
						c_data->deg_control[i][0] = true;
						c_data->deg_control[i][1] = false;

						//Correction in case of degeneration
						Matrix temp_coordinates(2);
						for (int index = 0; index < 2; index++)
							temp_coordinates(index, 0) = (*cNR1[i])(index, 0);
						//Transforming into degeneration basis
						temp_coordinates = transp(*c_data->P[i]) * temp_coordinates;
						//Fixing degenerated coordinates
						for (int index = 0; index < 2; index++)
							if (c_data->deg_control[i][index] == true)
								temp_coordinates(index, 0) = c_data->copy_deg_coordinates[i][index];
						//Transforming into original basis
						temp_coordinates = (*c_data->P[i]) * temp_coordinates;
						//Copying info into cNR1 vector
						for (int index = 0; index < 2; index++)
							(*cNR1[i])(index, 0) = temp_coordinates(index, 0);

						if (write_report)
							fprintf(f_TR_report[i], "Direct search for degenerated minimum for coordinate 1\n");

						//Degenerative operator
						c_data->MountDegenerativeOperator();

						converged1 = FindMinimumSolutionDegenerated(c_data, c_data->P_0[i], cNR1[i]);

						surf1->SplinePoint((*cNR1[i])(0, 0), GammaA);
						surf2->SplinePoint((*cNR1[i])(1, 0), GammaB);
						g_deg1 = norm(GammaA - GammaB);
						cNR1_deg(0, 0) = (*cNR1[i])(0, 0);
						cNR1_deg(1, 0) = (*cNR1[i])(1, 0);

						//Degenerated coordinates index (in the new basis)
						////Initially assumed as false
						c_data->deg_control[i][0] = false;
						c_data->deg_control[i][1] = true;

						//Correction in case of degeneration						
						for (int index = 0; index < 2; index++)
							temp_coordinates(index, 0) = (*cNR1[i])(index, 0);
						//Transforming into degeneration basis
						temp_coordinates = transp(*c_data->P[i]) * temp_coordinates;
						//Fixing degenerated coordinates
						for (int index = 0; index < 2; index++)
							if (c_data->deg_control[i][index] == true)
								temp_coordinates(index, 0) = c_data->copy_deg_coordinates[i][index];
						//Transforming into original basis
						temp_coordinates = (*c_data->P[i]) * temp_coordinates;
						//Copying info into cNR1 vector
						for (int index = 0; index < 2; index++)
							(*cNR1[i])(index, 0) = temp_coordinates(index, 0);

						if (write_report)
							fprintf(f_TR_report[i], "Direct search for degenerated minimum for coordinate 2\n");

						//Degenerative operator
						c_data->MountDegenerativeOperator();

						converged1 = FindMinimumSolutionDegenerated(c_data, c_data->P_0[i], cNR1[i]);

						surf1->SplinePoint((*cNR1[i])(0, 0), GammaA);
						surf2->SplinePoint((*cNR1[i])(1, 0), GammaB);
						g_deg2 = norm(GammaA - GammaB);
						cNR2_deg(0, 0) = (*cNR1[i])(0, 0);
						cNR2_deg(1, 0) = (*cNR1[i])(1, 0);

						if (deg1_converged == true && deg2_converged == true) {							
							if (g_deg1 <= g_deg2) {
								//Degenerated coordinates index (in the new basis)
								////Initially assumed as false
								c_data->deg_control[i][0] = true;
								c_data->deg_control[i][1] = false;

								//Correction in case of degeneration
								Matrix temp_coordinates(2);
								for (int index = 0; index < 2; index++)
									temp_coordinates(index, 0) = (*cNR1[i])(index, 0);
								//Transforming into degeneration basis
								temp_coordinates = transp(*c_data->P[i]) * temp_coordinates;
								//Fixing degenerated coordinates
								for (int index = 0; index < 2; index++)
									if (c_data->deg_control[i][index] == true)
										temp_coordinates(index, 0) = c_data->copy_deg_coordinates[i][index];
								//Transforming into original basis
								temp_coordinates = (*c_data->P[i]) * temp_coordinates;
								//Copying info into cNR1 vector
								for (int index = 0; index < 2; index++)
									(*cNR1[i])(index, 0) = temp_coordinates(index, 0);

								if (write_report)
									fprintf(f_TR_report[i], "Direct search for degenerated minimum for coordinate 1\n");

								//Degenerative operator
								c_data->MountDegenerativeOperator();
								(*cNR1[i])(0, 0) = cNR1_deg(0, 0);
								(*cNR1[i])(1, 0) = cNR1_deg(1, 0);
					}
							else {
						//Degenerated coordinates index (in the new basis)
						////Initially assumed as false
						c_data->deg_control[i][0] = false;
						c_data->deg_control[i][1] = true;

						//Correction in case of degeneration
						for (int index = 0; index < 2; index++)
							temp_coordinates(index, 0) = (*cNR1[i])(index, 0);
						//Transforming into degeneration basis
						temp_coordinates = transp(*c_data->P[i]) * temp_coordinates;
						//Fixing degenerated coordinates
						for (int index = 0; index < 2; index++)
							if (c_data->deg_control[i][index] == true)
								temp_coordinates(index, 0) = c_data->copy_deg_coordinates[i][index];
						//Transforming into original basis
						temp_coordinates = (*c_data->P[i]) * temp_coordinates;
						//Copying info into cNR1 vector
						for (int index = 0; index < 2; index++)
							(*cNR1[i])(index, 0) = temp_coordinates(index, 0);

						if (write_report)
							fprintf(f_TR_report[i], "Direct search for degenerated minimum for coordinate 2\n");

						//Degenerative operator
						c_data->MountDegenerativeOperator();
								(*cNR1[i])(0, 0) = cNR2_deg(0, 0);
								(*cNR1[i])(1, 0) = cNR2_deg(1, 0);
							}
						}

					}

				}

				if (converged1 == false)
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

				if (converged1 == true) {
					//Saves in the variables
					c_data->convective[i][0] = (*cNR1[i])(0, 0);
					c_data->convective[i][1] = (*cNR1[i])(1, 0);

					//Fills return value with 0, 2 or 4
					c_data->return_value[i] = VerifyConvectiveRange(*cNR1[i]);
					//Correction (for not creating a not-strong candidate here!):
					if (c_data->return_value[i] == 2)
						c_data->return_value[i] = 4;
					if (write_report)
						fprintf(f_TR_report[i], "SolveLCP return value is %d \n", c_data->return_value[i]);

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

int SplineElementPair::CharacterizeCriticalPoint(Matrix* solution)
{
	///////////////////////////////Pointers and variables to make access easier//////////////////////////////////////////
	SplineElement* surf1 = db.splines[spline1_ID - 1]->sp_element[surf1_ID];		//Pointer to the surface 1
	SplineElement* surf2 = db.splines[spline2_ID - 1]->sp_element[surf2_ID];		//Pointer to the surface 2

	//0 - minimum strict
	//1 - minimum not strict (intersection)
	//2 - transition minimum-saddle (just-touch)
	//3 - saddle 2 negative eigenvalues
	//4 - other
	Matrix xk(2);
	Matrix Gra(2);
	Matrix Hes(2, 2);
	Matrix P(2, 2);
	Matrix D(2, 2);
	for (int i = 0; i < 2; i++)
		xk(i, 0) = (*solution)(i, 0);
	double ob = ObjectivePhase1(xk);
	//Hessian of the objective function
	HessianPhase1(xk, Hes);
	//Calculating directions main and curvatures main of the objective function
	fulleigen1(Hes, P, D, tol_eig);

	double max_eig = -1e100;
	for (int i = 0; i < 2; i++)
	{
		if (D(i, i) > max_eig)
			max_eig = D(i, i);
	}
	double tol_intersect = max_eig * tol_convective * tol_convective;
	double tol_small = max_eig * tol_small_1;

	////minimum not strict (intersection)
	//if (ob < tol_intersect)
	//{
	//	if (write_report)
	//		fprintf(f_TR_report[seq_number], "Intersection found (tolerance %.6e). Eigenvalues are %.6e\t%.6e\n", tol_intersect, D(0, 0), D(1, 1));
	//	return 1;
	//}
	//minimum strict
	if (D(0, 0) >= tol_small && D(1, 1) >= tol_small)
	{
		if (write_report)
			fprintf(f_TR_report[seq_number], "Strict minimum found. Eigenvalues are %.6e\t%.6e\n", D(0, 0), D(1, 1));
		return 0;
	}
	if (write_report)
		fprintf(f_TR_report[seq_number], "Other critical point found. Eigenvalues are %.6e\t%.6e. Tolerance is %.6e.\n", D(0, 0), D(1, 1), tol_small);
	return 4;
}

int SplineElementPair::CharacterizeCriticalPointDegenerated(Matrix* solution, Matrix* P_0, bool print)
{
	///////////////////////////////Pointers and variables to make access easier//////////////////////////////////////////
	SplineElement* surf1 = db.splines[spline1_ID - 1]->sp_element[surf1_ID];		//Pointer to the surface 1
	SplineElement* surf2 = db.splines[spline2_ID - 1]->sp_element[surf2_ID];		//Pointer to the surface 2

	//0 - minimum strict
	//1 - minimum not strict (intersection)
	//2 - transition minimum-saddle (just-touch)
	//3 - saddle 2 negative eigenvalues
	//4 - other
	Matrix xk(2);
	Matrix Hes(2, 2);
	Matrix dHes(1, 1);
	Matrix P(1, 1);
	Matrix D(1, 1);
	for (int i = 0; i < 2; i++)
		xk(i, 0) = (*solution)(i, 0);
	double ob = ObjectivePhase1(xk);
	//Hessian of the objective function
	HessianPhase1(xk, Hes);
	dHes = transp(*P_0) * Hes * (*P_0);
	fulleigen1(dHes, P, D, tol_eig);

	double max_eig = -1e100;
	for (int i = 0; i < 1; i++)
	{
		if (D(i, i) > max_eig)
			max_eig = D(i, i);
	}
	double tol_intersect = max_eig * tol_convective * tol_convective;
	double tol_small = max_eig * tol_small_1;

	////minimum not strict (intersection)
	//if (ob < tol_intersect)
	//{
	//	if (write_report)
	//		fprintf(f_TR_report[seq_number], "Intersection found. Eigenvalues are %.6e\n", D(0, 0));
	//	return 1;
	//}
	//minimum strict
	if (D(0, 0) >= tol_small)
	{
		if (write_report)
			fprintf(f_TR_report[seq_number], "Strict minimum found. Eigenvalue is %.6e. Tolerance is %.6e.\n", D(0, 0), tol_small);
		return 0;
	}

	if (write_report)
		fprintf(f_TR_report[seq_number], "Other critical point found. Eigenvalues are %.6e\n", D(0, 0));
	return 4;
}

bool SplineElementPair::EndStepCheck(SPContactData* c_data)
{

	if (write_report)
	{
		for (int i = 0; i < c_data->n_solutions; i++)
		{
			OpenTRReport(i);
			fprintf(f_TR_report[i], "///////////////////////////////////////////////////////////////////////\n");
			fprintf(f_TR_report[i], "\nTime\t%.6f\tIteration\t%d\t", db.last_converged_time + db.current_time_step, db.current_iteration_number);
			fprintf(f_TR_report[i], "EndStepCheck\n");

			OpenDEGReport(i);
			fprintf(f_DEG_report[i], "\nTime\t%.6f\tIteration\t%d\t%d\t%d\t%d\t", db.last_converged_time + db.current_time_step, db.current_iteration_number, c_data->degenerated[i], c_data->deg_control[i][0], c_data->deg_control[i][1]);
		}
	}
	///////////////////////////////Pointers and variables to make access easier//////////////////////////////////////////
	SplineElement* surf1 = db.splines[spline1_ID - 1]->sp_element[surf1_ID];		//Pointer to the surface 1
	SplineElement* surf2 = db.splines[spline2_ID - 1]->sp_element[surf2_ID];		//Pointer to the surface 2

	//Return with problems: true
	//Return without problems: false
	bool c_problem = false;

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
				db.myprintf("LCP between surfaces %d and %d has presented problems. Code 1.\n", surf1_ID, surf2_ID);
				if (write_report)
				{
					fprintf(f_TR_report[i], "LCP between surfaces %d and %d has presented problems. Code 1.\n", surf1_ID, surf2_ID);
				}
				c_problem = true;
				//return true;
			}

			//Criterion numérico for indicate inversão of the normal vector of contact
			//double pen_tol = 1e-14;
			//if (c_data->g_n[i] < -pen_tol)

			//if (dot(*c_data->g[i], *c_data->g[i]) < tol_small_1)
			//{
			//	db.myprintf("LCP between surfaces %d and %d has presented problems. Code 5.\n", surf1_ID, surf2_ID);
			//	if (write_report)
			//	{
			//		fprintf(f_TR_report[i], "LCP between surfaces %d and %d has presented problems. Code 5.\n", surf1_ID, surf2_ID);
			//	}
			//	c_problem = true;
			//	//return true;
			//}

			if (dot(*c_data->n[i], *c_data->copy_n[i]) < -0.5)
			{
				db.myprintf("LCP between surfaces %d and %d has presented problems. Code 4.\n", surf1_ID, surf2_ID);
				if (write_report)
				{
					fprintf(f_TR_report[i], "LCP between surfaces %d and %d has presented problems. Code 4.\n", surf1_ID, surf2_ID);
				}
				c_problem = true;
				//return true;
			}

			//If there was muita change of range of the convective coordinates and indicativo of perda of the solution desired for c_bar
			if (c_data->return_value[i] == 3)
			{
				db.myprintf("LCP between surfaces %d and %d has presented problems. Code 3.\n", surf1_ID, surf2_ID);
				if (write_report)
				{
					fprintf(f_TR_report[i], "LCP between surfaces %d and %d has presented problems. Code 3.\n", surf1_ID, surf2_ID);
				}
				c_problem = true;
				//return true;
			}
		}
	}

	if (write_report)
	{
		for (int ip = 0; ip < c_data->n_solutions; ip++)
			fclose(f_TR_report[ip]);
		for (int ip = 0; ip < c_data->n_solutions; ip++)
			fclose(f_DEG_report[ip]);
	}

	return c_problem;
}

bool SplineElementPair::FindMinimumSolution(SPContactData* c_data, Matrix* solution, int& return_info)
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
	Matrix Hes(2, 2);
	Matrix Gra(2, 1);
	Matrix pGra(2, 1);
	Matrix pk(2);
	Matrix pb(2);
	Matrix pc(2);
	Matrix xk(2);
	Matrix P(2, 2);
	Matrix D(2, 2);
	Matrix cHes(2, 2);

	//Initialization of the method - initial guess
	for (int i = 0; i < 2; i++)
		xk(i, 0) = (*solution)(i, 0);

	//Stopping criterion
	HessianPhase1(xk, Hes);
	fulleigen1(Hes, P, D, tol_eig);
	double max_eig = -1e100;
	for (int i = 0; i < 2; i++)
	{
		if (D(i, i) > max_eig)
			max_eig = D(i, i);
	}
	double tol_ortho = tol_convective * abs(max_eig);
	double tol_small = tol_small_1;

	//Criterion for identify a intersection
	double tol_intersect = max_eig * tol_convective * tol_convective;

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
		fprintf(f_TR_report[seq_number], "%d\t%.12e\t%.12e\t%.6e\t%.6e\t%.6e\t%c\t%.6e\t%.6e\t%.6e\n", it, xk(0, 0), xk(1, 0), ob, error, Deltak, c, rhok, actual_reduction, predicted_reduction);
	/////////////////////////////////////////////////////BEGIN///////////////////////////////////////////////

	while ((error > tol_ortho || (norm(pk) > tol_convective && ob > tol_intersect)) && it <= max_it)
	{
		//Determination of the point of Cauchy
		double gragra = (transp(Gra) * Gra)(0, 0);
		double grahesgra = (transp(Gra) * Hes * Gra)(0, 0);
		double normgra = norm(Gra);
		Matrix pc;	//direction of the Cauchy point
		/////////////////Point of Cauchy///////////////////////
		if (grahesgra <= 0.0)
			pc = -(Deltak / normgra) * Gra;
		else
			pc = -gragra / grahesgra * Gra;
		double normpc = norm(pc);
		//Cauchy point outside the TR - use a fraction of it
		if ((normpc + tol_convective) >= Deltak)
		{
			pk = (Deltak / normpc) * pc; //Steep descent
			c = 'C';
		}
		else//Cauchy point inside the TR
		{
			//Calculating directions main and curvatures main of the objective function
			cHes = Hes;
			fulleigen1(cHes, P, D, tol_eig);
			//Writing gradient in the directions main
			pGra = transp(P) * Gra;
			//Determination of the smallest eigenvalue (min_eig)
			double min_eig = 1e100;
			for (int i = 0; i < 2; i++)
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
				for (int i = 0; i < 2; i++)
					pb(i, 0) = -pGra(i, 0) / (D(i, i) - (min_eig - abs(min_eig) * tol_ascent));
			}
			//If the smallest eigenvalue and larger that zero (tol_small) - direction of NR and chosen
			else
			{
				for (int i = 0; i < 2; i++)
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
				a1 = norm(pb - pc) * norm(pb - pc);
				b1 = 2 * (transp(pc) * (pb - pc))(0, 0);
				c1 = normpc * normpc - Deltak * Deltak;
				thetak = (-b1 + sqrt(b1 * b1 - 4 * a1 * c1)) / (2 * a1);
				//Determination of the path
				pk = pc + thetak * (pb - pc);
				c = 'D';
			}
		}

		//////////////////////////UPDATING SOLUTION////////////////////////////////
		//Calculation of rhok
		actual_reduction = ObjectivePhase1(xk) - ObjectivePhase1(xk + pk);
		predicted_reduction = -(transp(Gra) * pk + 0.5 * transp(pk) * Hes * pk)(0, 0);
		rhok = actual_reduction / predicted_reduction;

		if (abs(actual_reduction / ObjectivePhase1(xk)) < tol_ascent)
			rhok = 1.0;

		if (abs(predicted_reduction) < tol_small || abs(actual_reduction) < tol_small)
			rhok = 1.0;

		if (rhok >= 0)
		{
			//low reduction
			if (rhok < 0.25)
				Deltak = 0.25 * norm(pk);
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
			fprintf(f_TR_report[seq_number], "%d\t%.12e\t%.12e\t%.6e\t%.6e\t%.6e\t%c\t%.6e\t%.6e\t%.6e\n", it, xk(0, 0), xk(1, 0), ob, error, Deltak, c, rhok, actual_reduction, predicted_reduction);
	}
	if (write_report)
		fprintf(f_TR_report[seq_number], "Error\t%.6e\tTolerance\t%.6e\n", error, tol_ortho);

	//Return of the function
	if (error < tol_ortho)
	{
		for (int i = 0; i < 2; i++)
			(*solution)(i, 0) = xk(i, 0);
		return true;
	}
	else
		return false;
}

bool SplineElementPair::FindMinimumSolutionDegenerated(SPContactData* c_data, Matrix* P_0, Matrix* solution)
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

	Matrix Hes(2, 2);
	Matrix Gra(2, 1);
	Matrix xk(2);
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
	for (int i = 0; i < 2; i++)
		xk(i, 0) = (*solution)(i, 0);

	//Stopping criterion
	HessianPhase1(xk, Hes);
	deg_Hes = transp(*P_0) * Hes * (*P_0);
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
	deg_Gra = transp(*P_0) * Gra;
	deg_Hes = transp(*P_0) * Hes * (*P_0);
	//Error - forcing first input
	double error = tol_ortho + 1.0;
	//Initial guess report
	if (write_report)
		fprintf(f_TR_report[seq_number], "%d\t%.12e\t%.12e\t%.12e\t%.12e\t%.6e\t%.6e\t%.6e\t%c\t%.6e\t%.6e\t%.6e\n", it, xk(0, 0), xk(1, 0), ob, error, Deltak, c, rhok, actual_reduction, predicted_reduction);
	/////////////////////////////////////////////////////BEGIN///////////////////////////////////////////////
	while ((error > tol_ortho || norm((*P_0) * pk) > tol_convective) && it <= max_it)
	{
		//Determination of the point of Cauchy
		double gragra = (transp(deg_Gra) * deg_Gra)(0, 0);
		double grahesgra = (transp(deg_Gra) * deg_Hes * deg_Gra)(0, 0);
		double normgra = norm(deg_Gra);
		Matrix pc;	//direction of the Cauchy point
		/////////////////Point of Cauchy///////////////////////
		if (grahesgra <= 0.0)
			pc = -(Deltak / normgra) * deg_Gra;
		else
			pc = -gragra / grahesgra * deg_Gra;
		double normpc = norm(pc);
		//Cauchy point outside the TR - use a fraction of it
		if ((normpc + tol_convective) >= Deltak)
		{
			pk = (Deltak / normpc) * pc; //Steep descent
			c = 'C';
		}

		else//Cauchy point inside the TR
		{
			//Calculating directions main and curvatures main of the objective function
			cHes = deg_Hes;
			fulleigen1(cHes, P, D, tol_eig);
			//Writing gradient in the directions main
			pGra = transp(P) * deg_Gra;
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
					pb(i, 0) = -pGra(i, 0) / (D(i, i) - (min_eig - abs(min_eig) * tol_ascent));
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
				a = norm(pb - pc) * norm(pb - pc);
				b = 2 * (transp(pc) * (pb - pc))(0, 0);
				c = normpc * normpc - Deltak * Deltak;
				thetak = (-b + sqrt(b * b - 4 * a * c)) / (2 * a);
				//Determination of the path
				pk = pc + thetak * (pb - pc);
				c = 'D';
			}
		}

		//////////////////////////UPDATING SOLUTION////////////////////////////////
		//Calculation of rhok
		double actual_reduction = ObjectivePhase1(xk) - ObjectivePhase1(xk + (*P_0) * pk);
		double predicted_reduction = -(transp(deg_Gra) * pk + 0.5 * transp(pk) * deg_Hes * pk)(0, 0);
		rhok = actual_reduction / predicted_reduction;

		if (abs(actual_reduction / ObjectivePhase1(xk)) < tol_ascent)
			rhok = 1.0;

		if (abs(predicted_reduction) < tol_small || abs(actual_reduction) < tol_small)
			rhok = 1.0;

		if (rhok >= 0.0)
		{
			if (rhok < 0.25)//low reduction or even augmenting the objective function
				Deltak = 0.25 * norm(pk);//reduce TR
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
				xk = xk + (*P_0) * pk;
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
		deg_Gra = transp(*P_0) * Gra;
		deg_Hes = transp(*P_0) * Hes * (*P_0);
		//Error - norm of the gradient
		error = norm(deg_Gra);
		//Initial guess report
		if (write_report)
			fprintf(f_TR_report[seq_number], "%d\t%.12e\t%.12e\t%.6e\t%.6e\t%.6e\t%c\t%.6e\t%.6e\t%.6e\n", it, xk(0, 0), xk(1, 0), ob, error, Deltak, c, rhok, actual_reduction, predicted_reduction);
	}
	if (write_report)
		fprintf(f_TR_report[seq_number], "Error\t%.6e\tTolerance\t%.6e\n", error, tol_ortho);

	//Return of the function
	if (error < tol_ortho)
	{
		//Saves result in solution
		for (int i = 0; i < 2; i++)
			(*solution)(i, 0) = xk(i, 0);
		return true;
	}
	else
		return false;
}



//Calculates the inversa of the Hessian
void SplineElementPair::EvaluateInvertedHessian(SPContactData* c_data)
{
	for (int sol = 0; sol < c_data->n_solutions; sol++)
	{
		if (c_data->repeated[sol] == false && (c_data->return_value[sol] == 0 || c_data->return_value[sol] == 4))
		{
			//c_data->P_0[sol]->print();
			Matrix Hes(2, 2);
			Matrix xk(2);
			for (int i = 0; i < 2; i++)
				xk(i, 0) = c_data->convective[sol][i];
			HessianPhase1(xk, Hes);
			Matrix Hes_minor = transp(*c_data->P_0[sol]) * Hes * (*c_data->P_0[sol]);
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
			Matrix invHes = (*c_data->P_0[sol]) * P * D * transp(P) * transp(*c_data->P_0[sol]);
			for (int i = 0; i < 2; i++)
				for (int j = 0; j < 2; j++)
					c_data->invHessian[sol][i][j] = invHes(i, j);
		}
	}
}

void SplineElementPair::WriteConvectiveRange()
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
			convective_min(1, 0), convective_max(1, 0), convective_range(1, 0)/*,
		convective_min(2, 0), convective_max(2, 0), convective_range(2, 0),
			convective_min(3, 0), convective_max(3, 0), convective_range(3, 0)*/);
		fclose(f_TR_report[seq_number]);
	}
}