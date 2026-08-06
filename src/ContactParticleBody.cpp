#include "ContactParticleBody.h"
#include "PlatformCompat.h"

//#include "CollisionDetection.h"
#include "SSContactData.h"
#include "Matrix.h"
#include "Encoding.h"
#include "Particle.h"
#include "BodyGeometry.h"
#include "BoundingCylinder.h"
#include "BoundingTriangularBox.h"
#include "BoundingSphere.h"
#include "Boundary.h"
#include "Geometry.h"

#include "Database.h"
//Global variables
extern
Database db;

ContactParticleBody::ContactParticleBody()
{
}

ContactParticleBody::~ContactParticleBody()
{
}

void ContactParticleBody::Clear()
{
	
}

void ContactParticleBody::DefaultValues()
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
}

void ContactParticleBody::SaveConfiguration()
{
	//Cleaning contact pair
	prev_eligible = eligible;

	/*if (eligible == false)
	{
		Free();
	}*/
	//Zeroes the gap tangential if contact not and eligible
	if (eligible == false)
	{
		zeros(cd->g_t[0]);
		zeros(cd->copy_g_t[0]);
	}
	else
	{
		*cd->copy_g_t[0] = *cd->g_t[0];
	}

	cd->copy_g_n[0] = cd->g_n[0];
	*cd->copy_g[0] = *cd->g[0];
	*cd->copy_n[0] = *cd->n[0];
	cd->copy_return_value[0] = cd->return_value[0];
	cd->copy_characterization_index[0] = cd->characterization_index[0];
	cd->copy_convective[0][0] = cd->convective[0][0];
	cd->copy_convective[0][1] = cd->convective[0][1];
	cd->copy_convective[0][2] = cd->convective[0][2];
	cd->copy_convective[0][3] = cd->convective[0][3];
	cd->copy_degenerated[0] = cd->degenerated[0];
	cd->copy_stick[0] = cd->stick[0];
	cd->copy_convergedLCP[0] = cd->convergedLCP[0];
}

void ContactParticleBody::MountContact()
{
	if (eligible)
	{
		Alloc();
		//SetVariables();				//Sets variables for next evaluations
		EvaluateInvertedHessian();
		MountLocalContributions();		//Local contact contributions
	}
}

void ContactParticleBody::MountContactExplicit(double t)
{
	//ALL-Explicit
}

bool ContactParticleBody::NightOwlContact()
{
	//db.myprintf("\nNight owl contact test between body surfaces %d and %d\n", db.body_geometries[index1]->ptr_geom[sub_index1]->number, db.body_geometries[index2]->ptr_geom[sub_index2]->number);
	//Solves previous configuration - as a candidate do Night Owl contact, surely no previous evaluation of LCP and gap have been processed
	SolvePreviousContact();
	ProcessSurfacePair();
	//If contact is detected (eligible pair)
	if (eligible == true)
	{
		db.myprintf("\nBodyGeometry %d and BodyGeometry %d, %d and %d\n", index1, index2, sub_index1, sub_index2);
		db.myprintf("Gap value %.6e\n", cd[0].g_n[0]);
		db.myprintf("Prev Gap value %.6e\n", cd[0].copy_g_n[0]);

		return true;
	}
	else
	{
		//Error reported
		if (HaveErrors())
			return true;
	}
	return false;
}

void ContactParticleBody::EvaluateNormalGap()
{
	//Calculation of the function gap (scalar)
	SurfacePoints();
	//Gap cross
	*cd->g[0] = *GammaA - *GammaB;
	//Normal of the contact
	if (norm(*cd->g[0]) != 0.0)
		*cd->n[0] = (1.0 / norm(*cd->g[0]))*(*cd->g[0]);
	else
		zeros(cd->n[0]);
	cd->g_n[0] = +1.0*norm(*cd->g[0]);
	if (interface_1_flag && cd->return_value[0] == 0)
	{
		if (cd->g_n[0] <= *gnb)
			eligible = true;
	}
	if (interface_0_flag)
	{
		////Gap scalar - attempt of calculation with base in the normal previous (if the same exist)
		//bool previouscontact = false;
		//bool contact = false;
		//if (cd->copy_g_n[0] <= 0)
		//	previouscontact = true;
		//if (dot(*cd->n[0], *cd->copy_n[0]) >= 0)
		//	contact = previouscontact;
		//else
		//{
		//	if (previouscontact)
		//		contact = false;
		//	else
		//		contact = true;
		//}
		//if (contact)
		//	cd->g_n[0] = -1.0*norm(*cd->g[0]);
		//else
		//	cd->g_n[0] = +1.0*norm(*cd->g[0]);

		//if (cd->g_n[0] < 0.0 && cd->return_value[0] == 0)
		//	eligible = true;
	}
}

void ContactParticleBody::SolveLCP()
{
	if (write_report)
	{
		OpenTRReport();
		fprintf(f_TR_report, "///////////////////////////////////////////////////////////////////////\n");
		fprintf(f_TR_report, "\nTime\t%.6f\tIteration\t%d\t", db.last_converged_time + db.current_time_step, db.current_iteration_number);
		fprintf(f_TR_report, "SolveLCP\n");
	}

	//Return value:
	//0 - in range
	//4 - out of range
	//1 - problematic return

	bool converged1 = false;
	bool converged2 = false;
	int info = 0;
	Matrix cb(4);

	//Initial guess - in the case of divergence previous or first calculation
	//Tambem enters if in the step previous the solution encontrada era a intersection outside of the range - avoids solutions "viciadas" by anteriores inadequadas.
	if (cd->return_value[0] == 1 || cd->return_value[0] == 2 || (cd->copy_characterization_index[0] == 1 && cd->return_value[0] == 4) || cd->return_value[0] == 3 || cd->copy_convergedLCP[0]==false || cd->convergedLCP[0]==false)
	{
		InitialGuess();
		cb(0, 0) = cd->convective[0][0];
		cb(1, 0) = cd->convective[0][1];
		cb(2, 0) = cd->convective[0][2];
		cb(3, 0) = cd->convective[0][3];
	}
	//If not, takes solution previous in the convective. Note: in the case of first input after sucesso of convergence in step previous, convective = copy_convective
	else
	{
		cb(0, 0) = cd->convective[0][0];
		cb(1, 0) = cd->convective[0][1];
		cb(2, 0) = cd->convective[0][2];
		cb(3, 0) = cd->convective[0][3];
	}

	if (cd->degenerated[0] == true)
	{
		//Correction in case of degeneration - to ensure the correct degenerated convective coordinates
		Matrix temp_coordinates(4);
		for (int index = 0; index < 4; index++)
			temp_coordinates(index, 0) = cb(index, 0);
		//Transforming into degeneration basis
		temp_coordinates = transp(*cd->P[0])*temp_coordinates;
		//Fixing degenerated coordinates
		for (int index = 0; index < 4; index++)
			if (cd->deg_control[0][index] == true)
				temp_coordinates(index, 0) = cd->copy_deg_coordinates[0][index];
		//Transforming into original basis
		temp_coordinates = (*cd->P[0])*temp_coordinates;
		//Copying info into cNR1 vector
		for (int index = 0; index < 4; index++)
			cb(index, 0) = temp_coordinates(index, 0);
	}

	////////////////////////////////////////////////////////////////
	////////////////////////////LCP solution////////////////////////
	////////////////////////////////////////////////////////////////

	//Search always point of minimum - but not can be a point of intersection by conta of the method Barrier

	//bool intersection = false;
	if (interface_1_flag)
	{
		//No degeneration
		if (cd->degenerated[0] == false)
		{
			if (write_report)
				fprintf(f_TR_report, "Search for minimum\n");
			converged1 = FindMinimumSolution(cd, &cb, info);
			converged2 = converged1;
		}

		//Degeneration
		else
		{
			if (write_report)
				fprintf(f_TR_report, "Search for minimum\n");
			converged1 = FindMinimumSolutionDegenerated(cd, cd->P_0[0], &cb);
			converged2 = converged1;
		}

	}
	

	if (previous_evaluation == false)
	{
		//Next convective coordinates:
		cd->convective[0][0] = cb(0, 0);
		cd->convective[0][1] = cb(1, 0);
		cd->convective[0][2] = cb(2, 0);
		cd->convective[0][3] = cb(3, 0);
		//Fills return value with 0 or 4
		cd->return_value[0] = VerifyConvectiveRange(cb);
		cd->convergedLCP[0] = converged2;
		if (converged2 == true)
		{
			cd->characterization_index[0] = CharacterizeCriticalPoint(&cb);
		}
		else
		{
			if (cd->return_value[0] == 0)
				cd->return_value[0] = 1;
			//db.myprintf("\nUnconverged LCP between particle %d and surface %d.\n", pA->number, sB->number);
		}
		
		
		if (write_report)
			fprintf(f_TR_report, "SolveLCP return value is %d. LCP converged flag is %d. \n", cd->return_value[0], converged2);
	}
	else
	{
		//Current convective coordinates:
		cd->copy_convective[0][0] = cb(0, 0);
		cd->copy_convective[0][1] = cb(1, 0);
		cd->copy_convective[0][2] = cb(2, 0);
		cd->copy_convective[0][3] = cb(3, 0);
		//Fills return value with 0 or 4
		cd->copy_return_value[0] = VerifyConvectiveRange(cb);
		cd->copy_convergedLCP[0] = converged2;
		if (converged2 == true)
		{
			cd->copy_characterization_index[0] = CharacterizeCriticalPoint(&cb);
		}
		else
		{
			if (cd->copy_return_value[0] == 0)
				cd->copy_return_value[0] = 1;
			db.myprintf("\nUnconverged LCP between particle %d and surface %d.\n", pA->number, sB->number);
		}
			
		
		if (write_report)
			fprintf(f_TR_report, "SolveLCP return value is %d. LCP converged flag is %d. \n", cd->copy_return_value[0], converged2);
	}

	if (write_report)
		fclose(f_TR_report);
}

void ContactParticleBody::ProcessSurfacePair()
{
	eligible = false;
	if (invalid == false)
	{
		SolveLCP();
		EvaluateNormalGap();
	}
	
}


void ContactParticleBody::FinalProcessSurfacePairsExplicit(double t)
{
	eligible = false;
	if (invalid == false)
	{
		SolveLCP();
		EvaluateNormalGap();
	}
	//ALL - function to update tangential gap, etc. (based on AceGen code)
}

void ContactParticleBody::WriteVTK_XMLForces(FILE *f)
{
	//vectors for writing in the formato binary - using the function 'enconde'
	std::vector<float> float_vector;
	std::vector<int> int_vector;

	int count = 0;
	if (eligible)
		count = count + 2;

	if (count != 0)
	{
		//Opens Piece
		fprintf(f, "\t\t<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", count, count);
		//Opens Points
		fprintf(f, "\t\t\t<Points>\n");
		//Opens DataArray
		fprintf(f, "\t\t\t\t<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"binary\">\n");
		float_vector.clear();
		if (eligible)
		{
			float_vector.push_back((float)((*GammaA)(0, 0)));
			float_vector.push_back((float)((*GammaA)(1, 0)));
			float_vector.push_back((float)((*GammaA)(2, 0)));

			float_vector.push_back((float)((*GammaB)(0, 0)));
			float_vector.push_back((float)((*GammaB)(1, 0)));
			float_vector.push_back((float)((*GammaB)(2, 0)));
		}

		fprintf(f, encodeData<float>(float_vector).c_str());
		fprintf(f, "\n");
		//Closes DataArray
		fprintf(f, "\t\t\t\t</DataArray>\n");
		//Closes Points
		fprintf(f, "\t\t\t</Points>\n");
		//Opens Cells
		fprintf(f, "\t\t\t<Cells>\n");
		//Opens DataArray
		fprintf(f, "\t\t\t\t<DataArray type=\"Int32\" Name=\"connectivity\" format=\"binary\">\n");
		int_vector.clear();
		for (int cell = 0; cell < count; cell++)
			int_vector.push_back(cell);
		fprintf(f, encodeData(int_vector).c_str());
		fprintf(f, "\n");
		//Closes DataArray
		fprintf(f, "\t\t\t\t</DataArray>\n");
		//Opens DataArray
		int_vector.clear();
		fprintf(f, "\t\t\t\t<DataArray type=\"Int32\" Name=\"types\" format=\"binary\">\n");
		for (int cell = 0; cell < count; cell++)
			int_vector.push_back(1);
		fprintf(f, encodeData(int_vector).c_str());
		fprintf(f, "\n");
		//Closes DataArray
		fprintf(f, "\t\t\t\t</DataArray>\n");
		//Opens DataArray
		fprintf(f, "\t\t\t\t<DataArray type=\"Int32\" Name=\"offsets\" format=\"binary\">\n");

		int_vector.clear();
		for (int cell = 0; cell < count; cell++)
			int_vector.push_back(cell + 1);
		fprintf(f, encodeData(int_vector).c_str());
		fprintf(f, "\n");
		//Closes DataArray
		fprintf(f, "\t\t\t\t</DataArray>\n");
		//Closes Cells
		fprintf(f, "\t\t\t</Cells>\n");

		//Opens PointData
		fprintf(f, "\t\t\t<PointData Vectors = \"ContactForces\">\n");
		Matrix normal(3);
		float_vector.clear();
		fprintf(f, "\t\t\t\t<DataArray Name = \"Normal\" type = \"Float32\" NumberOfComponents= \"3\" format=\"binary\">\n");
		if (eligible)
		{

			float_vector.push_back(-(float)(fn[0]));
			float_vector.push_back(-(float)(fn[1]));
			float_vector.push_back(-(float)(fn[2]));

			float_vector.push_back(+(float)(fn[0]));
			float_vector.push_back(+(float)(fn[1]));
			float_vector.push_back(+(float)(fn[2]));
		}
		fprintf(f, encodeData(float_vector).c_str());
		fprintf(f, "\n");
		//Closes DataArray
		fprintf(f, "\t\t\t\t</DataArray>\n");

		float_vector.clear();
		fprintf(f, "\t\t\t\t<DataArray Name = \"Friction\" type = \"Float32\" NumberOfComponents= \"3\" format=\"binary\">\n");
		if (eligible)
		{
			float_vector.push_back(-(float)(ft[0]));
			float_vector.push_back(-(float)(ft[1]));
			float_vector.push_back(-(float)(ft[2]));

			float_vector.push_back(+(float)(ft[0]));
			float_vector.push_back(+(float)(ft[1]));
			float_vector.push_back(+(float)(ft[2]));
		}
		fprintf(f, encodeData(float_vector).c_str());
		fprintf(f, "\n");
		//Closes DataArray
		fprintf(f, "\t\t\t\t</DataArray>\n");

		//Closes PointData
		fprintf(f, "\t\t\t</PointData>\n");

		//Closes Piece
		fprintf(f, "\t\t</Piece>\n");
	}
}

//Calculates the inversa of the Hessian
void ContactParticleBody::EvaluateInvertedHessian()
{
	Matrix invHes(4, 4);
	if (cd->P_0[0]->getColumns() != 0)
	{
		//Tolerance for precision of eigenvalues extraction
		double tol_eig = 1e-14;
		Matrix Hes(4, 4);
		Matrix xk(4);
		for (int i = 0; i < 4; i++)
			xk(i, 0) = cd->convective[0][i];
		HessianPhase1(xk,Hes);
		//Hes.print();
		Matrix Hes_minor = transp(*cd->P_0[0])*Hes*(*cd->P_0[0]);
		int order_minor = cd->P_0[0]->getColumns();
		Matrix P(order_minor, order_minor);
		Matrix D(order_minor, order_minor);
		fulleigen1(Hes_minor, P, D, tol_eig);
		//Inversão of the Hessian
		for (int i = 0; i < order_minor; i++)
		{
			D(i, i) = 1.0 / D(i, i);
		}
		invHes = (*cd->P_0[0])*P*D*transp(P)*transp(*cd->P_0[0]);
	}

	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
			cd->invHessian[0][i][j] = invHes(i, j);
	//invHes.print();
}


//Optimization - determination of minimum
bool ContactParticleBody::FindMinimumSolution(SSContactData* c_data, Matrix* solution, int &return_info)
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
		fprintf(f_TR_report, "FindMinimumSolution\n");
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
		fprintf(f_TR_report, "%d\t%.12e\t%.12e\t%.12e\t%.12e\t%.6e\t%.6e\t%.6e\t%c\t%.6e\t%.6e\t%.6e\n", it, xk(0, 0), xk(1, 0), xk(2, 0), xk(3, 0), ob, error, Deltak, c, rhok, actual_reduction, predicted_reduction);
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
			fprintf(f_TR_report, "%d\t%.12e\t%.12e\t%.12e\t%.12e\t%.6e\t%.6e\t%.6e\t%c\t%.6e\t%.6e\t%.6e\n", it, xk(0, 0), xk(1, 0), xk(2, 0), xk(3, 0), ob, error, Deltak, c, rhok, actual_reduction, predicted_reduction);
	}
	if (write_report)
		fprintf(f_TR_report, "Error\t%.6e\tTolerance\t%.6e\n", error, tol_ortho);

	for (int i = 0; i < 4; i++)
		(*solution)(i, 0) = xk(i, 0);
	//Return of the function
	if (error < tol_ortho)
		return true;
	else
		return false;
}

//Optimization - determination of saddle
//bool ContactParticleBody::FindSaddleSolution(SSContactData* c_data, Matrix* solution, int &return_info, bool return_gap)
//{
//	//Data - trust region
//	double Deltamax = 10;			//maximum radius of the trust region allowed
//	double Deltak = 0.01;			//current radius of trust region
//	double etha = 0.15;				//value between 0 and 0.15 - indicates that the approximation and ruim and veta the increment
//	double rhok = 0.0;				//reason between the difference in the objective function and the difference of the approximation used
//	double actual_change = 0.0;
//	double predicted_change = 0.0;
//	double tol_small = tol_small_1;
//	int max_it = max_it_2;			//max iterations
//	//Variables internal of use geral
//	Matrix Hes(4, 4);
//	Matrix Gra(4, 1);
//	Matrix pGra(4, 1);
//	Matrix P(4, 4);
//	Matrix D(4, 4);
//	Matrix cHes(4, 4);
//	Matrix pk(4);
//	Matrix pb(4);
//	Matrix xk(4);
//	//Initialization of the method - initial guess - obtained of the solution previous - problem of minimum distance
//	for (int i = 0; i < 4; i++)
//		xk(i, 0) = (*solution)(i, 0);
//
//	//Error of the minimum distance
//	GradientPhase1(xk, Gra);
//	double ob = ObjectivePhase1(xk);
//	//Hessian of the objective function
//	HessianPhase1(xk, Hes);
//	//Calculating directions main and curvatures main of the objective function
//	cHes = Hes;
//	fulleigen1(cHes, P, D, tol_eig);
//
//	//Stopping criterion
//	double max_eig = -1e100;
//	for (int i = 0; i < 4; i++)
//	{
//		if (D(i, i) > max_eig)
//			max_eig = D(i, i);
//	}
//	double tol_ortho = tol_convective * abs(max_eig);
//
//	bool saddle_basin = false;		//variable boolean that indicates that the bacia of attraction of the saddle was localizada
//	int it = 1;
//	bool conv1 = false;
//	saddle_basin = false;
//
//	if (CharacterizeCriticalPoint(&xk) == 3)
//	{
//		if (write_report)
//			fprintf(f_TR_report, "FindSaddleSolution initial guess is already in saddle basin\n");
//		saddle_basin = true;//saddle basin
//	}
//	else
//	{
//		if (write_report)
//			fprintf(f_TR_report, "FindSaddleSolution initial guess is not in saddle basin\n");
//
//		Geometry* surf1 = db.body_geometries[index1]->ptr_geom[sub_index1];
//		Geometry* surf2 = db.body_geometries[index2]->ptr_geom[sub_index2];
//		Matrix nA(3);
//		surf1->NormalExt(xk(0, 0), xk(1, 0), &nA, true);
//		Matrix nB(3);
//		surf2->NormalExt(xk(2, 0), xk(3, 0), &nB, true);
//		bool normal_fixed = false;
//
//		//Degeneration local
//		Matrix P_0(4, 2);
//		for (int i = 0; i < 4; i++)
//			P_0(i, 0) = P(i, 0);//first eigenvector
//		for (int i = 0; i < 4; i++)
//			P_0(i, 1) = P(i, 1);//second eigenvector
//		return_info = 0;
//		//Search bi-direcional
//		conv1 = FindMinimumGapDegenerated(c_data, &P_0, &xk, return_info, normal_fixed, in the, nB);
//
//		///////
//		////Hessian of the objective function
//		//HessianPhase1(xk, Hes);
//		////Calculating directions main and curvatures main of the objective function
//		//cHes = Hes;
//		//fulleigen1(cHes, P, D, tol_eig);
//		////Degeneration local
//		//for (int i = 0; i < 4; i++)
//		//	P_0(i, 0) = P(i, 0);//first eigenvector
//		//for (int i = 0; i < 4; i++)
//		//	P_0(i, 1) = P(i, 1);//second eigenvector
//		//conv1 = FindMinimumGapDegenerated(c_data, &P_0, &xk, return_info, normal_fixed, in the, nB);
//		///////
//
//		///////
//		////Hessian of the objective function
//		//HessianPhase1(xk, Hes);
//		////Calculating directions main and curvatures main of the objective function
//		//cHes = Hes;
//		//fulleigen1(cHes, P, D, tol_eig);
//		////Degeneration local
//		//for (int i = 0; i < 4; i++)
//		//	P_0(i, 0) = P(i, 0);//first eigenvector
//		//for (int i = 0; i < 4; i++)
//		//	P_0(i, 1) = P(i, 1);//second eigenvector
//		//conv1 = FindMinimumGapDegenerated(c_data, &P_0, &xk, return_info, normal_fixed, in the, nB);
//		///////
//
//
//		if (conv1 == false)
//			return false;
//		else
//		{
//			if (return_gap == true)
//			{
//				//Saves result in solution
//				for (int i = 0; i < 4; i++)
//					(*solution)(i, 0) = xk(i, 0);
//				return true;
//			}
//		}
//	}
//	if (CharacterizeCriticalPoint(&xk) == 3)
//	{
//		if (write_report)
//			fprintf(f_TR_report, "After FindMinimumGap, initial guess is in saddle basin\n");
//	}
//	else
//	{
//		if (write_report)
//			fprintf(f_TR_report, "After FindMinimumGap, initial guess is not in saddle basin\n");
//	}
//
//	//Calculations of error and objective
//	GradientPhase1(xk, Gra);
//	//Error - forcing first input
//	double error = tol_ortho + 1.0;
//	ob = ObjectivePhase1(xk);
//	//Hessian
//	HessianPhase1(xk, Hes);
//	int flag_error = 0;
//	if (write_report)
//		fprintf(f_TR_report, "FindSaddleSolution\n");
//	char c = 'I';
//	if (write_report)
//		fprintf(f_TR_report, "%d\t%.12e\t%.12e\t%.12e\t%.12e\t%.6e\t%.6e\t%.6e\t%c\t%.6e\t%.6e\t%.6e\n", it, xk(0, 0), xk(1, 0), xk(2, 0), xk(3, 0), ob, error, Deltak, c, rhok, actual_change, predicted_change);
//	/////////////////////////////////Step of analysis of the bacia of attraction of the point of saddle/////////////////////////////////////
//	while ((error > tol_ortho || norm(pk) > tol_convective) && it <= max_it)
//	{
//
//		cHes = Hes;
//		pGra = Gra;
//		//Solves linear system
//		pb = fullsystem(cHes, -1.0*pGra, &flag_error);
//		//Limita TR
//		pk = pb;
//		//Limitando TR
//		double norminc = norm(pk);
//		if ((norminc + tol_convective) >= Deltak)
//			pk = (Deltak / norminc)*pk;
//		c = 'N';
//		//////////////////////////UPDATING SOLUTION////////////////////////////////
//		//Calculation of rhok
//		actual_change = ObjectivePhase1(xk) - ObjectivePhase1(xk + pk);
//		predicted_change = -(transp(Gra)*(pk)+0.5*transp(pk)*Hes*(pk))(0, 0);
//		rhok = actual_change / predicted_change;
//
//		if (abs(actual_change / ObjectivePhase1(xk)) < tol_ascent)
//			rhok = 1.0;
//
//		if (predicted_change < tol_small)
//			rhok = 1.0;
//
//		//update - case the model function represente bem the região. This avoids exit of the bacia of attraction of the saddle
//		if (rhok >= 0.9 && rhok <= 1.1)
//		{
//			xk = xk + pk;
//			Deltak = 2 * Deltak;
//			if (Deltak > Deltamax)
//				Deltak = Deltamax;
//		}
//		else
//			Deltak = Deltak / 4;
//		//increments iterations
//		it++;
//		//Calculations that will be used for check output of the loop - or in the next iteration, case not saia
//		//Objective
//		ob = ObjectivePhase1(xk);
//		//Gradient
//		GradientPhase1(xk, Gra);
//		//Hessian
//		HessianPhase1(xk, Hes);
//		error = norm(Gra);
//		if (write_report)
//			fprintf(f_TR_report, "%d\t%.12e\t%.12e\t%.12e\t%.12e\t%.6e\t%.6e\t%.6e\t%c\t%.6e\t%.6e\t%.6e\n", it, xk(0, 0), xk(1, 0), xk(2, 0), xk(3, 0), ob, error, Deltak, c, rhok, actual_change, predicted_change);
//	}
//	if (write_report)
//		fprintf(f_TR_report, "Error\t%.6e\tTolerance\t%.6e\n", error, tol_ortho);
//
//	////////////////////////////////Check final - if the solution and point of saddle/////////////////////////////////////
//	//Saves result in solution
//	for (int i = 0; i < 4; i++)
//		(*solution)(i, 0) = xk(i, 0);
//	//Return of the function
//	if (error < tol_ortho)
//	{
//		//Checks the solution type encontrada
//		int character = CharacterizeCriticalPoint(&xk);
//		if (character == 3)
//			return true;
//		else
//			return false;
//	}
//	else
//		return false;
//}

//Optimization - determination of minimum
bool ContactParticleBody::FindMinimumSolutionDegenerated(SSContactData* c_data, Matrix* P_0, Matrix* solution)
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
		fprintf(f_TR_report, "FindMinimumSolutionDegenerated\n");
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
		fprintf(f_TR_report, "%d\t%.12e\t%.12e\t%.12e\t%.12e\t%.6e\t%.6e\t%.6e\t%c\t%.6e\t%.6e\t%.6e\n", it, xk(0, 0), xk(1, 0), xk(2, 0), xk(3, 0), ob, error, Deltak, c, rhok, actual_reduction, predicted_reduction);
	/////////////////////////////////////////////////////BEGIN///////////////////////////////////////////////
	while ((error > tol_ortho || norm((*P_0)*pk) > tol_convective) && it <= max_it)
	{
		/*P_0->print();
		Gra.print();
		Hes.print();*/
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
			xk = xk + (*P_0)*pk;	//Increments iterations
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
			fprintf(f_TR_report, "%d\t%.12e\t%.12e\t%.12e\t%.12e\t%.6e\t%.6e\t%.6e\t%c\t%.6e\t%.6e\t%.6e\n", it, xk(0, 0), xk(1, 0), xk(2, 0), xk(3, 0), ob, error, Deltak, c, rhok, actual_reduction, predicted_reduction);
	}
	if (write_report)
		fprintf(f_TR_report, "Error\t%.6e\tTolerance\t%.6e\n", error, tol_ortho);

	//Saves result in solution
	for (int i = 0; i < 4; i++)
		(*solution)(i, 0) = xk(i, 0);
	//Return of the function
	if (error < tol_ortho)
		return true;
	else
		return false;
}

//Optimization - determination of minimum
//bool ContactParticleBody::FindMinimumGapDegenerated(SSContactData* c_data, Matrix* P_0, Matrix* solution, int &return_info, bool fixed_normals, Matrix& nA, Matrix& nB)
//{
//	//Data - trust region   
//	double Deltak = minimum_convective_range / 1000;		//current radius of trust region
//	double Deltamax = 100;									//maximum radius of the trust region allowed
//	double etha = 0.15;										//value between 0 and 0.15 - indicates that the approximation and ruim and veta the increment
//	double rhok = 0.0;										//reason between the difference in the objective function and the difference of the approximation used
//	double last_reduction = 0.0;
//	double reduction = 0.0;
//	double actual_reduction = 0;
//	double predicted_reduction = 0;
//	double tol_small = tol_small_1;
//	int max_it = max_it_2;			//Max iterations
//	Matrix Hes(4, 4);
//	Matrix Gra(4, 1);
//	Matrix xk(4);
//	int order = P_0->getColumns();
//	Matrix pk(order);
//	Matrix pb(order);
//	Matrix pc(order);
//	Matrix pGra(order, 1);
//	Matrix deg_Hes(order, order);
//	Matrix deg_Gra(order, 1);
//	Matrix P(order, order);
//	Matrix D(order, order);
//	Matrix cHes(order, order);
//
//	//	//Initialization of the method - initial guess
//	for (int i = 0; i < 4; i++)
//		xk(i, 0) = (*solution)(i, 0);
//
//	//Stopping criterion
//	HessianGap(xk, Hes, fixed_normals, nA, nB);
//	deg_Hes = transp(*P_0)*Hes*(*P_0);
//	fulleigen1(deg_Hes, P, D, tol_eig);
//	double max_eig = -1e100;
//	for (int i = 0; i < order; i++)
//	{
//		if (D(i, i) > max_eig)
//			max_eig = D(i, i);
//	}
//	double tol_ortho = tol_convective * abs(max_eig);
//
//	if (write_report)
//		fprintf(f_TR_report, "FindMinimumGapDegenerated\n");
//	char c = 'I';
//
//	int it = 1;
//	//Objective
//	double ob = Gap(xk, fixed_normals, nA, nB);
//	//Gradient
//	GradientGap(xk, Gra, fixed_normals, nA, nB);
//	//Hessian
//	HessianGap(xk, Hes, fixed_normals, nA, nB);
//	deg_Gra = transp(*P_0)*Gra;
//	deg_Hes = transp(*P_0)*Hes*(*P_0);
//	//Error - forcing first input
//	double error = tol_ortho + 1.0;
//	//Initial guess report
//	if (write_report)
//		fprintf(f_TR_report, "%d\t%.12e\t%.12e\t%.12e\t%.12e\t%.6e\t%.6e\t%.6e\t%c\t%.6e\t%.6e\t%.6e\n", it, xk(0, 0), xk(1, 0), xk(2, 0), xk(3, 0), ob, error, Deltak, c, rhok, actual_reduction, predicted_reduction);
//	/////////////////////////////////////////////////////BEGIN///////////////////////////////////////////////
//	while ((error > tol_ortho || norm((*P_0)*pk) > tol_convective) && it <= max_it)
//	{
//		//Determination of the point of Cauchy
//		double gragra = (transp(deg_Gra)*deg_Gra)(0, 0);
//		double grahesgra = (transp(deg_Gra)*deg_Hes*deg_Gra)(0, 0);
//		double normgra = norm(deg_Gra);
//		Matrix pc;	//direction of the Cauchy point
//		/////////////////Point of Cauchy///////////////////////
//		if (grahesgra <= 0.0)
//			pc = -(Deltak / normgra)*deg_Gra;
//		else
//			pc = -gragra / grahesgra * deg_Gra;
//		double normpc = norm(pc);
//		//Cauchy point outside the TR - use a fraction of it
//		if ((normpc + tol_convective) >= Deltak)
//		{
//			pk = (Deltak / normpc)*pc; //Steep descent
//			c = 'C';
//		}
//		else//Cauchy point inside the TR
//		{
//			//Calculating directions main and curvatures main of the objective function
//			cHes = deg_Hes;
//			fulleigen1(cHes, P, D, tol_eig);
//			//Writing gradient in the directions main
//			pGra = transp(P)*deg_Gra;
//			//Determination of the smallest eigenvalue (min_eig)
//			double min_eig = 1e100;
//			for (int i = 0; i < order; i++)
//			{
//				if (D(i, i) < min_eig)
//					min_eig = D(i, i);
//			}
//			//Construction of the direction of search
//			//Direction of search based in NR - modificada through the smallest eigenvalue
//			zeros(&pb);
//			//If the smallest eigenvalue and smaller or igual the zero (tol_small) - modifies the direction of NR for garantir direction descendente
//			if (min_eig < tol_small)
//			{
//				for (int i = 0; i < order; i++)
//					pb(i, 0) = -pGra(i, 0) / (D(i, i) - (min_eig - abs(min_eig)*tol_ascent));
//			}
//			//If the smallest eigenvalue and larger that zero (tol_small) - direction of NR and chosen
//			else
//			{
//				for (int i = 0; i < order; i++)
//					pb(i, 0) = -pGra(i, 0) / D(i, i);
//			}
//			//Writing direction of search in the coordinates originais
//			pb = P * pb;
//			double normpb = norm(pb);
//			double thetak;
//			//Newton point inside the TR -  use it
//			if (normpb <= Deltak)
//			{
//				pk = pb;
//				c = 'N';
//			}
//
//			else
//			{
//				////////////////////////////Dogleg path//////////////////////////
//				double a1, b1, c1;
//				a1 = norm(pb - pc)*norm(pb - pc);
//				b1 = 2 * (transp(pc)*(pb - pc))(0, 0);
//				c1 = normpc * normpc - Deltak * Deltak;
//				thetak = (-b1 + sqrt(b1*b1 - 4 * a1*c1)) / (2 * a1);
//				//Determination of the path
//				pk = pc + thetak * (pb - pc);
//				c = 'D';
//			}
//		}
//		//////////////////////////UPDATING SOLUTION////////////////////////////////
//		//Calculation of rhok
//		actual_reduction = Gap(xk, fixed_normals, nA, nB) - Gap(xk + (*P_0)*pk, fixed_normals, nA, nB);
//		predicted_reduction = -(transp(deg_Gra)*pk + 0.5*transp(pk)*deg_Hes*pk)(0, 0);
//		rhok = actual_reduction / predicted_reduction;
//
//		if (abs(predicted_reduction / Gap(xk, fixed_normals, in the, nB)) < tol_ascent)
//			rhok = 1.0;
//
//		if (abs(predicted_reduction) < tol_small || abs(actual_reduction) < tol_small)
//			rhok = 1.0;
//
//		if (rhok >= 0.0)
//		{
//			if (rhok < 0.25)//low reduction
//				Deltak = 0.25*norm(pk);//reduce TR
//			else
//			{
//				if (rhok > 0.75 && ((norm(pk) + tol_convective) >= Deltak && (norm(pk) - tol_convective) <= Deltak))//high reduction and testing the limits of the trust region
//				{
//					//augments the radius of TR
//					if (2.0 * Deltak < Deltamax)
//						Deltak = 2.0 * Deltak;
//					else
//						Deltak = Deltamax;
//				}
//			}
//			if (rhok >= etha)
//				xk = xk + (*P_0)*pk;
//			//error = actual_reduction;
//		}
//		else
//			Deltak = Deltak / 2.0;
//
//		//Increments iterations
//		it++;
//		//Objective
//		ob = Gap(xk, fixed_normals, nA, nB);
//		//Gradient
//		GradientGap(xk, Gra, fixed_normals, nA, nB);
//		//Hessian
//		HessianGap(xk, Hes, fixed_normals, nA, nB);
//		deg_Gra = transp(*P_0)*Gra;
//		deg_Hes = transp(*P_0)*Hes*(*P_0);
//		//Error - norm of the gradient
//		error = norm(deg_Gra);
//		//Initial guess report
//		if (write_report)
//			fprintf(f_TR_report, "%d\t%.12e\t%.12e\t%.12e\t%.12e\t%.6e\t%.6e\t%.6e\t%c\t%.6e\t%.6e\t%.6e\n", it, xk(0, 0), xk(1, 0), xk(2, 0), xk(3, 0), ob, error, Deltak, c, rhok, actual_reduction, predicted_reduction);
//	}
//	if (write_report)
//		fprintf(f_TR_report, "Error\t%.6e\tTolerance\t%.6e\n", error, tol_ortho);
//
//
//	//Return of the function
//	if (error < tol_ortho)
//	{
//		//Saves result in solution
//		for (int i = 0; i < 4; i++)
//			(*solution)(i, 0) = xk(i, 0);
//		return true;
//	}
//	else
//		return false;
//}

//Optimization - determination of saddle
//bool ContactParticleBody::FindSaddleSolutionDegenerated(SSContactData* c_data, Matrix* P_0, Matrix* solution, bool return_gap)
//{
//	//Function desenvolvida for degeneration of order 3
//	//Order 4 - call function FindSaddleSolution
//	//Order 2 - call function FindMinimumSolution (not and necessary determine the point of saddle)
//
//	//Data - trust region
//	double Deltamax = 10;			//maximum radius of the trust region allowed
//	double Deltak = 0.01;			//current radius of trust region
//	double etha = 0.15;				//value between 0 and 0.15 - indicates that the approximation and ruim and veta the increment
//	double rhok = 0.0;				//reason between the difference in the objective function and the difference of the approximation used
//	double actual_change = 0.0;
//	double predicted_change = 0.0;
//
//	double tol_small = tol_small_1;
//	int max_it = max_it_2;			//max iterations
//
//	//Variables internal of use geral
//	Matrix Hes(4, 4);
//	Matrix Gra(4, 1);
//	Matrix pGra(3, 1);
//
//	Matrix dHes(3, 3);
//	Matrix dGra(3, 1);
//
//	Matrix P(3, 3);
//	Matrix D(3, 3);
//	Matrix cHes(3, 3);
//
//	Matrix Hes2(4, 4);
//	Matrix Gra2(4, 1);
//
//	Matrix pk(3);
//	Matrix pb(3);
//	Matrix xk(4);
//
//	if (write_report)
//		fprintf(f_TR_report, "FindSaddleSolutionDegenerated\n");
//
//	//Initialization of the method - initial guess - obtained of the solution previous - problem of minimum distance
//	for (int i = 0; i < 4; i++)
//		xk(i, 0) = (*solution)(i, 0);
//
//	//Error of the minimum distance
//	GradientPhase1(xk, Gra);
//	double error = norm(Gra);
//	double ob = ObjectivePhase1(xk);
//
//	bool saddle_basin = false;		//variable boolean that indicates that the bacia of attraction of the saddle was localizada
//	int it = 1;
//
//	//Hessian of the objective function
//	HessianPhase1(xk, Hes);
//	dHes = transp(*P_0)*Hes*(*P_0);
//	//Calculating directions main and curvatures main of the objective function
//	cHes = dHes;
//	fulleigen1(cHes, P, D, tol_eig);
//
//	//Stopping criterion
//	double max_eig = -1e100;
//	for (int i = 0; i < 3; i++)
//	{
//		if (D(i, i) > max_eig)
//			max_eig = D(i, i);
//	}
//	double tol_ortho = tol_convective * abs(max_eig);
//
//	if (CharacterizeCriticalPointDegenerated(&xk, P_0) == 3)
//	{
//		if (write_report)
//			fprintf(f_TR_report, "FindSaddleSolution initial guess is already in saddle basin\n");
//		saddle_basin = true;//saddle basin
//	}
//	else
//	{
//		if (write_report)
//			fprintf(f_TR_report, "FindSaddleSolution initial guess is not in saddle basin\n");
//
//		//Pointers to surfaces
//		Geometry *sA = static_cast<Geometry*>(db.body_geometries[index1]->ptr_geom[sub_index1]);
//		Geometry *sB = static_cast<Geometry*>(db.body_geometries[index2]->ptr_geom[sub_index2]);
//
//		Matrix nA(3);
//		sA->NormalExt(xk(0, 0), xk(1, 0), &nA, true);
//		Matrix nB(3);
//		sB->NormalExt(xk(2, 0), xk(3, 0), &nB, true);
//		bool normal_fixed = false;
//
//		int info = 0;
//
//		//Degeneration local
//		//Direction of search formada through the eigenvector associated to the smallest eigenvalue
//		Matrix P_0local(3, 1);
//		for (int i = 0; i < 3; i++)
//			P_0local(i, 0) = P(i, 0);
//		Matrix composed_P_0 = (*P_0)*P_0local;
//		bool conv1 = FindMinimumGapDegenerated(c_data, &composed_P_0, &xk, info, normal_fixed, in the, nB);
//
//		if (conv1 == false)
//			return false;
//		//else
//		//{
//		//	if (return_gap == true)
//		//	{
//		//		//Saves result in solution
//		//		for (int i = 0; i < 4; i++)
//		//			(*solution)(i, 0) = xk(i, 0);
//		//		return true;
//		//	}
//
//		//	///////////////
//		//	if (sA->deg_u1 == false && sB->deg_u2 == false)
//		//		sA->NormalExt(xk(0, 0), xk(1, 0), &nA, true);
//		//	if (sB->deg_u1 == false && sB->deg_u2 == false)
//		//		sB->NormalExt(xk(2, 0), xk(3, 0), &nB, true);
//
//		//	normal_fixed = true;
//		//	info = 0;
//		//	//Hessian of the objective function
//		//	HessianPhase1(xk, Hes);
//		//	dHes = transp(*P_0)*Hes*(*P_0);
//		//	//Calculating directions main and curvatures main of the objective function
//		//	cHes = dHes;
//		//	fulleigen1(cHes, P, D, tol_eig);
//		//	//Degeneration local
//		//	//Direction of search formada through the eigenvector associated to the smallest eigenvalue
//
//		//	for (int i = 0; i < 3; i++)
//		//		P_0local(i, 0) = P(i, 0);
//		//	composed_P_0 = (*P_0)*P_0local;
//		//	conv1 = FindMinimumGapDegenerated(c_data, &composed_P_0, &xk, info, normal_fixed, in the, nB);
//
//		//	//////////
//
//
//		//}
//	}
//	/*if (CharacterizeCriticalPointDegenerated(&xk, P_0) == 3)
//	{
//		if (write_report)
//			fprintf(f_TR_report[seq_number], "After FindMinimumGap, initial guess is in saddle basin\n");
//	}
//	else
//	{
//		if (write_report)
//			fprintf(f_TR_report[seq_number], "After FindMinimumGap, initial guess is not in saddle basin\n");
//	}*/
//
//	//Calculations of error and objective
//	GradientPhase1(xk, Gra);
//	//Gradient degenerated
//	dGra = transp(*P_0)*Gra;
//	//Hessian
//	HessianPhase1(xk, Hes);
//	//Hessian degenerated
//	dHes = transp(*P_0)*Hes*(*P_0);
//	//Error - forcing first input
//	error = tol_ortho + 1.0;
//	ob = ObjectivePhase1(xk);
//	int flag_error = 0;
//	if (write_report)
//		fprintf(f_TR_report, "FindSaddleSolutionDegenerated\n");
//	char c = 'I';
//	if (write_report)
//		fprintf(f_TR_report, "%d\t%.12e\t%.12e\t%.12e\t%.12e\t%.6e\t%.6e\t%.6e\t%c\t%.6e\t%.6e\t%.6e\n", it, xk(0, 0), xk(1, 0), xk(2, 0), xk(3, 0), ob, error, Deltak, c, rhok, actual_change, predicted_change);
//	/////////////////////////////////Step of analysis of the bacia of attraction of the point of saddle/////////////////////////////////////
//	while ((error > tol_ortho || norm(pk) > tol_convective) && it <= max_it)
//	{
//		//Calculating directions main and curvatures main of the objective function
//		cHes = dHes;
//		//Solves linear system - pb
//		pb = fullsystem(cHes, -1.0*dGra, &flag_error);
//		//Writing direction of search in the coordinates originais
//		pk = (*P_0)*pb;
//		double normpb = norm(pb);
//		//Limita TR
//		double norminc = norm(pk);
//		if ((norminc + tol_convective) >= Deltak)
//			pk = (Deltak / norminc)*pk;
//		c = 'N';
//		//////////////////////////UPDATING SOLUTION////////////////////////////////
//		//Calculation of rhok
//		double actual_change = ObjectivePhase1(xk) - ObjectivePhase1(xk + pk);
//		double predicted_change = -(transp(Gra)*(pk)+0.5*transp(pk)*Hes*(pk))(0, 0);
//		rhok = actual_change / predicted_change;
//
//		if (abs(actual_change / ObjectivePhase1(xk)) < tol_ascent)
//			rhok = 1.0;
//
//		if (predicted_change < tol_small)
//			rhok = 1.0;
//
//		//update - case the model function represente bem the região. This avoids exit of the bacia of attraction of the saddle
//		if (rhok >= 0.9 && rhok <= 1.1)
//		{
//			xk = xk + pk;
//			Deltak = 2 * Deltak;
//			if (Deltak > Deltamax)
//				Deltak = Deltamax;
//		}
//		else
//			Deltak = Deltak / 4;
//		//increments iterations
//		it++;
//		//Calculations that will be used for check output of the loop - or in the next iteration, case not saia
//		//Objective
//		ob = ObjectivePhase1(xk);
//		//Gradient
//		GradientPhase1(xk, Gra);
//		//Gradient degenerated
//		dGra = transp(*P_0)*Gra;
//		//Hessian
//		HessianPhase1(xk, Hes);
//		//Hessian degenerated
//		dHes = transp(*P_0)*Hes*(*P_0);
//		//Error - norm of the gradient
//		error = norm(dGra);
//		if (write_report)
//			fprintf(f_TR_report, "%d\t%.12e\t%.12e\t%.12e\t%.12e\t%.6e\t%.6e\t%.6e\t%c\t%.6e\t%.6e\t%.6e\n", it, xk(0, 0), xk(1, 0), xk(2, 0), xk(3, 0), ob, error, Deltak, c, rhok, actual_change, predicted_change);
//	}
//	if (write_report)
//		fprintf(f_TR_report, "Error\t%.6e\tTolerance\t%.6e\n", error, tol_ortho);
//
//
//	////////////////////////////////Check final - if the solution and point of saddle/////////////////////////////////////
//	//Saves result in solution
//	for (int i = 0; i < 4; i++)
//		(*solution)(i, 0) = xk(i, 0);
//	//Return of the function
//	if (error < tol_ortho)
//	{
//		//Checks the solution type encontrada
//		int character = CharacterizeCriticalPointDegenerated(&xk, P_0);
//		if (character == 3)
//			return true;
//		else
//		{
//			//Tests adicionais - converged but the saddle not was encontrada
//			return false;
//		}
//	}
//	else
//		return false;
//}


int ContactParticleBody::CharacterizeCriticalPoint(Matrix* solution)
{
	///////////////////////////////Pointers and variables to make access easier//////////////////////////////////////////
	//Surface* surf1 = db.surfaces[surf1_ID - 1];		//Pointer to the surface 1
	//Surface* surf2 = db.surfaces[surf2_ID - 1];		//Pointer to the surface 2

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
			fprintf(f_TR_report, "Intersection found (tolerance %.6e). Eigenvalues are %.6e\t%.6e\t%.6e\t%.6e\n", tol_intersect, D(0, 0), D(1, 1), D(2, 2), D(3, 3));
		return 1;
	}
	//minimum strict
	if (D(0, 0) >= tol_small && D(1, 1) >= tol_small && D(2, 2) >= tol_small && D(3, 3) >= tol_small)
	{
		if (write_report)
			fprintf(f_TR_report, "Strict minimum found. Eigenvalues are %.6e\t%.6e\t%.6e\t%.6e\n", D(0, 0), D(1, 1), D(2, 2), D(3, 3));
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
			fprintf(f_TR_report, "Saddle found. Eigenvalues are %.6e\t%.6e\t%.6e\t%.6e\n", D(0, 0), D(1, 1), D(2, 2), D(3, 3));
		return 3;
	}
	if (write_report)
		fprintf(f_TR_report, "Other critical point found. Eigenvalues are %.6e\t%.6e\t%.6e\t%.6e. Tolerance is %.6e.\n", D(0, 0), D(1, 1), D(2, 2), D(3, 3), tol_small);
	return 4;
}

int ContactParticleBody::CharacterizeCriticalPointDegenerated(Matrix* solution, Matrix* P_0, bool print)
{
	///////////////////////////////Pointers and variables to make access easier//////////////////////////////////////////
	//Surface* surf1 = db.surfaces[surf1_ID - 1];		//Pointer to the surface 1
	//Surface* surf2 = db.surfaces[surf2_ID - 1];		//Pointer to the surface 2

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
			fprintf(f_TR_report, "Intersection found. Eigenvalues are %.6e\t%.6e\t%.6e\n", D(0, 0), D(1, 1), D(2, 2));
		return 1;
	}
	//minimum strict
	if (D(0, 0) >= tol_small && D(1, 1) >= tol_small && D(2, 2) >= tol_small)
	{
		if (write_report)
			fprintf(f_TR_report, "Strict minimum found. Eigenvalues are %.6e\t%.6e\t%.6e. Tolerance is %.6e.\n", D(0, 0), D(1, 1), D(2, 2), tol_small);
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
			fprintf(f_TR_report, "Saddle found. Eigenvalues are %.6e\t%.6e\t%.6e\n", D(0, 0), D(1, 1), D(2, 2));
		return 3;
	}
	if (write_report)
		fprintf(f_TR_report, "Other critical point found. Eigenvalues are %.6e\t%.6e\t%.6e\n", D(0, 0), D(1, 1), D(2, 2));
	return 4;
}

void ContactParticleBody::OpenTRReport()
{
	//Writing the name of the contact pair
	char pair_name[100];
	sprintf(pair_name, "TR_report_particle_%d_surf_%d", db.particles[index1]->number, db.body_geometries[index2]->ptr_geom[sub_index2]->number);
	strcpy(name, db.folder_name);	//folder of the job
	strcat(name, "TR/");			//directory TR
	GIRAFFE_MKDIR(name);					//creating directory TR
	strcat(name, pair_name);		//file name
	strcat(name, ".txt");			//creating file
	f_TR_report = fopen(name, "a");	//opens file
}

void ContactParticleBody::InitializeTRReport()
{
	if (write_report)
	{
		//Writing the name of the contact pair
		char pair_name[100];
		sprintf(pair_name, "TR_report_particle_%d_surf_%d", db.particles[index1]->number, db.body_geometries[index2]->ptr_geom[sub_index2]->number);
		strcpy(name, db.folder_name);	//folder of the job
		strcat(name, "TR/");			//directory TR
		GIRAFFE_MKDIR(name);					//creating directory TR
		strcat(name, pair_name);		//file name
		strcat(name, ".txt");			//creating file
		f_TR_report = fopen(name, "w");	//opens file

		fprintf(f_TR_report, "///////////////////////////////////////////////////////////////////////\n");
		fprintf(f_TR_report, "InitializeConvectiveRange\n");
		fprintf(f_TR_report, "Convective 1:\t%.6e\tto\t%.6e.\tRange:\t%.6e\n", (*convective_min)(0, 0), (*convective_max)(0, 0), (*convective_range)(0, 0));
		fprintf(f_TR_report, "Convective 2:\t%.6e\tto\t%.6e.\tRange:\t%.6e\n", (*convective_min)(1, 0), (*convective_max)(1, 0), (*convective_range)(1, 0));
		fprintf(f_TR_report, "Convective 3:\t%.6e\tto\t%.6e.\tRange:\t%.6e\n", (*convective_min)(2, 0), (*convective_max)(2, 0), (*convective_range)(2, 0));
		fprintf(f_TR_report, "Convective 4:\t%.6e\tto\t%.6e.\tRange:\t%.6e\n", (*convective_min)(3, 0), (*convective_max)(3, 0), (*convective_range)(3, 0));
		fclose(f_TR_report);

		fclose(f_TR_report);			//closes file	
	}
}
