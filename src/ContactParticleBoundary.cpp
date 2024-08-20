#include "ContactParticleBoundary.h"
#include "Database.h"

//Variáveis globais
extern
Database db;


ContactParticleBoundary::ContactParticleBoundary()
{
}


ContactParticleBoundary::~ContactParticleBoundary()
{
}

void ContactParticleBoundary::Clear()
{
	for (int i = 0; i < contact_pairs.size(); i++)
		contact_pairs[i]->SetUnnactive();
}

bool ContactParticleBoundary::HaveErrors()
{
	//Verifications of unconverged contacts (LCP or other issues)
	for (int i = 0; i < contact_pairs.size(); i++)
	{
		if (contact_pairs[i]->GetActive())
			if (contact_pairs[i]->HaveErrors())
				return true;
	}

	return false;
}

void ContactParticleBoundary::SaveConfiguration()
{
	//Cleaning contact pairs
	for (int i = 0; i < contact_pairs.size(); i++)
	{
		contact_pairs[i]->prev_eligible = contact_pairs[i]->eligible;

		if (contact_pairs[i]->eligible == false)
		{
			contact_pairs[i]->Free();
		}
		if (contact_pairs[i]->GetActive())
		{
			//Zera o gap tangencial se:
			//1 - gn for positivo (não há contato)
			//2 - copy_gn for positivo (é a primeira ocorrência de contato - não há como acumular gap tangencial pois o contato acaba de começar)
			//3 - return value é 2 - contato não estabelecido - e não é strong candidate, mas pode ocorrer em próximos instantes
			if (contact_pairs[i]->eligible == false)
			{
				zeros(contact_pairs[i]->cd->g_t[0]);
				zeros(contact_pairs[i]->cd->copy_g_t[0]);
				contact_pairs[i]->td->steps_count_impact = 0;
			}
			else
			{
				*contact_pairs[i]->cd->copy_g_t[0] = *contact_pairs[i]->cd->g_t[0];
				contact_pairs[i]->td->steps_count_impact++;
			}
				

			contact_pairs[i]->cd->copy_g_n[0] = contact_pairs[i]->cd->g_n[0];
			*contact_pairs[i]->cd->copy_g[0] = *contact_pairs[i]->cd->g[0];
			*contact_pairs[i]->cd->copy_n[0] = *contact_pairs[i]->cd->n[0];
			contact_pairs[i]->cd->copy_return_value[0] = contact_pairs[i]->cd->return_value[0];
			contact_pairs[i]->cd->copy_convective[0][0] = contact_pairs[i]->cd->convective[0][0];
			contact_pairs[i]->cd->copy_convective[0][1] = contact_pairs[i]->cd->convective[0][1];
			contact_pairs[i]->cd->copy_convective[0][2] = contact_pairs[i]->cd->convective[0][2];
			contact_pairs[i]->cd->copy_convective[0][3] = contact_pairs[i]->cd->convective[0][3];
			contact_pairs[i]->cd->copy_degenerated[0] = contact_pairs[i]->cd->degenerated[0];
			contact_pairs[i]->cd->copy_stick[0] = contact_pairs[i]->cd->stick[0];
		}
	}
}

void ContactParticleBoundary::MountContacts()
{
	//#pragma omp parallel
	{
		//#pragma omp for
		for (int i = 0; i < contact_pairs.size(); i++)
		{
			if (contact_pairs[i]->eligible)
			{
				contact_pairs[i]->Alloc();
				contact_pairs[i]->SetVariables();					//Sets variables for next evaluations
				contact_pairs[i]->EvaluateInvertedHessian();
				contact_pairs[i]->MountLocalContributions();		//Local contact contributions
			}
		}
	}
}

void ContactParticleBoundary::MountContactsExplicit(double t)
{
	for (int i = 0; i < contact_pairs.size(); i++)
	{
		if (contact_pairs[i]->eligible)
		{
			contact_pairs[i]->Alloc();
			contact_pairs[i]->SetVariablesExplicit(t);					//Sets variables for next evaluations
			contact_pairs[i]->EvaluateInvertedHessian();
			contact_pairs[i]->MountLocalContributionsExplicit(t);		//Local contact contributions
		}
	}
}

bool ContactParticleBoundary::NightOwlContact()
{
	ProcessSurfacePairs();
	//Checking if some of the contact pairs is eligible. If yes, return true.
	for (int i = 0; i < contact_pairs.size(); i++)
	{
		//If contact is detected (eligible pair)
		if (contact_pairs[i]->eligible == true)
		{
			db.myprintf("\nParticle %d and Boundary %d, %d and %d\n", contact_pairs[i]->index1, contact_pairs[i]->index2, contact_pairs[i]->sub_index1, contact_pairs[i]->sub_index2);
			db.myprintf("Gap value %.6e\n", contact_pairs[i]->cd[0].g_n[0]);
			db.myprintf("Prev Gap value %.6e\n", contact_pairs[i]->cd[0].copy_g_n[0]);
			//Bounding volume 1
			if (typeid(*db.particles[index1]->sub_bv[sub_index1]) == typeid(BoundingSphere))
				db.myprintf("BoundingSphere\n");
			if (typeid(*db.particles[index1]->sub_bv[sub_index1]) == typeid(BoundingCylinder))
			{
				db.myprintf("BoundingCylinder\n");
				BoundingCylinder* temp = static_cast<BoundingCylinder*>(db.particles[index1]->sub_bv[sub_index1]);
				temp->Report();
			}

			if (typeid(*db.particles[index1]->sub_bv[sub_index1]) == typeid(BoundingTriangularBox))
				db.myprintf("BoundingTriangularBox\n");
			
			//Bounding volume 2
			if (typeid(*db.boundaries[index2]->sub_bv[sub_index2]) == typeid(BoundingSphere))
				db.myprintf("BoundingSphere\n");
			if (typeid(*db.boundaries[index2]->sub_bv[sub_index2]) == typeid(BoundingCylinder))
			{
				db.myprintf("BoundingCylinder\n");
				BoundingCylinder* temp = static_cast<BoundingCylinder*>(db.particles[index2]->sub_bv[sub_index2]);
				temp->Report();
			}
			if (typeid(*db.boundaries[index2]->sub_bv[sub_index2]) == typeid(BoundingTriangularBox))
				db.myprintf("BoundingTriangularBox\n");

			return true;
		}
		else
		{
			//Error reported - inversion of normal - even not eligible, but "passing through"
			if (contact_pairs[i]->HaveErrors())
				return true;
		}

	}
	return false;
}

void ContactParticleBoundary::WriteVTK_XMLForces(FILE *f)
{
	//vetores para escrita no formato binário - usando a função 'enconde'
	std::vector<float> float_vector;
	std::vector<int> int_vector;

	int count = 0;
	bool action_only = true;
	int n;
	if (action_only)
		n = 1;
	else
		n = 2;

	for (int i = 0; i < contact_pairs.size(); i++)
		if (contact_pairs[i]->eligible)
			count = count + n;

	if (count != 0)
	{
		//Opens Piece
		fprintf(f, "\t\t<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", count, count);
		//Opens Points
		fprintf(f, "\t\t\t<Points>\n");
		//Opens DataArray
		fprintf(f, "\t\t\t\t<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"binary\">\n");
		float_vector.clear();
		for (int i = 0; i < contact_pairs.size(); i++)
		{
			if (contact_pairs[i]->eligible)
			{
				float_vector.push_back((float)((*contact_pairs[i]->GammaA)(0, 0)));
				float_vector.push_back((float)((*contact_pairs[i]->GammaA)(1, 0)));
				float_vector.push_back((float)((*contact_pairs[i]->GammaA)(2, 0)));
				if (!action_only)
				{
					float_vector.push_back((float)((*contact_pairs[i]->GammaB)(0, 0)));
					float_vector.push_back((float)((*contact_pairs[i]->GammaB)(1, 0)));
					float_vector.push_back((float)((*contact_pairs[i]->GammaB)(2, 0)));
				}
			}
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
		for (int i = 0; i < contact_pairs.size(); i++)
		{
			if (contact_pairs[i]->eligible)
			{

				float_vector.push_back(-(float)(contact_pairs[i]->fn[0]));
				float_vector.push_back(-(float)(contact_pairs[i]->fn[1]));
				float_vector.push_back(-(float)(contact_pairs[i]->fn[2]));
				if (!action_only)
				{
					float_vector.push_back(+(float)(contact_pairs[i]->fn[0]));
					float_vector.push_back(+(float)(contact_pairs[i]->fn[1]));
					float_vector.push_back(+(float)(contact_pairs[i]->fn[2]));
				}
			}
		}
		fprintf(f, encodeData(float_vector).c_str());
		fprintf(f, "\n");
		//Closes DataArray
		fprintf(f, "\t\t\t\t</DataArray>\n");

		float_vector.clear();
		fprintf(f, "\t\t\t\t<DataArray Name = \"Friction\" type = \"Float32\" NumberOfComponents= \"3\" format=\"binary\">\n");
		for (int i = 0; i < contact_pairs.size(); i++)
		{
			if (contact_pairs[i]->eligible)
			{
				float_vector.push_back(-(float)(contact_pairs[i]->ft[0]));
				float_vector.push_back(-(float)(contact_pairs[i]->ft[1]));
				float_vector.push_back(-(float)(contact_pairs[i]->ft[2]));
				if (!action_only)
				{
					float_vector.push_back(+(float)(contact_pairs[i]->ft[0]));
					float_vector.push_back(+(float)(contact_pairs[i]->ft[1]));
					float_vector.push_back(+(float)(contact_pairs[i]->ft[2]));
				}
			}
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

double ContactParticleBoundary::TimeStepControl(double kin)
{
	double return_step = db.solution[db.current_solution_number - 1]->end_time;				//valor alto
	for (int i = 0; i < contact_pairs.size(); i++)
	{
		if (typeid(*db.solution[db.current_solution_number - 1]) == typeid(Dynamic))
		{
			Dynamic* ptr_sol = static_cast<Dynamic*>(db.solution[db.current_solution_number - 1]);
			//Setting time step variables (impact control)
			contact_pairs[i]->td->n_steps_impact = ptr_sol->n_steps_impact;
			//*contact_pairs[i]->impactcontrol = true;
		}
		contact_pairs[i]->PredictorTimeStep(kin);
		double try_step = contact_pairs[i]->td->time_step_impact;
		//db.myprintf("try step %.6e\n", try_step);
		if (try_step < return_step)
			return_step = try_step;
	}
	return return_step;
}