#include "SuperNodalLoad.h"
#include "IO.h"


SuperNodalLoad::SuperNodalLoad()
{
	number = 0;
	table = NULL;
	mcode = NULL;
	n_times = 0;
	n_values = 1;

	localDOF = 0;
	super_node_ID = 0;
	load_value = 0.0;
}

SuperNodalLoad::~SuperNodalLoad()
{
	if (table != NULL)
		delete table;
	if (mcode != NULL)
		delete mcode;
}

//Reads input file
bool SuperNodalLoad::Read(FILE *f)
{
	char s[1000];
	fscanf(f, "%s", s);
	number = atoi(s);
	fscanf(f, "%s", s);
	if (!strcmp(s, "SuperNode"))
	{
		fscanf(f, "%s", s);
		super_node_ID = atoi(s);
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "LocalDOF"))
	{
		fscanf(f, "%s", s);
		localDOF = atoi(s);
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "NTimes"))
	{
		fscanf(f, "%s", s);
		n_times = atoi(s);
		//Allocating table
		table = new Table(n_times, n_values);
		TryComment(f);
		//Reads table data
		if (!table->Read(f))
			return false;
	}
	else
	{
		if (!strcmp(s, "MathCode"))
		{
			//Allocating MathCode
			mcode = new MathCode(n_values);
			//Reads MathCode
			if (!mcode->Read(f))
				return false;
		}
		else
			return false;
	}
	return true;
}

//Checking inconsistencies
bool SuperNodalLoad::Check()
{
	if (super_node_ID > db.number_super_nodes)
		return false;
	return true;
}

//Writes output file
void SuperNodalLoad::Write(FILE *f)
{
	fprintf(f, "SuperNodalLoad\t%d\tSuperNode\t%d\tLocalDOF\t%d\t", number, super_node_ID, localDOF);
	if (table != NULL)
	{
		fprintf(f, "NTimes\t%d\n", n_times);
		table->Write(f);
	}
	if (mcode != NULL)
	{
		fprintf(f, "MathCode\n");
		mcode->Write(f);
	}
}
//Writes VTK XML data for post-processing
void SuperNodalLoad::WriteVTK_XML(FILE *f)
{
	//If the local DOF is valid
	if (localDOF <= db.super_nodes[super_node_ID - 1]->n_DOFs)
	{
		int vertex = 0;
		//localDOF is a displacement DOF
		if (localDOF <= db.super_nodes[super_node_ID - 1]->n_displacement_DOFs)
			vertex = (localDOF - 1) / 3 + 1;
		//localDOF is a temperature DOF
		else
			vertex = (localDOF - db.super_nodes[super_node_ID - 1]->n_displacement_DOFs - 1) / 3 + 1;

		//vetores para escrita no formato binário - usando a função 'enconde'
		std::vector<float> float_vector;
		std::vector<int> int_vector;

		//Opens Piece
		fprintf(f, "\t\t<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", 1, 1);
		//Opens Points
		fprintf(f, "\t\t\t<Points>\n");
		//Opens DataArray
		fprintf(f, "\t\t\t\t<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"binary\">\n");
		float_vector.clear();
		//Preenchendo as coordenadas do vértice do super node
		float_vector.push_back((float)(db.super_nodes[super_node_ID - 1]->copy_coordinates[(vertex - 1) * 3]));
		float_vector.push_back((float)(db.super_nodes[super_node_ID - 1]->copy_coordinates[(vertex - 1) * 3 + 1]));
		float_vector.push_back((float)(db.super_nodes[super_node_ID - 1]->copy_coordinates[(vertex - 1) * 3 + 2]));
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
		int_vector.push_back(0);
		fprintf(f, encodeData(int_vector).c_str());
		fprintf(f, "\n");
		//Closes DataArray
		fprintf(f, "\t\t\t\t</DataArray>\n");
		//Opens DataArray
		int_vector.clear();
		fprintf(f, "\t\t\t\t<DataArray type=\"Int32\" Name=\"types\" format=\"binary\">\n");
		int_vector.push_back(1);
		fprintf(f, encodeData(int_vector).c_str());
		fprintf(f, "\n");
		//Closes DataArray
		fprintf(f, "\t\t\t\t</DataArray>\n");
		//Opens DataArray
		fprintf(f, "\t\t\t\t<DataArray type=\"Int32\" Name=\"offsets\" format=\"binary\">\n");
		int_vector.clear();
		int_vector.push_back(1);
		fprintf(f, encodeData(int_vector).c_str());
		fprintf(f, "\n");
		//Closes DataArray
		fprintf(f, "\t\t\t\t</DataArray>\n");
		//Closes Cells
		fprintf(f, "\t\t\t</Cells>\n");

		//Opens PointData
		fprintf(f, "\t\t\t<PointData Scalars = \"Loads\">\n");

		float_vector.clear();
		//Opens DataArray
		fprintf(f, "\t\t\t\t<DataArray Name = \"Load\" type = \"Float32\" format=\"binary\">\n");
		float_vector.push_back((float)(load_value));
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
//Pre-calculus
void SuperNodalLoad::PreCalc()
{

}

//Atualiza dados necessários e que sejam dependentes de DOFs ativos/inativos - chamado no início de cada solution step
void SuperNodalLoad::UpdateforSolutionStep()
{
	
}

void SuperNodalLoad::Mount()
{
	load_value = GetValueAt(db.last_converged_time + db.current_time_step, 0);
	//Global contributions
	int GL_lin;
	//If the local DOF is valid
	if (localDOF <= db.super_nodes[super_node_ID - 1]->n_DOFs)
		GL_lin = db.super_nodes[super_node_ID - 1]->DOFs[localDOF - 1];
	else
		GL_lin = 0;
	if (GL_lin != 0)
	{
		if (GL_lin > 0)	//Grau de liberdade livre
			db.global_P_A(GL_lin - 1, 0) += -1.0*load_value;
		else
			db.global_P_B(-GL_lin - 1, 0) += -1.0*load_value;
	}
}


void SuperNodalLoad::EvaluateExplicit(double t)
{

}
