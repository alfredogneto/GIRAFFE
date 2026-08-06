#include "SpringDashpot_1.h"

#include "Matrix.h"
#include "Node.h"
#include "Encoding.h"
#include "Dynamic.h"
#include"Database.h"

#define PI 3.1415926535897932384626433832795

//Global variables
extern
Database db;

SpringDashpot_1::SpringDashpot_1()
{
	strain_energy = 0.0;
	kinetic_energy = 0.0;
	potential_gravitational_energy = 0.0;

	VTK_type = 3;
	nDOFs = 6;
	material = 0;
	section = 0;
	n_nodes = 2;
	number = 0;
	nodes = new int[n_nodes];
	VTK_nodes = new int[n_nodes];
	VTK_nodes[0] = 0;
	VTK_nodes[1] = 1;
	
	DOFs = new int *[n_nodes];
	for (int i = 0; i < n_nodes;i++)
		DOFs[i] = new int[db.number_GLs_node];

	type_name = new char[20];//Name of the type of the element
	sprintf(type_name, "SpringDashpot_1");

	//Routine for activate the GLS of each node of the element
	for (int i = 0; i < n_nodes; i++)
	{
		for (int j = 0; j < db.number_GLs_node; j++)
		{
			DOFs[i][j] = 0;
		}
		DOFs[i][0] = 1;
		DOFs[i][1] = 1;
		DOFs[i][2] = 1;
	}

	//Zeroing coefficients of the element
	k = 0.0;
	c = 0.0;
	initial_distance = 0.0;

	gn = 0.0;
	elastic_force = 0.0;
	damping_force = 0.0;

	c_stiffness = Matrix(6, 6);									//Stiffness matrix
	c_damping = Matrix(6, 6);									//Damping matrix
	c_damping_modal = Matrix(6, 6);								//Damping matrix
	c_stiffness_force = Matrix(6, 1);							//Vector of internal forces elasticos
	c_damping_force = Matrix(6, 1);								//Vector of internal forces of damping
	I3 = Matrix(3, 3);
	I3(0, 0) = 1.0;
	I3(1, 1) = 1.0;
	I3(2, 2) = 1.0;

	z1 = Matrix(3, 1);
	z2 = Matrix(3, 1);
	xd1 = Matrix(3, 1);
	xd2 = Matrix(3, 1);

	z1z2 = Matrix(3, 1);		//distance current between nodes
	n = Matrix(3, 1);			//normal direction of the spring
	last_n = Matrix(3, 1);		//last normal direction converged of the spring
	first_evaluation = true;
	non = Matrix(3, 3);
	f = 0.0;
	C1 = Matrix(3, 3);
	C2 = Matrix(3, 3);
	C3 = Matrix(3, 3);
}

SpringDashpot_1::~SpringDashpot_1()
{
	delete[] nodes;
	delete[] VTK_nodes;
	if (DOFs != NULL)
	{
		for (int i = 0; i < n_nodes; i++)
			delete[] DOFs[i];
		delete[] DOFs;
	}
	delete[]type_name;
}

//Checks inconsistencies in the element for avoid errors of execution
bool SpringDashpot_1::Check()
{
	for (int i = 0; i < n_nodes; i++)
	{
		if (nodes[i] > db.number_nodes)
			return false;
	}
	return true;
}

bool SpringDashpot_1::Read(FILE *f)
{
	char s[1000];
	fscanf(f, "%s", s);
	number = atoi(s);

	fscanf(f, "%s", s);
	if (!strcmp(s, "Stiffness"))
	{
		fscanf(f, "%s", s);
		k = atof(s);
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "Damping"))
	{
		fscanf(f, "%s", s);
		c = atof(s);
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "Nodes"))
	{
		for (int n = 0; n < n_nodes; n++)
		{
			fscanf(f, "%s", s);
			nodes[n] = atoi(s);
		}
	}
	else
		return false;
	return true;
}

void SpringDashpot_1::Write(FILE *f)
{
	fprintf(f, "SpringDashpot_1\t%d\tStiffness\t%.6e\tDamping\t%.6e\tNodes\t%d\t%d\n",
		number,
		k,
		c,
		nodes[0],
		nodes[1]);
}
//Writes file of results
void SpringDashpot_1::WriteResults(FILE *f)
{
	fprintf(f, "SpringDashpot_1\t%d\tDeformation\t%.6e\tElasticForce\t%.6e\tDampingForce\t%.6e\n",
		number, gn, elastic_force, damping_force);
}

//Writes in the monitor of the element//Writes in the monitor of the element
void SpringDashpot_1::WriteMonitor(FILE *f, bool first_record, double time)
{
	//Header
	if (first_record == true)
		fprintf(f, "TIME\tDeformation\tElasticForce\tDampingForce\n");
	//Information the to be saved
	fprintf(f, "%.6e\t%.6e\t%.6e\t%.6e\n",
		time,gn,elastic_force,damping_force);
}

void SpringDashpot_1::WriteVTK_XMLBase(std::vector<float> *float_vector)
{
	//Prints the results of the element
	int res_element = 3;
	float_vector->push_back((float)gn);
	float_vector->push_back((float)elastic_force);
	float_vector->push_back((float)damping_force);
	//Prints values null for results that not fazem sentido for this element type
	for (int i = res_element; i < db.n_element_results; i++)
		float_vector->push_back(0.0);
}
void SpringDashpot_1::WriteVTK_XMLRender(FILE *f)
{
	//vectors for writing in the formato binary - using the function 'enconde'
	std::vector<float> float_vector;
	std::vector<int> int_vector;

	//Plotting of the spring, if the stiffness not for nula
	if (k != 0.0)
	{
		int n;								/*number of espiras*/
		float r, L, l;						/*radius of the espira, length of the spring without the extremidades and distance between espiras*/
		float xa, ya, za, xb, yb, zb;       /*Points of the extremidades*/
		float d;                            /*Length total*/
		float D;                            /*Length total of the spring relaxada*/
		float fr, fl;                       /*factors radius r/d and length of the extremidade c/d*/
		int numpontos, i, offsets;			/*number of points OF THE ESPIRAS*/
		int tamanho;						/*size of the vectors*/
		float *x, *y, *z;					/*points of the espiras*/
		double teta;						/*angle in rad*/
		float c_i;							/*length of the extremidades*/
		float *ro;							/*distance of the point the origin - change of coordinates*/
		float alfa, beta;                   /*angles rotation in the eixos x and y - change of coordinates*/
		float *xg, *yg, *zg;				/*points of the espiras in the system global*/

		//Assignment of parameters:
		n = 10;
		xa = (float)db.nodes[nodes[0] - 1]->copy_coordinates[0];
		ya = (float)db.nodes[nodes[0] - 1]->copy_coordinates[1];
		za = (float)db.nodes[nodes[0] - 1]->copy_coordinates[2];
		xb = (float)db.nodes[nodes[1] - 1]->copy_coordinates[0];
		yb = (float)db.nodes[nodes[1] - 1]->copy_coordinates[1];
		zb = (float)db.nodes[nodes[1] - 1]->copy_coordinates[2];
		fr = (float)0.1;
		fl = (float)0.2;
		D = (float)initial_distance;

		tamanho = (20 * n) + 6;
		x = new float[tamanho];
		y = new float[tamanho];
		z = new float[tamanho];
		ro = new float[tamanho];
		xg = new float[tamanho];
		yg = new float[tamanho];
		zg = new float[tamanho];

		numpontos = (20 * n) + 1;

		d = sqrt((xb - xa)*(xb - xa) + (yb - ya)*(yb - ya) + (zb - za)*(zb - za));
		r = fr*D;
		c_i = fl*D;
		L = d - 2 * c_i;
		l = L / n;


		teta = 0;
		for (i = 2; i<numpontos + 2; i++)
		{
			x[i] = r*(float)cos(teta);
			teta = teta + (PI / 10); /*18º in radianos*/
		}
		x[0] = x[1] = x[numpontos + 2] = x[numpontos + 3] = 0;


		teta = 0;
		for (i = 2; i<numpontos + 2; i++)
		{
			y[i] = r*(float)sin(teta);
			teta = teta + (PI / 10); /*18º in radianos*/
		}
		y[0] = y[1] = y[numpontos + 2] = y[numpontos + 3] = 0;



		for (i = 3; i<numpontos + 2; i++)
		{
			z[0] = 0;
			z[1] = c_i;
			z[2] = z[1];
			z[i] = z[i - 1] + (l / 20);
		}
		z[numpontos + 2] = z[numpontos + 1];
		z[numpontos + 3] = d;

		/*Change of coordinates*/
		alfa = atan2((ya - yb), (zb - za));
		beta = asin((xb - xa) / d);
		for (i = 0; i<numpontos + 4; i++)
		{
			ro[i] = sqrt((x[i] * x[i]) + (y[i] * y[i]) + (z[i] * z[i]));
			xg[i] = xa + x[i] * cos(beta) + z[i] * sin(beta);
			yg[i] = ya + x[i] * sin(alfa)*sin(beta) + y[i] * cos(alfa) - z[i] * sin(alfa)*cos(beta);
			zg[i] = za - x[i] * cos(alfa)*sin(beta) + y[i] * sin(alfa) + z[i] * cos(alfa)*cos(beta);
		}
		fprintf(f, "     <Piece NumberOfPoints = \"%d\" NumberOfCells = \"%d\">\n", numpontos + 4, numpontos + 3);
		fprintf(f, "         <Points>\n");
		fprintf(f, "             <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"binary\">\n");
		float_vector.clear();
		float_vector.push_back(xg[0]);
		float_vector.push_back(yg[0]);
		float_vector.push_back(zg[0]);
		float_vector.push_back(xg[1]);
		float_vector.push_back(yg[1]);
		float_vector.push_back(zg[1]);
		//fprintf(f, "                 %.6f %.6f %.6f\n", xg[0], yg[0], zg[0]);
		//fprintf(f, "                 %.6f %.6f %.6f\n", xg[1], yg[1], zg[1]);
		for (i = 2; i<numpontos + 2; i++)
		{
			float_vector.push_back(xg[i]);
			float_vector.push_back(yg[i]);
			float_vector.push_back(zg[i]);
			//fprintf(f, "                 %.6f %.6f %.6f\n", xg[i], yg[i], zg[i]);
		}
		float_vector.push_back(xg[numpontos + 2]);
		float_vector.push_back(yg[numpontos + 2]);
		float_vector.push_back(zg[numpontos + 2]);
		float_vector.push_back(xg[numpontos + 3]);
		float_vector.push_back(yg[numpontos + 3]);
		float_vector.push_back(zg[numpontos + 3]);
		fprintf(f, encodeData<float>(float_vector).c_str());
		fprintf(f, "\n");
		//fprintf(f, "                 %.6f %.6f %.6f\n", xg[numpontos + 2], yg[numpontos + 2], zg[numpontos + 2]);
		//fprintf(f, "                 %.6f %.6f %.6f\n", xg[numpontos + 3], yg[numpontos + 3], zg[numpontos + 3]);
		fprintf(f, "             </DataArray>\n");
		fprintf(f, "         </Points>\n");


		fprintf(f, "         <Cells>\n");
		fprintf(f, "             <DataArray type=\"Int32\" Name=\"connectivity\" format=\"binary\">\n");
		int_vector.clear();
		for (i = 0; i<numpontos + 3; i++)
		{
			int_vector.push_back(i);
			int_vector.push_back(i+1);
			//fprintf(f, "                  %d %d\n", i, i + 1);
		}
		fprintf(f, encodeData(int_vector).c_str());
		fprintf(f, "\n");
		fprintf(f, "             </DataArray>\n");
		fprintf(f, "             <DataArray type=\"Int32\" Name=\"types\" format=\"binary\">\n");
		int_vector.clear();
		for (i = 1; i<numpontos + 4; i++)
		{
			int_vector.push_back(3);
			//fprintf(f, "                 3\n");
		}
		fprintf(f, encodeData(int_vector).c_str());
		fprintf(f, "\n");
		fprintf(f, "             </DataArray>\n");
		fprintf(f, "             <DataArray type=\"Int32\" Name=\"offsets\" format=\"binary\">\n");
		offsets = 0;
		int_vector.clear();
		for (i = 1; i<numpontos + 4; i++)
		{
			offsets = offsets + 2;
			int_vector.push_back(offsets);
			//fprintf(f, "                 %d\n", offsets);
		}
		fprintf(f, encodeData(int_vector).c_str());
		fprintf(f, "\n");
		fprintf(f, "             </DataArray>\n");
		fprintf(f, "         </Cells>\n");

		/////////////////////////////////////////////////////////////////////
		//Opens CellData
		fprintf(f, "\t\t\t<CellData FieldData=\"ElementData\">\n");
		float_vector.clear();
		//Opens DataArray
		int n_cells = numpontos + 3;
		fprintf(f, "\t\t\t\t<DataArray Name=\"ElementResults\" type=\"Float32\" NumberOfComponents=\"%d\" format=\"binary\">\n", db.n_element_results);
		for (int cell = 0; cell < n_cells; cell++)
		{
			//Prints the results of the element
			int res_element = 3;
			float_vector.push_back((float)(gn));
			float_vector.push_back((float)(elastic_force));
			float_vector.push_back((float)(damping_force));
			//Prints values null for results that not fazem sentido for this element type
			for (int i = res_element; i < db.n_element_results; i++)
				float_vector.push_back(0.0);
		}
		fprintf(f, encodeData(float_vector).c_str());
		fprintf(f, "\n");
		//Closes DataArray
		fprintf(f, "\t\t\t\t</DataArray>\n");
		int_vector.clear();
		//Opens DataArray
		fprintf(f, "\t\t\t\t<DataArray Name=\"ElementProperties\" type=\"Int32\" NumberOfComponents=\"%d\" format=\"binary\">\n", 4);
		for (int cell = 0; cell < n_cells; cell++)
		{
			int_vector.push_back(5);		//Element ID
			int_vector.push_back(0);
			int_vector.push_back(0);
			int_vector.push_back(0);
		}
		fprintf(f, encodeData(int_vector).c_str());
		fprintf(f, "\n");
		//Closes DataArray
		fprintf(f, "\t\t\t\t</DataArray>\n");
		//Closes CellData
		fprintf(f, "\t\t\t</CellData>\n");
		/////////////////////////////////////////////////////////////////////

		fprintf(f, "     </Piece>\n");

		delete[]x;
		delete[]y;
		delete[]z;
		delete[]ro;
		delete[]xg;
		delete[]yg;
		delete[]zg;
	}

	//Plotting of the amortecedor

	if (c != 0)
	{
		float L;							/*Comprimentoda spring without the extremidades*/
		float xa, ya, za, xb, yb, zb;       /*Points of the extremidades*/
		int n;								/*number of segmentos of the circunferencia*/
		float d;                            /*Length total*/
		float D;                            /*Length total of the spring relaxada*/
		float r, R;							/*Raios of the cilindros smaller and larger*/
		float fR;							/*factor radius of the cilindro larger R/D*/
		float fl;							/*Length of the extremidade c/D*/
		float fd;							/*length cilindro of larger radius dp/D*/
		int numpontos, i, offsets;			/*number of points of the cilindros*/
		int tamanho;						/*size of the vectors*/
		float *x, *y, *z;					/*points of the espiras*/
		double teta;						/*angle in rad*/
		float c_i;							/*length of the extremidades*/
		float dp;							/*length cilindro of larger radius*/
		float *ro;							/*distance of the point the origin - change of coordinates*/
		float alfa, beta;                   /*angles rotation in the eixos x and y - change of coordinates*/
		float *xg, *yg, *zg;				/*points of the espiras in the system global*/

		//Atrubuição of the Parêmetros

		xa = (float)db.nodes[nodes[0] - 1]->copy_coordinates[0];
		ya = (float)db.nodes[nodes[0] - 1]->copy_coordinates[1];
		za = (float)db.nodes[nodes[0] - 1]->copy_coordinates[2];
		xb = (float)db.nodes[nodes[1] - 1]->copy_coordinates[0];
		yb = (float)db.nodes[nodes[1] - 1]->copy_coordinates[1];
		zb = (float)db.nodes[nodes[1] - 1]->copy_coordinates[2];
		D = (float)initial_distance;
		fl = (float)0.2;
		fd = (float)0.2;
		fR = (float)0.08;
		n = 20;

		tamanho = n * 4 + 6;
		x = new float[tamanho];
		y = new float[tamanho];
		z = new float[tamanho];
		ro = new float[tamanho];
		xg = new float[tamanho];
		yg = new float[tamanho];
		zg = new float[tamanho];

		d = sqrt((xb - xa)*(xb - xa) + (yb - ya)*(yb - ya) + (zb - za)*(zb - za));
		c_i = fl*D;
		dp = fd*D;
		R = fR*D;
		r = R / 3;				/*radius of the cilindro smaller*/
		L = d - 2 * c_i;
		numpontos = n * 4;

		teta = 0;
		for (i = 2; i<numpontos / 2 + 2; i++)				/*n repartições*/
		{
			x[i] = r*(float)cos(teta);
			teta = teta + (PI / (n / 2));
		}
		teta = 0;
		for (i = numpontos / 2 + 2; i<numpontos + 2; i++)	/*n repartições*/
		{
			x[i] = R*(float)cos(teta);
			teta = teta + (PI / (n / 2));
		}
		x[0] = x[1] = x[numpontos + 2] = x[numpontos + 3] = 0;


		teta = 0;
		for (i = 2; i<numpontos / 2 + 2; i++)
		{
			y[i] = r*(float)sin(teta);
			teta = teta + (PI / (n / 2));
		}
		teta = 0;
		for (i = numpontos / 2 + 2; i<numpontos + 2; i++)
		{
			y[i] = R*(float)sin(teta);
			teta = teta + (PI / (n / 2));
		}
		y[0] = y[1] = y[numpontos + 2] = y[numpontos + 3] = 0;


		for (i = 2; i<n + 2; i++)
		{
			z[i] = c_i;
			z[i + n] = z[i + (n * 2)] = c_i + L - dp;
			z[i + (n * 3)] = c_i + L;

		}
		z[0] = 0;
		z[1] = c_i;
		z[numpontos + 2] = c_i + L;
		z[numpontos + 3] = d;


		/*Change of coordinates*/

		alfa = atan2((ya - yb), (zb - za));
		beta = asin((xb - xa) / d);
		for (i = 0; i<numpontos + 4; i++)
		{
			ro[i] = sqrt((x[i] * x[i]) + (y[i] * y[i]) + (z[i] * z[i]));
			xg[i] = xa + x[i] * cos(beta) + z[i] * sin(beta);
			yg[i] = ya + x[i] * sin(alfa)*sin(beta) + y[i] * cos(alfa) - z[i] * sin(alfa)*cos(beta);
			zg[i] = za - x[i] * cos(alfa)*sin(beta) + y[i] * sin(alfa) + z[i] * cos(alfa)*cos(beta);
		}

		fprintf(f, "     <Piece NumberOfPoints = \"%d\" NumberOfCells = \"%d\">\n", numpontos + 4, (2 * n) + 5);
		fprintf(f, "         <Points>\n");
		fprintf(f, "             <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"binary\">\n");
		float_vector.clear();
		float_vector.push_back(xg[0]);
		float_vector.push_back(yg[0]);
		float_vector.push_back(zg[0]);
		float_vector.push_back(xg[1]);
		float_vector.push_back(yg[1]);
		float_vector.push_back(zg[1]);
		//fprintf(f, "                 %.6f %.6f %.6f\n", xg[0], yg[0], zg[0]);
		//fprintf(f, "                 %.6f %.6f %.6f\n", xg[1], yg[1], zg[1]);
		for (i = 2; i<numpontos + 2; i++)
		{
			float_vector.push_back(xg[i]);
			float_vector.push_back(yg[i]);
			float_vector.push_back(zg[i]);
			//fprintf(f, "                 %.6f %.6f %.6f\n", xg[i], yg[i], zg[i]);
		}
		//fprintf(f, "                 %.6f %.6f %.6f\n", xg[numpontos + 2], yg[numpontos + 2], zg[numpontos + 2]);
		//fprintf(f, "                 %.6f %.6f %.6f\n", xg[numpontos + 3], yg[numpontos + 3], zg[numpontos + 3]);
		float_vector.push_back(xg[numpontos + 2]);
		float_vector.push_back(yg[numpontos + 2]);
		float_vector.push_back(zg[numpontos + 2]);
		float_vector.push_back(xg[numpontos + 3]);
		float_vector.push_back(yg[numpontos + 3]);
		float_vector.push_back(zg[numpontos + 3]);
		fprintf(f, encodeData<float>(float_vector).c_str());
		fprintf(f, "\n");
		fprintf(f, "             </DataArray>\n");
		fprintf(f, "         </Points>\n");


		fprintf(f, "         <Cells>\n");
		fprintf(f, "             <DataArray type=\"Int32\" Name=\"connectivity\" format=\"binary\">\n");
		int_vector.clear();
		int_vector.push_back(0);
		int_vector.push_back(1);
		int_vector.push_back((4 * n) + 2);
		int_vector.push_back((4 * n) + 3);
		//fprintf(f, "                  %d %d\n", 0, 1);
		//fprintf(f, "                  %d %d\n		", (4 * n) + 2, (4 * n) + 3);
		
		for (i = 2; i <n + 2; i++)
		{
			int_vector.push_back(i);
			//fprintf(f, "%d ", i);
		}
		
		//fprintf(f, "\n		");
		for (i = (2 * n) + 2; i <(3 * n) + 2; i++)
		{
			int_vector.push_back(i);
			//fprintf(f, "%d ", i);
		}
		//fprintf(f, "\n		");
		for (i = (3 * n) + 2; i <(4 * n) + 2; i++)
		{
			int_vector.push_back(i);
			//fprintf(f, "%d ", i);
		}
	
		for (i = 2; i<n + 1; i++)
		{
			int_vector.push_back(n + i);
			int_vector.push_back(n + 1 + i);
			int_vector.push_back(i + 1);
			int_vector.push_back(i);
			//fprintf(f, "                  %d %d %d %d\n", n + i, n + 1 + i, i + 1, i);
		}
		int_vector.push_back(2 * n + 1);
		int_vector.push_back(n + 2);
		int_vector.push_back(2);
		int_vector.push_back(n + 1);
		//fprintf(f, "                  %d %d %d %d\n", 2 * n + 1, n + 2, 2, n + 1);
		for (i = 2 * n + 2; i<3 * n + 1; i++)
		{
			int_vector.push_back(n + i);
			int_vector.push_back(n + 1 + i);
			int_vector.push_back(i + 1);
			int_vector.push_back(i);
			//fprintf(f, "                  %d %d %d %d\n", n + i, n + 1 + i, i + 1, i);
		}
		int_vector.push_back(4 * n + 1);
		int_vector.push_back(3*n + 2);
		int_vector.push_back(2 * n + 2);
		int_vector.push_back(3*n + 1);
		//fprintf(f, "                  %d %d %d %d\n", 4 * n + 1, 3 * n + 2, 2 * n + 2, 3 * n + 1);

		fprintf(f, encodeData(int_vector).c_str());
		fprintf(f, "\n");

		fprintf(f, "             </DataArray>\n");
		fprintf(f, "             <DataArray type=\"Int32\" Name=\"types\" format=\"binary\">\n");
		int_vector.clear();
		int_vector.push_back(3);
		int_vector.push_back(3);
		//fprintf(f, "                 3\n");
		//fprintf(f, "                 3\n");
		for (i = 0; i< 2 * n + 3; i++)
		{
			int_vector.push_back(7);
			//fprintf(f, "                 7\n");
		}
		fprintf(f, encodeData(int_vector).c_str());
		fprintf(f, "\n");
		fprintf(f, "             </DataArray>\n");
		fprintf(f, "             <DataArray type=\"Int32\" Name=\"offsets\" format=\"binary\">\n");
		int_vector.clear();
		int_vector.push_back(2);
		int_vector.push_back(4);
		int_vector.push_back(n + 4);
		int_vector.push_back(2 * n + 4);
		int_vector.push_back(3 * n + 4);
		//fprintf(f, "                 %d\n", 2);
		//fprintf(f, "                 %d\n", 4);
		//fprintf(f, "                 %d\n", n + 4);
		//fprintf(f, "                 %d\n", 2 * n + 4);
		//fprintf(f, "                 %d\n", 3 * n + 4);
		offsets = 3 * n + 4;

		for (i = 0; i <2 * n; i++)
		{
			offsets = offsets + 4;
			int_vector.push_back(offsets);
			//fprintf(f, "                 %d\n", offsets);
		}
		fprintf(f, encodeData(int_vector).c_str());
		fprintf(f, "\n");
		fprintf(f, "             </DataArray>\n");
		fprintf(f, "         </Cells>\n");

		/////////////////////////////////////////////////////////////////////
		//Opens CellData
		fprintf(f, "\t\t\t<CellData FieldData=\"ElementData\">\n");
		float_vector.clear();
		//Opens DataArray
		int n_cells = (2 * n) + 5;
		fprintf(f, "\t\t\t\t<DataArray Name=\"ElementResults\" type=\"Float32\" NumberOfComponents=\"%d\" format=\"binary\">\n", db.n_element_results);
		for (int cell = 0; cell < n_cells; cell++)
		{
			//Prints the results of the element
			int res_element = 3;
			float_vector.push_back((float)(gn));
			float_vector.push_back((float)(elastic_force));
			float_vector.push_back((float)(damping_force));
			//Prints values null for results that not fazem sentido for this element type
			for (int i = res_element; i < db.n_element_results; i++)
				float_vector.push_back(0.0);
		}
		fprintf(f, encodeData(float_vector).c_str());
		fprintf(f, "\n");
		//Closes DataArray
		fprintf(f, "\t\t\t\t</DataArray>\n");
		int_vector.clear();
		//Opens DataArray
		fprintf(f, "\t\t\t\t<DataArray Name=\"ElementProperties\" type=\"Int32\" NumberOfComponents=\"%d\" format=\"binary\">\n", 4);
		for (int cell = 0; cell < n_cells; cell++)
		{
			int_vector.push_back(5);		//Element ID
			int_vector.push_back(0);
			int_vector.push_back(0);
			int_vector.push_back(0);
		}
		fprintf(f, encodeData(int_vector).c_str());
		fprintf(f, "\n");
		//Closes DataArray
		fprintf(f, "\t\t\t\t</DataArray>\n");
		//Closes CellData
		fprintf(f, "\t\t\t</CellData>\n");
		/////////////////////////////////////////////////////////////////////

		fprintf(f, "     </Piece>\n");

		delete[]x;
		delete[]y;
		delete[]z;
		delete[]ro;
		delete[]xg;
		delete[]yg;
		delete[]zg;
	}
}

//Assembles loadings associated to the element
void SpringDashpot_1::MountElementLoads()
{

}

//Assembles elements
void SpringDashpot_1::Mount()
{
	for (int ind = 0; ind < 3; ind++)
	{
		xd1(ind, 0) = db.nodes[nodes[0] - 1]->vel[ind];	//Velocity of the Node 1
		xd2(ind, 0) = db.nodes[nodes[1] - 1]->vel[ind];	//Velocity of the Node 2
		z1(ind, 0) = db.nodes[nodes[0] - 1]->copy_coordinates[ind] + db.nodes[nodes[0] - 1]->displacements[ind];	//Position of the node 1
		z2(ind, 0) = db.nodes[nodes[1] - 1]->copy_coordinates[ind] + db.nodes[nodes[1] - 1]->displacements[ind];	//Position of the node 2
	}
	
	//Calculation of the distance between the points
	z1z2 = z1 - z2;
	//Normal direction
	n = (1.0 / norm(z1z2))*z1z2;
	double norm_z1z2 = norm(z1z2);

	//Calculation of the gap - when not there is inversão of points or for the first time that calculates (until convergir some step of the modelo)
	if (dot(last_n, n) >= 0.0 || first_evaluation == true)
	{
		gn = norm_z1z2 - initial_distance;
		f = (gn / norm_z1z2);
	}
		
	else//there was inversão of the points
	{
		gn = -(norm_z1z2 + initial_distance);
		n = -1.0*n;//Inversão of the normal
		f = gn / (norm_z1z2 + initial_distance);
	}
	
	//Force in the spring
	elastic_force = k*gn;
	//Force in the amortecedor
	damping_force = c*dot(xd1 - xd2, n);

	//Matrices auxiliary for operator tangent
	non = dyadic(n, n);
	if (typeid(*db.solution[db.current_solution_number - 1]) == typeid(Dynamic))
	{
		Dynamic* ptr_sol = static_cast<Dynamic*>(db.solution[db.current_solution_number - 1]);

		C1 = ptr_sol->a4*c*non;
		C2 = (1.0 / norm(z1z2))*c*(dyadic(n, xd1 - xd2)*(I3 - non));
		C3 = (1.0 / norm(z1z2))*c*dot(xd1 - xd2, n)*(I3 - non);
	}
	
	
	for (int l = 0; l < 3; l++)
	{
		c_stiffness_force(l, 0) = elastic_force*n(l, 0);
		c_stiffness_force(l + 3, 0) = -elastic_force*n(l, 0);
		c_damping_force(l, 0) = damping_force*n(l, 0);
		c_damping_force(l + 3, 0) = -damping_force*n(l, 0);

		for (int m = 0; m < 3; m++)
		{
			c_stiffness(l, m) = k*(f*I3(l, m) + (1 - f)*non(l, m));
			c_stiffness(l + 3, m + 3) = k*(f*I3(l, m) + (1 - f)*non(l, m));
			c_stiffness(l + 3, m) = -k*(f*I3(l, m) + (1 - f)*non(l, m));
			c_stiffness(l, m + 3) = -k*(f*I3(l, m) + (1 - f)*non(l, m));
			
			c_damping(l, m) = C1(l, m) + C2(l, m) + C3(l, m);
			c_damping(l + 3, m + 3) = C1(l, m) + C2(l, m) + C3(l, m);
			c_damping(l + 3, m) = -C1(l, m) - C2(l, m) - C3(l, m);
			c_damping(l, m + 3) = -C1(l, m) - C2(l, m) - C3(l, m);
		}
	}
}
//Assembles transformation matrix of coordinates
void SpringDashpot_1::TransformMatrix()
{
	//DOES NOTHING
}
//Zeroes matrices local of the element
void SpringDashpot_1::Zeros()
{
	c_stiffness.clear();									
	c_damping.clear();									
	c_stiffness_force.clear();							
	c_damping_force.clear();
	c_damping_modal.clear();
	kinetic_energy = 0.0;
	strain_energy = 0.0;
	potential_gravitational_energy = 0.0;
}
//Fills the contribution of the element in the matrices global
void SpringDashpot_1::MountGlobal()
{
	int GL_global_1 = 0;
	int GL_global_2 = 0;
	double anterior = 0;
	for (int i = 0; i < 6; i++)
	{
		//Node 1
		if (i<3)
			GL_global_1 = db.nodes[nodes[0] - 1]->GLs[i];
		//Node 2
		else
			GL_global_1 = db.nodes[nodes[1] - 1]->GLs[i - 3];

		//Case the degree of freedom be free:
		if (GL_global_1 > 0)
		{
			anterior = db.global_P_A(GL_global_1 - 1, 0);
			db.global_P_A(GL_global_1 - 1, 0) = anterior + c_stiffness_force(i, 0);
			anterior = db.global_I_A(GL_global_1 - 1, 0);
			db.global_I_A(GL_global_1 - 1, 0) = anterior + c_stiffness_force(i, 0);
		}
		else
		{
			anterior = db.global_P_B(-GL_global_1 - 1, 0);
			db.global_P_B(-GL_global_1 - 1, 0) = anterior + c_stiffness_force(i, 0);
		}
		for (int j = 0; j < 6; j++)
		{
			//Node 1
			if (j<3)
				GL_global_2 = db.nodes[nodes[0] - 1]->GLs[j];
			//Node 2
			else
				GL_global_2 = db.nodes[nodes[1] - 1]->GLs[j - 3];

			//Case the degrees of freedom be both free (Matrix Kaa)
			if (GL_global_1 > 0 && GL_global_2 > 0)
				db.global_stiffness_AA.setValue(GL_global_1 - 1, GL_global_2 - 1, c_stiffness(i, j));
			//Case the degrees of freedom be both fixed (Matrix Kbb)
			if (GL_global_1 < 0 && GL_global_2 < 0)
				db.global_stiffness_BB.setValue(-GL_global_1 - 1, -GL_global_2 - 1, c_stiffness(i, j));
			//Case the degrees of freedom be free and fixed (Matrix Kab)
			if (GL_global_1 > 0 && GL_global_2 < 0)
				db.global_stiffness_AB.setValue(GL_global_1 - 1, -GL_global_2 - 1, c_stiffness(i, j));
			//Case the degrees of freedom be fixed and free (Matrix Kba)
			if (GL_global_1 < 0 && GL_global_2 > 0)
				db.global_stiffness_BA.setValue(-GL_global_1 - 1, GL_global_2 - 1, c_stiffness(i, j));
		}
	}
}
//Saves variables in the gauss points useful for description lagrangiana updated
void SpringDashpot_1::SaveLagrange()
{
	//Saves the direction converged of the spring (for check inversões in the step subsequent)
	last_n = n;
	first_evaluation = false;
}
//Pre-calculation of variables that and done a single time in the start
void SpringDashpot_1::PreCalc()
{
	z1(0, 0) = db.nodes[nodes[0] - 1]->ref_coordinates[0];
	z1(1, 0) = db.nodes[nodes[0] - 1]->ref_coordinates[1];
	z1(2, 0) = db.nodes[nodes[0] - 1]->ref_coordinates[2];

	z2(0, 0) = db.nodes[nodes[1] - 1]->ref_coordinates[0];
	z2(1, 0) = db.nodes[nodes[1] - 1]->ref_coordinates[1];
	z2(2, 0) = db.nodes[nodes[1] - 1]->ref_coordinates[2];

	initial_distance = norm(z1 - z2);
}

//Assembles the mass matrix
void SpringDashpot_1::MountMass()
{
	//DOES NOTHING
}

//Assembles the mass matrix
void SpringDashpot_1::MountMassModal()
{
	Zeros();
	//DOES NOTHING
}

//Assembles the damping matrix for performing of the analysis modal
void SpringDashpot_1::MountDampingModal()
{
	Zeros();
	//ALL
}

//Assembles the damping matrix
void SpringDashpot_1::MountDamping(bool update_rayleigh)
{
	//DOES NOTHING
}

//Assemblies - Newmark
void SpringDashpot_1::MountDyn()
{
	//Modificações of the dynamics
	c_stiffness = c_stiffness + c_damping;
	//Modificações of the dynamics in the internal forces - presença of the forces of damping
	c_stiffness_force = c_stiffness_force + c_damping_force;
}

//Assemblies for analysis modal - insertion of the mass matrix and damping in the stiffness matrix for subsequent assembly global
void SpringDashpot_1::MountDynModal()
{
	c_stiffness = c_damping_modal;
}

