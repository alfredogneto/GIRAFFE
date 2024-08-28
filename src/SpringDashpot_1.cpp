#include "SpringDashpot_1.h"

#include "Matrix.h"
#include "Node.h"
#include "Encoding.h"
#include "Dynamic.h"
#include"Database.h"

#define PI 3.1415926535897932384626433832795

//Variaveis globais
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

	type_name = new char[20];//Nome do tipo do elemento
	sprintf(type_name, "SpringDashpot_1");

	//Rotina para ativar os GLS de cada nó do elemento
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

	//Zerando coeficientes do elemento
	k = 0.0;
	c = 0.0;
	initial_distance = 0.0;

	gn = 0.0;
	elastic_force = 0.0;
	damping_force = 0.0;

	c_stiffness = Matrix(6, 6);									//Matriz de rigidez
	c_damping = Matrix(6, 6);									//Matriz de amortecimento
	c_damping_modal = Matrix(6, 6);								//Matriz de amortecimento
	c_stiffness_force = Matrix(6, 1);							//Vetor de esforços elasticos
	c_damping_force = Matrix(6, 1);								//Vetor de esforços de amortecimento
	I3 = Matrix(3, 3);
	I3(0, 0) = 1.0;
	I3(1, 1) = 1.0;
	I3(2, 2) = 1.0;

	z1 = Matrix(3, 1);
	z2 = Matrix(3, 1);
	xd1 = Matrix(3, 1);
	xd2 = Matrix(3, 1);

	z1z2 = Matrix(3, 1);		//distancia atual entre nós
	n = Matrix(3, 1);			//direção normal da mola
	last_n = Matrix(3, 1);		//ultima direção normal convergida da mola
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

//Checa inconsistências no elemento para evitar erros de execução
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
//Escreve arquivo de resultados
void SpringDashpot_1::WriteResults(FILE *f)
{
	fprintf(f, "SpringDashpot_1\t%d\tDeformation\t%.6e\tElasticForce\t%.6e\tDampingForce\t%.6e\n",
		number, gn, elastic_force, damping_force);
}

//Escreve no monitor do elemento//Escreve no monitor do elemento
void SpringDashpot_1::WriteMonitor(FILE *f, bool first_record, double time)
{
	//Cabeçalho
	if (first_record == true)
		fprintf(f, "TIME\tDeformation\tElasticForce\tDampingForce\n");
	//Informações a serem salvas
	fprintf(f, "%.6e\t%.6e\t%.6e\t%.6e\n",
		time,gn,elastic_force,damping_force);
}

void SpringDashpot_1::WriteVTK_XMLBase(std::vector<float> *float_vector)
{
	//Imprime os resultados do elemento
	int res_element = 3;
	float_vector->push_back((float)gn);
	float_vector->push_back((float)elastic_force);
	float_vector->push_back((float)damping_force);
	//Imprime valores nulos para resultados que não fazem sentido para esse tipo de elemento
	for (int i = res_element; i < db.n_element_results; i++)
		float_vector->push_back(0.0);
}
void SpringDashpot_1::WriteVTK_XMLRender(FILE *f)
{
	//vetores para escrita no formato binario - usando a função 'enconde'
	std::vector<float> float_vector;
	std::vector<int> int_vector;

	//Plotagem da mola, se a rigidez não for nula
	if (k != 0.0)
	{
		int n;								/*numero de espiras*/
		float r, L, l;						/*raio da espira, comprimento da mola sem as extremidades e distancia entre espiras*/
		float xa, ya, za, xb, yb, zb;       /*Pontos das extremidades*/
		float d;                            /*Comprimento total*/
		float D;                            /*Comprimento total da mola relaxada*/
		float fr, fl;                       /*fatores raio r/d e comprimento da extremidade c/d*/
		int numpontos, i, offsets;			/*numero de pontos DAS ESPIRAS*/
		int tamanho;						/*tamanho dos vetores*/
		float *x, *y, *z;					/*pontos das espiras*/
		double teta;						/*angulo em rad*/
		float c_i;							/*comprimento das extremidades*/
		float *ro;							/*distancia do ponto a origem - mudança de coordenadas*/
		float alfa, beta;                   /*angulos rotação nos eixos x e y - mudança de coordenadas*/
		float *xg, *yg, *zg;				/*pontos das espiras no sistema global*/

		//Atribuição de parametros:
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
			teta = teta + (PI / 10); /*18º em radianos*/
		}
		x[0] = x[1] = x[numpontos + 2] = x[numpontos + 3] = 0;


		teta = 0;
		for (i = 2; i<numpontos + 2; i++)
		{
			y[i] = r*(float)sin(teta);
			teta = teta + (PI / 10); /*18º em radianos*/
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

		/*Mudança de coordenadas*/
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
			//Imprime os resultados do elemento
			int res_element = 3;
			float_vector.push_back((float)(gn));
			float_vector.push_back((float)(elastic_force));
			float_vector.push_back((float)(damping_force));
			//Imprime valores nulos para resultados que não fazem sentido para esse tipo de elemento
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

	//Plotagem do amortecedor

	if (c != 0)
	{
		float L;							/*Comprimentoda mola sem as extremidades*/
		float xa, ya, za, xb, yb, zb;       /*Pontos das extremidades*/
		int n;								/*numero de segmentos da circunferencia*/
		float d;                            /*Comprimento total*/
		float D;                            /*Comprimento total da mola relaxada*/
		float r, R;							/*Raios dos cilindros menor e maior*/
		float fR;							/*fator raio do cilindro maior R/D*/
		float fl;							/*Comprimento da extremidade c/D*/
		float fd;							/*comprimento cilindro de maior raio dp/D*/
		int numpontos, i, offsets;			/*numero de pontos dos cilindros*/
		int tamanho;						/*tamanho dos vetores*/
		float *x, *y, *z;					/*pontos das espiras*/
		double teta;						/*angulo em rad*/
		float c_i;							/*comprimento das extremidades*/
		float dp;							/*comprimento cilindro de maior raio*/
		float *ro;							/*distancia do ponto a origem - mudança de coordenadas*/
		float alfa, beta;                   /*angulos rotação nos eixos x e y - mudança de coordenadas*/
		float *xg, *yg, *zg;				/*pontos das espiras no sistema global*/

		//Atrubuição dos Parêmetros

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
		r = R / 3;				/*raio do cilindro menor*/
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


		/*Mudança de coordenadas*/

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
			//Imprime os resultados do elemento
			int res_element = 3;
			float_vector.push_back((float)(gn));
			float_vector.push_back((float)(elastic_force));
			float_vector.push_back((float)(damping_force));
			//Imprime valores nulos para resultados que não fazem sentido para esse tipo de elemento
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

//Monta carregamentos associados ao elemento
void SpringDashpot_1::MountElementLoads()
{

}

//Monta elementos
void SpringDashpot_1::Mount()
{
	for (int ind = 0; ind < 3; ind++)
	{
		xd1(ind, 0) = db.nodes[nodes[0] - 1]->vel[ind];	//Velocidade do Nó 1
		xd2(ind, 0) = db.nodes[nodes[1] - 1]->vel[ind];	//Velocidade do Nó 2
		z1(ind, 0) = db.nodes[nodes[0] - 1]->copy_coordinates[ind] + db.nodes[nodes[0] - 1]->displacements[ind];	//Posição do nó 1
		z2(ind, 0) = db.nodes[nodes[1] - 1]->copy_coordinates[ind] + db.nodes[nodes[1] - 1]->displacements[ind];	//Posição do nó 2
	}
	
	//Calculo da distancia entre os pontos
	z1z2 = z1 - z2;
	//Direção normal
	n = (1.0 / norm(z1z2))*z1z2;
	double norm_z1z2 = norm(z1z2);

	//Calculo do gap - quando não houver inversão de pontos ou for a primeira vez que calcula (ate convergir algum passo do modelo)
	if (dot(last_n, n) >= 0.0 || first_evaluation == true)
	{
		gn = norm_z1z2 - initial_distance;
		f = (gn / norm_z1z2);
	}
		
	else//houve inversão dos pontos
	{
		gn = -(norm_z1z2 + initial_distance);
		n = -1.0*n;//Inversão da normal
		f = gn / (norm_z1z2 + initial_distance);
	}
	
	//Força na mola
	elastic_force = k*gn;
	//Força no amortecedor
	damping_force = c*dot(xd1 - xd2, n);

	//Matrizes auxiliares para operador tangente
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
//Monta matriz de transformação de coordenadas
void SpringDashpot_1::TransformMatrix()
{
	//DOES NOTHING
}
//Zera matrizes locais do elemento
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
//Preenche a contribuição do elemento nas matrizes globais
void SpringDashpot_1::MountGlobal()
{
	int GL_global_1 = 0;
	int GL_global_2 = 0;
	double anterior = 0;
	for (int i = 0; i < 6; i++)
	{
		//Nó 1
		if (i<3)
			GL_global_1 = db.nodes[nodes[0] - 1]->GLs[i];
		//Nó 2
		else
			GL_global_1 = db.nodes[nodes[1] - 1]->GLs[i - 3];

		//Caso o grau de liberdade seja livre:
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
			//Nó 1
			if (j<3)
				GL_global_2 = db.nodes[nodes[0] - 1]->GLs[j];
			//Nó 2
			else
				GL_global_2 = db.nodes[nodes[1] - 1]->GLs[j - 3];

			//Caso os graus de liberdade sejam ambos livres (Matriz Kaa)
			if (GL_global_1 > 0 && GL_global_2 > 0)
				db.global_stiffness_AA.setValue(GL_global_1 - 1, GL_global_2 - 1, c_stiffness(i, j));
			//Caso os graus de liberdade sejam ambos fixos (Matriz Kbb)
			if (GL_global_1 < 0 && GL_global_2 < 0)
				db.global_stiffness_BB.setValue(-GL_global_1 - 1, -GL_global_2 - 1, c_stiffness(i, j));
			//Caso os graus de liberdade sejam livre e fixo (Matriz Kab)
			if (GL_global_1 > 0 && GL_global_2 < 0)
				db.global_stiffness_AB.setValue(GL_global_1 - 1, -GL_global_2 - 1, c_stiffness(i, j));
			//Caso os graus de liberdade sejam fixo e livre (Matriz Kba)
			if (GL_global_1 < 0 && GL_global_2 > 0)
				db.global_stiffness_BA.setValue(-GL_global_1 - 1, GL_global_2 - 1, c_stiffness(i, j));
		}
	}
}
//Salva variaveis nos pontos de Gauss uteis para descrição lagrangiana atualizada
void SpringDashpot_1::SaveLagrange()
{
	//Salva a direção convergida da mola (para checar inversões no passo posterior)
	last_n = n;
	first_evaluation = false;
}
//Pre-calculo de variaveis que e feito uma unica vez no inicio
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

//Monta a matriz de massa
void SpringDashpot_1::MountMass()
{
	//DOES NOTHING
}

//Monta a matriz de massa
void SpringDashpot_1::MountMassModal()
{
	Zeros();
	//DOES NOTHING
}

//Monta a matriz de amortecimento para realização da analise modal
void SpringDashpot_1::MountDampingModal()
{
	Zeros();
	//TODO
}

//Monta a matriz de amortecimento
void SpringDashpot_1::MountDamping(bool update_rayleigh)
{
	//DOES NOTHING
}

//Montagens - Newmark
void SpringDashpot_1::MountDyn()
{
	//Modificações da dinamica
	c_stiffness = c_stiffness + c_damping;
	//Modificações da dinamica nos esforços - presença das forças de amortecimento
	c_stiffness_force = c_stiffness_force + c_damping_force;
}

//Montagens para analise modal - inserção da matriz de massa e amortecimento na matriz de rigidez para posterior montagem global
void SpringDashpot_1::MountDynModal()
{
	c_stiffness = c_damping_modal;
}

