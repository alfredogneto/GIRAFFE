#include "STLSurface.h"
#include <Eigen/dense>

#include "MatrixFloat.h"
#include "TriangularFace.h"
#include "Particle.h"
#include "VEMPolyhedron.h"
#include "Encoding.h"
#include "Database.h"
//Variaveis globais
extern
Database db;

STLSurface::STLSurface()
{
	G = new Matrix(3);
	J_O = new Matrix(3, 3);

	number = 0;
	n_CAD_points = 0;
	coord_double.clear();
	coord_float.clear();
	faces = NULL;
	radius = 0.0;
	volume = 0.0;
	vertice_factors = NULL;
	
	vertices.clear();
	edges.clear();
	tetras.clear();

	mesh_available = false;
}

STLSurface::~STLSurface()
{
	delete G;
	delete J_O;

	coord_double.clear();
	coord_float.clear();
	if (faces != NULL)
	{
		for (int i = 0; i < n_faces; i++)
			delete faces[i];
		delete[]faces;
	}

	vertices.clear();
	edges.clear();
	tetras.clear();
	if (vertice_factors != NULL)
		delete[]vertice_factors;
}

//Leitura do arquivo de entrada
bool STLSurface::Read(FILE *f)
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
void STLSurface::Write(FILE *f)
{
	fprintf(f, "STLSurface\t%d\t%s\n",
		number,
		file
	);
}

//Leitura do arquivo de malha (opcional)
bool STLSurface::ReadMeshFile()
{
	FILE *f1 = NULL;
	char name_file[500];
	strcpy(name_file, db.folder_name);
	strcat(name_file, "CAD/");
	strcat(name_file, file);
	//Correction on file termination - eliminating the ".stl" and replacing for ".msh"
	char *temp;
	temp = strrchr(name_file, '.');   //Get the pointer to the last occurrence to the character '.'
	*temp = '\0';					  //Replace token with null char
	strcat(name_file, ".msh");
	
	//Reading the file
	f1 = fopen(name_file, "r");
	if (f1 == NULL)
	{
		//Not available - return false
		return false;
	}
	else
	{
		char s[1000];

		fscanf(f1, "%s", s);
		if (!strcmp(s, "VERTICES"))
		{
			fscanf(f1, "%s", s);
			while (!strcmp(s, "Vertex"))
			{
				Vertex v1;
				//ID
				fscanf(f1, "%s", s);
				v1.ID = atoi(s) + 1;
				fscanf(f1, "%s", s);
				(*v1.coord_float)(0, 0) = (float)atof(s);
				(*v1.coord_double)(0, 0) = atof(s);
				fscanf(f1, "%s", s);
				(*v1.coord_float)(1, 0) = (float)atof(s);
				(*v1.coord_double)(1, 0) = atof(s);
				fscanf(f1, "%s", s);
				(*v1.coord_float)(2, 0) = (float)atof(s);
				(*v1.coord_double)(2, 0) = atof(s);

				vertices.push_back(v1);
				fscanf(f1, "%s", s);
			}
		}
		else
		{
			db.myprintf("\nError reading mesh file for CAD ID %d\n\n", number);
			fclose(f1);
			return false;
		}
			
		//fscanf(f1, "%s", s);
		if (!strcmp(s, "FACES"))
		{
			fscanf(f1, "%s", s);
			while (!strcmp(s, "Face"))
			{
				//ID
				fscanf(f1, "%s", s);
				int temp_ID = atoi(s) + 1;
				fscanf(f1, "%s", s);
				faces[temp_ID - 1]->verticesIDs[0] = atoi(s) + 1;
				fscanf(f1, "%s", s);
				faces[temp_ID - 1]->verticesIDs[1] = atoi(s) + 1;
				fscanf(f1, "%s", s);
				faces[temp_ID - 1]->verticesIDs[2] = atoi(s) + 1;
				fscanf(f1, "%s", s);
			}
		}
		else
		{
			db.myprintf("\nError reading mesh file for CAD ID %d\n\n", number);
			fclose(f1);
			return false;
		}

		//fscanf(f1, "%s", s);
		if (!strcmp(s, "TETRAHEDRA"))
		{
			fscanf(f1, "%s", s);
			while (!strcmp(s, "Tetrahedron"))
			{
				Tetrahedron t1;
				//ID
				fscanf(f1, "%s", s);
				t1.ID = atoi(s) + 1;
				fscanf(f1, "%s", s);
				t1.verticesIDs[0] = atoi(s) + 1;
				fscanf(f1, "%s", s);
				t1.verticesIDs[1] = atoi(s) + 1;
				fscanf(f1, "%s", s);
				t1.verticesIDs[2] = atoi(s) + 1;
				fscanf(f1, "%s", s);
				t1.verticesIDs[3] = atoi(s) + 1;

				t1.CAD_ID = number;
				tetras.push_back(t1); 
				fscanf(f1, "%s", s);
			}
		}
		else
		{
			db.myprintf("\nError reading mesh file for CAD ID %d\n\n", number);
			fclose(f1);
			return false;
		}

		//Loop on faces to set data
		for (int i = 0; i < n_faces; i++)
		{
			//Creation of 3 edges
			Edge e1;
			Edge e2;
			Edge e3;
			e1.verticesIDs[0] = faces[i]->verticesIDs[0];
			e1.verticesIDs[1] = faces[i]->verticesIDs[1];
			e2.verticesIDs[0] = faces[i]->verticesIDs[1];
			e2.verticesIDs[1] = faces[i]->verticesIDs[2];
			e3.verticesIDs[0] = faces[i]->verticesIDs[2];
			e3.verticesIDs[1] = faces[i]->verticesIDs[0];
			//IDs
			e1.ID = i * 3 + 1;
			e2.ID = i * 3 + 2;
			e3.ID = i * 3 + 3;

			faces[i]->edgesIDs[0] = e1.ID;
			faces[i]->edgesIDs[1] = e2.ID;
			faces[i]->edgesIDs[2] = e3.ID;

			//Adding data to the lists
			edges.push_back(e1);
			edges.push_back(e2);
			edges.push_back(e3);
		}
		fclose(f1);
		return true;
	}
	
}

void STLSurface::PreCalc()
{
	//Set number of faces (triangular only)
	n_faces = n_CAD_points / 3;
	//Alloc faces vector
	faces = new TriangularFace*[n_faces];
	for (int i = 0; i < n_faces; i++)
	{
		faces[i] = new TriangularFace();
		faces[i]->CAD_ID = number;
		faces[i]->ID = i + 1;
	}
	//Leitura do arquivo de malha (opcional)
	if (ReadMeshFile())
		mesh_available = true;
	else
	{
		//Creation of vertices, edges and a connected mesh for the surface
		CreateVerticesEdges();
		MergeVertices();
	}

	MergeEdges();
	OrganizeNumbering();
	SetNumberingInfo();

	//Loop on faces to set data
	for (int i = 0; i < n_faces; i++)
		faces[i]->PreCalc();
	EvaluateRadius();
	EvaluateVolume();
	EvaluateCentroid();
	EvaluateInertiaTensor();
	MarkConcaveEdges();
	PointNormalEdges();
	//Loop on edges to evaluate lengths
	for (int i = 0; i < edges.size(); i++)
	{
		edges[i].CAD_ID = number;
		edges[i].PreCalc();
	}

	//Checks if this surface is employed for some VEMPolyhedron. If this is the case, perform some special functions for generating a volumetric mesh
	bool VEMtrue = false;
	for (int part = 0; part < db.number_particles; part++)
	{
		if (typeid(*db.particles[part]) == typeid(VEMPolyhedron))
		{
			VEMPolyhedron* ptr = static_cast<VEMPolyhedron*>(db.particles[part]);
			if (ptr->CADDATA_ID == number)
				VEMtrue = true;
		}
	}
	if (VEMtrue && mesh_available == false)
	{
		GenerateTetraMesh();
	}
	EvaluateVerticeFactors();
	PrintSurfaceReport();

	////Heavy top
	//volume = 15;
	//(*G)(0, 0) = 0;
	//(*G)(1, 0) = 1.0;
	//(*G)(2, 0) = 0;
	//(*J_O).clear();
	//(*J_O)(0, 0) = 0.234375;
	//(*J_O)(1, 1) = 0.46875;
	//(*J_O)(2, 2) = 0.234375;
	//Matrix I3(3, 3);
	//I3(0, 0) = 1;
	//I3(1, 1) = 1;
	//I3(2, 2) = 1;
	//(*J_O) = (*J_O) + 15 * (dot(*G, *G)*(I3) - dyadic(*G, *G));
}

//Leitura do arquivo de CAD
bool  STLSurface::ReadCADFile()
{
	//Faz leitura do arquivo 'file'
	//atribuir valor para n_CAD_points

	FILE *f1 = NULL;
	char name_file[500];
	strcpy(name_file, db.folder_name);
	strcat(name_file, "CAD/");
	strcat(name_file, file);
	f1 = fopen(name_file, "r");
	if (f1 == NULL)
		return false;//Erro de leitura do arquivo de CAD

	//Leitura do CAD file
	bool scape = false;
	n_CAD_points = 0;
	char s[200];
	coord_float.clear();
	coord_double.clear();
	while (fscanf(f1, "%s", s) != EOF)
	{
		if (!strcmp(s, "vertex"))
		{
			n_CAD_points++;

			fscanf(f1, "%s", s);
			coord_double.push_back(atof(s));
			coord_float.push_back((float)atof(s));
			fscanf(f1, "%s", s);
			coord_double.push_back(atof(s));
			coord_float.push_back((float)atof(s));
			fscanf(f1, "%s", s);
			coord_double.push_back(atof(s));
			coord_float.push_back((float)atof(s));
		}
		if (!strcmp(s, "endsolid"))
			scape = true;
	}

	fclose(f1);
	return true;
}

void STLSurface::WriteVTK_XMLRender(FILE *f, Matrix& pos, Matrix& rot, int number, int matID)
{
	std::vector<int> vint;
	int c2 = 0, c3 = 0, k = 3, j = 0, c4 = 0;
	Matrix rP(3);
	Matrix xP(3);

	//Abre piece
	fprintf(f, "\n<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n<Points>\n<DataArray type = \"Float32\" NumberOfComponents = \"3\" format = \"binary\">\n", n_CAD_points, n_CAD_points / 3);
	std::vector<float> coord2;

	for (int point = 0; point < 3 * n_CAD_points; point += 3)
	{
		rP(0, 0) = coord_float[point];
		rP(1, 0) = coord_float[point + 1];
		rP(2, 0) = coord_float[point + 2];
		//rP = (pegar info do vector) 
		//Calcular a posição atual de cada ponto no sistema global de coordenadas
		xP = pos + rot * rP;
		//		coord2.clear();
		coord2.push_back((float)xP(0, 0));
		coord2.push_back((float)xP(1, 0));
		coord2.push_back((float)xP(2, 0));
	}
	fprintf(f, encodeData<float>(coord2).c_str());
	fprintf(f, "\n</DataArray>\n</Points>");
	fprintf(f, "\n<Cells><DataArray type=\"Int32\" Name=\"connectivity\" format=\"binary\">\n");

	while (c2 < n_CAD_points) {
		vint.push_back(j);
		vint.push_back(j + 1);
		vint.push_back(j + 2);
		c2 = c2 + 3;
		j = j + 3;
	}
	fprintf(f, encodeData<int>(vint).c_str());

	vint.clear();
	fprintf(f, "\n</DataArray>\n<DataArray type=\"Int32\" Name=\"offsets\" format=\"binary\">\n");
	while (c3 < n_CAD_points) {
		vint.push_back(k);
		k = k + 3;
		c3 = c3 + 3;
	}
	fprintf(f, encodeData<int>(vint).c_str());

	vint.clear();
	fprintf(f, "\n</DataArray>\n<DataArray type=\"Int32\" Name=\"types\" format=\"binary\">\n");
	while (c4 < n_CAD_points) {
		vint.push_back(5);
		c4 = c4 + 3;
	}
	fprintf(f, encodeData<int>(vint).c_str());

	vint.clear();
	fprintf(f, "</DataArray>\n</Cells>\n");

	//Opens CellData
	fprintf(f, "\t\t\t<CellData FieldData=\"ElementData\">\n");
	vint.clear();
	//Opens DataArray
	fprintf(f, "\t\t\t\t<DataArray Name=\"GeneralData\" type=\"Int32\" NumberOfComponents=\"%d\" format=\"binary\">\n", 4);
	for (int cell = 0; cell < (n_CAD_points/3); cell++)
	{
		vint.push_back(number);		//Surface ID
		vint.push_back(matID);		//Material ID
		vint.push_back(0);
		vint.push_back(0);
	}
	fprintf(f, encodeData<int>(vint).c_str());
	fprintf(f, "\n");
	//Closes DataArray
	fprintf(f, "\t\t\t\t</DataArray>\n");
	//Closes CellData
	fprintf(f, "\t\t\t</CellData>\n");

	
	fprintf(f, "</Piece>\n");
}

void STLSurface::EvaluateVolume()
{
	volume = 0.0;
	Matrix a(3), b(3), c(3), n(3);
	for (int i = 0; i < 3 * n_CAD_points; i += 9)
	{
		a(0, 0) = coord_double[i];
		a(1, 0) = coord_double[i + 1];
		a(2, 0) = coord_double[i + 2];
		b(0, 0) = coord_double[i + 3];
		b(1, 0) = coord_double[i + 4];
		b(2, 0) = coord_double[i + 5];
		c(0, 0) = coord_double[i + 6];
		c(1, 0) = coord_double[i + 7];
		c(2, 0) = coord_double[i + 8];
		n = cross(b - a, c - a);
		volume += dot(a, n) / 6;
	}
	//printf("Volume: %f\n", volume);
}
void STLSurface::EvaluateCentroid()
{
	Matrix a(3), b(3), ab(3), bc(3), ca(3), c(3), n(3);
	for (int i = 0; i < 3 * n_CAD_points; i += 9)
	{
		a(0, 0) = coord_double[i];
		a(1, 0) = coord_double[i + 1];
		a(2, 0) = coord_double[i + 2];
		b(0, 0) = coord_double[i + 3];
		b(1, 0) = coord_double[i + 4];
		b(2, 0) = coord_double[i + 5];
		c(0, 0) = coord_double[i + 6];
		c(1, 0) = coord_double[i + 7];
		c(2, 0) = coord_double[i + 8];
		n = cross(b - a, c - a);
		ab = a + b;
		bc = b + c;
		ca = c + a;
		for (int j = 0; j < 3; j++)
		{
			(*G)(j, 0) += n(j, 0)*(ab(j, 0)*ab(j, 0) + bc(j, 0)*bc(j, 0) + ca(j, 0)*ca(j, 0)) / (24 * 2 * volume);
		}
	}
}
void STLSurface::EvaluateInertiaTensor()
{
	//TODO
	//Inertia tensor with respect to the origin in local coordinate system
					//      | Jxx -Jxy -Jxz|     //
					//      |-Jxy  Jyy -Jyz|     //
					//      |-Jxz -Jyz  Jzz|     //

	Matrix coord(3);
	double zeta1 = 0.585410196624969, zeta2 = 0.138196601125011;
	double ve = 0.0, w = 0.25;
	for (int i = 0; i < 3 * n_CAD_points; i += 9)
	{
		Eigen::Matrix4d T;

		/*T << 1, 1, 1, 1,
			coord_double[i] - (*G)(0, 0), coord_double[i + 3] - (*G)(0, 0), coord_double[i + 6] - (*G)(0, 0), 0,
			coord_double[i + 1] - (*G)(1, 0), coord_double[i + 4] - (*G)(1, 0), coord_double[i + 7] - (*G)(1, 0), 0,
			coord_double[i + 2] - (*G)(2, 0), coord_double[i + 5] - (*G)(2, 0), coord_double[i + 8] - (*G)(2, 0), 0;*/
		T << 1, 1, 1, 1,
			coord_double[i], coord_double[i + 3], coord_double[i + 6], (*G)(0, 0),
			coord_double[i + 1], coord_double[i + 4], coord_double[i + 7], (*G)(1, 0),
			coord_double[i + 2], coord_double[i + 5], coord_double[i + 8], (*G)(2, 0);
		ve = (-1)*T.determinant() / 6;
		for (int j = 1; j < 4; j++)
		{
			coord(j - 1, 0) = zeta1 * T(j, 0) + zeta2 * T(j, 1) + zeta2 * T(j, 2) + zeta2 * T(j, 3);
		}
		(*J_O)(0, 0) += ve * w * (coord(1, 0) * coord(1, 0) + coord(2, 0) * coord(2, 0));
		(*J_O)(1, 1) += ve * w * (coord(0, 0) * coord(0, 0) + coord(2, 0) * coord(2, 0));
		(*J_O)(2, 2) += ve * w * (coord(0, 0) * coord(0, 0) + coord(1, 0) * coord(1, 0));
		(*J_O)(1, 0) += ve * w * (coord(0, 0) * coord(1, 0));
		(*J_O)(2, 0) += ve * w * (coord(0, 0) * coord(2, 0));
		(*J_O)(2, 1) += ve * w * (coord(1, 0) * coord(2, 0));
		for (int j = 1; j < 4; j++)
		{
			coord(j - 1, 0) = zeta2 * T(j, 0) + zeta1 * T(j, 1) + zeta2 * T(j, 2) + zeta2 * T(j, 3);
		}
		(*J_O)(0, 0) += ve * w * (coord(1, 0) * coord(1, 0) + coord(2, 0) * coord(2, 0));
		(*J_O)(1, 1) += ve * w * (coord(0, 0) * coord(0, 0) + coord(2, 0) * coord(2, 0));
		(*J_O)(2, 2) += ve * w * (coord(0, 0) * coord(0, 0) + coord(1, 0) * coord(1, 0));
		(*J_O)(1, 0) += ve * w * (coord(0, 0) * coord(1, 0));
		(*J_O)(2, 0) += ve * w * (coord(0, 0) * coord(2, 0));
		(*J_O)(2, 1) += ve * w * (coord(1, 0) * coord(2, 0));
		for (int j = 1; j < 4; j++)
		{
			coord(j - 1, 0) = zeta2 * T(j, 0) + zeta2 * T(j, 1) + zeta1 * T(j, 2) + zeta2 * T(j, 3);
		}
		(*J_O)(0, 0) += ve * w * (coord(1, 0) * coord(1, 0) + coord(2, 0) * coord(2, 0));
		(*J_O)(1, 1) += ve * w * (coord(0, 0) * coord(0, 0) + coord(2, 0) * coord(2, 0));
		(*J_O)(2, 2) += ve * w * (coord(0, 0) * coord(0, 0) + coord(1, 0) * coord(1, 0));
		(*J_O)(1, 0) += ve * w * (coord(0, 0) * coord(1, 0));
		(*J_O)(2, 0) += ve * w * (coord(0, 0) * coord(2, 0));
		(*J_O)(2, 1) += ve * w * (coord(1, 0) * coord(2, 0));
		for (int j = 1; j < 4; j++)
		{
			coord(j - 1, 0) = zeta2 * T(j, 0) + zeta2 * T(j, 1) + zeta2 * T(j, 2) + zeta1 * T(j, 3);
		}
		(*J_O)(0, 0) += ve * w * (coord(1, 0) * coord(1, 0) + coord(2, 0) * coord(2, 0));
		(*J_O)(1, 1) += ve * w * (coord(0, 0) * coord(0, 0) + coord(2, 0) * coord(2, 0));
		(*J_O)(2, 2) += ve * w * (coord(0, 0) * coord(0, 0) + coord(1, 0) * coord(1, 0));
		(*J_O)(1, 0) += ve * w * (coord(0, 0) * coord(1, 0));
		(*J_O)(2, 0) += ve * w * (coord(0, 0) * coord(2, 0));
		(*J_O)(2, 1) += ve * w * (coord(1, 0) * coord(2, 0));
	}
	(*J_O)(1, 0) = (-1)*(*J_O)(1, 0);
	(*J_O)(0, 1) = (*J_O)(1, 0);
	(*J_O)(2, 0) = (-1)*(*J_O)(2, 0);
	(*J_O)(0, 2) = (*J_O)(2, 0);
	(*J_O)(2, 1) = (-1)*(*J_O)(2, 1);
	(*J_O)(1, 2) = (*J_O)(2, 1);
	//J_O->print();
}
void STLSurface::EvaluateRadius()
{
	//Evaluates the particle radius (w/r to the origin)
	//Loop on vertices
	radius = 0.0;
	float temp_len = 0.0;
	for (int i = 0; i < 3 * n_CAD_points; i += 3)
	{
		temp_len = sqrt(coord_float[i] * coord_float[i] + coord_float[i + 1] * coord_float[i + 1] + coord_float[i + 2] * coord_float[i + 2]);
		if (temp_len > radius)
			radius = temp_len;
	}
}

void STLSurface::CreateVerticesEdges()
{
	//Loop on faces to set data
	for (int i = 0; i < n_faces; i++)
	{
		//Creation of 3 vertices
		Vertex v1;
		Vertex v2;
		Vertex v3;
		//IDs
		v1.ID = i * 3 + 1;
		v2.ID = i * 3 + 2;
		v3.ID = i * 3 + 3;
		//Creation of 3 edges
		Edge e1;
		Edge e2;
		Edge e3;
		e1.verticesIDs[0] = v1.ID;
		e1.verticesIDs[1] = v2.ID;
		e2.verticesIDs[0] = v2.ID;
		e2.verticesIDs[1] = v3.ID;
		e3.verticesIDs[0] = v3.ID;
		e3.verticesIDs[1] = v1.ID;
		//IDs
		e1.ID = i * 3 + 1;
		e2.ID = i * 3 + 2;
		e3.ID = i * 3 + 3;
		
		//Loop on x,y,z
		for (int j = 0; j < 3; j++)
		{
			(*v1.coord_float)(j, 0) = coord_float[i * 9 + j];
			(*v2.coord_float)(j, 0) = coord_float[i * 9 + j + 3];
			(*v3.coord_float)(j, 0) = coord_float[i * 9 + j + 6];

			(*v1.coord_double)(j, 0) = coord_double[i * 9 + j];
			(*v2.coord_double)(j, 0) = coord_double[i * 9 + j + 3];
			(*v3.coord_double)(j, 0) = coord_double[i * 9 + j + 6];
		}
		
		//Setting faces data
		faces[i]->verticesIDs[0] = v1.ID;
		faces[i]->verticesIDs[1] = v2.ID;
		faces[i]->verticesIDs[2] = v3.ID;
		faces[i]->edgesIDs[0] = e1.ID;
		faces[i]->edgesIDs[1] = e2.ID;
		faces[i]->edgesIDs[2] = e3.ID;

		//Adding data to the lists
		vertices.push_back(v1);
		vertices.push_back(v2);
		vertices.push_back(v3);
		edges.push_back(e1);
		edges.push_back(e2);
		edges.push_back(e3);
	}
}

void STLSurface::MergeVertices()
{
	//Searching for vertices with coincident location
	int i = 0;
	int j = i + 1;
	while (i < vertices.size())
	{
		j = i + 1;
		while (j < vertices.size())
		{
			if (vertices[i] == vertices[j])
			{
				ReplaceVertexID(j, vertices[i].ID);
				vertices.erase(vertices.begin() + j);
			}
			else
				j++;
		}
		i++;
	}
}

void STLSurface::MergeEdges()
{
	//Searching for equal edges
	int i = 0;
	int j = i + 1;
	while (i < edges.size())
	{
		j = i + 1;
		while (j < edges.size())
		{
			if (edges[i] == edges[j])
			{
				ReplaceEdgeID(j, edges[i].ID);
				edges.erase(edges.begin() + j);
			}
			else
				j++;
		}
		i++;
	}
}

void STLSurface::OrganizeNumbering()
{
	//Re-organize vertices and edges IDs, by imposing the same ID as the position of the vector. This is done for further avoiding searching on the vector for a given ID
	for (int ver = 0; ver < vertices.size(); ver++)
	{
		int cur_max_ID = n_faces * 3 + 1;
		//Looking for previous ID equal to the new one to be set
		for (int ver2 = 0; ver2 < vertices.size(); ver2++)
		{
			if (vertices[ver2].ID == ver + 1)
			{
				ReplaceVertexID(ver2, cur_max_ID);
				cur_max_ID++;
			}
		}
		ReplaceVertexID(ver, ver + 1);
	}
	for (int ed = 0; ed < edges.size(); ed++)
	{
		int cur_max_ID = n_faces * 3 + 1;
		//Looking for previous ID equal to the new one to be set
		for (int ed2 = 0; ed2 < edges.size(); ed2++)
		{
			if (edges[ed2].ID == ed + 1)
			{
				ReplaceEdgeID(ed2, cur_max_ID);
				cur_max_ID++;
			}
		}
		ReplaceEdgeID(ed, ed + 1);
	}
	vertices.shrink_to_fit();
	edges.shrink_to_fit();
}

void STLSurface::ReplaceVertexID(int index, int newID)
{
	int oldID = vertices[index].ID;
	//Replacing vertices[j] for vertices[i] on edges
	for (int ed = 0; ed < edges.size(); ed++)
	{
		if (edges[ed].verticesIDs[0] == oldID)
			edges[ed].verticesIDs[0] = newID;
		if (edges[ed].verticesIDs[1] == oldID)
			edges[ed].verticesIDs[1] = newID;
	}
	//Replacing vertices[j] for vertices[i] on faces
	for (int fac = 0; fac < n_faces; fac++)
	{
		if (faces[fac]->verticesIDs[0] == oldID)
			faces[fac]->verticesIDs[0] = newID;
		if (faces[fac]->verticesIDs[1] == oldID)
			faces[fac]->verticesIDs[1] = newID;
		if (faces[fac]->verticesIDs[2] == oldID)
			faces[fac]->verticesIDs[2] = newID;
	}
	vertices[index].ID = newID;
}

void STLSurface::ReplaceEdgeID(int index, int newID)
{
	int oldID = edges[index].ID;

	//Replacing edges[j] for edges[i] on faces
	for (int fac = 0; fac < n_faces; fac++)
	{
		if (faces[fac]->edgesIDs[0] == oldID)
			faces[fac]->edgesIDs[0] = newID;
		if (faces[fac]->edgesIDs[1] == oldID)
			faces[fac]->edgesIDs[1] = newID;
		if (faces[fac]->edgesIDs[2] == oldID)
			faces[fac]->edgesIDs[2] = newID;
	}
	edges[index].ID = newID;
}

void STLSurface::SetNumberingInfo()
{
	//For each vertex, saves to which faces it is connected to and which edges it is connected to
	for (int i = 0; i < vertices.size(); i++)
	{
		//Loops on faces
		for (int fac = 0; fac < n_faces; fac++)
		{
			//If the face contains that vertex
			if (faces[fac]->verticesIDs[0] == vertices[i].ID || faces[fac]->verticesIDs[1] == vertices[i].ID || faces[fac]->verticesIDs[2] == vertices[i].ID)
			{
				bool already_set = false;
				for (int ii = 0; ii < vertices[i].faceIDs.size(); ii++)
				{
					if (faces[fac]->ID == vertices[i].faceIDs[ii])
						already_set = true;
				}
				if (already_set == false)
					vertices[i].faceIDs.push_back(faces[fac]->ID);
			}
		}
		//Loops on edges
		for (int ed = 0; ed < edges.size(); ed++)
		{
			//If the edge contains that vertex
			if (edges[ed].verticesIDs[0] == vertices[i].ID || edges[ed].verticesIDs[1] == vertices[i].ID)
			{
				bool already_set = false;
				for (int ii = 0; ii < vertices[i].edgeIDs.size(); ii++)
				{
					if (edges[ed].ID == vertices[i].edgeIDs[ii])
						already_set = true;
				}
				if (already_set == false)
					vertices[i].edgeIDs.push_back(edges[ed].ID);
			}
		}
	}

	//For each edge, saves to which faces it is connected to
	for (int i = 0; i < edges.size(); i++)
	{
		//Loops on faces
		for (int fac = 0; fac < n_faces; fac++)
		{
			//If the face contains that edge
			if (faces[fac]->edgesIDs[0] == edges[i].ID || faces[fac]->edgesIDs[1] == edges[i].ID || faces[fac]->edgesIDs[2] == edges[i].ID)
			{
				bool already_set = false;
				for (int ii = 0; ii < edges[i].faceIDs.size(); ii++)
				{
					if (faces[fac]->ID == edges[i].faceIDs[ii])
						already_set = true;
				}
				if (already_set == false)
					edges[i].faceIDs.push_back(faces[fac]->ID);
			}
		}
	}
}

void STLSurface::WriteVTK_XMLMesh(MatrixFloat& pos, MatrixFloat& rot)
{
	std::vector<int> vint;
	MatrixFloat rP(3);
	MatrixFloat xP(3);

	FILE *f;
	char name[1000];
	strcpy(name, db.folder_name);
	strcat(name, "CAD/");
	strcat(name, file);
	//Correction on file termination - eliminating the ".stl" and replacing for "_report.txt"
	char *temp;
	temp = strrchr(name, '.');   //Get the pointer to the last occurrence to the character '.'
	*temp = '\0';					  //Replace token with null char
	strcat(name, "_vertices.vtu");
	f = fopen(name, "w");

	//Cabeçalho do arquivo VTK
	fprintf(f, "<?xml version=\"1.0\"?>\n");
	fprintf(f, "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n");
	fprintf(f, "<!--INPUT: %s-->\n", db.file_name);
	fprintf(f, "\t<UnstructuredGrid>\n");

	//Abre piece
	fprintf(f, "\n<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n<Points>\n<DataArray type = \"Float32\" NumberOfComponents = \"3\" format = \"binary\">\n", (int)vertices.size(), (int)vertices.size());
	std::vector<float> coord2;

	for (int i = 0; i < vertices.size(); i++)
	{
		rP(0, 0) = (*vertices[i].coord_float)(0, 0);
		rP(1, 0) = (*vertices[i].coord_float)(1, 0);
		rP(2, 0) = (*vertices[i].coord_float)(2, 0);
		//Calcular a posição atual de cada ponto no sistema global de coordenadas
		xP = pos + rot * rP;
		//		coord2.clear();
		coord2.push_back(xP(0, 0));
		coord2.push_back(xP(1, 0));
		coord2.push_back(xP(2, 0));
	}

	fprintf(f, encodeData<float>(coord2).c_str());
	fprintf(f, "\n</DataArray>\n</Points>");

	fprintf(f, "\n<Cells><DataArray type=\"Int32\" Name=\"connectivity\" format=\"binary\">\n");

	//1 - Vertices
	for (int i = 0; i < vertices.size(); i++)
		vint.push_back(vertices[i].ID - 1);
	
	fprintf(f, encodeData<int>(vint).c_str());

	vint.clear();
	fprintf(f, "\n</DataArray>\n<DataArray type=\"Int32\" Name=\"offsets\" format=\"binary\">\n");
	int cur_offset = 0;
	for (int i = 0; i < vertices.size(); i++)
	{
		cur_offset = cur_offset + 1;
		vint.push_back(cur_offset);
	}
	fprintf(f, encodeData<int>(vint).c_str());

	vint.clear();
	fprintf(f, "\n</DataArray>\n<DataArray type=\"Int32\" Name=\"types\" format=\"binary\">\n");
	for (int i = 0; i < vertices.size(); i++)
	{
		vint.push_back(1);
	}
	fprintf(f, encodeData<int>(vint).c_str());

	vint.clear();
	fprintf(f, "</DataArray>\n</Cells>\n");
	
	fprintf(f, "</Piece>\n");

	fprintf(f, "\t</UnstructuredGrid>\n");
	fprintf(f, "</VTKFile>\n");
	fclose(f);

	//Edges
	strcpy(name, db.folder_name);
	strcat(name, "CAD/");
	strcat(name, file);
	//Correction on file termination - eliminating the ".stl" and replacing for "_report.txt"
	temp = strrchr(name, '.');   //Get the pointer to the last occurrence to the character '.'
	*temp = '\0';					  //Replace token with null char
	strcat(name, "_edges.vtu");
	f = fopen(name, "w");

	//Cabeçalho do arquivo VTK
	fprintf(f, "<?xml version=\"1.0\"?>\n");
	fprintf(f, "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n");
	fprintf(f, "<!--INPUT: %s-->\n", db.file_name);
	fprintf(f, "\t<UnstructuredGrid>\n");

	//Abre piece
	fprintf(f, "\n<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n<Points>\n<DataArray type = \"Float32\" NumberOfComponents = \"3\" format = \"binary\">\n", (int)vertices.size(), (int)edges.size());
	coord2.clear();

	for (int i = 0; i < vertices.size(); i++)
	{
		rP(0, 0) = (*vertices[i].coord_float)(0, 0);
		rP(1, 0) = (*vertices[i].coord_float)(1, 0);
		rP(2, 0) = (*vertices[i].coord_float)(2, 0);
		//Calcular a posição atual de cada ponto no sistema global de coordenadas
		xP = pos + rot * rP;
		//		coord2.clear();
		coord2.push_back(xP(0, 0));
		coord2.push_back(xP(1, 0));
		coord2.push_back(xP(2, 0));
	}

	fprintf(f, encodeData<float>(coord2).c_str());
	fprintf(f, "\n</DataArray>\n</Points>");

	fprintf(f, "\n<Cells><DataArray type=\"Int32\" Name=\"connectivity\" format=\"binary\">\n");

	//2 - Edges
	for (int i = 0; i < edges.size(); i++)
	{
		vint.push_back(edges[i].verticesIDs[0] - 1);
		vint.push_back(edges[i].verticesIDs[1] - 1);
	}
	
	fprintf(f, encodeData<int>(vint).c_str());

	vint.clear();
	fprintf(f, "\n</DataArray>\n<DataArray type=\"Int32\" Name=\"offsets\" format=\"binary\">\n");
	cur_offset = 0;
	
	for (int i = 0; i < edges.size(); i++)
	{
		cur_offset = cur_offset + 2;
		vint.push_back(cur_offset);
	}
	
	fprintf(f, encodeData<int>(vint).c_str());

	vint.clear();
	fprintf(f, "\n</DataArray>\n<DataArray type=\"Int32\" Name=\"types\" format=\"binary\">\n");
	
	for (int i = 0; i < edges.size(); i++)
	{
		vint.push_back(3);
	}
	
	fprintf(f, encodeData<int>(vint).c_str());

	vint.clear();
	fprintf(f, "</DataArray>\n</Cells>\n");

	fprintf(f, "</Piece>\n");

	fprintf(f, "\t</UnstructuredGrid>\n");
	fprintf(f, "</VTKFile>\n");
	fclose(f);


	
	//Edges
	strcpy(name, db.folder_name);
	strcat(name, "CAD/");
	strcat(name, file);
	//Correction on file termination - eliminating the ".stl" and replacing for "_report.txt"
	temp = strrchr(name, '.');   //Get the pointer to the last occurrence to the character '.'
	*temp = '\0';					  //Replace token with null char
	strcat(name, "_faces.vtu");
	f = fopen(name, "w");

	//Cabeçalho do arquivo VTK
	fprintf(f, "<?xml version=\"1.0\"?>\n");
	fprintf(f, "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n");
	fprintf(f, "<!--INPUT: %s-->\n", db.file_name);
	fprintf(f, "\t<UnstructuredGrid>\n");

	//Abre piece
	fprintf(f, "\n<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n<Points>\n<DataArray type = \"Float32\" NumberOfComponents = \"3\" format = \"binary\">\n", (int)vertices.size(), n_faces);
	coord2.clear();

	for (int i = 0; i < vertices.size(); i++)
	{
		rP(0, 0) = (*vertices[i].coord_float)(0, 0);
		rP(1, 0) = (*vertices[i].coord_float)(1, 0);
		rP(2, 0) = (*vertices[i].coord_float)(2, 0);
		//Calcular a posição atual de cada ponto no sistema global de coordenadas
		xP = pos + rot * rP;
		//		coord2.clear();
		coord2.push_back(xP(0, 0));
		coord2.push_back(xP(1, 0));
		coord2.push_back(xP(2, 0));
	}

	fprintf(f, encodeData<float>(coord2).c_str());
	fprintf(f, "\n</DataArray>\n</Points>");

	fprintf(f, "\n<Cells><DataArray type=\"Int32\" Name=\"connectivity\" format=\"binary\">\n");

	//3 - Faces
	for (int i = 0; i < n_faces; i++)
	{
		vint.push_back(faces[i]->verticesIDs[0] - 1);
		vint.push_back(faces[i]->verticesIDs[1] - 1);
		vint.push_back(faces[i]->verticesIDs[2] - 1);
	}

	fprintf(f, encodeData<int>(vint).c_str());

	vint.clear();
	fprintf(f, "\n</DataArray>\n<DataArray type=\"Int32\" Name=\"offsets\" format=\"binary\">\n");
	cur_offset = 0;
	for (int i = 0; i < n_faces; i++)
	{
		cur_offset = cur_offset + 3;
		vint.push_back(cur_offset);
	}
	fprintf(f, encodeData<int>(vint).c_str());

	vint.clear();
	fprintf(f, "\n</DataArray>\n<DataArray type=\"Int32\" Name=\"types\" format=\"binary\">\n");
	for (int i = 0; i < n_faces; i++)
	{
		vint.push_back(5);
	}
	fprintf(f, encodeData<int>(vint).c_str());

	vint.clear();
	fprintf(f, "</DataArray>\n</Cells>\n");

	fprintf(f, "</Piece>\n");

	fprintf(f, "\t</UnstructuredGrid>\n");
	fprintf(f, "</VTKFile>\n");
	fclose(f);
}

void STLSurface::MarkConcaveEdges()
{
	//This function evaluates for each edge of the surface if it is a concave or a convex edge (for closed surfaces)
	//When the edge is no shared between faces, it is set as convex (concave_indicator = false)
	float epsilon = 1e-6f;
	for (int i = 0; i < edges.size(); i++)
	{
		if (edges[i].faceIDs.size() == 1)
			edges[i].concave_indicator = -1; //convex
		else
		{
			//face 1
			MatrixFloat C_face1 = *faces[edges[i].faceIDs[0] - 1]->centroid;
			//face 2
			MatrixFloat C_face2 = *faces[edges[i].faceIDs[1] - 1]->centroid;

			int v1ID_face2 = faces[edges[i].faceIDs[1] - 1]->verticesIDs[0];
			int v2ID_face2 = faces[edges[i].faceIDs[1] - 1]->verticesIDs[1];
			int v3ID_face2 = faces[edges[i].faceIDs[1] - 1]->verticesIDs[2];

			MatrixFloat t1_face2 = *vertices[v2ID_face2 - 1].coord_float - *vertices[v1ID_face2 - 1].coord_float;
			MatrixFloat t2_face2 = *vertices[v3ID_face2 - 1].coord_float - *vertices[v1ID_face2 - 1].coord_float;
			MatrixFloat n_face2 = (float)1.0 / norm(cross(t1_face2, t2_face2))*cross(t1_face2, t2_face2);
			
			MatrixFloat C1C2 = ((float)1.0/norm(C_face1 - C_face2))*(C_face1 - C_face2);

			//Concave or convex evaluation
			if (dot(C1C2, n_face2) < -epsilon)
				edges[i].concave_indicator = -1; //convex
			else
			{
				if (dot(C1C2, n_face2) > +epsilon)
					edges[i].concave_indicator = +1; //concave
				else
					edges[i].concave_indicator = 0; //flat
			}
				
		}
	}

	//Set has_concave_edge variables on faces
	for (int i = 0; i < n_faces; i++)
	{
		if (edges[faces[i]->edgesIDs[0] - 1].concave_indicator > 0)
			faces[i]->has_concave_edge = true;
		if (edges[faces[i]->edgesIDs[1] - 1].concave_indicator > 0)
			faces[i]->has_concave_edge = true;
		if (edges[faces[i]->edgesIDs[2] - 1].concave_indicator > 0)
			faces[i]->has_concave_edge = true;
	}
}

void STLSurface::PointNormalEdges()
{
	//This function associates pointers of edge normals with face normals
	for (int i = 0; i < edges.size(); i++)
	{
		if (edges[i].faceIDs.size() == 1)
			edges[i].n2 = NULL;
		else
		{
			//face 1
			edges[i].n1 = faces[edges[i].faceIDs[0] - 1]->normal;
			//face 2
			edges[i].n2 = faces[edges[i].faceIDs[1] - 1]->normal;
		}
	}
}

void STLSurface::PrintSurfaceReport()
{
	//PARAVIEW FILE
	MatrixFloat pos(3);
	MatrixFloat rot(3, 3);
	rot(0, 0) = 1.0;
	rot(1, 1) = 1.0;
	rot(2, 2) = 1.0;
	WriteVTK_XMLMesh(pos, rot);


	//REPORT
	FILE *f;
	char name[1000];
	strcpy(name, db.folder_name);
	strcat(name, "CAD/");
	strcat(name, file);
	//Correction on file termination - eliminating the ".stl" and replacing for "_report.txt"
	char *temp;
	temp = strrchr(name, '.');   //Get the pointer to the last occurrence to the character '.'
	*temp = '\0';					  //Replace token with null char
	strcat(name, "_report.txt");

	//Escrevendo o nome do arquivo
	f = fopen(name, "w");
	fprintf(f, "IDs on this report are zero-based to match the Paraview IDs\n\n");
	fprintf(f, "FACES\n");
	for (int i = 0; i < n_faces; i++)
		faces[i]->Print(f);
	fprintf(f, "EDGES\n");
	for (int i = 0; i < edges.size(); i++)
		edges[i].Print(f);
	fprintf(f, "VERTICES\n");
	for (int i = 0; i < vertices.size(); i++)
		vertices[i].Print(f);
	fprintf(f, "TETRAHEDRA\n");
	for (int i = 0; i < tetras.size(); i++)
		tetras[i].Print(f);
	
	fprintf(f, "\n\nInertia properties are based only on geometry\n\n");
	fprintf(f, "Volume\n");
	fprintf(f, "%.6e\n",volume);
	fprintf(f, "Surface Area\n");
	fprintf(f, "%.6e\n", total_ref_area);
	fprintf(f, "Barycenter\n");
	fprintf(f, "(%.6e, %.6e, %.6e)\n", (*G)(0, 0), (*G)(1, 0), (*G)(2, 0));
	fprintf(f, "Inertia matrix (evaluated as a rigid volume)\n");
	fprintf(f, "|%.6e\t%.6e\t%.6e|\n", (*J_O)(0, 0), (*J_O)(0, 1), (*J_O)(0, 2));
	fprintf(f, "|%.6e\t%.6e\t%.6e|\n", (*J_O)(1, 0), (*J_O)(1, 1), (*J_O)(1, 2));
	fprintf(f, "|%.6e\t%.6e\t%.6e|\n", (*J_O)(2, 0), (*J_O)(2, 1), (*J_O)(2, 2));

	/*fprintf(f, "Face Areas\n");
	for (int i=0;i<n_faces;i++)
		fprintf(f, "Face\t%d\t%.6e\n",faces[i]->ID-1,faces[i]->area);*/

	/*fprintf(f, "Inertia matrix (evaluated as a set of lumped volumes)\n");
	double jxx = 0.0;
	double jxy = 0.0;
	double jxz = 0.0;
	double jyy = 0.0;
	double jyz = 0.0;
	double jzz = 0.0;
	double vol2 = 0.0;
	for (int i = 0; i < (int)vertices.size(); i++)
	{
		jxx += vertice_factors[i] * volume*((*vertices[i].coord_double)(1, 0)*(*vertices[i].coord_double)(1, 0) + 
			(*vertices[i].coord_double)(2, 0)*(*vertices[i].coord_double)(2, 0));
		jyy += vertice_factors[i] * volume*((*vertices[i].coord_double)(0, 0)*(*vertices[i].coord_double)(0, 0) +
			(*vertices[i].coord_double)(2, 0)*(*vertices[i].coord_double)(2, 0));
		jzz += vertice_factors[i] * volume*((*vertices[i].coord_double)(1, 0)*(*vertices[i].coord_double)(1, 0) +
			(*vertices[i].coord_double)(0, 0)*(*vertices[i].coord_double)(0, 0));
		jxy += vertice_factors[i] * volume*(*vertices[i].coord_double)(0, 0)*(*vertices[i].coord_double)(1, 0);
		jxz += vertice_factors[i] * volume*(*vertices[i].coord_double)(0, 0)*(*vertices[i].coord_double)(2, 0);
		jyz += vertice_factors[i] * volume*(*vertices[i].coord_double)(1, 0)*(*vertices[i].coord_double)(2, 0);
		vol2 += vertice_factors[i] * volume;
	}
	fprintf(f, "|%.6e\t%.6e\t%.6e|\n", jxx, -jxy, -jxz);
	fprintf(f, "|%.6e\t%.6e\t%.6e|\n", -jxy, jyy, -jyz);
	fprintf(f, "|%.6e\t%.6e\t%.6e|\n", -jxz, -jyz, jzz);
	fprintf(f, "Volume (evaluated as a set of lumped volumes)\n");
	fprintf(f, "%.6e\n", vol2);
	fprintf(f, "Vertice factors\n");
	double sum = 0.0;
	for (int i = 0; i < (int)vertices.size(); i++)
	{
		sum += vertice_factors[i];
		fprintf(f, "%.6e\n", vertice_factors[i]);
	}
	fprintf(f, "Sum: %.6e\n", sum);*/
	fclose(f);
}

int STLSurface::GetVertexAssociatedwithBothEdges(int edge1, int edge2)
{
	int v1 = edges[edge1 - 1].verticesIDs[0];
	int v2 = edges[edge1 - 1].verticesIDs[1];
	int v3 = edges[edge2 - 1].verticesIDs[0];
	int v4 = edges[edge2 - 1].verticesIDs[1];

	if (v1 == v3)
		return v1;
	if (v1 == v4)
		return v1;
	if (v2 == v3)
		return v2;
	if (v2 == v4)
		return v2;
	//If no return yet called here, no common vertex exists
	return 0;
}

void STLSurface::GenerateTetraMesh()
{
	db.myprintf("No required mesh information of the CAD was provided\n");
	//Generates the mesh of tetrahedrons on the interior of the body (if this is desired)
	//Teste
	/*Tetrahedron t1;
	t1.CAD_ID = number;
	t1.ID = 1;
	t1.verticesIDs[0] = 1;
	t1.verticesIDs[1] = 2;
	t1.verticesIDs[2] = 3;
	t1.verticesIDs[3] = 4;
	tetras.push_back(t1);*/
	//Creation of 3 vertices
	/*Tetrahedron t1;
	t1.CAD_ID = number;
	t1.ID = 1;
	t1.verticesIDs[0] = 8;
	t1.verticesIDs[1] = 1;
	t1.verticesIDs[2] = 2;
	t1.verticesIDs[3] = 5;
	tetras.push_back(t1);
	t1.CAD_ID = number;
	t1.ID = 2;
	t1.verticesIDs[0] = 1;
	t1.verticesIDs[1] = 4;
	t1.verticesIDs[2] = 8;
	t1.verticesIDs[3] = 2;
	tetras.push_back(t1);
	t1.CAD_ID = number;
	t1.ID = 3;
	t1.verticesIDs[0] = 4;
	t1.verticesIDs[1] = 3;
	t1.verticesIDs[2] = 1;
	t1.verticesIDs[3] = 8;
	tetras.push_back(t1);
	t1.CAD_ID = number;
	t1.ID = 4;
	t1.verticesIDs[0] = 1;
	t1.verticesIDs[1] = 3;
	t1.verticesIDs[2] = 7;
	t1.verticesIDs[3] = 8;
	tetras.push_back(t1);
	t1.CAD_ID = number;
	t1.ID = 5;
	t1.verticesIDs[0] = 8;
	t1.verticesIDs[1] = 7;
	t1.verticesIDs[2] = 1;
	t1.verticesIDs[3] = 5;
	tetras.push_back(t1);
	t1.CAD_ID = number;
	t1.ID = 6;
	t1.verticesIDs[0] = 1;
	t1.verticesIDs[1] = 7;
	t1.verticesIDs[2] = 6;
	t1.verticesIDs[3] = 5;
	tetras.push_back(t1);*/
	//{8, 1, 2, 5}, {1, 4, 8, 2}, {4, 3, 1, 8}, {1, 3, 7, 8}, {8, 7, 1, 5}, {1, 7, 6, 5}
}

void STLSurface::EvaluateVerticeFactors()
{
	vertice_factors = new double[(int)vertices.size()];
	//Surface area evaluation
	total_ref_area = 0.0;
	for (int i = 0; i < n_faces; i++)
		total_ref_area += faces[i]->area;

	double temp_w = 0.0;
	for (int i = 0; i < (int)vertices.size(); i++)
	{
		temp_w = 0.0;
		//Loops on faces of the vertex under investigation
		for (int f = 0; f < (int)vertices[i].faceIDs.size(); f++)
		{
			int temp_face = vertices[i].faceIDs[f];
			//For each face, attributes the area contribution to the vertex
			if (faces[temp_face - 1]->verticesIDs[0] == vertices[i].ID)
				temp_w += faces[temp_face - 1]->vertice_weight[0] * faces[temp_face - 1]->area;
			if (faces[temp_face - 1]->verticesIDs[1] == vertices[i].ID)
				temp_w += faces[temp_face - 1]->vertice_weight[1] * faces[temp_face - 1]->area;
			if (faces[temp_face - 1]->verticesIDs[2] == vertices[i].ID)
				temp_w += faces[temp_face - 1]->vertice_weight[2] * faces[temp_face - 1]->area;

		}
		temp_w = temp_w / total_ref_area;
		vertice_factors[i] = temp_w;
	}
	/*for (int i = 0; i < (int)vertices.size(); i++)
		vertice_factors[i] = 0.0;
	double sum = 0.0;
	for (int i = 0; i < (int)n_faces; i++)
	{
		sum = 0.0;
		vertice_factors[faces[i]->verticesIDs[0] - 1] += faces[i]->vertice_weight[0] * faces[i]->area;
		vertice_factors[faces[i]->verticesIDs[1] - 1] += faces[i]->vertice_weight[1] * faces[i]->area;
		vertice_factors[faces[i]->verticesIDs[2] - 1] += faces[i]->vertice_weight[2] * faces[i]->area;
		sum = faces[i]->vertice_weight[0] + faces[i]->vertice_weight[1] + faces[i]->vertice_weight[2];
		printf("%lf\n", sum);
	}
	for (int i = 0; i < (int)vertices.size(); i++)
		vertice_factors[i] = vertice_factors[i] / total_ref_area;*/
}