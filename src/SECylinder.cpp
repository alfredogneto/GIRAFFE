#include "SECylinder.h"

#include "BoundingCylinder.h"
#include "MatrixFloat.h"
#include "Node.h"
#include "CoordinateSystem.h"
#include "GeneralContactSearch.h"
#include "Interface_1.h"
#include "PostFiles.h"
#include "Encoding.h"

#include"Database.h"

//extern
//FILE *fdebug;

//Variáveis globais
extern
Database db;

SECylinder::SECylinder()
{
	super_node = 0;

	nDOFs = 12;
	n_nodes = 2;

	Alloc();
	AllocSpecific();

	//Rotina para ativar os GLS de cada nó
	for (int i = 0; i < n_nodes; i++)
	{
		for (int j = 0; j < db.number_GLs_node; j++)
		{
			DOFs[i][j] = 0;
		}
		//Translação
		DOFs[i][0] = 1;
		DOFs[i][1] = 1;
		DOFs[i][2] = 1;
		//Rotação
		DOFs[i][3] = 1;
		DOFs[i][4] = 1;
		DOFs[i][5] = 1;
	}

	rad = 0;
}

void SECylinder::AllocSpecific()
{
	bv = new BoundingCylinder();
	sprintf(type_name, "SECylinder");
	bv_offset = 0.0f;

	//Variáveis para calcular estado atual (nas funcoes de bounding volume)
	xAf = new MatrixFloat(3);
	xBf = new MatrixFloat(3);

	a = new double;
	b = new double;
	n = new double;
	e = new double;

	*a = 0.0;
	*b = 0.0;
	*n = 0.0;
	*e = 0.0;

	flag_normal_int = false;

	x_Ai = new Matrix(3);
	x_Bi = new Matrix(3);
	Q_Ai = new Matrix(3, 3);
	Q_Bi = new Matrix(3, 3);
	
	d = new Matrix(12);
	dui = new Matrix(12);
	ddui = new Matrix(12);
	
	I3 = new Matrix(3, 3);

	(*I3)(0, 0) = 1.0;
	(*I3)(1, 1) = 1.0;
	(*I3)(2, 2) = 1.0;

	x_Ap = new Matrix(3);
	x_Bp = new Matrix(3);
	Q_Ap = new Matrix(3, 3);
	Q_Bp = new Matrix(3, 3);

	Q0A = new Matrix(3, 3);
	Q0B = new Matrix(3, 3);

	QAi = new double*[3];
	for (int i = 0; i < 3; i++)
		QAi[i] = new double[3];
	QBi = new double*[3];
	for (int i = 0; i < 3; i++)
		QBi[i] = new double[3];
}

void SECylinder::FreeSpecific()
{
	delete bv;

	delete xAf;
	delete xBf;

	delete x_Ai;
	delete x_Bi;
	delete Q_Ai;
	delete Q_Bi;

	delete Q0A;
	delete Q0B;
	delete d;
	delete dui;
	delete ddui;

	delete I3;

	delete x_Ap;
	delete x_Bp;
	delete Q_Ap;
	delete Q_Bp;

	delete a;
	delete b;
	delete n;
	delete e;

	for (int i = 0; i < 3; i++)
		delete[]QAi[i];
	delete[]QAi;
	for (int i = 0; i < 3; i++)
		delete[]QBi[i];
	delete[]QBi;
}

SECylinder::~SECylinder()
{
	Free();
	FreeSpecific();
}

bool SECylinder::Read(FILE* f)
{
	char s[1000];
	fscanf(f, "%s", s);
	number = atoi(s);
	
	fscanf(f, "%s", s);
	if (!strcmp(s, "Material"))
	{
		fscanf(f, "%s", s);
		material = atoi(s);
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "A"))
	{
		fscanf(f, "%s", s);
		*a = atof(s);
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "B"))
	{
		fscanf(f, "%s", s);
		*b = atof(s);
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "N"))
	{
		fscanf(f, "%s", s);
		*n = atof(s);
	}
	else
		return false;
	//Leitura dos sistemas de coordenadas para criação da superfície
	//Salva a posição (stream)
	fpos_t pos;
	fgetpos(f, &pos);
	//Caso 1 - um único sistema
	fscanf(f, "%s", s);
	if (!strcmp(s, "CS"))
	{
		fscanf(f, "%s", s);
		csA = atoi(s);
		csB = atoi(s);
	}
	else
	{
		fsetpos(f, &pos);
		//Caso 2 - dois sistemas
		fscanf(f, "%s", s);
		if (!strcmp(s, "CSA"))
		{
			fscanf(f, "%s", s);
			csA = atoi(s);
		}
		else
			return false;
		fscanf(f, "%s", s);
		if (!strcmp(s, "CSB"))
		{
			fscanf(f, "%s", s);
			csB = atoi(s);
		}
		else
			return false;
	}

	fscanf(f, "%s", s);
	if (!strcmp(s, "NormalInterior"))
		flag_normal_int = true;
	else
	{
		if (!strcmp(s, "NormalExterior"))
			flag_normal_int = false;
		else
			return false;
	}

	fscanf(f, "%s", s);
	if (!strcmp(s, "Nodes"))
	{
		fscanf(f, "%s", s);
		nodes[0] = atoi(s);
		fscanf(f, "%s", s);
		nodes[1] = atoi(s);
	}
	else
		return false;

	//Degeneration
	fgetpos(f, &pos);
	fscanf(f, "%s", s);
	if (!strcmp(s, "Degeneration"))
	{
		fgetpos(f, &pos);
		fscanf(f, "%s", s);
		if (!strcmp(s, "Coord1"))
		{
			fscanf(f, "%s", s);
			deg_u1 = true;
			deg_u1_value = atof(s);
		}
		else
			fsetpos(f, &pos);

		fgetpos(f, &pos);
		fscanf(f, "%s", s);
		if (!strcmp(s, "Coord2"))
		{
			fscanf(f, "%s", s);
			deg_u2 = true;
			deg_u2_value = atof(s);
		}
		else
			fsetpos(f, &pos);
	}
	else
		fsetpos(f, &pos);


	return true;
}

void SECylinder::Write(FILE *f)
{
	char s[20];
	if (flag_normal_int == true)
		sprintf(s, "NormalInterior");
	else
		sprintf(s, "NormalExterior");
	fprintf(f, "SECylinder\t%d\tMaterial\t%d\tA\t%.6e\tB\t%.6e\tN\t%.6e\tCSA\t%d\tCSB\t%d\t%s\tNodes\t%d\t%d\n", number, material, *a, *b, *n, csA, csB, s, nodes[0], nodes[1]);
	//PS: degeneration info is not printed
}

bool SECylinder::Check()
{
	for (int i = 0; i < n_nodes; i++)
	{
		if (nodes[i] > db.number_nodes)
			return false;
	}
	if (csA > db.number_CS)
		return false;
	if (csB > db.number_CS)
		return false;

	if (material > db.number_materials)
		return false;

	return true;
}
void SECylinder::PreCalc()
{
	//Apontando para posição que indica valor dos GLs globais
	for (int i = 0; i < 6; i++)
	{
		GLs[i] = &db.nodes[nodes[0] - 1]->GLs[i];
		GLs[i + 6] = &db.nodes[nodes[1] - 1]->GLs[i];
	}
	*e = 2.0 / *n;	//expoente

	//Transformação de coordenadas
	Matrix e1g(3);
	e1g(0, 0) = 1.0;
	Matrix e2g(3);
	e2g(1, 0) = 1.0;
	Matrix e3g(3);
	e3g(2, 0) = 1.0;
	Matrix e1l = *db.CS[csA - 1]->E1;
	Matrix e2l = *db.CS[csA - 1]->E2;
	Matrix e3l = *db.CS[csA - 1]->E3;

	//Salva a matriz de transformação de coordenadas (para orientar o plano da ST de acordo com a orientação de referência da superfície)
	(*Q0A)(0, 0) = dot(e1g, e1l);
	(*Q0A)(0, 1) = dot(e1g, e2l);
	(*Q0A)(0, 2) = dot(e1g, e3l);

	(*Q0A)(1, 0) = dot(e2g, e1l);
	(*Q0A)(1, 1) = dot(e2g, e2l);
	(*Q0A)(1, 2) = dot(e2g, e3l);

	(*Q0A)(2, 0) = dot(e3g, e1l);
	(*Q0A)(2, 1) = dot(e3g, e2l);
	(*Q0A)(2, 2) = dot(e3g, e3l);

	e1l = *db.CS[csB - 1]->E1;
	e2l = *db.CS[csB - 1]->E2;
	e3l = *db.CS[csB - 1]->E3;

	//Salva a matriz de transformação de coordenadas (para orientar o plano da ST de acordo com a orientação de referência da superfície)
	(*Q0B)(0, 0) = dot(e1g, e1l);
	(*Q0B)(0, 1) = dot(e1g, e2l);
	(*Q0B)(0, 2) = dot(e1g, e3l);

	(*Q0B)(1, 0) = dot(e2g, e1l);
	(*Q0B)(1, 1) = dot(e2g, e2l);
	(*Q0B)(1, 2) = dot(e2g, e3l);

	(*Q0B)(2, 0) = dot(e3g, e1l);
	(*Q0B)(2, 1) = dot(e3g, e2l);
	(*Q0B)(2, 2) = dot(e3g, e3l);

	//////////////////////////////////////Bounding Volumes////////////////////////////////////
	if (db.gcs_exist)
		inc_len_factor = db.gcs->inc_len_factor;

	//Computing the largest gnb available, from interface laws (to avoid smaller BVs of small geometric entities)
	largest_gnb = 0.0;
	for (int i = 0; i < db.number_contactinterfaces; i++)
	{
		Interface_1* ptr = static_cast<Interface_1*>(db.contactinterfaces[i]);
		if (ptr->gnb > largest_gnb)
			largest_gnb = (float)ptr->gnb;
	}
	rad = 0;
	if ((float)*a > rad)
		rad = (float)*a;
	if ((float)*b > rad)
		rad = (float)*b;
	bv_offset = (1.0f + inc_len_factor) * (rad + largest_gnb);	//offset a ser usado no BV do BodyGeometry

	BoundingCylinder* ptr_bv = static_cast<BoundingCylinder*>(bv);
	
	ptr_bv->radius = (1.0f + inc_len_factor) * (rad + largest_gnb);
	ptr_bv->ref_radius = (1.0f + inc_len_factor) * (rad + largest_gnb);
	//ptr_bv->size = ptr_bv->radius;
	ptr_bv->inc_len_factor = inc_len_factor;
	//Associated entity:
	//associated_ID is set on PreCalc of BodyGeometry class
	ptr_bv->associated_type = 'G';
	ptr_bv->associated_sub_ID = number;
	
	UpdateVariables();
	UpdateBoundingVolumes();
	SaveLagrange();

	//AceGen Mirrors
	xAi = x_Ai->getMatrix();
	xBi = x_Bi->getMatrix();
}

void SECylinder::UpdateVariables()
{
	Matrix alpha_A(3);
	Matrix alpha_B(3);
	Matrix u_A(3);
	Matrix u_B(3);

	for (int i = 0; i < 3; i++)
	{
		(*d)(i + 0, 0) = db.nodes[nodes[0] - 1]->displacements[i];
		(*d)(i + 3, 0) = db.nodes[nodes[0] - 1]->displacements[i + 3];
		(*d)(i + 6, 0) = db.nodes[nodes[1] - 1]->displacements[i];
		(*d)(i + 9, 0) = db.nodes[nodes[1] - 1]->displacements[i + 3];

		u_A(i, 0) = db.nodes[nodes[0] - 1]->displacements[i];
		u_B(i, 0) = db.nodes[nodes[1] - 1]->displacements[i];
		alpha_A(i, 0) = db.nodes[nodes[0] - 1]->displacements[i + 3];
		alpha_B(i, 0) = db.nodes[nodes[1] - 1]->displacements[i + 3];
	}
	//Q_AAp
	double alpha = norm(alpha_A);							//Valor escalar do parametro alpha
	Matrix A = skew(alpha_A);								//Matriz A
	double g = 4.0 / (4.0 + alpha * alpha);					//função g(alpha) - em algumas ref. tb. chamado de h(alpha)
	Matrix Qd = *I3 + g * (A + 0.5*(A*A));					//Tensor de rotação
	*Q_Ap = Qd * (*db.nodes[nodes[0] - 1]->Q) * (*Q0A);
	//Q_BAp
	alpha = norm(alpha_B);									//Valor escalar do parametro alpha
	A = skew(alpha_B);										//Matriz A
	g = 4.0 / (4.0 + alpha * alpha);						//função g(alpha) - em algumas ref. tb. chamado de h(alpha)
	Qd = *I3 + g * (A + 0.5*(A*A));							//Tensor de rotação
	*Q_Bp = Qd * (*db.nodes[nodes[1] - 1]->Q) * (*Q0B);
	//x_Ap
	*x_Ap = *x_Ai + u_A;
	//x_Bp
	*x_Bp = *x_Bi + u_B;
}

void SECylinder::UpdateBoundingVolumes()
{
	//Conversão para single precision
	*xAf = *x_Ap;
	*xBf = *x_Bp;
	
	MatrixFloat dir = (1.0f / norm(*xBf - *xAf)) * (*xBf - *xAf);
	//Updating the bounding volume
	BoundingCylinder* ptr_bv = static_cast<BoundingCylinder*>(bv);
	*ptr_bv->xb = *xAf - (rad + largest_gnb) * dir;
	*ptr_bv->xt = *xBf + (rad + largest_gnb) * dir;

	//Setting max and min (used in AABB of BodyGeometry)
	for (int i = 0; i < 3; i++)
	{
		if ((*ptr_bv->xb)(i, 0) > (*ptr_bv->xt)(i, 0))
		{
			max[i] = (*ptr_bv->xb)(i, 0);
			min[i] = (*ptr_bv->xt)(i, 0);
		}
		else
		{
			min[i] = (*ptr_bv->xb)(i, 0);
			max[i] = (*ptr_bv->xt)(i, 0);
		}
	}
	//Updating center and size (used in Verlet/LinkedCells schemes)
	float half_len =  (1.0f + inc_len_factor) * (0.5f * norm(*ptr_bv->xb - *ptr_bv->xt) );
	ptr_bv->size = sqrt(ptr_bv->radius * ptr_bv->radius + half_len * half_len);
	ptr_bv->x_center[0] = 0.5f * ((*ptr_bv->xb)(0, 0) + (*ptr_bv->xt)(0, 0));
	ptr_bv->x_center[1] = 0.5f * ((*ptr_bv->xb)(1, 0) + (*ptr_bv->xt)(1, 0));
	ptr_bv->x_center[2] = 0.5f * ((*ptr_bv->xb)(2, 0) + (*ptr_bv->xt)(2, 0));
	
}

void SECylinder::SaveLagrange()
{
	//Do not call UpdateVaribles because node data are updated first!
	*x_Ai = *x_Ap;
	*x_Bi = *x_Bp;
	*Q_Ai = *Q_Ap;
	*Q_Bi = *Q_Bp;

	//Computes dui, ddui - SaveLagrange from "Node" is called first!
	for (int i = 0; i < 3; i++)
	{
		(*dui)(i + 0, 0) = db.nodes[nodes[0] - 1]->copy_vel[i];
		(*dui)(i + 3, 0) = db.nodes[nodes[0] - 1]->copy_vel[i + 3];
		(*dui)(i + 6, 0) = db.nodes[nodes[1] - 1]->copy_vel[i];
		(*dui)(i + 9, 0) = db.nodes[nodes[1] - 1]->copy_vel[i + 3];

		(*ddui)(i + 0, 0) = db.nodes[nodes[0] - 1]->copy_accel[i];
		(*ddui)(i + 3, 0) = db.nodes[nodes[0] - 1]->copy_accel[i + 3];
		(*ddui)(i + 6, 0) = db.nodes[nodes[1] - 1]->copy_accel[i];
		(*ddui)(i + 9, 0) = db.nodes[nodes[1] - 1]->copy_accel[i + 3];

		(*x_Ai)(i, 0) = db.nodes[nodes[0] - 1]->copy_coordinates[i];
		(*x_Bi)(i, 0) = db.nodes[nodes[1] - 1]->copy_coordinates[i];

		(*d)(i + 0, 0) = 0.0;
		(*d)(i + 3, 0) = 0.0;
		(*d)(i + 6, 0) = 0.0;
		(*d)(i + 9, 0) = 0.0;
	}

	//Saving bounding volume
	bv->SaveConfiguration();

	//AceGen Mirrors
	Q_Ai->MatrixToPtr(QAi, 3);
	Q_Bi->MatrixToPtr(QBi, 3);
}
void SECylinder::WriteVTK_XMLRender(FILE *f)
{
	if (db.post_files->WriteFlexibleContactSurfaces_flag == true)
	{
		int n_circ = 24;
		double theta = 0;
		double phi;
		double** points;
		points = new double*[n_circ];
		for (int i = 0; i < n_circ; i++)
			points[i] = new double[2];
		for (int index = 0; index < n_circ; index++)
		{
			theta = (index * 2 * PI) / n_circ;
			phi = 1.0 / (pow((pow(abs(sin(theta)), *n) + pow(abs(cos(theta)), *n)), 1.0 / *n));
			points[index][0] = phi * *a*cos(theta);
			points[index][1] = phi * *b*sin(theta);
		}
		//vetores para escrita no formato binário - usando a função 'enconde'
		std::vector<float> float_vector;
		std::vector<int> int_vector;

		Matrix vec_P(3);
		//Número de pontos a serem gerados
		int n_points = n_circ * 2;
		//Número de células a serem geradas
		int n_cells = n_circ;
		//Opens Piece
		fprintf(f, "\t\t<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", n_points, n_cells);
		//Opens Points
		fprintf(f, "\t\t\t<Points>\n");
		//Opens DataArray
		fprintf(f, "\t\t\t\t<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"binary\">\n");
		float_vector.clear();
		//Preenchendo as coordenadas dos pontos
		for (int i = 0; i < 2; i++)//percorrendo os 2 nós das extremidades
		{
			for (int point = 0; point < n_circ; point++)//Percorre os nós que descrevem o perímetro da ST
			{
				//Posição de cada ponto P no plano xy (referência)
				vec_P(0, 0) = points[point][0];
				vec_P(1, 0) = points[point][1];
				vec_P(2, 0) = 0.0;
				if (i == 0)
					vec_P = (*Q_Ai)*vec_P;//Operando rotacionando para o sistema da barra
				else
					vec_P = (*Q_Bi)*vec_P;//Operando rotacionando para o sistema da barra
				for (int c = 0; c < 3; c++)
					vec_P(c, 0) += db.nodes[nodes[i] - 1]->ref_coordinates[c] + db.post_files->mag_factor*(db.nodes[nodes[i] - 1]->copy_coordinates[c] - db.nodes[nodes[i] - 1]->ref_coordinates[c]);//Translação - soma a posição do centro da barra (do ponto em questão)

				float_vector.push_back((float)(vec_P(0, 0)));
				float_vector.push_back((float)(vec_P(1, 0)));
				float_vector.push_back((float)(vec_P(2, 0)));
			}
		}
		fprintf(f, encodeData<float>(float_vector).c_str());
		fprintf(f, "\n");
		//Closes DataArray
		fprintf(f, "\t\t\t\t</DataArray>\n");
		//Closes Points
		fprintf(f, "\t\t\t</Points>\n");

		int nodes[4];
		//Opens Cells
		fprintf(f, "\t\t\t<Cells>\n");
		//Opens DataArray
		fprintf(f, "\t\t\t\t<DataArray type=\"Int32\" Name=\"connectivity\" format=\"binary\">\n");
		int_vector.clear();
		for (int cell = 0; cell < n_cells; cell++)
		{
			if (cell != n_cells - 1)
			{
				nodes[0] = cell;
				nodes[1] = cell + n_cells;
				nodes[2] = cell + n_cells + 1;
				nodes[3] = cell + 1;
			}
			else
			{
				nodes[0] = cell;
				nodes[1] = 2 * n_cells - 1;
				nodes[2] = n_cells;
				nodes[3] = 0;
			}
			int_vector.push_back(nodes[0]);
			int_vector.push_back(nodes[1]);
			int_vector.push_back(nodes[2]);
			int_vector.push_back(nodes[3]);
		}

		fprintf(f, encodeData(int_vector).c_str());
		fprintf(f, "\n");
		//Closes DataArray
		fprintf(f, "\t\t\t\t</DataArray>\n");
		//Opens DataArray
		int_vector.clear();
		fprintf(f, "\t\t\t\t<DataArray type=\"Int32\" Name=\"types\" format=\"binary\">\n");
		for (int cell = 0; cell < n_cells; cell++)
			int_vector.push_back(9);
		fprintf(f, encodeData(int_vector).c_str());
		fprintf(f, "\n");
		//Closes DataArray
		fprintf(f, "\t\t\t\t</DataArray>\n");
		//Opens DataArray
		fprintf(f, "\t\t\t\t<DataArray type=\"Int32\" Name=\"offsets\" format=\"binary\">\n");
		int cur_off = 0;
		int_vector.clear();
		for (int cell = 0; cell < n_cells; cell++)
		{
			cur_off += 4;
			int_vector.push_back(cur_off);
		}
		fprintf(f, encodeData(int_vector).c_str());
		fprintf(f, "\n");
		//Closes DataArray
		fprintf(f, "\t\t\t\t</DataArray>\n");
		//Closes Cells
		fprintf(f, "\t\t\t</Cells>\n");

		//Opens CellData
		fprintf(f, "\t\t\t<CellData FieldData=\"SurfaceData\">\n");
		int_vector.clear();
		//Opens DataArray
		fprintf(f, "\t\t\t\t<DataArray Name=\"SurfaceProperties\" type=\"Int32\" NumberOfComponents=\"%d\" format=\"binary\">\n", 1);
		for (int cell = 0; cell < n_cells; cell++)
		{
			int_vector.push_back(material);		//Material ID
		}
		fprintf(f, encodeData(int_vector).c_str());
		fprintf(f, "\n");
		//Closes DataArray
		fprintf(f, "\t\t\t\t</DataArray>\n");
		//Closes CellData
		fprintf(f, "\t\t\t</CellData>\n");


		//Closes Piece
		fprintf(f, "\t\t</Piece>\n");
		for (int i = 0; i < n_circ; i++)
			delete[] points[i];
		delete[]points;
	}
}

void SECylinder::SurfacePosition(Matrix* Gamma, double z, double th, bool next)
{
	double v[234];
	double c[2];
	c[0] = z;
	c[1] = th;
	double dof[12];
	if (next)
	{
		for (int i = 0; i < 12; i++)
			dof[i] = (*d)(i, 0);
	}
	else
	{
		for (int i = 0; i < 12; i++)
			dof[i] = 0.0;
	}
	double Gam[3];

	//AceGen code
	v[1] = c[0];
	v[2] = c[1];
	v[113] = sin(v[2]);
	v[111] = cos(v[2]);
	v[6] = dof[3];
	v[59] = (v[6] * v[6]);
	v[7] = dof[4];
	v[119] = v[7] / 2e0;
	v[57] = v[119] * v[6];
	v[52] = (v[7] * v[7]);
	v[8] = dof[5];
	v[64] = v[119] * v[8];
	v[62] = (v[6] * v[8]) / 2e0;
	v[53] = (v[8] * v[8]);
	v[121] = v[52] + v[53];
	v[12] = dof[9];
	v[75] = (v[12] * v[12]);
	v[13] = dof[10];
	v[120] = v[13] / 2e0;
	v[73] = v[12] * v[120];
	v[68] = (v[13] * v[13]);
	v[14] = dof[11];
	v[80] = v[120] * v[14];
	v[78] = (v[12] * v[14]) / 2e0;
	v[69] = (v[14] * v[14]);
	v[123] = v[68] + v[69];
	v[24] = QAi[0][0];
	v[25] = QAi[0][1];
	v[27] = QAi[1][0];
	v[28] = QAi[1][1];
	v[30] = QAi[2][0];
	v[31] = QAi[2][1];
	v[33] = QBi[0][0];
	v[34] = QBi[0][1];
	v[36] = QBi[1][0];
	v[37] = QBi[1][1];
	v[39] = QBi[2][0];
	v[40] = QBi[2][1];
	v[51] = 4e0 / (4e0 + v[121] + v[59]);
	v[122] = -0.5e0*v[51];
	v[54] = 1e0 + v[121] * v[122];
	v[55] = v[51] * (v[57] - v[8]);
	v[56] = v[51] * (v[62] + v[7]);
	v[58] = v[51] * (v[57] + v[8]);
	v[60] = 1e0 + v[122] * (v[53] + v[59]);
	v[61] = v[51] * (-v[6] + v[64]);
	v[63] = v[51] * (v[62] - v[7]);
	v[65] = v[51] * (v[6] + v[64]);
	v[66] = 1e0 - v[122] * (-v[52] - v[59]);
	v[67] = 4e0 / (4e0 + v[123] + v[75]);
	v[124] = -0.5e0*v[67];
	v[70] = 1e0 + v[123] * v[124];
	v[71] = v[67] * (-v[14] + v[73]);
	v[72] = v[67] * (v[13] + v[78]);
	v[74] = v[67] * (v[14] + v[73]);
	v[76] = 1e0 + v[124] * (v[69] + v[75]);
	v[77] = v[67] * (-v[12] + v[80]);
	v[79] = v[67] * (-v[13] + v[78]);
	v[81] = v[67] * (v[12] + v[80]);
	v[82] = 1e0 - v[124] * (-v[68] - v[75]);
	v[107] = (1e0 - v[1]) / 2e0;
	v[108] = (1e0 + v[1]) / 2e0;
	v[109] = 2e0 / (*e);
	v[110] = 1e0 / Power(Power(fabs(v[111]), v[109]) + Power(fabs(v[113]), v[109]), 1 / v[109]);
	v[112] = (*a)*v[110] * v[111];
	v[114] = (*b)*v[110] * v[113];
	Gam[0] = v[107] * (dof[0] + v[112] * (v[24] * v[54] + v[27] * v[55] + v[30] * v[56]) + v[114] * (v[25] * v[54] + v[28] * v[55]
		+ v[31] * v[56]) + xAi[0]) + v[108] * (dof[6] + v[112] * (v[33] * v[70] + v[36] * v[71] + v[39] * v[72]) + v[114] *
		(v[34] * v[70] + v[37] * v[71] + v[40] * v[72]) + xBi[0]);
	Gam[1] = v[107] * (dof[1] + v[112] * (v[24] * v[58] + v[27] * v[60] + v[30] * v[61]) + v[114] * (v[25] * v[58] + v[28] * v[60]
		+ v[31] * v[61]) + xAi[1]) + v[108] * (dof[7] + v[112] * (v[33] * v[74] + v[36] * v[76] + v[39] * v[77]) + v[114] *
		(v[34] * v[74] + v[37] * v[76] + v[40] * v[77]) + xBi[1]);
	Gam[2] = v[107] * (dof[2] + v[112] * (v[24] * v[63] + v[27] * v[65] + v[30] * v[66]) + v[114] * (v[25] * v[63] + v[28] * v[65]
		+ v[31] * v[66]) + xAi[2]) + v[108] * (dof[8] + v[112] * (v[33] * v[79] + v[36] * v[81] + v[39] * v[82]) + v[114] *
		(v[34] * v[79] + v[37] * v[81] + v[40] * v[82]) + xBi[2]);

	//Copying final result
	(*Gamma)(0, 0) = Gam[0];
	(*Gamma)(1, 0) = Gam[1];
	(*Gamma)(2, 0) = Gam[2];
}

void SECylinder::NormalExt(double z, double th, Matrix* n, bool next)
{
	double v[262];
	double c[2];
	c[0] = z;
	c[1] = th;
	double dof[12];
	if (next)
	{
		for (int i = 0; i < 12; i++)
			dof[i] = (*d)(i, 0);
	}
	else
	{
		for (int i = 0; i < 12; i++)
			dof[i] = 0.0;
	}
	double normal[3];
	bool* normalint = &flag_normal_int;
	
	//AceGen code
	v[1] = c[0];
	v[2] = c[1];
	v[119] = cos(v[2]);
	v[114] = sin(v[2]);
	v[123] = fabs(v[114]);
	v[121] = fabs(v[119]);
	v[6] = dof[3];
	v[60] = (v[6] * v[6]);
	v[7] = dof[4];
	v[148] = v[7] / 2e0;
	v[58] = v[148] * v[6];
	v[53] = (v[7] * v[7]);
	v[8] = dof[5];
	v[65] = v[148] * v[8];
	v[63] = (v[6] * v[8]) / 2e0;
	v[54] = (v[8] * v[8]);
	v[150] = v[53] + v[54];
	v[12] = dof[9];
	v[76] = (v[12] * v[12]);
	v[13] = dof[10];
	v[149] = v[13] / 2e0;
	v[74] = v[12] * v[149];
	v[69] = (v[13] * v[13]);
	v[14] = dof[11];
	v[81] = v[14] * v[149];
	v[79] = (v[12] * v[14]) / 2e0;
	v[70] = (v[14] * v[14]);
	v[152] = v[69] + v[70];
	v[15] = (*a);
	v[16] = (*b);
	v[24] = QAi[0][0];
	v[25] = QAi[0][1];
	v[27] = QAi[1][0];
	v[28] = QAi[1][1];
	v[30] = QAi[2][0];
	v[31] = QAi[2][1];
	v[33] = QBi[0][0];
	v[34] = QBi[0][1];
	v[36] = QBi[1][0];
	v[37] = QBi[1][1];
	v[39] = QBi[2][0];
	v[40] = QBi[2][1];
	v[52] = 4e0 / (4e0 + v[150] + v[60]);
	v[151] = -0.5e0*v[52];
	v[55] = 1e0 + v[150] * v[151];
	v[56] = v[52] * (v[58] - v[8]);
	v[57] = v[52] * (v[63] + v[7]);
	v[59] = v[52] * (v[58] + v[8]);
	v[61] = 1e0 + v[151] * (v[54] + v[60]);
	v[62] = v[52] * (-v[6] + v[65]);
	v[64] = v[52] * (v[63] - v[7]);
	v[66] = v[52] * (v[6] + v[65]);
	v[67] = 1e0 - v[151] * (-v[53] - v[60]);
	v[68] = 4e0 / (4e0 + v[152] + v[76]);
	v[153] = -0.5e0*v[68];
	v[71] = 1e0 + v[152] * v[153];
	v[72] = v[68] * (-v[14] + v[74]);
	v[73] = v[68] * (v[13] + v[79]);
	v[75] = v[68] * (v[14] + v[74]);
	v[77] = 1e0 + v[153] * (v[70] + v[76]);
	v[78] = v[68] * (-v[12] + v[81]);
	v[80] = v[68] * (-v[13] + v[79]);
	v[82] = v[68] * (v[12] + v[81]);
	v[83] = 1e0 - v[153] * (-v[69] - v[76]);
	v[84] = v[24] * v[55] + v[27] * v[56] + v[30] * v[57];
	v[85] = v[25] * v[55] + v[28] * v[56] + v[31] * v[57];
	v[87] = v[24] * v[59] + v[27] * v[61] + v[30] * v[62];
	v[88] = v[25] * v[59] + v[28] * v[61] + v[31] * v[62];
	v[90] = v[24] * v[64] + v[27] * v[66] + v[30] * v[67];
	v[91] = v[25] * v[64] + v[28] * v[66] + v[31] * v[67];
	v[93] = v[33] * v[71] + v[36] * v[72] + v[39] * v[73];
	v[94] = v[34] * v[71] + v[37] * v[72] + v[40] * v[73];
	v[96] = v[33] * v[75] + v[36] * v[77] + v[39] * v[78];
	v[97] = v[34] * v[75] + v[37] * v[77] + v[40] * v[78];
	v[99] = v[33] * v[80] + v[36] * v[82] + v[39] * v[83];
	v[100] = v[34] * v[80] + v[37] * v[82] + v[40] * v[83];
	v[108] = (1e0 - v[1]) / 2e0;
	v[109] = (1e0 + v[1]) / 2e0;
	v[110] = 2e0 / (*e);
	v[122] = -1e0 + v[110];
	v[120] = Power(v[121], v[110]) + Power(v[123], v[110]);
	v[124] = Power(v[120], -1e0 - 1e0 / v[110])*(-(v[119] * Power(v[123], v[122])*_copysign(1.e0, v[114]))
		+ v[114] * Power(v[121], v[122])*_copysign(1.e0, v[119]));
	v[111] = 1e0 / Power(v[120], 1 / v[110]);
	v[156] = v[111] * v[114];
	v[155] = v[111] * v[119];
	v[126] = (v[114] * v[124] + v[155])*v[16];
	v[125] = v[15] * (v[119] * v[124] - v[156]);
	v[129] = v[108] * (v[125] * v[90] + v[126] * v[91]) + v[109] * (v[100] * v[126] + v[125] * v[99]);
	v[128] = v[108] * (v[125] * v[87] + v[126] * v[88]) + v[109] * (v[125] * v[96] + v[126] * v[97]);
	v[127] = v[108] * (v[125] * v[84] + v[126] * v[85]) + v[109] * (v[125] * v[93] + v[126] * v[94]);
	v[113] = v[15] * v[155];
	v[115] = v[156] * v[16];
	v[138] = (-dof[2] + dof[8] + v[100] * v[115] - v[113] * v[90] - v[115] * v[91] + v[113] * v[99] - xAi[2] + xBi[2]) / 2e0;
	v[135] = (-dof[1] + dof[7] - v[113] * v[87] - v[115] * v[88] + v[113] * v[96] + v[115] * v[97] - xAi[1] + xBi[1]) / 2e0;
	v[132] = (-dof[0] + dof[6] - v[113] * v[84] - v[115] * v[85] + v[113] * v[93] + v[115] * v[94] - xAi[0] + xBi[0]) / 2e0;
	v[145] = -(v[128] * v[132]) + v[127] * v[135];
	v[143] = v[129] * v[132] - v[127] * v[138];
	v[140] = -(v[129] * v[135]) + v[128] * v[138];
	v[157] = ((*normalint) ? -1 : 1) / sqrt((v[140] * v[140]) + (v[143] * v[143]) + (v[145] * v[145]));
	normal[0] = v[140] * v[157];
	normal[1] = v[143] * v[157];
	normal[2] = v[145] * v[157];

	//Copying final result
	(*n)(0, 0) = normal[0];
	(*n)(1, 0) = normal[1];
	(*n)(2, 0) = normal[2];
}