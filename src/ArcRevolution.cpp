#include "ArcRevolution.h"

#include "BoundingBoxAxesOriented.h"
#include "MatrixFloat.h"
#include "ArcCirc.h"
#include "Interface_1.h"
#include "Node.h"
#include "CoordinateSystem.h"
#include "GeneralContactSearch.h"
#include "PostFiles.h"
#include "Encoding.h"

#include"Database.h"
//Variaveis globais
extern
Database db;


ArcRevolution::ArcRevolution()
{
	super_node = 0;

	nDOFs = 6;
	n_nodes = 1;

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

	x_min = 0.0f;
	x_max = 0.0f;
	y_min = 0.0f;
	y_max = 0.0f;
}

ArcRevolution::~ArcRevolution()
{
	Free();
	FreeSpecific();
}

void ArcRevolution::AllocSpecific()
{
	bv = new BoundingBoxAxesOriented();
	sprintf(type_name, "ArcRevolution");
	bv_offset = 0.0f;

	//Variaveis para calcular estado atual (nas funcoes de bounding volume)
	xAf = new MatrixFloat(3);

	flag_normal_int = false;

	x_Ai = new Matrix(3);
	Q_Ai = new Matrix(3, 3);

	Q0 = new Matrix(3, 3);
	Qf = new MatrixFloat(3, 3);
	local_half_widths = new MatrixFloat(3);

	d = new Matrix(6);
	dui = new Matrix(6);
	ddui = new Matrix(6);

	I3 = new Matrix(3, 3);

	(*I3)(0, 0) = 1.0;
	(*I3)(1, 1) = 1.0;
	(*I3)(2, 2) = 1.0;

	x_Ap = new Matrix(3);
	Q_Ap = new Matrix(3, 3);

	QAi = new double*[3];
	for (int i = 0; i < 3; i++)
		QAi[i] = new double[3];

	center_local = new Matrix(3);
	
}

void ArcRevolution::FreeSpecific()
{
	delete bv;

	delete xAf;

	delete x_Ai;
	delete Q_Ai;

	delete Q0;
	delete Qf;
	delete local_half_widths;

	delete d;
	delete dui;
	delete ddui;

	delete I3;

	delete x_Ap;
	delete Q_Ap;

	for (int i = 0; i < 3; i++)
		delete[]QAi[i];
	delete[]QAi;

	delete center_local;
}
bool ArcRevolution::Read(FILE *f)
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
	if (!strcmp(s, "Arc"))
	{
		fscanf(f, "%s", s);
		arc_ID = atoi(s);
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "CS"))
	{
		fscanf(f, "%s", s);
		cs = atoi(s);
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "Node"))
	{
		fscanf(f, "%s", s);
		nodes[0] = atoi(s);
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "Concave"))
		flag_normal_int = true;
	else
	{
		if (!strcmp(s, "Convex"))
			flag_normal_int = false;
		else
			return false;
	}

	fpos_t pos;
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
void ArcRevolution::Write(FILE *f)
{
	char s[20];
	if (flag_normal_int == true)
		sprintf(s, "Concave");
	else
		sprintf(s, "Convex");
	fprintf(f, "ArcRevolution\t%d\tArc\t%d\tCS\t%d\tNode\t%d\t%s\n",
		number,
		arc_ID,
		cs,
		nodes[0],
		s);
	//PS: degeneration info is not printed
}
bool ArcRevolution::Check()
{
	for (int i = 0; i < n_nodes; i++)
	{
		if (nodes[i] > db.number_nodes)
			return false;
	}
	if (cs > db.number_CS)
		return false;
	if (arc_ID > db.number_arcs)
		return false;

	if (material > db.number_materials)
		return false;

	return true;
}
void ArcRevolution::PreCalc()
{
	//Atribuindo valores do arco - ponteiros
	radius = &db.arcs[arc_ID - 1]->radius;
	c_point = &db.arcs[arc_ID - 1]->c_point;
	i_point = &db.arcs[arc_ID - 1]->i_point;
	f_point = &db.arcs[arc_ID - 1]->f_point;
	theta_i = &db.arcs[arc_ID - 1]->theta_i;
	theta_f = &db.arcs[arc_ID - 1]->theta_f;
	if ((*c_point)(0, 0) < (*radius))
	{
		theta_valid_min = -acos(-(*c_point)(0, 0) / (*radius));
		theta_valid_max = +acos(-(*c_point)(0, 0) / (*radius));
	}
	else
	{
		theta_valid_min = -PI;
		theta_valid_max = +PI;
	}

	//if (number == 14)
	//	int test = 1;

	//Apontando para posição que indica valor dos GLs globais
	for (int i = 0; i < 6; i++)
	{
		GLs[i] = &db.nodes[nodes[0] - 1]->GLs[i];
	}
	 
	//Transformação de coordenadas
	Matrix e1g(3);
	e1g(0, 0) = 1.0;
	Matrix e2g(3);
	e2g(1, 0) = 1.0;
	Matrix e3g(3);
	e3g(2, 0) = 1.0;
	Matrix e1l = *db.CS[cs - 1]->E1;
	Matrix e2l = *db.CS[cs - 1]->E2;
	Matrix e3l = *db.CS[cs - 1]->E3;

	//Salva a matriz de transformação de coordenadas (para orientar o plano da ST de acordo com a orientação de referência da superficie)
	(*Q0)(0, 0) = dot(e1g, e1l);
	(*Q0)(0, 1) = dot(e1g, e2l);
	(*Q0)(0, 2) = dot(e1g, e3l);

	(*Q0)(1, 0) = dot(e2g, e1l);
	(*Q0)(1, 1) = dot(e2g, e2l);
	(*Q0)(1, 2) = dot(e2g, e3l);

	(*Q0)(2, 0) = dot(e3g, e1l);
	(*Q0)(2, 1) = dot(e3g, e2l);
	(*Q0)(2, 2) = dot(e3g, e3l);

	//////////////////////////////////////Bounding Volumes////////////////////////////////////
	if (db.gcs_exist)
		inc_len_factor = db.gcs->inc_len_factor;

	//Computing the largest gnb available, from interface laws (to avoid smaller BVs of small geometric entities)
	largest_gnb = 0.0f;
	for (int i = 0; i < db.number_contactinterfaces; i++)
	{
		Interface_1* ptr = static_cast<Interface_1*>(db.contactinterfaces[i]);
		if (ptr->gnb > largest_gnb)
			largest_gnb = (float)ptr->gnb;
	}

	db.arcs[arc_ID - 1]->Center(*center_local);

	//Evaluation of local position for the cylinder extremes and its BVradius
	db.arcs[arc_ID - 1]->BoundingRectangle(&x_min, &x_max, &y_min, &y_max);
	BoundingBoxAxesOriented* ptr_bv = static_cast<BoundingBoxAxesOriented*>(bv);
	//Associated entity:
	//associated_ID is set on PreCalc of BodyGeometry class
	ptr_bv->associated_type = 'G';
	ptr_bv->associated_sub_ID = number;

	(*local_half_widths)(0, 0) = (1.0f + inc_len_factor) * x_max;
	(*local_half_widths)(1, 0) = (1.0f + inc_len_factor) * 0.5f * (y_max - y_min);
	(*local_half_widths)(2, 0) = (1.0f + inc_len_factor) * x_max;

	float max_dimension = (*local_half_widths)(0, 0);
	if ((*local_half_widths)(1, 0) > (*local_half_widths)(0, 0))
		max_dimension = (*local_half_widths)(1, 0);

	bv_offset = (1.0f + inc_len_factor) * (max_dimension + largest_gnb);	//offset a ser usado no AABB do BodyGeometry
	ptr_bv->inc_len_factor = inc_len_factor;

	

	UpdateVariables();
	UpdateBoundingVolumes();
	SaveLagrange();

	//AceGen Mirrors
	xAi = x_Ai->getMatrix();
}
void ArcRevolution::UpdateVariables()
{
	Matrix alpha_A(3);
	Matrix u_A(3);
	
	for (int i = 0; i < 3; i++)
	{
		(*d)(i + 0, 0) = db.nodes[nodes[0] - 1]->displacements[i];
		(*d)(i + 3, 0) = db.nodes[nodes[0] - 1]->displacements[i + 3];

		u_A(i, 0) = db.nodes[nodes[0] - 1]->displacements[i];
		alpha_A(i, 0) = db.nodes[nodes[0] - 1]->displacements[i + 3];
	}
	//Q_AAp
	double alpha = norm(alpha_A);							//Valor escalar do parametro alpha
	Matrix A = skew(alpha_A);								//Matriz A
	double g = 4.0 / (4.0 + alpha * alpha);					//função g(alpha) - em algumas ref. tb. chamado de h(alpha)
	Matrix Qd = *I3 + g * (A + 0.5*(A*A));					//Tensor de rotação
	*Q_Ap = Qd * (*db.nodes[nodes[0] - 1]->Q) * (*Q0);
	//x_Ap
	*x_Ap = *x_Ai + u_A;
}
void ArcRevolution::UpdateBoundingVolumes()
{
	//Conversão para single precision
	*xAf = *x_Ap;
	*Qf = *Q_Ap;
	
	//Updating the bounding volume
	BoundingBoxAxesOriented* ptr_bv = static_cast<BoundingBoxAxesOriented*>(bv);
	MatrixFloat local(3);	//centro da caixa no sistema local
	local(1, 0) = 0.5f * (y_min + y_max);
	*ptr_bv->center = *xAf + (*Qf) * local;	//centro da caixa no sistema global
	*ptr_bv->halfwidths = (*local_half_widths);	//dimensões da caixa no sistema local
	
	MatrixFloat e1(3);
	e1(0, 0) = 1.0f;
	MatrixFloat e2(3);
	e2(1, 0) = 1.0f;
	MatrixFloat e3(3);
	e3(2, 0) = 1.0f;
	//Orientação da caixa atual
	*ptr_bv->x_local = (*Qf) * e1;
	*ptr_bv->y_local = (*Qf) * e2;
	*ptr_bv->z_local = (*Qf) * e3;

	//ptr_bv->Report();
	//Setting max and min (used in AABB of BodyGeometry)
	for (int i = 0; i < 3; i++)
	{
		max[i] = (*ptr_bv->center)(i, 0);
		min[i] = (*ptr_bv->center)(i, 0);
	}

	//Updating center and size (used in Verlet/LinkedCells schemes)
	ptr_bv->size = (1.0f + inc_len_factor) * (norm(*ptr_bv->halfwidths) + largest_gnb);
	ptr_bv->x_center[0] = (*ptr_bv->center)(0, 0);
	ptr_bv->x_center[1] = (*ptr_bv->center)(1, 0);
	ptr_bv->x_center[2] = (*ptr_bv->center)(2, 0);

}
void ArcRevolution::SaveLagrange()
{
	//Do not call UpdateVaribles because node data are updated first!
	*x_Ai = *x_Ap;
	*Q_Ai = *Q_Ap;

	//Computes dui, ddui - SaveLagrange from "Node" is called first!
	for (int i = 0; i < 3; i++)
	{
		(*dui)(i + 0, 0) = db.nodes[nodes[0] - 1]->copy_vel[i];
		(*dui)(i + 3, 0) = db.nodes[nodes[0] - 1]->copy_vel[i + 3];
		
		(*ddui)(i + 0, 0) = db.nodes[nodes[0] - 1]->copy_accel[i];
		(*ddui)(i + 3, 0) = db.nodes[nodes[0] - 1]->copy_accel[i + 3];
		
		(*x_Ai)(i, 0) = db.nodes[nodes[0] - 1]->copy_coordinates[i];
		
		(*d)(i + 0, 0) = 0.0;
		(*d)(i + 3, 0) = 0.0;
	}

	//Saving bounding volume
	bv->SaveConfiguration();

	//AceGen Mirrors
	Q_Ai->MatrixToPtr(QAi, 3);
}

void ArcRevolution::WriteVTK_XMLRender(FILE *f)
{
	if (db.post_files->WriteRigidContactSurfaces_flag == true)
	{
		int n_circ1 = 24; //numero de divisões ao longo da revolução
		int n_circ2 = 24; //numero de divisões ao longo do arco
		//Numero de pontos a serem gerados
		int n_points = (2 * n_circ1 + 1) * (2 * n_circ2 + 1);
		//Numero de celulas a serem geradas
		int n_cells = n_circ1 * n_circ2;
		double** points;
		points = new double*[n_points];
		for (int i = 0; i < n_points; i++)
			points[i] = new double[3];
		int index = 0;
		
		double phi_f = 2 * PI;	//angulo de revolução
	
		//
		double theta_i = db.arcs[arc_ID - 1]->theta_i;
		double dtheta = db.arcs[arc_ID - 1]->AngularRange() / (2 * n_circ2);
		double dphi = (phi_f) / (2 * n_circ1);

		double phi = 0;
		double theta;

		//
		for (int m = 0; m < 2 * n_circ1 + 1; m++)
		{
			theta = theta_i;
			for (int n = 0; n < 2 * n_circ2 + 1; n++)
			{
				points[index][0] = ((*radius)*cos(theta) + (*c_point)(0, 0))*cos(phi);
				points[index][1] = (*radius)*sin(theta) + (*c_point)(1, 0);
				points[index][2] = -1.0*((*radius)*cos(theta) + (*c_point)(0, 0))*sin(phi);
				index += 1;
				theta += dtheta;
			}

			phi += dphi;
		}
		//vetores para escrita no formato binario - usando a função 'enconde'
		std::vector<float> float_vector;
		std::vector<int> int_vector;

		Matrix vec_P(3);

		//Opens Piece
		fprintf(f, "\t\t<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", n_points, n_cells);
		//Opens Points
		fprintf(f, "\t\t\t<Points>\n");
		//Opens DataArray
		fprintf(f, "\t\t\t\t<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"binary\">\n");
		float_vector.clear();
		//Preenchendo as coordenadas dos pontos
		Matrix xO(3);
		xO(0, 0) = db.nodes[nodes[0] - 1]->copy_coordinates[0];
		xO(1, 0) = db.nodes[nodes[0] - 1]->copy_coordinates[1];
		xO(2, 0) = db.nodes[nodes[0] - 1]->copy_coordinates[2];
		for (int point = 0; point < n_points; point++)//Percorre os nós que descrevem o perimetro da ST
		{
			//Posição de cada ponto P no plano xy (referência)
			vec_P(0, 0) = points[point][0];
			vec_P(1, 0) = points[point][1];// -2 * (*c_point)(1, 0);
			vec_P(2, 0) = points[point][2];
			vec_P = xO + (*Q_Ai)*vec_P;//Operando rotacionando para o sistema da barra e translação

			float_vector.push_back((float)(vec_P(0, 0)));
			float_vector.push_back((float)(vec_P(1, 0)));
			float_vector.push_back((float)(vec_P(2, 0)));

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
		//Linhas
		int nj = (2 * n_circ2 + 1);
		int_vector.clear();
		for (int i = 0; i < (2 * n_circ1 - 1); i = i + 2)
			for (int j = 0; j < (2 * n_circ2 - 1); j = j + 2)
			{
				int_vector.push_back(j + 0 + (i + 0)*nj);
				int_vector.push_back(j + 0 + (i + 2)*nj);
				int_vector.push_back(j + 2 + (i + 2)*nj);
				int_vector.push_back(j + 2 + (i + 0)*nj);
				int_vector.push_back(j + 0 + (i + 1)*nj);
				int_vector.push_back(j + 1 + (i + 2)*nj);
				int_vector.push_back(j + 2 + (i + 1)*nj);
				int_vector.push_back(j + 1 + (i + 0)*nj);
			}
		fprintf(f, encodeData(int_vector).c_str());
		fprintf(f, "\n");
		//Closes DataArray
		fprintf(f, "\t\t\t\t</DataArray>\n");
		//Opens DataArray
		fprintf(f, "\t\t\t\t<DataArray type=\"Int32\" Name=\"types\" format=\"binary\">\n");
		int_vector.clear();
		for (int i = 0; i < n_cells; i++)
		{
			int_vector.push_back(23);
			//fprintf(f, "\t\t\t\t\t%d\n", 23);
		}
		fprintf(f, encodeData(int_vector).c_str());
		fprintf(f, "\n");
		//Closes DataArray
		fprintf(f, "\t\t\t\t</DataArray>\n");
		//Opens DataArray
		fprintf(f, "\t\t\t\t<DataArray type=\"Int32\" Name=\"offsets\" format=\"binary\">\n");
		int_vector.clear();
		for (int i = 0; i < n_cells; i++)
		{
			int_vector.push_back((i + 1) * 8);
			//fprintf(f, "\t\t\t\t\t%d\n", (i + 1) * 8);
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
		for (int i = 0; i < n_points; i++)
			delete[] points[i];
		delete[]points;
	}
}

void ArcRevolution::SurfacePosition(Matrix* Gamma, double z, double th, bool next)
{
	double* rad = radius;
	double* cpoint = c_point->getMatrix();
	double v[200];
	double c[2];
	c[0] = z;
	c[1] = th;
	double dof[6];
	if (next)
	{
		for (int i = 0; i < 6; i++)
			dof[i] = (*d)(i, 0);
	}
	else
	{
		for (int i = 0; i < 6; i++)
			dof[i] = 0.0;
	}
	double Gam[3];

	//AceGen code
	v[1] = c[0];
	v[2] = c[1];
	v[6] = dof[3];
	v[41] = (v[6] * v[6]);
	v[7] = dof[4];
	v[69] = v[7] / 2e0;
	v[39] = v[6] * v[69];
	v[34] = (v[7] * v[7]);
	v[8] = dof[5];
	v[46] = v[69] * v[8];
	v[44] = (v[6] * v[8]) / 2e0;
	v[35] = (v[8] * v[8]);
	v[70] = v[34] + v[35];
	v[9] = (*rad);
	v[63] = cpoint[0] + v[9] * cos(v[2]);
	v[15] = QAi[0][0];
	v[16] = QAi[0][1];
	v[17] = QAi[0][2];
	v[18] = QAi[1][0];
	v[19] = QAi[1][1];
	v[20] = QAi[1][2];
	v[21] = QAi[2][0];
	v[22] = QAi[2][1];
	v[23] = QAi[2][2];
	v[33] = 4e0 / (4e0 + v[41] + v[70]);
	v[71] = -0.5e0*v[33];
	v[36] = 1e0 + v[70] * v[71];
	v[37] = v[33] * (v[39] - v[8]);
	v[38] = v[33] * (v[44] + v[7]);
	v[40] = v[33] * (v[39] + v[8]);
	v[42] = 1e0 + (v[35] + v[41])*v[71];
	v[43] = v[33] * (v[46] - v[6]);
	v[45] = v[33] * (v[44] - v[7]);
	v[47] = v[33] * (v[46] + v[6]);
	v[48] = 1e0 - (-v[34] - v[41])*v[71];
	v[61] = v[63] * cos(v[1]);
	v[62] = cpoint[1] + v[9] * sin(v[2]);
	v[64] = -(v[63] * sin(v[1]));
	Gam[0] = dof[0] + (v[15] * v[36] + v[18] * v[37] + v[21] * v[38])*v[61] + (v[16] * v[36] + v[19] * v[37] + v[22] * v[38]
		)*v[62] + (v[17] * v[36] + v[20] * v[37] + v[23] * v[38])*v[64] + xAi[0];
	Gam[1] = dof[1] + (v[15] * v[40] + v[18] * v[42] + v[21] * v[43])*v[61] + (v[16] * v[40] + v[19] * v[42] + v[22] * v[43]
		)*v[62] + (v[17] * v[40] + v[20] * v[42] + v[23] * v[43])*v[64] + xAi[1];
	Gam[2] = dof[2] + (v[15] * v[45] + v[18] * v[47] + v[21] * v[48])*v[61] + (v[16] * v[45] + v[19] * v[47] + v[22] * v[48]
		)*v[62] + (v[17] * v[45] + v[20] * v[47] + v[23] * v[48])*v[64] + xAi[2];

	//Copying final result
	(*Gamma)(0, 0) = Gam[0];
	(*Gamma)(1, 0) = Gam[1];
	(*Gamma)(2, 0) = Gam[2];

}

void ArcRevolution::NormalExt(double z, double th, Matrix* n, bool next)
{
	double v[210];
	double c[2];
	c[0] = z;
	c[1] = th;
	double* rad = radius;
	double* cpoint = c_point->getMatrix();

	double dof[6];
	if (next)
	{
		for (int i = 0; i < 6; i++)
			dof[i] = (*d)(i, 0);
	}
	else
	{
		for (int i = 0; i < 6; i++)
			dof[i] = 0.0;
	}
	double normal[3];
	bool* normalint = &flag_normal_int;

	//AceGen code
	v[1] = c[0];
	v[77] = sin(v[1]);
	v[75] = cos(v[1]);
	v[2] = c[1];
	v[6] = dof[3];
	v[42] = (v[6] * v[6]);
	v[7] = dof[4];
	v[92] = v[7] / 2e0;
	v[40] = v[6] * v[92];
	v[35] = (v[7] * v[7]);
	v[8] = dof[5];
	v[47] = v[8] * v[92];
	v[45] = (v[6] * v[8]) / 2e0;
	v[36] = (v[8] * v[8]);
	v[93] = v[35] + v[36];
	v[9] = (*rad);
	v[79] = v[9] * cos(v[2]);
	v[74] = v[9] * sin(v[2]);
	v[78] = v[74] * v[77];
	v[76] = -(v[74] * v[75]);
	v[64] = cpoint[0] + v[79];
	v[70] = v[64] * v[75];
	v[65] = -(v[64] * v[77]);
	v[15] = QAi[0][0];
	v[16] = QAi[0][1];
	v[17] = QAi[0][2];
	v[18] = QAi[1][0];
	v[19] = QAi[1][1];
	v[20] = QAi[1][2];
	v[21] = QAi[2][0];
	v[22] = QAi[2][1];
	v[23] = QAi[2][2];
	v[34] = 4e0 / (4e0 + v[42] + v[93]);
	v[94] = -0.5e0*v[34];
	v[37] = 1e0 + v[93] * v[94];
	v[38] = v[34] * (v[40] - v[8]);
	v[39] = v[34] * (v[45] + v[7]);
	v[41] = v[34] * (v[40] + v[8]);
	v[43] = 1e0 + (v[36] + v[42])*v[94];
	v[44] = v[34] * (v[47] - v[6]);
	v[46] = v[34] * (v[45] - v[7]);
	v[48] = v[34] * (v[47] + v[6]);
	v[49] = 1e0 - (-v[35] - v[42])*v[94];
	v[50] = v[15] * v[37] + v[18] * v[38] + v[21] * v[39];
	v[51] = v[16] * v[37] + v[19] * v[38] + v[22] * v[39];
	v[52] = v[17] * v[37] + v[20] * v[38] + v[23] * v[39];
	v[80] = v[50] * v[76] + v[52] * v[78] + v[51] * v[79];
	v[71] = v[50] * v[65] - v[52] * v[70];
	v[53] = v[15] * v[41] + v[18] * v[43] + v[21] * v[44];
	v[54] = v[16] * v[41] + v[19] * v[43] + v[22] * v[44];
	v[55] = v[17] * v[41] + v[20] * v[43] + v[23] * v[44];
	v[81] = v[53] * v[76] + v[55] * v[78] + v[54] * v[79];
	v[72] = v[53] * v[65] - v[55] * v[70];
	v[89] = -(v[72] * v[80]) + v[71] * v[81];
	v[56] = v[15] * v[46] + v[18] * v[48] + v[21] * v[49];
	v[57] = v[16] * v[46] + v[19] * v[48] + v[22] * v[49];
	v[58] = v[17] * v[46] + v[20] * v[48] + v[23] * v[49];
	v[82] = v[56] * v[76] + v[58] * v[78] + v[57] * v[79];
	v[73] = v[56] * v[65] - v[58] * v[70];
	v[87] = v[73] * v[80] - v[71] * v[82];
	v[63] = cpoint[1] + v[74];
	v[84] = -(v[73] * v[81]) + v[72] * v[82];
	v[95] = ((*normalint) ? -1 : 1) / sqrt((v[84] * v[84]) + (v[87] * v[87]) + (v[89] * v[89]));
	normal[0] = v[84] * v[95];
	normal[1] = v[87] * v[95];
	normal[2] = v[89] * v[95];

	//Copying final result
	(*n)(0, 0) = normal[0];
	(*n)(1, 0) = normal[1];
	(*n)(2, 0) = normal[2];
}
