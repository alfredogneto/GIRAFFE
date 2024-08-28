#include "RigidNURBS_1.h"

#include "NURBSSurface.h"
#include "PostFiles.h"
#include "Node.h"
#include "CoordinateSystem.h"

#include"Database.h"
//Variaveis globais
extern
Database db;

RigidNURBS_1::RigidNURBS_1()
{
	nDOFs = 6;
	pilot_node = 1;
	pilot_is_used = true;
	n_nodes = 1;
	invert_normal = false;
	CADData_ID = 0;
	nodes = new int[n_nodes];
	DOFs = new int *[n_nodes];
	for (int i = 0; i < n_nodes; i++)
		DOFs[i] = new int[db.number_GLs_node];
	
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

	GLs = new int*[nDOFs];
	for (int i = 0; i < nDOFs; i++)
		GLs[i] = NULL;

	Qi = Matrix(3, 3);
	Qip= Matrix(3, 3);
	xi = Matrix(3);
	xip= Matrix(3);
	I3 = Matrix(3, 3);
	I3(0, 0) = 1.0;
	I3(1, 1) = 1.0;
	I3(2, 2) = 1.0;

	InitializeDegeneration();
}

RigidNURBS_1::~RigidNURBS_1()
{
	delete[] nodes;
	delete[]GLs;
	if (DOFs != NULL)
	{
		for (int i = 0; i < n_nodes; i++)
			delete[] DOFs[i];
		delete[] DOFs;
	}

	FreeDegeneration();
	
}

void RigidNURBS_1::SetMinMaxRange()
{
	u1_min = surf->U_knot_vector[surf->U_order];
	u1_max = surf->U_knot_vector[surf->U_order + surf->U_dim];
	u1_range = u1_max - u1_min;

	u2_min = surf->V_knot_vector[surf->V_order];
	u2_max = surf->V_knot_vector[surf->V_order + surf->V_dim];
	u2_range = u2_max - u2_min;
}

//Normal exterior a superficie na posição escolhida
void RigidNURBS_1::NormalExt(double* zeta, double* theta, Matrix* n)
{
	Matrix** data;
	int d = 1;
	data = new Matrix*[d + 1];
	for (int i = 0; i < d + 1; i++)
	{
		data[i] = new Matrix[d + 1];
		for (int j = 0; j < d + 1; j++)
			data[i][j] = Matrix(3);
	}
	//material point (zeta,theta)
	surf->NURBSDerivatives(*zeta, *theta, data, d);
	Matrix t1 = 1.0 / norm(data[1][0])*data[1][0];
	Matrix t2 = 1.0 / norm(data[0][1])*data[0][1];
	*n = 1.0 / norm(cross(t1, t2))*cross(t1, t2);
	if (invert_normal == true)
		*n = -1.0* (*n);
	*n = Qip * (*n);
	for (int i = 0; i < d + 1; i++)
		delete[] data[i];
	delete[]data;

} 

//Obtem ponto da superficie
void RigidNURBS_1::SurfacePoint(double& zeta, double& theta, Matrix& point)
{
	Matrix p(3);
	surf->NURBSPoint(zeta, theta, p);
	point = xip + Qip * p;
}

//Atualiza bounding box
void RigidNURBS_1::UpdateBox()
{
	//Large value of bounding box - encompassing the whole model
	double diag = db.EvaluateBoundingBoxDiag();
	//Bounding box
	double x[8], y[8], z[8];
	x[0] = -diag;
	y[0] = -diag;
	z[0] = -diag;

	x[1] = +diag;
	y[1] = -diag;
	z[1] = -diag;

	x[2] = +diag;
	y[2] = +diag;
	z[2] = -diag;

	x[3] = -diag;
	y[3] = +diag;
	z[3] = -diag;
	
	x[4] = -diag;
	y[4] = -diag;
	z[4] = +diag;

	x[5] = +diag;
	y[5] = -diag;
	z[5] = +diag;

	x[6] = +diag;
	y[6] = +diag;
	z[6] = +diag;

	x[7] = -diag;
	y[7] = +diag;
	z[7] = +diag;

	//Setando os vertices
	box.SetVertices(x, y, z);
}

void RigidNURBS_1::WriteVTK_XMLRender(FILE *f)
{
	if (db.post_files->WriteRigidContactSurfaces_flag == true)
	{
		Matrix pilot(3);
		//Posição de cada ponto P no plano xy (referência)
		pilot(0, 0) = db.post_files->mag_factor*(db.nodes[pilot_node - 1]->copy_coordinates[0] - db.nodes[pilot_node - 1]->ref_coordinates[0]) + db.nodes[pilot_node - 1]->ref_coordinates[0];
		pilot(1, 0) = db.post_files->mag_factor*(db.nodes[pilot_node - 1]->copy_coordinates[1] - db.nodes[pilot_node - 1]->ref_coordinates[1]) + db.nodes[pilot_node - 1]->ref_coordinates[1];
		pilot(2, 0) = db.post_files->mag_factor*(db.nodes[pilot_node - 1]->copy_coordinates[2] - db.nodes[pilot_node - 1]->ref_coordinates[2]) + db.nodes[pilot_node - 1]->ref_coordinates[2];

		surf->WriteVTK_XMLRender(f, pilot, Qi, number);
	}
}

bool RigidNURBS_1::Read(FILE *f)
{
	char s[1000];
	fscanf(f, "%s", s);
	number = atoi(s);

	fscanf(f, "%s", s);
	if (!strcmp(s, "CS"))
	{
		fscanf(f, "%s", s);
		number_CS = atoi(s);
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "PilotNode"))
	{
		fscanf(f, "%s", s);
		pilot_node = atoi(s);
		nodes[0] = pilot_node;
	}
	else
		return false;

	fscanf(f, "%s", s);
	if (!strcmp(s, "CADData"))
	{
		fscanf(f, "%s", s);
		CADData_ID = atoi(s);
	}
	else
		return false;

	//Flag normal invertida
	//Salva a posição (stream)
	fpos_t pos;
	fgetpos(f, &pos);

	fscanf(f, "%s", s);
	if (!strcmp(s, "ExteriorNormalInverted"))
		invert_normal = true;
	else
	{
		fsetpos(f, &pos);
		invert_normal = false;
	}

	ReadCommon(f);

	return true;
}

void RigidNURBS_1::Write(FILE *f)
{
	char s[30];
	if (invert_normal == true)
		sprintf(s, "ExteriorNormalInverted");
	else
		sprintf(s, "");
	fprintf(f, "RigidNURBS_1\t%d\tCS\t%d\tPilotNode\t%d\tCADData\t%d\t%s\n", number, number_CS, pilot_node, CADData_ID, s);
}

//Checa inconsistências para evitar erros de execução
bool RigidNURBS_1::Check()
{
	if (pilot_node > db.number_nodes)
		return false;
	if (number_CS > db.number_CS)
		return false;
	if (CADData_ID > db.number_cad_data)
		return false;
	if (typeid(*db.cad_data[CADData_ID - 1]) != typeid(NURBSSurface))
		return false;
	return true;
}

//Realiza chute inicial para as variaveis zeta e theta
void RigidNURBS_1::InitialGuess(Matrix* xS, double** convective, int n_solutions)
{
	
}

void RigidNURBS_1::PreCalc()
{
	surf = static_cast<NURBSSurface*>(db.cad_data[CADData_ID - 1]);
	
	//Apontando para posição que indica valor dos GLs globais
	for (int i = 0; i < 6; i++)
		GLs[i] = &db.nodes[pilot_node - 1]->GLs[i];
	//Matriz para transformação de coordenadas  local-global
	Qi(0, 0) = (*db.CS[number_CS - 1]->E1)(0, 0);
	Qi(1, 0) = (*db.CS[number_CS - 1]->E1)(1, 0);
	Qi(2, 0) = (*db.CS[number_CS - 1]->E1)(2, 0);
	Qi(0, 1) = (*db.CS[number_CS - 1]->E2)(0, 0);
	Qi(1, 1) = (*db.CS[number_CS - 1]->E2)(1, 0);
	Qi(2, 1) = (*db.CS[number_CS - 1]->E2)(2, 0);
	Qi(0, 2) = (*db.CS[number_CS - 1]->E3)(0, 0);
	Qi(1, 2) = (*db.CS[number_CS - 1]->E3)(1, 0);
	Qi(2, 2) = (*db.CS[number_CS - 1]->E3)(2, 0);
	//Position - pilot node
	xi(0, 0) = db.nodes[pilot_node - 1]->ref_coordinates[0];
	xi(1, 0) = db.nodes[pilot_node - 1]->ref_coordinates[1];
	xi(2, 0) = db.nodes[pilot_node - 1]->ref_coordinates[2];

	DegenerationPreCalc();
}

//Retorna as coordenadas da superficie para um par (zeta,theta) - configuração anterior convergida
void RigidNURBS_1::Gamma_and_Triad(Matrix* G_p, Matrix* t1_p, Matrix* t2_p, Matrix* n_p, Matrix* G_i, Matrix* t1_i, Matrix* t2_i, Matrix* n_i, Matrix* G_ip, double* zi, double* thi, double* zp, double* thp)
{
	Matrix** data;
	int d = 1;
	data = new Matrix*[d + 1];
	for (int i = 0; i < d + 1; i++)
	{
		data[i] = new Matrix[d + 1];
		for (int j = 0; j < d + 1; j++)
			data[i][j] = Matrix(3);
	}
	//material point 'i'
	surf->NURBSDerivatives(*zi, *thi, data, d);
	Matrix Gi = data[0][0];
	Matrix t1i = 1.0 / norm(data[1][0])*data[1][0];
	Matrix t2i = 1.0 / norm(data[0][1])*data[0][1];
	Matrix ni = 1.0 / norm(cross(t1i, t2i))*cross(t1i, t2i);
	if (invert_normal)
		ni = -1.0*ni;
	//material point 'p'
	surf->NURBSDerivatives(*zp, *thp, data, d);
	Matrix Gp = data[0][0];
	Matrix t1p = 1.0 / norm(data[1][0])*data[1][0];
	Matrix t2p = 1.0 / norm(data[0][1])*data[0][1];
	Matrix np = 1.0 / norm(cross(t1p, t2p))*cross(t1i, t2p);
	if (invert_normal)
		np = -1.0*np;
	//Retornos
	*G_i = xi + Qi * Gi;
	*G_p = xip + Qip * Gp;
	*G_ip = xip + Qip * Gi;

	*t1_p = Qip * t1p;
	*t2_p = Qip * t2p;
	*n_p = Qip * np;

	*t1_i = Qi * t1i;
	*t2_i = Qi * t2i;
	*n_i = Qi * ni;

	for (int i = 0; i < d + 1; i++)
		delete[] data[i];
	delete[]data;
	return;
}

//Dado o ponto xS, calcula as coordenadas (zeta,theta) referentes a minima distancia
void RigidNURBS_1::FindMinimimumParameters(Matrix* xS, NSContactData* cd)
{
	
}
//Atualiza as variaveis internas da superficie, para pegarem info do pilot node para uso posterior com posição atualizada
void RigidNURBS_1::FillNodes()
{
	//Updating pilot node position and rotation tensor
	Matrix disp(3);
	disp(0, 0) = db.nodes[pilot_node - 1]->displacements[0];
	disp(1, 0) = db.nodes[pilot_node - 1]->displacements[1];
	disp(2, 0) = db.nodes[pilot_node - 1]->displacements[2];
	xip = xi + disp;
	Matrix alpha_1(3);
	alpha_1(0, 0) = db.nodes[pilot_node - 1]->displacements[3];
	alpha_1(1, 0) = db.nodes[pilot_node - 1]->displacements[4];
	alpha_1(2, 0) = db.nodes[pilot_node - 1]->displacements[5];
	double alpha_escalar = norm(alpha_1);
	Matrix A = skew(alpha_1);
	double g = 4.0 / (4.0 + alpha_escalar * alpha_escalar);
	Matrix Q = I3 + g * (A + 0.5*((A)*(A)));
	Qip = Q * Qi;
}

//Retorna coordenadas globais do ponto central da superficie a ser utilizado para calculos grosseiros de sua localização (pinball)
void RigidNURBS_1::CenterPoint(Matrix* center)
{
	double avg_u = surf->U_knot_vector[(surf->U_order + surf->U_dim) / 2];
	double avg_v = surf->V_knot_vector[(surf->V_order + surf->V_dim) / 2];
	Matrix local(3);
	surf->NURBSPoint(avg_u, avg_v, local);
	*center = xi + Qi * local;
}

//Salva vetores de configuração convergida
void RigidNURBS_1::SaveConfiguration()
{
	Qi = Qip;
	xi = xip;
}

//Calcula contribuições de contato entre esfera e superficie
void RigidNURBS_1::ContactSphereSurfaceSticking(double* Rc, double** Kc, double zetap, double thetap, double zetai, double thetai, double* gti, int node, double* epsn, double* epst, double* cn, double* ct, double* mu, double* radius)
{
	
}

//Calcula contribuições de contato entre esfera e superficie
void RigidNURBS_1::ContactSphereSurfaceSliding(double* Rc, double** Kc, double zetap, double thetap, double zetai, double thetai, double* gti, int node, double* epsn, double* epst, double* cn, double* ct, double* mu, double* radius)
{
	
}