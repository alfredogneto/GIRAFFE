#include "FlexibleArcExtrusion_1_RigidArcRevolution_1.h"

#include "Dynamic.h"
#include "FlexibleArcExtrusion_1.h"
#include "RigidArcRevolution_1.h"
#include "ArcCirc.h"
#include "Node.h"
#include "Matrix.h"
#include "SSContactData.h"

#include "Database.h"
//Variaveis globais
extern
Database db;
#define PI 3.1415926535897932384626433832795

FlexibleArcExtrusion_1_RigidArcRevolution_1::FlexibleArcExtrusion_1_RigidArcRevolution_1()
{
	DefaultValues();
}

FlexibleArcExtrusion_1_RigidArcRevolution_1::~FlexibleArcExtrusion_1_RigidArcRevolution_1()
{
	Free();
}

void FlexibleArcExtrusion_1_RigidArcRevolution_1::InitializeConvectiveRange()
{
	FlexibleArcExtrusion_1* surf1;		//Ponteiro para a superficie 1
	RigidArcRevolution_1* surf2;		//Ponteiro para a superficie 2
	surf1 = static_cast<FlexibleArcExtrusion_1*>(db.surfaces[surf1_ID - 1]);
	surf2 = static_cast<RigidArcRevolution_1*>(db.surfaces[surf2_ID - 1]);

	double thetai1 = db.arcs[surf1->arc_ID - 1]->theta_i;
	double thetaf1 = db.arcs[surf1->arc_ID - 1]->theta_f;
	double phii = 0.0;
	double phif = surf2->rev_Ang;
	double thetai2 = db.arcs[surf2->arc_ID - 1]->theta_i;
	double thetaf2 = db.arcs[surf2->arc_ID - 1]->theta_f;

	//Variaveis de medicao de range de coordenadas convectivas
	double range_zeta;
	double range_phi;
	double range_theta1;
	double range_theta2;

	
	range_zeta = 2.0;
	if (thetaf1 >= thetai1)
		range_theta1 = abs(thetaf1 - thetai1);
	else
		range_theta1 = 2 * PI - abs(thetaf1 - thetai1);
	range_phi = abs(phif - phii);
	if (thetaf2 >= thetai2)
		range_theta2 = abs(thetaf2 - thetai2);
	else
		range_theta2 = 2 * PI - abs(thetaf2 - thetai2);

	convective_min(0, 0) = -1;
	convective_max(0, 0) = +1;
	//Degeneration
	if (surf1->degeneration[0] == true)
	{
		convective_range(0, 0) = 0.0;
	}
	//No degeneration
	else
	{
		convective_range(0, 0) = range_zeta;
	}

	convective_min(1, 0) = thetai1;
	convective_max(1, 0) = thetaf1;
	//Degeneration
	if (surf1->degeneration[1] == true)
	{
		convective_range(1, 0) = 0.0;
	}
	//No degeneration
	else
	{
		convective_range(1, 0) = range_theta1;
	}

	convective_min(2, 0) = phii;
	convective_max(2, 0) = phif;
	//Degeneration
	if (surf2->degeneration[0] == true)
	{
		convective_range(2, 0) = 0.0;
	}
	//No degeneration
	else
	{
		convective_range(2, 0) = range_phi;
	}

	convective_min(3, 0) = thetai2;
	convective_max(3, 0) = thetaf2;
	//Degeneration
	if (surf2->degeneration[1] == true)
	{
		convective_range(3, 0) = 0.0;
	}
	//No degeneration
	else
	{
		convective_range(3, 0) = range_theta2;
	}
}

bool FlexibleArcExtrusion_1_RigidArcRevolution_1::SpecialLCP(Matrix& solution)
{
	///////////////////////////////Ponteiros e variaveis para facilitar acesso//////////////////////////////////////////
	FlexibleArcExtrusion_1* surf_1 = static_cast<FlexibleArcExtrusion_1*>(db.surfaces[surf1_ID - 1]);
	int arc1id = surf_1->arc_ID;
	RigidArcRevolution_1* surf_2 = static_cast<RigidArcRevolution_1*>(db.surfaces[surf2_ID - 1]);
	int arc2id = surf_2->arc_ID;
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	// Pontos dos arcos nos sistemas locais
	Matrix c1(3);	// Centro extrudada - trilho
	c1(0, 0) = db.arcs[arc1id - 1]->c_point(0, 0);
	c1(1, 0) = db.arcs[arc1id - 1]->c_point(1, 0);
	c1(2, 0) = 0.0;
	Matrix c2(3);	// Centro revolucionada - roda
	c2(0, 0) = db.arcs[arc2id - 1]->c_point(0, 0);
	c2(1, 0) = db.arcs[arc2id - 1]->c_point(1, 0);
	c2(2, 0) = 0.0;
	Matrix i1(3);	// Inicial extrudada - trilho
	i1(0, 0) = db.arcs[arc1id - 1]->i_point(0, 0);
	i1(1, 0) = db.arcs[arc1id - 1]->i_point(1, 0);
	i1(2, 0) = 0.0;
	Matrix f1(3);	// Final extrudada - trilho
	f1(0, 0) = db.arcs[arc1id - 1]->f_point(0, 0);
	f1(1, 0) = db.arcs[arc1id - 1]->f_point(1, 0);
	f1(2, 0) = 0.0;
	Matrix i2(3);	// Inicial revolucionada - roda
	i2(0, 0) = db.arcs[arc2id - 1]->i_point(0, 0);
	i2(1, 0) = db.arcs[arc2id - 1]->i_point(1, 0);
	i2(2, 0) = 0.0;
	Matrix f2(3);	// Final revolucionada - roda
	f2(0, 0) = db.arcs[arc2id - 1]->f_point(0, 0);
	f2(1, 0) = db.arcs[arc2id - 1]->f_point(1, 0);
	f2(2, 0) = 0.0;


	//	Matrizes de Rotação
	//	Extrudada - Trilho
	Matrix alpha_AA_c(3);
	Matrix xAA_c(3);
	for (int i = 0; i < 3; i++)
	{
		alpha_AA_c(i, 0) = db.nodes[surf_1->nodes[0] - 1]->displacements[i + 3];
		xAA_c(i, 0) = db.nodes[surf_1->nodes[0] - 1]->copy_coordinates[i] + db.nodes[surf_1->nodes[0] - 1]->displacements[i];
	}
	//	Revolucionada - Roda
	Matrix alpha_AB_c(3);
	Matrix xAB_c(3);
	for (int i = 0; i < 3; i++)
	{
		alpha_AB_c(i, 0) = db.nodes[surf_2->nodes[0] - 1]->displacements[i + 3];
		xAB_c(i, 0) = db.nodes[surf_2->nodes[0] - 1]->copy_coordinates[i] + db.nodes[surf_2->nodes[0] - 1]->displacements[i];
	}
	//	Q_AAi - Extrudada - Trilho
	Matrix I3(3, 3);
	I3(0, 0) = 1.0;
	I3(1, 1) = 1.0;
	I3(2, 2) = 1.0;
	double alpha = norm(alpha_AA_c);							//Valor escalar do parametro alpha
	Matrix A = skew(alpha_AA_c);								//Matriz A
	double g = 4.0 / (4.0 + alpha * alpha);					//função g(alpha) - em algumas ref. tb. chamado de h(alpha)
	Matrix Q_AAd = I3 + g * (A + 0.5*(A*A));						//Tensor de rotação
	Matrix Q1 = Q_AAd * (*surf_1->Q_AAi);
	//	Q_ABi - revolucionada - Roda
	alpha = norm(alpha_AB_c);							//Valor escalar do parametro alpha
	A = skew(alpha_AB_c);								//Matriz A
	g = 4.0 / (4.0 + alpha * alpha);					//função g(alpha) - em algumas ref. tb. chamado de h(alpha)
	Matrix Q_ABd = I3 + g * (A + 0.5*(A*A));						//Tensor de rotação
	Matrix Q2 = Q_ABd * (*surf_2->Q_ABi);

	xAA_c(0, 0) = 0.0;
	xAB_c(0, 0) = 0.0;
	// Sistema local para Global
	Matrix c1g(3);
	Matrix c2g(3);
	Matrix i1g(3);
	Matrix f1g(3);
	Matrix i2g(3);
	Matrix f2g(3);
	c1g = Q1 * c1 + xAA_c;
	c2g = Q2 * c2 + xAB_c;
	i1g = Q1 * i1 + xAA_c;
	f1g = Q1 * f1 + xAA_c;
	i2g = Q2 * i2 + xAB_c;
	f2g = Q2 * f2 + xAB_c;

	// Coordenadas convectivas que não vamos usar
	Matrix xk(4);
	xk(0, 0) = 0.0;		
	xk(2, 0) = 0.0;		

	// Valores de theta no sistema global para achar o ponto de interseção
	double theta_glob1;
	double theta_glob2;
	theta_glob1 = atan2((c2g(2, 0) - c1g(2, 0)), (c2g(1, 0) - c1g(1, 0)));		// Extrudada - Trilho
	theta_glob2 = atan2((c1g(2, 0) - c2g(2, 0)), (c1g(1, 0) - c2g(1, 0)));		// Revolucionada - Roda

	////	Verifica se o theta esta na direção correta ou oposta
	//if (dot(0.5*(i1g + f1g) - c1g, c2g - c1g) < 0)	// Extrudada - Trilho
	//{
	//	theta_glob1 = theta_glob1 + PI;
	//}
	//if (dot(0.5*(i2g + f2g) - c2g, c1g - c2g) < 0)	//Revolucionada - Roda
	//{
	//	theta_glob2 = theta_glob2 + PI;
	//}

	////	Pontos da interseção no sistema global
	//Matrix p1g(3);	// Extrudada - Trilho
	//p1g(0, 0) = 0;
	//p1g(1, 0) = c1g(1, 0) + (db.arcs[arc1id - 1]->radius)*cos(theta_glob1);
	//p1g(2, 0) = c1g(2, 0) + (db.arcs[arc1id - 1]->radius)*sin(theta_glob1);
	//Matrix p2g(3);	// Revolucionada - Roda
	//p2g(0, 0) = 0;
	//p2g(1, 0) = c2g(1, 0) + (db.arcs[arc2id - 1]->radius)*cos(theta_glob2);
	//p2g(2, 0) = c2g(2, 0) + (db.arcs[arc2id - 1]->radius)*sin(theta_glob2);


	//	Pontos da interseção no sistema global
	Matrix p1g(3);	// Extrudada - Trilho
	p1g(0, 0) = 0;
	p1g(1, 0) = c1g(1, 0) + (db.arcs[arc1id - 1]->radius)*cos(theta_glob1);
	p1g(2, 0) = c1g(2, 0) + (db.arcs[arc1id - 1]->radius)*sin(theta_glob1);
	Matrix p2g(3);	// Revolucionada - Roda
	p2g(0, 0) = 0;
	p2g(1, 0) = c2g(1, 0) + (db.arcs[arc2id - 1]->radius)*cos(theta_glob2);
	p2g(2, 0) = c2g(2, 0) + (db.arcs[arc2id - 1]->radius)*sin(theta_glob2);

	
	//	Verifica se o theta esta na direção correta ou oposta
	if (dot(0.5*(i1g + f1g) - c1g, p1g - c1g) < 0)	// Extrudada - Trilho
	{
		theta_glob1 = theta_glob1 + PI;

		p1g(1, 0) = c1g(1, 0) + (db.arcs[arc1id - 1]->radius)*cos(theta_glob1);
		p1g(2, 0) = c1g(2, 0) + (db.arcs[arc1id - 1]->radius)*sin(theta_glob1);
	}
	if (dot(0.5*(i2g + f2g) - c2g, p2g - c2g) < 0)	//Revolucionada - Roda
	{
		theta_glob2 = theta_glob2 + PI;

		p2g(1, 0) = c2g(1, 0) + (db.arcs[arc2id - 1]->radius)*cos(theta_glob2);
		p2g(2, 0) = c2g(2, 0) + (db.arcs[arc2id - 1]->radius)*sin(theta_glob2);
	}


	//	Pontos da interseção no sistema local
	Matrix p1(3);	// Extrudada - Trilho
	Matrix p2(3);	// Revolucionada - Roda
	p1 = transp(Q1)*(p1g - xAA_c);
	p2 = transp(Q2)*(p2g - xAB_c);

	// Coordenada convectivas (thetas no sistema local)
	xk(1, 0) = atan2(p1(1, 0) - c1(1, 0), p1(0, 0) - c1(0, 0));	// Extrudada - Trilho
	xk(3, 0) = atan2(p2(1, 0) - c2(1, 0), p2(0, 0) - c2(0, 0));	// Revolucionada - Roda


	////////////////////
	//xk.print();

	for (int i = 0; i < 4; i++)
		(solution)(i, 0) = xk(i, 0);

	return true;
}

//Verifica range de coordenadas convectivas
int FlexibleArcExtrusion_1_RigidArcRevolution_1::VerifyConvectiveRange(Matrix& mc)
{
	//Retornos:
	//0 - Range fisico da superficie
	//4 - Fora do range fisico da superficie - proximo
	//2 - Fora do range fisico da superficie - distante (nao strong) 
	FlexibleArcExtrusion_1* surf1;		//Ponteiro para a superficie 1
	RigidArcRevolution_1* surf2;		//Ponteiro para a superficie 2
	surf1 = static_cast<FlexibleArcExtrusion_1*>(db.surfaces[surf1_ID - 1]);
	surf2 = static_cast<RigidArcRevolution_1*>(db.surfaces[surf2_ID - 1]);
	double thetai1 = db.arcs[surf1->arc_ID - 1]->theta_i;
	double thetaf1 = db.arcs[surf1->arc_ID - 1]->theta_f;
	double phii = 0.0;
	double phif = surf2->rev_Ang;
	double thetai2 = db.arcs[surf2->arc_ID - 1]->theta_i;
	double thetaf2 = db.arcs[surf2->arc_ID - 1]->theta_f;

	//Variaveis de medicao de range de coordenadas convectivas
	double range_zeta;
	double range_phi;
	double range_theta1;
	double range_theta2;

	int return_vector[4];
	for (int i = 0; i < 4; i++)
		return_vector[i] = 0;
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//Teste 1 - zeta - superficie 1
	range_zeta = perc * convective_range(0, 0);
	if (abs(mc(0, 0)) <= 1.00)
		return_vector[0] = 0;
	else
	{
		if (abs(mc(0, 0)) <= (1.00+range_zeta))
			return_vector[0] = 4;
		else
			return_vector[0] = 2;
	}
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//Teste 2 - theta_1 - superficie 1
	range_theta1 = perc * convective_range(1, 0);
	if (thetaf1 >= thetai1)
	{
		if (ArcReduction(mc(1, 0)) >= thetai1 && ArcReduction(mc(1, 0)) <= thetaf1)
			return_vector[1] = 0;
		else
		{
			if (ArcReduction(mc(1, 0)) >= (thetai1 - range_theta1) && ArcReduction(mc(1, 0)) <= (thetaf1 + range_theta1))
				return_vector[1] = 4;
			else
				return_vector[1] = 2;
		}
	}
	else
	{
		//Verificacao do complemento
		if (ArcReduction(mc(1, 0)) >= thetaf1 && ArcReduction(mc(1, 0)) <= thetai1)
		{
			if (ArcReduction(mc(1, 0)) >= (thetaf1 + range_theta1) && ArcReduction(mc(1, 0)) <= (thetai1 - range_theta1))
				return_vector[1] = 2;
			else
				return_vector[1] = 4;
		}
		else
			return_vector[1] = 0;
	}
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//Teste 3 - phi - superficie 2
	range_phi = perc * convective_range(2, 0);
	if (ArcReduction2p(mc(2, 0)) >= phii && ArcReduction2p(mc(2, 0)) <= phif)
		return_vector[2] = 0;
	else
	{
		if (ArcReduction2p(mc(2, 0)) >= (phii - range_phi) && ArcReduction2p(mc(2, 0)) <= (phif + range_phi))
			return_vector[2] = 4;
		else
			return_vector[2] = 2;
	}
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//Teste 4 - theta_2 - superficie 2
	range_theta2 = perc * convective_range(3, 0);
	if (thetaf2 >= thetai2)
	{
		if (ArcReduction(mc(3, 0)) >= thetai2 && ArcReduction(mc(3, 0)) <= thetaf2)
			return_vector[3] = 0;
		else
		{
			if (ArcReduction(mc(3, 0)) >= (thetai2 - range_theta2) && ArcReduction(mc(3, 0)) <= (thetaf2 + range_theta2))
				return_vector[3] = 4;
			else
				return_vector[3] = 2;
		}
	}
	else
	{
		//Verificacao do complemento
		if (ArcReduction(mc(3, 0)) >= thetaf2 && ArcReduction(mc(3, 0)) <= thetai2)
		{
			if (ArcReduction(mc(3, 0)) >= (thetaf2 + range_theta2) && ArcReduction(mc(3, 0)) <= (thetai2 - range_theta2))
				return_vector[3] = 2;
			else
				return_vector[3] = 4;
		}
		else
			return_vector[3] = 0;
	}
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//Varredura final para retorno
	int return_value = 0;
	for (int i = 0; i < 4; i++)
	{
		if (return_value == 0)
			return_value = return_vector[i];
		else
		{
			if (return_value == 4 && return_vector[i] == 2)
				return_value = 2;
		}
	}
	
	return return_value;
}

//Calcula e rotorna o gap (com sinal)
double FlexibleArcExtrusion_1_RigidArcRevolution_1::Gap(Matrix& mc, bool fixed_normals, Matrix& nA, Matrix& nB)
{
	double v[2000];		//variavel temporaria - AceGen

	bool* fixnormal = &fixed_normals;
	double* normalA = nA.getMatrix();
	double* normalB = nB.getMatrix();

	///////////////////////////////Ponteiros e variaveis para facilitar acesso//////////////////////////////////////////
	FlexibleArcExtrusion_1* surf1;		//Ponteiro para a superficie 1
	RigidArcRevolution_1* surf2;		//Ponteiro para a superficie 2
	double* radA;
	double* radB;
	double* xfac;
	double* zfac;
	double* cpointA;
	double* cpointB;
	double* dA;
	double* dB;
	double* duiA;
	double* dduiA;
	double* duiB;
	double* dduiB;
	bool* normalintA;
	bool* normalintB;
	double* xAAi;
	double* xBAi;
	double* xABi;
	double** QAAi;
	double** QBAi;
	double** QABi;
	surf1 = static_cast<FlexibleArcExtrusion_1*>(db.surfaces[surf1_ID - 1]);
	surf2 = static_cast<RigidArcRevolution_1*>(db.surfaces[surf2_ID - 1]);
	radA = surf1->radius;
	radB = surf2->radius;
	xfac = &surf2->x_fac;
	zfac = &surf2->z_fac;
	cpointA = surf1->c_point->getMatrix();
	cpointB = surf2->c_point->getMatrix();
	dA = surf1->d_A->getMatrix();
	dB = surf2->d_B->getMatrix();
	duiA = surf1->dui_A->getMatrix();
	dduiA = surf1->ddui_A->getMatrix();
	duiB = surf2->dui_B->getMatrix();
	dduiB = surf2->ddui_B->getMatrix();
	xAAi = surf1->x_AAi->getMatrix();
	xBAi = surf1->x_BAi->getMatrix();
	xABi = surf2->x_ABi->getMatrix();
	normalintA = &surf1->flag_normal_int;
	normalintB = &surf2->flag_normal_int;
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	QAAi = new double*[3];
	QBAi = new double*[3];
	QABi = new double*[3];
	for (int i = 0; i < 3; i++)
	{
		QAAi[i] = new double[3];
		QBAi[i] = new double[3];
		QABi[i] = new double[3];
	}
	//Salvando variaveis locais para montagem de superficies
	surf1->Q_AAi->MatrixToPtr(QAAi, 3);
	surf1->Q_BAi->MatrixToPtr(QBAi, 3);
	surf2->Q_ABi->MatrixToPtr(QABi, 3);
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/*surf1->Q_AAi->print();
	surf1->Q_BAi->print();
	surf2->Q_ABi->print();*/
	double gap;
	double *c = mc.getMatrix();
	

	//AceGen
	
	v[215] = cos(c[2]);
	v[210] = sin(c[2]);
	v[250] = v[210] * (*zfac);
	v[93] = Power(dA[3], 2);
	v[91] = (dA[3] * dA[4]) / 2e0;
	v[86] = Power(dA[4], 2);
	v[98] = (dA[4] * dA[5]) / 2e0;
	v[96] = (dA[3] * dA[5]) / 2e0;
	v[87] = Power(dA[5], 2);
	v[246] = v[86] + v[87];
	v[112] = Power(dA[9], 2);
	v[110] = (dA[10] * dA[9]) / 2e0;
	v[105] = Power(dA[10], 2);
	v[117] = (dA[10] * dA[11]) / 2e0;
	v[115] = (dA[11] * dA[9]) / 2e0;
	v[106] = Power(dA[11], 2);
	v[247] = v[105] + v[106];
	v[193] = (*radA)*cos(c[1]);
	v[192] = (*radA)*sin(c[1]);
	v[162] = Power(dB[3], 2);
	v[160] = (dB[3] * dB[4]) / 2e0;
	v[155] = Power(dB[4], 2);
	v[167] = (dB[4] * dB[5]) / 2e0;
	v[165] = (dB[3] * dB[5]) / 2e0;
	v[156] = Power(dB[5], 2);
	v[248] = v[155] + v[156];
	v[209] = (*radB)*cos(c[3]);
	v[206] = (*radB)*sin(c[3]);
	v[187] = cpointB[0] + v[209];
	v[249] = v[187] * v[215];
	v[216] = -(v[187] * v[210] * (*xfac));
	v[208] = -(v[206] * v[215] * (*xfac));
	v[217] = -(v[249] * (*zfac));
	v[211] = v[206] * v[250];
	v[85] = 4e0 / (4e0 + v[246] + v[93]);
	v[88] = 1e0 - (v[246] * v[85]) / 2e0;
	v[89] = v[85] * (-dA[5] + v[91]);
	v[90] = v[85] * (dA[4] + v[96]);
	v[92] = v[85] * (dA[5] + v[91]);
	v[94] = 1e0 - (v[85] * (v[87] + v[93])) / 2e0;
	v[95] = v[85] * (-dA[3] + v[98]);
	v[97] = v[85] * (-dA[4] + v[96]);
	v[99] = v[85] * (dA[3] + v[98]);
	v[100] = 1e0 - (v[85] * (v[86] + v[93])) / 2e0;
	v[104] = 4e0 / (4e0 + v[112] + v[247]);
	v[107] = 1e0 - (v[104] * v[247]) / 2e0;
	v[108] = v[104] * (-dA[11] + v[110]);
	v[109] = v[104] * (dA[10] + v[115]);
	v[111] = v[104] * (dA[11] + v[110]);
	v[113] = 1e0 - (v[104] * (v[106] + v[112])) / 2e0;
	v[114] = v[104] * (-dA[9] + v[117]);
	v[116] = v[104] * (-dA[10] + v[115]);
	v[118] = v[104] * (dA[9] + v[117]);
	v[119] = 1e0 - (v[104] * (v[105] + v[112])) / 2e0;
	v[123] = QAAi[0][0] * v[88] + QAAi[1][0] * v[89] + QAAi[2][0] * v[90];
	v[124] = QAAi[0][1] * v[88] + QAAi[1][1] * v[89] + QAAi[2][1] * v[90];
	v[126] = QAAi[0][0] * v[92] + QAAi[1][0] * v[94] + QAAi[2][0] * v[95];
	v[127] = QAAi[0][1] * v[92] + QAAi[1][1] * v[94] + QAAi[2][1] * v[95];
	v[129] = QAAi[2][0] * v[100] + QAAi[0][0] * v[97] + QAAi[1][0] * v[99];
	v[130] = QAAi[2][1] * v[100] + QAAi[0][1] * v[97] + QAAi[1][1] * v[99];
	v[132] = QBAi[0][0] * v[107] + QBAi[1][0] * v[108] + QBAi[2][0] * v[109];
	v[133] = QBAi[0][1] * v[107] + QBAi[1][1] * v[108] + QBAi[2][1] * v[109];
	v[135] = QBAi[0][0] * v[111] + QBAi[1][0] * v[113] + QBAi[2][0] * v[114];
	v[136] = QBAi[0][1] * v[111] + QBAi[1][1] * v[113] + QBAi[2][1] * v[114];
	v[138] = QBAi[0][0] * v[116] + QBAi[1][0] * v[118] + QBAi[2][0] * v[119];
	v[139] = QBAi[0][1] * v[116] + QBAi[1][1] * v[118] + QBAi[2][1] * v[119];
	v[147] = (1e0 - c[0]) / 2e0;
	v[148] = (1e0 + c[0]) / 2e0;
	v[196] = v[147] * (-(v[129] * v[192]) + v[130] * v[193]) + v[148] * (-(v[138] * v[192]) + v[139] * v[193]);
	v[195] = v[147] * (-(v[126] * v[192]) + v[127] * v[193]) + v[148] * (-(v[135] * v[192]) + v[136] * v[193]);
	v[194] = v[147] * (-(v[123] * v[192]) + v[124] * v[193]) + v[148] * (-(v[132] * v[192]) + v[133] * v[193]);
	v[149] = cpointA[0] + v[193];
	v[150] = cpointA[1] + v[192];
	v[204] = dA[8] + v[138] * v[149] + v[139] * v[150] + xBAi[2];
	v[203] = dA[2] + v[129] * v[149] + v[130] * v[150] + xAAi[2];
	v[205] = (-v[203] + v[204]) / 2e0;
	v[201] = dA[7] + v[135] * v[149] + v[136] * v[150] + xBAi[1];
	v[200] = dA[1] + v[126] * v[149] + v[127] * v[150] + xAAi[1];
	v[202] = (-v[200] + v[201]) / 2e0;
	v[198] = dA[6] + v[132] * v[149] + v[133] * v[150] + xBAi[0];
	v[197] = dA[0] + v[123] * v[149] + v[124] * v[150] + xAAi[0];
	v[199] = (-v[197] + v[198]) / 2e0;
	v[154] = 4e0 / (4e0 + v[162] + v[248]);
	v[157] = 1e0 - (v[154] * v[248]) / 2e0;
	v[158] = v[154] * (-dB[5] + v[160]);
	v[159] = v[154] * (dB[4] + v[165]);
	v[161] = v[154] * (dB[5] + v[160]);
	v[163] = 1e0 - (v[154] * (v[156] + v[162])) / 2e0;
	v[164] = v[154] * (-dB[3] + v[167]);
	v[166] = v[154] * (-dB[4] + v[165]);
	v[168] = v[154] * (dB[3] + v[167]);
	v[169] = 1e0 - (v[154] * (v[155] + v[162])) / 2e0;
	v[173] = QABi[0][0] * v[157] + QABi[1][0] * v[158] + QABi[2][0] * v[159];
	v[174] = QABi[0][1] * v[157] + QABi[1][1] * v[158] + QABi[2][1] * v[159];
	v[175] = QABi[0][2] * v[157] + QABi[1][2] * v[158] + QABi[2][2] * v[159];
	v[218] = v[173] * v[216] + v[175] * v[217];
	v[212] = v[173] * v[208] + v[174] * v[209] + v[175] * v[211];
	v[176] = QABi[0][0] * v[161] + QABi[1][0] * v[163] + QABi[2][0] * v[164];
	v[177] = QABi[0][1] * v[161] + QABi[1][1] * v[163] + QABi[2][1] * v[164];
	v[178] = QABi[0][2] * v[161] + QABi[1][2] * v[163] + QABi[2][2] * v[164];
	v[219] = v[176] * v[216] + v[178] * v[217];
	v[213] = v[176] * v[208] + v[177] * v[209] + v[178] * v[211];
	v[179] = QABi[0][0] * v[166] + QABi[1][0] * v[168] + QABi[2][0] * v[169];
	v[180] = QABi[0][1] * v[166] + QABi[1][1] * v[168] + QABi[2][1] * v[169];
	v[181] = QABi[0][2] * v[166] + QABi[1][2] * v[168] + QABi[2][2] * v[169];
	v[220] = v[179] * v[216] + v[181] * v[217];
	v[214] = v[179] * v[208] + v[180] * v[209] + v[181] * v[211];
	v[185] = v[249] * (*xfac);
	v[186] = cpointB[1] + v[206];
	v[188] = -(v[187] * v[250]);
	v[242] = -dB[0] - v[173] * v[185] - v[174] * v[186] - v[175] * v[188] + v[147] * v[197] + v[148] * v[198] - xABi[0];
	v[243] = -dB[1] - v[176] * v[185] - v[177] * v[186] - v[178] * v[188] + v[147] * v[200] + v[148] * v[201] - xABi[1];
	v[244] = -dB[2] - v[179] * v[185] - v[180] * v[186] - v[181] * v[188] + v[147] * v[203] + v[148] * v[204] - xABi[2];
	v[223] = -(v[196] * v[202]) + v[195] * v[205];
	v[230] = v[214] * v[219] - v[213] * v[220];
	if ((*fixnormal)) {
		v[241] = 0.5e0*((-normalA[0] + normalB[0])*v[242] + (-normalA[1] + normalB[1])*v[243] + (-normalA[2]
			+ normalB[2])*v[244]);
	}
	else {
		v[235] = v[213] * v[218] - v[212] * v[219];
		v[233] = -(v[214] * v[218]) + v[212] * v[220];
		v[228] = -(v[195] * v[199]) + v[194] * v[202];
		v[226] = v[196] * v[199] - v[194] * v[205];
		v[241] = (-0.5e0*((*normalintA) ? -1 : 1)*(v[223] * v[242] + v[226] * v[243] + v[228] * v[244])) / sqrt(
			(v[223] * v[223]) + (v[226] * v[226]) + (v[228] * v[228])) + (0.5e0*((*normalintB) ? -1 : 1)*(v[230] * v[242]
				+ v[233] * v[243] + v[235] * v[244])) / sqrt((v[230] * v[230]) + (v[233] * v[233]) + (v[235] * v[235]));
	};
	(gap) = v[241];

	////////////////////////////////////////////////Desalocando memória//////////////////////////////////////////////////////
	for (int i = 0; i < 3; i++)
	{
		delete[] QAAi[i];
		delete[] QBAi[i];
		delete[] QABi[i];
	}
	delete[] QAAi;
	delete[] QBAi;
	delete[] QABi;

	return gap;
}

//Calcula o Gradiente do gap
void FlexibleArcExtrusion_1_RigidArcRevolution_1::GradientGap(Matrix& mc, Matrix& mGra, bool fixed_normals, Matrix& nA, Matrix& nB)
{
	double v[2000];		//variavel temporaria - AceGen

	bool* fixnormal = &fixed_normals;
	double* normalA = nA.getMatrix();
	double* normalB = nB.getMatrix();
	///////////////////////////////Ponteiros e variaveis para facilitar acesso//////////////////////////////////////////
	FlexibleArcExtrusion_1* surf1;		//Ponteiro para a superficie 1
	RigidArcRevolution_1* surf2;		//Ponteiro para a superficie 2
	double* radA;
	double* radB;
	double* xfac;
	double* zfac;
	double* cpointA;
	double* cpointB;
	double* dA;
	double* dB;
	double* duiA;
	double* dduiA;
	double* duiB;
	double* dduiB;
	bool* normalintA;
	bool* normalintB;
	double* xAAi;
	double* xBAi;
	double* xABi;
	double** QAAi;
	double** QBAi;
	double** QABi;
	surf1 = static_cast<FlexibleArcExtrusion_1*>(db.surfaces[surf1_ID - 1]);
	surf2 = static_cast<RigidArcRevolution_1*>(db.surfaces[surf2_ID - 1]);
	radA = surf1->radius;
	radB = surf2->radius;
	xfac = &surf2->x_fac;
	zfac = &surf2->z_fac;
	cpointA = surf1->c_point->getMatrix();
	cpointB = surf2->c_point->getMatrix();
	dA = surf1->d_A->getMatrix();
	dB = surf2->d_B->getMatrix();
	duiA = surf1->dui_A->getMatrix();
	dduiA = surf1->ddui_A->getMatrix();
	duiB = surf2->dui_B->getMatrix();
	dduiB = surf2->ddui_B->getMatrix();
	xAAi = surf1->x_AAi->getMatrix();
	xBAi = surf1->x_BAi->getMatrix();
	xABi = surf2->x_ABi->getMatrix();
	normalintA = &surf1->flag_normal_int;
	normalintB = &surf2->flag_normal_int;
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	QAAi = new double*[3];
	QBAi = new double*[3];
	QABi = new double*[3];
	for (int i = 0; i < 3; i++)
	{
		QAAi[i] = new double[3];
		QBAi[i] = new double[3];
		QABi[i] = new double[3];
	}
	//Salvando variaveis locais para montagem de superficies
	surf1->Q_AAi->MatrixToPtr(QAAi, 3);
	surf1->Q_BAi->MatrixToPtr(QBAi, 3);
	surf2->Q_ABi->MatrixToPtr(QABi, 3);
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	
	double *c = mc.getMatrix();
	double Gra[4];
	

	//AceGen
	int i246, i312, i313, b248;
	v[215] = cos(c[2]);
	v[317] = -(v[215] * (*xfac));
	v[210] = sin(c[2]);
	v[306] = v[210] * (*zfac);
	v[93] = Power(dA[3], 2);
	v[91] = (dA[3] * dA[4]) / 2e0;
	v[86] = Power(dA[4], 2);
	v[98] = (dA[4] * dA[5]) / 2e0;
	v[96] = (dA[3] * dA[5]) / 2e0;
	v[87] = Power(dA[5], 2);
	v[302] = v[86] + v[87];
	v[112] = Power(dA[9], 2);
	v[110] = (dA[10] * dA[9]) / 2e0;
	v[105] = Power(dA[10], 2);
	v[117] = (dA[10] * dA[11]) / 2e0;
	v[115] = (dA[11] * dA[9]) / 2e0;
	v[106] = Power(dA[11], 2);
	v[303] = v[105] + v[106];
	v[193] = (*radA)*cos(c[1]);
	v[192] = (*radA)*sin(c[1]);
	v[162] = Power(dB[3], 2);
	v[160] = (dB[3] * dB[4]) / 2e0;
	v[155] = Power(dB[4], 2);
	v[167] = (dB[4] * dB[5]) / 2e0;
	v[165] = (dB[3] * dB[5]) / 2e0;
	v[156] = Power(dB[5], 2);
	v[304] = v[155] + v[156];
	v[209] = (*radB)*cos(c[3]);
	v[206] = (*radB)*sin(c[3]);
	v[187] = cpointB[0] + v[209];
	v[305] = v[187] * v[215];
	v[216] = -(v[187] * v[210] * (*xfac));
	v[208] = v[206] * v[317];
	v[217] = -(v[305] * (*zfac));
	v[211] = v[206] * v[306];
	v[85] = 4e0 / (4e0 + v[302] + v[93]);
	v[88] = 1e0 - (v[302] * v[85]) / 2e0;
	v[89] = v[85] * (-dA[5] + v[91]);
	v[90] = v[85] * (dA[4] + v[96]);
	v[92] = v[85] * (dA[5] + v[91]);
	v[94] = 1e0 - (v[85] * (v[87] + v[93])) / 2e0;
	v[95] = v[85] * (-dA[3] + v[98]);
	v[97] = v[85] * (-dA[4] + v[96]);
	v[99] = v[85] * (dA[3] + v[98]);
	v[100] = 1e0 - (v[85] * (v[86] + v[93])) / 2e0;
	v[104] = 4e0 / (4e0 + v[112] + v[303]);
	v[107] = 1e0 - (v[104] * v[303]) / 2e0;
	v[108] = v[104] * (-dA[11] + v[110]);
	v[109] = v[104] * (dA[10] + v[115]);
	v[111] = v[104] * (dA[11] + v[110]);
	v[113] = 1e0 - (v[104] * (v[106] + v[112])) / 2e0;
	v[114] = v[104] * (-dA[9] + v[117]);
	v[116] = v[104] * (-dA[10] + v[115]);
	v[118] = v[104] * (dA[9] + v[117]);
	v[119] = 1e0 - (v[104] * (v[105] + v[112])) / 2e0;
	v[123] = QAAi[0][0] * v[88] + QAAi[1][0] * v[89] + QAAi[2][0] * v[90];
	v[124] = QAAi[0][1] * v[88] + QAAi[1][1] * v[89] + QAAi[2][1] * v[90];
	v[294] = -(v[123] * v[192]) + v[124] * v[193];
	v[126] = QAAi[0][0] * v[92] + QAAi[1][0] * v[94] + QAAi[2][0] * v[95];
	v[127] = QAAi[0][1] * v[92] + QAAi[1][1] * v[94] + QAAi[2][1] * v[95];
	v[295] = -(v[126] * v[192]) + v[127] * v[193];
	v[129] = QAAi[2][0] * v[100] + QAAi[0][0] * v[97] + QAAi[1][0] * v[99];
	v[130] = QAAi[2][1] * v[100] + QAAi[0][1] * v[97] + QAAi[1][1] * v[99];
	v[293] = -(v[129] * v[192]) + v[130] * v[193];
	v[132] = QBAi[0][0] * v[107] + QBAi[1][0] * v[108] + QBAi[2][0] * v[109];
	v[133] = QBAi[0][1] * v[107] + QBAi[1][1] * v[108] + QBAi[2][1] * v[109];
	v[297] = -(v[132] * v[192]) + v[133] * v[193];
	v[135] = QBAi[0][0] * v[111] + QBAi[1][0] * v[113] + QBAi[2][0] * v[114];
	v[136] = QBAi[0][1] * v[111] + QBAi[1][1] * v[113] + QBAi[2][1] * v[114];
	v[298] = -(v[135] * v[192]) + v[136] * v[193];
	v[138] = QBAi[0][0] * v[116] + QBAi[1][0] * v[118] + QBAi[2][0] * v[119];
	v[139] = QBAi[0][1] * v[116] + QBAi[1][1] * v[118] + QBAi[2][1] * v[119];
	v[296] = -(v[138] * v[192]) + v[139] * v[193];
	v[147] = (1e0 - c[0]) / 2e0;
	v[148] = (1e0 + c[0]) / 2e0;
	v[196] = v[147] * v[293] + v[148] * v[296];
	v[195] = v[147] * v[295] + v[148] * v[298];
	v[194] = v[147] * v[294] + v[148] * v[297];
	v[149] = cpointA[0] + v[193];
	v[150] = cpointA[1] + v[192];
	v[204] = dA[8] + v[138] * v[149] + v[139] * v[150] + xBAi[2];
	v[203] = dA[2] + v[129] * v[149] + v[130] * v[150] + xAAi[2];
	v[205] = (-v[203] + v[204]) / 2e0;
	v[201] = dA[7] + v[135] * v[149] + v[136] * v[150] + xBAi[1];
	v[200] = dA[1] + v[126] * v[149] + v[127] * v[150] + xAAi[1];
	v[202] = (-v[200] + v[201]) / 2e0;
	v[198] = dA[6] + v[132] * v[149] + v[133] * v[150] + xBAi[0];
	v[197] = dA[0] + v[123] * v[149] + v[124] * v[150] + xAAi[0];
	v[199] = (-v[197] + v[198]) / 2e0;
	v[228] = -(v[195] * v[199]) + v[194] * v[202];
	v[226] = v[196] * v[199] - v[194] * v[205];
	v[154] = 4e0 / (4e0 + v[162] + v[304]);
	v[157] = 1e0 - (v[154] * v[304]) / 2e0;
	v[158] = v[154] * (-dB[5] + v[160]);
	v[159] = v[154] * (dB[4] + v[165]);
	v[161] = v[154] * (dB[5] + v[160]);
	v[163] = 1e0 - (v[154] * (v[156] + v[162])) / 2e0;
	v[164] = v[154] * (-dB[3] + v[167]);
	v[166] = v[154] * (-dB[4] + v[165]);
	v[168] = v[154] * (dB[3] + v[167]);
	v[169] = 1e0 - (v[154] * (v[155] + v[162])) / 2e0;
	v[173] = QABi[0][0] * v[157] + QABi[1][0] * v[158] + QABi[2][0] * v[159];
	v[174] = QABi[0][1] * v[157] + QABi[1][1] * v[158] + QABi[2][1] * v[159];
	v[175] = QABi[0][2] * v[157] + QABi[1][2] * v[158] + QABi[2][2] * v[159];
	v[218] = v[173] * v[216] + v[175] * v[217];
	v[212] = v[173] * v[208] + v[174] * v[209] + v[175] * v[211];
	v[176] = QABi[0][0] * v[161] + QABi[1][0] * v[163] + QABi[2][0] * v[164];
	v[177] = QABi[0][1] * v[161] + QABi[1][1] * v[163] + QABi[2][1] * v[164];
	v[178] = QABi[0][2] * v[161] + QABi[1][2] * v[163] + QABi[2][2] * v[164];
	v[219] = v[176] * v[216] + v[178] * v[217];
	v[213] = v[176] * v[208] + v[177] * v[209] + v[178] * v[211];
	v[235] = v[213] * v[218] - v[212] * v[219];
	v[179] = QABi[0][0] * v[166] + QABi[1][0] * v[168] + QABi[2][0] * v[169];
	v[180] = QABi[0][1] * v[166] + QABi[1][1] * v[168] + QABi[2][1] * v[169];
	v[181] = QABi[0][2] * v[166] + QABi[1][2] * v[168] + QABi[2][2] * v[169];
	v[220] = v[179] * v[216] + v[181] * v[217];
	v[214] = v[179] * v[208] + v[180] * v[209] + v[181] * v[211];
	v[233] = -(v[214] * v[218]) + v[212] * v[220];
	v[185] = v[305] * (*xfac);
	v[186] = cpointB[1] + v[206];
	v[188] = -(v[187] * v[306]);
	v[242] = -dB[0] - v[173] * v[185] - v[174] * v[186] - v[175] * v[188] + v[147] * v[197] + v[148] * v[198] - xABi[0];
	v[243] = -dB[1] - v[176] * v[185] - v[177] * v[186] - v[178] * v[188] + v[147] * v[200] + v[148] * v[201] - xABi[1];
	v[244] = -dB[2] - v[179] * v[185] - v[180] * v[186] - v[181] * v[188] + v[147] * v[203] + v[148] * v[204] - xABi[2];
	v[223] = -(v[196] * v[202]) + v[195] * v[205];
	v[264] = (v[223] * v[223]) + (v[226] * v[226]) + (v[228] * v[228]);
	v[225] = 1e0 / sqrt(v[264]);
	v[230] = v[214] * v[219] - v[213] * v[220];
	v[262] = (v[230] * v[230]) + (v[233] * v[233]) + (v[235] * v[235]);
	v[232] = 1e0 / sqrt(v[262]);
	if ((*fixnormal)) {
		v[311] = -normalA[0] + normalB[0];
		v[310] = -normalA[1] + normalB[1];
		v[309] = -normalA[2] + normalB[2];
		v[241] = 0.5e0*(v[244] * v[309] + v[243] * v[310] + v[242] * v[311]);
	}
	else {
		i313 = ((*normalintA) ? -1 : 1);
		i312 = ((*normalintB) ? -1 : 1);
		v[308] = i313 * v[225];
		v[307] = i312 * v[232];
		v[314] = v[235] * v[307] - v[228] * v[308];
		v[315] = v[233] * v[307] - v[226] * v[308];
		v[316] = v[230] * v[307] - v[223] * v[308];
		v[241] = 0.5e0*(v[244] * v[314] + v[243] * v[315] + v[242] * v[316]);
	};
	b248 = (*fixnormal);
	if (b248) {
		v[249] = 0e0;
		v[250] = 0e0;
		v[251] = 0e0;
		v[252] = 0e0;
		v[253] = 0e0;
		v[254] = 0e0;
		v[255] = 0e0;
		v[256] = 0e0;
		v[257] = 0.5e0*v[309];
		v[258] = 0.5e0*v[310];
		v[259] = 0.5e0*v[311];
	}
	else {
		v[261] = -0.5e0*v[308];
		v[260] = 0.5e0*v[307];
		v[254] = 0.5e0*i312*(v[230] * v[242] + v[233] * v[243] + v[235] * v[244]);
		v[253] = v[242] * v[260];
		v[250] = -0.5e0*i313*(v[223] * v[242] + v[226] * v[243] + v[228] * v[244]);
		v[249] = v[242] * v[261];
		v[257] = 0.5e0*v[314];
		v[258] = 0.5e0*v[315];
		v[259] = 0.5e0*v[316];
		v[255] = v[243] * v[260];
		v[256] = 1e0*v[244] * v[260];
		v[251] = v[243] * v[261];
		v[252] = 1e0*v[244] * v[261];
	};
	v[265] = -((v[225] * v[250]) / v[264]);
	v[263] = -((v[232] * v[254]) / v[262]);
	v[253] = v[253] + v[230] * v[263];
	v[255] = v[255] + v[233] * v[263];
	v[256] = v[256] + v[235] * v[263];
	v[249] = v[249] + v[223] * v[265];
	v[251] = v[251] + v[226] * v[265];
	v[252] = v[252] + v[228] * v[265];
	v[268] = v[219] * v[253] - v[218] * v[255];
	v[269] = -(v[213] * v[253]) + v[212] * v[255];
	v[270] = -(v[220] * v[253]) + v[218] * v[256];
	v[271] = v[214] * v[253] - v[212] * v[256];
	v[272] = v[220] * v[255] - v[219] * v[256];
	v[273] = -(v[214] * v[255]) + v[213] * v[256];
	v[274] = v[181] * v[268] + v[178] * v[270] + v[175] * v[272];
	v[275] = v[179] * v[268] + v[176] * v[270] + v[173] * v[272];
	v[318] = (-(v[179] * v[257]) - v[176] * v[258] - v[173] * v[259])*(*xfac) - (v[181] * v[269] + v[178] * v[271]
		+ v[175] * v[273])*(*zfac);
	v[319] = -((v[179] * v[269] + v[176] * v[271] + v[173] * v[273])*(*xfac)) - (-(v[181] * v[257]) - v[178] * v[258]
		- v[175] * v[259])*(*zfac);
	v[291] = (-(v[195] * v[249]) + v[194] * v[251]) / 2e0;
	v[279] = -(v[202] * v[249]) + v[199] * v[251];
	v[285] = (-(v[196] * v[251]) + v[195] * v[252]) / 2e0;
	v[288] = (v[196] * v[249] - v[194] * v[252]) / 2e0;
	v[282] = -(v[205] * v[251]) + v[202] * v[252];
	v[283] = v[205] * v[249] - v[199] * v[252];
	v[284] = v[147] * v[259] + v[285];
	v[286] = v[148] * v[259] - v[285];
	v[287] = v[147] * v[258] + v[288];
	v[289] = v[148] * v[258] - v[288];
	v[290] = v[147] * v[257] + v[291];
	v[292] = v[148] * v[257] - v[291];
	v[368] = (-(v[203] * v[257]) + v[204] * v[257] - v[200] * v[258] + v[201] * v[258] - v[197] * v[259] + v[198] * v[259]
		- v[279] * v[293] - v[282] * v[294] - v[283] * v[295] + v[279] * v[296] + v[282] * v[297] + v[283] * v[298]) / 2e0;
	v[369] = -(v[192] * ((v[130] * v[147] + v[139] * v[148])*v[279] + (v[124] * v[147] + v[133] * v[148])*v[282] +
		(v[127] * v[147] + v[136] * v[148])*v[283] + v[123] * v[284] + v[132] * v[286] + v[126] * v[287] + v[135] * v[289]
		+ v[129] * v[290] + v[138] * v[292])) + v[193] * ((-(v[129] * v[147]) - v[138] * v[148])*v[279] + (-(v[123] * v[147])
			- v[132] * v[148])*v[282] + (-(v[126] * v[147]) - v[135] * v[148])*v[283] + v[124] * v[284] + v[133] * v[286]
			+ v[127] * v[287] + v[136] * v[289] + v[130] * v[290] + v[139] * v[292]);
	v[370] = -(v[210] * (v[187] * v[318] - v[206] * v[275] * (*xfac))) + v[215] * (v[187] * v[319] + v[206] * v[274] * (*zfac)
		);
	v[371] = v[209] * (-(v[180] * v[257]) - v[177] * v[258] - v[174] * v[259] + v[274] * v[306] + v[275] * v[317]) - v[206] *
		(v[180] * v[268] + v[177] * v[270] + v[174] * v[272] + v[215] * v[318] + v[210] * v[319]);
	
	for (i246 = 1; i246 <= 4; i246++) {
		Gra[i246 - 1] = v[367 + i246];
	};/* end for */

	for (int i = 0; i < 4; i++)
		mGra(i, 0) = Gra[i];

	////////////////////////////////////////////////Desalocando memória//////////////////////////////////////////////////////
	for (int i = 0; i < 3; i++)
	{
		delete[] QAAi[i];
		delete[] QBAi[i];
		delete[] QABi[i];
	}
	delete[] QAAi;
	delete[] QBAi;
	delete[] QABi;
}

//Calcula a Hessiana do gap
void FlexibleArcExtrusion_1_RigidArcRevolution_1::HessianGap(Matrix& mc, Matrix& mHes, bool fixed_normals, Matrix& nA, Matrix& nB)
{
	double v[2000];		//variavel temporaria - AceGen

	bool* fixnormal = &fixed_normals;
	double* normalA = nA.getMatrix();
	double* normalB = nB.getMatrix();
	///////////////////////////////Ponteiros e variaveis para facilitar acesso//////////////////////////////////////////
	FlexibleArcExtrusion_1* surf1;		//Ponteiro para a superficie 1
	RigidArcRevolution_1* surf2;		//Ponteiro para a superficie 2
	double* radA;
	double* radB;
	double* xfac;
	double* zfac;
	double* cpointA;
	double* cpointB;
	double* dA;
	double* dB;
	double* duiA;
	double* dduiA;
	double* duiB;
	double* dduiB;
	bool* normalintA;
	bool* normalintB;
	double* xAAi;
	double* xBAi;
	double* xABi;
	double** QAAi;
	double** QBAi;
	double** QABi;
	surf1 = static_cast<FlexibleArcExtrusion_1*>(db.surfaces[surf1_ID - 1]);
	surf2 = static_cast<RigidArcRevolution_1*>(db.surfaces[surf2_ID - 1]);
	radA = surf1->radius;
	radB = surf2->radius;
	xfac = &surf2->x_fac;
	zfac = &surf2->z_fac;
	cpointA = surf1->c_point->getMatrix();
	cpointB = surf2->c_point->getMatrix();
	dA = surf1->d_A->getMatrix();
	dB = surf2->d_B->getMatrix();
	duiA = surf1->dui_A->getMatrix();
	dduiA = surf1->ddui_A->getMatrix();
	duiB = surf2->dui_B->getMatrix();
	dduiB = surf2->ddui_B->getMatrix();
	xAAi = surf1->x_AAi->getMatrix();
	xBAi = surf1->x_BAi->getMatrix();
	xABi = surf2->x_ABi->getMatrix();
	normalintA = &surf1->flag_normal_int;
	normalintB = &surf2->flag_normal_int;
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	QAAi = new double*[3];
	QBAi = new double*[3];
	QABi = new double*[3];
	for (int i = 0; i < 3; i++)
	{
		QAAi[i] = new double[3];
		QBAi[i] = new double[3];
		QABi[i] = new double[3];
	}
	//Salvando variaveis locais para montagem de superficies
	surf1->Q_AAi->MatrixToPtr(QAAi, 3);
	surf1->Q_BAi->MatrixToPtr(QBAi, 3);
	surf2->Q_ABi->MatrixToPtr(QABi, 3);
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	
	double *c = mc.getMatrix();
	double Hes[4][4];
	
	//AceGen
	int i246, i302, i438, i441, i450, b248, b361, b462;
	v[215] = cos(c[2]);
	v[210] = sin(c[2]);
	v[93] = Power(dA[3], 2);
	v[91] = (dA[3] * dA[4]) / 2e0;
	v[86] = Power(dA[4], 2);
	v[98] = (dA[4] * dA[5]) / 2e0;
	v[96] = (dA[3] * dA[5]) / 2e0;
	v[87] = Power(dA[5], 2);
	v[423] = v[86] + v[87];
	v[112] = Power(dA[9], 2);
	v[110] = (dA[10] * dA[9]) / 2e0;
	v[105] = Power(dA[10], 2);
	v[117] = (dA[10] * dA[11]) / 2e0;
	v[115] = (dA[11] * dA[9]) / 2e0;
	v[106] = Power(dA[11], 2);
	v[424] = v[105] + v[106];
	v[193] = (*radA)*cos(c[1]);
	v[192] = (*radA)*sin(c[1]);
	v[162] = Power(dB[3], 2);
	v[160] = (dB[3] * dB[4]) / 2e0;
	v[155] = Power(dB[4], 2);
	v[167] = (dB[4] * dB[5]) / 2e0;
	v[165] = (dB[3] * dB[5]) / 2e0;
	v[156] = Power(dB[5], 2);
	v[425] = v[155] + v[156];
	v[209] = (*radB)*cos(c[3]);
	v[206] = (*radB)*sin(c[3]);
	v[315] = v[206] * v[215];
	v[635] = 0e0;
	v[636] = 0e0;
	v[637] = v[315];
	v[638] = v[209] * v[210];
	v[314] = v[206] * v[210];
	v[627] = 0e0;
	v[628] = 0e0;
	v[629] = v[314];
	v[630] = -(v[209] * v[215]);
	v[187] = cpointB[0] + v[209];
	v[427] = -(v[187] * v[210]);
	v[619] = 0e0;
	v[620] = 0e0;
	v[621] = -v[427];
	v[622] = v[315];
	v[426] = v[187] * v[215];
	v[611] = 0e0;
	v[612] = 0e0;
	v[613] = -v[426];
	v[614] = v[314];
	v[216] = v[427] * (*xfac);
	v[208] = -(v[315] * (*xfac));
	v[217] = -(v[426] * (*zfac));
	v[211] = v[314] * (*zfac);
	v[85] = 4e0 / (4e0 + v[423] + v[93]);
	v[88] = 1e0 - (v[423] * v[85]) / 2e0;
	v[89] = v[85] * (-dA[5] + v[91]);
	v[90] = v[85] * (dA[4] + v[96]);
	v[92] = v[85] * (dA[5] + v[91]);
	v[94] = 1e0 - (v[85] * (v[87] + v[93])) / 2e0;
	v[95] = v[85] * (-dA[3] + v[98]);
	v[97] = v[85] * (-dA[4] + v[96]);
	v[99] = v[85] * (dA[3] + v[98]);
	v[100] = 1e0 - (v[85] * (v[86] + v[93])) / 2e0;
	v[104] = 4e0 / (4e0 + v[112] + v[424]);
	v[107] = 1e0 - (v[104] * v[424]) / 2e0;
	v[108] = v[104] * (-dA[11] + v[110]);
	v[109] = v[104] * (dA[10] + v[115]);
	v[111] = v[104] * (dA[11] + v[110]);
	v[113] = 1e0 - (v[104] * (v[106] + v[112])) / 2e0;
	v[114] = v[104] * (-dA[9] + v[117]);
	v[116] = v[104] * (-dA[10] + v[115]);
	v[118] = v[104] * (dA[9] + v[117]);
	v[119] = 1e0 - (v[104] * (v[105] + v[112])) / 2e0;
	v[123] = QAAi[0][0] * v[88] + QAAi[1][0] * v[89] + QAAi[2][0] * v[90];
	v[124] = QAAi[0][1] * v[88] + QAAi[1][1] * v[89] + QAAi[2][1] * v[90];
	v[294] = -(v[123] * v[192]) + v[124] * v[193];
	v[126] = QAAi[0][0] * v[92] + QAAi[1][0] * v[94] + QAAi[2][0] * v[95];
	v[127] = QAAi[0][1] * v[92] + QAAi[1][1] * v[94] + QAAi[2][1] * v[95];
	v[295] = -(v[126] * v[192]) + v[127] * v[193];
	v[129] = QAAi[2][0] * v[100] + QAAi[0][0] * v[97] + QAAi[1][0] * v[99];
	v[130] = QAAi[2][1] * v[100] + QAAi[0][1] * v[97] + QAAi[1][1] * v[99];
	v[293] = -(v[129] * v[192]) + v[130] * v[193];
	v[132] = QBAi[0][0] * v[107] + QBAi[1][0] * v[108] + QBAi[2][0] * v[109];
	v[133] = QBAi[0][1] * v[107] + QBAi[1][1] * v[108] + QBAi[2][1] * v[109];
	v[297] = -(v[132] * v[192]) + v[133] * v[193];
	v[455] = -v[294] + v[297];
	v[135] = QBAi[0][0] * v[111] + QBAi[1][0] * v[113] + QBAi[2][0] * v[114];
	v[136] = QBAi[0][1] * v[111] + QBAi[1][1] * v[113] + QBAi[2][1] * v[114];
	v[298] = -(v[135] * v[192]) + v[136] * v[193];
	v[456] = -v[295] + v[298];
	v[138] = QBAi[0][0] * v[116] + QBAi[1][0] * v[118] + QBAi[2][0] * v[119];
	v[139] = QBAi[0][1] * v[116] + QBAi[1][1] * v[118] + QBAi[2][1] * v[119];
	v[296] = -(v[138] * v[192]) + v[139] * v[193];
	v[454] = -v[293] + v[296];
	v[147] = (1e0 - c[0]) / 2e0;
	v[148] = (1e0 + c[0]) / 2e0;
	v[336] = -(v[129] * v[147]) - v[138] * v[148];
	v[335] = v[130] * v[147] + v[139] * v[148];
	v[603] = -v[293] / 2e0 + v[296] / 2e0;
	v[604] = -(v[192] * v[335]) + v[193] * v[336];
	v[605] = 0e0;
	v[606] = 0e0;
	v[331] = -(v[123] * v[147]) - v[132] * v[148];
	v[330] = v[124] * v[147] + v[133] * v[148];
	v[599] = -v[294] / 2e0 + v[297] / 2e0;
	v[600] = -(v[192] * v[330]) + v[193] * v[331];
	v[601] = 0e0;
	v[602] = 0e0;
	v[329] = -(v[126] * v[147]) - v[135] * v[148];
	v[328] = v[127] * v[147] + v[136] * v[148];
	v[595] = -v[295] / 2e0 + v[298] / 2e0;
	v[596] = -(v[192] * v[328]) + v[193] * v[329];
	v[597] = 0e0;
	v[598] = 0e0;
	v[196] = v[147] * v[293] + v[148] * v[296];
	v[195] = v[147] * v[295] + v[148] * v[298];
	v[194] = v[147] * v[294] + v[148] * v[297];
	v[149] = cpointA[0] + v[193];
	v[150] = cpointA[1] + v[192];
	v[204] = dA[8] + v[138] * v[149] + v[139] * v[150] + xBAi[2];
	v[203] = dA[2] + v[129] * v[149] + v[130] * v[150] + xAAi[2];
	v[459] = -v[203] + v[204];
	v[205] = v[459] / 2e0;
	v[201] = dA[7] + v[135] * v[149] + v[136] * v[150] + xBAi[1];
	v[200] = dA[1] + v[126] * v[149] + v[127] * v[150] + xAAi[1];
	v[458] = -v[200] + v[201];
	v[202] = v[458] / 2e0;
	v[198] = dA[6] + v[132] * v[149] + v[133] * v[150] + xBAi[0];
	v[197] = dA[0] + v[123] * v[149] + v[124] * v[150] + xAAi[0];
	v[457] = -v[197] + v[198];
	v[199] = v[457] / 2e0;
	v[228] = -(v[195] * v[199]) + v[194] * v[202];
	v[226] = v[196] * v[199] - v[194] * v[205];
	v[154] = 4e0 / (4e0 + v[162] + v[425]);
	v[157] = 1e0 - (v[154] * v[425]) / 2e0;
	v[158] = v[154] * (-dB[5] + v[160]);
	v[159] = v[154] * (dB[4] + v[165]);
	v[161] = v[154] * (dB[5] + v[160]);
	v[163] = 1e0 - (v[154] * (v[156] + v[162])) / 2e0;
	v[164] = v[154] * (-dB[3] + v[167]);
	v[166] = v[154] * (-dB[4] + v[165]);
	v[168] = v[154] * (dB[3] + v[167]);
	v[169] = 1e0 - (v[154] * (v[155] + v[162])) / 2e0;
	v[173] = QABi[0][0] * v[157] + QABi[1][0] * v[158] + QABi[2][0] * v[159];
	v[174] = QABi[0][1] * v[157] + QABi[1][1] * v[158] + QABi[2][1] * v[159];
	v[326] = v[174] * v[209];
	v[175] = QABi[0][2] * v[157] + QABi[1][2] * v[158] + QABi[2][2] * v[159];
	v[218] = v[173] * v[216] + v[175] * v[217];
	v[212] = v[173] * v[208] + v[175] * v[211] + v[326];
	v[176] = QABi[0][0] * v[161] + QABi[1][0] * v[163] + QABi[2][0] * v[164];
	v[177] = QABi[0][1] * v[161] + QABi[1][1] * v[163] + QABi[2][1] * v[164];
	v[324] = v[177] * v[209];
	v[178] = QABi[0][2] * v[161] + QABi[1][2] * v[163] + QABi[2][2] * v[164];
	v[219] = v[176] * v[216] + v[178] * v[217];
	v[213] = v[176] * v[208] + v[178] * v[211] + v[324];
	v[235] = v[213] * v[218] - v[212] * v[219];
	v[179] = QABi[0][0] * v[166] + QABi[1][0] * v[168] + QABi[2][0] * v[169];
	v[180] = QABi[0][1] * v[166] + QABi[1][1] * v[168] + QABi[2][1] * v[169];
	v[322] = v[180] * v[209];
	v[181] = QABi[0][2] * v[166] + QABi[1][2] * v[168] + QABi[2][2] * v[169];
	v[220] = v[179] * v[216] + v[181] * v[217];
	v[214] = v[179] * v[208] + v[181] * v[211] + v[322];
	v[233] = -(v[214] * v[218]) + v[212] * v[220];
	v[185] = v[426] * (*xfac);
	v[186] = cpointB[1] + v[206];
	v[188] = v[427] * (*zfac);
	v[242] = -dB[0] - v[173] * v[185] - v[174] * v[186] - v[175] * v[188] + v[147] * v[197] + v[148] * v[198] - xABi[0];
	v[243] = -dB[1] - v[176] * v[185] - v[177] * v[186] - v[178] * v[188] + v[147] * v[200] + v[148] * v[201] - xABi[1];
	v[244] = -dB[2] - v[179] * v[185] - v[180] * v[186] - v[181] * v[188] + v[147] * v[203] + v[148] * v[204] - xABi[2];
	v[221] = ((*normalintA) ? -1 : 1);
	v[222] = ((*normalintB) ? -1 : 1);
	v[223] = -(v[196] * v[202]) + v[195] * v[205];
	v[264] = (v[223] * v[223]) + (v[226] * v[226]) + (v[228] * v[228]);
	v[477] = 1e0 / Power(v[264], 2);
	v[225] = 1e0 / sqrt(v[264]);
	v[230] = v[214] * v[219] - v[213] * v[220];
	v[262] = (v[230] * v[230]) + (v[233] * v[233]) + (v[235] * v[235]);
	v[476] = 1e0 / Power(v[262], 2);
	v[232] = 1e0 / sqrt(v[262]);
	if ((*fixnormal)) {
		v[432] = -normalA[0] + normalB[0];
		v[431] = -normalA[1] + normalB[1];
		v[430] = -normalA[2] + normalB[2];
		v[241] = 0.5e0*(v[244] * v[430] + v[243] * v[431] + v[242] * v[432]);
	}
	else {
		v[659] = v[205];
		v[660] = 0e0;
		v[661] = 0e0;
		v[662] = -v[322];
		v[663] = v[202];
		v[664] = 0e0;
		v[665] = 0e0;
		v[666] = -v[324];
		v[667] = v[199];
		v[668] = 0e0;
		v[669] = 0e0;
		v[670] = -v[326];
		v[429] = v[221] * v[225];
		v[428] = v[222] * v[232];
		v[433] = v[235] * v[428] - v[228] * v[429];
		v[434] = v[233] * v[428] - v[226] * v[429];
		v[435] = v[230] * v[428] - v[223] * v[429];
		v[241] = 0.5e0*(v[244] * v[433] + v[243] * v[434] + v[242] * v[435]);
	};
	b248 = (*fixnormal);
	if (b248) {
		v[249] = 0e0;
		v[250] = 0e0;
		v[251] = 0e0;
		v[252] = 0e0;
		v[253] = 0e0;
		v[254] = 0e0;
		v[255] = 0e0;
		v[256] = 0e0;
		v[257] = 0.5e0*v[430];
		v[258] = 0.5e0*v[431];
		v[259] = 0.5e0*v[432];
	}
	else {
		v[261] = -0.5e0*v[429];
		v[260] = 0.5e0*v[428];
		v[254] = 0.5e0*v[222] * (v[230] * v[242] + v[233] * v[243] + v[235] * v[244]);
		v[253] = v[242] * v[260];
		v[250] = -0.5e0*v[221] * (v[223] * v[242] + v[226] * v[243] + v[228] * v[244]);
		v[249] = v[242] * v[261];
		v[257] = 0.5e0*v[433];
		v[258] = 0.5e0*v[434];
		v[259] = 0.5e0*v[435];
		v[255] = v[243] * v[260];
		v[256] = 1e0*v[244] * v[260];
		v[251] = v[243] * v[261];
		v[252] = 1e0*v[244] * v[261];
	};
	v[449] = -(v[250] / v[264]);
	v[265] = v[225] * v[449];
	v[448] = -(v[254] / v[262]);
	v[263] = v[232] * v[448];
	v[253] = v[253] + v[230] * v[263];
	v[255] = v[255] + v[233] * v[263];
	v[256] = v[256] + v[235] * v[263];
	v[249] = v[249] + v[223] * v[265];
	v[251] = v[251] + v[226] * v[265];
	v[252] = v[252] + v[228] * v[265];
	v[266] = -(v[181] * v[257]) - v[178] * v[258] - v[175] * v[259];
	v[267] = -(v[179] * v[257]) - v[176] * v[258] - v[173] * v[259];
	v[268] = v[219] * v[253] - v[218] * v[255];
	v[269] = -(v[213] * v[253]) + v[212] * v[255];
	v[270] = -(v[220] * v[253]) + v[218] * v[256];
	v[271] = v[214] * v[253] - v[212] * v[256];
	v[272] = v[220] * v[255] - v[219] * v[256];
	v[273] = -(v[214] * v[255]) + v[213] * v[256];
	v[274] = v[181] * v[268] + v[178] * v[270] + v[175] * v[272];
	v[461] = v[274] * (*zfac);
	v[275] = v[179] * v[268] + v[176] * v[270] + v[173] * v[272];
	v[460] = -(v[275] * (*xfac));
	v[416] = -(v[180] * v[257]) - v[177] * v[258] - v[174] * v[259] + v[215] * v[460] + v[210] * v[461];
	v[276] = v[181] * v[269] + v[178] * v[271] + v[175] * v[273];
	v[436] = -(v[276] * (*zfac));
	v[413] = v[187] * v[436] + (v[187] * v[267] - v[206] * v[275])*(*xfac);
	v[412] = v[436] + v[267] * (*xfac);
	v[277] = v[179] * v[269] + v[176] * v[271] + v[173] * v[273];
	v[437] = -(v[277] * (*xfac));
	v[415] = v[180] * v[268] + v[177] * v[270] + v[174] * v[272] + (v[215] * v[267] - v[210] * v[277])*(*xfac) -
		(v[210] * v[266] + v[215] * v[276])*(*zfac);
	v[731] = 0e0;
	v[732] = 0e0;
	v[733] = 0e0;
	v[734] = -v[415];
	v[414] = v[437] - v[266] * (*zfac);
	v[482] = -(v[210] * v[412]) + v[215] * v[414];
	v[747] = 0e0;
	v[748] = 0e0;
	v[749] = -v[413];
	v[750] = -(v[206] * v[414]) + v[209] * v[461];
	v[411] = v[187] * v[437] + (-(v[187] * v[266]) + v[206] * v[274])*(*zfac);
	v[751] = 0e0;
	v[752] = 0e0;
	v[753] = v[411];
	v[754] = -(v[206] * v[412]) + v[209] * v[460];
	v[291] = (-(v[195] * v[249]) + v[194] * v[251]) / 2e0;
	v[279] = -(v[202] * v[249]) + v[199] * v[251];
	v[453] = v[192] * v[279];
	v[285] = (-(v[196] * v[251]) + v[195] * v[252]) / 2e0;
	v[288] = (v[196] * v[249] - v[194] * v[252]) / 2e0;
	v[282] = -(v[205] * v[251]) + v[202] * v[252];
	v[283] = v[205] * v[249] - v[199] * v[252];
	v[478] = v[193] * ((v[129] - v[138])*v[279] + (v[123] - v[132])*v[282] + (v[126] - v[135])*v[283]) + v[192] * (
		(v[124] - v[133])*v[282] + (v[127] - v[136])*v[283]) + v[130] * v[453] - v[139] * v[453] + v[257] * v[454]
		+ v[259] * v[455] + v[258] * v[456];
	v[284] = v[147] * v[259] + v[285];
	v[286] = v[148] * v[259] - v[285];
	v[287] = v[147] * v[258] + v[288];
	v[289] = v[148] * v[258] - v[288];
	v[290] = v[147] * v[257] + v[291];
	v[292] = v[148] * v[257] - v[291];
	v[410] = v[124] * v[284] + v[133] * v[286] + v[127] * v[287] + v[136] * v[289] + v[130] * v[290] + v[139] * v[292]
		+ v[283] * v[329] + v[282] * v[331] + v[279] * v[336];
	v[409] = v[123] * v[284] + v[132] * v[286] + v[126] * v[287] + v[135] * v[289] + v[129] * v[290] + v[138] * v[292]
		+ v[283] * v[328] + v[282] * v[330] + v[279] * v[335];
	v[559] = (-(v[203] * v[257]) + v[204] * v[257] - v[200] * v[258] + v[201] * v[258] - v[197] * v[259] + v[198] * v[259]
		- v[279] * v[293] - v[282] * v[294] - v[283] * v[295] + v[279] * v[296] + v[282] * v[297] + v[283] * v[298]) / 2e0;
	v[560] = -(v[192] * v[409]) + v[193] * v[410];
	v[561] = v[215] * v[411] - v[210] * v[413];
	v[562] = -(v[206] * v[415]) + v[209] * v[416];
	
	for (i246 = 1; i246 <= 4; i246++) {
		b462 = i246 == 4;
		i450 = i246 == 2;
		v[445] = (b462 ? -v[206] : 0e0);
		i441 = (i450 ? 1 : 0);
		v[442] = i441 / 2e0;
		i438 = (i246 == 1 ? 1 : 0);
		v[440] = -i438 / 2e0;
		v[439] = i438 / 2e0;
		v[407] = (i438*v[279]) / 2e0;
		v[404] = v[283] * v[439];
		v[401] = v[282] * v[439];
		v[397] = v[257] * v[440];
		v[393] = v[258] * v[440];
		v[388] = -(i438*v[259]) / 2e0;
		v[342] = v[594 + i246];
		v[341] = v[602 + i246];
		v[337] = v[598 + i246];
		v[447] = v[634 + i246] * (*zfac);
		v[446] = v[626 + i246] * (*xfac);
		v[311] = v[618 + i246];
		v[444] = v[311] * (*zfac);
		v[310] = v[610 + i246];
		v[443] = v[310] * (*xfac);
		v[306] = (i441*v[456]) / 2e0;
		v[307] = v[442] * v[455];
		v[333] = v[194] * v[306] - v[195] * v[307] + v[202] * v[337] - v[199] * v[342];
		v[308] = v[442] * v[454];
		v[344] = -(v[196] * v[306]) + v[195] * v[308] - v[202] * v[341] + v[205] * v[342];
		v[339] = v[196] * v[307] - v[194] * v[308] - v[205] * v[337] + v[199] * v[341];
		v[309] = v[173] * v[443] + v[175] * v[444];
		v[312] = v[176] * v[443] + v[178] * v[444];
		v[313] = v[179] * v[443] + v[181] * v[444];
		v[316] = v[174] * v[445] + v[173] * v[446] + v[175] * v[447];
		v[320] = v[177] * v[445] + v[176] * v[446] + v[178] * v[447];
		v[348] = v[213] * v[309] - v[212] * v[312] - v[219] * v[316] + v[218] * v[320];
		v[321] = v[180] * v[445] + v[179] * v[446] + v[181] * v[447];
		v[354] = v[214] * v[312] - v[213] * v[313] - v[220] * v[320] + v[219] * v[321];
		v[351] = -(v[214] * v[309]) + v[212] * v[313] + v[220] * v[316] - v[218] * v[321];
		v[332] = v[265] * v[333];
		v[338] = v[265] * v[339];
		v[343] = v[228] * v[333] + v[226] * v[339] + v[223] * v[344];
		v[345] = v[265] * v[344];
		v[347] = v[263] * v[348];
		v[350] = v[263] * v[351];
		v[353] = v[235] * v[348] + v[233] * v[351] + v[230] * v[354];
		v[355] = v[263] * v[354];
		v[358] = v[353] * v[448];
		v[360] = v[343] * v[449];
		b361 = (*fixnormal);
		if (b361) {
			v[362] = 0e0;
			v[363] = 0e0;
			v[364] = 0e0;
		}
		else {
			v[452] = v[311] * (*xfac);
			v[451] = -(v[310] * (*zfac));
			v[327] = (i450 ? v[194] : 0e0) + v[175] * v[451] + v[173] * v[452] + v[666 + i246];
			v[325] = (i450 ? v[195] : 0e0) + v[178] * v[451] + v[176] * v[452] + v[662 + i246];
			v[323] = (i450 ? v[196] : 0e0) + v[181] * v[451] + v[179] * v[452] + v[658 + i246];
			v[366] = (-0.5e0*v[353] * v[428]) / v[262];
			v[365] = (0.5e0*v[343] * v[429]) / v[264];
			v[345] = v[261] * v[327] + v[345] + v[242] * v[365];
			v[338] = v[261] * v[325] + v[338] + v[243] * v[365];
			v[332] = v[261] * v[323] + v[332] + 1e0*v[244] * v[365];
			v[355] = v[260] * v[327] + v[355] + v[242] * v[366];
			v[364] = 1e0*v[261] * v[333] + 1e0*v[260] * v[348] + 1e0*v[228] * v[365] + v[235] * v[366];
			v[363] = v[261] * v[339] + v[260] * v[351] + 1e0*v[226] * v[365] + 1e0*v[233] * v[366];
			v[362] = v[261] * v[344] + v[260] * v[354] + 1e0*v[223] * v[365] + 1e0*v[230] * v[366];
			v[350] = v[260] * v[325] + v[350] + 1e0*v[243] * v[366];
			v[347] = v[260] * v[323] + v[347] + 1e0*v[244] * v[366];
			v[358] = v[222] * (0.5e0*v[235] * v[323] + 0.5e0*v[233] * v[325] + 0.5e0*v[230] * v[327] + 0.5e0*v[244] * v[348]
				+ 0.5e0*v[243] * v[351] + 0.5e0*v[242] * v[354]) + v[358];
			v[360] = v[221] * (-0.5e0*v[228] * v[323] - 0.5e0*v[226] * v[325] - 0.5e0*v[223] * v[327] - 0.5e0*v[244] * v[333]
				- 0.5e0*v[243] * v[339] - 0.5e0*v[242] * v[344]) + v[360];
		};
		v[367] = -(v[232] * (-2e0*v[254] * v[353] + v[262] * v[358])*v[476]) / 2e0;
		v[355] = v[355] + 2e0*v[230] * v[367];
		v[350] = v[350] + 2e0*v[233] * v[367];
		v[347] = v[347] + 2e0*v[235] * v[367];
		v[368] = -(v[225] * (-2e0*v[250] * v[343] + v[264] * v[360])*v[477]) / 2e0;
		v[345] = v[345] + 2e0*v[223] * v[368];
		v[338] = v[338] + 2e0*v[226] * v[368];
		v[332] = v[332] + 2e0*v[228] * v[368];
		v[369] = -(v[175] * v[362]) - v[178] * v[363] - v[181] * v[364];
		v[370] = -(v[173] * v[362]) - v[176] * v[363] - v[179] * v[364];
		v[371] = -(v[255] * v[309]) + v[253] * v[312] - v[218] * v[350] + v[219] * v[355];
		v[372] = v[255] * v[316] - v[253] * v[320] + v[212] * v[350] - v[213] * v[355];
		v[373] = v[256] * v[309] - v[253] * v[313] + v[218] * v[347] - v[220] * v[355];
		v[374] = -(v[256] * v[316]) + v[253] * v[321] - v[212] * v[347] + v[214] * v[355];
		v[375] = -(v[256] * v[312]) + v[255] * v[313] - v[219] * v[347] + v[220] * v[350];
		v[376] = v[256] * v[320] - v[255] * v[321] + v[213] * v[347] - v[214] * v[350];
		v[377] = v[181] * v[371] + v[178] * v[373] + v[175] * v[375];
		v[378] = v[179] * v[371] + v[176] * v[373] + v[173] * v[375];
		v[379] = v[181] * v[372] + v[178] * v[374] + v[175] * v[376];
		v[380] = v[179] * v[372] + v[176] * v[374] + v[173] * v[376];
		v[398] = (v[251] * v[337] + v[194] * v[338] - v[249] * v[342] - v[195] * v[345]) / 2e0;
		v[382] = -(v[249] * v[306]) + v[251] * v[307] + v[199] * v[338] - v[202] * v[345];
		v[389] = (v[195] * v[332] - v[196] * v[338] - v[251] * v[341] + v[252] * v[342]) / 2e0;
		v[394] = (-(v[194] * v[332]) - v[252] * v[337] + v[249] * v[341] + v[196] * v[345]) / 2e0;
		v[385] = v[252] * v[306] - v[251] * v[308] + v[202] * v[332] - v[205] * v[338];
		v[386] = -(v[252] * v[307]) + v[249] * v[308] - v[199] * v[332] + v[205] * v[345];
		v[387] = v[147] * v[362] + v[388] + v[389];
		v[390] = v[148] * v[362] - v[388] - v[389];
		v[392] = v[147] * v[363] + v[393] + v[394];
		v[395] = v[148] * v[363] - v[393] - v[394];
		v[396] = v[147] * v[364] + v[397] + v[398];
		v[399] = v[148] * v[364] - v[397] - v[398];
		v[400] = v[148] * v[385] + v[401];
		v[402] = v[147] * v[385] - v[401];
		v[403] = v[148] * v[386] + v[404];
		v[405] = v[147] * v[386] - v[404];
		v[406] = v[148] * v[382] + v[407];
		v[408] = v[147] * v[382] - v[407];
		v[755] = ((i450 ? v[478] : 0e0) + v[382] * v[454] + v[385] * v[455] + v[386] * v[456] + v[362] * v[457] + v[363] * v[458]
			+ v[364] * v[459]) / 2e0;
		v[756] = v[193] * ((i450 ? -v[409] : 0e0) + v[124] * v[387] + v[133] * v[390] + v[127] * v[392] + v[136] * v[395]
			+ v[130] * v[396] + v[139] * v[399] - v[132] * v[400] - v[123] * v[402] - v[135] * v[403] - v[126] * v[405] - v[138] * v[406]
			- v[129] * v[408]) - v[192] * ((i450 ? v[410] : 0e0) + v[123] * v[387] + v[132] * v[390] + v[126] * v[392] + v[135] * v[395]
				+ v[129] * v[396] + v[138] * v[399] + v[133] * v[400] + v[124] * v[402] + v[136] * v[403] + v[127] * v[405] + v[139] * v[406]
				+ v[130] * v[408]);
		v[757] = v[215] * (v[746 + i246] - v[187] * v[380] * (*xfac) + (-(v[187] * v[369]) + v[206] * v[377])*(*zfac)) - v[210] *
			(v[750 + i246] + (v[187] * v[370] - v[206] * v[378])*(*xfac) - v[187] * v[379] * (*zfac));
		v[758] = v[209] * (-(v[174] * v[362]) - v[177] * v[363] - v[180] * v[364] + v[730 + i246] - v[215] * v[378] * (*xfac)
			+ v[210] * v[377] * (*zfac)) - v[206] * ((b462 ? v[416] : 0e0) + (i246 == 3 ? v[482] : 0e0) + v[180] * v[371] + v[177] * v[373]
				+ v[174] * v[375] + (v[215] * v[370] - v[210] * v[380])*(*xfac) - (v[210] * v[369] + v[215] * v[379])*(*zfac));
		
		for (i302 = i246; i302 <= 4; i302++) {
			v[418] = v[754 + i302];
			Hes[i246 - 1][i302 - 1] = v[418];
			if (i246 != i302) {
				Hes[i302 - 1][i246 - 1] = v[418];
			}
			else {
			};
		};/* end for */
	};/* end for */

	for (int i = 0; i < 4; i++)
	for (int j = 0; j < 4; j++)
		mHes(i, j) = Hes[i][j];

	////////////////////////////////////////////////Desalocando memória//////////////////////////////////////////////////////
	for (int i = 0; i < 3; i++)
	{
		delete[] QAAi[i];
		delete[] QBAi[i];
		delete[] QABi[i];
	}
	delete[] QAAi;
	delete[] QBAi;
	delete[] QABi;
}

//Chute inicial para coordenadas convectivas do par de superficies
void FlexibleArcExtrusion_1_RigidArcRevolution_1::InitialGuess(SSContactData* c_data)
{
	///////////////////////////////Ponteiros e variaveis para facilitar acesso//////////////////////////////////////////
	FlexibleArcExtrusion_1* surf1;		//Ponteiro para a superfÌcie 1
	RigidArcRevolution_1* surf2;		//Ponteiro para a superfÌcie 2
	surf1 = static_cast<FlexibleArcExtrusion_1*>(db.surfaces[surf1_ID - 1]);
	surf2 = static_cast<RigidArcRevolution_1*>(db.surfaces[surf2_ID - 1]);

	FlexibleArcExtrusion_1* surf_1 = static_cast<FlexibleArcExtrusion_1*>(db.surfaces[surf1_ID - 1]);
	int arc1id = surf_1->arc_ID;
	RigidArcRevolution_1* surf_2 = static_cast<RigidArcRevolution_1*>(db.surfaces[surf2_ID - 1]);
	int arc2id = surf_2->arc_ID;

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	Matrix xAAi = *surf1->x_AAi;
	Matrix xBAi = *surf1->x_BAi;
	Matrix xABi = *surf2->x_ABi;

	double radB = *surf2->radius;
	Matrix cpoint = *surf2->c_point;
	Matrix c(3);
	
	c(0, 0) = 0.0;
	c(1, 0) = 0.0;
	c(2, 0) = 0.0;

	Matrix Q = *surf2->Q_ABi;

	//Vetor c escrito no sistema global
	Matrix cglob;
	cglob = Q * c;

	//*** Inicio do chute inicial do zeta
	//Express„o obtida no Mathematica escrevendo um vetor parametrizado em theta para descrever a reta e projetando o ponto da superfÌcie 2 nessa reta. Impor ortogonalidade e encontrar o zeta que satisfaz a equaa„o. O zeta e obtido pela express„o abaixo
	//double zetaguess = ((xAAi(0, 0)*xAAi(0, 0)) + (xAAi(1, 0)*xAAi(1, 0)) + (xAAi(2, 0)*xAAi(2, 0)) - (2 * xAAi(0, 0)*xABi(0, 0)) - (2 * xAAi(1, 0)*xABi(1, 0)) - (2 * xAAi(2, 0)*xABi(2, 0)) + (2 * xABi(0, 0)*xBAi(0, 0)) - (xBAi(0, 0)*xBAi(0, 0)) + (2 * xABi(1, 0)*xBAi(1, 0)) - (xBAi(1, 0)*xBAi(1, 0)) + (2 * xABi(2, 0)*xBAi(2, 0)) - (xBAi(2, 0)*xBAi(2, 0))) / ((xAAi(0, 0)*xAAi(0, 0)) + (xAAi(1, 0)*xAAi(1, 0)) + (xAAi(2, 0)*xAAi(2, 0)) - (2 * xAAi(0, 0)*xBAi(0, 0)) + (xBAi(0, 0)*xBAi(0, 0)) - (2 * xAAi(1, 0)*xBAi(1, 0)) + (xBAi(1, 0)*xBAi(1, 0)) - (2 * xAAi(2, 0)*xBAi(2, 0)) + (xBAi(2, 0)*xBAi(2, 0)));

	//Express„o obtida no Mathematica escrevendo um vetor parametrizado em theta para descrever a reta e projetando o ponto da superfÌcie 2 + centro de curvatura do arco nessa reta. Impor ortogonalidade e encontrar o zeta que satisfaz a equaa„o. O zeta e obtido pela express„o abaixo
	double zetaguess = (-(2 * cglob(0, 0)*xAAi(0, 0)) + (xAAi(0, 0)*xAAi(0, 0)) - (2 * cglob(1, 0)*xAAi(1, 0)) + (xAAi(1, 0)*xAAi(1, 0)) - (2 * cglob(2, 0)*xAAi(2, 0)) + (xAAi(2, 0)*xAAi(2, 0)) - (2 * xAAi(0, 0)*xABi(0, 0)) - (2 * xAAi(1, 0)*xABi(1, 0)) - (2 * xAAi(2, 0)*xABi(2, 0)) + (2 * cglob(0, 0)*xBAi(0, 0)) + (2 * xABi(0, 0)*xBAi(0, 0)) - (xBAi(0, 0)*xBAi(0, 0)) + (2 * cglob(1, 0)*xBAi(1, 0)) + (2 * xABi(1, 0)*xBAi(1, 0)) - (xBAi(1, 0)*xBAi(1, 0)) - (2 * cglob(2, 0)*xBAi(2, 0)) + (2 * xABi(2, 0)*xBAi(2, 0)) - (xBAi(2, 0)*xBAi(2, 0))) / ((xAAi(0, 0)*xAAi(0, 0)) + (xAAi(1, 0)*xAAi(1, 0)) + (xAAi(2, 0)*xAAi(2, 0)) - (2 * xAAi(0, 0)*xBAi(0, 0)) + (xBAi(0, 0)*xBAi(0, 0)) - (2 * xAAi(1, 0)*xBAi(1, 0)) + (xBAi(1, 0)*xBAi(1, 0)) - (2 * xAAi(2, 0)*xBAi(2, 0)) + (xBAi(2, 0)*xBAi(2, 0)));


	//Vetor global que d· a posia„o do zetaguess
	Matrix zetaglob;
	zetaglob = 0.5*xAAi*(1 - zetaguess) + 0.5*xBAi*(1 + zetaguess);

	//*** Fim do chute inicial do zeta

	//Usando um vetor que vai do nÛ do RB ate o ponto dado por zetaglobal, porem, ao inves de pegar esse ponto sobre o eixo da viga, 
	//pegar esse ponto sobre um eixo paralelo que passa pelo centro de curvatura do arco extrudado
	//Centro de curvatura do arco extrudado
	Matrix cpoint_ext = *surf1->c_point;
	// Matriz rotaa„o/transformaa„o (nÛ A da viga, podia ser o nÛ B - pode causar problemas futuros fazer uma escolha aqui)
	Matrix Q_ext = *surf1->Q_AAi;
	//Escrever zetaglob no CS do arco extrudado
	Matrix zetaloc;
	zetaloc = transp(Q_ext)*zetaglob;
	//Transladar o ponto no plano x-y local
	Matrix zetaloct(3);
	zetaloct(0, 0) = zetaloc(0, 0) + cpoint_ext(0, 0);
	zetaloct(1, 0) = zetaloc(1, 0) + cpoint_ext(1, 0);
	zetaloct(2, 0) = zetaloc(2, 0);
	//Voltar para o CS global
	Matrix zetaglobt;
	zetaglobt = Q_ext * zetaloct;
	//Usar o zetaglobt ao inves de zetaglob para calcular vglob


	//*** Inicio do chite inicial do theta2 e do phi
	//Vetor global que vai do centro do corpo rÌgido ate o zetaguess
	Matrix vglob;
	vglob = zetaglobt - (xABi + cglob);

	//Vetor anterior escrito no sistema local do arco
	Matrix vloc;
	vloc = transp(Q)*vglob;

	Matrix versorvloc;
	versorvloc = (1 / norm(vloc))*vloc;

	if (versorvloc(1, 0) > 1.0)
		versorvloc(1, 0) = 1.0;
	if (versorvloc(1, 0) < -1.0)
		versorvloc(1, 0) = -1.0;

	double theta2guess;
	theta2guess = asin(versorvloc(1, 0));

	double cosphiguess;
	cosphiguess = (versorvloc(0, 0) / cos(theta2guess));

	if (cosphiguess > 1.0)
		cosphiguess = 1.0;
	if (cosphiguess < -1.0)
		cosphiguess = -1.0;

	double sinphiguess;
	sinphiguess = (-versorvloc(2, 0) / cos(theta2guess));

	if (sinphiguess > 1.0)
		sinphiguess = 1.0;
	if (sinphiguess < -1.0)
		sinphiguess = -1.0;

	double phiguess;

	if (sinphiguess >= 0 && cosphiguess >= 0)
		phiguess = asin(sinphiguess);
	if (sinphiguess >= 0 && cosphiguess < 0)
		phiguess = -asin(sinphiguess) + PI;
	if (sinphiguess < 0 && cosphiguess >= 0)
		phiguess = asin(sinphiguess);
	if (sinphiguess < 0 && cosphiguess < 0)
		phiguess = -asin(sinphiguess) - PI;

	//*** Fim do chute inicial do theta2 e do phi

	//*** Inicio do chute inicial do theta1

	Matrix chord = 0.5*(db.arcs[arc1id - 1]->f_point - db.arcs[arc1id - 1]->i_point);
	Matrix midarc = (db.arcs[arc1id - 1]->i_point + chord);

	double theta1guess;
	theta1guess = atan2(midarc(1, 0) - db.arcs[arc1id - 1]->c_point(1, 0), midarc(0, 0) - db.arcs[arc1id - 1]->c_point(0, 0));

	//*** Fim do chute inicial do theta1

	//Preenchendo as coordenadas convectivas:
	c_data->convective[0][0] = zetaguess;
	c_data->convective[0][1] = theta1guess;
	c_data->convective[0][2] = phiguess;
	c_data->convective[0][3] = theta2guess;

	//Inicial no q1 e final no q2
	if ((db.arcs[arc2id - 1]->theta_i >= 0 && db.arcs[arc2id - 1]->theta_i <= PI / 2) && (db.arcs[arc2id - 1]->theta_f >= PI / 2 && db.arcs[arc2id - 1]->theta_f <= PI))
	{
		if (theta2guess < 0)
			c_data->convective[0][3] = theta2guess + PI;
	}
	//Inicial e final no q2
	if ((db.arcs[arc2id - 1]->theta_i >= PI / 2 && db.arcs[arc2id - 1]->theta_i <= PI) && (db.arcs[arc2id - 1]->theta_f >= PI / 2 && db.arcs[arc2id - 1]->theta_f <= PI))
		c_data->convective[0][3] = theta2guess + PI;
	//Inicial no q2 e final no q3
	if ((db.arcs[arc2id - 1]->theta_i >= PI / 2 && db.arcs[arc2id - 1]->theta_i <= PI) && (db.arcs[arc2id - 1]->theta_f >= -PI && db.arcs[arc2id - 1]->theta_f <= -PI / 2))
		c_data->convective[0][3] = theta2guess + PI;
	//Inicial e final no q3
	if ((db.arcs[arc2id - 1]->theta_i >= -PI && db.arcs[arc2id - 1]->theta_i <= -PI / 2) && (db.arcs[arc2id - 1]->theta_f >= -PI && db.arcs[arc2id - 1]->theta_f <= -PI / 2))
		c_data->convective[0][3] = theta2guess + PI;
	//Inicial no q3 e final no q4
	if ((db.arcs[arc2id - 1]->theta_i >= -PI && db.arcs[arc2id - 1]->theta_i <= -PI / 2) && (db.arcs[arc2id - 1]->theta_f >= -PI / 2 && db.arcs[arc2id - 1]->theta_f <= 0))
	{
		if (theta2guess > 0)
			c_data->convective[0][3] = theta2guess + PI;
	}

	for (int ip = 1; ip < c_data->n_solutions; ip++)
	{
		//Preenchendo as coordenadas convectivas:
		c_data->convective[ip][0] = c_data->convective[0][0];
		c_data->convective[ip][1] = c_data->convective[0][1];
		c_data->convective[ip][2] = c_data->convective[0][2];
		c_data->convective[ip][3] = c_data->convective[0][3];
	}
}

//Calcula a função objetivo para um conjunto de coordenadas convectivas - Phase 1
double FlexibleArcExtrusion_1_RigidArcRevolution_1::ObjectivePhase1(Matrix& mc)
{
	double v[2000];		//variavel temporaria - AceGen

	///////////////////////////////Ponteiros e variaveis para facilitar acesso//////////////////////////////////////////
	FlexibleArcExtrusion_1* surf1;		//Ponteiro para a superficie 1
	RigidArcRevolution_1* surf2;		//Ponteiro para a superficie 2
	double* radA;
	double* radB;
	double* xfac;
	double* zfac;
	double* cpointA;
	double* cpointB;
	double* dA;
	double* dB;
	double* duiA;
	double* dduiA;
	double* duiB;
	double* dduiB;
	bool* normalintA;
	bool* normalintB;
	double* xAAi;
	double* xBAi;
	double* xABi;
	double** QAAi;
	double** QBAi;
	double** QABi;
	surf1 = static_cast<FlexibleArcExtrusion_1*>(db.surfaces[surf1_ID - 1]);
	surf2 = static_cast<RigidArcRevolution_1*>(db.surfaces[surf2_ID - 1]);
	radA = surf1->radius;
	radB = surf2->radius;
	xfac = &surf2->x_fac;
	zfac = &surf2->z_fac;
	cpointA = surf1->c_point->getMatrix();
	cpointB = surf2->c_point->getMatrix();
	dA = surf1->d_A->getMatrix();
	dB = surf2->d_B->getMatrix();
	duiA = surf1->dui_A->getMatrix();
	dduiA = surf1->ddui_A->getMatrix();
	duiB = surf2->dui_B->getMatrix();
	dduiB = surf2->ddui_B->getMatrix();
	xAAi = surf1->x_AAi->getMatrix();
	xBAi = surf1->x_BAi->getMatrix();
	xABi = surf2->x_ABi->getMatrix();
	normalintA = &surf1->flag_normal_int;
	normalintB = &surf2->flag_normal_int;
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	QAAi = new double*[3];
	QBAi = new double*[3];
	QABi = new double*[3];
	for (int i = 0; i < 3; i++)
	{
		QAAi[i] = new double[3];
		QBAi[i] = new double[3];
		QABi[i] = new double[3];
	}
	//Salvando variaveis locais para montagem de superficies
	surf1->Q_AAi->MatrixToPtr(QAAi, 3);
	surf1->Q_BAi->MatrixToPtr(QBAi, 3);
	surf2->Q_ABi->MatrixToPtr(QABi, 3);

	double *c = mc.getMatrix();
	double Ob;
	
	//Ace Gen
	v[93] = Power(dA[3], 2);
	v[91] = (dA[3] * dA[4]) / 2e0;
	v[86] = Power(dA[4], 2);
	v[98] = (dA[4] * dA[5]) / 2e0;
	v[96] = (dA[3] * dA[5]) / 2e0;
	v[87] = Power(dA[5], 2);
	v[246] = v[86] + v[87];
	v[112] = Power(dA[9], 2);
	v[110] = (dA[10] * dA[9]) / 2e0;
	v[105] = Power(dA[10], 2);
	v[117] = (dA[10] * dA[11]) / 2e0;
	v[115] = (dA[11] * dA[9]) / 2e0;
	v[106] = Power(dA[11], 2);
	v[247] = v[105] + v[106];
	v[162] = Power(dB[3], 2);
	v[160] = (dB[3] * dB[4]) / 2e0;
	v[155] = Power(dB[4], 2);
	v[167] = (dB[4] * dB[5]) / 2e0;
	v[165] = (dB[3] * dB[5]) / 2e0;
	v[156] = Power(dB[5], 2);
	v[248] = v[155] + v[156];
	v[187] = cpointB[0] + (*radB)*cos(c[3]);
	v[85] = 4e0 / (4e0 + v[246] + v[93]);
	v[88] = 1e0 - (v[246] * v[85]) / 2e0;
	v[89] = v[85] * (-dA[5] + v[91]);
	v[90] = v[85] * (dA[4] + v[96]);
	v[92] = v[85] * (dA[5] + v[91]);
	v[94] = 1e0 - (v[85] * (v[87] + v[93])) / 2e0;
	v[95] = v[85] * (-dA[3] + v[98]);
	v[97] = v[85] * (-dA[4] + v[96]);
	v[99] = v[85] * (dA[3] + v[98]);
	v[100] = 1e0 - (v[85] * (v[86] + v[93])) / 2e0;
	v[104] = 4e0 / (4e0 + v[112] + v[247]);
	v[107] = 1e0 - (v[104] * v[247]) / 2e0;
	v[108] = v[104] * (-dA[11] + v[110]);
	v[109] = v[104] * (dA[10] + v[115]);
	v[111] = v[104] * (dA[11] + v[110]);
	v[113] = 1e0 - (v[104] * (v[106] + v[112])) / 2e0;
	v[114] = v[104] * (-dA[9] + v[117]);
	v[116] = v[104] * (-dA[10] + v[115]);
	v[118] = v[104] * (dA[9] + v[117]);
	v[119] = 1e0 - (v[104] * (v[105] + v[112])) / 2e0;
	v[147] = (1e0 - c[0]) / 2e0;
	v[148] = (1e0 + c[0]) / 2e0;
	v[149] = cpointA[0] + (*radA)*cos(c[1]);
	v[150] = cpointA[1] + (*radA)*sin(c[1]);
	v[154] = 4e0 / (4e0 + v[162] + v[248]);
	v[157] = 1e0 - (v[154] * v[248]) / 2e0;
	v[158] = v[154] * (-dB[5] + v[160]);
	v[159] = v[154] * (dB[4] + v[165]);
	v[161] = v[154] * (dB[5] + v[160]);
	v[163] = 1e0 - (v[154] * (v[156] + v[162])) / 2e0;
	v[164] = v[154] * (-dB[3] + v[167]);
	v[166] = v[154] * (-dB[4] + v[165]);
	v[168] = v[154] * (dB[3] + v[167]);
	v[169] = 1e0 - (v[154] * (v[155] + v[162])) / 2e0;
	v[185] = v[187] * (*xfac)*cos(c[2]);
	v[186] = cpointB[1] + (*radB)*sin(c[3]);
	v[188] = -(v[187] * (*zfac)*sin(c[2]));
	(Ob) = (Power(-dB[0] - (QABi[0][0] * v[157] + QABi[1][0] * v[158] + QABi[2][0] * v[159])*v[185] -
		(QABi[0][1] * v[157] + QABi[1][1] * v[158] + QABi[2][1] * v[159])*v[186] - (QABi[0][2] * v[157] + QABi[1][2] * v[158]
			+ QABi[2][2] * v[159])*v[188] + v[147] * (dA[0] + v[149] * (QAAi[0][0] * v[88] + QAAi[1][0] * v[89] + QAAi[2][0] * v[90]
				) + v[150] * (QAAi[0][1] * v[88] + QAAi[1][1] * v[89] + QAAi[2][1] * v[90]) + xAAi[0]) - xABi[0] + v[148] * (dA[6] +
				(QBAi[0][0] * v[107] + QBAi[1][0] * v[108] + QBAi[2][0] * v[109])*v[149] + (QBAi[0][1] * v[107] + QBAi[1][1] * v[108]
					+ QBAi[2][1] * v[109])*v[150] + xBAi[0]), 2) + Power(-dB[1] - (QABi[0][0] * v[161] + QABi[1][0] * v[163]
						+ QABi[2][0] * v[164])*v[185] - (QABi[0][1] * v[161] + QABi[1][1] * v[163] + QABi[2][1] * v[164])*v[186] -
						(QABi[0][2] * v[161] + QABi[1][2] * v[163] + QABi[2][2] * v[164])*v[188] + v[147] * (dA[1] + v[149] *
						(QAAi[0][0] * v[92] + QAAi[1][0] * v[94] + QAAi[2][0] * v[95]) + v[150] * (QAAi[0][1] * v[92] + QAAi[1][1] * v[94]
							+ QAAi[2][1] * v[95]) + xAAi[1]) - xABi[1] + v[148] * (dA[7] + (QBAi[0][0] * v[111] + QBAi[1][0] * v[113]
								+ QBAi[2][0] * v[114])*v[149] + (QBAi[0][1] * v[111] + QBAi[1][1] * v[113] + QBAi[2][1] * v[114])*v[150] + xBAi[1])
						, 2) + Power(-dB[2] - (QABi[0][0] * v[166] + QABi[1][0] * v[168] + QABi[2][0] * v[169])*v[185] - (QABi[0][1] * v[166]
							+ QABi[1][1] * v[168] + QABi[2][1] * v[169])*v[186] - (QABi[0][2] * v[166] + QABi[1][2] * v[168] + QABi[2][2] * v[169]
								)*v[188] + v[147] * (dA[2] + v[149] * (QAAi[2][0] * v[100] + QAAi[0][0] * v[97] + QAAi[1][0] * v[99]) + v[150] *
								(QAAi[2][1] * v[100] + QAAi[0][1] * v[97] + QAAi[1][1] * v[99]) + xAAi[2]) - xABi[2] + v[148] * (dA[8] +
									(QBAi[0][0] * v[116] + QBAi[1][0] * v[118] + QBAi[2][0] * v[119])*v[149] + (QBAi[0][1] * v[116] + QBAi[1][1] * v[118]
										+ QBAi[2][1] * v[119])*v[150] + xBAi[2]), 2)) / 2e0;
	//

	////////////////////////////////////////////////Desalocando memória//////////////////////////////////////////////////////
	for (int i = 0; i < 3; i++)
	{
		delete[] QAAi[i];
		delete[] QBAi[i];
		delete[] QABi[i];
	}
	delete[] QAAi;
	delete[] QBAi;
	delete[] QABi;

	return Ob;
}

//Calcula o Gradiente da função objetivo - Phase 1
void FlexibleArcExtrusion_1_RigidArcRevolution_1::GradientPhase1(Matrix& mc, Matrix& mGra)
{
	double v[2000];		//variavel temporaria - AceGen

	///////////////////////////////Ponteiros e variaveis para facilitar acesso//////////////////////////////////////////
	FlexibleArcExtrusion_1* surf1;		//Ponteiro para a superficie 1
	RigidArcRevolution_1* surf2;		//Ponteiro para a superficie 2
	double* radA;
	double* radB;
	double* xfac;
	double* zfac;
	double* cpointA;
	double* cpointB;
	double* dA;
	double* dB;
	double* duiA;
	double* dduiA;
	double* duiB;
	double* dduiB;
	bool* normalintA;
	bool* normalintB;
	double* xAAi;
	double* xBAi;
	double* xABi;
	double** QAAi;
	double** QBAi;
	double** QABi;
	surf1 = static_cast<FlexibleArcExtrusion_1*>(db.surfaces[surf1_ID - 1]);
	surf2 = static_cast<RigidArcRevolution_1*>(db.surfaces[surf2_ID - 1]);
	radA = surf1->radius;
	radB = surf2->radius;
	xfac = &surf2->x_fac;
	zfac = &surf2->z_fac;
	cpointA = surf1->c_point->getMatrix();
	cpointB = surf2->c_point->getMatrix();
	dA = surf1->d_A->getMatrix();
	dB = surf2->d_B->getMatrix();
	duiA = surf1->dui_A->getMatrix();
	dduiA = surf1->ddui_A->getMatrix();
	duiB = surf2->dui_B->getMatrix();
	dduiB = surf2->ddui_B->getMatrix();
	xAAi = surf1->x_AAi->getMatrix();
	xBAi = surf1->x_BAi->getMatrix();
	xABi = surf2->x_ABi->getMatrix();
	normalintA = &surf1->flag_normal_int;
	normalintB = &surf2->flag_normal_int;
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	QAAi = new double*[3];
	QBAi = new double*[3];
	QABi = new double*[3];
	for (int i = 0; i < 3; i++)
	{
		QAAi[i] = new double[3];
		QBAi[i] = new double[3];
		QABi[i] = new double[3];
	}
	//Salvando variaveis locais para montagem de superficies
	surf1->Q_AAi->MatrixToPtr(QAAi, 3);
	surf1->Q_BAi->MatrixToPtr(QBAi, 3);
	surf2->Q_ABi->MatrixToPtr(QABi, 3);
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	double *c = mc.getMatrix();
	double Gra[4];

	//AceGen
	int i246;
	v[215] = cos(c[2]);
	v[258] = v[215] * (*xfac);
	v[210] = sin(c[2]);
	v[257] = -(v[210] * (*zfac));
	v[93] = Power(dA[3], 2);
	v[91] = (dA[3] * dA[4]) / 2e0;
	v[86] = Power(dA[4], 2);
	v[98] = (dA[4] * dA[5]) / 2e0;
	v[96] = (dA[3] * dA[5]) / 2e0;
	v[87] = Power(dA[5], 2);
	v[254] = v[86] + v[87];
	v[112] = Power(dA[9], 2);
	v[110] = (dA[10] * dA[9]) / 2e0;
	v[105] = Power(dA[10], 2);
	v[117] = (dA[10] * dA[11]) / 2e0;
	v[115] = (dA[11] * dA[9]) / 2e0;
	v[106] = Power(dA[11], 2);
	v[255] = v[105] + v[106];
	v[193] = (*radA)*cos(c[1]);
	v[192] = (*radA)*sin(c[1]);
	v[162] = Power(dB[3], 2);
	v[160] = (dB[3] * dB[4]) / 2e0;
	v[155] = Power(dB[4], 2);
	v[167] = (dB[4] * dB[5]) / 2e0;
	v[165] = (dB[3] * dB[5]) / 2e0;
	v[156] = Power(dB[5], 2);
	v[256] = v[155] + v[156];
	v[209] = (*radB)*cos(c[3]);
	v[206] = (*radB)*sin(c[3]);
	v[187] = cpointB[0] + v[209];
	v[85] = 4e0 / (4e0 + v[254] + v[93]);
	v[88] = 1e0 - (v[254] * v[85]) / 2e0;
	v[89] = v[85] * (-dA[5] + v[91]);
	v[90] = v[85] * (dA[4] + v[96]);
	v[92] = v[85] * (dA[5] + v[91]);
	v[94] = 1e0 - (v[85] * (v[87] + v[93])) / 2e0;
	v[95] = v[85] * (-dA[3] + v[98]);
	v[97] = v[85] * (-dA[4] + v[96]);
	v[99] = v[85] * (dA[3] + v[98]);
	v[100] = 1e0 - (v[85] * (v[86] + v[93])) / 2e0;
	v[104] = 4e0 / (4e0 + v[112] + v[255]);
	v[107] = 1e0 - (v[104] * v[255]) / 2e0;
	v[108] = v[104] * (-dA[11] + v[110]);
	v[109] = v[104] * (dA[10] + v[115]);
	v[111] = v[104] * (dA[11] + v[110]);
	v[113] = 1e0 - (v[104] * (v[106] + v[112])) / 2e0;
	v[114] = v[104] * (-dA[9] + v[117]);
	v[116] = v[104] * (-dA[10] + v[115]);
	v[118] = v[104] * (dA[9] + v[117]);
	v[119] = 1e0 - (v[104] * (v[105] + v[112])) / 2e0;
	v[123] = QAAi[0][0] * v[88] + QAAi[1][0] * v[89] + QAAi[2][0] * v[90];
	v[124] = QAAi[0][1] * v[88] + QAAi[1][1] * v[89] + QAAi[2][1] * v[90];
	v[126] = QAAi[0][0] * v[92] + QAAi[1][0] * v[94] + QAAi[2][0] * v[95];
	v[127] = QAAi[0][1] * v[92] + QAAi[1][1] * v[94] + QAAi[2][1] * v[95];
	v[129] = QAAi[2][0] * v[100] + QAAi[0][0] * v[97] + QAAi[1][0] * v[99];
	v[130] = QAAi[2][1] * v[100] + QAAi[0][1] * v[97] + QAAi[1][1] * v[99];
	v[132] = QBAi[0][0] * v[107] + QBAi[1][0] * v[108] + QBAi[2][0] * v[109];
	v[133] = QBAi[0][1] * v[107] + QBAi[1][1] * v[108] + QBAi[2][1] * v[109];
	v[135] = QBAi[0][0] * v[111] + QBAi[1][0] * v[113] + QBAi[2][0] * v[114];
	v[136] = QBAi[0][1] * v[111] + QBAi[1][1] * v[113] + QBAi[2][1] * v[114];
	v[138] = QBAi[0][0] * v[116] + QBAi[1][0] * v[118] + QBAi[2][0] * v[119];
	v[139] = QBAi[0][1] * v[116] + QBAi[1][1] * v[118] + QBAi[2][1] * v[119];
	v[147] = (1e0 - c[0]) / 2e0;
	v[148] = (1e0 + c[0]) / 2e0;
	v[149] = cpointA[0] + v[193];
	v[150] = cpointA[1] + v[192];
	v[204] = dA[8] + v[138] * v[149] + v[139] * v[150] + xBAi[2];
	v[203] = dA[2] + v[129] * v[149] + v[130] * v[150] + xAAi[2];
	v[201] = dA[7] + v[135] * v[149] + v[136] * v[150] + xBAi[1];
	v[200] = dA[1] + v[126] * v[149] + v[127] * v[150] + xAAi[1];
	v[198] = dA[6] + v[132] * v[149] + v[133] * v[150] + xBAi[0];
	v[197] = dA[0] + v[123] * v[149] + v[124] * v[150] + xAAi[0];
	v[154] = 4e0 / (4e0 + v[162] + v[256]);
	v[157] = 1e0 - (v[154] * v[256]) / 2e0;
	v[158] = v[154] * (-dB[5] + v[160]);
	v[159] = v[154] * (dB[4] + v[165]);
	v[161] = v[154] * (dB[5] + v[160]);
	v[163] = 1e0 - (v[154] * (v[156] + v[162])) / 2e0;
	v[164] = v[154] * (-dB[3] + v[167]);
	v[166] = v[154] * (-dB[4] + v[165]);
	v[168] = v[154] * (dB[3] + v[167]);
	v[169] = 1e0 - (v[154] * (v[155] + v[162])) / 2e0;
	v[173] = QABi[0][0] * v[157] + QABi[1][0] * v[158] + QABi[2][0] * v[159];
	v[174] = QABi[0][1] * v[157] + QABi[1][1] * v[158] + QABi[2][1] * v[159];
	v[175] = QABi[0][2] * v[157] + QABi[1][2] * v[158] + QABi[2][2] * v[159];
	v[176] = QABi[0][0] * v[161] + QABi[1][0] * v[163] + QABi[2][0] * v[164];
	v[177] = QABi[0][1] * v[161] + QABi[1][1] * v[163] + QABi[2][1] * v[164];
	v[178] = QABi[0][2] * v[161] + QABi[1][2] * v[163] + QABi[2][2] * v[164];
	v[179] = QABi[0][0] * v[166] + QABi[1][0] * v[168] + QABi[2][0] * v[169];
	v[180] = QABi[0][1] * v[166] + QABi[1][1] * v[168] + QABi[2][1] * v[169];
	v[181] = QABi[0][2] * v[166] + QABi[1][2] * v[168] + QABi[2][2] * v[169];
	v[185] = v[187] * v[258];
	v[186] = cpointB[1] + v[206];
	v[188] = v[187] * v[257];
	v[242] = -dB[0] - v[173] * v[185] - v[174] * v[186] - v[175] * v[188] + v[147] * v[197] + v[148] * v[198] - xABi[0];
	v[243] = -dB[1] - v[176] * v[185] - v[177] * v[186] - v[178] * v[188] + v[147] * v[200] + v[148] * v[201] - xABi[1];
	v[248] = -dB[2] - v[179] * v[185] - v[180] * v[186] - v[181] * v[188] + v[147] * v[203] + v[148] * v[204] - xABi[2];
	v[249] = -(v[175] * v[242]) - v[178] * v[243] - v[181] * v[248];
	v[250] = -(v[173] * v[242]) - v[176] * v[243] - v[179] * v[248];
	v[279] = (-(v[197] * v[242]) + v[198] * v[242] - v[200] * v[243] + v[201] * v[243] - v[203] * v[248] + v[204] * v[248])
		/ 2e0;
	v[280] = -(v[192] * ((v[123] * v[147] + v[132] * v[148])*v[242] + (v[126] * v[147] + v[135] * v[148])*v[243] +
		(v[129] * v[147] + v[138] * v[148])*v[248])) + v[193] * ((v[124] * v[147] + v[133] * v[148])*v[242] + (v[127] * v[147]
			+ v[136] * v[148])*v[243] + (v[130] * v[147] + v[139] * v[148])*v[248]);
	v[281] = v[187] * (-(v[210] * v[250] * (*xfac)) - v[215] * v[249] * (*zfac));
	v[282] = v[209] * (-(v[174] * v[242]) - v[177] * v[243] - v[180] * v[248]) - v[206] * (v[249] * v[257] + v[250] * v[258]);
	
	for (i246 = 1; i246 <= 4; i246++) {
		Gra[i246 - 1] = v[278 + i246];
	};/* end for */

	for (int i = 0; i < 4; i++)
		mGra(i, 0) = Gra[i];

	////////////////////////////////////////////////Desalocando memória//////////////////////////////////////////////////////
	for (int i = 0; i < 3; i++)
	{
		delete[] QAAi[i];
		delete[] QBAi[i];
		delete[] QABi[i];
	}
	delete[] QAAi;
	delete[] QBAi;
	delete[] QABi;

}

//Calcula a Hessiana da função objetivo - Phase 1
void FlexibleArcExtrusion_1_RigidArcRevolution_1::HessianPhase1(Matrix& mc, Matrix& mHes)
{
	double v[2000];		//variavel temporaria - AceGen

	///////////////////////////////Ponteiros e variaveis para facilitar acesso//////////////////////////////////////////
	FlexibleArcExtrusion_1* surf1;		//Ponteiro para a superficie 1
	RigidArcRevolution_1* surf2;		//Ponteiro para a superficie 2
	double* radA;
	double* radB;
	double* xfac;
	double* zfac;
	double* cpointA;
	double* cpointB;
	double* dA;
	double* dB;
	double* duiA;
	double* dduiA;
	double* duiB;
	double* dduiB;
	bool* normalintA;
	bool* normalintB;
	double* xAAi;
	double* xBAi;
	double* xABi;
	double** QAAi;
	double** QBAi;
	double** QABi;
	surf1 = static_cast<FlexibleArcExtrusion_1*>(db.surfaces[surf1_ID - 1]);
	surf2 = static_cast<RigidArcRevolution_1*>(db.surfaces[surf2_ID - 1]);
	radA = surf1->radius;
	radB = surf2->radius;
	xfac = &surf2->x_fac;
	zfac = &surf2->z_fac;
	cpointA = surf1->c_point->getMatrix();
	cpointB = surf2->c_point->getMatrix();
	dA = surf1->d_A->getMatrix();
	dB = surf2->d_B->getMatrix();
	duiA = surf1->dui_A->getMatrix();
	dduiA = surf1->ddui_A->getMatrix();
	duiB = surf2->dui_B->getMatrix();
	dduiB = surf2->ddui_B->getMatrix();
	xAAi = surf1->x_AAi->getMatrix();
	xBAi = surf1->x_BAi->getMatrix();
	xABi = surf2->x_ABi->getMatrix();
	normalintA = &surf1->flag_normal_int;
	normalintB = &surf2->flag_normal_int;
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	QAAi = new double*[3];
	QBAi = new double*[3];
	QABi = new double*[3];
	for (int i = 0; i < 3; i++)
	{
		QAAi[i] = new double[3];
		QBAi[i] = new double[3];
		QABi[i] = new double[3];
	}
	//Salvando variaveis locais para montagem de superficies
	surf1->Q_AAi->MatrixToPtr(QAAi, 3);
	surf1->Q_BAi->MatrixToPtr(QBAi, 3);
	surf2->Q_ABi->MatrixToPtr(QABi, 3);
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	double *c = mc.getMatrix();
	double Hes[4][4];

	//AceGen
	int i246, i254, i295, b297;
	v[215] = cos(c[2]);
	v[294] = v[215] * (*xfac);
	v[210] = sin(c[2]);
	v[303] = v[210] * (*xfac);
	v[293] = v[210] * (*zfac);
	v[93] = Power(dA[3], 2);
	v[91] = (dA[3] * dA[4]) / 2e0;
	v[86] = Power(dA[4], 2);
	v[98] = (dA[4] * dA[5]) / 2e0;
	v[96] = (dA[3] * dA[5]) / 2e0;
	v[87] = Power(dA[5], 2);
	v[289] = v[86] + v[87];
	v[112] = Power(dA[9], 2);
	v[110] = (dA[10] * dA[9]) / 2e0;
	v[105] = Power(dA[10], 2);
	v[117] = (dA[10] * dA[11]) / 2e0;
	v[115] = (dA[11] * dA[9]) / 2e0;
	v[106] = Power(dA[11], 2);
	v[290] = v[105] + v[106];
	v[193] = (*radA)*cos(c[1]);
	v[192] = (*radA)*sin(c[1]);
	v[162] = Power(dB[3], 2);
	v[160] = (dB[3] * dB[4]) / 2e0;
	v[155] = Power(dB[4], 2);
	v[167] = (dB[4] * dB[5]) / 2e0;
	v[165] = (dB[3] * dB[5]) / 2e0;
	v[156] = Power(dB[5], 2);
	v[291] = v[155] + v[156];
	v[209] = (*radB)*cos(c[3]);
	v[206] = (*radB)*sin(c[3]);
	v[187] = cpointB[0] + v[209];
	v[304] = v[187] * (*xfac);
	v[292] = v[187] * v[215];
	v[216] = -(v[187] * v[303]);
	v[344] = 0e0;
	v[345] = 0e0;
	v[346] = v[216];
	v[347] = -(v[206] * v[294]);
	v[217] = -(v[292] * (*zfac));
	v[348] = 0e0;
	v[349] = 0e0;
	v[350] = v[217];
	v[351] = v[206] * v[293];
	v[85] = 4e0 / (4e0 + v[289] + v[93]);
	v[88] = 1e0 - (v[289] * v[85]) / 2e0;
	v[89] = v[85] * (-dA[5] + v[91]);
	v[90] = v[85] * (dA[4] + v[96]);
	v[92] = v[85] * (dA[5] + v[91]);
	v[94] = 1e0 - (v[85] * (v[87] + v[93])) / 2e0;
	v[95] = v[85] * (-dA[3] + v[98]);
	v[97] = v[85] * (-dA[4] + v[96]);
	v[99] = v[85] * (dA[3] + v[98]);
	v[100] = 1e0 - (v[85] * (v[86] + v[93])) / 2e0;
	v[104] = 4e0 / (4e0 + v[112] + v[290]);
	v[107] = 1e0 - (v[104] * v[290]) / 2e0;
	v[108] = v[104] * (-dA[11] + v[110]);
	v[109] = v[104] * (dA[10] + v[115]);
	v[111] = v[104] * (dA[11] + v[110]);
	v[113] = 1e0 - (v[104] * (v[106] + v[112])) / 2e0;
	v[114] = v[104] * (-dA[9] + v[117]);
	v[116] = v[104] * (-dA[10] + v[115]);
	v[118] = v[104] * (dA[9] + v[117]);
	v[119] = 1e0 - (v[104] * (v[105] + v[112])) / 2e0;
	v[123] = QAAi[0][0] * v[88] + QAAi[1][0] * v[89] + QAAi[2][0] * v[90];
	v[299] = v[123] * v[192];
	v[124] = QAAi[0][1] * v[88] + QAAi[1][1] * v[89] + QAAi[2][1] * v[90];
	v[126] = QAAi[0][0] * v[92] + QAAi[1][0] * v[94] + QAAi[2][0] * v[95];
	v[127] = QAAi[0][1] * v[92] + QAAi[1][1] * v[94] + QAAi[2][1] * v[95];
	v[129] = QAAi[2][0] * v[100] + QAAi[0][0] * v[97] + QAAi[1][0] * v[99];
	v[130] = QAAi[2][1] * v[100] + QAAi[0][1] * v[97] + QAAi[1][1] * v[99];
	v[132] = QBAi[0][0] * v[107] + QBAi[1][0] * v[108] + QBAi[2][0] * v[109];
	v[298] = -(v[132] * v[192]);
	v[133] = QBAi[0][1] * v[107] + QBAi[1][1] * v[108] + QBAi[2][1] * v[109];
	v[135] = QBAi[0][0] * v[111] + QBAi[1][0] * v[113] + QBAi[2][0] * v[114];
	v[136] = QBAi[0][1] * v[111] + QBAi[1][1] * v[113] + QBAi[2][1] * v[114];
	v[138] = QBAi[0][0] * v[116] + QBAi[1][0] * v[118] + QBAi[2][0] * v[119];
	v[139] = QBAi[0][1] * v[116] + QBAi[1][1] * v[118] + QBAi[2][1] * v[119];
	v[147] = (1e0 - c[0]) / 2e0;
	v[148] = (1e0 + c[0]) / 2e0;
	v[149] = cpointA[0] + v[193];
	v[150] = cpointA[1] + v[192];
	v[204] = dA[8] + v[138] * v[149] + v[139] * v[150] + xBAi[2];
	v[203] = dA[2] + v[129] * v[149] + v[130] * v[150] + xAAi[2];
	v[300] = -v[203] + v[204];
	v[201] = dA[7] + v[135] * v[149] + v[136] * v[150] + xBAi[1];
	v[200] = dA[1] + v[126] * v[149] + v[127] * v[150] + xAAi[1];
	v[301] = -v[200] + v[201];
	v[198] = dA[6] + v[132] * v[149] + v[133] * v[150] + xBAi[0];
	v[197] = dA[0] + v[123] * v[149] + v[124] * v[150] + xAAi[0];
	v[302] = -v[197] + v[198];
	v[154] = 4e0 / (4e0 + v[162] + v[291]);
	v[157] = 1e0 - (v[154] * v[291]) / 2e0;
	v[158] = v[154] * (-dB[5] + v[160]);
	v[159] = v[154] * (dB[4] + v[165]);
	v[161] = v[154] * (dB[5] + v[160]);
	v[163] = 1e0 - (v[154] * (v[156] + v[162])) / 2e0;
	v[164] = v[154] * (-dB[3] + v[167]);
	v[166] = v[154] * (-dB[4] + v[165]);
	v[168] = v[154] * (dB[3] + v[167]);
	v[169] = 1e0 - (v[154] * (v[155] + v[162])) / 2e0;
	v[173] = QABi[0][0] * v[157] + QABi[1][0] * v[158] + QABi[2][0] * v[159];
	v[174] = QABi[0][1] * v[157] + QABi[1][1] * v[158] + QABi[2][1] * v[159];
	v[360] = v[302] / 2e0;
	v[361] = v[148] * (v[133] * v[193] + v[298]) + v[147] * (v[124] * v[193] - v[299]);
	v[362] = 0e0;
	v[363] = -(v[174] * v[209]);
	v[175] = QABi[0][2] * v[157] + QABi[1][2] * v[158] + QABi[2][2] * v[159];
	v[176] = QABi[0][0] * v[161] + QABi[1][0] * v[163] + QABi[2][0] * v[164];
	v[177] = QABi[0][1] * v[161] + QABi[1][1] * v[163] + QABi[2][1] * v[164];
	v[356] = v[301] / 2e0;
	v[357] = v[147] * (-(v[126] * v[192]) + v[127] * v[193]) + v[148] * (-(v[135] * v[192]) + v[136] * v[193]);
	v[358] = 0e0;
	v[359] = -(v[177] * v[209]);
	v[178] = QABi[0][2] * v[161] + QABi[1][2] * v[163] + QABi[2][2] * v[164];
	v[179] = QABi[0][0] * v[166] + QABi[1][0] * v[168] + QABi[2][0] * v[169];
	v[180] = QABi[0][1] * v[166] + QABi[1][1] * v[168] + QABi[2][1] * v[169];
	v[352] = v[300] / 2e0;
	v[353] = v[147] * (-(v[129] * v[192]) + v[130] * v[193]) + v[148] * (-(v[138] * v[192]) + v[139] * v[193]);
	v[354] = 0e0;
	v[355] = -(v[180] * v[209]);
	v[181] = QABi[0][2] * v[166] + QABi[1][2] * v[168] + QABi[2][2] * v[169];
	v[185] = v[292] * (*xfac);
	v[186] = cpointB[1] + v[206];
	v[188] = -(v[187] * v[293]);
	v[242] = -dB[0] - v[173] * v[185] - v[174] * v[186] - v[175] * v[188] + v[147] * v[197] + v[148] * v[198] - xABi[0];
	v[243] = -dB[1] - v[176] * v[185] - v[177] * v[186] - v[178] * v[188] + v[147] * v[200] + v[148] * v[201] - xABi[1];
	v[248] = -dB[2] - v[179] * v[185] - v[180] * v[186] - v[181] * v[188] + v[147] * v[203] + v[148] * v[204] - xABi[2];
	v[311] = v[192] * ((v[126] - v[135])*v[243] + (v[129] - v[138])*v[248]) + v[193] * ((-v[124] + v[133])*v[242] + (
		-v[127] + v[136])*v[243] + (-v[130] + v[139])*v[248]) + v[242] * v[298] + v[242] * v[299];
	v[282] = -(v[174] * v[242]) - v[177] * v[243] - v[180] * v[248];
	v[279] = (v[123] * v[147] + v[132] * v[148])*v[242] + (v[126] * v[147] + v[135] * v[148])*v[243] + (v[129] * v[147]
		+ v[138] * v[148])*v[248];
	v[278] = (v[124] * v[147] + v[133] * v[148])*v[242] + (v[127] * v[147] + v[136] * v[148])*v[243] + (v[130] * v[147]
		+ v[139] * v[148])*v[248];
	v[249] = -(v[175] * v[242]) - v[178] * v[243] - v[181] * v[248];
	v[305] = v[249] * (*zfac);
	v[250] = -(v[173] * v[242]) - v[176] * v[243] - v[179] * v[248];
	v[314] = v[206] * (v[250] * v[303] + v[215] * v[305]);
	v[420] = 0e0;
	v[421] = 0e0;
	v[422] = -(v[187] * v[305]);
	v[423] = -(v[206] * v[250] * (*xfac));
	v[416] = 0e0;
	v[417] = 0e0;
	v[418] = -(v[250] * v[304]);
	v[419] = v[206] * v[305];
	v[281] = -(v[249] * v[293]) + v[250] * v[294];
	v[315] = -(v[209] * v[281]) - v[206] * v[282];
	v[336] = (-(v[197] * v[242]) + v[198] * v[242] - v[200] * v[243] + v[201] * v[243] - v[203] * v[248] + v[204] * v[248])
		/ 2e0;
	v[337] = v[193] * v[278] - v[192] * v[279];
	v[338] = v[217] * v[249] + v[216] * v[250];
	v[339] = -(v[206] * v[281]) + v[209] * v[282];
	for (i246 = 1; i246 <= 4; i246++) {
		b297 = i246 == 2;
		i295 = (i246 == 1 ? 1 : 0);
		v[296] = -i295 / 2e0;
		v[275] = -(i295*v[242]) / 2e0;
		v[270] = v[243] * v[296];
		v[266] = v[248] * v[296];
		v[257] = v[343 + i246];
		v[258] = v[347 + i246];
		v[260] = -(v[179] * v[257]) - v[181] * v[258] + v[351 + i246];
		v[262] = -(v[176] * v[257]) - v[178] * v[258] + v[355 + i246];
		v[264] = -(v[173] * v[257]) - v[175] * v[258] + v[359 + i246];
		v[265] = v[147] * v[260] + v[266];
		v[267] = v[148] * v[260] - v[266];
		v[269] = v[147] * v[262] + v[270];
		v[271] = v[148] * v[262] - v[270];
		v[272] = -(v[181] * v[260]) - v[178] * v[262] - v[175] * v[264];
		v[273] = -(v[179] * v[260]) - v[176] * v[262] - v[173] * v[264];
		v[274] = v[147] * v[264] + v[275];
		v[276] = v[148] * v[264] - v[275];
		v[424] = ((b297 ? v[311] : 0e0) + v[260] * v[300] + v[262] * v[301] + v[264] * v[302]) / 2e0;
		v[425] = -(v[192] * ((b297 ? v[278] : 0e0) + v[129] * v[265] + v[138] * v[267] + v[126] * v[269] + v[135] * v[271]
			+ v[123] * v[274] + v[132] * v[276])) + v[193] * ((b297 ? -v[279] : 0e0) + v[130] * v[265] + v[139] * v[267]
				+ v[127] * v[269] + v[136] * v[271] + v[124] * v[274] + v[133] * v[276]);
		v[426] = -(v[210] * (v[273] * v[304] + v[419 + i246])) + v[215] * (v[415 + i246] - v[187] * v[272] * (*zfac));
		v[427] = (i246 == 3 ? v[314] : 0e0) + (i246 == 4 ? v[315] : 0e0) - v[209] * (v[180] * v[260] + v[177] * v[262] + v[174] * v[264]
			) + v[206] * (v[272] * v[293] - v[273] * v[294]);
		for (i254 = i246; i254 <= 4; i254++) {
			v[284] = v[423 + i254];
			Hes[i246 - 1][i254 - 1] = v[284];
			if (i246 != i254) {
				Hes[i254 - 1][i246 - 1] = v[284];
			}
			else {
			};
		};/* end for */
	};/* end for */

	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
			mHes(i, j) = Hes[i][j];

	////////////////////////////////////////////////Desalocando memória//////////////////////////////////////////////////////
	for (int i = 0; i < 3; i++)
	{
		delete[] QAAi[i];
		delete[] QBAi[i];
		delete[] QABi[i];
	}
	delete[] QAAi;
	delete[] QBAi;
	delete[] QABi;
}

void FlexibleArcExtrusion_1_RigidArcRevolution_1::ContactSS(bool *stick, bool *stickupdated, bool *previouscontact, double* Rc, double** Kc, double** invH, double* convective, double* copy_convective, double* gti, double* gtpupdated, double* epsn, double* epsn0, double* epst, double* cn, double* ct, double* mus, double* mud, double* fn, double* ft)
{
	double v[30000];		//variavel temporaria - AceGen
	double value = 0.0;
	double* a4;
	double* a5;
	double* a6;
	if (typeid(*db.solution[db.current_solution_number - 1]) == typeid(Dynamic))
	{
		Dynamic* ptr_sol = static_cast<Dynamic*>(db.solution[db.current_solution_number - 1]);
		a4 = &ptr_sol->a4;
		a5 = &ptr_sol->a5;
		a6 = &ptr_sol->a6;
	}
	else
	{
		a4 = &value;
		a5 = &value;
		a6 = &value;
	}
	//Zerando matrizes e vetores
	for (int i = 0; i < 18; i++)
	{
		Rc[i] = 0.0;
		for (int j = 0; j < 18; j++)
			Kc[i][j] = 0.0;
	}

	///////////////////////////////Ponteiros e variaveis para facilitar acesso//////////////////////////////////////////
	FlexibleArcExtrusion_1* surf1;		//Ponteiro para a superficie 1
	RigidArcRevolution_1* surf2;		//Ponteiro para a superficie 2
	double* radA;
	double* radB;
	double* xfac;
	double* zfac;
	double* cpointA;
	double* cpointB;
	double* dA;
	double* dB;
	double* duiA;
	double* dduiA;
	double* duiB;
	double* dduiB;
	bool* normalintA;
	bool* normalintB;
	double* xAAi;
	double* xBAi;
	double* xABi;
	double** QAAi;
	double** QBAi;
	double** QABi;

	surf1 = static_cast<FlexibleArcExtrusion_1*>(db.surfaces[surf1_ID - 1]);
	surf2 = static_cast<RigidArcRevolution_1*>(db.surfaces[surf2_ID - 1]);
	radA = surf1->radius;
	radB = surf2->radius;
	xfac = &surf2->x_fac;
	zfac = &surf2->z_fac;
	cpointA = surf1->c_point->getMatrix();
	cpointB = surf2->c_point->getMatrix();
	dA = surf1->d_A->getMatrix();
	dB = surf2->d_B->getMatrix();
	duiA = surf1->dui_A->getMatrix();
	dduiA = surf1->ddui_A->getMatrix();
	duiB = surf2->dui_B->getMatrix();
	dduiB = surf2->ddui_B->getMatrix();
	xAAi = surf1->x_AAi->getMatrix();
	xBAi = surf1->x_BAi->getMatrix();
	xABi = surf2->x_ABi->getMatrix();
	normalintA = &surf1->flag_normal_int;
	normalintB = &surf2->flag_normal_int;
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	QAAi = new double*[3];
	QBAi = new double*[3];
	QABi = new double*[3];
	for (int i = 0; i < 3; i++)
	{
		QAAi[i] = new double[3];
		QBAi[i] = new double[3];
		QABi[i] = new double[3];
	}
	//Salvando variaveis locais para montagem de superficies
	surf1->Q_AAi->MatrixToPtr(QAAi, 3);
	surf1->Q_BAi->MatrixToPtr(QBAi, 3);
	surf2->Q_ABi->MatrixToPtr(QABi, 3);
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	
	double* ci = copy_convective;
	double* cp = convective;


#pragma region AceGen
	double v01; double v010; double v011; double v012; double v013; double v014;
	double v015; double v016; double v017; double v018; double v019; double v02;
	double v020; double v021; double v022; double v023; double v024; double v025;
	double v026; double v027; double v03; double v04; double v05; double v06; double v07;
	double v08; double v09;
	int i1087, i1169, i1514, i1876, i3049, i3128, i3711, i3712, i3713, i3714, i3715, i3716
		, i3717, i3718, i3719, i3789, i3790, i3791, i3792, i3793, i3794, i3795, i3796, i3797, i3798
		, i3799, i3800, i3801, i3802, i3803, i3813, i3814, i3815, i4016, i4017, i4018, i4019, i4020
		, i4021, i4022, i4023, i4024, b909, b1017, b1051
		, b1519, b1520, b1521, b1565, b1608, b2238, b2283, b2366, b2367, b2384, b2547
		, b2578;
	v[3734] = -((*a4)*(*ct));
	v[1438] = cos(cp[2]);
	v[3704] = v[1438] * (*xfac);
	v[992] = sin(cp[2]);
	v[3633] = v[992] * (*zfac);
	v[427] = dA[3] / 2e0;
	v[425] = 2e0*dA[3];
	v[154] = Power(dA[3], 2);
	v[428] = 2e0*dA[4];
	v[5293] = 0e0;
	v[5294] = 0e0;
	v[5295] = 0e0;
	v[5296] = v[425];
	v[5297] = v[428];
	v[5298] = 0e0;
	v[5299] = 0e0;
	v[5300] = 0e0;
	v[5301] = 0e0;
	v[5302] = 0e0;
	v[5303] = 0e0;
	v[5304] = 0e0;
	v[5305] = 0e0;
	v[5306] = 0e0;
	v[5307] = 0e0;
	v[5308] = 0e0;
	v[5309] = 0e0;
	v[5310] = 0e0;
	v[5985] = 0e0;
	v[5986] = 0e0;
	v[5987] = 0e0;
	v[5988] = -v[425];
	v[5989] = -v[428];
	v[5990] = 0e0;
	v[5991] = 0e0;
	v[5992] = 0e0;
	v[5993] = 0e0;
	v[5994] = 0e0;
	v[5995] = 0e0;
	v[5996] = 0e0;
	v[5997] = 0e0;
	v[5998] = 0e0;
	v[5999] = 0e0;
	v[6000] = 0e0;
	v[6001] = 0e0;
	v[6002] = 0e0;
	v[426] = dA[4] / 2e0;
	v[4951] = 0e0;
	v[4952] = 0e0;
	v[4953] = 0e0;
	v[4954] = v[426];
	v[4955] = v[427];
	v[4956] = 0e0;
	v[4957] = 0e0;
	v[4958] = 0e0;
	v[4959] = 0e0;
	v[4960] = 0e0;
	v[4961] = 0e0;
	v[4962] = 0e0;
	v[4963] = 0e0;
	v[4964] = 0e0;
	v[4965] = 0e0;
	v[4966] = 0e0;
	v[4967] = 0e0;
	v[4968] = 0e0;
	v[152] = (dA[3] * dA[4]) / 2e0;
	v[147] = Power(dA[4], 2);
	v[483] = -v[147] - v[154];
	v[461] = dA[5] + v[152];
	v[451] = -dA[5] + v[152];
	v[430] = 2e0*dA[5];
	v[7901] = 0e0;
	v[7902] = 0e0;
	v[7903] = 0e0;
	v[7904] = -v[425];
	v[7905] = 0e0;
	v[7906] = -v[430];
	v[7907] = 0e0;
	v[7908] = 0e0;
	v[7909] = 0e0;
	v[7910] = 0e0;
	v[7911] = 0e0;
	v[7912] = 0e0;
	v[7913] = 0e0;
	v[7914] = 0e0;
	v[7915] = 0e0;
	v[7916] = 0e0;
	v[7917] = 0e0;
	v[7918] = 0e0;
	v[5221] = 0e0;
	v[5222] = 0e0;
	v[5223] = 0e0;
	v[5224] = v[425];
	v[5225] = 0e0;
	v[5226] = v[430];
	v[5227] = 0e0;
	v[5228] = 0e0;
	v[5229] = 0e0;
	v[5230] = 0e0;
	v[5231] = 0e0;
	v[5232] = 0e0;
	v[5233] = 0e0;
	v[5234] = 0e0;
	v[5235] = 0e0;
	v[5236] = 0e0;
	v[5237] = 0e0;
	v[5238] = 0e0;
	v[5149] = 0e0;
	v[5150] = 0e0;
	v[5151] = 0e0;
	v[5152] = 0e0;
	v[5153] = v[428];
	v[5154] = v[430];
	v[5155] = 0e0;
	v[5156] = 0e0;
	v[5157] = 0e0;
	v[5158] = 0e0;
	v[5159] = 0e0;
	v[5160] = 0e0;
	v[5161] = 0e0;
	v[5162] = 0e0;
	v[5163] = 0e0;
	v[5164] = 0e0;
	v[5165] = 0e0;
	v[5166] = 0e0;
	v[5949] = 0e0;
	v[5950] = 0e0;
	v[5951] = 0e0;
	v[5952] = 0e0;
	v[5953] = -v[428];
	v[5954] = -v[430];
	v[5955] = 0e0;
	v[5956] = 0e0;
	v[5957] = 0e0;
	v[5958] = 0e0;
	v[5959] = 0e0;
	v[5960] = 0e0;
	v[5961] = 0e0;
	v[5962] = 0e0;
	v[5963] = 0e0;
	v[5964] = 0e0;
	v[5965] = 0e0;
	v[5966] = 0e0;
	v[4933] = 0e0;
	v[4934] = 0e0;
	v[4935] = 0e0;
	v[4936] = v[425];
	v[4937] = v[428];
	v[4938] = v[430];
	v[4939] = 0e0;
	v[4940] = 0e0;
	v[4941] = 0e0;
	v[4942] = 0e0;
	v[4943] = 0e0;
	v[4944] = 0e0;
	v[4945] = 0e0;
	v[4946] = 0e0;
	v[4947] = 0e0;
	v[4948] = 0e0;
	v[4949] = 0e0;
	v[4950] = 0e0;
	v[429] = dA[5] / 2e0;
	v[4969] = 0e0;
	v[4970] = 0e0;
	v[4971] = 0e0;
	v[4972] = v[429];
	v[4973] = 0e0;
	v[4974] = v[427];
	v[4975] = 0e0;
	v[4976] = 0e0;
	v[4977] = 0e0;
	v[4978] = 0e0;
	v[4979] = 0e0;
	v[4980] = 0e0;
	v[4981] = 0e0;
	v[4982] = 0e0;
	v[4983] = 0e0;
	v[4984] = 0e0;
	v[4985] = 0e0;
	v[4986] = 0e0;
	v[4987] = 0e0;
	v[4988] = 0e0;
	v[4989] = 0e0;
	v[4990] = 0e0;
	v[4991] = v[429];
	v[4992] = v[426];
	v[4993] = 0e0;
	v[4994] = 0e0;
	v[4995] = 0e0;
	v[4996] = 0e0;
	v[4997] = 0e0;
	v[4998] = 0e0;
	v[4999] = 0e0;
	v[5000] = 0e0;
	v[5001] = 0e0;
	v[5002] = 0e0;
	v[5003] = 0e0;
	v[5004] = 0e0;
	v[159] = (dA[4] * dA[5]) / 2e0;
	v[478] = dA[3] + v[159];
	v[470] = -dA[3] + v[159];
	v[157] = (dA[3] * dA[5]) / 2e0;
	v[474] = -dA[4] + v[157];
	v[457] = dA[4] + v[157];
	v[148] = Power(dA[5], 2);
	v[1402] = 4e0 + v[147] + v[148] + v[154];
	v[4256] = 1e0 / Power(v[1402], 4);
	v[1850] = 1e0 / Power(v[1402], 3);
	v[1953] = -8e0*v[1850] * v[425];
	v[1951] = -8e0*v[1850] * v[430];
	v[1949] = 8e0*v[1850] * v[428];
	v[1224] = 1e0 / Power(v[1402], 2);
	v[4003] = 4e0*v[1224];
	v[465] = -v[148] - v[154];
	v[446] = -v[147] - v[148];
	v[3188] = -2e0*v[1224] * v[446];
	v[445] = 4e0*v[1224] * v[430];
	v[487] = -(v[445] * v[483]) / 2e0;
	v[444] = -4e0*v[1224] * v[428];
	v[467] = (v[444] * v[465]) / 2e0;
	v[443] = 4e0*v[1224] * v[425];
	v[447] = -(v[443] * v[446]) / 2e0;
	v[433] = dA[9] / 2e0;
	v[431] = 2e0*dA[9];
	v[173] = Power(dA[9], 2);
	v[434] = 2e0*dA[10];
	v[5455] = 0e0;
	v[5456] = 0e0;
	v[5457] = 0e0;
	v[5458] = 0e0;
	v[5459] = 0e0;
	v[5460] = 0e0;
	v[5461] = 0e0;
	v[5462] = 0e0;
	v[5463] = 0e0;
	v[5464] = v[431];
	v[5465] = v[434];
	v[5466] = 0e0;
	v[5467] = 0e0;
	v[5468] = 0e0;
	v[5469] = 0e0;
	v[5470] = 0e0;
	v[5471] = 0e0;
	v[5472] = 0e0;
	v[6093] = 0e0;
	v[6094] = 0e0;
	v[6095] = 0e0;
	v[6096] = 0e0;
	v[6097] = 0e0;
	v[6098] = 0e0;
	v[6099] = 0e0;
	v[6100] = 0e0;
	v[6101] = 0e0;
	v[6102] = -v[431];
	v[6103] = -v[434];
	v[6104] = 0e0;
	v[6105] = 0e0;
	v[6106] = 0e0;
	v[6107] = 0e0;
	v[6108] = 0e0;
	v[6109] = 0e0;
	v[6110] = 0e0;
	v[432] = dA[10] / 2e0;
	v[5023] = 0e0;
	v[5024] = 0e0;
	v[5025] = 0e0;
	v[5026] = 0e0;
	v[5027] = 0e0;
	v[5028] = 0e0;
	v[5029] = 0e0;
	v[5030] = 0e0;
	v[5031] = 0e0;
	v[5032] = v[432];
	v[5033] = v[433];
	v[5034] = 0e0;
	v[5035] = 0e0;
	v[5036] = 0e0;
	v[5037] = 0e0;
	v[5038] = 0e0;
	v[5039] = 0e0;
	v[5040] = 0e0;
	v[171] = (dA[10] * dA[9]) / 2e0;
	v[166] = Power(dA[10], 2);
	v[528] = -v[166] - v[173];
	v[506] = dA[11] + v[171];
	v[496] = -dA[11] + v[171];
	v[436] = 2e0*dA[11];
	v[7919] = 0e0;
	v[7920] = 0e0;
	v[7921] = 0e0;
	v[7922] = 0e0;
	v[7923] = 0e0;
	v[7924] = 0e0;
	v[7925] = 0e0;
	v[7926] = 0e0;
	v[7927] = 0e0;
	v[7928] = -v[431];
	v[7929] = 0e0;
	v[7930] = -v[436];
	v[7931] = 0e0;
	v[7932] = 0e0;
	v[7933] = 0e0;
	v[7934] = 0e0;
	v[7935] = 0e0;
	v[7936] = 0e0;
	v[5383] = 0e0;
	v[5384] = 0e0;
	v[5385] = 0e0;
	v[5386] = 0e0;
	v[5387] = 0e0;
	v[5388] = 0e0;
	v[5389] = 0e0;
	v[5390] = 0e0;
	v[5391] = 0e0;
	v[5392] = v[431];
	v[5393] = 0e0;
	v[5394] = v[436];
	v[5395] = 0e0;
	v[5396] = 0e0;
	v[5397] = 0e0;
	v[5398] = 0e0;
	v[5399] = 0e0;
	v[5400] = 0e0;
	v[5311] = 0e0;
	v[5312] = 0e0;
	v[5313] = 0e0;
	v[5314] = 0e0;
	v[5315] = 0e0;
	v[5316] = 0e0;
	v[5317] = 0e0;
	v[5318] = 0e0;
	v[5319] = 0e0;
	v[5320] = 0e0;
	v[5321] = v[434];
	v[5322] = v[436];
	v[5323] = 0e0;
	v[5324] = 0e0;
	v[5325] = 0e0;
	v[5326] = 0e0;
	v[5327] = 0e0;
	v[5328] = 0e0;
	v[6057] = 0e0;
	v[6058] = 0e0;
	v[6059] = 0e0;
	v[6060] = 0e0;
	v[6061] = 0e0;
	v[6062] = 0e0;
	v[6063] = 0e0;
	v[6064] = 0e0;
	v[6065] = 0e0;
	v[6066] = 0e0;
	v[6067] = -v[434];
	v[6068] = -v[436];
	v[6069] = 0e0;
	v[6070] = 0e0;
	v[6071] = 0e0;
	v[6072] = 0e0;
	v[6073] = 0e0;
	v[6074] = 0e0;
	v[5005] = 0e0;
	v[5006] = 0e0;
	v[5007] = 0e0;
	v[5008] = 0e0;
	v[5009] = 0e0;
	v[5010] = 0e0;
	v[5011] = 0e0;
	v[5012] = 0e0;
	v[5013] = 0e0;
	v[5014] = v[431];
	v[5015] = v[434];
	v[5016] = v[436];
	v[5017] = 0e0;
	v[5018] = 0e0;
	v[5019] = 0e0;
	v[5020] = 0e0;
	v[5021] = 0e0;
	v[5022] = 0e0;
	v[435] = dA[11] / 2e0;
	v[5041] = 0e0;
	v[5042] = 0e0;
	v[5043] = 0e0;
	v[5044] = 0e0;
	v[5045] = 0e0;
	v[5046] = 0e0;
	v[5047] = 0e0;
	v[5048] = 0e0;
	v[5049] = 0e0;
	v[5050] = v[435];
	v[5051] = 0e0;
	v[5052] = v[433];
	v[5053] = 0e0;
	v[5054] = 0e0;
	v[5055] = 0e0;
	v[5056] = 0e0;
	v[5057] = 0e0;
	v[5058] = 0e0;
	v[5059] = 0e0;
	v[5060] = 0e0;
	v[5061] = 0e0;
	v[5062] = 0e0;
	v[5063] = 0e0;
	v[5064] = 0e0;
	v[5065] = 0e0;
	v[5066] = 0e0;
	v[5067] = 0e0;
	v[5068] = 0e0;
	v[5069] = v[435];
	v[5070] = v[432];
	v[5071] = 0e0;
	v[5072] = 0e0;
	v[5073] = 0e0;
	v[5074] = 0e0;
	v[5075] = 0e0;
	v[5076] = 0e0;
	v[178] = (dA[10] * dA[11]) / 2e0;
	v[523] = dA[9] + v[178];
	v[515] = -dA[9] + v[178];
	v[176] = (dA[11] * dA[9]) / 2e0;
	v[519] = -dA[10] + v[176];
	v[502] = dA[10] + v[176];
	v[167] = Power(dA[11], 2);
	v[1412] = 4e0 + v[166] + v[167] + v[173];
	v[4250] = 1e0 / Power(v[1412], 4);
	v[1852] = 1e0 / Power(v[1412], 3);
	v[1980] = -8e0*v[1852] * v[431];
	v[1978] = -8e0*v[1852] * v[436];
	v[1976] = 8e0*v[1852] * v[434];
	v[1238] = 1e0 / Power(v[1412], 2);
	v[4007] = 4e0*v[1238];
	v[510] = -v[167] - v[173];
	v[491] = -v[166] - v[167];
	v[3199] = -2e0*v[1238] * v[491];
	v[490] = 4e0*v[1238] * v[436];
	v[532] = -(v[490] * v[528]) / 2e0;
	v[489] = -4e0*v[1238] * v[434];
	v[512] = (v[489] * v[510]) / 2e0;
	v[488] = 4e0*v[1238] * v[431];
	v[492] = -(v[488] * v[491]) / 2e0;
	v[288] = (*a4)*dA[3] + (*a6)*dduiA[3] + (*a5)*duiA[3];
	v[294] = (*a4)*dA[4] + (*a6)*dduiA[4] + (*a5)*duiA[4];
	v[296] = (*a4)*dA[5] + (*a6)*dduiA[5] + (*a5)*duiA[5];
	v[314] = (*a4)*dA[9] + (*a6)*dduiA[9] + (*a5)*duiA[9];
	v[320] = (*a4)*dA[10] + (*a6)*dduiA[10] + (*a5)*duiA[10];
	v[322] = (*a4)*dA[11] + (*a6)*dduiA[11] + (*a5)*duiA[11];
	v[1439] = (*radA)*cos(cp[1]);
	v[993] = (*radA)*sin(cp[1]);
	v[1336] = -((*epsn)*v[993]) / 2e0;
	v[1335] = ((*epsn)*v[1439]) / 2e0;
	v[439] = dB[3] / 2e0;
	v[437] = 2e0*dB[3];
	v[233] = Power(dB[3], 2);
	v[440] = 2e0*dB[4];
	v[5617] = 0e0;
	v[5618] = 0e0;
	v[5619] = 0e0;
	v[5620] = 0e0;
	v[5621] = 0e0;
	v[5622] = 0e0;
	v[5623] = 0e0;
	v[5624] = 0e0;
	v[5625] = 0e0;
	v[5626] = 0e0;
	v[5627] = 0e0;
	v[5628] = 0e0;
	v[5629] = 0e0;
	v[5630] = 0e0;
	v[5631] = 0e0;
	v[5632] = v[437];
	v[5633] = v[440];
	v[5634] = 0e0;
	v[6309] = 0e0;
	v[6310] = 0e0;
	v[6311] = 0e0;
	v[6312] = 0e0;
	v[6313] = 0e0;
	v[6314] = 0e0;
	v[6315] = 0e0;
	v[6316] = 0e0;
	v[6317] = 0e0;
	v[6318] = 0e0;
	v[6319] = 0e0;
	v[6320] = 0e0;
	v[6321] = 0e0;
	v[6322] = 0e0;
	v[6323] = 0e0;
	v[6324] = -v[437];
	v[6325] = -v[440];
	v[6326] = 0e0;
	v[438] = dB[4] / 2e0;
	v[5095] = 0e0;
	v[5096] = 0e0;
	v[5097] = 0e0;
	v[5098] = 0e0;
	v[5099] = 0e0;
	v[5100] = 0e0;
	v[5101] = 0e0;
	v[5102] = 0e0;
	v[5103] = 0e0;
	v[5104] = 0e0;
	v[5105] = 0e0;
	v[5106] = 0e0;
	v[5107] = 0e0;
	v[5108] = 0e0;
	v[5109] = 0e0;
	v[5110] = v[438];
	v[5111] = v[439];
	v[5112] = 0e0;
	v[231] = (dB[3] * dB[4]) / 2e0;
	v[226] = Power(dB[4], 2);
	v[699] = -v[226] - v[233];
	v[677] = dB[5] + v[231];
	v[667] = -dB[5] + v[231];
	v[442] = 2e0*dB[5];
	v[7937] = 0e0;
	v[7938] = 0e0;
	v[7939] = 0e0;
	v[7940] = 0e0;
	v[7941] = 0e0;
	v[7942] = 0e0;
	v[7943] = 0e0;
	v[7944] = 0e0;
	v[7945] = 0e0;
	v[7946] = 0e0;
	v[7947] = 0e0;
	v[7948] = 0e0;
	v[7949] = 0e0;
	v[7950] = 0e0;
	v[7951] = 0e0;
	v[7952] = -v[437];
	v[7953] = 0e0;
	v[7954] = -v[442];
	v[5545] = 0e0;
	v[5546] = 0e0;
	v[5547] = 0e0;
	v[5548] = 0e0;
	v[5549] = 0e0;
	v[5550] = 0e0;
	v[5551] = 0e0;
	v[5552] = 0e0;
	v[5553] = 0e0;
	v[5554] = 0e0;
	v[5555] = 0e0;
	v[5556] = 0e0;
	v[5557] = 0e0;
	v[5558] = 0e0;
	v[5559] = 0e0;
	v[5560] = v[437];
	v[5561] = 0e0;
	v[5562] = v[442];
	v[5473] = 0e0;
	v[5474] = 0e0;
	v[5475] = 0e0;
	v[5476] = 0e0;
	v[5477] = 0e0;
	v[5478] = 0e0;
	v[5479] = 0e0;
	v[5480] = 0e0;
	v[5481] = 0e0;
	v[5482] = 0e0;
	v[5483] = 0e0;
	v[5484] = 0e0;
	v[5485] = 0e0;
	v[5486] = 0e0;
	v[5487] = 0e0;
	v[5488] = 0e0;
	v[5489] = v[440];
	v[5490] = v[442];
	v[6273] = 0e0;
	v[6274] = 0e0;
	v[6275] = 0e0;
	v[6276] = 0e0;
	v[6277] = 0e0;
	v[6278] = 0e0;
	v[6279] = 0e0;
	v[6280] = 0e0;
	v[6281] = 0e0;
	v[6282] = 0e0;
	v[6283] = 0e0;
	v[6284] = 0e0;
	v[6285] = 0e0;
	v[6286] = 0e0;
	v[6287] = 0e0;
	v[6288] = 0e0;
	v[6289] = -v[440];
	v[6290] = -v[442];
	v[5077] = 0e0;
	v[5078] = 0e0;
	v[5079] = 0e0;
	v[5080] = 0e0;
	v[5081] = 0e0;
	v[5082] = 0e0;
	v[5083] = 0e0;
	v[5084] = 0e0;
	v[5085] = 0e0;
	v[5086] = 0e0;
	v[5087] = 0e0;
	v[5088] = 0e0;
	v[5089] = 0e0;
	v[5090] = 0e0;
	v[5091] = 0e0;
	v[5092] = v[437];
	v[5093] = v[440];
	v[5094] = v[442];
	v[441] = dB[5] / 2e0;
	v[5113] = 0e0;
	v[5114] = 0e0;
	v[5115] = 0e0;
	v[5116] = 0e0;
	v[5117] = 0e0;
	v[5118] = 0e0;
	v[5119] = 0e0;
	v[5120] = 0e0;
	v[5121] = 0e0;
	v[5122] = 0e0;
	v[5123] = 0e0;
	v[5124] = 0e0;
	v[5125] = 0e0;
	v[5126] = 0e0;
	v[5127] = 0e0;
	v[5128] = v[441];
	v[5129] = 0e0;
	v[5130] = v[439];
	v[5131] = 0e0;
	v[5132] = 0e0;
	v[5133] = 0e0;
	v[5134] = 0e0;
	v[5135] = 0e0;
	v[5136] = 0e0;
	v[5137] = 0e0;
	v[5138] = 0e0;
	v[5139] = 0e0;
	v[5140] = 0e0;
	v[5141] = 0e0;
	v[5142] = 0e0;
	v[5143] = 0e0;
	v[5144] = 0e0;
	v[5145] = 0e0;
	v[5146] = 0e0;
	v[5147] = v[441];
	v[5148] = v[438];
	v[238] = (dB[4] * dB[5]) / 2e0;
	v[694] = dB[3] + v[238];
	v[686] = -dB[3] + v[238];
	v[236] = (dB[3] * dB[5]) / 2e0;
	v[690] = -dB[4] + v[236];
	v[673] = dB[4] + v[236];
	v[227] = Power(dB[5], 2);
	v[1422] = 4e0 + v[226] + v[227] + v[233];
	v[4244] = 1e0 / Power(v[1422], 4);
	v[1854] = 1e0 / Power(v[1422], 3);
	v[2061] = -8e0*v[1854] * v[437];
	v[2059] = -8e0*v[1854] * v[442];
	v[2057] = 8e0*v[1854] * v[440];
	v[1252] = 1e0 / Power(v[1422], 2);
	v[4011] = 4e0*v[1252];
	v[681] = -v[227] - v[233];
	v[662] = -v[226] - v[227];
	v[3210] = -2e0*v[1252] * v[662];
	v[661] = 4e0*v[1252] * v[442];
	v[703] = -(v[661] * v[699]) / 2e0;
	v[660] = -4e0*v[1252] * v[440];
	v[683] = (v[660] * v[681]) / 2e0;
	v[659] = 4e0*v[1252] * v[437];
	v[663] = -(v[659] * v[662]) / 2e0;
	v[340] = (*a4)*dB[3] + (*a6)*dduiB[3] + (*a5)*duiB[3];
	v[346] = (*a4)*dB[4] + (*a6)*dduiB[4] + (*a5)*duiB[4];
	v[348] = (*a4)*dB[5] + (*a6)*dduiB[5] + (*a5)*duiB[5];
	v[1440] = (*radB)*cos(cp[3]);
	v[994] = (*radB)*sin(cp[3]);
	v[262] = cpointB[0] + v[1440];
	v[3918] = v[262] * (*zfac);
	v[3725] = v[262] * (*xfac);
	v[3632] = v[1438] * v[262];
	v[258] = cpointB[0] + (*radB)*cos(ci[3]);
	v[418] = -(v[3704] * v[994]);
	v[410] = -(v[3725] * v[992]);
	v[420] = v[3633] * v[994];
	v[412] = -(v[3632] * (*zfac));
	v[146] = 4e0 / v[1402];
	v[4001] = -v[146] / 2e0;
	v[485] = -(v[146] * v[428]) / 2e0;
	v[486] = (v[444] * v[483]) / 2e0 + v[485];
	v[482] = -(v[146] * v[425]) / 2e0;
	v[484] = v[482] - (v[443] * v[483]) / 2e0;
	v[479] = v[146] - v[443] * v[478];
	v[476] = -v[146] + v[444] * v[474];
	v[471] = -v[146] - v[443] * v[470];
	v[468] = -(v[146] * v[430]) / 2e0;
	v[469] = -(v[445] * v[465]) / 2e0 + v[468];
	v[466] = -(v[443] * v[465]) / 2e0 + v[482];
	v[464] = v[146] - v[445] * v[461];
	v[459] = v[146] + v[444] * v[457];
	v[456] = -(v[146] * v[429]);
	v[480] = -v[456] + v[444] * v[478];
	v[555] = QAAi[0][1] * v[476] + QAAi[1][1] * v[480] + QAAi[2][1] * v[486];
	v[552] = QAAi[0][0] * v[476] + QAAi[1][0] * v[480] + QAAi[2][0] * v[486];
	v[475] = -v[456] - v[443] * v[474];
	v[554] = QAAi[0][1] * v[475] + QAAi[1][1] * v[479] + QAAi[2][1] * v[484];
	v[551] = QAAi[0][0] * v[475] + QAAi[1][0] * v[479] + QAAi[2][0] * v[484];
	v[472] = -v[456] + v[444] * v[470];
	v[458] = -v[456] - v[443] * v[457];
	v[455] = -v[146] - v[445] * v[451];
	v[453] = -(v[146] * v[427]);
	v[477] = -v[453] - v[445] * v[474];
	v[463] = -v[453] + v[444] * v[461];
	v[546] = QAAi[0][1] * v[463] + QAAi[1][1] * v[467] + QAAi[2][1] * v[472];
	v[543] = QAAi[0][0] * v[463] + QAAi[1][0] * v[467] + QAAi[2][0] * v[472];
	v[460] = -v[453] - v[445] * v[457];
	v[454] = v[444] * v[451] - v[453];
	v[450] = v[146] * v[426];
	v[481] = v[450] - v[445] * v[478];
	v[556] = QAAi[0][1] * v[477] + QAAi[1][1] * v[481] + QAAi[2][1] * v[487];
	v[553] = QAAi[0][0] * v[477] + QAAi[1][0] * v[481] + QAAi[2][0] * v[487];
	v[473] = v[450] - v[445] * v[470];
	v[547] = QAAi[0][1] * v[464] + QAAi[1][1] * v[469] + QAAi[2][1] * v[473];
	v[544] = QAAi[0][0] * v[464] + QAAi[1][0] * v[469] + QAAi[2][0] * v[473];
	v[462] = v[450] - v[443] * v[461];
	v[545] = QAAi[0][1] * v[462] + QAAi[1][1] * v[466] + QAAi[2][1] * v[471];
	v[542] = QAAi[0][0] * v[462] + QAAi[1][0] * v[466] + QAAi[2][0] * v[471];
	v[452] = v[450] - v[443] * v[451];
	v[536] = QAAi[0][1] * v[447] + QAAi[1][1] * v[452] + QAAi[2][1] * v[458];
	v[533] = QAAi[0][0] * v[447] + QAAi[1][0] * v[452] + QAAi[2][0] * v[458];
	v[449] = -(v[445] * v[446]) / 2e0 + v[468];
	v[538] = QAAi[0][1] * v[449] + QAAi[1][1] * v[455] + QAAi[2][1] * v[460];
	v[535] = QAAi[0][0] * v[449] + QAAi[1][0] * v[455] + QAAi[2][0] * v[460];
	v[448] = (v[444] * v[446]) / 2e0 + v[485];
	v[537] = QAAi[0][1] * v[448] + QAAi[1][1] * v[454] + QAAi[2][1] * v[459];
	v[534] = QAAi[0][0] * v[448] + QAAi[1][0] * v[454] + QAAi[2][0] * v[459];
	v[285] = (v[146] * v[146]);
	v[4260] = 1e0 / Power(v[285], 3);
	v[3902] = QAAi[2][0] * v[285];
	v[3901] = QAAi[2][1] * v[285];
	v[3625] = v[296] / v[285];
	v[3624] = v[294] / v[285];
	v[3623] = v[288] / v[285];
	v[149] = 1e0 + (v[146] * v[446]) / 2e0;
	v[1717] = v[149] * v[3623];
	v[150] = v[146] * v[451];
	v[1719] = v[150] * v[3624];
	v[2315] = v[1717] + v[1719];
	v[151] = v[146] * v[457];
	v[4045] = v[151] / v[285];
	v[1720] = v[151] * v[3625];
	v[2320] = v[1717] + v[1720];
	v[2309] = v[1719] + v[2320];
	v[153] = v[146] * v[461];
	v[1724] = v[153] * v[3623];
	v[155] = 1e0 + (v[146] * v[465]) / 2e0;
	v[1712] = v[155] * v[3624];
	v[2311] = v[1712] + v[1724];
	v[156] = v[146] * v[470];
	v[4046] = v[156] / v[285];
	v[1713] = v[156] * v[3625];
	v[2321] = v[1712] + v[1713];
	v[2314] = v[1724] + v[2321];
	v[158] = v[146] * v[474];
	v[1728] = v[158] * v[3623];
	v[160] = v[146] * v[478];
	v[4047] = v[160] / v[285];
	v[1708] = v[160] * v[3624];
	v[161] = 1e0 + (v[146] * v[483]) / 2e0;
	v[1709] = v[161] * v[3625];
	v[2319] = v[1708] + v[1709] + v[1728];
	v[2316] = -v[1728] + v[2319];
	v[2310] = -v[1708] + v[2319];
	v[303] = -(v[450] * v[456]);
	v[3638] = v[303] - v[443];
	v[301] = v[453] * v[456];
	v[3636] = v[301] - v[444];
	v[292] = -(v[450] * v[453]);
	v[3634] = v[292] - v[445];
	v[165] = 4e0 / v[1412];
	v[4005] = -v[165] / 2e0;
	v[530] = -(v[165] * v[434]) / 2e0;
	v[531] = (v[489] * v[528]) / 2e0 + v[530];
	v[527] = -(v[165] * v[431]) / 2e0;
	v[529] = v[527] - (v[488] * v[528]) / 2e0;
	v[524] = v[165] - v[488] * v[523];
	v[521] = -v[165] + v[489] * v[519];
	v[516] = -v[165] - v[488] * v[515];
	v[513] = -(v[165] * v[436]) / 2e0;
	v[514] = -(v[490] * v[510]) / 2e0 + v[513];
	v[511] = -(v[488] * v[510]) / 2e0 + v[527];
	v[509] = v[165] - v[490] * v[506];
	v[504] = v[165] + v[489] * v[502];
	v[501] = -(v[165] * v[435]);
	v[525] = -v[501] + v[489] * v[523];
	v[582] = QBAi[0][1] * v[521] + QBAi[1][1] * v[525] + QBAi[2][1] * v[531];
	v[579] = QBAi[0][0] * v[521] + QBAi[1][0] * v[525] + QBAi[2][0] * v[531];
	v[520] = -v[501] - v[488] * v[519];
	v[581] = QBAi[0][1] * v[520] + QBAi[1][1] * v[524] + QBAi[2][1] * v[529];
	v[578] = QBAi[0][0] * v[520] + QBAi[1][0] * v[524] + QBAi[2][0] * v[529];
	v[517] = -v[501] + v[489] * v[515];
	v[503] = -v[501] - v[488] * v[502];
	v[500] = -v[165] - v[490] * v[496];
	v[498] = -(v[165] * v[433]);
	v[522] = -v[498] - v[490] * v[519];
	v[508] = -v[498] + v[489] * v[506];
	v[573] = QBAi[0][1] * v[508] + QBAi[1][1] * v[512] + QBAi[2][1] * v[517];
	v[570] = QBAi[0][0] * v[508] + QBAi[1][0] * v[512] + QBAi[2][0] * v[517];
	v[505] = -v[498] - v[490] * v[502];
	v[499] = v[489] * v[496] - v[498];
	v[495] = v[165] * v[432];
	v[526] = v[495] - v[490] * v[523];
	v[583] = QBAi[0][1] * v[522] + QBAi[1][1] * v[526] + QBAi[2][1] * v[532];
	v[580] = QBAi[0][0] * v[522] + QBAi[1][0] * v[526] + QBAi[2][0] * v[532];
	v[518] = v[495] - v[490] * v[515];
	v[574] = QBAi[0][1] * v[509] + QBAi[1][1] * v[514] + QBAi[2][1] * v[518];
	v[571] = QBAi[0][0] * v[509] + QBAi[1][0] * v[514] + QBAi[2][0] * v[518];
	v[507] = v[495] - v[488] * v[506];
	v[572] = QBAi[0][1] * v[507] + QBAi[1][1] * v[511] + QBAi[2][1] * v[516];
	v[569] = QBAi[0][0] * v[507] + QBAi[1][0] * v[511] + QBAi[2][0] * v[516];
	v[497] = v[495] - v[488] * v[496];
	v[563] = QBAi[0][1] * v[492] + QBAi[1][1] * v[497] + QBAi[2][1] * v[503];
	v[560] = QBAi[0][0] * v[492] + QBAi[1][0] * v[497] + QBAi[2][0] * v[503];
	v[494] = -(v[490] * v[491]) / 2e0 + v[513];
	v[565] = QBAi[0][1] * v[494] + QBAi[1][1] * v[500] + QBAi[2][1] * v[505];
	v[562] = QBAi[0][0] * v[494] + QBAi[1][0] * v[500] + QBAi[2][0] * v[505];
	v[493] = (v[489] * v[491]) / 2e0 + v[530];
	v[564] = QBAi[0][1] * v[493] + QBAi[1][1] * v[499] + QBAi[2][1] * v[504];
	v[561] = QBAi[0][0] * v[493] + QBAi[1][0] * v[499] + QBAi[2][0] * v[504];
	v[311] = (v[165] * v[165]);
	v[4254] = 1e0 / Power(v[311], 3);
	v[3896] = QBAi[2][0] * v[311];
	v[3895] = QBAi[2][1] * v[311];
	v[3628] = v[322] / v[311];
	v[3627] = v[320] / v[311];
	v[3626] = v[314] / v[311];
	v[168] = 1e0 + (v[165] * v[491]) / 2e0;
	v[1690] = v[168] * v[3626];
	v[169] = v[165] * v[496];
	v[1692] = v[169] * v[3627];
	v[2330] = v[1690] + v[1692];
	v[170] = v[165] * v[502];
	v[4042] = v[170] / v[311];
	v[1693] = v[170] * v[3628];
	v[2335] = v[1690] + v[1693];
	v[2324] = v[1692] + v[2335];
	v[172] = v[165] * v[506];
	v[1697] = v[172] * v[3626];
	v[174] = 1e0 + (v[165] * v[510]) / 2e0;
	v[1685] = v[174] * v[3627];
	v[2326] = v[1685] + v[1697];
	v[175] = v[165] * v[515];
	v[4043] = v[175] / v[311];
	v[1686] = v[175] * v[3628];
	v[2336] = v[1685] + v[1686];
	v[2329] = v[1697] + v[2336];
	v[177] = v[165] * v[519];
	v[1701] = v[177] * v[3626];
	v[179] = v[165] * v[523];
	v[4044] = v[179] / v[311];
	v[1681] = v[179] * v[3627];
	v[180] = 1e0 + (v[165] * v[528]) / 2e0;
	v[1682] = v[180] * v[3628];
	v[2334] = v[1681] + v[1682] + v[1701];
	v[2331] = -v[1701] + v[2334];
	v[2325] = -v[1681] + v[2334];
	v[329] = -(v[495] * v[501]);
	v[3643] = v[329] - v[488];
	v[327] = v[498] * v[501];
	v[3641] = v[327] - v[489];
	v[318] = -(v[495] * v[498]);
	v[3639] = v[318] - v[490];
	v[184] = QAAi[0][0] * v[149] + QAAi[1][0] * v[150] + QAAi[2][0] * v[151];
	v[185] = QAAi[0][1] * v[149] + QAAi[1][1] * v[150] + QAAi[2][1] * v[151];
	v[404] = v[1439] * v[185] - v[184] * v[993];
	v[187] = QAAi[0][0] * v[153] + QAAi[1][0] * v[155] + QAAi[2][0] * v[156];
	v[188] = QAAi[0][1] * v[153] + QAAi[1][1] * v[155] + QAAi[2][1] * v[156];
	v[402] = v[1439] * v[188] - v[187] * v[993];
	v[190] = QAAi[0][0] * v[158] + QAAi[1][0] * v[160] + QAAi[2][0] * v[161];
	v[191] = QAAi[0][1] * v[158] + QAAi[1][1] * v[160] + QAAi[2][1] * v[161];
	v[400] = v[1439] * v[191] - v[190] * v[993];
	v[193] = QBAi[0][0] * v[168] + QBAi[1][0] * v[169] + QBAi[2][0] * v[170];
	v[194] = QBAi[0][1] * v[168] + QBAi[1][1] * v[169] + QBAi[2][1] * v[170];
	v[403] = v[1439] * v[194] - v[193] * v[993];
	v[196] = QBAi[0][0] * v[172] + QBAi[1][0] * v[174] + QBAi[2][0] * v[175];
	v[197] = QBAi[0][1] * v[172] + QBAi[1][1] * v[174] + QBAi[2][1] * v[175];
	v[401] = v[1439] * v[197] - v[196] * v[993];
	v[199] = QBAi[0][0] * v[177] + QBAi[1][0] * v[179] + QBAi[2][0] * v[180];
	v[200] = QBAi[0][1] * v[177] + QBAi[1][1] * v[179] + QBAi[2][1] * v[180];
	v[399] = v[1439] * v[200] - v[199] * v[993];
	v[202] = dA[0] + xAAi[0];
	v[203] = dA[1] + xAAi[1];
	v[204] = dA[2] + xAAi[2];
	v[205] = dA[6] + xBAi[0];
	v[206] = dA[7] + xBAi[1];
	v[207] = dA[8] + xBAi[2];
	v[208] = (1e0 - ci[0]) / 2e0;
	v[209] = (1e0 + ci[0]) / 2e0;
	v[210] = (1e0 - cp[0]) / 2e0;
	v[211] = (1e0 + cp[0]) / 2e0;
	v[5635] = 0e0;
	v[5636] = 0e0;
	v[5637] = v[210];
	v[5638] = 0e0;
	v[5639] = 0e0;
	v[5640] = 0e0;
	v[5641] = 0e0;
	v[5642] = 0e0;
	v[5643] = v[211];
	v[5644] = 0e0;
	v[5645] = 0e0;
	v[5646] = 0e0;
	v[5647] = 0e0;
	v[5648] = 0e0;
	v[5649] = -1e0;
	v[5650] = 0e0;
	v[5651] = 0e0;
	v[5652] = 0e0;
	v[5653] = 0e0;
	v[5654] = v[210];
	v[5655] = 0e0;
	v[5656] = 0e0;
	v[5657] = 0e0;
	v[5658] = 0e0;
	v[5659] = 0e0;
	v[5660] = v[211];
	v[5661] = 0e0;
	v[5662] = 0e0;
	v[5663] = 0e0;
	v[5664] = 0e0;
	v[5665] = 0e0;
	v[5666] = -1e0;
	v[5667] = 0e0;
	v[5668] = 0e0;
	v[5669] = 0e0;
	v[5670] = 0e0;
	v[5671] = v[210];
	v[5672] = 0e0;
	v[5673] = 0e0;
	v[5674] = 0e0;
	v[5675] = 0e0;
	v[5676] = 0e0;
	v[5677] = v[211];
	v[5678] = 0e0;
	v[5679] = 0e0;
	v[5680] = 0e0;
	v[5681] = 0e0;
	v[5682] = 0e0;
	v[5683] = -1e0;
	v[5684] = 0e0;
	v[5685] = 0e0;
	v[5686] = 0e0;
	v[5687] = 0e0;
	v[5688] = 0e0;
	v[1276] = v[184] * v[210] + v[193] * v[211];
	v[1275] = v[185] * v[210] + v[194] * v[211];
	v[1273] = v[187] * v[210] + v[196] * v[211];
	v[1272] = v[188] * v[210] + v[197] * v[211];
	v[1269] = v[190] * v[210] + v[199] * v[211];
	v[1268] = v[191] * v[210] + v[200] * v[211];
	v[407] = v[211] * v[399] + v[210] * v[400];
	v[784] = v[211] * v[407];
	v[778] = v[210] * v[407];
	v[406] = v[211] * v[401] + v[210] * v[402];
	v[783] = v[211] * v[406];
	v[777] = v[210] * v[406];
	v[405] = v[211] * v[403] + v[210] * v[404];
	v[782] = v[211] * v[405];
	v[776] = v[210] * v[405];
	v[212] = cpointA[0] + (*radA)*cos(ci[1]);
	v[213] = cpointA[1] + (*radA)*sin(ci[1]);
	v[214] = cpointA[0] + v[1439];
	v[215] = cpointA[1] + v[993];
	v[1342] = QBAi[2][0] * v[214] + QBAi[2][1] * v[215];
	v[1341] = QBAi[1][0] * v[214] + QBAi[1][1] * v[215];
	v[1340] = QBAi[0][0] * v[214] + QBAi[0][1] * v[215];
	v[1339] = QAAi[2][0] * v[214] + QAAi[2][1] * v[215];
	v[1338] = QAAi[1][0] * v[214] + QAAi[1][1] * v[215];
	v[1337] = QAAi[0][0] * v[214] + QAAi[0][1] * v[215];
	v[634] = v[214] * v[535] + v[215] * v[538];
	v[643] = v[210] * v[634];
	v[633] = v[214] * v[534] + v[215] * v[537];
	v[642] = v[210] * v[633];
	v[632] = v[214] * v[533] + v[215] * v[536];
	v[641] = v[210] * v[632];
	v[631] = v[214] * v[562] + v[215] * v[565];
	v[646] = v[211] * v[631];
	v[630] = v[214] * v[561] + v[215] * v[564];
	v[645] = v[211] * v[630];
	v[629] = v[214] * v[560] + v[215] * v[563];
	v[644] = v[211] * v[629];
	v[622] = v[214] * v[544] + v[215] * v[547];
	v[649] = v[210] * v[622];
	v[621] = v[214] * v[543] + v[215] * v[546];
	v[648] = v[210] * v[621];
	v[620] = v[214] * v[542] + v[215] * v[545];
	v[647] = v[210] * v[620];
	v[619] = v[214] * v[571] + v[215] * v[574];
	v[652] = v[211] * v[619];
	v[618] = v[214] * v[570] + v[215] * v[573];
	v[651] = v[211] * v[618];
	v[617] = v[214] * v[569] + v[215] * v[572];
	v[650] = v[211] * v[617];
	v[610] = v[214] * v[553] + v[215] * v[556];
	v[655] = v[210] * v[610];
	v[609] = v[214] * v[552] + v[215] * v[555];
	v[654] = v[210] * v[609];
	v[608] = v[214] * v[551] + v[215] * v[554];
	v[653] = v[210] * v[608];
	v[607] = v[214] * v[580] + v[215] * v[583];
	v[658] = v[211] * v[607];
	v[606] = v[214] * v[579] + v[215] * v[582];
	v[657] = v[211] * v[606];
	v[605] = v[214] * v[578] + v[215] * v[581];
	v[656] = v[211] * v[605];
	v[394] = v[207] + v[199] * v[214] + v[200] * v[215];
	v[393] = v[204] + v[190] * v[214] + v[191] * v[215];
	v[3707] = -v[393] + v[394];
	v[395] = v[3707] / 2e0;
	v[391] = v[206] + v[196] * v[214] + v[197] * v[215];
	v[390] = v[203] + v[187] * v[214] + v[188] * v[215];
	v[3706] = -v[390] + v[391];
	v[392] = v[3706] / 2e0;
	v[388] = v[205] + v[193] * v[214] + v[194] * v[215];
	v[387] = v[202] + v[184] * v[214] + v[185] * v[215];
	v[3705] = -v[387] + v[388];
	v[389] = v[3705] / 2e0;
	v[225] = 4e0 / v[1422];
	v[4009] = -v[225] / 2e0;
	v[701] = -(v[225] * v[440]) / 2e0;
	v[702] = (v[660] * v[699]) / 2e0 + v[701];
	v[698] = -(v[225] * v[437]) / 2e0;
	v[700] = v[698] - (v[659] * v[699]) / 2e0;
	v[695] = v[225] - v[659] * v[694];
	v[692] = -v[225] + v[660] * v[690];
	v[687] = -v[225] - v[659] * v[686];
	v[684] = -(v[225] * v[442]) / 2e0;
	v[685] = -(v[661] * v[681]) / 2e0 + v[684];
	v[682] = -(v[659] * v[681]) / 2e0 + v[698];
	v[680] = v[225] - v[661] * v[677];
	v[675] = v[225] + v[660] * v[673];
	v[672] = -(v[225] * v[441]);
	v[696] = -v[672] + v[660] * v[694];
	v[741] = QABi[0][2] * v[692] + QABi[1][2] * v[696] + QABi[2][2] * v[702];
	v[738] = QABi[0][1] * v[692] + QABi[1][1] * v[696] + QABi[2][1] * v[702];
	v[735] = QABi[0][0] * v[692] + QABi[1][0] * v[696] + QABi[2][0] * v[702];
	v[691] = -v[672] - v[659] * v[690];
	v[740] = QABi[0][2] * v[691] + QABi[1][2] * v[695] + QABi[2][2] * v[700];
	v[737] = QABi[0][1] * v[691] + QABi[1][1] * v[695] + QABi[2][1] * v[700];
	v[734] = QABi[0][0] * v[691] + QABi[1][0] * v[695] + QABi[2][0] * v[700];
	v[688] = -v[672] + v[660] * v[686];
	v[674] = -v[672] - v[659] * v[673];
	v[671] = -v[225] - v[661] * v[667];
	v[669] = -(v[225] * v[439]);
	v[693] = -v[669] - v[661] * v[690];
	v[679] = -v[669] + v[660] * v[677];
	v[726] = QABi[0][2] * v[679] + QABi[1][2] * v[683] + QABi[2][2] * v[688];
	v[723] = QABi[0][1] * v[679] + QABi[1][1] * v[683] + QABi[2][1] * v[688];
	v[720] = QABi[0][0] * v[679] + QABi[1][0] * v[683] + QABi[2][0] * v[688];
	v[676] = -v[669] - v[661] * v[673];
	v[670] = v[660] * v[667] - v[669];
	v[666] = v[225] * v[438];
	v[697] = v[666] - v[661] * v[694];
	v[742] = QABi[0][2] * v[693] + QABi[1][2] * v[697] + QABi[2][2] * v[703];
	v[739] = QABi[0][1] * v[693] + QABi[1][1] * v[697] + QABi[2][1] * v[703];
	v[736] = QABi[0][0] * v[693] + QABi[1][0] * v[697] + QABi[2][0] * v[703];
	v[689] = v[666] - v[661] * v[686];
	v[727] = QABi[0][2] * v[680] + QABi[1][2] * v[685] + QABi[2][2] * v[689];
	v[724] = QABi[0][1] * v[680] + QABi[1][1] * v[685] + QABi[2][1] * v[689];
	v[721] = QABi[0][0] * v[680] + QABi[1][0] * v[685] + QABi[2][0] * v[689];
	v[678] = v[666] - v[659] * v[677];
	v[725] = QABi[0][2] * v[678] + QABi[1][2] * v[682] + QABi[2][2] * v[687];
	v[722] = QABi[0][1] * v[678] + QABi[1][1] * v[682] + QABi[2][1] * v[687];
	v[719] = QABi[0][0] * v[678] + QABi[1][0] * v[682] + QABi[2][0] * v[687];
	v[668] = v[666] - v[659] * v[667];
	v[710] = QABi[0][2] * v[663] + QABi[1][2] * v[668] + QABi[2][2] * v[674];
	v[707] = QABi[0][1] * v[663] + QABi[1][1] * v[668] + QABi[2][1] * v[674];
	v[704] = QABi[0][0] * v[663] + QABi[1][0] * v[668] + QABi[2][0] * v[674];
	v[665] = -(v[661] * v[662]) / 2e0 + v[684];
	v[712] = QABi[0][2] * v[665] + QABi[1][2] * v[671] + QABi[2][2] * v[676];
	v[709] = QABi[0][1] * v[665] + QABi[1][1] * v[671] + QABi[2][1] * v[676];
	v[706] = QABi[0][0] * v[665] + QABi[1][0] * v[671] + QABi[2][0] * v[676];
	v[664] = (v[660] * v[662]) / 2e0 + v[701];
	v[711] = QABi[0][2] * v[664] + QABi[1][2] * v[670] + QABi[2][2] * v[675];
	v[708] = QABi[0][1] * v[664] + QABi[1][1] * v[670] + QABi[2][1] * v[675];
	v[705] = QABi[0][0] * v[664] + QABi[1][0] * v[670] + QABi[2][0] * v[675];
	v[337] = (v[225] * v[225]);
	v[4248] = 1e0 / Power(v[337], 3);
	v[3631] = v[348] / v[337];
	v[3630] = v[346] / v[337];
	v[3629] = v[340] / v[337];
	v[228] = 1e0 + (v[225] * v[662]) / 2e0;
	v[1663] = v[228] * v[3629];
	v[229] = v[225] * v[667];
	v[1665] = v[229] * v[3630];
	v[2345] = v[1663] + v[1665];
	v[230] = v[225] * v[673];
	v[4039] = v[230] / v[337];
	v[1666] = v[230] * v[3631];
	v[2350] = v[1663] + v[1666];
	v[2339] = v[1665] + v[2350];
	v[232] = v[225] * v[677];
	v[1670] = v[232] * v[3629];
	v[234] = 1e0 + (v[225] * v[681]) / 2e0;
	v[1658] = v[234] * v[3630];
	v[2341] = v[1658] + v[1670];
	v[235] = v[225] * v[686];
	v[4040] = v[235] / v[337];
	v[1659] = v[235] * v[3631];
	v[2351] = v[1658] + v[1659];
	v[2344] = v[1670] + v[2351];
	v[237] = v[225] * v[690];
	v[1674] = v[237] * v[3629];
	v[239] = v[225] * v[694];
	v[4041] = v[239] / v[337];
	v[1654] = v[239] * v[3630];
	v[240] = 1e0 + (v[225] * v[699]) / 2e0;
	v[1655] = v[240] * v[3631];
	v[2349] = v[1654] + v[1655] + v[1674];
	v[2346] = -v[1674] + v[2349];
	v[2340] = -v[1654] + v[2349];
	v[355] = -(v[666] * v[672]);
	v[3648] = v[355] - v[659];
	v[353] = v[669] * v[672];
	v[3646] = v[353] - v[660];
	v[344] = -(v[666] * v[669]);
	v[3644] = v[344] - v[661];
	v[244] = QABi[0][0] * v[228] + QABi[1][0] * v[229] + QABi[2][0] * v[230];
	v[245] = QABi[0][1] * v[228] + QABi[1][1] * v[229] + QABi[2][1] * v[230];
	v[4025] = v[1440] * v[245];
	v[246] = QABi[0][2] * v[228] + QABi[1][2] * v[229] + QABi[2][2] * v[230];
	v[421] = v[4025] + v[244] * v[418] + v[246] * v[420];
	v[812] = -(v[211] * v[421]);
	v[806] = -(v[210] * v[421]);
	v[413] = v[244] * v[410] + v[246] * v[412];
	v[797] = -(v[211] * v[413]);
	v[791] = -(v[210] * v[413]);
	v[247] = QABi[0][0] * v[232] + QABi[1][0] * v[234] + QABi[2][0] * v[235];
	v[248] = QABi[0][1] * v[232] + QABi[1][1] * v[234] + QABi[2][1] * v[235];
	v[4026] = v[1440] * v[248];
	v[249] = QABi[0][2] * v[232] + QABi[1][2] * v[234] + QABi[2][2] * v[235];
	v[422] = v[4026] + v[247] * v[418] + v[249] * v[420];
	v[813] = -(v[211] * v[422]);
	v[807] = -(v[210] * v[422]);
	v[414] = v[247] * v[410] + v[249] * v[412];
	v[798] = -(v[211] * v[414]);
	v[792] = -(v[210] * v[414]);
	v[250] = QABi[0][0] * v[237] + QABi[1][0] * v[239] + QABi[2][0] * v[240];
	v[251] = QABi[0][1] * v[237] + QABi[1][1] * v[239] + QABi[2][1] * v[240];
	v[4028] = v[1440] * v[251];
	v[252] = QABi[0][2] * v[237] + QABi[1][2] * v[239] + QABi[2][2] * v[240];
	v[423] = v[4028] + v[250] * v[418] + v[252] * v[420];
	v[817] = -(v[421] * v[646]) - v[422] * v[652] - v[423] * v[658];
	v[816] = -(v[421] * v[645]) - v[422] * v[651] - v[423] * v[657];
	v[815] = -(v[421] * v[644]) - v[422] * v[650] - v[423] * v[656];
	v[814] = -(v[211] * v[423]);
	v[811] = -(v[421] * v[643]) - v[422] * v[649] - v[423] * v[655];
	v[810] = -(v[421] * v[642]) - v[422] * v[648] - v[423] * v[654];
	v[809] = -(v[421] * v[641]) - v[422] * v[647] - v[423] * v[653];
	v[808] = -(v[210] * v[423]);
	v[415] = v[250] * v[410] + v[252] * v[412];
	v[802] = -(v[413] * v[646]) - v[414] * v[652] - v[415] * v[658];
	v[801] = -(v[413] * v[645]) - v[414] * v[651] - v[415] * v[657];
	v[800] = -(v[413] * v[644]) - v[414] * v[650] - v[415] * v[656];
	v[799] = -(v[211] * v[415]);
	v[796] = -(v[413] * v[643]) - v[414] * v[649] - v[415] * v[655];
	v[795] = -(v[413] * v[642]) - v[414] * v[648] - v[415] * v[654];
	v[794] = -(v[413] * v[641]) - v[414] * v[647] - v[415] * v[653];
	v[793] = -(v[210] * v[415]);
	v[253] = dB[0] + xABi[0];
	v[254] = dB[1] + xABi[1];
	v[255] = dB[2] + xABi[2];
	v[256] = v[258] * (*xfac)*cos(ci[2]);
	v[257] = cpointB[1] + (*radB)*sin(ci[3]);
	v[259] = -(v[258] * (*zfac)*sin(ci[2]));
	v[260] = v[3632] * (*xfac);
	v[261] = cpointB[1] + v[994];
	v[263] = -(v[262] * v[3633]);
	v[757] = v[260] * v[736] + v[261] * v[739] + v[263] * v[742];
	v[756] = v[260] * v[735] + v[261] * v[738] + v[263] * v[741];
	v[755] = v[260] * v[734] + v[261] * v[737] + v[263] * v[740];
	v[754] = v[260] * v[721] + v[261] * v[724] + v[263] * v[727];
	v[753] = v[260] * v[720] + v[261] * v[723] + v[263] * v[726];
	v[752] = v[260] * v[719] + v[261] * v[722] + v[263] * v[725];
	v[751] = v[260] * v[706] + v[261] * v[709] + v[263] * v[712];
	v[790] = -(v[405] * v[751]) - v[406] * v[754] - v[407] * v[757];
	v[775] = -(v[389] * v[751]) - v[392] * v[754] - v[395] * v[757];
	v[750] = v[260] * v[705] + v[261] * v[708] + v[263] * v[711];
	v[789] = -(v[405] * v[750]) - v[406] * v[753] - v[407] * v[756];
	v[774] = -(v[389] * v[750]) - v[392] * v[753] - v[395] * v[756];
	v[749] = v[260] * v[704] + v[261] * v[707] + v[263] * v[710];
	v[788] = -(v[405] * v[749]) - v[406] * v[752] - v[407] * v[755];
	v[773] = -(v[389] * v[749]) - v[392] * v[752] - v[395] * v[755];
	v[273] = (*a4)*dA[0] + (*a6)*dduiA[0] + (*a5)*duiA[0];
	v[3653] = v[210] * v[273];
	v[274] = (*a4)*dA[1] + (*a6)*dduiA[1] + (*a5)*duiA[1];
	v[3656] = v[210] * v[274];
	v[275] = (*a4)*dA[2] + (*a6)*dduiA[2] + (*a5)*duiA[2];
	v[3677] = v[210] * v[275];
	v[276] = (*a4)*dA[6] + (*a6)*dduiA[6] + (*a5)*duiA[6];
	v[3652] = v[211] * v[276];
	v[277] = (*a4)*dA[7] + (*a6)*dduiA[7] + (*a5)*duiA[7];
	v[3655] = v[211] * v[277];
	v[278] = (*a4)*dA[8] + (*a6)*dduiA[8] + (*a5)*duiA[8];
	v[3676] = v[211] * v[278];
	v[279] = (*a4)*dB[0] + (*a6)*dduiB[0] + (*a5)*duiB[0];
	v[280] = (*a4)*dB[1] + (*a6)*dduiB[1] + (*a5)*duiB[1];
	v[281] = (*a4)*dB[2] + (*a6)*dduiB[2] + (*a5)*duiB[2];
	v[283] = v[301] + v[444];
	v[3940] = v[161] * v[283];
	v[3637] = v[283] / v[285];
	v[284] = v[292] + v[445];
	v[3936] = v[155] * v[284];
	v[3635] = v[284] / v[285];
	v[286] = v[285] + (v[453] * v[453]);
	v[3037] = v[286] / v[285];
	v[1861] = v[3636] / v[285];
	v[1859] = v[3634] / v[285];
	v[6597] = 0e0;
	v[6598] = 0e0;
	v[6599] = 0e0;
	v[6600] = 0e0;
	v[6601] = v[1859];
	v[6602] = v[1861];
	v[6603] = 0e0;
	v[6604] = 0e0;
	v[6605] = 0e0;
	v[6606] = 0e0;
	v[6607] = 0e0;
	v[6608] = 0e0;
	v[6609] = 0e0;
	v[6610] = 0e0;
	v[6611] = 0e0;
	v[6612] = 0e0;
	v[6613] = 0e0;
	v[6614] = 0e0;
	v[287] = v[1859] * v[294] + v[1861] * v[296] + v[288] * v[3037];
	v[289] = v[285] + (v[450] * v[450]);
	v[3939] = v[156] * v[289];
	v[3938] = v[153] * v[289];
	v[3934] = v[149] * v[3634];
	v[295] = v[303] + v[443];
	v[3941] = v[161] * v[295];
	v[3038] = v[289] / v[285];
	v[1860] = v[3638] / v[285];
	v[6615] = 0e0;
	v[6616] = 0e0;
	v[6617] = 0e0;
	v[6618] = v[3635];
	v[6619] = 0e0;
	v[6620] = v[1860];
	v[6621] = 0e0;
	v[6622] = 0e0;
	v[6623] = 0e0;
	v[6624] = 0e0;
	v[6625] = 0e0;
	v[6626] = 0e0;
	v[6627] = 0e0;
	v[6628] = 0e0;
	v[6629] = 0e0;
	v[6630] = 0e0;
	v[6631] = 0e0;
	v[6632] = 0e0;
	v[297] = v[1860] * v[296] + v[294] * v[3038] + v[288] * v[3635];
	v[298] = v[285] + (v[456] * v[456]);
	v[3933] = v[160] * v[298];
	v[3932] = v[158] * v[298];
	v[3935] = v[149] * v[3636];
	v[3937] = v[155] * v[3638];
	v[3039] = v[298] / v[285];
	v[1858] = v[295] / v[285];
	v[6633] = 0e0;
	v[6634] = 0e0;
	v[6635] = 0e0;
	v[6636] = v[3637];
	v[6637] = v[1858];
	v[6638] = 0e0;
	v[6639] = 0e0;
	v[6640] = 0e0;
	v[6641] = 0e0;
	v[6642] = 0e0;
	v[6643] = 0e0;
	v[6644] = 0e0;
	v[6645] = 0e0;
	v[6646] = 0e0;
	v[6647] = 0e0;
	v[6648] = 0e0;
	v[6649] = 0e0;
	v[6650] = 0e0;
	v[307] = v[1858] * v[294] + v[296] * v[3039] + v[288] * v[3637];
	v[309] = v[327] + v[489];
	v[3927] = v[180] * v[309];
	v[3642] = v[309] / v[311];
	v[310] = v[318] + v[490];
	v[3923] = v[174] * v[310];
	v[3640] = v[310] / v[311];
	v[312] = v[311] + (v[498] * v[498]);
	v[3040] = v[312] / v[311];
	v[1867] = v[3641] / v[311];
	v[1865] = v[3639] / v[311];
	v[6651] = 0e0;
	v[6652] = 0e0;
	v[6653] = 0e0;
	v[6654] = 0e0;
	v[6655] = 0e0;
	v[6656] = 0e0;
	v[6657] = 0e0;
	v[6658] = 0e0;
	v[6659] = 0e0;
	v[6660] = 0e0;
	v[6661] = v[1865];
	v[6662] = v[1867];
	v[6663] = 0e0;
	v[6664] = 0e0;
	v[6665] = 0e0;
	v[6666] = 0e0;
	v[6667] = 0e0;
	v[6668] = 0e0;
	v[313] = v[3040] * v[314] + v[1865] * v[320] + v[1867] * v[322];
	v[315] = v[311] + (v[495] * v[495]);
	v[3926] = v[175] * v[315];
	v[3925] = v[172] * v[315];
	v[3921] = v[168] * v[3639];
	v[321] = v[329] + v[488];
	v[3928] = v[180] * v[321];
	v[3041] = v[315] / v[311];
	v[1866] = v[3643] / v[311];
	v[6669] = 0e0;
	v[6670] = 0e0;
	v[6671] = 0e0;
	v[6672] = 0e0;
	v[6673] = 0e0;
	v[6674] = 0e0;
	v[6675] = 0e0;
	v[6676] = 0e0;
	v[6677] = 0e0;
	v[6678] = v[3640];
	v[6679] = 0e0;
	v[6680] = v[1866];
	v[6681] = 0e0;
	v[6682] = 0e0;
	v[6683] = 0e0;
	v[6684] = 0e0;
	v[6685] = 0e0;
	v[6686] = 0e0;
	v[323] = v[3041] * v[320] + v[1866] * v[322] + v[314] * v[3640];
	v[324] = v[311] + (v[501] * v[501]);
	v[3920] = v[179] * v[324];
	v[3919] = v[177] * v[324];
	v[3922] = v[168] * v[3641];
	v[3924] = v[174] * v[3643];
	v[3042] = v[324] / v[311];
	v[1864] = v[321] / v[311];
	v[6687] = 0e0;
	v[6688] = 0e0;
	v[6689] = 0e0;
	v[6690] = 0e0;
	v[6691] = 0e0;
	v[6692] = 0e0;
	v[6693] = 0e0;
	v[6694] = 0e0;
	v[6695] = 0e0;
	v[6696] = v[3642];
	v[6697] = v[1864];
	v[6698] = 0e0;
	v[6699] = 0e0;
	v[6700] = 0e0;
	v[6701] = 0e0;
	v[6702] = 0e0;
	v[6703] = 0e0;
	v[6704] = 0e0;
	v[333] = v[1864] * v[320] + v[3042] * v[322] + v[314] * v[3642];
	v[335] = v[353] + v[660];
	v[3913] = v[240] * v[335];
	v[3647] = v[335] / v[337];
	v[336] = v[344] + v[661];
	v[3909] = v[234] * v[336];
	v[3645] = v[336] / v[337];
	v[338] = v[337] + (v[669] * v[669]);
	v[3043] = v[338] / v[337];
	v[1873] = v[3646] / v[337];
	v[1871] = v[3644] / v[337];
	v[6705] = 0e0;
	v[6706] = 0e0;
	v[6707] = 0e0;
	v[6708] = 0e0;
	v[6709] = 0e0;
	v[6710] = 0e0;
	v[6711] = 0e0;
	v[6712] = 0e0;
	v[6713] = 0e0;
	v[6714] = 0e0;
	v[6715] = 0e0;
	v[6716] = 0e0;
	v[6717] = 0e0;
	v[6718] = 0e0;
	v[6719] = 0e0;
	v[6720] = 0e0;
	v[6721] = v[1871];
	v[6722] = v[1873];
	v[339] = v[3043] * v[340] + v[1871] * v[346] + v[1873] * v[348];
	v[341] = v[337] + (v[666] * v[666]);
	v[3912] = v[235] * v[341];
	v[3911] = v[232] * v[341];
	v[3907] = v[228] * v[3644];
	v[347] = v[355] + v[659];
	v[3914] = v[240] * v[347];
	v[3044] = v[341] / v[337];
	v[1872] = v[3648] / v[337];
	v[6723] = 0e0;
	v[6724] = 0e0;
	v[6725] = 0e0;
	v[6726] = 0e0;
	v[6727] = 0e0;
	v[6728] = 0e0;
	v[6729] = 0e0;
	v[6730] = 0e0;
	v[6731] = 0e0;
	v[6732] = 0e0;
	v[6733] = 0e0;
	v[6734] = 0e0;
	v[6735] = 0e0;
	v[6736] = 0e0;
	v[6737] = 0e0;
	v[6738] = v[3645];
	v[6739] = 0e0;
	v[6740] = v[1872];
	v[349] = v[3044] * v[346] + v[1872] * v[348] + v[340] * v[3645];
	v[350] = v[337] + (v[672] * v[672]);
	v[3906] = v[239] * v[350];
	v[3905] = v[237] * v[350];
	v[3908] = v[228] * v[3646];
	v[3910] = v[234] * v[3648];
	v[3045] = v[350] / v[337];
	v[1870] = v[347] / v[337];
	v[6741] = 0e0;
	v[6742] = 0e0;
	v[6743] = 0e0;
	v[6744] = 0e0;
	v[6745] = 0e0;
	v[6746] = 0e0;
	v[6747] = 0e0;
	v[6748] = 0e0;
	v[6749] = 0e0;
	v[6750] = 0e0;
	v[6751] = 0e0;
	v[6752] = 0e0;
	v[6753] = 0e0;
	v[6754] = 0e0;
	v[6755] = 0e0;
	v[6756] = v[3647];
	v[6757] = v[1870];
	v[6758] = 0e0;
	v[359] = v[1870] * v[346] + v[3045] * v[348] + v[340] * v[3647];
	v[3654] = -v[279] + v[287] * v[641] + v[297] * v[642] + v[307] * v[643] + v[313] * v[644] + v[323] * v[645] + v[333] * v[646]
		- v[339] * v[749] - v[349] * v[750] - v[359] * v[751];
	v[3650] = -v[281] + v[287] * v[653] + v[297] * v[654] + v[307] * v[655] + v[313] * v[656] + v[323] * v[657] + v[333] * v[658]
		- v[339] * v[755] - v[349] * v[756] - v[359] * v[757];
	v[3649] = -v[280] + v[287] * v[647] + v[297] * v[648] + v[307] * v[649] + v[313] * v[650] + v[323] * v[651] + v[333] * v[652]
		- v[339] * v[752] - v[349] * v[753] - v[359] * v[754];
	v[3294] = v[3652] + v[3653] + v[3654];
	v[3290] = v[3649] + v[3655] + v[3656];
	v[3287] = v[3650] + v[3676] + v[3677];
	v[3852] = v[3287] - v[3650];
	v[360] = -v[253] - v[244] * v[260] - v[245] * v[261] - v[246] * v[263] + v[210] * v[387] + v[211] * v[388];
	v[3721] = 2e0*v[360];
	v[764] = -v[360] / 2e0;
	v[765] = v[211] * v[389] - v[764];
	v[758] = v[210] * v[389] + v[764];
	v[361] = -v[254] - v[247] * v[260] - v[248] * v[261] - v[249] * v[263] + v[210] * v[390] + v[211] * v[391];
	v[3722] = 2e0*v[361];
	v[766] = -v[361] / 2e0;
	v[767] = v[211] * v[392] - v[766];
	v[759] = v[210] * v[392] + v[766];
	v[362] = -v[255] - v[250] * v[260] - v[251] * v[261] - v[252] * v[263] + v[210] * v[393] + v[211] * v[394];
	v[3723] = 2e0*v[362];
	v[4081] = v[1336] * (v[193] * v[3721] + v[196] * v[3722] + v[199] * v[3723]) + v[1335] * (v[194] * v[3721]
		+ v[197] * v[3722] + v[200] * v[3723]);
	v[4079] = v[1336] * (v[184] * v[3721] + v[187] * v[3722] + v[190] * v[3723]) + v[1335] * (v[185] * v[3721]
		+ v[188] * v[3722] + v[191] * v[3723]);
	v[1391] = 2e0*(v[1275] * v[360] + v[1272] * v[361] + v[1268] * v[362]);
	v[1389] = 2e0*(v[1276] * v[360] + v[1273] * v[361] + v[1269] * v[362]);
	v[1383] = -((*epsn)*(v[245] * v[360] + v[248] * v[361] + v[251] * v[362]));
	v[1644] = sqrt((v[360] * v[360]) + (v[361] * v[361]) + (v[362] * v[362]));
	v[2634] = 1e0 / Power(v[1644], 2);
	v[949] = -(QABi[2][0] * v[260]) - QABi[2][1] * v[261] - QABi[2][2] * v[263];
	v[950] = -(QABi[1][0] * v[260]) - QABi[1][1] * v[261] - QABi[1][2] * v[263];
	v[951] = -(QABi[0][0] * v[260]) - QABi[0][1] * v[261] - QABi[0][2] * v[263];
	v[952] = v[1342] * v[211];
	v[953] = v[1341] * v[211];
	v[954] = v[1340] * v[211];
	v[955] = v[1339] * v[210];
	v[956] = v[1338] * v[210];
	v[957] = v[1337] * v[210];
	v[820] = -(v[360] * (v[418] * v[706] + v[1440] * v[709] + v[420] * v[712])) - v[361] * (v[418] * v[721] + v[1440] * v[724]
		+ v[420] * v[727]) - v[362] * (v[418] * v[736] + v[1440] * v[739] + v[420] * v[742]) + v[421] * v[751] + v[422] * v[754]
		+ v[423] * v[757];
	v[819] = -(v[360] * (v[418] * v[705] + v[1440] * v[708] + v[420] * v[711])) - v[361] * (v[418] * v[720] + v[1440] * v[723]
		+ v[420] * v[726]) - v[362] * (v[418] * v[735] + v[1440] * v[738] + v[420] * v[741]) + v[421] * v[750] + v[422] * v[753]
		+ v[423] * v[756];
	v[818] = -(v[360] * (v[418] * v[704] + v[1440] * v[707] + v[420] * v[710])) - v[361] * (v[418] * v[719] + v[1440] * v[722]
		+ v[420] * v[725]) - v[362] * (v[418] * v[734] + v[1440] * v[737] + v[420] * v[740]) + v[421] * v[749] + v[422] * v[752]
		+ v[423] * v[755];
	v[805] = -(v[360] * (v[410] * v[706] + v[412] * v[712])) - v[361] * (v[410] * v[721] + v[412] * v[727]) - v[362] *
		(v[410] * v[736] + v[412] * v[742]) + v[413] * v[751] + v[414] * v[754] + v[415] * v[757];
	v[804] = -(v[360] * (v[410] * v[705] + v[412] * v[711])) - v[361] * (v[410] * v[720] + v[412] * v[726]) - v[362] *
		(v[410] * v[735] + v[412] * v[741]) + v[413] * v[750] + v[414] * v[753] + v[415] * v[756];
	v[803] = -(v[360] * (v[410] * v[704] + v[412] * v[710])) - v[361] * (v[410] * v[719] + v[412] * v[725]) - v[362] *
		(v[410] * v[734] + v[412] * v[740]) + v[413] * v[749] + v[414] * v[752] + v[415] * v[755];
	v[787] = v[405] * v[646] + v[406] * v[652] + v[407] * v[658] + v[211] * (v[360] * (v[1439] * v[565] - v[562] * v[993])
		+ v[361] * (v[1439] * v[574] - v[571] * v[993]) + v[362] * (v[1439] * v[583] - v[580] * v[993]));
	v[786] = v[405] * v[645] + v[406] * v[651] + v[407] * v[657] + v[211] * (v[360] * (v[1439] * v[564] - v[561] * v[993])
		+ v[361] * (v[1439] * v[573] - v[570] * v[993]) + v[362] * (v[1439] * v[582] - v[579] * v[993]));
	v[785] = v[405] * v[644] + v[406] * v[650] + v[407] * v[656] + v[211] * (v[360] * (v[1439] * v[563] - v[560] * v[993])
		+ v[361] * (v[1439] * v[572] - v[569] * v[993]) + v[362] * (v[1439] * v[581] - v[578] * v[993]));
	v[781] = v[405] * v[643] + v[406] * v[649] + v[407] * v[655] + v[210] * (v[360] * (v[1439] * v[538] - v[535] * v[993])
		+ v[361] * (v[1439] * v[547] - v[544] * v[993]) + v[362] * (v[1439] * v[556] - v[553] * v[993]));
	v[780] = v[405] * v[642] + v[406] * v[648] + v[407] * v[654] + v[210] * (v[360] * (v[1439] * v[537] - v[534] * v[993])
		+ v[361] * (v[1439] * v[546] - v[543] * v[993]) + v[362] * (v[1439] * v[555] - v[552] * v[993]));
	v[779] = v[405] * v[641] + v[406] * v[647] + v[407] * v[653] + v[210] * (v[360] * (v[1439] * v[536] - v[533] * v[993])
		+ v[361] * (v[1439] * v[545] - v[542] * v[993]) + v[362] * (v[1439] * v[554] - v[551] * v[993]));
	v[772] = (v[362] * v[607] + v[361] * v[619] + v[360] * v[631] + 2e0*v[389] * v[646] + 2e0*v[392] * v[652]
		+ 2e0*v[395] * v[658]) / 2e0;
	v[771] = (v[362] * v[606] + v[361] * v[618] + v[360] * v[630] + 2e0*v[389] * v[645] + 2e0*v[392] * v[651]
		+ 2e0*v[395] * v[657]) / 2e0;
	v[770] = (v[362] * v[605] + v[361] * v[617] + v[360] * v[629] + 2e0*v[389] * v[644] + 2e0*v[392] * v[650]
		+ 2e0*v[395] * v[656]) / 2e0;
	v[768] = -v[362] / 2e0;
	v[769] = v[211] * v[395] - v[768];
	v[763] = -(v[362] * v[610]) / 2e0 - (v[361] * v[622]) / 2e0 - (v[360] * v[634]) / 2e0 + v[389] * v[643] + v[392] * v[649]
		+ v[395] * v[655];
	v[762] = -(v[362] * v[609]) / 2e0 - (v[361] * v[621]) / 2e0 - (v[360] * v[633]) / 2e0 + v[389] * v[642] + v[392] * v[648]
		+ v[395] * v[654];
	v[761] = -(v[362] * v[608]) / 2e0 - (v[361] * v[620]) / 2e0 - (v[360] * v[632]) / 2e0 + v[389] * v[641] + v[392] * v[647]
		+ v[395] * v[653];
	v[760] = v[210] * v[395] + v[768];
	v[363] = -(QABi[0][0] * v[256]) - QABi[0][1] * v[257] - QABi[0][2] * v[259] + v[208] * (QAAi[0][0] * v[212]
		+ QAAi[0][1] * v[213] + xAAi[0]) - xABi[0] + v[209] * (QBAi[0][0] * v[212] + QBAi[0][1] * v[213] + xBAi[0]);
	v[364] = -(QABi[1][0] * v[256]) - QABi[1][1] * v[257] - QABi[1][2] * v[259] + v[208] * (QAAi[1][0] * v[212]
		+ QAAi[1][1] * v[213] + xAAi[1]) - xABi[1] + v[209] * (QBAi[1][0] * v[212] + QBAi[1][1] * v[213] + xBAi[1]);
	v[365] = -(QABi[2][0] * v[256]) - QABi[2][1] * v[257] - QABi[2][2] * v[259] + v[208] * (QAAi[2][0] * v[212]
		+ QAAi[2][1] * v[213] + xAAi[2]) - xABi[2] + v[209] * (QBAi[2][0] * v[212] + QBAi[2][1] * v[213] + xBAi[2]);
	if (v[1644] > 0.1e-7) { v01 = 1e0 / v[1644]; v02 = (-(v01 / v[1644])); v03 = (2e0*v01) / Power(v[1644], 2); }
	else {
		v01 = (12500000e0 / 3e0)*(24e0 - (-0.1e-7 + v[1644])*(0.24e10 - 2e0*(-1e0 + 100000000e0*v[1644])*
			(0.2399999997e10 - 0.1199999994e18*v[1644] - 0.3e17*(v[1644] * v[1644]))));
		v02 = (-50000000e0 / 3e0)*(0.3599999994e10 - 0.4799999982e18*v[1644] + 0.6e25*Power(v[1644], 3)
			+ 0.1799999982e26*(v[1644] * v[1644]));
		v03 = 0.1e17*(799999997e0 - 0.599999994e17*v[1644] - 0.3e17*(v[1644] * v[1644]));
	};
	v[372] = v03;
	v[373] = v02;
	v[374] = v01;
	v[375] = v[360] * v[374];
	v[3651] = -((*cn)*v[375]);
	v[3289] = v[3649] * v[375];
	v[3284] = v[3650] * v[375];
	v[2506] = v[359] * v[3651];
	v[2504] = v[349] * v[3651];
	v[2501] = v[339] * v[3651];
	v[2479] = -(v[375] * v[751]);
	v[2478] = -(v[375] * v[750]);
	v[2477] = -(v[375] * v[749]);
	v[2476] = v[375] * v[646];
	v[2475] = v[375] * v[645];
	v[2474] = v[375] * v[644];
	v[2473] = v[375] * v[643];
	v[2472] = v[375] * v[642];
	v[2471] = v[375] * v[641];
	v[2470] = -(v[279] * v[375]);
	v[982] = (v[375] * v[375]);
	v[376] = v[361] * v[374];
	v[3657] = v[375] * v[376];
	v[3296] = v[3654] * v[376];
	v[3285] = v[3650] * v[376];
	v[2499] = -(v[376] * v[754]);
	v[3675] = v[2479] + v[2499];
	v[2497] = -(v[376] * v[753]);
	v[3674] = v[2478] + v[2497];
	v[2495] = -(v[376] * v[752]);
	v[3673] = v[2477] + v[2495];
	v[2493] = v[376] * v[652];
	v[3672] = v[2476] + v[2493];
	v[2491] = v[376] * v[651];
	v[3671] = v[2475] + v[2491];
	v[2489] = v[376] * v[650];
	v[3670] = v[2474] + v[2489];
	v[2487] = v[376] * v[649];
	v[3669] = v[2473] + v[2487];
	v[2485] = v[376] * v[648];
	v[3668] = v[2472] + v[2485];
	v[2483] = v[376] * v[647];
	v[3667] = v[2471] + v[2483];
	v[2482] = -(v[280] * v[376]);
	v[3702] = v[2470] + v[2482] + v[3652] * v[375] + v[3653] * v[375] + v[3655] * v[376] + v[3656] * v[376];
	v[3286] = v[287] * v[3667] + v[297] * v[3668] + v[307] * v[3669] + v[313] * v[3670] + v[323] * v[3671] + v[333] * v[3672]
		+ v[339] * v[3673] + v[349] * v[3674] + v[359] * v[3675] + v[3702];
	v[986] = (v[376] * v[376]);
	v[985] = v[211] * v[3657];
	v[984] = v[210] * v[3657];
	v[377] = v[362] * v[374];
	v[3666] = -(v[377] * v[757]);
	v[3695] = v[2479] + v[3666];
	v[3686] = -v[2499] - v[3666];
	v[3665] = -(v[377] * v[756]);
	v[3694] = v[2478] + v[3665];
	v[3685] = -v[2497] - v[3665];
	v[3664] = -(v[377] * v[755]);
	v[3693] = v[2477] + v[3664];
	v[3684] = -v[2495] - v[3664];
	v[3663] = v[377] * v[658];
	v[3692] = v[2476] + v[3663];
	v[3683] = v[2493] + v[3663];
	v[3662] = v[377] * v[657];
	v[3691] = v[2475] + v[3662];
	v[3682] = v[2491] + v[3662];
	v[3661] = v[377] * v[656];
	v[3690] = v[2474] + v[3661];
	v[3681] = v[2489] + v[3661];
	v[3660] = v[377] * v[655];
	v[3689] = v[2473] + v[3660];
	v[3680] = v[2487] + v[3660];
	v[3659] = v[377] * v[654];
	v[3688] = v[2472] + v[3659];
	v[3679] = v[2485] + v[3659];
	v[3658] = v[377] * v[653];
	v[3687] = v[2471] + v[3658];
	v[3678] = v[2483] + v[3658];
	v[3297] = v[3294] * v[377];
	v[3292] = v[3290] * v[377];
	v[2631] = v[3687] * v[376] + v[647] * v[986];
	v[2630] = v[3678] * v[375] + v[641] * v[982];
	v[2627] = v[3688] * v[376] + v[648] * v[986];
	v[2626] = v[3679] * v[375] + v[642] * v[982];
	v[2623] = v[3689] * v[376] + v[649] * v[986];
	v[2622] = v[3680] * v[375] + v[643] * v[982];
	v[2619] = v[3690] * v[376] + v[650] * v[986];
	v[2618] = v[3681] * v[375] + v[644] * v[982];
	v[2615] = v[3691] * v[376] + v[651] * v[986];
	v[2614] = v[3682] * v[375] + v[645] * v[982];
	v[2611] = v[3692] * v[376] + v[652] * v[986];
	v[2610] = v[3683] * v[375] + v[646] * v[982];
	v[2607] = v[3693] * v[376] - v[752] * v[986];
	v[2606] = -(v[3684] * v[375]) - v[749] * v[982];
	v[2603] = v[3694] * v[376] - v[753] * v[986];
	v[2602] = -(v[3685] * v[375]) - v[750] * v[982];
	v[2599] = v[3695] * v[376] - v[754] * v[986];
	v[2598] = -(v[3686] * v[375]) - v[751] * v[982];
	v[990] = (v[377] * v[377]);
	v[3499] = (v[276] * v[375] + v[277] * v[376])*v[377] + v[278] * v[990];
	v[3498] = (v[273] * v[375] + v[274] * v[376])*v[377] + v[275] * v[990];
	v[4227] = -v[3498] + v[3499];
	v[2632] = v[3667] * v[377] + v[653] * v[990];
	v[2628] = v[3668] * v[377] + v[654] * v[990];
	v[2624] = v[3669] * v[377] + v[655] * v[990];
	v[2620] = v[3670] * v[377] + v[656] * v[990];
	v[2616] = v[3671] * v[377] + v[657] * v[990];
	v[2612] = v[3672] * v[377] + v[658] * v[990];
	v[2608] = v[3673] * v[377] - v[755] * v[990];
	v[2604] = v[3674] * v[377] - v[756] * v[990];
	v[2600] = v[3675] * v[377] - v[757] * v[990];
	v[2481] = (-v[281] + v[3676] + v[3677])*v[377];
	v[3701] = v[2470] + v[2481];
	v[3700] = v[2481] + v[2482];
	v[3295] = v[287] * v[3678] + v[297] * v[3679] + v[307] * v[3680] + v[313] * v[3681] + v[323] * v[3682] + v[333] * v[3683]
		- v[339] * v[3684] - v[349] * v[3685] - v[359] * v[3686] + v[3700];
	v[3291] = v[287] * v[3687] + v[297] * v[3688] + v[307] * v[3689] + v[313] * v[3690] + v[323] * v[3691] + v[333] * v[3692]
		+ v[339] * v[3693] + v[349] * v[3694] + v[359] * v[3695] + v[3701];
	v[378] = sqrt((v[363] * v[363]) + (v[364] * v[364]) + (v[365] * v[365]));
	if (v[378] > 0.1e-7) { v04 = 1e0 / v[378]; v05 = (-(v04 / v[378])); v06 = (2e0*v04) / Power(v[378], 2); }
	else {
		v04 = (12500000e0 / 3e0)*(24e0 - (-0.1e-7 + v[378])*(0.24e10 - 2e0*(-1e0 + 100000000e0*v[378])*(0.2399999997e10
			- 0.1199999994e18*v[378] - 0.3e17*(v[378] * v[378]))));
		v05 = (-50000000e0 / 3e0)*(0.3599999994e10 - 0.4799999982e18*v[378] + 0.6e25*Power(v[378], 3)
			+ 0.1799999982e26*(v[378] * v[378]));
		v06 = 0.1e17*(799999997e0 - 0.599999994e17*v[378] - 0.3e17*(v[378] * v[378]));
	};
	v[383] = v04;
	v[384] = v[363] * v[383];
	v[385] = v[364] * v[383];
	v[386] = v[365] * v[383];
	v[837] = -(invH[0][0] * v[758]) - invH[0][1] * v[776] - invH[0][2] * v[791] - invH[0][3] * v[806];
	v[838] = -(invH[0][0] * v[759]) - invH[0][1] * v[777] - invH[0][2] * v[792] - invH[0][3] * v[807];
	v[839] = -(invH[0][0] * v[760]) - invH[0][1] * v[778] - invH[0][2] * v[793] - invH[0][3] * v[808];
	v[840] = -(invH[0][0] * v[761]) - invH[0][1] * v[779] - invH[0][2] * v[794] - invH[0][3] * v[809];
	v[841] = -(invH[0][0] * v[762]) - invH[0][1] * v[780] - invH[0][2] * v[795] - invH[0][3] * v[810];
	v[842] = -(invH[0][0] * v[763]) - invH[0][1] * v[781] - invH[0][2] * v[796] - invH[0][3] * v[811];
	v[843] = -(invH[0][0] * v[765]) - invH[0][1] * v[782] - invH[0][2] * v[797] - invH[0][3] * v[812];
	v[844] = -(invH[0][0] * v[767]) - invH[0][1] * v[783] - invH[0][2] * v[798] - invH[0][3] * v[813];
	v[845] = -(invH[0][0] * v[769]) - invH[0][1] * v[784] - invH[0][2] * v[799] - invH[0][3] * v[814];
	v[846] = -(invH[0][0] * v[770]) - invH[0][1] * v[785] - invH[0][2] * v[800] - invH[0][3] * v[815];
	v[847] = -(invH[0][0] * v[771]) - invH[0][1] * v[786] - invH[0][2] * v[801] - invH[0][3] * v[816];
	v[848] = -(invH[0][0] * v[772]) - invH[0][1] * v[787] - invH[0][2] * v[802] - invH[0][3] * v[817];
	v[849] = invH[0][0] * v[389] + invH[0][1] * v[405] - invH[0][2] * v[413] - invH[0][3] * v[421];
	v[850] = invH[0][0] * v[392] + invH[0][1] * v[406] - invH[0][2] * v[414] - invH[0][3] * v[422];
	v[851] = invH[0][0] * v[395] + invH[0][1] * v[407] - invH[0][2] * v[415] - invH[0][3] * v[423];
	v[852] = -(invH[0][0] * v[773]) - invH[0][1] * v[788] - invH[0][2] * v[803] - invH[0][3] * v[818];
	v[853] = -(invH[0][0] * v[774]) - invH[0][1] * v[789] - invH[0][2] * v[804] - invH[0][3] * v[819];
	v[854] = -(invH[0][0] * v[775]) - invH[0][1] * v[790] - invH[0][2] * v[805] - invH[0][3] * v[820];
	v[1533] = -(v[273] * v[837]) - v[274] * v[838] - v[275] * v[839] - v[287] * v[840] - v[297] * v[841] - v[307] * v[842]
		- v[276] * v[843] - v[277] * v[844] - v[278] * v[845] - v[313] * v[846] - v[323] * v[847] - v[333] * v[848] - v[279] * v[849]
		- v[280] * v[850] - v[281] * v[851] - v[339] * v[852] - v[349] * v[853] - v[359] * v[854];
	v[4857] = v[837];
	v[4858] = v[838];
	v[4859] = v[839];
	v[4860] = v[840];
	v[4861] = v[841];
	v[4862] = v[842];
	v[4863] = v[843];
	v[4864] = v[844];
	v[4865] = v[845];
	v[4866] = v[846];
	v[4867] = v[847];
	v[4868] = v[848];
	v[4869] = v[849];
	v[4870] = v[850];
	v[4871] = v[851];
	v[4872] = v[852];
	v[4873] = v[853];
	v[4874] = v[854];
	v[855] = -(invH[1][0] * v[758]) - invH[1][1] * v[776] - invH[1][2] * v[791] - invH[1][3] * v[806];
	v[856] = -(invH[1][0] * v[759]) - invH[1][1] * v[777] - invH[1][2] * v[792] - invH[1][3] * v[807];
	v[857] = -(invH[1][0] * v[760]) - invH[1][1] * v[778] - invH[1][2] * v[793] - invH[1][3] * v[808];
	v[858] = -(invH[1][0] * v[761]) - invH[1][1] * v[779] - invH[1][2] * v[794] - invH[1][3] * v[809];
	v[859] = -(invH[1][0] * v[762]) - invH[1][1] * v[780] - invH[1][2] * v[795] - invH[1][3] * v[810];
	v[860] = -(invH[1][0] * v[763]) - invH[1][1] * v[781] - invH[1][2] * v[796] - invH[1][3] * v[811];
	v[861] = -(invH[1][0] * v[765]) - invH[1][1] * v[782] - invH[1][2] * v[797] - invH[1][3] * v[812];
	v[862] = -(invH[1][0] * v[767]) - invH[1][1] * v[783] - invH[1][2] * v[798] - invH[1][3] * v[813];
	v[863] = -(invH[1][0] * v[769]) - invH[1][1] * v[784] - invH[1][2] * v[799] - invH[1][3] * v[814];
	v[864] = -(invH[1][0] * v[770]) - invH[1][1] * v[785] - invH[1][2] * v[800] - invH[1][3] * v[815];
	v[865] = -(invH[1][0] * v[771]) - invH[1][1] * v[786] - invH[1][2] * v[801] - invH[1][3] * v[816];
	v[866] = -(invH[1][0] * v[772]) - invH[1][1] * v[787] - invH[1][2] * v[802] - invH[1][3] * v[817];
	v[867] = invH[1][0] * v[389] + invH[1][1] * v[405] - invH[1][2] * v[413] - invH[1][3] * v[421];
	v[868] = invH[1][0] * v[392] + invH[1][1] * v[406] - invH[1][2] * v[414] - invH[1][3] * v[422];
	v[869] = invH[1][0] * v[395] + invH[1][1] * v[407] - invH[1][2] * v[415] - invH[1][3] * v[423];
	v[870] = -(invH[1][0] * v[773]) - invH[1][1] * v[788] - invH[1][2] * v[803] - invH[1][3] * v[818];
	v[871] = -(invH[1][0] * v[774]) - invH[1][1] * v[789] - invH[1][2] * v[804] - invH[1][3] * v[819];
	v[872] = -(invH[1][0] * v[775]) - invH[1][1] * v[790] - invH[1][2] * v[805] - invH[1][3] * v[820];
	v[1535] = -(v[273] * v[855]) - v[274] * v[856] - v[275] * v[857] - v[287] * v[858] - v[297] * v[859] - v[307] * v[860]
		- v[276] * v[861] - v[277] * v[862] - v[278] * v[863] - v[313] * v[864] - v[323] * v[865] - v[333] * v[866] - v[279] * v[867]
		- v[280] * v[868] - v[281] * v[869] - v[339] * v[870] - v[349] * v[871] - v[359] * v[872];
	v[4875] = v[855];
	v[4876] = v[856];
	v[4877] = v[857];
	v[4878] = v[858];
	v[4879] = v[859];
	v[4880] = v[860];
	v[4881] = v[861];
	v[4882] = v[862];
	v[4883] = v[863];
	v[4884] = v[864];
	v[4885] = v[865];
	v[4886] = v[866];
	v[4887] = v[867];
	v[4888] = v[868];
	v[4889] = v[869];
	v[4890] = v[870];
	v[4891] = v[871];
	v[4892] = v[872];
	v[873] = -(invH[2][0] * v[758]) - invH[2][1] * v[776] - invH[2][2] * v[791] - invH[2][3] * v[806];
	v[874] = -(invH[2][0] * v[759]) - invH[2][1] * v[777] - invH[2][2] * v[792] - invH[2][3] * v[807];
	v[875] = -(invH[2][0] * v[760]) - invH[2][1] * v[778] - invH[2][2] * v[793] - invH[2][3] * v[808];
	v[876] = -(invH[2][0] * v[761]) - invH[2][1] * v[779] - invH[2][2] * v[794] - invH[2][3] * v[809];
	v[877] = -(invH[2][0] * v[762]) - invH[2][1] * v[780] - invH[2][2] * v[795] - invH[2][3] * v[810];
	v[878] = -(invH[2][0] * v[763]) - invH[2][1] * v[781] - invH[2][2] * v[796] - invH[2][3] * v[811];
	v[879] = -(invH[2][0] * v[765]) - invH[2][1] * v[782] - invH[2][2] * v[797] - invH[2][3] * v[812];
	v[880] = -(invH[2][0] * v[767]) - invH[2][1] * v[783] - invH[2][2] * v[798] - invH[2][3] * v[813];
	v[881] = -(invH[2][0] * v[769]) - invH[2][1] * v[784] - invH[2][2] * v[799] - invH[2][3] * v[814];
	v[882] = -(invH[2][0] * v[770]) - invH[2][1] * v[785] - invH[2][2] * v[800] - invH[2][3] * v[815];
	v[883] = -(invH[2][0] * v[771]) - invH[2][1] * v[786] - invH[2][2] * v[801] - invH[2][3] * v[816];
	v[884] = -(invH[2][0] * v[772]) - invH[2][1] * v[787] - invH[2][2] * v[802] - invH[2][3] * v[817];
	v[885] = invH[2][0] * v[389] + invH[2][1] * v[405] - invH[2][2] * v[413] - invH[2][3] * v[421];
	v[886] = invH[2][0] * v[392] + invH[2][1] * v[406] - invH[2][2] * v[414] - invH[2][3] * v[422];
	v[887] = invH[2][0] * v[395] + invH[2][1] * v[407] - invH[2][2] * v[415] - invH[2][3] * v[423];
	v[888] = -(invH[2][0] * v[773]) - invH[2][1] * v[788] - invH[2][2] * v[803] - invH[2][3] * v[818];
	v[889] = -(invH[2][0] * v[774]) - invH[2][1] * v[789] - invH[2][2] * v[804] - invH[2][3] * v[819];
	v[890] = -(invH[2][0] * v[775]) - invH[2][1] * v[790] - invH[2][2] * v[805] - invH[2][3] * v[820];
	v[1529] = v[273] * v[873] + v[274] * v[874] + v[275] * v[875] + v[287] * v[876] + v[297] * v[877] + v[307] * v[878]
		+ v[276] * v[879] + v[277] * v[880] + v[278] * v[881] + v[313] * v[882] + v[323] * v[883] + v[333] * v[884] + v[279] * v[885]
		+ v[280] * v[886] + v[281] * v[887] + v[339] * v[888] + v[349] * v[889] + v[359] * v[890];
	v[4893] = v[873];
	v[4894] = v[874];
	v[4895] = v[875];
	v[4896] = v[876];
	v[4897] = v[877];
	v[4898] = v[878];
	v[4899] = v[879];
	v[4900] = v[880];
	v[4901] = v[881];
	v[4902] = v[882];
	v[4903] = v[883];
	v[4904] = v[884];
	v[4905] = v[885];
	v[4906] = v[886];
	v[4907] = v[887];
	v[4908] = v[888];
	v[4909] = v[889];
	v[4910] = v[890];
	v[891] = -(invH[3][0] * v[758]) - invH[3][1] * v[776] - invH[3][2] * v[791] - invH[3][3] * v[806];
	v[1443] = v[395] * v[837] + v[407] * v[855] - v[415] * v[873] - v[423] * v[891];
	v[1442] = v[392] * v[837] + v[406] * v[855] - v[414] * v[873] - v[422] * v[891];
	v[1441] = v[389] * v[837] + v[405] * v[855] - v[413] * v[873] - v[421] * v[891];
	v[892] = -(invH[3][0] * v[759]) - invH[3][1] * v[777] - invH[3][2] * v[792] - invH[3][3] * v[807];
	v[1447] = v[395] * v[838] + v[407] * v[856] - v[415] * v[874] - v[423] * v[892];
	v[1446] = v[392] * v[838] + v[406] * v[856] - v[414] * v[874] - v[422] * v[892];
	v[1445] = v[389] * v[838] + v[405] * v[856] - v[413] * v[874] - v[421] * v[892];
	v[893] = -(invH[3][0] * v[760]) - invH[3][1] * v[778] - invH[3][2] * v[793] - invH[3][3] * v[808];
	v[1451] = v[395] * v[839] + v[407] * v[857] - v[415] * v[875] - v[423] * v[893];
	v[1450] = v[392] * v[839] + v[406] * v[857] - v[414] * v[875] - v[422] * v[893];
	v[1449] = v[389] * v[839] + v[405] * v[857] - v[413] * v[875] - v[421] * v[893];
	v[894] = -(invH[3][0] * v[761]) - invH[3][1] * v[779] - invH[3][2] * v[794] - invH[3][3] * v[809];
	v[1455] = v[395] * v[840] + v[407] * v[858] - v[415] * v[876] - v[423] * v[894];
	v[1454] = v[392] * v[840] + v[406] * v[858] - v[414] * v[876] - v[422] * v[894];
	v[1453] = v[389] * v[840] + v[405] * v[858] - v[413] * v[876] - v[421] * v[894];
	v[895] = -(invH[3][0] * v[762]) - invH[3][1] * v[780] - invH[3][2] * v[795] - invH[3][3] * v[810];
	v[1459] = v[395] * v[841] + v[407] * v[859] - v[415] * v[877] - v[423] * v[895];
	v[1458] = v[392] * v[841] + v[406] * v[859] - v[414] * v[877] - v[422] * v[895];
	v[1457] = v[389] * v[841] + v[405] * v[859] - v[413] * v[877] - v[421] * v[895];
	v[896] = -(invH[3][0] * v[763]) - invH[3][1] * v[781] - invH[3][2] * v[796] - invH[3][3] * v[811];
	v[1463] = v[395] * v[842] + v[407] * v[860] - v[415] * v[878] - v[423] * v[896];
	v[1462] = v[392] * v[842] + v[406] * v[860] - v[414] * v[878] - v[422] * v[896];
	v[1461] = v[389] * v[842] + v[405] * v[860] - v[413] * v[878] - v[421] * v[896];
	v[897] = -(invH[3][0] * v[765]) - invH[3][1] * v[782] - invH[3][2] * v[797] - invH[3][3] * v[812];
	v[1467] = v[395] * v[843] + v[407] * v[861] - v[415] * v[879] - v[423] * v[897];
	v[1466] = v[392] * v[843] + v[406] * v[861] - v[414] * v[879] - v[422] * v[897];
	v[1465] = v[389] * v[843] + v[405] * v[861] - v[413] * v[879] - v[421] * v[897];
	v[898] = -(invH[3][0] * v[767]) - invH[3][1] * v[783] - invH[3][2] * v[798] - invH[3][3] * v[813];
	v[1471] = v[395] * v[844] + v[407] * v[862] - v[415] * v[880] - v[423] * v[898];
	v[1470] = v[392] * v[844] + v[406] * v[862] - v[414] * v[880] - v[422] * v[898];
	v[1469] = v[389] * v[844] + v[405] * v[862] - v[413] * v[880] - v[421] * v[898];
	v[899] = -(invH[3][0] * v[769]) - invH[3][1] * v[784] - invH[3][2] * v[799] - invH[3][3] * v[814];
	v[1475] = v[395] * v[845] + v[407] * v[863] - v[415] * v[881] - v[423] * v[899];
	v[1474] = v[392] * v[845] + v[406] * v[863] - v[414] * v[881] - v[422] * v[899];
	v[1473] = v[389] * v[845] + v[405] * v[863] - v[413] * v[881] - v[421] * v[899];
	v[900] = -(invH[3][0] * v[770]) - invH[3][1] * v[785] - invH[3][2] * v[800] - invH[3][3] * v[815];
	v[1479] = v[395] * v[846] + v[407] * v[864] - v[415] * v[882] - v[423] * v[900];
	v[1478] = v[392] * v[846] + v[406] * v[864] - v[414] * v[882] - v[422] * v[900];
	v[1477] = v[389] * v[846] + v[405] * v[864] - v[413] * v[882] - v[421] * v[900];
	v[901] = -(invH[3][0] * v[771]) - invH[3][1] * v[786] - invH[3][2] * v[801] - invH[3][3] * v[816];
	v[1483] = v[395] * v[847] + v[407] * v[865] - v[415] * v[883] - v[423] * v[901];
	v[1482] = v[392] * v[847] + v[406] * v[865] - v[414] * v[883] - v[422] * v[901];
	v[1481] = v[389] * v[847] + v[405] * v[865] - v[413] * v[883] - v[421] * v[901];
	v[902] = -(invH[3][0] * v[772]) - invH[3][1] * v[787] - invH[3][2] * v[802] - invH[3][3] * v[817];
	v[1487] = v[395] * v[848] + v[407] * v[866] - v[415] * v[884] - v[423] * v[902];
	v[1486] = v[392] * v[848] + v[406] * v[866] - v[414] * v[884] - v[422] * v[902];
	v[1485] = v[389] * v[848] + v[405] * v[866] - v[413] * v[884] - v[421] * v[902];
	v[903] = invH[3][0] * v[389] + invH[3][1] * v[405] - invH[3][2] * v[413] - invH[3][3] * v[421];
	v[1491] = v[395] * v[849] + v[407] * v[867] - v[415] * v[885] - v[423] * v[903];
	v[1490] = v[392] * v[849] + v[406] * v[867] - v[414] * v[885] - v[422] * v[903];
	v[1489] = v[389] * v[849] + v[405] * v[867] - v[413] * v[885] - v[421] * v[903];
	v[904] = invH[3][0] * v[392] + invH[3][1] * v[406] - invH[3][2] * v[414] - invH[3][3] * v[422];
	v[1495] = v[395] * v[850] + v[407] * v[868] - v[415] * v[886] - v[423] * v[904];
	v[1494] = v[392] * v[850] + v[406] * v[868] - v[414] * v[886] - v[422] * v[904];
	v[1493] = v[389] * v[850] + v[405] * v[868] - v[413] * v[886] - v[421] * v[904];
	v[905] = invH[3][0] * v[395] + invH[3][1] * v[407] - invH[3][2] * v[415] - invH[3][3] * v[423];
	v[1499] = v[395] * v[851] + v[407] * v[869] - v[415] * v[887] - v[423] * v[905];
	v[1498] = v[392] * v[851] + v[406] * v[869] - v[414] * v[887] - v[422] * v[905];
	v[1497] = v[389] * v[851] + v[405] * v[869] - v[413] * v[887] - v[421] * v[905];
	v[906] = -(invH[3][0] * v[773]) - invH[3][1] * v[788] - invH[3][2] * v[803] - invH[3][3] * v[818];
	v[1503] = v[395] * v[852] + v[407] * v[870] - v[415] * v[888] - v[423] * v[906];
	v[1502] = v[392] * v[852] + v[406] * v[870] - v[414] * v[888] - v[422] * v[906];
	v[1501] = v[389] * v[852] + v[405] * v[870] - v[413] * v[888] - v[421] * v[906];
	v[907] = -(invH[3][0] * v[774]) - invH[3][1] * v[789] - invH[3][2] * v[804] - invH[3][3] * v[819];
	v[1507] = v[395] * v[853] + v[407] * v[871] - v[415] * v[889] - v[423] * v[907];
	v[1506] = v[392] * v[853] + v[406] * v[871] - v[414] * v[889] - v[422] * v[907];
	v[1505] = v[389] * v[853] + v[405] * v[871] - v[413] * v[889] - v[421] * v[907];
	v[908] = -(invH[3][0] * v[775]) - invH[3][1] * v[790] - invH[3][2] * v[805] - invH[3][3] * v[820];
	v[1531] = v[273] * v[891] + v[274] * v[892] + v[275] * v[893] + v[287] * v[894] + v[297] * v[895] + v[307] * v[896]
		+ v[276] * v[897] + v[277] * v[898] + v[278] * v[899] + v[313] * v[900] + v[323] * v[901] + v[333] * v[902] + v[279] * v[903]
		+ v[280] * v[904] + v[281] * v[905] + v[339] * v[906] + v[349] * v[907] + v[359] * v[908];
	v[1511] = v[395] * v[854] + v[407] * v[872] - v[415] * v[890] - v[423] * v[908];
	v[1510] = v[392] * v[854] + v[406] * v[872] - v[414] * v[890] - v[422] * v[908];
	v[1509] = v[389] * v[854] + v[405] * v[872] - v[413] * v[890] - v[421] * v[908];
	v[4839] = v[891];
	v[4840] = v[892];
	v[4841] = v[893];
	v[4842] = v[894];
	v[4843] = v[895];
	v[4844] = v[896];
	v[4845] = v[897];
	v[4846] = v[898];
	v[4847] = v[899];
	v[4848] = v[900];
	v[4849] = v[901];
	v[4850] = v[902];
	v[4851] = v[903];
	v[4852] = v[904];
	v[4853] = v[905];
	v[4854] = v[906];
	v[4855] = v[907];
	v[4856] = v[908];
	b909 = sqrt(Power(v[376] * v[384] - v[375] * v[385], 2) + Power(-(v[377] * v[384]) + v[375] * v[386], 2) + Power
	(v[377] * v[385] - v[376] * v[386], 2)) > 0.1e-7;
	if (b909) {
		v[911] = v[377] * v[385] - v[376] * v[386];
		v[912] = -(v[377] * v[384]) + v[375] * v[386];
		v[913] = v[376] * v[384] - v[375] * v[385];
		v[914] = sqrt((v[911] * v[911]) + (v[912] * v[912]) + (v[913] * v[913]));
		v[2245] = 1e0 / Power(v[914], 2);
		v[1637] = v[914];
		v[2256] = 1e0 - (v[1637] * v[1637]);
		v[4118] = 1e0 / Power(v[2256], 0.15e1);
		v[2251] = 1e0 / sqrt(v[2256]);
		v[1636] = asin(v[1637]) / 2e0;
		v[4117] = tan(v[1636]);
		v[2250] = 1e0 / Power(cos(v[1636]), 2);
		v[3808] = v[2250] * v[2251];
		v[916] = 2e0*tan(v[1636]);
		if (v[914] > 0.1e-7) { v07 = 1e0 / v[914]; v08 = (-(v07 / v[914])); v09 = (2e0*v07) / Power(v[914], 2); }
		else {
			v07 = (12500000e0 / 3e0)*(24e0 - (-0.1e-7 + v[914])*(0.24e10 - 2e0*(-1e0 + 100000000e0*v[914])*
				(0.2399999997e10 - 0.1199999994e18*v[914] - 0.3e17*(v[914] * v[914]))));
			v08 = (-50000000e0 / 3e0)*(0.3599999994e10 - 0.4799999982e18*v[914] + 0.6e25*Power(v[914], 3)
				+ 0.1799999982e26*(v[914] * v[914]));
			v09 = 0.1e17*(799999997e0 - 0.599999994e17*v[914] - 0.3e17*(v[914] * v[914]));
		};
		v[920] = v09;
		v[921] = v08;
		v[922] = v07;
		v[4116] = v[916] * v[921] + v[3808] * v[922];
		v[3696] = v[916] * v[922];
		v[923] = v[3696] * v[911];
		v[934] = (v[923] * v[923]);
		v[924] = v[3696] * v[912];
		v[932] = (v[923] * v[924]) / 2e0;
		v[927] = (v[924] * v[924]);
		v[1620] = -v[927] - v[934];
		v[925] = v[3696] * v[913];
		v[1615] = v[925] + v[932];
		v[1613] = -v[925] + v[932];
		v[939] = (v[924] * v[925]) / 2e0;
		v[1619] = v[923] + v[939];
		v[1617] = -v[923] + v[939];
		v[937] = (v[923] * v[925]) / 2e0;
		v[1618] = -v[924] + v[937];
		v[1614] = v[924] + v[937];
		v[928] = (v[925] * v[925]);
		v[1624] = 4e0 + v[927] + v[928] + v[934];
		v[4119] = 1e0 / Power(v[1624], 3);
		v[2579] = 1e0 / Power(v[1624], 2);
		v[1616] = -v[928] - v[934];
		v[1612] = -v[927] - v[928];
		v[926] = 4e0 / v[1624];
		v[929] = 1e0 + (v[1612] * v[926]) / 2e0;
		v[930] = v[1613] * v[926];
		v[931] = v[1614] * v[926];
		v[933] = v[1615] * v[926];
		v[935] = 1e0 + (v[1616] * v[926]) / 2e0;
		v[936] = v[1617] * v[926];
		v[938] = v[1618] * v[926];
		v[940] = v[1619] * v[926];
		v[941] = 1e0 + (v[1620] * v[926]) / 2e0;
	}
	else {
		v[929] = 1e0;
		v[930] = 0e0;
		v[931] = 0e0;
		v[933] = 0e0;
		v[935] = 1e0;
		v[936] = 0e0;
		v[938] = 0e0;
		v[940] = 0e0;
		v[941] = 1e0;
	};
	if ((*previouscontact)) {
		v[1572] = 1e0 - v[990];
		v[1570] = 1e0 - v[986];
		v[1568] = 1e0 - v[982];
		v[946] = v[208] * (v[204] + v[190] * v[212] + v[191] * v[213]) + v[209] * (v[207] + v[199] * v[212] + v[200] * v[213])
			- v[255] - v[250] * v[256] - v[251] * v[257] - v[252] * v[259] + gti[0] * v[938] + gti[1] * v[940] + gti[2] * v[941];
		v[3697] = v[377] * v[946];
		v[945] = v[208] * (v[203] + v[187] * v[212] + v[188] * v[213]) + v[209] * (v[206] + v[196] * v[212] + v[197] * v[213])
			- v[254] - v[247] * v[256] - v[248] * v[257] - v[249] * v[259] + gti[0] * v[933] + gti[1] * v[935] + gti[2] * v[936];
		v[3699] = v[376] * v[945];
		v[3740] = v[3697] + v[3699];
		v[944] = v[208] * (v[202] + v[184] * v[212] + v[185] * v[213]) + v[209] * (v[205] + v[193] * v[212] + v[194] * v[213])
			- v[253] - v[244] * v[256] - v[245] * v[257] - v[246] * v[259] + gti[0] * v[929] + gti[1] * v[930] + gti[2] * v[931];
		v[3698] = -(v[375] * v[944]);
		v[3739] = -v[3697] + v[3698];
		v[3738] = v[3698] - v[3699];
		v[943] = -(v[3740] * v[375]) + v[1568] * v[944];
		v[947] = v[3739] * v[376] + v[1570] * v[945];
		v[948] = v[3738] * v[377] + v[1572] * v[946];
	}
	else {
		v[943] = 0e0;
		v[947] = 0e0;
		v[948] = 0e0;
	};
	v[998] = (*epsn)*v[360];
	v[999] = (*epsn)*v[361];
	v[1000] = (*epsn)*v[362];
	v[5689] = 0e0;
	v[5690] = 0e0;
	v[5691] = 0e0;
	v[5692] = 0e0;
	v[5693] = 0e0;
	v[5694] = 0e0;
	v[5695] = v[998];
	v[5696] = v[999];
	v[5697] = v[1000];
	v[5698] = 0e0;
	v[5699] = 0e0;
	v[5700] = 0e0;
	v[5701] = 0e0;
	v[5702] = 0e0;
	v[5703] = 0e0;
	v[5704] = 0e0;
	v[5705] = 0e0;
	v[5706] = 0e0;
	v[5707] = v[998];
	v[5708] = v[999];
	v[5709] = v[1000];
	v[5710] = 0e0;
	v[5711] = 0e0;
	v[5712] = 0e0;
	v[5713] = 0e0;
	v[5714] = 0e0;
	v[5715] = 0e0;
	v[5716] = 0e0;
	v[5717] = 0e0;
	v[5718] = 0e0;
	v[5719] = 0e0;
	v[5720] = 0e0;
	v[5721] = 0e0;
	v[5722] = 0e0;
	v[5723] = 0e0;
	v[5724] = 0e0;
	v[1001] = (*cn)*(v[2630] * v[287] + v[2626] * v[297] + v[2622] * v[307] + v[2618] * v[313] + v[2614] * v[323]
		+ v[2610] * v[333] + v[2606] * v[339] + v[2602] * v[349] + v[2598] * v[359] + v[3700] * v[375] + (-v[279] + v[3294]
			- v[3654])*v[982] + v[274] * v[984] + v[277] * v[985]);
	v[1002] = (*cn)*(v[2631] * v[287] + v[2627] * v[297] + v[2623] * v[307] + v[2619] * v[313] + v[2615] * v[323]
		+ v[2611] * v[333] + v[2607] * v[339] + v[2603] * v[349] + v[2599] * v[359] + v[3701] * v[376] + v[273] * v[984]
		+ v[276] * v[985] + (-v[280] + v[3290] - v[3649])*v[986]);
	v[1003] = (*cn)*(v[2632] * v[287] + v[2628] * v[297] + v[2624] * v[307] + v[2620] * v[313] + v[2616] * v[323]
		+ v[2612] * v[333] + v[2608] * v[339] + v[2604] * v[349] + v[2600] * v[359] + v[3702] * v[377] + (-v[281] + v[3852]
			)*v[990]);
	v[4272] = v[1439] * (v[1001] * v[194] + v[1002] * v[197] + v[1003] * v[200]) - (v[1001] * v[193] + v[1002] * v[196]
		+ v[1003] * v[199])*v[993];
	v[7955] = 0e0;
	v[7956] = 0e0;
	v[7957] = 0e0;
	v[7958] = 0e0;
	v[7959] = 0e0;
	v[7960] = 0e0;
	v[7961] = v[1001];
	v[7962] = v[1002];
	v[7963] = v[1003];
	v[7964] = 0e0;
	v[7965] = 0e0;
	v[7966] = 0e0;
	v[7967] = 0e0;
	v[7968] = 0e0;
	v[7969] = 0e0;
	v[7970] = 0e0;
	v[7971] = 0e0;
	v[7972] = 0e0;
	v[4270] = v[1439] * (v[1001] * v[185] + v[1002] * v[188] + v[1003] * v[191]) - (v[1001] * v[184] + v[1002] * v[187]
		+ v[1003] * v[190])*v[993];
	v[7973] = v[1001];
	v[7974] = v[1002];
	v[7975] = v[1003];
	v[7976] = 0e0;
	v[7977] = 0e0;
	v[7978] = 0e0;
	v[7979] = 0e0;
	v[7980] = 0e0;
	v[7981] = 0e0;
	v[7982] = 0e0;
	v[7983] = 0e0;
	v[7984] = 0e0;
	v[7985] = 0e0;
	v[7986] = 0e0;
	v[7987] = 0e0;
	v[7988] = 0e0;
	v[7989] = 0e0;
	v[7990] = 0e0;
	v[3578] = v[1003] * v[1269] + v[1002] * v[1273] + v[1001] * v[1276];
	v[3577] = v[1003] * v[1268] + v[1002] * v[1272] + v[1001] * v[1275];
	v[3564] = -(v[1001] * v[245]) - v[1002] * v[248] - v[1003] * v[251];
	v[1004] = v[1001] + v[998];
	v[1005] = v[1002] + v[999];
	v[1006] = v[1000] + v[1003];
	v[2382] = (v[1004] * v[1004]) + (v[1005] * v[1005]) + (v[1006] * v[1006]);
	v[1007] = (*epst)*v[943];
	v[1008] = (*epst)*v[947];
	v[1009] = (*epst)*v[948];
	v[1013] = v[1007] - (*ct)*(v[1441] * v[273] + v[1445] * v[274] + v[1449] * v[275] + v[1465] * v[276] + v[1469] * v[277]
		+ v[1473] * v[278] + v[1489] * v[279] + v[1493] * v[280] + v[1497] * v[281] + v[1453] * v[287] + v[1457] * v[297]
		+ v[1461] * v[307] + v[1477] * v[313] + v[1481] * v[323] + v[1485] * v[333] + v[1501] * v[339] + v[1505] * v[349]
		+ v[1509] * v[359]);
	v[1014] = v[1008] - (*ct)*(v[1442] * v[273] + v[1446] * v[274] + v[1450] * v[275] + v[1466] * v[276] + v[1470] * v[277]
		+ v[1474] * v[278] + v[1490] * v[279] + v[1494] * v[280] + v[1498] * v[281] + v[1454] * v[287] + v[1458] * v[297]
		+ v[1462] * v[307] + v[1478] * v[313] + v[1482] * v[323] + v[1486] * v[333] + v[1502] * v[339] + v[1506] * v[349]
		+ v[1510] * v[359]);
	v[1015] = v[1009] - (*ct)*(v[1443] * v[273] + v[1447] * v[274] + v[1451] * v[275] + v[1467] * v[276] + v[1471] * v[277]
		+ v[1475] * v[278] + v[1491] * v[279] + v[1495] * v[280] + v[1499] * v[281] + v[1455] * v[287] + v[1459] * v[297]
		+ v[1463] * v[307] + v[1479] * v[313] + v[1483] * v[323] + v[1487] * v[333] + v[1503] * v[339] + v[1507] * v[349]
		+ v[1511] * v[359]);
	v[2379] = (v[1013] * v[1013]) + (v[1014] * v[1014]) + (v[1015] * v[1015]);
	if ((*stick)) {
		b1017 = sqrt((v[1013] * v[1013]) + (v[1014] * v[1014]) + (v[1015] * v[1015])) <= (*mus)*sqrt((v[1004] * v[1004]) +
			(v[1005] * v[1005]) + (v[1006] * v[1006]));
		if (b1017) {
			v[1019] = v[1013];
			v[1020] = v[1014];
			v[1021] = v[1015];
			v[1022] = 1e0;
		}
		else {
			v[2381] = 1e0 / sqrt(v[2379]);
			v[2383] = 1e0 / sqrt(v[2382]);
			v[1032] = sqrt(v[2382]);
			v[1023] = sqrt(v[2379]);
			if (v[1023] > 0.1e-5) {
				v010 = 1e0 / v[1023]; v011 = (-(v010 / v[1023])); v012 = (2e0*v010) / Power(v[1023], 2
				);
			}
			else {
				v010 = (24000000e0 - (-1e0 + 1000000e0*v[1023])*(71999994e0 - 0.71999982e14*v[1023] + 0.6e19*Power(v[1023]
					, 3) + 0.23999982e20*(v[1023] * v[1023]))) / 24e0;
				v011 = (-500000e0 / 3e0)*(35999994e0 - 0.47999982e14*v[1023] + 0.6e19*Power(v[1023], 3) + 0.17999982e20*
					(v[1023] * v[1023]));
				v012 = 0.1e13*(7999997e0 - 0.5999994e13*v[1023] - 0.3e13*(v[1023] * v[1023]));
			};
			v[1027] = v011;
			v[1028] = v010;
			v[2380] = (*mud)*v[1028] * v[1032];
			v[1019] = v[1013] * v[2380];
			v[1020] = v[1014] * v[2380];
			v[1021] = v[1015] * v[2380];
			v[1022] = 0e0;
		};
		if (sqrt((v[1007] * v[1007]) + (v[1008] * v[1008]) + (v[1009] * v[1009])) > (*mus)*sqrt((v[1004] * v[1004]) +
			(v[1005] * v[1005]) + (v[1006] * v[1006]))) {
			if ((*epst) > 0.1e-5) {
				v013 = 1e0 / (*epst); v014 = (-(v013 / (*epst))); v015 = (2e0*v013) / Power((*epst), 2
				);
			}
			else {
				v013 = (24000000e0 - (-1e0 + 1000000e0*(*epst))*(71999994e0 - 0.71999982e14*(*epst) + 0.23999982e20*Power(
					(*epst), 2) + 0.6e19*Power((*epst), 3))) / 24e0;
				v014 = (-500000e0 / 3e0)*(35999994e0 - 0.47999982e14*(*epst) + 0.17999982e20*Power((*epst), 2)
					+ 0.6e19*Power((*epst), 3));
				v015 = 0.1e13*(7999997e0 - 0.5999994e13*(*epst) - 0.3e13*Power((*epst), 2));
			};
			v[1041] = sqrt((v[1007] * v[1007]) + (v[1008] * v[1008]) + (v[1009] * v[1009]));
			if (v[1041] > 0.1e-5) {
				v016 = 1e0 / v[1041]; v017 = (-(v016 / v[1041])); v018 = (2e0*v016) / Power(v[1041], 2
				);
			}
			else {
				v016 = (24000000e0 - (-1e0 + 1000000e0*v[1041])*(71999994e0 - 0.71999982e14*v[1041] + 0.6e19*Power(v[1041]
					, 3) + 0.23999982e20*(v[1041] * v[1041]))) / 24e0;
				v017 = (-500000e0 / 3e0)*(35999994e0 - 0.47999982e14*v[1041] + 0.6e19*Power(v[1041], 3) + 0.17999982e20*
					(v[1041] * v[1041]));
				v018 = 0.1e13*(7999997e0 - 0.5999994e13*v[1041] - 0.3e13*(v[1041] * v[1041]));
			};
			v[1048] = -((*mud)*v013*v016*sqrt(v[2382]));
			v[1047] = v[1007] * v[1048] + v[943];
			v[1049] = v[1008] * v[1048] + v[947];
			v[1050] = v[1009] * v[1048] + v[948];
		}
		else {
			v[1047] = 0e0;
			v[1049] = 0e0;
			v[1050] = 0e0;
		};
	}
	else {
		b1051 = sqrt((v[1013] * v[1013]) + (v[1014] * v[1014]) + (v[1015] * v[1015])) <= (*mud)*sqrt((v[1004] * v[1004]) +
			(v[1005] * v[1005]) + (v[1006] * v[1006]));
		if (b1051) {
			v[1019] = v[1013];
			v[1020] = v[1014];
			v[1021] = v[1015];
			v[1022] = 1e0;
		}
		else {
			v[2389] = 1e0 / sqrt(v[2379]);
			v[2393] = 1e0 / sqrt(v[2382]);
			v[1062] = sqrt(v[2382]);
			v[3822] = (*mud)*v[1062];
			v[1053] = sqrt(v[2379]);
			if (v[1053] > 0.1e-5) {
				v019 = 1e0 / v[1053]; v020 = (-(v019 / v[1053])); v021 = (2e0*v019) / Power(v[1053], 2
				);
			}
			else {
				v019 = (24000000e0 - (-1e0 + 1000000e0*v[1053])*(71999994e0 - 0.71999982e14*v[1053] + 0.6e19*Power(v[1053]
					, 3) + 0.23999982e20*(v[1053] * v[1053]))) / 24e0;
				v020 = (-500000e0 / 3e0)*(35999994e0 - 0.47999982e14*v[1053] + 0.6e19*Power(v[1053], 3) + 0.17999982e20*
					(v[1053] * v[1053]));
				v021 = 0.1e13*(7999997e0 - 0.5999994e13*v[1053] - 0.3e13*(v[1053] * v[1053]));
			};
			v[1057] = v020;
			v[1058] = v019;
			v[3703] = (*mud)*v[1058] * v[1062];
			v[1019] = v[1013] * v[3703];
			v[1020] = v[1014] * v[3703];
			v[1021] = v[1015] * v[3703];
			v[1022] = 0e0;
		};
		if (sqrt((v[1007] * v[1007]) + (v[1008] * v[1008]) + (v[1009] * v[1009])) > (*mud)*sqrt((v[1004] * v[1004]) +
			(v[1005] * v[1005]) + (v[1006] * v[1006]))) {
			if ((*epst) > 0.1e-5) {
				v022 = 1e0 / (*epst); v023 = (-(v022 / (*epst))); v024 = (2e0*v022) / Power((*epst), 2
				);
			}
			else {
				v022 = (24000000e0 - (-1e0 + 1000000e0*(*epst))*(71999994e0 - 0.71999982e14*(*epst) + 0.23999982e20*Power(
					(*epst), 2) + 0.6e19*Power((*epst), 3))) / 24e0;
				v023 = (-500000e0 / 3e0)*(35999994e0 - 0.47999982e14*(*epst) + 0.17999982e20*Power((*epst), 2)
					+ 0.6e19*Power((*epst), 3));
				v024 = 0.1e13*(7999997e0 - 0.5999994e13*(*epst) - 0.3e13*Power((*epst), 2));
			};
			v[1071] = sqrt((v[1007] * v[1007]) + (v[1008] * v[1008]) + (v[1009] * v[1009]));
			if (v[1071] > 0.1e-5) {
				v025 = 1e0 / v[1071]; v026 = (-(v025 / v[1071])); v027 = (2e0*v025) / Power(v[1071], 2
				);
			}
			else {
				v025 = (24000000e0 - (-1e0 + 1000000e0*v[1071])*(71999994e0 - 0.71999982e14*v[1071] + 0.6e19*Power(v[1071]
					, 3) + 0.23999982e20*(v[1071] * v[1071]))) / 24e0;
				v026 = (-500000e0 / 3e0)*(35999994e0 - 0.47999982e14*v[1071] + 0.6e19*Power(v[1071], 3) + 0.17999982e20*
					(v[1071] * v[1071]));
				v027 = 0.1e13*(7999997e0 - 0.5999994e13*v[1071] - 0.3e13*(v[1071] * v[1071]));
			};
			v[1077] = -((*mud)*v022*v025*sqrt(v[2382]));
			v[1047] = v[1007] * v[1077] + v[943];
			v[1049] = v[1008] * v[1077] + v[947];
			v[1050] = v[1009] * v[1077] + v[948];
		}
		else {
			v[1047] = 0e0;
			v[1049] = 0e0;
			v[1050] = 0e0;
		};
	};
	fn[0] = v[1004];
	fn[1] = v[1005];
	fn[2] = v[1006];
	ft[0] = v[1019];
	ft[1] = v[1020];
	ft[2] = v[1021];
	(*stickupdated) = v[1022];
	gtpupdated[0] = -v[1047] + v[943];
	gtpupdated[1] = -v[1049] + v[947];
	gtpupdated[2] = -v[1050] + v[948];
	v[1089] = -((*epsn)*(v[246] * v[360] + v[249] * v[361] + v[252] * v[362]));
	v[1090] = -((*epsn)*(v[244] * v[360] + v[247] * v[361] + v[250] * v[362]));
	v[1382] = -(v[1089] * v[3633]) + v[1090] * v[3704];
	v[1091] = v[1000] * v[949];
	v[1092] = v[1000] * v[950];
	v[1156] = v[1092] * v[225];
	v[1093] = v[1000] * v[951];
	v[1162] = v[1093] * v[225];
	v[1094] = v[949] * v[999];
	v[1157] = v[1094] * v[225];
	v[1095] = v[950] * v[999];
	v[1159] = -(v[1095] * v[225]) / 2e0;
	v[1096] = v[951] * v[999];
	v[1165] = v[1096] * v[225];
	v[1097] = v[949] * v[998];
	v[1163] = v[1097] * v[225];
	v[1098] = v[950] * v[998];
	v[1166] = v[1098] * v[225];
	v[1099] = v[951] * v[998];
	v[1164] = -(v[1099] * v[225]) / 2e0;
	v[1100] = v[1000] * v[952];
	v[1101] = v[1000] * v[953];
	v[1145] = v[1101] * v[165];
	v[1102] = v[1000] * v[954];
	v[1151] = v[1102] * v[165];
	v[1103] = v[952] * v[999];
	v[1146] = v[1103] * v[165];
	v[1104] = v[953] * v[999];
	v[1148] = -(v[1104] * v[165]) / 2e0;
	v[1105] = v[954] * v[999];
	v[1154] = v[1105] * v[165];
	v[1106] = v[952] * v[998];
	v[1152] = v[1106] * v[165];
	v[1107] = v[953] * v[998];
	v[1155] = v[1107] * v[165];
	v[1108] = v[954] * v[998];
	v[1153] = -(v[1108] * v[165]) / 2e0;
	v[1109] = v[1000] * v[955];
	v[1110] = v[1000] * v[956];
	v[1134] = v[1110] * v[146];
	v[1111] = v[1000] * v[957];
	v[1140] = v[1111] * v[146];
	v[1112] = v[955] * v[999];
	v[1135] = v[1112] * v[146];
	v[1113] = v[956] * v[999];
	v[1137] = -(v[1113] * v[146]) / 2e0;
	v[1114] = v[957] * v[999];
	v[1143] = v[1114] * v[146];
	v[1115] = v[955] * v[998];
	v[1141] = v[1115] * v[146];
	v[1116] = v[956] * v[998];
	v[1144] = v[1116] * v[146];
	v[1117] = v[957] * v[998];
	v[1142] = -(v[1117] * v[146]) / 2e0;
	v[1118] = v[1090] * v[410] + v[1089] * v[412];
	v[1119] = v[1156] + v[1157];
	v[1425] = v[1119] / 2e0;
	v[1421] = (v[1162] + v[1163]) / 2e0;
	v[1121] = v[1165] + v[1166];
	v[1420] = v[1121] / 2e0;
	v[1122] = (v[1099] * v[662]) / 2e0 + v[1098] * v[667] + v[1097] * v[673] + v[1096] * v[677] + (v[1095] * v[681]) / 2e0
		+ v[1094] * v[686] + v[1093] * v[690] + v[1092] * v[694] + (v[1091] * v[699]) / 2e0;
	v[1428] = v[1159] + v[1164] - 4e0*v[1122] * v[1252];
	v[1424] = -v[1159] + v[1428] - (v[1091] * v[225]) / 2e0;
	v[1419] = v[1159] - v[1164] + v[1424];
	v[1123] = v[1336] * v[1389] + v[1335] * v[1391];
	v[1124] = -((*epsn)*(-(v[360] * v[3705]) - v[361] * v[3706] - v[362] * v[3707])) / 2e0;
	v[1125] = v[1145] + v[1146];
	v[1415] = v[1125] / 2e0;
	v[1411] = (v[1151] + v[1152]) / 2e0;
	v[1127] = v[1154] + v[1155];
	v[1410] = v[1127] / 2e0;
	v[1128] = (v[1108] * v[491]) / 2e0 + v[1107] * v[496] + v[1106] * v[502] + v[1105] * v[506] + (v[1104] * v[510]) / 2e0
		+ v[1103] * v[515] + v[1102] * v[519] + v[1101] * v[523] + (v[1100] * v[528]) / 2e0;
	v[1418] = v[1148] + v[1153] - 4e0*v[1128] * v[1238];
	v[1414] = -v[1148] + v[1418] - (v[1100] * v[165]) / 2e0;
	v[1409] = v[1148] - v[1153] + v[1414];
	v[1129] = v[1134] + v[1135];
	v[1405] = v[1129] / 2e0;
	v[1401] = (v[1140] + v[1141]) / 2e0;
	v[1131] = v[1143] + v[1144];
	v[1400] = v[1131] / 2e0;
	v[1132] = (v[1117] * v[446]) / 2e0 + v[1116] * v[451] + v[1115] * v[457] + v[1114] * v[461] + (v[1113] * v[465]) / 2e0
		+ v[1112] * v[470] + v[1111] * v[474] + v[1110] * v[478] + (v[1109] * v[483]) / 2e0;
	v[1408] = v[1137] + v[1142] - 4e0*v[1132] * v[1224];
	v[1404] = -v[1137] + v[1408] - (v[1109] * v[146]) / 2e0;
	v[1399] = v[1137] - v[1142] + v[1404];
	v[4821] = v[210] * v[998];
	v[4822] = v[210] * v[999];
	v[4823] = v[1000] * v[210];
	v[4824] = v[1134] - v[1135] + 2e0*dA[3] * v[1399] + dA[5] * v[1401] + v[1131] * v[426];
	v[4825] = (dA[5] * v[1129]) / 2e0 + (dA[3] * v[1131]) / 2e0 - v[1140] + v[1141] + 2e0*dA[4] * v[1404];
	v[4826] = v[1143] - v[1144] + dA[3] * v[1401] + 2e0*dA[5] * v[1408] + v[1129] * v[426];
	v[4827] = v[211] * v[998];
	v[4828] = v[211] * v[999];
	v[4829] = v[1000] * v[211];
	v[4830] = v[1145] - v[1146] + 2e0*dA[9] * v[1409] + dA[11] * v[1411] + v[1127] * v[432];
	v[4831] = (dA[11] * v[1125]) / 2e0 + (dA[9] * v[1127]) / 2e0 - v[1151] + v[1152] + 2e0*dA[10] * v[1414];
	v[4832] = v[1154] - v[1155] + dA[9] * v[1411] + 2e0*dA[11] * v[1418] + v[1125] * v[432];
	v[4833] = -v[998];
	v[4834] = -v[999];
	v[4835] = -v[1000];
	v[4836] = v[1156] - v[1157] + 2e0*dB[3] * v[1419] + dB[5] * v[1421] + v[1121] * v[438];
	v[4837] = (dB[5] * v[1119]) / 2e0 + (dB[3] * v[1121]) / 2e0 - v[1162] + v[1163] + 2e0*dB[4] * v[1424];
	v[4838] = v[1165] - v[1166] + dB[3] * v[1421] + 2e0*dB[5] * v[1428] + v[1119] * v[438];
	v[1133] = v[1383] * v[1440] - v[1382] * v[994];
	for (i1087 = 1; i1087 <= 18; i1087++) {
		i3719 = (i1087 == 16 ? 1 : 0);
		i3718 = (i1087 == 17 ? 1 : 0);
		i3717 = (i1087 == 18 ? 1 : 0);
		i3716 = (i1087 == 10 ? 1 : 0);
		i3715 = (i1087 == 11 ? 1 : 0);
		i3714 = (i1087 == 12 ? 1 : 0);
		i3713 = (i1087 == 4 ? 1 : 0);
		i3712 = (i1087 == 5 ? 1 : 0);
		i3711 = (i1087 == 6 ? 1 : 0);
		v[1267] = v[4874 + i1087];
		v[3710] = v[1267] * v[998];
		v[3709] = v[1267] * v[999];
		v[3708] = v[1000] * v[1267];
		v[1333] = v[211] * v[3708];
		v[1330] = v[210] * v[3708];
		v[1327] = v[211] * v[3709];
		v[1324] = v[210] * v[3709];
		v[1321] = v[211] * v[3710];
		v[1318] = v[210] * v[3710];
		v[1266] = v[4856 + i1087];
		v[1296] = -(v[1266] * v[998]) / 2e0;
		v[1288] = -(v[1266] * v[999]) / 2e0;
		v[1282] = -(v[1000] * v[1266]) / 2e0;
		v[1194] = v[4838 + i1087];
		v[3727] = v[1194] * v[994];
		v[3720] = -(v[1194] * v[1440]);
		v[1271] = (*epsn)*v[3720];
		v[1193] = v[4892 + i1087];
		v[3726] = v[1089] * v[1193];
		v[3724] = v[1090] * v[1193];
		v[1172] = v[4932 + i1087];
		v[1227] = -4e0*v[1172] * v[1224];
		v[1173] = v[4950 + i1087];
		v[1174] = v[4968 + i1087];
		v[1175] = v[4986 + i1087];
		v[1176] = v[5004 + i1087];
		v[1241] = -4e0*v[1176] * v[1238];
		v[1177] = v[5022 + i1087];
		v[1178] = v[5040 + i1087];
		v[1179] = v[5058 + i1087];
		v[1180] = v[5076 + i1087];
		v[1255] = -4e0*v[1180] * v[1252];
		v[1181] = v[5094 + i1087];
		v[1182] = v[5112 + i1087];
		v[1183] = v[5130 + i1087];
		v[1184] = v[5148 + i1087];
		v[1185] = v[5220 + i1087];
		v[1186] = v[5292 + i1087];
		v[1187] = v[5310 + i1087];
		v[1188] = v[5382 + i1087];
		v[1189] = v[5454 + i1087];
		v[1190] = v[5472 + i1087];
		v[1191] = v[5544 + i1087];
		v[1192] = v[5616 + i1087];
		v[1195] = v[1193] * v[410] + v[1194] * v[418];
		v[1196] = v[1193] * v[412] + v[1194] * v[420];
		v[1197] = -i3711 + v[1173];
		v[1199] = i3711 + v[1173];
		v[1200] = i3712 + v[1174];
		v[1202] = -i3712 + v[1174];
		v[1203] = -i3713 + v[1175];
		v[1205] = i3713 + v[1175];
		v[1206] = -i3714 + v[1177];
		v[1208] = i3714 + v[1177];
		v[1209] = i3715 + v[1178];
		v[1211] = -i3715 + v[1178];
		v[1212] = -i3716 + v[1179];
		v[1214] = i3716 + v[1179];
		v[1215] = -i3717 + v[1181];
		v[1217] = i3717 + v[1181];
		v[1218] = i3718 + v[1182];
		v[1220] = -i3718 + v[1182];
		v[1221] = -i3719 + v[1183];
		v[1223] = i3719 + v[1183];
		v[1225] = -(v[1184] * v[146]) / 2e0 + v[1172] * v[3188];
		v[1226] = v[1197] * v[146] + v[1227] * v[451];
		v[1228] = v[1200] * v[146] + v[1227] * v[457];
		v[1229] = v[1199] * v[146] + v[1227] * v[461];
		v[1230] = (-(v[1185] * v[146]) + v[1227] * v[465]) / 2e0;
		v[1231] = v[1203] * v[146] + v[1227] * v[470];
		v[1232] = v[1202] * v[146] + v[1227] * v[474];
		v[1233] = v[1000] * v[1232] + v[1225] * v[998] + v[1229] * v[999];
		v[1234] = v[1205] * v[146] + v[1227] * v[478];
		v[1235] = v[1000] * v[1234] + v[1226] * v[998] + v[1230] * v[999];
		v[1236] = (-(v[1186] * v[146]) + v[1227] * v[483]) / 2e0;
		v[1237] = v[1000] * v[1236] + v[1228] * v[998] + v[1231] * v[999];
		v[1239] = -(v[1187] * v[165]) / 2e0 + v[1176] * v[3199];
		v[1240] = v[1206] * v[165] + v[1241] * v[496];
		v[1242] = v[1209] * v[165] + v[1241] * v[502];
		v[1243] = v[1208] * v[165] + v[1241] * v[506];
		v[1244] = (-(v[1188] * v[165]) + v[1241] * v[510]) / 2e0;
		v[1245] = v[1212] * v[165] + v[1241] * v[515];
		v[1246] = v[1211] * v[165] + v[1241] * v[519];
		v[1247] = v[1000] * v[1246] + v[1239] * v[998] + v[1243] * v[999];
		v[1248] = v[1214] * v[165] + v[1241] * v[523];
		v[1249] = v[1000] * v[1248] + v[1240] * v[998] + v[1244] * v[999];
		v[1250] = (-(v[1189] * v[165]) + v[1241] * v[528]) / 2e0;
		v[1251] = v[1000] * v[1250] + v[1242] * v[998] + v[1245] * v[999];
		v[1253] = -(v[1190] * v[225]) / 2e0 + v[1180] * v[3210];
		v[1254] = v[1215] * v[225] + v[1255] * v[667];
		v[1256] = v[1218] * v[225] + v[1255] * v[673];
		v[1257] = v[1217] * v[225] + v[1255] * v[677];
		v[1258] = (-(v[1191] * v[225]) + v[1255] * v[681]) / 2e0;
		v[1259] = v[1221] * v[225] + v[1255] * v[686];
		v[1260] = v[1220] * v[225] + v[1255] * v[690];
		v[1261] = v[1000] * v[1260] + v[1253] * v[998] + v[1257] * v[999];
		v[1262] = v[1223] * v[225] + v[1255] * v[694];
		v[1263] = v[1000] * v[1262] + v[1254] * v[998] + v[1258] * v[999];
		v[1264] = (-(v[1192] * v[225]) + v[1255] * v[699]) / 2e0;
		v[1265] = v[1000] * v[1264] + v[1256] * v[998] + v[1259] * v[999];
		v[1270] = v[1271] * v[251] + (*epsn)*(-(v[1195] * v[250]) - v[1196] * v[252] + v[1266] * v[395] + v[1267] * v[407]
			+ v[5634 + i1087] + v[1264] * v[949] + v[1262] * v[950] + v[1260] * v[951] + v[1250] * v[952] + v[1248] * v[953]
			+ v[1246] * v[954] + v[1236] * v[955] + v[1234] * v[956] + v[1232] * v[957]);
		v[1274] = v[1271] * v[248] + (*epsn)*(-(v[1195] * v[247]) - v[1196] * v[249] + v[1266] * v[392] + v[1267] * v[406]
			+ v[5652 + i1087] + v[1259] * v[949] + v[1258] * v[950] + v[1257] * v[951] + v[1245] * v[952] + v[1244] * v[953]
			+ v[1243] * v[954] + v[1231] * v[955] + v[1230] * v[956] + v[1229] * v[957]);
		v[1277] = v[1271] * v[245] + (*epsn)*(-(v[1195] * v[244]) - v[1196] * v[246] + v[1266] * v[389] + v[1267] * v[405]
			+ v[5670 + i1087] + v[1256] * v[949] + v[1254] * v[950] + v[1253] * v[951] + v[1242] * v[952] + v[1240] * v[953]
			+ v[1239] * v[954] + v[1228] * v[955] + v[1226] * v[956] + v[1225] * v[957]);
		v[1278] = -(v[1000] * v[1196]) - v[1270] * v[263];
		v[1279] = -(v[1270] * v[261]) + v[1000] * v[3720];
		v[1280] = -(v[1000] * v[1195]) - v[1270] * v[260];
		v[1281] = v[1282] + v[1270] * v[210];
		v[1283] = -v[1282] + v[1270] * v[211];
		v[1284] = -(v[1274] * v[263]) - v[1196] * v[999];
		v[1285] = -(v[1274] * v[261]) + v[3720] * v[999];
		v[1286] = -(v[1274] * v[260]) - v[1195] * v[999];
		v[1287] = v[1288] + v[1274] * v[210];
		v[1289] = -v[1288] + v[1274] * v[211];
		v[1290] = -(QABi[0][2] * v[1261]) - QABi[1][2] * v[1263] - QABi[2][2] * v[1265] - v[1277] * v[246] - v[1274] * v[249]
			- v[1270] * v[252];
		v[1291] = -(QABi[0][0] * v[1261]) - QABi[1][0] * v[1263] - QABi[2][0] * v[1265] - v[1277] * v[244] - v[1274] * v[247]
			- v[1270] * v[250];
		v[1292] = -(v[1277] * v[263]) - v[1196] * v[998];
		v[1293] = -(v[1277] * v[261]) + v[3720] * v[998];
		v[1294] = -(v[1277] * v[260]) - v[1195] * v[998];
		v[1295] = v[1296] + v[1277] * v[210];
		v[1297] = -v[1296] + v[1277] * v[211];
		v[1298] = QABi[2][2] * v[1278] + QABi[2][1] * v[1279] + QABi[2][0] * v[1280];
		v[1299] = QABi[1][2] * v[1278] + QABi[1][1] * v[1279] + QABi[1][0] * v[1280];
		v[1300] = QABi[0][2] * v[1278] + QABi[0][1] * v[1279] + QABi[0][0] * v[1280];
		v[1301] = QABi[2][2] * v[1284] + QABi[2][1] * v[1285] + QABi[2][0] * v[1286];
		v[1302] = QABi[1][2] * v[1284] + QABi[1][1] * v[1285] + QABi[1][0] * v[1286];
		v[1303] = QABi[0][2] * v[1284] + QABi[0][1] * v[1285] + QABi[0][0] * v[1286];
		v[1304] = QABi[2][2] * v[1292] + QABi[2][1] * v[1293] + QABi[2][0] * v[1294];
		v[1305] = QABi[1][2] * v[1292] + QABi[1][1] * v[1293] + QABi[1][0] * v[1294];
		v[1306] = QABi[0][2] * v[1292] + QABi[0][1] * v[1293] + QABi[0][0] * v[1294];
		v[1307] = (v[1091] * v[1255] + v[1298] * v[225]) / 2e0;
		v[1308] = v[1092] * v[1255] + v[1299] * v[225];
		v[1309] = v[1093] * v[1255] + v[1300] * v[225];
		v[1310] = v[1094] * v[1255] + v[1301] * v[225];
		v[1312] = v[1096] * v[1255] + v[1303] * v[225];
		v[1313] = v[1097] * v[1255] + v[1304] * v[225];
		v[1314] = v[1098] * v[1255] + v[1305] * v[225];
		v[3733] = 8e0*v[1122] * v[1180] * v[1854] - 4e0*v[1252] * (-(v[1099] * v[1190]) / 2e0 - (v[1095] * v[1191]) / 2e0 -
			(v[1091] * v[1192]) / 2e0 + v[1098] * v[1215] + v[1096] * v[1217] + v[1097] * v[1218] + v[1093] * v[1220]
			+ v[1094] * v[1221] + v[1092] * v[1223] + (v[1306] * v[662]) / 2e0 + v[1305] * v[667] + v[1304] * v[673] + v[1303] * v[677]
			+ (v[1302] * v[681]) / 2e0 + v[1301] * v[686] + v[1300] * v[690] + v[1299] * v[694] + (v[1298] * v[699]) / 2e0);
		v[1423] = -(v[1095] * v[1255]) / 2e0 - (v[1302] * v[225]) / 2e0 + v[3733];
		v[1316] = (v[1099] * v[1255] + v[1306] * v[225]) / 2e0;
		v[1317] = v[1318] * v[1439] + v[1295] * v[215];
		v[1319] = v[1295] * v[214] - v[1318] * v[993];
		v[1320] = v[1321] * v[1439] + v[1297] * v[215];
		v[1322] = v[1297] * v[214] - v[1321] * v[993];
		v[1323] = v[1324] * v[1439] + v[1287] * v[215];
		v[1325] = v[1287] * v[214] - v[1324] * v[993];
		v[1326] = v[1327] * v[1439] + v[1289] * v[215];
		v[1328] = v[1289] * v[214] - v[1327] * v[993];
		v[1329] = v[1330] * v[1439] + v[1281] * v[215];
		v[1331] = v[1281] * v[214] - v[1330] * v[993];
		v[1332] = v[1333] * v[1439] + v[1283] * v[215];
		v[1334] = v[1283] * v[214] - v[1333] * v[993];
		v[1343] = (-(v[1233] * v[1337]) - v[1235] * v[1338] - v[1237] * v[1339] + v[1247] * v[1340] + v[1249] * v[1341]
			+ v[1251] * v[1342] - v[1277] * v[387] + v[1277] * v[388] - v[1274] * v[390] + v[1274] * v[391] - v[1270] * v[393]
			+ v[1270] * v[394] - v[1267] * v[4079] + v[1267] * v[4081] + v[5688 + i1087] - v[5706 + i1087]) / 2e0;
		v[1344] = QBAi[2][1] * v[1332] + QBAi[2][0] * v[1334];
		v[1345] = QBAi[1][1] * v[1332] + QBAi[1][0] * v[1334];
		v[1346] = QBAi[0][1] * v[1332] + QBAi[0][0] * v[1334];
		v[1347] = QBAi[2][1] * v[1326] + QBAi[2][0] * v[1328];
		v[1348] = QBAi[1][1] * v[1326] + QBAi[1][0] * v[1328];
		v[1349] = QBAi[0][1] * v[1326] + QBAi[0][0] * v[1328];
		v[1350] = QBAi[2][1] * v[1320] + QBAi[2][0] * v[1322];
		v[1351] = QBAi[1][1] * v[1320] + QBAi[1][0] * v[1322];
		v[1352] = QBAi[0][1] * v[1320] + QBAi[0][0] * v[1322];
		v[1353] = QAAi[2][1] * v[1329] + QAAi[2][0] * v[1331];
		v[1354] = QAAi[1][1] * v[1329] + QAAi[1][0] * v[1331];
		v[1355] = QAAi[0][1] * v[1329] + QAAi[0][0] * v[1331];
		v[1356] = QAAi[2][1] * v[1323] + QAAi[2][0] * v[1325];
		v[1357] = QAAi[1][1] * v[1323] + QAAi[1][0] * v[1325];
		v[1358] = QAAi[0][1] * v[1323] + QAAi[0][0] * v[1325];
		v[1359] = QAAi[2][1] * v[1317] + QAAi[2][0] * v[1319];
		v[1360] = QAAi[1][1] * v[1317] + QAAi[1][0] * v[1319];
		v[1361] = QAAi[0][1] * v[1317] + QAAi[0][0] * v[1319];
		v[1362] = (v[1100] * v[1241] + v[1344] * v[165]) / 2e0;
		v[1363] = v[1101] * v[1241] + v[1345] * v[165];
		v[1364] = v[1102] * v[1241] + v[1346] * v[165];
		v[1365] = v[1103] * v[1241] + v[1347] * v[165];
		v[1367] = v[1105] * v[1241] + v[1349] * v[165];
		v[1368] = v[1106] * v[1241] + v[1350] * v[165];
		v[1369] = v[1107] * v[1241] + v[1351] * v[165];
		v[3731] = 8e0*v[1128] * v[1176] * v[1852] - 4e0*v[1238] * (-(v[1108] * v[1187]) / 2e0 - (v[1104] * v[1188]) / 2e0 -
			(v[1100] * v[1189]) / 2e0 + v[1107] * v[1206] + v[1105] * v[1208] + v[1106] * v[1209] + v[1102] * v[1211]
			+ v[1103] * v[1212] + v[1101] * v[1214] + (v[1352] * v[491]) / 2e0 + v[1351] * v[496] + v[1350] * v[502] + v[1349] * v[506]
			+ (v[1348] * v[510]) / 2e0 + v[1347] * v[515] + v[1346] * v[519] + v[1345] * v[523] + (v[1344] * v[528]) / 2e0);
		v[1413] = -(v[1104] * v[1241]) / 2e0 - (v[1348] * v[165]) / 2e0 + v[3731];
		v[1371] = (v[1108] * v[1241] + v[1352] * v[165]) / 2e0;
		v[1372] = (v[1109] * v[1227] + v[1353] * v[146]) / 2e0;
		v[1373] = v[1110] * v[1227] + v[1354] * v[146];
		v[1374] = v[1111] * v[1227] + v[1355] * v[146];
		v[1375] = v[1112] * v[1227] + v[1356] * v[146];
		v[1377] = v[1114] * v[1227] + v[1358] * v[146];
		v[1378] = v[1115] * v[1227] + v[1359] * v[146];
		v[1379] = v[1116] * v[1227] + v[1360] * v[146];
		v[3729] = 8e0*v[1132] * v[1172] * v[1850] - 4e0*v[1224] * (-(v[1117] * v[1184]) / 2e0 - (v[1113] * v[1185]) / 2e0 -
			(v[1109] * v[1186]) / 2e0 + v[1116] * v[1197] + v[1114] * v[1199] + v[1115] * v[1200] + v[1111] * v[1202]
			+ v[1112] * v[1203] + v[1110] * v[1205] + (v[1361] * v[446]) / 2e0 + v[1360] * v[451] + v[1359] * v[457] + v[1358] * v[461]
			+ (v[1357] * v[465]) / 2e0 + v[1356] * v[470] + v[1355] * v[474] + v[1354] * v[478] + (v[1353] * v[483]) / 2e0);
		v[1403] = -(v[1113] * v[1227]) / 2e0 - (v[1357] * v[146]) / 2e0 + v[3729];
		v[1381] = (v[1117] * v[1227] + v[1361] * v[146]) / 2e0;
		v[1384] = -(v[1440] * (QABi[0][1] * v[1261] + QABi[1][1] * v[1263] + QABi[2][1] * v[1265] + v[1194] * v[1382]
			+ v[1277] * v[245] + v[1274] * v[248] + v[1270] * v[251])) + v[994] * (-(v[1194] * v[1383]) + (-(v[1291] * v[1438])
				+ v[3724] * v[992])*(*xfac) + (v[1438] * v[3726] + v[1290] * v[992])*(*zfac));
		v[3732] = (v[1309] + v[1313]) / 2e0;
		v[1386] = v[1308] + v[1310];
		v[1387] = v[1312] + v[1314];
		v[1388] = v[992] * (v[3726] * v[3918] + (-(v[1291] * v[262]) + v[1090] * v[3727])*(*xfac)) - v[1438] *
			(v[3724] * v[3725] + (v[1290] * v[262] - v[1089] * v[3727])*(*zfac));
		v[1390] = ((*epsn)*v[1267]) / 2e0;
		v[1392] = v[1439] * (-(v[1389] * v[1390]) + v[1295] * v[185] + v[1287] * v[188] + v[1281] * v[191] + v[1297] * v[194]
			+ v[1289] * v[197] + v[1283] * v[200] + (QAAi[0][1] * v[1233] + QAAi[1][1] * v[1235] + QAAi[2][1] * v[1237])*v[210] +
			(QBAi[0][1] * v[1247] + QBAi[1][1] * v[1249] + QBAi[2][1] * v[1251])*v[211]) - (v[1390] * v[1391] + v[1295] * v[184]
				+ v[1287] * v[187] + v[1281] * v[190] + v[1297] * v[193] + v[1289] * v[196] + v[1283] * v[199] + (QAAi[0][0] * v[1233]
					+ QAAi[1][0] * v[1235] + QAAi[2][0] * v[1237])*v[210] + (QBAi[0][0] * v[1247] + QBAi[1][0] * v[1249]
						+ QBAi[2][0] * v[1251])*v[211])*v[993];
		v[3730] = (v[1364] + v[1368]) / 2e0;
		v[1394] = v[1363] + v[1365];
		v[1395] = v[1367] + v[1369];
		v[3728] = (v[1374] + v[1378]) / 2e0;
		v[1397] = v[1373] + v[1375];
		v[1398] = v[1377] + v[1379];
		v[5797] = 0e0;
		v[5798] = 0e0;
		v[5799] = 0e0;
		v[5800] = 2e0*v[1399];
		v[5801] = v[1400];
		v[5802] = v[1401];
		v[5803] = 0e0;
		v[5804] = 0e0;
		v[5805] = 0e0;
		v[5806] = 0e0;
		v[5807] = 0e0;
		v[5808] = 0e0;
		v[5809] = 0e0;
		v[5810] = 0e0;
		v[5811] = 0e0;
		v[5812] = 0e0;
		v[5813] = 0e0;
		v[5814] = 0e0;
		v[5779] = 0e0;
		v[5780] = 0e0;
		v[5781] = 0e0;
		v[5782] = v[1400];
		v[5783] = 2e0*v[1404];
		v[5784] = v[1405];
		v[5785] = 0e0;
		v[5786] = 0e0;
		v[5787] = 0e0;
		v[5788] = 0e0;
		v[5789] = 0e0;
		v[5790] = 0e0;
		v[5791] = 0e0;
		v[5792] = 0e0;
		v[5793] = 0e0;
		v[5794] = 0e0;
		v[5795] = 0e0;
		v[5796] = 0e0;
		v[5761] = 0e0;
		v[5762] = 0e0;
		v[5763] = 0e0;
		v[5764] = v[1401];
		v[5765] = v[1405];
		v[5766] = 2e0*v[1408];
		v[5767] = 0e0;
		v[5768] = 0e0;
		v[5769] = 0e0;
		v[5770] = 0e0;
		v[5771] = 0e0;
		v[5772] = 0e0;
		v[5773] = 0e0;
		v[5774] = 0e0;
		v[5775] = 0e0;
		v[5776] = 0e0;
		v[5777] = 0e0;
		v[5778] = 0e0;
		v[5743] = 0e0;
		v[5744] = 0e0;
		v[5745] = 0e0;
		v[5746] = 0e0;
		v[5747] = 0e0;
		v[5748] = 0e0;
		v[5749] = 0e0;
		v[5750] = 0e0;
		v[5751] = 0e0;
		v[5752] = 2e0*v[1409];
		v[5753] = v[1410];
		v[5754] = v[1411];
		v[5755] = 0e0;
		v[5756] = 0e0;
		v[5757] = 0e0;
		v[5758] = 0e0;
		v[5759] = 0e0;
		v[5760] = 0e0;
		v[5725] = 0e0;
		v[5726] = 0e0;
		v[5727] = 0e0;
		v[5728] = 0e0;
		v[5729] = 0e0;
		v[5730] = 0e0;
		v[5731] = 0e0;
		v[5732] = 0e0;
		v[5733] = 0e0;
		v[5734] = v[1410];
		v[5735] = 2e0*v[1414];
		v[5736] = v[1415];
		v[5737] = 0e0;
		v[5738] = 0e0;
		v[5739] = 0e0;
		v[5740] = 0e0;
		v[5741] = 0e0;
		v[5742] = 0e0;
		v[5815] = 0e0;
		v[5816] = 0e0;
		v[5817] = 0e0;
		v[5818] = 0e0;
		v[5819] = 0e0;
		v[5820] = 0e0;
		v[5821] = 0e0;
		v[5822] = 0e0;
		v[5823] = 0e0;
		v[5824] = v[1411];
		v[5825] = v[1415];
		v[5826] = 2e0*v[1418];
		v[5827] = 0e0;
		v[5828] = 0e0;
		v[5829] = 0e0;
		v[5830] = 0e0;
		v[5831] = 0e0;
		v[5832] = 0e0;
		v[5833] = 0e0;
		v[5834] = 0e0;
		v[5835] = 0e0;
		v[5836] = 0e0;
		v[5837] = 0e0;
		v[5838] = 0e0;
		v[5839] = 0e0;
		v[5840] = 0e0;
		v[5841] = 0e0;
		v[5842] = 0e0;
		v[5843] = 0e0;
		v[5844] = 0e0;
		v[5845] = 0e0;
		v[5846] = 0e0;
		v[5847] = 0e0;
		v[5848] = 2e0*v[1419];
		v[5849] = v[1420];
		v[5850] = v[1421];
		v[5851] = 0e0;
		v[5852] = 0e0;
		v[5853] = 0e0;
		v[5854] = 0e0;
		v[5855] = 0e0;
		v[5856] = 0e0;
		v[5857] = 0e0;
		v[5858] = 0e0;
		v[5859] = 0e0;
		v[5860] = 0e0;
		v[5861] = 0e0;
		v[5862] = 0e0;
		v[5863] = 0e0;
		v[5864] = 0e0;
		v[5865] = 0e0;
		v[5866] = v[1420];
		v[5867] = 2e0*v[1424];
		v[5868] = v[1425];
		v[5869] = 0e0;
		v[5870] = 0e0;
		v[5871] = 0e0;
		v[5872] = 0e0;
		v[5873] = 0e0;
		v[5874] = 0e0;
		v[5875] = 0e0;
		v[5876] = 0e0;
		v[5877] = 0e0;
		v[5878] = 0e0;
		v[5879] = 0e0;
		v[5880] = 0e0;
		v[5881] = 0e0;
		v[5882] = 0e0;
		v[5883] = 0e0;
		v[5884] = v[1421];
		v[5885] = v[1425];
		v[5886] = 2e0*v[1428];
		v[5887] = v[1295];
		v[5888] = v[1287];
		v[5889] = v[1281];
		v[5890] = v[1373] - v[1375] + 2e0*dA[3] * (-v[1372] + v[1403]) + dA[5] * v[3728] + v[1398] * v[426] + v[5796 + i1087];
		v[5891] = -v[1374] + v[1378] + (dA[5] * v[1397]) / 2e0 + (dA[3] * v[1398]) / 2e0 + 2e0*dA[4] * (-v[1372] - v[1381]
			+ v[3729]) + v[5778 + i1087];
		v[5892] = v[1377] - v[1379] + 2e0*dA[5] * (-v[1381] + v[1403]) + dA[3] * v[3728] + v[1397] * v[426] + v[5760 + i1087];
		v[5893] = v[1297];
		v[5894] = v[1289];
		v[5895] = v[1283];
		v[5896] = v[1363] - v[1365] + 2e0*dA[9] * (-v[1362] + v[1413]) + dA[11] * v[3730] + v[1395] * v[432] + v[5742 + i1087];
		v[5897] = -v[1364] + v[1368] + (dA[11] * v[1394]) / 2e0 + (dA[9] * v[1395]) / 2e0 + 2e0*dA[10] * (-v[1362] - v[1371]
			+ v[3731]) + v[5724 + i1087];
		v[5898] = v[1367] - v[1369] + 2e0*dA[11] * (-v[1371] + v[1413]) + dA[9] * v[3730] + v[1394] * v[432] + v[5814 + i1087];
		v[5899] = -v[1277];
		v[5900] = -v[1274];
		v[5901] = -v[1270];
		v[5902] = v[1308] - v[1310] + 2e0*dB[3] * (-v[1307] + v[1423]) + dB[5] * v[3732] + v[1387] * v[438] + v[5832 + i1087];
		v[5903] = -v[1309] + v[1313] + (dB[5] * v[1386]) / 2e0 + (dB[3] * v[1387]) / 2e0 + 2e0*dB[4] * (-v[1307] - v[1316]
			+ v[3733]) + v[5850 + i1087];
		v[5904] = v[1312] - v[1314] + 2e0*dB[5] * (-v[1316] + v[1423]) + dB[3] * v[3732] + v[1386] * v[438] + v[5868 + i1087];
		Rc[i1087 - 1] += v[1118] * v[1193] + v[1133] * v[1194] + v[1124] * v[1266] + v[1123] * v[1267] + v[4820 + i1087];
		for (i1169 = i1087; i1169 <= 18; i1169++) {
			v[1433] = v[1384] * v[4838 + i1169] + v[1343] * v[4856 + i1169] + v[1392] * v[4874 + i1169] + v[1388] * v[4892 + i1169]
				+ v[5886 + i1169];
			Kc[i1087 - 1][i1169 - 1] += v[1433];
			if (i1087 != i1169) {
				Kc[i1169 - 1][i1087 - 1] += v[1433];
			}
			else {
			};
		};/* end for */
	};/* end for */
	v[1516] = 0e0;
	v[1517] = 0e0;
	v[1518] = 0e0;
	b1519 = (*stick);
	if (b1519) {
		b1520 = b1017;
		if (b1520) {
			v[1518] = 0e0;
			v[1517] = 0e0;
			v[1516] = 0e0;
		}
		else {
		};
	}
	else {
		b1521 = b1051;
		if (b1521) {
			v[1518] = 0e0;
			v[1517] = 0e0;
			v[1516] = 0e0;
		}
		else {
		};
	};
	v[3737] = (*ct)*v[1516];
	v[3832] = -v[1019] + v[1516] * v[3734];
	v[3736] = (*ct)*v[1517];
	v[3830] = -v[1020] + v[1517] * v[3734];
	v[3735] = (*ct)*v[1518];
	v[3828] = -v[1021] + v[1518] * v[3734];
	v[1525] = v[1529] * v[3735];
	v[1526] = v[1531] * v[3735];
	v[1775] = -(v[1533] * v[3735]) / 2e0;
	v[1528] = v[1535] * v[3735];
	v[3889] = v[1439] * v[1528];
	v[3762] = v[1528] * v[210];
	v[3761] = v[1528] * v[211];
	v[1530] = v[1529] * v[3736];
	v[1532] = v[1531] * v[3736];
	v[1772] = -(v[1533] * v[3736]) / 2e0;
	v[1536] = v[1535] * v[3736];
	v[3890] = v[1439] * v[1536];
	v[3760] = v[1536] * v[210];
	v[3759] = v[1536] * v[211];
	v[1537] = -((*ct)*(v[1509] * v[1516] + v[1510] * v[1517] + v[1511] * v[1518]));
	v[1751] = v[1537] * v[1873];
	v[1733] = v[1537] * v[3045];
	v[1538] = -((*ct)*(v[1505] * v[1516] + v[1506] * v[1517] + v[1507] * v[1518]));
	v[1740] = v[1538] * v[3044];
	v[3745] = v[1740] + v[1537] * v[1872];
	v[3744] = v[1733] + v[1538] * v[1870];
	v[1539] = -((*ct)*(v[1501] * v[1516] + v[1502] * v[1517] + v[1503] * v[1518]));
	v[1749] = v[1539] * v[3043];
	v[3998] = v[1749] + v[1751];
	v[3746] = v[1749] + v[1538] * v[1871];
	v[3989] = v[1751] + v[3746];
	v[1742] = v[1539] * v[3645];
	v[3995] = v[1742] + v[3745];
	v[3990] = v[1740] + v[1742];
	v[1736] = v[1539] * v[3647];
	v[3999] = v[1736] + v[3744];
	v[3991] = v[1733] + v[1736];
	v[1540] = -((*ct)*(v[1485] * v[1516] + v[1486] * v[1517] + v[1487] * v[1518]));
	v[1792] = v[1540] * v[1867];
	v[1778] = v[1540] * v[3042];
	v[1541] = -((*ct)*(v[1481] * v[1516] + v[1482] * v[1517] + v[1483] * v[1518]));
	v[1785] = v[1541] * v[3041];
	v[3764] = v[1785] + v[1540] * v[1866];
	v[3763] = v[1778] + v[1541] * v[1864];
	v[1542] = -((*ct)*(v[1477] * v[1516] + v[1478] * v[1517] + v[1479] * v[1518]));
	v[1790] = v[1542] * v[3040];
	v[3985] = v[1790] + v[1792];
	v[3765] = v[1790] + v[1541] * v[1865];
	v[3976] = v[1792] + v[3765];
	v[1787] = v[1542] * v[3640];
	v[3982] = v[1787] + v[3764];
	v[3977] = v[1785] + v[1787];
	v[1781] = v[1542] * v[3642];
	v[3986] = v[1781] + v[3763];
	v[3978] = v[1778] + v[1781];
	v[1543] = -((*ct)*(v[1461] * v[1516] + v[1462] * v[1517] + v[1463] * v[1518]));
	v[1810] = v[1543] * v[1861];
	v[1796] = v[1543] * v[3039];
	v[1544] = -((*ct)*(v[1457] * v[1516] + v[1458] * v[1517] + v[1459] * v[1518]));
	v[1803] = v[1544] * v[3038];
	v[3767] = v[1803] + v[1543] * v[1860];
	v[3766] = v[1796] + v[1544] * v[1858];
	v[1545] = -((*ct)*(v[1453] * v[1516] + v[1454] * v[1517] + v[1455] * v[1518]));
	v[1808] = v[1545] * v[3037];
	v[5927] = (*ct)*(-(v[1441] * v[1516]) - v[1442] * v[1517] - v[1443] * v[1518]);
	v[5928] = (*ct)*(-(v[1445] * v[1516]) - v[1446] * v[1517] - v[1447] * v[1518]);
	v[5929] = (*ct)*(-(v[1449] * v[1516]) - v[1450] * v[1517] - v[1451] * v[1518]);
	v[5930] = v[1808] + v[1544] * v[3635] + v[1543] * v[3637];
	v[5931] = v[1803] + v[1543] * v[1858] + v[1545] * v[1859];
	v[5932] = v[1796] + v[1544] * v[1860] + v[1545] * v[1861];
	v[5933] = (*ct)*(-(v[1465] * v[1516]) - v[1466] * v[1517] - v[1467] * v[1518]);
	v[5934] = (*ct)*(-(v[1469] * v[1516]) - v[1470] * v[1517] - v[1471] * v[1518]);
	v[5935] = (*ct)*(-(v[1473] * v[1516]) - v[1474] * v[1517] - v[1475] * v[1518]);
	v[5936] = v[1790] + v[1541] * v[3640] + v[1540] * v[3642];
	v[5937] = v[1785] + v[1540] * v[1864] + v[1542] * v[1865];
	v[5938] = v[1778] + v[1541] * v[1866] + v[1542] * v[1867];
	v[5939] = (*ct)*(-(v[1489] * v[1516]) - v[1490] * v[1517] - v[1491] * v[1518]);
	v[5940] = (*ct)*(-(v[1493] * v[1516]) - v[1494] * v[1517] - v[1495] * v[1518]);
	v[5941] = (*ct)*(-(v[1497] * v[1516]) - v[1498] * v[1517] - v[1499] * v[1518]);
	v[5942] = v[1749] + v[1538] * v[3645] + v[1537] * v[3647];
	v[5943] = v[1740] + v[1537] * v[1870] + v[1539] * v[1871];
	v[5944] = v[1733] + v[1538] * v[1872] + v[1539] * v[1873];
	v[3973] = v[1808] + v[1810];
	v[3768] = v[1808] + v[1544] * v[1859];
	v[3964] = v[1810] + v[3768];
	v[1805] = v[1545] * v[3635];
	v[3970] = v[1805] + v[3767];
	v[3965] = v[1803] + v[1805];
	v[1799] = v[1545] * v[3637];
	v[3974] = v[1799] + v[3766];
	v[3966] = v[1796] + v[1799];
	v[1546] = v[1529] * v[3737];
	v[1547] = v[1531] * v[3737];
	v[1769] = -(v[1533] * v[3737]) / 2e0;
	v[1549] = v[1535] * v[3737];
	v[3891] = v[1439] * v[1549];
	v[3758] = v[1549] * v[210];
	v[3757] = v[1549] * v[211];
	v[1550] = (*epst)*v[1518];
	v[2303] = -(v[1550] * v[377]);
	v[1551] = (*epst)*v[1517];
	v[2305] = -(v[1551] * v[376]);
	v[2302] = v[2303] + v[2305];
	v[1552] = (*epst)*v[1516];
	v[2306] = -(v[1552] * v[375]);
	v[2307] = v[2305] + v[2306];
	v[2304] = v[2303] + v[2306];
	v[1553] = 0e0;
	v[1554] = 0e0;
	v[1555] = 0e0;
	v[1556] = 0e0;
	v[1557] = 0e0;
	v[1558] = 0e0;
	v[1559] = 0e0;
	v[1560] = 0e0;
	v[1561] = 0e0;
	v[1562] = 0e0;
	v[1563] = 0e0;
	v[1564] = 0e0;
	b1565 = (*previouscontact);
	if (b1565) {
		v[1566] = -(v[1550] * v[946]);
		v[1567] = -(v[1551] * v[945]);
		v[1569] = v[1552] * v[1568] + v[2302] * v[375];
		v[1571] = v[1551] * v[1570] + v[2304] * v[376];
		v[1573] = v[1550] * v[1572] + v[2307] * v[377];
		v[1555] = v[1550] * v[3738] + v[2307] * v[946];
		v[1554] = v[1551] * v[3739] + v[2304] * v[945];
		v[1575] = -(v[1552] * v[944]);
		v[1553] = -(v[1552] * v[3740]) + v[2302] * v[944];
		v[1556] = gti[0] * v[1569];
		v[1557] = gti[1] * v[1569];
		v[1558] = gti[2] * v[1569];
		v[1578] = -v[1569];
		v[1579] = -(v[1569] * v[259]);
		v[1580] = -(v[1569] * v[257]);
		v[1581] = -(v[1569] * v[256]);
		v[1582] = v[1569] * v[209];
		v[1583] = v[1569] * v[208];
		v[1584] = v[1582] * v[213];
		v[1585] = v[1582] * v[212];
		v[1586] = v[1583] * v[213];
		v[1587] = v[1583] * v[212];
		v[1559] = gti[0] * v[1571];
		v[1560] = gti[1] * v[1571];
		v[1561] = gti[2] * v[1571];
		v[1588] = -v[1571];
		v[1589] = -(v[1571] * v[259]);
		v[1590] = -(v[1571] * v[257]);
		v[1591] = -(v[1571] * v[256]);
		v[1592] = v[1571] * v[209];
		v[1593] = v[1571] * v[208];
		v[1594] = v[1592] * v[213];
		v[1595] = v[1592] * v[212];
		v[1596] = v[1593] * v[213];
		v[1597] = v[1593] * v[212];
		v[1562] = gti[0] * v[1573];
		v[1563] = gti[1] * v[1573];
		v[1564] = gti[2] * v[1573];
		v[1598] = -v[1573];
		v[1599] = -(v[1573] * v[259]);
		v[1600] = -(v[1573] * v[257]);
		v[1601] = -(v[1573] * v[256]);
		v[1602] = v[1573] * v[209];
		v[1603] = v[1573] * v[208];
		v[1604] = v[1602] * v[213];
		v[1605] = v[1602] * v[212];
		v[1606] = v[1603] * v[213];
		v[1607] = v[1603] * v[212];
	}
	else {
		v[1587] = 0e0;
		v[1586] = 0e0;
		v[1597] = 0e0;
		v[1596] = 0e0;
		v[1607] = 0e0;
		v[1606] = 0e0;
		v[1585] = 0e0;
		v[1584] = 0e0;
		v[1595] = 0e0;
		v[1594] = 0e0;
		v[1605] = 0e0;
		v[1604] = 0e0;
		v[1583] = 0e0;
		v[1593] = 0e0;
		v[1603] = 0e0;
		v[1582] = 0e0;
		v[1592] = 0e0;
		v[1602] = 0e0;
		v[1581] = 0e0;
		v[1580] = 0e0;
		v[1579] = 0e0;
		v[1591] = 0e0;
		v[1590] = 0e0;
		v[1589] = 0e0;
		v[1601] = 0e0;
		v[1600] = 0e0;
		v[1599] = 0e0;
		v[1578] = 0e0;
		v[1588] = 0e0;
		v[1598] = 0e0;
		v[1575] = 0e0;
		v[1567] = 0e0;
		v[1566] = 0e0;
	};
	b1608 = b909;
	if (b1608) {
		v[1634] = -(v[1564] * v[926]) / 2e0;
		v[1633] = -(v[1560] * v[926]) / 2e0;
		v[1632] = v[1563] * v[926];
		v[1631] = v[1561] * v[926];
		v[1627] = v[1562] * v[926];
		v[1626] = v[1558] * v[926];
		v[1623] = v[1559] * v[926];
		v[1622] = v[1557] * v[926];
		v[1609] = v[1631] + v[1632];
		v[1610] = v[1626] + v[1627];
		v[1611] = v[1622] + v[1623];
		v[1621] = (v[1556] * v[1612]) / 2e0 + v[1557] * v[1613] + v[1558] * v[1614] + v[1559] * v[1615] + (v[1560] * v[1616])
			/ 2e0 + v[1561] * v[1617] + v[1562] * v[1618] + v[1563] * v[1619] + (v[1564] * v[1620]) / 2e0;
		v[2269] = (-4e0*v[1621]) / Power(v[1624], 2) + v[1633] + v[1634];
		v[2268] = -v[1633] + v[2269] - (v[1556] * v[926]) / 2e0;
		v[2267] = v[1633] - v[1634] + v[2268];
		v[1625] = (-2e0*v[1622] + 2e0*v[1623] + v[1610] * v[923] + v[1609] * v[924] + 4e0*v[2267] * v[925]) / 2e0;
		v[1630] = (2e0*v[1626] - 2e0*v[1627] + v[1611] * v[923] + 4e0*v[2268] * v[924] + v[1609] * v[925]) / 2e0;
		v[1635] = (-2e0*v[1631] + 2e0*v[1632] + 4e0*v[2269] * v[923] + v[1611] * v[924] + v[1610] * v[925]) / 2e0;
		v[3741] = v[1635] * v[911] + v[1630] * v[912] + v[1625] * v[913];
		v[2255] = v[3741] * v[922];
		v[2252] = v[3741] * v[916];
		v[1638] = v[2252] * v[921] + v[2255] / (Power(cos(v[1636]), 2)*sqrt(v[2256]));
		v[3812] = v[1638] / v[914];
		v[3742] = v[1638] / v[914];
		v[1639] = v[1625] * v[3696] + v[3742] * v[913];
		v[1641] = v[1630] * v[3696] + v[3742] * v[912];
		v[1642] = v[1635] * v[3696] + v[3742] * v[911];
		v[1553] = v[1553] - v[1639] * v[385] + v[1641] * v[386];
		v[1555] = v[1555] - v[1641] * v[384] + v[1642] * v[385];
		v[1554] = v[1554] + v[1639] * v[384] - v[1642] * v[386];
	}
	else {
	};
	v[1555] = v[1555] + 2e0*v[1566] * v[377];
	v[1554] = v[1554] + 2e0*v[1567] * v[376];
	v[1553] = v[1553] + 2e0*v[1575] * v[375];
	v[2635] = v[1553] * v[360] + v[1554] * v[361] + v[1555] * v[362];
	v[1643] = v[2635] * v[373];
	v[3743] = v[1643] / v[1644];
	v[1645] = v[1555] * v[374] + v[362] * v[3743];
	v[1647] = v[1554] * v[374] + v[361] * v[3743];
	v[1648] = v[1553] * v[374] + v[360] * v[3743];
	v[1598] = v[1598] - v[1645];
	v[1588] = v[1588] - v[1647];
	v[1578] = v[1578] - v[1648];
	v[1651] = v[1537] * v[2351];
	v[1652] = v[1537] * v[2350];
	v[1653] = v[1537] * v[2349];
	v[1656] = v[1538] * v[2346];
	v[1657] = v[1537] * v[346] + v[1538] * v[348];
	v[3747] = v[1657] * v[230];
	v[1660] = v[1538] * v[2344];
	v[1661] = v[1651] + v[1656];
	v[1662] = -v[1651] + v[1656];
	v[3959] = 8e0*v[1662];
	v[1664] = v[1538] * v[2345] + v[3747] / v[337];
	v[1667] = v[1539] * v[2339];
	v[1668] = v[1537] * v[340] + v[1539] * v[348];
	v[3748] = v[1668] * v[235];
	v[1669] = v[1538] * v[340] + v[1539] * v[346];
	v[3749] = v[1669] * v[239];
	v[4249] = -(v[3644] * v[3747]) - v[336] * v[3748] - v[335] * v[3749] - v[1537] * (v[348] * v[350] + v[340] * v[3905]
		+ v[346] * v[3906] + v[340] * v[3908] + v[346] * v[3910]) - v[1539] * (v[338] * (v[340] + v[229] * v[346] + v[230] * v[348])
			+ v[346] * v[3909] + v[348] * v[3913]) - v[1538] * (v[341] * v[346] + v[340] * v[3907] + v[340] * v[3911] + v[348] * v[3912]
				+ v[348] * v[3914]);
	v[1671] = v[1539] * v[2341] + v[3748] / v[337];
	v[1672] = v[1664] + v[1671];
	v[1673] = -v[1664] + v[1671];
	v[3957] = 8e0*v[1673];
	v[1675] = v[1539] * v[2340] + v[3749] / v[337];
	v[1676] = v[1652] + v[1675];
	v[1677] = -v[1652] + v[1675];
	v[4245] = v[1662] * v[437] - v[1677] * v[440] + v[1673] * v[442];
	v[3958] = -8e0*v[1677];
	v[7623] = 0e0;
	v[7624] = 0e0;
	v[7625] = 0e0;
	v[7626] = 0e0;
	v[7627] = 0e0;
	v[7628] = 0e0;
	v[7629] = 0e0;
	v[7630] = 0e0;
	v[7631] = 0e0;
	v[7632] = 0e0;
	v[7633] = 0e0;
	v[7634] = 0e0;
	v[7635] = 0e0;
	v[7636] = 0e0;
	v[7637] = 0e0;
	v[7638] = v[3959];
	v[7639] = v[3958];
	v[7640] = v[3957];
	v[1678] = v[1540] * v[2336];
	v[1679] = v[1540] * v[2335];
	v[1680] = v[1540] * v[2334];
	v[1683] = v[1541] * v[2331];
	v[1684] = v[1540] * v[320] + v[1541] * v[322];
	v[3769] = v[1684] * v[170];
	v[1687] = v[1541] * v[2329];
	v[1688] = v[1678] + v[1683];
	v[1689] = -v[1678] + v[1683];
	v[3953] = 8e0*v[1689];
	v[1691] = v[1541] * v[2330] + v[3769] / v[311];
	v[1694] = v[1542] * v[2324];
	v[1695] = v[1540] * v[314] + v[1542] * v[322];
	v[3770] = v[1695] * v[175];
	v[1696] = v[1541] * v[314] + v[1542] * v[320];
	v[3771] = v[1696] * v[179];
	v[4255] = -(v[3639] * v[3769]) - v[310] * v[3770] - v[309] * v[3771] - v[1540] * (v[322] * v[324] + v[314] * v[3919]
		+ v[320] * v[3920] + v[314] * v[3922] + v[320] * v[3924]) - v[1542] * (v[312] * (v[314] + v[169] * v[320] + v[170] * v[322])
			+ v[320] * v[3923] + v[322] * v[3927]) - v[1541] * (v[315] * v[320] + v[314] * v[3921] + v[314] * v[3925] + v[322] * v[3926]
				+ v[322] * v[3928]);
	v[1698] = v[1542] * v[2326] + v[3770] / v[311];
	v[1699] = v[1691] + v[1698];
	v[1700] = -v[1691] + v[1698];
	v[3951] = 8e0*v[1700];
	v[1702] = v[1542] * v[2325] + v[3771] / v[311];
	v[1703] = v[1679] + v[1702];
	v[1704] = -v[1679] + v[1702];
	v[4251] = v[1689] * v[431] - v[1704] * v[434] + v[1700] * v[436];
	v[3952] = -8e0*v[1704];
	v[7641] = 0e0;
	v[7642] = 0e0;
	v[7643] = 0e0;
	v[7644] = 0e0;
	v[7645] = 0e0;
	v[7646] = 0e0;
	v[7647] = 0e0;
	v[7648] = 0e0;
	v[7649] = 0e0;
	v[7650] = v[3953];
	v[7651] = v[3952];
	v[7652] = v[3951];
	v[7653] = 0e0;
	v[7654] = 0e0;
	v[7655] = 0e0;
	v[7656] = 0e0;
	v[7657] = 0e0;
	v[7658] = 0e0;
	v[1705] = v[1543] * v[2321];
	v[1706] = v[1543] * v[2320];
	v[1707] = v[1543] * v[2319];
	v[1710] = v[1544] * v[2316];
	v[1711] = v[1543] * v[294] + v[1544] * v[296];
	v[3779] = v[151] * v[1711];
	v[1714] = v[1544] * v[2314];
	v[1715] = v[1705] + v[1710];
	v[1716] = -v[1705] + v[1710];
	v[3947] = 8e0*v[1716];
	v[1718] = v[1544] * v[2315] + v[3779] / v[285];
	v[1721] = v[1545] * v[2309];
	v[1722] = v[1543] * v[288] + v[1545] * v[296];
	v[3780] = v[156] * v[1722];
	v[1723] = v[1544] * v[288] + v[1545] * v[294];
	v[3781] = v[160] * v[1723];
	v[4261] = -(v[3634] * v[3779]) - v[284] * v[3780] - v[283] * v[3781] - v[1543] * (v[296] * v[298] + v[288] * v[3932]
		+ v[294] * v[3933] + v[288] * v[3935] + v[294] * v[3937]) - v[1545] * (v[286] * (v[288] + v[150] * v[294] + v[151] * v[296])
			+ v[294] * v[3936] + v[296] * v[3940]) - v[1544] * (v[289] * v[294] + v[288] * v[3934] + v[288] * v[3938] + v[296] * v[3939]
				+ v[296] * v[3941]);
	v[1725] = v[1545] * v[2311] + v[3780] / v[285];
	v[1726] = v[1718] + v[1725];
	v[1727] = -v[1718] + v[1725];
	v[3945] = 8e0*v[1727];
	v[1729] = v[1545] * v[2310] + v[3781] / v[285];
	v[1730] = v[1706] + v[1729];
	v[1731] = -v[1706] + v[1729];
	v[4257] = v[1716] * v[425] - v[1731] * v[428] + v[1727] * v[430];
	v[3946] = -8e0*v[1731];
	v[7677] = 0e0;
	v[7678] = 0e0;
	v[7679] = 0e0;
	v[7680] = v[3947];
	v[7681] = v[3946];
	v[7682] = v[3945];
	v[7683] = 0e0;
	v[7684] = 0e0;
	v[7685] = 0e0;
	v[7686] = 0e0;
	v[7687] = 0e0;
	v[7688] = 0e0;
	v[7689] = 0e0;
	v[7690] = 0e0;
	v[7691] = 0e0;
	v[7692] = 0e0;
	v[7693] = 0e0;
	v[7694] = 0e0;
	v[1599] = v[1599] - v[1645] * v[263] + v[1525] * v[412] + v[1526] * v[420];
	v[1600] = v[1440] * v[1526] + v[1600] - v[1645] * v[261];
	v[1601] = v[1601] - v[1645] * v[260] + v[1525] * v[410] + v[1526] * v[418];
	v[1732] = QABi[2][2] * v[1599] + QABi[2][1] * v[1600] + QABi[2][0] * v[1601] + v[348] * v[3999];
	v[1735] = QABi[1][2] * v[1599] + QABi[1][1] * v[1600] + QABi[1][0] * v[1601] + v[1669] * v[3647] + v[346] * v[3744];
	v[1737] = QABi[0][2] * v[1599] + QABi[0][1] * v[1600] + QABi[0][0] * v[1601] + v[340] * v[3991];
	v[1589] = v[1589] - v[1647] * v[263] + v[1530] * v[412] + v[1532] * v[420];
	v[1590] = v[1440] * v[1532] + v[1590] - v[1647] * v[261];
	v[1591] = v[1591] - v[1647] * v[260] + v[1530] * v[410] + v[1532] * v[418];
	v[1738] = QABi[2][2] * v[1589] + QABi[2][1] * v[1590] + QABi[2][0] * v[1591] + v[1668] * v[3645] + v[348] * v[3745];
	v[1741] = QABi[1][2] * v[1589] + QABi[1][1] * v[1590] + QABi[1][0] * v[1591] + v[346] * v[3995];
	v[1743] = QABi[0][2] * v[1589] + QABi[0][1] * v[1590] + QABi[0][0] * v[1591] + v[340] * v[3990];
	v[1579] = v[1579] - v[1648] * v[263] + v[1546] * v[412] + v[1547] * v[420];
	v[1580] = v[1440] * v[1547] + v[1580] - v[1648] * v[261];
	v[1581] = v[1581] - v[1648] * v[260] + v[1546] * v[410] + v[1547] * v[418];
	v[1748] = QABi[2][2] * v[1579] + QABi[2][1] * v[1580] + QABi[2][0] * v[1581] + v[1657] * v[1871] + v[348] * v[3998];
	v[1750] = QABi[1][2] * v[1579] + QABi[1][1] * v[1580] + QABi[1][0] * v[1581] + v[346] * v[3746];
	v[1753] = QABi[0][2] * v[1579] + QABi[0][1] * v[1580] + QABi[0][0] * v[1581] + v[340] * v[3989];
	v[1754] = -(v[1672] * v[666]) + 2e0*v[1667] * v[669] + v[1676] * v[672];
	v[3962] = -v[1754] / 2e0;
	v[1755] = 2e0*v[1660] * v[666] - v[1672] * v[669] - v[1661] * v[672];
	v[3961] = v[1755] / 2e0;
	v[1756] = -(v[1661] * v[666]) + v[1676] * v[669] + 2e0*v[1653] * v[672];
	v[3960] = -v[1756] / 2e0;
	v[7569] = 0e0;
	v[7570] = 0e0;
	v[7571] = 0e0;
	v[7572] = 0e0;
	v[7573] = 0e0;
	v[7574] = 0e0;
	v[7575] = 0e0;
	v[7576] = 0e0;
	v[7577] = 0e0;
	v[7578] = 0e0;
	v[7579] = 0e0;
	v[7580] = 0e0;
	v[7581] = 0e0;
	v[7582] = 0e0;
	v[7583] = 0e0;
	v[7584] = v[3962];
	v[7585] = v[3961];
	v[7586] = v[3960];
	v[1757] = (v[1732] * v[225]) / 2e0;
	v[1758] = v[1735] * v[225];
	v[1759] = v[1737] * v[225];
	v[1760] = v[1738] * v[225];
	v[1761] = (v[1741] * v[225]) / 2e0;
	v[1762] = v[1743] * v[225];
	v[1763] = v[1748] * v[225];
	v[1764] = v[1750] * v[225];
	v[1766] = 1e0 / Power(v[337], 2);
	v[3756] = -(v[1766] * v[346]);
	v[3755] = -(v[1766] * v[348]);
	v[3754] = -(v[1766] * v[340]);
	v[3753] = -(v[1766] * v[335]);
	v[3752] = -(v[1537] * v[1766]);
	v[3751] = -(v[1766] * v[336]);
	v[3750] = -(v[1766] * v[3644]);
	v[2952] = -(v[1766] * (v[230] * v[338] + v[3913]));
	v[2951] = -(v[1766] * (v[229] * v[338] + v[3909]));
	v[2950] = -(v[1766] * (v[3907] + v[3911]));
	v[2949] = -(v[1766] * (v[3912] + v[3914]));
	v[2948] = -(v[1766] * (v[3905] + v[3908]));
	v[2947] = -(v[1766] * (v[3906] + v[3910]));
	v[2946] = -(v[1766] * v[338]);
	v[3917] = v[2946] * v[340];
	v[2945] = -(v[1766] * v[341]);
	v[3916] = v[2945] * v[346];
	v[3915] = -(v[1766] * v[348] * v[350]);
	v[2943] = -(v[1766] * v[3747]);
	v[2942] = -(v[1766] * v[3748]);
	v[2941] = -(v[1766] * v[3749]);
	v[2940] = v[1669] * v[3753];
	v[2939] = v[1668] * v[3751];
	v[2938] = v[1657] * v[3750];
	v[2785] = v[1538] * v[3750];
	v[2784] = v[3646] * v[3752];
	v[2782] = v[1539] * v[2946];
	v[4225] = v[2938] + (v[2782] + v[2784])*v[348];
	v[3882] = v[2782] + v[2785];
	v[4226] = v[2784] + v[3882];
	v[2778] = v[1539] * v[3751];
	v[2775] = v[1538] * v[2945];
	v[4224] = v[2775] + v[2778];
	v[2774] = v[3648] * v[3752];
	v[3881] = v[2774] + v[2775];
	v[4223] = v[2778] + v[3881];
	v[4222] = v[2939] + v[348] * v[3881];
	v[2770] = v[1539] * v[3753];
	v[2768] = -(v[1538] * v[1766] * v[347]);
	v[2767] = v[350] * v[3752];
	v[4221] = v[2767] + v[2770];
	v[3880] = v[2767] + v[2768];
	v[4220] = v[2940] + v[346] * v[3880];
	v[4219] = v[2770] + v[3880];
	v[2670] = v[237] * v[3754];
	v[2667] = v[232] * v[3754];
	v[2664] = v[230] * v[3755];
	v[2663] = v[229] * v[3756];
	v[2658] = v[235] * v[3755];
	v[2657] = v[234] * v[3756];
	v[4205] = v[2657] + v[2667];
	v[3861] = v[2657] + v[2658];
	v[4203] = v[2667] + v[3861];
	v[2655] = v[228] * v[3754];
	v[4201] = v[2655] + v[2664];
	v[3862] = v[2655] + v[2663];
	v[4204] = v[2664] + v[3862];
	v[2651] = v[240] * v[3755];
	v[4206] = v[2651] + v[2670];
	v[2650] = v[239] * v[3756];
	v[3859] = v[2650] + v[2651];
	v[4202] = v[2670] + v[3859];
	v[2348] = v[2948] * v[340] + v[2947] * v[346] + v[3915];
	v[2343] = v[2950] * v[340] + v[2949] * v[348] + v[3916];
	v[2338] = v[2951] * v[346] + v[2952] * v[348] + v[3917];
	v[2167] = v[230] * v[3750];
	v[2165] = v[235] * v[3751];
	v[2161] = v[239] * v[3753];
	v[2937] = v[1653] + v[1660] + v[1667] + v[1669] * v[2161] + v[1668] * v[2165] + v[1657] * v[2167] + v[1539] * v[2338]
		+ v[1538] * v[2343] + v[1537] * v[2348];
	v[1767] = (-(dB[3] * v[1754]) + dB[4] * v[1755] - dB[5] * v[1756] + 4e0*v[225] * v[2937] + v[1753] * v[662]
		+ 2e0*v[1750] * v[667] + 2e0*v[1748] * v[673] + 2e0*v[1743] * v[677] + v[1741] * v[681] + 2e0*v[1738] * v[686]
		+ 2e0*v[1737] * v[690] + 2e0*v[1735] * v[694] + v[1732] * v[699]) / 2e0;
	v[1855] = -v[1761] - 4e0*v[1252] * v[1767] + v[1677] * v[2057] + v[1673] * v[2059] + v[1662] * v[2061];
	v[3036] = v[1855] - (v[1753] * v[225]) / 2e0;
	v[3034] = -v[1757] + v[1761] + v[3036];
	v[3029] = -v[1757] + v[1855];
	v[1768] = v[1769] + v[1648] * v[210];
	v[1770] = -v[1769] + v[1648] * v[211];
	v[1583] = v[1583] + v[1768];
	v[1582] = v[1582] + v[1770];
	v[1771] = v[1772] + v[1647] * v[210];
	v[1773] = -v[1772] + v[1647] * v[211];
	v[1593] = v[1593] + v[1771];
	v[1592] = v[1592] + v[1773];
	v[1774] = v[1775] + v[1645] * v[210];
	v[1776] = -v[1775] + v[1645] * v[211];
	v[1603] = v[1603] + v[1774];
	v[1602] = v[1602] + v[1776];
	v[1584] = v[1584] + v[1770] * v[215] + v[1439] * v[3757];
	v[1585] = v[1585] + v[1770] * v[214] - v[3757] * v[993];
	v[1586] = v[1586] + v[1768] * v[215] + v[1439] * v[3758];
	v[1587] = v[1587] + v[1768] * v[214] - v[3758] * v[993];
	v[1594] = v[1594] + v[1773] * v[215] + v[1439] * v[3759];
	v[1595] = v[1595] + v[1773] * v[214] - v[3759] * v[993];
	v[1596] = v[1596] + v[1771] * v[215] + v[1439] * v[3760];
	v[1597] = v[1597] + v[1771] * v[214] - v[3760] * v[993];
	v[1604] = v[1604] + v[1776] * v[215] + v[1439] * v[3761];
	v[1605] = v[1605] + v[1776] * v[214] - v[3761] * v[993];
	v[1606] = v[1606] + v[1774] * v[215] + v[1439] * v[3762];
	v[1607] = v[1607] + v[1774] * v[214] - v[3762] * v[993];
	v[1777] = QBAi[2][1] * v[1604] + QBAi[2][0] * v[1605] + v[322] * v[3986];
	v[1780] = QBAi[1][1] * v[1604] + QBAi[1][0] * v[1605] + v[1696] * v[3642] + v[320] * v[3763];
	v[1782] = QBAi[0][1] * v[1604] + QBAi[0][0] * v[1605] + v[314] * v[3978];
	v[1783] = QBAi[2][1] * v[1594] + QBAi[2][0] * v[1595] + v[1695] * v[3640] + v[322] * v[3764];
	v[1786] = QBAi[1][1] * v[1594] + QBAi[1][0] * v[1595] + v[320] * v[3982];
	v[1788] = QBAi[0][1] * v[1594] + QBAi[0][0] * v[1595] + v[314] * v[3977];
	v[1789] = QBAi[2][1] * v[1584] + QBAi[2][0] * v[1585] + v[1684] * v[1865] + v[322] * v[3985];
	v[1791] = QBAi[1][1] * v[1584] + QBAi[1][0] * v[1585] + v[320] * v[3765];
	v[1794] = QBAi[0][1] * v[1584] + QBAi[0][0] * v[1585] + v[314] * v[3976];
	v[1795] = QAAi[2][1] * v[1606] + QAAi[2][0] * v[1607] + v[296] * v[3974];
	v[1798] = QAAi[1][1] * v[1606] + QAAi[1][0] * v[1607] + v[1723] * v[3637] + v[294] * v[3766];
	v[1800] = QAAi[0][1] * v[1606] + QAAi[0][0] * v[1607] + v[288] * v[3966];
	v[1801] = QAAi[2][1] * v[1596] + QAAi[2][0] * v[1597] + v[1722] * v[3635] + v[296] * v[3767];
	v[1804] = QAAi[1][1] * v[1596] + QAAi[1][0] * v[1597] + v[294] * v[3970];
	v[1806] = QAAi[0][1] * v[1596] + QAAi[0][0] * v[1597] + v[288] * v[3965];
	v[1807] = QAAi[2][1] * v[1586] + QAAi[2][0] * v[1587] + v[1711] * v[1859] + v[296] * v[3973];
	v[1809] = QAAi[1][1] * v[1586] + QAAi[1][0] * v[1587] + v[294] * v[3768];
	v[1812] = QAAi[0][1] * v[1586] + QAAi[0][0] * v[1587] + v[288] * v[3964];
	v[1813] = -(v[1699] * v[495]) + 2e0*v[1694] * v[498] + v[1703] * v[501];
	v[3956] = -v[1813] / 2e0;
	v[1814] = 2e0*v[1687] * v[495] - v[1699] * v[498] - v[1688] * v[501];
	v[3955] = v[1814] / 2e0;
	v[1815] = -(v[1688] * v[495]) + v[1703] * v[498] + 2e0*v[1680] * v[501];
	v[3954] = -v[1815] / 2e0;
	v[7587] = 0e0;
	v[7588] = 0e0;
	v[7589] = 0e0;
	v[7590] = 0e0;
	v[7591] = 0e0;
	v[7592] = 0e0;
	v[7593] = 0e0;
	v[7594] = 0e0;
	v[7595] = 0e0;
	v[7596] = v[3956];
	v[7597] = v[3955];
	v[7598] = v[3954];
	v[7599] = 0e0;
	v[7600] = 0e0;
	v[7601] = 0e0;
	v[7602] = 0e0;
	v[7603] = 0e0;
	v[7604] = 0e0;
	v[1816] = (v[165] * v[1777]) / 2e0;
	v[1817] = v[165] * v[1780];
	v[1818] = v[165] * v[1782];
	v[1819] = v[165] * v[1783];
	v[1820] = (v[165] * v[1786]) / 2e0;
	v[1821] = v[165] * v[1788];
	v[1822] = v[165] * v[1789];
	v[1823] = v[165] * v[1791];
	v[1825] = 1e0 / Power(v[311], 2);
	v[3778] = -(v[1825] * v[320]);
	v[3777] = -(v[1825] * v[322]);
	v[3776] = -(v[1825] * v[314]);
	v[3775] = -(v[1825] * v[309]);
	v[3774] = -(v[1540] * v[1825]);
	v[3773] = -(v[1825] * v[310]);
	v[3772] = -(v[1825] * v[3639]);
	v[2981] = -(v[1825] * (v[170] * v[312] + v[3927]));
	v[2980] = -(v[1825] * (v[169] * v[312] + v[3923]));
	v[2979] = -(v[1825] * (v[3921] + v[3925]));
	v[2978] = -(v[1825] * (v[3926] + v[3928]));
	v[2977] = -(v[1825] * (v[3919] + v[3922]));
	v[2976] = -(v[1825] * (v[3920] + v[3924]));
	v[2975] = -(v[1825] * v[312]);
	v[3931] = v[2975] * v[314];
	v[2974] = -(v[1825] * v[315]);
	v[3930] = v[2974] * v[320];
	v[3929] = -(v[1825] * v[322] * v[324]);
	v[2972] = -(v[1825] * v[3769]);
	v[2971] = -(v[1825] * v[3770]);
	v[2970] = -(v[1825] * v[3771]);
	v[2969] = v[1696] * v[3775];
	v[2968] = v[1695] * v[3773];
	v[2967] = v[1684] * v[3772];
	v[2882] = v[1541] * v[3772];
	v[2881] = v[3641] * v[3774];
	v[2879] = v[1542] * v[2975];
	v[4234] = v[2967] + (v[2879] + v[2881])*v[322];
	v[3898] = v[2879] + v[2882];
	v[4235] = v[2881] + v[3898];
	v[2875] = v[1542] * v[3773];
	v[2872] = v[1541] * v[2974];
	v[4233] = v[2872] + v[2875];
	v[2871] = v[3643] * v[3774];
	v[3894] = v[2871] + v[2872];
	v[4232] = v[2875] + v[3894];
	v[4231] = v[2968] + v[322] * v[3894];
	v[2867] = v[1542] * v[3775];
	v[2865] = -(v[1541] * v[1825] * v[321]);
	v[2864] = v[324] * v[3774];
	v[4230] = v[2864] + v[2867];
	v[3893] = v[2864] + v[2865];
	v[4229] = v[2969] + v[320] * v[3893];
	v[4228] = v[2867] + v[3893];
	v[2700] = v[177] * v[3776];
	v[2697] = v[172] * v[3776];
	v[2694] = v[170] * v[3777];
	v[2693] = v[169] * v[3778];
	v[2688] = v[175] * v[3777];
	v[2687] = v[174] * v[3778];
	v[4211] = v[2687] + v[2697];
	v[3869] = v[2687] + v[2688];
	v[4209] = v[2697] + v[3869];
	v[2685] = v[168] * v[3776];
	v[4207] = v[2685] + v[2694];
	v[3870] = v[2685] + v[2693];
	v[4210] = v[2694] + v[3870];
	v[2681] = v[180] * v[3777];
	v[4212] = v[2681] + v[2700];
	v[2680] = v[179] * v[3778];
	v[3867] = v[2680] + v[2681];
	v[4208] = v[2700] + v[3867];
	v[2333] = v[2977] * v[314] + v[2976] * v[320] + v[3929];
	v[2328] = v[2979] * v[314] + v[2978] * v[322] + v[3930];
	v[2323] = v[2980] * v[320] + v[2981] * v[322] + v[3931];
	v[2139] = v[170] * v[3772];
	v[2137] = v[175] * v[3773];
	v[2133] = v[179] * v[3775];
	v[2966] = v[1680] + v[1687] + v[1694] + v[1696] * v[2133] + v[1695] * v[2137] + v[1684] * v[2139] + v[1542] * v[2323]
		+ v[1541] * v[2328] + v[1540] * v[2333];
	v[1826] = (-(dA[9] * v[1813]) + dA[10] * v[1814] - dA[11] * v[1815] + 4e0*v[165] * v[2966] + v[1794] * v[491]
		+ 2e0*v[1791] * v[496] + 2e0*v[1789] * v[502] + 2e0*v[1788] * v[506] + v[1786] * v[510] + 2e0*v[1783] * v[515]
		+ 2e0*v[1782] * v[519] + 2e0*v[1780] * v[523] + v[1777] * v[528]) / 2e0;
	v[1853] = -v[1820] - 4e0*v[1238] * v[1826] + v[1704] * v[1976] + v[1700] * v[1978] + v[1689] * v[1980];
	v[3028] = -(v[165] * v[1794]) / 2e0 + v[1853];
	v[3026] = -v[1816] + v[1820] + v[3028];
	v[3021] = -v[1816] + v[1853];
	v[1827] = -(v[1726] * v[450]) + 2e0*v[1721] * v[453] + v[1730] * v[456];
	v[3950] = -v[1827] / 2e0;
	v[1828] = 2e0*v[1714] * v[450] - v[1726] * v[453] - v[1715] * v[456];
	v[3949] = v[1828] / 2e0;
	v[1829] = -(v[1715] * v[450]) + v[1730] * v[453] + 2e0*v[1707] * v[456];
	v[3948] = -v[1829] / 2e0;
	v[7605] = 0e0;
	v[7606] = 0e0;
	v[7607] = 0e0;
	v[7608] = v[3950];
	v[7609] = v[3949];
	v[7610] = v[3948];
	v[7611] = 0e0;
	v[7612] = 0e0;
	v[7613] = 0e0;
	v[7614] = 0e0;
	v[7615] = 0e0;
	v[7616] = 0e0;
	v[7617] = 0e0;
	v[7618] = 0e0;
	v[7619] = 0e0;
	v[7620] = 0e0;
	v[7621] = 0e0;
	v[7622] = 0e0;
	v[1830] = (v[146] * v[1795]) / 2e0;
	v[1831] = v[146] * v[1798];
	v[1832] = v[146] * v[1800];
	v[1833] = v[146] * v[1801];
	v[1834] = (v[146] * v[1804]) / 2e0;
	v[1835] = v[146] * v[1806];
	v[1836] = v[146] * v[1807];
	v[1837] = v[146] * v[1809];
	v[1839] = 1e0 / Power(v[285], 2);
	v[3788] = -(v[1839] * v[294]);
	v[3787] = -(v[1839] * v[296]);
	v[3786] = -(v[1839] * v[288]);
	v[3785] = -(v[1839] * v[283]);
	v[3784] = -(v[1543] * v[1839]);
	v[3783] = -(v[1839] * v[284]);
	v[3782] = -(v[1839] * v[3634]);
	v[3008] = -(v[1839] * (v[151] * v[286] + v[3940]));
	v[3007] = -(v[1839] * (v[150] * v[286] + v[3936]));
	v[3006] = -(v[1839] * (v[3934] + v[3938]));
	v[3005] = -(v[1839] * (v[3939] + v[3941]));
	v[3004] = -(v[1839] * (v[3932] + v[3935]));
	v[3003] = -(v[1839] * (v[3933] + v[3937]));
	v[3002] = -(v[1839] * v[286]);
	v[3944] = v[288] * v[3002];
	v[3001] = -(v[1839] * v[289]);
	v[3943] = v[294] * v[3001];
	v[3942] = -(v[1839] * v[296] * v[298]);
	v[2999] = -(v[1839] * v[3779]);
	v[2998] = -(v[1839] * v[3780]);
	v[2997] = -(v[1839] * v[3781]);
	v[2996] = v[1723] * v[3785];
	v[2995] = v[1722] * v[3783];
	v[2994] = v[1711] * v[3782];
	v[2915] = v[1544] * v[3782];
	v[2914] = v[3636] * v[3784];
	v[2912] = v[1545] * v[3002];
	v[4242] = (v[2912] + v[2914])*v[296] + v[2994];
	v[3904] = v[2912] + v[2915];
	v[4243] = v[2914] + v[3904];
	v[2908] = v[1545] * v[3783];
	v[2905] = v[1544] * v[3001];
	v[4241] = v[2905] + v[2908];
	v[2904] = v[3638] * v[3784];
	v[3900] = v[2904] + v[2905];
	v[4240] = v[2908] + v[3900];
	v[4239] = v[2995] + v[296] * v[3900];
	v[2900] = v[1545] * v[3785];
	v[2898] = -(v[1544] * v[1839] * v[295]);
	v[2897] = v[298] * v[3784];
	v[4238] = v[2897] + v[2900];
	v[3899] = v[2897] + v[2898];
	v[4237] = v[2996] + v[294] * v[3899];
	v[4236] = v[2900] + v[3899];
	v[2730] = v[158] * v[3786];
	v[2727] = v[153] * v[3786];
	v[2724] = v[151] * v[3787];
	v[2723] = v[150] * v[3788];
	v[2718] = v[156] * v[3787];
	v[2717] = v[155] * v[3788];
	v[4217] = v[2717] + v[2727];
	v[3877] = v[2717] + v[2718];
	v[4215] = v[2727] + v[3877];
	v[2715] = v[149] * v[3786];
	v[4213] = v[2715] + v[2724];
	v[3878] = v[2715] + v[2723];
	v[4216] = v[2724] + v[3878];
	v[2711] = v[161] * v[3787];
	v[4218] = v[2711] + v[2730];
	v[2710] = v[160] * v[3788];
	v[3875] = v[2710] + v[2711];
	v[4214] = v[2730] + v[3875];
	v[2318] = v[294] * v[3003] + v[288] * v[3004] + v[3942];
	v[2313] = v[296] * v[3005] + v[288] * v[3006] + v[3943];
	v[2308] = v[294] * v[3007] + v[296] * v[3008] + v[3944];
	v[2111] = v[151] * v[3782];
	v[2109] = v[156] * v[3783];
	v[2105] = v[160] * v[3785];
	v[2993] = v[1707] + v[1714] + v[1721] + v[1723] * v[2105] + v[1722] * v[2109] + v[1711] * v[2111] + v[1545] * v[2308]
		+ v[1544] * v[2313] + v[1543] * v[2318];
	v[1840] = (-(dA[3] * v[1827]) + dA[4] * v[1828] - dA[5] * v[1829] + 4e0*v[146] * v[2993] + v[1812] * v[446]
		+ 2e0*v[1809] * v[451] + 2e0*v[1807] * v[457] + 2e0*v[1806] * v[461] + v[1804] * v[465] + 2e0*v[1801] * v[470]
		+ 2e0*v[1800] * v[474] + 2e0*v[1798] * v[478] + v[1795] * v[483]) / 2e0;
	v[1851] = -v[1834] - 4e0*v[1224] * v[1840] + v[1731] * v[1949] + v[1727] * v[1951] + v[1716] * v[1953];
	v[3020] = -(v[146] * v[1812]) / 2e0 + v[1851];
	v[3018] = -v[1830] + v[1834] + v[3020];
	v[3013] = -v[1830] + v[1851];
	v[3031] = (v[1759] + v[1763]) / 2e0;
	v[1842] = v[1758] + v[1760];
	v[3035] = v[1842] / 2e0;
	v[1843] = v[1762] + v[1764];
	v[3030] = v[1843] / 2e0;
	v[3023] = (v[1818] + v[1822]) / 2e0;
	v[1845] = v[1817] + v[1819];
	v[3027] = v[1845] / 2e0;
	v[1846] = v[1821] + v[1823];
	v[3022] = v[1846] / 2e0;
	v[3015] = (v[1832] + v[1836]) / 2e0;
	v[1848] = v[1831] + v[1833];
	v[3019] = v[1848] / 2e0;
	v[1849] = v[1835] + v[1837];
	v[5909] = -(v[1019] * v[1441]) - v[1020] * v[1442] - v[1021] * v[1443] + v[1583];
	v[5910] = -(v[1019] * v[1445]) - v[1020] * v[1446] - v[1021] * v[1447] + v[1593];
	v[5911] = -(v[1019] * v[1449]) - v[1020] * v[1450] - v[1021] * v[1451] + v[1603];
	v[5912] = -(v[1019] * v[1453]) - v[1020] * v[1454] - v[1021] * v[1455] + v[1831] - v[1833] + 2e0*dA[3] * v[3013]
		+ dA[5] * v[3015] + v[1224] * v[3947] + v[146] * v[3950] + v[1849] * v[426];
	v[5913] = -(v[1019] * v[1457]) - v[1020] * v[1458] - v[1021] * v[1459] - v[1832] + v[1836] + (dA[5] * v[1848]) / 2e0 +
		(dA[3] * v[1849]) / 2e0 + 2e0*dA[4] * v[3018] + v[1224] * v[3946] + v[146] * v[3949];
	v[5914] = -(v[1019] * v[1461]) - v[1020] * v[1462] - v[1021] * v[1463] + v[1835] - v[1837] + dA[3] * v[3015]
		+ 2e0*dA[5] * v[3020] + v[1224] * v[3945] + v[146] * v[3948] + v[1848] * v[426];
	v[5915] = -(v[1019] * v[1465]) - v[1020] * v[1466] - v[1021] * v[1467] + v[1582];
	v[5916] = -(v[1019] * v[1469]) - v[1020] * v[1470] - v[1021] * v[1471] + v[1592];
	v[5917] = -(v[1019] * v[1473]) - v[1020] * v[1474] - v[1021] * v[1475] + v[1602];
	v[5918] = -(v[1019] * v[1477]) - v[1020] * v[1478] - v[1021] * v[1479] + v[1817] - v[1819] + 2e0*dA[9] * v[3021]
		+ dA[11] * v[3023] + v[1238] * v[3953] + v[165] * v[3956] + v[1846] * v[432];
	v[5919] = -(v[1019] * v[1481]) - v[1020] * v[1482] - v[1021] * v[1483] - v[1818] + v[1822] + (dA[11] * v[1845]) / 2e0 +
		(dA[9] * v[1846]) / 2e0 + 2e0*dA[10] * v[3026] + v[1238] * v[3952] + v[165] * v[3955];
	v[5920] = -(v[1019] * v[1485]) - v[1020] * v[1486] - v[1021] * v[1487] + v[1821] - v[1823] + dA[9] * v[3023]
		+ 2e0*dA[11] * v[3028] + v[1238] * v[3951] + v[165] * v[3954] + v[1845] * v[432];
	v[5921] = -(v[1019] * v[1489]) - v[1020] * v[1490] - v[1021] * v[1491] + v[1578];
	v[5922] = -(v[1019] * v[1493]) - v[1020] * v[1494] - v[1021] * v[1495] + v[1588];
	v[5923] = -(v[1019] * v[1497]) - v[1020] * v[1498] - v[1021] * v[1499] + v[1598];
	v[5924] = -(v[1019] * v[1501]) - v[1020] * v[1502] - v[1021] * v[1503] + v[1758] - v[1760] + 2e0*dB[3] * v[3029]
		+ dB[5] * v[3031] + v[1252] * v[3959] + v[225] * v[3962] + v[1843] * v[438];
	v[5925] = -(v[1019] * v[1505]) - v[1020] * v[1506] - v[1021] * v[1507] - v[1759] + v[1763] + (dB[5] * v[1842]) / 2e0 +
		(dB[3] * v[1843]) / 2e0 + 2e0*dB[4] * v[3034] + v[1252] * v[3958] + v[225] * v[3961];
	v[5926] = -(v[1019] * v[1509]) - v[1020] * v[1510] - v[1021] * v[1511] + v[1762] - v[1764] + dB[3] * v[3031]
		+ 2e0*dB[5] * v[3036] + v[1252] * v[3957] + v[225] * v[3960] + v[1842] * v[438];
	v[3014] = v[1849] / 2e0;
	for (i1514 = 1; i1514 <= 18; i1514++) {
		i3815 = (i1514 == 13 ? 1 : 0);
		i3814 = (i1514 == 14 ? 1 : 0);
		i3813 = (i1514 == 15 ? 1 : 0);
		i3803 = (i1514 == 1 ? 1 : 0);
		i3802 = (i1514 == 7 ? 1 : 0);
		i3801 = (i1514 == 2 ? 1 : 0);
		i3800 = (i1514 == 8 ? 1 : 0);
		i3799 = (i1514 == 3 ? 1 : 0);
		i3798 = (i1514 == 9 ? 1 : 0);
		i3797 = (i1514 == 17 ? 1 : 0);
		v[3856] = (*a4)*i3797;
		i3796 = (i1514 == 16 ? 1 : 0);
		v[3857] = (*a4)*i3796;
		i3795 = (i1514 == 18 ? 1 : 0);
		v[3858] = (*a4)*i3795;
		i3794 = (i1514 == 11 ? 1 : 0);
		v[3864] = (*a4)*i3794;
		i3793 = (i1514 == 10 ? 1 : 0);
		v[3865] = (*a4)*i3793;
		i3792 = (i1514 == 12 ? 1 : 0);
		v[3866] = (*a4)*i3792;
		i3791 = (i1514 == 5 ? 1 : 0);
		v[3872] = (*a4)*i3791;
		i3790 = (i1514 == 4 ? 1 : 0);
		v[3873] = (*a4)*i3790;
		i3789 = (i1514 == 6 ? 1 : 0);
		v[3874] = (*a4)*i3789;
		v[1879] = v[4950 + i1514];
		v[1881] = v[4986 + i1514];
		v[1883] = v[4968 + i1514];
		v[1885] = v[5022 + i1514];
		v[1887] = v[5058 + i1514];
		v[1889] = v[5040 + i1514];
		v[1891] = v[5094 + i1514];
		v[1893] = v[5130 + i1514];
		v[1895] = v[5112 + i1514];
		v[1897] = v[4932 + i1514];
		v[3016] = -8e0*v[1850] * v[1897];
		v[1946] = 2e0*v[1224] * v[1897];
		v[1966] = -4e0*v[146] * v[1946];
		v[3804] = v[1966] * v[285];
		v[1898] = v[5948 + i1514];
		v[1899] = v[5984 + i1514];
		v[1901] = v[5004 + i1514];
		v[3024] = -8e0*v[1852] * v[1901];
		v[1973] = 2e0*v[1238] * v[1901];
		v[1993] = -4e0*v[165] * v[1973];
		v[3805] = v[1993] * v[311];
		v[1902] = v[6056 + i1514];
		v[1903] = v[6092 + i1514];
		v[1905] = v[5076 + i1514];
		v[3032] = -8e0*v[1854] * v[1905];
		v[2054] = 2e0*v[1252] * v[1905];
		v[2074] = -4e0*v[2054] * v[225];
		v[3806] = v[2074] * v[337];
		v[1906] = v[6272 + i1514];
		v[1907] = v[6308 + i1514];
		v[1926] = -i3789 + v[1879];
		v[1927] = i3789 + v[1879];
		v[1928] = -i3790 + v[1881];
		v[1929] = i3790 + v[1881];
		v[1930] = i3791 + v[1883];
		v[1931] = -i3791 + v[1883];
		v[1932] = -i3792 + v[1885];
		v[1933] = i3792 + v[1885];
		v[1934] = -i3793 + v[1887];
		v[1935] = i3793 + v[1887];
		v[1936] = i3794 + v[1889];
		v[1937] = -i3794 + v[1889];
		v[1938] = -i3795 + v[1891];
		v[1939] = i3795 + v[1891];
		v[1940] = -i3796 + v[1893];
		v[1941] = i3796 + v[1893];
		v[1942] = i3797 + v[1895];
		v[1943] = -i3797 + v[1895];
		v[1944] = 2e0*dA[4] * i3791 - v[1897];
		v[1945] = -(i3789*v[146]) / 2e0 + dA[5] * v[1946];
		v[1947] = (i3791*v[146]) / 2e0 - dA[4] * v[1946];
		v[1948] = -(i3790*v[146]) / 2e0 + dA[3] * v[1946];
		v[1950] = -8e0*i3791*v[1224] + v[1897] * v[1949];
		v[1952] = 8e0*i3789*v[1224] + v[1897] * v[1951];
		v[1954] = 8e0*i3790*v[1224] + v[1897] * v[1953];
		v[1955] = (v[146] * v[1898]) / 2e0 - v[1946] * v[446];
		v[2000] = v[1955] * v[288];
		v[1956] = v[146] * v[1926] - 2e0*v[1946] * v[451];
		v[1957] = v[146] * v[1930] - 2e0*v[1946] * v[457];
		v[2001] = v[1957] * v[296];
		v[1958] = v[146] * v[1927] - 2e0*v[1946] * v[461];
		v[1959] = (v[146] * v[1944]) / 2e0 - v[1946] * v[465];
		v[2005] = v[1959] * v[294];
		v[1960] = v[146] * v[1928] - 2e0*v[1946] * v[470];
		v[2006] = v[1960] * v[296];
		v[1961] = v[146] * v[1931] - 2e0*v[1946] * v[474];
		v[2011] = v[1961] * v[288];
		v[1962] = v[146] * v[1929] - 2e0*v[1946] * v[478];
		v[1963] = (v[146] * v[1899]) / 2e0 - v[1946] * v[483];
		v[2009] = v[1963] * v[296];
		v[1964] = v[1966] + 2e0*v[1945] * v[456];
		v[1965] = -(v[1945] * v[450]) - v[1947] * v[456];
		v[1967] = v[1966] + 2e0*v[1947] * v[450];
		v[2127] = v[1544] * v[1967];
		v[1968] = v[1945] * v[453] + v[1948] * v[456];
		v[1969] = -(v[1948] * v[450]) - v[1947] * v[453];
		v[1970] = v[1966] + 2e0*v[1948] * v[453];
		v[2123] = v[1545] * v[1970];
		v[1971] = 2e0*dA[10] * i3794 - v[1901];
		v[1972] = -(i3792*v[165]) / 2e0 + dA[11] * v[1973];
		v[1974] = (i3794*v[165]) / 2e0 - dA[10] * v[1973];
		v[1975] = -(i3793*v[165]) / 2e0 + dA[9] * v[1973];
		v[1977] = -8e0*i3794*v[1238] + v[1901] * v[1976];
		v[1979] = 8e0*i3792*v[1238] + v[1901] * v[1978];
		v[1981] = 8e0*i3793*v[1238] + v[1901] * v[1980];
		v[1982] = (v[165] * v[1902]) / 2e0 - v[1973] * v[491];
		v[2015] = v[1982] * v[314];
		v[1983] = v[165] * v[1932] - 2e0*v[1973] * v[496];
		v[1984] = v[165] * v[1936] - 2e0*v[1973] * v[502];
		v[2016] = v[1984] * v[322];
		v[1985] = v[165] * v[1933] - 2e0*v[1973] * v[506];
		v[1986] = (v[165] * v[1971]) / 2e0 - v[1973] * v[510];
		v[2020] = v[1986] * v[320];
		v[1987] = v[165] * v[1934] - 2e0*v[1973] * v[515];
		v[2021] = v[1987] * v[322];
		v[1988] = v[165] * v[1937] - 2e0*v[1973] * v[519];
		v[2026] = v[1988] * v[314];
		v[1989] = v[165] * v[1935] - 2e0*v[1973] * v[523];
		v[1990] = (v[165] * v[1903]) / 2e0 - v[1973] * v[528];
		v[2024] = v[1990] * v[322];
		v[1991] = v[1993] + 2e0*v[1972] * v[501];
		v[1992] = -(v[1972] * v[495]) - v[1974] * v[501];
		v[1994] = v[1993] + 2e0*v[1974] * v[495];
		v[2155] = v[1541] * v[1994];
		v[1995] = v[1972] * v[498] + v[1975] * v[501];
		v[1996] = -(v[1975] * v[495]) - v[1974] * v[498];
		v[1997] = v[1993] + 2e0*v[1975] * v[498];
		v[2151] = v[1542] * v[1997];
		v[1998] = v[2000] + v[1956] * v[294];
		v[1999] = v[1998] + v[2001] + v[3873];
		v[2002] = v[2000] + v[2001];
		v[2003] = v[2005] + v[1958] * v[288];
		v[2004] = v[2003] + v[2006] + v[3872];
		v[2007] = v[2005] + v[2006];
		v[2008] = v[2009] + v[2011];
		v[2010] = v[2009] + v[1962] * v[294];
		v[2012] = v[2010] + v[2011] + v[3874];
		v[2013] = v[2015] + v[1983] * v[320];
		v[2014] = v[2013] + v[2016] + v[3865];
		v[2017] = v[2015] + v[2016];
		v[2018] = v[2020] + v[1985] * v[314];
		v[2019] = v[2018] + v[2021] + v[3864];
		v[2022] = v[2020] + v[2021];
		v[2023] = v[2024] + v[2026];
		v[2025] = v[2024] + v[1989] * v[320];
		v[2027] = v[2025] + v[2026] + v[3866];
		v[2028] = QAAi[0][0] * v[1961] + QAAi[1][0] * v[1962] + QAAi[2][0] * v[1963];
		v[2029] = QAAi[0][1] * v[1961] + QAAi[1][1] * v[1962] + QAAi[2][1] * v[1963];
		v[2030] = QBAi[0][0] * v[1988] + QBAi[1][0] * v[1989] + QBAi[2][0] * v[1990];
		v[2031] = QBAi[0][1] * v[1988] + QBAi[1][1] * v[1989] + QBAi[2][1] * v[1990];
		v[2032] = QAAi[0][0] * v[1958] + QAAi[1][0] * v[1959] + QAAi[2][0] * v[1960];
		v[2033] = QAAi[0][1] * v[1958] + QAAi[1][1] * v[1959] + QAAi[2][1] * v[1960];
		v[2034] = QBAi[0][0] * v[1985] + QBAi[1][0] * v[1986] + QBAi[2][0] * v[1987];
		v[2035] = QBAi[0][1] * v[1985] + QBAi[1][1] * v[1986] + QBAi[2][1] * v[1987];
		v[2036] = QAAi[0][0] * v[1955] + QAAi[1][0] * v[1956] + QAAi[2][0] * v[1957];
		v[2037] = QAAi[0][1] * v[1955] + QAAi[1][1] * v[1956] + QAAi[2][1] * v[1957];
		v[2038] = QBAi[0][0] * v[1982] + QBAi[1][0] * v[1983] + QBAi[2][0] * v[1984];
		v[2039] = QBAi[0][1] * v[1982] + QBAi[1][1] * v[1983] + QBAi[2][1] * v[1984];
		v[2040] = i3798 + v[2030] * v[214] + v[2031] * v[215];
		v[2041] = i3798;
		v[2042] = i3799 + v[2028] * v[214] + v[2029] * v[215];
		v[2043] = i3799;
		v[2044] = i3800 + v[2034] * v[214] + v[2035] * v[215];
		v[2045] = i3800;
		v[2046] = i3801 + v[2032] * v[214] + v[2033] * v[215];
		v[2047] = i3801;
		v[2048] = i3802 + v[2038] * v[214] + v[2039] * v[215];
		v[2049] = i3802;
		v[2050] = i3803 + v[2036] * v[214] + v[2037] * v[215];
		v[2051] = i3803;
		v[2052] = 2e0*dB[4] * i3797 - v[1905];
		v[2053] = dB[5] * v[2054] - (i3795*v[225]) / 2e0;
		v[2055] = -(dB[4] * v[2054]) + (i3797*v[225]) / 2e0;
		v[2056] = dB[3] * v[2054] - (i3796*v[225]) / 2e0;
		v[2058] = -8e0*i3797*v[1252] + v[1905] * v[2057];
		v[2060] = 8e0*i3795*v[1252] + v[1905] * v[2059];
		v[2062] = 8e0*i3796*v[1252] + v[1905] * v[2061];
		v[2063] = (v[1906] * v[225]) / 2e0 - v[2054] * v[662];
		v[2081] = v[2063] * v[340];
		v[2064] = v[1938] * v[225] - 2e0*v[2054] * v[667];
		v[2065] = v[1942] * v[225] - 2e0*v[2054] * v[673];
		v[2082] = v[2065] * v[348];
		v[2066] = v[1939] * v[225] - 2e0*v[2054] * v[677];
		v[2067] = (v[2052] * v[225]) / 2e0 - v[2054] * v[681];
		v[2089] = v[2067] * v[346];
		v[2068] = v[1940] * v[225] - 2e0*v[2054] * v[686];
		v[2090] = v[2068] * v[348];
		v[2069] = v[1943] * v[225] - 2e0*v[2054] * v[690];
		v[2098] = v[2069] * v[340];
		v[2070] = v[1941] * v[225] - 2e0*v[2054] * v[694];
		v[2071] = (v[1907] * v[225]) / 2e0 - v[2054] * v[699];
		v[2096] = v[2071] * v[348];
		v[2072] = v[2074] + 2e0*v[2053] * v[672];
		v[2073] = -(v[2053] * v[666]) - v[2055] * v[672];
		v[2075] = v[2074] + 2e0*v[2055] * v[666];
		v[2183] = v[1538] * v[2075];
		v[2076] = v[2053] * v[669] + v[2056] * v[672];
		v[2077] = -(v[2056] * v[666]) - v[2055] * v[669];
		v[2078] = v[2074] + 2e0*v[2056] * v[669];
		v[2179] = v[1539] * v[2078];
		v[2079] = v[2081] + v[2064] * v[346];
		v[2080] = v[2079] + v[2082] + v[3857];
		v[2083] = v[2081] + v[2082];
		v[2084] = QABi[0][0] * v[2063] + QABi[1][0] * v[2064] + QABi[2][0] * v[2065];
		v[2085] = QABi[0][1] * v[2063] + QABi[1][1] * v[2064] + QABi[2][1] * v[2065];
		v[2086] = QABi[0][2] * v[2063] + QABi[1][2] * v[2064] + QABi[2][2] * v[2065];
		v[2087] = v[2089] + v[2066] * v[340];
		v[2088] = v[2087] + v[2090] + v[3856];
		v[2091] = v[2089] + v[2090];
		v[2092] = QABi[0][0] * v[2066] + QABi[1][0] * v[2067] + QABi[2][0] * v[2068];
		v[2093] = QABi[0][1] * v[2066] + QABi[1][1] * v[2067] + QABi[2][1] * v[2068];
		v[2094] = QABi[0][2] * v[2066] + QABi[1][2] * v[2067] + QABi[2][2] * v[2068];
		v[2095] = v[2096] + v[2098];
		v[2097] = v[2096] + v[2070] * v[346];
		v[2099] = v[2097] + v[2098] + v[3858];
		v[2100] = QABi[0][0] * v[2069] + QABi[1][0] * v[2070] + QABi[2][0] * v[2071];
		v[2101] = QABi[0][1] * v[2069] + QABi[1][1] * v[2070] + QABi[2][1] * v[2071];
		v[2102] = QABi[0][2] * v[2069] + QABi[1][2] * v[2070] + QABi[2][2] * v[2071];
		v[2103] = v[1950] + v[1968];
		v[2121] = v[1545] * v[2103];
		v[2104] = -v[1950] + v[1968];
		v[2106] = (v[160] * v[2103] + v[1962] * v[283] + v[2105] * v[3804]) / v[285];
		v[2107] = v[1952] + v[1969];
		v[2129] = v[1545] * v[2107];
		v[2108] = -v[1952] + v[1969];
		v[2125] = v[1544] * v[2108];
		v[2110] = (v[156] * v[2107] + v[1960] * v[284] + v[2109] * v[3804]) / v[285];
		v[2112] = (v[151] * v[2108] + v[1957] * v[3634] + v[2111] * v[3804]) / v[285];
		v[2113] = v[2123] + v[2125];
		v[3971] = v[2113] / v[285];
		v[2114] = v[1954] + v[1965];
		v[2119] = v[1544] * v[2114];
		v[2115] = -v[1954] + v[1965];
		v[2116] = v[2127] + v[2129];
		v[3967] = v[2116] / v[285];
		v[2117] = v[1543] * v[1964] + v[2119] + v[2121];
		v[3975] = v[2117] / v[285];
		v[2120] = v[2117] - v[2121];
		v[2122] = v[2117] - v[2119];
		v[3968] = v[2122] / v[285];
		v[2124] = v[1543] * v[2104] + v[2123];
		v[2126] = v[2124] + v[2125];
		v[3969] = v[2126] / v[285];
		v[2128] = v[1543] * v[2115] + v[2127];
		v[2130] = v[2128] + v[2129];
		v[3972] = v[2130] / v[285];
		v[2131] = v[1977] + v[1995];
		v[2149] = v[1542] * v[2131];
		v[2132] = -v[1977] + v[1995];
		v[2134] = (v[179] * v[2131] + v[1989] * v[309] + v[2133] * v[3805]) / v[311];
		v[2135] = v[1979] + v[1996];
		v[2157] = v[1542] * v[2135];
		v[2136] = -v[1979] + v[1996];
		v[2153] = v[1541] * v[2136];
		v[2138] = (v[175] * v[2135] + v[1987] * v[310] + v[2137] * v[3805]) / v[311];
		v[2140] = (v[170] * v[2136] + v[1984] * v[3639] + v[2139] * v[3805]) / v[311];
		v[2141] = v[2151] + v[2153];
		v[3983] = v[2141] / v[311];
		v[2142] = v[1981] + v[1992];
		v[2147] = v[1541] * v[2142];
		v[2143] = -v[1981] + v[1992];
		v[2144] = v[2155] + v[2157];
		v[3979] = v[2144] / v[311];
		v[2145] = v[1540] * v[1991] + v[2147] + v[2149];
		v[3987] = v[2145] / v[311];
		v[2148] = v[2145] - v[2149];
		v[2150] = v[2145] - v[2147];
		v[3980] = v[2150] / v[311];
		v[2152] = v[1540] * v[2132] + v[2151];
		v[2154] = v[2152] + v[2153];
		v[3981] = v[2154] / v[311];
		v[2156] = v[1540] * v[2143] + v[2155];
		v[2158] = v[2156] + v[2157];
		v[3984] = v[2158] / v[311];
		v[2159] = v[2058] + v[2076];
		v[2177] = v[1539] * v[2159];
		v[2160] = -v[2058] + v[2076];
		v[2162] = (v[2159] * v[239] + v[2070] * v[335] + v[2161] * v[3806]) / v[337];
		v[2163] = v[2060] + v[2077];
		v[2185] = v[1539] * v[2163];
		v[2164] = -v[2060] + v[2077];
		v[2181] = v[1538] * v[2164];
		v[2166] = (v[2163] * v[235] + v[2068] * v[336] + v[2165] * v[3806]) / v[337];
		v[2168] = (v[2164] * v[230] + v[2065] * v[3644] + v[2167] * v[3806]) / v[337];
		v[2169] = v[2179] + v[2181];
		v[3996] = v[2169] / v[337];
		v[2170] = v[2062] + v[2073];
		v[2175] = v[1538] * v[2170];
		v[2171] = -v[2062] + v[2073];
		v[2172] = v[2183] + v[2185];
		v[3992] = v[2172] / v[337];
		v[2173] = v[1537] * v[2072] + v[2175] + v[2177];
		v[4000] = v[2173] / v[337];
		v[2176] = v[2173] - v[2177];
		v[2178] = v[2173] - v[2175];
		v[3993] = v[2178] / v[337];
		v[2180] = v[1537] * v[2160] + v[2179];
		v[2182] = v[2180] + v[2181];
		v[3994] = v[2182] / v[337];
		v[2184] = v[1537] * v[2171] + v[2183];
		v[2186] = v[2184] + v[2185];
		v[3997] = v[2186] / v[337];
		v[2187] = -i3815 + v[2050] * v[210] + v[2048] * v[211] - v[2084] * v[260] - v[2085] * v[261] - v[2086] * v[263];
		v[2189] = -i3814 + v[2046] * v[210] + v[2044] * v[211] - v[2092] * v[260] - v[2093] * v[261] - v[2094] * v[263];
		v[2191] = -i3813 + v[2042] * v[210] + v[2040] * v[211] - v[2100] * v[260] - v[2101] * v[261] - v[2102] * v[263];
		v[3854] = v[2187] * v[360] + v[2189] * v[361] + v[2191] * v[362];
		v[2193] = v[3854] / v[1644];
		v[3807] = v[2193] * v[373];
		v[2203] = v[2191] * v[374] + v[362] * v[3807];
		v[2199] = v[2189] * v[374] + v[361] * v[3807];
		v[2195] = v[2187] * v[374] + v[360] * v[3807];
		v[2196] = 2e0*v[1575] * v[2195];
		v[2197] = v[2195];
		v[2200] = 2e0*v[1567] * v[2199];
		v[2201] = v[2199];
		v[2204] = 2e0*v[1566] * v[2203];
		v[2205] = v[2203];
		v[2206] = 0e0;
		v[2207] = 0e0;
		v[2208] = 0e0;
		v[2209] = 0e0;
		v[2210] = 0e0;
		v[2211] = 0e0;
		v[2212] = 0e0;
		v[2213] = 0e0;
		v[2214] = 0e0;
		v[2215] = 0e0;
		v[2216] = 0e0;
		v[2217] = 0e0;
		v[2218] = 0e0;
		v[2219] = 0e0;
		v[2220] = 0e0;
		v[2221] = 0e0;
		v[2222] = 0e0;
		v[2223] = 0e0;
		v[2224] = 0e0;
		v[2225] = 0e0;
		v[2226] = 0e0;
		v[2227] = 0e0;
		v[2228] = 0e0;
		v[2229] = 0e0;
		v[2230] = 0e0;
		v[2231] = 0e0;
		v[2232] = 0e0;
		v[2233] = 0e0;
		v[2234] = 0e0;
		v[2235] = 0e0;
		v[2236] = 0e0;
		v[2237] = 0e0;
		b2238 = b909;
		if (b2238) {
			v[2241] = v[2205];
			v[2239] = v[2201];
			v[2240] = v[2205] * v[385] - v[2239] * v[386];
			v[2242] = -(v[2241] * v[384]) + v[2197] * v[386];
			v[2243] = v[2239] * v[384] - v[2197] * v[385];
			v[3811] = v[1635] * v[2240] + v[1630] * v[2242] + v[1625] * v[2243];
			v[3809] = v[2240] * v[911] + v[2242] * v[912] + v[2243] * v[913];
			v[2244] = v[3809] / v[914];
			v[3810] = v[2244] * v[4116];
			v[2254] = v[2244] * v[3808];
			v[2209] = -(v[1638] * v[2245] * v[3809]);
			v[2246] = v[2240] * v[3696] + v[3810] * v[911];
			v[2265] = 2e0*v[2246] * v[923];
			v[2248] = v[2242] * v[3696] + v[3810] * v[912];
			v[2262] = 2e0*v[2248] * v[924];
			v[2249] = v[2243] * v[3696] + v[3810] * v[913];
			v[2263] = 2e0*v[2249] * v[925];
			v[2212] = v[2254] * v[3741] + v[3811] * v[916];
			v[2211] = v[2244] * v[2252] * v[920];
			v[2210] = v[2244] * v[3741] * v[921] + v[3811] * v[922];
			v[2236] = 2e0*v[2254] * v[2255] * v[4117];
			v[2237] = v[1637] * v[2244] * v[2250] * v[2255] * v[4118];
			v[2208] = v[1625] * v[3810] + v[2243] * v[3812];
			v[2207] = v[1630] * v[3810] + v[2242] * v[3812];
			v[2206] = v[1635] * v[3810] + v[2240] * v[3812];
			v[2257] = (v[2248] * v[923] + v[2246] * v[924]) / 2e0;
			v[2258] = v[2262] + v[2265];
			v[2259] = v[2258] + v[2263];
			v[2260] = (v[2249] * v[923] + v[2246] * v[925]) / 2e0;
			v[2261] = (v[2249] * v[924] + v[2248] * v[925]) / 2e0;
			v[2264] = v[2262] + v[2263];
			v[2266] = v[2263] + v[2265];
			v[2215] = (v[1610] * v[2246] + v[1609] * v[2248] + 4e0*v[2249] * v[2267]) / 2e0;
			v[2214] = (v[1611] * v[2246] + v[1609] * v[2249] + 4e0*v[2248] * v[2268]) / 2e0;
			v[2213] = (v[1611] * v[2248] + v[1610] * v[2249] + 4e0*v[2246] * v[2269]) / 2e0;
			v[2270] = -4e0*v[2259] * v[2579];
			v[2235] = 8e0*v[1621] * v[2259] * v[4119];
			v[2226] = (v[1556] * v[2270]) / 2e0;
			v[2230] = (v[1560] * v[2270]) / 2e0;
			v[2228] = v[1558] * v[2270];
			v[2232] = v[1562] * v[2270];
			v[2231] = v[1561] * v[2270];
			v[2233] = v[1563] * v[2270];
			v[2227] = v[1557] * v[2270];
			v[2229] = v[1559] * v[2270];
			v[2234] = (v[1564] * v[2270]) / 2e0;
			v[2271] = -v[2249] + v[2257];
			v[2272] = v[2249] + v[2257];
			v[2273] = v[2248] + v[2260];
			v[2274] = -v[2248] + v[2260];
			v[2275] = -v[2246] + v[2261];
			v[2276] = v[2246] + v[2261];
			v[2218] = v[1613] * v[2270] + v[2271] * v[926];
			v[2220] = v[1615] * v[2270] + v[2272] * v[926];
			v[2219] = v[1614] * v[2270] + v[2273] * v[926];
			v[2223] = v[1618] * v[2270] + v[2274] * v[926];
			v[2217] = (v[1612] * v[2270] - v[2264] * v[926]) / 2e0;
			v[2222] = v[1617] * v[2270] + v[2275] * v[926];
			v[2224] = v[1619] * v[2270] + v[2276] * v[926];
			v[2221] = (v[1616] * v[2270] - v[2266] * v[926]) / 2e0;
			v[2225] = (v[1620] * v[2270] - v[2258] * v[926]) / 2e0;
			v[2216] = -(v[1564] * v[2258]) / 2e0 - (v[1556] * v[2264]) / 2e0 - (v[1560] * v[2266]) / 2e0 + v[1557] * v[2271]
				+ v[1559] * v[2272] + v[1558] * v[2273] + v[1562] * v[2274] + v[1561] * v[2275] + v[1563] * v[2276];
		}
		else {
		};
		v[2277] = 0e0;
		v[2278] = 0e0;
		v[2279] = 0e0;
		v[2280] = 0e0;
		v[2281] = 0e0;
		v[2282] = 0e0;
		b2283 = (*previouscontact);
		if (b2283) {
			v[2301] = v[2205];
			v[2295] = v[2201];
			v[2293] = v[2197];
			v[2291] = v[2218];
			v[2290] = v[2219];
			v[2288] = v[2221];
			v[2287] = v[2222];
			v[2285] = v[2224];
			v[2284] = v[2225];
			v[2299] = v[1552] * v[2293];
			v[2298] = v[1550] * v[2205];
			v[2296] = v[1551] * v[2201];
			v[2043] = v[2043] + v[2028] * v[212] + v[2029] * v[213];
			v[2041] = v[2041] + v[2030] * v[212] + v[2031] * v[213];
			v[2225] = 0e0;
			v[2224] = 0e0;
			v[2286] = -i3813 + v[2043] * v[208] + v[2041] * v[209] + gti[0] * v[2223] + gti[2] * v[2284] + gti[1] * v[2285]
				- v[2100] * v[256] - v[2101] * v[257] - v[2102] * v[259];
			v[3816] = -(v[2286] * v[377]) - v[2301] * v[946];
			v[2223] = 0e0;
			v[2047] = v[2047] + v[2032] * v[212] + v[2033] * v[213];
			v[2045] = v[2045] + v[2034] * v[212] + v[2035] * v[213];
			v[2222] = 0e0;
			v[2221] = 0e0;
			v[2289] = -i3814 + v[2047] * v[208] + v[2045] * v[209] + gti[0] * v[2220] + gti[2] * v[2287] + gti[1] * v[2288]
				- v[2092] * v[256] - v[2093] * v[257] - v[2094] * v[259];
			v[3818] = -(v[2289] * v[376]) - v[2295] * v[945];
			v[2220] = 0e0;
			v[2051] = v[2051] + v[2036] * v[212] + v[2037] * v[213];
			v[2049] = v[2049] + v[2038] * v[212] + v[2039] * v[213];
			v[2219] = 0e0;
			v[2218] = 0e0;
			v[2292] = -i3815 + v[2051] * v[208] + v[2049] * v[209] + gti[0] * v[2217] + gti[2] * v[2290] + gti[1] * v[2291]
				- v[2084] * v[256] - v[2085] * v[257] - v[2086] * v[259];
			v[3817] = -(v[2292] * v[375]) - v[2293] * v[944];
			v[2217] = 0e0;
			v[2197] = 0e0;
			v[2294] = v[2296] + v[2299];
			v[2201] = 0e0;
			v[2297] = v[2296] + v[2298];
			v[2300] = v[2298] + v[2299];
			v[2205] = 0e0;
			v[2282] = v[1550] * v[2286];
			v[2281] = v[1551] * v[2289];
			v[2280] = v[1552] * v[2292];
			v[2277] = v[2293] * v[2302] - (2e0*v[1552] * v[2195] + v[2297])*v[375];
			v[2196] = v[2196] + v[2292] * v[2302] + v[1552] * (v[3816] + v[3818]) - v[2297] * v[944];
			v[2278] = v[2295] * v[2304] - (2e0*v[1551] * v[2199] + v[2300])*v[376];
			v[2200] = v[2200] + v[2289] * v[2304] + v[1551] * (v[3816] + v[3817]) - v[2300] * v[945];
			v[2279] = v[2301] * v[2307] - (2e0*v[1550] * v[2203] + v[2294])*v[377];
			v[2204] = v[2204] + v[2286] * v[2307] + v[1550] * (v[3817] + v[3818]) - v[2294] * v[946];
		}
		else {
		};
		v[2312] = v[1966] * v[2308] + v[1970] * v[2309] + v[2103] * v[2310] + v[2107] * v[2311] + v[2106] * v[294]
			+ v[2110] * v[296] + v[1999] * v[3037] + v[2003] * v[3635] + v[2008] * v[3637] + (*a4)*v[6596 + i1514];
		v[2317] = v[1859] * v[1998] + v[1858] * v[2010] + v[1966] * v[2313] + v[1967] * v[2314] + v[2108] * v[2315]
			+ v[2114] * v[2316] + v[2106] * v[288] + v[2112] * v[296] + v[2004] * v[3038] + (*a4)*v[6614 + i1514];
		v[2322] = v[1861] * v[2002] + v[1860] * v[2007] + v[1966] * v[2318] + v[1964] * v[2319] + v[2104] * v[2320]
			+ v[2115] * v[2321] + v[2110] * v[288] + v[2112] * v[294] + v[2012] * v[3039] + (*a4)*v[6632 + i1514];
		v[2327] = v[1993] * v[2323] + v[1997] * v[2324] + v[2131] * v[2325] + v[2135] * v[2326] + v[2014] * v[3040]
			+ v[2134] * v[320] + v[2138] * v[322] + v[2018] * v[3640] + v[2023] * v[3642] + (*a4)*v[6650 + i1514];
		v[2332] = v[1865] * v[2013] + v[1864] * v[2025] + v[1993] * v[2328] + v[1994] * v[2329] + v[2136] * v[2330]
			+ v[2142] * v[2331] + v[2019] * v[3041] + v[2134] * v[314] + v[2140] * v[322] + (*a4)*v[6668 + i1514];
		v[2337] = v[1867] * v[2017] + v[1866] * v[2022] + v[1993] * v[2333] + v[1991] * v[2334] + v[2132] * v[2335]
			+ v[2143] * v[2336] + v[2027] * v[3042] + v[2138] * v[314] + v[2140] * v[320] + (*a4)*v[6686 + i1514];
		v[2342] = v[2074] * v[2338] + v[2078] * v[2339] + v[2159] * v[2340] + v[2163] * v[2341] + v[2080] * v[3043]
			+ v[2162] * v[346] + v[2166] * v[348] + v[2087] * v[3645] + v[2095] * v[3647] + (*a4)*v[6704 + i1514];
		v[2347] = v[1871] * v[2079] + v[1870] * v[2097] + v[2074] * v[2343] + v[2075] * v[2344] + v[2164] * v[2345]
			+ v[2170] * v[2346] + v[2088] * v[3044] + v[2162] * v[340] + v[2168] * v[348] + (*a4)*v[6722 + i1514];
		v[2352] = v[1873] * v[2083] + v[1872] * v[2091] + v[2074] * v[2348] + v[2072] * v[2349] + v[2160] * v[2350]
			+ v[2171] * v[2351] + v[2099] * v[3045] + v[2166] * v[340] + v[2168] * v[346] + (*a4)*v[6740 + i1514];
		v[2353] = (*ct)*(v[210] * (v[1439] * (v[1518] * v[2029] + v[1517] * v[2033] + v[1516] * v[2037]) - (v[1518] * v[2028]
			+ v[1517] * v[2032] + v[1516] * v[2036])*v[993]) + v[211] * (v[1439] * (v[1518] * v[2031] + v[1517] * v[2035]
				+ v[1516] * v[2039]) - (v[1518] * v[2030] + v[1517] * v[2034] + v[1516] * v[2038])*v[993]));
		v[2354] = ((*ct)*(v[1518] * (v[2040] - v[2042]) + v[1517] * (v[2044] - v[2046]) + v[1516] * (v[2048] - v[2050])))
			/ 2e0;
		v[2355] = (*ct)*(v[1516] * (v[1440] * v[2085] + v[2084] * v[418] + v[2086] * v[420]) + v[1517] * (v[1440] * v[2093]
			+ v[2092] * v[418] + v[2094] * v[420]) + v[1518] * (v[1440] * v[2101] + v[2100] * v[418] + v[2102] * v[420]));
		v[2356] = (*ct)*((v[1516] * v[2084] + v[1517] * v[2092] + v[1518] * v[2100])*v[410] + (v[1516] * v[2086]
			+ v[1517] * v[2094] + v[1518] * v[2102])*v[412]);
		v[2357] = -(i3803*v[1441]) - i3801 * v[1445] - i3799 * v[1449] - i3790 * v[1453] - i3791 * v[1457] - i3789 * v[1461]
			- i3802 * v[1465] - i3800 * v[1469] - i3798 * v[1473] - i3793 * v[1477] - i3794 * v[1481] - i3792 * v[1485] - i3815 * v[1489]
			- i3814 * v[1493] - i3813 * v[1497] - i3796 * v[1501] - i3797 * v[1505] - i3795 * v[1509];
		v[2372] = v[2357];
		v[2358] = -(i3803*v[1442]) - i3801 * v[1446] - i3799 * v[1450] - i3790 * v[1454] - i3791 * v[1458] - i3789 * v[1462]
			- i3802 * v[1466] - i3800 * v[1470] - i3798 * v[1474] - i3793 * v[1478] - i3794 * v[1482] - i3792 * v[1486] - i3815 * v[1490]
			- i3814 * v[1494] - i3813 * v[1498] - i3796 * v[1502] - i3797 * v[1506] - i3795 * v[1510];
		v[2370] = v[2358];
		v[2359] = -(i3803*v[1443]) - i3801 * v[1447] - i3799 * v[1451] - i3790 * v[1455] - i3791 * v[1459] - i3789 * v[1463]
			- i3802 * v[1467] - i3800 * v[1471] - i3798 * v[1475] - i3793 * v[1479] - i3794 * v[1483] - i3792 * v[1487] - i3815 * v[1491]
			- i3814 * v[1495] - i3813 * v[1499] - i3796 * v[1503] - i3797 * v[1507] - i3795 * v[1511];
		v[2368] = v[2359];
		v[2360] = 0e0;
		v[2361] = 0e0;
		v[2362] = 0e0;
		v[2363] = 0e0;
		v[2364] = 0e0;
		v[2365] = 0e0;
		b2366 = (*stick);
		if (b2366) {
			b2367 = b1017;
			if (b2367) {
				v[2365] = v[2359];
				v[2359] = 0e0;
				v[2364] = v[2358];
				v[2358] = 0e0;
				v[2363] = v[2357];
				v[2357] = 0e0;
			}
			else {
				v[3819] = (*mud)*(v[1015] * v[2368] + v[1014] * v[2370] + v[1013] * v[2372]);
				v[2359] = 0e0;
				v[2358] = 0e0;
				v[3821] = v[1028] * v[2383] * v[3819];
				v[2357] = 0e0;
				v[3820] = v[1027] * v[1032] * v[2381] * v[3819];
				v[2365] = v[2368] * v[2380] + v[1015] * v[3820];
				v[2364] = v[2370] * v[2380] + v[1014] * v[3820];
				v[2363] = v[2372] * v[2380] + v[1013] * v[3820];
				v[2362] = v[1006] * v[3821];
				v[2361] = v[1005] * v[3821];
				v[2360] = v[1004] * v[3821];
			};
		}
		else {
			b2384 = b1051;
			if (b2384) {
				v[2365] = v[2368];
				v[2359] = 0e0;
				v[2364] = v[2370];
				v[2358] = 0e0;
				v[2363] = v[2372];
				v[2357] = 0e0;
			}
			else {
				v[2390] = v[2372] * v[3822];
				v[2388] = v[2370] * v[3822];
				v[2387] = v[2368] * v[3822];
				v[3824] = (*mud)*v[1058] * (v[1015] * v[2368] + v[1014] * v[2370] + v[1013] * v[2372])*v[2393];
				v[3823] = v[1057] * v[2389] * (v[1015] * v[2387] + v[1014] * v[2388] + v[1013] * v[2390]);
				v[2365] = v[1058] * v[2387] + v[1015] * v[3823];
				v[2364] = v[1058] * v[2388] + v[1014] * v[3823];
				v[2363] = v[1058] * v[2390] + v[1013] * v[3823];
				v[2362] = v[1006] * v[3824];
				v[2361] = v[1005] * v[3824];
				v[2360] = v[1004] * v[3824];
			};
		};
		v[3835] = v[2360] * v[376];
		v[3833] = (*cn)*v[2360];
		v[2518] = -(v[2360] * v[3651]);
		v[2503] = -(v[3833] * v[982]);
		v[3825] = (*cn)*v[2361];
		v[2513] = v[376] * v[3825];
		v[2509] = -(v[3825] * v[986]);
		v[3834] = v[2361] * v[376] + v[2362] * v[377];
		v[3826] = (*cn)*v[2362];
		v[2531] = v[377] * v[3826];
		v[3988] = (v[2513] + v[2531])*v[375];
		v[3963] = -v[2503] + v[2531] * v[375];
		v[3839] = -v[2503] + v[3988];
		v[2515] = -(v[3826] * v[990]);
		v[3836] = v[2509] - v[2531] * v[376];
		v[3838] = -(v[2518] * v[376]) + v[3836];
		v[3831] = -((*ct)*v[2363]);
		v[3829] = -((*ct)*v[2364]);
		v[3827] = -((*ct)*v[2365]);
		v[2395] = -(i3795*v[1021]) - (*ct)*(v[1518] * v[2352] + v[2365] * v[359]);
		v[2396] = -(i3797*v[1021]) - (*ct)*(v[1518] * v[2347] + v[2365] * v[349]);
		v[2397] = -(i3796*v[1021]) - (*ct)*(v[1518] * v[2342] + v[2365] * v[339]);
		v[2398] = v[281] * v[3827] + i3813 * v[3828];
		v[2400] = v[280] * v[3827] + i3814 * v[3828];
		v[2401] = v[279] * v[3827] + i3815 * v[3828];
		v[2402] = -(i3792*v[1021]) - (*ct)*(v[1518] * v[2337] + v[2365] * v[333]);
		v[2403] = -(i3794*v[1021]) - (*ct)*(v[1518] * v[2332] + v[2365] * v[323]);
		v[2404] = -(i3793*v[1021]) - (*ct)*(v[1518] * v[2327] + v[2365] * v[313]);
		v[2406] = v[278] * v[3827] + i3798 * v[3828];
		v[2408] = v[277] * v[3827] + i3800 * v[3828];
		v[2410] = v[276] * v[3827] + i3802 * v[3828];
		v[2411] = -(i3789*v[1021]) - (*ct)*(v[1518] * v[2322] + v[2365] * v[307]);
		v[2412] = -(i3791*v[1021]) - (*ct)*(v[1518] * v[2317] + v[2365] * v[297]);
		v[2413] = -(i3790*v[1021]) - (*ct)*(v[1518] * v[2312] + v[2365] * v[287]);
		v[2415] = v[275] * v[3827] + i3799 * v[3828];
		v[2417] = v[274] * v[3827] + i3801 * v[3828];
		v[2419] = v[273] * v[3827] + i3803 * v[3828];
		v[2420] = -(i3795*v[1020]) - (*ct)*(v[1517] * v[2352] + v[2364] * v[359]);
		v[2421] = -(i3797*v[1020]) - (*ct)*(v[1517] * v[2347] + v[2364] * v[349]);
		v[2422] = -(i3796*v[1020]) - (*ct)*(v[1517] * v[2342] + v[2364] * v[339]);
		v[2423] = v[281] * v[3829] + i3813 * v[3830];
		v[2425] = v[280] * v[3829] + i3814 * v[3830];
		v[2426] = v[279] * v[3829] + i3815 * v[3830];
		v[2427] = -(i3792*v[1020]) - (*ct)*(v[1517] * v[2337] + v[2364] * v[333]);
		v[2428] = -(i3794*v[1020]) - (*ct)*(v[1517] * v[2332] + v[2364] * v[323]);
		v[2429] = -(i3793*v[1020]) - (*ct)*(v[1517] * v[2327] + v[2364] * v[313]);
		v[2430] = v[278] * v[3829] + i3798 * v[3830];
		v[2431] = v[277] * v[3829] + i3800 * v[3830];
		v[2432] = v[276] * v[3829] + i3802 * v[3830];
		v[2433] = -(i3789*v[1020]) - (*ct)*(v[1517] * v[2322] + v[2364] * v[307]);
		v[2434] = -(i3791*v[1020]) - (*ct)*(v[1517] * v[2317] + v[2364] * v[297]);
		v[2435] = -(i3790*v[1020]) - (*ct)*(v[1517] * v[2312] + v[2364] * v[287]);
		v[2436] = v[275] * v[3829] + i3799 * v[3830];
		v[2437] = v[274] * v[3829] + i3801 * v[3830];
		v[2438] = v[273] * v[3829] + i3803 * v[3830];
		v[2439] = -(i3795*v[1019]) - (*ct)*(v[1516] * v[2352] + v[2363] * v[359]);
		v[2440] = -(i3797*v[1019]) - (*ct)*(v[1516] * v[2347] + v[2363] * v[349]);
		v[2441] = -(i3796*v[1019]) - (*ct)*(v[1516] * v[2342] + v[2363] * v[339]);
		v[2442] = v[281] * v[3831] + i3813 * v[3832];
		v[2444] = v[280] * v[3831] + i3814 * v[3832];
		v[2445] = v[279] * v[3831] + i3815 * v[3832];
		v[2446] = -(i3792*v[1019]) - (*ct)*(v[1516] * v[2337] + v[2363] * v[333]);
		v[2447] = -(i3794*v[1019]) - (*ct)*(v[1516] * v[2332] + v[2363] * v[323]);
		v[2448] = -(i3793*v[1019]) - (*ct)*(v[1516] * v[2327] + v[2363] * v[313]);
		v[2449] = v[278] * v[3831] + i3798 * v[3832];
		v[2450] = v[277] * v[3831] + i3800 * v[3832];
		v[2451] = v[276] * v[3831] + i3802 * v[3832];
		v[2452] = -(i3789*v[1019]) - (*ct)*(v[1516] * v[2322] + v[2363] * v[307]);
		v[2453] = -(i3791*v[1019]) - (*ct)*(v[1516] * v[2317] + v[2363] * v[297]);
		v[2454] = -(i3790*v[1019]) - (*ct)*(v[1516] * v[2312] + v[2363] * v[287]);
		v[2455] = v[275] * v[3831] + i3799 * v[3832];
		v[2456] = v[274] * v[3831] + i3801 * v[3832];
		v[2457] = v[273] * v[3831] + i3803 * v[3832];
		v[2461] = v[3287] * v[3826];
		v[2462] = v[3290] * v[3825];
		v[2463] = v[2513] + v[2518];
		v[3892] = v[2463] * v[377];
		v[3837] = -v[2515] + v[3892];
		v[2204] = v[2204] + (*cn)*(v[2360] * v[3284] + v[2361] * v[3285] + v[2362] * v[3286]);
		v[2468] = (*cn)*(v[2361] * v[273] + v[2360] * v[274]);
		v[2469] = (*cn)*(v[2361] * v[276] + v[2360] * v[277]);
		v[3853] = v[210] * v[2468] + v[211] * v[2469];
		v[2200] = v[2200] + (*cn)*(v[2360] * v[3289] + v[2361] * v[3291] + v[2362] * v[3292]);
		v[2480] = v[3294] * v[3833];
		v[2196] = v[2196] + (*cn)*(v[2360] * v[3295] + v[2361] * v[3296] + v[2362] * v[3297]);
		v[2502] = v[2503] * v[339] + v[2501] * v[3834];
		v[2505] = v[2503] * v[349] + v[2504] * v[3834];
		v[2507] = v[2503] * v[359] + v[2506] * v[3834];
		v[2508] = v[2501] * v[3835] + v[339] * v[3836];
		v[2511] = v[2504] * v[3835] + v[349] * v[3836];
		v[2512] = v[2506] * v[3835] + v[359] * v[3836];
		v[2514] = v[2515] * v[339] + (v[2360] * v[2501] - v[2513] * v[339])*v[377];
		v[2516] = v[2515] * v[349] + (v[2360] * v[2504] - v[2513] * v[349])*v[377];
		v[2517] = v[2515] * v[359] + (v[2360] * v[2506] - v[2513] * v[359])*v[377];
		v[2519] = v[313] * v[3837];
		v[2520] = v[323] * v[3837];
		v[2521] = v[333] * v[3837];
		v[2522] = v[287] * v[3837];
		v[2523] = v[297] * v[3837];
		v[2524] = v[307] * v[3837];
		v[2525] = -(v[313] * v[3838]);
		v[2526] = -(v[323] * v[3838]);
		v[2527] = -(v[333] * v[3838]);
		v[2528] = -(v[287] * v[3838]);
		v[2529] = -(v[297] * v[3838]);
		v[2530] = -(v[307] * v[3838]);
		v[2532] = v[313] * v[3839];
		v[2533] = v[323] * v[3839];
		v[2534] = v[333] * v[3839];
		v[2535] = v[287] * v[3839];
		v[2536] = v[297] * v[3839];
		v[2537] = v[307] * v[3839];
		v[2538] = 0e0;
		v[2539] = 0e0;
		v[2540] = 0e0;
		v[2541] = 0e0;
		v[2542] = 0e0;
		v[2543] = 0e0;
		v[2544] = 0e0;
		v[2545] = 0e0;
		v[2546] = 0e0;
		b2547 = (*previouscontact);
		if (b2547) {
			v[2460] = (*epst)*v[2363];
			v[3842] = v[2460] * v[375];
			v[2459] = (*epst)*v[2364];
			v[3841] = v[2459] * v[376];
			v[3843] = v[3841] + v[3842];
			v[2458] = (*epst)*v[2365];
			v[3840] = v[2458] * v[377];
			v[3845] = v[3840] + v[3841];
			v[3844] = v[3840] + v[3842];
			v[2282] = v[2282] + v[2458] * v[946];
			v[2281] = v[2281] + v[2459] * v[945];
			v[2277] = v[2277] + v[1568] * v[2460] - v[375] * v[3845];
			v[2278] = v[2278] + v[1570] * v[2459] - v[376] * v[3844];
			v[2279] = v[2279] + v[1572] * v[2458] - v[377] * v[3843];
			v[2280] = v[2280] + v[2460] * v[944];
			v[2204] = v[2204] + v[2458] * v[3738] - v[3843] * v[946];
			v[2200] = v[2200] + v[2459] * v[3739] - v[3844] * v[945];
			v[2196] = v[2196] - v[2460] * v[3740] - v[3845] * v[944];
			v[2538] = gti[0] * v[2277];
			v[2539] = gti[1] * v[2277];
			v[2540] = gti[2] * v[2277];
			v[2548] = -v[2277];
			v[2549] = -(v[2277] * v[259]);
			v[2550] = -(v[2277] * v[257]);
			v[2551] = -(v[2277] * v[256]);
			v[2552] = v[209] * v[2277];
			v[2553] = v[208] * v[2277];
			v[2554] = v[213] * v[2552];
			v[2555] = v[212] * v[2552];
			v[2556] = v[213] * v[2553];
			v[2557] = v[212] * v[2553];
			v[2541] = gti[0] * v[2278];
			v[2542] = gti[1] * v[2278];
			v[2543] = gti[2] * v[2278];
			v[2558] = -v[2278];
			v[2559] = -(v[2278] * v[259]);
			v[2560] = -(v[2278] * v[257]);
			v[2561] = -(v[2278] * v[256]);
			v[2562] = v[209] * v[2278];
			v[2563] = v[208] * v[2278];
			v[2564] = v[213] * v[2562];
			v[2565] = v[212] * v[2562];
			v[2566] = v[213] * v[2563];
			v[2567] = v[212] * v[2563];
			v[2544] = gti[0] * v[2279];
			v[2545] = gti[1] * v[2279];
			v[2546] = gti[2] * v[2279];
			v[2568] = -v[2279];
			v[2569] = -(v[2279] * v[259]);
			v[2570] = -(v[2279] * v[257]);
			v[2571] = -(v[2279] * v[256]);
			v[2572] = v[209] * v[2279];
			v[2573] = v[208] * v[2279];
			v[2574] = v[213] * v[2572];
			v[2575] = v[212] * v[2572];
			v[2576] = v[213] * v[2573];
			v[2577] = v[212] * v[2573];
			v[2480] = -v[2280] + v[2480];
			v[2462] = -v[2281] + v[2462];
			v[2461] = -v[2282] + v[2461];
		}
		else {
			v[2557] = 0e0;
			v[2556] = 0e0;
			v[2567] = 0e0;
			v[2566] = 0e0;
			v[2577] = 0e0;
			v[2576] = 0e0;
			v[2555] = 0e0;
			v[2554] = 0e0;
			v[2565] = 0e0;
			v[2564] = 0e0;
			v[2575] = 0e0;
			v[2574] = 0e0;
			v[2553] = 0e0;
			v[2563] = 0e0;
			v[2573] = 0e0;
			v[2552] = 0e0;
			v[2562] = 0e0;
			v[2572] = 0e0;
			v[2551] = 0e0;
			v[2550] = 0e0;
			v[2549] = 0e0;
			v[2561] = 0e0;
			v[2560] = 0e0;
			v[2559] = 0e0;
			v[2571] = 0e0;
			v[2570] = 0e0;
			v[2569] = 0e0;
			v[2548] = 0e0;
			v[2558] = 0e0;
			v[2568] = 0e0;
		};
		b2578 = b909;
		if (b2578) {
			v[2234] = v[2234] + (v[2546] * v[926]) / 2e0;
			v[2233] = v[2233] + v[2545] * v[926];
			v[2232] = v[2232] + v[2544] * v[926];
			v[2231] = v[2231] + v[2543] * v[926];
			v[2230] = v[2230] + (v[2542] * v[926]) / 2e0;
			v[2229] = v[2229] + v[2541] * v[926];
			v[2228] = v[2228] + v[2540] * v[926];
			v[2227] = v[2227] + v[2539] * v[926];
			v[2216] = v[2216] + (v[1612] * v[2538]) / 2e0 + v[1613] * v[2539] + v[1614] * v[2540] + v[1615] * v[2541] +
				(v[1616] * v[2542]) / 2e0 + v[1617] * v[2543] + v[1618] * v[2544] + v[1619] * v[2545] + (v[1620] * v[2546]) / 2e0;
			v[2226] = v[2226] + (v[2538] * v[926]) / 2e0;
			v[2235] = v[2235] - 4e0*v[2216] * v[2579];
			v[3846] = v[2226] - v[2235];
			v[3848] = (v[2228] + v[2232]) / 2e0;
			v[3847] = (v[2231] + v[2233]) / 2e0;
			v[2215] = v[2215] - v[2227] + v[2229] + v[3848] * v[923] + v[3847] * v[924] - 2e0*(v[2230] + v[3846])*v[925];
			v[3849] = (v[2227] + v[2229]) / 2e0;
			v[2214] = v[2214] + v[2228] - v[2232] + v[3849] * v[923] - 2e0*(v[2234] + v[3846])*v[924] + v[3847] * v[925];
			v[2213] = v[2213] - v[2231] + v[2233] - 2e0*(v[2230] + v[2234] - v[2235])*v[923] + v[3849] * v[924]
				+ v[3848] * v[925];
			v[3850] = v[2213] * v[911] + v[2214] * v[912] + v[2215] * v[913];
			v[2212] = v[2212] + v[3850] * v[916];
			v[2210] = v[2210] + v[3850] * v[922];
			v[2211] = v[2211] + v[2212] * v[921];
			v[2236] = v[2236] + 2e0*v[2210] * v[2250];
			v[2237] = v[2237] + (v[2236] * v[2251]) / 2e0;
			v[2209] = v[2209] + v[2211] + v[2237];
			v[3851] = v[2209] / v[914];
			v[2208] = v[2208] + v[2215] * v[3696] + v[3851] * v[913];
			v[2207] = v[2207] + v[2214] * v[3696] + v[3851] * v[912];
			v[2206] = v[2206] + v[2213] * v[3696] + v[3851] * v[911];
			v[2196] = v[2196] - v[2208] * v[385] + v[2207] * v[386];
			v[2204] = v[2204] - v[2207] * v[384] + v[2206] * v[385];
			v[2200] = v[2200] + v[2208] * v[384] - v[2206] * v[386];
		}
		else {
		};
		v[2586] = -(v[2457] * v[873]) - v[2456] * v[874] - v[2455] * v[875] - v[2454] * v[876] - v[2453] * v[877]
			- v[2452] * v[878] - v[2451] * v[879] - v[2450] * v[880] - v[2449] * v[881] - v[2448] * v[882] - v[2447] * v[883]
			- v[2446] * v[884] - v[2445] * v[885] - v[2444] * v[886] - v[2442] * v[887] - v[2441] * v[888] - v[2440] * v[889]
			- v[2439] * v[890];
		v[2587] = -(v[2457] * v[891]) - v[2456] * v[892] - v[2455] * v[893] - v[2454] * v[894] - v[2453] * v[895]
			- v[2452] * v[896] - v[2451] * v[897] - v[2450] * v[898] - v[2449] * v[899] - v[2448] * v[900] - v[2447] * v[901]
			- v[2446] * v[902] - v[2445] * v[903] - v[2444] * v[904] - v[2442] * v[905] - v[2441] * v[906] - v[2440] * v[907]
			- v[2439] * v[908];
		v[2800] = (-(v[2457] * v[837]) - v[2456] * v[838] - v[2455] * v[839] - v[2454] * v[840] - v[2453] * v[841]
			- v[2452] * v[842] - v[2451] * v[843] - v[2450] * v[844] - v[2449] * v[845] - v[2448] * v[846] - v[2447] * v[847]
			- v[2446] * v[848] - v[2445] * v[849] - v[2444] * v[850] - v[2442] * v[851] - v[2441] * v[852] - v[2440] * v[853]
			- v[2439] * v[854]) / 2e0;
		v[2589] = v[2457] * v[855] + v[2456] * v[856] + v[2455] * v[857] + v[2454] * v[858] + v[2453] * v[859] + v[2452] * v[860]
			+ v[2451] * v[861] + v[2450] * v[862] + v[2449] * v[863] + v[2448] * v[864] + v[2447] * v[865] + v[2446] * v[866]
			+ v[2445] * v[867] + v[2444] * v[868] + v[2442] * v[869] + v[2441] * v[870] + v[2440] * v[871] + v[2439] * v[872];
		v[3884] = v[210] * v[2589];
		v[3883] = v[211] * v[2589];
		v[2590] = -(v[2438] * v[873]) - v[2437] * v[874] - v[2436] * v[875] - v[2435] * v[876] - v[2434] * v[877]
			- v[2433] * v[878] - v[2432] * v[879] - v[2431] * v[880] - v[2430] * v[881] - v[2429] * v[882] - v[2428] * v[883]
			- v[2427] * v[884] - v[2426] * v[885] - v[2425] * v[886] - v[2423] * v[887] - v[2422] * v[888] - v[2421] * v[889]
			- v[2420] * v[890];
		v[2591] = -(v[2438] * v[891]) - v[2437] * v[892] - v[2436] * v[893] - v[2435] * v[894] - v[2434] * v[895]
			- v[2433] * v[896] - v[2432] * v[897] - v[2431] * v[898] - v[2430] * v[899] - v[2429] * v[900] - v[2428] * v[901]
			- v[2427] * v[902] - v[2426] * v[903] - v[2425] * v[904] - v[2423] * v[905] - v[2422] * v[906] - v[2421] * v[907]
			- v[2420] * v[908];
		v[2803] = (-(v[2438] * v[837]) - v[2437] * v[838] - v[2436] * v[839] - v[2435] * v[840] - v[2434] * v[841]
			- v[2433] * v[842] - v[2432] * v[843] - v[2431] * v[844] - v[2430] * v[845] - v[2429] * v[846] - v[2428] * v[847]
			- v[2427] * v[848] - v[2426] * v[849] - v[2425] * v[850] - v[2423] * v[851] - v[2422] * v[852] - v[2421] * v[853]
			- v[2420] * v[854]) / 2e0;
		v[2593] = v[2438] * v[855] + v[2437] * v[856] + v[2436] * v[857] + v[2435] * v[858] + v[2434] * v[859] + v[2433] * v[860]
			+ v[2432] * v[861] + v[2431] * v[862] + v[2430] * v[863] + v[2429] * v[864] + v[2428] * v[865] + v[2427] * v[866]
			+ v[2426] * v[867] + v[2425] * v[868] + v[2423] * v[869] + v[2422] * v[870] + v[2421] * v[871] + v[2420] * v[872];
		v[3886] = v[210] * v[2593];
		v[3885] = v[211] * v[2593];
		v[2594] = -(v[2419] * v[873]) - v[2417] * v[874] - v[2415] * v[875] - v[2413] * v[876] - v[2412] * v[877]
			- v[2411] * v[878] - v[2410] * v[879] - v[2408] * v[880] - v[2406] * v[881] - v[2404] * v[882] - v[2403] * v[883]
			- v[2402] * v[884] - v[2401] * v[885] - v[2400] * v[886] - v[2398] * v[887] - v[2397] * v[888] - v[2396] * v[889]
			- v[2395] * v[890];
		v[2595] = -(v[2419] * v[891]) - v[2417] * v[892] - v[2415] * v[893] - v[2413] * v[894] - v[2412] * v[895]
			- v[2411] * v[896] - v[2410] * v[897] - v[2408] * v[898] - v[2406] * v[899] - v[2404] * v[900] - v[2403] * v[901]
			- v[2402] * v[902] - v[2401] * v[903] - v[2400] * v[904] - v[2398] * v[905] - v[2397] * v[906] - v[2396] * v[907]
			- v[2395] * v[908];
		v[2806] = (-(v[2419] * v[837]) - v[2417] * v[838] - v[2415] * v[839] - v[2413] * v[840] - v[2412] * v[841]
			- v[2411] * v[842] - v[2410] * v[843] - v[2408] * v[844] - v[2406] * v[845] - v[2404] * v[846] - v[2403] * v[847]
			- v[2402] * v[848] - v[2401] * v[849] - v[2400] * v[850] - v[2398] * v[851] - v[2397] * v[852] - v[2396] * v[853]
			- v[2395] * v[854]) / 2e0;
		v[2597] = v[2419] * v[855] + v[2417] * v[856] + v[2415] * v[857] + v[2413] * v[858] + v[2412] * v[859] + v[2411] * v[860]
			+ v[2410] * v[861] + v[2408] * v[862] + v[2406] * v[863] + v[2404] * v[864] + v[2403] * v[865] + v[2402] * v[866]
			+ v[2401] * v[867] + v[2400] * v[868] + v[2398] * v[869] + v[2397] * v[870] + v[2396] * v[871] + v[2395] * v[872];
		v[3888] = v[210] * v[2597];
		v[3887] = v[211] * v[2597];
		v[2601] = -((*ct)*(v[1509] * v[2363] + v[1510] * v[2364] + v[1511] * v[2365])) + (*cn)*(v[2360] * v[2598]
			+ v[2361] * v[2599] + v[2362] * v[2600]) - v[2354] * v[854] - v[2353] * v[872] + v[2356] * v[890] + v[2355] * v[908];
		v[2772] = v[2601] * v[3631];
		v[2605] = -((*ct)*(v[1505] * v[2363] + v[1506] * v[2364] + v[1507] * v[2365])) + (*cn)*(v[2360] * v[2602]
			+ v[2361] * v[2603] + v[2362] * v[2604]) - v[2354] * v[853] - v[2353] * v[871] + v[2356] * v[889] + v[2355] * v[907];
		v[3860] = v[2605] * v[337];
		v[2776] = v[2605] * v[3630];
		v[2609] = -((*ct)*(v[1501] * v[2363] + v[1502] * v[2364] + v[1503] * v[2365])) + (*cn)*(v[2360] * v[2606]
			+ v[2361] * v[2607] + v[2362] * v[2608]) - v[2354] * v[852] - v[2353] * v[870] + v[2356] * v[888] + v[2355] * v[906];
		v[3863] = v[2609] * v[337];
		v[2779] = v[2609] * v[3629];
		v[2613] = -((*ct)*(v[1485] * v[2363] + v[1486] * v[2364] + v[1487] * v[2365])) + (*cn)*(v[2360] * v[2610]
			+ v[2361] * v[2611] + v[2362] * v[2612]) - v[2354] * v[848] - v[2353] * v[866] + v[2356] * v[884] + v[2355] * v[902];
		v[2869] = v[2613] * v[3628];
		v[3897] = v[2869] * v[311];
		v[2617] = -((*ct)*(v[1481] * v[2363] + v[1482] * v[2364] + v[1483] * v[2365])) + (*cn)*(v[2360] * v[2614]
			+ v[2361] * v[2615] + v[2362] * v[2616]) - v[2354] * v[847] - v[2353] * v[865] + v[2356] * v[883] + v[2355] * v[901];
		v[3868] = v[2617] * v[311];
		v[2873] = v[2617] * v[3627];
		v[2621] = -((*ct)*(v[1477] * v[2363] + v[1478] * v[2364] + v[1479] * v[2365])) + (*cn)*(v[2360] * v[2618]
			+ v[2361] * v[2619] + v[2362] * v[2620]) - v[2354] * v[846] - v[2353] * v[864] + v[2356] * v[882] + v[2355] * v[900];
		v[3871] = v[2621] * v[311];
		v[2876] = v[2621] * v[3626];
		v[2625] = -((*ct)*(v[1461] * v[2363] + v[1462] * v[2364] + v[1463] * v[2365])) + (*cn)*(v[2360] * v[2622]
			+ v[2361] * v[2623] + v[2362] * v[2624]) - v[2354] * v[842] - v[2353] * v[860] + v[2356] * v[878] + v[2355] * v[896];
		v[2902] = v[2625] * v[3625];
		v[3903] = v[285] * v[2902];
		v[2629] = -((*ct)*(v[1457] * v[2363] + v[1458] * v[2364] + v[1459] * v[2365])) + (*cn)*(v[2360] * v[2626]
			+ v[2361] * v[2627] + v[2362] * v[2628]) - v[2354] * v[841] - v[2353] * v[859] + v[2356] * v[877] + v[2355] * v[895];
		v[3876] = v[2629] * v[285];
		v[2906] = v[2629] * v[3624];
		v[2633] = -((*ct)*(v[1453] * v[2363] + v[1454] * v[2364] + v[1455] * v[2365])) + (*cn)*(v[2360] * v[2630]
			+ v[2361] * v[2631] + v[2362] * v[2632]) - v[2354] * v[840] - v[2353] * v[858] + v[2356] * v[876] + v[2355] * v[894];
		v[7857] = (*ct)*(-(v[1441] * v[2363]) - v[1442] * v[2364] - v[1443] * v[2365]) + v[210] * v[3963] - v[2354] * v[837]
			- v[2353] * v[855] + v[2356] * v[873] + v[2355] * v[891] + v[3825] * v[984];
		v[7858] = (*ct)*(-(v[1445] * v[2363]) - v[1446] * v[2364] - v[1447] * v[2365]) - v[210] * v[3836] - v[2354] * v[838]
			- v[2353] * v[856] + v[2356] * v[874] + v[2355] * v[892] + v[3833] * v[984];
		v[7859] = (*ct)*(-(v[1449] * v[2363]) - v[1450] * v[2364] - v[1451] * v[2365]) + v[210] * v[3837] - v[2354] * v[839]
			- v[2353] * v[857] + v[2356] * v[875] + v[2355] * v[893];
		v[7860] = v[1544] * v[2106] + v[1543] * v[2110] + v[1966] * (v[2912] + v[1543] * v[3004] + v[1544] * v[3006])
			+ v[2633] * v[3037] + v[2629] * v[3635] + v[2625] * v[3637] + v[1955] * v[3964] + v[1958] * v[3965] + v[1961] * v[3966]
			+ v[153] * v[3967] + v[158] * v[3968] + v[149] * v[3969];
		v[7861] = v[1545] * v[2106] + v[1543] * v[2112] + v[1858] * v[2625] + v[1859] * v[2633] + v[1966] * (v[2905]
			+ v[1543] * v[3003] + v[1545] * v[3007]) + v[2629] * v[3038] + v[1962] * v[3766] + v[1956] * v[3768] + v[1959] * v[3970]
			+ v[150] * v[3971] + v[155] * v[3972] + v[2120] * v[4047];
		v[7862] = v[1545] * v[2110] + v[1544] * v[2112] + v[1860] * v[2629] + v[1861] * v[2633] + v[1966] * (v[2897]
			+ v[1544] * v[3005] + v[1545] * v[3008]) + v[2625] * v[3039] + v[1960] * v[3767] + v[1957] * v[3973] + v[1963] * v[3974]
			+ v[161] * v[3975] + v[2124] * v[4045] + v[2128] * v[4046];
		v[7863] = (*ct)*(-(v[1465] * v[2363]) - v[1466] * v[2364] - v[1467] * v[2365]) + v[211] * v[3963] - v[2354] * v[843]
			- v[2353] * v[861] + v[2356] * v[879] + v[2355] * v[897] + v[3825] * v[985];
		v[7864] = (*ct)*(-(v[1469] * v[2363]) - v[1470] * v[2364] - v[1471] * v[2365]) - v[211] * v[3836] - v[2354] * v[844]
			- v[2353] * v[862] + v[2356] * v[880] + v[2355] * v[898] + v[3833] * v[985];
		v[7865] = (*ct)*(-(v[1473] * v[2363]) - v[1474] * v[2364] - v[1475] * v[2365]) + v[211] * v[3837] - v[2354] * v[845]
			- v[2353] * v[863] + v[2356] * v[881] + v[2355] * v[899];
		v[7866] = v[1541] * v[2134] + v[1540] * v[2138] + v[1993] * (v[2879] + v[1540] * v[2977] + v[1541] * v[2979])
			+ v[2621] * v[3040] + v[2617] * v[3640] + v[2613] * v[3642] + v[1982] * v[3976] + v[1985] * v[3977] + v[1988] * v[3978]
			+ v[172] * v[3979] + v[177] * v[3980] + v[168] * v[3981];
		v[7867] = v[1542] * v[2134] + v[1540] * v[2140] + v[1864] * v[2613] + v[1865] * v[2621] + v[1993] * (v[2872]
			+ v[1540] * v[2976] + v[1542] * v[2980]) + v[2617] * v[3041] + v[1989] * v[3763] + v[1983] * v[3765] + v[1986] * v[3982]
			+ v[169] * v[3983] + v[174] * v[3984] + v[2148] * v[4044];
		v[7868] = v[1542] * v[2138] + v[1541] * v[2140] + v[1866] * v[2617] + v[1867] * v[2621] + v[1993] * (v[2864]
			+ v[1541] * v[2978] + v[1542] * v[2981]) + v[2613] * v[3042] + v[1987] * v[3764] + v[1984] * v[3985] + v[1990] * v[3986]
			+ v[180] * v[3987] + v[2152] * v[4042] + v[2156] * v[4043];
		v[7869] = (*ct)*(-(v[1489] * v[2363]) - v[1490] * v[2364] - v[1491] * v[2365]) + v[2503] - v[3988] - v[2354] * v[849]
			- v[2353] * v[867] + v[2356] * v[885] + v[2355] * v[903];
		v[7870] = (*ct)*(-(v[1493] * v[2363]) - v[1494] * v[2364] - v[1495] * v[2365]) + v[3838] - v[2354] * v[850]
			- v[2353] * v[868] + v[2356] * v[886] + v[2355] * v[904];
		v[7871] = (*ct)*(-(v[1497] * v[2363]) - v[1498] * v[2364] - v[1499] * v[2365]) + v[2515] - v[3892] - v[2354] * v[851]
			- v[2353] * v[869] + v[2356] * v[887] + v[2355] * v[905];
		v[7872] = v[1538] * v[2162] + v[1537] * v[2166] + v[2074] * (v[2782] + v[1537] * v[2948] + v[1538] * v[2950])
			+ v[2609] * v[3043] + v[2605] * v[3645] + v[2601] * v[3647] + v[2063] * v[3989] + v[2066] * v[3990] + v[2069] * v[3991]
			+ v[232] * v[3992] + v[237] * v[3993] + v[228] * v[3994];
		v[7873] = v[1539] * v[2162] + v[1537] * v[2168] + v[1870] * v[2601] + v[1871] * v[2609] + v[2074] * (v[2775]
			+ v[1537] * v[2947] + v[1539] * v[2951]) + v[2605] * v[3044] + v[2070] * v[3744] + v[2064] * v[3746] + v[2067] * v[3995]
			+ v[229] * v[3996] + v[234] * v[3997] + v[2176] * v[4041];
		v[7874] = v[1539] * v[2166] + v[1538] * v[2168] + v[1872] * v[2605] + v[1873] * v[2609] + v[2074] * (v[2767]
			+ v[1538] * v[2949] + v[1539] * v[2952]) + v[2601] * v[3045] + v[2068] * v[3745] + v[2065] * v[3998] + v[2071] * v[3999]
			+ v[240] * v[4000] + v[2180] * v[4039] + v[2184] * v[4040];
		v[3879] = v[2633] * v[285];
		v[2909] = v[2633] * v[3623];
		v[2204] = v[2204] + 2e0*v[2461] * v[377] + v[2463] * v[3852];
		v[2200] = v[2200] + 2e0*v[2462] * v[376] + v[375] * v[3853];
		v[2196] = v[2196] + 2e0*v[2480] * v[375] + v[376] * v[3853];
		v[3855] = (v[2193] * v[2635] * v[372] + (v[1553] * v[2187] + v[1554] * v[2189] + v[1555] * v[2191] + v[2196] * v[360]
			+ v[2200] * v[361] + v[2204] * v[362])*v[373] - v[1643] * v[2634] * v[3854]) / v[1644];
		v[2638] = (*epsn)*v[2362] + v[2204] * v[374] + v[2191] * v[3743] + v[1555] * v[3807] + v[362] * v[3855];
		v[2640] = (*epsn)*v[2361] + v[2200] * v[374] + v[2189] * v[3743] + v[1554] * v[3807] + v[361] * v[3855];
		v[2642] = (*epsn)*v[2360] + v[2196] * v[374] + v[2187] * v[3743] + v[1553] * v[3807] + v[360] * v[3855];
		v[2568] = v[2568] - v[2638];
		v[2558] = v[2558] - v[2640];
		v[2548] = v[2548] - v[2642];
		v[2643] = v[2601] * v[346] + v[1537] * v[3856];
		v[2644] = v[2601] * v[340] + v[1537] * v[3857];
		v[2645] = v[1659] * v[2601] + (v[234] * v[2643]) / v[337] + v[1537] * (v[2091] / v[337] + v[2074] * v[3861]);
		v[2646] = v[1666] * v[2601] + (v[228] * v[2644]) / v[337] + v[1537] * (v[2083] / v[337] + v[2074] * v[4201]);
		v[2647] = (v[239] * v[2643] + v[237] * v[2644] + v[1655] * v[2601] * v[337] + v[1537] * (v[2099] + v[3806] * v[4202]))
			/ v[337];
		v[2648] = v[2605] * v[348] + v[1538] * v[3858];
		v[2649] = v[2605] * v[340] + v[1538] * v[3857];
		v[2652] = v[1654] * v[2605] + (v[240] * v[2648]) / v[337] + v[1538] * (v[2097] / v[337] + v[2074] * v[3859]);
		v[2653] = v[2643] + v[2648];
		v[2654] = v[2645] + v[2652];
		v[2656] = (v[1657] * v[2065] + v[228] * v[2649] + v[230] * v[2653] + v[2943] * v[3806] + v[1665] * v[3860] + v[1538] *
			(v[2079] + v[3806] * v[3862])) / v[337];
		v[2659] = (v[235] * v[2648] + v[232] * v[2649] + v[1658] * v[3860] + v[1538] * (v[2088] + v[3806] * v[4203])) / v[337];
		v[2660] = v[2609] * v[346] + v[1539] * v[3856];
		v[2661] = v[2609] * v[348] + v[1539] * v[3858];
		v[2662] = v[2649] + v[2660];
		v[2665] = (v[229] * v[2660] + v[230] * v[2661] + v[1663] * v[3863] + v[1539] * (v[2080] + v[3806] * v[4204])) / v[337];
		v[2666] = v[2644] + v[2661];
		v[2668] = (v[1668] * v[2068] + v[234] * v[2660] + v[235] * v[2666] + v[2942] * v[3806] + v[1670] * v[3863] + v[1539] *
			(v[2087] + v[3806] * v[4205])) / v[337];
		v[2669] = v[2656] + v[2668];
		v[2671] = (v[1669] * v[2070] + v[240] * v[2661] + v[239] * v[2662] + v[2941] * v[3806] + v[1674] * v[3863] + v[1539] *
			(v[2095] + v[3806] * v[4206])) / v[337];
		v[2672] = v[2646] + v[2671];
		v[2673] = v[2613] * v[320] + v[1540] * v[3864];
		v[2674] = v[2613] * v[314] + v[1540] * v[3865];
		v[2675] = v[1686] * v[2613] + (v[174] * v[2673]) / v[311] + v[1540] * (v[2022] / v[311] + v[1993] * v[3869]);
		v[2676] = v[1693] * v[2613] + (v[168] * v[2674]) / v[311] + v[1540] * (v[2017] / v[311] + v[1993] * v[4207]);
		v[2677] = (v[179] * v[2673] + v[177] * v[2674] + v[1682] * v[2613] * v[311] + v[1540] * (v[2027] + v[3805] * v[4208]))
			/ v[311];
		v[2678] = v[2617] * v[322] + v[1541] * v[3866];
		v[2679] = v[2617] * v[314] + v[1541] * v[3865];
		v[2682] = v[1681] * v[2617] + (v[180] * v[2678]) / v[311] + v[1541] * (v[2025] / v[311] + v[1993] * v[3867]);
		v[2683] = v[2673] + v[2678];
		v[2684] = v[2675] + v[2682];
		v[2686] = (v[1684] * v[1984] + v[168] * v[2679] + v[170] * v[2683] + v[2972] * v[3805] + v[1692] * v[3868] + v[1541] *
			(v[2013] + v[3805] * v[3870])) / v[311];
		v[2689] = (v[175] * v[2678] + v[172] * v[2679] + v[1685] * v[3868] + v[1541] * (v[2019] + v[3805] * v[4209])) / v[311];
		v[2690] = v[2621] * v[320] + v[1542] * v[3864];
		v[2691] = v[2621] * v[322] + v[1542] * v[3866];
		v[2692] = v[2679] + v[2690];
		v[2695] = (v[169] * v[2690] + v[170] * v[2691] + v[1690] * v[3871] + v[1542] * (v[2014] + v[3805] * v[4210])) / v[311];
		v[2696] = v[2674] + v[2691];
		v[2698] = (v[1695] * v[1987] + v[174] * v[2690] + v[175] * v[2696] + v[2971] * v[3805] + v[1697] * v[3871] + v[1542] *
			(v[2018] + v[3805] * v[4211])) / v[311];
		v[2699] = v[2686] + v[2698];
		v[2701] = (v[1696] * v[1989] + v[180] * v[2691] + v[179] * v[2692] + v[2970] * v[3805] + v[1701] * v[3871] + v[1542] *
			(v[2023] + v[3805] * v[4212])) / v[311];
		v[2702] = v[2676] + v[2701];
		v[2703] = v[2625] * v[294] + v[1543] * v[3872];
		v[2704] = v[2625] * v[288] + v[1543] * v[3873];
		v[2705] = v[1713] * v[2625] + (v[155] * v[2703]) / v[285] + v[1543] * (v[2007] / v[285] + v[1966] * v[3877]);
		v[2706] = v[1720] * v[2625] + (v[149] * v[2704]) / v[285] + v[1543] * (v[2002] / v[285] + v[1966] * v[4213]);
		v[2707] = (v[160] * v[2703] + v[158] * v[2704] + v[1709] * v[2625] * v[285] + v[1543] * (v[2012] + v[3804] * v[4214]))
			/ v[285];
		v[2708] = v[2629] * v[296] + v[1544] * v[3874];
		v[2709] = v[2629] * v[288] + v[1544] * v[3873];
		v[2712] = v[1708] * v[2629] + (v[161] * v[2708]) / v[285] + v[1544] * (v[2010] / v[285] + v[1966] * v[3875]);
		v[2713] = v[2703] + v[2708];
		v[2714] = v[2705] + v[2712];
		v[2716] = (v[1711] * v[1957] + v[149] * v[2709] + v[151] * v[2713] + v[2999] * v[3804] + v[1719] * v[3876] + v[1544] *
			(v[1998] + v[3804] * v[3878])) / v[285];
		v[2719] = (v[156] * v[2708] + v[153] * v[2709] + v[1712] * v[3876] + v[1544] * (v[2004] + v[3804] * v[4215])) / v[285];
		v[2720] = v[2633] * v[294] + v[1545] * v[3872];
		v[2721] = v[2633] * v[296] + v[1545] * v[3874];
		v[2722] = v[2709] + v[2720];
		v[2725] = (v[150] * v[2720] + v[151] * v[2721] + v[1717] * v[3879] + v[1545] * (v[1999] + v[3804] * v[4216])) / v[285];
		v[2726] = v[2704] + v[2721];
		v[2728] = (v[1722] * v[1960] + v[155] * v[2720] + v[156] * v[2726] + v[2998] * v[3804] + v[1724] * v[3879] + v[1545] *
			(v[2003] + v[3804] * v[4217])) / v[285];
		v[2729] = v[2716] + v[2728];
		v[2731] = (v[1723] * v[1962] + v[161] * v[2721] + v[160] * v[2722] + v[2997] * v[3804] + v[1728] * v[3879] + v[1545] *
			(v[2008] + v[3804] * v[4218])) / v[285];
		v[2732] = v[2706] + v[2731];
		v[2733] = -(v[2502] * v[950]);
		v[2734] = -(v[2502] * v[949]);
		v[2735] = -(v[2502] * v[951]);
		v[2736] = -(v[2505] * v[951]);
		v[2737] = -(v[2505] * v[950]);
		v[2738] = -(v[2505] * v[949]);
		v[2739] = -(v[2507] * v[951]);
		v[2740] = -(v[2507] * v[949]);
		v[2741] = -(v[2507] * v[950]);
		v[2742] = -(v[2508] * v[951]);
		v[2743] = -(v[2508] * v[950]);
		v[2744] = -(v[2508] * v[949]);
		v[2745] = -(v[2511] * v[951]);
		v[2746] = -(v[2511] * v[949]);
		v[2747] = -(v[2511] * v[950]);
		v[2748] = -(v[2512] * v[949]);
		v[2749] = -(v[2512] * v[951]);
		v[2750] = -(v[2512] * v[950]);
		v[2751] = -(v[2514] * v[951]);
		v[2752] = -(v[2514] * v[950]);
		v[2753] = -(v[2514] * v[949]);
		v[2754] = -(v[2516] * v[950]);
		v[2755] = -(v[2516] * v[951]);
		v[2756] = -(v[2516] * v[949]);
		v[2757] = -(v[1648] * v[2086]) - v[1647] * v[2094] - v[1645] * v[2102] - v[252] * v[2638] - v[249] * v[2640]
			- v[246] * v[2642] + v[2502] * v[710] + v[2505] * v[711] + v[2507] * v[712] + v[2508] * v[725] + v[2511] * v[726]
			+ v[2512] * v[727] + v[2514] * v[740] + v[2516] * v[741] + v[2517] * v[742];
		v[2758] = -(v[1648] * v[2084]) - v[1647] * v[2092] - v[1645] * v[2100] - v[250] * v[2638] - v[247] * v[2640]
			- v[244] * v[2642] + v[2502] * v[704] + v[2505] * v[705] + v[2507] * v[706] + v[2508] * v[719] + v[2511] * v[720]
			+ v[2512] * v[721] + v[2514] * v[734] + v[2516] * v[735] + v[2517] * v[736];
		v[2759] = -(v[2517] * v[950]);
		v[2760] = -(v[2517] * v[951]);
		v[2761] = -(v[2517] * v[949]);
		v[2569] = v[2569] - v[263] * v[2638] + v[2594] * v[412] + v[2595] * v[420];
		v[2570] = v[2570] + v[1440] * v[2595] - v[261] * v[2638];
		v[2571] = v[2571] - v[260] * v[2638] + v[2594] * v[410] + v[2595] * v[418];
		v[2559] = v[2559] - v[263] * v[2640] + v[2590] * v[412] + v[2591] * v[420];
		v[2560] = v[2560] + v[1440] * v[2591] - v[261] * v[2640];
		v[2561] = v[2561] - v[260] * v[2640] + v[2590] * v[410] + v[2591] * v[418];
		v[2762] = v[1546] * v[2086] + v[1530] * v[2094] + v[1525] * v[2102] + v[246] * v[2586] + v[249] * v[2590]
			+ v[252] * v[2594];
		v[2763] = v[1546] * v[2084] + v[1530] * v[2092] + v[1525] * v[2100] + v[244] * v[2586] + v[247] * v[2590]
			+ v[250] * v[2594];
		v[2549] = v[2549] - v[263] * v[2642] + v[2586] * v[412] + v[2587] * v[420];
		v[2550] = v[2550] + v[1440] * v[2587] - v[261] * v[2642];
		v[2551] = v[2551] - v[260] * v[2642] + v[2586] * v[410] + v[2587] * v[418];
		v[2764] = v[1547] * v[2086] + v[1532] * v[2094] + v[1526] * v[2102] + v[246] * v[2587] + v[249] * v[2591]
			+ v[252] * v[2595];
		v[2765] = v[1547] * v[2084] + v[1532] * v[2092] + v[1526] * v[2100] + v[244] * v[2587] + v[247] * v[2591]
			+ v[250] * v[2595];
		v[2766] = QABi[2][2] * v[2569] + QABi[2][1] * v[2570] + QABi[2][0] * v[2571] + v[1870] * v[2648] + v[2772] * v[350]
			+ v[2661] * v[3647] + v[348] * (v[4000] + v[2074] * v[4219]);
		v[2769] = QABi[1][2] * v[2569] + QABi[1][1] * v[2570] + QABi[1][0] * v[2571] + v[2643] * v[3045] + (v[1669] * v[2159])
			/ v[337] + v[2776] * v[347] + v[2176] * v[3630] + v[2662] * v[3647] + v[2074] * v[4220];
		v[2771] = QABi[0][2] * v[2569] + QABi[0][1] * v[2570] + QABi[0][0] * v[2571] + v[2644] * v[3045] + v[2779] * v[335]
			+ v[340] * (v[3993] + v[2074] * v[4221]);
		v[2773] = QABi[2][2] * v[2559] + QABi[2][1] * v[2560] + QABi[2][0] * v[2561] + v[2648] * v[3044] + (v[1668] * v[2163])
			/ v[337] + v[2184] * v[3631] + v[2666] * v[3645] + v[2772] * v[3648] + v[2074] * v[4222];
		v[2777] = QABi[1][2] * v[2559] + QABi[1][1] * v[2560] + QABi[1][0] * v[2561] + v[1872] * v[2643] + v[2776] * v[341]
			+ v[2660] * v[3645] + v[346] * (v[3997] + v[2074] * v[4223]);
		v[2780] = QABi[0][2] * v[2559] + QABi[0][1] * v[2560] + QABi[0][0] * v[2561] + v[2649] * v[3044] + v[2779] * v[336]
			+ v[340] * (v[3992] + v[2074] * v[4224]);
		v[2781] = QABi[2][2] * v[2549] + QABi[2][1] * v[2550] + QABi[2][0] * v[2551] + v[1871] * v[2653] + v[2661] * v[3043] +
			(v[1657] * v[2164]) / v[337] + v[2180] * v[3631] + v[2772] * v[3646] + v[2074] * v[4225];
		v[2783] = QABi[1][2] * v[2549] + QABi[1][1] * v[2550] + QABi[1][0] * v[2551] + v[2660] * v[3043] + v[2776] * v[3644]
			+ v[346] * (v[2074] * v[3882] + v[3996]);
		v[2786] = QABi[0][2] * v[2549] + QABi[0][1] * v[2550] + QABi[0][0] * v[2551] + v[1873] * v[2644] + v[1871] * v[2649]
			+ v[2779] * v[338] + v[340] * (v[3994] + v[2074] * v[4226]);
		v[2787] = -(v[1661] * v[2053]) + 2e0*v[1660] * v[2055] - v[1672] * v[2056] + v[2733] + v[2742] + v[2748] + v[2759]
			+ 2e0*v[2659] * v[666] - v[2669] * v[669] - v[2654] * v[672];
		v[2788] = v[1676] * v[2053] - v[1672] * v[2055] + 2e0*v[1667] * v[2056] - v[2737] - v[2740] - v[2745] - v[2760]
			- v[2669] * v[666] + 2e0*v[2665] * v[669] + v[2672] * v[672];
		v[2789] = -2e0*v[1750] * v[2054] + v[225] * v[2783] - v[2733] * v[659] + v[2737] * v[660] - v[2741] * v[661];
		v[2790] = 2e0*v[1653] * v[2053] - v[1661] * v[2055] + v[1676] * v[2056] - v[2734] - v[2746] - v[2751] - v[2754]
			- v[2654] * v[666] + v[2672] * v[669] + 2e0*v[2647] * v[672];
		v[2791] = -2e0*v[1748] * v[2054] + v[225] * v[2781] - v[2734] * v[659] + v[2738] * v[660] - v[2740] * v[661];
		v[2792] = -2e0*v[1743] * v[2054] + v[225] * v[2780] - v[2742] * v[659] + v[2745] * v[660] - v[2749] * v[661];
		v[2793] = v[2739] + v[2750];
		v[2794] = -2e0*v[1738] * v[2054] + v[225] * v[2773] - v[2744] * v[659] + v[2746] * v[660] - v[2748] * v[661];
		v[2795] = -2e0*v[1737] * v[2054] + v[225] * v[2771] - v[2751] * v[659] + v[2755] * v[660] - v[2760] * v[661];
		v[2796] = -2e0*v[1735] * v[2054] + v[225] * v[2769] - v[2752] * v[659] + v[2754] * v[660] - v[2759] * v[661];
		v[2797] = v[2743] + v[2753];
		v[2798] = v[2736] + v[2756];
		v[2799] = v[210] * v[2642] + v[2800];
		v[2801] = v[211] * v[2642] - v[2800];
		v[2553] = v[2553] + v[2799];
		v[2552] = v[2552] + v[2801];
		v[2802] = v[210] * v[2640] + v[2803];
		v[2804] = v[211] * v[2640] - v[2803];
		v[2563] = v[2563] + v[2802];
		v[2562] = v[2562] + v[2804];
		v[2805] = v[210] * v[2638] + v[2806];
		v[2807] = v[211] * v[2638] - v[2806];
		v[2573] = v[2573] + v[2805];
		v[2572] = v[2572] + v[2807];
		v[2808] = v[2519] * v[954];
		v[2809] = v[2519] * v[953];
		v[2810] = v[2519] * v[952];
		v[2811] = v[2520] * v[953];
		v[2812] = v[2520] * v[954];
		v[2813] = v[2520] * v[952];
		v[2814] = v[2521] * v[953];
		v[2815] = v[2521] * v[954];
		v[2816] = v[2521] * v[952];
		v[2817] = v[2522] * v[957];
		v[2818] = v[2522] * v[956];
		v[2819] = v[2522] * v[955];
		v[2820] = v[2523] * v[956];
		v[2821] = v[2523] * v[957];
		v[2822] = v[2523] * v[955];
		v[2823] = v[2524] * v[956];
		v[2824] = v[2524] * v[957];
		v[2825] = v[2524] * v[955];
		v[2826] = v[2525] * v[954];
		v[2827] = v[2525] * v[953];
		v[2828] = v[2525] * v[952];
		v[2829] = v[2526] * v[954];
		v[2830] = v[2526] * v[952];
		v[2831] = v[2526] * v[953];
		v[2832] = v[2527] * v[952];
		v[2833] = v[2527] * v[954];
		v[2834] = v[2527] * v[953];
		v[2835] = v[2528] * v[957];
		v[2836] = v[2528] * v[956];
		v[2837] = v[2528] * v[955];
		v[2838] = v[2529] * v[957];
		v[2839] = v[2529] * v[955];
		v[2840] = v[2529] * v[956];
		v[2841] = v[2530] * v[955];
		v[2842] = v[2530] * v[957];
		v[2843] = v[2530] * v[956];
		v[2844] = v[2532] * v[953];
		v[2845] = v[2532] * v[952];
		v[2846] = v[2532] * v[954];
		v[2847] = v[2533] * v[954];
		v[2848] = v[2533] * v[953];
		v[2849] = v[2533] * v[952];
		v[2850] = v[2534] * v[954];
		v[2851] = v[2534] * v[952];
		v[2852] = v[2534] * v[953];
		v[2853] = v[2535] * v[956];
		v[2854] = v[2535] * v[955];
		v[2855] = v[2535] * v[957];
		v[2856] = v[2536] * v[957];
		v[2857] = v[2536] * v[956];
		v[2858] = v[2536] * v[955];
		v[2859] = v[2537] * v[957];
		v[2860] = v[2537] * v[955];
		v[2861] = v[2537] * v[956];
		v[2554] = v[2554] + v[215] * v[2801] + v[1439] * v[3883];
		v[2555] = v[2555] + v[214] * v[2801] - v[3883] * v[993];
		v[2556] = v[2556] + v[215] * v[2799] + v[1439] * v[3884];
		v[2557] = v[2557] + v[214] * v[2799] - v[3884] * v[993];
		v[2564] = v[2564] + v[215] * v[2804] + v[1439] * v[3885];
		v[2565] = v[2565] + v[214] * v[2804] - v[3885] * v[993];
		v[2566] = v[2566] + v[215] * v[2802] + v[1439] * v[3886];
		v[2567] = v[2567] + v[214] * v[2802] - v[3886] * v[993];
		v[2574] = v[2574] + v[215] * v[2807] + v[1439] * v[3887];
		v[2575] = v[2575] + v[214] * v[2807] - v[3887] * v[993];
		v[2576] = v[2576] + v[215] * v[2805] + v[1439] * v[3888];
		v[2577] = v[2577] + v[214] * v[2805] - v[3888] * v[993];
		v[2862] = (v[1645] * v[2040] - v[1645] * v[2042] + v[1647] * v[2044] - v[1647] * v[2046] + v[1648] * v[2048]
			- v[1648] * v[2050] + v[2503] * v[273] + v[2509] * v[274] - v[2503] * v[276] - v[2509] * v[277] - v[2468] * v[3657]
			+ v[2469] * v[3657] - v[2642] * v[387] + v[2642] * v[388] - v[2029] * v[3889] + v[2031] * v[3889] - v[2033] * v[3890]
			+ v[2035] * v[3890] - v[2037] * v[3891] + v[2039] * v[3891] - v[275] * v[3892] + v[278] * v[3892] - v[2640] * v[390]
			+ v[2640] * v[391] - v[2638] * v[393] + v[2638] * v[394] + v[2597] * v[399] - v[2597] * v[400] + v[2593] * v[401]
			- v[2593] * v[402] + v[2589] * v[403] - v[2589] * v[404] + v[3826] * v[4227] + v[2519] * v[605] + v[2520] * v[606]
			+ v[2521] * v[607] - v[2522] * v[608] - v[2523] * v[609] - v[2524] * v[610] + v[2525] * v[617] + v[2526] * v[618]
			+ v[2527] * v[619] - v[2528] * v[620] - v[2529] * v[621] - v[2530] * v[622] + v[2532] * v[629] + v[2533] * v[630]
			+ v[2534] * v[631] - v[2535] * v[632] - v[2536] * v[633] - v[2537] * v[634] + (v[1528] * (v[2028] - v[2030]) + v[1536] *
			(v[2032] - v[2034]) + v[1549] * (v[2036] - v[2038]))*v[993]) / 2e0;
		v[2863] = QBAi[2][1] * v[2574] + QBAi[2][0] * v[2575] + v[1864] * v[2678] + v[2869] * v[324] + v[2691] * v[3642]
			+ v[322] * (v[3987] + v[1993] * v[4228]);
		v[2866] = (v[1696] * v[2131] + v[2692] * v[309] + QBAi[1][1] * v[2574] * v[311] + QBAi[1][0] * v[2575] * v[311]
			+ v[2148] * v[320] + v[2873] * v[311] * v[321] + v[2673] * v[324] + v[3805] * v[4229]) / v[311];
		v[2868] = QBAi[0][1] * v[2574] + QBAi[0][0] * v[2575] + v[2674] * v[3042] + v[2876] * v[309] + v[314] * (v[3980]
			+ v[1993] * v[4230]);
		v[2870] = (v[1695] * v[2135] + v[2696] * v[310] + v[2678] * v[315] + v[2156] * v[322] + v[2564] * v[3895]
			+ v[2565] * v[3896] + v[3643] * v[3897] + v[3805] * v[4231]) / v[311];
		v[2874] = QBAi[1][1] * v[2564] + QBAi[1][0] * v[2565] + v[1866] * v[2673] + v[2873] * v[315] + v[2690] * v[3640]
			+ v[320] * (v[3984] + v[1993] * v[4232]);
		v[2877] = QBAi[0][1] * v[2564] + QBAi[0][0] * v[2565] + v[2679] * v[3041] + v[2876] * v[310] + v[314] * (v[3979]
			+ v[1993] * v[4233]);
		v[2878] = (v[1684] * v[2136] + v[2691] * v[312] + v[2152] * v[322] + v[2683] * v[3639] + v[2554] * v[3895]
			+ v[2555] * v[3896] + v[3641] * v[3897] + v[3805] * v[4234]) / v[311];
		v[2880] = QBAi[1][1] * v[2554] + QBAi[1][0] * v[2555] + v[2690] * v[3040] + v[2873] * v[3639] + v[320] *
			(v[1993] * v[3898] + v[3983]);
		v[2883] = QBAi[0][1] * v[2554] + QBAi[0][0] * v[2555] + v[1867] * v[2674] + v[1865] * v[2679] + v[2876] * v[312]
			+ v[314] * (v[3981] + v[1993] * v[4235]);
		v[2884] = -(v[1688] * v[1972]) + 2e0*v[1687] * v[1974] - v[1699] * v[1975] + v[2814] + v[2826] + v[2832] + v[2844]
			+ 2e0*v[2689] * v[495] - v[2699] * v[498] - v[2684] * v[501];
		v[2885] = v[1703] * v[1972] - v[1699] * v[1974] + 2e0*v[1694] * v[1975] - v[2815] - v[2829] - v[2848] - v[2851]
			- v[2699] * v[495] + 2e0*v[2695] * v[498] + v[2702] * v[501];
		v[2886] = -2e0*v[1791] * v[1973] + v[165] * v[2880] - v[2844] * v[488] + v[2848] * v[489] - v[2852] * v[490];
		v[2887] = 2e0*v[1680] * v[1972] - v[1688] * v[1974] + v[1703] * v[1975] - v[2808] - v[2811] - v[2830] - v[2845]
			- v[2684] * v[495] + v[2702] * v[498] + 2e0*v[2677] * v[501];
		v[2888] = -2e0*v[1789] * v[1973] + v[165] * v[2878] - v[2845] * v[488] + v[2849] * v[489] - v[2851] * v[490];
		v[2889] = -2e0*v[1788] * v[1973] + v[165] * v[2877] - v[2826] * v[488] + v[2829] * v[489] - v[2833] * v[490];
		v[2890] = v[2834] + v[2850];
		v[2891] = -2e0*v[1783] * v[1973] + v[165] * v[2870] - v[2828] * v[488] + v[2830] * v[489] - v[2832] * v[490];
		v[2892] = -2e0*v[1782] * v[1973] + v[165] * v[2868] - v[2808] * v[488] + v[2812] * v[489] - v[2815] * v[490];
		v[2893] = -2e0*v[1780] * v[1973] + v[165] * v[2866] - v[2809] * v[488] + v[2811] * v[489] - v[2814] * v[490];
		v[2894] = v[2810] + v[2827];
		v[2895] = v[2813] + v[2847];
		v[2896] = QAAi[2][1] * v[2576] + QAAi[2][0] * v[2577] + v[1858] * v[2708] + v[2902] * v[298] + v[2721] * v[3637]
			+ v[296] * (v[3975] + v[1966] * v[4236]);
		v[2899] = (v[1723] * v[2103] + v[2722] * v[283] + QAAi[1][1] * v[2576] * v[285] + QAAi[1][0] * v[2577] * v[285]
			+ v[2120] * v[294] + v[285] * v[2906] * v[295] + v[2703] * v[298] + v[3804] * v[4237]) / v[285];
		v[2901] = QAAi[0][1] * v[2576] + QAAi[0][0] * v[2577] + v[283] * v[2909] + v[2704] * v[3039] + v[288] * (v[3968]
			+ v[1966] * v[4238]);
		v[2903] = (v[1722] * v[2107] + v[2726] * v[284] + v[2708] * v[289] + v[2128] * v[296] + v[2566] * v[3901]
			+ v[2567] * v[3902] + v[3638] * v[3903] + v[3804] * v[4239]) / v[285];
		v[2907] = QAAi[1][1] * v[2566] + QAAi[1][0] * v[2567] + v[1860] * v[2703] + v[289] * v[2906] + v[2720] * v[3635]
			+ v[294] * (v[3972] + v[1966] * v[4240]);
		v[2910] = QAAi[0][1] * v[2566] + QAAi[0][0] * v[2567] + v[284] * v[2909] + v[2709] * v[3038] + v[288] * (v[3967]
			+ v[1966] * v[4241]);
		v[2911] = (v[1711] * v[2108] + v[2721] * v[286] + v[2124] * v[296] + v[2713] * v[3634] + v[2556] * v[3901]
			+ v[2557] * v[3902] + v[3636] * v[3903] + v[3804] * v[4242]) / v[285];
		v[2913] = QAAi[1][1] * v[2556] + QAAi[1][0] * v[2557] + v[2720] * v[3037] + v[2906] * v[3634] + v[294] *
			(v[1966] * v[3904] + v[3971]);
		v[2916] = QAAi[0][1] * v[2556] + QAAi[0][0] * v[2557] + v[1861] * v[2704] + v[1859] * v[2709] + v[286] * v[2909]
			+ v[288] * (v[3969] + v[1966] * v[4243]);
		v[2917] = -(v[1715] * v[1945]) + 2e0*v[1714] * v[1947] - v[1726] * v[1948] + v[2823] + v[2835] + v[2841] + v[2853]
			+ 2e0*v[2719] * v[450] - v[2729] * v[453] - v[2714] * v[456];
		v[2918] = v[1730] * v[1945] - v[1726] * v[1947] + 2e0*v[1721] * v[1948] - v[2824] - v[2838] - v[2857] - v[2860]
			- v[2729] * v[450] + 2e0*v[2725] * v[453] + v[2732] * v[456];
		v[2919] = -2e0*v[1809] * v[1946] + v[146] * v[2913] - v[2853] * v[443] + v[2857] * v[444] - v[2861] * v[445];
		v[2920] = 2e0*v[1707] * v[1945] - v[1715] * v[1947] + v[1730] * v[1948] - v[2817] - v[2820] - v[2839] - v[2854]
			- v[2714] * v[450] + v[2732] * v[453] + 2e0*v[2707] * v[456];
		v[2921] = -2e0*v[1807] * v[1946] + v[146] * v[2911] - v[2854] * v[443] + v[2858] * v[444] - v[2860] * v[445];
		v[2922] = -2e0*v[1806] * v[1946] + v[146] * v[2910] - v[2835] * v[443] + v[2838] * v[444] - v[2842] * v[445];
		v[2923] = v[2843] + v[2859];
		v[2924] = -2e0*v[1801] * v[1946] + v[146] * v[2903] - v[2837] * v[443] + v[2839] * v[444] - v[2841] * v[445];
		v[2925] = -2e0*v[1800] * v[1946] + v[146] * v[2901] - v[2817] * v[443] + v[2821] * v[444] - v[2824] * v[445];
		v[2926] = -2e0*v[1798] * v[1946] + v[146] * v[2899] - v[2818] * v[443] + v[2820] * v[444] - v[2823] * v[445];
		v[2927] = v[2819] + v[2836];
		v[2928] = v[2822] + v[2856];
		v[2929] = v[1440] * (-(v[1648] * v[2085]) - v[1647] * v[2093] - v[1645] * v[2101] - v[251] * v[2638] - v[248] * v[2640]
			- v[245] * v[2642] + v[2764] * v[3633] - v[2765] * v[3704] + v[2502] * v[707] + v[2505] * v[708] + v[2507] * v[709]
			+ v[2508] * v[722] + v[2511] * v[723] + v[2512] * v[724] + v[2514] * v[737] + v[2516] * v[738] + v[2517] * v[739])
			- v[994] * (v[1547] * v[2085] + v[1532] * v[2093] + v[1526] * v[2101] + v[245] * v[2587] + v[248] * v[2591]
				+ v[251] * v[2595] + (v[1438] * v[2758] - v[2763] * v[992])*(*xfac) - (v[1438] * v[2762] + v[2757] * v[992])*(*zfac)
				);
		v[2930] = (-2e0*v[2645] + 2e0*v[2652] - v[2735] * v[662] - 2e0*v[2733] * v[667] - 2e0*v[2734] * v[673]
			- 2e0*v[2742] * v[677] - v[2743] * v[681] - 2e0*v[2744] * v[686] - 2e0*v[2751] * v[690] - 2e0*v[2752] * v[694]
			- v[2753] * v[699]) / 2e0;
		v[2931] = (-2e0*v[1753] * v[2054] + v[225] * v[2786] - v[2735] * v[659] + v[2736] * v[660] - v[2739] * v[661]) / 2e0;
		v[2933] = -v[2646] + v[2671] + (v[2736] * v[662]) / 2e0 + v[2737] * v[667] + v[2738] * v[673] + v[2745] * v[677] +
			(v[2747] * v[681]) / 2e0 + v[2746] * v[686] + v[2755] * v[690] + v[2754] * v[694] + (v[2756] * v[699]) / 2e0;
		v[2934] = (-2e0*v[1741] * v[2054] + v[225] * v[2777] - v[2743] * v[659] + v[2747] * v[660] - v[2750] * v[661]) / 2e0;
		v[2935] = (-2e0*v[2656] + 2e0*v[2668] - v[2739] * v[662] - 2e0*v[2741] * v[667] - 2e0*v[2740] * v[673]
			- 2e0*v[2749] * v[677] - v[2750] * v[681] - 2e0*v[2748] * v[686] - 2e0*v[2760] * v[690] - 2e0*v[2759] * v[694]
			- v[2761] * v[699]) / 2e0;
		v[3033] = v[2057] * v[2933] - v[2934] + v[2059] * v[2935] + 24e0*v[1905] * v[4244] * v[4245] - 2e0*v[1252] *
			(v[1753] * v[1906] + v[1732] * v[1907] + 2e0*v[1750] * v[1938] + 2e0*v[1743] * v[1939] + 2e0*v[1738] * v[1940]
				+ 2e0*v[1735] * v[1941] + 2e0*v[1748] * v[1942] + 2e0*v[1737] * v[1943] + v[1741] * v[2052] + 2e0*v[2738]
				- 2e0*v[2741] - 2e0*v[2744] + 2e0*v[2749] + 2e0*v[2752] - 2e0*v[2755] + dB[4] * v[2787] - dB[3] * v[2788]
				- dB[5] * v[2790] - 8e0*v[2054] * v[2937] + 4e0*v[225] * (v[2647] + v[2176] * v[2650] + v[2173] * v[2651]
					+ v[2167] * v[2653] + v[2182] * v[2655] + v[2186] * v[2657] + v[2184] * v[2658] + v[2659] + v[2161] * v[2662]
					+ v[2169] * v[2663] + v[2180] * v[2664] + v[2665] + v[2165] * v[2666] + v[2172] * v[2667] + v[2178] * v[2670]
					+ v[2099] * v[2767] + v[2097] * v[2768] + v[2095] * v[2770] + v[2091] * v[2774] + v[2088] * v[2775] + v[2087] * v[2778]
					+ v[2080] * v[2782] + v[2083] * v[2784] + v[2079] * v[2785] + v[2065] * v[2938] + v[2068] * v[2939] + v[2070] * v[2940]
					+ v[2159] * v[2941] + v[2163] * v[2942] + v[2164] * v[2943] + v[2643] * v[2947] + v[2644] * v[2948] + v[2648] * v[2949]
					+ v[2649] * v[2950] + v[2660] * v[2951] + v[2661] * v[2952] + v[2601] * v[3915] + v[2605] * v[3916] + v[2609] * v[3917]
					- 2e0*v[2074] * v[4248] * v[4249]) - v[2797] * v[437] - v[2798] * v[440] - v[2793] * v[442] + v[2786] * v[662]
				+ 2e0*v[2783] * v[667] + 2e0*v[2781] * v[673] + 2e0*v[2780] * v[677] + v[2777] * v[681] + 2e0*v[2773] * v[686]
				+ 2e0*v[2771] * v[690] + 2e0*v[2769] * v[694] + v[2766] * v[699] + 2e0*v[7568 + i1514]) - 2e0*v[1854] * (
					-4e0*v[1767] * v[1905] + 4e0*v[2930] * v[437] + v[7622 + i1514]);
		v[4012] = v[3033] + (2e0*v[1732] * v[2054] - v[225] * v[2766] + v[2753] * v[659] - v[2756] * v[660] + v[2761] * v[661])
			/ 2e0;
		v[4010] = (v[2791] + v[2795]) / 2e0;
		v[2955] = v[2794] + v[2796];
		v[2956] = v[2789] + v[2792];
		v[2957] = v[992] * (v[2762] * v[3918] + (-(v[262] * v[2758]) + v[2765] * v[994])*(*xfac)) - v[1438] *
			(v[2763] * v[3725] + (v[262] * v[2757] - v[2764] * v[994])*(*zfac));
		v[2958] = v[1439] * (v[1774] * v[2029] + v[1776] * v[2031] + v[1771] * v[2033] + v[1773] * v[2035] + v[1768] * v[2037]
			+ v[1770] * v[2039] - v[1276] * v[2589] - v[1273] * v[2593] - v[1269] * v[2597] + v[185] * v[2799] + v[194] * v[2801]
			+ v[188] * v[2802] + v[197] * v[2804] + v[191] * v[2805] + v[200] * v[2807] + v[210] * (-(v[1528] * v[2028])
				- v[1536] * v[2032] - v[1549] * v[2036] + v[2535] * v[536] + v[2536] * v[537] + v[2537] * v[538] + v[2528] * v[545]
				+ v[2529] * v[546] + v[2530] * v[547] + v[2522] * v[554] + v[2523] * v[555] + v[2524] * v[556]) + v[211] * (-
				(v[1528] * v[2030]) - v[1536] * v[2034] - v[1549] * v[2038] + v[2532] * v[563] + v[2533] * v[564] + v[2534] * v[565]
					+ v[2525] * v[572] + v[2526] * v[573] + v[2527] * v[574] + v[2519] * v[581] + v[2520] * v[582] + v[2521] * v[583])) -
					(v[1774] * v[2028] + v[1776] * v[2030] + v[1771] * v[2032] + v[1773] * v[2034] + v[1768] * v[2036] + v[1770] * v[2038]
						+ v[1275] * v[2589] + v[1272] * v[2593] + v[1268] * v[2597] + v[184] * v[2799] + v[193] * v[2801] + v[187] * v[2802]
						+ v[196] * v[2804] + v[190] * v[2805] + v[199] * v[2807] + v[210] * (v[1528] * v[2029] + v[1536] * v[2033]
							+ v[1549] * v[2037] + v[2535] * v[533] + v[2536] * v[534] + v[2537] * v[535] + v[2528] * v[542] + v[2529] * v[543]
							+ v[2530] * v[544] + v[2522] * v[551] + v[2523] * v[552] + v[2524] * v[553]) + v[211] * (v[1528] * v[2031]
								+ v[1536] * v[2035] + v[1549] * v[2039] + v[2532] * v[560] + v[2533] * v[561] + v[2534] * v[562] + v[2525] * v[569]
								+ v[2526] * v[570] + v[2527] * v[571] + v[2519] * v[578] + v[2520] * v[579] + v[2521] * v[580]))*v[993];
		v[2959] = (-2e0*v[2675] + 2e0*v[2682] - v[2846] * v[491] - 2e0*v[2844] * v[496] - 2e0*v[2845] * v[502]
			- 2e0*v[2826] * v[506] - v[2827] * v[510] - 2e0*v[2828] * v[515] - 2e0*v[2808] * v[519] - 2e0*v[2809] * v[523]
			- v[2810] * v[528]) / 2e0;
		v[2960] = (-2e0*v[1794] * v[1973] + v[165] * v[2883] - v[2846] * v[488] + v[2847] * v[489] - v[2850] * v[490]) / 2e0;
		v[2962] = -v[2676] + v[2701] + (v[2847] * v[491]) / 2e0 + v[2848] * v[496] + v[2849] * v[502] + v[2829] * v[506] +
			(v[2831] * v[510]) / 2e0 + v[2830] * v[515] + v[2812] * v[519] + v[2811] * v[523] + (v[2813] * v[528]) / 2e0;
		v[2963] = (-2e0*v[1786] * v[1973] + v[165] * v[2874] - v[2827] * v[488] + v[2831] * v[489] - v[2834] * v[490]) / 2e0;
		v[2964] = (-2e0*v[2686] + 2e0*v[2698] - v[2850] * v[491] - 2e0*v[2852] * v[496] - 2e0*v[2851] * v[502]
			- 2e0*v[2833] * v[506] - v[2834] * v[510] - 2e0*v[2832] * v[515] - 2e0*v[2815] * v[519] - 2e0*v[2814] * v[523]
			- v[2816] * v[528]) / 2e0;
		v[3025] = v[1976] * v[2962] - v[2963] + v[1978] * v[2964] + 24e0*v[1901] * v[4250] * v[4251] - 2e0*v[1238] *
			(v[1794] * v[1902] + v[1777] * v[1903] + 2e0*v[1791] * v[1932] + 2e0*v[1788] * v[1933] + 2e0*v[1783] * v[1934]
				+ 2e0*v[1780] * v[1935] + 2e0*v[1789] * v[1936] + 2e0*v[1782] * v[1937] + v[1786] * v[1971] + 2e0*v[2809]
				- 2e0*v[2812] - 2e0*v[2828] + 2e0*v[2833] + 2e0*v[2849] - 2e0*v[2852] + dA[10] * v[2884] - dA[9] * v[2885]
				- dA[11] * v[2887] - 8e0*v[1973] * v[2966] + 4e0*v[165] * (v[2677] + v[2148] * v[2680] + v[2145] * v[2681]
					+ v[2139] * v[2683] + v[2154] * v[2685] + v[2158] * v[2687] + v[2156] * v[2688] + v[2689] + v[2133] * v[2692]
					+ v[2141] * v[2693] + v[2152] * v[2694] + v[2695] + v[2137] * v[2696] + v[2144] * v[2697] + v[2150] * v[2700]
					+ v[2027] * v[2864] + v[2025] * v[2865] + v[2023] * v[2867] + v[2022] * v[2871] + v[2019] * v[2872] + v[2018] * v[2875]
					+ v[2014] * v[2879] + v[2017] * v[2881] + v[2013] * v[2882] + v[1984] * v[2967] + v[1987] * v[2968] + v[1989] * v[2969]
					+ v[2131] * v[2970] + v[2135] * v[2971] + v[2136] * v[2972] + v[2673] * v[2976] + v[2674] * v[2977] + v[2678] * v[2978]
					+ v[2679] * v[2979] + v[2690] * v[2980] + v[2691] * v[2981] + v[2613] * v[3929] + v[2617] * v[3930] + v[2621] * v[3931]
					- 2e0*v[1993] * v[4254] * v[4255]) - v[2894] * v[431] - v[2895] * v[434] - v[2890] * v[436] + v[2883] * v[491]
				+ 2e0*v[2880] * v[496] + 2e0*v[2878] * v[502] + 2e0*v[2877] * v[506] + v[2874] * v[510] + 2e0*v[2870] * v[515]
				+ 2e0*v[2868] * v[519] + 2e0*v[2866] * v[523] + v[2863] * v[528] + 2e0*v[7586 + i1514]) - 2e0*v[1852] * (
					-4e0*v[1826] * v[1901] + 4e0*v[2959] * v[431] + v[7640 + i1514]);
		v[4008] = v[3025] + (2e0*v[1777] * v[1973] - v[165] * v[2863] + v[2810] * v[488] - v[2813] * v[489] + v[2816] * v[490])
			/ 2e0;
		v[4006] = (v[2888] + v[2892]) / 2e0;
		v[2984] = v[2891] + v[2893];
		v[2985] = v[2886] + v[2889];
		v[2986] = (-2e0*v[2705] + 2e0*v[2712] - v[2855] * v[446] - 2e0*v[2853] * v[451] - 2e0*v[2854] * v[457]
			- 2e0*v[2835] * v[461] - v[2836] * v[465] - 2e0*v[2837] * v[470] - 2e0*v[2817] * v[474] - 2e0*v[2818] * v[478]
			- v[2819] * v[483]) / 2e0;
		v[2987] = (-2e0*v[1812] * v[1946] + v[146] * v[2916] - v[2855] * v[443] + v[2856] * v[444] - v[2859] * v[445]) / 2e0;
		v[2989] = -v[2706] + v[2731] + (v[2856] * v[446]) / 2e0 + v[2857] * v[451] + v[2858] * v[457] + v[2838] * v[461] +
			(v[2840] * v[465]) / 2e0 + v[2839] * v[470] + v[2821] * v[474] + v[2820] * v[478] + (v[2822] * v[483]) / 2e0;
		v[2990] = (-2e0*v[1804] * v[1946] + v[146] * v[2907] - v[2836] * v[443] + v[2840] * v[444] - v[2843] * v[445]) / 2e0;
		v[2991] = (-2e0*v[2716] + 2e0*v[2728] - v[2859] * v[446] - 2e0*v[2861] * v[451] - 2e0*v[2860] * v[457]
			- 2e0*v[2842] * v[461] - v[2843] * v[465] - 2e0*v[2841] * v[470] - 2e0*v[2824] * v[474] - 2e0*v[2823] * v[478]
			- v[2825] * v[483]) / 2e0;
		v[3017] = v[1949] * v[2989] - v[2990] + v[1951] * v[2991] + 24e0*v[1897] * v[4256] * v[4257] - 2e0*v[1224] *
			(v[1812] * v[1898] + v[1795] * v[1899] + 2e0*v[1809] * v[1926] + 2e0*v[1806] * v[1927] + 2e0*v[1801] * v[1928]
				+ 2e0*v[1798] * v[1929] + 2e0*v[1807] * v[1930] + 2e0*v[1800] * v[1931] + v[1804] * v[1944] + 2e0*v[2818]
				- 2e0*v[2821] - 2e0*v[2837] + 2e0*v[2842] + 2e0*v[2858] - 2e0*v[2861] + dA[4] * v[2917] - dA[3] * v[2918]
				- dA[5] * v[2920] - 8e0*v[1946] * v[2993] - v[2927] * v[425] + 4e0*v[146] * (v[2707] + v[2120] * v[2710]
					+ v[2117] * v[2711] + v[2111] * v[2713] + v[2126] * v[2715] + v[2130] * v[2717] + v[2128] * v[2718] + v[2719]
					+ v[2105] * v[2722] + v[2113] * v[2723] + v[2124] * v[2724] + v[2725] + v[2109] * v[2726] + v[2116] * v[2727]
					+ v[2122] * v[2730] + v[2012] * v[2897] + v[2010] * v[2898] + v[2008] * v[2900] + v[2007] * v[2904] + v[2004] * v[2905]
					+ v[2003] * v[2908] + v[1999] * v[2912] + v[2002] * v[2914] + v[1998] * v[2915] + v[1957] * v[2994] + v[1960] * v[2995]
					+ v[1962] * v[2996] + v[2103] * v[2997] + v[2107] * v[2998] + v[2108] * v[2999] + v[2703] * v[3003] + v[2704] * v[3004]
					+ v[2708] * v[3005] + v[2709] * v[3006] + v[2720] * v[3007] + v[2721] * v[3008] + v[2625] * v[3942] + v[2629] * v[3943]
					+ v[2633] * v[3944] - 2e0*v[1966] * v[4260] * v[4261]) - v[2928] * v[428] - v[2923] * v[430] + v[2916] * v[446]
				+ 2e0*v[2913] * v[451] + 2e0*v[2911] * v[457] + 2e0*v[2910] * v[461] + v[2907] * v[465] + 2e0*v[2903] * v[470]
				+ 2e0*v[2901] * v[474] + 2e0*v[2899] * v[478] + v[2896] * v[483] + 2e0*v[7604 + i1514]) - 2e0*v[1850] * (
					-4e0*v[1840] * v[1897] + 4e0*v[2986] * v[425] + v[7676 + i1514]);
		v[4004] = v[3017] + (2e0*v[1795] * v[1946] - v[146] * v[2896] + v[2819] * v[443] - v[2822] * v[444] + v[2825] * v[445])
			/ 2e0;
		v[4002] = (v[2921] + v[2925]) / 2e0;
		v[3011] = v[2924] + v[2926];
		v[3012] = v[2919] + v[2922];
		v[7731] = 0e0;
		v[7732] = 0e0;
		v[7733] = 0e0;
		v[7734] = 2e0*v[3013];
		v[7735] = v[3014];
		v[7736] = v[3015];
		v[7737] = 0e0;
		v[7738] = 0e0;
		v[7739] = 0e0;
		v[7740] = 0e0;
		v[7741] = 0e0;
		v[7742] = 0e0;
		v[7743] = 0e0;
		v[7744] = 0e0;
		v[7745] = 0e0;
		v[7746] = 0e0;
		v[7747] = 0e0;
		v[7748] = 0e0;
		v[7713] = 0e0;
		v[7714] = 0e0;
		v[7715] = 0e0;
		v[7716] = v[3014];
		v[7717] = 2e0*v[3018];
		v[7718] = v[3019];
		v[7719] = 0e0;
		v[7720] = 0e0;
		v[7721] = 0e0;
		v[7722] = 0e0;
		v[7723] = 0e0;
		v[7724] = 0e0;
		v[7725] = 0e0;
		v[7726] = 0e0;
		v[7727] = 0e0;
		v[7728] = 0e0;
		v[7729] = 0e0;
		v[7730] = 0e0;
		v[7695] = 0e0;
		v[7696] = 0e0;
		v[7697] = 0e0;
		v[7698] = v[3015];
		v[7699] = v[3019];
		v[7700] = 2e0*v[3020];
		v[7701] = 0e0;
		v[7702] = 0e0;
		v[7703] = 0e0;
		v[7704] = 0e0;
		v[7705] = 0e0;
		v[7706] = 0e0;
		v[7707] = 0e0;
		v[7708] = 0e0;
		v[7709] = 0e0;
		v[7710] = 0e0;
		v[7711] = 0e0;
		v[7712] = 0e0;
		v[7659] = 0e0;
		v[7660] = 0e0;
		v[7661] = 0e0;
		v[7662] = 0e0;
		v[7663] = 0e0;
		v[7664] = 0e0;
		v[7665] = 0e0;
		v[7666] = 0e0;
		v[7667] = 0e0;
		v[7668] = 2e0*v[3021];
		v[7669] = v[3022];
		v[7670] = v[3023];
		v[7671] = 0e0;
		v[7672] = 0e0;
		v[7673] = 0e0;
		v[7674] = 0e0;
		v[7675] = 0e0;
		v[7676] = 0e0;
		v[7749] = 0e0;
		v[7750] = 0e0;
		v[7751] = 0e0;
		v[7752] = 0e0;
		v[7753] = 0e0;
		v[7754] = 0e0;
		v[7755] = 0e0;
		v[7756] = 0e0;
		v[7757] = 0e0;
		v[7758] = v[3022];
		v[7759] = 2e0*v[3026];
		v[7760] = v[3027];
		v[7761] = 0e0;
		v[7762] = 0e0;
		v[7763] = 0e0;
		v[7764] = 0e0;
		v[7765] = 0e0;
		v[7766] = 0e0;
		v[7767] = 0e0;
		v[7768] = 0e0;
		v[7769] = 0e0;
		v[7770] = 0e0;
		v[7771] = 0e0;
		v[7772] = 0e0;
		v[7773] = 0e0;
		v[7774] = 0e0;
		v[7775] = 0e0;
		v[7776] = v[3023];
		v[7777] = v[3027];
		v[7778] = 2e0*v[3028];
		v[7779] = 0e0;
		v[7780] = 0e0;
		v[7781] = 0e0;
		v[7782] = 0e0;
		v[7783] = 0e0;
		v[7784] = 0e0;
		v[7785] = 0e0;
		v[7786] = 0e0;
		v[7787] = 0e0;
		v[7788] = 0e0;
		v[7789] = 0e0;
		v[7790] = 0e0;
		v[7791] = 0e0;
		v[7792] = 0e0;
		v[7793] = 0e0;
		v[7794] = 0e0;
		v[7795] = 0e0;
		v[7796] = 0e0;
		v[7797] = 0e0;
		v[7798] = 0e0;
		v[7799] = 0e0;
		v[7800] = 2e0*v[3029];
		v[7801] = v[3030];
		v[7802] = v[3031];
		v[7803] = 0e0;
		v[7804] = 0e0;
		v[7805] = 0e0;
		v[7806] = 0e0;
		v[7807] = 0e0;
		v[7808] = 0e0;
		v[7809] = 0e0;
		v[7810] = 0e0;
		v[7811] = 0e0;
		v[7812] = 0e0;
		v[7813] = 0e0;
		v[7814] = 0e0;
		v[7815] = 0e0;
		v[7816] = 0e0;
		v[7817] = 0e0;
		v[7818] = v[3030];
		v[7819] = 2e0*v[3034];
		v[7820] = v[3035];
		v[7821] = 0e0;
		v[7822] = 0e0;
		v[7823] = 0e0;
		v[7824] = 0e0;
		v[7825] = 0e0;
		v[7826] = 0e0;
		v[7827] = 0e0;
		v[7828] = 0e0;
		v[7829] = 0e0;
		v[7830] = 0e0;
		v[7831] = 0e0;
		v[7832] = 0e0;
		v[7833] = 0e0;
		v[7834] = 0e0;
		v[7835] = 0e0;
		v[7836] = v[3031];
		v[7837] = v[3035];
		v[7838] = 2e0*v[3036];
		v[7839] = v[2553];
		v[7840] = v[2563];
		v[7841] = v[2573];
		v[7842] = v[1827] * v[1946] - v[2924] + v[2926] + v[2918] * v[4001] + dA[5] * v[4002] + 2e0*(v[1716] * v[3016]
			+ v[2927] * v[4001] + v[2986] * v[4003]) + 2e0*dA[3] * v[4004] + v[3012] * v[426] + v[7730 + i1514];
		v[7843] = -(v[1828] * v[1946]) + v[2921] - v[2925] + (dA[5] * v[3011]) / 2e0 + (dA[3] * v[3012]) / 2e0 - v[2917] * v[4001]
			+ 2e0*(-4e0*v[1224] * v[2989] - v[1731] * v[3016] + v[2928] * v[4001]) + 2e0*dA[4] * (-v[2987] + v[2990] + v[4004])
			+ v[7712 + i1514];
		v[7844] = v[1829] * v[1946] - v[2919] + v[2922] + 2e0*dA[5] * (-v[2987] + v[3017]) + v[2920] * v[4001] + dA[3] * v[4002]
			+ 2e0*(v[1727] * v[3016] + v[2923] * v[4001] + v[2991] * v[4003]) + v[3011] * v[426] + v[7694 + i1514];
		v[7845] = v[2552];
		v[7846] = v[2562];
		v[7847] = v[2572];
		v[7848] = v[1813] * v[1973] - v[2891] + v[2893] + v[2885] * v[4005] + dA[11] * v[4006] + 2e0*(v[1689] * v[3024]
			+ v[2894] * v[4005] + v[2959] * v[4007]) + 2e0*dA[9] * v[4008] + v[2985] * v[432] + v[7658 + i1514];
		v[7849] = -(v[1814] * v[1973]) + v[2888] - v[2892] + (dA[11] * v[2984]) / 2e0 + (dA[9] * v[2985]) / 2e0
			- v[2884] * v[4005] + 2e0*(-4e0*v[1238] * v[2962] - v[1704] * v[3024] + v[2895] * v[4005]) + 2e0*dA[10] * (-v[2960]
				+ v[2963] + v[4008]) + v[7748 + i1514];
		v[7850] = v[1815] * v[1973] - v[2886] + v[2889] + 2e0*dA[11] * (-v[2960] + v[3025]) + v[2887] * v[4005]
			+ dA[9] * v[4006] + 2e0*(v[1700] * v[3024] + v[2890] * v[4005] + v[2964] * v[4007]) + v[2984] * v[432] + v[7766
			+ i1514];
		v[7851] = v[2548];
		v[7852] = v[2558];
		v[7853] = v[2568];
		v[7854] = v[1754] * v[2054] - v[2794] + v[2796] + v[2788] * v[4009] + dB[5] * v[4010] + 2e0*(v[1662] * v[3032]
			+ v[2797] * v[4009] + v[2930] * v[4011]) + 2e0*dB[3] * v[4012] + v[2956] * v[438] + v[7784 + i1514];
		v[7855] = -(v[1755] * v[2054]) + v[2791] - v[2795] + (dB[5] * v[2955]) / 2e0 + (dB[3] * v[2956]) / 2e0 - v[2787] * v[4009]
			+ 2e0*(-4e0*v[1252] * v[2933] - v[1677] * v[3032] + v[2798] * v[4009]) + 2e0*dB[4] * (-v[2931] + v[2934] + v[4012])
			+ v[7802 + i1514];
		v[7856] = v[1756] * v[2054] - v[2789] + v[2792] + 2e0*dB[5] * (-v[2931] + v[3033]) + v[2790] * v[4009] + dB[3] * v[4010]
			+ 2e0*(v[1673] * v[3032] + v[2793] * v[4009] + v[2935] * v[4011]) + v[2955] * v[438] + v[7820 + i1514];
		Rc[i1514 - 1] += v[5908 + i1514] + (*a4)*v[5926 + i1514];
		for (i1876 = 1; i1876 <= 18; i1876++) {
			Kc[i1514 - 1][i1876 - 1] += v[2929] * v[4838 + i1876] + v[2862] * v[4856 + i1876] + v[2958] * v[4874 + i1876]
				+ v[2957] * v[4892 + i1876] + v[7838 + i1876] + (*a4)*v[7856 + i1876];
		};/* end for */
	};/* end for */
	v[3051] = v[1003] * v[949];
	v[3052] = v[1003] * v[950];
	v[3053] = v[1003] * v[951];
	v[3054] = v[1003] * v[952];
	v[3055] = v[1003] * v[953];
	v[3056] = v[1003] * v[954];
	v[3057] = v[1003] * v[955];
	v[3058] = v[1003] * v[956];
	v[3059] = v[1003] * v[957];
	v[3060] = v[1002] * v[949];
	v[3061] = v[1002] * v[950];
	v[3062] = v[1002] * v[951];
	v[3063] = v[1002] * v[952];
	v[3064] = v[1002] * v[953];
	v[3065] = v[1002] * v[954];
	v[3066] = v[1002] * v[955];
	v[3067] = v[1002] * v[956];
	v[3068] = v[1002] * v[957];
	v[3069] = -(v[1001] * v[246]) - v[1002] * v[249] - v[1003] * v[252];
	v[3070] = -(v[1001] * v[244]) - v[1002] * v[247] - v[1003] * v[250];
	v[3563] = -(v[3069] * v[3633]) + v[3070] * v[3704];
	v[3071] = v[1001] * v[949];
	v[3072] = v[1001] * v[950];
	v[3073] = v[1001] * v[951];
	v[3074] = v[1001] * v[952];
	v[3075] = v[1001] * v[953];
	v[3076] = v[1001] * v[954];
	v[3077] = v[1001] * v[955];
	v[3078] = v[1001] * v[956];
	v[3079] = v[1001] * v[957];
	v[3080] = v[3070] * v[410] + v[3069] * v[412];
	v[3081] = (v[225] * v[3051]) / 2e0;
	v[3082] = v[225] * v[3052];
	v[3083] = v[225] * v[3053];
	v[3084] = v[225] * v[3060];
	v[3085] = (v[225] * v[3061]) / 2e0;
	v[3086] = v[225] * v[3062];
	v[3087] = v[225] * v[3071];
	v[3088] = v[225] * v[3072];
	v[3089] = (v[3073] * v[662]) / 2e0 + v[3072] * v[667] + v[3071] * v[673] + v[3062] * v[677] + (v[3061] * v[681]) / 2e0
		+ v[3060] * v[686] + v[3053] * v[690] + v[3052] * v[694] + (v[3051] * v[699]) / 2e0;
	v[3614] = -v[3081] - v[3085] - 4e0*v[1252] * v[3089];
	v[3620] = -(v[225] * v[3073]) / 2e0 + v[3081] + v[3614];
	v[3618] = -v[3081] + v[3085] + v[3620];
	v[3091] = v[1439] * v[3577] - v[3578] * v[993];
	v[3092] = (-(v[1001] * v[387]) + v[1001] * v[388] - v[1002] * v[390] + v[1002] * v[391] - v[1003] * v[393]
		+ v[1003] * v[394]) / 2e0;
	v[3093] = (v[165] * v[3054]) / 2e0;
	v[3094] = v[165] * v[3055];
	v[3095] = v[165] * v[3056];
	v[3096] = v[165] * v[3063];
	v[3097] = (v[165] * v[3064]) / 2e0;
	v[3098] = v[165] * v[3065];
	v[3099] = v[165] * v[3074];
	v[3100] = v[165] * v[3075];
	v[3101] = (v[3076] * v[491]) / 2e0 + v[3075] * v[496] + v[3074] * v[502] + v[3065] * v[506] + (v[3064] * v[510]) / 2e0
		+ v[3063] * v[515] + v[3056] * v[519] + v[3055] * v[523] + (v[3054] * v[528]) / 2e0;
	v[3607] = -v[3093] - v[3097] - 4e0*v[1238] * v[3101];
	v[3613] = -(v[165] * v[3076]) / 2e0 + v[3093] + v[3607];
	v[3611] = -v[3093] + v[3097] + v[3613];
	v[3103] = (v[146] * v[3057]) / 2e0;
	v[3104] = v[146] * v[3058];
	v[3105] = v[146] * v[3059];
	v[3106] = v[146] * v[3066];
	v[3107] = (v[146] * v[3067]) / 2e0;
	v[3108] = v[146] * v[3068];
	v[3109] = v[146] * v[3077];
	v[3110] = v[146] * v[3078];
	v[3111] = (v[3079] * v[446]) / 2e0 + v[3078] * v[451] + v[3077] * v[457] + v[3068] * v[461] + (v[3067] * v[465]) / 2e0
		+ v[3066] * v[470] + v[3059] * v[474] + v[3058] * v[478] + (v[3057] * v[483]) / 2e0;
	v[3600] = -v[3103] - v[3107] - 4e0*v[1224] * v[3111];
	v[3606] = -(v[146] * v[3079]) / 2e0 + v[3103] + v[3600];
	v[3604] = -v[3103] + v[3107] + v[3606];
	v[3113] = v[1440] * v[3564] - v[3563] * v[994];
	v[3616] = (v[3083] + v[3087]) / 2e0;
	v[3115] = v[3082] + v[3084];
	v[3619] = v[3115] / 2e0;
	v[3116] = v[3086] + v[3088];
	v[3615] = v[3116] / 2e0;
	v[3609] = (v[3095] + v[3099]) / 2e0;
	v[3118] = v[3094] + v[3096];
	v[3612] = v[3118] / 2e0;
	v[3119] = v[3098] + v[3100];
	v[3608] = v[3119] / 2e0;
	v[3602] = (v[3105] + v[3109]) / 2e0;
	v[3121] = v[3104] + v[3106];
	v[3605] = v[3121] / 2e0;
	v[3122] = v[3108] + v[3110];
	v[7879] = v[1001] * v[210];
	v[7880] = v[1002] * v[210];
	v[7881] = v[1003] * v[210];
	v[7882] = v[3104] - v[3106] + 2e0*dA[3] * v[3600] + dA[5] * v[3602] + v[3122] * v[426];
	v[7883] = -v[3105] + v[3109] + (dA[5] * v[3121]) / 2e0 + (dA[3] * v[3122]) / 2e0 + 2e0*dA[4] * v[3604];
	v[7884] = v[3108] - v[3110] + dA[3] * v[3602] + 2e0*dA[5] * v[3606] + v[3121] * v[426];
	v[7885] = v[1001] * v[211];
	v[7886] = v[1002] * v[211];
	v[7887] = v[1003] * v[211];
	v[7888] = v[3094] - v[3096] + 2e0*dA[9] * v[3607] + dA[11] * v[3609] + v[3119] * v[432];
	v[7889] = -v[3095] + v[3099] + (dA[11] * v[3118]) / 2e0 + (dA[9] * v[3119]) / 2e0 + 2e0*dA[10] * v[3611];
	v[7890] = v[3098] - v[3100] + dA[9] * v[3609] + 2e0*dA[11] * v[3613] + v[3118] * v[432];
	v[7891] = -v[1001];
	v[7892] = -v[1002];
	v[7893] = -v[1003];
	v[7894] = v[3082] - v[3084] + 2e0*dB[3] * v[3614] + dB[5] * v[3616] + v[3116] * v[438];
	v[7895] = -v[3083] + v[3087] + (dB[5] * v[3115]) / 2e0 + (dB[3] * v[3116]) / 2e0 + 2e0*dB[4] * v[3618];
	v[7896] = v[3086] - v[3088] + dB[3] * v[3616] + 2e0*dB[5] * v[3620] + v[3115] * v[438];
	v[3601] = v[3122] / 2e0;
	for (i3049 = 1; i3049 <= 18; i3049++) {
		i4024 = (i3049 == 17 ? 1 : 0);
		i4023 = (i3049 == 16 ? 1 : 0);
		i4022 = (i3049 == 18 ? 1 : 0);
		i4021 = (i3049 == 11 ? 1 : 0);
		i4020 = (i3049 == 10 ? 1 : 0);
		i4019 = (i3049 == 12 ? 1 : 0);
		i4018 = (i3049 == 5 ? 1 : 0);
		i4017 = (i3049 == 4 ? 1 : 0);
		i4016 = (i3049 == 6 ? 1 : 0);
		v[3223] = v[4874 + i3049];
		v[4015] = v[1001] * v[3223];
		v[4014] = v[1002] * v[3223];
		v[4013] = v[1003] * v[3223];
		v[3442] = v[211] * v[4013];
		v[3439] = v[210] * v[4013];
		v[3436] = v[211] * v[4014];
		v[3433] = v[210] * v[4014];
		v[3430] = v[211] * v[4015];
		v[3427] = v[210] * v[4015];
		v[3222] = v[4856 + i3049];
		v[3319] = -(v[1001] * v[3222]) / 2e0;
		v[3313] = -(v[1002] * v[3222]) / 2e0;
		v[3307] = -(v[1003] * v[3222]) / 2e0;
		v[3158] = v[4838 + i3049];
		v[4059] = v[3158] * v[994];
		v[4038] = -(v[1440] * v[3158]);
		v[3157] = v[4892 + i3049];
		v[4060] = v[3069] * v[3157];
		v[4058] = v[3070] * v[3157];
		v[3130] = v[4950 + i3049];
		v[3131] = v[4986 + i3049];
		v[3132] = v[4968 + i3049];
		v[3133] = v[5022 + i3049];
		v[3134] = v[5058 + i3049];
		v[3135] = v[5040 + i3049];
		v[3136] = v[5094 + i3049];
		v[3137] = v[5130 + i3049];
		v[3138] = v[5112 + i3049];
		v[3140] = v[5948 + i3049];
		v[3141] = v[4932 + i3049];
		v[3191] = -4e0*v[1224] * v[3141];
		v[3142] = v[7900 + i3049];
		v[3144] = v[5984 + i3049];
		v[3146] = v[6056 + i3049];
		v[3147] = v[5004 + i3049];
		v[3202] = -4e0*v[1238] * v[3147];
		v[3148] = v[7918 + i3049];
		v[3150] = v[6092 + i3049];
		v[3152] = v[6272 + i3049];
		v[3153] = v[5076 + i3049];
		v[3213] = -4e0*v[1252] * v[3153];
		v[3154] = v[7936 + i3049];
		v[3156] = v[6308 + i3049];
		v[3159] = v[3157] * v[410] + v[3158] * v[418];
		v[3160] = v[3157] * v[412] + v[3158] * v[420];
		v[3161] = -i4016 + v[3130];
		v[3163] = i4016 + v[3130];
		v[3164] = -i4017 + v[3131];
		v[3166] = i4017 + v[3131];
		v[3167] = i4018 + v[3132];
		v[3169] = -i4018 + v[3132];
		v[3170] = -i4019 + v[3133];
		v[3172] = i4019 + v[3133];
		v[3173] = -i4020 + v[3134];
		v[3175] = i4020 + v[3134];
		v[3176] = i4021 + v[3135];
		v[3178] = -i4021 + v[3135];
		v[3179] = -i4022 + v[3136];
		v[3181] = i4022 + v[3136];
		v[3182] = -i4023 + v[3137];
		v[3184] = i4023 + v[3137];
		v[3185] = i4024 + v[3138];
		v[3187] = -i4024 + v[3138];
		v[3189] = (v[146] * v[3140]) / 2e0 + v[3141] * v[3188];
		v[3190] = v[146] * v[3161] + v[3191] * v[451];
		v[3192] = v[146] * v[3167] + v[3191] * v[457];
		v[3193] = v[146] * v[3163] + v[3191] * v[461];
		v[3194] = (v[146] * v[3142] + v[3191] * v[465]) / 2e0;
		v[3195] = v[146] * v[3164] + v[3191] * v[470];
		v[3196] = v[146] * v[3169] + v[3191] * v[474];
		v[3197] = v[146] * v[3166] + v[3191] * v[478];
		v[3198] = (v[146] * v[3144] + v[3191] * v[483]) / 2e0;
		v[3200] = (v[165] * v[3146]) / 2e0 + v[3147] * v[3199];
		v[3201] = v[165] * v[3170] + v[3202] * v[496];
		v[3203] = v[165] * v[3176] + v[3202] * v[502];
		v[3204] = v[165] * v[3172] + v[3202] * v[506];
		v[3205] = (v[165] * v[3148] + v[3202] * v[510]) / 2e0;
		v[3206] = v[165] * v[3173] + v[3202] * v[515];
		v[3207] = v[165] * v[3178] + v[3202] * v[519];
		v[3208] = v[165] * v[3175] + v[3202] * v[523];
		v[3209] = (v[165] * v[3150] + v[3202] * v[528]) / 2e0;
		v[3211] = (v[225] * v[3152]) / 2e0 + v[3153] * v[3210];
		v[3212] = v[225] * v[3179] + v[3213] * v[667];
		v[3214] = v[225] * v[3185] + v[3213] * v[673];
		v[3215] = v[225] * v[3181] + v[3213] * v[677];
		v[3216] = (v[225] * v[3154] + v[3213] * v[681]) / 2e0;
		v[3217] = v[225] * v[3182] + v[3213] * v[686];
		v[3218] = v[225] * v[3187] + v[3213] * v[690];
		v[3219] = v[225] * v[3184] + v[3213] * v[694];
		v[3220] = (v[225] * v[3156] + v[3213] * v[699]) / 2e0;
		v[3224] = -(v[244] * v[3159]) - v[246] * v[3160] + v[3222] * v[389] - v[3158] * v[4025] + v[3223] * v[405] + v[5670
			+ i3049] + v[3214] * v[949] + v[3212] * v[950] + v[3211] * v[951] + v[3203] * v[952] + v[3201] * v[953] + v[3200] * v[954]
			+ v[3192] * v[955] + v[3190] * v[956] + v[3189] * v[957];
		v[4062] = (*cn)*v[3224];
		v[4031] = v[3224] * v[376];
		v[3264] = -(v[3224] * v[3651]);
		v[3251] = -(v[4062] * v[982]);
		v[3226] = -(v[247] * v[3159]) - v[249] * v[3160] + v[3222] * v[392] - v[3158] * v[4026] + v[3223] * v[406] + v[5652
			+ i3049] + v[3217] * v[949] + v[3216] * v[950] + v[3215] * v[951] + v[3206] * v[952] + v[3205] * v[953] + v[3204] * v[954]
			+ v[3195] * v[955] + v[3194] * v[956] + v[3193] * v[957];
		v[4027] = (*cn)*v[3226];
		v[3259] = v[376] * v[4027];
		v[3255] = -(v[4027] * v[986]);
		v[3227] = v[1001] * v[3189] + v[1002] * v[3193] + v[1003] * v[3196];
		v[3228] = v[1001] * v[3190] + v[1002] * v[3194] + v[1003] * v[3197];
		v[3229] = v[1001] * v[3192] + v[1002] * v[3195] + v[1003] * v[3198];
		v[3230] = v[1001] * v[3200] + v[1002] * v[3204] + v[1003] * v[3207];
		v[3231] = v[1001] * v[3201] + v[1002] * v[3205] + v[1003] * v[3208];
		v[3232] = v[1001] * v[3203] + v[1002] * v[3206] + v[1003] * v[3209];
		v[3233] = v[1001] * v[3211] + v[1002] * v[3215] + v[1003] * v[3218];
		v[3234] = v[1001] * v[3212] + v[1002] * v[3216] + v[1003] * v[3219];
		v[3236] = -(v[250] * v[3159]) - v[252] * v[3160] + v[3222] * v[395] - v[3158] * v[4028] + v[3223] * v[407] + v[5634
			+ i3049] + v[3220] * v[949] + v[3219] * v[950] + v[3218] * v[951] + v[3209] * v[952] + v[3208] * v[953] + v[3207] * v[954]
			+ v[3198] * v[955] + v[3197] * v[956] + v[3196] * v[957];
		v[4030] = v[3226] * v[376] + v[3236] * v[377];
		v[4029] = (*cn)*v[3236];
		v[3277] = v[377] * v[4029];
		v[4061] = -v[3251] + v[3277] * v[375];
		v[4035] = -v[3251] + (v[3259] + v[3277])*v[375];
		v[3261] = -(v[4029] * v[990]);
		v[4032] = v[3255] - v[3277] * v[376];
		v[4034] = -(v[3264] * v[376]) + v[4032];
		v[3237] = v[1001] * v[3214] + v[1002] * v[3217] + v[1003] * v[3220];
		v[3238] = v[3259] + v[3264];
		v[4051] = v[3238] * v[377];
		v[4033] = -v[3261] + v[4051];
		v[3239] = (*cn)*(v[274] * v[3224] + v[273] * v[3226]);
		v[3240] = (*cn)*(v[277] * v[3224] + v[276] * v[3226]);
		v[4036] = v[210] * v[3239] + v[211] * v[3240];
		v[3241] = (*cn)*(v[2598] * v[3224] + v[2599] * v[3226] + v[2600] * v[3236]);
		v[3410] = v[1873] * v[3241];
		v[3396] = v[3045] * v[3241];
		v[3242] = (*cn)*(v[2602] * v[3224] + v[2603] * v[3226] + v[2604] * v[3236]);
		v[3403] = v[3044] * v[3242];
		v[4049] = v[1872] * v[3241] + v[3403];
		v[4048] = v[1870] * v[3242] + v[3396];
		v[3243] = (*cn)*(v[2606] * v[3224] + v[2607] * v[3226] + v[2608] * v[3236]);
		v[3408] = v[3043] * v[3243];
		v[4050] = v[1871] * v[3242] + v[3408];
		v[3405] = v[3243] * v[3645];
		v[3399] = v[3243] * v[3647];
		v[3244] = (*cn)*(v[2610] * v[3224] + v[2611] * v[3226] + v[2612] * v[3236]);
		v[3516] = v[1867] * v[3244];
		v[3502] = v[3042] * v[3244];
		v[3245] = (*cn)*(v[2614] * v[3224] + v[2615] * v[3226] + v[2616] * v[3236]);
		v[3509] = v[3041] * v[3245];
		v[4053] = v[1866] * v[3244] + v[3509];
		v[4052] = v[1864] * v[3245] + v[3502];
		v[3246] = (*cn)*(v[2618] * v[3224] + v[2619] * v[3226] + v[2620] * v[3236]);
		v[3514] = v[3040] * v[3246];
		v[4054] = v[1865] * v[3245] + v[3514];
		v[3511] = v[3246] * v[3640];
		v[3505] = v[3246] * v[3642];
		v[3247] = (*cn)*(v[2622] * v[3224] + v[2623] * v[3226] + v[2624] * v[3236]);
		v[3534] = v[1861] * v[3247];
		v[3520] = v[3039] * v[3247];
		v[3248] = (*cn)*(v[2626] * v[3224] + v[2627] * v[3226] + v[2628] * v[3236]);
		v[3527] = v[3038] * v[3248];
		v[4056] = v[1860] * v[3247] + v[3527];
		v[4055] = v[1858] * v[3248] + v[3520];
		v[3249] = (*cn)*(v[2630] * v[3224] + v[2631] * v[3226] + v[2632] * v[3236]);
		v[3532] = v[3037] * v[3249];
		v[8171] = v[210] * v[4061] + v[4027] * v[984];
		v[8172] = -(v[210] * v[4032]) + v[4062] * v[984];
		v[8173] = v[210] * v[4033];
		v[8174] = v[3532] + v[3248] * v[3635] + v[3247] * v[3637];
		v[8175] = v[1858] * v[3247] + v[1859] * v[3249] + v[3527];
		v[8176] = v[1860] * v[3248] + v[1861] * v[3249] + v[3520];
		v[8177] = v[211] * v[4061] + v[4027] * v[985];
		v[8178] = -(v[211] * v[4032]) + v[4062] * v[985];
		v[8179] = v[211] * v[4033];
		v[8180] = v[3514] + v[3245] * v[3640] + v[3244] * v[3642];
		v[8181] = v[1864] * v[3244] + v[1865] * v[3246] + v[3509];
		v[8182] = v[1866] * v[3245] + v[1867] * v[3246] + v[3502];
		v[8183] = -v[4035];
		v[8184] = v[4034];
		v[8185] = -v[4033];
		v[8186] = v[3408] + v[3242] * v[3645] + v[3241] * v[3647];
		v[8187] = v[1870] * v[3241] + v[1871] * v[3243] + v[3403];
		v[8188] = v[1872] * v[3242] + v[1873] * v[3243] + v[3396];
		v[4057] = v[1859] * v[3248] + v[3532];
		v[3529] = v[3249] * v[3635];
		v[3523] = v[3249] * v[3637];
		v[3250] = v[3251] * v[339] + v[2501] * v[4030];
		v[3252] = v[3251] * v[349] + v[2504] * v[4030];
		v[3253] = v[3251] * v[359] + v[2506] * v[4030];
		v[3254] = v[2501] * v[4031] + v[339] * v[4032];
		v[3257] = v[2504] * v[4031] + v[349] * v[4032];
		v[3258] = v[2506] * v[4031] + v[359] * v[4032];
		v[3260] = v[3261] * v[339] + (v[2501] * v[3224] - v[3259] * v[339])*v[377];
		v[3262] = v[3261] * v[349] + (v[2504] * v[3224] - v[3259] * v[349])*v[377];
		v[3263] = v[3261] * v[359] + (v[2506] * v[3224] - v[3259] * v[359])*v[377];
		v[3265] = v[313] * v[4033];
		v[3266] = v[323] * v[4033];
		v[3267] = v[333] * v[4033];
		v[3268] = v[287] * v[4033];
		v[3269] = v[297] * v[4033];
		v[3270] = v[307] * v[4033];
		v[3271] = -(v[313] * v[4034]);
		v[3272] = -(v[323] * v[4034]);
		v[3273] = -(v[333] * v[4034]);
		v[3274] = -(v[287] * v[4034]);
		v[3275] = -(v[297] * v[4034]);
		v[3276] = -(v[307] * v[4034]);
		v[3278] = v[313] * v[4035];
		v[3279] = v[323] * v[4035];
		v[3280] = v[333] * v[4035];
		v[3281] = v[287] * v[4035];
		v[3282] = v[297] * v[4035];
		v[3283] = v[307] * v[4035];
		v[3288] = (*cn)*(v[3224] * v[3284] + v[3226] * v[3285] + v[3236] * v[3286]) + 2e0*v[3277] * v[3287]
			+ v[3238] * v[3852];
		v[3293] = 2e0*v[3259] * v[3290] + (*cn)*(v[3224] * v[3289] + v[3226] * v[3291] + v[3236] * v[3292])
			+ v[375] * v[4036];
		v[3298] = 2e0*v[3264] * v[3294] + (*cn)*(v[3224] * v[3295] + v[3226] * v[3296] + v[3236] * v[3297])
			+ v[376] * v[4036];
		v[4037] = ((v[3298] * v[360] + v[3293] * v[361] + v[3288] * v[362])*v[373]) / v[1644];
		v[3300] = v[3288] * v[374] + v[362] * v[4037];
		v[3301] = v[3293] * v[374] + v[361] * v[4037];
		v[3302] = v[3298] * v[374] + v[360] * v[4037];
		v[3303] = -(v[1003] * v[3160]) - v[263] * v[3300];
		v[3304] = -(v[261] * v[3300]) + v[1003] * v[4038];
		v[3305] = -(v[1003] * v[3159]) - v[260] * v[3300];
		v[3306] = v[210] * v[3300] + v[3307];
		v[3308] = v[211] * v[3300] - v[3307];
		v[3309] = -(v[1002] * v[3160]) - v[263] * v[3301];
		v[3310] = -(v[261] * v[3301]) + v[1002] * v[4038];
		v[3311] = -(v[1002] * v[3159]) - v[260] * v[3301];
		v[3312] = v[210] * v[3301] + v[3313];
		v[3314] = v[211] * v[3301] - v[3313];
		v[3315] = -(v[1001] * v[3160]) - v[263] * v[3302];
		v[3316] = -(v[261] * v[3302]) + v[1001] * v[4038];
		v[3317] = -(v[1001] * v[3159]) - v[260] * v[3302];
		v[3318] = v[210] * v[3302] + v[3319];
		v[3320] = v[211] * v[3302] - v[3319];
		v[3321] = v[2351] * v[3241];
		v[3322] = v[2350] * v[3241];
		v[3323] = v[2349] * v[3241];
		v[3324] = v[2346] * v[3242];
		v[3325] = v[3241] * v[346] + v[3242] * v[348];
		v[3326] = v[2344] * v[3242];
		v[3327] = v[3321] + v[3324];
		v[3328] = v[2345] * v[3242] + v[3325] * v[4039];
		v[3329] = v[2339] * v[3243];
		v[3330] = v[3241] * v[340] + v[3243] * v[348];
		v[3331] = v[3242] * v[340] + v[3243] * v[346];
		v[3332] = v[2341] * v[3243] + v[3330] * v[4040];
		v[3333] = v[3328] + v[3332];
		v[3334] = v[2340] * v[3243] + v[3331] * v[4041];
		v[3335] = v[3322] + v[3334];
		v[3336] = v[2336] * v[3244];
		v[3337] = v[2335] * v[3244];
		v[3338] = v[2334] * v[3244];
		v[3339] = v[2331] * v[3245];
		v[3340] = v[320] * v[3244] + v[322] * v[3245];
		v[3341] = v[2329] * v[3245];
		v[3342] = v[3336] + v[3339];
		v[3343] = v[2330] * v[3245] + v[3340] * v[4042];
		v[3344] = v[2324] * v[3246];
		v[3345] = v[314] * v[3244] + v[322] * v[3246];
		v[3346] = v[314] * v[3245] + v[320] * v[3246];
		v[3347] = v[2326] * v[3246] + v[3345] * v[4043];
		v[3348] = v[3343] + v[3347];
		v[3349] = v[2325] * v[3246] + v[3346] * v[4044];
		v[3350] = v[3337] + v[3349];
		v[3351] = v[2321] * v[3247];
		v[3352] = v[2320] * v[3247];
		v[3353] = v[2319] * v[3247];
		v[3354] = v[2316] * v[3248];
		v[3355] = v[294] * v[3247] + v[296] * v[3248];
		v[3356] = v[2314] * v[3248];
		v[3357] = v[3351] + v[3354];
		v[3358] = v[2315] * v[3248] + v[3355] * v[4045];
		v[3359] = v[2309] * v[3249];
		v[3360] = v[288] * v[3247] + v[296] * v[3249];
		v[3361] = v[288] * v[3248] + v[294] * v[3249];
		v[3362] = v[2311] * v[3249] + v[3360] * v[4046];
		v[3363] = v[3358] + v[3362];
		v[3364] = v[2310] * v[3249] + v[3361] * v[4047];
		v[3365] = v[3352] + v[3364];
		v[3366] = -(v[3250] * v[950]);
		v[3367] = -(v[3250] * v[949]);
		v[3368] = -(v[3250] * v[951]);
		v[3369] = -(v[3252] * v[951]);
		v[3370] = -(v[3252] * v[950]);
		v[3371] = -(v[3252] * v[949]);
		v[3372] = -(v[3253] * v[951]);
		v[3373] = -(v[3253] * v[949]);
		v[3374] = -(v[3253] * v[950]);
		v[3375] = -(v[3254] * v[951]);
		v[3376] = -(v[3254] * v[950]);
		v[3377] = -(v[3254] * v[949]);
		v[3378] = -(v[3257] * v[951]);
		v[3379] = -(v[3257] * v[949]);
		v[3380] = -(v[3257] * v[950]);
		v[3381] = -(v[3258] * v[949]);
		v[3382] = -(v[3258] * v[951]);
		v[3383] = -(v[3258] * v[950]);
		v[3384] = -(v[3260] * v[951]);
		v[3385] = -(v[3260] * v[950]);
		v[3386] = -(v[3260] * v[949]);
		v[3387] = -(v[3262] * v[950]);
		v[3388] = -(v[3262] * v[951]);
		v[3389] = -(v[3262] * v[949]);
		v[3390] = -(QABi[0][2] * v[3233]) - QABi[1][2] * v[3234] - QABi[2][2] * v[3237] - v[252] * v[3300] - v[249] * v[3301]
			- v[246] * v[3302] + v[3250] * v[710] + v[3252] * v[711] + v[3253] * v[712] + v[3254] * v[725] + v[3257] * v[726]
			+ v[3258] * v[727] + v[3260] * v[740] + v[3262] * v[741] + v[3263] * v[742];
		v[3391] = -(QABi[0][0] * v[3233]) - QABi[1][0] * v[3234] - QABi[2][0] * v[3237] - v[250] * v[3300] - v[247] * v[3301]
			- v[244] * v[3302] + v[3250] * v[704] + v[3252] * v[705] + v[3253] * v[706] + v[3254] * v[719] + v[3257] * v[720]
			+ v[3258] * v[721] + v[3260] * v[734] + v[3262] * v[735] + v[3263] * v[736];
		v[3392] = -(v[3263] * v[950]);
		v[3393] = -(v[3263] * v[951]);
		v[3394] = -(v[3263] * v[949]);
		v[3395] = QABi[2][2] * v[3303] + QABi[2][1] * v[3304] + QABi[2][0] * v[3305] + v[348] * (v[3399] + v[4048]);
		v[3398] = QABi[1][2] * v[3303] + QABi[1][1] * v[3304] + QABi[1][0] * v[3305] + v[3331] * v[3647] + v[346] * v[4048];
		v[3400] = QABi[0][2] * v[3303] + QABi[0][1] * v[3304] + QABi[0][0] * v[3305] + (v[3396] + v[3399])*v[340];
		v[3401] = QABi[2][2] * v[3309] + QABi[2][1] * v[3310] + QABi[2][0] * v[3311] + v[3330] * v[3645] + v[348] * v[4049];
		v[3404] = QABi[1][2] * v[3309] + QABi[1][1] * v[3310] + QABi[1][0] * v[3311] + v[346] * (v[3405] + v[4049]);
		v[3406] = QABi[0][2] * v[3309] + QABi[0][1] * v[3310] + QABi[0][0] * v[3311] + v[340] * (v[3403] + v[3405]);
		v[3407] = QABi[2][2] * v[3315] + QABi[2][1] * v[3316] + QABi[2][0] * v[3317] + v[1871] * v[3325] + (v[3408] + v[3410]
			)*v[348];
		v[3409] = QABi[1][2] * v[3315] + QABi[1][1] * v[3316] + QABi[1][0] * v[3317] + v[346] * v[4050];
		v[3412] = QABi[0][2] * v[3315] + QABi[0][1] * v[3316] + QABi[0][0] * v[3317] + v[340] * (v[3410] + v[4050]);
		v[3413] = v[3366] + v[3375] + v[3381] + v[3392] + 2e0*v[3326] * v[666] - v[3333] * v[669] - v[3327] * v[672];
		v[3414] = -v[3370] - v[3373] - v[3378] - v[3393] - v[3333] * v[666] + 2e0*v[3329] * v[669] + v[3335] * v[672];
		v[3415] = v[3072] * v[3213] + v[225] * v[3409] - v[3366] * v[659] + v[3370] * v[660] - v[3374] * v[661];
		v[3416] = -v[3367] - v[3379] - v[3384] - v[3387] - v[3327] * v[666] + v[3335] * v[669] + 2e0*v[3323] * v[672];
		v[3417] = v[3071] * v[3213] + v[225] * v[3407] - v[3367] * v[659] + v[3371] * v[660] - v[3373] * v[661];
		v[3418] = v[3062] * v[3213] + v[225] * v[3406] - v[3375] * v[659] + v[3378] * v[660] - v[3382] * v[661];
		v[3419] = v[3372] + v[3383];
		v[3420] = v[3060] * v[3213] + v[225] * v[3401] - v[3377] * v[659] + v[3379] * v[660] - v[3381] * v[661];
		v[3421] = v[3053] * v[3213] + v[225] * v[3400] - v[3384] * v[659] + v[3388] * v[660] - v[3393] * v[661];
		v[3422] = v[3052] * v[3213] + v[225] * v[3398] - v[3385] * v[659] + v[3387] * v[660] - v[3392] * v[661];
		v[3423] = v[3376] + v[3386];
		v[3424] = v[3369] + v[3389];
		v[3426] = v[215] * v[3318] + v[1439] * v[3427];
		v[3428] = v[214] * v[3318] - v[3427] * v[993];
		v[3429] = v[215] * v[3320] + v[1439] * v[3430];
		v[3431] = v[214] * v[3320] - v[3430] * v[993];
		v[3432] = v[215] * v[3312] + v[1439] * v[3433];
		v[3434] = v[214] * v[3312] - v[3433] * v[993];
		v[3435] = v[215] * v[3314] + v[1439] * v[3436];
		v[3437] = v[214] * v[3314] - v[3436] * v[993];
		v[3438] = v[215] * v[3306] + v[1439] * v[3439];
		v[3440] = v[214] * v[3306] - v[3439] * v[993];
		v[3441] = v[215] * v[3308] + v[1439] * v[3442];
		v[3443] = v[214] * v[3308] - v[3442] * v[993];
		v[3444] = v[3265] * v[954];
		v[3445] = v[3265] * v[953];
		v[3446] = v[3265] * v[952];
		v[3447] = v[3266] * v[953];
		v[3448] = v[3266] * v[954];
		v[3449] = v[3266] * v[952];
		v[3450] = v[3267] * v[953];
		v[3451] = v[3267] * v[954];
		v[3452] = v[3267] * v[952];
		v[3453] = v[3268] * v[957];
		v[3454] = v[3268] * v[956];
		v[3455] = v[3268] * v[955];
		v[3456] = v[3269] * v[956];
		v[3457] = v[3269] * v[957];
		v[3458] = v[3269] * v[955];
		v[3459] = v[3270] * v[956];
		v[3460] = v[3270] * v[957];
		v[3461] = v[3270] * v[955];
		v[3462] = v[3271] * v[954];
		v[3463] = v[3271] * v[953];
		v[3464] = v[3271] * v[952];
		v[3465] = v[3272] * v[954];
		v[3466] = v[3272] * v[952];
		v[3467] = v[3272] * v[953];
		v[3468] = v[3273] * v[952];
		v[3469] = v[3273] * v[954];
		v[3470] = v[3273] * v[953];
		v[3471] = v[3274] * v[957];
		v[3472] = v[3274] * v[956];
		v[3473] = v[3274] * v[955];
		v[3474] = v[3275] * v[957];
		v[3475] = v[3275] * v[955];
		v[3476] = v[3275] * v[956];
		v[3477] = v[3276] * v[955];
		v[3478] = v[3276] * v[957];
		v[3479] = v[3276] * v[956];
		v[3480] = v[3278] * v[953];
		v[3481] = v[3278] * v[952];
		v[3482] = v[3278] * v[954];
		v[3483] = v[3279] * v[954];
		v[3484] = v[3279] * v[953];
		v[3485] = v[3279] * v[952];
		v[3486] = v[3280] * v[954];
		v[3487] = v[3280] * v[952];
		v[3488] = v[3280] * v[953];
		v[3489] = v[3281] * v[956];
		v[3490] = v[3281] * v[955];
		v[3491] = v[3281] * v[957];
		v[3492] = v[3282] * v[957];
		v[3493] = v[3282] * v[956];
		v[3494] = v[3282] * v[955];
		v[3495] = v[3283] * v[957];
		v[3496] = v[3283] * v[955];
		v[3497] = v[3283] * v[956];
		v[3500] = (-(v[1337] * v[3227]) - v[1338] * v[3228] - v[1339] * v[3229] + v[1340] * v[3230] + v[1341] * v[3231]
			+ v[1342] * v[3232] + v[273] * v[3251] - v[276] * v[3251] + v[274] * v[3255] - v[277] * v[3255] - v[3239] * v[3657]
			+ v[3240] * v[3657] - v[3302] * v[387] + v[3302] * v[388] - v[3301] * v[390] + v[3301] * v[391] - v[3300] * v[393]
			+ v[3300] * v[394] - v[3498] * v[4029] + v[3499] * v[4029] - v[275] * v[4051] + v[278] * v[4051] - v[3223] * v[4270]
			+ v[3223] * v[4272] + v[3265] * v[605] + v[3266] * v[606] + v[3267] * v[607] - v[3268] * v[608] - v[3269] * v[609]
			- v[3270] * v[610] + v[3271] * v[617] + v[3272] * v[618] + v[3273] * v[619] - v[3274] * v[620] - v[3275] * v[621]
			- v[3276] * v[622] + v[3278] * v[629] + v[3279] * v[630] + v[3280] * v[631] - v[3281] * v[632] - v[3282] * v[633]
			- v[3283] * v[634] + v[7954 + i3049] - v[7972 + i3049]) / 2e0;
		v[3501] = QBAi[2][1] * v[3441] + QBAi[2][0] * v[3443] + v[322] * (v[3505] + v[4052]);
		v[3504] = QBAi[1][1] * v[3441] + QBAi[1][0] * v[3443] + v[3346] * v[3642] + v[320] * v[4052];
		v[3506] = QBAi[0][1] * v[3441] + QBAi[0][0] * v[3443] + v[314] * (v[3502] + v[3505]);
		v[3507] = QBAi[2][1] * v[3435] + QBAi[2][0] * v[3437] + v[3345] * v[3640] + v[322] * v[4053];
		v[3510] = QBAi[1][1] * v[3435] + QBAi[1][0] * v[3437] + v[320] * (v[3511] + v[4053]);
		v[3512] = QBAi[0][1] * v[3435] + QBAi[0][0] * v[3437] + v[314] * (v[3509] + v[3511]);
		v[3513] = v[1865] * v[3340] + QBAi[2][1] * v[3429] + QBAi[2][0] * v[3431] + v[322] * (v[3514] + v[3516]);
		v[3515] = QBAi[1][1] * v[3429] + QBAi[1][0] * v[3431] + v[320] * v[4054];
		v[3518] = QBAi[0][1] * v[3429] + QBAi[0][0] * v[3431] + v[314] * (v[3516] + v[4054]);
		v[3519] = QAAi[2][1] * v[3438] + QAAi[2][0] * v[3440] + v[296] * (v[3523] + v[4055]);
		v[3522] = QAAi[1][1] * v[3438] + QAAi[1][0] * v[3440] + v[3361] * v[3637] + v[294] * v[4055];
		v[3524] = QAAi[0][1] * v[3438] + QAAi[0][0] * v[3440] + v[288] * (v[3520] + v[3523]);
		v[3525] = QAAi[2][1] * v[3432] + QAAi[2][0] * v[3434] + v[3360] * v[3635] + v[296] * v[4056];
		v[3528] = QAAi[1][1] * v[3432] + QAAi[1][0] * v[3434] + v[294] * (v[3529] + v[4056]);
		v[3530] = QAAi[0][1] * v[3432] + QAAi[0][0] * v[3434] + v[288] * (v[3527] + v[3529]);
		v[3531] = v[1859] * v[3355] + QAAi[2][1] * v[3426] + QAAi[2][0] * v[3428] + v[296] * (v[3532] + v[3534]);
		v[3533] = QAAi[1][1] * v[3426] + QAAi[1][0] * v[3428] + v[294] * v[4057];
		v[3536] = QAAi[0][1] * v[3426] + QAAi[0][0] * v[3428] + v[288] * (v[3534] + v[4057]);
		v[3537] = v[3450] + v[3462] + v[3468] + v[3480] + 2e0*v[3341] * v[495] - v[3348] * v[498] - v[3342] * v[501];
		v[3538] = -v[3451] - v[3465] - v[3484] - v[3487] - v[3348] * v[495] + 2e0*v[3344] * v[498] + v[3350] * v[501];
		v[3539] = v[3075] * v[3202] + v[165] * v[3515] - v[3480] * v[488] + v[3484] * v[489] - v[3488] * v[490];
		v[3540] = -v[3444] - v[3447] - v[3466] - v[3481] - v[3342] * v[495] + v[3350] * v[498] + 2e0*v[3338] * v[501];
		v[3541] = v[3074] * v[3202] + v[165] * v[3513] - v[3481] * v[488] + v[3485] * v[489] - v[3487] * v[490];
		v[3542] = v[3065] * v[3202] + v[165] * v[3512] - v[3462] * v[488] + v[3465] * v[489] - v[3469] * v[490];
		v[3543] = v[3470] + v[3486];
		v[3544] = v[3063] * v[3202] + v[165] * v[3507] - v[3464] * v[488] + v[3466] * v[489] - v[3468] * v[490];
		v[3545] = v[3056] * v[3202] + v[165] * v[3506] - v[3444] * v[488] + v[3448] * v[489] - v[3451] * v[490];
		v[3546] = v[3055] * v[3202] + v[165] * v[3504] - v[3445] * v[488] + v[3447] * v[489] - v[3450] * v[490];
		v[3547] = v[3446] + v[3463];
		v[3548] = v[3449] + v[3483];
		v[3550] = v[3459] + v[3471] + v[3477] + v[3489] + 2e0*v[3356] * v[450] - v[3363] * v[453] - v[3357] * v[456];
		v[3551] = -v[3460] - v[3474] - v[3493] - v[3496] - v[3363] * v[450] + 2e0*v[3359] * v[453] + v[3365] * v[456];
		v[3552] = v[3078] * v[3191] + v[146] * v[3533] - v[3489] * v[443] + v[3493] * v[444] - v[3497] * v[445];
		v[3553] = -v[3453] - v[3456] - v[3475] - v[3490] - v[3357] * v[450] + v[3365] * v[453] + 2e0*v[3353] * v[456];
		v[3554] = v[3077] * v[3191] + v[146] * v[3531] - v[3490] * v[443] + v[3494] * v[444] - v[3496] * v[445];
		v[3555] = v[3068] * v[3191] + v[146] * v[3530] - v[3471] * v[443] + v[3474] * v[444] - v[3478] * v[445];
		v[3556] = v[3479] + v[3495];
		v[3557] = v[3066] * v[3191] + v[146] * v[3525] - v[3473] * v[443] + v[3475] * v[444] - v[3477] * v[445];
		v[3558] = v[3059] * v[3191] + v[146] * v[3524] - v[3453] * v[443] + v[3457] * v[444] - v[3460] * v[445];
		v[3559] = v[3058] * v[3191] + v[146] * v[3522] - v[3454] * v[443] + v[3456] * v[444] - v[3459] * v[445];
		v[3560] = v[3455] + v[3472];
		v[3561] = v[3458] + v[3492];
		v[3565] = v[1440] * (-(QABi[0][1] * v[3233]) - QABi[1][1] * v[3234] - QABi[2][1] * v[3237] - v[251] * v[3300]
			- v[248] * v[3301] - v[245] * v[3302] - v[3158] * v[3563] + v[3250] * v[707] + v[3252] * v[708] + v[3253] * v[709]
			+ v[3254] * v[722] + v[3257] * v[723] + v[3258] * v[724] + v[3260] * v[737] + v[3262] * v[738] + v[3263] * v[739])
			+ v[994] * (-(v[3158] * v[3564]) + (-(v[1438] * v[3391]) + v[4058] * v[992])*(*xfac) + (v[1438] * v[4060]
				+ v[3390] * v[992])*(*zfac));
		v[3566] = (-2e0*v[3321] + 2e0*v[3324] - v[3368] * v[662] - 2e0*v[3366] * v[667] - 2e0*v[3367] * v[673]
			- 2e0*v[3375] * v[677] - v[3376] * v[681] - 2e0*v[3377] * v[686] - 2e0*v[3384] * v[690] - 2e0*v[3385] * v[694]
			- v[3386] * v[699]) / 2e0;
		v[3567] = (v[3073] * v[3213] + v[225] * v[3412] - v[3368] * v[659] + v[3369] * v[660] - v[3372] * v[661]) / 2e0;
		v[3569] = -v[3322] + v[3334] + (v[3369] * v[662]) / 2e0 + v[3370] * v[667] + v[3371] * v[673] + v[3378] * v[677] +
			(v[3380] * v[681]) / 2e0 + v[3379] * v[686] + v[3388] * v[690] + v[3387] * v[694] + (v[3389] * v[699]) / 2e0;
		v[3570] = (v[3061] * v[3213] + v[225] * v[3404] - v[3376] * v[659] + v[3380] * v[660] - v[3383] * v[661]) / 2e0;
		v[3571] = (-2e0*v[3328] + 2e0*v[3332] - v[3372] * v[662] - 2e0*v[3374] * v[667] - 2e0*v[3373] * v[673]
			- 2e0*v[3382] * v[677] - v[3383] * v[681] - 2e0*v[3381] * v[686] - 2e0*v[3393] * v[690] - 2e0*v[3392] * v[694]
			- v[3394] * v[699]) / 2e0;
		v[3617] = v[2057] * v[3569] - v[3570] + v[2059] * v[3571] + 8e0*v[1854] * (v[3089] * v[3153] - v[3566] * v[437])
			- 2e0*v[1252] * (v[3073] * v[3152] + v[3061] * v[3154] + v[3051] * v[3156] + 2e0*v[3072] * v[3179]
				+ 2e0*v[3062] * v[3181] + 2e0*v[3060] * v[3182] + 2e0*v[3052] * v[3184] + 2e0*v[3071] * v[3185]
				+ 2e0*v[3053] * v[3187] + 4e0*v[225] * (v[2348] * v[3241] + v[2343] * v[3242] + v[2338] * v[3243] + v[3323]
					+ v[2167] * v[3325] + v[3326] + v[3329] + v[2165] * v[3330] + v[2161] * v[3331]) + 2e0*v[3371] - 2e0*v[3374]
				- 2e0*v[3377] + 2e0*v[3382] + 2e0*v[3385] - 2e0*v[3388] + dB[4] * v[3413] - dB[3] * v[3414] - dB[5] * v[3416]
				- v[3423] * v[437] - v[3424] * v[440] - v[3419] * v[442] + v[3412] * v[662] + 2e0*v[3409] * v[667] + 2e0*v[3407] * v[673]
				+ 2e0*v[3406] * v[677] + v[3404] * v[681] + 2e0*v[3401] * v[686] + 2e0*v[3400] * v[690] + 2e0*v[3398] * v[694]
				+ v[3395] * v[699]);
		v[4068] = v[3617] + (-(v[3051] * v[3213]) - v[225] * v[3395] + v[3386] * v[659] - v[3389] * v[660] + v[3394] * v[661])
			/ 2e0;
		v[4067] = (v[3417] + v[3421]) / 2e0;
		v[3574] = v[3420] + v[3422];
		v[3575] = v[3415] + v[3418];
		v[3576] = v[992] * (v[3918] * v[4060] + (-(v[262] * v[3391]) + v[3070] * v[4059])*(*xfac)) - v[1438] *
			(v[3725] * v[4058] + (v[262] * v[3390] - v[3069] * v[4059])*(*zfac));
		v[3579] = v[1439] * (v[191] * v[3306] + v[200] * v[3308] + v[188] * v[3312] + v[197] * v[3314] + v[185] * v[3318]
			+ v[194] * v[3320] - v[3223] * v[3578] + v[210] * (QAAi[0][1] * v[3227] + QAAi[1][1] * v[3228] + QAAi[2][1] * v[3229]
				+ v[3281] * v[536] + v[3282] * v[537] + v[3283] * v[538] + v[3274] * v[545] + v[3275] * v[546] + v[3276] * v[547]
				+ v[3268] * v[554] + v[3269] * v[555] + v[3270] * v[556]) + v[211] * (QBAi[0][1] * v[3230] + QBAi[1][1] * v[3231]
					+ QBAi[2][1] * v[3232] + v[3278] * v[563] + v[3279] * v[564] + v[3280] * v[565] + v[3271] * v[572] + v[3272] * v[573]
					+ v[3273] * v[574] + v[3265] * v[581] + v[3266] * v[582] + v[3267] * v[583])) - (v[190] * v[3306] + v[199] * v[3308]
						+ v[187] * v[3312] + v[196] * v[3314] + v[184] * v[3318] + v[193] * v[3320] + v[3223] * v[3577] + v[210] *
						(QAAi[0][0] * v[3227] + QAAi[1][0] * v[3228] + QAAi[2][0] * v[3229] + v[3281] * v[533] + v[3282] * v[534]
							+ v[3283] * v[535] + v[3274] * v[542] + v[3275] * v[543] + v[3276] * v[544] + v[3268] * v[551] + v[3269] * v[552]
							+ v[3270] * v[553]) + v[211] * (QBAi[0][0] * v[3230] + QBAi[1][0] * v[3231] + QBAi[2][0] * v[3232] + v[3278] * v[560]
								+ v[3279] * v[561] + v[3280] * v[562] + v[3271] * v[569] + v[3272] * v[570] + v[3273] * v[571] + v[3265] * v[578]
								+ v[3266] * v[579] + v[3267] * v[580]))*v[993];
		v[3580] = (-2e0*v[3336] + 2e0*v[3339] - v[3482] * v[491] - 2e0*v[3480] * v[496] - 2e0*v[3481] * v[502]
			- 2e0*v[3462] * v[506] - v[3463] * v[510] - 2e0*v[3464] * v[515] - 2e0*v[3444] * v[519] - 2e0*v[3445] * v[523]
			- v[3446] * v[528]) / 2e0;
		v[3581] = (v[3076] * v[3202] + v[165] * v[3518] - v[3482] * v[488] + v[3483] * v[489] - v[3486] * v[490]) / 2e0;
		v[3583] = -v[3337] + v[3349] + (v[3483] * v[491]) / 2e0 + v[3484] * v[496] + v[3485] * v[502] + v[3465] * v[506] +
			(v[3467] * v[510]) / 2e0 + v[3466] * v[515] + v[3448] * v[519] + v[3447] * v[523] + (v[3449] * v[528]) / 2e0;
		v[3584] = (v[3064] * v[3202] + v[165] * v[3510] - v[3463] * v[488] + v[3467] * v[489] - v[3470] * v[490]) / 2e0;
		v[3585] = (-2e0*v[3343] + 2e0*v[3347] - v[3486] * v[491] - 2e0*v[3488] * v[496] - 2e0*v[3487] * v[502]
			- 2e0*v[3469] * v[506] - v[3470] * v[510] - 2e0*v[3468] * v[515] - 2e0*v[3451] * v[519] - 2e0*v[3450] * v[523]
			- v[3452] * v[528]) / 2e0;
		v[3610] = v[1976] * v[3583] - v[3584] + v[1978] * v[3585] + 8e0*v[1852] * (v[3101] * v[3147] - v[3580] * v[431])
			- 2e0*v[1238] * (v[3076] * v[3146] + v[3064] * v[3148] + v[3054] * v[3150] + 2e0*v[3075] * v[3170]
				+ 2e0*v[3065] * v[3172] + 2e0*v[3063] * v[3173] + 2e0*v[3055] * v[3175] + 2e0*v[3074] * v[3176]
				+ 2e0*v[3056] * v[3178] + 4e0*v[165] * (v[2333] * v[3244] + v[2328] * v[3245] + v[2323] * v[3246] + v[3338]
					+ v[2139] * v[3340] + v[3341] + v[3344] + v[2137] * v[3345] + v[2133] * v[3346]) + 2e0*v[3445] - 2e0*v[3448]
				- 2e0*v[3464] + 2e0*v[3469] + 2e0*v[3485] - 2e0*v[3488] + dA[10] * v[3537] - dA[9] * v[3538] - dA[11] * v[3540]
				- v[3547] * v[431] - v[3548] * v[434] - v[3543] * v[436] + v[3518] * v[491] + 2e0*v[3515] * v[496] + 2e0*v[3513] * v[502]
				+ 2e0*v[3512] * v[506] + v[3510] * v[510] + 2e0*v[3507] * v[515] + 2e0*v[3506] * v[519] + 2e0*v[3504] * v[523]
				+ v[3501] * v[528]);
		v[4066] = v[3610] + (-(v[3054] * v[3202]) - v[165] * v[3501] + v[3446] * v[488] - v[3449] * v[489] + v[3452] * v[490])
			/ 2e0;
		v[4065] = (v[3541] + v[3545]) / 2e0;
		v[3588] = v[3544] + v[3546];
		v[3589] = v[3539] + v[3542];
		v[3590] = (-2e0*v[3351] + 2e0*v[3354] - v[3491] * v[446] - 2e0*v[3489] * v[451] - 2e0*v[3490] * v[457]
			- 2e0*v[3471] * v[461] - v[3472] * v[465] - 2e0*v[3473] * v[470] - 2e0*v[3453] * v[474] - 2e0*v[3454] * v[478]
			- v[3455] * v[483]) / 2e0;
		v[3591] = (v[3079] * v[3191] + v[146] * v[3536] - v[3491] * v[443] + v[3492] * v[444] - v[3495] * v[445]) / 2e0;
		v[3593] = -v[3352] + v[3364] + (v[3492] * v[446]) / 2e0 + v[3493] * v[451] + v[3494] * v[457] + v[3474] * v[461] +
			(v[3476] * v[465]) / 2e0 + v[3475] * v[470] + v[3457] * v[474] + v[3456] * v[478] + (v[3458] * v[483]) / 2e0;
		v[3594] = (v[3067] * v[3191] + v[146] * v[3528] - v[3472] * v[443] + v[3476] * v[444] - v[3479] * v[445]) / 2e0;
		v[3595] = (-2e0*v[3358] + 2e0*v[3362] - v[3495] * v[446] - 2e0*v[3497] * v[451] - 2e0*v[3496] * v[457]
			- 2e0*v[3478] * v[461] - v[3479] * v[465] - 2e0*v[3477] * v[470] - 2e0*v[3460] * v[474] - 2e0*v[3459] * v[478]
			- v[3461] * v[483]) / 2e0;
		v[3603] = v[1949] * v[3593] - v[3594] + v[1951] * v[3595] + 8e0*v[1850] * (v[3111] * v[3141] - v[3590] * v[425])
			- 2e0*v[1224] * (v[3079] * v[3140] + v[3067] * v[3142] + v[3057] * v[3144] + 2e0*v[3078] * v[3161]
				+ 2e0*v[3068] * v[3163] + 2e0*v[3066] * v[3164] + 2e0*v[3058] * v[3166] + 2e0*v[3077] * v[3167]
				+ 2e0*v[3059] * v[3169] + 4e0*v[146] * (v[2318] * v[3247] + v[2313] * v[3248] + v[2308] * v[3249] + v[3353]
					+ v[2111] * v[3355] + v[3356] + v[3359] + v[2109] * v[3360] + v[2105] * v[3361]) + 2e0*v[3454] - 2e0*v[3457]
				- 2e0*v[3473] + 2e0*v[3478] + 2e0*v[3494] - 2e0*v[3497] + dA[4] * v[3550] - dA[3] * v[3551] - dA[5] * v[3553]
				- v[3560] * v[425] - v[3561] * v[428] - v[3556] * v[430] + v[3536] * v[446] + 2e0*v[3533] * v[451] + 2e0*v[3531] * v[457]
				+ 2e0*v[3530] * v[461] + v[3528] * v[465] + 2e0*v[3525] * v[470] + 2e0*v[3524] * v[474] + 2e0*v[3522] * v[478]
				+ v[3519] * v[483]);
		v[4064] = v[3603] + (-(v[3057] * v[3191]) - v[146] * v[3519] + v[3455] * v[443] - v[3458] * v[444] + v[3461] * v[445])
			/ 2e0;
		v[4063] = (v[3554] + v[3558]) / 2e0;
		v[3598] = v[3557] + v[3559];
		v[3599] = v[3552] + v[3555];
		v[8045] = 0e0;
		v[8046] = 0e0;
		v[8047] = 0e0;
		v[8048] = 2e0*v[3600];
		v[8049] = v[3601];
		v[8050] = v[3602];
		v[8051] = 0e0;
		v[8052] = 0e0;
		v[8053] = 0e0;
		v[8054] = 0e0;
		v[8055] = 0e0;
		v[8056] = 0e0;
		v[8057] = 0e0;
		v[8058] = 0e0;
		v[8059] = 0e0;
		v[8060] = 0e0;
		v[8061] = 0e0;
		v[8062] = 0e0;
		v[8027] = 0e0;
		v[8028] = 0e0;
		v[8029] = 0e0;
		v[8030] = v[3601];
		v[8031] = 2e0*v[3604];
		v[8032] = v[3605];
		v[8033] = 0e0;
		v[8034] = 0e0;
		v[8035] = 0e0;
		v[8036] = 0e0;
		v[8037] = 0e0;
		v[8038] = 0e0;
		v[8039] = 0e0;
		v[8040] = 0e0;
		v[8041] = 0e0;
		v[8042] = 0e0;
		v[8043] = 0e0;
		v[8044] = 0e0;
		v[8009] = 0e0;
		v[8010] = 0e0;
		v[8011] = 0e0;
		v[8012] = v[3602];
		v[8013] = v[3605];
		v[8014] = 2e0*v[3606];
		v[8015] = 0e0;
		v[8016] = 0e0;
		v[8017] = 0e0;
		v[8018] = 0e0;
		v[8019] = 0e0;
		v[8020] = 0e0;
		v[8021] = 0e0;
		v[8022] = 0e0;
		v[8023] = 0e0;
		v[8024] = 0e0;
		v[8025] = 0e0;
		v[8026] = 0e0;
		v[7991] = 0e0;
		v[7992] = 0e0;
		v[7993] = 0e0;
		v[7994] = 0e0;
		v[7995] = 0e0;
		v[7996] = 0e0;
		v[7997] = 0e0;
		v[7998] = 0e0;
		v[7999] = 0e0;
		v[8000] = 2e0*v[3607];
		v[8001] = v[3608];
		v[8002] = v[3609];
		v[8003] = 0e0;
		v[8004] = 0e0;
		v[8005] = 0e0;
		v[8006] = 0e0;
		v[8007] = 0e0;
		v[8008] = 0e0;
		v[8063] = 0e0;
		v[8064] = 0e0;
		v[8065] = 0e0;
		v[8066] = 0e0;
		v[8067] = 0e0;
		v[8068] = 0e0;
		v[8069] = 0e0;
		v[8070] = 0e0;
		v[8071] = 0e0;
		v[8072] = v[3608];
		v[8073] = 2e0*v[3611];
		v[8074] = v[3612];
		v[8075] = 0e0;
		v[8076] = 0e0;
		v[8077] = 0e0;
		v[8078] = 0e0;
		v[8079] = 0e0;
		v[8080] = 0e0;
		v[8081] = 0e0;
		v[8082] = 0e0;
		v[8083] = 0e0;
		v[8084] = 0e0;
		v[8085] = 0e0;
		v[8086] = 0e0;
		v[8087] = 0e0;
		v[8088] = 0e0;
		v[8089] = 0e0;
		v[8090] = v[3609];
		v[8091] = v[3612];
		v[8092] = 2e0*v[3613];
		v[8093] = 0e0;
		v[8094] = 0e0;
		v[8095] = 0e0;
		v[8096] = 0e0;
		v[8097] = 0e0;
		v[8098] = 0e0;
		v[8099] = 0e0;
		v[8100] = 0e0;
		v[8101] = 0e0;
		v[8102] = 0e0;
		v[8103] = 0e0;
		v[8104] = 0e0;
		v[8105] = 0e0;
		v[8106] = 0e0;
		v[8107] = 0e0;
		v[8108] = 0e0;
		v[8109] = 0e0;
		v[8110] = 0e0;
		v[8111] = 0e0;
		v[8112] = 0e0;
		v[8113] = 0e0;
		v[8114] = 2e0*v[3614];
		v[8115] = v[3615];
		v[8116] = v[3616];
		v[8117] = 0e0;
		v[8118] = 0e0;
		v[8119] = 0e0;
		v[8120] = 0e0;
		v[8121] = 0e0;
		v[8122] = 0e0;
		v[8123] = 0e0;
		v[8124] = 0e0;
		v[8125] = 0e0;
		v[8126] = 0e0;
		v[8127] = 0e0;
		v[8128] = 0e0;
		v[8129] = 0e0;
		v[8130] = 0e0;
		v[8131] = 0e0;
		v[8132] = v[3615];
		v[8133] = 2e0*v[3618];
		v[8134] = v[3619];
		v[8135] = 0e0;
		v[8136] = 0e0;
		v[8137] = 0e0;
		v[8138] = 0e0;
		v[8139] = 0e0;
		v[8140] = 0e0;
		v[8141] = 0e0;
		v[8142] = 0e0;
		v[8143] = 0e0;
		v[8144] = 0e0;
		v[8145] = 0e0;
		v[8146] = 0e0;
		v[8147] = 0e0;
		v[8148] = 0e0;
		v[8149] = 0e0;
		v[8150] = v[3616];
		v[8151] = v[3619];
		v[8152] = 2e0*v[3620];
		v[8153] = v[3318];
		v[8154] = v[3312];
		v[8155] = v[3306];
		v[8156] = -v[3557] + v[3559] + v[3551] * v[4001] + 2e0*(v[3560] * v[4001] + v[3590] * v[4003]) + dA[5] * v[4063]
			+ 2e0*dA[3] * v[4064] + v[3599] * v[426] + v[8044 + i3049];
		v[8157] = v[3554] - v[3558] + (dA[5] * v[3598]) / 2e0 + (dA[3] * v[3599]) / 2e0 - v[3550] * v[4001] + 2e0*(
			-4e0*v[1224] * v[3593] + v[3561] * v[4001]) + 2e0*dA[4] * (-v[3591] + v[3594] + v[4064]) + v[8026 + i3049];
		v[8158] = -v[3552] + v[3555] + 2e0*dA[5] * (-v[3591] + v[3603]) + v[3553] * v[4001] + 2e0*(v[3556] * v[4001]
			+ v[3595] * v[4003]) + dA[3] * v[4063] + v[3598] * v[426] + v[8008 + i3049];
		v[8159] = v[3320];
		v[8160] = v[3314];
		v[8161] = v[3308];
		v[8162] = -v[3544] + v[3546] + v[3538] * v[4005] + 2e0*(v[3547] * v[4005] + v[3580] * v[4007]) + dA[11] * v[4065]
			+ 2e0*dA[9] * v[4066] + v[3589] * v[432] + v[7990 + i3049];
		v[8163] = v[3541] - v[3545] + (dA[11] * v[3588]) / 2e0 + (dA[9] * v[3589]) / 2e0 - v[3537] * v[4005] + 2e0*(
			-4e0*v[1238] * v[3583] + v[3548] * v[4005]) + 2e0*dA[10] * (-v[3581] + v[3584] + v[4066]) + v[8062 + i3049];
		v[8164] = -v[3539] + v[3542] + 2e0*dA[11] * (-v[3581] + v[3610]) + v[3540] * v[4005] + 2e0*(v[3543] * v[4005]
			+ v[3585] * v[4007]) + dA[9] * v[4065] + v[3588] * v[432] + v[8080 + i3049];
		v[8165] = -v[3302];
		v[8166] = -v[3301];
		v[8167] = -v[3300];
		v[8168] = -v[3420] + v[3422] + v[3414] * v[4009] + 2e0*(v[3423] * v[4009] + v[3566] * v[4011]) + dB[5] * v[4067]
			+ 2e0*dB[3] * v[4068] + v[3575] * v[438] + v[8098 + i3049];
		v[8169] = v[3417] - v[3421] + (dB[5] * v[3574]) / 2e0 + (dB[3] * v[3575]) / 2e0 - v[3413] * v[4009] + 2e0*(
			-4e0*v[1252] * v[3569] + v[3424] * v[4009]) + 2e0*dB[4] * (-v[3567] + v[3570] + v[4068]) + v[8116 + i3049];
		v[8170] = -v[3415] + v[3418] + 2e0*dB[5] * (-v[3567] + v[3617]) + v[3416] * v[4009] + 2e0*(v[3419] * v[4009]
			+ v[3571] * v[4011]) + dB[3] * v[4067] + v[3574] * v[438] + v[8134 + i3049];
		Rc[i3049 - 1] += v[3080] * v[3157] + v[3113] * v[3158] + v[3092] * v[3222] + v[3091] * v[3223] + v[7878 + i3049];
		for (i3128 = 1; i3128 <= 18; i3128++) {
			Kc[i3049 - 1][i3128 - 1] += v[3565] * v[4838 + i3128] + v[3500] * v[4856 + i3128] + v[3579] * v[4874 + i3128]
				+ v[3576] * v[4892 + i3128] + v[8152 + i3128] + (*a4)*v[8170 + i3128];
		};/* end for */
	};/* end for */
#pragma endregion


}