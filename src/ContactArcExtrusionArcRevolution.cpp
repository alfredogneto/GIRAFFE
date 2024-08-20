#include "ContactArcExtrusionArcRevolution.h"
#include "Database.h"

//Variáveis globais
extern
Database db;

////////////////////////////////////////////////////////////////////
#define _CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include <crtdbg.h>
#ifdef _DEBUG
#define DBG_NEW new ( _NORMAL_BLOCK , __FILE__ , __LINE__ )
// Replace _NORMAL_BLOCK with _CLIENT_BLOCK if you want the
// allocations to be of _CLIENT_BLOCK type
#else
#define DBG_NEW new
#endif
////////////////////////////////////////////////////////////////////


ContactArcExtrusionArcRevolution::ContactArcExtrusionArcRevolution()
{
	index1 = 0;				//Body 1 - index
	index2 = 0;				//Body 2 - index
	sub_index1 = 0;			//Body 1 - sub_index
	sub_index2 = 0;			//Body 2 - sub_index
	invert = false;

	prev_active = false;
	cur_active = false;

	I3 = new Matrix(3, 3);
	(*I3)(0, 0) = 1.0;
	(*I3)(1, 1) = 1.0;
	(*I3)(2, 2) = 1.0;

	index1 = 0;				//index 1
	index2 = 0;				//index 2
	sub_index1 = 0;			//sub index 1
	sub_index2 = 0;			//sub index 2

	minimum_convective_range = 0.0;

	cd = NULL;
	Rc = NULL;
	Kc = NULL;

	alloc_control = false;

	GammaA = new Matrix(3);
	GammaB = new Matrix(3);

	eligible = false;
	prev_eligible = false;

	previous_evaluation = false;

	me = new double;

	interface_0_flag = false;
	interface_1_flag = false;
	inter_0 = NULL;
	inter_1 = NULL;

	DefaultValues();

	write_report = false;

	dA_zero = new double[12];
	for (int i = 0; i < 12; i++)
		dA_zero[i] = 0.0;
	dB_zero = new double[6];
	for (int i = 0; i < 6; i++)
		dB_zero[i] = 0.0;
}


ContactArcExtrusionArcRevolution::~ContactArcExtrusionArcRevolution()
{
	delete GammaA;
	delete GammaB;
	delete I3;
	delete me;
	delete[]dA_zero;
	delete[]dB_zero;
	Free();
}

void ContactArcExtrusionArcRevolution::Alloc()
{

	nDOF = sA->nDOFs + sB->nDOFs;

	if (alloc_control == false)
	{
		cd = new SSContactData();
		cd->n_solutions = 1;
		cd->Alloc();

		fn = DBG_NEW double[3];
		ft = DBG_NEW double[3];

		Rc = DBG_NEW double[nDOF];
		Kc = DBG_NEW double*[nDOF];
		for (int i = 0; i < nDOF; i++)
			Kc[i] = DBG_NEW double[nDOF];
		for (int ni = 0; ni < nDOF; ni++)
			for (int nj = 0; nj < nDOF; nj++)
				Kc[ni][nj] = 0.0;

		cAp = DBG_NEW double[2];
		cBp = DBG_NEW double[2];
		cAi = DBG_NEW double[2];
		cBi = DBG_NEW double[2];

		convective_range = new Matrix(4);
		convective_max = new Matrix(4);
		convective_min = new Matrix(4);

		write_report = db.execution_data->print_contact_report;

		alloc_control = true;
	}
}

void ContactArcExtrusionArcRevolution::Free()
{
	if (alloc_control == true)
	{
		delete cd;

		delete[] fn;
		delete[] ft;

		delete[] Rc;
		for (int i = 0; i < nDOF; i++)
			delete[]Kc[i];
		delete[]Kc;

		delete[]cAp;
		delete[]cBp;
		delete[]cAi;
		delete[]cBi;

		delete convective_range;
		delete convective_max;
		delete convective_min;

		alloc_control = false;
	}
}

//Chute inicial para coordenadas convectivas do par de superfícies
void ContactArcExtrusionArcRevolution::InitialGuess()
{
	//Local specific pointers
	ArcExtrusion *sA_local;
	ArcRevolution *sB_local;
	if (invert)
	{
		sA_local = static_cast<ArcExtrusion*>(db.body_geometries[index2]->ptr_geom[sub_index2]);
		sB_local = static_cast<ArcRevolution*>(db.body_geometries[index1]->ptr_geom[sub_index1]);
	}
	else
	{
		sA_local = static_cast<ArcExtrusion*>(db.body_geometries[index1]->ptr_geom[sub_index1]);
		sB_local = static_cast<ArcRevolution*>(db.body_geometries[index2]->ptr_geom[sub_index2]);
	}
	

	//Determinação de zeta_barra
	//b,t calculados a partir das posições nodais do arco extrudado (eq. reta). Centro do arco revolucionado é projetado nesta reta, estimando o zeta_barra
	Matrix b = 0.5*(*sA_local->x_Ap + *sA_local->x_Bp);
	Matrix t = 0.5*(*sA_local->x_Bp - *sA_local->x_Ap);
	Matrix center = *sB_local->x_Ap + (*sB_local->Q_Ap) * (*sB_local->center_local);
	double zeta_barra = (dot(t, center) - dot(t, b)) / (dot(t, t));
	Matrix x_barra = b + zeta_barra * t;
	
	//Determinação do phi_barra
	Matrix dir = x_barra - center;
	Matrix e1(3);
	Matrix e2(3);
	Matrix e3(3);
	//Obs: eixos já escritos no sistema local
	e1(0, 0) = 1.0;
	e2(1, 0) = 1.0;
	e3(2, 0) = 1.0;
	//Orientações atuais dos eixos e1, e2 e e3
	e1 = (*sB_local->Q_Ap)*e1;
	e2 = (*sB_local->Q_Ap)*e2;
	e3 = (*sB_local->Q_Ap)*e3;
	dir = dir - dot(dir, e2)*e2;	//projetando dir no plano e1e3
	dir = (1.0 / norm(dir))*dir;	//normalização de dir
	double e1_dir = dot(e1, dir);
	double e3_dir = dot(e3, dir);
	//Testes para evitar erros de arredondamento
	if (e1_dir > 1.0)
		e1_dir = 1.0;
	if (e1_dir < -1.0)
		e1_dir = -1.0;
	if (e3_dir > 1.0)
		e3_dir = 1.0;
	if (e3_dir < -1.0)
		e3_dir = -1.0;
	double phi_barra;
	if (e1_dir >= 0.0 && e3_dir <= 0.0)
		phi_barra = acos(e1_dir);
	else
	{
		if (e1_dir >= 0.0 && e3_dir > 0.0)
			phi_barra = 2 * PI - acos(e1_dir);
		else
		{
			if (e1_dir < 0.0 && e3_dir <= 0.0)
				phi_barra =  acos(e1_dir);
			else
				phi_barra = PI + acos(-e1_dir);
		}
	}

	/*if (sA_local->number == 19 && sB_local->number == 5)
	{
		center.print();
		x_barra.print();
		e1.print();
		e3.print();
	}*/

	//Preenchendo as coordenadas convectivas:
	cd->convective[0][0] = zeta_barra;
	cd->convective[0][1] = db.arcs[sA_local->arc_ID - 1]->CentralTheta();
	cd->convective[0][2] = phi_barra;
	cd->convective[0][3] = db.arcs[sB_local->arc_ID - 1]->CentralTheta();
}

void ContactArcExtrusionArcRevolution::SetVariables()
{
	//Local specific pointers
	ArcExtrusion *sA_local;
	ArcRevolution *sB_local;
	if (invert)
	{
		sA_local = static_cast<ArcExtrusion*>(db.body_geometries[index2]->ptr_geom[sub_index2]);
		sB_local = static_cast<ArcRevolution*>(db.body_geometries[index1]->ptr_geom[sub_index1]);
	}
	else
	{
		sA_local = static_cast<ArcExtrusion*>(db.body_geometries[index1]->ptr_geom[sub_index1]);
		sB_local = static_cast<ArcRevolution*>(db.body_geometries[index2]->ptr_geom[sub_index2]);
	}

	radA = &db.arcs[sA_local->arc_ID - 1]->radius;
	cpointA = db.arcs[sA_local->arc_ID - 1]->c_point.getMatrix();

	radB = &db.arcs[sB_local->arc_ID - 1]->radius;
	cpointB = db.arcs[sB_local->arc_ID - 1]->c_point.getMatrix();
	
	dA = sA->d->getMatrix();
	dB = sB->d->getMatrix();
	duiA = sA->dui->getMatrix();
	dduiA = sA->ddui->getMatrix();
	duiB = sB->dui->getMatrix();
	dduiB = sB->ddui->getMatrix();
	xAAi = sA_local->xAi;
	xBAi = sA_local->xBi;
	xABi = sB_local->xAi;
	
	normalintA = &sA_local->flag_normal_int;
	normalintB = &sB_local->flag_normal_int;
	QAAi = sA_local->QAi;
	QBAi = sA_local->QBi;
	QABi = sB_local->QAi;
	
	gti = cd->copy_g_t[0]->getMatrix();
	gtpupdated = cd->g_t[0]->getMatrix();
	stick = &cd->copy_stick[0];
	stickupdated = &cd->stick[0];
	invH = cd->invHessian[0];

	interfacelaw0 = &interface_0_flag;
}

void ContactArcExtrusionArcRevolution::PrintAceGenPointers()
{
	printf("%lf\n", *radA);
	printf("%lf\n", *radB);
	db.PrintPtr(cpointA, 2);
	db.PrintPtr(cpointB, 2);
	
	db.PrintPtr(dA, sA->nDOFs);
	db.PrintPtr(dB, sB->nDOFs);

	db.PrintPtr(duiA, sA->nDOFs);
	db.PrintPtr(duiB, sB->nDOFs);

	db.PrintPtr(dduiA, sA->nDOFs);
	db.PrintPtr(dduiB, sB->nDOFs);

	db.PrintPtr(xAAi, 3);
	db.PrintPtr(xBAi, 3);

	db.PrintPtr(xABi, 3);

	db.PrintPtr(QAAi, 3, 3);
	db.PrintPtr(QBAi, 3, 3);

	db.PrintPtr(QABi, 3, 3);

	db.PrintPtr(invH, 4, 4);

	cd->Plot();
}

void ContactArcExtrusionArcRevolution::Report()
{
	if (db.execution_data->print_contact_report)
	{
		db.myprintf("ContactArcExtrusionArcRevolution between Surface\t%d and Surface\t%d\n", sA->number, sB->number);
		Matrix cb(4);
		cb(0, 0) = cd->convective[0][0];
		cb(1, 0) = cd->convective[0][1];
		cb(2, 0) = cd->convective[0][2];
		cb(3, 0) = cd->convective[0][3];
		//int charac = CharacterizeCriticalPoint(&cb);
		//db.myprintf("Characterization %d\n", charac);
		cd->PlotSmallReport();
		if (eligible)
		{
			db.myprintf("\t\tThis contact is ELIGIBLE!\n");
			if (interface_1_flag)
			{
				if (cd->g_n[0] < *gnbb)
					db.myprintf("\t\tBarrier activated! ---------- %.1f %c of contact layer is active.\n", 100.0*(1.0 - cd->g_n[0] / (*gnb)), 37);
			}
		}
	}
}


void ContactArcExtrusionArcRevolution::InitializeConvectiveRange()
{
	//Local specific pointers
	ArcExtrusion *sA_local;
	ArcRevolution *sB_local;
	if (invert)
	{
		sA_local = static_cast<ArcExtrusion*>(db.body_geometries[index2]->ptr_geom[sub_index2]);
		sB_local = static_cast<ArcRevolution*>(db.body_geometries[index1]->ptr_geom[sub_index1]);
	}
	else
	{
		sA_local = static_cast<ArcExtrusion*>(db.body_geometries[index1]->ptr_geom[sub_index1]);
		sB_local = static_cast<ArcRevolution*>(db.body_geometries[index2]->ptr_geom[sub_index2]);
	}

	//Degeneration
	if (sA->deg_u1 == true)
	{
		(*convective_min)(0, 0) = sA->deg_u1_value;
		(*convective_max)(0, 0) = sA->deg_u1_value;
		(*convective_range)(0, 0) = 0.0;
	}
	//No degeneration
	else
	{
		(*convective_min)(0, 0) = -1.0;
		(*convective_max)(0, 0) = +1.0;
		(*convective_range)(0, 0) = +2.0;
	}

	//Degeneration
	if (sA->deg_u2 == true)
	{
		(*convective_min)(1, 0) = sA->deg_u2_value;
		(*convective_max)(1, 0) = sA->deg_u2_value;
		(*convective_range)(1, 0) = 0.0;
	}

	//No degeneration
	else
	{
		(*convective_min)(1, 0) = db.arcs[sA_local->arc_ID - 1]->theta_i;
		(*convective_max)(1, 0) = db.arcs[sA_local->arc_ID - 1]->theta_f;
		(*convective_range)(1, 0) = db.arcs[sA_local->arc_ID - 1]->AngularRange();
	}

	//Degeneration
	if (sB->deg_u1 == true)
	{
		(*convective_min)(2, 0) = sB->deg_u1_value;
		(*convective_max)(2, 0) = sB->deg_u1_value;
		(*convective_range)(2, 0) = 0.0;
	}
	//No degeneration
	else
	{
		(*convective_min)(2, 0) = 0;
		(*convective_max)(2, 0) = +2.0 * PI;
		(*convective_range)(2, 0) = + 2.0 * PI;
	}

	//Degeneration
	if (sB->deg_u2 == true)
	{
		(*convective_min)(3, 0) = sB->deg_u2_value;
		(*convective_max)(3, 0) = sB->deg_u2_value;
		(*convective_range)(3, 0) = 0.0;
	}

	//No degeneration
	else
	{
		(*convective_min)(3, 0) = db.arcs[sB_local->arc_ID - 1]->theta_i;
		(*convective_max)(3, 0) = db.arcs[sB_local->arc_ID - 1]->theta_f;
		(*convective_range)(3, 0) = db.arcs[sB_local->arc_ID - 1]->AngularRange();
	}

	//Setting the minimum convective range
	minimum_convective_range = 1e100;
	for (int i = 0; i < 4; i++)
	{
		if ((*convective_range)(i, 0) < minimum_convective_range && (*convective_range)(i, 0) != 0.0)
			minimum_convective_range = (*convective_range)(i, 0);
	}
}

//Verifica range de coordenadas convectivas
int ContactArcExtrusionArcRevolution::VerifyConvectiveRange(Matrix& mc)
{
	//Local specific pointers
	ArcExtrusion *sA_local;
	ArcRevolution *sB_local;
	if (invert)
	{
		sA_local = static_cast<ArcExtrusion*>(db.body_geometries[index2]->ptr_geom[sub_index2]);
		sB_local = static_cast<ArcRevolution*>(db.body_geometries[index1]->ptr_geom[sub_index1]);
	}
	else
	{
		sA_local = static_cast<ArcExtrusion*>(db.body_geometries[index1]->ptr_geom[sub_index1]);
		sB_local = static_cast<ArcRevolution*>(db.body_geometries[index2]->ptr_geom[sub_index2]);
	}

	//Retornos:
	//0 - Range fisico da superficie
	//3 - Fora do range físico da superfície e em região proibida de parâmetro
	//4 - Fora do range físico da superficie

	//if (sA_local->number == 7 && sB_local->number == 2)
	//	int test = 1;

	//Primeiro teste - região proibitiva de parâmetro theta - relacionada à revolução (theta_2)
	if (ArcReduction(mc(3, 0)) < sB_local->theta_valid_min || ArcReduction(mc(3, 0)) > sB_local->theta_valid_max)
		return 3;
	
	if (!(abs(mc(0, 0)) < 1.0 || sA->deg_u1))
		return 4;
	//if (!(abs(mc(2, 0) < 2 * PI && mc(2, 0) > 0.0) || sB->deg_u1))
	//	return 4;
	if (!(db.arcs[sA_local->arc_ID - 1]->InsideArc(mc(1, 0)) || sA->deg_u2))
		return 4;
	if (!(db.arcs[sB_local->arc_ID - 1]->InsideArc(mc(3, 0)) || sB->deg_u2))
		return 4;
	
	//If no return hapened until here, the point is on the valid range for all conditions
	return 0;
}

void ContactArcExtrusionArcRevolution::CompactReport()
{
	db.myprintf("Eligible %d\n", (int)eligible);
	db.myprintf("Gap %.6e\n", cd->g_n[0]);
}

//Calcula a função objetivo para um conjunto de coordenadas convectivas - Phase 1
double ContactArcExtrusionArcRevolution::ObjectivePhase1(Matrix& mc)
{
	double v[2000];
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
	v[185] = v[187] * cos(c[2]);
	v[186] = cpointB[1] + (*radB)*sin(c[3]);
	v[188] = -(v[187] * sin(c[2]));
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

	return Ob;
}

//Calcula o Gradiente da função objetivo - Phase 1
void ContactArcExtrusionArcRevolution::GradientPhase1(Matrix& mc, Matrix& mGra)
{
	double v[2000];
	double *c = mc.getMatrix();
	double Gra[4];

	//ACEGEN
	int i246;
	v[215] = cos(c[2]);
	v[258] = v[215];
	v[210] = sin(c[2]);
	v[257] = -(v[210]);
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
	v[281] = v[187] * (-(v[210] * v[250] ) - v[215] * v[249] );
	v[282] = v[209] * (-(v[174] * v[242]) - v[177] * v[243] - v[180] * v[248]) - v[206] * (v[249] * v[257] + v[250] * v[258]);

	for (i246 = 1; i246 <= 4; i246++) {
		Gra[i246 - 1] = v[278 + i246];
	};/* end for */

	for (int i = 0; i < 4; i++)
		mGra(i, 0) = Gra[i];
}

//Calcula a Hessiana da função objetivo - Phase 1
void ContactArcExtrusionArcRevolution::HessianPhase1(Matrix& mc, Matrix& mHes)
{
	double v[2000];
	double *c = mc.getMatrix();
	double Hes[4][4];

	//ACEGEN
	int i246, i254, i295, b287, b297;
	v[215] = cos(c[2]);
	v[294] = v[215] ;
	v[210] = sin(c[2]);
	v[303] = v[210] ;
	v[293] = v[210] ;
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
	v[304] = v[187] ;
	v[292] = v[187] * v[215];
	v[216] = -(v[187] * v[303]);
	v[344] = 0e0;
	v[345] = 0e0;
	v[346] = v[216];
	v[347] = -(v[206] * v[294]);
	v[217] = -(v[292] );
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
	v[185] = v[292] ;
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
	v[305] = v[249] ;
	v[250] = -(v[173] * v[242]) - v[176] * v[243] - v[179] * v[248];
	v[314] = v[206] * (v[250] * v[303] + v[215] * v[305]);
	v[420] = 0e0;
	v[421] = 0e0;
	v[422] = -(v[187] * v[305]);
	v[423] = -(v[206] * v[250] );
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
		v[426] = -(v[210] * (v[273] * v[304] + v[419 + i246])) + v[215] * (v[415 + i246] - v[187] * v[272] );
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
}

//Calcula e rotorna o gap (com sinal)
double ContactArcExtrusionArcRevolution::Gap(Matrix& mc, bool fixed_normals, Matrix& nA, Matrix& nB)
{
	double v[2000];		//variável temporária - AceGen
	double gap;
	double *c = mc.getMatrix();

	bool* fixnormal = &fixed_normals;
	double* normalA = nA.getMatrix();
	double* normalB = nB.getMatrix();

	//ACEGEN
	int b237;
	v[215] = cos(c[2]);
	v[210] = sin(c[2]);
	v[250] = v[210] ;
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
	v[216] = -(v[187] * v[210] );
	v[208] = -(v[206] * v[215] );
	v[217] = -(v[249] );
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
	v[185] = v[249] ;
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

	return gap;
}

//Calcula o Gradiente do gap
void ContactArcExtrusionArcRevolution::GradientGap(Matrix& mc, Matrix& mGra, bool fixed_normals, Matrix& nA, Matrix& nB)
{
	double v[2000];		//variável temporária - AceGen
	double *c = mc.getMatrix();
	double Gra[4];
	
	bool* fixnormal = &fixed_normals;
	double* normalA = nA.getMatrix();
	double* normalB = nB.getMatrix();

	//ACEGEN
	int i246, i312, i313, b237, b248;
	v[215] = cos(c[2]);
	v[317] = -(v[215] );
	v[210] = sin(c[2]);
	v[306] = v[210] ;
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
	v[216] = -(v[187] * v[210] );
	v[208] = v[206] * v[317];
	v[217] = -(v[305] );
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
	v[185] = v[305] ;
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
	v[318] = (-(v[179] * v[257]) - v[176] * v[258] - v[173] * v[259]) - (v[181] * v[269] + v[178] * v[271]
		+ v[175] * v[273]);
	v[319] = -((v[179] * v[269] + v[176] * v[271] + v[173] * v[273])) - (-(v[181] * v[257]) - v[178] * v[258]
		- v[175] * v[259]);
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
	v[370] = -(v[210] * (v[187] * v[318] - v[206] * v[275] )) + v[215] * (v[187] * v[319] + v[206] * v[274] 
		);
	v[371] = v[209] * (-(v[180] * v[257]) - v[177] * v[258] - v[174] * v[259] + v[274] * v[306] + v[275] * v[317]) - v[206] *
		(v[180] * v[268] + v[177] * v[270] + v[174] * v[272] + v[215] * v[318] + v[210] * v[319]);

	for (i246 = 1; i246 <= 4; i246++) {
		Gra[i246 - 1] = v[367 + i246];
	};/* end for */

	for (int i = 0; i < 4; i++)
		mGra(i, 0) = Gra[i];
}

//Calcula a Hessiana do gap
void ContactArcExtrusionArcRevolution::HessianGap(Matrix& mc, Matrix& mHes, bool fixed_normals, Matrix& nA, Matrix& nB)
{
	double v[2000];		//variável temporária - AceGen
	double *c = mc.getMatrix();
	double Hes[4][4];

	bool* fixnormal = &fixed_normals;
	double* normalA = nA.getMatrix();
	double* normalB = nB.getMatrix();

	//ACEGEN
	int i246, i302, i438, i441, i450, b237, b248, b361, b421, b462;
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
	v[216] = v[427] ;
	v[208] = -(v[315] );
	v[217] = -(v[426] );
	v[211] = v[314] ;
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
	v[185] = v[426];
	v[186] = cpointB[1] + v[206];
	v[188] = v[427] ;
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
	v[461] = v[274] ;
	v[275] = v[179] * v[268] + v[176] * v[270] + v[173] * v[272];
	v[460] = -(v[275] );
	v[416] = -(v[180] * v[257]) - v[177] * v[258] - v[174] * v[259] + v[215] * v[460] + v[210] * v[461];
	v[276] = v[181] * v[269] + v[178] * v[271] + v[175] * v[273];
	v[436] = -(v[276] );
	v[413] = v[187] * v[436] + (v[187] * v[267] - v[206] * v[275]);
	v[412] = v[436] + v[267] ;
	v[277] = v[179] * v[269] + v[176] * v[271] + v[173] * v[273];
	v[437] = -(v[277] );
	v[415] = v[180] * v[268] + v[177] * v[270] + v[174] * v[272] + (v[215] * v[267] - v[210] * v[277]) -
		(v[210] * v[266] + v[215] * v[276]);
	v[731] = 0e0;
	v[732] = 0e0;
	v[733] = 0e0;
	v[734] = -v[415];
	v[414] = v[437] - v[266] ;
	v[482] = -(v[210] * v[412]) + v[215] * v[414];
	v[747] = 0e0;
	v[748] = 0e0;
	v[749] = -v[413];
	v[750] = -(v[206] * v[414]) + v[209] * v[461];
	v[411] = v[187] * v[437] + (-(v[187] * v[266]) + v[206] * v[274]);
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
		v[447] = v[634 + i246] ;
		v[446] = v[626 + i246] ;
		v[311] = v[618 + i246];
		v[444] = v[311] ;
		v[310] = v[610 + i246];
		v[443] = v[310] ;
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
			v[452] = v[311] ;
			v[451] = -(v[310] );
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
		v[757] = v[215] * (v[746 + i246] - v[187] * v[380]  + (-(v[187] * v[369]) + v[206] * v[377])) - v[210] *
			(v[750 + i246] + (v[187] * v[370] - v[206] * v[378]) - v[187] * v[379] );
		v[758] = v[209] * (-(v[174] * v[362]) - v[177] * v[363] - v[180] * v[364] + v[730 + i246] - v[215] * v[378] 
			+ v[210] * v[377] ) - v[206] * ((b462 ? v[416] : 0e0) + (i246 == 3 ? v[482] : 0e0) + v[180] * v[371] + v[177] * v[373]
				+ v[174] * v[375] + (v[215] * v[370] - v[210] * v[380]) - (v[210] * v[369] + v[215] * v[379]));

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
}

void ContactArcExtrusionArcRevolution::MountLocalContributions()
{
	v = DBG_NEW double[21000];
	cAp[0] = cd->convective[0][0];
	cAp[1] = cd->convective[0][1];
	cBp[0] = cd->convective[0][2];
	cBp[1] = cd->convective[0][3];

	cAi[0] = cd->copy_convective[0][0];
	cAi[1] = cd->copy_convective[0][1];
	cBi[0] = cd->copy_convective[0][2];
	cBi[1] = cd->copy_convective[0][3];

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
	for (int i = 0; i < nDOF; i++)
	{
		Rc[i] = 0.0;
		for (int j = 0; j < nDOF; j++)
			Kc[i][j] = 0.0;
	}

	bool *previouscontact = &prev_eligible;

	//Avalia contribuições de contato

	//ACEGEN
	double v01; double v010; double v011; double v012; double v013; double v014;
	double v015; double v016; double v017; double v018; double v019; double v02;
	double v020; double v021; double v022; double v023; double v024; double v025;
	double v026; double v027; double v028; double v029; double v03; double v030;
	double v031; double v032; double v033; double v04; double v05; double v06;
	double v07; double v08; double v09;
	int i1149, i1284, i1687, i2526, i4521, i5016, i6660, i6661, i6662, i6663, i6664, i6665
		, i6673, i6674, i6675, i6957, i6958, i6959, i6960, i6961, i6962, i6963, i6990, i6991, i7003
		, i7004, b4, b5, b6, b382, b418, b419, b941, b942, b974, b1031, b1040, b1041, b1049, b1050
		, b1051, b1052, b1055, b1056, b1064, b1065, b1078, b1079, b1080, b1081, b1097, b1114, b1115
		, b1127, b1157, b1161, b1162, b1166, b1168, b1413, b1427, b1431, b1432, b1440, b1444, b1692
		, b1693, b1694, b1695, b1782, b1783, b1787, b1791, b1792, b1806, b1807, b2033, b2073, b2126
		, b2133, b3008, b3022, b3083, b3129, b3583, b3587, b3589, b3590, b3595, b3596, b3622, b3623
		, b3624, b3640, b3771, b3772, b3782, b3786, b3790, b3804, b3805, b3921, b3952, b4021, b4022
		, b4536, b4537, b4545, b4546, b4550, b4551, b4552, b4684, b4688, b5437, b5454, b5729, b5730
		, b5731, b5745, b5746, b5756, b5760, b5764, b5775, b5776, b5919, b5920, b6526, b6527, b6528
		, b6529, b6530, b6531, b6992, b7005;
	v[1] = gti[0];
	v[2] = gti[1];
	v[3] = gti[2];
	b4 = (*previouscontact);
	b5 = (*stick);
	b6 = (*interfacelaw0);
	v[7] = (*a4);
	v[8] = (*a5);
	v[9] = (*a6);
	v[10] = invH[0][0];
	v[11] = invH[0][1];
	v[12] = invH[0][2];
	v[13] = invH[0][3];
	v[14] = invH[1][0];
	v[15] = invH[1][1];
	v[16] = invH[1][2];
	v[17] = invH[1][3];
	v[18] = invH[2][0];
	v[19] = invH[2][1];
	v[20] = invH[2][2];
	v[21] = invH[2][3];
	v[22] = invH[3][0];
	v[23] = invH[3][1];
	v[24] = invH[3][2];
	v[25] = invH[3][3];
	v[26] = (*epsn1);
	v[27] = (*gnb);
	v[28] = (*gnbb);
	v[29] = (*n1);
	v[6506] = v[26] * v[29];
	v[1045] = -1e0 + v[29];
	v[3599] = -2e0 + v[1045];
	v[3585] = -1e0 + v[1045];
	v[3584] = v[1045] * v[6506];
	v[30] = (*n2);
	v[31] = (*zetan);
	v[6734] = 2e0*v[31];
	v[32] = (*epst);
	v[33] = (*ct);
	v[34] = (*mus);
	v[35] = (*mud);
	v[36] = (*me);
	v[6726] = v[31] * v[36];
	v[6426] = v[3584] * v[36];
	v[5754] = v[31] * v[6426];
	v[3794] = v[6426] / 2e0;
	v[37] = dA[0];
	v[38] = dA[1];
	v[39] = dA[2];
	v[40] = dA[3];
	v[459] = v[40] / 2e0;
	v[457] = 2e0*v[40];
	v[6874] = -0.5e0*v[457];
	v[167] = (v[40] * v[40]);
	v[41] = dA[4];
	v[460] = 2e0*v[41];
	v[9953] = 0e0;
	v[9954] = 0e0;
	v[9955] = 0e0;
	v[9956] = v[457];
	v[9957] = v[460];
	v[9958] = 0e0;
	v[9959] = 0e0;
	v[9960] = 0e0;
	v[9961] = 0e0;
	v[9962] = 0e0;
	v[9963] = 0e0;
	v[9964] = 0e0;
	v[9965] = 0e0;
	v[9966] = 0e0;
	v[9967] = 0e0;
	v[9968] = 0e0;
	v[9969] = 0e0;
	v[9970] = 0e0;
	v[6875] = -0.5e0*v[460];
	v[10447] = 0e0;
	v[10448] = 0e0;
	v[10449] = 0e0;
	v[10450] = -v[457];
	v[10451] = -v[460];
	v[10452] = 0e0;
	v[10453] = 0e0;
	v[10454] = 0e0;
	v[10455] = 0e0;
	v[10456] = 0e0;
	v[10457] = 0e0;
	v[10458] = 0e0;
	v[10459] = 0e0;
	v[10460] = 0e0;
	v[10461] = 0e0;
	v[10462] = 0e0;
	v[10463] = 0e0;
	v[10464] = 0e0;
	v[458] = v[41] / 2e0;
	v[9233] = 0e0;
	v[9234] = 0e0;
	v[9235] = 0e0;
	v[9236] = v[458];
	v[9237] = v[459];
	v[9238] = 0e0;
	v[9239] = 0e0;
	v[9240] = 0e0;
	v[9241] = 0e0;
	v[9242] = 0e0;
	v[9243] = 0e0;
	v[9244] = 0e0;
	v[9245] = 0e0;
	v[9246] = 0e0;
	v[9247] = 0e0;
	v[9248] = 0e0;
	v[9249] = 0e0;
	v[9250] = 0e0;
	v[165] = v[40] * v[458];
	v[160] = (v[41] * v[41]);
	v[515] = -v[160] - v[167];
	v[6623] = v[515] / 2e0;
	v[42] = dA[5];
	v[493] = v[165] + v[42];
	v[6641] = -2e0*v[493];
	v[483] = v[165] - v[42];
	v[6643] = -2e0*v[483];
	v[462] = 2e0*v[42];
	v[9719] = 0e0;
	v[9720] = 0e0;
	v[9721] = 0e0;
	v[9722] = v[457];
	v[9723] = 0e0;
	v[9724] = v[462];
	v[9725] = 0e0;
	v[9726] = 0e0;
	v[9727] = 0e0;
	v[9728] = 0e0;
	v[9729] = 0e0;
	v[9730] = 0e0;
	v[9731] = 0e0;
	v[9732] = 0e0;
	v[9733] = 0e0;
	v[9734] = 0e0;
	v[9735] = 0e0;
	v[9736] = 0e0;
	v[9485] = 0e0;
	v[9486] = 0e0;
	v[9487] = 0e0;
	v[9488] = 0e0;
	v[9489] = v[460];
	v[9490] = v[462];
	v[9491] = 0e0;
	v[9492] = 0e0;
	v[9493] = 0e0;
	v[9494] = 0e0;
	v[9495] = 0e0;
	v[9496] = 0e0;
	v[9497] = 0e0;
	v[9498] = 0e0;
	v[9499] = 0e0;
	v[9500] = 0e0;
	v[9501] = 0e0;
	v[9502] = 0e0;
	v[6873] = -0.5e0*v[462];
	v[10519] = 0e0;
	v[10520] = 0e0;
	v[10521] = 0e0;
	v[10522] = 0e0;
	v[10523] = -v[460];
	v[10524] = -v[462];
	v[10525] = 0e0;
	v[10526] = 0e0;
	v[10527] = 0e0;
	v[10528] = 0e0;
	v[10529] = 0e0;
	v[10530] = 0e0;
	v[10531] = 0e0;
	v[10532] = 0e0;
	v[10533] = 0e0;
	v[10534] = 0e0;
	v[10535] = 0e0;
	v[10536] = 0e0;
	v[9215] = 0e0;
	v[9216] = 0e0;
	v[9217] = 0e0;
	v[9218] = v[457];
	v[9219] = v[460];
	v[9220] = v[462];
	v[9221] = 0e0;
	v[9222] = 0e0;
	v[9223] = 0e0;
	v[9224] = 0e0;
	v[9225] = 0e0;
	v[9226] = 0e0;
	v[9227] = 0e0;
	v[9228] = 0e0;
	v[9229] = 0e0;
	v[9230] = 0e0;
	v[9231] = 0e0;
	v[9232] = 0e0;
	v[461] = v[42] / 2e0;
	v[9251] = 0e0;
	v[9252] = 0e0;
	v[9253] = 0e0;
	v[9254] = v[461];
	v[9255] = 0e0;
	v[9256] = v[459];
	v[9257] = 0e0;
	v[9258] = 0e0;
	v[9259] = 0e0;
	v[9260] = 0e0;
	v[9261] = 0e0;
	v[9262] = 0e0;
	v[9263] = 0e0;
	v[9264] = 0e0;
	v[9265] = 0e0;
	v[9266] = 0e0;
	v[9267] = 0e0;
	v[9268] = 0e0;
	v[9269] = 0e0;
	v[9270] = 0e0;
	v[9271] = 0e0;
	v[9272] = 0e0;
	v[9273] = v[461];
	v[9274] = v[458];
	v[9275] = 0e0;
	v[9276] = 0e0;
	v[9277] = 0e0;
	v[9278] = 0e0;
	v[9279] = 0e0;
	v[9280] = 0e0;
	v[9281] = 0e0;
	v[9282] = 0e0;
	v[9283] = 0e0;
	v[9284] = 0e0;
	v[9285] = 0e0;
	v[9286] = 0e0;
	v[172] = v[41] * v[461];
	v[510] = v[172] + v[40];
	v[6638] = -2e0*v[510];
	v[502] = v[172] - v[40];
	v[6640] = -2e0*v[502];
	v[170] = v[40] * v[461];
	v[506] = v[170] - v[41];
	v[6639] = -2e0*v[506];
	v[489] = v[170] + v[41];
	v[6642] = -2e0*v[489];
	v[161] = (v[42] * v[42]);
	v[1584] = 4e0 + v[160] + v[161] + v[167];
	v[7199] = 24e0 / Power(v[1584], 4);
	v[6427] = 8e0 / Power(v[1584], 3);
	v[2587] = -(v[457] * v[6427]);
	v[2585] = v[460] * v[6427];
	v[2582] = -(v[462] * v[6427]);
	v[1351] = 1e0 / (v[1584] * v[1584]);
	v[6428] = 4e0*v[1351];
	v[497] = -v[161] - v[167];
	v[6624] = v[497] / 2e0;
	v[478] = -v[160] - v[161];
	v[6625] = v[478] / 2e0;
	v[477] = v[462] * v[6428];
	v[6436] = -0.5e0*v[477];
	v[519] = v[515] * v[6436];
	v[476] = -(v[460] * v[6428]);
	v[6433] = v[476] / 2e0;
	v[499] = v[497] * v[6433];
	v[475] = v[457] * v[6428];
	v[6435] = -0.5e0*v[475];
	v[479] = v[478] * v[6435];
	v[43] = dA[6];
	v[44] = dA[7];
	v[45] = dA[8];
	v[46] = dA[9];
	v[465] = v[46] / 2e0;
	v[463] = 2e0*v[46];
	v[6860] = -0.5e0*v[463];
	v[186] = (v[46] * v[46]);
	v[47] = dA[10];
	v[466] = 2e0*v[47];
	v[10007] = 0e0;
	v[10008] = 0e0;
	v[10009] = 0e0;
	v[10010] = 0e0;
	v[10011] = 0e0;
	v[10012] = 0e0;
	v[10013] = 0e0;
	v[10014] = 0e0;
	v[10015] = 0e0;
	v[10016] = v[463];
	v[10017] = v[466];
	v[10018] = 0e0;
	v[10019] = 0e0;
	v[10020] = 0e0;
	v[10021] = 0e0;
	v[10022] = 0e0;
	v[10023] = 0e0;
	v[10024] = 0e0;
	v[6861] = -0.5e0*v[466];
	v[10555] = 0e0;
	v[10556] = 0e0;
	v[10557] = 0e0;
	v[10558] = 0e0;
	v[10559] = 0e0;
	v[10560] = 0e0;
	v[10561] = 0e0;
	v[10562] = 0e0;
	v[10563] = 0e0;
	v[10564] = -v[463];
	v[10565] = -v[466];
	v[10566] = 0e0;
	v[10567] = 0e0;
	v[10568] = 0e0;
	v[10569] = 0e0;
	v[10570] = 0e0;
	v[10571] = 0e0;
	v[10572] = 0e0;
	v[464] = v[47] / 2e0;
	v[9305] = 0e0;
	v[9306] = 0e0;
	v[9307] = 0e0;
	v[9308] = 0e0;
	v[9309] = 0e0;
	v[9310] = 0e0;
	v[9311] = 0e0;
	v[9312] = 0e0;
	v[9313] = 0e0;
	v[9314] = v[464];
	v[9315] = v[465];
	v[9316] = 0e0;
	v[9317] = 0e0;
	v[9318] = 0e0;
	v[9319] = 0e0;
	v[9320] = 0e0;
	v[9321] = 0e0;
	v[9322] = 0e0;
	v[184] = v[46] * v[464];
	v[179] = (v[47] * v[47]);
	v[560] = -v[179] - v[186];
	v[6605] = v[560] / 2e0;
	v[48] = dA[11];
	v[538] = v[184] + v[48];
	v[6635] = -2e0*v[538];
	v[528] = v[184] - v[48];
	v[6637] = -2e0*v[528];
	v[468] = 2e0*v[48];
	v[9773] = 0e0;
	v[9774] = 0e0;
	v[9775] = 0e0;
	v[9776] = 0e0;
	v[9777] = 0e0;
	v[9778] = 0e0;
	v[9779] = 0e0;
	v[9780] = 0e0;
	v[9781] = 0e0;
	v[9782] = v[463];
	v[9783] = 0e0;
	v[9784] = v[468];
	v[9785] = 0e0;
	v[9786] = 0e0;
	v[9787] = 0e0;
	v[9788] = 0e0;
	v[9789] = 0e0;
	v[9790] = 0e0;
	v[9539] = 0e0;
	v[9540] = 0e0;
	v[9541] = 0e0;
	v[9542] = 0e0;
	v[9543] = 0e0;
	v[9544] = 0e0;
	v[9545] = 0e0;
	v[9546] = 0e0;
	v[9547] = 0e0;
	v[9548] = 0e0;
	v[9549] = v[466];
	v[9550] = v[468];
	v[9551] = 0e0;
	v[9552] = 0e0;
	v[9553] = 0e0;
	v[9554] = 0e0;
	v[9555] = 0e0;
	v[9556] = 0e0;
	v[6859] = -0.5e0*v[468];
	v[10627] = 0e0;
	v[10628] = 0e0;
	v[10629] = 0e0;
	v[10630] = 0e0;
	v[10631] = 0e0;
	v[10632] = 0e0;
	v[10633] = 0e0;
	v[10634] = 0e0;
	v[10635] = 0e0;
	v[10636] = 0e0;
	v[10637] = -v[466];
	v[10638] = -v[468];
	v[10639] = 0e0;
	v[10640] = 0e0;
	v[10641] = 0e0;
	v[10642] = 0e0;
	v[10643] = 0e0;
	v[10644] = 0e0;
	v[9287] = 0e0;
	v[9288] = 0e0;
	v[9289] = 0e0;
	v[9290] = 0e0;
	v[9291] = 0e0;
	v[9292] = 0e0;
	v[9293] = 0e0;
	v[9294] = 0e0;
	v[9295] = 0e0;
	v[9296] = v[463];
	v[9297] = v[466];
	v[9298] = v[468];
	v[9299] = 0e0;
	v[9300] = 0e0;
	v[9301] = 0e0;
	v[9302] = 0e0;
	v[9303] = 0e0;
	v[9304] = 0e0;
	v[467] = v[48] / 2e0;
	v[9323] = 0e0;
	v[9324] = 0e0;
	v[9325] = 0e0;
	v[9326] = 0e0;
	v[9327] = 0e0;
	v[9328] = 0e0;
	v[9329] = 0e0;
	v[9330] = 0e0;
	v[9331] = 0e0;
	v[9332] = v[467];
	v[9333] = 0e0;
	v[9334] = v[465];
	v[9335] = 0e0;
	v[9336] = 0e0;
	v[9337] = 0e0;
	v[9338] = 0e0;
	v[9339] = 0e0;
	v[9340] = 0e0;
	v[9341] = 0e0;
	v[9342] = 0e0;
	v[9343] = 0e0;
	v[9344] = 0e0;
	v[9345] = 0e0;
	v[9346] = 0e0;
	v[9347] = 0e0;
	v[9348] = 0e0;
	v[9349] = 0e0;
	v[9350] = 0e0;
	v[9351] = v[467];
	v[9352] = v[464];
	v[9353] = 0e0;
	v[9354] = 0e0;
	v[9355] = 0e0;
	v[9356] = 0e0;
	v[9357] = 0e0;
	v[9358] = 0e0;
	v[191] = v[467] * v[47];
	v[555] = v[191] + v[46];
	v[6632] = -2e0*v[555];
	v[547] = v[191] - v[46];
	v[6634] = -2e0*v[547];
	v[189] = v[46] * v[467];
	v[551] = v[189] - v[47];
	v[6633] = -2e0*v[551];
	v[534] = v[189] + v[47];
	v[6636] = -2e0*v[534];
	v[180] = (v[48] * v[48]);
	v[1591] = 4e0 + v[179] + v[180] + v[186];
	v[7189] = 24e0 / Power(v[1591], 4);
	v[6429] = 8e0 / Power(v[1591], 3);
	v[2613] = -(v[463] * v[6429]);
	v[2611] = v[466] * v[6429];
	v[2608] = -(v[468] * v[6429]);
	v[1355] = 1e0 / (v[1591] * v[1591]);
	v[6430] = 4e0*v[1355];
	v[542] = -v[180] - v[186];
	v[6606] = v[542] / 2e0;
	v[523] = -v[179] - v[180];
	v[6607] = v[523] / 2e0;
	v[522] = v[468] * v[6430];
	v[6443] = -0.5e0*v[522];
	v[564] = v[560] * v[6443];
	v[521] = -(v[466] * v[6430]);
	v[6440] = v[521] / 2e0;
	v[544] = v[542] * v[6440];
	v[520] = v[463] * v[6430];
	v[6442] = -0.5e0*v[520];
	v[524] = v[523] * v[6442];
	v[301] = v[40] * v[7] + duiA[3] * v[8] + dduiA[3] * v[9];
	v[307] = v[41] * v[7] + duiA[4] * v[8] + dduiA[4] * v[9];
	v[309] = v[42] * v[7] + duiA[5] * v[8] + dduiA[5] * v[9];
	v[327] = v[46] * v[7] + duiA[9] * v[8] + dduiA[9] * v[9];
	v[333] = v[47] * v[7] + duiA[10] * v[8] + dduiA[10] * v[9];
	v[335] = v[48] * v[7] + duiA[11] * v[8] + dduiA[11] * v[9];
	v[73] = (*radA);
	v[74] = cpointA[0];
	v[75] = cpointA[1];
	v[76] = xAAi[0];
	v[77] = xAAi[1];
	v[78] = xAAi[2];
	v[79] = xBAi[0];
	v[80] = xBAi[1];
	v[81] = xBAi[2];
	v[82] = QAAi[0][0];
	v[83] = QAAi[0][1];
	v[85] = QAAi[1][0];
	v[86] = QAAi[1][1];
	v[88] = QAAi[2][0];
	v[89] = QAAi[2][1];
	v[91] = QBAi[0][0];
	v[92] = QBAi[0][1];
	v[94] = QBAi[1][0];
	v[95] = QBAi[1][1];
	v[97] = QBAi[2][0];
	v[98] = QBAi[2][1];
	v[100] = cAp[0];
	v[101] = cAp[1];
	v[1610] = v[73] * cos(v[101]);
	v[1024] = v[73] * sin(v[101]);
	v[102] = cAi[0];
	v[103] = cAi[1];
	v[104] = dB[0];
	v[105] = dB[1];
	v[106] = dB[2];
	v[107] = dB[3];
	v[471] = v[107] / 2e0;
	v[469] = 2e0*v[107];
	v[6838] = -0.5e0*v[469];
	v[246] = (v[107] * v[107]);
	v[108] = dB[4];
	v[472] = 2e0*v[108];
	v[10061] = 0e0;
	v[10062] = 0e0;
	v[10063] = 0e0;
	v[10064] = 0e0;
	v[10065] = 0e0;
	v[10066] = 0e0;
	v[10067] = 0e0;
	v[10068] = 0e0;
	v[10069] = 0e0;
	v[10070] = 0e0;
	v[10071] = 0e0;
	v[10072] = 0e0;
	v[10073] = 0e0;
	v[10074] = 0e0;
	v[10075] = 0e0;
	v[10076] = v[469];
	v[10077] = v[472];
	v[10078] = 0e0;
	v[6839] = -0.5e0*v[472];
	v[10663] = 0e0;
	v[10664] = 0e0;
	v[10665] = 0e0;
	v[10666] = 0e0;
	v[10667] = 0e0;
	v[10668] = 0e0;
	v[10669] = 0e0;
	v[10670] = 0e0;
	v[10671] = 0e0;
	v[10672] = 0e0;
	v[10673] = 0e0;
	v[10674] = 0e0;
	v[10675] = 0e0;
	v[10676] = 0e0;
	v[10677] = 0e0;
	v[10678] = -v[469];
	v[10679] = -v[472];
	v[10680] = 0e0;
	v[470] = v[108] / 2e0;
	v[9377] = 0e0;
	v[9378] = 0e0;
	v[9379] = 0e0;
	v[9380] = 0e0;
	v[9381] = 0e0;
	v[9382] = 0e0;
	v[9383] = 0e0;
	v[9384] = 0e0;
	v[9385] = 0e0;
	v[9386] = 0e0;
	v[9387] = 0e0;
	v[9388] = 0e0;
	v[9389] = 0e0;
	v[9390] = 0e0;
	v[9391] = 0e0;
	v[9392] = v[470];
	v[9393] = v[471];
	v[9394] = 0e0;
	v[244] = v[107] * v[470];
	v[239] = (v[108] * v[108]);
	v[731] = -v[239] - v[246];
	v[6575] = v[731] / 2e0;
	v[109] = dB[5];
	v[709] = v[109] + v[244];
	v[6629] = -2e0*v[709];
	v[699] = -v[109] + v[244];
	v[6627] = -2e0*v[699];
	v[474] = 2e0*v[109];
	v[9827] = 0e0;
	v[9828] = 0e0;
	v[9829] = 0e0;
	v[9830] = 0e0;
	v[9831] = 0e0;
	v[9832] = 0e0;
	v[9833] = 0e0;
	v[9834] = 0e0;
	v[9835] = 0e0;
	v[9836] = 0e0;
	v[9837] = 0e0;
	v[9838] = 0e0;
	v[9839] = 0e0;
	v[9840] = 0e0;
	v[9841] = 0e0;
	v[9842] = v[469];
	v[9843] = 0e0;
	v[9844] = v[474];
	v[9593] = 0e0;
	v[9594] = 0e0;
	v[9595] = 0e0;
	v[9596] = 0e0;
	v[9597] = 0e0;
	v[9598] = 0e0;
	v[9599] = 0e0;
	v[9600] = 0e0;
	v[9601] = 0e0;
	v[9602] = 0e0;
	v[9603] = 0e0;
	v[9604] = 0e0;
	v[9605] = 0e0;
	v[9606] = 0e0;
	v[9607] = 0e0;
	v[9608] = 0e0;
	v[9609] = v[472];
	v[9610] = v[474];
	v[6837] = -0.5e0*v[474];
	v[10735] = 0e0;
	v[10736] = 0e0;
	v[10737] = 0e0;
	v[10738] = 0e0;
	v[10739] = 0e0;
	v[10740] = 0e0;
	v[10741] = 0e0;
	v[10742] = 0e0;
	v[10743] = 0e0;
	v[10744] = 0e0;
	v[10745] = 0e0;
	v[10746] = 0e0;
	v[10747] = 0e0;
	v[10748] = 0e0;
	v[10749] = 0e0;
	v[10750] = 0e0;
	v[10751] = -v[472];
	v[10752] = -v[474];
	v[9359] = 0e0;
	v[9360] = 0e0;
	v[9361] = 0e0;
	v[9362] = 0e0;
	v[9363] = 0e0;
	v[9364] = 0e0;
	v[9365] = 0e0;
	v[9366] = 0e0;
	v[9367] = 0e0;
	v[9368] = 0e0;
	v[9369] = 0e0;
	v[9370] = 0e0;
	v[9371] = 0e0;
	v[9372] = 0e0;
	v[9373] = 0e0;
	v[9374] = v[469];
	v[9375] = v[472];
	v[9376] = v[474];
	v[473] = v[109] / 2e0;
	v[9395] = 0e0;
	v[9396] = 0e0;
	v[9397] = 0e0;
	v[9398] = 0e0;
	v[9399] = 0e0;
	v[9400] = 0e0;
	v[9401] = 0e0;
	v[9402] = 0e0;
	v[9403] = 0e0;
	v[9404] = 0e0;
	v[9405] = 0e0;
	v[9406] = 0e0;
	v[9407] = 0e0;
	v[9408] = 0e0;
	v[9409] = 0e0;
	v[9410] = v[473];
	v[9411] = 0e0;
	v[9412] = v[471];
	v[9413] = 0e0;
	v[9414] = 0e0;
	v[9415] = 0e0;
	v[9416] = 0e0;
	v[9417] = 0e0;
	v[9418] = 0e0;
	v[9419] = 0e0;
	v[9420] = 0e0;
	v[9421] = 0e0;
	v[9422] = 0e0;
	v[9423] = 0e0;
	v[9424] = 0e0;
	v[9425] = 0e0;
	v[9426] = 0e0;
	v[9427] = 0e0;
	v[9428] = 0e0;
	v[9429] = v[473];
	v[9430] = v[470];
	v[251] = v[108] * v[473];
	v[726] = v[107] + v[251];
	v[6630] = -2e0*v[726];
	v[718] = -v[107] + v[251];
	v[6628] = -2e0*v[718];
	v[249] = v[107] * v[473];
	v[722] = -v[108] + v[249];
	v[6631] = -2e0*v[722];
	v[705] = v[108] + v[249];
	v[6626] = -2e0*v[705];
	v[240] = (v[109] * v[109]);
	v[1598] = 4e0 + v[239] + v[240] + v[246];
	v[7173] = 24e0 / Power(v[1598], 4);
	v[6431] = 8e0 / Power(v[1598], 3);
	v[2639] = -(v[469] * v[6431]);
	v[2637] = v[472] * v[6431];
	v[2634] = -(v[474] * v[6431]);
	v[1359] = 1e0 / (v[1598] * v[1598]);
	v[6432] = 4e0*v[1359];
	v[713] = -v[240] - v[246];
	v[6576] = v[713] / 2e0;
	v[694] = -v[239] - v[240];
	v[6577] = v[694] / 2e0;
	v[693] = v[474] * v[6432];
	v[6450] = -0.5e0*v[693];
	v[735] = v[6450] * v[731];
	v[692] = -(v[472] * v[6432]);
	v[6447] = v[692] / 2e0;
	v[715] = v[6447] * v[713];
	v[691] = v[469] * v[6432];
	v[6449] = -0.5e0*v[691];
	v[695] = v[6449] * v[694];
	v[353] = v[107] * v[7] + duiB[3] * v[8] + dduiB[3] * v[9];
	v[359] = v[108] * v[7] + duiB[4] * v[8] + dduiB[4] * v[9];
	v[361] = v[109] * v[7] + duiB[5] * v[8] + dduiB[5] * v[9];
	v[122] = (*radB);
	v[123] = cpointB[0];
	v[124] = cpointB[1];
	v[125] = xABi[0];
	v[126] = xABi[1];
	v[127] = xABi[2];
	v[128] = QABi[0][0];
	v[129] = QABi[0][1];
	v[130] = QABi[0][2];
	v[131] = QABi[1][0];
	v[132] = QABi[1][1];
	v[133] = QABi[1][2];
	v[134] = QABi[2][0];
	v[135] = QABi[2][1];
	v[136] = QABi[2][2];
	v[137] = cBp[0];
	v[1611] = sin(v[137]);
	v[1025] = cos(v[137]);
	v[138] = cBp[1];
	v[1612] = v[122] * cos(v[138]);
	v[1026] = v[122] * sin(v[138]);
	v[451] = v[1026] * v[1611];
	v[449] = -(v[1025] * v[1026]);
	v[275] = v[123] + v[1612];
	v[1613] = v[1025] * v[275];
	v[1027] = -(v[1611] * v[275]);
	v[139] = cBi[0];
	v[140] = cBi[1];
	v[271] = v[123] + v[122] * cos(v[140]);
	v[159] = 4e0 / v[1584];
	v[6646] = 2e0*v[159];
	v[6434] = -0.5e0*v[159];
	v[517] = v[460] * v[6434];
	v[518] = v[517] + v[515] * v[6433];
	v[514] = v[457] * v[6434];
	v[516] = v[514] + v[515] * v[6435];
	v[511] = v[159] - v[475] * v[510];
	v[508] = -v[159] + v[476] * v[506];
	v[503] = -v[159] - v[475] * v[502];
	v[500] = v[462] * v[6434];
	v[501] = v[500] + v[497] * v[6436];
	v[498] = v[514] + v[497] * v[6435];
	v[496] = v[159] - v[477] * v[493];
	v[491] = v[159] + v[476] * v[489];
	v[488] = -(v[159] * v[461]);
	v[6653] = 2e0*v[488];
	v[512] = -v[488] + v[476] * v[510];
	v[587] = v[508] * v[83] + v[512] * v[86] + v[518] * v[89];
	v[584] = v[508] * v[82] + v[512] * v[85] + v[518] * v[88];
	v[507] = -v[488] - v[475] * v[506];
	v[586] = v[507] * v[83] + v[511] * v[86] + v[516] * v[89];
	v[583] = v[507] * v[82] + v[511] * v[85] + v[516] * v[88];
	v[504] = -v[488] + v[476] * v[502];
	v[490] = -v[488] - v[475] * v[489];
	v[487] = -v[159] - v[477] * v[483];
	v[485] = -(v[159] * v[459]);
	v[6654] = 2e0*v[485];
	v[509] = -v[485] - v[477] * v[506];
	v[495] = -v[485] + v[476] * v[493];
	v[578] = v[495] * v[83] + v[499] * v[86] + v[504] * v[89];
	v[575] = v[495] * v[82] + v[499] * v[85] + v[504] * v[88];
	v[492] = -v[485] - v[477] * v[489];
	v[486] = v[476] * v[483] - v[485];
	v[482] = v[159] * v[458];
	v[6655] = 2e0*v[482];
	v[513] = v[482] - v[477] * v[510];
	v[588] = v[509] * v[83] + v[513] * v[86] + v[519] * v[89];
	v[585] = v[509] * v[82] + v[513] * v[85] + v[519] * v[88];
	v[505] = v[482] - v[477] * v[502];
	v[579] = v[496] * v[83] + v[501] * v[86] + v[505] * v[89];
	v[576] = v[496] * v[82] + v[501] * v[85] + v[505] * v[88];
	v[494] = v[482] - v[475] * v[493];
	v[577] = v[494] * v[83] + v[498] * v[86] + v[503] * v[89];
	v[574] = v[494] * v[82] + v[498] * v[85] + v[503] * v[88];
	v[484] = v[482] - v[475] * v[483];
	v[568] = v[479] * v[83] + v[484] * v[86] + v[490] * v[89];
	v[565] = v[479] * v[82] + v[484] * v[85] + v[490] * v[88];
	v[481] = v[500] + v[478] * v[6436];
	v[570] = v[481] * v[83] + v[487] * v[86] + v[492] * v[89];
	v[567] = v[481] * v[82] + v[487] * v[85] + v[492] * v[88];
	v[480] = v[517] + v[478] * v[6433];
	v[569] = v[480] * v[83] + v[486] * v[86] + v[491] * v[89];
	v[566] = v[480] * v[82] + v[486] * v[85] + v[491] * v[88];
	v[298] = (v[159] * v[159]);
	v[6439] = v[309] / v[298];
	v[6438] = v[307] / v[298];
	v[6437] = v[301] / v[298];
	v[7200] = 2e0 / Power(v[298], 3);
	v[162] = 1e0 - v[478] * v[6434];
	v[7124] = v[162] / v[298];
	v[2224] = v[162] * v[6437];
	v[163] = v[159] * v[483];
	v[6620] = v[163] * v[307];
	v[2226] = v[163] * v[6438];
	v[3179] = v[2224] + v[2226];
	v[164] = v[159] * v[489];
	v[7217] = v[164] / v[298];
	v[6619] = v[164] * v[309];
	v[2227] = v[164] * v[6439];
	v[3184] = v[2224] + v[2227];
	v[3173] = v[2226] + v[3184];
	v[166] = v[159] * v[493];
	v[2231] = v[166] * v[6437];
	v[168] = 1e0 - v[497] * v[6434];
	v[7123] = v[168] / v[298];
	v[2220] = v[168] * v[6438];
	v[3175] = v[2220] + v[2231];
	v[169] = v[159] * v[502];
	v[7218] = v[169] / v[298];
	v[2221] = v[169] * v[6439];
	v[3185] = v[2220] + v[2221];
	v[3178] = v[2231] + v[3185];
	v[171] = v[159] * v[506];
	v[2234] = v[171] * v[6437];
	v[173] = v[159] * v[510];
	v[7213] = v[173] / v[298];
	v[2216] = v[173] * v[6438];
	v[174] = 1e0 - v[515] * v[6434];
	v[7127] = v[174] / v[298];
	v[2217] = v[174] * v[6439];
	v[7125] = v[2217] * v[298];
	v[3183] = v[2216] + v[2217] + v[2234];
	v[3180] = -v[2234] + v[3183];
	v[3174] = -v[2216] + v[3183];
	v[316] = -(v[482] * v[488]);
	v[6461] = v[316] - v[475];
	v[314] = v[485] * v[488];
	v[6459] = v[314] - v[476];
	v[305] = -(v[482] * v[485]);
	v[6456] = v[305] - v[477];
	v[178] = 4e0 / v[1591];
	v[6647] = 2e0*v[178];
	v[6441] = -0.5e0*v[178];
	v[562] = v[466] * v[6441];
	v[563] = v[562] + v[560] * v[6440];
	v[559] = v[463] * v[6441];
	v[561] = v[559] + v[560] * v[6442];
	v[556] = v[178] - v[520] * v[555];
	v[553] = -v[178] + v[521] * v[551];
	v[548] = -v[178] - v[520] * v[547];
	v[545] = v[468] * v[6441];
	v[546] = v[545] + v[542] * v[6443];
	v[543] = v[559] + v[542] * v[6442];
	v[541] = v[178] - v[522] * v[538];
	v[536] = v[178] + v[521] * v[534];
	v[533] = -(v[178] * v[467]);
	v[6657] = 2e0*v[533];
	v[557] = -v[533] + v[521] * v[555];
	v[614] = v[553] * v[92] + v[557] * v[95] + v[563] * v[98];
	v[611] = v[553] * v[91] + v[557] * v[94] + v[563] * v[97];
	v[552] = -v[533] - v[520] * v[551];
	v[613] = v[552] * v[92] + v[556] * v[95] + v[561] * v[98];
	v[610] = v[552] * v[91] + v[556] * v[94] + v[561] * v[97];
	v[549] = -v[533] + v[521] * v[547];
	v[535] = -v[533] - v[520] * v[534];
	v[532] = -v[178] - v[522] * v[528];
	v[530] = -(v[178] * v[465]);
	v[6658] = 2e0*v[530];
	v[554] = -v[530] - v[522] * v[551];
	v[540] = -v[530] + v[521] * v[538];
	v[605] = v[540] * v[92] + v[544] * v[95] + v[549] * v[98];
	v[602] = v[540] * v[91] + v[544] * v[94] + v[549] * v[97];
	v[537] = -v[530] - v[522] * v[534];
	v[531] = v[521] * v[528] - v[530];
	v[527] = v[178] * v[464];
	v[6659] = 2e0*v[527];
	v[558] = v[527] - v[522] * v[555];
	v[615] = v[554] * v[92] + v[558] * v[95] + v[564] * v[98];
	v[612] = v[554] * v[91] + v[558] * v[94] + v[564] * v[97];
	v[550] = v[527] - v[522] * v[547];
	v[606] = v[541] * v[92] + v[546] * v[95] + v[550] * v[98];
	v[603] = v[541] * v[91] + v[546] * v[94] + v[550] * v[97];
	v[539] = v[527] - v[520] * v[538];
	v[604] = v[539] * v[92] + v[543] * v[95] + v[548] * v[98];
	v[601] = v[539] * v[91] + v[543] * v[94] + v[548] * v[97];
	v[529] = v[527] - v[520] * v[528];
	v[595] = v[524] * v[92] + v[529] * v[95] + v[535] * v[98];
	v[592] = v[524] * v[91] + v[529] * v[94] + v[535] * v[97];
	v[526] = v[545] + v[523] * v[6443];
	v[597] = v[526] * v[92] + v[532] * v[95] + v[537] * v[98];
	v[594] = v[526] * v[91] + v[532] * v[94] + v[537] * v[97];
	v[525] = v[562] + v[523] * v[6440];
	v[596] = v[525] * v[92] + v[531] * v[95] + v[536] * v[98];
	v[593] = v[525] * v[91] + v[531] * v[94] + v[536] * v[97];
	v[324] = (v[178] * v[178]);
	v[6446] = v[335] / v[324];
	v[6445] = v[333] / v[324];
	v[6444] = v[327] / v[324];
	v[7190] = 2e0 / Power(v[324], 3);
	v[181] = 1e0 - v[523] * v[6441];
	v[7115] = v[181] / v[324];
	v[2200] = v[181] * v[6444];
	v[182] = v[178] * v[528];
	v[6602] = v[182] * v[333];
	v[2202] = v[182] * v[6445];
	v[3194] = v[2200] + v[2202];
	v[183] = v[178] * v[534];
	v[7229] = v[183] / v[324];
	v[6601] = v[183] * v[335];
	v[2203] = v[183] * v[6446];
	v[3199] = v[2200] + v[2203];
	v[3188] = v[2202] + v[3199];
	v[185] = v[178] * v[538];
	v[2207] = v[185] * v[6444];
	v[187] = 1e0 - v[542] * v[6441];
	v[7114] = v[187] / v[324];
	v[2196] = v[187] * v[6445];
	v[3190] = v[2196] + v[2207];
	v[188] = v[178] * v[547];
	v[7230] = v[188] / v[324];
	v[2197] = v[188] * v[6446];
	v[3200] = v[2196] + v[2197];
	v[3193] = v[2207] + v[3200];
	v[190] = v[178] * v[551];
	v[2210] = v[190] * v[6444];
	v[192] = v[178] * v[555];
	v[7225] = v[192] / v[324];
	v[2192] = v[192] * v[6445];
	v[193] = 1e0 - v[560] * v[6441];
	v[7118] = v[193] / v[324];
	v[2193] = v[193] * v[6446];
	v[7116] = v[2193] * v[324];
	v[3198] = v[2192] + v[2193] + v[2210];
	v[3195] = -v[2210] + v[3198];
	v[3189] = -v[2192] + v[3198];
	v[342] = -(v[527] * v[533]);
	v[6469] = v[342] - v[520];
	v[340] = v[530] * v[533];
	v[6467] = v[340] - v[521];
	v[331] = -(v[527] * v[530]);
	v[6464] = v[331] - v[522];
	v[197] = v[162] * v[82] + v[163] * v[85] + v[164] * v[88];
	v[198] = v[162] * v[83] + v[163] * v[86] + v[164] * v[89];
	v[437] = -(v[1024] * v[197]) + v[1610] * v[198];
	v[200] = v[166] * v[82] + v[168] * v[85] + v[169] * v[88];
	v[201] = v[166] * v[83] + v[168] * v[86] + v[169] * v[89];
	v[435] = -(v[1024] * v[200]) + v[1610] * v[201];
	v[203] = v[171] * v[82] + v[173] * v[85] + v[174] * v[88];
	v[204] = v[171] * v[83] + v[173] * v[86] + v[174] * v[89];
	v[433] = -(v[1024] * v[203]) + v[1610] * v[204];
	v[206] = v[181] * v[91] + v[182] * v[94] + v[183] * v[97];
	v[207] = v[181] * v[92] + v[182] * v[95] + v[183] * v[98];
	v[436] = -(v[1024] * v[206]) + v[1610] * v[207];
	v[209] = v[185] * v[91] + v[187] * v[94] + v[188] * v[97];
	v[210] = v[185] * v[92] + v[187] * v[95] + v[188] * v[98];
	v[434] = -(v[1024] * v[209]) + v[1610] * v[210];
	v[212] = v[190] * v[91] + v[192] * v[94] + v[193] * v[97];
	v[213] = v[190] * v[92] + v[192] * v[95] + v[193] * v[98];
	v[432] = -(v[1024] * v[212]) + v[1610] * v[213];
	v[215] = v[37] + v[76];
	v[216] = v[38] + v[77];
	v[217] = v[39] + v[78];
	v[218] = v[43] + v[79];
	v[219] = v[44] + v[80];
	v[220] = v[45] + v[81];
	v[221] = (1e0 - v[102]) / 2e0;
	v[222] = (1e0 + v[102]) / 2e0;
	v[223] = (1e0 - v[100]) / 2e0;
	v[633] = v[223] * (-(v[1024] * v[567]) + v[1610] * v[570]);
	v[632] = v[223] * (-(v[1024] * v[566]) + v[1610] * v[569]);
	v[631] = v[223] * (-(v[1024] * v[565]) + v[1610] * v[568]);
	v[627] = v[223] * (-(v[1024] * v[576]) + v[1610] * v[579]);
	v[626] = v[223] * (-(v[1024] * v[575]) + v[1610] * v[578]);
	v[625] = v[223] * (-(v[1024] * v[574]) + v[1610] * v[577]);
	v[621] = v[223] * (-(v[1024] * v[585]) + v[1610] * v[588]);
	v[620] = v[223] * (-(v[1024] * v[584]) + v[1610] * v[587]);
	v[619] = v[223] * (-(v[1024] * v[583]) + v[1610] * v[586]);
	v[224] = (1e0 + v[100]) / 2e0;
	v[12453] = 0e0;
	v[12454] = 0e0;
	v[12455] = -v[223];
	v[12456] = 0e0;
	v[12457] = 0e0;
	v[12458] = 0e0;
	v[12459] = 0e0;
	v[12460] = 0e0;
	v[12461] = -v[224];
	v[12462] = 0e0;
	v[12463] = 0e0;
	v[12464] = 0e0;
	v[12465] = 0e0;
	v[12466] = 0e0;
	v[12467] = 0e0;
	v[12468] = 0e0;
	v[12469] = 0e0;
	v[12470] = 0e0;
	v[12435] = 0e0;
	v[12436] = -v[223];
	v[12437] = 0e0;
	v[12438] = 0e0;
	v[12439] = 0e0;
	v[12440] = 0e0;
	v[12441] = 0e0;
	v[12442] = -v[224];
	v[12443] = 0e0;
	v[12444] = 0e0;
	v[12445] = 0e0;
	v[12446] = 0e0;
	v[12447] = 0e0;
	v[12448] = 0e0;
	v[12449] = 0e0;
	v[12450] = 0e0;
	v[12451] = 0e0;
	v[12452] = 0e0;
	v[12417] = -v[223];
	v[12418] = 0e0;
	v[12419] = 0e0;
	v[12420] = 0e0;
	v[12421] = 0e0;
	v[12422] = 0e0;
	v[12423] = -v[224];
	v[12424] = 0e0;
	v[12425] = 0e0;
	v[12426] = 0e0;
	v[12427] = 0e0;
	v[12428] = 0e0;
	v[12429] = 0e0;
	v[12430] = 0e0;
	v[12431] = 0e0;
	v[12432] = 0e0;
	v[12433] = 0e0;
	v[12434] = 0e0;
	v[1404] = v[204] * v[223] + v[213] * v[224];
	v[1403] = v[203] * v[223] + v[212] * v[224];
	v[1399] = v[201] * v[223] + v[210] * v[224];
	v[1398] = v[200] * v[223] + v[209] * v[224];
	v[1394] = v[198] * v[223] + v[207] * v[224];
	v[1393] = v[197] * v[223] + v[206] * v[224];
	v[636] = v[224] * (-(v[1024] * v[594]) + v[1610] * v[597]);
	v[635] = v[224] * (-(v[1024] * v[593]) + v[1610] * v[596]);
	v[634] = v[224] * (-(v[1024] * v[592]) + v[1610] * v[595]);
	v[630] = v[224] * (-(v[1024] * v[603]) + v[1610] * v[606]);
	v[629] = v[224] * (-(v[1024] * v[602]) + v[1610] * v[605]);
	v[628] = v[224] * (-(v[1024] * v[601]) + v[1610] * v[604]);
	v[624] = v[224] * (-(v[1024] * v[612]) + v[1610] * v[615]);
	v[623] = v[224] * (-(v[1024] * v[611]) + v[1610] * v[614]);
	v[622] = v[224] * (-(v[1024] * v[610]) + v[1610] * v[613]);
	v[440] = v[224] * v[432] + v[223] * v[433];
	v[816] = v[224] * v[440];
	v[810] = v[223] * v[440];
	v[439] = v[224] * v[434] + v[223] * v[435];
	v[815] = v[224] * v[439];
	v[809] = v[223] * v[439];
	v[438] = v[224] * v[436] + v[223] * v[437];
	v[814] = v[224] * v[438];
	v[808] = v[223] * v[438];
	v[225] = v[74] + v[73] * cos(v[103]);
	v[226] = v[75] + v[73] * sin(v[103]);
	v[227] = v[1610] + v[74];
	v[228] = v[1024] + v[75];
	v[1471] = v[227] * v[88] + v[228] * v[89];
	v[1470] = v[227] * v[85] + v[228] * v[86];
	v[1469] = v[227] * v[82] + v[228] * v[83];
	v[1467] = v[227] * v[97] + v[228] * v[98];
	v[1466] = v[227] * v[94] + v[228] * v[95];
	v[1465] = v[227] * v[91] + v[228] * v[92];
	v[666] = v[227] * v[567] + v[228] * v[570];
	v[675] = v[223] * v[666];
	v[665] = v[227] * v[566] + v[228] * v[569];
	v[674] = v[223] * v[665];
	v[664] = v[227] * v[565] + v[228] * v[568];
	v[673] = v[223] * v[664];
	v[663] = v[227] * v[594] + v[228] * v[597];
	v[678] = v[224] * v[663];
	v[662] = v[227] * v[593] + v[228] * v[596];
	v[677] = v[224] * v[662];
	v[661] = v[227] * v[592] + v[228] * v[595];
	v[676] = v[224] * v[661];
	v[654] = v[227] * v[576] + v[228] * v[579];
	v[681] = v[223] * v[654];
	v[653] = v[227] * v[575] + v[228] * v[578];
	v[680] = v[223] * v[653];
	v[652] = v[227] * v[574] + v[228] * v[577];
	v[679] = v[223] * v[652];
	v[651] = v[227] * v[603] + v[228] * v[606];
	v[684] = v[224] * v[651];
	v[650] = v[227] * v[602] + v[228] * v[605];
	v[683] = v[224] * v[650];
	v[649] = v[227] * v[601] + v[228] * v[604];
	v[682] = v[224] * v[649];
	v[642] = v[227] * v[585] + v[228] * v[588];
	v[687] = v[223] * v[642];
	v[641] = v[227] * v[584] + v[228] * v[587];
	v[686] = v[223] * v[641];
	v[640] = v[227] * v[583] + v[228] * v[586];
	v[685] = v[223] * v[640];
	v[639] = v[227] * v[612] + v[228] * v[615];
	v[690] = v[224] * v[639];
	v[638] = v[227] * v[611] + v[228] * v[614];
	v[689] = v[224] * v[638];
	v[637] = v[227] * v[610] + v[228] * v[613];
	v[688] = v[224] * v[637];
	v[427] = v[220] + v[212] * v[227] + v[213] * v[228];
	v[426] = v[217] + v[203] * v[227] + v[204] * v[228];
	v[428] = (-v[426] + v[427]) / 2e0;
	v[6480] = 2e0*v[428];
	v[424] = v[219] + v[209] * v[227] + v[210] * v[228];
	v[423] = v[216] + v[200] * v[227] + v[201] * v[228];
	v[425] = (-v[423] + v[424]) / 2e0;
	v[6479] = 2e0*v[425];
	v[421] = v[218] + v[206] * v[227] + v[207] * v[228];
	v[420] = v[215] + v[197] * v[227] + v[198] * v[228];
	v[422] = (-v[420] + v[421]) / 2e0;
	v[6478] = 2e0*v[422];
	v[238] = 4e0 / v[1598];
	v[6648] = 2e0*v[238];
	v[6448] = -0.5e0*v[238];
	v[733] = v[472] * v[6448];
	v[734] = v[6447] * v[731] + v[733];
	v[730] = v[469] * v[6448];
	v[732] = v[730] + v[6449] * v[731];
	v[727] = v[238] - v[691] * v[726];
	v[724] = -v[238] + v[692] * v[722];
	v[719] = -v[238] - v[691] * v[718];
	v[716] = v[474] * v[6448];
	v[717] = v[6450] * v[713] + v[716];
	v[714] = v[6449] * v[713] + v[730];
	v[712] = v[238] - v[693] * v[709];
	v[707] = v[238] + v[692] * v[705];
	v[704] = -(v[238] * v[473]);
	v[6667] = 2e0*v[704];
	v[728] = -v[704] + v[692] * v[726];
	v[773] = v[130] * v[724] + v[133] * v[728] + v[136] * v[734];
	v[770] = v[129] * v[724] + v[132] * v[728] + v[135] * v[734];
	v[767] = v[128] * v[724] + v[131] * v[728] + v[134] * v[734];
	v[723] = -v[704] - v[691] * v[722];
	v[772] = v[130] * v[723] + v[133] * v[727] + v[136] * v[732];
	v[769] = v[129] * v[723] + v[132] * v[727] + v[135] * v[732];
	v[766] = v[128] * v[723] + v[131] * v[727] + v[134] * v[732];
	v[720] = -v[704] + v[692] * v[718];
	v[706] = -v[704] - v[691] * v[705];
	v[703] = -v[238] - v[693] * v[699];
	v[701] = -(v[238] * v[471]);
	v[6668] = 2e0*v[701];
	v[725] = -v[701] - v[693] * v[722];
	v[711] = -v[701] + v[692] * v[709];
	v[758] = v[130] * v[711] + v[133] * v[715] + v[136] * v[720];
	v[755] = v[129] * v[711] + v[132] * v[715] + v[135] * v[720];
	v[752] = v[128] * v[711] + v[131] * v[715] + v[134] * v[720];
	v[708] = -v[701] - v[693] * v[705];
	v[702] = v[692] * v[699] - v[701];
	v[698] = v[238] * v[470];
	v[6669] = 2e0*v[698];
	v[729] = v[698] - v[693] * v[726];
	v[774] = v[130] * v[725] + v[133] * v[729] + v[136] * v[735];
	v[771] = v[129] * v[725] + v[132] * v[729] + v[135] * v[735];
	v[768] = v[128] * v[725] + v[131] * v[729] + v[134] * v[735];
	v[721] = v[698] - v[693] * v[718];
	v[759] = v[130] * v[712] + v[133] * v[717] + v[136] * v[721];
	v[756] = v[129] * v[712] + v[132] * v[717] + v[135] * v[721];
	v[753] = v[128] * v[712] + v[131] * v[717] + v[134] * v[721];
	v[710] = v[698] - v[691] * v[709];
	v[757] = v[130] * v[710] + v[133] * v[714] + v[136] * v[719];
	v[754] = v[129] * v[710] + v[132] * v[714] + v[135] * v[719];
	v[751] = v[128] * v[710] + v[131] * v[714] + v[134] * v[719];
	v[700] = v[698] - v[691] * v[699];
	v[742] = v[130] * v[695] + v[133] * v[700] + v[136] * v[706];
	v[739] = v[129] * v[695] + v[132] * v[700] + v[135] * v[706];
	v[736] = v[128] * v[695] + v[131] * v[700] + v[134] * v[706];
	v[697] = v[6450] * v[694] + v[716];
	v[744] = v[130] * v[697] + v[133] * v[703] + v[136] * v[708];
	v[741] = v[129] * v[697] + v[132] * v[703] + v[135] * v[708];
	v[738] = v[128] * v[697] + v[131] * v[703] + v[134] * v[708];
	v[696] = v[6447] * v[694] + v[733];
	v[743] = v[130] * v[696] + v[133] * v[702] + v[136] * v[707];
	v[740] = v[129] * v[696] + v[132] * v[702] + v[135] * v[707];
	v[737] = v[128] * v[696] + v[131] * v[702] + v[134] * v[707];
	v[350] = (v[238] * v[238]);
	v[6453] = v[361] / v[350];
	v[6452] = v[359] / v[350];
	v[6451] = v[353] / v[350];
	v[7174] = 2e0 / Power(v[350], 3);
	v[241] = 1e0 - v[6448] * v[694];
	v[7106] = v[241] / v[350];
	v[2176] = v[241] * v[6451];
	v[242] = v[238] * v[699];
	v[6572] = v[242] * v[359];
	v[2178] = v[242] * v[6452];
	v[3209] = v[2176] + v[2178];
	v[243] = v[238] * v[705];
	v[7241] = v[243] / v[350];
	v[6571] = v[243] * v[361];
	v[2179] = v[243] * v[6453];
	v[3214] = v[2176] + v[2179];
	v[3203] = v[2178] + v[3214];
	v[245] = v[238] * v[709];
	v[2183] = v[245] * v[6451];
	v[247] = 1e0 - v[6448] * v[713];
	v[7105] = v[247] / v[350];
	v[2172] = v[247] * v[6452];
	v[3205] = v[2172] + v[2183];
	v[248] = v[238] * v[718];
	v[7242] = v[248] / v[350];
	v[2173] = v[248] * v[6453];
	v[3215] = v[2172] + v[2173];
	v[3208] = v[2183] + v[3215];
	v[250] = v[238] * v[722];
	v[2186] = v[250] * v[6451];
	v[252] = v[238] * v[726];
	v[7237] = v[252] / v[350];
	v[2168] = v[252] * v[6452];
	v[253] = 1e0 - v[6448] * v[731];
	v[7109] = v[253] / v[350];
	v[2169] = v[253] * v[6453];
	v[7107] = v[2169] * v[350];
	v[3213] = v[2168] + v[2169] + v[2186];
	v[3210] = -v[2186] + v[3213];
	v[3204] = -v[2168] + v[3213];
	v[368] = -(v[698] * v[704]);
	v[6477] = v[368] - v[691];
	v[366] = v[701] * v[704];
	v[6475] = v[366] - v[692];
	v[357] = -(v[698] * v[701]);
	v[6472] = v[357] - v[693];
	v[257] = v[128] * v[241] + v[131] * v[242] + v[134] * v[243];
	v[258] = v[129] * v[241] + v[132] * v[242] + v[135] * v[243];
	v[259] = v[130] * v[241] + v[133] * v[242] + v[136] * v[243];
	v[453] = v[1612] * v[258] + v[257] * v[449] + v[259] * v[451];
	v[844] = -(v[224] * v[453]);
	v[838] = -(v[223] * v[453]);
	v[443] = v[1027] * v[257] - v[1613] * v[259];
	v[829] = -(v[224] * v[443]);
	v[823] = -(v[223] * v[443]);
	v[260] = v[128] * v[245] + v[131] * v[247] + v[134] * v[248];
	v[261] = v[129] * v[245] + v[132] * v[247] + v[135] * v[248];
	v[262] = v[130] * v[245] + v[133] * v[247] + v[136] * v[248];
	v[454] = v[1612] * v[261] + v[260] * v[449] + v[262] * v[451];
	v[845] = -(v[224] * v[454]);
	v[839] = -(v[223] * v[454]);
	v[444] = v[1027] * v[260] - v[1613] * v[262];
	v[830] = -(v[224] * v[444]);
	v[824] = -(v[223] * v[444]);
	v[263] = v[128] * v[250] + v[131] * v[252] + v[134] * v[253];
	v[264] = v[129] * v[250] + v[132] * v[252] + v[135] * v[253];
	v[265] = v[130] * v[250] + v[133] * v[252] + v[136] * v[253];
	v[455] = v[1612] * v[264] + v[263] * v[449] + v[265] * v[451];
	v[849] = -(v[453] * v[678]) - v[454] * v[684] - v[455] * v[690];
	v[848] = -(v[453] * v[677]) - v[454] * v[683] - v[455] * v[689];
	v[847] = -(v[453] * v[676]) - v[454] * v[682] - v[455] * v[688];
	v[846] = -(v[224] * v[455]);
	v[843] = -(v[453] * v[675]) - v[454] * v[681] - v[455] * v[687];
	v[842] = -(v[453] * v[674]) - v[454] * v[680] - v[455] * v[686];
	v[841] = -(v[453] * v[673]) - v[454] * v[679] - v[455] * v[685];
	v[840] = -(v[223] * v[455]);
	v[445] = v[1027] * v[263] - v[1613] * v[265];
	v[834] = -(v[443] * v[678]) - v[444] * v[684] - v[445] * v[690];
	v[833] = -(v[443] * v[677]) - v[444] * v[683] - v[445] * v[689];
	v[832] = -(v[443] * v[676]) - v[444] * v[682] - v[445] * v[688];
	v[831] = -(v[224] * v[445]);
	v[828] = -(v[443] * v[675]) - v[444] * v[681] - v[445] * v[687];
	v[827] = -(v[443] * v[674]) - v[444] * v[680] - v[445] * v[686];
	v[826] = -(v[443] * v[673]) - v[444] * v[679] - v[445] * v[685];
	v[825] = -(v[223] * v[445]);
	v[266] = v[104] + v[125];
	v[267] = v[105] + v[126];
	v[268] = v[106] + v[127];
	v[269] = v[271] * cos(v[139]);
	v[270] = v[124] + v[122] * sin(v[140]);
	v[272] = -(v[271] * sin(v[139]));
	v[274] = v[1026] + v[124];
	v[789] = v[1613] * v[768] + v[274] * v[771] + v[1027] * v[774];
	v[788] = v[1613] * v[767] + v[274] * v[770] + v[1027] * v[773];
	v[787] = v[1613] * v[766] + v[274] * v[769] + v[1027] * v[772];
	v[786] = v[1613] * v[753] + v[274] * v[756] + v[1027] * v[759];
	v[785] = v[1613] * v[752] + v[274] * v[755] + v[1027] * v[758];
	v[784] = v[1613] * v[751] + v[274] * v[754] + v[1027] * v[757];
	v[783] = v[1613] * v[738] + v[274] * v[741] + v[1027] * v[744];
	v[822] = -(v[438] * v[783]) - v[439] * v[786] - v[440] * v[789];
	v[807] = -(v[422] * v[783]) - v[425] * v[786] - v[428] * v[789];
	v[782] = v[1613] * v[737] + v[274] * v[740] + v[1027] * v[743];
	v[821] = -(v[438] * v[782]) - v[439] * v[785] - v[440] * v[788];
	v[806] = -(v[422] * v[782]) - v[425] * v[785] - v[428] * v[788];
	v[781] = v[1613] * v[736] + v[274] * v[739] + v[1027] * v[742];
	v[820] = -(v[438] * v[781]) - v[439] * v[784] - v[440] * v[787];
	v[805] = -(v[422] * v[781]) - v[425] * v[784] - v[428] * v[787];
	v[286] = v[37] * v[7] + duiA[0] * v[8] + dduiA[0] * v[9];
	v[3458] = v[223] * v[286];
	v[287] = v[38] * v[7] + duiA[1] * v[8] + dduiA[1] * v[9];
	v[3480] = v[223] * v[287];
	v[288] = v[39] * v[7] + duiA[2] * v[8] + dduiA[2] * v[9];
	v[3030] = v[223] * v[288];
	v[289] = v[43] * v[7] + duiA[6] * v[8] + dduiA[6] * v[9];
	v[7153] = v[286] - v[289];
	v[3459] = v[224] * v[289];
	v[5483] = v[3458] + v[3459];
	v[290] = v[44] * v[7] + duiA[7] * v[8] + dduiA[7] * v[9];
	v[7154] = v[287] - v[290];
	v[3481] = v[224] * v[290];
	v[5481] = v[3480] + v[3481];
	v[291] = v[45] * v[7] + duiA[8] * v[8] + dduiA[8] * v[9];
	v[7383] = v[288] - v[291];
	v[3031] = v[224] * v[291];
	v[5461] = v[3030] + v[3031];
	v[292] = v[104] * v[7] + duiB[0] * v[8] + dduiB[0] * v[9];
	v[4632] = -v[292] + v[5483];
	v[293] = v[105] * v[7] + duiB[1] * v[8] + dduiB[1] * v[9];
	v[4600] = -v[293] + v[5481];
	v[294] = v[106] * v[7] + duiB[2] * v[8] + dduiB[2] * v[9];
	v[4566] = -v[294] + v[5461];
	v[296] = v[314] + v[476];
	v[6460] = v[296] / v[298];
	v[6454] = v[174] * v[296];
	v[297] = v[305] + v[477];
	v[6457] = v[297] / v[298];
	v[6455] = v[168] * v[297];
	v[299] = v[298] + (v[485] * v[485]);
	v[6401] = -(v[309] * v[6454]) - v[307] * v[6455] - v[299] * (v[301] + v[6619] + v[6620]);
	v[4497] = v[299] / v[298];
	v[2506] = v[6459] / v[298];
	v[2504] = v[6456] / v[298];
	v[12237] = 0e0;
	v[12238] = 0e0;
	v[12239] = 0e0;
	v[12240] = 0e0;
	v[12241] = v[2504];
	v[12242] = v[2506];
	v[12243] = 0e0;
	v[12244] = 0e0;
	v[12245] = 0e0;
	v[12246] = 0e0;
	v[12247] = 0e0;
	v[12248] = 0e0;
	v[12249] = 0e0;
	v[12250] = 0e0;
	v[12251] = 0e0;
	v[12252] = 0e0;
	v[12253] = 0e0;
	v[12254] = 0e0;
	v[11113] = 0e0;
	v[11114] = 0e0;
	v[11115] = 0e0;
	v[11116] = 0e0;
	v[11117] = v[2504] * v[7];
	v[11118] = v[2506] * v[7];
	v[11119] = 0e0;
	v[11120] = 0e0;
	v[11121] = 0e0;
	v[11122] = 0e0;
	v[11123] = 0e0;
	v[11124] = 0e0;
	v[11125] = 0e0;
	v[11126] = 0e0;
	v[11127] = 0e0;
	v[11128] = 0e0;
	v[11129] = 0e0;
	v[11130] = 0e0;
	v[300] = v[2504] * v[307] + v[2506] * v[309] + v[301] * v[4497];
	v[302] = v[298] + (v[482] * v[482]);
	v[6608] = v[166] * v[302] + v[162] * v[6456];
	v[308] = v[316] + v[475];
	v[6458] = v[169] * v[302] + v[174] * v[308];
	v[6400] = -(v[302] * v[307]) - v[309] * v[6458] - v[301] * v[6608];
	v[4500] = v[302] / v[298];
	v[2505] = v[6461] / v[298];
	v[12255] = 0e0;
	v[12256] = 0e0;
	v[12257] = 0e0;
	v[12258] = v[6457];
	v[12259] = 0e0;
	v[12260] = v[2505];
	v[12261] = 0e0;
	v[12262] = 0e0;
	v[12263] = 0e0;
	v[12264] = 0e0;
	v[12265] = 0e0;
	v[12266] = 0e0;
	v[12267] = 0e0;
	v[12268] = 0e0;
	v[12269] = 0e0;
	v[12270] = 0e0;
	v[12271] = 0e0;
	v[12272] = 0e0;
	v[11149] = 0e0;
	v[11150] = 0e0;
	v[11151] = 0e0;
	v[11152] = v[6457] * v[7];
	v[11153] = 0e0;
	v[11154] = v[2505] * v[7];
	v[11155] = 0e0;
	v[11156] = 0e0;
	v[11157] = 0e0;
	v[11158] = 0e0;
	v[11159] = 0e0;
	v[11160] = 0e0;
	v[11161] = 0e0;
	v[11162] = 0e0;
	v[11163] = 0e0;
	v[11164] = 0e0;
	v[11165] = 0e0;
	v[11166] = 0e0;
	v[310] = v[2505] * v[309] + v[307] * v[4500] + v[301] * v[6457];
	v[311] = v[298] + (v[488] * v[488]);
	v[6609] = v[171] * v[311] + v[162] * v[6459];
	v[6610] = v[173] * v[311] + v[168] * v[6461];
	v[6399] = -(v[309] * v[311]) - v[301] * v[6609] - v[307] * v[6610];
	v[4502] = v[311] / v[298];
	v[2503] = v[308] / v[298];
	v[12273] = 0e0;
	v[12274] = 0e0;
	v[12275] = 0e0;
	v[12276] = v[6460];
	v[12277] = v[2503];
	v[12278] = 0e0;
	v[12279] = 0e0;
	v[12280] = 0e0;
	v[12281] = 0e0;
	v[12282] = 0e0;
	v[12283] = 0e0;
	v[12284] = 0e0;
	v[12285] = 0e0;
	v[12286] = 0e0;
	v[12287] = 0e0;
	v[12288] = 0e0;
	v[12289] = 0e0;
	v[12290] = 0e0;
	v[11185] = 0e0;
	v[11186] = 0e0;
	v[11187] = 0e0;
	v[11188] = v[6460] * v[7];
	v[11189] = v[2503] * v[7];
	v[11190] = 0e0;
	v[11191] = 0e0;
	v[11192] = 0e0;
	v[11193] = 0e0;
	v[11194] = 0e0;
	v[11195] = 0e0;
	v[11196] = 0e0;
	v[11197] = 0e0;
	v[11198] = 0e0;
	v[11199] = 0e0;
	v[11200] = 0e0;
	v[11201] = 0e0;
	v[11202] = 0e0;
	v[320] = v[2503] * v[307] + v[309] * v[4502] + v[301] * v[6460];
	v[322] = v[340] + v[521];
	v[6468] = v[322] / v[324];
	v[6462] = v[193] * v[322];
	v[323] = v[331] + v[522];
	v[6465] = v[323] / v[324];
	v[6463] = v[187] * v[323];
	v[325] = v[324] + (v[530] * v[530]);
	v[6378] = -(v[335] * v[6462]) - v[333] * v[6463] - v[325] * (v[327] + v[6601] + v[6602]);
	v[4504] = v[325] / v[324];
	v[2514] = v[6467] / v[324];
	v[2512] = v[6464] / v[324];
	v[12291] = 0e0;
	v[12292] = 0e0;
	v[12293] = 0e0;
	v[12294] = 0e0;
	v[12295] = 0e0;
	v[12296] = 0e0;
	v[12297] = 0e0;
	v[12298] = 0e0;
	v[12299] = 0e0;
	v[12300] = 0e0;
	v[12301] = v[2512];
	v[12302] = v[2514];
	v[12303] = 0e0;
	v[12304] = 0e0;
	v[12305] = 0e0;
	v[12306] = 0e0;
	v[12307] = 0e0;
	v[12308] = 0e0;
	v[11221] = 0e0;
	v[11222] = 0e0;
	v[11223] = 0e0;
	v[11224] = 0e0;
	v[11225] = 0e0;
	v[11226] = 0e0;
	v[11227] = 0e0;
	v[11228] = 0e0;
	v[11229] = 0e0;
	v[11230] = 0e0;
	v[11231] = v[2512] * v[7];
	v[11232] = v[2514] * v[7];
	v[11233] = 0e0;
	v[11234] = 0e0;
	v[11235] = 0e0;
	v[11236] = 0e0;
	v[11237] = 0e0;
	v[11238] = 0e0;
	v[326] = v[2512] * v[333] + v[2514] * v[335] + v[327] * v[4504];
	v[328] = v[324] + (v[527] * v[527]);
	v[6590] = v[185] * v[328] + v[181] * v[6464];
	v[334] = v[342] + v[520];
	v[6466] = v[188] * v[328] + v[193] * v[334];
	v[6377] = -(v[328] * v[333]) - v[335] * v[6466] - v[327] * v[6590];
	v[4507] = v[328] / v[324];
	v[2513] = v[6469] / v[324];
	v[12309] = 0e0;
	v[12310] = 0e0;
	v[12311] = 0e0;
	v[12312] = 0e0;
	v[12313] = 0e0;
	v[12314] = 0e0;
	v[12315] = 0e0;
	v[12316] = 0e0;
	v[12317] = 0e0;
	v[12318] = v[6465];
	v[12319] = 0e0;
	v[12320] = v[2513];
	v[12321] = 0e0;
	v[12322] = 0e0;
	v[12323] = 0e0;
	v[12324] = 0e0;
	v[12325] = 0e0;
	v[12326] = 0e0;
	v[11257] = 0e0;
	v[11258] = 0e0;
	v[11259] = 0e0;
	v[11260] = 0e0;
	v[11261] = 0e0;
	v[11262] = 0e0;
	v[11263] = 0e0;
	v[11264] = 0e0;
	v[11265] = 0e0;
	v[11266] = v[6465] * v[7];
	v[11267] = 0e0;
	v[11268] = v[2513] * v[7];
	v[11269] = 0e0;
	v[11270] = 0e0;
	v[11271] = 0e0;
	v[11272] = 0e0;
	v[11273] = 0e0;
	v[11274] = 0e0;
	v[336] = v[2513] * v[335] + v[333] * v[4507] + v[327] * v[6465];
	v[337] = v[324] + (v[533] * v[533]);
	v[6591] = v[190] * v[337] + v[181] * v[6467];
	v[6592] = v[192] * v[337] + v[187] * v[6469];
	v[6376] = -(v[335] * v[337]) - v[327] * v[6591] - v[333] * v[6592];
	v[4509] = v[337] / v[324];
	v[2511] = v[334] / v[324];
	v[12327] = 0e0;
	v[12328] = 0e0;
	v[12329] = 0e0;
	v[12330] = 0e0;
	v[12331] = 0e0;
	v[12332] = 0e0;
	v[12333] = 0e0;
	v[12334] = 0e0;
	v[12335] = 0e0;
	v[12336] = v[6468];
	v[12337] = v[2511];
	v[12338] = 0e0;
	v[12339] = 0e0;
	v[12340] = 0e0;
	v[12341] = 0e0;
	v[12342] = 0e0;
	v[12343] = 0e0;
	v[12344] = 0e0;
	v[11293] = 0e0;
	v[11294] = 0e0;
	v[11295] = 0e0;
	v[11296] = 0e0;
	v[11297] = 0e0;
	v[11298] = 0e0;
	v[11299] = 0e0;
	v[11300] = 0e0;
	v[11301] = 0e0;
	v[11302] = v[6468] * v[7];
	v[11303] = v[2511] * v[7];
	v[11304] = 0e0;
	v[11305] = 0e0;
	v[11306] = 0e0;
	v[11307] = 0e0;
	v[11308] = 0e0;
	v[11309] = 0e0;
	v[11310] = 0e0;
	v[346] = v[2511] * v[333] + v[335] * v[4509] + v[327] * v[6468];
	v[348] = v[366] + v[692];
	v[6476] = v[348] / v[350];
	v[6470] = v[253] * v[348];
	v[349] = v[357] + v[693];
	v[6473] = v[349] / v[350];
	v[6471] = v[247] * v[349];
	v[351] = v[350] + (v[701] * v[701]);
	v[6352] = -(v[361] * v[6470]) - v[359] * v[6471] - v[351] * (v[353] + v[6571] + v[6572]);
	v[4511] = v[351] / v[350];
	v[2522] = v[6475] / v[350];
	v[2520] = v[6472] / v[350];
	v[12345] = 0e0;
	v[12346] = 0e0;
	v[12347] = 0e0;
	v[12348] = 0e0;
	v[12349] = 0e0;
	v[12350] = 0e0;
	v[12351] = 0e0;
	v[12352] = 0e0;
	v[12353] = 0e0;
	v[12354] = 0e0;
	v[12355] = 0e0;
	v[12356] = 0e0;
	v[12357] = 0e0;
	v[12358] = 0e0;
	v[12359] = 0e0;
	v[12360] = 0e0;
	v[12361] = v[2520];
	v[12362] = v[2522];
	v[11329] = 0e0;
	v[11330] = 0e0;
	v[11331] = 0e0;
	v[11332] = 0e0;
	v[11333] = 0e0;
	v[11334] = 0e0;
	v[11335] = 0e0;
	v[11336] = 0e0;
	v[11337] = 0e0;
	v[11338] = 0e0;
	v[11339] = 0e0;
	v[11340] = 0e0;
	v[11341] = 0e0;
	v[11342] = 0e0;
	v[11343] = 0e0;
	v[11344] = 0e0;
	v[11345] = v[2520] * v[7];
	v[11346] = v[2522] * v[7];
	v[352] = v[2520] * v[359] + v[2522] * v[361] + v[353] * v[4511];
	v[354] = v[350] + (v[698] * v[698]);
	v[6560] = v[245] * v[354] + v[241] * v[6472];
	v[360] = v[368] + v[691];
	v[6474] = v[248] * v[354] + v[253] * v[360];
	v[6351] = -(v[354] * v[359]) - v[361] * v[6474] - v[353] * v[6560];
	v[4514] = v[354] / v[350];
	v[2521] = v[6477] / v[350];
	v[12363] = 0e0;
	v[12364] = 0e0;
	v[12365] = 0e0;
	v[12366] = 0e0;
	v[12367] = 0e0;
	v[12368] = 0e0;
	v[12369] = 0e0;
	v[12370] = 0e0;
	v[12371] = 0e0;
	v[12372] = 0e0;
	v[12373] = 0e0;
	v[12374] = 0e0;
	v[12375] = 0e0;
	v[12376] = 0e0;
	v[12377] = 0e0;
	v[12378] = v[6473];
	v[12379] = 0e0;
	v[12380] = v[2521];
	v[11365] = 0e0;
	v[11366] = 0e0;
	v[11367] = 0e0;
	v[11368] = 0e0;
	v[11369] = 0e0;
	v[11370] = 0e0;
	v[11371] = 0e0;
	v[11372] = 0e0;
	v[11373] = 0e0;
	v[11374] = 0e0;
	v[11375] = 0e0;
	v[11376] = 0e0;
	v[11377] = 0e0;
	v[11378] = 0e0;
	v[11379] = 0e0;
	v[11380] = v[6473] * v[7];
	v[11381] = 0e0;
	v[11382] = v[2521] * v[7];
	v[362] = v[2521] * v[361] + v[359] * v[4514] + v[353] * v[6473];
	v[363] = v[350] + (v[704] * v[704]);
	v[6561] = v[250] * v[363] + v[241] * v[6475];
	v[6562] = v[252] * v[363] + v[247] * v[6477];
	v[6350] = -(v[361] * v[363]) - v[353] * v[6561] - v[359] * v[6562];
	v[4516] = v[363] / v[350];
	v[2519] = v[360] / v[350];
	v[12381] = 0e0;
	v[12382] = 0e0;
	v[12383] = 0e0;
	v[12384] = 0e0;
	v[12385] = 0e0;
	v[12386] = 0e0;
	v[12387] = 0e0;
	v[12388] = 0e0;
	v[12389] = 0e0;
	v[12390] = 0e0;
	v[12391] = 0e0;
	v[12392] = 0e0;
	v[12393] = 0e0;
	v[12394] = 0e0;
	v[12395] = 0e0;
	v[12396] = v[6476];
	v[12397] = v[2519];
	v[12398] = 0e0;
	v[11401] = 0e0;
	v[11402] = 0e0;
	v[11403] = 0e0;
	v[11404] = 0e0;
	v[11405] = 0e0;
	v[11406] = 0e0;
	v[11407] = 0e0;
	v[11408] = 0e0;
	v[11409] = 0e0;
	v[11410] = 0e0;
	v[11411] = 0e0;
	v[11412] = 0e0;
	v[11413] = 0e0;
	v[11414] = 0e0;
	v[11415] = 0e0;
	v[11416] = v[6476] * v[7];
	v[11417] = v[2519] * v[7];
	v[11418] = 0e0;
	v[372] = v[2519] * v[359] + v[361] * v[4516] + v[353] * v[6476];
	v[3426] = -v[292] + v[300] * v[673] + v[310] * v[674] + v[320] * v[675] + v[326] * v[676] + v[336] * v[677] + v[346] * v[678]
		- v[352] * v[781] - v[362] * v[782] - v[372] * v[783];
	v[3301] = v[3426] + v[5483];
	v[3241] = -v[293] + v[300] * v[679] + v[310] * v[680] + v[320] * v[681] + v[326] * v[682] + v[336] * v[683] + v[346] * v[684]
		- v[352] * v[784] - v[362] * v[785] - v[372] * v[786];
	v[3401] = v[3241] + v[5481];
	v[3219] = -v[294] + v[300] * v[685] + v[310] * v[686] + v[320] * v[687] + v[326] * v[688] + v[336] * v[689] + v[346] * v[690]
		- v[352] * v[787] - v[362] * v[788] - v[372] * v[789];
	v[3506] = v[3219] + v[5461];
	v[373] = -(v[1613] * v[257]) - v[1027] * v[259] - v[266] - v[258] * v[274] + v[223] * v[420] + v[224] * v[421];
	v[796] = -0.5e0*v[373];
	v[797] = v[224] * v[422] - v[796];
	v[790] = v[223] * v[422] + v[796];
	v[374] = -(v[1613] * v[260]) - v[1027] * v[262] - v[267] - v[261] * v[274] + v[223] * v[423] + v[224] * v[424];
	v[798] = -0.5e0*v[374];
	v[799] = v[224] * v[425] - v[798];
	v[791] = v[223] * v[425] + v[798];
	v[375] = -(v[1613] * v[263]) - v[1027] * v[265] - v[268] - v[264] * v[274] + v[223] * v[426] + v[224] * v[427];
	v[981] = -(v[1027] * v[136]) - v[134] * v[1613] - v[135] * v[274];
	v[982] = -(v[1027] * v[133]) - v[131] * v[1613] - v[132] * v[274];
	v[983] = -(v[1027] * v[130]) - v[128] * v[1613] - v[129] * v[274];
	v[984] = v[1467] * v[224];
	v[985] = v[1466] * v[224];
	v[986] = v[1465] * v[224];
	v[987] = v[1471] * v[223];
	v[988] = v[1470] * v[223];
	v[989] = v[1469] * v[223];
	v[6913] = v[6432] / 2e0;
	v[6651] = -2e0*v[6432];
	v[6899] = v[6430] / 2e0;
	v[6650] = -2e0*v[6430];
	v[6885] = v[6428] / 2e0;
	v[6649] = -2e0*v[6428];
	v[852] = -(v[373] * (v[449] * v[738] + v[1612] * v[741] + v[451] * v[744])) - v[374] * (v[449] * v[753] + v[1612] * v[756]
		+ v[451] * v[759]) - v[375] * (v[449] * v[768] + v[1612] * v[771] + v[451] * v[774]) + v[453] * v[783] + v[454] * v[786]
		+ v[455] * v[789];
	v[851] = -(v[373] * (v[449] * v[737] + v[1612] * v[740] + v[451] * v[743])) - v[374] * (v[449] * v[752] + v[1612] * v[755]
		+ v[451] * v[758]) - v[375] * (v[449] * v[767] + v[1612] * v[770] + v[451] * v[773]) + v[453] * v[782] + v[454] * v[785]
		+ v[455] * v[788];
	v[850] = -(v[373] * (v[449] * v[736] + v[1612] * v[739] + v[451] * v[742])) - v[374] * (v[449] * v[751] + v[1612] * v[754]
		+ v[451] * v[757]) - v[375] * (v[449] * v[766] + v[1612] * v[769] + v[451] * v[772]) + v[453] * v[781] + v[454] * v[784]
		+ v[455] * v[787];
	v[837] = -(v[373] * (v[1027] * v[738] - v[1613] * v[744])) - v[374] * (v[1027] * v[753] - v[1613] * v[759]) - v[375] *
		(v[1027] * v[768] - v[1613] * v[774]) + v[443] * v[783] + v[444] * v[786] + v[445] * v[789];
	v[836] = -(v[373] * (v[1027] * v[737] - v[1613] * v[743])) - v[374] * (v[1027] * v[752] - v[1613] * v[758]) - v[375] *
		(v[1027] * v[767] - v[1613] * v[773]) + v[443] * v[782] + v[444] * v[785] + v[445] * v[788];
	v[835] = -(v[373] * (v[1027] * v[736] - v[1613] * v[742])) - v[374] * (v[1027] * v[751] - v[1613] * v[757]) - v[375] *
		(v[1027] * v[766] - v[1613] * v[772]) + v[443] * v[781] + v[444] * v[784] + v[445] * v[787];
	v[819] = v[375] * v[624] + v[374] * v[630] + v[373] * v[636] + v[438] * v[678] + v[439] * v[684] + v[440] * v[690];
	v[818] = v[375] * v[623] + v[374] * v[629] + v[373] * v[635] + v[438] * v[677] + v[439] * v[683] + v[440] * v[689];
	v[817] = v[375] * v[622] + v[374] * v[628] + v[373] * v[634] + v[438] * v[676] + v[439] * v[682] + v[440] * v[688];
	v[813] = v[375] * v[621] + v[374] * v[627] + v[373] * v[633] + v[438] * v[675] + v[439] * v[681] + v[440] * v[687];
	v[812] = v[375] * v[620] + v[374] * v[626] + v[373] * v[632] + v[438] * v[674] + v[439] * v[680] + v[440] * v[686];
	v[811] = v[375] * v[619] + v[374] * v[625] + v[373] * v[631] + v[438] * v[673] + v[439] * v[679] + v[440] * v[685];
	v[804] = (v[375] * v[639] + v[374] * v[651] + v[373] * v[663] + v[6478] * v[678] + v[6479] * v[684] + v[6480] * v[690])
		/ 2e0;
	v[803] = (v[375] * v[638] + v[374] * v[650] + v[373] * v[662] + v[6478] * v[677] + v[6479] * v[683] + v[6480] * v[689])
		/ 2e0;
	v[802] = (v[375] * v[637] + v[374] * v[649] + v[373] * v[661] + v[6478] * v[676] + v[6479] * v[682] + v[6480] * v[688])
		/ 2e0;
	v[800] = -0.5e0*v[375];
	v[801] = v[224] * v[428] - v[800];
	v[795] = v[422] * v[675] + v[425] * v[681] + v[428] * v[687] + v[666] * v[796] + v[654] * v[798] + v[642] * v[800];
	v[794] = v[422] * v[674] + v[425] * v[680] + v[428] * v[686] + v[665] * v[796] + v[653] * v[798] + v[641] * v[800];
	v[793] = v[422] * v[673] + v[425] * v[679] + v[428] * v[685] + v[664] * v[796] + v[652] * v[798] + v[640] * v[800];
	v[792] = v[223] * v[428] + v[800];
	v[401] = sqrt((v[373] * v[373]) + (v[374] * v[374]) + (v[375] * v[375]));
	v[1406] = 1e0 / (v[401] * v[401]);
	b1049 = v[401] < v[27];
	v[376] = -v[125] - v[128] * v[269] - v[129] * v[270] - v[130] * v[272] + v[221] * (v[76] + v[225] * v[82] + v[226] * v[83])
		+ v[222] * (v[79] + v[225] * v[91] + v[226] * v[92]);
	v[377] = -v[126] - v[131] * v[269] - v[132] * v[270] - v[133] * v[272] + v[221] * (v[77] + v[225] * v[85] + v[226] * v[86])
		+ v[222] * (v[80] + v[225] * v[94] + v[226] * v[95]);
	v[378] = -v[127] - v[134] * v[269] - v[135] * v[270] - v[136] * v[272] + v[221] * (v[78] + v[225] * v[88] + v[226] * v[89])
		+ v[222] * (v[81] + v[225] * v[97] + v[226] * v[98]);
	v[411] = sqrt((v[376] * v[376]) + (v[377] * v[377]) + (v[378] * v[378]));
	if (b6) {
		v[7477] = Power(v[401], v[3585]);
		v[7476] = Power(v[401], v[1045]);
		v[7341] = Power(v[401], v[1045]);
		v[7334] = Power(v[401], v[3599]);
		v[3780] = Power(v[401], v[3585]);
		v[7257] = Power(v[401], v[1045]);
		if (v[401] > 0.1e-7) { v01 = 1e0 / v[401]; v02 = (-(v01 / v[401])); v03 = (2e0*v01) / (v[401] * v[401]); }
		else {
			v01 = (12500000e0 / 3e0)*(24e0 - (-0.1e-7 + v[401])*(0.24e10 - 2e0*(-1e0 + 100000000e0*v[401])*
				(0.2399999997e10 - 0.1199999994e18*v[401] - 0.3e17*(v[401] * v[401]))));
			v02 = (-50000000e0 / 3e0)*(0.3599999994e10 - 0.4799999982e18*v[401] + 0.6e25*Power(v[401], 3)
				+ 0.1799999982e26*(v[401] * v[401]));
			v03 = 0.1e17*(799999997e0 - 0.599999994e17*v[401] - 0.3e17*(v[401] * v[401]));
		};
		v[386] = v03;
		v[387] = v02;
		v[388] = v01;
		v[389] = v[373] * v[388];
		v[390] = v[374] * v[388];
		v[391] = v[375] * v[388];
		if (v[411] > 0.1e-7) { v04 = 1e0 / v[411]; v05 = (-(v04 / v[411])); v06 = (2e0*v04) / (v[411] * v[411]); }
		else {
			v04 = (12500000e0 / 3e0)*(24e0 - (-0.1e-7 + v[411])*(0.24e10 - 2e0*(-1e0 + 100000000e0*v[411])*
				(0.2399999997e10 - 0.1199999994e18*v[411] - 0.3e17*(v[411] * v[411]))));
			v05 = (-50000000e0 / 3e0)*(0.3599999994e10 - 0.4799999982e18*v[411] + 0.6e25*Power(v[411], 3)
				+ 0.1799999982e26*(v[411] * v[411]));
			v06 = 0.1e17*(799999997e0 - 0.599999994e17*v[411] - 0.3e17*(v[411] * v[411]));
		};
		v[397] = v04;
		v[398] = v[376] * v[397];
		v[399] = v[377] * v[397];
		v[400] = v[378] * v[397];
	}
	else {
		if (v[401] > 0.1e-7) { v07 = 1e0 / v[401]; v08 = (-(v07 / v[401])); v09 = (2e0*v07) / (v[401] * v[401]); }
		else {
			v07 = (12500000e0 / 3e0)*(24e0 - (-0.1e-7 + v[401])*(0.24e10 - 2e0*(-1e0 + 100000000e0*v[401])*
				(0.2399999997e10 - 0.1199999994e18*v[401] - 0.3e17*(v[401] * v[401]))));
			v08 = (-50000000e0 / 3e0)*(0.3599999994e10 - 0.4799999982e18*v[401] + 0.6e25*Power(v[401], 3)
				+ 0.1799999982e26*(v[401] * v[401]));
			v09 = 0.1e17*(799999997e0 - 0.599999994e17*v[401] - 0.3e17*(v[401] * v[401]));
		};
		v[405] = v09;
		v[406] = v08;
		v[407] = v07;
		v[408] = v[373] * v[407];
		v[409] = v[374] * v[407];
		v[410] = v[375] * v[407];
		if (v[411] > 0.1e-7) { v010 = 1e0 / v[411]; v011 = (-(v010 / v[411])); v012 = (2e0*v010) / (v[411] * v[411]); }
		else {
			v010 = (12500000e0 / 3e0)*(24e0 - (-0.1e-7 + v[411])*(0.24e10 - 2e0*(-1e0 + 100000000e0*v[411])*
				(0.2399999997e10 - 0.1199999994e18*v[411] - 0.3e17*(v[411] * v[411]))));
			v011 = (-50000000e0 / 3e0)*(0.3599999994e10 - 0.4799999982e18*v[411] + 0.6e25*Power(v[411], 3)
				+ 0.1799999982e26*(v[411] * v[411]));
			v012 = 0.1e17*(799999997e0 - 0.599999994e17*v[411] - 0.3e17*(v[411] * v[411]));
		};
		v[417] = v010;
		v[398] = v[376] * v[417];
		v[399] = v[377] * v[417];
		v[400] = v[378] * v[417];
		b418 = v[398] * v[408] + v[399] * v[409] + v[400] * v[410] < 0e0;
		if (b418) {
			v[389] = -v[408];
			v[390] = -v[409];
			v[391] = -v[410];
		}
		else {
			v[389] = v[408];
			v[390] = v[409];
			v[391] = v[410];
		};
	};
	v[12399] = 0e0;
	v[12400] = 0e0;
	v[12401] = v[223] * v[391];
	v[12402] = 0e0;
	v[12403] = 0e0;
	v[12404] = 0e0;
	v[12405] = 0e0;
	v[12406] = 0e0;
	v[12407] = v[224] * v[391];
	v[12408] = 0e0;
	v[12409] = 0e0;
	v[12410] = 0e0;
	v[12411] = 0e0;
	v[12412] = 0e0;
	v[12413] = 0e0;
	v[12414] = 0e0;
	v[12415] = 0e0;
	v[12416] = 0e0;
	v[6731] = v[26] * v[391];
	v[6680] = 2e0*v[391];
	v[1965] = -(v[391] * v[789]);
	v[1963] = -(v[391] * v[788]);
	v[1961] = -(v[391] * v[787]);
	v[1959] = v[391] * v[690];
	v[1957] = v[391] * v[689];
	v[1955] = v[391] * v[688];
	v[1953] = v[391] * v[687];
	v[1951] = v[391] * v[686];
	v[1949] = v[391] * v[685];
	v[3826] = (-v[293] + v[3480] + v[3481])*v[391];
	v[3530] = v[3826] + (v[293] + v[3241])*v[391];
	v[3825] = (-v[292] + v[3458] + v[3459])*v[391];
	v[3552] = v[3825] + (v[292] + v[3426])*v[391];
	v[1022] = (v[391] * v[391]);
	v[1946] = (-v[294] + v[3030] + v[3031])*v[391];
	v[6732] = v[26] * v[390];
	v[6683] = 2e0*v[390];
	v[3425] = v[3426] * v[390];
	v[3337] = v[3219] * v[390];
	v[1964] = -(v[390] * v[786]);
	v[3283] = v[1964] + v[1965];
	v[1962] = -(v[390] * v[785]);
	v[3285] = v[1962] + v[1963];
	v[1960] = -(v[390] * v[784]);
	v[3287] = v[1960] + v[1961];
	v[1958] = v[390] * v[684];
	v[3289] = v[1958] + v[1959];
	v[1956] = v[390] * v[683];
	v[3291] = v[1956] + v[1957];
	v[1954] = v[390] * v[682];
	v[3293] = v[1954] + v[1955];
	v[1952] = v[390] * v[681];
	v[3295] = v[1952] + v[1953];
	v[1950] = v[390] * v[680];
	v[3297] = v[1950] + v[1951];
	v[1948] = v[390] * v[679];
	v[3299] = v[1948] + v[1949];
	v[1947] = -(v[293] * v[390]);
	v[4636] = v[1946] + v[1947];
	v[3262] = v[326] * v[3293] + v[320] * v[3295] + v[310] * v[3297] + v[300] * v[3299] + v[3291] * v[336] + v[3289] * v[346]
		+ v[3287] * v[352] + v[3285] * v[362] + v[3283] * v[372] + v[4636];
	v[1018] = (v[390] * v[390]);
	v[6818] = v[289] * v[389] + v[290] * v[390];
	v[6812] = v[286] * v[389] + v[287] * v[390];
	v[6733] = v[26] * v[389];
	v[6686] = 2e0*v[389];
	v[6546] = v[389] * v[391];
	v[6499] = v[389] * v[390];
	v[6497] = v[389] * v[673];
	v[6498] = v[1949] + v[6497];
	v[6495] = v[389] * v[674];
	v[6496] = v[1951] + v[6495];
	v[6493] = v[389] * v[675];
	v[6494] = v[1953] + v[6493];
	v[6491] = v[389] * v[676];
	v[6492] = v[1955] + v[6491];
	v[6489] = v[389] * v[677];
	v[6490] = v[1957] + v[6489];
	v[6487] = v[389] * v[678];
	v[6488] = v[1959] + v[6487];
	v[6485] = -(v[389] * v[781]);
	v[6486] = v[1961] + v[6485];
	v[6483] = -(v[389] * v[782]);
	v[6484] = v[1963] + v[6483];
	v[6481] = -(v[389] * v[783]);
	v[6482] = v[1965] + v[6481];
	v[3240] = v[3241] * v[389];
	v[3217] = v[3219] * v[389];
	v[1905] = v[390] * v[6498] + v[1018] * v[679];
	v[1903] = v[390] * v[6496] + v[1018] * v[680];
	v[1901] = v[390] * v[6494] + v[1018] * v[681];
	v[1899] = v[390] * v[6492] + v[1018] * v[682];
	v[1897] = v[390] * v[6490] + v[1018] * v[683];
	v[1895] = v[390] * v[6488] + v[1018] * v[684];
	v[1893] = v[390] * v[6486] - v[1018] * v[784];
	v[1891] = v[390] * v[6484] - v[1018] * v[785];
	v[1889] = v[390] * v[6482] - v[1018] * v[786];
	v[7449] = v[3283] + 2e0*v[6481];
	v[3482] = v[1964] + v[6481];
	v[7404] = 2e0*v[1965] + v[3482];
	v[7423] = 2e0*v[1964] + v[6482];
	v[7448] = v[3285] + 2e0*v[6483];
	v[3484] = v[1962] + v[6483];
	v[7403] = 2e0*v[1963] + v[3484];
	v[7422] = 2e0*v[1962] + v[6484];
	v[7447] = v[3287] + 2e0*v[6485];
	v[3486] = v[1960] + v[6485];
	v[7402] = 2e0*v[1961] + v[3486];
	v[7421] = 2e0*v[1960] + v[6486];
	v[7446] = v[3289] + 2e0*v[6487];
	v[3488] = v[1958] + v[6487];
	v[7401] = 2e0*v[1959] + v[3488];
	v[7420] = 2e0*v[1958] + v[6488];
	v[7445] = v[3291] + 2e0*v[6489];
	v[3490] = v[1956] + v[6489];
	v[7400] = 2e0*v[1957] + v[3490];
	v[7419] = 2e0*v[1956] + v[6490];
	v[7444] = v[3293] + 2e0*v[6491];
	v[3492] = v[1954] + v[6491];
	v[7399] = 2e0*v[1955] + v[3492];
	v[7418] = 2e0*v[1954] + v[6492];
	v[7443] = v[3295] + 2e0*v[6493];
	v[3494] = v[1952] + v[6493];
	v[7398] = 2e0*v[1953] + v[3494];
	v[7417] = 2e0*v[1952] + v[6494];
	v[7442] = v[3297] + 2e0*v[6495];
	v[3496] = v[1950] + v[6495];
	v[7397] = 2e0*v[1951] + v[3496];
	v[7416] = 2e0*v[1950] + v[6496];
	v[7441] = v[3299] + 2e0*v[6497];
	v[3498] = v[1948] + v[6497];
	v[7396] = 2e0*v[1949] + v[3498];
	v[7415] = 2e0*v[1948] + v[6498];
	v[1877] = -(v[292] * v[389]);
	v[4602] = v[1877] + v[1946];
	v[4572] = v[1877] + v[1947];
	v[3446] = v[346] * v[3488] + v[336] * v[3490] + v[326] * v[3492] + v[320] * v[3494] + v[310] * v[3496] + v[300] * v[3498]
		+ v[3486] * v[352] + v[3484] * v[362] + v[3482] * v[372] + v[4572] + v[223] * v[6812] + v[224] * v[6818];
	v[3362] = v[4602] + v[372] * v[6482] + v[362] * v[6484] + v[352] * v[6486] + v[346] * v[6488] + v[336] * v[6490]
		+ v[326] * v[6492] + v[320] * v[6494] + v[310] * v[6496] + v[300] * v[6498];
	v[1832] = v[3498] * v[391] + v[1022] * v[685];
	v[1830] = v[3496] * v[391] + v[1022] * v[686];
	v[1828] = v[3494] * v[391] + v[1022] * v[687];
	v[1826] = v[3492] * v[391] + v[1022] * v[688];
	v[1824] = v[3490] * v[391] + v[1022] * v[689];
	v[1822] = v[3488] * v[391] + v[1022] * v[690];
	v[1820] = v[3486] * v[391] - v[1022] * v[787];
	v[1818] = v[3484] * v[391] - v[1022] * v[788];
	v[1816] = v[3482] * v[391] - v[1022] * v[789];
	v[1017] = v[224] * v[6499];
	v[1016] = v[223] * v[6499];
	v[1014] = (v[389] * v[389]);
	v[1984] = v[3299] * v[389] + v[1014] * v[673];
	v[1982] = v[3297] * v[389] + v[1014] * v[674];
	v[1980] = v[3295] * v[389] + v[1014] * v[675];
	v[1978] = v[3293] * v[389] + v[1014] * v[676];
	v[1976] = v[3291] * v[389] + v[1014] * v[677];
	v[1974] = v[3289] * v[389] + v[1014] * v[678];
	v[1972] = v[3287] * v[389] - v[1014] * v[781];
	v[1970] = v[3285] * v[389] - v[1014] * v[782];
	v[1968] = v[3283] * v[389] - v[1014] * v[783];
	v[869] = -(v[10] * v[790]) - v[11] * v[808] - v[12] * v[823] - v[13] * v[838];
	v[870] = -(v[10] * v[791]) - v[11] * v[809] - v[12] * v[824] - v[13] * v[839];
	v[871] = -(v[10] * v[792]) - v[11] * v[810] - v[12] * v[825] - v[13] * v[840];
	v[872] = -(v[10] * v[793]) - v[11] * v[811] - v[12] * v[826] - v[13] * v[841];
	v[873] = -(v[10] * v[794]) - v[11] * v[812] - v[12] * v[827] - v[13] * v[842];
	v[874] = -(v[10] * v[795]) - v[11] * v[813] - v[12] * v[828] - v[13] * v[843];
	v[875] = -(v[10] * v[797]) - v[11] * v[814] - v[12] * v[829] - v[13] * v[844];
	v[876] = -(v[10] * v[799]) - v[11] * v[815] - v[12] * v[830] - v[13] * v[845];
	v[877] = -(v[10] * v[801]) - v[11] * v[816] - v[12] * v[831] - v[13] * v[846];
	v[878] = -(v[10] * v[802]) - v[11] * v[817] - v[12] * v[832] - v[13] * v[847];
	v[879] = -(v[10] * v[803]) - v[11] * v[818] - v[12] * v[833] - v[13] * v[848];
	v[880] = -(v[10] * v[804]) - v[11] * v[819] - v[12] * v[834] - v[13] * v[849];
	v[881] = v[10] * v[422] + v[11] * v[438] - v[12] * v[443] - v[13] * v[453];
	v[882] = v[10] * v[425] + v[11] * v[439] - v[12] * v[444] - v[13] * v[454];
	v[883] = v[10] * v[428] + v[11] * v[440] - v[12] * v[445] - v[13] * v[455];
	v[884] = -(v[10] * v[805]) - v[11] * v[820] - v[12] * v[835] - v[13] * v[850];
	v[885] = -(v[10] * v[806]) - v[11] * v[821] - v[12] * v[836] - v[13] * v[851];
	v[886] = -(v[10] * v[807]) - v[11] * v[822] - v[12] * v[837] - v[13] * v[852];
	v[6544] = (v[286] * v[869] + v[287] * v[870] + v[288] * v[871] + v[300] * v[872] + v[310] * v[873] + v[320] * v[874]
		+ v[289] * v[875] + v[290] * v[876] + v[291] * v[877] + v[326] * v[878] + v[336] * v[879] + v[346] * v[880] + v[292] * v[881]
		+ v[293] * v[882] + v[294] * v[883] + v[352] * v[884] + v[362] * v[885] + v[372] * v[886]) / 2e0;
	v[9157] = v[869];
	v[9158] = v[870];
	v[9159] = v[871];
	v[9160] = v[872];
	v[9161] = v[873];
	v[9162] = v[874];
	v[9163] = v[875];
	v[9164] = v[876];
	v[9165] = v[877];
	v[9166] = v[878];
	v[9167] = v[879];
	v[9168] = v[880];
	v[9169] = v[881];
	v[9170] = v[882];
	v[9171] = v[883];
	v[9172] = v[884];
	v[9173] = v[885];
	v[9174] = v[886];
	v[887] = -(v[14] * v[790]) - v[15] * v[808] - v[16] * v[823] - v[17] * v[838];
	v[888] = -(v[14] * v[791]) - v[15] * v[809] - v[16] * v[824] - v[17] * v[839];
	v[889] = -(v[14] * v[792]) - v[15] * v[810] - v[16] * v[825] - v[17] * v[840];
	v[890] = -(v[14] * v[793]) - v[15] * v[811] - v[16] * v[826] - v[17] * v[841];
	v[891] = -(v[14] * v[794]) - v[15] * v[812] - v[16] * v[827] - v[17] * v[842];
	v[892] = -(v[14] * v[795]) - v[15] * v[813] - v[16] * v[828] - v[17] * v[843];
	v[893] = -(v[14] * v[797]) - v[15] * v[814] - v[16] * v[829] - v[17] * v[844];
	v[894] = -(v[14] * v[799]) - v[15] * v[815] - v[16] * v[830] - v[17] * v[845];
	v[895] = -(v[14] * v[801]) - v[15] * v[816] - v[16] * v[831] - v[17] * v[846];
	v[896] = -(v[14] * v[802]) - v[15] * v[817] - v[16] * v[832] - v[17] * v[847];
	v[897] = -(v[14] * v[803]) - v[15] * v[818] - v[16] * v[833] - v[17] * v[848];
	v[898] = -(v[14] * v[804]) - v[15] * v[819] - v[16] * v[834] - v[17] * v[849];
	v[899] = v[14] * v[422] + v[15] * v[438] - v[16] * v[443] - v[17] * v[453];
	v[900] = v[14] * v[425] + v[15] * v[439] - v[16] * v[444] - v[17] * v[454];
	v[901] = v[14] * v[428] + v[15] * v[440] - v[16] * v[445] - v[17] * v[455];
	v[902] = -(v[14] * v[805]) - v[15] * v[820] - v[16] * v[835] - v[17] * v[850];
	v[903] = -(v[14] * v[806]) - v[15] * v[821] - v[16] * v[836] - v[17] * v[851];
	v[904] = -(v[14] * v[807]) - v[15] * v[822] - v[16] * v[837] - v[17] * v[852];
	v[1747] = -(v[286] * v[887]) - v[287] * v[888] - v[288] * v[889] - v[300] * v[890] - v[310] * v[891] - v[320] * v[892]
		- v[289] * v[893] - v[290] * v[894] - v[291] * v[895] - v[326] * v[896] - v[336] * v[897] - v[346] * v[898] - v[292] * v[899]
		- v[293] * v[900] - v[294] * v[901] - v[352] * v[902] - v[362] * v[903] - v[372] * v[904];
	v[9175] = v[887];
	v[9176] = v[888];
	v[9177] = v[889];
	v[9178] = v[890];
	v[9179] = v[891];
	v[9180] = v[892];
	v[9181] = v[893];
	v[9182] = v[894];
	v[9183] = v[895];
	v[9184] = v[896];
	v[9185] = v[897];
	v[9186] = v[898];
	v[9187] = v[899];
	v[9188] = v[900];
	v[9189] = v[901];
	v[9190] = v[902];
	v[9191] = v[903];
	v[9192] = v[904];
	v[905] = -(v[18] * v[790]) - v[19] * v[808] - v[20] * v[823] - v[21] * v[838];
	v[906] = -(v[18] * v[791]) - v[19] * v[809] - v[20] * v[824] - v[21] * v[839];
	v[907] = -(v[18] * v[792]) - v[19] * v[810] - v[20] * v[825] - v[21] * v[840];
	v[908] = -(v[18] * v[793]) - v[19] * v[811] - v[20] * v[826] - v[21] * v[841];
	v[909] = -(v[18] * v[794]) - v[19] * v[812] - v[20] * v[827] - v[21] * v[842];
	v[910] = -(v[18] * v[795]) - v[19] * v[813] - v[20] * v[828] - v[21] * v[843];
	v[911] = -(v[18] * v[797]) - v[19] * v[814] - v[20] * v[829] - v[21] * v[844];
	v[912] = -(v[18] * v[799]) - v[19] * v[815] - v[20] * v[830] - v[21] * v[845];
	v[913] = -(v[18] * v[801]) - v[19] * v[816] - v[20] * v[831] - v[21] * v[846];
	v[914] = -(v[18] * v[802]) - v[19] * v[817] - v[20] * v[832] - v[21] * v[847];
	v[915] = -(v[18] * v[803]) - v[19] * v[818] - v[20] * v[833] - v[21] * v[848];
	v[916] = -(v[18] * v[804]) - v[19] * v[819] - v[20] * v[834] - v[21] * v[849];
	v[917] = v[18] * v[422] + v[19] * v[438] - v[20] * v[443] - v[21] * v[453];
	v[918] = v[18] * v[425] + v[19] * v[439] - v[20] * v[444] - v[21] * v[454];
	v[919] = v[18] * v[428] + v[19] * v[440] - v[20] * v[445] - v[21] * v[455];
	v[920] = -(v[18] * v[805]) - v[19] * v[820] - v[20] * v[835] - v[21] * v[850];
	v[921] = -(v[18] * v[806]) - v[19] * v[821] - v[20] * v[836] - v[21] * v[851];
	v[922] = -(v[18] * v[807]) - v[19] * v[822] - v[20] * v[837] - v[21] * v[852];
	v[1741] = v[286] * v[905] + v[287] * v[906] + v[288] * v[907] + v[300] * v[908] + v[310] * v[909] + v[320] * v[910]
		+ v[289] * v[911] + v[290] * v[912] + v[291] * v[913] + v[326] * v[914] + v[336] * v[915] + v[346] * v[916] + v[292] * v[917]
		+ v[293] * v[918] + v[294] * v[919] + v[352] * v[920] + v[362] * v[921] + v[372] * v[922];
	v[9139] = v[905];
	v[9140] = v[906];
	v[9141] = v[907];
	v[9142] = v[908];
	v[9143] = v[909];
	v[9144] = v[910];
	v[9145] = v[911];
	v[9146] = v[912];
	v[9147] = v[913];
	v[9148] = v[914];
	v[9149] = v[915];
	v[9150] = v[916];
	v[9151] = v[917];
	v[9152] = v[918];
	v[9153] = v[919];
	v[9154] = v[920];
	v[9155] = v[921];
	v[9156] = v[922];
	v[923] = -(v[22] * v[790]) - v[23] * v[808] - v[24] * v[823] - v[25] * v[838];
	v[1616] = v[428] * v[869] + v[440] * v[887] - v[445] * v[905] - v[455] * v[923];
	v[1615] = v[425] * v[869] + v[439] * v[887] - v[444] * v[905] - v[454] * v[923];
	v[1614] = v[422] * v[869] + v[438] * v[887] - v[443] * v[905] - v[453] * v[923];
	v[924] = -(v[22] * v[791]) - v[23] * v[809] - v[24] * v[824] - v[25] * v[839];
	v[1620] = v[428] * v[870] + v[440] * v[888] - v[445] * v[906] - v[455] * v[924];
	v[1619] = v[425] * v[870] + v[439] * v[888] - v[444] * v[906] - v[454] * v[924];
	v[1618] = v[422] * v[870] + v[438] * v[888] - v[443] * v[906] - v[453] * v[924];
	v[925] = -(v[22] * v[792]) - v[23] * v[810] - v[24] * v[825] - v[25] * v[840];
	v[1624] = v[428] * v[871] + v[440] * v[889] - v[445] * v[907] - v[455] * v[925];
	v[1623] = v[425] * v[871] + v[439] * v[889] - v[444] * v[907] - v[454] * v[925];
	v[1622] = v[422] * v[871] + v[438] * v[889] - v[443] * v[907] - v[453] * v[925];
	v[926] = -(v[22] * v[793]) - v[23] * v[811] - v[24] * v[826] - v[25] * v[841];
	v[1628] = v[428] * v[872] + v[440] * v[890] - v[445] * v[908] - v[455] * v[926];
	v[1627] = v[425] * v[872] + v[439] * v[890] - v[444] * v[908] - v[454] * v[926];
	v[1626] = v[422] * v[872] + v[438] * v[890] - v[443] * v[908] - v[453] * v[926];
	v[927] = -(v[22] * v[794]) - v[23] * v[812] - v[24] * v[827] - v[25] * v[842];
	v[1632] = v[428] * v[873] + v[440] * v[891] - v[445] * v[909] - v[455] * v[927];
	v[1631] = v[425] * v[873] + v[439] * v[891] - v[444] * v[909] - v[454] * v[927];
	v[1630] = v[422] * v[873] + v[438] * v[891] - v[443] * v[909] - v[453] * v[927];
	v[928] = -(v[22] * v[795]) - v[23] * v[813] - v[24] * v[828] - v[25] * v[843];
	v[1636] = v[428] * v[874] + v[440] * v[892] - v[445] * v[910] - v[455] * v[928];
	v[1635] = v[425] * v[874] + v[439] * v[892] - v[444] * v[910] - v[454] * v[928];
	v[1634] = v[422] * v[874] + v[438] * v[892] - v[443] * v[910] - v[453] * v[928];
	v[929] = -(v[22] * v[797]) - v[23] * v[814] - v[24] * v[829] - v[25] * v[844];
	v[1640] = v[428] * v[875] + v[440] * v[893] - v[445] * v[911] - v[455] * v[929];
	v[1639] = v[425] * v[875] + v[439] * v[893] - v[444] * v[911] - v[454] * v[929];
	v[1638] = v[422] * v[875] + v[438] * v[893] - v[443] * v[911] - v[453] * v[929];
	v[930] = -(v[22] * v[799]) - v[23] * v[815] - v[24] * v[830] - v[25] * v[845];
	v[1644] = v[428] * v[876] + v[440] * v[894] - v[445] * v[912] - v[455] * v[930];
	v[1643] = v[425] * v[876] + v[439] * v[894] - v[444] * v[912] - v[454] * v[930];
	v[1642] = v[422] * v[876] + v[438] * v[894] - v[443] * v[912] - v[453] * v[930];
	v[931] = -(v[22] * v[801]) - v[23] * v[816] - v[24] * v[831] - v[25] * v[846];
	v[1648] = v[428] * v[877] + v[440] * v[895] - v[445] * v[913] - v[455] * v[931];
	v[1647] = v[425] * v[877] + v[439] * v[895] - v[444] * v[913] - v[454] * v[931];
	v[1646] = v[422] * v[877] + v[438] * v[895] - v[443] * v[913] - v[453] * v[931];
	v[932] = -(v[22] * v[802]) - v[23] * v[817] - v[24] * v[832] - v[25] * v[847];
	v[1652] = v[428] * v[878] + v[440] * v[896] - v[445] * v[914] - v[455] * v[932];
	v[1651] = v[425] * v[878] + v[439] * v[896] - v[444] * v[914] - v[454] * v[932];
	v[1650] = v[422] * v[878] + v[438] * v[896] - v[443] * v[914] - v[453] * v[932];
	v[933] = -(v[22] * v[803]) - v[23] * v[818] - v[24] * v[833] - v[25] * v[848];
	v[1656] = v[428] * v[879] + v[440] * v[897] - v[445] * v[915] - v[455] * v[933];
	v[1655] = v[425] * v[879] + v[439] * v[897] - v[444] * v[915] - v[454] * v[933];
	v[1654] = v[422] * v[879] + v[438] * v[897] - v[443] * v[915] - v[453] * v[933];
	v[934] = -(v[22] * v[804]) - v[23] * v[819] - v[24] * v[834] - v[25] * v[849];
	v[1660] = v[428] * v[880] + v[440] * v[898] - v[445] * v[916] - v[455] * v[934];
	v[1659] = v[425] * v[880] + v[439] * v[898] - v[444] * v[916] - v[454] * v[934];
	v[1658] = v[422] * v[880] + v[438] * v[898] - v[443] * v[916] - v[453] * v[934];
	v[935] = v[22] * v[422] + v[23] * v[438] - v[24] * v[443] - v[25] * v[453];
	v[1664] = v[428] * v[881] + v[440] * v[899] - v[445] * v[917] - v[455] * v[935];
	v[1663] = v[425] * v[881] + v[439] * v[899] - v[444] * v[917] - v[454] * v[935];
	v[1662] = v[422] * v[881] + v[438] * v[899] - v[443] * v[917] - v[453] * v[935];
	v[936] = v[22] * v[425] + v[23] * v[439] - v[24] * v[444] - v[25] * v[454];
	v[1668] = v[428] * v[882] + v[440] * v[900] - v[445] * v[918] - v[455] * v[936];
	v[1667] = v[425] * v[882] + v[439] * v[900] - v[444] * v[918] - v[454] * v[936];
	v[1666] = v[422] * v[882] + v[438] * v[900] - v[443] * v[918] - v[453] * v[936];
	v[937] = v[22] * v[428] + v[23] * v[440] - v[24] * v[445] - v[25] * v[455];
	v[1672] = v[428] * v[883] + v[440] * v[901] - v[445] * v[919] - v[455] * v[937];
	v[1671] = v[425] * v[883] + v[439] * v[901] - v[444] * v[919] - v[454] * v[937];
	v[1670] = v[422] * v[883] + v[438] * v[901] - v[443] * v[919] - v[453] * v[937];
	v[938] = -(v[22] * v[805]) - v[23] * v[820] - v[24] * v[835] - v[25] * v[850];
	v[1676] = v[428] * v[884] + v[440] * v[902] - v[445] * v[920] - v[455] * v[938];
	v[1675] = v[425] * v[884] + v[439] * v[902] - v[444] * v[920] - v[454] * v[938];
	v[1674] = v[422] * v[884] + v[438] * v[902] - v[443] * v[920] - v[453] * v[938];
	v[939] = -(v[22] * v[806]) - v[23] * v[821] - v[24] * v[836] - v[25] * v[851];
	v[1680] = v[428] * v[885] + v[440] * v[903] - v[445] * v[921] - v[455] * v[939];
	v[1679] = v[425] * v[885] + v[439] * v[903] - v[444] * v[921] - v[454] * v[939];
	v[1678] = v[422] * v[885] + v[438] * v[903] - v[443] * v[921] - v[453] * v[939];
	v[940] = -(v[22] * v[807]) - v[23] * v[822] - v[24] * v[837] - v[25] * v[852];
	v[1743] = v[286] * v[923] + v[287] * v[924] + v[288] * v[925] + v[300] * v[926] + v[310] * v[927] + v[320] * v[928]
		+ v[289] * v[929] + v[290] * v[930] + v[291] * v[931] + v[326] * v[932] + v[336] * v[933] + v[346] * v[934] + v[292] * v[935]
		+ v[293] * v[936] + v[294] * v[937] + v[352] * v[938] + v[362] * v[939] + v[372] * v[940];
	v[1684] = v[428] * v[886] + v[440] * v[904] - v[445] * v[922] - v[455] * v[940];
	v[1683] = v[425] * v[886] + v[439] * v[904] - v[444] * v[922] - v[454] * v[940];
	v[1682] = v[422] * v[886] + v[438] * v[904] - v[443] * v[922] - v[453] * v[940];
	v[9121] = v[923];
	v[9122] = v[924];
	v[9123] = v[925];
	v[9124] = v[926];
	v[9125] = v[927];
	v[9126] = v[928];
	v[9127] = v[929];
	v[9128] = v[930];
	v[9129] = v[931];
	v[9130] = v[932];
	v[9131] = v[933];
	v[9132] = v[934];
	v[9133] = v[935];
	v[9134] = v[936];
	v[9135] = v[937];
	v[9136] = v[938];
	v[9137] = v[939];
	v[9138] = v[940];
	b941 = sqrt(Power(v[390] * v[398] - v[389] * v[399], 2) + Power(-(v[391] * v[398]) + v[389] * v[400], 2) + Power
	(v[391] * v[399] - v[390] * v[400], 2)) > 0.1e-7;
	if (b941) {
		v[943] = v[391] * v[399] - v[390] * v[400];
		v[944] = -(v[391] * v[398]) + v[389] * v[400];
		v[945] = v[390] * v[398] - v[389] * v[399];
		v[946] = sqrt((v[943] * v[943]) + (v[944] * v[944]) + (v[945] * v[945]));
		v[3091] = 1e0 / (v[946] * v[946]);
		v[2120] = v[946];
		v[3102] = 1e0 - (v[2120] * v[2120]);
		v[7307] = 1e0 / Power(v[3102], 0.15e1);
		v[3097] = 1e0 / sqrt(v[3102]);
		v[2119] = asin(v[2120]) / 2e0;
		v[3096] = 1e0 / Power(cos(v[2119]), 2);
		v[6687] = v[3096] * v[3097];
		v[948] = 2e0*tan(v[2119]);
		if (v[946] > 0.1e-7) { v013 = 1e0 / v[946]; v014 = (-(v013 / v[946])); v015 = (2e0*v013) / (v[946] * v[946]); }
		else {
			v013 = (12500000e0 / 3e0)*(24e0 - (-0.1e-7 + v[946])*(0.24e10 - 2e0*(-1e0 + 100000000e0*v[946])*
				(0.2399999997e10 - 0.1199999994e18*v[946] - 0.3e17*(v[946] * v[946]))));
			v014 = (-50000000e0 / 3e0)*(0.3599999994e10 - 0.4799999982e18*v[946] + 0.6e25*Power(v[946], 3)
				+ 0.1799999982e26*(v[946] * v[946]));
			v015 = 0.1e17*(799999997e0 - 0.599999994e17*v[946] - 0.3e17*(v[946] * v[946]));
		};
		v[952] = v015;
		v[953] = v014;
		v[954] = v013;
		v[7306] = v[948] * v[953] + v[6687] * v[954];
		v[6500] = v[948] * v[954];
		v[955] = v[6500] * v[943];
		v[6748] = 2e0*v[955];
		v[6553] = v[955] / 2e0;
		v[966] = (v[955] * v[955]);
		v[956] = v[6500] * v[944];
		v[6501] = v[956] / 2e0;
		v[964] = v[6501] * v[955];
		v[959] = (v[956] * v[956]);
		v[2074] = -v[959] - v[966];
		v[957] = v[6500] * v[945];
		v[6745] = 2e0*v[957];
		v[2104] = -v[957] + v[964];
		v[2095] = v[957] + v[964];
		v[971] = v[6501] * v[957];
		v[2086] = -v[955] + v[971];
		v[2078] = v[955] + v[971];
		v[969] = v[6553] * v[957];
		v[2099] = v[956] + v[969];
		v[2082] = -v[956] + v[969];
		v[960] = (v[957] * v[957]);
		v[2114] = 4e0 + v[959] + v[960] + v[966];
		v[7308] = 1e0 / Power(v[2114], 3);
		v[6744] = -4e0 / (v[2114] * v[2114]);
		v[2109] = -v[959] - v[960];
		v[2091] = -v[960] - v[966];
		v[958] = 4e0 / v[2114];
		v[6502] = v[958] / 2e0;
		v[961] = 1e0 + v[2109] * v[6502];
		v[962] = v[2104] * v[958];
		v[963] = v[2099] * v[958];
		v[965] = v[2095] * v[958];
		v[967] = 1e0 + v[2091] * v[6502];
		v[968] = v[2086] * v[958];
		v[970] = v[2082] * v[958];
		v[972] = v[2078] * v[958];
		v[973] = 1e0 + v[2074] * v[6502];
	}
	else {
		v[961] = 1e0;
		v[962] = 0e0;
		v[963] = 0e0;
		v[965] = 0e0;
		v[967] = 1e0;
		v[968] = 0e0;
		v[970] = 0e0;
		v[972] = 0e0;
		v[973] = 1e0;
	};
	if (b4) {
		v[2038] = 1e0 - v[1022];
		v[2036] = 1e0 - v[1018];
		v[2034] = 1e0 - v[1014];
		v[978] = v[221] * (v[217] + v[203] * v[225] + v[204] * v[226]) + v[222] * (v[220] + v[212] * v[225] + v[213] * v[226])
			- v[268] - v[263] * v[269] - v[264] * v[270] - v[265] * v[272] + v[1] * v[970] + v[2] * v[972] + v[3] * v[973];
		v[6503] = v[391] * v[978];
		v[977] = v[221] * (v[216] + v[200] * v[225] + v[201] * v[226]) + v[222] * (v[219] + v[209] * v[225] + v[210] * v[226])
			- v[267] - v[260] * v[269] - v[261] * v[270] - v[262] * v[272] + v[1] * v[965] + v[2] * v[967] + v[3] * v[968];
		v[6505] = v[390] * v[977];
		v[6550] = v[6503] + v[6505];
		v[976] = v[221] * (v[215] + v[197] * v[225] + v[198] * v[226]) + v[222] * (v[218] + v[206] * v[225] + v[207] * v[226])
			- v[266] - v[257] * v[269] - v[258] * v[270] - v[259] * v[272] + v[1] * v[961] + v[2] * v[962] + v[3] * v[963];
		v[6504] = -(v[389] * v[976]);
		v[6552] = v[6504] - v[6505];
		v[6551] = -v[6503] + v[6504];
		v[975] = -(v[389] * v[6550]) + v[2034] * v[976];
		v[979] = v[390] * v[6551] + v[2036] * v[977];
		v[980] = v[391] * v[6552] + v[2038] * v[978];
	}
	else {
		v[975] = 0e0;
		v[979] = 0e0;
		v[980] = 0e0;
	};
	v[1015] = v[1016] * v[287] + v[1017] * v[290] + v[1984] * v[300] + v[1982] * v[310] + v[1980] * v[320] + v[1978] * v[326]
		+ v[1976] * v[336] + v[1974] * v[346] + v[1972] * v[352] + v[1970] * v[362] + v[1968] * v[372] + v[1014] * v[4632]
		+ v[389] * v[4636];
	v[1021] = v[1016] * v[286] + v[1017] * v[289] + v[1905] * v[300] + v[1903] * v[310] + v[1901] * v[320] + v[1899] * v[326]
		+ v[1897] * v[336] + v[1895] * v[346] + v[1893] * v[352] + v[1891] * v[362] + v[1889] * v[372] + v[1018] * v[4600]
		+ v[390] * v[4602];
	v[1023] = v[1832] * v[300] + v[1830] * v[310] + v[1828] * v[320] + v[1826] * v[326] + v[1824] * v[336] + v[1822] * v[346]
		+ v[1820] * v[352] + v[1818] * v[362] + v[1816] * v[372] + v[3825] * v[389] + v[3826] * v[390] + v[1022] * v[4566];
	if (b6) {
		v[1785] = sqrt(v[36] * v[6506] * Power(v[401], v[1045]));
		v[7333] = 1e0 / (v[1785] * v[1785]);
		v[1037] = 2e0*v[1785] * v[31];
		v[6507] = v[26] * Power(v[401], v[29]);
		v[1032] = v[389] * v[6507];
		v[1034] = v[390] * v[6507];
		v[1035] = v[391] * v[6507];
		v[1036] = v[1015] * v[1037];
		v[1038] = v[1021] * v[1037];
		v[1039] = v[1023] * v[1037];
		b1040 = (v[1032] + v[1036])*v[389] + (v[1034] + v[1038])*v[390] + (v[1035] + v[1039])*v[391] > 0e0;
		if (b1040) {
			v[1042] = v[1036];
			v[1043] = v[1038];
			v[1044] = v[1039];
		}
		else {
			v[1042] = -v[1032];
			v[1043] = -v[1034];
			v[1044] = -v[1035];
		};
	}
	else {
		v[1062] = -1e0 + v[30];
		v[3606] = -1e0 + v[1062];
		v[1047] = v[27] - v[28];
		v[1046] = -((v[6506] * Power(v[1047], v[1045])) / (v[30] * Power(v[28], v[1062])));
		v[6509] = v[1046] * v[30];
		v[3801] = v[1062] * v[6509];
		if (b1049) {
			b1051 = v[401] > v[28];
			if (b1051) {
				v[1057] = v[27] - v[401];
				v[6508] = -(v[26] * Power(v[1057], v[29]));
				v[1032] = v[389] * v[6508];
				v[1034] = v[390] * v[6508];
				v[1035] = v[391] * v[6508];
			}
			else {
				v[1054] = -(v[26] * Power(v[1047], v[29])) + v[1046] * (Power(v[28], v[30]) - Power(v[401], v[30]));
				v[1032] = v[1054] * v[389];
				v[1034] = v[1054] * v[390];
				v[1035] = v[1054] * v[391];
			};
		}
		else {
			v[1032] = 0e0;
			v[1034] = 0e0;
			v[1035] = 0e0;
		};
		if (b1049) {
			if (b1051) {
				v[1799] = sqrt(v[36] * v[6506] * Power(v[1057], v[1045]));
				v[1059] = v[1799] * v[6734];
				v[1058] = v[1015] * v[1059];
				v[1060] = v[1021] * v[1059];
				v[1061] = v[1023] * v[1059];
			}
			else {
				v[1805] = sqrt(-(v[36] * v[6509] * Power(v[401], v[1062])));
				v[1063] = v[1805] * v[6734];
				v[1058] = v[1015] * v[1063];
				v[1060] = v[1021] * v[1063];
				v[1061] = v[1023] * v[1063];
			};
			b1064 = v[401] < v[27] && (v[1032] + v[1058])*v[389] + (v[1034] + v[1060])*v[390] + (v[1035] + v[1061]
				)*v[391] < 0e0;
			if (b1064) {
				v[1042] = v[1058];
				v[1043] = v[1060];
				v[1044] = v[1061];
			}
			else {
				v[1042] = -v[1032];
				v[1043] = -v[1034];
				v[1044] = -v[1035];
			};
		}
		else {
			v[1042] = 0e0;
			v[1043] = 0e0;
			v[1044] = 0e0;
		};
	};
	v[1066] = v[1032] + v[1042];
	v[1067] = v[1034] + v[1043];
	v[1068] = v[1035] + v[1044];
	v[3638] = (v[1066] * v[1066]) + (v[1067] * v[1067]) + (v[1068] * v[1068]);
	v[1069] = v[32] * v[975];
	v[1070] = v[32] * v[979];
	v[1071] = v[32] * v[980];
	v[1075] = v[1069] - v[33] * (v[1614] * v[286] + v[1618] * v[287] + v[1622] * v[288] + v[1638] * v[289] + v[1642] * v[290]
		+ v[1646] * v[291] + v[1662] * v[292] + v[1666] * v[293] + v[1670] * v[294] + v[1626] * v[300] + v[1630] * v[310]
		+ v[1634] * v[320] + v[1650] * v[326] + v[1654] * v[336] + v[1658] * v[346] + v[1674] * v[352] + v[1678] * v[362]
		+ v[1682] * v[372]);
	v[1076] = v[1070] - v[33] * (v[1615] * v[286] + v[1619] * v[287] + v[1623] * v[288] + v[1639] * v[289] + v[1643] * v[290]
		+ v[1647] * v[291] + v[1663] * v[292] + v[1667] * v[293] + v[1671] * v[294] + v[1627] * v[300] + v[1631] * v[310]
		+ v[1635] * v[320] + v[1651] * v[326] + v[1655] * v[336] + v[1659] * v[346] + v[1675] * v[352] + v[1679] * v[362]
		+ v[1683] * v[372]);
	v[1077] = v[1071] - v[33] * (v[1616] * v[286] + v[1620] * v[287] + v[1624] * v[288] + v[1640] * v[289] + v[1644] * v[290]
		+ v[1648] * v[291] + v[1664] * v[292] + v[1668] * v[293] + v[1672] * v[294] + v[1628] * v[300] + v[1632] * v[310]
		+ v[1636] * v[320] + v[1652] * v[326] + v[1656] * v[336] + v[1660] * v[346] + v[1676] * v[352] + v[1680] * v[362]
		+ v[1684] * v[372]);
	v[3634] = (v[1075] * v[1075]) + (v[1076] * v[1076]) + (v[1077] * v[1077]);
	if (b1049) {
		if (b5) {
			b1080 = sqrt((v[1075] * v[1075]) + (v[1076] * v[1076]) + (v[1077] * v[1077])) <= v[34] * sqrt((v[1066] * v[1066]) +
				(v[1067] * v[1067]) + (v[1068] * v[1068]));
			if (b1080) {
				v[1082] = v[1075];
				v[1083] = v[1076];
				v[1084] = v[1077];
				v[1085] = 1e0;
			}
			else {
				v[6511] = v[35] * sqrt(v[3638]);
				v[1086] = sqrt(v[3634]);
				if (v[1086] > 0.1e-5) {
					v016 = 1e0 / v[1086]; v017 = (-(v016 / v[1086])); v018 = (2e0*v016) / (v[1086] * v[1086]
						);
				}
				else {
					v016 = (24000000e0 - (-1e0 + 1000000e0*v[1086])*(71999994e0 - 0.71999982e14*v[1086] + 0.6e19*Power(v[1086]
						, 3) + 0.23999982e20*(v[1086] * v[1086]))) / 24e0;
					v017 = (-500000e0 / 3e0)*(35999994e0 - 0.47999982e14*v[1086] + 0.6e19*Power(v[1086], 3) + 0.17999982e20*
						(v[1086] * v[1086]));
					v018 = 0.1e13*(7999997e0 - 0.5999994e13*v[1086] - 0.3e13*(v[1086] * v[1086]));
				};
				v[1090] = v017;
				v[1091] = v016;
				v[1092] = v[1075] * v[1091];
				v[1093] = v[1076] * v[1091];
				v[1094] = v[1077] * v[1091];
				v[1082] = v[1092] * v[6511];
				v[1083] = v[1093] * v[6511];
				v[1084] = v[1094] * v[6511];
				v[1085] = 0e0;
			};
			if (sqrt((v[1069] * v[1069]) + (v[1070] * v[1070]) + (v[1071] * v[1071])) > v[34] * sqrt((v[1066] * v[1066]) +
				(v[1067] * v[1067]) + (v[1068] * v[1068]))) {
				if (v[32] > 0.1e-5) { v019 = 1e0 / v[32]; v020 = (-(v019 / v[32])); v021 = (2e0*v019) / (v[32] * v[32]); }
				else {
					v019 = (24000000e0 - (-1e0 + 1000000e0*v[32])*(71999994e0 - 0.71999982e14*v[32] + 0.6e19*Power(v[32], 3)
						+ 0.23999982e20*(v[32] * v[32]))) / 24e0;
					v020 = (-500000e0 / 3e0)*(35999994e0 - 0.47999982e14*v[32] + 0.6e19*Power(v[32], 3) + 0.17999982e20*
						(v[32] * v[32]));
					v021 = 0.1e13*(7999997e0 - 0.5999994e13*v[32] - 0.3e13*(v[32] * v[32]));
				};
				v[1104] = sqrt((v[1069] * v[1069]) + (v[1070] * v[1070]) + (v[1071] * v[1071]));
				if (v[1104] > 0.1e-5) {
					v022 = 1e0 / v[1104]; v023 = (-(v022 / v[1104])); v024 = (2e0*v022) / (v[1104] * v[1104]
						);
				}
				else {
					v022 = (24000000e0 - (-1e0 + 1000000e0*v[1104])*(71999994e0 - 0.71999982e14*v[1104] + 0.6e19*Power(v[1104]
						, 3) + 0.23999982e20*(v[1104] * v[1104]))) / 24e0;
					v023 = (-500000e0 / 3e0)*(35999994e0 - 0.47999982e14*v[1104] + 0.6e19*Power(v[1104], 3) + 0.17999982e20*
						(v[1104] * v[1104]));
					v024 = 0.1e13*(7999997e0 - 0.5999994e13*v[1104] - 0.3e13*(v[1104] * v[1104]));
				};
				v[1111] = -(v019*v022*v[35] * sqrt(v[3638]));
				v[1110] = v[1069] * v[1111] + v[975];
				v[1112] = v[1070] * v[1111] + v[979];
				v[1113] = v[1071] * v[1111] + v[980];
			}
			else {
				v[1110] = 0e0;
				v[1112] = 0e0;
				v[1113] = 0e0;
			};
		}
		else {
			b1114 = sqrt((v[1075] * v[1075]) + (v[1076] * v[1076]) + (v[1077] * v[1077])) <= v[35] * sqrt((v[1066] * v[1066]) +
				(v[1067] * v[1067]) + (v[1068] * v[1068]));
			if (b1114) {
				v[1082] = v[1075];
				v[1083] = v[1076];
				v[1084] = v[1077];
				v[1085] = 1e0;
			}
			else {
				v[1125] = sqrt(v[3638]);
				v[6512] = v[1125] * v[35];
				v[1116] = sqrt(v[3634]);
				if (v[1116] > 0.1e-5) {
					v025 = 1e0 / v[1116]; v026 = (-(v025 / v[1116])); v027 = (2e0*v025) / (v[1116] * v[1116]
						);
				}
				else {
					v025 = (24000000e0 - (-1e0 + 1000000e0*v[1116])*(71999994e0 - 0.71999982e14*v[1116] + 0.6e19*Power(v[1116]
						, 3) + 0.23999982e20*(v[1116] * v[1116]))) / 24e0;
					v026 = (-500000e0 / 3e0)*(35999994e0 - 0.47999982e14*v[1116] + 0.6e19*Power(v[1116], 3) + 0.17999982e20*
						(v[1116] * v[1116]));
					v027 = 0.1e13*(7999997e0 - 0.5999994e13*v[1116] - 0.3e13*(v[1116] * v[1116]));
				};
				v[1120] = v026;
				v[1121] = v025;
				v[1122] = v[1075] * v[1121];
				v[1123] = v[1076] * v[1121];
				v[1124] = v[1077] * v[1121];
				v[1082] = v[1122] * v[6512];
				v[1083] = v[1123] * v[6512];
				v[1084] = v[1124] * v[6512];
				v[1085] = 0e0;
			};
			if (sqrt((v[1069] * v[1069]) + (v[1070] * v[1070]) + (v[1071] * v[1071])) > v[35] * sqrt((v[1066] * v[1066]) +
				(v[1067] * v[1067]) + (v[1068] * v[1068]))) {
				if (v[32] > 0.1e-5) { v028 = 1e0 / v[32]; v029 = (-(v028 / v[32])); v030 = (2e0*v028) / (v[32] * v[32]); }
				else {
					v028 = (24000000e0 - (-1e0 + 1000000e0*v[32])*(71999994e0 - 0.71999982e14*v[32] + 0.6e19*Power(v[32], 3)
						+ 0.23999982e20*(v[32] * v[32]))) / 24e0;
					v029 = (-500000e0 / 3e0)*(35999994e0 - 0.47999982e14*v[32] + 0.6e19*Power(v[32], 3) + 0.17999982e20*
						(v[32] * v[32]));
					v030 = 0.1e13*(7999997e0 - 0.5999994e13*v[32] - 0.3e13*(v[32] * v[32]));
				};
				v[1134] = sqrt((v[1069] * v[1069]) + (v[1070] * v[1070]) + (v[1071] * v[1071]));
				if (v[1134] > 0.1e-5) {
					v031 = 1e0 / v[1134]; v032 = (-(v031 / v[1134])); v033 = (2e0*v031) / (v[1134] * v[1134]
						);
				}
				else {
					v031 = (24000000e0 - (-1e0 + 1000000e0*v[1134])*(71999994e0 - 0.71999982e14*v[1134] + 0.6e19*Power(v[1134]
						, 3) + 0.23999982e20*(v[1134] * v[1134]))) / 24e0;
					v032 = (-500000e0 / 3e0)*(35999994e0 - 0.47999982e14*v[1134] + 0.6e19*Power(v[1134], 3) + 0.17999982e20*
						(v[1134] * v[1134]));
					v033 = 0.1e13*(7999997e0 - 0.5999994e13*v[1134] - 0.3e13*(v[1134] * v[1134]));
				};
				v[1140] = -(v028*v031*v[35] * sqrt(v[3638]));
				v[1110] = v[1069] * v[1140] + v[975];
				v[1112] = v[1070] * v[1140] + v[979];
				v[1113] = v[1071] * v[1140] + v[980];
			}
			else {
				v[1110] = 0e0;
				v[1112] = 0e0;
				v[1113] = 0e0;
			};
		};
	}
	else {
		v[1082] = 0e0;
		v[1083] = 0e0;
		v[1084] = 0e0;
	};
	fn[0] = v[1066];
	fn[1] = v[1067];
	fn[2] = v[1068];
	ft[0] = v[1082];
	ft[1] = v[1083];
	ft[2] = v[1084];
	(*stickupdated) = v[1085];
	gtpupdated[0] = -v[1110] + v[975];
	gtpupdated[1] = -v[1112] + v[979];
	gtpupdated[2] = -v[1113] + v[980];
	v[1151] = v[1035];
	v[1173] = v[1151];
	v[1152] = v[1034];
	v[1175] = v[1152];
	v[1153] = v[1032];
	v[1177] = v[1153];
	b1157 = b6;
	if (b1157) {
	}
	else {
		b1161 = b1049;
		if (b1161) {
			b1162 = b1051;
		}
		else {
		};
	};
	b1166 = b6;
	if (b1166) {
		v[1167] = 0e0;
	}
	else {
		b1168 = b418;
		if (b1168) {
			v[1169] = 0e0;
			v[1170] = 0e0;
			v[1171] = 0e0;
		}
		else {
			v[1171] = 0e0;
			v[1170] = 0e0;
			v[1169] = 0e0;
		};
		v[1151] = v[1173] + v[1169] * v[407];
		v[1152] = v[1175] + v[1170] * v[407];
		v[1176] = v[1171] * v[373] + v[1170] * v[374] + v[1169] * v[375];
		v[1153] = v[1177] + v[1171] * v[407];
		v[1167] = v[1176] * v[406];
	};
	v[6513] = v[1167] / v[401];
	v[1151] = v[1151] + v[375] * v[6513];
	v[1152] = v[1152] + v[374] * v[6513];
	v[1153] = v[1153] + v[373] * v[6513];
	v[7261] = -(v[1153] * v[198]) - v[1152] * v[201] - v[1151] * v[204];
	v[7260] = -(v[1153] * v[197]) - v[1152] * v[200] - v[1151] * v[203];
	v[7259] = v[1153] * v[207] + v[1152] * v[210] + v[1151] * v[213];
	v[7258] = v[1153] * v[206] + v[1152] * v[209] + v[1151] * v[212];
	v[1180] = v[1151] * v[981];
	v[6525] = v[1180] / 2e0;
	v[1181] = v[1151] * v[982];
	v[1272] = v[1181] * v[238];
	v[1182] = v[1151] * v[983];
	v[1277] = v[1182] * v[238];
	v[1187] = v[1151] * v[984];
	v[6534] = v[1187] / 2e0;
	v[1188] = v[1151] * v[985];
	v[1262] = v[1188] * v[178];
	v[1189] = v[1151] * v[986];
	v[1267] = v[1189] * v[178];
	v[1190] = v[1151] * v[987];
	v[6537] = v[1190] / 2e0;
	v[1191] = v[1151] * v[988];
	v[1252] = v[1191] * v[159];
	v[1192] = v[1151] * v[989];
	v[1257] = v[1192] * v[159];
	v[1196] = v[1151] * v[224];
	v[1197] = v[1151] * v[223];
	v[1199] = v[1152] * v[981];
	v[1273] = v[1199] * v[238];
	v[1200] = v[1152] * v[982];
	v[6524] = v[1200] / 2e0;
	v[1275] = v[1200] * v[6448];
	v[1201] = v[1152] * v[983];
	v[1280] = v[1201] * v[238];
	v[1206] = v[1152] * v[984];
	v[1263] = v[1206] * v[178];
	v[1207] = v[1152] * v[985];
	v[6533] = v[1207] / 2e0;
	v[1265] = v[1207] * v[6441];
	v[1208] = v[1152] * v[986];
	v[1270] = v[1208] * v[178];
	v[1209] = v[1152] * v[987];
	v[1253] = v[1209] * v[159];
	v[1210] = v[1152] * v[988];
	v[6536] = v[1210] / 2e0;
	v[1255] = v[1210] * v[6434];
	v[1211] = v[1152] * v[989];
	v[1260] = v[1211] * v[159];
	v[1215] = v[1152] * v[224];
	v[1216] = v[1152] * v[223];
	v[1217] = -(v[1153] * v[258]) - v[1152] * v[261] - v[1151] * v[264];
	v[1218] = v[1153] * v[981];
	v[1278] = v[1218] * v[238];
	v[1219] = v[1153] * v[982];
	v[1281] = v[1219] * v[238];
	v[1220] = v[1153] * v[983];
	v[6523] = v[1220] / 2e0;
	v[1279] = v[1220] * v[6448];
	v[1221] = v[1153] * v[1394] + v[1152] * v[1399] + v[1151] * v[1404];
	v[1222] = v[1153] * v[1393] + v[1152] * v[1398] + v[1151] * v[1403];
	v[1225] = v[1153] * v[984];
	v[1268] = v[1225] * v[178];
	v[1226] = v[1153] * v[985];
	v[1271] = v[1226] * v[178];
	v[1227] = v[1153] * v[986];
	v[6532] = v[1227] / 2e0;
	v[1269] = v[1227] * v[6441];
	v[1228] = v[1153] * v[987];
	v[1258] = v[1228] * v[159];
	v[1229] = v[1153] * v[988];
	v[1261] = v[1229] * v[159];
	v[1230] = v[1153] * v[989];
	v[6535] = v[1230] / 2e0;
	v[1259] = v[1230] * v[6434];
	v[1231] = -(v[1153] * v[257]) - v[1152] * v[260] - v[1151] * v[263];
	v[1232] = -(v[1153] * v[259]) - v[1152] * v[262] - v[1151] * v[265];
	v[1570] = -(v[1025] * v[1231]) + v[1232] * v[1611];
	v[1234] = v[1153] * v[224];
	v[1235] = v[1153] * v[223];
	v[1236] = v[1272] + v[1273];
	v[1237] = v[1277] + v[1278];
	v[1238] = v[1280] + v[1281];
	v[1239] = v[6523] * v[694] + v[1219] * v[699] + v[1218] * v[705] + v[1201] * v[709] + v[6524] * v[713] + v[1199] * v[718]
		+ v[1182] * v[722] + v[1181] * v[726] + v[6525] * v[731];
	v[1603] = v[1275] + v[1279] - v[1239] * v[6432];
	v[1600] = -v[1275] + v[1603] + v[1180] * v[6448];
	v[1597] = v[1275] - v[1279] + v[1600];
	v[1240] = -(v[1024] * v[1222]) + v[1221] * v[1610];
	v[1241] = (-(v[1153] * v[420]) + v[1153] * v[421] - v[1152] * v[423] + v[1152] * v[424] - v[1151] * v[426]
		+ v[1151] * v[427]) / 2e0;
	v[1242] = v[1262] + v[1263];
	v[1243] = v[1267] + v[1268];
	v[1244] = v[1270] + v[1271];
	v[1245] = v[1226] * v[528] + v[1225] * v[534] + v[1208] * v[538] + v[1206] * v[547] + v[1189] * v[551] + v[1188] * v[555]
		+ v[523] * v[6532] + v[542] * v[6533] + v[560] * v[6534];
	v[1596] = v[1265] + v[1269] - v[1245] * v[6430];
	v[1593] = -v[1265] + v[1596] + v[1187] * v[6441];
	v[1590] = v[1265] - v[1269] + v[1593];
	v[1246] = v[1252] + v[1253];
	v[1247] = v[1257] + v[1258];
	v[1248] = v[1260] + v[1261];
	v[1249] = v[1229] * v[483] + v[1228] * v[489] + v[1211] * v[493] + v[1209] * v[502] + v[1192] * v[506] + v[1191] * v[510]
		+ v[478] * v[6535] + v[497] * v[6536] + v[515] * v[6537];
	v[1589] = v[1255] + v[1259] - v[1249] * v[6428];
	v[1586] = -v[1255] + v[1589] + v[1190] * v[6434];
	v[1583] = v[1255] - v[1259] + v[1586];
	v[9103] = v[1235];
	v[9104] = v[1216];
	v[9105] = v[1197];
	v[9106] = v[1252] - v[1253] + v[1583] * v[457] + v[1248] * v[458] + v[1247] * v[461];
	v[9107] = -v[1257] + v[1258] + v[1248] * v[459] + v[1586] * v[460] + v[1246] * v[461];
	v[9108] = v[1260] - v[1261] + v[1246] * v[458] + v[1247] * v[459] + v[1589] * v[462];
	v[9109] = v[1234];
	v[9110] = v[1215];
	v[9111] = v[1196];
	v[9112] = v[1262] - v[1263] + v[1590] * v[463] + v[1244] * v[464] + v[1243] * v[467];
	v[9113] = -v[1267] + v[1268] + v[1244] * v[465] + v[1593] * v[466] + v[1242] * v[467];
	v[9114] = v[1270] - v[1271] + v[1242] * v[464] + v[1243] * v[465] + v[1596] * v[468];
	v[9115] = -v[1153];
	v[9116] = -v[1152];
	v[9117] = -v[1151];
	v[9118] = v[1272] - v[1273] + v[1597] * v[469] + v[1238] * v[470] + v[1237] * v[473];
	v[9119] = -v[1277] + v[1278] + v[1238] * v[471] + v[1600] * v[472] + v[1236] * v[473];
	v[9120] = v[1280] - v[1281] + v[1236] * v[470] + v[1237] * v[471] + v[1603] * v[474];
	v[1250] = v[1027] * v[1231] - v[1232] * v[1613];
	v[1251] = v[1026] * v[1570] + v[1217] * v[1612];
	for (i1149 = 1; i1149 <= 18; i1149++) {
		b6531 = i1149 == 1;
		b6530 = i1149 == 2;
		b6529 = i1149 == 3;
		b6528 = i1149 == 7;
		b6527 = i1149 == 8;
		b6526 = i1149 == 9;
		v[1349] = (i1149 == 16 ? 1 : 0);
		v[1346] = (i1149 == 17 ? 1 : 0);
		v[1343] = (i1149 == 18 ? 1 : 0);
		v[1340] = (i1149 == 10 ? 1 : 0);
		v[1337] = (i1149 == 11 ? 1 : 0);
		v[1334] = (i1149 == 12 ? 1 : 0);
		v[1331] = (i1149 == 4 ? 1 : 0);
		v[1328] = (i1149 == 5 ? 1 : 0);
		v[1325] = (i1149 == 6 ? 1 : 0);
		v[1309] = v[9174 + i1149];
		v[1307] = v[9156 + i1149];
		v[1302] = v[9120 + i1149];
		v[6538] = v[1026] * v[1302];
		v[1301] = v[9138 + i1149];
		v[1287] = v[9214 + i1149];
		v[6514] = -(v[1287] * v[6428]);
		v[1364] = -2e0*v[1287] * v[1351];
		v[6517] = 2e0*v[1364];
		v[1288] = v[9232 + i1149];
		v[1289] = v[9250 + i1149];
		v[1290] = v[9268 + i1149];
		v[1291] = v[9286 + i1149];
		v[6515] = -(v[1291] * v[6430]);
		v[1368] = -2e0*v[1291] * v[1355];
		v[6518] = 2e0*v[1368];
		v[1292] = v[9304 + i1149];
		v[1293] = v[9322 + i1149];
		v[1294] = v[9340 + i1149];
		v[1295] = v[9358 + i1149];
		v[6516] = -(v[1295] * v[6432]);
		v[1372] = -2e0*v[1295] * v[1359];
		v[6519] = 2e0*v[1372];
		v[1296] = v[9376 + i1149];
		v[1297] = v[9394 + i1149];
		v[1298] = v[9412 + i1149];
		v[1303] = -(v[1301] * v[1613]) + v[1302] * v[451];
		v[1304] = v[1027] * v[1301] + v[1302] * v[449];
		v[1305] = v[9484 + i1149];
		v[1306] = v[9538 + i1149];
		v[1308] = -0.5e0*v[1307];
		v[1493] = v[1153] * v[1308];
		v[1481] = v[1152] * v[1308];
		v[1463] = v[1151] * v[1308];
		v[1310] = -(v[1024] * v[1309]);
		v[1311] = v[1309] * v[1610];
		v[1312] = v[9592 + i1149];
		v[1313] = v[1302] * v[1612];
		v[1316] = v[9718 + i1149];
		v[1317] = v[9772 + i1149];
		v[1318] = v[9826 + i1149];
		v[1321] = v[9952 + i1149];
		v[1322] = v[10006 + i1149];
		v[1323] = v[10060 + i1149];
		v[1324] = v[1288] - v[1325];
		v[1326] = v[1288] + v[1325];
		v[1327] = v[1289] + v[1328];
		v[1329] = v[1289] - v[1328];
		v[1330] = v[1290] - v[1331];
		v[1332] = v[1290] + v[1331];
		v[1333] = v[1292] - v[1334];
		v[1335] = v[1292] + v[1334];
		v[1336] = v[1293] + v[1337];
		v[1338] = v[1293] - v[1337];
		v[1339] = v[1294] - v[1340];
		v[1341] = v[1294] + v[1340];
		v[1342] = v[1296] - v[1343];
		v[1344] = v[1296] + v[1343];
		v[1345] = v[1297] + v[1346];
		v[1347] = v[1297] - v[1346];
		v[1348] = v[1298] - v[1349];
		v[1350] = v[1298] + v[1349];
		v[1352] = v[1364] * v[478] + v[1305] * v[6434];
		v[1353] = v[1324] * v[159] + v[483] * v[6514];
		v[1354] = v[1327] * v[159] + v[489] * v[6514];
		v[1356] = v[1368] * v[523] + v[1306] * v[6441];
		v[1357] = v[1333] * v[178] + v[528] * v[6515];
		v[1358] = v[1336] * v[178] + v[534] * v[6515];
		v[1360] = v[1312] * v[6448] + v[1372] * v[694];
		v[1361] = v[1342] * v[238] + v[6516] * v[699];
		v[1362] = v[1345] * v[238] + v[6516] * v[705];
		v[1396] = (b6528 ? v[224] : 0e0) + (b6531 ? v[223] : 0e0) + (i1149 == 13 ? -1 : 0) + v[1310] * v[1393] + v[1311] * v[1394]
			- v[1304] * v[257] - v[1313] * v[258] - v[1303] * v[259] - v[1308] * v[6478] + v[1362] * v[981] + v[1361] * v[982]
			+ v[1360] * v[983] + v[1358] * v[984] + v[1357] * v[985] + v[1356] * v[986] + v[1354] * v[987] + v[1353] * v[988]
			+ v[1352] * v[989];
		v[1363] = v[1326] * v[159] + v[493] * v[6514];
		v[1365] = v[1364] * v[497] + v[1316] * v[6434];
		v[1366] = v[1330] * v[159] + v[502] * v[6517];
		v[1367] = v[1335] * v[178] + v[538] * v[6515];
		v[1369] = v[1368] * v[542] + v[1317] * v[6441];
		v[1370] = v[1339] * v[178] + v[547] * v[6518];
		v[1371] = v[1344] * v[238] + v[6516] * v[709];
		v[1373] = v[1318] * v[6448] + v[1372] * v[713];
		v[1374] = v[1348] * v[238] + v[6519] * v[718];
		v[1401] = (b6527 ? v[224] : 0e0) + (b6530 ? v[223] : 0e0) + (i1149 == 14 ? -1 : 0) + v[1310] * v[1398] + v[1311] * v[1399]
			- v[1304] * v[260] - v[1313] * v[261] - v[1303] * v[262] - v[1308] * v[6479] + v[1374] * v[981] + v[1373] * v[982]
			+ v[1371] * v[983] + v[1370] * v[984] + v[1369] * v[985] + v[1367] * v[986] + v[1366] * v[987] + v[1365] * v[988]
			+ v[1363] * v[989];
		v[1375] = v[1329] * v[159] + v[506] * v[6517];
		v[1376] = v[1153] * v[1352] + v[1152] * v[1363] + v[1151] * v[1375];
		v[1377] = v[1332] * v[159] + v[510] * v[6517];
		v[1378] = v[1153] * v[1353] + v[1152] * v[1365] + v[1151] * v[1377];
		v[1379] = v[1364] * v[515] + v[1321] * v[6434];
		v[1380] = v[1153] * v[1354] + v[1152] * v[1366] + v[1151] * v[1379];
		v[1381] = v[1338] * v[178] + v[551] * v[6518];
		v[1382] = v[1153] * v[1356] + v[1152] * v[1367] + v[1151] * v[1381];
		v[1383] = v[1341] * v[178] + v[555] * v[6518];
		v[1384] = v[1153] * v[1357] + v[1152] * v[1369] + v[1151] * v[1383];
		v[1385] = v[1368] * v[560] + v[1322] * v[6441];
		v[1386] = v[1153] * v[1358] + v[1152] * v[1370] + v[1151] * v[1385];
		v[1387] = v[1347] * v[238] + v[6519] * v[722];
		v[1388] = v[1153] * v[1360] + v[1152] * v[1371] + v[1151] * v[1387];
		v[1389] = v[1350] * v[238] + v[6519] * v[726];
		v[1390] = v[1153] * v[1361] + v[1152] * v[1373] + v[1151] * v[1389];
		v[1391] = v[1323] * v[6448] + v[1372] * v[731];
		v[1407] = (b6526 ? v[224] : 0e0) + (b6529 ? v[223] : 0e0) + (i1149 == 15 ? -1 : 0) + v[1310] * v[1403] + v[1311] * v[1404]
			- v[1304] * v[263] - v[1313] * v[264] - v[1303] * v[265] - v[1308] * v[6480] + v[1391] * v[981] + v[1389] * v[982]
			+ v[1387] * v[983] + v[1385] * v[984] + v[1383] * v[985] + v[1381] * v[986] + v[1379] * v[987] + v[1377] * v[988]
			+ v[1375] * v[989];
		v[1392] = v[1153] * v[1362] + v[1152] * v[1374] + v[1151] * v[1391];
		v[1395] = v[1396] * v[6513];
		v[1397] = v[1396];
		v[1418] = v[1397];
		v[1400] = v[1401] * v[6513];
		v[1402] = v[1401];
		v[6520] = v[1397] * v[373] + v[1402] * v[374] + v[1407] * v[375];
		v[1419] = v[1402];
		v[1405] = v[6520] / v[401];
		v[1408] = -(v[1167] * v[1406] * v[6520]);
		v[1439] = v[1408];
		v[1409] = v[1407] * v[6513];
		v[1410] = v[1407];
		v[1420] = v[1410];
		v[1411] = 0e0;
		v[1412] = 0e0;
		b1413 = b6;
		if (b1413) {
			v[1414] = 0e0;
			v[1415] = 0e0;
			v[1416] = 0e0;
		}
		else {
			v[1417] = v[1405] * v[406];
			v[1411] = v[1176] * v[1405] * v[405];
			v[1416] = v[1418];
			v[1397] = 0e0;
			v[1395] = v[1395] + v[1171] * v[1417];
			v[1415] = v[1419];
			v[1402] = 0e0;
			v[1400] = v[1400] + v[1170] * v[1417];
			v[1414] = v[1420];
			v[1412] = v[1171] * v[1418] + v[1170] * v[1419] + v[1169] * v[1420];
			v[1410] = 0e0;
			v[1409] = v[1409] + v[1169] * v[1417];
		};
		v[1453] = v[1395];
		v[1452] = v[1400];
		v[1451] = v[1409];
		v[1397] = v[1397] + v[1416];
		v[1421] = v[1397];
		v[1437] = v[1421];
		v[1402] = v[1402] + v[1415];
		v[1422] = v[1402];
		v[1435] = v[1422];
		v[1410] = v[1410] + v[1414];
		v[1423] = v[1410];
		v[1433] = v[1423];
		v[1424] = 0e0;
		v[1425] = 0e0;
		v[1426] = 0e0;
		b1427 = b6;
		if (b1427) {
			v[1428] = v[1423] * v[6731];
			v[1426] = v[1423] * v[6507];
			v[1423] = 0e0;
			v[1429] = v[1428] + v[1422] * v[6732];
			v[1425] = v[1422] * v[6507];
			v[1422] = 0e0;
			v[1430] = v[1429] + v[1421] * v[6733];
			v[1424] = v[1421] * v[6507];
			v[1421] = 0e0;
			v[1408] = v[1408] + v[1430] * v[29] * v[7257];
		}
		else {
			b1431 = b1049;
			if (b1431) {
				v[6521] = -(v[1437] * v[389]) - v[1435] * v[390] - v[1433] * v[391];
				b1432 = b1051;
				if (b1432) {
					v[1426] = v[1433] * v[6508];
					v[1423] = 0e0;
					v[1425] = v[1435] * v[6508];
					v[1422] = 0e0;
					v[1424] = v[1437] * v[6508];
					v[1421] = 0e0;
					v[1408] = v[1439] - v[6506] * v[6521] * Power(v[1057], v[1045]);
				}
				else {
					v[1426] = v[1054] * v[1433];
					v[1423] = 0e0;
					v[1425] = v[1054] * v[1435];
					v[1422] = 0e0;
					v[1424] = v[1054] * v[1437];
					v[1421] = 0e0;
					v[1408] = v[1439] + v[6509] * v[6521] * Power(v[401], v[1062]);
				};
			}
			else {
			};
		};
		v[1454] = v[1408];
		v[1449] = v[1424];
		v[1447] = v[1425];
		v[1445] = v[1426];
		b1440 = b6;
		if (b1440) {
			v[1441] = v[1426] * v[375];
			v[1409] = v[1409] + v[1426] * v[388];
			v[1426] = 0e0;
			v[1442] = v[1441] + v[1425] * v[374];
			v[1400] = v[1400] + v[1425] * v[388];
			v[1425] = 0e0;
			v[1443] = v[1442] + v[1424] * v[373];
			v[1395] = v[1395] + v[1424] * v[388];
			v[1424] = 0e0;
			v[1408] = v[1408] + v[1443] * v[387];
		}
		else {
			b1444 = b418;
			if (b1444) {
				v[1446] = -v[1445];
				v[1426] = 0e0;
				v[1448] = -v[1447];
				v[1425] = 0e0;
				v[1450] = -v[1449];
				v[1424] = 0e0;
			}
			else {
				v[1446] = v[1445];
				v[1448] = v[1447];
				v[1450] = v[1449];
			};
			v[1412] = v[1412] + v[1446] * v[375];
			v[1409] = v[1451] + v[1446] * v[407];
			v[1412] = v[1412] + v[1448] * v[374];
			v[1400] = v[1452] + v[1448] * v[407];
			v[1412] = v[1412] + v[1450] * v[373];
			v[1395] = v[1453] + v[1450] * v[407];
			v[1411] = v[1411] + v[1412] * v[406];
			v[1408] = v[1411] + v[1454];
		};
		v[6522] = v[1408] / v[401];
		v[1409] = v[1409] + v[375] * v[6522];
		v[1400] = v[1400] + v[374] * v[6522];
		v[1395] = v[1395] + v[373] * v[6522];
		v[1459] = -(v[1151] * v[1303]) - v[1027] * v[1409];
		v[1460] = -(v[1151] * v[1313]) - v[1409] * v[274];
		v[1461] = -(v[1151] * v[1304]) - v[1409] * v[1613];
		v[1462] = v[1463] + v[1409] * v[223];
		v[1464] = -v[1463] + v[1409] * v[224];
		v[1477] = -(v[1152] * v[1303]) - v[1027] * v[1400];
		v[1478] = -(v[1152] * v[1313]) - v[1400] * v[274];
		v[1479] = -(v[1152] * v[1304]) - v[1400] * v[1613];
		v[1480] = v[1481] + v[1400] * v[223];
		v[1482] = -v[1481] + v[1400] * v[224];
		v[1489] = -(v[1153] * v[1303]) - v[1027] * v[1395];
		v[1490] = -(v[1153] * v[1313]) - v[1395] * v[274];
		v[1491] = -(v[1153] * v[1304]) - v[1395] * v[1613];
		v[1492] = v[1493] + v[1395] * v[223];
		v[1494] = -v[1493] + v[1395] * v[224];
		v[1497] = -(v[1232] * v[1301]) - v[128] * v[1388] - v[131] * v[1390] - v[134] * v[1392] - v[1395] * v[257]
			- v[1400] * v[260] - v[1409] * v[263];
		v[1498] = v[1231] * v[1301] - v[130] * v[1388] - v[133] * v[1390] - v[136] * v[1392] - v[1395] * v[259] - v[1400] * v[262]
			- v[1409] * v[265];
		v[1500] = v[136] * v[1459] + v[135] * v[1460] + v[134] * v[1461];
		v[1501] = v[133] * v[1459] + v[132] * v[1460] + v[131] * v[1461];
		v[1502] = v[130] * v[1459] + v[129] * v[1460] + v[128] * v[1461];
		v[1503] = v[136] * v[1477] + v[135] * v[1478] + v[134] * v[1479];
		v[1504] = v[133] * v[1477] + v[132] * v[1478] + v[131] * v[1479];
		v[1505] = v[130] * v[1477] + v[129] * v[1478] + v[128] * v[1479];
		v[1506] = v[136] * v[1489] + v[135] * v[1490] + v[134] * v[1491];
		v[1507] = v[133] * v[1489] + v[132] * v[1490] + v[131] * v[1491];
		v[1508] = v[130] * v[1489] + v[129] * v[1490] + v[128] * v[1491];
		v[1509] = v[1180] * v[1372] - v[1500] * v[6448];
		v[1510] = v[1501] * v[238] + v[1181] * v[6519];
		v[1511] = v[1502] * v[238] + v[1182] * v[6519];
		v[1512] = v[1503] * v[238] + v[1199] * v[6519];
		v[1514] = v[1505] * v[238] + v[1201] * v[6519];
		v[1515] = v[1506] * v[238] + v[1218] * v[6519];
		v[1516] = v[1507] * v[238] + v[1219] * v[6519];
		v[6541] = v[1239] * v[1295] * v[6431] - v[6432] * (v[1219] * v[1342] + v[1201] * v[1344] + v[1218] * v[1345]
			+ v[1182] * v[1347] + v[1199] * v[1348] + v[1181] * v[1350] - v[1312] * v[6523] - v[1318] * v[6524] - v[1323] * v[6525]
			+ v[1500] * v[6575] + v[1504] * v[6576] + v[1508] * v[6577] + v[1507] * v[699] + v[1506] * v[705] + v[1505] * v[709]
			+ v[1503] * v[718] + v[1502] * v[722] + v[1501] * v[726]);
		v[1599] = -(v[1200] * v[1372]) + v[1504] * v[6448] + v[6541];
		v[1518] = v[1220] * v[1372] - v[1508] * v[6448];
		v[1519] = v[1235] * v[1311] + v[1492] * v[228];
		v[1520] = v[1235] * v[1310] + v[1492] * v[227];
		v[1521] = v[1234] * v[1311] + v[1494] * v[228];
		v[1522] = v[1234] * v[1310] + v[1494] * v[227];
		v[1523] = v[1216] * v[1311] + v[1480] * v[228];
		v[1524] = v[1216] * v[1310] + v[1480] * v[227];
		v[1525] = v[1215] * v[1311] + v[1482] * v[228];
		v[1526] = v[1215] * v[1310] + v[1482] * v[227];
		v[1527] = v[1197] * v[1311] + v[1462] * v[228];
		v[1528] = v[1197] * v[1310] + v[1462] * v[227];
		v[1529] = v[1196] * v[1311] + v[1464] * v[228];
		v[1530] = v[1196] * v[1310] + v[1464] * v[227];
		v[1531] = ((b6526 ? v[1151] : 0e0) + (b6527 ? v[1152] : 0e0) + (b6528 ? v[1153] : 0e0) + (b6529 ? -v[1151] : 0e0) + (b6530 ?
			-v[1152] : 0e0) + (b6531 ? -v[1153] : 0e0) + v[1382] * v[1465] + v[1384] * v[1466] + v[1386] * v[1467] - v[1376] * v[1469]
			- v[1378] * v[1470] - v[1380] * v[1471] + v[1395] * (-v[420] + v[421]) + v[1400] * (-v[423] + v[424]) + v[1409] * (
				-v[426] + v[427]) + v[1310] * (v[7258] + v[7260]) + v[1311] * (v[7259] + v[7261])) / 2e0;
		v[1532] = v[1530] * v[97] + v[1529] * v[98];
		v[1533] = v[1530] * v[94] + v[1529] * v[95];
		v[1534] = v[1530] * v[91] + v[1529] * v[92];
		v[1535] = v[1526] * v[97] + v[1525] * v[98];
		v[1536] = v[1526] * v[94] + v[1525] * v[95];
		v[1537] = v[1526] * v[91] + v[1525] * v[92];
		v[1538] = v[1522] * v[97] + v[1521] * v[98];
		v[1539] = v[1522] * v[94] + v[1521] * v[95];
		v[1540] = v[1522] * v[91] + v[1521] * v[92];
		v[1541] = v[1528] * v[88] + v[1527] * v[89];
		v[1542] = v[1528] * v[85] + v[1527] * v[86];
		v[1543] = v[1528] * v[82] + v[1527] * v[83];
		v[1544] = v[1524] * v[88] + v[1523] * v[89];
		v[1545] = v[1524] * v[85] + v[1523] * v[86];
		v[1546] = v[1524] * v[82] + v[1523] * v[83];
		v[1547] = v[1520] * v[88] + v[1519] * v[89];
		v[1548] = v[1520] * v[85] + v[1519] * v[86];
		v[1549] = v[1520] * v[82] + v[1519] * v[83];
		v[1550] = v[1187] * v[1368] - v[1532] * v[6441];
		v[1551] = v[1533] * v[178] + v[1188] * v[6518];
		v[1552] = v[1534] * v[178] + v[1189] * v[6518];
		v[1553] = v[1535] * v[178] + v[1206] * v[6518];
		v[1555] = v[1537] * v[178] + v[1208] * v[6518];
		v[1556] = v[1538] * v[178] + v[1225] * v[6518];
		v[1557] = v[1539] * v[178] + v[1226] * v[6518];
		v[6540] = v[1245] * v[1291] * v[6429] - v[6430] * (v[1226] * v[1333] + v[1208] * v[1335] + v[1225] * v[1336]
			+ v[1189] * v[1338] + v[1206] * v[1339] + v[1188] * v[1341] + v[1539] * v[528] + v[1538] * v[534] + v[1537] * v[538]
			+ v[1535] * v[547] + v[1534] * v[551] + v[1533] * v[555] - v[1306] * v[6532] - v[1317] * v[6533] - v[1322] * v[6534]
			+ v[1532] * v[6605] + v[1536] * v[6606] + v[1540] * v[6607]);
		v[1592] = -(v[1207] * v[1368]) + v[1536] * v[6441] + v[6540];
		v[1559] = v[1227] * v[1368] - v[1540] * v[6441];
		v[1560] = v[1190] * v[1364] - v[1541] * v[6434];
		v[1561] = v[1542] * v[159] + v[1191] * v[6517];
		v[1562] = v[1543] * v[159] + v[1192] * v[6517];
		v[1563] = v[1544] * v[159] + v[1209] * v[6517];
		v[1565] = v[1546] * v[159] + v[1211] * v[6517];
		v[1566] = v[1547] * v[159] + v[1228] * v[6517];
		v[1567] = v[1548] * v[159] + v[1229] * v[6517];
		v[6539] = v[1249] * v[1287] * v[6427] - v[6428] * (v[1229] * v[1324] + v[1211] * v[1326] + v[1228] * v[1327]
			+ v[1192] * v[1329] + v[1209] * v[1330] + v[1191] * v[1332] + v[1548] * v[483] + v[1547] * v[489] + v[1546] * v[493]
			+ v[1544] * v[502] + v[1543] * v[506] + v[1542] * v[510] - v[1305] * v[6535] - v[1316] * v[6536] - v[1321] * v[6537]
			+ v[1541] * v[6623] + v[1545] * v[6624] + v[1549] * v[6625]);
		v[1585] = -(v[1210] * v[1364]) + v[1545] * v[6434] + v[6539];
		v[1569] = v[1230] * v[1364] - v[1549] * v[6434];
		v[1571] = -(v[1026] * (v[1217] * v[1302] + v[1025] * v[1497] - v[1498] * v[1611])) - v[1612] * (v[129] * v[1388]
			+ v[132] * v[1390] + v[135] * v[1392] - v[1302] * v[1570] + v[1395] * v[258] + v[1400] * v[261] + v[1409] * v[264]);
		v[1572] = v[1611] * (-(v[1497] * v[275]) + v[1231] * v[6538]) + v[1025] * (-(v[1498] * v[275]) + v[1232] * v[6538]);
		v[1573] = v[1511] + v[1515];
		v[1574] = v[1510] + v[1512];
		v[1575] = v[1514] + v[1516];
		v[1576] = v[1024] * (-(v[1221] * v[1309]) - v[1492] * v[197] - v[1480] * v[200] - v[1462] * v[203] - v[1494] * v[206]
			- v[1482] * v[209] - v[1464] * v[212] + v[223] * (-(v[1376] * v[82]) - v[1378] * v[85] - v[1380] * v[88]) + v[224] * (-
			(v[1382] * v[91]) - v[1384] * v[94] - v[1386] * v[97])) + v[1610] * (-(v[1222] * v[1309]) + v[1492] * v[198]
				+ v[1480] * v[201] + v[1462] * v[204] + v[1494] * v[207] + v[1482] * v[210] + v[1464] * v[213] + v[223] * (v[1376] * v[83]
					+ v[1378] * v[86] + v[1380] * v[89]) + v[224] * (v[1382] * v[92] + v[1384] * v[95] + v[1386] * v[98]));
		v[1577] = v[1552] + v[1556];
		v[1578] = v[1551] + v[1553];
		v[1579] = v[1555] + v[1557];
		v[1580] = v[1562] + v[1566];
		v[1581] = v[1561] + v[1563];
		v[1582] = v[1565] + v[1567];
		v[10367] = 0e0;
		v[10368] = 0e0;
		v[10369] = 0e0;
		v[10370] = 0e0;
		v[10371] = v[1248];
		v[10372] = v[1247];
		v[10373] = 0e0;
		v[10374] = 0e0;
		v[10375] = 0e0;
		v[10376] = 0e0;
		v[10377] = 0e0;
		v[10378] = 0e0;
		v[10379] = 0e0;
		v[10380] = 0e0;
		v[10381] = 0e0;
		v[10382] = 0e0;
		v[10383] = 0e0;
		v[10384] = 0e0;
		v[10331] = 0e0;
		v[10332] = 0e0;
		v[10333] = 0e0;
		v[10334] = v[1248];
		v[10335] = 0e0;
		v[10336] = v[1246];
		v[10337] = 0e0;
		v[10338] = 0e0;
		v[10339] = 0e0;
		v[10340] = 0e0;
		v[10341] = 0e0;
		v[10342] = 0e0;
		v[10343] = 0e0;
		v[10344] = 0e0;
		v[10345] = 0e0;
		v[10346] = 0e0;
		v[10347] = 0e0;
		v[10348] = 0e0;
		v[10295] = 0e0;
		v[10296] = 0e0;
		v[10297] = 0e0;
		v[10298] = v[1247];
		v[10299] = v[1246];
		v[10300] = 0e0;
		v[10301] = 0e0;
		v[10302] = 0e0;
		v[10303] = 0e0;
		v[10304] = 0e0;
		v[10305] = 0e0;
		v[10306] = 0e0;
		v[10307] = 0e0;
		v[10308] = 0e0;
		v[10309] = 0e0;
		v[10310] = 0e0;
		v[10311] = 0e0;
		v[10312] = 0e0;
		v[10259] = 0e0;
		v[10260] = 0e0;
		v[10261] = 0e0;
		v[10262] = 0e0;
		v[10263] = 0e0;
		v[10264] = 0e0;
		v[10265] = 0e0;
		v[10266] = 0e0;
		v[10267] = 0e0;
		v[10268] = 0e0;
		v[10269] = v[1244];
		v[10270] = v[1243];
		v[10271] = 0e0;
		v[10272] = 0e0;
		v[10273] = 0e0;
		v[10274] = 0e0;
		v[10275] = 0e0;
		v[10276] = 0e0;
		v[10223] = 0e0;
		v[10224] = 0e0;
		v[10225] = 0e0;
		v[10226] = 0e0;
		v[10227] = 0e0;
		v[10228] = 0e0;
		v[10229] = 0e0;
		v[10230] = 0e0;
		v[10231] = 0e0;
		v[10232] = v[1244];
		v[10233] = 0e0;
		v[10234] = v[1242];
		v[10235] = 0e0;
		v[10236] = 0e0;
		v[10237] = 0e0;
		v[10238] = 0e0;
		v[10239] = 0e0;
		v[10240] = 0e0;
		v[10187] = 0e0;
		v[10188] = 0e0;
		v[10189] = 0e0;
		v[10190] = 0e0;
		v[10191] = 0e0;
		v[10192] = 0e0;
		v[10193] = 0e0;
		v[10194] = 0e0;
		v[10195] = 0e0;
		v[10196] = v[1243];
		v[10197] = v[1242];
		v[10198] = 0e0;
		v[10199] = 0e0;
		v[10200] = 0e0;
		v[10201] = 0e0;
		v[10202] = 0e0;
		v[10203] = 0e0;
		v[10204] = 0e0;
		v[10115] = 0e0;
		v[10116] = 0e0;
		v[10117] = 0e0;
		v[10118] = 0e0;
		v[10119] = 0e0;
		v[10120] = 0e0;
		v[10121] = 0e0;
		v[10122] = 0e0;
		v[10123] = 0e0;
		v[10124] = 0e0;
		v[10125] = 0e0;
		v[10126] = 0e0;
		v[10127] = 0e0;
		v[10128] = 0e0;
		v[10129] = 0e0;
		v[10130] = 0e0;
		v[10131] = v[1238];
		v[10132] = v[1237];
		v[10097] = 0e0;
		v[10098] = 0e0;
		v[10099] = 0e0;
		v[10100] = 0e0;
		v[10101] = 0e0;
		v[10102] = 0e0;
		v[10103] = 0e0;
		v[10104] = 0e0;
		v[10105] = 0e0;
		v[10106] = 0e0;
		v[10107] = 0e0;
		v[10108] = 0e0;
		v[10109] = 0e0;
		v[10110] = 0e0;
		v[10111] = 0e0;
		v[10112] = v[1238];
		v[10113] = 0e0;
		v[10114] = v[1236];
		v[10079] = 0e0;
		v[10080] = 0e0;
		v[10081] = 0e0;
		v[10082] = 0e0;
		v[10083] = 0e0;
		v[10084] = 0e0;
		v[10085] = 0e0;
		v[10086] = 0e0;
		v[10087] = 0e0;
		v[10088] = 0e0;
		v[10089] = 0e0;
		v[10090] = 0e0;
		v[10091] = 0e0;
		v[10092] = 0e0;
		v[10093] = 0e0;
		v[10094] = v[1237];
		v[10095] = v[1236];
		v[10096] = 0e0;
		v[10403] = v[1492];
		v[10404] = v[1480];
		v[10405] = v[1462];
		v[10406] = v[10366 + i1149] / 2e0 + v[1561] - v[1563] + 2e0*v[1331] * v[1583] + (-v[1560] + v[1585])*v[457]
			+ v[1582] * v[458] + v[1580] * v[461];
		v[10407] = v[10330 + i1149] / 2e0 - v[1562] + v[1566] + 2e0*v[1328] * v[1586] + v[1582] * v[459] + v[1581] * v[461]
			+ v[460] * (-v[1560] - v[1569] + v[6539]);
		v[10408] = v[10294 + i1149] / 2e0 + v[1565] - v[1567] + 2e0*v[1325] * v[1589] + v[1581] * v[458] + v[1580] * v[459] + (
			-v[1569] + v[1585])*v[462];
		v[10409] = v[1494];
		v[10410] = v[1482];
		v[10411] = v[1464];
		v[10412] = v[10258 + i1149] / 2e0 + v[1551] - v[1553] + 2e0*v[1340] * v[1590] + (-v[1550] + v[1592])*v[463]
			+ v[1579] * v[464] + v[1577] * v[467];
		v[10413] = v[10222 + i1149] / 2e0 - v[1552] + v[1556] + 2e0*v[1337] * v[1593] + v[1579] * v[465] + v[1578] * v[467]
			+ v[466] * (-v[1550] - v[1559] + v[6540]);
		v[10414] = v[10186 + i1149] / 2e0 + v[1555] - v[1557] + 2e0*v[1334] * v[1596] + v[1578] * v[464] + v[1577] * v[465] + (
			-v[1559] + v[1592])*v[468];
		v[10415] = -v[1395];
		v[10416] = -v[1400];
		v[10417] = -v[1409];
		v[10418] = v[10114 + i1149] / 2e0 + v[1510] - v[1512] + 2e0*v[1349] * v[1597] + (-v[1509] + v[1599])*v[469]
			+ v[1575] * v[470] + v[1573] * v[473];
		v[10419] = v[10096 + i1149] / 2e0 - v[1511] + v[1515] + 2e0*v[1346] * v[1600] + v[1575] * v[471] + v[1574] * v[473]
			+ v[472] * (-v[1509] - v[1518] + v[6541]);
		v[10420] = v[10078 + i1149] / 2e0 + v[1514] - v[1516] + 2e0*v[1343] * v[1603] + v[1574] * v[470] + v[1573] * v[471] + (
			-v[1518] + v[1599])*v[474];
		Rc[i1149 - 1] += v[1250] * v[1301] + v[1251] * v[1302] + v[1241] * v[1307] + v[1240] * v[1309] + v[9102 + i1149];
		for (i1284 = 1; i1284 <= 18; i1284++) {
			Kc[i1149 - 1][i1284 - 1] += v[10402 + i1284] + v[1571] * v[9120 + i1284] + v[1572] * v[9138 + i1284] + v[1531] * v[9156
				+ i1284] + v[1576] * v[9174 + i1284];
		};/* end for */
	};/* end for */
	v[1689] = 0e0;
	v[1690] = 0e0;
	v[1691] = 0e0;
	b1692 = b1049;
	if (b1692) {
		b1693 = b5;
		if (b1693) {
			b1694 = b1080;
			if (b1694) {
				v[1691] = 0e0;
				v[1690] = 0e0;
				v[1689] = 0e0;
			}
			else {
			};
		}
		else {
			b1695 = b1114;
			if (b1695) {
				v[1691] = 0e0;
				v[1690] = 0e0;
				v[1689] = 0e0;
			}
			else {
			};
		};
	}
	else {
	};
	v[6545] = v[1689] * v[33];
	v[6543] = v[1690] * v[33];
	v[6542] = v[1691] * v[33];
	v[1718] = v[1741] * v[6542];
	v[1719] = v[1743] * v[6542];
	v[2328] = v[6542] * v[6544];
	v[1721] = v[1747] * v[6542];
	v[6583] = v[1721] * v[223];
	v[6582] = v[1721] * v[224];
	v[1742] = v[1741] * v[6543];
	v[1744] = v[1743] * v[6543];
	v[2325] = v[6543] * v[6544];
	v[1748] = v[1747] * v[6543];
	v[6581] = v[1748] * v[223];
	v[6580] = v[1748] * v[224];
	v[1768] = v[1741] * v[6545];
	v[1769] = v[1743] * v[6545];
	v[2322] = v[6544] * v[6545];
	v[1771] = v[1747] * v[6545];
	v[6579] = v[1771] * v[223];
	v[6578] = v[1771] * v[224];
	v[1772] = v[1691] * v[32];
	v[3165] = -(v[1772] * v[391]);
	v[1773] = v[1690] * v[32];
	v[3166] = -(v[1773] * v[390]);
	v[3168] = v[3165] + v[3166];
	v[1774] = v[1689] * v[32];
	v[3161] = -(v[1774] * v[389]);
	v[3169] = v[3161] + v[3165];
	v[3167] = v[3161] + v[3166];
	v[1775] = 0e0;
	v[1776] = 0e0;
	v[1777] = 0e0;
	v[1778] = 0e0;
	v[1779] = 0e0;
	v[1780] = 0e0;
	v[1781] = 0e0;
	b1782 = b6;
	if (b1782) {
		b1783 = b1040;
		if (b1783) {
			v[1781] = 0e0;
			v[1780] = 0e0;
			v[1784] = 0e0;
			v[1779] = 0e0;
		}
		else {
			v[1784] = 0e0;
		};
		v[6725] = v[1784] * v[5754];
		v[1777] = 0e0;
		v[1776] = 0e0;
		v[1775] = 0e0;
		v[1778] = (v[6725] * Power(v[401], v[3585])) / v[1785];
	}
	else {
		v[1786] = 0e0;
		b1787 = b1049;
		if (b1787) {
			v[1788] = 0e0;
			v[1789] = 0e0;
			v[1790] = 0e0;
			b1791 = b1064;
			if (b1791) {
				v[1790] = 0e0;
				v[1789] = 0e0;
				v[1788] = 0e0;
			}
			else {
			};
			v[1804] = (v[1015] * v[1788] + v[1021] * v[1789] + v[1023] * v[1790])*v[6726];
			b1792 = b1051;
			if (b1792) {
				v[1781] = v[1059] * v[1790];
				v[1780] = v[1059] * v[1789];
				v[1779] = v[1059] * v[1788];
				v[1786] = (v[1804] * v[3584] * Power(v[1057], v[3585])) / v[1799];
			}
			else {
				v[1781] = v[1063] * v[1790];
				v[1780] = v[1063] * v[1789];
				v[1779] = v[1063] * v[1788];
				v[1778] = -((v[1804] * v[3801] * Power(v[401], v[3606])) / v[1805]);
			};
		}
		else {
		};
		b1806 = b1049;
		if (b1806) {
			b1807 = b1051;
			if (b1807) {
				v[1777] = 0e0;
				v[1776] = 0e0;
				v[1775] = 0e0;
				v[1778] = v[1778] - v[1786];
			}
			else {
			};
		}
		else {
		};
	};
	v[2143] = v[1778];
	v[6702] = v[1779] * v[389];
	v[6645] = v[1014] * v[1779];
	v[1996] = -(v[1779] * v[6499]);
	v[3339] = -(v[1018] * v[1780]);
	v[6720] = v[1781] * v[391];
	v[1840] = v[1781] * v[6546];
	v[1838] = v[390] * v[6720];
	v[6548] = v[1838] - v[1996] - v[3339];
	v[3342] = -v[1838] + v[3339];
	v[1775] = v[1775] + v[1781] * v[3552];
	v[1776] = v[1776] + v[1781] * v[3530];
	v[1815] = v[1781] * v[3506];
	v[1777] = v[1777] + v[1781] * v[3446];
	v[1775] = v[1775] + v[1780] * v[3425];
	v[1876] = v[1780] * v[3401];
	v[1776] = v[1776] + v[1780] * v[3362];
	v[1887] = v[1780] * v[390];
	v[6706] = -(v[1887] * v[391]) - v[1779] * v[6546];
	v[6549] = -(v[1022] * v[1781]) + v[6706];
	v[3340] = -(v[1887] * v[389]);
	v[6547] = v[1840] - v[3340] + v[6645];
	v[6644] = v[3340] + v[6547];
	v[1777] = v[1777] + v[1780] * v[3337];
	v[1943] = v[1779] * v[3301];
	v[1944] = v[1780] * v[286] + v[1779] * v[287];
	v[1945] = v[1780] * v[289] + v[1779] * v[290];
	v[3047] = v[1944] * v[223] + v[1945] * v[224];
	v[1775] = v[1775] + v[1779] * v[3262];
	v[1776] = v[1776] + v[1779] * v[3240];
	v[1966] = v[1887] + v[6702];
	v[1777] = v[1777] + v[1779] * v[3217];
	v[1969] = v[1781] * v[1816] + v[1780] * v[1889] + v[1779] * v[1968] - (v[1682] * v[1689] + v[1683] * v[1690]
		+ v[1684] * v[1691])*v[33];
	v[2298] = v[1969] * v[2522];
	v[2268] = v[1969] * v[4516];
	v[1971] = v[1781] * v[1818] + v[1780] * v[1891] + v[1779] * v[1970] - (v[1678] * v[1689] + v[1679] * v[1690]
		+ v[1680] * v[1691])*v[33];
	v[2281] = v[1971] * v[4514];
	v[6558] = v[2281] + v[1969] * v[2521];
	v[6557] = v[2268] + v[1971] * v[2519];
	v[1973] = v[1781] * v[1820] + v[1780] * v[1893] + v[1779] * v[1972] - (v[1674] * v[1689] + v[1675] * v[1690]
		+ v[1676] * v[1691])*v[33];
	v[2296] = v[1973] * v[4511];
	v[6919] = v[2296] + v[2298];
	v[6559] = v[2296] + v[1971] * v[2520];
	v[6907] = v[2298] + v[6559];
	v[2283] = v[1973] * v[6473];
	v[6915] = v[2283] + v[6558];
	v[6909] = v[2281] + v[2283];
	v[2271] = v[1973] * v[6476];
	v[6918] = v[2271] + v[6557];
	v[6908] = v[2268] + v[2271];
	v[1975] = v[1781] * v[1822] + v[1780] * v[1895] + v[1779] * v[1974] - (v[1658] * v[1689] + v[1659] * v[1690]
		+ v[1660] * v[1691])*v[33];
	v[2411] = v[1975] * v[2514];
	v[2391] = v[1975] * v[4509];
	v[1977] = v[1781] * v[1824] + v[1780] * v[1897] + v[1779] * v[1976] - (v[1654] * v[1689] + v[1655] * v[1690]
		+ v[1656] * v[1691])*v[33];
	v[2401] = v[1977] * v[4507];
	v[6585] = v[2401] + v[1975] * v[2513];
	v[6584] = v[2391] + v[1977] * v[2511];
	v[1979] = v[1781] * v[1826] + v[1780] * v[1899] + v[1779] * v[1978] - (v[1650] * v[1689] + v[1651] * v[1690]
		+ v[1652] * v[1691])*v[33];
	v[2409] = v[1979] * v[4504];
	v[6905] = v[2409] + v[2411];
	v[6586] = v[2409] + v[1977] * v[2512];
	v[6893] = v[2411] + v[6586];
	v[2403] = v[1979] * v[6465];
	v[6901] = v[2403] + v[6585];
	v[6895] = v[2401] + v[2403];
	v[2394] = v[1979] * v[6468];
	v[6904] = v[2394] + v[6584];
	v[6894] = v[2391] + v[2394];
	v[1981] = v[1781] * v[1828] + v[1780] * v[1901] + v[1779] * v[1980] - (v[1634] * v[1689] + v[1635] * v[1690]
		+ v[1636] * v[1691])*v[33];
	v[2438] = v[1981] * v[2506];
	v[2418] = v[1981] * v[4502];
	v[1983] = v[1781] * v[1830] + v[1780] * v[1903] + v[1779] * v[1982] - (v[1630] * v[1689] + v[1631] * v[1690]
		+ v[1632] * v[1691])*v[33];
	v[2428] = v[1983] * v[4500];
	v[6588] = v[2428] + v[1981] * v[2505];
	v[6587] = v[2418] + v[1983] * v[2503];
	v[1985] = v[1781] * v[1832] + v[1780] * v[1905] + v[1779] * v[1984] - (v[1626] * v[1689] + v[1627] * v[1690]
		+ v[1628] * v[1691])*v[33];
	v[2436] = v[1985] * v[4497];
	v[6891] = v[2436] + v[2438];
	v[6589] = v[2436] + v[1983] * v[2504];
	v[6879] = v[2438] + v[6589];
	v[2430] = v[1985] * v[6457];
	v[6887] = v[2430] + v[6588];
	v[6881] = v[2428] + v[2430];
	v[2421] = v[1985] * v[6460];
	v[6890] = v[2421] + v[6587];
	v[6880] = v[2418] + v[2421];
	v[1993] = -(v[352] * v[6547]);
	v[1994] = -(v[362] * v[6547]);
	v[1995] = -(v[372] * v[6547]);
	v[1997] = -(v[352] * v[6548]);
	v[1998] = -(v[362] * v[6548]);
	v[1999] = -(v[372] * v[6548]);
	v[2001] = v[352] * v[6549];
	v[2002] = v[362] * v[6549];
	v[2003] = v[372] * v[6549];
	v[2004] = -(v[326] * v[6549]);
	v[2005] = -(v[336] * v[6549]);
	v[2006] = -(v[346] * v[6549]);
	v[2007] = -(v[300] * v[6549]);
	v[2008] = -(v[310] * v[6549]);
	v[2009] = -(v[320] * v[6549]);
	v[2010] = v[326] * v[6548];
	v[2011] = v[336] * v[6548];
	v[2012] = v[346] * v[6548];
	v[2013] = v[300] * v[6548];
	v[2014] = v[310] * v[6548];
	v[2015] = v[320] * v[6548];
	v[2016] = v[326] * v[6547];
	v[2017] = v[336] * v[6547];
	v[2018] = v[346] * v[6547];
	v[2019] = v[300] * v[6547];
	v[2020] = v[310] * v[6547];
	v[2021] = v[320] * v[6547];
	v[2024] = 0e0;
	v[2025] = 0e0;
	v[2026] = 0e0;
	v[2027] = 0e0;
	v[2028] = 0e0;
	v[2029] = 0e0;
	v[2030] = 0e0;
	v[2031] = 0e0;
	v[2032] = 0e0;
	b2033 = b4;
	if (b2033) {
		v[1815] = v[1815] - v[1772] * v[978];
		v[1876] = v[1876] - v[1773] * v[977];
		v[2035] = v[1774] * v[2034] + v[3168] * v[389];
		v[2037] = v[1773] * v[2036] + v[3169] * v[390];
		v[2039] = v[1772] * v[2038] + v[3167] * v[391];
		v[1943] = v[1943] - v[1774] * v[976];
		v[1775] = v[1775] - v[1774] * v[6550] + v[3168] * v[976];
		v[1776] = v[1776] + v[1773] * v[6551] + v[3169] * v[977];
		v[1777] = v[1777] + v[1772] * v[6552] + v[3167] * v[978];
		v[2024] = v[1] * v[2035];
		v[2025] = v[2] * v[2035];
		v[2026] = v[2035] * v[3];
		v[2043] = -v[2035];
		v[2044] = -(v[2035] * v[272]);
		v[2045] = -(v[2035] * v[270]);
		v[2046] = -(v[2035] * v[269]);
		v[2047] = v[2035] * v[222];
		v[2048] = v[2035] * v[221];
		v[2049] = v[2047] * v[226];
		v[2050] = v[2047] * v[225];
		v[2051] = v[2048] * v[226];
		v[2052] = v[2048] * v[225];
		v[2027] = v[1] * v[2037];
		v[2028] = v[2] * v[2037];
		v[2029] = v[2037] * v[3];
		v[2053] = -v[2037];
		v[2054] = -(v[2037] * v[272]);
		v[2055] = -(v[2037] * v[270]);
		v[2056] = -(v[2037] * v[269]);
		v[2057] = v[2037] * v[222];
		v[2058] = v[2037] * v[221];
		v[2059] = v[2057] * v[226];
		v[2060] = v[2057] * v[225];
		v[2061] = v[2058] * v[226];
		v[2062] = v[2058] * v[225];
		v[2030] = v[1] * v[2039];
		v[2031] = v[2] * v[2039];
		v[2032] = v[2039] * v[3];
		v[2063] = -v[2039];
		v[2064] = -(v[2039] * v[272]);
		v[2065] = -(v[2039] * v[270]);
		v[2066] = -(v[2039] * v[269]);
		v[2067] = v[2039] * v[222];
		v[2068] = v[2039] * v[221];
		v[2069] = v[2067] * v[226];
		v[2070] = v[2067] * v[225];
		v[2071] = v[2068] * v[226];
		v[2072] = v[2068] * v[225];
	}
	else {
		v[2052] = 0e0;
		v[2051] = 0e0;
		v[2062] = 0e0;
		v[2061] = 0e0;
		v[2072] = 0e0;
		v[2071] = 0e0;
		v[2050] = 0e0;
		v[2049] = 0e0;
		v[2060] = 0e0;
		v[2059] = 0e0;
		v[2070] = 0e0;
		v[2069] = 0e0;
		v[2048] = 0e0;
		v[2058] = 0e0;
		v[2068] = 0e0;
		v[2047] = 0e0;
		v[2057] = 0e0;
		v[2067] = 0e0;
		v[2046] = 0e0;
		v[2045] = 0e0;
		v[2044] = 0e0;
		v[2056] = 0e0;
		v[2055] = 0e0;
		v[2054] = 0e0;
		v[2066] = 0e0;
		v[2065] = 0e0;
		v[2064] = 0e0;
		v[2043] = 0e0;
		v[2053] = 0e0;
		v[2063] = 0e0;
	};
	v[6693] = v[2024] / 2e0;
	v[6694] = v[2028] / 2e0;
	v[6695] = v[2032] / 2e0;
	b2073 = b941;
	if (b2073) {
		v[2107] = -(v[2025] * v[958]);
		v[2102] = v[2026] * v[958];
		v[2089] = v[2029] * v[958];
		v[2076] = -(v[2032] * v[6502]);
		v[2080] = v[2031] * v[958];
		v[2084] = v[2030] * v[958];
		v[2088] = v[2080] + v[2089];
		v[2093] = -(v[2028] * v[6502]);
		v[2097] = v[2027] * v[958];
		v[2101] = v[2084] + v[2102];
		v[2108] = v[2097] - v[2107];
		v[2110] = v[2031] * v[2078] + v[2030] * v[2082] + v[2029] * v[2086] + v[2027] * v[2095] + v[2026] * v[2099]
			+ v[2025] * v[2104] + v[2109] * v[6693] + v[2091] * v[6694] + v[2074] * v[6695];
		v[3115] = v[2076] + v[2093] - (4e0*v[2110]) / (v[2114] * v[2114]);
		v[6692] = 4e0*v[3115];
		v[3113] = -v[2076] + v[3115] - v[2024] * v[6502];
		v[6691] = 4e0*(v[2076] - v[2093] + v[3113]);
		v[2115] = v[2097] + v[2107] + v[2088] * v[6501] + v[2101] * v[6553] + 2e0*v[3113] * v[957];
		v[2117] = (-2e0*v[2084] + 2e0*v[2102] + v[2108] * v[955] + v[6691] * v[956] + v[2088] * v[957]) / 2e0;
		v[2118] = (2e0*v[2080] - 2e0*v[2089] + v[6692] * v[955] + v[2108] * v[956] + v[2101] * v[957]) / 2e0;
		v[6554] = v[2118] * v[943] + v[2117] * v[944] + v[2115] * v[945];
		v[3101] = v[6554] * v[954];
		v[3098] = v[6554] * v[948];
		v[2121] = v[3098] * v[953] + v[3101] / (Power(cos(v[2119]), 2)*sqrt(v[3102]));
		v[6555] = v[2121] / v[946];
		v[2122] = v[2115] * v[6500] + v[6555] * v[945];
		v[2124] = v[2117] * v[6500] + v[6555] * v[944];
		v[2125] = v[2118] * v[6500] + v[6555] * v[943];
		v[1775] = v[1775] - v[2122] * v[399] + v[2124] * v[400];
		v[1776] = v[1776] + v[2122] * v[398] - v[2125] * v[400];
		v[1777] = v[1777] - v[2124] * v[398] + v[2125] * v[399];
	}
	else {
	};
	v[1775] = v[1775] + v[1943] * v[6686];
	v[1775] = v[1775] + v[3047] * v[390];
	v[2138] = v[1775];
	v[1776] = v[1776] + v[3047] * v[389] + v[1876] * v[6683];
	v[2136] = v[1776];
	v[1777] = v[1777] + v[1966] * v[5461] + v[1815] * v[6680];
	v[2134] = v[1777];
	b2126 = b6;
	if (b2126) {
		v[2127] = v[1777] * v[375];
		v[2128] = v[1777] * v[388];
		v[1777] = 0e0;
		v[2129] = v[2127] + v[1776] * v[374];
		v[2130] = v[1776] * v[388];
		v[1776] = 0e0;
		v[2131] = v[2129] + v[1775] * v[373];
		v[2132] = v[1775] * v[388];
		v[1775] = 0e0;
		v[1778] = v[1778] + v[2131] * v[387];
	}
	else {
		b2133 = b418;
		if (b2133) {
			v[2135] = -v[2134];
			v[1777] = 0e0;
			v[2137] = -v[2136];
			v[1776] = 0e0;
			v[2139] = -v[2138];
			v[1775] = 0e0;
		}
		else {
			v[2135] = v[2134];
			v[2137] = v[2136];
			v[2139] = v[2138];
		};
		v[2128] = v[2135] * v[407];
		v[2130] = v[2137] * v[407];
		v[2142] = v[2139] * v[373] + v[2137] * v[374] + v[2135] * v[375];
		v[2132] = v[2139] * v[407];
		v[1778] = v[2143] + v[2142] * v[406];
	};
	v[6678] = -(v[1406] * v[1778]);
	v[6556] = v[1778] / v[401];
	v[2128] = v[2128] + v[375] * v[6556];
	v[2130] = v[2130] + v[374] * v[6556];
	v[2132] = v[2132] + v[373] * v[6556];
	v[2063] = v[2063] - v[2128];
	v[2064] = v[2064] - v[1027] * v[2128];
	v[2065] = v[2065] - v[2128] * v[274];
	v[2066] = v[2066] - v[1613] * v[2128];
	v[2053] = v[2053] - v[2130];
	v[2054] = v[2054] - v[1027] * v[2130];
	v[2055] = v[2055] - v[2130] * v[274];
	v[2056] = v[2056] - v[1613] * v[2130];
	v[2043] = v[2043] - v[2132];
	v[2044] = v[2044] - v[1027] * v[2132];
	v[2045] = v[2045] - v[2132] * v[274];
	v[2046] = v[2046] - v[1613] * v[2132];
	v[2165] = v[1969] * v[3215];
	v[2166] = v[1969] * v[3214];
	v[2167] = v[1969] * v[3213];
	v[2170] = v[1971] * v[3210];
	v[2171] = v[1969] * v[359] + v[1971] * v[361];
	v[6563] = v[2171] * v[243];
	v[2174] = v[1971] * v[3208];
	v[2175] = v[2165] + v[2170];
	v[2177] = v[1971] * v[3209] + v[6563] / v[350];
	v[2180] = v[1973] * v[3203];
	v[2181] = v[1969] * v[353] + v[1973] * v[361];
	v[6564] = v[2181] * v[248];
	v[2182] = v[1971] * v[353] + v[1973] * v[359];
	v[6565] = v[2182] * v[252];
	v[7373] = -(v[1969] * v[6350]) - v[1971] * v[6351] - v[1973] * v[6352] + v[6472] * v[6563] + v[349] * v[6564]
		+ v[348] * v[6565];
	v[2184] = v[1973] * v[3205] + v[6564] / v[350];
	v[2185] = v[2177] + v[2184];
	v[2187] = v[1973] * v[3204] + v[6565] / v[350];
	v[2188] = v[2166] + v[2187];
	v[2189] = v[1975] * v[3200];
	v[2190] = v[1975] * v[3199];
	v[2191] = v[1975] * v[3198];
	v[2194] = v[1977] * v[3195];
	v[2195] = v[1975] * v[333] + v[1977] * v[335];
	v[6593] = v[183] * v[2195];
	v[2198] = v[1977] * v[3193];
	v[2199] = v[2189] + v[2194];
	v[2201] = v[1977] * v[3194] + v[6593] / v[324];
	v[2204] = v[1979] * v[3188];
	v[2205] = v[1975] * v[327] + v[1979] * v[335];
	v[6594] = v[188] * v[2205];
	v[2206] = v[1977] * v[327] + v[1979] * v[333];
	v[6595] = v[192] * v[2206];
	v[7377] = -(v[1975] * v[6376]) - v[1977] * v[6377] - v[1979] * v[6378] + v[6464] * v[6593] + v[323] * v[6594]
		+ v[322] * v[6595];
	v[2208] = v[1979] * v[3190] + v[6594] / v[324];
	v[2209] = v[2201] + v[2208];
	v[2211] = v[1979] * v[3189] + v[6595] / v[324];
	v[2212] = v[2190] + v[2211];
	v[2213] = v[1981] * v[3185];
	v[2214] = v[1981] * v[3184];
	v[2215] = v[1981] * v[3183];
	v[2218] = v[1983] * v[3180];
	v[2219] = v[1981] * v[307] + v[1983] * v[309];
	v[6611] = v[164] * v[2219];
	v[2222] = v[1983] * v[3178];
	v[2223] = v[2213] + v[2218];
	v[2225] = v[1983] * v[3179] + v[6611] / v[298];
	v[2228] = v[1985] * v[3173];
	v[2229] = v[1981] * v[301] + v[1985] * v[309];
	v[6612] = v[169] * v[2229];
	v[2230] = v[1983] * v[301] + v[1985] * v[307];
	v[6613] = v[173] * v[2230];
	v[7381] = -(v[1981] * v[6399]) - v[1983] * v[6400] - v[1985] * v[6401] + v[6456] * v[6611] + v[297] * v[6612]
		+ v[296] * v[6613];
	v[2232] = v[1985] * v[3175] + v[6612] / v[298];
	v[2233] = v[2225] + v[2232];
	v[2235] = v[1985] * v[3174] + v[6613] / v[298];
	v[2236] = v[2214] + v[2235];
	v[2237] = -(v[1993] * v[982]);
	v[2238] = -(v[1993] * v[981]);
	v[2239] = -(v[1993] * v[983]);
	v[2240] = -(v[1994] * v[983]);
	v[2241] = -(v[1994] * v[982]);
	v[2242] = -(v[1994] * v[981]);
	v[2243] = -(v[1995] * v[983]);
	v[2244] = -(v[1995] * v[981]);
	v[2245] = -(v[1995] * v[982]);
	v[2246] = -(v[1997] * v[983]);
	v[2247] = -(v[1997] * v[982]);
	v[2248] = -(v[1997] * v[981]);
	v[2249] = -(v[1998] * v[983]);
	v[2250] = -(v[1998] * v[981]);
	v[2251] = -(v[1998] * v[982]);
	v[2252] = -(v[1999] * v[981]);
	v[2253] = -(v[1999] * v[983]);
	v[2254] = -(v[1999] * v[982]);
	v[2255] = -(v[2001] * v[983]);
	v[2256] = -(v[2001] * v[982]);
	v[2257] = -(v[2001] * v[981]);
	v[6835] = -2e0*v[2165] + 2e0*v[2170] + v[2238] * v[6626] + v[2237] * v[6627] + v[2248] * v[6628] + v[2246] * v[6629]
		+ v[2256] * v[6630] + v[2255] * v[6631] - v[2239] * v[694] - v[2247] * v[713] - v[2257] * v[731];
	v[2258] = -(v[2002] * v[982]);
	v[2259] = -(v[2002] * v[983]);
	v[2260] = -(v[2002] * v[981]);
	v[2261] = -(v[2003] * v[982]);
	v[2262] = -(v[2003] * v[983]);
	v[2263] = -(v[2003] * v[981]);
	v[6836] = -2e0*v[2177] + 2e0*v[2184] + v[2244] * v[6626] + v[2245] * v[6627] + v[2252] * v[6628] + v[2253] * v[6629]
		+ v[2261] * v[6630] + v[2262] * v[6631] - v[2243] * v[694] - v[2254] * v[713] - v[2263] * v[731];
	v[2064] = -(v[1613] * v[1718]) + v[2064] + v[1719] * v[451];
	v[2065] = v[1612] * v[1719] + v[2065];
	v[2066] = v[1027] * v[1718] + v[2066] + v[1719] * v[449];
	v[2276] = v[136] * v[2064] + v[135] * v[2065] + v[134] * v[2066] + v[361] * v[6918];
	v[2277] = v[133] * v[2064] + v[132] * v[2065] + v[131] * v[2066] + v[2182] * v[6476] + v[359] * v[6557];
	v[2278] = v[130] * v[2064] + v[129] * v[2065] + v[128] * v[2066] + v[353] * v[6908];
	v[2054] = -(v[1613] * v[1742]) + v[2054] + v[1744] * v[451];
	v[2055] = v[1612] * v[1744] + v[2055];
	v[2056] = v[1027] * v[1742] + v[2056] + v[1744] * v[449];
	v[2288] = v[136] * v[2054] + v[135] * v[2055] + v[134] * v[2056] + v[2181] * v[6473] + v[361] * v[6558];
	v[2289] = v[133] * v[2054] + v[132] * v[2055] + v[131] * v[2056] + v[359] * v[6915];
	v[2290] = v[130] * v[2054] + v[129] * v[2055] + v[128] * v[2056] + v[353] * v[6909];
	v[2044] = -(v[1613] * v[1768]) + v[2044] + v[1769] * v[451];
	v[2045] = v[1612] * v[1769] + v[2045];
	v[2046] = v[1027] * v[1768] + v[2046] + v[1769] * v[449];
	v[2304] = v[136] * v[2044] + v[135] * v[2045] + v[134] * v[2046] + v[2171] * v[2520] + v[361] * v[6919];
	v[2305] = v[133] * v[2044] + v[132] * v[2045] + v[131] * v[2046] + v[359] * v[6559];
	v[2306] = v[130] * v[2044] + v[129] * v[2045] + v[128] * v[2046] + v[353] * v[6907];
	v[2307] = v[2237] + v[2246] + v[2252] + v[2261] + v[2174] * v[6669] - v[2185] * v[701] - v[2175] * v[704];
	v[2308] = -v[2241] - v[2244] - v[2249] - v[2262] + v[2180] * v[6668] - v[2185] * v[698] + v[2188] * v[704];
	v[2309] = v[2305] * v[238] - v[2237] * v[691] + v[2241] * v[692] - v[2245] * v[693];
	v[2310] = -v[2238] - v[2250] - v[2255] - v[2258] + v[2167] * v[6667] - v[2175] * v[698] + v[2188] * v[701];
	v[2311] = v[2304] * v[238] - v[2238] * v[691] + v[2242] * v[692] - v[2244] * v[693];
	v[2312] = v[2290] * v[238] - v[2246] * v[691] + v[2249] * v[692] - v[2253] * v[693];
	v[2313] = v[2243] + v[2254];
	v[2314] = v[2288] * v[238] - v[2248] * v[691] + v[2250] * v[692] - v[2252] * v[693];
	v[2315] = v[2278] * v[238] - v[2255] * v[691] + v[2259] * v[692] - v[2262] * v[693];
	v[2316] = v[2277] * v[238] - v[2256] * v[691] + v[2258] * v[692] - v[2261] * v[693];
	v[2317] = v[2247] + v[2257];
	v[2318] = v[2240] + v[2260];
	v[11797] = 0e0;
	v[11798] = 0e0;
	v[11799] = 0e0;
	v[11800] = 0e0;
	v[11801] = 0e0;
	v[11802] = 0e0;
	v[11803] = 0e0;
	v[11804] = 0e0;
	v[11805] = 0e0;
	v[11806] = 0e0;
	v[11807] = 0e0;
	v[11808] = 0e0;
	v[11809] = 0e0;
	v[11810] = 0e0;
	v[11811] = 0e0;
	v[11812] = -0.5e0*v[2308] - v[2317];
	v[11813] = v[2307] / 2e0 - v[2318];
	v[11814] = -0.5e0*v[2310] - v[2313];
	v[2319] = 1e0 / (v[350] * v[350]);
	v[6924] = -(v[2319] * v[360]);
	v[6574] = -(v[2319] * v[359]);
	v[6573] = -(v[2319] * v[361]);
	v[6570] = -(v[2319] * v[353]);
	v[6569] = -(v[2319] * v[348]);
	v[6568] = -(v[1969] * v[2319]);
	v[6567] = -(v[2319] * v[349]);
	v[6566] = -(v[2319] * v[6472]);
	v[4437] = -(v[2319] * (v[243] * v[351] + v[6470]));
	v[4436] = -(v[2319] * (v[242] * v[351] + v[6471]));
	v[4435] = -(v[2319] * v[6560]);
	v[4434] = -(v[2319] * v[6474]);
	v[4433] = -(v[2319] * v[6561]);
	v[4432] = -(v[2319] * v[6562]);
	v[4431] = -(v[2319] * v[351]);
	v[6842] = v[353] * v[4431];
	v[4430] = -(v[2319] * v[354]);
	v[6841] = v[359] * v[4430];
	v[6840] = -(v[2319] * v[361] * v[363]);
	v[4428] = -(v[2319] * v[6563]);
	v[4427] = -(v[2319] * v[6564]);
	v[4426] = -(v[2319] * v[6565]);
	v[4425] = v[2171] * v[6566];
	v[4424] = v[2181] * v[6567];
	v[4423] = v[2182] * v[6569];
	v[4212] = v[1971] * v[6566];
	v[4211] = v[6475] * v[6568];
	v[4209] = v[1973] * v[4431];
	v[6805] = v[4209] + v[4212];
	v[7353] = v[4211] + v[6805];
	v[4195] = v[1973] * v[6567];
	v[4192] = v[1971] * v[4430];
	v[6804] = v[4192] + v[4195];
	v[4191] = v[6477] * v[6568];
	v[7352] = v[4191] + v[6804];
	v[4181] = v[1973] * v[6569];
	v[4179] = v[1971] * v[6924];
	v[4178] = v[363] * v[6568];
	v[6802] = v[4178] + v[4181];
	v[7351] = v[4179] + v[6802];
	v[4084] = v[250] * v[6570];
	v[4081] = v[245] * v[6570];
	v[4078] = -(v[2319] * v[6571]);
	v[4077] = -(v[2319] * v[6572]);
	v[4072] = v[248] * v[6573];
	v[4071] = v[247] * v[6574];
	v[5994] = v[4071] + v[4081];
	v[5986] = v[4072] + v[5994];
	v[5972] = v[4071] + v[4072];
	v[4069] = v[241] * v[6570];
	v[5991] = v[4069] + v[4077] + v[4078];
	v[5984] = -v[4078] + v[5991];
	v[5974] = -v[4077] + v[5991];
	v[4065] = v[253] * v[6573];
	v[5997] = v[4065] + v[4084];
	v[4064] = v[252] * v[6574];
	v[5980] = v[4064] + v[4065];
	v[5976] = v[4064] + v[5997];
	v[3212] = v[359] * v[4432] + v[353] * v[4433] + v[6840];
	v[3207] = v[361] * v[4434] + v[353] * v[4435] + v[6841];
	v[3202] = v[359] * v[4436] + v[361] * v[4437] + v[6842];
	v[2955] = v[243] * v[6566];
	v[2953] = v[248] * v[6567];
	v[2949] = v[252] * v[6569];
	v[4422] = v[2167] + v[2174] + v[2180] + v[2182] * v[2949] + v[2181] * v[2953] + v[2171] * v[2955] + v[1973] * v[3202]
		+ v[1971] * v[3207] + v[1969] * v[3212];
	v[2320] = v[2242] - v[2245] - v[2248] + v[2253] + v[2256] - v[2259] + v[2307] * v[470] - v[2308] * v[471]
		- v[2310] * v[473] + v[2276] * v[6575] + v[2289] * v[6576] + v[2306] * v[6577] + v[4422] * v[6648] + v[2313] * v[6837]
		+ v[2317] * v[6838] + v[2318] * v[6839] + v[2305] * v[699] + v[2304] * v[705] + v[2290] * v[709] + v[2288] * v[718]
		+ v[2278] * v[722] + v[2277] * v[726];
	v[2321] = v[2132] * v[223] + v[2322];
	v[2323] = v[2132] * v[224] - v[2322];
	v[2048] = v[2048] + v[2321];
	v[2047] = v[2047] + v[2323];
	v[2324] = v[2130] * v[223] + v[2325];
	v[2326] = v[2130] * v[224] - v[2325];
	v[2058] = v[2058] + v[2324];
	v[2057] = v[2057] + v[2326];
	v[2327] = v[2128] * v[223] + v[2328];
	v[2329] = v[2128] * v[224] - v[2328];
	v[2068] = v[2068] + v[2327];
	v[2067] = v[2067] + v[2329];
	v[2330] = v[2004] * v[986];
	v[2331] = v[2004] * v[985];
	v[2332] = v[2004] * v[984];
	v[2333] = v[2005] * v[985];
	v[2334] = v[2005] * v[986];
	v[2335] = v[2005] * v[984];
	v[2336] = v[2006] * v[985];
	v[2337] = v[2006] * v[986];
	v[2338] = v[2006] * v[984];
	v[2339] = v[2007] * v[989];
	v[2340] = v[2007] * v[988];
	v[2341] = v[2007] * v[987];
	v[2342] = v[2008] * v[988];
	v[2343] = v[2008] * v[989];
	v[2344] = v[2008] * v[987];
	v[2345] = v[2009] * v[988];
	v[2346] = v[2009] * v[989];
	v[2347] = v[2009] * v[987];
	v[2348] = v[2010] * v[986];
	v[2349] = v[2010] * v[985];
	v[2350] = v[2010] * v[984];
	v[2351] = v[2011] * v[986];
	v[2352] = v[2011] * v[984];
	v[2353] = v[2011] * v[985];
	v[2354] = v[2012] * v[984];
	v[2355] = v[2012] * v[986];
	v[2356] = v[2012] * v[985];
	v[2357] = v[2013] * v[989];
	v[2358] = v[2013] * v[988];
	v[2359] = v[2013] * v[987];
	v[2360] = v[2014] * v[989];
	v[2361] = v[2014] * v[987];
	v[2362] = v[2014] * v[988];
	v[2363] = v[2015] * v[987];
	v[2364] = v[2015] * v[989];
	v[2365] = v[2015] * v[988];
	v[2366] = v[2016] * v[985];
	v[2367] = v[2016] * v[984];
	v[2368] = v[2016] * v[986];
	v[6857] = -2e0*v[2189] + 2e0*v[2194] - v[2368] * v[523] - v[2349] * v[542] - v[2332] * v[560] + v[2331] * v[6632]
		+ v[2330] * v[6633] + v[2350] * v[6634] + v[2348] * v[6635] + v[2367] * v[6636] + v[2366] * v[6637];
	v[2369] = v[2017] * v[986];
	v[2370] = v[2017] * v[985];
	v[2371] = v[2017] * v[984];
	v[2372] = v[2018] * v[986];
	v[2373] = v[2018] * v[984];
	v[2374] = v[2018] * v[985];
	v[6858] = -2e0*v[2201] + 2e0*v[2208] - v[2372] * v[523] - v[2356] * v[542] - v[2338] * v[560] + v[2336] * v[6632]
		+ v[2337] * v[6633] + v[2354] * v[6634] + v[2355] * v[6635] + v[2373] * v[6636] + v[2374] * v[6637];
	v[2375] = v[2019] * v[988];
	v[2376] = v[2019] * v[987];
	v[2377] = v[2019] * v[989];
	v[6871] = -2e0*v[2213] + 2e0*v[2218] - v[2377] * v[478] - v[2358] * v[497] - v[2341] * v[515] + v[2340] * v[6638]
		+ v[2339] * v[6639] + v[2359] * v[6640] + v[2357] * v[6641] + v[2376] * v[6642] + v[2375] * v[6643];
	v[2378] = v[2020] * v[989];
	v[2379] = v[2020] * v[988];
	v[2380] = v[2020] * v[987];
	v[2381] = v[2021] * v[989];
	v[2382] = v[2021] * v[987];
	v[2383] = v[2021] * v[988];
	v[6872] = -2e0*v[2225] + 2e0*v[2232] - v[2381] * v[478] - v[2365] * v[497] - v[2347] * v[515] + v[2345] * v[6638]
		+ v[2346] * v[6639] + v[2363] * v[6640] + v[2364] * v[6641] + v[2382] * v[6642] + v[2383] * v[6643];
	v[2049] = v[2049] + v[228] * v[2323] + v[1610] * v[6578];
	v[2050] = v[2050] + v[227] * v[2323] - v[1024] * v[6578];
	v[2051] = v[2051] + v[228] * v[2321] + v[1610] * v[6579];
	v[2052] = v[2052] + v[227] * v[2321] - v[1024] * v[6579];
	v[2059] = v[2059] + v[228] * v[2326] + v[1610] * v[6580];
	v[2060] = v[2060] + v[227] * v[2326] - v[1024] * v[6580];
	v[2061] = v[2061] + v[228] * v[2324] + v[1610] * v[6581];
	v[2062] = v[2062] + v[227] * v[2324] - v[1024] * v[6581];
	v[2069] = v[2069] + v[228] * v[2329] + v[1610] * v[6582];
	v[2070] = v[2070] + v[227] * v[2329] - v[1024] * v[6582];
	v[2071] = v[2071] + v[228] * v[2327] + v[1610] * v[6583];
	v[2072] = v[2072] + v[227] * v[2327] - v[1024] * v[6583];
	v[2396] = v[335] * v[6904] + v[2070] * v[97] + v[2069] * v[98];
	v[2397] = v[2206] * v[6468] + v[333] * v[6584] + v[2070] * v[94] + v[2069] * v[95];
	v[2398] = v[327] * v[6894] + v[2070] * v[91] + v[2069] * v[92];
	v[2405] = v[2205] * v[6465] + v[335] * v[6585] + v[2060] * v[97] + v[2059] * v[98];
	v[2406] = v[333] * v[6901] + v[2060] * v[94] + v[2059] * v[95];
	v[2407] = v[327] * v[6895] + v[2060] * v[91] + v[2059] * v[92];
	v[2414] = v[2195] * v[2512] + v[335] * v[6905] + v[2050] * v[97] + v[2049] * v[98];
	v[2415] = v[333] * v[6586] + v[2050] * v[94] + v[2049] * v[95];
	v[2416] = v[327] * v[6893] + v[2050] * v[91] + v[2049] * v[92];
	v[2423] = v[309] * v[6890] + v[2072] * v[88] + v[2071] * v[89];
	v[2424] = v[2230] * v[6460] + v[307] * v[6587] + v[2072] * v[85] + v[2071] * v[86];
	v[2425] = v[301] * v[6880] + v[2072] * v[82] + v[2071] * v[83];
	v[2432] = v[2229] * v[6457] + v[309] * v[6588] + v[2062] * v[88] + v[2061] * v[89];
	v[2433] = v[307] * v[6887] + v[2062] * v[85] + v[2061] * v[86];
	v[2434] = v[301] * v[6881] + v[2062] * v[82] + v[2061] * v[83];
	v[2441] = v[2219] * v[2504] + v[309] * v[6891] + v[2052] * v[88] + v[2051] * v[89];
	v[2442] = v[307] * v[6589] + v[2052] * v[85] + v[2051] * v[86];
	v[2443] = v[301] * v[6879] + v[2052] * v[82] + v[2051] * v[83];
	v[2444] = v[2336] + v[2348] + v[2354] + v[2366] - v[2209] * v[530] - v[2199] * v[533] + v[2198] * v[6659];
	v[2445] = -v[2337] - v[2351] - v[2370] - v[2373] - v[2209] * v[527] + v[2212] * v[533] + v[2204] * v[6658];
	v[2446] = v[178] * v[2415] - v[2366] * v[520] + v[2370] * v[521] - v[2374] * v[522];
	v[2447] = -v[2330] - v[2333] - v[2352] - v[2367] - v[2199] * v[527] + v[2212] * v[530] + v[2191] * v[6657];
	v[2448] = v[178] * v[2414] - v[2367] * v[520] + v[2371] * v[521] - v[2373] * v[522];
	v[2449] = v[178] * v[2407] - v[2348] * v[520] + v[2351] * v[521] - v[2355] * v[522];
	v[2450] = v[2356] + v[2372];
	v[2451] = v[178] * v[2405] - v[2350] * v[520] + v[2352] * v[521] - v[2354] * v[522];
	v[2452] = v[178] * v[2398] - v[2330] * v[520] + v[2334] * v[521] - v[2337] * v[522];
	v[2453] = v[178] * v[2397] - v[2331] * v[520] + v[2333] * v[521] - v[2336] * v[522];
	v[2454] = v[2332] + v[2349];
	v[2455] = v[2335] + v[2369];
	v[11815] = 0e0;
	v[11816] = 0e0;
	v[11817] = 0e0;
	v[11818] = 0e0;
	v[11819] = 0e0;
	v[11820] = 0e0;
	v[11821] = 0e0;
	v[11822] = 0e0;
	v[11823] = 0e0;
	v[11824] = -0.5e0*v[2445] - v[2454];
	v[11825] = v[2444] / 2e0 - v[2455];
	v[11826] = -0.5e0*v[2447] - v[2450];
	v[11827] = 0e0;
	v[11828] = 0e0;
	v[11829] = 0e0;
	v[11830] = 0e0;
	v[11831] = 0e0;
	v[11832] = 0e0;
	v[2456] = 1e0 / (v[324] * v[324]);
	v[6926] = -(v[2456] * v[334]);
	v[6604] = -(v[2456] * v[333]);
	v[6603] = -(v[2456] * v[335]);
	v[6600] = -(v[2456] * v[327]);
	v[6599] = -(v[2456] * v[322]);
	v[6598] = -(v[1975] * v[2456]);
	v[6597] = -(v[2456] * v[323]);
	v[6596] = -(v[2456] * v[6464]);
	v[4465] = -(v[2456] * (v[183] * v[325] + v[6462]));
	v[4464] = -(v[2456] * (v[182] * v[325] + v[6463]));
	v[4463] = -(v[2456] * v[6590]);
	v[4462] = -(v[2456] * v[6466]);
	v[4461] = -(v[2456] * v[6591]);
	v[4460] = -(v[2456] * v[6592]);
	v[4459] = -(v[2456] * v[325]);
	v[6864] = v[327] * v[4459];
	v[4458] = -(v[2456] * v[328]);
	v[6863] = v[333] * v[4458];
	v[6862] = -(v[2456] * v[335] * v[337]);
	v[4456] = -(v[2456] * v[6593]);
	v[4455] = -(v[2456] * v[6594]);
	v[4454] = -(v[2456] * v[6595]);
	v[4453] = v[2195] * v[6596];
	v[4452] = v[2205] * v[6597];
	v[4451] = v[2206] * v[6599];
	v[4336] = v[1977] * v[6596];
	v[4335] = v[6467] * v[6598];
	v[4333] = v[1979] * v[4459];
	v[7360] = v[335] * (v[4333] + v[4335]) + v[4453];
	v[6825] = v[4333] + v[4336];
	v[7361] = v[4335] + v[6825];
	v[4326] = v[1979] * v[6597];
	v[4323] = v[1977] * v[4458];
	v[7359] = v[4323] + v[4326];
	v[4322] = v[6469] * v[6598];
	v[6824] = v[4322] + v[4323];
	v[7358] = v[4326] + v[6824];
	v[7357] = v[4452] + v[335] * v[6824];
	v[4315] = v[1979] * v[6599];
	v[4313] = v[1977] * v[6926];
	v[4312] = v[337] * v[6598];
	v[7356] = v[4312] + v[4315];
	v[6823] = v[4312] + v[4313];
	v[7355] = v[4451] + v[333] * v[6823];
	v[7354] = v[4315] + v[6823];
	v[4114] = v[190] * v[6600];
	v[4111] = v[185] * v[6600];
	v[4108] = -(v[2456] * v[6601]);
	v[4107] = -(v[2456] * v[6602]);
	v[4102] = v[188] * v[6603];
	v[4101] = v[187] * v[6604];
	v[6024] = v[4101] + v[4111];
	v[6016] = v[4102] + v[6024];
	v[6002] = v[4101] + v[4102];
	v[4099] = v[181] * v[6600];
	v[6021] = v[4099] + v[4107] + v[4108];
	v[6014] = -v[4108] + v[6021];
	v[6004] = -v[4107] + v[6021];
	v[4095] = v[193] * v[6603];
	v[6027] = v[4095] + v[4114];
	v[4094] = v[192] * v[6604];
	v[6010] = v[4094] + v[4095];
	v[6006] = v[4094] + v[6027];
	v[3197] = v[333] * v[4460] + v[327] * v[4461] + v[6862];
	v[3192] = v[335] * v[4462] + v[327] * v[4463] + v[6863];
	v[3187] = v[333] * v[4464] + v[335] * v[4465] + v[6864];
	v[2927] = v[183] * v[6596];
	v[2925] = v[188] * v[6597];
	v[2921] = v[192] * v[6599];
	v[4450] = v[2191] + v[2198] + v[2204] + v[2206] * v[2921] + v[2205] * v[2925] + v[2195] * v[2927] + v[1979] * v[3187]
		+ v[1977] * v[3192] + v[1975] * v[3197];
	v[2457] = v[2331] - v[2334] - v[2350] + v[2355] + v[2371] - v[2374] + v[2444] * v[464] - v[2445] * v[465]
		- v[2447] * v[467] + v[2415] * v[528] + v[2414] * v[534] + v[2407] * v[538] + v[2405] * v[547] + v[2398] * v[551]
		+ v[2397] * v[555] + v[2396] * v[6605] + v[2406] * v[6606] + v[2416] * v[6607] + v[4450] * v[6647] + v[2450] * v[6859]
		+ v[2454] * v[6860] + v[2455] * v[6861];
	v[2458] = v[2345] + v[2357] + v[2363] + v[2375] - v[2233] * v[485] - v[2223] * v[488] + v[2222] * v[6655];
	v[2459] = -v[2346] - v[2360] - v[2379] - v[2382] - v[2233] * v[482] + v[2236] * v[488] + v[2228] * v[6654];
	v[2460] = v[159] * v[2442] - v[2375] * v[475] + v[2379] * v[476] - v[2383] * v[477];
	v[2461] = -v[2339] - v[2342] - v[2361] - v[2376] - v[2223] * v[482] + v[2236] * v[485] + v[2215] * v[6653];
	v[2462] = v[159] * v[2441] - v[2376] * v[475] + v[2380] * v[476] - v[2382] * v[477];
	v[2463] = v[159] * v[2434] - v[2357] * v[475] + v[2360] * v[476] - v[2364] * v[477];
	v[2464] = v[2365] + v[2381];
	v[2465] = v[159] * v[2432] - v[2359] * v[475] + v[2361] * v[476] - v[2363] * v[477];
	v[2466] = v[159] * v[2425] - v[2339] * v[475] + v[2343] * v[476] - v[2346] * v[477];
	v[2467] = v[159] * v[2424] - v[2340] * v[475] + v[2342] * v[476] - v[2345] * v[477];
	v[2468] = v[2341] + v[2358];
	v[2469] = v[2344] + v[2378];
	v[11833] = 0e0;
	v[11834] = 0e0;
	v[11835] = 0e0;
	v[11836] = -0.5e0*v[2459] - v[2468];
	v[11837] = v[2458] / 2e0 - v[2469];
	v[11838] = -0.5e0*v[2461] - v[2464];
	v[11839] = 0e0;
	v[11840] = 0e0;
	v[11841] = 0e0;
	v[11842] = 0e0;
	v[11843] = 0e0;
	v[11844] = 0e0;
	v[11845] = 0e0;
	v[11846] = 0e0;
	v[11847] = 0e0;
	v[11848] = 0e0;
	v[11849] = 0e0;
	v[11850] = 0e0;
	v[2470] = 1e0 / (v[298] * v[298]);
	v[6928] = -(v[2470] * v[308]);
	v[6622] = -(v[2470] * v[307]);
	v[6621] = -(v[2470] * v[309]);
	v[6618] = -(v[2470] * v[301]);
	v[6617] = -(v[2470] * v[296]);
	v[6616] = -(v[1981] * v[2470]);
	v[6615] = -(v[2470] * v[297]);
	v[6614] = -(v[2470] * v[6456]);
	v[4492] = -(v[2470] * (v[164] * v[299] + v[6454]));
	v[4491] = -(v[2470] * (v[163] * v[299] + v[6455]));
	v[4490] = -(v[2470] * v[6608]);
	v[4489] = -(v[2470] * v[6458]);
	v[4488] = -(v[2470] * v[6609]);
	v[4487] = -(v[2470] * v[6610]);
	v[4486] = -(v[2470] * v[299]);
	v[6878] = v[301] * v[4486];
	v[4485] = -(v[2470] * v[302]);
	v[6877] = v[307] * v[4485];
	v[6876] = -(v[2470] * v[309] * v[311]);
	v[4483] = -(v[2470] * v[6611]);
	v[4482] = -(v[2470] * v[6612]);
	v[4481] = -(v[2470] * v[6613]);
	v[4480] = v[2219] * v[6614];
	v[4479] = v[2229] * v[6615];
	v[4478] = v[2230] * v[6617];
	v[4366] = v[1983] * v[6614];
	v[4365] = v[6459] * v[6616];
	v[4363] = v[1985] * v[4486];
	v[7368] = v[309] * (v[4363] + v[4365]) + v[4480];
	v[6828] = v[4363] + v[4366];
	v[7369] = v[4365] + v[6828];
	v[4356] = v[1985] * v[6615];
	v[4353] = v[1983] * v[4485];
	v[7367] = v[4353] + v[4356];
	v[4352] = v[6461] * v[6616];
	v[6827] = v[4352] + v[4353];
	v[7366] = v[4356] + v[6827];
	v[7365] = v[4479] + v[309] * v[6827];
	v[4345] = v[1985] * v[6617];
	v[4343] = v[1983] * v[6928];
	v[4342] = v[311] * v[6616];
	v[7364] = v[4342] + v[4345];
	v[6826] = v[4342] + v[4343];
	v[7363] = v[4478] + v[307] * v[6826];
	v[7362] = v[4345] + v[6826];
	v[4144] = v[171] * v[6618];
	v[4141] = v[166] * v[6618];
	v[4138] = -(v[2470] * v[6619]);
	v[4137] = -(v[2470] * v[6620]);
	v[4132] = v[169] * v[6621];
	v[4131] = v[168] * v[6622];
	v[6054] = v[4131] + v[4141];
	v[6046] = v[4132] + v[6054];
	v[6032] = v[4131] + v[4132];
	v[4129] = v[162] * v[6618];
	v[6051] = v[4129] + v[4137] + v[4138];
	v[6044] = -v[4138] + v[6051];
	v[6034] = -v[4137] + v[6051];
	v[4125] = v[174] * v[6621];
	v[6057] = v[4125] + v[4144];
	v[4124] = v[173] * v[6622];
	v[6040] = v[4124] + v[4125];
	v[6036] = v[4124] + v[6057];
	v[3182] = v[307] * v[4487] + v[301] * v[4488] + v[6876];
	v[3177] = v[309] * v[4489] + v[301] * v[4490] + v[6877];
	v[3172] = v[307] * v[4491] + v[309] * v[4492] + v[6878];
	v[2899] = v[164] * v[6614];
	v[2897] = v[169] * v[6615];
	v[2893] = v[173] * v[6617];
	v[4477] = v[2215] + v[2222] + v[2228] + v[2230] * v[2893] + v[2229] * v[2897] + v[2219] * v[2899] + v[1985] * v[3172]
		+ v[1983] * v[3177] + v[1981] * v[3182];
	v[2471] = v[2340] - v[2343] - v[2359] + v[2364] + v[2380] - v[2383] + v[2458] * v[458] - v[2459] * v[459]
		- v[2461] * v[461] + v[2442] * v[483] + v[2441] * v[489] + v[2434] * v[493] + v[2432] * v[502] + v[2425] * v[506]
		+ v[2424] * v[510] + v[2423] * v[6623] + v[2433] * v[6624] + v[2443] * v[6625] + v[4477] * v[6646] + v[2464] * v[6873]
		+ v[2468] * v[6874] + v[2469] * v[6875];
	v[2472] = v[6835] / 2e0;
	v[2474] = -v[2166] + v[2187] + v[2260] * v[6575] + v[2251] * v[6576] + v[2240] * v[6577] + v[2241] * v[699]
		+ v[2242] * v[705] + v[2249] * v[709] + v[2250] * v[718] + v[2259] * v[722] + v[2258] * v[726];
	v[11455] = 0e0;
	v[11456] = 0e0;
	v[11457] = 0e0;
	v[11458] = 0e0;
	v[11459] = 0e0;
	v[11460] = 0e0;
	v[11461] = 0e0;
	v[11462] = 0e0;
	v[11463] = 0e0;
	v[11464] = 0e0;
	v[11465] = 0e0;
	v[11466] = 0e0;
	v[11467] = 0e0;
	v[11468] = 0e0;
	v[11469] = 0e0;
	v[11470] = -v[6835];
	v[11471] = 2e0*v[2474];
	v[11472] = -v[6836];
	v[2475] = (v[2289] * v[238] - v[2247] * v[691] + v[2251] * v[692] - v[2254] * v[693]) / 2e0;
	v[2476] = v[6836] / 2e0;
	v[7371] = v[2472] * v[469] - v[2474] * v[472] + v[2476] * v[474];
	v[2518] = -v[2475] + v[2476] * v[2634] + v[2474] * v[2637] + v[2472] * v[2639] - v[2320] * v[6432];
	v[4517] = v[2518] + (-(v[2306] * v[238]) + v[2239] * v[691] - v[2240] * v[692] + v[2243] * v[693]) / 2e0;
	v[2477] = (v[2276] * v[238] - v[2257] * v[691] + v[2260] * v[692] - v[2263] * v[693]) / 2e0;
	v[4515] = v[2475] - v[2477] + v[4517];
	v[4512] = -v[2477] + v[2518];
	v[2478] = v[2311] + v[2315];
	v[2479] = v[2314] + v[2316];
	v[2480] = v[2309] + v[2312];
	v[2481] = v[6857] / 2e0;
	v[2483] = -v[2190] + v[2211] + v[2370] * v[528] + v[2371] * v[534] + v[2351] * v[538] + v[2352] * v[547] + v[2334] * v[551]
		+ v[2333] * v[555] + v[2335] * v[6605] + v[2353] * v[6606] + v[2369] * v[6607];
	v[11437] = 0e0;
	v[11438] = 0e0;
	v[11439] = 0e0;
	v[11440] = 0e0;
	v[11441] = 0e0;
	v[11442] = 0e0;
	v[11443] = 0e0;
	v[11444] = 0e0;
	v[11445] = 0e0;
	v[11446] = -v[6857];
	v[11447] = 2e0*v[2483];
	v[11448] = -v[6858];
	v[11449] = 0e0;
	v[11450] = 0e0;
	v[11451] = 0e0;
	v[11452] = 0e0;
	v[11453] = 0e0;
	v[11454] = 0e0;
	v[2484] = (v[178] * v[2406] - v[2349] * v[520] + v[2353] * v[521] - v[2356] * v[522]) / 2e0;
	v[2485] = v[6858] / 2e0;
	v[7375] = v[2481] * v[463] - v[2483] * v[466] + v[2485] * v[468];
	v[2510] = -v[2484] + v[2485] * v[2608] + v[2483] * v[2611] + v[2481] * v[2613] - v[2457] * v[6430];
	v[4510] = v[2510] + (-(v[178] * v[2416]) + v[2368] * v[520] - v[2369] * v[521] + v[2372] * v[522]) / 2e0;
	v[2486] = (v[178] * v[2396] - v[2332] * v[520] + v[2335] * v[521] - v[2338] * v[522]) / 2e0;
	v[4508] = v[2484] - v[2486] + v[4510];
	v[4505] = -v[2486] + v[2510];
	v[2487] = v[2448] + v[2452];
	v[2488] = v[2451] + v[2453];
	v[2489] = v[2446] + v[2449];
	v[2490] = v[6871] / 2e0;
	v[2492] = -v[2214] + v[2235] + v[2379] * v[483] + v[2380] * v[489] + v[2360] * v[493] + v[2361] * v[502] + v[2343] * v[506]
		+ v[2342] * v[510] + v[2344] * v[6623] + v[2362] * v[6624] + v[2378] * v[6625];
	v[11419] = 0e0;
	v[11420] = 0e0;
	v[11421] = 0e0;
	v[11422] = -v[6871];
	v[11423] = 2e0*v[2492];
	v[11424] = -v[6872];
	v[11425] = 0e0;
	v[11426] = 0e0;
	v[11427] = 0e0;
	v[11428] = 0e0;
	v[11429] = 0e0;
	v[11430] = 0e0;
	v[11431] = 0e0;
	v[11432] = 0e0;
	v[11433] = 0e0;
	v[11434] = 0e0;
	v[11435] = 0e0;
	v[11436] = 0e0;
	v[2493] = (v[159] * v[2433] - v[2358] * v[475] + v[2362] * v[476] - v[2365] * v[477]) / 2e0;
	v[2494] = v[6872] / 2e0;
	v[7379] = v[2490] * v[457] - v[2492] * v[460] + v[2494] * v[462];
	v[2502] = -v[2493] + v[2494] * v[2582] + v[2492] * v[2585] + v[2490] * v[2587] - v[2471] * v[6428];
	v[4503] = v[2502] + (-(v[159] * v[2443]) + v[2377] * v[475] - v[2378] * v[476] + v[2381] * v[477]) / 2e0;
	v[2495] = (v[159] * v[2423] - v[2341] * v[475] + v[2344] * v[476] - v[2347] * v[477]) / 2e0;
	v[4501] = v[2493] - v[2495] + v[4503];
	v[4498] = -v[2495] + v[2502];
	v[2496] = v[2462] + v[2466];
	v[2497] = v[2465] + v[2467];
	v[2498] = v[2460] + v[2463];
	v[10425] = -(v[1082] * v[1614]) - v[1083] * v[1615] - v[1084] * v[1616] + v[2048] + (v[1016] * v[1780] + (-
		(v[1614] * v[1689]) - v[1615] * v[1690] - v[1616] * v[1691])*v[33] + v[223] * v[6644])*v[7];
	v[10426] = -(v[1082] * v[1618]) - v[1083] * v[1619] - v[1084] * v[1620] + v[2058] + (v[1016] * v[1779] + (-
		(v[1618] * v[1689]) - v[1619] * v[1690] - v[1620] * v[1691])*v[33] - v[223] * v[3342])*v[7];
	v[10427] = -(v[1082] * v[1622]) - v[1083] * v[1623] - v[1084] * v[1624] + v[2068] + ((-(v[1622] * v[1689])
		- v[1623] * v[1690] - v[1624] * v[1691])*v[33] - v[223] * v[6549])*v[7];
	v[10428] = -(v[1082] * v[1626]) - v[1083] * v[1627] - v[1084] * v[1628] - v[2465] + v[2467] + v[4498] * v[457]
		+ v[2498] * v[458] + v[2496] * v[461] + v[2459] * v[6434] + 2e0*(v[2490] * v[6428] + v[2468] * v[6434]) + (v[2436]
			+ v[1983] * v[6457] + v[1981] * v[6460])*v[7];
	v[10429] = -(v[1082] * v[1630]) - v[1083] * v[1631] - v[1084] * v[1632] + v[2462] - v[2466] + v[2498] * v[459]
		+ v[4501] * v[460] + v[2497] * v[461] - v[2458] * v[6434] + 2e0*(-(v[2492] * v[6428]) + v[2469] * v[6434]) + (v[2428]
			+ v[1981] * v[2503] + v[1985] * v[2504])*v[7];
	v[10430] = -(v[1082] * v[1634]) - v[1083] * v[1635] - v[1084] * v[1636] - v[2460] + v[2463] + v[2497] * v[458]
		+ v[2496] * v[459] + v[4503] * v[462] + v[2461] * v[6434] + 2e0*(v[2494] * v[6428] + v[2464] * v[6434]) + (v[2418]
			+ v[1983] * v[2505] + v[1985] * v[2506])*v[7];
	v[10431] = -(v[1082] * v[1638]) - v[1083] * v[1639] - v[1084] * v[1640] + v[2047] + (v[1017] * v[1780] + (-
		(v[1638] * v[1689]) - v[1639] * v[1690] - v[1640] * v[1691])*v[33] + v[224] * v[6644])*v[7];
	v[10432] = -(v[1082] * v[1642]) - v[1083] * v[1643] - v[1084] * v[1644] + v[2057] + (v[1017] * v[1779] + (-
		(v[1642] * v[1689]) - v[1643] * v[1690] - v[1644] * v[1691])*v[33] - v[224] * v[3342])*v[7];
	v[10433] = -(v[1082] * v[1646]) - v[1083] * v[1647] - v[1084] * v[1648] + v[2067] + ((-(v[1646] * v[1689])
		- v[1647] * v[1690] - v[1648] * v[1691])*v[33] - v[224] * v[6549])*v[7];
	v[10434] = -(v[1082] * v[1650]) - v[1083] * v[1651] - v[1084] * v[1652] - v[2451] + v[2453] + v[4505] * v[463]
		+ v[2489] * v[464] + v[2487] * v[467] + v[2445] * v[6441] + 2e0*(v[2481] * v[6430] + v[2454] * v[6441]) + (v[2409]
			+ v[1977] * v[6465] + v[1975] * v[6468])*v[7];
	v[10435] = -(v[1082] * v[1654]) - v[1083] * v[1655] - v[1084] * v[1656] + v[2448] - v[2452] + v[2489] * v[465]
		+ v[4508] * v[466] + v[2488] * v[467] - v[2444] * v[6441] + 2e0*(-(v[2483] * v[6430]) + v[2455] * v[6441]) + (v[2401]
			+ v[1975] * v[2511] + v[1979] * v[2512])*v[7];
	v[10436] = -(v[1082] * v[1658]) - v[1083] * v[1659] - v[1084] * v[1660] - v[2446] + v[2449] + v[2488] * v[464]
		+ v[2487] * v[465] + v[4510] * v[468] + v[2447] * v[6441] + 2e0*(v[2485] * v[6430] + v[2450] * v[6441]) + (v[2391]
			+ v[1977] * v[2513] + v[1979] * v[2514])*v[7];
	v[10437] = -(v[1082] * v[1662]) - v[1083] * v[1663] - v[1084] * v[1664] + v[2043] + (-v[1840] + (-(v[1662] * v[1689])
		- v[1663] * v[1690] - v[1664] * v[1691])*v[33] + v[3340] - v[6645])*v[7];
	v[10438] = -(v[1082] * v[1666]) - v[1083] * v[1667] - v[1084] * v[1668] + v[2053] + (v[1996] + (-(v[1666] * v[1689])
		- v[1667] * v[1690] - v[1668] * v[1691])*v[33] + v[3342])*v[7];
	v[10439] = -(v[1082] * v[1670]) - v[1083] * v[1671] - v[1084] * v[1672] + v[2063] + ((-(v[1670] * v[1689])
		- v[1671] * v[1690] - v[1672] * v[1691])*v[33] + v[6549])*v[7];
	v[10440] = -(v[1082] * v[1674]) - v[1083] * v[1675] - v[1084] * v[1676] - v[2314] + v[2316] + v[4512] * v[469]
		+ v[2480] * v[470] + v[2478] * v[473] + v[2308] * v[6448] + 2e0*(v[2472] * v[6432] + v[2317] * v[6448]) + (v[2296]
			+ v[1971] * v[6473] + v[1969] * v[6476])*v[7];
	v[10441] = -(v[1082] * v[1678]) - v[1083] * v[1679] - v[1084] * v[1680] + v[2311] - v[2315] + v[2480] * v[471]
		+ v[4515] * v[472] + v[2479] * v[473] - v[2307] * v[6448] + 2e0*(-(v[2474] * v[6432]) + v[2318] * v[6448]) + (v[2281]
			+ v[1969] * v[2519] + v[1973] * v[2520])*v[7];
	v[10442] = -(v[1082] * v[1682]) - v[1083] * v[1683] - v[1084] * v[1684] - v[2309] + v[2312] + v[2479] * v[470]
		+ v[2478] * v[471] + v[4517] * v[474] + v[2310] * v[6448] + 2e0*(v[2476] * v[6432] + v[2313] * v[6448]) + (v[2268]
			+ v[1971] * v[2521] + v[1973] * v[2522])*v[7];
	for (i1687 = 1; i1687 <= 18; i1687++) {
		i6675 = (i1687 == 15 ? 1 : 0);
		i6674 = (i1687 == 14 ? 1 : 0);
		i6673 = (i1687 == 13 ? 1 : 0);
		i6665 = (i1687 == 1 ? 1 : 0);
		i6664 = (i1687 == 7 ? 1 : 0);
		i6663 = (i1687 == 2 ? 1 : 0);
		i6662 = (i1687 == 8 ? 1 : 0);
		i6661 = (i1687 == 3 ? 1 : 0);
		i6660 = (i1687 == 9 ? 1 : 0);
		v[2529] = v[9232 + i1687];
		v[2531] = v[9268 + i1687];
		v[2533] = v[9250 + i1687];
		v[2534] = v[10446 + i1687];
		v[2536] = v[9214 + i1687];
		v[2661] = -(v[2536] * v[6428]);
		v[2683] = v[2661] * v[6646];
		v[6670] = v[2683] * v[298];
		v[2590] = -0.5e0*v[2661];
		v[2537] = v[10518 + i1687];
		v[2539] = v[9304 + i1687];
		v[2541] = v[9340 + i1687];
		v[2543] = v[9322 + i1687];
		v[2544] = v[10554 + i1687];
		v[2546] = v[9286 + i1687];
		v[2699] = -(v[2546] * v[6430]);
		v[2721] = v[2699] * v[6647];
		v[6671] = v[2721] * v[324];
		v[2616] = -0.5e0*v[2699];
		v[2547] = v[10626 + i1687];
		v[2549] = v[9376 + i1687];
		v[2551] = v[9412 + i1687];
		v[2553] = v[9394 + i1687];
		v[2554] = v[10662 + i1687];
		v[2556] = v[9358 + i1687];
		v[2821] = -(v[2556] * v[6432]);
		v[2843] = v[2821] * v[6648];
		v[6803] = v[2843] * v[361];
		v[6801] = v[2843] * v[359];
		v[6672] = v[2843] * v[350];
		v[2642] = -0.5e0*v[2821];
		v[2557] = v[10734 + i1687];
		v[2567] = (i1687 == 18 ? 1 : 0);
		v[6786] = v[2567] * v[7];
		v[2568] = (i1687 == 17 ? 1 : 0);
		v[6784] = v[2568] * v[7];
		v[2569] = (i1687 == 16 ? 1 : 0);
		v[6785] = v[2569] * v[7];
		v[2570] = (i1687 == 12 ? 1 : 0);
		v[6791] = v[2570] * v[7];
		v[2571] = (i1687 == 11 ? 1 : 0);
		v[6789] = v[2571] * v[7];
		v[2572] = (i1687 == 10 ? 1 : 0);
		v[6790] = v[2572] * v[7];
		v[2573] = (i1687 == 6 ? 1 : 0);
		v[6796] = v[2573] * v[7];
		v[2574] = (i1687 == 5 ? 1 : 0);
		v[6794] = v[2574] * v[7];
		v[2575] = (i1687 == 4 ? 1 : 0);
		v[6795] = v[2575] * v[7];
		v[2576] = v[2529] + v[2573];
		v[6865] = 2e0*v[2576];
		v[2577] = v[2529] - v[2573];
		v[6866] = 2e0*v[2577];
		v[2578] = v[2531] + v[2575];
		v[6867] = 2e0*v[2578];
		v[2579] = v[2531] - v[2575];
		v[6868] = 2e0*v[2579];
		v[2580] = v[2533] - v[2574];
		v[6869] = 2e0*v[2580];
		v[2581] = v[2533] + v[2574];
		v[6870] = 2e0*v[2581];
		v[2583] = v[2536] * v[2582] - v[2573] * v[6649];
		v[2584] = -v[2536] + v[2574] * v[460];
		v[2586] = v[2536] * v[2585] + v[2574] * v[6649];
		v[2588] = v[2536] * v[2587] - v[2575] * v[6649];
		v[6652] = 2e0*(-(v[159] * v[2574]) + v[2590] * v[460]);
		v[2591] = -(v[159] * v[2575]) + v[2590] * v[457];
		v[2592] = -(v[159] * v[2573]) + v[2590] * v[462];
		v[2593] = -(v[2661] * v[461]) + v[2573] * v[6434];
		v[2594] = -(v[2661] * v[459]) + v[2575] * v[6434];
		v[2595] = v[2661] * v[458] - v[2574] * v[6434];
		v[2596] = -(v[2590] * v[515]) - v[2534] * v[6434];
		v[2741] = v[2596] * v[309];
		v[2597] = (-(v[2534] * v[477]) - v[2583] * v[515]) / 2e0;
		v[2598] = -(v[2590] * v[497]) - v[2584] * v[6434];
		v[2737] = v[2598] * v[307];
		v[2599] = (v[2584] * v[476] + v[2586] * v[497]) / 2e0;
		v[2600] = -(v[2590] * v[478]) - v[2537] * v[6434];
		v[2732] = v[2600] * v[301];
		v[2601] = (-(v[2537] * v[475]) - v[2588] * v[478]) / 2e0;
		v[2602] = v[2539] + v[2570];
		v[6851] = 2e0*v[2602];
		v[2603] = v[2539] - v[2570];
		v[6852] = 2e0*v[2603];
		v[2604] = v[2541] + v[2572];
		v[6853] = 2e0*v[2604];
		v[2605] = v[2541] - v[2572];
		v[6854] = 2e0*v[2605];
		v[2606] = v[2543] - v[2571];
		v[6855] = 2e0*v[2606];
		v[2607] = v[2543] + v[2571];
		v[6856] = 2e0*v[2607];
		v[2609] = v[2546] * v[2608] - v[2570] * v[6650];
		v[2610] = -v[2546] + v[2571] * v[466];
		v[2612] = v[2546] * v[2611] + v[2571] * v[6650];
		v[2614] = v[2546] * v[2613] - v[2572] * v[6650];
		v[6656] = 2e0*(-(v[178] * v[2571]) + v[2616] * v[466]);
		v[2617] = -(v[178] * v[2572]) + v[2616] * v[463];
		v[2618] = -(v[178] * v[2570]) + v[2616] * v[468];
		v[2619] = -(v[2699] * v[467]) + v[2570] * v[6441];
		v[2620] = -(v[2699] * v[465]) + v[2572] * v[6441];
		v[2621] = v[2699] * v[464] - v[2571] * v[6441];
		v[2622] = -(v[2616] * v[560]) - v[2544] * v[6441];
		v[2756] = v[2622] * v[335];
		v[2623] = (-(v[2544] * v[522]) - v[2609] * v[560]) / 2e0;
		v[2624] = -(v[2616] * v[542]) - v[2610] * v[6441];
		v[2752] = v[2624] * v[333];
		v[2625] = (v[2610] * v[521] + v[2612] * v[542]) / 2e0;
		v[2626] = -(v[2616] * v[523]) - v[2547] * v[6441];
		v[2747] = v[2626] * v[327];
		v[2627] = (-(v[2547] * v[520]) - v[2614] * v[523]) / 2e0;
		v[2628] = v[2549] + v[2567];
		v[6829] = 2e0*v[2628];
		v[2629] = v[2549] - v[2567];
		v[6830] = 2e0*v[2629];
		v[2630] = v[2551] + v[2569];
		v[6831] = 2e0*v[2630];
		v[2631] = v[2551] - v[2569];
		v[6832] = 2e0*v[2631];
		v[2632] = v[2553] - v[2568];
		v[6833] = 2e0*v[2632];
		v[2633] = v[2553] + v[2568];
		v[6834] = 2e0*v[2633];
		v[2635] = v[2556] * v[2634] - v[2567] * v[6651];
		v[2636] = -v[2556] + v[2568] * v[472];
		v[2638] = v[2556] * v[2637] + v[2568] * v[6651];
		v[2640] = v[2556] * v[2639] - v[2569] * v[6651];
		v[6666] = 2e0*(-(v[238] * v[2568]) + v[2642] * v[472]);
		v[2643] = -(v[238] * v[2569]) + v[2642] * v[469];
		v[2644] = -(v[238] * v[2567]) + v[2642] * v[474];
		v[2645] = -(v[2821] * v[473]) + v[2567] * v[6448];
		v[2646] = -(v[2821] * v[471]) + v[2569] * v[6448];
		v[2647] = v[2821] * v[470] - v[2568] * v[6448];
		v[2648] = -(v[2554] * v[6448]) - v[2642] * v[731];
		v[2869] = v[2648] * v[361];
		v[2649] = (-(v[2554] * v[693]) - v[2635] * v[731]) / 2e0;
		v[2650] = -(v[2636] * v[6448]) - v[2642] * v[713];
		v[2862] = v[2650] * v[359];
		v[2651] = (v[2636] * v[692] + v[2638] * v[713]) / 2e0;
		v[2652] = -(v[2557] * v[6448]) - v[2642] * v[694];
		v[2854] = v[2652] * v[353];
		v[2653] = (-(v[2557] * v[691]) - v[2640] * v[694]) / 2e0;
		v[2654] = (v[2537] * v[476] + v[2586] * v[478] + v[6652]) / 2e0;
		v[2655] = (v[2534] * v[476] + v[2586] * v[515] + v[6652]) / 2e0;
		v[2656] = v[2591] + v[2584] * v[6435] - v[2588] * v[6624];
		v[2657] = v[2591] + v[2534] * v[6435] - v[2588] * v[6623];
		v[2658] = v[159] * v[2578] + v[2661] * v[510];
		v[2659] = v[2661] - v[2578] * v[475] - v[2588] * v[510];
		v[2660] = v[159] * v[2580] + v[2661] * v[506];
		v[2743] = v[2660] * v[301];
		v[2662] = -v[2661] + v[2580] * v[476] + v[2586] * v[506];
		v[2663] = v[159] * v[2579] + v[2661] * v[502];
		v[2738] = v[2663] * v[309];
		v[2664] = -v[2661] - v[2579] * v[475] - v[2588] * v[502];
		v[2665] = v[2592] + v[2537] * v[6436] - v[2583] * v[6625];
		v[2666] = v[2592] + v[2584] * v[6436] - v[2583] * v[6624];
		v[2667] = v[159] * v[2576] + v[2661] * v[493];
		v[2668] = v[2661] - v[2576] * v[477] - v[2583] * v[493];
		v[2669] = v[159] * v[2581] + v[2661] * v[489];
		v[2733] = v[2669] * v[309];
		v[2670] = v[2661] + v[2581] * v[476] + v[2586] * v[489];
		v[2671] = -v[2593] - v[2581] * v[475] - v[2588] * v[489];
		v[2672] = -v[2593] + v[2579] * v[476] + v[2586] * v[502];
		v[2673] = -v[2593] + v[2578] * v[476] + v[2586] * v[510];
		v[2674] = -v[2593] - v[2580] * v[475] - v[2588] * v[506];
		v[2675] = v[2683] + v[2593] * v[6653];
		v[2676] = v[159] * v[2577] + v[2661] * v[483];
		v[2677] = -v[2661] - v[2577] * v[477] - v[2583] * v[483];
		v[2678] = -v[2594] - v[2581] * v[477] - v[2583] * v[489];
		v[2679] = -v[2594] + v[2577] * v[476] + v[2586] * v[483];
		v[2680] = -v[2594] + v[2576] * v[476] + v[2586] * v[493];
		v[2681] = -v[2594] - v[2580] * v[477] - v[2583] * v[506];
		v[2682] = v[2593] * v[485] + v[2594] * v[488];
		v[2684] = v[2683] + v[2594] * v[6654];
		v[2911] = v[1985] * v[2684];
		v[2685] = v[2595] - v[2577] * v[475] - v[2588] * v[483];
		v[2686] = v[2595] - v[2579] * v[477] - v[2583] * v[502];
		v[2687] = v[2595] - v[2576] * v[475] - v[2588] * v[493];
		v[2688] = v[2595] - v[2578] * v[477] - v[2583] * v[510];
		v[2689] = -(v[2594] * v[482]) - v[2595] * v[485];
		v[2690] = -(v[2593] * v[482]) - v[2595] * v[488];
		v[2691] = v[2683] + v[2595] * v[6655];
		v[2915] = v[1983] * v[2691];
		v[2692] = (v[2547] * v[521] + v[2612] * v[523] + v[6656]) / 2e0;
		v[2693] = (v[2544] * v[521] + v[2612] * v[560] + v[6656]) / 2e0;
		v[2694] = v[2617] + v[2610] * v[6442] - v[2614] * v[6606];
		v[2695] = v[2617] + v[2544] * v[6442] - v[2614] * v[6605];
		v[2696] = v[178] * v[2604] + v[2699] * v[555];
		v[2697] = v[2699] - v[2604] * v[520] - v[2614] * v[555];
		v[2698] = v[178] * v[2606] + v[2699] * v[551];
		v[2758] = v[2698] * v[327];
		v[2700] = -v[2699] + v[2606] * v[521] + v[2612] * v[551];
		v[2701] = v[178] * v[2605] + v[2699] * v[547];
		v[2753] = v[2701] * v[335];
		v[2702] = -v[2699] - v[2605] * v[520] - v[2614] * v[547];
		v[2703] = v[2618] + v[2547] * v[6443] - v[2609] * v[6607];
		v[2704] = v[2618] + v[2610] * v[6443] - v[2609] * v[6606];
		v[2705] = v[178] * v[2602] + v[2699] * v[538];
		v[2706] = v[2699] - v[2602] * v[522] - v[2609] * v[538];
		v[2707] = v[178] * v[2607] + v[2699] * v[534];
		v[2748] = v[2707] * v[335];
		v[2708] = v[2699] + v[2607] * v[521] + v[2612] * v[534];
		v[2709] = -v[2619] - v[2607] * v[520] - v[2614] * v[534];
		v[2710] = -v[2619] + v[2605] * v[521] + v[2612] * v[547];
		v[2711] = -v[2619] + v[2604] * v[521] + v[2612] * v[555];
		v[2712] = -v[2619] - v[2606] * v[520] - v[2614] * v[551];
		v[2713] = v[2721] + v[2619] * v[6657];
		v[2714] = v[178] * v[2603] + v[2699] * v[528];
		v[2715] = -v[2699] - v[2603] * v[522] - v[2609] * v[528];
		v[2716] = -v[2620] - v[2607] * v[522] - v[2609] * v[534];
		v[2717] = -v[2620] + v[2603] * v[521] + v[2612] * v[528];
		v[2718] = -v[2620] + v[2602] * v[521] + v[2612] * v[538];
		v[2719] = -v[2620] - v[2606] * v[522] - v[2609] * v[551];
		v[2720] = v[2619] * v[530] + v[2620] * v[533];
		v[2722] = v[2721] + v[2620] * v[6658];
		v[2939] = v[1979] * v[2722];
		v[2723] = v[2621] - v[2603] * v[520] - v[2614] * v[528];
		v[2724] = v[2621] - v[2605] * v[522] - v[2609] * v[547];
		v[2725] = v[2621] - v[2602] * v[520] - v[2614] * v[538];
		v[2726] = v[2621] - v[2604] * v[522] - v[2609] * v[555];
		v[2727] = -(v[2620] * v[527]) - v[2621] * v[530];
		v[2728] = -(v[2619] * v[527]) - v[2621] * v[533];
		v[2729] = v[2721] + v[2621] * v[6659];
		v[2943] = v[1977] * v[2729];
		v[2730] = v[2732] + v[2676] * v[307];
		v[2731] = v[2730] + v[2733] + v[6795];
		v[2734] = v[2732] + v[2733];
		v[2735] = v[2737] + v[2667] * v[301];
		v[2736] = v[2735] + v[2738] + v[6794];
		v[2739] = v[2737] + v[2738];
		v[2740] = v[2741] + v[2743];
		v[2742] = v[2741] + v[2658] * v[307];
		v[2744] = v[2742] + v[2743] + v[6796];
		v[2745] = v[2747] + v[2714] * v[333];
		v[2746] = v[2745] + v[2748] + v[6790];
		v[2749] = v[2747] + v[2748];
		v[2750] = v[2752] + v[2705] * v[327];
		v[2751] = v[2750] + v[2753] + v[6789];
		v[2754] = v[2752] + v[2753];
		v[2755] = v[2756] + v[2758];
		v[2757] = v[2756] + v[2696] * v[333];
		v[2759] = v[2757] + v[2758] + v[6791];
		v[2760] = i6665 * v[7];
		v[2761] = i6663 * v[7];
		v[2762] = i6661 * v[7];
		v[2763] = i6664 * v[7];
		v[2764] = i6662 * v[7];
		v[2765] = i6660 * v[7];
		v[6681] = v[223] * v[2762] + v[224] * v[2765];
		v[2766] = v[2660] * v[82] + v[2658] * v[85] + v[2596] * v[88];
		v[2767] = v[2660] * v[83] + v[2658] * v[86] + v[2596] * v[89];
		v[6843] = v[1721] * v[2767];
		v[2768] = v[2698] * v[91] + v[2696] * v[94] + v[2622] * v[97];
		v[2769] = v[2698] * v[92] + v[2696] * v[95] + v[2622] * v[98];
		v[6846] = v[1721] * v[2769];
		v[2770] = v[2667] * v[82] + v[2598] * v[85] + v[2663] * v[88];
		v[2771] = v[2667] * v[83] + v[2598] * v[86] + v[2663] * v[89];
		v[6844] = v[1748] * v[2771];
		v[2772] = v[2705] * v[91] + v[2624] * v[94] + v[2701] * v[97];
		v[2773] = v[2705] * v[92] + v[2624] * v[95] + v[2701] * v[98];
		v[6847] = v[1748] * v[2773];
		v[2774] = v[2600] * v[82] + v[2676] * v[85] + v[2669] * v[88];
		v[6849] = -(v[1721] * v[2766]) - v[1748] * v[2770] - v[1771] * v[2774];
		v[2775] = v[2600] * v[83] + v[2676] * v[86] + v[2669] * v[89];
		v[6845] = v[1771] * v[2775];
		v[2776] = v[2626] * v[91] + v[2714] * v[94] + v[2707] * v[97];
		v[6850] = -(v[1721] * v[2768]) - v[1748] * v[2772] - v[1771] * v[2776];
		v[2777] = v[2626] * v[92] + v[2714] * v[95] + v[2707] * v[98];
		v[6848] = v[1771] * v[2777];
		v[2778] = v[2678] * v[987] + v[2677] * v[988] + v[2665] * v[989];
		v[2779] = v[2670] * v[987] + v[2679] * v[988] + v[2654] * v[989];
		v[2780] = v[2671] * v[987] + v[2685] * v[988] + v[2601] * v[989];
		v[2781] = v[2716] * v[984] + v[2715] * v[985] + v[2703] * v[986];
		v[2782] = v[2708] * v[984] + v[2717] * v[985] + v[2692] * v[986];
		v[2783] = v[2709] * v[984] + v[2723] * v[985] + v[2627] * v[986];
		v[2784] = v[2686] * v[987] + v[2666] * v[988] + v[2668] * v[989];
		v[3326] = -(v[2784] * v[320]);
		v[2785] = v[2672] * v[987] + v[2599] * v[988] + v[2680] * v[989];
		v[3327] = -(v[2785] * v[310]);
		v[2786] = v[2664] * v[987] + v[2656] * v[988] + v[2687] * v[989];
		v[3328] = -(v[2786] * v[300]);
		v[2787] = v[2724] * v[984] + v[2704] * v[985] + v[2706] * v[986];
		v[3329] = -(v[2787] * v[346]);
		v[2788] = v[2710] * v[984] + v[2625] * v[985] + v[2718] * v[986];
		v[3330] = -(v[2788] * v[336]);
		v[2789] = v[2702] * v[984] + v[2694] * v[985] + v[2725] * v[986];
		v[3331] = -(v[2789] * v[326]);
		v[2790] = v[2597] * v[987] + v[2688] * v[988] + v[2681] * v[989];
		v[6711] = -(v[2790] * v[320]);
		v[2791] = v[2655] * v[987] + v[2673] * v[988] + v[2662] * v[989];
		v[6712] = -(v[2791] * v[310]);
		v[2792] = v[2009] * v[2597] + v[2008] * v[2655] + v[2007] * v[2657] + v[2013] * v[2664] + v[2020] * v[2670]
			+ v[2019] * v[2671] + v[2014] * v[2672] + v[2021] * v[2678] + v[2015] * v[2686];
		v[2793] = v[2014] * v[2599] + v[2013] * v[2656] + v[2007] * v[2659] + v[2015] * v[2666] + v[2008] * v[2673]
			+ v[2021] * v[2677] + v[2020] * v[2679] + v[2019] * v[2685] + v[2009] * v[2688];
		v[2794] = v[2657] * v[987] + v[2659] * v[988] + v[2674] * v[989];
		v[6713] = -(v[2794] * v[300]);
		v[2795] = v[2019] * v[2601] + v[2020] * v[2654] + v[2008] * v[2662] + v[2021] * v[2665] + v[2015] * v[2668]
			+ v[2007] * v[2674] + v[2014] * v[2680] + v[2009] * v[2681] + v[2013] * v[2687];
		v[2796] = v[2623] * v[984] + v[2726] * v[985] + v[2719] * v[986];
		v[6714] = -(v[2796] * v[346]);
		v[2797] = v[2693] * v[984] + v[2711] * v[985] + v[2700] * v[986];
		v[6715] = -(v[2797] * v[336]);
		v[2798] = v[2006] * v[2623] + v[2005] * v[2693] + v[2004] * v[2695] + v[2010] * v[2702] + v[2017] * v[2708]
			+ v[2016] * v[2709] + v[2011] * v[2710] + v[2018] * v[2716] + v[2012] * v[2724];
		v[2799] = v[2011] * v[2625] + v[2010] * v[2694] + v[2004] * v[2697] + v[2012] * v[2704] + v[2005] * v[2711]
			+ v[2018] * v[2715] + v[2017] * v[2717] + v[2016] * v[2723] + v[2006] * v[2726];
		v[2800] = v[2695] * v[984] + v[2697] * v[985] + v[2712] * v[986];
		v[6716] = -(v[2800] * v[326]);
		v[2801] = v[2016] * v[2627] + v[2017] * v[2692] + v[2005] * v[2700] + v[2018] * v[2703] + v[2012] * v[2706]
			+ v[2004] * v[2712] + v[2011] * v[2718] + v[2006] * v[2719] + v[2010] * v[2725];
		v[2802] = i6660 + v[227] * v[2768] + v[228] * v[2769];
		v[2803] = i6660;
		v[2804] = i6661 + v[227] * v[2766] + v[228] * v[2767];
		v[2805] = i6661;
		v[2806] = i6662 + v[227] * v[2772] + v[228] * v[2773];
		v[2807] = i6662;
		v[2808] = i6663 + v[227] * v[2770] + v[228] * v[2771];
		v[2809] = i6663;
		v[2810] = i6664 + v[227] * v[2776] + v[228] * v[2777];
		v[2811] = i6664;
		v[2812] = i6665 + v[227] * v[2774] + v[228] * v[2775];
		v[2813] = i6665;
		v[2814] = (v[6666] + v[2554] * v[692] + v[2638] * v[731]) / 2e0;
		v[2815] = (v[6666] + v[2557] * v[692] + v[2638] * v[694]) / 2e0;
		v[2816] = v[2643] + v[2554] * v[6449] - v[2640] * v[6575];
		v[2817] = v[2643] + v[2636] * v[6449] - v[2640] * v[6576];
		v[2818] = v[238] * v[2630] + v[2821] * v[726];
		v[2819] = v[2821] - v[2630] * v[691] - v[2640] * v[726];
		v[2820] = v[238] * v[2632] + v[2821] * v[722];
		v[2871] = v[2820] * v[353];
		v[2822] = -v[2821] + v[2632] * v[692] + v[2638] * v[722];
		v[2823] = v[238] * v[2631] + v[2821] * v[718];
		v[2863] = v[2823] * v[361];
		v[2824] = -v[2821] - v[2631] * v[691] - v[2640] * v[718];
		v[2825] = v[2644] + v[2636] * v[6450] - v[2635] * v[6576];
		v[2826] = v[2644] + v[2557] * v[6450] - v[2635] * v[6577];
		v[2827] = v[238] * v[2628] + v[2821] * v[709];
		v[2828] = v[2821] - v[2628] * v[693] - v[2635] * v[709];
		v[2829] = v[238] * v[2633] + v[2821] * v[705];
		v[2855] = v[2829] * v[361];
		v[2830] = v[2821] + v[2633] * v[692] + v[2638] * v[705];
		v[2831] = -v[2645] + v[2630] * v[692] + v[2638] * v[726];
		v[2832] = -v[2645] - v[2632] * v[691] - v[2640] * v[722];
		v[2833] = -v[2645] + v[2631] * v[692] + v[2638] * v[718];
		v[2834] = -v[2645] - v[2633] * v[691] - v[2640] * v[705];
		v[2835] = v[2843] + v[2645] * v[6667];
		v[2836] = v[238] * v[2629] + v[2821] * v[699];
		v[2837] = -v[2821] - v[2629] * v[693] - v[2635] * v[699];
		v[2838] = -v[2646] - v[2632] * v[693] - v[2635] * v[722];
		v[2839] = -v[2646] + v[2628] * v[692] + v[2638] * v[709];
		v[2840] = -v[2646] - v[2633] * v[693] - v[2635] * v[705];
		v[2841] = -v[2646] + v[2629] * v[692] + v[2638] * v[699];
		v[2842] = v[2645] * v[701] + v[2646] * v[704];
		v[2844] = v[2843] + v[2646] * v[6668];
		v[2967] = v[1973] * v[2844];
		v[2845] = v[2647] - v[2630] * v[693] - v[2635] * v[726];
		v[2846] = v[2647] - v[2631] * v[693] - v[2635] * v[718];
		v[2847] = v[2647] - v[2628] * v[691] - v[2640] * v[709];
		v[2848] = v[2647] - v[2629] * v[691] - v[2640] * v[699];
		v[2849] = -(v[2646] * v[698]) - v[2647] * v[701];
		v[2850] = -(v[2645] * v[698]) - v[2647] * v[704];
		v[2851] = v[2843] + v[2647] * v[6669];
		v[2971] = v[1971] * v[2851];
		v[2852] = v[2854] + v[2836] * v[359];
		v[2853] = v[2852] + v[2855] + v[6785];
		v[2856] = v[2854] + v[2855];
		v[2857] = v[128] * v[2652] + v[134] * v[2829] + v[131] * v[2836];
		v[2858] = v[129] * v[2652] + v[135] * v[2829] + v[132] * v[2836];
		v[2859] = v[130] * v[2652] + v[136] * v[2829] + v[133] * v[2836];
		v[2860] = v[2862] + v[2827] * v[353];
		v[2861] = v[2860] + v[2863] + v[6784];
		v[2864] = v[2862] + v[2863];
		v[2865] = v[131] * v[2650] + v[134] * v[2823] + v[128] * v[2827];
		v[2866] = v[132] * v[2650] + v[135] * v[2823] + v[129] * v[2827];
		v[2867] = v[133] * v[2650] + v[136] * v[2823] + v[130] * v[2827];
		v[2868] = v[2869] + v[2871];
		v[2870] = v[2869] + v[2818] * v[359];
		v[2872] = v[2870] + v[2871] + v[6786];
		v[2873] = v[134] * v[2648] + v[131] * v[2818] + v[128] * v[2820];
		v[2874] = v[135] * v[2648] + v[132] * v[2818] + v[129] * v[2820];
		v[2875] = v[136] * v[2648] + v[133] * v[2818] + v[130] * v[2820];
		v[2876] = i6673 * v[7];
		v[2877] = i6674 * v[7];
		v[2878] = i6675 * v[7];
		v[2879] = -(v[2649] * v[981]) - v[2845] * v[982] - v[2838] * v[983];
		v[6717] = v[2879] * v[372];
		v[2880] = -(v[2814] * v[981]) - v[2831] * v[982] - v[2822] * v[983];
		v[6718] = v[2880] * v[362];
		v[2881] = -(v[2816] * v[981]) - v[2819] * v[982] - v[2832] * v[983];
		v[6719] = v[2881] * v[352];
		v[6735] = -v[2878] + v[6681] - v[6711] - v[6712] - v[6713] - v[6714] - v[6715] - v[6716] - v[6717] - v[6718] - v[6719];
		v[2882] = -(v[2846] * v[981]) - v[2825] * v[982] - v[2828] * v[983];
		v[3332] = v[2882] * v[372];
		v[2883] = -(v[2833] * v[981]) - v[2651] * v[982] - v[2839] * v[983];
		v[3333] = v[2883] * v[362];
		v[2884] = -(v[2824] * v[981]) - v[2817] * v[982] - v[2847] * v[983];
		v[3334] = v[2884] * v[352];
		v[2885] = -(v[2840] * v[981]) - v[2837] * v[982] - v[2826] * v[983];
		v[2886] = -(v[2830] * v[981]) - v[2841] * v[982] - v[2815] * v[983];
		v[2887] = -(v[1993] * v[2653]) - v[1994] * v[2815] - v[2002] * v[2822] - v[1995] * v[2826] - v[1999] * v[2828]
			- v[2001] * v[2832] - v[2003] * v[2838] - v[1998] * v[2839] - v[1997] * v[2847];
		v[2888] = -(v[2003] * v[2649]) - v[2002] * v[2814] - v[2001] * v[2816] - v[1997] * v[2824] - v[1994] * v[2830]
			- v[1998] * v[2833] - v[1993] * v[2834] - v[1995] * v[2840] - v[1999] * v[2846];
		v[2889] = -(v[2834] * v[981]) - v[2848] * v[982] - v[2653] * v[983];
		v[6704] = -v[2876] + v[2780] * v[300] + v[2779] * v[310] + v[2778] * v[320] + v[2783] * v[326] + v[2782] * v[336]
			+ v[2781] * v[346] - v[2889] * v[352] - v[2886] * v[362] - v[2885] * v[372];
		v[6703] = v[223] * v[2760] + v[224] * v[2763] + v[6704];
		v[2890] = -(v[1998] * v[2651]) - v[1997] * v[2817] - v[2001] * v[2819] - v[1999] * v[2825] - v[2002] * v[2831]
			- v[1995] * v[2837] - v[1994] * v[2841] - v[2003] * v[2845] - v[1993] * v[2848];
		v[2891] = v[2586] + v[2682];
		v[2909] = v[1985] * v[2891];
		v[2892] = -v[2586] + v[2682];
		v[2894] = (v[173] * v[2891] + v[2658] * v[296] + v[2893] * v[6670]) / v[298];
		v[2895] = v[2583] + v[2689];
		v[2917] = v[1985] * v[2895];
		v[2896] = -v[2583] + v[2689];
		v[2913] = v[1983] * v[2896];
		v[2898] = (v[169] * v[2895] + v[2663] * v[297] + v[2897] * v[6670]) / v[298];
		v[2900] = (v[164] * v[2896] + v[2669] * v[6456] + v[2899] * v[6670]) / v[298];
		v[2901] = v[2911] + v[2913];
		v[6888] = v[2901] / v[298];
		v[2902] = v[2588] + v[2690];
		v[2907] = v[1983] * v[2902];
		v[2903] = -v[2588] + v[2690];
		v[2904] = v[2915] + v[2917];
		v[6882] = v[2904] / v[298];
		v[2905] = v[1981] * v[2675] + v[2907] + v[2909];
		v[6892] = v[2905] / v[298];
		v[2908] = v[2905] - v[2909];
		v[2910] = v[2905] - v[2907];
		v[6883] = v[2910] / v[298];
		v[2912] = v[1981] * v[2892] + v[2911];
		v[2914] = v[2912] + v[2913];
		v[6884] = v[2914] / v[298];
		v[2916] = v[1981] * v[2903] + v[2915];
		v[2918] = v[2916] + v[2917];
		v[6889] = v[2918] / v[298];
		v[2919] = v[2612] + v[2720];
		v[2937] = v[1979] * v[2919];
		v[2920] = -v[2612] + v[2720];
		v[2922] = (v[192] * v[2919] + v[2696] * v[322] + v[2921] * v[6671]) / v[324];
		v[2923] = v[2609] + v[2727];
		v[2945] = v[1979] * v[2923];
		v[2924] = -v[2609] + v[2727];
		v[2941] = v[1977] * v[2924];
		v[2926] = (v[188] * v[2923] + v[2701] * v[323] + v[2925] * v[6671]) / v[324];
		v[2928] = (v[183] * v[2924] + v[2707] * v[6464] + v[2927] * v[6671]) / v[324];
		v[2929] = v[2939] + v[2941];
		v[6902] = v[2929] / v[324];
		v[2930] = v[2614] + v[2728];
		v[2935] = v[1977] * v[2930];
		v[2931] = -v[2614] + v[2728];
		v[2932] = v[2943] + v[2945];
		v[6896] = v[2932] / v[324];
		v[2933] = v[1975] * v[2713] + v[2935] + v[2937];
		v[6906] = v[2933] / v[324];
		v[2936] = v[2933] - v[2937];
		v[2938] = v[2933] - v[2935];
		v[6897] = v[2938] / v[324];
		v[2940] = v[1975] * v[2920] + v[2939];
		v[2942] = v[2940] + v[2941];
		v[6898] = v[2942] / v[324];
		v[2944] = v[1975] * v[2931] + v[2943];
		v[2946] = v[2944] + v[2945];
		v[6903] = v[2946] / v[324];
		v[2947] = v[2638] + v[2842];
		v[2965] = v[1973] * v[2947];
		v[2948] = -v[2638] + v[2842];
		v[2950] = (v[252] * v[2947] + v[2818] * v[348] + v[2949] * v[6672]) / v[350];
		v[2951] = v[2635] + v[2849];
		v[2973] = v[1973] * v[2951];
		v[2952] = -v[2635] + v[2849];
		v[2969] = v[1971] * v[2952];
		v[2954] = (v[248] * v[2951] + v[2823] * v[349] + v[2953] * v[6672]) / v[350];
		v[2956] = (v[243] * v[2952] + v[2829] * v[6472] + v[2955] * v[6672]) / v[350];
		v[2957] = v[2967] + v[2969];
		v[6916] = v[2957] / v[350];
		v[2958] = v[2640] + v[2850];
		v[2963] = v[1971] * v[2958];
		v[2959] = -v[2640] + v[2850];
		v[2960] = v[2971] + v[2973];
		v[6910] = v[2960] / v[350];
		v[2961] = v[1969] * v[2835] + v[2963] + v[2965];
		v[6920] = v[2961] / v[350];
		v[2964] = v[2961] - v[2965];
		v[2966] = v[2961] - v[2963];
		v[6911] = v[2966] / v[350];
		v[2968] = v[1969] * v[2948] + v[2967];
		v[2970] = v[2968] + v[2969];
		v[6912] = v[2970] / v[350];
		v[2972] = v[1969] * v[2959] + v[2971];
		v[2974] = v[2972] + v[2973];
		v[6917] = v[2974] / v[350];
		v[2975] = v[224] * v[2810] + v[223] * v[2812] - v[1613] * v[2857];
		v[2976] = -(v[2132] * v[2857]) - v[1768] * v[2859] - v[1742] * v[2867] - v[1718] * v[2875];
		v[2975] = -(v[274] * v[2858]) + v[2975];
		v[2977] = -(v[2132] * v[2858]);
		v[2975] = -(v[1027] * v[2859]) + v[2975];
		v[2991] = -i6673 + v[2975];
		v[6676] = v[2991] * v[373];
		v[2978] = v[1768] * v[2857] - v[2132] * v[2859] + v[1742] * v[2865] + v[1718] * v[2873];
		v[2980] = v[224] * v[2806] + v[223] * v[2808] - v[1613] * v[2865];
		v[2981] = -(v[2130] * v[2865]) + v[2976];
		v[2980] = -(v[274] * v[2866]) + v[2980];
		v[2982] = -(v[2130] * v[2866]) + v[2977];
		v[2980] = -(v[1027] * v[2867]) + v[2980];
		v[2994] = -i6674 + v[2980];
		v[6677] = v[2994] * v[374];
		v[2983] = -(v[2130] * v[2867]) + v[2978];
		v[2985] = v[224] * v[2802] + v[223] * v[2804] - v[1613] * v[2873];
		v[2986] = -(v[2128] * v[2873]) + v[2981];
		v[2985] = -(v[274] * v[2874]) + v[2985];
		v[2987] = -(v[2128] * v[2874]) + v[2982];
		v[2985] = -(v[1027] * v[2875]) + v[2985];
		v[2996] = -i6675 + v[2985];
		v[6679] = v[2996] * v[375];
		v[2988] = -(v[2128] * v[2875]) + v[2983];
		v[2990] = v[6676] / v[401];
		v[2992] = v[6676] * v[6678];
		v[2993] = v[2991] * v[6556];
		v[3015] = v[2993];
		v[2975] = v[2991];
		v[2990] = v[2990] + v[6677] / v[401];
		v[2992] = v[2992] + v[6677] * v[6678];
		v[2995] = v[2994] * v[6556];
		v[3018] = v[2995];
		v[2980] = v[2994];
		v[2990] = v[2990] + v[6679] / v[401];
		v[3011] = v[2990];
		v[2992] = v[2992] + v[6678] * v[6679];
		v[3604] = v[2992];
		v[2997] = v[2996] * v[6556];
		v[3021] = v[2997];
		v[2985] = v[2996];
		v[2998] = 0e0;
		v[2999] = 0e0;
		v[3000] = 0e0;
		v[3001] = 0e0;
		v[3002] = 0e0;
		v[3003] = 0e0;
		v[3004] = 0e0;
		v[3005] = 0e0;
		v[3006] = 0e0;
		v[3007] = 0e0;
		b3008 = b6;
		if (b3008) {
			v[3009] = 0e0;
			v[3010] = v[2990] * v[387];
			v[2998] = v[2131] * v[2990] * v[386];
			v[3002] = v[2975] * v[388];
			v[2999] = v[2138] * v[2975];
			v[3002] = v[3002] + v[3010] * v[373];
			v[2993] = v[2993] + v[2138] * v[3010];
			v[3003] = v[2980] * v[388];
			v[2999] = v[2136] * v[2980] + v[2999];
			v[3003] = v[3003] + v[3010] * v[374];
			v[2995] = v[2995] + v[2136] * v[3010];
			v[3004] = v[2985] * v[388];
			v[2999] = v[2134] * v[2985] + v[2999];
			v[3004] = v[3004] + v[3010] * v[375];
			v[2997] = v[2997] + v[2134] * v[3010];
		}
		else {
			v[3012] = v[3011] * v[406];
			v[3009] = v[3011];
			v[3000] = v[2142] * v[3011] * v[405];
			v[2990] = 0e0;
			v[3001] = v[2139] * v[2975];
			v[3014] = v[3012] * v[373] + v[2975] * v[407];
			v[2993] = v[2139] * v[3012] + v[3015];
			v[3001] = v[2137] * v[2980] + v[3001];
			v[3017] = v[3012] * v[374] + v[2980] * v[407];
			v[2995] = v[2137] * v[3012] + v[3018];
			v[3001] = v[2135] * v[2985] + v[3001];
			v[3020] = v[3012] * v[375] + v[2985] * v[407];
			v[2997] = v[2135] * v[3012] + v[3021];
			b3022 = b418;
			if (b3022) {
				v[3002] = 0e0;
				v[3007] = -v[3014];
				v[3003] = 0e0;
				v[3006] = -v[3017];
				v[3004] = 0e0;
				v[3005] = -v[3020];
			}
			else {
				v[3007] = v[3014];
				v[3006] = v[3017];
				v[3005] = v[3020];
			};
		};
		v[4031] = v[2993];
		v[4030] = v[2995];
		v[4029] = v[2997];
		v[3004] = v[3004] + v[3005];
		v[6682] = v[1966] * v[3004];
		v[3029] = v[3004] * v[6680];
		v[3032] = v[3004] * v[5461] + v[391] * v[6681];
		v[3033] = 2e0*v[1815] * v[3004] + v[1966] * v[6681];
		v[3036] = v[2128] * v[2802] + v[2130] * v[2806] + v[2132] * v[2810] + v[1966] * (v[291] * v[3004] + v[2765] * v[391])
			+ v[1610] * (v[6846] + v[6847] + v[6848]) + v[1024] * v[6850];
		v[3037] = v[2128] * v[2804] + v[2130] * v[2808] + v[2132] * v[2812] + v[1966] * (v[288] * v[3004] + v[2762] * v[391])
			+ v[1610] * (v[6843] + v[6844] + v[6845]) + v[1024] * v[6849];
		v[3003] = v[3003] + v[3006];
		v[6684] = v[3003] * v[389];
		v[3038] = v[3003] * v[6683];
		v[3041] = v[3003] * v[3047];
		v[3042] = 2e0*v[1876] * v[3003];
		v[3002] = v[3002] + v[3007];
		v[6685] = v[3002] * v[390];
		v[3045] = v[224] * (v[6684] + v[6685]);
		v[3046] = v[223] * (v[6684] + v[6685]);
		v[3042] = v[3042] + v[3002] * v[3047];
		v[3050] = v[3002] * v[6686];
		v[3041] = 2e0*v[1943] * v[3002] + v[3041];
		v[3051] = 0e0;
		v[3052] = 0e0;
		v[3053] = 0e0;
		v[3054] = 0e0;
		v[3055] = 0e0;
		v[3056] = 0e0;
		v[3057] = 0e0;
		v[3058] = 0e0;
		v[3059] = 0e0;
		v[3060] = 0e0;
		v[3061] = 0e0;
		v[3062] = 0e0;
		v[3063] = 0e0;
		v[3064] = 0e0;
		v[3065] = 0e0;
		v[3066] = 0e0;
		v[3067] = 0e0;
		v[3068] = 0e0;
		v[3069] = 0e0;
		v[3070] = 0e0;
		v[3071] = 0e0;
		v[3072] = 0e0;
		v[3073] = 0e0;
		v[3074] = 0e0;
		v[3075] = 0e0;
		v[3076] = 0e0;
		v[3077] = 0e0;
		v[3078] = 0e0;
		v[3079] = 0e0;
		v[3080] = 0e0;
		v[3081] = 0e0;
		v[3082] = 0e0;
		b3083 = b941;
		if (b3083) {
			v[3084] = v[3004] * v[399];
			v[3085] = -(v[3004] * v[398]);
			v[3086] = v[3084] - v[3003] * v[400];
			v[3087] = v[3003] * v[398];
			v[3088] = v[3085] + v[3002] * v[400];
			v[3089] = v[3087] - v[3002] * v[399];
			v[6690] = v[2118] * v[3086] + v[2117] * v[3088] + v[2115] * v[3089];
			v[6688] = v[3086] * v[943] + v[3088] * v[944] + v[3089] * v[945];
			v[3090] = v[6688] / v[946];
			v[6689] = v[3090] * v[7306];
			v[3100] = v[3090] * v[6687];
			v[3054] = -(v[2121] * v[3091] * v[6688]);
			v[3092] = v[3086] * v[6500] + v[6689] * v[943];
			v[3107] = v[3092] * v[6748];
			v[3094] = v[3088] * v[6500] + v[6689] * v[944];
			v[3111] = 2e0*v[3094] * v[956];
			v[3095] = v[3089] * v[6500] + v[6689] * v[945];
			v[3108] = v[3095] * v[6745];
			v[3057] = v[3100] * v[6554] + v[6690] * v[948];
			v[3056] = v[3090] * v[3098] * v[952];
			v[3055] = v[3090] * v[6554] * v[953] + v[6690] * v[954];
			v[3081] = v[3100] * v[3101] * v[948];
			v[3082] = v[2120] * v[3090] * v[3096] * v[3101] * v[7307];
			v[3053] = v[3089] * v[6555] + v[2115] * v[6689];
			v[3052] = v[3088] * v[6555] + v[2117] * v[6689];
			v[3051] = v[3086] * v[6555] + v[2118] * v[6689];
			v[3103] = (v[3094] * v[955] + v[3092] * v[956]) / 2e0;
			v[3104] = v[3107] + v[3111];
			v[3105] = v[3104] + v[3108];
			v[3106] = (v[3095] * v[955] + v[3092] * v[957]) / 2e0;
			v[3109] = v[3107] + v[3108];
			v[3110] = (v[3095] * v[956] + v[3094] * v[957]) / 2e0;
			v[3112] = v[3108] + v[3111];
			v[3060] = (v[2101] * v[3092] + v[2088] * v[3094] + 4e0*v[3095] * v[3113]) / 2e0;
			v[3059] = (v[2108] * v[3092] + v[2088] * v[3095] + v[3094] * v[6691]) / 2e0;
			v[3058] = (v[2108] * v[3094] + v[2101] * v[3095] + v[3092] * v[6692]) / 2e0;
			v[3116] = v[3105] * v[6744];
			v[3080] = 8e0*v[2110] * v[3105] * v[7308];
			v[3079] = v[3116] * v[6693];
			v[3117] = v[3095] + v[3103];
			v[3118] = v[3095] - v[3103];
			v[3078] = v[2025] * v[3116];
			v[3119] = -v[3094] + v[3106];
			v[3120] = v[3094] + v[3106];
			v[3077] = v[2026] * v[3116];
			v[3065] = v[2095] * v[3116] + v[3117] * v[958];
			v[3076] = v[2027] * v[3116];
			v[3066] = (v[2091] * v[3116] - v[3109] * v[958]) / 2e0;
			v[3075] = v[3116] * v[6694];
			v[3121] = v[3092] + v[3110];
			v[3122] = -v[3092] + v[3110];
			v[3074] = v[2029] * v[3116];
			v[3068] = v[2082] * v[3116] + v[3119] * v[958];
			v[3073] = v[2030] * v[3116];
			v[3069] = v[2078] * v[3116] + v[3121] * v[958];
			v[3072] = v[2031] * v[3116];
			v[3070] = (v[2074] * v[3116] - v[3104] * v[958]) / 2e0;
			v[3071] = v[3116] * v[6695];
			v[3067] = v[2086] * v[3116] + v[3122] * v[958];
			v[3064] = v[2099] * v[3116] + v[3120] * v[958];
			v[3063] = v[2104] * v[3116] - v[3118] * v[958];
			v[3062] = (v[2109] * v[3116] - v[3112] * v[958]) / 2e0;
			v[3061] = v[2027] * v[3117] - v[2025] * v[3118] + v[2030] * v[3119] + v[2026] * v[3120] + v[2031] * v[3121]
				+ v[2029] * v[3122] - v[3112] * v[6693] - v[3109] * v[6694] - v[3104] * v[6695];
		}
		else {
		};
		v[3123] = 0e0;
		v[3124] = 0e0;
		v[3125] = 0e0;
		v[3126] = 0e0;
		v[3127] = 0e0;
		v[3128] = 0e0;
		b3129 = b4;
		if (b3129) {
			v[6698] = -(v[3002] * v[976]);
			v[6697] = -(v[3003] * v[977]);
			v[6696] = -(v[3004] * v[978]);
			v[3163] = v[1774] * v[3002];
			v[2805] = v[225] * v[2766] + v[2805];
			v[2805] = v[226] * v[2767] + v[2805];
			v[2803] = v[225] * v[2768] + v[2803];
			v[2803] = v[226] * v[2769] + v[2803];
			v[3132] = v[222] * v[2803] + v[221] * v[2805] - v[269] * v[2873];
			v[3133] = -(v[270] * v[2874]) + v[3132];
			v[3134] = -(v[272] * v[2875]) + v[3133];
			v[3136] = -i6675 + v[3] * v[3070] + v[3134];
			v[3070] = 0e0;
			v[3137] = v[2] * v[3069] + v[3136];
			v[3069] = 0e0;
			v[3138] = v[1] * v[3068] + v[3137];
			v[6700] = -(v[3138] * v[391]);
			v[3068] = 0e0;
			v[2809] = v[225] * v[2770] + v[2809];
			v[2809] = v[226] * v[2771] + v[2809];
			v[2807] = v[225] * v[2772] + v[2807];
			v[2807] = v[226] * v[2773] + v[2807];
			v[3141] = v[222] * v[2807] + v[221] * v[2809] - v[269] * v[2865];
			v[3142] = -(v[270] * v[2866]) + v[3141];
			v[3143] = -(v[272] * v[2867]) + v[3142];
			v[3145] = -i6674 + v[3] * v[3067] + v[3143];
			v[3067] = 0e0;
			v[3146] = v[2] * v[3066] + v[3145];
			v[3066] = 0e0;
			v[3147] = v[1] * v[3065] + v[3146];
			v[6699] = -(v[3147] * v[390]);
			v[3065] = 0e0;
			v[2813] = v[225] * v[2774] + v[2813];
			v[2813] = v[226] * v[2775] + v[2813];
			v[2811] = v[225] * v[2776] + v[2811];
			v[2811] = v[226] * v[2777] + v[2811];
			v[3150] = v[222] * v[2811] + v[221] * v[2813] - v[269] * v[2857];
			v[3151] = -(v[270] * v[2858]) + v[3150];
			v[3152] = -(v[272] * v[2859]) + v[3151];
			v[3154] = -i6673 + v[3] * v[3064] + v[3152];
			v[3064] = 0e0;
			v[3155] = v[2] * v[3063] + v[3154];
			v[3063] = 0e0;
			v[3156] = v[1] * v[3062] + v[3155];
			v[6701] = -(v[3156] * v[389]);
			v[3062] = 0e0;
			v[3157] = v[1772] * v[3004];
			v[3125] = v[3004] * v[3167];
			v[3041] = v[3041] + v[1774] * v[6696];
			v[3042] = v[3042] + v[1773] * v[6696];
			v[3159] = v[1773] * v[3003];
			v[3160] = v[3157] + v[3159];
			v[3124] = v[3003] * v[3169];
			v[3041] = v[3041] + v[1774] * v[6697];
			v[3033] = v[3033] + v[1772] * v[6697];
			v[3162] = v[3159] + v[3163];
			v[3164] = v[3157] + v[3163];
			v[3123] = v[3002] * v[3168];
			v[3042] = v[3042] + v[1773] * v[6698];
			v[3033] = v[3033] + v[1772] * v[6698];
			v[3123] = -(v[1774] * v[3050]) + v[3123];
			v[3128] = v[1772] * v[3138];
			v[3127] = v[1773] * v[3147];
			v[3126] = v[1774] * v[3156];
			v[3124] = -(v[1773] * v[3038]) + v[3124];
			v[3125] = -(v[1772] * v[3029]) + v[3125];
			v[3125] = v[3125] - v[3162] * v[391];
			v[3033] = v[3033] + v[3138] * v[3167] + v[1772] * (v[6699] + v[6701]) - v[3162] * v[978];
			v[3123] = v[3123] - v[3160] * v[389];
			v[3041] = v[3041] + v[3156] * v[3168] + v[1774] * (v[6699] + v[6700]) - v[3160] * v[976];
			v[3124] = v[3124] - v[3164] * v[390];
			v[3042] = v[3042] + v[3147] * v[3169] + v[1773] * (v[6700] + v[6701]) - v[3164] * v[977];
		}
		else {
		};
		v[6766] = v[1781] * v[3029];
		v[6710] = -(v[1780] * v[3038]);
		v[6709] = -(v[1887] * v[3002]);
		v[3243] = -(v[3003] * v[6702]);
		v[6708] = v[1781] * v[3004];
		v[6705] = v[3004] * v[3219];
		v[3221] = -(v[3004] * v[6702]);
		v[6707] = -(v[1887] * v[3004]) + v[3221];
		v[3170] = v[2878] + v[6711] + v[6712] + v[6713] + v[6714] + v[6715] + v[6716] + v[6717] + v[6718] + v[6719];
		v[6721] = -(v[3170] * v[391]);
		v[3171] = v[2877] + v[3326] + v[3327] + v[3328] + v[3329] + v[3330] + v[3331] + v[3332] + v[3333] + v[3334];
		v[6722] = v[3171] * v[390];
		v[3176] = v[11112 + i1687] + v[2894] * v[307] + v[2898] * v[309] + v[2683] * v[3172] + v[2684] * v[3173]
			+ v[2891] * v[3174] + v[2895] * v[3175] + v[2731] * v[4497] + v[2735] * v[6457] + v[2740] * v[6460];
		v[3181] = v[11148 + i1687] + v[2504] * v[2730] + v[2503] * v[2742] + v[2894] * v[301] + v[2900] * v[309]
			+ v[2683] * v[3177] + v[2691] * v[3178] + v[2896] * v[3179] + v[2902] * v[3180] + v[2736] * v[4500];
		v[3186] = v[11184 + i1687] + v[2506] * v[2734] + v[2505] * v[2739] + v[2898] * v[301] + v[2900] * v[307]
			+ v[2683] * v[3182] + v[2675] * v[3183] + v[2892] * v[3184] + v[2903] * v[3185] + v[2744] * v[4502];
		v[3191] = v[11220 + i1687] + v[2721] * v[3187] + v[2722] * v[3188] + v[2919] * v[3189] + v[2923] * v[3190]
			+ v[2922] * v[333] + v[2926] * v[335] + v[2746] * v[4504] + v[2750] * v[6465] + v[2755] * v[6468];
		v[3196] = v[11256 + i1687] + v[2512] * v[2745] + v[2511] * v[2757] + v[2721] * v[3192] + v[2729] * v[3193]
			+ v[2924] * v[3194] + v[2930] * v[3195] + v[2922] * v[327] + v[2928] * v[335] + v[2751] * v[4507];
		v[3201] = v[11292 + i1687] + v[2514] * v[2749] + v[2513] * v[2754] + v[2721] * v[3197] + v[2713] * v[3198]
			+ v[2920] * v[3199] + v[2931] * v[3200] + v[2926] * v[327] + v[2928] * v[333] + v[2759] * v[4509];
		v[3206] = v[11328 + i1687] + v[2843] * v[3202] + v[2844] * v[3203] + v[2947] * v[3204] + v[2951] * v[3205]
			+ v[2950] * v[359] + v[2954] * v[361] + v[2853] * v[4511] + v[2860] * v[6473] + v[2868] * v[6476];
		v[3211] = v[11364 + i1687] + v[2520] * v[2852] + v[2519] * v[2870] + v[2843] * v[3207] + v[2851] * v[3208]
			+ v[2952] * v[3209] + v[2958] * v[3210] + v[2950] * v[353] + v[2956] * v[361] + v[2861] * v[4514];
		v[3216] = v[11400 + i1687] + v[2522] * v[2856] + v[2521] * v[2864] + v[2843] * v[3212] + v[2835] * v[3213]
			+ v[2948] * v[3214] + v[2959] * v[3215] + v[2954] * v[353] + v[2956] * v[359] + v[2872] * v[4516];
		v[3218] = v[1016] * v[2761] + v[1017] * v[2764] + v[1984] * v[3176] + v[1982] * v[3181] + v[1980] * v[3186]
			+ v[1978] * v[3191] + v[1976] * v[3196] + v[1974] * v[3201] + v[1972] * v[3206] + v[1970] * v[3211] + v[1968] * v[3216]
			+ v[3004] * v[3217] + v[1014] * v[6703];
		v[3041] = v[3041] + v[1779] * v[6705];
		v[3218] = v[3218] + v[3003] * v[3240] + v[3032] * v[389];
		v[3041] = v[3041] + v[1779] * (v[3032] + v[3003] * v[3241]);
		v[3218] = v[3218] + v[3002] * v[3262];
		v[3263] = v[1779] * v[3002];
		v[3264] = v[300] * v[3263];
		v[3265] = v[310] * v[3263];
		v[3266] = v[320] * v[3263];
		v[3267] = v[326] * v[3263];
		v[3268] = v[3263] * v[336];
		v[3269] = v[3263] * v[346];
		v[3270] = v[3263] * v[352];
		v[3271] = v[3263] * v[362];
		v[3272] = v[3263] * v[372];
		v[3218] = v[290] * v[3045] + v[287] * v[3046] + v[3218] + v[3050] * v[3301];
		v[3311] = -(v[1779] * v[3050]);
		v[6723] = -(v[1779] * v[3170]) + v[1781] * v[6703];
		v[3335] = v[223] * v[2761] + v[224] * v[2764] - v[2877] - v[3326] - v[3327] - v[3328] - v[3329] - v[3330] - v[3331]
			- v[3332] - v[3333] - v[3334];
		v[6724] = v[3335] * v[390];
		v[3336] = v[3032] + v[389] * v[6704] + v[6721];
		v[3338] = v[1016] * v[2760] + v[1017] * v[2763] + v[1905] * v[3176] + v[1903] * v[3181] + v[1901] * v[3186]
			+ v[1899] * v[3191] + v[1897] * v[3196] + v[1895] * v[3201] + v[1893] * v[3206] + v[1891] * v[3211] + v[1889] * v[3216]
			+ v[1018] * v[3335] + v[3004] * v[3337];
		v[3042] = v[3042] + v[1780] * v[6705];
		v[3478] = v[372] * v[6708];
		v[3476] = v[362] * v[6708];
		v[3474] = v[352] * v[6708];
		v[3472] = v[346] * v[6708];
		v[3470] = v[336] * v[6708];
		v[3468] = v[326] * v[6708];
		v[3466] = v[320] * v[6708];
		v[3464] = v[310] * v[6708];
		v[3462] = v[300] * v[6708];
		v[3338] = v[3338] + v[3003] * v[3362] + v[3336] * v[390];
		v[3363] = v[1780] * v[3003];
		v[3364] = v[300] * v[3363];
		v[3365] = v[310] * v[3363];
		v[3366] = v[320] * v[3363];
		v[3367] = v[326] * v[3363];
		v[3368] = v[336] * v[3363];
		v[3369] = v[3363] * v[346];
		v[3370] = v[3363] * v[352];
		v[3371] = v[3363] * v[362];
		v[3372] = v[3363] * v[372];
		v[3373] = v[3263] + v[3363];
		v[3374] = v[3264] + v[3364];
		v[3375] = v[3265] + v[3365];
		v[3376] = v[3266] + v[3366];
		v[3377] = v[3267] + v[3367];
		v[3378] = v[3268] + v[3368];
		v[3379] = v[3269] + v[3369];
		v[3380] = v[3270] + v[3370];
		v[3381] = v[3271] + v[3371];
		v[3382] = v[3272] + v[3372];
		v[3338] = v[3338] + v[3038] * v[3401];
		v[3338] = v[3338] + v[3002] * v[3425];
		v[3042] = v[3042] + v[1780] * (v[3336] + v[3002] * v[3426]);
		v[3436] = v[3311] + v[6709];
		v[3338] = v[289] * v[3045] + v[286] * v[3046] + v[3338];
		v[3601] = v[3338];
		v[3447] = v[1832] * v[3176] + v[1830] * v[3181] + v[1828] * v[3186] + v[1826] * v[3191] + v[1824] * v[3196]
			+ v[1822] * v[3201] + v[1820] * v[3206] + v[1818] * v[3211] + v[1816] * v[3216] + v[3004] * v[3446] + v[1022] * v[6735];
		v[3448] = v[3363] + v[6708];
		v[3449] = v[3364] + v[3462];
		v[3450] = v[3365] + v[3464];
		v[3451] = v[3366] + v[3466];
		v[3452] = v[3367] + v[3468];
		v[3453] = v[3368] + v[3470];
		v[3454] = v[3369] + v[3472];
		v[3455] = v[3370] + v[3474];
		v[3456] = v[3371] + v[3476];
		v[3457] = v[3372] + v[3478];
		v[3041] = v[3041] + v[1887] * v[6704] + v[5483] * v[6708];
		v[3461] = v[3263] + v[6708];
		v[3463] = v[3264] + v[3462];
		v[3465] = v[3265] + v[3464];
		v[3467] = v[3266] + v[3466];
		v[3469] = v[3267] + v[3468];
		v[3471] = v[3268] + v[3470];
		v[3473] = v[3269] + v[3472];
		v[3475] = v[3270] + v[3474];
		v[3477] = v[3271] + v[3476];
		v[3479] = v[3272] + v[3478];
		v[3042] = v[3042] + v[5481] * v[6708];
		v[3447] = v[3447] + v[3029] * v[3506];
		v[3447] = v[3447] + v[3003] * v[3530];
		v[3531] = v[1781] * v[3003];
		v[6799] = v[3531] * v[391];
		v[6755] = -v[3243] - v[6710] + v[6799];
		v[3033] = v[3033] - v[1887] * v[3170] + v[3241] * v[3531];
		v[3447] = v[3447] + v[3002] * v[3552];
		v[3553] = v[1781] * v[3002];
		v[6752] = -(v[3553] * v[391]);
		v[3033] = v[3033] + v[3426] * v[3553];
		v[3447] = v[3447] + v[391] * (v[389] * v[6703] + v[6724]);
		v[3602] = v[3447];
		v[3042] = v[3042] - v[3171] * v[6702] + v[3335] * v[6720];
		v[3809] = v[3042];
		v[3218] = v[3218] - v[389] * (-v[6721] + v[6722]);
		v[3600] = v[3218];
		v[3041] = v[3041] - v[1779] * v[6722] + v[391] * v[6723];
		v[3811] = v[3041];
		v[3033] = v[3033] + v[389] * v[6723] + v[1781] * v[6724];
		v[3807] = v[3033];
		v[2990] = v[2990] + v[3009];
		v[3591] = v[2990];
		v[3574] = 0e0;
		v[3575] = 0e0;
		v[3576] = 0e0;
		v[3577] = 0e0;
		v[3578] = 0e0;
		v[3579] = 0e0;
		v[3580] = 0e0;
		v[3581] = 0e0;
		v[3582] = 0e0;
		b3583 = b6;
		if (b3583) {
			v[3586] = v[2990] * v[6725];
			v[3580] = -(v[3586] * v[3780] * v[7333]);
			v[2992] = v[2992] + (v[3585] * v[3586] * v[7334]) / v[1785];
			v[2990] = 0e0;
			v[3002] = 0e0;
			v[3003] = 0e0;
			v[3004] = 0e0;
			b3587 = b1040;
			if (b3587) {
				v[3218] = 0e0;
				v[3338] = 0e0;
				v[3447] = 0e0;
			}
			else {
			};
		}
		else {
			v[3588] = 0e0;
			b3589 = b1049;
			if (b3589) {
				b3590 = b1051;
				if (b3590) {
					v[2990] = v[3591];
					v[3588] = -v[3591];
					v[3002] = 0e0;
					v[3003] = 0e0;
					v[3004] = 0e0;
				}
				else {
				};
			}
			else {
			};
			b3595 = b1049;
			if (b3595) {
				v[3607] = v[1788] * v[3600] + v[1789] * v[3601] + v[1790] * v[3602];
				b3596 = b1051;
				if (b3596) {
					v[3795] = Power(v[1057], v[3585]);
					v[3598] = (v[3584] * v[3588]) / v[1799];
					v[3597] = v[3598] * v[3795];
					v[3581] = -((v[1804] * v[3597]) / v[1799]);
					v[3577] = v[1804] * v[3585] * v[3598] * Power(v[1057], v[3599]);
					v[3588] = 0e0;
					v[3218] = 0e0;
					v[3338] = 0e0;
					v[3578] = v[3607];
					v[3447] = 0e0;
				}
				else {
					v[3802] = Power(v[401], v[3606]);
					v[3605] = v[2990] * v[3801];
					v[3603] = v[3605] * v[3802];
					v[3597] = -(v[3603] / v[1805]);
					v[3582] = (v[1804] * v[3603]) / (v[1805] * v[1805]);
					v[2992] = v[3604] - (v[1804] * v[3605] * v[3606] * Power(v[401], -2e0 + v[1062])) / v[1805];
					v[2990] = 0e0;
					v[3218] = 0e0;
					v[3338] = 0e0;
					v[3579] = v[3607];
					v[3447] = 0e0;
				};
				v[3608] = v[3597] * v[6726];
				v[3574] = v[1788] * v[3608];
				v[3575] = v[1789] * v[3608];
				v[3576] = v[1790] * v[3608];
			}
			else {
			};
		};
		v[3800] = v[2992];
		v[3793] = v[3574];
		v[3792] = v[3575];
		v[3791] = v[3576];
		v[3609] = (v[223] * (-(v[1024] * (v[1691] * v[2766] + v[1690] * v[2770] + v[1689] * v[2774])) + v[1610] *
			(v[1691] * v[2767] + v[1690] * v[2771] + v[1689] * v[2775])) + v[224] * (-(v[1024] * (v[1691] * v[2768]
				+ v[1690] * v[2772] + v[1689] * v[2776])) + v[1610] * (v[1691] * v[2769] + v[1690] * v[2773] + v[1689] * v[2777]))
			)*v[33];
		v[3610] = ((v[1691] * (v[2802] - v[2804]) + v[1690] * (v[2806] - v[2808]) + v[1689] * (v[2810] - v[2812]))*v[33])
			/ 2e0;
		v[3611] = v[33] * (v[1689] * (v[1612] * v[2858] + v[2857] * v[449] + v[2859] * v[451]) + v[1690] * (v[1612] * v[2866]
			+ v[2865] * v[449] + v[2867] * v[451]) + v[1691] * (v[1612] * v[2874] + v[2873] * v[449] + v[2875] * v[451]));
		v[3612] = (v[1027] * (v[1689] * v[2857] + v[1690] * v[2865] + v[1691] * v[2873]) - v[1613] * (v[1689] * v[2859]
			+ v[1690] * v[2867] + v[1691] * v[2875]))*v[33];
		v[3613] = -(i6665*v[1614]) - i6663 * v[1618] - i6661 * v[1622] - i6664 * v[1638] - i6662 * v[1642] - i6660 * v[1646]
			- i6673 * v[1662] - i6674 * v[1666] - i6675 * v[1670] - v[1682] * v[2567] - v[1678] * v[2568] - v[1674] * v[2569]
			- v[1658] * v[2570] - v[1654] * v[2571] - v[1650] * v[2572] - v[1634] * v[2573] - v[1630] * v[2574] - v[1626] * v[2575];
		v[3629] = v[3613];
		v[3614] = -(i6665*v[1615]) - i6663 * v[1619] - i6661 * v[1623] - i6664 * v[1639] - i6662 * v[1643] - i6660 * v[1647]
			- i6673 * v[1663] - i6674 * v[1667] - i6675 * v[1671] - v[1683] * v[2567] - v[1679] * v[2568] - v[1675] * v[2569]
			- v[1659] * v[2570] - v[1655] * v[2571] - v[1651] * v[2572] - v[1635] * v[2573] - v[1631] * v[2574] - v[1627] * v[2575];
		v[3627] = v[3614];
		v[3615] = -(i6665*v[1616]) - i6663 * v[1620] - i6661 * v[1624] - i6664 * v[1640] - i6662 * v[1644] - i6660 * v[1648]
			- i6673 * v[1664] - i6674 * v[1668] - i6675 * v[1672] - v[1684] * v[2567] - v[1680] * v[2568] - v[1676] * v[2569]
			- v[1660] * v[2570] - v[1656] * v[2571] - v[1652] * v[2572] - v[1636] * v[2573] - v[1632] * v[2574] - v[1628] * v[2575];
		v[3625] = v[3615];
		v[3616] = 0e0;
		v[3617] = 0e0;
		v[3618] = 0e0;
		v[3619] = 0e0;
		v[3620] = 0e0;
		v[3621] = 0e0;
		b3622 = b1049;
		if (b3622) {
			b3623 = b5;
			if (b3623) {
				b3624 = b1080;
				if (b3624) {
					v[3621] = v[3615];
					v[3615] = 0e0;
					v[3620] = v[3614];
					v[3614] = 0e0;
					v[3619] = v[3613];
					v[3613] = 0e0;
				}
				else {
					v[3637] = v[3629] * v[6511];
					v[3635] = v[3627] * v[6511];
					v[3633] = v[3625] * v[6511];
					v[3615] = 0e0;
					v[3614] = 0e0;
					v[6728] = (v[35] * (v[1094] * v[3625] + v[1093] * v[3627] + v[1092] * v[3629])) / sqrt(v[3638]);
					v[3613] = 0e0;
					v[6727] = (v[1090] * (v[1077] * v[3633] + v[1076] * v[3635] + v[1075] * v[3637])) / sqrt(v[3634]);
					v[3621] = v[1091] * v[3633] + v[1077] * v[6727];
					v[3620] = v[1091] * v[3635] + v[1076] * v[6727];
					v[3619] = v[1091] * v[3637] + v[1075] * v[6727];
					v[3618] = v[1068] * v[6728];
					v[3617] = v[1067] * v[6728];
					v[3616] = v[1066] * v[6728];
				};
			}
			else {
				b3640 = b1114;
				if (b3640) {
					v[3621] = v[3625];
					v[3615] = 0e0;
					v[3620] = v[3627];
					v[3614] = 0e0;
					v[3619] = v[3629];
					v[3613] = 0e0;
				}
				else {
					v[3646] = v[1125] * v[35] * v[3629];
					v[3644] = v[1125] * v[35] * v[3627];
					v[3643] = v[1125] * v[35] * v[3625];
					v[3615] = 0e0;
					v[3614] = 0e0;
					v[6730] = (v[35] * (v[1124] * v[3625] + v[1123] * v[3627] + v[1122] * v[3629])) / sqrt(v[3638]);
					v[3613] = 0e0;
					v[6729] = (v[1120] * (v[1077] * v[3643] + v[1076] * v[3644] + v[1075] * v[3646])) / sqrt(v[3634]);
					v[3621] = v[1121] * v[3643] + v[1077] * v[6729];
					v[3620] = v[1121] * v[3644] + v[1076] * v[6729];
					v[3619] = v[1121] * v[3646] + v[1075] * v[6729];
					v[3618] = v[1068] * v[6730];
					v[3617] = v[1067] * v[6730];
					v[3616] = v[1066] * v[6730];
				};
			};
		}
		else {
		};
		v[3652] = -(v[1084] * v[2567]) - v[33] * (v[1691] * v[3216] + v[3621] * v[372]);
		v[3653] = -(v[1084] * v[2568]) - v[33] * (v[1691] * v[3211] + v[362] * v[3621]);
		v[3654] = -(v[1084] * v[2569]) - v[33] * (v[1691] * v[3206] + v[352] * v[3621]);
		v[3655] = -(i6675*v[1084]) - v[33] * (v[1691] * v[2878] + v[294] * v[3621]);
		v[3656] = -(i6674*v[1084]) - v[33] * (v[1691] * v[2877] + v[293] * v[3621]);
		v[3657] = -(i6673*v[1084]) - v[33] * (v[1691] * v[2876] + v[292] * v[3621]);
		v[3658] = -(v[1084] * v[2570]) - v[33] * (v[1691] * v[3201] + v[346] * v[3621]);
		v[3659] = -(v[1084] * v[2571]) - v[33] * (v[1691] * v[3196] + v[336] * v[3621]);
		v[3660] = -(v[1084] * v[2572]) - v[33] * (v[1691] * v[3191] + v[326] * v[3621]);
		v[3661] = -(i6660*v[1084]) - v[33] * (v[1691] * v[2765] + v[291] * v[3621]);
		v[3662] = -(i6662*v[1084]) - v[33] * (v[1691] * v[2764] + v[290] * v[3621]);
		v[3663] = -(i6664*v[1084]) - v[33] * (v[1691] * v[2763] + v[289] * v[3621]);
		v[3664] = -(v[1084] * v[2573]) - v[33] * (v[1691] * v[3186] + v[320] * v[3621]);
		v[3665] = -(v[1084] * v[2574]) - v[33] * (v[1691] * v[3181] + v[310] * v[3621]);
		v[3666] = -(v[1084] * v[2575]) - v[33] * (v[1691] * v[3176] + v[300] * v[3621]);
		v[3667] = -(i6661*v[1084]) - v[33] * (v[1691] * v[2762] + v[288] * v[3621]);
		v[3668] = -(i6663*v[1084]) - v[33] * (v[1691] * v[2761] + v[287] * v[3621]);
		v[3669] = -(i6665*v[1084]) - v[33] * (v[1691] * v[2760] + v[286] * v[3621]);
		v[3689] = -(v[1083] * v[2567]) - v[33] * (v[1690] * v[3216] + v[3620] * v[372]);
		v[3690] = -(v[1083] * v[2568]) - v[33] * (v[1690] * v[3211] + v[362] * v[3620]);
		v[3691] = -(v[1083] * v[2569]) - v[33] * (v[1690] * v[3206] + v[352] * v[3620]);
		v[3692] = -(i6675*v[1083]) - v[33] * (v[1690] * v[2878] + v[294] * v[3620]);
		v[3693] = -(i6674*v[1083]) - v[33] * (v[1690] * v[2877] + v[293] * v[3620]);
		v[3694] = -(i6673*v[1083]) - v[33] * (v[1690] * v[2876] + v[292] * v[3620]);
		v[3695] = -(v[1083] * v[2570]) - v[33] * (v[1690] * v[3201] + v[346] * v[3620]);
		v[3696] = -(v[1083] * v[2571]) - v[33] * (v[1690] * v[3196] + v[336] * v[3620]);
		v[3697] = -(v[1083] * v[2572]) - v[33] * (v[1690] * v[3191] + v[326] * v[3620]);
		v[3698] = -(i6660*v[1083]) - v[33] * (v[1690] * v[2765] + v[291] * v[3620]);
		v[3699] = -(i6662*v[1083]) - v[33] * (v[1690] * v[2764] + v[290] * v[3620]);
		v[3700] = -(i6664*v[1083]) - v[33] * (v[1690] * v[2763] + v[289] * v[3620]);
		v[3701] = -(v[1083] * v[2573]) - v[33] * (v[1690] * v[3186] + v[320] * v[3620]);
		v[3702] = -(v[1083] * v[2574]) - v[33] * (v[1690] * v[3181] + v[310] * v[3620]);
		v[3703] = -(v[1083] * v[2575]) - v[33] * (v[1690] * v[3176] + v[300] * v[3620]);
		v[3704] = -(i6661*v[1083]) - v[33] * (v[1690] * v[2762] + v[288] * v[3620]);
		v[3705] = -(i6663*v[1083]) - v[33] * (v[1690] * v[2761] + v[287] * v[3620]);
		v[3706] = -(i6665*v[1083]) - v[33] * (v[1690] * v[2760] + v[286] * v[3620]);
		v[3726] = -(v[1082] * v[2567]) - v[33] * (v[1689] * v[3216] + v[3619] * v[372]);
		v[3727] = -(v[1082] * v[2568]) - v[33] * (v[1689] * v[3211] + v[3619] * v[362]);
		v[3728] = -(v[1082] * v[2569]) - v[33] * (v[1689] * v[3206] + v[352] * v[3619]);
		v[3729] = -(i6675*v[1082]) - v[33] * (v[1689] * v[2878] + v[294] * v[3619]);
		v[3730] = -(i6674*v[1082]) - v[33] * (v[1689] * v[2877] + v[293] * v[3619]);
		v[3731] = -(i6673*v[1082]) - v[33] * (v[1689] * v[2876] + v[292] * v[3619]);
		v[3732] = -(v[1082] * v[2570]) - v[33] * (v[1689] * v[3201] + v[346] * v[3619]);
		v[3733] = -(v[1082] * v[2571]) - v[33] * (v[1689] * v[3196] + v[336] * v[3619]);
		v[3734] = -(v[1082] * v[2572]) - v[33] * (v[1689] * v[3191] + v[326] * v[3619]);
		v[3735] = -(i6660*v[1082]) - v[33] * (v[1689] * v[2765] + v[291] * v[3619]);
		v[3736] = -(i6662*v[1082]) - v[33] * (v[1689] * v[2764] + v[290] * v[3619]);
		v[3737] = -(i6664*v[1082]) - v[33] * (v[1689] * v[2763] + v[289] * v[3619]);
		v[3738] = -(v[1082] * v[2573]) - v[33] * (v[1689] * v[3186] + v[320] * v[3619]);
		v[3739] = -(v[1082] * v[2574]) - v[33] * (v[1689] * v[3181] + v[310] * v[3619]);
		v[3740] = -(v[1082] * v[2575]) - v[33] * (v[1689] * v[3176] + v[300] * v[3619]);
		v[3741] = -(i6661*v[1082]) - v[33] * (v[1689] * v[2762] + v[288] * v[3619]);
		v[3742] = -(i6663*v[1082]) - v[33] * (v[1689] * v[2761] + v[287] * v[3619]);
		v[3743] = -(i6665*v[1082]) - v[33] * (v[1689] * v[2760] + v[286] * v[3619]);
		v[3762] = v[32] * v[3621];
		v[3763] = v[32] * v[3620];
		v[3764] = v[32] * v[3619];
		v[3765] = v[3618];
		v[3766] = v[3618];
		v[3774] = v[3766];
		v[3767] = v[3617];
		v[3768] = v[3617];
		v[3775] = v[3768];
		v[3769] = v[3616];
		v[3770] = v[3616];
		v[3776] = v[3770];
		b3771 = b6;
		if (b3771) {
			b3772 = b1040;
			if (b3772) {
				v[3773] = v[1023] * v[3766];
				v[3576] = v[3576] + v[1037] * v[3766];
				v[3766] = 0e0;
				v[3773] = v[1021] * v[3768] + v[3773];
				v[3575] = v[3575] + v[1037] * v[3768];
				v[3768] = 0e0;
				v[3773] = v[1015] * v[3770] + v[3773];
				v[3574] = v[3574] + v[1037] * v[3770];
				v[3770] = 0e0;
			}
			else {
				v[3773] = 0e0;
				v[3765] = 0e0;
				v[3766] = 0e0;
				v[3767] = 0e0;
				v[3768] = 0e0;
				v[3769] = 0e0;
				v[3770] = 0e0;
			};
			v[3777] = v[3765] * v[6731];
			v[3033] = v[3033] + v[3765] * v[6507];
			v[3765] = 0e0;
			v[3778] = v[3777] + v[3767] * v[6732];
			v[3042] = v[3042] + v[3767] * v[6507];
			v[3767] = 0e0;
			v[3779] = v[3778] + v[3769] * v[6733];
			v[3041] = v[3041] + v[3769] * v[6507];
			v[3769] = 0e0;
			v[3580] = v[3580] + v[3773] * v[6734];
			v[2992] = v[2992] + (v[3580] * v[3780] * v[3794]) / v[1785] + v[29] * v[3779] * v[7341];
		}
		else {
			b3782 = b1049;
			if (b3782) {
				v[3783] = 0e0;
				v[3784] = 0e0;
				v[3785] = 0e0;
				b3786 = b1064;
				if (b3786) {
					v[3785] = v[3774];
					v[3766] = 0e0;
					v[3784] = v[3775];
					v[3768] = 0e0;
					v[3783] = v[3776];
					v[3770] = 0e0;
				}
				else {
					v[3765] = 0e0;
					v[3766] = 0e0;
					v[3767] = 0e0;
					v[3768] = 0e0;
					v[3769] = 0e0;
					v[3770] = 0e0;
				};
				v[3799] = v[1015] * v[3783];
				v[3798] = v[1021] * v[3784];
				v[3797] = v[1023] * v[3785];
				b3790 = b1051;
				if (b3790) {
					v[3578] = v[3578] + v[3797];
					v[3576] = v[1059] * v[3785] + v[3791];
					v[3578] = v[3578] + v[3798];
					v[3575] = v[1059] * v[3784] + v[3792];
					v[3578] = v[3578] + v[3799];
					v[3574] = v[1059] * v[3783] + v[3793];
					v[3581] = v[3581] + v[3578] * v[6734];
					v[3577] = v[3577] + (v[3581] * v[3794] * v[3795]) / v[1799];
				}
				else {
					v[3579] = v[3579] + v[3797];
					v[3576] = v[1063] * v[3785] + v[3791];
					v[3579] = v[3579] + v[3798];
					v[3575] = v[1063] * v[3784] + v[3792];
					v[3579] = v[3579] + v[3799];
					v[3574] = v[1063] * v[3783] + v[3793];
					v[3582] = v[3582] + v[3579] * v[6734];
					v[2992] = v[3800] - (v[3582] * v[36] * v[3801] * v[3802]) / (2e0*v[1805]);
				};
			}
			else {
			};
			v[3815] = v[2992];
			v[3814] = v[3769];
			v[3813] = v[3767];
			v[3812] = v[3765];
			b3804 = b1049;
			if (b3804) {
				b3805 = b1051;
				if (b3805) {
					v[3806] = -(v[3765] * v[6731]);
					v[3033] = v[3807] + v[3765] * v[6508];
					v[3765] = 0e0;
					v[3808] = v[3806] - v[3767] * v[6732];
					v[3042] = v[3809] + v[3767] * v[6508];
					v[3767] = 0e0;
					v[3810] = v[3808] - v[3769] * v[6733];
					v[3041] = v[3811] + v[3769] * v[6508];
					v[3769] = 0e0;
					v[3577] = v[3577] + v[29] * v[3810] * Power(v[1057], v[1045]);
					v[2992] = v[2992] - v[3577];
				}
				else {
					v[3033] = v[3807] + v[1054] * v[3812];
					v[3765] = 0e0;
					v[3042] = v[3809] + v[1054] * v[3813];
					v[3767] = 0e0;
					v[3041] = v[3811] + v[1054] * v[3814];
					v[3769] = 0e0;
					v[2992] = v[3815] - (v[3814] * v[389] + v[3813] * v[390] + v[3812] * v[391])*v[6509] * Power(v[401], v[1062]);
				};
			}
			else {
			};
		};
		v[4032] = v[2992];
		v[6813] = -(v[1014] * v[3574]);
		v[6814] = v[3311] + v[6813];
		v[6737] = -(v[3574] * v[389]);
		v[6815] = -(v[1018] * v[3575]);
		v[6816] = -v[6710] - v[6815];
		v[6736] = v[3575] * v[390];
		v[6817] = v[1022] * v[3576];
		v[6800] = -(v[3576] * v[391]);
		v[3816] = v[1781] * v[3216] + v[3576] * v[372];
		v[6763] = v[3816] * v[391];
		v[3817] = v[1781] * v[3211] + v[3576] * v[362];
		v[6762] = v[3817] * v[391];
		v[3818] = v[1781] * v[3206] + v[352] * v[3576];
		v[6761] = v[3818] * v[391];
		v[3819] = v[1781] * v[3201] + v[346] * v[3576];
		v[6760] = v[3819] * v[391];
		v[3820] = v[1781] * v[3196] + v[336] * v[3576];
		v[6759] = v[3820] * v[391];
		v[3821] = v[1781] * v[3191] + v[326] * v[3576];
		v[6758] = v[3821] * v[391];
		v[3822] = v[1781] * v[3186] + v[320] * v[3576];
		v[6757] = v[3822] * v[391];
		v[3823] = v[1781] * v[3181] + v[310] * v[3576];
		v[6756] = v[3823] * v[391];
		v[3824] = v[1781] * v[3176] + v[300] * v[3576];
		v[6754] = v[3824] * v[391];
		v[3041] = v[3041] + v[3576] * v[3825];
		v[3042] = v[3042] + v[3576] * v[3826];
		v[3827] = v[3576] * v[4566] + v[1781] * v[6735];
		v[3828] = v[3553] + v[3576] * v[389];
		v[6822] = -v[3311] + v[3828] * v[391] + v[389] * v[6708] - v[6813];
		v[3830] = v[3531] + v[3576] * v[390];
		v[6821] = v[3830] * v[391] + v[390] * v[6708] + v[6816];
		v[3033] = v[3033] + v[3576] * v[4572];
		v[3850] = v[1780] * v[3216] + v[3575] * v[372];
		v[6781] = v[3850] * v[390];
		v[3851] = v[1780] * v[3211] + v[3575] * v[362];
		v[6779] = v[3851] * v[390];
		v[3852] = v[1780] * v[3206] + v[352] * v[3575];
		v[6777] = v[3852] * v[390];
		v[3853] = v[1780] * v[3201] + v[346] * v[3575];
		v[6775] = v[3853] * v[390];
		v[3854] = v[1780] * v[3196] + v[336] * v[3575];
		v[6773] = v[3854] * v[390];
		v[3855] = v[1780] * v[3191] + v[326] * v[3575];
		v[6771] = v[3855] * v[390];
		v[3856] = v[1780] * v[3186] + v[320] * v[3575];
		v[6769] = v[3856] * v[390];
		v[3857] = v[1780] * v[3181] + v[310] * v[3575];
		v[6767] = v[3857] * v[390];
		v[3858] = v[1780] * v[3176] + v[300] * v[3575];
		v[6764] = v[3858] * v[390];
		v[3041] = v[3041] - v[292] * v[6736];
		v[3859] = v[1780] * v[3335] + v[3575] * v[4600];
		v[3042] = v[3042] + v[3575] * v[4602];
		v[3033] = v[3033] - v[294] * v[6736];
		v[3880] = v[1779] * v[3216] + v[3574] * v[372];
		v[6782] = v[3880] * v[389];
		v[3881] = v[1779] * v[3211] + v[3574] * v[362];
		v[6780] = v[3881] * v[389];
		v[3882] = v[1779] * v[3206] + v[352] * v[3574];
		v[6778] = v[3882] * v[389];
		v[3883] = v[1779] * v[3201] + v[346] * v[3574];
		v[6776] = v[3883] * v[389];
		v[3884] = v[1779] * v[3196] + v[336] * v[3574];
		v[6774] = v[3884] * v[389];
		v[3885] = v[1779] * v[3191] + v[326] * v[3574];
		v[6772] = v[3885] * v[389];
		v[3886] = v[1779] * v[3186] + v[320] * v[3574];
		v[6770] = v[3886] * v[389];
		v[3887] = v[1779] * v[3181] + v[310] * v[3574];
		v[6768] = v[3887] * v[389];
		v[3888] = v[1779] * v[3176] + v[300] * v[3574];
		v[6765] = v[3888] * v[389];
		v[3889] = v[3574] * v[4632] + v[1779] * v[6703];
		v[3890] = v[1780] * v[2760] + v[1779] * v[2761] + v[287] * v[3574] + v[286] * v[3575];
		v[3891] = v[1780] * v[2763] + v[1779] * v[2764] + v[290] * v[3574] + v[289] * v[3575];
		v[6753] = v[223] * v[3890] + v[224] * v[3891];
		v[3041] = v[3041] + v[3574] * v[4636];
		v[3042] = v[3042] + v[293] * v[6737];
		v[3033] = v[3033] + v[294] * v[6737];
		v[3912] = 0e0;
		v[3913] = 0e0;
		v[3914] = 0e0;
		v[3915] = 0e0;
		v[3916] = 0e0;
		v[3917] = 0e0;
		v[3918] = 0e0;
		v[3919] = 0e0;
		v[3920] = 0e0;
		b3921 = b4;
		if (b3921) {
			v[6740] = v[3764] * v[389];
			v[6739] = v[3763] * v[390];
			v[6743] = v[6739] + v[6740];
			v[6738] = v[3762] * v[391];
			v[6742] = v[6738] + v[6740];
			v[6741] = v[6738] + v[6739];
			v[3128] = v[3128] + v[3762] * v[978];
			v[3127] = v[3127] + v[3763] * v[977];
			v[3123] = v[3123] + v[2034] * v[3764] - v[389] * v[6741];
			v[3124] = v[3124] + v[2036] * v[3763] - v[390] * v[6742];
			v[3125] = v[3125] + v[2038] * v[3762] - v[391] * v[6743];
			v[3126] = v[3126] + v[3764] * v[976];
			v[3041] = v[3041] - v[3764] * v[6550] - v[6741] * v[976];
			v[3042] = v[3042] + v[3763] * v[6551] - v[6742] * v[977];
			v[3033] = v[3033] + v[3762] * v[6552] - v[6743] * v[978];
			v[3912] = v[1] * v[3123];
			v[3913] = v[2] * v[3123];
			v[3914] = v[3] * v[3123];
			v[3922] = -v[3123];
			v[3923] = -(v[272] * v[3123]);
			v[3924] = -(v[270] * v[3123]);
			v[3925] = -(v[269] * v[3123]);
			v[3926] = v[222] * v[3123];
			v[3927] = v[221] * v[3123];
			v[3928] = v[226] * v[3926];
			v[3929] = v[225] * v[3926];
			v[3930] = v[226] * v[3927];
			v[3931] = v[225] * v[3927];
			v[3915] = v[1] * v[3124];
			v[3916] = v[2] * v[3124];
			v[3917] = v[3] * v[3124];
			v[3932] = -v[3124];
			v[3933] = -(v[272] * v[3124]);
			v[3934] = -(v[270] * v[3124]);
			v[3935] = -(v[269] * v[3124]);
			v[3936] = v[222] * v[3124];
			v[3937] = v[221] * v[3124];
			v[3938] = v[226] * v[3936];
			v[3939] = v[225] * v[3936];
			v[3940] = v[226] * v[3937];
			v[3941] = v[225] * v[3937];
			v[3918] = v[1] * v[3125];
			v[3919] = v[2] * v[3125];
			v[3920] = v[3] * v[3125];
			v[3942] = -v[3125];
			v[3943] = -(v[272] * v[3125]);
			v[3944] = -(v[270] * v[3125]);
			v[3945] = -(v[269] * v[3125]);
			v[3946] = v[222] * v[3125];
			v[3947] = v[221] * v[3125];
			v[3948] = v[226] * v[3946];
			v[3949] = v[225] * v[3946];
			v[3950] = v[226] * v[3947];
			v[3951] = v[225] * v[3947];
			v[3889] = -v[3126] + v[3889];
			v[3859] = -v[3127] + v[3859];
			v[3827] = -v[3128] + v[3827];
		}
		else {
			v[3931] = 0e0;
			v[3930] = 0e0;
			v[3941] = 0e0;
			v[3940] = 0e0;
			v[3951] = 0e0;
			v[3950] = 0e0;
			v[3929] = 0e0;
			v[3928] = 0e0;
			v[3939] = 0e0;
			v[3938] = 0e0;
			v[3949] = 0e0;
			v[3948] = 0e0;
			v[3927] = 0e0;
			v[3937] = 0e0;
			v[3947] = 0e0;
			v[3926] = 0e0;
			v[3936] = 0e0;
			v[3946] = 0e0;
			v[3925] = 0e0;
			v[3924] = 0e0;
			v[3923] = 0e0;
			v[3935] = 0e0;
			v[3934] = 0e0;
			v[3933] = 0e0;
			v[3945] = 0e0;
			v[3944] = 0e0;
			v[3943] = 0e0;
			v[3922] = 0e0;
			v[3932] = 0e0;
			v[3942] = 0e0;
		};
		b3952 = b941;
		if (b3952) {
			v[3061] = v[3061] + (v[2074] * v[3920]) / 2e0;
			v[3071] = v[3071] + v[3920] * v[6502];
			v[3061] = v[3061] + v[2078] * v[3919];
			v[3072] = v[3072] + v[3919] * v[958];
			v[3061] = v[3061] + v[2082] * v[3918];
			v[3073] = v[3073] + v[3918] * v[958];
			v[3061] = v[3061] + v[2086] * v[3917];
			v[3074] = v[3074] + v[3917] * v[958];
			v[3061] = v[3061] + (v[2091] * v[3916]) / 2e0;
			v[3075] = v[3075] + v[3916] * v[6502];
			v[3061] = v[3061] + v[2095] * v[3915];
			v[3076] = v[3076] + v[3915] * v[958];
			v[3061] = v[3061] + v[2099] * v[3914];
			v[3077] = v[3077] + v[3914] * v[958];
			v[3061] = v[3061] + v[2104] * v[3913];
			v[3078] = v[3078] + v[3913] * v[958];
			v[3061] = v[3061] + (v[2109] * v[3912]) / 2e0;
			v[3079] = v[3079] + v[3912] * v[6502];
			v[3080] = v[3080] + v[3061] * v[6744];
			v[6749] = -v[3075] + v[3080];
			v[3059] = v[3059] - v[3073];
			v[3962] = v[3073] + v[3077];
			v[3059] = v[3059] + v[3077];
			v[3058] = v[3058] + v[3072] + (v[3962] * v[957]) / 2e0;
			v[3964] = v[3072] + v[3074];
			v[3058] = v[3058] - v[3074];
			v[3060] = v[3060] + v[3076] + v[3964] * v[6501] + v[3962] * v[6553] + v[6745] * (-v[3079] + v[6749]);
			v[3060] = v[3060] - v[3078];
			v[6746] = v[3060] * v[945];
			v[3966] = v[3076] + v[3078];
			v[3057] = v[3057] + v[6746] * v[948];
			v[3055] = v[3055] + v[6746] * v[954];
			v[3053] = v[3053] + v[3060] * v[6500];
			v[3059] = v[3059] + (v[3966] * v[955] - 4e0*(v[3071] + v[3079] - v[3080])*v[956] + v[3964] * v[957]) / 2e0;
			v[6747] = v[3059] * v[944];
			v[3057] = v[3057] + v[6747] * v[948];
			v[3055] = v[3055] + v[6747] * v[954];
			v[3052] = v[3052] + v[3059] * v[6500];
			v[3058] = v[3058] + v[3966] * v[6501] + v[6748] * (-v[3071] + v[6749]);
			v[6750] = v[3058] * v[943];
			v[3057] = v[3057] + v[6750] * v[948];
			v[3055] = v[3055] + v[6750] * v[954];
			v[3051] = v[3051] + v[3058] * v[6500];
			v[3056] = v[3056] + v[3057] * v[953];
			v[3054] = v[3054] + v[3056];
			v[3081] = v[3081] + 2e0*v[3055] * v[3096];
			v[3082] = v[3082] + (v[3081] * v[3097]) / 2e0;
			v[3054] = v[3054] + v[3082];
			v[6751] = v[3054] / v[946];
			v[3053] = v[3053] + v[6751] * v[945];
			v[3052] = v[3052] + v[6751] * v[944];
			v[3051] = v[3051] + v[6751] * v[943];
			v[3041] = v[3041] - v[3053] * v[399];
			v[3042] = v[3042] + v[3053] * v[398];
			v[3041] = v[3041] + v[3052] * v[400];
			v[3033] = v[3033] - v[3052] * v[398];
			v[3042] = v[3042] - v[3051] * v[400];
			v[3033] = v[3033] + v[3051] * v[399];
		}
		else {
		};
		v[3972] = -(v[3743] * v[905]) - v[3742] * v[906] - v[3741] * v[907] - v[3740] * v[908] - v[3739] * v[909]
			- v[3738] * v[910] - v[3737] * v[911] - v[3736] * v[912] - v[3735] * v[913] - v[3734] * v[914] - v[3733] * v[915]
			- v[3732] * v[916] - v[3731] * v[917] - v[3730] * v[918] - v[3729] * v[919] - v[3728] * v[920] - v[3727] * v[921]
			- v[3726] * v[922];
		v[3973] = -(v[3743] * v[923]) - v[3742] * v[924] - v[3741] * v[925] - v[3740] * v[926] - v[3739] * v[927]
			- v[3738] * v[928] - v[3737] * v[929] - v[3736] * v[930] - v[3735] * v[931] - v[3734] * v[932] - v[3733] * v[933]
			- v[3732] * v[934] - v[3731] * v[935] - v[3730] * v[936] - v[3729] * v[937] - v[3728] * v[938] - v[3727] * v[939]
			- v[3726] * v[940];
		v[4242] = (-(v[3743] * v[869]) - v[3742] * v[870] - v[3741] * v[871] - v[3740] * v[872] - v[3739] * v[873]
			- v[3738] * v[874] - v[3737] * v[875] - v[3736] * v[876] - v[3735] * v[877] - v[3734] * v[878] - v[3733] * v[879]
			- v[3732] * v[880] - v[3731] * v[881] - v[3730] * v[882] - v[3729] * v[883] - v[3728] * v[884] - v[3727] * v[885]
			- v[3726] * v[886]) / 2e0;
		v[3975] = v[3743] * v[887] + v[3742] * v[888] + v[3741] * v[889] + v[3740] * v[890] + v[3739] * v[891] + v[3738] * v[892]
			+ v[3737] * v[893] + v[3736] * v[894] + v[3735] * v[895] + v[3734] * v[896] + v[3733] * v[897] + v[3732] * v[898]
			+ v[3731] * v[899] + v[3730] * v[900] + v[3729] * v[901] + v[3728] * v[902] + v[3727] * v[903] + v[3726] * v[904];
		v[6807] = v[223] * v[3975];
		v[6806] = v[224] * v[3975];
		v[3976] = -(v[3706] * v[905]) - v[3705] * v[906] - v[3704] * v[907] - v[3703] * v[908] - v[3702] * v[909]
			- v[3701] * v[910] - v[3700] * v[911] - v[3699] * v[912] - v[3698] * v[913] - v[3697] * v[914] - v[3696] * v[915]
			- v[3695] * v[916] - v[3694] * v[917] - v[3693] * v[918] - v[3692] * v[919] - v[3691] * v[920] - v[3690] * v[921]
			- v[3689] * v[922];
		v[3977] = -(v[3706] * v[923]) - v[3705] * v[924] - v[3704] * v[925] - v[3703] * v[926] - v[3702] * v[927]
			- v[3701] * v[928] - v[3700] * v[929] - v[3699] * v[930] - v[3698] * v[931] - v[3697] * v[932] - v[3696] * v[933]
			- v[3695] * v[934] - v[3694] * v[935] - v[3693] * v[936] - v[3692] * v[937] - v[3691] * v[938] - v[3690] * v[939]
			- v[3689] * v[940];
		v[4245] = (-(v[3706] * v[869]) - v[3705] * v[870] - v[3704] * v[871] - v[3703] * v[872] - v[3702] * v[873]
			- v[3701] * v[874] - v[3700] * v[875] - v[3699] * v[876] - v[3698] * v[877] - v[3697] * v[878] - v[3696] * v[879]
			- v[3695] * v[880] - v[3694] * v[881] - v[3693] * v[882] - v[3692] * v[883] - v[3691] * v[884] - v[3690] * v[885]
			- v[3689] * v[886]) / 2e0;
		v[3979] = v[3706] * v[887] + v[3705] * v[888] + v[3704] * v[889] + v[3703] * v[890] + v[3702] * v[891] + v[3701] * v[892]
			+ v[3700] * v[893] + v[3699] * v[894] + v[3698] * v[895] + v[3697] * v[896] + v[3696] * v[897] + v[3695] * v[898]
			+ v[3694] * v[899] + v[3693] * v[900] + v[3692] * v[901] + v[3691] * v[902] + v[3690] * v[903] + v[3689] * v[904];
		v[6809] = v[223] * v[3979];
		v[6808] = v[224] * v[3979];
		v[3980] = -(v[3669] * v[905]) - v[3668] * v[906] - v[3667] * v[907] - v[3666] * v[908] - v[3665] * v[909]
			- v[3664] * v[910] - v[3663] * v[911] - v[3662] * v[912] - v[3661] * v[913] - v[3660] * v[914] - v[3659] * v[915]
			- v[3658] * v[916] - v[3657] * v[917] - v[3656] * v[918] - v[3655] * v[919] - v[3654] * v[920] - v[3653] * v[921]
			- v[3652] * v[922];
		v[3981] = -(v[3669] * v[923]) - v[3668] * v[924] - v[3667] * v[925] - v[3666] * v[926] - v[3665] * v[927]
			- v[3664] * v[928] - v[3663] * v[929] - v[3662] * v[930] - v[3661] * v[931] - v[3660] * v[932] - v[3659] * v[933]
			- v[3658] * v[934] - v[3657] * v[935] - v[3656] * v[936] - v[3655] * v[937] - v[3654] * v[938] - v[3653] * v[939]
			- v[3652] * v[940];
		v[4248] = (-(v[3669] * v[869]) - v[3668] * v[870] - v[3667] * v[871] - v[3666] * v[872] - v[3665] * v[873]
			- v[3664] * v[874] - v[3663] * v[875] - v[3662] * v[876] - v[3661] * v[877] - v[3660] * v[878] - v[3659] * v[879]
			- v[3658] * v[880] - v[3657] * v[881] - v[3656] * v[882] - v[3655] * v[883] - v[3654] * v[884] - v[3653] * v[885]
			- v[3652] * v[886]) / 2e0;
		v[3983] = v[3669] * v[887] + v[3668] * v[888] + v[3667] * v[889] + v[3666] * v[890] + v[3665] * v[891] + v[3664] * v[892]
			+ v[3663] * v[893] + v[3662] * v[894] + v[3661] * v[895] + v[3660] * v[896] + v[3659] * v[897] + v[3658] * v[898]
			+ v[3657] * v[899] + v[3656] * v[900] + v[3655] * v[901] + v[3654] * v[902] + v[3653] * v[903] + v[3652] * v[904];
		v[6811] = v[223] * v[3983];
		v[6810] = v[224] * v[3983];
		v[3984] = v[3263] * v[3283] + v[1968] * v[3574] + v[1889] * v[3575] + v[1816] * v[3576] + v[33] * (-(v[1682] * v[3619])
			- v[1683] * v[3620] - v[1684] * v[3621]) + v[3363] * v[6482] - v[2885] * v[6547] - v[2882] * v[6548] + v[2879] * v[6706]
			+ v[3482] * v[6708] + (v[3311] + v[6709])*v[783] + (v[3243] + v[6710])*v[786] + v[391] * (-(v[3553] * v[783])
				- v[3531] * v[786]) + v[6707] * v[789] + v[1781] * (-(v[1022] * v[2879]) - v[3029] * v[789]) - v[3610] * v[886]
			- v[3609] * v[904] + v[3612] * v[922] + v[3611] * v[940];
		v[4189] = v[3984] * v[6453];
		v[3985] = v[3263] * v[3285] + v[1970] * v[3574] + v[1891] * v[3575] + v[1818] * v[3576] + v[33] * (-(v[1678] * v[3619])
			- v[1679] * v[3620] - v[1680] * v[3621]) + v[3363] * v[6484] - v[2886] * v[6547] - v[2883] * v[6548] + v[2880] * v[6706]
			+ v[3484] * v[6708] + (v[3311] + v[6709])*v[782] + (v[3243] + v[6710])*v[785] + v[391] * (-(v[3553] * v[782])
				- v[3531] * v[785]) + v[6707] * v[788] + v[1781] * (-(v[1022] * v[2880]) - v[3029] * v[788]) - v[3610] * v[885]
			- v[3609] * v[903] + v[3612] * v[921] + v[3611] * v[939];
		v[6787] = v[350] * v[3985];
		v[4193] = v[3985] * v[6452];
		v[3986] = v[3263] * v[3287] + v[1972] * v[3574] + v[1893] * v[3575] + v[1820] * v[3576] + v[33] * (-(v[1674] * v[3619])
			- v[1675] * v[3620] - v[1676] * v[3621]) + v[3363] * v[6486] - v[2889] * v[6547] - v[2884] * v[6548] + v[2881] * v[6706]
			+ v[3486] * v[6708] + (v[3311] + v[6709])*v[781] + (v[3243] + v[6710])*v[784] + v[391] * (-(v[3553] * v[781])
				- v[3531] * v[784]) + v[6707] * v[787] + v[1781] * (-(v[1022] * v[2881]) - v[3029] * v[787]) - v[3610] * v[884]
			- v[3609] * v[902] + v[3612] * v[920] + v[3611] * v[938];
		v[6788] = v[350] * v[3986];
		v[4196] = v[3986] * v[6451];
		v[3987] = v[3263] * v[3289] + v[1974] * v[3574] + v[1895] * v[3575] + v[1822] * v[3576] + v[33] * (-(v[1658] * v[3619])
			- v[1659] * v[3620] - v[1660] * v[3621]) + v[3363] * v[6488] + v[2781] * v[6547] + v[2787] * v[6548] - v[2796] * v[6706]
			+ v[3488] * v[6708] + (-v[3311] - v[6709])*v[678] + (-v[3243] - v[6710])*v[684] + v[391] * (v[3553] * v[678]
				+ v[3531] * v[684]) - v[6707] * v[690] + v[1781] * (v[1022] * v[2796] + v[3029] * v[690]) - v[3610] * v[880]
			- v[3609] * v[898] + v[3612] * v[916] + v[3611] * v[934];
		v[4320] = v[3987] * v[6446];
		v[3988] = v[3263] * v[3291] + v[1976] * v[3574] + v[1897] * v[3575] + v[1824] * v[3576] + v[33] * (-(v[1654] * v[3619])
			- v[1655] * v[3620] - v[1656] * v[3621]) + v[3363] * v[6490] + v[2782] * v[6547] + v[2788] * v[6548] - v[2797] * v[6706]
			+ v[3490] * v[6708] + (-v[3311] - v[6709])*v[677] + (-v[3243] - v[6710])*v[683] + v[391] * (v[3553] * v[677]
				+ v[3531] * v[683]) - v[6707] * v[689] + v[1781] * (v[1022] * v[2797] + v[3029] * v[689]) - v[3610] * v[879]
			- v[3609] * v[897] + v[3612] * v[915] + v[3611] * v[933];
		v[6792] = v[324] * v[3988];
		v[4324] = v[3988] * v[6445];
		v[3989] = v[3263] * v[3293] + v[1978] * v[3574] + v[1899] * v[3575] + v[1826] * v[3576] + v[33] * (-(v[1650] * v[3619])
			- v[1651] * v[3620] - v[1652] * v[3621]) + v[3363] * v[6492] + v[2783] * v[6547] + v[2789] * v[6548] - v[2800] * v[6706]
			+ v[3492] * v[6708] + (-v[3311] - v[6709])*v[676] + (-v[3243] - v[6710])*v[682] + v[391] * (v[3553] * v[676]
				+ v[3531] * v[682]) - v[6707] * v[688] + v[1781] * (v[1022] * v[2800] + v[3029] * v[688]) - v[3610] * v[878]
			- v[3609] * v[896] + v[3612] * v[914] + v[3611] * v[932];
		v[6793] = v[324] * v[3989];
		v[4327] = v[3989] * v[6444];
		v[3990] = v[3263] * v[3295] + v[1980] * v[3574] + v[1901] * v[3575] + v[1828] * v[3576] + v[33] * (-(v[1634] * v[3619])
			- v[1635] * v[3620] - v[1636] * v[3621]) + v[3363] * v[6494] + v[2778] * v[6547] + v[2784] * v[6548] - v[2790] * v[6706]
			+ v[3494] * v[6708] + (-v[3311] - v[6709])*v[675] + (-v[3243] - v[6710])*v[681] + v[391] * (v[3553] * v[675]
				+ v[3531] * v[681]) - v[6707] * v[687] + v[1781] * (v[1022] * v[2790] + v[3029] * v[687]) - v[3610] * v[874]
			- v[3609] * v[892] + v[3612] * v[910] + v[3611] * v[928];
		v[4350] = v[3990] * v[6439];
		v[3991] = v[3263] * v[3297] + v[1982] * v[3574] + v[1903] * v[3575] + v[1830] * v[3576] + v[33] * (-(v[1630] * v[3619])
			- v[1631] * v[3620] - v[1632] * v[3621]) + v[3363] * v[6496] + v[2779] * v[6547] + v[2785] * v[6548] - v[2791] * v[6706]
			+ v[3496] * v[6708] + (-v[3311] - v[6709])*v[674] + (-v[3243] - v[6710])*v[680] + v[391] * (v[3553] * v[674]
				+ v[3531] * v[680]) - v[6707] * v[686] + v[1781] * (v[1022] * v[2791] + v[3029] * v[686]) - v[3610] * v[873]
			- v[3609] * v[891] + v[3612] * v[909] + v[3611] * v[927];
		v[6797] = v[298] * v[3991];
		v[4354] = v[3991] * v[6438];
		v[3992] = v[3263] * v[3299] + v[1984] * v[3574] + v[1905] * v[3575] + v[1832] * v[3576] + v[33] * (-(v[1626] * v[3619])
			- v[1627] * v[3620] - v[1628] * v[3621]) + v[3363] * v[6498] + v[2780] * v[6547] + v[2786] * v[6548] - v[2794] * v[6706]
			+ v[3498] * v[6708] + (-v[3311] - v[6709])*v[673] + (-v[3243] - v[6710])*v[679] + v[391] * (v[3553] * v[673]
				+ v[3531] * v[679]) - v[6707] * v[685] + v[1781] * (v[1022] * v[2794] + v[3029] * v[685]) - v[3610] * v[872]
			- v[3609] * v[890] + v[3612] * v[908] + v[3611] * v[926];
		v[6798] = v[298] * v[3992];
		v[4357] = v[3992] * v[6437];
		v[3889] = v[3889] + v[3888] * v[673] + v[3887] * v[674] + v[3886] * v[675] + v[3885] * v[676] + v[3884] * v[677]
			+ v[3883] * v[678] - v[3882] * v[781] - v[3881] * v[782] - v[3880] * v[783];
		v[3041] = v[3041] + v[3283] * v[3880] + v[3285] * v[3881] + v[3287] * v[3882] + v[3289] * v[3883] + v[3291] * v[3884]
			+ v[3293] * v[3885] + v[3295] * v[3886] + v[3297] * v[3887] + v[3299] * v[3888] + v[3889] * v[6686];
		v[3827] = v[3827] + v[3824] * v[685] + v[3823] * v[686] + v[3822] * v[687] + v[3821] * v[688] + v[3820] * v[689]
			+ v[3819] * v[690] - v[3818] * v[787] - v[3817] * v[788] - v[3816] * v[789];
		v[3993] = -(v[1014] * v[3880]) + v[372] * (v[3436] + v[6752]) - v[389] * (v[3457] + v[6763] + v[6781]);
		v[3994] = -(v[1014] * v[3881]) + v[362] * (v[3436] + v[6752]) - v[389] * (v[3456] + v[6762] + v[6779]);
		v[3995] = -(v[1014] * v[3882]) + v[352] * (v[3436] + v[6752]) - v[389] * (v[3455] + v[6761] + v[6777]);
		v[3996] = v[1014] * v[3883] + v[346] * (-v[3436] - v[6752]) + v[389] * (v[3454] + v[6760] + v[6775]);
		v[3997] = v[1014] * v[3884] + v[336] * (-v[3436] - v[6752]) + v[389] * (v[3453] + v[6759] + v[6773]);
		v[3998] = v[1014] * v[3885] + v[326] * (-v[3436] - v[6752]) + v[389] * (v[3452] + v[6758] + v[6771]);
		v[3999] = v[1014] * v[3886] + v[320] * (-v[3436] - v[6752]) + v[389] * (v[3451] + v[6757] + v[6769]);
		v[4000] = v[1014] * v[3887] + v[310] * (-v[3436] - v[6752]) + v[389] * (v[3450] + v[6756] + v[6767]);
		v[3041] = v[3041] - v[292] * v[3448] + v[3449] * v[673] + v[3450] * v[674] + v[3451] * v[675] + v[3452] * v[676]
			+ v[3453] * v[677] + v[3454] * v[678] - v[3455] * v[781] - v[3456] * v[782] - v[3457] * v[783] + v[391] * (v[3824] * v[673]
				+ v[3823] * v[674] + v[3822] * v[675] + v[3821] * v[676] + v[3820] * v[677] + v[3819] * v[678] - v[3818] * v[781]
				- v[3817] * v[782] - v[3816] * v[783]) + v[390] * (v[3858] * v[673] + v[3857] * v[674] + v[3856] * v[675] + v[6753]
					+ v[3855] * v[676] + v[3854] * v[677] + v[3853] * v[678] - v[3852] * v[781] - v[3851] * v[782] - v[3850] * v[783]);
		v[4027] = v[3041];
		v[3859] = v[3859] + v[3858] * v[679] + v[3857] * v[680] + v[3856] * v[681] + v[3855] * v[682] + v[3854] * v[683]
			+ v[3853] * v[684] - v[3852] * v[784] - v[3851] * v[785] - v[3850] * v[786];
		v[4001] = v[1014] * v[3888] + v[300] * (-v[3436] - v[6752]) + v[389] * (v[3449] + v[6754] + v[6764]);
		v[3042] = v[3042] + v[3850] * v[6482] + v[3851] * v[6484] + v[3852] * v[6486] + v[3853] * v[6488] + v[3854] * v[6490]
			+ v[3855] * v[6492] + v[3856] * v[6494] + v[3857] * v[6496] + v[3858] * v[6498] + v[3859] * v[6683] + v[391] *
			(v[3824] * v[679] + v[3823] * v[680] + v[3822] * v[681] + v[3821] * v[682] + v[3820] * v[683] + v[3819] * v[684]
				- v[3818] * v[784] - v[3817] * v[785] - v[3816] * v[786]) + v[389] * (v[6753] + v[3888] * v[679] + v[3887] * v[680]
					+ v[3886] * v[681] + v[3885] * v[682] + v[3884] * v[683] + v[3883] * v[684] - v[3882] * v[784] - v[3881] * v[785]
					- v[3880] * v[786]);
		v[4002] = v[1018] * v[3858] + v[300] * v[6755] + v[390] * (v[3463] + v[6754] + v[6765]);
		v[4003] = v[1018] * v[3857] + v[310] * v[6755] + v[390] * (v[3465] + v[6756] + v[6768]);
		v[4004] = v[1018] * v[3856] + v[320] * v[6755] + v[390] * (v[3467] + v[6757] + v[6770]);
		v[4005] = v[1018] * v[3855] + v[326] * v[6755] + v[390] * (v[3469] + v[6758] + v[6772]);
		v[4006] = v[1018] * v[3854] + v[336] * v[6755] + v[390] * (v[3471] + v[6759] + v[6774]);
		v[4007] = v[1018] * v[3853] + v[346] * v[6755] + v[390] * (v[3473] + v[6760] + v[6776]);
		v[4008] = -(v[1018] * v[3852]) - v[352] * v[6755] - v[390] * (v[3475] + v[6761] + v[6778]);
		v[4009] = -(v[1018] * v[3851]) - v[362] * v[6755] - v[390] * (v[3477] + v[6762] + v[6780]);
		v[3042] = v[3042] - v[293] * v[3461] + v[3463] * v[679] + v[3465] * v[680] + v[3467] * v[681] + v[3469] * v[682]
			+ v[3471] * v[683] + v[3473] * v[684] - v[3475] * v[784] - v[3477] * v[785] - v[3479] * v[786];
		v[4025] = v[3042];
		v[4010] = -(v[1018] * v[3850]) - v[372] * v[6755] - v[390] * (v[3479] + v[6763] + v[6782]);
		v[4011] = v[3373] + v[6736] - v[6737];
		v[6819] = -(v[391] * v[4011]);
		v[6820] = v[6766] + v[6817] - v[6819];
		v[3033] = v[3033] - v[294] * v[3373] + v[3482] * v[3816] + v[3484] * v[3817] + v[3486] * v[3818] + v[3488] * v[3819]
			+ v[3490] * v[3820] + v[3492] * v[3821] + v[3494] * v[3822] + v[3496] * v[3823] + v[3498] * v[3824] + v[4011] * v[5461]
			+ v[3827] * v[6680] + v[390] * (v[3858] * v[685] + v[3857] * v[686] + v[3856] * v[687] + v[3855] * v[688]
				+ v[3854] * v[689] + v[3853] * v[690] - v[3852] * v[787] - v[3851] * v[788] - v[3850] * v[789]) + v[389] *
				(v[3888] * v[685] + v[3887] * v[686] + v[3886] * v[687] + v[3885] * v[688] + v[3884] * v[689] + v[3883] * v[690]
					- v[3882] * v[787] - v[3881] * v[788] - v[3880] * v[789]);
		v[4012] = v[1022] * v[3824] + v[391] * (v[3374] + v[6764] + v[6765]) + v[300] * (-v[6707] + v[6766]);
		v[4013] = v[1022] * v[3823] + v[310] * (-v[6707] + v[6766]) + v[391] * (v[3375] + v[6767] + v[6768]);
		v[4014] = v[1022] * v[3822] + v[320] * (-v[6707] + v[6766]) + v[391] * (v[3376] + v[6769] + v[6770]);
		v[4015] = v[1022] * v[3821] + v[326] * (-v[6707] + v[6766]) + v[391] * (v[3377] + v[6771] + v[6772]);
		v[4016] = v[1022] * v[3820] + v[336] * (-v[6707] + v[6766]) + v[391] * (v[3378] + v[6773] + v[6774]);
		v[4017] = v[1022] * v[3819] + v[346] * (-v[6707] + v[6766]) + v[391] * (v[3379] + v[6775] + v[6776]);
		v[4018] = -(v[1022] * v[3818]) + v[352] * (v[6707] - v[6766]) + v[391] * (-v[3380] - v[6777] - v[6778]);
		v[4019] = -(v[1022] * v[3817]) + v[362] * (v[6707] - v[6766]) + v[391] * (-v[3381] - v[6779] - v[6780]);
		v[3033] = v[3033] + v[3830] * v[5481] + v[3828] * v[5483] + v[3374] * v[685] + v[3375] * v[686] + v[3376] * v[687]
			+ v[3377] * v[688] + v[3378] * v[689] + v[3379] * v[690] - v[3380] * v[787] - v[3381] * v[788] - v[3382] * v[789];
		v[4023] = v[3033];
		v[4020] = -(v[1022] * v[3816]) + v[372] * (v[6707] - v[6766]) + v[391] * (-v[3382] - v[6781] - v[6782]);
		b4021 = b6;
		if (b4021) {
			v[2999] = v[2999] + v[3033] * v[375];
			v[2997] = v[2997] + v[3033] * v[388];
			v[3033] = 0e0;
			v[2999] = v[2999] + v[3042] * v[374];
			v[2995] = v[2995] + v[3042] * v[388];
			v[3042] = 0e0;
			v[2999] = v[2999] + v[3041] * v[373];
			v[2993] = v[2993] + v[3041] * v[388];
			v[3041] = 0e0;
			v[2998] = v[2998] + v[2999] * v[387];
			v[2992] = v[2992] + v[2998];
		}
		else {
			b4022 = b418;
			if (b4022) {
				v[4024] = -v[4023];
				v[3033] = 0e0;
				v[4026] = -v[4025];
				v[3042] = 0e0;
				v[4028] = -v[4027];
				v[3041] = 0e0;
			}
			else {
				v[4024] = v[4023];
				v[4026] = v[4025];
				v[4028] = v[4027];
			};
			v[3001] = v[3001] + v[375] * v[4024];
			v[2997] = v[4029] + v[4024] * v[407];
			v[3001] = v[3001] + v[374] * v[4026];
			v[2995] = v[4030] + v[4026] * v[407];
			v[3001] = v[3001] + v[373] * v[4028];
			v[2993] = v[4031] + v[4028] * v[407];
			v[3000] = v[3000] + v[3001] * v[406];
			v[2992] = v[3000] + v[4032];
		};
		v[6783] = v[2992] / v[401];
		v[2997] = v[2997] + v[375] * v[6783];
		v[2995] = v[2995] + v[374] * v[6783];
		v[2993] = v[2993] + v[373] * v[6783];
		v[3942] = -v[2997] + v[3942];
		v[3943] = -(v[1027] * v[2997]) + v[3943];
		v[3944] = -(v[274] * v[2997]) + v[3944];
		v[3945] = -(v[1613] * v[2997]) + v[3945];
		v[3932] = -v[2995] + v[3932];
		v[3933] = -(v[1027] * v[2995]) + v[3933];
		v[3934] = -(v[274] * v[2995]) + v[3934];
		v[3935] = -(v[1613] * v[2995]) + v[3935];
		v[3922] = -v[2993] + v[3922];
		v[3923] = -(v[1027] * v[2993]) + v[3923];
		v[3924] = -(v[274] * v[2993]) + v[3924];
		v[3925] = -(v[1613] * v[2993]) + v[3925];
		v[4057] = v[359] * v[3984] + v[1969] * v[6784];
		v[4058] = v[353] * v[3984] + v[1969] * v[6785];
		v[4059] = v[2173] * v[3984] + v[1969] * (v[2864] / v[350] + v[2843] * v[5972]) + v[4057] * v[7105];
		v[4060] = v[2179] * v[3984] + v[1969] * (v[2856] / v[350] + v[2843] * v[5974]) + v[4058] * v[7106];
		v[4061] = (v[252] * v[4057] + v[250] * v[4058] + v[1969] * (v[2872] + v[5976] * v[6672]) + v[3984] * v[7107]) / v[350];
		v[4062] = v[361] * v[3985] + v[1971] * v[6786];
		v[4063] = v[353] * v[3985] + v[1971] * v[6785];
		v[4066] = v[2168] * v[3985] + v[1971] * (v[2870] / v[350] + v[2843] * v[5980]) + v[4062] * v[7109];
		v[4067] = v[4057] + v[4062];
		v[4068] = v[4059] + v[4066];
		v[4070] = (v[2171] * v[2829] + v[241] * v[4063] + v[243] * v[4067] + v[4428] * v[6672] + v[1971] * (v[2852]
			+ v[5984] * v[6672]) + v[2178] * v[6787]) / v[350];
		v[4073] = (v[248] * v[4062] + v[245] * v[4063] + v[1971] * (v[2861] + v[5986] * v[6672]) + v[2172] * v[6787]) / v[350];
		v[4074] = v[359] * v[3986] + v[1973] * v[6784];
		v[4075] = v[361] * v[3986] + v[1973] * v[6786];
		v[4076] = v[4063] + v[4074];
		v[4079] = (v[242] * v[4074] + v[243] * v[4075] + v[1973] * (v[2853] + v[5991] * v[6672]) + v[2176] * v[6788]) / v[350];
		v[4080] = v[4058] + v[4075];
		v[4082] = (v[2181] * v[2823] + v[247] * v[4074] + v[248] * v[4080] + v[4427] * v[6672] + v[1973] * (v[2860]
			+ v[5994] * v[6672]) + v[2183] * v[6788]) / v[350];
		v[4083] = v[4070] + v[4082];
		v[4085] = (v[2182] * v[2818] + v[253] * v[4075] + v[252] * v[4076] + v[4426] * v[6672] + v[1973] * (v[2868]
			+ v[5997] * v[6672]) + v[2186] * v[6788]) / v[350];
		v[4086] = v[4060] + v[4085];
		v[4087] = v[333] * v[3987] + v[1975] * v[6789];
		v[4088] = v[327] * v[3987] + v[1975] * v[6790];
		v[4089] = v[2197] * v[3987] + v[1975] * (v[2754] / v[324] + v[2721] * v[6002]) + v[4087] * v[7114];
		v[4090] = v[2203] * v[3987] + v[1975] * (v[2749] / v[324] + v[2721] * v[6004]) + v[4088] * v[7115];
		v[4091] = (v[192] * v[4087] + v[190] * v[4088] + v[1975] * (v[2759] + v[6006] * v[6671]) + v[3987] * v[7116]) / v[324];
		v[4092] = v[335] * v[3988] + v[1977] * v[6791];
		v[4093] = v[327] * v[3988] + v[1977] * v[6790];
		v[4096] = v[2192] * v[3988] + v[1977] * (v[2757] / v[324] + v[2721] * v[6010]) + v[4092] * v[7118];
		v[4097] = v[4087] + v[4092];
		v[4098] = v[4089] + v[4096];
		v[4100] = (v[2195] * v[2707] + v[181] * v[4093] + v[183] * v[4097] + v[4456] * v[6671] + v[1977] * (v[2745]
			+ v[6014] * v[6671]) + v[2202] * v[6792]) / v[324];
		v[4103] = (v[188] * v[4092] + v[185] * v[4093] + v[1977] * (v[2751] + v[6016] * v[6671]) + v[2196] * v[6792]) / v[324];
		v[4104] = v[333] * v[3989] + v[1979] * v[6789];
		v[4105] = v[335] * v[3989] + v[1979] * v[6791];
		v[4106] = v[4093] + v[4104];
		v[4109] = (v[182] * v[4104] + v[183] * v[4105] + v[1979] * (v[2746] + v[6021] * v[6671]) + v[2200] * v[6793]) / v[324];
		v[4110] = v[4088] + v[4105];
		v[4112] = (v[2205] * v[2701] + v[187] * v[4104] + v[188] * v[4110] + v[4455] * v[6671] + v[1979] * (v[2750]
			+ v[6024] * v[6671]) + v[2207] * v[6793]) / v[324];
		v[4113] = v[4100] + v[4112];
		v[4115] = (v[2206] * v[2696] + v[193] * v[4105] + v[192] * v[4106] + v[4454] * v[6671] + v[1979] * (v[2755]
			+ v[6027] * v[6671]) + v[2210] * v[6793]) / v[324];
		v[4116] = v[4090] + v[4115];
		v[4117] = v[307] * v[3990] + v[1981] * v[6794];
		v[4118] = v[301] * v[3990] + v[1981] * v[6795];
		v[4119] = v[2221] * v[3990] + v[1981] * (v[2739] / v[298] + v[2683] * v[6032]) + v[4117] * v[7123];
		v[4120] = v[2227] * v[3990] + v[1981] * (v[2734] / v[298] + v[2683] * v[6034]) + v[4118] * v[7124];
		v[4121] = (v[173] * v[4117] + v[171] * v[4118] + v[1981] * (v[2744] + v[6036] * v[6670]) + v[3990] * v[7125]) / v[298];
		v[4122] = v[309] * v[3991] + v[1983] * v[6796];
		v[4123] = v[301] * v[3991] + v[1983] * v[6795];
		v[4126] = v[2216] * v[3991] + v[1983] * (v[2742] / v[298] + v[2683] * v[6040]) + v[4122] * v[7127];
		v[4127] = v[4117] + v[4122];
		v[4128] = v[4119] + v[4126];
		v[4130] = (v[2219] * v[2669] + v[162] * v[4123] + v[164] * v[4127] + v[4483] * v[6670] + v[1983] * (v[2730]
			+ v[6044] * v[6670]) + v[2226] * v[6797]) / v[298];
		v[4133] = (v[169] * v[4122] + v[166] * v[4123] + v[1983] * (v[2736] + v[6046] * v[6670]) + v[2220] * v[6797]) / v[298];
		v[4134] = v[307] * v[3992] + v[1985] * v[6794];
		v[4135] = v[309] * v[3992] + v[1985] * v[6796];
		v[4136] = v[4123] + v[4134];
		v[4139] = (v[163] * v[4134] + v[164] * v[4135] + v[1985] * (v[2731] + v[6051] * v[6670]) + v[2224] * v[6798]) / v[298];
		v[4140] = v[4118] + v[4135];
		v[4142] = (v[2229] * v[2663] + v[168] * v[4134] + v[169] * v[4140] + v[4482] * v[6670] + v[1985] * (v[2735]
			+ v[6054] * v[6670]) + v[2231] * v[6798]) / v[298];
		v[4143] = v[4130] + v[4142];
		v[4145] = (v[2230] * v[2658] + v[174] * v[4135] + v[173] * v[4136] + v[4481] * v[6670] + v[1985] * (v[2740]
			+ v[6057] * v[6670]) + v[2234] * v[6798]) / v[298];
		v[4146] = v[4120] + v[4145];
		v[4147] = -(v[3995] * v[982]);
		v[4148] = -(v[3995] * v[981]);
		v[4149] = -(v[3995] * v[983]);
		v[4150] = -(v[3994] * v[983]);
		v[4151] = -(v[3994] * v[982]);
		v[4152] = -(v[3994] * v[981]);
		v[4153] = -(v[3993] * v[983]);
		v[4154] = -(v[3993] * v[981]);
		v[4155] = -(v[3993] * v[982]);
		v[4156] = -(v[4008] * v[983]);
		v[4157] = -(v[4008] * v[982]);
		v[4158] = -(v[4008] * v[981]);
		v[4159] = -(v[4009] * v[983]);
		v[4160] = -(v[4009] * v[981]);
		v[4161] = -(v[4009] * v[982]);
		v[4162] = -(v[4010] * v[981]);
		v[4163] = -(v[4010] * v[983]);
		v[4164] = -(v[4010] * v[982]);
		v[4165] = -(v[4018] * v[983]);
		v[4166] = -(v[4018] * v[982]);
		v[4167] = -(v[4018] * v[981]);
		v[4168] = -(v[4019] * v[982]);
		v[4169] = -(v[4019] * v[983]);
		v[4170] = -(v[4019] * v[981]);
		v[4171] = -(v[4020] * v[982]);
		v[4172] = -(v[4020] * v[983]);
		v[4173] = -(v[4020] * v[981]);
		v[3943] = v[3943] - v[1613] * v[3980] + v[3981] * v[451];
		v[3944] = v[3944] + v[1612] * v[3981];
		v[3945] = v[3945] + v[1027] * v[3980] + v[3981] * v[449];
		v[3933] = v[3933] - v[1613] * v[3976] + v[3977] * v[451];
		v[3934] = v[3934] + v[1612] * v[3977];
		v[3935] = v[3935] + v[1027] * v[3976] + v[3977] * v[449];
		v[4204] = -(v[130] * v[2887]) - v[136] * v[2888] - v[133] * v[2890] + v[2988] - v[259] * v[2993] - v[262] * v[2995]
			- v[265] * v[2997] + v[257] * v[3972] + v[260] * v[3976] + v[263] * v[3980] + v[3995] * v[742] + v[3994] * v[743]
			+ v[3993] * v[744] + v[4008] * v[757] + v[4009] * v[758] + v[4010] * v[759] + v[4018] * v[772] + v[4019] * v[773]
			+ v[4020] * v[774];
		v[4205] = -(v[128] * v[2887]) - v[134] * v[2888] - v[131] * v[2890] + v[2986] - v[257] * v[2993] - v[260] * v[2995]
			- v[263] * v[2997] - v[259] * v[3972] - v[262] * v[3976] - v[265] * v[3980] + v[3995] * v[736] + v[3994] * v[737]
			+ v[3993] * v[738] + v[4008] * v[751] + v[4009] * v[752] + v[4010] * v[753] + v[4018] * v[766] + v[4019] * v[767]
			+ v[4020] * v[768];
		v[3923] = v[3923] - v[1613] * v[3972] + v[3973] * v[451];
		v[3924] = v[3924] + v[1612] * v[3973];
		v[3925] = v[3925] + v[1027] * v[3972] + v[3973] * v[449];
		v[4206] = v[1769] * v[2857] + v[1744] * v[2865] + v[1719] * v[2873] + v[257] * v[3973] + v[260] * v[3977]
			+ v[263] * v[3981];
		v[4207] = v[1769] * v[2859] + v[1744] * v[2867] + v[1719] * v[2875] + v[259] * v[3973] + v[262] * v[3977]
			+ v[265] * v[3981];
		v[4220] = v[136] * v[3943] + v[135] * v[3944] + v[134] * v[3945] + v[2519] * v[4062] + v[363] * v[4189]
			+ v[4075] * v[6476] + v[361] * (v[6920] + v[2843] * v[7351]);
		v[4221] = v[133] * v[3943] + v[132] * v[3944] + v[131] * v[3945] + (v[2182] * v[2947] + v[2964] * v[359]
			+ v[363] * v[4057] + v[348] * v[4076]) / v[350] + v[360] * v[4193] + v[2843] * v[4423] + v[4178] * v[6801]
			+ v[4179] * v[6801];
		v[4222] = v[130] * v[3943] + v[129] * v[3944] + v[128] * v[3945] + v[348] * v[4196] + v[4058] * v[4516] + v[353] *
			(v[2843] * v[6802] + v[6911]);
		v[4223] = v[136] * v[3933] + v[135] * v[3934] + v[134] * v[3935] + (v[2181] * v[2951] + v[2972] * v[361]
			+ v[354] * v[4062] + v[349] * v[4080]) / v[350] + v[2843] * v[4424] + v[4189] * v[6477] + (v[4191] + v[4192])*v[6803];
		v[4224] = v[133] * v[3933] + v[132] * v[3934] + v[131] * v[3935] + v[2521] * v[4057] + v[354] * v[4193]
			+ v[4074] * v[6473] + v[359] * (v[6917] + v[2843] * v[7352]);
		v[4225] = v[130] * v[3933] + v[129] * v[3934] + v[128] * v[3935] + v[349] * v[4196] + v[4063] * v[4514] + v[353] *
			(v[2843] * v[6804] + v[6910]);
		v[4226] = v[136] * v[3923] + v[135] * v[3924] + v[134] * v[3925] + v[2843] * v[4425] + (v[2171] * v[2952]
			+ v[2968] * v[361] + v[351] * v[4075] + v[4067] * v[6472]) / v[350] + v[4189] * v[6475] + (v[4209] + v[4211])*v[6803];
		v[4227] = v[133] * v[3923] + v[132] * v[3924] + v[131] * v[3925] + v[4074] * v[4511] + v[4193] * v[6472] + v[359] *
			(v[2843] * v[6805] + v[6916]);
		v[4228] = v[130] * v[3923] + v[129] * v[3924] + v[128] * v[3925] + v[2522] * v[4058] + v[2520] * v[4063]
			+ v[351] * v[4196] + v[353] * (v[6912] + v[2843] * v[7353]);
		v[4229] = -(v[2175] * v[2645]) - v[2185] * v[2646] + 2e0*v[2174] * v[2647] + v[4147] + v[4156] + v[4162] + v[4171]
			+ v[4073] * v[6669] - v[4083] * v[701] - v[4068] * v[704];
		v[4230] = v[2188] * v[2645] + 2e0*v[2180] * v[2646] - v[2185] * v[2647] - v[4151] - v[4154] - v[4159] - v[4172]
			+ v[4079] * v[6668] - v[4083] * v[698] + v[4086] * v[704];
		v[4231] = -(v[2245] * v[2635]) + v[2241] * v[2638] - v[2237] * v[2640] + v[2305] * v[2821] + v[238] * v[4227]
			- v[4147] * v[691] + v[4151] * v[692] - v[4155] * v[693];
		v[4232] = 2e0*v[2167] * v[2645] + v[2188] * v[2646] - v[2175] * v[2647] - v[4148] - v[4160] - v[4165] - v[4168]
			+ v[4061] * v[6667] - v[4068] * v[698] + v[4086] * v[701];
		v[4233] = -(v[2244] * v[2635]) + v[2242] * v[2638] - v[2238] * v[2640] + v[2304] * v[2821] + v[238] * v[4226]
			- v[4148] * v[691] + v[4152] * v[692] - v[4154] * v[693];
		v[4234] = -(v[2253] * v[2635]) + v[2249] * v[2638] - v[2246] * v[2640] + v[2290] * v[2821] + v[238] * v[4225]
			- v[4156] * v[691] + v[4159] * v[692] - v[4163] * v[693];
		v[4235] = v[4153] + v[4164];
		v[4236] = -(v[2252] * v[2635]) + v[2250] * v[2638] - v[2248] * v[2640] + v[2288] * v[2821] + v[238] * v[4223]
			- v[4158] * v[691] + v[4160] * v[692] - v[4162] * v[693];
		v[4237] = -(v[2262] * v[2635]) + v[2259] * v[2638] - v[2255] * v[2640] + v[2278] * v[2821] + v[238] * v[4222]
			- v[4165] * v[691] + v[4169] * v[692] - v[4172] * v[693];
		v[4238] = -(v[2261] * v[2635]) + v[2258] * v[2638] - v[2256] * v[2640] + v[2277] * v[2821] + v[238] * v[4221]
			- v[4166] * v[691] + v[4168] * v[692] - v[4171] * v[693];
		v[4239] = v[4157] + v[4167];
		v[4240] = v[4150] + v[4170];
		v[4241] = v[223] * v[2993] + v[4242];
		v[4243] = v[224] * v[2993] - v[4242];
		v[3927] = v[3927] + v[4241];
		v[3926] = v[3926] + v[4243];
		v[4244] = v[223] * v[2995] + v[4245];
		v[4246] = v[224] * v[2995] - v[4245];
		v[3937] = v[3937] + v[4244];
		v[3936] = v[3936] + v[4246];
		v[4247] = v[223] * v[2997] + v[4248];
		v[4249] = v[224] * v[2997] - v[4248];
		v[3947] = v[3947] + v[4247];
		v[3946] = v[3946] + v[4249];
		v[4250] = v[4015] * v[986];
		v[4251] = v[4015] * v[985];
		v[4252] = v[4015] * v[984];
		v[4253] = v[4016] * v[985];
		v[4254] = v[4016] * v[986];
		v[4255] = v[4016] * v[984];
		v[4256] = v[4017] * v[985];
		v[4257] = v[4017] * v[986];
		v[4258] = v[4017] * v[984];
		v[4259] = v[4012] * v[989];
		v[4260] = v[4012] * v[988];
		v[4261] = v[4012] * v[987];
		v[4262] = v[4013] * v[988];
		v[4263] = v[4013] * v[989];
		v[4264] = v[4013] * v[987];
		v[4265] = v[4014] * v[988];
		v[4266] = v[4014] * v[989];
		v[4267] = v[4014] * v[987];
		v[4268] = v[4005] * v[986];
		v[4269] = v[4005] * v[985];
		v[4270] = v[4005] * v[984];
		v[4271] = v[4006] * v[986];
		v[4272] = v[4006] * v[984];
		v[4273] = v[4006] * v[985];
		v[4274] = v[4007] * v[984];
		v[4275] = v[4007] * v[986];
		v[4276] = v[4007] * v[985];
		v[4277] = v[4002] * v[989];
		v[4278] = v[4002] * v[988];
		v[4279] = v[4002] * v[987];
		v[4280] = v[4003] * v[989];
		v[4281] = v[4003] * v[987];
		v[4282] = v[4003] * v[988];
		v[4283] = v[4004] * v[987];
		v[4284] = v[4004] * v[989];
		v[4285] = v[4004] * v[988];
		v[4286] = v[3998] * v[985];
		v[4287] = v[3998] * v[984];
		v[4288] = v[3998] * v[986];
		v[4289] = v[3997] * v[986];
		v[4290] = v[3997] * v[985];
		v[4291] = v[3997] * v[984];
		v[4292] = v[3996] * v[986];
		v[4293] = v[3996] * v[984];
		v[4294] = v[3996] * v[985];
		v[4295] = v[4001] * v[988];
		v[4296] = v[4001] * v[987];
		v[4297] = v[4001] * v[989];
		v[4298] = v[4000] * v[989];
		v[4299] = v[4000] * v[988];
		v[4300] = v[4000] * v[987];
		v[4301] = v[3999] * v[989];
		v[4302] = v[3999] * v[987];
		v[4303] = v[3999] * v[988];
		v[3928] = v[3928] + v[228] * v[4243] + v[1610] * v[6806];
		v[3929] = v[3929] + v[227] * v[4243] - v[1024] * v[6806];
		v[3930] = v[3930] + v[228] * v[4241] + v[1610] * v[6807];
		v[3931] = v[3931] + v[227] * v[4241] - v[1024] * v[6807];
		v[3938] = v[3938] + v[228] * v[4246] + v[1610] * v[6808];
		v[3939] = v[3939] + v[227] * v[4246] - v[1024] * v[6808];
		v[3940] = v[3940] + v[228] * v[4244] + v[1610] * v[6809];
		v[3941] = v[3941] + v[227] * v[4244] - v[1024] * v[6809];
		v[3948] = v[3948] + v[228] * v[4249] + v[1610] * v[6810];
		v[3949] = v[3949] + v[227] * v[4249] - v[1024] * v[6810];
		v[3950] = v[3950] + v[228] * v[4247] + v[1610] * v[6811];
		v[3951] = v[3951] + v[227] * v[4247] - v[1024] * v[6811];
		v[4304] = (-(v[1471] * v[2792]) - v[1470] * v[2793] - v[1469] * v[2795] + v[1467] * v[2798] + v[1466] * v[2799]
			+ v[1465] * v[2801] + v[1781] * (-(v[1022] * v[2762]) - v[288] * v[3029]) + v[1781] * (v[1022] * v[2765]
				+ v[291] * v[3029]) + v[3036] - v[3037] + v[2761] * v[3342] - v[2764] * v[3342] + v[391] * (-(v[286] * v[3828])
					- v[287] * v[3830] - v[288] * v[4011]) + v[391] * (v[289] * v[3828] + v[290] * v[3830] + v[291] * v[4011])
			- v[2993] * v[420] + v[2993] * v[421] - v[2995] * v[423] + v[2995] * v[424] - v[2997] * v[426] + v[2997] * v[427]
			+ v[3983] * v[432] - v[3983] * v[433] + v[3979] * v[434] - v[3979] * v[435] + v[3975] * v[436] - v[3975] * v[437]
			+ v[4015] * v[637] + v[4016] * v[638] + v[4017] * v[639] - v[4012] * v[640] - v[4013] * v[641] - v[4014] * v[642]
			+ v[4005] * v[649] - v[3890] * v[6499] + v[3891] * v[6499] + v[4006] * v[650] + v[4007] * v[651] - v[4002] * v[652]
			- v[4003] * v[653] - v[4004] * v[654] + v[3998] * v[661] + v[3997] * v[662] + v[3996] * v[663] - v[4001] * v[664]
			- v[2760] * v[6644] + v[2763] * v[6644] - v[4000] * v[665] - v[3999] * v[666] + v[1944] * (-v[6684] - v[6685]) + v[1945] *
			(v[6684] + v[6685]) - v[6708] * v[6812] + v[286] * v[6814] - v[289] * v[6814] - v[287] * v[6816] + v[290] * v[6816]
			- v[288] * v[6817] + v[291] * v[6817] + v[6708] * v[6818]) / 2e0;
		v[4371] = v[2511] * v[4092] + v[337] * v[4320] + v[4105] * v[6468] + v[335] * (v[6906] + v[2721] * v[7354])
			+ v[3949] * v[97] + v[3948] * v[98];
		v[4372] = (v[2206] * v[2919]) / v[324] + v[334] * v[4324] + v[4087] * v[4509] + v[2936] * v[6445] + v[4106] * v[6468]
			+ v[2721] * v[7355] + v[3949] * v[94] + v[3948] * v[95];
		v[4373] = v[322] * v[4327] + v[4088] * v[4509] + v[327] * (v[6897] + v[2721] * v[7356]) + v[3949] * v[91]
			+ v[3948] * v[92];
		v[4374] = (v[2205] * v[2923]) / v[324] + v[4092] * v[4507] + v[2944] * v[6446] + v[4110] * v[6465] + v[4320] * v[6469]
			+ v[2721] * v[7357] + v[3939] * v[97] + v[3938] * v[98];
		v[4375] = v[2513] * v[4087] + v[328] * v[4324] + v[4104] * v[6465] + v[333] * (v[6903] + v[2721] * v[7358])
			+ v[3939] * v[94] + v[3938] * v[95];
		v[4376] = v[323] * v[4327] + v[4093] * v[4507] + v[327] * (v[6896] + v[2721] * v[7359]) + v[3939] * v[91]
			+ v[3938] * v[92];
		v[4377] = (v[2195] * v[2924]) / v[324] + v[2512] * v[4097] + v[4105] * v[4504] + v[2940] * v[6446] + v[4320] * v[6467]
			+ v[2721] * v[7360] + v[3929] * v[97] + v[3928] * v[98];
		v[4378] = v[4104] * v[4504] + v[4324] * v[6464] + v[333] * (v[2721] * v[6825] + v[6902]) + v[3929] * v[94]
			+ v[3928] * v[95];
		v[4379] = v[2514] * v[4088] + v[2512] * v[4093] + v[325] * v[4327] + v[327] * (v[6898] + v[2721] * v[7361])
			+ v[3929] * v[91] + v[3928] * v[92];
		v[4380] = -(v[2199] * v[2619]) - v[2209] * v[2620] + 2e0*v[2198] * v[2621] + v[4256] + v[4268] + v[4274] + v[4286]
			- v[4113] * v[530] - v[4098] * v[533] + v[4103] * v[6659];
		v[4381] = v[2212] * v[2619] + 2e0*v[2204] * v[2620] - v[2209] * v[2621] - v[4257] - v[4271] - v[4290] - v[4293]
			- v[4113] * v[527] + v[4116] * v[533] + v[4109] * v[6658];
		v[4382] = -(v[2374] * v[2609]) + v[2370] * v[2612] - v[2366] * v[2614] + v[2415] * v[2699] + v[178] * v[4378]
			- v[4286] * v[520] + v[4290] * v[521] - v[4294] * v[522];
		v[4383] = 2e0*v[2191] * v[2619] + v[2212] * v[2620] - v[2199] * v[2621] - v[4250] - v[4253] - v[4272] - v[4287]
			- v[4098] * v[527] + v[4116] * v[530] + v[4091] * v[6657];
		v[4384] = -(v[2373] * v[2609]) + v[2371] * v[2612] - v[2367] * v[2614] + v[2414] * v[2699] + v[178] * v[4377]
			- v[4287] * v[520] + v[4291] * v[521] - v[4293] * v[522];
		v[4385] = -(v[2355] * v[2609]) + v[2351] * v[2612] - v[2348] * v[2614] + v[2407] * v[2699] + v[178] * v[4376]
			- v[4268] * v[520] + v[4271] * v[521] - v[4275] * v[522];
		v[4386] = v[4276] + v[4292];
		v[4387] = -(v[2354] * v[2609]) + v[2352] * v[2612] - v[2350] * v[2614] + v[2405] * v[2699] + v[178] * v[4374]
			- v[4270] * v[520] + v[4272] * v[521] - v[4274] * v[522];
		v[4388] = -(v[2337] * v[2609]) + v[2334] * v[2612] - v[2330] * v[2614] + v[2398] * v[2699] + v[178] * v[4373]
			- v[4250] * v[520] + v[4254] * v[521] - v[4257] * v[522];
		v[4389] = -(v[2336] * v[2609]) + v[2333] * v[2612] - v[2331] * v[2614] + v[2397] * v[2699] + v[178] * v[4372]
			- v[4251] * v[520] + v[4253] * v[521] - v[4256] * v[522];
		v[4390] = v[4252] + v[4269];
		v[4391] = v[4255] + v[4289];
		v[4392] = v[2503] * v[4122] + v[311] * v[4350] + v[4135] * v[6460] + v[309] * (v[6892] + v[2683] * v[7362])
			+ v[3951] * v[88] + v[3950] * v[89];
		v[4393] = (v[2230] * v[2891]) / v[298] + v[308] * v[4354] + v[4117] * v[4502] + v[2908] * v[6438] + v[4136] * v[6460]
			+ v[2683] * v[7363] + v[3951] * v[85] + v[3950] * v[86];
		v[4394] = v[296] * v[4357] + v[4118] * v[4502] + v[301] * (v[6883] + v[2683] * v[7364]) + v[3951] * v[82]
			+ v[3950] * v[83];
		v[4395] = (v[2229] * v[2895]) / v[298] + v[4122] * v[4500] + v[2916] * v[6439] + v[4140] * v[6457] + v[4350] * v[6461]
			+ v[2683] * v[7365] + v[3941] * v[88] + v[3940] * v[89];
		v[4396] = v[2505] * v[4117] + v[302] * v[4354] + v[4134] * v[6457] + v[307] * (v[6889] + v[2683] * v[7366])
			+ v[3941] * v[85] + v[3940] * v[86];
		v[4397] = v[297] * v[4357] + v[4123] * v[4500] + v[301] * (v[6882] + v[2683] * v[7367]) + v[3941] * v[82]
			+ v[3940] * v[83];
		v[4398] = (v[2219] * v[2896]) / v[298] + v[2504] * v[4127] + v[4135] * v[4497] + v[2912] * v[6439] + v[4350] * v[6459]
			+ v[2683] * v[7368] + v[3931] * v[88] + v[3930] * v[89];
		v[4399] = v[4134] * v[4497] + v[4354] * v[6456] + v[307] * (v[2683] * v[6828] + v[6888]) + v[3931] * v[85]
			+ v[3930] * v[86];
		v[4400] = v[2506] * v[4118] + v[2504] * v[4123] + v[299] * v[4357] + v[301] * (v[6884] + v[2683] * v[7369])
			+ v[3931] * v[82] + v[3930] * v[83];
		v[4401] = -(v[2223] * v[2593]) - v[2233] * v[2594] + 2e0*v[2222] * v[2595] + v[4265] + v[4277] + v[4283] + v[4295]
			- v[4143] * v[485] - v[4128] * v[488] + v[4133] * v[6655];
		v[4402] = v[2236] * v[2593] + 2e0*v[2228] * v[2594] - v[2233] * v[2595] - v[4266] - v[4280] - v[4299] - v[4302]
			- v[4143] * v[482] + v[4146] * v[488] + v[4139] * v[6654];
		v[4403] = -(v[2383] * v[2583]) + v[2379] * v[2586] - v[2375] * v[2588] + v[2442] * v[2661] + v[159] * v[4399]
			- v[4295] * v[475] + v[4299] * v[476] - v[4303] * v[477];
		v[4404] = 2e0*v[2215] * v[2593] + v[2236] * v[2594] - v[2223] * v[2595] - v[4259] - v[4262] - v[4281] - v[4296]
			- v[4128] * v[482] + v[4146] * v[485] + v[4121] * v[6653];
		v[4405] = -(v[2382] * v[2583]) + v[2380] * v[2586] - v[2376] * v[2588] + v[2441] * v[2661] + v[159] * v[4398]
			- v[4296] * v[475] + v[4300] * v[476] - v[4302] * v[477];
		v[4406] = -(v[2364] * v[2583]) + v[2360] * v[2586] - v[2357] * v[2588] + v[2434] * v[2661] + v[159] * v[4397]
			- v[4277] * v[475] + v[4280] * v[476] - v[4284] * v[477];
		v[4407] = v[4285] + v[4301];
		v[4408] = -(v[2363] * v[2583]) + v[2361] * v[2586] - v[2359] * v[2588] + v[2432] * v[2661] + v[159] * v[4395]
			- v[4279] * v[475] + v[4281] * v[476] - v[4283] * v[477];
		v[4409] = -(v[2346] * v[2583]) + v[2343] * v[2586] - v[2339] * v[2588] + v[2425] * v[2661] + v[159] * v[4394]
			- v[4259] * v[475] + v[4263] * v[476] - v[4266] * v[477];
		v[4410] = -(v[2345] * v[2583]) + v[2342] * v[2586] - v[2340] * v[2588] + v[2424] * v[2661] + v[159] * v[4393]
			- v[4260] * v[475] + v[4262] * v[476] - v[4265] * v[477];
		v[4411] = v[4261] + v[4278];
		v[4412] = v[4264] + v[4298];
		v[4413] = -(v[1026] * (v[1769] * v[2858] + v[1744] * v[2866] + v[1719] * v[2874] + v[258] * v[3973] + v[261] * v[3977]
			+ v[264] * v[3981] - v[1611] * v[4204] + v[1025] * v[4205])) + v[1612] * (-(v[129] * v[2887]) - v[135] * v[2888]
				- v[132] * v[2890] + v[2987] - v[258] * v[2993] - v[261] * v[2995] - v[264] * v[2997] - v[1025] * v[4206]
				+ v[1611] * v[4207] + v[3995] * v[739] + v[3994] * v[740] + v[3993] * v[741] + v[4008] * v[754] + v[4009] * v[755]
				+ v[4010] * v[756] + v[4018] * v[769] + v[4019] * v[770] + v[4020] * v[771]);
		v[4414] = -(v[275] * (v[1025] * v[4204] + v[1611] * v[4205])) + v[1026] * (v[1611] * v[4206] + v[1025] * v[4207]);
		v[4415] = (-(v[2257] * v[2554]) - v[2239] * v[2557] - v[2247] * v[2636] - 2e0*v[4059] + 2e0*v[4066]
			+ v[4148] * v[6626] + v[4147] * v[6627] + v[4158] * v[6628] + v[4156] * v[6629] + v[4166] * v[6630] + v[4165] * v[6631]
			- v[2246] * v[6829] - v[2237] * v[6830] - v[2256] * v[6831] - v[2248] * v[6832] - v[2255] * v[6833] - v[2238] * v[6834]
			- v[4149] * v[694] - v[4157] * v[713] - v[4167] * v[731]) / 2e0;
		v[4416] = (-(v[2243] * v[2635]) + v[2240] * v[2638] - v[2239] * v[2640] + v[2306] * v[2821] + v[238] * v[4228]
			- v[4149] * v[691] + v[4150] * v[692] - v[4153] * v[693]) / 2e0;
		v[4418] = (v[2260] * v[2554] + v[2240] * v[2557] + v[2251] * v[2636] - 2e0*v[4060] + 2e0*v[4085] - v[4152] * v[6626]
			- v[4151] * v[6627] - v[4160] * v[6628] - v[4159] * v[6629] - v[4168] * v[6630] - v[4169] * v[6631] + v[2249] * v[6829]
			+ v[2241] * v[6830] + v[2258] * v[6831] + v[2250] * v[6832] + v[2259] * v[6833] + v[2242] * v[6834] + v[4150] * v[694]
			+ v[4161] * v[713] + v[4170] * v[731]) / 2e0;
		v[4419] = (-(v[2254] * v[2635]) + v[2251] * v[2638] - v[2247] * v[2640] + v[2289] * v[2821] + v[238] * v[4224]
			- v[4157] * v[691] + v[4161] * v[692] - v[4164] * v[693]) / 2e0;
		v[4420] = (-(v[2263] * v[2554]) - v[2243] * v[2557] - v[2254] * v[2636] - 2e0*v[4070] + 2e0*v[4082]
			+ v[4154] * v[6626] + v[4155] * v[6627] + v[4162] * v[6628] + v[4163] * v[6629] + v[4171] * v[6630] + v[4172] * v[6631]
			- v[2253] * v[6829] - v[2245] * v[6830] - v[2261] * v[6831] - v[2252] * v[6832] - v[2262] * v[6833] - v[2244] * v[6834]
			- v[4153] * v[694] - v[4164] * v[713] - v[4173] * v[731]) / 2e0;
		v[4513] = v[2637] * v[4418] - v[4419] + v[2634] * v[4420] + (v[11454 + i1687] + v[2320] * v[2556] - v[4415] * v[469]
			)*v[6431] + v[2556] * v[7173] * v[7371] - v[6432] * (v[11796 + i1687] + (v[2276] * v[2554]) / 2e0 + (v[2306] * v[2557])
				/ 2e0 + v[2290] * v[2628] + v[2305] * v[2629] + v[2277] * v[2630] + v[2288] * v[2631] + v[2278] * v[2632]
				+ v[2304] * v[2633] + (v[2289] * v[2636]) / 2e0 + v[4152] - v[4155] - v[4158] + v[4163] + v[4166] - v[4169]
				+ 2e0*v[2821] * v[4422] + v[4229] * v[470] - v[4230] * v[471] - v[4232] * v[473] + v[4220] * v[6575] + v[4224] * v[6576]
				+ v[4228] * v[6577] + v[4235] * v[6837] + v[4239] * v[6838] + v[4240] * v[6839] + v[4227] * v[699] + v[4226] * v[705]
				+ v[4225] * v[709] + v[4223] * v[718] + v[4222] * v[722] + v[4221] * v[726] + v[6648] * (v[4061] + v[2964] * v[4064]
					+ v[2961] * v[4065] + v[2955] * v[4067] + v[2970] * v[4069] + v[2974] * v[4071] + v[2972] * v[4072] + v[4073]
					+ v[2949] * v[4076] + v[2957] * v[4077] + v[2968] * v[4078] + v[4079] + v[2953] * v[4080] + v[2960] * v[4081]
					+ v[2966] * v[4084] + v[2872] * v[4178] + v[2870] * v[4179] + v[2868] * v[4181] + v[2864] * v[4191] + v[2861] * v[4192]
					+ v[2860] * v[4195] + v[2853] * v[4209] + v[2856] * v[4211] + v[2852] * v[4212] + v[2818] * v[4423] + v[2823] * v[4424]
					+ v[2829] * v[4425] + v[2947] * v[4426] + v[2951] * v[4427] + v[2952] * v[4428] + v[4057] * v[4432] + v[4058] * v[4433]
					+ v[4062] * v[4434] + v[4063] * v[4435] + v[4074] * v[4436] + v[4075] * v[4437] + v[3984] * v[6840] + v[3985] * v[6841]
					+ v[3986] * v[6842] + v[2843] * v[7174] * v[7373]));
		v[6914] = v[4513] + (v[2263] * v[2635] - v[2260] * v[2638] + v[2257] * v[2640] - v[2276] * v[2821] - v[238] * v[4220]
			+ v[4167] * v[691] - v[4170] * v[692] + v[4173] * v[693]) / 2e0;
		v[4439] = v[4233] + v[4237];
		v[4440] = v[4236] + v[4238];
		v[4441] = v[4231] + v[4234];
		v[4442] = v[1024] * (-(v[2327] * v[2766]) - v[2329] * v[2768] - v[2324] * v[2770] - v[2326] * v[2772]
			- v[2321] * v[2774] - v[2323] * v[2776] - v[1394] * v[3975] - v[1399] * v[3979] - v[1404] * v[3983] - v[197] * v[4241]
			- v[206] * v[4243] - v[200] * v[4244] - v[209] * v[4246] - v[203] * v[4247] - v[212] * v[4249] + v[223] * (-
			(v[4001] * v[565]) - v[4000] * v[566] - v[3999] * v[567] - v[4002] * v[574] - v[4003] * v[575] - v[4004] * v[576]
				- v[4012] * v[583] - v[4013] * v[584] - v[4014] * v[585] - v[6843] - v[6844] - v[6845] - v[2795] * v[82] - v[2793] * v[85]
				- v[2792] * v[88]) + v[224] * (-(v[3998] * v[592]) - v[3997] * v[593] - v[3996] * v[594] - v[4005] * v[601]
					- v[4006] * v[602] - v[4007] * v[603] - v[4015] * v[610] - v[4016] * v[611] - v[4017] * v[612] - v[6846] - v[6847]
					- v[6848] - v[2801] * v[91] - v[2799] * v[94] - v[2798] * v[97])) + v[1610] * (v[2327] * v[2767] + v[2329] * v[2769]
						+ v[2324] * v[2771] + v[2326] * v[2773] + v[2321] * v[2775] + v[2323] * v[2777] - v[1393] * v[3975] - v[1398] * v[3979]
						- v[1403] * v[3983] + v[198] * v[4241] + v[207] * v[4243] + v[201] * v[4244] + v[210] * v[4246] + v[204] * v[4247]
						+ v[213] * v[4249] + v[223] * (v[4001] * v[568] + v[4000] * v[569] + v[3999] * v[570] + v[4002] * v[577] + v[4003] * v[578]
							+ v[4004] * v[579] + v[4012] * v[586] + v[4013] * v[587] + v[4014] * v[588] + v[6849] + v[2795] * v[83] + v[2793] * v[86]
							+ v[2792] * v[89]) + v[224] * (v[3998] * v[595] + v[3997] * v[596] + v[3996] * v[597] + v[4005] * v[604] + v[4006] * v[605]
								+ v[4007] * v[606] + v[4015] * v[613] + v[4016] * v[614] + v[4017] * v[615] + v[6850] + v[2801] * v[92] + v[2799] * v[95]
								+ v[2798] * v[98]));
		v[4443] = (-(v[2332] * v[2544]) - v[2368] * v[2547] - v[2349] * v[2610] - 2e0*v[4089] + 2e0*v[4096] - v[4288] * v[523]
			- v[4269] * v[542] - v[4252] * v[560] + v[4251] * v[6632] + v[4250] * v[6633] + v[4270] * v[6634] + v[4268] * v[6635]
			+ v[4287] * v[6636] + v[4286] * v[6637] - v[2348] * v[6851] - v[2366] * v[6852] - v[2331] * v[6853] - v[2350] * v[6854]
			- v[2330] * v[6855] - v[2367] * v[6856]) / 2e0;
		v[4444] = (-(v[2372] * v[2609]) + v[2369] * v[2612] - v[2368] * v[2614] + v[2416] * v[2699] + v[178] * v[4379]
			- v[4288] * v[520] + v[4289] * v[521] - v[4292] * v[522]) / 2e0;
		v[4446] = (v[2335] * v[2544] + v[2369] * v[2547] + v[2353] * v[2610] - 2e0*v[4090] + 2e0*v[4115] + v[4289] * v[523]
			+ v[4273] * v[542] + v[4255] * v[560] - v[4253] * v[6632] - v[4254] * v[6633] - v[4272] * v[6634] - v[4271] * v[6635]
			- v[4291] * v[6636] - v[4290] * v[6637] + v[2351] * v[6851] + v[2370] * v[6852] + v[2333] * v[6853] + v[2352] * v[6854]
			+ v[2334] * v[6855] + v[2371] * v[6856]) / 2e0;
		v[4447] = (-(v[2356] * v[2609]) + v[2353] * v[2612] - v[2349] * v[2614] + v[2406] * v[2699] + v[178] * v[4375]
			- v[4269] * v[520] + v[4273] * v[521] - v[4276] * v[522]) / 2e0;
		v[4448] = (-(v[2338] * v[2544]) - v[2372] * v[2547] - v[2356] * v[2610] - 2e0*v[4100] + 2e0*v[4112] - v[4292] * v[523]
			- v[4276] * v[542] - v[4258] * v[560] + v[4256] * v[6632] + v[4257] * v[6633] + v[4274] * v[6634] + v[4275] * v[6635]
			+ v[4293] * v[6636] + v[4294] * v[6637] - v[2355] * v[6851] - v[2374] * v[6852] - v[2336] * v[6853] - v[2354] * v[6854]
			- v[2337] * v[6855] - v[2373] * v[6856]) / 2e0;
		v[4506] = v[2611] * v[4446] - v[4447] + v[2608] * v[4448] + (v[11436 + i1687] + v[2457] * v[2546] - v[4443] * v[463]
			)*v[6429] + v[2546] * v[7189] * v[7375] - v[6430] * (v[11814 + i1687] + (v[2396] * v[2544]) / 2e0 + (v[2416] * v[2547])
				/ 2e0 + v[2407] * v[2602] + v[2415] * v[2603] + v[2397] * v[2604] + v[2405] * v[2605] + v[2398] * v[2606]
				+ v[2414] * v[2607] + (v[2406] * v[2610]) / 2e0 + v[4251] - v[4254] - v[4270] + v[4275] + v[4291] - v[4294]
				+ 2e0*v[2699] * v[4450] + v[4380] * v[464] - v[4381] * v[465] - v[4383] * v[467] + v[4378] * v[528] + v[4377] * v[534]
				+ v[4376] * v[538] + v[4374] * v[547] + v[4373] * v[551] + v[4372] * v[555] + v[4371] * v[6605] + v[4375] * v[6606]
				+ v[4379] * v[6607] + v[4386] * v[6859] + v[4390] * v[6860] + v[4391] * v[6861] + v[6647] * (v[4091] + v[2936] * v[4094]
					+ v[2933] * v[4095] + v[2927] * v[4097] + v[2942] * v[4099] + v[2946] * v[4101] + v[2944] * v[4102] + v[4103]
					+ v[2921] * v[4106] + v[2929] * v[4107] + v[2940] * v[4108] + v[4109] + v[2925] * v[4110] + v[2932] * v[4111]
					+ v[2938] * v[4114] + v[2759] * v[4312] + v[2757] * v[4313] + v[2755] * v[4315] + v[2754] * v[4322] + v[2751] * v[4323]
					+ v[2750] * v[4326] + v[2746] * v[4333] + v[2749] * v[4335] + v[2745] * v[4336] + v[2696] * v[4451] + v[2701] * v[4452]
					+ v[2707] * v[4453] + v[2919] * v[4454] + v[2923] * v[4455] + v[2924] * v[4456] + v[4087] * v[4460] + v[4088] * v[4461]
					+ v[4092] * v[4462] + v[4093] * v[4463] + v[4104] * v[4464] + v[4105] * v[4465] + v[3987] * v[6862] + v[3988] * v[6863]
					+ v[3989] * v[6864] + v[2721] * v[7190] * v[7377]));
		v[6900] = v[4506] + (v[2338] * v[2609] - v[2335] * v[2612] + v[2332] * v[2614] - v[2396] * v[2699] - v[178] * v[4371]
			+ v[4252] * v[520] - v[4255] * v[521] + v[4258] * v[522]) / 2e0;
		v[4467] = v[4384] + v[4388];
		v[4468] = v[4387] + v[4389];
		v[4469] = v[4382] + v[4385];
		v[4470] = (-(v[2341] * v[2534]) - v[2377] * v[2537] - v[2358] * v[2584] - 2e0*v[4119] + 2e0*v[4126] - v[4297] * v[478]
			- v[4278] * v[497] - v[4261] * v[515] + v[4260] * v[6638] + v[4259] * v[6639] + v[4279] * v[6640] + v[4277] * v[6641]
			+ v[4296] * v[6642] + v[4295] * v[6643] - v[2357] * v[6865] - v[2375] * v[6866] - v[2340] * v[6867] - v[2359] * v[6868]
			- v[2339] * v[6869] - v[2376] * v[6870]) / 2e0;
		v[4471] = (-(v[2381] * v[2583]) + v[2378] * v[2586] - v[2377] * v[2588] + v[2443] * v[2661] + v[159] * v[4400]
			- v[4297] * v[475] + v[4298] * v[476] - v[4301] * v[477]) / 2e0;
		v[4473] = (v[2344] * v[2534] + v[2378] * v[2537] + v[2362] * v[2584] - 2e0*v[4120] + 2e0*v[4145] + v[4298] * v[478]
			+ v[4282] * v[497] + v[4264] * v[515] - v[4262] * v[6638] - v[4263] * v[6639] - v[4281] * v[6640] - v[4280] * v[6641]
			- v[4300] * v[6642] - v[4299] * v[6643] + v[2360] * v[6865] + v[2379] * v[6866] + v[2342] * v[6867] + v[2361] * v[6868]
			+ v[2343] * v[6869] + v[2380] * v[6870]) / 2e0;
		v[4474] = (-(v[2365] * v[2583]) + v[2362] * v[2586] - v[2358] * v[2588] + v[2433] * v[2661] + v[159] * v[4396]
			- v[4278] * v[475] + v[4282] * v[476] - v[4285] * v[477]) / 2e0;
		v[4475] = (-(v[2347] * v[2534]) - v[2381] * v[2537] - v[2365] * v[2584] - 2e0*v[4130] + 2e0*v[4142] - v[4301] * v[478]
			- v[4285] * v[497] - v[4267] * v[515] + v[4265] * v[6638] + v[4266] * v[6639] + v[4283] * v[6640] + v[4284] * v[6641]
			+ v[4302] * v[6642] + v[4303] * v[6643] - v[2364] * v[6865] - v[2383] * v[6866] - v[2345] * v[6867] - v[2363] * v[6868]
			- v[2346] * v[6869] - v[2382] * v[6870]) / 2e0;
		v[4499] = v[2585] * v[4473] - v[4474] + v[2582] * v[4475] + (v[11418 + i1687] + v[2471] * v[2536] - v[4470] * v[457]
			)*v[6427] + v[2536] * v[7199] * v[7379] - v[6428] * (v[11832 + i1687] + (v[2423] * v[2534]) / 2e0 + (v[2443] * v[2537])
				/ 2e0 + v[2434] * v[2576] + v[2442] * v[2577] + v[2424] * v[2578] + v[2432] * v[2579] + v[2425] * v[2580]
				+ v[2441] * v[2581] + (v[2433] * v[2584]) / 2e0 + v[4260] - v[4263] - v[4279] + v[4284] + v[4300] - v[4303]
				+ 2e0*v[2661] * v[4477] + v[4401] * v[458] - v[4402] * v[459] - v[4404] * v[461] + v[4399] * v[483] + v[4398] * v[489]
				+ v[4397] * v[493] + v[4395] * v[502] + v[4394] * v[506] + v[4393] * v[510] + v[4392] * v[6623] + v[4396] * v[6624]
				+ v[4400] * v[6625] + v[4407] * v[6873] + v[4411] * v[6874] + v[4412] * v[6875] + v[6646] * (v[4121] + v[2908] * v[4124]
					+ v[2905] * v[4125] + v[2899] * v[4127] + v[2914] * v[4129] + v[2918] * v[4131] + v[2916] * v[4132] + v[4133]
					+ v[2893] * v[4136] + v[2901] * v[4137] + v[2912] * v[4138] + v[4139] + v[2897] * v[4140] + v[2904] * v[4141]
					+ v[2910] * v[4144] + v[2744] * v[4342] + v[2742] * v[4343] + v[2740] * v[4345] + v[2739] * v[4352] + v[2736] * v[4353]
					+ v[2735] * v[4356] + v[2731] * v[4363] + v[2734] * v[4365] + v[2730] * v[4366] + v[2658] * v[4478] + v[2663] * v[4479]
					+ v[2669] * v[4480] + v[2891] * v[4481] + v[2895] * v[4482] + v[2896] * v[4483] + v[4117] * v[4487] + v[4118] * v[4488]
					+ v[4122] * v[4489] + v[4123] * v[4490] + v[4134] * v[4491] + v[4135] * v[4492] + v[3990] * v[6876] + v[3991] * v[6877]
					+ v[3992] * v[6878] + v[2683] * v[7200] * v[7381]));
		v[6886] = v[4499] + (v[2347] * v[2583] - v[2344] * v[2586] + v[2341] * v[2588] - v[2423] * v[2661] - v[159] * v[4392]
			+ v[4261] * v[475] - v[4264] * v[476] + v[4267] * v[477]) / 2e0;
		v[4494] = v[4405] + v[4409];
		v[4495] = v[4408] + v[4410];
		v[4496] = v[4403] + v[4406];
		v[12157] = 0e0;
		v[12158] = 0e0;
		v[12159] = 0e0;
		v[12160] = 0e0;
		v[12161] = v[2498];
		v[12162] = v[2496];
		v[12163] = 0e0;
		v[12164] = 0e0;
		v[12165] = 0e0;
		v[12166] = 0e0;
		v[12167] = 0e0;
		v[12168] = 0e0;
		v[12169] = 0e0;
		v[12170] = 0e0;
		v[12171] = 0e0;
		v[12172] = 0e0;
		v[12173] = 0e0;
		v[12174] = 0e0;
		v[12103] = 0e0;
		v[12104] = 0e0;
		v[12105] = 0e0;
		v[12106] = v[2498];
		v[12107] = 0e0;
		v[12108] = v[2497];
		v[12109] = 0e0;
		v[12110] = 0e0;
		v[12111] = 0e0;
		v[12112] = 0e0;
		v[12113] = 0e0;
		v[12114] = 0e0;
		v[12115] = 0e0;
		v[12116] = 0e0;
		v[12117] = 0e0;
		v[12118] = 0e0;
		v[12119] = 0e0;
		v[12120] = 0e0;
		v[12067] = 0e0;
		v[12068] = 0e0;
		v[12069] = 0e0;
		v[12070] = v[2496];
		v[12071] = v[2497];
		v[12072] = 0e0;
		v[12073] = 0e0;
		v[12074] = 0e0;
		v[12075] = 0e0;
		v[12076] = 0e0;
		v[12077] = 0e0;
		v[12078] = 0e0;
		v[12079] = 0e0;
		v[12080] = 0e0;
		v[12081] = 0e0;
		v[12082] = 0e0;
		v[12083] = 0e0;
		v[12084] = 0e0;
		v[12049] = 0e0;
		v[12050] = 0e0;
		v[12051] = 0e0;
		v[12052] = 0e0;
		v[12053] = 0e0;
		v[12054] = 0e0;
		v[12055] = 0e0;
		v[12056] = 0e0;
		v[12057] = 0e0;
		v[12058] = 0e0;
		v[12059] = v[2489];
		v[12060] = v[2487];
		v[12061] = 0e0;
		v[12062] = 0e0;
		v[12063] = 0e0;
		v[12064] = 0e0;
		v[12065] = 0e0;
		v[12066] = 0e0;
		v[11995] = 0e0;
		v[11996] = 0e0;
		v[11997] = 0e0;
		v[11998] = 0e0;
		v[11999] = 0e0;
		v[12000] = 0e0;
		v[12001] = 0e0;
		v[12002] = 0e0;
		v[12003] = 0e0;
		v[12004] = v[2489];
		v[12005] = 0e0;
		v[12006] = v[2488];
		v[12007] = 0e0;
		v[12008] = 0e0;
		v[12009] = 0e0;
		v[12010] = 0e0;
		v[12011] = 0e0;
		v[12012] = 0e0;
		v[11959] = 0e0;
		v[11960] = 0e0;
		v[11961] = 0e0;
		v[11962] = 0e0;
		v[11963] = 0e0;
		v[11964] = 0e0;
		v[11965] = 0e0;
		v[11966] = 0e0;
		v[11967] = 0e0;
		v[11968] = v[2487];
		v[11969] = v[2488];
		v[11970] = 0e0;
		v[11971] = 0e0;
		v[11972] = 0e0;
		v[11973] = 0e0;
		v[11974] = 0e0;
		v[11975] = 0e0;
		v[11976] = 0e0;
		v[11941] = 0e0;
		v[11942] = 0e0;
		v[11943] = 0e0;
		v[11944] = 0e0;
		v[11945] = 0e0;
		v[11946] = 0e0;
		v[11947] = 0e0;
		v[11948] = 0e0;
		v[11949] = 0e0;
		v[11950] = 0e0;
		v[11951] = 0e0;
		v[11952] = 0e0;
		v[11953] = 0e0;
		v[11954] = 0e0;
		v[11955] = 0e0;
		v[11956] = 0e0;
		v[11957] = v[2480];
		v[11958] = v[2478];
		v[11887] = 0e0;
		v[11888] = 0e0;
		v[11889] = 0e0;
		v[11890] = 0e0;
		v[11891] = 0e0;
		v[11892] = 0e0;
		v[11893] = 0e0;
		v[11894] = 0e0;
		v[11895] = 0e0;
		v[11896] = 0e0;
		v[11897] = 0e0;
		v[11898] = 0e0;
		v[11899] = 0e0;
		v[11900] = 0e0;
		v[11901] = 0e0;
		v[11902] = v[2480];
		v[11903] = 0e0;
		v[11904] = v[2479];
		v[11851] = 0e0;
		v[11852] = 0e0;
		v[11853] = 0e0;
		v[11854] = 0e0;
		v[11855] = 0e0;
		v[11856] = 0e0;
		v[11857] = 0e0;
		v[11858] = 0e0;
		v[11859] = 0e0;
		v[11860] = 0e0;
		v[11861] = 0e0;
		v[11862] = 0e0;
		v[11863] = 0e0;
		v[11864] = 0e0;
		v[11865] = 0e0;
		v[11866] = v[2478];
		v[11867] = v[2479];
		v[11868] = 0e0;
		v[12175] = v[3927] + v[7] * (v[1780] * v[3046] + v[1016] * v[3575] + v[33] * (-(v[1614] * v[3619]) - v[1615] * v[3620]
			- v[1616] * v[3621]) + v[223] * v[6822] - v[3610] * v[869] - v[3609] * v[887] + v[3612] * v[905] + v[3611] * v[923]);
		v[12176] = v[3937] + v[7] * (v[1779] * v[3046] + v[1016] * v[3574] + v[33] * (-(v[1618] * v[3619]) - v[1619] * v[3620]
			- v[1620] * v[3621]) + v[223] * v[6821] - v[3610] * v[870] - v[3609] * v[888] + v[3612] * v[906] + v[3611] * v[924]);
		v[12177] = v[3947] + v[7] * (v[33] * (-(v[1622] * v[3619]) - v[1623] * v[3620] - v[1624] * v[3621]) + v[223] * (v[6682]
			+ v[6820]) - v[3610] * v[871] - v[3609] * v[889] + v[3612] * v[907] + v[3611] * v[925]);
		v[12178] = (v[12156 + i1687] - v[2459] * v[2661] - v[159] * v[4402]) / 2e0 - v[4408] + v[4410] + v[4496] * v[458]
			+ v[4494] * v[461] + 2e0*(v[2575] * v[4498] + v[4470] * v[6428] + v[4411] * v[6434] + v[2536] * (-(v[2490] * v[6427])
				+ v[2468] * v[6885])) + v[457] * v[6886] + (v[1983] * v[2894] + v[1981] * v[2898] + v[2683] * (v[4363]
					+ v[1981] * v[4488] + v[1983] * v[4490]) + v[3992] * v[4497] + v[3991] * v[6457] + v[3990] * v[6460] + v[2600] * v[6879]
					+ v[2660] * v[6880] + v[2667] * v[6881] + v[166] * v[6882] + v[171] * v[6883] + v[162] * v[6884])*v[7];
		v[12179] = (v[12102 + i1687] + v[2458] * v[2661] + v[159] * v[4401]) / 2e0 + v[4405] - v[4409] + v[4496] * v[459]
			+ v[4495] * v[461] + 2e0*(v[2574] * v[4501] - v[4473] * v[6428] + v[4412] * v[6434] + v[2536] * (v[2492] * v[6427]
				+ v[2469] * v[6885])) + v[460] * (-v[4471] + v[4474] + v[6886]) + v[7] * (v[1985] * v[2894] + v[1981] * v[2900]
					+ v[2503] * v[3990] + v[2504] * v[3992] + v[2683] * (v[4353] + v[1981] * v[4487] + v[1985] * v[4491]) + v[3991] * v[4500]
					+ v[2658] * v[6587] + v[2676] * v[6589] + v[2598] * v[6887] + v[163] * v[6888] + v[168] * v[6889] + v[2908] * v[7213]);
		v[12180] = -v[4403] + (v[12066 + i1687] - v[2461] * v[2661] - v[159] * v[4404]) / 2e0 + v[4406] + v[4495] * v[458]
			+ v[4494] * v[459] + (-v[4471] + v[4499])*v[462] + 2e0*(v[2573] * v[4503] + v[4475] * v[6428] + v[4407] * v[6434]
				+ v[2536] * (-(v[2494] * v[6427]) + v[2464] * v[6885])) + v[7] * (v[1985] * v[2898] + v[1983] * v[2900]
					+ v[2505] * v[3991] + v[2506] * v[3992] + v[2683] * (v[4342] + v[1983] * v[4489] + v[1985] * v[4492]) + v[3990] * v[4502]
					+ v[2663] * v[6588] + v[2596] * v[6890] + v[2669] * v[6891] + v[174] * v[6892] + v[2912] * v[7217] + v[2916] * v[7218]);
		v[12181] = v[3926] + v[7] * (v[1780] * v[3045] + v[1017] * v[3575] + v[33] * (-(v[1638] * v[3619]) - v[1639] * v[3620]
			- v[1640] * v[3621]) + v[224] * v[6822] - v[3610] * v[875] - v[3609] * v[893] + v[3612] * v[911] + v[3611] * v[929]);
		v[12182] = v[3936] + v[7] * (v[1779] * v[3045] + v[1017] * v[3574] + v[33] * (-(v[1642] * v[3619]) - v[1643] * v[3620]
			- v[1644] * v[3621]) + v[224] * v[6821] - v[3610] * v[876] - v[3609] * v[894] + v[3612] * v[912] + v[3611] * v[930]);
		v[12183] = v[3946] + v[7] * (v[33] * (-(v[1646] * v[3619]) - v[1647] * v[3620] - v[1648] * v[3621]) + v[224] * (v[6682]
			+ v[6820]) - v[3610] * v[877] - v[3609] * v[895] + v[3612] * v[913] + v[3611] * v[931]);
		v[12184] = (v[12048 + i1687] - v[2445] * v[2699] - v[178] * v[4381]) / 2e0 - v[4387] + v[4389] + v[4469] * v[464]
			+ v[4467] * v[467] + 2e0*(v[2572] * v[4505] + v[4443] * v[6430] + v[4390] * v[6441] + v[2546] * (-(v[2481] * v[6429])
				+ v[2454] * v[6899])) + v[463] * v[6900] + (v[1977] * v[2922] + v[1975] * v[2926] + v[2721] * (v[4333]
					+ v[1975] * v[4461] + v[1977] * v[4463]) + v[3989] * v[4504] + v[3988] * v[6465] + v[3987] * v[6468] + v[2626] * v[6893]
					+ v[2698] * v[6894] + v[2705] * v[6895] + v[185] * v[6896] + v[190] * v[6897] + v[181] * v[6898])*v[7];
		v[12185] = (v[11994 + i1687] + v[2444] * v[2699] + v[178] * v[4380]) / 2e0 + v[4384] - v[4388] + v[4469] * v[465]
			+ v[4468] * v[467] + 2e0*(v[2571] * v[4508] - v[4446] * v[6430] + v[4391] * v[6441] + v[2546] * (v[2483] * v[6429]
				+ v[2455] * v[6899])) + v[466] * (-v[4444] + v[4447] + v[6900]) + v[7] * (v[1979] * v[2922] + v[1975] * v[2928]
					+ v[2511] * v[3987] + v[2512] * v[3989] + v[2721] * (v[4323] + v[1975] * v[4460] + v[1979] * v[4464]) + v[3988] * v[4507]
					+ v[2696] * v[6584] + v[2714] * v[6586] + v[2624] * v[6901] + v[182] * v[6902] + v[187] * v[6903] + v[2936] * v[7225]);
		v[12186] = -v[4382] + (v[11958 + i1687] - v[2447] * v[2699] - v[178] * v[4383]) / 2e0 + v[4385] + v[4468] * v[464]
			+ v[4467] * v[465] + (-v[4444] + v[4506])*v[468] + 2e0*(v[2570] * v[4510] + v[4448] * v[6430] + v[4386] * v[6441]
				+ v[2546] * (-(v[2485] * v[6429]) + v[2450] * v[6899])) + v[7] * (v[1979] * v[2926] + v[1977] * v[2928]
					+ v[2513] * v[3988] + v[2514] * v[3989] + v[2721] * (v[4312] + v[1977] * v[4462] + v[1979] * v[4465]) + v[3987] * v[4509]
					+ v[2701] * v[6585] + v[2622] * v[6904] + v[2707] * v[6905] + v[193] * v[6906] + v[2940] * v[7229] + v[2944] * v[7230]);
		v[12187] = v[3922] + v[7] * (v[3436] + v[33] * (-(v[1662] * v[3619]) - v[1663] * v[3620] - v[1664] * v[3621]) + v[6752]
			- v[389] * (v[3448] + v[6736] - v[6800]) + v[6813] - v[3610] * v[881] - v[3609] * v[899] + v[3612] * v[917]
			+ v[3611] * v[935]);
		v[12188] = v[3932] + v[7] * (v[3243] + v[33] * (-(v[1666] * v[3619]) - v[1667] * v[3620] - v[1668] * v[3621]) + v[6710]
			- v[6799] + v[390] * (-v[3461] + v[6737] + v[6800]) + v[6815] - v[3610] * v[882] - v[3609] * v[900] + v[3612] * v[918]
			+ v[3611] * v[936]);
		v[12189] = v[3942] + v[7] * (v[33] * (-(v[1670] * v[3619]) - v[1671] * v[3620] - v[1672] * v[3621]) + v[6707] - v[6766]
			- v[6817] + v[6819] - v[3610] * v[883] - v[3609] * v[901] + v[3612] * v[919] + v[3611] * v[937]);
		v[12190] = (v[11940 + i1687] - v[2308] * v[2821] - v[238] * v[4230]) / 2e0 - v[4236] + v[4238] + v[4441] * v[470]
			+ v[4439] * v[473] + 2e0*(v[2569] * v[4512] + v[4415] * v[6432] + v[4239] * v[6448] + v[2556] * (-(v[2472] * v[6431])
				+ v[2317] * v[6913])) + v[469] * v[6914] + (v[1971] * v[2950] + v[1969] * v[2954] + v[2843] * (v[4209]
					+ v[1969] * v[4433] + v[1971] * v[4435]) + v[3986] * v[4511] + v[3985] * v[6473] + v[3984] * v[6476] + v[2652] * v[6907]
					+ v[2820] * v[6908] + v[2827] * v[6909] + v[245] * v[6910] + v[250] * v[6911] + v[241] * v[6912])*v[7];
		v[12191] = (v[11886 + i1687] + v[2307] * v[2821] + v[238] * v[4229]) / 2e0 + v[4233] - v[4237] + v[4441] * v[471]
			+ v[4440] * v[473] + 2e0*(v[2568] * v[4515] - v[4418] * v[6432] + v[4240] * v[6448] + v[2556] * (v[2474] * v[6431]
				+ v[2318] * v[6913])) + v[472] * (-v[4416] + v[4419] + v[6914]) + v[7] * (v[1973] * v[2950] + v[1969] * v[2956]
					+ v[2519] * v[3984] + v[2520] * v[3986] + v[2843] * (v[4192] + v[1969] * v[4432] + v[1973] * v[4436]) + v[3985] * v[4514]
					+ v[2818] * v[6557] + v[2836] * v[6559] + v[2650] * v[6915] + v[242] * v[6916] + v[247] * v[6917] + v[2964] * v[7237]);
		v[12192] = -v[4231] + (v[11850 + i1687] - v[2310] * v[2821] - v[238] * v[4232]) / 2e0 + v[4234] + v[4440] * v[470]
			+ v[4439] * v[471] + (-v[4416] + v[4513])*v[474] + 2e0*(v[2567] * v[4517] + v[4420] * v[6432] + v[4235] * v[6448]
				+ v[2556] * (-(v[2476] * v[6431]) + v[2313] * v[6913])) + v[7] * (v[1973] * v[2954] + v[1971] * v[2956]
					+ v[2521] * v[3985] + v[2522] * v[3986] + v[2843] * (v[4178] + v[1971] * v[4434] + v[1973] * v[4437]) + v[3984] * v[4516]
					+ v[2823] * v[6558] + v[2648] * v[6918] + v[2829] * v[6919] + v[253] * v[6920] + v[2968] * v[7241] + v[2972] * v[7242]);
		Rc[i1687 - 1] += v[10424 + i1687];
		for (i2526 = 1; i2526 <= 18; i2526++) {
			Kc[i1687 - 1][i2526 - 1] += v[12174 + i2526] + v[4413] * v[9120 + i2526] + v[4414] * v[9138 + i2526] + v[4304] * v[9156
				+ i2526] + v[4442] * v[9174 + i2526];
		};/* end for */
	};/* end for */
	v[4523] = v[1044];
	v[4696] = v[4523];
	v[4524] = v[1043];
	v[4698] = v[4524];
	v[4525] = v[1042];
	v[4700] = v[4525];
	v[4529] = 0e0;
	v[4530] = 0e0;
	v[4531] = 0e0;
	v[4532] = 0e0;
	v[4533] = 0e0;
	v[4534] = 0e0;
	v[4535] = 0e0;
	b4536 = b6;
	if (b4536) {
		b4537 = b1040;
		v[4531] = 0e0;
		v[4530] = 0e0;
		v[4529] = 0e0;
		v[4532] = 0e0;
	}
	else {
		v[4544] = 0e0;
		b4545 = b1049;
		if (b4545) {
			b4546 = b1064;
			b4550 = b1051;
			if (b4550) {
				v[4535] = 0e0;
				v[4534] = 0e0;
				v[4533] = 0e0;
				v[4544] = 0e0;
			}
			else {
			};
		}
		else {
		};
		b4551 = b1049;
		if (b4551) {
			b4552 = b1051;
			if (b4552) {
				v[4531] = 0e0;
				v[4530] = 0e0;
				v[4529] = 0e0;
				v[4532] = -v[4544];
			}
			else {
			};
		}
		else {
		};
	};
	v[4701] = v[4532];
	v[6930] = v[1014] * v[4533];
	v[6922] = v[389] * v[4533];
	v[6929] = v[1018] * v[4534];
	v[6921] = v[390] * v[4534];
	v[4557] = v[372] * v[4535];
	v[5679] = -(v[391] * v[4557]);
	v[4558] = v[362] * v[4535];
	v[5676] = -(v[391] * v[4558]);
	v[4559] = v[352] * v[4535];
	v[5673] = -(v[391] * v[4559]);
	v[4560] = v[346] * v[4535];
	v[5688] = v[391] * v[4560];
	v[4561] = v[336] * v[4535];
	v[5685] = v[391] * v[4561];
	v[4562] = v[326] * v[4535];
	v[5682] = v[391] * v[4562];
	v[4563] = v[320] * v[4535];
	v[5695] = v[391] * v[4563];
	v[4564] = v[310] * v[4535];
	v[5692] = v[391] * v[4564];
	v[4565] = v[300] * v[4535];
	v[7025] = v[4565] * v[673] + v[4564] * v[674] + v[4563] * v[675] + v[4562] * v[676] + v[4561] * v[677] + v[4560] * v[678]
		- v[4559] * v[781] - v[4558] * v[782] - v[4557] * v[783];
	v[6983] = v[4565] * v[679] + v[4564] * v[680] + v[4563] * v[681] + v[4562] * v[682] + v[4561] * v[683] + v[4560] * v[684]
		- v[4559] * v[784] - v[4558] * v[785] - v[4557] * v[786];
	v[5641] = v[391] * v[4565];
	v[4529] = v[4529] + v[3825] * v[4535];
	v[4530] = v[4530] + v[3826] * v[4535];
	v[5515] = v[4535] * v[4566] + v[4565] * v[685] + v[4564] * v[686] + v[4563] * v[687] + v[4562] * v[688] + v[4561] * v[689]
		+ v[4560] * v[690] - v[4559] * v[787] - v[4558] * v[788] - v[4557] * v[789];
	v[4568] = v[389] * v[4535];
	v[4570] = v[390] * v[4535];
	v[4531] = v[4531] + v[4535] * v[4572];
	v[4582] = -(v[1022] * v[4535]);
	v[4583] = -(v[391] * v[4570]);
	v[6951] = v[4583] - v[6929];
	v[4584] = -(v[391] * v[4568]);
	v[6952] = v[4584] - v[6930];
	v[4591] = v[372] * v[4534];
	v[5680] = -(v[390] * v[4591]);
	v[7033] = -v[5679] - v[5680];
	v[4592] = v[362] * v[4534];
	v[5677] = -(v[390] * v[4592]);
	v[7032] = -v[5676] - v[5677];
	v[4593] = v[352] * v[4534];
	v[5674] = -(v[390] * v[4593]);
	v[7031] = -v[5673] - v[5674];
	v[4594] = v[346] * v[4534];
	v[5689] = v[390] * v[4594];
	v[7028] = v[5688] + v[5689];
	v[4595] = v[336] * v[4534];
	v[5686] = v[390] * v[4595];
	v[7027] = v[5685] + v[5686];
	v[4596] = v[326] * v[4534];
	v[5683] = v[390] * v[4596];
	v[7026] = v[5682] + v[5683];
	v[4597] = v[320] * v[4534];
	v[5696] = v[390] * v[4597];
	v[7030] = v[5695] + v[5696];
	v[4598] = v[310] * v[4534];
	v[5693] = v[390] * v[4598];
	v[7029] = v[5692] + v[5693];
	v[4599] = v[300] * v[4534];
	v[6972] = v[4599] * v[685] + v[4598] * v[686] + v[4597] * v[687] + v[4596] * v[688] + v[4595] * v[689] + v[4594] * v[690]
		- v[4593] * v[787] - v[4592] * v[788] - v[4591] * v[789];
	v[5642] = v[390] * v[4599];
	v[7022] = v[5641] + v[5642];
	v[4529] = v[4529] - v[292] * v[6921];
	v[5572] = v[4534] * v[4600] + v[4599] * v[679] + v[4598] * v[680] + v[4597] * v[681] + v[4596] * v[682] + v[4595] * v[683]
		+ v[4594] * v[684] - v[4593] * v[784] - v[4592] * v[785] - v[4591] * v[786];
	v[4530] = v[4530] + v[4534] * v[4602];
	v[4531] = v[4531] - v[294] * v[6921];
	v[4623] = v[372] * v[4533];
	v[5595] = -(v[389] * v[4623]);
	v[7452] = 2e0*v[5595] + v[5679] + v[5680];
	v[7426] = v[5595] + v[5679] - v[4591] * v[6683];
	v[7408] = v[5595] + v[5680] - v[4557] * v[6680];
	v[7042] = -v[5595] - v[5680];
	v[7038] = -v[5595] - v[5679];
	v[4624] = v[362] * v[4533];
	v[5593] = -(v[389] * v[4624]);
	v[7451] = 2e0*v[5593] + v[5676] + v[5677];
	v[7425] = v[5593] + v[5676] - v[4592] * v[6683];
	v[7407] = v[5593] + v[5677] - v[4558] * v[6680];
	v[7041] = -v[5593] - v[5677];
	v[7037] = -v[5593] - v[5676];
	v[4625] = v[352] * v[4533];
	v[5591] = -(v[389] * v[4625]);
	v[7450] = 2e0*v[5591] + v[5673] + v[5674];
	v[7424] = v[5591] + v[5673] - v[4593] * v[6683];
	v[7406] = v[5591] + v[5674] - v[4559] * v[6680];
	v[7040] = -v[5591] - v[5674];
	v[7036] = -v[5591] - v[5673];
	v[4626] = v[346] * v[4533];
	v[5601] = v[389] * v[4626];
	v[7455] = 2e0*v[5601] + v[7028];
	v[6986] = v[5601] + v[5688];
	v[7429] = 2e0*v[5689] + v[6986];
	v[6976] = v[5601] + v[5689];
	v[7411] = 2e0*v[5688] + v[6976];
	v[4627] = v[336] * v[4533];
	v[5599] = v[389] * v[4627];
	v[7454] = 2e0*v[5599] + v[7027];
	v[6985] = v[5599] + v[5685];
	v[7428] = 2e0*v[5686] + v[6985];
	v[6975] = v[5599] + v[5686];
	v[7410] = 2e0*v[5685] + v[6975];
	v[4628] = v[326] * v[4533];
	v[5597] = v[389] * v[4628];
	v[7453] = 2e0*v[5597] + v[7026];
	v[6984] = v[5597] + v[5682];
	v[7427] = 2e0*v[5683] + v[6984];
	v[6974] = v[5597] + v[5683];
	v[7409] = 2e0*v[5682] + v[6974];
	v[4629] = v[320] * v[4533];
	v[5607] = v[389] * v[4629];
	v[7458] = 2e0*v[5607] + v[7030];
	v[6989] = v[5607] + v[5695];
	v[7432] = 2e0*v[5696] + v[6989];
	v[6979] = v[5607] + v[5696];
	v[7414] = 2e0*v[5695] + v[6979];
	v[4630] = v[310] * v[4533];
	v[5605] = v[389] * v[4630];
	v[7457] = 2e0*v[5605] + v[7029];
	v[6988] = v[5605] + v[5692];
	v[7431] = 2e0*v[5693] + v[6988];
	v[6978] = v[5605] + v[5693];
	v[7413] = 2e0*v[5692] + v[6978];
	v[4631] = v[300] * v[4533];
	v[6970] = v[4631] * v[685] + v[4630] * v[686] + v[4629] * v[687] + v[4628] * v[688] + v[4627] * v[689] + v[4626] * v[690]
		- v[4625] * v[787] - v[4624] * v[788] - v[4623] * v[789];
	v[5603] = v[389] * v[4631];
	v[7456] = 2e0*v[5603] + v[7022];
	v[6987] = v[5603] + v[5641];
	v[7430] = 2e0*v[5642] + v[6987];
	v[6977] = v[5603] + v[5642];
	v[7412] = 2e0*v[5641] + v[6977];
	v[5643] = v[4533] * v[4632] + v[4631] * v[673] + v[4630] * v[674] + v[4629] * v[675] + v[4628] * v[676] + v[4627] * v[677]
		+ v[4626] * v[678] - v[4625] * v[781] - v[4624] * v[782] - v[4623] * v[783];
	v[4634] = v[287] * v[4533] + v[286] * v[4534];
	v[4635] = v[290] * v[4533] + v[289] * v[4534];
	v[7384] = v[4634] - v[4635];
	v[6931] = v[223] * v[4634] + v[224] * v[4635];
	v[7024] = v[4599] * v[673] + v[4598] * v[674] + v[4597] * v[675] + v[4596] * v[676] + v[4595] * v[677] + v[4594] * v[678]
		+ v[6931] - v[4593] * v[781] - v[4592] * v[782] - v[4591] * v[783];
	v[6982] = v[4631] * v[679] + v[4630] * v[680] + v[4629] * v[681] + v[4628] * v[682] + v[4627] * v[683] + v[4626] * v[684]
		+ v[6931] - v[4625] * v[784] - v[4624] * v[785] - v[4623] * v[786];
	v[4529] = v[4529] + v[4533] * v[4636];
	v[4530] = v[4530] - v[293] * v[6922];
	v[4637] = v[6921] + v[6922];
	v[12471] = 0e0;
	v[12472] = 0e0;
	v[12473] = v[223] * v[4637];
	v[12474] = 0e0;
	v[12475] = 0e0;
	v[12476] = 0e0;
	v[12477] = 0e0;
	v[12478] = 0e0;
	v[12479] = v[224] * v[4637];
	v[12480] = 0e0;
	v[12481] = 0e0;
	v[12482] = 0e0;
	v[12483] = 0e0;
	v[12484] = 0e0;
	v[12485] = 0e0;
	v[12486] = 0e0;
	v[12487] = 0e0;
	v[12488] = 0e0;
	v[6953] = -(v[391] * v[4637]);
	v[4531] = v[4531] - v[294] * v[6922];
	v[4639] = v[1968] * v[4533] + v[1889] * v[4534] + v[1816] * v[4535];
	v[6923] = -(v[2319] * v[4639]);
	v[6107] = v[6475] * v[6923];
	v[6097] = v[6477] * v[6923];
	v[6090] = v[363] * v[6923];
	v[4785] = v[2521] * v[4639];
	v[4781] = v[2522] * v[4639];
	v[4771] = v[4516] * v[4639];
	v[4640] = v[1970] * v[4533] + v[1891] * v[4534] + v[1818] * v[4535];
	v[6108] = v[4640] * v[6566];
	v[6098] = v[4430] * v[4640];
	v[6091] = v[4640] * v[6924];
	v[4789] = v[2519] * v[4640];
	v[7235] = v[4771] + v[4789];
	v[4775] = v[4514] * v[4640];
	v[7239] = v[4775] + v[4785];
	v[4641] = v[1972] * v[4533] + v[1893] * v[4534] + v[1820] * v[4535];
	v[6105] = v[4431] * v[4641];
	v[7134] = v[6105] + v[6108];
	v[7481] = v[6107] + v[7134];
	v[6101] = v[4641] * v[6567];
	v[7133] = v[6098] + v[6101];
	v[7480] = v[6097] + v[7133];
	v[6093] = v[4641] * v[6569];
	v[7131] = v[6090] + v[6093];
	v[7479] = v[6091] + v[7131];
	v[4779] = v[4511] * v[4641];
	v[7240] = v[4779] + v[4781];
	v[6936] = v[2520] * v[4640] + v[4779];
	v[7231] = v[4781] + v[6936];
	v[6935] = v[4775] + v[4641] * v[6473];
	v[7234] = v[4785] + v[6935];
	v[6934] = v[4771] + v[4641] * v[6476];
	v[7238] = v[4789] + v[6934];
	v[4642] = v[1974] * v[4533] + v[1895] * v[4534] + v[1822] * v[4535];
	v[6925] = -(v[2456] * v[4642]);
	v[6234] = v[6467] * v[6925];
	v[6224] = v[6469] * v[6925];
	v[6217] = v[337] * v[6925];
	v[4818] = v[2513] * v[4642];
	v[4814] = v[2514] * v[4642];
	v[4804] = v[4509] * v[4642];
	v[4643] = v[1976] * v[4533] + v[1897] * v[4534] + v[1824] * v[4535];
	v[6235] = v[4643] * v[6596];
	v[6225] = v[4458] * v[4643];
	v[7159] = v[6224] + v[6225];
	v[6218] = v[4643] * v[6926];
	v[7158] = v[6217] + v[6218];
	v[4822] = v[2511] * v[4643];
	v[7223] = v[4804] + v[4822];
	v[4808] = v[4507] * v[4643];
	v[7227] = v[4808] + v[4818];
	v[4644] = v[1978] * v[4533] + v[1899] * v[4534] + v[1826] * v[4535];
	v[6232] = v[4459] * v[4644];
	v[7160] = v[6232] + v[6235];
	v[7489] = v[6234] + v[7160];
	v[6228] = v[4644] * v[6597];
	v[7487] = v[6225] + v[6228];
	v[7486] = v[6228] + v[7159];
	v[6220] = v[4644] * v[6599];
	v[7484] = v[6217] + v[6220];
	v[7482] = v[6220] + v[7158];
	v[4812] = v[4504] * v[4644];
	v[7228] = v[4812] + v[4814];
	v[6942] = v[2512] * v[4643] + v[4812];
	v[7219] = v[4814] + v[6942];
	v[6941] = v[4808] + v[4644] * v[6465];
	v[7222] = v[4818] + v[6941];
	v[6940] = v[4804] + v[4644] * v[6468];
	v[7226] = v[4822] + v[6940];
	v[4645] = v[1980] * v[4533] + v[1901] * v[4534] + v[1828] * v[4535];
	v[6927] = -(v[2470] * v[4645]);
	v[6294] = v[6459] * v[6927];
	v[6284] = v[6461] * v[6927];
	v[6277] = v[311] * v[6927];
	v[4851] = v[2505] * v[4645];
	v[4847] = v[2506] * v[4645];
	v[4837] = v[4502] * v[4645];
	v[4646] = v[1982] * v[4533] + v[1903] * v[4534] + v[1830] * v[4535];
	v[6295] = v[4646] * v[6614];
	v[6285] = v[4485] * v[4646];
	v[7162] = v[6284] + v[6285];
	v[6278] = v[4646] * v[6928];
	v[7161] = v[6277] + v[6278];
	v[4855] = v[2503] * v[4646];
	v[7211] = v[4837] + v[4855];
	v[4841] = v[4500] * v[4646];
	v[7215] = v[4841] + v[4851];
	v[4647] = v[1984] * v[4533] + v[1905] * v[4534] + v[1832] * v[4535];
	v[6292] = v[4486] * v[4647];
	v[7163] = v[6292] + v[6295];
	v[7497] = v[6294] + v[7163];
	v[6288] = v[4647] * v[6615];
	v[7495] = v[6285] + v[6288];
	v[7494] = v[6288] + v[7162];
	v[6280] = v[4647] * v[6617];
	v[7492] = v[6277] + v[6280];
	v[7490] = v[6280] + v[7161];
	v[4845] = v[4497] * v[4647];
	v[7216] = v[4845] + v[4847];
	v[6948] = v[2504] * v[4646] + v[4845];
	v[7207] = v[4847] + v[6948];
	v[6947] = v[4841] + v[4647] * v[6457];
	v[7210] = v[4851] + v[6947];
	v[6946] = v[4837] + v[4647] * v[6460];
	v[7214] = v[4855] + v[6946];
	v[4648] = v[4582] + v[6953];
	v[12903] = -v[4584];
	v[12904] = -v[4583];
	v[12905] = -v[4648];
	v[12906] = 0e0;
	v[12907] = 0e0;
	v[12908] = 0e0;
	v[12909] = 0e0;
	v[12910] = 0e0;
	v[12911] = 0e0;
	v[12912] = 0e0;
	v[12913] = 0e0;
	v[12914] = 0e0;
	v[12915] = 0e0;
	v[12916] = 0e0;
	v[12917] = 0e0;
	v[12918] = 0e0;
	v[12919] = 0e0;
	v[12920] = 0e0;
	v[12885] = 0e0;
	v[12886] = 0e0;
	v[12887] = 0e0;
	v[12888] = 0e0;
	v[12889] = 0e0;
	v[12890] = 0e0;
	v[12891] = -v[4584];
	v[12892] = -v[4583];
	v[12893] = -v[4648];
	v[12894] = 0e0;
	v[12895] = 0e0;
	v[12896] = 0e0;
	v[12897] = 0e0;
	v[12898] = 0e0;
	v[12899] = 0e0;
	v[12900] = 0e0;
	v[12901] = 0e0;
	v[12902] = 0e0;
	v[12215] = v[1016] * v[4534] + v[223] * (-v[4584] + v[6930]);
	v[12216] = v[1016] * v[4533] + v[223] * (-v[4583] + v[6929]);
	v[12217] = v[223] * (-v[4582] - v[6953]);
	v[12218] = v[4845] + v[4646] * v[6457] + v[4645] * v[6460];
	v[12219] = v[2503] * v[4645] + v[2504] * v[4647] + v[4841];
	v[12220] = v[2505] * v[4646] + v[2506] * v[4647] + v[4837];
	v[12221] = v[1017] * v[4534] + v[224] * (-v[4584] + v[6930]);
	v[12222] = v[1017] * v[4533] + v[224] * (-v[4583] + v[6929]);
	v[12223] = v[224] * (-v[4582] - v[6953]);
	v[12224] = v[4812] + v[4643] * v[6465] + v[4642] * v[6468];
	v[12225] = v[2511] * v[4642] + v[2512] * v[4644] + v[4808];
	v[12226] = v[2513] * v[4643] + v[2514] * v[4644] + v[4804];
	v[12227] = -(v[389] * v[6921]) + v[6952];
	v[12228] = -(v[390] * v[6922]) + v[6951];
	v[12229] = v[4648];
	v[12230] = v[4779] + v[4640] * v[6473] + v[4639] * v[6476];
	v[12231] = v[2519] * v[4639] + v[2520] * v[4641] + v[4775];
	v[12232] = v[2521] * v[4640] + v[2522] * v[4641] + v[4771];
	v[4657] = -(v[1014] * v[4623]) - v[389] * v[7033];
	v[4658] = -(v[1018] * v[4591]) - v[390] * v[7038];
	v[4659] = -(v[1022] * v[4557]) - v[391] * v[7042];
	v[4660] = -(v[1014] * v[4624]) - v[389] * v[7032];
	v[4661] = -(v[1018] * v[4592]) - v[390] * v[7037];
	v[4662] = -(v[1022] * v[4558]) - v[391] * v[7041];
	v[4663] = -(v[1014] * v[4625]) - v[389] * v[7031];
	v[4664] = -(v[1018] * v[4593]) - v[390] * v[7036];
	v[4665] = -(v[1022] * v[4559]) - v[391] * v[7040];
	v[4666] = v[1022] * v[4560] + v[391] * v[6976];
	v[7137] = v[224] * v[4666];
	v[6993] = -(v[4666] * v[612]);
	v[4667] = v[1018] * v[4594] + v[390] * v[6986];
	v[7143] = v[224] * v[4667];
	v[6994] = -(v[4667] * v[603]);
	v[4668] = v[1014] * v[4626] + v[389] * v[7028];
	v[7149] = v[224] * v[4668];
	v[6995] = -(v[4668] * v[594]);
	v[4669] = v[1022] * v[4561] + v[391] * v[6975];
	v[7136] = v[224] * v[4669];
	v[6996] = -(v[4669] * v[611]);
	v[4670] = v[1018] * v[4595] + v[390] * v[6985];
	v[7142] = v[224] * v[4670];
	v[6997] = -(v[4670] * v[602]);
	v[4671] = v[1014] * v[4627] + v[389] * v[7027];
	v[7148] = v[224] * v[4671];
	v[6998] = -(v[4671] * v[593]);
	v[4672] = v[1022] * v[4562] + v[391] * v[6974];
	v[7135] = v[224] * v[4672];
	v[6999] = -(v[4672] * v[610]);
	v[4673] = v[1018] * v[4596] + v[390] * v[6984];
	v[7141] = v[224] * v[4673];
	v[7000] = -(v[4673] * v[601]);
	v[4674] = v[1014] * v[4628] + v[389] * v[7026];
	v[7147] = v[224] * v[4674];
	v[7002] = v[4674] * v[595] + v[4671] * v[596] + v[4668] * v[597] + v[4673] * v[604] + v[4670] * v[605] + v[4667] * v[606]
		+ v[4672] * v[613] + v[4669] * v[614] + v[4666] * v[615];
	v[7001] = -(v[4674] * v[592]);
	v[7436] = v[1024] * v[6993] - v[1024] * (-v[6994] - v[6995] - v[6996] - v[6997] - v[6998] - v[6999] - v[7000] - v[7001])
		+ v[1610] * v[7002];
	v[4675] = v[1022] * v[4563] + v[391] * v[6979];
	v[7140] = v[223] * v[4675];
	v[7006] = -(v[4675] * v[585]);
	v[4676] = v[1018] * v[4597] + v[390] * v[6989];
	v[7146] = v[223] * v[4676];
	v[7007] = -(v[4676] * v[576]);
	v[4677] = v[1014] * v[4629] + v[389] * v[7030];
	v[7152] = v[223] * v[4677];
	v[7008] = -(v[4677] * v[567]);
	v[4678] = v[1022] * v[4564] + v[391] * v[6978];
	v[7139] = v[223] * v[4678];
	v[7009] = -(v[4678] * v[584]);
	v[4679] = v[1018] * v[4598] + v[390] * v[6988];
	v[7145] = v[223] * v[4679];
	v[7010] = -(v[4679] * v[575]);
	v[4680] = v[1014] * v[4630] + v[389] * v[7029];
	v[7151] = v[223] * v[4680];
	v[7011] = -(v[4680] * v[566]);
	v[4529] = v[4529] + v[3283] * v[4623] + v[3285] * v[4624] + v[3287] * v[4625] + v[3289] * v[4626] + v[3291] * v[4627]
		+ v[3293] * v[4628] + v[3295] * v[4629] + v[3297] * v[4630] + v[3299] * v[4631] + v[5643] * v[6686] + v[390] * v[7024]
		+ v[391] * v[7025];
	v[4693] = v[4529];
	v[4681] = v[1022] * v[4565] + v[391] * v[6977];
	v[7138] = v[223] * v[4681];
	v[7012] = -(v[4681] * v[583]);
	v[4682] = v[1018] * v[4599] + v[390] * v[6987];
	v[7144] = v[223] * v[4682];
	v[7013] = -(v[4682] * v[574]);
	v[4683] = v[1014] * v[4631] + v[389] * v[7022];
	v[7150] = v[223] * v[4683];
	v[7015] = v[4683] * v[568] + v[4680] * v[569] + v[4677] * v[570] + v[4682] * v[577] + v[4679] * v[578] + v[4676] * v[579]
		+ v[4681] * v[586] + v[4678] * v[587] + v[4675] * v[588];
	v[7014] = -(v[4683] * v[565]);
	v[7440] = v[1024] * v[7006] - v[1024] * (-v[7007] - v[7008] - v[7009] - v[7010] - v[7011] - v[7012] - v[7013] - v[7014])
		+ v[1610] * v[7015];
	v[4530] = v[4530] + v[4591] * v[6482] + v[4592] * v[6484] + v[4593] * v[6486] + v[4594] * v[6488] + v[4595] * v[6490]
		+ v[4596] * v[6492] + v[4597] * v[6494] + v[4598] * v[6496] + v[4599] * v[6498] + v[5572] * v[6683] + v[389] * v[6982]
		+ v[391] * v[6983];
	v[4691] = v[4530];
	v[4531] = v[4531] + v[3482] * v[4557] + v[3484] * v[4558] + v[3486] * v[4559] + v[3488] * v[4560] + v[3490] * v[4561]
		+ v[3492] * v[4562] + v[3494] * v[4563] + v[3496] * v[4564] + v[3498] * v[4565] + v[4637] * v[5461] + v[4570] * v[5481]
		+ v[4568] * v[5483] + v[5515] * v[6680] + v[389] * v[6970] + v[390] * v[6972];
	v[4689] = v[4531];
	b4684 = b6;
	if (b4684) {
		v[4685] = v[375] * v[4531];
		v[4523] = v[4523] + v[388] * v[4531];
		v[4531] = 0e0;
		v[4686] = v[374] * v[4530] + v[4685];
		v[4524] = v[4524] + v[388] * v[4530];
		v[4530] = 0e0;
		v[4687] = v[373] * v[4529] + v[4686];
		v[4525] = v[4525] + v[388] * v[4529];
		v[4529] = 0e0;
		v[4532] = v[4532] + v[387] * v[4687];
	}
	else {
		b4688 = b418;
		if (b4688) {
			v[4690] = -v[4689];
			v[4531] = 0e0;
			v[4692] = -v[4691];
			v[4530] = 0e0;
			v[4694] = -v[4693];
			v[4529] = 0e0;
		}
		else {
			v[4690] = v[4689];
			v[4692] = v[4691];
			v[4694] = v[4693];
		};
		v[4523] = v[407] * v[4690] + v[4696];
		v[4524] = v[407] * v[4692] + v[4698];
		v[4699] = v[375] * v[4690] + v[374] * v[4692] + v[373] * v[4694];
		v[4525] = v[407] * v[4694] + v[4700];
		v[4532] = v[406] * v[4699] + v[4701];
	};
	v[6932] = v[4532] / v[401];
	v[4523] = v[4523] + v[375] * v[6932];
	v[4524] = v[4524] + v[374] * v[6932];
	v[4525] = v[4525] + v[373] * v[6932];
	v[7439] = v[204] * v[4523] + v[201] * v[4524] + v[198] * v[4525];
	v[7438] = v[203] * v[4523] + v[200] * v[4524] + v[197] * v[4525];
	v[7435] = v[213] * v[4523] + v[210] * v[4524] + v[207] * v[4525];
	v[7434] = v[212] * v[4523] + v[209] * v[4524] + v[206] * v[4525];
	v[4719] = v[224] * v[4523];
	v[4720] = v[223] * v[4523];
	v[4738] = v[224] * v[4524];
	v[4739] = v[223] * v[4524];
	v[6337] = -(v[264] * v[4523]) - v[261] * v[4524] - v[258] * v[4525] + v[4663] * v[739] + v[4660] * v[740]
		+ v[4657] * v[741] + v[4664] * v[754] + v[4661] * v[755] + v[4658] * v[756] + v[4665] * v[769] + v[4662] * v[770]
		+ v[4659] * v[771];
	v[6363] = v[1404] * v[4523] + v[1399] * v[4524] + v[1394] * v[4525] + v[224] * v[7002] + v[223] * v[7015];
	v[6364] = -(v[1403] * v[4523]) - v[1398] * v[4524] - v[1393] * v[4525] + v[224] * (v[6993] + v[6994] + v[6995] + v[6996]
		+ v[6997] + v[6998] + v[6999] + v[7000] + v[7001]) + v[223] * (v[7006] + v[7007] + v[7008] + v[7009] + v[7010] + v[7011]
			+ v[7012] + v[7013] + v[7014]);
	v[4757] = v[224] * v[4525];
	v[4758] = v[223] * v[4525];
	v[4759] = v[3215] * v[4639];
	v[4760] = v[3214] * v[4639];
	v[4761] = v[3213] * v[4639];
	v[4762] = v[3210] * v[4640];
	v[4763] = v[359] * v[4639] + v[361] * v[4640];
	v[6933] = v[243] * v[4763];
	v[6358] = -(v[2319] * v[6933]);
	v[6355] = v[4763] * v[6566];
	v[4764] = v[3208] * v[4640];
	v[4765] = v[4759] + v[4762];
	v[4766] = v[3209] * v[4640] + v[6933] / v[350];
	v[4767] = v[3203] * v[4641];
	v[4768] = v[353] * v[4639] + v[361] * v[4641];
	v[6937] = v[248] * v[4768];
	v[6357] = -(v[2319] * v[6937]);
	v[6354] = v[4768] * v[6567];
	v[4769] = v[353] * v[4640] + v[359] * v[4641];
	v[6938] = v[252] * v[4769];
	v[7501] = -(v[4639] * v[6350]) - v[4640] * v[6351] - v[4641] * v[6352] + v[6472] * v[6933] + v[349] * v[6937]
		+ v[348] * v[6938];
	v[6356] = -(v[2319] * v[6938]);
	v[6353] = v[4769] * v[6569];
	v[6348] = v[3212] * v[4639] + v[3207] * v[4640] + v[3202] * v[4641] + v[4761] + v[2955] * v[4763] + v[4764] + v[4767]
		+ v[2953] * v[4768] + v[2949] * v[4769];
	v[4770] = v[361] * v[7238] + v[4523] * v[981];
	v[4773] = v[353] * v[6934] + v[4523] * v[983];
	v[4774] = v[359] * v[7234] + v[4524] * v[982];
	v[4777] = v[353] * v[6935] + v[4524] * v[983];
	v[4778] = v[2520] * v[4763] + v[361] * v[7240] + v[4525] * v[981];
	v[4780] = v[359] * v[6936] + v[4525] * v[982];
	v[4783] = v[353] * v[7231] + v[4525] * v[983];
	v[4784] = v[3205] * v[4641] + v[6937] / v[350];
	v[4786] = v[4768] * v[6473] + v[361] * v[7239] + v[4524] * v[981];
	v[4787] = v[4766] + v[4784];
	v[4788] = v[3204] * v[4641] + v[6938] / v[350];
	v[4790] = v[4769] * v[6476] + v[359] * v[7235] + v[4523] * v[982];
	v[4791] = v[4760] + v[4788];
	v[4792] = v[3200] * v[4642];
	v[4793] = v[3199] * v[4642];
	v[4794] = v[3198] * v[4642];
	v[4795] = v[3195] * v[4643];
	v[4796] = v[333] * v[4642] + v[335] * v[4643];
	v[6939] = v[183] * v[4796];
	v[6384] = -(v[2456] * v[6939]);
	v[6381] = v[4796] * v[6596];
	v[7488] = v[335] * (v[6232] + v[6234]) + v[6381];
	v[4797] = v[3193] * v[4643];
	v[4798] = v[4792] + v[4795];
	v[4799] = v[3194] * v[4643] + v[6939] / v[324];
	v[4800] = v[3188] * v[4644];
	v[4801] = v[327] * v[4642] + v[335] * v[4644];
	v[6943] = v[188] * v[4801];
	v[6383] = -(v[2456] * v[6943]);
	v[6380] = v[4801] * v[6597];
	v[7485] = v[6380] + v[335] * v[7159];
	v[4802] = v[327] * v[4643] + v[333] * v[4644];
	v[6944] = v[192] * v[4802];
	v[7505] = -(v[4642] * v[6376]) - v[4643] * v[6377] - v[4644] * v[6378] + v[6464] * v[6939] + v[323] * v[6943]
		+ v[322] * v[6944];
	v[6382] = -(v[2456] * v[6944]);
	v[6379] = v[4802] * v[6599];
	v[7483] = v[6379] + v[333] * v[7158];
	v[6374] = v[3197] * v[4642] + v[3192] * v[4643] + v[3187] * v[4644] + v[4794] + v[2927] * v[4796] + v[4797] + v[4800]
		+ v[2925] * v[4801] + v[2921] * v[4802];
	v[4803] = v[335] * v[7226] + v[4523] * v[984];
	v[4806] = v[327] * v[6940] + v[4523] * v[986];
	v[4807] = v[333] * v[7222] + v[4524] * v[985];
	v[4810] = v[327] * v[6941] + v[4524] * v[986];
	v[4811] = v[2512] * v[4796] + v[335] * v[7228] + v[4525] * v[984];
	v[4813] = v[333] * v[6942] + v[4525] * v[985];
	v[4816] = v[327] * v[7219] + v[4525] * v[986];
	v[4817] = v[3190] * v[4644] + v[6943] / v[324];
	v[4819] = v[4801] * v[6465] + v[335] * v[7227] + v[4524] * v[984];
	v[4820] = v[4799] + v[4817];
	v[4821] = v[3189] * v[4644] + v[6944] / v[324];
	v[4823] = v[4802] * v[6468] + v[333] * v[7223] + v[4523] * v[985];
	v[4824] = v[4793] + v[4821];
	v[4825] = v[3185] * v[4645];
	v[4826] = v[3184] * v[4645];
	v[4827] = v[3183] * v[4645];
	v[4828] = v[3180] * v[4646];
	v[4829] = v[307] * v[4645] + v[309] * v[4646];
	v[6945] = v[164] * v[4829];
	v[6407] = -(v[2470] * v[6945]);
	v[6404] = v[4829] * v[6614];
	v[7496] = v[309] * (v[6292] + v[6294]) + v[6404];
	v[4830] = v[3178] * v[4646];
	v[4831] = v[4825] + v[4828];
	v[4832] = v[3179] * v[4646] + v[6945] / v[298];
	v[4833] = v[3173] * v[4647];
	v[4834] = v[301] * v[4645] + v[309] * v[4647];
	v[6949] = v[169] * v[4834];
	v[6406] = -(v[2470] * v[6949]);
	v[6403] = v[4834] * v[6615];
	v[7493] = v[6403] + v[309] * v[7162];
	v[4835] = v[301] * v[4646] + v[307] * v[4647];
	v[6950] = v[173] * v[4835];
	v[7509] = -(v[4645] * v[6399]) - v[4646] * v[6400] - v[4647] * v[6401] + v[6456] * v[6945] + v[297] * v[6949]
		+ v[296] * v[6950];
	v[6405] = -(v[2470] * v[6950]);
	v[6402] = v[4835] * v[6617];
	v[7491] = v[6402] + v[307] * v[7161];
	v[6397] = v[3182] * v[4645] + v[3177] * v[4646] + v[3172] * v[4647] + v[4827] + v[2899] * v[4829] + v[4830] + v[4833]
		+ v[2897] * v[4834] + v[2893] * v[4835];
	v[4836] = v[309] * v[7214] + v[4523] * v[987];
	v[4839] = v[301] * v[6946] + v[4523] * v[989];
	v[4840] = v[307] * v[7210] + v[4524] * v[988];
	v[4843] = v[301] * v[6947] + v[4524] * v[989];
	v[4844] = v[2504] * v[4829] + v[309] * v[7216] + v[4525] * v[987];
	v[4846] = v[307] * v[6948] + v[4525] * v[988];
	v[4849] = v[301] * v[7207] + v[4525] * v[989];
	v[4850] = v[3175] * v[4647] + v[6949] / v[298];
	v[4852] = v[4834] * v[6457] + v[309] * v[7215] + v[4524] * v[987];
	v[4853] = v[4832] + v[4850];
	v[4854] = v[3174] * v[4647] + v[6950] / v[298];
	v[4856] = v[4835] * v[6460] + v[307] * v[7211] + v[4523] * v[988];
	v[4857] = v[4826] + v[4854];
	v[4858] = -(v[4663] * v[982]);
	v[4859] = -(v[4663] * v[981]);
	v[4860] = -(v[4663] * v[983]);
	v[4861] = -(v[4660] * v[983]);
	v[4862] = -(v[4660] * v[982]);
	v[4863] = -(v[4660] * v[981]);
	v[4864] = -(v[4657] * v[983]);
	v[4865] = -(v[4657] * v[981]);
	v[4866] = -(v[4657] * v[982]);
	v[4867] = -(v[4664] * v[983]);
	v[4868] = -(v[4664] * v[982]);
	v[4869] = -(v[4664] * v[981]);
	v[4870] = -(v[4661] * v[983]);
	v[4871] = -(v[4661] * v[981]);
	v[4872] = -(v[4661] * v[982]);
	v[4873] = -(v[4658] * v[981]);
	v[4874] = -(v[4658] * v[983]);
	v[4875] = -(v[4658] * v[982]);
	v[4876] = -(v[4665] * v[983]);
	v[4877] = -(v[4665] * v[982]);
	v[4878] = -(v[4665] * v[981]);
	v[7171] = -2e0*v[4759] + 2e0*v[4762] + v[4859] * v[6626] + v[4858] * v[6627] + v[4869] * v[6628] + v[4867] * v[6629]
		+ v[4877] * v[6630] + v[4876] * v[6631] - v[4860] * v[694] - v[4868] * v[713] - v[4878] * v[731];
	v[4879] = -(v[4662] * v[982]);
	v[4880] = -(v[4662] * v[983]);
	v[4881] = -(v[4662] * v[981]);
	v[4882] = -(v[4659] * v[982]);
	v[4883] = -(v[4659] * v[983]);
	v[4884] = -(v[265] * v[4523]) - v[262] * v[4524] - v[259] * v[4525] + v[4663] * v[742] + v[4660] * v[743]
		+ v[4657] * v[744] + v[4664] * v[757] + v[4661] * v[758] + v[4658] * v[759] + v[4665] * v[772] + v[4662] * v[773]
		+ v[4659] * v[774];
	v[4885] = -(v[263] * v[4523]) - v[260] * v[4524] - v[257] * v[4525] + v[4663] * v[736] + v[4660] * v[737]
		+ v[4657] * v[738] + v[4664] * v[751] + v[4661] * v[752] + v[4658] * v[753] + v[4665] * v[766] + v[4662] * v[767]
		+ v[4659] * v[768];
	v[6336] = v[1611] * v[4884] - v[1025] * v[4885];
	v[4886] = -(v[4659] * v[981]);
	v[7172] = -2e0*v[4766] + 2e0*v[4784] + v[4865] * v[6626] + v[4866] * v[6627] + v[4873] * v[6628] + v[4874] * v[6629]
		+ v[4882] * v[6630] + v[4883] * v[6631] - v[4864] * v[694] - v[4875] * v[713] - v[4886] * v[731];
	v[4887] = v[4858] + v[4867] + v[4873] + v[4882] + v[4764] * v[6669] - v[4787] * v[701] - v[4765] * v[704];
	v[4888] = -v[4862] - v[4865] - v[4870] - v[4883] + v[4767] * v[6668] - v[4787] * v[698] + v[4791] * v[704];
	v[4889] = v[238] * v[4780] - v[4858] * v[691] + v[4862] * v[692] - v[4866] * v[693];
	v[4890] = -v[4859] - v[4871] - v[4876] - v[4879] + v[4761] * v[6667] - v[4765] * v[698] + v[4791] * v[701];
	v[4891] = v[238] * v[4778] - v[4859] * v[691] + v[4863] * v[692] - v[4865] * v[693];
	v[4892] = v[238] * v[4777] - v[4867] * v[691] + v[4870] * v[692] - v[4874] * v[693];
	v[4893] = v[4864] + v[4875];
	v[4894] = v[238] * v[4786] - v[4869] * v[691] + v[4871] * v[692] - v[4873] * v[693];
	v[4895] = v[238] * v[4773] - v[4876] * v[691] + v[4880] * v[692] - v[4883] * v[693];
	v[4896] = v[238] * v[4790] - v[4877] * v[691] + v[4879] * v[692] - v[4882] * v[693];
	v[4897] = v[4868] + v[4878];
	v[4898] = v[4861] + v[4881];
	v[12867] = 0e0;
	v[12868] = 0e0;
	v[12869] = 0e0;
	v[12870] = 0e0;
	v[12871] = 0e0;
	v[12872] = 0e0;
	v[12873] = 0e0;
	v[12874] = 0e0;
	v[12875] = 0e0;
	v[12876] = 0e0;
	v[12877] = 0e0;
	v[12878] = 0e0;
	v[12879] = 0e0;
	v[12880] = 0e0;
	v[12881] = 0e0;
	v[12882] = -0.5e0*v[4888] - v[4897];
	v[12883] = v[4887] / 2e0 - v[4898];
	v[12884] = -0.5e0*v[4890] - v[4893];
	v[4899] = v[4863] - v[4866] - v[4869] + v[4874] + v[4877] - v[4880] + v[470] * v[4887] - v[471] * v[4888]
		- v[473] * v[4890] + v[4770] * v[6575] + v[4774] * v[6576] + v[4783] * v[6577] + v[6348] * v[6648] + v[4893] * v[6837]
		+ v[4897] * v[6838] + v[4898] * v[6839] + v[4780] * v[699] + v[4778] * v[705] + v[4777] * v[709] + v[4786] * v[718]
		+ v[4773] * v[722] + v[4790] * v[726];
	v[4900] = v[4672] * v[986];
	v[4901] = v[4672] * v[985];
	v[4902] = v[4672] * v[984];
	v[4903] = v[4669] * v[985];
	v[4904] = v[4669] * v[986];
	v[4905] = v[4669] * v[984];
	v[4906] = v[4666] * v[985];
	v[4907] = v[4666] * v[986];
	v[4908] = v[4666] * v[984];
	v[4909] = v[4681] * v[989];
	v[4910] = v[4681] * v[988];
	v[4911] = v[4681] * v[987];
	v[4912] = v[4678] * v[988];
	v[4913] = v[4678] * v[989];
	v[4914] = v[4678] * v[987];
	v[4915] = v[4675] * v[988];
	v[4916] = v[4675] * v[989];
	v[4917] = v[4675] * v[987];
	v[4918] = v[4673] * v[986];
	v[4919] = v[4673] * v[985];
	v[4920] = v[4673] * v[984];
	v[4921] = v[4670] * v[986];
	v[4922] = v[4670] * v[984];
	v[4923] = v[4670] * v[985];
	v[4924] = v[4667] * v[984];
	v[4925] = v[4667] * v[986];
	v[4926] = v[4667] * v[985];
	v[4927] = v[4682] * v[989];
	v[4928] = v[4682] * v[988];
	v[4929] = v[4682] * v[987];
	v[4930] = v[4679] * v[989];
	v[4931] = v[4679] * v[987];
	v[4932] = v[4679] * v[988];
	v[4933] = v[4676] * v[987];
	v[4934] = v[4676] * v[989];
	v[4935] = v[4676] * v[988];
	v[4936] = v[4674] * v[985];
	v[4937] = v[4674] * v[984];
	v[4938] = v[4674] * v[986];
	v[7187] = -2e0*v[4792] + 2e0*v[4795] - v[4938] * v[523] - v[4919] * v[542] - v[4902] * v[560] + v[4901] * v[6632]
		+ v[4900] * v[6633] + v[4920] * v[6634] + v[4918] * v[6635] + v[4937] * v[6636] + v[4936] * v[6637];
	v[4939] = v[4671] * v[986];
	v[4940] = v[4671] * v[985];
	v[4941] = v[4671] * v[984];
	v[4942] = v[4668] * v[986];
	v[4943] = v[4668] * v[984];
	v[4944] = v[4668] * v[985];
	v[7188] = -2e0*v[4799] + 2e0*v[4817] - v[4942] * v[523] - v[4926] * v[542] - v[4908] * v[560] + v[4906] * v[6632]
		+ v[4907] * v[6633] + v[4924] * v[6634] + v[4925] * v[6635] + v[4943] * v[6636] + v[4944] * v[6637];
	v[4945] = v[4683] * v[988];
	v[4946] = v[4683] * v[987];
	v[4947] = v[4683] * v[989];
	v[7197] = -2e0*v[4825] + 2e0*v[4828] - v[478] * v[4947] - v[4928] * v[497] - v[4911] * v[515] + v[4910] * v[6638]
		+ v[4909] * v[6639] + v[4929] * v[6640] + v[4927] * v[6641] + v[4946] * v[6642] + v[4945] * v[6643];
	v[4948] = v[4680] * v[989];
	v[4949] = v[4680] * v[988];
	v[4950] = v[4680] * v[987];
	v[4951] = v[4677] * v[989];
	v[4952] = v[4677] * v[987];
	v[4953] = v[4677] * v[988];
	v[7198] = -2e0*v[4832] + 2e0*v[4850] - v[478] * v[4951] - v[4935] * v[497] - v[4917] * v[515] + v[4915] * v[6638]
		+ v[4916] * v[6639] + v[4933] * v[6640] + v[4934] * v[6641] + v[4952] * v[6642] + v[4953] * v[6643];
	v[4954] = v[1610] * v[6363] + v[1024] * v[6364];
	v[4955] = (-(v[426] * v[4523]) + v[427] * v[4523] - v[423] * v[4524] + v[424] * v[4524] - v[420] * v[4525]
		+ v[421] * v[4525] + v[288] * v[4648] - v[291] * v[4648] + v[4672] * v[637] + v[4669] * v[638] + v[4666] * v[639]
		- v[4681] * v[640] - v[4678] * v[641] - v[4675] * v[642] + v[4673] * v[649] - v[4634] * v[6499] + v[4635] * v[6499]
		+ v[4670] * v[650] + v[4667] * v[651] - v[4682] * v[652] - v[4679] * v[653] - v[4676] * v[654] + v[4674] * v[661]
		+ v[4671] * v[662] + v[4668] * v[663] - v[4683] * v[664] - v[4680] * v[665] - v[4677] * v[666] + v[287] * v[6951]
		- v[290] * v[6951] + v[286] * v[6952] - v[289] * v[6952]) / 2e0;
	v[4956] = v[4906] + v[4918] + v[4924] + v[4936] - v[4820] * v[530] - v[4798] * v[533] + v[4797] * v[6659];
	v[4957] = -v[4907] - v[4921] - v[4940] - v[4943] - v[4820] * v[527] + v[4824] * v[533] + v[4800] * v[6658];
	v[4958] = v[178] * v[4813] - v[4936] * v[520] + v[4940] * v[521] - v[4944] * v[522];
	v[4959] = -v[4900] - v[4903] - v[4922] - v[4937] - v[4798] * v[527] + v[4824] * v[530] + v[4794] * v[6657];
	v[4960] = v[178] * v[4811] - v[4937] * v[520] + v[4941] * v[521] - v[4943] * v[522];
	v[4961] = v[178] * v[4810] - v[4918] * v[520] + v[4921] * v[521] - v[4925] * v[522];
	v[4962] = v[4926] + v[4942];
	v[4963] = v[178] * v[4819] - v[4920] * v[520] + v[4922] * v[521] - v[4924] * v[522];
	v[4964] = v[178] * v[4806] - v[4900] * v[520] + v[4904] * v[521] - v[4907] * v[522];
	v[4965] = v[178] * v[4823] - v[4901] * v[520] + v[4903] * v[521] - v[4906] * v[522];
	v[4966] = v[4902] + v[4919];
	v[4967] = v[4905] + v[4939];
	v[12921] = 0e0;
	v[12922] = 0e0;
	v[12923] = 0e0;
	v[12924] = 0e0;
	v[12925] = 0e0;
	v[12926] = 0e0;
	v[12927] = 0e0;
	v[12928] = 0e0;
	v[12929] = 0e0;
	v[12930] = -0.5e0*v[4957] - v[4966];
	v[12931] = v[4956] / 2e0 - v[4967];
	v[12932] = -0.5e0*v[4959] - v[4962];
	v[12933] = 0e0;
	v[12934] = 0e0;
	v[12935] = 0e0;
	v[12936] = 0e0;
	v[12937] = 0e0;
	v[12938] = 0e0;
	v[4968] = v[4901] - v[4904] - v[4920] + v[4925] + v[4941] - v[4944] + v[464] * v[4956] - v[465] * v[4957]
		- v[467] * v[4959] + v[4813] * v[528] + v[4811] * v[534] + v[4810] * v[538] + v[4819] * v[547] + v[4806] * v[551]
		+ v[4823] * v[555] + v[4803] * v[6605] + v[4807] * v[6606] + v[4816] * v[6607] + v[6374] * v[6647] + v[4962] * v[6859]
		+ v[4966] * v[6860] + v[4967] * v[6861];
	v[4969] = -(v[485] * v[4853]) - v[4831] * v[488] + v[4915] + v[4927] + v[4933] + v[4945] + v[4830] * v[6655];
	v[4970] = -(v[482] * v[4853]) + v[4857] * v[488] - v[4916] - v[4930] - v[4949] - v[4952] + v[4833] * v[6654];
	v[4971] = v[159] * v[4846] - v[475] * v[4945] + v[476] * v[4949] - v[477] * v[4953];
	v[4972] = -(v[482] * v[4831]) + v[485] * v[4857] - v[4909] - v[4912] - v[4931] - v[4946] + v[4827] * v[6653];
	v[4973] = v[159] * v[4844] - v[475] * v[4946] + v[476] * v[4950] - v[477] * v[4952];
	v[4974] = v[159] * v[4843] - v[475] * v[4927] + v[476] * v[4930] - v[477] * v[4934];
	v[4975] = v[4935] + v[4951];
	v[4976] = v[159] * v[4852] - v[475] * v[4929] + v[476] * v[4931] - v[477] * v[4933];
	v[4977] = v[159] * v[4839] - v[475] * v[4909] + v[476] * v[4913] - v[477] * v[4916];
	v[4978] = v[159] * v[4856] - v[475] * v[4910] + v[476] * v[4912] - v[477] * v[4915];
	v[4979] = v[4911] + v[4928];
	v[4980] = v[4914] + v[4948];
	v[12939] = 0e0;
	v[12940] = 0e0;
	v[12941] = 0e0;
	v[12942] = -0.5e0*v[4970] - v[4979];
	v[12943] = v[4969] / 2e0 - v[4980];
	v[12944] = -0.5e0*v[4972] - v[4975];
	v[12945] = 0e0;
	v[12946] = 0e0;
	v[12947] = 0e0;
	v[12948] = 0e0;
	v[12949] = 0e0;
	v[12950] = 0e0;
	v[12951] = 0e0;
	v[12952] = 0e0;
	v[12953] = 0e0;
	v[12954] = 0e0;
	v[12955] = 0e0;
	v[12956] = 0e0;
	v[4981] = v[483] * v[4846] + v[4844] * v[489] + v[4910] - v[4913] - v[4929] + v[4843] * v[493] + v[4934] + v[4950]
		- v[4953] + v[458] * v[4969] - v[459] * v[4970] - v[461] * v[4972] + v[4852] * v[502] + v[4839] * v[506] + v[4856] * v[510]
		+ v[4836] * v[6623] + v[4840] * v[6624] + v[4849] * v[6625] + v[6397] * v[6646] + v[4975] * v[6873] + v[4979] * v[6874]
		+ v[4980] * v[6875];
	v[4982] = -(v[1613] * v[4884]) + v[1027] * v[4885];
	v[4983] = v[1026] * v[6336] + v[1612] * v[6337];
	v[4984] = v[7171] / 2e0;
	v[4986] = -v[4760] + v[4788] + v[4881] * v[6575] + v[4872] * v[6576] + v[4861] * v[6577] + v[4862] * v[699]
		+ v[4863] * v[705] + v[4870] * v[709] + v[4871] * v[718] + v[4880] * v[722] + v[4879] * v[726];
	v[12525] = 0e0;
	v[12526] = 0e0;
	v[12527] = 0e0;
	v[12528] = 0e0;
	v[12529] = 0e0;
	v[12530] = 0e0;
	v[12531] = 0e0;
	v[12532] = 0e0;
	v[12533] = 0e0;
	v[12534] = 0e0;
	v[12535] = 0e0;
	v[12536] = 0e0;
	v[12537] = 0e0;
	v[12538] = 0e0;
	v[12539] = 0e0;
	v[12540] = -v[7171];
	v[12541] = 2e0*v[4986];
	v[12542] = -v[7172];
	v[4987] = (v[238] * v[4774] - v[4868] * v[691] + v[4872] * v[692] - v[4875] * v[693]) / 2e0;
	v[4988] = v[7172] / 2e0;
	v[7499] = v[469] * v[4984] - v[472] * v[4986] + v[474] * v[4988];
	v[5013] = v[2639] * v[4984] + v[2637] * v[4986] - v[4987] + v[2634] * v[4988] - v[4899] * v[6432];
	v[6423] = v[5013] + (-(v[238] * v[4783]) + v[4860] * v[691] - v[4861] * v[692] + v[4864] * v[693]) / 2e0;
	v[4989] = (v[238] * v[4770] - v[4878] * v[691] + v[4881] * v[692] - v[4886] * v[693]) / 2e0;
	v[6422] = v[4987] - v[4989] + v[6423];
	v[6420] = -v[4989] + v[5013];
	v[4990] = v[4891] + v[4895];
	v[4991] = v[4894] + v[4896];
	v[4992] = v[4889] + v[4892];
	v[4993] = v[7187] / 2e0;
	v[4995] = -v[4793] + v[4821] + v[4940] * v[528] + v[4941] * v[534] + v[4921] * v[538] + v[4922] * v[547] + v[4904] * v[551]
		+ v[4903] * v[555] + v[4905] * v[6605] + v[4923] * v[6606] + v[4939] * v[6607];
	v[12507] = 0e0;
	v[12508] = 0e0;
	v[12509] = 0e0;
	v[12510] = 0e0;
	v[12511] = 0e0;
	v[12512] = 0e0;
	v[12513] = 0e0;
	v[12514] = 0e0;
	v[12515] = 0e0;
	v[12516] = -v[7187];
	v[12517] = 2e0*v[4995];
	v[12518] = -v[7188];
	v[12519] = 0e0;
	v[12520] = 0e0;
	v[12521] = 0e0;
	v[12522] = 0e0;
	v[12523] = 0e0;
	v[12524] = 0e0;
	v[4996] = (v[178] * v[4807] - v[4919] * v[520] + v[4923] * v[521] - v[4926] * v[522]) / 2e0;
	v[4997] = v[7188] / 2e0;
	v[7503] = v[463] * v[4993] - v[466] * v[4995] + v[468] * v[4997];
	v[5012] = v[2613] * v[4993] + v[2611] * v[4995] - v[4996] + v[2608] * v[4997] - v[4968] * v[6430];
	v[6419] = v[5012] + (-(v[178] * v[4816]) + v[4938] * v[520] - v[4939] * v[521] + v[4942] * v[522]) / 2e0;
	v[4998] = (v[178] * v[4803] - v[4902] * v[520] + v[4905] * v[521] - v[4908] * v[522]) / 2e0;
	v[6418] = v[4996] - v[4998] + v[6419];
	v[6416] = -v[4998] + v[5012];
	v[4999] = v[4960] + v[4964];
	v[5000] = v[4963] + v[4965];
	v[5001] = v[4958] + v[4961];
	v[5002] = v[7197] / 2e0;
	v[5004] = -v[4826] + v[4854] + v[493] * v[4930] + v[483] * v[4949] + v[489] * v[4950] + v[4931] * v[502] + v[4913] * v[506]
		+ v[4912] * v[510] + v[4914] * v[6623] + v[4932] * v[6624] + v[4948] * v[6625];
	v[12489] = 0e0;
	v[12490] = 0e0;
	v[12491] = 0e0;
	v[12492] = -v[7197];
	v[12493] = 2e0*v[5004];
	v[12494] = -v[7198];
	v[12495] = 0e0;
	v[12496] = 0e0;
	v[12497] = 0e0;
	v[12498] = 0e0;
	v[12499] = 0e0;
	v[12500] = 0e0;
	v[12501] = 0e0;
	v[12502] = 0e0;
	v[12503] = 0e0;
	v[12504] = 0e0;
	v[12505] = 0e0;
	v[12506] = 0e0;
	v[5005] = (v[159] * v[4840] - v[475] * v[4928] + v[476] * v[4932] - v[477] * v[4935]) / 2e0;
	v[5006] = v[7198] / 2e0;
	v[7507] = v[457] * v[5002] - v[460] * v[5004] + v[462] * v[5006];
	v[5011] = v[2587] * v[5002] + v[2585] * v[5004] - v[5005] + v[2582] * v[5006] - v[4981] * v[6428];
	v[6415] = (-(v[159] * v[4849]) + v[475] * v[4947] - v[476] * v[4948] + v[477] * v[4951]) / 2e0 + v[5011];
	v[5007] = (v[159] * v[4836] - v[475] * v[4911] + v[476] * v[4914] - v[477] * v[4917]) / 2e0;
	v[6414] = v[5005] - v[5007] + v[6415];
	v[6412] = -v[5007] + v[5011];
	v[5008] = v[4973] + v[4977];
	v[5009] = v[4976] + v[4978];
	v[5010] = v[4971] + v[4974];
	v[12197] = v[4758];
	v[12198] = v[4739];
	v[12199] = v[4720];
	v[12200] = -v[4976] + v[4978] + v[461] * v[5008] + v[458] * v[5010] + v[457] * v[6412] + v[4970] * v[6434] + 2e0*
		(v[5002] * v[6428] + v[4979] * v[6434]);
	v[12201] = v[4973] - v[4977] + v[461] * v[5009] + v[459] * v[5010] + v[460] * v[6414] - v[4969] * v[6434] + 2e0*(-
		(v[5004] * v[6428]) + v[4980] * v[6434]);
	v[12202] = -v[4971] + v[4974] + v[459] * v[5008] + v[458] * v[5009] + v[462] * v[6415] + v[4972] * v[6434] + 2e0*
		(v[5006] * v[6428] + v[4975] * v[6434]);
	v[12203] = v[4757];
	v[12204] = v[4738];
	v[12205] = v[4719];
	v[12206] = -v[4963] + v[4965] + v[467] * v[4999] + v[464] * v[5001] + v[463] * v[6416] + v[4957] * v[6441] + 2e0*
		(v[4993] * v[6430] + v[4966] * v[6441]);
	v[12207] = v[4960] - v[4964] + v[467] * v[5000] + v[465] * v[5001] + v[466] * v[6418] - v[4956] * v[6441] + 2e0*(-
		(v[4995] * v[6430]) + v[4967] * v[6441]);
	v[12208] = -v[4958] + v[4961] + v[465] * v[4999] + v[464] * v[5000] + v[468] * v[6419] + v[4959] * v[6441] + 2e0*
		(v[4997] * v[6430] + v[4962] * v[6441]);
	v[12209] = -v[4525];
	v[12210] = -v[4524];
	v[12211] = -v[4523];
	v[12212] = -v[4894] + v[4896] + v[473] * v[4990] + v[470] * v[4992] + v[469] * v[6420] + v[4888] * v[6448] + 2e0*
		(v[4984] * v[6432] + v[4897] * v[6448]);
	v[12213] = v[4891] - v[4895] + v[473] * v[4991] + v[471] * v[4992] + v[472] * v[6422] - v[4887] * v[6448] + 2e0*(-
		(v[4986] * v[6432]) + v[4898] * v[6448]);
	v[12214] = -v[4889] + v[4892] + v[471] * v[4990] + v[470] * v[4991] + v[474] * v[6423] + v[4890] * v[6448] + 2e0*
		(v[4988] * v[6432] + v[4893] * v[6448]);
	for (i4521 = 1; i4521 <= 18; i4521++) {
		b7005 = 3 == i4521;
		i7004 = i4521 == 2;
		i7003 = i4521 == 1;
		b6992 = 9 == i4521;
		i6991 = i4521 == 8;
		i6990 = i4521 == 7;
		i6963 = (i6990 ? 1 : 0);
		i6962 = (i7003 ? 1 : 0);
		i6961 = (i4521 == 13 ? 1 : 0);
		i6960 = (i6991 ? 1 : 0);
		i6959 = (i7004 ? 1 : 0);
		i6958 = (i4521 == 14 ? 1 : 0);
		i6957 = (i4521 == 15 ? 1 : 0);
		v[5124] = (i4521 == 17 ? 1 : 0);
		v[7103] = v[5124] * v[7];
		v[5121] = (i4521 == 16 ? 1 : 0);
		v[7104] = v[5121] * v[7];
		v[5118] = (i4521 == 18 ? 1 : 0);
		v[7108] = v[5118] * v[7];
		v[5098] = (i4521 == 11 ? 1 : 0);
		v[7112] = v[5098] * v[7];
		v[5095] = (i4521 == 10 ? 1 : 0);
		v[7113] = v[5095] * v[7];
		v[5092] = (i4521 == 12 ? 1 : 0);
		v[7117] = v[5092] * v[7];
		v[5072] = (i4521 == 5 ? 1 : 0);
		v[7121] = v[5072] * v[7];
		v[5069] = (i4521 == 4 ? 1 : 0);
		v[7122] = v[5069] * v[7];
		v[5066] = (i4521 == 6 ? 1 : 0);
		v[7126] = v[5066] * v[7];
		v[5050] = v[9174 + i4521];
		v[5048] = v[9156 + i4521];
		v[5043] = v[9138 + i4521];
		v[5042] = v[9120 + i4521];
		v[7164] = v[1026] * v[5042];
		v[5018] = v[9232 + i4521];
		v[5019] = v[9268 + i4521];
		v[5020] = v[9250 + i4521];
		v[5022] = v[10446 + i4521];
		v[5023] = v[9214 + i4521];
		v[5149] = -(v[5023] * v[6428]);
		v[5172] = v[5149] * v[6646];
		v[6964] = v[298] * v[5172];
		v[5079] = -0.5e0*v[5149];
		v[5025] = v[10518 + i4521];
		v[5026] = v[9304 + i4521];
		v[5027] = v[9340 + i4521];
		v[5028] = v[9322 + i4521];
		v[5030] = v[10554 + i4521];
		v[5031] = v[9286 + i4521];
		v[5187] = -(v[5031] * v[6430]);
		v[5210] = v[5187] * v[6647];
		v[6965] = v[324] * v[5210];
		v[5105] = -0.5e0*v[5187];
		v[5033] = v[10626 + i4521];
		v[5034] = v[9376 + i4521];
		v[5035] = v[9412 + i4521];
		v[5036] = v[9394 + i4521];
		v[5038] = v[10662 + i4521];
		v[5039] = v[9358 + i4521];
		v[5243] = -(v[5039] * v[6432]);
		v[5266] = v[5243] * v[6648];
		v[7132] = v[361] * v[5266];
		v[7130] = v[359] * v[5266];
		v[6966] = v[350] * v[5266];
		v[5131] = -0.5e0*v[5243];
		v[5041] = v[10734 + i4521];
		v[5044] = v[449] * v[5042] + v[1027] * v[5043];
		v[5045] = v[451] * v[5042] - v[1613] * v[5043];
		v[5049] = -0.5e0*v[5048];
		v[7081] = -(v[5049] * v[7154]);
		v[7072] = -(v[5049] * v[7153]);
		v[6973] = v[5049] * v[7383];
		v[6971] = v[5049] * v[7384];
		v[5963] = v[4525] * v[5049];
		v[5951] = v[4524] * v[5049];
		v[5939] = v[4523] * v[5049];
		v[5853] = v[4584] * v[5049];
		v[5852] = -(v[5049] * v[6930]);
		v[5850] = v[4583] * v[5049];
		v[5849] = -(v[5049] * v[6929]);
		v[5817] = v[5049] * v[6953];
		v[5816] = v[4582] * v[5049];
		v[5533] = -(v[5049] * v[6499]);
		v[5051] = -(v[1024] * v[5050]);
		v[5052] = v[1610] * v[5050];
		v[5053] = v[1612] * v[5042];
		v[5058] = i6962 * v[7];
		v[5059] = i6959 * v[7];
		v[5060] = i6963 * v[7];
		v[5061] = i6960 * v[7];
		v[5062] = i6961 * v[7];
		v[5063] = i6958 * v[7];
		v[7035] = -(v[4533] * v[5063]);
		v[5064] = i6957 * v[7];
		v[7039] = -(v[4533] * v[5064]);
		v[5065] = v[5018] + v[5066];
		v[7191] = 2e0*v[5065];
		v[5067] = v[5018] - v[5066];
		v[7192] = 2e0*v[5067];
		v[5068] = v[5019] + v[5069];
		v[7193] = 2e0*v[5068];
		v[5070] = v[5019] - v[5069];
		v[7194] = 2e0*v[5070];
		v[5071] = v[5020] - v[5072];
		v[7195] = 2e0*v[5071];
		v[5073] = v[5020] + v[5072];
		v[7196] = 2e0*v[5073];
		v[5074] = v[2582] * v[5023] - v[5066] * v[6649];
		v[5075] = -v[5023] + v[460] * v[5072];
		v[5076] = v[2585] * v[5023] + v[5072] * v[6649];
		v[5077] = v[2587] * v[5023] - v[5069] * v[6649];
		v[6954] = 2e0*(-(v[159] * v[5072]) + v[460] * v[5079]);
		v[5080] = -(v[159] * v[5069]) + v[457] * v[5079];
		v[5081] = -(v[159] * v[5066]) + v[462] * v[5079];
		v[5082] = -(v[461] * v[5149]) + v[5066] * v[6434];
		v[5083] = -(v[459] * v[5149]) + v[5069] * v[6434];
		v[5084] = v[458] * v[5149] - v[5072] * v[6434];
		v[5085] = -(v[5079] * v[515]) - v[5022] * v[6434];
		v[5301] = v[309] * v[5085];
		v[5086] = (-(v[477] * v[5022]) - v[5074] * v[515]) / 2e0;
		v[5087] = -(v[497] * v[5079]) - v[5075] * v[6434];
		v[5297] = v[307] * v[5087];
		v[5088] = (v[476] * v[5075] + v[497] * v[5076]) / 2e0;
		v[5089] = -(v[478] * v[5079]) - v[5025] * v[6434];
		v[5292] = v[301] * v[5089];
		v[5090] = (-(v[475] * v[5025]) - v[478] * v[5077]) / 2e0;
		v[5091] = v[5026] + v[5092];
		v[7181] = 2e0*v[5091];
		v[5093] = v[5026] - v[5092];
		v[7182] = 2e0*v[5093];
		v[5094] = v[5027] + v[5095];
		v[7183] = 2e0*v[5094];
		v[5096] = v[5027] - v[5095];
		v[7184] = 2e0*v[5096];
		v[5097] = v[5028] - v[5098];
		v[7185] = 2e0*v[5097];
		v[5099] = v[5028] + v[5098];
		v[7186] = 2e0*v[5099];
		v[5100] = v[2608] * v[5031] - v[5092] * v[6650];
		v[5101] = -v[5031] + v[466] * v[5098];
		v[5102] = v[2611] * v[5031] + v[5098] * v[6650];
		v[5103] = v[2613] * v[5031] - v[5095] * v[6650];
		v[6955] = 2e0*(-(v[178] * v[5098]) + v[466] * v[5105]);
		v[5106] = -(v[178] * v[5095]) + v[463] * v[5105];
		v[5107] = -(v[178] * v[5092]) + v[468] * v[5105];
		v[5108] = -(v[467] * v[5187]) + v[5092] * v[6441];
		v[5109] = -(v[465] * v[5187]) + v[5095] * v[6441];
		v[5110] = v[464] * v[5187] - v[5098] * v[6441];
		v[5111] = -(v[5105] * v[560]) - v[5030] * v[6441];
		v[5341] = v[335] * v[5111];
		v[5112] = (-(v[5030] * v[522]) - v[5100] * v[560]) / 2e0;
		v[5113] = -(v[5105] * v[542]) - v[5101] * v[6441];
		v[5337] = v[333] * v[5113];
		v[5114] = (v[5101] * v[521] + v[5102] * v[542]) / 2e0;
		v[5115] = -(v[5105] * v[523]) - v[5033] * v[6441];
		v[5332] = v[327] * v[5115];
		v[5116] = (-(v[5033] * v[520]) - v[5103] * v[523]) / 2e0;
		v[5117] = v[5034] + v[5118];
		v[7165] = 2e0*v[5117];
		v[5119] = v[5034] - v[5118];
		v[7166] = 2e0*v[5119];
		v[5120] = v[5035] + v[5121];
		v[7167] = 2e0*v[5120];
		v[5122] = v[5035] - v[5121];
		v[7168] = 2e0*v[5122];
		v[5123] = v[5036] - v[5124];
		v[7169] = 2e0*v[5123];
		v[5125] = v[5036] + v[5124];
		v[7170] = 2e0*v[5125];
		v[5126] = v[2634] * v[5039] - v[5118] * v[6651];
		v[5127] = -v[5039] + v[472] * v[5124];
		v[5128] = v[2637] * v[5039] + v[5124] * v[6651];
		v[5129] = v[2639] * v[5039] - v[5121] * v[6651];
		v[6956] = 2e0*(-(v[238] * v[5124]) + v[472] * v[5131]);
		v[5132] = -(v[238] * v[5121]) + v[469] * v[5131];
		v[5133] = -(v[238] * v[5118]) + v[474] * v[5131];
		v[5134] = -(v[473] * v[5243]) + v[5118] * v[6448];
		v[5135] = -(v[471] * v[5243]) + v[5121] * v[6448];
		v[5136] = v[470] * v[5243] - v[5124] * v[6448];
		v[5137] = -(v[5038] * v[6448]) - v[5131] * v[731];
		v[5381] = v[361] * v[5137];
		v[5138] = (-(v[5038] * v[693]) - v[5126] * v[731]) / 2e0;
		v[5139] = -(v[5127] * v[6448]) - v[5131] * v[713];
		v[5377] = v[359] * v[5139];
		v[5140] = (v[5127] * v[692] + v[5128] * v[713]) / 2e0;
		v[5141] = -(v[5041] * v[6448]) - v[5131] * v[694];
		v[5372] = v[353] * v[5141];
		v[5142] = (-(v[5041] * v[691]) - v[5129] * v[694]) / 2e0;
		v[5143] = (v[476] * v[5025] + v[478] * v[5076] + v[6954]) / 2e0;
		v[5144] = (v[476] * v[5022] + v[5076] * v[515] + v[6954]) / 2e0;
		v[5145] = v[5080] + v[5075] * v[6435] - v[5077] * v[6624];
		v[5146] = v[5080] + v[5022] * v[6435] - v[5077] * v[6623];
		v[5147] = -(v[475] * v[5068]) - v[5077] * v[510] + v[5149];
		v[5148] = v[159] * v[5068] + v[510] * v[5149];
		v[5150] = v[476] * v[5071] + v[506] * v[5076] - v[5149];
		v[5151] = v[159] * v[5071] + v[506] * v[5149];
		v[5303] = v[301] * v[5151];
		v[5152] = -(v[475] * v[5070]) - v[502] * v[5077] - v[5149];
		v[5153] = v[159] * v[5070] + v[502] * v[5149];
		v[5298] = v[309] * v[5153];
		v[5154] = v[5081] + v[5025] * v[6436] - v[5074] * v[6625];
		v[5155] = v[5081] + v[5075] * v[6436] - v[5074] * v[6624];
		v[5156] = -(v[477] * v[5065]) - v[493] * v[5074] + v[5149];
		v[5157] = v[159] * v[5065] + v[493] * v[5149];
		v[5158] = v[476] * v[5073] + v[489] * v[5076] + v[5149];
		v[5159] = v[159] * v[5073] + v[489] * v[5149];
		v[5293] = v[309] * v[5159];
		v[5160] = -(v[475] * v[5073]) - v[489] * v[5077] - v[5082];
		v[5161] = v[476] * v[5070] + v[502] * v[5076] - v[5082];
		v[5162] = v[476] * v[5068] - v[5082] + v[5076] * v[510];
		v[5163] = -(v[475] * v[5071]) - v[506] * v[5077] - v[5082];
		v[5164] = v[5172] + v[5082] * v[6653];
		v[5165] = -(v[477] * v[5067]) - v[483] * v[5074] - v[5149];
		v[5166] = v[159] * v[5067] + v[483] * v[5149];
		v[5167] = -(v[477] * v[5073]) - v[489] * v[5074] - v[5083];
		v[5168] = v[476] * v[5067] + v[483] * v[5076] - v[5083];
		v[5169] = v[476] * v[5065] + v[493] * v[5076] - v[5083];
		v[5170] = -(v[477] * v[5071]) - v[506] * v[5074] - v[5083];
		v[5171] = v[485] * v[5082] + v[488] * v[5083];
		v[5173] = v[5172] + v[5083] * v[6654];
		v[5316] = v[4647] * v[5173];
		v[5174] = -(v[475] * v[5067]) - v[483] * v[5077] + v[5084];
		v[5175] = -(v[477] * v[5070]) - v[502] * v[5074] + v[5084];
		v[5176] = -(v[475] * v[5065]) - v[493] * v[5077] + v[5084];
		v[5177] = -(v[477] * v[5068]) + v[5084] - v[5074] * v[510];
		v[5178] = -(v[482] * v[5083]) - v[485] * v[5084];
		v[5179] = -(v[482] * v[5082]) - v[488] * v[5084];
		v[5180] = v[5172] + v[5084] * v[6655];
		v[5320] = v[4646] * v[5180];
		v[5181] = (v[5033] * v[521] + v[5102] * v[523] + v[6955]) / 2e0;
		v[5182] = (v[5030] * v[521] + v[5102] * v[560] + v[6955]) / 2e0;
		v[5183] = v[5106] + v[5101] * v[6442] - v[5103] * v[6606];
		v[5184] = v[5106] + v[5030] * v[6442] - v[5103] * v[6605];
		v[5185] = v[5187] - v[5094] * v[520] - v[5103] * v[555];
		v[5186] = v[178] * v[5094] + v[5187] * v[555];
		v[5188] = -v[5187] + v[5097] * v[521] + v[5102] * v[551];
		v[5189] = v[178] * v[5097] + v[5187] * v[551];
		v[5343] = v[327] * v[5189];
		v[5190] = -v[5187] - v[5096] * v[520] - v[5103] * v[547];
		v[5191] = v[178] * v[5096] + v[5187] * v[547];
		v[5338] = v[335] * v[5191];
		v[5192] = v[5107] + v[5033] * v[6443] - v[5100] * v[6607];
		v[5193] = v[5107] + v[5101] * v[6443] - v[5100] * v[6606];
		v[5194] = v[5187] - v[5091] * v[522] - v[5100] * v[538];
		v[5195] = v[178] * v[5091] + v[5187] * v[538];
		v[5196] = v[5187] + v[5099] * v[521] + v[5102] * v[534];
		v[5197] = v[178] * v[5099] + v[5187] * v[534];
		v[5333] = v[335] * v[5197];
		v[5198] = -v[5108] - v[5099] * v[520] - v[5103] * v[534];
		v[5199] = -v[5108] + v[5096] * v[521] + v[5102] * v[547];
		v[5200] = -v[5108] + v[5094] * v[521] + v[5102] * v[555];
		v[5201] = -v[5108] - v[5097] * v[520] - v[5103] * v[551];
		v[5202] = v[5210] + v[5108] * v[6657];
		v[5203] = -v[5187] - v[5093] * v[522] - v[5100] * v[528];
		v[5204] = v[178] * v[5093] + v[5187] * v[528];
		v[5205] = -v[5109] - v[5099] * v[522] - v[5100] * v[534];
		v[5206] = -v[5109] + v[5093] * v[521] + v[5102] * v[528];
		v[5207] = -v[5109] + v[5091] * v[521] + v[5102] * v[538];
		v[5208] = -v[5109] - v[5097] * v[522] - v[5100] * v[551];
		v[5209] = v[5108] * v[530] + v[5109] * v[533];
		v[5211] = v[5210] + v[5109] * v[6658];
		v[5356] = v[4644] * v[5211];
		v[5212] = v[5110] - v[5093] * v[520] - v[5103] * v[528];
		v[5213] = v[5110] - v[5096] * v[522] - v[5100] * v[547];
		v[5214] = v[5110] - v[5091] * v[520] - v[5103] * v[538];
		v[5215] = v[5110] - v[5094] * v[522] - v[5100] * v[555];
		v[5216] = -(v[5109] * v[527]) - v[5110] * v[530];
		v[5217] = -(v[5108] * v[527]) - v[5110] * v[533];
		v[5218] = v[5210] + v[5110] * v[6659];
		v[5360] = v[4643] * v[5218];
		v[5219] = v[5050] * v[633] + v[5049] * v[666] + v[5167] * v[987] + v[5165] * v[988] + v[5154] * v[989];
		v[7066] = v[389] * v[5219];
		v[5220] = v[5050] * v[632] + v[5049] * v[665] + v[5158] * v[987] + v[5168] * v[988] + v[5143] * v[989];
		v[7069] = v[389] * v[5220];
		v[5221] = v[5050] * v[631] + v[5049] * v[664] + v[5160] * v[987] + v[5174] * v[988] + v[5090] * v[989];
		v[7023] = v[389] * v[5221];
		v[5222] = v[5050] * v[636] - v[5049] * v[663] + v[5205] * v[984] + v[5203] * v[985] + v[5192] * v[986];
		v[7057] = v[389] * v[5222];
		v[5223] = v[5050] * v[635] - v[5049] * v[662] + v[5196] * v[984] + v[5206] * v[985] + v[5181] * v[986];
		v[7060] = v[389] * v[5223];
		v[5224] = v[5050] * v[634] - v[5049] * v[661] + v[5198] * v[984] + v[5212] * v[985] + v[5116] * v[986];
		v[7063] = v[389] * v[5224];
		v[5225] = v[5050] * v[627] + v[5049] * v[654] + v[5175] * v[987] + v[5155] * v[988] + v[5156] * v[989];
		v[7064] = v[390] * v[5225];
		v[5226] = v[5050] * v[626] + v[5049] * v[653] + v[5161] * v[987] + v[5088] * v[988] + v[5169] * v[989];
		v[7067] = v[390] * v[5226];
		v[5227] = v[5050] * v[625] + v[5049] * v[652] + v[5152] * v[987] + v[5145] * v[988] + v[5176] * v[989];
		v[7019] = v[390] * v[5227];
		v[5228] = v[5050] * v[630] - v[5049] * v[651] + v[5213] * v[984] + v[5193] * v[985] + v[5194] * v[986];
		v[7055] = v[390] * v[5228];
		v[5229] = v[5050] * v[629] - v[5049] * v[650] + v[5199] * v[984] + v[5114] * v[985] + v[5207] * v[986];
		v[7058] = v[390] * v[5229];
		v[5230] = v[5050] * v[628] - v[5049] * v[649] + v[5190] * v[984] + v[5183] * v[985] + v[5214] * v[986];
		v[7061] = v[390] * v[5230];
		v[5231] = v[5050] * v[621] + v[5049] * v[642] + v[5086] * v[987] + v[5177] * v[988] + v[5170] * v[989];
		v[7065] = v[391] * v[5231];
		v[5232] = v[5050] * v[620] + v[5049] * v[641] + v[5144] * v[987] + v[5162] * v[988] + v[5150] * v[989];
		v[7068] = v[391] * v[5232];
		v[5233] = v[5050] * v[619] + v[5049] * v[640] + v[5146] * v[987] + v[5147] * v[988] + v[5163] * v[989];
		v[7017] = v[391] * v[5233];
		v[7021] = v[7017] + v[7019];
		v[5234] = v[5050] * v[624] - v[5049] * v[639] + v[5112] * v[984] + v[5215] * v[985] + v[5208] * v[986];
		v[7056] = v[391] * v[5234];
		v[5235] = v[5050] * v[623] - v[5049] * v[638] + v[5182] * v[984] + v[5200] * v[985] + v[5188] * v[986];
		v[7059] = v[391] * v[5235];
		v[5236] = v[5050] * v[622] - v[5049] * v[637] + v[5184] * v[984] + v[5185] * v[985] + v[5201] * v[986];
		v[7062] = v[391] * v[5236];
		v[5237] = (v[5038] * v[692] + v[6956] + v[5128] * v[731]) / 2e0;
		v[5238] = (v[5041] * v[692] + v[5128] * v[694] + v[6956]) / 2e0;
		v[5239] = v[5132] + v[5038] * v[6449] - v[5129] * v[6575];
		v[5240] = v[5132] + v[5127] * v[6449] - v[5129] * v[6576];
		v[5241] = v[5243] - v[5120] * v[691] - v[5129] * v[726];
		v[5242] = v[238] * v[5120] + v[5243] * v[726];
		v[5244] = -v[5243] + v[5123] * v[692] + v[5128] * v[722];
		v[5245] = v[238] * v[5123] + v[5243] * v[722];
		v[5423] = (b6992 ? v[224] : 0e0) + (b7005 ? v[223] : 0e0) - i6957 - v[263] * v[5044] - v[265] * v[5045] + v[1403] * v[5051]
			+ v[1404] * v[5052] - v[264] * v[5053] - v[5049] * v[6480] + v[5137] * v[981] + v[5242] * v[982] + v[5245] * v[983]
			+ v[5111] * v[984] + v[5186] * v[985] + v[5189] * v[986] + v[5085] * v[987] + v[5148] * v[988] + v[5151] * v[989];
		v[5383] = v[353] * v[5245];
		v[5246] = -v[5243] - v[5122] * v[691] - v[5129] * v[718];
		v[5247] = v[238] * v[5122] + v[5243] * v[718];
		v[5378] = v[361] * v[5247];
		v[5248] = v[5133] + v[5127] * v[6450] - v[5126] * v[6576];
		v[5249] = v[5133] + v[5041] * v[6450] - v[5126] * v[6577];
		v[5250] = v[5243] - v[5117] * v[693] - v[5126] * v[709];
		v[5251] = v[238] * v[5117] + v[5243] * v[709];
		v[5419] = -i6958 + i6959 * v[223] + i6960 * v[224] - v[260] * v[5044] - v[262] * v[5045] + v[1398] * v[5051]
			+ v[1399] * v[5052] - v[261] * v[5053] - v[5049] * v[6479] + v[5247] * v[981] + v[5139] * v[982] + v[5251] * v[983]
			+ v[5191] * v[984] + v[5113] * v[985] + v[5195] * v[986] + v[5153] * v[987] + v[5087] * v[988] + v[5157] * v[989];
		v[5252] = v[5243] + v[5125] * v[692] + v[5128] * v[705];
		v[5253] = v[238] * v[5125] + v[5243] * v[705];
		v[5373] = v[361] * v[5253];
		v[5254] = -v[5134] + v[5120] * v[692] + v[5128] * v[726];
		v[5255] = -v[5134] - v[5123] * v[691] - v[5129] * v[722];
		v[5256] = -v[5134] + v[5122] * v[692] + v[5128] * v[718];
		v[5257] = -v[5134] - v[5125] * v[691] - v[5129] * v[705];
		v[5258] = v[5266] + v[5134] * v[6667];
		v[5259] = -v[5243] - v[5119] * v[693] - v[5126] * v[699];
		v[5260] = v[238] * v[5119] + v[5243] * v[699];
		v[5415] = -i6961 + i6962 * v[223] + i6963 * v[224] - v[257] * v[5044] - v[259] * v[5045] + v[1393] * v[5051]
			+ v[1394] * v[5052] - v[258] * v[5053] - v[5049] * v[6478] + v[5253] * v[981] + v[5260] * v[982] + v[5141] * v[983]
			+ v[5197] * v[984] + v[5204] * v[985] + v[5115] * v[986] + v[5159] * v[987] + v[5166] * v[988] + v[5089] * v[989];
		v[5261] = -v[5135] - v[5123] * v[693] - v[5126] * v[722];
		v[5262] = -v[5135] + v[5117] * v[692] + v[5128] * v[709];
		v[5263] = -v[5135] - v[5125] * v[693] - v[5126] * v[705];
		v[5264] = -v[5135] + v[5119] * v[692] + v[5128] * v[699];
		v[5265] = v[5134] * v[701] + v[5135] * v[704];
		v[5267] = v[5266] + v[5135] * v[6668];
		v[5396] = v[4641] * v[5267];
		v[5268] = v[5136] - v[5120] * v[693] - v[5126] * v[726];
		v[5269] = v[5136] - v[5122] * v[693] - v[5126] * v[718];
		v[5270] = v[5136] - v[5117] * v[691] - v[5129] * v[709];
		v[5271] = v[5136] - v[5119] * v[691] - v[5129] * v[699];
		v[5272] = -(v[5135] * v[698]) - v[5136] * v[701];
		v[5273] = -(v[5134] * v[698]) - v[5136] * v[704];
		v[5274] = v[5266] + v[5136] * v[6669];
		v[5400] = v[4640] * v[5274];
		v[5275] = v[5044] * v[768] + v[5053] * v[771] + v[5045] * v[774] - v[5138] * v[981] - v[5268] * v[982] - v[5261] * v[983];
		v[7046] = -(v[391] * v[5275]);
		v[5276] = v[5044] * v[767] + v[5053] * v[770] + v[5045] * v[773] - v[5237] * v[981] - v[5254] * v[982] - v[5244] * v[983];
		v[7049] = -(v[391] * v[5276]);
		v[5277] = v[5044] * v[766] + v[5053] * v[769] + v[5045] * v[772] - v[5239] * v[981] - v[5241] * v[982] - v[5255] * v[983];
		v[7052] = -(v[391] * v[5277]);
		v[5278] = v[5044] * v[753] + v[5053] * v[756] + v[5045] * v[759] - v[5269] * v[981] - v[5248] * v[982] - v[5250] * v[983];
		v[7047] = -(v[390] * v[5278]);
		v[5279] = v[5044] * v[752] + v[5053] * v[755] + v[5045] * v[758] - v[5256] * v[981] - v[5140] * v[982] - v[5262] * v[983];
		v[7050] = -(v[390] * v[5279]);
		v[5280] = v[5044] * v[751] + v[5053] * v[754] + v[5045] * v[757] - v[5246] * v[981] - v[5240] * v[982] - v[5270] * v[983];
		v[7053] = -(v[390] * v[5280]);
		v[5281] = v[5044] * v[738] + v[5053] * v[741] + v[5045] * v[744] - v[5263] * v[981] - v[5259] * v[982] - v[5249] * v[983];
		v[7048] = -(v[389] * v[5281]);
		v[5282] = v[5044] * v[737] + v[5053] * v[740] + v[5045] * v[743] - v[5252] * v[981] - v[5264] * v[982] - v[5238] * v[983];
		v[7051] = -(v[389] * v[5282]);
		v[5283] = v[5044] * v[736] + v[5053] * v[739] + v[5045] * v[742] - v[5257] * v[981] - v[5271] * v[982] - v[5142] * v[983];
		v[7054] = -(v[389] * v[5283]);
		v[5284] = v[5076] + v[5171];
		v[5314] = v[4647] * v[5284];
		v[5285] = -v[5076] + v[5171];
		v[5286] = (v[296] * v[5148] + v[173] * v[5284] + v[2893] * v[6964]) / v[298];
		v[5287] = v[5074] + v[5178];
		v[5322] = v[4647] * v[5287];
		v[5288] = -v[5074] + v[5178];
		v[5318] = v[4646] * v[5288];
		v[5289] = (v[297] * v[5153] + v[169] * v[5287] + v[2897] * v[6964]) / v[298];
		v[5290] = v[307] * v[5166] + v[5292];
		v[5291] = v[5290] + v[5293] + v[7122];
		v[5294] = v[5292] + v[5293];
		v[5295] = v[301] * v[5157] + v[5297];
		v[5296] = v[5295] + v[5298] + v[7121];
		v[5299] = v[5297] + v[5298];
		v[5300] = v[5301] + v[5303];
		v[5302] = v[307] * v[5148] + v[5301];
		v[5304] = v[5302] + v[5303] + v[7126];
		v[5305] = (v[164] * v[5288] + v[5159] * v[6456] + v[2899] * v[6964]) / v[298];
		v[5306] = v[5316] + v[5318];
		v[7212] = v[5306] / v[298];
		v[5307] = v[5077] + v[5179];
		v[5312] = v[4646] * v[5307];
		v[5308] = -v[5077] + v[5179];
		v[5309] = v[5320] + v[5322];
		v[7208] = v[5309] / v[298];
		v[5310] = v[4645] * v[5164] + v[5312] + v[5314];
		v[5313] = v[5310] - v[5314];
		v[5315] = v[5310] - v[5312];
		v[7209] = v[5315] / v[298];
		v[5317] = v[4645] * v[5285] + v[5316];
		v[5319] = v[5317] + v[5318];
		v[5321] = v[4645] * v[5308] + v[5320];
		v[5323] = v[5321] + v[5322];
		v[5324] = v[5102] + v[5209];
		v[5354] = v[4644] * v[5324];
		v[5325] = -v[5102] + v[5209];
		v[5326] = (v[322] * v[5186] + v[192] * v[5324] + v[2921] * v[6965]) / v[324];
		v[5327] = v[5100] + v[5216];
		v[5362] = v[4644] * v[5327];
		v[5328] = -v[5100] + v[5216];
		v[5358] = v[4643] * v[5328];
		v[5329] = (v[323] * v[5191] + v[188] * v[5327] + v[2925] * v[6965]) / v[324];
		v[5330] = v[333] * v[5204] + v[5332];
		v[5331] = v[5330] + v[5333] + v[7113];
		v[5334] = v[5332] + v[5333];
		v[5335] = v[327] * v[5195] + v[5337];
		v[5336] = v[5335] + v[5338] + v[7112];
		v[5339] = v[5337] + v[5338];
		v[5340] = v[5341] + v[5343];
		v[5342] = v[333] * v[5186] + v[5341];
		v[5344] = v[5342] + v[5343] + v[7117];
		v[5345] = (v[183] * v[5328] + v[5197] * v[6464] + v[2927] * v[6965]) / v[324];
		v[5346] = v[5356] + v[5358];
		v[7224] = v[5346] / v[324];
		v[5347] = v[5103] + v[5217];
		v[5352] = v[4643] * v[5347];
		v[5348] = -v[5103] + v[5217];
		v[5349] = v[5360] + v[5362];
		v[7220] = v[5349] / v[324];
		v[5350] = v[4642] * v[5202] + v[5352] + v[5354];
		v[5353] = v[5350] - v[5354];
		v[5355] = v[5350] - v[5352];
		v[7221] = v[5355] / v[324];
		v[5357] = v[4642] * v[5325] + v[5356];
		v[5359] = v[5357] + v[5358];
		v[5361] = v[4642] * v[5348] + v[5360];
		v[5363] = v[5361] + v[5362];
		v[5364] = v[5128] + v[5265];
		v[5394] = v[4641] * v[5364];
		v[5365] = -v[5128] + v[5265];
		v[5366] = (v[348] * v[5242] + v[252] * v[5364] + v[2949] * v[6966]) / v[350];
		v[5367] = v[5126] + v[5272];
		v[5402] = v[4641] * v[5367];
		v[5368] = -v[5126] + v[5272];
		v[5398] = v[4640] * v[5368];
		v[5369] = (v[349] * v[5247] + v[248] * v[5367] + v[2953] * v[6966]) / v[350];
		v[5370] = v[359] * v[5260] + v[5372];
		v[5371] = v[5370] + v[5373] + v[7104];
		v[5374] = v[5372] + v[5373];
		v[5375] = v[353] * v[5251] + v[5377];
		v[5376] = v[5375] + v[5378] + v[7103];
		v[5379] = v[5377] + v[5378];
		v[5380] = v[5381] + v[5383];
		v[5382] = v[359] * v[5242] + v[5381];
		v[5384] = v[5382] + v[5383] + v[7108];
		v[5385] = (v[243] * v[5368] + v[5253] * v[6472] + v[2955] * v[6966]) / v[350];
		v[5386] = v[5396] + v[5398];
		v[7236] = v[5386] / v[350];
		v[5387] = v[5129] + v[5273];
		v[5392] = v[4640] * v[5387];
		v[5388] = -v[5129] + v[5273];
		v[5389] = v[5400] + v[5402];
		v[7232] = v[5389] / v[350];
		v[5390] = v[4639] * v[5258] + v[5392] + v[5394];
		v[5393] = v[5390] - v[5394];
		v[5395] = v[5390] - v[5392];
		v[7233] = v[5395] / v[350];
		v[5397] = v[4639] * v[5365] + v[5396];
		v[5399] = v[5397] + v[5398];
		v[5401] = v[4639] * v[5388] + v[5400];
		v[5403] = v[5401] + v[5402];
		v[5404] = v[4525] * v[5089] + v[4683] * v[5090] + v[4680] * v[5143] + v[4678] * v[5150] + v[4523] * v[5151]
			+ v[4677] * v[5154] + v[4676] * v[5156] + v[4524] * v[5157] + v[4681] * v[5163] + v[4679] * v[5169] + v[4675] * v[5170]
			+ v[4682] * v[5176];
		v[7175] = v[223] * v[5404];
		v[5405] = v[4524] * v[5087] + v[4679] * v[5088] + v[4682] * v[5145] + v[4681] * v[5147] + v[4523] * v[5148]
			+ v[4676] * v[5155] + v[4678] * v[5162] + v[4677] * v[5165] + v[4525] * v[5166] + v[4680] * v[5168] + v[4683] * v[5174]
			+ v[4675] * v[5177];
		v[7176] = v[223] * v[5405];
		v[5406] = v[4523] * v[5085] + v[4675] * v[5086] + v[4678] * v[5144] + v[4681] * v[5146] + v[4682] * v[5152]
			+ v[4524] * v[5153] + v[4680] * v[5158] + v[4525] * v[5159] + v[4683] * v[5160] + v[4679] * v[5161] + v[4677] * v[5167]
			+ v[4676] * v[5175];
		v[7177] = v[223] * v[5406];
		v[5407] = v[4525] * v[5115] + v[4674] * v[5116] + v[4671] * v[5181] + v[4669] * v[5188] + v[4523] * v[5189]
			+ v[4668] * v[5192] + v[4667] * v[5194] + v[4524] * v[5195] + v[4672] * v[5201] + v[4670] * v[5207] + v[4666] * v[5208]
			+ v[4673] * v[5214];
		v[7178] = v[224] * v[5407];
		v[5408] = v[4524] * v[5113] + v[4670] * v[5114] + v[4673] * v[5183] + v[4672] * v[5185] + v[4523] * v[5186]
			+ v[4667] * v[5193] + v[4669] * v[5200] + v[4668] * v[5203] + v[4525] * v[5204] + v[4671] * v[5206] + v[4674] * v[5212]
			+ v[4666] * v[5215];
		v[7179] = v[224] * v[5408];
		v[5409] = v[4523] * v[5111] + v[4666] * v[5112] + v[4669] * v[5182] + v[4672] * v[5184] + v[4673] * v[5190]
			+ v[4524] * v[5191] + v[4671] * v[5196] + v[4525] * v[5197] + v[4674] * v[5198] + v[4670] * v[5199] + v[4668] * v[5205]
			+ v[4667] * v[5213];
		v[7180] = v[224] * v[5409];
		v[5410] = v[4525] * v[5141] - v[4663] * v[5142] - v[4660] * v[5238] - v[4662] * v[5244] + v[4523] * v[5245]
			- v[4657] * v[5249] - v[4658] * v[5250] + v[4524] * v[5251] - v[4665] * v[5255] - v[4659] * v[5261] - v[4661] * v[5262]
			- v[4664] * v[5270];
		v[5411] = v[4524] * v[5139] - v[4661] * v[5140] - v[4664] * v[5240] - v[4665] * v[5241] + v[4523] * v[5242]
			- v[4658] * v[5248] - v[4662] * v[5254] - v[4657] * v[5259] + v[4525] * v[5260] - v[4660] * v[5264] - v[4659] * v[5268]
			- v[4663] * v[5271];
		v[5412] = v[4523] * v[5137] - v[4659] * v[5138] - v[4662] * v[5237] - v[4665] * v[5239] - v[4664] * v[5246]
			+ v[4524] * v[5247] - v[4660] * v[5252] + v[4525] * v[5253] - v[4661] * v[5256] - v[4663] * v[5257] - v[4657] * v[5263]
			- v[4658] * v[5269];
		v[5414] = v[5415] * v[6932];
		v[5447] = v[5414];
		v[5416] = v[5415];
		v[5445] = v[5416];
		v[5418] = v[5419] * v[6932];
		v[5450] = v[5418];
		v[5420] = v[5419];
		v[6967] = v[373] * v[5416] + v[374] * v[5420] + v[375] * v[5423];
		v[5448] = v[5420];
		v[5422] = v[6967] / v[401];
		v[5443] = v[5422];
		v[5424] = -(v[1406] * v[4532] * v[6967]);
		v[5774] = v[5424];
		v[5425] = v[5423] * v[6932];
		v[5453] = v[5425];
		v[5426] = v[5423];
		v[5451] = v[5426];
		v[5427] = 0e0;
		v[5428] = 0e0;
		v[5429] = 0e0;
		v[5430] = 0e0;
		v[5431] = 0e0;
		v[5432] = 0e0;
		v[5433] = 0e0;
		v[5434] = 0e0;
		v[5435] = 0e0;
		v[5436] = 0e0;
		b5437 = b6;
		if (b5437) {
			v[5438] = 0e0;
			v[5439] = 0e0;
			v[5440] = 0e0;
			v[5441] = 0e0;
			v[5442] = v[387] * v[5422];
			v[5427] = v[386] * v[4687] * v[5422];
			v[5431] = v[388] * v[5416];
			v[5428] = v[4693] * v[5416];
			v[5431] = v[5431] + v[373] * v[5442];
			v[5414] = v[5414] + v[4693] * v[5442];
			v[5432] = v[388] * v[5420];
			v[5428] = v[4691] * v[5420] + v[5428];
			v[5432] = v[5432] + v[374] * v[5442];
			v[5418] = v[5418] + v[4691] * v[5442];
			v[5433] = v[388] * v[5426];
			v[5428] = v[4689] * v[5426] + v[5428];
			v[5433] = v[5433] + v[375] * v[5442];
			v[5425] = v[5425] + v[4689] * v[5442];
		}
		else {
			v[5444] = v[406] * v[5443];
			v[5441] = v[5443];
			v[5429] = v[405] * v[4699] * v[5443];
			v[5422] = 0e0;
			v[5440] = v[5445];
			v[5416] = 0e0;
			v[5446] = v[373] * v[5444] + v[407] * v[5445];
			v[5414] = v[4694] * v[5444] + v[5447];
			v[5439] = v[5448];
			v[5420] = 0e0;
			v[5449] = v[374] * v[5444] + v[407] * v[5448];
			v[5418] = v[4692] * v[5444] + v[5450];
			v[5438] = v[5451];
			v[5430] = v[4694] * v[5445] + v[4692] * v[5448] + v[4690] * v[5451];
			v[5426] = 0e0;
			v[5452] = v[375] * v[5444] + v[407] * v[5451];
			v[5425] = v[4690] * v[5444] + v[5453];
			b5454 = b418;
			if (b5454) {
				v[5431] = 0e0;
				v[5436] = -v[5446];
				v[5432] = 0e0;
				v[5435] = -v[5449];
				v[5433] = 0e0;
				v[5434] = -v[5452];
			}
			else {
				v[5436] = v[5446];
				v[5435] = v[5449];
				v[5434] = v[5452];
			};
		};
		v[5929] = v[5414];
		v[5928] = v[5418];
		v[5927] = v[5425];
		v[5433] = v[5433] + v[5434];
		v[7034] = -(v[294] * v[5433]);
		v[6969] = v[390] * v[5433];
		v[6968] = v[389] * v[5433];
		v[5462] = v[5433] * v[5461] + v[391] * v[6973] + v[12398 + i4521] * v[7];
		v[5482] = v[5433] * v[5481];
		v[5484] = v[5433] * v[5483];
		v[5485] = v[5433] * v[6680];
		v[5810] = v[4535] * v[5485];
		v[5486] = v[5433] * v[7396];
		v[5487] = v[5433] * v[7397];
		v[5488] = v[5433] * v[7398];
		v[5489] = v[5433] * v[7399];
		v[5490] = v[5433] * v[7400];
		v[5491] = v[5433] * v[7401];
		v[5492] = v[5433] * v[7402];
		v[5493] = v[5433] * v[7403];
		v[5494] = v[5433] * v[7404];
		v[5495] = v[4565] * v[5433];
		v[5496] = v[4564] * v[5433];
		v[5497] = v[4563] * v[5433];
		v[5498] = v[4562] * v[5433];
		v[5499] = v[4561] * v[5433];
		v[5500] = v[4560] * v[5433];
		v[5501] = v[4559] * v[5433];
		v[5502] = v[4558] * v[5433];
		v[5503] = v[4557] * v[5433];
		v[5504] = v[5433] * v[6970] + v[390] * v[6971];
		v[5514] = v[389] * v[6971] + v[5433] * v[6972];
		v[5516] = 2e0*v[5433] * v[5515] + v[4637] * v[6973] + v[12470 + i4521] * v[7];
		v[5517] = v[4637] * v[5433];
		v[5518] = v[4570] * v[5433];
		v[5519] = v[4568] * v[5433];
		v[5523] = v[5433] * v[7406];
		v[5524] = v[5433] * v[7407];
		v[5525] = v[5433] * v[7408];
		v[5526] = v[5433] * v[7409];
		v[5527] = v[5433] * v[7410];
		v[5528] = v[5433] * v[7411];
		v[5529] = v[5433] * v[7412];
		v[5530] = v[5433] * v[7413];
		v[5531] = v[5433] * v[7414];
		v[5432] = v[5432] + v[5435];
		v[6981] = v[391] * v[5432];
		v[6980] = v[389] * v[5432];
		v[5589] = v[4591] * v[5432];
		v[5587] = v[4592] * v[5432];
		v[5585] = v[4593] * v[5432];
		v[5583] = v[4594] * v[5432];
		v[5581] = v[4595] * v[5432];
		v[5579] = v[4596] * v[5432];
		v[5577] = v[4597] * v[5432];
		v[5575] = v[4598] * v[5432];
		v[5573] = v[4599] * v[5432];
		v[5544] = v[5432] * v[6683];
		v[5812] = v[4534] * v[5544];
		v[5545] = v[685] * v[6969] + v[5432] * v[7415];
		v[5546] = v[686] * v[6969] + v[5432] * v[7416];
		v[5547] = v[687] * v[6969] + v[5432] * v[7417];
		v[5548] = v[688] * v[6969] + v[5432] * v[7418];
		v[5549] = v[689] * v[6969] + v[5432] * v[7419];
		v[5550] = v[690] * v[6969] + v[5432] * v[7420];
		v[5551] = v[5432] * v[7421] - v[6969] * v[787];
		v[5552] = v[5432] * v[7422] - v[6969] * v[788];
		v[5553] = v[5432] * v[7423] - v[6969] * v[789];
		v[5563] = v[5495] + v[5573];
		v[5564] = v[5496] + v[5575];
		v[5565] = v[5497] + v[5577];
		v[5566] = v[5498] + v[5579];
		v[5567] = v[5499] + v[5581];
		v[5568] = v[5500] + v[5583];
		v[5569] = v[5501] + v[5585];
		v[5570] = v[5502] + v[5587];
		v[5571] = v[5503] + v[5589];
		v[5504] = v[5504] + v[5432] * v[6982];
		v[5514] = v[5514] + 2e0*v[5432] * v[5572];
		v[5516] = v[5516] + v[5432] * v[6983];
		v[5592] = v[5432] * v[7424];
		v[5594] = v[5432] * v[7425];
		v[5596] = v[5432] * v[7426];
		v[5598] = v[5432] * v[7427];
		v[5600] = v[5432] * v[7428];
		v[5602] = v[5432] * v[7429];
		v[5604] = v[5432] * v[7430];
		v[5606] = v[5432] * v[7431];
		v[5608] = v[5432] * v[7432];
		v[5431] = v[5431] + v[5436];
		v[7020] = v[391] * v[5431];
		v[7018] = v[5431] * v[673] + v[7023];
		v[7016] = v[390] * v[5431];
		v[5671] = v[4623] * v[5431];
		v[5669] = v[4624] * v[5431];
		v[5667] = v[4625] * v[5431];
		v[5665] = v[4626] * v[5431];
		v[5663] = v[4627] * v[5431];
		v[5661] = v[4628] * v[5431];
		v[5659] = v[4629] * v[5431];
		v[5657] = v[4630] * v[5431];
		v[5655] = v[4631] * v[5431];
		v[5611] = v[5533] + v[224] * (v[6980] + v[7016]);
		v[5612] = -v[5533] + v[223] * (v[6980] + v[7016]);
		v[5613] = v[5431] * v[6686];
		v[5721] = v[4533] * v[5613];
		v[5614] = v[1014] * v[5221] + v[685] * v[6968] + v[679] * v[6980] + v[389] * v[7021] + v[5431] * v[7441];
		v[5615] = v[686] * v[6968] + v[680] * v[6980] + v[5431] * v[7442];
		v[5616] = v[687] * v[6968] + v[681] * v[6980] + v[5431] * v[7443];
		v[5617] = v[688] * v[6968] + v[682] * v[6980] + v[5431] * v[7444];
		v[5618] = v[689] * v[6968] + v[683] * v[6980] + v[5431] * v[7445];
		v[5619] = v[690] * v[6968] + v[684] * v[6980] + v[5431] * v[7446];
		v[5620] = v[5431] * v[7447] - v[6980] * v[784] - v[6968] * v[787];
		v[5621] = v[5431] * v[7448] - v[6980] * v[785] - v[6968] * v[788];
		v[5622] = v[5431] * v[7449] - v[6980] * v[786] - v[6968] * v[789];
		v[5504] = v[5504] + 2e0*v[5431] * v[5643] + v[4631] * v[7021] + v[5221] * v[7022];
		v[5644] = v[5495] + v[5655];
		v[5645] = v[5496] + v[5657];
		v[5646] = v[5497] + v[5659];
		v[5647] = v[5498] + v[5661];
		v[5648] = v[5499] + v[5663];
		v[5649] = v[5500] + v[5665];
		v[5650] = v[5501] + v[5667];
		v[5651] = v[5502] + v[5669];
		v[5652] = v[5503] + v[5671];
		v[5514] = v[5514] + v[5227] * v[6987] + v[4599] * (v[7017] + v[7023]) + v[5431] * v[7024];
		v[5656] = v[5573] + v[5655];
		v[5658] = v[5575] + v[5657];
		v[5660] = v[5577] + v[5659];
		v[5662] = v[5579] + v[5661];
		v[5664] = v[5581] + v[5663];
		v[5666] = v[5583] + v[5665];
		v[5668] = v[5585] + v[5667];
		v[5670] = v[5587] + v[5669];
		v[5672] = v[5589] + v[5671];
		v[5516] = v[5516] + v[5233] * v[6977] + v[4565] * (v[7019] + v[7023]) + v[5431] * v[7025];
		v[5675] = v[5431] * v[7450];
		v[5678] = v[5431] * v[7451];
		v[5681] = v[5431] * v[7452];
		v[5684] = v[5431] * v[7453];
		v[5687] = v[5431] * v[7454];
		v[5690] = v[5431] * v[7455];
		v[5691] = v[5431] * v[7456];
		v[5694] = v[5431] * v[7457];
		v[5697] = v[5431] * v[7458];
		v[5700] = v[3172] * v[5172] + v[3173] * v[5173] + v[3174] * v[5284] + v[307] * v[5286] + v[3175] * v[5287]
			+ v[309] * v[5289] + v[4497] * v[5291] + v[5295] * v[6457] + v[5300] * v[6460] + v[12236 + i4521] * v[7];
		v[5701] = v[3177] * v[5172] + v[3178] * v[5180] + v[301] * v[5286] + v[3179] * v[5288] + v[2504] * v[5290]
			+ v[4500] * v[5296] + v[2503] * v[5302] + v[309] * v[5305] + v[3180] * v[5307] + v[12254 + i4521] * v[7];
		v[5702] = v[3183] * v[5164] + v[3182] * v[5172] + v[3184] * v[5285] + v[301] * v[5289] + v[2506] * v[5294]
			+ v[2505] * v[5299] + v[4502] * v[5304] + v[307] * v[5305] + v[3185] * v[5308] + v[12272 + i4521] * v[7];
		v[5703] = v[3187] * v[5210] + v[3188] * v[5211] + v[3189] * v[5324] + v[333] * v[5326] + v[3190] * v[5327]
			+ v[335] * v[5329] + v[4504] * v[5331] + v[5335] * v[6465] + v[5340] * v[6468] + v[12290 + i4521] * v[7];
		v[5704] = v[3192] * v[5210] + v[3193] * v[5218] + v[327] * v[5326] + v[3194] * v[5328] + v[2512] * v[5330]
			+ v[4507] * v[5336] + v[2511] * v[5342] + v[335] * v[5345] + v[3195] * v[5347] + v[12308 + i4521] * v[7];
		v[5705] = v[3198] * v[5202] + v[3197] * v[5210] + v[3199] * v[5325] + v[327] * v[5329] + v[2514] * v[5334]
			+ v[2513] * v[5339] + v[4509] * v[5344] + v[333] * v[5345] + v[3200] * v[5348] + v[12326 + i4521] * v[7];
		v[5706] = v[3202] * v[5266] + v[3203] * v[5267] + v[3204] * v[5364] + v[359] * v[5366] + v[3205] * v[5367]
			+ v[361] * v[5369] + v[4511] * v[5371] + v[5375] * v[6473] + v[5380] * v[6476] + v[12344 + i4521] * v[7];
		v[5707] = v[3207] * v[5266] + v[3208] * v[5274] + v[353] * v[5366] + v[3209] * v[5368] + v[2520] * v[5370]
			+ v[4514] * v[5376] + v[2519] * v[5382] + v[361] * v[5385] + v[3210] * v[5387] + v[12362 + i4521] * v[7];
		v[5708] = v[3213] * v[5258] + v[3212] * v[5266] + v[3214] * v[5365] + v[353] * v[5369] + v[2522] * v[5374]
			+ v[2521] * v[5379] + v[4516] * v[5384] + v[359] * v[5385] + v[3215] * v[5388] + v[12380 + i4521] * v[7];
		v[5504] = v[5504] + v[5224] * v[7026] + v[5223] * v[7027] + v[5222] * v[7028] + v[5220] * v[7029] + v[5219] * v[7030]
			- v[5283] * v[7031] - v[5282] * v[7032] - v[5281] * v[7033] + v[4533] * v[7034] + v[390] * (v[4629] * v[5225]
				+ v[4630] * v[5226] + v[4626] * v[5228] + v[4627] * v[5229] + v[4628] * v[5230] - v[4623] * v[5278] - v[4624] * v[5279]
				- v[4625] * v[5280] + v[7035]) + v[391] * (v[4629] * v[5231] + v[4630] * v[5232] + v[4626] * v[5234] + v[4627] * v[5235]
					+ v[4628] * v[5236] - v[4623] * v[5275] - v[4624] * v[5276] - v[4625] * v[5277] + v[7039]);
		v[5504] = v[4533] * (-(v[293] * v[5432]) + v[5462]) + v[5504];
		v[5711] = v[4533] * v[5431];
		v[5713] = v[5519] + v[5721];
		v[5714] = v[5062] + v[12416 + i4521] * v[7] + v[7072];
		v[5715] = v[5063] + v[12434 + i4521] * v[7] + v[7081];
		v[5724] = v[4535] * v[5433];
		v[5718] = v[4534] * v[5432];
		v[5719] = v[5711] + v[5718];
		v[5720] = v[5518] + v[5812];
		v[5514] = v[391] * (v[4597] * v[5231] + v[4598] * v[5232] + v[4594] * v[5234] + v[4595] * v[5235] + v[4596] * v[5236]
			- v[4591] * v[5275] - v[4592] * v[5276] - v[4593] * v[5277]) + v[5514] + v[5230] * v[6984] + v[5229] * v[6985]
			+ v[5228] * v[6986] + v[5226] * v[6988] + v[5225] * v[6989] + v[4534] * (-(v[389] * v[5062]) - v[391] * v[5064]
				- v[292] * v[5431] + v[5462] + v[7034]) + v[389] * (v[4597] * v[5219] + v[4598] * v[5220] + v[4594] * v[5222]
					+ v[4595] * v[5223] + v[4596] * v[5224] - v[4591] * v[5281] - v[4592] * v[5282] - v[4593] * v[5283] + v[7035])
			- v[5280] * v[7036] - v[5279] * v[7037] - v[5278] * v[7038];
		v[5722] = -v[5721] - v[5431] * v[6921];
		v[5516] = v[390] * (v[4563] * v[5225] + v[4564] * v[5226] + v[4560] * v[5228] + v[4561] * v[5229] + v[4562] * v[5230]
			- v[4557] * v[5278] - v[4558] * v[5279] - v[4559] * v[5280]) + v[5516] - v[4568] * v[5714] - v[4570] * v[5715]
			- v[5064] * v[6921] + v[5236] * v[6974] + v[5235] * v[6975] + v[5234] * v[6976] + v[5232] * v[6978] + v[5231] * v[6979]
			+ v[389] * (v[4563] * v[5219] + v[4564] * v[5220] + v[4560] * v[5222] + v[4561] * v[5223] + v[4562] * v[5224]
				- v[4557] * v[5281] - v[4558] * v[5282] - v[4559] * v[5283] + v[7039]) - v[5277] * v[7040] - v[5276] * v[7041]
			- v[5275] * v[7042];
		v[5778] = v[5516];
		v[5723] = v[5718] + v[5724];
		v[5725] = v[5711] + v[5724];
		v[5514] = v[5514] + v[4535] * (v[5482] - v[391] * v[5715]);
		v[5780] = v[5514];
		v[5504] = v[5504] + v[4535] * (v[5484] - v[391] * v[5714]) - v[5062] * v[6921];
		v[5782] = v[5504];
		v[5726] = v[5517] + v[5810];
		v[5727] = v[4535] * v[5432];
		v[5728] = v[4535] * v[5431];
		v[5422] = v[5422] + v[5441];
		b5729 = b6;
		if (b5729) {
			v[5422] = 0e0;
			v[5431] = 0e0;
			v[5432] = 0e0;
			v[5433] = 0e0;
		}
		else {
			b5730 = b1049;
			if (b5730) {
				b5731 = b1051;
				if (b5731) {
					v[5422] = 0e0;
					v[5431] = 0e0;
					v[5432] = 0e0;
					v[5433] = 0e0;
				}
				else {
				};
			}
			else {
			};
		};
		v[5416] = v[5416] + v[5440];
		v[5736] = v[5416];
		v[5750] = v[5736];
		v[5420] = v[5420] + v[5439];
		v[5737] = v[5420];
		v[5749] = v[5737];
		v[5426] = v[5426] + v[5438];
		v[5738] = v[5426];
		v[5748] = v[5738];
		v[5739] = 0e0;
		v[5740] = 0e0;
		v[5741] = 0e0;
		v[5742] = 0e0;
		v[5743] = 0e0;
		v[5744] = 0e0;
		b5745 = b6;
		if (b5745) {
			b5746 = b1040;
			if (b5746) {
				v[5747] = v[1023] * v[5738];
				v[5741] = v[1037] * v[5738];
				v[5738] = 0e0;
				v[5747] = v[1021] * v[5737] + v[5747];
				v[5740] = v[1037] * v[5737];
				v[5737] = 0e0;
				v[5747] = v[1015] * v[5736] + v[5747];
				v[5739] = v[1037] * v[5736];
				v[5736] = 0e0;
			}
			else {
				v[5747] = 0e0;
				v[5744] = -v[5748];
				v[5738] = 0e0;
				v[5743] = -v[5749];
				v[5737] = 0e0;
				v[5742] = -v[5750];
				v[5736] = 0e0;
			};
			v[5751] = v[5744] * v[6731];
			v[5516] = v[5516] + v[5744] * v[6507];
			v[5744] = 0e0;
			v[5752] = v[5751] + v[5743] * v[6732];
			v[5514] = v[5514] + v[5743] * v[6507];
			v[5743] = 0e0;
			v[5753] = v[5752] + v[5742] * v[6733];
			v[5504] = v[5504] + v[5742] * v[6507];
			v[5742] = 0e0;
			v[5424] = v[5424] + v[29] * v[5753] * v[7476] + (v[5747] * v[5754] * v[7477]) / v[1785];
		}
		else {
			v[5755] = 0e0;
			b5756 = b1049;
			if (b5756) {
				v[5757] = 0e0;
				v[5758] = 0e0;
				v[5759] = 0e0;
				b5760 = b1064;
				if (b5760) {
					v[5759] = v[5748];
					v[5738] = 0e0;
					v[5758] = v[5749];
					v[5737] = 0e0;
					v[5757] = v[5750];
					v[5736] = 0e0;
				}
				else {
					v[5744] = -v[5748];
					v[5738] = 0e0;
					v[5743] = -v[5749];
					v[5737] = 0e0;
					v[5742] = -v[5750];
					v[5736] = 0e0;
				};
				v[5769] = v[1015] * v[5757] + v[1021] * v[5758] + v[1023] * v[5759];
				b5764 = b1051;
				if (b5764) {
					v[5741] = v[1059] * v[5759];
					v[5740] = v[1059] * v[5758];
					v[5739] = v[1059] * v[5757];
					v[5755] = (v[5754] * v[5769] * Power(v[1057], v[3585])) / v[1799];
				}
				else {
					v[5741] = v[1063] * v[5759];
					v[5740] = v[1063] * v[5758];
					v[5739] = v[1063] * v[5757];
					v[5424] = v[5774] - (v[3801] * v[5769] * v[6726] * Power(v[401], v[3606])) / v[1805];
				};
			}
			else {
			};
			v[5786] = v[5424];
			v[5785] = v[5742];
			v[5784] = v[5743];
			v[5783] = v[5744];
			b5775 = b1049;
			if (b5775) {
				b5776 = b1051;
				if (b5776) {
					v[5777] = -(v[5744] * v[6731]);
					v[5516] = v[5778] + v[5744] * v[6508];
					v[5744] = 0e0;
					v[5779] = v[5777] - v[5743] * v[6732];
					v[5514] = v[5780] + v[5743] * v[6508];
					v[5743] = 0e0;
					v[5781] = v[5779] - v[5742] * v[6733];
					v[5504] = v[5782] + v[5742] * v[6508];
					v[5742] = 0e0;
					v[5755] = v[5755] + v[29] * v[5781] * Power(v[1057], v[1045]);
					v[5424] = v[5424] - v[5755];
				}
				else {
					v[5516] = v[5778] + v[1054] * v[5783];
					v[5744] = 0e0;
					v[5514] = v[5780] + v[1054] * v[5784];
					v[5743] = 0e0;
					v[5504] = v[5782] + v[1054] * v[5785];
					v[5742] = 0e0;
					v[5424] = v[5786] - (v[391] * v[5783] + v[390] * v[5784] + v[389] * v[5785])*v[6509] * Power(v[401], v[1062]);
				};
			}
			else {
			};
		};
		v[5930] = v[5424];
		v[7071] = v[1014] * v[5739];
		v[7045] = -(v[389] * v[5739]);
		v[7070] = v[1018] * v[5740];
		v[7044] = v[390] * v[5740];
		v[7043] = v[1022] * v[5741];
		v[5787] = v[4535] * v[5708] + v[372] * v[5741];
		v[7083] = v[391] * v[5787];
		v[5788] = v[4535] * v[5707] + v[362] * v[5741];
		v[7080] = v[391] * v[5788];
		v[5789] = v[4535] * v[5706] + v[352] * v[5741];
		v[7079] = v[391] * v[5789];
		v[5790] = v[4535] * v[5705] + v[346] * v[5741];
		v[7078] = v[391] * v[5790];
		v[5791] = v[4535] * v[5704] + v[336] * v[5741];
		v[7077] = v[391] * v[5791];
		v[5792] = v[4535] * v[5703] + v[326] * v[5741];
		v[7076] = v[391] * v[5792];
		v[5793] = v[4535] * v[5702] + v[320] * v[5741];
		v[7075] = v[391] * v[5793];
		v[5794] = v[4535] * v[5701] + v[310] * v[5741];
		v[7074] = v[391] * v[5794];
		v[5795] = v[4535] * v[5700] + v[300] * v[5741];
		v[7073] = v[391] * v[5795];
		v[5504] = v[5504] + v[3825] * v[5741];
		v[5514] = v[5514] + v[3826] * v[5741];
		v[5797] = v[389] * v[5741];
		v[5799] = v[390] * v[5741];
		v[5516] = v[5516] + v[4572] * v[5741];
		v[5823] = v[4534] * v[5708] + v[372] * v[5740];
		v[7100] = v[390] * v[5823];
		v[5824] = v[4534] * v[5707] + v[362] * v[5740];
		v[7098] = v[390] * v[5824];
		v[5825] = v[4534] * v[5706] + v[352] * v[5740];
		v[7096] = v[390] * v[5825];
		v[5826] = v[4534] * v[5705] + v[346] * v[5740];
		v[7094] = v[390] * v[5826];
		v[5827] = v[4534] * v[5704] + v[336] * v[5740];
		v[7092] = v[390] * v[5827];
		v[5828] = v[4534] * v[5703] + v[326] * v[5740];
		v[7090] = v[390] * v[5828];
		v[5829] = v[4534] * v[5702] + v[320] * v[5740];
		v[7088] = v[390] * v[5829];
		v[5830] = v[4534] * v[5701] + v[310] * v[5740];
		v[7086] = v[390] * v[5830];
		v[5831] = v[4534] * v[5700] + v[300] * v[5740];
		v[7084] = v[390] * v[5831];
		v[5504] = v[5504] - v[292] * v[7044];
		v[5514] = v[5514] + v[4602] * v[5740];
		v[5516] = v[5516] - v[294] * v[7044];
		v[5857] = v[4533] * v[5708] + v[372] * v[5739];
		v[7101] = v[389] * v[5857];
		v[5858] = v[4533] * v[5707] + v[362] * v[5739];
		v[7099] = v[389] * v[5858];
		v[5859] = v[4533] * v[5706] + v[352] * v[5739];
		v[7097] = v[389] * v[5859];
		v[5860] = v[4533] * v[5705] + v[346] * v[5739];
		v[7095] = v[389] * v[5860];
		v[5861] = v[4533] * v[5704] + v[336] * v[5739];
		v[7093] = v[389] * v[5861];
		v[5862] = v[4533] * v[5703] + v[326] * v[5739];
		v[7091] = v[389] * v[5862];
		v[5863] = v[4533] * v[5702] + v[320] * v[5739];
		v[7089] = v[389] * v[5863];
		v[5864] = v[4533] * v[5701] + v[310] * v[5739];
		v[7087] = v[389] * v[5864];
		v[5865] = v[4533] * v[5700] + v[300] * v[5739];
		v[7085] = v[389] * v[5865];
		v[5867] = v[4534] * v[5058] + v[4533] * v[5059] + v[287] * v[5739] + v[286] * v[5740];
		v[5868] = v[4534] * v[5060] + v[4533] * v[5061] + v[290] * v[5739] + v[289] * v[5740];
		v[7082] = v[223] * v[5867] + v[224] * v[5868];
		v[5504] = v[5504] + v[4636] * v[5739];
		v[5514] = v[5514] + v[293] * v[7045];
		v[5516] = v[5516] + v[294] * v[7045];
		v[5871] = v[1968] * v[5739] + v[1889] * v[5740] + v[1816] * v[5741] + v[4533] * (-(v[1014] * v[5281]) + v[5622]
			+ v[389] * (v[7046] + v[7047])) + v[4534] * (-(v[1018] * v[5278]) + v[5553] + v[390] * (v[7046] + v[7048])
				- v[7016] * v[783]) + v[4535] * (-(v[1022] * v[5275]) + v[5494] + v[391] * (v[7047] + v[7048]) - v[7020] * v[783]
					- v[6981] * v[786]);
		v[6095] = v[5871] * v[6453];
		v[5872] = v[1970] * v[5739] + v[1891] * v[5740] + v[1818] * v[5741] + v[4533] * (-(v[1014] * v[5282]) + v[5621]
			+ v[389] * (v[7049] + v[7050])) + v[4534] * (-(v[1018] * v[5279]) + v[5552] + v[390] * (v[7049] + v[7051])
				- v[7016] * v[782]) + v[4535] * (-(v[1022] * v[5276]) + v[5493] + v[391] * (v[7050] + v[7051]) - v[7020] * v[782]
					- v[6981] * v[785]);
		v[7110] = v[350] * v[5872];
		v[6099] = v[5872] * v[6452];
		v[5873] = v[1972] * v[5739] + v[1893] * v[5740] + v[1820] * v[5741] + v[4533] * (-(v[1014] * v[5283]) + v[5620]
			+ v[389] * (v[7052] + v[7053])) + v[4534] * (-(v[1018] * v[5280]) + v[5551] + v[390] * (v[7052] + v[7054])
				- v[7016] * v[781]) + v[4535] * (-(v[1022] * v[5277]) + v[5492] + v[391] * (v[7053] + v[7054]) - v[7020] * v[781]
					- v[6981] * v[784]);
		v[7111] = v[350] * v[5873];
		v[6102] = v[5873] * v[6451];
		v[5874] = v[1974] * v[5739] + v[1895] * v[5740] + v[1822] * v[5741] + v[4533] * (v[1014] * v[5222] + v[5619] + v[389] *
			(v[7055] + v[7056])) + v[4535] * (v[1022] * v[5234] + v[5491] + v[684] * v[6981] + v[678] * v[7020] + v[391] * (v[7055]
				+ v[7057])) + v[4534] * (v[1018] * v[5228] + v[5550] + v[678] * v[7016] + v[390] * (v[7056] + v[7057]));
		v[6222] = v[5874] * v[6446];
		v[5875] = v[1976] * v[5739] + v[1897] * v[5740] + v[1824] * v[5741] + v[4533] * (v[1014] * v[5223] + v[5618] + v[389] *
			(v[7058] + v[7059])) + v[4535] * (v[1022] * v[5235] + v[5490] + v[683] * v[6981] + v[677] * v[7020] + v[391] * (v[7058]
				+ v[7060])) + v[4534] * (v[1018] * v[5229] + v[5549] + v[677] * v[7016] + v[390] * (v[7059] + v[7060]));
		v[7119] = v[324] * v[5875];
		v[6226] = v[5875] * v[6445];
		v[5876] = v[1978] * v[5739] + v[1899] * v[5740] + v[1826] * v[5741] + v[4533] * (v[1014] * v[5224] + v[5617] + v[389] *
			(v[7061] + v[7062])) + v[4535] * (v[1022] * v[5236] + v[5489] + v[682] * v[6981] + v[676] * v[7020] + v[391] * (v[7061]
				+ v[7063])) + v[4534] * (v[1018] * v[5230] + v[5548] + v[676] * v[7016] + v[390] * (v[7062] + v[7063]));
		v[7120] = v[324] * v[5876];
		v[6229] = v[5876] * v[6444];
		v[5877] = v[1980] * v[5739] + v[1901] * v[5740] + v[1828] * v[5741] + v[4533] * (v[1014] * v[5219] + v[5616] + v[389] *
			(v[7064] + v[7065])) + v[4535] * (v[1022] * v[5231] + v[5488] + v[681] * v[6981] + v[675] * v[7020] + v[391] * (v[7064]
				+ v[7066])) + v[4534] * (v[1018] * v[5225] + v[5547] + v[675] * v[7016] + v[390] * (v[7065] + v[7066]));
		v[6282] = v[5877] * v[6439];
		v[5878] = v[1982] * v[5739] + v[1903] * v[5740] + v[1830] * v[5741] + v[4533] * (v[1014] * v[5220] + v[5615] + v[389] *
			(v[7067] + v[7068])) + v[4535] * (v[1022] * v[5232] + v[5487] + v[680] * v[6981] + v[674] * v[7020] + v[391] * (v[7067]
				+ v[7069])) + v[4534] * (v[1018] * v[5226] + v[5546] + v[674] * v[7016] + v[390] * (v[7068] + v[7069]));
		v[7128] = v[298] * v[5878];
		v[6286] = v[5878] * v[6438];
		v[5879] = v[4533] * v[5614] + v[1984] * v[5739] + v[1905] * v[5740] + v[1832] * v[5741] + v[4534] * (v[1018] * v[5227]
			+ v[5545] + v[390] * (v[7017] + v[7018])) + v[4535] * (v[1022] * v[5233] + v[5486] + v[679] * v[6981] + v[391] * (v[7018]
				+ v[7019]));
		v[7129] = v[298] * v[5879];
		v[6289] = v[5879] * v[6437];
		v[5889] = v[5681] - v[1014] * v[5857] - v[389] * (v[5571] + v[7083] + v[7100]);
		v[5890] = v[5678] - v[1014] * v[5858] - v[389] * (v[5570] + v[7080] + v[7098]);
		v[5891] = v[5675] - v[1014] * v[5859] - v[389] * (v[5569] + v[7079] + v[7096]);
		v[5892] = v[5690] + v[1014] * v[5860] + v[389] * (v[5568] + v[7078] + v[7094]);
		v[5893] = v[5687] + v[1014] * v[5861] + v[389] * (v[5567] + v[7077] + v[7092]);
		v[5894] = v[5684] + v[1014] * v[5862] + v[389] * (v[5566] + v[7076] + v[7090]);
		v[5895] = v[5697] + v[1014] * v[5863] + v[389] * (v[5565] + v[7075] + v[7088]);
		v[5896] = v[5694] + v[1014] * v[5864] + v[389] * (v[5564] + v[7074] + v[7086]);
		v[5504] = v[5504] - v[292] * v[5723] + v[3283] * v[5857] + v[3285] * v[5858] + v[3287] * v[5859] + v[3289] * v[5860]
			+ v[3291] * v[5861] + v[3293] * v[5862] + v[3295] * v[5863] + v[3297] * v[5864] + v[3299] * v[5865] + v[5563] * v[673]
			+ v[5564] * v[674] + v[5565] * v[675] + v[5566] * v[676] + v[5567] * v[677] + v[5568] * v[678] - v[5569] * v[781]
			- v[5570] * v[782] - v[5571] * v[783] + v[391] * (v[5795] * v[673] + v[5794] * v[674] + v[5793] * v[675] + v[5792] * v[676]
				+ v[5791] * v[677] + v[5790] * v[678] - v[5789] * v[781] - v[5788] * v[782] - v[5787] * v[783]) + v[390] *
				(v[5831] * v[673] + v[5830] * v[674] + v[5829] * v[675] + v[5828] * v[676] + v[5827] * v[677] + v[5826] * v[678] + v[7082]
					- v[5825] * v[781] - v[5824] * v[782] - v[5823] * v[783]) + v[6686] * (v[4629] * v[5219] + v[4630] * v[5220]
						+ v[4631] * v[5221] + v[4626] * v[5222] + v[4627] * v[5223] + v[4628] * v[5224] - v[4623] * v[5281] - v[4624] * v[5282]
						- v[4625] * v[5283] + v[4632] * v[5739] + v[5865] * v[673] + v[5864] * v[674] + v[5863] * v[675] + v[5862] * v[676]
						+ v[5861] * v[677] + v[5860] * v[678] + v[4533] * (v[223] * v[5058] + v[224] * v[5060] - v[5062] - v[7072])
						- v[5859] * v[781] - v[5858] * v[782] - v[5857] * v[783]);
		v[5925] = v[5504];
		v[5897] = v[5691] + v[1014] * v[5865] + v[389] * (v[5563] + v[7073] + v[7084]);
		v[5898] = v[5604] + v[1018] * v[5831] + v[390] * (v[5644] + v[7073] + v[7085]);
		v[5899] = v[5606] + v[1018] * v[5830] + v[390] * (v[5645] + v[7074] + v[7087]);
		v[5900] = v[5608] + v[1018] * v[5829] + v[390] * (v[5646] + v[7075] + v[7089]);
		v[5901] = v[5598] + v[1018] * v[5828] + v[390] * (v[5647] + v[7076] + v[7091]);
		v[5902] = v[5600] + v[1018] * v[5827] + v[390] * (v[5648] + v[7077] + v[7093]);
		v[5903] = v[5602] + v[1018] * v[5826] + v[390] * (v[5649] + v[7078] + v[7095]);
		v[5904] = v[5592] - v[1018] * v[5825] - v[390] * (v[5650] + v[7079] + v[7097]);
		v[5905] = v[5594] - v[1018] * v[5824] - v[390] * (v[5651] + v[7080] + v[7099]);
		v[5514] = v[5514] - v[293] * v[5725] + v[5823] * v[6482] + v[5824] * v[6484] + v[5825] * v[6486] + v[5826] * v[6488]
			+ v[5827] * v[6490] + v[5828] * v[6492] + v[5829] * v[6494] + v[5830] * v[6496] + v[5831] * v[6498] + v[5644] * v[679]
			+ v[5645] * v[680] + v[5646] * v[681] + v[5647] * v[682] + v[5648] * v[683] + v[5649] * v[684] - v[5650] * v[784]
			- v[5651] * v[785] - v[5652] * v[786] + v[391] * (v[5795] * v[679] + v[5794] * v[680] + v[5793] * v[681] + v[5792] * v[682]
				+ v[5791] * v[683] + v[5790] * v[684] - v[5789] * v[784] - v[5788] * v[785] - v[5787] * v[786]) + v[6683] *
				(v[4597] * v[5225] + v[4598] * v[5226] + v[4599] * v[5227] + v[4594] * v[5228] + v[4595] * v[5229] + v[4596] * v[5230]
					- v[4591] * v[5278] - v[4592] * v[5279] - v[4593] * v[5280] + v[4600] * v[5740] + v[5831] * v[679] + v[5830] * v[680]
					+ v[5829] * v[681] + v[5828] * v[682] + v[5827] * v[683] + v[5826] * v[684] + v[4534] * (v[223] * v[5059]
						+ v[224] * v[5061] - v[5063] - v[7081]) - v[5825] * v[784] - v[5824] * v[785] - v[5823] * v[786]) + v[389] *
						(v[5865] * v[679] + v[5864] * v[680] + v[5863] * v[681] + v[5862] * v[682] + v[5861] * v[683] + v[5860] * v[684] + v[7082]
							- v[5859] * v[784] - v[5858] * v[785] - v[5857] * v[786]);
		v[5923] = v[5514];
		v[5906] = v[5596] - v[1018] * v[5823] - v[390] * (v[5652] + v[7083] + v[7101]);
		v[5907] = v[5719] + v[7044] - v[7045];
		v[7155] = v[391] * v[5907];
		v[7206] = v[5726] + v[7155];
		v[5908] = v[5728] + v[5797];
		v[7156] = v[391] * v[5908];
		v[7204] = v[5713] + v[7156];
		v[5909] = v[5727] + v[5799];
		v[7157] = v[391] * v[5909];
		v[7205] = v[5720] + v[7157];
		v[13299] = v[4534] * v[5612] + v[1016] * v[5740] - v[5852] - v[5853] + v[223] * (v[7071] + v[7204]);
		v[13300] = v[4533] * v[5612] + v[1016] * v[5739] - v[5849] - v[5850] + v[223] * (v[7070] + v[7205]);
		v[13301] = -v[5816] - v[5817] + v[223] * (v[7043] + v[7206]);
		v[13302] = v[4646] * v[5286] + v[4645] * v[5289] + v[4497] * v[5879] + v[5172] * (v[4488] * v[4645] + v[4490] * v[4646]
			+ v[6292]) + v[5878] * v[6457] + v[5877] * v[6460] + v[5151] * v[6946] + v[5157] * v[6947] + v[5319] * v[7124]
			+ v[5089] * v[7207] + v[166] * v[7208] + v[171] * v[7209];
		v[13303] = v[4647] * v[5286] + v[4645] * v[5305] + v[2503] * v[5877] + v[4500] * v[5878] + v[2504] * v[5879] + v[5172] *
			(v[4487] * v[4645] + v[4491] * v[4647] + v[6285]) + v[5166] * v[6948] + v[5323] * v[7123] + v[5087] * v[7210]
			+ v[5148] * v[7211] + v[163] * v[7212] + v[5313] * v[7213];
		v[13304] = v[4647] * v[5289] + v[4646] * v[5305] + v[4502] * v[5877] + v[2505] * v[5878] + v[2506] * v[5879] + v[5172] *
			(v[4489] * v[4646] + v[4492] * v[4647] + v[6277]) + v[5310] * v[7127] + v[5085] * v[7214] + v[5153] * v[7215]
			+ v[5159] * v[7216] + v[5317] * v[7217] + v[5321] * v[7218];
		v[13305] = v[4534] * v[5611] + v[1017] * v[5740] + v[5852] + v[5853] + v[224] * (v[7071] + v[7204]);
		v[13306] = v[4533] * v[5611] + v[1017] * v[5739] + v[5849] + v[5850] + v[224] * (v[7070] + v[7205]);
		v[13307] = v[5816] + v[5817] + v[224] * (v[7043] + v[7206]);
		v[13308] = v[4643] * v[5326] + v[4642] * v[5329] + v[4504] * v[5876] + v[5210] * (v[4461] * v[4642] + v[4463] * v[4643]
			+ v[6232]) + v[5875] * v[6465] + v[5874] * v[6468] + v[5189] * v[6940] + v[5195] * v[6941] + v[5359] * v[7115]
			+ v[5115] * v[7219] + v[185] * v[7220] + v[190] * v[7221];
		v[13309] = v[4644] * v[5326] + v[4642] * v[5345] + v[2511] * v[5874] + v[4507] * v[5875] + v[2512] * v[5876] + v[5210] *
			(v[4460] * v[4642] + v[4464] * v[4644] + v[6225]) + v[5204] * v[6942] + v[5363] * v[7114] + v[5113] * v[7222]
			+ v[5186] * v[7223] + v[182] * v[7224] + v[5353] * v[7225];
		v[13310] = v[4644] * v[5329] + v[4643] * v[5345] + v[4509] * v[5874] + v[2513] * v[5875] + v[2514] * v[5876] + v[5210] *
			(v[4462] * v[4643] + v[4465] * v[4644] + v[6217]) + v[5350] * v[7118] + v[5111] * v[7226] + v[5191] * v[7227]
			+ v[5197] * v[7228] + v[5357] * v[7229] + v[5361] * v[7230];
		v[13311] = v[5722] + v[391] * (-v[5728] - v[5797]) + v[389] * (-v[5723] - v[7044]) - v[7071];
		v[13312] = -(v[390] * v[5725]) + v[391] * (-v[5727] - v[5799]) - v[5812] - v[5739] * v[6499] - v[4533] * v[6980]
			- v[7070];
		v[13313] = -v[5517] - v[5810] - v[7043] + v[391] * (-v[5719] - v[7044] + v[7045]);
		v[13314] = v[4640] * v[5366] + v[4639] * v[5369] + v[4511] * v[5873] + v[5266] * (v[4433] * v[4639] + v[4435] * v[4640]
			+ v[6105]) + v[5872] * v[6473] + v[5871] * v[6476] + v[5245] * v[6934] + v[5251] * v[6935] + v[5399] * v[7106]
			+ v[5141] * v[7231] + v[245] * v[7232] + v[250] * v[7233];
		v[13315] = v[4641] * v[5366] + v[4639] * v[5385] + v[2519] * v[5871] + v[4514] * v[5872] + v[2520] * v[5873] + v[5266] *
			(v[4432] * v[4639] + v[4436] * v[4641] + v[6098]) + v[5260] * v[6936] + v[5403] * v[7105] + v[5139] * v[7234]
			+ v[5242] * v[7235] + v[242] * v[7236] + v[5393] * v[7237];
		v[13316] = v[4641] * v[5369] + v[4640] * v[5385] + v[4516] * v[5871] + v[2521] * v[5872] + v[2522] * v[5873] + v[5266] *
			(v[4434] * v[4640] + v[4437] * v[4641] + v[6090]) + v[5390] * v[7109] + v[5137] * v[7238] + v[5247] * v[7239]
			+ v[5253] * v[7240] + v[5397] * v[7241] + v[5401] * v[7242];
		v[5910] = v[5529] + v[1022] * v[5795] + v[391] * (v[5656] + v[7084] + v[7085]);
		v[5911] = v[5530] + v[1022] * v[5794] + v[391] * (v[5658] + v[7086] + v[7087]);
		v[5912] = v[5531] + v[1022] * v[5793] + v[391] * (v[5660] + v[7088] + v[7089]);
		v[5913] = v[5526] + v[1022] * v[5792] + v[391] * (v[5662] + v[7090] + v[7091]);
		v[5914] = v[5527] + v[1022] * v[5791] + v[391] * (v[5664] + v[7092] + v[7093]);
		v[5915] = v[5528] + v[1022] * v[5790] + v[391] * (v[5666] + v[7094] + v[7095]);
		v[5916] = v[5523] - v[1022] * v[5789] - v[391] * (v[5668] + v[7096] + v[7097]);
		v[5917] = v[5524] - v[1022] * v[5788] - v[391] * (v[5670] + v[7098] + v[7099]);
		v[5516] = v[5516] - v[294] * v[5719] - v[293] * v[5727] - v[292] * v[5728] + v[3482] * v[5787] + v[3484] * v[5788]
			+ v[3486] * v[5789] + v[3488] * v[5790] + v[3490] * v[5791] + v[3492] * v[5792] + v[3494] * v[5793] + v[3496] * v[5794]
			+ v[3498] * v[5795] + v[5461] * v[5907] + v[5483] * v[5908] + v[5481] * v[5909] + v[5656] * v[685] + v[5658] * v[686]
			+ v[5660] * v[687] + v[5662] * v[688] + v[5664] * v[689] + v[5666] * v[690] - v[5668] * v[787] - v[5670] * v[788]
			- v[5672] * v[789] + v[6680] * (v[4563] * v[5231] + v[4564] * v[5232] + v[4565] * v[5233] + v[4560] * v[5234]
				+ v[4561] * v[5235] + v[4562] * v[5236] - v[4557] * v[5275] - v[4558] * v[5276] - v[4559] * v[5277] + v[4566] * v[5741]
				+ v[5795] * v[685] + v[5794] * v[686] + v[5793] * v[687] + v[5792] * v[688] + v[5791] * v[689] + v[5790] * v[690]
				- v[4535] * (v[5064] - v[6973] + v[12452 + i4521] * v[7]) - v[5789] * v[787] - v[5788] * v[788] - v[5787] * v[789])
			+ v[390] * (v[5831] * v[685] + v[5830] * v[686] + v[5829] * v[687] + v[5828] * v[688] + v[5827] * v[689] + v[5826] * v[690]
				- v[5825] * v[787] - v[5824] * v[788] - v[5823] * v[789]) + v[389] * (v[5865] * v[685] + v[5864] * v[686]
					+ v[5863] * v[687] + v[5862] * v[688] + v[5861] * v[689] + v[5860] * v[690] - v[5859] * v[787] - v[5858] * v[788]
					- v[5857] * v[789]);
		v[5921] = v[5516];
		v[5918] = v[5525] - v[1022] * v[5787] - v[391] * (v[5672] + v[7100] + v[7101]);
		b5919 = b6;
		if (b5919) {
			v[5428] = v[5428] + v[375] * v[5516];
			v[5425] = v[5425] + v[388] * v[5516];
			v[5516] = 0e0;
			v[5428] = v[5428] + v[374] * v[5514];
			v[5418] = v[5418] + v[388] * v[5514];
			v[5514] = 0e0;
			v[5428] = v[5428] + v[373] * v[5504];
			v[5414] = v[5414] + v[388] * v[5504];
			v[5504] = 0e0;
			v[5427] = v[5427] + v[387] * v[5428];
			v[5424] = v[5424] + v[5427];
		}
		else {
			b5920 = b418;
			if (b5920) {
				v[5922] = -v[5921];
				v[5516] = 0e0;
				v[5924] = -v[5923];
				v[5514] = 0e0;
				v[5926] = -v[5925];
				v[5504] = 0e0;
			}
			else {
				v[5922] = v[5921];
				v[5924] = v[5923];
				v[5926] = v[5925];
			};
			v[5430] = v[5430] + v[375] * v[5922];
			v[5425] = v[407] * v[5922] + v[5927];
			v[5430] = v[5430] + v[374] * v[5924];
			v[5418] = v[407] * v[5924] + v[5928];
			v[5430] = v[5430] + v[373] * v[5926];
			v[5414] = v[407] * v[5926] + v[5929];
			v[5429] = v[5429] + v[406] * v[5430];
			v[5424] = v[5429] + v[5930];
		};
		v[7102] = v[5424] / v[401];
		v[5425] = v[5425] + v[375] * v[7102];
		v[5418] = v[5418] + v[374] * v[7102];
		v[5414] = v[5414] + v[373] * v[7102];
		v[5935] = -(v[4523] * v[5045]) - v[1027] * v[5425];
		v[5936] = -(v[4523] * v[5053]) - v[274] * v[5425];
		v[5937] = -(v[4523] * v[5044]) - v[1613] * v[5425];
		v[5938] = v[223] * v[5425] + v[5939];
		v[5940] = v[224] * v[5425] - v[5939];
		v[5947] = -(v[4524] * v[5045]) - v[1027] * v[5418];
		v[5948] = -(v[4524] * v[5053]) - v[274] * v[5418];
		v[5949] = -(v[4524] * v[5044]) - v[1613] * v[5418];
		v[5950] = v[223] * v[5418] + v[5951];
		v[5952] = v[224] * v[5418] - v[5951];
		v[5959] = -(v[4525] * v[5045]) - v[1027] * v[5414];
		v[5960] = -(v[4525] * v[5053]) - v[274] * v[5414];
		v[5961] = -(v[4525] * v[5044]) - v[1613] * v[5414];
		v[5962] = v[223] * v[5414] + v[5963];
		v[5964] = v[224] * v[5414] - v[5963];
		v[5970] = v[359] * v[5871] + v[4639] * v[7103];
		v[5971] = v[353] * v[5871] + v[4639] * v[7104];
		v[5973] = v[2173] * v[5871] + v[4639] * (v[5379] / v[350] + v[5266] * v[5972]) + v[5970] * v[7105];
		v[5975] = v[2179] * v[5871] + v[4639] * (v[5374] / v[350] + v[5266] * v[5974]) + v[5971] * v[7106];
		v[5977] = (v[252] * v[5970] + v[250] * v[5971] + v[4639] * (v[5384] + v[5976] * v[6966]) + v[5871] * v[7107]) / v[350];
		v[5978] = v[361] * v[5872] + v[4640] * v[7108];
		v[5979] = v[353] * v[5872] + v[4640] * v[7104];
		v[5981] = v[2168] * v[5872] + v[4640] * (v[5382] / v[350] + v[5266] * v[5980]) + v[5978] * v[7109];
		v[5982] = v[5970] + v[5978];
		v[5983] = v[5973] + v[5981];
		v[5985] = (v[4763] * v[5253] + v[241] * v[5979] + v[243] * v[5982] + v[6358] * v[6966] + v[4640] * (v[5370]
			+ v[5984] * v[6966]) + v[2178] * v[7110]) / v[350];
		v[5987] = (v[248] * v[5978] + v[245] * v[5979] + v[4640] * (v[5376] + v[5986] * v[6966]) + v[2172] * v[7110]) / v[350];
		v[5988] = v[359] * v[5873] + v[4641] * v[7103];
		v[5989] = v[361] * v[5873] + v[4641] * v[7108];
		v[5990] = v[5979] + v[5988];
		v[5992] = (v[242] * v[5988] + v[243] * v[5989] + v[4641] * (v[5371] + v[5991] * v[6966]) + v[2176] * v[7111]) / v[350];
		v[5993] = v[5971] + v[5989];
		v[5995] = (v[4768] * v[5247] + v[247] * v[5988] + v[248] * v[5993] + v[6357] * v[6966] + v[4641] * (v[5375]
			+ v[5994] * v[6966]) + v[2183] * v[7111]) / v[350];
		v[5996] = v[5985] + v[5995];
		v[5998] = (v[4769] * v[5242] + v[253] * v[5989] + v[252] * v[5990] + v[6356] * v[6966] + v[4641] * (v[5380]
			+ v[5997] * v[6966]) + v[2186] * v[7111]) / v[350];
		v[5999] = v[5975] + v[5998];
		v[6000] = v[333] * v[5874] + v[4642] * v[7112];
		v[6001] = v[327] * v[5874] + v[4642] * v[7113];
		v[6003] = v[2197] * v[5874] + v[4642] * (v[5339] / v[324] + v[5210] * v[6002]) + v[6000] * v[7114];
		v[6005] = v[2203] * v[5874] + v[4642] * (v[5334] / v[324] + v[5210] * v[6004]) + v[6001] * v[7115];
		v[6007] = (v[192] * v[6000] + v[190] * v[6001] + v[4642] * (v[5344] + v[6006] * v[6965]) + v[5874] * v[7116]) / v[324];
		v[6008] = v[335] * v[5875] + v[4643] * v[7117];
		v[6009] = v[327] * v[5875] + v[4643] * v[7113];
		v[6011] = v[2192] * v[5875] + v[4643] * (v[5342] / v[324] + v[5210] * v[6010]) + v[6008] * v[7118];
		v[6012] = v[6000] + v[6008];
		v[6013] = v[6003] + v[6011];
		v[6015] = (v[4796] * v[5197] + v[181] * v[6009] + v[183] * v[6012] + v[6384] * v[6965] + v[4643] * (v[5330]
			+ v[6014] * v[6965]) + v[2202] * v[7119]) / v[324];
		v[6017] = (v[188] * v[6008] + v[185] * v[6009] + v[4643] * (v[5336] + v[6016] * v[6965]) + v[2196] * v[7119]) / v[324];
		v[6018] = v[333] * v[5876] + v[4644] * v[7112];
		v[6019] = v[335] * v[5876] + v[4644] * v[7117];
		v[6020] = v[6009] + v[6018];
		v[6022] = (v[182] * v[6018] + v[183] * v[6019] + v[4644] * (v[5331] + v[6021] * v[6965]) + v[2200] * v[7120]) / v[324];
		v[6023] = v[6001] + v[6019];
		v[6025] = (v[4801] * v[5191] + v[187] * v[6018] + v[188] * v[6023] + v[6383] * v[6965] + v[4644] * (v[5335]
			+ v[6024] * v[6965]) + v[2207] * v[7120]) / v[324];
		v[6026] = v[6015] + v[6025];
		v[6028] = (v[4802] * v[5186] + v[193] * v[6019] + v[192] * v[6020] + v[6382] * v[6965] + v[4644] * (v[5340]
			+ v[6027] * v[6965]) + v[2210] * v[7120]) / v[324];
		v[6029] = v[6005] + v[6028];
		v[6030] = v[307] * v[5877] + v[4645] * v[7121];
		v[6031] = v[301] * v[5877] + v[4645] * v[7122];
		v[6033] = v[2221] * v[5877] + v[4645] * (v[5299] / v[298] + v[5172] * v[6032]) + v[6030] * v[7123];
		v[6035] = v[2227] * v[5877] + v[4645] * (v[5294] / v[298] + v[5172] * v[6034]) + v[6031] * v[7124];
		v[6037] = (v[173] * v[6030] + v[171] * v[6031] + v[4645] * (v[5304] + v[6036] * v[6964]) + v[5877] * v[7125]) / v[298];
		v[6038] = v[309] * v[5878] + v[4646] * v[7126];
		v[6039] = v[301] * v[5878] + v[4646] * v[7122];
		v[6041] = v[2216] * v[5878] + v[4646] * (v[5302] / v[298] + v[5172] * v[6040]) + v[6038] * v[7127];
		v[6042] = v[6030] + v[6038];
		v[6043] = v[6033] + v[6041];
		v[6045] = (v[4829] * v[5159] + v[162] * v[6039] + v[164] * v[6042] + v[6407] * v[6964] + v[4646] * (v[5290]
			+ v[6044] * v[6964]) + v[2226] * v[7128]) / v[298];
		v[6047] = (v[169] * v[6038] + v[166] * v[6039] + v[4646] * (v[5296] + v[6046] * v[6964]) + v[2220] * v[7128]) / v[298];
		v[6048] = v[307] * v[5879] + v[4647] * v[7121];
		v[6049] = v[309] * v[5879] + v[4647] * v[7126];
		v[6050] = v[6039] + v[6048];
		v[6052] = (v[163] * v[6048] + v[164] * v[6049] + v[4647] * (v[5291] + v[6051] * v[6964]) + v[2224] * v[7129]) / v[298];
		v[6053] = v[6031] + v[6049];
		v[6055] = (v[4834] * v[5153] + v[168] * v[6048] + v[169] * v[6053] + v[6406] * v[6964] + v[4647] * (v[5295]
			+ v[6054] * v[6964]) + v[2231] * v[7129]) / v[298];
		v[6056] = v[6045] + v[6055];
		v[6058] = (v[4835] * v[5148] + v[174] * v[6049] + v[173] * v[6050] + v[6405] * v[6964] + v[4647] * (v[5300]
			+ v[6057] * v[6964]) + v[2234] * v[7129]) / v[298];
		v[6059] = v[6035] + v[6058];
		v[6060] = v[4663] * v[5044] + v[1613] * v[5891];
		v[6061] = v[4663] * v[5053] + v[274] * v[5891];
		v[6062] = v[4663] * v[5045] + v[1027] * v[5891];
		v[6063] = v[4660] * v[5044] + v[1613] * v[5890];
		v[6064] = v[4660] * v[5053] + v[274] * v[5890];
		v[6065] = v[4660] * v[5045] + v[1027] * v[5890];
		v[6066] = v[4657] * v[5044] + v[1613] * v[5889];
		v[6067] = v[4657] * v[5053] + v[274] * v[5889];
		v[6068] = v[4657] * v[5045] + v[1027] * v[5889];
		v[6069] = v[4664] * v[5044] + v[1613] * v[5904];
		v[6070] = v[4664] * v[5053] + v[274] * v[5904];
		v[6071] = v[4664] * v[5045] + v[1027] * v[5904];
		v[6072] = v[4661] * v[5044] + v[1613] * v[5905];
		v[6073] = v[4661] * v[5053] + v[274] * v[5905];
		v[6074] = v[4661] * v[5045] + v[1027] * v[5905];
		v[6075] = v[4658] * v[5044] + v[1613] * v[5906];
		v[6076] = v[4658] * v[5053] + v[274] * v[5906];
		v[6077] = v[4658] * v[5045] + v[1027] * v[5906];
		v[6078] = v[4665] * v[5044] + v[1613] * v[5916];
		v[6079] = v[4665] * v[5053] + v[274] * v[5916];
		v[6080] = v[4665] * v[5045] + v[1027] * v[5916];
		v[6081] = v[4662] * v[5044] + v[1613] * v[5917];
		v[6082] = v[4662] * v[5053] + v[274] * v[5917];
		v[6083] = v[4662] * v[5045] + v[1027] * v[5917];
		v[6084] = v[4659] * v[5044] + v[1613] * v[5918];
		v[6085] = v[4659] * v[5053] + v[274] * v[5918];
		v[6086] = v[4659] * v[5045] + v[1027] * v[5918];
		v[6087] = v[4885] * v[5043] - v[130] * v[5410] - v[133] * v[5411] - v[136] * v[5412] - v[259] * v[5414] - v[262] * v[5418]
			- v[265] * v[5425] + v[5891] * v[742] + v[5890] * v[743] + v[5889] * v[744] + v[5904] * v[757] + v[5905] * v[758]
			+ v[5906] * v[759] + v[5916] * v[772] + v[5917] * v[773] + v[5918] * v[774];
		v[6088] = -(v[4884] * v[5043]) - v[128] * v[5410] - v[131] * v[5411] - v[134] * v[5412] - v[257] * v[5414]
			- v[260] * v[5418] - v[263] * v[5425] + v[5891] * v[736] + v[5890] * v[737] + v[5889] * v[738] + v[5904] * v[751]
			+ v[5905] * v[752] + v[5906] * v[753] + v[5916] * v[766] + v[5917] * v[767] + v[5918] * v[768];
		v[6089] = v[136] * v[5935] + v[135] * v[5936] + v[134] * v[5937] + v[2519] * v[5978] + v[363] * v[6095]
			+ v[5989] * v[6476] + v[361] * (v[5390] / v[350] + v[5266] * v[7479]);
		v[6092] = v[133] * v[5935] + v[132] * v[5936] + v[131] * v[5937] + (v[4769] * v[5364] + v[359] * v[5393]
			+ v[363] * v[5970] + v[348] * v[5990]) / v[350] + v[360] * v[6099] + v[5266] * v[6353] + v[6090] * v[7130]
			+ v[6091] * v[7130];
		v[6094] = v[130] * v[5935] + v[129] * v[5936] + v[128] * v[5937] + v[4516] * v[5971] + v[348] * v[6102] + v[353] *
			(v[5266] * v[7131] + v[7233]);
		v[6096] = v[136] * v[5947] + v[135] * v[5948] + v[134] * v[5949] + (v[4768] * v[5367] + v[361] * v[5401]
			+ v[354] * v[5978] + v[349] * v[5993]) / v[350] + v[5266] * v[6354] + v[6095] * v[6477] + (v[6097] + v[6098])*v[7132];
		v[6100] = v[133] * v[5947] + v[132] * v[5948] + v[131] * v[5949] + v[2521] * v[5970] + v[354] * v[6099]
			+ v[5988] * v[6473] + v[359] * (v[5403] / v[350] + v[5266] * v[7480]);
		v[6103] = v[130] * v[5947] + v[129] * v[5948] + v[128] * v[5949] + v[4514] * v[5979] + v[349] * v[6102] + v[353] *
			(v[5266] * v[7133] + v[7232]);
		v[6104] = v[136] * v[5959] + v[135] * v[5960] + v[134] * v[5961] + v[5266] * v[6355] + (v[4763] * v[5368]
			+ v[361] * v[5397] + v[351] * v[5989] + v[5982] * v[6472]) / v[350] + v[6095] * v[6475] + (v[6105] + v[6107])*v[7132];
		v[6106] = v[133] * v[5959] + v[132] * v[5960] + v[131] * v[5961] + v[4511] * v[5988] + v[6099] * v[6472] + v[359] *
			(v[5266] * v[7134] + v[7236]);
		v[6109] = v[130] * v[5959] + v[129] * v[5960] + v[128] * v[5961] + v[2522] * v[5971] + v[2520] * v[5979]
			+ v[351] * v[6102] + v[353] * (v[5399] / v[350] + v[5266] * v[7481]);
		v[6110] = v[128] * v[6063] + v[129] * v[6064] + v[130] * v[6065];
		v[6111] = v[131] * v[6063] + v[132] * v[6064] + v[133] * v[6065];
		v[6112] = v[134] * v[6063] + v[135] * v[6064] + v[136] * v[6065];
		v[6113] = v[128] * v[6066] + v[129] * v[6067] + v[130] * v[6068];
		v[6114] = v[134] * v[6066] + v[135] * v[6067] + v[136] * v[6068];
		v[6115] = v[131] * v[6066] + v[132] * v[6067] + v[133] * v[6068];
		v[6116] = v[131] * v[6060] + v[132] * v[6061] + v[133] * v[6062];
		v[6117] = v[134] * v[6060] + v[135] * v[6061] + v[136] * v[6062];
		v[6118] = v[128] * v[6060] + v[129] * v[6061] + v[130] * v[6062];
		v[6119] = v[128] * v[6069] + v[129] * v[6070] + v[130] * v[6071];
		v[6120] = v[131] * v[6069] + v[132] * v[6070] + v[133] * v[6071];
		v[6121] = v[134] * v[6069] + v[135] * v[6070] + v[136] * v[6071];
		v[6122] = v[134] * v[6075] + v[135] * v[6076] + v[136] * v[6077];
		v[6123] = v[128] * v[6075] + v[129] * v[6076] + v[130] * v[6077];
		v[6124] = v[131] * v[6075] + v[132] * v[6076] + v[133] * v[6077];
		v[6125] = v[131] * v[6084] + v[132] * v[6085] + v[133] * v[6086];
		v[6126] = v[128] * v[6084] + v[129] * v[6085] + v[130] * v[6086];
		v[6127] = v[134] * v[6084] + v[135] * v[6085] + v[136] * v[6086];
		v[6128] = -(v[4765] * v[5134]) - v[4787] * v[5135] + 2e0*v[4764] * v[5136] + v[6116] + v[6119] + v[6122] + v[6125]
			+ v[5987] * v[6669] - v[5996] * v[701] - v[5983] * v[704];
		v[6129] = v[128] * v[6072] + v[129] * v[6073] + v[130] * v[6074];
		v[6130] = v[134] * v[6072] + v[135] * v[6073] + v[136] * v[6074];
		v[6131] = v[131] * v[6072] + v[132] * v[6073] + v[133] * v[6074];
		v[6132] = v[4791] * v[5134] + 2e0*v[4767] * v[5135] - v[4787] * v[5136] - v[6111] - v[6114] - v[6126] - v[6129]
			+ v[5992] * v[6668] - v[5996] * v[698] + v[5999] * v[704];
		v[6133] = -(v[4866] * v[5126]) + v[4862] * v[5128] - v[4858] * v[5129] + v[4780] * v[5243] + v[238] * v[6106]
			- v[6116] * v[691] + v[6111] * v[692] - v[6115] * v[693];
		v[6134] = v[128] * v[6078] + v[129] * v[6079] + v[130] * v[6080];
		v[6135] = v[131] * v[6078] + v[132] * v[6079] + v[133] * v[6080];
		v[6136] = v[134] * v[6078] + v[135] * v[6079] + v[136] * v[6080];
		v[6137] = v[131] * v[6081] + v[132] * v[6082] + v[133] * v[6083];
		v[6138] = v[128] * v[6081] + v[129] * v[6082] + v[130] * v[6083];
		v[6139] = v[134] * v[6081] + v[135] * v[6082] + v[136] * v[6083];
		v[6140] = 2e0*v[4761] * v[5134] + v[4791] * v[5135] - v[4765] * v[5136] - v[6117] - v[6130] - v[6134] - v[6137]
			+ v[5977] * v[6667] - v[5983] * v[698] + v[5999] * v[701];
		v[6141] = -(v[4865] * v[5126]) + v[4863] * v[5128] - v[4859] * v[5129] + v[4778] * v[5243] + v[238] * v[6104]
			- v[6117] * v[691] + v[6112] * v[692] - v[6114] * v[693];
		v[6142] = -(v[4874] * v[5126]) + v[4870] * v[5128] - v[4867] * v[5129] + v[4777] * v[5243] + v[238] * v[6103]
			- v[6119] * v[691] + v[6129] * v[692] - v[6123] * v[693];
		v[6143] = v[6113] + v[6124];
		v[6144] = -(v[4873] * v[5126]) + v[4871] * v[5128] - v[4869] * v[5129] + v[4786] * v[5243] + v[238] * v[6096]
			- v[6121] * v[691] + v[6130] * v[692] - v[6122] * v[693];
		v[6145] = -(v[4883] * v[5126]) + v[4880] * v[5128] - v[4876] * v[5129] + v[4773] * v[5243] + v[238] * v[6094]
			- v[6134] * v[691] + v[6138] * v[692] - v[6126] * v[693];
		v[6146] = -(v[4882] * v[5126]) + v[4879] * v[5128] - v[4877] * v[5129] + v[4790] * v[5243] + v[238] * v[6092]
			- v[6135] * v[691] + v[6137] * v[692] - v[6125] * v[693];
		v[6147] = v[6120] + v[6136];
		v[6148] = v[6110] + v[6139];
		v[6149] = v[4758] * v[5052] + v[228] * v[5962];
		v[6150] = v[4758] * v[5051] + v[227] * v[5962];
		v[6151] = v[4757] * v[5052] + v[228] * v[5964];
		v[6152] = v[4757] * v[5051] + v[227] * v[5964];
		v[6153] = v[4739] * v[5052] + v[228] * v[5950];
		v[6154] = v[4739] * v[5051] + v[227] * v[5950];
		v[6155] = v[4738] * v[5052] + v[228] * v[5952];
		v[6156] = v[4738] * v[5051] + v[227] * v[5952];
		v[6157] = v[4720] * v[5052] + v[228] * v[5938];
		v[6158] = v[4720] * v[5051] + v[227] * v[5938];
		v[6159] = v[4719] * v[5052] + v[228] * v[5940];
		v[6160] = v[4719] * v[5051] + v[227] * v[5940];
		v[6161] = -(v[4672] * v[5049]) + v[224] * v[5913];
		v[6162] = v[227] * v[6161] + v[5051] * v[7135];
		v[6163] = v[228] * v[6161] + v[5052] * v[7135];
		v[6164] = -(v[4669] * v[5049]) + v[224] * v[5914];
		v[6165] = v[227] * v[6164] + v[5051] * v[7136];
		v[6166] = v[228] * v[6164] + v[5052] * v[7136];
		v[6167] = -(v[4666] * v[5049]) + v[224] * v[5915];
		v[6168] = v[227] * v[6167] + v[5051] * v[7137];
		v[6169] = v[228] * v[6167] + v[5052] * v[7137];
		v[6170] = v[4681] * v[5049] + v[223] * v[5910];
		v[6171] = v[227] * v[6170] + v[5051] * v[7138];
		v[6172] = v[228] * v[6170] + v[5052] * v[7138];
		v[6173] = v[4678] * v[5049] + v[223] * v[5911];
		v[6174] = v[227] * v[6173] + v[5051] * v[7139];
		v[6175] = v[228] * v[6173] + v[5052] * v[7139];
		v[6176] = v[4675] * v[5049] + v[223] * v[5912];
		v[6177] = v[227] * v[6176] + v[5051] * v[7140];
		v[6178] = v[228] * v[6176] + v[5052] * v[7140];
		v[6179] = -(v[4673] * v[5049]) + v[224] * v[5901];
		v[6180] = v[227] * v[6179] + v[5051] * v[7141];
		v[6181] = v[228] * v[6179] + v[5052] * v[7141];
		v[6182] = -(v[4670] * v[5049]) + v[224] * v[5902];
		v[6183] = v[227] * v[6182] + v[5051] * v[7142];
		v[6184] = v[228] * v[6182] + v[5052] * v[7142];
		v[6185] = -(v[4667] * v[5049]) + v[224] * v[5903];
		v[6186] = v[227] * v[6185] + v[5051] * v[7143];
		v[6187] = v[228] * v[6185] + v[5052] * v[7143];
		v[6188] = v[4682] * v[5049] + v[223] * v[5898];
		v[6189] = v[227] * v[6188] + v[5051] * v[7144];
		v[6190] = v[228] * v[6188] + v[5052] * v[7144];
		v[6191] = v[4679] * v[5049] + v[223] * v[5899];
		v[6192] = v[227] * v[6191] + v[5051] * v[7145];
		v[6193] = v[228] * v[6191] + v[5052] * v[7145];
		v[6194] = v[4676] * v[5049] + v[223] * v[5900];
		v[6195] = v[227] * v[6194] + v[5051] * v[7146];
		v[6196] = v[228] * v[6194] + v[5052] * v[7146];
		v[6197] = -(v[4674] * v[5049]) + v[224] * v[5894];
		v[6198] = v[227] * v[6197] + v[5051] * v[7147];
		v[6199] = v[228] * v[6197] + v[5052] * v[7147];
		v[6200] = -(v[4671] * v[5049]) + v[224] * v[5893];
		v[6201] = v[227] * v[6200] + v[5051] * v[7148];
		v[6202] = v[228] * v[6200] + v[5052] * v[7148];
		v[6203] = -(v[4668] * v[5049]) + v[224] * v[5892];
		v[6204] = v[227] * v[6203] + v[5051] * v[7149];
		v[6205] = v[228] * v[6203] + v[5052] * v[7149];
		v[6206] = v[4683] * v[5049] + v[223] * v[5897];
		v[6207] = v[227] * v[6206] + v[5051] * v[7150];
		v[6208] = v[228] * v[6206] + v[5052] * v[7150];
		v[6209] = v[4680] * v[5049] + v[223] * v[5896];
		v[6210] = v[227] * v[6209] + v[5051] * v[7151];
		v[6211] = v[228] * v[6209] + v[5052] * v[7151];
		v[6212] = v[4677] * v[5049] + v[223] * v[5895];
		v[6213] = v[227] * v[6212] + v[5051] * v[7152];
		v[6214] = v[228] * v[6212] + v[5052] * v[7152];
		v[6215] = ((b6992 ? v[4523] : 0e0) + (b7005 ? -v[4523] : 0e0) + (i6990 ? v[4525] : 0e0) + (i6991 ? v[4524] : 0e0) + (i7003 ?
			-v[4525] : 0e0) + (i7004 ? -v[4524] : 0e0) - v[1469] * v[5404] - v[1470] * v[5405] - v[1471] * v[5406] + v[1465] * v[5407]
			+ v[1466] * v[5408] + v[1467] * v[5409] + (-v[420] + v[421])*v[5414] + (-v[423] + v[424])*v[5418] + (-v[426] + v[427]
				)*v[5425] + v[5913] * v[637] + v[5914] * v[638] + v[5915] * v[639] - v[5910] * v[640] - v[5911] * v[641]
			- v[5912] * v[642] + v[5901] * v[649] + (-v[5867] + v[5868])*v[6499] + v[5902] * v[650] + v[5903] * v[651]
			- v[5898] * v[652] - v[5899] * v[653] - v[5900] * v[654] + v[5894] * v[661] + v[5893] * v[662] + v[5892] * v[663]
			- v[5897] * v[664] - v[5896] * v[665] - v[5895] * v[666] + (v[12884 + i4521] - v[12902 + i4521])*v[7] - (v[4634]
				- v[4635])*(v[6980] + v[7016]) - v[1014] * (v[4533] * (v[5058] - v[5060]) + v[5739] * v[7153]) - v[1018] * (v[4534] *
				(v[5059] - v[5061]) + v[5740] * v[7154]) - (v[288] - v[291])*(v[5726] + v[7043] + v[7155]) - (v[286] - v[289])*
					(v[5713] + v[7156]) - (v[287] - v[290])*(v[5720] + v[7157]) + v[5051] * (v[7434] - v[7438]) + v[5052] * (v[7435]
						- v[7439]) + v[5050] * (v[7436] - v[7440])) / 2e0;
		v[6216] = v[2511] * v[6008] + v[337] * v[6222] + v[6019] * v[6468] + v[335] * (v[5350] / v[324] + v[5210] * v[7482])
			+ v[6160] * v[97] + v[6159] * v[98];
		v[6219] = (v[4802] * v[5324]) / v[324] + v[4509] * v[6000] + v[334] * v[6226] + v[5353] * v[6445] + v[6020] * v[6468]
			+ v[5210] * v[7483] + v[6160] * v[94] + v[6159] * v[95];
		v[6221] = v[4509] * v[6001] + v[322] * v[6229] + v[327] * (v[7221] + v[5210] * v[7484]) + v[6160] * v[91]
			+ v[6159] * v[92];
		v[6223] = (v[4801] * v[5327]) / v[324] + v[4507] * v[6008] + v[5361] * v[6446] + v[6023] * v[6465] + v[6222] * v[6469]
			+ v[5210] * v[7485] + v[6156] * v[97] + v[6155] * v[98];
		v[6227] = v[2513] * v[6000] + v[328] * v[6226] + v[6018] * v[6465] + v[333] * (v[5363] / v[324] + v[5210] * v[7486])
			+ v[6156] * v[94] + v[6155] * v[95];
		v[6230] = v[4507] * v[6009] + v[323] * v[6229] + v[327] * (v[7220] + v[5210] * v[7487]) + v[6156] * v[91]
			+ v[6155] * v[92];
		v[6231] = (v[4796] * v[5328]) / v[324] + v[2512] * v[6012] + v[4504] * v[6019] + v[5357] * v[6446] + v[6222] * v[6467]
			+ v[5210] * v[7488] + v[6152] * v[97] + v[6151] * v[98];
		v[6233] = v[4504] * v[6018] + v[6226] * v[6464] + v[333] * (v[5210] * v[7160] + v[7224]) + v[6152] * v[94]
			+ v[6151] * v[95];
		v[6236] = v[2514] * v[6001] + v[2512] * v[6009] + v[325] * v[6229] + v[327] * (v[5359] / v[324] + v[5210] * v[7489])
			+ v[6152] * v[91] + v[6151] * v[92];
		v[6237] = v[6201] * v[91] + v[6202] * v[92];
		v[6238] = v[6201] * v[94] + v[6202] * v[95];
		v[6239] = v[6201] * v[97] + v[6202] * v[98];
		v[6240] = v[6204] * v[91] + v[6205] * v[92];
		v[6241] = v[6204] * v[97] + v[6205] * v[98];
		v[6242] = v[6204] * v[94] + v[6205] * v[95];
		v[6243] = v[6198] * v[94] + v[6199] * v[95];
		v[6244] = v[6198] * v[97] + v[6199] * v[98];
		v[6245] = v[6198] * v[91] + v[6199] * v[92];
		v[6246] = v[6180] * v[91] + v[6181] * v[92];
		v[6247] = v[6180] * v[94] + v[6181] * v[95];
		v[6248] = v[6180] * v[97] + v[6181] * v[98];
		v[6249] = v[6186] * v[97] + v[6187] * v[98];
		v[6250] = v[6186] * v[91] + v[6187] * v[92];
		v[6251] = v[6186] * v[94] + v[6187] * v[95];
		v[6252] = v[6168] * v[94] + v[6169] * v[95];
		v[6253] = v[6168] * v[91] + v[6169] * v[92];
		v[6254] = v[6168] * v[97] + v[6169] * v[98];
		v[6255] = -(v[4798] * v[5108]) - v[4820] * v[5109] + 2e0*v[4797] * v[5110] - v[533] * v[6013] - v[530] * v[6026]
			+ v[6243] + v[6246] + v[6249] + v[6252] + v[6017] * v[6659];
		v[6256] = v[6183] * v[91] + v[6184] * v[92];
		v[6257] = v[6183] * v[97] + v[6184] * v[98];
		v[6258] = v[6183] * v[94] + v[6184] * v[95];
		v[6259] = v[4824] * v[5108] + 2e0*v[4800] * v[5109] - v[4820] * v[5110] - v[527] * v[6026] + v[533] * v[6029] - v[6238]
			- v[6241] - v[6253] - v[6256] + v[6022] * v[6658];
		v[6260] = -(v[4944] * v[5100]) + v[4940] * v[5102] - v[4936] * v[5103] + v[4813] * v[5187] + v[178] * v[6233]
			+ v[521] * v[6238] - v[522] * v[6242] - v[520] * v[6243];
		v[6261] = v[6162] * v[91] + v[6163] * v[92];
		v[6262] = v[6162] * v[94] + v[6163] * v[95];
		v[6263] = v[6162] * v[97] + v[6163] * v[98];
		v[6264] = v[6165] * v[94] + v[6166] * v[95];
		v[6265] = v[6165] * v[91] + v[6166] * v[92];
		v[6266] = v[6165] * v[97] + v[6166] * v[98];
		v[6267] = 2e0*v[4794] * v[5108] + v[4824] * v[5109] - v[4798] * v[5110] - v[527] * v[6013] + v[530] * v[6029] - v[6244]
			- v[6257] - v[6261] - v[6264] + v[6007] * v[6657];
		v[6268] = -(v[4943] * v[5100]) + v[4941] * v[5102] - v[4937] * v[5103] + v[4811] * v[5187] + v[178] * v[6231]
			+ v[521] * v[6239] - v[522] * v[6241] - v[520] * v[6244];
		v[6269] = -(v[4925] * v[5100]) + v[4921] * v[5102] - v[4918] * v[5103] + v[4810] * v[5187] + v[178] * v[6230]
			- v[520] * v[6246] - v[522] * v[6250] + v[521] * v[6256];
		v[6270] = v[6240] + v[6251];
		v[6271] = -(v[4924] * v[5100]) + v[4922] * v[5102] - v[4920] * v[5103] + v[4819] * v[5187] + v[178] * v[6223]
			- v[520] * v[6248] - v[522] * v[6249] + v[521] * v[6257];
		v[6272] = -(v[4907] * v[5100]) + v[4904] * v[5102] - v[4900] * v[5103] + v[4806] * v[5187] + v[178] * v[6221]
			- v[522] * v[6253] - v[520] * v[6261] + v[521] * v[6265];
		v[6273] = -(v[4906] * v[5100]) + v[4903] * v[5102] - v[4901] * v[5103] + v[4823] * v[5187] + v[178] * v[6219]
			- v[522] * v[6252] - v[520] * v[6262] + v[521] * v[6264];
		v[6274] = v[6247] + v[6263];
		v[6275] = v[6237] + v[6266];
		v[6276] = v[2503] * v[6038] + v[311] * v[6282] + v[6049] * v[6460] + v[309] * (v[5310] / v[298] + v[5172] * v[7490])
			+ v[6158] * v[88] + v[6157] * v[89];
		v[6279] = (v[4835] * v[5284]) / v[298] + v[4502] * v[6030] + v[308] * v[6286] + v[5313] * v[6438] + v[6050] * v[6460]
			+ v[5172] * v[7491] + v[6158] * v[85] + v[6157] * v[86];
		v[6281] = v[4502] * v[6031] + v[296] * v[6289] + v[301] * (v[7209] + v[5172] * v[7492]) + v[6158] * v[82]
			+ v[6157] * v[83];
		v[6283] = (v[4834] * v[5287]) / v[298] + v[4500] * v[6038] + v[5321] * v[6439] + v[6053] * v[6457] + v[6282] * v[6461]
			+ v[5172] * v[7493] + v[6154] * v[88] + v[6153] * v[89];
		v[6287] = v[2505] * v[6030] + v[302] * v[6286] + v[6048] * v[6457] + v[307] * (v[5323] / v[298] + v[5172] * v[7494])
			+ v[6154] * v[85] + v[6153] * v[86];
		v[6290] = v[4500] * v[6039] + v[297] * v[6289] + v[301] * (v[7208] + v[5172] * v[7495]) + v[6154] * v[82]
			+ v[6153] * v[83];
		v[6291] = (v[4829] * v[5288]) / v[298] + v[2504] * v[6042] + v[4497] * v[6049] + v[5317] * v[6439] + v[6282] * v[6459]
			+ v[5172] * v[7496] + v[6150] * v[88] + v[6149] * v[89];
		v[6293] = v[4497] * v[6048] + v[6286] * v[6456] + v[307] * (v[5172] * v[7163] + v[7212]) + v[6150] * v[85]
			+ v[6149] * v[86];
		v[6296] = v[2506] * v[6031] + v[2504] * v[6039] + v[299] * v[6289] + v[301] * (v[5319] / v[298] + v[5172] * v[7497])
			+ v[6150] * v[82] + v[6149] * v[83];
		v[6297] = v[6210] * v[82] + v[6211] * v[83];
		v[6298] = v[6210] * v[85] + v[6211] * v[86];
		v[6299] = v[6210] * v[88] + v[6211] * v[89];
		v[6300] = v[6213] * v[82] + v[6214] * v[83];
		v[6301] = v[6213] * v[88] + v[6214] * v[89];
		v[6302] = v[6213] * v[85] + v[6214] * v[86];
		v[6303] = v[6207] * v[85] + v[6208] * v[86];
		v[6304] = v[6207] * v[88] + v[6208] * v[89];
		v[6305] = v[6207] * v[82] + v[6208] * v[83];
		v[6306] = v[6189] * v[82] + v[6190] * v[83];
		v[6307] = v[6189] * v[85] + v[6190] * v[86];
		v[6308] = v[6189] * v[88] + v[6190] * v[89];
		v[6309] = v[6195] * v[88] + v[6196] * v[89];
		v[6310] = v[6195] * v[82] + v[6196] * v[83];
		v[6311] = v[6195] * v[85] + v[6196] * v[86];
		v[6312] = v[6177] * v[85] + v[6178] * v[86];
		v[6313] = v[6177] * v[82] + v[6178] * v[83];
		v[6314] = v[6177] * v[88] + v[6178] * v[89];
		v[6315] = -(v[4831] * v[5082]) - v[4853] * v[5083] + 2e0*v[4830] * v[5084] - v[488] * v[6043] - v[485] * v[6056]
			+ v[6303] + v[6306] + v[6309] + v[6312] + v[6047] * v[6655];
		v[6316] = v[6192] * v[82] + v[6193] * v[83];
		v[6317] = v[6192] * v[88] + v[6193] * v[89];
		v[6318] = v[6192] * v[85] + v[6193] * v[86];
		v[6319] = v[4857] * v[5082] + 2e0*v[4833] * v[5083] - v[4853] * v[5084] - v[482] * v[6056] + v[488] * v[6059] - v[6298]
			- v[6301] - v[6313] - v[6316] + v[6052] * v[6654];
		v[6320] = -(v[4953] * v[5074]) + v[4949] * v[5076] - v[4945] * v[5077] + v[4846] * v[5149] + v[159] * v[6293]
			+ v[476] * v[6298] - v[477] * v[6302] - v[475] * v[6303];
		v[6321] = v[6171] * v[82] + v[6172] * v[83];
		v[6322] = v[6171] * v[85] + v[6172] * v[86];
		v[6323] = v[6171] * v[88] + v[6172] * v[89];
		v[6324] = v[6174] * v[85] + v[6175] * v[86];
		v[6325] = v[6174] * v[82] + v[6175] * v[83];
		v[6326] = v[6174] * v[88] + v[6175] * v[89];
		v[6327] = 2e0*v[4827] * v[5082] + v[4857] * v[5083] - v[4831] * v[5084] - v[482] * v[6043] + v[485] * v[6059] - v[6304]
			- v[6317] - v[6321] - v[6324] + v[6037] * v[6653];
		v[6328] = -(v[4952] * v[5074]) + v[4950] * v[5076] - v[4946] * v[5077] + v[4844] * v[5149] + v[159] * v[6291]
			+ v[476] * v[6299] - v[477] * v[6301] - v[475] * v[6304];
		v[6329] = -(v[4934] * v[5074]) + v[4930] * v[5076] - v[4927] * v[5077] + v[4843] * v[5149] + v[159] * v[6290]
			- v[475] * v[6306] - v[477] * v[6310] + v[476] * v[6316];
		v[6330] = v[6300] + v[6311];
		v[6331] = -(v[4933] * v[5074]) + v[4931] * v[5076] - v[4929] * v[5077] + v[4852] * v[5149] + v[159] * v[6283]
			- v[475] * v[6308] - v[477] * v[6309] + v[476] * v[6317];
		v[6332] = -(v[4916] * v[5074]) + v[4913] * v[5076] - v[4909] * v[5077] + v[4839] * v[5149] + v[159] * v[6281]
			- v[477] * v[6313] - v[475] * v[6321] + v[476] * v[6325];
		v[6333] = -(v[4915] * v[5074]) + v[4912] * v[5076] - v[4910] * v[5077] + v[4856] * v[5149] + v[159] * v[6279]
			- v[477] * v[6312] - v[475] * v[6322] + v[476] * v[6324];
		v[6334] = v[6307] + v[6323];
		v[6335] = v[6297] + v[6326];
		v[6338] = v[1026] * (v[1611] * v[6087] - v[1025] * v[6088] - v[5042] * v[6337]) + v[1612] * (-(v[129] * v[5410])
			- v[132] * v[5411] - v[135] * v[5412] - v[258] * v[5414] - v[261] * v[5418] - v[264] * v[5425] + v[5042] * v[6336]
			+ v[5891] * v[739] + v[5890] * v[740] + v[5889] * v[741] + v[5904] * v[754] + v[5905] * v[755] + v[5906] * v[756]
			+ v[5916] * v[769] + v[5917] * v[770] + v[5918] * v[771]);
		v[6339] = v[1025] * (-(v[275] * v[6087]) + v[4884] * v[7164]) + v[1611] * (-(v[275] * v[6088]) + v[4885] * v[7164]);
		v[6340] = (-(v[4878] * v[5038]) - v[4860] * v[5041] - v[4868] * v[5127] - 2e0*v[5973] + 2e0*v[5981]
			+ v[6117] * v[6626] + v[6116] * v[6627] + v[6121] * v[6628] + v[6119] * v[6629] + v[6135] * v[6630] + v[6134] * v[6631]
			- v[6118] * v[694] - v[6120] * v[713] - v[4867] * v[7165] - v[4858] * v[7166] - v[4877] * v[7167] - v[4869] * v[7168]
			- v[4876] * v[7169] - v[4859] * v[7170] - v[6136] * v[731]) / 2e0;
		v[6341] = (-(v[4864] * v[5126]) + v[4861] * v[5128] - v[4860] * v[5129] + v[4783] * v[5243] + v[238] * v[6109]
			- v[6118] * v[691] + v[6110] * v[692] - v[6113] * v[693]) / 2e0;
		v[6343] = (v[4881] * v[5038] + v[4861] * v[5041] + v[4872] * v[5127] - 2e0*v[5975] + 2e0*v[5998] - v[6112] * v[6626]
			- v[6111] * v[6627] - v[6130] * v[6628] - v[6129] * v[6629] - v[6137] * v[6630] - v[6138] * v[6631] + v[6110] * v[694]
			+ v[6131] * v[713] + v[4870] * v[7165] + v[4862] * v[7166] + v[4879] * v[7167] + v[4871] * v[7168] + v[4880] * v[7169]
			+ v[4863] * v[7170] + v[6139] * v[731]) / 2e0;
		v[6344] = (-(v[4875] * v[5126]) + v[4872] * v[5128] - v[4868] * v[5129] + v[4774] * v[5243] + v[238] * v[6100]
			- v[6120] * v[691] + v[6131] * v[692] - v[6124] * v[693]) / 2e0;
		v[6345] = (-(v[4886] * v[5038]) - v[4864] * v[5041] - v[4875] * v[5127] - 2e0*v[5985] + 2e0*v[5995]
			+ v[6114] * v[6626] + v[6115] * v[6627] + v[6122] * v[6628] + v[6123] * v[6629] + v[6125] * v[6630] + v[6126] * v[6631]
			- v[6113] * v[694] - v[6124] * v[713] - v[4874] * v[7165] - v[4866] * v[7166] - v[4882] * v[7167] - v[4873] * v[7168]
			- v[4883] * v[7169] - v[4865] * v[7170] - v[6127] * v[731]) / 2e0;
		v[6421] = v[2637] * v[6343] - v[6344] + v[2634] * v[6345] + (v[12524 + i4521] + v[4899] * v[5039] - v[469] * v[6340]
			)*v[6431] + v[5039] * v[7173] * v[7499] - v[6432] * (v[12866 + i4521] + (v[4770] * v[5038]) / 2e0 + (v[4783] * v[5041])
				/ 2e0 + v[4777] * v[5117] + v[4780] * v[5119] + v[4790] * v[5120] + v[4786] * v[5122] + v[4773] * v[5123]
				+ v[4778] * v[5125] + (v[4774] * v[5127]) / 2e0 + v[6112] - v[6115] - v[6121] + v[6123] + v[470] * v[6128]
				- v[471] * v[6132] + v[6135] - v[6138] - v[473] * v[6140] + 2e0*v[5243] * v[6348] + v[6089] * v[6575] + v[6100] * v[6576]
				+ v[6109] * v[6577] + v[6143] * v[6837] + v[6147] * v[6838] + v[6148] * v[6839] + v[6106] * v[699] + v[6104] * v[705]
				+ v[6103] * v[709] + v[6096] * v[718] + v[6094] * v[722] + v[6092] * v[726] + v[6648] * (v[4077] * v[5386]
					+ v[4081] * v[5389] + v[4065] * v[5390] + v[4064] * v[5393] + v[4084] * v[5395] + v[4078] * v[5397] + v[4069] * v[5399]
					+ v[4072] * v[5401] + v[4071] * v[5403] + v[4432] * v[5970] + v[4433] * v[5971] + v[5977] + v[4434] * v[5978]
					+ v[4435] * v[5979] + v[2955] * v[5982] + v[5987] + v[4436] * v[5988] + v[4437] * v[5989] + v[2949] * v[5990] + v[5992]
					+ v[2953] * v[5993] + v[5384] * v[6090] + v[5382] * v[6091] + v[5380] * v[6093] + v[5379] * v[6097] + v[5376] * v[6098]
					+ v[5375] * v[6101] + v[5371] * v[6105] + v[5374] * v[6107] + v[5370] * v[6108] + v[5242] * v[6353] + v[5247] * v[6354]
					+ v[5253] * v[6355] + v[5364] * v[6356] + v[5367] * v[6357] + v[5368] * v[6358] + v[5871] * v[6840] + v[5872] * v[6841]
					+ v[5873] * v[6842] + v[5266] * v[7174] * v[7501]));
		v[7203] = v[6421] + (v[4886] * v[5126] - v[4881] * v[5128] + v[4878] * v[5129] - v[4770] * v[5243] - v[238] * v[6089]
			+ v[6136] * v[691] - v[6139] * v[692] + v[6127] * v[693]) / 2e0;
		v[6360] = v[6141] + v[6145];
		v[6361] = v[6144] + v[6146];
		v[6362] = v[6133] + v[6142];
		v[6365] = -(v[1024] * (v[203] * v[5938] + v[212] * v[5940] + v[200] * v[5950] + v[209] * v[5952] + v[197] * v[5962]
			+ v[206] * v[5964] + v[610] * v[6161] + v[611] * v[6164] + v[612] * v[6167] + v[583] * v[6170] + v[584] * v[6173]
			+ v[585] * v[6176] + v[601] * v[6179] + v[602] * v[6182] + v[603] * v[6185] + v[574] * v[6188] + v[575] * v[6191]
			+ v[576] * v[6194] + v[592] * v[6197] + v[593] * v[6200] + v[594] * v[6203] + v[565] * v[6206] + v[566] * v[6209]
			+ v[567] * v[6212] + v[5050] * v[6363] + v[7175] * v[82] + v[7176] * v[85] + v[7177] * v[88] + v[7178] * v[91]
			+ v[7179] * v[94] + v[7180] * v[97])) + v[1610] * (v[204] * v[5938] + v[213] * v[5940] + v[201] * v[5950]
				+ v[210] * v[5952] + v[198] * v[5962] + v[207] * v[5964] + v[613] * v[6161] + v[614] * v[6164] + v[615] * v[6167]
				+ v[586] * v[6170] + v[587] * v[6173] + v[588] * v[6176] + v[604] * v[6179] + v[605] * v[6182] + v[606] * v[6185]
				+ v[577] * v[6188] + v[578] * v[6191] + v[579] * v[6194] + v[595] * v[6197] + v[596] * v[6200] + v[597] * v[6203]
				+ v[568] * v[6206] + v[569] * v[6209] + v[570] * v[6212] + v[5050] * v[6364] + v[7175] * v[83] + v[7176] * v[86]
				+ v[7177] * v[89] + v[7178] * v[92] + v[7179] * v[95] + v[7180] * v[98]);
		v[6366] = (-(v[4902] * v[5030]) - v[4938] * v[5033] - v[4919] * v[5101] - 2e0*v[6003] + 2e0*v[6011] - v[523] * v[6245]
			- v[542] * v[6247] - v[560] * v[6263] + v[6262] * v[6632] + v[6261] * v[6633] + v[6248] * v[6634] + v[6246] * v[6635]
			+ v[6244] * v[6636] + v[6243] * v[6637] - v[4918] * v[7181] - v[4936] * v[7182] - v[4901] * v[7183] - v[4920] * v[7184]
			- v[4900] * v[7185] - v[4937] * v[7186]) / 2e0;
		v[6367] = (-(v[4942] * v[5100]) + v[4939] * v[5102] - v[4938] * v[5103] + v[4816] * v[5187] + v[178] * v[6236]
			+ v[521] * v[6237] - v[522] * v[6240] - v[520] * v[6245]) / 2e0;
		v[6369] = (v[4905] * v[5030] + v[4939] * v[5033] + v[4923] * v[5101] - 2e0*v[6005] + 2e0*v[6028] + v[523] * v[6237]
			+ v[542] * v[6258] + v[560] * v[6266] - v[6264] * v[6632] - v[6265] * v[6633] - v[6257] * v[6634] - v[6256] * v[6635]
			- v[6239] * v[6636] - v[6238] * v[6637] + v[4921] * v[7181] + v[4940] * v[7182] + v[4903] * v[7183] + v[4922] * v[7184]
			+ v[4904] * v[7185] + v[4941] * v[7186]) / 2e0;
		v[6370] = (-(v[4926] * v[5100]) + v[4923] * v[5102] - v[4919] * v[5103] + v[4807] * v[5187] + v[178] * v[6227]
			- v[520] * v[6247] - v[522] * v[6251] + v[521] * v[6258]) / 2e0;
		v[6371] = (-(v[4908] * v[5030]) - v[4942] * v[5033] - v[4926] * v[5101] - 2e0*v[6015] + 2e0*v[6025] - v[523] * v[6240]
			- v[542] * v[6251] - v[560] * v[6254] + v[6252] * v[6632] + v[6253] * v[6633] + v[6249] * v[6634] + v[6250] * v[6635]
			+ v[6241] * v[6636] + v[6242] * v[6637] - v[4925] * v[7181] - v[4944] * v[7182] - v[4906] * v[7183] - v[4924] * v[7184]
			- v[4907] * v[7185] - v[4943] * v[7186]) / 2e0;
		v[6417] = v[2611] * v[6369] - v[6370] + v[2608] * v[6371] + (v[12506 + i4521] + v[4968] * v[5031] - v[463] * v[6366]
			)*v[6429] + v[5031] * v[7189] * v[7503] - v[6430] * (v[12920 + i4521] + (v[4803] * v[5030]) / 2e0 + (v[4816] * v[5033])
				/ 2e0 + v[4810] * v[5091] + v[4813] * v[5093] + v[4823] * v[5094] + v[4819] * v[5096] + v[4806] * v[5097]
				+ v[4811] * v[5099] + (v[4807] * v[5101]) / 2e0 + v[555] * v[6219] + v[551] * v[6221] + v[547] * v[6223] + v[538] * v[6230]
				+ v[534] * v[6231] + v[528] * v[6233] + v[6239] - v[6242] - v[6248] + v[6250] + v[464] * v[6255] - v[465] * v[6259]
				+ v[6262] - v[6265] - v[467] * v[6267] + 2e0*v[5187] * v[6374] + v[6216] * v[6605] + v[6227] * v[6606]
				+ v[6236] * v[6607] + v[6270] * v[6859] + v[6274] * v[6860] + v[6275] * v[6861] + v[6647] * (v[4107] * v[5346]
					+ v[4111] * v[5349] + v[4095] * v[5350] + v[4094] * v[5353] + v[4114] * v[5355] + v[4108] * v[5357] + v[4099] * v[5359]
					+ v[4102] * v[5361] + v[4101] * v[5363] + v[4460] * v[6000] + v[4461] * v[6001] + v[6007] + v[4462] * v[6008]
					+ v[4463] * v[6009] + v[2927] * v[6012] + v[6017] + v[4464] * v[6018] + v[4465] * v[6019] + v[2921] * v[6020] + v[6022]
					+ v[2925] * v[6023] + v[5344] * v[6217] + v[5342] * v[6218] + v[5340] * v[6220] + v[5339] * v[6224] + v[5336] * v[6225]
					+ v[5335] * v[6228] + v[5331] * v[6232] + v[5334] * v[6234] + v[5330] * v[6235] + v[5186] * v[6379] + v[5191] * v[6380]
					+ v[5197] * v[6381] + v[5324] * v[6382] + v[5327] * v[6383] + v[5328] * v[6384] + v[5874] * v[6862] + v[5875] * v[6863]
					+ v[5876] * v[6864] + v[5210] * v[7190] * v[7505]));
		v[7202] = (v[4908] * v[5100] - v[4905] * v[5102] + v[4902] * v[5103] - v[4803] * v[5187] - v[178] * v[6216]
			+ v[522] * v[6254] + v[520] * v[6263] - v[521] * v[6266]) / 2e0 + v[6417];
		v[6386] = v[6268] + v[6272];
		v[6387] = v[6271] + v[6273];
		v[6388] = v[6260] + v[6269];
		v[6389] = (-(v[4911] * v[5022]) - v[4947] * v[5025] - v[4928] * v[5075] - 2e0*v[6033] + 2e0*v[6041] - v[478] * v[6305]
			- v[497] * v[6307] - v[515] * v[6323] + v[6322] * v[6638] + v[6321] * v[6639] + v[6308] * v[6640] + v[6306] * v[6641]
			+ v[6304] * v[6642] + v[6303] * v[6643] - v[4927] * v[7191] - v[4945] * v[7192] - v[4910] * v[7193] - v[4929] * v[7194]
			- v[4909] * v[7195] - v[4946] * v[7196]) / 2e0;
		v[6390] = (-(v[4951] * v[5074]) + v[4948] * v[5076] - v[4947] * v[5077] + v[4849] * v[5149] + v[159] * v[6296]
			+ v[476] * v[6297] - v[477] * v[6300] - v[475] * v[6305]) / 2e0;
		v[6392] = (v[4914] * v[5022] + v[4948] * v[5025] + v[4932] * v[5075] - 2e0*v[6035] + 2e0*v[6058] + v[478] * v[6297]
			+ v[497] * v[6318] + v[515] * v[6326] - v[6324] * v[6638] - v[6325] * v[6639] - v[6317] * v[6640] - v[6316] * v[6641]
			- v[6299] * v[6642] - v[6298] * v[6643] + v[4930] * v[7191] + v[4949] * v[7192] + v[4912] * v[7193] + v[4931] * v[7194]
			+ v[4913] * v[7195] + v[4950] * v[7196]) / 2e0;
		v[6393] = (-(v[4935] * v[5074]) + v[4932] * v[5076] - v[4928] * v[5077] + v[4840] * v[5149] + v[159] * v[6287]
			- v[475] * v[6307] - v[477] * v[6311] + v[476] * v[6318]) / 2e0;
		v[6394] = (-(v[4917] * v[5022]) - v[4951] * v[5025] - v[4935] * v[5075] - 2e0*v[6045] + 2e0*v[6055] - v[478] * v[6300]
			- v[497] * v[6311] - v[515] * v[6314] + v[6312] * v[6638] + v[6313] * v[6639] + v[6309] * v[6640] + v[6310] * v[6641]
			+ v[6301] * v[6642] + v[6302] * v[6643] - v[4934] * v[7191] - v[4953] * v[7192] - v[4915] * v[7193] - v[4933] * v[7194]
			- v[4916] * v[7195] - v[4952] * v[7196]) / 2e0;
		v[6413] = v[2585] * v[6392] - v[6393] + v[2582] * v[6394] + (v[12488 + i4521] + v[4981] * v[5023] - v[457] * v[6389]
			)*v[6427] + v[5023] * v[7199] * v[7507] - v[6428] * (v[12938 + i4521] + (v[4836] * v[5022]) / 2e0 + (v[4849] * v[5025])
				/ 2e0 + v[4843] * v[5065] + v[4846] * v[5067] + v[4856] * v[5068] + v[4852] * v[5070] + v[4839] * v[5071]
				+ v[4844] * v[5073] + (v[4840] * v[5075]) / 2e0 + v[510] * v[6279] + v[506] * v[6281] + v[502] * v[6283] + v[493] * v[6290]
				+ v[489] * v[6291] + v[483] * v[6293] + v[6299] - v[6302] - v[6308] + v[6310] + v[458] * v[6315] - v[459] * v[6319]
				+ v[6322] - v[6325] - v[461] * v[6327] + 2e0*v[5149] * v[6397] + v[6276] * v[6623] + v[6287] * v[6624]
				+ v[6296] * v[6625] + v[6330] * v[6873] + v[6334] * v[6874] + v[6335] * v[6875] + v[6646] * (v[4137] * v[5306]
					+ v[4141] * v[5309] + v[4125] * v[5310] + v[4124] * v[5313] + v[4144] * v[5315] + v[4138] * v[5317] + v[4129] * v[5319]
					+ v[4132] * v[5321] + v[4131] * v[5323] + v[4487] * v[6030] + v[4488] * v[6031] + v[6037] + v[4489] * v[6038]
					+ v[4490] * v[6039] + v[2899] * v[6042] + v[6047] + v[4491] * v[6048] + v[4492] * v[6049] + v[2893] * v[6050] + v[6052]
					+ v[2897] * v[6053] + v[5304] * v[6277] + v[5302] * v[6278] + v[5300] * v[6280] + v[5299] * v[6284] + v[5296] * v[6285]
					+ v[5295] * v[6288] + v[5291] * v[6292] + v[5294] * v[6294] + v[5290] * v[6295] + v[5148] * v[6402] + v[5153] * v[6403]
					+ v[5159] * v[6404] + v[5284] * v[6405] + v[5287] * v[6406] + v[5288] * v[6407] + v[5877] * v[6876] + v[5878] * v[6877]
					+ v[5879] * v[6878] + v[5172] * v[7200] * v[7509]));
		v[7201] = (v[4917] * v[5074] - v[4914] * v[5076] + v[4911] * v[5077] - v[4836] * v[5149] - v[159] * v[6276]
			+ v[477] * v[6314] + v[475] * v[6323] - v[476] * v[6326]) / 2e0 + v[6413];
		v[6409] = v[6328] + v[6332];
		v[6410] = v[6331] + v[6333];
		v[6411] = v[6320] + v[6329];
		v[13263] = 0e0;
		v[13264] = 0e0;
		v[13265] = 0e0;
		v[13266] = 0e0;
		v[13267] = v[5010];
		v[13268] = v[5008];
		v[13269] = 0e0;
		v[13270] = 0e0;
		v[13271] = 0e0;
		v[13272] = 0e0;
		v[13273] = 0e0;
		v[13274] = 0e0;
		v[13275] = 0e0;
		v[13276] = 0e0;
		v[13277] = 0e0;
		v[13278] = 0e0;
		v[13279] = 0e0;
		v[13280] = 0e0;
		v[13209] = 0e0;
		v[13210] = 0e0;
		v[13211] = 0e0;
		v[13212] = v[5010];
		v[13213] = 0e0;
		v[13214] = v[5009];
		v[13215] = 0e0;
		v[13216] = 0e0;
		v[13217] = 0e0;
		v[13218] = 0e0;
		v[13219] = 0e0;
		v[13220] = 0e0;
		v[13221] = 0e0;
		v[13222] = 0e0;
		v[13223] = 0e0;
		v[13224] = 0e0;
		v[13225] = 0e0;
		v[13226] = 0e0;
		v[13173] = 0e0;
		v[13174] = 0e0;
		v[13175] = 0e0;
		v[13176] = v[5008];
		v[13177] = v[5009];
		v[13178] = 0e0;
		v[13179] = 0e0;
		v[13180] = 0e0;
		v[13181] = 0e0;
		v[13182] = 0e0;
		v[13183] = 0e0;
		v[13184] = 0e0;
		v[13185] = 0e0;
		v[13186] = 0e0;
		v[13187] = 0e0;
		v[13188] = 0e0;
		v[13189] = 0e0;
		v[13190] = 0e0;
		v[13155] = 0e0;
		v[13156] = 0e0;
		v[13157] = 0e0;
		v[13158] = 0e0;
		v[13159] = 0e0;
		v[13160] = 0e0;
		v[13161] = 0e0;
		v[13162] = 0e0;
		v[13163] = 0e0;
		v[13164] = 0e0;
		v[13165] = v[5001];
		v[13166] = v[4999];
		v[13167] = 0e0;
		v[13168] = 0e0;
		v[13169] = 0e0;
		v[13170] = 0e0;
		v[13171] = 0e0;
		v[13172] = 0e0;
		v[13101] = 0e0;
		v[13102] = 0e0;
		v[13103] = 0e0;
		v[13104] = 0e0;
		v[13105] = 0e0;
		v[13106] = 0e0;
		v[13107] = 0e0;
		v[13108] = 0e0;
		v[13109] = 0e0;
		v[13110] = v[5001];
		v[13111] = 0e0;
		v[13112] = v[5000];
		v[13113] = 0e0;
		v[13114] = 0e0;
		v[13115] = 0e0;
		v[13116] = 0e0;
		v[13117] = 0e0;
		v[13118] = 0e0;
		v[13065] = 0e0;
		v[13066] = 0e0;
		v[13067] = 0e0;
		v[13068] = 0e0;
		v[13069] = 0e0;
		v[13070] = 0e0;
		v[13071] = 0e0;
		v[13072] = 0e0;
		v[13073] = 0e0;
		v[13074] = v[4999];
		v[13075] = v[5000];
		v[13076] = 0e0;
		v[13077] = 0e0;
		v[13078] = 0e0;
		v[13079] = 0e0;
		v[13080] = 0e0;
		v[13081] = 0e0;
		v[13082] = 0e0;
		v[13047] = 0e0;
		v[13048] = 0e0;
		v[13049] = 0e0;
		v[13050] = 0e0;
		v[13051] = 0e0;
		v[13052] = 0e0;
		v[13053] = 0e0;
		v[13054] = 0e0;
		v[13055] = 0e0;
		v[13056] = 0e0;
		v[13057] = 0e0;
		v[13058] = 0e0;
		v[13059] = 0e0;
		v[13060] = 0e0;
		v[13061] = 0e0;
		v[13062] = 0e0;
		v[13063] = v[4992];
		v[13064] = v[4990];
		v[12993] = 0e0;
		v[12994] = 0e0;
		v[12995] = 0e0;
		v[12996] = 0e0;
		v[12997] = 0e0;
		v[12998] = 0e0;
		v[12999] = 0e0;
		v[13000] = 0e0;
		v[13001] = 0e0;
		v[13002] = 0e0;
		v[13003] = 0e0;
		v[13004] = 0e0;
		v[13005] = 0e0;
		v[13006] = 0e0;
		v[13007] = 0e0;
		v[13008] = v[4992];
		v[13009] = 0e0;
		v[13010] = v[4991];
		v[12957] = 0e0;
		v[12958] = 0e0;
		v[12959] = 0e0;
		v[12960] = 0e0;
		v[12961] = 0e0;
		v[12962] = 0e0;
		v[12963] = 0e0;
		v[12964] = 0e0;
		v[12965] = 0e0;
		v[12966] = 0e0;
		v[12967] = 0e0;
		v[12968] = 0e0;
		v[12969] = 0e0;
		v[12970] = 0e0;
		v[12971] = 0e0;
		v[12972] = v[4990];
		v[12973] = v[4991];
		v[12974] = 0e0;
		v[13281] = v[5962];
		v[13282] = v[5950];
		v[13283] = v[5938];
		v[13284] = (v[13262 + i4521] - v[4970] * v[5149] - v[159] * v[6319]) / 2e0 - v[6331] + v[6333] + v[461] * v[6409]
			+ v[458] * v[6411] + 2e0*(v[5069] * v[6412] + v[6389] * v[6428] + v[6334] * v[6434] + v[5023] * (-(v[5002] * v[6427])
				+ v[4979] * v[6885])) + v[457] * v[7201];
		v[13285] = (v[13208 + i4521] + v[4969] * v[5149] + v[159] * v[6315]) / 2e0 + v[6328] - v[6332] + v[461] * v[6410]
			+ v[459] * v[6411] + 2e0*(v[5072] * v[6414] - v[6392] * v[6428] + v[6335] * v[6434] + v[5023] * (v[5004] * v[6427]
				+ v[4980] * v[6885])) + v[460] * (-v[6390] + v[6393] + v[7201]);
		v[13286] = -v[6320] + (v[13172 + i4521] - v[4972] * v[5149] - v[159] * v[6327]) / 2e0 + v[6329] + v[459] * v[6409]
			+ v[458] * v[6410] + v[462] * (-v[6390] + v[6413]) + 2e0*(v[5066] * v[6415] + v[6394] * v[6428] + v[6330] * v[6434]
				+ v[5023] * (-(v[5006] * v[6427]) + v[4975] * v[6885]));
		v[13287] = v[5964];
		v[13288] = v[5952];
		v[13289] = v[5940];
		v[13290] = (v[13154 + i4521] - v[4957] * v[5187] - v[178] * v[6259]) / 2e0 - v[6271] + v[6273] + v[467] * v[6386]
			+ v[464] * v[6388] + 2e0*(v[5095] * v[6416] + v[6366] * v[6430] + v[6274] * v[6441] + v[5031] * (-(v[4993] * v[6429])
				+ v[4966] * v[6899])) + v[463] * v[7202];
		v[13291] = (v[13100 + i4521] + v[4956] * v[5187] + v[178] * v[6255]) / 2e0 + v[6268] - v[6272] + v[467] * v[6387]
			+ v[465] * v[6388] + 2e0*(v[5098] * v[6418] - v[6369] * v[6430] + v[6275] * v[6441] + v[5031] * (v[4995] * v[6429]
				+ v[4967] * v[6899])) + v[466] * (-v[6367] + v[6370] + v[7202]);
		v[13292] = -v[6260] + (v[13064 + i4521] - v[4959] * v[5187] - v[178] * v[6267]) / 2e0 + v[6269] + v[465] * v[6386]
			+ v[464] * v[6387] + v[468] * (-v[6367] + v[6417]) + 2e0*(v[5092] * v[6419] + v[6371] * v[6430] + v[6270] * v[6441]
				+ v[5031] * (-(v[4997] * v[6429]) + v[4962] * v[6899]));
		v[13293] = -v[5414];
		v[13294] = -v[5418];
		v[13295] = -v[5425];
		v[13296] = (v[13046 + i4521] - v[4888] * v[5243] - v[238] * v[6132]) / 2e0 - v[6144] + v[6146] + v[473] * v[6360]
			+ v[470] * v[6362] + 2e0*(v[5121] * v[6420] + v[6340] * v[6432] + v[6147] * v[6448] + v[5039] * (-(v[4984] * v[6431])
				+ v[4897] * v[6913])) + v[469] * v[7203];
		v[13297] = (v[12992 + i4521] + v[4887] * v[5243] + v[238] * v[6128]) / 2e0 + v[6141] - v[6145] + v[473] * v[6361]
			+ v[471] * v[6362] + 2e0*(v[5124] * v[6422] - v[6343] * v[6432] + v[6148] * v[6448] + v[5039] * (v[4986] * v[6431]
				+ v[4898] * v[6913])) + v[472] * (-v[6341] + v[6344] + v[7203]);
		v[13298] = -v[6133] + (v[12956 + i4521] - v[4890] * v[5243] - v[238] * v[6140]) / 2e0 + v[6142] + v[471] * v[6360]
			+ v[470] * v[6361] + v[474] * (-v[6341] + v[6421]) + 2e0*(v[5118] * v[6423] + v[6345] * v[6432] + v[6143] * v[6448]
				+ v[5039] * (-(v[4988] * v[6431]) + v[4893] * v[6913]));
		Rc[i4521 - 1] += v[12196 + i4521] + v[4983] * v[5042] + v[4982] * v[5043] + v[4955] * v[5048] + v[4954] * v[5050]
			+ v[12214 + i4521] * v[7];
		for (i5016 = 1; i5016 <= 18; i5016++) {
			Kc[i4521 - 1][i5016 - 1] += v[13280 + i5016] + v[13298 + i5016] * v[7] + v[6338] * v[9120 + i5016] + v[6339] * v[9138
				+ i5016] + v[6215] * v[9156 + i5016] + v[6365] * v[9174 + i5016];
		};/* end for */
	};/* end for */

	delete[]v;

	//db.PrintPtr(Rc, nDOF);
}

